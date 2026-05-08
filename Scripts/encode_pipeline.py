#!/usr/bin/env python3
"""ENCODE bulk-genomics pipeline.

End-to-end analysis surface for the major ENCODE bulk assays:

  ChIP-seq (TF + Histone), ATAC-seq, DNase-seq, Hi-C, capture Hi-C,
  ChIA-PET, RNA-seq, MNase-seq, FAIRE-seq, RAMPAGE, CAGE.

Subcommands

  retrieve         Search the ENCODE Portal for experiments matching
                   assay / biosample / target / assembly filters.
  manifest         Hydrate one or more accessions and emit a per-file
                   manifest (signal bigWig, peaks BED, alignments BAM,
                   contact matrices mcool/.hic, ChIA-PET clusters bedpe).
  download         Pull manifest files under a size cap, filter by
                   ``content_type`` / ``file_format``.
  describe         Per-accession plain-language report: assay, biosample,
                   target, ENCODE pipeline, replicates, file inventory,
                   quality metrics if present.
  analyze-peaks    Peak-level QC from a BED file (count, width, score
                   distribution, density per chromosome).
  super-enhancers  ROSE-style super-enhancer calling: stitch H3K27ac (or
                   any enhancer-mark) peaks within --stitching-distance,
                   optionally exclude TSS-proximal stretches, rank by
                   signal, find the inflection point, emit SE BED +
                   hockey-stick plot.
  integrate-ccre   Overlay peaks with SCREEN cCRE classes (PLS / pELS /
                   dELS / CTCF / DNase-H3K4me3) and emit per-class
                   counts + a stacked-bar overview plot.
  browser          Render an IGV-style multi-track SVG for a genomic
                   region from any combination of BED / bedpe / BedGraph
                   tracks.
  write-playbook   Emit ``Docs/Skills/ENCODE_PIPELINE_SKILLS.md``.

All endpoint URLs resolve through ``Scripts/_endpoints.py``; cCRE
downloads route through the existing SCREEN catalog. Outputs follow the
project layout:

  Data/IGVF/ENCODE/Metadata/         hydrated portal JSON
  Data/IGVF/ENCODE/Downloads/        per-experiment file bundles
  Data/Manifests/ENCODE/             dataset / per-file manifests
  Data/Cache/ENCODE/                 cCRE BED caches, etc.
  Docs/ENCODE/                       per-run reports + plots
  Docs/Logs/                         runtime logs
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint, host as _resolve_host

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "ENCODE"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "ENCODE"
METADATA_DIR = DATA_DIR / "IGVF" / "ENCODE" / "Metadata"
DOWNLOAD_DIR = DATA_DIR / "IGVF" / "ENCODE" / "Downloads"
CACHE_DIR = DATA_DIR / "Cache" / "ENCODE"

ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
WENGLAB_DL = _resolve_endpoint("wenglab_dl")

USER_AGENT = "IGVFdataAgent-ENCODE/0.1"


# --------------------------- Assay taxonomy ---------------------------------

# Each entry: assay group label -> ENCODE assay_title facets, file
# preferences, and post-processing hints. Keys are case-insensitive
# canonical names users pass to --assay.
ENCODE_ASSAYS: "dict[str, dict[str, Any]]" = {
    "ChIP-seq": {
        "facets":          ["TF ChIP-seq", "ChIP-seq"],
        "needs_target":    True,
        "primary_files":   ["bigWig", "bed narrowPeak", "bigBed narrowPeak"],
        "key_metrics":     ["FRiP", "NRF", "PBC1", "PBC2", "NSC", "RSC"],
    },
    "Histone ChIP-seq": {
        "facets":          ["Histone ChIP-seq"],
        "needs_target":    True,
        "primary_files":   ["bigWig", "bed broadPeak", "bed narrowPeak"],
        "key_metrics":     ["FRiP"],
    },
    "ATAC-seq": {
        "facets":          ["ATAC-seq"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed narrowPeak"],
        "key_metrics":     ["TSS enrichment", "FRiP"],
    },
    "DNase-seq": {
        "facets":          ["DNase-seq"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed narrowPeak"],
        "key_metrics":     ["FRiP", "SPOT score"],
    },
    "Hi-C": {
        "facets":          ["Hi-C", "intact Hi-C", "in situ Hi-C"],
        "needs_target":    False,
        "primary_files":   ["hic", "mcool", "bed bedpe", "bedpe"],
        "key_metrics":     ["compartment AB", "TADs", "loops"],
    },
    "capture Hi-C": {
        "facets":          ["capture Hi-C"],
        "needs_target":    False,
        "primary_files":   ["hic", "mcool", "bedpe"],
        "key_metrics":     ["loops", "interactions"],
    },
    "ChIA-PET": {
        "facets":          ["ChIA-PET"],
        "needs_target":    True,
        "primary_files":   ["bedpe", "bigBed bedpe", "bed"],
        "key_metrics":     ["valid pairs", "intra-loops"],
    },
    "RNA-seq": {
        "facets":          ["polyA plus RNA-seq", "total RNA-seq",
                             "small RNA-seq", "RNA-seq"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "tsv"],
    },
    "MNase-seq": {
        "facets":          ["MNase-seq"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed narrowPeak"],
    },
    "FAIRE-seq": {
        "facets":          ["FAIRE-seq"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed narrowPeak"],
    },
    "CAGE": {
        "facets":          ["CAGE"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed"],
    },
    "RAMPAGE": {
        "facets":          ["RAMPAGE"],
        "needs_target":    False,
        "primary_files":   ["bigWig", "bed"],
    },
}

HISTONE_MARK_NOTES = {
    "H3K27ac":  "active enhancer / super-enhancer marker",
    "H3K4me1":  "primed enhancer marker (often paired with H3K27ac)",
    "H3K4me3":  "active promoter marker",
    "H3K27me3": "polycomb-repressed (Polycomb Repressive Complex 2 mark)",
    "H3K9me3":  "constitutive heterochromatin",
    "H3K36me3": "gene body / active transcription elongation",
    "H3K9ac":   "active gene marker",
    "H3K79me3": "active gene body / DOT1L",
    "H4K20me1": "transcribed gene body",
    "H2AFZ":    "active promoter / enhancer (H2A.Z)",
    "H3K4me2":  "active promoter / enhancer mixture",
}

# SCREEN cCRE class palette (matches V4 categories used by the catalog).
CCRE_CLASSES = {
    "PLS":       {"label": "Promoter-like (PLS)",      "color": "#D55E00"},
    "pELS":      {"label": "Proximal enhancer (pELS)", "color": "#0072B2"},
    "dELS":      {"label": "Distal enhancer (dELS)",   "color": "#56B4E9"},
    "CTCF-only": {"label": "CTCF-only",                "color": "#7B5BA6"},
    "DNase-H3K4me3": {"label": "DNase + H3K4me3",      "color": "#009E73"},
    "Low-DNase": {"label": "Low DNase",                "color": "#999999"},
}


# --------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"encode_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (DATA_DIR, REPORT_DIR, PLOT_DIR, MANIFEST_DIR, METADATA_DIR,
              DOWNLOAD_DIR, CACHE_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# --------------------------- HTTP helpers ----------------------------------

def request_headers(json_only: bool = True) -> "dict[str, str]":
    h = {"User-Agent": USER_AGENT,
         "Accept": "application/json,*/*" if json_only else "*/*"}
    return h


def fetch_json(url: str, timeout: int = 60) -> "tuple[int, Any]":
    logging.info("GET %s", url)
    req = urllib.request.Request(url, headers=request_headers(),
                                    method="GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace") if e.fp else ""
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            data = {"http_error_body": body[:600]}
        return e.code, data
    except urllib.error.URLError as e:
        return 0, {"network_error": str(e.reason)}


def encode_get(path: str, **params) -> "tuple[int, Any]":
    url = ENCODE_BASE + path
    if params:
        url = (url + ("&" if "?" in url else "?")
               + urllib.parse.urlencode({k: v for k, v in params.items()
                                            if v is not None}, doseq=True))
    return fetch_json(url)


def download_file(url: str, dest: Path, max_bytes: Optional[int] = None) -> int:
    logging.info("Download %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                  "Accept": "*/*"})
    written = 0
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            written += len(chunk)
            if max_bytes and written > max_bytes:
                f.close()
                dest.unlink(missing_ok=True)
                raise RuntimeError(f"Exceeded size cap on {url}")
    return written


# --------------------------- Manifest writing ------------------------------

def write_csv(path: Path, rows: "list[dict]",
                cols: "Optional[list[str]]" = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --------------------------- Portal queries --------------------------------

def search_experiments(assay: str, *, biosample: Optional[str] = None,
                         target: Optional[str] = None,
                         assembly: Optional[str] = None,
                         status: str = "released",
                         limit: int = 50) -> "list[dict]":
    """Search ENCODE for experiments matching the requested assay group."""
    cfg = _assay_cfg(assay)
    rows: "list[dict]" = []
    for facet in cfg["facets"]:
        params: "dict[str, Any]" = {
            "type":         "Experiment",
            "format":       "json",
            "limit":        limit,
            "status":       status,
            "assay_title":  facet,
        }
        if biosample:
            params["biosample_ontology.term_name"] = biosample
        if target and cfg.get("needs_target"):
            params["target.label"] = target
        if assembly:
            params["assembly"] = assembly
        sc, data = encode_get("/search/", **params)
        if sc != 200 or not isinstance(data, dict):
            logging.warning("search %s -> %s", facet, sc)
            continue
        for g in (data.get("@graph") or []):
            if isinstance(g, dict):
                g["_assay_group"] = assay
                rows.append(g)
        time.sleep(0.1)
    return rows


def hydrate_experiment(accession: str) -> "dict":
    sc, data = encode_get(f"/experiments/{accession}/", format="json")
    if sc == 200 and isinstance(data, dict):
        return data
    return {"accession": accession, "_http_status": sc}


def files_for_experiment(meta: dict) -> "list[dict]":
    out: "list[dict]" = []
    for f in (meta.get("files") or []):
        if not isinstance(f, dict):
            continue
        href = f.get("href") or ""
        out.append({
            "experiment":     meta.get("accession", ""),
            "file_accession": f.get("accession", ""),
            "file_format":    f.get("file_format", ""),
            "file_type":      f.get("file_type", ""),
            "output_type":    f.get("output_type", ""),
            "output_category": f.get("output_category", ""),
            "assembly":       f.get("assembly", ""),
            "genome_annotation": f.get("genome_annotation", ""),
            "biological_replicates":
                ",".join(str(x) for x in (f.get("biological_replicates") or [])),
            "file_size_bytes": f.get("file_size") or 0,
            "file_size_gb":   round((f.get("file_size") or 0) / (1024 ** 3), 3),
            "md5sum":         f.get("md5sum", ""),
            "status":         f.get("status", ""),
            "url":            (ENCODE_BASE + href) if href else "",
        })
    return out


def _assay_cfg(assay: str) -> dict:
    if assay in ENCODE_ASSAYS:
        return ENCODE_ASSAYS[assay]
    # case-insensitive fallback
    for k, v in ENCODE_ASSAYS.items():
        if k.lower() == assay.lower():
            return v
    raise SystemExit(
        f"Unknown assay: {assay!r}. Known: " + ", ".join(ENCODE_ASSAYS.keys())
    )


def summarize_experiment(meta: dict) -> dict:
    bio = meta.get("biosample_ontology") or {}
    if isinstance(bio, list):
        bio = bio[0] if bio else {}
    target = meta.get("target") or {}
    if isinstance(target, list):
        target = target[0] if target else {}
    lab = meta.get("lab") or {}
    return {
        "accession":        meta.get("accession", ""),
        "assay_title":      meta.get("assay_title", ""),
        "assay_term_name":  meta.get("assay_term_name", ""),
        "biosample":        bio.get("term_name", "") if isinstance(bio, dict) else "",
        "biosample_classification":
            bio.get("classification", "") if isinstance(bio, dict) else "",
        "target":           (target.get("label", "")
                              if isinstance(target, dict) else ""),
        "assembly":         ", ".join(meta.get("assembly", []) or []),
        "lab":              (lab.get("title", "")
                              if isinstance(lab, dict) else str(lab)),
        "status":           meta.get("status", ""),
        "n_files":          len(meta.get("files") or []),
        "description":      (meta.get("description") or "")[:240]
                              .replace("\n", " "),
        "url":              ENCODE_BASE + (meta.get("@id") or ""),
        "_assay_group":     meta.get("_assay_group", ""),
    }


# --------------------------- Subcommands -----------------------------------

def cmd_retrieve(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    rows = search_experiments(
        args.assay, biosample=args.biosample, target=args.target,
        assembly=args.assembly, status=args.status, limit=args.limit,
    )
    if args.fetch_file_details:
        for r in rows:
            acc = r.get("accession")
            if not acc:
                continue
            meta = hydrate_experiment(acc)
            (METADATA_DIR / f"{acc}.json").write_text(
                json.dumps(meta, indent=2, default=str))
            r.update(meta)
            time.sleep(0.1)

    summaries = [summarize_experiment(r) for r in rows]
    label = safe_label(args.label or f"{args.assay}_{args.biosample or 'any'}")
    out = MANIFEST_DIR / f"{ts}_{label}_experiments.csv"
    cols = ["accession", "assay_title", "assay_term_name", "_assay_group",
            "biosample", "biosample_classification", "target", "assembly",
            "lab", "status", "n_files", "description", "url"]
    write_csv(out, summaries, cols)
    print(f"Manifest:    {out}")
    print(f"Experiments: {len(summaries)}")
    return out


def cmd_manifest(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    accs = [a.strip() for a in args.accessions.split(",") if a.strip()]
    label = safe_label(args.label or (accs[0] if accs else "manifest"))
    summaries: "list[dict]" = []
    files_rows: "list[dict]" = []
    for acc in accs:
        meta = hydrate_experiment(acc)
        (METADATA_DIR / f"{acc}.json").write_text(
            json.dumps(meta, indent=2, default=str))
        summaries.append(summarize_experiment(meta))
        files_rows.extend(files_for_experiment(meta))
    summary_path = MANIFEST_DIR / f"{ts}_{label}_summary.csv"
    files_path = MANIFEST_DIR / f"{ts}_{label}_files.csv"
    write_csv(summary_path, summaries)
    write_csv(files_path, files_rows,
              cols=["experiment", "file_accession", "file_format",
                     "file_type", "output_type", "output_category",
                     "assembly", "genome_annotation",
                     "biological_replicates", "file_size_gb",
                     "file_size_bytes", "md5sum", "status", "url"])
    print(f"Summary:        {summary_path}")
    print(f"File manifest:  {files_path}  ({len(files_rows)} files)")
    return files_path


def cmd_download(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    rows = list(csv.DictReader(open(args.manifest)))
    if not rows:
        raise SystemExit("Empty manifest")
    cap = int(args.max_download_gb * (1024 ** 3))
    only_kinds = {x.lower() for x in (args.only or [])}
    only_formats = {x.lower() for x in (args.formats or [])}
    used = 0
    log_rows: "list[dict]" = []
    for r in rows:
        if r.get("status") and r["status"].lower() not in {"released",
                                                              "in progress"}:
            continue
        kind = (r.get("output_type") or "").lower()
        fmt = (r.get("file_format") or "").lower()
        if only_kinds and kind not in only_kinds and \
           not any(k in kind for k in only_kinds):
            continue
        if only_formats and fmt not in only_formats:
            continue
        size = int(r.get("file_size_bytes") or 0)
        if used + size > cap:
            logging.info("Skipping %s — cap reached", r.get("file_accession"))
            continue
        exp = r.get("experiment") or "misc"
        out_dir = DOWNLOAD_DIR / safe_label(exp)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / safe_label(
            f"{r.get('file_accession','file')}.{fmt or 'bin'}"
        )
        if out_file.exists() and out_file.stat().st_size == size:
            written = size
        else:
            try:
                written = download_file(r["url"], out_file,
                                          max_bytes=size + (1 << 20))
            except Exception as e:
                logging.warning("Failed %s: %s", r.get("url"), e)
                continue
        used += written
        log_rows.append({**r, "downloaded_path": str(out_file),
                          "downloaded_bytes": written})
    log_path = MANIFEST_DIR / f"{timestamp()}_download_log.csv"
    write_csv(log_path, log_rows)
    print(f"Downloaded {len(log_rows)} files / {used / (1024**3):.2f} GB")
    print(f"Log: {log_path}")
    return log_path


# --------------------------- describe --------------------------------------

def cmd_describe(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    acc = args.accession
    meta_path = METADATA_DIR / f"{acc}.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
    else:
        meta = hydrate_experiment(acc)
        meta_path.write_text(json.dumps(meta, indent=2, default=str))

    summary = summarize_experiment(meta)
    files = files_for_experiment(meta)
    fmts = Counter(r["file_format"] for r in files if r.get("file_format"))
    output_types = Counter(r["output_type"] for r in files
                              if r.get("output_type"))
    histone = summary.get("target", "")
    histone_note = HISTONE_MARK_NOTES.get(histone, "")
    out_dir = REPORT_DIR / f"{ts}_describe_{safe_label(acc)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / f"{acc}_describe.md"
    lines = [
        f"# ENCODE experiment: `{acc}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Summary",
        "",
        f"- **Assay**: {summary.get('assay_title','')}  "
            f"({summary.get('assay_term_name','')})",
        f"- **Biosample**: {summary.get('biosample','')}  "
            f"_({summary.get('biosample_classification','')})_",
        (f"- **Target**: `{histone}` — {histone_note}"
         if histone and histone_note
         else (f"- **Target**: `{histone}`" if histone else "")),
        f"- **Assembly**: {summary.get('assembly','')}",
        f"- **Lab**: {summary.get('lab','')}",
        f"- **Status**: `{summary.get('status','')}`",
        f"- **URL**: {summary.get('url','')}",
        "",
        f"### Description",
        "",
        f"{(meta.get('description') or '_no description_').strip()}",
        "",
        f"## File inventory ({len(files)})",
        "",
        "| file_format | count |",
        "|---|---:|",
    ]
    for fmt, n in sorted(fmts.items(), key=lambda x: -x[1]):
        lines.append(f"| {fmt} | {n} |")
    lines += ["", "| output_type | count |", "|---|---:|"]
    for ot, n in sorted(output_types.items(), key=lambda x: -x[1]):
        lines.append(f"| {ot} | {n} |")

    # Quality metrics, if present in meta.
    qc = meta.get("quality_metric_summary") or meta.get("audit", None)
    if isinstance(qc, list) and qc:
        lines += ["", "## Quality metrics (truncated)", ""]
        for q in qc[:10]:
            lines.append(f"- `{q}`"[:300])

    lines += [
        "",
        "## Suggested follow-ups",
        "",
        "- `igvfagent encode download --manifest <files.csv> "
        "--only signal --formats bigWig --max-download-gb 1`",
        "- `igvfagent encode analyze-peaks --bed <peaks.bed>`  "
        "(after downloading the narrowPeak / broadPeak files)",
        "- `igvfagent encode integrate-ccre --bed <peaks.bed>`  "
        "(overlay against SCREEN cCREs)",
        "- `igvfagent encode super-enhancers --bed <H3K27ac_peaks.bed>`  "
        "(only for H3K27ac / Mediator / BRD4 ChIP-seq)",
        "- `igvfagent encode browser --region chr19:44903000-44912000 "
        "--track <peaks.bed:peaks>`  (multi-track SVG view)",
    ]
    report.write_text("\n".join(l for l in lines if l is not None))
    files_path = out_dir / f"{acc}_files.csv"
    write_csv(files_path, files)
    print(f"Report:        {report}")
    print(f"Files manifest:{files_path}")
    return report


# --------------------------- BED parsing -----------------------------------

def parse_bed(path: Path, max_rows: Optional[int] = None
                ) -> "tuple[list[dict], dict]":
    """Read a BED / narrowPeak / broadPeak file (plain or .gz) into a list
    of dicts. Returns (rows, stats) where stats contains chromosome
    counts and median width."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: "list[dict]" = []
    widths: "list[int]" = []
    chrom_counts: Counter = Counter()
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            chrom = parts[0]
            row = {"chrom": chrom, "start": start, "end": end}
            if len(parts) > 3: row["name"] = parts[3]
            if len(parts) > 4:
                try:
                    row["score"] = float(parts[4])
                except ValueError:
                    row["score"] = 0.0
            if len(parts) > 5: row["strand"] = parts[5]
            # narrowPeak: 7=signalValue, 8=pValue, 9=qValue, 10=peak
            if len(parts) > 6:
                try: row["signal"] = float(parts[6])
                except ValueError: pass
            if len(parts) > 7:
                try: row["pvalue"] = float(parts[7])
                except ValueError: pass
            if len(parts) > 8:
                try: row["qvalue"] = float(parts[8])
                except ValueError: pass
            rows.append(row)
            widths.append(end - start)
            chrom_counts[chrom] += 1
            if max_rows and len(rows) >= max_rows:
                break
    median_w = float(sorted(widths)[len(widths) // 2]) if widths else 0.0
    stats = {
        "n_peaks":       len(rows),
        "n_chromosomes": len(chrom_counts),
        "median_width":  median_w,
        "min_width":     min(widths) if widths else 0,
        "max_width":     max(widths) if widths else 0,
        "chrom_counts":  dict(chrom_counts),
    }
    return rows, stats


# --------------------------- analyze-peaks ----------------------------------

def cmd_analyze_peaks(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    bed = Path(args.bed)
    rows, stats = parse_bed(bed)
    if not rows:
        raise SystemExit(f"No peaks parsed from {bed}")
    ts = timestamp()
    label = safe_label(args.label or bed.stem)
    out_dir = REPORT_DIR / f"{ts}_analyze_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / "Plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Plots if matplotlib is available
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        widths = np.array([r["end"] - r["start"] for r in rows])
        scores = np.array([r.get("score", 0.0) for r in rows])

        fig, axes = plt.subplots(1, 3, figsize=(13, 3.5),
                                   constrained_layout=True)
        axes[0].hist(np.log10(widths + 1), bins=40, color="#0072B2")
        axes[0].set_xlabel("log10(peak width)"); axes[0].set_ylabel("Count")
        axes[0].set_title("Peak width distribution")

        if (scores > 0).any():
            axes[1].hist(scores, bins=40, color="#D55E00")
            axes[1].set_xlabel("Score"); axes[1].set_ylabel("Count")
            axes[1].set_title("Peak score distribution")
        else:
            axes[1].axis("off"); axes[1].set_title("(no score column)")

        chrom_counts = Counter(r["chrom"] for r in rows)
        chroms_sorted = sorted(chrom_counts.items(), key=lambda x: -x[1])[:24]
        axes[2].bar(range(len(chroms_sorted)),
                     [c[1] for c in chroms_sorted], color="#009E73")
        axes[2].set_xticks(range(len(chroms_sorted)))
        axes[2].set_xticklabels([c[0] for c in chroms_sorted],
                                  rotation=70, fontsize=7)
        axes[2].set_ylabel("Peak count")
        axes[2].set_title("Peaks per chromosome (top 24)")
        for ext in ("png", "svg"):
            fig.savefig(plot_dir / f"{label}_peak_qc.{ext}", dpi=200,
                        bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_peak_qc.md"
    lines = [
        f"# ENCODE peak QC: `{bed.name}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Summary",
        "",
        f"- Peaks: **{stats['n_peaks']:,}**",
        f"- Chromosomes covered: **{stats['n_chromosomes']}**",
        f"- Width median: **{stats['median_width']:.0f}** bp  "
            f"(min {stats['min_width']}, max {stats['max_width']})",
        "",
        "## Peaks per chromosome (top)",
        "",
        "| chrom | peaks |",
        "|---|---:|",
    ]
    for chrom, n in sorted(stats["chrom_counts"].items(),
                              key=lambda x: -x[1])[:24]:
        lines.append(f"| {chrom} | {n} |")
    lines += ["", f"![QC plots](Plots/{label}_peak_qc.png)", ""]
    report.write_text("\n".join(lines))
    print(f"Report: {report}")
    return report


# --------------------------- super-enhancers --------------------------------

def stitch_peaks(rows: "list[dict]",
                   stitching_distance: int = 12500,
                   exclude_tss: "Optional[list[tuple[str,int]]]" = None,
                   tss_distance: int = 2000) -> "list[dict]":
    """Cluster peaks within ``stitching_distance`` bp on the same chrom.

    Returns a list of stitched regions with summed signal/score and
    constituent peak count. Peaks within ``tss_distance`` of any TSS in
    ``exclude_tss`` are removed before stitching (ROSE option).
    """
    if exclude_tss:
        tss_by_chr: "dict[str, list[int]]" = defaultdict(list)
        for chrom, pos in exclude_tss:
            tss_by_chr[chrom].append(pos)
        for chrom in tss_by_chr:
            tss_by_chr[chrom].sort()
        kept = []
        for r in rows:
            ts_list = tss_by_chr.get(r["chrom"], [])
            mid = (r["start"] + r["end"]) // 2
            if ts_list:
                # binary search for nearest TSS
                import bisect
                i = bisect.bisect_left(ts_list, mid)
                near = []
                if i < len(ts_list): near.append(ts_list[i])
                if i > 0:            near.append(ts_list[i - 1])
                if any(abs(t - mid) <= tss_distance for t in near):
                    continue
            kept.append(r)
        rows = kept

    rows_sorted = sorted(rows, key=lambda r: (r["chrom"], r["start"]))
    stitched: "list[dict]" = []
    cur: Optional[dict] = None
    for r in rows_sorted:
        if cur is None or r["chrom"] != cur["chrom"] \
                or r["start"] - cur["end"] > stitching_distance:
            if cur is not None:
                stitched.append(cur)
            cur = {"chrom":      r["chrom"],
                    "start":      r["start"],
                    "end":        r["end"],
                    "score_sum":  r.get("signal", r.get("score", 0.0)),
                    "n_peaks":    1}
            continue
        cur["end"] = max(cur["end"], r["end"])
        cur["score_sum"] += r.get("signal", r.get("score", 0.0))
        cur["n_peaks"] += 1
    if cur is not None:
        stitched.append(cur)
    return stitched


def find_inflection(scores: "list[float]") -> int:
    """ROSE inflection: project onto the [0,1] x [0,1] diagonal and pick
    the index where the curve is furthest above the diagonal."""
    if not scores:
        return 0
    sorted_scores = sorted(scores, reverse=True)
    n = len(sorted_scores)
    smax = sorted_scores[0] or 1.0
    # Normalize x from rank 0..1, y from score/max
    best_idx = 0
    best_dist = 0.0
    for i, s in enumerate(sorted_scores):
        x = i / max(n - 1, 1)
        y = s / smax
        # Distance ABOVE the diagonal (y - x)
        d = y - x
        # Want the LAST high point — break ties on x
        if d >= best_dist:
            best_dist = d; best_idx = i
    # The inflection is roughly where the slope = 1 from below; we
    # take the index from the high end of the curve where y - x first
    # falls below a threshold.
    # Simpler fallback: the index 1/3 of the way down — keep the
    # geometric estimate as the SE call.
    return best_idx + 1   # number of regions ABOVE the inflection


def cmd_super_enhancers(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    bed = Path(args.bed)
    rows, _ = parse_bed(bed)
    if not rows:
        raise SystemExit(f"No peaks in {bed}")

    excl: "Optional[list[tuple[str,int]]]" = None
    if args.tss_bed:
        tss_rows, _ = parse_bed(Path(args.tss_bed))
        excl = [(r["chrom"], r["start"]) for r in tss_rows]

    stitched = stitch_peaks(
        rows,
        stitching_distance=args.stitching_distance,
        exclude_tss=excl,
        tss_distance=args.tss_distance,
    )
    if not stitched:
        raise SystemExit("Stitching produced no regions.")

    stitched.sort(key=lambda r: r["score_sum"], reverse=True)
    n_se = find_inflection([r["score_sum"] for r in stitched])
    for i, r in enumerate(stitched):
        r["rank"] = i + 1
        r["is_super_enhancer"] = i < n_se

    ts = timestamp()
    label = safe_label(args.label or bed.stem + "_SE")
    out_dir = REPORT_DIR / f"{ts}_SE_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    se_path = out_dir / f"{label}_super_enhancers.bed"
    typical_path = out_dir / f"{label}_typical_enhancers.bed"
    all_path = out_dir / f"{label}_stitched_all.csv"
    with se_path.open("w") as f:
        for r in stitched[:n_se]:
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t"
                     f"SE_{r['rank']}\t{r['score_sum']:.3f}\t.\n")
    with typical_path.open("w") as f:
        for r in stitched[n_se:]:
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t"
                     f"TE_{r['rank']}\t{r['score_sum']:.3f}\t.\n")
    write_csv(all_path, stitched,
              cols=["rank", "is_super_enhancer", "chrom", "start", "end",
                     "score_sum", "n_peaks"])

    # Hockey-stick plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        ranks = np.arange(1, len(stitched) + 1)
        scores = np.array([r["score_sum"] for r in stitched])
        # rank from low to high so the "hockey stick" turns up at the right
        order = np.argsort(scores)
        plt.figure(figsize=(7, 5))
        plt.plot(np.arange(len(scores))[::-1], scores[order], lw=1.5,
                  color="#0072B2")
        plt.axvline(len(scores) - n_se - 0.5, color="#D55E00", lw=1.0,
                     linestyle="--",
                     label=f"inflection: {n_se} super-enhancers")
        plt.xlabel("Stitched-enhancer rank (low → high)")
        plt.ylabel("Summed signal (score)")
        plt.title(f"Super-enhancer ranking — {label}")
        plt.legend()
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_hockey_stick.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_SE_report.md"
    lines = [
        f"# Super-enhancer call: `{bed.name}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Parameters",
        "",
        f"- Stitching distance: **{args.stitching_distance:,} bp**",
        (f"- TSS exclusion: peaks within **{args.tss_distance} bp** of any "
         f"TSS in `{args.tss_bed}` are removed before stitching"
         if args.tss_bed else "- TSS exclusion: _none_"),
        "",
        "## Result",
        "",
        f"- Stitched regions: **{len(stitched):,}**",
        f"- Super-enhancers: **{n_se:,}** (above the inflection point)",
        f"- Typical enhancers: **{len(stitched) - n_se:,}**",
        "",
        "## Top 25 super-enhancers",
        "",
        "| rank | chrom | start | end | score_sum | n_peaks |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for r in stitched[:25]:
        lines.append(f"| {r['rank']} | {r['chrom']} | {r['start']} | "
                      f"{r['end']} | {r['score_sum']:.2f} | {r['n_peaks']} |")
    lines += [
        "",
        "## Outputs",
        "",
        f"- Super-enhancer BED: `{se_path.relative_to(ROOT)}`",
        f"- Typical-enhancer BED: `{typical_path.relative_to(ROOT)}`",
        f"- Full stitched CSV: `{all_path.relative_to(ROOT)}`",
        f"- Hockey-stick plot: `Plots/{label}_hockey_stick.png` "
            f"+ SVG variant",
        "",
        "## Notes",
        "",
        "ROSE-style stitching merges nearby peaks (default 12.5 kb) and "
        "ranks the resulting regions by total signal. The 'inflection' "
        "is the geometric break in the ranked score curve — regions above "
        "it are typically called super-enhancers in the literature.",
        "If you ran this on an H3K4me3-only or polycomb mark, the call "
        "is *not* meaningful: super-enhancers are most informative for "
        "H3K27ac, BRD4, MED1, MED12, P300, and similar enhancer marks.",
    ]
    report.write_text("\n".join(lines))
    print(f"Report:           {report}")
    print(f"Super-enhancers:  {se_path}  ({n_se} regions)")
    return report


# --------------------------- integrate-ccre ---------------------------------

def _ensure_screen_ccres(cclass: Optional[str] = None) -> Path:
    """Download the SCREEN cCRE BED if not cached and return its path."""
    name = "GRCh38-cCREs.bed" if not cclass else f"GRCh38-cCREs.{cclass}.bed"
    cache = CACHE_DIR / name
    if cache.exists() and cache.stat().st_size > 0:
        return cache
    url = f"{WENGLAB_DL}/Registry-V4/{name}"
    download_file(url, cache)
    return cache


def overlap_peaks_with_ccres(peaks: "list[dict]",
                                ccre_bed: Path) -> "tuple[list[dict], dict]":
    """Naive O(n+m) chrom-bucketed overlap. Each peak is assigned the
    cCRE class of its first overlapping cCRE; peaks with no overlap get
    class 'none'."""
    ccre_by_chr: "dict[str, list[tuple[int, int, str]]]" = defaultdict(list)
    with open(ccre_bed) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            cclass = parts[5] if len(parts) > 5 else (
                parts[3].split(",", 1)[0] if "," in parts[3] else parts[3]
            )
            ccre_by_chr[parts[0]].append((start, end, cclass))
    for chrom in ccre_by_chr:
        ccre_by_chr[chrom].sort()

    annotated: "list[dict]" = []
    class_counts: Counter = Counter()
    for p in peaks:
        bucket = ccre_by_chr.get(p["chrom"], [])
        # binary search by start
        import bisect
        starts = [x[0] for x in bucket]
        i = bisect.bisect_left(starts, p["end"])
        cclass = "none"
        # Walk back at most a few entries for overlap
        for j in range(max(0, i - 200), min(len(bucket), i + 1)):
            cs, ce, cls_ = bucket[j]
            if ce >= p["start"] and cs <= p["end"]:
                cclass = cls_; break
        class_counts[cclass] += 1
        annotated.append({**p, "ccre_class": cclass})
    return annotated, dict(class_counts)


def cmd_integrate_ccre(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    bed = Path(args.bed)
    peaks, _ = parse_bed(bed)
    if not peaks:
        raise SystemExit(f"No peaks in {bed}")
    ccre_bed = (Path(args.ccre_bed) if args.ccre_bed
                else _ensure_screen_ccres(args.ccre_class))

    annotated, counts = overlap_peaks_with_ccres(peaks, ccre_bed)
    ts = timestamp()
    label = safe_label(args.label or bed.stem + "_ccre")
    out_dir = REPORT_DIR / f"{ts}_ccre_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{label}_peaks_with_ccre.csv"
    write_csv(csv_path, annotated,
              cols=["chrom", "start", "end", "name", "score", "signal",
                     "pvalue", "qvalue", "ccre_class"])

    # Stacked-bar plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        order = list(CCRE_CLASSES) + ["none"]
        labels = [CCRE_CLASSES.get(o, {"label": o})["label"] for o in order]
        values = [counts.get(o, 0) for o in order]
        colors = [CCRE_CLASSES.get(o, {"color": "#aaaaaa"})["color"]
                   for o in order]
        plt.figure(figsize=(8, 4))
        plt.barh(range(len(order)), values, color=colors)
        plt.yticks(range(len(order)), labels)
        plt.xlabel("Peaks")
        plt.title(f"cCRE class composition — {label}")
        plt.tight_layout()
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_ccre_classes.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_ccre_report.md"
    total = len(annotated)
    overlap_n = total - counts.get("none", 0)
    lines = [
        f"# Peaks ↔ SCREEN cCRE overlap: `{bed.name}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"- cCRE BED used: `{ccre_bed.name}`",
        f"- Peaks: **{total:,}**",
        f"- Peaks overlapping any cCRE: **{overlap_n:,}** "
            f"({overlap_n / total * 100:.1f}%)",
        "",
        "## Class distribution",
        "",
        "| Class | Peaks | % |",
        "|---|---:|---:|",
    ]
    for cls in (list(CCRE_CLASSES) + ["none"]):
        n = counts.get(cls, 0)
        pct = n / total * 100 if total else 0.0
        label_ = CCRE_CLASSES.get(cls, {"label": cls})["label"]
        lines.append(f"| {label_} | {n:,} | {pct:.2f} |")
    lines += [
        "",
        f"![cCRE composition](Plots/{label}_ccre_classes.png)",
        "",
        f"Annotated peak CSV: `{csv_path.relative_to(ROOT)}`",
    ]
    report.write_text("\n".join(lines))
    print(f"Report:    {report}")
    print(f"Annotated: {csv_path}")
    return report


# --------------------------- browser SVG ------------------------------------

def _parse_track_arg(spec: str) -> "tuple[str, Path]":
    """``label:path`` -> (label, Path). If no ':', label = filename."""
    if ":" in spec:
        label, _, path = spec.partition(":")
    else:
        label, path = Path(spec).stem, spec
    return label, Path(path)


def _bed_in_region(bed: Path, chrom: str, start: int, end: int
                    ) -> "list[tuple[int, int, str, float]]":
    rows, _ = parse_bed(bed)
    out = []
    for r in rows:
        if r["chrom"] != chrom and r["chrom"] != f"chr{chrom.replace('chr','')}":
            continue
        if r["end"] < start or r["start"] > end:
            continue
        out.append((max(r["start"], start), min(r["end"], end),
                     r.get("name", ""), r.get("score", 0.0)))
    return out


def cmd_browser(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    m = re.match(r"^(chr\w+):(\d+)-(\d+)$", args.region)
    if not m:
        raise SystemExit("--region must look like chr19:44903000-44912000")
    chrom = m.group(1); start = int(m.group(2)); end = int(m.group(3))
    span = max(end - start, 1)

    tracks: "list[dict]" = []
    for spec in args.track or []:
        label, path = _parse_track_arg(spec)
        if not path.exists():
            raise SystemExit(f"Track file not found: {path}")
        tracks.append({"label": label, "path": path,
                        "rows": _bed_in_region(path, chrom, start, end)})
    if args.with_ccre:
        ccre_bed = (Path(args.ccre_bed) if args.ccre_bed
                    else _ensure_screen_ccres())
        tracks.append({"label": "SCREEN cCREs", "path": ccre_bed,
                        "rows": _bed_in_region(ccre_bed, chrom, start, end)})

    width = args.width
    track_h = 28
    margin_top = 50
    margin_left = 90
    margin_right = 30
    plot_w = width - margin_left - margin_right
    height = margin_top + track_h * (len(tracks) + 1) + 40

    def x_of(pos: int) -> float:
        return margin_left + (pos - start) / span * plot_w

    parts: "list[str]" = [
        f'<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        f'font-family="Helvetica, Arial, sans-serif" font-size="11">',
        f'<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_left}" y="22" font-size="13" font-weight="bold">'
        f'{args.region}  ·  {span:,} bp  ·  {len(tracks)} track(s)</text>',
    ]

    # Coordinate axis
    n_ticks = 6
    parts.append(f'<line x1="{margin_left}" y1="{margin_top - 8}" '
                  f'x2="{width - margin_right}" y2="{margin_top - 8}" '
                  f'stroke="black" stroke-width="0.5"/>')
    for i in range(n_ticks + 1):
        pos = start + (i / n_ticks) * span
        x = x_of(int(pos))
        parts.append(f'<line x1="{x}" y1="{margin_top - 12}" '
                      f'x2="{x}" y2="{margin_top - 4}" stroke="black" '
                      f'stroke-width="0.5"/>')
        parts.append(f'<text x="{x}" y="{margin_top - 16}" '
                      f'text-anchor="middle" font-size="9">'
                      f'{int(pos):,}</text>')

    # Track strips
    for ti, tr in enumerate(tracks):
        y_top = margin_top + ti * track_h
        y_mid = y_top + track_h / 2
        # Label
        parts.append(f'<text x="{margin_left - 6}" y="{y_mid + 4}" '
                      f'text-anchor="end" font-size="10" font-weight="bold">'
                      f'{tr["label"]}</text>')
        # Track baseline
        parts.append(f'<line x1="{margin_left}" y1="{y_mid}" '
                      f'x2="{width - margin_right}" y2="{y_mid}" '
                      f'stroke="#dddddd" stroke-width="0.5"/>')
        for s, e, name, score in tr["rows"]:
            x1 = x_of(s); x2 = x_of(e)
            w = max(x2 - x1, 1.0)
            color = _color_for_track(tr["label"], name)
            parts.append(f'<rect x="{x1:.2f}" y="{y_top + 6}" '
                          f'width="{w:.2f}" height="{track_h - 12}" '
                          f'fill="{color}" fill-opacity="0.85" '
                          f'stroke="{color}" stroke-width="0.4"/>')
        if not tr["rows"]:
            parts.append(f'<text x="{(margin_left + width - margin_right)/2}" '
                          f'y="{y_mid + 4}" text-anchor="middle" '
                          f'font-size="9" fill="#888">'
                          f'(no features in region)</text>')

    parts.append('</svg>')
    ts = timestamp()
    label = safe_label(args.label or f"{chrom}_{start}_{end}")
    out = PLOT_DIR / f"{ts}_browser_{label}.svg"
    out.write_text("\n".join(parts))
    print(f"Browser SVG: {out}")
    return out


def _color_for_track(track_label: str, name: str) -> str:
    if "ccre" in track_label.lower() or "SCREEN" in track_label:
        for cls, info in CCRE_CLASSES.items():
            if cls.lower() in (name or "").lower():
                return info["color"]
        return "#888888"
    return "#0072B2"


# --------------------------- write-playbook ---------------------------------

def cmd_write_playbook(_args) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "ENCODE_PIPELINE_SKILLS.md"
    lines = [
        "# Skill: ENCODE bulk-genomics pipeline",
        "",
        "End-to-end retrieval, description, peak QC, super-enhancer "
        "calling, cCRE integration, and browser-style SVG visualization "
        "for the major ENCODE bulk assays:",
        "",
        ", ".join(f"`{a}`" for a in ENCODE_ASSAYS),
        "",
        "## Subcommands",
        "",
        "### `retrieve` — find experiments by assay / biosample / target",
        "",
        "```bash",
        "igvfagent encode retrieve --assay 'Histone ChIP-seq' \\",
        "    --target H3K27ac --biosample K562 --assembly GRCh38 \\",
        "    --limit 50 --label k562_h3k27ac",
        "igvfagent encode retrieve --assay ATAC-seq --biosample 'liver'",
        "igvfagent encode retrieve --assay Hi-C --biosample GM12878",
        "```",
        "",
        "### `manifest` — per-file inventory for one or more accessions",
        "",
        "```bash",
        "igvfagent encode manifest --accessions ENCSR000DUB --label demo",
        "```",
        "",
        "### `download` — pull files under a size cap",
        "",
        "```bash",
        "igvfagent encode download \\",
        "    --manifest Data/Manifests/ENCODE/<files.csv> \\",
        "    --max-download-gb 5 --formats bed bigWig",
        "```",
        "",
        "### `describe` — plain-language report for one experiment",
        "",
        "```bash",
        "igvfagent encode describe --accession ENCSR000DUB",
        "```",
        "",
        "### `analyze-peaks` — peak QC stats from a BED file",
        "",
        "```bash",
        "igvfagent encode analyze-peaks --bed peaks.bed.gz",
        "```",
        "",
        "Outputs: peak count, width distribution, score distribution, "
        "per-chromosome counts, plus a 3-panel QC PNG/SVG.",
        "",
        "### `super-enhancers` — ROSE-style super-enhancer call",
        "",
        "```bash",
        "igvfagent encode super-enhancers --bed h3k27ac_peaks.bed \\",
        "    --stitching-distance 12500 --tss-bed tss.bed --tss-distance 2000 \\",
        "    --label k562_h3k27ac",
        "```",
        "",
        "Stitches enhancer-mark peaks within `--stitching-distance` (default "
        "12.5 kb), optionally excludes TSS-proximal peaks given a TSS BED, "
        "ranks the stitched regions by summed signal, and finds the "
        "geometric inflection point. Emits two BED files (super- and typical-"
        "enhancers) plus a hockey-stick PNG/SVG.",
        "",
        "Reliable for H3K27ac, BRD4, MED1, MED12, P300; not meaningful for "
        "polycomb / heterochromatin marks.",
        "",
        "### `integrate-ccre` — overlay peaks with SCREEN cCREs",
        "",
        "```bash",
        "igvfagent encode integrate-ccre --bed peaks.bed",
        "igvfagent encode integrate-ccre --bed peaks.bed \\",
        "    --ccre-bed Data/Cache/cCRELinkage/<custom.bed>",
        "```",
        "",
        "Auto-downloads the SCREEN V4 cCRE BED on first use (cached in "
        "`Data/Cache/ENCODE/`). Annotates each peak with its cCRE class "
        "(PLS / pELS / dELS / CTCF-only / DNase-H3K4me3) and emits a "
        "stacked-bar overview.",
        "",
        "### `browser` — IGV-style multi-track SVG",
        "",
        "```bash",
        "igvfagent encode browser \\",
        "    --region chr19:44903000-44912000 \\",
        "    --track 'H3K27ac peaks:peaks.bed' \\",
        "    --track 'ATAC peaks:atac_peaks.bed' \\",
        "    --with-ccre \\",
        "    --label apoe_locus",
        "```",
        "",
        "Each `--track LABEL:PATH` adds a horizontal track strip; "
        "`--with-ccre` overlays SCREEN cCREs colored by class.",
        "",
        "## How this chains with other skills",
        "",
        "- After `retrieve` / `manifest` / `download`, hand peak BEDs to "
        "`analyze-peaks` and `integrate-ccre`. For H3K27ac (or any "
        "enhancer mark), pipe to `super-enhancers` next.",
        "- `browser` accepts any combination of BED tracks — including "
        "outputs from `splitseq_pipeline` / `multiome_10x_pipeline` "
        "ATAC peak BEDs, advanced-variant-analysis cCRE overlaps, etc.",
        "- The agent runtime exposes `encode_retrieve`, `encode_describe`, "
        "`encode_super_enhancers`, `encode_integrate_ccre`, and "
        "`encode_browser` as tools so a single `igvfagent ask` can drive "
        "the whole pipeline end-to-end.",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI -------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="ENCODE bulk-genomics retrieval, analysis, and "
                    "visualization pipeline."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("retrieve",
                        help="Search ENCODE for experiments.")
    s.add_argument("--assay", required=True,
                    help="One of: " + ", ".join(ENCODE_ASSAYS.keys()))
    s.add_argument("--biosample", default=None)
    s.add_argument("--target", default=None,
                    help="ChIP-seq / Histone target (e.g. H3K27ac, CTCF).")
    s.add_argument("--assembly", default=None)
    s.add_argument("--status", default="released")
    s.add_argument("--limit", type=int, default=50)
    s.add_argument("--fetch-file-details", action="store_true")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_retrieve)

    s = sub.add_parser("manifest",
                        help="Per-file inventory for one or more accessions.")
    s.add_argument("--accessions", required=True,
                    help="Comma-separated experiment accessions.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_manifest)

    s = sub.add_parser("download",
                        help="Download files referenced by a manifest.")
    s.add_argument("--manifest", required=True)
    s.add_argument("--max-download-gb", type=float, default=2.0)
    s.add_argument("--only", nargs="*", default=None,
                    help="Filter by output_type substrings.")
    s.add_argument("--formats", nargs="*", default=None,
                    help="Filter by file_format (e.g. 'bigWig' 'bed').")
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("describe",
                        help="Plain-language report for one experiment.")
    s.add_argument("--accession", required=True)
    s.set_defaults(func=cmd_describe)

    s = sub.add_parser("analyze-peaks",
                        help="Peak QC from a BED file.")
    s.add_argument("--bed", required=True)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_analyze_peaks)

    s = sub.add_parser("super-enhancers",
                        help="ROSE-style super-enhancer call from a BED.")
    s.add_argument("--bed", required=True)
    s.add_argument("--stitching-distance", type=int, default=12500)
    s.add_argument("--tss-bed", default=None)
    s.add_argument("--tss-distance", type=int, default=2000)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_super_enhancers)

    s = sub.add_parser("integrate-ccre",
                        help="Overlay BED peaks with SCREEN cCREs.")
    s.add_argument("--bed", required=True)
    s.add_argument("--ccre-bed", default=None,
                    help="Optional local cCRE BED. Auto-downloads SCREEN V4.")
    s.add_argument("--ccre-class", default=None,
                    help="Restrict download to PLS / pELS / dELS / CTCF.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_integrate_ccre)

    s = sub.add_parser("browser",
                        help="Render IGV-style multi-track SVG.")
    s.add_argument("--region", required=True,
                    help="chr19:44903000-44912000")
    s.add_argument("--track", action="append", default=None,
                    help="LABEL:PATH for each track. Repeatable.")
    s.add_argument("--with-ccre", action="store_true",
                    help="Add a SCREEN cCRE track (auto-downloads).")
    s.add_argument("--ccre-bed", default=None)
    s.add_argument("--width", type=int, default=1000)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_browser)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/ENCODE_PIPELINE_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
