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
    is_gz = str(path).endswith(".gz")
    if not is_gz:
        try:
            with open(path, "rb") as probe:
                is_gz = probe.read(2) == b"\x1f\x8b"
        except Exception:
            pass
    opener = gzip.open if is_gz else open
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


# --------------------------- bigWig: FRiP + TSS heatmap --------------------

def _open_bigwig(path: Path):
    try:
        import pyBigWig
    except ImportError as e:
        raise RuntimeError(
            "bigWig analysis needs the `pyBigWig` package. "
            "Install with: pip install 'igvfagent[hic]'"
        ) from e
    return pyBigWig.open(str(path))


def cmd_bigwig_frip(args: argparse.Namespace) -> Path:
    """Compute FRiP (fraction of reads / signal in peaks) from a bigWig
    signal file + a peak BED. Reports per-chromosome FRiP plus a global
    FRiP and a stacked histogram of in-peak vs out-of-peak signal."""
    setup_logging(); mkdirs()
    bw = _open_bigwig(Path(args.bigwig))
    peaks, _ = parse_bed(Path(args.bed))
    if not peaks:
        raise SystemExit(f"No peaks parsed from {args.bed}")

    # Total signal per chromosome (use stats(... type='sum') over each chrom).
    total_per_chrom: "dict[str, float]" = {}
    for chrom, length in bw.chroms().items():
        s = bw.stats(chrom, 0, length, type="sum")
        total_per_chrom[chrom] = float(s[0]) if s and s[0] is not None else 0.0
    in_peak_per_chrom: "dict[str, float]" = defaultdict(float)
    n_per_chrom: "dict[str, int]" = defaultdict(int)
    for p in peaks:
        s = bw.stats(p["chrom"], p["start"], p["end"], type="sum")
        if s and s[0] is not None:
            in_peak_per_chrom[p["chrom"]] += float(s[0])
        n_per_chrom[p["chrom"]] += 1
    bw.close()

    rows: "list[dict]" = []
    g_in = g_total = 0.0
    for chrom, total in total_per_chrom.items():
        if total <= 0:
            continue
        in_p = in_peak_per_chrom.get(chrom, 0.0)
        rows.append({
            "chrom":         chrom,
            "n_peaks":       n_per_chrom.get(chrom, 0),
            "total_signal":  total,
            "in_peak_signal": in_p,
            "frip":          (in_p / total) if total else 0.0,
        })
        g_in += in_p; g_total += total
    rows.sort(key=lambda r: -r["total_signal"])
    global_frip = g_in / g_total if g_total else 0.0

    ts = timestamp()
    label = safe_label(args.label or Path(args.bigwig).stem + "_FRiP")
    out_dir = REPORT_DIR / f"{ts}_FRiP_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / f"{label}_per_chromosome.csv", rows,
              cols=["chrom", "n_peaks", "total_signal",
                     "in_peak_signal", "frip"])

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        top = rows[:24]
        x = np.arange(len(top))
        in_p = np.array([r["in_peak_signal"] for r in top])
        out_p = np.array([r["total_signal"] - r["in_peak_signal"]
                            for r in top])
        plt.figure(figsize=(10, 4))
        plt.bar(x, out_p, color="#cccccc", label="out of peaks")
        plt.bar(x, in_p, bottom=out_p, color="#D55E00",
                  label=f"in peaks (FRiP={global_frip:.3f})")
        plt.xticks(x, [r["chrom"] for r in top], rotation=70, fontsize=7)
        plt.ylabel("Signal sum")
        plt.title(f"FRiP — {label}  (global = {global_frip:.3f})")
        plt.legend(loc="upper right")
        plt.tight_layout()
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_frip.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_frip_report.md"
    report.write_text(
        f"# FRiP report — `{Path(args.bigwig).name}`\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Peaks: **{len(peaks):,}**\n"
        f"- Global FRiP: **{global_frip:.4f}**\n"
        f"  (in-peak signal / total signal across all chromosomes)\n\n"
        f"FRiP > 0.01 is the ENCODE minimum; > 0.05 is good. Histone "
        f"ChIP-seq tends to score lower than TF ChIP-seq because the "
        f"signal is broader.\n\n"
        f"![FRiP per chromosome](Plots/{label}_frip.png)\n"
    )
    print(f"Global FRiP: {global_frip:.4f}")
    print(f"Report:      {report}")
    return report


def cmd_bigwig_tss_heatmap(args: argparse.Namespace) -> Path:
    """Aggregate bigWig signal over a window centered on each anchor in
    the BED (TSS BED, peak summits, etc.) and render a sorted heatmap +
    average meta-profile."""
    setup_logging(); mkdirs()
    bw = _open_bigwig(Path(args.bigwig))
    anchors, _ = parse_bed(Path(args.anchor_bed))
    if not anchors:
        raise SystemExit(f"No anchors in {args.anchor_bed}")

    # Down-sample to keep the heatmap manageable
    if args.max_anchors and len(anchors) > args.max_anchors:
        step = max(1, len(anchors) // args.max_anchors)
        anchors = anchors[::step][: args.max_anchors]

    half = args.window // 2
    bins = args.bins
    bin_w = max(args.window // bins, 1)
    try:
        import numpy as np
    except ImportError:
        raise RuntimeError("Heatmap needs numpy. pip install 'igvfagent[analysis]'")

    matrix = np.zeros((len(anchors), bins), dtype=float)
    for i, a in enumerate(anchors):
        center = (a["start"] + a["end"]) // 2
        s = max(center - half, 0); e = center + half
        try:
            row = bw.stats(a["chrom"], s, e, type="mean", nBins=bins)
        except Exception:
            row = [0.0] * bins
        if row:
            matrix[i] = [float(r) if r is not None else 0.0 for r in row]
    bw.close()

    # Sort rows by total signal (descending) for the classic heatmap look.
    order = np.argsort(-matrix.sum(axis=1))
    matrix = matrix[order]

    ts = timestamp()
    label = safe_label(args.label or Path(args.bigwig).stem + "_TSS")
    out_dir = REPORT_DIR / f"{ts}_signal_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{label}_matrix.npy", matrix)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        vmax = float(np.nanpercentile(matrix, 99))
        fig = plt.figure(figsize=(7, 9), constrained_layout=True)
        gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 4])
        ax_top = fig.add_subplot(gs[0])
        ax_top.plot(np.linspace(-half, half, bins),
                     np.nanmean(matrix, axis=0), color="#0072B2", lw=1.5)
        ax_top.axvline(0, color="#888", lw=0.5)
        ax_top.set_ylabel("Mean signal")
        ax_top.set_title(f"{label} — meta-profile")
        ax_bot = fig.add_subplot(gs[1])
        ax_bot.imshow(matrix, aspect="auto", cmap="OrRd", vmin=0,
                       vmax=vmax,
                       extent=[-half, half, len(matrix), 0])
        ax_bot.set_xlabel("Distance from anchor center (bp)")
        ax_bot.set_ylabel(f"{len(matrix):,} anchors (sorted by signal)")
        for ext in ("png", "svg"):
            fig.savefig(plot_dir / f"{label}_heatmap.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_signal_report.md"
    report.write_text(
        f"# Anchor-centered signal heatmap — `{Path(args.bigwig).name}`\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Anchors: **{len(matrix):,}** (window ± {half} bp, {bins} bins)\n"
        f"- Source bigWig: `{args.bigwig}`\n"
        f"- Anchor BED: `{args.anchor_bed}`\n\n"
        f"![meta-profile + heatmap](Plots/{label}_heatmap.png)\n"
    )
    print(f"Heatmap report: {report}")
    return report


# --------------------------- Hi-C matrix + insulation ----------------------

def _load_hic_region(path: Path, region: str, resolution: int,
                       balance: bool = True):
    """Load a square contact matrix for ``region`` from .mcool or .hic."""
    chrom, start, end = _parse_region_str(region)
    suffix = path.suffix.lower()
    try:
        if suffix in (".mcool", ".cool"):
            import cooler
            uri = (f"{path}::resolutions/{resolution}"
                   if suffix == ".mcool" else str(path))
            c = cooler.Cooler(uri)
            mat = c.matrix(balance=balance).fetch(region)
            return mat, chrom, start, end
        if suffix == ".hic":
            import hicstraw  # type: ignore
            hic = hicstraw.HiCFile(str(path))
            norm = "KR" if balance else "NONE"
            mzd = hic.getMatrixZoomData(
                chrom.replace("chr", ""), chrom.replace("chr", ""),
                "observed", norm, "BP", resolution,
            )
            mat = mzd.getRecordsAsMatrix(start, end, start, end)
            return mat, chrom, start, end
        raise SystemExit(f"Unsupported Hi-C format: {suffix}")
    except ImportError as e:
        raise RuntimeError(
            "Hi-C analysis needs the `cooler` (.mcool) or `hic-straw` "
            "(.hic) package. Install with: pip install 'igvfagent[hic]'"
        ) from e


def _parse_region_str(region: str) -> "tuple[str, int, int]":
    m = re.match(r"^(chr\w+):(\d+)-(\d+)$", region)
    if not m:
        raise SystemExit("Region format: chr19:44900000-45100000")
    return m.group(1), int(m.group(2)), int(m.group(3))


def cmd_hic_matrix(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    mat, chrom, start, end = _load_hic_region(
        Path(args.input), args.region, args.resolution, balance=args.balance,
    )
    import numpy as np
    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = 0.0

    ts = timestamp()
    label = safe_label(args.label or
                         f"{Path(args.input).stem}_{chrom}_{start}_{end}_"
                         f"{args.resolution}")
    out_dir = REPORT_DIR / f"{ts}_hic_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"{label}_matrix.npy", mat)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        vmax = float(np.nanpercentile(mat, 99))
        plt.figure(figsize=(7, 7))
        plt.imshow(np.log1p(mat), cmap="YlOrRd", vmin=0,
                     vmax=np.log1p(vmax), origin="upper",
                     extent=[start, end, end, start])
        plt.title(f"Hi-C contacts — {chrom}:{start:,}-{end:,}  "
                    f"({args.resolution} bp bins)")
        plt.xlabel("Position"); plt.ylabel("Position")
        plt.colorbar(label="log1p(contacts)")
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_heatmap.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_hic_report.md"
    report.write_text(
        f"# Hi-C contact matrix — `{Path(args.input).name}`\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Region: `{chrom}:{start:,}-{end:,}`\n"
        f"- Resolution: **{args.resolution:,}** bp\n"
        f"- Balanced: `{args.balance}`\n"
        f"- Matrix shape: `{mat.shape}`\n"
        f"- Total contacts (sum, post-balance): "
        f"**{float(np.nansum(mat)):.2f}**\n\n"
        f"![contact heatmap](Plots/{label}_heatmap.png)\n"
        f"\nRaw matrix saved as `{label}_matrix.npy` "
        f"for downstream analysis (e.g. cooler-balance / "
        f"insulation-score chaining).\n"
    )
    print(f"Hi-C report: {report}")
    return report


def cmd_hic_insulation(args: argparse.Namespace) -> Path:
    """Crane-et-al-style insulation score across a region.

    For each diagonal bin i, score = mean(M[i-w:i, i:i+w]) where w is the
    sliding window in bins. Local minima of the insulation signal mark
    candidate TAD boundaries.
    """
    setup_logging(); mkdirs()
    mat, chrom, start, end = _load_hic_region(
        Path(args.input), args.region, args.resolution, balance=args.balance,
    )
    import numpy as np
    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = 0.0
    n = mat.shape[0]
    w_bins = max(args.window // args.resolution, 2)
    if n < 2 * w_bins + 4:
        raise SystemExit("Region too small for the requested window. "
                         "Use a wider --region or smaller --window.")
    ins = np.full(n, np.nan)
    for i in range(w_bins, n - w_bins):
        sub = mat[i - w_bins: i, i + 1: i + w_bins + 1]
        ins[i] = float(np.nanmean(sub)) if sub.size else np.nan
    # log-mean-normalize so the score is comparable across regions
    valid = ins[np.isfinite(ins) & (ins > 0)]
    if len(valid):
        log_mean = float(np.log2(valid.mean()))
        ins_norm = np.log2(np.clip(ins, 1e-9, None)) - log_mean
    else:
        ins_norm = ins

    # Local minima as candidate TAD boundaries
    boundaries: "list[int]" = []
    for i in range(1, n - 1):
        if (np.isfinite(ins_norm[i]) and
                ins_norm[i] < ins_norm[i - 1] and
                ins_norm[i] < ins_norm[i + 1] and
                ins_norm[i] < args.boundary_threshold):
            boundaries.append(i)

    ts = timestamp()
    label = safe_label(args.label or f"insulation_{chrom}_{start}_{end}")
    out_dir = REPORT_DIR / f"{ts}_ins_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Write boundaries BED
    bed_path = out_dir / f"{label}_boundaries.bed"
    with bed_path.open("w") as f:
        for b in boundaries:
            pos = start + b * args.resolution
            f.write(f"{chrom}\t{pos}\t{pos + args.resolution}\t"
                     f"boundary_{b}\t{ins_norm[b]:.3f}\t.\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = np.arange(n) * args.resolution + start
        plt.figure(figsize=(11, 3.5))
        plt.plot(x, ins_norm, color="#0072B2", lw=1.0)
        plt.axhline(args.boundary_threshold, color="#D55E00", ls="--",
                      lw=0.8)
        for b in boundaries:
            plt.axvline(start + b * args.resolution, color="#D55E00",
                          alpha=0.4, lw=0.4)
        plt.xlabel("Position"); plt.ylabel("log2 insulation (centered)")
        plt.title(f"Insulation score — {chrom}:{start:,}-{end:,}  "
                    f"(window {args.window:,} bp)")
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_insulation.{ext}", dpi=200,
                          bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_insulation_report.md"
    report.write_text(
        f"# Hi-C insulation score — `{Path(args.input).name}`\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Region: `{chrom}:{start:,}-{end:,}`\n"
        f"- Resolution: **{args.resolution:,}** bp; window "
        f"**{args.window:,}** bp\n"
        f"- Candidate TAD boundaries (below threshold "
        f"`{args.boundary_threshold}`): **{len(boundaries)}**\n"
        f"- Boundaries BED: `{bed_path.relative_to(ROOT)}`\n\n"
        f"![insulation score](Plots/{label}_insulation.png)\n\n"
        f"Boundaries are local minima of the centered log2 insulation "
        f"score below the configured threshold. They approximate TAD "
        f"boundaries (Crane et al. 2015). For loop-level calls "
        f"(HiCCUPS / Mustache / Peakachu), use a dedicated tool and "
        f"feed the resulting bedpe back through `loops-analyze`.\n"
    )
    print(f"Insulation report: {report}")
    print(f"Boundaries BED:    {bed_path}")
    return report


# --------------------------- loops-analyze (Hi-C + ChIA-PET) ---------------

def _parse_bedpe(path: Path) -> "list[dict]":
    rows: "list[dict]" = []
    # Detect gzip by magic bytes as well as by extension — ENCODE files
    # are sometimes downloaded with the extension stripped.
    is_gz = str(path).endswith(".gz")
    if not is_gz:
        try:
            with open(path, "rb") as probe:
                is_gz = probe.read(2) == b"\x1f\x8b"
        except Exception:
            pass
    opener = gzip.open if is_gz else open
    with opener(path, "rt") as f:
        for line in f:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                s1, e1, s2, e2 = int(p[1]), int(p[2]), int(p[4]), int(p[5])
            except ValueError:
                continue
            rows.append({
                "chrom1": p[0], "start1": s1, "end1": e1,
                "chrom2": p[3], "start2": s2, "end2": e2,
                "name":   p[6] if len(p) > 6 else "",
                "score":  float(p[7]) if (len(p) > 7
                            and re.match(r"^-?\d+(\.\d+)?$", p[7] or "")) else 0.0,
            })
    return rows


def cmd_loops_analyze(args: argparse.Namespace) -> Path:
    """Analyze a .bedpe loop / interaction file (Hi-C, ChIA-PET, capture
    Hi-C). Reports per-anchor count, intra- vs inter-chromosomal
    breakdown, loop length distribution, and (if peak BEDs are
    provided) anchor ↔ peak overlap stats."""
    setup_logging(); mkdirs()
    loops = _parse_bedpe(Path(args.bedpe))
    if not loops:
        raise SystemExit(f"No loops parsed from {args.bedpe}")

    intra = [L for L in loops if L["chrom1"] == L["chrom2"]]
    inter = [L for L in loops if L["chrom1"] != L["chrom2"]]
    lengths = [(((L["start2"] + L["end2"]) // 2) -
                ((L["start1"] + L["end1"]) // 2)) for L in intra]
    lengths = [abs(x) for x in lengths]

    # Peak overlap (anchor sets)
    peak_overlap_stats: "dict[str, dict]" = {}
    for spec in (args.peaks or []):
        label, _, path = spec.partition(":")
        if not path:
            label, path = Path(spec).stem, spec
        peaks, _ = parse_bed(Path(path))
        anchor1_overlap = 0
        anchor2_overlap = 0
        both = 0
        peak_by_chr: "dict[str, list[tuple[int,int]]]" = defaultdict(list)
        for p in peaks:
            peak_by_chr[p["chrom"]].append((p["start"], p["end"]))
        for chrom in peak_by_chr:
            peak_by_chr[chrom].sort()

        def overlaps(chrom, s, e):
            from bisect import bisect_left
            arr = peak_by_chr.get(chrom, [])
            if not arr:
                return False
            i = bisect_left([a[0] for a in arr], e)
            for j in range(max(0, i - 100), min(len(arr), i + 1)):
                if arr[j][1] >= s and arr[j][0] <= e:
                    return True
            return False

        for L in loops:
            a1 = overlaps(L["chrom1"], L["start1"], L["end1"])
            a2 = overlaps(L["chrom2"], L["start2"], L["end2"])
            anchor1_overlap += int(a1)
            anchor2_overlap += int(a2)
            both += int(a1 and a2)
        peak_overlap_stats[label] = {
            "n_loops": len(loops),
            "anchor1_overlap":  anchor1_overlap,
            "anchor2_overlap":  anchor2_overlap,
            "both_anchors":     both,
            "either_anchor":    sum(1 for L in loops
                                     if overlaps(L["chrom1"], L["start1"], L["end1"])
                                     or overlaps(L["chrom2"], L["start2"], L["end2"])),
            "peak_count":       len(peaks),
        }

    ts = timestamp()
    label = safe_label(args.label or Path(args.bedpe).stem + "_loops")
    out_dir = REPORT_DIR / f"{ts}_loops_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        if lengths:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4),
                                       constrained_layout=True)
            axes[0].hist(np.log10(np.array(lengths) + 1), bins=40,
                          color="#0072B2")
            axes[0].set_xlabel("log10(loop length, bp)")
            axes[0].set_ylabel("Count")
            axes[0].set_title("Loop length distribution (intra-chrom)")
            chr_counts = Counter(L["chrom1"] for L in intra)
            top = sorted(chr_counts.items(), key=lambda x: -x[1])[:24]
            axes[1].bar(range(len(top)), [c[1] for c in top],
                          color="#009E73")
            axes[1].set_xticks(range(len(top)))
            axes[1].set_xticklabels([c[0] for c in top], rotation=70,
                                       fontsize=7)
            axes[1].set_ylabel("Loops"); axes[1].set_title("Loops per chromosome")
            for ext in ("png", "svg"):
                fig.savefig(plot_dir / f"{label}_loop_qc.{ext}", dpi=200,
                              bbox_inches="tight", facecolor="white")
            plt.close(fig)
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_loops_report.md"
    lines = [
        f"# Loop / interaction analysis — `{Path(args.bedpe).name}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"- Total loops: **{len(loops):,}**",
        f"- Intra-chromosomal: **{len(intra):,}**",
        f"- Inter-chromosomal: **{len(inter):,}**",
    ]
    if lengths:
        lengths_arr = sorted(lengths)
        lines += [
            f"- Loop length median: **{lengths_arr[len(lengths_arr)//2]:,}** bp",
            f"  (min {min(lengths_arr):,}, max {max(lengths_arr):,})",
        ]
    if peak_overlap_stats:
        lines += ["", "## Anchor ↔ peak overlap", "",
                   "| peak set | peaks | anchor1 | anchor2 | both | either |",
                   "|---|---:|---:|---:|---:|---:|"]
        for label_, st in peak_overlap_stats.items():
            lines.append(
                f"| {label_} | {st['peak_count']:,} | "
                f"{st['anchor1_overlap']:,} | {st['anchor2_overlap']:,} | "
                f"{st['both_anchors']:,} | {st['either_anchor']:,} |"
            )
    lines += ["", f"![loop QC](Plots/{label}_loop_qc.png)", ""]
    report.write_text("\n".join(lines))
    print(f"Loops report: {report}")
    return report


# --------------------------- Motif enrichment ------------------------------

# Compact JASPAR-style PFM database (counts; treated as PFMs after
# pseudocount + normalisation in score_motif). Public-domain (JASPAR
# matrices are CC0).
_INLINE_PFMS = {
    "CTCF (MA0139.1)": [
        # A  C  G  T
        [ 87, 167, 281,  56],
        [ 99, 174, 132, 122],
        [125, 127, 207,  45],
        [ 60,  38, 387,  42],
        [ 40,  61, 388,  21],
        [114,  36, 290,  92],
        [197, 133, 117,  85],
        [ 91, 147, 273,  30],
        [114, 178,  74, 175],
        [ 30,  20, 460,  30],
        [277,  68, 191,  60],
        [ 23, 477,   1,  41],
        [ 14, 446,  43,  39],
        [ 40, 146,   8, 343],
        [ 32, 350,  37, 109],
        [ 62, 211,  72, 136],
        [122, 245, 102,  45],
        [109, 285,  26, 102],
        [ 86,  64, 314, 136],
    ],
    "AP-1 / FOS-JUN (MA0099.3)": [
        [ 18,  17, 165,  39],
        [ 23, 169,  18,  29],
        [219,  10,   2,   8],
        [  4,   7,  10, 218],
        [ 17,  29,   7, 186],
        [  3,   1, 235,   0],
        [148,  64,  18,   9],
        [ 41,  20, 128,  50],
        [  5,  62,  71, 101],
        [ 24,  56,  53, 106],
        [ 46,  18,  98,  77],
    ],
    "GATA1 (MA0035.4)": [
        [ 38,  14, 199,  10],
        [256,   1,   3,   1],
        [  1, 256,   0,   4],
        [  4,   0,   0, 257],
        [259,   0,   0,   2],
        [108,  16,  56,  81],
        [193,  12,  47,   9],
        [196,  19,  29,  17],
    ],
    "ETS / ELF1 (MA0473.4)": [
        [108,  73,  43,  37],
        [ 32,  18,   3, 208],
        [ 16,  20,   1, 224],
        [  5,  27, 224,   5],
        [  0,   1, 260,   0],
        [261,   0,   0,   0],
        [  0,   3, 257,   1],
        [ 41, 144,  20,  56],
        [ 87,  81,  25,  68],
    ],
    "NFKB / RELA (MA0107.1)": [
        [ 48,  34,  84,  10],
        [  6,   3, 161,   6],
        [  3,   3, 163,   7],
        [  1, 159,   3,  13],
        [ 90,  10,   2,  74],
        [  4,  15,   3, 154],
        [110,   5,  10,  51],
        [128,   3,  40,   5],
        [134,   3,  35,   4],
        [ 61,  34,  25,  56],
    ],
    "STAT1 (MA0137.3)": [
        [186,  24,  28,  62],
        [  3, 261,  31,   5],
        [ 20,   2,   0, 278],
        [288,   0,   0,  12],
        [ 23,   0,   1, 276],
        [ 24, 110, 117,  49],
        [299,   0,   0,   1],
        [ 27, 153,  54,  66],
        [ 86,  87,  44,  83],
    ],
    "FOXA1 (MA0148.4)": [
        [ 38,  31, 167,  38],
        [ 35,  10,  58, 171],
        [  4,   1,   0, 269],
        [  3, 268,   2,   1],
        [223,   1,   1,  49],
        [187,  16,  61,  10],
        [186,  36,  19,  33],
        [ 31, 150,  20,  73],
    ],
    "TP53 (MA0106.3)": [
        [ 15, 157,  30,  18],
        [120,   2,  87,  11],
        [  4, 200,   2,  14],
        [  6, 122,  87,   5],
        [ 73,  24,  16, 107],
        [ 57,  15, 135,  13],
        [  6, 112,  13,  89],
        [ 13, 168,  19,  20],
        [125,  20,  60,  15],
        [  1, 208,   1,  10],
    ],
    "MYC (MA0147.3)": [
        [ 38,  60,  40, 142],
        [  0, 277,   1,   2],
        [277,   0,   3,   0],
        [  3, 262,   2,  13],
        [  0,   3, 274,   3],
        [  0, 266,   3,  11],
        [ 19,  41, 197,  23],
        [ 23, 158,  46,  53],
    ],
    "SP1 / KLF (MA0079.5)": [
        [ 48, 206,  20,   6],
        [  4, 260,   1,  15],
        [  3, 263,   0,  14],
        [  2,  24, 247,   7],
        [  1, 270,   0,   9],
        [  2, 263,   1,  14],
        [ 33, 187,  28,  32],
        [ 46, 166,  43,  25],
    ],
}

_DNA_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3, "N": -1, "a": 0, "c": 1,
              "g": 2, "t": 3, "n": -1}


def _pfm_to_pwm(pfm, pseudocount: float = 0.5):
    import math
    total_counts = [sum(row) + 4 * pseudocount for row in pfm]
    pwm = []
    for row, total in zip(pfm, total_counts):
        log_odds = []
        for c in row:
            p = (c + pseudocount) / total
            log_odds.append(math.log2(p / 0.25))
        pwm.append(log_odds)
    return pwm


def _score_pwm(pwm, seq: str) -> float:
    best = -1e9
    L = len(pwm)
    for i in range(0, len(seq) - L + 1):
        s = 0.0
        ok = True
        for j in range(L):
            idx = _DNA_INDEX.get(seq[i + j], -1)
            if idx < 0:
                ok = False; break
            s += pwm[j][idx]
        if ok and s > best:
            best = s
    return best


def _revcomp(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A",
             "a": "t", "c": "g", "g": "c", "t": "a", "N": "N", "n": "n"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def _shuffle_seq_dinuc(seq: str, rng) -> str:
    """Mononucleotide-preserving shuffle (lighter than dinuc; sufficient
    for first-pass enrichment)."""
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def cmd_motif_enrichment(args: argparse.Namespace) -> Path:
    """Score peak sequences against a curated TF motif PFM set, compare
    to mononucleotide-shuffled controls, and report per-motif
    enrichment + a bar plot.

    Inputs:
      --bed        peaks BED
      --genome     genome FASTA (gzip OK if pyfaidx supports it)
      --top        cap on peak count for speed
      --score-cutoff   PWM log-odds cutoff to count a 'hit'
    """
    setup_logging(); mkdirs()
    try:
        from pyfaidx import Fasta
    except ImportError as e:
        raise RuntimeError(
            "Motif enrichment needs `pyfaidx`. "
            "Install with: pip install 'igvfagent[motif]'"
        ) from e
    import random
    rng = random.Random(args.seed)

    peaks, _ = parse_bed(Path(args.bed))
    if not peaks:
        raise SystemExit(f"No peaks in {args.bed}")
    if args.top and len(peaks) > args.top:
        peaks = peaks[: args.top]

    fa = Fasta(args.genome, sequence_always_upper=True, as_raw=True)
    seqs: "list[str]" = []
    for p in peaks:
        try:
            s = str(fa[p["chrom"]][p["start"]: p["end"]])
        except KeyError:
            continue
        if "N" in s and s.count("N") > 0.3 * len(s):
            continue
        seqs.append(s)
    if not seqs:
        raise SystemExit("No usable peak sequences (check chromosome "
                         "naming + genome FASTA).")

    pwm_cache = {name: _pfm_to_pwm(pfm) for name, pfm in _INLINE_PFMS.items()}
    cutoff = args.score_cutoff

    fg_hits = {name: 0 for name in pwm_cache}
    bg_hits = {name: 0 for name in pwm_cache}
    fg_total = len(seqs)
    bg_seqs = [_shuffle_seq_dinuc(s, rng) for s in seqs]

    for name, pwm in pwm_cache.items():
        for s in seqs:
            top_score = max(_score_pwm(pwm, s), _score_pwm(pwm, _revcomp(s)))
            if top_score >= cutoff:
                fg_hits[name] += 1
        for s in bg_seqs:
            top_score = max(_score_pwm(pwm, s), _score_pwm(pwm, _revcomp(s)))
            if top_score >= cutoff:
                bg_hits[name] += 1

    rows: "list[dict]" = []
    for name in pwm_cache:
        fg = fg_hits[name]; bg = bg_hits[name]
        # Fisher's exact (right-tailed) approximation via log odds + p-est
        a, b = fg, fg_total - fg
        c, d = bg, fg_total - bg
        # avoid div by 0
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        # quick chi-sq stat
        n = a + b + c + d
        expected_a = (a + b) * (a + c) / n
        chi = ((a - expected_a) ** 2) / max(expected_a, 0.5)
        rows.append({
            "motif":   name,
            "fg_hits": fg, "bg_hits": bg,
            "fg_rate": fg / fg_total,
            "bg_rate": bg / fg_total,
            "log2_odds": math.log2(max(odds, 1e-9)),
            "chi_sq":  chi,
        })
    rows.sort(key=lambda r: -r["log2_odds"])

    ts = timestamp()
    label = safe_label(args.label or Path(args.bed).stem + "_motifs")
    out_dir = REPORT_DIR / f"{ts}_motif_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / f"{label}_motif_enrichment.csv", rows,
              cols=["motif", "fg_hits", "bg_hits", "fg_rate", "bg_rate",
                     "log2_odds", "chi_sq"])

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        names = [r["motif"] for r in rows]
        odds = np.array([r["log2_odds"] for r in rows])
        plt.figure(figsize=(8, 0.35 * len(rows) + 1))
        colors = ["#D55E00" if x > 0 else "#0072B2" for x in odds]
        plt.barh(range(len(rows)), odds, color=colors)
        plt.yticks(range(len(rows)), names, fontsize=8)
        plt.gca().invert_yaxis()
        plt.axvline(0, color="black", lw=0.5)
        plt.xlabel("log2(odds ratio) vs shuffled background")
        plt.title(f"Motif enrichment — {label}")
        plt.tight_layout()
        for ext in ("png", "svg"):
            plt.savefig(plot_dir / f"{label}_motif_enrichment.{ext}",
                          dpi=200, bbox_inches="tight", facecolor="white")
        plt.close()
    except Exception as e:
        logging.warning("Plot skipped (%s)", e)

    report = out_dir / f"{label}_motif_report.md"
    report.write_text(
        f"# Motif enrichment — `{Path(args.bed).name}`\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Peaks scanned: **{len(seqs):,}**\n"
        f"- Genome FASTA: `{args.genome}`\n"
        f"- PWM cutoff (log2 odds): **{cutoff}**\n"
        f"- PFM database: bundled JASPAR-derived "
        f"({len(_INLINE_PFMS)} motifs)\n\n"
        + "| motif | fg | bg | fg/bg | log2-odds |\n"
        + "|---|---:|---:|---:|---:|\n"
        + "\n".join(
            f"| {r['motif']} | {r['fg_hits']} | {r['bg_hits']} | "
            f"{(r['fg_rate'] / max(r['bg_rate'], 1e-9)):.2f} | "
            f"{r['log2_odds']:+.3f} |"
            for r in rows
        )
        + f"\n\n![motif enrichment](Plots/{label}_motif_enrichment.png)\n"
        + "\n\nNotes: this is a lightweight first-pass scanner over a "
          "curated set of common TF motifs (CTCF, AP-1, GATA1, ETS, "
          "NFkB, STAT1, FOXA1, TP53, MYC, SP1/KLF). For a publication-"
          "grade motif analysis use HOMER (`findMotifsGenome.pl`) or "
          "MEME/MEME-ChIP with their full motif libraries.\n"
    )
    print(f"Motif report: {report}")
    return report


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
        "### `bigwig-frip` — fraction-of-signal-in-peaks (needs `[hic]` extras)",
        "",
        "```bash",
        "pip install 'igvfagent[hic]'",
        "igvfagent encode bigwig-frip --bigwig signal.bw --bed peaks.bed",
        "```",
        "",
        "Returns the global FRiP plus per-chromosome breakdown. ENCODE's "
        "FRiP minimum is 0.01 for ChIP-seq; > 0.05 is good. Histone "
        "ChIP-seq tends to score lower because the signal is broader.",
        "",
        "### `bigwig-tss-heatmap` — anchor-centered signal heatmap (needs `[hic]`)",
        "",
        "```bash",
        "igvfagent encode bigwig-tss-heatmap \\",
        "    --bigwig signal.bw --anchor-bed tss.bed \\",
        "    --window 4000 --bins 200 --max-anchors 5000 --label tss_h3k27ac",
        "```",
        "",
        "Aggregates bigWig signal over `--window` bp centered on each "
        "anchor (TSS BED, peak summits, etc.), sorts rows by total "
        "signal, and renders a heatmap + meta-profile in one figure.",
        "",
        "### `hic-matrix` — Hi-C contact heatmap (needs `[hic]`)",
        "",
        "```bash",
        "pip install 'igvfagent[hic]'",
        "igvfagent encode hic-matrix \\",
        "    --input ENCFF000XYZ.mcool --region chr19:44900000-45100000 \\",
        "    --resolution 10000 --balance --label apoe_hic",
        "igvfagent encode hic-matrix \\",
        "    --input ENCFF111ABC.hic --region chr19:44900000-45100000 \\",
        "    --resolution 10000",
        "```",
        "",
        "Loads the contact matrix for a region from `.mcool` (cooler) or "
        "`.hic` (hic-straw), saves the raw matrix as `.npy`, and renders "
        "a `log1p`-scaled contact heatmap PNG/SVG.",
        "",
        "### `hic-insulation` — Crane-style TAD-boundary score (needs `[hic]`)",
        "",
        "```bash",
        "igvfagent encode hic-insulation \\",
        "    --input ENCFF000XYZ.mcool --region chr19:44000000-46000000 \\",
        "    --resolution 10000 --window 200000 --balance \\",
        "    --boundary-threshold -0.3 --label apoe_ins",
        "```",
        "",
        "Computes the insulation score along the diagonal, identifies "
        "local minima below `--boundary-threshold` as candidate TAD "
        "boundaries, and emits a BED of boundaries plus a 1-D track "
        "plot. For loop-level calls, run a dedicated tool (HiCCUPS / "
        "Mustache / Peakachu) and feed the resulting bedpe into "
        "`loops-analyze`.",
        "",
        "### `loops-analyze` — Hi-C / ChIA-PET / capture Hi-C loops",
        "",
        "```bash",
        "igvfagent encode loops-analyze --bedpe loops.bedpe \\",
        "    --peaks 'CTCF peaks:ctcf.bed' --peaks 'H3K27ac peaks:h3k27ac.bed' \\",
        "    --label gm12878_loops",
        "```",
        "",
        "Parses a `.bedpe` and reports total / intra / inter / median "
        "loop length plus optional anchor ↔ peak overlap stats per "
        "input peak set. Works on outputs from any loop caller (Mustache, "
        "HiCCUPS, MAPS, Fit-Hi-C, ChIA-PET clusters, etc.).",
        "",
        "### `motif-enrichment` — TF motif enrichment in peak sequences",
        "",
        "```bash",
        "pip install 'igvfagent[motif]'",
        "# Download a genome FASTA from UCSC (one-time):",
        "#   curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O",
        "igvfagent encode motif-enrichment \\",
        "    --bed peaks.bed --genome hg38.fa.gz \\",
        "    --top 2000 --score-cutoff 8 --label k562_h3k27ac_motifs",
        "```",
        "",
        "Scans peak sequences against a curated, bundled JASPAR-derived "
        "PFM set (CTCF, AP-1, GATA1, ETS, NFkB, STAT1, FOXA1, TP53, "
        "MYC, SP1/KLF) and reports log2 odds enrichment vs "
        "mononucleotide-shuffled background. For full TF coverage use "
        "HOMER `findMotifsGenome.pl` or MEME-ChIP.",
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

    # ---- bigWig signal subcommands ----
    s = sub.add_parser("bigwig-frip",
                        help="FRiP from a bigWig + peak BED.")
    s.add_argument("--bigwig", required=True)
    s.add_argument("--bed", required=True)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_bigwig_frip)

    s = sub.add_parser("bigwig-tss-heatmap",
                        help="bigWig signal heatmap centered on anchors "
                             "(TSS BED, peak summits, etc.).")
    s.add_argument("--bigwig", required=True)
    s.add_argument("--anchor-bed", required=True)
    s.add_argument("--window", type=int, default=4000,
                    help="Window size in bp (centered on each anchor).")
    s.add_argument("--bins", type=int, default=200)
    s.add_argument("--max-anchors", type=int, default=5000)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_bigwig_tss_heatmap)

    # ---- Hi-C subcommands ----
    s = sub.add_parser("hic-matrix",
                        help="Render a contact heatmap for a region from "
                             ".mcool / .hic.")
    s.add_argument("--input", required=True,
                    help="Path to .mcool or .hic file.")
    s.add_argument("--region", required=True,
                    help="chr19:44900000-45100000")
    s.add_argument("--resolution", type=int, default=10000)
    s.add_argument("--balance", action="store_true",
                    help=".mcool: use ICE-balanced matrix; .hic: KR.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_hic_matrix)

    s = sub.add_parser("hic-insulation",
                        help="Crane-style insulation score across a region.")
    s.add_argument("--input", required=True)
    s.add_argument("--region", required=True)
    s.add_argument("--resolution", type=int, default=10000)
    s.add_argument("--window", type=int, default=200000,
                    help="Insulation window in bp.")
    s.add_argument("--balance", action="store_true")
    s.add_argument("--boundary-threshold", type=float, default=-0.3,
                    help="Local-minimum threshold for boundary calls.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_hic_insulation)

    # ---- bedpe loop / interaction analysis ----
    s = sub.add_parser("loops-analyze",
                        help="Analyze a .bedpe loops / interactions file "
                             "(Hi-C, ChIA-PET, capture Hi-C). Optionally "
                             "overlay anchors with peak BEDs.")
    s.add_argument("--bedpe", required=True)
    s.add_argument("--peaks", action="append", default=None,
                    help="Optional `LABEL:PATH` peak BED for anchor "
                         "overlap. Repeatable.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_loops_analyze)

    # ---- motif enrichment ----
    s = sub.add_parser("motif-enrichment",
                        help="JASPAR-derived TF motif enrichment in peak "
                             "sequences (foreground vs shuffled "
                             "background).")
    s.add_argument("--bed", required=True)
    s.add_argument("--genome", required=True,
                    help="Indexed genome FASTA (.fa / .fa.gz). Use UCSC "
                         "GRCh38: hg38.fa.gz")
    s.add_argument("--top", type=int, default=2000)
    s.add_argument("--score-cutoff", type=float, default=8.0,
                    help="PWM log2-odds cutoff for a 'hit' (default 8).")
    s.add_argument("--seed", type=int, default=42)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_motif_enrichment)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/ENCODE_PIPELINE_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
