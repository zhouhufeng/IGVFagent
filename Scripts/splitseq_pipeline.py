#!/usr/bin/env python3
"""SPLiT-seq pipeline for IGVF Portal data.

Discovers, retrieves, processes, and analyzes Parse-Biosciences-style
combinatorial-barcoding single-nucleus RNA-seq datasets (Parse SPLiT-seq /
SPLiT-seq) hosted on the IGVF Portal. The pipeline is opinionated about
SPLiT-seq's quirks: data come in **multiplexed sub-pools** with tens of
donors per pool, files are typically distributed as ``tar.gz`` bundles of
MatrixMarket sparse-count matrices, and the analytical value is in
**per-strain / per-donor / per-pool** comparisons.

Subcommands

  retrieve         Discover SPLiT-seq AnalysisSets in the IGVF Portal,
                   write a manifest, optionally download files.
  manifest         Hydrate metadata for one accession (or a list of
                   accessions) and write a per-pool / per-donor manifest.
  download         Download sparse-matrix and cell-annotation files for a
                   manifest with a configurable size cap.
  process-local    Load downloaded matrices into per-pool AnnData objects
                   and concatenate; persist as ``.h5ad``.
  analyze          End-to-end: load -> QC -> normalize -> integrate -> UMAP
                   -> cluster -> annotate -> per-strain composition.
  plot             Produce publication-style multi-panel figures from a
                   processed AnnData.
  compare-strains  Per-cell-type strain effect: pseudobulk + Wilcoxon DEG
                   between mouse founder strains, with volcano plots.
  demux-script     Emit a ready-to-run shell script for genotype-based
                   donor demultiplexing (souporcell / vireo). The skill
                   does not run the demultiplexer itself; it scaffolds it.
  write-playbook   Emit ``Docs/Skills/SPLITSEQ_SKILLS.md``.

Outputs follow the project layout:

  Data/IGVF/SPLiTseq/Metadata/         hydrated portal JSON
  Data/IGVF/SPLiTseq/Downloads/        per-AnalysisSet bundles
  Data/Manifests/SPLiTseq/             dataset / pool / donor manifests
  Data/Cache/SPLiTseq/                 in-process AnnData (.h5ad)
  Docs/SPLiTseq/                       per-run reports + plots
  Docs/Logs/                           runtime logs
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import io
import json
import logging
import os
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "SPLiTseq"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "SPLiTseq"
METADATA_DIR = DATA_DIR / "IGVF" / "SPLiTseq" / "Metadata"
DOWNLOAD_DIR = DATA_DIR / "IGVF" / "SPLiTseq" / "Downloads"
H5AD_CACHE_DIR = DATA_DIR / "Cache" / "SPLiTseq"

PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
PORTAL_PUBLIC_BASE = _resolve_endpoint("portal", "IGVF_PORTAL_PUBLIC_BASE")

# Canonical Collaborative Cross / Diversity Outbred founder strains.
FOUNDERS_8CUBE = (
    "A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ",
    "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ",
)
FOUNDER_COLORS = {
    "A/J":          "#D55E00",
    "C57BL/6J":     "#000000",
    "129S1/SvImJ":  "#0072B2",
    "NOD/ShiLtJ":   "#E69F00",
    "NZO/HlLtJ":    "#009E73",
    "CAST/EiJ":     "#7B5BA6",
    "PWK/PhJ":      "#A14586",
    "WSB/EiJ":      "#56B4E9",
}

# Mouse marker gene panels per tissue. Used by `analyze` when no shipped
# cell annotation file is available.
MOUSE_MARKERS: dict[str, dict[str, list[str]]] = {
    "gonad": {
        "Germ cell":          ["Ddx4", "Dazl", "Sycp3", "Stra8"],
        "Sertoli":            ["Sox9", "Dhh", "Amh", "Wt1"],
        "Leydig":             ["Cyp17a1", "Star", "Cyp11a1", "Hsd3b1"],
        "Granulosa":          ["Foxl2", "Inha", "Inhba", "Cyp19a1"],
        "Theca":              ["Cyp17a1", "Lhcgr"],
        "Stromal":            ["Pdgfra", "Tcf21", "Nr2f2"],
        "Endothelial":        ["Pecam1", "Cdh5", "Kdr"],
        "Macrophage":         ["Cd68", "Csf1r", "Adgre1"],
    },
    "adrenal": {
        "Cortex zG":          ["Cyp11b2", "Dab2"],
        "Cortex zF":           ["Cyp11b1", "Cyp17a1"],
        "Cortex zR":           ["Cyb5a", "Sult2a1"],
        "Medulla":             ["Th", "Pnmt", "Chga", "Dbh"],
        "Capsule/Stromal":     ["Acta2", "Myh11", "Pdgfra"],
        "Endothelial":         ["Pecam1", "Cdh5"],
        "Macrophage":          ["Cd68", "Csf1r"],
    },
    "brain": {
        "Excitatory neuron":   ["Slc17a7", "Slc17a6", "Camk2a"],
        "Inhibitory neuron":   ["Gad1", "Gad2", "Slc32a1"],
        "Astrocyte":           ["Aqp4", "Gfap", "Slc1a2", "Aldh1l1"],
        "Oligodendrocyte":     ["Plp1", "Mog", "Mobp", "Mbp"],
        "OPC":                 ["Pdgfra", "Cspg4"],
        "Microglia":           ["Csf1r", "Cx3cr1", "P2ry12"],
        "Endothelial":         ["Pecam1", "Cldn5"],
        "Pericyte":            ["Pdgfrb", "Rgs5"],
    },
    "liver": {
        "Hepatocyte":          ["Alb", "Cyp3a11", "Hnf4a"],
        "Cholangiocyte":       ["Krt19", "Sox9", "Epcam"],
        "Endothelial":         ["Pecam1", "Stab2"],
        "Stellate":            ["Lrat", "Pdgfrb"],
        "Kupffer":             ["Clec4f", "Csf1r"],
    },
    "heart": {
        "Cardiomyocyte":       ["Tnnt2", "Myh6", "Myh7"],
        "Fibroblast":          ["Pdgfra", "Tcf21", "Col1a1"],
        "Endothelial":         ["Pecam1", "Cdh5"],
        "Smooth muscle":       ["Acta2", "Myh11"],
        "Macrophage":          ["Cd68", "Csf1r"],
    },
    "kidney": {
        "Proximal tubule":     ["Slc34a1", "Lrp2"],
        "Distal tubule":       ["Slc12a3", "Wnk1"],
        "Loop of Henle":       ["Slc12a1", "Umod"],
        "Collecting duct":     ["Aqp2", "Slc4a1"],
        "Podocyte":            ["Nphs1", "Nphs2"],
        "Endothelial":         ["Pecam1", "Cdh5"],
        "Stromal":             ["Pdgfra"],
    },
    "muscle": {
        "Type I myofiber":     ["Myh7"],
        "Type II myofiber":    ["Myh1", "Myh4"],
        "Satellite cell":      ["Pax7", "Myf5"],
        "Endothelial":         ["Pecam1"],
        "Fibroblast":          ["Pdgfra", "Tcf21"],
    },
}

# Conservative SPLiT-seq QC defaults; tunable from the CLI.
DEFAULT_QC = {
    "min_genes":     500,    # genes detected per cell
    "max_genes":     8000,
    "min_counts":    1000,   # UMIs per cell
    "max_counts":    50000,
    "max_pct_mito":  10.0,   # %
    "min_cells":     3,      # gene must be expressed in at least N cells
}


# ----------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"splitseq_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
              DOWNLOAD_DIR, H5AD_CACHE_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".")
                   else "_" for ch in label)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# ------------------------------- HTTP helpers --------------------------------

def request_headers() -> dict[str, str]:
    h = {"Accept": "application/json,*/*",
         "User-Agent": "IGVFdataAgent/0.1 splitseq-pipeline"}
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        h["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return h


def fetch_json(url: str) -> tuple[int, Any]:
    logging.info("GET %s", url)
    req = urllib.request.Request(url, headers=request_headers(), method="GET")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            content = resp.read()
            if "json" in resp.headers.get("Content-Type", ""):
                return resp.status, json.loads(content)
            return resp.status, {"text_response": content.decode(errors="replace"),
                                 "url": url}
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace")
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            data = {"http_error_body": body, "url": url}
        return e.code, data
    except urllib.error.URLError as e:
        return 0, {"network_error": str(e.reason), "url": url}


def build_url(base: str, path: str, params: dict | None = None) -> str:
    p = path if path.startswith("/") else f"/{path}"
    url = f"{base}{p}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def download_file(url: str, dest: Path, max_bytes: int | None = None) -> int:
    """Stream a URL to disk. Returns bytes written; raises on truncation."""
    logging.info("Download %s -> %s", url, dest)
    req = urllib.request.Request(url, headers=request_headers(), method="GET")
    written = 0
    with urllib.request.urlopen(req, timeout=120) as resp, dest.open("wb") as f:
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


# ------------------------------- Portal queries ------------------------------

def search_splitseq_analysissets(limit: int = 50,
                                  lab: str | None = None,
                                  sample_type: str | None = None,
                                  status: str = "released",
                                  taxa: str | None = None) -> list[dict]:
    """Search IGVF AnalysisSets for Parse SPLiT-seq."""
    params: dict[str, Any] = {
        "type": "AnalysisSet",
        "preferred_assay_titles": "Parse SPLiT-seq",
        "format": "json",
        "limit": limit,
        "status": status,
    }
    if lab:
        params["lab.title"] = lab
    if sample_type:
        params["samples.sample_terms.term_name"] = sample_type
    if taxa:
        params["donors.taxa"] = taxa
    url = build_url(PORTAL_API_BASE, "/search/", params)
    status_code, data = fetch_json(url)
    if status_code != 200:
        logging.warning("Portal search HTTP %s", status_code)
        return []
    if isinstance(data, dict):
        return [g for g in data.get("@graph", []) if isinstance(g, dict)]
    return []


def hydrate_analysisset(accession: str) -> dict:
    url = build_url(PORTAL_API_BASE, f"/analysis-sets/{accession}/",
                    {"format": "json"})
    status_code, data = fetch_json(url)
    if status_code == 200 and isinstance(data, dict):
        return data
    logging.warning("Hydrate %s -> HTTP %s", accession, status_code)
    return {"accession": accession, "_http_status": status_code}


def hydrate_files_for_analysisset(meta: dict) -> list[dict]:
    """Inline files in the hydrated AnalysisSet JSON; resolve any @id-only refs."""
    files = meta.get("files") or []
    out = []
    for f in files:
        if isinstance(f, dict) and "accession" in f:
            out.append(simplify_file(f))
            continue
        if isinstance(f, str):
            url = build_url(PORTAL_API_BASE, f"{f}", {"format": "json"})
            sc, data = fetch_json(url)
            if sc == 200 and isinstance(data, dict):
                out.append(simplify_file(data))
    return out


def simplify_file(f: dict) -> dict:
    href = f.get("href") or ""
    public_url = (PORTAL_API_BASE + href) if href else ""
    return {
        "accession":         f.get("accession", ""),
        "content_type":      f.get("content_type", ""),
        "file_format":       f.get("file_format", ""),
        "assembly":          f.get("assembly", ""),
        "controlled_access": f.get("controlled_access", ""),
        "file_size_bytes":   f.get("file_size", ""),
        "file_size_gb":      round((f.get("file_size") or 0) / (1024**3), 3),
        "md5sum":            f.get("md5sum", ""),
        "s3_uri":            f.get("s3_uri", ""),
        "url":               public_url,
        "status":            f.get("status", ""),
        "summary":           f.get("summary", "") or f.get("description", ""),
    }


# ------------------------------- Manifest writing ----------------------------

def write_csv(path: Path, rows: list[dict], cols: list[str] | None = None):
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def summarize_dataset(meta: dict) -> dict:
    samples = meta.get("samples") or []
    donors = meta.get("donors") or []
    sample_acc = []
    for s in samples:
        if isinstance(s, dict):
            sample_acc.append(s.get("accession", ""))
        elif isinstance(s, str):
            sample_acc.append(s.rsplit("/", 2)[-2] if s.endswith("/") else s)
    donor_acc = []
    donor_strains: list[str] = []
    for d in donors:
        if isinstance(d, dict):
            donor_acc.append(d.get("accession", ""))
            for s in (d.get("strain_background") or []):
                donor_strains.append(s)
            if d.get("strain_background_term") and not d.get("strain_background"):
                donor_strains.append(d.get("strain_background_term"))
    lab = meta.get("lab") or {}
    return {
        "accession":              meta.get("accession", ""),
        "file_set_type":          meta.get("file_set_type", ""),
        "preferred_assay_titles": ", ".join(meta.get("preferred_assay_titles", []) or []),
        "assay_titles":           ", ".join(meta.get("assay_titles", []) or []),
        "summary":                meta.get("summary", ""),
        "description":            (meta.get("description") or "").replace("\n", " "),
        "lab":                    (lab.get("title") if isinstance(lab, dict) else str(lab)),
        "status":                 meta.get("status", ""),
        "release_timestamp":      meta.get("release_timestamp", ""),
        "n_samples":              len(samples),
        "samples":                "; ".join(sample_acc[:10]),
        "n_donors":               len(donors),
        "donor_accessions":       "; ".join(donor_acc[:10]),
        "founder_strains":        "; ".join(sorted(set(donor_strains))[:8]),
        "n_input_file_sets":      len(meta.get("input_file_sets") or []),
        "url":                    (PORTAL_PUBLIC_BASE
                                    + "/analysis-sets/" + meta.get("accession", "")),
    }


def per_pool_rows(meta: dict, files: list[dict]) -> list[dict]:
    """Cross-product rows per (file × any plate / GEM-well metadata in
    aliases). For SPLiT-seq the alias often encodes ``..._SubpoolN_...``
    or ``Plate_XX_GY``; we surface that here."""
    aliases = meta.get("aliases") or []
    pool_label = ""
    for a in aliases:
        m = re.search(r"(Subpool[_-]?\d+|plate[_-]?\d+|GEM[_-]?\w+)", a, re.I)
        if m:
            pool_label = m.group(0)
            break
    out = []
    for f in files:
        out.append({
            **f,
            "analysis_set_accession": meta.get("accession", ""),
            "pool_label":             pool_label,
            "alias":                  aliases[0] if aliases else "",
        })
    return out


# ------------------------------- Subcommands ---------------------------------

def cmd_retrieve(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    hits = search_splitseq_analysissets(
        limit=args.limit, lab=args.lab,
        sample_type=args.sample_type, status=args.status, taxa=args.taxa,
    )
    if args.fetch_file_details:
        for h in hits:
            acc = h.get("accession")
            if not acc:
                continue
            meta = hydrate_analysisset(acc)
            (METADATA_DIR / f"{acc}.json").write_text(json.dumps(meta, indent=2))
            time.sleep(0.15)
    summaries = []
    for h in hits:
        acc = h.get("accession")
        meta_path = METADATA_DIR / f"{acc}.json"
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else h
        summaries.append(summarize_dataset(meta))
    cols = ["accession", "file_set_type", "preferred_assay_titles",
            "assay_titles", "summary", "lab", "status", "release_timestamp",
            "n_samples", "samples", "n_donors", "donor_accessions",
            "founder_strains", "n_input_file_sets", "url", "description"]
    label = safe_label(args.label or "splitseq_corpus")
    manifest = MANIFEST_DIR / f"{ts}_{label}_dataset_manifest.csv"
    write_csv(manifest, summaries, cols)
    logging.info("Wrote: %s", manifest)
    print(f"Manifest: {manifest}")
    print(f"Datasets: {len(summaries)}")
    return manifest


def cmd_manifest(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    accs = [a.strip() for a in args.accessions.split(",") if a.strip()]
    label = safe_label(args.label or accs[0] if accs else "manifest")
    summaries: list[dict] = []
    pool_rows: list[dict] = []
    for acc in accs:
        meta = hydrate_analysisset(acc)
        (METADATA_DIR / f"{acc}.json").write_text(json.dumps(meta, indent=2))
        files = hydrate_files_for_analysisset(meta)
        summaries.append(summarize_dataset(meta))
        pool_rows.extend(per_pool_rows(meta, files))
    summary_path = MANIFEST_DIR / f"{ts}_{label}_dataset_summary.csv"
    files_path   = MANIFEST_DIR / f"{ts}_{label}_files_per_pool.csv"
    write_csv(summary_path, summaries)
    write_csv(files_path, pool_rows,
              cols=["analysis_set_accession", "pool_label", "alias",
                     "accession", "content_type", "file_format",
                     "file_size_gb", "file_size_bytes", "status",
                     "controlled_access", "url", "summary"])
    print(f"Summary manifest:  {summary_path}")
    print(f"Files-per-pool:    {files_path}")
    print(f"Datasets summarized: {len(summaries)}")
    return files_path


def cmd_download(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    rows = list(csv.DictReader(open(args.manifest)))
    if not rows:
        raise SystemExit("Empty manifest")
    cap = int(args.max_download_gb * (1024**3))
    only_kinds = set(args.only or [])
    download_log: list[dict] = []
    used = 0
    for r in rows:
        kind = (r.get("content_type") or "").lower()
        if only_kinds and kind not in only_kinds:
            continue
        size = int(r.get("file_size_bytes") or 0)
        if used + size > cap:
            logging.info("Skipping %s — would exceed size cap", r.get("accession"))
            continue
        ds_acc = r.get("analysis_set_accession") or "misc"
        out_dir = DOWNLOAD_DIR / safe_label(ds_acc)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / safe_label(
            (r.get("accession") or "file") + "." + (r.get("file_format") or "bin")
        )
        if out_file.exists() and out_file.stat().st_size == size:
            logging.info("Cached %s", out_file)
            written = size
        else:
            try:
                written = download_file(r["url"], out_file, max_bytes=size + (1 << 20))
            except Exception as e:
                logging.warning("Failed %s: %s", r.get("url"), e)
                continue
        used += written
        download_log.append({**r, "downloaded_path": str(out_file),
                              "downloaded_bytes": written})
    log_path = MANIFEST_DIR / f"{timestamp()}_download_log.csv"
    write_csv(log_path, download_log)
    print(f"Downloaded {len(download_log)} files / {used / (1024**3):.2f} GB")
    print(f"Log: {log_path}")
    return log_path


# ------------------------------- Local processing ----------------------------

def find_dataset_bundles() -> dict[str, Path]:
    """Return mapping of analysis-set accession -> local download dir."""
    out = {}
    if not DOWNLOAD_DIR.exists():
        return out
    for d in sorted(DOWNLOAD_DIR.iterdir()):
        if d.is_dir() and d.name.startswith("IGVFDS"):
            out[d.name] = d
    return out


def find_mtx_bundle(directory: Path) -> tuple[Path, Path, Path] | None:
    """Locate (mtx, features, barcodes) inside a downloaded folder, untarring
    a single tar.gz if found.

    Handles both the typical Parse SPLiT-seq layout
    (``BANN1632_PUT_counts.mtx`` + features.tsv.gz + barcodes.tsv.gz) and
    a generic ``matrix.mtx`` 10x-style layout.
    """
    tars = sorted(list(directory.glob("*.tar.gz")) + list(directory.glob("*.tar")))
    for t in tars:
        try:
            with tarfile.open(t) as tf:
                tf.extractall(directory)
        except Exception as e:
            logging.warning("Could not untar %s: %s", t, e)
    mtx = next(iter(sorted(directory.glob("*.mtx"))), None)
    if not mtx:
        # 10x-style?
        for sub in directory.glob("**/matrix.mtx*"):
            mtx = sub; break
    if not mtx:
        return None
    parent = mtx.parent
    feats = next(iter(sorted(parent.glob("*features.tsv*")) +
                       sorted(parent.glob("genes.tsv*"))), None)
    bcs = next(iter(sorted(parent.glob("*barcodes.tsv*"))), None)
    if not feats or not bcs:
        return None
    return mtx, feats, bcs


def load_mtx_to_anndata(mtx: Path, feats: Path, bcs: Path,
                         sample_id: str) -> "ad.AnnData":  # noqa: F821
    import anndata as ad  # heavy import deferred
    import pandas as pd
    import scipy.io
    logging.info("Loading mtx %s for %s", mtx.name, sample_id)
    raw = scipy.io.mmread(str(mtx)).tocsr()
    n_features, n_barcodes = raw.shape
    # mtx is features × barcodes; transpose to cells × genes
    raw = raw.T.tocsr()
    bc_df = pd.read_csv(bcs, sep="\t", header=None, names=["barcode"])
    fe_df = pd.read_csv(feats, sep="\t", header=None)
    if fe_df.shape[1] >= 2:
        fe_df.columns = ["gene_id", "gene"] + list(fe_df.columns[2:])
    else:
        fe_df.columns = ["gene"]
    bc_df["barcode_full"] = sample_id + ":" + bc_df["barcode"]
    bc_df = bc_df.set_index("barcode_full")
    var = fe_df.set_index("gene") if "gene" in fe_df.columns else fe_df
    # collapse duplicate gene symbols
    if var.index.duplicated().any():
        keep = ~var.index.duplicated(keep="first")
        raw = raw[:, keep.to_numpy()]
        var = var.loc[keep]
    a = ad.AnnData(X=raw, obs=bc_df, var=var)
    a.obs["sample_id"] = sample_id
    return a


def attach_annotation_tsv(a: "ad.AnnData",  # noqa: F821
                            ann_tsv: Path) -> "ad.AnnData":  # noqa: F821
    import pandas as pd
    df = pd.read_csv(ann_tsv, sep="\t")
    if "cell_barcode" in df.columns:
        df["bc_full"] = df["analysis_set_accession"].astype(str) + ":" + df["cell_barcode"]
        df = df.set_index("bc_full")
    elif "barcode" in df.columns:
        df = df.set_index("barcode")
    common = a.obs_names.intersection(df.index)
    logging.info("Annotation join: matrix=%d, ann=%d, common=%d",
                  a.n_obs, len(df), len(common))
    if len(common) == 0:
        return a
    a = a[common, :].copy()
    df = df.loc[common]
    for c in df.columns:
        if c not in a.obs.columns:
            a.obs[c] = df[c].values
    return a


def cmd_process_local(args: argparse.Namespace) -> Path:
    """Load downloaded SPLiT-seq matrices into a single AnnData."""
    setup_logging(); mkdirs()
    import anndata as ad
    bundles = find_dataset_bundles()
    if args.accessions:
        keep = set(args.accessions.split(","))
        bundles = {k: v for k, v in bundles.items() if k in keep}
    if not bundles:
        raise SystemExit("No downloaded SPLiT-seq bundles found. "
                         "Run `download` first.")
    adatas = []
    for acc, dirpath in bundles.items():
        triple = find_mtx_bundle(dirpath)
        if not triple:
            logging.warning("Skipping %s — no usable matrix bundle", acc)
            continue
        mtx, feats, bcs = triple
        a = load_mtx_to_anndata(mtx, feats, bcs, acc)
        ann = next(iter(dirpath.glob("*cell*annotation*tsv*"))
                    if dirpath.exists() else iter(()), None)
        if not ann:
            ann = next(iter(dirpath.glob("*.tsv.gz")), None)
        if ann and ann.suffix.lower() in (".gz", ".tsv"):
            a = attach_annotation_tsv(a, ann)
        adatas.append(a)
    if not adatas:
        raise SystemExit("No matrices loaded.")
    if len(adatas) > 1:
        full = ad.concat(adatas, join="outer", label="dataset",
                          keys=list(bundles.keys()))
    else:
        full = adatas[0]
    label = safe_label(args.label or "splitseq_pool")
    out = H5AD_CACHE_DIR / f"{label}.h5ad"
    full.write_h5ad(out)
    print(f"Wrote AnnData: {out} (cells={full.n_obs}, genes={full.n_vars})")
    print(f"obs cols: {list(full.obs.columns)[:25]}")
    return out


# ------------------------------- Analysis ------------------------------------

def run_qc_normalize(a, qc=DEFAULT_QC):
    import scanpy as sc
    a.var["mt"] = a.var_names.str.lower().str.startswith("mt-") | \
                   a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True,
                                percent_top=None, log1p=False)
    pre = a.n_obs
    # Stash pre-filter distributions so cmd_plot can render the
    # pre-vs-post QC violin overlay (Chipeyown / Fowler-lab style).
    import numpy as np
    a.uns["qc_pre"] = {
        "n_genes":       np.asarray(a.obs["n_genes_by_counts"], dtype=float),
        "total_counts":  np.asarray(a.obs["total_counts"], dtype=float),
        "pct_counts_mt": np.asarray(a.obs["pct_counts_mt"], dtype=float),
        "thresholds":    {k: float(v) for k, v in qc.items()},
        "n_pre":         int(pre),
    }
    keep = (
        (a.obs["n_genes_by_counts"] >= qc["min_genes"]) &
        (a.obs["n_genes_by_counts"] <= qc["max_genes"]) &
        (a.obs["total_counts"]      >= qc["min_counts"]) &
        (a.obs["total_counts"]      <= qc["max_counts"]) &
        (a.obs["pct_counts_mt"]     <= qc["max_pct_mito"])
    )
    a = a[keep].copy()
    sc.pp.filter_genes(a, min_cells=qc["min_cells"])
    logging.info("QC kept %d / %d cells", a.n_obs, pre)
    a.uns["qc_pre"]["n_post"] = int(a.n_obs)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2500, flavor="seurat_v3"
                                  if "counts" in a.layers else "seurat",
                                  layer="counts" if "counts" in a.layers else None)
    # Keep the full log-normalized matrix in .raw (for marker scoring), then
    # restrict to HVGs BEFORE scaling. sc.pp.scale densifies its input, so
    # scaling all genes on a large atlas (e.g. Rosenberg's 156k×27k) needs
    # tens of GB; scaling only the 2,500 HVGs keeps peak memory ~100× smaller.
    a.raw = a
    a = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50)
    return a


def integrate_pools(a, batch_key="dataset"):
    import scanpy as sc
    if a.obs[batch_key].nunique() < 2:
        sc.pp.neighbors(a, n_pcs=40)
        return a
    try:
        import scanpy.external as sce
        sce.pp.harmony_integrate(a, key=batch_key)
        sc.pp.neighbors(a, n_pcs=40, use_rep="X_pca_harmony")
        a.uns["integration"] = "harmony"
    except Exception as e:
        logging.info("Harmony unavailable (%s); using uncorrected PCA", e)
        sc.pp.neighbors(a, n_pcs=40)
        a.uns["integration"] = "pca"
    return a


def cluster_and_embed(a, resolution: float = 0.6):
    import scanpy as sc
    sc.tl.umap(a)
    sc.tl.leiden(a, resolution=resolution, flavor="igraph", n_iterations=2,
                  directed=False)
    return a


def annotate_with_markers(a, tissue: str):
    """Score each Leiden cluster against the configured marker panel and
    assign the best-scoring cell-type label."""
    import scanpy as sc
    import numpy as np
    panel = MOUSE_MARKERS.get(tissue, {})
    if not panel:
        a.obs["cell_type_marker"] = a.obs["leiden"].astype(str)
        return a
    scores = {}
    # After QC we restrict a.X to HVGs but keep all genes in a.raw, so score
    # markers against raw (marker genes are frequently not among the HVGs).
    use_raw = a.raw is not None
    marker_vocab = set(a.raw.var_names if use_raw else a.var_names)
    for ct, genes in panel.items():
        present = [g for g in genes if g in marker_vocab]
        if not present:
            continue
        sc.tl.score_genes(a, present, score_name=f"_score_{ct}", use_raw=use_raw)
        scores[ct] = a.obs[f"_score_{ct}"].to_numpy()
    if not scores:
        a.obs["cell_type_marker"] = a.obs["leiden"].astype(str)
        return a
    cell_types = list(scores.keys())
    score_mat = np.stack([scores[c] for c in cell_types], axis=1)
    cluster_call = {}
    for cl, idx in a.obs.groupby("leiden", observed=True).indices.items():
        cluster_means = score_mat[idx].mean(axis=0)
        cluster_call[cl] = cell_types[int(cluster_means.argmax())]
    a.obs["cell_type_marker"] = a.obs["leiden"].map(cluster_call).astype("category")
    return a


def cmd_analyze(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    import anndata as ad
    a = ad.read_h5ad(args.input)
    qc = DEFAULT_QC.copy()
    qc.update({k: getattr(args, k) for k in DEFAULT_QC
               if getattr(args, k, None) is not None})
    a = run_qc_normalize(a, qc=qc)
    if getattr(args, "detect_doublets", False):
        info = _run_scrublet(a)
        if info:
            n_dbl = int(sum(info["predicted_doublet"]))
            logging.info("Scrublet: %d / %d (%.1f%%) cells flagged as doublets",
                          n_dbl, a.n_obs, 100 * n_dbl / max(a.n_obs, 1))
    bk = args.batch_key if args.batch_key in a.obs.columns else "sample_id"
    a = integrate_pools(a, batch_key=bk)
    a = cluster_and_embed(a, resolution=args.resolution)
    if "cell_name" in a.obs.columns:
        a.obs["cell_type"] = a.obs["cell_name"]
    elif args.tissue:
        a = annotate_with_markers(a, args.tissue)
        a.obs["cell_type"] = a.obs["cell_type_marker"]
    else:
        a.obs["cell_type"] = a.obs["leiden"].astype(str)
    out = H5AD_CACHE_DIR / f"{safe_label(args.label or 'splitseq_analysis')}_processed.h5ad"
    a.write_h5ad(out)
    print(f"Wrote analyzed AnnData: {out}")
    print(f"Cells: {a.n_obs}, genes: {a.n_vars}, "
          f"cell types: {a.obs['cell_type'].nunique()}")
    return out


# ------------------------------- Plotting ------------------------------------

def _setup_pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300,
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.linewidth": 0.8, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


# ---------------------------------------------------------------------------
# SPLiT-seq QC helpers — knee plot, pre/post violins, per-Rd1 well distribution,
# optional doublet detection. Closes the HIGH-priority gaps surfaced by the
# Chipeyown SDAT audit: a canonical SPLiT-seq toolkit ships these as
# minimum-viable QC; the IGVFagent skill previously did not.
# ---------------------------------------------------------------------------

_SPLITSEQ_BARCODE_DIR = Path(__file__).resolve().parents[1] / \
                             "Data" / "Reference" / "SPLiTseq"


def _load_splitseq_barcodes(rd: int) -> "dict[str, str]":
    """Load the vendored SPLiT-seq Rd{1,2,3} whitelist as {seq → well_id}.

    Source: Chipeyown SDAT toolkit (MIT) — see
    ``Data/Reference/SPLiTseq/README.md`` for attribution.
    """
    path = _SPLITSEQ_BARCODE_DIR / f"barcodes_rd{rd}.tsv"
    if not path.exists():
        return {}
    out: "dict[str, str]" = {}
    with path.open() as f:
        next(f)  # header
        for line in f:
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) >= 4:
                out[cols[3].upper()] = cols[0]
    return out


def _knee_plot(ax, total_counts):
    """Barcode-rank vs total UMI counts (Drop-seq / 10x knee plot)."""
    import numpy as np
    counts = np.sort(np.asarray(total_counts, dtype=float))[::-1]
    if len(counts) == 0:
        ax.text(0.5, 0.5, "(no cells)", ha="center", va="center",
                  transform=ax.transAxes)
        return
    ax.loglog(range(1, len(counts) + 1), counts, "-", color="#1F77B4", lw=1.2)
    ax.set_xlabel("Cell barcode rank (log)")
    ax.set_ylabel("Total UMI counts (log)")
    ax.set_title("Knee plot — barcode rank vs total UMI",
                  loc="left", weight="bold", fontsize=10)
    ax.grid(alpha=0.25, linestyle=":")


def _pre_post_violin(ax, pre, post, label, threshold=None):
    """Side-by-side violin: distribution before vs after QC filtering."""
    import numpy as np
    parts = ax.violinplot([pre, post], positions=[0, 1], widths=0.8,
                            showmedians=True)
    for body, color in zip(parts["bodies"], ("#FCB16E", "#3498DB")):
        body.set_facecolor(color); body.set_alpha(0.7)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["pre-filter", "post-filter"],
                                                 fontsize=8)
    ax.set_ylabel(label, fontsize=9)
    if threshold is not None:
        ax.axhline(threshold, color="red", ls=":", lw=0.6,
                     label=f"threshold = {threshold:g}")
        ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.25, axis="y", linestyle=":")


def _rd1_well_distribution(a):
    """Tally cells per Rd1 well by matching the first 8 bp of each cell
    barcode against the vendored Rd1 whitelist.

    Returns a dict {well_id → n_cells} for matched cells, plus a
    metadata dict with no-match count etc.
    """
    rd1_map = _load_splitseq_barcodes(1)
    if not rd1_map:
        return {}, {"error": "Rd1 whitelist not available"}
    matched: "dict[str, int]" = {}
    n_match = 0
    n_total = a.n_obs
    for bc in a.obs_names:
        # Cell barcodes in IGVF SPLiT-seq matrices typically look like
        # "<Rd1><Rd2><Rd3>" (24 bp) or with dashes. Try a few shapes.
        token = (str(bc).replace("-", "").replace("_", "")[:8]).upper()
        if token in rd1_map:
            well = rd1_map[token]
            matched[well] = matched.get(well, 0) + 1
            n_match += 1
    return matched, {"n_total": n_total, "n_matched": n_match,
                       "n_unmatched": n_total - n_match}


def _per_well_heatmap(ax, matched):
    """8 × 12 plate heatmap of cell counts per Rd1 well."""
    import numpy as np
    rows = "ABCDEFGH"
    cols = list(range(1, 13))
    M = np.zeros((8, 12), dtype=float)
    for well, n in matched.items():
        if not well or len(well) < 2:
            continue
        r = well[0]
        if r not in rows:
            continue
        try:
            c = int(well[1:])
        except ValueError:
            continue
        if 1 <= c <= 12:
            M[rows.index(r), c - 1] = n
    im = ax.imshow(M, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(12)); ax.set_xticklabels(cols, fontsize=7)
    ax.set_yticks(range(8));  ax.set_yticklabels(list(rows), fontsize=7)
    ax.set_title("Cells per Rd1 well (8 × 12 plate)",
                  loc="left", weight="bold", fontsize=10)
    for i in range(8):
        for j in range(12):
            if M[i, j] > 0:
                ax.text(j, i, f"{int(M[i, j]):,}", ha="center", va="center",
                          fontsize=6,
                          color="white" if M[i, j] > M.max() / 2 else "black")
    return im


def _run_scrublet(a) -> "Optional[dict]":
    """Run Scrublet doublet detection; soft-skip if scrublet not installed.

    Returns ``{'doublet_score': np.ndarray, 'predicted_doublet': np.ndarray}``
    or None if Scrublet is unavailable / matrix is too small.
    """
    try:
        import scrublet as scr  # type: ignore
    except Exception:
        try:
            # Newer scanpy ships sc.pp.scrublet
            import scanpy as sc  # type: ignore
            sc.pp.scrublet(a, batch_key=("dataset"
                                            if "dataset" in a.obs.columns
                                            else None))
            return {
                "doublet_score":     a.obs["doublet_score"].values.copy(),
                "predicted_doublet": a.obs["predicted_doublet"].values.copy(),
            }
        except Exception as e:
            logging.info("Doublet detection skipped: %s", e)
            return None
    try:
        import numpy as np
        X = (a.layers["counts"] if "counts" in a.layers else a.X)
        scrub = scr.Scrublet(X)
        scores, calls = scrub.scrub_doublets(verbose=False)
        a.obs["doublet_score"] = scores
        a.obs["predicted_doublet"] = calls
        return {"doublet_score": scores, "predicted_doublet": calls}
    except Exception as e:
        logging.warning("Scrublet failed: %s", e)
        return None


def cmd_plot(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.gridspec import GridSpec
    _setup_pub_style()
    a = ad.read_h5ad(args.input)
    ts = timestamp()
    label = safe_label(args.label or "splitseq")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / "Plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    cell_types = list(a.obs["cell_type"].astype(str).unique())
    palette = make_palette(cell_types)

    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.45,
                  top=0.93, bottom=0.07, left=0.05, right=0.95)

    # a: UMAP by cell type
    ax = fig.add_subplot(gs[0, 0])
    _umap(ax, a, color_by="cell_type", palette=palette,
          title="a  UMAP by cell type", legend=True)

    # b: UMAP by batch / pool
    ax = fig.add_subplot(gs[0, 1])
    if "dataset" in a.obs.columns:
        bk = "dataset"
    elif "sample_id" in a.obs.columns:
        bk = "sample_id"
    else:
        bk = "leiden"
    _umap(ax, a, color_by=bk, title=f"b  UMAP by {bk}",
          legend=False, alpha=0.6)

    # c: UMAP by strain (if available)
    ax = fig.add_subplot(gs[0, 2])
    strain_key = next((k for k in ("strain", "strain_background", "DonorStrain")
                       if k in a.obs.columns), None)
    if strain_key:
        _umap(ax, a, color_by=strain_key,
              palette={s: FOUNDER_COLORS.get(s, "#999")
                       for s in a.obs[strain_key].astype(str).unique()},
              title=f"c  UMAP by {strain_key}", legend=True)
    else:
        ax.text(0.5, 0.5, "No strain column found\n"
                          "(populate `strain` after demux)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.axis("off"); ax.set_title("c  Strain UMAP", loc="left", weight="bold")

    # d: marker dot plot
    ax = fig.add_subplot(gs[1, 0:2])
    panel = MOUSE_MARKERS.get(args.tissue or "", {})
    if panel:
        _dotplot(ax, a, panel, palette=palette)
        ax.set_title("d  Marker expression", loc="left", weight="bold")
    else:
        ax.text(0.5, 0.5, "Pass --tissue (gonad / adrenal / brain / "
                          "liver / heart / kidney / muscle) to render the "
                          "marker dot plot.", ha="center", va="center",
                transform=ax.transAxes, fontsize=9)
        ax.axis("off")

    # e: composition stacked bar by batch / strain
    ax = fig.add_subplot(gs[1, 2])
    by = strain_key or bk
    _composition_bar(ax, a, by=by, palette=palette,
                      title=f"e  Cell-type composition by {by}")

    fig.suptitle("SPLiT-seq pipeline overview", fontsize=12, weight="bold",
                  y=0.97)
    main_path = plot_dir / f"{ts}_{label}_main.png"
    fig.savefig(main_path, bbox_inches="tight", facecolor="white")
    fig.savefig(main_path.with_suffix(".svg"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # ---------- QC panel: knee + pre/post violins + per-Rd1 well -----------
    # Adds the Chipeyown-style QC outputs the audit flagged as HIGH gaps:
    # universal knee plot, pre-vs-post-filter violins, per-Rd1-well sample
    # distribution heatmap, optional doublet detection.
    qc_path = plot_dir / f"{ts}_{label}_qc.png"
    qc_pre = a.uns.get("qc_pre", {}) if hasattr(a, "uns") else {}
    fig2 = plt.figure(figsize=(15, 9))
    gs2 = GridSpec(2, 3, figure=fig2, hspace=0.40, wspace=0.45,
                     top=0.92, bottom=0.08, left=0.06, right=0.96)

    # (a) Knee plot
    ax = fig2.add_subplot(gs2[0, 0])
    _knee_plot(ax, np.asarray(a.obs.get("total_counts",
                                            a.X.sum(axis=1).A1
                                            if hasattr(a.X, "A1")
                                            else np.asarray(a.X.sum(axis=1)).ravel())))

    # (b) (c) (d): pre vs post QC violins for the 3 canonical metrics
    metrics = [
        ("n_genes",       "Genes / cell",      "min_genes"),
        ("total_counts",  "Total UMI / cell",  "min_counts"),
        ("pct_counts_mt", "% mito",            "max_pct_mito"),
    ]
    for i, (key, lab, thr_key) in enumerate(metrics):
        ax = fig2.add_subplot(gs2[0, i] if i > 0 else gs2[0, i])
        if i == 0:
            continue
        pre = (qc_pre.get(key, []) if isinstance(qc_pre, dict)
                else []) or []
        post_col = ("n_genes_by_counts" if key == "n_genes" else key)
        post = (np.asarray(a.obs[post_col], dtype=float)
                 if post_col in a.obs.columns else [])
        if len(pre) and len(post):
            thr = qc_pre.get("thresholds", {}).get(thr_key)
            _pre_post_violin(ax, pre, post, lab, threshold=thr)
            ax.set_title(f"{chr(0x62 + i)}  {lab} — pre vs post filter",
                          loc="left", weight="bold", fontsize=10)
        else:
            ax.text(0.5, 0.5,
                     f"Run `process-local` to populate `qc_pre`\n"
                     f"for the pre-vs-post overlay.",
                     ha="center", va="center", transform=ax.transAxes,
                     fontsize=8)
            ax.axis("off")

    # (e) Per-Rd1-well plate heatmap (SPLiT-seq native sample multiplexing)
    ax = fig2.add_subplot(gs2[1, 0:2])
    matched, info = _rd1_well_distribution(a)
    if matched:
        im = _per_well_heatmap(ax, matched)
        fig2.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="cells")
        ax.set_title(
            f"e  Cells per Rd1 well "
            f"({info['n_matched']:,}/{info['n_total']:,} cell-barcodes "
            f"matched the 96-well whitelist)",
            loc="left", weight="bold", fontsize=10,
        )
        # Write a TSV the user can join with their sample sheet
        well_tsv = out_dir / f"{ts}_{label}_rd1_wells.tsv"
        with well_tsv.open("w") as f:
            f.write("well_id\tn_cells\n")
            for w in sorted(matched, key=lambda w: (w[0], int(w[1:]))):
                f.write(f"{w}\t{matched[w]}\n")
        logging.info("Per-well counts: %s", well_tsv)
    else:
        ax.text(0.5, 0.5,
                  f"Rd1 well distribution unavailable\n"
                  f"({info.get('error', 'no cell-barcode prefix match')}).\n\n"
                  f"Vendor the Chipeyown SDAT 96×96×96 barcodes at\n"
                  f"Data/Reference/SPLiTseq/barcodes_rd1.tsv to enable.",
                  ha="center", va="center", transform=ax.transAxes,
                  fontsize=8)
        ax.axis("off")

    # (f) Doublet panel — optional Scrublet
    ax = fig2.add_subplot(gs2[1, 2])
    if "predicted_doublet" in a.obs.columns:
        n_dbl = int((a.obs["predicted_doublet"].astype(bool)).sum())
        pct = 100 * n_dbl / max(a.n_obs, 1)
        ax.hist(np.asarray(a.obs["doublet_score"], dtype=float),
                  bins=40, color="#7570B3", alpha=0.85,
                  edgecolor="white")
        ax.set_xlabel("Doublet score (Scrublet)")
        ax.set_ylabel("Cells")
        ax.set_title(f"f  Doublet detection — "
                       f"{n_dbl:,} called ({pct:.1f}%)",
                       loc="left", weight="bold", fontsize=10)
        ax.grid(alpha=0.25, linestyle=":")
    else:
        ax.text(0.5, 0.5,
                  "Run `sc.pp.scrublet(adata)` or install scrublet\n"
                  "to populate predicted_doublet / doublet_score.",
                  ha="center", va="center", transform=ax.transAxes,
                  fontsize=8)
        ax.axis("off")

    fig2.suptitle("SPLiT-seq QC panel", fontsize=12, weight="bold", y=0.96)
    fig2.savefig(qc_path, bbox_inches="tight", facecolor="white")
    fig2.savefig(qc_path.with_suffix(".svg"), bbox_inches="tight",
                   facecolor="white")
    plt.close(fig2)

    # write report
    report = out_dir / f"{ts}_{label}_report.md"
    counts = a.obs["cell_type"].value_counts()
    lines = [f"# SPLiT-seq Run Report: `{label}`", "",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "",
             f"Cells: **{a.n_obs:,}**, genes: **{a.n_vars:,}**",
             f"Cell types: {len(counts)}",
             f"Integration: {a.uns.get('integration','n/a')}",
             "",
             "## Cell-type counts", "",
             "| Cell type | n cells | % |",
             "|---|---:|---:|"]
    for ct, n in counts.items():
        lines.append(f"| {ct} | {n:,} | {n / a.n_obs * 100:.2f} |")
    lines += ["", f"![Main figure]({main_path.relative_to(out_dir)})", "",
                f"![QC panel]({qc_path.relative_to(out_dir)})", ""]
    if matched:
        lines += ["## Per-Rd1-well sample distribution",
                    f"_(see `{(out_dir / (ts + '_' + label + '_rd1_wells.tsv')).name}` "
                    f"for the full TSV; matched {info['n_matched']:,} of "
                    f"{info['n_total']:,} cells to the 96-well whitelist)_", ""]
    report.write_text("\n".join(lines))
    print(f"Plot:    {main_path}")
    print(f"QC:      {qc_path}")
    print(f"Report:  {report}")
    return main_path


def _umap(ax, a, color_by, palette=None, title="", legend=True,
            alpha=0.85, s=2.0):
    import numpy as np
    xy = a.obsm.get("X_umap")
    if xy is None:
        ax.axis("off"); return
    vals = a.obs[color_by].astype(str).to_numpy()
    cats = sorted(set(vals))
    if palette is None:
        palette = make_palette(cats)
    for c in cats:
        m = vals == c
        ax.scatter(xy[m, 0], xy[m, 1], s=s, alpha=alpha,
                    c=palette.get(c, "#888"), linewidths=0,
                    label=c, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(title, loc="left", weight="bold")
    if legend and len(cats) <= 20:
        leg = ax.legend(markerscale=2.5, loc="center left",
                         bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
        for h in leg.legend_handles:
            h.set_alpha(1.0)


def _dotplot(ax, a, panel, palette=None, size_max=120):
    import matplotlib as mpl
    import numpy as np
    flat = []
    breaks = []
    for grp, genes in panel.items():
        for g in genes:
            if g in a.var_names:
                flat.append(g)
        breaks.append((grp, len(flat)))
    if not flat:
        ax.text(0.5, 0.5, "No panel markers found in matrix",
                ha="center", va="center", transform=ax.transAxes); return
    types = sorted([t for t in a.obs["cell_type"].astype(str).unique()])
    means = np.zeros((len(types), len(flat)))
    pcts = np.zeros((len(types), len(flat)))
    var_idx = {g: i for i, g in enumerate(a.var_names)}
    X = a.X
    for i, t in enumerate(types):
        mask = (a.obs["cell_type"].astype(str) == t).to_numpy()
        if not mask.any():
            continue
        sub = X[mask]
        for j, g in enumerate(flat):
            col = sub[:, var_idx[g]]
            arr = col.toarray().ravel() if hasattr(col, "toarray") else \
                   np.asarray(col).ravel()
            means[i, j] = float(arr.mean())
            pcts[i, j] = float((arr > 0).mean() * 100)
    norm = mpl.colors.Normalize(vmin=0,
                                  vmax=float(np.nanpercentile(means, 99) + 1e-6))
    cmap = mpl.colormaps["Reds"]
    for i in range(len(types)):
        for j in range(len(flat)):
            ax.scatter(j, i, s=pcts[i, j] / 100 * size_max + 4,
                        c=[cmap(norm(means[i, j]))],
                        edgecolors="black", linewidths=0.3, zorder=3)
    ax.set_xticks(np.arange(len(flat)))
    ax.set_xticklabels(flat, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(types))); ax.set_yticklabels(types, fontsize=8)
    ax.set_xlim(-0.6, len(flat) - 0.4); ax.set_ylim(len(types) - 0.4, -0.6)
    ax.tick_params(axis="both", which="both", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = 0
    trans = ax.get_xaxis_transform()
    for grp, end in breaks:
        if end == prev:
            continue
        if prev > 0:
            ax.axvline(prev - 0.5, color="#cccccc", lw=0.5, zorder=1)
        ax.plot([prev - 0.4, end - 0.6], [1.02, 1.02],
                 color="black", lw=0.8, transform=trans, clip_on=False)
        ax.text((prev + end - 1) / 2, 1.05, grp, ha="center", va="bottom",
                fontsize=7, fontweight="bold", transform=trans, clip_on=False)
        prev = end
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = ax.figure.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Mean log-expression", fontsize=7)
    cb.ax.tick_params(labelsize=6)


def _composition_bar(ax, a, by, palette, title):
    import numpy as np
    import pandas as pd
    df = a.obs[[by, "cell_type"]].astype(str)
    ct = pd.crosstab(df[by], df["cell_type"])
    frac = ct.div(ct.sum(axis=1), axis=0)
    bottom = np.zeros(frac.shape[0])
    x = np.arange(frac.shape[0])
    for c in frac.columns:
        ax.bar(x, frac[c], bottom=bottom, color=palette.get(c, "#888"),
                width=0.8, linewidth=0)
        bottom += frac[c].to_numpy()
    ax.set_xticks(x)
    ax.set_xticklabels(frac.index, rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0, 1); ax.set_ylabel("Fraction"); ax.set_xlabel(by)
    ax.set_title(title, loc="left", weight="bold")


def make_palette(categories: list[str]) -> dict[str, str]:
    base = ["#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
            "#CC79A7", "#7B5BA6", "#F0E442", "#999999", "#A14586",
            "#73C03E", "#B5651D", "#8DA0CB", "#FFD92F", "#444444",
            "#4E79A7", "#F28E2B", "#76B7B2", "#59A14F", "#EDC948"]
    return {c: base[i % len(base)] for i, c in enumerate(sorted(categories))}


# ------------------------------- Strain comparison ---------------------------

def cmd_compare_strains(args: argparse.Namespace) -> Path:
    """Per cell type, run a Wilcoxon DEG between a target strain and the rest."""
    setup_logging(); mkdirs()
    import anndata as ad
    import scanpy as sc
    a = ad.read_h5ad(args.input)
    strain_key = next((k for k in ("strain", "strain_background", "DonorStrain")
                       if k in a.obs.columns), None)
    if not strain_key:
        raise SystemExit("No strain column on the AnnData. Run demultiplexing "
                         "and add a `strain` obs column first.")
    ts = timestamp()
    label = safe_label(args.label or f"strain_compare_{strain_key}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict] = []
    for ct in sorted(a.obs["cell_type"].astype(str).unique()):
        sub = a[a.obs["cell_type"].astype(str) == ct].copy()
        if sub.obs[strain_key].nunique() < 2 or sub.n_obs < args.min_cells:
            continue
        sc.tl.rank_genes_groups(sub, strain_key, method="wilcoxon",
                                  n_genes=200)
        df = sc.get.rank_genes_groups_df(sub, group=None)
        df["cell_type"] = ct
        rows = df[df["pvals_adj"] < args.padj_cutoff].sort_values("logfoldchanges",
                                                                    key=abs,
                                                                    ascending=False)
        rows.to_csv(out_dir / f"DEG_{safe_label(ct)}.csv", index=False)
        summary_rows.append({
            "cell_type": ct,
            "n_cells": sub.n_obs,
            "n_strains": sub.obs[strain_key].nunique(),
            "n_DEGs": len(rows),
        })
    write_csv(out_dir / "compare_strains_summary.csv", summary_rows)
    print(f"Wrote per-cell-type DEG tables in {out_dir}")
    return out_dir


# ------------------------------- Demultiplex script --------------------------

def cmd_demux_script(args: argparse.Namespace) -> Path:
    """Emit a runnable shell script for genotype demultiplexing.

    The skill does not run the demultiplexer (souporcell or vireo) itself —
    these need a BAM, a reference genome, and per-donor genotype VCFs that
    live outside the IGVF agent's scope. The emitted script wires up the
    common-case command lines so the user can run them locally.
    """
    setup_logging(); mkdirs()
    ts = timestamp()
    out = REPORT_DIR / f"{ts}_demux_{safe_label(args.tool)}.sh"
    bam = args.bam or "<path-to-cellranger-bam>"
    ref = args.reference or "<path-to-genome.fa>"
    bcs = args.barcodes or "<path-to-filtered-barcodes.tsv>"
    vcf = args.vcf or "<path-to-donor-genotypes.vcf.gz>"
    n_donors = args.n_donors or 8
    if args.tool == "souporcell":
        body = f"""#!/usr/bin/env bash
# Souporcell genotype-free donor demultiplexing for SPLiT-seq.
# https://github.com/wheaton5/souporcell
set -euo pipefail
BAM={bam}
REF={ref}
BCS={bcs}
N_DONORS={n_donors}
OUT_DIR=./souporcell_output
mkdir -p "$OUT_DIR"
souporcell_pipeline.py \\
    --bam "$BAM" \\
    --barcodes "$BCS" \\
    --fasta "$REF" \\
    --threads 16 \\
    --out_dir "$OUT_DIR" \\
    --clusters $N_DONORS
echo "Souporcell finished. Use clusters.tsv to assign donor IDs."
"""
    else:
        body = f"""#!/usr/bin/env bash
# Vireo genotype-aware donor demultiplexing for SPLiT-seq.
# https://github.com/single-cell-genetics/vireo
set -euo pipefail
BAM={bam}
BCS={bcs}
VCF={vcf}
N_DONORS={n_donors}
CELLSNP_OUT=./cellsnp_output
VIREO_OUT=./vireo_output
mkdir -p "$CELLSNP_OUT" "$VIREO_OUT"
cellsnp-lite -s "$BAM" -b "$BCS" -O "$CELLSNP_OUT" \\
    -R "$VCF" -p 16 --minMAF 0.1 --minCOUNT 20 --gzip
vireo -c "$CELLSNP_OUT" -d "$VCF" -o "$VIREO_OUT" -N $N_DONORS \\
    --extraDonor=0
echo "Vireo finished. Use donor_ids.tsv to assign donor IDs."
"""
    out.write_text(body)
    out.chmod(0o755)
    print(f"Wrote: {out}")
    print("Edit the variables at the top before running.")
    return out


# ------------------------------- Playbook ------------------------------------

def cmd_write_playbook(_args) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "SPLITSEQ_SKILLS.md"
    lines = [
        "# Skill: Parse SPLiT-seq pipeline (IGVF)",
        "",
        "Reads, processes, and analyzes Parse-Biosciences-style "
        "combinatorial-barcoding single-nucleus RNA-seq datasets from the "
        "IGVF Portal. Designed around SPLiT-seq's structural quirks: "
        "multiplexed sub-pools, tens of donors per pool, mtx tarballs, "
        "and the analytical value living in **per-strain / per-donor / "
        "per-pool** comparisons.",
        "",
        "## Subcommands",
        "",
        "### `retrieve` — discover SPLiT-seq AnalysisSets",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py retrieve \\",
        "    --limit 100 --fetch-file-details --label mortazavi_8cube",
        "python3 Scripts/splitseq_pipeline.py retrieve \\",
        "    --limit 50 --lab 'Ali Mortazavi, UCI' --taxa 'Mus musculus'",
        "python3 Scripts/splitseq_pipeline.py retrieve \\",
        "    --limit 20 --sample-type gonad",
        "```",
        "",
        "### `manifest` — per-pool / per-donor manifest for one or more accessions",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py manifest \\",
        "    --accessions IGVFDS3222WCZH,IGVFDS6290WNNH \\",
        "    --label 8cube_mouse_demo",
        "```",
        "",
        "### `download` — pull files for a manifest under a size cap",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py download \\",
        "    --manifest Data/Manifests/SPLiTseq/<files_per_pool>.csv \\",
        "    --max-download-gb 5 \\",
        "    --only 'sparse gene count matrix' 'cell annotations'",
        "```",
        "",
        "### `process-local` — load downloaded mtx bundles into AnnData",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py process-local \\",
        "    --label 8cube_mouse_demo",
        "```",
        "",
        "Concatenates every downloaded AnalysisSet under "
        "`Data/IGVF/SPLiTseq/Downloads/` into a single `.h5ad` under "
        "`Data/Cache/SPLiTseq/`. Each cell carries a `dataset` and "
        "`sample_id` column for downstream batch correction.",
        "",
        "### `analyze` — full pipeline (QC → norm → integrate → UMAP → cluster → label)",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py analyze \\",
        "    --input Data/Cache/SPLiTseq/8cube_mouse_demo.h5ad \\",
        "    --tissue gonad \\",
        "    --resolution 0.6 \\",
        "    --label 8cube_mouse_demo",
        "```",
        "",
        "Tunable QC thresholds via `--min-genes`, `--max-genes`, `--min-counts`, "
        "`--max-counts`, `--max-pct-mito`. Pool integration uses Harmony if "
        "`scanpy.external` is available, else falls back to uncorrected PCA.",
        "",
        "Auto-annotation panels currently shipped: `gonad`, `adrenal`, "
        "`brain`, `liver`, `heart`, `kidney`, `muscle`. Pass `--tissue` to "
        "use one. If the AnnData already carries a `cell_name` obs column "
        "from an IGVF cell-annotation TSV (as in some principal "
        "AnalysisSets), that label is used directly.",
        "",
        "### `plot` — publication-style multi-panel figure",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py plot \\",
        "    --input Data/Cache/SPLiTseq/8cube_mouse_demo_processed.h5ad \\",
        "    --tissue gonad --label 8cube_mouse_demo",
        "```",
        "",
        "Panels: UMAP by cell type, UMAP by pool/batch, UMAP by strain "
        "(if the `strain` obs column is populated), marker dot plot, "
        "stacked composition bar.",
        "",
        "### `compare-strains` — per-cell-type DEG across the 8 founders",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py compare-strains \\",
        "    --input Data/Cache/SPLiTseq/8cube_mouse_demo_processed.h5ad \\",
        "    --label 8cube_strain_DEG --padj-cutoff 0.05 --min-cells 50",
        "```",
        "",
        "Requires a `strain` (or `strain_background` / `DonorStrain`) "
        "column on `obs`. Populate it after running donor demultiplexing.",
        "",
        "### `demux-script` — emit a runnable demux pipeline",
        "",
        "```bash",
        "python3 Scripts/splitseq_pipeline.py demux-script --tool souporcell --n-donors 8",
        "python3 Scripts/splitseq_pipeline.py demux-script --tool vireo --n-donors 8",
        "```",
        "",
        "The skill does NOT run the demultiplexer — these tools need BAMs "
        "and (for vireo) per-donor VCFs that live outside the agent. The "
        "emitted script wires up the common-case command lines so you can "
        "run them locally and feed the resulting `donor_ids.tsv` back as a "
        "`strain` obs column.",
        "",
        "## Outputs",
        "",
        "- `Data/IGVF/SPLiTseq/Metadata/<accession>.json` — hydrated portal JSON",
        "- `Data/IGVF/SPLiTseq/Downloads/<accession>/` — per-AnalysisSet bundles",
        "- `Data/Manifests/SPLiTseq/` — dataset + per-pool manifests, download log",
        "- `Data/Cache/SPLiTseq/<label>.h5ad` — processed AnnData",
        "- `Docs/SPLiTseq/<timestamp>_<label>/` — reports + plots",
        "",
        "## How this chains with other skills",
        "",
        "- After `retrieve`, hand the manifest to `download` and then "
        "`process-local`.",
        "- After `analyze` produces a cell-type-labeled `.h5ad`, hand a "
        "gene list (e.g. top markers) to `Scripts/reference_skill.py "
        "validate` to check prior literature.",
        "- For workflow scaffolding before retrieval, run "
        "`Scripts/reference_skill.py design --data-type parse_split_seq` "
        "to see the curated workflow plus cognate published studies.",
        "- The bundled mouse marker panels (gonad / adrenal / brain / "
        "liver / heart / kidney / muscle) match the tissues covered by "
        "the Mortazavi 8-cube founder atlas, so a one-shot run through "
        "`retrieve → download → process-local → analyze → plot` produces "
        "publication-style figures with no manual annotation.",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI ---------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="SPLiT-seq pipeline for IGVF Portal data."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("retrieve",
                        help="Search SPLiT-seq AnalysisSets and write a manifest.")
    s.add_argument("--limit", type=int, default=50)
    s.add_argument("--fetch-file-details", action="store_true")
    s.add_argument("--lab", default=None)
    s.add_argument("--sample-type", default=None,
                    help="biosample term name (e.g. 'gonad', 'adrenal gland').")
    s.add_argument("--status", default="released")
    s.add_argument("--taxa", default=None,
                    help="e.g. 'Mus musculus' or 'Homo sapiens'.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_retrieve)

    s = sub.add_parser("manifest",
                        help="Hydrate metadata for one or more accessions.")
    s.add_argument("--accessions", required=True,
                    help="Comma-separated AnalysisSet accessions.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_manifest)

    s = sub.add_parser("download",
                        help="Download files referenced by a manifest.")
    s.add_argument("--manifest", required=True)
    s.add_argument("--max-download-gb", type=float, default=5.0)
    s.add_argument("--only", nargs="*", default=None,
                    help="Filter by content_type (e.g. 'sparse gene count matrix').")
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("process-local",
                        help="Load downloaded mtx bundles into AnnData.")
    s.add_argument("--accessions", default=None,
                    help="Comma-separated whitelist; default = all downloaded.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_process_local)

    s = sub.add_parser("analyze",
                        help="QC -> normalize -> integrate -> UMAP -> cluster -> annotate.")
    s.add_argument("--input", required=True)
    s.add_argument("--tissue",
                    choices=sorted(MOUSE_MARKERS.keys()),
                    default=None)
    s.add_argument("--resolution", type=float, default=0.6)
    s.add_argument("--batch-key", default="dataset")
    s.add_argument("--min-genes", type=int, default=None)
    s.add_argument("--max-genes", type=int, default=None)
    s.add_argument("--min-counts", type=int, default=None)
    s.add_argument("--max-counts", type=int, default=None)
    s.add_argument("--max-pct-mito", type=float, default=None)
    s.add_argument("--min-cells", type=int, default=None)
    s.add_argument("--detect-doublets", action="store_true",
                    help="Run Scrublet (or sc.pp.scrublet) and annotate "
                         "predicted_doublet / doublet_score on obs.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_analyze)

    s = sub.add_parser("plot",
                        help="Publication-style multi-panel figure.")
    s.add_argument("--input", required=True)
    s.add_argument("--tissue", choices=sorted(MOUSE_MARKERS.keys()),
                    default=None)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_plot)

    s = sub.add_parser("compare-strains",
                        help="Per-cell-type strain DEG (Wilcoxon).")
    s.add_argument("--input", required=True)
    s.add_argument("--padj-cutoff", type=float, default=0.05)
    s.add_argument("--min-cells", type=int, default=50)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_compare_strains)

    s = sub.add_parser("demux-script",
                        help="Emit a souporcell or vireo runnable script.")
    s.add_argument("--tool", choices=["souporcell", "vireo"], default="souporcell")
    s.add_argument("--n-donors", type=int, default=8)
    s.add_argument("--bam", default=None)
    s.add_argument("--reference", default=None)
    s.add_argument("--barcodes", default=None)
    s.add_argument("--vcf", default=None)
    s.set_defaults(func=cmd_demux_script)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/SPLITSEQ_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
