#!/usr/bin/env python3
"""Download and analyze ten scRNA-seq files from the IGVF Portal.

Selection:
  - All 6 publicly released "gene quantifications" TSV files
    (content_type=gene quantifications)
  - 4 smallest released "sparse gene count matrix" h5ad files
    (file_format=h5ad)

Outputs are written under repository-relative paths:
  Data/SingleCell/Portal10/<accession>.<ext>            raw downloads
  Data/Manifests/SingleCell/portal10_selection.csv      file manifest
  Data/Manifests/SingleCell/portal10_download_log.csv   download log
  Docs/SingleCell/Plots/Portal10/*.png                  per-file plots
  Docs/SingleCell/Portal10_report.md                    summary report
  Docs/Logs/portal_scrna_10_*.log                       runtime log
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
DOWNLOAD_DIR = DATA_DIR / "SingleCell" / "Portal10"
MANIFEST_DIR = DATA_DIR / "Manifests" / "SingleCell"
PLOT_DIR = DOCS_DIR / "SingleCell" / "Plots" / "Portal10"
REPORT_PATH = DOCS_DIR / "SingleCell" / "Portal10_report.md"
SELECTION_CSV = MANIFEST_DIR / "portal10_selection.csv"
DOWNLOAD_LOG_CSV = MANIFEST_DIR / "portal10_download_log.csv"
ANALYSIS_JSON = DATA_DIR / "SingleCell" / "Portal10" / "analysis_summary.json"

PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://api.data.igvf.org").rstrip("/")
USER_AGENT = "IGVFagent-portal-scrna-10/0.1"


_PALETTE = [
    "#2f6f9f", "#c44e52", "#55a868", "#8172b2", "#ccb974", "#64b5cd",
    "#8c613c", "#dc7ec0", "#4c72b0", "#dd8452", "#937860", "#6acc64",
]


def palette(index: int) -> str:
    return _PALETTE[index % len(_PALETTE)]


def rel_path(path: Path | str) -> str:
    """Return path as POSIX-style string relative to repo root.

    Stored in committed outputs to keep the repo location private and portable.
    """
    p = Path(path).resolve()
    try:
        return p.relative_to(ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def abs_from_rel(rel: str) -> Path:
    """Resolve a stored repo-relative path back to an absolute Path."""
    p = Path(rel)
    return p if p.is_absolute() else (ROOT / p)
MANIFEST_FIELDS = [
    "accession",
    "kind",
    "file_format",
    "content_type",
    "file_size_bytes",
    "assembly",
    "transcriptome_annotation",
    "lab",
    "creation_timestamp",
    "href",
    "download_url",
    "local_path",
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"portal_scrna_10_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def http_get_json(url: str, timeout: int = 60) -> Any:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    logging.info("GET %s", url)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def search(params: dict[str, Any]) -> list[dict[str, Any]]:
    query = dict(params)
    query.setdefault("format", "json")
    query.setdefault("frame", "object")
    url = f"{PORTAL_BASE}/search/?{urllib.parse.urlencode(query, doseq=True)}"
    payload = http_get_json(url)
    return [row for row in payload.get("@graph", []) if isinstance(row, dict)]


def file_record(row: dict[str, Any], kind: str) -> dict[str, Any]:
    href = row.get("href") or ""
    download_url = f"{PORTAL_BASE}{href}" if href else ""
    local_abs = (
        DOWNLOAD_DIR
        / f"{row.get('accession','UNKNOWN')}.{(row.get('file_format') or 'bin').replace('h5ad','h5ad')}"
    )
    return {
        "accession": row.get("accession", ""),
        "kind": kind,
        "file_format": row.get("file_format", ""),
        "content_type": row.get("content_type", ""),
        "file_size_bytes": row.get("file_size") or 0,
        "assembly": row.get("assembly", ""),
        "transcriptome_annotation": row.get("transcriptome_annotation", ""),
        "lab": row.get("lab", ""),
        "creation_timestamp": row.get("creation_timestamp", ""),
        "href": href,
        "download_url": download_url,
        "local_path": rel_path(local_abs),
    }


def select_files(h5ad_count: int = 4) -> list[dict[str, Any]]:
    gene_quant_rows = search(
        {"type": "File", "content_type": "gene quantifications", "limit": "100"}
    )
    h5ad_rows = search({"type": "File", "file_format": "h5ad", "limit": "200"})
    h5ad_rows = [r for r in h5ad_rows if r.get("file_size") and r.get("href")]
    h5ad_rows.sort(key=lambda r: r["file_size"])

    records: list[dict[str, Any]] = []
    for row in gene_quant_rows:
        records.append(file_record(row, "gene_quant_tsv"))
    for row in h5ad_rows[:h5ad_count]:
        records.append(file_record(row, "sparse_h5ad"))
    return records


def write_selection(records: list[dict[str, Any]]) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    with SELECTION_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_FIELDS, extrasaction="ignore"
        )
        writer.writeheader()
        for row in records:
            writer.writerow(row)
    logging.info("Wrote selection manifest: %s", SELECTION_CSV)
    return SELECTION_CSV


def read_selection() -> list[dict[str, Any]]:
    if not SELECTION_CSV.exists():
        raise SystemExit(
            f"selection manifest missing: {SELECTION_CSV}. Run 'select' first."
        )
    with SELECTION_CSV.open(encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def download_one(url: str, target: Path, force: bool = False) -> dict[str, Any]:
    target.parent.mkdir(parents=True, exist_ok=True)
    rel = rel_path(target)
    if target.exists() and not force and target.stat().st_size > 0:
        return {"status": "exists", "bytes": target.stat().st_size, "url": url, "path": rel}
    request = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"}
    )
    started = time.time()
    try:
        with urllib.request.urlopen(request, timeout=300) as response, target.open("wb") as fh:
            total = 0
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                total += len(chunk)
        elapsed = time.time() - started
        logging.info(
            "downloaded %s (%.2f MB in %.1fs)",
            target.name,
            total / 1024 / 1024,
            elapsed,
        )
        return {"status": "downloaded", "bytes": total, "url": url, "path": rel}
    except urllib.error.HTTPError as exc:
        return {"status": f"http_error_{exc.code}", "bytes": 0, "url": url, "path": rel}
    except urllib.error.URLError as exc:
        return {"status": f"network_error_{exc.reason}", "bytes": 0, "url": url, "path": rel}


def download_selection(force: bool = False) -> Path:
    records = read_selection()
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for rec in records:
        url = rec["download_url"]
        target = abs_from_rel(rec["local_path"])
        result = download_one(url, target, force=force)
        result["accession"] = rec["accession"]
        result["kind"] = rec["kind"]
        rows.append(result)
        print(f"{rec['accession']:18s} {rec['kind']:18s} {result['status']:>14s} {result['bytes']:>12,} bytes")
    with DOWNLOAD_LOG_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["accession", "kind", "status", "bytes", "url", "path"]
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    logging.info("Wrote download log: %s", DOWNLOAD_LOG_CSV)
    return DOWNLOAD_LOG_CSV


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def _import_analysis_libs():
    import numpy as np  # type: ignore
    import pandas as pd  # type: ignore
    import scipy.sparse as sp  # type: ignore
    import scipy.sparse.linalg as spla  # type: ignore
    import matplotlib  # type: ignore

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    return np, pd, sp, spla, plt


def analyze_starsolo_gene_quant(path, accession, rows, cols, metadata_lines, plt, np):
    """Schema: geneID, unstranded, first_read_strand, second_read_strand (STARsolo per-gene aggregate)."""
    gene_col = "geneID" if "geneID" in cols else "gene_id"

    def to_int(v):
        try:
            return int(float(v))
        except (TypeError, ValueError):
            return 0

    counts = np.array([to_int(r.get("unstranded")) for r in rows], dtype=float)
    detected = counts > 0
    log_counts = np.log10(counts[detected] + 1) if detected.any() else np.array([])

    plots: list[str] = []
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    if log_counts.size:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(log_counts, bins=60, color="#2f6f9f", edgecolor="white")
        ax.set_xlabel("log10(unstranded counts + 1) over detected genes")
        ax.set_ylabel("genes")
        ax.set_title(f"{accession}: gene-level expression (n_detected={int(detected.sum()):,}/{len(rows):,})")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_log_counts_hist.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    # Strand consistency (if present)
    if "first_read_strand" in cols and "second_read_strand" in cols:
        first = np.array([to_int(r.get("first_read_strand")) for r in rows], dtype=float)
        second = np.array([to_int(r.get("second_read_strand")) for r in rows], dtype=float)
        mask = (first + second) > 0
        if mask.any():
            ratio = first[mask] / (first[mask] + second[mask])
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(ratio, bins=40, color="#55a868", edgecolor="white")
            ax.set_xlabel("first_read_strand / (first + second)")
            ax.set_ylabel("genes")
            ax.set_title(f"{accession}: strand bias (genes with reads={int(mask.sum()):,})")
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_strand_bias.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(str(out))

    # Top expressed genes
    order = np.argsort(-counts)[:20]
    top_labels = [str(rows[i].get(gene_col, "?")) for i in order]
    top_values = counts[order]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.barh(range(len(top_labels))[::-1], top_values, color="#dd8452")
    ax.set_yticks(range(len(top_labels))[::-1])
    ax.set_yticklabels(top_labels, fontsize=7)
    ax.set_xlabel("unstranded counts")
    ax.set_title(f"{accession}: top-20 genes by counts")
    fig.tight_layout()
    out = PLOT_DIR / f"{accession}_top_genes.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    plots.append(str(out))

    # Saturation curve (Lorenz-style cumulative share)
    sorted_counts = np.sort(counts)[::-1]
    total = sorted_counts.sum()
    if total > 0:
        cumshare = np.cumsum(sorted_counts) / total
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(np.arange(1, len(cumshare) + 1), cumshare, color="#c44e52")
        ax.set_xscale("log")
        ax.set_xlabel("gene rank (log scale)")
        ax.set_ylabel("cumulative share of total counts")
        ax.set_title(f"{accession}: count concentration over genes")
        ax.axhline(0.5, ls="--", color="#888")
        ax.axhline(0.9, ls="--", color="#888")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_count_concentration.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    return {
        "accession": accession,
        "kind": "gene_quant_tsv",
        "schema": "starsolo_unstranded",
        "input": rel_path(path),
        "metadata": metadata_lines[:8],
        "n_rows": len(rows),
        "columns_detected": {"gene": gene_col, "counts": "unstranded"},
        "counts_summary": {
            "n_genes": int(len(rows)),
            "n_detected": int(detected.sum()),
            "total": float(counts.sum()),
            "max": float(counts.max()) if counts.size else 0.0,
            "median_detected": float(np.median(counts[detected])) if detected.any() else 0.0,
        },
        "top_genes_by_counts": [
            {"gene": label, "unstranded": float(value)}
            for label, value in zip(top_labels, top_values)
        ],
        "plots": [rel_path(p) for p in plots],
    }


def analyze_gene_quant_tsv(path: Path, accession: str, plt, np):
    """Analyze a scE2G/IGVF gene quantification TSV (gzipped or plain)."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    is_gzip = magic == b"\x1f\x8b"
    opener = gzip.open if is_gzip else open
    rows: list[dict[str, str]] = []
    metadata_lines: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                metadata_lines.append(line[1:].strip())
                continue
            if header is None:
                header = line.split("\t")
                continue
            rows.append(dict(zip(header, line.split("\t"))))

    def to_float(value: Any) -> float | None:
        if value is None or value == "":
            return None
        s = str(value).strip()
        if not s or s.upper() in {"NA", "NAN", "NULL", "NONE"}:
            return None
        try:
            v = float(s)
        except ValueError:
            return None
        if v != v or v == float("inf") or v == float("-inf"):
            return None
        return v

    rna_col_candidates = ("RNA_meanLogNorm", "RNA_log", "log_tpm", "logTPM")
    tpm_col_candidates = ("RNA_pseudobulkTPM", "TPM", "tpm")
    atac_col_candidates = ("normalizedATAC_prom", "promoter_atac", "atac_promoter")
    pct_col_candidates = ("RNA_percentCellsDetected", "percent_detected")
    gene_col_candidates = ("GeneSymbol", "gene", "gene_symbol", "Gene", "name")

    def pick(cols, candidates):
        for c in candidates:
            if c in cols:
                return c
        return None

    if not rows:
        return {
            "accession": accession,
            "kind": "gene_quant_tsv",
            "input": rel_path(path),
            "n_rows": 0,
            "metadata": metadata_lines[:8],
            "plots": [],
        }

    cols = list(rows[0].keys())

    # Detect STARsolo-style raw count schema: geneID + unstranded/first/second
    if "unstranded" in cols and ("geneID" in cols or "gene_id" in cols):
        return analyze_starsolo_gene_quant(path, accession, rows, cols, metadata_lines, plt, np)

    gene_col = pick(cols, gene_col_candidates)
    rna_col = pick(cols, rna_col_candidates)
    tpm_col = pick(cols, tpm_col_candidates)
    atac_col = pick(cols, atac_col_candidates)
    pct_col = pick(cols, pct_col_candidates)

    rna = np.array([to_float(r.get(rna_col)) for r in rows], dtype=object) if rna_col else None
    tpm = np.array([to_float(r.get(tpm_col)) for r in rows], dtype=object) if tpm_col else None
    atac = np.array([to_float(r.get(atac_col)) for r in rows], dtype=object) if atac_col else None
    pct = np.array([to_float(r.get(pct_col)) for r in rows], dtype=object) if pct_col else None

    def to_num(values):
        if values is None:
            return np.array([])
        out = np.array([v for v in values if v is not None], dtype=float)
        return out

    rna_v = to_num(rna)
    tpm_v = to_num(tpm)
    atac_v = to_num(atac)
    pct_v = to_num(pct)

    plots: list[str] = []
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # Plot 1: RNA distribution
    if rna_v.size > 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(rna_v, bins=50, color="#2f6f9f", edgecolor="white")
        ax.set_xlabel(rna_col)
        ax.set_ylabel("genes")
        ax.set_title(f"{accession}: gene RNA distribution\n(n={rna_v.size:,})")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_rna_hist.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    # Plot 2: TPM (log scale) distribution
    if tpm_v.size > 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        positive = tpm_v[tpm_v > 0]
        if positive.size > 0:
            ax.hist(np.log10(positive), bins=50, color="#55a868", edgecolor="white")
            ax.set_xlabel(f"log10({tpm_col} + 0)")
        else:
            ax.hist(tpm_v, bins=50, color="#55a868", edgecolor="white")
            ax.set_xlabel(tpm_col)
        ax.set_ylabel("genes")
        ax.set_title(f"{accession}: pseudobulk TPM distribution\n(n>0={positive.size:,})")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_tpm_hist.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    # Plot 3: ATAC promoter vs RNA scatter (only when both columns are present)
    if rna_col and atac_col:
        pairs = [
            (to_float(r.get(atac_col)), to_float(r.get(rna_col)), r.get(gene_col, ""))
            for r in rows
        ]
        pairs = [(a, b, g) for a, b, g in pairs if a is not None and b is not None]
        if pairs:
            xs = np.array([p[0] for p in pairs])
            ys = np.array([p[1] for p in pairs])
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(xs, ys, s=4, alpha=0.4, color="#c44e52", edgecolors="none")
            ax.set_xlabel(atac_col)
            ax.set_ylabel(rna_col)
            ax.set_title(f"{accession}: promoter ATAC vs RNA (n={len(pairs):,})")
            # Annotate top-15 by RNA
            ranked = sorted(pairs, key=lambda p: p[1], reverse=True)[:15]
            for a, b, g in ranked:
                ax.annotate(g, (a, b), fontsize=7, alpha=0.85)
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_atac_vs_rna.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(str(out))

    # Plot 4: % cells detected
    if pct_v.size > 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(pct_v, bins=50, color="#8172b2", edgecolor="white")
        ax.set_xlabel(pct_col)
        ax.set_ylabel("genes")
        ax.set_title(f"{accession}: detection fraction across cells")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_pct_detected.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    top_genes = []
    if rna_col:
        valued = [(r.get(gene_col, ""), to_float(r.get(rna_col))) for r in rows]
        valued = [(g, v) for g, v in valued if v is not None]
        valued.sort(key=lambda x: x[1], reverse=True)
        top_genes = [{"gene": g, rna_col: v} for g, v in valued[:15]]

    def stats(arr):
        if arr.size == 0:
            return None
        return {
            "n": int(arr.size),
            "min": float(arr.min()),
            "median": float(np.median(arr)),
            "mean": float(arr.mean()),
            "max": float(arr.max()),
        }

    return {
        "accession": accession,
        "kind": "gene_quant_tsv",
        "input": rel_path(path),
        "metadata": metadata_lines[:8],
        "n_rows": len(rows),
        "columns_detected": {
            "gene": gene_col,
            "rna": rna_col,
            "tpm": tpm_col,
            "promoter_atac": atac_col,
            "percent_detected": pct_col,
        },
        "rna_summary": stats(rna_v),
        "tpm_summary": stats(tpm_v),
        "promoter_atac_summary": stats(atac_v),
        "percent_detected_summary": stats(pct_v),
        "top_genes_by_rna": top_genes,
        "plots": [rel_path(p) for p in plots],
    }


# Tirosh et al. Science 2016 cell-cycle markers (human SYMBOLs).
TIROSH_S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]
TIROSH_G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


def _gene_index_by_symbol(adata, var_index, symbols_target: list[str]) -> list[int]:
    """Find variable indices that match any of the target gene symbols.

    Looks first in adata.var.gene_name (scNT-seq2 dynast files), then in
    adata.var.index (Ensembl IDs — won't match symbols, returns []).
    """
    candidate_cols = ["gene_name", "gene_symbol", "feature_name", "symbol", "Symbol"]
    series = None
    for col in candidate_cols:
        if col in adata.var.columns:
            series = adata.var[col].astype(str)
            break
    if series is None:
        return []
    upper = series.str.upper()
    targets = {s.upper() for s in symbols_target}
    return [i for i, s in enumerate(upper) if s in targets]


def publication_style_analyses(
    *, accession, adata, np, sp, plt,
    Xn, gene_means, gene_var, hvg_indices,
    cluster_labels, umap_xy, cell_idx,
    symbols, var_index, counts_per_cell,
) -> dict[str, Any]:
    """HVG plot, cluster markers (Welch t-test), dot plot, heatmap,
    cell-cycle scoring, and per-cell kinetic-fraction plots."""
    plots: list[Path] = []
    summary: dict[str, Any] = {}

    if Xn is None or cluster_labels is None or hvg_indices is None:
        return {"plots": plots, "summary": summary}

    # ----- HVG mean-variance scatter ---------------------------------
    try:
        fig, ax = plt.subplots(figsize=(6, 5))
        x_log = np.log10(np.asarray(gene_means) + 1e-6)
        y_log = np.log10(np.asarray(gene_var) + 1e-6)
        is_hvg = np.zeros_like(gene_var, dtype=bool)
        is_hvg[hvg_indices] = True
        ax.scatter(x_log[~is_hvg], y_log[~is_hvg], s=3, color="#bbbbbb", alpha=0.5, edgecolors="none", label="non-HVG")
        ax.scatter(x_log[is_hvg], y_log[is_hvg], s=4, color="#c44e52", alpha=0.7, edgecolors="none", label=f"HVG (n={int(is_hvg.sum())})")
        ax.set_xlabel("log10(mean log-norm expression + 1e-6)")
        ax.set_ylabel("log10(variance + 1e-6)")
        ax.set_title(f"{accession}: highly variable gene selection")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_hvg_mean_variance.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(out)
    except Exception as exc:
        logging.warning("HVG plot failed for %s: %s", accession, exc)

    # ----- Cluster markers via Welch t-test on log-normalized HVGs ---
    Xv = Xn[:, hvg_indices].toarray() if sp.issparse(Xn) else np.asarray(Xn)[:, hvg_indices]
    n_cells, n_hvgs = Xv.shape
    cluster_ids = np.asarray(cluster_labels)
    unique_c = np.unique(cluster_ids)

    cluster_markers: dict[int, list[dict[str, Any]]] = {}
    chosen_marker_indices: list[int] = []
    chosen_marker_clusters: list[int] = []
    chosen_marker_labels: list[str] = []
    seen_marker_idx: set[int] = set()
    n_top = 5

    for c in unique_c:
        mask = cluster_ids == c
        n1 = int(mask.sum())
        n2 = n_cells - n1
        if n1 < 3 or n2 < 3:
            continue
        m1 = Xv[mask].mean(axis=0)
        m2 = Xv[~mask].mean(axis=0)
        v1 = Xv[mask].var(axis=0, ddof=1)
        v2 = Xv[~mask].var(axis=0, ddof=1)
        denom = np.sqrt(v1 / max(n1, 2) + v2 / max(n2, 2))
        denom = np.where(denom > 0, denom, np.nan)
        t = (m1 - m2) / denom
        # Restrict to upregulated in cluster
        t = np.where(m1 > m2, t, -np.inf)
        order = np.argsort(-t)[:50]
        order = [int(i) for i in order if np.isfinite(t[i])][:n_top]
        cluster_markers[int(c)] = []
        for i in order:
            gene_idx = int(hvg_indices[i])
            label = (
                str(symbols[gene_idx]) if gene_idx < len(symbols) else str(var_index[gene_idx])
            )
            cluster_markers[int(c)].append(
                {"gene": label, "t_stat": float(t[i]), "mean_in": float(m1[i]), "mean_out": float(m2[i])}
            )
            if gene_idx not in seen_marker_idx:
                chosen_marker_indices.append(gene_idx)
                chosen_marker_clusters.append(int(c))
                chosen_marker_labels.append(label)
                seen_marker_idx.add(gene_idx)

    summary["cluster_markers"] = cluster_markers

    if chosen_marker_indices:
        marker_idx_arr = np.asarray(chosen_marker_indices)
        # ----- Marker dot plot ---------------------------------------
        try:
            Xm = Xn[:, marker_idx_arr].toarray() if sp.issparse(Xn) else np.asarray(Xn)[:, marker_idx_arr]
            n_clusters_actual = len(unique_c)
            mean_mat = np.zeros((n_clusters_actual, len(marker_idx_arr)))
            pct_mat = np.zeros((n_clusters_actual, len(marker_idx_arr)))
            for ci, c in enumerate(unique_c):
                mask = cluster_ids == c
                if mask.any():
                    block = Xm[mask]
                    mean_mat[ci] = block.mean(axis=0)
                    pct_mat[ci] = (block > 0).mean(axis=0)

            fig, ax = plt.subplots(
                figsize=(max(6, 0.35 * len(marker_idx_arr) + 2), max(3, 0.45 * n_clusters_actual + 1.5))
            )
            xs, ys, sz, cl = [], [], [], []
            for ci in range(n_clusters_actual):
                for gi in range(len(marker_idx_arr)):
                    xs.append(gi)
                    ys.append(ci)
                    sz.append(20 + 200 * pct_mat[ci, gi])
                    cl.append(mean_mat[ci, gi])
            sc = ax.scatter(xs, ys, s=sz, c=cl, cmap="Reds", edgecolors="#444", linewidths=0.4)
            ax.set_xticks(range(len(marker_idx_arr)))
            ax.set_xticklabels(chosen_marker_labels, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(n_clusters_actual))
            ax.set_yticklabels([f"c{int(c)}" for c in unique_c])
            ax.set_xlabel("marker gene")
            ax.set_ylabel("KMeans cluster")
            ax.set_title(f"{accession}: top-{n_top} markers/cluster (dot size = % expr)")
            fig.colorbar(sc, ax=ax, label="mean log-norm")
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_markers_dotplot.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(out)
        except Exception as exc:
            logging.warning("dot plot failed for %s: %s", accession, exc)

        # ----- Marker heatmap (cells × markers) ----------------------
        try:
            rng = np.random.default_rng(1)
            cells_per_cluster = 40
            sample_rows: list[int] = []
            row_cluster: list[int] = []
            for c in unique_c:
                pos = np.where(cluster_ids == c)[0]
                if pos.size == 0:
                    continue
                pick = rng.choice(pos, size=min(cells_per_cluster, pos.size), replace=False)
                sample_rows.extend(pick.tolist())
                row_cluster.extend([int(c)] * pick.size)
            sample_rows_arr = np.asarray(sample_rows)
            heat = (
                Xn[sample_rows_arr][:, marker_idx_arr].toarray()
                if sp.issparse(Xn)
                else np.asarray(Xn)[sample_rows_arr][:, marker_idx_arr]
            )
            # z-score per column for visual contrast
            mean = heat.mean(axis=0)
            std = heat.std(axis=0)
            std = np.where(std > 0, std, 1.0)
            heat_z = np.clip((heat - mean) / std, -3, 3)

            fig, ax = plt.subplots(
                figsize=(max(6, 0.35 * heat_z.shape[1] + 2), max(4, 0.04 * heat_z.shape[0] + 2))
            )
            im = ax.imshow(heat_z, aspect="auto", cmap="RdBu_r", interpolation="nearest")
            ax.set_xticks(range(heat_z.shape[1]))
            ax.set_xticklabels(chosen_marker_labels, rotation=60, ha="right", fontsize=7)
            ax.set_yticks([])
            # Side bar of cluster identity
            row_cluster_arr = np.asarray(row_cluster).reshape(-1, 1)
            cax2 = ax.twinx()
            cax2.set_yticks([])
            ax.set_xlabel("marker gene")
            ax.set_ylabel(f"cells (n={heat_z.shape[0]:,}, sorted by cluster)")
            ax.set_title(f"{accession}: marker heatmap (z-scored log-norm)")
            fig.colorbar(im, ax=ax, label="z(log-norm)")
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_markers_heatmap.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(out)
        except Exception as exc:
            logging.warning("heatmap failed for %s: %s", accession, exc)

    # ----- Cell-cycle scoring (Tirosh) — only when symbols available -
    cc_summary: dict[str, Any] = {}
    s_idx = _gene_index_by_symbol(adata, var_index, TIROSH_S_GENES)
    g2m_idx = _gene_index_by_symbol(adata, var_index, TIROSH_G2M_GENES)
    cc_summary["n_S_genes_matched"] = len(s_idx)
    cc_summary["n_G2M_genes_matched"] = len(g2m_idx)
    if (
        umap_xy is not None
        and len(s_idx) >= 5 and len(g2m_idx) >= 5
        and Xn is not None
    ):
        try:
            S_mat = Xn[:, np.asarray(s_idx)].toarray() if sp.issparse(Xn) else np.asarray(Xn)[:, np.asarray(s_idx)]
            G_mat = Xn[:, np.asarray(g2m_idx)].toarray() if sp.issparse(Xn) else np.asarray(Xn)[:, np.asarray(g2m_idx)]
            s_score = S_mat.mean(axis=1)
            g_score = G_mat.mean(axis=1)
            phase = np.where(
                (s_score < 0.05) & (g_score < 0.05),
                "G1",
                np.where(s_score > g_score, "S", "G2M"),
            )

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            axes[0].scatter(s_score, g_score, s=6, alpha=0.6, edgecolors="none", c="#2f6f9f")
            axes[0].set_xlabel("S score (mean log-norm of Tirosh S genes)")
            axes[0].set_ylabel("G2/M score")
            axes[0].set_title(f"{accession}: cell-cycle score (per cell)")
            sc2 = axes[1].scatter(
                umap_xy[:, 0], umap_xy[:, 1], s=5, alpha=0.7,
                c=g_score - s_score, cmap="coolwarm", vmin=-0.3, vmax=0.3, edgecolors="none",
            )
            axes[1].set_xlabel("UMAP-1")
            axes[1].set_ylabel("UMAP-2")
            axes[1].set_title("UMAP coloured by G2/M − S")
            fig.colorbar(sc2, ax=axes[1], label="G2/M − S")
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_cellcycle_score.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(out)

            phase_counts = {p: int((phase == p).sum()) for p in np.unique(phase)}
            cc_summary["phase_counts"] = phase_counts
            cc_summary["s_score_mean"] = float(s_score.mean())
            cc_summary["g2m_score_mean"] = float(g_score.mean())
        except Exception as exc:
            logging.warning("cell cycle plot failed for %s: %s", accession, exc)
    summary["cell_cycle"] = cc_summary

    # ----- Layer-aware kinetic fraction -----------------------------
    layers = set(adata.layers.keys())
    kinetic_summary: dict[str, Any] = {"layers_seen": sorted(layers)}
    fraction = None
    fraction_label = None

    def _layer_csr(name):
        m = adata.layers[name][cell_idx]
        return sp.csr_matrix(m) if not sp.issparse(m) else m.tocsr()

    try:
        if {"labeled_TC", "total"}.issubset(layers):
            num = np.asarray(_layer_csr("labeled_TC").sum(axis=1)).ravel()
            den = np.asarray(_layer_csr("total").sum(axis=1)).ravel()
            den = np.where(den > 0, den, 1.0)
            fraction = num / den
            fraction_label = "labeled fraction (newly synthesized RNA)"
            kinetic_summary["mode"] = "scNT-seq2 labeled/total"
        elif {"nascent", "mature"}.issubset(layers):
            n = np.asarray(_layer_csr("nascent").sum(axis=1)).ravel()
            m = np.asarray(_layer_csr("mature").sum(axis=1)).ravel()
            den = n + m
            den = np.where(den > 0, den, 1.0)
            fraction = n / den
            fraction_label = "nascent / (nascent + mature)"
            kinetic_summary["mode"] = "kb-python nascent/mature"
        elif {"unspliced", "spliced"}.issubset(layers):
            u = np.asarray(_layer_csr("unspliced").sum(axis=1)).ravel()
            s = np.asarray(_layer_csr("spliced").sum(axis=1)).ravel()
            den = u + s
            den = np.where(den > 0, den, 1.0)
            fraction = u / den
            fraction_label = "unspliced / (unspliced + spliced)"
            kinetic_summary["mode"] = "spliced/unspliced"
    except Exception as exc:
        logging.warning("kinetic layer extraction failed for %s: %s", accession, exc)

    if fraction is not None:
        try:
            fraction = np.clip(fraction, 0.0, 1.0)
            kinetic_summary["mean"] = float(fraction.mean())
            kinetic_summary["median"] = float(np.median(fraction))
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            axes[0].hist(fraction, bins=40, color="#55a868", edgecolor="white")
            axes[0].set_xlabel(fraction_label)
            axes[0].set_ylabel("cells")
            axes[0].set_title(f"{accession}: per-cell kinetic fraction\n({kinetic_summary['mode']})")
            if umap_xy is not None and umap_xy.shape[0] == fraction.shape[0]:
                sc3 = axes[1].scatter(
                    umap_xy[:, 0], umap_xy[:, 1], s=5, alpha=0.7,
                    c=fraction, cmap="viridis", vmin=0, vmax=min(1.0, float(fraction.max() if fraction.size else 1.0)),
                    edgecolors="none",
                )
                axes[1].set_xlabel("UMAP-1")
                axes[1].set_ylabel("UMAP-2")
                axes[1].set_title("UMAP coloured by kinetic fraction")
                fig.colorbar(sc3, ax=axes[1], label=fraction_label)
            else:
                axes[1].axis("off")
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_kinetic_fraction.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            plots.append(out)
        except Exception as exc:
            logging.warning("kinetic plot failed for %s: %s", accession, exc)

    summary["kinetic"] = kinetic_summary

    return {"plots": [str(p) for p in plots], "summary": summary}


def analyze_h5ad(path: Path, accession: str, np, sp, spla, plt) -> dict[str, Any]:
    import anndata  # type: ignore

    logging.info("loading h5ad: %s", path)
    adata = anndata.read_h5ad(path, backed=None)
    logging.info("loaded shape=%s", adata.shape)
    n_obs, n_vars = adata.shape
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    # Per-cell stats
    counts_per_cell = np.asarray(X.sum(axis=1)).ravel()
    binarized = X.copy()
    binarized.data = (binarized.data > 0).astype(np.float32)
    genes_per_cell = np.asarray(binarized.sum(axis=1)).ravel()

    # Per-gene stats
    counts_per_gene = np.asarray(X.sum(axis=0)).ravel()
    cells_per_gene = np.asarray(binarized.sum(axis=0)).ravel()

    # Mitochondrial fraction (best-effort)
    var_index = adata.var.index.astype(str)
    gene_symbol_col = None
    for col in ("gene_symbol", "feature_name", "symbol", "Symbol"):
        if col in adata.var.columns:
            gene_symbol_col = col
            break
    symbols = (
        adata.var[gene_symbol_col].astype(str)
        if gene_symbol_col
        else var_index
    )
    is_mito = np.array(
        [s.upper().startswith(("MT-", "MT.")) or s.upper().startswith("MT_") for s in symbols]
    )
    mito_counts = (
        np.asarray(X[:, is_mito].sum(axis=1)).ravel() if is_mito.any() else np.zeros(n_obs)
    )
    pct_mito = np.where(counts_per_cell > 0, 100 * mito_counts / counts_per_cell, 0.0)

    # Cell QC filter: keep cells with min_genes detected. For large unfiltered
    # matrices (e.g., raw barcodes), a stricter threshold drops empty droplets;
    # for already-filtered matrices, fall back to a low threshold so most cells
    # survive.
    min_genes_strict = 200
    pass_strict = int((genes_per_cell >= min_genes_strict).sum())
    if pass_strict >= 100:
        keep_cells = genes_per_cell >= min_genes_strict
        applied_min_genes = min_genes_strict
    else:
        relaxed = max(50, int(np.percentile(genes_per_cell, 10)))
        keep_cells = genes_per_cell >= relaxed
        applied_min_genes = relaxed
    n_kept = int(keep_cells.sum())
    logging.info(
        "%s: cell filter min_genes>=%d kept %d/%d cells",
        accession, applied_min_genes, n_kept, n_obs,
    )

    # Sub-sample to a manageable size for PCA/UMAP/tSNE on big inputs.
    rng = np.random.default_rng(0)
    max_cells_for_embedding = 8000
    cell_idx_pool = np.where(keep_cells)[0]
    if cell_idx_pool.size > max_cells_for_embedding:
        idx = rng.choice(cell_idx_pool, size=max_cells_for_embedding, replace=False)
        idx.sort()
    else:
        idx = cell_idx_pool

    embedding_info: dict[str, Any] = {
        "applied_min_genes": int(applied_min_genes),
        "n_cells_kept_after_filter": n_kept,
        "n_cells_used_for_embedding": int(idx.size),
    }

    Xs = X[idx] if idx.size else X[:0]
    cs = counts_per_cell[idx] if idx.size else np.array([])
    gs = genes_per_cell[idx] if idx.size else np.array([])
    ms = pct_mito[idx] if idx.size else np.array([])

    # Compute embeddings only when we have enough cells.
    pc_scores = None
    cluster_labels = None
    tsne_xy = None
    umap_xy = None
    Xn = None
    gene_means = None
    gene_var = None
    hvg_indices = None
    if Xs.shape[0] >= 50:
        target_sum = 1e4
        row_sum = np.asarray(Xs.sum(axis=1)).ravel()
        row_sum[row_sum == 0] = 1.0
        scale = target_sum / row_sum
        Xn = Xs.multiply(scale[:, None]).tocsr()
        Xn.data = np.log1p(Xn.data).astype(np.float32)

        detected = np.asarray((Xn > 0).sum(axis=0)).ravel()
        top_n_genes = min(2000, int((detected >= 10).sum()))
        if top_n_genes < 50:
            top_n_genes = min(2000, int((detected >= max(3, Xs.shape[0] // 200)).sum()))
        if top_n_genes >= 50:
            gene_means = np.asarray(Xn.mean(axis=0)).ravel()
            sq = Xn.copy()
            sq.data = sq.data ** 2
            gene_var = np.asarray(sq.mean(axis=0)).ravel() - gene_means ** 2
            gene_var = np.clip(gene_var, 0, None)
            keep = np.argsort(-gene_var)[:top_n_genes]
            hvg_indices = keep
            Xv = Xn[:, keep]
            embedding_info["n_hvgs"] = int(top_n_genes)

            n_pcs = max(2, min(50, Xv.shape[0] - 1, Xv.shape[1] - 1))
            try:
                U, S, _ = spla.svds(Xv.astype(np.float32), k=n_pcs)
                order = np.argsort(-S)
                U = U[:, order]
                S = S[order]
                pc_scores = (U * S).astype(np.float32)
                embedding_info["n_pcs"] = int(n_pcs)
            except Exception as exc:
                logging.warning("svds failed for %s: %s", accession, exc)

            if pc_scores is not None and pc_scores.shape[0] >= 50:
                from sklearn.cluster import KMeans  # type: ignore
                from sklearn.manifold import TSNE  # type: ignore
                import umap  # type: ignore

                k = max(2, min(8, pc_scores.shape[0] // 80 + 2))
                try:
                    cluster_labels = KMeans(
                        n_clusters=k, n_init=10, random_state=0
                    ).fit_predict(pc_scores)
                    embedding_info["n_clusters"] = int(k)
                except Exception as exc:
                    logging.warning("KMeans failed for %s: %s", accession, exc)

                try:
                    tsne_perplexity = max(5, min(30, (pc_scores.shape[0] - 1) // 4))
                    tsne_xy = TSNE(
                        n_components=2,
                        perplexity=tsne_perplexity,
                        init="pca",
                        learning_rate="auto",
                        random_state=0,
                    ).fit_transform(pc_scores)
                    embedding_info["tsne_perplexity"] = int(tsne_perplexity)
                except Exception as exc:
                    logging.warning("TSNE failed for %s: %s", accession, exc)

                try:
                    n_neighbors = max(5, min(15, pc_scores.shape[0] - 1))
                    umap_xy = umap.UMAP(
                        n_components=2,
                        n_neighbors=n_neighbors,
                        min_dist=0.3,
                        random_state=0,
                    ).fit_transform(pc_scores)
                    embedding_info["umap_n_neighbors"] = int(n_neighbors)
                except Exception as exc:
                    logging.warning("UMAP failed for %s: %s", accession, exc)

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plots: list[str] = []

    # Plot 1: histogram total counts per cell
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.log10(counts_per_cell + 1), bins=60, color="#2f6f9f", edgecolor="white")
    ax.set_xlabel("log10(total UMI counts + 1) per cell")
    ax.set_ylabel("cells")
    ax.set_title(f"{accession}: counts per cell  (n_cells={n_obs:,})")
    fig.tight_layout()
    out = PLOT_DIR / f"{accession}_counts_per_cell.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    plots.append(str(out))

    # Plot 2: histogram genes per cell
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(genes_per_cell, bins=60, color="#55a868", edgecolor="white")
    ax.set_xlabel("genes detected per cell")
    ax.set_ylabel("cells")
    ax.set_title(f"{accession}: genes per cell")
    fig.tight_layout()
    out = PLOT_DIR / f"{accession}_genes_per_cell.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    plots.append(str(out))

    # Plot 3: total counts vs genes, colored by % mito (filtered/sampled cells)
    if cs.size:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc = ax.scatter(
            np.log10(cs + 1),
            gs,
            c=ms,
            s=4,
            alpha=0.6,
            cmap="magma",
            vmin=0,
            vmax=min(20, float(ms.max()) if ms.size else 20),
        )
        ax.set_xlabel("log10(counts + 1)")
        ax.set_ylabel("genes detected")
        ax.set_title(f"{accession}: per-cell QC ({len(idx):,} sampled)")
        fig.colorbar(sc, ax=ax, label="% mitochondrial")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_qc_counts_genes_mito.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    # Plot 4: PCA scatter (PC1/PC2 from truncated SVD on log-normalized HVGs)
    if pc_scores is not None and pc_scores.shape[1] >= 2:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc_p = ax.scatter(
            pc_scores[:, -1], pc_scores[:, -2],
            c=np.log10(cs + 1), s=4, alpha=0.7, cmap="viridis", edgecolors="none",
        )
        ax.set_xlabel("PC1 (truncated SVD on log-normalized HVGs)")
        ax.set_ylabel("PC2")
        ax.set_title(f"{accession}: PCA ({len(idx):,} cells)")
        fig.colorbar(sc_p, ax=ax, label="log10(counts+1)")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_pca.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        plots.append(str(out))

    # Helper to render a 2D embedding twice (cluster + log counts)
    def _embedding_pair(xy, name: str):
        outs: list[Path] = []
        if xy is None:
            return outs
        # By cluster
        if cluster_labels is not None:
            fig, ax = plt.subplots(figsize=(6, 5))
            n_c = int(cluster_labels.max()) + 1
            for ci in range(n_c):
                mask = cluster_labels == ci
                ax.scatter(
                    xy[mask, 0], xy[mask, 1],
                    s=5, alpha=0.7, label=f"c{ci} (n={int(mask.sum())})",
                    color=palette(ci), edgecolors="none",
                )
            ax.set_xlabel(f"{name}-1")
            ax.set_ylabel(f"{name}-2")
            ax.set_title(f"{accession}: {name} by KMeans cluster ({len(idx):,} cells)")
            ax.legend(fontsize=7, loc="best", frameon=False)
            fig.tight_layout()
            out = PLOT_DIR / f"{accession}_{name.lower()}_clusters.png"
            fig.savefig(out, dpi=130)
            plt.close(fig)
            outs.append(out)
        # By log counts
        fig, ax = plt.subplots(figsize=(6, 5))
        sc_e = ax.scatter(
            xy[:, 0], xy[:, 1],
            c=np.log10(cs + 1), s=5, alpha=0.7,
            cmap="viridis", edgecolors="none",
        )
        ax.set_xlabel(f"{name}-1")
        ax.set_ylabel(f"{name}-2")
        ax.set_title(f"{accession}: {name} coloured by log counts")
        fig.colorbar(sc_e, ax=ax, label="log10(counts+1)")
        fig.tight_layout()
        out = PLOT_DIR / f"{accession}_{name.lower()}_counts.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        outs.append(out)
        return outs

    for path_obj in _embedding_pair(tsne_xy, "tSNE"):
        plots.append(str(path_obj))
    for path_obj in _embedding_pair(umap_xy, "UMAP"):
        plots.append(str(path_obj))

    # Plot 5: top expressed genes bar chart
    top_idx = np.argsort(-counts_per_gene)[:20]
    top_labels = [str(symbols[i]) if i < len(symbols) else str(var_index[i]) for i in top_idx]
    top_values = counts_per_gene[top_idx].astype(float)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.barh(range(len(top_labels))[::-1], top_values, color="#dd8452")
    ax.set_yticks(range(len(top_labels))[::-1])
    ax.set_yticklabels(top_labels, fontsize=8)
    ax.set_xlabel("total UMI counts")
    ax.set_title(f"{accession}: top-20 expressed genes")
    fig.tight_layout()
    out = PLOT_DIR / f"{accession}_top_genes.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    plots.append(str(out))

    # ------------------------------------------------------------------
    # Publication-style downstream analyses
    # ------------------------------------------------------------------
    pub = publication_style_analyses(
        accession=accession,
        adata=adata,
        np=np,
        sp=sp,
        plt=plt,
        Xn=Xn,
        gene_means=gene_means,
        gene_var=gene_var,
        hvg_indices=hvg_indices,
        cluster_labels=cluster_labels,
        umap_xy=umap_xy,
        cell_idx=idx,
        symbols=symbols,
        var_index=var_index,
        counts_per_cell=counts_per_cell,
    )
    plots.extend(pub["plots"])

    return {
        "accession": accession,
        "kind": "sparse_h5ad",
        "input": rel_path(path),
        "n_cells": int(n_obs),
        "n_genes": int(n_vars),
        "var_columns": list(adata.var.columns)[:12],
        "obs_columns": list(adata.obs.columns)[:12],
        "counts_per_cell": {
            "median": float(np.median(counts_per_cell)),
            "min": float(counts_per_cell.min()),
            "max": float(counts_per_cell.max()),
            "mean": float(counts_per_cell.mean()),
        },
        "genes_per_cell": {
            "median": float(np.median(genes_per_cell)),
            "min": float(genes_per_cell.min()),
            "max": float(genes_per_cell.max()),
            "mean": float(genes_per_cell.mean()),
        },
        "percent_mito": {
            "n_mito_genes": int(is_mito.sum()),
            "median": float(np.median(pct_mito)),
            "max": float(pct_mito.max()) if pct_mito.size else 0.0,
        },
        "top_expressed_genes": [
            {"gene": label, "total_counts": float(value)}
            for label, value in zip(top_labels, top_values)
        ],
        "embedding": embedding_info,
        "publication_style": pub["summary"],
        "plots": [rel_path(p) for p in plots],
    }


def write_report(records: list[dict[str, Any]], analyses: list[dict[str, Any]]) -> Path:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# IGVF Portal scRNA-seq — Ten-File Analysis")
    lines.append("")
    lines.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    lines.append("")
    lines.append("## Selection")
    lines.append("")
    lines.append(
        "| # | accession | kind | format | content_type | size (MB) | assembly | annotation |"
    )
    lines.append("|---|-----------|------|--------|--------------|-----------|----------|------------|")
    for i, rec in enumerate(records, 1):
        size_mb = float(rec.get("file_size_bytes") or 0) / 1024 / 1024
        lines.append(
            "| {i} | {acc} | {kind} | {fmt} | {ct} | {sz:.2f} | {asm} | {ann} |".format(
                i=i,
                acc=rec.get("accession", ""),
                kind=rec.get("kind", ""),
                fmt=rec.get("file_format", ""),
                ct=rec.get("content_type", ""),
                sz=size_mb,
                asm=rec.get("assembly", ""),
                ann=rec.get("transcriptome_annotation", ""),
            )
        )
    lines.append("")
    lines.append(f"Selection manifest: `{SELECTION_CSV.relative_to(ROOT)}`  ")
    lines.append(f"Download log: `{DOWNLOAD_LOG_CSV.relative_to(ROOT)}`")
    lines.append("")

    for analysis in analyses:
        acc = analysis["accession"]
        lines.append(f"## {acc} — {analysis['kind']}")
        lines.append("")
        lines.append(f"Input: `{rel_path(abs_from_rel(analysis['input']))}`")
        lines.append("")
        if analysis["kind"] == "gene_quant_tsv":
            lines.append(f"- gene rows parsed: {analysis['n_rows']:,}")
            cols = analysis.get("columns_detected", {})
            lines.append(f"- detected columns: {cols}")
            if analysis.get("schema") == "starsolo_unstranded":
                cs = analysis.get("counts_summary", {})
                lines.append(
                    f"- counts: n_genes={cs.get('n_genes',0):,} n_detected={cs.get('n_detected',0):,} "
                    f"total={cs.get('total',0):,.0f} max={cs.get('max',0):,.0f} median(detected)={cs.get('median_detected',0):.0f}"
                )
                top = analysis.get("top_genes_by_counts") or []
                if top:
                    lines.append("")
                    lines.append("Top genes (by unstranded counts):")
                    for g in top[:10]:
                        lines.append(f"  - {g['gene']}: {g['unstranded']:.0f}")
            else:
                for key in (
                    "rna_summary",
                    "tpm_summary",
                    "promoter_atac_summary",
                    "percent_detected_summary",
                ):
                    stats = analysis.get(key)
                    if stats:
                        lines.append(
                            f"- {key}: n={stats['n']:,} min={stats['min']:.3g} median={stats['median']:.3g} mean={stats['mean']:.3g} max={stats['max']:.3g}"
                        )
                top = analysis.get("top_genes_by_rna") or []
                if top:
                    rna_col = analysis["columns_detected"].get("rna") or "rna"
                    lines.append("")
                    lines.append("Top genes (by " + rna_col + "):")
                    for g in top[:10]:
                        val = g.get(rna_col)
                        lines.append(f"  - {g['gene']}: {val:.3f}")
        else:
            if "error" in analysis:
                lines.append(f"- ERROR: {analysis['error']}")
                lines.append("")
                continue
            lines.append(f"- cells: {analysis['n_cells']:,}")
            lines.append(f"- genes: {analysis['n_genes']:,}")
            cpc = analysis["counts_per_cell"]
            gpc = analysis["genes_per_cell"]
            mito = analysis["percent_mito"]
            lines.append(
                f"- counts/cell: median={cpc['median']:.0f} mean={cpc['mean']:.0f} min={cpc['min']:.0f} max={cpc['max']:.0f}"
            )
            lines.append(
                f"- genes/cell: median={gpc['median']:.0f} mean={gpc['mean']:.0f} min={gpc['min']:.0f} max={gpc['max']:.0f}"
            )
            lines.append(
                f"- mito: n_mito_genes={mito['n_mito_genes']}, median %mito={mito['median']:.2f}, max %mito={mito['max']:.2f}"
            )
            emb = analysis.get("embedding") or {}
            if emb:
                lines.append(
                    "- embedding: filtered cells (min_genes>="
                    f"{emb.get('applied_min_genes','?')}) = {emb.get('n_cells_kept_after_filter','?'):,}; "
                    f"used {emb.get('n_cells_used_for_embedding','?'):,} for tSNE/UMAP on "
                    f"{emb.get('n_pcs','?')} PCs of {emb.get('n_hvgs','?')} HVGs; "
                    f"k={emb.get('n_clusters','?')} KMeans"
                )
            pub = analysis.get("publication_style") or {}
            kin = pub.get("kinetic") or {}
            if kin.get("mode"):
                lines.append(
                    f"- kinetic fraction ({kin['mode']}): mean={kin.get('mean',0):.3f}, "
                    f"median={kin.get('median',0):.3f}; layers seen: "
                    f"{', '.join(kin.get('layers_seen', [])[:8])}"
                )
            cc = pub.get("cell_cycle") or {}
            if cc.get("phase_counts"):
                pc_str = ", ".join(f"{p}={n}" for p, n in cc["phase_counts"].items())
                lines.append(
                    f"- cell-cycle (Tirosh): S genes matched={cc.get('n_S_genes_matched',0)}, "
                    f"G2M matched={cc.get('n_G2M_genes_matched',0)}, "
                    f"phases: {pc_str}"
                )
            elif cc.get("n_S_genes_matched", 0) == 0 and cc.get("n_G2M_genes_matched", 0) == 0:
                lines.append("- cell-cycle: skipped (no gene_name in var; only Ensembl IDs)")
            cm = pub.get("cluster_markers") or {}
            if cm:
                lines.append("")
                lines.append("Top cluster markers (Welch t-test, top 5 per cluster, log-norm HVGs):")
                for cid in sorted(cm.keys()):
                    names = ", ".join(g["gene"] for g in cm[cid][:5])
                    lines.append(f"  - c{cid}: {names}")
            lines.append("")
            lines.append("Top expressed genes:")
            for g in analysis["top_expressed_genes"][:10]:
                lines.append(f"  - {g['gene']}: {g['total_counts']:.0f}")
        lines.append("")
        lines.append("Plots:")
        for plot in analysis.get("plots", []):
            plot_abs = abs_from_rel(plot)
            try:
                rel = plot_abs.resolve().relative_to(REPORT_PATH.parent.resolve()).as_posix()
            except ValueError:
                rel = rel_path(plot_abs)
            lines.append(f"- ![{Path(plot).stem}]({rel})")
        lines.append("")

    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote report: %s", REPORT_PATH)
    return REPORT_PATH


def analyze_selection() -> Path:
    np, pd, sp, spla, plt = _import_analysis_libs()
    records = read_selection()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    analyses: list[dict[str, Any]] = []
    for rec in records:
        local = abs_from_rel(rec["local_path"])
        if not local.exists() or local.stat().st_size == 0:
            logging.warning("skipping (missing): %s", local)
            continue
        if rec["kind"] == "gene_quant_tsv":
            try:
                summary = analyze_gene_quant_tsv(local, rec["accession"], plt, np)
            except Exception as exc:
                logging.exception("gene_quant_tsv analysis failed for %s", rec["accession"])
                summary = {
                    "accession": rec["accession"],
                    "kind": rec["kind"],
                    "input": rel_path(local),
                    "error": str(exc),
                    "plots": [],
                }
        else:
            try:
                summary = analyze_h5ad(local, rec["accession"], np, sp, spla, plt)
            except Exception as exc:
                logging.exception("h5ad analysis failed for %s", rec["accession"])
                summary = {
                    "accession": rec["accession"],
                    "kind": rec["kind"],
                    "input": rel_path(local),
                    "error": str(exc),
                    "plots": [],
                }
        analyses.append(summary)
        plots_str = ", ".join(Path(p).name for p in summary.get("plots", []))
        print(f"analyzed {rec['accession']:18s} {rec['kind']:18s} plots: {plots_str}")

    ANALYSIS_JSON.parent.mkdir(parents=True, exist_ok=True)
    ANALYSIS_JSON.write_text(json.dumps(analyses, indent=2, default=str), encoding="utf-8")
    logging.info("Wrote analysis JSON: %s", ANALYSIS_JSON)
    write_report(records, analyses)
    return REPORT_PATH


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sel = sub.add_parser("select", help="Query portal and write the 10-file selection manifest.")
    sel.add_argument("--h5ad-count", type=int, default=4, help="Number of small h5ad files to include.")

    dl = sub.add_parser("download", help="Download the selected files into Data/SingleCell/Portal10/")
    dl.add_argument("--force", action="store_true", help="Re-download even if files already exist.")

    sub.add_parser("analyze", help="Analyze downloaded files, write plots and report.")
    sub.add_parser("run", help="Run select -> download -> analyze in one shot.")

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "select":
        records = select_files(h5ad_count=args.h5ad_count)
        write_selection(records)
        for rec in records:
            print(
                f"{rec['accession']:18s} {rec['kind']:18s} "
                f"{(rec.get('file_size_bytes') or 0) / 1024 / 1024:8.2f} MB "
                f"{rec['download_url']}"
            )
        return 0
    if args.command == "download":
        download_selection(force=args.force)
        return 0
    if args.command == "analyze":
        analyze_selection()
        return 0
    if args.command == "run":
        records = select_files()
        write_selection(records)
        download_selection(force=False)
        analyze_selection()
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
