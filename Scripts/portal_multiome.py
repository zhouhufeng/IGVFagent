#!/usr/bin/env python3
"""Download and analyze ten 10x Multiome AnalysisSets from the IGVF Portal.

Each AnalysisSet ships a quartet of files:
  - sparse gene count matrix       (RNA, Cell Ranger ARC tar.gz triplet)
  - annotated sparse peak count matrix  (ATAC, .rds — R only, skipped here)
  - fragments                       (BED.gz — used for ATAC QC on a subset)
  - cell annotations                (TSV.gz — IGVF cell-type labels per barcode)

We download:
  - all 10 sets: RNA tar + cell annotations
  - 3 smallest sets:  + fragments BED for ATAC QC and nucleosomal-banding plot

Outputs (paths repo-relative; no absolute paths leaked):
  Data/SingleCell/Multiome10/<acc>/...                raw downloads
  Data/Manifests/Multiome10/portal_multiome_selection.csv
  Data/Manifests/Multiome10/portal_multiome_download_log.csv
  Data/SingleCell/Multiome10/analysis_summary.json    machine-readable summary
  Docs/SingleCell/Plots/Multiome10/<acc>/*.png        per-set plots
  Docs/SingleCell/Multiome10_report.md                final report
  Docs/Logs/portal_multiome_*.log                     runtime log
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import os
import sys
import tarfile
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
DOWNLOAD_DIR = DATA_DIR / "SingleCell" / "Multiome10"
MANIFEST_DIR = DATA_DIR / "Manifests" / "Multiome10"
PLOT_DIR = DOCS_DIR / "SingleCell" / "Plots" / "Multiome10"
REPORT_PATH = DOCS_DIR / "SingleCell" / "Multiome10_report.md"
SURVEY_PATH = DOCS_DIR / "SingleCell" / "Multiome10_survey.md"
SELECTION_CSV = MANIFEST_DIR / "portal_multiome_selection.csv"
DOWNLOAD_LOG = MANIFEST_DIR / "portal_multiome_download_log.csv"
ANALYSIS_JSON = DOWNLOAD_DIR / "analysis_summary.json"

PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://api.data.igvf.org").rstrip("/")
USER_AGENT = "IGVFagent-portal-multiome/0.1"

# Standard cell-cycle marker symbols (Tirosh et al. Science 2016).
TIROSH_S = [
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
    "CDCA7","DTL","PRIM1","UHRF1","CENPU","HELLS","RFC2","RPA2","NASP",
    "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
    "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
    "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8",
]
TIROSH_G2M = [
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","PIMREG","SMC4","CCNB2",
    "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
    "KIF20B","HJURP","CDCA3","JPT1","CDC20","TTK","CDC25C","KIF2C","RANGAP1",
    "NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1",
    "ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA",
]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"portal_multiome_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def rel_path(path: Path | str) -> str:
    p = Path(path).resolve()
    try:
        return p.relative_to(ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def abs_from_rel(rel: str) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else (ROOT / p)


def palette(i: int) -> str:
    cols = [
        "#2f6f9f","#c44e52","#55a868","#8172b2","#ccb974","#64b5cd",
        "#8c613c","#dc7ec0","#4c72b0","#dd8452","#937860","#6acc64",
    ]
    return cols[i % len(cols)]


def http_get_json(url: str, timeout: int = 60) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    logging.info("GET %s", url)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def latest_files_manifest() -> Path | None:
    """Return the most recent files manifest produced by multiome_10x_pipeline."""
    cands = sorted(MANIFEST_DIR.glob("*_igvf_10x_multiome_files.csv"))
    if not cands:
        # Fall back to any pre-existing manifest dir from the original skill
        cands = sorted((DATA_DIR / "Manifests" / "Multiome10x").glob("*_igvf_10x_multiome_files.csv"))
    return cands[-1] if cands else None


def select_sets(n_sets: int = 10, n_with_fragments: int = 3) -> list[dict[str, Any]]:
    """Pick the n_sets smallest AnalysisSets (by tar RNA + cell annotations).

    The first n_with_fragments of those also get their fragments BED
    downloaded for ATAC QC.
    """
    manifest = latest_files_manifest()
    if manifest is None:
        raise SystemExit(
            "No 10x multiome files manifest found. Run "
            "`Scripts/multiome_10x_pipeline.py retrieve --count 30 --fetch-file-details` first."
        )
    logging.info("Reading manifest: %s", manifest)
    rows = list(csv.DictReader(manifest.open()))
    by_set: dict[str, dict[str, dict[str, str]]] = {}
    for r in rows:
        ct = r.get("content_type", "")
        by_set.setdefault(r["analysis_set_accession"], {})[ct] = r

    candidates = []
    for acc, fmap in by_set.items():
        rna = fmap.get("sparse gene count matrix")
        ann = fmap.get("cell annotations")
        frag = fmap.get("fragments")
        if not (rna and ann):
            continue
        if rna.get("file_format") != "tar":  # restrict to readable tar triplet
            continue
        rna_sz = int(float(rna.get("file_size_bytes") or 0))
        ann_sz = int(float(ann.get("file_size_bytes") or 0))
        frag_sz = int(float((frag or {}).get("file_size_bytes") or 0))
        candidates.append({
            "accession": acc,
            "rna": rna,
            "ann": ann,
            "frag": frag,
            "small_size": rna_sz + ann_sz,
            "frag_size": frag_sz,
        })
    candidates.sort(key=lambda x: x["small_size"])
    chosen = candidates[:n_sets]
    chosen.sort(key=lambda x: x["small_size"])
    for i, item in enumerate(chosen):
        item["include_fragments"] = i < n_with_fragments
    return chosen


def write_selection(items: list[dict[str, Any]]) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "accession", "kind", "file_accession", "file_format", "content_type",
        "file_size_bytes", "download_url", "local_path",
    ]

    def kind_path(set_acc: str, kind: str, file_acc: str, fmt: str) -> Path:
        suffix_map = {"rna_tar": "tar.gz", "annotations": "tsv.gz", "fragments": "tsv.gz"}
        suffix = suffix_map.get(kind, fmt)
        return DOWNLOAD_DIR / set_acc / f"{file_acc}.{suffix}"

    out_rows = []
    for item in items:
        for kind, file_dict in (("rna_tar", item["rna"]), ("annotations", item["ann"]), ("fragments", item["frag"] if item["include_fragments"] else None)):
            if file_dict is None:
                continue
            local = kind_path(item["accession"], kind, file_dict["file_accession"], file_dict["file_format"])
            out_rows.append({
                "accession": item["accession"],
                "kind": kind,
                "file_accession": file_dict["file_accession"],
                "file_format": file_dict["file_format"],
                "content_type": file_dict["content_type"],
                "file_size_bytes": int(float(file_dict.get("file_size_bytes") or 0)),
                "download_url": file_dict.get("download_url", ""),
                "local_path": rel_path(local),
            })

    with SELECTION_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    logging.info("Wrote selection: %s (%d rows)", SELECTION_CSV, len(out_rows))
    return SELECTION_CSV


def read_selection() -> list[dict[str, Any]]:
    if not SELECTION_CSV.exists():
        raise SystemExit(f"selection missing: {SELECTION_CSV}. Run 'select' first.")
    with SELECTION_CSV.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download_one(url: str, target: Path, force: bool = False) -> dict[str, Any]:
    target.parent.mkdir(parents=True, exist_ok=True)
    rel = rel_path(target)
    if target.exists() and not force and target.stat().st_size > 0:
        return {"status": "exists", "bytes": target.stat().st_size, "url": url, "path": rel}
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    started = time.time()
    try:
        with urllib.request.urlopen(req, timeout=600) as r, target.open("wb") as fh:
            total = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                total += len(chunk)
        elapsed = time.time() - started
        logging.info("downloaded %s (%.2f MB in %.1fs)", target.name, total/1024/1024, elapsed)
        return {"status": "downloaded", "bytes": total, "url": url, "path": rel}
    except urllib.error.HTTPError as exc:
        return {"status": f"http_error_{exc.code}", "bytes": 0, "url": url, "path": rel}
    except urllib.error.URLError as exc:
        return {"status": f"network_error_{exc.reason}", "bytes": 0, "url": url, "path": rel}


def download_selection(force: bool = False) -> Path:
    rows = read_selection()
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    log_rows = []
    for r in rows:
        result = download_one(r["download_url"], abs_from_rel(r["local_path"]), force=force)
        log_rows.append({
            "accession": r["accession"], "kind": r["kind"],
            "file_accession": r["file_accession"], "status": result["status"],
            "bytes": result["bytes"], "url": result["url"], "path": result["path"],
        })
        print(
            f"{r['accession']:18s} {r['kind']:14s} {r['file_accession']:14s} "
            f"{result['status']:>14s} {result['bytes']:>13,} bytes"
        )
    with DOWNLOAD_LOG.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log_rows[0].keys()))
        w.writeheader()
        for row in log_rows:
            w.writerow(row)
    logging.info("Wrote download log: %s", DOWNLOAD_LOG)
    return DOWNLOAD_LOG


# ---------------------------------------------------------------------------
# Cell Ranger ARC tar.gz reader
# ---------------------------------------------------------------------------


def read_tenx_tar(tar_path: Path, np, sp):
    """Parse a Cell Ranger ARC tar.gz containing matrix.mtx[.gz], features.tsv[.gz],
    barcodes.tsv[.gz]. Returns (X csr matrix, var DataFrame, obs DataFrame)."""
    import scipy.io  # type: ignore
    import pandas as pd  # type: ignore

    barcodes = features = matrix = None
    with tarfile.open(tar_path, "r:*") as tar:
        members = tar.getmembers()

        def _extract(name_keys):
            for m in members:
                base = Path(m.name).name.lower()
                if any(k in base for k in name_keys):
                    f = tar.extractfile(m)
                    if f is None:
                        continue
                    raw = f.read()
                    if base.endswith(".gz"):
                        raw = gzip.decompress(raw)
                    return raw, m.name
            return None, None

        # 10x Cell Ranger ARC default names + IGVF DACC '<sample>_<tissue>_*' style
        b_raw, b_name = _extract(("barcodes.tsv",))
        f_raw, f_name = _extract(("features.tsv", "genes.tsv"))
        m_raw, m_name = _extract(("matrix.mtx", "counts.mtx"))

    if b_raw is None or f_raw is None or m_raw is None:
        raise RuntimeError(f"Missing 10x triplet inside {tar_path.name}")

    barcodes = [line.strip() for line in b_raw.decode().splitlines() if line.strip()]
    feature_rows = [line.split("\t") for line in f_raw.decode().splitlines() if line.strip()]
    n_feat_cols = max((len(r) for r in feature_rows), default=1)
    if n_feat_cols == 1:
        # IGVF DACC pipelines often write features.tsv with a single column.
        # Detect Ensembl IDs vs gene symbols heuristically.
        sample = [r[0] for r in feature_rows[:200] if r]
        looks_ensembl = sum(1 for s in sample if s.startswith(("ENSG", "ENSMUSG"))) > len(sample) * 0.5
        features_cols = ["gene_id"] if looks_ensembl else ["gene_name"]
    else:
        features_cols = ["gene_id", "gene_name", "feature_type"][:n_feat_cols]
    var = pd.DataFrame(feature_rows, columns=features_cols)

    X = scipy.io.mmread(io.BytesIO(m_raw)).tocsr()
    if X.shape[0] == len(var) and X.shape[1] == len(barcodes):
        # 10x stores genes × cells; transpose to cells × genes
        X = X.T.tocsr()

    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    return X, var, obs, {"barcodes_file": b_name, "features_file": f_name, "matrix_file": m_name}


def read_annotation_tsv(path: Path):
    """Read IGVF cell annotations TSV (gzipped) → pandas DataFrame indexed by barcode."""
    import pandas as pd  # type: ignore
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    # IGVF cell-annotations TSVs commonly have 'barcode' or 'cell_barcode' column
    for col in ("barcode", "cell_barcode", "cell_id", "cellID", "Cell"):
        if col in df.columns:
            df = df.set_index(col)
            break
    return df


# ---------------------------------------------------------------------------
# Per-set analyses
# ---------------------------------------------------------------------------


def _import_libs():
    import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla  # type: ignore
    import matplotlib  # type: ignore
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    import pandas as pd  # type: ignore
    return np, sp, spla, plt, pd


def analyze_one_set(item: dict[str, Any], np, sp, spla, plt, pd) -> dict[str, Any]:
    set_acc = item["accession"]
    set_dir = DOWNLOAD_DIR / set_acc
    plot_dir = PLOT_DIR / set_acc
    plot_dir.mkdir(parents=True, exist_ok=True)
    plots: list[Path] = []

    summary: dict[str, Any] = {"accession": set_acc}

    # Locate downloaded files for this set from the selection CSV
    rows = read_selection()
    files = {r["kind"]: abs_from_rel(r["local_path"]) for r in rows if r["accession"] == set_acc}

    # ---- RNA: load tar triplet ------------------------------------------
    rna_path = files.get("rna_tar")
    if rna_path is None or not rna_path.exists():
        return {"accession": set_acc, "error": "no rna tar", "plots": []}
    logging.info("[%s] loading RNA from %s", set_acc, rna_path)
    X, var, obs, file_meta = read_tenx_tar(rna_path, np, sp)
    n_cells, n_genes = X.shape
    summary["rna"] = {"n_cells": int(n_cells), "n_genes": int(n_genes)}

    # ---- Per-cell QC ----------------------------------------------------
    counts_per_cell = np.asarray(X.sum(axis=1)).ravel()
    binarized = X.copy()
    binarized.data = (binarized.data > 0).astype(np.float32)
    genes_per_cell = np.asarray(binarized.sum(axis=1)).ravel()
    counts_per_gene = np.asarray(X.sum(axis=0)).ravel()

    gene_names = var["gene_name"].astype(str) if "gene_name" in var.columns else var.index.astype(str)
    is_mito = np.array([s.upper().startswith(("MT-", "MT.", "MT_")) for s in gene_names])
    mito_counts = np.asarray(X[:, is_mito].sum(axis=1)).ravel() if is_mito.any() else np.zeros(n_cells)
    pct_mito = np.where(counts_per_cell > 0, 100 * mito_counts / counts_per_cell, 0.0)

    summary["rna"]["counts_per_cell_median"] = float(np.median(counts_per_cell))
    summary["rna"]["genes_per_cell_median"] = float(np.median(genes_per_cell))
    summary["rna"]["pct_mito_median"] = float(np.median(pct_mito))

    # QC plot 1: counts vs genes coloured by % mito (sub-sampled for speed)
    rng = np.random.default_rng(0)
    sub_idx = (
        rng.choice(n_cells, size=min(8000, n_cells), replace=False) if n_cells > 8000 else np.arange(n_cells)
    )
    sub_idx.sort()
    fig, ax = plt.subplots(figsize=(6, 5))
    sc_ = ax.scatter(
        np.log10(counts_per_cell[sub_idx] + 1), genes_per_cell[sub_idx],
        c=pct_mito[sub_idx], cmap="magma", s=4, alpha=0.6, vmin=0, vmax=20,
    )
    ax.set_xlabel("log10(counts + 1)"); ax.set_ylabel("genes detected")
    ax.set_title(f"{set_acc}: per-cell QC ({len(sub_idx):,} sampled)")
    fig.colorbar(sc_, ax=ax, label="% mitochondrial")
    fig.tight_layout()
    out = plot_dir / "rna_qc_counts_genes_mito.png"
    fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

    # ---- Cell filter + embedding (RNA) ---------------------------------
    keep = genes_per_cell >= 200
    if keep.sum() < 100:
        keep = genes_per_cell >= max(50, int(np.percentile(genes_per_cell, 10)))
    pool = np.where(keep)[0]
    rng = np.random.default_rng(0)
    if pool.size > 5000:
        idx = rng.choice(pool, size=5000, replace=False); idx.sort()
    else:
        idx = pool
    summary["rna"]["n_cells_kept"] = int(keep.sum())
    summary["rna"]["n_cells_embedded"] = int(idx.size)

    pc_scores = cluster_labels = umap_xy = tsne_xy = None
    Xn = None; hvg_idx = None
    if idx.size >= 50:
        Xs = X[idx]
        row_sum = np.asarray(Xs.sum(axis=1)).ravel()
        row_sum[row_sum == 0] = 1.0
        Xn = Xs.multiply((1e4 / row_sum)[:, None]).tocsr()
        Xn.data = np.log1p(Xn.data).astype(np.float32)

        detected = np.asarray((Xn > 0).sum(axis=0)).ravel()
        n_hvg = min(2000, int((detected >= 10).sum()))
        if n_hvg < 50:
            n_hvg = min(2000, int((detected >= 3).sum()))
        if n_hvg >= 50:
            gm = np.asarray(Xn.mean(axis=0)).ravel()
            sq = Xn.copy(); sq.data = sq.data ** 2
            gv = np.asarray(sq.mean(axis=0)).ravel() - gm ** 2
            gv = np.clip(gv, 0, None)
            hvg_idx = np.argsort(-gv)[:n_hvg]
            Xv = Xn[:, hvg_idx]
            n_pc = max(2, min(50, Xv.shape[0] - 1, Xv.shape[1] - 1))
            try:
                U, S, _ = spla.svds(Xv.astype(np.float32), k=n_pc)
                order = np.argsort(-S)
                pc_scores = (U[:, order] * S[order]).astype(np.float32)
            except Exception as exc:
                logging.warning("[%s] svds failed: %s", set_acc, exc)

        if pc_scores is not None and pc_scores.shape[0] >= 50:
            from sklearn.cluster import KMeans  # type: ignore
            from sklearn.manifold import TSNE  # type: ignore
            import umap  # type: ignore
            try:
                k = max(2, min(8, pc_scores.shape[0] // 80 + 2))
                cluster_labels = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(pc_scores)
                tsne_xy = TSNE(
                    n_components=2,
                    perplexity=max(5, min(30, (pc_scores.shape[0] - 1) // 4)),
                    init="pca", learning_rate="auto", random_state=0,
                ).fit_transform(pc_scores)
                umap_xy = umap.UMAP(
                    n_components=2,
                    n_neighbors=max(5, min(15, pc_scores.shape[0] - 1)),
                    min_dist=0.3, random_state=0,
                ).fit_transform(pc_scores)
                summary["rna"]["n_clusters"] = int(k)
            except Exception as exc:
                logging.warning("[%s] embedding failed: %s", set_acc, exc)

    # ---- UMAP plots ------------------------------------------------------
    cs = counts_per_cell[idx] if idx.size else np.array([])
    if umap_xy is not None and cluster_labels is not None:
        fig, ax = plt.subplots(figsize=(6, 5))
        for ci in range(int(cluster_labels.max()) + 1):
            m = cluster_labels == ci
            ax.scatter(umap_xy[m, 0], umap_xy[m, 1], s=6, alpha=0.7,
                       label=f"c{ci} (n={int(m.sum())})", color=palette(ci), edgecolors="none")
        ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
        ax.set_title(f"{set_acc}: RNA UMAP by KMeans cluster ({len(idx):,} cells)")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        out = plot_dir / "rna_umap_clusters.png"
        fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

        fig, ax = plt.subplots(figsize=(6, 5))
        scc = ax.scatter(umap_xy[:, 0], umap_xy[:, 1], s=5, alpha=0.7,
                          c=np.log10(cs + 1), cmap="viridis", edgecolors="none")
        ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
        ax.set_title(f"{set_acc}: RNA UMAP coloured by log counts")
        fig.colorbar(scc, ax=ax, label="log10(counts+1)")
        fig.tight_layout()
        out = plot_dir / "rna_umap_counts.png"
        fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

    # ---- Cell annotations + UMAP coloured by IGVF cell-type ------------
    ann_path = files.get("annotations")
    cell_type_summary: dict[str, Any] = {}
    if ann_path is not None and ann_path.exists():
        try:
            ann = read_annotation_tsv(ann_path)
            ann_cols = [c for c in ann.columns]
            cell_type_summary["columns"] = ann_cols[:20]
            label_col = None
            for cand in ("cell_type", "celltype", "Cell_type", "cluster", "subclass", "subclass_label", "leiden", "annotation"):
                if cand in ann.columns:
                    label_col = cand
                    break
            if label_col is None and ann_cols:
                # heuristic: first column with low cardinality strings
                for c in ann_cols:
                    series = ann[c].astype(str)
                    nunique = series.nunique()
                    if 2 <= nunique <= 80:
                        label_col = c
                        break
            cell_type_summary["label_column"] = label_col

            if label_col is not None and umap_xy is not None:
                # Match annotation barcodes to RNA matrix barcodes
                rna_barcodes = obs.index.astype(str)
                ann_index = ann.index.astype(str)
                # Some IGVF annotation barcodes carry a sample suffix; try a quick alignment
                shared = ann_index.intersection(rna_barcodes)
                if shared.size < 0.1 * rna_barcodes.size:
                    # try stripping common suffixes
                    rna_b_set = set(rna_barcodes)
                    def _strip(s):
                        return s.split("-")[0]
                    candidate = ann_index.map(_strip)
                    if any(c in rna_b_set for c in candidate[:200]):
                        ann.index = candidate
                        ann_index = ann.index
                        shared = ann_index.intersection(rna_barcodes)
                cell_type_summary["barcodes_overlap"] = int(shared.size)

                if shared.size >= 50:
                    # Build a label vector aligned to embedded cells (idx)
                    rna_b_to_pos = {b: i for i, b in enumerate(rna_barcodes)}
                    label_full = pd.Series([None] * n_cells, index=rna_barcodes)
                    label_full.loc[shared] = ann.loc[shared, label_col].astype(str).values
                    labels_emb = label_full.iloc[idx].values
                    unique_labels, counts = np.unique([l for l in labels_emb if l is not None], return_counts=True)
                    cell_type_summary["unique_labels"] = [
                        {"label": str(l), "count": int(c)} for l, c in zip(unique_labels, counts)
                    ][:25]

                    # plot UMAP coloured by cell type label (top 12 by count)
                    top_labels = list(unique_labels[np.argsort(-counts)][:12])
                    fig, ax = plt.subplots(figsize=(7, 5))
                    other = np.array([(l not in top_labels) for l in labels_emb])
                    if other.any():
                        ax.scatter(umap_xy[other, 0], umap_xy[other, 1], s=5, alpha=0.4,
                                   color="#dddddd", edgecolors="none", label="other/none")
                    for i, lab in enumerate(top_labels):
                        m = np.array([str(l) == str(lab) for l in labels_emb])
                        if m.any():
                            ax.scatter(umap_xy[m, 0], umap_xy[m, 1], s=6, alpha=0.85,
                                       color=palette(i), edgecolors="none",
                                       label=f"{lab} (n={int(m.sum())})")
                    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
                    ax.set_title(f"{set_acc}: UMAP by IGVF cell-type label ({label_col})")
                    ax.legend(fontsize=6, loc="best", frameon=False)
                    fig.tight_layout()
                    out = plot_dir / "rna_umap_celltype.png"
                    fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)
        except Exception as exc:
            logging.warning("[%s] annotation processing failed: %s", set_acc, exc)
    summary["cell_annotations"] = cell_type_summary

    # ---- Cluster markers via Welch t-test ------------------------------
    cluster_markers_summary: list[dict[str, Any]] = []
    chosen_marker_idx: list[int] = []
    chosen_marker_labels: list[str] = []
    if Xn is not None and cluster_labels is not None and hvg_idx is not None:
        try:
            Xv = Xn[:, hvg_idx].toarray()
            for c in np.unique(cluster_labels):
                mask = cluster_labels == c
                n1, n2 = int(mask.sum()), int((~mask).sum())
                if n1 < 3 or n2 < 3:
                    continue
                m1 = Xv[mask].mean(axis=0); m2 = Xv[~mask].mean(axis=0)
                v1 = Xv[mask].var(axis=0, ddof=1); v2 = Xv[~mask].var(axis=0, ddof=1)
                denom = np.sqrt(v1 / max(n1, 2) + v2 / max(n2, 2))
                denom = np.where(denom > 0, denom, np.nan)
                t = np.where(m1 > m2, (m1 - m2) / denom, -np.inf)
                top = np.argsort(-t)[:5]
                top = [int(i) for i in top if np.isfinite(t[i])]
                cluster_markers = []
                for i in top:
                    gid = int(hvg_idx[i])
                    label = str(gene_names[gid])
                    cluster_markers.append({
                        "gene": label, "t": float(t[i]),
                        "mean_in": float(m1[i]), "mean_out": float(m2[i]),
                    })
                    if gid not in chosen_marker_idx:
                        chosen_marker_idx.append(gid)
                        chosen_marker_labels.append(label)
                cluster_markers_summary.append({"cluster": int(c), "markers": cluster_markers})
        except Exception as exc:
            logging.warning("[%s] markers failed: %s", set_acc, exc)
    summary["cluster_markers"] = cluster_markers_summary

    # ---- Marker dot plot -----------------------------------------------
    if Xn is not None and cluster_labels is not None and chosen_marker_idx:
        try:
            Xm = Xn[:, np.asarray(chosen_marker_idx)].toarray()
            unique_c = np.unique(cluster_labels)
            mean_mat = np.zeros((len(unique_c), len(chosen_marker_idx)))
            pct_mat = np.zeros_like(mean_mat)
            for ci, c in enumerate(unique_c):
                m = cluster_labels == c
                if m.any():
                    mean_mat[ci] = Xm[m].mean(axis=0)
                    pct_mat[ci] = (Xm[m] > 0).mean(axis=0)
            fig, ax = plt.subplots(figsize=(max(7, 0.32 * len(chosen_marker_idx) + 2),
                                            max(3.5, 0.45 * len(unique_c) + 1.5)))
            xs, ys, sz, cl = [], [], [], []
            for ci in range(len(unique_c)):
                for gi in range(len(chosen_marker_idx)):
                    xs.append(gi); ys.append(ci)
                    sz.append(20 + 200 * pct_mat[ci, gi])
                    cl.append(mean_mat[ci, gi])
            sc2 = ax.scatter(xs, ys, s=sz, c=cl, cmap="Reds", edgecolors="#333", linewidths=0.4)
            ax.set_xticks(range(len(chosen_marker_idx)))
            ax.set_xticklabels(chosen_marker_labels, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(len(unique_c)))
            ax.set_yticklabels([f"c{int(c)}" for c in unique_c])
            ax.set_xlabel("marker"); ax.set_ylabel("cluster")
            ax.set_title(f"{set_acc}: cluster marker dot plot")
            fig.colorbar(sc2, ax=ax, label="mean log-norm")
            fig.tight_layout()
            out = plot_dir / "rna_markers_dotplot.png"
            fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)
        except Exception as exc:
            logging.warning("[%s] dot plot failed: %s", set_acc, exc)

    # ---- Cell-cycle score ---------------------------------------------
    cc: dict[str, Any] = {}
    if Xn is not None and umap_xy is not None and "gene_name" in var.columns:
        try:
            sym = var["gene_name"].astype(str).str.upper()
            s_idx = np.where(sym.isin([s.upper() for s in TIROSH_S]))[0]
            g_idx = np.where(sym.isin([s.upper() for s in TIROSH_G2M]))[0]
            cc["s_genes_matched"] = int(s_idx.size)
            cc["g2m_genes_matched"] = int(g_idx.size)
            if s_idx.size >= 5 and g_idx.size >= 5:
                S_mat = Xn[:, s_idx].toarray()
                G_mat = Xn[:, g_idx].toarray()
                s_score = S_mat.mean(axis=1); g_score = G_mat.mean(axis=1)
                phase = np.where((s_score < 0.05) & (g_score < 0.05), "G1",
                                 np.where(s_score > g_score, "S", "G2M"))
                fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
                axes[0].scatter(s_score, g_score, s=5, alpha=0.6, edgecolors="none", c="#2f6f9f")
                axes[0].set_xlabel("S score"); axes[0].set_ylabel("G2M score")
                axes[0].set_title(f"{set_acc}: cell-cycle score")
                sc3 = axes[1].scatter(umap_xy[:, 0], umap_xy[:, 1], s=5, alpha=0.7,
                                      c=g_score - s_score, cmap="coolwarm", vmin=-0.3, vmax=0.3,
                                      edgecolors="none")
                axes[1].set_xlabel("UMAP-1"); axes[1].set_ylabel("UMAP-2")
                axes[1].set_title("UMAP: G2M − S")
                fig.colorbar(sc3, ax=axes[1], label="G2M − S")
                fig.tight_layout()
                out = plot_dir / "rna_cellcycle.png"
                fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)
                vals, counts = np.unique(phase, return_counts=True)
                cc["phase_counts"] = {str(v): int(c) for v, c in zip(vals, counts)}
        except Exception as exc:
            logging.warning("[%s] cell-cycle failed: %s", set_acc, exc)
    summary["cell_cycle"] = cc

    # ---- ATAC fragments QC (only if downloaded) ------------------------
    frag_path = files.get("fragments")
    atac: dict[str, Any] = {"fragments_present": bool(frag_path and frag_path.exists())}
    if frag_path and frag_path.exists():
        try:
            sizes, frag_per_barcode = sample_fragments(frag_path, max_lines=2_000_000)
            atac["fragments_sampled"] = int(sum(frag_per_barcode.values()))
            atac["unique_barcodes_sampled"] = int(len(frag_per_barcode))
            # Plot 1: fragment size distribution (nucleosomal banding)
            fig, ax = plt.subplots(figsize=(7, 4.5))
            sizes_arr = np.asarray(sizes)
            sizes_arr = sizes_arr[(sizes_arr > 0) & (sizes_arr < 1000)]
            ax.hist(sizes_arr, bins=200, color="#2f6f9f", edgecolor="none")
            ax.set_xlabel("fragment size (bp)")
            ax.set_ylabel("fragments (sampled)")
            ax.set_title(f"{set_acc}: fragment size distribution (n={sizes_arr.size:,} sampled)")
            for x in (147, 294, 441):
                ax.axvline(x, color="#c44e52", linestyle="--", linewidth=0.7, alpha=0.6)
            fig.tight_layout()
            out = plot_dir / "atac_fragment_size.png"
            fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

            # Plot 2: fragments per barcode rank plot (knee)
            counts = np.array(sorted(frag_per_barcode.values(), reverse=True))
            fig, ax = plt.subplots(figsize=(6, 4.5))
            ax.plot(np.arange(1, counts.size + 1), counts, color="#55a868")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel("barcode rank (log)")
            ax.set_ylabel("fragments per barcode (log)")
            ax.set_title(f"{set_acc}: barcode-rank fragment counts (sampled)")
            fig.tight_layout()
            out = plot_dir / "atac_barcode_rank.png"
            fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

            # Joint plot: ATAC fragments vs RNA counts per cell (matched barcodes)
            if pc_scores is not None and idx.size > 0:
                rna_b = obs.index.astype(str).values
                rna_b_emb = rna_b[idx]
                rna_total_emb = counts_per_cell[idx]
                frag_lookup = {b: c for b, c in frag_per_barcode.items()}
                atac_emb = np.array([frag_lookup.get(str(b), 0) for b in rna_b_emb], dtype=float)
                # If almost no barcode overlap, try stripping suffix
                if (atac_emb > 0).sum() < 0.1 * rna_b_emb.size:
                    atac_emb = np.array([frag_lookup.get(str(b).split("-")[0], 0) for b in rna_b_emb], dtype=float)
                atac["matched_cells_with_fragments"] = int((atac_emb > 0).sum())
                if (atac_emb > 0).any():
                    fig, ax = plt.subplots(figsize=(6, 5))
                    ax.scatter(np.log10(rna_total_emb + 1), np.log10(atac_emb + 1),
                               s=5, alpha=0.5, edgecolors="none", color="#8172b2")
                    ax.set_xlabel("log10(RNA UMI counts + 1)")
                    ax.set_ylabel("log10(ATAC fragments + 1)")
                    ax.set_title(f"{set_acc}: per-cell RNA vs ATAC depth")
                    fig.tight_layout()
                    out = plot_dir / "joint_rna_vs_atac_depth.png"
                    fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)

                if umap_xy is not None and (atac_emb > 0).any():
                    fig, ax = plt.subplots(figsize=(6, 5))
                    sc4 = ax.scatter(umap_xy[:, 0], umap_xy[:, 1], s=5, alpha=0.7,
                                     c=np.log10(atac_emb + 1), cmap="plasma", edgecolors="none")
                    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
                    ax.set_title(f"{set_acc}: RNA UMAP coloured by log10(ATAC frag + 1)")
                    fig.colorbar(sc4, ax=ax, label="log10(ATAC frag + 1)")
                    fig.tight_layout()
                    out = plot_dir / "joint_umap_atac_overlay.png"
                    fig.savefig(out, dpi=130); plt.close(fig); plots.append(out)
        except Exception as exc:
            logging.warning("[%s] fragments QC failed: %s", set_acc, exc)
    summary["atac"] = atac

    summary["plots"] = [rel_path(p) for p in plots]
    return summary


def sample_fragments(path: Path, max_lines: int = 1_000_000):
    """Stream a 10x ATAC fragments BED.gz; return (fragment_sizes, fragments_per_barcode).

    Skips comment headers ("#"). 10x format columns: chrom start end barcode dup_count.
    """
    sizes: list[int] = []
    counts: dict[str, int] = {}
    n_data = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            sizes.append(end - start)
            barcode = parts[3]
            counts[barcode] = counts.get(barcode, 0) + 1
            n_data += 1
            if n_data >= max_lines:
                break
    return sizes, counts


# ---------------------------------------------------------------------------
# Top-level analyze + report
# ---------------------------------------------------------------------------


def analyze_all() -> Path:
    np, sp, spla, plt, pd = _import_libs()
    rows = read_selection()
    sets: dict[str, dict[str, Any]] = {}
    for r in rows:
        sets.setdefault(r["accession"], {"accession": r["accession"], "include_fragments": False})
        if r["kind"] == "fragments":
            sets[r["accession"]]["include_fragments"] = True
    items = list(sets.values())
    items.sort(key=lambda x: x["accession"])
    analyses = []
    for item in items:
        try:
            analyses.append(analyze_one_set(item, np, sp, spla, plt, pd))
        except Exception as exc:
            logging.exception("analysis failed for %s", item["accession"])
            analyses.append({"accession": item["accession"], "error": str(exc), "plots": []})
        n_plots = len(analyses[-1].get("plots", []))
        print(f"analyzed {item['accession']:18s} plots={n_plots}")

    ANALYSIS_JSON.parent.mkdir(parents=True, exist_ok=True)
    ANALYSIS_JSON.write_text(json.dumps(analyses, indent=2, default=str), encoding="utf-8")
    write_report(analyses)
    return REPORT_PATH


def write_report(analyses: list[dict[str, Any]]) -> Path:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# IGVF Portal 10x Multiome — Ten-Set Analysis",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Selection manifest: `{rel_path(SELECTION_CSV)}`  ",
        f"Download log: `{rel_path(DOWNLOAD_LOG)}`",
        "",
    ]
    for a in analyses:
        acc = a["accession"]
        lines.append(f"## {acc}")
        lines.append("")
        if "error" in a:
            lines.append(f"- ERROR: {a['error']}")
            lines.append("")
            continue
        rna = a.get("rna", {})
        lines.append(f"- RNA cells × genes: {rna.get('n_cells', '?'):,} × {rna.get('n_genes', '?'):,}")
        lines.append(
            f"- RNA QC medians: counts={rna.get('counts_per_cell_median',0):.0f}, "
            f"genes={rna.get('genes_per_cell_median',0):.0f}, %mito={rna.get('pct_mito_median',0):.2f}"
        )
        lines.append(
            f"- RNA cells kept (n_genes ≥ 200): {rna.get('n_cells_kept', '?'):,}; "
            f"embedded: {rna.get('n_cells_embedded', '?'):,}; clusters: {rna.get('n_clusters', '?')}"
        )
        ann = a.get("cell_annotations", {})
        if ann.get("label_column"):
            uniq = ann.get("unique_labels", [])[:8]
            uniq_str = ", ".join(f"{u['label']}:{u['count']}" for u in uniq)
            lines.append(
                f"- IGVF cell annotations column: `{ann['label_column']}`; "
                f"barcodes overlap = {ann.get('barcodes_overlap','?'):,}; top labels: {uniq_str}"
            )
        cc = a.get("cell_cycle", {})
        if cc.get("phase_counts"):
            pc_str = ", ".join(f"{p}={n}" for p, n in cc["phase_counts"].items())
            lines.append(
                f"- cell-cycle (Tirosh): S matched={cc.get('s_genes_matched',0)}, "
                f"G2M matched={cc.get('g2m_genes_matched',0)}, phases: {pc_str}"
            )
        atac = a.get("atac", {})
        if atac.get("fragments_present"):
            lines.append(
                f"- ATAC fragments: sampled {atac.get('fragments_sampled',0):,} fragments "
                f"across {atac.get('unique_barcodes_sampled',0):,} barcodes; "
                f"matched cells = {atac.get('matched_cells_with_fragments','?')}"
            )
        else:
            lines.append("- ATAC fragments: not downloaded for this set")
        cm = a.get("cluster_markers") or []
        if cm:
            lines.append("")
            lines.append("Top RNA cluster markers (Welch t-test, top 5 / cluster):")
            for entry in cm:
                names = ", ".join(g["gene"] for g in entry["markers"][:5])
                lines.append(f"  - c{entry['cluster']}: {names}")
        lines.append("")
        lines.append("Plots:")
        for p in a.get("plots", []):
            try:
                rel = abs_from_rel(p).resolve().relative_to(REPORT_PATH.parent.resolve()).as_posix()
            except ValueError:
                rel = p
            lines.append(f"- ![{Path(p).stem}]({rel})")
        lines.append("")
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote report: %s", REPORT_PATH)
    return REPORT_PATH


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sel = sub.add_parser("select", help="Pick 10 small AnalysisSets and write the selection manifest.")
    sel.add_argument("--n-sets", type=int, default=10)
    sel.add_argument("--n-with-fragments", type=int, default=3)

    dl = sub.add_parser("download", help="Download the selected files.")
    dl.add_argument("--force", action="store_true")

    sub.add_parser("analyze", help="Analyze the downloaded files and write report.")
    sub.add_parser("run", help="select -> download -> analyze.")

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "select":
        items = select_sets(n_sets=args.n_sets, n_with_fragments=args.n_with_fragments)
        write_selection(items)
        for it in items:
            tag = " +frag" if it["include_fragments"] else ""
            small_mb = it["small_size"] / 1024 / 1024
            print(f"  {it['accession']:18s} small={small_mb:6.1f}MB{tag}")
        return 0
    if args.command == "download":
        download_selection(force=args.force)
        return 0
    if args.command == "analyze":
        analyze_all()
        return 0
    if args.command == "run":
        items = select_sets()
        write_selection(items)
        download_selection(force=False)
        analyze_all()
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
