"""Single-cell analysis skill — UMAP / t-SNE / clustering / markers.

Closes the gap where IGVFagent could discover and download single-cell
datasets but had to hand off to "use Scanpy / Seurat" for the actual
visualization. This skill drives the standard Scanpy workflow end-to-
end with one CLI:

    raw counts  →  QC + filter  →  log-normalize + HVG
                →  PCA  →  k-NN graph + UMAP + t-SNE
                →  Leiden clustering  →  rank-genes-groups (markers)
                →  publication-style figures + report

Input formats auto-detected:
    .h5ad                   anndata native
    .h5                     10x Genomics CellRanger HDF5
    .mtx / .mtx.gz          10x sparse matrix (with companion barcodes.tsv
                            and features.tsv)
    .csv / .tsv             genes × cells (or cells × genes — set --transpose)

Outputs land under ``Docs/SingleCell/<ts>_<label>/`` with:
    processed.h5ad          full anndata snapshot (resumable)
    Plots/qc_violins.png
    Plots/pca_variance.png
    Plots/umap_clusters.png
    Plots/umap_sample.png   (if sample col present)
    Plots/tsne_clusters.png
    Plots/umap_top_markers.png
    Plots/heatmap_top_markers.png
    markers.csv
    qc_summary.json
    report.md

Subcommands
-----------
    qc                  Load + QC + write violins, basic stats
    normalize           Log-normalize + select HVGs
    pca                 PCA on HVGs
    umap                Build k-NN graph + UMAP + save figure
    tsne                t-SNE + save figure
    cluster             Leiden (or Louvain fallback)
    markers             rank-genes-groups + Top-N CSV
    pipeline            Full end-to-end (QC→norm→PCA→UMAP→tSNE→cluster→markers)
    plot-embedding      Re-plot an existing UMAP/tSNE colored by gene or obs
    write-playbook      Write Docs/Skills/SINGLECELL_ANALYSIS_SKILLS.md
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "SingleCell"
DOCS_DIR = ROOT / "Docs" / "SingleCell"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "SINGLECELL_ANALYSIS_SKILLS.md"

logger = logging.getLogger("singlecell_analysis")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in (s or "run"))


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"sc_analyze_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def mkdirs() -> None:
    for d in (DATA_DIR, DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def _new_run_dir(label: str) -> Path:
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}"
    (out / "Plots").mkdir(parents=True, exist_ok=True)
    return out


def _scanpy():
    """Soft import scanpy + matplotlib. Raise an actionable error on missing
    deps so the agent's failure message tells the user exactly what to do."""
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import scanpy as sc  # type: ignore
        import anndata as ad  # type: ignore
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "single-cell analysis requires scanpy + anndata + numpy + "
            "pandas + matplotlib in this environment. Install with:\n"
            "  pip install scanpy 'anndata>=0.10' umap-learn leidenalg "
            "python-igraph matplotlib\n"
            f"Original error: {exc}"
        )
    return sc, ad, np, pd


# ---------------------------------------------------------------------------
# Loader — auto-detect format
# ---------------------------------------------------------------------------


def load_counts(path: Path, *, transpose: bool = False):
    """Load any of .h5ad / 10x .h5 / .mtx / .csv / .tsv  -> AnnData.

    For .mtx input, expects companion ``barcodes.tsv`` and
    ``features.tsv`` (or ``genes.tsv``) in the same directory.
    """
    sc, ad, np, pd = _scanpy()
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = "".join(path.suffixes).lower()

    if suffix.endswith(".h5ad"):
        logger.info("Reading anndata: %s", path)
        adata = sc.read_h5ad(path)
    elif suffix.endswith(".h5"):
        logger.info("Reading 10x HDF5: %s", path)
        adata = sc.read_10x_h5(str(path))
        adata.var_names_make_unique()
    elif suffix.endswith(".mtx") or suffix.endswith(".mtx.gz"):
        logger.info("Reading 10x mtx (dir=%s)", path.parent)
        adata = sc.read_mtx(str(path)).T  # mtx is usually genes × cells
        # Look for barcodes + features in the same dir
        for bn in ("barcodes.tsv", "barcodes.tsv.gz"):
            bp = path.parent / bn
            if bp.exists():
                opener = gzip.open if bp.suffix == ".gz" else open
                with opener(bp, "rt") as f:
                    adata.obs_names = [ln.strip() for ln in f if ln.strip()]
                break
        for fn in ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"):
            fp = path.parent / fn
            if fp.exists():
                opener = gzip.open if fp.suffix == ".gz" else open
                with opener(fp, "rt") as f:
                    names = [ln.strip().split("\t")[1] if "\t" in ln else
                                ln.strip() for ln in f if ln.strip()]
                adata.var_names = names
                adata.var_names_make_unique()
                break
    elif suffix.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        sep = "," if ".csv" in suffix else "\t"
        logger.info("Reading text matrix (sep=%r): %s", sep, path)
        df = pd.read_csv(path, sep=sep, index_col=0,
                          compression="infer")
        if transpose:
            df = df.T
        # Heuristic: cells should be rows; if rows > cols, assume genes are rows
        if df.shape[0] < df.shape[1]:
            df = df.T
        adata = ad.AnnData(df.values.astype("float32"))
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)
    else:
        raise ValueError(f"Unsupported input format: {suffix}")

    logger.info("Loaded AnnData: n_obs=%d  n_vars=%d", adata.n_obs, adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------


def qc(adata, *, min_genes: int = 200, min_cells: int = 3,
        max_mito: float = 20.0, mito_prefix: str = "MT-",
        out: Path) -> dict:
    sc, ad, np, pd = _scanpy()
    import matplotlib.pyplot as plt  # type: ignore

    n_before = adata.n_obs
    g_before = adata.n_vars
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # Mito QC
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
                                  percent_top=None, log1p=False, inplace=True)
    n_kept = (adata.obs["pct_counts_mt"] < max_mito).sum()
    adata._inplace_subset_obs(adata.obs["pct_counts_mt"] < max_mito)
    n_after = adata.n_obs

    summary = {
        "n_obs_before": n_before,
        "n_obs_after_filter_cells": n_kept,
        "n_obs_after_mito_filter": n_after,
        "n_vars_before": g_before,
        "n_vars_after": adata.n_vars,
        "median_n_genes": float(np.median(adata.obs["n_genes_by_counts"])),
        "median_total_counts": float(np.median(adata.obs["total_counts"])),
        "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"])),
        "min_genes": min_genes,
        "min_cells": min_cells,
        "max_mito": max_mito,
        "mito_genes_detected": int(adata.var["mt"].sum()),
    }

    # QC violins
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    for ax, col, color, ttl in zip(
        axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["#1F77B4", "#2CA02C", "#D62728"],
        ["Genes / cell", "Counts / cell", "% mito"]):
        vals = adata.obs[col].values
        ax.violinplot([vals], showmedians=True, widths=0.7)
        ax.scatter(np.ones_like(vals) + np.random.uniform(-0.15, 0.15, len(vals)),
                    vals, s=2, alpha=0.25, color=color)
        ax.set_title(ttl); ax.set_xticks([]); ax.grid(alpha=0.25, linestyle=":")
    fig.suptitle(f"QC violins  ·  n_cells = {adata.n_obs:,}  ·  n_genes = {adata.n_vars:,}",
                  fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "Plots" / "qc_violins.png", dpi=150)
    plt.close(fig)
    return summary


# ---------------------------------------------------------------------------
# Normalize + HVG
# ---------------------------------------------------------------------------


def normalize(adata, *, n_hvg: int = 2000, target_sum: float = 1e4):
    sc, _, _, _ = _scanpy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    # Keep raw for marker DE later
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    return adata


# ---------------------------------------------------------------------------
# PCA + neighbors + UMAP + tSNE
# ---------------------------------------------------------------------------


def pca(adata, *, n_comps: int = 50, out: Path):
    sc, _, np, _ = _scanpy()
    import matplotlib.pyplot as plt  # type: ignore
    n_comps_eff = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps_eff, svd_solver="arpack")
    # Variance elbow
    var = adata.uns["pca"]["variance_ratio"]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(range(1, len(var) + 1), var, marker="o", ms=3, color="#08519C")
    ax.set_xlabel("PC"); ax.set_ylabel("Variance explained")
    ax.set_title(f"PCA variance ratio (top {n_comps_eff} PCs)")
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out / "Plots" / "pca_variance.png", dpi=150)
    plt.close(fig)
    return adata


def umap(adata, *, n_neighbors: int = 15, n_pcs: int = 40,
         out: Optional[Path] = None):
    sc, _, _, _ = _scanpy()
    n_pcs_eff = min(n_pcs, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_eff)
    sc.tl.umap(adata)
    return adata


def tsne(adata, *, n_pcs: int = 30):
    sc, _, _, _ = _scanpy()
    n_pcs_eff = min(n_pcs, adata.obsm["X_pca"].shape[1])
    sc.tl.tsne(adata, n_pcs=n_pcs_eff)
    return adata


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------


def cluster(adata, *, resolution: float = 1.0):
    sc, _, _, _ = _scanpy()
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                      flavor="igraph", n_iterations=2, directed=False)
    except Exception as e:
        logger.warning("Leiden failed (%s) — falling back to Louvain", e)
        sc.tl.louvain(adata, resolution=resolution, key_added="leiden")
    return adata


# ---------------------------------------------------------------------------
# Markers + figures
# ---------------------------------------------------------------------------


def markers(adata, *, n_top: int = 25, group_key: str = "leiden",
             out: Path):
    sc, _, np, pd = _scanpy()
    import matplotlib.pyplot as plt  # type: ignore
    sc.tl.rank_genes_groups(adata, groupby=group_key, method="wilcoxon",
                              use_raw=True)
    # Build a long DataFrame
    res = adata.uns["rank_genes_groups"]
    groups = res["names"].dtype.names
    rows = []
    for g in groups:
        for i in range(min(n_top, len(res["names"][g]))):
            rows.append({
                "group": g,
                "rank": i + 1,
                "gene": res["names"][g][i],
                "logfc": float(res["logfoldchanges"][g][i]),
                "pval": float(res["pvals"][g][i]),
                "padj": float(res["pvals_adj"][g][i]),
                "score": float(res["scores"][g][i]),
            })
    df = pd.DataFrame(rows)
    df.to_csv(out / "markers.csv", index=False)

    # Heatmap of top-3 markers per cluster
    top3 = (df[df["rank"] <= 3].groupby("group")["gene"]
              .apply(list).reset_index())
    genes = []
    for _, r in top3.iterrows():
        for g in r["gene"]:
            if g not in genes:
                genes.append(g)
    if genes:
        try:
            sc.pl.heatmap(adata, var_names=genes, groupby=group_key,
                           use_raw=True, swap_axes=False, show=False,
                           save=False, dendrogram=False, standard_scale="var")
            import matplotlib.pyplot as plt  # noqa
            plt.gcf().savefig(out / "Plots" / "heatmap_top_markers.png",
                                dpi=150, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            logger.warning("Marker heatmap skipped: %s", e)
    return df


def plot_umap(adata, *, color: "list[str]", out: Path, suffix: str = "",
              ncols: int = 2, palette: Optional[str] = None) -> Optional[Path]:
    sc, _, _, _ = _scanpy()
    import matplotlib.pyplot as plt  # type: ignore
    if "X_umap" not in adata.obsm:
        return None
    # Subset color list to keys that exist in obs or var_names
    keep = []
    for c in color:
        if c in adata.obs.columns or c in adata.var_names \
            or (adata.raw is not None and c in adata.raw.var_names):
            keep.append(c)
    if not keep:
        return None
    sc.pl.umap(adata, color=keep, ncols=ncols, show=False,
                use_raw=True, frameon=False, palette=palette)
    fn = out / "Plots" / f"umap{suffix}.png"
    plt.gcf().savefig(fn, dpi=150, bbox_inches="tight")
    plt.close("all")
    return fn


def plot_tsne(adata, *, color: "list[str]", out: Path, suffix: str = "",
               ncols: int = 2) -> Optional[Path]:
    sc, _, _, _ = _scanpy()
    import matplotlib.pyplot as plt  # type: ignore
    if "X_tsne" not in adata.obsm:
        return None
    keep = [c for c in color
             if c in adata.obs.columns or c in adata.var_names
             or (adata.raw is not None and c in adata.raw.var_names)]
    if not keep:
        return None
    sc.pl.tsne(adata, color=keep, ncols=ncols, show=False,
                use_raw=True, frameon=False)
    fn = out / "Plots" / f"tsne{suffix}.png"
    plt.gcf().savefig(fn, dpi=150, bbox_inches="tight")
    plt.close("all")
    return fn


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def run_pipeline(*, input_path: Path, label: str,
                  min_genes: int = 200, min_cells: int = 3,
                  max_mito: float = 20.0, mito_prefix: str = "MT-",
                  n_hvg: int = 2000, n_pcs: int = 50,
                  n_neighbors: int = 15, resolution: float = 1.0,
                  do_tsne: bool = True, n_markers: int = 25,
                  sample_col: Optional[str] = None,
                  highlight_genes: "Optional[list[str]]" = None,
                  transpose: bool = False) -> Path:
    sc, ad, np, pd = _scanpy()
    mkdirs()
    out = _new_run_dir(label)
    summary: "dict[str, Any]" = {"label": label, "input": str(input_path),
                                    "started_at": time.strftime("%FT%T")}

    logger.info("==> Load")
    adata = load_counts(input_path, transpose=transpose)
    summary["loaded_n_obs"] = adata.n_obs
    summary["loaded_n_vars"] = adata.n_vars

    logger.info("==> QC + filter")
    summary["qc"] = qc(adata, min_genes=min_genes, min_cells=min_cells,
                         max_mito=max_mito, mito_prefix=mito_prefix, out=out)

    logger.info("==> Normalize + HVG (top %d)", n_hvg)
    adata = normalize(adata, n_hvg=n_hvg)

    logger.info("==> PCA (n_comps=%d)", n_pcs)
    adata = pca(adata, n_comps=n_pcs, out=out)

    logger.info("==> Neighbors + UMAP")
    adata = umap(adata, n_neighbors=n_neighbors, n_pcs=min(40, n_pcs), out=out)

    if do_tsne:
        logger.info("==> t-SNE")
        try:
            adata = tsne(adata, n_pcs=min(30, n_pcs))
        except Exception as e:
            logger.warning("t-SNE failed: %s", e)

    logger.info("==> Leiden clustering (res=%.2f)", resolution)
    adata = cluster(adata, resolution=resolution)
    summary["n_clusters"] = int(adata.obs["leiden"].nunique())

    logger.info("==> Marker genes (top %d / cluster)", n_markers)
    markers_df = markers(adata, n_top=n_markers, group_key="leiden", out=out)
    summary["markers_total"] = len(markers_df)

    # UMAP / tSNE figures
    color_list = ["leiden"]
    if sample_col and sample_col in adata.obs.columns:
        color_list.append(sample_col)
    if highlight_genes:
        color_list.extend(highlight_genes)
    plot_umap(adata, color=["leiden"], out=out, suffix="_clusters")
    if sample_col and sample_col in adata.obs.columns:
        plot_umap(adata, color=[sample_col], out=out, suffix="_sample")
    if highlight_genes:
        plot_umap(adata, color=highlight_genes, out=out,
                    suffix="_top_markers", ncols=2)
    if do_tsne and "X_tsne" in adata.obsm:
        plot_tsne(adata, color=["leiden"], out=out, suffix="_clusters")

    # Auto-overlay top 4 marker genes (cluster-specific) on UMAP
    if not highlight_genes and len(markers_df) > 0:
        top_per_cluster = (markers_df[markers_df["rank"] == 1]
                              .sort_values("score", ascending=False)
                              .head(4)["gene"].tolist())
        if top_per_cluster:
            plot_umap(adata, color=top_per_cluster, out=out,
                        suffix="_auto_top_markers", ncols=2)
            summary["auto_top_markers"] = top_per_cluster

    # Save processed h5ad for downstream re-use
    h5ad_path = out / "processed.h5ad"
    adata.write_h5ad(h5ad_path, compression="gzip")
    summary["processed_h5ad"] = str(h5ad_path)
    summary["ended_at"] = time.strftime("%FT%T")
    (out / "qc_summary.json").write_text(json.dumps(summary, indent=2,
                                                       default=str))

    # Markdown report
    md = ["# Single-cell analysis — " + label, "",
          f"- Input: `{input_path}`",
          f"- Loaded: **{summary['loaded_n_obs']:,} cells × {summary['loaded_n_vars']:,} genes**",
          f"- After filter: **{summary['qc']['n_obs_after_mito_filter']:,} cells × {summary['qc']['n_vars_after']:,} genes**",
          f"- HVGs retained: **{n_hvg}**",
          f"- PCs computed: **{min(n_pcs, summary['qc']['n_vars_after'])}**",
          f"- Clusters found (Leiden res={resolution}): **{summary['n_clusters']}**",
          f"- Markers (top {n_markers}/cluster): **{summary['markers_total']:,}** rows",
          "",
          "## Plots",
          "- `Plots/qc_violins.png` — QC violins (genes / counts / mito %)",
          "- `Plots/pca_variance.png` — PCA scree",
          "- `Plots/umap_clusters.png` — UMAP coloured by Leiden cluster"]
    if sample_col:
        md.append(f"- `Plots/umap_sample.png` — UMAP coloured by `{sample_col}`")
    if "auto_top_markers" in summary:
        md.append("- `Plots/umap_auto_top_markers.png` — UMAP overlaid with "
                  f"top markers ({', '.join(summary['auto_top_markers'])})")
    if do_tsne and "X_tsne" in adata.obsm:
        md.append("- `Plots/tsne_clusters.png` — t-SNE coloured by Leiden cluster")
    md.append("- `Plots/heatmap_top_markers.png` — top-3 marker heatmap")
    md += ["", f"## Saved AnnData", f"- `{h5ad_path}`  (compressed h5ad — "
            f"resumable by passing this file as `--input` to any sub-step)"]
    (out / "report.md").write_text("\n".join(md) + "\n")
    logger.info("Wrote %s", out / "report.md")
    return out


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_qc(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("qc_" + (args.label or "run"))
    out = _new_run_dir(args.label or "qc")
    adata = load_counts(Path(args.input), transpose=args.transpose)
    summary = qc(adata, min_genes=args.min_genes, min_cells=args.min_cells,
                  max_mito=args.max_mito, mito_prefix=args.mito_prefix, out=out)
    (out / "qc_summary.json").write_text(json.dumps(summary, indent=2,
                                                       default=str))
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out}")
    print(json.dumps(summary, indent=2, default=str))
    return 0


def cmd_normalize(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("normalize")
    out = _new_run_dir(args.label or "normalize")
    adata = load_counts(Path(args.input))
    adata = normalize(adata, n_hvg=args.n_hvg)
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out / 'processed.h5ad'}")
    return 0


def cmd_pca(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("pca")
    out = _new_run_dir(args.label or "pca")
    adata = load_counts(Path(args.input))
    adata = pca(adata, n_comps=args.n_pcs, out=out)
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out}")
    return 0


def cmd_umap(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("umap")
    out = _new_run_dir(args.label or "umap")
    adata = load_counts(Path(args.input))
    if "X_pca" not in adata.obsm:
        adata = pca(adata, n_comps=args.n_pcs, out=out)
    adata = umap(adata, n_neighbors=args.n_neighbors,
                    n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]),
                    out=out)
    color = ["leiden"] if "leiden" in adata.obs.columns else []
    if not color and "sample" in adata.obs.columns:
        color = ["sample"]
    if not color:
        color = [adata.obs.columns[0]] if len(adata.obs.columns) else []
    plot_umap(adata, color=color or ["n_genes_by_counts"], out=out,
                suffix="_" + (color[0] if color else "qc"))
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out}")
    return 0


def cmd_tsne(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("tsne")
    out = _new_run_dir(args.label or "tsne")
    adata = load_counts(Path(args.input))
    if "X_pca" not in adata.obsm:
        adata = pca(adata, n_comps=args.n_pcs, out=out)
    adata = tsne(adata, n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]))
    plot_tsne(adata, color=(["leiden"] if "leiden" in adata.obs.columns
                              else [adata.obs.columns[0]]),
                out=out, suffix="_clusters")
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out}")
    return 0


def cmd_cluster(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("cluster")
    out = _new_run_dir(args.label or "cluster")
    adata = load_counts(Path(args.input))
    if "neighbors" not in adata.uns:
        if "X_pca" not in adata.obsm:
            adata = pca(adata, n_comps=args.n_pcs, out=out)
        adata = umap(adata, n_neighbors=args.n_neighbors,
                        n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]),
                        out=out)
    adata = cluster(adata, resolution=args.resolution)
    plot_umap(adata, color=["leiden"], out=out, suffix="_clusters")
    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    print(f"Output: {out}  ·  clusters={adata.obs['leiden'].nunique()}")
    return 0


def cmd_markers(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("markers")
    out = _new_run_dir(args.label or "markers")
    adata = load_counts(Path(args.input))
    if args.group_key not in adata.obs.columns:
        raise SystemExit(f"--group-key {args.group_key!r} not in obs. "
                         f"Run `cluster` first.")
    df = markers(adata, n_top=args.n_top, group_key=args.group_key, out=out)
    print(f"Output: {out}/markers.csv  ·  rows={len(df)}")
    return 0


def cmd_pipeline(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("pipeline_" + (args.label or "run"))
    highlight = [g.strip() for g in (args.highlight_genes or "").split(",")
                  if g.strip()]
    out = run_pipeline(
        input_path=Path(args.input),
        label=args.label or "pipeline",
        min_genes=args.min_genes, min_cells=args.min_cells,
        max_mito=args.max_mito, mito_prefix=args.mito_prefix,
        n_hvg=args.n_hvg, n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors, resolution=args.resolution,
        do_tsne=not args.skip_tsne, n_markers=args.n_markers,
        sample_col=args.sample_col,
        highlight_genes=highlight or None,
        transpose=args.transpose,
    )
    print(f"Output: {out}")
    return 0


def cmd_plot_embedding(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("plot_embedding")
    out = _new_run_dir(args.label or "plot")
    adata = load_counts(Path(args.input))
    color = [c.strip() for c in args.color.split(",") if c.strip()]
    paths = []
    if args.embedding in ("umap", "both"):
        p = plot_umap(adata, color=color, out=out,
                       suffix="_" + "_".join(color)[:60], ncols=args.ncols)
        if p: paths.append(str(p))
    if args.embedding in ("tsne", "both"):
        p = plot_tsne(adata, color=color, out=out,
                       suffix="_" + "_".join(color)[:60], ncols=args.ncols)
        if p: paths.append(str(p))
    print("\n".join(paths) if paths else "(no embedding present in input)")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: Single-cell analysis (UMAP / t-SNE / Leiden / markers)

End-to-end Scanpy-driven single-cell workflow. Closes the gap where
IGVFagent could discover and download counts matrices but had to hand
off to "run Scanpy or Seurat" for the actual analysis.

## Subcommands

### qc
```
igvfagent sc-analyze qc --input counts.h5ad --label demo \\
    --min-genes 200 --min-cells 3 --max-mito 20
```
Loads any of `.h5ad` / 10x `.h5` / `.mtx` / `.csv` / `.tsv`, filters
cells & genes, computes mito %, drops high-mito cells, writes QC
violins, and saves the cleaned `processed.h5ad`.

### normalize
```
igvfagent sc-analyze normalize --input processed.h5ad --n-hvg 2000
```
Total-count normalize, log1p, select top-N highly-variable genes,
scale to unit variance. Saves raw counts in `.raw` for downstream
marker DE.

### pca / umap / tsne / cluster / markers
```
igvfagent sc-analyze pca       --input processed.h5ad --n-pcs 50
igvfagent sc-analyze umap      --input processed.h5ad --n-neighbors 15
igvfagent sc-analyze tsne      --input processed.h5ad --n-pcs 30
igvfagent sc-analyze cluster   --input processed.h5ad --resolution 1.0
igvfagent sc-analyze markers   --input processed.h5ad --n-top 25
```
Each step accepts the `processed.h5ad` from any previous step and
saves a new `processed.h5ad` carrying the additional fields
(`X_pca`, `X_umap`, `X_tsne`, `leiden`, `rank_genes_groups`).

### pipeline (recommended one-shot)
```
igvfagent sc-analyze pipeline --input counts.h5ad --label k562_demo \\
    --min-genes 200 --min-cells 3 --max-mito 20 \\
    --n-hvg 2000 --n-pcs 50 --resolution 1.0 \\
    --sample-col sample \\
    --highlight-genes APOE,TREM2,LDLR,GFAP
```
Drives QC → normalize → PCA → neighbors → UMAP → t-SNE → Leiden →
marker DE → figures → markdown report in one go. Saves the
`processed.h5ad` so any later step can resume.

### plot-embedding
```
igvfagent sc-analyze plot-embedding --input processed.h5ad \\
    --embedding umap --color leiden,APOE,TREM2 --label apoe_view
```
Re-render UMAP or t-SNE coloured by any combination of obs columns
and genes. Useful for asking the agent "show me APOE expression
overlaid on the UMAP we just made."

## Inputs

| Format | Detection rule |
|---|---|
| AnnData `.h5ad` | suffix `.h5ad` |
| 10x CellRanger | suffix `.h5` |
| 10x sparse | suffix `.mtx` / `.mtx.gz` (looks for `barcodes.tsv` + `features.tsv` in same dir) |
| Text matrix | `.csv` / `.tsv` / `.txt` (and `.gz` variants) — auto-orients to cells×genes |

## Outputs

```
Docs/SingleCell/<ts>_<label>/
    processed.h5ad
    qc_summary.json
    markers.csv
    Plots/
        qc_violins.png
        pca_variance.png
        umap_clusters.png
        umap_sample.png            (if --sample-col given)
        umap_auto_top_markers.png  (auto-picks 4 cluster markers)
        tsne_clusters.png
        heatmap_top_markers.png
    report.md
```

## Cross-skill chaining

- `singlecell` / `multiome` / `splitseq` — discover datasets, then
  feed their downloaded counts straight into `sc-analyze pipeline`.
- `geo download` → `sc-analyze pipeline` — wide net retrieval from
  GEO, then end-to-end analysis.
- `kg` / `proteomics kg-visualize` — after `sc-analyze markers`,
  drop the top marker genes into the KG visualizer to see their
  protein-protein interaction context.

## Dependencies

`scanpy>=1.10`, `anndata>=0.10`, `umap-learn`, `scikit-learn`,
`leidenalg`, `python-igraph`, `matplotlib`. Installed in the
project venv. If you're using an external Python, install with:

```
pip install scanpy 'anndata>=0.10' umap-learn leidenalg python-igraph
```
"""


def cmd_write_playbook(_a) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_TEXT)
    print(f"Wrote: {PLAYBOOK_PATH}")
    return 0


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------


def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(
        prog="sc-analyze",
        description="Single-cell analysis: QC, PCA, UMAP, t-SNE, "
                    "Leiden clustering, marker-gene DE, publication "
                    "figures.")
    sub = p.add_subparsers(dest="cmd")

    def _common_io(s):
        s.add_argument("--input", required=True,
                        help="Input counts (.h5ad / .h5 / .mtx / .csv / .tsv).")
        s.add_argument("--label", default=None)
        s.add_argument("--transpose", action="store_true",
                        help="Force transpose on .csv/.tsv input.")

    s = sub.add_parser("qc", help="Load + QC + write violins.")
    _common_io(s)
    s.add_argument("--min-genes", type=int, default=200)
    s.add_argument("--min-cells", type=int, default=3)
    s.add_argument("--max-mito", type=float, default=20.0)
    s.add_argument("--mito-prefix", default="MT-")
    s.set_defaults(func=cmd_qc)

    s = sub.add_parser("normalize", help="Total-count normalize + HVG.")
    _common_io(s)
    s.add_argument("--n-hvg", type=int, default=2000)
    s.set_defaults(func=cmd_normalize)

    s = sub.add_parser("pca", help="PCA on HVGs.")
    _common_io(s)
    s.add_argument("--n-pcs", type=int, default=50)
    s.set_defaults(func=cmd_pca)

    s = sub.add_parser("umap", help="k-NN + UMAP + figure.")
    _common_io(s)
    s.add_argument("--n-pcs", type=int, default=40)
    s.add_argument("--n-neighbors", type=int, default=15)
    s.set_defaults(func=cmd_umap)

    s = sub.add_parser("tsne", help="t-SNE + figure.")
    _common_io(s)
    s.add_argument("--n-pcs", type=int, default=30)
    s.set_defaults(func=cmd_tsne)

    s = sub.add_parser("cluster", help="Leiden clustering on the k-NN graph.")
    _common_io(s)
    s.add_argument("--resolution", type=float, default=1.0)
    s.add_argument("--n-pcs", type=int, default=40)
    s.add_argument("--n-neighbors", type=int, default=15)
    s.set_defaults(func=cmd_cluster)

    s = sub.add_parser("markers", help="Rank genes per cluster (Wilcoxon).")
    _common_io(s)
    s.add_argument("--n-top", type=int, default=25)
    s.add_argument("--group-key", default="leiden")
    s.set_defaults(func=cmd_markers)

    s = sub.add_parser("pipeline",
                        help="Full QC→norm→PCA→UMAP→t-SNE→Leiden→markers.")
    _common_io(s)
    s.add_argument("--min-genes", type=int, default=200)
    s.add_argument("--min-cells", type=int, default=3)
    s.add_argument("--max-mito", type=float, default=20.0)
    s.add_argument("--mito-prefix", default="MT-")
    s.add_argument("--n-hvg", type=int, default=2000)
    s.add_argument("--n-pcs", type=int, default=50)
    s.add_argument("--n-neighbors", type=int, default=15)
    s.add_argument("--resolution", type=float, default=1.0)
    s.add_argument("--skip-tsne", action="store_true")
    s.add_argument("--n-markers", type=int, default=25)
    s.add_argument("--sample-col", default=None,
                    help="Obs column to color UMAP by (e.g. sample, batch).")
    s.add_argument("--highlight-genes", default=None,
                    help="Comma-separated gene symbols to overlay on UMAP.")
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("plot-embedding",
                        help="Re-plot UMAP/t-SNE colored by gene or obs.")
    _common_io(s)
    s.add_argument("--embedding", choices=["umap", "tsne", "both"],
                    default="umap")
    s.add_argument("--color", required=True,
                    help="Comma-separated obs columns and/or gene symbols.")
    s.add_argument("--ncols", type=int, default=2)
    s.set_defaults(func=cmd_plot_embedding)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/SINGLECELL_ANALYSIS_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
