#!/usr/bin/env python3
"""IGVF TF Perturb-seq analysis stages, rewritten from the consortium package.

The IGVF TF Perturb-seq pillar project (https://github.com/IGVF/tf_perturb_seq,
package ``src/tf_perturb_seq``, no licence file) runs five stages on every
dataset: portal staging, the Nextflow CRISPR pipeline on GCP, local QC,
energy distance, and gene programs (cNMF). Stages 1 and 2 are infrastructure
(GCS uploads, Nextflow on Google Batch) and stay upstream. Everything the
package computes locally after the pipeline has produced
``inference_mudata.h5mu`` is rewritten here so it runs inside IGVFagent:

  qc-gene        per-cell UMIs / genes / mito% summary, overall and per batch
  qc-guide       guide UMIs, guides per cell, cells per guide, per-guide capture
  qc-target      intended-target knockdown: strong / significant fractions and
                 the AUROC / AUPRC of targeting vs non-targeting guides
  calibrate      empirical p-values for PerTurbo per-element results against
                 the non-targeting null (eCDF or t-fit), BH on targeting tests,
                 cis / direct-target annotation, the four calibrated tables
                 that the WG3 disease/GWAS analysis (`tf-perturb run`) consumes
  pathways       per-element Fisher over-representation of significant DEGs
                 in GMT gene sets, BH per element
  filter-cells   drop cells carrying more than N assigned guides
  filter-guides  drop outlier guides (BH on the energy-distance outlier tables)
                 and every cell that received one
  edist-prep     normalise / log1p / scale / PCA(50), guide -> cells map and the
                 annotation table the energy-distance steps read
  edist-filter   Step 1 of the energy-distance pipeline: DISCO permutation test
                 among sibling guides, hypergeometric outlier ranking, K-means
                 outliers among non-targeting guides
  edist          Step 2: per-target energy distance vs random non-targeting
                 backgrounds with permutation p-values (pval_edist_full.csv)
  edist-summary  cross-dataset roll-up of pval_edist_full.csv files
  edist-validate schema / value checks on an energy-distance output folder
  cnmf-export    inference_mudata.h5mu -> PerturbNMF single-modality .h5ad
  cnmf-validate  file / shape / value checks on a torch-cNMF run folder
  pipeline-summary  one row per dataset from the QC metric tables

Provenance and relationship, module by module:

  * QC, calibration, pathways, the filters, cNMF export, the validators and
    the cross-dataset summaries are PORTS of ``src/tf_perturb_seq`` (source
    consulted, logic rewritten; column names and output file names kept so
    downstream steps and the DACC schemas line up).
  * The energy-distance steps are a CLEAN-ROOM reimplementation of
    https://github.com/Chikara-Takeuchi/energy_dist_pipeline (no licence
    file; GPU torch code). The statistic is the same -- for cell sets A and
    B in PCA space, E = 2 mean||a-b||^2 - mean||a-a'||^2 - mean||b-b'||^2 --
    and so is the design: `num_of_bg` random non-targeting backgrounds of
    `non_target_pick` cells, `permute_per_bg` label permutations each,
    p = fraction of permuted E >= observed, then distance_mean / pval_mean,
    pval_mean_log = -log10(pval_mean + 0.1 / (permute_per_bg * num_of_bg)).
    It runs on numpy/scipy, so it is CPU-only and meant for the sizes a
    laptop or the hosted instance can hold; the upstream container is the
    right tool for 300k-cell datasets on a GPU node.
  * cNMF / PerturbNMF factorisation itself (torch-cNMF, GPU) and the Nextflow
    CRISPR pipeline are not reimplemented; `cnmf-export` prepares their input
    and `cnmf-validate` checks their output.

Where this rewrite differs, stated so nobody mistakes the two:

  * `cnmf-export` selects highly variable genes with scanpy's Seurat-flavour
    dispersion on normalised counts, not cnmf's `get_highvar_genes_sparse`
    on TPM; the cell filter it drives (drop cells with zero counts in the HVG
    set) is the same rule. Pass `--no-hvg-filter` to skip it.
  * BH adjustment is computed here (statsmodels' `fdr_bh` procedure) so
    statsmodels is not required. `t-fit` calibration needs scipy.
  * `edist-filter` keeps upstream's tests (DISCO F-statistic with label
    permutations; hypergeometric right tail on how often a guide's pairs fall
    in the top half of intra-target distances; K-means k=2 on the
    non-targeting distance matrix) and its output columns (`pval_outlier`),
    with the subset-combination scheme for targets with more than
    `threshold_gRNA_num` guides simplified to one-vs-rest comparisons.

mudata is optional throughout: an .h5mu is HDF5 with AnnData groups under
``mod/``, read through h5py + anndata when the package is absent; writing a
filtered .h5mu needs mudata, otherwise one .h5ad per modality is written.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import os
import pickle
import re
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
OUT_ROOT = ROOT / "Docs" / "TFPerturbSeq"
LOG_DIR = ROOT / "Docs" / "Logs"

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
log = logging.getLogger("tf_perturb_stages")


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"tf_perturb_stages_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def run_dir(label: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _np():
    import numpy as np  # type: ignore
    return np


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("tf-perturb stages need pandas/numpy/scipy: pip install 'igvfagent[analysis]'") from e


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    return path


def bh_adjust(pvals) -> Any:
    """Benjamini-Hochberg step-up (statsmodels fdr_bh); NaNs stay NaN."""
    np = _np()
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def column_stats(series) -> dict:
    return {"median": float(series.median()), "mean": float(series.mean()), "std": float(series.std()),
            "min": float(series.min()), "max": float(series.max()),
            "q25": float(series.quantile(0.25)), "q75": float(series.quantile(0.75))}


def _to_dense_col(mat, j: int):
    np = _np()
    col = mat[:, j]
    return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()


def _row_sums(mat):
    np = _np()
    return np.asarray(mat.sum(axis=1)).ravel()


def _col_sums(mat):
    np = _np()
    return np.asarray(mat.sum(axis=0)).ravel()


# ---------------------------------------------------------------------------
# MuData access with or without the mudata package
# ---------------------------------------------------------------------------

class Modalities:
    """The two AnnData modalities of an inference_mudata.h5mu plus top-level obs/uns."""

    def __init__(self, gene, guide, obs, uns: dict, source: Path, backend: str) -> None:
        self.gene, self.guide, self.obs, self.uns, self.source, self.backend = gene, guide, obs, uns, source, backend


def read_mudata(path: Path, need_uns: bool = True) -> Modalities:
    try:
        import mudata as md  # type: ignore
        m = md.read_h5mu(str(path))
        gene_key = "gene" if "gene" in m.mod else ("rna" if "rna" in m.mod else list(m.mod)[0])
        guide_key = "guide" if "guide" in m.mod else ("gRNA" if "gRNA" in m.mod else list(m.mod)[-1])
        return Modalities(m.mod[gene_key], m.mod[guide_key], m.obs, dict(m.uns), path, "mudata")
    except ImportError:
        pass
    import h5py  # type: ignore
    try:
        from anndata.io import read_elem  # type: ignore
    except ImportError:
        from anndata.experimental import read_elem  # type: ignore
    with h5py.File(str(path), "r") as h:
        mods = list(h["mod"].keys())
        gene_key = "gene" if "gene" in mods else ("rna" if "rna" in mods else mods[0])
        guide_key = "guide" if "guide" in mods else ("gRNA" if "gRNA" in mods else mods[-1])
        gene = read_elem(h[f"mod/{gene_key}"])
        guide = read_elem(h[f"mod/{guide_key}"])
        obs = read_elem(h["obs"]) if "obs" in h else gene.obs
        uns = read_elem(h["uns"]) if (need_uns and "uns" in h) else {}
    return Modalities(gene, guide, obs, uns, path, "h5py")


def uns_table(mods: Modalities, key: str):
    """A DataFrame stored in uns (mudata writes DataFrames as groups)."""
    pd = _pd()
    v = mods.uns.get(key)
    if v is None:
        raise SystemExit(f"'{key}' not found in mdata.uns; available: {sorted(mods.uns)}")
    if isinstance(v, dict):
        return pd.DataFrame({k: (list(x) if not hasattr(x, "shape") else x) for k, x in v.items()})
    return pd.DataFrame(v)


def guide_var_normalised(guide_var, non_targeting_label: str = "non_targeting"):
    """guide.var with the columns the QC needs, inferred when a newer schema lacks them."""
    np = _np()
    gv = guide_var.copy()
    if "guide_id" not in gv.columns:
        gv["guide_id"] = gv.index.astype(str)
    if "targeting" not in gv.columns:
        if "type" in gv.columns:
            gv["targeting"] = ~gv["type"].astype(str).str.lower().str.replace(r"[-\s]", "_", regex=True).eq(
                non_targeting_label.lower().replace("-", "_"))
        else:
            gv["targeting"] = True
    if "gene_name" not in gv.columns:
        gv["gene_name"] = gv["guide_id"].astype(str).str.split("#").str[0]
        gv.loc[~gv["targeting"].astype(bool), "gene_name"] = np.nan
    if "label" not in gv.columns:
        if "type" in gv.columns:
            gv["label"] = gv["type"].astype(str).str.replace("-", "_", regex=False)
        else:
            gv["label"] = np.where(gv["targeting"].astype(bool), "targeting", non_targeting_label)
    return gv


def is_non_targeting(guide_var, non_targeting_label: str = "non_targeting"):
    if "type" in guide_var.columns:
        norm = guide_var["type"].astype(str).str.lower().str.replace(r"[-\s]", "_", regex=True)
        return norm == non_targeting_label.lower().replace("-", "_").replace(" ", "_")
    if "targeting" in guide_var.columns:
        return ~guide_var["targeting"].astype(bool)
    raise SystemExit("cannot identify non-targeting guides: guide.var has neither 'type' nor 'targeting'")


# ---------------------------------------------------------------------------
# Stage 3a: gene mapping QC
# ---------------------------------------------------------------------------

def gene_metrics(gene, batch_label: str, umi_col: str, genes_col: str, mito_col: str) -> dict:
    m = {"batch": batch_label, "n_cells": int(gene.n_obs)}
    for col, prefix in ((umi_col, "umi"), (genes_col, "genes"), (mito_col, "mito")):
        if col in gene.obs.columns:
            for k, v in column_stats(gene.obs[col].astype(float)).items():
                m[f"{prefix}_{k}"] = v
    return m


def cmd_qc_gene(args) -> int:
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    out = run_dir(args.label or "qc_gene")
    mods = read_mudata(Path(args.mudata), need_uns=False)
    gene = mods.gene
    umi_col, genes_col, mito_col, batch_col = args.umi_col, args.genes_col, args.mito_col, args.batch_col
    # Fall back to what the pipeline actually writes.
    if umi_col not in gene.obs.columns and "total_counts" in gene.obs.columns:
        umi_col = "total_counts"
    if umi_col not in gene.obs.columns:
        gene.obs[umi_col] = _row_sums(gene.X)
    if genes_col not in gene.obs.columns:
        if "log1p_n_genes_by_counts" in gene.obs.columns:
            gene.obs[genes_col] = np.expm1(gene.obs["log1p_n_genes_by_counts"].astype(float)).round().astype(int)
        else:
            X = gene.X
            gene.obs[genes_col] = np.asarray((X > 0).sum(axis=1)).ravel()
    if mito_col not in gene.obs.columns and "mt" in gene.var.columns:
        mt = gene.var["mt"].astype(bool).to_numpy()
        tot = _row_sums(gene.X)
        gene.obs[mito_col] = 100.0 * _row_sums(gene.X[:, mt]) / np.where(tot > 0, tot, 1)
    rows = [gene_metrics(gene, "all", umi_col, genes_col, mito_col)]
    if batch_col in gene.obs.columns:
        for b in sorted(gene.obs[batch_col].astype(str).unique()):
            rows.append(gene_metrics(gene[gene.obs[batch_col].astype(str) == b], b, umi_col, genes_col, mito_col))
    df = pd.DataFrame(rows)
    metrics_path = out / f"{args.prefix}_metrics.tsv"
    df.to_csv(metrics_path, sep="\t", index=False)
    figs = []
    plt = _plt()
    if plt and not args.no_plots:
        fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
        for ax, (col, ttl) in zip(axes, ((umi_col, "UMIs per cell"), (genes_col, "genes per cell"), (mito_col, "mito %"))):
            if col in gene.obs.columns:
                v = gene.obs[col].astype(float)
                ax.hist(np.log10(v + 1) if col != mito_col else v, bins=60, color=BLUE)
                _style(ax, ttl, "log10(x + 1)" if col != mito_col else "%", "cells")
        figs.append(_save(fig, out / f"{args.prefix}_histograms.png"))
        fig, ax = plt.subplots(figsize=(6, 4))
        v = np.sort(gene.obs[umi_col].astype(float).to_numpy())[::-1]
        ax.plot(np.arange(1, len(v) + 1), v, lw=1.2, color=BLUE)
        ax.set_xscale("log")
        ax.set_yscale("log")
        _style(ax, "Knee plot: gene UMIs per cell", "cell rank", "UMIs")
        figs.append(_save(fig, out / f"{args.prefix}_knee.png"))
    r = rows[0]
    print(f"Cells: {r['n_cells']:,}; median UMIs {r.get('umi_median', float('nan')):,.0f}; median genes "
          f"{r.get('genes_median', float('nan')):,.0f}; median mito {r.get('mito_median', float('nan')):.2f}%  "
          f"(batches: {len(rows) - 1})")
    print(f"TSV: {metrics_path}")
    for f in figs:
        print(f"Figure: {f}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Stage 3b: guide mapping QC
# ---------------------------------------------------------------------------

def guide_assignment_counts(guide, layer: str = "guide_assignment"):
    if layer not in guide.layers:
        raise SystemExit(f"layer {layer!r} not in guide.layers ({list(guide.layers.keys())})")
    A = guide.layers[layer]
    is_assigned = A > 0
    guide.obs["n_guides_per_cell"] = _row_sums(is_assigned)
    guide.var["n_cells_per_guide"] = _col_sums(is_assigned)
    return is_assigned


def guide_metrics(guide, batch_label: str, umi_col: str, include_per_guide: bool) -> dict:
    m = {"batch": batch_label, "n_cells": int(guide.n_obs)}
    if umi_col in guide.obs.columns:
        for k, v in column_stats(guide.obs[umi_col].astype(float)).items():
            m[f"guide_umi_{k}"] = v
    g = guide.obs["n_guides_per_cell"].astype(float)
    m.update({"guides_per_cell_mean": float(g.mean()), "guides_per_cell_std": float(g.std()),
              "guides_per_cell_min": float(g.min()), "guides_per_cell_max": float(g.max()),
              "guides_per_cell_median": float(g.median()),
              "n_cells_with_guide": int((g > 0).sum()), "n_cells_exactly_1_guide": int((g == 1).sum()),
              "frac_cells_with_guide": float((g > 0).mean()) if len(g) else 0.0})
    if include_per_guide:
        c = guide.var["n_cells_per_guide"].astype(float)
        m.update({"n_guides_total": int(guide.n_vars), "cells_per_guide_median": float(c.median()),
                  "cells_per_guide_mean": float(c.mean()), "cells_per_guide_std": float(c.std()),
                  "cells_per_guide_min": float(c.min()), "cells_per_guide_max": float(c.max())})
    return m


def per_guide_capture(guide, label_col: str = "label"):
    pd = _pd()
    np = _np()
    X = guide.X
    A = guide.layers.get("guide_assignment")
    n_cells = guide.n_obs
    rows = []
    for i, gid in enumerate(guide.var_names):
        counts = _to_dense_col(X, i)
        det = counts > 0
        n_det = int(det.sum())
        n_asg = int((_to_dense_col(A, i) > 0).sum()) if A is not None else float("nan")
        vals = counts[det]
        rows.append({"guide_id": gid, "n_cells_detected": n_det, "frac_cells_detected": n_det / n_cells,
                     "n_cells_assigned": n_asg, "frac_cells_assigned": (n_asg / n_cells) if A is not None else float("nan"),
                     "total_umi": float(vals.sum()) if n_det else 0.0,
                     "mean_umi": float(vals.mean()) if n_det else float("nan"),
                     "median_umi": float(np.median(vals)) if n_det else float("nan"),
                     "std_umi": float(vals.std()) if n_det else float("nan"),
                     "max_umi": float(vals.max()) if n_det else 0.0})
    df = pd.DataFrame(rows)
    df["total_cells"] = n_cells
    for col in (label_col, "gene_name", "type"):
        if col in guide.var.columns:
            df[col] = df["guide_id"].map(guide.var[col].to_dict())
    return df


def cmd_qc_guide(args) -> int:
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    out = run_dir(args.label or "qc_guide")
    mods = read_mudata(Path(args.mudata), need_uns=False)
    guide = mods.guide
    guide_assignment_counts(guide, args.assignment_layer)
    umi_col = args.umi_col if args.umi_col in guide.obs.columns else None
    if umi_col is None:
        guide.obs[args.umi_col] = _row_sums(guide.X)
        umi_col = args.umi_col
    rows = [guide_metrics(guide, "all", umi_col, True)]
    if args.batch_col in guide.obs.columns:
        for b in sorted(guide.obs[args.batch_col].astype(str).unique()):
            rows.append(guide_metrics(guide[guide.obs[args.batch_col].astype(str) == b], b, umi_col, False))
    df = pd.DataFrame(rows)
    metrics_path = out / f"{args.prefix}_metrics.tsv"
    df.to_csv(metrics_path, sep="\t", index=False)
    cap = per_guide_capture(guide)
    cap_path = out / f"{args.prefix}_per_guide_capture.tsv"
    cap.to_csv(cap_path, sep="\t", index=False)
    figs = []
    plt = _plt()
    if plt and not args.no_plots:
        fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
        axes[0].hist(np.log10(guide.obs[umi_col].astype(float) + 1), bins=60, color=BLUE)
        _style(axes[0], "guide UMIs per cell", "log10(UMIs + 1)", "cells")
        g = guide.obs["n_guides_per_cell"].astype(int)
        axes[1].hist(g, bins=range(0, int(g.max()) + 2), color=BLUE)
        _style(axes[1], "guides assigned per cell", "guides", "cells")
        axes[2].hist(np.log10(guide.var["n_cells_per_guide"].astype(float) + 1), bins=50, color=BLUE)
        _style(axes[2], "cells per guide", "log10(cells + 1)", "guides")
        figs.append(_save(fig, out / f"{args.prefix}_histograms.png"))
    r = rows[0]
    print(f"Cells: {r['n_cells']:,}; with >= 1 guide {r['n_cells_with_guide']:,} ({r['frac_cells_with_guide']:.1%}); "
          f"guides/cell median {r['guides_per_cell_median']:.0f}, mean {r['guides_per_cell_mean']:.2f}; "
          f"guides {r['n_guides_total']:,}, cells/guide median {r['cells_per_guide_median']:.0f}")
    print(f"TSV: {metrics_path}")
    print(f"TSV: {cap_path}")
    for f in figs:
        print(f"Figure: {f}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Stage 3c: intended-target knockdown QC
# ---------------------------------------------------------------------------

def intended_target_table(trans, guide_var, log2fc_col: str, pvalue_col: str):
    """Trans tests where the guide's gene_name equals the target parsed from guide_id."""
    id_map = guide_var.drop_duplicates("intended_target_name").set_index("intended_target_name")["gene_name"]
    res = trans.copy()
    res["gene_name"] = res["gene_id"].map(id_map)
    res["target_name"] = res["guide_id"].astype(str).str.split("#").str[0]
    intended = res[res["gene_name"] == res["target_name"]].drop_duplicates(["guide_id", "gene_id"])
    merged = guide_var.reset_index(drop=True).merge(intended[["guide_id", "gene_id", log2fc_col, pvalue_col]],
                                                    on="guide_id", how="left")
    return merged


def knockdown_metrics(intended, log2fc_col: str, pvalue_col: str, fc_threshold: float, pval_threshold: float) -> dict:
    np = _np()
    valid = intended.dropna(subset=[log2fc_col, pvalue_col])
    n = len(valid)
    thr = math.log2(fc_threshold)
    strong = valid[log2fc_col] <= thr
    sig = valid[pvalue_col] < pval_threshold
    return {"n_guides_total": int(len(intended)), "n_guides_tested": int(n), "fc_threshold": fc_threshold,
            "log2fc_threshold": thr, "pval_threshold": pval_threshold,
            "n_strong_knockdowns": int(strong.sum()), "n_significant": int(sig.sum()),
            "n_strong_and_significant": int((strong & sig).sum()),
            "frac_strong_knockdowns": float(strong.mean()) if n else 0.0,
            "frac_significant": float(sig.mean()) if n else 0.0,
            "frac_strong_and_significant": float((strong & sig).mean()) if n else 0.0,
            "median_log2fc": float(valid[log2fc_col].median()) if n else float("nan"),
            "mean_log2fc": float(valid[log2fc_col].mean()) if n else float("nan")}


def evaluation_table(trans, guide_var, pvalue_col: str, non_targeting_label: str, seed: int = 42):
    """Positives: targeting guide vs its own target. Negatives: non-targeting guides vs the
    intended-target genes, downsampled to the number of positives (CRISPR_Pipeline's
    evaluate_controls.py convention)."""
    pd = _pd()
    id_map = guide_var.drop_duplicates("intended_target_name").set_index("intended_target_name")["gene_name"]
    res = trans.copy()
    res["gene_name"] = res["gene_id"].map(id_map)
    res["target_name"] = res["guide_id"].astype(str).str.split("#").str[0]
    res = res.merge(guide_var.reset_index(drop=True)[["guide_id", "label"]].drop_duplicates(), on="guide_id", how="left")
    targets = set(guide_var["gene_name"].dropna().unique())
    pos = res[(res["gene_name"] == res["target_name"]) & (res["label"] != non_targeting_label) & res[pvalue_col].notna()].copy()
    neg = res[(res["label"] == non_targeting_label) & res["gene_name"].isin(targets) & res[pvalue_col].notna()].copy()
    if len(pos) == 0 or len(neg) == 0:
        return pd.DataFrame()
    if len(neg) > len(pos):
        neg = neg.sample(n=len(pos), random_state=seed)
    pos["direct_target"], neg["direct_target"] = 1, 0
    ev = pd.concat([pos, neg], ignore_index=True)[["guide_id", "gene_id", pvalue_col, "direct_target"]]
    return ev.rename(columns={pvalue_col: "p_value"})


def roc_pr(labels, scores) -> "tuple[float, float, dict]":
    """AUROC and AUPRC (trapezoid, sklearn convention) with tie-aware thresholds."""
    np = _np()
    y = np.asarray(labels, dtype=bool)
    s = np.asarray(scores, dtype=float)
    order = np.argsort(-s, kind="mergesort")
    y, s = y[order], s[order]
    P, N = int(y.sum()), int((~y).sum())
    if P == 0 or N == 0:
        return float("nan"), float("nan"), {}
    distinct = np.r_[np.where(np.diff(s))[0], len(s) - 1]
    tps = np.cumsum(y)[distinct]
    fps = np.cumsum(~y)[distinct]
    tpr = np.r_[0, tps / P]
    fpr = np.r_[0, fps / N]
    auroc = float(np.trapz(tpr, fpr))
    prec = tps / (tps + fps)
    rec = tps / P
    prec = np.r_[1.0, prec]
    rec = np.r_[0.0, rec]
    auprc = float(np.trapz(prec, rec))
    return auroc, auprc, {"fpr": fpr.tolist(), "tpr": tpr.tolist(), "precision": prec.tolist(), "recall": rec.tolist()}


def cmd_qc_target(args) -> int:
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    out = run_dir(args.label or "qc_target")
    mods = read_mudata(Path(args.mudata), need_uns=True)
    gv = guide_var_normalised(mods.guide.var, args.non_targeting_label)
    trans = uns_table(mods, args.results_key)
    for c in ("guide_id", "gene_id", args.log2fc_col, args.pvalue_col):
        if c not in trans.columns:
            raise SystemExit(f"{args.results_key} lacks column {c!r}; has {list(trans.columns)[:12]}")
    intended = intended_target_table(trans, gv, args.log2fc_col, args.pvalue_col)
    metrics = knockdown_metrics(intended, args.log2fc_col, args.pvalue_col, args.fc_threshold, args.pval_threshold)
    ev = evaluation_table(trans, gv, args.pvalue_col, args.non_targeting_label)
    curves: dict = {}
    if len(ev):
        auroc, auprc, curves = roc_pr(ev["direct_target"].to_numpy() == 1, 1.0 - ev["p_value"].to_numpy())
        metrics.update({"auroc": auroc, "auprc": auprc, "n_eval_positives": int((ev["direct_target"] == 1).sum()),
                        "n_eval_negatives": int((ev["direct_target"] == 0).sum())})
    else:
        metrics.update({"auroc": float("nan"), "auprc": float("nan"), "n_eval_positives": 0, "n_eval_negatives": 0})
    metrics_path = out / f"{args.prefix}_metrics.tsv"
    pd.DataFrame([metrics]).to_csv(metrics_path, sep="\t", index=False)
    cols = [c for c in ("guide_id", "gene_name", "label", args.log2fc_col, args.pvalue_col) if c in intended.columns]
    results_path = out / f"{args.prefix}_results.tsv"
    intended[cols].to_csv(results_path, sep="\t", index=False)
    figs = []
    plt = _plt()
    if plt and not args.no_plots:
        valid = intended.dropna(subset=[args.log2fc_col, args.pvalue_col])
        fig, axes = plt.subplots(1, 3 if curves else 1, figsize=(14 if curves else 5, 4), squeeze=False)
        ax = axes[0, 0]
        ax.scatter(valid[args.log2fc_col], -np.log10(valid[args.pvalue_col].clip(lower=1e-300)), s=6, color=BLUE, alpha=0.6, lw=0)
        ax.axvline(metrics["log2fc_threshold"], color=ORANGE, ls="--", lw=0.8)
        ax.axhline(-math.log10(args.pval_threshold), color=ORANGE, ls="--", lw=0.8)
        _style(ax, "Intended-target knockdown", "log2 fold change", "-log10 p")
        if curves:
            axes[0, 1].plot(curves["fpr"], curves["tpr"], color=BLUE, lw=1.6)
            axes[0, 1].plot([0, 1], [0, 1], color=AXIS, ls="--", lw=0.8)
            _style(axes[0, 1], f"ROC (AUROC {metrics['auroc']:.3f})", "false positive rate", "true positive rate")
            axes[0, 2].plot(curves["recall"], curves["precision"], color=BLUE, lw=1.6)
            _style(axes[0, 2], f"PR (AUPRC {metrics['auprc']:.3f})", "recall", "precision")
        figs.append(_save(fig, out / f"{args.prefix}_plots.png"))
    print(f"Intended-target tests: {metrics['n_guides_tested']:,} of {metrics['n_guides_total']:,} guides; strong knockdown "
          f"(FC <= {args.fc_threshold}) {metrics['frac_strong_knockdowns']:.1%}; significant (p < {args.pval_threshold}) "
          f"{metrics['frac_significant']:.1%}; median log2FC {metrics['median_log2fc']:.2f}")
    if len(ev):
        print(f"Targeting vs non-targeting: AUROC {metrics['auroc']:.3f}, AUPRC {metrics['auprc']:.3f} "
              f"({metrics['n_eval_positives']} positives / {metrics['n_eval_negatives']} negatives)")
    print(f"TSV: {metrics_path}")
    print(f"TSV: {results_path}")
    for f in figs:
        print(f"Figure: {f}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Calibration of PerTurbo per-element results against the non-targeting null
# ---------------------------------------------------------------------------

def empirical_pvals_ecdf(null_stats, test_stats, two_sided: bool = True):
    """p = (r + 1) / (B + 1), r = number of null statistics at least as extreme."""
    np = _np()
    null = np.asarray(null_stats, dtype=float)
    null = null[np.isfinite(null)]
    test = np.asarray(test_stats, dtype=float)
    if two_sided:
        null, test = np.abs(null), np.abs(test)
    ns = np.sort(null)
    B = ns.size
    r = B - np.searchsorted(ns, test, side="right")
    return (r + 1.0) / (B + 1.0)


def empirical_pvals_tfit(null_stats, test_stats, two_sided: bool = True, winsor: Optional[float] = 0.01):
    """Fit a t distribution (location 0) to winsorised null z-values; p from its survival function."""
    np = _np()
    from scipy import stats  # type: ignore
    null = np.asarray(null_stats, dtype=float)
    null = null[np.isfinite(null)]
    if winsor is not None:
        lo, hi = np.quantile(null, [winsor, 1.0 - winsor])
        null = np.clip(null, lo, hi)
    df_hat, _, scale_hat = stats.t.fit(null, floc=0.0)
    z = np.asarray(test_stats, dtype=float) / scale_hat
    return 2.0 * stats.t.sf(np.abs(z), df_hat) if two_sided else stats.t.sf(z, df_hat)


def annotate_cis(trans, gene_var, cis_window: int):
    """Same chromosome and element-to-gene midpoint distance within the window."""
    pd = _pd()
    if "intended_target_chr" not in trans.columns:
        return pd.Series(False, index=trans.index)
    coords = None
    for chr_col, s_col, e_col in (("gene_chr", "gene_start", "gene_end"), ("chromosome", "start", "end"),
                                  ("chr", "start", "end"), ("seqname", "start", "end")):
        if chr_col in gene_var.columns and s_col in gene_var.columns:
            coords = gene_var[[chr_col, s_col, e_col]].copy()
            coords.columns = ["gene_chr", "gene_start", "gene_end"]
            break
    if coords is None:
        return pd.Series(False, index=trans.index)
    coords["gene_id"] = gene_var.index.astype(str).values
    m = trans[["gene_id", "intended_target_chr", "intended_target_start", "intended_target_end"]].reset_index(drop=True).merge(
        coords.reset_index(drop=True), on="gene_id", how="left")
    same = m["intended_target_chr"].astype(str) == m["gene_chr"].astype(str)
    dist = ((m["intended_target_start"].astype(float) + m["intended_target_end"].astype(float)) / 2
            - (m["gene_start"].astype(float) + m["gene_end"].astype(float)) / 2).abs()
    is_cis = (same & (dist <= cis_window)).fillna(False)
    is_cis.index = trans.index
    return is_cis


CALIBRATED_COLS = ["element_id", "element_symbol", "element_label", "tested_gene_id", "tested_gene_symbol", "n_cells",
                   "log2fc", "log2fc_se", "is_cis", "is_direct_target", "posterior_pval", "empirical_pval", "empirical_pval_adj"]


def cmd_calibrate(args) -> int:
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    out = run_dir(args.label or f"calibrate_{args.prefix}")
    trans = pd.read_csv(args.trans_results, sep="\t", low_memory=False)
    mods = read_mudata(Path(args.mudata), need_uns=False)
    gv, gene_var = mods.guide.var, mods.gene.var
    if "log2_fc_std" not in trans.columns:
        print("ERROR: trans results lack `log2_fc_std`; calibration needs z = log2_fc / log2_fc_std "
              "(sceptre_v11-style outputs without standard errors cannot be calibrated)", file=sys.stderr)
        return 2
    ntc_mask_var = is_non_targeting(gv, args.non_targeting_label)
    ntc_elements = set(gv.loc[ntc_mask_var, "intended_target_name"].astype(str))
    # an element is null if its guides are non-targeting in the guide metadata, or if its name itself
    # carries the non-targeting label (per-element outputs often name them non-targeting / non-targeting_<k> / NTC...)
    norm_label = args.non_targeting_label.lower().replace("-", "_").replace(" ", "_")
    names_norm = trans["intended_target_name"].astype(str).str.lower().str.replace("-", "_", regex=False).str.replace(" ", "_", regex=False)
    is_ntc = (trans["intended_target_name"].astype(str).isin(ntc_elements) | names_norm.str.startswith(norm_label)
              | names_norm.str.startswith("non_targeting") | names_norm.str.startswith("ntc") | names_norm.eq("safe_targeting"))
    valid = trans["log2_fc_std"].notna() & (trans["log2_fc_std"] > 0)
    trans["z_value"] = np.where(valid, trans["log2_fc"] / trans["log2_fc_std"], np.nan)
    null = trans.loc[is_ntc, "z_value"].dropna().to_numpy()
    if null.size == 0:
        print("ERROR: no non-targeting z-values to build the null from", file=sys.stderr)
        return 2
    fn = empirical_pvals_ecdf if args.method == "ecdf" else empirical_pvals_tfit
    trans["empirical_pval"] = fn(null, trans["z_value"].to_numpy(), two_sided=not args.one_sided)
    trans.loc[trans["z_value"].isna(), "empirical_pval"] = np.nan
    trans["empirical_pval_adj"] = np.nan
    disc = (~is_ntc) & trans["empirical_pval"].notna()
    trans.loc[disc, "empirical_pval_adj"] = bh_adjust(trans.loc[disc, "empirical_pval"].to_numpy())
    # symbols, labels, cells
    sym_map = gene_var["symbol"].to_dict() if "symbol" in gene_var.columns else {}
    el_sym = {}
    if "gene_name" in gv.columns:
        el_sym = gv.drop_duplicates("intended_target_name").set_index("intended_target_name")["gene_name"].to_dict()
        if not any(isinstance(v, str) and v and not v.startswith("ENSG") for v in el_sym.values()):
            el_sym = {}
    if not el_sym:
        el_sym = sym_map
    el_type = gv.drop_duplicates("intended_target_name").set_index("intended_target_name")["type"].to_dict() if "type" in gv.columns else {}
    trans["element_symbol"] = trans["intended_target_name"].map(el_sym)
    trans["element_label"] = trans["intended_target_name"].map(el_type)
    trans["tested_gene_symbol"] = trans["gene_id"].map(sym_map)
    # cells per element from the assignment layer
    if "guide_assignment" in mods.guide.layers:
        A = mods.guide.layers["guide_assignment"] > 0
        per_guide = _col_sums(A)
        cells_by_elem = pd.Series(per_guide, index=gv["intended_target_name"].astype(str).values).groupby(level=0).sum()
        trans["n_cells"] = trans["intended_target_name"].astype(str).map(cells_by_elem)
    trans["is_direct_target"] = trans["gene_id"].astype(str) == trans["intended_target_name"].astype(str)
    trans["is_cis"] = annotate_cis(trans, gene_var, args.cis_window)
    trans.loc[trans["is_direct_target"], "is_cis"] = True
    trans["element_id"] = trans["intended_target_name"].astype(str)
    if {"intended_target_chr", "intended_target_start", "intended_target_end"} <= set(trans.columns):
        st = pd.to_numeric(trans["intended_target_start"], errors="coerce")
        en = pd.to_numeric(trans["intended_target_end"], errors="coerce")
        has = st.notna() & en.notna() & trans["intended_target_chr"].astype(str).str.len().gt(0) & (trans["intended_target_chr"].astype(str) != "nan")
        trans.loc[has, "element_id"] = (trans.loc[has, "intended_target_name"].astype(str) + "|" + trans.loc[has, "intended_target_chr"].astype(str)
                                        + ":" + st[has].astype(int).astype(str) + "-" + en[has].astype(int).astype(str))
    outdf = trans.rename(columns={"gene_id": "tested_gene_id", "log2_fc": "log2fc", "log2_fc_std": "log2fc_se", "p_value": "posterior_pval"})
    outdf = outdf[[c for c in CALIBRATED_COLS if c in outdf.columns]]
    paths = {}
    for name, sel in (("all_results", slice(None)), ("direct_target_results", outdf["is_direct_target"]),
                      ("cis_results", outdf["is_cis"]), ("trans_results", ~outdf["is_cis"])):
        p = out / f"{args.prefix}_calibrated_{name}.tsv"
        (outdf if isinstance(sel, slice) else outdf[sel]).to_csv(p, sep="\t", index=False)
        paths[name] = p
    n10, n05 = int((outdf["empirical_pval_adj"] < 0.10).sum()), int((outdf["empirical_pval_adj"] < 0.05).sum())
    print(f"Calibration ({args.method}): {len(outdf):,} tests; null from {null.size:,} non-targeting z-values; "
          f"direct-target {int(outdf['is_direct_target'].sum()):,}, cis {int(outdf['is_cis'].sum()):,}, trans {int((~outdf['is_cis']).sum()):,}; "
          f"FDR<0.10: {n10:,}; FDR<0.05: {n05:,}")
    for p in paths.values():
        print(f"TSV: {p}")
    print(f"Calibrated prefix: {out / args.prefix}_calibrated_   (pass to `tf-perturb run --calibrated-prefix`)")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Per-element pathway over-representation
# ---------------------------------------------------------------------------

def load_gmt(path: Path) -> dict:
    sets = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = {"name": parts[1] or parts[0], "genes": {g for g in parts[2:] if g}}
    return sets


def fisher_over(deg: set, pw: set, bg: set) -> "tuple[int, float, float]":
    from scipy.stats import fisher_exact  # type: ignore
    a = len(deg & pw)
    b = len(deg - pw)
    c = len((pw & bg) - deg)
    d = len(bg - pw - deg)
    _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    exp = len(deg) * len(pw & bg) / len(bg) if bg else 0.0
    return a, (a / exp if exp > 0 else float("inf")), float(p)


PATHWAY_COLS = ["element_id", "element_symbol", "pathway_id", "pathway_name", "pathway_source", "n_degs_in_pathway",
                "n_degs_total", "pathway_size", "background_size", "fold_enrichment", "pvalue", "pvalue_adj"]


def cmd_pathways(args) -> int:
    pd = _pd()
    log_path = setup_logging()
    out = run_dir(args.label or f"pathways_{args.prefix}")
    df = pd.read_csv(args.calibrated, sep="\t", low_memory=False)
    sets = {k: v for k, v in load_gmt(Path(args.gmt)).items() if args.min_pathway_size <= len(v["genes"]) <= args.max_pathway_size}
    bg = set(df[args.gene_col].dropna().astype(str))
    for v in sets.values():
        v["genes"] &= bg
    sets = {k: v for k, v in sets.items() if v["genes"]}
    sig = df[df[args.pval_col].notna() & (df[args.pval_col] < args.fdr)]
    degs = sig.groupby("element_id")[args.gene_col].apply(lambda s: set(s.dropna().astype(str))).to_dict()
    symbols = df.drop_duplicates("element_id").set_index("element_id")["element_symbol"].to_dict()
    tested = {e: g for e, g in degs.items() if len(g) >= args.min_genes}
    rows = []
    for eid, genes in tested.items():
        er = []
        for pid, info in sets.items():
            n, fold, p = fisher_over(genes, info["genes"], bg)
            if n < 1:
                continue
            er.append({"element_id": eid, "element_symbol": symbols.get(eid, ""), "pathway_id": pid, "pathway_name": info["name"],
                       "n_degs_in_pathway": n, "n_degs_total": len(genes), "pathway_size": len(info["genes"]),
                       "background_size": len(bg), "fold_enrichment": fold, "pvalue": p})
        if er:
            ed = pd.DataFrame(er)
            ed["pvalue_adj"] = bh_adjust(ed["pvalue"].to_numpy())
            rows.append(ed)
    res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=PATHWAY_COLS)
    res["pathway_source"] = Path(args.gmt).name.split(".")[0]
    res = res[PATHWAY_COLS]
    p = out / f"{args.prefix}_calibrated_pathway_results.tsv"
    res.to_csv(p, sep="\t", index=False)
    n_sig = int((res["pvalue_adj"] < 0.05).sum()) if len(res) else 0
    print(f"Pathways: {len(tested)} elements with >= {args.min_genes} DEGs (FDR < {args.fdr}) x {len(sets)} gene sets "
          f"-> {len(res):,} tests, {n_sig:,} at FDR < 0.05 in {res.loc[res['pvalue_adj'] < 0.05, 'element_id'].nunique() if len(res) else 0} elements")
    print(f"TSV: {p}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# MuData filters
# ---------------------------------------------------------------------------

def guide_counts_h5py(path: Path):
    """Assigned guides per cell straight from the HDF5 layer; no full load."""
    np = _np()
    import h5py  # type: ignore
    with h5py.File(str(path), "r") as h:
        mods = list(h["mod"].keys())
        gk = "guide" if "guide" in mods else mods[-1]
        layer = h[f"mod/{gk}/layers/guide_assignment"]
        enc = layer.attrs.get("encoding-type", "")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        if "csr" in enc.lower():
            return np.diff(layer["indptr"][:]).astype(np.int32)
        if "csc" in enc.lower():
            shape = tuple(layer.attrs["shape"])
            return np.bincount(layer["indices"][:], minlength=shape[0]).astype(np.int32)
        data = layer["data"] if isinstance(layer, h5py.Group) else layer
        counts = np.zeros(data.shape[0], dtype=np.int32)
        for s in range(0, data.shape[0], 50_000):
            counts[s:s + 50_000] = (data[s:s + 50_000] != 0).sum(axis=1)
        return counts


def write_filtered(mods: Modalities, keep_cells, out_path: Path, keep_guides=None) -> "list[Path]":
    """Subset both modalities; write .h5mu with mudata, else one .h5ad per modality."""
    np = _np()
    gene = mods.gene[keep_cells].copy()
    guide = mods.guide[keep_cells].copy()
    if keep_guides is not None:
        guide = guide[:, keep_guides].copy()
    if "guide_assignment" in guide.layers:
        guide.obs["num_expressed_guides"] = _row_sums(guide.layers["guide_assignment"] > 0)
    guide.obs["total_guide_umis"] = _row_sums(guide.X)
    try:
        import mudata as md  # type: ignore
        m = md.MuData({"gene": gene, "guide": guide})
        for c in mods.obs.columns:
            if c not in m.obs.columns:
                m.obs[c] = mods.obs.loc[m.obs_names, c].values if set(m.obs_names) <= set(mods.obs.index) else np.nan
        m.update()
        m.write_h5mu(str(out_path))
        return [out_path]
    except ImportError:
        base = out_path.with_suffix("")
        p1, p2 = Path(f"{base}.gene.h5ad"), Path(f"{base}.guide.h5ad")
        gene.write_h5ad(str(p1))
        guide.write_h5ad(str(p2))
        return [p1, p2]


def cmd_filter_cells(args) -> int:
    np = _np()
    log_path = setup_logging()
    counts = guide_counts_h5py(Path(args.mudata))
    keep = counts <= args.max_guides
    mods = read_mudata(Path(args.mudata), need_uns=False)
    written = write_filtered(mods, keep, Path(args.out))
    print(f"Cells: {len(counts):,}; with >= 1 guide {int((counts > 0).sum()):,}; guides/cell median {int(np.median(counts))}, "
          f"max {int(counts.max())}; removed {int((~keep).sum()):,} with > {args.max_guides} guides; kept {int(keep.sum()):,}")
    for p in written:
        print(f"Wrote: {p}")
    print(f"Log: {log_path}")
    return 0


def outlier_guides(non_targeting_table: Path, targeting_table: Path, fdr: float) -> "tuple[set, dict]":
    pd = _pd()
    found, info = set(), {}
    for name, path in (("non_targeting", non_targeting_table), ("targeting", targeting_table)):
        df = pd.read_csv(path, index_col=0)
        if "pval_outlier" not in df.columns:
            raise SystemExit(f"{path}: expected a 'pval_outlier' column")
        adj = bh_adjust(df["pval_outlier"].astype(float).to_numpy())
        rej = [str(i) for i, a in zip(df.index, adj) if a == a and a < fdr]
        found |= set(rej)
        info[name] = (len(rej), len(df))
    return found, info


def cmd_filter_guides(args) -> int:
    np = _np()
    log_path = setup_logging()
    mods = read_mudata(Path(args.mudata), need_uns=False)
    guide = mods.guide
    out_ids, info = outlier_guides(Path(args.non_targeting_outliers), Path(args.targeting_outliers), args.fdr)
    names = guide.var_names.astype(str)
    mask_var = names.isin(out_ids)
    A = guide.layers["guide_assignment"]
    cells_out = _row_sums(A[:, np.where(mask_var)[0]] > 0) > 0 if mask_var.any() else np.zeros(guide.n_obs, dtype=bool)
    written = write_filtered(mods, ~cells_out, Path(args.out), keep_guides=~mask_var)
    print(f"Outlier guides at FDR {args.fdr}: non-targeting {info['non_targeting'][0]}/{info['non_targeting'][1]}, "
          f"targeting {info['targeting'][0]}/{info['targeting'][1]}; {int(mask_var.sum())} present in the object; "
          f"cells removed {int(cells_out.sum()):,} of {guide.n_obs:,}; guides kept {int((~mask_var).sum()):,}")
    for p in written:
        print(f"Wrote: {p}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Energy distance: preprocessing
# ---------------------------------------------------------------------------

def pca_embedding(gene, n_comps: int = 50, seed: int = 0):
    """normalize_total -> log1p -> scale -> PCA, scanpy when present, numpy otherwise."""
    np = _np()
    try:
        import scanpy as sc  # type: ignore
        ad = gene.copy()
        sc.pp.filter_genes(ad, min_counts=1)
        sc.pp.normalize_total(ad)
        sc.pp.log1p(ad)
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=min(n_comps, min(ad.shape) - 1), random_state=seed)
        return ad.obsm["X_pca"]
    except ImportError:
        X = gene.X.toarray() if hasattr(gene.X, "toarray") else np.asarray(gene.X, dtype=float)
        X = X[:, X.sum(axis=0) > 0]
        tot = X.sum(axis=1, keepdims=True)
        X = np.log1p(X / np.where(tot > 0, tot, 1) * np.median(tot))
        X = (X - X.mean(axis=0)) / np.where(X.std(axis=0) > 0, X.std(axis=0), 1)
        X = np.clip(X, -10, 10)
        k = min(n_comps, min(X.shape) - 1)
        rng = np.random.default_rng(seed)
        Q = rng.standard_normal((X.shape[1], k + 10))
        for _ in range(4):
            Q, _ = np.linalg.qr(X.T @ (X @ Q))
        B = X @ Q
        U, S, Vt = np.linalg.svd(B, full_matrices=False)
        return (U[:, :k] * S[:k])


def target_label(row) -> str:
    if str(row.get("type", "")).lower().replace("_", "-") == "non-targeting" or not row.get("targeting", True):
        return "non-targeting"
    if all(k in row and row[k] == row[k] and str(row[k]) not in ("", "nan") for k in ("intended_target_chr", "intended_target_start", "intended_target_end")):
        return f"{row['intended_target_name']}|{row['intended_target_chr']}:{int(float(row['intended_target_start']))}-{int(float(row['intended_target_end']))}"
    return str(row["intended_target_name"])


def cmd_edist_prep(args) -> int:
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    out = Path(args.out_dir) if args.out_dir else run_dir(args.label or "edist_prep")
    out.mkdir(parents=True, exist_ok=True)
    mods = read_mudata(Path(args.mudata), need_uns=False)
    pcs = pca_embedding(mods.gene, args.n_comps)
    pca_df = pd.DataFrame(pcs, index=mods.gene.obs_names.astype(str))
    pca_df.to_pickle(out / "pca_dataframe.pickle")
    guide = mods.guide
    A = guide.layers["guide_assignment"] if "guide_assignment" in guide.layers else guide.X
    A = A.tocsc() if hasattr(A, "tocsc") else A
    grna = {}
    for j, g in enumerate(guide.var_names.astype(str)):
        col = _to_dense_col(A, j)
        grna[g] = list(guide.obs_names.astype(str)[col > 0])
    with open(out / "gRNA_dictionary.pickle", "wb") as fh:
        pickle.dump(grna, fh)
    gv = guide_var_normalised(guide.var)
    keep = ["guide_id"] + [c for c in ("intended_target_name", "type", "spacer", "targeting", "intended_target_chr",
                                      "intended_target_start", "intended_target_end") if c in gv.columns]
    ann = gv[keep].copy()
    ann["target"] = [target_label(r) for _, r in ann.iterrows()]
    ann.to_csv(out / "annotation_table.csv", index=False)
    print(f"PCA: {pca_df.shape[0]:,} cells x {pca_df.shape[1]} PCs; guides {len(grna):,}; targets "
          f"{ann.loc[ann['target'] != 'non-targeting', 'target'].nunique():,} + non-targeting")
    for f in ("pca_dataframe.pickle", "gRNA_dictionary.pickle", "annotation_table.csv"):
        print(f"Wrote: {out / f}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Energy distance: the statistic, the DISCO test, the permutation test
# ---------------------------------------------------------------------------

def _sqdist(A, B):
    from scipy.spatial.distance import cdist  # type: ignore
    return cdist(A, B, metric="sqeuclidean")


def energy_distance(A, B) -> float:
    """E = 2 mean||a-b||^2 - mean||a-a'||^2 - mean||b-b'||^2 (squared Euclidean, all pairs incl. self)."""
    if len(A) == 0 or len(B) == 0:
        return float("nan")
    return float(2 * _sqdist(A, B).mean() - _sqdist(A, A).mean() - _sqdist(B, B).mean())


def permutation_edist(A, B, n_perm: int, rng) -> "tuple[float, Any]":
    """Observed E(A,B) and E over label permutations of the pooled cells, from one distance matrix."""
    np = _np()
    pooled = np.vstack([A, B])
    D = _sqdist(pooled, pooled)
    n1 = len(A)
    n = len(pooled)

    def e_from(idx1, idx2):
        return 2 * D[np.ix_(idx1, idx2)].mean() - D[np.ix_(idx1, idx1)].mean() - D[np.ix_(idx2, idx2)].mean()

    obs = e_from(np.arange(n1), np.arange(n1, n))
    perms = np.empty(n_perm)
    for i in range(n_perm):
        p = rng.permutation(n)
        perms[i] = e_from(p[:n1], p[n1:])
    return float(obs), perms


def disco_stat(groups: "list") -> float:
    """Distance-components F statistic: between-group vs within-group energy (Rizzo & Szekely)."""
    np = _np()
    k = len(groups)
    pooled = np.vstack(groups)
    N = len(pooled)
    D = _sqdist(pooled, pooled) ** 0.5
    sizes = [len(g) for g in groups]
    bounds = np.cumsum([0] + sizes)
    total = D.sum() / (2 * N)
    within = sum(D[bounds[i]:bounds[i + 1], bounds[i]:bounds[i + 1]].sum() / (2 * sizes[i]) for i in range(k))
    between = total - within
    if within <= 0 or k < 2 or N - k <= 0:
        return float("nan")
    return float((between / (k - 1)) / (within / (N - k)))


def cmd_edist_filter(args) -> int:
    """Step 1: outlier guides among siblings (DISCO + hypergeometric ranks) and among non-targeting (K-means)."""
    pd = _pd()
    np = _np()
    from scipy.stats import hypergeom  # type: ignore
    log_path = setup_logging()
    src = Path(args.prep_dir)
    out = Path(args.out_dir) if args.out_dir else src
    out.mkdir(parents=True, exist_ok=True)
    pca = pd.read_pickle(src / "pca_dataframe.pickle")
    grna = pickle.load(open(src / "gRNA_dictionary.pickle", "rb"))
    ann = pd.read_csv(src / "annotation_table.csv")
    rng = np.random.default_rng(args.seed)

    def cells(g):
        ids = [c for c in grna.get(g, []) if c in pca.index]
        if len(ids) > args.max_cells_per_guide:
            ids = list(rng.choice(ids, args.max_cells_per_guide, replace=False))
        return pca.loc[ids].to_numpy()

    # --- targeting guides, per target
    t_rows = []
    for target, sub in ann[ann["target"] != "non-targeting"].groupby("target"):
        guides = [g for g in sub["guide_id"].astype(str) if len(grna.get(g, [])) >= args.min_cells]
        if len(guides) < 2:
            for g in sub["guide_id"].astype(str):
                t_rows.append({"gRNA_name": g, "target": target, "n_cells": len(grna.get(g, [])), "disco_pvalue": float("nan"),
                               "pval_outlier": 1.0, "note": "fewer than 2 guides with enough cells"})
            continue
        groups = [cells(g) for g in guides]
        f_obs = disco_stat(groups)
        pooled = np.vstack(groups)
        sizes = [len(g) for g in groups]
        f_perm = []
        for _ in range(args.disco_permutations):
            p = rng.permutation(len(pooled))
            bounds = np.cumsum([0] + sizes)
            f_perm.append(disco_stat([pooled[p[bounds[i]:bounds[i + 1]]] for i in range(len(sizes))]))
        disco_p = float(np.mean(np.asarray(f_perm) > f_obs)) if f_obs == f_obs else float("nan")
        # pairwise energy distances among sibling guides; hypergeometric right tail on how often a
        # guide's pairs fall in the top half of all intra-target distances
        pairs = list(combinations(range(len(guides)), 2))
        dist = {(i, j): energy_distance(groups[i], groups[j]) for i, j in pairs}
        order = sorted(pairs, key=lambda ij: -dist[ij])
        n_top = max(1, int(round(args.significance_fraction * len(pairs))))
        top = set(order[:n_top])
        for gi, g in enumerate(guides):
            mine = [ij for ij in pairs if gi in ij]
            x = sum(1 for ij in mine if ij in top)
            if disco_p == disco_p and disco_p <= args.disco_alpha and len(guides) > 2:
                pval = float(hypergeom.sf(x - 1, len(pairs), n_top, len(mine)))
            else:
                pval = 1.0
            t_rows.append({"gRNA_name": g, "target": target, "n_cells": len(grna.get(g, [])), "disco_pvalue": disco_p,
                           "n_pairs": len(mine), "n_pairs_in_top": x, "pval_outlier": pval, "note": ""})
    tdf = pd.DataFrame(t_rows).set_index("gRNA_name") if t_rows else pd.DataFrame(columns=["pval_outlier"])
    tdf.to_csv(out / "targeting_outlier_table.csv")
    # --- non-targeting guides: pairwise E, K-means k=2 on the distance matrix, minority cluster = outliers
    nt = [g for g in ann.loc[ann["target"] == "non-targeting", "guide_id"].astype(str) if len(grna.get(g, [])) >= args.min_cells]
    nt_rows = []
    if len(nt) >= 3:
        mats = [cells(g) for g in nt]
        M = np.zeros((len(nt), len(nt)))
        for i, j in combinations(range(len(nt)), 2):
            M[i, j] = M[j, i] = energy_distance(mats[i], mats[j])
        try:
            from sklearn.cluster import KMeans  # type: ignore
            lab = KMeans(n_clusters=2, n_init=10, random_state=args.seed).fit_predict(M)
        except ImportError:
            med = np.median(M.mean(axis=1))
            lab = (M.mean(axis=1) > med).astype(int)
        major = np.bincount(lab).argmax()
        for g, l, row in zip(nt, lab, M):
            nt_rows.append({"gRNA_name": g, "n_cells": len(grna.get(g, [])), "mean_edist_to_others": float(row.sum() / max(1, len(nt) - 1)),
                            "cluster": int(l), "pval_outlier": 0.0 if l != major else 1.0})
    ndf = pd.DataFrame(nt_rows).set_index("gRNA_name") if nt_rows else pd.DataFrame(columns=["pval_outlier"])
    ndf.to_csv(out / "non_targeting_outlier_table.csv")
    n_t_out = int((bh_adjust(tdf["pval_outlier"].to_numpy()) < args.fdr).sum()) if len(tdf) else 0
    n_nt_out = int((ndf["pval_outlier"] < 0.5).sum()) if len(ndf) else 0
    print(f"Targeting guides: {len(tdf)} in {tdf['target'].nunique() if len(tdf) else 0} targets; DISCO significant targets "
          f"{int((tdf.drop_duplicates('target')['disco_pvalue'] <= args.disco_alpha).sum()) if len(tdf) else 0}; "
          f"outlier guides at FDR {args.fdr}: {n_t_out}. Non-targeting guides: {len(ndf)}, K-means outliers {n_nt_out}.")
    print(f"CSV: {out / 'targeting_outlier_table.csv'}")
    print(f"CSV: {out / 'non_targeting_outlier_table.csv'}")
    print(f"Log: {log_path}")
    return 0


def cmd_edist(args) -> int:
    """Step 2: per-target energy distance vs random non-targeting backgrounds, permutation p-values."""
    pd = _pd()
    np = _np()
    log_path = setup_logging()
    src = Path(args.prep_dir)
    out = Path(args.out_dir) if args.out_dir else src
    out.mkdir(parents=True, exist_ok=True)
    pca = pd.read_pickle(src / "pca_dataframe.pickle")
    grna = pickle.load(open(src / "gRNA_dictionary.pickle", "rb"))
    ann = pd.read_csv(src / "annotation_table.csv")
    rng = np.random.default_rng(args.seed)
    excluded = set()
    for name in ("targeting_outlier_table.csv", "non_targeting_outlier_table.csv"):
        p = src / name
        if p.is_file() and not args.keep_outliers:
            t = pd.read_csv(p, index_col=0)
            if "pval_outlier" in t.columns and len(t):
                adj = bh_adjust(t["pval_outlier"].astype(float).to_numpy()) if "targeting" == name.split("_")[0] else t["pval_outlier"].to_numpy()
                excluded |= {str(i) for i, a in zip(t.index, adj) if a == a and a < (args.fdr if name.startswith("targeting") else 0.5)}
    nt_guides = [g for g in ann.loc[ann["target"] == "non-targeting", "guide_id"].astype(str) if g not in excluded]
    nt_cells = sorted({c for g in nt_guides for c in grna.get(g, []) if c in pca.index})
    if len(nt_cells) < 10:
        print("ERROR: fewer than 10 non-targeting cells", file=sys.stderr)
        return 2
    n_pick = min(args.non_target_pick, len(nt_cells))
    backgrounds = [rng.choice(nt_cells, n_pick, replace=False) for _ in range(args.num_bg)]
    d_tmp = 0.1 / (args.permutations * args.num_bg)
    rows = {}
    targets = ann[ann["target"] != "non-targeting"].groupby("target")
    for target, sub in targets:
        guides = [g for g in sub["guide_id"].astype(str) if g not in excluded]
        cells = sorted({c for g in guides for c in grna.get(g, []) if c in pca.index})
        ttype = str(sub["type"].iloc[0]) if "type" in sub.columns else "targeting"
        if len(cells) < args.min_cells:
            continue
        if len(cells) > args.target_cell_max:
            cells = list(rng.choice(cells, args.target_cell_max, replace=False))
        A = pca.loc[cells].to_numpy()
        row = {"cell_count": len(cells), "type": ttype}
        ds, ps = [], []
        for b, bg in enumerate(backgrounds):
            bg_cells = [c for c in bg if c not in set(cells)]
            B = pca.loc[bg_cells].to_numpy()
            obs, perms = permutation_edist(A, B, args.permutations, rng)
            p = float((perms >= obs).sum() / args.permutations)
            row[f"distance_{b}"], row[f"pval_{b}"] = obs, p
            ds.append(obs)
            ps.append(p)
        row["distance_mean"] = float(np.mean(ds))
        row["pval_mean"] = float(np.mean(ps))
        row["pval_mean_log"] = float(-np.log10(row["pval_mean"] + d_tmp))
        row["distance_mean_log"] = float(np.log10(row["distance_mean"])) if row["distance_mean"] > 0 else float("nan")
        rows[target] = row
    df = pd.DataFrame.from_dict(rows, orient="index")
    df.index.name = "target"
    p = out / "pval_edist_full.csv"
    df.to_csv(p)
    n_sig = int((df["pval_mean"] < 0.05).sum()) if len(df) else 0
    nc = df[df["type"].astype(str).str.contains("negative", case=False)] if len(df) else df
    above = int(((~df.index.isin(nc.index)) & (df["distance_mean"] > nc["distance_mean"].max())).sum()) if len(nc) else None
    print(f"Energy distance: {len(df)} targets vs {args.num_bg} backgrounds of {n_pick} non-targeting cells "
          f"({len(nt_cells):,} available; {len(excluded)} outlier guides excluded), {args.permutations} permutations each; "
          f"pval_mean < 0.05: {n_sig}" + (f"; targets above the negative-control max distance: {above}" if above is not None else ""))
    print(f"CSV: {p}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Energy distance: cross-dataset summary and validation
# ---------------------------------------------------------------------------

EDIST_REQUIRED = {"cell_count", "type", "distance_mean", "pval_mean", "pval_mean_log", "distance_mean_log"}
EDIST_TYPES = {"targeting", "target", "positive control", "negative control", "non-targeting"}


def edist_summary_row(dataset: str, df) -> dict:
    row: dict = {"dataset_id": dataset, "n_targets_total": int(len(df))}
    if "type" in df.columns:
        for t in ("targeting", "negative control", "positive control"):
            row[f"n_targets_{t.replace(' ', '_')}"] = int((df["type"] == t).sum())
    if "distance_mean" in df.columns:
        row["distance_mean_median_all"] = float(df["distance_mean"].median())
        if "type" in df.columns:
            for t in ("targeting", "negative control", "positive control"):
                s = df.loc[df["type"] == t, "distance_mean"]
                if len(s):
                    row[f"distance_mean_median_{t.replace(' ', '_')}"] = float(s.median())
            nc = df.loc[df["type"] == "negative control", "distance_mean"]
            if len(nc):
                row["nc_distance_mean_max"] = float(nc.max())
                row["n_targeting_above_nc_max"] = int(((df["type"] == "targeting") & (df["distance_mean"] > nc.max())).sum())
    if "pval_mean" in df.columns:
        row["n_pval_eq_0"] = int((df["pval_mean"] == 0).sum())
        row["n_pval_lt_0p05"] = int((df["pval_mean"] < 0.05).sum())
        if "type" in df.columns:
            row["n_nc_pval_eq_0"] = int(((df["type"] == "negative control") & (df["pval_mean"] == 0)).sum())
    if "cell_count" in df.columns:
        row["cell_count_median"] = float(df["cell_count"].median())
        row["cell_count_min"] = int(df["cell_count"].min())
    return row


def cmd_edist_summary(args) -> int:
    pd = _pd()
    log_path = setup_logging()
    out = run_dir(args.label or "edist_summary")
    rows = []
    for spec in args.inputs:
        name, _, p = spec.partition("=")
        if not p:
            name, p = Path(spec).parent.name, spec
        path = Path(p)
        if path.is_dir():
            path = path / "pval_edist_full.csv"
        rows.append(edist_summary_row(name, pd.read_csv(path, index_col=0)))
    df = pd.DataFrame(rows)
    p = out / "cross_dataset_edistance_summary.tsv"
    df.to_csv(p, sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"TSV: {p}")
    print(f"Log: {log_path}")
    return 0


def cmd_edist_validate(args) -> int:
    pd = _pd()
    src = Path(args.dir)
    problems, notes = [], []
    for name in ("pval_edist_full.csv", "targeting_outlier_table.csv", "non_targeting_outlier_table.csv"):
        p = src / name
        if not p.is_file() or p.stat().st_size == 0:
            problems.append(f"{name}: missing or empty")
            continue
        df = pd.read_csv(p, index_col=0)
        if name == "pval_edist_full.csv":
            miss = EDIST_REQUIRED - set(df.columns)
            if miss:
                problems.append(f"{name}: missing columns {sorted(miss)}")
            if not any(c.startswith("distance_") and c[9:].isdigit() for c in df.columns):
                problems.append(f"{name}: no per-background distance_<i> columns")
            if "cell_count" in df.columns and (df["cell_count"] <= 0).any():
                problems.append(f"{name}: cell_count <= 0 rows")
            if "pval_mean" in df.columns and ((df["pval_mean"] < 0) | (df["pval_mean"] > 1)).any():
                problems.append(f"{name}: pval_mean outside [0, 1]")
            if "type" in df.columns:
                unknown = set(df["type"].dropna().astype(str)) - EDIST_TYPES
                if unknown:
                    notes.append(f"{name}: unexpected type values {sorted(unknown)}")
            notes.append(f"{name}: {len(df)} targets, {len(df.columns)} columns")
        else:
            if "pval_outlier" not in df.columns:
                problems.append(f"{name}: no pval_outlier column")
            notes.append(f"{name}: {len(df)} guides")
    for n in notes:
        print(f"  ok    {n}")
    for p in problems:
        print(f"  FAIL  {p}")
    print(f"Validation: {len(problems)} problem(s)")
    return 1 if problems else 0


# ---------------------------------------------------------------------------
# cNMF export and validation
# ---------------------------------------------------------------------------

def cmd_cnmf_export(args) -> int:
    pd = _pd()
    np = _np()
    import anndata as ad  # type: ignore
    import scipy.sparse as sp  # type: ignore
    log_path = setup_logging()
    mods = read_mudata(Path(args.mudata), need_uns=False)
    g, gd = mods.gene, mods.guide
    var = g.var.copy()
    var.index = var[args.symbol_col].astype(str).values
    var.index.name = args.symbol_col
    keep = ~var.index.duplicated()
    a = ad.AnnData(X=g.X[:, np.where(keep)[0]], obs=g.obs.copy(), var=var[keep])
    A = gd.layers["guide_assignment"]
    a.obsm["guide_assignment"] = A if sp.issparse(A) else sp.csr_matrix(A)
    a.uns["guide_names"] = np.array(gd.var[args.guide_id_col].astype(str).tolist()) if args.guide_id_col in gd.var.columns else np.array(gd.var_names.astype(str))
    a.uns["guide_targets"] = np.array(gd.var[args.target_col].astype(str).tolist())
    dropped = 0
    if not args.no_hvg_filter:
        try:
            import scanpy as sc  # type: ignore
            tmp = a.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.pp.highly_variable_genes(tmp, n_top_genes=min(args.num_highvar_genes, tmp.n_vars - 1), flavor="seurat")
            hvg = tmp.var["highly_variable"].to_numpy()
            keep_cells = _row_sums(a.X[:, np.where(hvg)[0]]) > 0
            dropped = int((~keep_cells).sum())
            a = a[keep_cells].copy()
        except ImportError:
            print("note: scanpy not installed; HVG zero-count cell filter skipped")
    Path(args.out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"PerturbNMF h5ad: {a.shape[0]:,} cells x {a.shape[1]:,} genes (symbols; {int((~keep).sum())} duplicate symbols dropped); "
          f"{len(a.uns['guide_names'])} guides, {pd.Series(a.uns['guide_targets']).nunique()} targets; "
          f"{dropped} cells dropped for zero counts in the {args.num_highvar_genes} HVGs"
          + ("" if "batch" in a.obs.columns else "; WARNING no 'batch' column in obs"))
    print(f"Wrote: {args.out_h5ad}")
    print(f"Log: {log_path}")
    return 0


def cmd_cnmf_validate(args) -> int:
    pd = _pd()
    src = Path(args.dir)
    k = args.selected_k
    problems, notes = [], []
    if not src.is_dir():
        print(f"FAIL  {src} is not a directory")
        return 1
    files = {p.name: p for p in src.rglob("*") if p.is_file()}
    spectra = [p for n, p in files.items() if re.search(rf"gene_spectra_score.*k_?{k}\b", n) or re.search(rf"k_?{k}.*gene_spectra_score", n)]
    usages = [p for n, p in files.items() if re.search(rf"usages.*k_?{k}\b", n) or re.search(rf"k_?{k}.*usages", n)]
    if not spectra:
        problems.append(f"no gene_spectra_score table for k={k}")
    if not usages:
        problems.append(f"no usages table for k={k}")
    for label, group, axis in (("gene_spectra_score", spectra, "rows"), ("usages", usages, "columns")):
        for p in group[:1]:
            try:
                df = pd.read_csv(p, sep="\t" if p.suffix in (".tsv", ".txt") else ",", index_col=0)
            except Exception as exc:
                problems.append(f"{p.name}: unreadable ({exc})")
                continue
            n = df.shape[0] if axis == "rows" else df.shape[1]
            if n != k:
                problems.append(f"{p.name}: {axis} = {n}, expected k = {k}")
            else:
                notes.append(f"{p.name}: {df.shape[0]} x {df.shape[1]}")
            vals = df.select_dtypes("number").to_numpy()
            if label == "usages" and (vals < 0).any():
                problems.append(f"{p.name}: negative usages")
            if not pd.notna(vals).all():
                problems.append(f"{p.name}: NaN values")
    sweep = sorted({int(m.group(1)) for n in files for m in [re.search(r"k_?(\d+)", n)] if m})
    notes.append(f"k values present: {sweep}" if sweep else "no k-sweep files found")
    for n in notes:
        print(f"  ok    {n}")
    for p in problems:
        print(f"  FAIL  {p}")
    print(f"Validation: {len(problems)} problem(s)")
    return 1 if problems else 0


# ---------------------------------------------------------------------------
# Cross-dataset pipeline summary from QC metric tables
# ---------------------------------------------------------------------------

def pipeline_summary_row(dataset: str, tables: dict) -> dict:
    row: dict = {"dataset_id": dataset}
    gm = tables.get("gene")
    if gm is not None and "batch" in gm.columns and (gm["batch"] == "all").any():
        r = gm[gm["batch"] == "all"].iloc[0]
        row.update({"n_cells": int(r.get("n_cells", 0)), "gene_umi_median": float(r.get("umi_median", 0)),
                    "mito_pct_median": float(r.get("mito_median", 0))})
    gu = tables.get("guide")
    if gu is not None and "batch" in gu.columns and (gu["batch"] == "all").any():
        r = gu[gu["batch"] == "all"].iloc[0]
        row.update({"guide_umi_median": float(r.get("guide_umi_median", 0)), "guides_per_cell_mean": float(r.get("guides_per_cell_mean", 0)),
                    "frac_cells_with_guide": float(r.get("frac_cells_with_guide", 0)), "n_guides_total": int(r.get("n_guides_total", 0))})
    it = tables.get("intended_target")
    if it is not None and len(it):
        r = it.iloc[0]
        row.update({"intended_n_guides_tested": int(r.get("n_guides_tested", 0)), "intended_n_significant": int(r.get("n_significant", 0)),
                    "intended_frac_significant": float(r.get("frac_significant", 0)), "intended_median_log2fc": float(r.get("median_log2fc", 0)),
                    "intended_auroc": float(r.get("auroc", float("nan"))), "intended_auprc": float(r.get("auprc", float("nan")))})
    return row


def cmd_pipeline_summary(args) -> int:
    pd = _pd()
    log_path = setup_logging()
    out = run_dir(args.label or "pipeline_summary")
    rows = []
    for spec in args.datasets:
        name, _, d = spec.partition("=")
        if not d:
            name, d = Path(spec).name, spec
        root = Path(d)
        tables = {}
        for key, pat in (("gene", "*gene*metrics.tsv"), ("guide", "*guide*metrics.tsv"), ("intended_target", "*intended_target*metrics.tsv")):
            hits = sorted(root.rglob(pat))
            hits = [h for h in hits if "per_guide" not in h.name]
            if hits:
                tables[key] = pd.read_csv(hits[-1], sep="\t")
        rows.append(pipeline_summary_row(name, tables))
    df = pd.DataFrame(rows)
    p = out / "cross_dataset_pipeline_summary.tsv"
    df.to_csv(p, sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"TSV: {p}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# argparse wiring, attached to the `tf-perturb` CLI by tf_perturb_seq_skill
# ---------------------------------------------------------------------------

def add_stage_parsers(sub) -> None:
    def common(p, default_label: str):
        p.add_argument("--label", default="", help=f"Run-directory label (default {default_label}).")
        p.add_argument("--no-plots", action="store_true")

    p = sub.add_parser("qc-gene", help="Stage 3: gene-expression mapping QC on inference_mudata.h5mu.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--batch-col", default="batch")
    p.add_argument("--umi-col", default="total_gene_umis")
    p.add_argument("--genes-col", default="num_expressed_genes")
    p.add_argument("--mito-col", default="percent_mito")
    p.add_argument("--prefix", default="gene")
    common(p, "qc_gene")
    p.set_defaults(func=cmd_qc_gene)

    p = sub.add_parser("qc-guide", help="Stage 3: guide mapping QC (UMIs, guides/cell, cells/guide, per-guide capture).")
    p.add_argument("--mudata", required=True)
    p.add_argument("--batch-col", default="batch")
    p.add_argument("--umi-col", default="total_guide_umis")
    p.add_argument("--assignment-layer", default="guide_assignment")
    p.add_argument("--prefix", default="guide")
    common(p, "qc_guide")
    p.set_defaults(func=cmd_qc_guide)

    p = sub.add_parser("qc-target", help="Stage 3: intended-target knockdown efficiency and targeting-vs-NT AUROC/AUPRC.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--results-key", default="trans_per_guide_results")
    p.add_argument("--log2fc-col", default="log2_fc")
    p.add_argument("--pvalue-col", default="p_value")
    p.add_argument("--fc-threshold", type=float, default=0.4, help="Fold change for a strong knockdown (0.4 = 60%% down).")
    p.add_argument("--pval-threshold", type=float, default=0.05)
    p.add_argument("--non-targeting-label", default="non_targeting")
    p.add_argument("--prefix", default="intended_target")
    common(p, "qc_target")
    p.set_defaults(func=cmd_qc_target)

    p = sub.add_parser("calibrate", help="Empirical p-values for PerTurbo per-element results against the non-targeting null.")
    p.add_argument("--trans-results", required=True, help="perturbo_trans_per_element_output.tsv(.gz)")
    p.add_argument("--mudata", required=True)
    p.add_argument("--prefix", required=True, help="Output prefix, e.g. <dataset>_<run>")
    p.add_argument("--method", choices=("ecdf", "t-fit"), default="t-fit")
    p.add_argument("--cis-window", type=int, default=100_000)
    p.add_argument("--non-targeting-label", default="non_targeting")
    p.add_argument("--one-sided", action="store_true")
    p.add_argument("--label", default="")
    p.set_defaults(func=cmd_calibrate)

    p = sub.add_parser("pathways", help="Per-element Fisher over-representation of significant DEGs in GMT gene sets.")
    p.add_argument("--calibrated", required=True, help="<prefix>_calibrated_trans_results.tsv (or all_results)")
    p.add_argument("--gmt", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--fdr", type=float, default=0.10)
    p.add_argument("--min-genes", type=int, default=3)
    p.add_argument("--min-pathway-size", type=int, default=10)
    p.add_argument("--max-pathway-size", type=int, default=500)
    p.add_argument("--gene-col", default="tested_gene_symbol")
    p.add_argument("--pval-col", default="empirical_pval_adj")
    p.add_argument("--label", default="")
    p.set_defaults(func=cmd_pathways)

    p = sub.add_parser("filter-cells", help="Drop cells with more than N assigned guides.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--out", required=True, help="Output .h5mu (one .h5ad per modality without the mudata package).")
    p.add_argument("--max-guides", type=int, default=15)
    p.set_defaults(func=cmd_filter_cells)

    p = sub.add_parser("filter-guides", help="Drop outlier guides (BH on outlier tables) and the cells that received them.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--targeting-outliers", required=True, help="targeting_outlier_table.csv")
    p.add_argument("--non-targeting-outliers", required=True, help="non_targeting_outlier_table.csv")
    p.add_argument("--fdr", type=float, default=0.05)
    p.set_defaults(func=cmd_filter_guides)

    p = sub.add_parser("edist-prep", help="Energy distance step 0: PCA, guide->cells map, annotation table.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--out-dir", help="Default: a new Docs/TFPerturbSeq run directory.")
    p.add_argument("--n-comps", type=int, default=50)
    p.add_argument("--label", default="")
    p.set_defaults(func=cmd_edist_prep)

    p = sub.add_parser("edist-filter", help="Energy distance step 1: outlier guides (DISCO + hypergeometric; K-means for non-targeting).")
    p.add_argument("--prep-dir", required=True)
    p.add_argument("--out-dir")
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--max-cells-per-guide", type=int, default=1000)
    p.add_argument("--disco-permutations", type=int, default=200)
    p.add_argument("--disco-alpha", type=float, default=0.05)
    p.add_argument("--significance-fraction", type=float, default=0.5)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_edist_filter)

    p = sub.add_parser("edist", help="Energy distance step 2: per-target E vs random non-targeting backgrounds, permutation p-values.")
    p.add_argument("--prep-dir", required=True)
    p.add_argument("--out-dir")
    p.add_argument("--num-bg", type=int, default=20)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--non-target-pick", type=int, default=2000)
    p.add_argument("--target-cell-max", type=int, default=2000)
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--keep-outliers", action="store_true", help="Ignore outlier tables from edist-filter.")
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_edist)

    p = sub.add_parser("edist-summary", help="Cross-dataset roll-up of pval_edist_full.csv files.")
    p.add_argument("--inputs", action="append", required=True, metavar="NAME=DIR_OR_CSV")
    p.add_argument("--label", default="")
    p.set_defaults(func=cmd_edist_summary)

    p = sub.add_parser("edist-validate", help="Schema and value checks on an energy-distance output folder.")
    p.add_argument("--dir", required=True)
    p.set_defaults(func=cmd_edist_validate)

    p = sub.add_parser("cnmf-export", help="inference_mudata.h5mu -> PerturbNMF single-modality .h5ad.")
    p.add_argument("--mudata", required=True)
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--guide-id-col", default="guide_id")
    p.add_argument("--target-col", default="gene_name")
    p.add_argument("--symbol-col", default="symbol")
    p.add_argument("--num-highvar-genes", type=int, default=2000)
    p.add_argument("--no-hvg-filter", action="store_true")
    p.set_defaults(func=cmd_cnmf_export)

    p = sub.add_parser("cnmf-validate", help="File / shape / value checks on a torch-cNMF run folder.")
    p.add_argument("--dir", required=True)
    p.add_argument("--selected-k", type=int, required=True)
    p.set_defaults(func=cmd_cnmf_validate)

    p = sub.add_parser("pipeline-summary", help="One row per dataset from QC metric tables.")
    p.add_argument("--datasets", action="append", required=True, metavar="NAME=DIR")
    p.add_argument("--label", default="")
    p.set_defaults(func=cmd_pipeline_summary)


STAGE_NAMES = ("qc-gene", "qc-guide", "qc-target", "calibrate", "pathways", "filter-cells", "filter-guides",
               "edist-prep", "edist-filter", "edist", "edist-summary", "edist-validate", "cnmf-export",
               "cnmf-validate", "pipeline-summary")


# ---------------------------------------------------------------------------
# Self-test: a synthetic inference_mudata.h5mu with planted effects, every stage run
# ---------------------------------------------------------------------------

def make_synthetic_h5mu(dirpath: Path, seed: int = 11) -> dict:
    """600 cells x 200 genes; 6 targets x 3 guides + 6 non-targeting + 1 shifted non-targeting outlier.

    Planted: TF1 cells have genes 0-19 knocked down (gene 0 is TF1's own target);
    the outlier NT guide's cells are shifted in genes 100-119; 20 cells carry 3
    guides; 4 cells carry 20 guides (for filter-cells).
    """
    np = _np()
    pd = _pd()
    import h5py  # type: ignore
    import scipy.sparse as sp  # type: ignore
    try:
        from anndata.io import write_elem  # type: ignore
    except ImportError:
        from anndata.experimental import write_elem  # type: ignore
    rng = np.random.default_rng(seed)
    n_cells, n_genes = 600, 200
    genes = [f"ENSG{100000 + i:08d}" for i in range(n_genes)]
    symbols = [f"TF{i + 1}" if i < 6 else f"G{i}" for i in range(n_genes)]
    targets = [f"TF{i + 1}" for i in range(6)]
    target_ids = genes[:6]
    guide_ids = [f"{t}#g{k}" for t in targets for k in range(3)] + [f"NT#{k}" for k in range(6)] + ["NT#outlier"]
    n_guides = len(guide_ids)
    # one guide per cell, round-robin; extra multi-guide cells
    assign = np.zeros((n_cells, n_guides), dtype=np.int8)
    for c in range(n_cells):
        assign[c, c % n_guides] = 1
    for c in range(0, 20):
        assign[c, (c + 5) % n_guides] = 1
        assign[c, (c + 9) % n_guides] = 1
    for c in range(20, 24):
        assign[c, :20] = 1
    base = rng.poisson(lam=np.tile(rng.gamma(2.0, 1.5, size=n_genes), (n_cells, 1)))
    X = base.astype(np.float64)
    tf1_cells = np.where(assign[:, 0:3].sum(axis=1) > 0)[0]
    X[np.ix_(tf1_cells, np.arange(0, 20))] = rng.poisson(0.15 * np.maximum(base[np.ix_(tf1_cells, np.arange(0, 20))], 1))
    out_cells = np.where(assign[:, n_guides - 1] > 0)[0]
    X[np.ix_(out_cells, np.arange(100, 120))] += rng.poisson(12, size=(len(out_cells), 20))
    guide_umi = assign * rng.poisson(30, size=assign.shape) + rng.poisson(0.05, size=assign.shape)
    batch = np.array([f"IGVFDS{c % 2}" for c in range(n_cells)])
    cells = [f"cell{c}" for c in range(n_cells)]
    gene_var = pd.DataFrame({"symbol": symbols, "mt": [False] * n_genes, "gene_chr": ["chr1"] * n_genes,
                             "gene_start": [1_000_000 + i * 100_000 for i in range(n_genes)],
                             "gene_end": [1_000_000 + i * 100_000 + 20_000 for i in range(n_genes)]}, index=genes)
    gene_obs = pd.DataFrame({"batch": batch, "total_gene_umis": X.sum(axis=1), "num_expressed_genes": (X > 0).sum(axis=1),
                             "percent_mito": rng.uniform(1, 8, size=n_cells)}, index=cells)
    gv = pd.DataFrame({
        "guide_id": guide_ids,
        "targeting": [not g.startswith("NT#") for g in guide_ids],
        "type": ["targeting" if not g.startswith("NT#") else "non-targeting" for g in guide_ids],
        "intended_target_name": [target_ids[targets.index(g.split("#")[0])] if not g.startswith("NT#") else "non-targeting" for g in guide_ids],
        "gene_name": [g.split("#")[0] if not g.startswith("NT#") else "non-targeting" for g in guide_ids],
        "intended_target_chr": ["chr1" if not g.startswith("NT#") else "" for g in guide_ids],
        "intended_target_start": [(1_000_000 + targets.index(g.split("#")[0]) * 100_000 - 500) if not g.startswith("NT#") else np.nan for g in guide_ids],
        "intended_target_end": [(1_000_000 + targets.index(g.split("#")[0]) * 100_000 + 500) if not g.startswith("NT#") else np.nan for g in guide_ids],
        "spacer": ["ACGT" * 5] * n_guides,
    }, index=guide_ids)
    guide_obs = pd.DataFrame({"batch": batch, "total_guide_umis": guide_umi.sum(axis=1),
                              "num_expressed_guides": assign.sum(axis=1)}, index=cells)
    # trans_per_guide_results: every guide x the 6 target genes
    rows = []
    for g in guide_ids:
        for t, tid in zip(targets, target_ids):
            own = (not g.startswith("NT#")) and g.split("#")[0] == t
            lfc = -2.0 + rng.normal(0, 0.1) if own else rng.normal(0, 0.1)
            p = 10 ** -rng.uniform(5, 9) if own else rng.uniform(0.02, 1.0)
            rows.append({"guide_id": g, "gene_id": tid, "log2_fc": lfc, "p_value": p})
    trans_guide = pd.DataFrame(rows)
    path = dirpath / "inference_mudata.h5mu"
    with h5py.File(str(path), "w") as h:
        write_elem(h, "obs", pd.DataFrame({"batch": batch}, index=cells))
        write_elem(h, "mod/gene/X", sp.csr_matrix(X))
        write_elem(h, "mod/gene/obs", gene_obs)
        write_elem(h, "mod/gene/var", gene_var)
        write_elem(h, "mod/guide/X", sp.csr_matrix(guide_umi.astype(np.float64)))
        write_elem(h, "mod/guide/obs", guide_obs)
        write_elem(h, "mod/guide/var", gv)
        write_elem(h, "mod/guide/layers/guide_assignment", sp.csr_matrix(assign.astype(np.float64)))
        write_elem(h, "uns/trans_per_guide_results", trans_guide)
        h["mod/gene"].attrs["encoding-type"] = "anndata"
        h["mod/guide"].attrs["encoding-type"] = "anndata"
        h["mod/gene"].attrs["encoding-version"] = "0.1.0"
        h["mod/guide"].attrs["encoding-version"] = "0.1.0"
    # PerTurbo per-element trans output: element x all genes
    elems = [(tid, "chr1", 1_000_000 + i * 100_000 - 500, 1_000_000 + i * 100_000 + 500) for i, tid in enumerate(target_ids)]
    elems += [(f"non-targeting_{k}", "", np.nan, np.nan) for k in range(12)]
    prows = []
    for tid, ch, st, en in elems:
        for gi, gid in enumerate(genes):
            planted = tid == target_ids[0] and gi < 20
            lfc = (-1.5 + rng.normal(0, 0.1)) if planted else rng.normal(0, 0.2)
            se = 0.2
            prows.append({"intended_target_name": tid, "intended_target_chr": ch, "intended_target_start": st,
                          "intended_target_end": en, "gene_id": gid, "log2_fc": lfc, "log2_fc_std": se,
                          "p_value": float(2 * (1 - 0.5 * (1 + math.erf(abs(lfc / se) / math.sqrt(2)))))})
    ppath = dirpath / "perturbo_trans_per_element_output.tsv"
    pd.DataFrame(prows).to_csv(ppath, sep="\t", index=False)
    gmt = dirpath / "sets.gmt"
    with open(gmt, "w") as fh:
        fh.write("PLANTED\tgenes 0-19\t" + "\t".join(symbols[:20]) + "\n")
        for k in range(8):
            fh.write(f"RANDOM{k}\trandom set\t" + "\t".join(rng.choice(symbols[20:], 25, replace=False)) + "\n")
    return {"h5mu": path, "perturbo": ppath, "gmt": gmt, "n_cells": n_cells, "n_genes": n_genes, "n_guides": n_guides,
            "tf1_cells": int(len(tf1_cells)), "outlier_cells": int(len(out_cells))}


def stages_selftest(no_plots: bool = True) -> int:
    import tempfile
    pd = _pd()
    np = _np()
    checks: "list[tuple[bool, str]]" = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def ns(**kw):
        return argparse.Namespace(**kw)

    def latest(label):
        return sorted(OUT_ROOT.glob(f"*_{label}"))[-1]

    made: "list[Path]" = []
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        syn = make_synthetic_h5mu(d)
        h5 = str(syn["h5mu"])
        print("\nstage 3 QC")
        rc = cmd_qc_gene(ns(mudata=h5, batch_col="batch", umi_col="total_gene_umis", genes_col="num_expressed_genes",
                            mito_col="percent_mito", prefix="gene", label="st_qc_gene", no_plots=no_plots))
        m = pd.read_csv(latest("st_qc_gene") / "gene_metrics.tsv", sep="\t"); made.append(latest("st_qc_gene"))
        check(rc == 0 and int(m.loc[m.batch == "all", "n_cells"].iloc[0]) == syn["n_cells"] and len(m) == 3,
              "qc-gene: all + 2 batches, cell count exact")
        rc = cmd_qc_guide(ns(mudata=h5, batch_col="batch", umi_col="total_guide_umis", assignment_layer="guide_assignment",
                             prefix="guide", label="st_qc_guide", no_plots=no_plots))
        g = pd.read_csv(latest("st_qc_guide") / "guide_metrics.tsv", sep="\t"); made.append(latest("st_qc_guide"))
        r = g[g.batch == "all"].iloc[0]
        cap = pd.read_csv(latest("st_qc_guide") / "guide_per_guide_capture.tsv", sep="\t")
        check(rc == 0 and int(r.n_guides_total) == syn["n_guides"] and float(r.frac_cells_with_guide) == 1.0
              and int(r.guides_per_cell_max) == 21 and len(cap) == syn["n_guides"],
              f"qc-guide: {int(r.n_guides_total)} guides, every cell assigned, max 21 guides/cell, per-guide capture table")
        rc = cmd_qc_target(ns(mudata=h5, results_key="trans_per_guide_results", log2fc_col="log2_fc", pvalue_col="p_value",
                              fc_threshold=0.4, pval_threshold=0.05, non_targeting_label="non_targeting", prefix="intended_target",
                              label="st_qc_target", no_plots=no_plots))
        t = pd.read_csv(latest("st_qc_target") / "intended_target_metrics.tsv", sep="\t").iloc[0]; made.append(latest("st_qc_target"))
        check(rc == 0 and int(t.n_guides_tested) == 18 and t.frac_strong_knockdowns == 1.0 and t.auroc > 0.95 and t.auprc > 0.95,
              f"qc-target: 18 targeting guides tested, all strong knockdowns, AUROC {t.auroc:.3f} / AUPRC {t.auprc:.3f}")

        print("\ncalibration and pathways")
        rc = cmd_calibrate(ns(trans_results=str(syn["perturbo"]), mudata=h5, prefix="syn", method="ecdf", cis_window=100_000,
                              non_targeting_label="non_targeting", one_sided=False, label="st_calibrate"))
        cal_dir = latest("st_calibrate"); made.append(cal_dir)
        allr = pd.read_csv(cal_dir / "syn_calibrated_all_results.tsv", sep="\t")
        tf1 = allr[allr.element_symbol == "TF1"]
        planted = tf1[tf1.tested_gene_id.isin([f"ENSG{100000 + i:08d}" for i in range(20)])]
        others = allr[(allr.element_label == "targeting") & ~allr.tested_gene_id.isin(planted.tested_gene_id)]
        check(rc == 0 and set(CALIBRATED_COLS) <= set(allr.columns), "calibrate: output has the upstream column set")
        check((planted.empirical_pval_adj < 0.05).mean() > 0.9 and (others.empirical_pval_adj < 0.05).mean() < 0.05,
              f"calibrate: planted TF1 effects {(planted.empirical_pval_adj < 0.05).mean():.0%} significant, others "
              f"{(others.empirical_pval_adj < 0.05).mean():.1%}")
        direct = pd.read_csv(cal_dir / "syn_calibrated_direct_target_results.tsv", sep="\t")
        cis = pd.read_csv(cal_dir / "syn_calibrated_cis_results.tsv", sep="\t")
        # genes sit 100 kb apart; midpoints: the upstream neighbour is 90 kb away (cis), the downstream one 110 kb (trans)
        check(len(direct) == 6 and direct.is_direct_target.all() and len(cis) == 11
              and int((cis.element_symbol == "TF1").sum()) == 1 and int((cis.element_symbol == "TF2").sum()) == 2,
              f"calibrate: 6 direct-target rows; cis = direct target + the neighbour whose midpoint is within 100 kb ({len(cis)} rows)")
        check(allr.loc[allr.element_label != "targeting", "empirical_pval_adj"].isna().all(),
              "calibrate: non-targeting rows are not in the BH discovery set")
        rc = cmd_calibrate(ns(trans_results=str(syn["perturbo"]), mudata=h5, prefix="synt", method="t-fit", cis_window=100_000,
                              non_targeting_label="non_targeting", one_sided=False, label="st_calibrate_t"))
        made.append(latest("st_calibrate_t"))
        check(rc == 0, "calibrate: t-fit method runs")
        rc = cmd_pathways(ns(calibrated=str(cal_dir / "syn_calibrated_all_results.tsv"), gmt=str(syn["gmt"]), prefix="syn", fdr=0.05,
                             min_genes=3, min_pathway_size=5, max_pathway_size=500, gene_col="tested_gene_symbol",
                             pval_col="empirical_pval_adj", label="st_pathways"))
        pw = pd.read_csv(latest("st_pathways") / "syn_calibrated_pathway_results.tsv", sep="\t"); made.append(latest("st_pathways"))
        top = pw.sort_values("pvalue").iloc[0] if len(pw) else None
        check(rc == 0 and top is not None and top.pathway_id == "PLANTED" and top.element_symbol == "TF1" and top.pvalue_adj < 1e-6,
              "pathways: the planted gene set is the top hit for TF1")

        print("\nfilters")
        rc = cmd_filter_cells(ns(mudata=h5, out=str(d / "filtered_cells.h5mu"), max_guides=15))
        outs = list(d.glob("filtered_cells*"))
        check(rc == 0 and outs, f"filter-cells: wrote {[p.name for p in outs]}")
        m2 = read_mudata(outs[0] if outs[0].suffix == ".h5mu" else d / "filtered_cells.gene.h5ad") if outs[0].suffix == ".h5mu" else None
        if m2 is not None:
            check(m2.gene.n_obs == syn["n_cells"] - 4, f"filter-cells: 4 cells with 21 guides removed ({m2.gene.n_obs} left)")
        else:
            import anndata as ad  # type: ignore
            a = ad.read_h5ad(str(d / "filtered_cells.gene.h5ad"))
            check(a.n_obs == syn["n_cells"] - 4, f"filter-cells: 4 cells with 21 guides removed ({a.n_obs} left; per-modality h5ad)")

        print("\nenergy distance")
        prep = d / "edist"
        rc = cmd_edist_prep(ns(mudata=h5, out_dir=str(prep), n_comps=20, label=""))
        check(rc == 0 and (prep / "pca_dataframe.pickle").is_file() and (prep / "annotation_table.csv").is_file(),
              "edist-prep: PCA, guide dictionary and annotation table written")
        rc = cmd_edist_filter(ns(prep_dir=str(prep), out_dir=None, min_cells=10, max_cells_per_guide=500, disco_permutations=60,
                                 disco_alpha=0.05, significance_fraction=0.5, fdr=0.05, seed=0))
        nt = pd.read_csv(prep / "non_targeting_outlier_table.csv", index_col=0)
        tg = pd.read_csv(prep / "targeting_outlier_table.csv", index_col=0)
        check(rc == 0 and nt.loc["NT#outlier", "pval_outlier"] == 0.0 and (nt.drop("NT#outlier")["pval_outlier"] == 1.0).all(),
              "edist-filter: the shifted non-targeting guide is the only K-means outlier")
        check(len(tg) == 18 and "pval_outlier" in tg.columns and tg["target"].nunique() == 6,
              "edist-filter: 18 sibling guides scored over 6 targets")
        rc = cmd_edist(ns(prep_dir=str(prep), out_dir=None, num_bg=3, permutations=60, non_target_pick=120, target_cell_max=200,
                          min_cells=10, fdr=0.05, keep_outliers=False, seed=0))
        ed = pd.read_csv(prep / "pval_edist_full.csv", index_col=0)
        tf1_row = ed[ed.index.str.startswith("ENSG00100000|")].iloc[0]
        rest = ed[~ed.index.str.startswith("ENSG00100000|")]
        check(rc == 0 and len(ed) == 6 and set(EDIST_REQUIRED) <= set(ed.columns) and "distance_2" in ed.columns,
              "edist: six targets, upstream column set incl. per-background distance_i / pval_i")
        check(tf1_row.pval_mean < 0.05 and tf1_row.distance_mean > rest.distance_mean.max(),
              f"edist: TF1 (planted) pval_mean {tf1_row.pval_mean:.3f}, distance {tf1_row.distance_mean:.2f} > others max {rest.distance_mean.max():.2f}")
        check((rest.pval_mean > 0.05).mean() >= 0.8, f"edist: unperturbed targets mostly non-significant ({(rest.pval_mean > 0.05).mean():.0%})")
        rc = cmd_edist_validate(ns(dir=str(prep)))
        check(rc == 0, "edist-validate: passes on the run folder")
        rc = cmd_edist_summary(ns(inputs=[f"syn={prep}"], label="st_edist_summary"))
        es = pd.read_csv(latest("st_edist_summary") / "cross_dataset_edistance_summary.tsv", sep="\t"); made.append(latest("st_edist_summary"))
        check(rc == 0 and int(es.n_targets_total.iloc[0]) == 6, "edist-summary: one row, six targets")
        rc = cmd_filter_guides(ns(mudata=h5, out=str(d / "filtered_guides.h5mu"), targeting_outliers=str(prep / "targeting_outlier_table.csv"),
                                  non_targeting_outliers=str(prep / "non_targeting_outlier_table.csv"), fdr=0.05))
        outs = list(d.glob("filtered_guides*"))
        check(rc == 0 and outs, "filter-guides: runs on the outlier tables")
        try:
            gmod = read_mudata(d / "filtered_guides.h5mu", need_uns=False).guide if (d / "filtered_guides.h5mu").is_file() else None
        except Exception:
            gmod = None
        if gmod is None:
            import anndata as ad  # type: ignore
            gmod = ad.read_h5ad(str(d / "filtered_guides.guide.h5ad"))
        check("NT#outlier" not in set(gmod.var_names) and gmod.n_obs <= syn["n_cells"] - syn["outlier_cells"] + 24,
              f"filter-guides: outlier guide and its cells removed ({gmod.n_vars} guides, {gmod.n_obs} cells)")

        print("\ncNMF export and summaries")
        rc = cmd_cnmf_export(ns(mudata=h5, out_h5ad=str(d / "perturbnmf.h5ad"), guide_id_col="guide_id", target_col="gene_name",
                                symbol_col="symbol", num_highvar_genes=50, no_hvg_filter=False))
        import anndata as ad  # type: ignore
        a = ad.read_h5ad(str(d / "perturbnmf.h5ad"))
        check(rc == 0 and "guide_assignment" in a.obsm and len(a.uns["guide_names"]) == syn["n_guides"] and a.var_names[0] == "TF1",
              f"cnmf-export: symbols as var_names, guide_assignment in obsm, {a.n_obs} cells kept")
        cn = d / "cnmf_run"; cn.mkdir()
        pd.DataFrame(np.random.rand(5, 30)).to_csv(cn / "run.gene_spectra_score.k_5.dt_2_0.txt", sep="\t")
        pd.DataFrame(np.random.rand(40, 5)).to_csv(cn / "run.usages.k_5.dt_2_0.consensus.txt", sep="\t")
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            v_ok, v_bad = cmd_cnmf_validate(ns(dir=str(cn), selected_k=5)), cmd_cnmf_validate(ns(dir=str(cn), selected_k=7))
        check(v_ok == 0 and v_bad == 1,
              "cnmf-validate: passes for the right k and fails for the wrong one")
        qc_root = d / "qc_bundle"; qc_root.mkdir()
        import shutil
        for src in (latest("st_qc_gene") / "gene_metrics.tsv", latest("st_qc_guide") / "guide_metrics.tsv",
                    latest("st_qc_target") / "intended_target_metrics.tsv"):
            shutil.copy(src, qc_root / src.name)
        rc = cmd_pipeline_summary(ns(datasets=[f"syn={qc_root}"], label="st_pipeline_summary"))
        ps = pd.read_csv(latest("st_pipeline_summary") / "cross_dataset_pipeline_summary.tsv", sep="\t"); made.append(latest("st_pipeline_summary"))
        check(rc == 0 and int(ps.n_cells.iloc[0]) == syn["n_cells"] and ps.intended_auroc.iloc[0] > 0.95,
              "pipeline-summary: one row joining the three QC tables")
    import shutil
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print("stages selftest: all checks pass" if ok else f"stages selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1
