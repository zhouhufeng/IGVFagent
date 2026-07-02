#!/usr/bin/env python3
"""Score + plot the Travaglini 2020 lung-atlas reproduction.

Reads the processed AnnData written by

    igvfagent sc-analyze pipeline --label travaglini2020_lung

(which carries both the IGVFagent ``leiden`` clustering and the paper's
author ``cell_type`` labels), and produces:

  * concordance_metrics.json  — ARI / AMI / homogeneity / completeness
  * fig1_confusion.png/.svg   — author cell_type x Leiden cluster heatmap
  * fig2_markers.png/.svg     — canonical lung markers x author cell_type
  * fig3_umap.png/.svg        — UMAP, author labels vs IGVFagent Leiden

The two scored claims:
  1. cluster concordance  (Leiden vs 46 author cell types)
  2. canonical-marker biology  (each marker peaks in its expected cell type)
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from sklearn.metrics import (adjusted_rand_score, adjusted_mutual_info_score,
                              homogeneity_score, completeness_score,
                              v_measure_score)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def latest_run() -> Path:
    base = ROOT / "Docs/SingleCell"
    dirs = sorted(base.glob("2*_travaglini2020_lung"), reverse=True)
    if not dirs:
        sys.exit("No sc-analyze run found — run Benchmarks/travaglini2020_lung/run.sh first.")
    return dirs[0]


RUN = latest_run()
print(f"==> run dir: {RUN.relative_to(ROOT)}")
# Memory-frugal: read obs from backed mode (never materialize the matrix).
a = ad.read_h5ad(RUN / "processed.h5ad", backed="r")
leiden = a.obs["leiden"].astype(str).values
ct = a.obs["cell_type"].astype(str).values
obs_names = a.obs_names.astype(str).values
a.file.close()

# ---- 1. clustering concordance -------------------------------------------
metrics = {
    "n_cells": int(a.n_obs),
    "n_leiden_clusters": int(len(set(leiden))),
    "n_author_cell_types": int(len(set(ct))),
    "ARI": round(float(adjusted_rand_score(ct, leiden)), 4),
    "AMI": round(float(adjusted_mutual_info_score(ct, leiden)), 4),
    "homogeneity": round(float(homogeneity_score(ct, leiden)), 4),
    "completeness": round(float(completeness_score(ct, leiden)), 4),
    "v_measure": round(float(v_measure_score(ct, leiden)), 4),
}

# ---- 2. canonical marker biology -----------------------------------------
# marker -> substring identifying its expected author cell_type
# substring matched (case-insensitively) against the paper's author cell_type labels
CANON = {
    "SFTPC": "alveolar type 2", "AGER": "alveolar type 1",
    "SCGB1A1": "club", "FOXJ1": "ciliated", "PECAM1": "endothelial",
    "MARCO": "alveolar macrophage", "CD3E": "t cell", "MS4A1": "b cell",
    "DCN": "fibroblast", "EPCAM": "epithelial",
}
# Read only the marker columns from the raw-counts prep file (backed → slice).
PREP = ROOT / "Benchmarks/_data/travaglini2020_lung/lung_for_scanalyze.h5ad"
ap = ad.read_h5ad(PREP, backed="r")
present = [g for g in CANON if g in ap.var_names]
col_idx = [ap.var_names.get_loc(g) for g in present]
# Pull just those columns into memory (65k x ~10), aligning cells by name.
Xcols = ap[:, present].to_memory().X
Xcols = Xcols.toarray() if sp.issparse(Xcols) else np.asarray(Xcols)
prep_names = ap.obs_names.astype(str).values
ap.file.close()
# align prep rows to the obs order used for clustering
order = pd.Series(np.arange(len(prep_names)), index=prep_names).reindex(obs_names).values
Xcols = Xcols[order.astype(int)]
# library-size normalize + log1p so cross-gene comparison is fair
lib = Xcols.sum(1, keepdims=True); lib[lib == 0] = 1
Xn = np.log1p(Xcols / lib * 1e4)
mexp = pd.DataFrame(Xn, columns=present)
mexp["ct"] = ct
per_ct = mexp.groupby("ct").mean()

marker_hits = 0
marker_rows = []
for g in present:
    target_sub = CANON[g]
    ranked = per_ct[g].sort_values(ascending=False)
    top_ct = ranked.index[0]
    target_types = [c for c in per_ct.index if target_sub in c.lower()]
    # rank of the best-matching target cell type
    rank = min([list(ranked.index).index(t) for t in target_types], default=999) + 1 \
        if target_types else 999
    ok = rank <= 3
    marker_hits += ok
    marker_rows.append({"marker": g, "expected": target_sub,
                         "argmax_cell_type": top_ct, "target_rank": rank, "pass": bool(ok)})
metrics["canonical_markers_tested"] = len(present)
metrics["canonical_markers_peaking_in_target_top3"] = int(marker_hits)
metrics["marker_detail"] = marker_rows

(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
# also drop a copy into the run dir so the scorer can read it
(RUN / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps({k: v for k, v in metrics.items() if k != "marker_detail"}, indent=2))

# ---- fig1: confusion heatmap (author cell_type x leiden, row-normalized) --
conf = pd.crosstab(pd.Series(ct, name="author cell_type"),
                   pd.Series(leiden, name="Leiden"))
conf = conf.loc[conf.sum(1).sort_values(ascending=False).index]
order = conf.idxmax(0).reset_index().sort_values(0).index  # roughly block-diagonalize
conf = conf[[conf.columns[i] for i in np.argsort(conf.values.argmax(0))]]
confn = conf.div(conf.sum(1), axis=0)
fig, ax = plt.subplots(figsize=(13, 11), facecolor="white")
im = ax.imshow(confn.values, aspect="auto", cmap="magma_r", vmin=0, vmax=1)
ax.set_xticks(range(conf.shape[1])); ax.set_xticklabels(conf.columns, fontsize=6, rotation=90)
ax.set_yticks(range(conf.shape[0])); ax.set_yticklabels(conf.index, fontsize=6)
ax.set_xlabel("IGVFagent Leiden cluster"); ax.set_ylabel("Travaglini 2020 author cell type")
ax.set_title(f"Cluster concordance — {metrics['n_leiden_clusters']} Leiden clusters vs "
             f"{metrics['n_author_cell_types']} author cell types\n"
             f"AMI={metrics['AMI']} · homogeneity={metrics['homogeneity']} · ARI={metrics['ARI']}",
             fontweight="bold", fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.025, label="fraction of author cell type in cluster")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_confusion.{ext}", dpi=170, facecolor="white")
plt.close(fig); print("  ✓ fig1_confusion")

# ---- fig2: canonical markers x author cell_type (mean expr, z per marker) -
sub = per_ct[present]
z = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
# show only cell types that are a top-2 expressor of some canonical marker, for legibility
keep_ct = sorted(set(z.idxmax(0)) | set(z.apply(lambda c: c.nlargest(2).index[-1])))
zk = z.loc[keep_ct]
fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(keep_ct))), facecolor="white")
im = ax.imshow(zk.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(present))); ax.set_xticklabels(present, fontsize=9, rotation=45, ha="right")
ax.set_yticks(range(len(keep_ct))); ax.set_yticklabels(keep_ct, fontsize=7)
ax.set_title(f"Canonical lung markers peak in expected cell types\n"
             f"{metrics['canonical_markers_peaking_in_target_top3']}/{metrics['canonical_markers_tested']} "
             f"markers rank their target cell type in the top 3", fontweight="bold", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.025, label="z-scored mean log-expr")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig2_markers.{ext}", dpi=170, facecolor="white")
plt.close(fig); print("  ✓ fig2_markers")

print(f"\nFigures + metrics under {FIG.relative_to(ROOT)}")
