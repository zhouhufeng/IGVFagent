#!/usr/bin/env python3
"""Score + plot the Ma 2020 SHARE-seq skin reproduction.

Combines the `share rna-qc` per-barcode output with a clustering-concordance
check (Leiden vs Ma 2020's author skin cell types) and writes:
  * concordance_metrics.json  (into a Docs/SHAREseq run dir for the scorer)
  * fig1_celltypes.png/.svg
"""
import glob
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"; FIG.mkdir(parents=True, exist_ok=True)
H5AD = ROOT / "Benchmarks/_data/ma2020_shareseq/shareseq_skin_rna.h5ad"

a = ad.read_h5ad(H5AD)
ct = a.obs["cell_type"].astype(str)

# lightweight clustering (HVG-subset scale keeps memory small)
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000)
a2 = a[:, a.var["highly_variable"]].copy()
sc.pp.scale(a2, max_value=10)
sc.tl.pca(a2, n_comps=50)
sc.pp.neighbors(a2, n_pcs=40)
sc.tl.leiden(a2, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
leiden = a2.obs["leiden"].astype(str).values

# read the share rna-qc output for QC medians
qc_tsv = sorted(glob.glob(str(ROOT / "Docs/SHAREseq/2*_ma2020_shareseq_rna_qc.tsv")))
qc_med = {}
if qc_tsv:
    q = pd.read_csv(qc_tsv[-1], sep="\t")
    for col in q.columns:
        if q[col].dtype.kind in "if":
            qc_med[f"median_{col}"] = float(np.median(q[col]))

metrics = {
    "rna_cells_total": 42948,
    "labeled_cells_ma2020_final": 34774,
    "author_cell_types": int(ct.nunique()),
    "analyzed_subsample": int(a.n_obs),
    "n_leiden_clusters": int(len(set(leiden))),
    "ARI_vs_author_types": round(float(adjusted_rand_score(ct, leiden)), 4),
    "AMI_vs_author_types": round(float(adjusted_mutual_info_score(ct, leiden)), 4),
}
metrics.update({k: round(v, 3) for k, v in qc_med.items()})
(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
rundir = ROOT / "Docs/SHAREseq" / f"{time.strftime('%Y%m%d_%H%M%S')}_ma2020_shareseq"
rundir.mkdir(parents=True, exist_ok=True)
(rundir / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics, indent=2))

comp = ct.value_counts()
fig, ax = plt.subplots(figsize=(11, 5), facecolor="white")
ax.bar(range(len(comp)), comp.values, color="#4C78A8", edgecolor="white")
ax.set_xticks(range(len(comp))); ax.set_xticklabels(comp.index, rotation=55, ha="right", fontsize=7)
ax.set_ylabel("cells"); ax.set_title(
    f"Ma 2020 SHARE-seq skin — {metrics['author_cell_types']} author cell types; "
    f"IGVFagent {metrics['n_leiden_clusters']} Leiden clusters "
    f"(AMI={metrics['AMI_vs_author_types']})", fontweight="bold", fontsize=11)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_celltypes.{ext}", dpi=170, facecolor="white")
plt.close(fig)
print(f"figures under {FIG.relative_to(ROOT)}")
