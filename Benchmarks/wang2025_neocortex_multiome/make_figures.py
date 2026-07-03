#!/usr/bin/env python3
"""Score + plot the Wang 2025 developing-neocortex reproduction (sc-analyze)."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.metrics import (adjusted_rand_score, adjusted_mutual_info_score,
                             homogeneity_score, completeness_score, v_measure_score)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"; FIG.mkdir(parents=True, exist_ok=True)


def latest_run() -> Path:
    dirs = sorted((ROOT / "Docs/SingleCell").glob("2*_wang2025_neocortex_multiome"),
                  reverse=True)
    if not dirs:
        sys.exit("No sc-analyze run — run Benchmarks/wang2025_neocortex_multiome/run.sh first.")
    return dirs[0]


RUN = latest_run()
a = ad.read_h5ad(RUN / "processed.h5ad", backed="r")
leiden = a.obs["leiden"].astype(str).values
ct = a.obs["cell_type"].astype(str).values
a.file.close()

metrics = {
    "full_atlas_nuclei": 232328,
    "analyzed_subsample": int(len(ct)),
    "author_cell_types": int(len(set(ct))),
    "n_leiden_clusters": int(len(set(leiden))),
    "ARI_vs_author_types": round(float(adjusted_rand_score(ct, leiden)), 4),
    "AMI_vs_author_types": round(float(adjusted_mutual_info_score(ct, leiden)), 4),
    "homogeneity": round(float(homogeneity_score(ct, leiden)), 4),
    "completeness": round(float(completeness_score(ct, leiden)), 4),
    "v_measure": round(float(v_measure_score(ct, leiden)), 4),
}
(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
(RUN / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics, indent=2))

conf = pd.crosstab(pd.Series(ct, name="author cell type"),
                   pd.Series(leiden, name="Leiden"))
conf = conf.loc[conf.sum(1).sort_values(ascending=False).index]
conf = conf[[conf.columns[i] for i in np.argsort(conf.values.argmax(0))]]
confn = conf.div(conf.sum(1), axis=0)
fig, ax = plt.subplots(figsize=(13, 9), facecolor="white")
im = ax.imshow(confn.values, aspect="auto", cmap="magma_r", vmin=0, vmax=1)
ax.set_xticks(range(conf.shape[1])); ax.set_xticklabels(conf.columns, fontsize=5, rotation=90)
ax.set_yticks(range(conf.shape[0])); ax.set_yticklabels(conf.index, fontsize=7)
ax.set_xlabel("IGVFagent Leiden cluster"); ax.set_ylabel("Wang 2025 author cell type")
ax.set_title(f"Wang 2025 neocortex multiome — {metrics['n_leiden_clusters']} Leiden clusters "
             f"vs {metrics['author_cell_types']} author cell types\n"
             f"AMI={metrics['AMI_vs_author_types']} · homogeneity={metrics['homogeneity']} "
             f"(20k-cell subsample of {metrics['full_atlas_nuclei']:,})",
             fontweight="bold", fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.025, label="fraction of author type in cluster")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_confusion.{ext}", dpi=160, facecolor="white")
plt.close(fig)
print(f"figures under {FIG.relative_to(ROOT)}")
