#!/usr/bin/env python3
"""Score + plot the Rosenberg 2018 SPLiT-seq CNS reproduction.

Reads the processed AnnData from `igvfagent splitseq analyze` and writes:
  * concordance_metrics.json  — full-atlas N, subsample N, clusters, CNS cell types
  * fig1_celltypes.png/.svg   — recovered CNS cell-type composition
"""
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
PROC = ROOT / "Data/Cache/SPLiTseq/rosenberg2018_splitseq_processed.h5ad"

if not PROC.is_file():
    sys.exit(f"missing {PROC} — run Benchmarks/rosenberg2018_splitseq/run.sh first")

a = ad.read_h5ad(PROC)
ct = a.obs["cell_type"].value_counts()
# canonical CNS lineages Rosenberg 2018 resolves
CNS = {"Excitatory neuron", "Inhibitory neuron", "Astrocyte",
       "Oligodendrocyte", "OPC", "Microglia", "Endothelial", "Pericyte"}
recovered = sorted(set(ct.index) & CNS)

metrics = {
    "full_atlas_nuclei": int(a.obs["full_atlas_n"].iloc[0]) if "full_atlas_n" in a.obs else None,
    "analyzed_nuclei": int(a.n_obs),
    "genes": int(a.n_vars),
    "n_leiden_clusters": int(a.obs["leiden"].nunique()),
    "n_cell_types": int(ct.shape[0]),
    "cns_lineages_recovered": recovered,
    "n_cns_lineages_recovered": len(recovered),
    "cell_type_counts": ct.to_dict(),
}
(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
# also drop into a Docs/SPLiTseq run dir so concordance.py discovers it
rundir = ROOT / "Docs/SPLiTseq" / f"{time.strftime('%Y%m%d_%H%M%S')}_rosenberg2018_splitseq"
rundir.mkdir(parents=True, exist_ok=True)
(rundir / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps({k: v for k, v in metrics.items() if k != "cell_type_counts"},
                 indent=2))

fig, ax = plt.subplots(figsize=(9, 5), facecolor="white")
ax.bar(ct.index, ct.values, color="#4C78A8", edgecolor="white")
ax.set_ylabel("nuclei"); ax.set_xticklabels(ct.index, rotation=35, ha="right", fontsize=9)
ax.set_title(f"Rosenberg 2018 SPLiT-seq CNS — {metrics['n_cell_types']} cell types / "
             f"{metrics['n_leiden_clusters']} Leiden clusters\n"
             f"(subsample of the {metrics['full_atlas_nuclei']:,}-nucleus atlas)",
             fontweight="bold", fontsize=11)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_celltypes.{ext}", dpi=170, facecolor="white")
plt.close(fig)
print(f"figures under {FIG.relative_to(ROOT)}")
