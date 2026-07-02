#!/usr/bin/env python3
"""Score + plot the Trevino 2021 multiome peak->gene reproduction.

Reads the peak2gene TSV written by `igvfagent multiome peak2gene` on the
cortical-lineage panel and produces:

  * concordance_metrics.json  — link counts, positive fraction, per-gene table
  * fig1_peak2gene.png/.svg   — per-gene significant cis-links + best-peak corr
  * fig2_distance.png/.svg    — correlation vs distance-to-TSS (cis enrichment)
"""
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
D = ROOT / "Benchmarks/_data/trevino2021_cortex_multiome"


def latest_tsv() -> Path:
    hits = sorted((ROOT / "Docs/Multiome10x").glob("2*trevino2021_cortex_multiome_peak2gene.tsv"),
                  reverse=True)
    if not hits:
        sys.exit("No peak2gene TSV — run Benchmarks/trevino2021_cortex_multiome/run.sh first.")
    return hits[0]


TSV = latest_tsv()
df = pd.read_csv(TSV, sep="\t")
ens2sym = json.load(open(D / "panel_tss_grch38.json"))
i2s = {v["ensembl"]: s for s, v in ens2sym.items()}
df["symbol"] = df["gene"].map(i2s)

sig = df[df["padj"] < 0.05]
pos = sig[sig["correlation"] > 0]

per_gene = (pos.groupby("symbol")
              .agg(n_sig_pos=("correlation", "size"),
                   best_corr=("correlation", "max"),
                   best_dist=("distance", lambda s: int(s.loc[s.abs().idxmin()])))
              .sort_values("n_sig_pos", ascending=False))

metrics = {
    "panel_genes": int(df["gene"].nunique()),
    "peaks_tested": int(df["peak"].nunique()),
    "candidate_pairs": int(len(df)),
    "significant_pairs": int(len(sig)),
    "significant_positive_pairs": int(len(pos)),
    "positive_fraction_of_significant": round(len(pos) / max(1, len(sig)), 4),
    "median_abs_distance_sig_pos_bp": int(pos["distance"].abs().median()),
    "genes_with_significant_positive_link": int(per_gene.shape[0]),
    "per_gene": per_gene.reset_index().to_dict(orient="records"),
}
(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
# also drop into the skill output dir so concordance.py (flat-file convention) finds it
(ROOT / "Docs/Multiome10x/concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps({k: v for k, v in metrics.items() if k != "per_gene"}, indent=2))

# fig1: per-gene significant positive links + best correlation
g = per_gene.head(26)
fig, ax = plt.subplots(figsize=(11, 6), facecolor="white")
bars = ax.bar(range(len(g)), g["n_sig_pos"], color="#4C78A8", edgecolor="white")
ax.set_xticks(range(len(g))); ax.set_xticklabels(g.index, rotation=60, ha="right", fontsize=8)
ax.set_ylabel("significant positive cis peak→gene links (padj<0.05)")
ax2 = ax.twinx()
ax2.plot(range(len(g)), g["best_corr"], "o-", color="#E45756", ms=5, lw=1)
ax2.set_ylabel("best peak→gene correlation", color="#E45756")
ax2.tick_params(axis="y", labelcolor="#E45756")
ax.set_title(f"Trevino 2021 multiome — cis peak→gene linkage on a {metrics['panel_genes']}-gene cortical panel\n"
             f"{metrics['significant_positive_pairs']:,} significant positive links "
             f"({100*metrics['positive_fraction_of_significant']:.0f}% of all significant)",
             fontweight="bold", fontsize=11)
for s in ("top",): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_peak2gene.{ext}", dpi=170, facecolor="white")
plt.close(fig); print("  ✓ fig1_peak2gene")

# fig2: correlation vs distance to TSS (cis enrichment near 0)
fig, ax = plt.subplots(figsize=(9, 5), facecolor="white")
ax.scatter(df["distance"] / 1000, df["correlation"], s=4, alpha=0.15, color="#999999", label="all tested")
ax.scatter(pos["distance"] / 1000, pos["correlation"], s=8, alpha=0.5, color="#E45756",
           label="significant positive")
ax.axvline(0, color="#1F2933", ls=":", lw=0.8)
ax.set_xlabel("peak midpoint − gene TSS (kb)")
ax.set_ylabel("Spearman correlation (ATAC peak vs gene RNA)")
ax.set_title("Cis peak→gene correlation concentrates near the TSS", fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig2_distance.{ext}", dpi=170, facecolor="white")
plt.close(fig); print("  ✓ fig2_distance")
print(f"\nFigures + metrics under {FIG.relative_to(ROOT)}")
