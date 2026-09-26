#!/usr/bin/env python3
"""Figures for the LLCB Python-port causal network (edges.csv)."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROC = ROOT / "Data/Weinstock2024/processed"
FIG_DIR = ROOT / "Benchmarks/weinstock2024_cd4_crispr/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"
PAPER_POINTS = [(0.020, 350), (0.025, 211), (0.030, 151)]  # paper Fig 2/Table S1

edges = pd.read_csv(PROC / "edges.csv")

# ----- Fig 4: edge count vs. magnitude threshold -----
thresholds = np.round(np.arange(0.01, 0.16, 0.005), 3)
counts = [(edges["estimate"].abs() > t).sum() for t in thresholds]

fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
ax.plot(thresholds, counts, color=COL_PRIMARY, lw=2, marker="o", ms=3, label="this port (GSE271788, own VST/PCs)")
px, py = zip(*PAPER_POINTS)
ax.plot(px, py, color=COL_HIGHLIGHT, lw=0, marker="D", ms=8, zorder=5,
        label="paper Fig 2/Table S1 (350, 211, 151)")
for x, y in PAPER_POINTS:
    ax.annotate(f"{y}", (x, y), textcoords="offset points", xytext=(8, -3),
                color=COL_HIGHLIGHT, fontsize=9, fontweight="bold")
ax.set_xlabel(r"edge-calling threshold on $|\hat\beta|$")
ax.set_ylabel("recovered directed edges (of 6,972 possible pairs)")
ax.set_title("Python-port LLCB: recovered edge count vs. threshold\n"
              "same 84-gene panel and thresholds as the paper — our network runs "
              "~3–4× denser at any given |β| cutoff",
              fontweight="bold", fontsize=10)
ax.legend(fontsize=9, frameon=False)
ax.grid(ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig4_llcb_edge_threshold.{ext}", dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig4_llcb_edge_threshold")

# ----- Fig 5: top |estimate| edges, Th17/IL2/JAK-STAT axis genes highlighted -----
AXIS_GENES = {"STAT5A", "STAT5B", "IL2RA", "JAK3", "RORC", "IRF4", "RELA", "STAT3", "STAT1", "STAT2", "KMT2A"}
top = edges.reindex(edges["estimate"].abs().sort_values(ascending=False).index).head(25)
labels = [f"{r} → {c}" for r, c in zip(top["row"], top["col"])]
colors = [
    COL_HIGHLIGHT if (r in AXIS_GENES or c in AXIS_GENES) else COL_PRIMARY
    for r, c in zip(top["row"], top["col"])
]

fig, ax = plt.subplots(figsize=(8, 8), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, top["estimate"], color=colors, edgecolor="white", linewidth=0.5)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel(r"posterior mean direct effect $\hat\beta$")
ax.set_title("Top 25 recovered edges by |effect size|\n"
              "★ = touches a paper-highlighted Th17/IL2/JAK-STAT axis gene",
              fontweight="bold", fontsize=10)
ax.grid(axis="x", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig5_llcb_top_edges.{ext}", dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig5_llcb_top_edges")

print(f"\nFigures saved under {FIG_DIR}")
