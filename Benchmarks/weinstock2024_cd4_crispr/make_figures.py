#!/usr/bin/env python3
"""Regenerate figures embedded in Benchmarks/weinstock2024_cd4_crispr/README.md.

Reads the Perturbation Catalogue search JSON (CRISPR-screen modality)
plus the Weinstock GSE171674 GEO report to render publication panels.
"""
import csv
import glob
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/weinstock2024_cd4_crispr/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

search_json = sorted(
    (ROOT / "Data/Perturbation/Searches").glob("*_modality_crispr-screen_KMT2A.json")
)
if not search_json:
    sys.exit("No KMT2A crispr-screen search JSON found.")
d = json.load(search_json[-1].open())
total = d.get("total_datasets_count", 0)
facets = d.get("facet_counts", {}) or {}


def topn(field, n=10):
    rows = facets.get(field, [])
    if not isinstance(rows, list):
        return [], []
    rows = sorted(rows, key=lambda r: -r["count"])[:n]
    return [r["value"] for r in rows], [r["count"] for r in rows]


print(f"Total CRISPR-screen datasets w/ KMT2A: {total:,}")

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"

# ----- Fig 1: top tissues -----
HIGHLIGHT_T = {"lymphoid tissue", "blood", "bone marrow"}
labels, counts = topn("dataset_tissues", 12)
colors = [COL_HIGHLIGHT if l in HIGHLIGHT_T else COL_PRIMARY for l in labels]
fig, ax = plt.subplots(figsize=(9, 5), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    suffix = "  ★ T-cell relevant" if l in HIGHLIGHT_T else ""
    ax.text(c + 2, y[i], f"  n={c}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l in HIGHLIGHT_T else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("CRISPR-screen datasets (Perturbation Catalogue)")
ax.set_title(f"Top tissues across all {total:,} CRISPR-screen datasets\n"
              "Lymphoid tissue + blood + bone marrow (★) are the T-cell context",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
ax.set_xlim(0, max(counts) * 1.25)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_top_tissues.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_top_tissues")

# ----- Fig 2: top cell types -----
HIGHLIGHT_C = {"t cell", "b cell", "plasma cell"}
labels, counts = topn("dataset_cell_types", 10)
colors = [COL_HIGHLIGHT if l in HIGHLIGHT_C else COL_PRIMARY for l in labels]
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    suffix = "  ★ adaptive-immune lineage" if l in HIGHLIGHT_C else ""
    ax.text(c + 0.2, y[i], f"  n={c}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l in HIGHLIGHT_C else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("CRISPR-screen datasets")
ax.set_title("Top cell types — adaptive-immune lineages (★) "
              "are the Weinstock 2024 context",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
ax.set_xlim(0, max(counts) * 1.3)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_cell_types.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_cell_types")

# ----- Fig 3: perturbation type breakdown (knockout dominant) -----
labels, counts = topn("dataset_perturbation_types", 5)
fig, ax = plt.subplots(figsize=(7, 4), facecolor="white")
bars = ax.bar(labels, counts, color=COL_PRIMARY, edgecolor="white",
                linewidth=0.5)
for b, c in zip(bars, counts):
    pct = c / total * 100
    ax.text(b.get_x() + b.get_width() / 2,
            c + max(counts) * 0.02,
            f"{c:,}\n({pct:.1f} %)",
            ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("Dataset count")
ax.set_title("Perturbation modality — CRISPR-KO (CRISPRn) dominates\n"
              "matching Weinstock 2024's 84-gene knockout design",
              fontweight="bold", fontsize=11)
ax.set_ylim(0, max(counts) * 1.15)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig3_perturbation_types.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig3_perturbation_types")

print(f"\nFigures saved under {FIG_DIR}")
