#!/usr/bin/env python3
"""Regenerate the figures embedded in Benchmarks/yao2024_encode4_crispri/README.md.

Pulls live data from www.encodeproject.org via IGVFagent's encode_pipeline
(through the on-disk manifest CSV the run.sh writes) and renders three
publication-grade panels under figures/.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/yao2024_encode4_crispri/figures"
MANIFEST_DIR = ROOT / "Data/Manifests/ENCODE"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# Find the most recent yao2024_encode4 manifest
candidates = sorted(MANIFEST_DIR.glob("*_yao2024_encode4_experiments.csv"))
if not candidates:
    sys.exit("No Yao 2024 manifest found under Data/Manifests/ENCODE/ — "
              "run `bash Benchmarks/yao2024_encode4_crispri/run.sh` first.")
manifest = candidates[-1]
rows = list(csv.DictReader(manifest.open()))
print(f"Loaded {len(rows):,} CRISPR-screen FCEs from {manifest.name}")

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"
COL_GREEN = "#7CA663"

# ----- Fig 1: screen count by biosample -----
bios = Counter(r["biosample"] for r in rows).most_common(12)
labels = [b for b, _ in bios]
counts = [c for _, c in bios]
HL = {"K562", "HepG2", "Jurkat"}
colors = [COL_HIGHLIGHT if l in HL else COL_PRIMARY for l in labels]
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    suffix = "  ★ Yao 2024 cell line" if l in HL else ""
    ax.text(c + 2, y[i], f"  n={c}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l in HL else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("ENCODE CRISPR-screen FCEs")
ax.set_title("ENCODE FunctionalCharacterizationExperiment "
              "CRISPR screens by biosample\n"
              "K562 / HepG2 / Jurkat (★) = Yao 2024's three cell lines",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
ax.set_xlim(0, max(counts) * 1.20)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_screens_by_biosample.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_screens_by_biosample")

# ----- Fig 2: screen count by assay_title (readout) -----
assays = Counter(r["assay_title"] for r in rows).most_common()
labels = [a for a, _ in assays]
counts = [c for _, c in assays]
fig, ax = plt.subplots(figsize=(9, 3.8), facecolor="white")
bars = ax.bar(labels, counts, color=COL_PRIMARY,
                edgecolor="white", linewidth=0.5)
# Highlight Flow-FISH (paper Fig 1 headline readout)
for b, l, c in zip(bars, labels, counts):
    if l == "Flow-FISH CRISPR screen":
        b.set_color(COL_HIGHLIGHT)
    ax.text(b.get_x() + b.get_width() / 2, c + 4, f"{c}",
            ha="center", fontsize=10,
            fontweight="bold" if l == "Flow-FISH CRISPR screen" else "normal")
ax.set_ylabel("ENCODE FCE count")
ax.set_title("CRISPR-screen readout breakdown — Flow-FISH dominates\n"
              "(matches Yao 2024 Fig 1: Flow-FISH is the paper's primary readout)",
              fontweight="bold", fontsize=11)
ax.set_ylim(0, max(counts) * 1.15)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for t in ax.get_xticklabels():
    t.set_rotation(15)
    t.set_ha("right")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_screens_by_assay.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_screens_by_assay")

# ----- Fig 3: top depositing labs (Engreitz dominance) -----
labs = Counter(r["lab"] for r in rows).most_common(8)
labels = [l for l, _ in labs]
counts = [c for _, c in labs]
ENGREITZ = "Jesse Engreitz, Stanford"
engreitz_n = next((c for l, c in labs if l == ENGREITZ), 0)
engreitz_pct = round(engreitz_n / len(rows) * 100)
colors = [COL_HIGHLIGHT if l == ENGREITZ else COL_PRIMARY for l in labels]
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    suffix = "  ★ Yao 2024 senior-author lab" if l == ENGREITZ else ""
    ax.text(c + 2, y[i], f"  n={c}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l == ENGREITZ else "normal")
ax.set_yticks(y)
ax.set_yticklabels([l.split(",")[0] for l in labels], fontsize=10)
ax.set_xlabel("ENCODE CRISPR-screen FCEs deposited")
ax.set_title("Depositing labs — Engreitz (Yao 2024 senior author) is dominant\n"
              f"{engreitz_n} / {len(rows)} ({engreitz_pct} %) of all ENCODE "
              "CRISPR-screen FCEs",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
ax.set_xlim(0, max(counts) * 1.25)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig3_top_labs.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig3_top_labs")

print(f"\nFigures saved under {FIG_DIR}")
