#!/usr/bin/env python3
"""Regenerate figures embedded in Benchmarks/mitra2024_scarlink/README.md.

Reads the IGVF Portal multiome manifest written by `run.sh` and produces
publication-grade matplotlib panels.
"""
import csv
import glob
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/mitra2024_scarlink/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Find the most recent mitra2024_scarlink analysis-set manifest
candidates = sorted((ROOT / "Data/Manifests/Multiome10x").glob(
    "*_mitra2024_scarlink_analysis_sets.csv"))
file_candidates = sorted((ROOT / "Data/Manifests/Multiome10x").glob(
    "*_mitra2024_scarlink_files.csv"))
if not candidates or not file_candidates:
    sys.exit("No Mitra 2024 manifest found — run "
              "`bash Benchmarks/mitra2024_scarlink/run.sh` first.")

rows = list(csv.DictReader(candidates[-1].open()))
files = list(csv.DictReader(file_candidates[-1].open()))
print(f"Loaded {len(rows)} AnalysisSets, {len(files)} files from "
       f"{candidates[-1].name}")

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"

# ----- Fig 1: AnalysisSets by sample-term (brain region) -----
samples = Counter(r["sample_terms"] for r in rows)
labels = list(samples.keys())
counts = list(samples.values())
sizes = [float(r["total_file_size_gb"]) for r in rows]

fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
labels_with_size = []
for r in rows:
    s = r["sample_terms"].replace(";", " /").strip()
    labels_with_size.append(f"{s}\n{float(r['total_file_size_gb']):.2f} GB")

y = np.arange(len(rows))[::-1]
sizes_gb = [float(r["total_file_size_gb"]) for r in rows]
ax.barh(y, sizes_gb, color=COL_PRIMARY, edgecolor="white", linewidth=0.5)
for i, (l, s) in enumerate(zip(labels_with_size, sizes_gb)):
    ax.text(s + 0.05, y[i], f"  {s:.2f} GB", va="center", fontsize=9)
ax.set_yticks(y)
ax.set_yticklabels(labels_with_size, fontsize=9)
ax.set_xlabel("Total file size (GB)")
ax.set_title("IGVF Portal 10x multiome AnalysisSets selected for Mitra-style analysis\n"
              "5 / 505 principal-analysis multiome AnalysisSets (Corces/Gladstone)",
              fontweight="bold", fontsize=11)
ax.set_xlim(0, max(sizes_gb) * 1.30)
ax.grid(axis="x", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_analysis_sets.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_analysis_sets")

# ----- Fig 2: file content type breakdown -----
# The auto-report says content_types = {annotated cell by peak matrix: 5,
# cell by gene matrix: 5, cell annotations: 5, fragments: 5}
content_types = Counter(f["content_type"] for f in files
                          if f.get("content_type"))
labels = [k for k, _ in content_types.most_common()]
counts = [v for _, v in content_types.most_common()]
fig, ax = plt.subplots(figsize=(8, 4), facecolor="white")
bars = ax.bar(labels, counts, color=COL_PRIMARY, edgecolor="white",
                linewidth=0.5)
for b, c in zip(bars, counts):
    ax.text(b.get_x() + b.get_width() / 2, c + 0.1, f"{c}",
            ha="center", fontsize=11, fontweight="bold")
ax.set_ylabel("File count across 5 AnalysisSets")
ax.set_title("Content-type completeness — every AnalysisSet ships "
              "the full SCARlink-style payload\n"
              "(peak matrix + gene matrix + cell annotations + ATAC fragments)",
              fontweight="bold", fontsize=11)
ax.set_ylim(0, max(counts) * 1.20)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for t in ax.get_xticklabels():
    t.set_rotation(15); t.set_ha("right")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_content_types.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_content_types")

# ----- Fig 3: Donut — Portal scale + Corces share -----
fig, ax = plt.subplots(figsize=(7, 5), facecolor="white")
PORTAL_TOTAL = 505     # from the live API summary
CORCES_SELECTED = 5
other = PORTAL_TOTAL - CORCES_SELECTED
sizes = [CORCES_SELECTED, other]
colors = [COL_HIGHLIGHT, COL_PRIMARY]
labels = [f"Corces / Gladstone\nbenchmark slice ({CORCES_SELECTED})",
            f"Remaining IGVF multiome\nAnalysisSets ({other:,})"]
wedges, texts = ax.pie(sizes, labels=labels, colors=colors,
                          startangle=90, wedgeprops=dict(width=0.4),
                          textprops=dict(fontsize=10))
ax.set_title(f"Portal scale: {PORTAL_TOTAL} released public "
              "10x-multiome AnalysisSets\n"
              "(IGVFagent enumerates the full universe, "
              "this benchmark slices 5 for plotting)",
              fontweight="bold", fontsize=11)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig3_portal_scale.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig3_portal_scale")

print(f"\nFigures saved under {FIG_DIR}")
