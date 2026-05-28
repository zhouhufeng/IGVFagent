#!/usr/bin/env python3
"""Regenerate figures embedded in Benchmarks/zheng2024_invivo_perturbseq/README.md.

Reads the GSE249416 GEO file manifest written by `run.sh` and renders
panels summarising the published Zheng 2024 in-vivo Perturb-seq cohort.
"""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/zheng2024_invivo_perturbseq/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

files_csv = sorted((ROOT / "Data/Manifests/GEO").glob("*_GSE249416_files.csv"))
if not files_csv:
    sys.exit("No GSE249416 manifest found.")
files = list(csv.DictReader(files_csv[-1].open()))
print(f"Loaded {len(files)} files from {files_csv[-1].name}")

# Read the GEO report markdown to extract the platform / sample fields
report = sorted((ROOT / "Docs/GEO").glob("*_GSE249416_geo_report.md"))[-1].read_text()


def grep_line(label):
    m = re.search(rf"\*\*{re.escape(label)}\*\*: (.+)", report)
    return m.group(1) if m else None


n_samples = int(grep_line("n_samples") or 0)
platforms = (grep_line("platforms") or "").split(";")
pubmed = grep_line("pubmed_id")
print(f"Samples: {n_samples}, platforms: {platforms}, pubmed: {pubmed}")

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"
COL_GREEN = "#7CA663"

# ----- Fig 1: file category breakdown -----
cats = Counter(f["category"] for f in files)
labels = list(cats.keys())
counts = list(cats.values())
fig, ax = plt.subplots(figsize=(7, 3.8), facecolor="white")
bars = ax.bar(labels, counts, color=COL_PRIMARY, edgecolor="white",
                linewidth=0.5)
for b, c in zip(bars, counts):
    ax.text(b.get_x() + b.get_width() / 2, c + 0.1, f"{c}",
            ha="center", fontsize=11, fontweight="bold")
ax.set_ylabel("File count")
ax.set_title("GSE249416 file inventory (Zheng 2024)\n"
              "matrix · R Seurat objects · qs serialised perturb data · soft",
              fontweight="bold", fontsize=11)
ax.set_ylim(0, max(counts) * 1.20)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_file_categories.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_file_categories")

# ----- Fig 2: file-name → content-class breakdown -----
def classify(fname):
    fn = fname.lower()
    if "perturb_sg" in fn:
        return "Perturb sgRNA assignments\n(.qs.gz)"
    if "perturb" in fn:
        return "Perturb-seq main object\n(.qs.gz)"
    if "aav_all" in fn or "aav_ctxobj" in fn:
        return "AAV titration controls\n(.Robj.gz)"
    if "3p5p" in fn:
        return "3'/5' library comparison\n(.Robj.gz)"
    if "soft" in fn:
        return "GEO SOFT metadata"
    if "matrix" in fn:
        return "GEO series matrix"
    return "other"


cls = Counter(classify(f["name"]) for f in files)
labels = [k for k, _ in cls.most_common()]
counts = [v for _, v in cls.most_common()]
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
y = np.arange(len(labels))[::-1]
HL_KEYS = {k for k in labels if "Perturb" in k}
colors = [COL_HIGHLIGHT if l in HL_KEYS else COL_PRIMARY for l in labels]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    ax.text(c + 0.05, y[i], f"  n={c}", va="center", fontsize=10,
            fontweight="bold" if l in HL_KEYS else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("File count in GSE249416 supplementary")
ax.set_title("Supplementary-file classes — Foxg1 in-utero Perturb-seq (★) "
              "is the published analysis target",
              fontweight="bold", fontsize=11)
ax.set_xlim(0, max(counts) * 1.4)
ax.grid(axis="x", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_content_classes.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_content_classes")

print(f"\nFigures saved under {FIG_DIR}")
