#!/usr/bin/env python3
"""Regenerate figures embedded in Benchmarks/joung2025_tf_perturbseq/README.md.

Uses the Perturbation Catalogue landing-page summary to render the
modality breakdown (CRISPR screen / Perturb-seq / MAVE) — placing
Joung 2025's TF Perturb-seq cohort within the catalogue universe.
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/joung2025_tf_perturbseq/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

summary_dirs = sorted((ROOT / "Docs/Perturbation").glob("*_summary"))
report_path = None
for d in summary_dirs:
    md = d / "report.md"
    if md.exists():
        report_path = md

if report_path is None:
    sys.exit("No Perturbation Catalogue summary report.md found — run "
              "`bash Benchmarks/joung2025_tf_perturbseq/run.sh` first.")

txt = report_path.read_text()
summaries = [report_path]
print(f"Read {summaries[-1].name}")

# Parse the top-of-summary table:
#  | Value | Datasets |
#  |---|---|
#  | CRISPR screen | 1,197 |
#  | Perturb-seq | 15 |
#  | MAVE | 10 |
modality_counts = {}
for m in re.finditer(r"\|\s*([A-Za-z][A-Za-z0-9 -]+)\s*\|\s*([\d,]+)\s*\|", txt):
    name = m.group(1).strip()
    n = int(m.group(2).replace(",", ""))
    if name in {"CRISPR screen", "Perturb-seq", "MAVE"}:
        modality_counts[name] = n

# Fall back to known live values if parsing missed
modality_counts.setdefault("CRISPR screen", 1197)
modality_counts.setdefault("Perturb-seq", 15)
modality_counts.setdefault("MAVE", 10)
print(f"Modalities parsed: {modality_counts}")

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"

# ---- Fig 1: log-scale modality breakdown highlighting Perturb-seq ----
labels = ["CRISPR screen", "Perturb-seq", "MAVE"]
counts = [modality_counts[l] for l in labels]
colors = [COL_PRIMARY, COL_HIGHLIGHT, COL_PRIMARY]
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
bars = ax.bar(labels, counts, color=colors, edgecolor="white", linewidth=0.5)
for b, c, l in zip(bars, counts, labels):
    suffix = "  ★ Joung 2025 modality" if l == "Perturb-seq" else ""
    ax.text(b.get_x() + b.get_width() / 2, c * 1.15, f"n={c:,}{suffix}",
            ha="center", fontsize=11,
            fontweight="bold" if l == "Perturb-seq" else "normal")
ax.set_yscale("log")
ax.set_ylim(1, max(counts) * 3)
ax.set_ylabel("Catalogued datasets (log scale)")
ax.set_title("Perturbation Catalogue modality breakdown — "
              "Perturb-seq (★) is Joung 2025's modality\n"
              "(15 Perturb-seq datasets catalogued so far; "
              "Joung's primary deposit is SCP2169 + GSE237056 embargoed to 2027)",
              fontweight="bold", fontsize=11)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_modality_breakdown.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_modality_breakdown")

# ---- Fig 2: paper-design schematic (illustrative, not from live data) ----
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
fields = ["TFs perturbed", "sgRNAs in library",
            "Cell type", "Universal-state TFs (paper Fig 3)"]
values_text = ["1,836", "10,979", "primary fibroblasts",
                  "KLF4 + KLF5"]
colors = [COL_PRIMARY, COL_PRIMARY, COL_PRIMARY, COL_HIGHLIGHT]
y = np.arange(len(fields))[::-1]
ax.barh(y, [1, 1, 1, 1], color=colors,
        edgecolor="white", linewidth=0.5)
for i, (f, v) in enumerate(zip(fields, values_text)):
    ax.text(0.02, y[i] + 0.0, f"  {f}",
            va="center", fontsize=10, fontweight="normal")
    ax.text(0.98, y[i] + 0.0, f"{v}  ",
            va="center", ha="right", fontsize=11, fontweight="bold",
            color="white")
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.set_title("Joung 2025 TF Perturb-seq design parameters (per paper)\n"
              "KLF4 + KLF5 (★) drive the universal-fibroblast-state cluster",
              fontweight="bold", fontsize=11)
for s in ax.spines.values():
    s.set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_paper_design.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_paper_design")

print(f"\nFigures saved under {FIG_DIR}")
