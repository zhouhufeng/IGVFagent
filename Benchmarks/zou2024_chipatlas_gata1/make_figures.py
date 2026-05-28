#!/usr/bin/env python3
"""Regenerate the figures embedded in Benchmarks/zou2024_chipatlas_gata1/README.md.

Pulls live data from chip-atlas.org/data/chip_antigen and saves
PNG + SVG under figures/.
"""
import json, sys
import urllib.request, urllib.parse
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "Benchmarks/zou2024_chipatlas_gata1/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def fetch_antigens(genome, ag_class, cl_class):
    url = ("https://chip-atlas.org/data/chip_antigen?"
            f"genome={genome}&agClass={urllib.parse.quote(ag_class)}"
            f"&clClass={urllib.parse.quote(cl_class)}")
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.loads(r.read())


tf_data = fetch_antigens("hg38", "TFs and others", "Blood")
items = [d for d in tf_data if d.get("id") not in ("-", "All")
          and d.get("count", 0) > 0]
items.sort(key=lambda d: -d["count"])
print(f"Loaded {len(items):,} TFs from chip-atlas.org (hg38, Blood)")

hist_data = fetch_antigens("hg38", "Histone", "Blood")
hist_items = [d for d in hist_data if d.get("id") not in ("-", "All")
                and d.get("count", 0) > 0]
hist_items.sort(key=lambda d: -d["count"])

COL_DEFAULT, COL_HIGHLIGHT = "#5C8DAA", "#C77F49"
HIGHLIGHT = {"GATA1", "TAL1", "RUNX1", "KLF1", "LMO2", "GATA2",
              "FLI1", "NFE2", "MYB"}

# Fig 1: top 25
top25 = items[:25]
fig, ax = plt.subplots(figsize=(10, 6.5), facecolor="white")
y = np.arange(len(top25))[::-1]
labels = [d["id"] for d in top25]
counts = [d["count"] for d in top25]
colors = [COL_HIGHLIGHT if l in HIGHLIGHT else COL_DEFAULT for l in labels]
bars = ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for b, c, l in zip(bars, counts, labels):
    suffix = "  ★ hematopoietic" if l in HIGHLIGHT else ""
    ax.text(c + 3, b.get_y() + b.get_height()/2, f"  n={c:,}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l in HIGHLIGHT else "normal")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("ChIP-seq experiments (hg38, Blood)")
ax.set_title("Top 25 TFs in ChIP-Atlas for hg38 Blood — "
              "IGVFagent chipatlas list-antigens\n"
              "Canonical hematopoietic TFs (★) cluster at the top, "
              "as Zou 2024 predicts",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.set_xlim(0, max(counts) * 1.20)
fig.tight_layout()
for ext in ("png","svg"):
    fig.savefig(FIG_DIR / f"fig1_top25_tfs_blood.{ext}", dpi=200,
                  facecolor="white")
plt.close(fig); print("  ✓ fig1")

# Fig 2: hematopoietic TF panel
gata = [d for d in items if d["id"].startswith("GATA") or d["id"] in HIGHLIGHT]
gata.sort(key=lambda d: -d["count"])
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
ts = [d["id"] for d in gata]; cs = [d["count"] for d in gata]
colors = [COL_HIGHLIGHT if t == "GATA1" else COL_DEFAULT for t in ts]
bars = ax.bar(ts, cs, color=colors, edgecolor="white", linewidth=0.5)
for b, c, t in zip(bars, cs, ts):
    weight = "bold" if t == "GATA1" else "normal"
    ax.text(b.get_x()+b.get_width()/2, c+2, f"{c:,}",
            ha="center", fontsize=10, fontweight=weight)
ax.set_ylabel("ChIP-seq experiments (hg38, Blood)")
ax.set_title("Canonical hematopoietic TFs in ChIP-Atlas",
              fontweight="bold")
ax.set_ylim(0, max(cs)*1.15)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
for t in ax.get_xticklabels():
    t.set_rotation(35); t.set_ha("right")
fig.tight_layout()
for ext in ("png","svg"):
    fig.savefig(FIG_DIR / f"fig2_hematopoietic_tfs.{ext}", dpi=200,
                  facecolor="white")
plt.close(fig); print("  ✓ fig2")

# Fig 3: TFs vs histones
fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
top_tfs = items[:15]; top_hists = hist_items[:15]
x_tf = np.arange(len(top_tfs))
ax.bar(x_tf - 0.2, [d["count"] for d in top_tfs], 0.4,
        color=COL_DEFAULT, label="TFs")
ax.bar(x_tf + 0.2, [d["count"] for d in top_hists], 0.4,
        color="#7CA663", label="Histone marks")
ax.set_xticks(x_tf)
ax.set_xticklabels([f"{t['id']}\n{h['id']}"
                      for t,h in zip(top_tfs, top_hists)],
                     fontsize=8)
ax.set_ylabel("ChIP-seq experiments (hg38, Blood)")
ax.set_title("Top 15 antigens — TFs vs Histones (hg38 Blood)",
              fontweight="bold")
ax.legend()
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png","svg"):
    fig.savefig(FIG_DIR / f"fig3_tfs_vs_histones.{ext}", dpi=200,
                  facecolor="white")
plt.close(fig); print("  ✓ fig3")
print(f"\nFigures saved under {FIG_DIR}")
