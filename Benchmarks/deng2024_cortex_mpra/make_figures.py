#!/usr/bin/env python3
"""Regenerate figures embedded in Benchmarks/deng2024_cortex_mpra/README.md.

Two source layers:
  1. IGVF Portal MPRA-class manifest (mpra portal-manifest output).
  2. Paper-headline statistics from Deng 2024 (Methods + Fig 1)
     used to render an illustrative element-fate panel.
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
FIG_DIR = ROOT / "Benchmarks/deng2024_cortex_mpra/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"
COL_GREEN = "#7CA663"
COL_RED = "#C9605B"

# ---- Source 1: IGVF Portal MPRA-class manifest ----
manifests = sorted(
    (ROOT / "Data/Manifests/MPRA").glob("*_deng2024_cortex_mpra_portal_many_manifest.csv")
)
if not manifests:
    print("[warn] No MPRA portal manifest; only paper-headline figures will render.")
    rows = []
else:
    rows = list(csv.DictReader(manifests[-1].open()))
    print(f"Loaded {len(rows)} IGVF MPRA-class entries from {manifests[-1].name}")

# ---- Fig 1: Paper headline — element fate distribution (Deng Fig 1) ----
N_TOTAL = 102767
N_ACTIVE = 46802
N_INACTIVE = N_TOTAL - N_ACTIVE
N_VARIANTS = 164  # variants that significantly alter activity
N_SCZ = 100  # approximation per paper Discussion
N_BP = 64    # approximation per paper Discussion

# Element-level activity stack
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
bars = ax.barh(
    ["Elements tested"], [N_TOTAL],
    color="#dddddd", edgecolor="white"
)
# Active fraction
ax.barh(
    ["Elements tested"], [N_ACTIVE],
    color=COL_HIGHLIGHT, edgecolor="white",
    label=f"Active enhancers ({N_ACTIVE:,} = {N_ACTIVE/N_TOTAL*100:.1f} %)"
)
ax.text(N_TOTAL * 0.50, 0,
        f"  {N_ACTIVE:,} active  /  {N_INACTIVE:,} inactive",
        va="center", ha="center", fontsize=11, fontweight="bold",
        color="white")
ax.set_xlabel("Element count")
ax.set_yticks([])
ax.set_title("Deng 2024 cortex lentiMPRA — published element-fate breakdown\n"
              f"({N_TOTAL:,} elements tested → {N_ACTIVE:,} active = "
              f"{N_ACTIVE/N_TOTAL*100:.1f} % active fraction)",
              fontweight="bold", fontsize=11)
ax.set_xlim(0, N_TOTAL * 1.05)
ax.grid(axis="x", ls=":", alpha=0.4)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.legend(loc="lower right", frameon=False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_paper_headline.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_paper_headline")

# ---- Fig 2: variant-impact panel (paper Fig 5) ----
labels = ["Variants altering\nMPRA activity",
            "SCZ-GWAS\nimplicated",
            "BP-GWAS\nimplicated"]
counts = [N_VARIANTS, N_SCZ, N_BP]
colors = [COL_HIGHLIGHT, COL_RED, COL_GREEN]
fig, ax = plt.subplots(figsize=(8, 4), facecolor="white")
bars = ax.bar(labels, counts, color=colors,
                edgecolor="white", linewidth=0.5)
for b, c in zip(bars, counts):
    ax.text(b.get_x() + b.get_width() / 2, c + 3,
            f"{c}",
            ha="center", fontsize=12, fontweight="bold")
ax.set_ylabel("Variant count")
ax.set_title("Deng 2024 — psychiatric-disorder variant calls (paper Fig 5)\n"
              "(SCZ + BP counts are approximations from Discussion)",
              fontweight="bold", fontsize=11)
ax.set_ylim(0, max(counts) * 1.25)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_variant_impact.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig2_variant_impact")

# ---- Fig 3 (NEW): Synapse NeuREs file inventory (live, via new synapse skill) ----
synapse_manifests = sorted(
    list((ROOT / "Data/Manifests/Synapse").glob("*deng2024*mpra*children.csv"))
    + list((ROOT / "Data/Manifests/Synapse").glob("*deng2024*capstone*children.csv"))
)
if synapse_manifests:
    syn_rows = list(csv.DictReader(synapse_manifests[-1].open()))
    print(f"Loaded {len(syn_rows)} Synapse MPRA_CapstoneII children")

    def bucket(name):
        n = name.lower()
        if "primary" in n and "dna" in n: return "Primary cortex — DNA"
        if "primary" in n and "rna" in n: return "Primary cortex — RNA"
        if "organoid" in n and "dna" in n: return "Cerebral organoid — DNA"
        if "organoid" in n and "rna" in n: return "Cerebral organoid — RNA"
        if (n.startswith(("a", "s"))
              and (n.endswith(".fq.gz") or n.endswith(".fastq.gz"))):
            return "Per-donor bulk RNA"
        return "Other"

    counter = Counter(bucket(r["name"]) for r in syn_rows)
    labels = [k for k, _ in counter.most_common()]
    counts = [v for _, v in counter.most_common()]
    fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="white")
    y = np.arange(len(labels))[::-1]
    color_map = {
        "Primary cortex — DNA": COL_GREEN,
        "Primary cortex — RNA": COL_HIGHLIGHT,
        "Cerebral organoid — DNA": "#7FB069",
        "Cerebral organoid — RNA": "#D08B5B",
        "Per-donor bulk RNA": COL_PRIMARY,
        "Other": "#aaaaaa",
    }
    colors = [color_map.get(l, COL_PRIMARY) for l in labels]
    ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
    for i, (l, c) in enumerate(zip(labels, counts)):
        ax.text(c + 0.3, y[i], f"  n={c}",
                va="center", fontsize=10, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("File count")
    ax.set_title(f"Deng 2024 lentiMPRA Synapse deposit "
                  f"({len(syn_rows)} files in MPRA_CapstoneII, syn51090452)\n"
                  "Paired primary cortex + cerebral organoid DNA/RNA design",
                  fontweight="bold", fontsize=11)
    ax.set_xlim(0, max(counts) * 1.20)
    ax.grid(axis="x", ls=":", alpha=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG_DIR / f"fig4_synapse_inventory.{ext}",
                      dpi=200, facecolor="white")
    plt.close(fig)
    print("  ✓ fig4_synapse_inventory")

# ---- Fig 3: IGVF Portal MPRA-class entries (live) ----
if rows:
    fig, ax = plt.subplots(figsize=(9, 4), facecolor="white")
    descs = [(r["accession"], (r.get("description") or "")[:60] or "<no desc>")
              for r in rows]
    descs = descs[:15]
    y = np.arange(len(descs))[::-1]
    ax.barh(y, [1] * len(descs), color=COL_PRIMARY,
            edgecolor="white", linewidth=0.5)
    for i, (acc, desc) in enumerate(descs):
        ax.text(0.02, y[i], f"  {acc}",
                va="center", fontsize=8, fontweight="bold",
                color="white")
        ax.text(0.20, y[i], f"  {desc}",
                va="center", fontsize=8, color="white")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax.set_title(f"IGVF Portal MPRA / construct-library-set entries "
                  f"({len(rows)} reachable via `mpra portal-manifest`)",
                  fontweight="bold", fontsize=11)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG_DIR / f"fig3_portal_mpra_entries.{ext}",
                      dpi=200, facecolor="white")
    plt.close(fig)
    print("  ✓ fig3_portal_mpra_entries")

print(f"\nFigures saved under {FIG_DIR}")
