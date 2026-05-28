#!/usr/bin/env python3
"""Regenerate figures for Benchmarks/martyn2025_variant_flowfish/README.md.

Reads two live data sources:
  1. IGVF Portal MeasurementSet manifest (flowfish pull-portal output)
  2. End-to-end synthetic-pipeline element-score TSV (flowfish score-elements)

And renders three publication-grade panels under figures/.
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
FIG_DIR = ROOT / "Benchmarks/martyn2025_variant_flowfish/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COL_PRIMARY = "#5C8DAA"
COL_HIGHLIGHT = "#C77F49"
COL_GREEN = "#7CA663"

# ---- Source 1: IGVF MeasurementSet manifest ----
manifests = sorted((ROOT / "Docs/FlowFISH").glob("*_martyn2025*_portal.tsv"))
if not manifests:
    sys.exit("No flowfish pull-portal TSV found — run "
              "`bash Benchmarks/martyn2025_variant_flowfish/run.sh` first.")
m = manifests[-1]
rows = list(csv.DictReader(m.open(), delimiter="\t"))
print(f"Loaded {len(rows)} MeasurementSets from {m.name}")

# Tag each row as flowfish-relevant or not
def is_flowfish(r):
    s = (r.get("assay_titles") or "" + " "
         + (r.get("preferred_assay_titles") or "")).lower()
    return ("flow" in s) or ("facs" in s) or (
        "Variant-EFFECTS" in (r.get("preferred_assay_titles") or "")
    ) or ("variant-effects" in s)


ff_subset = [r for r in rows if is_flowfish(r)]
print(f"  Flow-FISH / FACS / Variant-EFFECTS subset: {len(ff_subset)}")

# ---- Fig 1: preferred_assay_title breakdown across all IGVF MeasurementSets ----
prefs = Counter(r["preferred_assay_titles"] for r in rows
                  if r["preferred_assay_titles"])
labels = [k for k, _ in prefs.most_common(15)]
counts = [v for _, v in prefs.most_common(15)]
HL = {"CRISPR FlowFISH screen", "CRISPR FACS screen", "Variant-EFFECTS",
       "VAMP-seq (MultiSTEP)"}
colors = [COL_HIGHLIGHT if l in HL else COL_PRIMARY for l in labels]
fig, ax = plt.subplots(figsize=(9.5, 5.5), facecolor="white")
y = np.arange(len(labels))[::-1]
ax.barh(y, counts, color=colors, edgecolor="white", linewidth=0.5)
for i, (l, c) in enumerate(zip(labels, counts)):
    suffix = "  ★ Flow-FISH-family" if l in HL else ""
    ax.text(c + 1, y[i], f"  n={c}{suffix}",
            va="center", fontsize=9,
            fontweight="bold" if l in HL else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("IGVF MeasurementSets in this 500-row page")
ax.set_title("IGVF Portal MeasurementSet preferred-title mix\n"
              "Flow-FISH-family assays (★) = the Martyn 2025 reproducibility surface\n"
              "(portal-wide total: 5,780 MeasurementSets; this page = 500)",
              fontweight="bold", fontsize=11)
ax.grid(axis="x", ls=":", alpha=0.4)
ax.set_xlim(0, max(counts) * 1.30)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_assay_mix.{ext}",
                  dpi=200, facecolor="white")
plt.close(fig)
print("  ✓ fig1_assay_mix")

# ---- Source 2: end-to-end pipeline output ----
pipeline_files = sorted(
    (ROOT / "Docs/FlowFISH").glob("*martyn2025_pipeline*FullEnhancerScore*.tsv")
)
if not pipeline_files:
    print("[skip] No pipeline-output TSV; run the simulate + estimate-effects + "
          "real-space + score-elements chain.")
else:
    p = pipeline_files[-1]
    elems = list(csv.DictReader(p.open(), delimiter="\t"))
    print(f"Loaded {len(elems)} element-level scores from {p.name}")
    print(f"  columns: {list(elems[0].keys())[:10]}")

    # Fig 2: pipeline element calls (Significant vs Regulated vs others)
    sig_col = next((c for c in elems[0].keys() if "significant" in c.lower()), None)
    reg_col = next((c for c in elems[0].keys() if "regulated" in c.lower()), None)
    if sig_col and reg_col:
        n_total = len(elems)
        n_sig = sum(1 for r in elems if r[sig_col].lower() in ("true", "1"))
        n_reg = sum(1 for r in elems if r[reg_col].lower() in ("true", "1"))
        labels = ["Elements tested", "Significant\n(MWU FDR<0.05)",
                    "Regulated\n(Sig AND |1−mean|≥0.1)"]
        counts = [n_total, n_sig, n_reg]
        colors = [COL_PRIMARY, COL_GREEN, COL_HIGHLIGHT]
        fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
        bars = ax.bar(labels, counts, color=colors,
                        edgecolor="white", linewidth=0.5)
        for b, c in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, c + 0.3,
                    f"{c}",
                    ha="center", fontsize=12, fontweight="bold")
        ax.set_ylabel("Element count")
        ax.set_title("End-to-end FlowFISH pipeline on synthetic data\n"
                      "estimate-effects → real-space → score-elements",
                      fontweight="bold", fontsize=11)
        ax.set_ylim(0, max(counts) * 1.20)
        ax.grid(axis="y", ls=":", alpha=0.4)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(FIG_DIR / f"fig2_pipeline_calls.{ext}",
                          dpi=200, facecolor="white")
        plt.close(fig)
        print("  ✓ fig2_pipeline_calls")

    # Fig 3: per-element effect distribution (mean-based)
    mean_col = next((c for c in elems[0].keys()
                       if "mean" in c.lower() and "p" not in c.lower()), None)
    if mean_col:
        means = [float(r[mean_col]) for r in elems]
        fig, ax = plt.subplots(figsize=(8, 4), facecolor="white")
        # Color by significance if available
        if sig_col:
            sigs = [r[sig_col].lower() in ("true", "1") for r in elems]
            sig_means = [m for m, s in zip(means, sigs) if s]
            ns_means = [m for m, s in zip(means, sigs) if not s]
            ax.hist([ns_means, sig_means], bins=15,
                    color=[COL_PRIMARY, COL_HIGHLIGHT], stacked=True,
                    edgecolor="white",
                    label=["non-significant", "Significant (FDR<0.05)"])
            ax.legend()
        else:
            ax.hist(means, bins=15, color=COL_PRIMARY,
                    edgecolor="white")
        ax.axvline(1.0, ls="--", color="black", alpha=0.5)
        ax.set_xlabel("Per-element effect (1 = wild-type)")
        ax.set_ylabel("Element count")
        ax.set_title("Element-effect distribution from clean-room pipeline\n"
                      "Knockdowns push elements toward 0 (loss of GoI signal)",
                      fontweight="bold", fontsize=11)
        ax.grid(axis="y", ls=":", alpha=0.4)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(FIG_DIR / f"fig3_effect_distribution.{ext}",
                          dpi=200, facecolor="white")
        plt.close(fig)
        print("  ✓ fig3_effect_distribution")

print(f"\nFigures saved under {FIG_DIR}")
