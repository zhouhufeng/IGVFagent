#!/usr/bin/env python3
"""Figures for the scE2G crowdsourced-feature benchmark.

Reads the merged run under Docs/scE2G/*_k562_crowdsourced_feature_benchmark
(benchmark_summary.tsv) and the per-batch pr_curves.tsv files, and writes
PNG + SVG under Benchmarks/sce2g_crowdsourced_features_k562/figures/.

  fig1  ranked AUPRC, every predictor above the distance baseline plus the
        baseline and the four reference scores, with bootstrap intervals
  fig2  precision-recall curves: released scE2G vs the best single features
  fig3  best single feature per contributed table (the "family" view)
"""
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DOCS = os.path.join(ROOT, "Docs", "scE2G")

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def latest(pattern):
    hits = sorted(glob.glob(os.path.join(DOCS, pattern)))
    if not hits:
        raise SystemExit(f"no run matches {pattern}; run run.sh first")
    return hits[-1]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def style(ax, title, xlabel="", ylabel=""):
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def save(fig, name):
    fig.patch.set_facecolor(SURFACE)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", dpi=130)
    plt.close(fig)
    print("wrote", name)


merged = latest("*_k562_crowdsourced_feature_benchmark")
rows = read_tsv(os.path.join(merged, "benchmark_summary.tsv"))
for r in rows:
    for k in ("auprc", "auprc_ci95_low", "auprc_ci95_high", "precision_at_70_recall", "auprc_negated"):
        r[k] = float(r[k])
    r["n_pairs"], r["n_positive"] = int(r["n_pairs"]), int(r["n_positive"])
baseline_prec = rows[0]["n_positive"] / rows[0]["n_pairs"]
dist = next(r["auprc"] for r in rows if r["pred_id"] == "baseline:E2G_Distance")

# --- fig1: the references, the distance baseline, and the 12 strongest single features
feats = [r for r in rows if ":" in r["pred_id"] and not r["pred_id"].startswith("baseline:")]
feats = sorted(feats, key=lambda r: -r["auprc"])[:12]
keep = [r for r in rows if ":" not in r["pred_id"] or r["pred_id"].startswith("baseline:")] + feats
keep = sorted(keep, key=lambda r: r["auprc"])
fig, ax = plt.subplots(figsize=(9, 0.32 * len(keep) + 1.5))
names = [r["pred_id"] for r in keep]
vals = [r["auprc"] for r in keep]
colors = [ORANGE if n.startswith("baseline:") else ("#4a3aa7" if ":" not in n else BLUE) for n in names]
ax.barh(names, vals, color=colors, height=0.72)
ax.errorbar(vals, range(len(vals)),
            xerr=[[max(0, r["auprc"] - r["auprc_ci95_low"]) for r in keep],
                  [max(0, r["auprc_ci95_high"] - r["auprc"]) for r in keep]],
            fmt="none", ecolor=INK2, elinewidth=0.8, capsize=2)
ax.axvline(baseline_prec, color=AXIS, ls="--", lw=0.8)
for i, v in enumerate(vals):
    ax.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=8, color=INK)
style(ax, "AUPRC on 10,356 K562 CRISPR pairs: 12 strongest single crowdsourced features (blue) vs the\n"
          "released scE2G / rE2G model scores (violet) and distance (orange); dashed = random", "AUPRC")
save(fig, "fig1_auprc_ranked")

# --- fig2: PR curves, released scE2G vs best single features
curves = {}
for run in ("*_k562_sce2g_reference", "*_k562_features_A", "*_k562_features_B", "*_k562_features_D_epcot"):
    for r in read_tsv(os.path.join(latest(run), "pr_curves.tsv")):
        curves.setdefault(r["pred_id"], []).append((float(r["recall"]), float(r["precision"])))
want = ["scE2G", "ABC", "Pinloop:Pinloop", "Signac:Signac_Score", "SCENT:SCENT_beta", "EPCOT:EP300 signal at E",
        "ChromHMM:E2G_Distance"]
auprc = {r["pred_id"]: r["auprc"] for r in rows}
fig, ax = plt.subplots(figsize=(7, 5.2))
for k, pid in enumerate(want):
    c = curves.get(pid)
    if not c:
        continue
    label = "distance baseline" if pid.endswith("E2G_Distance") else pid
    a = auprc.get("baseline:E2G_Distance" if pid.endswith("E2G_Distance") else pid, float("nan"))
    ax.plot([p[0] for p in c], [p[1] for p in c], lw=1.6, color=SERIES[k % len(SERIES)],
            label=f"{label} (AUPRC {a:.3f})")
ax.axhline(baseline_prec, color=AXIS, ls="--", lw=0.8)
ax.axvline(0.7, color=INK2, ls=":", lw=0.8)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(fontsize=8, frameon=False, loc="upper right")
style(ax, "Precision-recall on K562 CRISPR pairs: released scE2G vs the strongest single features",
      "recall", "precision")
save(fig, "fig2_pr_curves")

# --- fig3: best feature per contributed table
best = {}
for r in rows:
    if ":" not in r["pred_id"] or r["pred_id"].startswith("baseline:"):
        continue
    table = r["pred_id"].split(":")[0]
    if table not in best or r["auprc"] > best[table]["auprc"]:
        best[table] = r
items = sorted(best.values(), key=lambda r: r["auprc"])
fig, ax = plt.subplots(figsize=(9, 0.36 * len(items) + 1.5))
ax.barh([f"{r['pred_id'].split(':')[0]}  ({r['pred_id'].split(':', 1)[1]})" for r in items],
        [r["auprc"] for r in items], color=BLUE, height=0.72)
ax.axvline(dist, color=ORANGE, ls="--", lw=1.0)
ax.axvline(baseline_prec, color=AXIS, ls="--", lw=0.8)
for i, r in enumerate(items):
    ax.text(r["auprc"] + 0.004, i, f"{r['auprc']:.3f}", va="center", fontsize=8, color=INK)
style(ax, "Best single feature per contributed table (orange dashed = distance baseline; grey = random)", "AUPRC")
save(fig, "fig3_best_feature_per_table")
