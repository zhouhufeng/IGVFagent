#!/usr/bin/env python3
"""Score + plot the deMULTIplex2 / Stoeckius cell-hashing reproduction.

Reads the classification written by `igvfagent multiseq demultiplex` on the
Stoeckius 2018 8-donor PBMC HTO matrix (bundled with deMULTIplex2) and writes:

  * concordance_metrics.json  — singlet/multiplet/negative rates, HTO balance
  * fig1_classification.png/.svg — droplet-type breakdown + per-HTO singlet bar
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def latest_run() -> Path:
    dirs = sorted((ROOT / "Docs/MultiSeq").glob("2*_demultiplex2_stoeckius"),
                  reverse=True)
    if not dirs:
        sys.exit("No multiseq run — run Benchmarks/demultiplex2_stoeckius/run.sh first.")
    return dirs[0]


RUN = latest_run()
summary = json.loads((RUN / "summary.json").read_text())
cls = pd.read_csv(RUN / "classifications.csv", index_col=0)

by_type = cls["droplet_type"].value_counts().to_dict()
n = int(cls.shape[0])
singlets = cls[cls["droplet_type"] == "singlet"]
hto_counts = singlets["barcode_assign"].value_counts().sort_index()
n_groups = int(hto_counts.shape[0])
# balance = min/max singlet count across HTOs (1.0 = perfectly even 8-donor pool)
balance = float(hto_counts.min() / hto_counts.max()) if n_groups else 0.0

metrics = {
    "n_cells": n,
    "n_tags": int(summary.get("n_tags", 0)),
    "n_singlet_hto_groups": n_groups,
    "singlet_rate": round(by_type.get("singlet", 0) / n, 4),
    "multiplet_rate": round(by_type.get("multiplet", 0) / n, 4),
    "negative_rate": round(by_type.get("negative", 0) / n, 4),
    "hto_balance_min_over_max": round(balance, 3),
    "failed_tags": summary.get("failed_tags", []),
}
(FIG / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
(RUN / "concordance_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics, indent=2))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor="white")
types = ["singlet", "multiplet", "negative"]
vals = [by_type.get(t, 0) for t in types]
cols = ["#4C78A8", "#E45756", "#BAB0AC"]
ax1.bar(types, vals, color=cols, edgecolor="white")
for i, v in enumerate(vals):
    ax1.text(i, v + n * 0.01, f"{v:,}\n{100*v/n:.1f}%", ha="center",
             fontsize=10, fontweight="bold")
ax1.set_ylabel("cells"); ax1.set_ylim(0, max(vals) * 1.18)
ax1.set_title("Droplet-type classification", fontweight="bold")
for s in ("top", "right"):
    ax1.spines[s].set_visible(False)

ax2.bar(hto_counts.index, hto_counts.values, color="#54A24B", edgecolor="white")
ax2.set_ylabel("singlet cells"); ax2.set_xlabel("HTO")
ax2.set_title(f"Per-HTO singlets — {n_groups}/8 groups, "
              f"balance {balance:.2f}", fontweight="bold")
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)
fig.suptitle("deMULTIplex2 reproduction — Stoeckius 2018 8-donor PBMC "
             f"({n:,} cells)", fontweight="bold")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG / f"fig1_classification.{ext}", dpi=170, facecolor="white")
plt.close(fig)
print(f"figures under {FIG.relative_to(ROOT)}")
