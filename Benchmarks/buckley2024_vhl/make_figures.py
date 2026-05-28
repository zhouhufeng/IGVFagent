#!/usr/bin/env python3
"""Regenerate the figures embedded in Benchmarks/buckley2024_vhl/README.md."""
import csv, re, sys
from collections import Counter
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "Data/Proteomics/Sources/MaveDB/urn-mavedb-00000675-a-1.csv"
FIG_DIR = ROOT / "Benchmarks/buckley2024_vhl/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

if not CSV.is_file():
    sys.exit(f"MaveDB CSV not found at {CSV}\n"
              f"  Run first: bash Benchmarks/buckley2024_vhl/run.sh")

_HGVS_NT_RE = re.compile(
    r"^[^:]+:c\."
    r"(?:(?P<star>\*\d+)|(?P<sign>-?)(?P<pos>\d+)(?P<off>[\+\-]\d+)?)"
    r"(?P<ref>[ACGT])>(?P<alt>[ACGT])$"
)

LOF1_THR, LOF2_THR, INT_THR = -1.26, -0.4, -0.22
COL = {"LOF1": "#5C8DAA", "LOF2": "#7CA663", "Intermediate": "#D9A33E",
       "Neutral": "#A3A3A3"}

rows = []
with CSV.open() as fh:
    for r in csv.DictReader(fh):
        try:
            score = float(r.get("score", "NA"))
        except (ValueError, TypeError):
            continue
        m = _HGVS_NT_RE.match(r.get("hgvs_nt", ""))
        if not m: continue
        if m.group("star"):
            region, pos = "3UTR", int(m.group("star")[1:])
        elif m.group("sign") == "-":
            region, pos = "5UTR", -int(m.group("pos"))
        elif m.group("off"):
            region, pos = "intronic", int(m.group("pos"))
        else:
            region, pos = "CDS", int(m.group("pos"))
        rows.append({"score": score, "region": region,
                      "cds_pos": pos if region == "CDS" else None})

print(f"Loaded {len(rows):,} variants")
scores = np.array([r["score"] for r in rows])
n_lof1 = int(np.sum(scores < LOF1_THR))
n_lof2 = int(np.sum((LOF1_THR <= scores) & (scores < LOF2_THR)))
n_int  = int(np.sum((LOF2_THR <= scores) & (scores < INT_THR)))
n_neut = int(np.sum(scores >= INT_THR))

# Fig 1: score distribution
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
ax.hist(np.clip(scores, -2.5, 1.5), bins=80, color=COL["Neutral"],
         edgecolor="white", linewidth=0.4)
for thr, lbl, col in [(LOF1_THR, f"LOF1 (n={n_lof1:,})", COL["LOF1"]),
                       (LOF2_THR, f"LOF2 (n={n_lof2:,})", COL["LOF2"]),
                       (INT_THR,  f"Inter (n={n_int:,})", COL["Intermediate"])]:
    ax.axvline(thr, color=col, ls="--", lw=1.5, label=f"{lbl} threshold @ {thr}")
ax.set_xlabel("VHL functional score (D20 vs D6 depletion)")
ax.set_ylabel("Variant count")
ax.set_title(f"Buckley 2024 VHL SGE — depletion-score distribution\n"
              f"{len(rows):,} SNVs (MaveDB urn:mavedb:00000675-a-1)",
              fontweight="bold")
ax.legend(loc="upper left", fontsize=9)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig1_score_distribution.{ext}", dpi=200,
                  facecolor="white")
plt.close(fig); print("  ✓ fig1")

# Fig 2: 4-bucket bar
fig, ax = plt.subplots(figsize=(7, 4.5), facecolor="white")
buckets = ["LOF1", "LOF2", "Intermediate", "Neutral"]
vals = [n_lof1, n_lof2, n_int, n_neut]
bars = ax.bar(buckets, vals, color=[COL[b] for b in buckets],
                edgecolor="white", linewidth=0.5)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 25,
            f"{v:,}\n({100 * v / sum(vals):.1f}%)",
            ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("Variant count"); ax.set_ylim(0, max(vals) * 1.18)
ax.set_title("Buckley 2024 VHL — 4-bucket functional classification\n"
              "IGVFagent reproduction (paper thresholds)",
              fontweight="bold")
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig2_buckets.{ext}", dpi=200, facecolor="white")
plt.close(fig); print("  ✓ fig2")

# Fig 3: score along cDNA
fig, ax = plt.subplots(figsize=(11, 4.5), facecolor="white")
xs = np.array([r["cds_pos"] for r in rows if r["region"] == "CDS"])
ys = np.array([r["score"]   for r in rows if r["region"] == "CDS"])
classes = np.where(ys < LOF1_THR, "LOF1",
                   np.where(ys < LOF2_THR, "LOF2",
                            np.where(ys < INT_THR, "Intermediate", "Neutral")))
for cls in buckets:
    mask = classes == cls
    ax.scatter(xs[mask], ys[mask], s=8, color=COL[cls], alpha=0.7,
                edgecolors="none", label=f"{cls} (n={mask.sum():,})")
for name, s, e, c in [("β-domain (HIF1α binding)", 1, 300, "#FFE9C7"),
                       ("linker", 300, 465, "#F5F5F5"),
                       ("α-domain (Elongin C binding)", 465, 642, "#D4E4F0")]:
    ax.axvspan(s, e, color=c, alpha=0.5, zorder=0)
    ax.text((s + e) / 2, 0.97, name, ha="center", va="top",
            fontsize=8, fontweight="bold", color="#1F2933",
            transform=ax.get_xaxis_transform())
for thr in (LOF1_THR, LOF2_THR, INT_THR):
    ax.axhline(thr, color="#1F2933", ls=":", lw=0.7, alpha=0.4)
ax.set_xlabel("cDNA position (nt)"); ax.set_ylabel("Functional score")
ax.set_title("VHL functional score along the cDNA, colored by paper bucket",
              fontweight="bold")
ax.legend(loc="lower right", fontsize=8, ncol=2)
ax.set_xlim(0, max(xs) * 1.02)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(FIG_DIR / f"fig3_score_along_cdna.{ext}", dpi=200,
                  facecolor="white")
plt.close(fig); print("  ✓ fig3")
print(f"\nFigures saved under {FIG_DIR}")
