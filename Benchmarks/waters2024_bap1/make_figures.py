#!/usr/bin/env python3
"""Regenerate the figures embedded in Benchmarks/waters2024_bap1/README.md.

Run after a fresh `bash Benchmarks/waters2024_bap1/run.sh` has put the
MaveDB CSV at the expected path. Saves PNG + SVG under figures/.
"""
import csv, re, sys
from collections import Counter
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "Data/Proteomics/Sources/MaveDB/urn-mavedb-00000662-0-1.csv"
FIG_DIR = ROOT / "Benchmarks/waters2024_bap1/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

if not CSV.is_file():
    sys.exit(f"MaveDB CSV not found at {CSV}\n"
              f"  Run first: bash Benchmarks/waters2024_bap1/run.sh\n")

_HGVS_NT_RE = re.compile(
    r"^[^:]+:c\."
    r"(?:(?P<star>\*\d+)|(?P<sign>-?)(?P<pos>\d+)(?P<off>[\+\-]\d+)?)"
    r"(?P<ref>[ACGT])>(?P<alt>[ACGT])$"
)

COL_LOF, COL_GOF, COL_NEUTRAL = "#5C8DAA", "#C77F49", "#A3A3A3"
COL_PAPER, COL_OURS = "#1F2933", "#5C8DAA"

rows = []
with CSV.open() as fh:
    for r in csv.DictReader(fh):
        try:
            z21 = float(r.get("processed_Z_D4_D21", "NA"))
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
        cls = "LOF" if z21 < -2.5 else "GOF" if z21 > 2.5 else "Neutral"
        rows.append({"z": z21, "region": region,
                      "cds_pos": pos if region == "CDS" else None, "class": cls})

print(f"Loaded {len(rows):,} variants")
n_lof = sum(r["class"]=="LOF" for r in rows)
n_gof = sum(r["class"]=="GOF" for r in rows)

# Fig 1: z-score distribution
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
ax.hist(np.clip([r["z"] for r in rows], -6, 6), bins=80,
         color=COL_NEUTRAL, edgecolor="white", linewidth=0.4)
ax.axvline(-2.5, color=COL_LOF, ls="--", lw=1.5,
            label=f"LOF threshold (z < −2.5): n={n_lof:,}")
ax.axvline(+2.5, color=COL_GOF, ls="--", lw=1.5,
            label=f"GOF threshold (z > +2.5): n={n_gof:,}")
ax.set_xlabel("Z-score (processed_Z_D4_D21)"); ax.set_ylabel("Variant count")
ax.set_title("Waters 2024 BAP1 SGE — depletion z-score distribution\n"
              "18,108 SNVs (MaveDB urn:mavedb:00000662-0-1)", fontweight="bold")
ax.legend(); ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png","svg"): fig.savefig(FIG_DIR/f"fig1_zscore_distribution.{ext}", dpi=200, facecolor="white")
plt.close(fig); print("  ✓ fig1")

# Fig 2: concordance
fig, ax = plt.subplots(figsize=(7, 4.5), facecolor="white")
ours, paper = [n_lof, n_gof], [5665, 531]
x = np.arange(2); w = 0.34
ax.bar(x - w/2, paper, w, color=COL_PAPER, label="Waters 2024 (paper)", edgecolor="white")
ax.bar(x + w/2, ours, w, color=COL_OURS, label="IGVFagent reproduction", edgecolor="white")
for i,(p,o) in enumerate(zip(paper, ours)):
    ax.text(x[i]-w/2, p+80, f"{p:,}", ha="center", color=COL_PAPER)
    d = o-p; pct = 100*d/p
    ax.text(x[i]+w/2, o+80, f"{o:,}\nΔ={d:+,} ({pct:+.1f}%)",
            ha="center", color=COL_OURS, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(["LOF","GOF"])
ax.set_ylabel("Variant count"); ax.set_ylim(0, max(ours+paper)*1.18)
ax.set_title("BAP1 SGE concordance: IGVFagent vs Waters 2024", fontweight="bold")
ax.legend(); ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png","svg"): fig.savefig(FIG_DIR/f"fig2_concordance.{ext}", dpi=200, facecolor="white")
plt.close(fig); print("  ✓ fig2")

# Fig 3: region distribution
fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="white")
regions = ["5UTR","CDS","intronic","3UTR"]
counts = {r: Counter() for r in regions}
for r in rows: counts[r["region"]][r["class"]] += 1
bottoms = np.zeros(len(regions))
for cls,col in [("LOF",COL_LOF),("GOF",COL_GOF),("Neutral",COL_NEUTRAL)]:
    vals = [counts[r][cls] for r in regions]
    ax.bar(regions, vals, bottom=bottoms, color=col, edgecolor="white", label=cls)
    bottoms += np.array(vals)
for i,r in enumerate(regions):
    ax.text(i, sum(counts[r].values())+80, f"n={sum(counts[r].values()):,}",
            ha="center", fontweight="bold")
ax.set_ylabel("Variant count")
ax.set_title("Variant distribution by transcript region\n"
              "IGVFagent SGE skill (parse_hgvsc_full) on Waters 2024 BAP1", fontweight="bold")
ax.legend(); ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png","svg"): fig.savefig(FIG_DIR/f"fig3_region_distribution.{ext}", dpi=200, facecolor="white")
plt.close(fig); print("  ✓ fig3")

# Fig 4: LOF along cDNA with paper-domain overlay
fig, ax = plt.subplots(figsize=(11, 4.5), facecolor="white")
lof_pos = np.array([r["cds_pos"] for r in rows if r["region"]=="CDS" and r["class"]=="LOF"])
all_pos = np.array([r["cds_pos"] for r in rows if r["region"]=="CDS"])
bins = np.arange(0, 2300, 50)
lof_h, _ = np.histogram(lof_pos, bins=bins)
all_h, _ = np.histogram(all_pos, bins=bins)
lof_frac = np.divide(lof_h, all_h, out=np.zeros_like(lof_h, dtype=float), where=all_h>0)
centers = (bins[:-1] + bins[1:]) / 2
ax.bar(centers, lof_frac, width=45, color=COL_LOF, edgecolor="white", alpha=0.85)
domains = [("UCH catalytic domain (paper: most depleted)", 1, 720, "#FFE9C7"),
            ("BBD / HCF-1 binding", 1095, 1155, "#E7F4D8"),
            ("ULD", 1785, 2163, "#D4E4F0")]
for name, s, e, c in domains:
    ax.axvspan(s, e, color=c, alpha=0.5, zorder=0)
    ax.text((s+e)/2, 0.94, name, ha="center", va="top", fontsize=8,
            fontweight="bold", color="#1F2933",
            transform=ax.get_xaxis_transform())
ax.set_xlabel("cDNA position (nt)"); ax.set_ylabel("Fraction LOF (z < −2.5)")
ax.set_title("LOF fraction along BAP1 cDNA — paper claims UCH domain "
              "(exons 1–9, c.~1–720) is most depleted", fontweight="bold")
ax.set_xlim(0, 2300); ax.set_ylim(0, 1.0)
ax.grid(axis="y", ls=":", alpha=0.4)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png","svg"): fig.savefig(FIG_DIR/f"fig4_lof_along_cdna.{ext}", dpi=200, facecolor="white")
plt.close(fig); print("  ✓ fig4")
print(f"\nFigures saved under {FIG_DIR}")
