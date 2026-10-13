#!/usr/bin/env python3
"""Overview figure: IGVFagent reproduces three published papers.

Laid out like Paper2Agent's Fig. 3:
  a  paper -> IGVFagent (ported methods, deposited data, authors' code) -> agent
  b  the chain of `igvfagent bench` steps the igvf-replicate-paper skill runs,
     and the report it returns
  c  results per paper: a redrawn panel, every quoted number vs the paper,
     and analysis coverage from `igvfagent bench score`

All numbers are the ones recorded in each benchmark's README
(Benchmarks/<paper-id>/README.md). Redrawn panels are read from
Benchmarks/<paper-id>/figures/.

    python3 Benchmarks/make_reproduction_figure.py
writes Benchmarks/figures/reproduction_overview.{png,pdf,svg}.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "figures" / "reproduction_overview"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

INK = "#1F2933"
MUTED = "#616E7C"
TOOL = "#FBE3D3"      # ported methods
DATA = "#DDF1DC"      # deposited data
CODE = "#DCEBF7"      # authors' code
USER = "#F6DDF0"      # user bubble, as in Paper2Agent
AGENT = "#E6EEF8"
PARTIAL = "#F4A261"
BLOCKED = "#BDBDBD"
MONO = "DejaVu Sans Mono"

# (paper id, header, colour, panel file, panel caption)
PAPERS = [
    ("liu2025_kidney_multiome", "Liu et al. 2025, Science", "#2E7EBE",
     "figures/fig3d_basa_vs_ase.png",
     "Fig. 3D: chromatin effect vs expression effect\n(R = 0.520, paper 0.51)"),
    ("gschwind2026_encyclopedia_enhancer", "Gschwind et al. 2026, Nature",
     "#D1495B", None,
     "Fig. 2b,c: CRISPR benchmark AUPRC\n(ours from Supp. Table 3 vs paper)"),
    ("pan2026_esig_interaction", "Pan et al. 2026, Nat Methods", "#2A9D8F",
     "figures/fig1c_roc_disease.png",
     "Fig. 1c: eSIG-Net finds disease mutations\n(ROC-AUC 0.908, paper 0.91)"),
]

# Every number each README compares with the paper: (label, paper, ours).
NUMBERS = {
    "liu2025_kidney_multiome": [
        ("GWAS variants, P<5e-8", 115726, 115725),
        ("non-MHC variants", 99595, 102520),
        ("independent loci", 1026, 1020),
        ("ASE genes, tubule", 8573, 8573),
        ("ASE genes, glomerulus", 8725, 8725),
        ("ASE genes, union", 10398, 10398),
        ("bASA peaks", 19083, 19083),
        ("bASA SNPs", 523700, 523700),
        ("bASA vs ASE R", 0.51, 0.520),
        ("snASA SNPs", 25766, 25766),
        ("snASA in open chromatin (%)", 93.2, 93.20),
        ("PT snASA vs bASA rho", 0.90, 0.895),
        ("Open4Gene links", 142452, 142452),
        ("Open4Gene peaks", 125699, 125699),
        ("median peaks/gene", 8, 4),
        ("GWAS variants in peaks", 7137, 7344),
        ("coding variants", 1363, 1363),
        ("coding genes", 782, 782),
        ("convergent genes", 601, 601),
        ("validated (%)", 93.8, 93.7),
    ],
    "gschwind2026_encyclopedia_enhancer": [
        ("biosamples", 1458, 1458),
        ("cell types / tissues", 369, 369),
        ("elements per biosample", 40210, 40209.5),
        ("links per biosample", 63221, 63221.0),
        ("predicted links, total", 92176227, 92176227),
        ("genes per enhancer", 1, 1),
        ("enhancers per gene", 3, 3),
        ("AUPRC, rE2G", 0.66, 0.6622),
        ("AUPRC, ABC", 0.56, 0.5649),
        ("AUPRC, rE2G Extended", 0.74, 0.7373),
        ("AUPRC, ABC (Extended set)", 0.61, 0.6128),
    ],
    "pan2026_esig_interaction": [
        ("accuracy, disease", 0.85, 0.848),
        ("ROC-AUC, disease", 0.91, 0.908),
        ("AP, disease", 0.86, 0.856),
        ("accuracy, population", 0.90, 0.907),
        ("ROC-AUC, population", 0.87, 0.866),
        ("AP, population", 0.71, 0.703),
        ("ROC-AUC vs structure", 0.91, 0.908),
        ("AP vs structure", 0.86, 0.856),
        ("FoldDock mutant acc.", 0.625, 0.625),
        ("FoldDock consistency (%)", 94.4, 94.4),
        ("FoldDock 17.6% (%)", 17.6, 17.58),
        ("FoldDock WT acc.", 0.922, 0.949),
    ],
}

# `igvfagent bench score` coverage: reproduced, partial / not reproduced,
# blocked or pending.
COVERAGE = {
    "liu2025_kidney_multiome": (11, 3, 2),
    "gschwind2026_encyclopedia_enhancer": (5, 0, 0),
    "pan2026_esig_interaction": (11, 1, 4),
}

# Fig. 2b,c AUPRC (Supplementary Table 3): (model, paper, ours).
GSCHWIND_AUPRC = [
    ("ENCODE-rE2G", 0.66, 0.6622),
    ("ABC", 0.56, 0.5649),
    ("rE2G\nExtended", 0.74, 0.7373),
    ("ABC\n(Ext. set)", 0.61, 0.6128),
]


# ---------------------------------------------------------------- helpers
def box(fig, x, y, w, h, fc, ec="none", lw=0.6, ls="-", r=0.008, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       transform=fig.transFigure, fc=fc, ec=ec, lw=lw, ls=ls,
                       zorder=z)
    fig.add_artist(p)
    return p


def arrow(fig, x0, y0, x1, y1, color=MUTED):
    fig.add_artist(FancyArrowPatch(
        (x0, y0), (x1, y1), transform=fig.transFigure, arrowstyle="-|>",
        mutation_scale=8, lw=1.0, color=color, zorder=3))


def text(fig, x, y, s, **kw):
    kw.setdefault("color", INK)
    kw.setdefault("va", "center")
    return fig.text(x, y, s, **kw)


def panel_label(fig, x, y, s):
    text(fig, x, y, s, fontsize=10, weight="bold", va="top")


def logo(fig, x, y, size):
    path = ROOT / "Docs" / "Figures" / "logo-mark.png"
    if path.is_file():
        ax = fig.add_axes([x, y, size, size * fig.get_figwidth()
                           / fig.get_figheight()], zorder=10)
        ax.imshow(mpimg.imread(path))
        ax.axis("off")


# ---------------------------------------------------------------- panel a
def panel_a(fig):
    panel_label(fig, 0.01, 0.995, "a")
    top, bot = 0.975, 0.745

    # Papers, stacked like a pile of PDFs.
    for i, (_, _, color, *_rest) in enumerate(reversed(PAPERS)):
        k = 2 - i                      # 2 = furthest back
        x, y = 0.022 + 0.010 * k, 0.800 + 0.012 * k
        box(fig, x, y, 0.130, 0.130, "white", ec="#9AA5B1", r=0.003, z=2 + i)
        box(fig, x, y + 0.113, 0.130, 0.017, color, r=0.003, z=2 + i)
    text(fig, 0.028, 0.922, "Science", fontsize=5.5, color="white",
         weight="bold", zorder=7)
    text(fig, 0.028, 0.880, "Kidney multiome-based\ngenetic scorecard reveals\n"
         "convergent coding and\nregulatory variants", fontsize=4.7,
         weight="bold", zorder=7)
    text(fig, 0.028, 0.824, "Liu H, … Susztak K.\nScience 2025\n"
         "doi:10.1126/science.adp4753", fontsize=4.5, color=MUTED, zorder=7)
    text(fig, 0.100, 0.775, "Title, DOI, PMID\nor citation", fontsize=6,
         ha="center", color=MUTED)

    arrow(fig, 0.182, 0.86, 0.232, 0.86)
    text(fig, 0.207, 0.885, "igvf-replicate-\npaper", fontsize=5.6, ha="center",
         color=MUTED, style="italic")

    # IGVF Agent server box.
    x0, x1 = 0.236, 0.630
    box(fig, x0, bot, x1 - x0, top - bot, "white", ec="#7B8794", ls=(0, (3, 2)),
        r=0.012)
    logo(fig, x0 + 0.012, top - 0.033, 0.024)
    text(fig, x0 + 0.042, top - 0.018, "IGVF Agent", fontsize=8, weight="bold")

    rows = [
        ("liu2025-open4gene: hurdle negative-binomial peak→gene test\n"
         "(port of R/Open4Gene.R; zero-model β and P identical)"),
        ("encode-re2g: ENCODE-rE2G enhancer→gene model, 9 pretrained\n"
         "(max |Δscore| 7.8e-16 vs the authors' K562 test)"),
        ("esig_net_predict: eSIG-Net inference on the released checkpoint\n"
         "(bit-identical to the authors' test.sh)"),
    ]
    ry = top - 0.060
    for s in rows:
        box(fig, x0 + 0.012, ry - 0.022, x1 - x0 - 0.024, 0.040, TOOL, r=0.006)
        text(fig, x0 + 0.020, ry - 0.002, s, fontsize=5.5)
        ry -= 0.047

    yb = bot + 0.034
    bw = (x1 - x0 - 0.036) / 2
    box(fig, x0 + 0.012, yb - 0.022, bw, 0.040, DATA, r=0.006)
    text(fig, x0 + 0.012 + bw / 2, yb - 0.002,
         "Deposited data: Figshare,\nSupp. Tables, Source Data",
         fontsize=5.5, ha="center")
    box(fig, x0 + 0.024 + bw, yb - 0.022, bw, 0.040, CODE, r=0.006)
    text(fig, x0 + 0.024 + bw + bw / 2, yb - 0.002,
         "Authors' code at pinned\ncommits, run as reference", fontsize=5.5,
         ha="center")

    # Legend chips under the server box.
    lx = x0 + 0.02
    for lab, c in [("Ported methods (CLI tools)", TOOL),
                   ("Data", DATA), ("Reference code", CODE)]:
        box(fig, lx, bot - 0.022, 0.010, 0.010, c, r=0.002)
        text(fig, lx + 0.014, bot - 0.017, lab, fontsize=5.5)
        lx += 0.014 + 0.0058 * len(lab) + 0.02

    arrow(fig, 0.636, 0.87, 0.680, 0.87)
    text(fig, 0.658, 0.900, "Connect\nto an\nagent", fontsize=5.6, ha="center",
         color=MUTED, style="italic")

    # Agent chat.
    ax0, ax1 = 0.685, 0.990
    box(fig, ax0, bot, ax1 - ax0, top - bot, "#F7F9FB", ec="#CBD2D9", r=0.012)
    logo(fig, (ax0 + ax1) / 2 - 0.030, top - 0.030, 0.018)
    text(fig, (ax0 + ax1) / 2 - 0.008, top - 0.016, "IGVF Agent", fontsize=6.5,
         weight="bold")
    box(fig, ax0 + 0.035, top - 0.078, ax1 - ax0 - 0.045, 0.040, USER, r=0.008)
    text(fig, ax0 + 0.042, top - 0.058,
         "User: Replicate Liu et al. 2025,\nScience (doi:10.1126/science.adp4753)",
         fontsize=5.5)
    box(fig, ax0 + 0.010, bot + 0.010, ax1 - ax0 - 0.020, 0.130, AGENT, r=0.008)
    text(fig, ax0 + 0.017, bot + 0.133,
         "Agent: 11 of 16 analyses reproduced.\nFig. 3D re-drawn: R = 0.520 "
         "(paper 0.51).\n2 blocked: individual-level data\n(controlled access).",
         fontsize=5.5, va="top")
    img = HERE / PAPERS[0][0] / PAPERS[0][3]
    if img.is_file():
        ax = fig.add_axes([ax0 + 0.085, bot + 0.013, 0.13, 0.066], zorder=10)
        ax.imshow(mpimg.imread(img))
        ax.axis("off")


# ---------------------------------------------------------------- panel b
def panel_b(fig):
    panel_label(fig, 0.01, 0.712, "b")
    top, bot = 0.700, 0.495

    # The chain of tools, as the skill runs them.
    x0, x1 = 0.030, 0.565
    box(fig, x0, bot, x1 - x0, top - bot, "#1F2933", r=0.010)
    for i, c in enumerate(["#FF5F56", "#FFBD2E", "#27C93F"]):
        fig.add_artist(plt.Circle((x0 + 0.014 + i * 0.013, top - 0.012), 0.004,
                                  transform=fig.transFigure, color=c, zorder=3))
    text(fig, x0 + 0.06, top - 0.012,
         "/igvf-replicate-paper  Liu et al. 2025 Science", fontsize=5.6,
         family=MONO, color="#E4E7EB")
    steps = [
        ("bench resolve", "PMID 39913582 · PMC12013656"),
        ("bench harvest", "Figshare 26299093 · 2 repos pinned"),
        ("bench plan --seed", "16 analyses, one per figure panel"),
        ("paper-code pipeline", "authors' scripts = reference"),
        ("port register", "liu2025-open4gene, liu2025-gwas-loci"),
        ("bench verify-port", "port vs R pscl::hurdle: P identical"),
        ("bench run", "recompute Figs. 1–6 from tables"),
        ("figures redraw", "Fig. 3D redrawn, side by side"),
        ("bench score", "reproduction: incomplete, 11/16"),
    ]
    y = top - 0.032
    for i, (cmd, out) in enumerate(steps, 1):
        text(fig, x0 + 0.012, y, f"{i}", fontsize=5.6, family=MONO,
             color="#7B8794")
        text(fig, x0 + 0.026, y, f"igvfagent {cmd}", fontsize=5.6, family=MONO,
             color="#9FE5A5")
        text(fig, x0 + 0.245, y, f"→ {out}", fontsize=5.2, family=MONO,
             color="#E4E7EB")
        y -= 0.0182
    text(fig, x0 + 0.012, bot + 0.010,
         "loop 4–9 until every analysis is reproduced or blocked",
         fontsize=5.4, family=MONO, color="#F7C948", va="bottom")
    text(fig, (x0 + x1) / 2, bot - 0.012, "Chain of tools", fontsize=6.5,
         ha="center", color=MUTED)

    arrow(fig, 0.570, 0.60, 0.592, 0.60)

    # The report it returns.
    rx0, rx1 = 0.597, 0.990
    box(fig, rx0, bot, rx1 - rx0, top - bot, "white", ec="#CBD2D9", r=0.010)
    text(fig, rx0 + 0.012, top - 0.014, "Output: Benchmarks/liu2025_kidney_"
         "multiome/README.md", fontsize=5.6, weight="bold")
    report = [
        ("Coverage", "11/16 analyses reproduced; 43/46 checks pass"),
        ("Ports", "liu2025-open4gene, liu2025-gwas-loci, each\n"
                  "verified against the authors' own code"),
        ("Exact", "ASE genes 8,573 / 8,725 / 10,398; bASA peaks\n"
                  "19,083; Open4Gene links 142,452; coding 1,363"),
        ("Differs", "median peaks per gene 4 vs 8 (author\n"
                    "re-checking); non-MHC variants +2.9%"),
        ("Blocked", "RASQUAL, snASA, CARMA: individual-level\n"
                    "genotypes, controlled access"),
        ("Asks", "a table of exact data needed from the authors,\n"
                 "ordered by how much each item unlocks"),
    ]
    y = top - 0.036
    for k, v in report:
        text(fig, rx0 + 0.012, y, k, fontsize=5.6, weight="bold", va="top",
             color="#3E4C59")
        text(fig, rx0 + 0.082, y, v, fontsize=5.5, va="top")
        y -= 0.0285 if "\n" in v else 0.018
    text(fig, (rx0 + rx1) / 2, bot - 0.012, "Output", fontsize=6.5,
         ha="center", color=MUTED)


# ---------------------------------------------------------------- panel c
def gschwind_bars(ax, color):
    import numpy as np
    x = np.arange(len(GSCHWIND_AUPRC))
    paper = [p for _, p, _ in GSCHWIND_AUPRC]
    ours = [o for _, _, o in GSCHWIND_AUPRC]
    ax.bar(x - 0.19, paper, 0.36, color="white", ec=color, lw=0.8,
           label="paper")
    ax.bar(x + 0.19, ours, 0.36, color=color, label="IGVF Agent")
    for xi, o in zip(x, ours):
        ax.text(xi + 0.19, o + 0.015, f"{o:.2f}", ha="center", fontsize=5)
    ax.set_xticks(x, [m for m, *_ in GSCHWIND_AUPRC], fontsize=5.2)
    ax.set_ylim(0, 0.98)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylabel("AUPRC (CRISPR benchmark)", fontsize=6)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2,
              handlelength=1.2, columnspacing=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_c(fig):
    panel_label(fig, 0.01, 0.458, "c")
    box(fig, 0.30, 0.435, 0.40, 0.022, USER, r=0.006)
    text(fig, 0.50, 0.446, "User: Replicate this paper with IGVF Agent",
         fontsize=6.5, ha="center")

    left0, colw, gap = 0.075, 0.255, 0.075
    for i, (pid, header, color, img, caption) in enumerate(PAPERS):
        x = left0 + i * (colw + gap)
        text(fig, x - 0.04, 0.418, header, fontsize=7, weight="bold",
             color=color)

        # Row 1: a redrawn panel.
        ax = fig.add_axes([x, 0.300, colw, 0.098] if not img
                          else [x, 0.272, colw, 0.112])
        if img:
            ax.imshow(mpimg.imread(HERE / pid / img))
            ax.axis("off")
        else:
            gschwind_bars(ax, color)
        text(fig, x + colw / 2, 0.400 + 0.002, caption, fontsize=5.6, ha="center",
             color=INK)

        # Row 2: every quoted number, ours vs paper.
        rows = NUMBERS[pid]
        n_ok = sum(abs(o - p) / abs(p) <= 0.05 for _, p, o in rows)
        h = 0.0082 * max(len(r) for r in NUMBERS.values())
        ax = fig.add_axes([x + 0.085, 0.240 - h, colw - 0.085, h])
        ax.axvspan(-5, 5, color="#E4E7EB", lw=0)
        ax.axvline(0, color="#9AA5B1", lw=0.5)
        for j, (lab, p, o) in enumerate(rows):
            d = 100 * (o - p) / p
            dd = max(min(d, 60), -60)
            out = abs(d) > 5
            c = color if not out else PARTIAL
            if lab == "median peaks/gene":
                c = "#C0392B"
            ax.plot(dd, j, "o", ms=2.8, color=c, mec="none")
            if out:
                ax.text(dd + 5, j, f"{d:+.0f}%", ha="left", va="center",
                        fontsize=5, color=c)
        ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=5)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xlim(-65, 65)
        ax.set_xticks([-50, 0, 50])
        ax.tick_params(axis="x", labelsize=5.5)
        ax.set_xlabel("IGVF Agent vs paper (%)", fontsize=5.8, labelpad=1)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        text(fig, x + colw / 2, 0.255,
             f"{n_ok} of {len(rows)} numbers within 5% of the paper",
             fontsize=6, ha="center", weight="bold")

        # Row 3: analysis coverage.
        rep, part, blk = COVERAGE[pid]
        tot = rep + part + blk
        ax = fig.add_axes([x, 0.000, colw, 0.018])
        left = 0
        for val, c in [(rep, color), (part, PARTIAL), (blk, BLOCKED)]:
            if val:
                ax.barh(0, val, left=left, color=c, height=1)
                ax.text(left + val / 2, 0, str(val), ha="center", va="center",
                        fontsize=6, color="white" if c == color else INK,
                        weight="bold")
            left += val
        ax.set_xlim(0, 16)
        ax.set_ylim(-0.5, 0.5)
        ax.axis("off")
        text(fig, x + colw / 2, 0.030,
             f"{rep} of {tot} analyses reproduced", fontsize=6, ha="center",
             weight="bold")

    # Shared legend for the coverage bars.
    lx = 0.20
    for lab, c in [("reproduced (paper colour)", "#52606D"),
                   ("partial / not reproduced", PARTIAL),
                   ("blocked or pending", BLOCKED)]:
        box(fig, lx, -0.030, 0.012, 0.009, c, r=0.001)
        text(fig, lx + 0.016, -0.0255, lab, fontsize=5.8)
        lx += 0.21


def main():
    fig = plt.figure(figsize=(7.2, 10.0), facecolor="white")
    panel_a(fig)
    panel_b(fig)
    panel_c(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    print(OUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
