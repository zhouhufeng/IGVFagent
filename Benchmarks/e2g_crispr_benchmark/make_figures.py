#!/usr/bin/env python3
"""Figures for the E2G CRISPR benchmark (Part A: K562 vs CRISPR ground truth;
Part B: coronary-artery scE2G cross-cell-type characterization)."""
import gzip, json, os
from collections import Counter, defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# Published / third-party reference numbers (see README for provenance)
REF = {"ENCODE-rE2G (published, K562)": 0.634, "ABC (published, K562)": 0.612}

METHOD_LABEL = {
    "rE2G_ext:FullModel.Score": "ENCODE-rE2G extended",
    "rE2G_ext:FullModel_minus_EP300.Score": "rE2G ext (−EP300)",
    "rE2G_base:Full.Score": "ENCODE-rE2G (base)",
    "Score.ignoreTPM": "scE2G (ignoreTPM)",
    "Score": "scE2G (Score)",
    "ARC.E2G.Score": "ARC-E2G",
    "ABC.Score": "ABC (baseline)",
}
COLORS = {"rE2G_ext:FullModel.Score": "#7b2d8b",
          "rE2G_ext:FullModel_minus_EP300.Score": "#a86bc4",
          "rE2G_base:Full.Score": "#b8002e",
          "Score": "#1b6ca8", "Score.ignoreTPM": "#0aa3a3",
          "ARC.E2G.Score": "#e07b39", "ABC.Score": "#888888"}


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", dpi=130)
    plt.close(fig)
    print("wrote", name)


def _load_all_methods():
    """Merge scE2G-table methods + rE2G base/extended into one ordered dict."""
    RES = os.path.join(HERE, "results")
    d = json.load(open(os.path.join(RES, "k562_scE2G_vs_crispr.json")))
    methods = dict(d["methods"])
    meta = {"n_positive": d["n_positive"], "n_pairs": d["n_pairs"],
            "match_rate": d["match_rate"]}
    for tag, fn in (("rE2G_base", "k562_rE2G_base_vs_crispr.json"),
                    ("rE2G_ext", "k562_rE2G_extended_vs_crispr.json")):
        p = os.path.join(RES, fn)
        if os.path.exists(p):
            dd = json.load(open(p))
            for m, v in dd["methods"].items():
                methods[f"{tag}:{m}"] = v
    # order by the label table, keeping only known methods
    ordered = {k: methods[k] for k in METHOD_LABEL if k in methods}
    return ordered, meta


def part_a():
    methods, d = _load_all_methods()

    # fig1 — PR curves
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    for m, v in methods.items():
        c = v["curve"]
        ax.plot(c["recall"], c["precision"], label=f"{METHOD_LABEL[m]} (AUPRC={v['auprc']:.3f})",
                color=COLORS[m], lw=2)
    ax.axvline(0.70, ls="--", color="grey", lw=1, alpha=.7)
    ax.text(0.70, 0.02, " 70% recall", color="grey", fontsize=8, rotation=90, va="bottom")
    base = d["n_positive"] / d["n_pairs"]
    ax.axhline(base, ls=":", color="red", lw=1, alpha=.6)
    ax.text(0.01, base + .01, f"random baseline = {base:.3f}", color="red", fontsize=8)
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_title(f"K562 E2G vs ENCODE CRISPR benchmark\n{d['n_positive']} positives / "
                 f"{d['n_pairs']} pairs, {d['match_rate']*100:.0f}% element-overlap")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.legend(fontsize=8, loc="upper right")
    save(fig, "fig1_k562_pr_curves")

    # fig2 — AUPRC bars vs published reference
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    labels = [METHOD_LABEL[m] for m in methods]
    vals = [methods[m]["auprc"] for m in methods]
    cols = [COLORS[m] for m in methods]
    bars = ax.bar(labels, vals, color=cols)
    for b, val in zip(bars, vals):
        ax.text(b.get_x() + b.get_width()/2, val + .008, f"{val:.3f}", ha="center", fontsize=9)
    for name, y in REF.items():
        ax.axhline(y, ls="--", lw=1.2, alpha=.8,
                   color="#b8002e" if "rE2G" in name else "#555")
        ax.text(len(labels)-0.4, y + .004, name, fontsize=7.5, ha="right",
                color="#b8002e" if "rE2G" in name else "#555")
    ax.set_ylabel("AUPRC"); ax.set_ylim(0, 0.75)
    ax.set_title("AUPRC on K562 CRISPR benchmark — IGVFagent vs published references")
    plt.xticks(rotation=15, ha="right")
    save(fig, "fig2_k562_auprc_bars")

    # fig3 — precision@70% recall
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    vals = [methods[m]["precision_at_70_recall"] for m in methods]
    bars = ax.bar(labels, vals, color=cols)
    for b, val in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, val + .006, f"{val:.3f}", ha="center", fontsize=9)
    ax.set_ylabel("Precision at 70% recall"); ax.set_ylim(0, 0.6)
    ax.set_title("Precision at 70% recall — K562 CRISPR benchmark")
    plt.xticks(rotation=15, ha="right")
    save(fig, "fig3_k562_precision_at_70")


def part_b():
    path = os.path.join(ROOT, "Data", "E2G_benchmark", "scE2G_coronary_artery",
                        "IGVFFI8252JBBA.tsv.gz")
    f = gzip.open(path, "rt")
    hdr = None
    for line in f:
        if line.startswith("#"): continue
        hdr = line.rstrip("\n").split("\t"); break
    ix = {c: i for i, c in enumerate(hdr)}
    ct = Counter(); ec = Counter(); scores = []; dist = []
    ct_score = defaultdict(list)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        cell = p[ix["CellType"]].replace(" (coronary artery)", "")
        ct[cell] += 1
        ec[p[ix["ElementClass"]]] += 1
        try:
            s = float(p[ix["Score"]]); scores.append(s); ct_score[cell].append(s)
        except ValueError:
            pass
        try:
            mid = (int(p[ix["ElementStart"]]) + int(p[ix["ElementEnd"]])) / 2
            dist.append(abs(mid - int(p[ix["GeneTSS"]])))
        except ValueError:
            pass

    # fig4 — predictions per cell type
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    items = ct.most_common()
    ax.barh([k for k, _ in items][::-1], [v for _, v in items][::-1], color="#2c7fb8")
    ax.set_xlabel("E2G predictions"); ax.set_title("scE2G coronary-artery predictions per cell type (IGVFDS7543BJBA)")
    for i, (k, v) in enumerate(items[::-1]):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    save(fig, "fig4_coronary_per_celltype")

    # fig5 — score distribution + element class
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    a1.hist(scores, bins=50, color="#41ab5d", edgecolor="white")
    a1.axvline(0.5, ls="--", color="grey"); a1.set_xlabel("scE2G Score"); a1.set_ylabel("pairs")
    a1.set_title(f"Score distribution (n={len(scores):,})")
    eitems = ec.most_common()
    a2.bar([k for k, _ in eitems], [v for _, v in eitems],
           color=["#756bb1", "#9e9ac8", "#cbc9e2"][:len(eitems)])
    for i, (k, v) in enumerate(eitems):
        a2.text(i, v, f"{100*v/sum(ec.values()):.0f}%", ha="center", va="bottom", fontsize=9)
    a2.set_ylabel("pairs"); a2.set_title("Element class")
    save(fig, "fig5_coronary_score_elementclass")

    # fig6 — distance to TSS decay
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    kb = [x/1000 for x in dist if x <= 300000]
    ax.hist(kb, bins=60, color="#dd8452", edgecolor="white")
    med = sorted(dist)[len(dist)//2] / 1000
    ax.axvline(med, ls="--", color="black"); ax.text(med, ax.get_ylim()[1]*.9, f" median {med:.0f} kb", fontsize=9)
    ax.set_xlabel("|element midpoint − TSS| (kb, ≤300 kb shown)"); ax.set_ylabel("pairs")
    ax.set_title("Enhancer→gene distance decay (coronary artery)")
    save(fig, "fig6_coronary_distance_decay")


if __name__ == "__main__":
    part_a()
    part_b()
    print("all figures ->", FIG)
