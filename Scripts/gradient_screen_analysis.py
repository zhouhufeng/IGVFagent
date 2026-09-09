#!/usr/bin/env python3
"""Sorted expression-bin CRISPR screens: score constructs across bins A-F.

WHAT THIS IS FOR. 554 MeasurementSets on the Portal are FACS screens sorted
into LETTERED expression bins rather than two tails, and nothing could
analyse them. Two families, both Engreitz lab:

  398 sets  Variant-EFFECTS / prime editing, endogenous or exogenous allelic
            sequencing, e.g. PPIF promoter variants in THP-1 and Jurkat
  156 sets  CRISPR FlowFISH ("CRUDO"), gRNA sequencing, e.g. KITLG, MYC,
            FAM3C with an auxin time course

crispr_screen_analysis handles the OTHER shape -- a bottom20%/top20% tail
sort -- and cannot be stretched to cover this one. A tail sort asks "is this
construct enriched in the low tail or the high tail"; a six-bin gradient
asks "what is this construct's mean expression", which is a different
statistic computed from all bins at once.

WHY THE ALIAS IS PARSED BY TOKEN, NOT BY POSITION. Twice now a positional
regex has silently found a fraction of a screen: `(?P<side>bottom|top)`
matched none of the Sherwood `Bot20` bins, and `_Rep(\\d+)_` matched none of
these `-BioRep3-FFrep3-` ones. The real aliases vary in separator and in
token order:

  220228_8merInsertion-THP1-Amp_PPIF_promoter-BioRep3-FFrep3-BinF_AmpliconSequencing
  measurement_CRUDO_KITLG-Auxin6hrs-FF4-BinE-PCR2_S95
  measurement_CRUDO_MYC-Auxin0hrs-FF1-BinE-PCR4_Rep1_S23
  measurement_CRUDO_FAM3C_Auxin6hrs-FF2-BinB-PCR3-Rep1_S200

So each axis is found wherever it sits, and the SERIES is whatever remains
once the recognised axes are removed. A convention this tool has never seen
still yields a series and a bin, instead of yielding nothing.

    igvfagent gradient-screen discover IGVFDS8710ZSOZ
    igvfagent gradient-screen analyze  IGVFDS8710ZSOZ --label kitlg

WHAT THE SCORE MEANS. Each construct gets a mean-bin score per sort:

    score = sum_b (freq_b * rank_b) / sum_b freq_b        rank A=1 .. F=6

freq_b is the construct's share of bin b's assigned reads, so bins sequenced
to different depths contribute on equal terms. The score is that construct's
centre of mass along the expression axis: higher = the cells carrying it
sorted into higher-expression bins. Scores are then centred on the
non-targeting controls where the library has them, so the units are "bins
away from a construct that does nothing", and compared across sorts with the
same moderated t used by the tail-sort screens (Scripts/_stats.py).
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import raw_data_pipeline as rp                              # noqa: E402
from _stats import benjamini_hochberg, moderated_t, _percentile  # noqa: E402

ROOT = rp.ROOT
OUT_DIR = ROOT / "Docs" / "GradientScreen"
LOG_DIR = ROOT / "Docs" / "Logs"

# Bin letters are the expression axis. A=lowest .. F=highest, which is the
# order the Engreitz aliases use and the order the sort gates are numbered.
BIN_RANK = {c: i + 1 for i, c in enumerate("ABCDEFGH")}

# Each axis is looked for anywhere in the alias. Order matters only within a
# pair: BioRep must be tried before Rep, and FFrep before FF, or the shorter
# name matches the longer token's tail and the number comes out wrong.
_AXES = (
    ("bin",      re.compile(r"Bin([A-Z])(?![a-z])")),
    ("bio_rep",  re.compile(r"BioRep(\d+)|(?<![A-Za-z])Rep(\d+)", re.I)),
    ("flow_rep", re.compile(r"FFrep(\d+)|(?<![A-Za-z])FF(\d+)", re.I)),
    ("pcr_rep",  re.compile(r"PCR(\d+)", re.I)),
)
# Trailing sample index from the sequencing submission: not an experimental
# axis, but it must be stripped or it lands in the series and splits a screen
# into one series per sample.
_SAMPLE_SUFFIX = re.compile(r"[-_]S\d+$")
# Readout suffixes that describe how the bin was sequenced, not which bin.
_READOUT_SUFFIX = re.compile(
    r"[-_](?:AmpliconSequencing|GuideSequencing|gRNASequencing)$", re.I)


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    p = LOG_DIR / f"gradient_screen_{time.strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)-6s %(message)s",
        handlers=[logging.FileHandler(p), logging.StreamHandler(sys.stdout)])
    return p


def parse_bin_alias(alias: str) -> "Optional[dict]":
    """Pull the experimental axes out of a submitter alias.

    Returns {series, bin, rank, bio_rep, flow_rep, pcr_rep} or None when no
    lettered bin is present -- which is the one thing this shape cannot do
    without.
    """
    if not alias:
        return None
    a = alias.split(":", 1)[-1].strip()
    a = _SAMPLE_SUFFIX.sub("", a)
    a = _READOUT_SUFFIX.sub("", a)

    found: "dict[str, Optional[str]]" = {}
    series = a
    for name, rx in _AXES:
        m = rx.search(series)
        if not m:
            found[name] = None
            continue
        found[name] = next(g for g in m.groups() if g is not None)
        # Remove the token so it cannot end up in the series, and close the
        # gap left behind so "A-BinC-B" does not become "A--B".
        series = (series[:m.start()] + series[m.end():])
    if found["bin"] is None:
        return None
    series = re.sub(r"[-_]{2,}", "-", series).strip("-_ ")
    return {
        "series": series,
        "bin": found["bin"].upper(),
        "rank": BIN_RANK.get(found["bin"].upper(), 0),
        "bio_rep": int(found["bio_rep"]) if found["bio_rep"] else None,
        "flow_rep": int(found["flow_rep"]) if found["flow_rep"] else None,
        "pcr_rep": int(found["pcr_rep"]) if found["pcr_rep"] else None,
    }


def _aliases(fs: dict) -> "list[str]":
    return [str(x) for x in (fs.get("aliases") or [])]


def _parse_any(fs: dict) -> "Optional[dict]":
    """First alias of a set that yields a bin. A CRUDO set carries one alias
    per PCR replicate, all agreeing on series and bin."""
    for a in _aliases(fs):
        p = parse_bin_alias(a)
        if p:
            p["alias"] = a.split(":", 1)[-1]
            return p
    return None


def discover_screen(accession: str) -> dict:
    """Every sorted bin of the screen the given accession belongs to.

    Siblings are found through the shared ConstructLibrarySet and then
    filtered to the same series. The library alone is not enough: the PPIF
    guide library is used by THP-1 AND Jurkat experiments, so searching by
    library returned 175 sets spanning cell types and conditions that must
    not be pooled.
    """
    st, me = rp.portal_json(f"/measurement-sets/{accession}/?format=json")
    if not me:
        return {"error": f"{accession} was not returned by the Portal "
                         f"(status {st})."}
    mine = _parse_any(me)
    if not mine:
        return {"error": f"{accession} has no lettered-bin alias, so its "
                         f"sorted screen cannot be identified. Aliases were: "
                         f"{_aliases(me) or 'none'}"}

    libs = [(x.get("accession") if isinstance(x, dict) else str(x).rstrip("/").split("/")[-1])
            for x in (me.get("construct_library_sets") or [])]
    rows: "dict[str, dict]" = {}
    for lib in [l for l in libs if l]:
        st, o = rp.portal_json(
            "/search/?type=MeasurementSet&construct_library_sets.accession="
            f"{lib}&limit=all&field=accession&field=aliases&field=lab"
            "&field=preferred_assay_titles&field=crispr_screen_readout"
            "&format=json")
        for r in (o or {}).get("@graph", []):
            rows[r["accession"]] = r
    rows.setdefault(accession, me)

    bins = []
    for acc, r in rows.items():
        p = _parse_any(r)
        if p and p["series"] == mine["series"]:
            bins.append({"accession": acc, **p})
    # Sort by the axes so output order is the experiment's order, not the
    # Portal's.
    bins.sort(key=lambda b: (b["bio_rep"] or 0, b["flow_rep"] or 0, b["rank"]))
    sorts = sorted({(b["bio_rep"], b["flow_rep"]) for b in bins},
                   key=lambda t: (t[0] or 0, t[1] or 0))
    return {
        "series": mine["series"], "query_set": accession,
        "assay": ", ".join(me.get("preferred_assay_titles") or []) or "?",
        "readout": me.get("crispr_screen_readout") or "?",
        "lab": (me.get("lab") or {}).get("title", "?")
                if isinstance(me.get("lab"), dict) else "?",
        "libraries": [l for l in libs if l],
        "bins": bins,
        "bin_letters": sorted({b["bin"] for b in bins}),
        "sorts": sorts,
    }


# ─── Scoring ────────────────────────────────────────────────────────────────

def sort_scores(per_bin: "dict[str, Counter]", bins: "list[dict]",
                 min_count: int) -> "tuple[dict[str, float], dict[str, float], dict[str, int]]":
    """Mean-bin score for every construct in ONE sort.

    Returns (score, counting_variance, total_reads) keyed by construct.

    freq_b is the construct's share of bin b's assigned reads, so a bin
    sequenced twice as deeply does not count twice as much. The score is the
    frequency-weighted mean bin rank -- the construct's centre of mass along
    the expression axis.

    counting_variance is the multinomial sampling variance of that weighted
    mean, sum_b p_b (rank_b - score)^2 / N, with N the construct's total
    reads. It is the floor on how precisely this construct's score could
    have been measured at the depth it actually got: a construct seen 20
    times cannot have a score good to 0.01 of a bin however tidily its
    replicates agree.
    """
    totals = {b["accession"]: max(sum(per_bin.get(b["accession"], Counter()).values()), 1)
              for b in bins}
    freq: "dict[str, dict[int, float]]" = defaultdict(dict)
    reads: "dict[str, int]" = defaultdict(int)
    for b in bins:
        c = per_bin.get(b["accession"])
        if not c:
            continue
        t = totals[b["accession"]]
        for gid, n in c.items():
            freq[gid][b["rank"]] = freq[gid].get(b["rank"], 0.0) + n / t
            reads[gid] += n

    score: "dict[str, float]" = {}
    cvar: "dict[str, float]" = {}
    for gid, by_rank in freq.items():
        if reads[gid] < min_count:
            continue
        w = sum(by_rank.values())
        if w <= 0:
            continue
        s = sum(r * f for r, f in by_rank.items()) / w
        score[gid] = s
        p = {r: f / w for r, f in by_rank.items()}
        var = sum(pi * (r - s) ** 2 for r, pi in p.items())
        cvar[gid] = var / max(reads[gid], 1)
    return score, cvar, dict(reads)


def score_gradient(per_sort: "dict[tuple, dict]", target_of: "dict[str, str]",
                    type_of: "dict[str, str]") -> "tuple[list[dict], dict]":
    """Combine per-sort mean-bin scores into per-construct effects.

    Each sort is centred on its own non-targeting controls before the sorts
    are combined. The gates are re-drawn for every sort, so an uncentred
    score carries the gate placement of that day as if it were biology;
    subtracting the controls of the same sort removes it and puts every score
    in the same units -- bins away from a construct that does nothing.

    Without controls the sort's own median score is used instead, and the
    output says so, because then a whole-library shift cannot be told from a
    gate shift.
    """
    centred: "dict[str, list[float]]" = defaultdict(list)
    cvars: "dict[str, list[float]]" = defaultdict(list)
    reads: "dict[str, int]" = defaultdict(int)
    centring = []
    for key, d in sorted(per_sort.items(), key=lambda kv: (kv[0][0] or 0, kv[0][1] or 0)):
        score, cvar, rd = d["score"], d["cvar"], d["reads"]
        ctrl = [s for g, s in score.items()
                if any(w in (type_of.get(g, "") or "").lower()
                       for w in ("non-targeting", "nontargeting", "control"))]
        if len(ctrl) >= 5:
            base = sum(ctrl) / len(ctrl)
            centring.append({"sort": list(key), "on": "controls",
                              "n": len(ctrl), "baseline": round(base, 4)})
        else:
            vals = sorted(score.values())
            base = _percentile(vals, 0.5) if vals else 0.0
            centring.append({"sort": list(key), "on": "library median",
                              "n": len(ctrl), "baseline": round(base, 4)})
        for g, s in score.items():
            centred[g].append(s - base)
            cvars[g].append(cvar.get(g, 0.0))
            reads[g] += rd.get(g, 0)

    # A screen-wide floor, from the constructs with enough sorts to estimate
    # a spread at all -- the same device the tail-sort screens use.
    sds = []
    for g, vals in centred.items():
        if len(vals) >= 3:
            m = sum(vals) / len(vals)
            sds.append(math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1)))
    sd_floor = _percentile([x for x in sds if x > 0], 0.10) if sds else 0.0

    rows = []
    for g, vals in centred.items():
        k = len(vals)
        sd_count = math.sqrt(sum(cvars[g]) / k) / math.sqrt(k) if k else 0.0
        st = moderated_t(vals, sd_floor, sd_count)
        rows.append({
            "construct": g, "target": target_of.get(g, g),
            "construct_type": type_of.get(g, ""),
            "n_sorts": k, "reads": reads[g],
            "delta_bins": round(st["mean"], 4),
            "sd": None if st["sd"] is None else round(st["sd"], 4),
            "sd_counting": round(sd_count, 4),
            "sd_used": None if st["sd_used"] is None else round(st["sd_used"], 4),
            "t": round(st["t"], 4), "df": st["df"], "p_value": st["p_value"],
        })
    for r, q in zip(rows, benjamini_hochberg([r["p_value"] for r in rows])):
        r["fdr"] = q
    rows.sort(key=lambda r: (r["p_value"], -abs(r["delta_bins"])))
    return rows, {"sd_floor": round(sd_floor, 4), "n_ref_sd": len(sds),
                  "centring": centring}


# ─── QC and plots ───────────────────────────────────────────────────────────

def concordance(per_sort: "dict[tuple, dict]") -> "list[dict]":
    """Pearson r between every pair of sorts, on their shared constructs.

    This is the QC that matters for a gradient screen: if two sorts of the
    same library disagree about which constructs sit high, no per-construct
    score from them means anything, however small its p-value.
    """
    keys = sorted(per_sort, key=lambda t: (t[0] or 0, t[1] or 0))
    out = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            sa, sb = per_sort[a]["score"], per_sort[b]["score"]
            shared = sorted(set(sa) & set(sb))
            if len(shared) < 3:
                out.append({"a": list(a), "b": list(b), "n": len(shared),
                             "r": None})
                continue
            xs = [sa[g] for g in shared]
            ys = [sb[g] for g in shared]
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            sxx = sum((x - mx) ** 2 for x in xs)
            syy = sum((y - my) ** 2 for y in ys)
            sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            r = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else None
            out.append({"a": list(a), "b": list(b), "n": len(shared),
                         "r": None if r is None else round(r, 4)})
    return out


def make_plots(out: Path, rows: "list[dict]", qc: "list[dict]",
                per_sort: "dict[tuple, dict]", conc: "list[dict]",
                bin_letters: "list[str]", series: str) -> "list[Path]":
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return []
    out.mkdir(parents=True, exist_ok=True)
    made = []

    # --- QC panel: depth, assignment rate, library coverage per library ---
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    labels = [f"{q['bin']}{'' if q['flow_rep'] is None else '/f'+str(q['flow_rep'])}"
              for q in qc]
    ax[0].bar(range(len(qc)), [q["reads_assigned"] for q in qc], color="#4C72B0")
    ax[0].set_title("Reads assigned per library")
    ax[0].set_ylabel("reads")
    ax[1].bar(range(len(qc)), [100 * q["assignment_rate"] for q in qc],
              color="#55A868")
    ax[1].set_title("Assignment rate")
    ax[1].set_ylabel("% of reads")
    ax[1].axhline(2, ls="--", lw=1, c="k")
    ax[2].bar(range(len(qc)), [100 * q["library_coverage"] for q in qc],
              color="#8172B2")
    ax[2].set_title("Library coverage")
    ax[2].set_ylabel("% of constructs seen")
    for a in ax:
        a.set_xticks(range(len(qc)))
        try:
            a.set_xticklabels(labels, rotation=90, fontsize=6)
        except Exception:
            pass
    fig.suptitle(f"QC — {series}", fontsize=10)
    fig.tight_layout()
    p = out / "qc_panel.png"
    fig.savefig(p, dpi=140); plt.close(fig); made.append(p)

    # --- effect distribution, controls overlaid ---
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    eff = np.array([r["delta_bins"] for r in rows])
    is_ctrl = np.array([any(w in (r["construct_type"] or "").lower()
                             for w in ("non-targeting", "nontargeting", "control"))
                         for r in rows])
    ax[0].hist(eff[~is_ctrl], bins=60, color="#4C72B0", alpha=.85,
               label=f"constructs (n={int((~is_ctrl).sum())})")
    if is_ctrl.any():
        ax[0].hist(eff[is_ctrl], bins=30, color="#C44E52", alpha=.85,
                   label=f"controls (n={int(is_ctrl.sum())})")
    ax[0].axvline(0, ls=":", lw=1, c="grey")
    ax[0].set_xlabel("delta mean bin (bins away from control)")
    ax[0].set_ylabel("constructs")
    ax[0].set_title("Effect size distribution")
    ax[0].legend(fontsize=8)

    fdr = np.array([max(r["fdr"], 1e-12) for r in rows])
    ax[1].scatter(eff[~is_ctrl], -np.log10(fdr[~is_ctrl]), s=12,
                  c="#4C72B0", alpha=.6, label="constructs")
    if is_ctrl.any():
        ax[1].scatter(eff[is_ctrl], -np.log10(fdr[is_ctrl]), s=26,
                      c="#C44E52", alpha=.9, label="controls")
    ax[1].axhline(-math.log10(0.05), ls="--", lw=1, c="k")
    ax[1].axvline(0, ls=":", lw=1, c="grey")
    ax[1].set_xlabel("delta mean bin")
    ax[1].set_ylabel("-log10 FDR")
    ax[1].set_title("positive = sorted into HIGHER expression bins")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    p = out / "effects.png"
    fig.savefig(p, dpi=140); plt.close(fig); made.append(p)

    # --- sort-to-sort concordance ---
    rs = [c["r"] for c in conc if c["r"] is not None]
    if rs:
        fig, ax = plt.subplots(figsize=(max(5, len(conc) * 0.32), 3.8))
        ax.bar(range(len(rs)), rs, color="#937860")
        ax.axhline(0, lw=1, c="k")
        ax.set_ylim(min(-0.1, min(rs) - 0.05), 1.0)
        ax.set_ylabel("Pearson r")
        ax.set_title("Score agreement between sorts (each bar = one pair)")
        ax.set_xticks([])
        fig.tight_layout()
        p = out / "sort_concordance.png"
        fig.savefig(p, dpi=140); plt.close(fig); made.append(p)

    # --- bin profile of the strongest constructs ---
    top = [r for r in rows if r["n_sorts"] >= 2][:20]
    if top and per_sort:
        ranks = sorted({b for d in per_sort.values() for b in d["profile_ranks"]})
        mat = []
        for r in top:
            prof = []
            for rk in ranks:
                vals = [d["profile"].get(r["construct"], {}).get(rk)
                        for d in per_sort.values()]
                vals = [v for v in vals if v is not None]
                prof.append(sum(vals) / len(vals) if vals else 0.0)
            tot = sum(prof) or 1.0
            mat.append([v / tot for v in prof])
        fig, ax = plt.subplots(figsize=(1.1 * len(ranks) + 4, 0.34 * len(top) + 1.6))
        im = ax.imshow(mat, aspect="auto", cmap="magma")
        ax.set_xticks(range(len(ranks)))
        ax.set_xticklabels([bin_letters[i] if i < len(bin_letters) else str(rk)
                             for i, rk in enumerate(ranks)])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([f"{r['target'][:30]}" for r in top], fontsize=6)
        ax.set_xlabel("expression bin")
        ax.set_title("Bin profile of the 20 most significant constructs")
        fig.colorbar(im, ax=ax, label="share of construct's reads")
        fig.tight_layout()
        p = out / "bin_profiles.png"
        fig.savefig(p, dpi=140); plt.close(fig); made.append(p)
    return made


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_discover(args: argparse.Namespace) -> int:
    d = discover_screen(args.accession)
    if "error" in d:
        print(d["error"])
        return 2
    print(f"Screen:     {d['series']}")
    print(f"Assay:      {d['assay']}   Readout: {d['readout']}")
    print(f"Lab:        {d['lab']}")
    print(f"Library:    {', '.join(d['libraries']) or 'none linked'}")
    print(f"Bins:       {''.join(d['bin_letters'])}  "
          f"(A = lowest expression)")
    print(f"Sorts:      {len(d['sorts'])}  "
          f"[(bio_rep, flow_rep): {', '.join(str(s) for s in d['sorts'])}]")
    print(f"\n{len(d['bins'])} measurement sets:")
    for b in d["bins"]:
        mark = "  <- queried" if b["accession"] == args.accession else ""
        rep = (f"bio{b['bio_rep']}/" if b["bio_rep"] is not None else "")
        print(f"  {b['accession']}  {rep}flow{b['flow_rep']}  "
              f"Bin{b['bin']}{mark}")
    exp = len(d["sorts"]) * len(d["bin_letters"])
    if len(d["bins"]) != exp:
        print(f"\nNOTE: {len(d['bins'])} sets for {len(d['sorts'])} sorts x "
              f"{len(d['bin_letters'])} bins = {exp} expected. Some sorts are "
              f"missing bins; scoring uses whichever bins each sort has.")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    setup_logging()
    d = discover_screen(args.accession)
    if "error" in d:
        print(d["error"])
        return 2
    bins = d["bins"]
    print(f"Screen:     {d['series']}  ({len(bins)} sets, "
          f"{len(d['sorts'])} sorts, bins {''.join(d['bin_letters'])})")
    print(f"Readout:    {d['readout']}")

    lib = rp.find_guide_library(args.accession)
    if not lib["resolved"]:
        print(f"No construct library: {lib['why']}")
        return 2
    gf = lib["guide_files"][0]
    idx = rp.load_guide_index(gf["id"] or gf["accession"])
    if not idx["resolved"]:
        print(f"Construct table {gf['accession']}: {idx['why']}")
        return 2
    print(f"Library:    {gf['accession']}  ({idx['n_rows']:,} constructs)")

    cal = rp.calibrate_key(idx, args.accession, args.calibrate_reads)
    if cal["tested"]:
        print(f"Key calibration on {cal['n_reads']:,} reads "
              f"({cal['read_len']} bp):")
        for t in sorted(cal["tested"], key=lambda t: -t["rate"]):
            mark = ("  <- used" if cal["chosen"]
                    and t["column"] == cal["chosen"]["column"]
                    and (not cal["read_type"]
                         or t["read_type"] in cal["read_type"]) else "")
            print(f"    {t['read_type']:4} {t['column']:24} "
                  f"in reads {t['rate']:>6.1%}  ambig {t['ambiguous']:>5.1%}"
                  f"{mark}")
    if not cal["chosen"]:
        print("\nNo candidate sequence column appears in these reads. Either "
              "they are not the construct amplicon, or the library table "
              "describes a different assay. Not guessing.")
        return 3
    key = cal["chosen"]
    matcher = rp.build_matcher(key["seq_to_guide"])
    print(f"Counting by: {key['column']}  ({key['distinct']:,} distinct)"
          f"   mates: {cal['read_type'] or 'all'}")

    sorts = d["sorts"][:args.max_sorts] if args.max_sorts else d["sorts"]
    use = [b for b in bins if (b["bio_rep"], b["flow_rep"]) in set(sorts)]
    if len(use) < len(bins):
        print(f"Counting:    {len(use)} of {len(bins)} sets "
              f"(--max-sorts {args.max_sorts})")

    per_bin: "dict[str, Counter]" = {}
    qc: "list[dict]" = []
    for b in use:
        c, scanned, assigned, cached = rp.count_guides(
            b["accession"], matcher, key["column"], args.max_reads,
            reuse=not args.recount, read_types=cal.get("read_type"))
        per_bin[b["accession"]] = c
        rate = assigned / scanned if scanned else 0.0
        qc.append({"accession": b["accession"], "bin": b["bin"],
                    "rank": b["rank"], "bio_rep": b["bio_rep"],
                    "flow_rep": b["flow_rep"], "reads_scanned": scanned,
                    "reads_assigned": assigned,
                    "assignment_rate": round(rate, 4),
                    "constructs_seen": len(c),
                    "library_coverage": round(len(c) / max(key["distinct"], 1), 4),
                    "from_cache": cached})
        rep = (f"bio{b['bio_rep']}/" if b["bio_rep"] is not None else "")
        print(f"  {rep}flow{b['flow_rep']} Bin{b['bin']}  {scanned:>8,} reads, "
              f"{assigned:>8,} assigned ({rate:>5.1%}), {len(c):>5} constructs"
              f"{'  (cached)' if cached else ''}")

    mean_rate = (sum(q["assignment_rate"] for q in qc) / len(qc)) if qc else 0.0
    if mean_rate < 0.02:
        print(f"\n  Only {mean_rate:.1%} of reads carry a known "
              f"{key['column']}. Too few to score; refusing to report "
              f"numbers built on it.")
        return 3

    # One set of scores per sort.
    per_sort: "dict[tuple, dict]" = {}
    for s in sorts:
        sb = [b for b in use if (b["bio_rep"], b["flow_rep"]) == s]
        if len(sb) < 2:
            print(f"  NOTE: sort {s} has {len(sb)} bin(s); a mean-bin score "
                  f"needs at least 2. Excluded.")
            continue
        score, cvar, reads = sort_scores(per_bin, sb, args.min_count)
        prof: "dict[str, dict[int, float]]" = defaultdict(dict)
        for b in sb:
            tot = max(sum(per_bin.get(b["accession"], Counter()).values()), 1)
            for g, n in per_bin.get(b["accession"], Counter()).items():
                prof[g][b["rank"]] = n / tot
        per_sort[s] = {"score": score, "cvar": cvar, "reads": reads,
                        "profile": prof,
                        "profile_ranks": sorted({b["rank"] for b in sb})}
    if not per_sort:
        print("\nNo sort had two or more bins; nothing to score.")
        return 3

    rows, mod = score_gradient(per_sort, idx["target"], idx["type"])
    conc = concordance(per_sort)

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_{d['series'][:40]}"
    out = OUT_DIR / label
    out.mkdir(parents=True, exist_ok=True)
    cols = ["construct", "target", "construct_type", "n_sorts", "reads",
            "delta_bins", "sd", "sd_counting", "sd_used", "t", "df",
            "p_value", "fdr"]
    with (out / "construct_effects.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with (out / "qc_per_library.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc[0]), delimiter="\t")
        w.writeheader()
        for q in qc:
            w.writerow(q)

    ctrls = [r for r in rows
             if any(w in (r["construct_type"] or "").lower()
                     for w in ("non-targeting", "nontargeting", "control"))]
    rs = [c["r"] for c in conc if c["r"] is not None]
    summary = {
        "screen": d["series"], "query_set": args.accession,
        "assay": d["assay"], "readout": d["readout"],
        "sets_in_screen": len(bins), "sets_counted": len(use),
        "sorts_scored": len(per_sort),
        "bins": d["bin_letters"],
        "library": gf["accession"],
        "counting_key": key["column"],
        "read_mates_counted": cal["read_type"] or "all",
        "mean_assignment_rate": round(mean_rate, 4),
        "constructs_scored": len(rows),
        "controls_scored": len(ctrls),
        "control_null_usable": len(ctrls) >= 20,
        "significant_fdr_0.05": sum(1 for r in rows if r["fdr"] < 0.05),
        "score": "frequency-weighted mean bin rank, A=1 .. "
                 f"{d['bin_letters'][-1] if d['bin_letters'] else '?'}"
                 f"={len(d['bin_letters'])}",
        "direction": "positive delta_bins = sorted into HIGHER expression bins",
        "test": "two-sided Student t on per-sort centred scores, df = k-1",
        "sd_floor": mod["sd_floor"],
        "centring": mod["centring"],
        "sort_concordance_r": {"n_pairs": len(rs),
                                "median": round(_percentile(sorted(rs), 0.5), 4)
                                if rs else None,
                                "min": round(min(rs), 4) if rs else None},
        "min_count": args.min_count,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    (out / "sort_concordance.json").write_text(json.dumps(conc, indent=2))
    plots = make_plots(out, rows, qc, per_sort, conc, d["bin_letters"],
                        d["series"])

    print()
    for k, v in summary.items():
        if k in ("centring", "sort_concordance_r"):
            print(f"  {k}: {json.dumps(v)[:150]}")
        else:
            print(f"  {k}: {v}")
    if len(ctrls) < 20:
        print(f"\n  CAVEAT: {len(ctrls)} control construct(s) scored, too few "
              f"for an empirical null, and each sort was centred on the "
              f"library median instead where controls were short. The FDR "
              f"then measures agreement between sorts, not a calibrated "
              f"false-discovery rate.")
    if rs and min(rs) < 0.2:
        print(f"\n  CAVEAT: the worst sort-to-sort agreement is r={min(rs):.2f}. "
              f"Sorts that disagree about which constructs sit high cannot "
              f"support per-construct claims; see sort_concordance.json.")
    print(f"\nTop 15 by significance:")
    print(f"  {'target':38} {'dbins':>7} {'k':>3} {'reads':>8} {'fdr':>9}")
    for r in rows[:15]:
        print(f"  {r['target'][:38]:38} {r['delta_bins']:>7.3f} "
              f"{r['n_sorts']:>3} {r['reads']:>8,} {r['fdr']:>9.2e}")
    print(f"\nEffects: {out / 'construct_effects.tsv'}")
    print(f"QC:      {out / 'qc_per_library.tsv'}")
    for p in plots:
        print(f"Plot: {p}")
    print(f"Output: {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gradient_screen_analysis",
        description="FACS screens sorted into lettered expression bins.")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("discover", help="Find every bin and sort of the screen.")
    d.add_argument("accession")
    a = sub.add_parser("analyze", help="Count constructs, score across bins.")
    a.add_argument("accession")
    a.add_argument("--min-count", type=int, default=20,
                    help="Minimum total reads for a construct within a sort.")
    a.add_argument("--max-reads", type=int, default=None,
                    help="Cap reads scanned per library.")
    a.add_argument("--max-sorts", type=int, default=None,
                    help="Score only the first N sorts (for a quick look).")
    a.add_argument("--calibrate-reads", type=int, default=20000)
    a.add_argument("--recount", action="store_true",
                    help="Ignore cached counts and rescan the FASTQs.")
    a.add_argument("--label")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    return {"discover": cmd_discover, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
