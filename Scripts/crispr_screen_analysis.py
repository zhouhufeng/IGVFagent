#!/usr/bin/env python3
"""CRISPR FACS screens: reprocess a whole screen from raw reads.

A single sorted bin is not analysable on its own. IGVFDS6464SOVZ is
"18loci_uptake_Rep1_bottom20" -- one tail of one replicate of a base-editing
screen in HepG2 sorted on LDL-C uptake. The measurement IS the comparison
between bins, so counting guides in that one library gives library
composition and nothing about biology. Asked to "analyse IGVFDS6464SOVZ",
the honest move is to find its 15 siblings and analyse the screen.

The pipeline mirrors what the lab's own AnalysisSet publishes, which is the
canonical shape for this assay:

  1. guide quantifications        reads per guide, per bin
  2. differential quantifications enrichment between the sorted tails
  3. variant effects              guides aggregated onto the variant each
                                  one installs, combined across replicates

Direction matters and is stated everywhere it is computed: a guide enriched
in the LOW-uptake tail installs a variant that REDUCES uptake, so the score
is log2(low / high) and negative means the variant raises uptake.

Discovery works off the submitter's alias
(`richard-sherwood:18loci_uptake_Rep1_bottom20_ms`), which encodes series,
replicate and bin. That is a convention, not a schema, so `screen discover`
prints what it matched for checking before any compute is spent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import raw_data_pipeline as rp                                # noqa: E402

ROOT = rp.ROOT
OUT_DIR = ROOT / "Docs" / "CrisprScreen"
LOG_DIR = ROOT / "Docs" / "Logs"

# Bin naming is NOT consistent, even within one lab. Both of these are
# Sherwood CRISPR FACS screens:
#   richard-sherwood:18loci_uptake_Rep1_bottom20_ms
#   richard-sherwood:0426_LDLR137-219_repo_Rep1_Bot20_ms
# A pattern accepting only "bottom" matched every Top bin of the LDLR screen
# and none of its Bot bins, so discovery returned 8 of 20 libraries, no
# replicate had a complete pair, and the run produced nothing. Hence: both
# spellings, plus low/high, plus the unsorted Bulk bin that screen also has.
_ALIAS = re.compile(
    r"(?P<series>[A-Za-z0-9_.\-]+?)_Rep(?P<rep>\d+)_"
    r"(?:(?P<side>bottom|bot|low|top|high)(?P<pct>\d+)|(?P<bulk>bulk))",
    re.I)

# Normalise the spellings onto one internal vocabulary.
_SIDE_ALIASES = {"bottom": "bottom", "bot": "bottom", "low": "bottom",
                 "top": "top", "high": "top"}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"crispr_screen_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(log),
                                  logging.StreamHandler(sys.stdout)])
    return log


def _bin_label(b: dict) -> str:
    """Human label for a bin. `bulk` has no percentage, so do not print one."""
    side = b["side"] if isinstance(b, dict) else b
    pct = b.get("pct", 0) if isinstance(b, dict) else 0
    return "bulk (unsorted)" if side == "bulk" else f"{side}{pct}%"


def _alias_of(fs: dict) -> str:
    return " ".join(fs.get("aliases") or [])


def parse_alias(alias: str) -> "Optional[dict]":
    m = _ALIAS.search(alias)
    if not m:
        return None
    series = m.group("series").split(":")[-1]
    if m.group("bulk"):
        # The unsorted population. Not a tail, but the natural denominator
        # when a screen provides one -- and it must not be mistaken for one.
        return {"series": series, "rep": int(m.group("rep")),
                "side": "bulk", "pct": 0}
    return {"series": series, "rep": int(m.group("rep")),
            "side": _SIDE_ALIASES[m.group("side").lower()],
            "pct": int(m.group("pct"))}


def discover_screen(accession: str) -> dict:
    """Every sorted bin belonging to the same screen as `accession`."""
    st, own = rp.portal_json(
        f"/search/?type=FileSet&accession={accession}&format=json")
    rows = (own or {}).get("@graph") or []
    if not rows:
        return {"error": f"{accession} not found"}
    st, full = rp.portal_json(f"{rows[0].get('@id')}?format=json")
    meta = parse_alias(_alias_of(full or {}))
    if not meta:
        return {"error": (f"{accession} has no Rep/bin alias, so its screen "
                          f"cannot be identified. Alias was: "
                          f"{_alias_of(full or {}) or 'none'}")}
    lab = ((full or {}).get("lab") or {}).get("title", "")
    assay = ((full or {}).get("preferred_assay_titles") or [""])[0]
    # limit=all, not a page. With limit=200 against 300 CRISPR FACS screens
    # this returned 10 of the screen's 16 bins and, worse, silently omitted
    # the top20 partner of the very set being queried -- so a replicate
    # would have been dropped with nothing to show it had been.
    st, d = rp.portal_json(
        f"/search/?type=MeasurementSet&limit=all&format=json"
        f"&preferred_assay_titles={assay.replace(' ', '+')}")
    bins = []
    for r in (d or {}).get("@graph") or []:
        m = parse_alias(_alias_of(r))
        if m and m["series"] == meta["series"]:
            bins.append({"accession": r.get("accession"),
                          "alias": _alias_of(r), **m})
    bins.sort(key=lambda b: (b["rep"], b["side"], b["pct"]))
    return {"series": meta["series"], "lab": lab, "assay": assay,
            "query_set": accession, "bins": bins,
            "replicates": sorted({b["rep"] for b in bins}),
            "bin_kinds": sorted(
                ({"side": s, "pct": p} for s, p in
                 {(b["side"], b["pct"]) for b in bins}),
                key=lambda d: (d["side"], d["pct"])),
            "sides": sorted({f"{b['side']}{b['pct']}" for b in bins})}


# ─── Guide counting ─────────────────────────────────────────────────────────

def count_guides(accession: str, matcher: dict,
                  max_reads: Optional[int]) -> "tuple[Counter, int, int]":
    """Reads per construct for one library.

    Returns (counts, reads_scanned, reads_assigned) so the caller can report
    the assignment rate. A rate near zero means the key is wrong for this
    library -- worth seeing, rather than inferring from an all-zero table.

    `lengths` must be longest-first: RT templates nest (a 9 bp template is a
    prefix of longer ones), and taking the first match found would hand the
    read to the least specific construct.
    """
    comp = str.maketrans("ACGTN", "TGCAN")
    counts: "Counter[str]" = Counter()
    scanned = assigned = 0
    rp.FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
    for f in rp.list_files(accession):
        if str(f.get("file_format", "")).lower() != "fastq":
            continue
        name = Path(str(f.get("href") or f.get("accession"))).name
        dest = rp.FASTQ_CACHE / name
        if not dest.exists():
            rp.portal_download(f["href"], dest)
        op = gzip.open if dest.read_bytes()[:2] == b"\x1f\x8b" else open
        with op(dest, "rt", errors="replace") as fh:        # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                if max_reads and scanned >= max_reads:
                    break
                scanned += 1
                seq = line.strip().upper()
                hit = None
                for s_ in (seq, seq.translate(comp)[::-1]):
                    hit = rp.match_read(s_, matcher)
                    if hit:
                        break
                if hit:
                    counts[hit] += 1
                    assigned += 1
    return counts, scanned, assigned


# ─── Statistics ─────────────────────────────────────────────────────────────

def _norm_sf(z: float) -> float:
    """Two-sided normal tail probability, via erfc. No scipy needed."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def benjamini_hochberg(pvals: "list[float]") -> "list[float]":
    """BH-adjusted p-values, order preserved."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        val = min(prev, pvals[idx] * n / i)
        adj[idx] = val
        prev = val
    return adj


def score_screen(per_bin: "dict[str, Counter]", bins: "list[dict]",
                  guide_target: "dict[str, str]", guide_type: "dict[str, str]",
                  low_pct: int, min_count: int) -> "list[dict]":
    """Per-variant effect from guide enrichment, low tail vs high tail.

    One log2 ratio per replicate, then combined across replicates. Combining
    replicate ratios rather than pooling raw counts is deliberate: pooling
    lets the deepest-sequenced replicate dominate, and the replicate spread
    is the only estimate of variability available for a p-value.
    """
    reps = sorted({b["rep"] for b in bins})
    usable, skipped = [], []
    for rep in reps:
        has_lo = any(b["rep"] == rep and b["side"] == "bottom"
                     and b["pct"] == low_pct for b in bins)
        has_hi = any(b["rep"] == rep and b["side"] == "top"
                     and b["pct"] == low_pct for b in bins)
        (usable if (has_lo and has_hi) else skipped).append(rep)
    if skipped:
        print(f"  NOTE: replicate(s) {skipped} lack a complete "
              f"bottom{low_pct}%/top{low_pct}% pair and are excluded. "
              f"Scoring uses replicates {usable}.")
    # guide -> [log2 ratio per replicate]
    ratios: "dict[str, list[float]]" = defaultdict(list)
    for rep in usable:
        lo = next((b for b in bins if b["rep"] == rep and b["side"] == "bottom"
                    and b["pct"] == low_pct), None)
        hi = next((b for b in bins if b["rep"] == rep and b["side"] == "top"
                    and b["pct"] == low_pct), None)
        if not lo or not hi:
            continue
        clo, chi = per_bin.get(lo["accession"]), per_bin.get(hi["accession"])
        if not clo or not chi:
            continue
        tlo, thi = max(sum(clo.values()), 1), max(sum(chi.values()), 1)
        for g in set(clo) | set(chi):
            a, b_ = clo.get(g, 0), chi.get(g, 0)
            if a + b_ < min_count:
                continue
            # +0.5 so a guide absent from one tail still yields a finite,
            # bounded ratio instead of being dropped or becoming infinite.
            ratios[g].append(math.log2(((a + 0.5) / tlo) / ((b_ + 0.5) / thi)))

    # Aggregate guides onto the variant each installs.
    by_target: "dict[str, list[float]]" = defaultdict(list)
    tgt_type: "dict[str, str]" = {}
    for g, vals in ratios.items():
        t = guide_target.get(g)
        if not t:
            continue
        by_target[t].extend(vals)
        tgt_type.setdefault(t, guide_type.get(g, ""))

    rows = []
    for t, vals in by_target.items():
        k = len(vals)
        mean = sum(vals) / k
        if k > 1:
            var = sum((v - mean) ** 2 for v in vals) / (k - 1)
            sd = math.sqrt(var)
            se = sd / math.sqrt(k)
        else:
            sd = se = float("nan")
        z = mean / se if se and se > 0 else 0.0
        rows.append({"target": t, "target_type": tgt_type.get(t, ""),
                      "n_guide_obs": k, "mean_log2_low_over_high": round(mean, 4),
                      "sd": None if math.isnan(sd) else round(sd, 4),
                      "z": round(z, 4),
                      "p_value": _norm_sf(z) if se and se > 0 else 1.0})
    ps = [r["p_value"] for r in rows]
    for r, q in zip(rows, benjamini_hochberg(ps)):
        r["fdr"] = q
    rows.sort(key=lambda r: r["p_value"])
    return rows


# ─── Plots ──────────────────────────────────────────────────────────────────

def make_plots(out: Path, rows: "list[dict]", per_bin: "dict[str, Counter]",
                bins: "list[dict]", low_pct: int) -> "list[Path]":
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return []
    out.mkdir(parents=True, exist_ok=True)
    made = []

    eff = np.array([r["mean_log2_low_over_high"] for r in rows])
    fdr = np.array([max(r["fdr"], 1e-12) for r in rows])
    is_ctrl = np.array([("control" in (r["target_type"] or "").lower())
                        for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].scatter(eff[~is_ctrl], -np.log10(fdr[~is_ctrl]), s=14,
                  c="#4C72B0", label="variants", alpha=.7)
    if is_ctrl.any():
        ax[0].scatter(eff[is_ctrl], -np.log10(fdr[is_ctrl]), s=26,
                      c="#C44E52", label="controls", alpha=.9)
    ax[0].axhline(-math.log10(0.05), ls="--", lw=1, c="k")
    ax[0].axvline(0, ls=":", lw=1, c="grey")
    ax[0].set_xlabel(f"mean log2( bottom{low_pct}% / top{low_pct}% )")
    ax[0].set_ylabel("-log10 FDR")
    ax[0].set_title("Variant effects on LDL-C uptake\n"
                    "positive = enriched in LOW uptake = variant reduces uptake")
    ax[0].legend(fontsize=8)

    # Replicate agreement: the honest check on whether any of this is signal.
    reps = sorted({b["rep"] for b in bins})
    pairs = []
    for rep in reps:
        lo = next((b for b in bins if b["rep"] == rep
                    and b["side"] == "bottom" and b["pct"] == low_pct), None)
        hi = next((b for b in bins if b["rep"] == rep
                    and b["side"] == "top" and b["pct"] == low_pct), None)
        if lo and hi and per_bin.get(lo["accession"]) and per_bin.get(hi["accession"]):
            pairs.append((rep, per_bin[lo["accession"]], per_bin[hi["accession"]]))
    if len(pairs) >= 2:
        (r1, l1, h1), (r2, l2, h2) = pairs[0], pairs[1]
        common = (set(l1) | set(h1)) & (set(l2) | set(h2))
        t1l, t1h = max(sum(l1.values()), 1), max(sum(h1.values()), 1)
        t2l, t2h = max(sum(l2.values()), 1), max(sum(h2.values()), 1)
        x = [math.log2(((l1.get(g, 0) + .5) / t1l) / ((h1.get(g, 0) + .5) / t1h))
             for g in common]
        y = [math.log2(((l2.get(g, 0) + .5) / t2l) / ((h2.get(g, 0) + .5) / t2h))
             for g in common]
        ax[1].scatter(x, y, s=8, alpha=.4, c="#55A868")
        if len(x) > 2:
            r = float(np.corrcoef(x, y)[0, 1])
            ax[1].set_title(f"guide log2 ratio, Rep{r1} vs Rep{r2}  (r = {r:.2f})")
        ax[1].set_xlabel(f"Rep{r1}")
        ax[1].set_ylabel(f"Rep{r2}")
        ax[1].axhline(0, ls=":", lw=1, c="grey")
        ax[1].axvline(0, ls=":", lw=1, c="grey")
    fig.tight_layout()
    p = out / "screen_effects.png"
    fig.savefig(p, dpi=140); plt.close(fig); made.append(p)

    # Library depth per bin -- a bin that failed to sequence invalidates its
    # replicate, and that must be visible rather than buried in a summary.
    fig, ax = plt.subplots(figsize=(11, 3.6))
    labels = [f"R{b['rep']} " + ("bulk" if b["side"] == "bulk"
               else f"{b['side'][:3]}{b['pct']}") for b in bins]
    vals = [sum(per_bin.get(b["accession"], Counter()).values()) for b in bins]
    ax.bar(labels, vals, color="#8172B2")
    ax.set_ylabel("reads assigned to a guide")
    ax.set_title("library depth by sorted bin")
    ax.tick_params(axis="x", rotation=60)
    fig.tight_layout()
    p = out / "screen_depth.png"
    fig.savefig(p, dpi=140); plt.close(fig); made.append(p)
    return made


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_discover(args: argparse.Namespace) -> int:
    scr = discover_screen(args.accession)
    if "error" in scr:
        print(scr["error"])
        return 2
    print(f"Screen:     {scr['series']}")
    print(f"Assay:      {scr['assay']}   Lab: {scr['lab']}")
    print(f"Replicates: {scr['replicates']}")
    print(f"Bins:       {', '.join(_bin_label(b) for b in scr['bin_kinds'])}")
    print(f"\n{len(scr['bins'])} libraries:")
    for b in scr["bins"]:
        mark = "  <- queried" if b["accession"] == args.accession else ""
        print(f"  {b['accession']}  Rep{b['rep']} {_bin_label(b)}{mark}")
    return 0


def counted_probe(scr: dict, tail: int) -> str:
    """A library to calibrate the key against: prefer a bin the scoring uses."""
    for b in scr["bins"]:
        if b["side"] in ("bottom", "top") and b["pct"] == tail:
            return b["accession"]
    return scr["bins"][0]["accession"]


def cmd_analyze(args: argparse.Namespace) -> int:
    scr = discover_screen(args.accession)
    if "error" in scr:
        print(scr["error"])
        return 2
    bins = scr["bins"]
    print(f"Screen:     {scr['series']}  ({len(bins)} libraries, "
          f"replicates {scr['replicates']})")

    lib = rp.find_guide_library(args.accession)
    if not lib["resolved"]:
        print(f"No guide table: {lib['why']}")
        return 2
    gf = lib["guide_files"][0]
    idx = rp.load_guide_index(gf["id"] or gf["accession"])
    if not idx["resolved"]:
        print(f"Guide table {gf['accession']}: {idx['why']}")
        return 2
    meta = {"target": idx["target"], "type": idx["type"]}
    print(f"Guide table: {gf['accession']}  ({idx['n_rows']:,} constructs)")

    cal = rp.calibrate_key(idx, (args.key_from or counted_probe(scr, args.tail)),
                        args.calibrate_reads)
    if cal["why"]:
        print(f"  {cal['why']}")
    if cal["tested"]:
        print(f"Key calibration on {cal['n_reads']:,} reads "
              f"({cal['read_len']} bp):")
        print(f"    {'column':24} {'separates':>11} {'in reads':>9} {'ambig':>7}")
        for t in sorted(cal["tested"], key=lambda t: -t["rate"]):
            mark = ("  <- used" if cal["chosen"]
                    and t["column"] == cal["chosen"]["column"] else "")
            print(f"    {t['column']:24} {t['fraction']:>10.1%} "
                  f"{t['rate']:>9.1%} {t['ambiguous']:>7.1%}{mark}")
    if not cal["chosen"]:
        print("\n  No candidate sequence column appears in the reads at all. "
              "Either these FASTQs are not the construct amplicon, or the "
              "library table describes a different assay. Not guessing.")
        return 3
    key = cal["chosen"]
    matcher, lengths = rp.build_matcher(key["seq_to_guide"]), key["lengths"]
    print(f"Counting by: {key['column']}  ({key['distinct']:,} distinct "
          f"sequences, {min(lengths)}-{max(lengths)} bp)")
    if key["distinct"] < idx["n_rows"] * 0.9:
        print(f"  WARNING: {key['column']} does not separate all constructs; "
              f"{idx['n_rows'] - key['distinct']:,} share a sequence with "
              f"another and their reads cannot be told apart.")

    # Counting a bin means downloading its whole FASTQ and scanning every
    # read for a spacer, so count only the bins the comparison will use.
    # This screen has 20 libraries but a bottom20/top20 comparison needs 8;
    # counting the 40% tails and the bulk bins as well was 2.5x the download
    # and the CPU for numbers nothing then read.
    counted = bins if args.all_bins else [
        b for b in bins if b["side"] in ("bottom", "top") and b["pct"] == args.tail]
    if not counted:
        have = ", ".join(_bin_label(b) for b in scr["bin_kinds"])
        print(f"\nNo bin matches --tail {args.tail}. This screen has: {have}")
        return 3
    if len(counted) < len(bins):
        skip = len(bins) - len(counted)
        print(f"Counting:    {len(counted)} of {len(bins)} libraries "
              f"(the {skip} not used by a bottom{args.tail}%/top{args.tail}% "
              f"comparison are skipped; --all-bins counts them too)")

    per_bin: "dict[str, Counter]" = {}
    rates = []
    for b in counted:
        c, scanned, assigned = count_guides(
            b["accession"], matcher, args.max_reads)
        per_bin[b["accession"]] = c
        rate = assigned / scanned if scanned else 0.0
        rates.append(rate)
        print(f"  Rep{b['rep']} {_bin_label(b):<16}  {scanned:>9,} reads, "
              f"{assigned:>8,} assigned ({rate:>5.1%}), "
              f"{len(c):>5} of {key['distinct']:,} constructs seen")
    mean_rate = sum(rates) / len(rates) if rates else 0.0
    if mean_rate < 0.02:
        print(f"\n  Only {mean_rate:.1%} of reads carry a known "
              f"{key['column']}. That is too few to score: the key is probably "
              f"wrong for this library, or these FASTQs are not the guide "
              f"amplicon. Refusing to report numbers built on it.")
        return 3

    rows = score_screen(per_bin, counted, meta["target"], meta["type"],
                         args.tail, args.min_count)
    if not rows:
        have = ", ".join(_bin_label(b) for b in scr["bin_kinds"])
        print(f"\nNo target could be scored. Scoring needs both bottom"
              f"{args.tail}% and top{args.tail}% for at least one replicate; "
              f"this screen has: {have}")
        return 3

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_{scr['series']}"
    out = OUT_DIR / label
    out.mkdir(parents=True, exist_ok=True)
    cols = ["target", "target_type", "n_guide_obs",
            "mean_log2_low_over_high", "sd", "z", "p_value", "fdr"]
    with (out / "variant_effects.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with (out / "guide_counts_by_bin.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["guide_id", "target"] +
                    [f"R{b['rep']}_{b['side']}{b['pct']}" for b in counted])
        for gid in sorted(idx["target"]):
            row = [per_bin[b["accession"]].get(gid, 0) for b in counted]
            if any(row):
                w.writerow([gid, meta["target"].get(gid, "")] + row)

    sig = [r for r in rows if r["fdr"] < 0.05]
    ctrls = [r for r in rows if "control" in (r["target_type"] or "").lower()]
    summary = {
        "screen": scr["series"], "query_set": args.accession,
        "libraries_in_screen": len(bins), "libraries_counted": len(counted),
        "bins_present": [_bin_label(b) for b in scr["bin_kinds"]],
        "replicates": scr["replicates"],
        "tail_compared": f"bottom{args.tail}% vs top{args.tail}%",
        "direction": "positive score = enriched in LOW uptake = variant reduces uptake",
        "guide_table": gf["accession"],
        "counting_key": key["column"],
        "constructs_in_library": idx["n_rows"],
        "constructs_distinguishable": key["distinct"],
        "read_length": cal["read_len"],
        "mean_read_assignment_rate": round(mean_rate, 4),
        "targets_scored": len(rows),
        "significant_fdr_0.05": len(sig),
        "controls_scored": len(ctrls),
        "controls_significant": sum(1 for r in ctrls if r["fdr"] < 0.05),
        "min_count": args.min_count,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    plots = make_plots(out, rows, per_bin, counted, args.tail)

    print()
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\nTop 10 by significance:")
    print(f"  {'target':34} {'type':18} {'log2':>7} {'n':>4} {'fdr':>9}")
    for r in rows[:10]:
        print(f"  {r['target'][:34]:34} {(r['target_type'] or '')[:18]:18} "
              f"{r['mean_log2_low_over_high']:>7.2f} {r['n_guide_obs']:>4} "
              f"{r['fdr']:>9.2e}")
    print(f"\nVariant effects: {out / 'variant_effects.tsv'}")
    for p in plots:
        print(f"Plot: {p}")
    print(f"Output: {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="crispr_screen_analysis",
        description="Reprocess a CRISPR FACS screen from raw reads.")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("discover", help="Find every sorted bin of the screen.")
    d.add_argument("accession")
    a = sub.add_parser("analyze", help="Count, compare tails, score variants.")
    a.add_argument("accession")
    a.add_argument("--tail", type=int, default=20,
                    help="Which tail pair to compare (20 or 40).")
    a.add_argument("--min-count", type=int, default=10,
                    help="Minimum combined reads for a guide in a replicate.")
    a.add_argument("--max-reads", type=int, default=None)
    a.add_argument("--key-from", help="Library accession to calibrate the "
                    "counting key against (default: a bin being scored).")
    a.add_argument("--calibrate-reads", type=int, default=20000,
                    help="Reads sampled to choose the counting key.")
    a.add_argument("--all-bins", action="store_true",
                    help="Count every library, including bins the tail "
                         "comparison does not use (slower; for QC).")
    a.add_argument("--label")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"discover": cmd_discover, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
