#!/usr/bin/env python3
"""CRISPR FACS screens: reprocess a whole screen from raw reads.

A single sorted bin is not analysable on its own. IGVFDS6464SOVZ is
"18loci_uptake_Rep1_bottom20" -- one tail of one replicate of a base-editing
screen in HepG2 sorted on LDL-C uptake. The measurement IS the comparison
between bins, so counting guides in that one library gives library
composition and nothing about biology. Asked to "analyse IGVFDS6464SOVZ",
the honest move is to find its 15 siblings and analyse the screen.

Two things about these screens refuse to be assumed, and both are settled by
measurement rather than convention:

  Bin names. The same lab writes both "_Rep1_bottom20_ms" and
  "_Rep1_Bot20_ms", and some screens add an unsorted "_Bulk_ms".

  What identifies a construct. A CRISPR-KO library has one unique spacer per
  guide. A prime-editing library does not -- the LDLR library's 1,741 pegRNAs
  share 52 spacers, because the spacer only sets the nick site and the variant
  lives in the RT template. Counting that by spacer yields effect sizes and
  FDRs for variants whose reads were never told apart. So the counting key is
  chosen by testing candidate columns against real reads; see
  raw_data_pipeline.load_guide_index and calibrate_key.

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
import hashlib
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
_LN2 = math.log(2.0)
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

COUNT_CACHE = rp.FASTQ_CACHE.parent / "_counts"


def _count_cache_path(accession: str, matcher: dict, key: str,
                       max_reads: Optional[int],
                       read_types: "Optional[list[str]]" = None) -> Path:
    """Where the counts for exactly this library-and-key combination live.

    The library itself is part of the key: a construct table that gains or
    loses a sequence must not silently reuse counts made against the old one.
    """
    sig = hashlib.sha256()
    sig.update(f"{key}|{max_reads}|{','.join(read_types or ['*'])}|".encode())
    for sq in sorted(matcher["pref"]):
        sig.update(sq.encode())
    for v in matcher["pref"].values():
        sig.update(str(len(v)).encode())
    return COUNT_CACHE / f"{accession}.{key}.{sig.hexdigest()[:16]}.json"


def count_guides(accession: str, matcher: dict, key: str,
                  max_reads: Optional[int], reuse: bool = True,
                  read_types: "Optional[list[str]]" = None
                  ) -> "tuple[Counter, int, int, bool]":
    """Reads per construct for one library.

    Returns (counts, reads_scanned, reads_assigned, from_cache). The scanned
    and assigned totals let the caller report an assignment rate: a rate near
    zero means the key is wrong for this library -- worth seeing, rather than
    left to be inferred from an all-zero table.

    Counts are cached per (library, counting key, read cap, construct set).
    Scanning a screen's eight libraries takes about twenty minutes, and every
    re-analysis at a different --tail or --min-count was repeating all of it
    to reach numbers that had not changed. FASTQ downloads were already
    cached; the counting was not.
    """
    cache = _count_cache_path(accession, matcher, key, max_reads, read_types)
    if reuse and cache.exists():
        try:
            d = json.loads(cache.read_text())
            return (Counter(d["counts"]), int(d["scanned"]),
                    int(d["assigned"]), True)
        except (ValueError, KeyError, OSError):
            # A truncated or hand-edited cache file is not worth a crash, and
            # not worth trusting either.
            logging.getLogger(__name__).warning(
                "ignoring unreadable count cache %s", cache)

    comp = str.maketrans("ACGTN", "TGCAN")
    counts: "Counter[str]" = Counter()
    scanned = assigned = 0
    rp.FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
    for f in rp.biological_fastqs(accession, read_types):
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
    try:
        COUNT_CACHE.mkdir(parents=True, exist_ok=True)
        tmp = cache.with_suffix(".json.part")
        tmp.write_text(json.dumps({"accession": accession, "key": key,
                                    "max_reads": max_reads,
                                    "scanned": scanned, "assigned": assigned,
                                    "counts": dict(counts)}))
        tmp.replace(cache)      # atomic: never leave a half-written cache
    except OSError as e:
        logging.getLogger(__name__).warning("could not cache counts: %s", e)
    return counts, scanned, assigned, False


# ─── Statistics ─────────────────────────────────────────────────────────────

def _norm_sf(z: float) -> float:
    """Two-sided normal tail probability, via erfc. No scipy needed."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def _t_sf(t: float, df: int) -> float:
    """Two-sided Student-t tail probability.

    The normal is the wrong reference here and it is not a subtle error. A
    variant with two replicate observations that happen to agree to within
    0.016 log2 units gets se = 0.011 and z = -215, which the normal reports
    as p = 0 exactly -- infinite confidence from two numbers. The t
    distribution with df = k-1 is the textbook correction: at df = 1 it is
    Cauchy, whose tails are heavy enough that |t| = 215 is p ~ 3e-3 rather
    than 0.

    scipy is present in both the local environment and the deployed image,
    but the fallback matters: without it this would silently revert to the
    normal and to implausible p-values, so it degrades to a documented
    df-scaled approximation instead of pretending nothing changed.
    """
    if df < 1 or not math.isfinite(t):
        return 1.0
    try:
        from scipy import stats
        return float(2.0 * stats.t.sf(abs(t), df))
    except ImportError:
        pass
    if df == 1:                       # Cauchy, exactly
        return 2.0 * (0.5 - math.atan(abs(t)) / math.pi)
    if df == 2:                       # closed form
        return 1.0 - abs(t) / math.sqrt(2.0 + t * t)
    # Otherwise the normal on a variance-inflated statistic. Approximate,
    # and conservative in the direction that matters: it does not turn a
    # small sample into certainty.
    return _norm_sf(abs(t) / math.sqrt(df / max(df - 2.0, 1.0)))


def _percentile(xs: "list[float]", q: float) -> float:
    """Linear-interpolated percentile. Small helper, avoids a numpy import
    in a function that otherwise only needs the standard library."""
    if not xs:
        return 0.0
    ys = sorted(xs)
    if len(ys) == 1:
        return ys[0]
    i = (len(ys) - 1) * q
    lo, hi = int(math.floor(i)), int(math.ceil(i))
    return ys[lo] + (ys[hi] - ys[lo]) * (i - lo)


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
                  low_pct: int, min_count: int) -> "tuple[list[dict], dict]":
    """Per-variant effect from construct enrichment, low tail vs high tail.

    Returns (rows, moderation) where moderation records the variance floor
    applied and how many targets it was estimated from.

    One log2 ratio per replicate, then combined across replicates. Combining
    replicate ratios rather than pooling raw counts is deliberate: pooling
    lets the deepest-sequenced replicate dominate, and the replicate spread
    is the only estimate of variability available for a p-value.

    That spread is a weak estimate at four replicates and a bad one at two,
    so the p-value needs two corrections that the first version lacked. It
    used a normal reference and no variance floor, and reported the top hit
    of this screen -- two observations agreeing to within 0.016 log2 units --
    at p = 0 exactly, with 41 targets under FDR 0.05. Both are artefacts of
    treating a two-point standard error as if it were known. The test is now
    Student t with df = k-1, and each sd is floored at the 10th percentile of
    the sds seen across targets with k >= 3.
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
    # guide -> [(log2 ratio, counting-noise variance of that ratio), ...]
    ratios: "dict[str, list[tuple[float, float]]]" = defaultdict(list)
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
            lr = math.log2(((a + 0.5) / tlo) / ((b_ + 0.5) / thi))
            # Counting noise alone puts a floor under how precise this ratio
            # can be. By the delta method, a log2 count has variance
            # 1/(n ln2^2), so the log2 ratio of two counts has variance
            # (1/a + 1/b)/ln2^2. A variant observed 50 times cannot be
            # measured to 0.001 no matter how well its replicates agree.
            pvar = (1.0 / (a + 0.5) + 1.0 / (b_ + 0.5)) / (_LN2 ** 2)
            ratios[g].append((lr, pvar))

    # Aggregate guides onto the variant each installs.
    by_target: "dict[str, list[tuple[float, float]]]" = defaultdict(list)
    tgt_type: "dict[str, str]" = {}
    for g, obs in ratios.items():
        t = guide_target.get(g)
        if not t:
            continue
        by_target[t].extend(obs)
        tgt_type.setdefault(t, guide_type.get(g, ""))

    stats_by_target = {}
    for t, obs in by_target.items():
        k = len(obs)
        vals = [v for v, _pv in obs]
        mean = sum(vals) / k
        if k > 1:
            var = sum((v - mean) ** 2 for v in vals) / (k - 1)
            sd = math.sqrt(var)
        else:
            sd = float("nan")
        # The counting-noise sd this target's own read depths imply for the
        # mean of k observations, converted back to an sd on the same scale
        # as the empirical sd so the two can be compared directly.
        sd_count = math.sqrt(sum(pv for _v, pv in obs) / k)
        stats_by_target[t] = (k, mean, sd, sd_count)

    # Variance moderation. Even under the t distribution, a target whose few
    # observations happen to agree almost exactly gets an sd near zero and a
    # t statistic set by luck rather than by effect size. Floor each sd at a
    # low percentile of the sds actually observed across the screen, taken
    # from targets with k >= 3 because their sd is the better estimate. This
    # is the idea behind limma's variance moderation, at its simplest: the
    # screen's own spread is a prior on how quiet a target can plausibly be.
    ref_sds = [sd for k, _m, sd, _sc in stats_by_target.values()
               if k >= 3 and not math.isnan(sd) and sd > 0]
    sd_floor = _percentile(ref_sds, 0.10) if ref_sds else 0.0

    rows = []
    for t, (k, mean, sd, sd_count) in stats_by_target.items():
        if k > 1:
            # Three lower bounds, whichever binds hardest. The empirical
            # floor needs k >= 3 targets to exist at all, and a screen with
            # only two replicates has none -- there the counting-noise floor
            # is the whole protection, which is why it is not optional.
            sd_used = max(sd, sd_floor, sd_count)
            se = sd_used / math.sqrt(k)
            tstat = mean / se if se > 0 else 0.0
            # df = k-1, the honest degrees of freedom for k observations.
            pval = _t_sf(tstat, k - 1) if se > 0 else 1.0
        else:
            sd_used = se = float("nan")
            tstat, pval = 0.0, 1.0
        rows.append({"target": t, "target_type": tgt_type.get(t, ""),
                      "n_guide_obs": k, "mean_log2_low_over_high": round(mean, 4),
                      "sd": None if math.isnan(sd) else round(sd, 4),
                      "sd_counting": round(sd_count, 4),
                      "sd_used": None if math.isnan(sd_used) else round(sd_used, 4),
                      "t": round(tstat, 4), "df": max(k - 1, 0),
                      "p_value": pval})
    ps = [r["p_value"] for r in rows]
    for r, q in zip(rows, benjamini_hochberg(ps)):
        r["fdr"] = q
    # Rank by p, then by effect size, so ties among the many k=4 targets that
    # share a p-value are not ordered by dictionary insertion.
    rows.sort(key=lambda r: (r["p_value"],
                              -abs(r["mean_log2_low_over_high"])))
    return rows, {"sd_floor": round(sd_floor, 4), "n_ref_sd": len(ref_sds)}


# ─── Plots ──────────────────────────────────────────────────────────────────

def is_control(target_type: "Optional[str]") -> bool:
    """Whether a construct is a control, across the spellings in use.

    This metadata says "non-targeting"; other libraries say "control". A
    predicate matching only "control" found none of the LDLR library's, so
    the volcano drew no control points and the summary reported none.
    """
    t = (target_type or "").lower()
    return any(w in t for w in ("control", "non-targeting", "nontargeting"))


def make_plots(out: Path, rows: "list[dict]", per_bin: "dict[str, Counter]",
                bins: "list[dict]", low_pct: int,
                screen: str = "") -> "list[Path]":
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
    is_ctrl = np.array([is_control(r["target_type"]) for r in rows])
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
    # The sorted phenotype is whatever this screen sorted on, which the tool
    # cannot know -- it was hardcoded to "LDL-C uptake" from the first screen
    # it was written against, and would have mislabelled every other one.
    ax[0].set_title(f"Variant effects{' — ' + screen if screen else ''}\n"
                    f"positive = enriched in the LOW bin = variant reduces "
                    f"the sorted phenotype")
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
        c, scanned, assigned, cached = count_guides(
            b["accession"], matcher, key["column"], args.max_reads,
            reuse=not args.recount, read_types=cal.get("read_type"))
        per_bin[b["accession"]] = c
        rate = assigned / scanned if scanned else 0.0
        rates.append(rate)
        print(f"  Rep{b['rep']} {_bin_label(b):<16}  {scanned:>9,} reads, "
              f"{assigned:>8,} assigned ({rate:>5.1%}), "
              f"{len(c):>5} of {key['distinct']:,} constructs seen"
              f"{'  (cached)' if cached else ''}")
    mean_rate = sum(rates) / len(rates) if rates else 0.0
    if mean_rate < 0.02:
        print(f"\n  Only {mean_rate:.1%} of reads carry a known "
              f"{key['column']}. That is too few to score: the key is probably "
              f"wrong for this library, or these FASTQs are not the guide "
              f"amplicon. Refusing to report numbers built on it.")
        return 3

    rows, moderation = score_screen(per_bin, counted, meta["target"],
                                     meta["type"], args.tail, args.min_count)
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
            "mean_log2_low_over_high", "sd", "sd_counting", "sd_used",
            "t", "df", "p_value", "fdr"]
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
    # Both spellings: this metadata uses "non-targeting", others "control".
    ctrls = [r for r in rows if is_control(r["target_type"])]
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
        "control_null_usable": len(ctrls) >= 20,
        "controls_significant": sum(1 for r in ctrls if r["fdr"] < 0.05),
        "min_count": args.min_count,
        "test": "two-sided Student t on per-replicate log2 ratios, df = k-1",
        "sd_floor": moderation["sd_floor"],
        "sd_floor_from_n_targets": moderation["n_ref_sd"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    plots = make_plots(out, rows, per_bin, counted, args.tail,
                        screen=scr['series'])

    print()
    for k, v in summary.items():
        print(f"  {k}: {v}")
    # A screen with a handful of non-targeting constructs cannot support an
    # empirical null, and saying so matters: the FDR below then rests entirely
    # on the spread between replicates. The LDLR library carries exactly one
    # non-targeting construct in 1,741, which reads as a control column being
    # present without being usable.
    if len(ctrls) < 20:
        print(f"\n  CAVEAT: only {len(ctrls)} non-targeting control(s) were "
              f"scored, too few for an empirical null. The FDR above comes "
              f"from the spread across replicates alone, so it measures "
              f"consistency between replicates, not the false-discovery rate "
              f"you would get by calibrating against constructs known to do "
              f"nothing. Treat the ranking as more trustworthy than the "
              f"absolute q-values.")
    print(f"\nTop 10 by significance:")
    print(f"  {'target':34} {'log2':>7} {'n':>3} {'sd':>7} {'t':>8} {'fdr':>9}")
    for r in rows[:10]:
        print(f"  {r['target'][:34]:34} "
              f"{r['mean_log2_low_over_high']:>7.2f} {r['n_guide_obs']:>3} "
              f"{(r['sd_used'] if r['sd_used'] is not None else float('nan')):>7.3f} "
              f"{r['t']:>8.1f} {r['fdr']:>9.2e}")
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
    a.add_argument("--recount", action="store_true",
                    help="Ignore cached counts and rescan the FASTQs.")
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
