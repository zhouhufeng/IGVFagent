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

# richard-sherwood:18loci_uptake_Rep1_bottom20_ms
_ALIAS = re.compile(
    r"(?P<series>[A-Za-z0-9_.\-]+?)_Rep(?P<rep>\d+)_(?P<side>bottom|top)"
    r"(?P<pct>\d+)", re.I)


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"crispr_screen_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(log),
                                  logging.StreamHandler(sys.stdout)])
    return log


def _alias_of(fs: dict) -> str:
    return " ".join(fs.get("aliases") or [])


def parse_alias(alias: str) -> "Optional[dict]":
    m = _ALIAS.search(alias)
    if not m:
        return None
    series = m.group("series").split(":")[-1]
    return {"series": series, "rep": int(m.group("rep")),
            "side": m.group("side").lower(), "pct": int(m.group("pct"))}


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
            "sides": sorted({f"{b['side']}{b['pct']}" for b in bins})}


# ─── Guide counting ─────────────────────────────────────────────────────────

def count_guides(accession: str, spacer_to_id: "dict[str, str]",
                  lengths: "list[int]", max_reads: Optional[int]) -> Counter:
    """Reads per guide for one library."""
    comp = str.maketrans("ACGTN", "TGCAN")
    counts: "Counter[str]" = Counter()
    n = 0
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
                if max_reads and n >= max_reads:
                    break
                n += 1
                seq = line.strip().upper()
                for s in (seq, seq.translate(comp)[::-1]):
                    hit = None
                    for L in lengths:
                        for off in range(0, len(s) - L + 1):
                            gid = spacer_to_id.get(s[off:off + L])
                            if gid:
                                hit = gid
                                break
                        if hit:
                            break
                    if hit:
                        counts[hit] += 1
                        break
    return counts


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
    labels = [f"R{b['rep']} {b['side'][:3]}{b['pct']}" for b in bins]
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
    print(f"Bins:       {', '.join(scr['sides'])}")
    print(f"\n{len(scr['bins'])} libraries:")
    for b in scr["bins"]:
        mark = "  <- queried" if b["accession"] == args.accession else ""
        print(f"  {b['accession']}  Rep{b['rep']} {b['side']}{b['pct']}%{mark}")
    return 0


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
    guides = rp.load_guide_table(gf["id"] or gf["accession"])
    if not guides:
        print(f"Guide table {gf['accession']} parsed to 0 guides.")
        return 2
    print(f"Guide table: {gf['accession']}  ({len(guides):,} guides)")

    # The table's own columns say what each guide targets and whether it is a
    # control -- far better than parsing the guide name.
    meta = _guide_metadata(gf["id"] or gf["accession"])
    spacer_to_id: "dict[str, str]" = {}
    for gid, sp in guides:
        spacer_to_id.setdefault(sp, gid)
    lengths = sorted({len(sp) for sp in spacer_to_id})

    per_bin: "dict[str, Counter]" = {}
    for b in bins:
        c = count_guides(b["accession"], spacer_to_id, lengths, args.max_reads)
        per_bin[b["accession"]] = c
        print(f"  Rep{b['rep']} {b['side']}{b['pct']:<3}  "
              f"{sum(c.values()):>8,} reads assigned, {len(c):>5} guides seen")

    rows = score_screen(per_bin, bins, meta["target"], meta["type"],
                         args.tail, args.min_count)
    if not rows:
        print(f"\nNo target could be scored. Are both bottom{args.tail}% and "
              f"top{args.tail}% present for at least one replicate?")
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
                    [f"R{b['rep']}_{b['side']}{b['pct']}" for b in bins])
        for gid, _sp in guides:
            row = [per_bin[b["accession"]].get(gid, 0) for b in bins]
            if any(row):
                w.writerow([gid, meta["target"].get(gid, "")] + row)

    sig = [r for r in rows if r["fdr"] < 0.05]
    ctrls = [r for r in rows if "control" in (r["target_type"] or "").lower()]
    summary = {
        "screen": scr["series"], "query_set": args.accession,
        "libraries": len(bins), "replicates": scr["replicates"],
        "tail_compared": f"bottom{args.tail}% vs top{args.tail}%",
        "direction": "positive score = enriched in LOW uptake = variant reduces uptake",
        "guide_table": gf["accession"],
        "targets_scored": len(rows),
        "significant_fdr_0.05": len(sig),
        "controls_scored": len(ctrls),
        "controls_significant": sum(1 for r in ctrls if r["fdr"] < 0.05),
        "min_count": args.min_count,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    plots = make_plots(out, rows, per_bin, bins, args.tail)

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


def _guide_metadata(file_id: str) -> dict:
    """guide_id -> intended target and guide type, from the library table."""
    st, obj = rp.portal_json(
        f"{file_id}?format=json" if file_id.startswith("/")
        else f"/tabular-files/{file_id}/?format=json")
    href = (obj or {}).get("href")
    out = {"target": {}, "type": {}}
    if not href:
        return out
    dest = rp.REF_DIR / "guides" / Path(href).name
    if not dest.exists():
        rp.portal_download(href, dest)
    text = rp.read_text_maybe_gzip(dest)
    lines = text.splitlines()
    if not lines:
        return out
    delim = max(("\t", ",", ";"), key=lambda d: len(lines[0].split(d)))
    hdr = [h.strip().lower() for h in lines[0].split(delim)]

    def col(*names):
        for n in names:
            if n in hdr:
                return hdr.index(n)
        return None

    gi = col("guide_id", "guide", "name", "id")
    ti = col("intended_target_name", "target", "genomic_element")
    ty = col("type", "guide_type", "targeting")
    for ln in lines[1:]:
        parts = ln.split(delim)
        if gi is None or len(parts) <= gi:
            continue
        gid = parts[gi].strip()
        if ti is not None and len(parts) > ti:
            out["target"][gid] = parts[ti].strip() or gid
        else:
            out["target"][gid] = gid
        if ty is not None and len(parts) > ty:
            out["type"][gid] = parts[ty].strip()
    return out


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
    a.add_argument("--label")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"discover": cmd_discover, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
