#!/usr/bin/env python3
"""Activity-by-Contact (ABC) enhancer-gene prediction.

Clean-room implementation of the ABC model (Fulco et al., Nat Genet 2019;
Nasser et al., Nature 2021; broadinstitute/ABC-Enhancer-Gene-Prediction, MIT).
No source copied. The upstream pipeline is Snakemake around macs2, samtools,
bedtools and juicer; none of those are installed here, and peak calling from
BAM is genuinely out of scope. What IS in scope, and what this does, is the
model itself -- starting from candidate elements and signal tracks, which
ENCODE publishes and IGVFagent already fetches.

THE MODEL, in full:

    ABC(E,G) = Activity(E) x Contact(E,G) / SUM over elements E' within the
               window of Activity(E') x Contact(E',G)

    Activity(E) = geometric mean of the accessibility and H3K27ac signal in E

It is a SHARE, not a score: the denominator makes every gene's predictions sum
to 1 across its neighbourhood. That is why an element can be strong in
absolute signal and still score low -- if it sits among stronger neighbours,
its share of that gene's regulatory input is small. Reading ABC as "how active
is this enhancer" gets this backwards.

CONTACT. With Hi-C, the observed contact frequency. Without it, the power law
Fulco 2019 fit genome-wide, contact ~ distance^gamma with gamma about -0.87
(the `--gamma` default). Nasser 2021 showed the power-law version performs
close to the Hi-C version, which is what makes this useful without a Hi-C map
-- and IGVFagent has no Hi-C for most biosamples.

WHAT THIS IS NOT. It does not call peaks (bring candidate elements), does not
handle BAMs (bring bigWigs), and does not reproduce the upstream's exact
quantile normalisation. Concordant, not bit-identical.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "ABC"

# Fulco 2019's genome-wide power-law fit.
DEFAULT_GAMMA = -0.87
# Elements further than this from a TSS are not considered for that gene.
DEFAULT_WINDOW = 5_000_000
# Nasser 2021's threshold for a "positive" prediction.
DEFAULT_THRESHOLD = 0.02


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def read_bed(path: str, name_col: int = 3):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                rec = {"chrom": f[0], "start": int(f[1]), "end": int(f[2])}
            except ValueError:
                continue
            rec["name"] = f[name_col] if len(f) > name_col else \
                f"{rec['chrom']}:{rec['start']}-{rec['end']}"
            if len(f) > 5 and f[5] in ("+", "-"):
                rec["strand"] = f[5]
            out.append(rec)
    return out


def signal_in(bw, chrom: str, start: int, end: int) -> float:
    """Mean signal over an interval; 0 when the track lacks the contig."""
    try:
        if chrom not in bw.chroms():
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt in bw.chroms():
                chrom = alt
            else:
                return 0.0
        clen = bw.chroms()[chrom]
        s, e = max(0, start), min(end, clen)
        if e <= s:
            return 0.0
        v = bw.stats(chrom, s, e, type="mean", exact=True)[0]
        return float(v) if v is not None else 0.0
    except Exception:                                       # noqa: BLE001
        return 0.0


def cmd_score(args) -> int:
    setup_logging()
    import numpy as np
    import pyBigWig

    elements = read_bed(args.elements)
    genes = read_bed(args.genes)
    if not elements or not genes:
        raise SystemExit("need non-empty --elements and --genes BED files")
    logging.info("%d candidate elements, %d genes", len(elements), len(genes))

    atac = pyBigWig.open(args.atac)
    h3k = pyBigWig.open(args.h3k27ac) if args.h3k27ac else None

    # Activity = geometric mean of the two marks. Geometric, not arithmetic,
    # so an element must be BOTH accessible and acetylated to score highly --
    # an arithmetic mean lets one huge signal carry an element that the other
    # assay says is inert.
    for e in elements:
        a = signal_in(atac, e["chrom"], e["start"], e["end"])
        if h3k is not None:
            k = signal_in(h3k, e["chrom"], e["start"], e["end"])
            e["activity"] = math.sqrt(max(a, 0.0) * max(k, 0.0))
        else:
            e["activity"] = max(a, 0.0)
        e["atac"] = a
        e["h3k27ac"] = k if h3k is not None else float("nan")
    atac.close()
    if h3k is not None:
        h3k.close()

    by_chrom: "dict[str, list]" = {}
    for e in elements:
        by_chrom.setdefault(e["chrom"], []).append(e)

    hic = None
    if args.hic:
        hic = {}
        with open(args.hic) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 3:
                    try:
                        hic[(f[0], int(f[1]))] = float(f[2])
                    except ValueError:
                        pass
        logging.info("loaded %d Hi-C bins", len(hic))

    rows = []
    for g in genes:
        tss = g["start"] if g.get("strand") != "-" else g["end"]
        cands = [e for e in by_chrom.get(g["chrom"], [])
                 if abs(((e["start"] + e["end"]) // 2) - tss) <= args.window]
        if not cands:
            continue
        num = []
        for e in cands:
            mid = (e["start"] + e["end"]) // 2
            dist = max(abs(mid - tss), args.min_distance)
            if hic is not None:
                contact = hic.get((g["chrom"], (mid // args.hic_resolution)
                                    * args.hic_resolution), 0.0)
                if contact <= 0:
                    contact = dist ** args.gamma
            else:
                contact = dist ** args.gamma
            num.append((e, e["activity"] * contact, dist, contact))
        denom = sum(v for _, v, _, _ in num)
        if denom <= 0:
            continue
        for e, v, dist, contact in num:
            score = v / denom
            if score < args.min_report:
                continue
            rows.append({
                "gene": g["name"], "chrom": g["chrom"], "tss": tss,
                "element": e["name"], "element_start": e["start"],
                "element_end": e["end"], "distance": dist,
                "activity": round(e["activity"], 6),
                "atac": round(e["atac"], 6),
                "h3k27ac": ("" if e["h3k27ac"] != e["h3k27ac"]
                            else round(e["h3k27ac"], 6)),
                "contact": f"{contact:.6g}",
                "abc_score": round(score, 6),
                "predicted": "TRUE" if score >= args.threshold else "FALSE",
            })

    rows.sort(key=lambda r: -r["abc_score"])
    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_abc"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "abc_predictions.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]) if rows
                           else ["gene", "element", "abc_score"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    npos = sum(1 for r in rows if r["predicted"] == "TRUE")
    summary = {"n_elements": len(elements), "n_genes": len(genes),
               "n_pairs": len(rows), "n_predicted": npos,
               "threshold": args.threshold, "gamma": args.gamma,
               "window": args.window,
               "contact": "Hi-C" if hic is not None else
                          f"power law distance^{args.gamma}",
               "activity": ("sqrt(ATAC x H3K27ac)" if args.h3k27ac
                            else "ATAC only (no H3K27ac given)")}
    js = out / "abc_summary.json"
    js.write_text(json.dumps(summary, indent=2))

    print(f"Report:  {tsv}")
    print(f"Summary: {js}")
    print(f"  {len(rows):,} element-gene pairs, {npos:,} at ABC >= {args.threshold}")
    print(f"  contact: {summary['contact']}  ·  activity: {summary['activity']}")
    for r in rows[:8]:
        print(f"    {r['gene']:12} {r['element']:22} ABC={r['abc_score']:.4f} "
              f"d={r['distance']:,}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent abc",
        description="Activity-by-Contact enhancer-gene prediction from "
                     "candidate elements + signal tracks.")
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("score", help="Candidate elements + signal -> ABC scores.")
    s.add_argument("--elements", required=True, help="Candidate element BED.")
    s.add_argument("--genes", required=True, help="Gene/TSS BED (name in col 4).")
    s.add_argument("--atac", required=True, help="ATAC/DNase bigWig.")
    s.add_argument("--h3k27ac", help="H3K27ac bigWig. Without it activity is "
                                      "accessibility alone, which is weaker.")
    s.add_argument("--hic", help="Optional 3-col contact file: chrom bin value.")
    s.add_argument("--hic-resolution", type=int, default=5000)
    s.add_argument("--gamma", type=float, default=DEFAULT_GAMMA)
    s.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    s.add_argument("--min-distance", type=int, default=1000)
    s.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    s.add_argument("--min-report", type=float, default=0.001)
    s.add_argument("--label")
    s.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"score": cmd_score}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
