#!/usr/bin/env python3
"""scNT-seq: split new from old RNA by 4sU T>C conversions.

Clean-room implementation of the conversion-counting method behind scNT-seq
(Qiu et al., Nat Methods 2020; hongjie7/scNT-seq_pipeline). No source ported.

THE CHEMISTRY, because the analysis only makes sense from it. Cells are fed
4-thiouridine, which is incorporated into RNA made DURING the labelling
window. TimeLapse/SLAM chemistry then converts 4sU to a cytosine analogue, so
newly made transcripts carry T>C mismatches against the reference while
pre-existing RNA does not. Counting those conversions per read splits new from
old RNA in the same library.

NOT THE SAME AS kb's `nac` WORKFLOW, which IGVFagent already has. `nac` calls
a transcript nascent from INTRON content -- unspliced means recently
transcribed. This calls it new from CHEMICAL LABELLING. They answer different
questions: an intronless transcript made an hour ago is new here and mature
there. Both are "nascent" in the vocabulary, which is exactly why they are
separate tools.

STRAND MATTERS AND IS THE EASY THING TO GET WRONG. The conversion is T>C on
the transcribed strand. A read aligned to the reverse strand shows it as A>G.
Counting only T>C would silently halve the signal and bias it by gene
orientation, so reverse-strand reads are scored on A>G.

BACKGROUND. Sequencing error and SNPs also produce T>C. A single conversion
is therefore weak evidence; the default requires `--min-conversions 2`, and
the per-cell background rate from all OTHER substitution types is reported so
a caller can see whether the threshold is doing any work on their data.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "scNTseq"


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def count_read_conversions(read, min_qual: int):
    """(n_labelling_conversions, n_other_substitutions) for one read.

    Uses the reference base from the MD tag via get_aligned_pairs(with_seq),
    so no FASTA is needed. Low-quality bases are skipped: a miscalled base is
    exactly what this must not count as a labelling event.
    """
    conv = other = 0
    try:
        pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
    except ValueError:
        return 0, 0                      # no MD tag on this read
    seq = read.query_sequence
    qual = read.query_qualities
    if seq is None:
        return 0, 0
    # On a reverse-strand read the T>C conversion presents as A>G.
    want = ("A", "G") if read.is_reverse else ("T", "C")
    for qpos, _rpos, rbase in pairs:
        if qpos is None or rbase is None:
            continue
        rb = rbase.upper()
        qb = seq[qpos].upper()
        if rb == qb:
            continue
        if qual is not None and qpos < len(qual) and qual[qpos] < min_qual:
            continue
        if rb == want[0] and qb == want[1]:
            conv += 1
        else:
            other += 1
    return conv, other


def cmd_count(args) -> int:
    setup_logging()
    import pysam

    bam = pysam.AlignmentFile(args.bam, "rb")
    per_cell_new: Counter = Counter()
    per_cell_old: Counter = Counter()
    per_cell_conv: Counter = Counter()
    per_cell_other: Counter = Counter()
    per_gene: "dict[str, list[int]]" = defaultdict(lambda: [0, 0])
    n_reads = n_used = 0
    t0 = time.time()

    for read in bam.fetch(until_eof=True):
        n_reads += 1
        if args.max_reads and n_reads > args.max_reads:
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < args.min_mapq:
            continue
        cb = read.get_tag(args.cell_tag) if read.has_tag(args.cell_tag) else None
        if args.require_cell and cb is None:
            continue
        gene = read.get_tag(args.gene_tag) if read.has_tag(args.gene_tag) else None
        conv, other = count_read_conversions(read, args.min_base_qual)
        n_used += 1
        cell = cb or "_nocell"
        per_cell_conv[cell] += conv
        per_cell_other[cell] += other
        is_new = conv >= args.min_conversions
        if is_new:
            per_cell_new[cell] += 1
        else:
            per_cell_old[cell] += 1
        if gene:
            per_gene[gene][0 if is_new else 1] += 1

    bam.close()
    logging.info("%d reads seen, %d used, %.1fs", n_reads, n_used,
                 time.time() - t0)
    if n_used == 0:
        raise SystemExit(
            "no usable reads. Check --min-mapq, and that the BAM carries MD "
            "tags (conversions are read from them) and the cell tag "
            f"{args.cell_tag!r}.")

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_scnt"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)

    cells = sorted(set(per_cell_new) | set(per_cell_old))
    ctsv = out / "per_cell_new_old.tsv"
    with ctsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cell", "new_reads", "old_reads", "new_fraction",
                    "labelling_conversions", "other_substitutions",
                    "background_ratio"])
        for c in cells:
            nn, no = per_cell_new[c], per_cell_old[c]
            tot = nn + no
            # The labelling signal is only meaningful relative to the rate of
            # every OTHER substitution type, which is pure noise.
            bg = (per_cell_conv[c] / per_cell_other[c]) if per_cell_other[c] else float("nan")
            w.writerow([c, nn, no, f"{nn / tot:.4f}" if tot else "",
                        per_cell_conv[c], per_cell_other[c],
                        "" if bg != bg else f"{bg:.3f}"])

    gtsv = out / "per_gene_new_old.tsv"
    with gtsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "new_reads", "old_reads", "new_fraction"])
        for g, (nn, no) in sorted(per_gene.items(), key=lambda kv: -sum(kv[1])):
            tot = nn + no
            w.writerow([g, nn, no, f"{nn / tot:.4f}" if tot else ""])

    tot_new = sum(per_cell_new.values())
    tot_old = sum(per_cell_old.values())
    tot_conv = sum(per_cell_conv.values())
    tot_other = sum(per_cell_other.values())
    summary = {
        "reads_seen": n_reads, "reads_used": n_used,
        "new_reads": tot_new, "old_reads": tot_old,
        "new_fraction": round(tot_new / (tot_new + tot_old), 4)
        if (tot_new + tot_old) else None,
        "labelling_conversions": tot_conv,
        "other_substitutions": tot_other,
        "signal_to_background": round(tot_conv / tot_other, 3) if tot_other else None,
        "min_conversions": args.min_conversions,
        "cells": len(cells), "genes": len(per_gene),
    }
    js = out / "scnt_summary.json"
    js.write_text(json.dumps(summary, indent=2))
    print(f"Report:  {ctsv}")
    print(f"Wrote:   {gtsv}")
    print(f"Summary: {js}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    if summary["signal_to_background"] is not None and \
            summary["signal_to_background"] < 1.5:
        print("  WARNING: labelling conversions are barely above the rate of "
              "other substitutions. Either labelling did not work, or these "
              "'new' calls are sequencing error.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent scnt-seq",
        description="scNT-seq: split new from old RNA by 4sU T>C conversions.")
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("count", help="BAM -> per-cell / per-gene new vs old.")
    c.add_argument("--bam", required=True, help="Aligned BAM with MD tags.")
    c.add_argument("--cell-tag", default="CB")
    c.add_argument("--gene-tag", default="GX")
    c.add_argument("--min-conversions", type=int, default=2,
                   help="Conversions needed to call a read NEW. 1 is weak: "
                        "sequencing error and SNPs also make T>C.")
    c.add_argument("--min-base-qual", type=int, default=20)
    c.add_argument("--min-mapq", type=int, default=10)
    c.add_argument("--require-cell", action="store_true")
    c.add_argument("--max-reads", type=int, default=0)
    c.add_argument("--label")
    c.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"count": cmd_count}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
