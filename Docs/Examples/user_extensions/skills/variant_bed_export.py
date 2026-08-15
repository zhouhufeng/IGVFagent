"""Export a variant list (chrom/pos TSV or CSV) to a sorted BED file.

User-skill template. Drop this file in ``~/.igvfagent/skills/`` or
``UserExtensions/skills/`` and it becomes a first-class subcommand:

    igvfagent variant-bed-export my_variants.tsv --out variants.bed

Conventions worth copying:

  * module docstring — its first line is the description shown by
    ``igvfagent --help`` and ``igvfagent extensions``
  * a ``main()`` entrypoint with its own argparse parser (the CLI
    dispatcher hands over ``sys.argv``; return an int exit code)
  * announce every artefact on stdout as ``Report: <path>`` — that is
    the machine contract the agent runtime parses to chain steps, and
    the localstore harvester picks the files up into the local KG/DB
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert a chrom/pos variant table to a sorted BED file.")
    parser.add_argument("table", help="TSV/CSV with chrom + pos columns "
                                       "(header row required).")
    parser.add_argument("--out", default=None,
                        help="Output BED path (default: <table>.bed).")
    parser.add_argument("--chrom-col", default="chrom")
    parser.add_argument("--pos-col", default="pos")
    args = parser.parse_args()

    src = Path(args.table)
    if not src.is_file():
        sys.stderr.write(f"no such file: {src}\n")
        return 2
    out = Path(args.out) if args.out else src.with_suffix(".bed")

    delim = "," if src.suffix.lower() == ".csv" else "\t"
    rows = []
    with src.open(newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter=delim):
            try:
                chrom = rec[args.chrom_col].strip()
                pos = int(rec[args.pos_col])
            except (KeyError, ValueError):
                continue
            rows.append((chrom, pos - 1, pos))
    if not rows:
        sys.stderr.write(f"no usable {args.chrom_col}/{args.pos_col} rows "
                         f"in {src}\n")
        return 1

    rows.sort()
    with out.open("w", encoding="utf-8") as fh:
        for chrom, start, end in rows:
            fh.write(f"{chrom}\t{start}\t{end}\n")

    print(f"{len(rows)} variants -> BED")
    print(f"Report: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
