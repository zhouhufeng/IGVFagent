#!/usr/bin/env python3
"""Split an IGVF 'reporter experiment barcode' file into pipeline inputs.

The IGVF standard packs every replicate into one wide file
(barcode, oligo_name, dna_count_<n>, rna_count_<n>, ...). MPRAsnakeflow's
count stage consumes one narrow file per replicate plus an assignment, so
this splits the wide file back out.

A barcode not observed in a replicate has BOTH count cells empty; such a
row is omitted from that replicate rather than written as a zero, which
is what the upstream per-replicate count files look like.
"""
import gzip
import re
import sys
from pathlib import Path


def main(barcode_file: str, out_dir: str) -> int:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    op = gzip.open if str(barcode_file).endswith(".gz") else open
    with op(barcode_file, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        reps = sorted({m.group(1) for c in cols
                       if (m := re.match(r"dna_count_(.+)", c))})
        if not reps:
            sys.exit(f"No dna_count_<replicate> columns in {barcode_file}")
        handles = {r: open(out / f"rep{r}.tsv", "w") for r in reps}
        assign = open(out / "assignment.tsv", "w")
        kept = {r: 0 for r in reps}
        n = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            bc, oligo = f[idx["barcode"]], f[idx["oligo_name"]]
            assign.write(f"{bc}\t{oligo}\texact\t1/1\n")
            n += 1
            for r in reps:
                d = f[idx[f"dna_count_{r}"]].strip()
                v = f[idx[f"rna_count_{r}"]].strip()
                if d == "" and v == "":
                    continue
                handles[r].write(f"{bc}\t{d or 0}\t{v or 0}\n")
                kept[r] += 1
    for h in handles.values():
        h.close()
    assign.close()
    print(f"barcodes: {n:,}")
    print(f"replicates: {', '.join(reps)}")
    for r in reps:
        print(f"  rep{r}: {kept[r]:,} observed barcodes")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: prepare.py <reporter_experiment_barcode.tsv.gz> <out_dir>")
    sys.exit(main(sys.argv[1], sys.argv[2]))
