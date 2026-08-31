#!/usr/bin/env python3
"""Score an IGVFagent master table against the published IGVF artefact.

Both files carry one row per (replicate, oligo_name). This compares every
shared numeric column value-by-value and writes a summary.json that
Benchmarks/concordance.py scores.
"""
import gzip
import json
import sys
from pathlib import Path

COLS = ["dna_counts", "rna_counts", "n_bc",
        "dna_normalized", "rna_normalized", "log2FoldChange"]


def load(path: str) -> dict:
    op = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(cols, line.rstrip("\n").split("\t")))
            out[(d["replicate"], d["oligo_name"])] = d
    return out


def main(mine_path: str, theirs_path: str, out_json: str) -> int:
    mine, theirs = load(mine_path), load(theirs_path)
    km, kt = set(mine), set(theirs)
    shared = km & kt
    per_col = {}
    max_abs = 0.0
    n_diff_total = 0
    for c in COLS:
        mx, nd = 0.0, 0
        for k in shared:
            try:
                a, b = float(mine[k][c]), float(theirs[k][c])
            except (KeyError, ValueError):
                continue
            d = abs(a - b)
            mx = max(mx, d)
            if d > 1e-9:
                nd += 1
        per_col[c] = {"max_abs_diff": mx, "rows_differing": nd}
        max_abs = max(max_abs, mx)
        n_diff_total += nd

    summary = {
        "rows_mine": len(mine),
        "rows_published": len(theirs),
        "rows_shared": len(shared),
        "keys_identical": int(km == kt),
        "only_mine": len(km - kt),
        "only_published": len(kt - km),
        "values_compared": len(shared) * len(COLS),
        "max_abs_diff": max_abs,
        "rows_differing_total": n_diff_total,
        "per_column": per_col,
    }
    Path(out_json).write_text(json.dumps(summary, indent=2))
    print(f"rows: mine={len(mine):,} published={len(theirs):,} "
          f"identical_keys={km == kt}")
    print(f"values compared: {summary['values_compared']:,}")
    print(f"max |difference|: {max_abs:g}")
    for c, s in per_col.items():
        verdict = "EXACT" if s["rows_differing"] == 0 else "DIFFERS"
        print(f"  {c:17s} max={s['max_abs_diff']:<12g} differing="
              f"{s['rows_differing']:>7,}  {verdict}")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: compare.py <mine.tsv.gz> <published.tsv.gz> <summary.json>")
    sys.exit(main(*sys.argv[1:]))
