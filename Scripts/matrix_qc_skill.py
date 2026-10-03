#!/usr/bin/env python3
"""Honest summaries and threshold aggregation for single-cell matrices.

Two gaps found by a hosted test on released IGVF SC-islet data (2026-09-15),
both of which pushed the agent into authoring its own tools badly.

**A sampled statistic was reported as a whole-matrix one.** An agent-authored
inspector computed ``X_max`` from ``X[:50]`` and returned it under a name that
claims to describe the matrix. On IGVFFI6543TDVW it reported 14; the true
maximum is 3,724, at row 21,043 / column 34,645. Nothing in the output said
"sampled", so neither the model nor the reader could tell. Here, a full
reduction is the default, sampling must be asked for by name, and a sampled
result is labelled ``scope: sampled`` with the row count that produced it --
the field cannot be read as a whole-matrix figure by accident.

**A two-column threshold count needed no new data and could not be done.** The
per-barcode QC table was correct; the remaining step was "count rows meeting
two numeric conditions and take medians". That failed through a
warehouse-query path that wanted an initialised warehouse, and two attempts to
author a small aggregator were rejected as overlapping with a GSEA tool on the
shared words "stats" and "tsv". ``threshold`` is that operation, first-class,
reading a plain table and needing nothing initialised.

Memory is bounded on purpose: the full reduction walks the sparse matrix in row
blocks and never materialises a dense array, because the matrices this is aimed
at are ~200k x ~60k and a dense copy is ~100 GB.

    igvfagent matrix-qc summary --h5ad matrix.h5ad
    igvfagent matrix-qc summary --h5ad matrix.h5ad --sample-rows 50
    igvfagent matrix-qc threshold --table qc.tsv \\
        --where "total_counts>=1000,n_genes_by_counts>=200" \\
        --median total_counts,n_genes_by_counts
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "MatrixQC"
log = logging.getLogger("matrix-qc")

DEFAULT_BLOCK = 20_000


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")


# --------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------


def _full_x_stats(X, n_obs: int, block: int) -> dict:
    """Reduce the whole matrix in row blocks, densifying one block at a time."""
    import numpy as np

    out = {
        "scope": "full", "rows_examined": int(n_obs),
        "X_min": math.inf, "X_max": -math.inf,
        "X_argmax_row": -1, "X_argmax_col": -1,
        "nnz": 0, "X_has_negatives": False,
        "X_sum": 0.0,
    }
    integer_like = True
    for start in range(0, n_obs, block):
        stop = min(start + block, n_obs)
        chunk = X[start:stop]
        dense = chunk.toarray() if hasattr(chunk, "toarray") else np.asarray(chunk)
        if dense.size == 0:
            continue
        cmin = float(dense.min())
        cmax = float(dense.max())
        if cmin < out["X_min"]:
            out["X_min"] = cmin
        if cmax > out["X_max"]:
            out["X_max"] = cmax
            flat = int(np.argmax(dense))
            out["X_argmax_row"] = int(start + flat // dense.shape[1])
            out["X_argmax_col"] = int(flat % dense.shape[1])
        out["nnz"] += int((dense != 0).sum())
        out["X_sum"] += float(dense.sum())
        if cmin < 0:
            out["X_has_negatives"] = True
        if integer_like:
            nz = dense[dense != 0]
            if nz.size and not bool(np.all(np.abs(nz - np.round(nz)) < 1e-6)):
                integer_like = False
    out["X_integer_valued"] = integer_like
    if not math.isfinite(out["X_min"]):
        out["X_min"] = 0.0
    if not math.isfinite(out["X_max"]):
        out["X_max"] = 0.0
    return out


def _sampled_x_stats(X, n_obs: int, rows: int) -> dict:
    import numpy as np

    stop = min(rows, n_obs)
    chunk = X[:stop]
    dense = chunk.toarray() if hasattr(chunk, "toarray") else np.asarray(chunk)
    # Every field says "sample". There is no key here that a reader could
    # mistake for a property of the whole matrix -- that confusion is the bug
    # this command exists to avoid repeating.
    return {
        "scope": "sampled",
        "rows_examined": int(stop),
        "rows_total": int(n_obs),
        "sample_X_min": float(dense.min()) if dense.size else 0.0,
        "sample_X_max": float(dense.max()) if dense.size else 0.0,
        "sample_nnz": int((dense != 0).sum()),
        "warning": (f"These are statistics of the FIRST {stop} rows, not of "
                    f"the {n_obs}-row matrix. Do not report sample_X_max as "
                    f"the matrix maximum. Re-run without --sample-rows for "
                    f"the full reduction."),
    }


def cmd_summary(args) -> int:
    try:
        import anndata
    except ImportError:
        print("ERROR: anndata is required for this command.", file=sys.stderr)
        return 2

    path = Path(args.h5ad)
    if not path.exists():
        print(f"ERROR: no such file: {path}", file=sys.stderr)
        return 2

    A = anndata.read_h5ad(str(path), backed="r" if args.backed else None)
    result = {
        "file": str(path),
        "n_obs": int(A.n_obs),
        "n_vars": int(A.n_vars),
        "obs_columns": list(map(str, A.obs.columns)),
        "var_columns": list(map(str, A.var.columns)),
        "obs_index_name": str(A.obs.index.name or ""),
        "var_index_name": str(A.var.index.name or ""),
        "layers": list(map(str, A.layers.keys())),
        "obsm": list(map(str, A.obsm.keys())),
        "has_raw": A.raw is not None,
        "X_dtype": str(getattr(A.X, "dtype", "")),
    }
    X = A.X
    if args.sample_rows:
        result["X"] = _sampled_x_stats(X, A.n_obs, args.sample_rows)
    elif args.no_x_stats:
        result["X"] = {"scope": "not computed",
                       "note": "re-run without --no-x-stats for matrix values"}
    else:
        result["X"] = _full_x_stats(X, A.n_obs, args.block)

    out_dir = Path(args.out_dir) if args.out_dir else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{path.stem}_summary.json"
    out.write_text(json.dumps(result, indent=2))

    x = result["X"]
    print(f"{path.name}")
    print(f"  shape            : {result['n_obs']:,} obs x {result['n_vars']:,} vars")
    # Spelled out rather than nested in the f-string: an empty obs/var column
    # list is a real and confusing state (the index is not an annotation
    # column, and reading it as one is its own bug), so it gets said plainly.
    obs_cols = ", ".join(result["obs_columns"]) or (
        "(none — the index is %r, which is not an annotation column)"
        % (result["obs_index_name"] or "unnamed"))
    var_cols = ", ".join(result["var_columns"]) or (
        "(none — the index is %r, which is not an annotation column)"
        % (result["var_index_name"] or "unnamed"))
    print(f"  obs columns      : {obs_cols}")
    print(f"  var columns      : {var_cols}")
    print(f"  layers           : {', '.join(result['layers']) or '(none)'}")
    print(f"  raw object       : {'present' if result['has_raw'] else 'absent'}")
    print(f"  X dtype          : {result['X_dtype']}")
    print(f"  X scope          : {x['scope'].upper()}"
          + (f" ({x['rows_examined']:,} of {x.get('rows_total', result['n_obs']):,} rows)"
             if x["scope"] == "sampled" else
             f" ({x.get('rows_examined', 0):,} rows)"))
    if x["scope"] == "full":
        print(f"  X min / max      : {x['X_min']:g} / {x['X_max']:g}")
        print(f"  X argmax         : row {x['X_argmax_row']:,}, "
              f"col {x['X_argmax_col']:,}")
        print(f"  non-zero entries : {x['nnz']:,}")
        print(f"  integer-valued   : {x['X_integer_valued']}")
    elif x["scope"] == "sampled":
        print(f"  sample X min/max : {x['sample_X_min']:g} / {x['sample_X_max']:g}")
        print(f"\n  {x['warning']}")
    print(f"\n  written: {out}")
    ls.record_analysis("matrix-qc", subcommand="summary", label=path.stem,
                       inputs=[str(path)], outputs=[str(out)])
    return 0


# --------------------------------------------------------------------------
# threshold
# --------------------------------------------------------------------------

_OPS = ("<=", ">=", "!=", "<", ">", "=")


def parse_where(spec: str) -> "list[tuple[str, str, float]]":
    """`total_counts>=1000,n_genes_by_counts>=200` -> [(col, op, value), …]"""
    out = []
    for clause in (spec or "").split(","):
        clause = clause.strip()
        if not clause:
            continue
        for op in _OPS:
            if op in clause:
                col, _, val = clause.partition(op)
                try:
                    out.append((col.strip(), op, float(val.strip())))
                except ValueError:
                    raise SystemExit(f"ERROR: {clause!r} — {val!r} is not a number")
                break
        else:
            raise SystemExit(f"ERROR: {clause!r} has no comparison operator "
                             f"(one of {', '.join(_OPS)})")
    return out


def _passes(row, conditions) -> bool:
    for col, op, want in conditions:
        try:
            v = float(row[col])
        except (KeyError, TypeError, ValueError):
            return False
        if op == ">=" and not v >= want:
            return False
        if op == "<=" and not v <= want:
            return False
        if op == ">" and not v > want:
            return False
        if op == "<" and not v < want:
            return False
        if op == "=" and not v == want:
            return False
        if op == "!=" and not v != want:
            return False
    return True


def _median(values) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2


def cmd_threshold(args) -> int:
    path = Path(args.table)
    if not path.exists():
        print(f"ERROR: no such table: {path}", file=sys.stderr)
        return 2
    conditions = parse_where(args.where)
    medians = [c.strip() for c in (args.median or "").split(",") if c.strip()]

    delim = "\t" if path.suffix.lower() in (".tsv", ".tab", ".txt") else ","
    if args.delimiter:
        delim = {"tab": "\t", "comma": ",", "\\t": "\t"}.get(args.delimiter,
                                                             args.delimiter)

    kept: "dict[str, list[float]]" = {c: [] for c in medians}
    before: "dict[str, list[float]]" = {c: [] for c in medians}
    total = n_pass = 0
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        missing = [c for c, _, _ in conditions if c not in (reader.fieldnames or [])]
        missing += [c for c in medians if c not in (reader.fieldnames or [])]
        if missing:
            print(f"ERROR: column(s) not in {path.name}: "
                  f"{', '.join(sorted(set(missing)))}\n"
                  f"available: {', '.join(reader.fieldnames or [])}",
                  file=sys.stderr)
            return 2
        for row in reader:
            total += 1
            for c in medians:
                try:
                    before[c].append(float(row[c]))
                except (TypeError, ValueError):
                    pass
            if _passes(row, conditions):
                n_pass += 1
                for c in medians:
                    try:
                        kept[c].append(float(row[c]))
                    except (TypeError, ValueError):
                        pass

    result = {
        "table": str(path),
        "rule": args.where,
        "rows_total": total,
        "rows_passing": n_pass,
        "fraction_passing": round(n_pass / total, 6) if total else 0.0,
        "medians_before": {c: _median(before[c]) for c in medians},
        "medians_after": {c: _median(kept[c]) for c in medians},
    }
    out_dir = Path(args.out_dir) if args.out_dir else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{path.stem}_threshold.json"
    out.write_text(json.dumps(result, indent=2, default=str))

    print(f"{path.name}")
    print(f"  rule             : {args.where}")
    print(f"  rows             : {total:,}")
    print(f"  passing          : {n_pass:,}  ({result['fraction_passing']:.4%})")
    for c in medians:
        print(f"  median {c:<22} before {_median(before[c]):g}   "
              f"after {_median(kept[c]):g}")
    print(f"\n  written: {out}")
    print("\n  Rows passing a numeric rule are QC-passing barcodes. They are "
          "not\n  validated cells — that requires a cell-calling method, not a "
          "threshold.")
    ls.record_analysis("matrix-qc", subcommand="threshold", label=path.stem,
                       inputs=[str(path)], outputs=[str(out)])
    return 0


# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent matrix-qc",
        description="Whole-matrix summaries that say what they measured, and "
                    "threshold aggregation that needs nothing initialised.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("summary", help="Structure and X statistics of an h5ad.")
    s.add_argument("--h5ad", required=True)
    s.add_argument("--out-dir")
    s.add_argument("--sample-rows", type=int, default=0,
                   help="Examine only the first N rows. The result is then "
                        "labelled sampled, with every field named accordingly.")
    s.add_argument("--no-x-stats", action="store_true",
                   help="Structure only; do not read matrix values.")
    s.add_argument("--block", type=int, default=DEFAULT_BLOCK,
                   help="Rows densified at a time during the full reduction.")
    s.add_argument("--backed", action="store_true",
                   help="Open backed on disk rather than into memory.")
    s.set_defaults(func=cmd_summary)

    s = sub.add_parser("threshold",
                       help="Count rows meeting numeric conditions, with "
                            "medians before and after.")
    s.add_argument("--table", required=True, help="TSV or CSV.")
    s.add_argument("--where", required=True,
                   help='e.g. "total_counts>=1000,n_genes_by_counts>=200"')
    s.add_argument("--median", default="",
                   help="Comma-separated columns to take medians of.")
    s.add_argument("--delimiter")
    s.add_argument("--out-dir")
    s.set_defaults(func=cmd_threshold)

    return p


def main(argv=None) -> int:
    setup_logging()
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
