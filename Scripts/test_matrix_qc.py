"""Threshold aggregation and scope labelling, asserted.

Both come from a hosted test on released IGVF SC-islet data (2026-09-15):

* An agent-authored inspector computed X_max from the first 50 rows and
  returned it under a name that reads as a whole-matrix statistic. It said 14;
  the true maximum was 3,724. The mistake was invisible in the output.
* The per-barcode QC table was correct, and "count rows meeting two numeric
  conditions, and take medians" still could not be completed.

    python3 Scripts/test_matrix_qc.py
"""
from __future__ import annotations

import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matrix_qc_skill as M  # noqa: E402

FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="igvf-mqc-"))

    print("\nrule parsing")
    check("two conditions", M.parse_where("a>=1000,b>=200"),
          [("a", ">=", 1000.0), ("b", ">=", 200.0)])
    check("mixed operators", M.parse_where("a<5,b!=0"),
          [("a", "<", 5.0), ("b", "!=", 0.0)])

    print("\nthreshold aggregation against an independent count")
    import random
    random.seed(7)
    rows = [(f"BC{i}", random.randint(1000, 8000), random.randint(200, 4000))
            for i in range(1858)]
    rows += [(f"BX{i}", random.randint(0, 999), random.randint(0, 199))
             for i in range(20000)]
    random.shuffle(rows)
    table = tmp / "qc.tsv"
    with open(table, "w") as fh:
        fh.write("barcode\ttotal_counts\tn_genes_by_counts\n")
        for b, t, g in rows:
            fh.write(f"{b}\t{t}\t{g}\n")

    passing = [(t, g) for _, t, g in rows if t >= 1000 and g >= 200]

    import json
    out = tmp / "out"
    args = M.build_parser().parse_args(
        ["threshold", "--table", str(table),
         "--where", "total_counts>=1000,n_genes_by_counts>=200",
         "--median", "total_counts,n_genes_by_counts",
         "--out-dir", str(out)])
    args.func(args)
    res = json.loads((out / "qc_threshold.json").read_text())

    check("rows counted", res["rows_total"], len(rows))
    check("rows passing matches an independent count",
          res["rows_passing"], len(passing))
    check("passing count is the planted 1,858", res["rows_passing"], 1858)
    check("median total_counts after filtering",
          res["medians_after"]["total_counts"],
          statistics.median(x[0] for x in passing))
    check("median n_genes after filtering",
          res["medians_after"]["n_genes_by_counts"],
          statistics.median(x[1] for x in passing))

    print("\na missing column is an error, not a silent zero")
    rc = subprocess.run(
        [sys.executable, str(HERE / "matrix_qc_skill.py"), "threshold",
         "--table", str(table), "--where", "nonexistent>=1"],
        capture_output=True, text=True).returncode
    check("exits non-zero", rc != 0, True)

    print("\nscope labelling (needs anndata; skipped without it)")
    try:
        import anndata  # noqa: F401
        import numpy as np
        import scipy.sparse as sp
    except ImportError:
        print("  SKIP  anndata/scipy not installed here")
        print()
        if FAILED:
            print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
            return 1
        print("all checks passed")
        return 0

    X = sp.random(500, 300, density=0.05, format="csr", random_state=0)
    X.data = np.round(X.data * 10) + 1
    X = X.tolil()
    X[400, 250] = 3724.0          # deliberately outside any 50-row sample
    A = anndata.AnnData(X=X.tocsr())
    h5 = tmp / "probe.h5ad"
    A.write_h5ad(str(h5))

    full = M._full_x_stats(A.X, A.n_obs, 128)
    check("the full reduction finds the planted maximum", full["X_max"], 3724.0)
    check("...and says where it is",
          (full["X_argmax_row"], full["X_argmax_col"]), (400, 250))
    check("...and labels itself full", full["scope"], "full")

    sampled = M._sampled_x_stats(A.X, A.n_obs, 50)
    check("a 50-row sample misses it", sampled["sample_X_max"] < 3724.0, True)
    check("...and labels itself sampled", sampled["scope"], "sampled")
    check("...and no key reads as a whole-matrix statistic",
          [k for k in sampled if k.startswith("X_")], [])
    check("...and carries an explicit warning",
          "not of the" in sampled["warning"], True)

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
