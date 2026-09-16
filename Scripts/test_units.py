"""Impossible statistics, caught before an answer ships.

From a hosted retest: the agent wrote "p_adj capped at 240". An adjusted P
value cannot exceed 1. The tool was right — it returned neg_log10_pvalue, for
which 240 is ordinary — and only the final prose renamed the field. Nothing
about "240" looks wrong until you know which field it belongs to, so a reader
has no chance; the name together with the value is contradictory on its face.

    python3 Scripts/test_units.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _units as U  # noqa: E402

FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def main() -> int:
    print("\nthe reported failure, in the phrasing it actually appeared in")
    check("p_adj capped at 240 is caught",
          bool(U.scan_text("The GATA3 trans summary shows p_adj capped at 240 "
                           "for the top hit.")), True)
    check("...naming the field it probably was",
          "neg_log10_pvalue" in U.scan_text("p_adj capped at 240")[0], True)

    print("\nother impossible combinations")
    for text in ("fdr = 12.5", "correlation = 1.8", "auc: 95",
                 "probability of 7.2", "mapping_rate was 87",
                 "p_adj of 240", "neg_log10_pvalue = -3"):
        check(f"caught: {text}", bool(U.scan_text(text)), True)

    print("\nmarkdown tables are read through their header, not by cell")
    tbl = ("| gene | p_adj | neg_log10_pvalue |\n|---|---|---|\n"
           "| GATA3 | 240 | 240 |\n| SOX9 | 0.003 | 2.5 |")
    hits = U.scan_text(tbl)
    check("the p_adj cell is caught", len(hits), 1)
    check("...and the neg_log10_pvalue cell of the same value is not",
          "neg_log10" in hits[0].split(":")[0], False)
    check("a number under a column with no declared range is ignored",
          U.scan_text("| gene | score |\n|---|---|\n| GATA3 | 240 |"), [])

    print("\nlegitimate values stay quiet — a false alarm is worse than a miss")
    for text in ("neg_log10_pvalue = 240 for the top hit",
                 "p_adj = 0.0031 after BH correction",
                 "total_counts = 3724 in that barcode",
                 "n_genes_by_counts = 1857",
                 "the correlation was 0.82",
                 "percent mitochondrial = 12.4",
                 "X_max = 3724",
                 "score = 0.9999999981",
                 "auc: 0.93",
                 "1,858 barcodes passed the filter",
                 "211,901 barcodes and 62,757 genes"):
        check(f"quiet: {text[:46]}", U.scan_text(text), [])

    print("\nsingle-value checks")
    check("p_adj 240 rejected", bool(U.check_value("p_adj", 240)), True)
    check("p_adj 0.04 accepted", U.check_value("p_adj", 0.04), None)
    check("unknown field ignored", U.check_value("widgets", 999), None)
    check("non-numeric ignored", U.check_value("p_adj", "n/a"), None)
    check("NaN ignored", U.check_value("p_adj", float("nan")), None)

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
