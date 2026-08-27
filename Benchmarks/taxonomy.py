#!/usr/bin/env python3
"""Classify what each benchmark actually establishes.

Reviewers of this work objected, correctly, to describing the whole suite as
"reproduction": some benchmarks recover a published quantitative result, while
others confirm that records were retrieved, rows counted, or files written.
Those are different claims and should not share a word.

The classification here is derived from the assertions in each
``expected.json`` rather than assigned by hand, so it cannot drift from what
the suite actually checks. Each benchmark is assigned the **strongest** class
among its checks:

  A  quantitative reproduction  a derived scientific quantity is asserted
                                against the published value (AMI, AUPRC,
                                homogeneity, singlet rate, aggregated Z/P …)
  B  method validation          agreement with a reference implementation is
                                asserted (e.g. correlation against the R
                                original), rather than against a paper figure
  C  retrieval and enumeration  counts and identities of retrieved records
                                (row counts, cell counts, gene identity)
  D  artefact generation        only that expected files were produced

Class D is a smoke test. Class C shows the agent reached the right data. Only
A and B say an implementation computes the right answer, and B says it against
a reference rather than a publication.

Usage::

    python3 Benchmarks/taxonomy.py            # table
    python3 Benchmarks/taxonomy.py --json     # machine-readable
    python3 Benchmarks/taxonomy.py --markdown # for the manuscript
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

BENCH = Path(__file__).resolve().parent

# Quantities that are derived scientific measurements rather than tallies.
# Anything matching is class A; a numeric check on anything else is a count.
_METRIC_RE = re.compile(
    r"(auprc|auroc|precision|recall|\bami\b|homogeneity|correlation|corr\b"
    r"|singlet_rate|balance|_frac_|fraction|rate\b|mean\b|median"
    r"|aggregated_[zp]|omega|max_abs_diff|match_rate|concordance)",
    re.IGNORECASE)

# Agreement against a reference implementation rather than a published figure.
_REFIMPL_RE = re.compile(r"(validation_vs_|_vs_r\b|vs_reference|reference_)",
                         re.IGNORECASE)

CLASS_LABEL = {
    "A": "quantitative reproduction",
    "B": "method validation",
    "C": "retrieval and enumeration",
    "D": "artefact generation",
}


def _check_class(check: dict) -> str:
    ctype = (check.get("type") or "").lower()
    path = str(check.get("path") or check.get("filename")
               or check.get("glob") or "")
    if ctype == "artefact":
        return "D"
    if ctype in ("range", "in_set", "csv_value_present"):
        if _REFIMPL_RE.search(path):
            return "B"
        if _METRIC_RE.search(path):
            return "A"
        return "C"          # numeric, but a count or an identity
    if ctype in ("row_count_tsv", "csv_row_count"):
        return "C"
    return "D"


def classify(bench_dir: Path) -> "dict | None":
    exp = bench_dir / "expected.json"
    if not exp.is_file():
        return None
    try:
        meta = json.loads(exp.read_text())
    except (OSError, json.JSONDecodeError):
        return None

    checks = meta.get("checks") or []
    per = [_check_class(c) for c in checks]
    best = min(per, default="D")          # "A" < "B" < "C" < "D" lexically
    counts = {k: per.count(k) for k in "ABCD" if per.count(k)}

    run = bench_dir / "run.sh"
    text = run.read_text(errors="replace") if run.is_file() else ""
    return {
        "benchmark": bench_dir.name,
        "class": best,
        "class_label": CLASS_LABEL[best],
        "n_checks": len(checks),
        "checks_by_class": counts,
        "synthetic": bool(re.search(r"\bsimulate\b", text)),
        "paper": meta.get("paper", {}),
        "strongest_assertions": [
            (c.get("path") or c.get("filename") or c.get("glob") or c.get("name"))
            for c, k in zip(checks, per) if k == best
        ][:4],
    }


def all_benchmarks() -> "list[dict]":
    out = []
    for d in sorted(BENCH.iterdir()):
        if not d.is_dir():
            continue
        c = classify(d)
        if c:
            out.append(c)
    return out


def _summary(rows) -> dict:
    s = {k: 0 for k in "ABCD"}
    for r in rows:
        s[r["class"]] += 1
    return s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--markdown", action="store_true")
    a = ap.parse_args(argv)

    rows = all_benchmarks()
    summ = _summary(rows)

    if a.json:
        print(json.dumps({"n_benchmarks": len(rows), "summary": summ,
                          "classes": CLASS_LABEL, "benchmarks": rows}, indent=2))
        return 0

    if a.markdown:
        print("| benchmark | class | what it establishes | checks |")
        print("|---|---|---|---|")
        for r in sorted(rows, key=lambda r: (r["class"], r["benchmark"])):
            print(f"| `{r['benchmark']}` | {r['class']} | "
                  f"{r['class_label']} | {r['n_checks']} |")
        print()
        chk = {k: 0 for k in "ABCD"}
        for r in rows:
            for k, n in r["checks_by_class"].items():
                chk[k] += n
        print(f"**{len(rows)} benchmarks**, "
              + ", ".join(f"{summ[k]} {CLASS_LABEL[k]}" for k in "ABCD" if summ[k])
              + f"; {sum(chk.values())} machine-checked criteria "
              + "(" + ", ".join(f"{chk[k]} class-{k}" for k in "ABCD" if chk[k]) + ").")
        return 0

    print(f"{'BENCHMARK':<32} {'CLASS':<6} {'CHECKS':>6}  ESTABLISHES")
    for r in sorted(rows, key=lambda r: (r["class"], r["benchmark"])):
        syn = " (includes synthetic)" if r["synthetic"] else ""
        print(f"{r['benchmark']:<32} {r['class']:<6} {r['n_checks']:>6}  "
              f"{r['class_label']}{syn}")
    print()
    # Per-CHECK totals as well as per-benchmark, because a benchmark takes
    # its strongest class: class-B checks are real but invisible in the
    # per-benchmark tally whenever the same benchmark also has a class-A check.
    chk = {k: 0 for k in "ABCD"}
    for r in rows:
        for k, n in r["checks_by_class"].items():
            chk[k] += n
    for k in "ABCD":
        if summ[k] or chk[k]:
            print(f"  {k}  {CLASS_LABEL[k]:<28} "
                  f"{summ[k]:>2} benchmark(s)   {chk[k]:>3} check(s)")
    print(f"\n  {len(rows)} benchmarks, "
          f"{sum(r['n_checks'] for r in rows)} machine-checked criteria")
    print("\n  Only classes A and B establish that an implementation computes "
          "the right answer.\n  All classes execute fixed command sequences "
          "with no model in the loop.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
