#!/usr/bin/env python3
"""Check IGVFagent's complexity estimates against the paper's stated numbers.

Rosen et al. (2025), "Complexity and sequencing depth estimation", reports
median assigned barcodes and median pairwise Lincoln index for three
libraries. Each is recomputed here from the published IGVF barcode file and
compared to the number printed in the paper.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Scripts"))
from mpra_snakeflow_skill import (  # noqa: E402
    read_wide_barcode_file, complexity_metrics)
from mpralib_skill import run_outliers  # noqa: E402

# (label, IGVF file accession, paper's median observed, paper's median Lincoln)
CLAIMS = [
    ("8K-neurons", "IGVFFI0032GUOI", 1_444_480, 1_523_572),
    ("80K-neurons", "IGVFFI8345QIJJ", 5_459_247, 6_243_618),
]


# Outlier replicate-consistency, Figure 4B. The paper reports
# flagged-in-all / mean-flagged-per-replicate, NOT over the union.
# (label, accession, paper percentage)
CONSISTENCY_CLAIMS = [
    ("12K-cardiop", "IGVFFI4856GPLO", 7.8),
]


def _consistency(path) -> float:
    res = run_outliers(str(path), "global", times_zscore=3.0,
                       times_activity=5.0)
    flagged = res["_flagged"]
    reps = res["replicates"]
    inter = set(flagged[reps[0]])
    for r in reps[1:]:
        inter &= flagged[r]
    mean_per_rep = sum(len(flagged[r]) for r in reps) / len(reps)
    return 100.0 * len(inter) / mean_per_rep if mean_per_rep else 0.0


def main(data_dir: str, out_json: str) -> int:
    d = Path(data_dir)
    results = []
    exact = 0
    for label, acc, exp_obs, exp_lin in CLAIMS:
        path = d / f"{acc}.tsv.gz"
        if not path.exists():
            print(f"  {label}: input {path} missing — skipped")
            continue
        observed, shared, reps, n = read_wide_barcode_file(path)
        m = complexity_metrics(observed, shared, reps)
        got_obs = int(m["median_observed_barcodes"])
        got_lin = round(m["median_lincoln_index"])
        ok = (got_obs == exp_obs) and (got_lin == exp_lin)
        exact += int(ok)
        results.append({
            "dataset": label, "file": acc,
            "median_observed_barcodes": got_obs,
            "paper_median_observed_barcodes": exp_obs,
            "median_lincoln_index": got_lin,
            "paper_median_lincoln_index": exp_lin,
            "percent_missing": round(m["percent_missing"], 2),
            "matches_paper": int(ok),
        })
        mark = "EXACT" if ok else "DIFFERS"
        print(f"  {label:14s} observed {got_obs:>12,} (paper {exp_obs:>12,})  "
              f"Lincoln {got_lin:>12,} (paper {exp_lin:>12,})  {mark}")

    for label, acc, exp_pct in CONSISTENCY_CLAIMS:
        path = d / f"{acc}.tsv.gz"
        if not path.exists():
            print(f"  {label}: input {path} missing — skipped")
            continue
        got = round(_consistency(path), 1)
        ok = abs(got - exp_pct) < 0.05
        exact += int(ok)
        results.append({
            "dataset": label, "file": acc,
            "metric": "global outlier replicate consistency (%)",
            "value": got, "paper_value": exp_pct,
            "matches_paper": int(ok),
        })
        print(f"  {label:14s} outlier consistency {got:>5.1f}% "
              f"(paper {exp_pct:>5.1f}%)  {'EXACT' if ok else 'DIFFERS'}")

    summary = {"claims_checked": len(results),
               "claims_exact": exact,
               "claims": results}
    Path(out_json).write_text(json.dumps(summary, indent=2))
    print(f"\n  {exact}/{len(results)} published complexity claims reproduced exactly")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: claims.py <data_dir> <out_json>")
    sys.exit(main(*sys.argv[1:]))
