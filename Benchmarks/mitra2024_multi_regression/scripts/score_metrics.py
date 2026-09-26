#!/usr/bin/env python3
"""Score the Mitra 2024 (SCARlink) BMMC peak->gene reproduction.

Reads the peak2gene TSV written by `igvfagent multiome peak2gene` on the
paper's exact BMMC donor/cell-type subset (site1_donor1/2/3,
site2_donor1/4/5, site3_donor6/7/10, site4_donor9; HSC, MK/E prog,
Proerythroblast, Erythroblast, Normoblast) and writes
mitra2024_multi_regression_concordance_metrics.json.

These are IGVFagent's OWN measured quantities from a generic per-peak
Pearson-correlation test, not a reimplementation of SCARlink's Poisson
regression model. They are NOT directly comparable to the paper's
model-specific gene counts (e.g. "1,655 genes for BMMC") — see the
benchmark README's Honest caveats section.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
LABEL = "mitra2024_multi_regression"

DONORS = {
    "site1_donor1_multiome", "site1_donor2_multiome", "site1_donor3_multiome",
    "site2_donor1_multiome", "site2_donor4_multiome", "site2_donor5_multiome",
    "site3_donor10_multiome", "site3_donor6_multiome", "site3_donor7_multiome",
    "site4_donor9_multiome",
}
CELL_TYPES = {"HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"}


def latest_tsv() -> Path:
    hits = sorted((ROOT / "Docs/Multiome10x").glob(f"2*{LABEL}_peak2gene.tsv"), reverse=True)
    if not hits:
        sys.exit(f"No peak2gene TSV for {LABEL} — run Benchmarks/{LABEL}/run.sh first.")
    return hits[0]


TSV = latest_tsv()
df = pd.read_csv(TSV, sep="\t")

sig = df[df["padj"] < 0.05]
pos = sig[sig["correlation"] > 0]

metrics = {
    "source_tsv": str(TSV.relative_to(ROOT)),
    "n_donors_matched": len(DONORS),
    "n_cell_types_matched": len(CELL_TYPES),
    "genes_tested": int(df["gene"].nunique()),
    "peaks_tested": int(df["peak"].nunique()),
    "candidate_pairs": int(len(df)),
    "significant_pairs": int(len(sig)),
    "significant_positive_pairs": int(len(pos)),
    "positive_fraction_of_significant": round(len(pos) / max(1, len(sig)), 4),
    "genes_with_significant_link": int(sig["gene"].nunique()),
    "genes_with_significant_positive_link": int(pos["gene"].nunique()),
    "median_abs_correlation_all_pairs": round(float(df["correlation"].abs().median()), 4),
}
out = ROOT / "Docs/Multiome10x" / f"{LABEL}_concordance_metrics.json"
out.write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics, indent=2))
print("Wrote", out.relative_to(ROOT))
