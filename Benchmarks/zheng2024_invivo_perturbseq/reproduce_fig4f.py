#!/usr/bin/env python3
"""Reproduce Zheng et al. 2024 (Cell) Fig 4F -- Foxg1 knockdown cell-type
proportion shifts -- from the paper's own deposited Seurat object metadata
(GSE249416_Perturb_all.qs, extracted by extract_metadata.R).

Paper's own method (STAR Methods, "Perturbation-associated analysis"):
  "Statistics for these pairwise composition comparisons were computed using
   the propeller.ttest function from speckle (R package v0.99.7). The
   proportions were first transformed using arcsin square root
   transformation, and the batch (10x channel) was additionally considered
   as another fixed effect to the linear models."

`speckle` isn't installed here, so this is a from-scratch Python
re-implementation of that recipe (arcsin-sqrt proportions + OLS with channel
as a fixed effect, gRNA-vs-NonTarget2 contrast) -- not a call into the
original R package. It is expected to approximate direction and
order-of-magnitude, not reproduce the paper's exact p-values.

Headline claim being checked (Results, "Distinct cell type proportion changes
by transcription factor perturbation in vivo"):
  "Cells perturbed by Foxg1-gRNA1 had a 9.9-fold reduction in L6-IT neurons
   (FDR=9.6x10^-6), accompanied by a 2.0-fold increase of upper layer
   projection neurons (FDR=0.014)."

Usage:
    python3 reproduce_fig4f.py <metadata_csv> <out_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


def arcsin_sqrt(p):
    return np.arcsin(np.sqrt(np.clip(p, 0, 1)))


def cell_type_proportions(meta, group_col="gRNA", batch_col="channel", type_col="cell_type"):
    tab = meta.groupby([group_col, batch_col, type_col]).size().rename("n").reset_index()
    totals = meta.groupby([group_col, batch_col]).size().rename("total").reset_index()
    tab = tab.merge(totals, on=[group_col, batch_col])
    tab["prop"] = tab["n"] / tab["total"]
    tab["asin_sqrt_prop"] = arcsin_sqrt(tab["prop"])
    return tab


def propeller_like_ttest(tab, test_group, ref_group, cell_type,
                          group_col="gRNA", batch_col="channel", type_col="cell_type"):
    sub = tab[(tab[type_col] == cell_type) & (tab[group_col].isin([test_group, ref_group]))].copy()
    if sub[group_col].nunique() < 2 or len(sub) < 4:
        return None
    sub["is_test"] = (sub[group_col] == test_group).astype(int)
    model = smf.ols(f"asin_sqrt_prop ~ is_test + C({batch_col})", data=sub).fit()
    if "is_test" not in model.params.index:
        return None
    prop_test = sub.loc[sub["is_test"] == 1, "prop"].mean()
    prop_ref = sub.loc[sub["is_test"] == 0, "prop"].mean()
    fold = (prop_test / prop_ref) if prop_ref > 0 else np.nan
    return dict(cell_type=cell_type, test_group=test_group, ref_group=ref_group,
                beta_asin_sqrt=model.params["is_test"], pval=model.pvalues["is_test"],
                prop_test=prop_test, prop_ref=prop_ref, fold=fold, n=len(sub))


def main():
    meta_path, out_dir = sys.argv[1], sys.argv[2]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(meta_path, index_col=0)
    n_cells_total = len(df)
    n_channels = df["orig.ident"].nunique()

    # QC-passed singlets only (paper's Methods: "cells assigned to that
    # gRNA/perturbation identity ... excluding cell barcodes that occur in
    # multiple 10x channels"; Keep==1 is the object's own post-QC flag).
    singlet = df["assignment"].notna() & (~df["assignment"].astype(str).str.contains(","))
    qc = (df["Keep"] == 1) & (df["lowQC"] == False)  # noqa: E712
    sub = df[singlet & qc].rename(columns={"assignment": "gRNA", "orig.ident": "channel", "CellType": "cell_type"}).copy()
    # paper excludes cell types with <200 cells from the composition test
    keep_types = sub["cell_type"].value_counts()
    keep_types = keep_types[keep_types >= 200].index
    sub = sub[sub["cell_type"].isin(keep_types)]
    n_cells_tested = len(sub)

    tab = cell_type_proportions(sub)
    tab.to_csv(out_dir / "proportion_table.csv", index=False)

    ref = "NonTarget2"
    results = []
    for g in sorted(sub["gRNA"].unique()):
        if g == ref:
            continue
        for ct in sorted(sub["cell_type"].unique()):
            r = propeller_like_ttest(tab, g, ref, ct)
            if r:
                results.append(r)
    res = pd.DataFrame(results)
    res["FDR"] = multipletests(res["pval"], method="fdr_bh")[1]
    res = res.sort_values("FDR")
    res.to_csv(out_dir / "cell_type_proportion_tests.csv", index=False)

    def lookup(gRNA, cell_type):
        row = res[(res["test_group"] == gRNA) & (res["cell_type"] == cell_type)]
        return row.iloc[0] if len(row) else None

    l6it = lookup("Foxg1_1", "Excit_L6IT")
    upper = lookup("Foxg1_1", "Excit_Upper")

    summary = {
        "paper": {
            "doi": "10.1016/j.cell.2024.04.050",
            "pmid": "38772369",
            "figure": "Fig 4F / Results text",
            "quote": "Cells perturbed by Foxg1-gRNA1 had a 9.9-fold reduction in L6-IT "
                     "neurons (FDR=9.6x10-6), accompanied by a 2.0-fold increase of upper "
                     "layer projection neurons (FDR=0.014).",
        },
        "n_cells_total": int(n_cells_total),
        "n_channels": int(n_channels),
        "n_cells_qc_singlet_tested": int(n_cells_tested),
        "foxg1_g1_l6it_reduction_rate": float(1 / l6it["fold"]) if l6it is not None else None,
        "foxg1_g1_l6it_reduction_pval": float(l6it["pval"]) if l6it is not None else None,
        "foxg1_g1_l6it_reduction_fdr": float(l6it["FDR"]) if l6it is not None else None,
        "foxg1_g1_upper_increase_rate": float(upper["fold"]) if upper is not None else None,
        "foxg1_g1_upper_increase_pval": float(upper["pval"]) if upper is not None else None,
        "foxg1_g1_upper_increase_fdr": float(upper["FDR"]) if upper is not None else None,
        "method_note": "arcsin-sqrt OLS with channel as a fixed effect, gRNA vs NonTarget2 "
                       "contrast -- a from-scratch re-implementation of the paper's described "
                       "propeller.ttest recipe (speckle package), not a call into speckle "
                       "itself. Expect directional + order-of-magnitude agreement, not exact "
                       "p-value/FDR reproduction.",
    }
    (out_dir / "concordance_metrics.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
