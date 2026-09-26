"""Build the interventionGraph-shaped input table the LLCB model expects:
one row per sample, one column per one of the paper's 84 KO'd genes (VST
expression, donor + a handful of expression PCs regressed out), plus
`donor` and `intervention` columns. Mirrors Weinstock's own
`regress_out_covariates` + `zero_intervened_nodes` steps in
`RNAseq-perturbation-CD4-pipeline/R/functions.R` and `data.jl`, with one
substitution: they regress out a curated `covariates.tsv` of unwanted-
variation covariates; we don't have that file, so we estimate our own
surrogate covariates as the top expression PCs computed from the AAVS1
control samples (the standard approach when a curated covariate table
isn't available -- and the approach their own README recommends verbatim
for anyone else applying this method: "we recommend estimating expression
PCs and regressing out those which correspond to unwanted global sources
of variation").

Input:  Data/Weinstock2024/processed/{vst.csv.gz,sample_meta.csv}
        Data/Weinstock2024/raw/gene_symbol_to_ensembl.json
Output: Data/Weinstock2024/processed/network_input.csv
        Data/Weinstock2024/processed/pc_variance_explained.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROC = ROOT / "Data/Weinstock2024/processed"
RAW = ROOT / "Data/Weinstock2024/raw"

N_PCS = 10  # number of control-derived expression PCs to regress out
N_HVG = 2000  # genes used to estimate the unwanted-variation PCs


def main() -> None:
    vst = pd.read_csv(PROC / "vst.csv.gz", index_col=0)  # genes x samples
    meta = pd.read_csv(PROC / "sample_meta.csv", index_col="sample", dtype={"donor": str})
    gene_map = json.loads((RAW / "gene_symbol_to_ensembl.json").read_text())

    vst = vst[meta.index]  # order columns to match meta

    # --- 1. Estimate unwanted-variation PCs from the AAVS1 control samples ---
    control_samples = meta.index[meta.is_control]
    hvg = vst.loc[:, control_samples].var(axis=1).sort_values(ascending=False).index[:N_HVG]
    ctrl_mat = vst.loc[hvg, control_samples].T.values  # controls x genes
    ctrl_mat_centered = ctrl_mat - ctrl_mat.mean(axis=0, keepdims=True)

    u, s, vt = np.linalg.svd(ctrl_mat_centered, full_matrices=False)
    var_explained = (s**2) / np.sum(s**2)
    print("top PC variance explained (controls, HVG):", np.round(var_explained[:N_PCS], 3))

    loadings = vt[:N_PCS].T  # genes(hvg) x N_PCS
    all_mat = vst.loc[hvg, :].T.values  # all_samples x genes(hvg)
    all_mat_centered = all_mat - ctrl_mat.mean(axis=0, keepdims=True)
    sample_pcs = all_mat_centered @ loadings  # all_samples x N_PCS

    json.dump(
        {"variance_explained": var_explained[:N_PCS].tolist(), "n_hvg": N_HVG},
        open(PROC / "pc_variance_explained.json", "w"),
        indent=1,
    )

    # --- 2. Regress donor + these PCs out of every one of the 84 target genes ---
    donor_dummies = pd.get_dummies(meta["donor"], prefix="donor", drop_first=True).astype(float)
    design = np.column_stack(
        [np.ones(len(meta)), donor_dummies.values, sample_pcs]
    )  # samples x (1 + n_donor_dummies + N_PCS)

    target_genes = sorted(gene_map.keys())
    residualized = {}
    for gene in target_genes:
        ens_ids = [e for e in gene_map[gene] if e in vst.index]
        if not ens_ids:
            print(f"WARNING: {gene} not found in VST matrix, skipping")
            continue
        y = vst.loc[ens_ids].mean(axis=0).values  # collapse multi-mapped ensembl ids
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        resid = y - design @ beta
        residualized[gene] = resid + y.mean()  # add back the mean, matching functions.R

    expr = pd.DataFrame(residualized, index=meta.index)
    print(f"residualized expression matrix: {expr.shape}")

    # --- 3. Zero out each sample's own knocked-out gene (sentinel used downstream) ---
    for gene in expr.columns:
        rows = meta.index[meta["ko"] == gene]
        expr.loc[rows, gene] = 0.0

    # --- 4. Assemble the final interventionGraph-shaped table ---
    out = expr.copy()
    out.insert(0, "intervention", meta["ko"])
    out.insert(0, "donor", meta["donor"])
    out.to_csv(PROC / "network_input.csv", index=False)
    print(f"wrote {PROC / 'network_input.csv'}  shape={out.shape}")
    print(out["intervention"].value_counts().head())


if __name__ == "__main__":
    main()
