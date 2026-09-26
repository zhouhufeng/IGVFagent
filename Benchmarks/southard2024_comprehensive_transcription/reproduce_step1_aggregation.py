#!/usr/bin/env python3
"""Faithful port of the authors' own notebook:
"Step 1 - Aggregation of cellranger outputs and guide thresholding for Hs27
experiment.ipynb" (norman-lab-msk/TFs_CRISPRa, commit 3637f77), run against
the RAW cellranger outputs deposited on Zenodo (10.5281/zenodo.15213597),
not the already-processed summary files used elsewhere in this benchmark.

This is a port, not a verbatim run, for two reasons stated plainly:

1. The notebook reads a single merged `cellranger aggr` output
   (EXPERIMENT/*/outs/count/filtered_feature_bc_matrix.h5). The Zenodo
   deposit only has the 16 per-lane `cellranger count` outputs (this is
   what the authors actually published), not the merged aggr output. This
   script concatenates the 16 lanes itself (barcode-suffixed per lane, like
   `cellranger aggr` does) instead of running `cellranger aggr` (not
   available in this environment). No depth-equalisation/subsampling is
   applied (`cellranger aggr`'s default), so absolute UMI counts here are
   not expected to be bit-identical to a real `cellranger aggr` run --
   documented, not hidden.
2. Cells 25-27 of the notebook (exploratory `sns.distplot` calls on
   individual guides, not used by any downstream cell) are omitted; they
   are dead ends for the pipeline, not part of the reproduction.

Every other step -- guide vs. GEX feature split, per-cell guide-UMI
stacking, >5 UMI threshold, dominant-guide assignment, gene mean-expression
filter -- reproduces the notebook's own cells line-for-line.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CR_DIR = ROOT / "Data" / "Benchmarks" / "southard2024_comprehensive_transcription" / "hs27_cellranger_count" / "cellranger_count"
OUT_DIR = ROOT / "Data" / "Benchmarks" / "southard2024_comprehensive_transcription" / "step1_reproduction"
REFERENCE_H5 = (ROOT / "Data" / "Benchmarks" / "southard2024_comprehensive_transcription"
                 / "fibroblast_CRISPRa_aggr_total_guide_umis.h5")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lanes = sorted(p.name for p in CR_DIR.iterdir() if p.is_dir())
    print(f"lanes found: {lanes}", file=sys.stderr)

    # --- cell 6: import + concatenate each lane, tagging obs['dataset'] ---
    adata_list = []
    for lane in lanes:
        h5 = CR_DIR / lane / "outs" / "filtered_feature_bc_matrix.h5"
        a = sc.read_10x_h5(h5, gex_only=False)
        a.var_names_make_unique()
        a.obs["dataset"] = lane
        adata_list.append(a)
        print(f"  {lane}: {a.shape}", file=sys.stderr)

    # cellranger aggr suffixes barcodes "-1", "-2", ... per input in the
    # order given; reproduce that instead of an arbitrary anndata.concat key.
    for i, a in enumerate(adata_list, start=1):
        a.obs_names = [bc.split("-")[0] + f"-{i}" for bc in a.obs_names]
    adata = sc.concat(adata_list, join="outer", index_unique=None, merge="same")
    print(f"concatenated: {adata.shape}", file=sys.stderr)

    # --- cells 13-14: split guide-capture vs gene-expression features ---
    adata_guides = adata[:, adata.var["feature_types"].isin(["CRISPR Guide Capture"])].copy()
    adata_gex = adata[:, adata.var["feature_types"].isin(["Gene Expression"])].copy()

    # --- cell 17: guide var index -> guide_identity / guide_target ---
    adata_guides.var = (adata_guides.var.rename(columns={"gene_ids": "guide_identity"})
                         .reset_index().set_index("guide_identity")
                         .rename(columns={"index": "guide_target"}))
    adata_guides.var["guide_target"] = adata_guides.var["guide_target"].map(lambda x: x.split("-")[0])

    # --- cell 19: chunked stack of nonzero guide UMIs per cell ---
    chunk_size = 1000
    stacked_list = []
    n = adata_guides.n_obs
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        df_chunk = adata_guides[start:end, :].to_df()
        m = df_chunk > 0
        masked_chunk = df_chunk.mask(~m)
        stacked_list.append(masked_chunk.stack())
        if start % 50000 == 0:
            print(f"  stacking {start}/{n}", file=sys.stderr)
    guide_umis = pd.concat(stacked_list)
    guide_umis.index.names = ["cell_barcode", "guide_identity"]

    # --- cells 30-31: write raw + threshold>5 guide UMI tables ---
    guide_umis.to_hdf(OUT_DIR / "aggr_total_guide_umis.h5", key="guide_umis", mode="w")
    filtered_guide_umis_raw = guide_umis[guide_umis > 5]
    filtered_guide_umis_raw.to_hdf(OUT_DIR / "aggr_thres5_guide_umis.h5", key="guide_umis", mode="w")

    # --- cells 36-42: dominant guide per cell + thresholded guide list ---
    filtered_guide_umis = (filtered_guide_umis_raw.reset_index()
                            .rename(columns={0: "guide_umi_count"}))
    filtered_guide_umis["guide_umi_count"] = filtered_guide_umis["guide_umi_count"].astype(int)

    cell_identities = (filtered_guide_umis.sort_values("guide_umi_count", ascending=False)
                       .groupby("cell_barcode").first())
    cell_identities["guide_umi_count"] = cell_identities["guide_umi_count"].astype(int)

    all_guides = filtered_guide_umis.groupby("cell_barcode")["guide_identity"].apply(lambda x: "|".join(x))
    all_guide_umi = filtered_guide_umis.groupby("cell_barcode")["guide_umi_count"].apply(
        lambda x: "|".join(x.astype(str)))
    num_cells = filtered_guide_umis.groupby("cell_barcode")["guide_identity"].count()

    cell_identities["thresholded_features"] = all_guides
    cell_identities["thresholded_guide_umi"] = all_guide_umi
    cell_identities["num_cells"] = num_cells

    # --- cells 43-46: merge into GEX, drop unassigned cells, gene filter ---
    adata_gex.var = (adata_gex.var.rename(columns={"gene_ids": "gene_id"})
                      .reset_index().set_index("gene_id").rename(columns={"index": "gene_name"}))
    adata_gex.obs["UMI_count"] = np.asarray(adata_gex.X.sum(axis=1)).flatten()
    adata_gex.obs = adata_gex.obs.merge(cell_identities, left_index=True, right_index=True, how="left")
    adata_gex = adata_gex[~adata_gex.obs["thresholded_features"].isnull()].copy()

    adata_gex.var["mean"] = np.asarray(adata_gex.X.mean(axis=0)).flatten()
    adata_gex = adata_gex[:, adata_gex.var["mean"] >= 0.05].copy()
    adata_gex.var["in_matrix"] = True

    adata_gex.write_h5ad(OUT_DIR / "cellranger_aggr_singlets_and_multiplets_5umi_thresh.h5ad")

    # --- independent cross-check against the already-verified Zenodo deposit ---
    guide_names_repro = set(guide_umis.index.get_level_values("guide_identity").astype(str))
    ref_guides = None
    if REFERENCE_H5.exists():
        import h5py
        with h5py.File(REFERENCE_H5, "r") as f:
            ref_guides = set(g.decode() for g in f["guide_umis/index_level1"][:])

    out = {
        "lanes_used": lanes,
        "n_lanes": len(lanes),
        "reproduction_method": "own concatenation of 16 per-lane cellranger count outputs (no cellranger aggr binary available); no depth-equalisation applied",
        "n_cells_concatenated_pre_filter": int(adata.n_obs),
        "n_guides_in_library": len(guide_names_repro),
        "n_cells_with_assigned_guide": int(adata_gex.n_obs),
        "n_genes_post_filter": int(adata_gex.n_vars),
        "cross_check_vs_zenodo_deposit": {
            "reference_file": str(REFERENCE_H5.relative_to(ROOT)) if REFERENCE_H5.exists() else None,
            "n_guides_in_reference": len(ref_guides) if ref_guides is not None else None,
            "guide_sets_identical": (guide_names_repro == ref_guides) if ref_guides is not None else None,
            "guide_sets_jaccard": (len(guide_names_repro & ref_guides) / len(guide_names_repro | ref_guides)
                                    if ref_guides is not None else None),
        },
    }
    (OUT_DIR / "step1_reproduction_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
