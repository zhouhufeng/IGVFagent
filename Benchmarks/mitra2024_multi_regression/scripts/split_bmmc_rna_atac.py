#!/usr/bin/env python
"""Split the GSE194122 NeurIPS-2021 multiome BMMC h5ad into RNA/ATAC
h5ads restricted to the donors and cell types the Mitra 2024 (SCARlink)
paper says it used, for `igvfagent multiome peak2gene`.

Paper (Data availability): "We used BMMC samples labeled as site1_donor1,
site1_donor2, site1_donor3, site2_donor1, site2_donor4, site2_donor5,
site3_donor10, site3_donor6, site3_donor7 and site4_donor9 and the cell
types HSC, MK/E progenitor, proerythroblast, erythroblast and normoblast."
"""
import anndata as ad

SRC = "Data/IGVF/GEO/Downloads/GSE194122/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
OUT_DIR = "Benchmarks/_data/mitra2024_multi_regression"

DONORS = {
    "site1_donor1_multiome", "site1_donor2_multiome", "site1_donor3_multiome",
    "site2_donor1_multiome", "site2_donor4_multiome", "site2_donor5_multiome",
    "site3_donor10_multiome", "site3_donor6_multiome", "site3_donor7_multiome",
    "site4_donor9_multiome",
}
CELL_TYPES = {"HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"}

import os
os.makedirs(OUT_DIR, exist_ok=True)

print("Reading", SRC)
adata = ad.read_h5ad(SRC)
print("Full shape:", adata.shape)

mask = adata.obs["Samplename"].isin(DONORS) & adata.obs["cell_type"].isin(CELL_TYPES)
sub = adata[mask].copy()
print("Filtered to paper's donors/cell types:", sub.shape)
print(sub.obs["Samplename"].value_counts())
print(sub.obs["cell_type"].value_counts())

gex = sub[:, sub.var["feature_types"] == "GEX"].copy()
atac = sub[:, sub.var["feature_types"] == "ATAC"].copy()
print("RNA (GEX):", gex.shape, " ATAC:", atac.shape)

# peak2gene's parse_peak wants "chr:start-end"; NeurIPS var_names are "chr-start-end".
new_names = []
for n in atac.var_names:
    chrom, start, end = n.rsplit("-", 2)
    new_names.append(f"{chrom}:{start}-{end}")
atac.var_names = new_names

gex.write_h5ad(f"{OUT_DIR}/mitra2024_bmmc_rna.h5ad")
atac.write_h5ad(f"{OUT_DIR}/mitra2024_bmmc_atac.h5ad")
print("Wrote", f"{OUT_DIR}/mitra2024_bmmc_rna.h5ad", "and", f"{OUT_DIR}/mitra2024_bmmc_atac.h5ad")
