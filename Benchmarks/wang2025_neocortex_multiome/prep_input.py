#!/usr/bin/env python3
"""Prepare the Wang 2025 developing-neocortex multiome (CELLxGENE) for sc-analyze.

The CELLxGENE deposit is the RNA side of the 10x-multiome atlas (232,328 nuclei;
the ATAC peaks live only on the token-gated Dryad copy). We memory-safely
subsample, pull raw counts + gene symbols (the internalized CELLxGENE handling),
and carry the author cell-type labels for concordance scoring.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import anndata as ad

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "Benchmarks/_data/wang2025_neocortex_multiome"
SRC = D / "wang_multiome.h5ad"
N_SUB = 20000
SEED = 0

if not SRC.is_file():
    sys.exit(f"missing {SRC} — download the CELLxGENE h5ad first (see run.sh)")

print("==> opening backed + subsampling")
a = ad.read_h5ad(SRC, backed="r")
rng = np.random.default_rng(SEED)
idx = np.sort(rng.choice(a.n_obs, size=min(N_SUB, a.n_obs), replace=False))
sub = a[idx].to_memory()              # loads only the subsample (with .raw)
a.file.close()

# raw counts -> X, gene symbols -> var_names (internalized CELLxGENE handling)
counts = sub.raw.X if sub.raw is not None else sub.X
if not sp.issparse(counts):
    counts = sp.csr_matrix(counts)
var = sub.raw.var if sub.raw is not None else sub.var
symbols = (var["feature_name"].astype(str).values
           if "feature_name" in var.columns else var.index.astype(str).values)

obs_keep = [c for c in ("cell_type", "Class", "Subclass", "Type_updated",
                        "Region", "Sample_ID") if c in sub.obs.columns]
clean = ad.AnnData(X=counts.astype(np.float32), obs=sub.obs[obs_keep].copy())
clean.var_names = symbols[: counts.shape[1]]
clean.var_names_make_unique()
clean.obs["sample_id"] = clean.obs.get("Sample_ID",
                                       np.array(["wang2025"] * clean.n_obs)).astype(str).values

out = D / "wang_rna_20k.h5ad"
clean.write_h5ad(out, compression="gzip")
print(f"==> wrote {out.name}: {clean.n_obs:,} cells × {clean.n_vars:,} genes; "
      f"{clean.obs['cell_type'].nunique()} author cell types "
      f"(full atlas 232,328 nuclei)")
