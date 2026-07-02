#!/usr/bin/env python3
"""Prepare the Travaglini 2020 Human Lung Cell Atlas h5ad for `sc-analyze`.

The CELLxGENE-distributed h5ad stores *normalized* log values in ``.X`` and
*raw counts* in ``.raw.X``; ``var_names`` are Ensembl gene IDs. IGVFagent's
``sc-analyze`` skill expects raw counts in ``.X`` and gene **symbols** as
``var_names`` (so the ``MT-`` mito prefix and canonical marker overlays
resolve). This step performs that standard CELLxGENE-to-Scanpy conversion and
nothing else — no filtering, no subsetting — so the downstream pipeline is the
thing under test.

Input : Benchmarks/_data/travaglini2020_lung/lung_atlas.h5ad   (CELLxGENE)
Output: Benchmarks/_data/travaglini2020_lung/lung_for_scanalyze.h5ad
"""
from pathlib import Path
import anndata as ad
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "Benchmarks/_data/travaglini2020_lung/lung_atlas.h5ad"
DST = ROOT / "Benchmarks/_data/travaglini2020_lung/lung_for_scanalyze.h5ad"

if not SRC.is_file():
    raise SystemExit(
        f"Missing {SRC}\n  Download first (see run.sh): the CELLxGENE asset for "
        f"DOI 10.1038/s41586-020-2922-4 (Krasnow Lab Human Lung Cell Atlas, 10X).")

print(f"==> Reading {SRC.name}")
a = ad.read_h5ad(SRC)

# raw counts -> X (raw spans the same 24,769 genes here)
counts = a.raw.X if a.raw is not None else a.X
if not sp.issparse(counts):
    counts = sp.csr_matrix(counts)

var = a.raw.var.copy() if a.raw is not None else a.var.copy()
symbols = var["feature_name"].astype(str).values if "feature_name" in var.columns \
    else var.index.astype(str).values

obs_keep = [c for c in ["cell_type", "free_annotation", "compartment",
                          "donor_id", "sample", "region", "location"]
            if c in a.obs.columns]

clean = ad.AnnData(X=counts.astype(np.float32),
                    obs=a.obs[obs_keep].copy())
clean.var_names = symbols
clean.var_names_make_unique()

# sanity: are these integer counts?
chk = clean.X[:500]
chk = chk.toarray() if sp.issparse(chk) else np.asarray(chk)
assert np.allclose(chk, np.round(chk)), "X is not integer counts after prep"

print(f"    {clean.n_obs:,} cells x {clean.n_vars:,} genes")
print(f"    author cell_type labels: {clean.obs['cell_type'].nunique()}")
print(f"==> Writing {DST.name}")
clean.write_h5ad(DST, compression="gzip")
print("done")
