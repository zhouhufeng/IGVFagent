#!/usr/bin/env python3
"""Build a SPLiT-seq AnnData from Rosenberg 2018's main CNS atlas MAT-file.

GSE110823 ships each sample as a MATLAB v5 ``.mat`` (DGE = cells x genes sparse,
plus ``barcodes`` / ``genes`` / ``sample_type``). We read the headline
150k-nuclei CNS atlas, subsample to a memory-safe number of cells (deterministic
seed), and write an h5ad for ``igvfagent splitseq analyze``.
"""
import gzip
import shutil
from pathlib import Path

import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "Benchmarks/_data/rosenberg2018_splitseq"
TAR = D / "GSE110823_RAW.tar"
MEMBER = "GSM3017261_150000_CNS_nuclei.mat.gz"
N_SUB = 12000   # memory-safe subsample of the 156,049-nucleus atlas (dense scale fits ~4 GB)
SEED = 0

matgz = D / MEMBER
mat = D / MEMBER[:-3]
if not mat.is_file():
    if not matgz.is_file():
        import tarfile
        with tarfile.open(TAR) as t:
            t.extract(MEMBER, path=D)
    with gzip.open(matgz, "rb") as fi, open(mat, "wb") as fo:
        shutil.copyfileobj(fi, fo)

def _strlist(arr):
    return [str(x).strip() for x in np.asarray(arr).ravel()]

print(f"==> loadmat {mat.name}")
m = sio.loadmat(str(mat))
dge = m["DGE"]
genes = _strlist(m["genes"])
barcodes = _strlist(m["barcodes"])
sample_type = _strlist(m["sample_type"])
# Orient to cells x genes by matching axis lengths to the label vectors.
if dge.shape == (len(genes), len(barcodes)):
    dge = dge.T.tocsr()
elif dge.shape != (len(barcodes), len(genes)):
    # fall back: assume rows are the shorter labelled axis
    if dge.shape[0] == len(genes):
        dge = dge.T.tocsr()
print(f"    full atlas: {dge.shape[0]:,} nuclei x {dge.shape[1]:,} genes "
      f"(barcodes={len(barcodes)}, genes={len(genes)})")

# deterministic subsample of cells (guard against label/matrix length drift)
n = min(dge.shape[0], len(barcodes))
rng = np.random.default_rng(SEED)
idx = np.sort(rng.choice(n, size=min(N_SUB, n), replace=False))
X = sp.csr_matrix(dge[idx, :])
del dge

a = ad.AnnData(X=X)
a.obs_names = [barcodes[i] for i in idx]
a.var_names = genes[:X.shape[1]]
a.var_names_make_unique()
a.obs["sample_type"] = [sample_type[i] if i < len(sample_type) else "NA"
                        for i in idx]
# splitseq analyze uses a sample_id batch key for its (single-batch) Harmony step
a.obs["sample_id"] = a.obs["sample_type"].astype(str).values
a.obs["full_atlas_n"] = n
print(f"    subsample: {a.n_obs:,} nuclei x {a.n_vars:,} genes (seed {SEED})")
out = D / "rosenberg_cns_30k.h5ad"
a.write_h5ad(out, compression="gzip")
print(f"==> wrote {out.name}")
