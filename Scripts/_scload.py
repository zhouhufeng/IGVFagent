"""Shared single-cell loaders — internalized from the reproducibility benchmarks.

Public single-cell datasets ship in a handful of awkward formats that the
benchmark suite had to handle one-off (CELLxGENE h5ad with normalized X + Ensembl
ids; dense gene-by-cell TSV; MatrixMarket triplets; MATLAB v5 DGE). This module
folds those handlers into one reusable place so any skill (`sc-analyze`,
`splitseq`, `multiome`, `share`) — and future analyses — can ingest them via a
single import instead of re-deriving the parsing each time.

All loaders return an ``AnnData`` oriented **cells × genes** with a sparse CSR
``.X`` of raw counts, ready for the skills' QC → normalize → cluster chains.

Memory: the dense-TSV and MatrixMarket paths stream, so a 40k-cell × 25k-gene
matrix loads in well under 1 GB regardless of the on-disk representation.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.sparse as sp


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _anndata():
    import anndata as ad
    return ad


def dense_gene_by_cell_tsv(path, *, sep: str = "\t",
                           first_col_is_gene: bool = True):
    """Stream a dense **genes × cells** TSV → cells × genes sparse AnnData.

    Row 0 is a header of cell barcodes (optionally with a leading label cell);
    each subsequent row is ``<gene>  v1 v2 …``. Used by SHARE-seq (Ma 2020) and
    the multiome RNA matrices (Trevino GSE162170).
    """
    ad = _anndata()
    rows, genes = [], []
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split(sep)
        barcodes = header[1:] if first_col_is_gene else header
        for ln in fh:
            gid, _, rest = ln.partition(sep)
            genes.append(gid.strip())
            # Sparsify each gene row immediately — never materialize the full
            # dense genes×cells array (which would be many GB for a real atlas).
            rows.append(sp.csr_matrix(
                np.fromstring(rest, sep=sep, dtype=np.float32)))
    mat = sp.vstack(rows).T.tocsr()             # cells × genes, stays sparse
    a = ad.AnnData(X=mat)
    a.obs_names = [b.strip() for b in barcodes]
    a.var_names = genes
    a.var_names_make_unique()
    return a


def matrixmarket(mtx_path, *, barcodes_path=None, features_path=None,
                 features_are_rows: bool = True):
    """MatrixMarket (features × cells or cells × features) → cells × genes AnnData."""
    ad = _anndata()
    from scipy.io import mmread
    m = sp.csr_matrix(mmread(str(mtx_path)))
    # mmread gives rows × cols as stored; orient to cells × features
    if features_are_rows:
        m = m.T.tocsr()
    a = ad.AnnData(X=m)
    if barcodes_path:
        a.obs_names = [l.strip() for l in _open(barcodes_path)][: a.n_obs]
    if features_path:
        feats = [l.split("\t")[0].strip() for l in _open(features_path)][: a.n_vars]
        a.var_names = feats
        a.var_names_make_unique()
    return a


def matlab_dge(mat_path, *, dge_key: str = "DGE", genes_key: str = "genes",
               barcodes_key: str = "barcodes"):
    """MATLAB v5 ``.mat`` DGE (Rosenberg 2018 SPLiT-seq) → cells × genes AnnData.

    ``mat_path`` may be a ``.mat`` or ``.mat.gz`` (auto-decompressed to a temp
    sibling). Orientation is inferred by matching axis lengths to the label
    vectors.
    """
    ad = _anndata()
    import scipy.io as sio
    p = Path(mat_path)
    if p.suffix == ".gz":
        raw = p.with_suffix("")
        if not raw.is_file():
            import shutil
            with gzip.open(p, "rb") as fi, open(raw, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        p = raw
    m = sio.loadmat(str(p))
    dge = m[dge_key]
    genes = [str(x).strip() for x in np.asarray(m[genes_key]).ravel()]
    barcodes = [str(x).strip() for x in np.asarray(m[barcodes_key]).ravel()]
    if dge.shape == (len(genes), len(barcodes)):
        dge = dge.T
    elif dge.shape != (len(barcodes), len(genes)) and dge.shape[0] == len(genes):
        dge = dge.T
    a = ad.AnnData(X=sp.csr_matrix(dge))
    n = min(a.n_obs, len(barcodes))
    a.obs_names = barcodes[:a.n_obs] if len(barcodes) >= a.n_obs else \
        barcodes + [f"cell_{i}" for i in range(a.n_obs - len(barcodes))]
    a.var_names = genes[:a.n_vars]
    a.var_names_make_unique()
    return a


def cellxgene_h5ad(path):
    """CELLxGENE h5ad → cells × genes AnnData with **raw counts** in X and gene
    **symbols** as var_names (Travaglini 2020 lung). Reads raw counts from
    ``.raw`` and ``feature_name`` symbols, dropping the normalized ``.X``.
    """
    ad = _anndata()
    src = ad.read_h5ad(str(path))
    counts = src.raw.X if src.raw is not None else src.X
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    var = src.raw.var if src.raw is not None else src.var
    symbols = (var["feature_name"].astype(str).values
               if "feature_name" in var.columns else var.index.astype(str).values)
    obs_keep = [c for c in ("cell_type", "free_annotation", "compartment")
                if c in src.obs.columns]
    a = ad.AnnData(X=counts.astype(np.float32), obs=src.obs[obs_keep].copy())
    a.var_names = symbols
    a.var_names_make_unique()
    return a


def subsample_cells(adata, n: int, *, seed: int = 0):
    """Deterministic cell subsample (memory-safe cap for laptop-scale runs)."""
    if adata.n_obs <= n:
        return adata
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False))
    return adata[idx].copy()


def attach_labels(adata, tsv_path, *, key_col: str, label_col: str,
                  obs_name: str = "cell_type", sep: str = "\t"):
    """Join an external per-cell label table (e.g. author cell types) into obs."""
    import pandas as pd
    df = pd.read_csv(tsv_path, sep=sep)
    df = df.set_index(df[key_col].astype(str))
    adata.obs[obs_name] = df[label_col].reindex(adata.obs_names).values
    return adata


__all__ = ["dense_gene_by_cell_tsv", "matrixmarket", "matlab_dge",
           "cellxgene_h5ad", "subsample_cells", "attach_labels"]
