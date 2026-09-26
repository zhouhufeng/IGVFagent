# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/buenrostrolab/FigR @ 094f5aa036aa
# (R/utils.R) for ma2020_shareseq. Unreviewed; provenance in Scripts/ported/registry.json.
#!/usr/bin/env python3
"""Chromatin potential (Ma et al., Cell 2020, Fig. 5H; STAR Methods).

No code for this analysis was released; this follows the STAR Methods text:

  "we first smoothed DORC scores (chromatin space) and corresponding gene
   expression (RNA space) over a k-nearest neighbor graph (k-NN, k = 50),
   calculated using normalized ATAC topics from cisTopic. Next, we calculated
   another k-NN (k = 10), between the smoothed chromatin profile of a given
   cell (C_atac,i), and the smoothed gene expression profile of each cell
   (C_rna,j). We then calculated the distance (D_i,j) between the C_atac,i and
   the average of C_rna,j in chromatin space."

Choices where the text is silent (recorded in summary.json):
  * ATAC low-dimensional space: LSI (TF-IDF, truncated SVD, first component
    dropped, L2-normalised) instead of cisTopic topics.
  * Chromatin and RNA profiles are z-scored per DORC gene before the
    cross-modal k-NN, so both live on one scale.
  * The arrow of cell i points from its position to the mean position of its
    10 RNA-neighbours in the ATAC embedding; its length is D_i / max(D).

The flow test counts, per source cell type, where the arrow's RNA neighbours
sit: a progenitor->differentiated flow means progenitor cells' neighbours are
more often differentiated than differentiated cells' neighbours are
progenitors.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

def _read_mtx(path):
    """MatrixMarket (.mtx/.mtx.gz, as deposited on GEO) or a scipy .npz cache."""
    path = str(path)
    if path.endswith(".npz"):
        return sp.load_npz(path)
    return scipy.io.mmread(gzip.open(path) if path.endswith(".gz") else path)



def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)



def read_rna_rows(path, genes):
    """Genes x cells dense TSV, keeping only ``genes`` (streams the file)."""
    want = set(genes)
    parts = [c[c.index.isin(want)] for c in pd.read_csv(path, sep="\t", index_col=0, chunksize=2000)]
    r = pd.concat(parts)
    r.columns = [c.replace(",", ".") for c in r.columns]
    return r

def lsi(counts_cells_x_peaks, n_comp=30, seed=0):
    X = sp.csr_matrix(counts_cells_x_peaks, dtype=np.float64)
    X.data[:] = 1.0  # binarise
    tf = sp.diags(1.0 / np.maximum(np.asarray(X.sum(1)).ravel(), 1)) @ X
    idf = np.log1p(X.shape[0] / np.maximum(np.asarray(X.sum(0)).ravel(), 1))
    tfidf = (tf @ sp.diags(idf)).tocsr()
    tfidf.data = np.log1p(tfidf.data * 1e4)
    Z = TruncatedSVD(n_components=n_comp + 1, random_state=seed).fit_transform(tfidf)[:, 1:]
    Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    return Z


def smooth(mat_cells_x_genes, nn_idx):
    return mat_cells_x_genes[nn_idx].mean(axis=1)


def zscore(m):
    sd = m.std(axis=0)
    sd[sd == 0] = 1.0
    return (m - m.mean(axis=0)) / sd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--atac-mtx", required=True, help="peaks x cells MatrixMarket")
    ap.add_argument("--atac-barcodes", required=True)
    ap.add_argument("--dorc-scores", required=True, help="genes x cells TSV from shareseq_dorc dorc-scores")
    ap.add_argument("--dorc-genes", help="restrict to these DORC genes (one per line)")
    ap.add_argument("--rna", required=True, help="genes x cells dense TSV (.gz)")
    ap.add_argument("--celltypes", required=True, help="TSV with rna.bc and celltype columns")
    ap.add_argument("--types", required=True, help="comma-separated cell types to include")
    ap.add_argument("--progenitor", required=True, help="comma-separated progenitor types (e.g. TAC-1,TAC-2)")
    ap.add_argument("--differentiated", required=True, help="comma-separated differentiated types")
    ap.add_argument("--pseudotime", help="TSV cell -> pseudotime (e.g. from shareseq-dorc-residuals) for the direction test")
    ap.add_argument("--k-smooth", type=int, default=50)
    ap.add_argument("--k-cross", type=int, default=10)
    ap.add_argument("--n-comp", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ct = pd.read_csv(a.celltypes, sep="\t").drop_duplicates("rna.bc").set_index("rna.bc")["celltype"]
    types = a.types.split(",")
    with _open(a.atac_barcodes) as fh:
        abcs = [l.strip() for l in fh if l.strip()]
    dorc = pd.read_csv(a.dorc_scores, sep="\t", index_col=0)
    if a.dorc_genes:
        with _open(a.dorc_genes) as fh:
            keep = [g.strip() for g in fh if g.strip()]
        dorc = dorc.loc[[g for g in keep if g in dorc.index]]
    rna = read_rna_rows(a.rna, dorc.index)
    genes = [g for g in dorc.index if g in rna.index]
    cells = [b for b in abcs if b in ct.index and ct[b] in types and b in rna.columns and b in dorc.columns]
    print(f"cells {len(cells):,}  DORC genes {len(genes)}")

    A = _read_mtx(a.atac_mtx).tocsc()
    col = {b: i for i, b in enumerate(abcs)}
    A = A[:, [col[b] for b in cells]].T.tocsr()
    Z = lsi(A, a.n_comp, a.seed)

    nn = NearestNeighbors(n_neighbors=a.k_smooth).fit(Z)
    idx = nn.kneighbors(Z, return_distance=False)

    R = rna[cells]
    R = R.div(R.sum(axis=0).replace(0, 1), axis=1) * R.sum(axis=0).mean()
    R = np.log1p(R.loc[genes].to_numpy().T)
    C = dorc.loc[genes, cells].to_numpy().T
    Cs = zscore(smooth(C, idx))
    Rs = zscore(smooth(R, idx))

    nn2 = NearestNeighbors(n_neighbors=a.k_cross).fit(Rs)
    d2, j2 = nn2.kneighbors(Cs)
    target = Z[j2].mean(axis=1)
    D = np.linalg.norm(target - Z, axis=1)
    length = D / D.max()

    lab = np.asarray([ct[b] for b in cells])
    prog, diff = set(a.progenitor.split(",")), set(a.differentiated.split(","))
    nb = lab[j2]
    is_prog = np.isin(lab, list(prog))
    is_diff = np.isin(lab, list(diff))
    nb_diff = np.isin(nb, list(diff)).mean(axis=1)
    nb_prog = np.isin(nb, list(prog)).mean(axis=1)
    fwd = {}
    if a.pseudotime:
        pt = pd.read_csv(a.pseudotime, sep="\t", index_col=0).iloc[:, 0].reindex(cells).to_numpy(float)
        ok = np.isfinite(pt)
        delta = np.nanmean(pt[j2], axis=1) - pt
        # baseline: same kNN done RNA->RNA (self excluded), i.e. no chromatin lead
        _, j_rr = NearestNeighbors(n_neighbors=a.k_cross + 1).fit(Rs).kneighbors(Rs)
        base = np.nanmean(pt[j_rr[:, 1:]], axis=1) - pt
        fwd = {
            "forward_fraction_all": float(np.mean(delta[ok] > 0)),
            "forward_fraction_progenitor": float(np.mean(delta[ok & is_prog] > 0)),
            "mean_pseudotime_shift_all": float(np.nanmean(delta)),
            "mean_pseudotime_shift_progenitor": float(np.nanmean(delta[is_prog])),
            "baseline_rna_rna_forward_fraction_all": float(np.mean(base[ok] > 0)),
            "baseline_rna_rna_mean_shift_all": float(np.nanmean(base)),
        }
    per_type = {t: {"n": int((lab == t).sum()),
                    "mean_frac_neighbours_differentiated": float(nb_diff[lab == t].mean()),
                    "mean_frac_neighbours_progenitor": float(nb_prog[lab == t].mean()),
                    "mean_arrow_length": float(length[lab == t].mean())} for t in types if (lab == t).any()}
    summary = {
        "n_cells": len(cells), "n_dorc_genes": len(genes),
        "embedding": f"LSI ({a.n_comp} comps, first dropped) in place of cisTopic",
        "k_smooth": a.k_smooth, "k_cross": a.k_cross,
        "prog_to_diff_frac": float(nb_diff[is_prog].mean()),
        "diff_to_prog_frac": float(nb_prog[is_diff].mean()),
        "flow_prog_to_diff_ratio": float(nb_diff[is_prog].mean() / max(nb_prog[is_diff].mean(), 1e-9)),
        **fwd,
        "per_type": per_type,
    }
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell": cells, "celltype": lab, "arrow_length": length,
                  "frac_nb_differentiated": nb_diff, "frac_nb_progenitor": nb_prog}
                 ).to_csv(out / "chromatin_potential_cells.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "per_type"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
