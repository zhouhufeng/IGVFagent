# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/buenrostrolab/FigR @ 094f5aa036aa
# (R/DORCs.R) for ma2020_shareseq. Unreviewed; provenance in Scripts/ported/registry.json.
#!/usr/bin/env python3
"""Chromatin-minus-expression residuals of DORC genes over pseudotime
(Ma et al., Cell 2020, Fig. 4B-C; STAR Methods "Pseudotime inference",
"Residual analysis"). No code was released; this follows the Methods text:

  * Pseudotime from scATAC: TAC, IRS and hair-shaft cells; 10 normalised
    ATAC components -> Palantir diffusion maps (n_components=10), then
    Palantir (num_waypoints=1000, knn=30). cisTopic topics are replaced by
    LSI components (TF-IDF + SVD, first component dropped).
  * Lineages: the paper picked cells above manually chosen branch-probability
    cut-offs; here a lineage is the progenitor cells plus one terminal type
    (by the authors' cell-type labels), which is the reproducible analogue.
  * "Both DORC scores and gene expression were smoothed over pseudotime with
    local polynomial regression fitting (loess) separately, then min-max
    normalized. The residual for each gene was calculated by subtracting
    normalized gene expression from normalized DORC scores."
The reported number is the share of DORC genes whose mean residual is
positive (chromatin ahead of expression), per lineage and pooled.
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
from statsmodels.nonparametric.smoothers_lowess import lowess

def _read_mtx(path):
    """MatrixMarket (.mtx/.mtx.gz, as deposited on GEO) or a scipy .npz cache."""
    path = str(path)
    if path.endswith(".npz"):
        return sp.load_npz(path)
    return scipy.io.mmread(gzip.open(path) if path.endswith(".gz") else path)


sys.path.insert(0, str(Path(__file__).resolve().parent))
from shareseq_chromatin_potential import lsi, read_rna_rows  # noqa: E402


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def minmax(v):
    lo, hi = np.nanmin(v), np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else np.zeros_like(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--atac-mtx", required=True)
    ap.add_argument("--atac-barcodes", required=True)
    ap.add_argument("--dorc-scores", required=True)
    ap.add_argument("--dorc-genes")
    ap.add_argument("--rna", required=True)
    ap.add_argument("--celltypes", required=True)
    ap.add_argument("--progenitor", required=True, help="comma-separated progenitor types")
    ap.add_argument("--terminals", required=True, help="comma-separated terminal types, one lineage each")
    ap.add_argument("--n-comp", type=int, default=10)
    ap.add_argument("--frac", type=float, default=0.3, help="lowess span")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    import palantir
    import scanpy as sc

    ct = pd.read_csv(a.celltypes, sep="\t").drop_duplicates("rna.bc").set_index("rna.bc")["celltype"]
    prog = a.progenitor.split(",")
    terms = a.terminals.split(",")
    with _open(a.atac_barcodes) as fh:
        abcs = [l.strip() for l in fh if l.strip()]
    dorc = pd.read_csv(a.dorc_scores, sep="\t", index_col=0)
    if a.dorc_genes:
        with _open(a.dorc_genes) as fh:
            keep = [g.strip() for g in fh if g.strip()]
        dorc = dorc.loc[[g for g in keep if g in dorc.index]]
    rna = read_rna_rows(a.rna, dorc.index)
    genes = [g for g in dorc.index if g in rna.index]
    cells = [b for b in abcs if b in ct.index and ct[b] in prog + terms and b in rna.columns]

    A = _read_mtx(a.atac_mtx).tocsc()
    col = {b: i for i, b in enumerate(abcs)}
    Z = lsi(A[:, [col[b] for b in cells]].T.tocsr(), a.n_comp, a.seed)
    pca = pd.DataFrame(Z, index=cells)
    dm = palantir.utils.run_diffusion_maps(pca, n_components=10)
    ms = palantir.utils.determine_multiscale_space(dm)
    lab = pd.Series([ct[b] for b in cells], index=cells)
    # root (the paper does not say how it was picked): the progenitor cell
    # farthest, in the multiscale diffusion space, from its nearest terminal
    # cell-type centroid
    cent = np.stack([ms[lab.to_numpy() == t].mean(axis=0) for t in terms])
    P = ms[lab.isin(prog).to_numpy()]
    dmin = np.min(np.linalg.norm(P.to_numpy()[:, None, :] - cent[None], axis=2), axis=1)
    root = P.index[int(np.argmax(dmin))]
    np.random.seed(a.seed)
    pr = palantir.core.run_palantir(ms, root, num_waypoints=1000, knn=30, use_early_cell_as_start=True)
    pt = pr.pseudotime

    R = rna[cells]
    R = np.log1p(R.div(R.sum(axis=0).replace(0, 1), axis=1) * R.sum(axis=0).mean())
    res = {}
    rows = []
    for t in terms:
        lin = lab.index[lab.isin(prog + [t])]
        order = pt[lin].sort_values()
        x = order.to_numpy()
        pos = 0
        for g in genes:
            c = dorc.loc[g, order.index].to_numpy(float)
            r = R.loc[g, order.index].to_numpy(float)
            cs = minmax(lowess(c, x, frac=a.frac, return_sorted=False))
            rs = minmax(lowess(r, x, frac=a.frac, return_sorted=False))
            m = float(np.mean(cs - rs))
            rows.append({"lineage": t, "gene": g, "mean_residual": m})
            pos += m > 0
        res[t] = {"n_cells": int(len(lin)), "n_genes": len(genes), "pct_positive": 100.0 * pos / max(len(genes), 1)}
    df = pd.DataFrame(rows)
    pooled = df.groupby("gene")["mean_residual"].mean()
    summary = {"root_cell": root, "lineages": res,
               "positive_residual_fraction_pct": float(100 * (pooled > 0).mean()),
               "positive_residual_fraction_pct_gene_lineage_pairs": float(100 * (df.mean_residual > 0).mean()),
               "embedding": f"LSI ({a.n_comp} comps) in place of cisTopic topics",
               "lineage_definition": "progenitor + one terminal cell type (author labels)"}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "residuals.tsv", sep="\t", index=False)
    pt.rename("pseudotime").to_csv(out / "pseudotime.tsv", sep="\t")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
