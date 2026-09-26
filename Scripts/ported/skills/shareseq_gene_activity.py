# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/satijalab/seurat @ b56d19493937
# (R/preprocessing.R) for ma2020_shareseq. Unreviewed; provenance in Scripts/ported/registry.json.
#!/usr/bin/env python3
"""Peak-matrix gene activity scores (Seurat v3 CreateGeneActivityMatrix).

Ma et al., Cell 2020 (STAR Methods "Computational pairing") used Seurat v3
to compute scATAC gene activity before CCA label transfer. Seurat v3's
CreateGeneActivityMatrix assigns each peak to the nearest gene whose body
extended by 2 kb upstream it overlaps (GenomicRanges distanceToNearest,
distance 0, first hit on ties) and sums peak counts per gene. Ported here
so the activity matrix is reproducible without the retired v3 function.
"""
from __future__ import annotations

import argparse
import gzip
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

def _read_mtx(path):
    """MatrixMarket (.mtx/.mtx.gz, as deposited on GEO) or a scipy .npz cache."""
    path = str(path)
    if path.endswith(".npz"):
        return sp.load_npz(path)
    return scipy.io.mmread(gzip.open(path) if path.endswith(".gz") else path)



def read_genes(gtf, seq_levels):
    rows = []
    with (gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            if chrom.replace("chr", "") not in seq_levels:
                continue
            attr = dict(a.strip().split(" ", 1) for a in f[8].strip().split(";") if a.strip())
            name = attr.get("gene_name", attr.get("gene_id", "")).strip('"')
            rows.append((chrom, int(f[3]), int(f[4]), f[6], name))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--atac-mtx", required=True)
    ap.add_argument("--atac-barcodes", required=True)
    ap.add_argument("--peaks", required=True, help="BED (0-based) in matrix row order")
    ap.add_argument("--gtf", required=True, help="gene annotation GTF (gene rows)")
    ap.add_argument("--upstream", type=int, default=2000)
    ap.add_argument("--seq-levels", default=",".join([str(i) for i in range(1, 20)] + ["X", "Y"]))
    ap.add_argument("--out", required=True, help="genes x cells MatrixMarket (.mtx.gz); writes .genes/.barcodes alongside")
    a = ap.parse_args(argv)

    m = _read_mtx(a.atac_mtx).tocsr()
    pk = pd.read_csv(a.peaks, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    g = read_genes(a.gtf, set(a.seq_levels.split(",")))
    # Extend(upstream=2000): strand-aware, gene body included
    g["lo"] = np.where(g.strand == "-", g.start, g.start - a.upstream)
    g["hi"] = np.where(g.strand == "-", g.end + a.upstream, g.end)
    # Seurat v3 builds GRanges straight from "chr:start-end" names, so the BED
    # start is taken as a 1-based start (0 is bumped to 1).
    ps, pe = np.maximum(pk.start.to_numpy(), 1), pk.end.to_numpy()
    gene_of_peak = np.full(len(pk), -1)
    for c, gi in g.groupby("chrom").groups.items():
        gg = g.loc[gi]
        lo, hi = gg.lo.to_numpy(), gg.hi.to_numpy()
        order = np.argsort(lo, kind="stable")
        lo_s, hi_s, ids = lo[order], hi[order], np.asarray(gi)[order]
        pidx = np.flatnonzero(pk.chrom.to_numpy() == c)
        for p in pidx:
            k = np.searchsorted(lo_s, pe[p], side="right")
            hit = np.flatnonzero(hi_s[:k] >= ps[p])
            if len(hit):
                gene_of_peak[p] = ids[hit].min()  # first subject hit
    keep = gene_of_peak >= 0
    names = g.gene.to_numpy()[gene_of_peak[keep]]
    uniq, inv = np.unique(names, return_inverse=True)
    agg = sp.csr_matrix((np.ones(keep.sum()), (inv, np.flatnonzero(keep))), shape=(len(uniq), len(pk)))
    act = (agg @ m).tocoo()
    scipy.io.mmwrite(a.out.replace(".gz", ""), act)
    import shutil, os
    with open(a.out.replace(".gz", ""), "rb") as fi, gzip.open(a.out if a.out.endswith(".gz") else a.out + ".gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(a.out.replace(".gz", ""))
    base = a.out.replace(".mtx.gz", "").replace(".gz", "")
    pd.Series(uniq).to_csv(base + ".genes.txt", index=False, header=False)
    with (gzip.open(a.atac_barcodes, "rt") if a.atac_barcodes.endswith(".gz") else open(a.atac_barcodes)) as fh:
        open(base + ".barcodes.txt", "w").write(fh.read())
    print(f"gene activity: {len(uniq)} genes x {m.shape[1]} cells; {keep.sum()} of {len(pk)} peaks assigned")
    return 0


if __name__ == "__main__":
    sys.exit(main())
