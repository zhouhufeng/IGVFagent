# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/buenrostrolab/FigR @ 094f5aa036aa
# (R/DORCs.R, R/utils.R) for ma2020_shareseq. Unreviewed; provenance in Scripts/ported/registry.json.
#!/usr/bin/env python3
"""SHARE-seq peak-gene associations and DORCs (Ma et al., Cell 2020).

Python port of buenrostrolab/FigR R/DORCs.R (runGenePeakcorr, PeakGeneCor,
chunkCore, getDORCScores) and R/utils.R (centerCounts), the Buenrostro lab's
released implementation of the Ma 2020 peak-gene / DORC framework
(STAR Methods, "Peak-gene cis-association and DORC identification").

Algorithm (faithful to FigR):
  1. ATAC counts are centred per cell: count / mean(count over all peaks).
  2. Peaks and genes with zero total signal are dropped.
  3. Peak summits (GRanges resize(width=1, fix="center")) overlapping a TSS
     window (GRanges flank(TSS, width=W, both=TRUE)) define the tested pairs.
  4. Spearman correlation (average ranks for ties, as R cor()) between each
     peak and its gene across cells; the same for the 100 chromVAR background
     peaks of every tested peak (GC- and accessibility-matched).
  5. One-sided z-test: p = 1 - Phi((rObs - mean(rBg)) / sd(rBg)), sd with n-1.
  6. Keep rObs > 0; for peaks mapping to several genes keep the max rObs
     (FigR, keepMultiMappingPeaks=FALSE). The paper text instead keeps the
     smallest p (``--multimap min_p``); ``--multimap all`` keeps every pair.
  7. DORC genes: > cutoff significant peaks (10 for skin, 5 for GM12878).

Background peaks: pass the chromVAR matrix (``--bg-matrix``, 1-based peak
indices as returned by chromVAR::getBackgroundPeaks) for an exact match to
the R reference. Without it, ``bg-peaks`` ports getBackgroundPeaks (same
whitening, 50x50 bin grid, Gaussian bin kernel w=0.1), which cannot reproduce
R's random stream, so p-values then agree only in distribution.
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
from scipy.stats import norm, rankdata

def _read_mtx(path):
    """MatrixMarket (.mtx/.mtx.gz, as deposited on GEO) or a scipy .npz cache."""
    path = str(path)
    if path.endswith(".npz"):
        return sp.load_npz(path)
    return scipy.io.mmread(gzip.open(path) if path.endswith(".gz") else path)


try:
    import numba
    from numba import njit, prange
except ImportError:  # pragma: no cover
    numba = None


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_atac(mtx, barcodes, peaks):
    """Peaks x cells CSR counts, cell barcodes, and peak table (0-based BED)."""
    m = _read_mtx(mtx)
    m = sp.csr_matrix(m, dtype=np.float64)
    with _open(barcodes) as fh:
        bcs = [l.strip() for l in fh if l.strip()]
    pk = pd.read_csv(peaks, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"])
    if m.shape != (len(pk), len(bcs)):
        sys.exit(f"ATAC shape {m.shape} != peaks {len(pk)} x barcodes {len(bcs)}")
    return m, bcs, pk


def load_rna(path, sep_from=",", sep_to="."):
    """Genes x cells counts from a dense TSV (first column = gene) or h5ad."""
    if str(path).endswith(".h5ad"):
        import anndata
        a = anndata.read_h5ad(path)
        X = sp.csr_matrix(a.X.T, dtype=np.float64)
        return X, list(a.var_names), [b.replace(sep_from, sep_to) for b in a.obs_names]
    genes, rows = [], []
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = [b.replace(sep_from, sep_to) for b in header[1:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes.append(f[0])
            v = np.asarray(f[1:], dtype=np.float64)
            rows.append(sp.csr_matrix(v))
    return sp.vstack(rows).tocsr(), genes, cells


def load_tss(path):
    """TSS table exported from FigR (chrom, start, end, strand, gene_name; 1-based)."""
    t = pd.read_csv(path, sep="\t")
    need = {"chrom", "start", "end", "strand", "gene_name"}
    if not need.issubset(t.columns):
        sys.exit(f"TSS table needs columns {sorted(need)}")
    return t


# --------------------------------------------------------------------------
# FigR building blocks
# --------------------------------------------------------------------------
def center_counts(m):
    """FigR centerCounts: divide each cell (column) by its mean over features."""
    means = np.asarray(m.mean(axis=0)).ravel()
    means[means == 0] = np.nan
    c = m.tocsc(copy=True)
    c = c @ sp.diags(1.0 / means)
    c.data[~np.isfinite(c.data)] = 0.0
    return c.tocsr()


def summits_1based(pk):
    """resize(GRanges(start+1, end), width=1, fix='center') -> 1-based summit."""
    start1 = pk["start"].to_numpy() + 1
    width = pk["end"].to_numpy() - pk["start"].to_numpy()
    return start1 + (width - 1) // 2


def tss_windows(tss, pad):
    """GenomicRanges::flank(x, width=pad, both=TRUE) on width-1 TSS ranges."""
    s = tss["start"].to_numpy()
    e = tss["end"].to_numpy()
    minus = (tss["strand"] == "-").to_numpy()
    lo = np.where(minus, e - pad + 1, s - pad)
    hi = np.where(minus, e + pad, s + pad - 1)
    return lo, hi


def find_pairs(tss, pad, pk, peak_keep):
    """findOverlaps(query=TSS windows, subject=peak summits), both 1-based closed."""
    lo, hi = tss_windows(tss, pad)
    summit = summits_1based(pk)
    genes_i, peaks_i = [], []
    kept = np.flatnonzero(peak_keep)
    by_chr = {}
    for idx in kept:
        by_chr.setdefault(pk["chrom"].iat[idx], []).append(idx)
    for c in by_chr:
        arr = np.asarray(by_chr[c])
        order = np.argsort(summit[arr], kind="stable")
        by_chr[c] = (arr[order], summit[arr][order])
    chroms = tss["chrom"].to_numpy()
    for g in range(len(tss)):
        hit = by_chr.get(chroms[g])
        if hit is None:
            continue
        idx, pos = hit
        a = np.searchsorted(pos, lo[g], side="left")
        b = np.searchsorted(pos, hi[g], side="right")
        if b > a:
            sel = np.sort(idx[a:b])  # findOverlaps orders subject hits by index
            genes_i.extend([g] * len(sel))
            peaks_i.extend(sel.tolist())
    return np.asarray(genes_i, dtype=np.int64), np.asarray(peaks_i, dtype=np.int64)


def sparse_rank_rows(m):
    """Per row: average ranks over all columns, returned as (rank - zero_rank)
    on the non-zeros (CSR, same pattern) plus the row sum of squared
    deviations from the mean rank (n+1)/2."""
    n = m.shape[1]
    m = m.tocsr()
    data = np.empty_like(m.data)
    ss = np.empty(m.shape[0])
    mu = (n + 1) / 2.0
    for i in range(m.shape[0]):
        a, b = m.indptr[i], m.indptr[i + 1]
        v = m.data[a:b]
        nnz = b - a
        nz0 = n - nnz
        z0 = (nz0 + 1) / 2.0
        if nnz:
            r = rankdata(v, method="average") + nz0
            data[a:b] = r - z0
            ss[i] = nz0 * (z0 - mu) ** 2 + np.sum((r - mu) ** 2)
        else:
            ss[i] = 0.0
    out = sp.csr_matrix((data, m.indices.copy(), m.indptr.copy()), shape=m.shape)
    return out, ss


def dense_rank_rows(m):
    """Genes: centred average ranks (dense) and sum of squares."""
    n = m.shape[1]
    r = np.empty(m.shape, dtype=np.float32)  # half-integer ranks are exact in float32
    ss = np.empty(m.shape[0])
    for i in range(m.shape[0]):
        row = m[i].toarray().ravel() if sp.issparse(m) else np.asarray(m[i]).ravel()
        v = rankdata(row, method="average") - (n + 1) / 2.0
        r[i] = v
        ss[i] = float(np.sum(v * v))
    return r, ss


if numba is not None:
    @njit(parallel=True, cache=False)
    def _pair_cor(indptr, indices, data, ssx, G, ssy, gi, pi):
        out = np.empty(len(gi))
        for k in prange(len(gi)):
            p = pi[k]
            g = gi[k]
            acc = 0.0
            for j in range(indptr[p], indptr[p + 1]):
                acc += data[j] * G[g, indices[j]]
            den = np.sqrt(ssx[p] * ssy[g])
            out[k] = acc / den if den > 0 else np.nan
        return out
else:
    def _pair_cor(indptr, indices, data, ssx, G, ssy, gi, pi):
        out = np.empty(len(gi))
        for k in range(len(gi)):
            p, g = pi[k], gi[k]
            a, b = indptr[p], indptr[p + 1]
            acc = float(np.dot(data[a:b], G[g, indices[a:b]]))
            den = np.sqrt(ssx[p] * ssy[g])
            out[k] = acc / den if den > 0 else np.nan
        return out


# --------------------------------------------------------------------------
# chromVAR getBackgroundPeaks port
# --------------------------------------------------------------------------
def gc_content(pk, fasta):
    """Per-peak GC fraction from an indexed FASTA (chromVAR addGCBias analogue)."""
    try:
        import pysam
    except ImportError:
        sys.exit("pysam is required to compute GC from --genome-fasta")
    fa = pysam.FastaFile(fasta)
    gc = np.empty(len(pk))
    for i, (c, s, e) in enumerate(pk[["chrom", "start", "end"]].itertuples(index=False)):
        seq = fa.fetch(c, int(s), int(e)).upper()
        acgt = sum(seq.count(x) for x in "ACGT")
        gc[i] = (seq.count("G") + seq.count("C")) / acgt if acgt else np.nan
    return gc


def background_peaks(counts_per_peak, bias, niterations=100, w=0.1, bs=50, seed=123):
    """chromVAR::getBackgroundPeaks: 0-based background peak index matrix."""
    rng = np.random.default_rng(seed)
    intensity = np.log10(counts_per_peak)
    X = np.column_stack([intensity, bias])
    L = np.linalg.cholesky(np.cov(X, rowvar=False))  # lower = t(chol())
    T = np.linalg.solve(L, X.T).T
    b1 = np.linspace(T[:, 0].min(), T[:, 0].max(), bs)
    b2 = np.linspace(T[:, 1].min(), T[:, 1].max(), bs)
    bins = np.column_stack([np.repeat(b1, bs), np.tile(b2, bs)])
    d = np.sqrt(((bins[:, None, :] - bins[None, :, :]) ** 2).sum(-1))
    bin_p = norm.pdf(d, 0, w)
    from scipy.spatial import cKDTree
    member = cKDTree(bins).query(T, k=1)[1]
    dens = np.bincount(member, minlength=bs * bs).astype(float)
    peaks_in = [np.flatnonzero(member == j) for j in range(bs * bs)]
    out = np.empty((len(T), niterations), dtype=np.int64)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(bs * bs):
            ix = peaks_in[i]
            if not len(ix):
                continue
            p = np.where(dens > 0, bin_p[i] / dens, 0.0)
            p = p * dens  # per-bin mass = sum over its peaks of p/dens
            p /= p.sum()
            chosen = rng.choice(bs * bs, size=len(ix) * niterations, p=p)
            picks = np.array([peaks_in[j][rng.integers(len(peaks_in[j]))] for j in chosen])
            out[ix] = picks.reshape(len(ix), niterations)
    return out


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def _paired_inputs(args):
    atac, abcs, pk = load_atac(args.atac_mtx, args.atac_barcodes, args.peaks)
    rna, genes, rbcs = load_rna(args.rna)
    common = [b for b in abcs if b in set(rbcs)]
    if args.cells:
        with _open(args.cells) as fh:
            keep = {l.strip().split("\t")[0] for l in fh if l.strip()}
        common = [b for b in common if b in keep]
    if not common:
        sys.exit("no shared ATAC/RNA barcodes")
    ai = {b: i for i, b in enumerate(abcs)}
    ri = {b: i for i, b in enumerate(rbcs)}
    atac = atac[:, [ai[b] for b in common]].tocsr()
    rna = rna[:, [ri[b] for b in common]].tocsr()
    return atac, pk, rna, genes, common


def cmd_peakgene(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    atac, pk, rna, genes, cells = _paired_inputs(args)
    n_cells = len(cells)
    print(f"cells {n_cells:,}  peaks {atac.shape[0]:,}  genes {rna.shape[0]:,}")

    # RNA library-size normalisation (Seurat LogNormalize; ranks are invariant
    # to the scale factor and the log, so counts/total is sufficient).
    tot = np.asarray(rna.sum(axis=0)).ravel()
    tot[tot == 0] = 1.0
    rna = (rna @ sp.diags(1.0 / tot)).tocsr()

    raw_peak_sum = np.asarray(atac.sum(axis=1)).ravel()
    atac_c = center_counts(atac)
    peak_keep = raw_peak_sum != 0
    gene_sum = np.asarray(rna.sum(axis=1)).ravel()
    gene_keep = {g for g, s in zip(genes, gene_sum) if s != 0}

    tss = load_tss(args.tss)
    tss = tss[tss["gene_name"].isin(gene_keep)]
    tss = tss.drop_duplicates("gene_name", keep="first")  # TSSg[genesToKeep]
    if args.genes:
        with _open(args.genes) as fh:
            want = {l.strip() for l in fh if l.strip()}
        tss = tss[tss["gene_name"].isin(want)]  # FigR geneList
    gpos = {g: i for i, g in enumerate(genes)}
    tss = tss.reset_index(drop=True)
    gi, pi = find_pairs(tss, args.window, pk, peak_keep)
    print(f"pairs {len(gi):,}  genes with a pair {len(np.unique(gi)):,}")

    rna_rows = np.asarray([gpos[g] for g in tss["gene_name"]])
    G, ssy = dense_rank_rows(rna[rna_rows])
    X, ssx = sparse_rank_rows(atac_c)
    ind = (X.indptr.astype(np.int64), X.indices.astype(np.int64), X.data)

    r_obs = _pair_cor(*ind, ssx, G, ssy, gi, pi)

    if args.bg_matrix:
        bg = pd.read_csv(args.bg_matrix, sep="\t", header=None).to_numpy(np.int64) - 1
    else:
        gc = (pd.read_csv(args.gc, sep="\t", header=None).iloc[:, -1].to_numpy()
              if args.gc else gc_content(pk, args.genome_fasta))
        # chromVAR runs on the zero-filtered SE, so index into kept peaks
        kept = np.flatnonzero(peak_keep)
        bg_k = background_peaks(raw_peak_sum[kept], gc[kept], args.n_bg, seed=args.seed)
        bg = np.full((len(pk), args.n_bg), -1, dtype=np.int64)
        bg[kept] = kept[bg_k]
    if bg.shape[0] != len(pk):
        sys.exit(f"background matrix rows {bg.shape[0]} != peaks {len(pk)}")
    n_bg = bg.shape[1]
    r_bg = np.empty((len(gi), n_bg))
    for k in range(n_bg):
        r_bg[:, k] = _pair_cor(*ind, ssx, G, ssy, gi, bg[pi, k])

    tab = pd.DataFrame({"Gene": tss["gene_name"].to_numpy()[gi], "Peak": pi + 1, "rObs": r_obs})
    tab["rBgMean"] = r_bg.mean(axis=1)
    tab["rBgSD"] = r_bg.std(axis=1, ddof=1)
    tab["pvalZ"] = 1 - norm.cdf(tab["rObs"], loc=tab["rBgMean"], scale=tab["rBgSD"])
    tab["PeakRanges"] = (pk["chrom"].to_numpy()[pi] + ":" + (pk["start"].to_numpy()[pi] + 1).astype(str)
                         + "-" + pk["end"].to_numpy()[pi].astype(str))
    tab.to_csv(out / "all_pairs.tsv.gz", sep="\t", index=False, float_format="%.10g")

    return _summarize(tab, args, out, n_cells)


def _summarize(tab, args, out, n_cells):
    """Filter (FigR: rObs > 0, per-peak dedup), write the association table,
    DORC ranking and summary. Pre-dedup statistics are reported too: the
    paper's 63,110 associations and 83.9% single-gene figure count every
    significant peak-gene pair before collapsing multi-gene peaks."""
    res = tab
    if args.pos_only:
        res = res[res["rObs"] > 0]
    if args.multimap == "max_cor":
        res = res[res["rObs"] == res.groupby("Peak")["rObs"].transform("max")]
    elif args.multimap == "min_p":
        res = res[res["pvalZ"] == res.groupby("Peak")["pvalZ"].transform("min")]
    res = res[["Peak", "PeakRanges", "Gene", "rObs", "pvalZ"]]
    res.to_csv(out / "gene_peak_cor.tsv.gz", sep="\t", index=False, float_format="%.10g")

    sig = res[res["pvalZ"] <= args.p_cut]
    per_gene = sig.groupby("Gene").size().sort_values(ascending=False)
    dorcs = per_gene[per_gene > args.dorc_cutoff]
    pre = tab[(tab["pvalZ"] <= args.p_cut) & ((tab["rObs"] > 0) if args.pos_only else True)]
    pre_gene = pre.groupby("Gene").size().sort_values(ascending=False)
    genes_per_peak = pre.groupby("Peak")["Gene"].nunique()
    assoc_single = float((pre["Peak"].map(genes_per_peak) == 1).mean() * 100) if len(pre) else float("nan")
    per_gene.rename("n_sig_peaks").to_csv(out / "dorc_rank.tsv", sep="\t", header=True)
    pre_gene.rename("n_sig_peaks").to_csv(out / "dorc_rank_before_dedup.tsv", sep="\t", header=True)
    nan = float("nan")
    summary = {
        "n_cells": n_cells,
        "n_pairs_tested": int(len(tab)),
        "window_bp": args.window,
        "p_cut": args.p_cut,
        "multimap": args.multimap,
        "n_sig_associations": int(len(sig)),
        "n_genes_with_sig": int(len(per_gene)),
        "sig_per_gene_mean": float(per_gene.mean()) if len(per_gene) else nan,
        "dorc_cutoff": args.dorc_cutoff,
        "n_dorcs": int(len(dorcs)),
        "n_sig_associations_before_dedup": int(len(pre)),
        "n_genes_with_sig_before_dedup": int(len(pre_gene)),
        "before_dedup_sig_per_gene_mean": float(pre_gene.mean()) if len(pre_gene) else nan,
        "n_dorcs_before_dedup": int((pre_gene > args.dorc_cutoff).sum()),
        "single_gene_assoc_fraction_pct": assoc_single,
        "peaks_4plus_genes_fraction_pct": float((genes_per_peak >= 4).mean() * 100) if len(genes_per_peak) else nan,
        "top_dorcs": dorcs.head(30).to_dict(),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "top_dorcs"}, indent=2))
    return 0


def cmd_summarize(args):
    """Re-derive the filtered table and summary from an existing all_pairs table."""
    out = Path(args.out)
    tab = pd.read_csv(out / "all_pairs.tsv.gz", sep="\t")
    prev = json.loads((out / "summary.json").read_text()) if (out / "summary.json").is_file() else {}
    args.window = prev.get("window_bp", args.window)
    return _summarize(tab, args, out, prev.get("n_cells"))


def cmd_dorc_scores(args):
    """FigR getDORCScores: per cell, sum of centred counts over each gene's
    significant peaks."""
    atac, abcs, pk = load_atac(args.atac_mtx, args.atac_barcodes, args.peaks)
    tab = pd.read_csv(args.gene_peak, sep="\t")
    tab = tab[tab["pvalZ"] <= args.p_cut]
    if args.genes:
        with _open(args.genes) as fh:
            want = {l.strip() for l in fh if l.strip()}
        tab = tab[tab["Gene"].isin(want)]  # FigR geneList
    ac = center_counts(atac)
    genes = sorted(tab["Gene"].unique())
    rows = []
    for g in genes:
        pks = np.unique(tab.loc[tab["Gene"] == g, "Peak"].to_numpy()) - 1
        rows.append(sp.csr_matrix(ac[pks].sum(axis=0)))
    M = sp.vstack(rows).tocsr()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(M.toarray(), index=genes, columns=abcs)
    df.to_csv(out, sep="\t", float_format="%.6g")
    print(f"DORC scores: {M.shape[0]} genes x {M.shape[1]} cells -> {out}")
    return 0


def cmd_bg_peaks(args):
    atac, abcs, pk = load_atac(args.atac_mtx, args.atac_barcodes, args.peaks)
    s = np.asarray(atac.sum(axis=1)).ravel()
    keep = np.flatnonzero(s > 0)
    gc = (pd.read_csv(args.gc, sep="\t", header=None).iloc[:, -1].to_numpy()
          if args.gc else gc_content(pk, args.genome_fasta))
    bg = background_peaks(s[keep], gc[keep], args.n_bg, seed=args.seed)
    full = np.zeros((len(pk), args.n_bg), dtype=np.int64)
    full[keep] = keep[bg] + 1
    pd.DataFrame(full).to_csv(args.out, sep="\t", header=False, index=False)
    print(f"background peaks: {full.shape} -> {args.out}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="shareseq-dorc", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--atac-mtx", required=True, help="peaks x cells MatrixMarket (.gz ok)")
        p.add_argument("--atac-barcodes", required=True)
        p.add_argument("--peaks", required=True, help="BED of peaks, same order as the matrix rows")

    p = sub.add_parser("peakgene", help="peak-gene Spearman correlations with background z-test + DORC calls")
    common(p)
    p.add_argument("--rna", required=True, help="genes x cells dense TSV (.gz) or cells x genes h5ad")
    p.add_argument("--tss", required=True, help="TSS table: chrom start end strand gene_name (1-based)")
    p.add_argument("--cells", help="optional barcode list restricting the cells")
    p.add_argument("--genes", help="optional gene list (FigR geneList) restricting the tested genes")
    p.add_argument("--window", type=int, default=50000, help="bp padded either side of each TSS")
    p.add_argument("--bg-matrix", help="chromVAR background peaks (peaks x n_bg, 1-based) for exact R parity")
    p.add_argument("--gc", help="per-peak GC table (last column), in peak order")
    p.add_argument("--genome-fasta", help="indexed FASTA to compute GC if --gc is absent")
    p.add_argument("--n-bg", type=int, default=100)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--pos-only", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--multimap", choices=["max_cor", "min_p", "all"], default="max_cor")
    p.add_argument("--p-cut", type=float, default=0.05)
    p.add_argument("--dorc-cutoff", type=int, default=10, help="DORC = more than this many significant peaks")
    p.add_argument("--out", required=True, help="output directory")
    p.set_defaults(func=cmd_peakgene)

    p = sub.add_parser("summarize", help="re-filter an existing peakgene output (all_pairs.tsv.gz)")
    p.add_argument("--out", required=True, help="peakgene output directory")
    p.add_argument("--window", type=int, default=50000)
    p.add_argument("--pos-only", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--multimap", choices=["max_cor", "min_p", "all"], default="max_cor")
    p.add_argument("--p-cut", type=float, default=0.05)
    p.add_argument("--dorc-cutoff", type=int, default=10)
    p.set_defaults(func=cmd_summarize)

    p = sub.add_parser("dorc-scores", help="per-cell DORC scores from significant pairs")
    common(p)
    p.add_argument("--gene-peak", required=True, help="gene_peak_cor.tsv(.gz) from peakgene")
    p.add_argument("--genes", help="optional gene list (e.g. the DORC genes)")
    p.add_argument("--p-cut", type=float, default=0.05)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_dorc_scores)

    p = sub.add_parser("bg-peaks", help="chromVAR-style GC/accessibility-matched background peaks")
    common(p)
    p.add_argument("--gc")
    p.add_argument("--genome-fasta")
    p.add_argument("--n-bg", type=int, default=100)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_bg_peaks)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
