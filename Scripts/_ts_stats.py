"""Numerical core for the Tabula Sapiens skill.

Clean-room implementations of the statistics the Tabula Sapiens 2.0 paper
(Tabula Sapiens Consortium, *Cell* 2026) delegates to external packages.
Each is written from the published description and, where the upstream
source is permissive, cross-checked against its documented behaviour --
never copied.

  ``tau``                 tspex's tissue-specificity statistic, the
                          backbone of Figures 2-3.
  ``decontx``             DecontX's variational-EM ambient-RNA model
                          (celda). Validated against the ``decontXcounts``
                          layer the published h5ads already carry.
  ``consensus_nmf``       cNMF's consensus non-negative matrix
                          factorisation (Kotliar 2019), used for the
                          senescence gene modules in Figure 4.
  ``aucell``              SCENIC's rank-based regulon activity score.
  ``grn_importance``      GRNBoost2-style co-expression importance, the
                          first stage of the SCENIC workflow.
  ``negative_binomial_de``  pseudobulk NB differential expression, an
                          Apache-2 stand-in for the GPL edgeR path.
  ``benjamini_hochberg``  shared FDR control.

Why this module exists apart from ``tabula_sapiens_skill``: it is pure
array math with no I/O, no argparse and no project paths, so it can be
unit-tested directly and reused by any other skill.

License: Apache-2.0. tspex (MIT), celda/DecontX (MIT), cNMF (MIT),
pySCENIC (GPL-3) and GSEApy (MIT) are referenced for their published
algorithms only; none is imported or vendored. pySCENIC in particular is
GPL-3 and is therefore reimplemented rather than depended on, matching
this codebase's no-GPL-runtime boundary.
"""

from __future__ import annotations

from typing import Any, Optional


def _np():
    import numpy as np
    return np


# ---------------------------------------------------------------------------
# Specificity
# ---------------------------------------------------------------------------

def tau(expression: Any, *, axis: int = 1) -> Any:
    """Tissue/cell-type specificity statistic (Yanai 2005).

    ``tau = sum(1 - x_i / max(x)) / (N - 1)`` over the ``N`` cell types of
    each gene's mean-expression vector. 0 means perfectly ubiquitous, 1
    means confined to a single cell type.

    The paper computes this on the mean of the **log-normalised**
    expression per cell type, using the droplet subset only, and splits
    specific from non-specific at 0.85.

    Parameters
    ----------
    expression
        2-D array of mean expression, genes along the axis that is *not*
        ``axis``.
    axis
        Axis holding the cell types (default 1, i.e. genes x cell types).

    Returns
    -------
    1-D array of tau, one per gene. Genes whose maximum is 0 -- expressed
    nowhere -- yield NaN rather than a spurious 0 or 1: they carry no
    specificity information at all, and the paper handles them separately
    (SHOX and ZBED1 are exactly this case).
    """
    np = _np()
    x = np.asarray(expression, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"expected a 2-D matrix, got shape {x.shape}")
    if axis == 0:
        x = x.T
    n = x.shape[1]
    if n < 2:
        raise ValueError("tau needs at least two cell types")

    mx = x.max(axis=1)
    out = np.full(x.shape[0], np.nan)
    ok = mx > 0
    if not ok.any():
        return out
    normed = x[ok] / mx[ok][:, None]
    out[ok] = (1.0 - normed).sum(axis=1) / (n - 1)
    return out


# ---------------------------------------------------------------------------
# Ambient RNA
# ---------------------------------------------------------------------------

def decontx(
    counts: Any,
    clusters: Any,
    *,
    max_iter: int = 200,
    tol: float = 1e-3,
    delta: "tuple[float, float]" = (10.0, 10.0),
    seed: int = 0,
) -> "dict[str, Any]":
    """Estimate and remove ambient ("soup") RNA contamination.

    A clean-room implementation of DecontX (Yang 2020, *Genome Biol*).
    Each cell's counts are modelled as a two-component mixture: its own
    cluster's expression profile, and a background profile built from
    every *other* cluster. Variational EM alternates between

      * **E step** -- per gene, split the observed counts between the
        native and contaminating components in proportion to their
        current expected rates;
      * **M step** -- re-estimate each cluster's native profile from the
        native share, and each cell's contamination fraction from a Beta
        posterior with pseudocounts ``delta``.

    The background is deliberately built from other clusters rather than
    from empty droplets: the published h5ads have no empty-droplet
    matrix, and the paper's own DecontX run used the same
    cluster-contrast formulation.

    Parameters
    ----------
    counts
        ``(n_cells, n_genes)`` raw count matrix, dense or scipy-sparse.
    clusters
        Per-cell cluster label. Any hashable dtype.
    max_iter, tol
        EM iteration cap and convergence threshold on the mean change in
        contamination fraction.
    delta
        Beta prior pseudocounts ``(native, contamination)``.
    seed
        Unused by the deterministic EM; accepted so callers can pass one
        uniformly.

    Returns
    -------
    ``{"contamination": (n_cells,), "decontaminated": (n_cells, n_genes),
    "n_iter": int, "converged": bool}``
    """
    np = _np()
    X = _as_dense(counts)
    n_cells, n_genes = X.shape
    labels = np.asarray(clusters)
    if labels.shape[0] != n_cells:
        raise ValueError(
            f"clusters has {labels.shape[0]} entries for {n_cells} cells")

    uniq, idx = np.unique(labels, return_inverse=True)
    k = uniq.size
    if k < 2:
        # With one cluster there is no contrast to separate native from
        # ambient, so decontamination is undefined. Return the input
        # untouched rather than inventing a correction.
        return {"contamination": np.zeros(n_cells),
                "decontaminated": X.copy(), "n_iter": 0, "converged": True}

    cell_totals = X.sum(axis=1)
    eps = 1e-10

    # Cluster profiles, rows summing to 1.
    phi = np.zeros((k, n_genes))
    for j in range(k):
        phi[j] = X[idx == j].sum(axis=0)
    phi /= np.maximum(phi.sum(axis=1, keepdims=True), eps)

    # Background for cluster j = every other cluster pooled. Using the
    # complement rather than the global total is what gives the model a
    # contrast: a gene high in j and nowhere else must be native to j.
    total = X.sum(axis=0)
    eta = np.zeros((k, n_genes))
    for j in range(k):
        eta[j] = np.maximum(total - X[idx == j].sum(axis=0), 0.0)
    eta /= np.maximum(eta.sum(axis=1, keepdims=True), eps)

    contamination = np.full(n_cells, 0.5)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        native_rate = (1.0 - contamination)[:, None] * phi[idx]
        contam_rate = contamination[:, None] * eta[idx]
        denom = native_rate + contam_rate + eps
        p_native = native_rate / denom

        native_counts = X * p_native
        contam_counts = X - native_counts

        new_contam = ((contam_counts.sum(axis=1) + delta[1]) /
                      (cell_totals + delta[0] + delta[1]))

        for j in range(k):
            m = idx == j
            if m.any():
                s = native_counts[m].sum(axis=0)
                phi[j] = s / max(s.sum(), eps)

        shift = float(np.abs(new_contam - contamination).mean())
        contamination = new_contam
        if shift < tol:
            converged = True
            break

    native_rate = (1.0 - contamination)[:, None] * phi[idx]
    contam_rate = contamination[:, None] * eta[idx]
    p_native = native_rate / (native_rate + contam_rate + eps)
    decon = X * p_native

    return {"contamination": contamination, "decontaminated": decon,
            "n_iter": it, "converged": converged}


def _as_dense(m: Any) -> Any:
    np = _np()
    if hasattr(m, "toarray"):
        return np.asarray(m.toarray(), dtype=float)
    return np.asarray(m, dtype=float)


# ---------------------------------------------------------------------------
# Consensus NMF
# ---------------------------------------------------------------------------

def consensus_nmf(
    X: Any,
    k: int,
    *,
    n_runs: int = 10,
    density_threshold: float = 0.5,
    top_genes: int = 100,
    seed: int = 0,
    max_iter: int = 400,
) -> "dict[str, Any]":
    """Consensus NMF gene-expression programs (Kotliar 2019).

    Runs NMF ``n_runs`` times from different seeds, pools all
    ``n_runs * k`` spectra, drops outlier spectra whose mean distance to
    their nearest neighbours exceeds ``density_threshold`` (in units of
    the median such distance), then clusters the survivors into ``k``
    consensus programs by k-means on the L2-normalised spectra.

    The filtering step is what makes it *consensus*: a spectrum that only
    one restart found is an artefact of that restart, and averaging it in
    would blur a real program.

    Parameters
    ----------
    X
        ``(n_cells, n_genes)`` non-negative matrix.
    k
        Number of programs. The paper selects K=34 for the senescence
        modules by trading reconstruction error against stability.
    n_runs
        NMF restarts (the paper uses 10).
    density_threshold
        Outlier cut, relative to the median nearest-neighbour distance.
    top_genes
        Genes reported per program.
    seed, max_iter
        RNG seed and per-run NMF iteration cap.

    Returns
    -------
    ``{"spectra": (k, n_genes), "usage": (n_cells, k),
    "top_gene_idx": (k, top_genes), "n_spectra_kept": int,
    "n_spectra_total": int}``
    """
    np = _np()
    from sklearn.decomposition import NMF
    from sklearn.cluster import KMeans

    A = _as_dense(X)
    A[A < 0] = 0.0
    n_cells, n_genes = A.shape
    if k < 1:
        raise ValueError("k must be >= 1")
    if k > min(A.shape):
        raise ValueError(f"k={k} exceeds the data rank {min(A.shape)}")

    spectra = []
    for r in range(n_runs):
        model = NMF(n_components=k, init="nndsvdar", random_state=seed + r,
                    max_iter=max_iter, tol=1e-4)
        model.fit(A)
        H = model.components_
        norms = np.linalg.norm(H, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        spectra.append(H / norms)
    S = np.vstack(spectra)

    # ── density filtering ────────────────────────────────────────────
    n_neighbors = max(1, min(n_runs - 1, S.shape[0] - 1))
    sq = (S ** 2).sum(axis=1)
    d2 = np.maximum(sq[:, None] + sq[None, :] - 2 * (S @ S.T), 0.0)
    np.fill_diagonal(d2, np.inf)
    knn = np.sort(d2, axis=1)[:, :n_neighbors]
    local = np.sqrt(knn).mean(axis=1)
    med = np.median(local)
    keep = local <= med * (1.0 + density_threshold) if med > 0 else np.ones(
        S.shape[0], dtype=bool)
    if keep.sum() < k:
        # Never filter below what the clustering needs.
        keep = np.ones(S.shape[0], dtype=bool)

    Sf = S[keep]
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(Sf)
    consensus = np.zeros((k, n_genes))
    for j in range(k):
        m = km.labels_ == j
        consensus[j] = Sf[m].mean(axis=0) if m.any() else Sf.mean(axis=0)
    norms = np.linalg.norm(consensus, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    consensus /= norms

    # Usage: non-negative least squares of each cell on the programs.
    usage = _nnls_rows(A, consensus)
    order = np.argsort(-consensus, axis=1)[:, :top_genes]

    return {"spectra": consensus, "usage": usage, "top_gene_idx": order,
            "n_spectra_kept": int(keep.sum()), "n_spectra_total": int(S.shape[0])}


def _nnls_rows(A: Any, H: Any, *, iters: int = 200) -> Any:
    """Multiplicative-update NNLS for ``A ~ W @ H`` with H fixed."""
    np = _np()
    n, k = A.shape[0], H.shape[0]
    W = np.full((n, k), 1.0 / k)
    HHt = H @ H.T
    AHt = A @ H.T
    for _ in range(iters):
        denom = W @ HHt
        denom[denom < 1e-12] = 1e-12
        W = W * (AHt / denom)
    return W


# ---------------------------------------------------------------------------
# Regulon activity (SCENIC)
# ---------------------------------------------------------------------------

def grn_importance(
    X: Any,
    tf_idx: Any,
    *,
    target_idx: Optional[Any] = None,
    n_estimators: int = 50,
    max_depth: int = 3,
    seed: int = 0,
    top_k: int = 50,
) -> "dict[int, list[tuple[int, float]]]":
    """GRNBoost2-style co-expression importance, TF -> targets.

    For each target gene, fit a gradient-boosted regression on the TF
    expression matrix and read off feature importances; a TF with high
    importance for a target is a candidate regulator. This is the first
    stage of the SCENIC workflow, before motif pruning.

    Parameters
    ----------
    X
        ``(n_cells, n_genes)`` log-normalised expression.
    tf_idx
        Column indices of transcription factors.
    target_idx
        Column indices to model as targets (default: all genes).
    n_estimators, max_depth, seed
        Booster settings. Defaults are far smaller than GRNBoost2's
        because this runs per cell type over up to 1,600 TFs.
    top_k
        Targets retained per TF, by importance.

    Returns
    -------
    ``{tf_column_index: [(target_column_index, importance), ...]}``,
    each list sorted by descending importance and truncated to ``top_k``.
    """
    np = _np()
    from sklearn.ensemble import GradientBoostingRegressor

    A = _as_dense(X)
    tfs = np.asarray(tf_idx, dtype=int)
    targets = (np.arange(A.shape[1]) if target_idx is None
               else np.asarray(target_idx, dtype=int))
    TFX = A[:, tfs]

    links: "dict[int, list[tuple[int, float]]]" = {int(t): [] for t in tfs}
    for t in targets:
        y = A[:, t]
        if y.std() == 0:
            continue
        # Never let a TF predict itself.
        mask = tfs != t
        if not mask.any():
            continue
        model = GradientBoostingRegressor(
            n_estimators=n_estimators, max_depth=max_depth,
            random_state=seed, subsample=0.9)
        model.fit(TFX[:, mask], y)
        imp = model.feature_importances_
        for tf, w in zip(tfs[mask], imp):
            if w > 0:
                links[int(tf)].append((int(t), float(w)))

    for tf in links:
        links[tf].sort(key=lambda p: -p[1])
        del links[tf][top_k:]
    return links


def aucell(
    X: Any,
    gene_sets: "dict[str, list[int]]",
    *,
    auc_max_rank_frac: float = 0.05,
) -> "tuple[Any, list[str]]":
    """Rank-based regulon activity score (SCENIC's AUCell).

    Each cell's genes are ranked by expression; a regulon's score is the
    area under the recovery curve of its member genes within the top
    ``auc_max_rank_frac`` of that ranking, normalised to [0, 1]. Being
    rank-based, it is invariant to per-cell normalisation -- which is why
    SCENIC uses it instead of mean expression.

    Parameters
    ----------
    X
        ``(n_cells, n_genes)`` expression.
    gene_sets
        Regulon name -> member column indices.
    auc_max_rank_frac
        Fraction of the ranking that counts as "recovered".

    Returns
    -------
    ``(scores, names)`` with ``scores`` shaped ``(n_cells, n_regulons)``.
    """
    np = _np()
    A = _as_dense(X)
    n_cells, n_genes = A.shape
    max_rank = max(1, int(round(n_genes * auc_max_rank_frac)))

    # Descending rank of every gene in every cell. argsort twice gives
    # the rank directly and avoids a per-cell Python loop.
    order = np.argsort(-A, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(n_cells)[:, None]
    ranks[rows, order] = np.arange(n_genes)[None, :]

    names = list(gene_sets)
    out = np.zeros((n_cells, len(names)))
    for j, name in enumerate(names):
        members = np.asarray(gene_sets[name], dtype=int)
        if members.size == 0:
            continue
        r = ranks[:, members]
        # Recovery curve area, computed in closed form: each member
        # inside the cutoff contributes (max_rank - rank) steps.
        inside = r < max_rank
        area = np.where(inside, max_rank - r, 0).sum(axis=1)
        out[:, j] = area / float(max_rank * members.size)
    return out, names


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------

def negative_binomial_de(
    counts_a: Any,
    counts_b: Any,
    *,
    size_factors_a: Optional[Any] = None,
    size_factors_b: Optional[Any] = None,
    min_total: int = 10,
) -> "dict[str, Any]":
    """Pseudobulk negative-binomial differential expression.

    An Apache-2 stand-in for the edgeR path the upstream sex-difference
    analysis uses (edgeR is GPL, so it cannot be a runtime dependency
    here). Per gene: normalise by library size, estimate a common
    dispersion by method of moments, and test the log fold-change with a
    Wald statistic on the NB variance.

    This will not reproduce edgeR's exact p-values -- edgeR shrinks
    dispersions empirically across genes and uses an exact/QL test -- so
    treat it as concordant in direction and ranking rather than
    numerically identical.

    Parameters
    ----------
    counts_a, counts_b
        ``(n_samples, n_genes)`` pseudobulk count matrices for the two
        groups. Samples are donors or donor-tissue pseudobulks, never
        individual cells: cells within a donor are not independent, and
        treating them as such is what produces implausibly tiny p-values.
    size_factors_a, size_factors_b
        Per-sample normalisation. Defaults to library size / mean.
    min_total
        Genes with fewer total counts across both groups are skipped and
        returned as NaN.

    Returns
    -------
    ``{"log2fc", "pvalue", "padj", "mean_a", "mean_b", "tested"}``
    """
    np = _np()
    from scipy import stats

    A = _as_dense(counts_a)
    B = _as_dense(counts_b)
    if A.shape[1] != B.shape[1]:
        raise ValueError(f"gene count mismatch: {A.shape[1]} vs {B.shape[1]}")
    n_genes = A.shape[1]

    def _sf(M, given):
        if given is not None:
            return np.asarray(given, dtype=float)
        lib = M.sum(axis=1)
        m = lib.mean()
        return lib / m if m > 0 else np.ones(M.shape[0])

    sa, sb = _sf(A, size_factors_a), _sf(B, size_factors_b)
    sa[sa <= 0] = 1.0
    sb[sb <= 0] = 1.0
    na = A / sa[:, None]
    nb = B / sb[:, None]

    mean_a, mean_b = na.mean(axis=0), nb.mean(axis=0)
    log2fc = np.full(n_genes, np.nan)
    pval = np.full(n_genes, np.nan)
    tested = np.zeros(n_genes, dtype=bool)

    totals = A.sum(axis=0) + B.sum(axis=0)
    for g in range(n_genes):
        if totals[g] < min_total:
            continue
        xa, xb = na[:, g], nb[:, g]
        ma, mb = xa.mean(), xb.mean()
        if ma <= 0 and mb <= 0:
            continue
        pooled = np.concatenate([xa, xb])
        mu = pooled.mean()
        var = pooled.var(ddof=1) if pooled.size > 1 else 0.0
        # Method-of-moments dispersion: var = mu + phi*mu^2.
        phi = max((var - mu) / (mu ** 2), 0.0) if mu > 0 else 0.0
        se2 = 0.0
        for x, m in ((xa, ma), (xb, mb)):
            v = m + phi * m ** 2
            se2 += v / max(len(x), 1) / max(m, 1e-8) ** 2
        se = np.sqrt(se2)
        lfc = np.log2((mb + 1e-8) / (ma + 1e-8))
        log2fc[g] = lfc
        if se > 0:
            z = np.log((mb + 1e-8) / (ma + 1e-8)) / se
            pval[g] = 2.0 * stats.norm.sf(abs(z))
            tested[g] = True

    padj = benjamini_hochberg(pval)
    return {"log2fc": log2fc, "pvalue": pval, "padj": padj,
            "mean_a": mean_a, "mean_b": mean_b, "tested": tested}


def benjamini_hochberg(pvalues: Any) -> Any:
    """BH-adjusted p-values; NaN in, NaN out, positions preserved."""
    np = _np()
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * n / (np.arange(n) + 1)
    # Enforce monotonicity from the largest p downward.
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out
