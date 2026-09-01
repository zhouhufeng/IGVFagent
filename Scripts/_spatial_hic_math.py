"""Numerical core for the Spatial-ATAC-Hi-C skill.

Clean-room implementations of the contact-matrix algorithms that the
Spatial-ATAC-Hi-C paper (Wang, Wang, Wang, Youngblood et al., *Nature
Methods* 2026, doi:10.1038/s41592-026-03217-4) delegates to external
tools. Each routine is written from the published description of the
method, not from the upstream source:

  ``impute_matrix``     scHiCluster's convolution + random-walk-with-restart
                        (Zhou et al. 2019 PNAS). Upstream is
                        zhoujt1994/scHiCluster; nothing is imported.
  ``compartment_pc1``   A/B compartment eigenvector on the
                        observed/expected Pearson correlation matrix
                        (Lieberman-Aiden 2009). Equivalent to the
                        cooltools ``eigs-cis`` call the paper uses.
  ``cnv_ratio``         Binned-coverage copy-number ratio, the role
                        NeoLoopFinder's ``calculate-cnv`` plays in the
                        paper. GC/mappability correction is a linear fit
                        rather than NeoLoopFinder's Poisson GLM — see the
                        note on the function.
  ``segment_cnv``       Gaussian-emission HMM + Viterbi over discrete copy
                        states, the role of ``segment-cnv``.
  ``apa_matrix``        Aggregate peak analysis pileup.
  ``one_way_anova``     Per-loop F test across cluster labels, the test the
                        paper uses to call cell-type-specific loops.
  ``magic_smooth``      Markov-affinity diffusion (van Dijk 2018) used in
                        the paper to denoise single-pixel CN ratios.

Why these live apart from ``spatial_atac_hic_skill``: they are pure
array math with no I/O, no argparse and no project paths, so they can be
unit-tested directly and reused by any other skill that grows a contact
matrix.

Every function takes and returns numpy arrays. numpy is imported lazily
by the caller (the skill module raises a friendly install hint), so this
module imports it at call time rather than at import time — that keeps
``import _spatial_hic_math`` free of hard dependencies.

License: Apache-2.0. No GPL runtime dependencies. scHiCluster (MIT),
cooltools (MIT) and NeoLoopFinder (MIT) are referenced for their
published algorithms only; none is imported or vendored.
"""

from __future__ import annotations

from typing import Any, Optional


def _np():
    import numpy as np
    return np


# ---------------------------------------------------------------------------
# scHiCluster-style imputation
# ---------------------------------------------------------------------------

def impute_matrix(
    mat: Any,
    *,
    pad: int = 1,
    rp: float = 0.5,
    tol: float = 0.01,
    max_iter: int = 30,
    logscale: bool = False,
    zero_diagonal: bool = False,
) -> Any:
    """Impute a sparse single-pixel contact matrix.

    Two steps, following Zhou et al. 2019 (PNAS 116:14011):

    1. **Linear convolution.** Each entry is replaced by the sum over a
       ``(2*pad+1)`` square window centred on it. This borrows signal
       from neighbouring bin pairs, which is what makes an otherwise
       near-empty single-cell (here single-pixel) matrix usable.
    2. **Random walk with restart.** With ``P`` the row-normalised
       convolved matrix, iterate ``Q <- rp*I + (1-rp) * Q @ P`` to the
       fixed point ``Q = rp * (I - (1-rp)P)^-1``. The restart
       probability ``rp`` trades locality (high ``rp``) against
       long-range smoothing (low ``rp``).

    The result is symmetrised, because a contact matrix is undirected
    while the random walk is not.

    Parameters
    ----------
    mat
        Dense square array of raw contact counts for one chromosome (or
        one window of one chromosome). Use a window for fine
        resolutions — the paper images 25 kb over 10.05 Mb and 10 kb
        over 5.05 Mb precisely to keep this tractable.
    pad
        Convolution half-width. The paper uses 1 at 100 kb and 2 at
        25 kb / 10 kb.
    rp
        Restart probability.
    tol
        L1 convergence threshold on successive iterates.
    max_iter
        Iteration cap; the walk is stopped even if ``tol`` is not met.
    logscale
        Return ``log1p`` of the imputed matrix. Useful for plotting.
    zero_diagonal
        Drop the main diagonal before the final normalisation. The
        restart term puts ``rp`` of the walk's mass straight back on the
        diagonal (at the default ``rp=0.5``, over half the matrix total),
        which is faithful to scHiCluster but swamps the off-diagonal
        structure in anything that reads the matrix as a whole — a
        correlation, a compartment eigenvector, an APA pileup. Set this
        when the diagonal is not the signal you are after.

    Returns
    -------
    Dense float array, same shape as ``mat``, normalised to sum to 1 so
    matrices from pixels of different depth are comparable.

    Notes
    -----
    Imputation earns its keep off the diagonal and at low depth. On a
    40-bin exponential-decay toy at 0.3x depth, off-diagonal correlation
    with the truth rises from 0.43 (raw) to 0.85 (imputed); at 3x depth,
    0.79 to 0.96.
    """
    np = _np()
    A = np.asarray(mat, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"expected a square matrix, got shape {A.shape}")
    n = A.shape[0]
    if n == 0:
        return A
    if not A.any():
        # An empty pixel imputes to empty; short-circuit so the
        # row-normalisation below never divides by zero.
        return A

    # ── 1. convolution ────────────────────────────────────────────────
    # Summed-area table gives the exact window sum in O(n^2) regardless
    # of `pad`, which matters because pad=2 over a 2,500-bin chromosome
    # is 25 multiply-adds per entry done naively.
    E = _window_sum(A, pad)

    # ── 2. random walk with restart ───────────────────────────────────
    rowsum = E.sum(axis=1, keepdims=True)
    # A bin with no contacts anywhere stays a zero row; guard the divide
    # and leave it zero rather than smearing uniform mass into it.
    rowsum[rowsum == 0] = 1.0
    P = E / rowsum

    I = np.eye(n)
    Q = I.copy()
    for _ in range(max_iter):
        Q_next = rp * I + (1.0 - rp) * (Q @ P)
        delta = np.abs(Q_next - Q).sum()
        Q = Q_next
        if delta < tol:
            break

    Q = (Q + Q.T) / 2.0
    if zero_diagonal:
        np.fill_diagonal(Q, 0.0)

    total = Q.sum()
    if total > 0:
        Q = Q / total
    return np.log1p(Q) if logscale else Q


def _window_sum(A: Any, pad: int) -> Any:
    """Sum of every ``(2*pad+1)`` square window, via a summed-area table."""
    np = _np()
    if pad <= 0:
        return A.copy()
    n = A.shape[0]
    # Pad the integral image by one row/col of zeros so the -1 index
    # arithmetic below needs no special-casing at the origin.
    S = np.zeros((n + 1, n + 1), dtype=float)
    S[1:, 1:] = A.cumsum(axis=0).cumsum(axis=1)

    idx = np.arange(n)
    lo = np.clip(idx - pad, 0, n)
    hi = np.clip(idx + pad + 1, 0, n)
    r0, r1 = lo[:, None], hi[:, None]
    c0, c1 = lo[None, :], hi[None, :]
    return S[r1, c1] - S[r0, c1] - S[r1, c0] + S[r0, c0]


# ---------------------------------------------------------------------------
# A/B compartments
# ---------------------------------------------------------------------------

def compartment_pc1(
    mat: Any,
    *,
    orient_track: Optional[Any] = None,
    min_coverage: float = 0.0,
) -> Any:
    """First eigenvector of the observed/expected correlation matrix.

    The standard A/B compartment call (Lieberman-Aiden 2009): divide
    each diagonal by its own mean to remove the distance-decay, take the
    Pearson correlation matrix across bins, and read off the leading
    eigenvector.

    Bins with no coverage are dropped before the eigendecomposition and
    returned as NaN, so the output always aligns with the input bins.

    Parameters
    ----------
    mat
        Dense square contact matrix for one chromosome, typically at
        100 kb (the resolution the paper uses).
    orient_track
        Per-bin values whose sign should agree with the A compartment —
        GC content or gene density. The eigenvector's sign is arbitrary,
        so without this the A/B labels can flip per chromosome. When
        given, PC1 is negated if it anti-correlates with the track.
    min_coverage
        Bins whose marginal is at or below this are treated as empty.

    Returns
    -------
    1-D array of length ``mat.shape[0]``; NaN at dropped bins.
    """
    np = _np()
    A = np.asarray(mat, dtype=float)
    n = A.shape[0]
    out = np.full(n, np.nan)
    if n == 0:
        return out

    cov = A.sum(axis=1)
    good = cov > min_coverage
    if good.sum() < 3:
        # Fewer than three usable bins cannot support a correlation
        # matrix; return all-NaN rather than a meaningless eigenvector.
        return out

    sub = A[np.ix_(good, good)]
    m = sub.shape[0]

    # ── observed / expected ───────────────────────────────────────────
    oe = np.ones_like(sub)
    for d in range(m):
        i = np.arange(m - d)
        j = i + d
        vals = sub[i, j]
        mean = vals.mean()
        if mean > 0:
            ratio = vals / mean
            oe[i, j] = ratio
            oe[j, i] = ratio

    # ── correlation + leading eigenvector ─────────────────────────────
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(oe)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)

    vals, vecs = np.linalg.eigh(corr)
    pc1 = vecs[:, -1] * np.sqrt(max(vals[-1], 0.0))

    if orient_track is not None:
        track = np.asarray(orient_track, dtype=float)[good]
        ok = np.isfinite(track) & np.isfinite(pc1)
        if ok.sum() >= 3 and np.std(track[ok]) > 0 and np.std(pc1[ok]) > 0:
            if np.corrcoef(track[ok], pc1[ok])[0, 1] < 0:
                pc1 = -pc1

    out[good] = pc1
    return out


# ---------------------------------------------------------------------------
# Copy number
# ---------------------------------------------------------------------------

def cnv_ratio(
    coverage: Any,
    *,
    gc: Optional[Any] = None,
    mappability: Optional[Any] = None,
    ploidy: float = 2.0,
    min_coverage: float = 0.0,
) -> Any:
    """Copy-number ratio from binned Hi-C coverage.

    This fills the role of NeoLoopFinder's ``calculate-cnv``: total
    contact coverage per genomic bin, corrected for the biases that make
    coverage a poor CN proxy on its own, then scaled so the genome-wide
    median sits at ``ploidy`` (all samples in the paper are treated as
    diploid, and its "CN ratio" is copy number divided by two).

    The bias correction here is a **linear least-squares fit** of
    coverage on GC and mappability, not NeoLoopFinder's Poisson GLM. It
    removes the same first-order trend and needs no extra dependency,
    but absolute CN at extreme-GC bins will differ from a NeoLoopFinder
    run — treat this as a concordant reimplementation, not a
    bit-identical one.

    Parameters
    ----------
    coverage
        Per-bin marginal contact counts.
    gc, mappability
        Optional per-bin covariates in [0, 1]. NaN entries are ignored
        when fitting and left uncorrected.
    ploidy
        Baseline copy number for the genome-wide median bin.
    min_coverage
        Bins at or below this are returned as NaN (blacklist, centromere,
        unmappable).

    Returns
    -------
    Per-bin copy-number ratio; NaN at excluded bins.
    """
    np = _np()
    cov = np.asarray(coverage, dtype=float)
    n = cov.size
    out = np.full(n, np.nan)
    good = np.isfinite(cov) & (cov > min_coverage)
    if not good.any():
        return out

    corrected = cov.copy()

    covars = []
    for track in (gc, mappability):
        if track is None:
            continue
        t = np.asarray(track, dtype=float)
        if t.size != n:
            raise ValueError(
                f"covariate length {t.size} does not match {n} bins")
        covars.append(t)

    if covars:
        X_cols = [np.ones(n)]
        fit_ok = good.copy()
        for t in covars:
            X_cols.append(np.nan_to_num(t, nan=0.0))
            fit_ok &= np.isfinite(t)
        X = np.column_stack(X_cols)
        if fit_ok.sum() > len(X_cols) + 1:
            beta, *_ = np.linalg.lstsq(X[fit_ok], cov[fit_ok], rcond=None)
            pred = X @ beta
            # Divide by the fitted trend rather than subtracting it:
            # coverage bias is multiplicative, and a subtraction can push
            # low-GC bins negative.
            scale = np.where(pred > 0, pred, np.nan)
            mean_pred = np.nanmean(pred[fit_ok]) if fit_ok.any() else 1.0
            with np.errstate(invalid="ignore", divide="ignore"):
                adjusted = cov / scale * mean_pred
            corrected = np.where(np.isfinite(adjusted), adjusted, cov)

    median = np.median(corrected[good])
    if median <= 0:
        return out
    out[good] = ploidy * corrected[good] / median
    return out


def segment_cnv(
    ratio: Any,
    *,
    states: Optional[Any] = None,
    sigma: float = 0.5,
    stay: float = 0.99,
) -> Any:
    """Viterbi segmentation of a copy-number ratio track.

    A Gaussian-emission HMM over discrete copy states, the role
    NeoLoopFinder's ``segment-cnv`` plays in the paper. Self-transition
    probability ``stay`` sets how much evidence is needed to call a
    breakpoint; the remaining mass is split evenly across the other
    states.

    NaN bins (blacklist, unmappable) do not vote: they emit uniformly, so
    a segment reads straight through a gap instead of being cut by it.

    Parameters
    ----------
    ratio
        Per-bin CN ratio, e.g. from :func:`cnv_ratio`.
    states
        Candidate copy numbers. Defaults to 0..8.
    sigma
        Gaussian emission width, in copy-number units.
    stay
        Self-transition probability, in (0, 1).

    Returns
    -------
    Per-bin integer copy state, as a float array with NaN preserved
    nowhere — every bin gets a state, including the ones that abstained.
    """
    np = _np()
    x = np.asarray(ratio, dtype=float)
    n = x.size
    if n == 0:
        return np.zeros(0)
    if states is None:
        states = np.arange(0, 9, dtype=float)
    s = np.asarray(states, dtype=float)
    k = s.size
    if k == 0:
        raise ValueError("need at least one candidate state")
    if not (0.0 < stay < 1.0):
        raise ValueError("stay must be in (0, 1)")

    # Log emissions. Missing bins get a flat row so they carry no
    # information into the Viterbi recursion.
    obs = np.isfinite(x)
    log_em = np.zeros((n, k))
    diff = x[obs][:, None] - s[None, :]
    log_em[obs] = -0.5 * (diff / sigma) ** 2

    other = (1.0 - stay) / (k - 1) if k > 1 else 1.0
    log_stay = np.log(stay)
    log_other = np.log(other) if other > 0 else -np.inf
    log_trans = np.full((k, k), log_other)
    np.fill_diagonal(log_trans, log_stay)

    delta = log_em[0].copy()
    back = np.zeros((n, k), dtype=int)
    for t in range(1, n):
        scores = delta[:, None] + log_trans
        best = scores.argmax(axis=0)
        back[t] = best
        delta = scores[best, np.arange(k)] + log_em[t]

    path = np.zeros(n, dtype=int)
    path[-1] = int(delta.argmax())
    for t in range(n - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return s[path]


# ---------------------------------------------------------------------------
# Loops
# ---------------------------------------------------------------------------

def apa_matrix(
    matrices: Any,
    anchors: Any,
    *,
    flank: int = 5,
) -> Any:
    """Aggregate peak analysis pileup.

    Stacks a ``(2*flank+1)`` square window centred on every anchor pair
    and averages. A real loop set produces a hot centre pixel; a shuffled
    one does not, which is why the paper shows an APA panel beside every
    cell-type-specific loop call.

    Parameters
    ----------
    matrices
        Mapping of chromosome name -> dense contact matrix.
    anchors
        Iterable of ``(chrom, bin_i, bin_j)`` triples, already binned at
        the matrices' resolution.
    flank
        Window half-width in bins.

    Returns
    -------
    ``(2*flank+1, 2*flank+1)`` mean pileup. All-NaN if no anchor fits
    inside its matrix with the full flank.
    """
    np = _np()
    w = 2 * flank + 1
    acc = np.zeros((w, w))
    count = 0
    for chrom, i, j in anchors:
        mat = matrices.get(chrom)
        if mat is None:
            continue
        n = mat.shape[0]
        i0, i1 = i - flank, i + flank + 1
        j0, j1 = j - flank, j + flank + 1
        # Skip rather than pad: a partially-truncated window would bias
        # the mean toward whichever corner survived.
        if i0 < 0 or j0 < 0 or i1 > n or j1 > n:
            continue
        acc += mat[i0:i1, j0:j1]
        count += 1
    if count == 0:
        return np.full((w, w), np.nan)
    return acc / count


def one_way_anova(groups: Any) -> "tuple[float, float]":
    """One-way ANOVA F statistic and p-value across groups.

    The test the paper applies per loop across cell types (loops with
    ``P < 0.05`` are called cell-type-specific). Implemented directly so
    the skill does not pull scipy in for a single F distribution — the
    survival function is taken from scipy when available and falls back
    to an incomplete-beta evaluation otherwise.

    Parameters
    ----------
    groups
        Sequence of 1-D arrays, one per group. Groups with fewer than
        one observation are dropped.

    Returns
    -------
    ``(F, p)``. ``(nan, nan)`` when fewer than two groups survive or
    there is no within-group variance to test against.
    """
    np = _np()
    arrs = [np.asarray(g, dtype=float) for g in groups]
    arrs = [a[np.isfinite(a)] for a in arrs]
    arrs = [a for a in arrs if a.size > 0]
    k = len(arrs)
    if k < 2:
        return float("nan"), float("nan")

    n_total = sum(a.size for a in arrs)
    if n_total <= k:
        return float("nan"), float("nan")

    grand = np.concatenate(arrs).mean()
    ss_between = sum(a.size * (a.mean() - grand) ** 2 for a in arrs)
    ss_within = sum(((a - a.mean()) ** 2).sum() for a in arrs)

    df_between = k - 1
    df_within = n_total - k
    if ss_within <= 0:
        # Zero within-group variance: the groups are either identical
        # (no effect) or perfectly separated (infinite F). Report the
        # separation honestly instead of dividing by zero.
        return (float("inf"), 0.0) if ss_between > 0 else (float("nan"), float("nan"))

    f = (ss_between / df_between) / (ss_within / df_within)
    return float(f), float(_f_sf(f, df_between, df_within))


def _f_sf(f: float, d1: int, d2: int) -> float:
    """P(F_{d1,d2} > f)."""
    if not (f == f) or f <= 0:
        return float("nan") if f != f else 1.0
    try:
        from scipy.stats import f as _fdist
        return float(_fdist.sf(f, d1, d2))
    except Exception:
        pass
    # Fallback: the F survival function is a regularised incomplete beta,
    # P(F > f) = I_{d2/(d2+d1 f)}(d2/2, d1/2).
    x = d2 / (d2 + d1 * f)
    return _betainc(d2 / 2.0, d1 / 2.0, x)


def _betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta I_x(a, b) via a continued fraction."""
    import math
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1 - x) - lbeta)
    # The continued fraction converges fast only for x < (a+1)/(a+b+2);
    # outside that, evaluate the symmetric case and complement.
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _betainc(b, a, 1 - x)

    f_, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30:
            c = 1e-30
        f_ *= c * d
        if abs(1.0 - c * d) < 1e-10:
            break
    return front * (f_ - 1.0) / a


# ---------------------------------------------------------------------------
# Spatial smoothing
# ---------------------------------------------------------------------------

def magic_smooth(
    X: Any,
    *,
    n_pca: int = 20,
    k: int = 10,
    t: int = 3,
) -> Any:
    """Markov-affinity graph imputation (MAGIC, van Dijk 2018).

    The paper uses MAGIC to denoise single-pixel CN ratios before
    plotting them in tissue space: neighbouring pixels in feature space
    share information, which turns a speckled per-pixel map into the
    spatial domains the tissue actually has.

    Three steps, as published: reduce with PCA, build an adaptive-kernel
    affinity graph, row-normalise it to a Markov matrix ``P``, and
    return ``P^t @ X``.

    Parameters
    ----------
    X
        ``(n_pixels, n_features)`` matrix, e.g. pixels by 5 Mb CN bins.
    n_pca
        Components used to build the graph. Capped at the data rank.
    k
        Neighbours defining the adaptive bandwidth.
    t
        Diffusion time; higher is smoother.

    Returns
    -------
    Smoothed array, same shape as ``X``.
    """
    np = _np()
    A = np.asarray(X, dtype=float)
    if A.ndim != 2:
        raise ValueError(f"expected a 2-D matrix, got shape {A.shape}")
    n = A.shape[0]
    if n < 2 or t <= 0:
        return A.copy()

    filled = np.nan_to_num(A, nan=0.0)

    # ── PCA ───────────────────────────────────────────────────────────
    centred = filled - filled.mean(axis=0, keepdims=True)
    comps = int(min(n_pca, min(centred.shape)))
    if comps >= 1:
        # economy SVD; we only need the left factor scaled by singulars
        U, S, _ = np.linalg.svd(centred, full_matrices=False)
        Y = U[:, :comps] * S[:comps]
    else:
        Y = centred

    # ── adaptive-bandwidth affinity ───────────────────────────────────
    sq = (Y ** 2).sum(axis=1)
    d2 = np.maximum(sq[:, None] + sq[None, :] - 2 * (Y @ Y.T), 0.0)
    kk = int(min(max(k, 1), n - 1))
    # Per-pixel bandwidth = distance to the k-th neighbour, so dense and
    # sparse regions of the tissue are smoothed on their own scale.
    knn_d2 = np.partition(d2, kk, axis=1)[:, kk]
    sigma = np.sqrt(np.maximum(knn_d2, 1e-12))
    W = np.exp(-d2 / (sigma[:, None] * sigma[None, :]))
    W = (W + W.T) / 2.0

    rowsum = W.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    P = W / rowsum

    out = filled
    for _ in range(int(t)):
        out = P @ out
    return out
