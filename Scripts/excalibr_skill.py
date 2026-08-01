"""exCALIBR — functional-assay calibration to ACMG/AMP evidence (IGVFagent port).

A clean-room Python reimplementation of **exCALIBR**
(github.com/rosstewart/exCALIBR, MIT, R. Stewart, Northeastern), which
implements the gene-based assay-calibration method of

    Zeiberg D, Tejura M, McEwen AE, Fayer S, Pejaver V, Rubin AF,
    Starita LM, Fowler DM, O'Donnell-Luria A, Radivojac P.
    "Gene-based calibration of high-throughput functional assays for
    clinical variant classification." bioRxiv 2025.04.29.651326.

The method turns a raw multiplexed-assay score (MAVE / VAMP-seq / SGE /
cell-fitness …) into **calibrated ACMG/AMP evidence strengths** (PS3 / BS3
at supporting / moderate / strong / very-strong), so an assay score can be
used directly in clinical variant classification instead of as a bare number.

Chain, faithful to upstream::

  1. Label variants into four samples from the scoreset table —
     0 P/LP (ClinVar), 1 B/LB (ClinVar), 2 population (gnomAD), 3 synonymous.
  2. Fit a **constrained skew-normal mixture** by EM: K components whose
     (skew, loc, scale) are *shared* across samples, each sample carrying its
     own mixing weights. The constraint forces the density ratio of adjacent
     components to be monotone (so the resulting LR+ is monotone in score);
     it is enforced parameter-by-parameter by binary search inside every
     M-step.
  3. **Bootstrap** the whole fit (per-sample resampling; best-of-N fits per
     bootstrap chosen on held-out log-likelihood).
  4. Estimate the **prior** P(pathogenic | population) per fit by EM against
     the population sample; take the median across bootstraps.
  5. Build **LR+(score) = f_P(score) / f_B(score)** per bootstrap, take the
     conservative percentile envelope across bootstraps.
  6. Solve for **Tavtigian's C (O_PVSt)** at that prior (Tavtigian 2018
     Bayesian ACMG framework), giving evidence thresholds C^(points/8), and
     convert the LR+ curve into **score ranges per evidence point value**.
  7. Optionally choose between a 2- and 3-component model by a paired
     Wilcoxon / 5th-percentile test on bootstrap validation likelihoods.

Reimplemented from the algorithm and the papers; no upstream source lines
copied. Deviations from upstream, all deliberate:

  * **Deterministic**: every random draw (bootstrap split, k-means seed,
    method-of-moments cut points, skew-sign table index) is derived from an
    explicit seed, so a rerun reproduces the calibration bit-for-bit. This
    matches IGVFagent's cross-backend consistency contract.
  * **stdlib parallelism** (`concurrent.futures`) instead of joblib; SLURM
    generation dropped in favour of a resumable in-process run — bootstraps
    are appended to a JSON-lines ledger with a progress heartbeat, so one
    long call finishes the job and `--resume` picks up where it stopped.
  * Empty sample categories are dropped at load time (with a name→column
    map) rather than index-shifted downstream.

CLI::

    igvfagent calibrate thresholds --prior 0.1
    igvfagent calibrate prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
    igvfagent calibrate run --table scores.csv --name MSH2_Jia_2021 \\
        --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
    igvfagent calibrate assign --calibration MSH2_2c_calibration.json \\
        --scores my_variants.csv
    igvfagent calibrate selftest

License: Apache-2.0 (upstream exCALIBR is MIT; method credited above).
Heavy deps imported lazily: numpy, scipy, matplotlib, scikit-learn
(k-means only — a numpy fallback is used when absent).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

warnings.filterwarnings("ignore")

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOC_DIR = ROOT / "Docs" / "Calibration"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "ASSAY_CALIBRATION_SKILL.md"

# Sample roster. Indices are part of the method's contract (upstream fixes
# them) — 0 P/LP, 1 B/LB, 2 population, 3 synonymous.
SAMPLE_NAMES = ["Pathogenic/Likely Pathogenic", "Benign/Likely Benign",
                "population", "Synonymous"]
PATHOGENIC, BENIGN, POPULATION, SYNONYMOUS = 0, 1, 2, 3

DEFAULT_POINT_VALUES = [1, 2, 3, 4, 5, 6, 7, 8]

# ACMG/AMP evidence strength per Tavtigian point value (8 points = very strong).
STRENGTH_BY_POINTS = {1: "supporting", 2: "moderate", 3: "moderate",
                      4: "strong", 5: "strong", 6: "strong", 7: "strong",
                      8: "very strong"}

CLINVAR_PATHOGENIC = {"Pathogenic", "Likely pathogenic",
                      "Pathogenic/Likely pathogenic"}
CLINVAR_BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}
# ClinVar review-status strings by star count (used for the quality gate).
ZERO_STAR = {"", "nan", "no assertion criteria provided", "no assertion provided",
             "no interpretation for the single variant",
             "no classification provided", "-"}
ONE_STAR = {"criteria provided, single submitter",
            "criteria provided, conflicting interpretations",
            "criteria provided, conflicting classifications"}
TWO_STAR = {"criteria provided, multiple submitters, no conflicts"}
THREE_STAR = {"reviewed by expert panel"}

_LOG_DENSITY_FLOOR = -7.0   # upstream's log-density support cutoff
_GRID = 1000                # constraint-check grid resolution
_SCORE_GRID = 10000         # calibration score grid resolution


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"calibrate_{label}_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def _lazy():
    """numpy + scipy.stats + scipy.special.logsumexp."""
    try:
        import numpy as np
        import scipy.stats as sps
        from scipy.special import logsumexp
        return np, sps, logsumexp
    except ImportError as exc:  # pragma: no cover
        sys.exit(f"calibrate needs numpy + scipy: {exc}")


# ─── Skew-normal mixture densities ──────────────────────────────────────────
#
# Canonical parameters are scipy's (a, loc, scale). The EM updates are closed
# form only in the "alternate" parameterisation (loc, Delta, Gamma) of
# Azzalini's stochastic representation, so we convert back and forth.

def canonical_to_alternate(a: float, loc: float, scale: float) -> Tuple[float, float, float]:
    delta = a / math.sqrt(1.0 + a * a)
    Delta = scale * delta
    Gamma = scale * scale - Delta * Delta
    return float(loc), float(Delta), float(Gamma)


def alternate_to_canonical(loc: float, Delta: float, Gamma: float) -> Tuple[float, float, float]:
    if Gamma <= 0:
        raise ZeroDivisionError(f"non-positive Gamma {Gamma}")
    a = math.copysign(math.sqrt(Delta * Delta / Gamma), Delta)
    if math.isinf(a) or math.isnan(a):
        raise ZeroDivisionError(f"invalid skew from Delta={Delta}, Gamma={Gamma}")
    scale = math.sqrt(Gamma + Delta * Delta)
    return float(a), float(loc), float(scale)


def log_joint_densities(np, sps, x, params, weights):
    """log(w_k * pdf_k(x)) for each component — shape (K, N)."""
    w = np.asarray(weights, dtype=float)
    log_pdfs = np.array([sps.skewnorm.logpdf(x, a, loc, scale)
                         for (a, loc, scale) in params])
    with np.errstate(divide="ignore"):
        log_w = np.log(w)
    log_w[w == 0] = -np.inf
    return log_w[:, None] + log_pdfs


def mixture_logpdf(np, sps, logsumexp, x, params, weights):
    """log of the mixture density — shape (N,)."""
    return logsumexp(log_joint_densities(np, sps, x, params, weights), axis=0)


def component_posteriors(np, sps, logsumexp, x, params, weights):
    """P(component k | x) under one sample's weights — shape (K, N)."""
    num = log_joint_densities(np, sps, np.asarray(x).ravel(), params, weights)
    den = logsumexp(num, axis=0)
    P = np.exp(num - den[None])
    P[np.isnan(P)] = 0.0
    return P


def total_loglik(np, sps, logsumexp, x, indicators, params, weights) -> float:
    """Sum over samples of the mixture log-likelihood under that sample's weights."""
    if params is None or weights is None:
        return -np.inf
    total = 0.0
    for s in range(indicators.shape[1]):
        mask = indicators[:, s]
        if not mask.any():
            continue
        total += float(logsumexp(
            log_joint_densities(np, sps, x[mask], params, weights[s]), axis=0).sum())
    return total


# ─── Monotone density-ratio constraint ──────────────────────────────────────

def constraint_violated(np, sps, params, xlims, tolerance: float = 0.0) -> bool:
    """True if any adjacent component pair has a non-decreasing log density ratio.

    Only the region where *both* components carry appreciable density
    (log pdf > -7) is examined, matching upstream: outside it the ratio is
    numerically meaningless.
    """
    grid = np.linspace(xlims[0], xlims[1], _GRID)
    try:
        log_pdfs = [sps.skewnorm.logpdf(grid, *p) for p in params]
    except (TypeError, ValueError):
        return True
    for i in range(len(log_pdfs) - 1):
        keep = (log_pdfs[i] > _LOG_DENSITY_FLOOR) & (log_pdfs[i + 1] > _LOG_DENSITY_FLOOR)
        if keep.sum() < 2:
            continue
        if np.any(np.diff(log_pdfs[i][keep] - log_pdfs[i + 1][keep]) > tolerance):
            return True
    return False


def _shrink_to_feasible(np, sps, params, xlims, max_trials: int = 300):
    """Shrink all scales by 5 % until the monotonicity constraint holds.

    Concentrating every component narrows the overlap regions where the log
    ratio can turn around; upstream uses the same scale-shrink repair.
    """
    params = sorted([tuple(float(v) for v in p) for p in params], key=lambda p: p[1])
    params = [(a, loc, max(scale, 1e-6)) for (a, loc, scale) in params]
    trial = 0
    while constraint_violated(np, sps, params, xlims) and trial < max_trials:
        params = [(a, loc, scale * 0.95) for (a, loc, scale) in params]
        trial += len(params)
        if min(p[2] for p in params) < 1e-6:
            break
    if constraint_violated(np, sps, params, xlims):
        return None
    return params


# ─── Initialisation ─────────────────────────────────────────────────────────

def _kmeans_1d(np, x, k, rng, iters: int = 100):
    """1-D k-means (numpy fallback when scikit-learn is unavailable)."""
    centres = rng.choice(x, size=k, replace=False).astype(float)
    for _ in range(iters):
        assign = np.abs(x[:, None] - centres[None, :]).argmin(axis=1)
        new = np.array([x[assign == j].mean() if (assign == j).any() else centres[j]
                        for j in range(k)])
        if np.allclose(new, centres):
            break
        centres = new
    return np.abs(x[:, None] - centres[None, :]).argmin(axis=1)


def kmeans_init(np, sps, x, k, rng, constrained: bool, attempts: int = 100):
    """Cluster the scores, fit a near-symmetric normal per cluster, repair."""
    for _ in range(attempts):
        try:
            from sklearn.cluster import KMeans
            km = KMeans(n_clusters=k, init="random", n_init=1,
                        random_state=int(rng.integers(0, 2 ** 31 - 1)))
            assign = km.fit_predict(x.reshape(-1, 1))
        except ImportError:
            assign = _kmeans_1d(np, x, k, rng)
        params = []
        ok = True
        for j in range(k):
            xj = x[assign == j]
            if len(xj) < 2:
                ok = False
                break
            loc, scale = float(xj.mean()), float(xj.std())
            if not np.isfinite(scale) or scale <= 0:
                ok = False
                break
            params.append((float(rng.uniform(-0.25, 0.25)), loc, scale))
        if not ok:
            continue
        params = sorted(params, key=lambda p: p[1])
        if constrained:
            params = _shrink_to_feasible(np, sps, params, (x.min(), x.max()))
            if params is None:
                continue
        return params
    return None


def _sn_method_of_moments(np, sps, x, rng):
    """Skew-normal method-of-moments estimator for one slice of the scores."""
    m1, m2, m3 = float(np.mean(x)), float(np.var(x)), float(sps.skew(x))
    if m2 < 1e-10:
        return None
    a1 = math.sqrt(2.0 / math.pi)
    b1 = (4.0 / math.pi - 1.0) / a1
    if m3 == 0:
        return float(rng.uniform(-0.25, 0.25)), m1, max(math.sqrt(m2), 1e-6)
    delta = math.copysign(1.0, m3) / math.sqrt(a1 ** 2 + m2 * (b1 / abs(m3)) ** (2.0 / 3.0))
    if math.isnan(delta) or abs(delta) >= 0.99:
        return float(rng.uniform(-0.25, 0.25)), m1, max(math.sqrt(m2), 1e-6)
    denom = 1.0 - a1 ** 2 * delta ** 2
    if denom <= 1e-10:
        return m3, m1, max(math.sqrt(m2), 1e-6)
    sigma = max(math.sqrt(m2 / denom), 1e-6)
    mu = m1 - a1 * delta * sigma
    if not all(map(math.isfinite, (mu, sigma, m3))):
        return None
    return m3, mu, sigma


def mom_init(np, sps, x, k, rng, constrained: bool, lambda_index: int,
             attempts: int = 200):
    """Slice the scores at (jittered) percentiles, method-of-moments each slice.

    The skew of component j is set to the j-th entry of the ±1 sign table
    selected by `lambda_index`, and every component is given the widest
    slice's scale — upstream's way of seeding a mixture that can satisfy the
    monotone-ratio constraint.
    """
    import itertools
    signs = list(itertools.product([-1, 1], repeat=k))[lambda_index % (2 ** k)]
    min_samples = max(10, int(0.05 * len(x)))
    q25, q75 = np.percentile(x, [25, 75])
    iqr = float(q75 - q25)
    for _ in range(attempts):
        if rng.random() < 0.7:
            base = np.linspace(0, 100, k + 1)[1:-1]
            jitter = rng.normal(0, max(iqr, 1e-9) * 0.1, len(base))
            cuts = np.percentile(x, np.sort(np.clip(base + jitter, 1, 99)))
        else:
            lo, hi = np.percentile(x, [5, 95])
            cuts = np.sort(rng.uniform(lo, hi, k - 1))
        params, ok = [], True
        for j in range(k):
            if j == 0:
                xj = x[x <= cuts[0]]
            elif j == k - 1:
                xj = x[x > cuts[-1]]
            else:
                xj = x[(x > cuts[j - 1]) & (x <= cuts[j])]
            if len(xj) < min_samples:
                ok = False
                break
            est = _sn_method_of_moments(np, sps, xj, rng)
            if est is None:
                ok = False
                break
            params.append((float(signs[j]), est[1], est[2]))
        if not ok:
            continue
        widest = max(p[2] for p in params)
        params = [(a, loc, widest) for (a, loc, _) in params]
        if constrained:
            params = _shrink_to_feasible(np, sps, params, (x.min(), x.max()))
            if params is None:
                continue
        return params
    return None


# ─── EM ─────────────────────────────────────────────────────────────────────

def _trunc_norm_moments(np, sps, mu, sigma):
    """First two moments of N(mu, sigma) truncated to the positive half-line."""
    z = mu / sigma
    cdf = sps.norm.cdf(z)
    pdf = sps.norm.pdf(z)
    ratio = np.zeros_like(pdf)
    zero = cdf == 0
    ratio[~zero] = pdf[~zero] / cdf[~zero]
    ratio[zero] = np.abs(z[zero])
    m1 = mu + sigma * ratio
    m2 = mu ** 2 + sigma ** 2 + sigma * mu * ratio
    return m1, m2


def _component_moments(np, sps, x, params_k):
    a, loc, scale = params_k
    delta = a / math.sqrt(1.0 + a * a)
    mu = delta / scale * (x - loc)
    sigma = math.sqrt(max(1.0 - delta * delta, 1e-300))
    return _trunc_norm_moments(np, sps, mu, sigma)


def sample_weights(np, sps, logsumexp, x, indicators, params, current_weights):
    """M-step for the per-sample mixing weights."""
    updated = np.zeros_like(current_weights)
    for s in range(current_weights.shape[0]):
        xs = x[indicators[:, s]]
        if len(xs) == 0:
            updated[s] = current_weights[s]
            continue
        posts = component_posteriors(np, sps, logsumexp, xs, params, current_weights[s])
        w = posts.mean(axis=1)
        if np.isnan(w).any():
            raise ValueError("NaN in weight update")
        updated[s] = w
    return updated


def _binary_search(np, sps, candidate, current_params, comp_i, param_i, xlims,
                   tol: float = 1e-4):
    """Largest step toward `candidate` that keeps the constraint satisfied.

    `current_params` is feasible, so the current value of the parameter is a
    valid lower bound; bisect between it and the unconstrained M-step value.
    """
    if constraint_violated(np, sps, current_params, xlims):
        raise ValueError("constraint already violated before binary search")
    alt = [list(canonical_to_alternate(*p)) for p in current_params]
    lower = alt[comp_i][param_i]
    upper = float(candidate)
    while abs(upper - lower) > tol:
        mid = (upper + lower) / 2.0
        trial = [list(p) for p in alt]
        trial[comp_i][param_i] = mid
        try:
            canon = [alternate_to_canonical(*p) for p in trial]
        except ZeroDivisionError:
            upper = mid
            continue
        if constraint_violated(np, sps, canon, xlims):
            upper = mid
        else:
            lower = mid
    return lower


def em_iteration(np, sps, logsumexp, x, indicators, params, weights,
                 constrained: bool, xlims):
    """One constrained EM iteration. Components are updated in order, each
    parameter (loc → Delta → Gamma) projected back into the feasible set."""
    if constrained and constraint_violated(np, sps, params, xlims):
        raise ValueError("constraint violated at start of EM iteration")
    K = len(params)
    # E-step: responsibilities under each observation's own sample weights.
    resp = np.zeros((K, len(x)))
    for s in range(indicators.shape[1]):
        mask = indicators[:, s]
        if not mask.any():
            continue
        resp[:, mask] = component_posteriors(np, sps, logsumexp, x[mask], params,
                                             weights[s])

    updated: List[Any] = [None] * K
    for k in range(K):
        r = resp[k]
        r_sum = r.sum()
        if r_sum <= 0 or not np.isfinite(r_sum):
            raise ValueError(f"degenerate responsibility for component {k}")
        v, w = _component_moments(np, sps, x, params[k])
        _, Delta_old, _ = canonical_to_alternate(*params[k])

        # loc update
        loc_cand = float(((x - v * Delta_old) * r).sum() / r_sum)
        bsearch = [updated[j] if j < k else params[j] for j in range(K)]
        loc_new = (_binary_search(np, sps, loc_cand, bsearch, k, 0, xlims)
                   if constrained else loc_cand)

        # Delta update (given the new loc)
        denom = float((w * r).sum())
        if denom == 0:
            raise ValueError("zero Delta denominator")
        Delta_cand = float((v * (x - loc_new) * r).sum() / denom)
        if constrained:
            bsearch = []
            for j in range(K):
                if j < k:
                    bsearch.append(updated[j])
                elif j > k:
                    bsearch.append(params[j])
                else:
                    _, D, G = canonical_to_alternate(*params[j])
                    bsearch.append(alternate_to_canonical(loc_new, D, G))
            Delta_new = _binary_search(np, sps, Delta_cand, bsearch, k, 1, xlims)
        else:
            Delta_new = Delta_cand

        # Gamma update (given the new loc + Delta)
        g = ((x - loc_new) ** 2 - 2 * Delta_new * v * (x - loc_new)
             + Delta_new ** 2 * w)
        Gamma_cand = float((g * r).sum() / r_sum)
        if constrained:
            bsearch = []
            for j in range(K):
                if j < k:
                    bsearch.append(updated[j])
                elif j > k:
                    bsearch.append(params[j])
                else:
                    _, _, G = canonical_to_alternate(*params[j])
                    bsearch.append(alternate_to_canonical(loc_new, Delta_new, G))
            Gamma_new = _binary_search(np, sps, Gamma_cand, bsearch, k, 2, xlims)
        else:
            Gamma_new = Gamma_cand

        updated[k] = alternate_to_canonical(loc_new, Delta_new, Gamma_new)
        if constrained and constraint_violated(
                np, sps, [*updated[:k + 1], *params[k + 1:]], xlims):
            raise ValueError(f"constraint violated after updating component {k}")

    return updated, sample_weights(np, sps, logsumexp, x, indicators, updated, weights)


def single_fit(np, sps, logsumexp, x, indicators, k: int, constrained: bool,
               init_method: str, seed: int, lambda_index: int = 0,
               max_em_iters: int = 10000, tol: float = 1e-8,
               initial_params=None, initial_weights=None) -> Dict[str, Any]:
    """One constrained skew-normal mixture fit. Never raises: a failed fit
    comes back with `loglik = -inf` so the caller can simply take the best."""
    MIN_SCALE = 1e-100
    rng = np.random.default_rng(seed)
    xlims = (float(x.min()), float(x.max()))
    S = indicators.shape[1]

    if initial_params is not None and initial_weights is not None:
        params = [tuple(map(float, p)) for p in initial_params]
        W = np.array(initial_weights, dtype=float)
    else:
        params = None
        if init_method == "method_of_moments":
            params = mom_init(np, sps, x, k, rng, constrained, lambda_index)
        if params is None:
            params = kmeans_init(np, sps, x, k, rng, constrained)
        if params is None:
            return {"component_params": None, "weights": None,
                    "loglik": -np.inf, "xlims": xlims, "n_iters": 0,
                    "status": "init-failed"}
        W = sample_weights(np, sps, logsumexp, x, indicators, params,
                           np.ones((S, k)) / k)

    lls = [total_loglik(np, sps, logsumexp, x, indicators, params, W) / len(x)]
    status = "converged"
    n_iters = 0
    try:
        for i in range(max_em_iters):
            params, W = em_iteration(np, sps, logsumexp, x, indicators, params, W,
                                     constrained, xlims)
            params = [(a, loc, max(scale, MIN_SCALE)) for (a, loc, scale) in params]
            if np.isnan(W).any() or not np.isfinite(np.array(params)).all():
                raise ValueError(f"NaN parameters at iteration {i}")
            lls.append(total_loglik(np, sps, logsumexp, x, indicators, params, W) / len(x))
            n_iters = i + 1
            if i > 1 and lls[-1] < lls[-2] and (lls[-2] - lls[-1]) > 1e-13:
                raise ValueError(f"likelihood decreased at iteration {i}")
            # Upstream takes one EM step before entering its loop and only then
            # allows early stopping, so the earliest exit is after three steps;
            # match that so both implementations stop at the same iterate.
            if i >= 2 and abs(lls[-1] - lls[-2]) / max(abs(lls[-2]), 1e-300) < tol:
                break
        else:
            status = "max-iters"
        if constrained and constraint_violated(np, sps, params, xlims):
            raise ValueError("final parameters violate the density constraint")
    except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
        return {"component_params": params, "weights": W, "loglik": -np.inf,
                "xlims": xlims, "n_iters": n_iters, "status": f"failed: {exc}"}

    return {"component_params": [tuple(map(float, p)) for p in params],
            "weights": W, "loglik": float(lls[-1]), "xlims": xlims,
            "n_iters": n_iters, "status": status}


# ─── Bootstrap ──────────────────────────────────────────────────────────────

def per_sample_bootstrap(np, indicators, seed: int):
    """Resample each sample with replacement; the out-of-bag rows are the
    validation set. Retried until every resampled sample leaves some out."""
    rng = np.random.default_rng(seed)
    train, val = [], []
    for s in range(indicators.shape[1]):
        idx = np.where(indicators[:, s])[0]
        if len(idx) == 0:
            continue
        if len(idx) == 1:
            train.append(idx)
            continue
        oob = np.array([], dtype=int)
        picked = idx
        for _ in range(100):
            picked = rng.choice(idx, size=len(idx), replace=True)
            oob = np.setdiff1d(idx, picked)
            if len(oob):
                break
        if not len(oob):
            raise ValueError("failed to generate a bootstrap split")
        train.append(picked)
        val.append(oob)
    return np.concatenate(train), (np.concatenate(val) if val
                                   else np.array([], dtype=int))


def _bootstrap_worker(payload):
    """Run every fit of one bootstrap iteration; return the best per K.

    Module-level and self-contained so it can be sent to a process pool.
    """
    np, sps, logsumexp = _lazy()
    scores = np.asarray(payload["scores"], dtype=float)
    indicators = np.asarray(payload["indicators"], dtype=bool)
    b = payload["bootstrap"]
    n_fits = payload["fits_per_bootstrap"]
    out = {"bootstrap": b}
    try:
        tr, va = per_sample_bootstrap(np, indicators, b)
    except ValueError as exc:
        return {"bootstrap": b, "error": str(exc)}
    x_tr, ind_tr = scores[tr], indicators[tr]
    x_va, ind_va = scores[va], indicators[va]
    # Deterministic per-bootstrap choice of initialisation strategy.
    pick = np.random.default_rng(b)
    methods = pick.choice(["kmeans", "method_of_moments"], size=n_fits)
    for k in payload["components"]:
        best, best_val = None, -np.inf
        for i in range(n_fits):
            fit = single_fit(np, sps, logsumexp, x_tr, ind_tr, k,
                             payload["constrained"], str(methods[i]),
                             seed=b * 100003 + i * 97 + k,
                             lambda_index=i % (2 ** k),
                             max_em_iters=payload["max_em_iters"])
            if not np.isfinite(fit["loglik"]):
                continue
            val_ll = (total_loglik(np, sps, logsumexp, x_va, ind_va,
                                   fit["component_params"], fit["weights"])
                      / len(x_va)) if len(x_va) else fit["loglik"]
            if np.isfinite(val_ll) and val_ll > best_val:
                best, best_val = fit, float(val_ll)
        if best is not None:
            out[f"{k}c"] = {
                "component_params": [list(p) for p in best["component_params"]],
                "weights": np.asarray(best["weights"]).tolist(),
                "xlims": list(best["xlims"]),
                "train_ll": best["loglik"],
                "val_ll": best_val,
                "n_iters": best["n_iters"],
            }
    return out


def run_bootstraps(np, scores, indicators, components: Sequence[int],
                   n_bootstraps: int, fits_per_bootstrap: int, jobs: int,
                   ledger: Optional[Path], constrained: bool = True,
                   max_em_iters: int = 10000,
                   heartbeat: int = 10) -> List[Dict[str, Any]]:
    """Fit every bootstrap, appending to a resumable JSON-lines ledger.

    One call runs to completion regardless of duration; progress is logged
    every `heartbeat` bootstraps and each finished bootstrap is durable, so
    `--resume` restarts from the ledger instead of from zero.
    """
    import concurrent.futures as cf

    done: Dict[int, Dict[str, Any]] = {}
    if ledger is not None and ledger.exists():
        with ledger.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                done[int(rec["bootstrap"])] = rec
        if done:
            logging.info("resuming: %d/%d bootstraps already in %s",
                         len(done), n_bootstraps, ledger)

    todo = [b for b in range(n_bootstraps) if b not in done]
    payloads = [{"scores": np.asarray(scores).tolist(),
                 "indicators": np.asarray(indicators).tolist(),
                 "bootstrap": b, "components": list(components),
                 "fits_per_bootstrap": fits_per_bootstrap,
                 "constrained": constrained, "max_em_iters": max_em_iters}
                for b in todo]

    t0 = time.time()
    completed = 0
    fh = ledger.open("a") if ledger is not None else None
    try:
        if jobs == 1 or len(payloads) <= 1:
            results_iter = (_bootstrap_worker(p) for p in payloads)
            for rec in results_iter:
                done[int(rec["bootstrap"])] = rec
                completed += 1
                if fh:
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                if completed % heartbeat == 0 or completed == len(payloads):
                    _heartbeat(completed, len(payloads), t0)
        else:
            with cf.ProcessPoolExecutor(max_workers=jobs) as pool:
                futures = {pool.submit(_bootstrap_worker, p): p["bootstrap"]
                           for p in payloads}
                for fut in cf.as_completed(futures):
                    try:
                        rec = fut.result()
                    except Exception as exc:      # keep the run alive
                        rec = {"bootstrap": futures[fut], "error": str(exc)}
                    done[int(rec["bootstrap"])] = rec
                    completed += 1
                    if fh:
                        fh.write(json.dumps(rec) + "\n")
                        fh.flush()
                    if completed % heartbeat == 0 or completed == len(payloads):
                        _heartbeat(completed, len(payloads), t0)
    finally:
        if fh:
            fh.close()
    return [done[b] for b in sorted(done)]


def _heartbeat(done: int, total: int, t0: float) -> None:
    elapsed = time.time() - t0
    rate = done / elapsed if elapsed > 0 else 0.0
    eta = (total - done) / rate if rate > 0 else float("nan")
    logging.info("bootstraps %d/%d (%.1f%%) — %.1fs elapsed, ~%.0fs remaining",
                 done, total, 100.0 * done / max(total, 1), elapsed, eta)


# ─── Prior estimation ───────────────────────────────────────────────────────

def benign_weights(np, weights, idx: Dict[str, Optional[int]], benign_method: str):
    """Effective benign mixing weights under the chosen benign definition."""
    b, s = idx.get("benign"), idx.get("synonymous")
    if benign_method == "avg" and b is not None and s is not None:
        return (np.asarray(weights[b], dtype=float)
                + np.asarray(weights[s], dtype=float)) / 2.0
    if benign_method == "synonymous" and s is not None:
        return np.asarray(weights[s], dtype=float)
    if b is not None:
        return np.asarray(weights[b], dtype=float)
    return np.asarray(weights[s], dtype=float)


def fit_prior(np, sps, logsumexp, pop_scores, params, weights,
              idx: Dict[str, Optional[int]], benign_method: str,
              tol: float = 1e-6, max_steps: int = 10000) -> float:
    """EM estimate of P(pathogenic) among the population sample.

    Standard case (both P and B densities available) is the Saerens/Latinne
    mixture-proportion EM; with only one labelled class it degenerates to the
    positive/negative-unlabelled variant.
    """
    p_idx, b_idx, s_idx = idx.get("pathogenic"), idx.get("benign"), idx.get("synonymous")
    if len(pop_scores) == 0:
        return float("nan")
    pop_density = np.exp(mixture_logpdf(np, sps, logsumexp, pop_scores, params,
                                        weights[idx["population"]]))
    fp = (np.exp(mixture_logpdf(np, sps, logsumexp, pop_scores, params,
                                weights[p_idx])) if p_idx is not None else None)
    fb = (np.exp(mixture_logpdf(np, sps, logsumexp, pop_scores, params,
                                benign_weights(np, weights, idx, benign_method)))
          if (b_idx is not None or s_idx is not None) else None)

    if fp is not None and fb is not None:
        mode, alpha = "standard", 0.5
    elif fp is not None:
        mode, alpha = "positive_unlabeled", 0.1
    elif fb is not None:
        mode, alpha = "negative_unlabeled", 0.9
    else:
        raise ValueError("need a pathogenic or a benign sample to estimate the prior")

    prev = alpha
    for _ in range(max_steps):
        with np.errstate(divide="ignore", invalid="ignore", over="ignore",
                         under="ignore"):
            if mode == "standard":
                post = 1.0 / (1.0 + (1.0 - alpha) / alpha * fb / fp)
            elif mode == "positive_unlabeled":
                post = alpha * fp / pop_density
            else:
                post = alpha * fb / pop_density
        post = np.clip(post, 0.0, 1.0)
        new = float(np.nanmean(post))
        converged = abs(new - prev) < tol
        prev, alpha = alpha, new
        if converged or alpha <= 0 or alpha >= 1:
            break
    if mode == "negative_unlabeled":
        alpha = 1.0 - alpha
    if not (0.0 < alpha < 1.0) or math.isnan(alpha):
        return float("nan")
    return float(alpha)


def class_log_densities(np, sps, logsumexp, grid, params, weights,
                        idx: Dict[str, Optional[int]], benign_method: str,
                        prior: float):
    """(log f_P, log f_B) on `grid`. A missing class is unmixed out of the
    population density with the prior: f_pop = α f_P + (1-α) f_B."""
    if math.isnan(prior) or not (0.0 < prior < 1.0):
        return None, None
    log_pop = mixture_logpdf(np, sps, logsumexp, grid, params,
                             weights[idx["population"]])
    pop = np.exp(log_pop)
    have_p = idx.get("pathogenic") is not None
    have_b = idx.get("benign") is not None or idx.get("synonymous") is not None

    if have_p:
        log_fp = mixture_logpdf(np, sps, logsumexp, grid, params,
                                weights[idx["pathogenic"]])
    if have_b:
        log_fb = mixture_logpdf(np, sps, logsumexp, grid, params,
                                benign_weights(np, weights, idx, benign_method))
    if not have_p:
        fp = np.maximum((pop - (1 - prior) * np.exp(log_fb)) / prior, pop * 1e-10)
        log_fp = np.log(fp)
    if not have_b:
        fb = np.maximum((pop - prior * np.exp(log_fp)) / (1 - prior), pop * 1e-10)
        log_fb = np.log(fb)
    return log_fp, log_fb


# ─── Tavtigian constant + evidence thresholds ───────────────────────────────

def _posterior_vec(np, lr_plus, prior: float):
    """Posterior probability of pathogenicity from a positive likelihood ratio."""
    return lr_plus * prior / ((lr_plus - 1) * prior + 1)


def tavtigian_constant(np, prior: float, original: bool = False,
                       strict: bool = False, c_max: int = 30000) -> int:
    """Solve for C = O_PVSt: the odds of pathogenicity of one very-strong line
    of evidence that best reproduces the ACMG/AMP combining rules at `prior`.

    Tavtigian et al., Genet Med 2018: each evidence strength contributes
    C^(1/8), C^(1/4), C^(1/2), C^1 (supporting → very strong). We score every
    candidate C by how many of the P / LP / B / LB rule combinations land
    outside their posterior bands, and take the argmin (upstream's criterion).
    """
    fracs = np.array([2.0 ** -3, 2.0 ** -2, 2.0 ** -1, 1.0])
    C = np.arange(1, c_max + 1, dtype=float)

    def post(counts):
        return _posterior_vec(np, C ** float(np.dot(counts, fracs)), prior)

    def post_b(counts):
        return _posterior_vec(np, C ** float(np.dot(counts, -fracs)), prior)

    # Pathogenic rules (counts of supporting / moderate / strong / very strong).
    path = [[0, 0, 1, 1], [0, 2, 0, 1], [1, 1, 0, 1], [2, 0, 0, 1],
            [0, 3, 1, 0], [2, 2, 1, 0], [4, 1, 1, 0]]
    path.append([0, 0, 2, 0] if original else [0, 1, 0, 1])
    lp = [[0, 1, 1, 0], [2, 0, 1, 0], [0, 3, 0, 0], [2, 2, 0, 0], [4, 1, 0, 0]]
    lp.append([0, 1, 0, 1] if original else [0, 0, 2, 0])

    path_post = np.round(np.stack([post(c) for c in path], axis=1), 3)
    lp_post = np.round(np.stack([post(c) for c in lp], axis=1), 3)
    ben_post = np.round(np.stack([post_b([0, 0, 2, 0])], axis=1), 3)
    lb_post = np.round(np.stack([post_b([1, 0, 1, 0]), post_b([2, 0, 0, 0])],
                                axis=1), 3)

    fails = (path_post < 0.99).sum(axis=1)
    lp_mask = lp_post < 0.90
    if strict:
        lp_mask |= lp_post > 0.99
    fails = fails + lp_mask.sum(axis=1) + (ben_post > 0.01).sum(axis=1)
    lb_mask = lb_post > 0.10
    if strict:
        lb_mask |= lb_post < 0.01
    fails = fails + lb_mask.sum(axis=1)
    return int(C[int(np.argmin(fails))])


def thresholds_from_prior(np, prior: float, point_values: Sequence[int],
                          **kwargs):
    """LR+ thresholds for each evidence point value: C^(points / max_points)."""
    C = tavtigian_constant(np, prior, **kwargs)
    pts = np.asarray(point_values, dtype=float)
    lr_p = C ** (pts / len(point_values))
    return lr_p, 1.0 / lr_p, C


# ─── LR+ → score ranges per evidence point value ────────────────────────────

def _assign_points(lr: float, tau, point_values, benign: bool) -> int:
    """Highest point value whose log-LR+ threshold band contains `lr`."""
    for i, t in enumerate(tau):
        if benign:
            if lr <= t and (i == len(tau) - 1 or lr > tau[i + 1]):
                return int(point_values[i])
        else:
            if lr >= t and (i == len(tau) - 1 or lr < tau[i + 1]):
                return int(point_values[i])
    return 0


def point_ranges_from_lr(np, scores, log_lr, tau, point_values, benign: bool):
    """Walk the score grid, recording contiguous ranges per point value."""
    pts = [-int(p) for p in point_values] if benign else [int(p) for p in point_values]
    ranges: Dict[int, List[List[float]]] = {p: [] for p in pts}
    open_at, open_point = float("nan"), float("nan")
    for s, lr in zip(scores, log_lr):
        p = _assign_points(lr, tau, pts, benign)
        if p != open_point:
            if open_point == open_point and open_point != 0:   # not NaN, not 0
                ranges[int(open_point)].append(sorted([float(open_at), float(s)]))
            open_at, open_point = s, p
    if open_point == open_point and open_point != 0:
        ranges[int(open_point)].append(sorted([float(open_at), float(scores[-1])]))
    return ranges


def calculate_score_ranges(np, log_lr_low, log_lr_high, prior, scores,
                           point_values):
    """Pathogenic ranges from the conservative (low) LR+ curve, benign ranges
    from the conservative (high) curve, both against the same prior's C."""
    lr_p, lr_b, C = thresholds_from_prior(np, prior, point_values)
    p = point_ranges_from_lr(np, scores, log_lr_low, np.log(lr_p), point_values, False)
    b = point_ranges_from_lr(np, scores, log_lr_high, np.log(lr_b), point_values, True)
    return p, b, C


def conservative_point_ranges(np, per_fit_p, per_fit_b, point_values,
                              score_range, flipped: bool):
    """Thresholds from the 5th-percentile of each bootstrap's *own* boundary.

    The `--no-median-prior` alternative to the LR+ envelope: for each evidence
    strength take the least-extreme boundary that 95 % of bootstraps reach, so a
    strength is only granted where nearly every bootstrap grants it. A bootstrap
    that never reaches a strength contributes ±inf, which drags the percentile
    the conservative way.
    """
    p_pct, b_pct = (5, 95) if not flipped else (95, 5)
    p_edge = max if not flipped else min
    b_edge = min if not flipped else max
    p_inf = -np.inf if not flipped else np.inf
    b_inf = np.inf if not flipped else -np.inf
    top = max(point_values)

    thresholds: Dict[int, float] = {}
    for p in point_values:
        thresholds[int(p)] = float(np.nanpercentile(
            [p_edge(r[int(p)]) if len(r.get(int(p), [])) else p_inf
             for r in per_fit_p], p_pct))
        thresholds[int(-p)] = float(np.nanpercentile(
            [b_edge(r[int(-p)]) if len(r.get(int(-p), [])) else b_inf
             for r in per_fit_b], b_pct))

    lo_lim, hi_lim = float(score_range[0]), float(score_range[-1])
    ranges: Dict[int, List[List[float]]] = {}
    for p, thr in thresholds.items():
        if math.isnan(thr) or math.isinf(thr):
            ranges[p] = []
            continue
        pathogenic_side = (p > 0 and not flipped) or (p < 0 and flipped)
        neighbour = thresholds.get(p + 1 if p > 0 else p - 1, float("nan"))
        if pathogenic_side:
            inner = lo_lim if math.isnan(neighbour) else neighbour
            ranges[p] = [sorted([float(inner), thr])] if abs(p) != top \
                else [sorted([lo_lim, thr])]
        else:
            inner = hi_lim if math.isnan(neighbour) else neighbour
            ranges[p] = [sorted([thr, float(inner)])] if abs(p) != top \
                else [sorted([thr, hi_lim])]
    return ranges


def enforce_monotonicity(point_ranges, point_values, flipped: bool) -> None:
    """Keep one range per point value, on the side the assay's scale implies.

    Upstream's 'liberal' rule: for pathogenic evidence keep the outermost
    range in the pathogenic direction (and mirror for benign), which makes the
    thresholds monotone in score without discarding evidence.
    """
    for p in point_values:
        if point_ranges.get(p):
            point_ranges[p] = ([point_ranges[p][-1]] if not flipped
                               else [point_ranges[p][0]])
        if point_ranges.get(-p):
            point_ranges[-p] = ([point_ranges[-p][0]] if not flipped
                                else [point_ranges[-p][-1]])


def extend_to_xlims(point_ranges, point_values, score_range, flipped: bool) -> None:
    """The strongest non-empty strength on each side runs to the score limit."""
    left, right = float(score_range[0]), float(score_range[-1])
    for p in point_values:
        if point_ranges.get(p):
            if not any(point_ranges.get(p + j) for j in range(1, max(point_values) + 1)):
                point_ranges[p] = ([[left, point_ranges[p][-1][-1]]] if not flipped
                                   else [[point_ranges[p][0][0], right]])
        if point_ranges.get(-p):
            if not any(point_ranges.get(-p - j) for j in range(1, max(point_values) + 1)):
                point_ranges[-p] = ([[left, point_ranges[-p][-1][-1]]] if flipped
                                    else [[point_ranges[-p][0][0], right]])


# ─── Model selection ────────────────────────────────────────────────────────

def paired_bootstrap_test(np, val_low, val_high, k_low: int, k_high: int,
                          alpha: float = 0.05) -> Dict[str, Any]:
    """Paired Wilcoxon signed-rank on per-bootstrap validation likelihoods,
    plus the conservative 5th-percentile rule (95 % of bootstraps improve)."""
    from scipy import stats
    d = np.asarray(val_high, dtype=float) - np.asarray(val_low, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) < 2:
        return {"selected_k": k_low, "conservative_k": k_low, "p_value": float("nan"),
                "n_samples": int(len(d)), "method": "insufficient-pairs",
                "k_low": k_low, "k_high": k_high}
    stat, p = stats.wilcoxon(d, alternative="greater")
    fifth = float(np.percentile(d, 5))
    lo, hi = np.percentile(d, [2.5, 97.5])
    return {"selected_k": int(k_high if p < alpha else k_low),
            "conservative_k": int(k_high if fifth > 0 else k_low),
            "p_value": float(p), "statistic": float(stat),
            "mean_diff": float(d.mean()), "median_diff": float(np.median(d)),
            "std_diff": float(d.std()), "fifth_percentile": fifth,
            "ci_95": [float(lo), float(hi)], "n_samples": int(len(d)),
            "method": "wilcoxon_paired", "k_low": int(k_low), "k_high": int(k_high)}


# ─── Scoreset loading ───────────────────────────────────────────────────────

class Scoreset:
    """Scores plus a one-hot sample-membership matrix over non-empty samples."""

    def __init__(self, np, scores, assignments, names: Sequence[str],
                 ids: Optional[Sequence[str]] = None):
        scores = np.asarray(scores, dtype=float)
        assignments = np.asarray(assignments, dtype=bool)
        keep = assignments.any(axis=1) & ~np.isnan(scores)
        self.scores = scores[keep]
        assignments = assignments[keep]
        self.ids = ([i for i, k in zip(ids, keep) if k] if ids is not None
                    else [str(i) for i in range(len(self.scores))])
        # Drop empty sample columns; keep a name → column map.
        nonempty = [j for j in range(assignments.shape[1]) if assignments[:, j].any()]
        self.assignments = assignments[:, nonempty]
        self.names = [names[j] for j in nonempty]
        self.counts = self.assignments.sum(axis=0).tolist()
        self.index = {}
        for key, name in (("pathogenic", SAMPLE_NAMES[PATHOGENIC]),
                          ("benign", SAMPLE_NAMES[BENIGN]),
                          ("population", SAMPLE_NAMES[POPULATION]),
                          ("synonymous", SAMPLE_NAMES[SYNONYMOUS])):
            self.index[key] = self.names.index(name) if name in self.names else None

    def one_hot(self, np, seed: int = 0):
        """Assign each multiply-labelled variant to exactly one sample.

        A variant can be both ClinVar-labelled and in gnomAD; the mixture
        needs disjoint samples, so ties are broken deterministically while
        guaranteeing no sample is emptied.
        """
        rng = np.random.default_rng(seed)
        A = self.assignments
        for _ in range(1000):
            out = np.zeros_like(A)
            for i in range(A.shape[0]):
                cols = np.where(A[i])[0]
                if len(cols):
                    out[i, rng.choice(cols)] = True
            if out.any(axis=0).all():
                return out
        raise ValueError("could not build disjoint samples (a sample stays empty)")

    def __repr__(self):
        return (f"<Scoreset n={len(self.scores)} "
                + ", ".join(f"{n}={c}" for n, c in zip(self.names, self.counts)) + ">")


def _stars_sufficient(review_status: str, min_star: int) -> bool:
    rs = (review_status or "").strip().lower()
    if min_star <= 0:
        return True
    bad = set(ZERO_STAR)
    if min_star >= 2:
        bad |= ONE_STAR
    if min_star >= 3:
        bad |= TWO_STAR
    if min_star >= 4:
        bad |= THREE_STAR
    return rs not in {b.lower() for b in bad}


def _open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def load_pillar(np, path: Path, dataset: Optional[str] = None,
                clinvar_release: str = "2025", min_clinvar_star: int = 1,
                population_type: str = "gnomAD") -> Scoreset:
    """Load an IGVF / Pillar-project-format scoreset table.

    The 75-column format carries the assay score (`auth_reported_score`) plus
    everything needed to label it: ClinVar significance + review status for a
    given release, gnomAD MAF, consequence, SpliceAI scores and a QC flag.
    Labelling reproduces upstream: rows flagged `*` are dropped; unless the
    assay itself measures splicing, splice-region variants and anything with a
    SpliceAI delta ≥ 0.2 are dropped; a variant that is synonymous goes to the
    synonymous sample exclusively; otherwise gnomAD membership, ClinVar P/LP
    and ClinVar B/LB are independent memberships.
    """
    sig_col = f"clinvar_sig_{clinvar_release}"
    star_col = f"clinvar_star_{clinvar_release}"
    rows: List[Dict[str, str]] = []
    with _open_maybe_gz(path) as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in ("auth_reported_score", sig_col, star_col)
                   if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"{path} is not a Pillar-format scoreset (missing {missing}). "
                     "Use --table for a plain score/sample CSV.")
        for row in reader:
            if dataset and (row.get("Dataset") or "") != dataset:
                continue
            rows.append(row)
    if not rows:
        sys.exit(f"no rows for dataset={dataset!r} in {path}")

    detects_splice = (rows[0].get("splice_measure") or "").strip().lower() == "yes"

    def _f(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    # Group by variant ID: one score per variant, memberships unioned.
    by_id: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        if (row.get("Flag") or "").strip() == "*":
            continue
        score = _f(row.get("auth_reported_score"))
        if math.isnan(score):
            continue
        cons = (row.get("simplified_consequence") or "").strip().lower()
        if not detects_splice:
            if cons in ("splice region", "splice_site_variant"):
                continue
            # Upstream requires all four SpliceAI deltas to be present and
            # < 0.2: a missing annotation cannot rule out a splice effect, so
            # the variant is dropped rather than assumed splice-neutral.
            ds = [_f(row.get(c)) for c in ("spliceAI_DS_AG", "spliceAI_DS_AL",
                                           "spliceAI_DS_DG", "spliceAI_DS_DL")]
            if any(math.isnan(v) or v >= 0.2 for v in ds):
                continue
        vid = row.get("ID") or f"row{len(by_id)}"
        sig = (row.get(sig_col) or "").strip()
        stars = row.get(star_col) or ""
        good = _stars_sufficient(stars, min_clinvar_star)
        maf = _f(row.get("gnomad_MAF"))
        is_snv = (len(str(row.get("ref_allele") or "")) == 1
                  and len(str(row.get("alt_allele") or "")) == 1)
        rec = by_id.setdefault(vid, {"score": score, "syn": False, "pop": False,
                                     "path": False, "ben": False})
        rec["syn"] |= cons in ("synonymous", "synonymous_variant")
        if population_type == "all_variants":
            rec["pop"] = True
        elif population_type in ("all_nsSNV", "gnomAD_nsSNV"):
            rec["pop"] |= is_snv and (population_type == "all_nsSNV"
                                      or not math.isnan(maf))
        else:                                    # 'gnomAD' (default)
            rec["pop"] |= not math.isnan(maf)
        rec["path"] |= good and sig in CLINVAR_PATHOGENIC
        rec["ben"] |= good and sig in CLINVAR_BENIGN

    ids, scores, assign = [], [], []
    for vid, rec in by_id.items():
        row = [False, False, False, False]
        if rec["syn"]:
            row[SYNONYMOUS] = True
        else:
            row[POPULATION] = rec["pop"]
            row[PATHOGENIC] = rec["path"]
            row[BENIGN] = rec["ben"]
        if any(row):
            ids.append(vid)
            scores.append(rec["score"])
            assign.append(row)
    if not scores:
        sys.exit("no labelled variants after filtering (check --clinvar-release "
                 "/ --min-clinvar-star)")
    return Scoreset(np, scores, assign, SAMPLE_NAMES, ids)


def load_table(np, path: Path, dataset: Optional[str] = None) -> Scoreset:
    """Load a plain `score,sample[,id]` CSV.

    `sample` is either the fixed index (0 P/LP, 1 B/LB, 2 population,
    3 synonymous) or one of the sample names.

    A variant may legitimately belong to several samples (ClinVar-labelled
    *and* seen in gnomAD). Rows that repeat an `ID` are therefore folded into
    one variant with unioned memberships — so the long-format table written by
    `calibrate prepare` round-trips to exactly the scoreset it came from,
    instead of double-counting those variants as independent observations.
    """
    name_to_idx = {n.lower(): i for i, n in enumerate(SAMPLE_NAMES)}
    name_to_idx["gnomad"] = POPULATION
    name_to_idx["p/lp"] = PATHOGENIC
    name_to_idx["b/lb"] = BENIGN
    by_id: Dict[str, Dict[str, Any]] = {}
    ids, scores, assign = [], [], []
    with _open_maybe_gz(path) as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        score_col = cols.get("score") or cols.get("scores")
        sample_col = cols.get("sample") or cols.get("sample_assignments")
        if not score_col or not sample_col:
            sys.exit(f"{path} needs 'score' and 'sample' columns "
                     f"(found {reader.fieldnames})")
        for n, row in enumerate(reader):
            if dataset and (row.get("Dataset") or row.get("dataset") or "") != dataset:
                continue
            try:
                score = float(row[score_col])
            except (TypeError, ValueError):
                continue
            raw = (row[sample_col] or "").strip()
            if raw.lstrip("-").isdigit():
                idx = int(raw)
            else:
                idx = name_to_idx.get(raw.lower(), -1)
            if not (0 <= idx < 4):
                continue
            vid = row.get("ID") or row.get("id") or f"row{n}"
            rec = by_id.setdefault(vid, {"score": score, "one": [False] * 4})
            rec["one"][idx] = True
    for vid, rec in by_id.items():
        ids.append(vid)
        scores.append(rec["score"])
        assign.append(rec["one"])
    if not scores:
        sys.exit(f"no usable rows in {path}")
    return Scoreset(np, scores, assign, SAMPLE_NAMES, ids)


def load_mapped_with_clinvar(np, mapped: Path, clinvar_tsv: Optional[Path],
                             min_clinvar_star: int = 1,
                             gnomad_tsv: Optional[Path] = None) -> Scoreset:
    """Label a `mavedb map-scoreset` output (chr/pos/ref/alt + score) by
    joining ClinVar's `variant_summary.txt.gz` on GRCh38 coordinates.

    This is the bridge from IGVFagent's own MAVE mapping chain into
    calibration: `mavedb map-scoreset` → this → `calibrate run`.
    """
    clinvar: Dict[Tuple[str, int, str, str], Tuple[str, str]] = {}
    if clinvar_tsv is not None:
        with _open_maybe_gz(clinvar_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if (row.get("Assembly") or "") != "GRCh38":
                    continue
                try:
                    pos = int(row.get("PositionVCF") or 0)
                except ValueError:
                    continue
                key = (str(row.get("Chromosome") or "").replace("chr", ""), pos,
                       (row.get("ReferenceAlleleVCF") or "").upper(),
                       (row.get("AlternateAlleleVCF") or "").upper())
                clinvar[key] = (row.get("ClinicalSignificance") or "",
                                row.get("ReviewStatus") or "")
        logging.info("loaded %d GRCh38 ClinVar records", len(clinvar))

    gnomad: set = set()
    if gnomad_tsv is not None:
        with _open_maybe_gz(gnomad_tsv) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    gnomad.add((f[0].replace("chr", ""), int(f[1]),
                                f[2].upper(), f[3].upper()))
        logging.info("loaded %d gnomAD positions", len(gnomad))

    ids, scores, assign = [], [], []
    with _open_maybe_gz(mapped) as fh:
        sniff = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sniff.splitlines()[0] else ","
        reader = csv.DictReader(fh, delimiter=delim)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        need = [c for c in ("chr", "pos", "ref", "alt") if c not in cols]
        score_col = next((cols[c] for c in ("score", "auth_reported_score",
                                            "mave_score") if c in cols), None)
        if need or not score_col:
            sys.exit(f"{mapped} needs chr/pos/ref/alt + a score column "
                     f"(missing {need or 'score'})")
        cons_col = next((cols[c] for c in ("consequence", "simplified_consequence")
                         if c in cols), None)
        for n, row in enumerate(reader):
            try:
                score = float(row[score_col])
                key = (str(row[cols["chr"]]).replace("chr", ""),
                       int(float(row[cols["pos"]])),
                       (row[cols["ref"]] or "").upper(),
                       (row[cols["alt"]] or "").upper())
            except (TypeError, ValueError):
                continue
            one = [False] * 4
            cons = (row.get(cons_col) or "").strip().lower() if cons_col else ""
            if cons in ("synonymous", "synonymous_variant"):
                one[SYNONYMOUS] = True
            else:
                sig, review = clinvar.get(key, ("", ""))
                good = _stars_sufficient(review, min_clinvar_star)
                one[PATHOGENIC] = good and sig in CLINVAR_PATHOGENIC
                one[BENIGN] = good and sig in CLINVAR_BENIGN
                one[POPULATION] = (key in gnomad) if gnomad else False
            if any(one):
                ids.append(row.get("ID") or row.get("id") or f"row{n}")
                scores.append(score)
                assign.append(one)
    if not scores:
        sys.exit("no variants could be labelled — supply --clinvar-tsv "
                 "(variant_summary.txt.gz) and/or --gnomad-tsv")
    return Scoreset(np, scores, assign, SAMPLE_NAMES, ids)


# ─── Calibration from bootstrap fits ────────────────────────────────────────

def calibrate_from_fits(np, sps, logsumexp, scoreset: Scoreset,
                        fits: List[Dict[str, Any]], benign_method: str,
                        point_values: Sequence[int],
                        use_median_prior: bool = True) -> Dict[str, Any]:
    """Aggregate bootstrap fits into one calibration (prior, C, point ranges)."""
    idx = dict(scoreset.index)
    if idx.get("population") is None:
        raise ValueError("a population (gnomAD) sample is required")
    if idx.get("pathogenic") is None and idx.get("benign") is None \
            and idx.get("synonymous") is None:
        raise ValueError("need a pathogenic or benign/synonymous sample")
    if idx.get("synonymous") is None and benign_method in ("avg", "synonymous"):
        logging.warning("no synonymous sample — benign_method %s -> benign",
                        benign_method)
        benign_method = "benign"
    elif idx.get("benign") is None and benign_method in ("avg", "benign"):
        logging.warning("no benign sample — benign_method %s -> synonymous",
                        benign_method)
        benign_method = "synonymous"

    pop_scores = scoreset.scores[scoreset.assignments[:, idx["population"]]]

    priors, kept = [], []
    for fit in fits:
        params = [tuple(p) for p in fit["component_params"]]
        W = np.asarray(fit["weights"], dtype=float)
        try:
            a = fit_prior(np, sps, logsumexp, pop_scores, params, W, idx,
                          benign_method)
        except (ValueError, FloatingPointError):
            a = float("nan")
        if not math.isnan(a) and 0.0 < a < 1.0:
            priors.append(a)
            kept.append(fit)
    if not kept:
        raise ValueError("no bootstrap fit produced a valid prior")
    priors = np.asarray(priors, dtype=float)
    prior = float(np.nanmedian(priors))
    logging.info("valid priors %d/%d — median prior %.6g",
                 len(kept), len(fits), prior)

    grid = np.linspace(scoreset.scores.min(), scoreset.scores.max(), _SCORE_GRID)
    log_lr = np.full((len(kept), len(grid)), np.nan)
    Cs = []
    per_fit_p: List[Dict[int, Any]] = []
    per_fit_b: List[Dict[int, Any]] = []
    for i, (fit, a) in enumerate(zip(kept, priors)):
        params = [tuple(p) for p in fit["component_params"]]
        W = np.asarray(fit["weights"], dtype=float)
        log_fp, log_fb = class_log_densities(np, sps, logsumexp, grid, params, W,
                                            idx, benign_method, prior)
        if log_fp is None:
            continue
        # A bootstrap only speaks for the score span it actually saw.
        xmin, xmax = fit["xlims"]
        inside = (grid >= xmin) & (grid <= xmax)
        lr = np.full(len(grid), np.nan)
        lr[inside] = log_fp[inside] - log_fb[inside]
        log_lr[i] = lr
        # This bootstrap's own boundaries, at its own prior — the input to the
        # conservative (--no-median-prior) threshold scheme.
        rp, rb, Ci = calculate_score_ranges(np, lr[inside], lr[inside], float(a),
                                           grid[inside], point_values)
        per_fit_p.append({int(k): np.asarray(v, dtype=float).reshape(-1)
                          for k, v in rp.items()})
        per_fit_b.append({int(k): np.asarray(v, dtype=float).reshape(-1)
                          for k, v in rb.items()})
        Cs.append(int(Ci))

    covered = ~np.isnan(log_lr).all(axis=0)
    if not covered.any():
        raise ValueError("no score position is covered by any bootstrap")

    flipped = _detect_flip(np, scoreset, idx, benign_method)
    lo = np.nanpercentile(log_lr[:, covered], 5, axis=0)
    hi = np.nanpercentile(log_lr[:, covered], 95, axis=0)
    scores_covered = grid[covered]
    if use_median_prior:
        ranges_p, ranges_b, C = calculate_score_ranges(
            np, lo, hi, prior, scores_covered, point_values)
        point_ranges = {**ranges_p, **ranges_b}
    else:
        C = tavtigian_constant(np, prior)
        point_ranges = conservative_point_ranges(
            np, per_fit_p, per_fit_b, point_values, scores_covered, flipped)
    enforce_monotonicity(point_ranges, list(point_values), flipped)
    extend_to_xlims(point_ranges, list(point_values), scores_covered, flipped)

    return {
        "prior": prior,
        "prior_iqr": [float(np.nanpercentile(priors, 25)),
                      float(np.nanpercentile(priors, 75))],
        "tavtigian_C": int(C),
        "C_bootstrap_5_95": ([float(np.nanpercentile(Cs, 5)),
                              float(np.nanpercentile(Cs, 95))] if Cs else None),
        "lr_thresholds": {int(p): float(C ** (p / len(point_values)))
                          for p in point_values},
        "point_ranges": {int(k): [[float(a), float(b)] for a, b in v]
                         for k, v in point_ranges.items()},
        "scoreset_flipped": bool(flipped),
        "benign_method": benign_method,
        "n_valid_fits": len(kept),
        "score_range": [float(scores_covered[0]), float(scores_covered[-1])],
        "sample_counts": dict(zip(scoreset.names, scoreset.counts)),
        "_grid": scores_covered,
        "_log_lr_low": lo,
        "_log_lr_high": hi,
        "_fits": kept,
    }


def _detect_flip(np, scoreset: Scoreset, idx, benign_method: str) -> bool:
    """True when higher assay scores mean *more* pathogenic (e.g. depletion
    assays reported as loss-of-function magnitude)."""
    if idx.get("pathogenic") is not None:
        p_mean = scoreset.scores[scoreset.assignments[:, idx["pathogenic"]]].mean()
    else:
        p_mean = scoreset.scores[scoreset.assignments[:, idx["population"]]].mean()
    b, s = idx.get("benign"), idx.get("synonymous")
    if benign_method == "avg" and b is not None and s is not None:
        b_mean = (scoreset.scores[scoreset.assignments[:, b]].mean()
                  + scoreset.scores[scoreset.assignments[:, s]].mean()) / 2
    elif benign_method == "synonymous" and s is not None:
        b_mean = scoreset.scores[scoreset.assignments[:, s]].mean()
    elif b is not None:
        b_mean = scoreset.scores[scoreset.assignments[:, b]].mean()
    elif s is not None:
        b_mean = scoreset.scores[scoreset.assignments[:, s]].mean()
    else:
        b_mean = scoreset.scores[scoreset.assignments[:, idx["population"]]].mean()
    return bool(p_mean > b_mean)


def assign_points(np, scores, point_ranges: Dict[int, List[List[float]]]):
    """Evidence points for each score from a calibration's point ranges."""
    scores = np.asarray(scores, dtype=float)
    pts = np.zeros(len(scores), dtype=int)
    for p, ranges in sorted(point_ranges.items(), key=lambda kv: abs(int(kv[0]))):
        for lo, hi in ranges:
            pts[(scores >= lo) & (scores <= hi)] = int(p)
    return pts


def evidence_label(points: int) -> str:
    """ACMG/AMP criterion + strength for a signed point value."""
    if points == 0:
        return "no evidence"
    code = "PS3" if points > 0 else "BS3"
    return f"{code} {STRENGTH_BY_POINTS.get(min(abs(points), 8), 'supporting')}"


# ─── Figure ─────────────────────────────────────────────────────────────────

def make_figure(np, scoreset: Scoreset, cal: Dict[str, Any], name: str,
                outdir: Path) -> List[Path]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib unavailable — skipping figure")
        return []
    _, sps, logsumexp = _lazy()

    grid = np.asarray(cal["_grid"], dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                             gridspec_kw={"height_ratios": [1.1, 1.0]})
    colours = {"Pathogenic/Likely Pathogenic": "#c1272d",
               "Benign/Likely Benign": "#0b6e4f",
               "population": "#4a4a4a", "Synonymous": "#2b6cb0"}

    ax = axes[0]
    for j, sname in enumerate(scoreset.names):
        vals = scoreset.scores[scoreset.assignments[:, j]]
        ax.hist(vals, bins=60, density=True, alpha=0.35,
                color=colours.get(sname, "#888"),
                label=f"{sname} (n={len(vals)})")
    # Median-bootstrap mixture densities for the P and B classes.
    fits = cal["_fits"]
    idx = scoreset.index
    if fits:
        mid = fits[len(fits) // 2]
        params = [tuple(p) for p in mid["component_params"]]
        W = np.asarray(mid["weights"], dtype=float)
        for key, colour, label in (("pathogenic", "#c1272d", "fitted $f_P$"),
                                   ("benign", "#0b6e4f", "fitted $f_B$")):
            if idx.get(key) is not None:
                dens = np.exp(mixture_logpdf(np, sps, logsumexp, grid, params,
                                             W[idx[key]]))
                ax.plot(grid, dens, color=colour, lw=1.8, label=label)
    ax.set_ylabel("density")
    ax.set_title(f"{name} — assay-score distributions and fitted mixture "
                 f"({cal['n_valid_fits']} bootstraps)")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    ax.plot(grid, cal["_log_lr_low"], color="#c1272d", lw=1.4,
            label="conservative log LR$^+$ (5th pct)")
    ax.plot(grid, cal["_log_lr_high"], color="#0b6e4f", lw=1.4,
            label="conservative log LR$^+$ (95th pct)")
    ax.axhline(0.0, color="#999", lw=0.8, ls=":")
    C = cal["tavtigian_C"]
    for p, style in ((1, ":"), (2, "-."), (4, "--"), (8, "-")):
        tau = math.log(C ** (p / 8))
        ax.axhline(tau, color="#c1272d", lw=0.7, ls=style, alpha=0.7)
        ax.axhline(-tau, color="#0b6e4f", lw=0.7, ls=style, alpha=0.7)
        ax.text(grid[-1], tau, f" {STRENGTH_BY_POINTS[p]}", fontsize=6,
                color="#c1272d", va="center")
    for p, ranges in cal["point_ranges"].items():
        for lo, hi in ranges:
            ax.axvspan(lo, hi, color="#c1272d" if int(p) > 0 else "#0b6e4f",
                       alpha=0.06 + 0.02 * min(abs(int(p)), 8))
    ax.set_xlabel("assay score")
    ax.set_ylabel("log LR$^+$")
    ax.set_title(f"Evidence thresholds — prior {cal['prior']:.4g}, "
                 f"Tavtigian C = {C}")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "svg"):
        p = outdir / f"{name}_calibration.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_thresholds(args) -> int:
    np, _, _ = _lazy()
    pts = [int(p) for p in args.point_values.split(",")] if args.point_values \
        else DEFAULT_POINT_VALUES
    lr_p, lr_b, C = thresholds_from_prior(np, args.prior, pts,
                                          original=args.original,
                                          strict=args.strict)
    print(f"prior (P(pathogenic) in the reference population) : {args.prior:g}")
    print(f"Tavtigian constant C = O_PVSt                      : {C}")
    print()
    print(f"{'points':>6}  {'strength':<12} {'LR+ (PS3, pathogenic)':>22}"
          f" {'LR+ (BS3, benign)':>19}")
    for p, a, b in zip(pts, lr_p, lr_b):
        print(f"{p:>6}  {STRENGTH_BY_POINTS.get(p, ''):<12} {a:>22.4g} {b:>19.4g}")
    if args.out:
        Path(args.out).write_text(json.dumps(
            {"prior": args.prior, "tavtigian_C": int(C),
             "point_values": pts,
             "lr_plus_pathogenic": [float(v) for v in lr_p],
             "lr_plus_benign": [float(v) for v in lr_b]}, indent=2))
        print(f"\nWrote {args.out}")
    return 0


def cmd_prepare(args) -> int:
    np, _, _ = _lazy()
    setup_logging("prepare")
    if args.pillar:
        ss = load_pillar(np, Path(args.pillar), dataset=args.name,
                         clinvar_release=args.clinvar_release,
                         min_clinvar_star=args.min_clinvar_star,
                         population_type=args.population_type)
    elif args.table:
        ss = load_table(np, Path(args.table), dataset=args.name)
    elif args.mapped:
        ss = load_mapped_with_clinvar(
            np, Path(args.mapped),
            Path(args.clinvar_tsv) if args.clinvar_tsv else None,
            min_clinvar_star=args.min_clinvar_star,
            gnomad_tsv=Path(args.gnomad_tsv) if args.gnomad_tsv else None)
    else:
        sys.exit("pick one of --pillar / --table / --mapped")

    label = args.name or "scoreset"
    outdir = DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}"
    outdir.mkdir(parents=True, exist_ok=True)
    out = Path(args.out) if args.out else outdir / f"{label}_calibration_input.csv"
    name_by_col = {j: n for j, n in enumerate(ss.names)}
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ID", "score", "sample", "sample_name"])
        for i, (vid, score) in enumerate(zip(ss.ids, ss.scores)):
            for j in np.where(ss.assignments[i])[0]:
                w.writerow([vid, f"{score:.10g}", SAMPLE_NAMES.index(name_by_col[j]),
                            name_by_col[j]])
    summary = {"dataset": label, "n_variants": int(len(ss.scores)),
               "samples": dict(zip(ss.names, [int(c) for c in ss.counts])),
               "score_min": float(ss.scores.min()),
               "score_max": float(ss.scores.max()),
               "input_table": str(out)}
    (outdir / f"{label}_prepare_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Prepared {summary['n_variants']} labelled variants -> {out}")
    for n, c in summary["samples"].items():
        print(f"  {n:<32} {c}")
    if summary["samples"].get("population", 0) == 0:
        print("  ! no population sample — calibration needs gnomAD membership")
    return 0


def cmd_run(args) -> int:
    np, sps, logsumexp = _lazy()
    setup_logging("run")
    label = args.name or "scoreset"
    if args.pillar:
        ss = load_pillar(np, Path(args.pillar), dataset=args.name,
                         clinvar_release=args.clinvar_release,
                         min_clinvar_star=args.min_clinvar_star,
                         population_type=args.population_type)
    elif args.table:
        ss = load_table(np, Path(args.table), dataset=None)
    else:
        sys.exit("pick one of --pillar / --table")
    logging.info("%s", ss)
    if len(ss.names) < 3:
        sys.exit(f"insufficient samples: {len(ss.names)} < 3 "
                 f"({dict(zip(ss.names, ss.counts))})")

    components = sorted({int(c) for c in args.components})
    if args.output_dir:
        outdir = Path(args.output_dir)
    elif args.resume:
        # Resume needs a stable directory: reuse the newest previous run for
        # this name (resolved *before* a new timestamped dir would shadow it).
        prior_runs = sorted(p for p in DOC_DIR.glob(f"*_{label}") if p.is_dir())
        outdir = (prior_runs[-1] if prior_runs
                  else DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}")
        if prior_runs:
            logging.info("resuming into %s", outdir)
    else:
        outdir = DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}"
    outdir.mkdir(parents=True, exist_ok=True)
    ledger = outdir / f"{label}_bootstrap_fits.jsonl"

    one_hot = ss.one_hot(np, seed=args.seed)
    logging.info("disjoint sample counts: %s",
                 dict(zip(ss.names, one_hot.sum(axis=0).tolist())))
    jobs = args.jobs if args.jobs > 0 else (os.cpu_count() or 1)
    logging.info("fitting %d bootstraps x %d fits, components %s, %d workers",
                 args.n_bootstraps, args.fits_per_bootstrap, components, jobs)

    records = run_bootstraps(np, ss.scores, one_hot, components,
                             args.n_bootstraps, args.fits_per_bootstrap, jobs,
                             ledger if not args.no_ledger else None,
                             max_em_iters=args.max_em_iters)
    errors = [r for r in records if "error" in r]
    if errors:
        logging.warning("%d bootstraps errored (e.g. %s)", len(errors),
                        errors[0].get("error"))

    # Model selection between exactly two component counts.
    selection = None
    chosen = components
    if len(components) == 2 and not args.no_auto_select:
        k_low, k_high = components
        pairs = [(r[f"{k_low}c"]["val_ll"], r[f"{k_high}c"]["val_ll"])
                 for r in records
                 if f"{k_low}c" in r and f"{k_high}c" in r]
        if pairs:
            selection = paired_bootstrap_test(
                np, [p[0] for p in pairs], [p[1] for p in pairs], k_low, k_high,
                alpha=args.selection_alpha)
            key = "conservative_k" if not args.no_conservative else "selected_k"
            chosen = [selection[key]]
            (outdir / f"{label}_model_selection.json").write_text(
                json.dumps(selection, indent=2))
            logging.info("model selection (%s): %dc (p=%.4g, 5th pct %+.3g)",
                         key, chosen[0], selection["p_value"],
                         selection["fifth_percentile"])

    written = []
    for k in chosen:
        fits = [r[f"{k}c"] for r in records if f"{k}c" in r]
        if not fits:
            logging.warning("no valid %dc fits", k)
            continue
        cal = calibrate_from_fits(np, sps, logsumexp, ss, fits,
                                  args.benign_method, DEFAULT_POINT_VALUES,
                                  use_median_prior=not args.no_median_prior)
        figs = make_figure(np, ss, cal, f"{label}_{k}c", outdir)
        public = {"dataset": label, "component": f"{k}c",
                  **{key: v for key, v in cal.items() if not key.startswith("_")},
                  "evidence_summary": {
                      str(p): evidence_label(int(p))
                      for p, r in sorted(cal["point_ranges"].items(),
                                         key=lambda kv: -int(kv[0])) if r},
                  "config": {"n_bootstraps": args.n_bootstraps,
                             "fits_per_bootstrap": args.fits_per_bootstrap,
                             "components": k, "benign_method": args.benign_method,
                             "use_median_prior": not args.no_median_prior,
                             "clinvar_release": args.clinvar_release,
                             "min_clinvar_star": args.min_clinvar_star,
                             "seed": args.seed},
                  "figures": [str(p) for p in figs],
                  "method": "exCALIBR / Zeiberg et al. bioRxiv 2025.04.29.651326"}
        path = outdir / f"{label}_{k}c_calibration.json"
        path.write_text(json.dumps(public, indent=2))
        written.append(path)
        print(f"\n{label} {k}c — prior {cal['prior']:.4g}, C = {cal['tavtigian_C']}, "
              f"{cal['n_valid_fits']} valid bootstraps"
              f"{' (flipped scale)' if cal['scoreset_flipped'] else ''}")
        for p in sorted(cal["point_ranges"], key=lambda v: -int(v)):
            for lo, hi in cal["point_ranges"][p]:
                print(f"  {evidence_label(int(p)):<22} "
                      f"score in [{lo:.4g}, {hi:.4g}]  ({int(p):+d} points)")
        print(f"  calibration -> {path}")
        for f in figs:
            print(f"  figure      -> {f}")
    if not written:
        sys.exit("no calibration produced — see the log for failed fits")
    return 0


def cmd_assign(args) -> int:
    np, _, _ = _lazy()
    cal = json.loads(Path(args.calibration).read_text())
    ranges = {int(k): v for k, v in cal["point_ranges"].items()}
    rows, scores = [], []
    with _open_maybe_gz(Path(args.scores)) as fh:
        sniff = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sniff else ","
        reader = csv.DictReader(fh, delimiter=delim)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        score_col = next((cols[c] for c in ("score", "auth_reported_score",
                                            "mave_score") if c in cols), None)
        if not score_col:
            sys.exit(f"{args.scores} has no score column ({reader.fieldnames})")
        for row in reader:
            try:
                scores.append(float(row[score_col]))
            except (TypeError, ValueError):
                continue
            rows.append(row)
    pts = assign_points(np, scores, ranges)
    out = Path(args.out) if args.out else Path(args.scores).with_suffix(".calibrated.tsv")
    extra = [c for c in (rows[0].keys() if rows else []) if c != score_col]
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["score", "evidence_points", "acmg_evidence"] + extra)
        for row, s, p in zip(rows, scores, pts):
            w.writerow([f"{s:.10g}", int(p), evidence_label(int(p))]
                       + [row.get(c, "") for c in extra])
    from collections import Counter
    counts = Counter(evidence_label(int(p)) for p in pts)
    print(f"Assigned ACMG evidence to {len(pts)} scores using "
          f"{cal.get('dataset')} {cal.get('component')} "
          f"(prior {cal.get('prior'):.4g}, C={cal.get('tavtigian_C')})")
    for label, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {label:<24} {n}")
    print(f"-> {out}")
    return 0


def cmd_selftest(args) -> int:
    """Validate the port: published Tavtigian constants, then recovery of a
    known calibration on a synthetic scoreset."""
    np, sps, logsumexp = _lazy()
    ok = True

    # 1. Tavtigian's Bayesian framework — published LR+ per evidence strength
    #    at a 10 % prior: supporting 2.08, moderate 4.33, strong 18.7, VS 350.
    C = tavtigian_constant(np, 0.1)
    got = [round(float(C ** (p / 8)), 2) for p in (1, 2, 4, 8)]
    want = [2.08, 4.32, 18.66, 348.0]
    print(f"[1] Tavtigian C at prior 0.10 : C={C}  LR+ (sup/mod/str/VS) = {got}")
    close = all(abs(g - w) <= max(0.02, 0.01 * w) for g, w in zip(got, want))
    print(f"    expected ~{want} (Tavtigian 2018 published 2.08 / 4.33 / 18.7 / 350)"
          f" -> {'PASS' if close else 'FAIL'}")
    ok &= close

    # C must increase as the prior falls (weaker prior needs stronger evidence).
    priors = [0.1, 0.05, 0.02]
    Cs = [tavtigian_constant(np, p) for p in priors]
    mono = all(b > a for a, b in zip(Cs, Cs[1:]))
    print(f"[2] C monotone in prior       : {list(zip(priors, Cs))} -> "
          f"{'PASS' if mono else 'FAIL'}")
    ok &= mono

    # 3. Synthetic scoreset with a known 20 % pathogenic fraction in the
    #    population: pathogenic scores low, benign high, population a mixture.
    rng = np.random.default_rng(0)
    true_prior = 0.2
    n_path, n_ben, n_pop, n_syn = 120, 150, 600, 100
    path = rng.normal(-2.0, 0.7, n_path)
    ben = rng.normal(0.0, 0.7, n_ben)
    n_pop_path = int(round(true_prior * n_pop))
    pop = np.concatenate([rng.normal(-2.0, 0.7, n_pop_path),
                          rng.normal(0.0, 0.7, n_pop - n_pop_path)])
    syn = rng.normal(0.0, 0.7, n_syn)
    scores = np.concatenate([path, ben, pop, syn])
    assign = np.zeros((len(scores), 4), dtype=bool)
    o = 0
    for col, n in ((PATHOGENIC, n_path), (BENIGN, n_ben), (POPULATION, n_pop),
                   (SYNONYMOUS, n_syn)):
        assign[o:o + n, col] = True
        o += n
    ss = Scoreset(np, scores, assign, SAMPLE_NAMES)
    one_hot = ss.one_hot(np, seed=0)
    records = run_bootstraps(np, ss.scores, one_hot, [2],
                             args.n_bootstraps, args.fits_per_bootstrap,
                             args.jobs if args.jobs > 0 else (os.cpu_count() or 1),
                             None, heartbeat=5)
    fits = [r["2c"] for r in records if "2c" in r]
    print(f"[3] synthetic 2c fits         : {len(fits)}/{args.n_bootstraps} bootstraps")
    if not fits:
        print("    -> FAIL (no fits)")
        return 1
    cal = calibrate_from_fits(np, sps, logsumexp, ss, fits, "avg",
                              DEFAULT_POINT_VALUES)
    err = abs(cal["prior"] - true_prior)
    print(f"[4] recovered prior           : {cal['prior']:.4f} "
          f"(truth {true_prior}) -> {'PASS' if err < 0.07 else 'FAIL'}")
    ok &= err < 0.07

    # Pathogenic evidence must sit at low scores, benign at high, and the
    # strongest available strengths must be reached on both sides.
    pos = [p for p, r in cal["point_ranges"].items() if int(p) > 0 and r]
    neg = [p for p, r in cal["point_ranges"].items() if int(p) < 0 and r]
    sides_ok = bool(pos) and bool(neg)
    if sides_ok:
        p_hi = max(cal["point_ranges"][max(pos, key=lambda v: int(v))],
                   key=lambda r: r[1])
        b_hi = max(cal["point_ranges"][min(neg, key=lambda v: int(v))],
                   key=lambda r: r[1])
        sides_ok = p_hi[0] < b_hi[1] and p_hi[1] <= b_hi[1]
        print(f"[5] evidence orientation      : PS3 up to {p_hi[1]:.3g}, "
              f"BS3 from {b_hi[0]:.3g} -> {'PASS' if sides_ok else 'FAIL'}")
    else:
        print("[5] evidence orientation      : FAIL (one side has no evidence)")
    ok &= sides_ok

    scored = assign_points(np, [-3.0, -2.0, 0.0, 1.0], cal["point_ranges"])
    print(f"[6] assign_points sanity      : {[evidence_label(int(p)) for p in scored]}")
    print(f"\nselftest: {'ALL PASS' if ok else 'FAILURES PRESENT'}")
    return 0 if ok else 1


def cmd_write_playbook(args) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(__doc__ or "")
    print(f"Wrote {PLAYBOOK_PATH}")
    return 0


# ─── CLI ────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent calibrate",
        description="Calibrate functional-assay scores to ACMG/AMP evidence "
                    "strengths (clean-room port of exCALIBR / Zeiberg et al.).")
    sub = p.add_subparsers(dest="cmd")

    sp = sub.add_parser("thresholds",
                        help="Tavtigian constant C + LR+ evidence thresholds "
                             "for a prior (no fitting).")
    sp.add_argument("--prior", type=float, required=True,
                    help="P(pathogenic) in the reference population, e.g. 0.1")
    sp.add_argument("--point-values", default=None,
                    help="comma list (default 1..8)")
    sp.add_argument("--original", action="store_true",
                    help="use the original ACMG/AMP combining rules")
    sp.add_argument("--strict", action="store_true",
                    help="also penalise posteriors above the LP/LB bands")
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_thresholds)

    sp = sub.add_parser("prepare",
                        help="Build a score/sample calibration input table.")
    src = sp.add_mutually_exclusive_group(required=True)
    src.add_argument("--pillar", help="IGVF / Pillar-format scoreset CSV")
    src.add_argument("--table", help="plain CSV with score + sample columns")
    src.add_argument("--mapped", help="mavedb map-scoreset output (chr/pos/ref/alt)")
    sp.add_argument("--name", default=None, help="dataset name (also filters --pillar)")
    sp.add_argument("--clinvar-release", default="2025", choices=["2025", "2018"])
    sp.add_argument("--min-clinvar-star", type=int, default=1)
    sp.add_argument("--population-type", default="gnomAD",
                    choices=["gnomAD", "gnomAD_nsSNV", "all_variants", "all_nsSNV"])
    sp.add_argument("--clinvar-tsv", default=None,
                    help="ClinVar variant_summary.txt.gz (with --mapped)")
    sp.add_argument("--gnomad-tsv", default=None,
                    help="chr/pos/ref/alt table of gnomAD variants (with --mapped)")
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_prepare)

    sp = sub.add_parser("run", help="Full bootstrap calibration pipeline.")
    src = sp.add_mutually_exclusive_group(required=True)
    src.add_argument("--pillar", help="IGVF / Pillar-format scoreset CSV")
    src.add_argument("--table", help="score/sample CSV (from `calibrate prepare`)")
    sp.add_argument("--name", default=None, help="dataset name for outputs")
    sp.add_argument("--components", nargs="+", type=int, default=[2, 3],
                    choices=[2, 3, 4],
                    help="mixture component counts to fit (default 2 3)")
    sp.add_argument("--n-bootstraps", type=int, default=1000)
    sp.add_argument("--fits-per-bootstrap", type=int, default=100)
    sp.add_argument("--max-em-iters", type=int, default=10000)
    sp.add_argument("--benign-method", default="avg",
                    choices=["avg", "benign", "synonymous"])
    sp.add_argument("--no-median-prior", action="store_true",
                    help="conservative percentile thresholds instead of the "
                         "median-prior envelope")
    sp.add_argument("--no-auto-select", action="store_true")
    sp.add_argument("--no-conservative", action="store_true",
                    help="select on the Wilcoxon p-value instead of the 5th "
                         "percentile rule")
    sp.add_argument("--selection-alpha", type=float, default=0.05)
    sp.add_argument("--clinvar-release", default="2025", choices=["2025", "2018"])
    sp.add_argument("--min-clinvar-star", type=int, default=1)
    sp.add_argument("--population-type", default="gnomAD",
                    choices=["gnomAD", "gnomAD_nsSNV", "all_variants", "all_nsSNV"])
    sp.add_argument("--jobs", type=int, default=-1, help="-1 = all CPUs")
    sp.add_argument("--seed", type=int, default=0)
    sp.add_argument("--resume", action="store_true",
                    help="continue the newest run for this name from its ledger")
    sp.add_argument("--no-ledger", action="store_true",
                    help="do not write the resumable bootstrap ledger")
    sp.add_argument("--output-dir", default=None)
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("assign",
                        help="Apply a calibration to a score table -> ACMG "
                             "evidence per variant.")
    sp.add_argument("--calibration", required=True, help="*_calibration.json")
    sp.add_argument("--scores", required=True, help="CSV/TSV with a score column")
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_assign)

    sp = sub.add_parser("selftest",
                        help="Validate against published Tavtigian constants "
                             "and a synthetic scoreset with a known prior.")
    sp.add_argument("--n-bootstraps", type=int, default=20)
    sp.add_argument("--fits-per-bootstrap", type=int, default=8)
    sp.add_argument("--jobs", type=int, default=-1)
    sp.set_defaults(func=cmd_selftest)

    sp = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/ASSAY_CALIBRATION_SKILL.md")
    sp.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
