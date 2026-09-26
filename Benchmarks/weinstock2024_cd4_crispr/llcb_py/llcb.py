"""A from-source Python port of the causal-network-inference core of
Weinstock et al. 2024's LLCB method (github.com/weinstockj/LLCB, Julia +
Turing.jl). Ports `graph.jl`'s `interventionGraph` / `estimate_total_effects`
/ `get_cyclic_matrices`, and replaces `turing_models.jl`'s
`joint_cyclic_model` + Pathfinder-based approximate inference with a
per-node closed-form/MAP solve (justified below).

Why per-block instead of porting Pathfinder/NUTS verbatim
-------------------------------------------------------------
`get_cyclic_matrices` builds one global system `T*beta = t`, but T is
block-diagonal: each of the `nv` blocks depends only on the `nv-1` beta
entries describing what feeds into one target gene. Of `joint_cyclic_model`'s
two soft regularization terms on the full adjacency matrix W, the per-column
sparsity penalty (`sum(abs(W_col)) ~ Normal(0, nv*0.1)`) is separable per
block too, since column v of W *is* beta_v -- `fit_block_map` carries it.
Only the spectral-radius stability penalty needs the full W jointly, and is
dropped (see its docstring). This turns a ~7,000-parameter coupled MCMC/
Pathfinder problem into `nv` independent (nv-1)-dimensional MAP solves,
each cheap enough to do in closed form or a short L-BFGS run.

A first version of this port used a single free noise variance per block,
re-estimated from the residual by an evidence-style fixed point (as in
`sklearn.linear_model.BayesianRidge`). That degenerates here: each block's
`T_v beta_v = t_v` is exactly square (nv-1 equations, nv-1 unknowns), so an
unconstrained per-block noise estimate always collapses toward zero (a
square system can always be fit exactly), which silently switches off the
`beta ~ Normal(0, 1)` shrinkage prior entirely. `estimate_total_effects` and
`fit_block_bayesian_ridge` instead use each total-effect estimate's own
known OLS standard error as fixed (not re-estimated) heteroskedastic
measurement noise -- a well-posed weighted ridge with no free noise
parameter to collapse.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm


def build_intervention_matrix(data: pd.DataFrame, targets: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Port of `graph.jl`'s `interventionGraph` constructor: returns
    (x, interventions), each samples x nv."""
    nv = len(targets)
    x = data[targets].to_numpy(dtype=float)
    interventions = np.zeros((len(data), nv), dtype=np.int8)
    intervention_col = data["intervention"].to_numpy()
    for j, gene in enumerate(targets):
        rows = intervention_col == gene
        interventions[rows, j] = 1
        assert np.allclose(x[rows, j], 0.0), f"{gene}: intervened samples must have expression zeroed"
    return x, interventions


def estimate_total_effects(x: np.ndarray, interventions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Port of `graph.jl`'s `estimate_total_effects` (log_normalize=False,
    center=True, scale=True, matching the paper's own quickstart call).

    Also returns each total-effect estimate's OLS standard error. `graph.jl`
    discards this (it only keeps the point estimate), but `t_v` in the
    downstream block system is a vector of *estimated* total effects, each
    with its own known sampling uncertainty -- treating that as the
    likelihood noise (see `fit_block_bayesian_ridge`) is the natural way to
    turn `beta ~ Normal(0, 1); t ~ Normal(T*beta, sigma_x^2)` into a
    well-posed problem when T is square (see module docstring)."""
    n, nv = x.shape
    control_samples = interventions.sum(axis=1) == 0

    ctrl_mean = x[control_samples].mean(axis=0)
    xc = x - ctrl_mean
    ctrl_std = xc[control_samples].std(axis=0, ddof=1)
    xs = xc / ctrl_std

    experiment_obs = [interventions[:, i].astype(bool) | control_samples for i in range(nv)]

    total_effects = np.zeros((nv, nv))
    total_effects_se = np.full((nv, nv), np.nan)
    for i in range(nv):
        obs = experiment_obs[i]
        xi = xs[obs, i]
        n_obs = obs.sum()
        design = np.column_stack([np.ones(n_obs), xi])
        xtx_inv = np.linalg.inv(design.T @ design)
        for j in range(nv):
            if i == j:
                continue
            yj = xs[obs, j]
            beta_hat, *_ = np.linalg.lstsq(design, yj, rcond=None)
            resid = yj - design @ beta_hat
            dof = max(n_obs - 2, 1)
            resid_var = float(resid @ resid) / dof
            se = np.sqrt(resid_var * xtx_inv[1, 1])
            total_effects[i, j] = beta_hat[1]
            total_effects_se[i, j] = se
    return total_effects, total_effects_se


def get_cyclic_blocks(total_effects: np.ndarray, total_effects_se: np.ndarray) -> list[dict]:
    """Port of `graph.jl`'s `get_cyclic_matrices`, kept as independent
    per-node blocks rather than assembled into one sparse global system
    (see module docstring)."""
    nv = total_effects.shape[0]
    blocks = []
    for observed in range(nv):
        idx = [i for i in range(nv) if i != observed]
        Tv = total_effects[np.ix_(idx, idx)].copy()
        np.fill_diagonal(Tv, 1.0)
        tv = total_effects[idx, observed].copy()
        tv_se = total_effects_se[idx, observed].copy()
        blocks.append({"observed": observed, "parents": idx, "T": Tv, "t": tv, "t_se": tv_se})
    return blocks


def fit_block_bayesian_ridge(
    T: np.ndarray, t: np.ndarray, t_se: np.ndarray, prior_precision: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form posterior mean/cov for beta_v ~ N(0, 1/prior_precision I),
    t_v ~ N(T_v beta_v, diag(t_se^2)) -- i.e. a Bayesian *weighted* ridge
    where each entry of t_v is weighted by its own known OLS uncertainty
    from `estimate_total_effects`, rather than a single free noise scalar
    re-estimated from the fit (which degenerates to zero on this square,
    exactly-determined system -- see module docstring)."""
    w = 1.0 / np.clip(t_se, 1e-6, None) ** 2
    I = np.eye(T.shape[1])
    precision = (T.T * w) @ T + prior_precision * I
    cov = np.linalg.inv(precision)
    mean = cov @ ((T.T * w) @ t)
    return mean, cov


def fit_block_map(
    T: np.ndarray,
    t: np.ndarray,
    t_se: np.ndarray,
    nv: int,
    prior_precision: float = 1.0,
    col_penalty_scale: float | None = None,
    eps: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """MAP estimate of the same block system, additionally carrying
    `joint_cyclic_model`'s per-column sparsity penalty
    (`sum(loglikelihood(Normal(0, nv*0.1), sum(abs(W_col))))` in
    `turing_models.jl`) -- separable per block, since column v of W *is*
    beta_v. The spectral-radius stability penalty is not block-separable
    (it depends on the full adjacency matrix) and is not included; see the
    module docstring for why this is a reasonable simplification.

    Warm-started from, and falls back to, the closed-form ridge above.
    """
    if col_penalty_scale is None:
        col_penalty_scale = nv * 0.1
    w = 1.0 / np.clip(t_se, 1e-6, None) ** 2
    beta0, _ = fit_block_bayesian_ridge(T, t, t_se, prior_precision)

    def smooth_abs(b):
        return np.sqrt(b**2 + eps)

    def nll_and_grad(beta):
        r = t - T @ beta
        wr = w * r
        nll = 0.5 * np.dot(r, wr)
        grad = -T.T @ wr

        nll += 0.5 * prior_precision * np.dot(beta, beta)
        grad += prior_precision * beta

        s = np.sum(smooth_abs(beta))
        nll += 0.5 * (s / col_penalty_scale) ** 2
        grad += (s / col_penalty_scale**2) * (beta / smooth_abs(beta))

        return nll, grad

    result = minimize(nll_and_grad, beta0, jac=True, method="L-BFGS-B")
    beta = result.x

    # Laplace covariance from the (Gauss-Newton) Hessian at the MAP.
    s = np.sum(smooth_abs(beta))
    d2s = eps / smooth_abs(beta) ** 3  # d/dbeta_i of beta_i/smooth_abs(beta_i)
    hessian = (
        (T.T * w) @ T
        + prior_precision * np.eye(T.shape[1])
        + np.outer(beta / smooth_abs(beta), beta / smooth_abs(beta)) / col_penalty_scale**2
        + np.diag((s / col_penalty_scale**2) * d2s)
    )
    cov = np.linalg.inv(hessian)
    return beta, cov


def fit_llcb(data: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    """End-to-end port: data (donor, intervention, <targets>) -> long-format
    edge table (row=regulator, col=target, estimate, std, ci_low, ci_high),
    matching the schema of `infer_cyclic_edges.jl`'s `parse_cyclic_chain`."""
    nv = len(targets)
    x, interventions = build_intervention_matrix(data, targets)
    total_effects, total_effects_se = estimate_total_effects(x, interventions)
    blocks = get_cyclic_blocks(total_effects, total_effects_se)

    rows = []
    for block in blocks:
        observed = block["observed"]
        mean, cov = fit_block_map(block["T"], block["t"], block["t_se"], nv)
        sd = np.sqrt(np.diag(cov))
        p_pos = 1 - norm.cdf(0.0, loc=mean, scale=sd)  # P(beta_k > 0) under the Laplace posterior
        lsfr = np.minimum(p_pos, 1 - p_pos)  # local sign false rate
        for k, parent in enumerate(block["parents"]):
            rows.append(
                {
                    "row": targets[parent],
                    "col": targets[observed],
                    "estimate": mean[k],
                    "std": sd[k],
                    "2.5%": mean[k] - 1.96 * sd[k],
                    "97.5%": mean[k] + 1.96 * sd[k],
                    "lsfr": lsfr[k],
                }
            )
    return pd.DataFrame(rows)
