# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/hbliu/Open4Gene @ 6e7f36aa80e6
# (R/Open4Gene.R) for liu2025_kidney_multiome. Unreviewed; provenance in Scripts/ported/registry.json.
"""Open4Gene peak-to-gene hurdle test, numerically faithful to the R package.

Port of hbliu/Open4Gene R/Open4Gene.R (Open4Gene + AssociationTest, Liu et
al., Science 2025). For each peak-gene pair it fits

    hurdle(RNA ~ ATAC + covariates, link = "logit", dist = "negbin")

exactly as pscl 1.5.9 does, so estimates agree with the authors' R output to
its printed precision:

  * starting values: glm.fit (IRLS) Poisson on all cells for the count part,
    binomial/logit on I(RNA > 0) for the zero part, theta = 1;
  * each part maximised separately with R's BFGS (optim's C routine vmmin,
    reltol = .Machine$double.eps^(1/1.6), maxit = 10000) on pscl's analytic
    log-likelihoods and gradients;
  * standard errors from optim's finite-difference Hessian (optimhess,
    ndeps = 1e-3), AIC/BIC from the summed log-likelihood.

Open4Gene details kept: covariates are coerced with as.integer() (truncation)
before fitting; a pair is tested only when some cells have RNA == 0, at least
MinNum.Cells have RNA > 0 and at least MinNum.Cells have ATAC > 0; Spearman
uses cor.test(exact = FALSE); betas/SEs/z/AIC/rho are rounded to 6 decimals
and p-values to 6 significant digits.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, gammaln

EPS = np.finfo(float).eps
RELTOL = EPS ** (1 / 1.6)
COLS = ["Peak", "Gene", "Celltype", "TotalCellNum", "ExpRessCellNum", "OpenCellNum",
        "hurdle.Res.zero.beta", "hurdle.Res.zero.se", "hurdle.Res.zero.z", "hurdle.Res.zero.p",
        "hurdle.Res.count.beta", "hurdle.Res.count.se", "hurdle.Res.count.z", "hurdle.Res.count.p",
        "hurdle.AIC", "hurdle.BIC", "spearman.rho", "spearman.p"]


# ---------------------------------------------------------------- R numerics
THRESH, MTHRESH, INVEPS = 30.0, -30.0, 1 / EPS


def logit_linkinv(eta):
    """C_logit_linkinv: exp clamped at |eta| > 30, then tmp / (1 + tmp)."""
    with np.errstate(over="ignore"):
        tmp = np.where(eta < MTHRESH, EPS, np.where(eta > THRESH, INVEPS, np.exp(eta)))
    return tmp / (1 + tmp)


def logit_mu_eta(eta):
    """C_logit_mu_eta: DBL_EPSILON outside |eta| <= 30."""
    with np.errstate(over="ignore", invalid="ignore"):
        opexp = 1 + np.exp(eta)
        v = np.exp(eta) / (opexp * opexp)
    return np.where((eta > THRESH) | (eta < MTHRESH), EPS, v)


def log_dnbinom0(theta, mu):
    """dnbinom(0, size = theta, mu = mu, log = TRUE) as R's dnbinom_mu computes it."""
    return theta * np.where(theta < mu, np.log(theta / (theta + mu)), np.log1p(-mu / (theta + mu)))


def glm_fit(X, y, family, maxit=25, epsilon=1e-8):
    """R's glm.fit IRLS (canonical links), returns coefficients."""
    if family == "poisson":
        mu = y + 0.1
        eta = np.log(mu)
        linkinv = mu_eta = lambda e: np.maximum(np.exp(e), EPS)   # poisson()$linkinv / mu.eta
        var = lambda m: m

        def dev(m):
            with np.errstate(divide="ignore", invalid="ignore"):
                t = np.where(y > 0, y * np.log(y / m), 0.0)
            return float(np.sum(2 * (t - (y - m))))
    else:  # binomial, logit
        mu = (y + 0.5) / 2.0
        eta = np.log(mu / (1 - mu))
        linkinv, mu_eta = logit_linkinv, logit_mu_eta
        var = lambda m: m * (1 - m)

        def dev(m):
            with np.errstate(divide="ignore", invalid="ignore"):
                t = np.where(y > 0, y * np.log(y / m), 0.0) + \
                    np.where(y < 1, (1 - y) * np.log((1 - y) / (1 - m)), 0.0)
            return float(np.sum(2 * t))
    mu = linkinv(eta)
    devold = dev(mu)
    coef = None
    for _ in range(maxit):
        me = mu_eta(eta)
        z = eta + (y - mu) / me
        w = np.sqrt(me ** 2 / var(mu))
        good = me != 0
        coef_new, *_ = np.linalg.lstsq(X[good] * w[good, None], z[good] * w[good], rcond=None)
        eta_new = X @ coef_new
        mu_new = linkinv(eta_new)
        dv = dev(mu_new)
        # step-halving when the deviance is not finite (R: "inner loop 1")
        k = 0
        while not np.isfinite(dv) and coef is not None and k < 25:
            coef_new = (coef_new + coef) / 2
            eta_new = X @ coef_new
            mu_new = linkinv(eta_new)
            dv = dev(mu_new)
            k += 1
        coef, eta, mu = coef_new, eta_new, mu_new
        if abs(dv - devold) / (abs(dv) + 0.1) < epsilon:
            break
        devold = dv
    return coef


def vmmin(fn, gr, b, maxit=10000, reltol=RELTOL, abstol=-np.inf):
    """R's variable-metric minimiser (src/appl/optim.c vmmin), minimising fn."""
    stepredn, acctol, reltest = 0.2, 0.0001, 10.0
    b = np.array(b, dtype=float)
    n = b.size
    f = fn(b)
    if not np.isfinite(f):
        raise FloatingPointError("initial value in 'vmmin' is not finite")
    Fmin = f
    g = gr(b)
    gradcount, it = 1, 1
    ilast = gradcount
    B = np.eye(n)
    while True:
        if ilast == gradcount:
            B = np.eye(n)
        X = b.copy()
        c = g.copy()
        t = -(B @ g)
        gradproj = float(t @ g)
        count = 0
        if gradproj < 0.0:
            steplength = 1.0
            accpoint = False
            while True:
                b = X + steplength * t
                count = int(np.sum(reltest + X == reltest + b))
                if count < n:
                    f = fn(b)
                    accpoint = bool(np.isfinite(f) and f <= Fmin + gradproj * steplength * acctol)
                    if not accpoint:
                        steplength *= stepredn
                if count == n or accpoint:
                    break
            enough = (f > abstol) and abs(f - Fmin) > reltol * (abs(Fmin) + reltol)
            if not enough:
                count = n
                Fmin = f
            if count < n:
                Fmin = f
                g = gr(b)
                gradcount += 1
                it += 1
                t = steplength * t
                c = g - c
                D1 = float(t @ c)
                if D1 > 0:
                    Xv = B @ c
                    D2 = 1.0 + float(Xv @ c) / D1
                    B = B + (D2 * np.outer(t, t) - np.outer(Xv, t) - np.outer(t, Xv)) / D1
                else:
                    ilast = gradcount
            else:
                if ilast < gradcount:
                    count = 0
                    ilast = gradcount
        else:
            count = 0
            if ilast == gradcount:
                count = n
            else:
                ilast = gradcount
        if it >= maxit:
            break
        if gradcount - ilast > 2 * n:
            ilast = gradcount
        if count == n and ilast == gradcount:
            break
    return b, Fmin


def optimhess(gr_ll, par, ndeps=1e-3):
    """optim(hessian = TRUE): central differences of the log-lik gradient."""
    k = par.size
    H = np.empty((k, k))
    for i in range(k):
        d = par.copy()
        d[i] += ndeps
        g1 = gr_ll(d)
        d[i] -= 2 * ndeps
        g2 = gr_ll(d)
        with np.errstate(invalid="ignore"):
            H[i] = (g1 - g2) / (2 * ndeps)
    return 0.5 * (H + H.T)


def _neg_inv(H):
    """-solve(H); pscl's tryCatch turns a failed solve() into an NA matrix."""
    try:
        with np.errstate(invalid="ignore"):
            return -np.linalg.inv(H)
    except np.linalg.LinAlgError:
        return np.full(H.shape, np.nan)


def _rsum(v):
    """R's sum(): sequential accumulation in 80-bit long double."""
    return float(np.cumsum(np.asarray(v, dtype=np.longdouble))[-1]) if np.size(v) else 0.0


def _colsums(M):
    """R's colSums(): per-column sequential long-double accumulation."""
    return np.cumsum(np.asarray(M, dtype=np.longdouble), axis=0)[-1].astype(float)


def _matvec(X, b):
    """X %*% b via reference-BLAS dgemv: y = 0; y += b[j] * X[:, j] column by column."""
    y = np.zeros(X.shape[0])
    for j in range(X.shape[1]):
        y = y + b[j] * X[:, j]
    return y


# ------------------------------------------------------- pscl hurdle (negbin)
_SFERR_HALVES = np.array([
    0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
    0.0413406959554092940938221, 0.03316287351993628748511048, 0.02767792568499833914878929,
    0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
    0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
    0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
    0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
    0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
    0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
    0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
    0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690])
_S0, _S1, _S2, _S3, _S4 = 1 / 12, 1 / 360, 1 / 1260, 1 / 1680, 1 / 1188
_LN_SQRT_2PI = 0.918938533204672741780329736406
_LN_2PI = 1.837877066409345483560659472811


def _stirlerr(n):
    """R 4.3 nmath stirlerr(), vectorised."""
    n = np.asarray(n, float)
    out = np.empty_like(n)
    small = n <= 15.0
    nn2 = n + n
    half = small & (nn2 == np.floor(nn2))
    out[half] = _SFERR_HALVES[nn2[half].astype(int)]
    o = small & ~half
    out[o] = gammaln(n[o] + 1.0) - (n[o] + 0.5) * np.log(n[o]) + n[o] - _LN_SQRT_2PI
    big = ~small
    nb = n[big]
    q = nb * nb
    r = np.where(nb > 500, (_S0 - _S1 / q) / nb,
        np.where(nb > 80, (_S0 - (_S1 - _S2 / q) / q) / nb,
        np.where(nb > 35, (_S0 - (_S1 - (_S2 - _S3 / q) / q) / q) / nb,
                 (_S0 - (_S1 - (_S2 - (_S3 - _S4 / q) / q) / q) / q) / nb)))
    out[big] = r
    return out


def _bd0(x, np_):
    """R 4.3 nmath bd0(): x log(x/np) + np - x, Taylor series near x == np."""
    x, np_ = np.broadcast_arrays(np.asarray(x, float), np.asarray(np_, float))
    out = x * np.log(x / np_) + np_ - x
    near = np.abs(x - np_) < 0.1 * (x + np_)
    if near.any():
        xn, pn = x[near], np_[near]
        v = (xn - pn) / (xn + pn)
        s = (xn - pn) * v
        ej = 2 * xn * v
        v = v * v
        done = np.abs(s) < np.finfo(float).tiny
        for j in range(1, 1000):
            if done.all():
                break
            ej = np.where(done, ej, ej * v)
            s_new = np.where(done, s, s + ej / (2 * j + 1))
            done = done | (s_new == s)
            s = s_new
        out[near] = s
    return out


def _nb_logpmf(y, theta, mu):
    """R 4.3 dnbinom_mu(y, size = theta, mu, log = TRUE) for y > 0."""
    y = np.asarray(y, float)
    mu = np.broadcast_to(np.asarray(mu, float), y.shape)
    out = np.empty(y.shape)
    mm = y < 1e-10 * theta                                   # MM's formula
    if mm.any():
        m = mu[mm]
        p = np.where(
            theta < m, np.log(theta / (1 + theta / m)), np.log(m / (1 + m / theta)))
        yy = y[mm]
        out[mm] = yy * p - m - gammaln(yy + 1) + np.log1p(yy * (yy - 1) / (2 * theta))
    o = ~mm
    if o.any():
        yy, m = y[o], mu[o]
        p = np.where(yy < theta, np.log1p(-yy / (theta + yy)), np.log(theta / (theta + yy)))
        n = yy + theta
        pp, qq = theta / (theta + m), m / (theta + m)
        xs = np.full(yy.shape, float(theta))
        lc = (_stirlerr(n) - _stirlerr(xs) - _stirlerr(n - xs)
              - _bd0(xs, n * pp) - _bd0(n - xs, n * qq))
        lf = _LN_2PI + np.log(xs) + np.log1p(-xs / n)
        out[o] = p + (lc - 0.5 * lf)
    return out


def hurdle_negbin(X, y):
    """pscl::hurdle(dist = "negbin", zero.dist = "binomial", link = "logit")."""
    kx = X.shape[1]
    y1 = y > 0
    X1, Y1 = X[y1], y[y1]
    zbin = y1.astype(float)

    def ll_count(p):
        # R returns NaN/Inf here instead of raising; vmmin then rejects the step.
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            return _ll_count(p)

    def _ll_count(p):
        mu = np.exp(_matvec(X, p[:kx]))[y1]
        th = np.exp(p[kx])
        l0 = log_dnbinom0(th, mu)
        return _rsum(_nb_logpmf(Y1, th, mu)) - _rsum(np.log(1 - np.exp(l0)))  # as pscl

    def gr_count(p):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            return _gr_count(p)

    def _gr_count(p):
        eta = _matvec(X, p[:kx])[y1]
        mu = np.exp(eta)
        th = np.exp(p[kx])
        l0 = log_dnbinom0(th, mu)                           # pnbinom(0, log.p = TRUE)
        logratio = l0 - np.log(-np.expm1(l0))               # - pnbinom(0, lower.tail = FALSE, log.p = TRUE)
        v = (Y1 - mu * (Y1 + th) / (mu + th)) - np.exp(logratio + np.log(th) - np.log(mu + th) + eta)
        r1 = _colsums(v[:, None] * X1)
        r2 = _rsum(digamma(Y1 + th) - digamma(th) + np.log(th) - np.log(mu + th) + 1
                   - (Y1 + th) / (mu + th)
                   + np.exp(logratio) * (np.log(th) - np.log(mu + th) + 1 - th / (mu + th))) * th
        return np.append(r1, r2)

    def ll_zero(p):
        mu = logit_linkinv(_matvec(X, p))
        with np.errstate(divide="ignore"):
            return _rsum(np.log(1 - mu[~y1])) + _rsum(np.log(mu[y1]))

    def gr_zero(p):
        eta = _matvec(X, p)
        mu = logit_linkinv(eta)
        me = logit_mu_eta(eta)
        with np.errstate(divide="ignore"):
            w = np.where(y1, 1 / mu, -1 / (1 - mu))
        return _colsums((w * me)[:, None] * X)

    start_count = glm_fit(X, y.astype(float), "poisson")
    start_zero = glm_fit(X, zbin, "binomial")
    pc, fc = vmmin(lambda p: -ll_count(p), lambda p: -gr_count(p), np.append(start_count, 0.0))
    pz, fz = vmmin(lambda p: -ll_zero(p), lambda p: -gr_zero(p), start_zero)
    vc_c = _neg_inv(optimhess(gr_count, pc))
    vc_z = _neg_inv(optimhess(gr_zero, pz))
    loglik = -fc - fz
    npar = kx + 1 + kx
    n = y.size

    def row(par, vc, j):
        beta = par[j]
        se = math.sqrt(vc[j, j]) if vc[j, j] > 0 else float("nan")   # sqrt(NA / negative) -> NaN
        z = beta / se
        return beta, se, z, 2 * stats.norm.sf(abs(z))
    return {"zero": row(pz, vc_z, 1), "count": row(pc, vc_c, 1),
            "AIC": -2 * loglik + 2 * npar, "BIC": -2 * loglik + math.log(n) * npar}


def spearman(a, b):
    """cor.test(method = "spearman", exact = FALSE)."""
    ra, rb = stats.rankdata(a), stats.rankdata(b)
    r = float(np.corrcoef(ra, rb)[0, 1])
    n = a.size
    tstat = r / math.sqrt((1 - r * r) / (n - 2))
    p = 2 * min(stats.t.cdf(tstat, n - 2), stats.t.sf(tstat, n - 2))
    return r, p


def _signif(x, digits=6):
    if not np.isfinite(x) or x == 0:
        return x
    return float(f"{x:.{digits}g}")


# ---------------------------------------------------------------- Open4Gene
def association_test(dm, gene, peak, celltype, covars):
    X = np.column_stack([np.ones(len(dm)), dm["ATAC"].to_numpy(float)]
                        + [dm[c].to_numpy(float) for c in covars])
    y = dm["RNA"].to_numpy(float)
    h = hurdle_negbin(X, y)
    rho, sp = spearman(dm["ATAC"].to_numpy(float), y)
    zb, zse, zz, zp = h["zero"]
    cb, cse, cz, cp = h["count"]
    r6 = lambda v: round(v, 6) if np.isfinite(v) else v
    return [peak, gene, celltype, len(dm), int((y > 0).sum()), int((dm["ATAC"] > 0).sum()),
            r6(zb), r6(zse), r6(zz), _signif(zp), r6(cb), r6(cse), r6(cz), _signif(cp),
            r6(h["AIC"]), h["BIC"], r6(rho), _signif(sp)]


def _read_matrix(path, names_path):
    if str(path).endswith(".mtx"):
        from scipy.io import mmread
        M = mmread(path).tocsr()
        names = pd.read_csv(names_path, header=None)[0].astype(str).tolist()
        return M, {n: i for i, n in enumerate(names)}
    df = pd.read_csv(path, index_col=0)
    from scipy.sparse import csr_matrix
    return csr_matrix(df.to_numpy()), {n: i for i, n in enumerate(df.index.astype(str))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rna", required=True, help="genes x cells counts (.mtx or dense .csv)")
    ap.add_argument("--rna-genes", help="gene names for .mtx rows")
    ap.add_argument("--atac", required=True, help="peaks x cells counts (.mtx or dense .csv)")
    ap.add_argument("--atac-peaks", help="peak names for .mtx rows")
    ap.add_argument("--meta", required=True, help="cell metadata CSV, cell ids in the first column")
    ap.add_argument("--pairs", required=True, help="CSV with Peak,Gene columns")
    ap.add_argument("--covariates", default="", help="comma-separated metadata columns")
    ap.add_argument("--celltype", default="All", help="'All', 'Each' or one cell type")
    ap.add_argument("--celltype-col", default="Cell_Type")
    ap.add_argument("--binary", action="store_true", help="binarise ATAC (> 0 -> 1)")
    ap.add_argument("--min-cells", type=int, default=5, help="MinNum.Cells")
    ap.add_argument("--out", required=True, help="output TSV")
    a = ap.parse_args(argv)

    covars = [c for c in a.covariates.split(",") if c]
    rna, gidx = _read_matrix(a.rna, a.rna_genes)
    atac, pidx = _read_matrix(a.atac, a.atac_peaks)
    meta = pd.read_csv(a.meta, index_col=0)
    for c in covars:                                   # as.integer(): truncation toward 0
        meta[c] = np.trunc(pd.to_numeric(meta[c])).astype(int)
    pairs = pd.read_csv(a.pairs)
    pairs.columns = ["Peak", "Gene"] + list(pairs.columns[2:])
    pairs = pairs[pairs.Peak.isin(pidx) & pairs.Gene.isin(gidx)]
    ctypes = meta[a.celltype_col].astype(str).to_numpy()
    rows = []
    for peak, gene in zip(pairs.Peak, pairs.Gene):
        dm = meta[covars].copy()
        dm["RNA"] = np.asarray(rna[gidx[gene]].todense()).ravel().astype(int)
        at = np.asarray(atac[pidx[peak]].todense()).ravel().astype(float)
        dm["ATAC"] = np.where(at > 0, 1.0, 0.0) if a.binary else at
        if a.celltype == "All":
            groups = [("All", np.ones(len(dm), bool))]
        elif a.celltype == "Each":
            groups = [(ct, ctypes == ct) for ct in pd.unique(ctypes)]
        else:
            groups = [(a.celltype, ctypes == a.celltype)]
        for ct, m in groups:
            sub = dm[m]
            if ((sub.RNA == 0).sum() > 0 and (sub.RNA > 0).sum() >= a.min_cells
                    and (sub.ATAC > 0).sum() >= a.min_cells):
                rows.append(association_test(sub, gene, peak, ct, covars))
    res = pd.DataFrame(rows, columns=COLS)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    summary = {"n_pairs_tested": int(len(res)),
               "n_zero_p_lt_0.05": int((res["hurdle.Res.zero.p"] < 0.05).sum()),
               "method": "port of hbliu/Open4Gene R/Open4Gene.R (pscl hurdle negbin/logit)"}
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
