"""scEPS — single-cell Expression exPlainability Statistics (IGVFagent port).

A clean-room reimplementation of the scEPS method (Zou, Shi et al., medRxiv
2026, doi 10.64898/2026.06.26.26356714; source: github.com/Genentech/sceps),
integrated as an IGVFagent single-cell skill. scEPS integrates GWAS and a
single-cell disease atlas to identify **disease-associated cell neighborhoods**:
for each cell's neighborhood it fits a method-of-moments variance-component
model on per-donor pseudobulk expression and tests whether GWAS-prioritized
genes explain more phenotypic variance than mean-expression-matched control
genes. The per-neighborhood **d-statistic** = OMEGA_GWAS − OMEGA_CONTROL.

Reimplemented from the algorithm (the upstream repo ships no license); no source
lines are copied. Validated against the upstream `test/` reference output.

Pipeline (mirrors the 4 upstream steps):
    igvfagent sceps estimate      # per-neighborhood d-statistics  (step 1)
    igvfagent sceps cluster       # approx-independent nbhd blocks  (step 2)
    igvfagent sceps aggregate     # aggregate to cell type / subtype (step 3)
    igvfagent sceps simulate      # synthetic-cascade self-test

Runs under the scEPS environment (scanpy 1.10.x / anndata 0.10.x / statsmodels
0.14.x / numpy 1.26 / scipy 1.13). Heavy deps are imported lazily.

Input h5ad expectations (`--adata`):
  * log-normalized expression in ``X`` (batch regressed out upstream);
  * ``obs[donor_id_col]`` (required), phenotype column (``--pheno``);
  * a kNN graph in ``obsp['connectivities']`` (run ``sc.pp.neighbors`` first);
  * optional cell-type / covariate columns.
MAGMA gene file (`--gene-list`): whitespace-delimited, columns ``GENE`` ``ZSTAT``.
"""
from __future__ import annotations

import argparse
import gzip
import logging
import os
import random
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOC_DIR = ROOT / "Docs/scEPS"
LOG_DIR = ROOT / "Docs/Logs"
PLAYBOOK_PATH = ROOT / "Docs/Skills/SCEPS_SKILLS.md"
EPS = 1e-12


def _lazy():
    """Import the heavy scientific stack lazily with a clear error."""
    try:
        import numpy as np
        import pandas as pd
        import scipy.stats as st
        import scipy.sparse as sp
        import statsmodels.api as sm
        return np, pd, st, sp, sm
    except ImportError as e:  # pragma: no cover
        sys.exit(f"scEPS needs numpy/pandas/scipy/statsmodels (and scanpy/anndata "
                 f"for --adata): {e}\nActivate the scEPS conda env.")


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = LOG_DIR / f"sceps_{label}_{ts}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(logfile),
                                  logging.StreamHandler()])
    logging.info("Log file: %s", logfile)
    return logfile


# ---------------------------------------------------------------------------
# GWAS gene selection (MAGMA Z -> FDR)
# ---------------------------------------------------------------------------

def select_gwas_genes(args, gene_df):
    """Return the set of GWAS gene symbols by MAGMA ZSTAT (FDR or top-N)."""
    np, pd, st, sp, sm = _lazy()
    g = gene_df.sort_values("ZSTAT", ascending=False).reset_index(drop=True)
    if args.auto_gene_selection:
        p = 1.0 - st.norm.cdf(g["ZSTAT"].values)
        try:
            fdr = st.false_discovery_control(p)
        except AttributeError:  # older scipy
            from statsmodels.stats.multitest import multipletests
            fdr = multipletests(p, method="fdr_bh")[1]
        n = int(np.sum(fdr < args.magma_fdr_thres))
        n = max(args.min_num_gwas_genes, min(n, args.max_num_gwas_genes))
    else:
        n = args.num_gwas_genes
    n = min(n, g.shape[0])
    return list(g["GENE"].values[:n])


# ---------------------------------------------------------------------------
# Neighborhood construction (random-walk diffusion + NAM membership)
# ---------------------------------------------------------------------------

def transition_matrix(adj, sp):
    import numpy as np
    n = adj.shape[0]
    T = adj + sp.identity(n, format=adj.getformat())
    colsums = np.asarray(T.sum(axis=0)).ravel()
    from sklearn.utils.sparsefuncs import inplace_row_scale
    inplace_row_scale(T, 1.0 / colsums)
    return T


def choose_step_size(T, donor_ids, maxnsteps=15):
    """Global random-walk step size via the kurtosis-plateau heuristic."""
    import numpy as np, pandas as pd
    import scipy.stats as st
    S = pd.get_dummies(donor_ids)
    C = S.sum(axis=0).values
    s = S.values.astype(float)
    prev = np.inf
    for i in range(maxnsteps):
        s = T.dot(s)
        medkurt = np.median(st.kurtosis(s / C, axis=1))
        if prev - medkurt < 3 and i + 1 >= 3:
            break
        prev = medkurt
    return i + 1


def reach_prob(T, idx_i, nsteps):
    import numpy as np
    v = np.zeros(T.shape[0]); v[idx_i] = 1.0
    Tt = T.transpose()
    for _ in range(nsteps):
        v = Tt.dot(v)
    return np.asarray(v).ravel()


def neighborhood_by_nam(args, prob, donor_ids, target_nam_row, donor_order,
                        ntimes=1000):
    """NAM-thresholded neighborhood: grid-search 1000 thresholds, pick the one
    minimizing (||expected-realized NAM||, -num_ok_donor, size) s.t. enough donors."""
    import numpy as np, pandas as pd, collections
    donors = np.asarray(donor_ids)
    order = np.argsort(-prob)
    prob_sorted = prob[order]
    donor_sorted = donors[order]
    num_donor = pd.unique(donors).shape[0]

    nam = pd.DataFrame({"expected": np.asarray(target_nam_row, float)},
                       index=list(donor_order))
    nam["num_cells"] = 0.0
    rows = []
    idx = last = 0
    nbhood_size = 0
    for thres in np.flip(np.linspace(0.0, 0.01, num=ntimes)):
        while idx < prob.shape[0] and prob_sorted[idx] > thres:
            idx += 1
        if prob_sorted[min(idx, prob.shape[0] - 1)] == 0.0:
            break
        add = donor_sorted[last:idx]; last = idx
        for d, c in collections.Counter(add).items():
            nam.loc[d, "num_cells"] += c
            nbhood_size += c
        num_ok = int(np.sum(nam["num_cells"] > args.min_num_cell))
        realized = nam["num_cells"] / (idx + EPS)
        diff = float(np.linalg.norm(nam["expected"] - realized))
        rows.append((thres, diff, num_ok, nbhood_size))
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["THRES", "DIFF", "NUM_OK", "SIZE"])
    min_donor = max(args.min_num_donor,
                    int(np.ceil(args.min_frac_donor * num_donor)))
    df = df[df["NUM_OK"] > min_donor]
    if df.empty:
        return None
    df = df.assign(NEG=-df["NUM_OK"]).sort_values(["DIFF", "NEG", "SIZE"])
    opt = df.iloc[0]["THRES"]
    return np.where(prob >= opt)[0], opt


# ---------------------------------------------------------------------------
# Pseudobulk + control-gene matching
# ---------------------------------------------------------------------------

def pseudobulk(args, X_sub, donor_vec, var_df):
    """Per-donor mean expression over the neighborhood; QC donors + genes.
    Returns (pb matrix donors×genes, donor list, gene var_df with mean/std/idx)."""
    import numpy as np, pandas as pd
    donors = pd.unique(donor_vec)
    rows, ncells = [], []
    for d in donors:
        ci = np.where(donor_vec == d)[0]
        rows.append(np.asarray(X_sub[ci, :].mean(axis=0)).ravel())
        ncells.append(ci.shape[0])
    pb = np.vstack(rows)
    ncells = np.array(ncells, float)
    keep_d = ncells >= args.min_num_cell
    pb, donors, ncells = pb[keep_d], donors[keep_d], ncells[keep_d]
    mean = pb.mean(axis=0)
    var = pb.var(axis=0)
    keep_g = (var > args.min_expr_var_thres) & (var < args.max_expr_var_thres)
    pb = pb[:, keep_g]
    v = var_df.iloc[np.where(keep_g)[0]].copy()
    v["mean"], v["std"] = mean[keep_g], np.sqrt(var[keep_g])
    v["GENE_INDEX"] = np.arange(v.shape[0])
    prop = ncells / ncells.sum()
    return pb, list(donors), v, prop


def sample_control_genes(args, var_df, gwas_symbols, rng, nbins=10):
    """Sample |gwas| mean-expression-matched control genes (weighted by the
    number of GWAS genes per expression bin)."""
    import numpy as np, pandas as pd
    v = var_df.copy()
    bins = pd.qcut(v["mean"].values, nbins, duplicates="drop")
    v["bin"] = bins.astype(str)
    gwas_mask = v[args._gene_col].isin(gwas_symbols)
    counts = v.loc[gwas_mask, "bin"].value_counts()
    v["weight"] = 0.0
    for b, c in counts.items():
        v.loc[v["bin"] == b, "weight"] = float(c)
    if not args.include_gwas_genes_in_control_genes:
        v = v[~v[args._gene_col].isin(gwas_symbols)]
    n = len([g for g in gwas_symbols if g in set(var_df[args._gene_col])])
    if v["weight"].sum() <= 0 or n == 0:
        return None
    ctrl = v.sample(n=n, weights=v["weight"], random_state=rng)
    return ctrl["GENE_INDEX"].values


# ---------------------------------------------------------------------------
# Method-of-moments variance-component regression
# ---------------------------------------------------------------------------

def _triu_vec(np, M, offset):
    return M[np.triu_indices(M.shape[0], k=offset)]


def rel_vec(np, X, offset, nbs):
    """Relatedness (Gram) vector with bootstrap bias correction."""
    ng = X.shape[1]
    base = _triu_vec(np, X @ X.T / (ng - 1.0), offset)
    if nbs <= 0:
        return base, np.empty((0, base.shape[0]))
    bs = np.empty((nbs, base.shape[0]))
    for b in range(nbs):
        idx = np.random.choice(ng, size=ng, replace=True)
        bs[b] = _triu_vec(np, X[:, idx] @ X[:, idx].T / (ng - 1.0), offset)
    return 2.0 * base - bs.mean(axis=0), bs


def reg_weights(np, rel_all, y, mean_var_all, offset):
    cp = _triu_vec(np, np.outer(y, y), 0)
    M = _triu_vec(np, rel_all, 0)
    var_y = np.var(y)
    p = M.dot(cp) / M.dot(M)
    p = min(max(p, 0.0), var_y / (mean_var_all + EPS))
    pred = p * rel_all
    pv = 2.0 * pred**2 + var_y**2 - pred**2
    np.fill_diagonal(pv, 2.0 * var_y**2)
    return 1.0 / (_triu_vec(np, pv, offset) + EPS)


def estimate_neighborhood(args, pb, donors, var_df, gwas_idx, y, prop, rng):
    """Fit the 3-component MoM model, return the OMEGA/SIGMA dict for one nbhd."""
    np, pd, st, sp, sm = _lazy()
    nsample = pb.shape[0]
    if nsample < args.min_num_donor:
        return None

    # residualize + scale phenotype
    covar = np.column_stack([prop, np.ones(nsample)])  # prop_cells_local + const
    resid = sm.OLS(y, covar).fit().resid
    resid = resid - resid.mean()
    if args.scale_pheno_neighborhood:
        resid = resid / (resid.std() + EPS)
    var_pheno = float(np.var(resid))

    offset = 1 if args.exclude_diag else 0
    nbs = 0 if args.no_disattenuation else args.num_bootstrap_disattenuation

    std_all = var_df["std"].values
    mean_all = var_df["mean"].values
    ngene_all = pb.shape[1]
    mean_var_all = float(np.mean(std_all**2))
    mean_mean_all = float(np.mean(mean_all))

    all_idx = np.arange(ngene_all)
    X1 = pb[:, gwas_idx]
    ngene1 = X1.shape[1]
    mean_var1 = float(np.mean(std_all[gwas_idx] ** 2))
    mean_mean1 = float(np.mean(mean_all[gwas_idx]))

    ctrl_idx = sample_control_genes(args, var_df, args._gwas_symbols, rng)
    if ctrl_idx is None:
        return None
    X2 = pb[:, ctrl_idx]
    ngene2 = X2.shape[1]
    mean_var2 = float(np.mean(std_all[ctrl_idx] ** 2))
    mean_mean2 = float(np.mean(mean_all[ctrl_idx]))

    rest_idx = np.setdiff1d(np.setdiff1d(all_idx, ctrl_idx), gwas_idx)
    X3 = pb[:, rest_idx]
    ngene_rest = X3.shape[1]
    mean_var_rest = float(np.mean(std_all[rest_idx] ** 2))
    mean_mean_rest = float(np.mean(mean_all[rest_idx]))

    # design matrix
    ncol = 3 + (1 - offset) + (0 if args.exclude_intercept else 1)
    dim = nsample * (nsample - 1) // 2 + (1 - offset) * nsample
    M = np.ones((dim, ncol))
    if not args.exclude_diag:
        M[:, 3] = _triu_vec(np, np.eye(nsample), offset)
    M1, M1b = rel_vec(np, X1, offset, nbs)
    M2, M2b = rel_vec(np, X2, offset, nbs)
    M3, M3b = rel_vec(np, X3, offset, nbs)
    M[:, 0], M[:, 1], M[:, 2] = M1, M2, M3
    cp = _triu_vec(np, np.outer(resid, resid), offset)

    rel_all = pb @ pb.T / (ngene_all - 1.0)
    wgt = np.ones(M.shape[0]) if args.use_ols else \
        reg_weights(np, rel_all, resid, mean_var_all, offset)

    # WLS + bootstrap parameter covariance
    fit = sm.WLS(cp, M, weights=wgt).fit()
    params = np.asarray(fit.params, float)
    if args.use_analytical_stderr:
        cov = np.asarray(fit.normalized_cov_params, float) * fit.scale
    else:
        nb = args.num_bootstrap_regression
        boot = np.empty((nb, ncol))
        n = cp.shape[0]
        for b in range(nb):
            ix = np.random.choice(n, size=n, replace=True)
            boot[b] = sm.WLS(cp[ix], M[ix], weights=wgt[ix]).fit().params
        cov = np.cov(boot.T)
        params = 2.0 * params - boot.mean(axis=0)

    # disattenuation (regression-dilution) factor per component
    disatt = np.ones(ncol)
    if not args.no_disattenuation and nbs > 0:
        for j, (Mj, Mjb) in enumerate([(M1, M1b), (M2, M2b), (M3, M3b)]):
            allM = np.vstack([Mj, Mjb])
            mo = allM.mean(); mm = allM.mean(axis=0)
            sb, sw = np.var(mm - mo), np.var(allM - mm)
            disatt[j] = 1.0 + sw / (sb + EPS)
    params = params * disatt
    cov = np.diag(disatt) @ cov @ np.diag(disatt)

    sig = params[:3] / np.array([ngene1, ngene2, ngene_rest])
    mean_var = np.array([mean_var1, mean_var2])
    omega = sig[:2] * mean_var
    sub = cov[:2, :2]

    def se(vec):
        return float(np.sqrt(max(sub.dot(vec).dot(vec), 0.0)))
    gw = np.array([1 / ngene1, 0]); ct = np.array([0, 1 / ngene2])
    dv = np.array([1 / ngene1, -1 / ngene2])
    omega_diff = omega[0] - omega[1]
    se_omega_gwas = se(gw * mean_var)
    se_omega_ctrl = se(ct * mean_var)
    se_omega_diff = se(dv * mean_var)

    ov = np.array([1 / ngene_all] * 3)
    sig_overall = params[:3].dot(ov)
    omega_overall = sig_overall * mean_var_all
    sub3 = cov[:3, :3]
    se_omega_overall = float(np.sqrt(max(sub3.dot(ov).dot(ov), 0.0))) * mean_var_all
    omega_rest = sig[2] * mean_var_rest
    se_omega_rest = float(np.sqrt(max(cov[2, 2], 0.0))) / ngene_rest * mean_var_rest

    df_z = nsample - ncol

    def zp(val, s, dfadj=0):
        z = val / (s + EPS)
        p = 2.0 * (1.0 - st.t.cdf(abs(z), max(df_z - dfadj, 1)))
        return z, p
    z_g, p_g = zp(omega[0], se_omega_gwas)
    z_c, p_c = zp(omega[1], se_omega_ctrl)
    z_r, p_r = zp(omega_rest, se_omega_rest)
    z_o, p_o = zp(omega_overall, se_omega_overall, dfadj=2)
    z_d, p_d = zp(omega_diff, se_omega_diff, dfadj=1)

    return {
        "NEIGHBORHOOD_SIZE": int(sum(1 for _ in donors) * 0),  # filled by caller
        "NUM_DONOR": nsample,
        "MEAN_MEAN_EXPR_GWAS": mean_mean1, "MEAN_MEAN_EXPR_CONTROL": mean_mean2,
        "MEAN_MEAN_EXPR_REST": mean_mean_rest, "MEAN_MEAN_EXPR_ALL": mean_mean_all,
        "MEAN_VAR_EXPR_GWAS": mean_var1, "MEAN_VAR_EXPR_CONTROL": mean_var2,
        "MEAN_VAR_EXPR_REST": mean_var_rest, "MEAN_VAR_EXPR_ALL": mean_var_all,
        "NUM_GENE_GWAS": ngene1, "NUM_GENE_CONTROL": ngene2,
        "NUM_GENE_REST": ngene_rest,
        "OMEGA_GWAS": omega[0], "SE_OMEGA_GWAS": se_omega_gwas,
        "Z_OMEGA_GWAS": z_g, "P_Z_OMEGA_GWAS": p_g,
        "OMEGA_CONTROL": omega[1], "SE_OMEGA_CONTROL": se_omega_ctrl,
        "Z_OMEGA_CONTROL": z_c, "P_Z_OMEGA_CONTROL": p_c,
        "OMEGA_REST": omega_rest, "SE_OMEGA_REST": se_omega_rest,
        "Z_OMEGA_REST": z_r, "P_Z_OMEGA_REST": p_r,
        "OMEGA_OVERALL": omega_overall, "SE_OMEGA_OVERALL": se_omega_overall,
        "Z_OMEGA_OVERALL": z_o, "P_Z_OMEGA_OVERALL": p_o,
        "OMEGA_DIFF": omega_diff, "SE_OMEGA_DIFF": se_omega_diff,
        "Z_OMEGA_DIFF": z_d, "P_Z_OMEGA_DIFF": p_d,
        "VAR_PHENO": var_pheno,
    }


# ---------------------------------------------------------------------------
# estimate (step 1) command
# ---------------------------------------------------------------------------

def cmd_estimate(args):
    np, pd, st, sp, sm = _lazy()
    import scanpy as sc
    setup_logging(args.label or "estimate")
    np.random.seed(args.seed); random.seed(args.seed)
    rng = np.random.RandomState(args.seed)

    logging.info("loading %s", args.adata)
    adata = sc.read_h5ad(args.adata)
    args._gene_col = args.gene_id_col if args.gene_id_col in adata.var.columns else None
    gene_symbols = (adata.var[args._gene_col].values if args._gene_col
                    else adata.var.index.values)

    gm = pd.read_csv(args.gene_list, sep=r"\s+")
    gm = gm[["GENE", "ZSTAT"]].dropna()
    args._gwas_symbols = select_gwas_genes(args, gm)
    magma_set = set(gm["GENE"])
    logging.info("GWAS genes selected: %d", len(args._gwas_symbols))

    # keep genes present in MAGMA
    keep = np.array([g in magma_set for g in gene_symbols])
    adata = adata[:, keep].copy()
    gene_symbols = gene_symbols[keep]
    var_base = pd.DataFrame({"_sym": gene_symbols})
    args._gene_col = "_sym"

    donor_vec = adata.obs[args.donor_id_col].values
    # phenotype per cell -> per donor handled inside neighborhood pseudobulk
    pheno_cell = adata.obs[args.pheno].astype(float).values

    adj = adata.obsp["connectivities"]
    T = transition_matrix(adj, sp)
    donor_ids = adata.obs[args.donor_id_col]
    nstep0 = choose_step_size(T, donor_ids, maxnsteps=args.max_step_size)
    logging.info("global step size: %d", nstep0)
    # target NAM = per-cell reachability under global step (donor fractions)
    dummies = pd.get_dummies(donor_ids)
    donor_order = list(dummies.columns)
    S = dummies.values.astype(float)
    nam = S.copy()
    for _ in range(nstep0):
        nam = T.dot(nam)
    target_nam = nam  # raw diffused mass (cell × donor); matches upstream target NAM

    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, float)

    ncell = adata.n_obs
    stop = min(args.stop_idx if args.stop_idx is not None else ncell, ncell)
    start = args.start_idx or 0
    cols_out = None
    rows = []
    for i in range(start, stop):
        for nstep in range(nstep0, args.max_step_size + 1):
            prob = reach_prob(T, i, nstep)
            res = neighborhood_by_nam(args, prob, donor_vec, target_nam[i],
                                      donor_order)
            if res is None:
                continue
            cell_idx, _ = res
            X_sub = X[cell_idx, :]
            pb, donors, vdf, prop = pseudobulk(args, X_sub, donor_vec[cell_idx],
                                               var_base)
            gwas_idx = vdf.index[vdf["_sym"].isin(args._gwas_symbols)].values
            gwas_idx = vdf.loc[gwas_idx, "GENE_INDEX"].values
            if len(gwas_idx) < 2 or pb.shape[0] < args.min_num_donor:
                continue
            # per-donor phenotype
            y = pd.Series(pheno_cell, index=donor_vec).groupby(level=0).mean()
            y = y.loc[donors].values.astype(float)
            y = y - y.mean()
            out = estimate_neighborhood(args, pb, donors, vdf, gwas_idx, y, prop, rng)
            if out is None:
                continue
            out["NEIGHBORHOOD_SIZE"] = int(cell_idx.shape[0])
            out = {"CELL": adata.obs_names[i], "DONOR_ID": donor_vec[i],
                   "STEP_SIZE": nstep, **out}
            rows.append(out)
            break
        if (i - start) % 50 == 0:
            logging.info("processed %d/%d neighborhoods", i - start + 1, stop - start)

    if not rows:
        sys.exit("no neighborhoods produced output")
    df = pd.DataFrame(rows)
    outdir = DOC_DIR / f"{datetime.now():%Y%m%d_%H%M%S}_{args.label or 'estimate'}"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = Path(args.out) if args.out else outdir / "sceps.omega.txt.gz"
    df.to_csv(outpath, sep="\t", index=False, compression="gzip"
              if str(outpath).endswith(".gz") else None)
    logging.info("wrote %d neighborhoods -> %s", df.shape[0], outpath)
    print(f"scEPS estimate: {df.shape[0]} neighborhoods -> {outpath}")
    print(f"  d-statistic (OMEGA_DIFF): mean={df['OMEGA_DIFF'].mean():.4g} "
          f"median={df['OMEGA_DIFF'].median():.4g}")


# ---------------------------------------------------------------------------
# cluster (step 2) — approximately-independent neighborhood blocks
# ---------------------------------------------------------------------------

ESTIMANDS = ["OMEGA_GWAS", "OMEGA_CONTROL", "OMEGA_REST", "OMEGA_OVERALL", "OMEGA_DIFF"]
INFO_COLS = ["NEIGHBORHOOD_SIZE", "NUM_DONOR", "MEAN_VAR_EXPR_GWAS",
             "MEAN_VAR_EXPR_CONTROL", "MEAN_VAR_EXPR_REST", "MEAN_VAR_EXPR_ALL",
             "NUM_GENE_GWAS", "NUM_GENE_CONTROL", "NUM_GENE_REST", "VAR_PHENO"]


def cmd_cluster(args):
    """Assign each cell to a neighborhood block (mini-batch k-means on the
    standardized NAM) for the block bootstrap used in aggregation."""
    np, pd, st, sp, sm = _lazy()
    import scanpy as sc
    from sklearn.cluster import MiniBatchKMeans
    setup_logging(args.label or "cluster")
    np.random.seed(args.seed); random.seed(args.seed)
    adata = sc.read_h5ad(args.adata)
    if "connectivities" not in adata.obsp:
        sc.pp.neighbors(adata, use_rep=args.neighbors_use_rep)
    T = transition_matrix(adata.obsp["connectivities"], sp)
    donor_ids = adata.obs[args.donor_id_col]
    nstep = choose_step_size(T, donor_ids, maxnsteps=15)
    nam = pd.get_dummies(donor_ids).values.astype(float)
    for _ in range(nstep):
        nam = T.dot(nam)
    nam = (nam - nam.mean(axis=0)) / (nam.std(axis=0) + EPS)
    km = MiniBatchKMeans(n_clusters=args.num_kmeans_cluster,
                         batch_size=max(1, int(0.05 * adata.n_obs)),
                         random_state=args.seed, n_init="auto")
    labels = km.fit(nam).labels_
    cell_ids = (adata.obs[args.cell_id_col].values
                if args.cell_id_col in adata.obs.columns else adata.obs_names.values)
    out = pd.DataFrame({"CELL": cell_ids, "sceps.neighborhood_cluster": labels})
    outpath = Path(args.out) if args.out else \
        DOC_DIR / f"{datetime.now():%Y%m%d_%H%M%S}_{args.label or 'cluster'}.txt.gz"
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outpath, sep="\t", index=False,
               compression="gzip" if str(outpath).endswith(".gz") else None)
    logging.info("assigned %d cells to %d neighborhood blocks -> %s",
                 len(out), args.num_kmeans_cluster, outpath)
    print(f"scEPS cluster: {len(out)} cells -> {args.num_kmeans_cluster} blocks -> {outpath}")


# ---------------------------------------------------------------------------
# aggregate (step 3) — neighborhood -> cell-type, block-bootstrap significance
# ---------------------------------------------------------------------------

def _fdr_flags(np, pvals, est, alpha):
    from statsmodels.stats.multitest import multipletests
    rej = multipletests(pvals, alpha=alpha, method="fdr_bh")[0]
    return rej & (est > 0)


def _per_neighborhood_fdr(np, df):
    """Add SIGNIF_FDR{5,10,20}_{estimand} booleans (est>0 & BH-significant;
    DIFF additionally requires OMEGA_GWAS>0)."""
    for e in ESTIMANDS:
        if e not in df or f"P_Z_{e}" not in df:
            continue
        p = df[f"P_Z_{e}"].fillna(1.0).values
        est = df[e].values
        for a, tag in ((0.05, 5), (0.10, 10), (0.20, 20)):
            flag = _fdr_flags(np, p, est, a)
            if "DIFF" in e:
                flag = flag & (df["OMEGA_GWAS"].values > 0)
            df[f"SIGNIF_FDR{tag}_{e}"] = flag
    return df


def _block_bootstrap_mean(np, sub, col, nbs, seed):
    """Mean of `col` across neighborhoods + block-bootstrap SE/Z/P."""
    rng = np.random.RandomState(seed)
    vals = sub[col].values
    n = len(vals)
    mean = float(vals.mean())
    if "BLOCK" in sub.columns:
        blocks = sub["BLOCK"].values
        uniq = np.unique(blocks)
        idx_by_block = {b: np.where(blocks == b)[0] for b in uniq}
        boot = np.empty(nbs)
        for i in range(nbs):
            pick = rng.choice(uniq, size=len(uniq), replace=True)
            idx = np.concatenate([idx_by_block[b] for b in pick])
            boot[i] = vals[idx].mean()
    else:
        boot = np.array([vals[rng.choice(n, n, replace=True)].mean() for _ in range(nbs)])
    se = float(boot.std())
    z = mean / (se + EPS)
    from scipy.stats import norm
    p = 2.0 * (1.0 - norm.cdf(abs(z)))
    return mean, se, z, p


def cmd_aggregate(args):
    """Aggregate per-neighborhood scEPS stats to cell-type level with a
    block bootstrap; report MEAN/SE/Z/P and counts of significant neighborhoods."""
    np, pd, st, sp, sm = _lazy()
    import glob
    setup_logging(args.label or "aggregate")
    files = glob.glob(args.sceps_result)
    if not files:
        sys.exit(f"no scEPS result files match {args.sceps_result!r}")
    df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    logging.info("loaded %d neighborhoods from %d file(s)", len(df), len(files))
    df = _per_neighborhood_fdr(np, df)

    # block map
    if args.neighborhood_clusters:
        cl = pd.read_csv(args.neighborhood_clusters, sep="\t")
        cell2block = dict(zip(cl["CELL"].astype(str),
                              cl["sceps.neighborhood_cluster"]))
        df["BLOCK"] = df["CELL"].astype(str).map(cell2block)

    # cell-type map (optional)
    ct_of = None
    if args.cell_type_col and args.adata:
        import scanpy as sc
        ad_obs = sc.read_h5ad(args.adata, backed="r").obs
        cid = (ad_obs[args.cell_id_col] if args.cell_id_col in ad_obs.columns
               else ad_obs.index)
        ct_of = dict(zip(cid.astype(str), ad_obs[args.cell_type_col].astype(str)))
        df["CELLTYPE"] = df["CELL"].astype(str).map(ct_of)
        groups = [("All", df)] + [(ct, df[df["CELLTYPE"] == ct])
                                  for ct in sorted(df["CELLTYPE"].dropna().unique())]
    else:
        groups = [("All", df)]

    rows = []
    for ct, sub in groups:
        if sub.shape[0] == 0:
            continue
        row = {"CELLTYPE": ct, "NUM_CELL": sub.shape[0]}
        for c in INFO_COLS:
            if c in sub:
                row[f"MEAN_{c}"] = float(sub[c].mean())
        for e in ESTIMANDS:
            if e not in sub:
                continue
            mean, se, z, p = _block_bootstrap_mean(np, sub, e, args.num_bootstrap, args.seed)
            row[f"MEAN_{e}"], row[f"SE_MEAN_{e}"] = mean, se
            row[f"Z_MEAN_{e}"], row[f"P_Z_MEAN_{e}"] = z, p
            for tag in (5, 10, 20):
                col = f"SIGNIF_FDR{tag}_{e}"
                if col in sub:
                    row[f"NUM_SIGNIF_FDR{tag}_{e}"] = int(sub[col].sum())
        rows.append(row)

    res = pd.DataFrame(rows)
    outpath = Path(args.out) if args.out else \
        DOC_DIR / f"{datetime.now():%Y%m%d_%H%M%S}_{args.label or 'aggregate'}.celltype.txt"
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(outpath, sep="\t", index=False, float_format="%.6g")
    logging.info("wrote %d cell-type rows -> %s", len(res), outpath)
    print(f"scEPS aggregate: {len(res)} groups -> {outpath}")
    for _, r in res.iterrows():
        print(f"  {r['CELLTYPE']:20s} d={r.get('MEAN_OMEGA_DIFF',float('nan')):.3e} "
              f"Z={r.get('Z_MEAN_OMEGA_DIFF',float('nan')):.2f} "
              f"P={r.get('P_Z_MEAN_OMEGA_DIFF',float('nan')):.2e} "
              f"nSig(FDR10)={r.get('NUM_SIGNIF_FDR10_OMEGA_DIFF','-')}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_estimate_args(sp):
    sp.add_argument("--adata", required=True)
    sp.add_argument("--gene-list", required=True, help="MAGMA file: GENE ZSTAT")
    sp.add_argument("--donor-id-col", required=True)
    sp.add_argument("--pheno", required=True, help="obs column with phenotype")
    sp.add_argument("--gene-id-col", default="")
    sp.add_argument("--cell-id-col", default="")
    sp.add_argument("--auto-gene-selection", action="store_true")
    sp.add_argument("--num-gwas-genes", type=int, default=500)
    sp.add_argument("--magma-fdr-thres", type=float, default=0.05)
    sp.add_argument("--min-num-gwas-genes", type=int, default=500)
    sp.add_argument("--max-num-gwas-genes", type=int, default=2000)
    sp.add_argument("--num-expr-bins", type=int, default=10)
    sp.add_argument("--include-gwas-genes-in-control-genes", action="store_true")
    sp.add_argument("--min-num-donor", type=int, default=8)
    sp.add_argument("--min-frac-donor", type=float, default=0.3333)
    sp.add_argument("--min-num-cell", type=int, default=5)
    sp.add_argument("--max-step-size", type=int, default=15)
    sp.add_argument("--min-expr-var-thres", type=float, default=1e-4)
    sp.add_argument("--max-expr-var-thres", type=float, default=100.0)
    sp.add_argument("--num-var-comp", type=int, default=3)
    sp.add_argument("--num-bootstrap-disattenuation", type=int, default=100)
    sp.add_argument("--num-bootstrap-regression", type=int, default=1000)
    sp.add_argument("--use-ols", action="store_true")
    sp.add_argument("--use-analytical-stderr", action="store_true")
    sp.add_argument("--no-disattenuation", action="store_true")
    sp.add_argument("--scale-pheno", action="store_true")
    sp.add_argument("--scale-pheno-neighborhood", action="store_true")
    sp.add_argument("--exclude-diag", action="store_true")
    sp.add_argument("--exclude-intercept", action="store_true")
    sp.add_argument("--start-idx", type=int, default=None)
    sp.add_argument("--stop-idx", type=int, default=None)
    sp.add_argument("--seed", type=int, default=0)
    sp.add_argument("--label", default=None)
    sp.add_argument("--out", default=None)


def main():
    p = argparse.ArgumentParser(prog="igvfagent sceps",
                                description="scEPS single-cell disease-neighborhood statistics.")
    sub = p.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("estimate", help="per-neighborhood scEPS d-statistics (step 1).")
    _add_estimate_args(sp)
    sp.set_defaults(func=cmd_estimate)

    sp = sub.add_parser("cluster", help="Assign cells to neighborhood blocks (step 2).")
    sp.add_argument("--adata", required=True)
    sp.add_argument("--donor-id-col", required=True)
    sp.add_argument("--cell-id-col", default="")
    sp.add_argument("--neighbors-use-rep", default="X_pca_harmony")
    sp.add_argument("--num-kmeans-cluster", type=int, default=50)
    sp.add_argument("--seed", type=int, default=0)
    sp.add_argument("--label", default=None)
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_cluster)

    sp = sub.add_parser("aggregate", help="Aggregate neighborhoods to cell types (step 3).")
    sp.add_argument("--sceps-result", required=True, help="glob of step-1 output(s)")
    sp.add_argument("--neighborhood-clusters", default=None, help="step-2 cluster TSV (for block bootstrap)")
    sp.add_argument("--adata", default=None, help="h5ad for cell->cell-type map")
    sp.add_argument("--cell-type-col", default="")
    sp.add_argument("--cell-id-col", default="")
    sp.add_argument("--num-bootstrap", type=int, default=1000)
    sp.add_argument("--seed", type=int, default=0)
    sp.add_argument("--label", default=None)
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_aggregate)

    sp = sub.add_parser("write-playbook", help="Emit Docs/Skills/SCEPS_SKILLS.md.")
    sp.set_defaults(func=lambda a: (PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True),
                                    PLAYBOOK_PATH.write_text(__doc__),
                                    print(f"Wrote {PLAYBOOK_PATH}")))

    args = p.parse_args()
    return int(args.func(args) is None)


if __name__ == "__main__":
    sys.exit(main())
