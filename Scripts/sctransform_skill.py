#!/usr/bin/env python3
"""SCTransform: variance-stabilising normalisation by regularised NB regression.

Clean-room Python implementation of the method in Hafemeister & Satija,
Genome Biology 2019 (Seurat's `SCTransform`, satijalab/seurat, MIT). No source
ported: Seurat is R and there is no R runtime here.

WHY IT EXISTS. Log-normalisation divides by total counts and takes a log. That
assumes every gene's expression scales the same way with sequencing depth, and
it does not: low-expressed genes stay depth-dependent after log-normalising,
so depth leaks into the PCA and clusters partly reflect library size rather
than biology.

THE METHOD.

  1. For each gene, fit a negative-binomial regression of its counts on
     sequencing depth:  log(mu_g) = b0_g + b1_g * log10(total UMI).
  2. REGULARISE. Per-gene estimates are noisy, especially theta for a barely
     detected gene. So the fitted parameters are smoothed across genes as a
     function of mean expression -- genes at similar abundance should have
     similar depth dependence, and pooling them is what makes the fit stable.
     That regularisation step IS the contribution; an unregularised
     gene-by-gene NB fit is a different, worse method.
  3. Pearson residuals against the regularised model:

         z = (y - mu) / sqrt(mu + mu^2 / theta)

     and clip them to +/- sqrt(n/30), Seurat's default, so one outlier cell
     cannot dominate the downstream PCA.

WHAT DIFFERS FROM SEURAT. The regularisation here is a Gaussian kernel
smoother over log10(gene mean) rather than Seurat's `ksmooth` bandwidth
heuristic, theta is estimated by statsmodels rather than MASS::theta.ml, and
genes are fitted directly instead of on Seurat's speed-up subsample. Expect
concordant residuals and the same downstream structure -- not bit-identical
numbers.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "SCTransform"


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _fit_gene(y, log_umi):
    """(b0, b1, theta) for one gene, or None when it cannot be fitted."""
    import numpy as np
    import statsmodels.api as sm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            X = np.column_stack([np.ones(len(y)), log_umi])
            res = sm.NegativeBinomial(y, X, loglike_method="nb2").fit(disp=0)
            b0, b1 = float(res.params[0]), float(res.params[1])
            alpha = float(res.params[-1])
            if not np.isfinite(alpha) or alpha <= 0:
                return None
            theta = 1.0 / alpha            # statsmodels alpha = 1/theta
            if not (np.isfinite(b0) and np.isfinite(b1)):
                return None
            return b0, b1, theta
        except Exception:                                   # noqa: BLE001
            return None


def _kernel_smooth(x, y, bandwidth: float):
    """Gaussian-kernel smooth of y over x. This is the regularisation."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out = np.empty_like(y)
    for i, xi in enumerate(x):
        w = np.exp(-0.5 * ((x - xi) / bandwidth) ** 2)
        s = w.sum()
        out[i] = float((w * y).sum() / s) if s > 0 else y[i]
    return out


def cmd_run(args) -> int:
    setup_logging()
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    adata = ad.read_h5ad(args.input)
    X = sp.csc_matrix(adata.X) if not sp.issparse(adata.X) else adata.X.tocsc()
    n_cells, n_genes = X.shape
    total = np.asarray(X.sum(axis=1)).ravel()
    keep_cells = total > 0
    if not keep_cells.all():
        logging.info("dropping %d cells with zero counts", int((~keep_cells).sum()))
        X = X[keep_cells]
        adata = adata[keep_cells].copy()
        total = total[keep_cells]
        n_cells = X.shape[0]
    log_umi = np.log10(total)

    detected = np.asarray((X > 0).sum(axis=0)).ravel()
    gene_idx = np.where(detected >= args.min_cells)[0]
    logging.info("%d cells, %d genes (%d detected in >= %d cells)",
                 n_cells, n_genes, len(gene_idx), args.min_cells)
    if len(gene_idx) == 0:
        raise SystemExit("no gene passes --min-cells")

    gene_mean = np.asarray(X[:, gene_idx].mean(axis=0)).ravel()
    log_mean = np.log10(np.maximum(gene_mean, 1e-8))

    # Step 1: per-gene NB fits.
    t0 = time.time()
    fits, ok_idx = [], []
    for k, j in enumerate(gene_idx):
        y = np.asarray(X[:, j].todense()).ravel()
        f = _fit_gene(y, log_umi)
        if f is not None:
            fits.append(f)
            ok_idx.append(k)
        if (k + 1) % 500 == 0:
            logging.info("  fitted %d/%d genes", k + 1, len(gene_idx))
    if not fits:
        raise SystemExit("no gene could be fitted")
    fits = np.asarray(fits, dtype=float)
    ok_idx = np.asarray(ok_idx)
    logging.info("fitted %d/%d genes in %.1fs", len(fits), len(gene_idx),
                 time.time() - t0)

    # Step 2: regularise against mean expression. This is the whole point.
    xs = log_mean[ok_idx]
    order = np.argsort(xs)
    xs_o = xs[order]
    reg = np.empty_like(fits)
    for c in range(3):
        vals = fits[order, c]
        if c == 2:                       # theta is smoothed in log space
            vals = np.log10(np.maximum(vals, 1e-8))
        sm_vals = _kernel_smooth(xs_o, vals, args.bandwidth)
        if c == 2:
            sm_vals = 10.0 ** sm_vals
        reg[order, c] = sm_vals

    # Step 3: Pearson residuals against the REGULARISED model.
    clip = float(np.sqrt(n_cells / 30.0))
    logging.info("clipping residuals to +/- %.2f", clip)
    resid = np.zeros((n_cells, len(ok_idx)), dtype=np.float32)
    names = []
    for r, (k, (b0, b1, theta)) in enumerate(zip(ok_idx, reg)):
        j = gene_idx[k]
        y = np.asarray(X[:, j].todense()).ravel()
        mu = np.exp(b0 + b1 * log_umi)
        mu = np.maximum(mu, 1e-8)
        z = (y - mu) / np.sqrt(mu + mu * mu / max(theta, 1e-8))
        resid[:, r] = np.clip(z, -clip, clip)
        names.append(str(adata.var_names[j]))

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_sct"
    out = Path(args.out) if args.out else (OUT_DIR / f"{label}.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    res = ad.AnnData(resid, obs=adata.obs.copy(),
                     var=adata.var.loc[names].copy()
                     if all(n in adata.var_names for n in names)
                     else None)
    res.var_names = names
    res.uns["sctransform"] = {
        "n_cells": int(n_cells), "n_genes_fitted": int(len(ok_idx)),
        "clip": clip, "bandwidth": args.bandwidth,
        "method": "regularised NB regression, Pearson residuals "
                   "(Hafemeister & Satija 2019)",
    }
    res.write_h5ad(out)

    # Residual variance ranks genes: the most variable AFTER depth is removed
    # is what SCTransform offers in place of a log-normalised HVG list.
    rv = resid.var(axis=0)
    top = [names[i] for i in np.argsort(-rv)[:args.n_top]]
    summary = {"output": str(out), "n_cells": int(n_cells),
               "n_genes_fitted": int(len(ok_idx)), "clip": round(clip, 3),
               "top_variable_genes": top[:20]}
    js = OUT_DIR / f"{label}_summary.json"
    js.parent.mkdir(parents=True, exist_ok=True)
    js.write_text(json.dumps(summary, indent=2))
    print(f"Wrote:   {out}")
    print(f"Summary: {js}")
    print(f"  {n_cells} cells x {len(ok_idx)} genes, residuals clipped to +/-{clip:.2f}")
    print(f"  top variable: {', '.join(top[:8])}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent sctransform",
        description="SCTransform: variance-stabilising normalisation by "
                     "regularised negative-binomial regression.")
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="Counts .h5ad -> Pearson residual .h5ad.")
    r.add_argument("--input", required=True, help="Raw-count .h5ad.")
    r.add_argument("--min-cells", type=int, default=5)
    r.add_argument("--bandwidth", type=float, default=0.3,
                   help="Kernel width over log10(gene mean) for regularisation.")
    r.add_argument("--n-top", type=int, default=3000)
    r.add_argument("--label")
    r.add_argument("--out")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"run": cmd_run}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
