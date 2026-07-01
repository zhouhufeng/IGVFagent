"""Open4Gene — hurdle-model peak-to-gene linkage (IGVFagent port).

A clean-room Python reimplementation of Open4Gene (Liu et al., Science 2025,
PMID 39913582; R package github.com/hbliu/Open4Gene), integrated as an
IGVFagent single-cell skill. Open4Gene links open-chromatin peaks (snATAC) to
gene expression (snRNA) across cells of a multiome dataset using a two-component
**hurdle model** that accounts for the excess zeros of single-nucleus RNA:

  * zero component  — logistic regression of  I(RNA>0) ~ ATAC + covariates
                      (does the peak predict *whether* the gene is expressed?);
  * count component — zero-truncated negative-binomial regression of the
                      positive RNA counts ~ ATAC + covariates
                      (does the peak predict *how much*?).

The ATAC coefficient in each component (β, se, z, p), plus a Spearman
correlation, are reported per peak-gene pair — optionally per cell type.

Reimplemented from the algorithm (upstream is R/pscl::hurdle); no source lines
copied. The two hurdle components are separable MLEs, fit here with statsmodels
Logit + TruncatedLFNegativeBinomialP — matching pscl's negbin/logit hurdle.

CLI::

    igvfagent open4gene link --rna RNA.mtx --rna-genes genes.txt \\
        --atac ATAC.mtx --atac-peaks peaks.txt --meta Meta.csv \\
        --pairs PeakGene.csv --covariates lognCount_RNA,percent.mt \\
        --celltype All --celltype-col Cell_Type --min-cells 5 --out links.tsv

Runs under a scientific-Python env (numpy/scipy/pandas/statsmodels; scipy.io for
Matrix Market). Heavy deps imported lazily.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DOC_DIR = ROOT / "Docs/Open4Gene"
LOG_DIR = ROOT / "Docs/Logs"
PLAYBOOK_PATH = ROOT / "Docs/Skills/OPEN4GENE_SKILLS.md"

OUT_COLS = ["Peak", "Gene", "Celltype", "TotalCellNum", "ExpressCellNum",
            "OpenCellNum",
            "hurdle.res.zero.beta", "hurdle.res.zero.se", "hurdle.res.zero.z",
            "hurdle.res.zero.p",
            "hurdle.res.count.beta", "hurdle.res.count.se", "hurdle.res.count.z",
            "hurdle.res.count.p",
            "hurdle.AIC", "hurdle.BIC", "spearman.rho", "spearman.p"]


def setup_logging(label):
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    lf = LOG_DIR / f"open4gene_{label}_{ts}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(lf), logging.StreamHandler()])
    return lf


def _lazy():
    try:
        import numpy as np, pandas as pd, scipy.stats as st
        import statsmodels.api as sm
        from statsmodels.discrete.truncated_model import TruncatedLFNegativeBinomialP
        return np, pd, st, sm, TruncatedLFNegativeBinomialP
    except ImportError as e:  # pragma: no cover
        sys.exit(f"open4gene needs numpy/pandas/scipy/statsmodels: {e}")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load_matrix(path, names_path):
    """Load a features×cells matrix (Matrix Market or dense csv) + row names."""
    import numpy as np, pandas as pd
    from scipy.io import mmread
    if str(path).endswith((".mtx", ".mtx.gz")):
        M = mmread(path).tocsr()
    else:
        df = pd.read_csv(path, index_col=0)
        M = df.values
        names_path = names_path or df.index
    names = (pd.read_csv(names_path, header=None)[0].values
             if isinstance(names_path, str) else np.asarray(names_path))
    return M, np.asarray(names)


# ---------------------------------------------------------------------------
# Hurdle model (two separable MLEs)
# ---------------------------------------------------------------------------

def fit_hurdle(np, pd, sm, TruncNB, y, exog):
    """Fit the negbin/logit hurdle. exog includes intercept + ATAC + covariates
    (ATAC is column index 1). Returns dict of zero/count ATAC stats + AIC/BIC."""
    pos = y > 0
    # --- zero component: logit of I(y>0) ---
    zero = sm.Logit(pos.astype(float), exog).fit(disp=0)
    # --- count component: zero-truncated NB on positive counts ---
    count = TruncNB(y[pos], exog[pos], truncation=0).fit(disp=0)
    # ATAC is column 1 (0 = const)
    zb, zse = zero.params[1], zero.bse[1]
    cb, cse = count.params[1], count.bse[1]
    zz, zp = zero.tvalues[1], zero.pvalues[1]
    cz, cp = count.tvalues[1], count.pvalues[1]
    # hurdle log-likelihood = zero-part LL + count-part LL; k = both param sets
    ll = zero.llf + count.llf
    k = exog.shape[1] * 2 + 1  # both linear predictors + NB alpha
    n = y.shape[0]
    aic = -2 * ll + 2 * k
    bic = -2 * ll + k * np.log(n)
    return {"hurdle.res.zero.beta": zb, "hurdle.res.zero.se": zse,
            "hurdle.res.zero.z": zz, "hurdle.res.zero.p": zp,
            "hurdle.res.count.beta": cb, "hurdle.res.count.se": cse,
            "hurdle.res.count.z": cz, "hurdle.res.count.p": cp,
            "hurdle.AIC": aic, "hurdle.BIC": bic}


def association_test(np, pd, st, sm, TruncNB, dm, gene, peak, celltype, covars):
    """Run the hurdle + Spearman for one peak-gene pair on data frame `dm`."""
    y = dm["RNA"].values.astype(float)
    exog = sm.add_constant(dm[["ATAC"] + covars].values.astype(float))
    rho, sp = st.spearmanr(dm["ATAC"].values, dm["RNA"].values)
    try:
        h = fit_hurdle(np, pd, sm, TruncNB, y, exog)
    except Exception as e:
        logging.warning("hurdle fit failed for %s~%s (%s): %s", peak, gene, celltype, e)
        return None
    row = {"Peak": peak, "Gene": gene, "Celltype": celltype,
           "TotalCellNum": dm.shape[0], "ExpressCellNum": int((y > 0).sum()),
           "OpenCellNum": int((dm["ATAC"].values > 0).sum()),
           "spearman.rho": rho, "spearman.p": sp, **h}
    return {k: row.get(k) for k in OUT_COLS}


# ---------------------------------------------------------------------------
# link command
# ---------------------------------------------------------------------------

def cmd_link(args):
    np, pd, st, sm, TruncNB = _lazy()
    setup_logging(args.label or "link")
    covars = [c for c in args.covariates.split(",") if c] if args.covariates else []

    logging.info("loading matrices")
    RNA, genes = _load_matrix(args.rna, args.rna_genes)
    ATAC, peaks = _load_matrix(args.atac, args.atac_peaks)
    meta = pd.read_csv(args.meta, index_col=0)
    for c in covars:
        meta[c] = meta[c].astype(int)
    pairs = pd.read_csv(args.pairs)
    pairs.columns = ["Peak", "Gene"] + list(pairs.columns[2:])
    gene_ix = {g: i for i, g in enumerate(genes)}
    peak_ix = {p: i for i, p in enumerate(peaks)}
    pairs = pairs[pairs["Peak"].isin(peak_ix) & pairs["Gene"].isin(gene_ix)]
    pairs = pairs.drop_duplicates(["Peak", "Gene"]).reset_index(drop=True)
    logging.info("%d peak-gene pairs to test (%s)", pairs.shape[0], args.celltype)

    ctcol = args.celltype_col
    rows = []
    for n, (peak, gene) in enumerate(zip(pairs["Peak"], pairs["Gene"])):
        atac_vec = np.asarray(ATAC[peak_ix[peak], :].todense()).ravel() \
            if hasattr(ATAC, "todense") else ATAC[peak_ix[peak], :]
        rna_vec = np.asarray(RNA[gene_ix[gene], :].todense()).ravel() \
            if hasattr(RNA, "todense") else RNA[gene_ix[gene], :]
        if args.binary:
            atac_vec = (atac_vec > 0).astype(float)
        dm = meta.copy()
        dm["ATAC"] = atac_vec
        dm["RNA"] = rna_vec.astype(int)

        def _ok(d):
            return ((d["RNA"] == 0).sum() > 0 and (d["RNA"] > 0).sum() >= args.min_cells
                    and (d["ATAC"] > 0).sum() >= args.min_cells)

        if args.celltype == "All":
            if _ok(dm):
                r = association_test(np, pd, st, sm, TruncNB, dm, gene, peak, "All", covars)
                if r: rows.append(r)
        elif args.celltype == "Each":
            for ct in dm[ctcol].dropna().unique():
                sub = dm[dm[ctcol] == ct]
                if _ok(sub):
                    r = association_test(np, pd, st, sm, TruncNB, sub, gene, peak, ct, covars)
                    if r: rows.append(r)
        else:
            sub = dm[dm[ctcol] == args.celltype]
            if _ok(sub):
                r = association_test(np, pd, st, sm, TruncNB, sub, gene, peak, args.celltype, covars)
                if r: rows.append(r)
        if (n + 1) % 20 == 0:
            logging.info("processed %d/%d pairs", n + 1, pairs.shape[0])

    if not rows:
        sys.exit("no peak-gene pairs produced output (check --min-cells / inputs)")
    df = pd.DataFrame(rows)[OUT_COLS]
    outdir = DOC_DIR / f"{datetime.now():%Y%m%d_%H%M%S}_{args.label or 'link'}"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = Path(args.out) if args.out else outdir / "open4gene.links.tsv"
    df.to_csv(outpath, sep="\t", index=False)
    logging.info("wrote %d links -> %s", df.shape[0], outpath)
    print(f"Open4Gene: {df.shape[0]} peak-gene links -> {outpath}")
    sig = (df["hurdle.res.zero.p"] < 0.05).sum()
    print(f"  significant zero-component links (p<0.05): {sig}")


def main():
    p = argparse.ArgumentParser(prog="igvfagent open4gene",
                                description="Hurdle-model peak-to-gene linkage (Open4Gene port).")
    sub = p.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("link", help="Peak-to-gene hurdle-model linkage test.")
    sp.add_argument("--rna", required=True, help="features×cells matrix (.mtx or dense .csv)")
    sp.add_argument("--rna-genes", default=None, help="gene names (one per line)")
    sp.add_argument("--atac", required=True)
    sp.add_argument("--atac-peaks", default=None)
    sp.add_argument("--meta", required=True, help="cell metadata CSV (cells in rows)")
    sp.add_argument("--pairs", required=True, help="peak-gene pairs CSV (Peak,Gene)")
    sp.add_argument("--covariates", default="", help="comma list of meta columns")
    sp.add_argument("--celltype", default="All", help="'All' | 'Each' | a cell-type value")
    sp.add_argument("--celltype-col", default="Cell_Type")
    sp.add_argument("--min-cells", type=int, default=5)
    sp.add_argument("--binary", action="store_true", help="binarize ATAC (>0 -> 1)")
    sp.add_argument("--label", default=None)
    sp.add_argument("--out", default=None)
    sp.set_defaults(func=cmd_link)

    sp = sub.add_parser("write-playbook", help="Emit Docs/Skills/OPEN4GENE_SKILLS.md.")
    sp.set_defaults(func=lambda a: (PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True),
                                    PLAYBOOK_PATH.write_text(__doc__),
                                    print(f"Wrote {PLAYBOOK_PATH}")))
    args = p.parse_args()
    return int(args.func(args) is None)


if __name__ == "__main__":
    sys.exit(main())
