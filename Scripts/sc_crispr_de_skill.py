#!/usr/bin/env python3
"""Single-cell CRISPR differential expression (per-guide NB GLM + rank aggregation).

Clean-room Python reimplementation of the method in **sc-crispr-de**
(Gersbach Lab Bioinformatics, MIT — https://github.com/Gersbachlab-Bioinformatics/sc-crispr-de).
No source was copied: that pipeline is R (Seurat, MASS::glm.nb, brglm2,
GenomicRanges) driven by a SLURM array, and neither R nor SLURM exists in this
container, so a port was never an option. What is reproduced is the METHOD.

The method, in one sentence: every guide is tested independently by comparing
the cells that carry it against the cells with NO detected guide, using a
negative-binomial GLM per gene, and the per-guide p-values are then aggregated
to a gene-level score.

Three stages, matching the upstream shape:

  prepare    10x sgRNA + gene-expression matrices -> filtered AnnData
             checkpoint. (Upstream writes .RData; this writes .h5ad, which is
             the same idea in a format the rest of IGVFagent can already read.)
  test       per-guide NB GLM. Upstream distributes one SLURM array task per
             guide; there is no scheduler here, so guides run across a process
             pool instead. Same unit of work, same independence.
  aggregate  per-guide p-values -> gene-level significance.

WHERE THIS DIFFERS FROM UPSTREAM, stated because a reimplementation that
quietly diverges is worse than one that does not exist:

  * Dispersion. MASS::glm.nb estimates theta by ML with an IRLS/profile
    alternation. statsmodels' discrete NegativeBinomial estimates alpha
    (= 1/theta) by direct ML. Both are ML estimates of the same parameter and
    agree closely, but the optimiser is not the same one and the last digits
    will differ.
  * Bias reduction. Upstream offers --apply-bias-reduction via brglm2
    (Firth-type penalised likelihood). There is no brglm2 equivalent here, so
    that option is NOT implemented rather than silently approximated. Guides
    with complete separation are reported as `separated` instead.
  * Aggregation. Upstream calls the external FRACTEL package. This implements
    alpha-RRA (robust rank aggregation, Kolde et al. 2012, as used by MAGeCK)
    with null calibration against non-targeting guides. That is the same
    FAMILY of method and the same input/output contract, but it is not
    FRACTEL and will not reproduce its numbers exactly. Column names are
    prefixed `rra_`, not `FRACTEL_`, so no downstream reader can mistake one
    for the other.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint          # noqa: E402,F401
import _localstore as ls                                     # noqa: E402

ROOT = ls.ROOT
DATA_DIR = ROOT / "Data" / "scCRISPRde"
REPORT_DIR = ROOT / "Docs" / "scCRISPRde"

# A guide with fewer cells than this on either side of the contrast cannot
# support a GLM; upstream exposes the same idea as
# --min-cells-per-gene-and-group.
DEFAULT_MIN_CELLS = 5


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _write_json(path: Path, obj: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))
    return path


# ─── stage 1: prepare ───────────────────────────────────────────────────────

def load_matrices(gex: str, sgrna: Optional[str] = None):
    """10x directories (or .h5ad) -> one AnnData with guide calls in .obs.

    Accepts either a CellRanger `filtered_feature_bc_matrix` directory or an
    .h5ad. When the sgRNA library is a separate matrix it is read and each
    cell's guides are called from nonzero counts; when guide features live in
    the same matrix (CRISPR Guide Capture), they are split out by feature type.
    """
    import anndata as ad
    import scanpy as sc

    def _read(p: str):
        path = Path(p)
        if path.is_dir():
            return sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        return ad.read_h5ad(path)

    adata = _read(gex)
    adata.var_names_make_unique()

    if sgrna:
        g = _read(sgrna)
        g.var_names_make_unique()
        shared = adata.obs_names.intersection(g.obs_names)
        if len(shared) == 0:
            raise SystemExit(
                "No cell barcodes are shared between the expression matrix and "
                "the sgRNA matrix. They must come from the same run; check "
                "whether one carries a -1 suffix and the other does not.")
        adata = adata[shared].copy()
        g = g[shared].copy()
        guide_names = list(map(str, g.var_names))
        gm = g.X
    else:
        # Guide features inside the same matrix.
        ft = None
        for col in ("feature_types", "feature_type"):
            if col in adata.var.columns:
                ft = adata.var[col].astype(str)
                break
        if ft is None:
            raise SystemExit(
                "No --sgrna matrix given and the expression matrix has no "
                "feature_types column, so guide features cannot be separated. "
                "Pass the sgRNA matrix explicitly.")
        mask = ft.str.contains("CRISPR", case=False, na=False)
        if not mask.any():
            raise SystemExit("No CRISPR guide features found in feature_types.")
        guide_names = list(map(str, adata.var_names[mask.values]))
        gm = adata[:, mask.values].X
        adata = adata[:, ~mask.values].copy()

    import numpy as np
    import scipy.sparse as sp
    gm = sp.csr_matrix(gm) if not sp.issparse(gm) else gm.tocsr()

    n_guides = np.asarray((gm > 0).sum(axis=1)).ravel()
    adata.obs["n_guides"] = n_guides
    # The contrast the method defines: cells carrying exactly one guide are
    # the test groups; cells with NO detected guide are the shared background.
    # Multi-guide cells belong to neither and are kept only so a caller can see
    # how many there were.
    top = np.full(gm.shape[0], "", dtype=object)
    single = np.where(n_guides == 1)[0]
    if len(single):
        idx = np.asarray(gm[single].argmax(axis=1)).ravel()
        for row, j in zip(single, idx):
            top[row] = guide_names[j]
    adata.obs["guide"] = top
    adata.obs["group"] = np.where(n_guides == 0, "background",
                                   np.where(n_guides == 1, "single", "multiplet"))
    adata.obs["nCount_RNA"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["nFeature_RNA"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata.uns["sc_crispr_de"] = {"n_guide_features": len(guide_names)}
    return adata


def cmd_prepare(args) -> int:
    setup_logging()
    import numpy as np
    adata = load_matrices(args.gex, args.sgrna)

    keep = (adata.obs["nFeature_RNA"] >= args.min_genes) & \
           (adata.obs["nCount_RNA"] >= args.min_counts)
    logging.info("cells %d -> %d after QC (min_genes=%d, min_counts=%d)",
                 adata.n_obs, int(keep.sum()), args.min_genes, args.min_counts)
    adata = adata[keep.values].copy()

    gene_cells = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    adata = adata[:, gene_cells >= args.min_cells_per_gene].copy()

    counts = adata.obs["group"].value_counts().to_dict()
    n_bg = int(counts.get("background", 0))
    if n_bg < args.min_cells:
        logging.warning(
            "only %d no-guide background cells (min %d). Every guide is tested "
            "AGAINST this group, so a thin background weakens every result, "
            "not just one.", n_bg, args.min_cells)

    out = Path(args.out) if args.out else (DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{args.label}.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    summary = {"checkpoint": str(out), "n_cells": int(adata.n_obs),
               "n_genes": int(adata.n_vars), "groups": counts,
               "n_guides_seen": int((adata.obs["guide"] != "").sum()),
               "distinct_guides": int(adata.obs.loc[adata.obs["guide"] != "", "guide"].nunique())}
    js = _write_json(REPORT_DIR / f"{Path(out).stem}_prepare.json", summary)
    print(f"Checkpoint: {out}")
    print(f"Summary:    {js}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return 0


# ─── stage 2: per-guide NB GLM ──────────────────────────────────────────────

def _fit_one(y, design):
    """NB GLM for one gene. Returns (coef, pvalue, status)."""
    import numpy as np
    import statsmodels.api as sm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            # alpha (=1/theta) estimated by ML, the same parameter glm.nb
            # profiles for. disp=0 keeps the optimiser quiet.
            res = sm.NegativeBinomial(y, design, loglike_method="nb2").fit(disp=0)
            coef = float(res.params[1])
            pval = float(res.pvalues[1])
            if not np.isfinite(coef) or not np.isfinite(pval):
                return coef, 1.0, "nonfinite"
            return coef, pval, "ok"
        except Exception as e:                              # noqa: BLE001
            return float("nan"), 1.0, f"failed:{type(e).__name__}"


def _test_guide(payload):
    """One guide against background, every gene. Runs in a worker process."""
    import numpy as np
    import scipy.sparse as sp

    (guide, counts_npz, gene_names, treat_idx, bg_idx,
     latent, min_cells) = payload
    X = sp.load_npz(counts_npz).tocsc()
    rows = np.concatenate([treat_idx, bg_idx])
    is_treat = np.concatenate([np.ones(len(treat_idx)), np.zeros(len(bg_idx))])

    design = [np.ones(len(rows)), is_treat]
    for col in latent:
        v = np.asarray(col)[rows].astype(float)
        # log1p on library size: the covariate upstream passes is nCount_RNA,
        # and on the log scale it behaves as the offset a count model wants.
        design.append(np.log1p(v))
    design = np.column_stack(design)

    out = []
    for j, gene in enumerate(gene_names):
        y = np.asarray(X[:, j].todense()).ravel()[rows]
        n_t = int((y[is_treat == 1] > 0).sum())
        n_b = int((y[is_treat == 0] > 0).sum())
        if n_t < min_cells or n_b < min_cells:
            continue
        if y.max() == 0:
            continue
        if n_t == len(treat_idx) and n_b == 0:
            out.append((guide, gene, float("nan"), 1.0, "separated", n_t, n_b))
            continue
        coef, pval, status = _fit_one(y, design)
        out.append((guide, gene, coef, pval, status, n_t, n_b))
    return out


def cmd_test(args) -> int:
    setup_logging()
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    adata = ad.read_h5ad(args.checkpoint)
    genes = list(map(str, adata.var_names))
    if args.gene_whitelist:
        wanted = {l.strip() for l in Path(args.gene_whitelist).read_text().split() if l.strip()}
        keep = [i for i, g in enumerate(genes) if g in wanted]
        if not keep:
            raise SystemExit(f"No gene in {args.gene_whitelist} matched the matrix.")
        adata = adata[:, keep].copy()
        genes = list(map(str, adata.var_names))
    logging.info("testing %d genes", len(genes))

    guides = sorted(set(adata.obs.loc[adata.obs["guide"] != "", "guide"]))
    if args.guides:
        chosen = {g.strip() for g in args.guides.split(",") if g.strip()}
        guides = [g for g in guides if g in chosen]
    if args.max_guides:
        guides = guides[: args.max_guides]
    bg_idx = np.where(adata.obs["group"].values == "background")[0]
    if len(bg_idx) < args.min_cells:
        raise SystemExit(
            f"only {len(bg_idx)} no-guide background cells; every guide is "
            f"tested against this group, so there is nothing to test against.")
    logging.info("%d guides vs %d background cells", len(guides), len(bg_idx))

    # The count matrix goes to workers as a file, not pickled per task: a
    # sparse matrix of this size copied into every worker is how a pool of 8
    # turns into an out-of-memory kill.
    work = DATA_DIR / "_work"
    work.mkdir(parents=True, exist_ok=True)
    npz = work / f"{Path(args.checkpoint).stem}_counts.npz"
    if not npz.exists():
        sp.save_npz(npz, sp.csr_matrix(adata.X))

    latent = []
    for name in (args.latentvar or "").split(","):
        name = name.strip()
        if name and name in adata.obs.columns:
            latent.append(adata.obs[name].values)

    payloads = []
    for g in guides:
        t_idx = np.where(adata.obs["guide"].values == g)[0]
        if len(t_idx) < args.min_cells:
            continue
        payloads.append((g, str(npz), genes, t_idx, bg_idx, latent, args.min_cells))
    logging.info("%d guides have >= %d cells", len(payloads), args.min_cells)

    rows = []
    t0 = time.time()
    workers = max(1, min(args.workers, len(payloads) or 1))
    if workers == 1 or len(payloads) <= 1:
        for p in payloads:
            rows.extend(_test_guide(p))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_test_guide, p): p[0] for p in payloads}
            for i, fut in enumerate(as_completed(futs), 1):
                try:
                    rows.extend(fut.result())
                except Exception as e:                      # noqa: BLE001
                    logging.warning("guide %s failed: %s", futs[fut], e)
                if i % 10 == 0:
                    logging.info("  %d/%d guides", i, len(payloads))

    import csv
    out = Path(args.out) if args.out else (
        DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{args.label}_per_guide.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["guide", "gene", "log2fc", "pvalue", "status",
                    "n_cells_treat_expr", "n_cells_bg_expr"])
        for guide, gene, coef, pval, status, n_t, n_b in rows:
            lfc = "" if coef != coef else f"{coef / math.log(2):.6g}"
            w.writerow([guide, gene, lfc, f"{pval:.6g}", status, n_t, n_b])
    ok = sum(1 for r in rows if r[4] == "ok")
    print(f"Per-guide results: {out}")
    print(f"  tests: {len(rows):,} ({ok:,} converged) over {len(payloads)} guides "
          f"in {time.time() - t0:.1f}s")
    return 0


# ─── stage 3: aggregation ───────────────────────────────────────────────────

def alpha_rra(ranks: "list[float]", alpha: float = 0.25) -> float:
    """alpha-RRA rho for one gene's normalised guide ranks (Kolde 2012).

    Under the null the k normalised ranks are uniform, so the i-th smallest is
    Beta(i, n-i+1). rho is the smallest of those tail probabilities. `alpha`
    restricts the scan to ranks that are actually good, which is what stops one
    mediocre guide out of ten from dragging a real hit down.
    """
    from scipy.stats import beta
    r = sorted(x for x in ranks if x == x)
    if not r:
        return 1.0
    n = len(r)
    rho = 1.0
    for i, ri in enumerate(r, start=1):
        if ri > alpha:
            break
        rho = min(rho, float(beta.cdf(ri, i, n - i + 1)))
    if rho == 1.0 and r:                    # nothing under alpha: use the best
        rho = float(beta.cdf(r[0], 1, n))
    return rho


def cmd_aggregate(args) -> int:
    setup_logging()
    import csv
    import numpy as np

    rows = list(csv.DictReader(Path(args.per_guide).open(), delimiter="\t"))
    rows = [r for r in rows if r.get("status") == "ok"]
    if not rows:
        raise SystemExit(f"{args.per_guide} has no converged tests to aggregate.")

    # Guide -> target gene. Without a map, the guide id itself is the target,
    # which is the convention when guides are named <gene>_<n>.
    gmap = {}
    if args.guide_map:
        for line in Path(args.guide_map).read_text().splitlines():
            parts = line.replace(",", "\t").split("\t")
            if len(parts) >= 2 and parts[0].strip():
                gmap[parts[0].strip()] = parts[1].strip()

    def target_of(guide: str) -> str:
        if guide in gmap:
            return gmap[guide]
        return guide.rsplit("_", 1)[0] if "_" in guide else guide

    ntc = {g.strip() for g in (args.nontargeting or "").split(",") if g.strip()}

    # Rank every test once, globally: RRA needs each guide's standing among
    # ALL tests, not within its own gene.
    rows.sort(key=lambda r: float(r["pvalue"]))
    n = len(rows)
    for i, r in enumerate(rows, start=1):
        r["_nrank"] = i / n

    by_target: "dict[str, list[dict]]" = {}
    null_ranks: "list[float]" = []
    for r in rows:
        tgt = target_of(r["guide"])
        is_ntc = (r["guide"] in ntc or tgt in ntc
                  or any(s in r["guide"].lower() for s in ("non-targeting", "nontargeting", "_ntc", "safe-harbor")))
        if is_ntc:
            null_ranks.append(r["_nrank"])
            continue
        by_target.setdefault(f"{tgt}\t{r['gene']}", []).append(r)

    out_rows = []
    for key, rs in by_target.items():
        tgt, gene = key.split("\t")
        rho = alpha_rra([x["_nrank"] for x in rs], alpha=args.alpha)
        lfcs = [float(x["log2fc"]) for x in rs if x.get("log2fc") not in ("", None)]
        out_rows.append({"target": tgt, "gene": gene, "n_guides": len(rs),
                         "rra_rho": rho,
                         "rra_effect_size": float(np.median(lfcs)) if lfcs else float("nan")})

    # Calibrate against the non-targeting guides when there are enough of them.
    # An uncalibrated rho is a tail probability under an assumed-uniform null;
    # the NTCs measure what the null ACTUALLY looks like in this experiment.
    if len(null_ranks) >= args.min_ntc:
        sizes = {}
        for o in out_rows:
            sizes.setdefault(o["n_guides"], []).append(o)
        rng = np.random.default_rng(args.seed)
        for k, group in sizes.items():
            draws = np.array([alpha_rra(list(rng.choice(null_ranks, size=k, replace=True)),
                                         alpha=args.alpha)
                              for _ in range(args.null_draws)])
            for o in group:
                o["rra_pval"] = float((np.sum(draws <= o["rra_rho"]) + 1) / (len(draws) + 1))
        calib = f"empirical null from {len(null_ranks)} non-targeting tests"
    else:
        for o in out_rows:
            o["rra_pval"] = o["rra_rho"]
        calib = (f"UNCALIBRATED — only {len(null_ranks)} non-targeting tests "
                 f"(need {args.min_ntc}); rho reported as the p-value")
        logging.warning("aggregation %s", calib)

    # Benjamini-Hochberg.
    out_rows.sort(key=lambda o: o["rra_pval"])
    m = len(out_rows)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        q = out_rows[i]["rra_pval"] * m / (i + 1)
        prev = min(prev, q, 1.0)
        out_rows[i]["rra_fdr"] = prev

    out = Path(args.out) if args.out else (
        DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{args.label}_gene_level.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["target", "gene", "n_guides", "rra_rho",
                                       "rra_pval", "rra_fdr", "rra_effect_size"])
        w.writeheader()
        for o in out_rows:
            w.writerow(o)
    sig = sum(1 for o in out_rows if o["rra_fdr"] < 0.05)
    js = _write_json(REPORT_DIR / f"{Path(out).stem}_aggregate.json",
                     {"results": str(out), "n_pairs": m, "significant_fdr05": sig,
                      "calibration": calib, "alpha": args.alpha})
    print(f"Gene-level results: {out}")
    print(f"Summary:            {js}")
    print(f"  {m:,} target-gene pairs, {sig:,} at FDR < 0.05")
    print(f"  calibration: {calib}")
    return 0


def cmd_pipeline(args) -> int:
    rc = cmd_prepare(args)
    if rc:
        return rc
    import glob
    ck = sorted(glob.glob(str(DATA_DIR / f"*_{args.label}.h5ad")))[-1]
    args.checkpoint = ck
    args.out = None
    rc = cmd_test(args)
    if rc:
        return rc
    pg = sorted(glob.glob(str(DATA_DIR / f"*_{args.label}_per_guide.tsv")))[-1]
    args.per_guide = pg
    args.out = None
    return cmd_aggregate(args)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent sc-crispr-de",
        description="Single-cell CRISPR differential expression: per-guide "
                     "negative-binomial GLM against no-guide cells, then rank "
                     "aggregation to gene level.")
    sub = p.add_subparsers(dest="command", required=True)

    def _common(sp):
        sp.add_argument("--label", default="sccrispr")
        sp.add_argument("--out")

    pr = sub.add_parser("prepare", help="10x matrices -> filtered checkpoint.")
    pr.add_argument("--gex", required=True, help="10x dir or .h5ad of gene expression.")
    pr.add_argument("--sgrna", help="10x dir or .h5ad of sgRNA counts (omit if guides are in --gex).")
    pr.add_argument("--min-genes", type=int, default=200)
    pr.add_argument("--min-counts", type=int, default=500)
    pr.add_argument("--min-cells-per-gene", type=int, default=3)
    pr.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS)
    _common(pr)

    te = sub.add_parser("test", help="Per-guide NB GLM against no-guide cells.")
    te.add_argument("--checkpoint", required=True)
    te.add_argument("--latentvar", default="nCount_RNA",
                    help="Comma list of .obs covariates (log1p-transformed).")
    te.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS)
    te.add_argument("--gene-whitelist")
    te.add_argument("--guides", help="Comma list; default every guide.")
    te.add_argument("--max-guides", type=int, default=0)
    te.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    _common(te)

    ag = sub.add_parser("aggregate", help="Per-guide p-values -> gene level (alpha-RRA).")
    ag.add_argument("--per-guide", required=True)
    ag.add_argument("--guide-map", help="guide<TAB>target lines; default splits <gene>_<n>.")
    ag.add_argument("--nontargeting", help="Comma list of NTC guide ids or targets.")
    ag.add_argument("--alpha", type=float, default=0.25)
    ag.add_argument("--null-draws", type=int, default=2000)
    ag.add_argument("--min-ntc", type=int, default=20)
    ag.add_argument("--seed", type=int, default=0)
    _common(ag)

    pl = sub.add_parser("pipeline", help="prepare -> test -> aggregate.")
    pl.add_argument("--gex", required=True)
    pl.add_argument("--sgrna")
    pl.add_argument("--min-genes", type=int, default=200)
    pl.add_argument("--min-counts", type=int, default=500)
    pl.add_argument("--min-cells-per-gene", type=int, default=3)
    pl.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS)
    pl.add_argument("--latentvar", default="nCount_RNA")
    pl.add_argument("--gene-whitelist")
    pl.add_argument("--guides")
    pl.add_argument("--max-guides", type=int, default=0)
    pl.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    pl.add_argument("--guide-map")
    pl.add_argument("--nontargeting")
    pl.add_argument("--alpha", type=float, default=0.25)
    pl.add_argument("--null-draws", type=int, default=2000)
    pl.add_argument("--min-ntc", type=int, default=20)
    pl.add_argument("--seed", type=int, default=0)
    _common(pl)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"prepare": cmd_prepare, "test": cmd_test,
            "aggregate": cmd_aggregate, "pipeline": cmd_pipeline}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
