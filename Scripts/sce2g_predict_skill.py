#!/usr/bin/env python3
"""scE2G-style enhancer-gene prediction from paired single-cell ATAC + RNA.

Clean-room implementation of the FEATURE COMPUTATION behind scE2G (Sheth et
al.; EngreitzLab/scE2G, MIT). No source ported.

WHAT THIS IS, AND IS NOT. scE2G's published predictor is a model TRAINED on
CRISPR-validated enhancer-gene pairs; its scientific asset is the fitted
weights, not the code. Those weights are not reproduced here, so this does not
emit "the scE2G score" and does not pretend to. What it computes is the
feature set that model reads, plus a transparent combined score:

  * kendall   rank correlation between a peak's accessibility and a gene's
              expression ACROSS CELLS. This is the single-cell signal ABC
              cannot see -- ABC uses bulk activity and distance, so it cannot
              tell a peak that co-varies with the gene from one that merely
              sits nearby and is busy.
  * abc       activity x contact / sum, from the ABC model (abc_skill).
  * distance  TSS to peak midpoint.

Kendall rather than Pearson deliberately: single-cell counts are sparse,
zero-inflated and far from normal, so a rank statistic is the honest choice
and is what scE2G uses.

IGVFagent already CONSUMED scE2G predictions (sce2g_kg_pull ingests them from
the Catalog, and benchmark #12 scores ENCODE-rE2G). This is the other
direction: generating candidate links from your own paired data.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "scE2G"
DEFAULT_WINDOW = 1_000_000


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def read_bed(path: str):
    out = []
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            rec = {"chrom": f[0], "start": int(f[1]), "end": int(f[2])}
        except ValueError:
            continue
        rec["name"] = f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}"
        if len(f) > 5 and f[5] in ("+", "-"):
            rec["strand"] = f[5]
        out.append(rec)
    return out


def metacells(A, R, k: int, seed: int):
    """Pool cells into k metacells.

    Single-cell ATAC is close to binary and desperately sparse: a correlation
    across raw cells is mostly measuring dropout. Pooling is what makes the
    rank statistic informative, and it is what scE2G does too. Pooling is by
    RNA-profile similarity so that metacells are cell-state coherent rather
    than arbitrary.
    """
    import numpy as np
    from sklearn.cluster import KMeans
    n = R.shape[0]
    k = max(2, min(k, n))
    Xr = np.log1p(R / np.maximum(R.sum(axis=1, keepdims=True), 1) * 1e4)
    if Xr.shape[1] > 50:
        Xc = Xr - Xr.mean(axis=0)
        _, _, vt = np.linalg.svd(Xc, full_matrices=False)
        Xr = Xc @ vt[:min(50, vt.shape[0])].T
    lab = KMeans(n_clusters=k, n_init=4, random_state=seed).fit_predict(Xr)
    Am = np.zeros((k, A.shape[1])); Rm = np.zeros((k, R.shape[1]))
    for i in range(k):
        m = lab == i
        if m.any():
            Am[i] = np.asarray(A[m].sum(axis=0)).ravel()
            Rm[i] = np.asarray(R[m].sum(axis=0)).ravel()
    # Depth-normalise, or the correlation is dominated by metacell size.
    Am = Am / np.maximum(Am.sum(axis=1, keepdims=True), 1) * 1e4
    Rm = Rm / np.maximum(Rm.sum(axis=1, keepdims=True), 1) * 1e4
    return Am, Rm


def cmd_predict(args) -> int:
    setup_logging()
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    from scipy.stats import kendalltau

    rna = ad.read_h5ad(args.rna)
    atac = ad.read_h5ad(args.atac)
    shared = rna.obs_names.intersection(atac.obs_names)
    if len(shared) < 10:
        raise SystemExit(
            f"only {len(shared)} shared cell barcodes between --rna and "
            f"--atac. They must be the SAME cells (multiome or matched).")
    rna = rna[shared].copy()
    atac = atac[shared].copy()
    logging.info("%d shared cells · %d genes · %d peaks",
                 len(shared), rna.n_vars, atac.n_vars)

    R = np.asarray(rna.X.todense() if sp.issparse(rna.X) else rna.X, dtype=float)
    A = np.asarray(atac.X.todense() if sp.issparse(atac.X) else atac.X, dtype=float)
    Am, Rm = metacells(A, R, args.metacells, args.seed)
    logging.info("pooled into %d metacells", Am.shape[0])

    peaks = read_bed(args.peaks)
    genes = read_bed(args.genes)
    # Peaks in the BED must line up with the ATAC matrix columns.
    if len(peaks) != atac.n_vars:
        logging.warning("--peaks has %d rows but the ATAC matrix has %d "
                        "columns; matching by name", len(peaks), atac.n_vars)
    avar = {str(v): i for i, v in enumerate(atac.var_names)}
    gvar = {str(v): i for i, v in enumerate(rna.var_names)}

    rows = []
    for g in genes:
        gi = gvar.get(g["name"])
        if gi is None:
            continue
        tss = g["start"] if g.get("strand") != "-" else g["end"]
        y = Rm[:, gi]
        if np.all(y == y[0]):
            continue
        for p in peaks:
            if p["chrom"] != g["chrom"]:
                continue
            mid = (p["start"] + p["end"]) // 2
            dist = abs(mid - tss)
            if dist > args.window:
                continue
            pi = avar.get(p["name"])
            if pi is None:
                continue
            x = Am[:, pi]
            if np.all(x == x[0]):
                continue
            tau, pv = kendalltau(x, y)
            if tau != tau:
                continue
            activity = float(x.mean())
            contact = max(dist, 1000) ** args.gamma
            rows.append({"gene": g["name"], "peak": p["name"],
                         "chrom": p["chrom"], "peak_start": p["start"],
                         "peak_end": p["end"], "distance": dist,
                         "kendall_tau": round(float(tau), 5),
                         "kendall_p": f"{pv:.4g}",
                         "atac_activity": round(activity, 5),
                         "abc_numerator": float(activity * contact)})

    # ABC share, per gene, over the same candidate set.
    by_gene: "dict[str, list]" = {}
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)
    for g, rs in by_gene.items():
        tot = sum(r["abc_numerator"] for r in rs) or 1.0
        for r in rs:
            r["abc"] = round(r["abc_numerator"] / tot, 6)
            del r["abc_numerator"]
            # A transparent combination, NOT the trained scE2G model: the
            # positive part of the correlation scaled by the ABC share.
            # Stated plainly so nobody reports it as an scE2G score.
            r["combined_score"] = round(max(r["kendall_tau"], 0.0) * r["abc"], 6)

    rows.sort(key=lambda r: -r["combined_score"])
    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_sce2g"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "sce2g_features.tsv"
    with tsv.open("w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow(r)
    js = out / "sce2g_summary.json"
    js.write_text(json.dumps({
        "pairs": len(rows), "genes": len(by_gene), "metacells": int(Am.shape[0]),
        "window": args.window, "gamma": args.gamma,
        "note": "features + a transparent combined score; NOT the trained "
                "scE2G model's output",
    }, indent=2))
    print(f"Report:  {tsv}")
    print(f"Summary: {js}")
    print(f"  {len(rows):,} peak-gene pairs over {len(by_gene)} genes")
    for r in rows[:8]:
        print(f"    {r['gene']:10} {r['peak']:22} tau={r['kendall_tau']:+.3f} "
              f"abc={r['abc']:.3f} score={r['combined_score']:.4f} d={r['distance']:,}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent sce2g-predict",
        description="scE2G-style enhancer-gene features from paired "
                     "single-cell ATAC + RNA.")
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("predict", help="Paired ATAC+RNA -> peak-gene features.")
    q.add_argument("--rna", required=True, help="Gene x cell .h5ad.")
    q.add_argument("--atac", required=True, help="Peak x cell .h5ad, same cells.")
    q.add_argument("--peaks", required=True, help="Peak BED; names match ATAC var.")
    q.add_argument("--genes", required=True, help="Gene TSS BED; names match RNA var.")
    q.add_argument("--metacells", type=int, default=50)
    q.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    q.add_argument("--gamma", type=float, default=-0.87)
    q.add_argument("--seed", type=int, default=0)
    q.add_argument("--label")
    q.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"predict": cmd_predict}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
