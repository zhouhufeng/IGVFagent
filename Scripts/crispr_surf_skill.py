#!/usr/bin/env python3
"""CRISPR-SURF: deconvolve a tiling-screen signal into regulatory regions.

Clean-room implementation of the METHOD described for CRISPR-SURF (Hsu et al.,
Nat Commun 2018; pinellolab/CRISPR-SURF, AGPL-3.0). No source was read or
ported: the upstream package is not on PyPI and its source does not parse
under Python 3 at all ("TabError: inconsistent use of tabs and spaces"), so
installing it the way BEAN is installed was not an option.

THE PROBLEM. In a tiling screen each guide perturbs a WINDOW, not a point:
Cas9 cuts locally, dCas9-KRAB spreads heterochromatin over hundreds of bases.
So a single functional element makes every guide within that window look
active, and the measured signal is the true regulatory profile SMEARED by a
perturbation kernel. Reading peaks straight off per-guide scores therefore
reports a wide blur where the element is narrow, and cannot separate two
nearby elements at all.

THE MODEL. Observed guide scores y are the true profile beta convolved with
that kernel:

    y = A beta + noise,    A[i, j] = kernel(guide_i position - bin_j position)

Deconvolution is solving that inverse problem for beta. It is ill-posed --
many profiles explain the same smeared signal -- so it is regularised with an
L1 penalty (`lam`), which prefers the explanation using the fewest active
bins. That is what turns a blur back into discrete elements.

SIGNIFICANCE comes from the negative-control guides, not from a parametric
assumption: they are resampled to build an empirical null for beta, which is
what makes a claim of significance mean something on this screen rather than
on an idealised one. A screen with no labelled negative controls therefore
gets an uncalibrated result, and the output says so rather than quietly
reporting a p-value that assumes a null nobody measured.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "CRISPRsurf"
DATA_DIR = ROOT / "Data" / "CRISPRsurf"

# Characteristic perturbation range in bp, by nuclease. Cas9 makes a
# double-strand break and its effect is local; CRISPRi/CRISPRa spread
# chromatin marks much further, which is why one number does not serve both.
RANGE_BY_NUCLEASE = {"cas9": 20, "cpf1": 20, "crispri": 250, "crispra": 250}


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def load_guides(path: str, score_cols: "Optional[list[str]]" = None):
    """Guide table -> (positions, scores matrix, classes, chrom)."""
    import numpy as np
    rows = list(csv.DictReader(open(path, newline="")))
    if not rows:
        raise SystemExit(f"{path} is empty.")
    cols = rows[0].keys()

    def pick(*names):
        for n in names:
            for c in cols:
                if c.lower() == n:
                    return c
        return None

    c_chr = pick("chr", "chrom", "chromosome")
    c_start = pick("start", "sgrna_start")
    c_stop = pick("stop", "end", "sgrna_stop")
    c_class = pick("class", "sgrna_type", "classification", "type")
    if not (c_chr and c_start):
        raise SystemExit(
            f"{path} needs chromosome and start columns; found: {list(cols)}")

    if score_cols:
        scols = [c for c in score_cols if c in cols]
        if not scols:
            raise SystemExit(f"none of {score_cols} are columns in {path}")
    else:
        # Any column that looks like a replicate log2FC.
        scols = [c for c in cols
                 if any(k in c.lower() for k in ("log2fc", "lfc", "score",
                                                   "rep", "enrich"))
                 and c not in (c_chr, c_start, c_stop, c_class)]
        if not scols:
            raise SystemExit(
                f"no score columns found in {path}. Pass --score-cols "
                f"explicitly; columns are: {list(cols)}")

    pos, mat, cls, chrom = [], [], [], None
    for r in rows:
        try:
            s = int(float(r[c_start]))
        except (TypeError, ValueError):
            continue
        e = int(float(r[c_stop])) if c_stop and r.get(c_stop) else s + 20
        vals = []
        for c in scols:
            try:
                vals.append(float(r[c]))
            except (TypeError, ValueError):
                vals.append(float("nan"))
        if all(v != v for v in vals):
            continue
        chrom = chrom or r[c_chr]
        pos.append((s + e) // 2)                 # guide midpoint
        mat.append(vals)
        cls.append((r.get(c_class) or "observation").strip().lower())
    if not pos:
        raise SystemExit(f"{path} yielded no usable guides.")
    order = np.argsort(np.asarray(pos))
    return (np.asarray(pos)[order], np.asarray(mat, dtype=float)[order],
            [cls[i] for i in order], chrom, scols)


def _is_negative(c: str) -> bool:
    return any(k in c for k in ("negative", "non-targeting", "nontargeting",
                                 "ntc", "neg_control", "negative_control"))


def build_design(pos, bins, rng: int):
    """A[i, j] = perturbation kernel weight of guide i on bin j.

    A triangular kernel over +/- rng: weight 1 at the guide, falling linearly
    to 0 at the edge of its range. The exact shape matters less than its
    WIDTH, which is the assay's characteristic perturbation length -- that is
    the thing being deconvolved out.
    """
    import numpy as np
    d = np.abs(pos[:, None] - bins[None, :])
    A = 1.0 - d / float(rng)
    np.clip(A, 0.0, None, out=A)
    return A


def deconvolve(pos, y, bins, rng: int, lam: float):
    """Solve y ~ A beta with an L1 penalty. Returns beta."""
    import numpy as np
    from sklearn.linear_model import Lasso
    A = build_design(pos, bins, rng)
    # fit_intercept=False: a genome-wide offset is not a regulatory element,
    # and letting the model absorb signal into one hides real effects.
    m = Lasso(alpha=lam, fit_intercept=False, max_iter=10000, positive=False)
    m.fit(A, y)
    return np.asarray(m.coef_, dtype=float)


def cmd_deconvolve(args) -> int:
    setup_logging()
    import numpy as np

    pos, mat, cls, chrom, scols = load_guides(
        args.guides, args.score_cols.split(",") if args.score_cols else None)
    logging.info("%d guides on %s, %d score column(s): %s",
                 len(pos), chrom, len(scols), ", ".join(scols))

    rng = args.range or RANGE_BY_NUCLEASE.get(args.nuclease.lower(), 20)
    logging.info("perturbation range %d bp (%s)", rng, args.nuclease)

    y = np.nanmean(mat, axis=1)
    bins = np.arange(pos.min() - rng, pos.max() + rng + 1, args.bin_size)
    logging.info("%d bins at %d bp", len(bins), args.bin_size)

    beta = deconvolve(pos, y, bins, rng, args.lam)
    logging.info("beta: %d nonzero of %d bins", int((beta != 0).sum()), len(beta))

    # Empirical null from the negative controls: deconvolve reshuffled
    # negative-control scores at the same positions, so the null carries this
    # screen's own noise, kernel and guide spacing.
    neg = np.asarray([i for i, c in enumerate(cls) if _is_negative(c)])
    calib = ""
    if len(neg) >= args.min_neg:
        rs = np.random.default_rng(args.seed)
        null = []
        for _ in range(args.null_draws):
            yy = y.copy()
            yy[:] = rs.choice(y[neg], size=len(y), replace=True)
            null.append(np.abs(deconvolve(pos, yy, bins, rng, args.lam)))
        null = np.concatenate(null)
        pvals = np.array([(np.sum(null >= abs(b)) + 1) / (null.size + 1)
                          for b in beta])
        calib = f"empirical null from {len(neg)} negative-control guides"
    else:
        pvals = np.ones(len(beta))
        calib = (f"UNCALIBRATED — {len(neg)} negative-control guides "
                 f"(need {args.min_neg}); no p-values computed")
        logging.warning(calib)

    # BH across bins.
    order = np.argsort(pvals)
    m = len(pvals)
    fdr = np.ones(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * m / (rank + 1), 1.0)
        fdr[i] = prev

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_surf"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)

    bg = out / "deconvolved_scores.bedgraph"
    with bg.open("w") as fh:
        for b, v in zip(bins, beta):
            fh.write(f"{chrom}\t{b}\t{b + args.bin_size}\t{v:.6g}\n")

    tbl = out / "deconvolved_bins.tsv"
    with tbl.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "start", "end", "beta", "pvalue", "fdr"])
        for b, v, p, q in zip(bins, beta, pvals, fdr):
            w.writerow([chrom, b, b + args.bin_size, f"{v:.6g}",
                        f"{p:.6g}", f"{q:.6g}"])

    # Merge adjacent significant bins into regions.
    sig = (fdr < args.fdr) & (np.abs(beta) > 0)
    regions, i = [], 0
    while i < len(bins):
        if sig[i]:
            j = i
            while j + 1 < len(bins) and sig[j + 1]:
                j += 1
            seg = beta[i:j + 1]
            regions.append({"chrom": chrom, "start": int(bins[i]),
                            "end": int(bins[j]) + args.bin_size,
                            "n_bins": j - i + 1,
                            "peak_beta": float(seg[np.argmax(np.abs(seg))]),
                            "min_fdr": float(fdr[i:j + 1].min())})
            i = j + 1
        else:
            i += 1
    reg = out / "significant_regions.csv"
    with reg.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chrom", "start", "end", "n_bins",
                                            "peak_beta", "min_fdr"])
        w.writeheader()
        for r in regions:
            w.writerow(r)

    summary = {"guides": int(len(pos)), "chrom": chrom, "range_bp": rng,
               "bin_size": args.bin_size, "lam": args.lam,
               "n_bins": int(len(bins)), "nonzero_bins": int((beta != 0).sum()),
               "negative_controls": int(len(neg)),
               "calibration": calib, "n_significant_regions": len(regions),
               "fdr_threshold": args.fdr}
    js = out / "surf_summary.json"
    js.write_text(json.dumps(summary, indent=2))

    print(f"Report:   {reg}")
    print(f"Wrote:    {bg}")
    print(f"Wrote:    {tbl}")
    print(f"Summary:  {js}")
    print(f"  {len(regions)} significant region(s) at FDR < {args.fdr}")
    print(f"  calibration: {calib}")
    for r in regions[:10]:
        print(f"    {r['chrom']}:{r['start']}-{r['end']}  "
              f"beta={r['peak_beta']:.3g}  FDR={r['min_fdr']:.3g}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent crispr-surf",
        description="Deconvolve a CRISPR tiling screen into regulatory "
                     "regions (CRISPR-SURF method).")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("deconvolve", help="Guide scores -> regulatory regions.")
    d.add_argument("--guides", required=True,
                   help="CSV/TSV with chrom,start[,stop],score cols[,class].")
    d.add_argument("--score-cols", help="Comma list of replicate score columns.")
    d.add_argument("--nuclease", default="cas9",
                   choices=sorted(RANGE_BY_NUCLEASE))
    d.add_argument("--range", type=int,
                   help="Perturbation range in bp; overrides --nuclease.")
    d.add_argument("--bin-size", type=int, default=10)
    d.add_argument("--lam", type=float, default=0.1,
                   help="L1 strength. Higher = sparser, fewer regions.")
    d.add_argument("--fdr", type=float, default=0.05)
    d.add_argument("--null-draws", type=int, default=40)
    d.add_argument("--min-neg", type=int, default=10)
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--label")
    d.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"deconvolve": cmd_deconvolve}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
