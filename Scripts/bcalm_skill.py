#!/usr/bin/env python3
"""BCalm-style barcode-level linear modelling for MPRA.

Clean-room implementation of the approach in BCalm (kircherlab/BCalm), which
extends mpralm to model MPRA activity at BARCODE level. No source ported:
BCalm is R, built on limma, and there is no R runtime here.

WHY BARCODE LEVEL. Every element in an MPRA is represented by many barcodes.
The usual analysis -- including IGVFagent's own mpra_activity -- sums those
barcodes into one count per element and tests that. Summing throws away the
spread AMONG barcodes, which is the best available evidence about how noisy
an element's measurement actually is. An element whose 40 barcodes agree and
one whose 40 barcodes disagree wildly get the same standard error once
summed, and the second is then reported with unearned confidence.

THE MODEL. Per barcode b of element e in replicate r, the activity is
log2((RNA + 1) / (DNA + 1)). Those are averaged per element with a variance
estimated from the barcode spread, and the per-element variances are then
MODERATED by empirical Bayes -- limma's central idea: an element's own
variance estimate is noisy when it has few barcodes, so it is shrunk toward
the trend across all elements.

    posterior_var = (d0 * s0^2 + d * s^2) / (d0 + d)

d0 and s0^2 come from the observed spread of the per-element variances
(method of moments on the scaled-F, rather than limma's full fitFDist). A
moderated t follows, then BH.

WHAT DIFFERS FROM UPSTREAM. limma's fitFDist uses a trigamma-based
Newton solve; this uses method of moments, which agrees closely in the middle
and less so in the extreme tails. voom's precision weights are not applied.
Concordant ranking, not identical p-values.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "BCalm"


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def load_barcode_table(path: str, oligo_col, bc_col):
    delim = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    rows = list(csv.DictReader(open(path), delimiter=delim))
    if not rows:
        raise SystemExit(f"{path} is empty.")
    cols = list(rows[0].keys())

    def pick(given, *names):
        if given:
            if given not in cols:
                raise SystemExit(f"column {given!r} not in {cols}")
            return given
        for n in names:
            for c in cols:
                if c.lower() == n:
                    return c
        return None

    c_ol = pick(oligo_col, "oligo", "element", "name", "insert", "sequence_id")
    c_bc = pick(bc_col, "barcode", "bc")
    if not c_ol:
        raise SystemExit(f"no oligo/element column found in {cols}")
    dna = [c for c in cols if c.lower().startswith("dna")]
    rna = [c for c in cols if c.lower().startswith("rna")]
    if not dna or not rna:
        raise SystemExit(
            f"need DNA* and RNA* count columns; found {cols}. Columns must be "
            f"named DNA_rep1/RNA_rep1 etc.")
    if len(dna) != len(rna):
        raise SystemExit(f"{len(dna)} DNA columns vs {len(rna)} RNA columns; "
                          f"they must pair up by replicate.")
    return rows, c_ol, c_bc, sorted(dna), sorted(rna)


def _fit_f_moments(s2, d):
    """(d0, s0^2) by method of moments on log per-element variances."""
    import numpy as np
    from scipy.special import polygamma
    s2 = np.asarray([v for v in s2 if v == v and v > 0], dtype=float)
    if len(s2) < 3:
        return 0.0, float(np.mean(s2)) if len(s2) else 1.0
    z = np.log(s2)
    ez = z.mean()
    # Var(log s^2) = trigamma(d/2) + trigamma(d0/2); solve for d0.
    target = z.var(ddof=1) - polygamma(1, d / 2.0)
    if target <= 0:
        return float("inf"), float(np.exp(ez + polygamma(0, d / 2.0) - math.log(d / 2.0)))
    lo, hi = 1e-4, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if polygamma(1, mid / 2.0) > target:
            lo = mid
        else:
            hi = mid
    d0 = math.sqrt(lo * hi)
    s0_2 = float(np.exp(ez + polygamma(0, d / 2.0) - math.log(d / 2.0)
                        - polygamma(0, d0 / 2.0) + math.log(d0 / 2.0)))
    return d0, s0_2


def cmd_activity(args) -> int:
    setup_logging()
    import numpy as np
    from scipy import stats

    rows, c_ol, c_bc, dna_cols, rna_cols = load_barcode_table(
        args.counts, args.oligo_col, args.barcode_col)
    logging.info("%d barcode rows · %d replicate pair(s)", len(rows), len(dna_cols))

    per_el: "dict[str, list[float]]" = defaultdict(list)
    n_bc: "dict[str, set]" = defaultdict(set)
    for r in rows:
        el = (r.get(c_ol) or "").strip()
        if not el:
            continue
        bc = (r.get(c_bc) or "").strip() if c_bc else None
        for dc, rc in zip(dna_cols, rna_cols):
            try:
                d = float(r[dc]); n = float(r[rc])
            except (TypeError, ValueError, KeyError):
                continue
            if d < args.min_dna:
                continue
            # +1 on both: a barcode with zero RNA is informative (silent), not
            # undefined, and dropping it biases activity upward.
            per_el[el].append(math.log2((n + 1.0) / (d + 1.0)))
            if bc:
                n_bc[el].add(bc)

    els = [e for e, v in per_el.items() if len(v) >= args.min_obs]
    if not els:
        raise SystemExit(f"no element has >= {args.min_obs} usable "
                          f"barcode x replicate observations")
    logging.info("%d elements with >= %d observations", len(els), args.min_obs)

    means, vars_, ns = [], [], []
    for e in els:
        v = np.asarray(per_el[e], dtype=float)
        means.append(float(v.mean()))
        vars_.append(float(v.var(ddof=1)) if len(v) > 1 else float("nan"))
        ns.append(len(v))
    means = np.asarray(means); vars_ = np.asarray(vars_); ns = np.asarray(ns)

    d_res = float(np.median(ns) - 1)
    d0, s0_2 = _fit_f_moments(vars_, max(d_res, 1.0))
    logging.info("empirical Bayes: d0=%.3g  s0^2=%.4g  (residual df %.0f)",
                 d0, s0_2, d_res)

    # Shrink each element's variance toward the shared prior.
    post = np.empty(len(els))
    for i, (s2, n) in enumerate(zip(vars_, ns)):
        d = max(n - 1, 1)
        if s2 != s2:
            post[i] = s0_2
        elif math.isinf(d0):
            post[i] = s0_2
        else:
            post[i] = (d0 * s0_2 + d * s2) / (d0 + d)
    se = np.sqrt(post / np.maximum(ns, 1))
    tstat = means / np.maximum(se, 1e-12)
    df = ns - 1 + (d0 if not math.isinf(d0) else 0.0)
    pvals = 2.0 * stats.t.sf(np.abs(tstat), np.maximum(df, 1.0))

    order = np.argsort(pvals)
    m = len(pvals)
    fdr = np.ones(m); prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * m / (rank + 1), 1.0)
        fdr[i] = prev

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_bcalm"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "bcalm_activity.tsv"
    idx = np.argsort(pvals)
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["element", "n_obs", "n_barcodes", "log2FC",
                    "raw_var", "moderated_var", "se", "t", "pvalue", "fdr"])
        for i in idx:
            e = els[i]
            w.writerow([e, int(ns[i]), len(n_bc.get(e, ())) or "",
                        f"{means[i]:.6g}",
                        "" if vars_[i] != vars_[i] else f"{vars_[i]:.6g}",
                        f"{post[i]:.6g}", f"{se[i]:.6g}",
                        f"{tstat[i]:.4f}", f"{pvals[i]:.4g}", f"{fdr[i]:.4g}"])

    nsig = int((fdr < args.fdr).sum())
    summary = {"elements": len(els), "observations": int(ns.sum()),
               "replicates": len(dna_cols),
               "eb_d0": None if math.isinf(d0) else round(d0, 4),
               "eb_s0_squared": round(s0_2, 6),
               "significant_fdr": nsig, "fdr_threshold": args.fdr}
    js = out / "bcalm_summary.json"
    js.write_text(json.dumps(summary, indent=2))
    print(f"Report:  {tsv}")
    print(f"Summary: {js}")
    print(f"  {len(els)} elements, {nsig} at FDR < {args.fdr}")
    print(f"  empirical Bayes prior: d0={summary['eb_d0']} s0^2={summary['eb_s0_squared']}")
    for i in idx[:8]:
        print(f"    {els[i]:16} log2FC={means[i]:+.3f} n={ns[i]:3d} "
              f"FDR={fdr[i]:.3g}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent bcalm",
        description="Barcode-level MPRA activity with empirical-Bayes "
                     "moderated statistics (BCalm approach).")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("activity", help="Barcode counts -> per-element activity.")
    a.add_argument("--counts", required=True,
                   help="Barcode-level table with DNA*/RNA* count columns.")
    a.add_argument("--oligo-col")
    a.add_argument("--barcode-col")
    a.add_argument("--min-dna", type=float, default=1.0,
                   help="Skip observations below this DNA count.")
    a.add_argument("--min-obs", type=int, default=3,
                   help="Minimum barcode x replicate observations per element.")
    a.add_argument("--fdr", type=float, default=0.05)
    a.add_argument("--label")
    a.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"activity": cmd_activity}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
