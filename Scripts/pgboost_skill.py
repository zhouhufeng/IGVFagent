#!/usr/bin/env python3
"""pgBoost: combine peak-gene link evidence into one calibrated probability.

Clean-room implementation of the approach in pgBoost
(elizabethdorans/pgboost, MIT). No source ported -- pgBoost is R. Unlike the
deep-learning tools on this list, its training data SHIPS with the upstream
repo, so the model can be TRAINED here rather than needing weights nobody
published.

THE PROBLEM IT SOLVES. Several methods score a peak-gene link and they
disagree: Signac, SCENT and Cicero each correlate accessibility with
expression differently, and distance is informative on its own. Picking one is
arbitrary; averaging them ignores that they are not equally trustworthy.
pgBoost learns the combination from fine-mapped eQTL links -- pairs with
independent genetic evidence of a real regulatory relationship -- and returns
one probability.

LEAVE-ONE-CHROMOSOME-OUT, and it is not a detail. Peak-gene links on the same
chromosome share LD structure and often the same eQTL signal, so an ordinary
random CV split puts near-duplicates of a training pair into the test set and
reports an accuracy the model will not reproduce genome-wide. Holding out a
whole chromosome at a time is what makes the estimate honest, and predictions
for each chromosome come from a model that never saw it.

FEATURES ARE BROUGHT, NOT COMPUTED. Signac / SCENT / Cicero are R packages
and are not installed. pgBoost takes their scores as input columns, which is
also how upstream works; `sce2g_predict` in this repo produces a comparable
Kendall correlation plus an ABC share, and those columns can be used directly.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "pgBoost"


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _open(p: str):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def load_table(path: str):
    delim = "," if str(path).rstrip(".gz").endswith(".csv") else "\t"
    rows = list(csv.DictReader(_open(path), delimiter=delim))
    if not rows:
        raise SystemExit(f"{path} is empty.")
    return rows


def _numeric(rows, cols):
    import numpy as np
    X = np.full((len(rows), len(cols)), np.nan)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            v = r.get(c)
            try:
                X[i, j] = float(v)
            except (TypeError, ValueError):
                pass
    return X


def pick_chrom(rows) -> str:
    for c in rows[0]:
        if c.lower() in ("chr", "chrom", "chromosome", "peak_chr"):
            return c
    raise SystemExit(
        "no chromosome column found; leave-one-chromosome-out needs one "
        f"(looked for chr/chrom/chromosome). Columns: {list(rows[0])}")


def cmd_train(args) -> int:
    setup_logging()
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score, average_precision_score

    train = load_table(args.training_file)
    label_col = args.label_col
    if label_col not in train[0]:
        cand = [c for c in train[0] if c.lower() in
                ("label", "y", "positive", "is_positive", "truth")]
        if not cand:
            raise SystemExit(
                f"no label column {label_col!r} in the training file. "
                f"Columns: {list(train[0])}")
        label_col = cand[0]
    preds = ([c.strip() for c in args.predictors.split(",")] if args.predictors
             else [c for c in train[0]
                   if c != label_col and _looks_numeric(train, c)])
    if not preds:
        raise SystemExit("no numeric predictor columns found")
    logging.info("%d training rows · %d predictors: %s",
                 len(train), len(preds), ", ".join(preds))

    chrom_col = pick_chrom(train)
    y = np.array([1 if str(r[label_col]).strip() in ("1", "True", "true", "TRUE")
                  else 0 for r in train])
    X = _numeric(train, preds)
    chroms = np.array([str(r[chrom_col]) for r in train])
    if y.sum() == 0 or y.sum() == len(y):
        raise SystemExit(f"label column {label_col!r} is constant; nothing to learn")
    logging.info("%d positives / %d negatives over %d chromosomes",
                 int(y.sum()), int((1 - y).sum()), len(set(chroms)))

    # Leave-one-chromosome-out: each fold's model never sees its test
    # chromosome, so the score is not inflated by LD-shared near-duplicates.
    oof = np.full(len(y), np.nan)
    for c in sorted(set(chroms)):
        te = chroms == c
        tr = ~te
        if y[tr].sum() == 0 or (1 - y[tr]).sum() == 0:
            continue
        m = HistGradientBoostingClassifier(
            max_iter=args.n_estimators, learning_rate=args.learning_rate,
            max_depth=args.max_depth if args.max_depth > 0 else None,
            random_state=args.seed)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(oof)
    auc = roc_auc_score(y[ok], oof[ok]) if len(set(y[ok])) > 1 else float("nan")
    auprc = average_precision_score(y[ok], oof[ok]) if len(set(y[ok])) > 1 else float("nan")
    logging.info("leave-one-chromosome-out: AUROC %.4f · AUPRC %.4f", auc, auprc)

    final = HistGradientBoostingClassifier(
        max_iter=args.n_estimators, learning_rate=args.learning_rate,
        max_depth=args.max_depth if args.max_depth > 0 else None,
        random_state=args.seed)
    final.fit(X, y)

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_pgboost"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    import pickle
    mp = out / "pgboost_model.pkl"
    with mp.open("wb") as fh:
        pickle.dump({"model": final, "predictors": preds}, fh)
    js = out / "pgboost_training.json"
    js.write_text(json.dumps({
        "training_rows": len(train), "positives": int(y.sum()),
        "predictors": preds, "chromosomes": sorted(set(chroms)),
        "loco_auroc": None if auc != auc else round(float(auc), 4),
        "loco_auprc": None if auprc != auprc else round(float(auprc), 4),
        "model": str(mp),
    }, indent=2))
    print(f"Model:   {mp}")
    print(f"Summary: {js}")
    print(f"  leave-one-chromosome-out AUROC={auc:.4f}  AUPRC={auprc:.4f}")
    print(f"  {int(y.sum())} positives / {len(y)} rows over "
          f"{len(set(chroms))} chromosomes")
    return 0


def _looks_numeric(rows, col, need: int = 5) -> bool:
    n = 0
    for r in rows[:200]:
        try:
            float(r[col]); n += 1
        except (TypeError, ValueError, KeyError):
            pass
    return n >= min(need, len(rows))


def cmd_predict(args) -> int:
    setup_logging()
    import pickle
    import numpy as np
    with open(args.model, "rb") as fh:
        bundle = pickle.load(fh)
    model, preds = bundle["model"], bundle["predictors"]
    rows = load_table(args.data_file)
    missing = [c for c in preds if c not in rows[0]]
    if missing:
        raise SystemExit(
            f"the model needs predictor columns {missing}, which are not in "
            f"{args.data_file}. Columns present: {list(rows[0])}")
    X = _numeric(rows, preds)
    p = model.predict_proba(X)[:, 1]
    order = np.argsort(np.argsort(p))
    pct = 100.0 * order / max(len(p) - 1, 1)

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_pgboost_pred"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "pgboost_predictions.tsv"
    cols = list(rows[0]) + ["pgBoost_probability", "pgBoost_percentile"]
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for r, pr, pc in zip(rows, p, pct):
            r = dict(r)
            r["pgBoost_probability"] = f"{pr:.6f}"
            r["pgBoost_percentile"] = f"{pc:.2f}"
            w.writerow(r)
    print(f"Report: {tsv}")
    print(f"  {len(rows):,} links scored · {int((p >= 0.5).sum()):,} at p >= 0.5")
    top = np.argsort(-p)[:8]
    for i in top:
        ident = " ".join(str(rows[i].get(k, "")) for k in list(rows[0])[:2])
        print(f"    {ident:34} p={p[i]:.4f}  pct={pct[i]:.1f}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent pgboost",
        description="Combine peak-gene link evidence into one calibrated "
                     "probability (gradient boosting, leave-one-chromosome-out).")
    sub = p.add_subparsers(dest="command", required=True)
    t = sub.add_parser("train", help="Labelled links -> model + LOCO metrics.")
    t.add_argument("--training-file", required=True)
    t.add_argument("--label-col", default="label")
    t.add_argument("--predictors", help="Comma list; numeric columns if omitted.")
    t.add_argument("--n-estimators", type=int, default=200)
    t.add_argument("--learning-rate", type=float, default=0.1)
    t.add_argument("--max-depth", type=int, default=0)
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--label"); t.add_argument("--out-dir")
    q = sub.add_parser("predict", help="Candidate links -> probability + percentile.")
    q.add_argument("--data-file", required=True)
    q.add_argument("--model", required=True)
    q.add_argument("--label"); q.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"train": cmd_train, "predict": cmd_predict}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
