# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/Gartner-Lab/deMULTIplex2 @ de48333b89f6
# (R/benchmarking.R, R/classify.R, R/em.R) for zhu2024_demultiplex2_robust. Unreviewed; provenance in Scripts/ported/registry.json.
"""deMULTIplex2 benchmark against ground-truth sample identities.

Port of Gartner-Lab/deMULTIplex2 R/benchmarking.R:
  * ``benchmark_demultiplex2``: run ``demultiplexTags`` with the benchmark
    defaults (init.cos.cut 0.5, max.iter 30, max.cell.fit 1000, prob.cut 0.5,
    quantile fit window 0.05-0.95, seed 1; each can be overridden as the
    per-dataset scripts do) and relabel calls as tag / Negative
    (0 positive tags) / Multiplet (>1);
  * ``confusion_stats``: per true sample, tp / fp / fn / tn, precision, recall
    and F = tp / (tp + (fp + fn) / 2) (0 when undefined), averaged over the
    samples of the tag mapping; plus the doublet-recovery rates.
The classifier itself is IGVFagent's line-by-line port of demultiplexTags /
fit.em / m.step / e.step / MASS::glm.nb (``multiseq demultiplex``).

Several captures can be given; each is demultiplexed on its own and the calls
are pooled before scoring, as the paper does for the lung cell-line and BAL
datasets ("bc_cbn"). Cells are scored on the truth table's cells: a truth cell
missing from the calls is 'NA', which is neither a tag nor Multiplet.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "Benchmarks").is_dir())


def confusion_stats(call_label: pd.Series, true_label: pd.Series, tag_mapping: pd.DataFrame,
                    call_multiplet: str = "Multiplet", true_multiplet: str = "doublet") -> dict:
    """R/benchmarking.R::confusion_stats, same arithmetic and NA handling."""
    call = call_label.reindex(true_label.index).astype(object)
    call = call.where(call.notna(), "NA").astype(str).values
    true = true_label.astype(object).values
    tag_stats = {}
    for tb in tag_mapping["true_label"]:
        tags_tb = set(tag_mapping.loc[tag_mapping["true_label"] == tb, "tag"])
        called = np.isin(call, list(tags_tb))
        is_tb = np.array([t == tb for t in true])
        tp = int(np.sum(called & is_tb))
        fp = int(np.sum(called)) - tp
        fn = int(np.sum(~called & is_tb))
        tn = int(np.sum(~called & ~is_tb))
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f = tp / (tp + 0.5 * (fp + fn)) if tp + fp + fn else 0.0
        tag_stats[tb] = dict(tp=tp, fp=fp, fn=fn, tn=tn, precision=precision,
                             recall=recall, f_score=f)
    avg = {k: float(np.mean([s[k] for s in tag_stats.values()]))
           for k in ("precision", "recall", "f_score")}
    all_tags = set(tag_mapping["tag"])
    is_dbl = np.array([t == true_multiplet for t in true])
    n_dbl = int(is_dbl.sum())
    in_tags = np.isin(call, list(all_tags))

    def rate(mask):
        return float(np.sum(is_dbl & mask) / n_dbl) if n_dbl else float("nan")

    doublet = {"recall": rate(call == call_multiplet),
               "doublet_called_singlet": rate(in_tags),
               "doublet_called_negative": rate(~in_tags & (call != call_multiplet))}
    return {"singlet_avg_stats": avg, "doublet_avg_stats": doublet, "tag_stats": tag_stats}


def benchmark_calls(tag_mtx: pd.DataFrame, *, max_cell_fit: float = 1000, max_iter: int = 30,
                    min_quantile_fit: float = 0.05, max_quantile_fit: float = 0.95,
                    seed: int = 1) -> pd.Series:
    """benchmark_demultiplex2's classifier call and relabelling."""
    from igvfagent.multiseq_analysis_skill import demultiplex_tags
    res = demultiplex_tags(tag_mtx, init_cos_cut=0.5, max_iter=max_iter, converge_threshold=1e-3,
                           prob_cut=0.5, min_cell_fit=10, max_cell_fit=max_cell_fit,
                           min_quantile_fit=min_quantile_fit, max_quantile_fit=max_quantile_fit,
                           residual_type="rqr", seed=seed)
    c = res["classifications"]
    calls = c["barcode_assign"].astype(str).copy()
    calls[c["barcode_count"] == 0] = "Negative"
    calls[c["barcode_count"] > 1] = "Multiplet"
    # R drops zero-count cells before classifying; they have no call there.
    return calls[c["total_tag_umi"] > 0]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="demultiplex2-benchmark", description=__doc__.split("\n")[0])
    ap.add_argument("--tags", nargs="+", required=True,
                    help="cells x tags count CSV(s), first column = cell id; one per capture")
    ap.add_argument("--truth", required=True, help="CSV with columns cell,truth")
    ap.add_argument("--tag-mapping", required=True,
                    help="JSON {tag: true_label} or CSV with columns tag,true_label")
    ap.add_argument("--true-multiplet", default="doublet")
    ap.add_argument("--max-cell-fit", type=float, default=1000)
    ap.add_argument("--max-iter", type=int, default=30)
    ap.add_argument("--min-quantile-fit", type=float, default=0.05)
    ap.add_argument("--max-quantile-fit", type=float, default=0.95,
                    help="0.9 for the paper's Gaublomme run")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--label", default="demultiplex2_benchmark")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args(argv)

    if a.tag_mapping.endswith(".json"):
        m = json.loads(Path(a.tag_mapping).read_text())
        m = m.get("tag_mapping", m)
        tag_mapping = pd.DataFrame({"tag": list(m), "true_label": list(m.values())})
    else:
        tag_mapping = pd.read_csv(a.tag_mapping)[["tag", "true_label"]]
    calls = pd.concat([benchmark_calls(pd.read_csv(p, index_col=0), max_cell_fit=a.max_cell_fit,
                                       max_iter=a.max_iter, min_quantile_fit=a.min_quantile_fit,
                                       max_quantile_fit=a.max_quantile_fit, seed=a.seed)
                       for p in a.tags])
    truth = pd.read_csv(a.truth)
    true_label = pd.Series(truth["truth"].values, index=truth["cell"].astype(str))
    cs = confusion_stats(calls, true_label, tag_mapping, true_multiplet=a.true_multiplet)

    out = Path(a.outdir) if a.outdir else ROOT / "Docs" / "MultiSeq" / f"{time.strftime('%Y%m%d_%H%M%S')}_{a.label}"
    out.mkdir(parents=True, exist_ok=True)
    calls.rename("call").rename_axis("cell").to_csv(out / "calls.csv")
    summary = {"label": a.label, "n_truth_cells": int(len(true_label)), "n_called_cells": int(len(calls)),
               **cs["singlet_avg_stats"], "doublet": cs["doublet_avg_stats"],
               "per_tag": cs["tag_stats"], "max_cell_fit": a.max_cell_fit, "max_iter": a.max_iter,
               "quantile_fit": [a.min_quantile_fit, a.max_quantile_fit],
               "inputs": {"tags": a.tags, "truth": a.truth, "tag_mapping": a.tag_mapping}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"F = {summary['f_score']:.4f}  precision = {summary['precision']:.4f}  "
          f"recall = {summary['recall']:.4f}")
    print(f"Output: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
