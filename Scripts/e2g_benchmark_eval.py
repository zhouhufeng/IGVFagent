#!/usr/bin/env python3
"""Benchmark enhancer-to-gene (E2G) predictions against the ENCODE CRISPR
ground-truth element-gene pairs (Gschwind et al.).

Given a CRISPR benchmark TSV (with `Regulated` labels) and one or more prediction
tables that carry per-(element, gene) scores, this maps each CRISPR-tested pair to
its predicted score (interval overlap of the perturbed element with a predicted
element for the *same gene*, taking the max score), then computes a
precision-recall curve, AUPRC, and precision at 70% recall for every score column.

Reusable across scE2G / ENCODE-rE2G / ABC style tables. No pandas; scipy/sklearn
used only if present (manual fallback otherwise).
"""
from __future__ import annotations
import argparse, gzip, json, sys
from bisect import bisect_left
from collections import defaultdict


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_ground_truth(path):
    """Return list of dicts: chrom, start, end, gene, regulated(bool)."""
    with _open(path) as f:
        # skip leading comment lines
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        ix = {c: i for i, c in enumerate(header)}
        need = ["chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Regulated"]
        for c in need:
            if c not in ix:
                sys.exit(f"ground truth missing column {c!r}; have {header}")
        out = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            try:
                out.append({
                    "chrom": p[ix["chrom"]],
                    "start": int(p[ix["chromStart"]]),
                    "end": int(p[ix["chromEnd"]]),
                    "gene": p[ix["measuredGeneSymbol"]],
                    "regulated": p[ix["Regulated"]].strip().upper() in ("TRUE", "1", "T"),
                })
            except ValueError:
                continue
    return out


def build_prediction_index(path, gene_col, chr_col, start_col, end_col,
                           score_cols, gene_whitelist):
    """gene -> chrom -> sorted list of (start, end, {score_col: float}).
    Only rows whose gene is in gene_whitelist are kept (keeps memory bounded)."""
    idx = defaultdict(lambda: defaultdict(list))
    kept = 0
    with _open(path) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        h = {c: i for i, c in enumerate(header)}
        for c in [gene_col, chr_col, start_col, end_col] + score_cols:
            if c not in h:
                sys.exit(f"prediction table missing column {c!r}; have {header}")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            g = p[h[gene_col]]
            if g not in gene_whitelist:
                continue
            try:
                st, en = int(p[h[start_col]]), int(p[h[end_col]])
            except ValueError:
                continue
            sc = {}
            for c in score_cols:
                try:
                    sc[c] = float(p[h[c]])
                except ValueError:
                    sc[c] = 0.0
            idx[g][p[h[chr_col]]].append((st, en, sc))
            kept += 1
    # sort each interval list by start for overlap search
    for g in idx:
        for c in idx[g]:
            idx[g][c].sort(key=lambda t: t[0])
    return idx, kept


def max_overlap_score(idx, gene, chrom, start, end, score_cols):
    """Max score per column over predicted elements overlapping [start,end] for gene."""
    best = {c: 0.0 for c in score_cols}
    chrom_map = idx.get(gene)
    if not chrom_map:
        return best, False
    intervals = chrom_map.get(chrom)
    if not intervals:
        return best, False
    starts = [t[0] for t in intervals]
    # candidate window: any interval starting at or before `end`
    hi = bisect_left(starts, end + 1)
    matched = False
    for i in range(hi):
        st, en, sc = intervals[i]
        if en >= start and st <= end:  # overlap
            matched = True
            for c in score_cols:
                if sc[c] > best[c]:
                    best[c] = sc[c]
    return best, matched


def pr_metrics(labels, scores):
    """Return dict with auprc, precision@70recall, n_pos, n, and the PR curve.

    Uses a dependency-free step-integration of average precision (identical
    convention to sklearn.average_precision_score: sum of precision * delta-recall),
    keeping this in line with the repo's pure-stdlib concordance tooling."""
    ap, p70, curve = _manual_pr(labels, scores)
    return {
        "auprc": round(ap, 4),
        "precision_at_70_recall": round(p70, 4),
        "n_pairs": len(labels),
        "n_positive": int(sum(labels)),
        "curve": curve,
    }


def _manual_pr(labels, scores):
    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    P = sum(labels)
    if P == 0:
        return 0.0, 0.0, {"precision": [], "recall": []}
    tp = fp = 0
    prev_rec = 0.0
    ap = 0.0
    pcs, rcs = [], []
    for i in order:
        if labels[i]:
            tp += 1
        else:
            fp += 1
        prec = tp / (tp + fp)
        rec = tp / P
        ap += prec * (rec - prev_rec)
        prev_rec = rec
        pcs.append(prec)
        rcs.append(rec)
    p70 = max((pcs[i] for i in range(len(rcs)) if rcs[i] >= 0.70), default=0.0)
    step = max(1, len(pcs) // 200)
    return ap, p70, {"precision": pcs[::step], "recall": rcs[::step]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ground-truth", required=True)
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--gene-col", default="GeneSymbol")
    ap.add_argument("--chr-col", default="ElementChr")
    ap.add_argument("--start-col", default="ElementStart")
    ap.add_argument("--end-col", default="ElementEnd")
    ap.add_argument("--score-cols", default="Score,ABC.Score,ARC.E2G.Score,Score.ignoreTPM")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    score_cols = [c.strip() for c in args.score_cols.split(",") if c.strip()]
    gt = read_ground_truth(args.ground_truth)
    genes = {r["gene"] for r in gt}
    sys.stderr.write(f"[gt] {len(gt)} pairs, {sum(r['regulated'] for r in gt)} positive, "
                     f"{len(genes)} unique genes\n")

    idx, kept = build_prediction_index(
        args.predictions, args.gene_col, args.chr_col, args.start_col,
        args.end_col, score_cols, genes)
    sys.stderr.write(f"[pred] kept {kept} prediction rows for benchmark genes\n")

    labels = [1 if r["regulated"] else 0 for r in gt]
    scored = {c: [] for c in score_cols}
    n_matched = 0
    for r in gt:
        best, matched = max_overlap_score(idx, r["gene"], r["chrom"],
                                          r["start"], r["end"], score_cols)
        n_matched += matched
        for c in score_cols:
            scored[c].append(best[c])

    results = {
        "ground_truth": args.ground_truth,
        "predictions": args.predictions,
        "n_pairs": len(gt),
        "n_positive": int(sum(labels)),
        "n_matched_to_prediction": n_matched,
        "match_rate": round(n_matched / len(gt), 4),
        "methods": {},
    }
    for c in score_cols:
        m = pr_metrics(labels, scored[c])
        results["methods"][c] = m
        sys.stderr.write(f"[metric] {c:20s} AUPRC={m['auprc']:.4f} "
                         f"P@70%recall={m['precision_at_70_recall']:.4f}\n")

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    sys.stderr.write(f"[out] wrote {args.out}\n")


if __name__ == "__main__":
    main()
