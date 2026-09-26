#!/usr/bin/env python3
"""Score the two optional follow-up real-data extensions to the Martyn 2025
Variant-EFFECTS benchmark (see OPERATIONS.md "Extending this benchmark"):

1. PPIF splice-site element — 3 untransfected-control replicate files, each
   holding a subset of the paper's "three edits disrupting the splice donor
   motif at the first 5' splice site of the PPIF gene" (paper's own wording;
   Fig. 1e/f). Combines the 3 replicates by variant ID and derives a percent
   expression change per edit from log2_fold_change, to compare against the
   paper's stated "-80% to -52%" range for these edits.

2. lentiMPRA vs. endogenous PPIF promoter — joins the 41 endogenous
   Variant-EFFECTS PPIF-promoter variants (real_data/PPIF_promoter_GRCh38.tsv)
   against the lentiMPRA reporter-variant file by the shared `variant` ID
   (identical NC_000010.11:pos:ref:alt encoding in both files) and computes
   the Pearson correlation, to compare against the paper's own stated
   "Pearson's r=0.54" between the two assays.

Reads/writes summary.json in the same timestamped run directory
analyze_real_data.py wrote to, adding two new top-level keys without
touching the existing "elements" / "ppif_enhancer_tss_distance_bp" keys.
"""
import csv
import json
import sys
from pathlib import Path

SPLICE_FILES = [
    "PPIF_splice_repA_GRCh38.tsv",
    "PPIF_splice_repB_GRCh38.tsv",
    "PPIF_splice_repC_GRCh38.tsv",
]
LENTIMPRA_FILE = "PPIF_promoter_lentiMPRA_GRCh38.tsv"
ENDOGENOUS_PROMOTER_FILE = "PPIF_promoter_GRCh38.tsv"


def score_splice_site(run_dir: Path) -> dict:
    by_variant: dict[str, list[float]] = {}
    n_reps_present = 0
    for fname in SPLICE_FILES:
        p = run_dir / fname
        if not p.is_file():
            print(f"missing {p}", file=sys.stderr)
            continue
        n_reps_present += 1
        with p.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                by_variant.setdefault(row["variant"], []).append(
                    float(row["log2_fold_change"])
                )

    variants = []
    for variant, l2fcs in sorted(by_variant.items()):
        mean_l2fc = sum(l2fcs) / len(l2fcs)
        fraction_retained = 2 ** mean_l2fc  # fraction of wild-type PPIF expression remaining
        pct_change = (fraction_retained - 1) * 100
        variants.append({
            "variant": variant,
            "n_replicates_observed": len(l2fcs),
            "mean_log2_fold_change": round(mean_l2fc, 4),
            "fraction_expression_retained": round(fraction_retained, 4),
            "pct_expression_change": round(pct_change, 2),
        })

    pct_changes = [v["pct_expression_change"] for v in variants]
    retained = [v["fraction_expression_retained"] for v in variants]
    return {
        "n_replicate_files": n_reps_present,
        "n_distinct_edits": len(variants),
        "variants": variants,
        "pct_expression_change_min": round(min(pct_changes), 2) if pct_changes else None,
        "pct_expression_change_max": round(max(pct_changes), 2) if pct_changes else None,
        "fraction_expression_retained_min": round(min(retained), 4) if retained else None,
        "fraction_expression_retained_max": round(max(retained), 4) if retained else None,
        "paper_claim": "three edits disrupting the splice donor motif at the "
                       "first 5' splice site of PPIF, effects from -80% to "
                       "-52% for the single edits (Fig. 1e/f)",
    }


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sx = sum((a - mx) ** 2 for a in xs) ** 0.5
    sy = sum((b - my) ** 2 for b in ys) ** 0.5
    return cov / (sx * sy)


def score_lentimpra_vs_endogenous(run_dir: Path) -> dict:
    endo_path = run_dir / ENDOGENOUS_PROMOTER_FILE
    mpra_path = run_dir / LENTIMPRA_FILE
    endo_effect: dict[str, float] = {}
    with endo_path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            endo_effect[row["variant"]] = float(row["effect_size"])
    mpra_logfc: dict[str, float] = {}
    with mpra_path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mpra_logfc[row["variant"]] = float(row["logFC"])

    common = sorted(set(endo_effect) & set(mpra_logfc))
    x = [endo_effect[v] for v in common]
    y = [mpra_logfc[v] for v in common]
    r = pearson(x, y) if len(common) >= 2 else None
    return {
        "n_endogenous_variants": len(endo_effect),
        "n_lentimpra_variants": len(mpra_logfc),
        "n_matched_by_variant_id": len(common),
        "correlation_endogenous_vs_lentimpra": round(r, 4) if r is not None else None,
        "paper_claim": "variant effects were positively correlated between "
                       "the two assays (Pearson's r=0.54); systematic "
                       "differences reflect known limitations of MPRA "
                       "(genomic-context specificity)",
    }


def main() -> int:
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    summary_path = run_dir / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}

    summary["ppif_splice_site"] = score_splice_site(run_dir)
    print(f"PPIF splice site: {summary['ppif_splice_site']}")

    summary["lentimpra_vs_endogenous"] = score_lentimpra_vs_endogenous(run_dir)
    print(f"lentiMPRA vs endogenous: {summary['lentimpra_vs_endogenous']}")

    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
