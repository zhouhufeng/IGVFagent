#!/usr/bin/env python3
"""Score the real Martyn 2025 Variant-EFFECTS variant-effect tables IGVFagent
downloaded from the IGVF Portal (GRCh38 tabular-files, uniform-pipeline
output — not raw reads, not synthetic data).

Reads the three per-element TSVs and writes summary.json under the same
timestamped run directory, with per-element variant/significance counts and
effect-size ranges, plus the PPIF enhancer-to-promoter distance as an
independent cross-check against the paper's stated "~60.5 kb upstream"
enhancer position.
"""
import csv
import json
import math
import sys
from pathlib import Path

FDR_THRESH = -math.log10(0.05)

ELEMENTS = {
    "PPIF_promoter": "PPIF_promoter_GRCh38.tsv",
    "PPIF_enhancer": "PPIF_enhancer_GRCh38.tsv",
    "IL2RA_promoter": "IL2RA_promoter_GRCh38.tsv",
}


def score_one(path: Path) -> dict:
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    n = len(rows)
    effects = [float(r["effect_size"]) for r in rows]
    sig = [r for r in rows if float(r["fdr_nlog10"]) > FDR_THRESH]
    positions = [int(r["pos"]) for r in rows]
    return {
        "n_variants": n,
        "n_significant_fdr05": len(sig),
        "fraction_significant": round(len(sig) / n, 4),
        "effect_size_min": round(min(effects), 4),
        "effect_size_max": round(max(effects), 4),
        "effect_size_mean": round(sum(effects) / n, 4),
        "pos_min": min(positions),
        "pos_max": max(positions),
        "chrom": rows[0]["chr"],
    }


def main() -> int:
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    summary: dict = {"elements": {}}
    for name, fname in ELEMENTS.items():
        p = run_dir / fname
        if not p.is_file():
            print(f"missing {p}", file=sys.stderr)
            return 1
        summary["elements"][name] = score_one(p)
        print(f"{name}: {summary['elements'][name]}")

    ppif_promoter_tss = summary["elements"]["PPIF_promoter"]["pos_min"]
    ppif_enhancer_pos = summary["elements"]["PPIF_enhancer"]["pos_min"]
    dist = abs(ppif_promoter_tss - ppif_enhancer_pos)
    summary["ppif_enhancer_tss_distance_bp"] = dist
    print(f"PPIF enhancer-TSS distance: {dist:,} bp "
          f"(paper states ~60.5 kb upstream)")

    out = run_dir / "summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
