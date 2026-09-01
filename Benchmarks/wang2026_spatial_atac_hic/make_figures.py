#!/usr/bin/env python3
"""Score the Spatial-ATAC-Hi-C chain against the planted ground truth.

Reads the artefacts written by ``run.sh`` and emits
``concordance_metrics.json`` into the QC run directory, plus a summary
figure. ``Benchmarks/concordance.py`` then checks those metrics against
``expected.json``.

Every metric here is a recovery statistic: what the pipeline reported
divided by, or compared against, what ``prep_input.py`` planted.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "Benchmarks" / "_data" / "wang2026_spatial_atac_hic"
DOCS = ROOT / "Docs" / "SpatialATACHiC"
LABEL = "wang2026_spatial_atac_hic"


def latest(suffix: str) -> Path:
    """Most recent run directory for one stage of the chain."""
    dirs = sorted((p for p in DOCS.glob(f"2*_{LABEL}_{suffix}") if p.is_dir()),
                  key=lambda p: p.name, reverse=True)
    if not dirs:
        raise SystemExit(f"no run directory Docs/SpatialATACHiC/*_{LABEL}_{suffix}")
    return dirs[0]


def read_json(p: Path):
    return json.loads(p.read_text())


def read_tsv(p: Path) -> "list[dict]":
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _f(v, default=float("nan")) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main() -> int:
    truth = read_json(DATA / "truth.json")
    m: "dict[str, object]" = {"benchmark": LABEL, "truth_source": "prep_input.py"}

    # ── 1. demultiplex ────────────────────────────────────────────────
    demux = read_json(latest("demux") / "demux_summary.json")
    m["demux_assigned_fraction"] = round(_f(demux["assigned_fraction"]), 6)
    m["pixels_recovered"] = int(demux["pixels_with_data"])
    m["pixels_expected"] = int(truth["n_pixels"])

    # ── 2. contact QC vs the paper's reported bands ───────────────────
    qc_dir = latest("qc")
    qc = read_json(qc_dir / "qc_summary.json")
    m["median_cis_fraction"] = round(_f(qc["median_cis_fraction"]), 6)
    m["median_long_range_ratio"] = round(_f(qc["median_long_range_ratio"]), 6)
    m["median_total_contacts"] = _f(qc["median_total_contacts"])
    # Score against what the generator actually produced (measured
    # independently in prep_input.py), not against its aspirational
    # target — otherwise a miscalibrated generator is charged to the skill.
    real = truth["realised"]
    m["cis_fraction_abs_error"] = round(
        abs(_f(qc["median_cis_fraction"]) - _f(real["cis_fraction"])), 6)
    m["long_range_ratio_abs_error"] = round(
        abs(_f(qc["median_long_range_ratio"]) - _f(real["long_range_ratio"])), 6)

    # ── 3. TSS enrichment ─────────────────────────────────────────────
    if "tss_enrichment" in qc:
        planted = _f(real["tss_enrichment"])
        got = _f(qc["tss_enrichment"])
        m["tss_enrichment"] = round(got, 4)
        m["tss_enrichment_planted"] = round(planted, 4)
        m["tss_enrichment_recovery"] = round(got / planted, 4) if planted else None

    # ── 4. copy number: pseudobulk gain, and clone separation ─────────
    cnv_dir = latest("cnv")
    rows = read_tsv(cnv_dir / "cnv_pseudobulk.tsv")
    g = truth["gain"]
    inside = [_f(r["cn_ratio"]) for r in rows
              if r["chrom"] == g["chrom"] and r["cn_ratio"] != ""
              and g["start"] <= int(r["start"]) < g["end"]]
    outside = [_f(r["cn_ratio"]) for r in rows
               if r["cn_ratio"] != "" and
               not (r["chrom"] == g["chrom"]
                    and g["start"] <= int(r["start"]) < g["end"])]
    mean_in = sum(inside) / len(inside) if inside else float("nan")
    mean_out = sum(outside) / len(outside) if outside else float("nan")
    m["cnv_background_ratio"] = round(mean_out, 4)
    observed = mean_in / mean_out if mean_out else float("nan")
    m["cnv_gain_fold_observed"] = round(observed, 4) if observed == observed else None
    m["cnv_gain_fold_planted"] = round(_f(real["gain_fold"]), 4)
    m["cnv_gain_fold_recovery"] = (
        round(observed / _f(real["gain_fold"]), 4)
        if observed == observed and _f(real["gain_fold"]) else None)

    # Per-pixel: does the gain localise to the planted clone?
    pp = cnv_dir / "cnv_per_pixel.tsv"
    if pp.is_file():
        lines = pp.read_text().splitlines()
        hdr = lines[0].split("\t")
        cols = [i for i, h in enumerate(hdr)
                if h.startswith(f"{g['chrom']}:")
                and g["start"] <= int(h.split(":")[1]) < g["end"]]
        by_clone: "dict[str, list[float]]" = {}
        for ln in lines[1:]:
            f = ln.split("\t")
            clone = truth["clone_of_pixel"].get(f[0])
            if clone is None:
                continue
            vals = [_f(f[i]) for i in cols if i < len(f) and f[i] != ""]
            vals = [v for v in vals if v == v]
            if vals:
                by_clone.setdefault(clone, []).append(sum(vals) / len(vals))
        carrier = g["clone"]
        other = next((k for k in by_clone if k != carrier), None)
        if carrier in by_clone and other:
            a = sum(by_clone[carrier]) / len(by_clone[carrier])
            b = sum(by_clone[other]) / len(by_clone[other])
            m["cnv_clone_separation"] = round(a / b, 4) if b else None
            m["cnv_carrier_clone_mean"] = round(a, 4)
            m["cnv_other_clone_mean"] = round(b, 4)

    # ── 5. loops: the planted one specific, the null one not ──────────
    loops_dir = latest("loops")
    lrows = read_tsv(loops_dir / "cluster_specific_loops.tsv")
    lsum = read_json(loops_dir / "loops_summary.json")
    planted_name = truth["loop"]["name"]
    hit = next((r for r in lrows if r["name"] == planted_name), None)
    null = next((r for r in lrows if r["name"] == truth["null_loop"]), None)
    if hit:
        m["planted_loop_p_value"] = _f(hit["p_value"])
        m["planted_loop_neglog10_p"] = (
            round(-math.log10(max(_f(hit["p_value"]), 1e-300)), 3)
            if hit["p_value"] else None)
        m["planted_loop_called_specific"] = hit["specific"] == "yes"
        m["planted_loop_top_cluster"] = hit["top_cluster"]
        m["planted_loop_clone_correct"] = hit["top_cluster"] == truth["loop"]["clone"]
    if null:
        m["null_loop_called_specific"] = null["specific"] == "yes"
    m["n_cluster_specific_loops"] = int(lsum.get("n_cluster_specific_loops", -1))

    # ── 6. compartments: do the called A/B blocks match the plant? ────
    comp_dir = latest("compartment")
    crows = read_tsv(comp_dir / "compartments.tsv")
    block = int(truth["compartment"]["block"])
    cchrom = truth["compartment"]["chrom"]
    agree = total = 0
    for r in crows:
        if r["chrom"] != cchrom or not r["compartment"]:
            continue
        planted_a = (int(r["start"]) // block) % 2 == 0
        total += 1
        if (r["compartment"] == "A") == planted_a:
            agree += 1
    if total:
        acc = agree / total
        # PC1's sign is only fixed up to orientation; report the better of
        # the two labellings so a globally flipped track is not scored as
        # a failure to find structure.
        m["compartment_accuracy"] = round(max(acc, 1 - acc), 4)
        m["compartment_bins_scored"] = total

    # ── write ─────────────────────────────────────────────────────────
    out = qc_dir / "concordance_metrics.json"
    out.write_text(json.dumps(m, indent=2, sort_keys=True))
    print(f"Wrote {out}")
    for k, v in sorted(m.items()):
        print(f"  {k}: {v}")

    _figure(m, qc_dir, truth)
    return 0


def _figure(m: dict, out_dir: Path, truth: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"(figure skipped: {exc})")
        return

    panels = [
        ("CN gain fold", m.get("cnv_gain_fold_observed"),
         m.get("cnv_gain_fold_planted")),
        ("CN clone sep.", m.get("cnv_clone_separation"),
         m.get("cnv_gain_fold_planted")),
        ("TSS enrichment", m.get("tss_enrichment"),
         m.get("tss_enrichment_planted")),
        ("cis fraction", m.get("median_cis_fraction"),
         _f(truth["realised"]["cis_fraction"])),
        ("long-range ratio", m.get("median_long_range_ratio"),
         _f(truth["realised"]["long_range_ratio"])),
        ("compartment acc.", m.get("compartment_accuracy"), 1.0),
    ]
    panels = [(n, o, e) for n, o, e in panels
              if isinstance(o, (int, float)) and o == o]
    if not panels:
        return

    fig, ax = plt.subplots(figsize=(8.2, 4.0), dpi=200)
    x = range(len(panels))
    ax.bar([i - 0.19 for i in x], [p[2] for p in panels], width=0.38,
           label="planted", color="#9aa5b1")
    ax.bar([i + 0.19 for i in x], [p[1] for p in panels], width=0.38,
           label="recovered", color="#2f6f9f")
    ax.set_xticks(list(x))
    ax.set_xticklabels([p[0] for p in panels], rotation=20, ha="right",
                       fontsize=9)
    ax.set_yscale("log")
    ax.set_ylabel("value (log scale)")
    ax.set_title("Spatial-ATAC-Hi-C: planted vs recovered "
                 "(Wang 2026 benchmark)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    png = out_dir / "Plots" / "benchmark_recovery.png"
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png)
    fig.savefig(png.with_suffix(".svg"))
    plt.close(fig)
    print(f"Figure: {png}")

    # Also drop a copy next to the benchmark README, which embeds it. The
    # run directory is timestamped and gitignored; this copy is the one
    # that survives in the repo.
    repo_fig = HERE / "figures"
    repo_fig.mkdir(exist_ok=True)
    for ext in (".png", ".svg"):
        (repo_fig / f"benchmark_recovery{ext}").write_bytes(
            png.with_suffix(ext).read_bytes())
    (repo_fig / "concordance_metrics.json").write_text(
        json.dumps(m, indent=2, sort_keys=True))


if __name__ == "__main__":
    sys.exit(main())
