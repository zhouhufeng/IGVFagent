#!/usr/bin/env python3
"""Score the Tabula Sapiens reproduction against the paper.

Reads whatever `run.sh` produced and writes ``concordance_metrics.json``
into the overview run directory, which ``Benchmarks/concordance.py``
then checks against ``expected.json``.

Metrics come in two kinds:

  * **Exact counts** the paper states outright -- cells, donors,
    tissues, cell types, populations, the droplet/FACS split. These are
    asserted equal, because a reproduction that gets 1,136,217 cells has
    a bug, not a rounding difference.
  * **Boolean claims** -- the age span, the sex split, the age-group
    breakdown, the multi-organ donors, and the biology checks (FOXP3 in
    T cells, germ-cell TFs in spermatogenic cells). These collapse a
    stated sentence in the paper to a single pass/fail.

Whatever tier of `run.sh` ran, only the metrics it could produce
appear. A metadata-tier run therefore scores **16/20 partial**: the four
Figure 2 checks (TFs found in the atlas, the ubiquitous-TF tau bound,
FOXP3 in T cells, germ-cell TFs in spermatogenic cells) need the 57 GB
atlas and report as failures until `run.sh --full` has run. That is the
honest outcome -- they are genuinely scoreable, just not from metadata
alone, so they are asserted rather than excused.

Separately, four quantities that need the FULL atlas to be comparable at
all -- the 890/745 tau split, 48,114 senescent cells, 3,792 SAGs, 17
pathways -- carry ``confirmed: false`` in ``expected.json`` and render as
NOT SCORED with their provenance.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOCS = ROOT / "Docs" / "TabulaSapiens"
LABEL = "quake2026_tabula_sapiens"

# Paper Fig 2B / 2F claims, checked by peak cell type.
TCELL_TOKENS = ("t cell", "thymocyte", "lymphocyte")
GERM_TOKENS = ("sperm", "germ", "oocyte")
UBIQUITOUS_TFS = ["NFAT5", "NCOA1", "FOXJ3", "FOXK2", "ATF4", "JUN",
                  "FOS", "STAT1"]


def latest(suffix: str) -> "Path | None":
    dirs = sorted((p for p in DOCS.glob(f"2*_{LABEL}_{suffix}") if p.is_dir()),
                  key=lambda p: p.name, reverse=True)
    return dirs[0] if dirs else None


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
    ov = latest("overview")
    if ov is None:
        raise SystemExit(
            f"no overview run under {DOCS}/2*_{LABEL}_overview — run run.sh")
    m: "dict[str, object]" = {"benchmark": LABEL}

    # ── dataset-level counts (Figure 1) ───────────────────────────────
    s = read_json(ov / "overview_summary.json")
    m["cells_total"] = int(s["cells"])
    m["donors"] = int(s["donors"])
    m["tissues"] = int(s["tissues"])
    m["fine_cell_types"] = int(s["fine_cell_types"])
    m["populations"] = int(s["populations"])
    by_method = s.get("by_method") or {}
    m["cells_droplet"] = int(by_method.get("10X", 0))
    m["cells_facs"] = int(by_method.get("smartseq", 0))

    # ── demographic claims, as booleans ───────────────────────────────
    m["age_span_ok"] = bool(_f(s.get("age_min")) == 22
                            and _f(s.get("age_max")) == 74)
    sex = {k.lower(): v for k, v in (s.get("sex") or {}).items()}
    m["sex_split_ok"] = bool(sex.get("male") == 11 and sex.get("female") == 13)
    ag = s.get("age_groups") or {}
    m["age_groups_ok"] = bool(ag.get("<40") == 7 and ag.get("40-59") == 11
                              and ag.get(">=60") == 6)

    donors = read_tsv(ov / "donors.tsv")
    tissue_counts = sorted((int(d["tissues"]) for d in donors), reverse=True)
    m["max_tissues_per_donor"] = tissue_counts[0] if tissue_counts else 0
    # The paper's nine NEW donors include one each with 20, 18 and 15
    # tissues; the full atlas also carries v1 donors, so require presence
    # rather than an exact multiset.
    m["multiorgan_donors_ok"] = bool({20, 18, 15} <= set(tissue_counts))

    # ── droplet-subset counts, the tau denominator ────────────────────
    # `tabula status` recomputes these; recompute here so the benchmark
    # does not depend on parsing console output.
    try:
        sys.path.insert(0, str(ROOT / "Scripts"))
        import _ts_atlas as A  # type: ignore
        meta = A.cell_metadata()
        dr = meta[meta["method"] == "10X"]
        m["fine_cell_types_droplet"] = int(dr["cell_ontology_class"].nunique())
        m["broad_cell_types_droplet"] = int(dr["broad_cell_class"].nunique())
    except Exception as exc:                       # pragma: no cover
        print(f"(droplet-subset counts skipped: {exc})")

    # ── Human TF database ─────────────────────────────────────────────
    try:
        import humantfs_skill as H  # type: ignore
        m["htf_tf_count"] = len(H.tf_symbols())
    except Exception as exc:                       # pragma: no cover
        print(f"(TF database check skipped: {exc})")

    # ── Figure 2 (atlas tier) ─────────────────────────────────────────
    spec = latest("tf_specificity")
    if spec is not None:
        ts = read_json(spec / "tf_specificity_summary.json")
        m["tf_in_atlas"] = int(ts["n_tf"])
        m["tf_specific"] = int(ts["n_specific"])
        m["tf_non_specific"] = int(ts["n_non_specific"])
        m["tf_cell_types"] = int(ts["n_cell_types"])
        zero = set(ts.get("zero_expression_genes") or [])
        m["zero_expression_genes"] = sorted(zero)
        # The paper's sharpest single claim: exactly SHOX and ZBED1.
        m["zero_expression_exactly_shox_zbed1"] = bool(zero == {"SHOX", "ZBED1"})

        rows = read_tsv(spec / "tf_tau.tsv")
        tau_of = {r["gene"]: _f(r["tau"]) for r in rows if r["tau"] != ""}
        peak_of = {r["gene"]: (r["peak_cell_type"] or "").lower() for r in rows}
        present = [g for g in UBIQUITOUS_TFS if g in tau_of]
        m["ubiquitous_tf_checked"] = len(present)
        m["ubiquitous_tf_max_tau"] = (round(max(tau_of[g] for g in present), 4)
                                      if present else None)
        m["ubiquitous_tf_all_below_threshold"] = bool(
            present and all(tau_of[g] <= 0.85 for g in present))
        m["foxp3_peak_is_tcell"] = bool(
            any(t in peak_of.get("FOXP3", "") for t in TCELL_TOKENS))
        germ = [g for g in ("SOX30", "DMRT1", "SALL4", "TCFL5", "NKX1-2")
                if g in peak_of]
        m["germ_tf_peak_ok"] = bool(
            germ and any(any(t in peak_of[g] for t in GERM_TOKENS)
                         for g in germ))

    # ── Figure 3 (atlas tier) ─────────────────────────────────────────
    enr = latest("tf_enrichment")
    if enr is not None:
        es = read_json(enr / "tf_enrichment_summary.json")
        m["enrichment_terms"] = int(es["n_terms"])
        m["enrichment_groups"] = len(es.get("terms_by_function_group") or {})

    # ── Figure 4 (atlas tier) ─────────────────────────────────────────
    sen = latest("senescence")
    if sen is not None:
        ss = read_json(sen / "senescence_summary.json")
        m["senescent_cells"] = int(ss["senescent_cells"])
        m["senescent_fraction"] = round(_f(ss["senescent_fraction"]), 5)
        m["senescence_tissues"] = int(ss["n_tissues"])
        m["senescence_donors"] = int(ss["n_donors"])
        m["senescence_cells_dropped"] = int(ss["cells_dropped_by_balance"])
        top = [t.lower() for t in (ss.get("highest_burden_tissues") or [])[:3]]
        m["senescence_top_tissues"] = top
        m["senescence_top_tissues_ok"] = bool(
            {"eye", "bladder", "tongue"} == set(top))
        ages = ss.get("median_fraction_by_age_group") or {}
        young, mid, old = (_f(ages.get("<40")), _f(ages.get("40-59")),
                           _f(ages.get(">=60")))
        # The paper's claim is directional, not a magnitude: "a small
        # increase ... from medium and old donors as compared with young".
        m["senescence_rises_with_age"] = bool(
            young == young and mid == mid and old == old
            and mid > young and old > young)

    out = ov / "concordance_metrics.json"
    out.write_text(json.dumps(m, indent=2, sort_keys=True))
    print(f"Wrote {out}")
    for k, v in sorted(m.items()):
        print(f"  {k}: {v}")

    fig_dir = HERE / "figures"
    fig_dir.mkdir(exist_ok=True)
    (fig_dir / "concordance_metrics.json").write_text(
        json.dumps(m, indent=2, sort_keys=True))
    for src in (ov / "Plots").glob("*.png"):
        (fig_dir / src.name).write_bytes(src.read_bytes())
    for stage in ("tf_specificity", "tf_enrichment", "senescence"):
        d = latest(stage)
        if d is None:
            continue
        for src in (d / "Plots").glob("*.png"):
            (fig_dir / src.name).write_bytes(src.read_bytes())
    return 0


if __name__ == "__main__":
    sys.exit(main())
