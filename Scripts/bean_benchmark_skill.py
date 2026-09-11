#!/usr/bin/env python3
"""Reproduce the crispr-bean paper's results, and say plainly where we differ.

Ryu et al., "Joint genotypic and phenotypic outcome modeling improves base
editing variant effect quantification", Nature Genetics 56, 925-937 (2024).
https://doi.org/10.1038/s41588-024-01726-6

WHY THIS IS A SEPARATE SKILL FROM `bean`. `igvfagent bean` analyses an IGVF
screen. This replicates a published analysis, and the two need opposite
defaults: an analysis wants the best available method, a replication wants
the method the paper used, on the paper's own data, scored by the paper's
own metric. Mixing them is how a "replication" quietly becomes a different
experiment that happens to agree.

The first thing this skill established is that the IGVF-deposited Sherwood
screens are NOT the paper's screens:

    IGVFDS6464SOVZ (18loci_uptake)   8,192 guides, ~12 splice-control genes,
                                     no labelled non-targeting controls
    paper, LDL-C GWAS library        3,455 guides, 6 splice-control genes
                                     (LDLR, MYLIP, ACAT2, SREBF2, HNF4A, LSS),
                                     100 non-targeting negative controls

Only LDLR and HNF4A overlap. Scoring the paper's numbers on the IGVF data
would compare different experiments and agree or disagree for the wrong
reason, so this skill works from the paper's own deposit at Zenodo
(10.5281/zenodo.10139794), which publishes the BEAN screen objects directly
-- the exact input `bean run` consumes.

PUBLISHED VALUES ARE DATA, NOT PROSE. Every number the paper reports lives
in PAPER_CLAIMS below with its source (figure, page or sentence). A
replication that quotes its target from memory is not checkable; this way a
mismatch is a diff, and adding a claim is adding a row.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base_editing_screen as bes                               # noqa: E402
from _stats import _percentile                                   # noqa: E402

logger = logging.getLogger("bean_benchmark")

ZENODO_RECORD = "10139794"
ZENODO_DOI = "10.5281/zenodo.10139794"
PAPER = ("Ryu et al., Nat Genet 56:925-937 (2024), "
          "doi:10.1038/s41588-024-01726-6")

DATA_DIR = Path(os.environ.get(
    "IGVF_BEAN_BENCH_DIR",
    str(Path(os.environ.get("IGVF_DATA_ROOT", "/workspace"))
        / "Data" / "Benchmarks" / "BEANpaper")))
OUT_DIR = Path(os.environ.get("IGVF_DATA_ROOT", "/workspace")) / "Docs" / "BEANbenchmark"

# The paper's two screens, by the file name it deposited them under.
SCREENS = {
    "ldlvar": {
        "file": "bean_count_LDLvar_annotated.h5ad",
        "label": "LDL-C GWAS library",
        "bytes": 394_500_000,
    },
    "ldlrcds": {
        "file": "bean_count_LDLRCDS_annotated_0.1_0.3.h5ad",
        "label": "LDLR CDS tiling library",
        "bytes": 555_900_000,
    },
}

# ── the paper's own numbers ────────────────────────────────────────────────
# `tol` is how close a reproduction has to be to count as agreeing. It is NOT
# a statistical statement -- it is a judgement about what the paper's own
# precision supports. A value quoted to two decimals cannot be reproduced to
# three, and a count of significant variants depends on a threshold applied
# to a stochastic fit, so it moves between runs of BEAN itself.
PAPER_CLAIMS = [
    {"id": "ldlvar.n_guides", "screen": "ldlvar", "kind": "count",
     "what": "gRNA species in the library", "value": 3455, "tol": 0,
     "source": "p.926: 'totaling 3,455 gRNA species'"},
    {"id": "ldlvar.n_variants", "screen": "ldlvar", "kind": "count",
     "what": "LDL-C associated variants targeted", "value": 583, "tol": 0,
     "source": "p.926: 'targets 583 variants associated with LDL-C levels'"},
    {"id": "ldlvar.n_negctrl", "screen": "ldlvar", "kind": "count",
     "what": "non-targeting negative control gRNAs", "value": 100, "tol": 0,
     "source": "p.926: '100 nontargeting negative control gRNA species'"},
    {"id": "ldlvar.replicate_rho", "screen": "ldlvar", "kind": "rho",
     "what": "median Spearman rho of gRNA counts across replicates",
     "value": 0.84, "tol": 0.05,
     "source": "p.928: 'median Spearman rho = 0.84 for the LDL-C GWAS library'"},
    {"id": "ldlvar.n_significant", "screen": "ldlvar", "kind": "count",
     "what": "variants with 95% CI excluding 0", "value": 54, "tol": 8,
     "source": "p.929: 'BEAN identified 54 variants that significantly alter "
                "LDL-C uptake'"},
    {"id": "ldlvar.auprc_bean", "screen": "ldlvar", "kind": "auprc",
     "what": "AUPRC, LDLR+MYLIP splice variants vs negative controls (BEAN)",
     "value": 0.90, "tol": 0.05,
     "source": "p.929, Fig.3b: 'mean AUPRC = 0.90 across 15 two-replicate "
                "subsamples'"},
    {"id": "ldlvar.auprc_bean_reporter", "screen": "ldlvar", "kind": "auprc",
     "what": "AUPRC, BEAN-Reporter (no accessibility)", "value": 0.87,
     "tol": 0.05, "source": "p.929: 'BEAN-Reporter (mean AUPRC = 0.87)'"},
    {"id": "ldlvar.auprc_bean_uniform", "screen": "ldlvar", "kind": "auprc",
     "what": "AUPRC, BEAN-Uniform (no reporter, uniform editing)",
     "value": 0.85, "tol": 0.05,
     "source": "p.929: 'BEAN-Uniform (mean AUPRC = 0.85)'"},
    {"id": "ldlvar.mean_edit_fraction", "screen": "ldlvar", "kind": "frac",
     "what": "average per-gRNA variant edit fraction", "value": 0.340,
     "tol": 0.05, "source": "p.929: 'an average edit fraction of 34.0%'"},
    {"id": "ldlvar.median_max_edit", "screen": "ldlvar", "kind": "frac",
     "what": "median maximal editing across the 5 gRNAs per variant",
     "value": 0.604, "tol": 0.05,
     "source": "p.929: 'median maximal editing of 60.4%'"},

    {"id": "ldlrcds.n_guides", "screen": "ldlrcds", "kind": "count",
     "what": "gRNA species in the tiling library", "value": 7500, "tol": 0,
     "source": "p.926: 'a total of 7,500 gRNA species'"},
    {"id": "ldlrcds.n_negctrl", "screen": "ldlrcds", "kind": "count",
     "what": "non-targeting negative control gRNAs", "value": 150, "tol": 0,
     "source": "p.926: '150 nontargeting negative control gRNA species'"},
    {"id": "ldlrcds.replicate_rho", "screen": "ldlrcds", "kind": "rho",
     "what": "median Spearman rho across replicates", "value": 0.88,
     "tol": 0.05, "source": "p.928: 'rho = 0.88 for the LDLR tiling library'"},
    {"id": "ldlrcds.n_variants", "screen": "ldlrcds", "kind": "count",
     "what": "distinct variants assessed", "value": 2182, "tol": 0,
     "source": "p.931: 'A total of 2,182 distinct variants were assessed'"},
    {"id": "ldlrcds.n_missense", "screen": "ldlrcds", "kind": "count",
     "what": "missense coding variants among them", "value": 874, "tol": 0,
     "source": "p.931: 'including 874 missense coding variants'"},
    {"id": "ldlrcds.n_significant", "screen": "ldlrcds", "kind": "count",
     "what": "variants with z < -1.96", "value": 145, "tol": 15,
     "source": "p.931: 'BEAN assigned significant z scores (<-1.96 ...) to "
                "145 variants'"},
    {"id": "ldlrcds.n_decreasing", "screen": "ldlrcds", "kind": "count",
     "what": "of those, variants decreasing LDL-C uptake", "value": 131,
     "tol": 15, "source": "p.931: '131 of which decreased LDL-C uptake'"},
    {"id": "ldlrcds.auprc_bean", "screen": "ldlrcds", "kind": "auprc",
     "what": "AUPRC, variant classification (BEAN)", "value": 0.88,
     "tol": 0.05, "source": "p.931: 'a high AUPRC of 0.88'"},
]

# Claims this skill cannot test, recorded so the gap is explicit rather than
# quietly missing from a report.
NOT_TESTABLE = [
    ("ldlrcds.ukb_spearman", "Spearman 0.40 / Pearson 0.45 vs UK Biobank "
      "patient LDL-C", "UKB is controlled access (application required)"),
    ("ldlrcds.fuse_ukb", "BEAN-FUSE Spearman 0.50 / Pearson 0.51 vs UKB",
      "needs both UKB access and the FUSE model"),
    ("ldlvar.individual_grna", "Spearman of BEAN effect vs individually "
      "tested gRNA LFC (26 gRNAs, 6 replicates)",
      "individual validation experiments are not in the deposit"),
]


def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


# ── metrics ────────────────────────────────────────────────────────────────

def spearman(a: "list[float]", b: "list[float]") -> float:
    """Rank correlation, ties averaged. Written out rather than imported so
    this skill keeps `bean`'s zero-heavy-dependency property."""
    if len(a) != len(b) or len(a) < 3:
        return float("nan")

    def ranks(xs):
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        r = [0.0] * len(xs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = ranks(a), ranks(b)
    n = len(ra)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((y - mb) ** 2 for y in rb))
    return num / (da * db) if da and db else float("nan")


def auprc(scores: "list[float]", labels: "list[int]") -> float:
    """Area under the precision-recall curve, by the trapezoid rule.

    The paper scores positive-control splice variants against non-targeting
    negative controls, so the positive class is rare and AUPRC is the right
    summary -- AUROC would look flattering on the same data.

    Higher score = more positive. Interpolation is the standard
    step-and-trapezoid over recall; ties are handled by advancing through the
    whole tied block before recording a point, or a run of equal scores would
    manufacture a staircase that inflates the area.
    """
    if not scores or len(scores) != len(labels):
        return float("nan")
    P = sum(1 for l in labels if l)
    if P == 0 or P == len(labels):
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = fp = 0
    prev_recall = 0.0
    area = 0.0
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for k in range(i, j + 1):
            if labels[order[k]]:
                tp += 1
            else:
                fp += 1
        recall = tp / P
        precision = tp / (tp + fp)
        area += (recall - prev_recall) * precision
        prev_recall = recall
        i = j + 1
    return area


# ── the paper's data ───────────────────────────────────────────────────────

def fetch(screen: str, *, force: bool = False) -> dict:
    """Download one screen object from the paper's Zenodo deposit."""
    spec = SCREENS[screen]
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    dest = DATA_DIR / spec["file"]
    if dest.exists() and not force:
        return {"screen": screen, "path": str(dest), "bytes": dest.stat().st_size,
                 "already_present": True}
    url = (f"https://zenodo.org/records/{ZENODO_RECORD}/files/"
           f"{spec['file']}?download=1")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=3600) as r, open(tmp, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.rename(dest)
    return {"screen": screen, "path": str(dest), "bytes": dest.stat().st_size,
             "already_present": False, "doi": ZENODO_DOI}


def describe(screen: str) -> dict:
    """What the deposited object actually contains.

    Run before anything is scored. A replication that starts by computing a
    metric cannot tell "the method disagrees" from "this is not the data the
    paper described", and those need completely different responses.
    """
    spec = SCREENS[screen]
    path = DATA_DIR / spec["file"]
    if not path.exists():
        return {"error": f"{path} not present — run `fetch {screen}` first"}
    try:
        import anndata
    except ImportError:
        return {"error": "anndata is not installed in this environment; "
                          "the BEAN venv has it (see Deploy/install-bean.sh)"}
    ad = anndata.read_h5ad(path)
    guides = ad.shape[0] if ad.shape else 0
    out = {"screen": screen, "label": spec["label"], "path": str(path),
            "n_guides": int(guides), "n_samples": int(ad.shape[1]),
            "layers": sorted(ad.layers.keys()),
            "guide_columns": sorted(ad.obs.columns.astype(str))[:40],
            "sample_columns": sorted(ad.var.columns.astype(str))[:40],
            "uns_keys": sorted(map(str, ad.uns.keys()))[:20]}
    # Guide classes, however this object spells them.
    for col in ("type", "group", "target_group", "Group", "class"):
        if col in ad.obs.columns:
            vc = ad.obs[col].value_counts().to_dict()
            out["guide_classes"] = {str(k): int(v) for k, v in vc.items()}
            out["guide_class_column"] = col
            break
    for col in ("target", "Target", "target_id"):
        if col in ad.obs.columns:
            out["n_distinct_targets"] = int(ad.obs[col].nunique())
            out["target_column"] = col
            break
    if "replicate" in ad.var.columns:
        out["replicates"] = sorted(map(str, ad.var["replicate"].unique()))
    if "condition" in ad.var.columns:
        out["conditions"] = sorted(map(str, ad.var["condition"].unique()))
    return out


def compare(measured: dict) -> "list[dict]":
    """Line the measured values up against PAPER_CLAIMS."""
    rows = []
    for c in PAPER_CLAIMS:
        got = measured.get(c["id"])
        if got is None:
            rows.append({**c, "measured": None, "status": "not measured"})
            continue
        delta = got - c["value"]
        ok = abs(delta) <= c["tol"]
        rows.append({**c, "measured": got, "delta": delta,
                      "status": "agrees" if ok else "DIFFERS"})
    return rows


def render(rows: "list[dict]", measured: dict) -> str:
    w = max((len(r["what"]) for r in rows), default = 20)
    out = [f"Benchmark against {PAPER}", f"  data: {ZENODO_DOI}", ""]
    for r in rows:
        got = r["measured"]
        gs = "—" if got is None else (f"{got:,}" if r["kind"] == "count"
                                       else f"{got:.3f}")
        vs = (f"{r['value']:,}" if r["kind"] == "count"
              else f"{r['value']:.3f}")
        mark = {"agrees": "ok  ", "DIFFERS": "DIFF", "not measured": "--  "}[r["status"]]
        out.append(f"  {mark} {r['what']:<{w}}  paper {vs:>8}   ours {gs:>8}")
    n_ok = sum(1 for r in rows if r["status"] == "agrees")
    n_diff = sum(1 for r in rows if r["status"] == "DIFFERS")
    n_na = sum(1 for r in rows if r["status"] == "not measured")
    out += ["", f"  {n_ok} agree, {n_diff} differ, {n_na} not yet measured",
             "", "  Claims this benchmark cannot test:"]
    for cid, what, why in NOT_TESTABLE:
        out.append(f"    - {what}  ({why})")
    return "\n".join(out)


# ── commands ───────────────────────────────────────────────────────────────

def cmd_fetch(args) -> int:
    setup_logging()
    names = list(SCREENS) if args.screen == "all" else [args.screen]
    for n in names:
        out = fetch(n, force=args.force)
        state = "already present" if out.get("already_present") else "downloaded"
        print(f"  {n:9} {state}: {out['bytes'] / 1e6:,.1f} MB  {out['path']}")
    return 0


def cmd_describe(args) -> int:
    setup_logging()
    out = describe(args.screen)
    if "error" in out:
        print(f"  {out['error']}")
        return 2
    print(f"{out['label']}  ({out['screen']})")
    print(f"  {out['n_guides']:,} guides x {out['n_samples']} samples")
    print(f"  layers:  {out['layers']}")
    if out.get("guide_classes"):
        print(f"  classes ({out['guide_class_column']}): {out['guide_classes']}")
    if out.get("n_distinct_targets"):
        print(f"  distinct targets: {out['n_distinct_targets']:,}")
    for k in ("replicates", "conditions"):
        if out.get(k):
            print(f"  {k}: {out[k]}")
    print(f"  guide columns:  {out['guide_columns']}")
    print(f"  sample columns: {out['sample_columns']}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"describe_{args.screen}.json").write_text(json.dumps(out, indent=2))
    print(f"\n  -> {OUT_DIR / f'describe_{args.screen}.json'}")
    return 0


def cmd_report(args) -> int:
    setup_logging()
    measured = {}
    p = OUT_DIR / "measured.json"
    if p.exists():
        measured = json.loads(p.read_text())
    rows = compare(measured)
    print(render(rows, measured))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "benchmark_report.json").write_text(
        json.dumps({"paper": PAPER, "doi": ZENODO_DOI, "rows": rows}, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent bean-benchmark",
        description="Reproduce the crispr-bean paper (Ryu et al. 2024) on its "
                     "own deposited data, and report where we differ.")
    sub = p.add_subparsers(dest="command", required=True)

    f = sub.add_parser("fetch", help="Download the paper's screen objects "
                                      "from Zenodo.")
    f.add_argument("screen", choices=list(SCREENS) + ["all"], default="all",
                   nargs="?")
    f.add_argument("--force", action="store_true")

    d = sub.add_parser("describe", help="What a deposited screen actually "
                                         "contains — run this first.")
    d.add_argument("screen", choices=list(SCREENS))

    sub.add_parser("report", help="Measured vs published, claim by claim.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"fetch": cmd_fetch, "describe": cmd_describe,
            "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
