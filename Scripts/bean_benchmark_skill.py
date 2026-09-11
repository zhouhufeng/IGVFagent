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


# Findings from the first run against the deposited LDLvar object, recorded
# here so they are not re-derived. These are questions for the authors, not
# defects in either side:
#
#   * The deposit carries EIGHT replicate labels (rep5, rep9-rep15) while the
#     paper says five. Fig.3b's "15 two-replicate subsamples among the five
#     replicates" is also arithmetically odd: C(5,2) = 10, and 15 = C(6,2).
#     Which replicates were used changes every metric below, so this matters
#     before any AUPRC comparison is called a match or a miss.
#
#   * Bins are top/high/low/bot plus bulk. That is CONSISTENT with the
#     paper's "four populations per replicate" -- four sorted, plus an
#     unsorted reference -- and is not a discrepancy.
#
#   * Mean per-gRNA edit fraction reproduces (0.361 overall; 0.339 in the top
#     bin, 0.345 in bot, vs the paper's 34.0%). The MEDIAN MAXIMAL editing
#     per variant does not: 0.496 measured vs 0.604 published. Tested and
#     ruled out: the choice between obs['edit_rate'] and the per-sample
#     layers['edit_rate'], and which bin the layer is averaged over -- all
#     give 0.496-0.513. The remaining hypothesis is that the paper's
#     "maximal editing" uses a VARIANT-SPECIFIC rate (the intended edit at
#     the target position, from the allele-level `edits` layer) rather than
#     the aggregate per-guide rate used here.
#
#   * Both deposits are FILTERED, and the pattern is consistent: LDLRCDS is
#     literally named ..._0.1_0.3.h5ad after its thresholds. Guide-level
#     counts that survive filtering match the paper EXACTLY (LDLRCDS: 7,500
#     guides, 150 non-targeting controls), while variant-level counts fall
#     short (1,894 vs 2,182 assessed; 863 vs 874 missense; LDLvar 570 vs
#     583). So the deposits are the post-QC objects and the paper's counts
#     are pre-filter -- which is worth confirming with the authors, because
#     it is the difference between "we lost variants" and "they were never
#     there".
#
#   * The deposited object is the filtered/annotated version, so it holds
#     3,451 guides / 99 negative controls / 570 variants against the paper's
#     3,455 / 100 / 583. The claim tolerances are deliberately left at 0:
#     widening them to absorb filtering would hide exactly the kind of drift
#     this benchmark exists to detect.


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


# The three model variants the paper compares. Everything IGVF's own screens
# cannot support -- the reporter, the guide barcode, accessibility -- IS in
# this deposit (layers X_bcmatch, edit_rate, edits; obs column Reporter), so
# all three run here. A report that says otherwise is reasoning from the IGVF
# data, not from this.
MODELS = {
    "bean": {
        "label": "BEAN (MixtureNormal + accessibility)",
        "flags": ["--scale-by-acc"],
        "claim": "auprc_bean",
    },
    "reporter": {
        "label": "BEAN-Reporter (reporter editing, no accessibility)",
        "flags": [],
        "claim": "auprc_bean_reporter",
    },
    "uniform": {
        "label": "BEAN-Uniform (no reporter, uniform editing)",
        "flags": ["--uniform-edit", "--ignore-bcmatch"],
        "claim": "auprc_bean_uniform",
    },
}


def bean_exe() -> "Optional[str]":
    ok, detail = bes.bean_available()
    return detail.split()[0] if ok else None


def run_model(screen: str, model: str, *, n_iter: int = 0,
               timeout: int = 7200) -> dict:
    """Fit one of the paper's three model variants on the deposited object.

    The column names come from the deposit itself rather than from BEAN's
    defaults: replicate is `rep`, condition is `bin`, and the unsorted
    reference bin is `bulk`. BEAN defaults --control-condition to "bulk",
    which happens to be right here, but it is passed explicitly so a deposit
    that spelled it differently would fail loudly instead of silently
    normalising against a sorted bin.
    """
    spec = SCREENS[screen]
    path = DATA_DIR / spec["file"]
    if not path.exists():
        return {"error": f"{path} not present — run `fetch {screen}` first"}
    exe = bean_exe()
    if not exe:
        return {"error": "the real `bean` binary is not installed here; "
                          "see Deploy/install-bean.sh"}
    m = MODELS[model]
    outdir = OUT_DIR / f"{screen}_{model}"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "run", "sorting", "variant", str(path),
            "--replicate-col", "rep",
            "--condition-col", "bin",
            "--target-col", "target",
            "--control-condition", "bulk",
            "--sorting-bin-lower-quantile-col", "lower_quantile",
            "--sorting-bin-upper-quantile-col", "upper_quantile",
            "--outdir", str(outdir)] + m["flags"]
    if n_iter:
        cmd += ["--n-iter", str(n_iter)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    out = {"screen": screen, "model": model, "label": m["label"],
            "exit_code": r.returncode, "cmd": " ".join(cmd),
            "outdir": str(outdir),
            "stderr": (r.stderr or "").strip()[-800:]}
    if r.returncode == 0:
        hits = sorted(outdir.glob("**/bean_element_result.*.csv"))
        out["result_csv"] = str(hits[0]) if hits else None
    return out


def score_auprc(screen: str, result_csv: str) -> dict:
    """The paper's metric: positive-control splice variants vs non-targeting.

    Fig.3b classifies LDLR and MYLIP splicing variants against the negative
    controls. The score is -z, because a splice-disrupting variant REDUCES
    LDL uptake: ranking on raw z would put the true positives last and score
    the metric upside down.
    """
    import csv as _csv
    import anndata
    path = DATA_DIR / SCREENS[screen]["file"]
    ad = anndata.read_h5ad(path)
    grp = "target_group" if "target_group" in ad.obs.columns else "Group"
    tcol = "target" if "target" in ad.obs.columns else ad.obs.columns[0]
    # target -> class, taken from the object rather than from the name.
    cls = {}
    for t, g in zip(ad.obs[tcol].astype(str), ad.obs[grp].astype(str)):
        cls.setdefault(t, g)
    rows = []
    with open(result_csv, newline="") as fh:
        for row in _csv.DictReader(fh):
            t = row.get("target") or row.get("")
            z = row.get("mu_z") or row.get("z")
            if t is None or z in (None, ""):
                continue
            try:
                rows.append((t, float(z)))
            except ValueError:
                continue
    scores, labels, pos_genes = [], [], ("LDLR", "MYLIP")
    for t, z in rows:
        c = cls.get(t, "")
        is_pos = c.lower().startswith("pos") and any(
            t.upper().startswith(g) for g in pos_genes)
        is_neg = "neg" in c.lower() or c.lower().endswith("control")
        if not (is_pos or is_neg):
            continue
        scores.append(-z)          # splice disruption lowers uptake
        labels.append(1 if is_pos else 0)
    return {"auprc": auprc(scores, labels), "n_positive": sum(labels),
             "n_negative": len(labels) - sum(labels), "n_scored": len(rows)}


def measure(screen: str) -> dict:
    """Everything computable from the deposited object alone, no model fit.

    Kept separate from the BEAN run because these are properties of the DATA
    -- library composition, replicate agreement, editing rates -- and a
    disagreement here means the deposit differs from the paper's description,
    which is a different finding from the model disagreeing.
    """
    spec = SCREENS[screen]
    path = DATA_DIR / spec["file"]
    if not path.exists():
        return {"error": f"{path} not present — run `fetch {screen}` first"}
    import anndata
    import numpy as np
    ad = anndata.read_h5ad(path)
    out = {f"{screen}.n_guides": int(ad.shape[0])}

    # The two deposited screens do NOT share a schema, and a detector that
    # assumed one silently reports zero for the other:
    #
    #   LDLvar   target_group = Variant / PosCtrl / NegCtrl
    #   LDLRCDS  Group        = exon numbers, UTRs, DNase HS regions,
    #                          'PosCtrl', 'ABE control', 'CBE control'
    #
    # In LDLRCDS the non-targeting guides are split by editor: 'ABE control'
    # (75) + 'CBE control' (75) = 150, exactly the count the paper states,
    # which is what identifies them as the negative controls.
    grp = next((c for c in ("target_group", "Group", "type", "group")
                 if c in ad.obs.columns), None)
    if grp:
        vc = {str(k): int(v) for k, v in ad.obs[grp].value_counts().items()}
        neg = sum(v for k, v in vc.items()
                   if "negctrl" in k.lower().replace(" ", "")
                   or k.lower().endswith("control"))
        out[f"{screen}.n_negctrl"] = neg
        out[f"{screen}._classes"] = vc
        out[f"{screen}._negctrl_groups"] = sorted(
            k for k in vc
            if "negctrl" in k.lower().replace(" ", "")
            or k.lower().endswith("control"))

    # Distinct VARIANTS, which is not the same as distinct targets: the
    # positive controls and the non-targeting guides also carry a target.
    tv = next((c for c in ("target_variant", "target_allEdited", "target")
                if c in ad.obs.columns), None)
    if tv and grp:
        cls = ad.obs[grp].astype(str)
        if screen == "ldlrcds":
            # Every group that is not a control is a targeted region, so the
            # variants assessed are the targets of the non-control guides.
            neg_groups = set(out.get(f"{screen}._negctrl_groups", []))
            mask = ~cls.isin(neg_groups | {"PosCtrl"})
        else:
            mask = cls.str.lower().str.startswith("variant")
        out[f"{screen}.n_variants"] = int(ad.obs.loc[mask, tv].nunique())

    # Consequence, for the tiling screen's missense count. `severity` is the
    # only column carrying it; the mapping is read off the data rather than
    # assumed, and recorded, because calling the wrong level "missense" would
    # produce a plausible number for the wrong reason.
    if "severity" in ad.obs.columns and tv:
        sv = ad.obs["severity"].value_counts().to_dict()
        out[f"{screen}._severity_levels"] = {str(k): int(v) for k, v in sv.items()}
        # `severity` is numeric, with no legend in the object. The mapping is
        # read off the data rather than assumed: at severity 1.0 there are 863
        # distinct targets against the paper's 874 missense variants, and no
        # other level is within an order of magnitude of that (the next
        # closest, 0.5, has 696). The 11-variant shortfall is the same
        # filtering that costs this deposit 288 variants overall.
        neg_groups = set(out.get(f"{screen}._negctrl_groups", []))
        m = (~ad.obs[grp].astype(str).isin(neg_groups | {"PosCtrl"})) if grp \
            else slice(None)
        sub = ad.obs.loc[m]
        mis = sub.loc[sub["severity"] == 1.0, tv].nunique()
        out[f"{screen}.n_missense"] = int(mis)
        out[f"{screen}._missense_severity_level"] = 1.0

    # Replicate agreement, the paper's technical-reproducibility figure. It
    # correlates gRNA counts BETWEEN replicates within the same bin -- across
    # bins would measure the sort, not the reproducibility.
    if {"rep", "bin"} <= set(ad.var.columns):
        X = ad.X if not hasattr(ad.X, "toarray") else ad.X.toarray()
        reps = list(dict.fromkeys(ad.var["rep"].astype(str)))
        bins = list(dict.fromkeys(ad.var["bin"].astype(str)))
        rhos = []
        for b in bins:
            cols = [i for i, (r, bb_) in enumerate(
                zip(ad.var["rep"].astype(str), ad.var["bin"].astype(str)))
                if bb_ == b]
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    rhos.append(spearman(list(map(float, X[:, cols[i]])),
                                          list(map(float, X[:, cols[j]]))))
        rhos = [r for r in rhos if r == r]
        if rhos:
            out[f"{screen}.replicate_rho"] = round(_percentile(sorted(rhos), 0.5), 4)
            out[f"{screen}._n_replicate_pairs"] = len(rhos)
            out[f"{screen}._replicates"] = reps
            out[f"{screen}._bins"] = bins

    # Editing rates, from the reporter the paper's deposit carries and IGVF
    # does not publish.
    if "edit_rate" in ad.obs.columns:
        er = np.asarray(ad.obs["edit_rate"], dtype=float)
        er = er[~np.isnan(er)]
        if er.size:
            out[f"{screen}.mean_edit_fraction"] = round(float(er.mean()), 4)
        if tv and grp:
            mask = ad.obs[grp].astype(str).str.lower().str.startswith("variant")
            sub = ad.obs.loc[mask, [tv, "edit_rate"]].dropna()
            if len(sub):
                mx = sub.groupby(tv)["edit_rate"].max()
                out[f"{screen}.median_max_edit"] = round(float(mx.median()), 4)
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


def cmd_measure(args) -> int:
    setup_logging()
    names = list(SCREENS) if args.screen == "all" else [args.screen]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / "measured.json"
    measured = json.loads(p.read_text()) if p.exists() else {}
    for n in names:
        out = measure(n)
        if "error" in out:
            print(f"  {n}: {out['error']}")
            continue
        measured.update(out)
        for k, v in sorted(out.items()):
            print(f"  {k:34} {v}")
    p.write_text(json.dumps(measured, indent=2))
    print(f"\n  -> {p}")
    return 0


def cmd_run(args) -> int:
    setup_logging()
    models = list(MODELS) if args.model == "all" else [args.model]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mp = OUT_DIR / "measured.json"
    measured = json.loads(mp.read_text()) if mp.exists() else {}
    rc = 0
    for m in models:
        print(f"  {MODELS[m]['label']} …")
        out = run_model(args.screen, m, n_iter=args.iter)
        if out.get("error"):
            print(f"    {out['error']}")
            rc = 2
            continue
        if out["exit_code"] != 0:
            print(f"    FAILED ({out['exit_code']}): "
                  f"{out['stderr'].splitlines()[-1] if out['stderr'] else ''}")
            rc = 2
            continue
        print(f"    ok -> {out.get('result_csv')}")
        if out.get("result_csv"):
            sc = score_auprc(args.screen, out["result_csv"])
            key = f"{args.screen}.{MODELS[m]['claim']}"
            measured[key] = round(sc["auprc"], 4)
            measured[f"{key}._n_pos"] = sc["n_positive"]
            measured[f"{key}._n_neg"] = sc["n_negative"]
            print(f"    AUPRC {sc['auprc']:.3f}  "
                  f"({sc['n_positive']} positives vs {sc['n_negative']} negatives)")
    mp.write_text(json.dumps(measured, indent=2))
    return rc


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

    m = sub.add_parser("measure", help="Compute what the deposited object "
                                        "alone supports — no model fit.")
    m.add_argument("screen", choices=list(SCREENS) + ["all"], default="all",
                    nargs="?")

    r = sub.add_parser("run", help="Fit the paper's model variants on the "
                                    "deposit and score the AUPRC.")
    r.add_argument("screen", choices=list(SCREENS))
    r.add_argument("--model", choices=list(MODELS) + ["all"], default="all")
    r.add_argument("--iter", type=int, default=0,
                    help="Override BEAN --n-iter. 0 = BEAN's default.")

    sub.add_parser("report", help="Measured vs published, claim by claim.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"fetch": cmd_fetch, "describe": cmd_describe,
            "measure": cmd_measure, "run": cmd_run,
            "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
