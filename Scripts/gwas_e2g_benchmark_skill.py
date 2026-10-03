#!/usr/bin/env python3
"""GWAS variant benchmark for enhancer-gene predictors (port of EngreitzLab/GWAS_E2G_benchmarking).

Port of https://github.com/EngreitzLab/GWAS_E2G_benchmarking (MIT; Snakemake +
R + bedtools + csvtk), the Engreitz-lab pipeline used in the scE2G and
ENCODE-rE2G papers to benchmark enhancer-gene (E2G) predictions against
fine-mapped UK Biobank GWAS variants.  It asks two questions: (1) variant
overlap -- how many fine-mapped noncoding GWAS variants for a trait fall in
the predicted enhancers of a biosample, and how enriched are they relative to
common 1000G SNPs; (2) gene linking -- how precisely do the predicted
enhancer-gene links recover the silver-standard causal gene of each credible
set (Weeks et al. 2023), alone and intersected with PoPS.  Every rule and R
script of the workflow was read and re-derived in Python (pandas + numpy +
scipy); no code was copied.  Relationship: port.  Pinned upstream commit
c43bee694c6abfa1cf54477fc433ba8e370ded18 (2025-10-27).

Definitions, exactly as upstream computes them
  variants          per trait (variant key -> variant file with chr, start, end,
                    rsid, pip, CredibleSet, trait): pip > thresholdPIP (strict,
                    default 0.1), chromosome in chrSizes, overlapping the
                    distal-noncoding partition (ABC / AllPeaks / Other /
                    OtherIntron).  Trait groups (plus ALL = every trait in the
                    key) concatenate their members, de-duplicate on
                    (chr, start, end, rsid, CredibleSet) and take the group
                    name as trait.  nVariantsTotal = distinct (chr, start, end)
                    per trait.
  background        1000G SNPs (chr, start, end, rsid) on the same partition;
                    nCommonVariantsTotal = their count.
  predictions       chr, start, end, TargetGene, <score_col> ('#' stripped from
                    the header, rows with blanks dropped), TargetGene in the
                    gene universe (4th column of the TSS reference), score and
                    threshold negated for inverse predictors; thresholded at
                    score >= threshold.  Biosample groups (plus ALL = every
                    biosample of the method; a configured group is used only
                    when all its members have predictions) concatenate the
                    members' thresholded predictions (gene linking) and merge
                    them (variant overlap, enhancer set size).
  enhancer set size bpEnhancers = bp of the merged thresholded elements.
  overlap           nVariantsOverlappingEnhancers = distinct variant locations of
                    the trait inside thresholded elements;
                    nCommonVariantsOverlappingEnhancers = distinct background
                    SNPs inside them; recall = nVO / nVariantsTotal;
                    enrichment = recall / (nCVO / nCommonVariantsTotal).
                    Trait x biosample rows exist only when nVO > 0.
  helpful_math.R    z = qnorm(thresholdPval / 2, lower = FALSE);
                    SE_log_enr = sqrt((n1-x1)/(x1 n1) + (n2-x2)/(x2 n2)) (see the
                    quirk below); CI_enr = exp(log(enr) -/+ z SE);
                    p_enr = phyper(x1, n1, n2, x1 + x2, lower = FALSE);
                    p_adjust_enr = Bonferroni over the table (dplyr semantics:
                    within the remaining groups after a grouped summarize);
                    recall_adjust = (x + 2)/(n + 4), SE_recall =
                    sqrt(ra (1 - ra)/(n + 4)), CI_recall = recall -/+ z SE_recall
                    (Agresti-Coull width centred on the raw recall, as upstream);
                    precision CI the same with TP and predicted positives.
  threshold span    per method, over the full (unthresholded) variant
                    intersections of the biosample x trait pairs in the
                    comparisons table: quantiles (R type 7) of the per-variant
                    max score at probs seq(0, 1, length = n) with n grown from
                    nThresholdSteps (25) until that many distinct thresholds or
                    the number of distinct scores, plus the method threshold,
                    plus an even grid min..max of nThresholdSteps; boolean
                    methods use n = 2.
  curves            enrichmentRecallAcrossThresholds.forBiosample<b>: the overlap
                    statistics at every span threshold for each comparison
                    biosample (groups = union of members); aggregated
                    "all_matched" = summed counts per (threshold, method) over a
                    comparison's biosample.trait keys; plotted points keep
                    recall > 0.005, p_adjust_enr < thresholdPval,
                    nVariantsTotal > 25 and nCVO / nCVT > 0.001.
  gene linking      gene prioritisation table (CredibleSet, Disease -> trait,
                    TargetGene, truth, POPS.Score, PromoterDistanceToBestSNP):
                    nCredibleSetsTotal = distinct credible sets per trait; truth
                    genes = rows with truth TRUE.  Per biosample x trait with any
                    overlapping prediction: restrict to truth credible sets, max
                    score per (CredibleSet, TargetGene), keep genes with
                    rank(-max, ties.method = "min") <= numPredGenes (2);
                    nCredibleSetsOverlappingEnhancers = distinct (CS, truth gene),
                    ...AnyGene = distinct (CS, TargetGene), ...CorrectGene = those
                    with TargetGene == truth; recall = correct / total, precision
                    = correct / any (NaN -> 0).  intersectPoPS = the same after
                    keeping only genes among the top numPoPSGenes (2) by
                    POPS.Score (row_number after sorting, file order breaks ties).
  baselines         PoPS: top numPoPSGenes by POPS.Score per credible set;
                    distanceToTSS: top numPoPSGenes by PromoterDistanceToBestSNP
                    (ascending); nCredibleSetGenePredictions(Correct) per trait;
                    the table is restricted to the gene universe.
  comparisons       variant overlap: rows whose biosample.trait key is in the
                    comparison, counts summed per method; gene linking: rows of
                    the comparison's biosamples (a group expands to its members,
                    ALL to every row) and traits (a trait group expands to its
                    members), de-duplicated, counts summed per method x
                    intersectPoPS; baselines per trait.
  metric ranges     per method min/max of enrichment_overlap (p_adjust_enr <
                    thresholdPval), recall_overlap (nVariantsTotal > 25,
                    nCVO/nCVT > 0.001) and precision/recall_linking (+/- PoPS).
  heatmaps          biosample x trait matrices (NA -> 0); rows and columns
                    ordered by hclust(dist(1 - cor(M)), "ward.D2"), kept in input
                    order when all distances are 0; colour limits 0..max
                    (plotFixedScale: 75th percentile of the per-method maxima
                    for enrichment, max for the others).

Upstream quirks, corrected by default and reproduced with --upstream-compat
  * helpful_math.R: SE_log_enr = sqrt(a) + b (the sqrt closes after the first
    term); the port uses sqrt(a + b), the textbook formula the code cites.
  * generate_quantile_threshold_span.R: `group_by(-c(predScore))` groups by the
    negated score, so the quantiles are taken over DISTINCT score values
    rather than per-variant maxima; the port uses per-variant maxima.
  * link_variants_genes.smk: individual biosamples are linked with their
    UNTHRESHOLDED predictions (enhancerPredictions.variantIntersection) while
    biosample groups use thresholded ones; the port thresholds both.
  * evaluate_gene_linking.R: the gene-universe filter of the prioritisation
    table is overwritten by a second fread (the baseline script keeps it); the
    port filters in both.
  * plot_thresholded_performance_comparison.R: biosample ALL in the gene-linking
    comparison also sums the biosample-GROUP rows (double counting); the port
    uses the individual biosamples.
  Other deliberate differences (not switchable): the group / biosample and
  trait-group / trait name-clash checks in utils.smk iterate over dict items and
  never fire -- they are enforced here; figures are PNG and the six per-method
  PDFs are drawn as two panels (overlap, linking); dendrogram leaf order comes
  from scipy and can differ from R's hclust order for equal-height merges.

Subcommands
  setup            fetch the small upstream resources (variant key, gene
                   prioritisation table ~9.5 MB, partition ~50 MB, TSS
                   reference, chromosome sizes) from the pinned commit into
                   Data/GWASE2G/resources; --traits also fetches per-trait
                   UK Biobank SuSiE variant lists (all 94 = 247 MB);
                   --synapse pulls the 1000G background SNPs (syn52264319).
  validate-config  resolve config.yml / flags, check methods, biosample and
                   trait groups and the comparisons table (ALL allowed), list
                   the biosamples and groups available per method.
  variants         process_variants.smk: PIP / chromosome / partition filters,
                   trait groups + ALL, background SNP filter; cached for run.
  baseline         evaluate_baseline_predictors.R: PoPS and distance-to-TSS
                   gene-linking precision / recall per trait.
  run              the whole workflow for a methods table + predictions table
                   + comparisons table: per method thresholded enrichment /
                   recall (trait x biosample), background overlap + enhancer
                   set size, threshold span, enrichment-recall curves, gene
                   linking; baselines, colour palette, metric ranges, heatmap
                   matrices, comparison tables (curves + thresholded
                   performance), figures, report.md, summary.json.
  plot             visualization.smk from an existing run directory (metric
                   ranges, heatmaps, curves, performance comparisons).
  selftest         synthetic genome with planted enhancers, credible sets,
                   silver-standard genes and PoPS scores; runs every
                   subcommand and asserts the definitions above.

Output: Docs/GWASE2G/<timestamp>_<label>/ laid out like the upstream results
directory (<method>/variant_overlap, <method>/gene_linking, baseline/,
plots/, reference_configs/).

Usage:
    igvfagent gwas-e2g setup [--traits RBC MCV Lym] [--synapse]
    igvfagent gwas-e2g validate-config --config config/config.yml
    igvfagent gwas-e2g variants --config config/config.yml --label ukbb
    igvfagent gwas-e2g baseline --gene-prioritization-table UKBiobank.ABCGene.anyabc.tsv
    igvfagent gwas-e2g run --config config/config.yml --label sce2g_blood
    igvfagent gwas-e2g run --methods-table methods.tsv --predictions-table preds.tsv --comparisons-table comparisons.tsv --variant-key key.tsv --bg-variants bg.bed --trait-group RBC_traits=MCH,MCV,Hb
    igvfagent gwas-e2g plot --run-dir Docs/GWASE2G/<run>
    igvfagent gwas-e2g selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eqtl_enrichment_skill import bonferroni, enrichment_stats, merged_bp, overlap_join  # noqa: E402

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "GWASE2G"
DATA_ROOT = ROOT / "Data" / "GWASE2G"

UPSTREAM_REPO = "EngreitzLab/GWAS_E2G_benchmarking"
UPSTREAM_COMMIT = "c43bee694c6abfa1cf54477fc433ba8e370ded18"
UPSTREAM_DATE = "2025-10-27T19:23:46Z"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
RESOURCE_FILES = {
    "resources/UKBB_variant_key.tsv": "variant key: 94 UK Biobank traits -> SuSiE fine-mapped variant lists",
    "resources/UKBiobank.ABCGene.anyabc.tsv": "gene prioritisation table: credible set x candidate gene, truth, POPS.Score, distance (~9.5 MB)",
    "resources/genome_annotation/PartitionCombined.bed": "hg38 partition ABC / Other / OtherIntron / TSS / CDS / UTR / splice (~50 MB)",
    "resources/genome_annotation/CollapsedGeneBound.hg38.intGENCODEv43.TSS500bp.bed": "TSS reference and gene universe",
    "resources/genome_annotation/GRCh38_main.chrom.sizes.tsv": "chromosome sizes (main chromosomes)",
}
DEFAULT_RES = {
    "variantKey": "resources/UKBB_variant_key.tsv",
    "genePrioritizationTable": "resources/UKBiobank.ABCGene.anyabc.tsv",
    "partition": "resources/genome_annotation/PartitionCombined.bed",
    "TSS": "resources/genome_annotation/CollapsedGeneBound.hg38.intGENCODEv43.TSS500bp.bed",
    "chrSizes": "resources/genome_annotation/GRCh38_main.chrom.sizes.tsv",
}
VARIANT_LIST_PATH = "resources/191010_UKBB_SuSiE_hg38_liftover/{trait}/variant.list.txt"
SYNAPSE_BG = "syn52264319"
DISTAL_CATEGORIES = {"ABC", "AllPeaks", "Other", "OtherIntron"}
DEFAULTS: "dict[str, Any]" = {"nThresholdSteps": 25, "thresholdPIP": 0.1, "thresholdPval": 0.05, "numPoPSGenes": 2,
                              "numPredGenes": 2, "plotFixedScale": False, "biosampleGroups": {}, "traitGroups": {},
                              "methods": []}
PATH_KEYS = ["methodsTable", "predictionsTable", "comparisonsTable", "variantKey", "genePrioritizationTable",
             "bgVariants", "partition", "TSS", "chrSizes"]
VAR_COLS = ["chr", "start", "end", "rsid", "CredibleSet", "trait"]
BASELINE_COLORS = [("PoPS", "Polygenic priority score (PoPS)", "#1c2a43"), ("distanceToTSS", "Distance to TSS", "#c5cad7")]
FILL_COLORS = ["#E495A5", "#BDAB66", "#65BC8C", "#55B8D0", "#C29DDE", "#D7A186", "#8FB67C", "#6FB3C4"]  # HCL 'Set 2'-like
ENR_OVL_COLORS = ["#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d"]
RECALL_OVL_COLORS = ["#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"]
PRECISION_COLORS = ["#edf8fb", "#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c"]
RECALL_LINK_COLORS = ["#ffffcc", "#c2e699", "#78c679", "#31a354", "#006837"]
INK, INK2, AXIS, SURFACE = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"

log = logging.getLogger("gwas_e2g")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"gwas_e2g_benchmark_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=8)
    ax.set_ylabel(ylabel, color=INK2, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=7)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path, figures: "list[str]") -> Path:
    fig.patch.set_facecolor(SURFACE)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    figures.append(str(path))
    return path


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        f = float(x)
        if abs(f) >= 1000:
            return f"{f:,.0f}"
        return f"{f:.{nd}g}" if 0 < abs(f) < 1e-3 else f"{f:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "t", "1", "yes", "y"}


def _num(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def write_tsv(df, path: Path, quiet: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)
    if not quiet:
        print(f"TSV: {path}")
    return path


def _header_line(path: Path) -> "tuple[list[str], int]":
    """Header columns and lines to skip; '#' removed from the header (awk 'NR==1{sub(/^#*/, "")}')."""
    opener = gzip.open if str(path).endswith(".gz") else open
    skip = 0
    with opener(path, "rt") as fh:
        for line in fh:
            skip += 1
            if not line.strip():
                continue
            if line.startswith("#") and "\t" not in line:
                continue
            return line.lstrip("#").rstrip("\n").rstrip("\r").split("\t"), skip
    return [], skip


def read_tsv(path: Path, usecols=None, dtype=None, **kw):
    pd = _pd()
    cols, skip = _header_line(path)
    return pd.read_csv(path, sep="\t", names=cols, skiprows=skip, usecols=usecols, dtype=dtype, low_memory=False, **kw)


def _chrom_key(c: str):
    s = str(c).replace("chr", "")
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


# ---------------------------------------------------------------------------
# Interval helpers (bedtools-free)
# ---------------------------------------------------------------------------

def merge_intervals(df) -> "dict[str, tuple[Any, Any]]":
    """bedtools merge per chromosome -> {chr: (starts, ends)} of disjoint, sorted intervals."""
    np = _np()
    out = {}
    for chrom, sub in df.groupby("chr", sort=False):
        s = sub["start"].to_numpy(dtype=np.int64); e = sub["end"].to_numpy(dtype=np.int64)
        o = np.argsort(s, kind="mergesort"); s, e = s[o], e[o]
        ce = np.maximum.accumulate(e)
        new = np.ones(s.size, dtype=bool)
        new[1:] = s[1:] > ce[:-1]
        grp = np.cumsum(new) - 1
        ms = s[new]
        me = np.zeros(ms.size, dtype=np.int64)
        np.maximum.at(me, grp, e)
        out[str(chrom)] = (ms, me)
    return out


def overlaps_union(df, merged: "dict[str, tuple[Any, Any]]"):
    """Boolean mask: row interval [start, end) overlaps the union of `merged` by >= 1 bp."""
    np = _np()
    mask = np.zeros(len(df), dtype=bool)
    if len(df) == 0:
        return mask
    chrs = df["chr"].astype(str).to_numpy()
    st = df["start"].to_numpy(dtype=np.int64); en = np.maximum(df["end"].to_numpy(dtype=np.int64), st + 1)
    for chrom in np.unique(chrs):
        if chrom not in merged:
            continue
        ms, me = merged[chrom]
        idx = np.flatnonzero(chrs == chrom)
        j = np.searchsorted(me, st[idx], side="right")
        ok = j < ms.size
        hit = np.zeros(idx.size, dtype=bool)
        hit[ok] = ms[j[ok]] < en[idx][ok]
        mask[idx] = hit
    return mask


# ---------------------------------------------------------------------------
# Configuration (config.yml + TSV tables), validation
# ---------------------------------------------------------------------------

def _yaml_scalar(v: str):
    v = v.strip()
    if not v:
        return None
    if (v[0] == v[-1]) and v[0] in "\"'" and len(v) >= 2:
        return v[1:-1]
    if v.startswith("[") and v.endswith("]"):
        inner = v[1:-1].strip()
        return [_yaml_scalar(x) for x in inner.split(",")] if inner else []
    if v.lower() in {"true", "false"}:
        return v.lower() == "true"
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


def _strip_comment(line: str) -> str:
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch
        elif ch == "#":
            break
        out.append(ch)
    return "".join(out).rstrip()


def parse_config_text(text: str) -> dict:
    """The subset of YAML the upstream config.yml uses (scalars, inline lists, one level of nested
    mappings, '- item' lists); PyYAML is used when it is installed."""
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text) or {}
    except ImportError:
        pass
    cfg: dict = {}
    parent: Optional[str] = None
    for raw in text.splitlines():
        line = _strip_comment(raw)
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        body = line.strip()
        if indent == 0:
            if ":" not in body:
                continue
            k, v = body.split(":", 1)
            k = k.strip()
            if v.strip():
                cfg[k] = _yaml_scalar(v); parent = None
            else:
                cfg[k] = None; parent = k
        elif parent is not None:
            if body.startswith("- "):
                if not isinstance(cfg[parent], list):
                    cfg[parent] = []
                cfg[parent].append(_yaml_scalar(body[2:]))
            elif ":" in body:
                if not isinstance(cfg[parent], dict):
                    cfg[parent] = {}
                k, v = body.split(":", 1)
                cfg[parent][k.strip()] = _yaml_scalar(v)
    return cfg


def _parse_groups(items: "Optional[list[str]]") -> "dict[str, list[str]]":
    out: "dict[str, list[str]]" = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"group spec must be NAME=member1,member2: {it}")
        k, v = it.split("=", 1)
        out[k.strip()] = [x.strip() for x in v.split(",") if x.strip()]
    return out


def _resolve(p, bases: "list[Path]") -> Optional[Path]:
    if p is None or str(p).strip() == "":
        return None
    q = Path(str(p)).expanduser()
    if q.is_absolute():
        return q
    for b in bases:
        if (b / q).exists():
            return (b / q).resolve()
    return (bases[0] / q) if bases else q


def load_config(args: argparse.Namespace) -> dict:
    """config.yml (optional) + command-line overrides -> resolved config dict."""
    cfg = json.loads(json.dumps(DEFAULTS))
    bases: "list[Path]" = []
    if getattr(args, "config", None):
        cp = Path(args.config).expanduser().resolve()
        if not cp.is_file():
            raise SystemExit(f"config not found: {cp}")
        text = cp.read_text()
        loaded = json.loads(text) if cp.suffix == ".json" else parse_config_text(text)
        for k, v in (loaded or {}).items():
            if v is not None:
                cfg[k] = v
        if cfg.get("baseDir") and Path(str(cfg["baseDir"])).is_dir():
            bases.append(Path(str(cfg["baseDir"])))
        bases += [cp.parent, cp.parent.parent]
    bases += [Path.cwd(), DATA_ROOT]
    over = {"methodsTable": "methods_table", "predictionsTable": "predictions_table", "comparisonsTable": "comparisons_table",
            "variantKey": "variant_key", "genePrioritizationTable": "gene_prioritization_table", "bgVariants": "bg_variants",
            "partition": "partition", "TSS": "tss", "chrSizes": "chr_sizes", "nThresholdSteps": "n_threshold_steps",
            "thresholdPIP": "threshold_pip", "thresholdPval": "threshold_pval", "numPoPSGenes": "num_pops_genes",
            "numPredGenes": "num_pred_genes"}
    for k, a in over.items():
        v = getattr(args, a, None)
        if v is not None:
            cfg[k] = v
    if getattr(args, "methods", None):
        cfg["methods"] = list(args.methods)
    if getattr(args, "plot_fixed_scale", False):
        cfg["plotFixedScale"] = True
    cfg["biosampleGroups"] = dict(cfg.get("biosampleGroups") or {})
    cfg["traitGroups"] = dict(cfg.get("traitGroups") or {})
    cfg["biosampleGroups"].update(_parse_groups(getattr(args, "biosample_group", None)))
    cfg["traitGroups"].update(_parse_groups(getattr(args, "trait_group", None)))
    for k in ("biosampleGroups", "traitGroups"):
        cfg[k] = {str(g): [str(x) for x in (v or [])] for g, v in cfg[k].items()}
    for k in PATH_KEYS:
        if cfg.get(k) in (None, "") and k in DEFAULT_RES:
            cfg[k] = str(DATA_ROOT / DEFAULT_RES[k])
        r = _resolve(cfg.get(k), bases)
        cfg[k] = str(r) if r else None
    for k in ("nThresholdSteps", "numPoPSGenes", "numPredGenes"):
        cfg[k] = int(cfg[k])
    for k in ("thresholdPIP", "thresholdPval"):
        cfg[k] = float(cfg[k])
    cfg["plotFixedScale"] = _truthy(cfg.get("plotFixedScale"))
    cfg["_bases"] = [str(b) for b in bases]
    return cfg


def read_variant_key(cfg: dict):
    pd = _pd()
    kp = Path(cfg["variantKey"])
    key = pd.read_csv(kp, sep="\t", dtype=str).dropna(subset=["trait"])
    bases = [kp.parent, kp.parent.parent] + [Path(b) for b in cfg["_bases"]]
    key["path"] = [str(_resolve(v, bases)) for v in key["variant_file"]]
    return key


class Method:
    def __init__(self, row: dict, biosamples: "list[str]", files: "list[str]", groups: "list[str]"):
        self.name = str(row["method"])
        self.boolean = _truthy(row.get("boolean", False))
        self.inverse = _truthy(row.get("inverse_predictor", False))
        self.long = str(row.get("pred_name_long") or self.name)
        thr = float(row["threshold"])
        self.threshold = -thr if self.inverse else thr
        self.score_col = str(row["score_col"])
        c = str(row.get("color") or "").strip()
        self.color = c if c and c.lower() not in {"nan", "none"} else ""
        self.biosamples = biosamples
        self.files = dict(zip(biosamples, files))
        self.groups = groups          # "ALL" + configured groups fully represented

    def as_dict(self) -> dict:
        return {"method": self.name, "boolean": self.boolean, "inverse_predictor": self.inverse, "pred_name_long": self.long,
                "threshold": self.threshold, "score_col": self.score_col, "color": self.color, "biosamples": self.biosamples,
                "predFiles": [self.files[b] for b in self.biosamples], "biosampleGroups": self.groups}


def read_tables(cfg: dict):
    pd = _pd()
    mt = pd.read_csv(cfg["methodsTable"], sep="\t", dtype=str).fillna("")
    pt = pd.read_csv(cfg["predictionsTable"], sep="\t", dtype=str)
    pt = pt.dropna(subset=["biosample"])
    pt["biosample"] = pt["biosample"].str.strip()
    comp = pd.read_csv(cfg["comparisonsTable"], sep="\t", dtype=str) if cfg.get("comparisonsTable") else \
        pd.DataFrame({"name": ["all"], "biosample": ["ALL"], "trait": ["ALL"]})
    comp = comp.dropna(how="all")
    for c in ("name", "biosample", "trait"):
        if c in comp.columns:
            comp[c] = comp[c].astype(str).str.strip()
    return mt, pt, comp


def build_methods(cfg: dict, mt, pt) -> "list[Method]":
    names = cfg["methods"] or [m for m in mt["method"] if m in pt.columns]
    rows = {r["method"]: r for r in mt.to_dict("records")}
    bases = [Path(cfg["predictionsTable"]).parent] + [Path(b) for b in cfg["_bases"]]
    out = []
    for m in names:
        if m not in rows or m not in pt.columns:
            continue
        sub = pt[["biosample", m]].dropna()
        sub = sub[sub[m].astype(str).str.strip() != ""]
        bs = [str(b) for b in sub["biosample"]]
        files = [str(_resolve(str(f).strip(), bases)) for f in sub[m]]
        groups = ["ALL"] + [g for g, mem in cfg["biosampleGroups"].items() if all(b in bs for b in mem)]
        try:
            out.append(Method(rows[m], bs, files, groups))
        except (TypeError, ValueError) as e:
            raise SystemExit(f"method {m}: bad methods-table row ({e})") from e
    return out


def validate(cfg: dict, check_files: bool = True) -> "tuple[list[str], list[str], Any]":
    """utils.smk validate_biosample_groups / validate_trait_groups + table / file checks."""
    errors: "list[str]" = []
    warnings: "list[str]" = []
    for k in ("methodsTable", "predictionsTable", "variantKey"):
        if not cfg.get(k) or not Path(cfg[k]).is_file():
            errors.append(f"{k} not found: {cfg.get(k)}")
    if errors:
        return errors, warnings, None
    mt, pt, comp = read_tables(cfg)
    for c in ("method", "boolean", "inverse_predictor", "pred_name_long", "threshold", "score_col"):
        if c not in mt.columns:
            errors.append(f"methods table lacks column {c}")
    if "biosample" not in pt.columns:
        errors.append("predictions table lacks column biosample")
    for c in ("name", "biosample", "trait"):
        if c not in comp.columns:
            errors.append(f"comparisons table lacks column {c}")
    if errors:
        return errors, warnings, None
    for m in cfg["methods"]:
        if m not in set(mt["method"]):
            errors.append(f"method {m} is not in the methods table")
        if m not in pt.columns:
            errors.append(f"method {m} has no column in the predictions table")
    bios = set(pt["biosample"])
    key = read_variant_key(cfg)
    traits = set(key["trait"])
    for g, mem in cfg["biosampleGroups"].items():
        if g in bios or g == "ALL":
            errors.append(f"Biosample group and biosample {g} cannot share the same name.")
        for b in mem:
            if b not in bios:
                errors.append(f"Biosample {b} in group {g} is not defined in predictions config.")
    for g, mem in cfg["traitGroups"].items():
        if g in traits:
            errors.append(f"Trait group and trait {g} cannot share the same name.")
        if g == "ALL":
            continue
        for t in mem:
            if t not in traits:
                errors.append(f"Trait {t} in group {g} is not defined in variant key.")
    for b in comp["biosample"]:
        if b not in cfg["biosampleGroups"] and b not in bios and b != "ALL":
            errors.append(f"Biosample {b} in comparisons config is not defined.")
    for t in comp["trait"]:
        if t not in cfg["traitGroups"] and t not in traits and t != "ALL":
            errors.append(f"Trait {t} in comparisons config is not defined.")
    methods = build_methods(cfg, mt, pt)
    if check_files:
        for m in methods:
            for b, f in m.files.items():
                if not Path(f).is_file():
                    errors.append(f"prediction file missing ({m.name}, {b}): {f}")
        for t, p in zip(key["trait"], key["path"]):
            if not Path(p).is_file():
                warnings.append(f"variant file missing for trait {t}: {p}")
        for k in ("partition", "TSS", "chrSizes", "bgVariants", "genePrioritizationTable"):
            if not cfg.get(k) or not Path(cfg[k]).is_file():
                (errors if k != "genePrioritizationTable" else warnings).append(f"{k} not found: {cfg.get(k)}")
    return errors, warnings, (mt, pt, comp, key, methods)


def cmd_validate(args: argparse.Namespace) -> int:
    cfg = load_config(args)
    errors, warnings, res = validate(cfg, check_files=not args.no_file_checks)
    for k in PATH_KEYS:
        print(f"  {k:<24} {cfg.get(k)}")
    if res is not None:
        mt, pt, comp, key, methods = res
        print(f"traits in variant key: {len(key)}; trait groups: {', '.join(cfg['traitGroups']) or '-'} (+ ALL)")
        print(f"comparisons: {', '.join(comp['name'].unique())}")
        for m in methods:
            print(f"  method {m.name}: threshold {m.threshold:g}{' (inverted)' if m.inverse else ''}; "
                  f"{len(m.biosamples)} biosamples; groups {', '.join(m.groups)}")
    for w in warnings:
        print("  WARN ", w)
    for e in errors:
        print("  ERROR", e)
    print("config valid" if not errors else f"config invalid: {len(errors)} error(s)")
    return 0 if not errors else 1


# ---------------------------------------------------------------------------
# setup
# ---------------------------------------------------------------------------

def _fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=180) as r, open(tmp, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.replace(dest)


def cmd_setup(args: argparse.Namespace) -> int:
    dest_root = Path(args.dest).resolve() if args.dest else DATA_ROOT
    todo = dict(RESOURCE_FILES)
    if args.traits:
        key_path = dest_root / "resources" / "UKBB_variant_key.tsv"
        if not key_path.is_file():
            _fetch(RAW + "resources/UKBB_variant_key.tsv", key_path)
        all_traits = [ln.split("\t")[0] for ln in key_path.read_text().splitlines()[1:] if ln.strip()]
        want = all_traits if [t.lower() for t in args.traits] == ["all"] else args.traits
        unknown = [t for t in want if t not in all_traits]
        if unknown:
            print(f"unknown traits (not in the UKBB variant key): {unknown}", file=sys.stderr)
        for t in want:
            if t in all_traits:
                todo[VARIANT_LIST_PATH.format(trait=t)] = f"SuSiE fine-mapped variants, {t}"
        if len(want) > 20:
            print(f"note: {len(want)} per-trait variant lists (~247 MB for all 94 traits)")
    got, skipped, failed = [], [], []
    for rel in todo:
        dest = dest_root / rel
        if dest.is_file() and dest.stat().st_size > 0 and not args.force:
            skipped.append(rel)
            continue
        try:
            _fetch(RAW + rel, dest)
            got.append(rel)
            print(f"Wrote: {dest}")
        except Exception as e:  # noqa: BLE001
            failed.append(f"{rel}: {e}")
    print(f"Resources in {dest_root}: fetched {len(got)}, already present {len(skipped)}, failed {len(failed)} "
          f"(pinned {UPSTREAM_REPO}@{UPSTREAM_COMMIT[:10]})")
    for f in failed:
        print("  FAILED", f, file=sys.stderr)
    if args.synapse:
        tok = os.environ.get("SYNAPSE_AUTH_TOKEN") or os.environ.get("SYNAPSE_PAT")
        if not tok:
            print("SYNAPSE_AUTH_TOKEN is not set; skipping the background-SNP download (syn52264319)", file=sys.stderr)
        else:
            cmd = [sys.executable, "-m", "igvfagent.cli", "synapse", "download", "--syn", SYNAPSE_BG, "--out-dir", str(dest_root)]
            print("Synapse", SYNAPSE_BG, "(all.bg.SNPs.hg38.baseline.v1.1.bed.sorted, ~10M 1000G SNPs)")
            if subprocess.call(cmd) != 0:
                failed.append(f"{SYNAPSE_BG} download failed")
    else:
        print(f"background SNPs: `igvfagent gwas-e2g setup --synapse` (Synapse {SYNAPSE_BG}) or pass --bg-variants")
    return 1 if failed else 0


# ---------------------------------------------------------------------------
# Variants (process_variants.smk)
# ---------------------------------------------------------------------------

def read_chr_sizes(path: str) -> "list[str]":
    with open(path) as fh:
        return [ln.split("\t")[0].strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def load_partition_distal(path: str, chroms: "list[str]"):
    pd = _pd()
    part = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "category"],
                       dtype={"chr": str, "category": str}, comment="#")
    part = part[part["chr"].isin(set(chroms)) & part["category"].isin(DISTAL_CATEGORIES)]
    return part[["chr", "start", "end"]].reset_index(drop=True)


def load_gene_universe(path: str) -> "set[str]":
    pd = _pd()
    tss = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[3], names=["name"], dtype=str)
    return set(tss["name"].dropna())


def _sort_bed(df, chroms: "list[str]"):
    order = {c: i for i, c in enumerate(chroms)}
    df = df.assign(_o=df["chr"].map(order).fillna(len(order)))
    return df.sort_values(["_o", "start", "end"], kind="mergesort").drop(columns="_o").reset_index(drop=True)


def process_variants(cfg: dict, key=None) -> dict:
    """filter_GWAS_variants + merge_GWAS_variants_trait_group + combine_GWAS_variants + background filter."""
    pd = _pd()
    chroms = read_chr_sizes(cfg["chrSizes"])
    part = merge_intervals(load_partition_distal(cfg["partition"], chroms))
    key = read_variant_key(cfg) if key is None else key
    per_trait, stats = {}, []
    for t, p in zip(key["trait"], key["path"]):
        if not Path(p).is_file():
            log.warning("variant file missing for %s: %s", t, p)
            stats.append({"trait": t, "n_rows": 0, "n_pip": 0, "n_distal_noncoding": 0, "file": p, "missing": True})
            continue
        cols = set(_header_line(Path(p))[0])
        need = ["chr", "start", "end", "rsid", "pip", "CredibleSet"] + (["trait"] if "trait" in cols else [])
        miss = [c for c in need if c not in cols]
        if miss:
            raise SystemExit(f"{p}: missing columns {miss}")
        df = read_tsv(Path(p), usecols=need, dtype={"chr": str, "rsid": str, "CredibleSet": str, "trait": str})
        if "trait" not in df.columns:
            df["trait"] = t
        n0 = len(df)
        df["pip"] = pd.to_numeric(df["pip"], errors="coerce")
        df = df[df["pip"] > cfg["thresholdPIP"]]
        n1 = len(df)
        df = df[df["chr"].isin(set(chroms))].copy()
        df["start"] = df["start"].astype("int64"); df["end"] = df["end"].astype("int64")
        df = df[overlaps_union(df, part)]
        df = df[VAR_COLS]
        per_trait[t] = df
        stats.append({"trait": t, "n_rows": n0, "n_pip": n1, "n_distal_noncoding": len(df), "file": p, "missing": False})
    groups = dict(cfg["traitGroups"])
    if "ALL" not in groups or not groups["ALL"]:
        groups["ALL"] = list(key["trait"])
    pieces = list(per_trait.values())
    for g, mem in groups.items():
        sub = [per_trait[t][VAR_COLS[:5]] for t in mem if t in per_trait]
        if not sub:
            continue
        gdf = pd.concat(sub, ignore_index=True).drop_duplicates()
        gdf["trait"] = g
        pieces.append(gdf)
    V = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=VAR_COLS)
    V = _sort_bed(V, chroms)
    bg = None
    if cfg.get("bgVariants") and Path(cfg["bgVariants"]).is_file():
        bg = pd.read_csv(cfg["bgVariants"], sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "rsid"],
                         dtype={"chr": str, "rsid": str}, comment="#")
        n_bg0 = len(bg)
        bg = bg[bg["chr"].isin(set(chroms))]
        bg = bg[overlaps_union(bg, part)].drop_duplicates(["chr", "start", "end", "rsid"])
        bg = _sort_bed(bg, chroms)
    else:
        n_bg0 = 0
    return {"variants": V, "bg": bg, "stats": pd.DataFrame(stats), "trait_groups": groups, "n_bg_input": n_bg0}


def write_variants_dir(res: dict, d: Path, cfg: dict) -> None:
    vd = d / "variants"
    write_tsv(res["variants"], vd / "filteredGWASVariants.merged.sorted.tsv.gz")
    write_tsv(res["stats"], vd / "variantFilterStats.tsv")
    if res["bg"] is not None:
        res["bg"].to_csv(vd / "bgVariants.distalNoncoding.bed.gz", sep="\t", header=False, index=False, compression="gzip")
        print(f"Wrote: {vd / 'bgVariants.distalNoncoding.bed.gz'}")
        (vd / "bgVariants.distalNoncoding.count.txt").write_text(f"{len(res['bg'])}\n")
    (vd / "variants_meta.json").write_text(json.dumps({"thresholdPIP": cfg["thresholdPIP"], "traitGroups": res["trait_groups"],
                                                        "partition": cfg["partition"], "chrSizes": cfg["chrSizes"],
                                                        "n_bg_input": res["n_bg_input"]}, indent=2))


def load_variants_dir(d: Path) -> dict:
    pd = _pd()
    vd = d / "variants" if (d / "variants").is_dir() else d
    V = pd.read_csv(vd / "filteredGWASVariants.merged.sorted.tsv.gz", sep="\t", dtype={"chr": str, "rsid": str, "CredibleSet": str, "trait": str})
    bgp = vd / "bgVariants.distalNoncoding.bed.gz"
    bg = pd.read_csv(bgp, sep="\t", header=None, names=["chr", "start", "end", "rsid"], dtype={"chr": str, "rsid": str}) if bgp.is_file() else None
    meta = json.loads((vd / "variants_meta.json").read_text()) if (vd / "variants_meta.json").is_file() else {}
    stats = pd.read_csv(vd / "variantFilterStats.tsv", sep="\t") if (vd / "variantFilterStats.tsv").is_file() else pd.DataFrame()
    return {"variants": V, "bg": bg, "stats": stats, "trait_groups": meta.get("traitGroups", {}), "n_bg_input": meta.get("n_bg_input", 0)}


def n_variants_per_trait(V):
    return V.drop_duplicates(["chr", "start", "end", "trait"]).groupby("trait").size().rename("nVariantsTotal").reset_index()


def cmd_variants(args: argparse.Namespace) -> int:
    cfg = load_config(args)
    for k in ("variantKey", "partition", "chrSizes"):
        if not cfg.get(k) or not Path(cfg[k]).is_file():
            raise SystemExit(f"{k} not found: {cfg.get(k)} (run `igvfagent gwas-e2g setup` or pass it)")
    res = process_variants(cfg)
    out = run_dir(args.label)
    write_variants_dir(res, out, cfg)
    nv = n_variants_per_trait(res["variants"])
    n_bg = 0 if res["bg"] is None else len(res["bg"])
    summ = {"label": args.label, "thresholdPIP": cfg["thresholdPIP"], "n_traits": int((~res["stats"]["missing"]).sum()) if len(res["stats"]) else 0,
            "trait_groups": res["trait_groups"], "n_variant_rows": int(len(res["variants"])),
            "nVariantsTotal": {r.trait: int(r.nVariantsTotal) for r in nv.itertuples()},
            "background_snps_input": int(res["n_bg_input"]), "background_snps_distal_noncoding": int(n_bg), "variants_dir": str(out)}
    (out / "summary.json").write_text(json.dumps(summ, indent=2))
    rep = [f"# GWAS variant processing — {args.label}", "",
           f"PIP > {cfg['thresholdPIP']}, distal-noncoding partition ({', '.join(sorted(DISTAL_CATEGORIES))}), "
           f"chromosomes from {Path(cfg['chrSizes']).name}.", "",
           md_table(["trait", "rows", "PIP pass", "distal noncoding"],
                    [(r["trait"], r["n_rows"], r["n_pip"], r["n_distal_noncoding"]) for r in res["stats"].to_dict("records")], 120), "",
           "Trait groups (merged, de-duplicated): " + "; ".join(f"{g} ({len(m)})" for g, m in res["trait_groups"].items()), "",
           f"Background SNPs: {res['n_bg_input']} -> {n_bg} distal noncoding.", "",
           f"Reuse with `igvfagent gwas-e2g run --variants-dir {out}`."]
    (out / "report.md").write_text("\n".join(rep) + "\n")
    print(f"JSON: {out / 'summary.json'}")
    print(f"Report: {out / 'report.md'}")
    return 0


# ---------------------------------------------------------------------------
# Statistics (helpful_math.R)
# ---------------------------------------------------------------------------

def _z(alpha: float) -> float:
    from scipy.stats import norm  # type: ignore
    return float(norm.isf(alpha / 2))


def add_enrichment_statistics(df, alpha: float, compat: bool = False, groups: "Optional[list[str]]" = None):
    """SE_log_enr, CI_enr_low/high, p_enr, p_adjust_enr (Bonferroni over the table, or within `groups`)."""
    np = _np()
    df = df.copy()
    if len(df) == 0:
        for c in ("SE_log_enr", "CI_enr_low", "CI_enr_high", "p_enr", "p_adjust_enr"):
            df[c] = []
        return df
    enr, lo, hi, se, p = enrichment_stats(df["nVariantsOverlappingEnhancers"], df["nVariantsTotal"],
                                          df["nCommonVariantsOverlappingEnhancers"], df["nCommonVariantsTotal"], alpha, compat)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = df["enrichment"].to_numpy(dtype=float)
        z = _z(alpha)
        df["SE_log_enr"] = se
        df["CI_enr_low"] = np.exp(np.log(e) - z * se)
        df["CI_enr_high"] = np.exp(np.log(e) + z * se)
    df["p_enr"] = p
    if groups:
        df["p_adjust_enr"] = df.groupby(groups, sort=False, dropna=False)["p_enr"].transform(lambda s: bonferroni(s.to_numpy()))
    else:
        df["p_adjust_enr"] = bonferroni(df["p_enr"].to_numpy())
    return df


def add_recall_overlap_statistics(df, alpha: float):
    np = _np()
    df = df.copy()
    z = _z(alpha)
    n = df["nVariantsTotal"].astype(float); x = df["nVariantsOverlappingEnhancers"].astype(float)
    df["recall_adjust"] = (x + 2) / (n + 4)
    df["SE_recall"] = np.sqrt(df["recall_adjust"] * (1 - df["recall_adjust"]) / (n + 4))
    df["CI_recall_low"] = df["recall"] - z * df["SE_recall"]
    df["CI_recall_high"] = df["recall"] + z * df["SE_recall"]
    return df


def _finish_overlap(df, alpha: float, compat: bool, groups: "Optional[list[str]]" = None):
    np = _np()
    with np.errstate(divide="ignore", invalid="ignore"):
        df["recall"] = df["nVariantsOverlappingEnhancers"] / df["nVariantsTotal"]
        df["enrichment"] = df["recall"] / (df["nCommonVariantsOverlappingEnhancers"] / df["nCommonVariantsTotal"])
    df = add_enrichment_statistics(df, alpha, compat, groups)
    return add_recall_overlap_statistics(df, alpha)


def summarize_grouped_enrichment_recall(df, keys: "list[str]", alpha: float, compat: bool = False):
    """summarize_grouped_enrichment_recall: summed counts per `keys`; as in dplyr, the result stays grouped by
    keys[:-1], so the Bonferroni correction runs within those groups."""
    cnt = ["nVariantsTotal", "nVariantsOverlappingEnhancers", "nCommonVariantsTotal", "nCommonVariantsOverlappingEnhancers"]
    g = df.groupby(keys, sort=True, dropna=False)[cnt].sum().reset_index()
    return _finish_overlap(g, alpha, compat, keys[:-1] if len(keys) > 1 else None)


def _tp_cols(kind: str) -> "tuple[str, str]":
    if kind == "enhancers":
        return "nCredibleSetsOverlappingEnhancersCorrectGene", "nCredibleSetsOverlappingEnhancersAnyGene"
    return "nCredibleSetGenePredictionsCorrect", "nCredibleSetGenePredictions"


def add_linking_statistics(df, alpha: float, kind: str):
    """add_recall_linking_statistics + add_precision_linking_statistics."""
    np = _np()
    df = df.copy()
    z = _z(alpha)
    tp, pp = _tp_cols(kind)
    n = df["nCredibleSetsTotal"].astype(float); x = df[tp].astype(float); m = df[pp].astype(float)
    df["recall_adjust"] = (x + 2) / (n + 4)
    df["SE_recall"] = np.sqrt(df["recall_adjust"] * (1 - df["recall_adjust"]) / (n + 4))
    df["CI_recall_low"] = df["recall"] - z * df["SE_recall"]
    df["CI_recall_high"] = df["recall"] + z * df["SE_recall"]
    df["precision_adjust"] = (x + 2) / (m + 4)
    df["SE_precision"] = np.sqrt(df["precision_adjust"] * (1 - df["precision_adjust"]) / (m + 4))
    df["CI_precision_low"] = df["precision"] - z * df["SE_precision"]
    df["CI_precision_high"] = df["precision"] + z * df["SE_precision"]
    return df


def summarize_grouped_precision_recall(df, keys: "list[str]", alpha: float, kind: str):
    np = _np()
    tp, pp = _tp_cols(kind)
    g = df.groupby(keys, sort=True, dropna=False)[["nCredibleSetsTotal", tp, pp]].sum().reset_index()
    with np.errstate(divide="ignore", invalid="ignore"):
        g["recall"] = g[tp] / g["nCredibleSetsTotal"]
        g["precision"] = g[tp] / g[pp]
    return add_linking_statistics(g, alpha, kind)


# ---------------------------------------------------------------------------
# Predictions and per-biosample intersections
# ---------------------------------------------------------------------------

def read_predictions(path: str, score_col: str, universe: "set[str]", invert: bool, biosample: str, chroms: "set[str]"):
    """process_predictions.smk + process_predictions.R."""
    pd = _pd()
    cols = _header_line(Path(path))[0]
    need = ["chr", "start", "end", "TargetGene", score_col]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(f"{path}: missing columns {missing}; available: {cols[:20]}")
    df = read_tsv(Path(path), usecols=need, dtype={"chr": str, "TargetGene": str})
    df = df.rename(columns={score_col: "score"})
    for c in ("score", "start", "end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end", "TargetGene", "score"])
    df = df[df["chr"].isin(chroms)]
    df["start"] = df["start"].astype("int64"); df["end"] = df["end"].astype("int64")
    df = df[df["TargetGene"].isin(universe)]
    if invert:
        df["score"] = -df["score"]
    df["biosample"] = biosample
    return df[["chr", "start", "end", "biosample", "TargetGene", "score"]].reset_index(drop=True)


class BiosampleData:
    """All intersections one prediction file needs: variant locations x predictions (full), per-location and
    per-background-SNP best score, thresholded intervals and enhancer set size."""

    def __init__(self, name: str, pred, locs, bg, threshold: float):
        np = _np()
        self.name = name
        self.n_pred = len(pred)
        thr = pred[pred["score"] >= threshold]
        self.thr_intervals = thr[["chr", "start", "end"]]
        self.bp = merged_bp(self.thr_intervals) if len(thr) else 0
        j = overlap_join(locs[["chr", "start", "end", "loc_id"]], pred[["chr", "start", "end", "TargetGene", "score"]])
        self.loc_int = j[["loc_id", "TargetGene", "score"]].reset_index(drop=True) if len(j) else \
            _pd().DataFrame({"loc_id": np.array([], dtype=np.int64), "TargetGene": [], "score": np.array([], dtype=float)})
        self.loc_best = self.loc_int.groupby("loc_id")["score"].max()
        if bg is not None and len(bg) and len(pred):
            el = pred.groupby(["chr", "start", "end"], sort=False, as_index=False)["score"].max()
            cand = bg[overlaps_union(bg, merge_intervals(el))]
            jb = overlap_join(cand[["chr", "start", "end", "bg_id"]], el)
            self.bg_best = jb.groupby("bg_id")["score"].max() if len(jb) else _pd().Series(dtype=float)
        else:
            self.bg_best = _pd().Series(dtype=float)


def _union_best(series_list):
    pd = _pd()
    s = [x for x in series_list if len(x)]
    if not s:
        return pd.Series(dtype=float)
    return pd.concat(s).groupby(level=0).max()


def expand_biosample(method: Method, name: str, cfg: dict) -> "list[str]":
    """get_single_biosamples."""
    if name in method.biosamples:
        return [name]
    if name == "ALL":
        return list(method.biosamples)
    return list(cfg["biosampleGroups"].get(name, []))


def comparison_pairs(method: Method, comp, cfg: dict) -> "list[tuple[str, str]]":
    """get_comparison_biosamples_per_method: (individual biosample, trait) pairs the comparisons need."""
    out = []
    for b, t in zip(comp["biosample"], comp["trait"]):
        if b in method.biosamples:
            out.append((b, t))
        elif b in method.groups:
            out += [(x, t) for x in expand_biosample(method, b, cfg)]
    return out


def comparison_names(method: Method, comp) -> "list[str]":
    """get_biosample_names_per_method (unique, in table order)."""
    avail = set(method.biosamples) | set(method.groups)
    seen, out = set(), []
    for b in comp["biosample"]:
        if b in avail and b not in seen:
            seen.add(b); out.append(b)
    return out


def quantile_threshold_span(scores, n_steps_target: int, boolean: bool, provided: float, compat: bool = False) -> "list[float]":
    """generate_quantile_threshold_span.R.  `scores` are per-variant maxima (default) or, with compat, the raw
    scores reduced to distinct values as upstream's group_by(-c(predScore)) does."""
    np = _np()
    if boolean:
        n_steps_target = 2
    s = np.asarray(scores, dtype=float)
    s = s[np.isfinite(s)]
    if compat:
        s = np.unique(s)
    if s.size == 0:
        return [float(provided)]
    n_distinct = np.unique(s).size
    even = np.linspace(s.min(), s.max(), n_steps_target)

    def q(n):
        return np.unique(np.concatenate([np.quantile(s, np.linspace(0, 1, n)), [provided]]))
    steps = n_steps_target
    thr = q(steps)
    while thr.size < n_steps_target and steps < n_distinct:
        steps += 1
        thr = q(steps)
    return sorted(set(np.concatenate([thr, even]).tolist()))


def overlap_counts(loc_best, bg_best, threshold: float, loc_trait):
    """Per trait distinct variant locations with best score >= threshold; background SNP count."""
    hit = set(loc_best.index[loc_best.to_numpy() >= threshold])
    sub = loc_trait[loc_trait["loc_id"].isin(hit)]
    counts = sub.groupby("trait").size().rename("nVariantsOverlappingEnhancers").reset_index()
    return counts, int((bg_best.to_numpy() >= threshold).sum())


# ---------------------------------------------------------------------------
# Gene linking (evaluate_gene_linking.R, evaluate_baseline_predictors.R)
# ---------------------------------------------------------------------------

def load_gene_prioritization(path: str, universe: "Optional[set[str]]"):
    pd = _pd()
    gp = pd.read_csv(path, sep="\t", low_memory=False)
    need = ["CredibleSet", "Disease", "TargetGene", "truth", "POPS.Score", "PromoterDistanceToBestSNP"]
    miss = [c for c in need if c not in gp.columns and not (c == "Disease" and "trait" in gp.columns)]
    if miss:
        raise SystemExit(f"{path}: gene prioritisation table lacks {miss}")
    gp = gp.rename(columns={"Disease": "trait"})
    gp["truth"] = gp["truth"].map(_truthy) if gp["truth"].dtype != bool else gp["truth"]
    gp["CredibleSet"] = gp["CredibleSet"].astype(str); gp["trait"] = gp["trait"].astype(str); gp["TargetGene"] = gp["TargetGene"].astype(str)
    if universe is not None:
        gp = gp[gp["TargetGene"].isin(universe)]
    return gp.reset_index(drop=True)


def _truth_sets(gp):
    return gp.loc[gp["truth"], ["CredibleSet", "trait", "TargetGene"]].rename(columns={"TargetGene": "TruthGene"}).drop_duplicates()


def _top_by(gp, col: str, n: int, ascending: bool):
    """arrange + group_by(CredibleSet, trait) + row_number() <= n (stable: file order breaks ties; NA last)."""
    s = gp.sort_values(col, ascending=ascending, kind="mergesort", na_position="last")
    r = s.groupby(["CredibleSet", "trait"], sort=False).cumcount() + 1
    return s[r <= n]


def gene_linking_counts(rows, cs_truth, pops_genes, num_pred: int):
    """One trait x biosample: (noPoPS counts, intPoPS counts) from rows CredibleSet, TargetGene, score."""
    pd = _pd()
    sub = rows[rows["CredibleSet"].isin(set(cs_truth["CredibleSet"]))]
    if len(sub):
        mx = sub.groupby(["CredibleSet", "TargetGene"], as_index=False)["score"].max().rename(columns={"score": "maxPredScore"})
        mx["predRank"] = mx.groupby("CredibleSet")["maxPredScore"].rank(ascending=False, method="min")
        mx = mx[mx["predRank"] <= num_pred]
        pred_only = mx.merge(cs_truth, on="CredibleSet", how="left")
    else:
        pred_only = pd.DataFrame(columns=["CredibleSet", "TargetGene", "maxPredScore", "predRank", "trait", "TruthGene"])

    def counts(df):
        a = df[["CredibleSet", "TruthGene"]].drop_duplicates().shape[0]
        b = df[["CredibleSet", "TargetGene"]].drop_duplicates().shape[0]
        c = df.loc[df["TargetGene"] == df["TruthGene"], ["CredibleSet", "TargetGene", "TruthGene"]].drop_duplicates().shape[0]
        return a, b, c
    pp = pred_only.merge(pops_genes, on=["CredibleSet", "trait"], how="left")
    pp = pp[pp["TargetGene"] == pp["PoPSGene"]]
    return counts(pred_only), counts(pp)


def evaluate_baseline(gp, num_pops: int, alpha: float):
    """evaluate_baseline_predictors.R -> precisionRecall.byTrait (PoPS, distanceToTSS)."""
    pd = _pd()
    cs_total = gp.groupby("trait")["CredibleSet"].nunique().rename("nCredibleSetsTotal").reset_index()
    truth = _truth_sets(gp)
    pops = _top_by(gp[["CredibleSet", "trait", "TargetGene", "POPS.Score"]], "POPS.Score", num_pops, False).merge(truth, on=["CredibleSet", "trait"], how="left")
    dist = gp[["CredibleSet", "trait", "TargetGene", "PromoterDistanceToBestSNP"]].sort_values(
        ["trait", "CredibleSet", "PromoterDistanceToBestSNP"], kind="mergesort", na_position="last")
    dist = dist[dist.groupby(["CredibleSet", "trait"], sort=False).cumcount() < num_pops]
    dist = dist.drop(columns="PromoterDistanceToBestSNP").merge(truth, on=["CredibleSet", "trait"], how="left")
    rows = []
    for name, df in (("PoPS", pops), ("distanceToTSS", dist)):
        for t in gp["trait"].unique():
            d = df[df["trait"] == t]
            rows.append({"trait": t, "method": name,
                         "nCredibleSetsWithPredictions": d[["CredibleSet", "TruthGene"]].drop_duplicates().shape[0],
                         "nCredibleSetGenePredictions": d[["CredibleSet", "TargetGene"]].drop_duplicates().shape[0],
                         "nCredibleSetGenePredictionsCorrect": d[d["TargetGene"] == d["TruthGene"]].drop_duplicates().shape[0]})
    res = pd.DataFrame(rows).merge(cs_total, on="trait", how="left")
    np = _np()
    with np.errstate(divide="ignore", invalid="ignore"):
        res["recall"] = res["nCredibleSetGenePredictionsCorrect"] / res["nCredibleSetsTotal"]
        res["precision"] = res["nCredibleSetGenePredictionsCorrect"] / res["nCredibleSetGenePredictions"]
    return add_linking_statistics(res, alpha, "baseline")


def cmd_baseline(args: argparse.Namespace) -> int:
    cfg = load_config(args)
    if not cfg.get("genePrioritizationTable") or not Path(cfg["genePrioritizationTable"]).is_file():
        raise SystemExit(f"gene prioritisation table not found: {cfg.get('genePrioritizationTable')}")
    universe = load_gene_universe(cfg["TSS"]) if cfg.get("TSS") and Path(cfg["TSS"]).is_file() else None
    gp = load_gene_prioritization(cfg["genePrioritizationTable"], universe)
    res = evaluate_baseline(gp, cfg["numPoPSGenes"], cfg["thresholdPval"])
    out = run_dir(args.label)
    write_tsv(res, out / "baseline" / "gene_linking" / "precisionRecall.byTrait.tsv.gz")
    agg = summarize_grouped_precision_recall(res, ["method"], cfg["thresholdPval"], "baseline")
    write_tsv(agg, out / "baseline" / "gene_linking" / "precisionRecall.allTraits.tsv")
    summ = {"label": args.label, "numPoPSGenes": cfg["numPoPSGenes"], "gene_universe_filter": universe is not None,
            "n_traits": int(gp["trait"].nunique()), "n_credible_sets": int(gp[["CredibleSet", "trait"]].drop_duplicates().shape[0]),
            "all_traits": {r["method"]: {"precision": r["precision"], "recall": r["recall"], "CI_precision": [r["CI_precision_low"], r["CI_precision_high"]],
                                         "CI_recall": [r["CI_recall_low"], r["CI_recall_high"]]} for r in agg.to_dict("records")}}
    (out / "summary.json").write_text(json.dumps(summ, indent=2, default=float))
    rep = [f"# Baseline gene linking — {args.label}", "",
           f"Top {cfg['numPoPSGenes']} genes per credible set by PoPS score and by promoter distance to the best SNP, "
           f"against the silver-standard truth genes ({summ['n_credible_sets']} credible sets, {summ['n_traits']} traits"
           f"{', gene universe applied' if universe is not None else ''}).", "",
           md_table(["method", "precision", "recall", "CS total", "predictions", "correct"],
                    [(r["method"], _fmt(r["precision"]), _fmt(r["recall"]), r["nCredibleSetsTotal"], r["nCredibleSetGenePredictions"],
                      r["nCredibleSetGenePredictionsCorrect"]) for r in agg.to_dict("records")])]
    (out / "report.md").write_text("\n".join(rep) + "\n")
    print(f"JSON: {out / 'summary.json'}")
    print(f"Report: {out / 'report.md'}")
    return 0


# ---------------------------------------------------------------------------
# Per-method workflow
# ---------------------------------------------------------------------------

def per_method(method: Method, cfg: dict, comp, V, locs, loc_trait, bg, n_bg_total: int, nvt, universe: "set[str]",
               chroms: "set[str]", gp, gp_raw, out: Path, compat: bool) -> dict:
    pd, np = _pd(), _np()
    alpha = cfg["thresholdPval"]
    trait_groups = set(cfg["traitGroups"]) | {"ALL"}
    mdir = out / method.name
    data: "dict[str, BiosampleData]" = {}
    for b in method.biosamples:
        pred = read_predictions(method.files[b], method.score_col, universe, method.inverse, b, chroms)
        data[b] = BiosampleData(b, pred, locs, bg, method.threshold)
        log.info("%s %s: %d predictions, %d bp thresholded", method.name, b, data[b].n_pred, data[b].bp)

    def combined(name):
        mem = expand_biosample(method, name, cfg)
        return mem, _union_best([data[m].loc_best for m in mem]), _union_best([data[m].bg_best for m in mem])

    # background overlap + enhancer set size, thresholded enrichment / recall
    stats_rows, thr_rows = [], []
    for name in method.biosamples + method.groups:
        mem, lb, bb = combined(name)
        if name in method.biosamples:
            bp = data[name].bp
        else:
            iv = pd.concat([data[m].thr_intervals for m in mem], ignore_index=True)
            bp = merged_bp(iv) if len(iv) else 0
        counts, nbg = overlap_counts(lb, bb, method.threshold, loc_trait)
        stats_rows.append({"biosample": name, "bp_enhancers": int(bp), "bg_variant_count": nbg})
        counts["biosample"] = name
        thr_rows.append(counts)
    stats = pd.DataFrame(stats_rows)
    write_tsv(stats, mdir / "variant_overlap" / "bgOverlap.enhancerSetSize.thresholdedPredictions.tsv", quiet=True)
    thr = pd.concat(thr_rows, ignore_index=True) if thr_rows else pd.DataFrame(columns=["trait", "nVariantsOverlappingEnhancers", "biosample"])
    thr = thr.merge(nvt, on="trait", how="left").merge(
        stats.rename(columns={"bp_enhancers": "bpEnhancers", "bg_variant_count": "nCommonVariantsOverlappingEnhancers"}), on="biosample", how="left")
    thr["nCommonVariantsTotal"] = n_bg_total
    with np.errstate(divide="ignore", invalid="ignore"):
        thr["recall"] = thr["nVariantsOverlappingEnhancers"] / thr["nVariantsTotal"]
        thr["enrichment"] = thr["recall"] / (thr["nCommonVariantsOverlappingEnhancers"] / thr["nCommonVariantsTotal"])
    thr["method"] = method.name
    thr["biosampleGroup"] = thr["biosample"].isin(set(method.groups))
    thr["traitGroup"] = thr["trait"].isin(trait_groups)
    thr = add_recall_overlap_statistics(add_enrichment_statistics(thr, alpha, compat), alpha)
    write_tsv(thr, mdir / "variant_overlap" / "enrichmentRecall.thresholded.traitByBiosample.tsv.gz", quiet=True)

    # threshold span over the comparisons' (biosample, trait) pairs, full intersections
    pairs = comparison_pairs(method, comp, cfg)
    sc = []
    for b in sorted({p[0] for p in pairs}):
        tr = {t for bb, t in pairs if bb == b}
        ids = set(loc_trait.loc[loc_trait["trait"].isin(tr), "loc_id"])
        li = data[b].loc_int
        sc.append(li[li["loc_id"].isin(ids)][["loc_id", "score"]])
    sc = pd.concat(sc, ignore_index=True) if sc else pd.DataFrame(columns=["loc_id", "score"])
    scores = sc["score"].to_numpy(dtype=float) if compat else sc.groupby("loc_id")["score"].max().to_numpy(dtype=float)
    span = quantile_threshold_span(scores, cfg["nThresholdSteps"], method.boolean, method.threshold, compat)
    write_tsv(pd.DataFrame({"threshold": span}), mdir / "thresholdSpan.tsv", quiet=True)

    # enrichment-recall curves per comparison biosample name
    curves = {}
    span_arr = np.asarray(span, dtype=float)
    for name in comparison_names(method, comp):
        mem, lb, bb = combined(name)
        m = loc_trait.merge(pd.DataFrame({"loc_id": lb.index.to_numpy(), "best": lb.to_numpy(dtype=float)}), on="loc_id", how="inner")
        bgs = np.sort(bb.to_numpy(dtype=float))
        bg_counts = bgs.size - np.searchsorted(bgs, span_arr, side="left")
        rows = []
        for t, sub in m.groupby("trait", sort=True):
            s = np.sort(sub["best"].to_numpy(dtype=float))
            c = s.size - np.searchsorted(s, span_arr, side="left")
            for thv, cv, bgc in zip(span_arr, c, bg_counts):
                if cv > 0:
                    rows.append((t, int(cv), int(bgc), float(thv)))
        cdf = pd.DataFrame(rows, columns=["trait", "nVariantsOverlappingEnhancers", "nCommonVariantsOverlappingEnhancers", "threshold"])
        cdf = cdf.sort_values(["threshold", "trait"], kind="mergesort").reset_index(drop=True)
        cdf = cdf.merge(nvt, on="trait", how="left")
        cdf["nCommonVariantsTotal"] = n_bg_total
        with np.errstate(divide="ignore", invalid="ignore"):
            cdf["recall"] = cdf["nVariantsOverlappingEnhancers"] / cdf["nVariantsTotal"]
            cdf["enrichment"] = cdf["recall"] / (cdf["nCommonVariantsOverlappingEnhancers"] / cdf["nCommonVariantsTotal"])
        cdf["biosample"] = name
        cdf["method"] = method.name
        cdf["biosampleGroup"] = len(mem) > 1
        cdf["traitGroup"] = cdf["trait"].isin(trait_groups)
        cdf = add_recall_overlap_statistics(add_enrichment_statistics(cdf, alpha, compat), alpha)
        write_tsv(cdf, mdir / "variant_overlap" / f"enrichmentRecallAcrossThresholds.forBiosample{name}.tsv.gz", quiet=True)
        curves[name] = cdf

    # gene linking
    link = pd.DataFrame()
    if gp is not None:
        gpl = gp_raw if compat else gp
        cs_total = gpl.groupby("trait")["CredibleSet"].nunique().rename("nCredibleSetsTotal").reset_index()
        truth = _truth_sets(gpl)
        pops = _top_by(gpl[["CredibleSet", "trait", "TargetGene", "POPS.Score"]], "POPS.Score", cfg["numPoPSGenes"], False)
        pops = pops[["CredibleSet", "trait", "TargetGene"]].rename(columns={"TargetGene": "PoPSGene"})
        gp_traits = set(gpl["trait"])
        vl = V[V["trait"].isin(gp_traits)][["loc_id", "CredibleSet", "trait"]]
        rows = []
        for name in method.biosamples + method.groups:
            is_group = name not in method.biosamples
            mem = expand_biosample(method, name, cfg)
            parts = []
            for mb in mem:
                li = data[mb].loc_int
                if is_group or not compat:
                    li = li[li["score"] >= method.threshold]
                parts.append(li)
            li = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["loc_id", "TargetGene", "score"])
            vi = vl.merge(li, on="loc_id", how="inner")
            for t, sub in vi.groupby("trait", sort=True):
                cs_t = truth[truth["trait"] == t]
                (a0, b0, c0), (a1, b1, c1) = gene_linking_counts(sub, cs_t, pops[pops["trait"] == t], cfg["numPredGenes"])
                for ip, (a, b_, c) in ((False, (a0, b0, c0)), (True, (a1, b1, c1))):
                    rows.append({"trait": t, "biosample": name, "group": is_group, "intersectPoPS": ip,
                                 "nCredibleSetsOverlappingEnhancers": a, "nCredibleSetsOverlappingEnhancersAnyGene": b_,
                                 "nCredibleSetsOverlappingEnhancersCorrectGene": c})
        cols = ["trait", "biosample", "group", "intersectPoPS", "nCredibleSetsOverlappingEnhancers",
                "nCredibleSetsOverlappingEnhancersAnyGene", "nCredibleSetsOverlappingEnhancersCorrectGene"]
        link = pd.DataFrame(rows, columns=cols).merge(cs_total, on="trait", how="left")
        link["nCredibleSetsTotal"] = link["nCredibleSetsTotal"].fillna(0).astype(int)
        with np.errstate(divide="ignore", invalid="ignore"):
            link["recall"] = (link["nCredibleSetsOverlappingEnhancersCorrectGene"] / link["nCredibleSetsTotal"]).fillna(0)
            link["precision"] = (link["nCredibleSetsOverlappingEnhancersCorrectGene"] / link["nCredibleSetsOverlappingEnhancersAnyGene"]).fillna(0)
        link = add_linking_statistics(link, alpha, "enhancers").fillna(0)
        link["method"] = method.name
        write_tsv(link, mdir / "gene_linking" / "precisionRecall.traitByBiosample.tsv.gz", quiet=True)
    print(f"Wrote: {mdir}")
    return {"thresholded": thr, "stats": stats, "span": span, "curves": curves, "linking": link}


# ---------------------------------------------------------------------------
# Visualisation-stage tables (colour palette, metric ranges, comparisons)
# ---------------------------------------------------------------------------

def color_palette(methods: "list[dict]"):
    pd = _pd()
    rows, k = [], 0
    for m in methods:
        c = m.get("color") or ""
        if not c:
            c = FILL_COLORS[k % len(FILL_COLORS)]; k += 1
        rows.append({"method": m["method"], "pred_name_long": m["pred_name_long"], "hex": c})
    rows += [{"method": a, "pred_name_long": b, "hex": c} for a, b, c in BASELINE_COLORS]
    return pd.DataFrame(rows)


def metric_ranges(thr_tables: "dict[str, Any]", link_tables: "dict[str, Any]", alpha: float):
    """gather_metric_ranges.R."""
    pd, np = _pd(), _np()
    rows = []

    def mm(method, metric, vals):
        vals = np.asarray(vals, dtype=float)
        for lab, f in (("min", np.min), ("max", np.max)):
            rows.append({"method": method, "metric": metric, "min_max": lab, "value": float(f(vals)) if vals.size else float("nan")})
    for m, df in thr_tables.items():
        d = df[np.isfinite(df["enrichment"]) & np.isfinite(df["recall"]) & (df["nVariantsTotal"] > 25)
               & (df["nCommonVariantsOverlappingEnhancers"] / df["nCommonVariantsTotal"] > 0.001)]
        mm(m, "enrichment_overlap", d.loc[d["p_adjust_enr"] < alpha, "enrichment"])
        mm(m, "recall_overlap", d["recall"])
        lk = link_tables.get(m)
        if lk is not None and len(lk):
            lk = lk[np.isfinite(lk["precision"]) & np.isfinite(lk["recall"]) & (lk["nCredibleSetsTotal"] > 0)]
            for metric in ("precision", "recall"):
                mm(m, f"{metric}_linking_intersectPoPS", lk.loc[lk["intersectPoPS"].astype(bool), metric])
                mm(m, f"{metric}_linking", lk.loc[~lk["intersectPoPS"].astype(bool), metric])
    return pd.DataFrame(rows, columns=["method", "metric", "min_max", "value"])


def expand_traits(traits: "Iterable[str]", trait_groups: "dict[str, list[str]]") -> "list[str]":
    out = []
    for t in traits:
        out += trait_groups[t] if t in trait_groups else [t]
    return list(dict.fromkeys(x for x in out if x))


def n_equals_variant_overlap(df, V, traits: "Iterable[str]", trait_groups: dict) -> str:
    tt = expand_traits(traits, trait_groups)
    n_var = V[V["trait"].isin(tt)][["chr", "start", "end"]].drop_duplicates().shape[0]
    n_pairs = int(df[["biosample", "trait", "nVariantsTotal"]].drop_duplicates()["nVariantsTotal"].sum()) if len(df) else 0
    return f"{n_var} variants from {len(tt)} traits\n({n_pairs} variant-biosample pairs)"


def n_equals_gene_linking(df) -> str:
    if not len(df):
        return "0 credible sets"
    cs = df[["trait", "nCredibleSetsTotal"]].drop_duplicates()
    pairs = df[["biosample", "trait", "nCredibleSetsTotal"]].drop_duplicates()
    return (f"{int(cs['nCredibleSetsTotal'].sum())} credible sets from {cs['trait'].nunique()} traits tested\n"
            f"({int(pairs['nCredibleSetsTotal'].sum())} credible set-biosample pairs)")


def curve_plot_filter(df, cp, alpha: float):
    """add_plotting_params (plot_enrichment_recall_curves.R)."""
    d = df[(df["recall"] > 0.005) & (df["p_adjust_enr"] < alpha) & df["enrichment"].notna() & (df["nVariantsTotal"] > 25)
           & (df["nCommonVariantsOverlappingEnhancers"] / df["nCommonVariantsTotal"] > 0.001)]
    d = d.merge(cp, on="method", how="left")
    n = d.groupby(["method", "biosample", "trait"]).size().rename("n_points").reset_index()
    return d.merge(n, on=["method", "biosample", "trait"], how="left")


def comparison_curves(name: str, comp, curves_by_method: "dict[str, dict]", cp, alpha: float, compat: bool):
    pd = _pd()
    sub = comp[comp["name"] == name]
    keys = set(sub["biosample"] + "." + sub["trait"])
    parts = []
    for m, cur in curves_by_method.items():
        for b, df in cur.items():
            if b in set(sub["biosample"]) and len(df):
                parts.append(df)
    if not parts:
        return pd.DataFrame(), pd.DataFrame()
    sep = pd.concat(parts, ignore_index=True)
    sep["key"] = sep["biosample"] + "." + sep["trait"]
    sep = sep[sep["key"].isin(keys)].drop(columns=["biosampleGroup", "traitGroup"])
    if not len(sep):
        return pd.DataFrame(), pd.DataFrame()
    agg = summarize_grouped_enrichment_recall(sep, ["threshold", "method"], alpha, compat)
    agg["trait"] = "all_matched"; agg["biosample"] = "all_matched"; agg["key"] = "all_matched.all_matched"
    return sep, agg


def filter_linking_rows(df, sub, trait_groups: dict, biosample_groups: dict, kind: str, compat: bool):
    """filter_linking_data (plot_thresholded_performance_comparison.R)."""
    pd = _pd()
    parts = []
    for b, t in zip(sub["biosample"], sub["trait"]):
        d = df
        if kind == "enhancers":
            if b == "ALL":
                d = df if compat else df[~df["group"].astype(bool)]
            elif b in biosample_groups:
                d = df[df["biosample"].isin(biosample_groups[b])]
            else:
                d = df[df["biosample"] == b]
        d = d[d["trait"].isin(trait_groups[t])] if t in trait_groups else d[d["trait"] == t]
        parts.append(d)
    return pd.concat(parts, ignore_index=True).drop_duplicates() if parts else df.iloc[0:0]


def performance_comparison(name: str, comp, thr_tables, link_tables, baseline, cp, trait_groups, biosample_groups,
                           alpha: float, compat: bool):
    pd = _pd()
    sub = comp[comp["name"] == name]
    keys = set(sub["biosample"] + "." + sub["trait"])
    ov = pd.concat([df for df in thr_tables.values()], ignore_index=True) if thr_tables else pd.DataFrame()
    ov_in = pd.DataFrame()
    ov_s = pd.DataFrame()
    if len(ov):
        ov_in = ov[(ov["biosample"] + "." + ov["trait"]).isin(keys)].merge(cp, on="method", how="left")
        if len(ov_in):
            ov_s = summarize_grouped_enrichment_recall(ov_in, ["method", "pred_name_long", "hex"], alpha, compat)
    lk_in = pd.DataFrame(); lk_s = pd.DataFrame()
    lk = pd.concat([df for df in link_tables.values() if df is not None and len(df)], ignore_index=True) \
        if any(df is not None and len(df) for df in link_tables.values()) else pd.DataFrame()
    if len(lk):
        lk_in = filter_linking_rows(lk, sub, trait_groups, biosample_groups, "enhancers", compat).merge(cp, on="method", how="left")
        if len(lk_in):
            lk_s = summarize_grouped_precision_recall(lk_in, ["method", "pred_name_long", "intersectPoPS"], alpha, "enhancers")
            lk_s = lk_s.drop(columns=list(_tp_cols("enhancers")))
    if baseline is not None and len(baseline):
        bl_in = filter_linking_rows(baseline, sub, trait_groups, biosample_groups, "baseline", compat).merge(cp, on="method", how="left")
        if len(bl_in):
            bl_s = summarize_grouped_precision_recall(bl_in, ["method", "pred_name_long"], alpha, "baseline")
            bl_s = bl_s.drop(columns=list(_tp_cols("baseline")))
            lk_s = pd.concat([lk_s, bl_s], ignore_index=True)
    return ov_in, ov_s, lk_in, lk_s


def cluster_order_upstream(M) -> "tuple[list, list]":
    """cluster_traits_biosamples: hclust(dist(1 - cor(M)), 'ward.D2') for traits (columns) and biosamples (rows)."""
    np = _np()
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage  # type: ignore
        from scipy.spatial.distance import pdist  # type: ignore
    except Exception:
        return list(M.index), list(M.columns)

    def order(X, labels):
        # X: variables in columns -> cor between columns
        if len(labels) < 3:
            return list(labels)
        with np.errstate(invalid="ignore", divide="ignore"):
            C = np.corrcoef(X, rowvar=False)
        D = pdist(1 - C)
        D = np.nan_to_num(D, nan=0.0)
        if D.sum() <= 0:
            return list(labels)
        return [labels[i] for i in leaves_list(linkage(D, method="ward"))]
    X = M.to_numpy(dtype=float)
    return order(X.T, list(M.index)), order(X, list(M.columns))


def heatmap_matrix(df, metric: str, keep_na: bool = False):
    """Clustered biosample x trait matrix (NA -> 0 for clustering; keep_na returns the NA cells for drawing)."""
    raw = df.pivot_table(index="biosample", columns="trait", values=metric, aggfunc="first")
    M = raw.fillna(0.0)
    ro, co = cluster_order_upstream(M)
    return raw.loc[ro, co] if keep_na else M.loc[ro, co]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _cmap(colors):
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
    return LinearSegmentedColormap.from_list("u", colors)


def _draw_heat(ax, M, colors, vmax, title, sig=None):
    np = _np()
    cm = _cmap(colors)
    cm.set_bad("#ffffff")
    im = ax.imshow(np.ma.masked_invalid(M.to_numpy(dtype=float)), aspect="auto", cmap=cm, vmin=0, vmax=max(vmax, 1e-9), interpolation="nearest")
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=6)
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, fontsize=6, rotation=60, ha="right")
    if sig is not None:
        yy, xx = np.nonzero(sig.loc[M.index, M.columns].to_numpy())
        ax.scatter(xx, yy, marker="*", s=6, color="white", linewidths=0)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK)
    import matplotlib.pyplot as plt  # type: ignore
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02).ax.tick_params(labelsize=6)


def plot_method_heatmaps(method: str, thr, link, ranges, cfg: dict, pdir: Path, figures: "list[str]", mats_out: Path) -> None:
    plt = _plt()
    np = _np()
    alpha = cfg["thresholdPval"]
    res = thr[np.isfinite(thr["enrichment"]) & np.isfinite(thr["recall"]) & (thr["nVariantsTotal"] > 20)].copy()
    if len(res):
        if cfg["plotFixedScale"] and len(ranges):
            e = ranges[(ranges["metric"] == "enrichment_overlap") & (ranges["min_max"] == "max")]["value"].dropna()
            r = ranges[(ranges["metric"] == "recall_overlap") & (ranges["min_max"] == "max")]["value"].dropna()
            max_enr = float(np.quantile(e, 0.75)) if len(e) else 1.0
            max_rec = float(r.max()) if len(r) else 1.0
        else:
            t = res[(res["nVariantsTotal"] > 25) & (res["nCommonVariantsOverlappingEnhancers"] / res["nCommonVariantsTotal"] > 0.001)]
            max_enr = float(t["enrichment"].max()) if len(t) else float(res["enrichment"].max())
            max_rec = float(t["recall"].max()) if len(t) else float(res["recall"].max())
        ME = heatmap_matrix(res, "enrichment", True); MR = heatmap_matrix(res, "recall", True)
        write_tsv(ME.fillna(0.0).reset_index(), mats_out / f"{method}.enrichmentOverlapHeatmap.matrix.tsv", quiet=True)
        write_tsv(MR.fillna(0.0).reset_index(), mats_out / f"{method}.recallOverlapHeatmap.matrix.tsv", quiet=True)
        if plt is not None:
            sig = res.assign(s=(res["p_adjust_enr"] < alpha).astype(float)).pivot_table(index="biosample", columns="trait", values="s", aggfunc="first")
            sig = sig.reindex(index=ME.index, columns=ME.columns).fillna(0.0) > 0.5
            h = max(3.0, 0.22 * ME.shape[0] + 1.5)
            fig, axes = plt.subplots(1, 2, figsize=(max(9, 0.5 * ME.shape[1] + 5), h))
            _draw_heat(axes[0], ME, ENR_OVL_COLORS, max_enr, f"{method}: enrichment (* p_adj < {alpha})", sig)
            _draw_heat(axes[1], MR, RECALL_OVL_COLORS, max_rec, f"{method}: recall (variant overlap)")
            fig.tight_layout()
            _save(fig, pdir / "overlapHeatmaps" / f"{method}.overlapHeatmaps.png", figures)
    if link is not None and len(link):
        lk = link[np.isfinite(link["precision"]) & np.isfinite(link["recall"]) & (link["nCredibleSetsTotal"] > 0)].copy()
        lk["intersectPoPS"] = lk["intersectPoPS"].astype(bool)
        if not len(lk):
            return
        mats = {}
        for metric in ("precision", "recall"):
            for ip in (False, True):
                d = lk[lk["intersectPoPS"] == ip]
                if not len(d):
                    continue
                if cfg["plotFixedScale"] and len(ranges):
                    mname = f"{metric}_linking" + ("_intersectPoPS" if ip else "")
                    v = ranges[(ranges["metric"] == mname) & (ranges["min_max"] == "max")]["value"].dropna()
                    vmax = float(v.max()) if len(v) else 1.0
                else:
                    vmax = float(d[metric].max())
                M = heatmap_matrix(d, metric, True)
                tag = ("intPoPS." if ip else "") + f"{metric}Linking"
                write_tsv(M.fillna(0.0).reset_index(), mats_out / f"{method}.{tag}Heatmap.matrix.tsv", quiet=True)
                mats[(metric, ip)] = (M, vmax)
        if plt is not None and mats:
            n0 = max(M.shape[0] for M, _ in mats.values()); n1 = max(M.shape[1] for M, _ in mats.values())
            fig, axes = plt.subplots(2, 2, figsize=(max(9, 0.5 * n1 + 5), max(5, 0.44 * n0 + 3)))
            for i, metric in enumerate(("precision", "recall")):
                for j, ip in enumerate((False, True)):
                    ax = axes[i][j]
                    if (metric, ip) not in mats:
                        ax.axis("off"); continue
                    M, vmax = mats[(metric, ip)]
                    _draw_heat(ax, M, PRECISION_COLORS if metric == "precision" else RECALL_LINK_COLORS, vmax,
                               f"{method}: {metric} linking" + (" (int. PoPS)" if ip else ""))
            fig.tight_layout()
            _save(fig, pdir / "linkingHeatmaps" / f"{method}.linkingHeatmaps.png", figures)


def _curve_panel(ax, df, cpd: dict, title: str, xlab: str) -> None:
    np = _np()
    if not len(df):
        ax.text(0.5, 0.5, "no points pass the filters", ha="center", va="center", fontsize=7, color=INK2, transform=ax.transAxes)
        _style(ax, title, xlab, "Enrichment (GWAS vs common variants)")
        return
    max_recall = df["recall"].max() if df["recall"].max() > 0 else 1
    abs_max = df["enrichment"].max()
    if abs_max <= 0:
        max_enr = 10
    else:
        q90 = float(np.quantile(df["enrichment"].dropna(), 0.9))
        max_enr = q90 if abs_max - q90 > 10 else abs_max
    for m, d in df.groupby("method", sort=False):
        col = cpd.get(m, INK2)
        d = d.sort_values("threshold")
        lab = str(d["pred_name_long"].iloc[0])
        line = d[d["n_points"] > 2]
        if len(line):
            ax.plot(line["recall"], line["enrichment"], color=col, lw=1.2, label=lab)
        binary = d[(d["n_points"] <= 2) & (d["threshold"] != 0)]
        if len(binary):
            ax.scatter(binary["recall"], binary["enrichment"], color=col, s=18, label=None if len(line) else lab)
        at = d[np.isclose(d["threshold"], d["score_threshold"])]
        if len(at):
            ax.scatter(at["recall"], at["enrichment"], color=col, s=22, zorder=3)
            ax.vlines(at["recall"], at["CI_enr_low"], at["CI_enr_high"], color=col, lw=1)
    ax.set_xlim(0, max_recall * 1.02); ax.set_ylim(0, max_enr * 1.05)
    _style(ax, title, xlab, "Enrichment (GWAS vs common variants)")


def plot_comparison_curves(name: str, sep, agg, cp, V, comp, trait_groups, pdir: Path, figures: "list[str]") -> None:
    plt = _plt()
    if plt is None:
        return
    cpd = dict(zip(cp["method"], cp["hex"]))
    sub = comp[comp["name"] == name]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    lab = n_equals_variant_overlap(sep, V, sub["trait"].unique(), trait_groups)
    _curve_panel(ax, agg, cpd, f"{name}: all matched pairs", "Recall (fraction of variants in predicted enhancers)\n" + lab)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, pdir / "enrichmentRecallCurves" / f"{name}.combinedCurves.png", figures)
    keys = list(dict.fromkeys(sub["biosample"] + "." + sub["trait"]))
    n = len(keys)
    ncol = 1 if n == 1 else 2
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), squeeze=False)
    for i, k in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        d = sep[sep["key"] == k]
        b, t = k.split(".", 1)
        _curve_panel(ax, d, cpd, f"Biosample: {b}\nTrait: {t}", "Recall\n" + n_equals_variant_overlap(d, V, [t], trait_groups))
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.tight_layout()
    _save(fig, pdir / "enrichmentRecallCurves" / f"{name}.individualCurves.png", figures)


def plot_performance(name: str, ov_in, ov_s, lk_in, lk_s, V, trait_groups, pdir: Path, figures: "list[str]", cp=None) -> None:
    plt = _plt()
    np = _np()
    hexes = dict(zip(cp["method"], cp["hex"])) if cp is not None else {}
    if plt is None or (not len(ov_s) and not len(lk_s)):
        return
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    ax = axes[0]
    if len(ov_s):
        vals = np.sort(np.concatenate([ov_s["CI_enr_high"].to_numpy(dtype=float), ov_s["enrichment"].to_numpy(dtype=float)]))[::-1]
        vals = vals[np.isfinite(vals)]
        emax = (vals[1] * 1.05 if vals.size > 1 and vals[0] > 1.2 * vals[1] else (vals[0] if vals.size else 1))
        for r in ov_s.to_dict("records"):
            c = r.get("hex") or INK2
            ax.hlines(r["enrichment"], max(r["CI_recall_low"], 0), r["CI_recall_high"], color=c, lw=1)
            ax.vlines(r["recall"], max(r["CI_enr_low"], 0), r["CI_enr_high"], color=c, lw=1)
            ax.scatter([r["recall"]], [r["enrichment"]], color=c, s=40, label=r["pred_name_long"])
        ax.set_xlim(0, max(ov_s["CI_recall_high"].max(), 1e-3)); ax.set_ylim(0, emax)
    lab = n_equals_variant_overlap(ov_in, V, ov_in["trait"].unique() if len(ov_in) else [], trait_groups)
    _style(ax, "Variant overlap", "Recall (fraction of variants in predicted enhancers)\n" + lab, "Enrichment\n(GWAS vs common variants)")
    en = lk_s[lk_s["intersectPoPS"].notna()] if len(lk_s) and "intersectPoPS" in lk_s.columns else lk_s.iloc[0:0]
    for k, withp in ((1, False), (2, True)):
        ax = axes[k]
        if len(en):
            x_max = min(1, en["CI_recall_high"].max()); y_max = min(1, en["CI_precision_high"].max())
            ip = en["intersectPoPS"].astype(bool)
            other, cur = en[ip != withp], en[ip == withp]
            for r in other.to_dict("records"):
                ax.scatter([r["recall"]], [r["precision"]], color=hexes.get(r["method"], INK2), s=40, alpha=0.25)
            for r in cur.to_dict("records"):
                c = hexes.get(r["method"], INK2)
                ax.scatter([r["recall"]], [r["precision"]], color=c, s=40, label=r["pred_name_long"])
                ax.hlines(r["precision"], max(r["CI_recall_low"], 0), r["CI_recall_high"], color=c, lw=1)
                ax.vlines(r["recall"], max(r["CI_precision_low"], 0), min(r["CI_precision_high"], 1), color=c, lw=1)
            ax.set_xlim(0, max(x_max, 1e-3)); ax.set_ylim(0, max(y_max, 1e-3))
        _style(ax, "Linking variants to known genes\n" + ("(E2G method int. PoPS)" if withp else "(E2G method)"),
               "Recall (credible sets linked to target gene)\n" + n_equals_gene_linking(lk_in),
               "Precision (correct credible set-gene links)")
    h, l_ = axes[0].get_legend_handles_labels()
    if h:
        axes[2].legend(h, l_, fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    _save(fig, pdir / "thresholdedPerformanceComparison" / f"{name}.scatter.png", figures)


# ---------------------------------------------------------------------------
# Visualisation stage (shared by run and plot)
# ---------------------------------------------------------------------------

def visualize(out: Path, cfg: dict, methods: "list[dict]", comp, V, thr_tables, link_tables, curves_by_method, baseline,
              compat: bool, no_plots: bool) -> dict:
    pd = _pd()
    alpha = cfg["thresholdPval"]
    pdir = out / "plots"
    figures: "list[str]" = []
    trait_groups = dict(cfg["traitGroups"])
    trait_groups.setdefault("ALL", sorted(set(V["trait"]) - set(cfg["traitGroups"]) - {"ALL"}))
    cp = color_palette(methods)
    thr_by = {m["method"]: m["threshold"] for m in methods}
    cp["score_threshold"] = cp["method"].map(thr_by)
    write_tsv(cp.drop(columns="score_threshold"), pdir / "colorPalette.tsv", quiet=True)
    ranges = metric_ranges(thr_tables, link_tables, alpha)
    write_tsv(ranges, pdir / "metricRanges.tsv")
    mats = pdir / "heatmapMatrices"
    for m in thr_tables:
        if no_plots:
            plot_method_heatmaps_tables_only(m, thr_tables[m], link_tables.get(m), mats)
        else:
            plot_method_heatmaps(m, thr_tables[m], link_tables.get(m), ranges, cfg, pdir, figures, mats)
    comparisons = {}
    for name in comp["name"].unique():
        sep, agg = comparison_curves(name, comp, curves_by_method, cp, alpha, compat)
        vals = pd.DataFrame()
        if len(sep):
            s2 = curve_plot_filter(sep, cp, alpha)
            a2 = curve_plot_filter(agg, cp, alpha)
            vals = pd.concat([s2, a2], ignore_index=True)
            write_tsv(vals, pdir / "enrichmentRecallCurves" / f"{name}.values.tsv.gz")
            if not no_plots:
                plot_comparison_curves(name, s2, a2, cp, V, comp, trait_groups, pdir, figures)
        ov_in, ov_s, lk_in, lk_s = performance_comparison(name, comp, thr_tables, link_tables, baseline, cp, trait_groups,
                                                          cfg["biosampleGroups"], alpha, compat)
        if len(ov_s):
            write_tsv(ov_s, pdir / "thresholdedPerformanceComparison" / f"{name}.variantOverlap.tsv")
        if len(lk_s):
            write_tsv(lk_s, pdir / "thresholdedPerformanceComparison" / f"{name}.geneLinking.tsv")
        if not no_plots:
            plot_performance(name, ov_in, ov_s, lk_in, lk_s, V, trait_groups, pdir, figures, cp)
        comparisons[name] = {"overlap": ov_s, "linking": lk_s, "curve_values": vals,
                             "pairs": comp[comp["name"] == name][["biosample", "trait"]].values.tolist()}
    return {"figures": figures, "ranges": ranges, "comparisons": comparisons, "palette": cp}


def plot_method_heatmaps_tables_only(method: str, thr, link, mats_out: Path) -> None:
    np = _np()
    res = thr[np.isfinite(thr["enrichment"]) & np.isfinite(thr["recall"]) & (thr["nVariantsTotal"] > 20)]
    if len(res):
        write_tsv(heatmap_matrix(res, "enrichment").reset_index(), mats_out / f"{method}.enrichmentOverlapHeatmap.matrix.tsv", quiet=True)
        write_tsv(heatmap_matrix(res, "recall").reset_index(), mats_out / f"{method}.recallOverlapHeatmap.matrix.tsv", quiet=True)
    if link is not None and len(link):
        lk = link[(link["nCredibleSetsTotal"] > 0)]
        for metric in ("precision", "recall"):
            for ip in (False, True):
                d = lk[lk["intersectPoPS"].astype(bool) == ip]
                if len(d):
                    tag = ("intPoPS." if ip else "") + f"{metric}Linking"
                    write_tsv(heatmap_matrix(d, metric).reset_index(), mats_out / f"{method}.{tag}Heatmap.matrix.tsv", quiet=True)


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

def _row_dict(r: dict, keys: "list[str]") -> dict:
    out = {}
    for k in keys:
        v = r.get(k)
        if isinstance(v, (bool,)) or v is None:
            out[k] = v
        else:
            try:
                f = float(v)
                out[k] = f if math.isfinite(f) else None
            except (TypeError, ValueError):
                out[k] = str(v)
    return out


OV_KEYS = ["nVariantsOverlappingEnhancers", "nVariantsTotal", "nCommonVariantsOverlappingEnhancers", "nCommonVariantsTotal",
           "recall", "enrichment", "CI_enr_low", "CI_enr_high", "p_enr", "p_adjust_enr", "CI_recall_low", "CI_recall_high"]
LK_KEYS = ["nCredibleSetsTotal", "precision", "recall", "CI_precision_low", "CI_precision_high", "CI_recall_low", "CI_recall_high"]


def save_reference_configs(out: Path, cfg: dict, methods: "list[Method]", comp) -> None:
    rd = out / "reference_configs"
    rd.mkdir(parents=True, exist_ok=True)
    tg = dict(cfg["traitGroups"])
    lines = [f"METHODS:  {[m.name for m in methods]}", f"BIOSAMPLE_GROUPS:  {cfg['biosampleGroups']}", f"TRAIT_GROUPS:  {tg}",
             f"NUM_POPS_GENES:  {cfg['numPoPSGenes']}", f"NUM_PRED_GENES: {cfg['numPredGenes']}",
             f"THRESHOLD_PIP:  {cfg['thresholdPIP']}", f"THRESHOLD_PVAL:  {cfg['thresholdPval']}", f"N_THRESHOLD_STEPS:  {cfg['nThresholdSteps']}"]
    (rd / "config_params.tsv").write_text("\n".join(lines) + "\n")
    write_tsv(comp, rd / "comparisons.tsv", quiet=True)
    pd = _pd()
    for fname, groups in (("trait_groups.tsv", tg), ("biosample_groups.tsv", cfg["biosampleGroups"])):
        if groups:
            n = max(len(v) for v in groups.values())
            pd.DataFrame({g: list(v) + [None] * (n - len(v)) for g, v in groups.items()}).to_csv(rd / fname, sep="\t", index=False)
    (rd / "methods_config.json").write_text(json.dumps([m.as_dict() for m in methods], indent=2))
    (rd / "run_config.json").write_text(json.dumps({k: v for k, v in cfg.items()}, indent=2, default=str))


def cmd_run(args: argparse.Namespace) -> int:
    pd = _pd()
    t0 = time.time()
    cfg = load_config(args)
    compat = bool(args.upstream_compat)
    errors, warnings, res = validate(cfg, check_files=True)
    if args.variants_dir:
        errors = [e for e in errors if not e.startswith(("partition", "bgVariants", "chrSizes"))]
    if errors:
        for e in errors:
            print("  ERROR", e, file=sys.stderr)
        raise SystemExit(f"config invalid ({len(errors)} errors); see `igvfagent gwas-e2g validate-config`")
    mt, pt, comp, key, methods = res
    for w in warnings:
        print("  WARN ", w)
    out = run_dir(args.label)
    print(f"Run dir: {out}")
    # variants
    if args.variants_dir:
        vres = load_variants_dir(Path(args.variants_dir))
        if not cfg["traitGroups"] and vres["trait_groups"]:
            cfg["traitGroups"] = {g: m for g, m in vres["trait_groups"].items() if g != "ALL"}
    else:
        vres = process_variants(cfg, key)
        write_variants_dir(vres, out, cfg)
    V = vres["variants"].copy()
    bg = vres["bg"]
    if bg is None:
        raise SystemExit("background SNPs are required (--bg-variants or a --variants-dir that has them)")
    bg = bg.reset_index(drop=True).copy()
    bg["bg_id"] = range(len(bg))
    n_bg_total = len(bg)
    V["loc_id"] = V.groupby(["chr", "start", "end"], sort=False).ngroup()
    locs = V.drop_duplicates("loc_id")[["chr", "start", "end", "loc_id"]].reset_index(drop=True)
    loc_trait = V[["loc_id", "trait"]].drop_duplicates()
    nvt = n_variants_per_trait(V)
    chroms = set(read_chr_sizes(cfg["chrSizes"])) if cfg.get("chrSizes") and Path(cfg["chrSizes"]).is_file() else set(V["chr"]) | set(bg["chr"])
    universe = load_gene_universe(cfg["TSS"])
    gp = gp_raw = baseline = None
    if cfg.get("genePrioritizationTable") and Path(cfg["genePrioritizationTable"]).is_file():
        gp_raw = load_gene_prioritization(cfg["genePrioritizationTable"], None)
        gp = gp_raw[gp_raw["TargetGene"].isin(universe)].reset_index(drop=True)
        baseline = evaluate_baseline(gp, cfg["numPoPSGenes"], cfg["thresholdPval"])
        write_tsv(baseline, out / "baseline" / "gene_linking" / "precisionRecall.byTrait.tsv.gz")
    save_reference_configs(out, cfg, methods, comp)
    results = {}
    for m in methods:
        print(f"method {m.name}: {len(m.biosamples)} biosamples, groups {', '.join(m.groups)}")
        results[m.name] = per_method(m, cfg, comp, V, locs, loc_trait, bg, n_bg_total, nvt, universe, chroms, gp, gp_raw, out, compat)
    thr_tables = {k: v["thresholded"] for k, v in results.items()}
    link_tables = {k: v["linking"] for k, v in results.items()}
    curves = {k: v["curves"] for k, v in results.items()}
    vis = visualize(out, cfg, [m.as_dict() for m in methods], comp, V, thr_tables, link_tables, curves, baseline, compat, args.no_plots)
    summ = build_summary(args.label, cfg, methods, results, vis, nvt, n_bg_total, baseline, compat, time.time() - t0)
    (out / "summary.json").write_text(json.dumps(summ, indent=2, default=str))
    (out / "report.md").write_text(build_report(args.label, cfg, methods, results, vis, nvt, n_bg_total, vres, compat))
    print(f"JSON: {out / 'summary.json'}")
    print(f"Report: {out / 'report.md'}")
    return 0


def build_summary(label, cfg, methods, results, vis, nvt, n_bg_total, baseline, compat, secs) -> dict:
    summ: dict = {"label": label, "upstream": f"{UPSTREAM_REPO}@{UPSTREAM_COMMIT}", "upstream_compat": compat,
                  "params": {k: cfg[k] for k in ("nThresholdSteps", "thresholdPIP", "thresholdPval", "numPoPSGenes", "numPredGenes", "plotFixedScale")},
                  "trait_groups": cfg["traitGroups"], "biosample_groups": cfg["biosampleGroups"],
                  "nVariantsTotal": {r.trait: int(r.nVariantsTotal) for r in nvt.itertuples()}, "nCommonVariantsTotal": int(n_bg_total),
                  "methods": {}, "comparisons": {}, "figures": vis["figures"], "seconds": round(secs, 1)}
    for m in methods:
        r = results[m.name]
        thr = r["thresholded"]
        st = r["stats"].set_index("biosample")
        allrow = thr[(thr["biosample"] == "ALL") & (thr["trait"] == "ALL")]
        summ["methods"][m.name] = {"threshold": m.threshold, "inverse_predictor": m.inverse, "boolean": m.boolean,
                                   "biosamples": m.biosamples, "biosample_groups": m.groups, "n_thresholds": len(r["span"]),
                                   "bp_enhancers": {b: int(v) for b, v in st["bp_enhancers"].items()},
                                   "bg_variant_count": {b: int(v) for b, v in st["bg_variant_count"].items()},
                                   "ALL_x_ALL": _row_dict(allrow.iloc[0].to_dict(), OV_KEYS) if len(allrow) else None,
                                   "n_trait_biosample_rows": int(len(thr)), "n_linking_rows": int(len(r["linking"]))}
    for name, c in vis["comparisons"].items():
        d: dict = {"pairs": c["pairs"], "overlap": {}, "linking": {}}
        for rr in c["overlap"].to_dict("records") if len(c["overlap"]) else []:
            d["overlap"][rr["method"]] = _row_dict(rr, OV_KEYS)
        for rr in c["linking"].to_dict("records") if len(c["linking"]) else []:
            ip = rr.get("intersectPoPS")
            tag = rr["method"] if ip is None or (isinstance(ip, float) and math.isnan(ip)) else f"{rr['method']}{'.intPoPS' if bool(ip) else ''}"
            d["linking"][tag] = _row_dict(rr, LK_KEYS)
        summ["comparisons"][name] = d
    if baseline is not None:
        summ["baseline_rows"] = int(len(baseline))
    return summ


def build_report(label, cfg, methods, results, vis, nvt, n_bg_total, vres, compat) -> str:
    L = [f"# GWAS E2G benchmark — {label}", "",
         f"Port of [{UPSTREAM_REPO}](https://github.com/{UPSTREAM_REPO}) (commit `{UPSTREAM_COMMIT[:10]}`). "
         f"Fine-mapped variants with PIP > {cfg['thresholdPIP']} in distal noncoding sequence vs {n_bg_total:,} common background SNPs; "
         f"gene linking against the silver-standard credible-set genes (top {cfg['numPredGenes']} predicted genes; PoPS top {cfg['numPoPSGenes']}). "
         f"{'**--upstream-compat**: upstream quirks reproduced.' if compat else 'Upstream quirks corrected (see the module docstring); rerun with --upstream-compat to reproduce them.'}",
         "", "## Variants", "",
         md_table(["trait", "nVariantsTotal"], [(r.trait, r.nVariantsTotal) for r in nvt.itertuples()], 40), "",
         "Trait groups: " + ("; ".join(f"{g} = {', '.join(m)}" for g, m in cfg["traitGroups"].items()) or "none") + " (+ ALL).",
         "", "## Methods", "",
         md_table(["method", "threshold", "biosamples", "groups", "thresholds in span", "ALL x ALL recall", "ALL x ALL enrichment"],
                  [(m.name, _fmt(m.threshold), len(m.biosamples), ", ".join(m.groups), len(results[m.name]["span"]),
                    _fmt((results[m.name]["thresholded"].query("biosample == 'ALL' and trait == 'ALL'")["recall"].tolist() or [float("nan")])[0]),
                    _fmt((results[m.name]["thresholded"].query("biosample == 'ALL' and trait == 'ALL'")["enrichment"].tolist() or [float("nan")])[0]))
                   for m in methods]), ""]
    for name, c in vis["comparisons"].items():
        L += [f"## Comparison `{name}`", "", "Pairs: " + ", ".join(f"{b} x {t}" for b, t in c["pairs"]), ""]
        if len(c["overlap"]):
            L += ["Variant overlap (counts summed over the pairs):", "",
                  md_table(["method", "recall", "enrichment", "95% CI", "p_adj", "variants in enh / total"],
                           [(r["method"], _fmt(r["recall"]), _fmt(r["enrichment"]), f"{_fmt(r['CI_enr_low'])}–{_fmt(r['CI_enr_high'])}",
                             _fmt(r["p_adjust_enr"]), f"{int(r['nVariantsOverlappingEnhancers'])}/{int(r['nVariantsTotal'])}")
                            for r in c["overlap"].sort_values("enrichment", ascending=False).to_dict("records")]), ""]
        if len(c["linking"]):
            def ipl(v):
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    return "baseline"
                return "int. PoPS" if bool(v) else "E2G only"
            L += ["Gene linking:", "",
                  md_table(["method", "mode", "precision", "recall", "credible sets"],
                           [(r["method"], ipl(r.get("intersectPoPS")), _fmt(r["precision"]), _fmt(r["recall"]), int(r["nCredibleSetsTotal"]))
                            for r in c["linking"].to_dict("records")]), ""]
    if vis["figures"]:
        L += ["## Figures", ""] + [f"- `plots/{Path(f).parent.name}/{Path(f).name}`" for f in vis["figures"]] + [""]
    L += ["## Outputs", "", "Per method: `<method>/variant_overlap/enrichmentRecall.thresholded.traitByBiosample.tsv.gz`, "
          "`bgOverlap.enhancerSetSize.thresholdedPredictions.tsv`, `enrichmentRecallAcrossThresholds.forBiosample<b>.tsv.gz`, "
          "`thresholdSpan.tsv`, `gene_linking/precisionRecall.traitByBiosample.tsv.gz`; `baseline/gene_linking/precisionRecall.byTrait.tsv.gz`; "
          "`plots/` (colour palette, metric ranges, heatmap matrices, comparison tables); `reference_configs/`.", ""]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# plot (visualisation stage from a run directory)
# ---------------------------------------------------------------------------

def cmd_plot(args: argparse.Namespace) -> int:
    pd = _pd()
    rd = Path(args.run_dir).resolve()
    cfg = json.loads((rd / "reference_configs" / "run_config.json").read_text())
    if args.plot_fixed_scale:
        cfg["plotFixedScale"] = True
    methods = json.loads((rd / "reference_configs" / "methods_config.json").read_text())
    comp = pd.read_csv(rd / "reference_configs" / "comparisons.tsv", sep="\t", dtype=str)
    vdir = rd / "variants" if (rd / "variants").is_dir() else (Path(args.variants_dir) if args.variants_dir else None)
    if vdir is None:
        raise SystemExit("run dir has no variants/ (it was run with --variants-dir); pass --variants-dir")
    V = load_variants_dir(vdir)["variants"]
    thr_tables, link_tables, curves = {}, {}, {}
    for m in methods:
        name = m["method"]
        p = rd / name / "variant_overlap" / "enrichmentRecall.thresholded.traitByBiosample.tsv.gz"
        if not p.is_file():
            continue
        thr_tables[name] = pd.read_csv(p, sep="\t")
        lp = rd / name / "gene_linking" / "precisionRecall.traitByBiosample.tsv.gz"
        link_tables[name] = pd.read_csv(lp, sep="\t") if lp.is_file() else None
        curves[name] = {}
        for f in sorted((rd / name / "variant_overlap").glob("enrichmentRecallAcrossThresholds.forBiosample*.tsv.gz")):
            b = f.name[len("enrichmentRecallAcrossThresholds.forBiosample"):-len(".tsv.gz")]
            curves[name][b] = pd.read_csv(f, sep="\t")
    bp = rd / "baseline" / "gene_linking" / "precisionRecall.byTrait.tsv.gz"
    baseline = pd.read_csv(bp, sep="\t") if bp.is_file() else None
    compat = bool(args.upstream_compat)
    vis = visualize(rd, cfg, methods, comp, V, thr_tables, link_tables, curves, baseline, compat, args.no_plots)
    summ = {"run_dir": str(rd), "figures": vis["figures"], "comparisons": list(vis["comparisons"]), "plotFixedScale": cfg["plotFixedScale"]}
    (rd / "plots" / "summary.json").write_text(json.dumps(summ, indent=2))
    rep = [f"# GWAS E2G benchmark — re-plotted", "", f"Run: `{rd}`", "", f"{len(vis['figures'])} figures; comparisons: "
           + ", ".join(vis["comparisons"]), ""] + [f"- `{f}`" for f in vis["figures"]]
    (rd / "plots" / "report.md").write_text("\n".join(rep) + "\n")
    print(f"JSON: {rd / 'plots' / 'summary.json'}")
    print(f"Report: {rd / 'plots' / 'report.md'}")
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def synthetic_world(d: Path, seed: int = 11) -> dict:
    """Two 2 Mb chromosomes, 60 genes, three traits x 15 credible sets with planted enhancers, a
    silver-standard gene table with PoPS scores and distances, four predictors x three biosamples."""
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    size = 2_000_000
    chroms = ["chr1", "chr2"]
    (d / "chrom.sizes.tsv").write_text("".join(f"{c}\t{size}\n" for c in chroms))
    genes = {}
    for c in chroms:
        for i in range(30):
            genes[(c, i)] = (f"G{c[-1]}{i:02d}", 5000 + 40000 * i)
    tss = pd.DataFrame([(c, t - 250, t + 250, n, 0, "+", f"ENSG{k:05d}", "protein_coding") for k, ((c, i), (n, t)) in enumerate(genes.items())],
                       columns=["chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"])
    with open(d / "tss.bed", "w") as fh:
        fh.write("#chr\tstart\tend\tname\tscore\tstrand\tEnsembl_ID\tgene_type\n")
        tss.to_csv(fh, sep="\t", header=False, index=False)
    special = [(c, t - 250, t + 250, "TSS-500bp") for (c, _), (_, t) in genes.items()]
    special += [("chr1", 1_500_000, 1_700_000, "CDS"), ("chr2", 1_300_000, 1_310_000, "UTR")]
    rows = []
    cyc = ["Other", "OtherIntron", "ABC"]
    for c in chroms:
        sp = sorted([s for s in special if s[0] == c], key=lambda x: x[1])
        pos, k = 0, 0
        for _, s, e, cat in sp + [(c, size, size, None)]:
            while pos < s:
                nxt = min(s, pos + 50_000)
                rows.append((c, pos, nxt, cyc[k % 3])); k += 1; pos = nxt
            if cat:
                rows.append((c, s, e, cat)); pos = e
    pd.DataFrame(rows).to_csv(d / "partition.bed", sep="\t", header=False, index=False)

    def slot(c, j):
        return 20_000 + 40_000 * j
    traits = {"T1": ("chr1", 0), "T2": ("chr1", 15), "T3": ("chr2", 0)}
    var_rows = {t: [] for t in traits}
    gp_rows = []
    cs_info = {}
    for ti, (t, (c, off)) in enumerate(traits.items()):
        for i in range(15):
            j = off + i
            s = slot(c, j)
            cs = f"{c}:{s}-{s + 500}-{ti + 1}"
            truth, decoy = genes[(c, j)][0], genes[(c, (j + 1) % 30)][0]
            others = [genes[(c, (j + k) % 30)][0] for k in (2, 3, 4)]
            cs_info[(t, i)] = (c, s, cs, truth, decoy)
            rs = f"rs{t}_{i}"
            var_rows[t] += [(c, s + 100, s + 101, rs + "_1", 0.4, cs, t), (c, s + 300, s + 301, rs + "_2", 0.3, cs, t),
                            ("chr2", 1_500_000 + (ti * 15 + i) * 2000 + 777, 1_500_000 + (ti * 15 + i) * 2000 + 778, rs + "_3", 0.2, cs, t),
                            (c, s + 200, s + 201, rs + "_4", 0.1, cs, t),
                            ("chr1", 1_500_000 + (ti * 15 + i) * 3000 + 55, 1_500_000 + (ti * 15 + i) * 3000 + 56, rs + "_5", 0.5, cs, t)]
            pops = {others[0]: 5.0, truth: 4.0, decoy: 1.0, others[1]: 0.5, others[2]: 0.2}
            dist = {truth: 1000 if i < 8 else 30000, decoy: 20000, others[0]: 25000, others[1]: 40000, others[2]: 50000}
            for g in [truth, decoy] + others:
                gp_rows.append({"TruthDistanceRank": 1, "CredibleSet": cs, "Disease": t, "TargetGene": g, "DistanceRank": 1,
                                "PromoterDistanceToBestSNP": dist[g], "POPS.Score": pops[g], "truth": "TRUE" if g == truth else "FALSE", "junk": "x"})
            if t == "T1" and i < 3:
                gp_rows.append({"TruthDistanceRank": 1, "CredibleSet": cs, "Disease": t, "TargetGene": "NOTAGENE", "DistanceRank": 1,
                                "PromoterDistanceToBestSNP": 999999, "POPS.Score": 10.0, "truth": "FALSE", "junk": "x"})
    # T2 credible sets 0-4 carry one extra variant at T1's CS i v1 location (shared rsid)
    for i in range(5):
        c, s1, _, _, _ = cs_info[("T1", i)]
        _, _, cs2, _, _ = cs_info[("T2", i)]
        var_rows["T2"].append((c, s1 + 100, s1 + 101, f"rsT1_{i}_1", 0.25, cs2, "T2"))
    key_rows = []
    for t, rws in var_rows.items():
        p = d / "variants" / t / "variant.list.txt"
        p.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rws, columns=["chr", "start", "end", "rsid", "pip", "CredibleSet", "trait"])
        df.insert(4, "beta", 0.01)
        df.to_csv(p, sep="\t", index=False)
        key_rows.append((t, f"variants/{t}/variant.list.txt"))
    pd.DataFrame(key_rows, columns=["trait", "variant_file"]).to_csv(d / "variant_key.tsv", sep="\t", index=False)
    pd.DataFrame(gp_rows).to_csv(d / "gene_prioritization.tsv", sep="\t", index=False)
    # background SNPs
    bgc = rng.choice(chroms, 40_000)
    bgp = rng.integers(0, size - 1, 40_000)
    bg = pd.DataFrame({"chr": bgc, "start": bgp, "end": bgp + 1, "rsid": [f"bg{k}" for k in range(40_000)]})
    bg.to_csv(d / "bg.bed", sep="\t", header=False, index=False)

    # predictions: element-gene rows with a "good" score; other methods re-score the same rows
    def cs_rows(t, strong_upto):
        out = []
        for i in range(15):
            c, s, _, truth, decoy = cs_info[(t, i)]
            hi = i < strong_upto
            out += [(c, s, s + 500, truth, 0.8 if hi else 0.4), (c, s, s + 500, decoy, 0.6 if hi else 0.3)]
        return out

    def background(n):
        c = rng.choice(chroms, n)
        s = rng.integers(1_200_000, 1_450_000, n)
        g = rng.choice([v[0] for v in genes.values()], n)
        return [(c[k], int(s[k]), int(s[k]) + 500, g[k], float(rng.uniform(0, 0.45))) for k in range(n)]
    base = {"B1": cs_rows("T1", 10) + cs_rows("T2", 15) + background(200),
            "B2": cs_rows("T3", 15) + background(200),
            "B3": background(200)}
    tss_of = {n: t for (_, _), (n, t) in genes.items()}
    methods_cfg = {"good": ("score", False, False, 0.5, "#792374", "Good predictor"),
                   "random": ("score", False, False, 0.5, "", "Random scores"),
                   "distance": ("distance", False, True, 30000, "#96a0b3", "Distance to TSS"),
                   "binary": ("HC_binary", True, False, 1, "#e0dcca", "Binary high-confidence")}
    pred_table = {"biosample": ["B1", "B2", "B3"]}
    for m, (col, boolean, inv, thr, color, long) in methods_cfg.items():
        pred_table[m] = []
        for b, rws in base.items():
            df = pd.DataFrame(rws, columns=["chr", "start", "end", "TargetGene", "good"])
            if m == "good":
                df[col] = df["good"]
            elif m == "random":
                df[col] = rng.uniform(0, 1, len(df))
            elif m == "distance":
                df[col] = ((df["start"] + df["end"]) / 2 - df["TargetGene"].map(tss_of)).abs()
            else:
                df[col] = (df["good"] >= 0.5).astype(int)
            df = df.drop(columns="good")
            # an out-of-universe gene with a top score on the reserved region (must be dropped by the gene universe)
            df.loc[len(df)] = ["chr2", 1_499_000, 1_600_000, "NOTAGENE", 1 if m == "binary" else (0.0 if m == "distance" else 0.99)]
            df["CellType"] = b
            if m == "binary" and b == "B3":
                pred_table[m].append("")
                continue
            p = d / "pred" / m / f"{b}.tsv.gz"
            p.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(p, "wt") as fh:
                fh.write("#" + "\t".join(df.columns) + "\n")
                df.to_csv(fh, sep="\t", header=False, index=False)
            pred_table[m].append(str(p))
    pd.DataFrame(pred_table).to_csv(d / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame([(m, str(v[1]).upper(), str(v[2]).upper(), v[5], v[3], v[0], v[4]) for m, v in methods_cfg.items()],
                 columns=["method", "boolean", "inverse_predictor", "pred_name_long", "threshold", "score_col", "color"]
                 ).to_csv(d / "methods.tsv", sep="\t", index=False)
    pd.DataFrame([("matched", "B1", "T12"), ("matched", "B2", "T3"), ("pooled", "B12", "T1"), ("all", "ALL", "ALL")],
                 columns=["name", "biosample", "trait"]).to_csv(d / "comparisons.tsv", sep="\t", index=False)
    (d / "config.yml").write_text(
        "## CONFIG\nresults: \"results/\"\nmethodsTable: \"methods.tsv\"\npredictionsTable: \"predictions.tsv\"\n"
        "methods: [good, random, distance, binary]\ncomparisonsTable: \"comparisons.tsv\" # comparisons\n\n"
        "biosampleGroups: # groups\n   B12: [B1, B2]\n\ntraitGroups:\n  T12: [T1, T2]\n\n"
        "nThresholdSteps: 25 # steps\nthresholdPIP: 0.1\nthresholdPval: 0.05\nnumPoPSGenes: 2\nnumPredGenes: 2\nplotFixedScale: False\n"
        "variantKey: \"variant_key.tsv\"\ngenePrioritizationTable: \"gene_prioritization.tsv\"\nbgVariants: \"bg.bed\"\n"
        "chrSizes: \"chrom.sizes.tsv\"\npartition: \"partition.bed\"\nTSS: \"tss.bed\"\n")
    return {"config": d / "config.yml", "bg": bg, "cs_info": cs_info, "base": base}


def _args(**kw) -> argparse.Namespace:
    base = dict(config=None, methods_table=None, predictions_table=None, comparisons_table=None, variant_key=None,
                gene_prioritization_table=None, bg_variants=None, partition=None, tss=None, chr_sizes=None,
                n_threshold_steps=None, threshold_pip=None, threshold_pval=None, num_pops_genes=None, num_pred_genes=None,
                methods=None, plot_fixed_scale=False, biosample_group=None, trait_group=None, label="selftest",
                upstream_compat=False, no_plots=True, variants_dir=None, no_file_checks=False)
    base.update(kw)
    return argparse.Namespace(**base)


def cmd_selftest(args: argparse.Namespace) -> int:
    import shutil
    import tempfile
    pd, np = _pd(), _np()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made: "list[Path]" = []
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_world(d)
        cfgp = str(W["config"])

        print("\nhelpers")
        mi = merge_intervals(pd.DataFrame({"chr": ["c"] * 3, "start": [0, 50, 200], "end": [100, 120, 300]}))
        q = pd.DataFrame({"chr": ["c", "c", "c", "d"], "start": [119, 120, 150, 5], "end": [120, 121, 250, 6]})
        check(list(mi["c"][0]) == [0, 200] and list(mi["c"][1]) == [120, 300] and list(overlaps_union(q, mi)) == [True, False, True, False],
              "merge_intervals / overlaps_union: half-open union overlap")
        e, lo, hi, se, p = enrichment_stats([50], [100], [100], [10_000], 0.05)
        _, _, _, sec, _ = enrichment_stats([50], [100], [100], [10_000], 0.05, compat=True)
        check(abs(se[0] - math.sqrt(50 / 5000 + 9900 / 1e6)) < 1e-12 and abs(sec[0] - (math.sqrt(50 / 5000) + 9900 / 1e6)) < 1e-12,
              "helpful_math: SE_log_enr = sqrt(a + b); --upstream-compat reproduces sqrt(a) + b")
        t = pd.DataFrame({"nVariantsOverlappingEnhancers": [10], "nVariantsTotal": [40]}); t["recall"] = 0.25
        t = add_recall_overlap_statistics(t, 0.05)
        ra = 12 / 44
        check(abs(t["SE_recall"][0] - math.sqrt(ra * (1 - ra) / 44)) < 1e-12 and abs(t["CI_recall_low"][0] - (0.25 - 1.959964 * t["SE_recall"][0])) < 1e-5,
              "helpful_math: Agresti-Coull recall SE ((x+2)/(n+4)), CI centred on the raw recall")
        span = quantile_threshold_span([0.1, 0.2, 0.2, 0.3, 0.9], 5, False, 0.55)
        check(0.55 in span and span == sorted(span) and 0.1 in span and 0.9 in span and len(span) >= 5,
              f"threshold span: quantiles + provided threshold + even grid ({len(span)} values)")
        check(quantile_threshold_span([0, 0, 1, 1, 1], 25, True, 1) == [0.0, 1.0], "threshold span: boolean predictor -> n = 2 -> {0, 1}")
        check(quantile_threshold_span([0.1, 0.1, 0.1, 0.9], 3, False, 0.5, compat=True) == sorted({0.1, 0.5, 0.9}),
              "threshold span --upstream-compat: quantiles over distinct score values")
        rows = pd.DataFrame({"CredibleSet": ["cs1"] * 4, "TargetGene": ["A", "B", "C", "A"], "score": [0.9, 0.9, 0.5, 0.1]})
        (a0, b0, c0), _ = gene_linking_counts(rows, pd.DataFrame({"CredibleSet": ["cs1"], "trait": ["T"], "TruthGene": ["C"]}),
                                               pd.DataFrame({"CredibleSet": ["cs1"], "trait": ["T"], "PoPSGene": ["A"]}), 2)
        check((a0, b0, c0) == (1, 2, 0), "gene linking: rank(-max, ties = 'min') keeps both tied top genes and drops rank 3")
        yml = parse_config_text(Path(cfgp).read_text())
        check(yml["methods"] == ["good", "random", "distance", "binary"] and yml["biosampleGroups"] == {"B12": ["B1", "B2"]}
              and yml["plotFixedScale"] is False and yml["thresholdPIP"] == 0.1 and yml["comparisonsTable"] == "comparisons.tsv",
              "config.yml parser: inline lists, nested groups, booleans, numbers, trailing comments")

        print("\nvalidate-config")
        check(cmd_validate(_args(config=cfgp)) == 0, "validate-config: synthetic config valid")
        bad = d / "comparisons_bad.tsv"
        bad.write_text("name\tbiosample\ttrait\nx\tB9\tT1\ny\tB1\tNOPE\n")
        errs, _, _ = validate(load_config(_args(config=cfgp, comparisons_table=str(bad), biosample_group=["B1=B2"])))
        check(any("B9" in e for e in errs) and any("NOPE" in e for e in errs) and any("cannot share the same name" in e for e in errs),
              "validate-config: undefined comparison biosample / trait and group-name clash are errors")

        print("\nvariants")
        rc = cmd_variants(_args(config=cfgp, label="st_gwas_variants"))
        vrun = sorted(OUT_ROOT.glob("*_st_gwas_variants"))[-1]; made.append(vrun)
        vs = json.loads((vrun / "summary.json").read_text())
        nv = vs["nVariantsTotal"]
        check(rc == 0 and nv == {"T1": 45, "T2": 50, "T3": 45, "T12": 90, "ALL": 135},
              f"variants: PIP > 0.1 strict, CDS dropped, trait group T12 and ALL de-duplicated by location {nv}")
        bgw = W["bg"]
        n_in_cds = int(((bgw["chr"] == "chr1") & (bgw["start"] >= 1_500_000) & (bgw["start"] < 1_700_000)).sum())
        check(vs["background_snps_distal_noncoding"] < 40_000 - n_in_cds + 1 and vs["background_snps_distal_noncoding"] > 0.9 * 40_000,
              f"variants: background SNPs in CDS / UTR / TSS removed ({vs['background_snps_distal_noncoding']} of 40000)")

        print("\nbaseline")
        rc = cmd_baseline(_args(config=cfgp, label="st_gwas_baseline"))
        brun = sorted(OUT_ROOT.glob("*_st_gwas_baseline"))[-1]; made.append(brun)
        bl = pd.read_csv(brun / "baseline" / "gene_linking" / "precisionRecall.byTrait.tsv.gz", sep="\t")
        pops = bl[bl["method"] == "PoPS"].set_index("trait"); dist = bl[bl["method"] == "distanceToTSS"].set_index("trait")
        check(rc == 0 and (pops["precision"] == 0.5).all() and (pops["recall"] == 1.0).all(),
              "baseline PoPS: top 2 by POPS.Score hold the truth gene (out-of-universe gene removed) -> precision 0.5, recall 1")
        check(np.allclose(dist["recall"], 8 / 15) and np.allclose(dist["precision"], 8 / 30),
              "baseline distanceToTSS: truth nearest in 8 of 15 credible sets -> recall 8/15")

        print("\nrun")
        rc = cmd_run(_args(config=cfgp, label="st_gwas_run", no_plots=args.no_plots))
        run = sorted(OUT_ROOT.glob("*_st_gwas_run"))[-1]; made.append(run)
        summ = json.loads((run / "summary.json").read_text())
        check(rc == 0 and (run / "report.md").is_file() and set(summ["methods"]) == {"good", "random", "distance", "binary"},
              "run: report + summary for 4 methods")
        thr = pd.read_csv(run / "good" / "variant_overlap" / "enrichmentRecall.thresholded.traitByBiosample.tsv.gz", sep="\t")
        T = thr.set_index(["biosample", "trait"])
        expect = {("B1", "T1"): 20, ("B1", "T2"): 35, ("B1", "T12"): 50, ("B1", "ALL"): 50, ("B2", "T3"): 30,
                  ("B12", "T1"): 20, ("B12", "T3"): 30, ("ALL", "ALL"): 80}
        got = {k: int(T.loc[k, "nVariantsOverlappingEnhancers"]) if k in T.index else None for k in expect}
        check(got == expect, f"run: nVariantsOverlappingEnhancers at threshold (strong enhancers only, groups = union) {got}")
        check(("B1", "T3") not in T.index and ("B3", "T1") not in T.index, "run: trait x biosample rows only where variants overlap (swapped pairs absent)")
        check(abs(T.loc[("B1", "T1"), "recall"] - 20 / 45) < 1e-12 and abs(T.loc[("B1", "T2"), "recall"] - 0.7) < 1e-12,
              "run: recall = overlapping / nVariantsTotal (the reserved-region variant and the NOTAGENE element are ignored)")
        st = pd.read_csv(run / "good" / "variant_overlap" / "bgOverlap.enhancerSetSize.thresholdedPredictions.tsv", sep="\t").set_index("biosample")
        els = pd.DataFrame([r for r in W["base"]["B1"] if r[4] >= 0.5], columns=["chr", "start", "end", "g", "s"]).drop_duplicates(["chr", "start", "end"])
        hit = np.zeros(len(bgw), dtype=bool)
        for r in els.itertuples():
            hit |= (bgw["chr"].to_numpy() == r.chr) & (bgw["start"].to_numpy() >= r.start) & (bgw["start"].to_numpy() < r.end)
        check(st.loc["B1", "bp_enhancers"] == 25 * 500 and st.loc["B1", "bg_variant_count"] == int(hit.sum())
              and st.loc["B12", "bp_enhancers"] == 40 * 500 and st.loc["B3", "bp_enhancers"] == 0,
              f"run: enhancer set size (merged bp) and background overlap match brute force ({int(hit.sum())} SNPs in B1)")
        n_bg = summ["nCommonVariantsTotal"]
        r = T.loc[("B1", "T1")]
        check(abs(r["enrichment"] - (20 / 45) / (r["nCommonVariantsOverlappingEnhancers"] / n_bg)) < 1e-9 and r["enrichment"] > 50 and r["p_adjust_enr"] < 1e-10,
              f"run: enrichment = recall / background fraction ({r['enrichment']:.0f}x, p_adj {r['p_adjust_enr']:.1e})")
        check(abs(r["p_adjust_enr"] - min(1.0, r["p_enr"] * len(thr))) < 1e-300 + 1e-12 * r["p_adjust_enr"], "run: Bonferroni over the method's table")
        rnd = pd.read_csv(run / "random" / "variant_overlap" / "enrichmentRecall.thresholded.traitByBiosample.tsv.gz", sep="\t").set_index(["biosample", "trait"])
        check(rnd.loc[("B1", "T1"), "enrichment"] < r["enrichment"], f"run: random scores enrich less ({rnd.loc[('B1', 'T1'), 'enrichment']:.0f}x)")
        check(summ["methods"]["distance"]["threshold"] == -30000 and summ["methods"]["binary"]["n_thresholds"] == 2
              and summ["methods"]["good"]["n_thresholds"] >= 25, "run: inverse threshold negated; binary span {0, 1}; continuous span >= nThresholdSteps")
        check(summ["methods"]["binary"]["biosample_groups"] == ["ALL", "B12"] and "B3" not in summ["methods"]["binary"]["biosamples"],
              "run: blank prediction entries skipped; groups used when all members are present")
        span = pd.read_csv(run / "good" / "thresholdSpan.tsv", sep="\t")["threshold"].tolist()
        check(0.5 in span and span == sorted(span), "run: threshold span includes the method threshold")
        cur = pd.read_csv(run / "good" / "variant_overlap" / "enrichmentRecallAcrossThresholds.forBiosampleB1.tsv.gz", sep="\t")
        at = cur[np.isclose(cur["threshold"], 0.5) & (cur["trait"] == "T1")].iloc[0]
        mono = all(g.sort_values("threshold")["nVariantsOverlappingEnhancers"].is_monotonic_decreasing for _, g in cur.groupby("trait"))
        check(at["nVariantsOverlappingEnhancers"] == 20 and at["nCommonVariantsOverlappingEnhancers"] == T.loc[("B1", "T1"), "nCommonVariantsOverlappingEnhancers"] and mono,
              "run: curve at the method threshold equals the thresholded table; counts non-increasing in the threshold")
        check(all((run / "good" / "variant_overlap" / f"enrichmentRecallAcrossThresholds.forBiosample{b}.tsv.gz").is_file() for b in ("B1", "B2", "B12", "ALL")),
              "run: one curve table per comparison biosample (incl. groups and ALL)")
        lk = pd.read_csv(run / "good" / "gene_linking" / "precisionRecall.traitByBiosample.tsv.gz", sep="\t")
        L = lk.set_index(["biosample", "trait", "intersectPoPS"])
        l1, l1p, l2 = L.loc[("B1", "T1", False)], L.loc[("B1", "T1", True)], L.loc[("B1", "T2", False)]
        check(l1["nCredibleSetsOverlappingEnhancersAnyGene"] == 20 and l1["nCredibleSetsOverlappingEnhancersCorrectGene"] == 10
              and abs(l1["recall"] - 10 / 15) < 1e-12 and abs(l1["precision"] - 0.5) < 1e-12,
              "gene linking: top-2 genes of 10 strong credible sets -> precision 0.5, recall 10/15 (weak ones below threshold)")
        check(l1p["precision"] == 1.0 and abs(l1p["recall"] - 10 / 15) < 1e-12, "gene linking int. PoPS: decoys removed -> precision 1")
        check(l2["nCredibleSetsOverlappingEnhancersAnyGene"] == 30 and l2["nCredibleSetsOverlappingEnhancersCorrectGene"] == 15,
              "gene linking: tied top scores from a shared enhancer both kept, rank-3 genes dropped (T2)")
        check(("B12", "T1", False) in L.index and bool(L.loc[("B12", "T1", False), "group"]), "gene linking: biosample-group rows from concatenated predictions")
        pc = pd.read_csv(run / "plots" / "thresholdedPerformanceComparison" / "matched.variantOverlap.tsv", sep="\t").set_index("method")
        g = pc.loc["good"]
        check(g["nVariantsOverlappingEnhancers"] == 80 and g["nVariantsTotal"] == 135 and abs(g["recall"] - 80 / 135) < 1e-12,
              "comparison 'matched': B1 x T12 + B2 x T3 counts summed (recall 80/135)")
        pl = pd.read_csv(run / "plots" / "thresholdedPerformanceComparison" / "matched.geneLinking.tsv", sep="\t")
        gl = pl[(pl["method"] == "good") & (pl["intersectPoPS"].astype(str).str.lower() == "false")].iloc[0]
        check(abs(gl["recall"] - 40 / 45) < 1e-12 and abs(gl["precision"] - 0.5) < 1e-12 and set(pl["method"]) >= {"PoPS", "distanceToTSS"},
              "comparison gene linking: trait group expanded to T1 + T2, baselines appended")
        vals = pd.read_csv(run / "plots" / "enrichmentRecallCurves" / "matched.values.tsv.gz", sep="\t")
        check(set(vals["key"]) >= {"B1.T12", "B2.T3", "all_matched.all_matched"} and (vals["recall"] > 0.005).all() and (vals["p_adjust_enr"] < 0.05).all(),
              "curves values: per-pair and all_matched rows, plotting filters applied")
        agg = vals[(vals["key"] == "all_matched.all_matched") & (vals["method"] == "good") & np.isclose(vals["threshold"], 0.5)]
        check(len(agg) == 1 and agg.iloc[0]["nVariantsOverlappingEnhancers"] == 80, "curves: all_matched sums counts over the comparison's pairs")
        rng_ = pd.read_csv(run / "plots" / "metricRanges.tsv", sep="\t")
        check(set(rng_["metric"]) == {"enrichment_overlap", "recall_overlap", "precision_linking", "recall_linking",
                                      "precision_linking_intersectPoPS", "recall_linking_intersectPoPS"}, "metric ranges: 6 metrics x min / max per method")
        hm = pd.read_csv(run / "plots" / "heatmapMatrices" / "good.enrichmentOverlapHeatmap.matrix.tsv", sep="\t")
        check(set(hm["biosample"]) == {"B1", "B2", "B12", "ALL"} and "T12" in hm.columns, "heatmap matrix: biosample x trait, clustered order")
        cp = pd.read_csv(run / "plots" / "colorPalette.tsv", sep="\t").set_index("method")
        check(cp.loc["good", "hex"] == "#792374" and cp.loc["random", "hex"].startswith("#") and "PoPS" in cp.index, "colour palette: configured, filled, baselines")
        check((run / "reference_configs" / "methods_config.json").is_file() and (run / "reference_configs" / "trait_groups.tsv").is_file(),
              "reference configs saved")
        if not args.no_plots:
            check(len(summ["figures"]) >= 8, f"run: {len(summ['figures'])} figures")

        print("\nplot")
        rc = cmd_plot(argparse.Namespace(run_dir=str(run), variants_dir=None, upstream_compat=False, no_plots=True, plot_fixed_scale=True))
        check(rc == 0 and (run / "plots" / "report.md").is_file(), "plot: visualisation stage re-run from the run directory (fixed scale)")

        print("\nrun --upstream-compat (reusing --variants-dir)")
        rc = cmd_run(_args(config=cfgp, label="st_gwas_compat", methods=["good"], upstream_compat=True, variants_dir=str(vrun)))
        run2 = sorted(OUT_ROOT.glob("*_st_gwas_compat"))[-1]; made.append(run2)
        lk2 = pd.read_csv(run2 / "good" / "gene_linking" / "precisionRecall.traitByBiosample.tsv.gz", sep="\t").set_index(["biosample", "trait", "intersectPoPS"])
        c1, c1p = lk2.loc[("B1", "T1", False)], lk2.loc[("B1", "T1", True)]
        check(rc == 0 and c1["recall"] == 1.0 and c1["nCredibleSetsOverlappingEnhancersAnyGene"] == 30,
              "compat: individual biosamples linked with unthresholded predictions (weak credible sets counted)")
        check(abs(c1p["recall"] - 12 / 15) < 1e-12 and lk2.loc[("B12", "T1", False), "recall"] == l1["recall"],
              "compat: gene table not restricted to the gene universe (NOTAGENE displaces the truth from PoPS top 2); groups still thresholded")
        t2 = pd.read_csv(run2 / "good" / "variant_overlap" / "enrichmentRecall.thresholded.traitByBiosample.tsv.gz", sep="\t").set_index(["biosample", "trait"])
        check(t2.loc[("B1", "T1"), "SE_log_enr"] != T.loc[("B1", "T1"), "SE_log_enr"] and t2.loc[("B1", "T1"), "enrichment"] == T.loc[("B1", "T1"), "enrichment"],
              "compat: SE formula quirk reproduced; point estimates unchanged")
        pl2 = pd.read_csv(run2 / "plots" / "thresholdedPerformanceComparison" / "all.geneLinking.tsv", sep="\t")
        pl1 = pd.read_csv(run / "plots" / "thresholdedPerformanceComparison" / "all.geneLinking.tsv", sep="\t")
        n2 = pl2[(pl2["method"] == "good")]["nCredibleSetsTotal"].max(); n1 = pl1[(pl1["method"] == "good")]["nCredibleSetsTotal"].max()
        check(n2 > n1, f"compat: comparison biosample ALL also sums biosample-group rows ({n2} vs {n1} credible set-biosample pairs)")

        print("\nsetup (offline)")
        check(all(RAW.startswith("https://raw.githubusercontent.com/") and UPSTREAM_COMMIT in RAW for _ in [0])
              and "resources/UKBiobank.ABCGene.anyabc.tsv" in RESOURCE_FILES and len(RESOURCE_FILES) == 5,
              "setup: five resources pinned to the upstream commit")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_config_args(p, full: bool = True) -> None:
    p.add_argument("--config", help="upstream-style config.yml (or .json); flags below override it")
    p.add_argument("--variant-key", help="TSV trait, variant_file (default: Data/GWASE2G/resources/UKBB_variant_key.tsv)")
    p.add_argument("--bg-variants", help="background SNP bed: chr, start, end, rsid (no header)")
    p.add_argument("--partition", help="PartitionCombined.bed")
    p.add_argument("--chr-sizes", help="chromosome sizes TSV")
    p.add_argument("--tss", help="TSS reference bed (4th column = gene universe)")
    p.add_argument("--gene-prioritization-table", help="UKBiobank.ABCGene.anyabc.tsv-style silver standard")
    p.add_argument("--threshold-pip", type=float, help="variants kept when pip > this (default 0.1)")
    p.add_argument("--trait-group", nargs="+", action="extend", help="NAME=trait1,trait2 (repeatable)")
    p.add_argument("--threshold-pval", type=float, help="alpha for CIs and significance (default 0.05)")
    p.add_argument("--num-pops-genes", type=int, help="PoPS rank cut-off (default 2)")
    if full:
        p.add_argument("--methods-table", help="method, boolean, inverse_predictor, pred_name_long, threshold, score_col, color")
        p.add_argument("--predictions-table", help="biosample + one column of prediction paths per method")
        p.add_argument("--comparisons-table", help="name, biosample[group|ALL], trait[group|ALL]")
        p.add_argument("--methods", nargs="+", help="methods to benchmark (default: config or all in the tables)")
        p.add_argument("--biosample-group", nargs="+", action="extend", help="NAME=biosample1,biosample2 (repeatable)")
        p.add_argument("--n-threshold-steps", type=int, help="enrichment-recall curve steps (default 25)")
        p.add_argument("--num-pred-genes", type=int, help="predicted-gene rank cut-off (default 2)")
        p.add_argument("--plot-fixed-scale", action="store_true", help="same heatmap colour scale across methods")


def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent gwas-e2g", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch the small upstream resources from the pinned commit")
    p.add_argument("--traits", nargs="+", help="also fetch these UK Biobank per-trait variant lists ('all' = 94 traits, ~247 MB)")
    p.add_argument("--synapse", action="store_true", help="also download the background SNPs syn52264319 (needs SYNAPSE_AUTH_TOKEN)")
    p.add_argument("--dest", help="destination root (default Data/GWASE2G)")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("validate-config", help="check config, groups and comparisons; list biosamples per method")
    _add_config_args(p)
    p.add_argument("--no-file-checks", action="store_true", help="skip existence checks of prediction / variant files")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("variants", help="PIP + distal-noncoding filters, trait groups, background SNPs (cached for run)")
    _add_config_args(p, full=False)
    p.add_argument("--label", default="gwas_variants")
    p.set_defaults(func=cmd_variants)

    p = sub.add_parser("baseline", help="PoPS and distance-to-TSS gene-linking baselines")
    p.add_argument("--config")
    p.add_argument("--gene-prioritization-table")
    p.add_argument("--tss")
    p.add_argument("--num-pops-genes", type=int)
    p.add_argument("--threshold-pval", type=float)
    p.add_argument("--label", default="gwas_baseline")
    p.set_defaults(func=cmd_baseline)

    p = sub.add_parser("run", help="the whole benchmark")
    _add_config_args(p)
    p.add_argument("--variants-dir", help="output of `variants` (skips variant processing)")
    p.add_argument("--upstream-compat", action="store_true", help="reproduce the upstream quirks (SE formula, span, linking)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="gwas_e2g")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("plot", help="visualisation stage from an existing run directory")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--variants-dir", help="needed when the run used --variants-dir")
    p.add_argument("--plot-fixed-scale", action="store_true")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser("selftest", help="synthetic genome, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
