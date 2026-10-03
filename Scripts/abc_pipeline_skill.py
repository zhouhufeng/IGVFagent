#!/usr/bin/env python3
"""Full Activity-by-Contact enhancer-gene pipeline: peaks -> candidate regions -> neighborhoods -> predictions -> thresholds -> QC (port of broadinstitute/ABC-Enhancer-Gene-Prediction).

Port of https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction (MIT, Copyright (c) 2020 Broad
Institute), pinned at commit 92ac50360231a6bcfd654f3147839d078be445d9 (main, 2026-08-24).  Relationship:
port -- every Snakemake rule (workflow/Snakefile, workflow/rules/*.smk), every script in workflow/scripts
(makeCandidateRegions, peaks, neighborhoods, run.neighborhoods, predict, predictor, hic, filter_predictions,
getVariantOverlap, metrics, grabMetrics, compute_powerlaw_fit_from_hic, makeAverageHiC, extract_avg_hic,
juicebox_dump, make_bedgraph_from_HiC, tools), config/config.yaml and tests/ was read and the algorithms were
re-derived in Python (pandas + numpy + scipy); no code was copied.  bedtools / pyranges / tabix are replaced by
vectorised interval arithmetic; macs2 and juicer_tools are optional binaries (run when on PATH, otherwise the
exact command is printed and the step exits 0); pysam (BAM), pyBigWig (bigWig) and hicstraw (.hic) are
optional imports.  This is the whole pipeline; `igvfagent abc score` (Scripts/abc_skill.py) remains the small
clean-room model-only scorer.

Definitions reproduced from upstream (defaults = config/config.yaml)
  peak calling      macs2 callpeak -f BED (tagAlign) | AUTO -g hs -p 0.1 -n macs2 --shift -75 --extsize 150
                    --nomodel --keep-dup all --call-summits; peaks kept on chrom-sizes chromosomes, sorted in
                    chrom-sizes order (macs2_peaks.narrowPeak.sorted).  Without macs2, --python-fallback runs a
                    MACS-like caller (tags shifted/extended identically, Poisson p against max(lambda_bg,
                    lambda_10kb), one summit per peak) -- an approximation, labelled as such.
  candidate regions reads counted in every peak (each accessibility file; averaged if several; chrX/chrY
                    counts doubled), peaks merged (overlapping or book-ended, max count), the nStrongestPeaks
                    (150000) with most reads kept, their narrowPeak rows -> summit (start + col 10)
                    +/- peakExtendFromSummit (250) clipped to the chromosome, merged, regions touching the
                    blocklist dropped, the TSS include-list (rows on known chromosomes) added, merged.
                    --ignore-summits: whole peak +/- extend, widened to minPeakWidth (500).
  read counting     tagAlign / fragments / BAM: reads overlapping the region by >= 1 bp; bigWig / bedGraph:
                    sum of per-base signal.  Totals: BAM mapped reads in the index; tagAlign lines matching
                    chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY; bigWig sum(length x mean).  fragments are split
                    into two tags (start..mid, mid+1..end) as in docs/usage/scATAC.rst.
  feature columns   <F>.<file>.readCount, .RPM = 1e6 count / total, .readCount.quantile and .RPM.quantile =
                    rank / n, .RPKM = 1e3 RPM / width, .RPKM.quantile; per feature <F>.RPM (mean over files),
                    <F>.RPM.quantile, <F>.RPKM, <F>.RPKM.quantile.  Bedgraphs <Enhancers|Genes|Genes.TSS1kb>.
                    <F>.<file>.CountReads.bedgraph.
  genes             BED6 + Ensembl_ID + gene_type (Ensembl IDs required); TSS = start (+) / end (-);
                    GeneList.TSS1kb = TSS +/- 500 (windows leaving the chromosome dropped); counts over the gene
                    body and the TSS window (suffix .TSS1Kb); PromoterActivityQuantile = pct-rank of
                    (1e-4 + H3K27ac.RPKM.quantile.TSS1Kb) x (1e-4 + <ACC>.RPKM.quantile.TSS1Kb) (accessibility
                    term alone without H3K27ac); Expression = mean of expression tables (NaN without);
                    is_ue = in UbiquitouslyExpressedGenes.
  element class     promoter if overlapping TSS +/- 500, else genic if overlapping a gene body, else
                    intergenic; name = class|chr:start-end; promoterSymbol, genicSymbol.
  quantile norm     rank method: promoters and non-promoters separately, q = rank / n_class, normalized =
                    linear interpolation (with extrapolation) of the K562 reference value at rank (1 - q) n_class,
                    clipped at 0; ATAC uses the reference DHS column if the reference has no ATAC.RPM;
                    --qnorm-method quantile interpolates on q.  No reference -> normalized = RPM.
  activity          activity_base = sqrt(normalized_h3K27ac x normalized_<acc>) with H3K27ac, else
                    normalized_<acc>; activity_base_no_qnorm from raw RPMs.
  expressed genes   isExpressed = Expression >= 1 OR (Expression NaN AND PromoterActivityQuantile >= 0.30).
  pairs             every element overlapping [TSS - 5 Mb, TSS + 5 Mb] (clipped) with |midpoint - TSS| < 5 Mb;
                    isSelfPromoter = class promoter AND start - 500 < TSS < end + 500.
  power law         powerlaw(d) = exp(scale - gamma ln(max(d, 5000) + 1)), gamma 1.024238616787792, scale
                    5.9594510043736655; reference gamma 0.87, reference scale -4.80 + 11.63 x 0.87.
  Hi-C              hic: .hic via hicstraw (observed, SCALE, BP), diagonal = max neighbouring bin within 5 kb,
                    divided by the mean row sum.  juicebox: <dir>/<chr>/<chr>.{KR,INTERSCALE}observed.gz + norm
                    (VC fallback), symmetrised, doubly stochastic unless row sums already 1 (+/-0.001), nonzero
                    diagonal = max(left, right) x tss_hic_contribution / 100, upper triangle,
                    |bin1 - bin2| <= window / resolution.  bedpe: <chr>/<chr>.bedpe.gz, both orientations.
                    avg: <chr>/<chr>.bed.gz (x1 x2 hic_contact).  Genes whose self-promoter rows sum to
                    hic_contact < 0.01 get power-law contact (min distance = resolution); NaN -> 0.
                    hic_contact_pl_scaled = hic_contact x powerlaw_contact_reference / powerlaw_contact
                    (--scale-hic-using-powerlaw, on in the config); hic_pseudocount = min(powerlaw_contact,
                    powerlaw(5000)); hic_contact_pl_scaled_adj = scaled + pseudocount.
  ABC               ABC.Score.Numerator = activity_base_enh x (hic_contact_pl_scaled_adj | powerlaw_contact);
                    ABC.Score = numerator / sum over (TargetGene, TargetGeneTSS); self-promoters -> 1;
                    powerlaw.Score likewise with powerlaw_contact.  Written %.6f, NaN, to
                    EnhancerPredictionsAllPutative.tsv.gz (TargetGeneIsExpressed) and
                    ...NonExpressedGenes.tsv.gz; ForVariantOverlap.shrunk150bp = (score > 0.015 & non-promoter)
                    | (promoter & score > 0.1), distance <= 2 Mb, 150 bp trimmed from both ends.
  thresholds        explicit threshold, else reference/abc_thresholds.tsv by (accessibility, H3K27ac?, Hi-C
                    type powerlaw | avg | intact_hic), else 0.02; kept = score > threshold, promoters removed
                    except self-promoters (include_self_promoter True); EnhancerPredictionsFull_threshold{t}
                    _self_promoter.tsv / .bedpe.gz, EnhancerPredictions_... (slim), GenePredictionStats_... .
  QC                QCSummary: MedianEnhPerGene, StdEnhPerGene, MedianGenePerEnh, StdGenePerEnh,
                    MeanEnhPerGene, MeanGenePerEnh, MedianEGDist, MeanEGDist, StdEGDist, NumEnhancersPerChrom,
                    MeanEnhancersPerChrom, NumEnhancers, EG10th, EG90th, NumPeaks, MedWidth, MeanWidth,
                    StdWidth, NumCandidate, MedWidthCandidate, MeanWidthCandidate, StdWidthCandidate,
                    counts{Enhancers,GeneTSS,Genes}_<F>; EnhancerPerGene.tsv, GenesPerEnhancer.tsv,
                    EnhancerGenePairsPerChrom.txt; distribution and power-law plots.
  power-law fit     ln(sum contact at distance / n entries at distance == resolution + 1e-6) regressed on
                    ln(distance + 1e-6) over 5 kb..1 Mb; hic_gamma = -slope, hic_scale = intercept.
  average Hi-C      per cell type KR matrix made doubly stochastic (no diagonal correction), rescaled by
                    powerlaw_ref / powerlaw_celltype (ref gamma 0.876, scale 5.41), outer-joined (absent = 0,
                    NaN kept), mean of non-NaN, NaN when fewer than 3 cell types; zero entries dropped.

Deliberate deviations (upstream behaviour reproducible with --upstream-compat)
  * UbiquitouslyExpressedGenes.txt is read with a header upstream, silently dropping its first gene (AARS);
    the port reads it headerless.
  * determine_threshold maps HiC_type avg to "avg_hic", which never matches abc_thresholds.tsv ("avg"), so
    average-Hi-C runs fall back to 0.02; the port matches "avg" (0.016 / 0.012 ...).
  * compute_powerlaw_fit_from_hic, hic_type avg: distance is |bin1 - bin2| (bins, not bp) and the window is
    not applied; the port uses bp and the 5 kb..1 Mb window.
  Fixed without a flag (upstream crashes or is dead code): run_qnorm without separate promoters indexes
  qnorm["enh_class" == "any"]; compute_powerlaw_fit passes an unknown interpolate_nan argument; the
  --hic_is_doubly_stochastic flag is never forwarded; make_bedgraph_from_HiC imports a HiC class hic.py does
  not define (rebuilt here on the juicebox loader); --ignore-summits widening can go below 0 (clipped).
  QC NumPeaks describes the candidate regions (the rule passes them as --macs_peaks), unless --narrowpeak is
  given.  Upstream's expected test outputs predate the current predict.py column names; this port writes the
  current ones (activity_base_enh, activity_base_squared_enh, normalized_<acc>_enh, normalized_<acc>_prom).

Subcommands
  setup              fetch reference/ (thresholds, K562 qnorm reference, ubiquitous genes, hg38/hg19 gene
                     bounds, TSS500bp, chrom sizes, blocklists) from the pinned commit (--dry-run lists them).
  call-peaks         macs2 wrapper (+ narrowPeak filter/sort); --python-fallback when macs2 is absent.
  fragments-to-tagalign  10x fragments -> tagAlign.gz (two tags per fragment), as the scATAC walkthrough.
  candidate-regions  makeCandidateRegions.py.
  neighborhoods      run.neighborhoods.py: GeneList.txt, EnhancerList.txt, count bedgraphs.
  predict            predict.py: AllPutative tables + variant-overlap file (powerlaw / hic / juicebox / bedpe / avg).
  filter             filter_predictions.py with automatic threshold choice.
  variant-overlap    getVariantOverlap.py on an AllPutative file.
  qc                 grabMetrics.py + metrics.py: QCSummary, per-gene/per-enhancer tables, plots.
  powerlaw-fit       compute_powerlaw_fit_from_hic.py.
  average-hic        makeAverageHiC.py.
  split-avg-hic      extract_avg_hic.py (genome-wide average Hi-C bed.gz -> AvgHiC/<chr>/<chr>.bed.gz).
  juicebox-dump      juicebox_dump.py (juicer_tools dump KR observed / norm per chromosome).
  hic-bedgraph       make_bedgraph_from_HiC.py (per-gene Hi-C row bedgraphs).
  compare            tests/test_full_abc_run.py comparison of two prediction files on the tested columns.
  run                the Snakefile for one biosample or a biosamples table (upstream's config columns).
  selftest           synthetic chromosome with planted enhancers, a power-law Hi-C map with a planted loop;
                     every subcommand, every definition asserted.

Usage:
    igvfagent abc-pipeline setup [--dry-run]
    igvfagent abc-pipeline run --biosample K562 --dhs dnase.bam --h3k27ac h3k27ac.bam --genes CollapsedGeneBounds.hg38.bed --tss CollapsedGeneBounds.hg38.TSS500bp.bed --chrom-sizes GRCh38_EBV.no_alt.chrom.sizes.tsv --blocklist GRCh38_unified_blacklist.bed --label k562
    igvfagent abc-pipeline run --biosamples-table config_biosamples.tsv --label batch
    igvfagent abc-pipeline candidate-regions --narrowpeak macs2_peaks.narrowPeak.sorted --accessibility atac.tagAlign.gz --chrom-sizes sizes.tsv --includelist tss.bed --blocklist blocklist.bed
    igvfagent abc-pipeline predict --enhancers EnhancerList.txt --genes GeneList.txt --chrom-sizes sizes.tsv --accessibility-feature DHS --hic-file ENCFF621AIY.hic --hic-type hic --hic-resolution 5000
    igvfagent abc-pipeline filter --pred-file EnhancerPredictionsAllPutative.tsv.gz --accessibility-feature DHS --has-h3k27ac --hic-type hic
    igvfagent abc-pipeline powerlaw-fit --hic-dir juicebox_dir --hic-type juicebox
    igvfagent abc-pipeline selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "ABCPipeline"
REF_DIR = ROOT / "Data" / "ABCPipeline" / "reference"

UPSTREAM_REPO = "broadinstitute/ABC-Enhancer-Gene-Prediction"
UPSTREAM_COMMIT = "92ac50360231a6bcfd654f3147839d078be445d9"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
REFERENCE_FILES = {
    "reference/abc_thresholds.tsv": "ABC score thresholds at 70% CRISPR recall",
    "reference/EnhancersQNormRef.K562.txt": "K562 quantile-normalisation reference (DHS.RPM, H3K27ac.RPM)",
    "reference/UbiquitouslyExpressedGenes.txt": "ubiquitously expressed genes",
    "reference/hg38/CollapsedGeneBounds.hg38.bed": "hg38 gene bounds (BED6 + Ensembl_ID + gene_type)",
    "reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed": "hg38 TSS +/- 250 bp include-list",
    "reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv": "hg38 chromosome sizes (no alt)",
    "reference/hg38/GRCh38_unified_blacklist.bed": "hg38 blocklist",
    "reference/hg19/CollapsedGeneBounds.hg19.bed": "hg19 gene bounds",
    "reference/hg19/CollapsedGeneBounds.hg19.TSS500bp.bed": "hg19 TSS include-list",
    "reference/hg19/chrom_sizes.tsv": "hg19 chromosome sizes",
    "reference/hg19/wgEncodeHg19ConsensusSignalArtifactRegions.bed": "hg19 blocklist",
}

# config/config.yaml
DEFAULTS = {
    "macs_pval": 0.1, "macs_genome_size": "hs", "peak_extend": 250, "n_strongest": 150000, "min_peak_width": 500,
    "hic_gamma": 1.024238616787792, "hic_scale": 5.9594510043736655, "hic_pseudocount_distance": 5000,
    "hic_gamma_reference": 0.87, "window": 5_000_000, "tss_slop": 500, "tss_hic_contribution": 100.0,
    "expression_cutoff": 1.0, "promoter_activity_quantile_cutoff": 0.30, "score_column": "ABC.Score",
    "default_threshold": 0.02,
}
# reference/abc_thresholds.tsv at the pinned commit
ABC_THRESHOLDS = [
    ("DHS", True, "intact_hic", 0.027), ("DHS", False, "intact_hic", 0.024), ("DHS", True, "avg", 0.016),
    ("DHS", True, "powerlaw", 0.017), ("DHS", False, "avg", 0.016), ("ATAC", True, "intact_hic", 0.025),
    ("DHS", False, "powerlaw", 0.016), ("ATAC", True, "avg", 0.016), ("ATAC", True, "powerlaw", 0.016),
    ("ATAC", False, "intact_hic", 0.021), ("ATAC", False, "avg", 0.012), ("ATAC", False, "powerlaw", 0.013),
]
TAGALIGN_TOTAL_RE = re.compile(r"chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY")
GENOME_SIZES = {"hs": 2.7e9, "mm": 1.87e9, "ce": 9e7, "dm": 1.2e8}
BIN = 1 << 16

BLUE, ORANGE, GREEN, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"

log = logging.getLogger("abc_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"abc_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir_for(args, default_label: str) -> Path:
    if getattr(args, "outdir", None):
        d = Path(args.outdir)
        d.mkdir(parents=True, exist_ok=True)
        return d
    return run_dir(getattr(args, "label", None) or default_label)


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("abc-pipeline needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


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
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path, pdf=None) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def emit(kind: str, path) -> None:
    print(f"{kind}: {path}")


def md_table(headers: "List[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 40) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 4) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        if isinstance(x, (int,)) and not isinstance(x, bool):
            return f"{x:,}"
        return f"{float(x):.{nd}g}"
    except (TypeError, ValueError):
        return str(x)


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_report(d: Path, title: str, sections: "List[Tuple[str, str]]", summary: dict) -> None:
    lines = [f"# {title}", "", f"_abc-pipeline, port of {UPSTREAM_REPO} @ {UPSTREAM_COMMIT[:10]}_", ""]
    for head, body in sections:
        lines += [f"## {head}", "", body, ""]
    (d / "report.md").write_text("\n".join(lines))
    (d / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2))
    emit("Report", d / "report.md")
    emit("JSON", d / "summary.json")


def write_tsv(df, path: Path, float_format: "Optional[str]" = "%.6f", header=True, na_rep="", quiet=False, index=False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, header=header, float_format=float_format, na_rep=na_rep,
              compression="gzip" if str(path).endswith(".gz") else None)
    if not quiet:
        emit("TSV", path)
    return path


def write_params(d: Path, name: str, params: dict) -> Path:
    p = Path(d) / name
    p.write_text("".join(f"{k} {v}\n" for k, v in params.items()))
    return p


def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "t", "1", "yes", "y"}


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v)) or str(v).strip() in ("", "nan", "None")


# ---------------------------------------------------------------------------
# BED arithmetic (bedtools sort / merge / slop / intersect replacements)
# ---------------------------------------------------------------------------

def read_chrom_sizes(path) -> "Dict[str, int]":
    sizes: "Dict[str, int]" = {}
    with _open_text(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and not line.startswith("#"):
                sizes[f[0]] = int(f[1])
    return sizes


def chrom_sizes_bed(sizes: "Dict[str, int]"):
    pd = _pd()
    return pd.DataFrame({"chr": list(sizes), "start": 0, "end": list(sizes.values())})


BED_EXTRA = ["name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]


def read_bed(path, extra_cols: "Optional[List[str]]" = None):
    """neighborhoods.read_bed: BED3 + named extra columns, 'track' line skipped, '#' comments, empty columns dropped."""
    pd = _pd()
    extra = BED_EXTRA if extra_cols is None else extra_cols
    with _open_text(path) as fh:
        first = fh.readline()
    skip = 1 if "track" in first else 0
    names = ["chr", "start", "end"] + list(extra)
    df = pd.read_csv(path, sep="\t", names=names, header=None, skiprows=skip, comment="#", low_memory=False)
    df = df.dropna(axis=1, how="all")
    if len(df):
        df["chr"] = df["chr"].astype(str)
        df["start"] = df["start"].astype("int64")
        df["end"] = df["end"].astype("int64")
    return df


def sort_bed(df, order: "Optional[Dict[str, int]]" = None, cols=("chr", "start")):
    """bedtools sort (-faidx): chromosome order of the sizes file (lexicographic without), then start."""
    if len(df) == 0:
        return df.reset_index(drop=True)
    c, s = cols
    if order:
        rank = {k: i for i, k in enumerate(order)}
        key = df[c].map(lambda x: rank.get(x, len(rank)))
    else:
        key = df[c]
    tmp = df.assign(_k=key.values)
    return tmp.sort_values(["_k", s], kind="mergesort").drop(columns="_k").reset_index(drop=True)


def merge_intervals(df, order=None, value_col: "Optional[str]" = None, how: str = "max"):
    """bedtools merge [-c 4 -o max]: overlapping or book-ended intervals are merged."""
    pd, np = _pd(), _np()
    if len(df) == 0:
        return pd.DataFrame(columns=["chr", "start", "end"] + ([value_col] if value_col else []))
    df = sort_bed(df, order)
    out = []
    for chrom, sub in df.groupby("chr", sort=False):
        s = sub["start"].to_numpy(np.int64)
        e = sub["end"].to_numpy(np.int64)
        # running max of ends; a new block starts where start > running max end of everything before
        cmax = np.maximum.accumulate(e)
        new = np.ones(len(s), dtype=bool)
        new[1:] = s[1:] > cmax[:-1]
        gid = np.cumsum(new) - 1
        g = pd.DataFrame({"g": gid, "start": s, "end": e})
        if value_col:
            g["v"] = sub[value_col].to_numpy()
        agg = {"start": "min", "end": "max"}
        if value_col:
            agg["v"] = how
        m = g.groupby("g", sort=True).agg(agg)
        piece = pd.DataFrame({"chr": chrom, "start": m["start"].to_numpy(), "end": m["end"].to_numpy()})
        if value_col:
            piece[value_col] = m["v"].to_numpy()
        out.append(piece)
    return pd.concat(out, ignore_index=True)


def clip_to_genome(df, sizes: "Dict[str, int]"):
    np = _np()
    df = df.copy()
    df["start"] = np.maximum(df["start"].to_numpy(np.int64), 0)
    ends = df["chr"].map(sizes).to_numpy(dtype=float)
    df["end"] = np.minimum(df["end"].to_numpy(np.int64), np.where(np.isnan(ends), np.inf, ends)).astype(np.int64)
    return df


def overlaps_any(a, b):
    """Boolean mask over rows of `a`: overlaps (>= 1 bp, half-open) any interval of `b` (bedtools intersect -u)."""
    np = _np()
    mask = np.zeros(len(a), dtype=bool)
    if len(a) == 0 or len(b) == 0:
        return mask
    a = a.reset_index(drop=True)
    for chrom, bsub in b.groupby("chr", sort=False):
        ai = np.flatnonzero((a["chr"] == chrom).to_numpy())
        if ai.size == 0:
            continue
        m = merge_intervals(bsub[["chr", "start", "end"]])
        bs, be = m["start"].to_numpy(np.int64), m["end"].to_numpy(np.int64)
        s = a["start"].to_numpy(np.int64)[ai]
        e = np.maximum(a["end"].to_numpy(np.int64)[ai], s + 1)
        # first merged block with end > s; overlap iff its start < e
        k = np.searchsorted(be, s, side="right")
        ok = k < bs.size
        hit = np.zeros(ai.size, dtype=bool)
        hit[ok] = bs[k[ok]] < e[ok]
        mask[ai] = hit
    return mask


def on_genome(df, sizes: "Dict[str, int]"):
    """bedtools intersect -u -b chrom_sizes.bed: rows on a listed chromosome that start inside it."""
    np = _np()
    size = df["chr"].map(sizes)
    return (size.notna() & (df["start"].to_numpy(np.int64) < size.fillna(0).to_numpy()) & (df["end"].to_numpy(np.int64) > 0)).to_numpy()


# ---------------------------------------------------------------------------
# Read counting (bedtools coverage / pysam / pyBigWig replacements)
# ---------------------------------------------------------------------------

_READ_CACHE: "Dict[str, Any]" = {}


def file_kind(path) -> str:
    name = os.path.basename(str(path))
    low = name.lower()
    if name.endswith(".bam"):
        return "bam"
    if "tagAlign" in name or "tagalign" in low:
        return "tagalign"
    if "fragments" in low:
        return "fragments"
    if low.endswith((".bw", ".bigwig")):
        return "bigwig"
    if low.endswith((".bedgraph", ".bedgraph.gz", ".bg", ".bg.gz")):
        return "bedgraph"
    raise ValueError(f"File {path} name format doesn't match bam, tagAlign, fragments, bigWig or bedGraph")


def fragments_to_tags(df):
    """docs/usage/scATAC.rst: mid = int((s+e)/2); tags (s, mid, +) and (mid+1, e, -)."""
    pd, np = _pd(), _np()
    s = df["start"].to_numpy(np.int64)
    e = df["end"].to_numpy(np.int64)
    mid = (s + e) // 2
    n = len(df)
    tags = pd.DataFrame({
        "chr": np.repeat(df["chr"].to_numpy(), 2),
        "start": np.column_stack([s, mid + 1]).ravel(),
        "end": np.column_stack([mid, e]).ravel(),
        "strand": np.tile(np.array(["+", "-"]), n),
    })
    return tags


def load_tags(path):
    """tagAlign / fragments / BAM -> per-chromosome arrays (starts, ends, strand '+'?); cached per file."""
    pd, np = _pd(), _np()
    key = str(Path(path).resolve())
    st = os.stat(key)
    ck = (key, st.st_mtime, st.st_size)
    if ck in _READ_CACHE:
        return _READ_CACHE[ck]
    kind = file_kind(path)
    if kind == "bam":
        try:
            import pysam  # type: ignore
        except ImportError as e:
            raise SystemExit("BAM input needs pysam (pip install pysam), or convert to tagAlign") from e
        chroms, starts, ends, plus = [], [], [], []
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_unmapped:
                    continue
                chroms.append(r.reference_name); starts.append(r.reference_start); ends.append(r.reference_end)
                plus.append(not r.is_reverse)
        df = pd.DataFrame({"chr": chroms, "start": starts, "end": ends, "strand": np.where(plus, "+", "-")})
    else:
        raw = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=None, low_memory=False)
        raw = raw.rename(columns={0: "chr", 1: "start", 2: "end"})
        raw["chr"] = raw["chr"].astype(str)
        if kind == "fragments":
            df = fragments_to_tags(raw)
        else:
            strand = raw[5] if 5 in raw.columns else "+"
            df = pd.DataFrame({"chr": raw["chr"], "start": raw["start"].astype(np.int64),
                               "end": raw["end"].astype(np.int64), "strand": strand})
    per = {}
    for chrom, sub in df.groupby("chr", sort=False):
        s = sub["start"].to_numpy(np.int64)
        e = sub["end"].to_numpy(np.int64)
        p = (sub["strand"].to_numpy() != "-")
        per[chrom] = {"start": s, "end": e, "plus": p, "s_sorted": np.sort(s), "e_sorted": np.sort(e)}
    total_regex = int(sum(len(v["start"]) for c, v in per.items() if TAGALIGN_TOTAL_RE.search(c)))
    out = {"per": per, "kind": kind, "n": int(len(df)), "total_regex": total_regex}
    _READ_CACHE[ck] = out
    return out


def _load_bedgraph(path):
    pd, np = _pd(), _np()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3], names=["chr", "start", "end", "v"])
    df = df[~df["chr"].astype(str).str.startswith(("track", "browser"))]
    per = {}
    for chrom, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values("start")
        s = sub["start"].to_numpy(np.int64); e = sub["end"].to_numpy(np.int64); v = sub["v"].to_numpy(float)
        cum = np.concatenate([[0.0], np.cumsum(v * (e - s))])
        per[str(chrom)] = (s, e, v, cum)
    return per


def _bedgraph_integral(track, x):
    np = _np()
    s, e, v, cum = track
    k = np.searchsorted(s, x, side="left")          # intervals with start < x
    total = cum[k]
    last = k - 1
    has = last >= 0
    over = np.zeros_like(x, dtype=float)
    li = last[has]
    over[has] = np.clip(e[li] - x[has], 0, None) * v[li]
    return total - over


def count_reads_in_regions(target, regions):
    """Counts per region (bedtools coverage -counts / pysam count / pyBigWig sum), X/Y doubled afterwards by caller."""
    np = _np()
    kind = file_kind(target)
    counts = np.zeros(len(regions), dtype=float)
    chr_arr = regions["chr"].astype(str).to_numpy()
    st = regions["start"].to_numpy(np.int64)
    en = regions["end"].to_numpy(np.int64)
    if kind in ("tagalign", "fragments", "bam"):
        per = load_tags(target)["per"]
        for chrom in np.unique(chr_arr):
            idx = np.flatnonzero(chr_arr == chrom)
            t = per.get(chrom)
            if t is None:
                continue
            counts[idx] = (np.searchsorted(t["s_sorted"], en[idx], side="left")
                           - np.searchsorted(t["e_sorted"], st[idx], side="right"))
        return counts
    if kind == "bedgraph":
        per = _load_bedgraph(target)
        for chrom in np.unique(chr_arr):
            idx = np.flatnonzero(chr_arr == chrom)
            t = per.get(chrom)
            if t is None:
                continue
            counts[idx] = _bedgraph_integral(t, en[idx]) - _bedgraph_integral(t, st[idx])
        return counts
    try:
        import pyBigWig  # type: ignore
    except ImportError as e:
        raise SystemExit("bigWig input needs pyBigWig (pip install pyBigWig); a bedGraph works without it") from e
    bw = pyBigWig.open(str(target))
    chroms = bw.chroms()
    for i, (c, s, e) in enumerate(zip(chr_arr, st, en)):
        if c in chroms:
            counts[i] = bw.stats(c, int(s), int(e), type="sum", exact=True)[0] or 0
    bw.close()
    return counts


def count_total(target) -> float:
    kind = file_kind(target)
    if kind in ("tagalign", "fragments"):
        n = load_tags(target)["total_regex"]
        if n <= 0:
            raise ValueError(f"no reads on chr1-22/X/Y in {target}")
        return float(n)
    if kind == "bam":
        try:
            import pysam  # type: ignore
            with pysam.AlignmentFile(str(target), "rb") as bam:
                n = bam.mapped
        except (ImportError, ValueError):
            n = load_tags(target)["n"]
        if not n > 0:
            raise ValueError("Error counting BAM file: count <= 0")
        return float(n)
    if kind == "bedgraph":
        per = _load_bedgraph(target)
        return float(sum(t[3][-1] for t in per.values()))
    import pyBigWig  # type: ignore
    bw = pyBigWig.open(str(target))
    tot = sum(l * (bw.stats(c, 0, l, "mean", exact=True)[0] or 0) for c, l in bw.chroms().items())
    bw.close()
    return float(tot)


def double_sex_chrom_counts(df, col: str):
    """neighborhoods.double_sex_chrom_counts: chromosomes whose name ends in X or Y count twice."""
    m = df["chr"].astype(str).str[-1:].isin(["X", "Y"])
    df.loc[m, col] = df.loc[m, col] * 2
    return df


def run_count_reads(target, regions, output: "Optional[Path]" = None):
    """regions (chr,start,end) -> DataFrame chr,start,end,count written as a headerless bedgraph."""
    np = _np()
    out = regions[["chr", "start", "end"]].copy().reset_index(drop=True)
    c = count_reads_in_regions(target, out)
    if file_kind(target) in ("tagalign", "fragments", "bam"):
        c = c.astype(np.int64)
    out["count"] = c
    out = double_sex_chrom_counts(out, "count")
    if output is not None:
        out.to_csv(output, sep="\t", header=False, index=False)
    return out


# ---------------------------------------------------------------------------
# Peak calling (rules/macs2.smk)
# ---------------------------------------------------------------------------

NARROWPEAK_COLS = ["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]


def macs2_command(accessibility: "List[str]", outdir: Path, pval: float = 0.1, genome_size: str = "hs", name: str = "macs2") -> "List[str]":
    fmt = "BED" if any("tagAlign" in os.path.basename(a) for a in accessibility) else "AUTO"
    return ["macs2", "callpeak", "-f", fmt, "-g", str(genome_size), "-p", str(pval), "-n", name, "--shift", "-75",
            "--extsize", "150", "--nomodel", "--keep-dup", "all", "--call-summits", "--outdir", str(outdir), "-t", *accessibility]


def read_narrowpeak(path):
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    df = df.iloc[:, :10]
    df.columns = NARROWPEAK_COLS[: df.shape[1]]
    df["chr"] = df["chr"].astype(str)
    return df


def sort_narrowpeak(narrowpeak, sizes: "Dict[str, int]", out: Path) -> Path:
    """rule sort_narrowpeaks: keep peaks on the chrom-sizes chromosomes, sort in that order."""
    df = read_narrowpeak(narrowpeak)
    df = df[on_genome(df, sizes)]
    df = sort_bed(df, sizes)
    df.to_csv(out, sep="\t", header=False, index=False)
    return out


def python_call_peaks(accessibility: "List[str]", sizes: "Dict[str, int]", outdir: Path, pval: float = 0.1,
                      genome_size="hs", name: str = "macs2", bin_size: int = 10, llocal: int = 10000,
                      shift: int = -75, extsize: int = 150) -> Path:
    """MACS-like fallback: tags shifted by `shift` and extended to `extsize` (both strands centred on the cut
    site for -75/150), pileup at bin centres, Poisson -log10 p against max(lambda_bg, lambda_llocal), runs of
    significant bins (gaps <= tag length merged, width >= extsize) become peaks, one summit per peak."""
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    gsize = GENOME_SIZES.get(str(genome_size), None) or float(genome_size)
    tagsets = [load_tags(a)["per"] for a in accessibility]
    n_total = sum(len(v["start"]) for ts in tagsets for v in ts.values())
    lam_bg = n_total * extsize / gsize
    thr = -math.log10(pval)
    rows = []
    for chrom, size in sizes.items():
        lo_l, hi_l, five_l, tlen = [], [], [], []
        for ts in tagsets:
            t = ts.get(chrom)
            if t is None:
                continue
            five = np.where(t["plus"], t["start"], t["end"])
            lo = np.where(t["plus"], five + shift, five - shift - extsize)
            lo_l.append(lo); hi_l.append(lo + extsize); five_l.append(five); tlen.append(t["end"] - t["start"])
        if not lo_l:
            continue
        lo = np.sort(np.concatenate(lo_l)); hi = np.sort(np.concatenate(hi_l)); five = np.sort(np.concatenate(five_l))
        max_gap = int(np.median(np.concatenate(tlen))) if tlen else 50
        centers = np.arange(bin_size // 2, size, bin_size)
        pile = np.searchsorted(lo, centers, side="right") - np.searchsorted(hi, centers, side="right")
        nwin = np.searchsorted(five, centers + llocal // 2) - np.searchsorted(five, centers - llocal // 2)
        lam = np.maximum(lam_bg, nwin * extsize / llocal)
        mlp = np.zeros(centers.size)
        pos = pile > 0
        mlp[pos] = -stats.poisson.logsf(pile[pos] - 1, lam[pos]) / math.log(10)
        sig = mlp >= thr
        if not sig.any():
            continue
        idx = np.flatnonzero(sig)
        breaks = np.flatnonzero(np.diff(idx) * bin_size > max_gap + bin_size)
        starts = np.concatenate([[0], breaks + 1]); ends = np.concatenate([breaks, [idx.size - 1]])
        for a, b in zip(starts, ends):
            i0, i1 = idx[a], idx[b]
            ps, pe = int(i0 * bin_size), int(min(size, (i1 + 1) * bin_size))
            if pe - ps < extsize:
                continue
            seg = pile[i0:i1 + 1]
            top = np.flatnonzero(seg == seg.max())
            k = i0 + int(top[len(top) // 2])
            rows.append((chrom, ps, pe, float(pile[k] / lam[k]), float(mlp[k]), int(centers[k]) - ps))
    df = pd.DataFrame(rows, columns=["chr", "start", "end", "signalValue", "pValue", "peak"])
    if len(df):
        p = np.power(10.0, -df["pValue"].to_numpy())
        order = np.argsort(p)
        q = np.empty_like(p)
        ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
        q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
        df["qValue"] = -np.log10(np.clip(q, 1e-300, 1))
    else:
        df["qValue"] = []
    df.insert(3, "name", [f"{name}_peak_{i + 1}" for i in range(len(df))])
    df.insert(4, "score", (np.minimum(1000, 10 * df["qValue"].to_numpy())).astype(int) if len(df) else [])
    df.insert(5, "strand", ".")
    df = df[NARROWPEAK_COLS]
    out = Path(outdir) / f"{name}_peaks.narrowPeak"
    df.to_csv(out, sep="\t", header=False, index=False, float_format="%.5f")
    summ = pd.DataFrame({"chr": df["chr"], "s": df["start"] + df["peak"], "e": df["start"] + df["peak"] + 1,
                         "name": df["name"], "score": df["pValue"]})
    summ.to_csv(Path(outdir) / f"{name}_summits.bed", sep="\t", header=False, index=False, float_format="%.5f")
    return out


def call_peaks(accessibility: "List[str]", sizes_path, outdir: Path, pval: float, genome_size: str,
               python_fallback: bool = False, force_python: bool = False) -> "Tuple[Optional[Path], str]":
    """Returns (sorted narrowPeak or None, method)."""
    outdir.mkdir(parents=True, exist_ok=True)
    sizes = read_chrom_sizes(sizes_path)
    cmd = macs2_command(accessibility, outdir, pval, genome_size)
    raw = outdir / "macs2_peaks.narrowPeak"
    method = "macs2"
    if shutil.which("macs2") and not force_python:
        print("Running: " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    elif python_fallback or force_python:
        print("macs2 not on PATH; running the Python MACS-like fallback (approximation). The macs2 command would be:")
        print("  " + " ".join(cmd))
        python_call_peaks(accessibility, sizes, outdir, pval, genome_size)
        method = "python-fallback"
    else:
        print("macs2 is not on PATH. Command that would run:")
        print("  " + " ".join(cmd))
        print("Install MACS2 (pip install MACS2) or rerun with --python-fallback.")
        return None, "not-run"
    out = outdir / "macs2_peaks.narrowPeak.sorted"
    sort_narrowpeak(raw, sizes, out)
    emit("Wrote", raw)
    emit("Wrote", out)
    return out, method


# ---------------------------------------------------------------------------
# Candidate regions (makeCandidateRegions.py / peaks.py)
# ---------------------------------------------------------------------------

def count_reads_over_peaks(narrowpeak, accessibility: "List[str]", outdir: Path):
    """peaks.get_read_counts: per-file <peaks>.<file>.Counts.bed, averaged when several files."""
    pd = _pd()
    peaks = read_narrowpeak(narrowpeak)[["chr", "start", "end"]]
    base = os.path.basename(str(narrowpeak))
    outs = []
    for acc in accessibility:
        o = Path(outdir) / f"{base}.{os.path.basename(str(acc))}.Counts.bed"
        outs.append(run_count_reads(acc, peaks, o))
    if len(outs) == 1:
        return outs[0], Path(outdir) / f"{base}.{os.path.basename(str(accessibility[0]))}.Counts.bed"
    avg = outs[0].copy()
    avg["count"] = sum(o["count"].astype(float) for o in outs) / len(outs)
    p = Path(outdir) / f"{base}.averageAccessibility.Counts.bed"
    avg.to_csv(p, sep="\t", header=False, index=False)
    return avg, p


def top_peaks_by_count(counts, n: int, order):
    """bedtools sort | bedtools merge -c 4 -o max | sort -nr -k 4 | head -n N."""
    np = _np()
    m = merge_intervals(counts, order, value_col="count", how="max")
    if len(m) == 0:
        return m
    line = m["chr"].astype(str) + "\t" + m["start"].astype(str) + "\t" + m["end"].astype(str) + "\t" + m["count"].astype(str)
    # GNU sort -nr -k4: numeric on field 4 descending, ties broken by the whole line, also reversed
    o = np.lexsort((line.to_numpy(), m["count"].to_numpy(dtype=float)))[::-1]
    return m.iloc[o[: int(n)]].reset_index(drop=True)


def make_candidate_regions(narrowpeak, accessibility: "List[str]", sizes_path, outdir: Path, includelist=None,
                           blocklist=None, n_strongest: int = 150000, peak_extend: int = 250,
                           ignore_summits: bool = False, min_peak_width: int = 500) -> "Dict[str, Any]":
    pd, np = _pd(), _np()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sizes = read_chrom_sizes(sizes_path)
    write_params(outdir, "params.txt", {
        "narrowPeak": narrowpeak, "accessibility": list(map(str, accessibility)), "chrom_sizes": sizes_path,
        "outDir": outdir, "nStrongestPeaks": n_strongest, "peakExtendFromSummit": peak_extend,
        "ignoreSummits": ignore_summits, "minPeakWidth": min_peak_width,
        "regions_includelist": includelist or "", "regions_blocklist": blocklist or ""})
    outfile = outdir / (os.path.basename(str(narrowpeak)) + ".candidateRegions.bed")
    counts, counts_path = count_reads_over_peaks(narrowpeak, accessibility, outdir)
    top = top_peaks_by_count(counts, n_strongest, sizes)
    peaks = read_narrowpeak(narrowpeak)
    keep = overlaps_any(peaks, top)
    sel = peaks[keep].copy()
    if not ignore_summits:
        summit = sel["start"].to_numpy(np.int64) + sel["peak"].to_numpy(np.int64)
        reg = pd.DataFrame({"chr": sel["chr"].to_numpy(), "start": summit - peak_extend, "end": summit + peak_extend})
        reg = clip_to_genome(reg, sizes)
    else:
        reg = pd.DataFrame({"chr": sel["chr"].to_numpy(), "start": sel["start"].to_numpy(np.int64) - peak_extend,
                            "end": sel["end"].to_numpy(np.int64) + peak_extend})
        reg = clip_to_genome(reg, sizes)
        width = (reg["end"] - reg["start"]).to_numpy()
        short = width < min_peak_width
        pad = ((min_peak_width - width) // 2).astype(np.int64)
        reg.loc[short, "start"] = np.maximum(reg.loc[short, "start"].to_numpy() - pad[short], 0)
        reg.loc[short, "end"] = reg.loc[short, "end"].to_numpy() + pad[short]
    reg = merge_intervals(reg, sizes)
    n_before_block = len(reg)
    if blocklist:
        bl = read_bed(blocklist, extra_cols=[])
        reg = reg[~overlaps_any(reg, bl)]
    n_blocked = n_before_block - len(reg)
    n_include = 0
    if includelist:
        inc = read_bed(includelist, extra_cols=["name", "score", "strand", "x7", "x8"])[["chr", "start", "end"]]
        inc = inc[on_genome(inc, sizes)]
        n_include = len(inc)
        reg = pd.concat([reg[["chr", "start", "end"]], inc], ignore_index=True)
    reg = merge_intervals(reg, sizes)[["chr", "start", "end"]]
    reg.to_csv(outfile, sep="\t", header=False, index=False)
    emit("Wrote", counts_path)
    emit("Wrote", outfile)
    widths = (reg["end"] - reg["start"]).to_numpy()
    return {"candidate_regions": str(outfile), "counts_file": str(counts_path), "n_peaks": int(len(peaks)),
            "n_merged_peaks": int(len(merge_intervals(counts, sizes))), "n_top_peaks": int(len(top)),
            "n_peak_rows_selected": int(keep.sum()), "n_blocklisted": int(n_blocked), "n_includelist_rows": int(n_include),
            "n_candidate_regions": int(len(reg)), "median_width": float(np.median(widths)) if len(widths) else float("nan"),
            "total_bp": int(widths.sum())}


# ---------------------------------------------------------------------------
# Neighborhoods (run.neighborhoods.py / neighborhoods.py)
# ---------------------------------------------------------------------------

def read_gene_bed_file(path):
    """BED6 + Ensembl_ID + gene_type; every Ensembl_ID must contain 'EN'."""
    pd = _pd()
    with _open_text(path) as fh:
        first = fh.readline()
    skip = 1 if "track" in first else 0
    cols = ["chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"]
    df = pd.read_csv(path, sep="\t", names=cols, header=None, skiprows=skip, comment="#", low_memory=False)
    ens = df["Ensembl_ID"]
    if ens.isna().all() or not ens.astype(str).str.contains("EN").all():
        raise SystemExit(f"Gene file {path} doesn't follow the correct format with Ensembl info (BED6 + Ensembl_ID + gene_type)")
    df["chr"] = df["chr"].astype(str)
    return df


def get_tss(bed):
    tss = bed["start"].copy()
    minus = bed["strand"] == "-"
    tss.loc[minus] = bed.loc[minus, "end"]
    return tss


def process_gene_bed(bed, name_cols: str = "symbol", main_name: str = "symbol", sizes: "Optional[Dict[str, int]]" = None,
                     fail_on_nonunique: bool = True):
    pd = _pd()
    bed = bed.drop(columns=[c for c in ["thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"] if c in bed.columns])
    names_list = name_cols.split(",")
    if main_name not in names_list:
        raise SystemExit(f"primary gene identifier {main_name} not in {name_cols}")
    names = bed["name"].astype(str).str.split(";", expand=True)
    if names.shape[1] != len(names_list):
        raise SystemExit(f"gene names have {names.shape[1]} ';'-fields but --gene-name-annotations lists {len(names_list)}")
    names.columns = names_list
    bed = pd.concat([bed.reset_index(drop=True), names.reset_index(drop=True)], axis=1)
    bed["name"] = bed[main_name]
    bed["tss"] = get_tss(bed)
    bed = bed.drop_duplicates()
    if sizes is not None:
        bed = bed[bed["chr"].astype(str).isin(set(sizes))]
    if fail_on_nonunique and bed["name"].nunique() != len(bed):
        raise SystemExit("Gene IDs are not unique! Failing. Please ensure unique identifiers are passed to --genes")
    return bed.reset_index(drop=True)


def processed_genes_file(genes_path, sizes: "Dict[str, int]", out: Path) -> Path:
    """rule create_neighborhoods: bedtools intersect -u with chrom sizes | bedtools sort -faidx | uniq."""
    lines = []
    with _open_text(genes_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                c, s, e = f[0], int(f[1]), int(f[2])
            except (IndexError, ValueError):
                continue
            if c in sizes and s < sizes[c] and e > 0:
                lines.append((c, s, line.rstrip("\n")))
    rank = {k: i for i, k in enumerate(sizes)}
    lines.sort(key=lambda x: (rank[x[0]], x[1]))
    out_lines, prev = [], None
    for _, _, l in lines:
        if l != prev:
            out_lines.append(l)
        prev = l
    Path(out).write_text("\n".join(out_lines) + ("\n" if out_lines else ""))
    return Path(out)


def count_single_feature(df, regions_file_df, feature: str, target, directory: Path, filebase: str,
                         skip_rpkm_quantile: bool = False):
    np = _np()
    orig_n = len(df)
    feature_name = feature + "." + os.path.basename(str(target))
    outfile = Path(directory) / f"{filebase}.{feature_name}.CountReads.bedgraph"
    counts = run_count_reads(target, regions_file_df, outfile)
    total = count_total(target)
    rc = feature_name + ".readCount"
    counts = counts.rename(columns={"count": rc}).drop_duplicates()
    counts["chr"] = counts["chr"].astype(str)
    df = df.merge(counts, on=["chr", "start", "end"], how="inner")
    if len(df) != orig_n:
        raise RuntimeError("Dimension mismatch while merging counts")
    n = float(len(df))
    df[feature_name + ".RPM"] = 1e6 * df[rc] / float(total)
    if not skip_rpkm_quantile:
        df[rc + ".quantile"] = df[rc].rank() / n
        df[feature_name + ".RPM.quantile"] = df[feature_name + ".RPM"].rank() / n
        df[feature_name + ".RPKM"] = 1e3 * df[feature_name + ".RPM"] / (df["end"] - df["start"]).astype(float)
        df[feature_name + ".RPKM.quantile"] = df[feature_name + ".RPKM"].rank() / n
    return df[~df.duplicated()]


def average_features(df, feature: str, files: "List[str]", skip_rpkm_quantile: bool = False):
    n = float(len(df))
    df[feature + ".RPM"] = df[[f"{feature}.{os.path.basename(str(f))}.RPM" for f in files]].mean(axis=1)
    if not skip_rpkm_quantile:
        df[feature + ".RPM.quantile"] = df[feature + ".RPM"].rank() / n
        df[feature + ".RPKM"] = df[[f"{feature}.{os.path.basename(str(f))}.RPKM" for f in files]].mean(axis=1)
        df[feature + ".RPKM.quantile"] = df[feature + ".RPKM"].rank() / n
    return df


def count_features_for_bed(df, regions_df, features: "Dict[str, List[str]]", directory: Path, filebase: str,
                           skip_rpkm_quantile: bool = False):
    for feature, files in features.items():
        t0 = time.time()
        for f in files:
            df = count_single_feature(df, regions_df, feature, f, directory, filebase, skip_rpkm_quantile)
        df = average_features(df, feature.replace("feature_", ""), files, skip_rpkm_quantile)
        log.info("feature %s over %s in %.1fs", feature, filebase, time.time() - t0)
    return df


def load_genes(genes_path, ue_file, sizes, outdir: Path, expression_tables: "List[str]", gene_id_names: str,
               primary_id: str, cell_type, class_gene_file=None, upstream_compat: bool = False):
    pd, np = _pd(), _np()
    bed = read_gene_bed_file(genes_path)
    genes = process_gene_bed(bed, gene_id_names, primary_id, sizes)
    genes[["chr", "start", "end", "name", "score", "strand"]].to_csv(Path(outdir) / "GeneList.bed", sep="\t", index=False, header=False)
    if expression_tables:
        names = []
        for t in expression_tables:
            nm = os.path.basename(str(t)) + ".Expression"
            expr = pd.read_csv(t, sep="\t", names=[primary_id, nm], header=None)
            expr[nm] = pd.to_numeric(expr[nm], errors="coerce")
            expr = expr.groupby(primary_id).max()
            genes = genes.merge(expr, how="left", right_index=True, left_on="symbol" if "symbol" in genes.columns else primary_id)
            names.append(nm)
        genes["Expression"] = genes[names].mean(axis=1)
        genes["Expression.quantile"] = genes["Expression"].rank(method="average", na_option="top", ascending=True, pct=True)
    else:
        genes["Expression"] = np.nan
    if ue_file:
        ubiq = pd.read_csv(ue_file, sep="\t", header=0 if upstream_compat else None)
        genes["is_ue"] = genes["name"].isin(ubiq.iloc[:, 0].astype(str).tolist())
    genes["cellType"] = cell_type
    if class_gene_file:
        cg = read_bed(class_gene_file)
        genes_for_class = process_gene_bed(cg, gene_id_names, primary_id, sizes, fail_on_nonunique=False)
    else:
        genes_for_class = genes
    return genes, genes_for_class


def make_tss_region_file(genes, outdir: Path, sizes: "Dict[str, int]", tss_slop: int = 500):
    """GeneList.TSS1kb.bed: TSS +/- 500; windows leaving the chromosome are removed (pyranges genome_bounds)."""
    np = _np()
    t = genes.loc[:, ["chr", "start", "end", "name", "score", "strand"]].copy()
    t["start"] = genes["tss"].to_numpy(np.int64) - tss_slop
    t["end"] = genes["tss"].to_numpy(np.int64) + tss_slop
    size = t["chr"].map(sizes)
    t = t[(t["start"] >= 0) & (t["end"] <= size)]
    sort_bed(t, sizes).to_csv(Path(outdir) / "GeneList.TSS1kb.bed", sep="\t", index=False, header=False)
    return t.reset_index(drop=True)


def annotate_genes_with_features(genes, sizes, features, outdir: Path, default_feature: str):
    tss1kb = make_tss_region_file(genes, outdir, sizes)
    genes = count_features_for_bed(genes, genes[["chr", "start", "end"]], features, outdir, "Genes")
    tsscounts = count_features_for_bed(tss1kb, tss1kb[["chr", "start", "end"]], features, outdir, "Genes.TSS1kb")
    tsscounts = tsscounts.drop(columns=["chr", "start", "end", "score", "strand"])
    merged = genes.merge(tsscounts, on="name", suffixes=["", ".TSS1Kb"])
    acc = default_feature + ".RPKM.quantile.TSS1Kb"
    if "H3K27ac.RPKM.quantile.TSS1Kb" in merged.columns:
        val = (0.0001 + merged["H3K27ac.RPKM.quantile.TSS1Kb"]) * (0.0001 + merged[acc])
    else:
        val = 0.0001 + merged[acc]
    merged["PromoterActivityQuantile"] = val.rank(method="average", na_option="top", ascending=True, pct=True)
    write_tsv(merged, Path(outdir) / "GeneList.txt")
    return merged


def assign_enhancer_classes(enh, genes, sizes, tss_slop: int = 500):
    """promoter (overlaps TSS +/- slop, clipped) > genic (overlaps gene body) > intergenic."""
    from eqtl_enrichment_skill import overlap_join  # noqa: E402
    np = _np()
    enh = enh.copy().reset_index(drop=True)
    enh["class"] = "intergenic"
    enh["uid"] = np.arange(len(enh))
    sym = "symbol" if "symbol" in genes.columns else "name"
    tss = genes[["chr", "tss", sym]].copy()
    tss["start"] = np.maximum(genes["tss"].to_numpy(np.int64) - tss_slop, 0)
    tss["end"] = np.minimum(genes["tss"].to_numpy(np.int64) + tss_slop, genes["chr"].map(sizes).fillna(np.inf).to_numpy())
    tss = tss[["chr", "start", "end", sym]].rename(columns={sym: "symbol"})
    body = genes[["chr", "start", "end", sym]].rename(columns={sym: "symbol"})
    base = enh[["chr", "start", "end", "uid"]]

    def symbols(targets):
        j = overlap_join(base, targets.astype({"start": "int64", "end": "int64"}))
        if len(j) == 0:
            return {}
        return j.groupby("uid")["symbol"].agg(lambda x: ",".join(sorted(set(map(str, x))))).to_dict()

    genic, prom = symbols(body), symbols(tss)
    enh.loc[enh["uid"].isin(list(genic)), "class"] = "genic"
    enh.loc[enh["uid"].isin(list(prom)), "class"] = "promoter"
    enh["isPromoterElement"] = enh["class"] == "promoter"
    enh["isGenicElement"] = enh["class"] == "genic"
    enh["isIntergenicElement"] = enh["class"] == "intergenic"
    enh["promoterSymbol"] = enh["uid"].map(prom).fillna("")
    enh["genicSymbol"] = enh["uid"].map(genic).fillna("")
    enh = enh.drop(columns="uid")
    enh["name"] = enh["class"] + "|" + enh["chr"].astype(str) + ":" + enh["start"].astype(str) + "-" + enh["end"].astype(str)
    return enh


def _interp(x_ref, y_ref, x):
    from scipy import interpolate  # type: ignore
    f = interpolate.interp1d(x_ref, y_ref, kind="linear", fill_value="extrapolate")
    return f(x)


def run_qnorm(df, qnorm, qnorm_method: str = "rank", separate_promoters: bool = True):
    pd, np = _pd(), _np()
    col_dict = {"DHS.RPM": "normalized_dhs", "ATAC.RPM": "normalized_atac", "H3K27ac.RPM": "normalized_h3K27ac"}
    if qnorm is None:
        for c, n in col_dict.items():
            if c in df.columns:
                df[n] = df[c]
        return df
    ref = pd.read_csv(qnorm, sep="\t") if not hasattr(qnorm, "columns") else qnorm.copy()
    for col in [c for c in ["DHS.RPM", "ATAC.RPM", "H3K27ac.RPM"] if c in df.columns]:
        if col == "ATAC.RPM" and "ATAC.RPM" not in ref.columns:
            ref["ATAC.RPM"] = ref["DHS.RPM"]
        if not separate_promoters:
            r = ref.loc[ref["enh_class"] == "any"]
            n = len(df)
            q = df[col].rank() / n
            x = (1 - q) * n if qnorm_method == "rank" else q
            df[col_dict[col]] = np.clip(_interp(r["rank" if qnorm_method == "rank" else "quantile"], r[col], x), 0, None)
            continue
        for enh_class in ["promoter", "nonpromoter"]:
            r = ref.loc[ref["enh_class"] == enh_class]
            is_prom = df["class"].isin(["tss", "promoter"])
            idx = df.index[is_prom] if enh_class == "promoter" else df.index[~is_prom]
            qcol = col + enh_class + ".quantile"
            df.loc[idx, qcol] = df.loc[idx, col].rank() / len(idx)
            if len(idx) == 0:
                continue
            if qnorm_method == "rank":
                vals = _interp(r["rank"], r[col], (1 - df.loc[idx, qcol]) * len(idx))
            else:
                vals = _interp(r["quantile"], r[col], df.loc[idx, qcol])
            df.loc[idx, col_dict[col]] = np.clip(vals, 0, None)
    return df


def compute_activity(df, access: str):
    np = _np()
    if access not in ("DHS", "ATAC"):
        raise SystemExit("At least one of ATAC or DHS must be provided!")
    norm = "normalized_dhs" if access == "DHS" else "normalized_atac"
    raw = access + ".RPM"
    if "H3K27ac.RPM" in df.columns:
        df["activity_base"] = np.sqrt(df["normalized_h3K27ac"] * df[norm])
        df["activity_base_no_qnorm"] = np.sqrt(df["H3K27ac.RPM"] * df[raw])
    else:
        df["activity_base"] = df[norm]
        df["activity_base_no_qnorm"] = df[raw]
    return df


def determine_accessibility_feature(default: "Optional[str]", atac, dhs) -> str:
    if default:
        return default
    if atac and dhs:
        raise SystemExit("Both DHS and ATAC have been provided. Must set one to be the default accessibility feature!")
    if atac:
        return "ATAC"
    if dhs:
        return "DHS"
    raise SystemExit("At least one of ATAC or DHS must be provided!")


def get_features(h3k27ac=None, atac=None, dhs=None, supplementary=None) -> "Dict[str, List[str]]":
    pd = _pd()
    feats: "Dict[str, List[str]]" = {}
    for name, v in (("H3K27ac", h3k27ac), ("ATAC", atac), ("DHS", dhs)):
        if v:
            feats[name] = v if isinstance(v, list) else str(v).split(",")
    if supplementary:
        supp = pd.read_csv(supplementary, sep="\t")
        for _, row in supp.iterrows():
            feats[row["feature_name"]] = str(row["file"]).split(",")
    return feats


def run_neighborhoods(candidate_regions, genes_path, sizes_path, outdir: Path, dhs=None, atac=None, h3k27ac=None,
                      default_feature=None, ubiquitous=None, expression_tables=None, qnorm=None, qnorm_method="rank",
                      separate_promoters=True, cell_type=None, gene_name_annotations="symbol", primary_gene_identifier="symbol",
                      genes_for_class_assignment=None, tss_slop_for_class_assignment=500, skip_gene_counts=False,
                      skip_rpkm_quantile=False, supplementary_features=None, enhancer_class_override=None,
                      upstream_compat=False) -> "Dict[str, Any]":
    pd, np = _pd(), _np()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sizes = read_chrom_sizes(sizes_path)
    access = determine_accessibility_feature(default_feature, atac, dhs)
    features = get_features(h3k27ac, atac, dhs, supplementary_features)
    genes, genes_cls = load_genes(genes_path, ubiquitous, sizes, outdir, expression_tables or [], gene_name_annotations,
                                  primary_gene_identifier, cell_type, genes_for_class_assignment, upstream_compat)
    if not skip_gene_counts:
        genes = annotate_genes_with_features(genes, sizes, features, outdir, access)
    enh = read_bed(candidate_regions, extra_cols=[])
    enh = enh[["chr", "start", "end"]]
    enh = count_features_for_bed(enh, enh[["chr", "start", "end"]], features, outdir, "Enhancers", skip_rpkm_quantile)
    if cell_type is not None:
        enh["cellType"] = cell_type
    enh = assign_enhancer_classes(enh, genes_cls, sizes, tss_slop_for_class_assignment)
    if enhancer_class_override:
        ov = read_bed(enhancer_class_override, extra_cols=["class"])
        m = enh.merge(ov[["chr", "start", "end", "class"]], on=["chr", "start", "end"], how="left", suffixes=("", "_override"))
        enh["class"] = m["class_override"].fillna(m["class"]).to_numpy()
    enh = run_qnorm(enh, qnorm, qnorm_method, separate_promoters)
    enh = compute_activity(enh, access)
    enh[["chr", "start", "end", "name"]].to_csv(outdir / "EnhancerList.bed", sep="\t", index=False, header=False)
    write_tsv(enh, outdir / "EnhancerList.txt")
    cls = enh["class"].value_counts().to_dict()
    return {"enhancer_list": str(outdir / "EnhancerList.txt"), "gene_list": str(outdir / "GeneList.txt"),
            "n_enhancers": int(len(enh)), "n_genes": int(len(genes)), "classes": {k: int(v) for k, v in cls.items()},
            "accessibility_feature": access, "features": {k: [os.path.basename(str(x)) for x in v] for k, v in features.items()},
            "qnorm": str(qnorm) if qnorm is not None and not hasattr(qnorm, "columns") else ("in-memory" if qnorm is not None else None)}


# ---------------------------------------------------------------------------
# Hi-C (hic.py / predictor.py)
# ---------------------------------------------------------------------------

def get_powerlaw_at_distance(distances, gamma: float, scale: float, min_distance: "Optional[float]" = 5000):
    np = _np()
    if not (gamma > 0 and scale > 0):
        raise ValueError("power-law gamma and scale must be positive")
    d = np.clip(np.asarray(distances, dtype=float), min_distance if min_distance is not None else -np.inf, np.inf)
    return np.exp(scale + -1 * gamma * np.log(d + 1))


def _exists_nonempty(p: Path) -> bool:
    if not p.exists():
        return False
    return p.stat().st_size > (100 if str(p).endswith("gz") else 0)


def get_hic_file(chromosome: str, hic_dir, allow_vc: bool = True, hic_type: str = "juicebox"):
    d = Path(hic_dir) / chromosome
    if hic_type == "juicebox":
        for ft in ("KR", "INTERSCALE"):
            for ext in (".gz", ""):
                f = d / f"{chromosome}.{ft}observed{ext}"
                if _exists_nonempty(f):
                    return f, d / f"{chromosome}.{ft}norm{ext}", False
        if allow_vc:
            for ext in (".gz", ""):
                f = d / f"{chromosome}.VCobserved{ext}"
                if _exists_nonempty(f):
                    log.info("using VC-normalised Hi-C %s", f)
                    return f, d / f"{chromosome}.VCnorm{ext}", True
        raise FileNotFoundError(f"Could not find KR, INTERSCALE or VC normalized hic files for {chromosome} in {hic_dir}")
    if hic_type == "bedpe":
        return d / f"{chromosome}.bedpe.gz", None, None
    if hic_type == "avg":
        return d / f"{chromosome}.bed.gz", None, None
    raise ValueError(f"unknown hic_type {hic_type}")


def hic_to_sparse(filename, norm_file, resolution: int, hic_is_doubly_stochastic: bool = False):
    pd, np = _pd(), _np()
    import scipy.sparse as ssp  # type: ignore
    hic = pd.read_csv(filename, sep="\t", names=["bin1", "bin2", "hic_contact"], header=None)
    if not np.all(hic["bin1"].to_numpy() <= hic["bin2"].to_numpy()):
        raise ValueError("juicebox dump must be upper triangular (bin1 <= bin2)")
    size = len(pd.read_csv(norm_file, header=None))
    row = np.floor(hic["bin1"].to_numpy() / resolution).astype(int)
    col = np.floor(hic["bin2"].to_numpy() / resolution).astype(int)
    dat = hic["hic_contact"].to_numpy(dtype=float)
    if not hic_is_doubly_stochastic:
        off = row != col
        row, col, dat = np.hstack([row, col[off]]), np.hstack([col, row[off]]), np.hstack([dat, dat[off]])
    return ssp.csr_matrix((dat, (row, col)), (size, size))


def process_hic(mat, hic_is_vc: bool, resolution: int, tss_hic_contribution: float, window: float, min_window: float = 0,
                hic_is_doubly_stochastic: bool = False, apply_diagonal_bin_correction: bool = True):
    pd, np = _pd(), _np()
    import scipy.sparse as ssp  # type: ignore
    mat = mat.tocsr()
    if not hic_is_doubly_stochastic and not hic_is_vc:
        mat.data = np.nan_to_num(mat.data, copy=False)      # upstream zeroes NaNs in place
        sums = np.asarray(mat.sum(axis=0)).ravel()
        sums = sums[~np.isnan(sums)]
        mean_sum = float(np.mean(sums[sums > 0]))
        if abs(mean_sum - 1) >= 0.001:
            mat = mat.multiply(1 / mean_sum).tocsr()
    if apply_diagonal_bin_correction:
        n = mat.shape[0]
        diag = mat.diagonal()
        nz = np.flatnonzero(diag != 0)
        if nz.size:
            left = np.where(nz > 0, np.asarray(mat[nz, np.maximum(nz - 1, 0)]).ravel(), -np.inf)
            right = np.where(nz < n - 1, np.asarray(mat[nz, np.minimum(nz + 1, n - 1)]).ravel(), -np.inf)
            new = np.maximum(left, right) * tss_hic_contribution / 100
            mat = mat.tolil()
            mat[nz, nz] = new
            mat = mat.tocsr()
    if not hic_is_vc:
        mat = ssp.triu(mat)
    else:
        colsum = np.asarray(mat.sum(axis=0)).ravel()
        colsum[colsum == 0] = 1
        mat = ssp.diags(1.0 / colsum) @ mat
    coo = mat.tocoo()
    df = pd.DataFrame({"bin1": coo.row, "bin2": coo.col, "hic_contact": coo.data})
    gap = (df["bin1"] - df["bin2"]).abs()
    return df[(gap <= window / resolution) & (gap >= min_window / resolution)].reset_index(drop=True)


def load_hic_juicebox(hic_file, hic_norm_file, hic_is_vc, hic_resolution, tss_hic_contribution=100.0, window=5e6,
                      min_window=0, hic_is_doubly_stochastic=False, apply_diagonal_bin_correction=True):
    m = hic_to_sparse(hic_file, hic_norm_file, hic_resolution)
    return process_hic(m, bool(hic_is_vc), hic_resolution, tss_hic_contribution, window, min_window,
                       hic_is_doubly_stochastic, apply_diagonal_bin_correction)


def load_hic_bedpe(hic_file):
    pd = _pd()
    return pd.read_csv(hic_file, sep="\t", names=["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "hic_contact"])


def load_hic_avg(hic_file, hic_resolution: int):
    pd, np = _pd(), _np()
    h = pd.read_csv(hic_file, sep="\t", names=["x1", "x2", "hic_contact"], usecols=[0, 1, 2],
                    dtype={"x1": np.int64, "x2": np.int64, "hic_contact": np.float64})
    h["x1"] = np.floor(h["x1"] / hic_resolution).astype(int)
    h["x2"] = np.floor(h["x2"] / hic_resolution).astype(int)
    return h.rename(columns={"x1": "bin1", "x2": "bin2"})


def hic_records_to_contacts(records, hic_resolution: int):
    """predictor.add_hic_from_hic_file on one full-chromosome query: bins = floor(pos / res); diagonal =
    max of neighbouring bins within 5 kb (from the raw records); counts / mean row sum."""
    pd, np = _pd(), _np()
    df = pd.DataFrame(records, columns=["binX", "binY", "counts"])
    raw = df.copy()
    bx, by, v = raw["binX"].to_numpy(), raw["binY"].to_numpy(), raw["counts"].to_numpy(float)
    s = pd.Series(v).groupby(bx).sum().add(pd.Series(np.where(bx != by, v, 0.0)).groupby(by).sum(), fill_value=0)
    row_mean = float(s.mean())
    df["binX"] = np.floor(df["binX"] / hic_resolution).astype(int)
    df["binY"] = np.floor(df["binY"] / hic_resolution).astype(int)
    df = df.groupby(["binX", "binY"], as_index=False)["counts"].max()
    look = dict(zip(zip(df["binX"], df["binY"]), df["counts"]))
    search = math.ceil(5000 / hic_resolution) if hic_resolution < 5000 else 1
    diag = df["binX"] == df["binY"]
    new = []
    for b in df.loc[diag, "binX"]:
        mx = 0.0
        for i in range(1, search + 1):
            for key in ((b - i, b), (b, b + i)):
                if key in look:
                    mx = max(mx, look[key])
        new.append(mx)
    df.loc[diag, "counts"] = new
    df["counts"] = df["counts"] / row_mean
    return df, row_mean


def add_hic_from_hic_file(pred, hic_file, chromosome: str, hic_resolution: int):
    np = _np()
    try:
        import hicstraw  # type: ignore
    except ImportError as e:
        raise SystemExit("hic_type 'hic' streams a .hic file and needs hicstraw (pip install hic-straw); "
                         "alternatively dump it with `abc-pipeline juicebox-dump` and use --hic-type juicebox") from e
    hic = hicstraw.HiCFile(str(hic_file))
    names = [c.name for c in hic.getChromosomes()]
    chrom = chromosome if (len(names) > 1 and names[1].startswith("chr")) else chromosome[3:]
    size = next(c.length for c in hic.getChromosomes() if c.name == chrom)
    mzd = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", hic_resolution)
    recs = [[r.binX, r.binY, r.counts] for r in mzd.getRecords(0, size, 0, size)]
    return merge_hic_records(pred, recs, hic_resolution)


def merge_hic_records(pred, records, hic_resolution: int):
    np = _np()
    df, _ = hic_records_to_contacts(records, hic_resolution)
    eb = np.floor(pred["enh_midpoint"] / hic_resolution).astype(int)
    tb = np.floor(pred["TargetGeneTSS"] / hic_resolution).astype(int)
    pred = pred.assign(binX=np.minimum(eb, tb), binY=np.maximum(eb, tb))
    pred = pred.merge(df, how="left", on=["binX", "binY"])
    pred = pred.drop(columns=["binX", "binY", "enh_idx", "gene_idx", "enh_midpoint"], errors="ignore")
    return pred.rename(columns={"counts": "hic_contact"})


def add_hic_from_directory(chromosome: str, enh, genes, pred, hic_dir, opts):
    pd, np = _pd(), _np()
    hic_file, hic_norm, hic_is_vc = get_hic_file(chromosome, hic_dir, hic_type=opts.hic_type)
    if opts.hic_type == "bedpe":
        from eqtl_enrichment_skill import overlap_join  # noqa: E402
        H = load_hic_bedpe(hic_file)
        H["hic_idx"] = np.arange(len(H))
        e = pd.DataFrame({"chr": enh["chr"].to_numpy(), "start": np.floor(enh["enh_midpoint"]).astype(np.int64),
                          "enh_idx": enh["enh_idx"].to_numpy()})
        e["end"] = e["start"] + 1
        g = pd.DataFrame({"chr": genes["chr"].to_numpy(), "start": genes["TargetGeneTSS"].astype(np.int64).to_numpy(),
                          "gene_idx": genes["gene_idx"].to_numpy()})
        g["end"] = g["start"] + 1
        h1 = H.rename(columns={"chr1": "chr", "x1": "start", "x2": "end"})[["chr", "start", "end", "hic_idx", "hic_contact"]]
        h2 = H.rename(columns={"chr2": "chr", "y1": "start", "y2": "end"})[["chr", "start", "end", "hic_idx", "hic_contact"]]
        eh1 = overlap_join(e, h1)[["enh_idx", "hic_idx", "hic_contact"]]
        gh2 = overlap_join(g, h2)[["gene_idx", "hic_idx"]]
        eh2 = overlap_join(e, h2)[["enh_idx", "hic_idx", "hic_contact"]]
        gh1 = overlap_join(g, h1)[["gene_idx", "hic_idx"]]
        ovl = pd.concat([eh1.merge(gh2, on="hic_idx"), eh2.merge(gh1, on="hic_idx")]).drop_duplicates()
        ovl = ovl.astype({"enh_idx": "int64", "gene_idx": "int64"})
        pred = pred.merge(ovl, on=["enh_idx", "gene_idx"], how="left")
    else:
        if opts.hic_type == "juicebox":
            H = load_hic_juicebox(hic_file, hic_norm, hic_is_vc, opts.hic_resolution, opts.tss_hic_contribution,
                                  opts.window, 0, getattr(opts, "hic_is_doubly_stochastic", False))
        else:
            H = load_hic_avg(hic_file, opts.hic_resolution)
        pred["enh_bin"] = np.floor(pred["enh_midpoint"] / opts.hic_resolution).astype(int)
        pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / opts.hic_resolution).astype(int)
        if not hic_is_vc:
            pred["bin1"] = np.minimum(pred["enh_bin"], pred["tss_bin"])
            pred["bin2"] = np.maximum(pred["enh_bin"], pred["tss_bin"])
            pred = pred.merge(H, how="left", on=["bin1", "bin2"])
        else:
            pred = pred.merge(H, how="left", left_on=["tss_bin", "enh_bin"], right_on=["bin1", "bin2"])
    return pred.drop(columns=["x1", "x2", "y1", "y2", "bin1", "bin2", "enh_idx", "gene_idx", "hic_idx", "enh_midpoint",
                              "tss_bin", "enh_bin"], errors="ignore")


# ---------------------------------------------------------------------------
# Predictions (predict.py / predictor.py / getVariantOverlap.py)
# ---------------------------------------------------------------------------

def determine_expressed_genes(genes, expression_cutoff: float, activity_quantile_cutoff: float):
    np = _np()
    e = genes["Expression"].astype(float)
    genes["isExpressed"] = np.logical_or(e >= expression_cutoff,
                                         np.logical_and(np.isnan(e), genes["PromoterActivityQuantile"] >= activity_quantile_cutoff))
    return genes


def make_pred_table(enh, genes, window: int, sizes: "Dict[str, int]"):
    """pyranges join of elements with [TSS - window, TSS + window] (clipped), then |midpoint - TSS| < window."""
    pd, np = _pd(), _np()
    enh = enh.copy()
    genes = genes.copy()
    enh["enh_midpoint"] = (enh["start"] + enh["end"]) / 2
    enh["enh_idx"] = enh.index
    genes["gene_idx"] = genes.index
    if len(enh) == 0 or len(genes) == 0:
        cols = list(enh.columns) + [c for c in genes.columns if c != "chr"] + ["distance"]
        return pd.DataFrame(columns=cols), enh, genes
    order = np.argsort(enh["start"].to_numpy(), kind="mergesort")
    es = enh["start"].to_numpy(np.int64)[order]
    ee = enh["end"].to_numpy(np.int64)[order]
    maxlen = int((ee - es).max())
    size = sizes.get(str(genes["chr"].iloc[0]), np.inf)
    tss = genes["TargetGeneTSS"].to_numpy(np.int64)
    lo = np.maximum(tss - window, 0)
    hi = np.minimum(tss + window, size)
    i0 = np.searchsorted(es, lo - maxlen, side="left")
    i1 = np.searchsorted(es, hi, side="left")
    n = np.maximum(i1 - i0, 0)
    gi = np.repeat(np.arange(len(genes)), n)
    offs = np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n)
    ei = np.repeat(i0, n) + offs
    keep = ee[ei] > lo[gi]
    gi, ei = gi[keep], order[ei[keep]]
    srt = np.lexsort((gi, ei))
    gi, ei = gi[srt], ei[srt]
    left = enh.iloc[ei].reset_index(drop=True)
    right = genes.iloc[gi].drop(columns=["chr"]).reset_index(drop=True)
    pred = pd.concat([left, right], axis=1)
    pred["distance"] = (pred["enh_midpoint"] - pred["TargetGeneTSS"]).abs()
    pred = pred.loc[pred["distance"] < window].reset_index(drop=True)
    return pred, enh, genes


def annotate_predictions(pred, tss_slop: int = 500):
    np = _np()
    pred["isSelfPromoter"] = np.logical_and.reduce((pred["class"] == "promoter",
                                                    pred["start"] - tss_slop < pred["TargetGeneTSS"],
                                                    pred["end"] + tss_slop > pred["TargetGeneTSS"]))
    return pred


def add_powerlaw_to_predictions(pred, opts):
    pred["powerlaw_contact"] = get_powerlaw_at_distance(pred["distance"].to_numpy(), opts.hic_gamma, opts.hic_scale)
    scale_ref = -4.80 + 11.63 * opts.hic_gamma_reference
    pred["powerlaw_contact_reference"] = get_powerlaw_at_distance(pred["distance"].to_numpy(), opts.hic_gamma_reference, scale_ref)
    return pred


def qc_hic(pred, gamma: float, scale: float, resolution, threshold: float = 0.01):
    summ = pred.loc[pred["isSelfPromoter"]].groupby("TargetGene").agg({"hic_contact": "sum"})
    bad = summ.index[summ["hic_contact"] < threshold]
    aff = pred["TargetGene"].isin(bad)
    if aff.any():
        pred.loc[aff, "hic_contact"] = get_powerlaw_at_distance(pred.loc[aff, "distance"].to_numpy(), gamma, scale,
                                                                min_distance=resolution)
    return pred


def scale_hic_with_powerlaw(pred, scale: bool):
    if not scale:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"]
    else:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"] * (pred["powerlaw_contact_reference"] / pred["powerlaw_contact"])
    return pred


def add_hic_pseudocount(pred, opts):
    np = _np()
    pc = float(get_powerlaw_at_distance(opts.hic_pseudocount_distance, opts.hic_gamma, opts.hic_scale, opts.hic_pseudocount_distance))
    pred["hic_pseudocount"] = np.minimum(pred["powerlaw_contact"].to_numpy(), pc)
    pred["hic_contact_pl_scaled_adj"] = pred["hic_contact_pl_scaled"] + pred["hic_pseudocount"]
    return pred


def compute_score(pred, product_terms, prefix: str, adjust_self_promoters: bool = True):
    np = _np()
    pred[prefix + ".Score.Numerator"] = np.column_stack(product_terms).prod(axis=1)
    denom = pred.groupby(["TargetGene", "TargetGeneTSS"])[prefix + ".Score.Numerator"].transform("sum")
    pred[prefix + ".Score"] = pred[prefix + ".Score.Numerator"] / denom
    if adjust_self_promoters:
        pred.loc[pred["isSelfPromoter"].astype(bool), prefix + ".Score"] = 1
    return pred


def make_predictions(chromosome: str, enh, genes, opts, sizes, hic_records=None):
    pred, enh2, genes2 = make_pred_table(enh, genes, opts.window, sizes)
    pred = annotate_predictions(pred, opts.tss_slop)
    pred = add_powerlaw_to_predictions(pred, opts)
    if opts.hic_file or hic_records is not None:
        if hic_records is not None:
            pred = merge_hic_records(pred, hic_records, opts.hic_resolution)
        elif opts.hic_type == "hic":
            pred = add_hic_from_hic_file(pred, opts.hic_file, chromosome, opts.hic_resolution)
        else:
            pred = add_hic_from_directory(chromosome, enh2, genes2, pred, opts.hic_file, opts)
        pred = qc_hic(pred, opts.hic_gamma, opts.hic_scale, opts.hic_resolution)
        pred["hic_contact"] = pred["hic_contact"].fillna(0)
        pred = scale_hic_with_powerlaw(pred, opts.scale_hic_using_powerlaw)
        pred = add_hic_pseudocount(pred, opts)
        pred = compute_score(pred, [pred["activity_base_enh"], pred["hic_contact_pl_scaled_adj"]], "ABC")
    else:
        pred = compute_score(pred, [pred["activity_base_enh"], pred["powerlaw_contact"]], "ABC")
    pred = compute_score(pred, [pred["activity_base_enh"], pred["powerlaw_contact"]], "powerlaw")
    return pred


def shrink_regions(df, bp: int, sizes: "Optional[Dict[str, int]]" = None):
    """bedtools slop -b -bp: trim both ends; intervals shorter than 2 bp collapse to their midpoint."""
    np = _np()
    df = df.copy()
    s = df["start"].to_numpy(np.int64) + bp
    e = df["end"].to_numpy(np.int64) - bp
    bad = s > e
    mid = (df["start"].to_numpy(np.int64) + df["end"].to_numpy(np.int64)) // 2
    s[bad] = mid[bad]; e[bad] = mid[bad]
    df["start"], df["end"] = s, e
    return df


def variant_overlap(all_putative, outdir: Path, score_column: str = "ABC.Score"):
    """getVariantOverlap.test_variant_overlap."""
    np = _np()
    sc = all_putative[score_column]
    prom = all_putative["class"] == "promoter"
    v = all_putative[((sc > 0.015) & ~prom) | (prom & (sc > 0.1))]
    v = v.dropna(subset=[score_column])
    v = v.loc[v["distance"] <= 2_000_000]
    path = Path(outdir) / "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz"
    write_tsv(shrink_regions(v, 150), path, na_rep="NaN")
    return path, int(len(v))


def predict(enhancers_path, genes_path, sizes_path, outdir: Path, opts, hic_records_by_chrom=None) -> "Dict[str, Any]":
    pd, np = _pd(), _np()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if opts.hic_file and opts.hic_type in ("hic", "juicebox") and not opts.hic_resolution:
        raise SystemExit("HiC resolution must be provided if hic_type is hic or juicebox")
    if not opts.hic_file and hic_records_by_chrom is None:
        print("WARNING: Hi-C not provided. Model will only compute ABC score using powerlaw!")
    write_params(outdir, "parameters.predict.txt", {k: v for k, v in vars(opts).items() if not k.startswith("_") and k != "func"})
    acc = opts.accessibility_feature
    if acc not in ("ATAC", "DHS"):
        raise SystemExit("The feature has to be either ATAC or DHS!")
    genes = pd.read_csv(genes_path, sep="\t")
    genes["chr"] = genes["chr"].astype(str)
    genes = determine_expressed_genes(genes, opts.expression_cutoff, opts.promoter_activity_quantile_cutoff)
    enh_full = pd.read_csv(enhancers_path, sep="\t")
    enh_full["chr"] = enh_full["chr"].astype(str)
    norm = f"normalized_{acc.lower()}"
    gcols = ["chr", "symbol", "tss", "Expression", "PromoterActivityQuantile", "isExpressed", "Ensembl_ID", f"{acc}.RPKM.quantile.TSS1Kb"]
    gnew = ["chr", "TargetGene", "TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile",
            "TargetGeneIsExpressed", "TargetGeneEnsembl_ID", f"{norm}_prom"]
    if "H3K27ac.RPKM.quantile.TSS1Kb" in genes.columns:
        gcols.append("H3K27ac.RPKM.quantile.TSS1Kb"); gnew.append("normalized_h3k27ac_prom")
    if "symbol" not in genes.columns:
        genes["symbol"] = genes["name"]
    genes = genes.loc[:, gcols]
    genes.columns = gnew
    enh = enh_full.loc[:, ["chr", "start", "end", "name", "class", "activity_base"]].copy()
    enh["activity_base_enh"] = enh_full["activity_base"]
    enh["activity_base_squared_enh"] = enh["activity_base_enh"] ** 2
    enh[f"{norm}_enh"] = enh_full[norm]
    if "normalized_h3K27ac" in enh_full.columns:
        enh["normalized_h3k27ac_enh"] = enh_full["normalized_h3K27ac"]
    sizes = read_chrom_sizes(sizes_path)
    if opts.chromosomes == "all":
        chroms = set(genes["chr"]) & set(enh["chr"])
        if not opts.include_chrY:
            chroms.discard("chrY")
        chroms = sorted(chroms)
    else:
        chroms = opts.chromosomes.split(",")
    parts = []
    for c in chroms:
        t0 = time.time()
        recs = hic_records_by_chrom.get(c) if hic_records_by_chrom else None
        this = make_predictions(c, enh.loc[enh["chr"] == c], genes.loc[genes["chr"] == c], opts, sizes,
                                hic_records=recs if hic_records_by_chrom else None)
        parts.append(this)
        log.info("chromosome %s: %d pairs in %.1fs", c, len(this), time.time() - t0)
    allp = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    allp["CellType"] = opts.cellType
    has_hic = bool(opts.hic_file) or hic_records_by_chrom is not None
    if has_hic:
        allp["hic_contact_squared"] = allp["hic_contact"] ** 2
    expr = allp["TargetGeneIsExpressed"].astype(bool) if len(allp) else []
    f1 = outdir / "EnhancerPredictionsAllPutative.tsv.gz"
    f2 = outdir / "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"
    write_tsv(allp.loc[expr] if len(allp) else allp, f1, na_rep="NaN")
    write_tsv(allp.loc[~expr] if len(allp) else allp, f2, na_rep="NaN")
    vo, nvo = variant_overlap(allp, outdir, opts.score_column)
    return {"all_putative": str(f1), "all_putative_nonexpressed": str(f2), "variant_overlap": str(vo),
            "n_pairs": int(len(allp)), "n_pairs_expressed": int(np.sum(expr)) if len(allp) else 0,
            "n_variant_overlap": nvo, "chromosomes": chroms, "contact": ("Hi-C " + str(opts.hic_type)) if has_hic else "powerlaw",
            "n_genes": int(allp["TargetGene"].nunique()) if len(allp) else 0}


def predict_opts(**kw) -> argparse.Namespace:
    base = dict(score_column=DEFAULTS["score_column"], accessibility_feature="DHS", cellType=None, hic_file=None,
                hic_resolution=None, hic_pseudocount_distance=DEFAULTS["hic_pseudocount_distance"], hic_type="hic",
                hic_is_doubly_stochastic=False, scale_hic_using_powerlaw=True, hic_gamma=DEFAULTS["hic_gamma"],
                hic_scale=DEFAULTS["hic_scale"], hic_gamma_reference=DEFAULTS["hic_gamma_reference"],
                expression_cutoff=DEFAULTS["expression_cutoff"],
                promoter_activity_quantile_cutoff=DEFAULTS["promoter_activity_quantile_cutoff"], window=DEFAULTS["window"],
                tss_hic_contribution=DEFAULTS["tss_hic_contribution"], tss_slop=DEFAULTS["tss_slop"], chromosomes="all",
                include_chrY=False)
    base.update(kw)
    return argparse.Namespace(**base)


# ---------------------------------------------------------------------------
# Thresholds and filtering (rules/utils.smk, filter_predictions.py, predictor.make_gene_prediction_stats)
# ---------------------------------------------------------------------------

def load_thresholds(path=None):
    pd = _pd()
    if path:
        t = pd.read_csv(path, sep="\t")
        t["has_h3k27ac"] = t["has_h3k27ac"].map(lambda v: _truthy(v) if not isinstance(v, bool) else v)
        return t
    return pd.DataFrame(ABC_THRESHOLDS, columns=["accessibility", "has_h3k27ac", "hic_type", "threshold"])


def determine_threshold(accessibility: str, has_h3k27ac: bool, hic_type: "Optional[str]", explicit: "Optional[float]" = None,
                        thresholds_path=None, upstream_compat: bool = False) -> "Tuple[float, str]":
    if explicit is not None:
        return float(explicit), "explicit"
    if _blank(hic_type):
        ht = "powerlaw"
    elif hic_type == "avg":
        ht = "avg_hic" if upstream_compat else "avg"
    elif hic_type == "hic":
        ht = "intact_hic"
    else:
        ht = str(hic_type)
    t = load_thresholds(thresholds_path)
    m = t[(t["accessibility"] == accessibility) & (t["has_h3k27ac"] == bool(has_h3k27ac)) & (t["hic_type"] == ht)]
    if len(m) == 0:
        return DEFAULTS["default_threshold"], f"default (no abc_thresholds row for {accessibility}/{bool(has_h3k27ac)}/{ht})"
    return float(m.iloc[0]["threshold"]), f"abc_thresholds.tsv ({accessibility}, H3K27ac={bool(has_h3k27ac)}, {ht})"


def filtered_file_format(threshold: float, include_self_promoter: bool, only_expressed_genes: bool) -> str:
    flags = []
    if include_self_promoter:
        flags.append("self_promoter")
    if only_expressed_genes:
        flags.append("only_expr_genes")
    sep = "_" if flags else ""
    return f"threshold{threshold}{sep}{'__'.join(flags)}"


def write_connections_bedpe(pred, outfile: Path, score_column: str):
    pd = _pd()
    pred = pred.drop_duplicates()
    t = pd.DataFrame({
        "chr1": pred["chr"], "x1": pred["start"], "x2": pred["end"], "chr2": pred["chr"],
        "y1": pred["TargetGeneTSS"], "y2": pred["TargetGeneTSS"],
        "name": pred["TargetGene"].astype(str) + "|" + pred["chr"].astype(str) + ":" + pred["start"].astype(str) + "-" + pred["end"].astype(str),
        "score": pred[score_column], "strand1": ".", "strand2": "."})
    t.to_csv(outfile, header=False, index=False, sep="\t", compression="gzip" if str(outfile).endswith(".gz") else None)
    return outfile


def make_gene_prediction_stats(pred, score_column: str, threshold: float, output_file: Path):
    np = _np()
    keys = ["chr", "TargetGene", "TargetGeneTSS"]
    s1 = pred.groupby(keys).agg(geneIsExpressed=("TargetGeneIsExpressed", "first"),
                                geneFailed=(score_column, lambda x: bool(np.all(np.isnan(x.to_numpy(float))))),
                                nEnhancersConsidered=("name", "count"))
    s2 = pred.loc[pred["class"] != "promoter"].groupby(keys).agg(nDistalEnhancersPredicted=(score_column, lambda x: int((x > threshold).sum())))
    s1 = s1.merge(s2, left_index=True, right_index=True)
    s1.to_csv(output_file, sep="\t", index=True)
    return s1


def remove_promoters(df, keep_self_promoters: bool):
    if keep_self_promoters:
        return df[(df["class"] != "promoter") | df["isSelfPromoter"].astype(bool)]
    return df[df["class"] != "promoter"]


def filter_predictions(pred_file, pred_nonexpressed_file, outdir: Path, score_column: str, threshold: float,
                       include_self_promoter: bool = True, only_expressed_genes: bool = False) -> "Dict[str, Any]":
    pd = _pd()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    allp = pd.read_csv(pred_file, sep="\t")
    if not only_expressed_genes and pred_nonexpressed_file and Path(pred_nonexpressed_file).exists():
        allp = pd.concat([allp, pd.read_csv(pred_nonexpressed_file, sep="\t")], ignore_index=True)
    filt = allp[allp[score_column] > threshold]
    filt = remove_promoters(filt, include_self_promoter)
    fmt = filtered_file_format(threshold, include_self_promoter, only_expressed_genes)
    full = outdir / f"EnhancerPredictionsFull_{fmt}.tsv"
    slim = outdir / f"EnhancerPredictions_{fmt}.tsv"
    bedpe = outdir / f"EnhancerPredictionsFull_{fmt}.bedpe.gz"
    stats = outdir / f"GenePredictionStats_{fmt}.tsv"
    write_tsv(filt, full)
    write_tsv(filt[["chr", "start", "end", "name", "TargetGene", "TargetGeneTSS", "CellType", score_column]], slim)
    write_connections_bedpe(filt, bedpe, score_column)
    emit("Wrote", bedpe)
    st = make_gene_prediction_stats(filt, score_column, threshold, stats)
    emit("TSV", stats)
    distal = filt[filt["class"] != "promoter"]
    return {"threshold": threshold, "file_format": fmt, "full": str(full), "slim": str(slim), "bedpe": str(bedpe),
            "gene_stats": str(stats), "n_input_pairs": int(len(allp)), "n_kept": int(len(filt)),
            "n_distal_kept": int(len(distal)), "n_self_promoters_kept": int(filt["isSelfPromoter"].astype(bool).sum()),
            "n_genes_with_distal": int(distal["TargetGene"].nunique()), "n_gene_stats_rows": int(len(st))}


# ---------------------------------------------------------------------------
# QC (grabMetrics.py / metrics.py)
# ---------------------------------------------------------------------------

def _stats3(s):
    return float(s.mean()), float(s.median()), float(s.std())


def _hist(plt, series, title, xlabel, path, pdf, stat="count", color=BLUE):
    import numpy as np  # type: ignore
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    vals = np.asarray(series, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size:
        ax.hist(vals, bins=50, color=color, alpha=0.85, density=(stat == "density"))
        mean, med = float(np.mean(vals)), float(np.median(vals))
        ax.axvline(mean, color="#e34948", lw=1.2, label=f"Mean={round(mean, 3)}")
        ax.axvline(med, color=GREEN, lw=1.2, label=f"Median={round(med, 3)}")
        ax.legend(fontsize=7, frameon=False)
    _style(ax, title, xlabel, stat.capitalize())
    return _save(fig, path, pdf)


def generate_qc(preds_file, candidate_regions, neighborhood_dir, sizes_path, outdir: Path, qc_name: str,
                hic_gamma: float, hic_scale: float, narrowpeak=None, plots: bool = True) -> "Dict[str, Any]":
    pd, np = _pd(), _np()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chrom_order = list(read_chrom_sizes(sizes_path))
    pred = pd.read_csv(preds_file, sep="\t")
    plt = _plt() if plots else None
    pdf = None
    figs = []
    if plt is not None:
        from matplotlib.backends.backend_pdf import PdfPages  # type: ignore
        pdf = PdfPages(str(outdir / f"QCPlots_{qc_name}.pdf"))
    m: "Dict[str, Any]" = {}
    epg = pred.groupby("TargetGene").size()
    epg.to_csv(outdir / "EnhancerPerGene.tsv", sep="\t")
    gpe = pred[["chr", "start", "end"]].groupby(["chr", "start", "end"]).size()
    gpe.to_csv(outdir / "GenesPerEnhancer.tsv", sep="\t")
    per_chr = pred.groupby("chr").size()
    per_chr = per_chr.reindex([c for c in chrom_order if c in per_chr.index])
    per_chr.to_csv(outdir / "EnhancerGenePairsPerChrom.txt", sep="\t")
    for p in ("EnhancerPerGene.tsv", "GenesPerEnhancer.tsv", "EnhancerGenePairsPerChrom.txt"):
        emit("TSV", outdir / p)
    gmean, gmed, gstd = _stats3(epg)
    emean, emed, estd = _stats3(gpe)
    cmean, cmed, _ = _stats3(per_chr)
    dist = pred["distance"]
    dist = dist[dist > 0]
    m["MedianEnhPerGene"] = gmed; m["StdEnhPerGene"] = gstd
    m["MedianGenePerEnh"] = emed; m["StdGenePerEnh"] = estd
    m["MeanEnhPerGene"] = gmean; m["MeanGenePerEnh"] = emean
    m["MedianEGDist"] = float(np.median(dist)) if len(dist) else float("nan")
    m["MeanEGDist"] = float(np.mean(dist)) if len(dist) else float("nan")
    m["StdEGDist"] = float(np.std(dist)) if len(dist) else float("nan")
    m["NumEnhancersPerChrom"] = cmed; m["MeanEnhancersPerChrom"] = cmean
    m["NumEnhancers"] = int(len(pred[["chr", "start", "end"]].drop_duplicates()))
    m["EG10th"] = float(np.percentile(dist, 10)) if len(dist) else float("nan")
    m["EG90th"] = float(np.percentile(dist, 90)) if len(dist) else float("nan")
    cand = pd.read_csv(candidate_regions, sep="\t", header=None, comment="#")
    peaks = pd.read_csv(narrowpeak, sep="\t", header=None, comment="#") if narrowpeak else cand
    pw = peaks[2] - peaks[1]
    cw = cand[2] - cand[1]
    m["NumPeaks"] = int(len(pw)); m["MedWidth"] = float(pw.median()); m["MeanWidth"] = float(pw.mean()); m["StdWidth"] = float(pw.std())
    m["NumCandidate"] = int(len(cw)); m["MedWidthCandidate"] = float(cw.median())
    m["MeanWidthCandidate"] = float(cw.mean()); m["StdWidthCandidate"] = float(cw.std())
    if plt is not None:
        figs.append(_hist(plt, epg, "Number of Enhancers Per Gene", "Enhancers per Gene", outdir / "qc_enhancers_per_gene.png", pdf))
        figs.append(_hist(plt, gpe, "Number Of Genes Per Enhancer", "Genes Per Enhancer", outdir / "qc_genes_per_enhancer.png", pdf))
        fig, ax = plt.subplots(figsize=(5.2, max(2.0, 0.25 * len(per_chr) + 1)))
        ax.barh(list(per_chr.index.astype(str)), per_chr.to_numpy(), color=BLUE)
        _style(ax, "Enhancer-Gene Pairs Per Chromosome", "Number of E-G pairs", "Chromosome")
        figs.append(_save(fig, outdir / "qc_pairs_per_chrom.png", pdf))
        figs.append(_hist(plt, np.log10(dist.to_numpy(float)), "Enhancer-Gene Distance", "Log10 Distance", outdir / "qc_eg_distance.png", pdf, stat="density"))
        figs.append(_hist(plt, cw, "WidthOfCandidateRegions", "Width", outdir / "qc_candidate_width.png", pdf, color=ORANGE))
        if "hic_contact" in pred.columns:
            h = pred.loc[(pred["distance"] > 10000) & (pred["distance"] < 1_000_000)]
            h = h.sample(min(10000, len(h)), random_state=0) if len(h) else h
            h = h[(h["distance"] > 0) & (h["hic_contact"] > 0)]
            if len(h):
                fig, ax = plt.subplots(figsize=(5.2, 3.6))
                lx, ly = np.log(h["distance"].to_numpy(float)), np.log(h["hic_contact"].to_numpy(float))
                ax.scatter(lx, ly, s=5, color=BLUE, alpha=0.4, lw=0)
                xs = np.sort(lx)
                ax.plot(xs, hic_scale - hic_gamma * xs, color="#e34948", lw=1.4, label="Fitted Powerlaw Fit")
                ax.legend(fontsize=7, frameon=False)
                _style(ax, "E-G Pair HiC Powerlaw Fit", "natural log (distance)", "natural log (hic_contact)")
                figs.append(_save(fig, outdir / "qc_hic_powerlaw.png", pdf))
    if pdf is not None:
        pdf.close()
        emit("Wrote", outdir / f"QCPlots_{qc_name}.pdf")
    if neighborhood_dir:
        nd = Path(neighborhood_dir)
        for feat in ("DHS", "H3K27ac", "ATAC"):
            x = sorted(nd.glob(f"Enhancers.{feat}.*CountReads.bedgraph"))
            y = sorted(nd.glob(f"Genes.TSS1kb.{feat}.*CountReads.bedgraph"))
            z = sorted(nd.glob(f"Genes.{feat}.*CountReads.bedgraph"))
            if x and y and z:
                rd = lambda p: pd.read_csv(p, sep="\t", header=None).iloc[:, 3].sum()  # noqa: E731
                m[f"countsEnhancers_{feat}"] = float(rd(x[0]))
                m[f"countsGeneTSS_{feat}"] = float(rd(y[0]))
                m[f"countsGenes_{feat}"] = float(rd(z[0]))
    summ = outdir / f"QCSummary_{qc_name}.tsv"
    with open(summ, "w") as fh:
        for k, v in m.items():
            fh.write(f"{k}\t{v}\n")
    emit("TSV", summ)
    return {"qc_summary": str(summ), "metrics": m, "figures": [str(f) for f in figs]}


# ---------------------------------------------------------------------------
# Hi-C utilities (compute_powerlaw_fit_from_hic, makeAverageHiC, extract_avg_hic, juicebox_dump, make_bedgraph_from_HiC)
# ---------------------------------------------------------------------------

def load_hic_for_powerlaw(chromosomes, hic_dir, hic_type: str, resolution: int, min_window: int, max_window: int,
                          upstream_compat: bool = False):
    pd, np = _pd(), _np()
    parts = []
    for c in chromosomes:
        try:
            f, norm, vc = get_hic_file(c, hic_dir, allow_vc=False, hic_type=hic_type)
            if hic_type == "juicebox":
                d = load_hic_juicebox(f, norm, vc, resolution, 100, max_window, min_window)
                d["dist_for_fit"] = (d["bin1"] - d["bin2"]).abs() * resolution
            elif hic_type == "bedpe":
                d = load_hic_bedpe(f)
                raw = ((d["x2"] + d["x1"]) / 2 - (d["y2"] + d["y1"]) / 2).abs()
                d["dist_for_fit"] = (raw // resolution) * resolution
                d = d[(d["dist_for_fit"] >= min_window) & (d["dist_for_fit"] <= max_window)]
            elif hic_type == "avg":
                d = load_hic_avg(f, resolution)
                if upstream_compat:
                    d["dist_for_fit"] = ((d["bin1"] - d["bin2"]) / resolution).abs() * resolution
                    d["dist_for_fit"] = d["dist_for_fit"].astype(int)
                else:
                    d["dist_for_fit"] = (d["bin1"] - d["bin2"]).abs() * resolution
                    d = d[(d["dist_for_fit"] >= min_window) & (d["dist_for_fit"] <= max_window)]
            else:
                raise ValueError("invalid --hic-type")
            parts.append(d)
        except (FileNotFoundError, OSError, ValueError) as e:
            print(f"skipping {c}: {e}")
    if not parts:
        raise SystemExit(f"no Hi-C data loaded from {hic_dir}")
    return pd.concat(parts, ignore_index=True)


def do_powerlaw_fit(H, resolution: int):
    np = _np()
    from scipy import stats  # type: ignore
    summ = H.groupby("dist_for_fit").agg({"hic_contact": "sum"})
    n_bins = int((H["dist_for_fit"] == resolution).sum())
    summ["hic_contact"] = summ["hic_contact"] / n_bins
    pc = 1e-6
    res = stats.linregress(np.log(summ.index.to_numpy(float) + pc), np.log(summ["hic_contact"].to_numpy(float) + pc))
    mv = H.groupby("dist_for_fit").agg({"hic_contact": ["mean", "var"]})
    mv.columns = ["mean", "var"]
    return float(res.slope), float(res.intercept), mv, float(res.rvalue ** 2)


def process_celltype_for_average(ct_dir: Path, chromosome: str, resolution: int, ref_scale: float, ref_gamma: float):
    """makeAverageHiC.process_chr: DS KR matrix, no diagonal correction, rescaled to the reference power law."""
    pd, np = _pd(), _np()
    f, norm, vc = get_hic_file(chromosome, ct_dir, allow_vc=True)
    if vc:
        return None
    plf = ct_dir / "powerlaw" / "hic.powerlaw.txt"
    if not plf.exists():
        plf = ct_dir / "powerlaw" / "hic.powerlaw.tsv"
    pl = pd.read_csv(plf, sep="\t")
    if "pl_gamma" in pl.columns:
        gamma, scale = float(pl["pl_gamma"].iloc[0]), float(pl["pl_scale"].iloc[0])      # upstream: negative gamma
    else:
        gamma, scale = -float(pl["hic_gamma"].iloc[0]), float(pl["hic_scale"].iloc[0])
    h = load_hic_juicebox(f, norm, False, resolution, float("nan"), float("inf"), 0, apply_diagonal_bin_correction=False)
    h = h.rename(columns={"hic_contact": "hic_kr"})
    dists = (h["bin2"] - h["bin1"]) * resolution
    h["hic_kr"] = h["hic_kr"] * (get_powerlaw_at_distance(dists, -ref_gamma, ref_scale) / get_powerlaw_at_distance(dists, -gamma, scale))
    return h


def make_average_hic(celltype_dirs: "List[Path]", chromosome: str, outdir: Path, resolution: int = 5000,
                     ref_scale: float = 5.41, ref_gamma: float = -0.876, min_cell_types: int = 3):
    pd, np = _pd(), _np()
    special = np.inf
    hs = []
    for d in celltype_dirs:
        h = process_celltype_for_average(Path(d), chromosome, resolution, ref_scale, ref_gamma)
        if h is None:
            continue
        h.loc[np.isnan(h["hic_kr"]), "hic_kr"] = special
        hs.append(h.set_index(["bin1", "bin2"])["hic_kr"])
    if not hs:
        raise SystemExit("no KR-normalised cell types found")
    allh = pd.concat(hs, axis=1, join="outer")
    allh = allh.fillna(0).replace(special, np.nan)
    avg = allh.mean(axis=1)
    good = allh.shape[1] - allh.isna().sum(axis=1)
    out = pd.DataFrame({"bin1": allh.index.get_level_values(0) * resolution, "bin2": allh.index.get_level_values(1) * resolution,
                        "avg_hic": avg.to_numpy()})
    out.loc[good.to_numpy() < min_cell_types, "avg_hic"] = np.nan
    out = out[(out["avg_hic"] > 0) | out["avg_hic"].isna()]
    p = Path(outdir) / chromosome / f"{chromosome}.avg.gz"
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(p, sep="\t", header=False, index=False, compression="gzip", na_rep="nan")
    return p, out


def split_avg_hic(avg_bed, outdir: Path) -> "Dict[str, Path]":
    base = Path(outdir) / "AvgHiC"
    base.mkdir(parents=True, exist_ok=True)
    handles: "Dict[str, Any]" = {}
    paths: "Dict[str, Path]" = {}
    with _open_text(avg_bed) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, rest = line.split("\t", 1)
            if chrom not in handles:
                (base / chrom).mkdir(parents=True, exist_ok=True)
                paths[chrom] = base / chrom / f"{chrom}.bed.gz"
                handles[chrom] = gzip.open(paths[chrom], "wt")
            handles[chrom].write(rest)
    for h in handles.values():
        h.close()
    return paths


def juicebox_dump_commands(juicebox: str, hic_file: str, outdir: Path, resolution: int, chromosomes: "List[str]",
                           include_raw: bool = False) -> "List[str]":
    cmds = []
    for c in chromosomes:
        c = str(c).replace("chr", "")
        od = f"{outdir}/chr{c}/"
        cmds.append(f"mkdir -p {od}")
        cmds.append(f"{juicebox} dump observed KR {hic_file} {c} {c} BP {resolution} {od}chr{c}.KRobserved")
        cmds.append(f"gzip {od}chr{c}.KRobserved")
        cmds.append(f"{juicebox} dump norm KR {hic_file} {c} BP {resolution} {od}chr{c}.KRnorm")
        cmds.append(f"gzip {od}chr{c}.KRnorm")
        if include_raw:
            cmds.append(f"{juicebox} dump observed NONE {hic_file} {c} {c} BP {resolution} {od}chr{c}.RAWobserved")
            cmds.append(f"gzip {od}chr{c}.RAWobserved")
    return cmds


def hic_bedgraphs(genes_path, hic_dir, outdir: Path, resolution: int = 5000, window: int = 5_000_000, kr_cutoff: float = 0.1,
                  gene_name_annotations="symbol", primary_gene_identifier="symbol", overwrite: bool = False) -> "List[Path]":
    """make_bedgraph_from_HiC: per gene, the Hi-C row at the TSS bin within `window`; bins whose KR norm is below
    kr_cutoff are NaN, NaNs interpolated linearly, leading/trailing NaN -> 0."""
    pd, np = _pd(), _np()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    genes = process_gene_bed(read_bed(genes_path), gene_name_annotations, primary_gene_identifier, fail_on_nonunique=False)
    written = []
    cache: "Dict[str, Any]" = {}
    for _, g in genes.iterrows():
        c = str(g["chr"])
        if c not in cache:
            try:
                f, norm, vc = get_hic_file(c, hic_dir, allow_vc=True)
            except FileNotFoundError:
                cache[c] = None
                print(f"No HiC data for {c}")
                continue
            m = hic_to_sparse(f, norm, resolution).tocsr()
            nv = pd.read_csv(norm, header=None)[0].to_numpy(float)
            cache[c] = (m, nv)
        if cache[c] is None:
            continue
        m, nv = cache[c]
        fn = outdir / f"{g['name'] or 'UNK'}_{c}_{int(g['tss'])}.bg.gz"
        if fn.exists() and not overwrite:
            continue
        tb = int(g["tss"]) // resolution
        if tb >= m.shape[0]:
            continue
        row = np.asarray(m[tb].todense()).ravel().astype(float)
        bad = ~(nv >= kr_cutoff)
        row[bad[: row.size]] = np.nan
        if nv[tb] < kr_cutoff:
            row[:] = np.nan
        idx = np.arange(row.size)
        keep = np.abs(idx * resolution - int(g["tss"])) < window
        df = pd.DataFrame({"chr": c, "start": idx[keep] * resolution, "end": (idx[keep] + 1) * resolution, "v": row[keep]})
        df["v"] = df["v"].interpolate().fillna(0)
        df.to_csv(fn, sep="\t", header=False, index=False, compression="gzip")
        written.append(fn)
    return written


def compare_prediction_files(test_file, expected_file, columns: "Optional[Dict[str, str]]" = None, atol: float = 1e-6):
    """tests/test_full_abc_run.py: the tested columns, rows matched on (chr, start, end, TargetGene)."""
    pd, np = _pd(), _np()
    cols = columns or {"chr": "str", "start": "int64", "end": "int64", "name": "str", "class": "str", "TargetGene": "str",
                       "ABC.Score.Numerator": "float64", "ABC.Score": "float64", "powerlaw.Score": "float64"}
    a = pd.read_csv(test_file, sep="\t")
    b = pd.read_csv(expected_file, sep="\t")
    use = [c for c in cols if c in a.columns and c in b.columns]
    key = [c for c in ["chr", "start", "end", "TargetGene"] if c in use]
    a = a[use].astype({c: cols[c] for c in use}); b = b[use].astype({c: cols[c] for c in use})
    m = a.merge(b, on=key, how="outer", suffixes=("_test", "_expected"), indicator=True)
    res = {"n_test": int(len(a)), "n_expected": int(len(b)), "n_matched": int((m["_merge"] == "both").sum()),
           "n_only_test": int((m["_merge"] == "left_only").sum()), "n_only_expected": int((m["_merge"] == "right_only").sum()),
           "columns": {}}
    both = m[m["_merge"] == "both"]
    ok = res["n_only_test"] == 0 and res["n_only_expected"] == 0
    for c in use:
        if c in key:
            continue
        x, y = both[c + "_test"], both[c + "_expected"]
        if cols[c] == "float64":
            d = np.abs(x.to_numpy(float) - y.to_numpy(float))
            nanmis = int((np.isnan(x.to_numpy(float)) != np.isnan(y.to_numpy(float))).sum())
            mx = float(np.nanmax(d)) if np.isfinite(d).any() else 0.0
            ncl = int((d > atol).sum())
            r = float(np.corrcoef(x.fillna(0), y.fillna(0))[0, 1]) if len(x) > 2 and x.std() > 0 and y.std() > 0 else float("nan")
            res["columns"][c] = {"max_abs_diff": mx, "n_differ": ncl, "nan_mismatch": nanmis, "pearson": r}
            ok = ok and ncl == 0 and nanmis == 0
        else:
            nd = int((x.astype(str) != y.astype(str)).sum())
            res["columns"][c] = {"n_differ": nd}
            ok = ok and nd == 0
    res["identical_within_atol"] = bool(ok)
    return res, both


def fragments_to_tagalign(fragments, out: Path) -> "Tuple[Path, int]":
    pd = _pd()
    raw = pd.read_csv(fragments, sep="\t", header=None, comment="#", usecols=[0, 1, 2], names=["chr", "start", "end"])
    tags = fragments_to_tags(raw)
    tags["name"] = "N"; tags["score"] = 1000
    tags = tags[["chr", "start", "end", "name", "score", "strand"]]
    tags = tags.sort_values(["chr", "start", "end"], kind="mergesort")
    tags.to_csv(out, sep="\t", header=False, index=False, compression="gzip")
    if shutil.which("tabix") and shutil.which("bgzip"):
        print("note: rewrite with bgzip + tabix -p bed for indexed access")
    return Path(out), int(len(tags))


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def ref_default(name: str) -> "Optional[Path]":
    p = REF_DIR / name
    return p if p.exists() else None


def resolve_qnorm(args):
    if getattr(args, "no_qnorm", False):
        return None
    if getattr(args, "qnorm", None):
        return args.qnorm
    p = ref_default("EnhancersQNormRef.K562.txt")
    if p is None:
        print("WARNING: no --qnorm and no reference/EnhancersQNormRef.K562.txt (run `abc-pipeline setup`); "
              "activity is not quantile-normalised, so the calibrated abc_thresholds do not strictly apply")
    return p


def cmd_setup(args) -> int:
    dest = Path(args.dest) if args.dest else REF_DIR
    rows = []
    for rel, desc in REFERENCE_FILES.items():
        url = RAW + rel
        out = dest / rel.replace("reference/", "", 1)
        if args.dry_run:
            print(f"would fetch {url} -> {out}")
            rows.append((rel, desc, "dry-run"))
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists() and not args.force:
            rows.append((rel, desc, "present"))
            continue
        with urllib.request.urlopen(url, timeout=120) as r:
            out.write_bytes(r.read())
        emit("Wrote", out)
        rows.append((rel, desc, f"{out.stat().st_size:,} bytes"))
    print(md_table(["file", "what", "status"], rows))
    return 0


def cmd_call_peaks(args) -> int:
    d = out_dir_for(args, "peaks")
    np_path, method = call_peaks(args.accessibility, args.chrom_sizes, d, args.pval, args.genome_size,
                                 python_fallback=args.python_fallback)
    summary = {"method": method, "narrowpeak_sorted": np_path, "accessibility": args.accessibility,
               "command": " ".join(macs2_command(args.accessibility, d, args.pval, args.genome_size))}
    if np_path is not None:
        summary["n_peaks"] = int(len(read_narrowpeak(np_path)))
    write_report(d, "ABC peak calling", [("Method", f"`{method}`" + (" (MACS-like approximation)" if method == "python-fallback" else "")),
                                         ("Command", "```\n" + summary["command"] + "\n```")], summary)
    return 0


def cmd_fragments(args) -> int:
    out = Path(args.out)
    p, n = fragments_to_tagalign(args.fragments, out)
    emit("Wrote", p)
    print(f"  {n:,} tags from {n // 2:,} fragments")
    return 0


def cmd_candidate_regions(args) -> int:
    d = out_dir_for(args, "candidate_regions")
    s = make_candidate_regions(args.narrowpeak, args.accessibility, args.chrom_sizes, d, args.includelist, args.blocklist,
                               args.n_strongest_peaks, args.peak_extend_from_summit, args.ignore_summits, args.min_peak_width)
    write_report(d, "ABC candidate regions", [("Counts", md_table(["metric", "value"], [(k, _fmt(v)) for k, v in s.items() if not isinstance(v, str)]))], s)
    return 0


def cmd_neighborhoods(args) -> int:
    d = out_dir_for(args, "neighborhoods")
    s = run_neighborhoods(args.candidate_regions, args.genes, args.chrom_sizes, d, dhs=args.dhs, atac=args.atac,
                          h3k27ac=args.h3k27ac, default_feature=args.default_accessibility_feature,
                          ubiquitous=args.ubiquitously_expressed_genes or ref_default("UbiquitouslyExpressedGenes.txt"),
                          expression_tables=(args.expression_table.split(",") if args.expression_table else []),
                          qnorm=resolve_qnorm(args), qnorm_method=args.qnorm_method, separate_promoters=not args.no_separate_promoters,
                          cell_type=args.cell_type, gene_name_annotations=args.gene_name_annotations,
                          primary_gene_identifier=args.primary_gene_identifier, genes_for_class_assignment=args.genes_for_class_assignment,
                          tss_slop_for_class_assignment=args.tss_slop_for_class_assignment, skip_gene_counts=args.skip_gene_counts,
                          skip_rpkm_quantile=args.skip_rpkm_quantile, supplementary_features=args.supplementary_features,
                          enhancer_class_override=args.enhancer_class_override, upstream_compat=args.upstream_compat)
    write_report(d, "ABC neighborhoods", [("Elements", md_table(["class", "n"], sorted(s["classes"].items()))),
                                          ("Features", md_table(["feature", "files"], [(k, ", ".join(v)) for k, v in s["features"].items()]))], s)
    return 0


def _predict_opts_from_args(args):
    return predict_opts(score_column=args.score_column, accessibility_feature=args.accessibility_feature, cellType=args.cell_type,
                        hic_file=args.hic_file, hic_resolution=args.hic_resolution, hic_pseudocount_distance=args.hic_pseudocount_distance,
                        hic_type=args.hic_type, hic_is_doubly_stochastic=args.hic_is_doubly_stochastic,
                        scale_hic_using_powerlaw=not args.no_scale_hic_using_powerlaw, hic_gamma=args.hic_gamma, hic_scale=args.hic_scale,
                        hic_gamma_reference=args.hic_gamma_reference, expression_cutoff=args.expression_cutoff,
                        promoter_activity_quantile_cutoff=args.promoter_activity_quantile_cutoff, window=args.window,
                        tss_hic_contribution=args.tss_hic_contribution, tss_slop=args.tss_slop, chromosomes=args.chromosomes,
                        include_chrY=args.include_chry)


def predictions_figure(all_putative_path, d: Path, score_col: str = "ABC.Score", threshold: "Optional[float]" = None):
    plt = _plt()
    if plt is None:
        return None
    pd, np = _pd(), _np()
    p = pd.read_csv(all_putative_path, sep="\t", usecols=lambda c: c in ("distance", score_col, "class", "isSelfPromoter"))
    p = p[~p["isSelfPromoter"].astype(bool) & (p[score_col] > 0)]
    if len(p) == 0:
        return None
    p = p.sample(min(len(p), 20000), random_state=0)
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.scatter(np.log10(p["distance"].clip(lower=1)), np.log10(p[score_col]), s=4, lw=0, alpha=0.4, color=BLUE)
    if threshold:
        ax.axhline(np.log10(threshold), color=ORANGE, lw=1.2, label=f"threshold {threshold}")
        ax.legend(fontsize=7, frameon=False)
    _style(ax, f"{score_col} vs distance (non-self-promoter pairs)", "log10 distance (bp)", f"log10 {score_col}")
    return _save(fig, d / "abc_vs_distance.png")


def cmd_predict(args) -> int:
    d = out_dir_for(args, "predict")
    opts = _predict_opts_from_args(args)
    s = predict(args.enhancers, args.genes, args.chrom_sizes, d, opts)
    if not args.no_plots:
        predictions_figure(s["all_putative"], d, args.score_column)
    write_report(d, "ABC predictions", [("Pairs", md_table(["metric", "value"], [(k, _fmt(v) if not isinstance(v, (str, list)) else v) for k, v in s.items()]))], s)
    return 0


def cmd_filter(args) -> int:
    d = out_dir_for(args, "filter")
    thr, why = determine_threshold(args.accessibility_feature, args.has_h3k27ac, args.hic_type, args.threshold,
                                   args.thresholds_table, args.upstream_compat)
    print(f"threshold {thr} ({why})")
    nonexp = args.pred_nonexpressed_file
    if nonexp is None:
        guess = Path(str(args.pred_file).replace("EnhancerPredictionsAllPutative.tsv.gz", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"))
        nonexp = guess if guess.exists() and str(guess) != str(args.pred_file) else None
    s = filter_predictions(args.pred_file, nonexp, d, args.score_column, thr, not args.exclude_self_promoter, args.only_expressed_genes)
    s["threshold_source"] = why
    write_report(d, "ABC thresholded predictions", [("Threshold", f"{thr} -- {why}"),
                                                    ("Kept", md_table(["metric", "value"], [(k, _fmt(v)) for k, v in s.items() if isinstance(v, (int, float))]))], s)
    return 0


def cmd_variant_overlap(args) -> int:
    pd = _pd()
    d = out_dir_for(args, "variant_overlap")
    p, n = variant_overlap(pd.read_csv(args.all_putative, sep="\t"), d, args.score_column)
    write_report(d, "ABC variant-overlap predictions", [("Rows", str(n))], {"file": p, "n_rows": n})
    return 0


def cmd_qc(args) -> int:
    d = out_dir_for(args, "qc")
    name = args.qc_name or Path(args.preds_file).name.replace("EnhancerPredictionsFull_", "").replace(".tsv", "")
    s = generate_qc(args.preds_file, args.candidate_regions, args.neighborhood_dir, args.chrom_sizes, d, name, args.hic_gamma,
                    args.hic_scale, args.narrowpeak, plots=not args.no_plots)
    write_report(d, "ABC QC", [("QCSummary", md_table(["metric", "value"], [(k, _fmt(v)) for k, v in s["metrics"].items()], 60))], s)
    return 0


def cmd_powerlaw_fit(args) -> int:
    pd = _pd()
    d = out_dir_for(args, "powerlaw_fit")
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"] if args.chr == "all" else args.chr.split(",")
    H = load_hic_for_powerlaw(chroms, args.hic_dir, args.hic_type, args.hic_resolution, args.min_window, args.max_window, args.upstream_compat)
    slope, intercept, mv, r2 = do_powerlaw_fit(H, args.hic_resolution)
    res = pd.DataFrame({"resolution": [args.hic_resolution], "maxWindow": [args.max_window], "minWindow": [args.min_window],
                        "hic_gamma": [-slope], "hic_scale": [intercept]})
    write_tsv(res, d / "hic.powerlaw.tsv", float_format=None)
    mv.to_csv(d / "hic.mean_var.tsv", sep="\t", index=True, header=True)
    emit("TSV", d / "hic.mean_var.tsv")
    if not args.no_plots and _plt() is not None:
        plt, np = _plt(), _np()
        fig, ax = plt.subplots(figsize=(5.2, 3.6))
        x = mv.index.to_numpy(float); y = mv["mean"].to_numpy(float)
        ok = (x > 0) & (y > 0)
        ax.scatter(np.log(x[ok]), np.log(y[ok]), s=8, color=BLUE, lw=0)
        ax.plot(np.log(x[ok]), intercept + slope * np.log(x[ok]), color="#e34948", lw=1.3, label=f"gamma {-slope:.3f}")
        ax.legend(fontsize=7, frameon=False)
        _style(ax, "Hi-C contact vs distance", "ln distance", "ln mean contact")
        _save(fig, d / "powerlaw_fit.png")
    s = {"hic_gamma": -slope, "hic_scale": intercept, "r2": r2, "n_entries": int(len(H)), "chromosomes": sorted(set(chroms))}
    write_report(d, "Hi-C power-law fit", [("Fit", md_table(["hic_gamma", "hic_scale", "r2"], [(_fmt(-slope, 6), _fmt(intercept, 6), _fmt(r2))]))], s)
    return 0


def cmd_average_hic(args) -> int:
    d = out_dir_for(args, "average_hic")
    if args.celltype_dirs:
        dirs = [Path(x) for x in args.celltype_dirs]
    else:
        dirs = [Path(args.basedir) / ct / "5kb_resolution_intra" for ct in args.celltypes.split(",")]
    write_params(d, "params.txt", {"celltypes": args.celltypes, "chromosome": args.chromosome, "basedir": args.basedir,
                                   "resolution": args.resolution, "ref_scale": args.ref_scale, "ref_gamma": args.ref_gamma,
                                   "min_cell_types_required": args.min_cell_types_required})
    p, out = make_average_hic(dirs, args.chromosome, d, args.resolution, args.ref_scale, args.ref_gamma, args.min_cell_types_required)
    emit("Wrote", p)
    s = {"file": p, "n_entries": int(len(out)), "n_nan": int(out["avg_hic"].isna().sum()), "n_cell_types": len(dirs)}
    write_report(d, "Average Hi-C", [("Entries", md_table(["metric", "value"], [(k, _fmt(v)) for k, v in s.items() if k != "file"]))], s)
    return 0


def cmd_split_avg_hic(args) -> int:
    paths = split_avg_hic(args.avg_hic_bed_file, Path(args.output_dir))
    for p in paths.values():
        emit("Wrote", p)
    return 0


def cmd_juicebox_dump(args) -> int:
    chroms = list(range(1, 23)) + ["X"] if args.chromosomes == "all" else args.chromosomes.split(",")
    cmds = juicebox_dump_commands(args.juicebox, args.hic_file, Path(args.outdir), args.resolution, chroms, args.include_raw)
    if args.skip_gzip:
        cmds = [c for c in cmds if not c.startswith("gzip ")]
    exe = args.juicebox.split()[0]
    if not shutil.which(exe) and not (exe == "java" and shutil.which("java")):
        print(f"{exe} not found; commands that would run:")
        for c in cmds:
            print("  " + c)
        return 0
    for c in cmds:
        print("Running: " + c)
        subprocess.run(c, shell=True, check=True)
    return 0


def cmd_hic_bedgraph(args) -> int:
    files = hic_bedgraphs(args.genes, args.hic_dir, Path(args.outdir), args.resolution, args.window, args.kr_cutoff,
                          args.gene_name_annotations, args.primary_gene_identifier, args.overwrite)
    for f in files[:20]:
        emit("Wrote", f)
    print(f"  {len(files)} bedgraphs")
    return 0


def cmd_compare(args) -> int:
    d = out_dir_for(args, "compare")
    res, both = compare_prediction_files(args.test, args.expected, atol=args.atol)
    rows = [(c, _fmt(v.get("max_abs_diff", float("nan"))), v.get("n_differ"), _fmt(v.get("pearson", float("nan")))) for c, v in res["columns"].items()]
    if not args.no_plots and _plt() is not None and "ABC.Score_test" in both.columns and len(both):
        plt = _plt()
        fig, ax = plt.subplots(figsize=(4.2, 4.0))
        ax.scatter(both["ABC.Score_expected"], both["ABC.Score_test"], s=5, lw=0, alpha=0.5, color=BLUE)
        lim = [0, float(max(both["ABC.Score_expected"].max(), both["ABC.Score_test"].max()))]
        ax.plot(lim, lim, color=AXIS, lw=1)
        _style(ax, "ABC.Score: test vs expected", "expected", "test")
        _save(fig, d / "compare_abc_score.png")
    verdict = "identical within atol" if res["identical_within_atol"] else "DIFFERENT"
    print(f"  {res['n_matched']:,} matched pairs; {verdict}")
    write_report(d, "ABC prediction comparison", [("Verdict", verdict),
                                                  ("Rows", f"test {res['n_test']:,}, expected {res['n_expected']:,}, matched {res['n_matched']:,}, "
                                                           f"only-test {res['n_only_test']:,}, only-expected {res['n_only_expected']:,}"),
                                                  ("Columns", md_table(["column", "max |diff|", "n differ", "pearson"], rows))], res)
    return 0 if (res["identical_within_atol"] or not args.fail_on_diff) else 1


# ---------------------------------------------------------------------------
# Whole pipeline (Snakefile)
# ---------------------------------------------------------------------------

BIOSAMPLE_COLS = ["biosample", "DHS", "ATAC", "H3K27ac", "default_accessibility_feature", "HiC_file", "HiC_type", "HiC_resolution", "alt_TSS", "alt_genes"]


def validate_biosample(row: dict) -> None:
    """utils.smk _validate_accessibility_feature / _validate_hic_info."""
    if not _blank(row.get("DHS")) and not _blank(row.get("ATAC")):
        raise SystemExit(f"{row['biosample']}: Can only specify one of DHS or ATAC for accessibility")
    if _blank(row.get("DHS")) and _blank(row.get("ATAC")):
        raise SystemExit(f"{row['biosample']}: Must provide either DHS or ATAC accessibility file")
    if not _blank(row.get("HiC_file")):
        if _blank(row.get("HiC_type")) or _blank(row.get("HiC_resolution")):
            raise SystemExit(f"{row['biosample']}: Must provide HiC type and resolution with file")
        if int(float(row["HiC_resolution"])) != 5000:
            raise SystemExit(f"{row['biosample']}: Only 5kb resolution supported at the moment")


def run_biosample(row: dict, args, base: Path) -> "Dict[str, Any]":
    np = _np()
    validate_biosample(row)
    b = str(row["biosample"])
    bdir = base / safe_label(b)
    peaks_dir, nb_dir, pred_dir, met_dir = bdir / "Peaks", bdir / "Neighborhoods", bdir / "Predictions", bdir / "Metrics"
    for x in (peaks_dir, nb_dir, pred_dir, met_dir):
        x.mkdir(parents=True, exist_ok=True)
    acc_feature = "DHS" if not _blank(row.get("DHS")) else "ATAC"
    acc_files = str(row["DHS"] if acc_feature == "DHS" else row["ATAC"]).split(",")
    default_feature = row.get("default_accessibility_feature") if not _blank(row.get("default_accessibility_feature")) else acc_feature
    h3 = None if _blank(row.get("H3K27ac")) else str(row["H3K27ac"]).split(",")
    tss = row.get("alt_TSS") if not _blank(row.get("alt_TSS")) else (args.tss or ref_default("hg38/CollapsedGeneBounds.hg38.TSS500bp.bed"))
    genes = row.get("alt_genes") if not _blank(row.get("alt_genes")) else (args.genes or ref_default("hg38/CollapsedGeneBounds.hg38.bed"))
    if genes is None:
        raise SystemExit("no gene file: pass --genes or run `abc-pipeline setup`")
    sizes_path = args.chrom_sizes or ref_default("hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv")
    if sizes_path is None:
        raise SystemExit("no chromosome sizes: pass --chrom-sizes or run `abc-pipeline setup`")
    blocklist = None if args.no_blocklist else (args.blocklist or ref_default("hg38/GRCh38_unified_blacklist.bed"))
    sizes = read_chrom_sizes(sizes_path)
    steps: "Dict[str, Any]" = {"biosample": b}
    # peaks
    if args.narrowpeak:
        np_sorted = sort_narrowpeak(args.narrowpeak, sizes, peaks_dir / "macs2_peaks.narrowPeak.sorted")
        steps["peaks"] = {"method": "provided", "file": str(np_sorted)}
    else:
        np_sorted, method = call_peaks(acc_files, sizes_path, peaks_dir, args.macs_pval, args.genome_size, python_fallback=args.python_fallback)
        steps["peaks"] = {"method": method, "file": str(np_sorted) if np_sorted else None}
        if np_sorted is None:
            steps["stopped"] = "macs2 not available (rerun with --python-fallback or pass --narrowpeak)"
            return steps
    steps["n_peaks"] = int(len(read_narrowpeak(np_sorted)))
    steps["candidate_regions"] = make_candidate_regions(np_sorted, acc_files, sizes_path, peaks_dir, tss, blocklist,
                                                        args.n_strongest_peaks, args.peak_extend_from_summit)
    pg = processed_genes_file(genes, sizes, bdir / "processed_genes_file.bed")
    steps["neighborhoods"] = run_neighborhoods(steps["candidate_regions"]["candidate_regions"], pg, sizes_path, nb_dir,
                                               dhs=acc_files if acc_feature == "DHS" else None, atac=acc_files if acc_feature == "ATAC" else None,
                                               h3k27ac=h3, default_feature=default_feature,
                                               ubiquitous=args.ubiquitously_expressed_genes or ref_default("UbiquitouslyExpressedGenes.txt"),
                                               qnorm=resolve_qnorm(args), cell_type=None, upstream_compat=args.upstream_compat)
    hic_file = None if _blank(row.get("HiC_file")) else str(row["HiC_file"])
    hic_type = None if _blank(row.get("HiC_type")) else str(row["HiC_type"])
    hic_res = None if _blank(row.get("HiC_resolution")) else int(float(row["HiC_resolution"]))
    score_col = args.score_column if hic_file else "powerlaw.Score"      # predictions.smk: no Hi-C -> --score_column powerlaw.Score
    opts = predict_opts(score_column=score_col, accessibility_feature=default_feature, cellType=b, hic_file=hic_file,
                        hic_type=hic_type or "hic", hic_resolution=hic_res, hic_gamma=args.hic_gamma, hic_scale=args.hic_scale,
                        hic_pseudocount_distance=args.hic_pseudocount_distance, scale_hic_using_powerlaw=not args.no_scale_hic_using_powerlaw)
    steps["predict"] = predict(nb_dir / "EnhancerList.txt", nb_dir / "GeneList.txt", sizes_path, pred_dir, opts)
    thr, why = determine_threshold(default_feature, bool(h3), hic_type, args.threshold, args.thresholds_table, args.upstream_compat)
    steps["threshold"] = {"value": thr, "source": why}
    steps["filter"] = filter_predictions(steps["predict"]["all_putative"], steps["predict"]["all_putative_nonexpressed"], pred_dir,
                                         args.score_column, thr, not args.exclude_self_promoter, args.only_expressed_genes)
    steps["qc"] = generate_qc(steps["filter"]["full"], steps["candidate_regions"]["candidate_regions"], nb_dir, sizes_path, met_dir,
                              steps["filter"]["file_format"], args.hic_gamma, args.hic_scale, plots=not args.no_plots)
    if not args.no_plots:
        predictions_figure(steps["predict"]["all_putative"], pred_dir, args.score_column, thr)
    return steps


def cmd_run(args) -> int:
    pd = _pd()
    d = out_dir_for(args, "abc_run")
    if args.biosamples_table:
        t = pd.read_csv(args.biosamples_table, sep="\t", dtype=str).fillna("")
        rows = t.to_dict("records")
    else:
        if not args.biosample or not (args.dhs or args.atac):
            raise SystemExit("pass --biosamples-table, or --biosample with --dhs or --atac")
        rows = [{"biosample": args.biosample, "DHS": ",".join(args.dhs or []), "ATAC": ",".join(args.atac or []),
                 "H3K27ac": ",".join(args.h3k27ac or []), "default_accessibility_feature": args.default_accessibility_feature or "",
                 "HiC_file": args.hic_file or "", "HiC_type": args.hic_type or "", "HiC_resolution": str(args.hic_resolution or ""),
                 "alt_TSS": "", "alt_genes": ""}]
    results = [run_biosample(r, args, d) for r in rows]
    table = []
    for r in results:
        if "stopped" in r:
            table.append((r["biosample"], r["peaks"]["method"], "-", "-", "-", "-", r["stopped"]))
            continue
        m = r["qc"]["metrics"]
        table.append((r["biosample"], r["peaks"]["method"], f"{r['candidate_regions']['n_candidate_regions']:,}",
                       f"{r['predict']['n_pairs']:,}", f"{r['threshold']['value']} ({r['threshold']['source']})",
                       f"{r['filter']['n_distal_kept']:,}", _fmt(m.get("MedianEnhPerGene"))))
    sections = [("Biosamples", md_table(["biosample", "peaks", "candidates", "all pairs", "threshold", "distal kept", "median enh/gene"], table)),
                ("Layout", "`<biosample>/Peaks`, `Neighborhoods`, `Predictions`, `Metrics` as in the upstream results/ tree.")]
    write_report(d, "ABC pipeline run", sections, {"biosamples": results})
    return 0 if all("stopped" not in r for r in results) or not args.python_fallback else 1


# ---------------------------------------------------------------------------
# Self-test: synthetic chromosome, planted enhancers, power-law Hi-C with a planted loop
# ---------------------------------------------------------------------------

SYN_SIZES = {"chr1": 2_000_000, "chrX": 400_000}
SYN_GENES = [  # name, chrom, tss, strand
    ("GENEA", "chr1", 150_000, "+"), ("GENEB", "chr1", 400_000, "-"), ("GENEC", "chr1", 600_000, "+"),
    ("GLOOP", "chr1", 900_000, "+"), ("GENEE", "chr1", 1_200_000, "-"), ("GNULL", "chr1", 1_600_000, "+"),
    ("GENEG", "chr1", 1_850_000, "+"), ("GENEX", "chrX", 200_000, "+"),
]
SYN_ENH = [  # name, chrom, centre, DHS reads, target gene
    ("E_A", "chr1", 180_000, 400, "GENEA"), ("E_B", "chr1", 370_000, 300, "GENEB"), ("E_C", "chr1", 610_000, 300, "GENEC"),
    ("E_LOOP", "chr1", 1_302_500, 90, "GLOOP"), ("E_E", "chr1", 1_170_000, 400, "GENEE"), ("E_G", "chr1", 1_880_000, 300, "GENEG"),
    ("E_X", "chrX", 230_000, 200, "GENEX"),
]
SYN_GAMMA, SYN_SCALE, SYN_RES = 1.0, 5.5, 5000


def _syn_tags(rng, sites, n_bg, sd, sizes, read_len=36):
    pd, np = _pd(), _np()
    rows = []
    for chrom, centre, n in sites:
        pos = np.round(rng.normal(centre, sd, n)).astype(int)
        rows.append(pd.DataFrame({"chr": chrom, "cut": pos}))
    tot = sum(sizes.values())
    for chrom, size in sizes.items():
        k = int(n_bg * size / tot)
        rows.append(pd.DataFrame({"chr": chrom, "cut": rng.integers(read_len, size - read_len, k)}))
    df = pd.concat(rows, ignore_index=True)
    plus = rng.random(len(df)) < 0.5
    start = np.where(plus, df["cut"], df["cut"] - read_len)
    df = pd.DataFrame({"chr": df["chr"], "start": start, "end": start + read_len, "name": "N", "score": 1000, "strand": np.where(plus, "+", "-")})
    rank = {c: i for i, c in enumerate(sizes)}
    df["_r"] = df["chr"].map(rank)
    return df.sort_values(["_r", "start"]).drop(columns="_r").reset_index(drop=True)


def _syn_hic_matrix(size: int, loops=(), zero_bins=(), rng=None):
    np = _np()
    n = size // SYN_RES
    i, j = np.triu_indices(n)
    d = (j - i) * SYN_RES
    v = 10.0 * np.exp(SYN_SCALE - SYN_GAMMA * np.log(np.maximum(d, 5000) + 1))
    if rng is not None:
        v = v * np.exp(rng.normal(0, 0.03, v.size))
    v[i == j] = v[i == j] * 5            # inflated diagonal, corrected by the pipeline
    for a, b, f in loops:
        v[(i == min(a, b)) & (j == max(a, b))] *= f
    for z in zero_bins:
        v[(i == z) | (j == z)] = 0.0
    return i, j, v, n


def synthetic_world(d: Path, seed: int = 11) -> dict:
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    W: "Dict[str, Any]" = {}
    sizes_p = d / "chrom.sizes.tsv"
    sizes_p.write_text("".join(f"{c}\t{s}\n" for c, s in SYN_SIZES.items()))
    W["sizes"] = sizes_p
    genes = []
    for k, (name, chrom, tss, strand) in enumerate(SYN_GENES):
        s, e = (tss, tss + 20_000) if strand == "+" else (tss - 20_000, tss)
        genes.append((chrom, s, e, name, 0, strand, f"ENSG{k:011d}", "protein_coding"))
    gdf = pd.DataFrame(genes)
    gpath = d / "genes.bed"
    with open(gpath, "w") as fh:
        fh.write("#chr\tstart\tend\tname\tscore\tstrand\tEnsembl_ID\tgene_type\n")
        gdf.to_csv(fh, sep="\t", header=False, index=False)
    W["genes"] = gpath
    tss = pd.DataFrame([(c, t - 250, t + 250, n, 0, st, f"ENSG{k:011d}", "protein_coding")
                        for k, (n, c, t, st) in enumerate(SYN_GENES)] + [("chrUn_x", 10, 510, "UNPLACED", 0, "+", "ENSG99999999999", "x")])
    tss.to_csv(d / "tss500.bed", sep="\t", header=False, index=False)
    W["tss"] = d / "tss500.bed"
    (d / "blocklist.bed").write_text("chr1\t1000000\t1000400\n")      # a background-only stretch
    W["blocklist"] = d / "blocklist.bed"
    (d / "ubiq.txt").write_text("GENEA\nGENEE\n")
    W["ubiq"] = d / "ubiq.txt"
    dhs_sites = [(c, x, n) for _, c, x, n, _ in SYN_ENH] + [(c, t, 200) for _, c, t, _ in SYN_GENES]
    dhs_sites += [("chr1", x, 25) for x in (60_000, 280_000, 500_000, 750_000, 1_000_200, 1_450_000, 1_720_000, 1_960_000)]
    h3_sites = []
    for c, x, n in dhs_sites:
        h3_sites += [(c, x - 250, int(0.3 * n)), (c, x + 250, int(0.3 * n))]
    W["dhs_df"] = _syn_tags(rng, dhs_sites, 30_000, 40, SYN_SIZES)
    W["h3_df"] = _syn_tags(rng, h3_sites, 20_000, 150, SYN_SIZES)
    W["dhs"] = d / "DHS.sample.tagAlign.gz"
    W["h3"] = d / "H3K27ac.sample.tagAlign.gz"
    W["dhs_df"].to_csv(W["dhs"], sep="\t", header=False, index=False, compression="gzip")
    W["h3_df"].to_csv(W["h3"], sep="\t", header=False, index=False, compression="gzip")
    # synthetic quantile-normalisation reference in the upstream format
    qs = np.round(np.arange(0.01, 1.0, 0.01), 5)
    ref = []
    for cls, n, lo, hi in (("promoter", 10, 20.0, 80.0), ("nonpromoter", 150, 0.2, 60.0), ("any", 160, 0.2, 70.0)):
        for q in qs:
            ref.append((cls, q, round((1 - q) * n, 5), round(lo + 1.3 * hi * q ** 6, 5), round(lo + hi * q ** 6, 5)))
    pd.DataFrame(ref, columns=["enh_class", "quantile", "rank", "H3K27ac.RPM", "DHS.RPM"]).to_csv(d / "qnorm_ref.tsv", sep="\t", index=False)
    W["qnorm"] = d / "qnorm_ref.tsv"
    # Hi-C: power law, noisy, inflated diagonal, a planted GLOOP <-> E_LOOP loop; chrX promoter of GENEX uncovered
    hic = d / "hic_juicebox"
    avg = d / "hic_avg"
    bedpe = d / "hic_bedpe"
    loop = (900_000 // SYN_RES, 1_302_500 // SYN_RES, 40.0)
    W["loop_bins"] = loop
    W["mats"] = {}
    for c, size in SYN_SIZES.items():
        zero = [200_000 // SYN_RES] if c == "chrX" else []
        i, j, v, n = _syn_hic_matrix(size, [loop] if c == "chr1" else [], zero, rng)
        W["mats"][c] = (i, j, v, n)
        for base in (hic, avg, bedpe):
            (base / c).mkdir(parents=True, exist_ok=True)
        keep = v > 0
        pd.DataFrame({"a": i[keep] * SYN_RES, "b": j[keep] * SYN_RES, "v": v[keep]}).to_csv(
            hic / c / f"{c}.KRobserved.gz", sep="\t", header=False, index=False, compression="gzip")
        pd.DataFrame({"x": np.ones(n)}).to_csv(hic / c / f"{c}.KRnorm.gz", header=False, index=False, compression="gzip")
        pd.DataFrame({"a": i[keep] * SYN_RES, "b": j[keep] * SYN_RES, "v": v[keep]}).to_csv(
            avg / c / f"{c}.bed.gz", sep="\t", header=False, index=False, compression="gzip")
        pd.DataFrame({"c1": c, "x1": i[keep] * SYN_RES, "x2": (i[keep] + 1) * SYN_RES, "c2": c, "y1": j[keep] * SYN_RES,
                      "y2": (j[keep] + 1) * SYN_RES, "name": ".", "v": v[keep]}).to_csv(
            bedpe / c / f"{c}.bedpe.gz", sep="\t", header=False, index=False, compression="gzip")
    W["hic"], W["avg"], W["bedpe"] = hic, avg, bedpe
    return W


def _write_bam(tags, sizes, path: Path) -> "Optional[Path]":
    try:
        import pysam  # type: ignore
    except ImportError:
        return None
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": s} for c, s in sizes.items()]}
    tid = {c: k for k, c in enumerate(sizes)}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for k, (c, s, e, strand) in enumerate(zip(tags["chr"], tags["start"], tags["end"], tags["strand"])):
            a = pysam.AlignedSegment()
            a.query_name = f"r{k}"
            a.query_sequence = "A" * int(e - s)
            a.flag = 16 if strand == "-" else 0
            a.reference_id = tid[c]
            a.reference_start = int(s)
            a.mapping_quality = 60
            a.cigar = [(0, int(e - s))]
            a.query_qualities = pysam.qualitystring_to_array("I" * int(e - s))
            out.write(a)
    pysam.index(str(path))
    return path


def cmd_selftest(args) -> int:
    import tempfile
    pd, np = _pd(), _np()
    checks: "List[Tuple[bool, str]]" = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    before = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
    np_flag = ["--no-plots"] if args.no_plots else []
    t_start = time.time()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_world(d)
        sizes = read_chrom_sizes(W["sizes"])

        print("\ninterval engine")
        m = merge_intervals(pd.DataFrame({"chr": ["chr1"] * 4, "start": [0, 100, 250, 400], "end": [100, 200, 300, 500], "count": [1, 5, 2, 3]}),
                            sizes, value_col="count")
        check(m[["start", "end", "count"]].values.tolist() == [[0, 200, 5], [250, 300, 2], [400, 500, 3]],
              "merge_intervals: book-ended intervals merge (bedtools merge -c 4 -o max)")
        a = pd.DataFrame({"chr": ["chr1"] * 3, "start": [10, 100, 200], "end": [20, 150, 300]})
        b = pd.DataFrame({"chr": ["chr1", "chr1"], "start": [20, 149], "end": [30, 151]})
        check(overlaps_any(a, b).tolist() == [False, True, False], "overlaps_any: half-open overlap (end 20 does not touch start 20)")
        tags = W["dhs_df"]
        reg = pd.DataFrame({"chr": ["chr1", "chr1", "chrX"], "start": [179_800, 1_000_000, 229_800], "end": [180_200, 1_000_400, 230_200]})
        c = run_count_reads(W["dhs"], reg)
        brute = [int(((tags["chr"] == r.chr) & (tags["start"] < r.end) & (tags["end"] > r.start)).sum()) for r in reg.itertuples()]
        check(c["count"].tolist()[:2] == brute[:2] and c["count"].iloc[2] == 2 * brute[2],
              f"tagAlign counting = brute-force overlap {brute[:2]}; chrX doubled ({brute[2]} -> {c['count'].iloc[2]})")
        check(count_total(W["dhs"]) == len(tags), "tagAlign total = lines matching chr[1-9]|...|chrX|chrY")
        bam = _write_bam(tags, SYN_SIZES, d / "DHS.sample.bam")
        if bam is not None:
            cb = run_count_reads(bam, reg)
            check(cb["count"].tolist() == c["count"].tolist() and count_total(bam) == len(tags),
                  "BAM counting via pysam = tagAlign counting; total = mapped reads")
        else:
            print("  skip  BAM counting (pysam not installed)")
        bg = d / "signal.bedgraph"
        bg.write_text("chr1\t100\t200\t2.0\nchr1\t300\t400\t1.0\n")
        cg = count_reads_in_regions(bg, pd.DataFrame({"chr": ["chr1", "chr1"], "start": [150, 0], "end": [350, 1000]}))
        check(np.allclose(cg, [2.0 * 50 + 1.0 * 50, 300.0]) and count_total(bg) == 300.0, "bedGraph signal = per-base sum (bigWig 'sum' semantics)")
        fr = d / "atac_fragments.tsv.gz"
        pd.DataFrame([("chr1", 100, 201, "BC1", 1), ("chr1", 500, 600, "BC2", 2)]).to_csv(fr, sep="\t", header=False, index=False, compression="gzip")
        rc = main(["fragments-to-tagalign", "--fragments", str(fr), "--out", str(d / "frag.tagAlign.gz")])
        ft = pd.read_csv(d / "frag.tagAlign.gz", sep="\t", header=None)
        check(rc == 0 and ft[[1, 2, 5]].values.tolist() == [[100, 150, "+"], [151, 201, "-"], [500, 550, "+"], [551, 600, "-"]],
              "fragments-to-tagalign: (start, mid, +) and (mid+1, end, -) per fragment")

        print("\ncandidate regions (exact upstream logic on a hand-built narrowPeak)")
        npk = d / "hand.narrowPeak"
        pd.DataFrame([("chr1", 1000, 1400, "p1", 0, ".", 1, 1, 1, 200), ("chr1", 1300, 1600, "p2", 0, ".", 1, 1, 1, 100),
                      ("chr1", 5000, 5300, "p3", 0, ".", 1, 1, 1, 150), ("chr1", 9000, 9200, "p4", 0, ".", 1, 1, 1, 100),
                      ("chr1", 20000, 20300, "p5", 0, ".", 1, 1, 1, 150), ("chr1", 1_999_900, 2_000_000, "p6", 0, ".", 1, 1, 1, 50)]
                     ).to_csv(npk, sep="\t", header=False, index=False)
        hand = []
        for centre, n in ((1200, 50), (5150, 40), (9100, 30), (20150, 5), (1_999_950, 20)):
            hand += [("chr1", centre - 18 + k % 5, centre + 18 + k % 5, "N", 0, "+") for k in range(n)]
        pd.DataFrame(hand).to_csv(d / "hand.tagAlign.gz", sep="\t", header=False, index=False, compression="gzip")
        (d / "hand_block.bed").write_text("chr1\t9300\t9310\n")
        (d / "hand_inc.bed").write_text("#chr\tstart\tend\n chr1\t30000\t30500\nchrZ\t1\t500\n".replace(" chr1", "chr1"))
        s = make_candidate_regions(npk, [d / "hand.tagAlign.gz"], W["sizes"], d / "hand_out", d / "hand_inc.bed", d / "hand_block.bed", n_strongest=4)
        cr = pd.read_csv(s["candidate_regions"], sep="\t", header=None).values.tolist()
        check(cr == [["chr1", 950, 1650], ["chr1", 4900, 5400], ["chr1", 30000, 30500], ["chr1", 1_999_700, 2_000_000]],
              "top-N by reads, summit +/- 250 clipped to chrom end, merge, blocklist removal, include-list added (chrZ row dropped)")
        cnt = pd.read_csv(s["counts_file"], sep="\t", header=None)
        check(cnt[3].tolist() == [50, 0, 40, 30, 5, 20], "Counts.bed: reads per narrowPeak row; merged p1+p2 ranked by max")
        s2 = make_candidate_regions(npk, [d / "hand.tagAlign.gz"], W["sizes"], d / "hand_out2", None, None, n_strongest=2, ignore_summits=True)
        cr2 = pd.read_csv(s2["candidate_regions"], sep="\t", header=None).values.tolist()
        check(cr2 == [["chr1", 750, 1850], ["chr1", 4750, 5550]], "--ignore-summits: peak +/- 250, n_strongest 2")
        s3 = make_candidate_regions(npk, [d / "hand.tagAlign.gz", d / "hand.tagAlign.gz"], W["sizes"], d / "hand_out3", None, None, n_strongest=1)
        check(Path(s3["counts_file"]).name.endswith("averageAccessibility.Counts.bed") and s3["n_top_peaks"] == 1,
              "several accessibility files -> averageAccessibility.Counts.bed")

        print("\npeak calling")
        pk = d / "peaks"
        rc = main(["call-peaks", "--accessibility", str(W["dhs"]), "--chrom-sizes", str(W["sizes"]), "--outdir", str(pk / "nomacs")])
        check(rc == 0 and not (pk / "nomacs" / "macs2_peaks.narrowPeak").exists() or shutil.which("macs2"),
              "call-peaks without macs2: prints the macs2 command, exits 0")
        cmd = macs2_command([str(W["dhs"])], pk, 0.1, "hs")
        check(cmd[:4] == ["macs2", "callpeak", "-f", "BED"] and "--call-summits" in cmd and cmd[cmd.index("--shift") + 1] == "-75",
              "macs2 command matches rules/macs2.smk")
        npf, method = call_peaks([str(W["dhs"])], W["sizes"], pk / "py", 0.1, "2.4e6", force_python=True)
        peaks = read_narrowpeak(npf)
        planted = pd.DataFrame([(c, x - 1, x + 1) for _, c, x, _, _ in SYN_ENH], columns=["chr", "start", "end"])
        hit = overlaps_any(planted, peaks)
        check(method == "python-fallback" and hit.all(), f"python fallback peak caller recovers all {len(planted)} planted enhancers ({len(peaks)} peaks)")
        rk = peaks["chr"].map({c_: k_ for k_, c_ in enumerate(SYN_SIZES)})
        check((np.diff(rk.to_numpy() * 10**9 + peaks["start"].to_numpy()) >= 0).all() and len(peaks.columns) == 10,
              "narrowPeak.sorted: 10 columns, chrom-sizes order then start")

        print("\nneighborhoods")
        cand = make_candidate_regions(npf, [W["dhs"]], W["sizes"], pk / "py", W["tss"], W["blocklist"])
        crd = pd.read_csv(cand["candidate_regions"], sep="\t", header=None, names=["chr", "start", "end"])
        check(not overlaps_any(crd, read_bed(W["blocklist"], extra_cols=[])).any() and overlaps_any(planted, crd).all(),
              f"{len(crd)} candidate regions: every planted enhancer present, none touching the blocklist")
        nb = d / "nb"
        rc = main(["neighborhoods", "--candidate-regions", cand["candidate_regions"], "--genes", str(W["genes"]), "--chrom-sizes", str(W["sizes"]),
                   "--dhs", str(W["dhs"]), "--h3k27ac", str(W["h3"]), "--default-accessibility-feature", "DHS",
                   "--ubiquitously-expressed-genes", str(W["ubiq"]), "--qnorm", str(W["qnorm"]), "--outdir", str(nb)])
        gl = pd.read_csv(nb / "GeneList.txt", sep="\t")
        el = pd.read_csv(nb / "EnhancerList.txt", sep="\t")
        check(rc == 0 and (nb / "Genes.TSS1kb.DHS.DHS.sample.tagAlign.gz.CountReads.bedgraph").exists()
              and (nb / "Enhancers.H3K27ac.H3K27ac.sample.tagAlign.gz.CountReads.bedgraph").exists(),
              "count bedgraphs named <Enhancers|Genes|Genes.TSS1kb>.<F>.<file>.CountReads.bedgraph")
        gb = gl.set_index("name")
        check(gb.loc["GENEB", "tss"] == 400_000 and gb.loc["GENEA", "tss"] == 150_000, "TSS = end for - strand, start for +")
        check(bool(gb.loc["GENEA", "is_ue"]) and not bool(gb.loc["GENEB", "is_ue"]), "is_ue: first line of the ubiquitous list kept (upstream drops it)")
        pa = ((1e-4 + gl["H3K27ac.RPKM.quantile.TSS1Kb"]) * (1e-4 + gl["DHS.RPKM.quantile.TSS1Kb"])).rank(pct=True)
        check(np.allclose(pa, gl["PromoterActivityQuantile"]), "PromoterActivityQuantile = pct-rank((1e-4 + H3K27ac q)(1e-4 + DHS q))")
        need = ["DHS.RPM", "DHS.RPM.quantile", "DHS.RPKM", "DHS.RPKM.quantile", "H3K27ac.RPM", "class", "isPromoterElement",
                "promoterSymbol", "genicSymbol", "name", "normalized_dhs", "normalized_h3K27ac", "DHS.RPMpromoter.quantile",
                "DHS.RPMnonpromoter.quantile", "activity_base", "activity_base_no_qnorm"]
        check(all(c in el.columns for c in need), "EnhancerList.txt carries the upstream columns")
        tot = count_total(W["dhs"])
        rpm = 1e6 * el["DHS.DHS.sample.tagAlign.gz.readCount"] / tot
        check(np.allclose(rpm, el["DHS.RPM"], atol=1e-6), "RPM = 1e6 x readCount / total")
        cls = {n: el.loc[overlaps_any(el, planted[planted.index == k]), "class"].iloc[0] for k, (n, *_r) in enumerate(SYN_ENH)}
        prom_el = el[el["class"] == "promoter"]
        check(cls["E_C"] == "genic" and cls["E_A"] == "intergenic" and len(prom_el) >= len(SYN_GENES),
              f"classes: E_C genic (inside GENEC), E_A intergenic, {len(prom_el)} promoter elements")
        check(el["name"].iloc[0] == f"{el['class'].iloc[0]}|{el['chr'].iloc[0]}:{el['start'].iloc[0]}-{el['end'].iloc[0]}", "name = class|chr:start-end")
        ref = pd.read_csv(W["qnorm"], sep="\t")
        npro = el[el["class"] != "promoter"]
        r = ref[ref["enh_class"] == "nonpromoter"].sort_values("rank")
        k = npro["DHS.RPM"].idxmax()
        x = (1 - npro.loc[k, "DHS.RPMnonpromoter.quantile"]) * len(npro)
        exp_v = max(0.0, float(_interp(r["rank"], r["DHS.RPM"], [x])[0]))
        check(abs(el.loc[k, "normalized_dhs"] - exp_v) < 1e-5, f"qnorm rank method: top non-promoter -> ref value at rank (1-q)n ({exp_v:.3f})")
        check(np.allclose(el["activity_base"], np.sqrt(el["normalized_h3K27ac"] * el["normalized_dhs"])), "activity_base = sqrt(normalized_h3K27ac x normalized_dhs)")
        el_nq = run_qnorm(el.copy(), None)
        check(np.allclose(el_nq["normalized_dhs"], el["DHS.RPM"]), "no qnorm reference: normalized_dhs = DHS.RPM")
        gu = load_genes(W["genes"], W["ubiq"], sizes, d, [], "symbol", "symbol", None, upstream_compat=True)[0].set_index("name")
        check(not bool(gu.loc["GENEA", "is_ue"]), "--upstream-compat reproduces the dropped first ubiquitous gene")
        expr = d / "expr.tsv"
        expr.write_text("GENEA\t5\nGENEB\t0.2\n")
        ge = load_genes(W["genes"], None, sizes, d, [str(expr)], "symbol", "symbol", None)[0].set_index("name")
        check(ge.loc["GENEA", "Expression"] == 5 and np.isnan(ge.loc["GENEC", "Expression"]), "expression table merged by symbol (NaN when absent)")

        print("\npredictions: power law")
        pl = d / "pred_pl"
        rc = main(["predict", "--enhancers", str(nb / "EnhancerList.txt"), "--genes", str(nb / "GeneList.txt"), "--chrom-sizes", str(W["sizes"]),
                   "--accessibility-feature", "DHS", "--cell-type", "SYN", "--score-column", "powerlaw.Score", "--outdir", str(pl)] + np_flag)
        P = pd.concat([pd.read_csv(pl / "EnhancerPredictionsAllPutative.tsv.gz", sep="\t"),
                       pd.read_csv(pl / "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz", sep="\t")], ignore_index=True)
        exp_cols = ["chr", "start", "end", "name", "class", "activity_base", "activity_base_enh", "activity_base_squared_enh",
                    "normalized_dhs_enh", "normalized_h3k27ac_enh", "enh_midpoint", "enh_idx", "TargetGene", "TargetGeneTSS",
                    "TargetGeneExpression", "TargetGenePromoterActivityQuantile", "TargetGeneIsExpressed", "TargetGeneEnsembl_ID",
                    "normalized_dhs_prom", "normalized_h3k27ac_prom", "gene_idx", "distance", "isSelfPromoter", "powerlaw_contact",
                    "powerlaw_contact_reference", "ABC.Score.Numerator", "ABC.Score", "powerlaw.Score.Numerator", "powerlaw.Score", "CellType"]
        check(rc == 0 and list(P.columns) == exp_cols, "AllPutative columns (power-law run) = upstream predict.py output")
        mid = (P["start"] + P["end"]) / 2
        check(np.allclose(P["distance"], (mid - P["TargetGeneTSS"]).abs()) and (P["distance"] < 5e6).all(), "distance = |midpoint - TSS| < 5 Mb")
        plc = np.exp(DEFAULTS["hic_scale"] - DEFAULTS["hic_gamma"] * np.log(np.maximum(P["distance"], 5000) + 1))
        check(np.allclose(P["powerlaw_contact"], plc, rtol=1e-5, atol=1e-6), "powerlaw_contact = exp(scale - gamma ln(max(d, 5000) + 1))")
        num = P["activity_base_enh"] * P["powerlaw_contact"]
        tol6 = 1e-6 * (1 + P["activity_base_enh"] + P["powerlaw_contact"])     # %.6f rounding of both factors
        check((np.abs(P["ABC.Score.Numerator"] - num) <= tol6).all(), "ABC.Score.Numerator = activity x powerlaw contact")
        den = P.groupby(["TargetGene", "TargetGeneTSS"])["ABC.Score.Numerator"].transform("sum")
        ns = ~P["isSelfPromoter"]
        check(np.allclose(P.loc[ns, "ABC.Score"], (P["ABC.Score.Numerator"] / den)[ns], rtol=1e-3, atol=1e-6)
              and (P.loc[P["isSelfPromoter"], "ABC.Score"] == 1).all(), "ABC.Score = numerator / per-gene sum; self-promoters = 1")
        selfp = P[P["isSelfPromoter"]]
        check(((selfp["start"] - 500 < selfp["TargetGeneTSS"]) & (selfp["end"] + 500 > selfp["TargetGeneTSS"]) & (selfp["class"] == "promoter")).all()
              and set(selfp["TargetGene"]) == {g[0] for g in SYN_GENES}, "every gene has its self-promoter (class promoter, TSS within 500 bp)")
        dist = P[P["class"] != "promoter"]
        for name, _c, x, _n, gene in SYN_ENH:
            if name == "E_LOOP":
                continue
            sub = dist[dist["TargetGene"] == gene].sort_values("ABC.Score", ascending=False)
            top = sub.iloc[0]
            check(top["start"] <= x <= top["end"], f"planted {name} is the top distal element for {gene} (ABC {top['ABC.Score']:.3f})")
        best_null = dist[(dist["TargetGene"] == "GNULL") & (dist["distance"] < 250_000)]["ABC.Score"].max()
        planted_min = min(dist[(dist["TargetGene"] == g) & (dist["start"] <= x) & (dist["end"] >= x)]["ABC.Score"].max()
                          for n_, _c, x, _n, g in SYN_ENH if n_ not in ("E_LOOP",))
        check(best_null < 0.5 * planted_min, f"null gene GNULL: best distal ABC {best_null:.4f} < half the weakest planted ({planted_min:.3f})")
        expd = P.groupby("TargetGene")["TargetGeneIsExpressed"].first()
        gq = gl.set_index("symbol")["PromoterActivityQuantile"]
        check(all(bool(expd[g]) == bool(gq[g] >= 0.30) for g in expd.index), "TargetGeneIsExpressed = PromoterActivityQuantile >= 0.30 (no expression data)")
        vo = pd.read_csv(pl / "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz", sep="\t")
        okv = (((vo["powerlaw.Score"] > 0.015) & (vo["class"] != "promoter")) | ((vo["class"] == "promoter") & (vo["powerlaw.Score"] > 0.1))).all()
        m_vo = vo.merge(P[["name", "TargetGene", "start", "end"]], on=["name", "TargetGene"], suffixes=("", "_orig"))
        check(okv and len(vo) > 0 and (m_vo["start"] == m_vo["start_orig"] + 150).all() and (m_vo["end"] == m_vo["end_orig"] - 150).all(),
              f"variant-overlap file: score rule, <= 2 Mb, 150 bp trimmed ({len(vo)} rows)")

        print("\npredictions: Hi-C (juicebox, avg, bedpe, .hic records)")
        runs = {}
        for ht, hdir in (("juicebox", W["hic"]), ("avg", W["avg"]), ("bedpe", W["bedpe"])):
            od = d / f"pred_{ht}"
            rc = main(["predict", "--enhancers", str(nb / "EnhancerList.txt"), "--genes", str(nb / "GeneList.txt"), "--chrom-sizes", str(W["sizes"]),
                       "--accessibility-feature", "DHS", "--cell-type", "SYN", "--hic-file", str(hdir), "--hic-type", ht,
                       "--hic-resolution", "5000", "--outdir", str(od), "--no-plots"])
            runs[ht] = pd.concat([pd.read_csv(od / "EnhancerPredictionsAllPutative.tsv.gz", sep="\t"),
                                  pd.read_csv(od / "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz", sep="\t")], ignore_index=True)
            check(rc == 0, f"predict --hic-type {ht}")
        J = runs["juicebox"]
        hic_cols = [c for c in exp_cols if c not in ("enh_midpoint", "enh_idx", "gene_idx", "ABC.Score.Numerator", "ABC.Score",
                                                     "powerlaw.Score.Numerator", "powerlaw.Score", "CellType")]
        hic_cols += ["hic_contact", "hic_contact_pl_scaled", "hic_pseudocount", "hic_contact_pl_scaled_adj", "ABC.Score.Numerator", "ABC.Score",
                     "powerlaw.Score.Numerator", "powerlaw.Score", "CellType", "hic_contact_squared"]
        check(list(J.columns) == hic_cols, "AllPutative columns (Hi-C run): hic_contact ... hic_contact_pl_scaled_adj, hic_contact_squared; idx columns dropped")
        i, j, v, n = W["mats"]["chr1"]
        full = np.zeros((n, n)); full[i, j] = v; full[j, i] = v
        dsm = full / full.sum(axis=0).mean()
        dg = np.arange(n)
        left = np.r_[-np.inf, dsm[dg[1:], dg[1:] - 1]]; right = np.r_[dsm[dg[:-1], dg[:-1] + 1], -np.inf]
        dsm[dg, dg] = np.maximum(left, right)
        J1 = J[(J["chr"] == "chr1")]
        eb = np.floor(((J1["start"] + J1["end"]) / 2) / 5000).astype(int).to_numpy()
        tb = np.floor(J1["TargetGeneTSS"] / 5000).astype(int).to_numpy()
        check(np.allclose(J1["hic_contact"], dsm[eb, tb], rtol=1e-4, atol=1e-6),
              "juicebox: contact = doubly-stochastic matrix entry, diagonal = max(left, right) neighbour")
        pl_ref = get_powerlaw_at_distance(J["distance"], 0.87, -4.80 + 11.63 * 0.87)
        ratio = pl_ref / J["powerlaw_contact"]
        check((np.abs(J["hic_contact_pl_scaled"] - J["hic_contact"] * ratio) <= 1e-6 * (1 + ratio + J["hic_contact"] / J["powerlaw_contact"] * 1e2)).all(),
              "hic_contact_pl_scaled = hic x powerlaw_reference / powerlaw")
        pcv = float(get_powerlaw_at_distance(5000, DEFAULTS["hic_gamma"], DEFAULTS["hic_scale"]))
        check(np.allclose(J["hic_pseudocount"], np.minimum(J["powerlaw_contact"], pcv), atol=1e-6)
              and np.allclose(J["hic_contact_pl_scaled_adj"], J["hic_contact_pl_scaled"] + J["hic_pseudocount"], atol=2e-6),
              "hic_pseudocount = min(powerlaw, powerlaw(5 kb)); adj = scaled + pseudocount")
        gx = J[J["TargetGene"] == "GENEX"]
        check(np.allclose(gx["hic_contact"], get_powerlaw_at_distance(gx["distance"], DEFAULTS["hic_gamma"], DEFAULTS["hic_scale"], 5000), rtol=1e-4, atol=1e-6),
              "GENEX promoter has no Hi-C coverage -> all its pairs fall back to the power law (qc_hic)")
        lp = lambda df: df[(df["TargetGene"] == "GLOOP") & (df["start"] <= 1_302_500) & (df["end"] >= 1_302_500)]["ABC.Score"].max()  # noqa: E731
        s_pl, s_j = lp(P), lp(J)
        check(s_j > 3 * s_pl and s_j > 0.02, f"planted loop: E_LOOP -> GLOOP ABC {s_pl:.4f} (power law) -> {s_j:.4f} (Hi-C)")
        A, B = runs["avg"], runs["bedpe"]
        k_ = ["chr", "start", "end", "TargetGene"]
        mab = A.merge(B, on=k_, suffixes=("_a", "_b"))
        check(len(mab) == len(A) == len(B) and np.allclose(mab["hic_contact_a"], mab["hic_contact_b"]) and np.allclose(mab["ABC.Score_a"], mab["ABC.Score_b"]),
              "avg and bedpe loaders give identical contacts and ABC scores for the same matrix")
        recs = {c: [[a_ * SYN_RES, b_ * SYN_RES, x_] for a_, b_, x_ in zip(*W["mats"][c][:3]) if x_ > 0] for c in SYN_SIZES}
        opts = predict_opts(accessibility_feature="DHS", cellType="SYN", hic_type="hic", hic_resolution=5000)
        predict(nb / "EnhancerList.txt", nb / "GeneList.txt", W["sizes"], d / "pred_hicrec", opts, hic_records_by_chrom=recs)
        R = pd.read_csv(d / "pred_hicrec" / "EnhancerPredictionsAllPutative.tsv.gz", sep="\t")
        mr = R.merge(J, on=k_, suffixes=("_r", "_j"))
        check(len(mr) == len(R) and np.allclose(mr["hic_contact_r"], mr["hic_contact_j"], rtol=1e-4, atol=1e-6),
              ".hic streaming path (records -> diagonal fill -> / mean row sum) = juicebox doubly-stochastic path")
        pv = pd.DataFrame({"chr": ["c"] * 3, "TargetGene": ["g", "g", "h"], "TargetGeneTSS": [1, 1, 5], "isSelfPromoter": [True, False, False],
                           "Numer": [2.0, 1.0, 3.0]})
        cs = compute_score(pv, [pv["Numer"]], "T")
        check(cs["T.Score"].tolist() == [1.0, 1 / 3, 1.0], "compute_score: share of the per-(gene, TSS) sum, self-promoter forced to 1")

        print("\nthresholds and filtering")
        check(determine_threshold("DHS", True, None)[0] == 0.017 and determine_threshold("DHS", True, "hic")[0] == 0.027
              and determine_threshold("ATAC", False, None)[0] == 0.013, "abc_thresholds: DHS+H3K27ac powerlaw 0.017, intact Hi-C 0.027, ATAC powerlaw 0.013")
        check(determine_threshold("DHS", True, "avg")[0] == 0.016 and determine_threshold("DHS", True, "avg", upstream_compat=True)[0] == 0.02,
              "avg Hi-C -> 0.016 (upstream's 'avg_hic' lookup falls back to 0.02 with --upstream-compat)")
        check(determine_threshold("DHS", True, "juicebox")[0] == 0.02 and determine_threshold("DHS", True, "hic", explicit=0.05)[0] == 0.05,
              "no table row -> 0.02; explicit threshold wins")
        fo = d / "filt"
        rc = main(["filter", "--pred-file", str(d / "pred_juicebox" / "EnhancerPredictionsAllPutative.tsv.gz"), "--accessibility-feature", "DHS",
                   "--has-h3k27ac", "--hic-type", "hic", "--outdir", str(fo)])
        full_f = fo / "EnhancerPredictionsFull_threshold0.027_self_promoter.tsv"
        F = pd.read_csv(full_f, sep="\t")
        check(rc == 0 and full_f.exists() and (fo / "EnhancerPredictionsFull_threshold0.027_self_promoter.bedpe.gz").exists()
              and (fo / "GenePredictionStats_threshold0.027_self_promoter.tsv").exists(), "file names EnhancerPredictionsFull_threshold0.027_self_promoter.*")
        check((F["ABC.Score"] > 0.027).all() and ((F["class"] != "promoter") | F["isSelfPromoter"]).all() and len(F) == int(
            ((J["ABC.Score"] > 0.027) & ((J["class"] != "promoter") | J["isSelfPromoter"])).sum()),
            f"kept = score > threshold, promoters only as self-promoters ({len(F)} rows, non-expressed genes included)")
        S = pd.read_csv(fo / "EnhancerPredictions_threshold0.027_self_promoter.tsv", sep="\t")
        check(list(S.columns) == ["chr", "start", "end", "name", "TargetGene", "TargetGeneTSS", "CellType", "ABC.Score"], "slim file columns")
        bp_ = pd.read_csv(fo / "EnhancerPredictionsFull_threshold0.027_self_promoter.bedpe.gz", sep="\t", header=None)
        check(bp_.shape[1] == 10 and bp_[6].iloc[0] == f"{F['TargetGene'].iloc[0]}|{F['chr'].iloc[0]}:{F['start'].iloc[0]}-{F['end'].iloc[0]}",
              "bedpe: chr1 x1 x2 chr2 y1 y2 TargetGene|chr:start-end score . .")
        GS = pd.read_csv(fo / "GenePredictionStats_threshold0.027_self_promoter.tsv", sep="\t")
        check(list(GS.columns) == ["chr", "TargetGene", "TargetGeneTSS", "geneIsExpressed", "geneFailed", "nEnhancersConsidered", "nDistalEnhancersPredicted"]
              and "GLOOP" in set(F.loc[F["class"] != "promoter", "TargetGene"]), "GenePredictionStats columns; GLOOP's looped enhancer passes")
        fo2 = d / "filt2"
        main(["filter", "--pred-file", str(d / "pred_juicebox" / "EnhancerPredictionsAllPutative.tsv.gz"), "--threshold", "0.05",
              "--exclude-self-promoter", "--only-expressed-genes", "--outdir", str(fo2)])
        F2 = pd.read_csv(fo2 / "EnhancerPredictionsFull_threshold0.05_only_expr_genes.tsv", sep="\t")
        check((F2["class"] != "promoter").all() and F2["TargetGeneIsExpressed"].all(), "exclude self-promoters / only expressed genes -> _only_expr_genes")

        print("\nQC")
        qd = d / "qc"
        rc = main(["qc", "--preds-file", str(full_f), "--candidate-regions", cand["candidate_regions"], "--neighborhood-dir", str(nb),
                   "--chrom-sizes", str(W["sizes"]), "--outdir", str(qd)] + np_flag)
        Q = dict(pd.read_csv(qd / "QCSummary_threshold0.027_self_promoter.tsv", sep="\t", header=None).values.tolist())
        keys = ["MedianEnhPerGene", "StdEnhPerGene", "MedianGenePerEnh", "StdGenePerEnh", "MeanEnhPerGene", "MeanGenePerEnh", "MedianEGDist",
                "MeanEGDist", "StdEGDist", "NumEnhancersPerChrom", "MeanEnhancersPerChrom", "NumEnhancers", "EG10th", "EG90th", "NumPeaks",
                "MedWidth", "MeanWidth", "StdWidth", "NumCandidate", "MedWidthCandidate", "MeanWidthCandidate", "StdWidthCandidate",
                "countsEnhancers_DHS", "countsGeneTSS_DHS", "countsGenes_DHS", "countsEnhancers_H3K27ac", "countsGeneTSS_H3K27ac", "countsGenes_H3K27ac"]
        check(rc == 0 and list(Q) == keys, "QCSummary keys and order = metrics.py")
        check(int(Q["NumCandidate"]) == len(crd) and Q["NumPeaks"] == Q["NumCandidate"] and int(Q["NumEnhancers"]) == len(F[["chr", "start", "end"]].drop_duplicates()),
              "NumCandidate = candidate regions; NumPeaks mirrors it (upstream passes candidates as --macs_peaks)")
        eb_ = pd.read_csv(nb / "Enhancers.DHS.DHS.sample.tagAlign.gz.CountReads.bedgraph", sep="\t", header=None)[3].sum()
        check(float(Q["countsEnhancers_DHS"]) == float(eb_) and abs(float(Q["MedianEnhPerGene"]) - F.groupby("TargetGene").size().median()) < 1e-9,
              "countsEnhancers_DHS = sum of the Enhancers bedgraph; MedianEnhPerGene recomputed")
        check((qd / "EnhancerPerGene.tsv").exists() and (qd / "GenesPerEnhancer.tsv").exists() and (qd / "EnhancerGenePairsPerChrom.txt").exists(),
              "EnhancerPerGene.tsv, GenesPerEnhancer.tsv, EnhancerGenePairsPerChrom.txt")

        print("\nHi-C utilities")
        fd = d / "fit"
        rc = main(["powerlaw-fit", "--hic-dir", str(W["hic"]), "--hic-type", "juicebox", "--chr", "chr1", "--max-window", "200000", "--outdir", str(fd), "--no-plots"])
        fit = pd.read_csv(fd / "hic.powerlaw.tsv", sep="\t")
        check(rc == 0 and list(fit.columns) == ["resolution", "maxWindow", "minWindow", "hic_gamma", "hic_scale"]
              and abs(fit["hic_gamma"].iloc[0] - SYN_GAMMA) < 0.05, f"powerlaw-fit recovers planted gamma {SYN_GAMMA} ({fit['hic_gamma'].iloc[0]:.3f})")
        fa = d / "fit_avg"
        main(["powerlaw-fit", "--hic-dir", str(W["avg"]), "--hic-type", "avg", "--chr", "chr1", "--max-window", "200000", "--outdir", str(fa), "--no-plots"])
        ga = pd.read_csv(fa / "hic.powerlaw.tsv", sep="\t")["hic_gamma"].iloc[0]
        H_c = load_hic_for_powerlaw(["chr1"], W["avg"], "avg", 5000, 5000, 200000, upstream_compat=True)
        check(abs(ga - SYN_GAMMA) < 0.05 and H_c["dist_for_fit"].max() < 400, f"avg fit in bp ({ga:.3f}); --upstream-compat keeps upstream's distance-in-bins")
        base = d / "avgbase"
        for ct in ("A", "B", "C"):
            cd = base / ct / "5kb_resolution_intra"
            shutil.copytree(W["hic"] / "chr1", cd / "chr1")
            (cd / "powerlaw").mkdir(parents=True)
            pd.DataFrame({"pl_gamma": [-SYN_GAMMA], "pl_scale": [SYN_SCALE]}).to_csv(cd / "powerlaw" / "hic.powerlaw.txt", sep="\t", index=False)
        ad = d / "avgout"
        rc = main(["average-hic", "--celltypes", "A,B,C", "--chromosome", "chr1", "--basedir", str(base), "--outdir", str(ad)])
        AV = pd.read_csv(ad / "chr1" / "chr1.avg.gz", sep="\t", header=None, names=["bin1", "bin2", "avg"])
        ds_raw = full / full.sum(axis=0).mean()
        dd = (AV["bin2"] - AV["bin1"]).to_numpy()
        expv = ds_raw[(AV["bin1"] // 5000).to_numpy(), (AV["bin2"] // 5000).to_numpy()] * (
            get_powerlaw_at_distance(dd, 0.876, 5.41) / get_powerlaw_at_distance(dd, SYN_GAMMA, SYN_SCALE))
        check(rc == 0 and np.allclose(AV["avg"], expv, rtol=1e-4), "average-hic: mean of DS matrices rescaled to the reference power law")
        p4, o4 = make_average_hic([base / ct / "5kb_resolution_intra" for ct in "AB"], "chr1", d / "avg2", min_cell_types=3)
        check(o4["avg_hic"].isna().all(), "average-hic: fewer than min_cell_types_required -> NaN")
        gw = d / "genomewide_avg.bed.gz"
        with gzip.open(gw, "wt") as fh:
            fh.write("#chr\tx1\tx2\tv\n")
            for c in SYN_SIZES:
                with gzip.open(W["avg"] / c / f"{c}.bed.gz", "rt") as src:
                    for line in src:
                        fh.write(f"{c}\t{line}")
        rc = main(["split-avg-hic", "--avg-hic-bed-file", str(gw), "--output-dir", str(d / "split")])
        sp = pd.read_csv(d / "split" / "AvgHiC" / "chr1" / "chr1.bed.gz", sep="\t", header=None)
        orig = pd.read_csv(W["avg"] / "chr1" / "chr1.bed.gz", sep="\t", header=None)
        check(rc == 0 and sp.shape == orig.shape and np.allclose(sp[2], orig[2]), "split-avg-hic: AvgHiC/<chr>/<chr>.bed.gz round-trips")
        cmds = juicebox_dump_commands("java -jar juicer_tools.jar", "x.hic", Path("out"), 5000, ["1", "X"])
        check(cmds[1] == "java -jar juicer_tools.jar dump observed KR x.hic 1 1 BP 5000 out/chr1/chr1.KRobserved" and len(cmds) == 10,
              "juicebox-dump: KR observed + norm commands per chromosome")
        rc = main(["juicebox-dump", "--hic-file", "x.hic", "--juicebox", "juicer_tools_not_installed", "--outdir", str(d / "jd"), "--chromosomes", "1"])
        check(rc == 0, "juicebox-dump without juicer: prints commands, exits 0")
        rc = main(["hic-bedgraph", "--genes", str(W["genes"]), "--hic-dir", str(W["hic"]), "--outdir", str(d / "bgs")])
        bgl = pd.read_csv(d / "bgs" / "GLOOP_chr1_900000.bg.gz", sep="\t", header=None)
        lb = bgl[(bgl[1] == 1_300_000)][3].iloc[0]
        nbv = bgl[(bgl[1] == 1_290_000) | (bgl[1] == 1_310_000)][3].max()
        check(rc == 0 and lb > 5 * nbv, "hic-bedgraph: GLOOP's row peaks at the planted loop bin")

        print("\ncompare")
        rc = main(["compare", "--test", str(d / "pred_avg" / "EnhancerPredictionsAllPutative.tsv.gz"),
                   "--expected", str(d / "pred_bedpe" / "EnhancerPredictionsAllPutative.tsv.gz"), "--outdir", str(d / "cmp1"), "--fail-on-diff", "--no-plots"])
        check(rc == 0, "compare: avg vs bedpe predictions identical on the upstream test columns")
        res, _ = compare_prediction_files(pl / "EnhancerPredictionsAllPutative.tsv.gz", d / "pred_juicebox" / "EnhancerPredictionsAllPutative.tsv.gz")
        check(not res["identical_within_atol"] and res["columns"]["powerlaw.Score"]["n_differ"] == 0,
              "compare: power law vs Hi-C differ in ABC.Score, agree on powerlaw.Score")

        print("\nsetup")
        rc = main(["setup", "--dry-run", "--dest", str(d / "ref")])
        check(rc == 0 and not (d / "ref").exists(), f"setup --dry-run lists {len(REFERENCE_FILES)} reference files at the pinned commit")

        print("\nrun (whole Snakefile on the synthetic genome)")
        rd = d / "run"
        rc = main(["run", "--biosample", "SYN", "--dhs", str(W["dhs"]), "--h3k27ac", str(W["h3"]), "--genes", str(W["genes"]), "--tss", str(W["tss"]),
                   "--chrom-sizes", str(W["sizes"]), "--blocklist", str(W["blocklist"]), "--qnorm", str(W["qnorm"]), "--hic-file", str(W["hic"]),
                   "--hic-type", "juicebox", "--hic-resolution", "5000", "--python-fallback", "--genome-size", "2.4e6", "--outdir", str(rd)] + np_flag)
        rf = rd / "SYN" / "Predictions" / "EnhancerPredictionsFull_threshold0.02_self_promoter.tsv"
        summ = json.loads((rd / "summary.json").read_text())
        check(rc == 0 and rf.exists() and (rd / "report.md").exists() and (rd / "SYN" / "Metrics" / "QCSummary_threshold0.02_self_promoter.tsv").exists(),
              "run: Peaks/Neighborhoods/Predictions/Metrics tree; juicebox has no table row -> threshold0.02")
        RF = pd.read_csv(rf, sep="\t")
        got = {g for g in RF.loc[RF["class"] != "promoter", "TargetGene"]}
        check({e[4] for e in SYN_ENH} <= got, f"run: every planted target gene has a thresholded distal enhancer ({sorted(got)})")
        check(summ["biosamples"][0]["peaks"]["method"] == "python-fallback", "run: peaks from the fallback caller when macs2 is absent")
        bt = d / "biosamples.tsv"
        pd.DataFrame([dict(zip(BIOSAMPLE_COLS, ["SYN_ATAC", "", str(W["dhs"]), "", "ATAC", "", "", "", "", ""]))]).to_csv(bt, sep="\t", index=False)
        rd2 = d / "run2"
        rc = main(["run", "--biosamples-table", str(bt), "--genes", str(W["genes"]), "--tss", str(W["tss"]), "--chrom-sizes", str(W["sizes"]),
                   "--no-blocklist", "--no-qnorm", "--python-fallback", "--genome-size", "2.4e6", "--outdir", str(rd2), "--no-plots"])
        check(rc == 0 and (rd2 / "SYN_ATAC" / "Predictions" / "EnhancerPredictionsFull_threshold0.013_self_promoter.tsv").exists(),
              "run --biosamples-table: ATAC-only power-law row -> threshold 0.013")
        try:
            validate_biosample({"biosample": "bad", "DHS": "a.bam", "ATAC": "b.tagAlign.gz"})
            bad = False
        except SystemExit:
            bad = True
        try:
            validate_biosample({"biosample": "bad2", "DHS": "a.bam", "HiC_file": "x.hic", "HiC_type": "hic", "HiC_resolution": "10000"})
            bad2 = False
        except SystemExit:
            bad2 = True
        check(bad and bad2, "biosample validation: DHS and ATAC both set, or Hi-C resolution != 5000, are rejected")
    after = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
    for p in after - before:
        shutil.rmtree(p, ignore_errors=True)
    n_fail = sum(1 for ok, _ in checks if not ok)
    print(f"\n{len(checks) - n_fail}/{len(checks)} checks passed in {time.time() - t_start:.0f}s")
    if n_fail:
        print("selftest: FAILED")
        return 1
    print("selftest: all checks pass")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_out(p, label: str) -> None:
    p.add_argument("--label", default=label, help="run label (output Docs/ABCPipeline/<timestamp>_<label>)")
    p.add_argument("--outdir", help="write here instead of a timestamped run directory")


def _add_predict_args(p) -> None:
    p.add_argument("--score-column", default=DEFAULTS["score_column"], help="score used for the variant-overlap file")
    p.add_argument("--hic-gamma", type=float, default=DEFAULTS["hic_gamma"])
    p.add_argument("--hic-scale", type=float, default=DEFAULTS["hic_scale"])
    p.add_argument("--hic-pseudocount-distance", type=int, default=DEFAULTS["hic_pseudocount_distance"])
    p.add_argument("--no-scale-hic-using-powerlaw", action="store_true", help="do not rescale Hi-C to the reference power law (config default: rescale)")


def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent abc-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch the upstream reference/ files from the pinned commit")
    p.add_argument("--dest", help=f"destination (default {REF_DIR})")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("call-peaks", help="macs2 callpeak (rules/macs2.smk) + narrowPeak sort; --python-fallback without macs2")
    p.add_argument("--accessibility", nargs="+", required=True, help="DNase BAM / ATAC tagAlign.gz file(s)")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--pval", type=float, default=DEFAULTS["macs_pval"])
    p.add_argument("--genome-size", default=DEFAULTS["macs_genome_size"], help="hs, mm, ce, dm or a number")
    p.add_argument("--python-fallback", action="store_true", help="MACS-like Python caller when macs2 is not on PATH")
    _add_out(p, "peaks")
    p.set_defaults(func=cmd_call_peaks)

    p = sub.add_parser("fragments-to-tagalign", help="10x fragments.tsv.gz -> tagAlign.gz (two tags per fragment)")
    p.add_argument("--fragments", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_fragments)

    p = sub.add_parser("candidate-regions", help="makeCandidateRegions.py")
    p.add_argument("--narrowpeak", required=True, help="macs2 narrowPeak with summits (column 10)")
    p.add_argument("--accessibility", nargs="+", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--includelist", help="regions forced in (TSS500bp file); overrides the blocklist")
    p.add_argument("--blocklist", help="regions excluded")
    p.add_argument("--n-strongest-peaks", type=int, default=DEFAULTS["n_strongest"])
    p.add_argument("--peak-extend-from-summit", type=int, default=DEFAULTS["peak_extend"])
    p.add_argument("--ignore-summits", action="store_true")
    p.add_argument("--min-peak-width", type=int, default=DEFAULTS["min_peak_width"])
    _add_out(p, "candidate_regions")
    p.set_defaults(func=cmd_candidate_regions)

    p = sub.add_parser("neighborhoods", help="run.neighborhoods.py: EnhancerList.txt + GeneList.txt")
    p.add_argument("--candidate-regions", required=True)
    p.add_argument("--genes", required=True, help="BED6 + Ensembl_ID + gene_type")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--dhs", nargs="+")
    p.add_argument("--atac", nargs="+")
    p.add_argument("--h3k27ac", nargs="+")
    p.add_argument("--default-accessibility-feature", choices=["DHS", "ATAC"])
    p.add_argument("--ubiquitously-expressed-genes")
    p.add_argument("--expression-table", help="comma-separated gene<TAB>value tables")
    p.add_argument("--qnorm", help="quantile-normalisation reference (default: reference/EnhancersQNormRef.K562.txt if set up)")
    p.add_argument("--no-qnorm", action="store_true")
    p.add_argument("--qnorm-method", choices=["rank", "quantile"], default="rank")
    p.add_argument("--no-separate-promoters", action="store_true", help="qnorm all elements against the 'any' reference rows")
    p.add_argument("--cell-type")
    p.add_argument("--gene-name-annotations", default="symbol")
    p.add_argument("--primary-gene-identifier", default="symbol")
    p.add_argument("--genes-for-class-assignment")
    p.add_argument("--tss-slop-for-class-assignment", type=int, default=500)
    p.add_argument("--skip-gene-counts", action="store_true")
    p.add_argument("--skip-rpkm-quantile", action="store_true")
    p.add_argument("--supplementary-features", help="TSV feature_name, file")
    p.add_argument("--enhancer-class-override", help="BED chr start end class")
    p.add_argument("--upstream-compat", action="store_true", help="read the ubiquitous-genes list with a header, as upstream")
    _add_out(p, "neighborhoods")
    p.set_defaults(func=cmd_neighborhoods)

    p = sub.add_parser("predict", help="predict.py: EnhancerPredictionsAllPutative*.tsv.gz")
    p.add_argument("--enhancers", required=True, help="EnhancerList.txt")
    p.add_argument("--genes", required=True, help="GeneList.txt")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--accessibility-feature", choices=["DHS", "ATAC"], required=True)
    p.add_argument("--cell-type")
    p.add_argument("--hic-file", help=".hic file / URL (hic) or directory (juicebox, bedpe, avg)")
    p.add_argument("--hic-type", default="hic", choices=["hic", "juicebox", "bedpe", "avg"])
    p.add_argument("--hic-resolution", type=int)
    p.add_argument("--hic-is-doubly-stochastic", action="store_true")
    p.add_argument("--hic-gamma-reference", type=float, default=DEFAULTS["hic_gamma_reference"])
    p.add_argument("--expression-cutoff", type=float, default=DEFAULTS["expression_cutoff"])
    p.add_argument("--promoter-activity-quantile-cutoff", type=float, default=DEFAULTS["promoter_activity_quantile_cutoff"])
    p.add_argument("--window", type=int, default=DEFAULTS["window"])
    p.add_argument("--tss-hic-contribution", type=float, default=DEFAULTS["tss_hic_contribution"])
    p.add_argument("--tss-slop", type=int, default=DEFAULTS["tss_slop"])
    p.add_argument("--chromosomes", default="all")
    p.add_argument("--include-chry", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    _add_predict_args(p)
    _add_out(p, "predict")
    p.set_defaults(func=cmd_predict)

    def add_threshold_args(p):
        p.add_argument("--threshold", type=float, help="explicit threshold (default: abc_thresholds lookup, else 0.02)")
        p.add_argument("--thresholds-table", help="abc_thresholds.tsv (default: the pinned table, built in)")
        p.add_argument("--exclude-self-promoter", action="store_true", help="drop self-promoters too (config include_self_promoter: True)")
        p.add_argument("--only-expressed-genes", action="store_true")
        p.add_argument("--upstream-compat", action="store_true", help="upstream's avg -> 'avg_hic' threshold lookup (and headed ubiquitous list in run)")

    p = sub.add_parser("filter", help="filter_predictions.py with automatic threshold")
    p.add_argument("--pred-file", required=True, help="EnhancerPredictionsAllPutative.tsv.gz")
    p.add_argument("--pred-nonexpressed-file", help="default: the NonExpressedGenes file next to --pred-file")
    p.add_argument("--score-column", default=DEFAULTS["score_column"])
    p.add_argument("--accessibility-feature", default="DHS", choices=["DHS", "ATAC"])
    p.add_argument("--has-h3k27ac", action="store_true")
    p.add_argument("--hic-type", help="hic, avg, juicebox, bedpe; omit for power law")
    add_threshold_args(p)
    _add_out(p, "filter")
    p.set_defaults(func=cmd_filter)

    p = sub.add_parser("variant-overlap", help="getVariantOverlap.py")
    p.add_argument("--all-putative", required=True)
    p.add_argument("--score-column", default=DEFAULTS["score_column"])
    _add_out(p, "variant_overlap")
    p.set_defaults(func=cmd_variant_overlap)

    p = sub.add_parser("qc", help="grabMetrics.py + metrics.py")
    p.add_argument("--preds-file", required=True, help="EnhancerPredictionsFull_*.tsv")
    p.add_argument("--candidate-regions", required=True)
    p.add_argument("--neighborhood-dir")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--narrowpeak", help="report peak widths from real peaks (upstream reports the candidate regions twice)")
    p.add_argument("--qc-name", help="suffix of QCSummary_<name>.tsv (default from the predictions file name)")
    p.add_argument("--hic-gamma", type=float, default=DEFAULTS["hic_gamma"])
    p.add_argument("--hic-scale", type=float, default=DEFAULTS["hic_scale"])
    p.add_argument("--no-plots", action="store_true")
    _add_out(p, "qc")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("powerlaw-fit", help="compute_powerlaw_fit_from_hic.py")
    p.add_argument("--hic-dir", required=True)
    p.add_argument("--hic-type", default="juicebox", choices=["juicebox", "bedpe", "avg"])
    p.add_argument("--hic-resolution", type=int, default=5000)
    p.add_argument("--min-window", type=int, default=5000)
    p.add_argument("--max-window", type=int, default=1_000_000)
    p.add_argument("--chr", default="all")
    p.add_argument("--upstream-compat", action="store_true", help="avg: distance in bins and no window, as upstream")
    p.add_argument("--no-plots", action="store_true")
    _add_out(p, "powerlaw_fit")
    p.set_defaults(func=cmd_powerlaw_fit)

    p = sub.add_parser("average-hic", help="makeAverageHiC.py")
    p.add_argument("--celltypes", default="", help="comma-separated cell types under --basedir/<ct>/5kb_resolution_intra")
    p.add_argument("--basedir", default=".")
    p.add_argument("--celltype-dirs", nargs="+", help="explicit juicebox dirs (each with powerlaw/hic.powerlaw.txt)")
    p.add_argument("--chromosome", required=True)
    p.add_argument("--resolution", type=int, default=5000)
    p.add_argument("--ref-scale", type=float, default=5.41)
    p.add_argument("--ref-gamma", type=float, default=-0.876)
    p.add_argument("--min-cell-types-required", type=int, default=3)
    _add_out(p, "average_hic")
    p.set_defaults(func=cmd_average_hic)

    p = sub.add_parser("split-avg-hic", help="extract_avg_hic.py")
    p.add_argument("--avg-hic-bed-file", required=True)
    p.add_argument("--output-dir", default=".")
    p.set_defaults(func=cmd_split_avg_hic)

    p = sub.add_parser("juicebox-dump", help="juicebox_dump.py (juicer_tools dump; prints the commands without it)")
    p.add_argument("--hic-file", required=True)
    p.add_argument("--juicebox", default="java -jar juicer_tools.jar")
    p.add_argument("--resolution", type=int, default=5000)
    p.add_argument("--outdir", default=".")
    p.add_argument("--include-raw", action="store_true")
    p.add_argument("--chromosomes", default="all")
    p.add_argument("--skip-gzip", action="store_true")
    p.set_defaults(func=cmd_juicebox_dump)

    p = sub.add_parser("hic-bedgraph", help="make_bedgraph_from_HiC.py")
    p.add_argument("--genes", required=True)
    p.add_argument("--hic-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--resolution", type=int, default=5000)
    p.add_argument("--window", type=int, default=5_000_000)
    p.add_argument("--kr-cutoff", type=float, default=0.1)
    p.add_argument("--gene-name-annotations", default="symbol")
    p.add_argument("--primary-gene-identifier", default="symbol")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=cmd_hic_bedgraph)

    p = sub.add_parser("compare", help="compare two prediction files on the upstream test columns")
    p.add_argument("--test", required=True)
    p.add_argument("--expected", required=True)
    p.add_argument("--atol", type=float, default=1e-6)
    p.add_argument("--fail-on-diff", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    _add_out(p, "compare")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("run", help="the whole Snakefile for one biosample or a biosamples table")
    p.add_argument("--biosamples-table", help="TSV: biosample DHS ATAC H3K27ac default_accessibility_feature HiC_file HiC_type HiC_resolution alt_TSS alt_genes")
    p.add_argument("--biosample")
    p.add_argument("--dhs", nargs="+")
    p.add_argument("--atac", nargs="+")
    p.add_argument("--h3k27ac", nargs="+")
    p.add_argument("--default-accessibility-feature", choices=["DHS", "ATAC"])
    p.add_argument("--hic-file")
    p.add_argument("--hic-type", choices=["hic", "juicebox", "bedpe", "avg"])
    p.add_argument("--hic-resolution", type=int)
    p.add_argument("--genes", help="default reference/hg38/CollapsedGeneBounds.hg38.bed")
    p.add_argument("--tss", help="TSS include-list (default reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed)")
    p.add_argument("--chrom-sizes")
    p.add_argument("--blocklist")
    p.add_argument("--no-blocklist", action="store_true")
    p.add_argument("--ubiquitously-expressed-genes")
    p.add_argument("--qnorm")
    p.add_argument("--no-qnorm", action="store_true")
    p.add_argument("--narrowpeak", help="skip peak calling and use this narrowPeak")
    p.add_argument("--python-fallback", action="store_true")
    p.add_argument("--macs-pval", type=float, default=DEFAULTS["macs_pval"])
    p.add_argument("--genome-size", default=DEFAULTS["macs_genome_size"])
    p.add_argument("--n-strongest-peaks", type=int, default=DEFAULTS["n_strongest"])
    p.add_argument("--peak-extend-from-summit", type=int, default=DEFAULTS["peak_extend"])
    add_threshold_args(p)
    _add_predict_args(p)
    p.add_argument("--no-plots", action="store_true")
    _add_out(p, "abc_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic chromosome, planted enhancers and loop, every subcommand asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
