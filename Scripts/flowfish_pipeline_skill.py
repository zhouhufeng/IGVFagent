#!/usr/bin/env python3
"""CRISPRi-FlowFISH analysis workflow, FASTQ to scored enhancers (port of EngreitzLab/CRISPRi-FlowFISH).

Port of https://github.com/EngreitzLab/CRISPRi-FlowFISH (MIT, Copyright (c) 2021
Engreitz Lab; Snakemake + Python + R + bowtie + bedtools), the Engreitz-lab
workflow behind the CRISPRi-FlowFISH screens of Fulco et al. 2019 (Nat Genet)
and Nasser et al. 2021 (Nature).  Every rule in workflow/Snakefile and
workflow/rules/*.smk and every script in workflow/scripts/ was read and its
method re-derived in Python (pandas + numpy + scipy); no code was copied.
Relationship: port.  Pinned commit 44008b8842b1ffd82874fc9df669dfb8bb70b2de
(2025-07-25).  The older `igvfagent flowfish` skill is a partial clean-room
re-derivation of three steps; this module covers the whole workflow.

Definitions, as upstream computes them
  experiment keys   ExperimentIDPCRRep     = '-'.join(key cols) + '-Rep' + '-'.join(rep cols) + '-PCR' + PCRRep
                    ExperimentIDReplicates = '-'.join(key cols) + '-Rep' + '-'.join(rep cols)
                    ExperimentID           = '-'.join(key cols)
                    required columns SampleID Bin PCRRep GeneSymbol qPCRGene; SampleID unique;
                    Bin == 'Water' rejected; bins = unique Bin minus All / Neg / blank.
  fastq discovery   {fastqdir}/{SampleID}_*_R1_*.fastq.gz, else {SampleID}*.fastq.gz; >1 match is an error.
  map_reads         bowtie -v0 -p1 --all -x {index} --un unaligned.fastq (0 mismatches, all
                    alignments); nReadsAligned = alignment lines, nReadsUnaligned = unaligned
                    FASTQ lines / 4.  Without bowtie the port maps exactly the same way in
                    Python (read, or its reverse complement, contained in a MappingSequence).
  count_reads       cut -f3 | sort | uniq -c : alignments per reference (OligoID).
  count tables      per ExperimentIDPCRRep and per ExperimentIDReplicates: outer-join sample
                    counts on OligoID (missing = 0), sum samples per Bin (bins first, then the
                    other Bin values e.g. All), drop all-zero columns; bin_freq = count / column sum.
  flat table        GuideCounts.flat.tsv.gz: count, OligoID, SampleID.
  replicate corr.   per ExperimentIDReplicates and Bin: Pearson r between PCR-replicate count
                    columns (inner join on OligoID), lower triangle -> ReplicateCorrelation.tsv.
  sort params       BigFoot_noTotal (the loader estimate_effect_sizes.R uses): tsv Bin Mean Min
                    Max Count; log10 of Mean/Min/Max; sd seed 0.5; mu seed = log10 Min of the
                    3rd bin; total 0.  BigFoot (csv, Barcode + 'Total' row, StdDev): sd seed =
                    log10(1 + (StdDev/Mean)^2) of Total.  Astrios / Astrios_nototal (tsv Barcode
                    Mean Bounds '(lo, hi)' Count, Std.Dev.): sd seed from row 1, mu seed = log10
                    lower bound of the 4th bin.
  estimate_effect_sizes.R
                    sumreads / numobsbins on raw counts; bin counts rescaled to sorted cells
                    (count / column sum * sorted-cell Count); WeightedAvg = sum(cells * bin
                    mid-point) / sumcells in log10, with bin A at its upper bound and F at its
                    lower bound; per guide a log-normal MLE on the binned cells plus a
                    "seventh bin" of cells outside every gate: one EM step (MLE with the 7th
                    bin = 0, sd in [0,1]) gives p_in and 7th = sum(cells)/p_in*(1-p_in), then
                    the final MLE with sd in [0.1,1]; mean in [log10 10, log10 1e6]; both
                    started at the sort-param seeds; negative log-likelihood
                    -sum(n_b log p_b), p_b <= 0 -> 1e-10; L-BFGS-B as R optim (central
                    differences, ndeps 1e-3, lmm 5, factr 1e7, maxit 100) and stats4::mle's
                    solve(hessian) -- any error (e.g. singular Hessian for an unobserved
                    guide) falls back to method 'weighted_avg' (logMean = WeightedAvg,
                    logSD = sd seed).
  convert_to_real_space.py
                    mleAvg = 10^(logMean + logSD^2/2), mleSD = sqrt((10^(logSD^2)-1) mleAvg^2);
                    guide filter 18 <= len(GuideSequence) <= 21, OffTargetScore >= 50, #G <= 10,
                    sumcells >= minsum (0.05 in the rule), numobsbins >= minbins (4) (falls back
                    to GuideSequenceMinusG when GuideSequence is missing); mleAvg / median(NC),
                    WeightedAvg / median(NC); clamp to [0, 5]; mleAvg / mean(NC) again.
  CalculateTilingStatistic.R (window 10 as the Snakefile sets, maxSpan 750, score mleAvg)
                    sliding windows of consecutive non-control guides (sorted by chr, start),
                    windows spanning > maxSpan bp skipped; mean = window mean - control mean.
                    The helper calculateSlidingWindowStatistic lives in the private
                    LanderLab-EP-Prediction repo, so its test columns are reconstructed.
  FlowFISHtssKD.py  best (lowest mean) window overlapping TSS +/- 500 bp (TSS_Override wins);
                    FlowFISH_at_TSS = 1 + mean; Background = (FF - qPCR) / (1 - FF).
  normalize_flowfish_to_qpcr.py
                    slope = (1 - qPCR)/(1 - FF); mleAvg -> slope (x - FF) + qPCR clipped to
                    [0, 5]; mleSD -> slope x; qPCR NA -> no transform.
  collapse          non-control guides with a start, overlapped with the enhancer list
                    (cut -f1-3, sorted; >= 1 bp half-open overlap), mleAvg collapsed per peak.
  ScoreEnhancers.py n >= minGuides (3); mean; Mann-Whitney U vs controls (scipy 1.5 default:
                    one-sided, normal approx, continuity + tie correction) and Student t
                    (equal_var=True); BH FDR each; Significant = fdr.ttest < 0.05 and
                    |mean - mean.ctrl| >= mineffect (0); Regulated = fdr.utest < 0.05 and
                    mean < mean.ctrl; PeakCallingSummary per experiment.
  power             the dead `if False: # POWER` block, with its indentation bug fixed:
                    ES' = ES/(Background+1); NC scores - ES' sampled (median n guides, 400
                    reps) and BH-tested with the real peaks; power = P(fdr.ttest < FDR).
  ScreenData / KnownEnhancers
                    ScreenData.txt exactly; KnownEnhancers.FlowFISH.txt is a reconstruction
                    (FormatCRISPRiScreensForPredictions.R is external): gene TSS, distance,
                    fraction change, FDR, Significant, Regulated, cell type, reference.
  PlotGuideCounts.R (written for allele counts; adapted to guide counts): PCR-replicate
                    correlations of %Reads per ExperimentIDReplicates x Bin, per-guide CV,
                    per-bin frequency relative to the average bin.

Deviations (reproduced with --upstream-compat):
  * estimate_effect_sizes.R reads sorted-cell counts with sort.params$Count[Bin], which under
    the pinned R 3.6 indexes by factor level code, not by bin name (right only when the
    file lists exactly the bins in alphabetical order); the port matches by name.
  * the MLE pairs count columns with gates positionally; the port aligns them by name.
  * BigFoot / Astrios formats are only reachable with --sort-params-format; BigFoot returns no
    mu seed upstream (every MLE would error), the port uses log10 Min of the 3rd bin.

Subcommands
  samplesheet            validate, discover FASTQs, add experiment keys (SampleList.snakemake.tsv)
  build-index            guide FASTA from the design (+ bowtie-build if on PATH)
  map-reads              bowtie (or exact Python mapper) -> alignments, counts, alignment_stats
  count-reads            alignment file -> {SampleID}.count.txt
  count-tables           per PCR rep / experimental rep bin counts + frequencies + flat table
  replicate-correlations PCR-replicate correlations per bin
  estimate-effects       estimate_effect_sizes.R (weighted average + MLE)
  real-space             convert_to_real_space.py
  windows                CalculateTilingStatistic.R
  tss-kd                 FlowFISHtssKD.py -> ScreenInfo.txt
  normalize-qpcr         normalize_flowfish_to_qpcr.py -> scaled.txt / .bedgraph
  collapse               enhancers x guides -> collapse.bed
  score-enhancers        ScoreEnhancers.py -> FullEnhancerScore, peaks bed, PeakCallingSummary
  power                  power analysis for the peak test
  format-screen          ScreenData.txt + KnownEnhancers.FlowFISH.txt
  guide-count-plots      PlotGuideCounts.R (guide-count version)
  run                    the whole Snakefile (all_input) from a sample sheet
  selftest               synthetic screen with planted enhancers through every subcommand

Output: Docs/FlowFISHPipeline/<timestamp>_<label>/ (report.md, summary.json, results/...).

Usage:
    igvfagent flowfish-pipeline run --sample-sheet SampleSheet.tsv --design design.txt --sortparams-dir sortParams --fastq-dir fastq --experiment-keycols CellLine --replicate-keycols FlowFISHRep --qpcr qPCR.txt --genelist GeneList.txt --enhancers EnhancerList.bed --label ppif
    igvfagent flowfish-pipeline estimate-effects --counts X.bin_counts.txt --sort-params B1_S1.txt --label X
    igvfagent flowfish-pipeline score-enhancers --collapsed X.collapse.bed --scaled X.scaled.txt --expt-name X
    igvfagent flowfish-pipeline selftest --no-plots
"""
from __future__ import annotations

import argparse
import glob
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
from collections import Counter
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "FlowFISHPipeline"

UPSTREAM_REPO = "EngreitzLab/CRISPRi-FlowFISH"
UPSTREAM_COMMIT = "44008b8842b1ffd82874fc9df669dfb8bb70b2de"

REQUIRED_COLS = ["SampleID", "Bin", "PCRRep", "GeneSymbol", "qPCRGene"]
MAXMEAN = 1_000_000.0
MINMEAN = 10.0
NEG = "negative_control"
BED_COLS = ["chr", "start", "end", "name", "score", "strand", "GuideSequence", "target", "OffTargetScore",
            "OligoID", "WeightedAvg", "mleAvg", "mleSD", "sumcells"]
SCORE_COLS = ["chr", "start", "end", "mean", "sem", "n", "mean.ctrl", "sem.ctrl", "n.ctrl", "p.utest", "fdr.utest",
              "p.ttest", "fdr.ttest", "CustomTargetGeneTSS", "Significant", "Regulated"]
SORT_FORMATS = ("bigfoot_nototal", "bigfoot", "astrios", "astrios_nototal")

INK, INK2, AXIS, SURFACE = "#1f2328", "#57606a", "#8c959f", "#ffffff"
C1, C2, C3 = "#2f6fb0", "#c2562b", "#8c959f"

log = logging.getLogger("flowfish_pipeline")
_LOG_PATH: Optional[Path] = None


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"flowfish_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers, rows, max_rows: int = 40) -> str:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(_fmt(c) for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 4) -> str:
    if x is None or type(x).__name__ in ("bool", "bool_"):
        return str(bool(x)) if x is not None else "None"
    if isinstance(x, (int,)):
        return str(x)
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x).replace("|", "\\|")
    if not math.isfinite(f):
        return "nan"
    if f == int(f) and abs(f) < 1e12 and not isinstance(x, float):
        return str(int(f))
    return f"{f:.{nd}g}"


def _jsonable(o):
    np = _np()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        f = float(o)
        return f if math.isfinite(f) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return str(o)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def emit(kind: str, path: Path) -> Path:
    print(f"{kind}: {path}")
    return path


def write_tsv(df, path: Path, kind: str = "TSV", **kw) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    kw.setdefault("index", False)
    df.to_csv(path, sep="\t", compression="gzip" if str(path).endswith(".gz") else None, **kw)
    return emit(kind, path)


def finish(out: Path, title: str, summary: dict, body: str) -> int:
    summary = dict(summary)
    summary.setdefault("upstream", {"repo": UPSTREAM_REPO, "ref": UPSTREAM_COMMIT})
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=2, default=_jsonable), encoding="utf-8")
    emit("JSON", out / "summary.json")
    (out / "report.md").write_text(f"# {title}\n\n{body.strip()}\n\n_Port of {UPSTREAM_REPO} @ {UPSTREAM_COMMIT[:7]}._\n",
                                   encoding="utf-8")
    emit("Report", out / "report.md")
    return 0


def read_table(path, **kw):
    pd = _pd()
    return pd.read_csv(path, sep="\t", **kw)


def r_make_names(name: str) -> str:
    """R read.delim(check.names=TRUE) column mangling for the characters that matter here."""
    return re.sub(r"[^A-Za-z0-9._]", ".", str(name))


def parse_cols(s) -> "list[str]":
    if s is None:
        return []
    if isinstance(s, (list, tuple)):
        s = ",".join(s)
    return [c.strip() for c in str(s).split(",") if c.strip()]


def bh_fdr(p):
    """statsmodels.stats.multitest.fdrcorrection (Benjamini-Hochberg), NaN-free input."""
    np = _np()
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ps = p[order]
    raw = ps * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(raw[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.empty(n)
    out[order] = adj
    return out


def mannwhitney_legacy(x, y) -> float:
    """scipy<=1.6 mannwhitneyu(x, y) with alternative=None (the pinned scipy 1.5.2 default):
    one-sided p = norm.sf(|z|), z = (max(U1,U2) - n1 n2/2 - 0.5) / sd, tie-corrected."""
    np = _np()
    from scipy import stats  # type: ignore
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = x.size, y.size
    ranked = stats.rankdata(np.concatenate([x, y]))
    u1 = n1 * n2 + n1 * (n1 + 1) / 2.0 - ranked[:n1].sum()
    u2 = n1 * n2 - u1
    t = stats.tiecorrect(ranked)
    if t == 0:
        return float("nan")
    sd = math.sqrt(t * n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (max(u1, u2) - (n1 * n2 / 2.0 + 0.5)) / sd
    return float(stats.norm.sf(abs(z)))


def mannwhitney(x, y, mode: str = "legacy") -> float:
    from scipy import stats  # type: ignore
    if mode == "legacy":
        return mannwhitney_legacy(x, y)
    try:
        return float(stats.mannwhitneyu(x, y, alternative="two-sided")[1])
    except ValueError:
        return float("nan")


# ---------------------------------------------------------------------------
# Sample sheet (rules/common.smk)
# ---------------------------------------------------------------------------

def validate_sample_sheet(ss, key_cols, rep_cols) -> "list[str]":
    msgs = []
    for col in REQUIRED_COLS:
        if col not in ss.columns:
            raise ValueError("Missing required column in sample sheet: " + col)
    if not ss["SampleID"].is_unique:
        raise ValueError("SampleID column in samplesheet must not contain duplicates.")
    for col in key_cols:
        if col not in ss.columns:
            raise ValueError("Missing column in sample sheet that is provided in experiment_keycols in the config file: " + col)
    for col in rep_cols:
        if col not in ss.columns:
            raise ValueError("Missing column in sample sheet that is provided in replicate_keycols in the config file: " + col)
    if (ss["Bin"] == "Water").any():
        raise ValueError("Found 'Water' value in Bin column. Please instead enter a blank in this slot.")
    if key_cols:
        for _, row in ss[key_cols].drop_duplicates().iterrows():
            msgs.append("\t".join(str(v) for v in row.values))
    return msgs


def find_fastq_files(ss, fastqdir) -> Any:
    col = "fastqR1"
    if col in ss.columns:
        return ss
    ss = ss.copy()
    ss[col] = ""
    if not fastqdir:
        return ss
    for i in ss.index:
        s = ss.at[i, "SampleID"]
        files = glob.glob("{}_*_R1_*.fastq.gz".format(os.path.join(fastqdir, s)))
        if len(files) < 1:
            files = glob.glob("{}*.fastq.gz".format(os.path.join(fastqdir, s)))
        if len(files) > 1:
            raise ValueError("Found more than one FASTQ file for sample :" + s)
        if len(files) == 1:
            ss.at[i, col] = files[0]
    return ss


def add_experiment_names(ss, key_cols, rep_cols):
    for c in ("ExperimentIDReplicates", "ExperimentID", "ExperimentIDPCRRep"):
        if c in ss.columns:
            log.warning("ExperimentID columns found and will be overwritten in the sample sheet")
            ss = ss.drop(columns=[c])

    def j(row, cols):
        return "-".join(str(v) for v in list(row[cols])) if cols else ""

    s = ss[key_cols + rep_cols + ["PCRRep"]].drop_duplicates()
    s["ExperimentIDPCRRep"] = [j(r, key_cols) + "-Rep" + j(r, rep_cols) + "-PCR" + str(r["PCRRep"]) for _, r in s.iterrows()]
    ss = ss.merge(s, how="left")
    s = ss[key_cols + rep_cols].drop_duplicates() if (key_cols + rep_cols) else None
    if s is not None:
        s["ExperimentIDReplicates"] = [j(r, key_cols) + "-Rep" + j(r, rep_cols) for _, r in s.iterrows()]
        ss = ss.merge(s, how="left")
    else:
        ss["ExperimentIDReplicates"] = "-Rep"
    if key_cols:
        s = ss[key_cols].drop_duplicates()
        s["ExperimentID"] = ["-".join(str(v) for v in r.values) for _, r in s.iterrows()]
        ss = ss.merge(s, how="left")
    else:
        ss["ExperimentID"] = ""
    return ss


def add_outputs(ss):
    ss["ExperimentIDPCRRep_BinCounts"] = ["results/byPCRRep/{}.bin_counts.txt".format(e) for e in ss["ExperimentIDPCRRep"]]
    ss["ExperimentIDReplicates_BinCounts"] = ["results/byExperimentRep/{}.bin_counts.txt".format(e) for e in ss["ExperimentIDReplicates"]]
    return ss


def load_sample_sheet(path, key_cols, rep_cols, fastqdir=None):
    pd = _pd()
    ss = pd.read_csv(path, sep="\t", dtype=str)
    msgs = validate_sample_sheet(ss, key_cols, rep_cols)
    ss = find_fastq_files(ss.reset_index(drop=True), fastqdir)
    ss = add_experiment_names(ss, key_cols, rep_cols)
    ss = add_outputs(ss)
    ss.index = ss["SampleID"].values
    return ss, msgs


def get_bin_list(ss) -> "list[str]":
    b = ss["Bin"].drop_duplicates()
    b = b[(b != "All") & (b != "Neg") & b.notnull()]
    return [str(x) for x in b]


def unique_for_expt(ss, expt: str, col: str, bins) -> str:
    cur = ss.loc[(ss["ExperimentIDReplicates"] == expt) & ss["Bin"].isin(bins)]
    vals = cur[col].unique() if col in cur.columns else []
    if len(vals) != 1:
        raise ValueError(f"Found zero or more than one possible {col} for {expt}. Correct the samplesheet and rerun.")
    return str(vals[0])


def sortparams_file(ss, expt: str, bins, sortparams_dir) -> Path:
    batch = unique_for_expt(ss, expt, "Batch", bins)
    sample = unique_for_expt(ss, expt, "Sample", bins)
    return Path(sortparams_dir) / f"{batch}_{sample}.txt"


# ---------------------------------------------------------------------------
# Mapping + counting (Snakefile map_reads / count_reads / aggregate_alignment_stats)
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def iter_fastq(path):
    with _open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            yield h.rstrip("\n")[1:].split()[0], seq, q


def write_guide_fasta(design, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for _, r in design.iterrows():
            if isinstance(r.get("MappingSequence"), str) and r["MappingSequence"]:
                fh.write(f">{r['OligoID']}\n{r['MappingSequence']}\n")
    return out


def bowtie_command(fastq, index, sam_gz, unaligned) -> str:
    return (f"zless {fastq} | bowtie -v0 -p1 --all -x {index} --un {unaligned} - | gzip > {sam_gz}")


class ExactMapper:
    """bowtie -v0 --all semantics: every (reference, offset, strand) where the whole read matches exactly."""

    def __init__(self, design):
        self.refs = [(str(r["OligoID"]), str(r["MappingSequence"]).upper()) for _, r in design.iterrows()
                     if isinstance(r.get("MappingSequence"), str) and r["MappingSequence"]]
        self._idx = {}

    def _index(self, L: int):
        if L not in self._idx:
            d = {}
            for oid, seq in self.refs:
                for off in range(0, len(seq) - L + 1):
                    d.setdefault(seq[off:off + L], []).append((oid, off, "+"))
                rc = revcomp(seq)
                for off in range(0, len(seq) - L + 1):
                    d.setdefault(rc[off:off + L], []).append((oid, len(seq) - L - off, "-"))
            self._idx[L] = d
        return self._idx[L]

    def hits(self, seq: str):
        seq = seq.upper()
        if not seq:
            return []
        return self._index(len(seq)).get(seq, [])


def map_fastq_python(fastq, design, sam_gz: Path, unaligned_gz: Path) -> "tuple[int, int]":
    mapper = ExactMapper(design)
    n_al = n_un = 0
    sam_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(sam_gz, "wt") as out, gzip.open(unaligned_gz, "wt") as un:
        for name, seq, q in iter_fastq(fastq):
            hits = mapper.hits(seq)
            if not hits:
                n_un += 1
                un.write(f"@{name}\n{seq}\n+\n{q}\n")
                continue
            for oid, off, strand in hits:
                n_al += 1
                out.write(f"{name}\t{strand}\t{oid}\t{off}\t{seq}\t{q}\t{len(hits) - 1}\t\n")
    return n_al, n_un


def count_reads(aln_path):
    """cut -f3 | sort | uniq -c on a bowtie (default format) or SAM alignment file."""
    pd = _pd()
    c = Counter()
    with _open(aln_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[2] == "*":
                continue
            c[f[2]] += 1
    keys = sorted(c)
    return pd.DataFrame({"count": [c[k] for k in keys], "OligoID": keys})


def write_count_file(df, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df[["count", "OligoID"]].to_csv(path, sep="\t", index=False, header=False)
    return path


def map_sample(sample_id: str, fastq, out_dir: Path, design=None, index=None, force_python: bool = False) -> dict:
    """map_reads + count_reads for one sample; returns alignment stats row (or the command if no mapper)."""
    sam = out_dir / "samFiles" / f"{sample_id}.sam.gz"
    un_fq = out_dir / "samFiles" / f"{sample_id}.unaligned.fastq"
    stats_path = out_dir / "samFiles" / f"{sample_id}.alignment_stats.txt"
    sam.parent.mkdir(parents=True, exist_ok=True)
    cmd = bowtie_command(fastq, index or "<guide_index>", sam, un_fq)
    mapper = None
    if index and shutil.which("bowtie") and not force_python:
        mapper = "bowtie"
        subprocess.run(["bash", "-c", f"set -o pipefail; {cmd}"], check=True)
        n_al = sum(1 for _ in _open(sam))
        n_un = sum(1 for _ in open(un_fq)) // 4
        subprocess.run(["gzip", "-f", str(un_fq)], check=True)
    elif design is not None:
        mapper = "python-exact"
        n_al, n_un = map_fastq_python(fastq, design, sam, Path(str(un_fq) + ".gz"))
    else:
        return {"SampleID": sample_id, "mapper": None, "command": cmd}
    stats_path.write_text(f"SampleID\tnReadsAligned\tnReadsUnaligned\n{sample_id}\t{n_al}\t{n_un}\n")
    counts = count_reads(sam)
    cpath = write_count_file(counts, out_dir / "counts" / f"{sample_id}.count.txt")
    return {"SampleID": sample_id, "nReadsAligned": n_al, "nReadsUnaligned": n_un, "mapper": mapper,
            "count_file": str(cpath), "n_guides": int(len(counts))}


# ---------------------------------------------------------------------------
# Count tables (rules/make_count_tables.smk)
# ---------------------------------------------------------------------------

def read_count_file(path, name: str):
    pd, np = _pd(), _np()
    cur = pd.read_csv(path, sep="\t", names=[name, "OligoID"], dtype={"OligoID": str})
    cur[name] = cur[name].astype(np.int64)
    return cur


def make_count_table(ss, group_col: str, group_id: str, bins, counts_dir: Path):
    pd = _pd()
    cur = ss.loc[ss[group_col] == group_id]
    tbls = []
    for _, row in cur.iterrows():
        f = Path(counts_dir) / f"{row['SampleID']}.count.txt"
        if f.exists():
            tbls.append(read_count_file(f, row["SampleID"]))
    extra = sorted({str(b) for b in cur["Bin"].unique() if isinstance(b, str)} - set(bins))
    bin_list = list(bins) + extra
    if tbls:
        tbl = tbls.pop()
        for t in tbls:
            tbl = tbl.merge(t, on="OligoID", how="outer")
        tbl = tbl.set_index("OligoID").fillna(0).astype(int)
        for b in bin_list:
            samples = cur.loc[cur["Bin"] == b, "SampleID"]
            samples = [s for s in samples if s in tbl.columns]
            tbl[b] = tbl[samples].sum(axis=1).values if samples else 0
        tbl = tbl[bin_list]
    else:
        tbl = pd.DataFrame({"OligoID": []})
        for b in bins:
            tbl[b] = []
        tbl = tbl.set_index("OligoID")
    tbl = tbl.loc[:, tbl.sum(axis=0) != 0]
    tbl.index.name = "OligoID"
    freq = tbl.div(tbl.sum(axis=0), axis=1)
    return tbl, freq


def write_count_table(tbl, freq, counts_out: Path, freq_out: Path):
    counts_out.parent.mkdir(parents=True, exist_ok=True)
    tbl.to_csv(counts_out, sep="\t")
    freq.to_csv(freq_out, sep="\t", float_format="%.6f")
    emit("TSV", counts_out)
    emit("TSV", freq_out)


def make_flat_table(ss, counts_dir: Path):
    pd = _pd()
    parts = []
    for _, row in ss.iterrows():
        f = Path(counts_dir) / f"{row['SampleID']}.count.txt"
        if f.exists():
            cur = read_count_file(f, "count")
            cur["SampleID"] = row["SampleID"]
            parts.append(cur)
    if not parts:
        return pd.DataFrame(columns=["count", "OligoID", "SampleID"])
    return pd.concat(parts, axis="index", ignore_index=True)


def all_count_tables(ss, counts_dir: Path, results: Path) -> dict:
    bins = get_bin_list(ss)
    out = {"byPCRRep": [], "byExperimentRep": []}
    for col, sub in (("ExperimentIDPCRRep", "byPCRRep"), ("ExperimentIDReplicates", "byExperimentRep")):
        for e in ss[col].drop_duplicates():
            tbl, freq = make_count_table(ss, col, e, bins, counts_dir)
            write_count_table(tbl, freq, results / sub / f"{e}.bin_counts.txt", results / sub / f"{e}.bin_freq.txt")
            out[sub].append({"id": e, "n_guides": int(len(tbl)), "bins": list(tbl.columns), "reads": int(tbl.values.sum())})
    flat = make_flat_table(ss, counts_dir)
    write_tsv(flat, results / "summary" / "GuideCounts.flat.tsv.gz")
    out["flat_rows"] = int(len(flat))
    return out


def replicate_correlations(ss, base: Path, groupcol: str = "ExperimentIDReplicates",
                           filecol: str = "ExperimentIDPCRRep_BinCounts"):
    """plotReplicateCorrelations.R: per group and Bin, Pearson r among PCR-replicate count columns."""
    pd, np = _pd(), _np()
    rows = []
    for group in ss[groupcol].unique():
        cur = ss[ss[groupcol] == group]
        files = [f for f in cur[filecol].unique() if (Path(base) / f).exists()]
        data = {f: pd.read_csv(Path(base) / f, sep="\t", dtype={"OligoID": str}) for f in files}
        if len(data) < 2:
            continue
        for b in cur["Bin"].unique():
            parts = []
            for f, d in data.items():
                if b in d.columns:
                    parts.append(d[["OligoID", b]].rename(columns={b: os.path.basename(f)}))
            if len(parts) < 2:
                continue
            m = parts[0]
            for p in parts[1:]:
                m = m.merge(p, on="OligoID")
            X = m.drop(columns=["OligoID"])
            C = X.corr().to_numpy()
            names = list(X.columns)
            for j in range(len(names)):
                for i in range(j + 1, len(names)):
                    rows.append({"Bin": b, "R": C[i, j], "Group": group, "RepA": names[j], "RepB": names[i]})
    return pd.DataFrame(rows, columns=["Bin", "R", "Group", "RepA", "RepB"])


def plot_replicate_correlations(cor, path: Path):
    plt = _plt()
    if plt is None or cor.empty:
        return None
    bins = list(dict.fromkeys(cor["Bin"]))
    fig, ax = plt.subplots(figsize=(max(3, 0.6 * len(bins) + 1.5), 2.6))
    ax.boxplot([cor.loc[cor["Bin"] == b, "R"].dropna().values for b in bins], labels=bins)
    ax.set_ylim(0, 1)
    _style(ax, "PCR replicate correlations", "Bin", "Pearson R")
    return _save(fig, path)


# ---------------------------------------------------------------------------
# Sort params loaders (estimate_effect_sizes.R)
# ---------------------------------------------------------------------------

def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _factor_counts(labels_all, counts_all, selected_labels):
    """R 3.6: Count[factor] indexes by the integer level code (levels = sorted unique labels)."""
    vals = [str(v) for v in labels_all if isinstance(v, str) or not (isinstance(v, float) and math.isnan(v))]
    numeric = all(re.fullmatch(r"-?\d+(\.\d+)?", v) for v in vals)
    out = []
    if numeric:  # read.delim reads a numeric Bin column as numbers, so Count[Bin] is positional by value
        for lab in selected_labels:
            k = int(float(lab))
            out.append(float(counts_all[k - 1]) if 1 <= k <= len(counts_all) else float("nan"))
        return out
    levels = sorted(set(vals))
    for lab in selected_labels:
        k = levels.index(str(lab)) + 1
        out.append(float(counts_all[k - 1]) if k <= len(counts_all) else float("nan"))
    return out


def _std_dev_col(cols) -> Optional[str]:
    for c in ("Std.Dev.", "StdDev", "Std.Dev", "SD", "Std..Dev."):
        if c in cols:
            return c
    return None


def load_sort_params(path, bin_names, fmt: str = "bigfoot_nototal", upstream_compat: bool = False,
                     total_binname: str = "Total") -> dict:
    """Returns {'bins': DataFrame[name, mean, lowerBound, upperBound, count] (log10, file order),
    'totalCount', 'sd_seed', 'mu_seed', 'format'} like the R loaders."""
    pd, np = _pd(), _np()
    fmt = fmt.lower()
    if fmt not in SORT_FORMATS:
        raise ValueError(f"unknown sort params format {fmt}; one of {SORT_FORMATS}")
    sp = pd.read_csv(path, sep="," if fmt == "bigfoot" else "\t")
    sp.columns = [r_make_names(c) for c in sp.columns]
    bin_names = [str(b) for b in bin_names]
    if fmt in ("astrios", "astrios_nototal"):
        req = ["Mean", "Bounds", "Barcode", "Count"]
        if not all(c in sp.columns for c in req):
            raise ValueError(f"Sort parameters file did not have the needed : {req}")
        if total_binname not in set(sp["Barcode"].astype(str)):
            raise ValueError(f"Barcode column in sort parameters file needs an entry called '{total_binname}'")
        rows = sp[sp["Barcode"].astype(str) != total_binname]
        names = [str(b).replace("-", ".") for b in rows["Barcode"]]
        bounds = []
        for s in rows["Bounds"]:
            lo, hi = re.sub(r"[() ]", "", str(s)).split(",")
            bounds.append((float(lo), float(hi)))
        bounds = np.array(bounds)
        total = float(sp.loc[sp["Barcode"].astype(str) == total_binname, "Count"].iloc[0]) if fmt == "astrios" else 0.0
        sdc = _std_dev_col(sp.columns)
        seed = math.log10(1 + (float(sp.iloc[0][sdc]) / float(sp.iloc[0]["Mean"])) ** 2) if sdc else 0.5
        if sum(1 for b in sp["Barcode"].astype(str) if b in names) != len(names):
            raise ValueError("Sort parameters file did not have all the bins")
        bins = pd.DataFrame({"name": names, "mean": np.log10(rows["Mean"].astype(float).values),
                             "lowerBound": np.log10(bounds[:, 0]), "upperBound": np.log10(bounds[:, 1]),
                             "count": rows["Count"].astype(float).values}, index=names)
        mu_seed = float(np.log10(bounds[3, 0])) if len(bounds) >= 4 else float("nan")
        return {"bins": bins, "totalCount": total, "sd_seed": seed, "mu_seed": mu_seed, "format": fmt}
    if fmt == "bigfoot":
        req = ["Mean", "Min", "Max", "Barcode", "Count"]
        if not all(c in sp.columns for c in req):
            raise ValueError(f"Sort parameters file did not have the needed : {req}")
        lab = sp["Barcode"].astype(str)
        if total_binname not in set(lab):
            raise ValueError(f"Barcode column in sort parameters file needs an entry called '{total_binname}'")
        if sum(1 for b in lab if b in bin_names) != len(bin_names):
            raise ValueError("Sort parameters file did not have all the bins")
        sel = sp[lab.isin(bin_names)]
        tot = sp[lab == total_binname].iloc[0]
        sdc = _std_dev_col(sp.columns)
        seed = math.log10(1 + (_num(tot[sdc]) / _num(tot["Mean"])) ** 2) if sdc else 0.5
        names = list(sel["Barcode"].astype(str))
        bins = pd.DataFrame({"name": names, "mean": np.log10(sel["Mean"].map(_num).values),
                             "lowerBound": np.log10(sel["Min"].map(_num).values),
                             "upperBound": np.log10(sel["Max"].map(_num).values),
                             "count": sel["Count"].map(_num).values}, index=names)
        mu_seed = float("nan") if upstream_compat else (float(bins["lowerBound"].iloc[2]) if len(bins) >= 3 else float("nan"))
        return {"bins": bins, "totalCount": _num(tot["Count"]), "sd_seed": seed, "mu_seed": mu_seed, "format": fmt}
    # bigfoot_nototal: what estimate_effect_sizes.R calls
    req = ["Mean", "Min", "Max", "Bin", "Count"]
    if not all(c in sp.columns for c in req):
        raise ValueError(f"Sort parameters file did not have the needed : {req}")
    lab = sp["Bin"].astype(str)
    sel = sp[lab.isin(bin_names)]
    names = list(sel["Bin"].astype(str))
    if upstream_compat:
        counts = _factor_counts(list(sp["Bin"]), list(sp["Count"].astype(float)), names)
    else:
        counts = list(sel["Count"].astype(float))
    bins = pd.DataFrame({"name": names, "mean": np.log10(sel["Mean"].astype(float).values),
                         "lowerBound": np.log10(sel["Min"].astype(float).values),
                         "upperBound": np.log10(sel["Max"].astype(float).values),
                         "count": counts}, index=names)
    mu_seed = float(bins["lowerBound"].iloc[2]) if len(bins) >= 3 else float("nan")
    return {"bins": bins, "totalCount": 0.0, "sd_seed": 0.5, "mu_seed": mu_seed, "format": fmt}


# ---------------------------------------------------------------------------
# estimate_effect_sizes.R: weighted average + binned log-normal MLE
# ---------------------------------------------------------------------------

class MLEError(Exception):
    pass


def _pnorm(q, mu: float, sd: float):
    np = _np()
    from scipy.special import ndtr  # type: ignore
    if sd == 0:
        return (q >= mu).astype(float)
    return ndtr((q - mu) / sd)


def nll_binned(mu: float, sd: float, obs, lb, ub) -> float:
    """The R `ll`: bins 1..n from the gates, bin n+1 = 1 - sum; p <= 0 -> 1e-10; -sum(obs * log p)."""
    np = _np()
    if sd < 0:
        return 1e10
    pe = _pnorm(ub, mu, sd) - _pnorm(lb, mu, sd)
    pe = np.append(pe, 1.0 - pe.sum())
    pe[pe <= 0] = 1e-10
    if pe.size != obs.size:
        raise MLEError("Bin counts and bin boundaries don't have matching dimensions")
    with np.errstate(invalid="ignore"):
        return float(-(np.log(pe) * obs).sum())


def _r_grad(fn, p, lower, upper, ndeps=1e-3, bounded=True):
    """R optim's numerical gradient (fmingr): central differences, clipped to the box when bounded."""
    np = _np()
    g = np.zeros(p.size)
    for i in range(p.size):
        x = p.copy()
        eps_up = eps_dn = ndeps
        t = p[i] + ndeps
        if bounded and t > upper[i]:
            t = upper[i]
            eps_up = t - p[i]
        x[i] = t
        v1 = fn(x)
        t = p[i] - ndeps
        if bounded and t < lower[i]:
            t = lower[i]
            eps_dn = p[i] - t
        x[i] = t
        v2 = fn(x)
        den = eps_up + eps_dn
        g[i] = (v1 - v2) / den if den > 0 else 0.0
    return g


def r_mle(fn2, start, lower, upper):
    """stats4::mle(method='L-BFGS-B') == optim(..., hessian=TRUE) followed by solve(hessian)."""
    np = _np()
    from scipy.optimize import minimize  # type: ignore
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)

    def fn(x):
        v = fn2(float(x[0]), float(x[1]))
        if not math.isfinite(v):
            raise MLEError("L-BFGS-B needs finite values of 'fn'")
        return v

    x0 = np.clip(np.asarray(start, float), lower, upper)
    res = minimize(fn, x0, jac=lambda x: _r_grad(fn, x, lower, upper), method="L-BFGS-B",
                   bounds=list(zip(lower, upper)),
                   options={"maxcor": 5, "ftol": 1e7 * np.finfo(float).eps, "gtol": 0.0, "maxiter": 100})
    par = np.asarray(res.x, float)
    # optimhess: unbounded central differences of the numerical gradient
    H = np.zeros((2, 2))
    for i in range(2):
        a = par.copy(); a[i] += 1e-3
        b = par.copy(); b[i] -= 1e-3
        H[i] = (_r_grad(fn, a, lower, upper, bounded=False) - _r_grad(fn, b, lower, upper, bounded=False)) / 2e-3
    H = 0.5 * (H + H.T)
    if not np.all(np.isfinite(H)):
        raise MLEError("non-finite Hessian")
    with np.errstate(divide="ignore"):
        cond = np.linalg.cond(H, 1)
    if not math.isfinite(cond) or 1.0 / cond < np.finfo(float).eps:
        raise MLEError("system is computationally singular (Hessian)")
    return float(par[0]), float(par[1])


def normal_mle_guide(mu_seed: float, sd_seed: float, counts, lb, ub,
                     minmean: float = MINMEAN, maxmean: float = MAXMEAN) -> "tuple[float, float, str]":
    """getNormalMLE with the EM estimate of the seventh (ungated) bin."""
    np = _np()
    counts = np.asarray(counts, float)
    if not (math.isfinite(mu_seed) and math.isfinite(sd_seed)):
        raise MLEError("missing seed")
    o = np.append(counts, 0.0)
    MU, SI = r_mle(lambda m, s: nll_binned(m, s, o, lb, ub), [mu_seed, sd_seed],
                   [math.log10(minmean), 0.0], [math.log10(maxmean), 1.0])
    ps = _pnorm(ub, MU, SI) - _pnorm(lb, MU, SI)
    with np.errstate(divide="ignore", invalid="ignore"):
        seventh = counts.sum() / ps.sum() * (1 - ps.sum())
    o = np.append(counts, seventh)
    MU, SI = r_mle(lambda m, s: nll_binned(m, s, o, lb, ub), [mu_seed, sd_seed],
                   [math.log10(minmean), 0.1], [math.log10(maxmean), 1.0])
    return MU, SI, "EM"


def estimate_effect_sizes(counts, sp: dict, upstream_compat: bool = False, mu_seed=None, sd_seed=None):
    """estimate_effect_sizes.R on a bin_counts table; returns (raw_effects DataFrame, info dict)."""
    pd, np = _pd(), _np()
    if len(counts) == 0:
        raise ValueError("Unable to merge design doc and counts file (empty counts table)")
    mS = counts.copy()
    mS.columns = ["OligoID" if c == "OligoID" else r_make_names(c) for c in mS.columns]
    bin_names = [c for c in mS.columns if c not in ("OligoID", "All")]
    bins = sp["bins"]
    raw = mS[bin_names].astype(float)
    mS["sumreads"] = raw.sum(axis=1)
    mS["numobsbins"] = (raw > 0).sum(axis=1)
    for b in bin_names:
        cnt = float(bins.loc[b, "count"]) if b in bins.index else float("nan")
        with np.errstate(divide="ignore", invalid="ignore"):
            mS[b] = raw[b] / raw[b].sum() * cnt
    num = mS.select_dtypes(include=[np.number]).columns
    mS[num] = mS[num].fillna(0)
    # weighted average (bin mid-points; A at its upper gate, F at its lower gate)
    mid = 0.5 * (bins["lowerBound"] + bins["upperBound"])
    if "A" in mid.index:
        mid["A"] = bins.loc["A", "upperBound"]
    if "F" in mid.index:
        mid["F"] = bins.loc["F", "lowerBound"]
    mS["sumcells"] = mS[bin_names].sum(axis=1)
    missing = [b for b in bin_names if b not in bins.index]
    wa_bins = bin_names if upstream_compat else [b for b in bin_names if b in bins.index]
    wmean = np.array([mid[b] if b in mid.index else np.nan for b in wa_bins])
    with np.errstate(invalid="ignore", divide="ignore"):
        wsum = (mS[wa_bins].to_numpy() * wmean).sum(axis=1)
        mS["WeightedAvg"] = np.where(mS["sumcells"].to_numpy() == 0, 0.0, wsum / mS["sumcells"].to_numpy())
    if "All" in mS.columns:
        mS["input.fraction"] = mS["All"] / mS["All"].sum()
    mu0 = sp["mu_seed"] if mu_seed is None else float(mu_seed)
    sd0 = sp["sd_seed"] if sd_seed is None else float(sd_seed)
    if upstream_compat:
        use = bin_names
        lb = bins["lowerBound"].to_numpy(float)
        ub = bins["upperBound"].to_numpy(float)
    else:
        use = [b for b in bin_names if b in bins.index]
        lb = bins.loc[use, "lowerBound"].to_numpy(float)
        ub = bins.loc[use, "upperBound"].to_numpy(float)
    C = mS[use].to_numpy(float)
    means, sds, methods, errors = [], [], [], 0
    for i in range(len(mS)):
        try:
            m, s, meth = normal_mle_guide(mu0, sd0, C[i], lb, ub)
        except (MLEError, ValueError, ZeroDivisionError, FloatingPointError, np.linalg.LinAlgError) as e:
            errors += 1
            log.info("MLE errored out at given initialization for guide %d, using weighted average instead: %s", i + 1, e)
            m, s, meth = float(mS["WeightedAvg"].iloc[i]), sd0, "weighted_avg"
        means.append(m); sds.append(s); methods.append(meth)
    mS["logMean"] = means
    mS["logSD"] = sds
    mS["method"] = methods
    info = {"bins": bin_names, "bins_missing_sort_params": missing, "mu_seed": mu0, "sd_seed": sd0,
            "n_guides": int(len(mS)), "n_mle": int(sum(1 for x in methods if x == "EM")), "n_weighted_avg_fallback": errors,
            "sort_params_format": sp.get("format")}
    return mS, info


# ---------------------------------------------------------------------------
# convert_to_real_space.py
# ---------------------------------------------------------------------------

def _is_valid_guide(row, minlen, maxlen, minofftarget, maxg, minsum, minbins) -> bool:
    def test(seq):
        return (len(seq) >= minlen and len(seq) <= maxlen and row["OffTargetScore"] >= minofftarget
                and seq.count("G") <= maxg and row["sumcells"] >= minsum and row["numobsbins"] >= minbins)
    try:
        return bool(test(row["GuideSequence"]))
    except (TypeError, AttributeError, KeyError):
        try:
            return bool(test(row["GuideSequenceMinusG"]))
        except (TypeError, AttributeError, KeyError):
            return False


def convert_to_real_space(mle, design, clamp: float = 5.0, background: float = 0.0, maxg: int = 10, maxlen: int = 21,
                          minlen: int = 18, minofftarget: float = 50, minsum: float = 0.0, minbins: int = 0):
    np = _np()
    lines = []
    mle = mle.copy()
    mle["OligoID"] = mle["OligoID"].astype(str)
    design = design.copy()
    design["OligoID"] = design["OligoID"].astype(str)
    mle = mle.merge(design, on="OligoID")
    lines.append(f"Input guides: {len(mle)}")
    lines.append(f"Input negative control guides: {int((mle['target'] == NEG).sum())}")
    mle["logMean"] = mle["logMean"].astype(float)
    mle["logSD"] = mle["logSD"].astype(float)
    mle["mleAvg"] = np.power(10, mle["logMean"] + np.square(mle["logSD"]) / 2)
    mle["mleSD"] = np.sqrt((np.power(10, np.square(mle["logSD"])) - 1) * np.square(mle["mleAvg"]))
    mle["mleAvg"] = mle["mleAvg"] - background
    keep = mle.apply(lambda r: _is_valid_guide(r, minlen, maxlen, minofftarget, maxg, minsum, minbins), axis=1) \
        if len(mle) else mle.index == mle.index
    mle = mle.loc[keep].copy()
    nc = mle["target"] == NEG
    lines.append(f"Filtered guides: {len(mle)}")
    lines.append(f"Filtered negative control guides: {int(nc.sum())}")
    if nc.sum() == 0:
        raise ValueError("no negative-control guides passed the filters")
    nc_med = float(np.median(mle.loc[nc, "mleAvg"]))
    lines.append(f"Negative control mean: {nc_med}")
    lines.append(f"All guides mean: {mle['mleAvg'].mean()}")
    mle["mleAvg"] = mle["mleAvg"] / nc_med
    mle["WeightedAvg"] = mle["WeightedAvg"] / float(np.median(mle.loc[nc, "WeightedAvg"]))
    lines.append(f"Negative control mean after scale: {mle.loc[nc, 'mleAvg'].mean()}")
    mle.loc[mle["mleAvg"] > clamp, "mleAvg"] = clamp
    mle.loc[mle["mleAvg"] < 0, "mleAvg"] = 0.0
    nc2 = float(mle.loc[nc, "mleAvg"].mean())
    lines.append(f"Negative control mean after clip: {nc2}")
    mle["mleAvg"] = mle["mleAvg"] / nc2
    lines.append(f"Negative control mean after final rescale: {mle.loc[nc, 'mleAvg'].mean()}")
    lines.append(f"All guides mean after final rescale: {mle['mleAvg'].mean()}")
    for c in BED_COLS:
        if c not in mle.columns:
            mle[c] = np.nan
    info = {"input_guides": None, "negative_control_median_raw": nc_med, "nc_mean_after_clip": nc2,
            "guides_passing": int(len(mle)), "nc_passing": int(nc.sum()), "log": lines}
    return mle[BED_COLS].reset_index(drop=True), info


def write_bedgraph(df, path: Path, score: str) -> Path:
    np = _np()
    b = df.loc[~df["chr"].isna(), ["chr", "start", "end", score]].copy()
    b["start"] = b["start"].astype(float).astype(np.int64)
    b["end"] = b["end"].astype(float).astype(np.int64)
    b.to_csv(path, sep="\t", index=False, header=False)
    return emit("Wrote", path)


# ---------------------------------------------------------------------------
# CalculateTilingStatistic.R (sliding guide windows)
# ---------------------------------------------------------------------------

def sliding_window_statistic(real, window: int = 10, max_span: float = 750, score_col: str = "mleAvg",
                             utest_mode: str = "two-sided"):
    """Windows of `window` consecutive non-control guides (by chr, start) spanning <= max_span bp;
    mean = window mean - control mean (so 20% remaining at a control mean of 1 is -0.8)."""
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    ctrl = real.loc[real["target"] == NEG, score_col].astype(float).dropna().to_numpy()
    if ctrl.size == 0:
        raise ValueError("no negative control guides for the window statistic")
    g = real[(real["target"] != NEG) & real["chr"].notna() & real["start"].notna()].copy()
    g["start"] = g["start"].astype(float)
    g["end"] = g["end"].astype(float)
    g = g.sort_values(["chr", "start", "end"], kind="mergesort")
    cm = float(ctrl.mean())
    rows = []
    for chrom, sub in g.groupby("chr", sort=True):
        s = sub["start"].to_numpy(); e = sub["end"].to_numpy(); v = sub[score_col].astype(float).to_numpy()
        ids = sub["OligoID"].astype(str).to_numpy() if "OligoID" in sub.columns else np.arange(len(sub)).astype(str)
        for i in range(0, len(sub) - window + 1):
            lo, hi = s[i:i + window].min(), e[i:i + window].max()
            if hi - lo > max_span:
                continue
            x = v[i:i + window]
            p_t = float(stats.ttest_ind(x, ctrl, equal_var=False)[1]) if window > 1 else float("nan")
            rows.append({"chr": chrom, "start": int(lo), "end": int(hi), "mean": float(x.mean() - cm),
                         "sd": float(x.std(ddof=1)) if window > 1 else float("nan"),
                         "sem": float(x.std(ddof=1) / math.sqrt(window)) if window > 1 else float("nan"),
                         "n": int(window), "span": int(hi - lo), "firstGuide": ids[i], "lastGuide": ids[i + window - 1],
                         "p.utest": mannwhitney(x, ctrl, utest_mode), "p.ttest": p_t})
    cols = ["chr", "start", "end", "mean", "sd", "sem", "n", "span", "firstGuide", "lastGuide", "p.utest", "p.ttest"]
    res = pd.DataFrame(rows, columns=cols)
    for c in ("utest", "ttest"):
        p = res[f"p.{c}"].to_numpy(float)
        f = np.full(p.size, np.nan)
        ok = np.isfinite(p)
        f[ok] = bh_fdr(p[ok])
        res[f"fdr.{c}"] = f
    return res


# ---------------------------------------------------------------------------
# FlowFISHtssKD.py + normalize_flowfish_to_qpcr.py
# ---------------------------------------------------------------------------

def flowfish_tss_kd(gene: str, windows, qpcr, genes, slop: int = 500):
    pd, np = _pd(), _np()
    if "-" in gene and "Screen" in qpcr.columns:
        q = qpcr[qpcr["Screen"] == gene]
    else:
        q = qpcr.loc[qpcr["qPCRGene"].astype(str) == gene]
        q = q.iloc[[0]] if len(q) else q
    if len(q) != 1:
        raise ValueError("Screen (e.g. GATA1-1) must appear exactly once in qPCR file: " + gene)
    q = q.reset_index(drop=True).copy()
    gcols = [c for c in ["name", "tss", "chr", "start", "end", "strand"] if c in genes.columns]
    data = pd.merge(q.drop(columns=[c for c in ["tss", "chr", "start", "end", "strand"] if c in q.columns]),
                    genes[gcols], on="name")
    if data.empty:
        raise ValueError(f"qPCR 'name' {q['name'].iloc[0]} not found in the gene list")
    chrom = data["chr"].iloc[0]
    t_start, t_end = float(data["tss"].iloc[0]) - slop, float(data["tss"].iloc[0]) + slop
    override = False
    if "TSS_Override" in q.columns and not q["TSS_Override"].isnull().any():
        override = True
        t_start = int(float(q["TSS_Override"].iloc[0]) - slop)
        t_end = int(float(q["TSS_Override"].iloc[0]) + slop)
    prom = windows[(windows["chr"] == chrom) & (windows["start"] < t_end) & (windows["end"] > t_start)]
    if prom.empty:
        raise ValueError(f"no guide windows overlap the TSS +/- {slop} bp of {gene}")
    best = prom.sort_values(by="mean", ascending=True, kind="mergesort").iloc[0]
    q["FlowFISH_at_TSS"] = 1 + float(best["mean"])
    qp = pd.to_numeric(q["TSS_qPCR"], errors="coerce")
    q["TSS_qPCR"] = qp
    q["Background"] = (q["FlowFISH_at_TSS"] - qp) / (1 - q["FlowFISH_at_TSS"])
    q["chr"] = chrom
    q["start"] = int(t_start)
    q["end"] = int(t_end)
    return q, {"n_windows_at_tss": int(len(prom)), "best_window": f"{best['chr']}:{int(best['start'])}-{int(best['end'])}",
               "tss_override": override}


def flowfish_to_qpcr(info, data, clamp: float = 5.0, columns=("mleAvg", "mleSD")):
    np = _np()
    si = info.iloc[0]
    lo_ff = float(si["FlowFISH_at_TSS"])
    lo_q = _num(si.get("TSS_qPCR"))
    adjusted = True
    if math.isnan(lo_q):
        lo_q = lo_ff
        adjusted = False
    if math.isnan(lo_ff):  # no TSS window at all: identity transform
        lo_ff = lo_q = 0.0
        adjusted = False
    if not (0 <= lo_q < 1):
        raise ValueError(f"TSS_qPCR must be in [0, 1): {lo_q}")
    slope = (1 - lo_q) / (1 - lo_ff)
    out = data.copy()
    for c in columns:
        v = out[c].astype(float)
        if c.startswith("mleAvg"):
            v = (slope * (v - lo_ff) + lo_q).clip(lower=0, upper=clamp)
        elif c.startswith("mleSD"):
            v = slope * v
        else:
            raise ValueError(f"Don't know how to normalize column {c}")
        out[c] = v
    out["start"] = ["" if (x is None or (isinstance(x, float) and np.isnan(x)) or x == "") else str(int(float(x))) for x in out["start"]]
    out["end"] = ["" if (x is None or (isinstance(x, float) and np.isnan(x)) or x == "") else str(int(float(x))) for x in out["end"]]
    return out, {"slope": slope, "FlowFISH_at_TSS": lo_ff, "TSS_qPCR": lo_q, "adjusted": adjusted}


def no_qpcr_screen_info(gene: str, windows, genes=None, slop: int = 500):
    """ScreenInfo when no qPCR table is given: TSS_qPCR = NA (so normalize is the identity)."""
    pd = _pd()
    q = pd.DataFrame({"qPCRGene": [gene], "name": [gene], "TSS_qPCR": [float("nan")]})
    if genes is not None:
        try:
            res, meta = flowfish_tss_kd(gene, windows, q, genes, slop)
            return res, meta
        except ValueError as e:
            log.info("no TSS window: %s", e)
    q["FlowFISH_at_TSS"] = float("nan")
    q["Background"] = float("nan")
    return q, {"n_windows_at_tss": 0}


# ---------------------------------------------------------------------------
# score_peaks.smk: collapse + ScoreEnhancers.py + power + formatting
# ---------------------------------------------------------------------------

_NEG_WORD = re.compile(r"(?<![A-Za-z0-9_])negative_control(?![A-Za-z0-9_])")


def read_enhancers(path):
    """make_enhancers: cut -f1-3 | bedtools sort."""
    pd = _pd()
    e = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"], comment="#",
                    dtype={0: str})
    e = e[pd.to_numeric(e["start"], errors="coerce").notna()].copy()
    e["start"] = e["start"].astype(int); e["end"] = e["end"].astype(int)
    return e.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)


def collapse_to_peaks(scaled, enhancers):
    """tail -n+2 scaled | grep -vw negative_control | awk '$2!="NA"&&$2!="NaN"' | bedtools map -o collapse -c 12."""
    pd, np = _pd(), _np()
    from eqtl_enrichment_skill import overlap_join  # type: ignore
    s = scaled.copy()
    line = s.astype(str).agg("\t".join, axis=1)
    s = s[~line.str.contains(_NEG_WORD)]
    st = s["start"].astype(str)
    s = s[(st != "") & (st != "NA") & (st != "NaN") & (st != "nan") & s["chr"].notna()].copy()
    s["start"] = s["start"].astype(float).astype(int)
    s["end"] = s["end"].astype(float).astype(int)
    s = s.sort_values(["chr", "start", "end"], kind="mergesort").reset_index(drop=True)
    s["_ord"] = np.arange(len(s))
    b = s[["chr", "start", "end", "mleAvg", "_ord"]].rename(columns={"chr": "gchr", "start": "gstart", "end": "gend"})
    e = enhancers.reset_index(drop=True).copy()
    e["_eid"] = np.arange(len(e))
    j = overlap_join(e, b, a_cols=("chr", "start", "end"), b_cols=("gchr", "gstart", "gend"))
    rows = []
    if len(j):
        j = j.sort_values(["_eid", "_ord"])
        for eid, sub in j.groupby("_eid", sort=True):
            r = e.iloc[int(eid)]
            rows.append({"chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]),
                         "scores": ",".join(repr(float(v)) for v in sub["mleAvg"])})
    return pd.DataFrame(rows, columns=["chr", "start", "end", "scores"])


def score_enhancers(collapsed, scaled, min_guides: int = 3, fdr: float = 0.05, min_effect: float = 0.0,
                    expt_name: str = "", utest_mode: str = "legacy"):
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    data = collapsed.copy()
    data["scores"] = data["scores"].apply(lambda y: [float(x) for x in str(y).split(",")])
    nc = scaled[scaled["target"] == NEG]["mleAvg"].astype(float).to_numpy()
    data["n"] = data["scores"].apply(len)
    data = data[data["n"] >= min_guides].copy()
    data["mean"] = data["scores"].apply(lambda x: float(np.mean(x)))
    data["p.utest"] = data["scores"].apply(lambda x: mannwhitney(x, nc, utest_mode))
    data["fdr.utest"] = bh_fdr(data["p.utest"].to_numpy(float)) if len(data) else []
    data["p.ttest"] = data["scores"].apply(lambda x: float(stats.ttest_ind(x, nc, equal_var=True)[1]))
    data["fdr.ttest"] = bh_fdr(data["p.ttest"].to_numpy(float)) if len(data) else []
    data["sem"] = data["scores"].apply(lambda x: float(stats.sem(x)))
    nc_mean = float(np.mean(nc)) if nc.size else float("nan")
    es = (data["mean"] - nc_mean).abs()
    sigu = (data["fdr.utest"] < fdr) & (es >= min_effect)
    data["Regulated"] = sigu & (data["mean"] < nc_mean)
    data["Significant"] = (data["fdr.ttest"] < fdr) & (es >= min_effect)
    data["mean.ctrl"] = nc_mean
    data["n.ctrl"] = int(nc.size)
    data["sem.ctrl"] = float(stats.sem(nc)) if nc.size > 1 else float("nan")
    data["CustomTargetGeneTSS"] = ""
    tw = data[data["n"] >= min_guides]
    full = tw[SCORE_COLS].reset_index(drop=True)
    peaks = tw[tw["Significant"]][["chr", "start", "end"]].reset_index(drop=True)
    n = tw["n"]
    summary = pd.DataFrame({"Screen": [expt_name], "GuidePerDHS.min": [n.min() if len(n) else np.nan],
                            "GuidePerDHS.max": [n.max() if len(n) else np.nan],
                            "GuidePerDHS.median": [float(np.median(n)) if len(n) else np.nan],
                            "GuidePerDHS.MINGUIDES": [min_guides], "mean.ctrl": [nc_mean], "n.ctrl": [int(nc.size)],
                            "sem.ctrl": [float(stats.sem(nc)) if nc.size > 1 else np.nan],
                            "sig.elements.t": [int(tw["Significant"].sum())], "all.elements": [int(len(tw))]})
    return full, peaks, summary, data


def power_analysis(data, scaled, background: float = 0.0, fdr: float = 0.05, reps: int = 400,
                   effect_sizes=(0.05, 0.10, 0.15, 0.2, 0.25, 0.3), seed: int = 1):
    """The POWER block of ScoreEnhancers.py (indentation bug fixed so every rep is tested)."""
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    rng = np.random.default_rng(seed)
    nc = scaled[scaled["target"] == NEG]["mleAvg"].astype(float).to_numpy()
    bg = 0.0 if (background is None or not math.isfinite(float(background))) else float(background)
    n_fake = int(np.median(data["n"])) if len(data) else 0
    real_p = [float(stats.ttest_ind(x, nc, equal_var=True)[1]) for x in data["scores"]]
    out = []
    for es in effect_sizes:
        es_adj = 1 - (bg + 1 - es) / (bg + 1)
        shift = nc - es_adj
        hits = 0
        if n_fake >= 2 and n_fake <= nc.size:
            for _ in range(reps):
                sample = rng.choice(shift, size=n_fake, replace=False)
                p0 = float(stats.ttest_ind(sample, nc, equal_var=True)[1])
                f = bh_fdr(np.array([p0] + real_p))
                hits += int(f[0] < fdr)
        out.append({"EffectSize": es, "EffectSizeAdjusted": es_adj, "Power": hits / reps if n_fake >= 2 else float("nan"),
                    "nGuides": n_fake})
    return pd.DataFrame(out)


def screen_data(gene: str, score_path) -> str:
    return "Screen\tScreenData\tRNAReadoutMethod\tReference\n" + f"{gene}\t{score_path}\tFlowFISH Screen\tThisStudy\n"


def known_enhancers(full, gene: str, genes=None, cell_line: str = "", fdr: float = 0.05, screen: str = "",
                    reference: str = "ThisStudy"):
    """Reconstruction of FormatCRISPRiScreensForPredictions.R's KnownEnhancers columns."""
    pd, np = _pd(), _np()
    tss = chr_tss = None
    if genes is not None and len(genes):
        g = genes[genes["name"].astype(str) == gene]
        if g.empty and "symbol" in genes.columns:
            g = genes[genes["symbol"].astype(str) == gene]
        if len(g):
            tss = float(g["tss"].iloc[0]); chr_tss = g["chr"].iloc[0]
    k = pd.DataFrame({"chr": full["chr"], "start": full["start"].astype(int), "end": full["end"].astype(int)})
    k["name"] = [f"{gene}|{c}:{s}-{e}" for c, s, e in zip(k["chr"], k["start"], k["end"])]
    k["Gene"] = gene
    k["GeneTSS"] = tss if tss is not None else np.nan
    mid = (k["start"] + k["end"]) / 2.0
    k["distance"] = (mid - tss).abs() if tss is not None else np.nan
    if tss is not None:
        k.loc[k["chr"] != chr_tss, "distance"] = np.nan
    mc = full["mean.ctrl"].astype(float)
    k["Fraction.change.in.gene.expr"] = full["mean"].astype(float) / mc - 1
    k["p.value"] = full["p.ttest"]
    k["Adjusted.p.value"] = full["fdr.ttest"]
    k["Significant"] = full["fdr.ttest"].astype(float) < fdr
    k["Regulated"] = full["Regulated"].astype(bool) & k["Significant"]
    k["nGuides"] = full["n"]
    k["CellType"] = cell_line
    k["Screen"] = screen
    k["RNAReadoutMethod"] = "FlowFISH Screen"
    k["Reference"] = reference
    return k


# ---------------------------------------------------------------------------
# PlotGuideCounts.R (guide-count version)
# ---------------------------------------------------------------------------

def guide_count_qc(flat, ss, bins):
    pd, np = _pd(), _np()
    f = flat.copy()
    f["%Reads"] = f["count"] / f.groupby("SampleID")["count"].transform("sum") * 100
    meta = ss[["SampleID", "ExperimentIDReplicates", "Bin"]].copy()
    meta["Grouping"] = meta["ExperimentIDReplicates"].astype(str) + meta["Bin"].astype(str)
    cor_rows, cv_rows = [], []
    for grp, cur in meta.groupby("Grouping", sort=False):
        sub = f[f["SampleID"].isin(cur["SampleID"])]
        if sub.empty:
            continue
        M = sub.pivot_table(index="OligoID", columns="SampleID", values="%Reads", aggfunc="sum", fill_value=0)
        if M.shape[1] > 1:
            C = np.corrcoef(M.to_numpy().T)
            for i in range(M.shape[1]):
                for j in range(i + 1, M.shape[1]):
                    cor_rows.append({"Grouping": grp, "row": M.columns[i], "column": M.columns[j], "cor": C[i, j]})
        if len(cur) >= 2 and M.shape[1] >= 2:
            X = M.to_numpy()
            mu = X.mean(axis=1); sd = X.std(axis=1, ddof=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                cv = sd / mu * 100
            cv_rows.append(pd.DataFrame({"Grouping": grp, "OligoID": M.index, "mean": mu, "sd": sd, "CV": cv}))
    cors = pd.DataFrame(cor_rows, columns=["Grouping", "row", "column", "cor"])
    cvs = pd.concat(cv_rows, ignore_index=True) if cv_rows else pd.DataFrame(columns=["Grouping", "OligoID", "mean", "sd", "CV"])
    fm = f.merge(ss[["SampleID", "ExperimentIDReplicates", "Bin"]], on="SampleID")
    fm = fm[fm["Bin"].isin(bins)]
    grp = fm.groupby(["ExperimentIDReplicates", "OligoID", "Bin"], as_index=False)["%Reads"].mean() \
        .rename(columns={"%Reads": "GroupedFrequency"})
    if len(grp):
        grp["GroupedFrequencyAvgBin"] = grp.groupby(["ExperimentIDReplicates", "OligoID"])["GroupedFrequency"].transform("mean")
        grp["GroupedFreqRelativeToAverageBin"] = grp["GroupedFrequency"] / grp["GroupedFrequencyAvgBin"]
    return cors, cvs, grp


def plot_guide_counts(cors, cvs, grouped, out: Path, prefix: str) -> "list[Path]":
    plt = _plt()
    figs = []
    if plt is None:
        return figs
    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    if len(cors):
        axes[0].hist(cors["cor"].dropna(), bins=30, color=C1)
    _style(axes[0], "PCR replicate correlation", "Pearson correlation", "pairs")
    if len(cvs):
        ok = cvs["mean"] > 0
        axes[1].scatter(cvs.loc[ok, "mean"], cvs.loc[ok, "CV"], s=4, alpha=0.5, color=C1)
        axes[1].set_xscale("log")
    _style(axes[1], "Guide CV across PCR reps", "Guide frequency (%)", "CV (%)")
    figs.append(_save(fig, out / f"{prefix}.replicateCorrelations.png"))
    if len(grouped):
        for e, sub in list(grouped.groupby("ExperimentIDReplicates"))[:4]:
            piv = sub.pivot_table(index="OligoID", columns="Bin", values="GroupedFreqRelativeToAverageBin")
            top = piv.loc[sub.groupby("OligoID")["GroupedFrequencyAvgBin"].first().sort_values(ascending=False).index[:30]]
            fig, ax = plt.subplots(figsize=(8, 3))
            im = ax.imshow(top.T.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=0, vmax=2)
            ax.set_yticks(range(top.shape[1])); ax.set_yticklabels(top.columns)
            ax.set_xticks([])
            fig.colorbar(im, ax=ax, label="freq / avg bin")
            _style(ax, f"{e}: bin frequencies of the 30 most abundant guides", "guide", "Bin")
            figs.append(_save(fig, out / f"{prefix}.binBarplots.{safe_label(e)}.png"))
    return figs


def plot_effects(scaled, full, out: Path, name: str) -> "list[Path]":
    plt = _plt()
    if plt is None:
        return []
    figs = []
    s = scaled.copy()
    s = s[s["target"] != NEG]
    s = s[s["start"].astype(str) != ""]
    if len(s):
        fig, ax = plt.subplots(figsize=(8, 2.8))
        x = s["start"].astype(float)
        ax.scatter(x, s["mleAvg"].astype(float), s=5, color=C3)
        if len(full):
            for _, r in full.iterrows():
                ax.hlines(r["mean"], r["start"], r["end"], color=C2 if r["Significant"] else C1, lw=3)
        ax.axhline(1, color=AXIS, lw=0.8, ls="--")
        _style(ax, f"{name}: guide effects (scaled) and peak means (orange = significant)", "position", "fraction remaining")
        figs.append(_save(fig, out / f"{name}.effects.png"))
    return figs


# ---------------------------------------------------------------------------
# Whole workflow
# ---------------------------------------------------------------------------

def process_experiment(expt: str, count_path: Path, sp_path: Path, design, results: Path, *, gene_q: str, gene_sym: str,
                       qpcr=None, genes=None, enhancers=None, cell_line: str = "", sp_format: str = "bigfoot_nototal",
                       upstream_compat: bool = False, window: int = 10, max_span: float = 750, minsum: float = 0.05,
                       minbins: int = 4, clamp: float = 5.0, min_guides: int = 3, fdr: float = 0.05,
                       min_effect: float = 0.0, utest_mode: str = "legacy", do_power: bool = False,
                       plots: bool = True) -> dict:
    pd = _pd()
    d = results / "byExperimentRep"
    d.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(count_path, sep="\t", dtype={"OligoID": str})
    bins_in = [c for c in counts.columns if c not in ("OligoID", "All")]
    sp = load_sort_params(sp_path, bins_in, sp_format, upstream_compat)
    raw, einfo = estimate_effect_sizes(counts, sp, upstream_compat)
    write_tsv(raw, d / f"{expt}.raw_effects.txt")
    (d / f"{expt}.mle_log.txt").write_text(json.dumps(_clean(einfo), default=_jsonable) + "\n")
    real, rinfo = convert_to_real_space(raw, design, clamp=clamp, minsum=minsum, minbins=minbins)
    write_tsv(real, d / f"{expt}.real_space.txt")
    write_bedgraph(real, d / f"{expt}.real_space.bedgraph", "mleAvg")
    write_bedgraph(real, d / f"{expt}.real_space.bedgraph.WeightedAvg.bedgraph", "WeightedAvg")
    (d / f"{expt}.convert.log").write_text("\n".join(rinfo["log"]) + "\n")
    win = sliding_window_statistic(real, window, max_span, "mleAvg")
    write_tsv(win, d / f"{expt}.windows.txt")
    wb = win[["chr", "start", "end"]].copy(); wb["score"] = -win["mean"]
    wb.to_csv(d / f"{expt}.windows.txt.mean.bedgraph", sep="\t", index=False, header=False)
    if qpcr is not None and genes is not None:
        info, tinfo = flowfish_tss_kd(gene_q, win, qpcr, genes)
    else:
        info, tinfo = no_qpcr_screen_info(gene_q, win, genes)
    write_tsv(info, d / f"{expt}.ScreenInfo.txt")
    scaled, ninfo = flowfish_to_qpcr(info, real, clamp)
    write_tsv(scaled, d / f"{expt}.scaled.txt")
    write_bedgraph(scaled[scaled["start"] != ""], d / f"{expt}.scaled.bedgraph", "mleAvg")
    out = {"ExperimentIDReplicates": expt, "estimate": einfo, "real_space": {k: v for k, v in rinfo.items() if k != "log"},
           "n_windows": int(len(win)), "tss": tinfo, "qpcr": ninfo}
    if enhancers is None:
        return out
    scaled_rt = pd.read_csv(d / f"{expt}.scaled.txt", sep="\t", keep_default_na=True)
    scaled_rt["start"] = scaled_rt["start"].map(lambda x: "" if pd.isna(x) else str(int(x)))
    col = collapse_to_peaks(scaled_rt, enhancers)
    col.to_csv(d / f"{expt}.collapse.bed", sep="\t", index=False, header=False)
    emit("Wrote", d / f"{expt}.collapse.bed")
    full, peaks, summ, data = score_enhancers(col, scaled_rt, min_guides, fdr, min_effect, expt, utest_mode)
    write_tsv(full, d / f"{expt}.FullEnhancerScore.txt")
    peaks.to_csv(d / f"{expt}.peaks.tfdr{fdr}.bed", sep="\t", index=False, header=False)
    emit("Wrote", d / f"{expt}.peaks.tfdr{fdr}.bed")
    write_tsv(summ, d / f"{expt}.PeakCallingSummary.txt")
    (d / f"{expt}.ScreenData.txt").write_text(screen_data(gene_sym, d / f"{expt}.FullEnhancerScore.txt"))
    emit("Wrote", d / f"{expt}.ScreenData.txt")
    ke = known_enhancers(full, gene_sym, genes, cell_line, fdr, expt)
    write_tsv(ke, d / f"{expt}.KnownEnhancers.FlowFISH.txt")
    if do_power:
        pw = power_analysis(data, scaled_rt, background=_num(info["Background"].iloc[0]) if "Background" in info else 0.0,
                            fdr=fdr, reps=100)
        write_tsv(pw, d / f"{expt}.BH-Power-fdr{fdr}.txt")
        out["power"] = pw.to_dict(orient="records")
    if plots:
        plot_effects(scaled_rt, full, d, expt)
    out.update({"n_peaks_tested": int(len(full)), "n_significant": int(full["Significant"].sum()),
                "peaks": full[["chr", "start", "end", "mean", "n", "fdr.ttest", "Significant", "Regulated"]].to_dict(orient="records")})
    return out


def cmd_run(args) -> int:
    pd = _pd()
    key_cols, rep_cols = parse_cols(args.experiment_keycols), parse_cols(args.replicate_keycols)
    out = run_dir(args.label)
    results = out / "results"
    ss, msgs = load_sample_sheet(args.sample_sheet, key_cols, rep_cols, args.fastq_dir)
    ss.to_csv(out / "SampleList.snakemake.tsv", sep="\t", index=False)
    emit("TSV", out / "SampleList.snakemake.tsv")
    bins = get_bin_list(ss)
    design = pd.read_csv(args.design, sep="\t", dtype={"OligoID": str})
    counts_dir = results / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    align = []
    commands = []
    for _, r in ss.iterrows():
        sid = r["SampleID"]
        pre = Path(args.counts_dir) / f"{sid}.count.txt" if args.counts_dir else None
        if pre is not None and pre.exists():
            shutil.copy(pre, counts_dir / f"{sid}.count.txt")
            continue
        fq = r.get("fastqR1", "")
        if not isinstance(fq, str) or not fq:
            log.warning("no FASTQ and no count file for %s", sid)
            continue
        st = map_sample(sid, fq, results, design=design, index=args.bowtie_index, force_python=args.python_mapper)
        if st.get("mapper") is None:
            commands.append(st["command"])
        else:
            align.append(st)
    if align:
        al = pd.DataFrame(align)[["SampleID", "nReadsAligned", "nReadsUnaligned"]]
        write_tsv(al, results / "summary" / "alignment_stats.tsv")
    ct = all_count_tables(ss, counts_dir, results)
    cor = replicate_correlations(ss, out)
    write_tsv(cor, results / "summary" / "ReplicateCorrelation.pdf.tsv")
    plots = not args.no_plots
    if plots:
        plot_replicate_correlations(cor, results / "summary" / "ReplicateCorrelation.png")
    flat = make_flat_table(ss, counts_dir)
    cors, cvs, grouped = guide_count_qc(flat, ss, bins)
    write_tsv(cors, results / "summary" / "GuideCounts.replicateCorrelations.tsv")
    write_tsv(cvs, results / "summary" / "GuideCounts.guideCV.tsv")
    write_tsv(grouped, results / "summary" / "GuideCounts.binFrequencies.tsv")
    if plots:
        plot_guide_counts(cors, cvs, grouped, results / "summary", "GuideCounts")
    qpcr = pd.read_csv(args.qpcr, sep="\t") if args.qpcr else None
    genes = pd.read_csv(args.genelist, sep="\t") if args.genelist else None
    enh = read_enhancers(args.enhancers) if args.enhancers else None
    if enh is not None:
        write_tsv(enh, out / "config" / "enhancers_from_neighborhoods.bed", header=False)
    experiments = []
    expts = ss.loc[ss["Bin"].isin(bins), "ExperimentIDReplicates"].unique()
    for e in expts:
        cpath = results / "byExperimentRep" / f"{e}.bin_counts.txt"
        spath = sortparams_file(ss, e, bins, args.sortparams_dir)
        res = process_experiment(e, cpath, spath, design, results, gene_q=unique_for_expt(ss, e, "qPCRGene", bins),
                                 gene_sym=unique_for_expt(ss, e, "GeneSymbol", bins), qpcr=qpcr, genes=genes,
                                 enhancers=enh, cell_line=args.cell_line, sp_format=args.sort_params_format,
                                 upstream_compat=args.upstream_compat, window=args.window, max_span=args.max_span,
                                 minsum=args.minsum, minbins=args.minbins, clamp=args.clamp, min_guides=args.min_guides,
                                 fdr=args.fdr, min_effect=args.min_effect, utest_mode=args.utest_mode,
                                 do_power=args.power, plots=plots)
        experiments.append(res)
    if enh is not None and experiments:
        pcs = pd.concat([pd.read_csv(results / "byExperimentRep" / f"{e['ExperimentIDReplicates']}.PeakCallingSummary.txt", sep="\t")
                         for e in experiments], ignore_index=True)
        write_tsv(pcs, results / "summary" / "PeakCallingSummary.tsv")
    else:
        pcs = pd.DataFrame()
    summary = {"samples": int(len(ss)), "bins": bins, "experiment_keys": msgs, "mapping_commands_not_run": commands,
               "alignment": align, "count_tables": ct, "replicate_correlation_median": float(cor["R"].median()) if len(cor) else None,
               "experiments": experiments}
    body = [f"Samples: {len(ss)}; bins: {', '.join(bins)}; experimental replicates: {len(expts)}."]
    if commands:
        body.append("bowtie is not on PATH and no design was given for the Python mapper; commands that would run:\n\n```\n"
                    + "\n".join(commands) + "\n```")
    if align:
        body.append("## Alignment\n\n" + md_table(["SampleID", "aligned", "unaligned", "mapper"],
                                                  [[a["SampleID"], a["nReadsAligned"], a["nReadsUnaligned"], a["mapper"]] for a in align]))
    if len(cor):
        body.append("## PCR replicate correlations\n\n" + md_table(["Bin", "median R", "n pairs"],
                    [[b, g["R"].median(), len(g)] for b, g in cor.groupby("Bin", sort=False)]))
    for e in experiments:
        body.append(f"## {e['ExperimentIDReplicates']}\n\nMLE for {e['estimate']['n_mle']} of {e['estimate']['n_guides']} guides "
                    f"({e['estimate']['n_weighted_avg_fallback']} weighted-average fallbacks); {e['real_space']['guides_passing']} "
                    f"guides pass filters; {e['n_windows']} windows; qPCR slope {_fmt(e['qpcr']['slope'])}.")
        if e.get("peaks"):
            body.append(md_table(["chr", "start", "end", "mean", "n", "fdr.ttest", "Significant", "Regulated"],
                                 [[p["chr"], p["start"], p["end"], p["mean"], p["n"], p["fdr.ttest"], p["Significant"], p["Regulated"]]
                                  for p in e["peaks"]]))
    if len(pcs):
        body.append("## PeakCallingSummary\n\n" + md_table(list(pcs.columns), pcs.values.tolist()))
    return finish(out, f"CRISPRi-FlowFISH pipeline: {args.label}", summary, "\n\n".join(body))


# ---------------------------------------------------------------------------
# Step subcommands
# ---------------------------------------------------------------------------

def _name(args, path) -> str:
    n = getattr(args, "name", None)
    if n:
        return n
    b = os.path.basename(str(path))
    for suf in (".bin_counts.txt", ".raw_effects.txt", ".real_space.txt", ".windows.txt", ".scaled.txt", ".collapse.bed",
                ".FullEnhancerScore.txt", ".txt", ".tsv"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b


def cmd_samplesheet(args) -> int:
    out = run_dir(args.label)
    ss, msgs = load_sample_sheet(args.sample_sheet, parse_cols(args.experiment_keycols), parse_cols(args.replicate_keycols),
                                 args.fastq_dir)
    write_tsv(ss, out / "SampleList.snakemake.tsv")
    bins = get_bin_list(ss)
    s = {"samples": int(len(ss)), "bins": bins, "experiments": msgs,
         "ExperimentIDPCRRep": sorted(ss["ExperimentIDPCRRep"].unique().tolist()),
         "ExperimentIDReplicates": sorted(ss["ExperimentIDReplicates"].unique().tolist()),
         "ExperimentID": sorted(ss["ExperimentID"].unique().tolist()),
         "samples_without_fastq": ss.loc[ss["fastqR1"].fillna("") == "", "SampleID"].tolist()}
    body = (f"Validated {len(ss)} samples; bins {', '.join(bins)}.\n\n"
            + md_table(["SampleID", "Bin", "PCRRep", "ExperimentIDPCRRep", "ExperimentIDReplicates", "fastqR1"],
                       ss[["SampleID", "Bin", "PCRRep", "ExperimentIDPCRRep", "ExperimentIDReplicates", "fastqR1"]].values.tolist()))
    return finish(out, "Sample sheet", s, body)


def cmd_build_index(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    design = pd.read_csv(args.design, sep="\t", dtype={"OligoID": str})
    fa = write_guide_fasta(design, out / "guides.fa")
    emit("Wrote", fa)
    cmd = f"bowtie-build {fa} {out / 'guides'}"
    ran = False
    if shutil.which("bowtie-build"):
        subprocess.run(cmd.split(), check=True)
        ran = True
    else:
        print(f"bowtie-build not on PATH; would run: {cmd}")
    return finish(out, "Guide index", {"fasta": str(fa), "n_guides": int(sum(1 for l in open(fa) if l.startswith(">"))),
                                       "command": cmd, "ran": ran},
                  f"Guide FASTA `{fa}`.\n\n```\n{cmd}\n```\n\n" + ("Index built." if ran else "bowtie-build not on PATH: not run."))


def cmd_map_reads(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    design = pd.read_csv(args.design, sep="\t", dtype={"OligoID": str}) if args.design else None
    if not (args.index and shutil.which("bowtie")) and design is None:
        cmd = bowtie_command(args.fastq, args.index or "<guide_index>", out / "samFiles" / f"{args.sample_id}.sam.gz",
                             out / "samFiles" / f"{args.sample_id}.unaligned.fastq")
        print(f"bowtie not on PATH (or no --index) and no --design for the Python mapper; would run:\n{cmd}")
        return finish(out, "Map reads", {"command": cmd, "ran": False}, f"Not run.\n\n```\n{cmd}\n```")
    st = map_sample(args.sample_id, args.fastq, out / "results", design=design, index=args.index,
                    force_python=args.python_mapper)
    emit("TSV", out / "results" / "samFiles" / f"{args.sample_id}.alignment_stats.txt")
    emit("Wrote", Path(st["count_file"]))
    return finish(out, f"Map reads: {args.sample_id}", st,
                  md_table(["SampleID", "aligned", "unaligned", "mapper", "guides"],
                           [[st["SampleID"], st["nReadsAligned"], st["nReadsUnaligned"], st["mapper"], st["n_guides"]]]))


def cmd_count_reads(args) -> int:
    out = run_dir(args.label)
    c = count_reads(args.alignments)
    p = write_count_file(c, out / f"{args.sample_id}.count.txt")
    emit("Wrote", p)
    return finish(out, f"Count reads: {args.sample_id}", {"n_guides": int(len(c)), "n_alignments": int(c["count"].sum())},
                  f"{int(c['count'].sum())} alignments over {len(c)} guides.")


def cmd_count_tables(args) -> int:
    out = run_dir(args.label)
    ss, _ = load_sample_sheet(args.sample_sheet, parse_cols(args.experiment_keycols), parse_cols(args.replicate_keycols), None)
    ss.to_csv(out / "SampleList.snakemake.tsv", sep="\t", index=False)
    ct = all_count_tables(ss, Path(args.counts_dir), out / "results")
    cor = replicate_correlations(ss, out)
    write_tsv(cor, out / "results" / "summary" / "ReplicateCorrelation.pdf.tsv")
    if not args.no_plots:
        plot_replicate_correlations(cor, out / "results" / "summary" / "ReplicateCorrelation.png")
    rows = [[k, x["id"], x["n_guides"], ", ".join(map(str, x["bins"])), x["reads"]] for k in ("byPCRRep", "byExperimentRep") for x in ct[k]]
    return finish(out, "Count tables", {"count_tables": ct, "replicate_correlations": cor.to_dict(orient="records")},
                  md_table(["level", "id", "guides", "bins", "reads"], rows))


def cmd_replicate_correlations(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    ss = pd.read_csv(args.sample_sheet, sep="\t", dtype=str)
    cor = replicate_correlations(ss, Path(args.base_dir), args.groupcol, args.filecol)
    write_tsv(cor, out / "ReplicateCorrelation.tsv")
    if not args.no_plots:
        plot_replicate_correlations(cor, out / "ReplicateCorrelation.png")
    return finish(out, "Replicate correlations", {"n_pairs": int(len(cor)), "median_R": float(cor["R"].median()) if len(cor) else None},
                  md_table(["Bin", "median R", "pairs"], [[b, g["R"].median(), len(g)] for b, g in cor.groupby("Bin", sort=False)]))


def cmd_estimate_effects(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    counts = pd.read_csv(args.counts, sep="\t", dtype={"OligoID": str})
    name = _name(args, args.counts)
    sp = load_sort_params(args.sort_params, [c for c in counts.columns if c not in ("OligoID", "All")],
                          args.sort_params_format, args.upstream_compat)
    raw, info = estimate_effect_sizes(counts, sp, args.upstream_compat, args.mu_seed, args.sd_seed)
    write_tsv(raw, out / f"{name}.raw_effects.txt")
    b = sp["bins"]
    body = (f"{info['n_mle']} of {info['n_guides']} guides by MLE; {info['n_weighted_avg_fallback']} weighted-average fallbacks. "
            f"Seeds mu={_fmt(info['mu_seed'])}, sd={_fmt(info['sd_seed'])}.\n\n"
            + md_table(["bin", "log10 mean", "lower", "upper", "sorted cells"], b[["name", "mean", "lowerBound", "upperBound", "count"]].values.tolist()))
    return finish(out, f"Effect sizes: {name}", info, body)


def cmd_real_space(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = _name(args, args.mle)
    mle = pd.read_csv(args.mle, sep="\t", dtype={"OligoID": str})
    design = pd.read_csv(args.design, sep="\t", dtype={"OligoID": str})
    real, info = convert_to_real_space(mle, design, args.clamp, args.background, args.maxg, args.maxlen, args.minlen,
                                       args.minofftarget, args.minsum, args.minbins)
    write_tsv(real, out / f"{name}.real_space.txt")
    write_bedgraph(real, out / f"{name}.real_space.bedgraph", "mleAvg")
    write_bedgraph(real, out / f"{name}.real_space.bedgraph.WeightedAvg.bedgraph", "WeightedAvg")
    (out / f"{name}.convert.log").write_text("\n".join(info["log"]) + "\n")
    return finish(out, f"Real space: {name}", {k: v for k, v in info.items()}, "\n".join(f"- {l}" for l in info["log"]))


def cmd_windows(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = _name(args, args.input)
    real = pd.read_csv(args.input, sep="\t", dtype={"OligoID": str})
    win = sliding_window_statistic(real, args.window, args.max_span, args.score_column)
    write_tsv(win, out / f"{name}.windows.txt")
    wb = win[["chr", "start", "end"]].copy(); wb["score"] = -win["mean"]
    wb.to_csv(out / f"{name}.windows.txt.mean.bedgraph", sep="\t", index=False, header=False)
    emit("Wrote", out / f"{name}.windows.txt.mean.bedgraph")
    top = win.sort_values("mean").head(10)
    return finish(out, f"Guide windows: {name}", {"n_windows": int(len(win)), "window": args.window, "max_span": args.max_span},
                  f"{len(win)} windows of {args.window} guides (span <= {args.max_span} bp).\n\n"
                  + md_table(["chr", "start", "end", "mean", "fdr.ttest"], top[["chr", "start", "end", "mean", "fdr.ttest"]].values.tolist()))


def cmd_tss_kd(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    win = pd.read_csv(args.windows, sep="\t")
    genes = pd.read_csv(args.genelist, sep="\t")
    if args.qpcr:
        info, meta = flowfish_tss_kd(args.gene, win, pd.read_csv(args.qpcr, sep="\t"), genes, args.slop)
    else:
        info, meta = no_qpcr_screen_info(args.gene, win, genes, args.slop)
    name = _name(args, args.windows)
    write_tsv(info, out / f"{name}.ScreenInfo.txt")
    return finish(out, f"TSS knockdown: {args.gene}", {**meta, "row": info.iloc[0].to_dict()},
                  md_table(list(info.columns), info.values.tolist()))


def cmd_normalize_qpcr(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = _name(args, args.input)
    info = pd.read_csv(args.info, sep="\t")
    data = pd.read_csv(args.input, sep="\t", dtype={"OligoID": str})
    scaled, meta = flowfish_to_qpcr(info, data, args.clamp)
    write_tsv(scaled, out / f"{name}.scaled.txt")
    write_bedgraph(scaled[scaled["start"] != ""], out / f"{name}.scaled.bedgraph", "mleAvg")
    return finish(out, f"qPCR normalisation: {name}", meta,
                  f"slope = (1 - qPCR)/(1 - FF) = {_fmt(meta['slope'])} (FF {_fmt(meta['FlowFISH_at_TSS'])}, qPCR {_fmt(meta['TSS_qPCR'])}).")


def _read_scaled(path):
    pd = _pd()
    s = pd.read_csv(path, sep="\t", dtype={"OligoID": str})
    s["start"] = s["start"].map(lambda x: "" if pd.isna(x) else str(int(x)))
    return s


def cmd_collapse(args) -> int:
    out = run_dir(args.label)
    name = _name(args, args.scaled)
    enh = read_enhancers(args.enhancers)
    col = collapse_to_peaks(_read_scaled(args.scaled), enh)
    p = out / f"{name}.collapse.bed"
    col.to_csv(p, sep="\t", index=False, header=False)
    emit("Wrote", p)
    return finish(out, f"Collapse to peaks: {name}", {"n_enhancers": int(len(enh)), "n_with_guides": int(len(col))},
                  f"{len(col)} of {len(enh)} elements overlap at least one guide.")


def cmd_score_enhancers(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = args.expt_name or _name(args, args.collapsed)
    col = pd.read_csv(args.collapsed, sep="\t", header=None, names=["chr", "start", "end", "scores"], dtype={"scores": str})
    scaled = _read_scaled(args.scaled)
    full, peaks, summ, data = score_enhancers(col, scaled, args.min_guides, args.fdr, args.min_effect, name, args.utest_mode)
    write_tsv(full, out / f"{name}.FullEnhancerScore.txt")
    peaks.to_csv(out / f"{name}.peaks.tfdr{args.fdr}.bed", sep="\t", index=False, header=False)
    emit("Wrote", out / f"{name}.peaks.tfdr{args.fdr}.bed")
    write_tsv(summ, out / f"{name}.PeakCallingSummary.txt")
    if not args.no_plots:
        plot_effects(scaled, full, out, name)
    return finish(out, f"Enhancer scores: {name}", {"n_tested": int(len(full)), "n_significant": int(full["Significant"].sum()),
                                                    "summary": summ.iloc[0].to_dict()},
                  md_table(["chr", "start", "end", "mean", "n", "p.ttest", "fdr.ttest", "Significant", "Regulated"],
                           full[["chr", "start", "end", "mean", "n", "p.ttest", "fdr.ttest", "Significant", "Regulated"]].values.tolist()))


def cmd_power(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = args.expt_name or _name(args, args.collapsed)
    col = pd.read_csv(args.collapsed, sep="\t", header=None, names=["chr", "start", "end", "scores"], dtype={"scores": str})
    scaled = _read_scaled(args.scaled)
    _, _, _, data = score_enhancers(col, scaled, args.min_guides, args.fdr, 0.0, name)
    bg = 0.0
    if args.info:
        inf = pd.read_csv(args.info, sep="\t")
        if len(inf) != 1:
            raise ValueError("Each Screen must be present exactly once")
        bg = _num(inf["Background"].iloc[0]) if "Background" in inf else 0.0
    pw = power_analysis(data, scaled, bg, args.fdr, args.reps, args.effect_sizes, args.seed)
    write_tsv(pw, out / f"{name}.BH-Power-fdr{args.fdr}.txt")
    return finish(out, f"Power: {name}", {"background": bg, "power": pw.to_dict(orient="records")},
                  md_table(list(pw.columns), pw.values.tolist()))


def cmd_format_screen(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    name = args.name or _name(args, args.score)
    (out / f"{name}.ScreenData.txt").write_text(screen_data(args.gene, os.path.abspath(args.score)))
    emit("Wrote", out / f"{name}.ScreenData.txt")
    full = pd.read_csv(args.score, sep="\t")
    genes = pd.read_csv(args.genelist, sep="\t") if args.genelist else None
    ke = known_enhancers(full, args.gene, genes, args.cell_line, args.fdr, name)
    write_tsv(ke, out / f"{name}.KnownEnhancers.FlowFISH.txt")
    return finish(out, f"KnownEnhancers: {name}", {"n_elements": int(len(ke)), "n_significant": int(ke["Significant"].sum())},
                  md_table(["name", "distance", "Fraction.change.in.gene.expr", "Adjusted.p.value", "Significant"],
                           ke[["name", "distance", "Fraction.change.in.gene.expr", "Adjusted.p.value", "Significant"]].values.tolist()))


def cmd_guide_count_plots(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    flat = pd.read_csv(args.guide_counts, sep="\t", dtype={"OligoID": str})
    ss = pd.read_csv(args.sample_sheet, sep="\t", dtype=str)
    bins = [b for b in ss["Bin"].dropna().unique() if b not in ("All", "Neg", "")]
    cors, cvs, grouped = guide_count_qc(flat, ss, bins)
    write_tsv(cors, out / "GuideCounts.replicateCorrelations.tsv")
    write_tsv(cvs, out / "GuideCounts.guideCV.tsv")
    write_tsv(grouped, out / "GuideCounts.binFrequencies.tsv")
    if not args.no_plots:
        plot_guide_counts(cors, cvs, grouped, out, "GuideCounts")
    return finish(out, "Guide count QC", {"n_pairs": int(len(cors)), "median_cor": float(cors["cor"].median()) if len(cors) else None,
                                          "median_cv": float(cvs["CV"].median()) if len(cvs) else None},
                  f"{len(cors)} PCR-replicate pairs, median r {_fmt(cors['cor'].median() if len(cors) else float('nan'))}.")


# ---------------------------------------------------------------------------
# Self-test: synthetic screen with planted enhancers
# ---------------------------------------------------------------------------

GATES = {"A": (1.5, 2.6), "B": (2.62, 2.8), "C": (2.82, 2.95), "D": (2.97, 3.1), "E": (3.12, 3.28), "F": (3.3, 4.5)}


def synthetic_screen(d: Path, seed: int = 11, depth: int = 8000) -> dict:
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    peaks = [("TSS", 99_750, 100_250, 0.3), ("E1", 130_000, 130_500, 0.5), ("E2", 160_000, 160_500, 0.72),
             ("N1", 180_000, 180_500, 1.0), ("N2", 200_000, 200_500, 1.0), ("N3", 220_000, 220_500, 1.0)]
    used = set()

    def seq20():
        while True:
            s = "".join(rng.choice(list("ACGT"), 20))
            if s not in used and s.count("G") <= 8:
                used.add(s); return s

    rows = []
    for pname, s0, e0, eff in peaks:
        for k in range(16):
            st = s0 + 10 + 29 * k
            sq = seq20()
            rows.append({"chr": "chr1", "start": st, "end": st + 20, "name": f"{pname}_g{k}", "score": 0, "strand": "+",
                         "GuideSequence": sq, "GuideSequenceMinusG": sq, "MappingSequence": sq, "OffTargetScore": 80,
                         "target": pname, "subpool": "tile", "OligoID": f"{pname}_g{k}", "_eff": eff})
    for k in range(2):  # chr2 element with too few guides
        sq = seq20()
        rows.append({"chr": "chr2", "start": 5000 + 30 * k, "end": 5020 + 30 * k, "name": f"FEW_g{k}", "score": 0, "strand": "+",
                     "GuideSequence": sq, "GuideSequenceMinusG": sq, "MappingSequence": sq, "OffTargetScore": 80,
                     "target": "FEW", "subpool": "tile", "OligoID": f"FEW_g{k}", "_eff": 0.3})
    for k, (flt, gs, ot) in enumerate([("lowOT", None, 20), ("lowOT2", None, 30), ("short", "ACGTACGTACGTACG", 80)]):
        sq = seq20()
        rows.append({"chr": "chr1", "start": 130_020 + 7 * k, "end": 130_040 + 7 * k, "name": f"BAD_{flt}", "score": 0, "strand": "+",
                     "GuideSequence": gs or sq, "GuideSequenceMinusG": sq, "MappingSequence": sq, "OffTargetScore": ot,
                     "target": "E1", "subpool": "tile", "OligoID": f"BAD_{flt}", "_eff": 3.0})
    for k in range(80):
        sq = seq20()
        rows.append({"chr": None, "start": None, "end": None, "name": f"NC_{k}", "score": 0, "strand": None,
                     "GuideSequence": (None if k == 0 else sq), "GuideSequenceMinusG": sq, "MappingSequence": sq,
                     "OffTargetScore": 80, "target": NEG, "subpool": "ctrl", "OligoID": f"NC_{k}", "_eff": 1.0})
    guides = pd.DataFrame(rows)
    design = guides.drop(columns=["_eff"])
    design_path = d / "design.txt"
    design.to_csv(design_path, sep="\t", index=False)
    eff = guides["_eff"].to_numpy()
    mu = 3.0 + np.log10(eff)
    names = list(GATES)
    lo = np.array([GATES[b][0] for b in names]); hi = np.array([GATES[b][1] for b in names])
    fq_dir = d / "fastq"; fq_dir.mkdir()
    sp_dir = d / "sortParams"; sp_dir.mkdir()
    ss_rows = []
    oligos = guides["OligoID"].to_numpy()
    seqs = guides["MappingSequence"].to_numpy()
    junk = "N" * 20

    def write_fastq(path, counts):
        with gzip.open(path, "wt", compresslevel=1) as fh:
            i = 0
            for g in np.flatnonzero(counts):
                for _ in range(int(counts[g])):
                    fh.write(f"@r{i}\n{seqs[g]}\n+\n{'I' * 20}\n"); i += 1
            for _ in range(int(0.02 * counts.sum())):
                fh.write(f"@r{i}\n{junk}\n+\n{'I' * 20}\n"); i += 1

    for rep in (1, 2):
        abund = rng.lognormal(0, 0.3, len(guides))
        ncell = rng.poisson(300 * abund)
        cells = np.zeros((len(guides), 6))
        tot_expr = np.zeros(6)
        for g in range(len(guides)):
            x = rng.normal(mu[g], 0.25, ncell[g])
            for b in range(6):
                m = (x >= lo[b]) & (x < hi[b])
                cells[g, b] = m.sum()
                tot_expr[b] += (10 ** x[m]).sum()
        sort_counts = cells.sum(axis=0)
        sp = pd.DataFrame({"Bin": names, "Mean": tot_expr / sort_counts, "Min": 10 ** lo, "Max": 10 ** hi, "Count": sort_counts.astype(int)})
        sp.to_csv(sp_dir / f"B1_S{rep}.txt", sep="\t", index=False)
        for b in range(6):
            for pcr in (1, 2):
                sid = f"R{rep}_{names[b]}_P{pcr}"
                reads = rng.multinomial(depth, cells[:, b] / cells[:, b].sum())
                write_fastq(fq_dir / f"{sid}_S1_R1_001.fastq.gz", reads)
                ss_rows.append({"SampleID": sid, "Bin": names[b], "PCRRep": str(pcr), "qPCRGene": "GENE1", "GeneSymbol": "GENE1",
                                "CellLine": "K562", "FlowFISHRep": str(rep), "Batch": "B1", "Sample": f"S{rep}"})
        sid = f"R{rep}_All_P1"
        write_fastq(fq_dir / f"{sid}_S1_R1_001.fastq.gz", rng.multinomial(depth, ncell / ncell.sum()))
        ss_rows.append({"SampleID": sid, "Bin": "All", "PCRRep": "1", "qPCRGene": "GENE1", "GeneSymbol": "GENE1",
                        "CellLine": "K562", "FlowFISHRep": str(rep), "Batch": "B1", "Sample": f"S{rep}"})
    ss = pd.DataFrame(ss_rows)
    ss_path = d / "SampleSheet.tsv"
    ss.to_csv(ss_path, sep="\t", index=False)
    genes = pd.DataFrame({"chr": ["chr1", "chr1"], "start": [100_000, 400_000], "end": [120_000, 410_000], "name": ["GENE1", "GENE2"],
                          "score": [0, 0], "strand": ["+", "-"], "symbol": ["GENE1", "GENE2"], "tss": [100_000, 410_000]})
    genes.to_csv(d / "GeneList.txt", sep="\t", index=False)
    pd.DataFrame({"qPCRGene": ["GENE1"], "name": ["GENE1"], "TSS_qPCR": [0.15]}).to_csv(d / "qPCR.txt", sep="\t", index=False)
    enh = [("chr1", p[1], p[2]) for p in peaks] + [("chr1", 300_000, 300_500), ("chr2", 4990, 5100)]
    with open(d / "EnhancerList.bed", "w") as fh:
        for c, s, e in reversed(enh):
            fh.write(f"{c}\t{s}\t{e}\tpeak\t0\n")
    return {"design": design_path, "sample_sheet": ss_path, "fastq": fq_dir, "sortparams": sp_dir, "genes": d / "GeneList.txt",
            "qpcr": d / "qPCR.txt", "enhancers": d / "EnhancerList.bed", "peaks": peaks, "guides": guides}


def cmd_selftest(args) -> int:
    import tempfile
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made = []

    def newest(label):
        p = sorted(OUT_ROOT.glob(f"*_{label}"))[-1]
        made.append(p)
        return p

    t0 = time.time()
    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            print("\nstatistics helpers")
            p = np.array([0.01, 0.04, 0.03, 0.2, 0.5])
            from statsmodels.stats.multitest import fdrcorrection  # type: ignore
            check(np.allclose(bh_fdr(p), fdrcorrection(p)[1]), "BH FDR equals statsmodels fdrcorrection")
            x, y = np.array([0.5, 0.6, 0.55, 0.7]), np.linspace(0.8, 1.2, 30)
            p2 = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")[1]
            check(abs(mannwhitney_legacy(x, y) - p2 / 2) < 1e-12, "legacy Mann-Whitney p = half the asymptotic two-sided p (scipy 1.5 default)")
            print("\nsample sheet + experiment keys")
            W = synthetic_screen(d)
            ss, msgs = load_sample_sheet(W["sample_sheet"], ["CellLine"], ["FlowFISHRep"], str(W["fastq"]))
            check(set(ss["ExperimentIDReplicates"]) == {"K562-Rep1", "K562-Rep2"}, "ExperimentIDReplicates = key-Rep<rep>")
            check("K562-Rep1-PCR2" in set(ss["ExperimentIDPCRRep"]) and set(ss["ExperimentID"]) == {"K562"}, "ExperimentIDPCRRep / ExperimentID")
            check(get_bin_list(ss) == list("ABCDEF"), "bin list excludes All")
            check(ss["fastqR1"].str.endswith("_S1_R1_001.fastq.gz").all(), "FASTQ discovery {SampleID}_*_R1_*.fastq.gz")
            bad = pd.read_csv(W["sample_sheet"], sep="\t", dtype=str)
            bad.loc[0, "Bin"] = "Water"
            bad.to_csv(d / "bad.tsv", sep="\t", index=False)
            try:
                load_sample_sheet(d / "bad.tsv", ["CellLine"], ["FlowFISHRep"]); check(False, "Water bin rejected")
            except ValueError:
                check(True, "Water bin rejected")
            bad = pd.read_csv(W["sample_sheet"], sep="\t", dtype=str).drop(columns=["qPCRGene"])
            bad.to_csv(d / "bad2.tsv", sep="\t", index=False)
            try:
                load_sample_sheet(d / "bad2.tsv", ["CellLine"], ["FlowFISHRep"]); check(False, "missing required column rejected")
            except ValueError:
                check(True, "missing required column rejected")

            print("\nmapping + counting")
            design = pd.read_csv(W["design"], sep="\t", dtype={"OligoID": str})
            m = ExactMapper(design)
            s0 = design["MappingSequence"].iloc[3]
            check(m.hits(s0)[0][0] == design["OligoID"].iloc[3] and m.hits(revcomp(s0))[0][2] == "-", "exact mapper: forward and reverse-complement hits")
            check(m.hits(s0[:-1] + ("A" if s0[-1] != "A" else "C")) == [], "exact mapper: one mismatch does not align (-v0)")
            sam = d / "x.sam"
            sam.write_text("@HD\tVN:1.0\nr1\t0\tg1\t1\nr2\t0\tg1\t1\nr3\t4\t*\t0\nr4\t0\tg2\t1\n")
            c = count_reads(sam)
            check(c.set_index("OligoID")["count"].to_dict() == {"g1": 2, "g2": 1}, "count_reads = uniq -c of column 3 (SAM header / unmapped skipped)")
            rc = main(["map-reads", "--fastq", str(ss["fastqR1"].iloc[0]), "--sample-id", "X", "--label", "st_ffp_map0"])
            newest("st_ffp_map0")
            check(rc == 0, "map-reads without bowtie or design prints the bowtie command and exits 0")

            print("\nsort params loaders")
            sp = load_sort_params(W["sortparams"] / "B1_S1.txt", list("ABCDEF"))
            check(abs(sp["mu_seed"] - 2.82) < 1e-9 and sp["sd_seed"] == 0.5 and sp["totalCount"] == 0, "BigFoot_noTotal: mu seed = log10 Min of bin 3, sd seed 0.5")
            spc = load_sort_params(W["sortparams"] / "B1_S1.txt", list("ABCDEF"), upstream_compat=True)
            check(np.allclose(sp["bins"]["count"], spc["bins"]["count"]), "alphabetical sort params: R factor-code Count == by-name Count")
            scr = pd.read_csv(W["sortparams"] / "B1_S1.txt", sep="\t").iloc[[5, 0, 1, 2, 3, 4]]
            scr.to_csv(d / "scr.txt", sep="\t", index=False)
            a = load_sort_params(d / "scr.txt", list("ABCDEF"))["bins"]["count"]
            b = load_sort_params(d / "scr.txt", list("ABCDEF"), upstream_compat=True)["bins"]["count"]
            check(a["F"] == scr["Count"].iloc[0] and b["F"] != a["F"], "--upstream-compat reproduces R 3.6 factor-code Count indexing on a reordered file")
            ast = pd.DataFrame({"Barcode": ["A", "B", "C", "D", "Total"], "Mean": [100, 300, 1000, 3000, 900],
                                "Std.Dev.": [50, 100, 300, 900, 1200], "Bounds": ["(50, 200)", "(200, 500)", "(500, 2000)", "(2000, 9000)", "(1, 1e5)"],
                                "Count": [10, 20, 30, 40, 1000]})
            ast.to_csv(d / "ast.txt", sep="\t", index=False)
            sa = load_sort_params(d / "ast.txt", list("ABCD"), "astrios")
            check(sa["totalCount"] == 1000 and abs(sa["mu_seed"] - math.log10(2000)) < 1e-12
                  and abs(sa["sd_seed"] - math.log10(1 + 0.25)) < 1e-12, "Astrios: Total count, mu seed = 4th bin lower bound, sd seed from row 1")
            check(load_sort_params(d / "ast.txt", list("ABCD"), "astrios_nototal")["totalCount"] == 0, "Astrios_nototal: total 0")
            bf = pd.DataFrame({"Barcode": ["A", "B", "C", "Total"], "Mean": [100, 300, 1000, 900], "Min": [50, 200, 500, 1],
                               "Max": [200, 500, 2000, 1e5], "Count": [10, 20, 30, 1000], "StdDev": [1, 1, 1, 900]})
            bf.to_csv(d / "bf.csv", index=False)
            sb = load_sort_params(d / "bf.csv", list("ABC"), "bigfoot")
            check(sb["totalCount"] == 1000 and abs(sb["sd_seed"] - math.log10(2)) < 1e-12, "BigFoot (csv): Total count and sd seed log10(1+(SD/Mean)^2)")

            print("\nMLE")
            lb = np.array([GATES[k][0] for k in "ABCDEF"]); ub = np.array([GATES[k][1] for k in "ABCDEF"])
            truth = 2.95
            pr = stats.norm.cdf(ub, truth, 0.25) - stats.norm.cdf(lb, truth, 0.25)
            mu_hat, sd_hat, meth = normal_mle_guide(2.9, 0.5, pr * 5000, lb, ub)
            check(abs(mu_hat - truth) < 0.01 and abs(sd_hat - 0.25) < 0.02 and meth == "EM", f"MLE with EM seventh bin recovers mu={truth}, sd=0.25 (got {mu_hat:.3f}, {sd_hat:.3f})")
            try:
                normal_mle_guide(2.9, 0.5, np.zeros(6), lb, ub); check(False, "unobserved guide -> singular Hessian error")
            except MLEError:
                check(True, "unobserved guide -> singular Hessian error (falls back to weighted_avg)")

            print("\nfull workflow (run)")
            rc = main(["run", "--sample-sheet", str(W["sample_sheet"]), "--design", str(W["design"]), "--sortparams-dir", str(W["sortparams"]),
                       "--fastq-dir", str(W["fastq"]), "--experiment-keycols", "CellLine", "--replicate-keycols", "FlowFISHRep",
                       "--qpcr", str(W["qpcr"]), "--genelist", str(W["genes"]), "--enhancers", str(W["enhancers"]),
                       "--cell-line", "K562", "--python-mapper", "--power", "--label", "st_ffp_run"] + (["--no-plots"] if args.no_plots else []))
            R = newest("st_ffp_run")
            res = R / "results"
            check(rc == 0 and (R / "report.md").exists() and (R / "summary.json").exists(), "run writes report.md and summary.json")
            al = pd.read_csv(res / "summary" / "alignment_stats.tsv", sep="\t")
            check(len(al) == 26 and (al["nReadsAligned"] == 8000).all() and (al["nReadsUnaligned"] == 160).all(), "alignment_stats: 8000 aligned, 2% junk unaligned per sample")
            pc = pd.read_csv(res / "byPCRRep" / "K562-Rep1-PCR1.bin_counts.txt", sep="\t")
            er = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.bin_counts.txt", sep="\t")
            check(list(pc.columns) == ["OligoID"] + list("ABCDEF") + ["All"] and list(er.columns) == ["OligoID"] + list("ABCDEF") + ["All"],
                  "count tables: bins then All, OligoID index")
            pc2 = pd.read_csv(res / "byPCRRep" / "K562-Rep1-PCR2.bin_counts.txt", sep="\t")
            mm = pc.set_index("OligoID")[list("ABCDEF")].add(pc2.set_index("OligoID")[list("ABCDEF")], fill_value=0)
            check(np.allclose(mm.loc[er["OligoID"]].to_numpy(), er[list("ABCDEF")].to_numpy()), "experimental-rep counts = sum of PCR reps per bin")
            fr = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.bin_freq.txt", sep="\t")
            check(np.allclose(fr[list("ABCDEF")].sum(), 1, atol=1e-4), "bin_freq columns sum to 1")
            flat = pd.read_csv(res / "summary" / "GuideCounts.flat.tsv.gz", sep="\t")
            check(list(flat.columns) == ["count", "OligoID", "SampleID"] and flat["count"].sum() == 26 * 8000, "flat count table")
            cor = pd.read_csv(res / "summary" / "ReplicateCorrelation.pdf.tsv", sep="\t")
            check(len(cor) == 12 and cor["R"].min() > 0.8, f"PCR replicate correlations per bin (min R {cor['R'].min():.2f})")
            raw = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.raw_effects.txt", sep="\t")
            need = ["sumreads", "numobsbins", "sumcells", "WeightedAvg", "input.fraction", "logMean", "logSD", "method"]
            check(all(c in raw.columns for c in need) and (raw["method"] == "EM").mean() > 0.95, "raw_effects columns; MLE converges for >95% of guides")
            spp = pd.read_csv(W["sortparams"] / "B1_S1.txt", sep="\t")
            check(np.allclose(raw[list("ABCDEF")].sum().to_numpy(), spp["Count"].to_numpy()), "read counts rescaled to sorted-cell counts per bin")
            gi = W["guides"].set_index("OligoID")["_eff"]
            raw["eff"] = raw["OligoID"].map(gi)
            rho_mle = stats.spearmanr(raw["logMean"], np.log10(raw["eff"]))[0]
            r_wa = stats.spearmanr(raw["WeightedAvg"], np.log10(raw["eff"]))[0]
            check(rho_mle > 0.7 and r_wa > 0.7, f"logMean and WeightedAvg track planted log effects (rho {rho_mle:.2f}, {r_wa:.2f})")
            real = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.real_space.txt", sep="\t")
            check(list(real.columns) == BED_COLS, "real_space.txt columns as upstream")
            check(not real["OligoID"].isin(["BAD_lowOT", "BAD_lowOT2", "BAD_short"]).any(), "guide filters drop low off-target and short guides")
            check("NC_0" in set(real["OligoID"]), "missing GuideSequence falls back to GuideSequenceMinusG")
            ncm = real.loc[real["target"] == NEG, "mleAvg"].mean()
            check(abs(ncm - 1) < 1e-9 and real["mleAvg"].between(0, 5.0 / 0.9).all(), "negative controls re-scaled to mean 1 after clamp")
            real["eff"] = real["OligoID"].map(gi)
            tss_med = real.loc[real["eff"] == 0.3, "mleAvg"].median()
            e1_med = real.loc[(real["eff"] == 0.5), "mleAvg"].median()
            check(0.2 < tss_med < 0.42 and 0.38 < e1_med < 0.65, f"real-space effects near planted (TSS {tss_med:.2f} ~0.3, E1 {e1_med:.2f} ~0.5)")
            win = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.windows.txt", sep="\t")
            check(len(win) > 0 and (win["span"] <= 750).all() and (win["n"] == 10).all(), "10-guide windows with span <= 750")
            check(win.sort_values("mean").iloc[0]["start"] < 100_250, "lowest window sits at the TSS")
            info = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.ScreenInfo.txt", sep="\t")
            ff = info["FlowFISH_at_TSS"].iloc[0]
            best = win[(win["chr"] == "chr1") & (win["start"] < 100_500) & (win["end"] > 99_500)]["mean"].min()
            check(abs(ff - (1 + best)) < 1e-12 and abs(info["Background"].iloc[0] - (ff - 0.15) / (1 - ff)) < 1e-12,
                  "ScreenInfo: FlowFISH_at_TSS = 1 + best window mean; Background = (FF-qPCR)/(1-FF)")
            sc = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.scaled.txt", sep="\t", dtype={"OligoID": str})
            slope = (1 - 0.15) / (1 - ff)
            rr = real.set_index("OligoID")["mleAvg"]
            exp_sc = (slope * (rr.loc[sc["OligoID"]].to_numpy() - ff) + 0.15).clip(0, 5)
            check(np.allclose(sc["mleAvg"].to_numpy(), exp_sc), "scaled = slope (x - FF) + qPCR, clipped [0,5]")
            colb = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.collapse.bed", sep="\t", header=None)
            check(len(colb) == 7 and not (colb[1] == 300_000).any(), "collapse: 6 tiled peaks + chr2 element; guide-less element dropped")
            full = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.FullEnhancerScore.txt", sep="\t")
            check(list(full.columns) == SCORE_COLS and len(full) == 6, "FullEnhancerScore columns; chr2 element (2 guides < 3) dropped")
            sig = dict(zip(full["start"], full["Significant"]))
            reg = dict(zip(full["start"], full["Regulated"]))
            check(sig[99_750] and sig[130_000] and sig[160_000] and reg[130_000], "planted TSS, E1 (50%) and E2 (28%) peaks significant and regulated")
            check(not any(sig[s] for s in (180_000, 200_000, 220_000)), "null peaks not significant")
            ph = full[full["start"] == 130_000].iloc[0]
            scores = [float(v) for v in colb[colb[1] == 130_000][3].iloc[0].split(",")]
            ncv = sc.loc[sc["target"] == NEG, "mleAvg"].to_numpy()
            check(abs(ph["p.ttest"] - stats.ttest_ind(scores, ncv, equal_var=True)[1]) < 1e-12 and ph["n"] == len(scores) == 16,
                  "p.ttest = Student t of the collapsed guides vs controls")
            pk = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.peaks.tfdr0.05.bed", sep="\t", header=None)
            check(len(pk) == int(full["Significant"].sum()), "peaks.tfdr0.05.bed lists the significant peaks")
            pcs = pd.read_csv(res / "summary" / "PeakCallingSummary.tsv", sep="\t")
            check(len(pcs) == 2 and (pcs["all.elements"] == 6).all() and (pcs["GuidePerDHS.MINGUIDES"] == 3).all(), "PeakCallingSummary.tsv per replicate")
            ke = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.KnownEnhancers.FlowFISH.txt", sep="\t")
            e1 = ke[ke["start"] == 130_000].iloc[0]
            check(e1["GeneTSS"] == 100_000 and e1["distance"] == 30_250 and e1["Fraction.change.in.gene.expr"] < -0.3, "KnownEnhancers: TSS, distance, fraction change")
            sd_txt = (res / "byExperimentRep" / "K562-Rep1.ScreenData.txt").read_text().splitlines()
            check(sd_txt[0] == "Screen\tScreenData\tRNAReadoutMethod\tReference" and sd_txt[1].startswith("GENE1\t"), "ScreenData.txt format")
            pw = pd.read_csv(res / "byExperimentRep" / "K562-Rep1.BH-Power-fdr0.05.txt", sep="\t")
            check(pw["Power"].is_monotonic_increasing or pw["Power"].iloc[-1] >= pw["Power"].iloc[0], "power grows with effect size")
            gcq = pd.read_csv(res / "summary" / "GuideCounts.replicateCorrelations.tsv", sep="\t")
            check(len(gcq) == 12 and gcq["cor"].min() > 0.8, "PlotGuideCounts: PCR-replicate %Reads correlations")

            print("\nstep subcommands on the run outputs")
            bc = res / "byExperimentRep" / "K562-Rep2.bin_counts.txt"
            rc = main(["estimate-effects", "--counts", str(bc), "--sort-params", str(W["sortparams"] / "B1_S2.txt"), "--label", "st_ffp_est"])
            E = newest("st_ffp_est")
            a = pd.read_csv(E / "K562-Rep2.raw_effects.txt", sep="\t")
            b2 = pd.read_csv(res / "byExperimentRep" / "K562-Rep2.raw_effects.txt", sep="\t")
            check(rc == 0 and np.allclose(a["logMean"], b2["logMean"]), "estimate-effects reproduces the run")
            rc = main(["real-space", "--mle", str(E / "K562-Rep2.raw_effects.txt"), "--design", str(W["design"]), "--minsum", "0.05",
                       "--minbins", "4", "--label", "st_ffp_real"])
            Rs = newest("st_ffp_real")
            rc2 = main(["windows", "--input", str(Rs / "K562-Rep2.real_space.txt"), "--label", "st_ffp_win"])
            Wn = newest("st_ffp_win")
            rc3 = main(["tss-kd", "--gene", "GENE1", "--windows", str(Wn / "K562-Rep2.windows.txt"), "--qpcr", str(W["qpcr"]),
                        "--genelist", str(W["genes"]), "--label", "st_ffp_tss"])
            T = newest("st_ffp_tss")
            rc4 = main(["normalize-qpcr", "--info", str(T / "K562-Rep2.ScreenInfo.txt"), "--input", str(Rs / "K562-Rep2.real_space.txt"),
                        "--label", "st_ffp_norm"])
            N = newest("st_ffp_norm")
            s2 = pd.read_csv(N / "K562-Rep2.scaled.txt", sep="\t")
            s2r = pd.read_csv(res / "byExperimentRep" / "K562-Rep2.scaled.txt", sep="\t")
            check(rc == rc2 == rc3 == rc4 == 0 and np.allclose(s2["mleAvg"], s2r["mleAvg"]), "real-space -> windows -> tss-kd -> normalize-qpcr reproduce the run")
            rc = main(["collapse", "--scaled", str(N / "K562-Rep2.scaled.txt"), "--enhancers", str(W["enhancers"]), "--label", "st_ffp_col"])
            Cc = newest("st_ffp_col")
            rc2 = main(["score-enhancers", "--collapsed", str(Cc / "K562-Rep2.collapse.bed"), "--scaled", str(N / "K562-Rep2.scaled.txt"),
                        "--expt-name", "K562-Rep2", "--label", "st_ffp_score", "--no-plots"])
            S = newest("st_ffp_score")
            f2 = pd.read_csv(S / "K562-Rep2.FullEnhancerScore.txt", sep="\t")
            f2r = pd.read_csv(res / "byExperimentRep" / "K562-Rep2.FullEnhancerScore.txt", sep="\t")
            check(rc == rc2 == 0 and np.allclose(f2["fdr.ttest"], f2r["fdr.ttest"]) and (f2["Significant"] == f2r["Significant"]).all(),
                  "collapse -> score-enhancers reproduce the run")
            rc = main(["power", "--collapsed", str(Cc / "K562-Rep2.collapse.bed"), "--scaled", str(N / "K562-Rep2.scaled.txt"),
                       "--info", str(T / "K562-Rep2.ScreenInfo.txt"), "--reps", "50", "--label", "st_ffp_pow"])
            P = newest("st_ffp_pow")
            pw = pd.read_csv(P / "K562-Rep2.BH-Power-fdr0.05.txt", sep="\t")
            check(rc == 0 and len(pw) == 6 and pw["Power"].iloc[-1] >= pw["Power"].iloc[0], "power subcommand")
            rc = main(["format-screen", "--score", str(S / "K562-Rep2.FullEnhancerScore.txt"), "--gene", "GENE1", "--genelist", str(W["genes"]),
                       "--cell-line", "K562", "--label", "st_ffp_fmt"])
            newest("st_ffp_fmt")
            rc2 = main(["guide-count-plots", "--guide-counts", str(res / "summary" / "GuideCounts.flat.tsv.gz"),
                        "--sample-sheet", str(R / "SampleList.snakemake.tsv"), "--label", "st_ffp_gc"] + (["--no-plots"] if args.no_plots else []))
            newest("st_ffp_gc")
            rc3 = main(["count-tables", "--sample-sheet", str(W["sample_sheet"]), "--counts-dir", str(res / "counts"),
                        "--experiment-keycols", "CellLine", "--replicate-keycols", "FlowFISHRep", "--label", "st_ffp_ct", "--no-plots"])
            Ct = newest("st_ffp_ct")
            e3 = pd.read_csv(Ct / "results" / "byExperimentRep" / "K562-Rep1.bin_counts.txt", sep="\t")
            rc4 = main(["samplesheet", "--sample-sheet", str(W["sample_sheet"]), "--experiment-keycols", "CellLine",
                        "--replicate-keycols", "FlowFISHRep", "--fastq-dir", str(W["fastq"]), "--label", "st_ffp_ss"])
            newest("st_ffp_ss")
            rc5 = main(["count-reads", "--alignments", str(res / "samFiles" / "R1_A_P1.sam.gz"), "--sample-id", "R1_A_P1", "--label", "st_ffp_cr"])
            Cr = newest("st_ffp_cr")
            c1 = (Cr / "R1_A_P1.count.txt").read_text() == (res / "counts" / "R1_A_P1.count.txt").read_text()
            rc6 = main(["build-index", "--design", str(W["design"]), "--label", "st_ffp_idx"])
            newest("st_ffp_idx")
            rc7 = main(["replicate-correlations", "--sample-sheet", str(R / "SampleList.snakemake.tsv"), "--base-dir", str(R),
                        "--label", "st_ffp_rc", "--no-plots"])
            newest("st_ffp_rc")
            fq = str(ss["fastqR1"].iloc[0])
            rc8 = main(["map-reads", "--fastq", fq, "--sample-id", "R1_A_P1", "--design", str(W["design"]), "--label", "st_ffp_map"])
            Mp = newest("st_ffp_map")
            c2 = (Mp / "results" / "counts" / "R1_A_P1.count.txt").read_text() == (res / "counts" / "R1_A_P1.count.txt").read_text()
            check(all(r == 0 for r in (rc, rc2, rc3, rc4, rc5, rc6, rc7, rc8)) and e3.equals(er) and c1 and c2,
                  "format-screen, guide-count-plots, count-tables, samplesheet, count-reads, build-index, replicate-correlations, map-reads")
            rc = main(["estimate-effects", "--counts", str(bc), "--sort-params", str(W["sortparams"] / "B1_S2.txt"), "--upstream-compat",
                       "--label", "st_ffp_estc"])
            Ec = newest("st_ffp_estc")
            ac = pd.read_csv(Ec / "K562-Rep2.raw_effects.txt", sep="\t")
            check(rc == 0 and np.allclose(ac["logMean"], a["logMean"]), "--upstream-compat is identical when the sort params are alphabetical")
    finally:
        for p in made:
            shutil.rmtree(p, ignore_errors=True)
    ok = bool(checks) and all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.0f}s)")
    print("selftest: all checks pass" if ok else f"selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent flowfish-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def keys(p, sheet=True):
        if sheet:
            p.add_argument("--sample-sheet", required=True, help="tab-delimited sample sheet (SampleID Bin PCRRep GeneSymbol qPCRGene + keys)")
        p.add_argument("--experiment-keycols", default="", help="comma-separated experiment key columns (config experiment_keycols)")
        p.add_argument("--replicate-keycols", default="", help="comma-separated replicate key columns (config replicate_keycols)")

    def est_opts(p):
        p.add_argument("--sort-params-format", default="bigfoot_nototal", choices=SORT_FORMATS)
        p.add_argument("--upstream-compat", action="store_true", help="R 3.6 factor-code Count indexing and positional bin pairing")

    p = sub.add_parser("samplesheet", help="validate the sample sheet, find FASTQs, add experiment keys")
    keys(p); p.add_argument("--fastq-dir"); p.add_argument("--label", default="samplesheet")
    p.set_defaults(func=cmd_samplesheet)

    p = sub.add_parser("build-index", help="guide FASTA (+ bowtie-build when on PATH)")
    p.add_argument("--design", required=True); p.add_argument("--label", default="guide_index")
    p.set_defaults(func=cmd_build_index)

    p = sub.add_parser("map-reads", help="bowtie -v0 --all (or exact Python mapper) + count_reads + alignment stats")
    p.add_argument("--fastq", required=True); p.add_argument("--sample-id", required=True)
    p.add_argument("--index", help="bowtie index prefix"); p.add_argument("--design", help="guide design (Python mapper)")
    p.add_argument("--python-mapper", action="store_true", help="use the Python mapper even if bowtie is on PATH")
    p.add_argument("--label", default="map_reads")
    p.set_defaults(func=cmd_map_reads)

    p = sub.add_parser("count-reads", help="alignment file -> {SampleID}.count.txt")
    p.add_argument("--alignments", required=True); p.add_argument("--sample-id", required=True); p.add_argument("--label", default="count_reads")
    p.set_defaults(func=cmd_count_reads)

    p = sub.add_parser("count-tables", help="bin count / frequency tables per PCR rep and experimental rep + flat table")
    keys(p); p.add_argument("--counts-dir", required=True, help="directory of {SampleID}.count.txt")
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--label", default="count_tables")
    p.set_defaults(func=cmd_count_tables)

    p = sub.add_parser("replicate-correlations", help="plotReplicateCorrelations.R")
    p.add_argument("--sample-sheet", required=True, help="SampleList.snakemake.tsv")
    p.add_argument("--base-dir", required=True, help="directory the count-file paths are relative to")
    p.add_argument("--groupcol", default="ExperimentIDReplicates"); p.add_argument("--filecol", default="ExperimentIDPCRRep_BinCounts")
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--label", default="replicate_correlations")
    p.set_defaults(func=cmd_replicate_correlations)

    p = sub.add_parser("estimate-effects", help="estimate_effect_sizes.R: weighted average + binned MLE")
    p.add_argument("--counts", required=True); p.add_argument("--sort-params", required=True)
    est_opts(p)
    p.add_argument("--mu-seed", type=float); p.add_argument("--sd-seed", type=float)
    p.add_argument("--name"); p.add_argument("--label", default="estimate_effects")
    p.set_defaults(func=cmd_estimate_effects)

    p = sub.add_parser("real-space", help="convert_to_real_space.py")
    p.add_argument("--mle", required=True); p.add_argument("--design", required=True)
    p.add_argument("--clamp", type=float, default=5.0); p.add_argument("--background", type=float, default=0.0)
    p.add_argument("--maxg", type=int, default=10); p.add_argument("--maxlen", type=int, default=21); p.add_argument("--minlen", type=int, default=18)
    p.add_argument("--minofftarget", type=float, default=50); p.add_argument("--minsum", type=float, default=0.0)
    p.add_argument("--minbins", type=int, default=0)
    p.add_argument("--name"); p.add_argument("--label", default="real_space")
    p.set_defaults(func=cmd_real_space)

    p = sub.add_parser("windows", help="CalculateTilingStatistic.R sliding guide windows")
    p.add_argument("--input", required=True); p.add_argument("--window", type=int, default=10); p.add_argument("--max-span", type=float, default=750)
    p.add_argument("--score-column", default="mleAvg"); p.add_argument("--name"); p.add_argument("--label", default="windows")
    p.set_defaults(func=cmd_windows)

    p = sub.add_parser("tss-kd", help="FlowFISHtssKD.py -> ScreenInfo.txt")
    p.add_argument("--gene", required=True); p.add_argument("--windows", required=True); p.add_argument("--qpcr")
    p.add_argument("--genelist", required=True); p.add_argument("--slop", type=int, default=500)
    p.add_argument("--name"); p.add_argument("--label", default="tss_kd")
    p.set_defaults(func=cmd_tss_kd)

    p = sub.add_parser("normalize-qpcr", help="normalize_flowfish_to_qpcr.py")
    p.add_argument("--info", required=True); p.add_argument("--input", required=True, help="real_space.txt")
    p.add_argument("--clamp", type=float, default=5.0); p.add_argument("--name"); p.add_argument("--label", default="normalize_qpcr")
    p.set_defaults(func=cmd_normalize_qpcr)

    p = sub.add_parser("collapse", help="enhancer list x guides -> collapse.bed")
    p.add_argument("--scaled", required=True); p.add_argument("--enhancers", required=True)
    p.add_argument("--name"); p.add_argument("--label", default="collapse")
    p.set_defaults(func=cmd_collapse)

    def score_opts(p):
        p.add_argument("--collapsed", required=True); p.add_argument("--scaled", required=True)
        p.add_argument("--expt-name", default=""); p.add_argument("--min-guides", type=int, default=3)
        p.add_argument("--fdr", type=float, default=0.05)

    p = sub.add_parser("score-enhancers", help="ScoreEnhancers.py")
    score_opts(p); p.add_argument("--min-effect", type=float, default=0.0)
    p.add_argument("--utest-mode", default="legacy", choices=["legacy", "two-sided"])
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--label", default="score_enhancers")
    p.set_defaults(func=cmd_score_enhancers)

    p = sub.add_parser("power", help="power analysis of the peak t-test")
    score_opts(p); p.add_argument("--info", help="ScreenInfo.txt (Background)")
    p.add_argument("--reps", type=int, default=400); p.add_argument("--effect-sizes", type=float, nargs="+", default=[0.05, 0.10, 0.15, 0.2, 0.25, 0.3])
    p.add_argument("--seed", type=int, default=1); p.add_argument("--label", default="power")
    p.set_defaults(func=cmd_power)

    p = sub.add_parser("format-screen", help="ScreenData.txt + KnownEnhancers.FlowFISH.txt")
    p.add_argument("--score", required=True, help="FullEnhancerScore.txt"); p.add_argument("--gene", required=True)
    p.add_argument("--genelist"); p.add_argument("--cell-line", default=""); p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--name"); p.add_argument("--label", default="format_screen")
    p.set_defaults(func=cmd_format_screen)

    p = sub.add_parser("guide-count-plots", help="PlotGuideCounts.R (guide-count version)")
    p.add_argument("--guide-counts", required=True, help="GuideCounts.flat.tsv.gz"); p.add_argument("--sample-sheet", required=True)
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--label", default="guide_counts")
    p.set_defaults(func=cmd_guide_count_plots)

    p = sub.add_parser("run", help="the whole workflow from a sample sheet")
    keys(p)
    p.add_argument("--design", required=True); p.add_argument("--sortparams-dir", required=True, help="{Batch}_{Sample}.txt files")
    p.add_argument("--fastq-dir"); p.add_argument("--counts-dir", help="pre-computed {SampleID}.count.txt (skips mapping)")
    p.add_argument("--bowtie-index"); p.add_argument("--python-mapper", action="store_true")
    p.add_argument("--qpcr"); p.add_argument("--genelist"); p.add_argument("--enhancers"); p.add_argument("--cell-line", default="")
    est_opts(p)
    p.add_argument("--window", type=int, default=10); p.add_argument("--max-span", type=float, default=750)
    p.add_argument("--minsum", type=float, default=0.05); p.add_argument("--minbins", type=int, default=4)
    p.add_argument("--clamp", type=float, default=5.0); p.add_argument("--min-guides", type=int, default=3)
    p.add_argument("--fdr", type=float, default=0.05); p.add_argument("--min-effect", type=float, default=0.0)
    p.add_argument("--utest-mode", default="legacy", choices=["legacy", "two-sided"])
    p.add_argument("--power", action="store_true", help="also run the power analysis per replicate")
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--label", default="flowfish")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic screen with planted enhancers through every subcommand")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    global _LOG_PATH
    if _LOG_PATH is None:
        _LOG_PATH = setup_logging()
        print(f"Log: {_LOG_PATH}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
