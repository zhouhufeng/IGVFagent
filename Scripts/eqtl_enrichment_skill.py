#!/usr/bin/env python3
"""eQTL enrichment benchmark for enhancer-gene predictors (port of EngreitzLab/eQTLEnrichment).

Port of https://github.com/EngreitzLab/eQTLEnrichment (MIT; Snakemake + R +
bedtools), the Engreitz-lab "eQTL benchmarking pipeline" used in the
ENCODE-rE2G and scE2G papers to ask two questions of an enhancer-gene model:
(i) are its predicted enhancers enriched for fine-mapped eQTL variants, and
(ii) does it link those enhancers to the right eGene?  Every rule of the
workflow was read and re-derived in Python (pandas + numpy + scipy); nothing
was copied.  Relationship: port.

Definitions, exactly as upstream computes them
  distal noncoding   eQTL variants (PIP >= thresholdPIP, default 0.5) and the
                     background 1000G SNPs are both restricted to partition
                     categories ABC / AllPeaks / Other / OtherIntron, which
                     removes coding sequence, UTRs, splice sites and promoters
                     (TSS +/- 250 bp) of protein-coding genes.
  gene universe      eGenes and prediction TargetGenes are restricted to the
                     TSS reference (CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp).
  enrichment         (fraction of eQTL variant LOCATIONS in a tissue overlapping
                     predicted enhancers with score >= t) / (fraction of
                     background SNPs overlapping them); log risk-ratio CI,
                     hypergeometric p (phyper upper tail), Bonferroni over the
                     table.
  recall (total)     fraction of eVariant-eGene PAIRS whose variant overlaps a
                     predicted enhancer with score >= t.
  recall (linking)   ... overlapping an enhancer linked to the correct eGene.
  correctGene.ifOverlap  linking / total.
  distance bins      |mean TSS centre of the eGene - variant start|, bins
                     0-d1, d1-d2, ..., plus "all" (0-30 Mb); variants beyond
                     the last bin only count in "all".
  threshold span     quantiles of the matched, correctly-linked scores
                     (nThresholdSteps, grown until that many distinct values
                     exist), the method's own threshold, and an even grid from
                     min to max; binary predictors use {0, 1}.
  inverse predictor  scores and the threshold are negated before everything.
  aggregate          "all_matches" rows sum numerators and denominators over
                     the matched (eQTL tissue, prediction biosample) pairs.
  enrichment at recall  per tissue and target recall, the threshold row whose
                     recall (linking) is closest to the target if within 10%
                     of it; pairwise two-sided z-tests on the log enrichment
                     difference, Bonferroni.
  enhancer set size  bp covered by merged elements with score >= threshold
                     that do not overlap the TSS reference.

Two upstream quirks are corrected by default and reproduced with
--upstream-compat so numbers can be compared against a pipeline run:
  * counts_to_enrichment_*.R:  SE(log RR) is written as
    sqrt(((n1-x1)/x1)/n1) + ((n2-x2)/x2)/n2 -- the square root closes after
    the first term.  The port uses sqrt(a + b), the textbook formula the code
    cites.
  * count_matrix_distance.R / counts_to_enrichment_distance.R:  in the
    by-distance tables the background-SNP overlap count ignores the score
    threshold, and for the "all" bin the eQTL numerator does too, so the
    enrichment heatmaps show "any candidate element" rather than "predicted
    enhancer at the threshold".  The port thresholds both sides.

Subcommands
  setup         fetch the seven resource files (partition, TSS reference,
                gene bounds, chromosome sizes) from the pinned upstream commit
                into Data/eQTLEnrichment/resources; --synapse also pulls the
                background SNPs (syn52264319) and the GTEx SuSiE fine-mapping
                release (syn52264297) with `igvfagent synapse download`.
  prepare-gtex  raw GTEx 19-column fine-mapping table -> the 7-column eQTL
                input: method (SUSIE), credible-set members only, Ensembl ->
                HGNC through CollapsedGeneBounds.hg38.tsv, and eGene median
                TPM > 1 in its tissue from the GTEx .gct.
  variants      PIP filter, distal-noncoding filter, TSS distance bins, gene
                universe, per-tissue counts, background SNP filter; cached so
                several `run`s can reuse it (--variants-dir).
  run           the whole benchmark for a methods table + predictions table
                (upstream's config formats): per-method enrichment / recall
                tables across thresholds and by distance, enrichment-recall
                curves per matched pair and aggregated, enrichment at target
                recalls with pairwise tests, enhancer set sizes, heatmap
                matrices, figures, report.md and summary.json.
  selftest      synthetic genome with planted enhancers, eQTLs and background
                SNPs; runs prepare-gtex, variants and run on four predictors
                (good, random, inverse distance, binary) and asserts every
                table and definition.

Output: Docs/eQTLEnrichment/<timestamp>_<label>/.  pandas, numpy and scipy
required; matplotlib optional (figures skipped without it).

Usage:
    igvfagent eqtl-enrich setup [--synapse]
    igvfagent eqtl-enrich prepare-gtex --raw GTEx_30tissues_release1.tsv.gz --expression GTEx_median_tpm.gct.gz --out eqtl.tsv.gz
    igvfagent eqtl-enrich variants --eqtl eqtl.tsv.gz --bg-variants all.bg.SNPs.hg38.bed --label gtex
    igvfagent eqtl-enrich run --methods-table methods.tsv --predictions-table predictions.tsv --eqtl eqtl.tsv.gz --bg-variants bg.bed --label k562_blood
    igvfagent eqtl-enrich selftest --no-plots
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import math
import os
import subprocess
import sys
import time
import urllib.request
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "eQTLEnrichment"
DATA_ROOT = ROOT / "Data" / "eQTLEnrichment"
RES_DIR = DATA_ROOT / "resources"

UPSTREAM_REPO = "EngreitzLab/eQTLEnrichment"
UPSTREAM_COMMIT = "b04a2a1aaf4c5a2a811347289981e217579d239b"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/resources/genome_annotation/"
RESOURCE_FILES = {
    "PartitionCombined.bed": "partition of hg38 into ABC / AllPeaks / TSS / CDS / UTR / splice / intron / other",
    "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6": "TSS reference (500 bp windows) and gene universe",
    "CollapsedGeneBounds.hg38.TSS500bp.bed": "TSS reference, older gene set",
    "CollapsedGeneBounds.hg38.bed": "gene bounds (bed6)",
    "CollapsedGeneBounds.hg38.tsv": "gene bounds with Ensembl IDs and gene types",
    "GRCh38_EBV.chrom.sizes.tsv": "chromosome sizes incl. EBV",
    "GRCh38_main.chrom.sizes.tsv": "chromosome sizes, main chromosomes",
}
SYNAPSE_RESOURCES = {
    "syn52264319": "all.bg.SNPs.hg38.baseline.v1.1.bed.sorted (1000G background SNPs, LDSC baseline v1.1)",
    "syn52264297": "GTEx_30tissues_release1.tsv.gz (SuSiE / DAP-G fine-mapped GTEx v8 eQTLs, hg38)",
}
DISTAL_CATEGORIES = {"ABC", "AllPeaks", "Other", "OtherIntron"}
EQTL_COLS = ["chr", "start", "end", "variant_id", "gene", "tissue", "pip"]
EQTL_INPUT_ALIASES = {"varID_hg38": "variant_id", "variant_id": "variant_id", "gene_hgnc": "gene", "gene": "gene",
                      "chr": "chr", "start": "start", "end": "end", "tissue": "tissue", "pip": "pip", "PIP": "pip"}
METHOD_COLS = ["method", "boolean", "inverse_predictor", "pred_name_long", "threshold", "score_col", "color"]
GTEX_RAW_COLS = ["chr", "start", "end", "variant_id", "varID_hg38", "allele1", "allele2", "cohort", "method", "tissue",
                 "gene_ensembl", "maf", "beta_marginal", "se_marginal", "z", "pip", "cs_id", "beta_posterior", "sd_posterior"]
ALL_BIN_KB = 30000
BIN = 1 << 16
EMPTY_FILE_BYTES = 100

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948",
          "#792374", "#006479", "#0096a0", "#d3a9ce", "#96a0b3", "#5e5948", "#916953", "#888364"]
DIST_COLORS = ["#003648", "#006479", "#0096a0", "#49bcbc", "#96ced3", "#cae5ee"]

log = logging.getLogger("eqtl_enrich")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"eqtl_enrichment_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "t", "1", "yes", "y"}


def write_tsv(df, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)
    print(f"TSV: {path}")
    return path


def _header_line(path: Path) -> "tuple[list[str], int]":
    """Header columns and the number of lines to skip to reach the data.

    Leading '#' comment lines without a tab (IGVF portal metadata such as '# Source: scE2G') are skipped; a
    '#'-prefixed tab-separated line is the header with the '#' removed (upstream: awk 'NR==1{sub(/^#*/, "")}').
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    skip = 0
    with opener(path, "rt") as fh:
        for line in fh:
            skip += 1
            if not line.strip():
                continue
            if line.startswith("#") and "\t" not in line:
                continue
            return line.lstrip("#").rstrip("\n").split("\t"), skip
    return [], skip


def read_table_header(path: Path) -> "list[str]":
    return _header_line(path)[0]


def read_tsv(path: Path, usecols: "Optional[list[str]]" = None, header: bool = True, names: "Optional[list[str]]" = None,
             dtype=None, **kw):
    """TSV with an optional '#'-prefixed header; comment lines above the header are skipped."""
    pd = _pd()
    if header:
        cols, skip = _header_line(path)
        return pd.read_csv(path, sep="\t", names=cols, skiprows=skip, usecols=usecols, dtype=dtype, low_memory=False, **kw)
    return pd.read_csv(path, sep="\t", header=None, names=names, usecols=usecols, dtype=dtype, comment="#",
                       low_memory=False, **kw)


# ---------------------------------------------------------------------------
# Interval overlap (bedtools intersect -wa -wb, half-open BED semantics)
# ---------------------------------------------------------------------------

def overlap_join(a, b, a_cols=("chr", "start", "end"), b_cols=("chr", "start", "end"), suffixes=("", "_b")):
    """Cross rows of `a` and `b` that overlap by >= 1 bp: same chromosome, a.start < b.end and a.end > b.start.

    Binned join per chromosome (bins of 65,536 bp; intervals are replicated across the bins they
    span, then de-duplicated), so it is vectorised with pandas merges and needs no bedtools.
    Returns a DataFrame with every column of `a` followed by every column of `b` (suffixed on clash).
    """
    pd, np = _pd(), _np()
    if len(a) == 0 or len(b) == 0:
        cols = list(a.columns) + [c if c not in a.columns else c + suffixes[1] for c in b.columns]
        return pd.DataFrame(columns=cols)
    ac, as_, ae = a_cols
    bc, bs, be = b_cols
    a = a.reset_index(drop=True)
    b = b.reset_index(drop=True)
    pieces = []
    for chrom in sorted(set(a[ac].unique()) & set(b[bc].unique())):
        ai = np.flatnonzero((a[ac] == chrom).to_numpy())
        bi = np.flatnonzero((b[bc] == chrom).to_numpy())
        if ai.size == 0 or bi.size == 0:
            continue
        a_start = a[as_].to_numpy()[ai].astype(np.int64)
        a_end = np.maximum(a[ae].to_numpy()[ai].astype(np.int64), a_start + 1)
        b_start = b[bs].to_numpy()[bi].astype(np.int64)
        b_end = np.maximum(b[be].to_numpy()[bi].astype(np.int64), b_start + 1)

        def explode(idx, start, end):
            lo, hi = start // BIN, (end - 1) // BIN
            n = hi - lo + 1
            rep = np.repeat(np.arange(idx.size), n)
            offs = np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n)
            return pd.DataFrame({"bin": np.repeat(lo, n) + offs, "i": rep})

        ea, eb = explode(ai, a_start, a_end), explode(bi, b_start, b_end)
        m = ea.merge(eb, on="bin", suffixes=("_a", "_b"))
        if m.empty:
            continue
        ia, ib = m["i_a"].to_numpy(), m["i_b"].to_numpy()
        keep = (a_start[ia] < b_end[ib]) & (a_end[ia] > b_start[ib])
        ia, ib = ia[keep], ib[keep]
        if ia.size == 0:
            continue
        pairs = pd.DataFrame({"ia": ai[ia], "ib": bi[ib]}).drop_duplicates()
        pieces.append(pairs)
    if not pieces:
        cols = list(a.columns) + [c if c not in a.columns else c + suffixes[1] for c in b.columns]
        return pd.DataFrame(columns=cols)
    pairs = pd.concat(pieces, ignore_index=True)
    left = a.iloc[pairs["ia"].to_numpy()].reset_index(drop=True)
    right = b.iloc[pairs["ib"].to_numpy()].reset_index(drop=True)
    right.columns = [c if c not in left.columns else c + suffixes[1] for c in right.columns]
    return pd.concat([left, right], axis=1)


def merged_bp(df, chr_col="chr", start_col="start", end_col="end") -> int:
    """Total bases covered by the union of intervals (bedtools merge, then sum of lengths)."""
    total = 0
    for _, sub in df.groupby(chr_col, sort=False):
        s = sub[start_col].to_numpy(); e = sub[end_col].to_numpy()
        order = s.argsort()
        s, e = s[order], e[order]
        cur_s, cur_e = None, None
        for x, y in zip(s, e):
            if cur_e is None or x > cur_e:
                if cur_e is not None:
                    total += cur_e - cur_s
                cur_s, cur_e = x, y
            else:
                cur_e = max(cur_e, y)
        if cur_e is not None:
            total += cur_e - cur_s
    return int(total)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def enrichment_stats(x1, n1, x2, n2, alpha: float = 0.05, compat: bool = False):
    """Risk ratio (enrichment), CI on the log scale, hypergeometric p, as counts_to_enrichment_*.R.

    x1 / n1: eQTL variants overlapping / total; x2 / n2: background SNPs overlapping / total.
    p = phyper(x1, n1, n2, x1 + x2, lower.tail = FALSE) = P(X > x1) for X ~ Hypergeom(n1 + n2, n1, x1 + x2).
    """
    np = _np()
    from scipy.stats import hypergeom, norm  # type: ignore
    x1 = np.asarray(x1, dtype=float); n1 = np.asarray(n1, dtype=float)
    x2 = np.asarray(x2, dtype=float); n2 = np.asarray(n2, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        enr = (x1 / n1) / (x2 / n2)
        a = (n1 - x1) / (x1 * n1)
        b = (n2 - x2) / (x2 * n2)
        se = np.sqrt(a) + b if compat else np.sqrt(a + b)
        z = norm.isf(alpha / 2)
        lo = np.exp(np.log(enr) - z * se)
        hi = np.exp(np.log(enr) + z * se)
    ok = np.isfinite(x1) & np.isfinite(x2) & (n1 > 0) & (n2 > 0)
    p = np.full(x1.shape, np.nan)
    p[ok] = hypergeom.sf(x1[ok], (n1 + n2)[ok], n1[ok], (x1 + x2)[ok])
    return enr, lo, hi, se, p


def bonferroni(p):
    np = _np()
    p = np.asarray(p, dtype=float)
    m = int(np.isfinite(p).sum())
    return np.minimum(p * max(m, 1), 1.0)


def counts_at_thresholds(max_scores, thresholds) -> "list[int]":
    """Number of items whose best score is >= t, for every t (items with no overlap are absent)."""
    np = _np()
    s = np.sort(np.asarray(max_scores, dtype=float))
    return [int(s.size - np.searchsorted(s, t, side="left")) for t in thresholds]


# ---------------------------------------------------------------------------
# setup
# ---------------------------------------------------------------------------

def cmd_setup(args: argparse.Namespace) -> int:
    RES_DIR.mkdir(parents=True, exist_ok=True)
    got, skipped, failed = [], [], []
    for name in RESOURCE_FILES:
        dest = RES_DIR / name
        if dest.is_file() and dest.stat().st_size > 0 and not args.force:
            skipped.append(name)
            continue
        url = RAW + name
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            got.append(name)
            print(f"Wrote: {dest}")
        except Exception as e:  # noqa: BLE001
            failed.append(f"{name}: {e}")
    print(f"Resources in {RES_DIR}: fetched {len(got)}, already present {len(skipped)}, failed {len(failed)}")
    for f in failed:
        print("  FAILED", f, file=sys.stderr)
    if args.synapse:
        tok = os.environ.get("SYNAPSE_AUTH_TOKEN") or os.environ.get("SYNAPSE_PAT")
        if not tok:
            print("SYNAPSE_AUTH_TOKEN is not set; skipping the Synapse downloads", file=sys.stderr)
        else:
            for syn, desc in SYNAPSE_RESOURCES.items():
                print(f"Synapse {syn}: {desc}")
                cmd = [sys.executable, "-m", "igvfagent.cli", "synapse", "download", "--syn", syn, "--out-dir", str(DATA_ROOT)]
                rc = subprocess.call(cmd)
                if rc != 0:
                    failed.append(f"{syn} download exit {rc}")
    return 1 if failed else 0


def resource(name: str, override: "Optional[str]" = None) -> Path:
    if override:
        p = Path(override)
        if not p.is_file():
            raise SystemExit(f"resource not found: {p}")
        return p
    p = RES_DIR / name
    if not p.is_file():
        raise SystemExit(f"resource {name} missing; run `igvfagent eqtl-enrich setup` or pass the path explicitly")
    return p


# ---------------------------------------------------------------------------
# prepare-gtex: raw fine-mapping release -> 7-column eQTL input
# ---------------------------------------------------------------------------

def _norm_tissue(s: str) -> str:
    """'Brain - Cortex' (.gct header) and 'Brain_Cortex' (fine-mapping table) -> 'braincortex'."""
    return "".join(ch for ch in str(s) if ch.isalnum()).lower()


def read_gct_median_tpm(path: Path):
    """GTEx median-TPM .gct -> long table gene_symbol, tissue (GTEx underscore style), TPM."""
    pd = _pd()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline()
        skip = 2 if first.startswith("#1.2") or first.strip().startswith("#") else 0
    df = pd.read_csv(path, sep="\t", skiprows=skip)
    cols = list(df.columns)
    sym = cols[1]
    tissues = cols[2:]
    long = df.melt(id_vars=[sym], value_vars=tissues, var_name="tissue", value_name="TPM").rename(columns={sym: "gene"})
    long["tissue_key"] = long["tissue"].map(_norm_tissue)
    return long[["gene", "tissue_key", "TPM"]]


def cmd_prepare_gtex(args: argparse.Namespace) -> int:
    pd = _pd()
    raw = Path(args.raw)
    df = read_tsv(raw, header=False, names=GTEX_RAW_COLS)
    n0 = len(df)
    df = df[df["method"].astype(str) == args.method]
    if not args.keep_all_cs:
        df = df[df["cs_id"].astype(str) != "-1"]
    n1 = len(df)
    genes = read_tsv(resource("CollapsedGeneBounds.hg38.tsv", args.gene_table))
    genes = genes[["name", "Ensembl_ID"]].rename(columns={"name": "gene_hgnc", "Ensembl_ID": "ens"})
    df["ens"] = df["gene_ensembl"].astype(str).str[:15]
    df = df.merge(genes, on="ens", how="inner")
    n2 = len(df)
    out = df[["chr", "start", "end", "varID_hg38", "gene_hgnc", "tissue", "pip"]].drop_duplicates()
    n3 = len(out)
    n4 = None
    if args.expression:
        expr = read_gct_median_tpm(Path(args.expression))
        out = out.assign(tissue_key=out["tissue"].map(_norm_tissue))
        out = out.merge(expr, left_on=["gene_hgnc", "tissue_key"], right_on=["gene", "tissue_key"], how="left")
        missing = out["TPM"].isna().sum()
        out = out[out["TPM"] > args.tpm].drop(columns=["gene", "tissue_key", "TPM"])
        n4 = len(out)
        if missing:
            print(f"  note: {missing:,} rows had no expression value for their tissue (tissue names did not match the .gct) and were dropped")
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, sep="\t", index=False, compression="gzip" if str(dest).endswith(".gz") else None)
    print(f"GTEx variants: {n0:,} rows -> {n1:,} {args.method}{'' if args.keep_all_cs else ' in a credible set'} -> "
          f"{n2:,} in the gene universe -> {n3:,} unique" + (f" -> {n4:,} with eGene TPM > {args.tpm}" if n4 is not None else ""))
    print(f"Tissues: {out['tissue'].nunique()}; eVariant-eGene pairs at PIP >= 0.5: {int((out['pip'] >= 0.5).sum()):,}")
    print(f"Wrote: {dest}")
    return 0


# ---------------------------------------------------------------------------
# variants: PIP + distal noncoding + distance bins + gene universe + background
# ---------------------------------------------------------------------------

def load_partition_distal(path: Path):
    pd = _pd()
    part = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "category"],
                       dtype={"chr": str, "category": str})
    part = part[part["category"].isin(DISTAL_CATEGORIES)]
    return part[["chr", "start", "end"]].reset_index(drop=True)


def load_tss(path: Path):
    pd = _pd()
    tss = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3],
                      names=["chr", "start", "end", "gene"], dtype={"chr": str, "gene": str})
    tss["center"] = (tss["start"] + tss["end"]) / 2.0
    return tss


def distance_bins(distances_kb: "list[int]") -> "tuple[list[int], list[int]]":
    """upstream: distances_min = [0] + d[:-1] + [0]; distances_max = d + [30000]."""
    d = sorted(int(x) for x in distances_kb)
    return [0] + d[:-1] + [0], d + [ALL_BIN_KB]


def assign_distance_group(distance_bp, distances_kb: "list[int]"):
    """Group label = upper edge in kb of the bin (0-d1, d1-d2, ..., d_last-30000]; 0 if beyond 30 Mb."""
    np = _np()
    edges = [0] + sorted(int(x) for x in distances_kb) + [ALL_BIN_KB]
    grp = np.zeros(len(distance_bp), dtype=int)
    d = np.asarray(distance_bp, dtype=float)
    for lo, hi in zip(edges[:-1], edges[1:]):
        grp[(d > lo * 1000) & (d <= hi * 1000)] = hi
    return grp


def read_eqtl_input(path: Path):
    pd = _pd()
    cols = read_table_header(path)
    rename = {c: EQTL_INPUT_ALIASES[c] for c in cols if c in EQTL_INPUT_ALIASES}
    need = {"chr", "start", "end", "variant_id", "gene", "tissue", "pip"}
    if not need <= set(rename.values()):
        raise SystemExit(f"eQTL file needs columns chr,start,end,varID_hg38|variant_id,gene_hgnc|gene,tissue,pip; found {cols}")
    df = read_tsv(path, usecols=list(rename), dtype={"chr": str})
    df = df.rename(columns=rename)[EQTL_COLS]
    df["pip"] = pd.to_numeric(df["pip"], errors="coerce")
    return df


def prepare_variants(eqtl_path: Path, bg_path: Path, partition_path: Path, tss_path: Path, pip_threshold: float,
                     distances_kb: "list[int]", out: Path) -> dict:
    """The upstream rules filter_all_variants, bg_variant_count, add_distance_to_variants,
    filter_variants_to_gene_universe and get_variants_per_GTEx_tissue, on local files."""
    pd = _pd()
    out.mkdir(parents=True, exist_ok=True)
    part = load_partition_distal(partition_path)
    tss = load_tss(tss_path)
    universe = set(tss["gene"].astype(str))

    eq = read_eqtl_input(eqtl_path)
    n_raw = len(eq)
    eq = eq[eq["pip"] >= pip_threshold].drop_duplicates()
    n_pip = len(eq)
    hit = overlap_join(eq, part, suffixes=("", "_part"))
    eq = hit[EQTL_COLS].drop_duplicates().reset_index(drop=True)
    n_distal = len(eq)
    # distance to the eGene's mean TSS centre; variants whose gene has no TSS are dropped (inner join upstream)
    centers = tss.groupby("gene", sort=False)["center"].mean()
    eq = eq[eq["gene"].isin(centers.index)].copy()
    eq["distance_bin"] = assign_distance_group((centers.reindex(eq["gene"]).to_numpy() - eq["start"].to_numpy()).__abs__(),
                                               distances_kb)
    eq = eq[eq["gene"].isin(universe)].reset_index(drop=True)
    n_final = len(eq)
    var_path = write_tsv(eq, out / "eqtl_variants.filtered.tsv.gz")

    bg = pd.read_csv(bg_path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3],
                     names=["chr", "start", "end", "rsid"], dtype={"chr": str, "rsid": str})
    n_bg_raw = len(bg)
    bg = overlap_join(bg, part, suffixes=("", "_part"))[["chr", "start", "end", "rsid"]].drop_duplicates().reset_index(drop=True)
    n_bg = len(bg)
    bg_path_out = write_tsv(bg, out / "background_snps.distal_noncoding.bed.gz")
    (out / "background_snp_count.txt").write_text(f"{n_bg}\n")

    # per-tissue counts of distinct variant LOCATIONS, overall and per bin
    dmin, dmax = distance_bins(distances_kb)
    rows = []
    for tissue, sub in eq.groupby("tissue", sort=True):
        loc = sub[["chr", "start", "end"]].drop_duplicates()
        r = {"tissue": tissue, "n_all": len(loc)}
        for i, mx in enumerate(dmax, start=1):
            r[f"n_bin{i}"] = len(sub[sub["distance_bin"] == mx][["chr", "start", "end"]].drop_duplicates())
        rows.append(r)
    per_tissue = pd.DataFrame(rows)
    write_tsv(per_tissue, out / "variants_per_tissue.tsv")
    summary = {"eqtl_rows": n_raw, "pip_threshold": pip_threshold, "after_pip": n_pip, "after_distal_noncoding": n_distal,
               "after_gene_universe": n_final, "tissues": int(eq["tissue"].nunique()),
               "unique_evariants": int(len(eq[["chr", "start", "end"]].drop_duplicates())),
               "unique_evariant_egene_pairs": int(len(eq[["chr", "start", "end", "gene"]].drop_duplicates())),
               "background_snps": n_bg_raw, "background_snps_distal_noncoding": n_bg,
               "distances_kb": list(distances_kb), "distances_min": dmin, "distances_max": dmax,
               "variants_file": str(var_path), "background_file": str(bg_path_out)}
    (out / "variants_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"eQTL variants: {n_raw:,} rows -> PIP >= {pip_threshold}: {n_pip:,} -> distal noncoding: {n_distal:,} -> "
          f"gene universe: {n_final:,} ({summary['unique_evariants']:,} eVariants, {summary['unique_evariant_egene_pairs']:,} "
          f"eVariant-eGene pairs, {summary['tissues']} tissues)")
    print(f"Background SNPs: {n_bg_raw:,} -> distal noncoding: {n_bg:,}")
    return summary


def load_variants_dir(d: Path) -> "tuple[Any, Any, Any, dict]":
    pd = _pd()
    summ = json.loads((d / "variants_summary.json").read_text())
    eq = pd.read_csv(d / "eqtl_variants.filtered.tsv.gz", sep="\t", dtype={"chr": str, "variant_id": str, "gene": str, "tissue": str})
    bg = pd.read_csv(d / "background_snps.distal_noncoding.bed.gz", sep="\t", dtype={"chr": str, "rsid": str})
    per_tissue = pd.read_csv(d / "variants_per_tissue.tsv", sep="\t")
    return eq, bg, per_tissue, summ


def cmd_variants(args: argparse.Namespace) -> int:
    out = run_dir(args.label or "variants")
    prepare_variants(Path(args.eqtl), Path(args.bg_variants), resource("PartitionCombined.bed", args.partition),
                     resource("CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6", args.tss), args.threshold_pip,
                     args.distances, out)
    print(f"Variants dir: {out}")
    return 0


# ---------------------------------------------------------------------------
# predictions
# ---------------------------------------------------------------------------

def read_predictions(path: Path, score_col: str, universe: "set[str]", invert: bool, biosample: str):
    """process_predictions: chr,start,end,TargetGene,<score_col>; rows with blanks dropped; TargetGene in the
    gene universe; score negated for inverse predictors; biosample column added."""
    pd = _pd()
    cols = read_table_header(path)
    need = ["chr", "start", "end", "TargetGene", score_col]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(f"{path}: missing columns {missing}; available: {cols[:20]}...")
    df = read_tsv(path, usecols=need, dtype={"chr": str, "TargetGene": str})
    df = df.rename(columns={score_col: "score"})
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end", "TargetGene", "score"])
    df["start"] = df["start"].astype("int64"); df["end"] = df["end"].astype("int64")
    df = df[df["TargetGene"].isin(universe)]
    if invert:
        df["score"] = -df["score"]
    df["biosample"] = biosample
    return df[["chr", "start", "end", "biosample", "TargetGene", "score"]].reset_index(drop=True)


def element_max_scores(pred):
    """Unique elements with the best score over their genes (what a variant 'overlaps' at threshold t)."""
    return pred.groupby(["chr", "start", "end"], sort=False, as_index=False)["score"].max()


def threshold_span(var_int, matches: "list[tuple[str, str]]", n_steps: int, binary: bool, provided: float):
    """generate_quantile_threshold_span.R."""
    np = _np()
    if binary:
        return [0.0, 1.0]
    keys = {f"{t}.{b}" for t, b in matches}
    sub = var_int[(var_int["tissue"] + "." + var_int["biosample"]).isin(keys) & (var_int["gene"] == var_int["TargetGene"])]
    sub = sub[["chr", "start", "end", "gene", "score"]].drop_duplicates()
    scores = sub["score"].to_numpy(dtype=float)
    scores = scores[np.isfinite(scores)]
    if scores.size == 0:
        return sorted({float(provided)})
    n_distinct = np.unique(scores).size
    even = np.linspace(scores.min(), scores.max(), n_steps)
    steps = n_steps
    thr = np.unique(np.concatenate([np.quantile(scores, np.linspace(0, 1, steps)), [provided]]))
    while thr.size < n_steps and steps < n_distinct:
        steps += 1
        thr = np.unique(np.concatenate([np.quantile(scores, np.linspace(0, 1, steps)), [provided]]))
    return sorted(set(np.concatenate([thr, even]).tolist()))


# ---------------------------------------------------------------------------
# the benchmark
# ---------------------------------------------------------------------------

class Method:
    def __init__(self, row: dict, pred_files: "dict[str, Path]", matches: "list[tuple[str, str]]"):
        self.name = str(row["method"])
        self.boolean = _truthy(row.get("boolean", False))
        self.inverse = _truthy(row.get("inverse_predictor", False))
        self.long = str(row.get("pred_name_long") or self.name)
        thr = float(row["threshold"])
        self.threshold = -thr if self.inverse else thr
        self.score_col = str(row["score_col"])
        self.color = str(row.get("color") or "")
        self.gene_universe = str(row.get("geneUniverse") or "") or None
        self.pred_files = pred_files            # biosample -> path
        self.matches = matches                  # (tissue, biosample)
        self.biosamples = list(pred_files)


def read_configs(methods_table: Path, predictions_table: Path, methods: "Optional[list[str]]") -> "list[Method]":
    pd = _pd()
    mt = pd.read_csv(methods_table, sep="\t", dtype=str).fillna("")
    pt = pd.read_csv(predictions_table, sep="\t", dtype=str)
    pt = pt.dropna(subset=["biosample"])
    if methods:
        mt = mt[mt["method"].isin(methods)]
        missing = set(methods) - set(mt["method"])
        if missing:
            raise SystemExit(f"methods not in the methods table: {sorted(missing)}")
    out = []
    for _, row in mt.iterrows():
        m = row["method"]
        if m not in pt.columns:
            raise SystemExit(f"predictions table has no column for method {m}")
        sub = pt[["biosample", m] + (["GTExTissue"] if "GTExTissue" in pt.columns else [])].dropna(subset=[m])
        sub = sub[sub[m].astype(str).str.strip() != ""]
        files = {str(r["biosample"]).strip(): Path(str(r[m]).strip()) for _, r in sub.iterrows()}
        matches = []
        if "GTExTissue" in sub.columns:
            for _, r in sub.dropna(subset=["GTExTissue"]).iterrows():
                for t in str(r["GTExTissue"]).split(","):
                    if t.strip():
                        matches.append((t.strip(), str(r["biosample"]).strip()))
        out.append(Method(row.to_dict(), files, matches))
    return out


def per_method(method: Method, eq, bg, per_tissue, n_bg_total: int, dmin, dmax, tss, universe_default: "set[str]",
               n_steps: int, alpha: float, compat: bool, out: Path) -> dict:
    """All per-method rules: intersections, threshold span, count matrices, enrichment and recall tables,
    enrichment-recall curves, enhancer set sizes."""
    pd, np = _pd(), _np()
    mdir = out / method.name
    mdir.mkdir(parents=True, exist_ok=True)
    universe = universe_default
    if method.gene_universe:
        p = Path(method.gene_universe)
        if not p.is_absolute():
            p = ROOT / p
        if p.is_file():
            universe = set(load_tss(p)["gene"].astype(str))
    tissues = list(per_tissue["tissue"])
    n_loc = dict(zip(per_tissue["tissue"], per_tissue["n_all"]))
    n_loc_bin = {mx: dict(zip(per_tissue["tissue"], per_tissue[f"n_bin{i}"])) for i, mx in enumerate(dmax, start=1)}
    # pair totals per tissue (recall denominators), overall and per bin
    pairs = eq[["chr", "start", "end", "gene", "tissue", "distance_bin"]].drop_duplicates(["chr", "start", "end", "gene", "tissue"])
    n_pairs = pairs.groupby("tissue").size().to_dict()
    n_pairs_bin = {mx: pairs[pairs["distance_bin"] == mx].groupby("tissue").size().to_dict() for mx in dmax}

    var_ints, bg_max, elements, n_pred_rows = {}, {}, {}, {}
    for bs, path in method.pred_files.items():
        if not path.is_file():
            raise SystemExit(f"{method.name}/{bs}: prediction file not found: {path}")
        pred = read_predictions(path, method.score_col, universe, method.inverse, bs)
        n_pred_rows[bs] = len(pred)
        el = element_max_scores(pred)
        elements[bs] = el
        vi = overlap_join(eq, pred[["chr", "start", "end", "TargetGene", "score"]], suffixes=("", "_enh"))
        var_ints[bs] = vi
        # background: best element score per rsID
        bi = overlap_join(bg, el, suffixes=("", "_enh"))
        bg_max[bs] = bi.groupby("rsid")["score"].max().to_numpy() if len(bi) else np.array([], dtype=float)
        log.info("%s/%s: %d prediction rows, %d elements, %d variant overlaps, %d background SNPs overlapping",
                 method.name, bs, len(pred), len(el), len(vi), bg_max[bs].size)
    all_int = pd.concat(var_ints.values(), ignore_index=True) if var_ints else pd.DataFrame(columns=EQTL_COLS + ["distance_bin", "TargetGene", "score", "biosample"])
    all_int["biosample"] = all_int["biosample"] if "biosample" in all_int.columns else ""
    for bs, vi in var_ints.items():
        vi["biosample"] = bs
    all_int = pd.concat(var_ints.values(), ignore_index=True) if var_ints else all_int

    thresholds = threshold_span(all_int, method.matches, n_steps, method.boolean, method.threshold)
    pd.DataFrame({"threshold": thresholds}).to_csv(mdir / "threshold_span.tsv", sep="\t", index=False, header=False)

    # --- per biosample x tissue: best score per variant location, per pair (any gene), per pair (correct gene)
    def best_scores(vi, dist_bin=None):
        sub = vi if dist_bin is None else vi[vi["distance_bin"] == dist_bin]
        loc = sub.groupby(["tissue", "chr", "start", "end"], sort=False)["score"].max().reset_index()
        pair = sub.groupby(["tissue", "chr", "start", "end", "gene"], sort=False)["score"].max().reset_index()
        link = sub[sub["gene"] == sub["TargetGene"]].groupby(["tissue", "chr", "start", "end", "gene"], sort=False)["score"].max().reset_index()
        return loc, pair, link

    match_keys = {f"{t}.{b}" for t, b in method.matches}
    enr_rows, rec_rows = [], []
    bg_counts = {bs: counts_at_thresholds(bg_max[bs], thresholds) for bs in method.biosamples}
    for bs in method.biosamples:
        loc, pair, link = best_scores(var_ints[bs])
        for tissue in tissues:
            l_s = loc.loc[loc["tissue"] == tissue, "score"].to_numpy()
            p_s = pair.loc[pair["tissue"] == tissue, "score"].to_numpy()
            k_s = link.loc[link["tissue"] == tissue, "score"].to_numpy()
            c_loc = counts_at_thresholds(l_s, thresholds)
            c_pair = counts_at_thresholds(p_s, thresholds)
            c_link = counts_at_thresholds(k_s, thresholds)
            for i, t in enumerate(thresholds):
                enr_rows.append({"Biosample": bs, "threshold": t, "GTExTissue": tissue, "nVariantsOverlappingEnhancers": c_loc[i],
                                 "nVariantsGTExTissue": int(n_loc.get(tissue, 0)), "nCommonVariantsOverlappingEnhancers": bg_counts[bs][i],
                                 "nCommonVariants": n_bg_total})
                rec_rows.append({"Biosample": bs, "threshold": t, "GTExTissue": tissue, "nVariantGenePairsOverlappingEnhancers": c_pair[i],
                                 "nVariantsOverlappingEnhancersCorrectGene": c_link[i], "total.variants": int(n_pairs.get(tissue, 0))})
    enr = pd.DataFrame(enr_rows)
    rec = pd.DataFrame(rec_rows)
    # aggregate over matched pairs
    if len(enr):
        enr["key"] = enr["GTExTissue"] + "." + enr["Biosample"]
        agg = enr[enr["key"].isin(match_keys)].groupby("threshold", as_index=False)[
            ["nVariantsGTExTissue", "nVariantsOverlappingEnhancers", "nCommonVariants", "nCommonVariantsOverlappingEnhancers"]].sum()
        agg["GTExTissue"] = "all_matches"; agg["Biosample"] = "all_matches"
        enr = pd.concat([enr.drop(columns=["key"]), agg], ignore_index=True)
        e, lo, hi, se, p = enrichment_stats(enr["nVariantsOverlappingEnhancers"], enr["nVariantsGTExTissue"],
                                            enr["nCommonVariantsOverlappingEnhancers"], enr["nCommonVariants"], alpha, compat)
        enr["enrichment"], enr["CI_enr_low"], enr["CI_enr_high"], enr["SE_log_enr"], enr["p_adjust_enr"] = e, lo, hi, se, bonferroni(p)
        enr["method"] = method.name
        rec["key"] = rec["GTExTissue"] + "." + rec["Biosample"]
        ragg = rec[rec["key"].isin(match_keys)].groupby("threshold", as_index=False)[
            ["nVariantGenePairsOverlappingEnhancers", "nVariantsOverlappingEnhancersCorrectGene", "total.variants"]].sum()
        ragg["GTExTissue"] = "all_matches"; ragg["Biosample"] = "all_matches"
        rec = pd.concat([rec.drop(columns=["key"]), ragg], ignore_index=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            rec["recall.total"] = rec["nVariantGenePairsOverlappingEnhancers"] / rec["total.variants"]
            rec["recall.linking"] = rec["nVariantsOverlappingEnhancersCorrectGene"] / rec["total.variants"]
            rec["correctGene.ifOverlap"] = rec["nVariantsOverlappingEnhancersCorrectGene"] / rec["nVariantGenePairsOverlappingEnhancers"]
        rec["method"] = method.name
    write_tsv(enr, mdir / "enrichment_across_thresholds.tsv.gz")
    write_tsv(rec, mdir / "recall_across_thresholds.tsv.gz")

    # --- by distance, at the method's own threshold
    thr = method.threshold
    enr_d_rows, rec_d_rows = [], []
    for bs in method.biosamples:
        vi = var_ints[bs]
        bg_all = int(bg_max[bs].size) if compat else int((bg_max[bs] >= thr).sum())
        for mn, mx in zip(dmin, dmax):
            is_all = mx == ALL_BIN_KB
            sub = vi if is_all else vi[vi["distance_bin"] == mx]
            thresholded = sub[sub["score"] >= thr]
            # upstream quirk: the "all" bin counts every overlap regardless of score
            enr_src = sub if (compat and is_all) else thresholded
            for tissue in tissues:
                s_e = enr_src[enr_src["tissue"] == tissue]
                s_r = thresholded[thresholded["tissue"] == tissue]
                n_over = len(s_e[["chr", "start", "end"]].drop_duplicates())
                n_pair_over = len(s_r[["chr", "start", "end", "gene"]].drop_duplicates())
                n_link = len(s_r[s_r["gene"] == s_r["TargetGene"]][["chr", "start", "end", "gene"]].drop_duplicates())
                enr_d_rows.append({"Biosample": bs, "GTExTissue": tissue, "nVariantsOverlappingEnhancers": n_over,
                                   "nVariantsGTExTissue": int(n_loc.get(tissue, 0)) if is_all else int(n_loc_bin[mx].get(tissue, 0)),
                                   "nCommonVariantsOverlappingEnhancers": bg_all, "nCommonVariants": n_bg_total,
                                   "distance_min": mn, "distance_max": mx})
                rec_d_rows.append({"GTExTissue": tissue, "nVariantGenePairsOverlappingEnhancers": n_pair_over,
                                   "nVariantsOverlappingEnhancersCorrectGene": n_link, "Biosample": bs, "distance_min": mn, "distance_max": mx,
                                   "total.variants": int(n_pairs.get(tissue, 0)) if is_all else int(n_pairs_bin[mx].get(tissue, 0))})
    enr_d = pd.DataFrame(enr_d_rows)
    rec_d = pd.DataFrame(rec_d_rows)
    if len(enr_d):
        parts = []
        for (mn, mx), sub in enr_d.groupby(["distance_min", "distance_max"], sort=False):
            sub = sub.copy()
            e, lo, hi, se, p = enrichment_stats(sub["nVariantsOverlappingEnhancers"], sub["nVariantsGTExTissue"],
                                                sub["nCommonVariantsOverlappingEnhancers"], sub["nCommonVariants"], alpha, compat)
            sub["enrichment"], sub["CI_enr_low"], sub["CI_enr_high"], sub["SE_log_enr"], sub["p_adjust_enr"] = e, lo, hi, se, bonferroni(p)
            parts.append(sub)
        enr_d = pd.concat(parts, ignore_index=True)
        enr_d["method"] = method.name
        with np.errstate(divide="ignore", invalid="ignore"):
            rec_d["recall.total"] = rec_d["nVariantGenePairsOverlappingEnhancers"] / rec_d["total.variants"]
            rec_d["recall.linking"] = rec_d["nVariantsOverlappingEnhancersCorrectGene"] / rec_d["total.variants"]
            rec_d["correctGene.ifOverlap"] = rec_d["nVariantsOverlappingEnhancersCorrectGene"] / rec_d["nVariantGenePairsOverlappingEnhancers"]
        rec_d["method"] = method.name
    write_tsv(enr_d, mdir / "enrichment_by_distance.tsv")
    write_tsv(rec_d, mdir / "recall_by_distance.tsv")

    # --- enrichment-recall curves: one per matched pair + AllMatches
    er_parts = []
    if len(enr) and len(rec):
        e_cols = ["GTExTissue", "Biosample", "threshold", "enrichment", "CI_enr_low", "CI_enr_high", "SE_log_enr", "p_adjust_enr"]
        for tissue, bs in method.matches + [("all_matches", "all_matches")]:
            r = rec[(rec["GTExTissue"] == tissue) & (rec["Biosample"] == bs)]
            e = enr[(enr["GTExTissue"] == tissue) & (enr["Biosample"] == bs)][e_cols]
            if r.empty:
                continue
            df = r.merge(e, on=["threshold", "GTExTissue", "Biosample"], how="left")
            df["curve"] = "AllMatches" if tissue == "all_matches" else f"GTExTissue{tissue}.Biosample{bs}"
            df["nPoints"] = len(df)
            er_parts.append(df)
    er = pd.concat(er_parts, ignore_index=True) if er_parts else pd.DataFrame()
    write_tsv(er, mdir / "enrichment_recall_curves.tsv")

    # --- enhancer set sizes at the threshold, elements overlapping the TSS reference excluded
    size_rows = []
    for bs in method.biosamples:
        el = elements[bs]
        el = el[el["score"] >= thr][["chr", "start", "end"]].drop_duplicates()
        if len(el):
            hit = overlap_join(el.reset_index(drop=True), tss[["chr", "start", "end"]], suffixes=("", "_tss"))
            hit_keys = set(zip(hit["chr"], hit["start"], hit["end"])) if len(hit) else set()
            keep = el[[k not in hit_keys for k in zip(el["chr"], el["start"], el["end"])]]
            bp = merged_bp(keep)
        else:
            bp = 0
        size_rows.append({"biosample": bs, "bp": bp, "n_elements_at_threshold": int(len(el))})
    sizes = pd.DataFrame(size_rows)
    write_tsv(sizes, mdir / "enhancer_set_sizes.tsv")

    return {"method": method.name, "thresholds": thresholds, "n_thresholds": len(thresholds), "threshold_used": thr,
            "biosamples": method.biosamples, "matches": method.matches, "n_prediction_rows": n_pred_rows,
            "n_elements": {bs: int(len(el)) for bs, el in elements.items()},
            "n_variant_overlaps": {bs: int(len(vi)) for bs, vi in var_ints.items()},
            "n_background_overlapping": {bs: int(v.size) for bs, v in bg_max.items()},
            "enrichment_across_thresholds": enr, "recall_across_thresholds": rec, "enrichment_by_distance": enr_d,
            "recall_by_distance": rec_d, "er_curves": er, "enhancer_set_sizes": sizes}


def enrichment_at_recall(er_tables: "list[Any]", target: float, alpha: float):
    """plot_enrichment_with_ci_at_recall.R: closest recall.linking within 10% of the target per curve; pairwise z-tests."""
    pd, np = _pd(), _np()
    from scipy.stats import norm  # type: ignore
    picked = []
    for df in er_tables:
        df = df.dropna(subset=["recall.linking", "enrichment"])
        df = df[np.isfinite(df["recall.linking"])]
        if df.empty:
            continue
        diff = (df["recall.linking"] - target).abs()
        if diff.min() <= target * 0.1:
            picked.append(df.loc[[diff.idxmin()]])
    if not picked:
        return pd.DataFrame(), pd.DataFrame()
    sel = pd.concat(picked, ignore_index=True)
    sel["key"] = sel["method"] + " (" + sel["Biosample"].astype(str) + ")"
    comp_rows = []
    for a, b in combinations(range(len(sel)), 2):
        d = math.log(sel.loc[a, "enrichment"]) - math.log(sel.loc[b, "enrichment"])
        se = math.sqrt(sel.loc[a, "SE_log_enr"] ** 2 + sel.loc[b, "SE_log_enr"] ** 2)
        z = d / se if se > 0 else float("nan")
        p = 2 * norm.cdf(-abs(z)) if math.isfinite(z) else float("nan")
        comp_rows.append({"group1": sel.loc[a, "key"], "group2": sel.loc[b, "key"], "p": p})
    comp = pd.DataFrame(comp_rows)
    if len(comp):
        comp["p_adjust"] = bonferroni(comp["p"])
        comp["significant"] = comp["p_adjust"] < alpha
    sel = sel.sort_values("enrichment", ascending=False).reset_index(drop=True)
    return sel, comp


def cluster_order(M) -> "tuple[list, list]":
    """Row and column orders from Ward clustering of 1 - Pearson correlation (plot_*_heatmap.R)."""
    np = _np()
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage  # type: ignore
        from scipy.spatial.distance import pdist  # type: ignore
    except Exception:
        return list(M.index), list(M.columns)

    def order(X, labels):
        if X.shape[0] < 3:
            return list(labels)
        with np.errstate(invalid="ignore"):
            C = np.corrcoef(X)
        C = np.nan_to_num(C, nan=0.0)
        D = 1 - C
        try:
            Z = linkage(pdist(D), method="ward")
            return [labels[i] for i in leaves_list(Z)]
        except Exception:
            return list(labels)
    X = M.to_numpy(dtype=float)
    return order(X, list(M.index)), order(X.T, list(M.columns))


def make_figures(results: "list[dict]", methods: "list[Method]", per_tissue, dmin, dmax, matched_tissues: "list[str]",
                 er_by_tissue: "dict[str, Any]", ear: "dict[tuple[str, float], Any]", all_enr, all_pred, heat: "dict[str, tuple]",
                 fig_dir: Path, alpha: float) -> "list[Path]":
    plt = _plt()
    if plt is None:
        print("matplotlib not available; figures skipped")
        return []
    pd, np = _pd(), _np()
    fig_dir.mkdir(parents=True, exist_ok=True)
    figs = []
    color = {m.name: (m.color if m.color.startswith("#") else SERIES[i % len(SERIES)]) for i, m in enumerate(methods)}

    # 1. variants per tissue by distance bin
    bins = [(mn, mx) for mn, mx in zip(dmin, dmax) if mx != ALL_BIN_KB]
    if len(per_tissue):
        pt = per_tissue.sort_values("n_all")
        fig, ax = plt.subplots(figsize=(7, max(3.5, 0.28 * len(pt) + 1)))
        left = np.zeros(len(pt))
        labels = [f"{mn}-{mx} kb" for mn, mx in bins] + [f"> {bins[-1][1]} kb" if bins else "all"]
        cols = [f"n_bin{i}" for i in range(1, len(dmax) + 1)]
        for k, (col, lab) in enumerate(zip(cols, labels)):
            vals = pt[col].to_numpy(dtype=float)
            ax.barh(pt["tissue"], vals, left=left, color=DIST_COLORS[k % len(DIST_COLORS)], label=lab)
            left += vals
        ax.axvline(50, ls="--", color="#96a0b3", lw=0.8)
        ax.legend(fontsize=7, frameon=False, title="eVariant-eGene distance", title_fontsize=7)
        _style(ax, "Fine-mapped eQTL variants above the PIP threshold per tissue (distal noncoding)", "number of variant locations")
        figs.append(_save(fig, fig_dir / "fig1_variants_per_tissue.png"))

    # 2. thresholded performance by distance (matched pairs): enrichment, recall total, correct gene if overlap
    if all_enr is not None and len(all_enr) and all_pred is not None and len(all_pred):
        enr_plot = all_enr[(all_enr["nVariantsGTExTissue"] > 20) & np.isfinite(all_enr["enrichment"]) & (all_enr["nVariantsOverlappingEnhancers"] >= 5)]
        pred_plot = all_pred[all_pred["total.variants"] > 20]
        fig, axes = plt.subplots(1, 3, figsize=(13, max(3.5, 0.55 * len(dmax) * len(methods) / 2 + 1.5)), sharey=True)
        dist_labels = [("All variants" if mx == ALL_BIN_KB else f"{mn}-{mx} kb") for mn, mx in zip(dmin, dmax)]
        panels = [(enr_plot, "enrichment", "Enrichment (eQTLs vs common variants)"),
                  (pred_plot, "recall.total", "Variants overlapping predicted enhancers"),
                  (pred_plot, "correctGene.ifOverlap", "Linked to the correct gene, given overlap")]
        nm = len(methods)
        for ax, (df, col, title) in zip(axes, panels):
            for j, m in enumerate(methods):
                data = []
                for mn, mx in zip(dmin, dmax):
                    sub = df[(df["method"] == m.name) & (df["distance_max"] == mx) & (df["distance_min"] == mn)][col]
                    data.append(sub[np.isfinite(sub)].to_numpy())
                pos = np.arange(len(dmax)) * (nm + 1) + j
                bp = ax.boxplot(data, positions=pos, vert=False, widths=0.8, patch_artist=True, showfliers=True,
                                flierprops={"markersize": 2}, medianprops={"color": INK})
                for patch in bp["boxes"]:
                    patch.set_facecolor(color[m.name]); patch.set_alpha(0.85)
            ax.set_yticks(np.arange(len(dmax)) * (nm + 1) + (nm - 1) / 2)
            ax.set_yticklabels(dist_labels)
            _style(ax, title, col)
        handles = [plt.Rectangle((0, 0), 1, 1, color=color[m.name]) for m in methods]
        axes[-1].legend(handles, [m.long for m in methods], fontsize=7, frameon=False, loc="lower right")
        fig.suptitle("Matched eQTL tissue / prediction biosample pairs at each method's threshold, stratified by eVariant-eGene distance",
                     fontsize=10, x=0.01, ha="left")
        figs.append(_save(fig, fig_dir / "fig2_thresholded_performance_by_distance.png"))

        # 3. scatter: recall (linking) vs log10 enrichment for the all-variants bin
        e_all = enr_plot[enr_plot["distance_max"] == ALL_BIN_KB]
        p_all = pred_plot[pred_plot["distance_max"] == ALL_BIN_KB]
        both = e_all.merge(p_all, on=["method", "GTExTissue", "Biosample", "distance_min", "distance_max"], how="inner")
        if len(both):
            fig, ax = plt.subplots(figsize=(6.5, 4.5))
            for m in methods:
                sub = both[both["method"] == m.name]
                ax.scatter(sub["recall.linking"], np.log10(sub["enrichment"]), s=28, alpha=0.8, color=color[m.name], label=m.long)
            ax.legend(fontsize=7, frameon=False)
            _style(ax, "Matched pairs at threshold: recall (linking) vs enrichment", "recall (variants overlapping an enhancer linked to the eGene)",
                   "log10 enrichment (eQTLs vs common variants)")
            figs.append(_save(fig, fig_dir / "fig3_matched_recall_vs_enrichment.png"))

    # 4. enrichment-recall curves per tissue + AllMatches
    for tissue, df in er_by_tissue.items():
        if df is None or df.empty:
            continue
        df = df[(df["total.variants"] > 20) & (df["recall.total"] > 0.001)]
        if df.empty:
            continue
        fig, ax = plt.subplots(figsize=(7, 4.6))
        ymax = 0
        for k, (curve, sub) in enumerate(df.groupby(["method", "curve"], sort=False)):
            m = next((mm for mm in methods if mm.name == curve[0]), None)
            sub = sub.sort_values("threshold")
            sub_plot = sub[sub["recall.linking"] > 0]
            if sub_plot.empty:
                continue
            c = color.get(curve[0], SERIES[k % len(SERIES)])
            label = f"{m.long if m else curve[0]}" + ("" if tissue == "AllMatches" else f" ({sub['Biosample'].iloc[0]})")
            if int(sub["nPoints"].iloc[0]) == 2:
                pt = sub_plot[sub_plot["threshold"] == 1]
                ax.scatter(pt["recall.linking"], pt["enrichment"], color=c, s=40, label=label)
            else:
                ax.plot(sub_plot["recall.linking"], sub_plot["enrichment"], color=c, lw=1.4, label=label)
            if m is not None:
                at = sub[np.isclose(sub["threshold"], m.threshold)]
                if len(at):
                    ax.errorbar(at["recall.linking"], at["enrichment"], yerr=[(at["enrichment"] - at["CI_enr_low"]).clip(lower=0),
                                                                              (at["CI_enr_high"] - at["enrichment"]).clip(lower=0)],
                                fmt="o", color=c, ms=5, capsize=2, lw=1)
            ymax = max(ymax, float(np.nanmax(sub_plot["enrichment"])) if len(sub_plot) else 0)
        ax.set_ylim(0, min(50, ymax * 1.05 if ymax > 0 else 1))
        ax.legend(fontsize=7, frameon=False, ncol=1 + len(df["curve"].unique()) // 12)
        n_var = df["total.variants"].mean()
        _style(ax, f"Enrichment-recall: eQTL tissue {tissue}" if tissue != "AllMatches" else "Enrichment-recall aggregated over all matched pairs",
               f"recall (variants overlapping a prediction linked to the eGene); {n_var:.0f} variant" + ("-biosample pairs" if tissue == "AllMatches" else "s"),
               "enrichment (eQTLs vs common variants)")
        figs.append(_save(fig, fig_dir / f"fig4_enrichment_recall.{safe_label(tissue)}.png"))

    # 5. enrichment at recall
    for (tissue, target), sel in ear.items():
        if sel is None or sel.empty or len(sel) < 2:
            continue
        fig, ax = plt.subplots(figsize=(max(4, 0.8 * len(sel) + 2), 4))
        cols = [color.get(m, SERIES[i % len(SERIES)]) for i, m in enumerate(sel["method"])]
        ax.bar(range(len(sel)), sel["enrichment"], color=cols)
        ax.errorbar(range(len(sel)), sel["enrichment"], yerr=[(sel["enrichment"] - sel["CI_enr_low"]).clip(lower=0),
                                                            (sel["CI_enr_high"] - sel["enrichment"]).clip(lower=0)], fmt="none", ecolor=INK2, capsize=3)
        ax.set_xticks(range(len(sel)))
        ax.set_xticklabels([f"{k}\n(recall {r:.3f})" for k, r in zip(sel["key"], sel["recall.linking"])], fontsize=7, rotation=30, ha="right")
        _style(ax, f"Enrichment at recall {target} — {tissue}", "", "enrichment (eQTLs vs common variants)")
        figs.append(_save(fig, fig_dir / f"fig5_enrichment_at_recall{target}.{safe_label(tissue)}.png"))

    # 6. heatmaps
    for name, (Me, Mr, sizes, nvar) in heat.items():
        for kind, M, lims in (("enrichment", Me, None), ("recall", Mr, (0, 0.25))):
            if M is None or M.empty:
                continue
            ro, co = cluster_order(M)
            M = M.loc[ro, co]
            fig, ax = plt.subplots(figsize=(max(5, 0.32 * M.shape[1] + 2), max(3.5, 0.3 * M.shape[0] + 1.5)))
            if lims is None:
                vals = M.to_numpy(dtype=float)
                vmax = max(2.0, float(np.nanquantile(vals[np.isfinite(vals)], 0.9)) if np.isfinite(vals).any() else 2.0)
                lims = (0, round(vmax, 1))
            im = ax.imshow(M.to_numpy(dtype=float), aspect="auto", cmap="PuBuGn", vmin=lims[0], vmax=lims[1])
            ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=7)
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label=kind)
            ax.set_title(f"{name}: {kind} per eQTL tissue (columns) x prediction biosample (rows)", loc="left", fontsize=10, fontweight="bold")
            figs.append(_save(fig, fig_dir / f"fig6_{kind}_heatmap.{safe_label(name)}.png"))
    return figs


def cmd_run(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    t0 = time.time()
    out = run_dir(args.label or "eqtl_benchmark")
    methods = read_configs(Path(args.methods_table), Path(args.predictions_table), args.methods)
    if not methods:
        raise SystemExit("no methods to run")
    tss_path = resource("CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6", args.tss)
    tss = load_tss(tss_path)
    universe = set(tss["gene"].astype(str))
    # variants
    if args.variants_dir:
        vdir = Path(args.variants_dir)
        eq, bg, per_tissue, vsumm = load_variants_dir(vdir)
        print(f"Variants: reused {vdir}")
    else:
        if not (args.eqtl and args.bg_variants):
            raise SystemExit("--eqtl and --bg-variants are required unless --variants-dir is given")
        vdir = out / "variants"
        vsumm = prepare_variants(Path(args.eqtl), Path(args.bg_variants), resource("PartitionCombined.bed", args.partition), tss_path,
                                 args.threshold_pip, args.distances, vdir)
        eq, bg, per_tissue, vsumm = load_variants_dir(vdir)
    dmin, dmax = vsumm["distances_min"], vsumm["distances_max"]
    n_bg_total = int(vsumm["background_snps_distal_noncoding"])
    matched_tissues = sorted({t for m in methods for t, _ in m.matches})
    matched_eq = eq[eq["tissue"].isin(matched_tissues)]
    counts_txt = (f"Number of unique eVariants: {len(matched_eq[['chr', 'start', 'end']].drop_duplicates())}\n"
                  f"Number of unique eVariant-eGene pairs: {len(matched_eq[['chr', 'start', 'end', 'gene']].drop_duplicates())}\n")
    (out / "unique_variant_counts.matched_tissues.txt").write_text(counts_txt)

    results = []
    for m in methods:
        print(f"\n== {m.name} ({m.long}): {len(m.biosamples)} biosample(s), {len(m.matches)} matched pair(s), threshold {m.threshold:g}"
              f"{' (negated: inverse predictor)' if m.inverse else ''}{' [binary]' if m.boolean else ''}")
        r = per_method(m, eq, bg, per_tissue, n_bg_total, dmin, dmax, tss, universe, args.n_threshold_steps, args.threshold_pval,
                       args.upstream_compat, out)
        results.append(r)
        er = r["er_curves"]
        if len(er):
            am = er[er["curve"] == "AllMatches"]
            at = am[np.isclose(am["threshold"], m.threshold)]
            if len(at):
                print(f"   all matches at threshold: enrichment {_fmt(at['enrichment'].iloc[0])} "
                      f"[{_fmt(at['CI_enr_low'].iloc[0])}, {_fmt(at['CI_enr_high'].iloc[0])}], recall total {_fmt(at['recall.total'].iloc[0])}, "
                      f"recall linking {_fmt(at['recall.linking'].iloc[0])}, {int(at['total.variants'].iloc[0]):,} variant-biosample pairs")

    # cross-method tables
    def matched(df, m):
        keys = {f"{t}.{b}" for t, b in m.matches}
        return df[(df["GTExTissue"] + "." + df["Biosample"]).isin(keys)]
    all_enr = pd.concat([matched(r["enrichment_by_distance"], m) for r, m in zip(results, methods) if len(r["enrichment_by_distance"])],
                        ignore_index=True) if results else pd.DataFrame()
    all_pred = pd.concat([matched(r["recall_by_distance"], m) for r, m in zip(results, methods) if len(r["recall_by_distance"])],
                         ignore_index=True) if results else pd.DataFrame()
    write_tsv(all_enr, out / "all_matched_enrichments.tsv")
    write_tsv(all_pred, out / "all_matched_prediction_metrics.tsv")

    er_by_tissue: "dict[str, Any]" = {}
    for tissue in matched_tissues + ["AllMatches"]:
        parts = []
        for r in results:
            er = r["er_curves"]
            if not len(er):
                continue
            parts.append(er[er["curve"] == "AllMatches"] if tissue == "AllMatches" else er[er["curve"].str.startswith(f"GTExTissue{tissue}.Biosample")])
        df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        er_by_tissue[tissue] = df
        if len(df):
            write_tsv(df, out / f"enrichment_recall.{safe_label(tissue)}.tsv")

    ear: "dict[tuple[str, float], Any]" = {}
    ear_summary = []
    for tissue in matched_tissues + ["AllMatches"]:
        df = er_by_tissue.get(tissue)
        if df is None or df.empty:
            continue
        curves = [sub for _, sub in df.groupby(["method", "curve"], sort=False)]
        for target in args.recalls:
            sel, comp = enrichment_at_recall(curves, float(target), args.threshold_pval)
            ear[(tissue, float(target))] = sel
            if len(sel):
                write_tsv(sel, out / f"enrichment_at_recall{target}.{safe_label(tissue)}.tsv")
                write_tsv(comp, out / f"pairwise_comparisons.recall{target}.{safe_label(tissue)}.tsv")
                for _, row in sel.iterrows():
                    ear_summary.append({"tissue": tissue, "target_recall": float(target), "method": row["method"], "biosample": row["Biosample"],
                                        "recall_linking": float(row["recall.linking"]), "enrichment": float(row["enrichment"]),
                                        "ci_low": float(row["CI_enr_low"]), "ci_high": float(row["CI_enr_high"])})

    heat: "dict[str, tuple]" = {}
    for r, m in zip(results, methods):
        e = r["enrichment_by_distance"]
        e = e[(e["distance_max"] == ALL_BIN_KB) & (e["nVariantsGTExTissue"] > 20) & np.isfinite(e["enrichment"])] if len(e) else e
        Me = e.pivot_table(index="Biosample", columns="GTExTissue", values="enrichment") if len(e) else pd.DataFrame()
        rc = r["recall_by_distance"]
        rc = rc[(rc["distance_max"] == ALL_BIN_KB) & (rc["total.variants"] > 20)] if len(rc) else rc
        Mr = rc.pivot_table(index="Biosample", columns="GTExTissue", values="recall.linking") if len(rc) else pd.DataFrame()
        if len(Mr):
            Mr = Mr.loc[:, Mr.nunique(axis=0) > 1] if Mr.shape[0] > 1 else Mr   # drop tissues constant across biosamples
        heat[m.name] = (Me, Mr, r["enhancer_set_sizes"], per_tissue)
        if len(Me):
            Me.to_csv(out / f"heatmap_enrichment.{safe_label(m.name)}.tsv", sep="\t")
        if len(Mr):
            Mr.to_csv(out / f"heatmap_recall.{safe_label(m.name)}.tsv", sep="\t")

    figs = [] if args.no_plots else make_figures(results, methods, per_tissue, dmin, dmax, matched_tissues, er_by_tissue, ear, all_enr, all_pred,
                                                 heat, out / "figures", args.threshold_pval)

    # summary + report
    method_summaries = {}
    for r, m in zip(results, methods):
        er = r["er_curves"]
        am = er[er["curve"] == "AllMatches"] if len(er) else er
        at = am[np.isclose(am["threshold"], m.threshold)] if len(am) else am
        me = matched(r["enrichment_by_distance"], m) if len(r["enrichment_by_distance"]) else r["enrichment_by_distance"]
        me = me[me["distance_max"] == ALL_BIN_KB] if len(me) else me
        mp = matched(r["recall_by_distance"], m) if len(r["recall_by_distance"]) else r["recall_by_distance"]
        mp = mp[mp["distance_max"] == ALL_BIN_KB] if len(mp) else mp
        method_summaries[m.name] = {
            "pred_name_long": m.long, "threshold": m.threshold, "inverse_predictor": m.inverse, "boolean": m.boolean,
            "biosamples": m.biosamples, "matches": m.matches, "n_thresholds": r["n_thresholds"],
            "n_prediction_rows": r["n_prediction_rows"], "n_elements": r["n_elements"],
            "enhancer_set_bp": {row["biosample"]: int(row["bp"]) for _, row in r["enhancer_set_sizes"].iterrows()},
            "all_matches_at_threshold": ({"enrichment": float(at["enrichment"].iloc[0]), "ci_low": float(at["CI_enr_low"].iloc[0]),
                                          "ci_high": float(at["CI_enr_high"].iloc[0]), "p_adjust": float(at["p_adjust_enr"].iloc[0]),
                                          "recall_total": float(at["recall.total"].iloc[0]), "recall_linking": float(at["recall.linking"].iloc[0]),
                                          "correct_gene_if_overlap": float(at["correctGene.ifOverlap"].iloc[0]),
                                          "total_variants": int(at["total.variants"].iloc[0])} if len(at) else None),
            "matched_pairs_at_threshold": {"n_pairs": int(len(me)),
                                           "median_enrichment": float(me["enrichment"].median()) if len(me) else None,
                                           "median_recall_linking": float(mp["recall.linking"].median()) if len(mp) else None,
                                           "median_recall_total": float(mp["recall.total"].median()) if len(mp) else None},
        }
    summary = {"label": args.label, "run_dir": str(out), "upstream": f"{UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}", "upstream_compat": bool(args.upstream_compat),
               "parameters": {"threshold_pip": vsumm["pip_threshold"], "threshold_pval": args.threshold_pval, "n_threshold_steps": args.n_threshold_steps,
                              "distances_kb": vsumm["distances_kb"], "recalls": [float(x) for x in args.recalls]},
               "variants": {k: v for k, v in vsumm.items() if k not in ("variants_file", "background_file")},
               "matched_tissues": matched_tissues, "methods": method_summaries, "enrichment_at_recall": ear_summary,
               "figures": [str(f) for f in figs], "seconds": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"JSON: {out / 'summary.json'}")

    lines = [f"# eQTL enrichment benchmark: {args.label or 'run'}", "",
             f"Port of [{UPSTREAM_REPO}](https://github.com/{UPSTREAM_REPO}) @ {UPSTREAM_COMMIT[:7]}"
             + (" — **upstream-compat mode** (reproduces the SE and by-distance counting quirks)" if args.upstream_compat else ""), "",
             "## Variants", "",
             f"- eQTL rows {vsumm['eqtl_rows']:,} → PIP ≥ {vsumm['pip_threshold']}: {vsumm['after_pip']:,} → distal noncoding: "
             f"{vsumm['after_distal_noncoding']:,} → gene universe: {vsumm['after_gene_universe']:,} rows "
             f"({vsumm['unique_evariants']:,} eVariants, {vsumm['unique_evariant_egene_pairs']:,} eVariant–eGene pairs, {vsumm['tissues']} tissues)",
             f"- background SNPs {vsumm['background_snps']:,} → distal noncoding: {vsumm['background_snps_distal_noncoding']:,}",
             f"- matched eQTL tissues: {', '.join(matched_tissues) or 'none'}; " + counts_txt.strip().replace("\n", "; "), "",
             "## Aggregated over matched pairs, at each method's threshold", "",
             md_table(["method", "threshold", "enrichment [95% CI]", "recall (total)", "recall (linking)", "correct gene | overlap", "pairs", "n thresholds"],
                      [[m.long, f"{m.threshold:g}",
                        (f"{_fmt(s['all_matches_at_threshold']['enrichment'])} [{_fmt(s['all_matches_at_threshold']['ci_low'])}, {_fmt(s['all_matches_at_threshold']['ci_high'])}]"
                         if s["all_matches_at_threshold"] else "—"),
                        _fmt(s["all_matches_at_threshold"]["recall_total"]) if s["all_matches_at_threshold"] else "—",
                        _fmt(s["all_matches_at_threshold"]["recall_linking"]) if s["all_matches_at_threshold"] else "—",
                        _fmt(s["all_matches_at_threshold"]["correct_gene_if_overlap"]) if s["all_matches_at_threshold"] else "—",
                        s["matched_pairs_at_threshold"]["n_pairs"], s["n_thresholds"]]
                       for m, s in ((m, method_summaries[m.name]) for m in methods)]), ""]
    if ear_summary:
        lines += ["## Enrichment at target recall (linking)", "",
                  md_table(["tissue", "target", "method", "biosample", "recall", "enrichment [95% CI]"],
                           [[e["tissue"], e["target_recall"], e["method"], e["biosample"], _fmt(e["recall_linking"]),
                             f"{_fmt(e['enrichment'])} [{_fmt(e['ci_low'])}, {_fmt(e['ci_high'])}]"] for e in ear_summary]), ""]
    if len(all_enr):
        lines += ["## Matched pairs by eVariant–eGene distance (median over pairs with > 20 variants)", ""]
        rows = []
        for m in methods:
            for mn, mx in zip(dmin, dmax):
                e = all_enr[(all_enr["method"] == m.name) & (all_enr["distance_max"] == mx) & (all_enr["distance_min"] == mn) & (all_enr["nVariantsGTExTissue"] > 20)]
                p = all_pred[(all_pred["method"] == m.name) & (all_pred["distance_max"] == mx) & (all_pred["distance_min"] == mn) & (all_pred["total.variants"] > 20)]
                rows.append([m.long, "all" if mx == ALL_BIN_KB else f"{mn}-{mx} kb", len(e), _fmt(e["enrichment"].median()) if len(e) else "—",
                             _fmt(p["recall.total"].median()) if len(p) else "—", _fmt(p["recall.linking"].median()) if len(p) else "—",
                             _fmt(p["correctGene.ifOverlap"].median()) if len(p) else "—"])
        lines += [md_table(["method", "distance", "pairs", "enrichment", "recall (total)", "recall (linking)", "correct gene | overlap"], rows), ""]
    lines += ["## Enhancer set sizes at threshold (bp, elements overlapping the TSS reference excluded)", "",
              md_table(["method", "biosample", "bp", "elements"],
                       [[m.name, row["biosample"], f"{int(row['bp']):,}", int(row["n_elements_at_threshold"])]
                        for r, m in zip(results, methods) for _, row in r["enhancer_set_sizes"].iterrows()]), "",
              "## Definitions", "",
              "- enrichment = (eQTL variant locations overlapping predicted enhancers with score ≥ t / all eQTL variant locations in the tissue) ÷ "
              "(background SNPs overlapping / all background SNPs); both sets restricted to distal noncoding regions (partition ABC/AllPeaks/Other/OtherIntron). "
              "CI from the log risk ratio; p from the hypergeometric upper tail, Bonferroni over the table.",
              "- recall (total) = eVariant–eGene pairs whose variant overlaps a predicted enhancer ÷ all pairs in the tissue; recall (linking) requires the enhancer to be linked to the eGene.",
              "- thresholds: quantiles of the matched, correctly-linked scores plus the method's threshold and an even grid; binary predictors use {0, 1}; inverse predictors are negated.",
              "- all_matches rows sum numerators and denominators over the matched (tissue, biosample) pairs before dividing.", "",
              "## Files", "",
              "- `<method>/enrichment_across_thresholds.tsv.gz`, `<method>/recall_across_thresholds.tsv.gz`, `<method>/enrichment_by_distance.tsv`, "
              "`<method>/recall_by_distance.tsv`, `<method>/enrichment_recall_curves.tsv`, `<method>/threshold_span.tsv`, `<method>/enhancer_set_sizes.tsv`",
              "- `all_matched_enrichments.tsv`, `all_matched_prediction_metrics.tsv`, `enrichment_recall.<tissue>.tsv`, `enrichment_at_recall<r>.<tissue>.tsv`, "
              "`pairwise_comparisons.recall<r>.<tissue>.tsv`, `heatmap_enrichment.<method>.tsv`, `heatmap_recall.<method>.tsv`, `variants/`, `figures/`"]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Report: {out / 'report.md'}")
    print(f"Done in {time.time() - t0:.1f}s -> {out}")
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def synthetic_world(d: Path, seed: int = 7) -> dict:
    """Two 6-Mb chromosomes, 60 genes, 3 true enhancers per gene; eQTLs in tissues A and B sitting in the true
    enhancers of tissue-specific gene sets; background SNPs uniform; four predictors over two biosamples."""
    np, pd = _np(), _pd()
    rng = np.random.default_rng(seed)
    chroms = {"chr1": 6_000_000, "chr2": 6_000_000}
    genes = []
    for ci, (c, L) in enumerate(chroms.items()):
        for gi in range(30):
            tss = 200_000 + gi * 190_000 + int(rng.integers(0, 20_000))
            genes.append({"chr": c, "tss": tss, "gene": f"G{ci * 30 + gi:03d}", "strand": "+" if gi % 2 == 0 else "-"})
    gdf = pd.DataFrame(genes)
    tss_bed = pd.DataFrame({"chr": gdf["chr"], "start": gdf["tss"] - 250, "end": gdf["tss"] + 250, "gene": gdf["gene"], "score": 0, "strand": gdf["strand"]})
    tss_path = d / "tss.bed6"; tss_bed.to_csv(tss_path, sep="\t", header=False, index=False)
    gene_tsv = pd.DataFrame({"chr": gdf["chr"], "start": gdf["tss"] - 1000, "end": gdf["tss"] + 5000, "name": gdf["gene"], "score": 0,
                             "strand": gdf["strand"], "Ensembl_ID": [f"ENSG{100000 + i:011d}" for i in range(len(gdf))], "gene_type": "protein_coding"})
    gene_tsv_path = d / "genes.tsv"; gene_tsv.to_csv(gene_tsv_path, sep="\t", index=False)
    # partition: promoters/CDS around each TSS (excluded), everything else "Other"
    part_rows = []
    for c, L in chroms.items():
        sub = gdf[gdf["chr"] == c].sort_values("tss")
        pos = 0
        for tss in sub["tss"]:
            a, b = tss - 250, tss + 2_000
            if a > pos:
                part_rows.append((c, pos, a, "Other"))
            part_rows.append((c, a, tss + 250, "TSS-500bp"))
            part_rows.append((c, tss + 250, b, "CDS"))
            pos = b
        part_rows.append((c, pos, L, "Other"))
    part_path = d / "partition.bed"
    pd.DataFrame(part_rows).to_csv(part_path, sep="\t", header=False, index=False)
    # true enhancers: 3 per gene at 8-120 kb, 400 bp wide
    enh = []
    for _, g in gdf.iterrows():
        for k in range(3):
            off = int(rng.integers(8_000, 120_000)) * (1 if rng.random() < 0.5 else -1)
            s = max(3000, g["tss"] + off)
            enh.append({"chr": g["chr"], "start": s, "end": s + 400, "gene": g["gene"], "k": k})
    edf = pd.DataFrame(enh)
    # tissue programmes: tissue A -> genes with even index, tissue B -> odd; tissue C shares A's
    gA = set(gdf["gene"][::2]); gB = set(gdf["gene"][1::2])
    # eQTLs: for each tissue gene, 4 variants: 3 inside true enhancers (2 with high PIP), one in a promoter (excluded)
    eq_rows = []
    for tissue, gs in (("Tissue_A", gA), ("Tissue_B", gB), ("Tissue_C", gA)):
        for g in sorted(gs):
            e = edf[edf["gene"] == g]
            for j, (_, row) in enumerate(e.iterrows()):
                pos = int(rng.integers(row["start"], row["end"]))
                pip = 0.9 if j < 2 else 0.2
                eq_rows.append({"chr": row["chr"], "start": pos, "end": pos + 1, "varID_hg38": f"{row['chr']}_{pos}_A_G", "gene_hgnc": g,
                                "tissue": tissue, "pip": pip})
            tss = int(gdf.loc[gdf["gene"] == g, "tss"].iloc[0])
            eq_rows.append({"chr": row["chr"], "start": tss + 10, "end": tss + 11, "varID_hg38": f"{row['chr']}_{tss + 10}_C_T", "gene_hgnc": g,
                            "tissue": tissue, "pip": 0.95})
            # one distal variant with PIP >= 0.5 far from any enhancer (in "Other"), linked to the gene
            far = tss + 2_000_000 if tss + 2_000_000 < 5_990_000 else tss - 2_000_000
            eq_rows.append({"chr": row["chr"], "start": far, "end": far + 1, "varID_hg38": f"{row['chr']}_{far}_G_A", "gene_hgnc": g,
                            "tissue": tissue, "pip": 0.6})
    eq = pd.DataFrame(eq_rows)
    eq_path = d / "eqtl.tsv.gz"; eq.to_csv(eq_path, sep="\t", index=False)
    # background SNPs: 40,000 uniform
    bg_rows = []
    for c, L in chroms.items():
        pos = np.sort(rng.integers(0, L, size=20_000))
        bg_rows.append(pd.DataFrame({"chr": c, "start": pos, "end": pos + 1, "rsid": [f"rs{c[3:]}_{i}" for i in range(pos.size)]}))
    bg = pd.concat(bg_rows, ignore_index=True)
    bg_path = d / "bg.bed"; bg.to_csv(bg_path, sep="\t", header=False, index=False)
    # predictors per biosample: bioA knows tissue A's enhancers, bioB tissue B's
    preds = {}
    for bs, gs in (("bioA", gA), ("bioB", gB)):
        true = edf[edf["gene"].isin(gs)].copy()
        rows = []
        for _, r in true.iterrows():
            rows.append({"chr": r["chr"], "start": r["start"], "end": r["end"], "TargetGene": r["gene"], "good": float(rng.uniform(0.6, 1.0)),
                         "random": float(rng.uniform(0, 1)), "distance": float(abs((r["start"] + 200) - gdf.loc[gdf["gene"] == r["gene"], "tss"].iloc[0])),
                         "binary": 1})
            # a wrong-gene link from the same element (nearest other gene), lower score
            other = gdf[(gdf["chr"] == r["chr"]) & (gdf["gene"] != r["gene"])].iloc[int(rng.integers(0, 29))]
            rows.append({"chr": r["chr"], "start": r["start"], "end": r["end"], "TargetGene": other["gene"], "good": float(rng.uniform(0.0, 0.4)),
                         "random": float(rng.uniform(0, 1)), "distance": float(abs((r["start"] + 200) - other["tss"])), "binary": 0})
        # random decoy elements
        for c, L in chroms.items():
            for _ in range(400):
                s = int(rng.integers(3000, L - 3000))
                g = gdf[gdf["chr"] == c].iloc[int(rng.integers(0, 30))]
                rows.append({"chr": c, "start": s, "end": s + 400, "TargetGene": g["gene"], "good": float(rng.uniform(0.0, 0.5)),
                             "random": float(rng.uniform(0, 1)), "distance": float(abs(s + 200 - g["tss"])), "binary": int(rng.random() < 0.3)})
        df = pd.DataFrame(rows)
        p = d / f"pred_{bs}.tsv.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("# synthetic predictions\n")
            df.to_csv(fh, sep="\t", index=False)
        preds[bs] = p
    methods = pd.DataFrame([
        {"method": "good", "boolean": "FALSE", "inverse_predictor": "FALSE", "pred_name_long": "Good predictor", "threshold": 0.6, "score_col": "good", "color": "#792374"},
        {"method": "random", "boolean": "FALSE", "inverse_predictor": "FALSE", "pred_name_long": "Random scores", "threshold": 0.5, "score_col": "random", "color": "#96a0b3"},
        {"method": "distance", "boolean": "FALSE", "inverse_predictor": "TRUE", "pred_name_long": "Distance to TSS (inverse)", "threshold": 60000, "score_col": "distance", "color": "#0096a0"},
        {"method": "binary", "boolean": "TRUE", "inverse_predictor": "FALSE", "pred_name_long": "Binary calls", "threshold": 1, "score_col": "binary", "color": "#e0dcca"},
    ])
    mpath = d / "methods.tsv"; methods.to_csv(mpath, sep="\t", index=False)
    ptab = pd.DataFrame([{"biosample": "bioA", "good": str(preds["bioA"]), "random": str(preds["bioA"]), "distance": str(preds["bioA"]),
                          "binary": str(preds["bioA"]), "GTExTissue": "Tissue_A,Tissue_C"},
                         {"biosample": "bioB", "good": str(preds["bioB"]), "random": str(preds["bioB"]), "distance": str(preds["bioB"]),
                          "binary": str(preds["bioB"]), "GTExTissue": "Tissue_B"}])
    ppath = d / "predictions.tsv"; ptab.to_csv(ppath, sep="\t", index=False)
    # raw GTEx-style file + gct for prepare-gtex
    raw_rows = []
    for _, r in eq.head(40).iterrows():
        ens = gene_tsv.loc[gene_tsv["name"] == r["gene_hgnc"], "Ensembl_ID"].iloc[0] + ".3"
        for method, cs in (("SUSIE", 1), ("SUSIE", -1), ("DAP-G", 1)):
            raw_rows.append([r["chr"], r["start"], r["end"], r["varID_hg38"], r["varID_hg38"], "A", "G", "GTEx", method, r["tissue"], ens,
                             0.2, 0.5, 0.1, 5.0, r["pip"], cs, 0.4, 0.1])
    raw_path = d / "gtex_raw.tsv.gz"
    pd.DataFrame(raw_rows).to_csv(raw_path, sep="\t", header=False, index=False)
    gct = d / "median_tpm.gct.gz"
    with gzip.open(gct, "wt") as fh:
        fh.write("#1.2\n")
        fh.write(f"{len(gene_tsv)}\t3\n")
        fh.write("Name\tDescription\tTissue - A\tTissue_B\tTissue C\n")
        for i, (_, g) in enumerate(gene_tsv.iterrows()):
            tpm_a = 5.0 if g["name"] in gA else 0.1
            tpm_b = 5.0 if g["name"] in gB else 0.1
            fh.write(f"{g['Ensembl_ID']}\t{g['name']}\t{tpm_a}\t{tpm_b}\t{tpm_a}\n")
    return {"tss": tss_path, "genes": gene_tsv_path, "partition": part_path, "eqtl": eq_path, "bg": bg_path, "methods": mpath,
            "predictions": ppath, "raw": raw_path, "gct": gct, "n_genes": len(gdf), "gA": gA, "gB": gB, "n_bg": len(bg), "eq": eq,
            "true_enhancers": edf}


def cmd_selftest(args: argparse.Namespace) -> int:
    import shutil
    import tempfile
    pd, np = _pd(), _np()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made = []
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_world(d)
        print("\noverlap engine")
        a = pd.DataFrame({"chr": ["chr1", "chr1", "chr2", "chr1"], "start": [100, 65_530, 5, 200_000], "end": [101, 65_531, 6, 200_001]})
        b = pd.DataFrame({"chr": ["chr1", "chr1", "chr2", "chr1"], "start": [50, 65_000, 0, 100], "end": [150, 70_000, 5, 200_000], "name": ["x", "y", "z", "w"]})
        j = overlap_join(a, b)
        check(sorted(j["name"]) == ["w", "w", "x", "y"], "overlap_join: bin-spanning interval found, half-open ends respected (chr2 5 not in [0,5); 200000 not in [100,200000))")
        check(merged_bp(pd.DataFrame({"chr": ["chr1"] * 3, "start": [0, 50, 200], "end": [100, 120, 300]})) == 220, "merged_bp: union of overlapping intervals")

        print("\nstatistics")
        e, lo, hi, se, p = enrichment_stats([50], [100], [100], [10_000], 0.05)
        check(abs(e[0] - 50.0) < 1e-9 and lo[0] < 50 < hi[0] and p[0] < 1e-30, f"enrichment_stats: RR 50, CI [{lo[0]:.1f}, {hi[0]:.1f}], hypergeometric p {p[0]:.1e}")
        _, _, _, se_c, _ = enrichment_stats([50], [100], [100], [10_000], 0.05, compat=True)
        check(abs(se[0] - math.sqrt(50 / (50 * 100) + 9900 / (100 * 10_000))) < 1e-12 and se_c[0] != se[0],
              "enrichment_stats: corrected SE = sqrt(a + b); --upstream-compat reproduces sqrt(a) + b")
        check(counts_at_thresholds([0.1, 0.5, 0.5, 0.9], [0.0, 0.5, 0.95]) == [4, 3, 0], "counts_at_thresholds: score >= t, inclusive")
        thr = threshold_span(pd.DataFrame({"tissue": ["T"] * 6, "biosample": ["B"] * 6, "gene": ["g"] * 6, "TargetGene": ["g"] * 6,
                                           "chr": ["c"] * 6, "start": range(6), "end": range(1, 7), "score": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]}),
                             [("T", "B")], 5, False, 0.33)
        check(0.33 in thr and thr == sorted(thr) and len(thr) >= 5, f"threshold_span: quantiles + provided threshold + even grid ({len(thr)} thresholds)")
        check(threshold_span(pd.DataFrame(columns=["tissue", "biosample", "gene", "TargetGene", "chr", "start", "end", "score"]), [], 5, True, 1) == [0.0, 1.0],
              "threshold_span: binary predictors use {0, 1}")

        print("\nprepare-gtex")
        out_eq = d / "gtex_prepared.tsv.gz"
        rc = cmd_prepare_gtex(argparse.Namespace(raw=str(W["raw"]), method="SUSIE", keep_all_cs=False, gene_table=str(W["genes"]),
                                                 expression=str(W["gct"]), tpm=1.0, out=str(out_eq)))
        prep = pd.read_csv(out_eq, sep="\t")
        check(rc == 0 and list(prep.columns) == ["chr", "start", "end", "varID_hg38", "gene_hgnc", "tissue", "pip"], "prepare-gtex: 7-column output")
        raw_head = W["eq"].head(40)
        expected = raw_head[[(g in W["gA"]) if t in ("Tissue_A", "Tissue_C") else (g in W["gB"]) for g, t in zip(raw_head["gene_hgnc"], raw_head["tissue"])]]
        check(len(prep) == len(expected), f"prepare-gtex: SUSIE + credible set + gene universe + TPM > 1 in tissue keeps {len(prep)} of 40 (DAP-G and cs -1 rows dropped; .gct 'Tissue - A' matched to 'Tissue_A')")

        print("\nvariants")
        vdir = d / "variants"
        vs = prepare_variants(W["eqtl"], W["bg"], W["partition"], W["tss"], 0.5, [10, 100, 250, 1000], vdir)
        eq, bg, per_tissue, _ = load_variants_dir(vdir)
        n_genes_A, n_genes_B = len(W["gA"]), len(W["gB"])
        # per tissue gene: 2 high-PIP enhancer variants + 1 far variant kept; promoter variant and PIP 0.2 dropped
        check(vs["after_pip"] == (n_genes_A * 2 + n_genes_B) * 4 and vs["after_gene_universe"] == (n_genes_A * 2 + n_genes_B) * 3,
              f"variants: PIP filter keeps 4 per gene-tissue, distal-noncoding drops the promoter variant -> 3 per gene-tissue ({vs['after_gene_universe']} rows)")
        check(vs["background_snps_distal_noncoding"] < vs["background_snps"] and vs["background_snps_distal_noncoding"] > 0.9 * vs["background_snps"],
              f"variants: background SNPs in promoters/CDS removed ({vs['background_snps']} -> {vs['background_snps_distal_noncoding']})")
        check(set(eq["distance_bin"].unique()) <= {10, 100, 250, 1000, 30000} and (eq["distance_bin"] == 30000).sum() == (n_genes_A * 2 + n_genes_B),
              "variants: distance bins from the mean TSS centre; the far variants land in the 1000-30000 kb group")
        check(list(per_tissue.columns) == ["tissue", "n_all", "n_bin1", "n_bin2", "n_bin3", "n_bin4", "n_bin5"] and len(per_tissue) == 3,
              "variants: per-tissue location counts overall and per bin")

        print("\nrun")
        rc = cmd_run(argparse.Namespace(label="st_eqtl", methods_table=str(W["methods"]), predictions_table=str(W["predictions"]), methods=None,
                                        eqtl=None, bg_variants=None, variants_dir=str(vdir), partition=str(W["partition"]), tss=str(W["tss"]),
                                        threshold_pip=0.5, distances=[10, 100, 250, 1000], recalls=[0.3, 0.6], threshold_pval=0.05,
                                        n_threshold_steps=20, upstream_compat=False, no_plots=args.no_plots))
        run = sorted(OUT_ROOT.glob("*_st_eqtl"))[-1]; made.append(run)
        summ = json.loads((run / "summary.json").read_text())
        check(rc == 0 and (run / "report.md").is_file() and set(summ["methods"]) == {"good", "random", "distance", "binary"}, "run: report + summary for 4 methods")
        good = summ["methods"]["good"]["all_matches_at_threshold"]
        rnd = summ["methods"]["random"]["all_matches_at_threshold"]
        check(good["enrichment"] > 20 and good["p_adjust"] < 1e-6, f"run: good predictor aggregated enrichment {good['enrichment']:.1f} (p_adj {good['p_adjust']:.1e})")
        check(abs(good["recall_linking"] - 2 / 3) < 0.02 and abs(good["recall_total"] - 2 / 3) < 0.02,
              f"run: good predictor recall linking {good['recall_linking']:.3f} = 2/3 (the far variant per gene is unreachable)")
        check(rnd["enrichment"] < good["enrichment"] and rnd["recall_linking"] < good["recall_linking"],
              f"run: random scores below the good predictor (enrichment {rnd['enrichment']:.1f}, recall linking {rnd['recall_linking']:.3f})")
        # unmatched pair (bioA predictions vs TissueB eQTLs) should not be enriched like the matched one
        enr_d = pd.read_csv(run / "good" / "enrichment_by_distance.tsv", sep="\t")
        allbin = enr_d[enr_d["distance_max"] == 30000].set_index(["Biosample", "GTExTissue"])["enrichment"]
        check(allbin[("bioA", "Tissue_A")] > 5 * max(allbin[("bioA", "Tissue_B")], 1e-9) and allbin[("bioB", "Tissue_B")] > 5 * max(allbin[("bioB", "Tissue_A")], 1e-9),
              f"run: cell-type specificity — matched pairs enriched ({allbin[('bioA', 'Tissue_A')]:.1f}, {allbin[('bioB', 'Tissue_B')]:.1f}), swapped pairs not "
              f"({allbin[('bioA', 'Tissue_B')]:.2f}, {allbin[('bioB', 'Tissue_A')]:.2f})")
        # aggregate = sums over matched pairs
        enr_t = pd.read_csv(run / "good" / "enrichment_across_thresholds.tsv.gz", sep="\t")
        t0 = 0.6
        at = enr_t[np.isclose(enr_t["threshold"], t0)]
        m = at[(at["GTExTissue"] + "." + at["Biosample"]).isin({"Tissue_A.bioA", "Tissue_C.bioA", "Tissue_B.bioB"})]
        agg = at[at["GTExTissue"] == "all_matches"].iloc[0]
        exp_enr = (m["nVariantsOverlappingEnhancers"].sum() / m["nVariantsGTExTissue"].sum()) / (m["nCommonVariantsOverlappingEnhancers"].sum() / m["nCommonVariants"].sum())
        check(abs(agg["enrichment"] - exp_enr) < 1e-9 and agg["nVariantsGTExTissue"] == m["nVariantsGTExTissue"].sum(),
              "run: all_matches enrichment = ratio of summed counts over the 3 matched pairs")
        # threshold span contains the method threshold; binary has 2 thresholds
        check(summ["methods"]["binary"]["n_thresholds"] == 2 and summ["methods"]["good"]["n_thresholds"] >= 20, "run: binary span {0,1}; continuous span >= nThresholdSteps")
        # inverse predictor: negated threshold, and being close to a TSS scores well
        dist = summ["methods"]["distance"]
        check(dist["threshold"] == -60000 and 0.15 < dist["all_matches_at_threshold"]["recall_total"] < 0.5,
              f"run: inverse predictor negated (threshold {dist['threshold']}); recall total {dist['all_matches_at_threshold']['recall_total']:.3f} within 60 kb")
        rec_d = pd.read_csv(run / "good" / "recall_by_distance.tsv", sep="\t")
        r = rec_d[(rec_d["Biosample"] == "bioA") & (rec_d["GTExTissue"] == "Tissue_A")].set_index("distance_max")
        check(r.loc[30000, "recall.total"] < 1.0 and r.loc[100, "total.variants"] + r.loc[250, "total.variants"] + r.loc[1000, "total.variants"] + r.loc[10, "total.variants"] < r.loc[30000, "total.variants"],
              "run: by-distance recall table — bins partition the near variants; far variants only in 'all'")
        er = pd.read_csv(run / "enrichment_recall.AllMatches.tsv", sep="\t")
        g = er[er["method"] == "good"].sort_values("threshold")
        check(len(g) == summ["methods"]["good"]["n_thresholds"] and g["recall.linking"].is_monotonic_decreasing,
              "run: enrichment-recall curve has one row per threshold with recall non-increasing in the threshold")
        ear_files = sorted(run.glob("enrichment_at_recall0.6.*.tsv"))
        check(len(ear_files) >= 1, f"run: enrichment-at-recall tables written ({len(ear_files)}) with pairwise comparisons")
        comp = sorted(run.glob("pairwise_comparisons.recall0.6.*.tsv"))
        if comp:
            c = pd.read_csv(comp[0], sep="\t")
            check({"group1", "group2", "p", "p_adjust", "significant"} <= set(c.columns), "run: pairwise log-RR z-tests, Bonferroni")
        sizes = pd.read_csv(run / "good" / "enhancer_set_sizes.tsv", sep="\t")
        check((sizes["bp"] > 0).all() and (sizes["bp"] <= sizes["n_elements_at_threshold"] * 400).all(), "run: enhancer set sizes = merged bp at threshold")
        check((run / "heatmap_enrichment.good.tsv").is_file() and (run / "all_matched_enrichments.tsv").is_file(), "run: heatmap matrix and all-matched tables")
        if not args.no_plots:
            check(len(summ["figures"]) >= 5, f"run: {len(summ['figures'])} figures")
        rc = cmd_run(argparse.Namespace(label="st_eqtl_compat", methods_table=str(W["methods"]), predictions_table=str(W["predictions"]), methods=["good"],
                                        eqtl=None, bg_variants=None, variants_dir=str(vdir), partition=str(W["partition"]), tss=str(W["tss"]),
                                        threshold_pip=0.5, distances=[10, 100, 250, 1000], recalls=[0.3], threshold_pval=0.05,
                                        n_threshold_steps=20, upstream_compat=True, no_plots=True))
        run2 = sorted(OUT_ROOT.glob("*_st_eqtl_compat"))[-1]; made.append(run2)
        e2 = pd.read_csv(run2 / "good" / "enrichment_by_distance.tsv", sep="\t")
        e2 = e2[(e2["distance_max"] == 30000) & (e2["Biosample"] == "bioA") & (e2["GTExTissue"] == "Tissue_A")].iloc[0]
        e1 = enr_d[(enr_d["distance_max"] == 30000) & (enr_d["Biosample"] == "bioA") & (enr_d["GTExTissue"] == "Tissue_A")].iloc[0]
        check(rc == 0 and e2["nCommonVariantsOverlappingEnhancers"] >= e1["nCommonVariantsOverlappingEnhancers"] and e2["nVariantsOverlappingEnhancers"] >= e1["nVariantsOverlappingEnhancers"],
              "run --upstream-compat: unthresholded 'all'-bin and background counts reproduced")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent eqtl-enrich", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch the upstream resource files (and, with --synapse, the background SNPs + GTEx release)")
    p.add_argument("--synapse", action="store_true", help="also download syn52264319 and syn52264297 (needs SYNAPSE_AUTH_TOKEN)")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("prepare-gtex", help="raw GTEx fine-mapping table -> 7-column eQTL input")
    p.add_argument("--raw", required=True, help="GTEx_30tissues_release1.tsv.gz (19 columns, no header)")
    p.add_argument("--method", default="SUSIE")
    p.add_argument("--keep-all-cs", action="store_true", help="keep variants outside a credible set (cs_id == -1)")
    p.add_argument("--gene-table", help="CollapsedGeneBounds.hg38.tsv (default: resources)")
    p.add_argument("--expression", help="GTEx median-TPM .gct(.gz); rows with eGene TPM <= --tpm in the tissue are dropped")
    p.add_argument("--tpm", type=float, default=1.0)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_prepare_gtex)

    def add_variant_args(p):
        p.add_argument("--partition", help="PartitionCombined.bed (default: resources)")
        p.add_argument("--tss", help="TSS reference bed6 (default: resources)")
        p.add_argument("--threshold-pip", type=float, default=0.5)
        p.add_argument("--distances", type=int, nargs="+", default=[10, 100, 250, 1000], help="distance bin edges in kb")

    p = sub.add_parser("variants", help="filter eQTL variants and background SNPs once, for reuse by run --variants-dir")
    p.add_argument("--eqtl", required=True, help="eQTL variants TSV(.gz): chr,start,end,varID_hg38,gene_hgnc,tissue,pip")
    p.add_argument("--bg-variants", required=True, help="background SNPs bed: chr,start,end,rsid (no header)")
    add_variant_args(p)
    p.add_argument("--label", default="variants")
    p.set_defaults(func=cmd_variants)

    p = sub.add_parser("run", help="the benchmark")
    p.add_argument("--methods-table", required=True, help="method,boolean,inverse_predictor,pred_name_long,threshold,score_col,color[,geneUniverse]")
    p.add_argument("--predictions-table", required=True, help="biosample,<method columns with prediction paths>,GTExTissue (comma-separated matches)")
    p.add_argument("--methods", nargs="+", help="subset of methods to run")
    p.add_argument("--eqtl")
    p.add_argument("--bg-variants")
    p.add_argument("--variants-dir", help="output of `variants` (skips the variant filtering)")
    add_variant_args(p)
    p.add_argument("--recalls", type=float, nargs="*", default=[0.03, 0.05, 0.1], help="target recall (linking) values for enrichment-at-recall")
    p.add_argument("--threshold-pval", type=float, default=0.05)
    p.add_argument("--n-threshold-steps", type=int, default=50)
    p.add_argument("--upstream-compat", action="store_true", help="reproduce the upstream SE formula and by-distance counting quirks")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="eqtl_benchmark")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic genome, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
