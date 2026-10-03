#!/usr/bin/env python3
"""IGVF single-cell CRISPR screen (Perturb-seq) pipeline, end to end (port of pinellolab/CRISPR_Pipeline).

Port of https://github.com/pinellolab/CRISPR_Pipeline (MIT; Nextflow DSL2 +
Python + R), the IGVF consortium pipeline that turns Perturb-seq FASTQs and a
guide design table into ``inference_mudata.h5mu`` (gene / guide [/ hashing]
modalities, a ``guide_assignment`` layer, cis and trans per-guide and
per-element test results) plus QC tables, evaluation plots and an HTML
dashboard.  Every Nextflow module, subworkflow and ``bin/`` script (61 Python,
2 R) was read; the analysis steps are rewritten here in Python (numpy / scipy
/ pandas / anndata / h5py, sklearn for ROC/PR).  Relationship: port (source
consulted, logic re-derived, no code copied).  External binaries (kb /
kallisto / bustools, seqspec, GMM-demux, cleanser, SCEPTRE in R, PerTurbo in
torch) are OPTIONAL: when one is on PATH (or --upstream-bin points at the
upstream scripts) the exact upstream command runs; otherwise the command is
printed and a Python implementation of the method runs instead.

Definitions reproduced from upstream (bin/*.py, bin/*.R, nextflow.config)
  seqspec -> kb -x   leaf regions walked from each read's primer (pos strand
                     forward, neg strand backward), max_len per region, clipped
                     to the read; string "bc:umi:feature", each "file,start,stop";
                     feature types cdna/gdna/protein/tag/sgrna_target/crispr.
                     With a spacer tag the feature part becomes "file,0,0".
  seqSpecCheck       per read x orientation: score = 3*HitRatio + 2*PosPurity
                     + 1*FlankPurity + 1*(1-Gini) over the first max_reads reads;
                     12-bp upstream flank; winner = argmax.
  feature reference  guide_features.txt = [spacer, guide_id] (reverse
                     complemented / spacer-tag prefixed on request);
                     hashing_table.txt = [sequence, hash_id]; kite 1-mismatch
                     variants, collisions dropped.
  knee filter        barcode rank vs log1p(UMI); knee1 = max distance to the
                     chord; knee2 = same inside the concave window found from
                     the 2nd derivative of UnivariateSpline(s=n) of the curve
                     up to knee1 (10% trimmed ends); keep UMI >= count at knee.
  preprocessing      QC_barcode_filter knee2 (else min_genes 500), filter_genes
                     min_cells=10, pct_counts_mt < 15, mt = symbol prefix MT-.
  MuData             guide metadata merged on guide_id (one-to-one); type
                     filled from targeting; controls (targeting False or type in
                     safe-targeting / non-targeting / negative control) bucketed
                     as "non-targeting|k", k = rank(guide_id) // median guides
                     per targeting element + 1, chr "non-targeting", start=end=k;
                     intended_target_key "coords::name::chr::start::end" or
                     "name_only::name"; moi high/low or mean guides/cell > 1.5;
                     obs num_expressed_guides / total_guide_umis; gene coords
                     from the GTF; barcodes intersected across modalities.
  guide assignment   per batch; SCEPTRE mixture (Poisson GLM mixture, EM with
                     n_em_rep restarts, posterior >= 0.8); cleanser (posterior
                     >= threshold); thresholding (UMI >= t); then genes kept if
                     expressed in > int(0.05 * n_cells) cells; DUAL_GUIDE
                     collapse of 2 guides on one element into "g1|g2".
  pairs to test      guide x GTF gene on the guide chromosome with
                     |gene start - guide_start| <= 1,000,000 (strategy default /
                     by_distance), user pairs, or all_by_all; strategy default
                     subsets to tested genes, tested + control guides and
                     cells with any guide UMI (cis MuData).
  inference          SCEPTRE (union -> per element, singleton -> per guide;
                     complement control group, side both, skew-normal
                     resampling approximation, auto formula log(n_umis) +
                     log(n_nonzero) + batch) and PerTurbo (NB, batch,
                     library size total_gene_umis, covariate
                     log1p_total_guide_umis centred; loc/scale/ln2 -> log2_fc /
                     log2_fc_std).  Default method: cis = SCEPTRE + PerTurbo on
                     the pairs, trans = PerTurbo on all pairs; uns keys
                     cis_/trans_ per_guide_results / per_element_results.
  QC (additional_qc) gene/guide metric tables (median/mean/std/min/max/q25/q75
                     overall and per batch); intended target: gene_id ==
                     intended_target_name, strong knockdown log2FC <= log2(0.4),
                     p < 0.05, AUROC/AUPRC (score 1-p, AUPRC = trapezoid under
                     the PR curve) of direct-target tests vs NT guides tested
                     on the same genes, negatives downsampled (seed 42); trans:
                     BH over all tests, per-guide significant counts excluding
                     the intended target, validated-link AUROC.
  evaluate_controls  same positives/negatives on trans_per_guide_results with
                     p_value; volcano, PR/ROC, bar plot; IGV bedpe (element ->
                     gene) and bedgraph (element on its own gene).
  TF benchmark       ENCODE peaks (quantile >= 0) overlapping +/-w bp of the
                     mean TSS, w in 500/1000/2500/5000/10000; significant =
                     p_value <= 0.001 (already-corrected flag as upstream);
                     one-sided Fisher (greater) of significant vs peak target.

Deliberate deviations (each documented at the function; --upstream-compat
restores the upstream behaviour where it is a plain bug)
  * TF benchmark TSS: upstream compares a GTF strand '+' to 1, so every gene's
    TSS becomes its END; the port uses start for '+' genes (--upstream-compat
    reproduces the END-for-all behaviour).
  * mitochondrial prefix: upstream always uses 'MT-' (its mouse prefix variable
    is unused); the port matches 'MT-' case-insensitively for non-human
    references (--upstream-compat: 'MT-' case-sensitive).
  * evaluate_controls samples as many negatives as positives with replacement
    False and crashes when there are fewer; the port samples min(n_pos, n_neg).
  * Python engines are method ports, not the packages: the SCEPTRE slot is a
    NB score statistic calibrated by label permutations (complement control
    group) with a moment-fitted skew-normal; the PerTurbo slot is a NB
    maximum-likelihood fold change on the covariate-adjusted null with a Wald
    test (no guide-efficacy latent); GMM-demux, scrublet and CLEANSER have
    Python fallbacks.  Results carry the upstream column names and
    uns['inference_engine'] says which engine produced them.

Subcommands
  setup            fetch the 5 ENCODE ChIP-seq BED files and example seqspecs
                   from the pinned upstream commit into Data/CRISPRPipeline
  samplesheet      validate a pipeline samplesheet; covariate table, formula,
                   per-modality FASTQ groups; --from-portal maps a portal TSV
  portal-samplesheet  analysis-set accession -> per-sample TSV (IGVF portal)
  portal-download  download/verify (md5) FASTQ, seqspec, onlist, guide design
  interface-config config table -> pipeline_input.config
  seqspec          seqspec YAML -> kb technology string (parsed_seqSpec.txt)
  seqspec-check    guide/hash position scan of R1/R2 FASTQs
  feature-ref      guide / hashing feature tables, kite mismatch FASTA, kb ref
  map              kb count wrapper (rna / guide / hash) or --engine python
                   feature counting for guide / hash libraries
  concat           <batch>_ks_*_out dirs -> one AnnData with batch covariates
  preprocess       knee barcode filter, QC metrics, mito / gene filters
  create-mudata    RNA + guide (+ hashing) -> mudata.h5mu
  hashing          hashing filter + demultiplex (GMM-demux or fallback)
  doublets         scrublet (or fallback) doublet removal on the MuData
  assign-guides    SCEPTRE mixture / cleanser / threshold assignment
  nt-targets       control-guide intended targets (same / median strategy)
  pairs            pairs_to_test + prepare_inference (cis subset)
  chunk            gene-chunked MuData files + manifest (SCEPTRE / PerTurbo)
  inference        sceptre / perturbo / sceptre,perturbo / default
  merge-results    merge SCEPTRE + PerTurbo, cis + trans, chunk results,
                   export per_guide_output / per_element_output
  qc               gene, guide, intended-target and trans QC + HTML report
  evaluate         evaluate_controls, IGV tracks, volcano and network plots
  expression       per-cell expression of a gene: guide cells vs NT-only cells
  tf-benchmark     ENCODE ChIP-seq TF target enrichment of trans results
  dashboard        dashboard plots + self-contained HTML dashboard
  run              count matrices -> ... -> dashboard in one run directory
  nextflow         print (or run) the upstream Nextflow command
  selftest         synthetic screen with planted effects; every subcommand

Output: Docs/CRISPRPipeline/<timestamp>_<label>/.

Usage:
    igvfagent crispr-pipeline seqspec --yaml guide_seqspec.yml --modalities guide --whitelist 737K.txt
    igvfagent crispr-pipeline seqspec-check --read1 g_R1.fq.gz --read2 g_R2.fq.gz --metadata guide_metadata.tsv
    igvfagent crispr-pipeline map --modality guide --fastqs g_R1.fq.gz g_R2.fq.gz --features guide_metadata.tsv --technology 0,0,16:0,16,28:1,0,20 --barcodes 737K.txt --batch B1 --engine python
    igvfagent crispr-pipeline run --rna B1_ks_transcripts_out B2_ks_transcripts_out --guide B1_ks_guide_out B2_ks_guide_out --guide-metadata guide_metadata.tsv --gtf gencode.gtf.gz --label screen1
    igvfagent crispr-pipeline inference --mudata mudata_inference_input.h5mu --method default --label cis_trans
    igvfagent crispr-pipeline qc --mudata inference_mudata.h5mu --label qc
    igvfagent crispr-pipeline evaluate --mudata inference_mudata.h5mu --gtf gencode.gtf.gz
    igvfagent crispr-pipeline tf-benchmark --mudata inference_mudata.h5mu --gtf gencode.gtf.gz
    igvfagent crispr-pipeline selftest --no-plots
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import io
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
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRPipeline"
DATA_ROOT = ROOT / "Data" / "CRISPRPipeline"

UPSTREAM_REPO = "pinellolab/CRISPR_Pipeline"
UPSTREAM_COMMIT = "5c2945402928fe10efdfda1bc497bc5de8ed2cc5"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
PORTAL = "https://api.data.igvf.org"

# nextflow.config defaults
DEFAULTS = {
    "QC_min_genes_per_cell": 500, "QC_min_cells_per_gene": 0.05, "QC_pct_mito": 15.0, "QC_barcode_filter": "knee2",
    "Multiplicity_of_infection": "high", "GUIDE_ASSIGNMENT_method": "sceptre", "GUIDE_ASSIGNMENT_capture_method": "CROP-seq",
    "GUIDE_ASSIGNMENT_cleanser_probability_threshold": 1.0, "GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold": 0.8,
    "GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep": 5, "INFERENCE_method": "default", "INFERENCE_target_guide_pairing_strategy": "default",
    "INFERENCE_max_target_distance_bp": 1_000_000, "INFERENCE_SCEPTRE_side": "both",
    "INFERENCE_SCEPTRE_resampling_approximation": "skew_normal", "INFERENCE_SCEPTRE_control_group": "complement",
    "INFERENCE_SCEPTRE_GENE_CHUNK_SIZE": 1000, "INFERENCE_SCEPTRE_MAX_MATRIX_ENTRIES": 2147483647,
    "INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK": 8000, "NETWORK_central_nodes_num": 1,
}
CONTROL_TYPES = {"safe-targeting", "non-targeting", "negative control"}
NT_LABEL = "non-targeting"
RESULT_KEYS = ["trans_per_guide_results", "per_guide_results", "cis_per_guide_results", "trans_test_results", "test_results"]
LOG2FC_CANDIDATES = ["log2_fc", "perturbo_log2_fc", "sceptre_log2_fc"]
PVALUE_CANDIDATES = ["p_value", "perturbo_p_value", "sceptre_p_value", "q_value"]
FEATURE_REGION_TYPES = {"cdna", "gdna", "protein", "tag", "sgrna_target", "crispr"}
PEAK_META_DICT = {  # tf_benchmark.py
    "ENSG00000166949": "ENCFF580MYL.bed.gz",  # SMAD3
    "ENSG00000081059": "ENCFF032CCV.bed.gz",  # TCF7
    "ENSG00000141646": "ENCFF126SGU.bed.gz",  # SMAD4
    "ENSG00000149948": "ENCFF614IHQ.bed.gz",  # HMGA2
    "ENSG00000072364": "ENCFF835WQS.bed.gz",  # AFF4
}
ENSEMBL_GENE_NAME_MAP = {"ENSG00000166949": "SMAD3", "ENSG00000081059": "TCF7", "ENSG00000141646": "SMAD4",
                         "ENSG00000149948": "HMGA2", "ENSG00000072364": "AFF4"}
PROMOTER_WINDOWS = [500, 1000, 2500, 5000, 10000]
SETUP_FILES = [f"encode_bed_files/{v}" for v in PEAK_META_DICT.values()] + [
    "example_data/yaml_files/rna_seqspec.yml", "example_data/yaml_files/guide_seqspec.yml",
    "example_data/yaml_files/hash_seqspec.yml", "example_data/guide_metadata.tsv", "example_data/hash_metadata.tsv"]

BLUE, ORANGE, RED, GREEN, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#d62728", "#1baf7a", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
GREY = "#8a8984"
log = logging.getLogger("crispr_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    k = 1
    while d.exists():
        k += 1
        d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_{k}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("crispr-pipeline needs pandas + numpy + scipy + anndata + h5py: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


def _sp():
    import scipy.sparse as sp  # type: ignore
    return sp


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
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), bbox_inches="tight", dpi=110)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(_fmt(c) if isinstance(c, float) else str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        if isinstance(x, (int,)) and not isinstance(x, bool):
            return str(x)
        v = float(x)
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
            return f"{v:.2e}"
        return f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _truthy(v) -> bool:
    np = _np()
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    try:
        import pandas as pd  # type: ignore
        if pd.isna(v):
            return False
    except (TypeError, ValueError):
        pass
    return str(v).strip().lower() in {"true", "1", "t", "yes", "y"}


def _jsonable(x):
    np = _np()
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if math.isfinite(v) else None
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, Path):
        return str(x)
    return x


def write_tsv(df, path: Path, sep: str = "\t", quiet: bool = False, **kw) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=kw.pop("index", False), compression="gzip" if str(path).endswith(".gz") else None, **kw)
    if not quiet:
        print(("CSV: " if sep == "," else "TSV: ") + str(path))
    return path


def write_json(obj, path: Path, quiet: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(obj), indent=2, default=str))
    if not quiet:
        print(f"JSON: {path}")
    return path


def finish(out: Path, title: str, summary: dict, sections: "list[str]") -> None:
    rep = out / "report.md"
    rep.write_text(f"# {title}\n\n" + "\n\n".join(s for s in sections if s) + "\n")
    print(f"Report: {rep}")
    write_json(summary, out / "summary.json")


def _open_text(path: Path):
    p = str(path)
    with open(p, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(p, "rt") if magic == b"\x1f\x8b" else open(p, "rt")


def which(binary: str) -> Optional[str]:
    return shutil.which(binary)


def external(cmd: "list[str]", binary: str, cwd: Optional[Path] = None, run: bool = True) -> "tuple[bool, str]":
    """Run an upstream command if its binary exists; otherwise print it and return (False, cmd)."""
    text = " ".join(str(c) for c in cmd)
    if run and which(binary):
        print(f"Running: {text}")
        log.info("run %s", text)
        res = subprocess.run([str(c) for c in cmd], cwd=str(cwd) if cwd else None, capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"{binary} failed ({res.returncode}): {res.stderr[-2000:]}")
        return True, text
    print(f"[{binary} not on PATH] command that would run:\n  {text}")
    return False, text


# ---------------------------------------------------------------------------
# AnnData / MuData I/O (mudata optional: an .h5mu is HDF5 with AnnData groups under mod/)
# ---------------------------------------------------------------------------

def _elem_io():
    try:
        from anndata.io import read_elem, write_elem  # type: ignore
    except ImportError:
        from anndata.experimental import read_elem, write_elem  # type: ignore
    return read_elem, write_elem


class MData:
    """The modalities of a MuData file (dict of AnnData) plus the global obs and uns."""

    def __init__(self, mod: dict, obs=None, uns: Optional[dict] = None) -> None:
        self.mod = dict(mod)
        first = next(iter(self.mod.values()))
        pd = _pd()
        self.obs = obs if obs is not None else pd.DataFrame(index=first.obs_names.copy())
        self.uns = dict(uns or {})

    @property
    def n_obs(self) -> int:
        return int(next(iter(self.mod.values())).n_obs)

    def __getitem__(self, key):
        return self.mod[key]

    def subset_cells(self, names) -> "MData":
        names = list(names)
        mod = {k: a[names, :].copy() for k, a in self.mod.items()}
        obs = self.obs.loc[names].copy() if set(names) <= set(self.obs.index) else None
        return MData(mod, obs, dict(self.uns))


def _sanitize_df(df):
    """Make a DataFrame writable by anndata: nullable ints -> float, mixed objects -> str."""
    pd, np = _pd(), _np()
    out = df.copy()
    for c in out.columns:
        s = out[c]
        if isinstance(s.dtype, pd.CategoricalDtype):
            cats = s.cat.categories
            if cats.dtype == object and any(not isinstance(v, str) for v in cats):
                out[c] = s.astype(str).astype("category")
            continue
        if str(s.dtype) in ("Int64", "Int32", "Int16", "Int8", "boolean", "Float64"):
            out[c] = s.astype("float64") if str(s.dtype) != "boolean" else s.astype(object).map(lambda v: bool(v) if pd.notna(v) else False)
        elif s.dtype == object:
            vals = s.dropna()
            if len(vals) and all(isinstance(v, (bool, np.bool_)) for v in vals) and not s.isna().any():
                out[c] = s.astype(bool)
            elif len(vals) and all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_)) for v in vals):
                out[c] = pd.to_numeric(s, errors="coerce").astype(float)
            else:
                out[c] = s.map(lambda v: "" if (v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA) else str(v))
    out.columns = [str(c) for c in out.columns]
    if out.index.dtype != object:
        out.index = out.index.astype(str)
    return out


def _sanitize_adata(a):
    b = a.copy()
    b.obs = _sanitize_df(b.obs)
    b.var = _sanitize_df(b.var)
    b.uns = _sanitize_uns(dict(b.uns))
    return b


def _sanitize_uns(uns: dict) -> dict:
    pd, np = _pd(), _np()
    out = {}
    for k, v in uns.items():
        if isinstance(v, pd.DataFrame):
            out[k] = _sanitize_df(v.reset_index(drop=True))
        elif isinstance(v, dict):
            try:
                out[k] = _sanitize_df(pd.DataFrame(v))
            except ValueError:
                out[k] = _sanitize_uns(v)
        elif isinstance(v, (list, tuple)):
            out[k] = np.array([str(x) for x in v], dtype=object) if any(isinstance(x, str) for x in v) else np.asarray(v)
        elif isinstance(v, str):
            out[k] = np.array([v], dtype=object)
        else:
            out[k] = v
    return out


def write_h5ad(a, path: Path, quiet: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    _sanitize_adata(a).write_h5ad(str(path))
    if not quiet:
        print(f"Wrote: {path}")
    return path


def read_h5ad(path: Path):
    import anndata as ad  # type: ignore
    return ad.read_h5ad(str(path))


def read_h5mu(path: Path) -> MData:
    try:
        import mudata as md  # type: ignore
        m = md.read_h5mu(str(path))
        return MData({k: m.mod[k] for k in m.mod}, m.obs.copy(), dict(m.uns))
    except ImportError:
        pass
    import h5py  # type: ignore
    read_elem, _ = _elem_io()
    with h5py.File(str(path), "r") as h:
        order = [k.decode() if isinstance(k, bytes) else str(k) for k in h["mod"].attrs.get("mod-order", list(h["mod"].keys()))]
        mods = {k: read_elem(h[f"mod/{k}"]) for k in order if k in h["mod"]}
        obs = read_elem(h["obs"]) if "obs" in h else None
        uns = read_elem(h["uns"]) if "uns" in h else {}
    return MData(mods, obs, dict(uns))


def write_h5mu(m: MData, path: Path, quiet: bool = False) -> Path:
    pd, np = _pd(), _np()
    path.parent.mkdir(parents=True, exist_ok=True)
    mods = {k: _sanitize_adata(a) for k, a in m.mod.items()}
    obs = _sanitize_df(m.obs) if m.obs is not None else pd.DataFrame(index=next(iter(mods.values())).obs_names)
    uns = _sanitize_uns(m.uns)
    try:
        import mudata as md  # type: ignore
        mm = md.MuData(mods)
        common = [c for c in obs.columns if c not in mm.obs.columns]
        for c in common:
            mm.obs[c] = obs.reindex(mm.obs_names)[c].values
        mm.uns = uns
        mm.write(str(path))
    except ImportError:
        import h5py  # type: ignore
        _, write_elem = _elem_io()
        names = list(mods)
        cells = mods[names[0]].obs_names
        with h5py.File(str(path), "w") as h:
            h.attrs["encoding-type"] = "MuData"
            h.attrs["encoding-version"] = "0.1.0"
            h.attrs["encoder"] = "igvfagent"
            g = h.create_group("mod")
            g.attrs["mod-order"] = np.array(names, dtype=object).astype("S")
            for k in names:
                write_elem(g, k, mods[k])
            write_elem(h, "obs", obs.reindex(cells))
            var_names = []
            for k in names:
                var_names.extend(list(mods[k].var_names))
            write_elem(h, "var", pd.DataFrame(index=pd.Index(var_names).astype(str)))
            om, vm = h.create_group("obsmap"), h.create_group("varmap")
            off = 0
            for k in names:
                pos = {c: i for i, c in enumerate(mods[k].obs_names)}
                om.create_dataset(k, data=np.array([pos.get(c, -1) + 1 for c in cells], dtype=np.uint32))
                vmap = np.zeros(len(var_names), dtype=np.uint32)
                n = mods[k].n_vars
                vmap[off:off + n] = np.arange(1, n + 1)
                off += n
                vm.create_dataset(k, data=vmap)
            write_elem(h, "uns", uns)
    if not quiet:
        print(f"Wrote: {path}")
    return path


def uns_df(m: MData, key: str):
    pd = _pd()
    v = m.uns.get(key)
    if v is None:
        return None
    if isinstance(v, pd.DataFrame):
        return v.copy()
    if isinstance(v, dict):
        return pd.DataFrame({k: (list(x) if not hasattr(x, "shape") else x) for k, x in v.items()})
    return pd.DataFrame(v)


def uns_str(m: MData, key: str, default: str = "") -> str:
    v = m.uns.get(key)
    if v is None:
        return default
    if hasattr(v, "__len__") and not isinstance(v, str):
        return str(list(v)[0]) if len(v) else default
    return str(v)


def dense_col(X, j: int):
    np, sp = _np(), _sp()
    col = X[:, j]
    return np.asarray(col.toarray()).ravel() if sp.issparse(col) else np.asarray(col).ravel()


def row_sums(X):
    np = _np()
    return np.asarray(X.sum(axis=1)).ravel()


def col_sums(X):
    np = _np()
    return np.asarray(X.sum(axis=0)).ravel()


def nnz_rows(X):
    np, sp = _np(), _sp()
    return np.asarray((X > 0).sum(axis=1)).ravel() if sp.issparse(X) else (np.asarray(X) > 0).sum(axis=1)


def as_csr(X):
    sp, np = _sp(), _np()
    return X.tocsr().astype(np.float64) if sp.issparse(X) else sp.csr_matrix(np.asarray(X, dtype=np.float64))


# ---------------------------------------------------------------------------
# GTF (gtfparse is not required)
# ---------------------------------------------------------------------------

def read_gtf(path: Path, features: "Optional[set]" = None):
    """GTF -> DataFrame (seqname, feature, start, end, strand, gene_id, gene_name, gene_type).

    Upstream reads the whole GTF with gtfparse and takes the first row per gene_id / gene_name; in a
    GENCODE-ordered GTF that is the 'gene' line, which is what is kept here by default (features={'gene'}).
    """
    pd = _pd()
    features = {"gene"} if features is None else features
    rows = []
    with _open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or (features and f[2] not in features):
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6], attrs.get("gene_id", ""), attrs.get("gene_name", attrs.get("gene_id", "")),
                         attrs.get("gene_type", attrs.get("gene_biotype", ""))))
    df = pd.DataFrame(rows, columns=["seqname", "feature", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])
    df["gene_id2"] = df["gene_id"].str.split(".").str[0]
    return df


# ---------------------------------------------------------------------------
# Guide metadata canonicalisation (intended_target_key_utils.py)
# ---------------------------------------------------------------------------

def _bool_series(values, index):
    pd = _pd()
    if values is None:
        return pd.Series(False, index=index)
    return pd.Series([_truthy(v) for v in values], index=index)


def _num_or_na(v):
    pd = _pd()
    try:
        if pd.isna(v):
            return pd.NA
    except (TypeError, ValueError):
        pass
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return pd.NA


def _str_or_na(v):
    pd = _pd()
    try:
        if pd.isna(v):
            return pd.NA
    except (TypeError, ValueError):
        pass
    return str(v)


def _norm_text(values, index):
    pd = _pd()
    if values is None:
        return pd.Series("", index=index)
    s = pd.Series(list(values), index=index).astype(object)
    return s.where(s.notna(), "").astype(str).str.strip().str.lower()


def build_target_key(name, chrom, start, end) -> str:
    pd = _pd()
    has = not (pd.isna(chrom) or pd.isna(start) or pd.isna(end))
    return f"coords::{name}::{chrom}::{int(start)}::{int(end)}" if has else f"name_only::{name}"


def control_mask(df):
    targeting = _bool_series(df["targeting"] if "targeting" in df.columns else None, df.index)
    gtype = _norm_text(df["type"] if "type" in df.columns else None, df.index)
    names = _norm_text(df["intended_target_name"] if "intended_target_name" in df.columns else None, df.index)
    return (~targeting) | gtype.isin(CONTROL_TYPES) | names.eq(NT_LABEL) | names.str.startswith(NT_LABEL + "|")


def annotate_intended_target_groups(guide_var, non_targeting_label: str = NT_LABEL):
    """Coordinate-aware target keys; controls bucketed as non-targeting|k (k-th block of median size)."""
    pd, np = _pd(), _np()
    df = guide_var.copy()
    if "guide_id" not in df.columns:
        df["guide_id"] = df.index.astype(str)
    if "intended_target_name" not in df.columns:
        targ = _bool_series(df["targeting"] if "targeting" in df.columns else None, df.index)
        df["intended_target_name"] = np.where(targ, df["guide_id"].astype(str), non_targeting_label)
    for c in ("intended_target_chr", "intended_target_start", "intended_target_end"):
        if c not in df.columns:
            df[c] = pd.NA
    df["guide_id"] = df["guide_id"].astype(str)
    df["intended_target_name"] = df["intended_target_name"].astype(object).map(_str_or_na)
    df["intended_target_chr"] = df["intended_target_chr"].astype(object).map(_str_or_na)
    df["intended_target_start"] = df["intended_target_start"].astype(object).map(_num_or_na).astype("Int64")
    df["intended_target_end"] = df["intended_target_end"].astype(object).map(_num_or_na).astype("Int64")
    cm = control_mask(df)
    n_ctrl = int(cm.sum())
    if n_ctrl:
        tdf = df.loc[~cm]
        if len(tdf):
            keys = [build_target_key(n, c, s, e) for n, c, s, e in zip(tdf["intended_target_name"], tdf["intended_target_chr"],
                                                                        tdf["intended_target_start"], tdf["intended_target_end"])]
            median = max(1, int(np.median(pd.Series(keys).value_counts().values)))
        else:
            median = 1
        idx = df.loc[cm, "guide_id"].astype(str).sort_values(kind="stable").index
        groups = np.arange(n_ctrl) // median + 1
        df.loc[idx, "intended_target_name"] = [f"{non_targeting_label}|{g}" for g in groups]
        df.loc[idx, "intended_target_chr"] = non_targeting_label
        df.loc[idx, "intended_target_start"] = pd.Series(groups, index=idx, dtype="Int64")
        df.loc[idx, "intended_target_end"] = pd.Series(groups, index=idx, dtype="Int64")
    if (df["intended_target_name"] == non_targeting_label).any():
        raise ValueError("found exact 'non-targeting' after bucketing")
    df["intended_target_key"] = [build_target_key(n, c, s, e) for n, c, s, e in zip(
        df["intended_target_name"], df["intended_target_chr"], df["intended_target_start"], df["intended_target_end"])]
    return df


def enrich_pairs_with_target_metadata(pairs, guide_var):
    pd = _pd()
    lookup = guide_var.reset_index(drop=True)[["guide_id", "intended_target_name", "intended_target_chr", "intended_target_start",
                                               "intended_target_end", "intended_target_key"]].drop_duplicates("guide_id")
    out = pairs.copy()
    for c in ("intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end", "intended_target_key"):
        if c not in out.columns:
            out[c] = pd.NA
    out["guide_id"] = out["guide_id"].astype(str)
    lookup = lookup.assign(guide_id=lookup["guide_id"].astype(str))
    mg = out.merge(lookup, on="guide_id", how="left", suffixes=("", "_from_guide"))
    for c in ("intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"):
        mg[c] = mg[c].where(mg[c].notna(), mg[f"{c}_from_guide"])
    mg["intended_target_name"] = mg["intended_target_name"].astype(object).map(_str_or_na)
    mg["intended_target_chr"] = mg["intended_target_chr"].astype(object).map(_str_or_na)
    mg["intended_target_start"] = mg["intended_target_start"].astype(object).map(_num_or_na).astype("Int64")
    mg["intended_target_end"] = mg["intended_target_end"].astype(object).map(_num_or_na).astype("Int64")
    mg["intended_target_key"] = [build_target_key(n, c, s, e) for n, c, s, e in zip(
        mg["intended_target_name"], mg["intended_target_chr"], mg["intended_target_start"], mg["intended_target_end"])]
    return mg.drop(columns=[c for c in mg.columns if c.endswith("_from_guide")])


def target_lookup(guide_var):
    cols = ["intended_target_key", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]
    return guide_var[cols].drop_duplicates("intended_target_key").reset_index(drop=True)


def guide_targeting(guide_var, non_targeting_name: str = NT_LABEL):
    """intended_target.py::_coerce_targeting_column -- bool, falling back to intended_target_name != NT."""
    pd = _pd()
    if "targeting" in guide_var.columns:
        def tb(v):
            if isinstance(v, bool):
                return v
            t = str(v).strip().upper()
            return True if t == "TRUE" else (False if t == "FALSE" else None)
        t = guide_var["targeting"].map(tb)
        if "intended_target_name" in guide_var.columns:
            fb = guide_var["intended_target_name"].astype(str) != non_targeting_name
            t = t.where(t.notna(), fb)
        return t.astype(bool)
    return (guide_var["intended_target_name"].astype(str) != non_targeting_name)


def guide_meta(guide_var, non_targeting_name: str = NT_LABEL):
    gm = guide_var.reset_index(drop=True) if "guide_id" in guide_var.columns else guide_var.reset_index().rename(columns={"index": "guide_id"})
    gm = gm.copy()
    gm["targeting"] = guide_targeting(gm, non_targeting_name).values
    cols = ["guide_id", "intended_target_name", "targeting"] + (["gene_name"] if "gene_name" in gm.columns else [])
    return gm[cols].drop_duplicates("guide_id")


def resolve_results_key(m: MData, requested: str = "auto") -> str:
    if requested != "auto" and requested in m.uns:
        return requested
    for k in RESULT_KEYS:
        if k in m.uns:
            return k
    raise SystemExit(f"no inference results in mdata.uns (have {sorted(m.uns)})")


def resolve_metric_columns(res, log2fc_col: str = "auto", pvalue_col: str = "auto") -> "tuple[str, str]":
    def pick(req, cands):
        if req != "auto" and req in res.columns:
            return req
        for c in cands:
            if c in res.columns and res[c].notna().any():
                return c
        raise SystemExit(f"cannot resolve column among {cands}; have {list(res.columns)}")
    return pick(log2fc_col, LOG2FC_CANDIDATES), pick(pvalue_col, PVALUE_CANDIDATES)


def bh_adjust(p):
    np = _np()
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def roc_pr(labels, scores) -> "tuple[float, float, dict]":
    """sklearn roc_curve / precision_recall_curve; AUPRC = auc(recall, precision) as upstream."""
    np = _np()
    from sklearn.metrics import auc, precision_recall_curve, roc_curve  # type: ignore
    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return float("nan"), float("nan"), {}
    fpr, tpr, _ = roc_curve(labels, scores)
    pre, rec, _ = precision_recall_curve(labels, scores)
    return float(auc(fpr, tpr)), float(auc(rec, pre)), {"fpr": fpr, "tpr": tpr, "precision": pre, "recall": rec}


# ---------------------------------------------------------------------------
# setup: pinned upstream resources (small: 5 ENCODE BED files ~830 KB + example seqspecs)
# ---------------------------------------------------------------------------

def cmd_setup(args) -> int:
    dest = Path(args.dest) if args.dest else DATA_ROOT
    got = []
    for rel in SETUP_FILES:
        out = dest / rel
        if out.is_file() and out.stat().st_size > 0 and not args.force:
            print(f"  have  {out}")
            got.append(str(out))
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        url = RAW + rel
        try:
            with urllib.request.urlopen(url, timeout=60) as r, open(out, "wb") as fh:
                fh.write(r.read())
            print(f"Wrote: {out}")
            got.append(str(out))
        except Exception as e:  # pragma: no cover - network
            print(f"  FAIL  {url}: {e}")
    print(f"setup: {len(got)}/{len(SETUP_FILES)} files under {dest} (upstream {UPSTREAM_REPO}@{UPSTREAM_COMMIT[:10]})")
    return 0 if len(got) == len(SETUP_FILES) else 1


# ---------------------------------------------------------------------------
# samplesheet / covariates (PIPELINE_INITIALISATION, prepare_covariate, parse_covariate, prepare_formula,
# update_samplesheet, process_reads, process_batches)
# ---------------------------------------------------------------------------

SAMPLESHEET_COLS = ["R1_path", "R2_path", "file_modality", "measurement_sets", "sequencing_run", "lane", "seqspec",
                    "barcode_onlist", "guide_design", "barcode_hashtag_map"]
PORTAL_MODALITY = {"scRNA sequencing": "scRNA", "gRNA sequencing": "gRNA", "cell hashing barcode sequencing": "hash"}


def portal_to_samplesheet(df):
    """update_samplesheet.py: portal per-sample TSV -> pipeline CSV columns, modality names mapped."""
    cols = [c for c in ["R1_path", "R2_path", "file_modality", "measurement_sets", "sequencing_run", "lane", "flowcell_id", "seqspec",
                        "barcode_onlist", "guide_design", "barcode_hashtag_map"] if c in df.columns]
    out = df[cols].copy()
    out["file_modality"] = out["file_modality"].replace(PORTAL_MODALITY)
    for c in ("local_R1_path", "local_R2_path", "local_seqspec", "local_barcode_onlist", "local_guide_design", "local_barcode_hashtag_map"):
        if c in df.columns:
            out[c[len("local_"):]] = df[c].where(df[c].notna(), out[c[len("local_"):]])
    return out


def parse_samplesheet(df) -> dict:
    """Group FASTQs as PIPELINE_INITIALISATION does: id = sample_<file_modality>_<measurement_sets>."""
    pd = _pd()
    missing = [c for c in ("R1_path", "file_modality", "measurement_sets") if c not in df.columns]
    if missing:
        raise SystemExit(f"samplesheet is missing columns: {missing}")
    groups: "dict[str, dict]" = {}
    for _, r in df.iterrows():
        mid = f"sample_{r['file_modality']}_{r['measurement_sets']}"
        single = not (isinstance(r.get("R2_path"), str) and r.get("R2_path"))
        g = groups.setdefault(mid, {"id": mid, "modality": str(r["file_modality"]).lower(), "measurement_sets": str(r["measurement_sets"]),
                                    "single_end": single, "fastqs": [], "meta": {c: r.get(c) for c in SAMPLESHEET_COLS if c in df.columns}})
        if g["single_end"] != single:
            raise SystemExit(f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: {mid}")
        g["fastqs"].extend([r["R1_path"]] + ([] if single else [r["R2_path"]]))
    mods = {"scrna": [], "grna": [], "hash": []}
    for g in groups.values():
        mods.setdefault(g["modality"], []).append(g)

    def first(mod, col):
        vals = [g["meta"].get(col) for g in mods.get(mod, []) if isinstance(g["meta"].get(col), str) and g["meta"].get(col)]
        return vals[0] if vals else None
    batches = sorted({g["measurement_sets"] for g in groups.values()})
    cov = pd.DataFrame({"batch": batches})
    formula = " + ".join([c for c in cov.columns if cov[c].nunique() > 1])
    return {"groups": groups, "modalities": {k: [g["id"] for g in v] for k, v in mods.items()},
            "rna_seqspec": first("scrna", "seqspec"), "guide_seqspec": first("grna", "seqspec"), "hash_seqspec": first("hash", "seqspec"),
            "barcode_onlist": first("scrna", "barcode_onlist"), "guide_design": first("grna", "guide_design"),
            "barcode_hashtag_map": first("hash", "barcode_hashtag_map"), "covariates": cov, "formula": formula,
            "hashing": bool(mods.get("hash"))}


def cmd_samplesheet(args) -> int:
    pd = _pd()
    p = Path(args.samplesheet)
    sep = "\t" if p.suffix in (".tsv", ".txt") or (p.suffix == ".gz" and ".tsv" in p.name) else ","
    df = pd.read_csv(p, sep=sep, dtype=str)
    out = run_dir(args.label)
    if args.from_portal or "file_set" in df.columns and "R1_md5sum" in df.columns and not df["file_modality"].isin(["scRNA", "gRNA", "hash"]).all():
        df = portal_to_samplesheet(df)
        write_tsv(df, out / "updated_samplesheet.csv", sep=",")
    info = parse_samplesheet(df)
    write_tsv(info["covariates"], out / "parse_covariate.csv", sep=",")
    (out / "cov_string.txt").write_text(info["formula"])
    rows = [[g["id"], g["modality"], g["measurement_sets"], "single" if g["single_end"] else "paired", len(g["fastqs"])] for g in info["groups"].values()]
    batches_txt = "\n".join(f"{i} {' '.join(g['fastqs'])}" for i, g in enumerate(info["groups"].values()))
    (out / "fastq_batches.txt").write_text(batches_txt + "\n")
    print(f"Wrote: {out / 'fastq_batches.txt'}")
    summ = {k: v for k, v in info.items() if k not in ("groups", "covariates")}
    summ["n_samples"] = len(info["groups"])
    finish(out, "CRISPR pipeline samplesheet", summ, [
        md_table(["id", "modality", "measurement_sets", "ends", "n_fastqs"], rows),
        f"Covariate formula: `{info['formula'] or '(none: one batch)'}`; hashing: {info['hashing']}",
        f"RNA seqspec: `{info['rna_seqspec']}`; guide seqspec: `{info['guide_seqspec']}`; hash seqspec: `{info['hash_seqspec']}`",
        f"Barcode onlist: `{info['barcode_onlist']}`; guide design: `{info['guide_design']}`; hashtag map: `{info['barcode_hashtag_map']}`"])
    return 0


def interface_config(config_table, mapping_tabs=("scRNA", "Guides", "Hash")) -> str:
    """parse_interface_configuration.py: config table (tab_name, batch_name, read1, read2, variable, variable_value) -> config text."""
    pd = _pd()
    seqs, tests, hashing = {}, [], "false"
    keymap = {"scRNA": "fastq_files_rna", "Guides": "fastq_files_guide", "Hash": "fastq_files_hashing"}
    for tab in mapping_tabs:
        sub = config_table[config_table["tab_name"] == tab]
        if sub.empty:
            continue
        grouped = [" ".join(f"{a} {b}" for a, b in zip(x["read1"], x["read2"])) for _, x in sub.groupby("batch_name")]
        seqs[tab] = f"{keymap.get(tab, 'fastq_files')} = [\n    " + ",\n    ".join(f'"{s}"' for s in grouped) + "\n    ]"
        if tab in ("Guides", "Hash"):
            d = sub.dropna(subset=["read1", "read2"])
            r1k, r2k = ("test_guide_fastq_r1", "test_guide_fastq_r2") if tab == "Guides" else ("test_hashing_fastq_r1", "test_hashing_fastq_r2")
            tests.append(f"{r1k} = {[str(x).strip() for x in d['read1']]}\n    {r2k} = {[str(x).strip() for x in d['read2']]}")
        if tab == "Hash":
            hashing = "true"
    nonseq = config_table[config_table["variable"] != "sequence"][["variable", "variable_value"]]
    lines, has_dist = [], False
    for _, r in nonseq.iterrows():
        name, val = r["variable"], r["variable_value"]
        if pd.isna(name):
            continue
        has_dist |= name == "distance_from_center"
        sval = str(val)
        try:
            v = float(sval) if "." in sval else int(sval)
        except ValueError:
            v = f"'{sval}'"
        lines.append(f"{name} = {v}")
    if not has_dist:
        lines.append("distance_from_center = 1000000")
    bn = config_table[config_table["batch_name"].notna()]["batch_name"].map(lambda x: str(x).split(", "))
    maxs = bn.map(len).max() if len(bn) else 1
    covs = pd.DataFrame(bn.tolist(), columns=["batch"] + [f"cov{i}" for i in range(1, maxs)]).drop_duplicates() if len(bn) else pd.DataFrame(columns=["batch"])
    batch = bn.map(lambda x: x[0]) if len(bn) else pd.Series(dtype=str)
    cov_list = f"params.covariate_list = [\n    batch: {batch.unique().tolist()},\n"
    for c, v in covs.drop(columns=["batch"]).to_dict(orient="list").items():
        cov_list += f"    {c}: {v},\n"
    cov_list = cov_list.rstrip(",\n") + "\n]"
    seq_txt = "".join(f"{s}\n    " for s in seqs.values())
    return (f"params {{\n    DATASET_HASHING = '{hashing}'\n\n    " + "\n    ".join(lines) + f"\n\n    {seq_txt}\n\n    "
            + "\n\n    ".join(tests) + f"\n    batch={batch.unique().tolist()}\n}}\n\n{cov_list}")


def cmd_interface_config(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    txt = interface_config(pd.read_csv(args.config_table), tuple(args.mapping_tabs))
    p = out / "pipeline_input.config"
    p.write_text(txt)
    print(f"Wrote: {p}")
    finish(out, "Pipeline interface configuration", {"config": str(p)}, ["```groovy\n" + txt[:3000] + "\n```"])
    return 0


# ---------------------------------------------------------------------------
# IGVF portal: per-sample TSV + downloads (download_pipeline/bin/generate_per_sample.py, download_igvf.py)
# ---------------------------------------------------------------------------

def _portal_fetch(path: str) -> dict:  # pragma: no cover - network
    url = path if path.startswith("http") else PORTAL + path
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    key, sec = os.environ.get("IGVF_API_KEY"), os.environ.get("IGVF_SECRET_KEY")
    if key and sec:
        req.add_header("Authorization", "Basic " + base64.b64encode(f"{key}:{sec}".encode()).decode())
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def generate_per_sample(accession: str, fetch: Callable[[str], dict], hash_seqspec=None, rna_seqspec=None, sgrna_seqspec=None) -> "list[dict]":
    """Analysis set -> one row per R1/R2 pair, with the portal's validity rules (see generate_per_sample.py)."""
    raw_fetch = fetch
    fetch = lambda path: raw_fetch(re.sub(r"/+@@", "/@@", path))  # noqa: E731  (portal @ids end in '/')
    aset = fetch(f"/analysis-sets/{accession}/@@object?format=json")
    cls = aset.get("construct_library_sets", [])
    if len(cls) != 1:
        raise SystemExit("Datasets with zero or multiple guide libraries are not currently supported by the pipeline.")
    clo = fetch(f"{cls[0]}/@@embedded?format=json")
    guides = [f["accession"] for f in clo.get("integrated_content_files", [])
              if f.get("content_type") == "guide RNA sequences" and f.get("status") in ("in progress", "preview", "released")]
    if len(guides) != 1:
        raise SystemExit("exactly one guide RNA sequences file is required on the construct library set")
    guide_design = guides[0]
    ms_props: "dict[str, dict]" = {}
    inputs = sorted([s for s in aset.get("input_file_sets", []) if s.startswith("/measurement-sets/") or s.startswith("/auxiliary-sets/")], reverse=True)
    rows = []
    for fs in inputs:
        obj = fetch(f"{fs}/@@object?format=json")
        ftype = obj["file_set_type"]
        modality = "scRNA sequencing" if ftype == "experimental data" else ftype
        if obj["@id"].startswith("/measurement-sets/"):
            mset = obj["accession"]
            onl, meth, strand = obj.get("onlist_files", []), obj.get("onlist_method", ""), obj.get("strand_specificity", "")
            if meth and meth != "no combination":
                raise SystemExit(f"onlist_method {meth} is not supported by the pipeline")
            ms_props[obj["@id"]] = {"barcode_onlist": onl, "onlist_method": meth, "strand_specificity": strand}
        else:
            links = obj.get("measurement_sets", [])
            pr = ms_props[links[0]]
            onl, meth, strand = pr["barcode_onlist"], pr["onlist_method"], pr["strand_specificity"]
            mset = ", ".join(x.split("/")[-2] for x in links if x in inputs)
        hmap = ""
        if modality == "cell hashing barcode sequencing":
            hmap = obj.get("barcode_map", "")
            if not hmap:
                raise SystemExit(f"Missing barcode_map on auxiliary set {obj['@id']}")
            hmap = hmap.split("/")[-2]
        index: "dict[tuple, dict]" = {}
        for f in obj.get("files", []):
            fo = fetch(f"{f}/@@object?format=json")
            if (not fo["@id"].startswith("/sequence-files/") or fo.get("status") in ("deleted", "revoked")
                    or fo.get("content_type") != "reads" or fo.get("illumina_read_type") not in ("R1", "R2")):
                continue
            key = (fo.get("sequencing_run"), fo.get("lane"), fo.get("flowcell_id"), fo.get("index"))
            index.setdefault(key, {})[fo["illumina_read_type"]] = fo
        for key, reads in index.items():
            r1, r2 = reads.get("R1"), reads.get("R2")
            if not r1 or not r2:
                continue
            seqspec = ""
            for s in r1.get("seqspecs", []):
                so = fetch(f"{s}/@@object?format=json")
                if so.get("upload_status") == "validated" and so.get("status") in ("in progress", "preview", "released"):
                    seqspec = s.split("/")[-2]
            if not seqspec:
                fb = {"scrna sequencing": rna_seqspec, "grna sequencing": sgrna_seqspec, "cell hashing barcode sequencing": hash_seqspec}
                seqspec = fb.get(modality.lower().strip()) or ""
                if not seqspec:
                    raise SystemExit(f"Missing seqspec for modality {modality} (R1 {r1['@id']}); pass --rna/--sgrna/--hash-seqspec")
            if not strand or not onl or not meth:
                raise SystemExit(f"measurement set of {obj['@id']} lacks strand_specificity / onlist_files / onlist_method")
            rows.append({"R1_path": r1["@id"].split("/")[-2], "R1_md5sum": r1.get("content_md5sum"), "R2_path": r2["@id"].split("/")[-2],
                         "R2_md5sum": r2.get("content_md5sum"), "file_modality": modality, "file_set": obj["accession"], "measurement_sets": mset,
                         "sequencing_run": key[0], "lane": key[1], "flowcell_id": key[2], "index": key[3], "seqspec": seqspec,
                         "barcode_onlist": onl[0].split("/")[-2] if onl else "", "onlist_method": meth, "strand_specificity": strand,
                         "guide_design": guide_design, "barcode_hashtag_map": hmap})
    return rows


def cmd_portal_samplesheet(args) -> int:
    pd = _pd()
    out = run_dir(args.label or args.accession)
    rows = generate_per_sample(args.accession, _portal_fetch, args.hash_seqspec, args.rna_seqspec, args.sgrna_seqspec)
    df = pd.DataFrame(rows)
    write_tsv(df, out / "per_sample.tsv")
    write_tsv(portal_to_samplesheet(df), out / "samplesheet.csv", sep=",")
    finish(out, f"Portal samplesheet {args.accession}", {"accession": args.accession, "n_rows": len(df),
                                                         "modalities": df["file_modality"].value_counts().to_dict() if len(df) else {}},
           [md_table(["R1", "R2", "modality", "measurement_sets", "seqspec"], df[["R1_path", "R2_path", "file_modality", "measurement_sets", "seqspec"]].values.tolist())])
    return 0


def portal_download_plan(df, file_types: str = "all", output_dir: Path = Path("downloads")) -> "list[dict]":
    """download_igvf.py routes: sequence-files .fastq.gz, configuration-files .yaml.gz, tabular-files .tsv.gz."""
    pd = _pd()
    plan = []
    for i, r in df.iterrows():
        for col in df.columns:
            v = r[col]
            if not isinstance(v, str) or not v.startswith("IGVFFI"):
                continue
            if col in ("R1_path", "R2_path") and file_types in ("fastq", "all"):
                # download_development/download_igvf.py: R1_path -> R1_md5sum (download_pipeline/ looks up R1_path_md5sum, never found)
                route, ext, md5 = "sequence-files", "fastq.gz", r.get(f"{col.replace('_path', '')}_md5sum")
            elif col == "seqspec" and file_types in ("other", "all"):
                route, ext, md5 = "configuration-files", "yaml.gz", None
            elif col in ("barcode_onlist", "guide_design", "barcode_hashtag_map") and file_types in ("other", "all"):
                route, ext, md5 = "tabular-files", "tsv.gz", None
            else:
                continue
            plan.append({"row": i, "column": col, "accession": v, "url": f"{PORTAL}/{route}/{v}/@@download/{v}.{ext}",
                         "path": str(Path(output_dir) / f"{v}.{ext}"), "md5": md5 if isinstance(md5, str) and not pd.isna(md5) else None})
    return plan


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_portal_download(args) -> int:  # network only with --execute
    pd = _pd()
    df = pd.read_csv(args.sample, sep="\t", dtype=str)
    out = run_dir(args.label)
    dest = Path(args.output_dir) if args.output_dir else out / "downloads"
    plan = portal_download_plan(df, args.file_types, dest)
    local = df.copy()
    done = 0
    for p in plan:
        local.at[p["row"], f"local_{p['column']}"] = str(Path(p["path"]).resolve())
        if not args.execute:
            print(f"  would download {p['url']} -> {p['path']}")
            continue
        path = Path(p["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.is_file() and (p["md5"] is None or _md5(path) == p["md5"]):
            done += 1
            continue
        req = urllib.request.Request(p["url"])
        key, sec = os.environ.get("IGVF_API_KEY"), os.environ.get("IGVF_SECRET_KEY")
        if key and sec:
            req.add_header("Authorization", "Basic " + base64.b64encode(f"{key}:{sec}".encode()).decode())
        with urllib.request.urlopen(req, timeout=600) as r, open(path, "wb") as fh:
            shutil.copyfileobj(r, fh, 1 << 20)
        if p["md5"] and _md5(path) != p["md5"]:
            print(f"  FAIL  md5 mismatch for {p['accession']}")
            continue
        done += 1
    write_tsv(local, out / "per_sample_local.tsv")
    finish(out, "Portal download", {"n_files": len(plan), "downloaded": done, "executed": bool(args.execute)},
           [f"{len(plan)} files planned ({'downloaded' if args.execute else 'dry run: pass --execute to download'}). FASTQs can be large."])
    return 0


# ---------------------------------------------------------------------------
# seqspec: minimal YAML reader + `seqspec index -t kb` (parsing_guide_metadata.py, extract_parsed_seqspec.py)
# ---------------------------------------------------------------------------

def _yaml_scalar(s: str):
    s = s.strip()
    if s in ("", "null", "~", "Null", "NULL"):
        return None
    if s in ("true", "True"):
        return True
    if s in ("false", "False"):
        return False
    if (s[0] == s[-1]) and s[0] in "'\"" and len(s) >= 2:
        return s[1:-1]
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        return [_yaml_scalar(x) for x in inner.split(",")] if inner else []
    if re.fullmatch(r"[-+]?\d+", s):
        return int(s)
    if re.fullmatch(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", s):
        return float(s)
    return s


def mini_yaml(text: str):
    """Block-style YAML subset used by seqspec files (mappings, lists, !Tags ignored, scalars)."""
    lines = []
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        ind = len(raw) - len(raw.lstrip(" "))
        content = raw.strip()
        if ind == 0 and re.fullmatch(r"!\w+", content):
            continue
        lines.append([ind, content])

    def strip_tag(s: str) -> str:
        return re.sub(r"^!\w+\s*", "", s)

    def parse(i: int, ind: int):
        if i >= len(lines):
            return None, i
        if lines[i][1].startswith("- ") or lines[i][1] == "-":
            out = []
            while i < len(lines) and lines[i][0] == ind and (lines[i][1].startswith("- ") or lines[i][1] == "-"):
                item = strip_tag(lines[i][1][1:].strip())
                if not item:
                    if i + 1 < len(lines) and lines[i + 1][0] > ind:
                        v, i = parse(i + 1, lines[i + 1][0])
                    else:
                        v, i = None, i + 1
                    out.append(v)
                elif re.match(r"^[^'\"\[][^:]*:(\s|$)", item):
                    lines[i] = [ind + 2, item]
                    v, i = parse(i, ind + 2)
                    out.append(v)
                else:
                    out.append(_yaml_scalar(item))
                    i += 1
            return out, i
        out = {}
        while i < len(lines) and lines[i][0] == ind and not lines[i][1].startswith("- "):
            key, _, rest = lines[i][1].partition(":")
            rest = strip_tag(rest.strip())
            if rest == "":
                if i + 1 < len(lines) and (lines[i + 1][0] > ind or (lines[i + 1][0] == ind and lines[i + 1][1].startswith("- "))):
                    v, i = parse(i + 1, lines[i + 1][0])
                else:
                    v, i = None, i + 1
            else:
                v, i = _yaml_scalar(rest), i + 1
            out[key.strip().strip("'\"")] = v
        return out, i

    val, _ = parse(0, lines[0][0]) if lines else ({}, 0)
    return val


def load_seqspec(path: Path) -> dict:
    with _open_text(path) as fh:
        text = fh.read()
    try:
        import yaml  # type: ignore

        class _L(yaml.SafeLoader):
            pass

        def _any(loader, suffix, node):
            if isinstance(node, yaml.MappingNode):
                return loader.construct_mapping(node, deep=True)
            if isinstance(node, yaml.SequenceNode):
                return loader.construct_sequence(node, deep=True)
            return loader.construct_scalar(node)
        _L.add_multi_constructor("!", _any)
        return yaml.load(text, Loader=_L)
    except ImportError:
        return mini_yaml(text)


def _leaves(region: dict) -> "list[dict]":
    subs = region.get("regions") or []
    if not subs:
        return [region]
    out = []
    for r in subs:
        out.extend(_leaves(r))
    return out


def seqspec_index_kb(spec: dict, modality: str) -> "tuple[str, list[dict]]":
    """Emulate `seqspec index -m <modality> -t kb`: 'bcs:umi:feature' with 'file,start,stop' triplets."""
    reads = [r for r in (spec.get("sequence_spec") or []) if str(r.get("modality", "")).lower() == modality.lower()]
    if not reads:
        reads = list(spec.get("sequence_spec") or [])
    libs = spec.get("library_spec") or []
    lib = next((l for l in libs if str(l.get("region_id", "")).lower() == modality.lower() or str(l.get("region_type", "")).lower() == modality.lower()), None)
    if lib is None:
        lib = libs[0] if len(libs) == 1 else {"regions": libs}
    leaves = _leaves(lib)
    ids = [str(l.get("region_id")) for l in leaves]
    bcs, umis, feats, detail = [], [], [], []
    for fidx, rd in enumerate(reads):
        primer = str(rd.get("primer_id"))
        if primer not in ids:
            raise SystemExit(f"read {rd.get('read_id')}: primer {primer} not among library regions {ids}")
        p = ids.index(primer)
        seq = leaves[p + 1:] if str(rd.get("strand", "pos")) == "pos" else leaves[:p][::-1]
        rlen = int(rd.get("max_len") or 0)
        pos = 0
        for reg in seq:
            if pos >= rlen:
                break
            L = int(reg.get("max_len") or reg.get("min_len") or 0)
            start, stop = pos, min(pos + L, rlen)
            pos += L
            t = str(reg.get("region_type") or "").lower()
            trip = f"{fidx},{start},{stop}"
            detail.append({"read": rd.get("read_id"), "file_index": fidx, "region_id": reg.get("region_id"), "region_type": t,
                           "start": start, "stop": stop})
            if t == "barcode":
                bcs.append(trip)
            elif t == "umi":
                umis.append(trip)
            elif t in FEATURE_REGION_TYPES:
                feats.append(trip)
    if not umis:
        umis = ["-1,-1,-1"]
    if not bcs:
        bcs = ["-1"]
    return ",".join(bcs) + ":" + ",".join(umis) + ":" + ",".join(feats), detail


def spacer_chemistry(chem: str) -> str:
    """mappingGuide/main.nf: with a spacer tag, the feature triplet becomes '<file>,0,0' (search the whole read)."""
    parts = chem.split(":")
    if len(parts) < 3:
        return chem
    t = parts[2].split(",")
    parts[2] = f"{t[0] or '1'},0,0"
    return ":".join(parts)


def cmd_seqspec(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    spec = load_seqspec(Path(args.yaml))
    rows, details = [], []
    for m in args.modalities:
        rep, det = seqspec_index_kb(spec, m)
        chk = None
        if which("seqspec"):
            reads = ",".join(str(r.get("read_id")) for r in spec.get("sequence_spec") or [])
            res = subprocess.run(["seqspec", "index", "-m", m, "-t", "kb", args.yaml, "-i", reads], capture_output=True, text=True)
            chk = res.stdout.strip() if res.returncode == 0 else None
        rows.append({"modality": m, "representation": rep, "barcode_whitelist": args.whitelist or "", "seqspec_file": args.yaml,
                     "representation_spacer": spacer_chemistry(rep), "seqspec_binary": chk})
        details.extend(dict(d, modality=m) for d in det)
    df = pd.DataFrame(rows)
    for m in args.modalities:
        write_tsv(df[df["modality"] == m][["modality", "representation", "barcode_whitelist", "seqspec_file"]], out / f"{m}_parsed_seqSpec.txt")
    write_tsv(pd.DataFrame(details), out / "seqspec_regions.tsv")
    finish(out, "seqspec -> kb technology", {"representations": {r["modality"]: r["representation"] for r in rows},
                                             "assay_id": spec.get("assay_id"), "seqspec_version": spec.get("seqspec_version")},
           [md_table(["modality", "kb -x", "with spacer tag", "seqspec binary"], [[r["modality"], r["representation"], r["representation_spacer"], r["seqspec_binary"]] for r in rows]),
            md_table(["read", "region", "type", "start", "stop"], [[d["read"], d["region_id"], d["region_type"], d["start"], d["stop"]] for d in details])])
    return 0


# ---------------------------------------------------------------------------
# seqSpecCheck.py: where do guides / hashtags sit in R1 / R2?
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_sequences(path: Path, max_reads: int) -> "list[str]":
    seqs = []
    name = str(path).lower()
    with _open_text(path) as fh:
        first = fh.readline()
        if not first:
            return seqs
        if first.startswith(">") or name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
            cur = ""
            for line in [first] + list(fh):
                line = line.strip()
                if line.startswith(">"):
                    if cur:
                        seqs.append(cur)
                        if len(seqs) >= max_reads:
                            return seqs
                    cur = ""
                else:
                    cur += line
            if cur and len(seqs) < max_reads:
                seqs.append(cur)
            return seqs
        header = first
        while header and len(seqs) < max_reads:
            s = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            seqs.append(s)
            header = fh.readline()
    return seqs


def gini(a) -> float:
    np = _np()
    a = np.sort(np.asarray(a, dtype=float))
    n = a.shape[0]
    idx = np.arange(1, n + 1)
    return float(np.sum((2 * idx - n - 1) * a) / (n * np.sum(a))) if n and np.sum(a) > 0 else 0.0


def analyze_guides_in_reads(reads, guides):
    positions, upstream, hits = [], {}, Counter({g: 0 for g in guides})
    for seq in reads:
        for g in guides:
            i = seq.find(g)
            if i != -1:
                positions.append(i)
                hits[g] += 1
                frag = seq[max(0, i - 12):i]
                upstream.setdefault(i, []).append("-" * (12 - len(frag)) + frag)
                break
    return positions, upstream, hits


def seqcheck_score(total, positions, upstream, hits) -> "tuple[float, float, float, float, float]":
    if not positions:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    hit_ratio = len(positions) / total
    mode_pos, mode_count = Counter(positions).most_common(1)[0]
    pos_purity = mode_count / len(positions)
    fl = upstream.get(mode_pos, [])
    flank_purity = Counter(fl).most_common(1)[0][1] / len(fl) if fl else 0.0
    g = gini(list(hits.values()))
    return 3 * hit_ratio + 2 * pos_purity + flank_purity + (1 - g), hit_ratio, pos_purity, flank_purity, g


def seqspec_check(read1s, read2s, features, max_reads: int = 10000):
    guides_fwd = [str(x).upper() for x in dict.fromkeys(features)]
    guides_rev = [revcomp(g) for g in guides_fwd]
    stats, best = [], []
    for r1, r2 in zip(read1s, read2s):
        s1, s2 = read_sequences(Path(r1), max_reads), read_sequences(Path(r2), max_reads)
        res = {}
        for lab, seqs, gl in (("R1_Fwd", s1, guides_fwd), ("R1_Rev", s1, guides_rev), ("R2_Fwd", s2, guides_fwd), ("R2_Rev", s2, guides_rev)):
            p, u, h = analyze_guides_in_reads(seqs, gl)
            sc = seqcheck_score(max(len(seqs), 1), p, u, h)
            mode = Counter(p).most_common(1)[0][0] if p else None
            flank = Counter(u.get(mode, [])).most_common(1)[0][0] if p else ""
            res[lab] = (sc, p, flank, mode)
        winner = max(res, key=lambda k: res[k][0][0])
        bn1, bn2 = os.path.basename(str(r1)).split(".")[0], os.path.basename(str(r2)).split(".")[0]
        best.append({"Sample_R1": bn1, "Sample_R2": bn2, "Best_Config": winner, "Score": res[winner][0][0],
                     "Equation": "3*HitRatio + 2*PosPurity + 1*FlankPurity + 1*(1-Gini)", "ModePosition": res[winner][3],
                     "TopFlank": res[winner][2], "reverse_complement_guides": winner.endswith("_Rev")})
        for lab, (sc, p, flank, mode) in res.items():
            stats.append({"Sample": bn1, "Config": lab, "TotalHits": len(p), "HitRatio": sc[1], "PosPurity": sc[2], "FlankPurity": sc[3],
                          "Gini": sc[4], "FinalScore": sc[0], "IsWinner": lab == winner, "ModePosition": mode, "TopFlank": flank})
    return stats, best


def cmd_seqspec_check(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    meta_p = Path(args.metadata)
    meta = pd.read_csv(meta_p, sep="," if meta_p.suffix == ".csv" else "\t")
    col = next((c for c in ["spacer", "sequence", "guide", "seq"] if c in meta.columns), None)
    if not col:
        raise SystemExit("No guide column found in metadata")
    stats, best = seqspec_check(args.read1, args.read2, meta[col].dropna().astype(str).tolist(), args.max_reads)
    st, be = pd.DataFrame(stats), pd.DataFrame(best)
    write_tsv(st, out / "position_table.csv", sep=",")
    write_tsv(be, out / "best_config.csv", sep=",")
    figs = []
    plt = _plt()
    if plt and not args.no_plots and len(st):
        fig, ax = plt.subplots(figsize=(7, 3))
        piv = st.pivot_table(index="Sample", columns="Config", values="FinalScore")
        piv.plot(kind="bar", ax=ax, color=[BLUE, ORANGE, GREEN, RED][: piv.shape[1]])
        _style(ax, "seqSpecCheck score per read orientation", "", "score")
        figs.append(str(_save(fig, out / "seqSpec_plots" / "seqSpec_check_plots.png")))
    finish(out, "seqSpecCheck", {"best": best, "figures": figs}, [
        md_table(["sample", "config", "hits", "hit ratio", "pos purity", "flank purity", "gini", "score", "winner", "mode pos", "flank"],
                 st[["Sample", "Config", "TotalHits", "HitRatio", "PosPurity", "FlankPurity", "Gini", "FinalScore", "IsWinner", "ModePosition", "TopFlank"]].values.tolist()),
        "Recommendation: " + "; ".join(f"{b['Sample_R1']}: {b['Best_Config']} at position {b['ModePosition']} (flank `{b['TopFlank']}`) -> "
                                       f"reverse_complement_guides={str(b['reverse_complement_guides']).lower()}" for b in best)])
    return 0


# ---------------------------------------------------------------------------
# Feature references (guide_table.py, hashing_table.py, kb ref --workflow kite)
# ---------------------------------------------------------------------------

def guide_feature_table(meta, rev_comp: bool = False, spacer: str = ""):
    m = meta.copy()
    if rev_comp:
        m["spacer"] = m["spacer"].astype(str).map(revcomp)
    if spacer:
        m["spacer"] = spacer + m["spacer"].astype(str)
    return m[["spacer", "guide_id"]]


def kite_mismatch_map(features: "dict[str, str]") -> "dict[str, str]":
    """seq -> feature for exact sequences and every 1-substitution variant not shared by two features (kite)."""
    exact = {s.upper(): f for f, s in features.items()}
    var: "dict[str, set]" = defaultdict(set)
    for f, s in features.items():
        s = s.upper()
        for i, b in enumerate(s):
            for nb in "ACGT":
                if nb != b:
                    var[s[:i] + nb + s[i + 1:]].add(f)
    out = dict(exact)
    for v, fs in var.items():
        if v not in exact and len(fs) == 1:
            out[v] = next(iter(fs))
    return out


def write_kite_reference(features: "dict[str, str]", out: Path, prefix: str) -> "tuple[Path, Path]":
    fa, t2g = out / f"{prefix}_mismatch.fa", out / ("t2guide.txt" if prefix == "guide" else f"t2g_{prefix}.txt")
    mm = kite_mismatch_map(features)
    with open(fa, "w") as fh, open(t2g, "w") as th:
        for k, (seq, f) in enumerate(sorted(mm.items(), key=lambda x: (x[1], x[0]))):
            name = f if features[f].upper() == seq else f"{f}-{k}"
            fh.write(f">{name}\n{seq}\n")
            th.write(f"{name}\t{f}\t{f}\n")
    return fa, t2g


def cmd_feature_ref(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    meta = pd.read_csv(args.table, sep="\t")
    if args.kind == "guide":
        tab = guide_feature_table(meta, args.rev_comp, args.spacer or "")
        tp = out / "guide_features.txt"
        feats = dict(zip(tab["guide_id"].astype(str), tab["spacer"].astype(str)))
    else:
        tab = meta[["sequence", "hash_id"]]
        tp = out / "hashing_table.txt"
        feats = dict(zip(tab["hash_id"].astype(str), tab["sequence"].astype(str)))
    tab.to_csv(tp, sep="\t", header=False, index=False)
    print(f"Wrote: {tp}")
    fa, t2g = write_kite_reference(feats, out, "guide" if args.kind == "guide" else "hashing")
    print(f"Wrote: {fa}")
    idx = out / ("guide_index.idx" if args.kind == "guide" else "hashing_index.idx")
    cmd = ["kb", "ref", "-i", str(idx), "-f1", str(fa), "-g", str(t2g)] + (["-k", "31"] if args.spacer else []) + ["--workflow", "kite", str(tp)]
    ran, text = external(cmd, "kb", out, run=not args.no_kb)
    finish(out, f"{args.kind} feature reference", {"features": len(feats), "kite_sequences": sum(1 for _ in open(fa)) // 2, "kb_ref": text, "kb_ran": ran},
           [f"{len(feats)} features; kite 1-mismatch FASTA `{fa.name}`; `{text}` ({'ran' if ran else 'kb not on PATH'})"])
    return 0


# ---------------------------------------------------------------------------
# Mapping: kb count wrappers (mappingscRNA / mappingGuide / mappingHashing) + Python feature counting
# ---------------------------------------------------------------------------

def parse_technology(tech: str) -> "tuple[list, list, list]":
    parts = tech.split(":")

    def trips(s):
        v = [int(x) for x in s.split(",")] if s and s != "-1" else []
        return [tuple(v[i:i + 3]) for i in range(0, len(v) - 2, 3) if v[i] >= 0]
    return trips(parts[0]), trips(parts[1]) if len(parts) > 1 else [], trips(parts[2]) if len(parts) > 2 else []


def _fastq_records(path: Path):
    with _open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            yield s


def _extract(reads: "list[str]", trips) -> str:
    return "".join(reads[f][s:(e if e > 0 else None)] for f, s, e in trips)


def count_features(fastqs: "list[str]", technology: str, features: "dict[str, str]", onlist: "Optional[list[str]]" = None,
                   spacer: str = "", max_reads: int = 0):
    """kite-style feature counting in Python: barcode (onlist exact, unique 1-mismatch when the list is small),
    UMI, feature by exact / 1-mismatch lookup (whole-read scan when the feature triplet is 'f,0,0')."""
    np, sp, pd = _np(), _sp(), _pd()
    import anndata as ad  # type: ignore
    bc_t, umi_t, ft_t = parse_technology(technology)
    nfiles = 1 + max([t[0] for t in bc_t + umi_t + ft_t] or [0])
    groups = [fastqs[i:i + nfiles] for i in range(0, len(fastqs), nfiles)]
    mm = kite_mismatch_map(features)
    lens = sorted({len(s) for s in features.values()})
    scan = any(s == 0 and e == 0 for _, s, e in ft_t) or not ft_t
    onl = set(onlist) if onlist else None
    bc_fix: "dict[str, str]" = {}
    if onl is not None and len(onl) <= 100_000:
        amb: "dict[str, set]" = defaultdict(set)
        for b in onl:
            for i, c in enumerate(b):
                for nb in "ACGT":
                    if nb != c:
                        amb[b[:i] + nb + b[i + 1:]].add(b)
        bc_fix = {v: next(iter(s)) for v, s in amb.items() if len(s) == 1 and v not in onl}
    umis: "dict[tuple, set]" = defaultdict(set)
    reads_per_bc: Counter = Counter()
    n_reads = n_mapped = n_bc_on = 0
    for grp in groups:
        its = [_fastq_records(Path(p)) for p in grp]
        for recs in zip(*its):
            n_reads += 1
            if max_reads and n_reads > max_reads:
                break
            recs = list(recs)
            bc = _extract(recs, bc_t) if bc_t else "NOBC"
            if onl is not None:
                if bc not in onl:
                    bc = bc_fix.get(bc)
                    if bc is None:
                        continue
                n_bc_on += 1
            umi = _extract(recs, umi_t) if umi_t else str(n_reads)
            feat = None
            if scan:
                src = recs[ft_t[0][0] if ft_t else nfiles - 1]
                if spacer:
                    j = src.find(spacer)
                    if j >= 0:
                        for L in lens:
                            feat = mm.get(src[j + len(spacer):j + len(spacer) + L])
                            if feat:
                                break
                if feat is None:
                    for L in lens:
                        for i in range(0, len(src) - L + 1):
                            feat = mm.get(src[i:i + L])
                            if feat:
                                break
                        if feat:
                            break
            else:
                seq = _extract(recs, ft_t)
                feat = mm.get(seq)
                if feat is None:
                    for L in lens:
                        feat = mm.get(seq[:L])
                        if feat:
                            break
            reads_per_bc[bc] += 1
            if feat is None:
                continue
            n_mapped += 1
            umis[(bc, feat)].add(umi)
    names = list(features)
    fidx = {f: i for i, f in enumerate(names)}
    bcs = sorted({b for b, _ in umis})
    bidx = {b: i for i, b in enumerate(bcs)}
    rows = [bidx[b] for b, _ in umis]
    cols = [fidx[f] for _, f in umis]
    vals = [len(u) for u in umis.values()]
    X = sp.csr_matrix((np.asarray(vals, dtype=np.float32), (rows, cols)), shape=(len(bcs), len(names)))
    var = pd.DataFrame(index=pd.Index(names, name="gene_id"))
    a = ad.AnnData(X, obs=pd.DataFrame(index=pd.Index(bcs, name="barcode")), var=var)
    n_umis = int(sum(vals))
    inspect = {"numRecords": n_reads, "numReads": n_reads, "numBarcodes": len(reads_per_bc), "numUMIs": n_umis,
               "numBarcodeUMIs": len(bcs), "meanReadsPerBarcode": n_reads / max(len(reads_per_bc), 1),
               "meanUMIsPerBarcode": n_umis / max(len(bcs), 1), "numBarcodesOnOnlist": len(bcs),
               "percentageBarcodesOnOnlist": 100.0 * len(bcs) / max(len(reads_per_bc), 1),
               "numReadsOnOnlist": n_bc_on, "percentageReadsOnOnlist": 100.0 * n_bc_on / max(n_reads, 1)}
    run_info = {"n_targets": len(names), "n_processed": n_reads, "n_pseudoaligned": n_mapped, "n_unique": n_mapped,
                "p_pseudoaligned": 100.0 * n_mapped / max(n_reads, 1), "p_unique": 100.0 * n_mapped / max(n_reads, 1),
                "call": f"igvfagent crispr-pipeline map --engine python -x {technology}", "start_time": time.strftime("%c")}
    return a, inspect, run_info


def write_kb_like(a, inspect: dict, run_info: dict, outdir: Path) -> Path:
    cu = outdir / "counts_unfiltered"
    cu.mkdir(parents=True, exist_ok=True)
    write_h5ad(a, cu / "adata.h5ad", quiet=True)
    (cu / "cells_x_genes.genes.names.txt").write_text("\n".join(map(str, a.var_names)) + "\n")
    (cu / "cells_x_genes.barcodes.txt").write_text("\n".join(map(str, a.obs_names)) + "\n")
    write_json(inspect, outdir / "inspect.json", quiet=True)
    write_json(run_info, outdir / "run_info.json", quiet=True)
    print(f"Wrote: {outdir}")
    return outdir


MAP_SUFFIX = {"rna": "ks_transcripts_out", "guide": "ks_guide_out", "hash": "ks_hashing_out"}


def kb_count_command(modality: str, fastqs, index, t2g, barcodes, chem, outdir, threads: int = 4, is_10x3v3: bool = False) -> "list[str]":
    workflow = []
    if modality in ("guide", "hash"):
        workflow = ["--workflow", "kite:10xFB" if is_10x3v3 else "kite"]
        if is_10x3v3:
            chem = "10XV3"
    return ["kb", "count", "-i", str(index), "-g", str(t2g), "--verbose", "-w", str(barcodes)] + workflow + [
        "--h5ad", "-x", chem, "-o", str(outdir), "-t", str(threads)] + [str(f) for f in fastqs] + ["--overwrite"]


def cmd_map(args) -> int:
    pd = _pd()
    out = Path(args.outdir) if args.outdir else run_dir(args.label or f"map_{args.modality}_{args.batch}")
    chem = args.technology
    if not chem and args.parsed_seqspec:
        chem = "".join(pd.read_csv(args.parsed_seqspec, sep="\t")["representation"].astype(str))
    if not chem and not args.is_10x3v3:
        raise SystemExit("--technology or --parsed-seqspec (from `seqspec`) is required")
    if args.modality == "guide" and args.spacer and len(args.spacer) > 1 and not args.is_10x3v3:
        chem = spacer_chemistry(chem)
    ks = out / f"{args.batch}_{MAP_SUFFIX[args.modality]}"
    engine = args.engine
    if engine == "auto":
        engine = "kb" if which("kb") and args.index else ("python" if args.modality in ("guide", "hash") else "kb")
    summ = {"modality": args.modality, "batch": args.batch, "technology": chem, "engine": engine, "outdir": str(ks)}
    if engine == "kb":
        cmd = kb_count_command(args.modality, args.fastqs, args.index or "INDEX.idx", args.t2g or "t2g.txt", args.barcodes or "barcodes.txt",
                               chem or "10XV3", ks, args.threads, args.is_10x3v3)
        ran, text = external(cmd, "kb")
        summ.update(kb_ran=ran, command=text)
        if not ran and args.modality == "rna":
            print("RNA pseudoalignment needs kb/kallisto; no Python fallback for the transcriptome (use `concat` on existing kb outputs).")
    else:
        if args.modality == "rna":
            raise SystemExit("--engine python counts guide / hash feature libraries only; RNA needs kb (kallisto | bustools)")
        meta = pd.read_csv(args.features, sep="\t")
        if args.modality == "guide":
            tab = guide_feature_table(meta, args.rev_comp, "")
            feats = dict(zip(tab["guide_id"].astype(str), tab["spacer"].astype(str)))
        else:
            feats = dict(zip(meta["hash_id"].astype(str), meta["sequence"].astype(str)))
        onlist = None
        if args.barcodes:
            with _open_text(Path(args.barcodes)) as fh:
                onlist = [l.strip().split("\t")[0] for l in fh if l.strip()]
        a, insp, ri = count_features(args.fastqs, chem, feats, onlist, args.spacer or "", args.max_reads)
        write_kb_like(a, insp, ri, ks)
        summ.update(n_barcodes=int(a.n_obs), n_features=int(a.n_vars), n_umis=int(a.X.sum()), p_pseudoaligned=ri["p_pseudoaligned"])
    if not args.outdir:
        finish(out, f"Mapping {args.modality} {args.batch}", summ, [md_table(["key", "value"], [[k, v] for k, v in summ.items()])])
    else:
        print(json.dumps(_jsonable(summ)))
    return 0


# ---------------------------------------------------------------------------
# anndata_concat.py / hashing_concat.py
# ---------------------------------------------------------------------------

def concat_ks_dirs(dirs: "list[str]", covariates=None):
    """Sort by basename; batch = '<batch>' from '<batch>_ks_*'; join covariate columns; outer join, index_unique '_'."""
    import anndata as ad  # type: ignore
    pd = _pd()
    adatas, var_name = [], None
    for i, d in enumerate(sorted(dirs, key=lambda x: os.path.basename(os.path.normpath(x)))):
        p = Path(d)
        h5 = p / "counts_unfiltered" / "adata.h5ad" if p.is_dir() else p
        a = read_h5ad(h5)
        if i == 0 and a.var_names.name is not None:
            var_name = a.var_names.name
        m = re.search(r"(.+)_ks_", os.path.basename(os.path.normpath(d)))
        if m:
            a.obs["batch"] = m.group(1)
            if covariates is not None and len(covariates.columns) > 1:
                cov = covariates.astype(str).set_index(covariates.columns[0])
                a.obs = a.obs.join(cov, on="batch")
        adatas.append(a)
    comb = ad.concat(adatas, join="outer", index_unique="_") if len(adatas) > 1 else adatas[0]
    if len(adatas) == 1:
        comb = comb.copy()
        comb.obs_names = [f"{n}_0" for n in comb.obs_names]
    if var_name:
        comb.var_names.name = var_name
    comb.X = as_csr(comb.X)
    return comb


def cmd_concat(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    cov = pd.read_csv(args.covariates) if args.covariates else None
    a = concat_ks_dirs(args.inputs, cov)
    p = write_h5ad(a, out / (args.output or "concatenated_adata.h5ad"))
    finish(out, "Concatenated AnnData", {"n_obs": int(a.n_obs), "n_vars": int(a.n_vars), "batches": a.obs["batch"].value_counts().to_dict() if "batch" in a.obs else {},
                                         "output": str(p)}, [f"{a.n_obs} barcodes x {a.n_vars} features from {len(args.inputs)} inputs"])
    return 0


# ---------------------------------------------------------------------------
# preprocess_adata.py: knee barcode filter + QC metrics + filters
# ---------------------------------------------------------------------------

def _knee_vectors(x, y):
    np = _np()
    from scipy.interpolate import UnivariateSpline  # type: ignore
    d2 = UnivariateSpline(x, y, s=len(x)).derivative(n=2)(x)
    ten = max(1, round(len(x) * 0.1))
    mid = d2[ten:-ten] if len(d2) > 2 * ten else d2
    if np.all(mid >= 0) or np.all(mid <= 0):
        return x, y
    amin = int(np.argmin(d2))
    c1 = np.where(d2[:amin + 1] >= 0)[0]
    if len(c1) == 0:
        return x, y
    c2 = np.where(d2[amin:] >= 0)[0]
    if len(c2) == 0:
        return x, y
    e1, e2 = c1[-1], amin + c2[0]
    if e1 >= e2:
        return x, y
    return x[e1:e2 + 1], y[e1:e2 + 1]


def elbow_knee_finder(x, y, mode: str = "basic"):
    np = _np()
    if mode == "advanced":
        if len(np.unique(x)) < 4:
            return None
        x, y = _knee_vectors(x, y)
    if len(x) == 0 or x[0] == x[-1]:
        return None
    slope = (y[-1] - y[0]) / (x[-1] - x[0])
    icpt = y[0] - slope * x[0]
    dist = np.abs(slope * x - y + icpt) / np.sqrt(slope ** 2 + 1)
    k = int(np.argmax(dist))
    return np.array([x[k], y[k]])


def knee_points(counts):
    """Barcode rank vs log1p(count): (knee1, knee2, sorted counts)."""
    np = _np()
    s = np.sort(np.asarray(counts, dtype=float))[::-1]
    x, y = np.arange(1, len(s) + 1, dtype=float), np.log1p(s)
    p1 = elbow_knee_finder(x, y, "basic")
    p2 = None
    if p1 is not None:
        end = max(1, min(len(x), int(round(p1[0]))))
        p2 = elbow_knee_finder(x[:end], y[:end], "advanced")
    return p1, p2, s


def barcode_filter_threshold(counts, method: str) -> "tuple[Optional[float], Optional[int]]":
    p1, p2, s = knee_points(counts)
    pt = p1 if method == "knee" else p2
    if pt is None:
        return None, None
    rank = max(1, min(len(s), int(round(pt[0]))))
    return float(s[rank - 1]), rank


def _gene_symbols(a):
    for c in ("symbol", "gene_name", "gene_symbol"):
        if c in a.var.columns:
            return a.var[c].astype(str).values
    return a.var_names.astype(str).values


def preprocess_rna(a, gene_names: "Optional[list[str]]", min_genes: int, pct_mito: float, reference: str, barcode_filter: str,
                   figdir: Optional[Path] = None, compat: bool = False, plots: bool = True) -> "tuple[Any, Any, dict]":
    np, pd = _np(), _pd()
    a = a.copy()
    a.X = as_csr(a.X)
    if gene_names is not None:
        if len(gene_names) != a.n_vars:
            raise SystemExit("The number of gene names does not match the number of variables in adata_rna")
        a.var["symbol"] = list(gene_names)
    elif "symbol" not in a.var.columns:
        a.var["symbol"] = _gene_symbols(a)
    a.var_names = pd.Index(a.var_names.astype(str)).str.split(".").str[0]
    a.var_names_make_unique()
    counts = row_sums(a.X)
    info = {"n_barcodes_in": int(a.n_obs), "barcode_filter": barcode_filter}
    p1, p2, s = knee_points(counts)
    figs = []
    plt = _plt() if plots else None
    if plt and figdir:
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(np.arange(1, len(s) + 1), np.log1p(s), color=BLUE, lw=1.2)
        for pt, col, lab in ((p1, RED, "Knee 1"), (p2, ORANGE, "Knee 2")):
            if pt is not None:
                ax.axvline(pt[0], color=col, ls="--", lw=1, label=lab)
        ax.set_xscale("log")
        ax.legend(frameon=False, fontsize=8)
        _style(ax, "Barcode rank plot", "barcode rank", "log1p UMI counts")
        figs.append(str(_save(fig, figdir / "knee_plot_scRNA.png")))
    if barcode_filter != "none":
        thr, rank = barcode_filter_threshold(counts, barcode_filter)
        info.update(knee_rank=rank, count_threshold=thr)
        if thr is not None:
            a = a[counts >= thr, :].copy()
    info["n_after_barcode_filter"] = int(a.n_obs)
    if "batch" in a.obs.columns:
        a.obs["batch_number"] = pd.factorize(a.obs["batch"])[0] + 1
    sym = a.var["symbol"].astype(str)
    if compat or reference == "human":
        a.var["mt"] = sym.str.startswith("MT-").values
    else:
        a.var["mt"] = sym.str.upper().str.startswith("MT-").values
    a.var["ribo"] = sym.str.startswith(("RPS", "RPL")).values
    a.X = a.X.astype(np.float32)
    try:
        import scanpy as sc  # type: ignore
        pt = [x for x in (50, 100, 200, 500) if x < a.n_vars] or None
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=pt)
    except ImportError:
        tot = row_sums(a.X)
        a.obs["n_genes_by_counts"] = nnz_rows(a.X)
        a.obs["total_counts"] = tot
        for q in ("mt", "ribo"):
            sub = row_sums(a.X[:, a.var[q].values])
            a.obs[f"total_counts_{q}"] = sub
            a.obs[f"pct_counts_{q}"] = np.where(tot > 0, 100 * sub / np.maximum(tot, 1), 0)
    if plt and figdir and a.n_obs:
        fig, axes = plt.subplots(1, 3, figsize=(9, 3))
        for ax, c in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
            ax.violinplot(a.obs[c].values, showmedians=True)
            _style(ax, c, "", "")
        figs.append(str(_save(fig, figdir / "violinplot_scrna.png")))
        fig, ax = plt.subplots(figsize=(5, 4))
        sc_ = ax.scatter(a.obs["total_counts"], a.obs["n_genes_by_counts"], c=a.obs["pct_counts_mt"], s=3, cmap="viridis")
        fig.colorbar(sc_, ax=ax, label="pct_counts_mt")
        _style(ax, "total counts vs genes", "total_counts", "n_genes_by_counts")
        figs.append(str(_save(fig, figdir / "scatterplot_scrna.png")))
    if barcode_filter == "none":
        ng = nnz_rows(a.X)
        a = a[ng >= min_genes, :].copy()
        a.obs["n_genes"] = nnz_rows(a.X)
    info["n_after_min_genes"] = int(a.n_obs)
    ncell = np.asarray((a.X > 0).sum(axis=0)).ravel()
    a = a[:, ncell >= 10].copy()
    a.var["n_cells"] = np.asarray((a.X > 0).sum(axis=0)).ravel()
    unfiltered = a.copy()
    a = a[a.obs["pct_counts_mt"].values < pct_mito, :].copy()
    info.update(n_after_mito=int(a.n_obs), n_genes=int(a.n_vars), figures=figs)
    return a, unfiltered, info


def cmd_preprocess(args) -> int:
    out = run_dir(args.label)
    a = read_h5ad(Path(args.adata))
    names = None
    if args.gene_names:
        gp = Path(args.gene_names)
        gp = gp / "counts_unfiltered" / "cells_x_genes.genes.names.txt" if gp.is_dir() else gp
        names = [l.strip() for l in open(gp) if l.strip()]
    f, _, info = preprocess_rna(a, names, args.min_genes, args.pct_mito, args.reference, args.barcode_filter, out / "figures",
                                args.upstream_compat, not args.no_plots)
    p = write_h5ad(f, out / "filtered_anndata.h5ad")
    info["output"] = str(p)
    finish(out, "scRNA preprocessing", info, [md_table(["step", "cells"], [["input barcodes", info["n_barcodes_in"]],
                                                                            [f"barcode filter ({args.barcode_filter})", info["n_after_barcode_filter"]],
                                                                            [f"min genes ({args.min_genes})", info["n_after_min_genes"]],
                                                                            [f"pct_counts_mt < {args.pct_mito}", info["n_after_mito"]]])])
    return 0


# ---------------------------------------------------------------------------
# create_mdata.py
# ---------------------------------------------------------------------------

def merge_guide_metadata(guide_var, meta):
    overl = [c for c in meta.columns if c != "guide_id" and c in guide_var.columns]
    gv = guide_var.drop(columns=overl)
    gv = gv.assign(guide_id=gv["guide_id"].astype(str))
    meta = meta.assign(guide_id=meta["guide_id"].astype(str))
    m = gv.merge(meta, on="guide_id", how="left", validate="one_to_one")
    if m["guide_id"].isna().any():
        raise SystemExit("guide_id is missing after merging guide metadata")
    return m


def read_guide_metadata(path: Path):
    pd = _pd()
    meta = pd.read_csv(path, sep="\t")
    if "guide_id" in meta.columns:
        meta = meta.dropna(subset=["guide_id"])
    if "spacer" in meta.columns:
        meta = meta.dropna(subset=["spacer"])
    for c in ("guide_start", "guide_end", "intended_target_start", "intended_target_end"):
        if c in meta.columns:
            meta[c] = pd.to_numeric(meta[c], errors="coerce")
    if not meta["guide_id"].is_unique:
        raise SystemExit("guide_id in guide metadata is not unique.")
    return meta


def create_mudata(rna, guide, meta, gtf=None, moi: str = "high", capture_method: str = "CROP-seq", hashing=None,
                  figdir: Optional[Path] = None, plots: bool = True) -> MData:
    np, pd = _np(), _pd()
    guide = guide.copy()
    gv = guide.var.reset_index()
    gv = gv.rename(columns={"gene_id": "guide_id", "feature_id": "guide_id"})
    if "guide_id" not in gv.columns:
        gv = gv.rename(columns={gv.columns[0]: "guide_id"})
    if len(meta) != len(gv):
        print(f"The numbers of guide_id differ: {len(gv)} in guide anndata, {len(meta)} in guide metadata.")
    gv = merge_guide_metadata(gv, meta)
    targ = gv["targeting"].map(_truthy) if "targeting" in gv.columns else pd.Series(True, index=gv.index)
    if "type" not in gv.columns:
        gv["type"] = np.where(targ, "targeting", "non-targeting")
    else:
        miss = gv["type"].isna()
        gv.loc[miss, "type"] = np.where(targ[miss], "targeting", "non-targeting")
    gv = annotate_intended_target_groups(gv)
    gv.index = pd.Index(gv["guide_id"].astype(str))
    guide.var = gv
    guide.var_names = pd.Index(gv["guide_id"].astype(str))
    guide.X = as_csr(guide.X)
    guide.uns["capture_method"] = np.array([capture_method], dtype=object)
    if moi in ("high", "low"):
        guide.uns["moi"] = np.array([moi], dtype=object)
    else:
        avg = float(np.mean(nnz_rows(guide.X)))
        guide.uns["moi"] = np.array(["high" if avg > 1.5 else "low"], dtype=object)
    guide.obs["num_expressed_guides"] = nnz_rows(guide.X)
    guide.obs["total_guide_umis"] = row_sums(guide.X)
    if hashing is not None and "batch" in guide.obs.columns:
        guide.obs["batch_number"] = pd.factorize(guide.obs["batch"])[0] + 1
    rna = rna.copy()
    if gtf is not None and len(gtf):
        g = gtf.drop_duplicates("gene_id2").set_index("gene_id2")[["seqname", "start", "end"]].rename(
            columns={"seqname": "gene_chr", "start": "gene_start", "end": "gene_end"})
        rna.var = rna.var.drop(columns=[c for c in ("gene_chr", "gene_start", "gene_end") if c in rna.var.columns]).join(g)
    rna.obs = rna.obs.rename(columns={"n_genes_by_counts": "n_counts", "pct_counts_mt": "percent_mito", "n_genes": "num_expressed_genes",
                                      "total_counts": "total_gene_umis"})
    if "num_expressed_genes" not in rna.obs.columns and "n_counts" in rna.obs.columns:
        rna.obs["num_expressed_genes"] = rna.obs["n_counts"].values  # n_genes_by_counts; upstream leaves it absent after a knee filter
    plt = _plt() if plots else None
    if plt and figdir:
        s = np.sort(row_sums(guide.X))[::-1]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(np.arange(len(s)), np.log1p(s), color=ORANGE, lw=1.2)
        _style(ax, "Knee Plot (guide UMIs)", "Barcode Index", "Log of UMI Counts")
        _save(fig, figdir / "knee_plot_guide.png")
    gset = set(guide.obs_names)
    cells = [c for c in rna.obs_names if c in gset]
    if hashing is not None:
        hset = set(hashing.obs_names)
        cells = [c for c in cells if c in hset]
    mods = {"gene": rna[cells, :].copy(), "guide": guide[cells, :].copy()}
    if hashing is not None:
        mods["hashing"] = hashing[cells, :].copy()
    common = [c for c in mods["guide"].obs.columns if c in set(mods["gene"].obs.columns)]
    if hashing is not None:
        common = [c for c in common if c in set(mods["hashing"].obs.columns)]
    m = MData(mods, mods["guide"].obs.loc[:, common].copy())
    m.uns["mudata_created_by"] = "igvfagent crispr-pipeline create-mudata"
    return m


def cmd_create_mudata(args) -> int:
    out = run_dir(args.label)
    gtf = read_gtf(Path(args.gtf)) if args.gtf else None
    hashing = read_h5ad(Path(args.hashing)) if args.hashing else None
    m = create_mudata(read_h5ad(Path(args.rna)), read_h5ad(Path(args.guide)), read_guide_metadata(Path(args.guide_metadata)), gtf,
                      args.moi, args.capture_method, hashing, out / "figures", not args.no_plots)
    p = write_h5mu(m, out / "mudata.h5mu")
    gv = m["guide"].var
    finish(out, "MuData", {"n_cells": m.n_obs, "n_genes": int(m["gene"].n_vars), "n_guides": int(m["guide"].n_vars), "moi": uns_str(MData({"g": m["guide"]}, uns=dict(m["guide"].uns)), "moi"),
                           "n_control_guides": int(control_mask(gv).sum()), "output": str(p)},
           [md_table(["guide_id", "type", "intended_target_name", "intended_target_key"], gv[["guide_id", "type", "intended_target_name", "intended_target_key"]].head(20).values.tolist())])
    return 0


# ---------------------------------------------------------------------------
# Hashing: filter_hashing.py, demultiplex (GMM-demux) + demultiplex_filter.py, hashing_concat.py
# ---------------------------------------------------------------------------

def gmm_demux_fallback(a) -> "tuple[Any, Any]":
    """Per-HTO two-component Gaussian mixture on log1p counts -> positive calls; cluster = set of positive HTOs.

    Returns (report: index barcode, column Cluster_id; config: cluster_id, hto_type) in GMM-demux's FULL format.
    """
    np, pd = _np(), _pd()
    from sklearn.mixture import GaussianMixture  # type: ignore
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    pos = np.zeros(X.shape, dtype=bool)
    for j in range(X.shape[1]):
        x = np.log1p(X[:, j]).reshape(-1, 1)
        if np.unique(x).size < 3:
            pos[:, j] = X[:, j] > 0
            continue
        gm = GaussianMixture(2, random_state=0).fit(x)
        hi = int(np.argmax(gm.means_.ravel()))
        pos[:, j] = gm.predict_proba(x)[:, hi] > 0.5
    names = list(map(str, a.var_names))
    labels = ["-".join(names[j] for j in np.where(r)[0]) or "negative" for r in pos]
    uniq = ["negative"] + sorted(set(labels) - {"negative"}, key=lambda s: (s.count("-"), s))
    cid = {l: i for i, l in enumerate(uniq)}
    report = pd.DataFrame({"Cluster_id": [cid[l] for l in labels]}, index=a.obs_names)
    config = pd.DataFrame({"cluster_id": list(range(len(uniq))), "hto_type": uniq})
    return report, config


def apply_demux(a, report, config):
    """demultiplex_filter.py: hto_type / hto_type_split ('multiplets' when several HTOs); drop negative + multiplets."""
    a = a.copy()
    config = config.copy()
    config.columns = ["cluster_id", "hto_type"]
    config["hto_type"] = config["hto_type"].astype(str).str.strip()
    a.obs["cluster_id"] = report["Cluster_id"].reindex(a.obs_names).values
    config["hto_type_split"] = config["hto_type"].str.split("-").str.join(",").map(lambda x: "multiplets" if "," in x else x)
    a.obs = a.obs.merge(config, how="left", on="cluster_id").set_index(a.obs.index)
    filt = a[(a.obs["hto_type"] != "negative").values & (a.obs["hto_type_split"] != "multiplets").values].copy()
    return filt, a


def demultiplex_batch(a, workdir: Path):
    pd = _pd()
    if which("GMM-demux"):
        d = workdir / "demuxfile"
        d.mkdir(parents=True, exist_ok=True)
        import scipy.io  # type: ignore
        with gzip.open(d / "barcodes.tsv.gz", "wt") as fh:
            fh.write("\n".join(a.obs_names) + "\n")
        with gzip.open(d / "features.tsv.gz", "wt") as fh:
            fh.write("\n".join(a.var_names) + "\n")
        buf = io.BytesIO()
        scipy.io.mmwrite(buf, as_csr(a.X).T)
        with gzip.open(d / "matrix.mtx.gz", "wb") as fh:
            fh.write(buf.getvalue())
        external(["GMM-demux", str(d), ",".join(a.var_names), "-f", str(workdir / "FULL")], "GMM-demux", workdir)
        report = pd.read_csv(workdir / "FULL" / "GMM_full.csv", index_col=0)
        config = pd.read_csv(workdir / "FULL" / "GMM_full.config", header=None)
        engine = "GMM-demux"
    else:
        report, config = gmm_demux_fallback(a)
        engine = "python-gmm"
    f, u = apply_demux(a, report, config)
    return f, u, engine


def run_hashing(hash_adata, rna_filtered, workdir: Path) -> "tuple[Any, Any, dict]":
    import anndata as ad  # type: ignore
    keep = [c for c in hash_adata.obs_names if c in set(rna_filtered.obs_names)]
    h = hash_adata[keep, :].copy()
    fs, us, engines = [], [], set()
    for b in sorted(h.obs["batch"].astype(str).unique()) if "batch" in h.obs.columns else ["all"]:
        sub = h[h.obs["batch"].astype(str).values == b].copy() if b != "all" else h
        f, u, eng = demultiplex_batch(sub, workdir / f"demux_{safe_label(b)}")
        fs.append(f)
        us.append(u)
        engines.add(eng)
    filt = ad.concat(fs, join="outer") if len(fs) > 1 else fs[0]
    unf = ad.concat(us, join="outer") if len(us) > 1 else us[0]
    info = {"n_intersecting": len(keep), "n_singlets": int(filt.n_obs), "engine": sorted(engines),
            "hto_type_split": unf.obs["hto_type_split"].value_counts().to_dict()}
    return filt, unf, info


def cmd_hashing(args) -> int:
    out = run_dir(args.label)
    f, u, info = run_hashing(read_h5ad(Path(args.hashing)), read_h5ad(Path(args.rna)), out)
    write_h5ad(f, out / "concatenated_hashing_demux.h5ad")
    write_h5ad(u, out / "concatenated_unfiltered_hashing_demux.h5ad")
    finish(out, "Hashing demultiplexing", info, [md_table(["hto_type_split", "cells"], sorted(info["hto_type_split"].items()))])
    return 0


# ---------------------------------------------------------------------------
# doublets.py (scrublet) with a Python fallback
# ---------------------------------------------------------------------------

def threshold_minimum(values, nbins: int = 256) -> float:
    """skimage-style threshold_minimum: smooth the histogram until bimodal, return the valley (Otsu if it never is)."""
    np = _np()
    v = np.asarray(values, dtype=float)
    hist, edges = np.histogram(v, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    h = hist.astype(float)
    for _ in range(10000):
        h = np.convolve(h, np.ones(3) / 3, mode="same")
        peaks = [i for i in range(1, len(h) - 1) if h[i - 1] < h[i] > h[i + 1]]
        if len(peaks) == 2:
            lo, hi = peaks
            return float(centers[lo + int(np.argmin(h[lo:hi + 1]))])
        if len(peaks) < 2:
            break
    w = np.cumsum(hist)
    mu = np.cumsum(hist * centers)
    tot, mt = w[-1], mu[-1]
    between = np.where((w > 0) & (w < tot), (mt * w / tot - mu) ** 2 / np.maximum(w * (tot - w), 1e-12), 0)
    return float(centers[int(np.argmax(between))])


def scrublet_fallback(X, expected_rate: float = 0.1, sim_ratio: float = 2.0, n_pcs: int = 30, seed: int = 42):
    """Scrublet's method: simulate doublets from random cell pairs, PCA, kNN, Bayesian doublet score, auto threshold."""
    np, sp = _np(), _sp()
    from sklearn.decomposition import PCA  # type: ignore
    from sklearn.neighbors import NearestNeighbors  # type: ignore
    rng = np.random.default_rng(seed)
    X = as_csr(X)
    n = X.shape[0]
    n_sim = int(round(sim_ratio * n))
    i1, i2 = rng.integers(0, n, n_sim), rng.integers(0, n, n_sim)
    S = X[i1] + X[i2]
    keep = np.asarray(((X >= 3).sum(axis=0))).ravel() >= 3
    if keep.sum() < 10:
        keep = np.asarray((X > 0).sum(axis=0)).ravel() > 0
    Xo, Xs = X[:, keep], S[:, keep]
    tgt = float(np.mean(row_sums(Xo))) or 1.0

    def norm(M):
        tot = row_sums(M)
        return np.log1p(sp.diags(tgt / np.maximum(tot, 1)).dot(M).toarray())
    No, Ns = norm(Xo), norm(Xs)
    var = No.var(axis=0)
    top = np.argsort(var)[::-1][: min(2000, No.shape[1])]
    mu, sd = No[:, top].mean(axis=0), No[:, top].std(axis=0) + 1e-8
    Zo, Zs = (No[:, top] - mu) / sd, (Ns[:, top] - mu) / sd
    pca = PCA(n_components=min(n_pcs, Zo.shape[1] - 1, n - 1), random_state=seed).fit(Zo)
    Po, Ps = pca.transform(Zo), pca.transform(Zs)
    P = np.vstack([Po, Ps])
    is_sim = np.r_[np.zeros(n, bool), np.ones(n_sim, bool)]
    k = max(1, int(round(0.5 * np.sqrt(n))))
    r = n_sim / n
    k_adj = int(round(k * (1 + r)))
    nn = NearestNeighbors(n_neighbors=min(k_adj + 1, len(P))).fit(P)
    _, idx = nn.kneighbors(P)
    nsim = is_sim[idx[:, 1:]].sum(axis=1)
    q = (nsim + 1) / (k_adj + 2)
    rho = expected_rate
    ld = q * rho / r / ((1 - rho) - q * (1 - rho - rho / r))
    obs_s, sim_s = ld[:n], ld[n:]
    thr = threshold_minimum(sim_s)
    return obs_s, obs_s > thr, thr, sim_s


def remove_doublets(m: MData, figdir: Optional[Path] = None, plots: bool = True) -> "tuple[MData, dict]":
    np = _np()
    g = m["gene"].copy()
    try:
        import scrublet as scr  # type: ignore
        s = scr.Scrublet(g.X, random_state=42)
        scores, pred = s.scrub_doublets()
        thr, sim, engine = float(s.threshold_), s.doublet_scores_sim_, "scrublet"
    except ImportError:
        scores, pred, thr, sim = scrublet_fallback(g.X)
        engine = "python-scrublet"
    pred = np.asarray(pred, dtype=bool)
    g.obs["doublet_scores"], g.obs["predicted_doublets"] = scores, pred
    g.obs["doublet_info"] = g.obs["predicted_doublets"].astype(str)
    plt = _plt() if plots else None
    if plt and figdir:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3))
        axes[0].hist(scores, bins=50, color=BLUE)
        axes[0].axvline(thr, color=RED, ls="--")
        _style(axes[0], "Observed transcriptomes", "doublet score", "cells")
        axes[1].hist(sim, bins=50, color=ORANGE)
        axes[1].axvline(thr, color=RED, ls="--")
        _style(axes[1], "Simulated doublets", "doublet score", "")
        _save(fig, figdir / "doublets_batch.png")
    g = g[~pred, :].copy()
    cells = [c for c in g.obs_names if c in set(m["guide"].obs_names)]
    mods = {"gene": g[cells, :].copy(), "guide": m["guide"][cells, :].copy()}
    common = [c for c in mods["guide"].obs.columns if c in set(mods["gene"].obs.columns)]
    out = MData(mods, mods["guide"].obs.loc[:, common].copy(), dict(m.uns))
    return out, {"engine": engine, "threshold": thr, "n_doublets": int(pred.sum()), "n_cells_after": out.n_obs}


def cmd_doublets(args) -> int:
    out = run_dir(args.label)
    m, info = remove_doublets(read_h5mu(Path(args.mudata)), out / "figures", not args.no_plots)
    write_h5mu(m, out / "mdata_doublets.h5mu")
    finish(out, "Doublet removal", info, [f"{info['engine']}: {info['n_doublets']} predicted doublets removed (threshold {info['threshold']:.3f}); {info['n_cells_after']} cells remain"])
    return 0


# ---------------------------------------------------------------------------
# Guide assignment (prepare_assignment, assign_grnas_sceptre.R, cleanser, mudata_concat, collapse_guides)
# ---------------------------------------------------------------------------

def poisson_irls(X, y, offset=None, n_iter: int = 25, w=None):
    """Poisson GLM by IRLS (log link)."""
    np = _np()
    n, p = X.shape
    off = np.zeros(n) if offset is None else offset
    w0 = np.ones(n) if w is None else w
    beta = np.zeros(p)
    beta[0] = np.log(max((w0 * y).sum() / max((w0 * np.exp(off)).sum(), 1e-12), 1e-8))
    for _ in range(n_iter):
        eta = X @ beta + off
        mu = np.exp(np.clip(eta, -30, 30))
        z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        W = w0 * mu
        XtW = X.T * W
        try:
            nb = np.linalg.solve(XtW @ X + 1e-8 * np.eye(p), XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(nb - beta)) < 1e-8:
            beta = nb
            break
        beta = nb
    return beta, np.exp(np.clip(X @ beta + off, -30, 30))


def sceptre_mixture_assign(y, covars, threshold: float = 0.8, n_em_rep: int = 5, seed: int = 0, n_iter: int = 100):
    """SCEPTRE 'mixture' gRNA assignment: Poisson GLM on covariates for the pilot means, then a two-component Poisson
    mixture with those log-means as offsets (rates exp(o + g0) / exp(o + g1), weight pi) fitted by EM with
    n_em_rep random restarts; a cell is assigned when P(perturbed) >= threshold and the count is > 0."""
    np = _np()
    from scipy.special import gammaln  # type: ignore
    y = np.asarray(y, dtype=float)
    if y.sum() == 0 or (y > 0).sum() < 2:
        return np.zeros(len(y), bool), np.zeros(len(y))
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in covars])
    _, mu = poisson_irls(X, y)
    o = np.log(np.maximum(mu, 1e-12))
    eo = np.exp(o)
    lgy = gammaln(y + 1)
    rng = np.random.default_rng(seed)
    best, best_ll = None, -np.inf
    for rep in range(max(1, n_em_rep)):
        pi = rng.uniform(0.01, 0.3)
        g0 = np.log(max(y.sum() / eo.sum(), 1e-8)) - rng.uniform(0.5, 2)
        g1 = g0 + rng.uniform(2, 6)
        ll_old = -np.inf
        for _ in range(n_iter):
            l0 = np.log(1 - pi) + y * (o + g0) - eo * np.exp(g0) - lgy
            l1 = np.log(pi) + y * (o + g1) - eo * np.exp(g1) - lgy
            mx = np.maximum(l0, l1)
            lse = mx + np.log(np.exp(l0 - mx) + np.exp(l1 - mx))
            r1 = np.exp(l1 - lse)
            ll = lse.sum()
            pi = float(np.clip(r1.mean(), 1e-6, 1 - 1e-6))
            g0 = np.log(max(((1 - r1) * y).sum(), 1e-10) / max(((1 - r1) * eo).sum(), 1e-12))
            g1 = np.log(max((r1 * y).sum(), 1e-10) / max((r1 * eo).sum(), 1e-12))
            if abs(ll - ll_old) < 1e-6 * max(1.0, abs(ll)):
                break
            ll_old = ll
        if g1 < g0:
            r1 = 1 - r1
        if ll > best_ll:
            best_ll, best = ll, r1
    post = best
    return (post >= threshold) & (y > 0), post


def cleanser_fallback_assign(y, total, threshold: float = 0.99, n_iter: int = 200):
    """CLEANSER-style: ambient Poisson (rate a * cell's total guide UMIs) vs captured-guide negative binomial, EM."""
    np = _np()
    from scipy.special import gammaln  # type: ignore
    y = np.asarray(y, float)
    T = np.maximum(np.asarray(total, float), 1.0)
    if y.sum() == 0:
        return np.zeros(len(y), bool), np.zeros(len(y))
    pos = y[y > 0]
    cut = np.quantile(pos, 0.5) if len(pos) else 1
    r = (y >= max(cut, 2)).astype(float)
    lgy = gammaln(y + 1)
    for _ in range(n_iter):
        pi = float(np.clip(r.mean(), 1e-6, 1 - 1e-6))
        a = max(((1 - r) * y).sum(), 1e-10) / max(((1 - r) * T).sum(), 1e-10)
        mu = max((r * y).sum() / max(r.sum(), 1e-10), 1e-6)
        var = (r * (y - mu) ** 2).sum() / max(r.sum(), 1e-10)
        th = mu ** 2 / (var - mu) if var > mu * 1.01 else 1e3
        lam = a * T
        lb = np.log(1 - pi) + y * np.log(lam) - lam - lgy
        ls = (np.log(pi) + gammaln(y + th) - gammaln(th) - lgy + th * np.log(th / (th + mu)) + y * np.log(mu / (th + mu)))
        mx = np.maximum(lb, ls)
        rn = np.exp(ls - (mx + np.log(np.exp(lb - mx) + np.exp(ls - mx))))
        if np.max(np.abs(rn - r)) < 1e-7:
            r = rn
            break
        r = rn
    return (r >= min(threshold, 0.999)) & (y > 0), r


def _external_assignment(m: MData, rows, method: str, threshold: Optional[float], n_em_rep: int, capture_method: str,
                         upstream_bin: Optional[str], workdir: Path):
    """guide_assignment_cleanser / guide_assignment_sceptre modules on one batch; None when the tool is unavailable."""
    import scipy.io  # type: ignore
    workdir.mkdir(parents=True, exist_ok=True)
    sub = m.subset_cells(list(m["guide"].obs_names[rows]))
    inp = write_h5mu(sub, workdir / "batch_mudata.h5mu", quiet=True)
    if method == "cleanser":
        cmd = ["cleanser", "-i", str(inp), "--posteriors-output", str(workdir / "batch_mudata_output.h5mu"), "--modality", "guide",
               f"--{capture_method}", "--output-layer", "guide_assignment"] + (["-t", str(threshold)] if threshold else [])
        ran, _ = external(cmd, "cleanser", workdir)
        if not ran:
            return None
        return as_csr(read_h5mu(workdir / "batch_mudata_output.h5mu")["guide"].layers["guide_assignment"])
    script = Path(upstream_bin) / "assign_grnas_sceptre.R" if upstream_bin else Path("assign_grnas_sceptre.R")
    ran, _ = external(["Rscript", str(script), str(inp), str(threshold if threshold is not None else 0.8), str(n_em_rep)], "Rscript", workdir,
                      run=bool(upstream_bin) and script.is_file())
    if not ran:
        return None
    return as_csr(scipy.io.mmread(str(workdir / "guide_assignment.mtx")).T)  # add_guide_assignment.py: transpose to cells x guides


def assign_guides(m: MData, method: str = "sceptre", threshold: Optional[float] = None, n_em_rep: int = 5, umi_threshold: int = 5,
                  batch_key: str = "batch", capture_method: str = "CROP-seq", seed: int = 0, engine: str = "python",
                  upstream_bin: Optional[str] = None, workdir: Optional[Path] = None) -> "tuple[MData, dict]":
    """prepare_assignment.py splits by batch; each batch is assigned independently, then re-assembled."""
    np, sp, pd = _np(), _sp(), _pd()
    g = m["guide"]
    C = as_csr(g.X).tocsc()
    n, k = C.shape
    obs = m.obs if batch_key in m.obs.columns else g.obs
    batches = obs[batch_key].astype(str).values if batch_key in obs.columns else np.array(["all"] * n)
    tot_guide = row_sums(g.X)
    gene_tot = m["gene"].obs["total_gene_umis"].values if "total_gene_umis" in m["gene"].obs.columns else row_sums(m["gene"].X)
    A = sp.lil_matrix((n, k), dtype=np.float64)
    used = set()
    for b in sorted(set(batches)):
        rows = np.where(batches == b)[0]
        if engine == "external" and method in ("cleanser", "sceptre"):
            ext = _external_assignment(m, rows, method, threshold, n_em_rep, capture_method, upstream_bin, (workdir or Path(".")) / f"assign_{safe_label(b)}")
            if ext is not None:
                ext = ext.tocoo()
                for i, j, v in zip(ext.row, ext.col, ext.data):
                    if v > 0:
                        A[rows[i], j] = 1.0
                used.add(f"external-{method}")
                continue
            print(f"external {method} unavailable; using the Python {method} assignment")
        used.add({"sceptre": "python-sceptre-mixture", "cleanser": "python-cleanser", "threshold": "threshold"}.get(method, method))
        for j in range(k):
            y = np.asarray(C[rows, j].toarray()).ravel()
            if y.sum() == 0:
                continue
            if method == "sceptre":
                thr = 0.8 if threshold is None else threshold
                call, _ = sceptre_mixture_assign(y, [np.log(np.maximum(gene_tot[rows], 1)), np.log1p(tot_guide[rows])], thr, n_em_rep, seed + j)
            elif method == "cleanser":
                call, _ = cleanser_fallback_assign(y, tot_guide[rows], 1.0 if threshold is None else threshold)
            elif method == "threshold":
                call = y >= umi_threshold
            else:
                raise SystemExit(f"Invalid GUIDE_ASSIGNMENT_method: {method}")
            for i in rows[call]:
                A[i, j] = 1.0
    g = g.copy()
    g.layers["guide_assignment"] = A.tocsr()
    out = MData({**m.mod, "guide": g}, m.obs.copy(), dict(m.uns))
    npc = nnz_rows(g.layers["guide_assignment"])
    return out, {"method": method, "engine": sorted(used), "n_batches": len(set(batches)), "cells_with_guide": int((npc > 0).sum()),
                 "mean_guides_per_cell": float(npc.mean()), "n_assignments": int(npc.sum())}


def filter_genes_by_cells(m: MData, fraction: float) -> MData:
    """mudata_concat.py: keep genes expressed in MORE than int(n_cells * fraction) cells."""
    np = _np()
    g = m["gene"]
    keep = np.asarray((g.X > 0).sum(axis=0)).ravel() > int(g.n_obs * fraction)
    return MData({**m.mod, "gene": g[:, keep].copy()}, m.obs, m.uns)


def collapse_guides(m: MData, max_elements: int = 1, max_guides: int = 2, min_guides: int = 2) -> "tuple[MData, dict]":
    """collapse_guides.py (DUAL_GUIDE): cells with exactly 2 guides on one element -> one 'g1|g2' pseudo-guide."""
    np, pd, sp = _np(), _pd(), _sp()
    import anndata as ad  # type: ignore
    g = m["guide"]
    A = as_csr(g.layers["guide_assignment"])
    gpc = row_sums(A)
    gbe = pd.get_dummies(g.var["intended_target_name"].astype(str)).astype(float)
    cbe = (A @ gbe.values) > 0
    epc = cbe.sum(axis=1)
    keep = (epc <= max_elements) & (gpc >= min_guides) & (gpc <= max_guides)
    sub = g[keep, :]
    Asub = as_csr(sub.layers["guide_assignment"])
    ci, gi = Asub.nonzero()
    bc = pd.DataFrame({"cell_bc": sub.obs_names[ci]})
    vdf = sub.var.iloc[gi].reset_index(drop=True)
    bv = pd.concat([bc, vdf], axis=1)
    grp = bv.sort_values("guide_id").groupby(["cell_bc", "intended_target_name"], observed=True)
    agg = {"guide_id": lambda x: "|".join(map(str, x))}
    if "spacer" in bv.columns:
        agg["spacer"] = lambda x: "|".join(map(str, x))
    joined = grp.agg(agg)
    rest = [c for c in ["intended_target_chr", "intended_target_start", "intended_target_end", "intended_target_key", "targeting", "type", "pam",
                        "gene_name"] if c in bv.columns]
    first = grp[rest].agg("first")
    ded = pd.concat([joined, first], axis=1).reset_index().set_index("cell_bc")
    ded["guide_chr"], ded["guide_start"], ded["guide_end"] = ded["intended_target_chr"], ded["intended_target_start"], ded["intended_target_end"]
    dm = pd.get_dummies(ded["guide_id"], dtype=float)
    var_new = ded.groupby("guide_id").first()
    var_new["guide_id"] = var_new.index
    cells = list(dm.index)
    na = ad.AnnData(sp.csr_matrix(dm.values), obs=g.obs.loc[cells].copy(), var=var_new.loc[dm.columns].copy())
    na.layers["guide_assignment"] = na.X.copy()
    na.uns = dict(g.uns)
    mods = {k: (na if k == "guide" else v[cells, :].copy()) for k, v in m.mod.items()}
    out = MData(mods, m.obs.loc[cells].copy() if set(cells) <= set(m.obs.index) else None, dict(m.uns))
    return out, {"n_cells_collapsed": len(cells), "n_assignments": int(len(bv)), "n_combinations": int(dm.shape[1])}


def cmd_assign_guides(args) -> int:
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    m, info = assign_guides(m, args.method, args.threshold, args.n_em_rep, args.umi_threshold, args.batch_key, args.capture_method,
                            engine=args.engine, upstream_bin=args.upstream_bin, workdir=out / "work")
    n0 = m["gene"].n_vars
    m = filter_genes_by_cells(m, args.min_cells_fraction)
    info.update(genes_before_filter=int(n0), genes_after_filter=int(m["gene"].n_vars))
    if args.dual_guide:
        m, ci = collapse_guides(m)
        info["collapse"] = ci
    write_h5mu(m, out / "concat_mudata.h5mu")
    finish(out, "Guide assignment", info, [md_table(["key", "value"], [[k, v] for k, v in info.items()])])
    return 0


def nt_intended_targets(guide_var, strategy: str = "same"):
    """create_nt_intended_targets.py: controls -> one 'non-targeting' element (same) or median-size groups (median)."""
    np = _np()
    gv = guide_var.copy()
    ctrl = (~gv["targeting"].map(_truthy)) | gv["type"].isin(list(CONTROL_TYPES))
    per = gv[~ctrl].groupby("intended_target_name", observed=True).size()
    if strategy == "median":
        med = int(per.median()) if len(per) else 1
        n = int(ctrl.sum())
        names = [f"non-targeting|{i + 1}" for i in range(int(np.ceil(n / max(med, 1)))) for _ in range(med)][:n]
    elif strategy == "same":
        names = ["non-targeting"] * int(ctrl.sum())
    else:
        raise SystemExit("strategy must be 'median' or 'same'")
    gv["intended_target_name"] = gv["intended_target_name"].astype(str)
    gv.loc[ctrl, "intended_target_name"] = names
    gv["intended_target_name"] = gv["intended_target_name"].astype("category")
    return gv


def cmd_nt_targets(args) -> int:
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    g = m["guide"].copy()
    g.var = nt_intended_targets(g.var, args.strategy)
    m = MData({**m.mod, "guide": g}, m.obs, m.uns)
    write_h5mu(m, out / "mudata_nt_targets.h5mu")
    vc = g.var["intended_target_name"].astype(str).value_counts()
    finish(out, "Control-guide intended targets", {"strategy": args.strategy, "n_elements": int(len(vc))},
           [md_table(["intended_target_name", "guides"], [[k, v] for k, v in vc.head(30).items()])])
    return 0


# ---------------------------------------------------------------------------
# pairs_to_test: create_pairs_to_test.py + prepare_inference.py
# ---------------------------------------------------------------------------

def create_pairs_to_test(guide_var, gene_ids, gtf, limit: Optional[int] = 1_000_000):
    """Per guide, GTF genes on guide_chr with |gene start - guide_start| <= limit (limit None: every gene in the MuData)."""
    pd, np = _pd(), _np()
    gd = annotate_intended_target_groups(guide_var)
    cols = ["guide_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end", "intended_target_key"]
    out = []
    if limit is None:
        genes = list(map(str, gene_ids))
        for _, r in gd.iterrows():
            d = pd.DataFrame({"guide_id": r["guide_id"], "gene_id": genes})
            for c in cols[1:]:
                d[c] = r[c]
            out.append(d)
    else:
        g = gtf.drop_duplicates("gene_name")
        by_chr = {c: sub for c, sub in g.groupby("seqname")}
        for _, r in gd.iterrows():
            gc, gs = r.get("guide_chr"), r.get("guide_start")
            if not isinstance(gc, str) or gc not in by_chr or pd.isna(gs):
                continue
            sub = by_chr[gc]
            sub = sub[np.abs(sub["start"].values - float(gs)) <= limit]
            if not len(sub):
                continue
            d = pd.DataFrame({"guide_id": r["guide_id"], "gene_id": sub["gene_id"].values})
            for c in cols[1:]:
                d[c] = r[c]
            out.append(d)
    if not out:
        return pd.DataFrame(columns=["guide_id", "gene_name"] + cols[1:] + ["pair_type"])
    df = pd.concat(out, ignore_index=True)
    df["pair_type"] = "discovery"
    df["gene_id"] = df["gene_id"].astype(str).str.split(".").str[0]
    df = df.rename(columns={"gene_id": "gene_name"})
    return df[["guide_id", "gene_name", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end",
               "intended_target_key", "pair_type"]]


def prepare_inference(m: MData, pairs, subset_for_cis: bool = False) -> "tuple[MData, dict]":
    np, pd = _np(), _pd()
    p = pairs.copy()
    if "gene_name" not in p.columns and "gene_id" in p.columns:
        p["gene_name"] = p["gene_id"]
    missing = [c for c in ("guide_id", "gene_name") if c not in p.columns]
    if missing:
        raise SystemExit(f"pairs file is missing columns: {missing}")
    if "pair_type" not in p.columns:
        p["pair_type"] = "discovery"
    p["guide_id"], p["gene_name"] = p["guide_id"].astype(str), p["gene_name"].astype(str)
    g = m["guide"].copy()
    g.var = annotate_intended_target_groups(g.var)
    genes, guides = set(map(str, m["gene"].var_names)), set(g.var["guide_id"].astype(str))
    sub = p[p["gene_name"].isin(genes) & p["guide_id"].isin(guides)].copy()
    if sub.empty:
        raise SystemExit("The subset of pairs_to_test is empty after filtering. Please check your input data.")
    sub = enrich_pairs_with_target_metadata(sub, g.var)
    sub = sub.rename(columns={"gene_name": "gene_id"})
    mods = {**m.mod, "guide": g}
    info = {"n_pairs_in": int(len(p)), "n_pairs": int(len(sub)), "n_genes_tested": int(sub["gene_id"].nunique()),
            "n_guides_tested": int(sub["guide_id"].nunique())}
    uns = dict(m.uns)
    uns["pairs_to_test"] = sub
    if subset_for_cis:
        gmask = m["gene"].var_names.isin(sub["gene_id"].unique())
        names = g.var["intended_target_name"].astype(str)
        ctrl = names.str.startswith("non-targeting|") | names.eq("non-targeting")
        gdm = g.var["guide_id"].astype(str).isin(sub["guide_id"].unique()).values | ctrl.values
        rs = m["gene"][:, gmask]
        gs = g[:, gdm]
        cells = row_sums(gs.X) > 0
        if not cells.any():
            raise SystemExit("No targeted cells found in the guide modality subset.")
        mods = {"gene": rs[cells].copy(), "guide": gs[cells].copy()}
        for k, v in m.mod.items():
            if k not in mods:
                mods[k] = v[cells].copy()
        obs = m.obs.iloc[np.where(cells)[0]].copy() if len(m.obs) == len(cells) else None
        out = MData(mods, obs, uns)
        info.update(cis_cells=int(cells.sum()), cis_genes=int(gmask.sum()), cis_guides=int(gdm.sum()))
        return out, info
    return MData(mods, m.obs.copy(), uns), info


def cmd_pairs(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    strategy = args.strategy
    if strategy == "predefined_pairs":
        if not args.pairs:
            raise SystemExit("--pairs is required for predefined_pairs")
        pairs = pd.read_csv(args.pairs)
    elif strategy in ("default", "by_distance", "all"):
        gtf = read_gtf(Path(args.gtf)) if args.gtf else None
        if gtf is None and strategy != "all":
            raise SystemExit("--gtf is required for distance-based pairs")
        pairs = create_pairs_to_test(m["guide"].var, m["gene"].var_names, gtf, None if strategy == "all" or args.limit == -1 else args.limit)
    else:
        raise SystemExit(f"Invalid INFERENCE_target_guide_pairing_strategy: {strategy}")
    write_tsv(pairs, out / "pairs_to_test.csv", sep=",")
    m2, info = prepare_inference(m, pairs, subset_for_cis=(strategy == "default"))
    write_h5mu(m2, out / "mudata_inference_input.h5mu")
    info["strategy"] = strategy
    finish(out, "Pairs to test", info, [md_table(["key", "value"], [[k, v] for k, v in info.items()])])
    return 0


# ---------------------------------------------------------------------------
# chunk_mudata.py (PerTurbo, balanced) and chunk_mudata_sceptre.py (auto / off / force)
# ---------------------------------------------------------------------------

def plan_chunk_ranges(n_genes: int, chunk_size: int) -> "list[tuple[int, int]]":
    if n_genes <= 0:
        return []
    if chunk_size <= 0 or n_genes <= chunk_size:
        return [(0, n_genes)]
    nc = math.ceil(n_genes / chunk_size)
    base, rem = divmod(n_genes, nc)
    out, s = [], 0
    for i in range(nc):
        e = s + base + (1 if i < rem else 0)
        out.append((s, e))
        s = e
    return out


def sceptre_chunk_plan(n_cells: int, n_genes: int, chunk_size: int, mode: str = "auto", threshold: int = 2147483647) -> "tuple[bool, list]":
    entries = int(n_cells) * int(n_genes)
    should = True if mode == "force" else (False if mode == "off" else entries > threshold)
    if not should or n_genes <= chunk_size:
        return False, [(0, n_genes)]
    return True, [(i * chunk_size, min((i + 1) * chunk_size, n_genes)) for i in range(math.ceil(n_genes / chunk_size))]


def cmd_chunk(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    ng, nc = int(m["gene"].n_vars), m.n_obs
    pairs = uns_df(m, "pairs_to_test")
    if args.mode == "sceptre":
        chunked, ranges = sceptre_chunk_plan(nc, ng, args.chunk_size, args.chunk_mode, args.auto_threshold_entries)
    else:
        chunked, ranges = len(plan_chunk_ranges(ng, args.chunk_size)) > 1, plan_chunk_ranges(ng, args.chunk_size)
    rows = []
    for i, (s, e) in enumerate(ranges):
        gene = m["gene"][:, s:e].copy()
        guide = m["guide"]
        uns = dict(m.uns)
        if pairs is not None:
            cp = pairs[pairs["gene_id"].isin(set(gene.var_names))]
            uns["pairs_to_test"] = cp
            if not args.keep_all_guides and args.mode == "perturbo":
                guide = guide[:, guide.var["guide_id"].isin(set(cp["guide_id"]))]
        f = out / "chunks" / (f"chunk.{i:03d}.h5mu" if args.mode == "sceptre" else f"chunk.{i:02d}.h5mu")
        write_h5mu(MData({"gene": gene, "guide": guide.copy()}, m.obs.copy(), uns), f, quiet=True)
        rows.append({"chunk_id": i, "chunk_file": f.name, "gene_start": s, "gene_end": e - 1, "chunk_gene_count": e - s,
                     "chunk_mode": "chunked" if chunked else "single", "mode": args.chunk_mode if args.mode == "sceptre" else "balanced",
                     "chunked": chunked, "n_cells": nc, "n_genes": ng, "entries": nc * ng, "threshold": args.auto_threshold_entries,
                     "chunk_size": args.chunk_size, "keep_all_guides": bool(args.keep_all_guides)})
    man = pd.DataFrame(rows)
    write_tsv(man, out / "chunk_manifest.tsv")
    finish(out, "MuData chunks", {"n_chunks": len(rows), "chunked": chunked}, [md_table(list(man.columns[:6]), man.iloc[:, :6].values.tolist())])
    return 0


# ---------------------------------------------------------------------------
# Inference engines (inference_sceptre.R / perturbo_inference.py; Python ports of the methods)
# ---------------------------------------------------------------------------

def nb_theta_ml(y, mu) -> float:
    np = _np()
    from scipy.optimize import minimize_scalar  # type: ignore
    from scipy.special import gammaln  # type: ignore
    mu = np.maximum(mu, 1e-10)

    def nll(lt):
        th = math.exp(lt)
        return -float(np.sum(gammaln(y + th) - gammaln(th) + th * np.log(th / (th + mu)) + y * np.log(mu / (th + mu))))
    if np.var(y) <= np.mean(y) * 1.0001:
        return 1e4
    r = minimize_scalar(nll, bounds=(-4.0, 9.0), method="bounded", options={"xatol": 1e-3})
    return float(math.exp(r.x))


def fit_nb_null(y, X, offset=None, n_iter: int = 3):
    """Covariate-only NB GLM (log link): Poisson IRLS start, ML theta, then NB IRLS refinements."""
    np = _np()
    off = np.zeros(len(y)) if offset is None else offset
    beta, mu = poisson_irls(X, y, off, n_iter=15)
    theta = nb_theta_ml(y, mu)
    for _ in range(n_iter):
        eta = X @ beta + off
        mu = np.exp(np.clip(eta, -30, 30))
        W = mu / (1 + mu / theta)
        z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        XtW = X.T * W
        try:
            beta = np.linalg.solve(XtW @ X + 1e-8 * np.eye(X.shape[1]), XtW @ z)
        except np.linalg.LinAlgError:
            break
    mu = np.exp(np.clip(X @ beta + off, -30, 30))
    return mu, nb_theta_ml(y, mu)


def skewnorm_moment_pvalue(z: float, null, side: str = "both") -> float:
    """SCEPTRE's resampling approximation: skew-normal fitted by moments to the resampled statistics."""
    np = _np()
    from scipy.stats import skew, skewnorm  # type: ignore
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if len(null) < 10 or not np.isfinite(z):
        return float("nan")
    m, s = float(null.mean()), float(null.std(ddof=1))
    if s <= 0:
        return 1.0
    g = float(np.clip(skew(null), -0.99, 0.99))
    a23 = abs(g) ** (2 / 3)
    delta = math.copysign(math.sqrt(math.pi / 2 * a23 / (a23 + ((4 - math.pi) / 2) ** (2 / 3))), g) if g != 0 else 0.0
    alpha = delta / math.sqrt(1 - delta ** 2)
    omega = s / math.sqrt(1 - 2 * delta ** 2 / math.pi)
    xi = m - omega * delta * math.sqrt(2 / math.pi)
    lo, hi = float(skewnorm.cdf(z, alpha, xi, omega)), float(skewnorm.sf(z, alpha, xi, omega))
    if side == "left":
        return lo
    if side == "right":
        return hi
    return float(min(1.0, 2 * min(lo, hi)))


def nb_mle_fold_change(y, mu0, theta: float, n_iter: int = 30) -> "tuple[float, float]":
    """MLE of beta in mu = mu0 * exp(beta) over the treated cells (NB, theta fixed); returns (beta, SE)."""
    np = _np()
    b = 0.0
    for _ in range(n_iter):
        mu = mu0 * math.exp(b)
        U = float(np.sum((y - mu) * theta / (theta + mu)))
        I = float(np.sum(mu * theta / (theta + mu)))
        if I <= 0:
            break
        step = U / I
        b = float(np.clip(b + step, -15, 15))
        if abs(step) < 1e-8:
            break
    mu = mu0 * math.exp(b)
    I = float(np.sum(mu * theta / (theta + mu)))
    return b, (1 / math.sqrt(I) if I > 0 else float("inf"))


class InferenceData:
    """Everything the tests need from a MuData: counts, covariates, the assignment matrix and the guide/element units."""

    def __init__(self, m: MData, batch_key: str = "batch"):
        np, pd, sp = _np(), _pd(), _sp()
        self.m = m
        g = m["gene"]
        self.Y = as_csr(g.X).tocsc()
        self.n = g.n_obs
        self.gene_ids = [str(x) for x in g.var_names]
        self.gidx = {x: i for i, x in enumerate(self.gene_ids)}
        lib = g.obs["total_gene_umis"].values.astype(float) if "total_gene_umis" in g.obs.columns else row_sums(g.X)
        lib = np.maximum(lib, 1.0)
        nnz = g.obs["num_expressed_genes"].values.astype(float) if "num_expressed_genes" in g.obs.columns else nnz_rows(g.X).astype(float)
        obs = m.obs if batch_key in m.obs.columns else g.obs
        batch = obs[batch_key].astype(str).values if batch_key in obs.columns else np.array(["b"] * self.n)
        bd = pd.get_dummies(pd.Series(batch)).values.astype(float)[:, 1:] if len(set(batch)) > 1 else np.zeros((self.n, 0))
        gu = m["guide"].obs["total_guide_umis"].values.astype(float) if "total_guide_umis" in m["guide"].obs.columns else row_sums(m["guide"].X)
        lgu = np.log1p(gu)
        self.X_sceptre = np.column_stack([np.ones(self.n), np.log(lib), np.log(np.maximum(nnz, 1)), bd])
        self.X_perturbo = np.column_stack([np.ones(self.n), bd, lgu - lgu.mean()])
        self.off_perturbo = np.log(lib)
        gv = annotate_intended_target_groups(m["guide"].var)
        self.guide_var = gv
        A = m["guide"].layers["guide_assignment"] if "guide_assignment" in m["guide"].layers else (m["guide"].X > 0)
        A = as_csr(A).tocsc()
        self.guide_ids = [str(x) for x in gv["guide_id"]]
        self.treated_guide = {gid: np.asarray(A[:, j].toarray()).ravel() > 0 for j, gid in enumerate(self.guide_ids)}
        self.treated_elem = {}
        for key, sub in gv.groupby("intended_target_key"):
            mask = np.zeros(self.n, bool)
            for gid in sub["guide_id"].astype(str):
                mask |= self.treated_guide[gid]
            self.treated_elem[str(key)] = mask
        self.lookup = target_lookup(gv).set_index("intended_target_key")
        self._null: "dict[tuple, tuple]" = {}

    def null(self, gene: str, slot: str):
        k = (gene, slot)
        if k not in self._null:
            np = _np()
            y = np.asarray(self.Y[:, self.gidx[gene]].toarray()).ravel()
            if slot == "sceptre":
                mu, th = fit_nb_null(y, self.X_sceptre)
            else:
                mu, th = fit_nb_null(y, self.X_perturbo, self.off_perturbo)
            self._null = {kk: v for kk, v in self._null.items() if kk[0] == gene}  # keep one gene cached
            self._null[k] = (y, mu, th)
        return self._null[k]


def _pair_list(d: InferenceData, unit: str, pairs=None, test_all: bool = False) -> "list[tuple[str, str]]":
    units = d.guide_ids if unit == "guide" else list(d.treated_elem)
    if test_all or pairs is None or not len(pairs):
        return [(u, g) for g in d.gene_ids for u in units]
    col = "guide_id" if unit == "guide" else "intended_target_key"
    p = pairs
    if col not in p.columns:
        p = enrich_pairs_with_target_metadata(p.rename(columns={"gene_id": "gene_name"}), d.guide_var).rename(columns={"gene_name": "gene_id"})
    p = p[p["gene_id"].astype(str).isin(d.gidx) & p[col].astype(str).isin(set(units))]
    return list(dict.fromkeys(zip(p[col].astype(str), p["gene_id"].astype(str))))


def run_tests(d: InferenceData, unit: str, slot: str, pairs=None, test_all: bool = False, n_resamples: int = 500, side: str = "both",
              seed: int = 0):
    """One unit type; see run_tests_multi."""
    return run_tests_multi(d, (unit,), slot, pairs, test_all, n_resamples, side, seed)[unit]


def run_tests_multi(d: InferenceData, units=("guide", "element"), slot: str = "sceptre", pairs=None, test_all: bool = False,
                    n_resamples: int = 500, side: str = "both", seed: int = 0, perm_budget: int = 50_000_000) -> dict:
    """slot 'sceptre': NB score statistic vs label permutations (complement control), skew-normal p-value;
       slot 'perturbo': NB MLE fold change on the covariate-adjusted null, Wald p-value, log2_fc_std.
    One covariate-only null fit per gene serves every unit type. Permutation index sets are drawn from a per-unit
    seed (crc32 of the unit id), so results do not depend on test order or on evictions from the bounded cache."""
    import zlib
    np, pd = _np(), _pd()
    from scipy.stats import norm  # type: ignore
    by_gene: "dict[str, list]" = defaultdict(list)
    for ut in units:
        for u, g in _pair_list(d, ut, pairs, test_all):
            by_gene[g].append((ut, u))
    perm_cache: "dict[str, Any]" = {}
    cache_entries = 0
    rows: "dict[str, list]" = {ut: [] for ut in units}
    for gene in d.gene_ids:
        if gene not in by_gene:
            continue
        y, mu, th = d.null(gene, slot)
        if slot == "sceptre":
            r = (y - mu) / (1 + mu / th)
            w = mu / (1 + mu / th)
        for ut, u in by_gene[gene]:
            tr = d.treated_guide if ut == "guide" else d.treated_elem
            T = np.where(tr[u])[0]
            if len(T) == 0:
                rows[ut].append((gene, u, float("nan"), float("nan"), float("nan")))
                continue
            if slot == "sceptre":
                ck = f"{ut}:{u}"
                if ck not in perm_cache:
                    rng = np.random.default_rng(seed + zlib.crc32(ck.encode()))
                    P = np.stack([rng.choice(d.n, size=len(T), replace=False) for _ in range(n_resamples)]).astype(np.int32)
                    while perm_cache and cache_entries + P.size > perm_budget:
                        old = next(iter(perm_cache))
                        cache_entries -= perm_cache.pop(old).size
                    perm_cache[ck] = P
                    cache_entries += P.size
                P = perm_cache[ck]
                z = float(r[T].sum() / math.sqrt(max(w[T].sum(), 1e-12)))
                zn = r[P].sum(axis=1) / np.sqrt(np.maximum(w[P].sum(axis=1), 1e-12))
                p = skewnorm_moment_pvalue(z, zn, side)
                lfc = math.log2((y[T].sum() + 0.5) / (mu[T].sum() + 0.5))
                rows[ut].append((gene, u, lfc, p, float("nan")))
            else:
                b, se = nb_mle_fold_change(y[T], mu[T], th)
                p = float(2 * norm.sf(abs(b / se))) if math.isfinite(se) and se > 0 else float("nan")
                rows[ut].append((gene, u, b / math.log(2), p, se / math.log(2)))
    out = {}
    for ut in units:
        ucol = "guide_id" if ut == "guide" else "intended_target_key"
        df = pd.DataFrame(rows[ut], columns=["gene_id", ucol, "log2_fc", "p_value", "log2_fc_std"])
        if ut == "element":
            df = df.merge(d.lookup.reset_index(), on="intended_target_key", how="left")
            cols = ["gene_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end", "log2_fc"]
        else:
            cols = ["gene_id", "guide_id", "log2_fc"]
        cols += (["log2_fc_std", "p_value"] if slot == "perturbo" else ["p_value"])
        sort = ["gene_id", "guide_id"] if ut == "guide" else ["gene_id", "intended_target_name"]
        out[ut] = df[cols].sort_values(sort).reset_index(drop=True)
    return out


def external_inference(slot: str, unit: str, mudata_path: Path, workdir: Path, test_all: bool, upstream_bin: Optional[str],
                       sceptre_args: dict) -> "Optional[Any]":
    """Run the upstream SCEPTRE (Rscript) / PerTurbo (python) scripts when available; else print the commands."""
    pd = _pd()
    workdir.mkdir(parents=True, exist_ok=True)
    ub = Path(upstream_bin) if upstream_bin else None
    if slot == "sceptre":
        (workdir / "args.txt").write_text("\n".join([str(mudata_path), sceptre_args.get("side", "both"), "union",
                                                      sceptre_args.get("resampling_approximation", "skew_normal"), "complement", "default", "4"]) + "\n")
        script = ub / "inference_sceptre.R" if ub else Path("inference_sceptre.R")
        ok = ub is not None and script.is_file()
        ran, _ = external(["Rscript", str(script), "args.txt"], "Rscript", workdir, run=ok)
        if not ran:
            return None
        f = workdir / ("per_guide_output.tsv.gz" if unit == "guide" else "per_element_output.tsv.gz")
        return pd.read_csv(f, sep="\t")
    out = workdir / f"perturbo_{unit}.tsv.gz"
    script = ub / "perturbo_inference.py" if ub else Path("perturbo_inference.py")
    cmd = [sys.executable, str(script), str(mudata_path), str(out), "--batch_size", "4096", "--num_workers", "0", "--efficiency_mode", "scaled",
           "--inference_type", unit] + (["--test_all_pairs"] if test_all else [])
    ok = ub is not None and script.is_file()
    try:
        import perturbo  # type: ignore  # noqa: F401
    except ImportError:
        ok = False
    ran, _ = external(cmd, sys.executable, workdir, run=ok)
    return pd.read_csv(out, sep="\t") if ran else None


def merge_method_results(sg, se, pg, pe):
    """merge_method_results.py: outer merges; element keys include coordinates when both tables have them."""
    pd = _pd()
    sg = sg.rename(columns={"log2_fc": "sceptre_log2_fc", "p_value": "sceptre_p_value"})
    se = se.rename(columns={"log2_fc": "sceptre_log2_fc", "p_value": "sceptre_p_value"})
    pg = pg.rename(columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"})
    pe = pe.rename(columns={"log2_fc": "perturbo_log2_fc", "p_value": "perturbo_p_value"})
    guide = pd.merge(sg[["gene_id", "guide_id", "sceptre_log2_fc", "sceptre_p_value"]], pg[["gene_id", "guide_id", "perturbo_log2_fc", "perturbo_p_value"]],
                     on=["gene_id", "guide_id"], how="outer")
    full = ["gene_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]
    keys = full if all(c in se.columns for c in full) and all(c in pe.columns for c in full) else ["gene_id", "intended_target_name"]

    def keyed(df, vals):
        d = df[keys + vals].copy()
        mk = []
        for c in keys:
            d[f"__merge_{c}"] = d[c].map(lambda v: "__NA__" if pd.isna(v) else (str(int(v)) if isinstance(v, float) and float(v).is_integer() else str(v)))
            mk.append(f"__merge_{c}")
        return d, mk
    a, lk = keyed(se, ["sceptre_log2_fc", "sceptre_p_value"])
    b, rk = keyed(pe, ["perturbo_log2_fc", "perturbo_p_value"])
    el = pd.merge(a, b, left_on=lk, right_on=rk, how="outer", suffixes=("", "_perturbo"))
    for c in keys:
        if f"{c}_perturbo" in el.columns:
            el[c] = el[c].combine_first(el[f"{c}_perturbo"])
            el = el.drop(columns=[f"{c}_perturbo"])
    el = el.drop(columns=[c for c in el.columns if c.startswith("__merge_")])
    el = el[keys + ["sceptre_log2_fc", "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value"]]
    return guide, el


def export_outputs(m: MData, per_guide, method_cols: bool = True) -> "tuple[Any, Any]":
    """export_output_single.py: per-guide + per-element tables with avg target expression and cell number."""
    np, pd = _np(), _pd()
    gv = m["guide"].var.reset_index(drop=True)
    res = per_guide.copy()
    if "intended_target_name" not in res.columns:
        res = res.merge(gv[["guide_id", "intended_target_name"]].astype({"guide_id": str}), on="guide_id", how="left")
    res = res.merge(gv[[c for c in ["guide_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"] if c in gv.columns]],
                    on=["intended_target_name", "guide_id"], how="left")
    for c in ("log2_fc", "p_value"):
        if c in res.columns and f"sceptre_{c}" not in res.columns and f"perturbo_{c}" not in res.columns:
            res = res.rename(columns={c: f"perturbo_{c}"})
    for c in ("sceptre_log2_fc", "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value"):
        if c not in res.columns:
            res[c] = np.nan
    X = as_csr(m["gene"].X)
    sym = _gene_symbols(m["gene"])
    pos: "dict[str, list]" = defaultdict(list)
    for i, s in enumerate(sym):
        pos[str(s)].append(i)
    avg = {t: (float(X[:, pos[t]].mean()) if t in pos else np.nan) for t in res["intended_target_name"].dropna().unique()}
    res["avg_gene_expression"] = res["intended_target_name"].map(avg)
    res["cell_number"] = m.n_obs
    cols = ["intended_target_name", "guide_id", "intended_target_chr", "intended_target_start", "intended_target_end", "gene_id",
            "sceptre_log2_fc", "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value", "cell_number", "avg_gene_expression"]
    pg = res.reindex(columns=cols).rename(columns={"guide_id": "guide_id(s)"})
    pe = pg.groupby(["gene_id", "intended_target_name"], as_index=False, observed=True).agg({
        "guide_id(s)": lambda x: ",".join(x.dropna().astype(str)), "intended_target_chr": "first", "intended_target_start": "first",
        "intended_target_end": "first", "sceptre_log2_fc": "first", "sceptre_p_value": "first", "perturbo_log2_fc": "first",
        "perturbo_p_value": "first", "cell_number": "first", "avg_gene_expression": "first"})
    return pg, pe[[c if c != "guide_id" else "guide_id(s)" for c in cols]]


def run_inference(m_cis: MData, method: str = "default", m_trans: Optional[MData] = None, engine: str = "python", n_resamples: int = 500,
                  side: str = "both", workdir: Optional[Path] = None, upstream_bin: Optional[str] = None, seed: int = 0) -> "tuple[MData, dict, dict]":
    """Returns (inference MuData, result tables by uns key, summary)."""
    np, pd = _np(), _pd()
    workdir = workdir or Path(".")
    pairs = uns_df(m_cis, "pairs_to_test")
    tables: "dict[str, Any]" = {}
    t0 = time.time()

    def slot_tables(m: MData, slot: str, test_all: bool, tag: str):
        if engine == "external":
            p = write_h5mu(m, workdir / f"{tag}_input.h5mu", quiet=True)
            g = external_inference(slot, "guide", p, workdir / f"{tag}_{slot}", test_all, upstream_bin, {"side": side})
            e = external_inference(slot, "element", p, workdir / f"{tag}_{slot}", test_all, upstream_bin, {"side": side})
            if g is not None and e is not None:
                return g, e, f"external-{slot}"
            print(f"external {slot} unavailable; falling back to the Python {slot} engine")
        d = InferenceData(m)
        pp = None if test_all else uns_df(m, "pairs_to_test")
        res = run_tests_multi(d, ("guide", "element"), slot, pp, test_all, n_resamples, side, seed)
        return res["guide"], res["element"], f"python-{slot}"

    engines = []
    if method in ("sceptre", "perturbo"):
        g, e, en = slot_tables(m_cis, method, pairs is None, "single")
        engines.append(en)
        tables = {"per_guide_results": g, "per_element_results": e}
        base = m_cis
    elif method == "sceptre,perturbo":
        sg, se, e1 = slot_tables(m_cis, "sceptre", pairs is None, "single")
        pg, pe, e2 = slot_tables(m_cis, "perturbo", pairs is None, "single")
        engines += [e1, e2]
        mg, me = merge_method_results(sg, se, pg, pe)
        tables = {"per_guide_results": mg, "per_element_results": me}
        base = m_cis
    elif method == "default":
        if pairs is None:
            raise SystemExit("INFERENCE_method='default' requires pairs_to_test (strategy 'default')")
        sg, se, e1 = slot_tables(m_cis, "sceptre", False, "cis")
        pg, pe, e2 = slot_tables(m_cis, "perturbo", False, "cis")
        mg, me = merge_method_results(sg, se, pg, pe)
        mt = m_trans if m_trans is not None else m_cis
        tg, te, e3 = slot_tables(mt, "perturbo", True, "trans")
        engines += [e1, e2, e3]
        tables = {"cis_per_guide_results": mg, "cis_per_element_results": me, "trans_per_guide_results": tg, "trans_per_element_results": te}
        base = mt
    else:
        raise SystemExit(f"Invalid INFERENCE_method: {method}. Valid options: sceptre, perturbo, sceptre,perturbo, default")
    uns = {k: v for k, v in base.uns.items() if k not in ("per_guide_results", "per_element_results")}
    if "pairs_to_test" not in uns and pairs is not None:
        uns["pairs_to_test"] = pairs
    uns.update(tables)
    uns["inference_engine"] = np.array(sorted(set(engines)), dtype=object)
    uns["inference_method"] = np.array([method], dtype=object)
    out = MData(dict(base.mod), base.obs.copy(), uns)
    summ = {"method": method, "engines": sorted(set(engines)), "seconds": round(time.time() - t0, 1),
            "tables": {k: int(len(v)) for k, v in tables.items()}, "n_resamples": n_resamples}
    return out, tables, summ


def cmd_inference(args) -> int:
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    mt = read_h5mu(Path(args.trans_mudata)) if args.trans_mudata else None
    im, tables, summ = run_inference(m, args.method, mt, args.engine, args.n_resamples, args.side, out / "work", args.upstream_bin, args.seed)
    write_h5mu(im, out / "inference_mudata.h5mu")
    for k, v in tables.items():
        write_tsv(v, out / f"{k}.tsv.gz")
    pgk = "cis_per_guide_results" if "cis_per_guide_results" in tables else "per_guide_results"
    pg, pe = export_outputs(im, tables[pgk])
    write_tsv(pg, out / "per_guide_output.tsv")
    write_tsv(pe, out / "per_element_output.tsv")
    secs = [md_table(["table", "rows"], [[k, len(v)] for k, v in tables.items()])]
    for k, v in tables.items():
        pc = next((c for c in ("sceptre_p_value", "perturbo_p_value", "p_value") if c in v.columns), None)
        if pc:
            top = v.sort_values(pc).head(10)
            secs.append(f"### {k} (top 10 by {pc})\n\n" + md_table(list(top.columns), top.values.tolist()))
    finish(out, f"Inference ({args.method})", summ, secs)
    return 0


def merge_chunk_tables(files: "list[str]", keys: "list[str]"):
    pd = _pd()
    df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    dup = df.duplicated(keys, keep=False)
    if dup.any():
        print(f"Warning: duplicate {keys} keys in chunk results ({int(dup.sum())} rows)")
    return df.sort_values(keys).reset_index(drop=True)


def cmd_merge_results(args) -> int:
    pd = _pd()
    out = run_dir(args.label)
    base = read_h5mu(Path(args.base_mudata)) if args.base_mudata else None
    uns_add: "dict[str, Any]" = {}
    rd = lambda p: pd.read_csv(p, sep="\t")  # noqa: E731
    if args.mode == "methods":
        g, e = merge_method_results(rd(args.sceptre_per_guide), rd(args.sceptre_per_element), rd(args.perturbo_per_guide), rd(args.perturbo_per_element))
        uns_add = {"per_guide_results": g, "per_element_results": e}
        write_tsv(g, out / "per_guide_output.tsv.gz")
        write_tsv(e, out / "per_element_output.tsv.gz")
    elif args.mode == "cis-trans":
        uns_add = {"cis_per_guide_results": rd(args.cis_per_guide), "cis_per_element_results": rd(args.cis_per_element),
                   "trans_per_guide_results": rd(args.trans_per_guide), "trans_per_element_results": rd(args.trans_per_element)}
        for k, v in uns_add.items():
            write_tsv(v, out / f"{k}.tsv.gz")
    elif args.mode == "chunks":
        g = merge_chunk_tables(args.per_guide_files, ["gene_id", "guide_id"])
        e = merge_chunk_tables(args.per_element_files, ["gene_id", "intended_target_name"])
        uns_add = {"per_guide_results": g, "per_element_results": e}
        if args.chunk_manifest:
            uns_add["sceptre_chunk_manifest"] = pd.read_csv(args.chunk_manifest, sep="\t")
        write_tsv(g, out / "sceptre_per_guide_output.tsv.gz")
        write_tsv(e, out / "sceptre_per_element_output.tsv.gz")
    elif args.mode == "legacy":
        cis, trans = read_h5mu(Path(args.cis_mudata)), read_h5mu(Path(args.trans_mudata))
        base = trans
        uns_add = {"cis_test_results": uns_df(cis, "test_results"), "trans_test_results": uns_df(trans, "test_results")}
        base.uns.pop("test_results", None)
    elif args.mode == "add-test-results":
        tr = pd.read_csv(args.test_results_csv).rename(columns={"log2_fc": "sceptre_log2_fc", "p_value": "sceptre_p_value"})
        uns_add = {"test_results": tr}
    elif args.mode == "export":
        key = resolve_results_key(base, args.results_key)
        pg, pe = export_outputs(base, uns_df(base, key))
        write_tsv(pg, out / "per_guide_output.tsv")
        write_tsv(pe, out / "per_element_output.tsv")
    if base is not None and args.mode != "export":
        if args.mode == "cis-trans":
            base.uns.pop("per_guide_results", None)
            base.uns.pop("per_element_results", None)
        base.uns.update(uns_add)
        write_h5mu(base, out / "inference_mudata.h5mu")
    finish(out, f"Merge results ({args.mode})", {"mode": args.mode, "tables": {k: int(len(v)) for k, v in uns_add.items() if v is not None}},
           [md_table(["uns key", "rows"], [[k, len(v)] for k, v in uns_add.items() if v is not None])])
    return 0


# ---------------------------------------------------------------------------
# Additional QC: mapping_gene.py, mapping_guide.py, intended_target.py, trans.py, generate_report.py
# ---------------------------------------------------------------------------

def column_stats(s) -> dict:
    return {"median": s.median(), "mean": s.mean(), "std": s.std(), "min": s.min(), "max": s.max(), "q25": s.quantile(0.25), "q75": s.quantile(0.75)}


def gene_metrics(a, label: str, umi_col="total_gene_umis", genes_col="num_expressed_genes", mito_col="percent_mito") -> dict:
    m = {"batch": label, "n_cells": int(a.n_obs)}
    for col, pre in ((umi_col, "umi"), (genes_col, "genes"), (mito_col, "mito")):
        if col in a.obs.columns:
            for k, v in column_stats(a.obs[col].astype(float)).items():
                m[f"{pre}_{k}"] = v
    return m


def guide_metrics(g, label: str, include_per_guide: bool, umi_col: str = "total_guide_umis") -> dict:
    m = {"batch": label, "n_cells": int(g.n_obs)}
    if umi_col in g.obs.columns:
        for k, v in column_stats(g.obs[umi_col].astype(float)).items():
            m[f"guide_umi_{k}"] = v
    gpc = g.obs["n_guides_per_cell"]
    m.update(guides_per_cell_mean=gpc.mean(), guides_per_cell_std=gpc.std(), guides_per_cell_min=gpc.min(), guides_per_cell_max=gpc.max(),
             guides_per_cell_median=gpc.median(), n_cells_with_guide=int((gpc > 0).sum()), n_cells_exactly_1_guide=int((gpc == 1).sum()),
             frac_cells_with_guide=float((gpc > 0).sum() / g.n_obs) if g.n_obs else 0.0)
    if include_per_guide:
        cpg = g.var["n_cells_per_guide"]
        m.update(n_guides_total=int(g.n_vars), cells_per_guide_median=cpg.median(), cells_per_guide_mean=cpg.mean(), cells_per_guide_std=cpg.std(),
                 cells_per_guide_min=cpg.min(), cells_per_guide_max=cpg.max())
    return m


def qc_gene(m: MData, outdir: Path, batch_col: str = "batch", plots: bool = True):
    pd, np = _pd(), _np()
    g = m["gene"]
    rows = [gene_metrics(g, "all")]
    if batch_col in g.obs.columns:
        for b in sorted(g.obs[batch_col].astype(str).unique()):
            rows.append(gene_metrics(g[g.obs[batch_col].astype(str).values == b], str(b)))
    df = pd.DataFrame(rows)
    write_tsv(df, outdir / "gene_metrics.tsv")
    plt = _plt() if plots else None
    if plt:
        s = np.sort(row_sums(g.X))[::-1]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(np.log1p(s), np.arange(len(s)), color=BLUE)
        _style(ax, f"Knee plot (N={g.n_obs} cells)", "Log(UMI counts + 1)", "Barcode rank")
        _save(fig, outdir / "gene_knee_plot.png")
        cols = [c for c in ("total_gene_umis", "num_expressed_genes", "percent_mito") if c in g.obs.columns]
        if cols:
            fig, axes = plt.subplots(len(cols), 1, figsize=(7, 2.4 * len(cols)))
            axes = np.atleast_1d(axes)
            for ax, c in zip(axes, cols):
                if batch_col in g.obs.columns:
                    for k, b in enumerate(sorted(g.obs[batch_col].astype(str).unique())):
                        v = g.obs.loc[g.obs[batch_col].astype(str) == b, c]
                        ax.hist(v, bins=60, histtype="step", density=True, color=[BLUE, ORANGE, GREEN, RED][k % 4], label=f"{b}: {v.median():.0f}")
                    ax.legend(frameon=False, fontsize=7)
                else:
                    ax.hist(g.obs[c], bins=60, color=BLUE)
                _style(ax, c, "", "")
            _save(fig, outdir / "gene_histograms_by_batch.png")
    return df


def qc_guide(m: MData, outdir: Path, batch_col: str = "batch", plots: bool = True):
    pd, np = _pd(), _np()
    g = m["guide"].copy()
    if "guide_assignment" not in g.layers:
        raise SystemExit("Layer 'guide_assignment' not found in guide.layers")
    A = as_csr(g.layers["guide_assignment"]) > 0
    g.obs["n_guides_per_cell"] = np.asarray(A.sum(axis=1)).ravel()
    g.var["n_cells_per_guide"] = np.asarray(A.sum(axis=0)).ravel()
    rows = [guide_metrics(g, "all", True)]
    if batch_col in g.obs.columns:
        for b in sorted(g.obs[batch_col].astype(str).unique()):
            rows.append(guide_metrics(g[g.obs[batch_col].astype(str).values == b], str(b), False))
    df = pd.DataFrame(rows)
    write_tsv(df, outdir / "guide_metrics.tsv")
    plt = _plt() if plots else None
    if plt:
        fig, axes = plt.subplots(3, 1, figsize=(7, 7))
        if "total_guide_umis" in g.obs.columns:
            axes[0].hist(g.obs["total_guide_umis"], bins=80, color=BLUE)
        _style(axes[0], "Guide UMI counts per cell", "total_guide_umis", "")
        axes[1].hist(g.obs["n_guides_per_cell"], bins=30, color=ORANGE)
        _style(axes[1], f"Guides assigned per cell (mean {g.obs['n_guides_per_cell'].mean():.2f})", "n_guides_per_cell", "")
        lab = g.var["type"].astype(str) if "type" in g.var.columns else pd.Series("guide", index=g.var.index)
        for k, t in enumerate(sorted(lab.unique())):
            axes[2].hist(g.var.loc[lab == t, "n_cells_per_guide"], bins=30, histtype="step", color=[BLUE, ORANGE, GREEN][k % 3], label=t)
        axes[2].legend(frameon=False, fontsize=7)
        _style(axes[2], "Cells per guide", "n_cells_per_guide", "")
        _save(fig, outdir / "guide_histograms.png")
    return df


def intended_target_qc(res, guide_var, log2fc_col: str, pvalue_col: str, fc_threshold: float = 0.4, pval_threshold: float = 0.05,
                       non_targeting_name: str = NT_LABEL, seed: int = 42, match_symbol: Optional[dict] = None):
    """intended_target.py: filter to gene_id == intended_target_name, knockdown metrics, AUROC/AUPRC vs NT guides."""
    np, pd = _np(), _pd()
    gm = guide_meta(guide_var, non_targeting_name)
    imap = gm.set_index("guide_id")["intended_target_name"].to_dict()
    tmap = gm.set_index("guide_id")["targeting"].to_dict()
    r = res.copy()
    r["intended_target_name"] = r["guide_id"].map(imap)
    r["targeting"] = r["guide_id"].map(tmap)
    gene_key = r["gene_id"].astype(str)
    if match_symbol:
        gene_key = gene_key.map(lambda x: match_symbol.get(x, x))
    direct = (r["gene_id"].astype(str) == r["intended_target_name"].astype(str)) | (gene_key == r["intended_target_name"].astype(str))
    it = r[direct].drop_duplicates(["guide_id", "gene_id"])
    intended = gm.merge(it[["guide_id", "gene_id", log2fc_col, pvalue_col]], on="guide_id", how="left")
    valid = intended.dropna(subset=[log2fc_col, pvalue_col])
    n = len(valid)
    l2 = math.log2(fc_threshold)
    strong = valid[log2fc_col] <= l2
    sig = valid[pvalue_col] < pval_threshold
    metrics = {"n_guides_total": int(len(intended)), "n_guides_tested": int(n), "fc_threshold": fc_threshold, "log2fc_threshold": l2,
               "pval_threshold": pval_threshold, "n_strong_knockdowns": int(strong.sum()), "n_significant": int(sig.sum()),
               "n_strong_and_significant": int((strong & sig).sum()), "frac_strong_knockdowns": float(strong.sum() / n) if n else 0.0,
               "frac_significant": float(sig.sum() / n) if n else 0.0, "frac_strong_and_significant": float((strong & sig).sum() / n) if n else 0.0,
               "median_log2fc": float(valid[log2fc_col].median()) if n else float("nan"), "mean_log2fc": float(valid[log2fc_col].mean()) if n else float("nan")}
    pos = r[direct & r[pvalue_col].notna()].copy()
    pos["direct_target"] = 1
    targets = set(r.loc[r["targeting"] == True, "intended_target_name"].dropna().astype(str))  # noqa: E712
    neg = r[(r["targeting"] == False) & (r["gene_id"].astype(str).isin(targets) | gene_key.isin(targets)) & r[pvalue_col].notna()].copy()  # noqa: E712
    neg["direct_target"] = 0
    curves = {}
    if len(pos) and len(neg):
        if len(neg) > len(pos):
            neg = neg.sample(n=len(pos), random_state=seed)
        ev = pd.concat([pos, neg], ignore_index=True)
        auroc, auprc, curves = roc_pr(ev["direct_target"].values, 1 - ev[pvalue_col].values)
        metrics.update(auroc=auroc, auprc=auprc, n_eval_positives=int(len(pos)), n_eval_negatives=int(len(neg)))
    else:
        metrics.update(auroc=float("nan"), auprc=float("nan"), n_eval_positives=0, n_eval_negatives=0)
    return intended, metrics, curves


def trans_qc(res, guide_var, log2fc_col: str, pvalue_col: str, pval_threshold: float = 0.05, fdr_method: str = "fdr_bh",
             non_targeting_name: str = NT_LABEL, validated_links=None, seed: int = 42):
    """trans.py: BH (or bonferroni / none) over all tests, per-guide trans counts excluding the intended target."""
    np, pd = _np(), _pd()
    r = res.copy()
    p = r[pvalue_col].values.astype(float)
    if fdr_method == "none":
        r["p_value_adj"] = p
    elif fdr_method == "bonferroni":
        r["p_value_adj"] = np.minimum(p * np.isfinite(p).sum(), 1.0)
    else:
        r["p_value_adj"] = bh_adjust(p)
    gm = guide_meta(guide_var, non_targeting_name)
    imap = gm.set_index("guide_id")["intended_target_name"].to_dict()
    r["intended_target_name"] = r["guide_id"].map(imap)
    rr = r[r["gene_id"] != r["intended_target_name"]].dropna(subset=[log2fc_col, "p_value_adj"])
    sig = rr["p_value_adj"] < pval_threshold
    grp = rr.assign(sig=sig, up=sig & (rr[log2fc_col] > 0), down=sig & (rr[log2fc_col] < 0), absfc=rr[log2fc_col].abs()).groupby("guide_id")
    per = pd.DataFrame({"n_tests": grp.size(), "n_significant_trans": grp["sig"].sum(), "n_upregulated": grp["up"].sum(), "n_downregulated": grp["down"].sum(),
                        "median_log2fc": grp[log2fc_col].median(), "mean_log2fc": grp[log2fc_col].mean(), "max_abs_log2fc": grp["absfc"].max()}).reset_index()
    per["frac_significant"] = per["n_significant_trans"] / per["n_tests"].where(per["n_tests"] > 0, 1)
    per = per.merge(gm, on="guide_id", how="left")
    t, nt = per[per["targeting"] == True], per[per["targeting"] == False]  # noqa: E712
    med = lambda s: float(s.median()) if len(s) else float("nan")  # noqa: E731
    metrics = {"n_guides_tested": int(len(per)), "n_targeting_guides": int(len(t)), "n_non_targeting_guides": int(len(nt)),
               "median_significant_per_guide_targeting": med(t["n_significant_trans"]), "mean_significant_per_guide_targeting": float(t["n_significant_trans"].mean()) if len(t) else float("nan"),
               "median_significant_per_guide_nt": med(nt["n_significant_trans"]), "mean_significant_per_guide_nt": float(nt["n_significant_trans"].mean()) if len(nt) else float("nan"),
               "total_significant_tests": int(per["n_significant_trans"].sum()), "median_genome_log2fc_targeting": med(t["median_log2fc"]),
               "median_genome_log2fc_nt": med(nt["median_log2fc"]), "fdr_method": fdr_method, "auroc": float("nan"), "auprc": float("nan"),
               "n_validated_links": 0, "n_eval_positives": 0, "n_eval_negatives": 0}
    curves = {}
    if validated_links is not None and len(validated_links) and "gene_name" in guide_var.columns:
        gvr = guide_var.reset_index(drop=True)
        id2sym = gvr.drop_duplicates("intended_target_name").set_index("intended_target_name")["gene_name"]
        tgt_sym = gm.set_index("guide_id")["gene_name"].to_dict() if "gene_name" in gm.columns else {}
        vl = validated_links[validated_links["guide_target"].isin(set(guide_var["gene_name"].dropna()))]
        pairs = set(zip(vl["guide_target"], vl["gene"]))
        genes = set(vl["gene"])
        r2 = r.copy()
        r2["gene_name"] = r2["gene_id"].map(id2sym)
        r2["target_name"] = r2["guide_id"].map(tgt_sym)
        r2["targeting"] = r2["guide_id"].map(gm.set_index("guide_id")["targeting"].to_dict())
        isv = [(a, b) in pairs for a, b in zip(r2["target_name"], r2["gene_name"])]
        posd = r2[pd.Series(isv, index=r2.index) & (r2["targeting"] == True) & r2[pvalue_col].notna()]  # noqa: E712
        negd = r2[(r2["targeting"] == False) & r2["gene_name"].isin(genes) & r2[pvalue_col].notna()]  # noqa: E712
        metrics["n_validated_links"] = int(len(vl))
        if len(posd) and len(negd):
            if len(negd) > len(posd):
                negd = negd.sample(n=len(posd), random_state=seed)
            lab = np.r_[np.ones(len(posd)), np.zeros(len(negd))]
            sc = 1 - np.r_[posd[pvalue_col].values, negd[pvalue_col].values]
            metrics["auroc"], metrics["auprc"], curves = roc_pr(lab, sc)
            metrics.update(n_eval_positives=int(len(posd)), n_eval_negatives=int(len(negd)))
    return r, per, metrics, curves


def _plot_rocpr(curves, auroc, auprc, path: Path, title: str = ""):
    plt = _plt()
    if not plt or not curves:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    axes[0].plot(curves["fpr"], curves["tpr"], color=BLUE, lw=2)
    axes[0].plot([0, 1], [0, 1], "--", color=GREY, lw=1)
    _style(axes[0], f"ROC (AUROC = {auroc:.3f})", "False Positive Rate", "True Positive Rate")
    axes[1].plot(curves["recall"], curves["precision"], color=RED, lw=2)
    axes[1].axhline(0.5, ls="--", color=GREY, lw=1)
    _style(axes[1], f"PR (AUPRC = {auprc:.3f})", "Recall", "Precision")
    if title:
        fig.suptitle(title, fontsize=9)
    return _save(fig, path)


def _plot_volcano(df, lfc, pcol, path: Path, title: str, highlight=None, thr: float = 0.05, n_label: int = 10, label_col: str = "guide_id"):
    np = _np()
    plt = _plt()
    if not plt:
        return None
    d = df.dropna(subset=[lfc, pcol])
    if not len(d):
        return None
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    y = -np.log10(d[pcol].astype(float).clip(lower=1e-300))
    col = np.where(highlight.reindex(d.index).fillna(False).values, RED, GREY) if highlight is not None else BLUE
    ax.scatter(d[lfc], y, s=8, c=col, alpha=0.75, edgecolors="none")
    ax.axhline(-math.log10(thr), color=GREY, ls="--", lw=1)
    ax.axvline(0, color=AXIS, lw=1)
    if label_col in d.columns:
        for _, row in d.nsmallest(n_label, pcol).iterrows():
            ax.annotate(str(row[label_col]), (row[lfc], -math.log10(max(row[pcol], 1e-300))), fontsize=6, color=INK2)
    _style(ax, title, "log2 fold change", "-log10(p-value)")
    return _save(fig, path)


QC_REPORT_CSS = """<style>
:root{--bg:#fcfcfb;--fg:#0b0b0b;--muted:#52514e;--line:#d0cfca;--accent:#2a78d6}
@media (prefers-color-scheme: dark){:root{--bg:#161615;--fg:#f2f1ed;--muted:#b8b6af;--line:#3a3935;--accent:#6aa7ef}}
body{background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1100px;margin:0 auto;padding:16px}
h1{font-size:22px}h2{font-size:17px;border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:28px}h3{font-size:14px;color:var(--muted)}
table{border-collapse:collapse;font-size:12px;margin:8px 0;display:block;overflow-x:auto}th,td{border:1px solid var(--line);padding:3px 7px;text-align:right}
th{background:rgba(127,127,127,.08)}td:first-child,th:first-child{text-align:left}img{max-width:100%;border:1px solid var(--line);margin:6px 0;background:#fff}
.desc{color:var(--muted)}.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:12px 0}.tabs a{color:var(--accent);text-decoration:none;border:1px solid var(--line);padding:3px 9px;border-radius:12px}
.flow span{display:block}.flow-step{border:1px solid var(--line);padding:6px 10px;margin:4px 0;border-radius:6px}.flow-arrow{text-align:center;color:var(--muted)}
</style>"""


def img_tag(path: Path, caption: str = "") -> str:
    if not path or not Path(path).is_file():
        return ""
    b64 = base64.b64encode(Path(path).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{b64}" alt="{caption}"><figcaption class="desc">{caption}</figcaption></figure>'


def table_html(df, max_rows: int = 50) -> str:
    if df is None or not len(df):
        return "<p class='desc'>(empty)</p>"
    d = df.head(max_rows).copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: _fmt(v))
    return d.to_html(index=False, escape=True, border=0)


def qc_report_html(sample: str, sections: "list[tuple[str, str, list]]") -> str:
    parts = [f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
             f"<title>QC Report</title>{QC_REPORT_CSS}</head><body><h1>QC Summary Report</h1><div class='desc'>{sample} · generated {time.strftime('%B %d, %Y')} "
             f"by igvfagent crispr-pipeline (port of {UPSTREAM_REPO})</div>"]
    for title, desc, items in sections:
        parts.append(f"<h2>{title}</h2><p class='desc'>{desc}</p>")
        for it in items:
            parts.append(it)
    parts.append("</body></html>")
    return "".join(parts)


def run_qc(m: MData, out: Path, results_key: str = "auto", log2fc_col: str = "auto", pvalue_col: str = "auto", fc_threshold: float = 0.4,
           pval_threshold: float = 0.05, fdr_method: str = "fdr_bh", validated_links: Optional[str] = None, batch_col: str = "batch",
           plots: bool = True, match_symbol="auto", sample: str = "Sample") -> dict:
    pd = _pd()
    d_gene, d_guide, d_it, d_tr = out / "gene", out / "guide", out / "intended_target", out / "trans"
    for d in (d_gene, d_guide, d_it, d_tr):
        d.mkdir(parents=True, exist_ok=True)
    gm_df = qc_gene(m, d_gene, batch_col, plots)
    gd_df = qc_guide(m, d_guide, batch_col, plots)
    summ: "dict[str, Any]" = {"gene": gm_df.iloc[0].to_dict(), "guide": gd_df.iloc[0].to_dict()}
    has = any(k in m.uns for k in RESULT_KEYS)
    sections = [("1. Gene Expression Mapping QC", "UMI counts, detected genes and mitochondrial content per cell, overall and per batch.",
                 [table_html(gm_df), img_tag(d_gene / "gene_histograms_by_batch.png", "distributions by batch"), img_tag(d_gene / "gene_knee_plot.png", "knee plot")]),
                ("2. Guide Mapping QC", "Guide UMIs, guides assigned per cell (guide_assignment layer) and cells per guide.",
                 [table_html(gd_df), img_tag(d_guide / "guide_histograms.png", "guide distributions")])]
    if has:
        key = resolve_results_key(m, results_key)
        res = uns_df(m, key)
        lfc, pc = resolve_metric_columns(res, log2fc_col, pvalue_col)
        sym_all = dict(zip(map(str, m["gene"].var_names), _gene_symbols(m["gene"])))
        if match_symbol == "auto":
            names = set(m["guide"].var["intended_target_name"].astype(str))
            match_symbol = not (names & set(sym_all)) and bool(names & set(sym_all.values()))
        sym = sym_all if match_symbol else None
        intended, it_metrics, curves = intended_target_qc(res, m["guide"].var, lfc, pc, fc_threshold, pval_threshold, match_symbol=sym)
        write_tsv(pd.DataFrame([it_metrics]), d_it / "intended_target_metrics.tsv")
        write_tsv(intended[[c for c in ["guide_id", "intended_target_name", "gene_name", "targeting", lfc, pc] if c in intended.columns]], d_it / "intended_target_results.tsv")
        if plots:
            _plot_rocpr(curves, it_metrics["auroc"], it_metrics["auprc"], d_it / "intended_target_roc_pr_curves.png")
            _plot_volcano(intended, lfc, pc, d_it / "intended_target_volcano.png", f"Intended target knockdown (N={intended[pc].notna().sum()} guides)", thr=pval_threshold)
        vl = pd.read_csv(validated_links, sep="\t") if validated_links else None
        tr_all, per, tr_metrics, tcurves = trans_qc(res, m["guide"].var, lfc, pc, pval_threshold, fdr_method, validated_links=vl)
        write_tsv(pd.DataFrame([tr_metrics]), d_tr / "trans_metrics.tsv")
        write_tsv(per, d_tr / "trans_per_guide_summary.tsv")
        sym_map = dict(zip(map(str, m["gene"].var_names), _gene_symbols(m["gene"])))
        ann = tr_all.rename(columns={"gene_id": "tested_gene_id"})
        ann["tested_gene_symbol"] = ann["tested_gene_id"].astype(str).map(sym_map)
        ann["guide_target_id"] = ann["intended_target_name"]
        ann["significant"] = ann["p_value_adj"] < pval_threshold
        cols = [c for c in ["guide_id", "guide_target_id", "tested_gene_id", "tested_gene_symbol", lfc, pc, "p_value_adj", "significant"] if c in ann.columns]
        write_tsv(ann[cols], d_tr / "trans_results.tsv.gz")
        write_tsv(ann.loc[ann["significant"], cols], d_tr / "trans_significant_results.tsv")
        if plots:
            _plot_volcano(tr_all.sample(n=min(len(tr_all), 10000), random_state=0), lfc, "p_value_adj", d_tr / "trans_volcano.png", "Trans effects (subsampled)", thr=pval_threshold)
            _plot_rocpr(tcurves, tr_metrics["auroc"], tr_metrics["auprc"], d_tr / "trans_roc_pr_curves.png")
        summ.update(results_key=key, log2fc_col=lfc, pvalue_col=pc, intended_target=it_metrics, trans=tr_metrics)
        sections += [("3. Intended Target Inference QC", f"{key}: tests where gene_id == intended_target_name; strong knockdown log2FC <= log2({fc_threshold}); "
                                                         "AUROC/AUPRC vs non-targeting guides tested on the same genes (balanced, seed 42).",
                      [table_html(pd.DataFrame([it_metrics])), img_tag(d_it / "intended_target_volcano.png", "volcano"),
                       img_tag(d_it / "intended_target_roc_pr_curves.png", "ROC / PR")]),
                     ("4. Trans QC", f"{fdr_method} over all tests; significant = adjusted p < {pval_threshold}; intended target excluded.",
                      [table_html(pd.DataFrame([tr_metrics])), table_html(per.sort_values("n_significant_trans", ascending=False), 20), img_tag(d_tr / "trans_volcano.png", "trans volcano")])]
    html = out / "qc_report.html"
    html.write_text(qc_report_html(sample, sections))
    print(f"Wrote: {html}")
    summ["html"] = str(html)
    return summ


def cmd_qc(args) -> int:
    out = run_dir(args.label)
    s = run_qc(read_h5mu(Path(args.mudata)), out, args.results_key, args.log2fc_col, args.pvalue_col, args.fc_threshold, args.pval_threshold,
               args.fdr_method, args.validated_links, args.batch_col, not args.no_plots, args.match_symbol, args.label)
    secs = [f"QC report: `{s['html']}`"]
    if "intended_target" in s:
        it = s["intended_target"]
        secs.append(md_table(["metric", "value"], [[k, it[k]] for k in ("n_guides_tested", "n_strong_knockdowns", "n_significant", "median_log2fc", "auroc", "auprc")]))
        tr = s["trans"]
        secs.append(md_table(["metric", "value"], [[k, tr[k]] for k in ("n_targeting_guides", "median_significant_per_guide_targeting", "median_significant_per_guide_nt", "total_significant_tests")]))
    finish(out, "Additional QC", s, secs)
    return 0


# ---------------------------------------------------------------------------
# Evaluation: evaluate_controls.py, igv.py, volcano_plot.py, network_plot(_undefined).py, select_nodes.py, utlis.py
# ---------------------------------------------------------------------------

def evaluate_controls(m: MData, out: Path, key: str = "trans_per_guide_results", pcol: str = "p_value", lfc: str = "log2_fc",
                      plots: bool = True, seed: int = 42, symbol_map: Optional[dict] = None) -> dict:
    np, pd = _np(), _pd()
    res = uns_df(m, key)
    if res is None:
        raise SystemExit(f"{key} not found in mdata.uns")
    if pcol not in res.columns:
        lfc, pcol = resolve_metric_columns(res)
    gv = m["guide"].var.copy()
    gv["targeting"] = gv["targeting"].map(lambda x: True if (x is True or str(x).upper() == "TRUE") else (False if (x is False or str(x).upper() == "FALSE") else x))
    imap = gv.set_index("guide_id")["intended_target_name"].astype(str).to_dict()
    tmap = gv.set_index("guide_id")["targeting"].to_dict()
    r = res.copy()
    r["intended_target_name"] = r["guide_id"].map(imap)
    gk = r["gene_id"].astype(str).map(lambda x: symbol_map.get(x, x)) if symbol_map else r["gene_id"].astype(str)
    direct = (r["gene_id"].astype(str) == r["intended_target_name"]) | (gk == r["intended_target_name"])
    pos = r[direct].drop_duplicates().copy()
    pos["direct_target"] = 1
    r["targeting_genes"] = r["guide_id"].map(tmap)
    targets = set(r.loc[r["targeting_genes"] == True, "intended_target_name"].dropna())  # noqa: E712
    ntc = r[r["targeting_genes"] == False]  # noqa: E712
    matched = ntc[ntc["gene_id"].astype(str).isin(targets) | gk.reindex(ntc.index).isin(targets)]
    warn = ""
    if len(matched):
        ntc = matched
    else:
        warn = "No non-targeting guides were tested against the intended targets; using a random sample of non-targeting tests instead."
    ntc = ntc.copy()
    ntc["direct_target"] = 0
    n = min(len(pos), len(ntc))
    ev = pd.concat([pos, ntc.sample(n=n, random_state=seed)]) if n else pos
    ev = ev.dropna(subset=[pcol])
    auroc, auprc, curves = roc_pr(ev["direct_target"].values, 1 - ev[pcol].values) if n else (float("nan"), float("nan"), {})
    write_tsv(ev[[c for c in ["guide_id", "gene_id", "intended_target_name", lfc, pcol, "direct_target"] if c in ev.columns]], out / "controls_evaluation_table.tsv")
    if plots and _plt():
        plt = _plt()
        fig, ax = plt.subplots(figsize=(5, 5.5))
        d = ev.dropna(subset=[lfc, pcol])
        ax.scatter(d[lfc], -np.log10(d[pcol].replace(0, np.nan)), c=np.where(d["direct_target"] == 1, RED, GREY), s=16, alpha=0.75, edgecolors="none")
        for x in (1, -1):
            ax.axvline(x, color=INK2, ls="--", lw=1, alpha=0.6)
        ax.axhline(-math.log10(0.05), color=INK2, ls="--", lw=1, alpha=0.6)
        _style(ax, "Volcano Plot (red: guides on direct target, grey: non-targeting)", "log2 Fold Change", "-log10(p-value)")
        _save(fig, out / "trans_perturbo_volcano_plot.png")
        _plot_rocpr(curves, auroc, auprc, out / "trans_perturbo_precision_recall_roc.png", warn)
        fig, ax = plt.subplots(figsize=(4, 3))
        cnt = ev.groupby("direct_target")["guide_id"].count()
        ax.bar([str(i) for i in cnt.index], cnt.values, color=[GREY, RED][: len(cnt)])
        _style(ax, "Direct targets vs control guides", "direct_target", "Number of guides")
        _save(fig, out / "trans_perturbo_barplot_direct_vs_control.png")
    return {"results_key": key, "auroc": auroc, "auprc": auprc, "n_positives": int((ev["direct_target"] == 1).sum()),
            "n_negatives": int((ev["direct_target"] == 0).sum()), "warning": warn}


def igv_tracks(m: MData, gtf, key: str, method: Optional[str]) -> "tuple[Any, Any]":
    """igv.py: bedgraph for elements tested on their own gene (promoter), bedpe element -> gene otherwise."""
    np, pd = _np(), _pd()
    lfc, pc = ("log2_fc", "p_value") if not method else (f"{method}_log2_fc", f"{method}_p_value")
    coords = {}
    gvar = m["gene"].var
    if {"gene_chr", "gene_start", "gene_end"} <= set(gvar.columns):
        for idx, r in gvar.iterrows():
            if pd.notna(r["gene_start"]) and pd.notna(r["gene_end"]):
                coords[str(idx)] = [r["gene_chr"], r["gene_start"], r["gene_end"]]
    for _, r in m["guide"].var.iterrows():
        n = str(r["intended_target_name"])
        if n in coords or n == NT_LABEL:
            continue
        coords[n] = [r.get("intended_target_chr"), r.get("intended_target_start"), r.get("intended_target_end")]
    if gtf is not None:
        names = gtf[["gene_id2", "gene_name"]].drop_duplicates().rename(columns={"gene_id2": "_gid"})
    else:
        names = pd.DataFrame({"_gid": list(map(str, m["gene"].var_names)), "gene_name": _gene_symbols(m["gene"])})
    res = uns_df(m, key).merge(names, left_on="gene_id", right_on="_gid", how="left").dropna(subset=[lfc, pc])
    bp, bg = defaultdict(list), defaultdict(list)
    for _, r in res.iterrows():
        t = str(r["intended_target_name"])
        if t == str(r["gene_name"]) or t == str(r["gene_id"]):
            if t in coords:
                c = coords[t]
                for k, v in zip(("chr", "start", "end", "p_value", "log2_fc"), (c[0], c[1], c[2], r[pc], r[lfc])):
                    bg[k].append(v)
        elif t in coords and str(r["gene_id"]) in coords:
            s, g = coords[t], coords[str(r["gene_id"])]
            for k, v in zip(("chr1", "start1", "end1", "chr2", "start2", "end2", "p_value", "log2_fc"), (s[0], s[1], s[2], g[0], g[1], g[2], r[pc], r[lfc])):
                bp[k].append(v)
    return pd.DataFrame(bp), pd.DataFrame(bg)


def select_central_nodes(res, num: int, source: str, target: str, weight: str, min_weight: Optional[float]) -> "list[str]":
    """network_plot_undefined.py: nodes ranked by weighted degree (|log2FC|), restricted to intended targets."""
    r = res.drop_duplicates()
    if min_weight is not None:
        r = r[r[weight].abs() >= min_weight]
    deg: Counter = Counter()
    for s, t, w in zip(r[source], r[target], r[weight]):
        if w == w:
            deg[s] += abs(w)
            deg[t] += abs(w)
    names = set(r["intended_target_name"])
    return [n for n, _ in sorted(deg.items(), key=lambda x: -x[1]) if n in names][:num]


def plot_network(res, node: str, weight: str, pcol: str, path: Path, source: str = "intended_target_name", target: str = "gene_id",
                 min_weight: Optional[float] = 0.1, title: str = ""):
    np = _np()
    plt = _plt()
    if not plt:
        return None
    r = res.drop_duplicates()
    if min_weight is not None:
        r = r[r[weight].abs() >= min_weight]
    r = r[r[source] == node].dropna(subset=[weight, pcol])
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.axis("off")
    if not len(r):
        ax.text(0.5, 0.5, f"No network found for\n{node}", ha="center", va="center")
    else:
        r = r.nsmallest(40, pcol)
        n = len(r)
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        xs, ys = np.cos(ang), np.sin(ang)
        vmax = max(float(r[weight].abs().max()), 1e-6)
        cmap = plt.cm.coolwarm
        for x, y_, (_, row) in zip(xs, ys, r.iterrows()):
            ax.plot([0, x], [0, y_], color=cmap(0.5 + row[weight] / (2 * vmax)), lw=1.5)
            ax.scatter([x], [y_], s=30 + 40 * min(-math.log10(max(row[pcol], 1e-300)), 10), color="#9cc9ee", zorder=3)
            ax.text(x * 1.12, y_ * 1.12, f"{row[target]}\nFC:{row[weight]:.2f}", fontsize=5, ha="center", va="center")
        ax.scatter([0], [0], s=300, color=ORANGE, zorder=4)
        ax.text(0, 0, node, fontsize=7, ha="center", va="center")
    ax.set_title(title or f"Network - {node}", fontsize=9, loc="left")
    return _save(fig, path)


def perturbation_expression(m: MData, guide_id: str, gene_id: str, layer: Optional[str] = None, non_targeting_label: str = NT_LABEL):
    """utlis.py::get_perturbation_expression: expression of gene_id in cells with guide_id vs cells with only NT guides."""
    np, pd = _np(), _pd()
    g, e = m["guide"], m["gene"]
    if guide_id not in set(map(str, g.var_names)):
        raise SystemExit(f"guide_id '{guide_id}' not found in guide.var_names")
    if gene_id not in set(map(str, e.var_names)):
        raise SystemExit(f"gene_id '{gene_id}' not found in gene.var_names")
    X = e.layers[layer] if layer else e.X
    expr = dense_col(as_csr(X), list(map(str, e.var_names)).index(gene_id))
    A = as_csr(g.layers["guide_assignment"])
    names = list(map(str, g.var_names))
    targ = guide_targeting(g.var, non_targeting_label).values
    with_g = dense_col(A, names.index(guide_id)) > 0
    nt = row_sums(A[:, np.where(~targ)[0]]) > 0
    tg = row_sums(A[:, np.where(targ)[0]]) > 0
    ctrl = nt & ~tg
    df = pd.DataFrame({"expression": np.r_[expr[with_g], expr[ctrl]], "group": [guide_id] * int(with_g.sum()) + [non_targeting_label] * int(ctrl.sum())},
                      index=pd.Index(list(e.obs_names[with_g]) + list(e.obs_names[ctrl]), name="barcode"))
    return df


def run_evaluate(m: MData, out: Path, gtf=None, plots: bool = True, num_nodes: int = 1, min_weight: float = 0.1, central_nodes=None) -> dict:
    pd = _pd()
    ev_dir = out / "evaluation_output"
    ev_dir.mkdir(parents=True, exist_ok=True)
    sym_all = dict(zip(map(str, m["gene"].var_names), _gene_symbols(m["gene"])))
    names = set(m["guide"].var["intended_target_name"].astype(str))
    sym = sym_all if (not names & set(sym_all) and names & set(sym_all.values())) else None
    summ: "dict[str, Any]" = {"figures": []}
    key = "trans_per_guide_results" if "trans_per_guide_results" in m.uns else resolve_results_key(m)
    summ["controls"] = evaluate_controls(m, ev_dir, key, plots=plots, symbol_map=sym)
    configs = [("cis_per_element_results", "cis"), ("trans_per_element_results", "trans")] if "cis_per_element_results" in m.uns else \
        [(k, None) for k in ("per_element_results", "test_results") if k in m.uns]
    tracks = {}
    for key, typ in configs:
        res = uns_df(m, key)
        methods = [None] if {"log2_fc", "p_value"} <= set(res.columns) else [x for x in ("sceptre", "perturbo") if f"{x}_log2_fc" in res.columns and res[f"{x}_log2_fc"].notna().any()]
        for meth in methods:
            bp, bg = igv_tracks(m, gtf, key, meth)
            pre = f"{typ}_" if typ else ""
            nm = meth or "evaluation"
            write_tsv(bp, ev_dir / f"{pre}{nm}.bedpe", header=False, quiet=True)
            write_tsv(bg, ev_dir / f"{pre}{nm}.bedgraph", header=False, quiet=True)
            print(f"Wrote: {ev_dir / f'{pre}{nm}.bedpe'}")
            tracks[f"{pre}{nm}"] = {"bedpe": len(bp), "bedgraph": len(bg)}
            lfc, pc = ("log2_fc", "p_value") if meth is None else (f"{meth}_log2_fc", f"{meth}_p_value")
            if plots:
                f = _plot_volcano(res, lfc, pc, ev_dir / f"{pre}{(meth + '_') if meth else ''}volcano_plot.png", f"{pre}{(meth or 'Differential Expression').capitalize()} Volcano Plot",
                                  label_col="intended_target_name")
                if f:
                    summ["figures"].append(str(f))
                nodes = central_nodes or select_central_nodes(res, num_nodes, "intended_target_name", "gene_id", lfc, min_weight)
                for node in nodes:
                    f = plot_network(res, node, lfc, pc, ev_dir / f"{pre}{(meth + '_') if meth else ''}network_{safe_label(node)}.png", min_weight=min_weight,
                                     title=f"{pre}{(meth or 'analysis').capitalize()} Network - {node}")
                    if f:
                        summ["figures"].append(str(f))
    summ["igv_tracks"] = tracks
    res0 = uns_df(m, key)
    pc0 = next((c for c in ("p_value", "perturbo_p_value", "sceptre_p_value") if c in res0.columns), None)
    if pc0 and "intended_target_name" in res0.columns:
        summ["top_nodes_by_p"] = list(res0.sort_values(pc0)["intended_target_name"].drop_duplicates().head(max(1, num_nodes)))
    return summ


def cmd_evaluate(args) -> int:
    out = run_dir(args.label)
    gtf = read_gtf(Path(args.gtf)) if args.gtf else None
    s = run_evaluate(read_h5mu(Path(args.mudata)), out, gtf, not args.no_plots, args.num_nodes, args.min_weight, args.central_nodes)
    c = s["controls"]
    finish(out, "Evaluation", s, [f"evaluate_controls on `{c['results_key']}`: AUROC {_fmt(c['auroc'])}, AUPRC {_fmt(c['auprc'])} "
                                  f"({c['n_positives']} direct-target tests vs {c['n_negatives']} non-targeting). {c['warning']}",
                                  md_table(["track", "bedpe rows", "bedgraph rows"], [[k, v["bedpe"], v["bedgraph"]] for k, v in s["igv_tracks"].items()])])
    return 0


def cmd_expression(args) -> int:
    out = run_dir(args.label)
    df = perturbation_expression(read_h5mu(Path(args.mudata)), args.guide_id, args.gene_id, args.layer)
    write_tsv(df.reset_index(), out / "perturbation_expression.tsv")
    st = df.groupby("group")["expression"].agg(["count", "mean", "median"]).reset_index()
    plt = _plt()
    if plt and not args.no_plots:
        fig, ax = plt.subplots(figsize=(4, 3.5))
        groups = list(st["group"])
        ax.violinplot([df.loc[df["group"] == g_, "expression"].values for g_ in groups], showmedians=True)
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups, fontsize=7)
        _style(ax, f"{args.gene_id}", "", "expression")
        _save(fig, out / "perturbation_expression.png")
    finish(out, f"Expression of {args.gene_id}: {args.guide_id} vs non-targeting-only cells", {"groups": st.to_dict(orient="records")},
           [md_table(list(st.columns), st.values.tolist())])
    return 0


# ---------------------------------------------------------------------------
# TF benchmark: tf_benchmark.py + tf_enrichment.py (pyranges / pybiomart not required)
# ---------------------------------------------------------------------------

def read_peaks(path: Path, threshold: float = 0.0):
    """ENCODE narrowPeak / bed.gz: keep peaks with Score >= the threshold quantile (upstream chipseq_threshold)."""
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df = df.rename(columns={0: "chr", 1: "start", 2: "end", 4: "score"})
    if "score" not in df.columns:
        df["score"] = 0.0
    df = df[df["score"] >= df["score"].quantile(threshold)]
    return df[["chr", "start", "end", "score"]]


def promoter_windows(gtf, width: int, genes=None, compat: bool = False):
    """Mean TSS per gene (strand-aware; --upstream-compat uses END for every gene), +/- width."""
    pd, np = _pd(), _np()
    g = gtf[gtf["feature"] == "gene"] if "feature" in gtf.columns else gtf
    g = g.drop_duplicates("gene_id2")
    if genes is not None:
        g = g[g["gene_id2"].isin(set(genes))]
    valid = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    g = g[g["seqname"].isin(valid)]
    tss = g["end"] if compat else np.where(g["strand"] == "+", g["start"], g["end"])
    t = pd.DataFrame({"gene_id": g["gene_id2"].values, "chr": g["seqname"].values, "TSS": np.asarray(tss)})
    ok = t.groupby("gene_id")["chr"].nunique() == 1
    t = t[t["gene_id"].isin(ok[ok].index)].groupby(["gene_id", "chr"], as_index=False)["TSS"].mean()
    t["TSS"] = t["TSS"].round().astype(int)
    t["start"], t["end"] = t["TSS"] - width, t["TSS"] + width
    return t


def tf_enrichment(results, peaks: "dict[str, Any]", gtf, width: int, genes=None, element_col: str = "intended_target_name", gene_col: str = "gene_id",
                  pval_col: str = "p_value", fdr_corrected: bool = True, cutoff: float = 0.001, compat: bool = False):
    """Fisher (greater) of significant tests vs ChIP-seq promoter targets, per TF (tf_enrichment.run_tf_enrichment)."""
    pd, np = _pd(), _np()
    from scipy.stats import fisher_exact  # type: ignore
    from eqtl_enrichment_skill import overlap_join  # type: ignore
    prom = promoter_windows(gtf, width, genes, compat)
    rows, tables = [], []
    for tf, pk in peaks.items():
        if tf not in set(results[element_col].astype(str)):
            print(f"Warning: TF {tf} not found in results (upstream raises); skipped")
            continue
        j = overlap_join(prom[["chr", "start", "end", "gene_id"]], pk[["chr", "start", "end"]])
        tg = set(j["gene_id"])
        sub = results[results[element_col].astype(str) == tf].copy()
        sub = sub[sub[gene_col].astype(str).isin(set(prom["gene_id"]))]
        sub["target"] = sub[gene_col].astype(str).isin(tg)
        if fdr_corrected:
            sub["significant"] = sub[pval_col] <= cutoff
        else:
            sub["significant"] = bh_adjust(sub[pval_col].values) <= cutoff
        t, s = sub["target"].values.astype(bool), sub["significant"].values.astype(bool)
        if len(set(t)) < 2 or len(set(s)) < 2:
            orv, p = float("nan"), float("nan")
        else:
            orv, p = fisher_exact([[int((t & s).sum()), int((~t & s).sum())], [int((t & ~s).sum()), int((~t & ~s).sum())]], alternative="greater")
        rows.append({"TF": tf, "Method": sub["method"].iloc[0] if "method" in sub.columns and len(sub) else "PerTurbo", "num_rejections": int(s.sum()),
                     "num_targets": int(t.sum()), "num_overlapping": int((t & s).sum()), "odds_ratio": float(orv), "pvalue": float(p),
                     "recall": float((t & s).sum() / t.sum()) if t.sum() else float("nan")})
        tables.append(sub.assign(TF=tf))
    return pd.DataFrame(rows), (pd.concat(tables, ignore_index=True) if tables else pd.DataFrame())


def run_tf_benchmark(m: MData, gtf, peak_files: "dict[str, str]", out: Path, windows=PROMOTER_WINDOWS, cutoff: float = 0.001,
                     compat: bool = False, plots: bool = True, name_map: Optional[dict] = None) -> dict:
    pd, np = _pd(), _np()
    tables = out / "benchmark_tables"
    tables.mkdir(parents=True, exist_ok=True)
    res = uns_df(m, "trans_per_element_results")
    if res is None:
        raise SystemExit("Missing trans_per_element_results in MuData.uns")
    res = res.drop_duplicates().copy()
    if "method" not in res.columns:
        res["method"] = "PerTurbo"
    name_map = name_map or ENSEMBL_GENE_NAME_MAP
    write_tsv(pd.DataFrame([{"tf_id": k, "tf_name": name_map.get(k, k), "bed_file": v} for k, v in peak_files.items()]), tables / "tf_peak_mapping.tsv", quiet=True)
    write_tsv(res, tables / "trans_per_element_results_used.tsv", quiet=True)
    peaks = {k: read_peaks(Path(v), 0.0) for k, v in peak_files.items()}
    pcol = "p_value" if "p_value" in res.columns else resolve_metric_columns(res)[1]
    runs = []
    for w in windows:
        er, _ = tf_enrichment(res, peaks, gtf, w, list(map(str, m["gene"].var_names)), pval_col=pcol, cutoff=cutoff, compat=compat)
        if len(er):
            runs.append(er.assign(TF_display=er["TF"].map(name_map).fillna(er["TF"]), promoter_window_width=w))
    if not runs:
        raise SystemExit("No enrichment results were produced. Check inputs/paths.")
    allr = pd.concat(runs, ignore_index=True)
    allr["log_pvalue_clipped"] = np.log10(allr["pvalue"].astype(float).clip(lower=1e-100, upper=1.0))
    write_tsv(allr, tables / "enrichment_all.tsv")
    order = allr.groupby("TF_display")["pvalue"].min().sort_values().index.tolist()
    write_tsv(pd.DataFrame({"TF": order}), tables / "tf_order.tsv", quiet=True)
    fig_p = None
    plt = _plt() if plots else None
    if plt:
        fig, axes = plt.subplots(1, 3, figsize=(12, 0.6 * len(order) + 2), sharey=True)
        cols = ["#003648", "#006479", "#0096a0", "#49bcbc", "#96ced3"]
        h = 0.8 / max(len(windows), 1)
        for k, w in enumerate(windows):
            sub = allr[allr["promoter_window_width"] == w].set_index("TF_display").reindex(order)
            ypos = np.arange(len(order)) + k * h
            for ax, c in zip(axes, ["num_rejections", "odds_ratio", "log_pvalue_clipped"]):
                v = sub[c].astype(float).replace([np.inf, -np.inf], np.nan)
                if c == "odds_ratio":  # an infinite odds ratio (no significant non-target) is drawn at 1.2x the largest finite one
                    fin = allr[c].astype(float).replace([np.inf, -np.inf], np.nan).max()
                    v = sub[c].astype(float).replace(np.inf, (fin if fin == fin else 10.0) * 1.2)
                ax.barh(ypos, v.fillna(0).values, height=h, color=cols[k % len(cols)], label=str(w))
        axes[0].set_yticks(np.arange(len(order)) + 0.4 - h / 2)
        axes[0].set_yticklabels(order)
        _style(axes[0], "N downstream genes", "", "Transcription factor")
        axes[1].axvline(1, ls="--", color=GREY)
        _style(axes[1], "Odds ratio", "", "")
        axes[2].axvline(math.log10(0.05), ls="--", color=GREY)
        _style(axes[2], "log10(p-value)", "", "")
        axes[2].legend(title="Promoter window (bp)", fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        fig_p = str(_save(fig, out / "tf_benchmark.png"))
    best = allr.sort_values("pvalue").drop_duplicates("TF").to_dict(orient="records")
    return {"n_tfs": int(allr["TF"].nunique()), "windows": list(windows), "best_per_tf": best, "figure": fig_p}


def cmd_tf_benchmark(args) -> int:
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    gtf = read_gtf(Path(args.gtf))
    peaks = {}
    if args.peaks:
        for spec in args.peaks:
            k, _, v = spec.partition("=")
            peaks[k] = v
    else:
        d = Path(args.encode_bed_dir) if args.encode_bed_dir else DATA_ROOT / "encode_bed_files"
        for tf, fn in PEAK_META_DICT.items():
            if not (d / fn).is_file():
                raise SystemExit(f"Missing ENCODE BED file: {d / fn} (run `igvfagent crispr-pipeline setup`)")
            peaks[tf] = str(d / fn)
    s = run_tf_benchmark(m, gtf, peaks, out, args.windows, args.cutoff, args.upstream_compat, not args.no_plots)
    finish(out, "TF benchmark (ENCODE ChIP-seq targets)", s, [md_table(["TF", "window", "rejections", "targets", "overlap", "odds ratio", "p", "recall"],
                                                                     [[r["TF_display"], r["promoter_window_width"], r["num_rejections"], r["num_targets"], r["num_overlapping"],
                                                                       r["odds_ratio"], r["pvalue"], r["recall"]] for r in s["best_per_tf"]])])
    return 0


# ---------------------------------------------------------------------------
# Dashboard: process_json.py, create_dashboard_plots(.py/_HASHING.py), create_dashboard_df(.py/_HASHING.py), create_dashboard.py
# ---------------------------------------------------------------------------

def human_format(num) -> str:
    try:
        num = float(num)
    except (TypeError, ValueError):
        return str(num)
    mag = 0
    while abs(num) >= 1000 and mag < 5:
        mag += 1
        num /= 1000.0
    return f"{num:.1f}{['', 'K', 'M', 'G', 'T', 'P'][mag]}"


def collect_kb_json(dirs: "list[Path]") -> "list[dict]":
    """process_json.py: <x>_ks_transcripts_out -> trans-<x>, _ks_guide_out -> guide-<x>, _ks_hashing_out -> hashing-<x>."""
    out = []
    pats = [(re.compile(r"(.+)_ks_transcripts_out"), "trans", "scRNA"), (re.compile(r"(.+)_ks_guide_out"), "guide", "Guide"),
            (re.compile(r"(.+)_ks_hashing_out"), "hashing", "Hashing")]
    for d in dirs:
        for pat, pre, mod in pats:
            mm = pat.match(Path(d).name)
            if mm and (Path(d) / "inspect.json").is_file() and (Path(d) / "run_info.json").is_file():
                comb = {**json.loads((Path(d) / "inspect.json").read_text()), **json.loads((Path(d) / "run_info.json").read_text())}
                comb.pop("start_time", None)
                comb.pop("call", None)
                out.append({"prefix": f"{pre}-{mm.group(1)}", "modality": mod, "sample": mm.group(1), "table": comb})
    return out


def dashboard_plots(m: MData, figdir: Path, hashing_unfiltered=None) -> "list[str]":
    np, pd = _np(), _pd()
    plt = _plt()
    if not plt:
        return []
    figs = []
    X, G = as_csr(m["gene"].X), as_csr(m["guide"].X)

    def bars(cuts, vals, title, xlabel, name):
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.bar([str(c) for c in cuts], vals, color=BLUE)
        for i, v in enumerate(vals):
            ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
        _style(ax, title, xlabel, "Number of cells")
        figs.append(str(_save(fig, figdir / name)))
    tot = row_sums(X)
    c1 = [200, 500, 1000, 2000, 5000]
    bars(c1, [int((tot > c).sum()) for c in c1], "scRNA barcodes by total UMI threshold", "Total UMI threshold", "scRNA_barcodes_UMI_thresholds.png")
    c2 = [200, 500, 1000, 1500, 2000]
    ng = nnz_rows(X)
    bars(c2, [int((ng > c).sum()) for c in c2], "scRNA barcodes with detected genes > threshold", "Total detected genes threshold", "scRNA_barcodes_detected_genes_thresholds.png")
    gt = row_sums(G)
    c3 = [0, 5, 10, 50, 100, 1000]
    bars(c3, [int((gt > c).sum()) for c in c3], "Cells by total guide UMI threshold", "Guides Total UMI Threshold", "guides_UMI_thresholds.png")
    A = as_csr(m["guide"].layers["guide_assignment"])
    freq = pd.DataFrame({"sgRNA": list(map(str, m["guide"].var_names)), "Frequency": col_sums(A)}).sort_values("Frequency")
    write_tsv(freq, figdir / "sgRNA_frequencies.csv", sep=",", quiet=True)
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.bar(range(len(freq)), freq["Frequency"].values, color=ORANGE, width=1.0)
    _style(ax, "Cells per sgRNA (guide_assignment)", "sgRNA (sorted)", "cells")
    figs.append(str(_save(fig, figdir / "guides_hist_num_sgRNA.png")))
    for vals, title, name in ((row_sums(A), "Histogram of Guides per Cell", "guides_per_cell_histogram.png"),
                              (np.sort(col_sums(A)), "Histogram of Cells per Guide", "cells_per_guide_histogram.png")):
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.hist(vals, bins=50, color="#9cc9ee", ec="#4b7fb3")
        _style(ax, title, "", "")
        figs.append(str(_save(fig, figdir / name)))
    if hashing_unfiltered is not None and "hto_type_split" in hashing_unfiltered.obs.columns:
        cnt = hashing_unfiltered.obs["hto_type_split"].astype(str).value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.bar(cnt.index, cnt.values, color="#9cc9ee")
        ax.set_yscale("log")
        _style(ax, "Number of Cells per HTO Type", "HTO Type", "cells (log)")
        figs.append(str(_save(fig, figdir / "cells_per_hto_barplot.png")))
        H = hashing_unfiltered.X.toarray() if hasattr(hashing_unfiltered.X, "toarray") else np.asarray(hashing_unfiltered.X)
        clr = np.vstack([np.log1p(r / np.exp(np.log1p(r[r > 0]).sum() / len(r))) if (r > 0).any() else r for r in H])
        from sklearn.decomposition import PCA  # type: ignore
        from sklearn.preprocessing import StandardScaler  # type: ignore
        Z = StandardScaler().fit_transform(clr)
        emb = PCA(n_components=min(2, Z.shape[1]), random_state=42).fit_transform(Z)
        try:
            from umap import UMAP  # type: ignore
            emb = UMAP(n_neighbors=10, min_dist=0.1, random_state=42).fit_transform(PCA(n_components=min(10, Z.shape[1]), random_state=42).fit_transform(Z))
            lab = "UMAP"
        except Exception:
            lab = "PC"
        fig, ax = plt.subplots(figsize=(5, 4.5))
        types = sorted(hashing_unfiltered.obs["hto_type_split"].astype(str).unique())
        cm = plt.cm.tab20
        for k, t in enumerate(types):
            mk = hashing_unfiltered.obs["hto_type_split"].astype(str).values == t
            ax.scatter(emb[mk, 0], emb[mk, 1] if emb.shape[1] > 1 else np.zeros(mk.sum()), s=4, color=cm(k % 20), label=t)
        ax.legend(fontsize=6, frameon=False, markerscale=3)
        _style(ax, f"HTO CLR embedding ({lab})", f"{lab}1", f"{lab}2")
        figs.append(str(_save(fig, figdir / "umap_hto.png")))
    return figs


def flow_html(sample_counts, concat_count, filter_count, filtered_count, guide_intersection, final_count, qc: dict, hashing: Optional[dict] = None) -> str:
    def step(title, count=None, removed=None, note=None, subitems=None):
        p = [f"<b>{title}</b>"] + [f"<span>{s}</span>" for s in (subitems or [])]
        if count is not None:
            p.append(f"<span>Cells: {count:,}</span>")
        if removed is not None:
            p.append(f"<span>Removed: {removed:,}</span>")
        if note:
            p.append(f"<span class='desc'>{note}</span>")
        return "<div class='flow-step'>" + "".join(p) + "</div>"
    rem = lambda a, b: None if a is None or b is None else max(int(a) - int(b), 0)  # noqa: E731
    steps = []
    tot = sum(c for _, c in sample_counts) if sample_counts else None
    steps.append(step("Unfiltered scRNA barcodes (per sample)", subitems=[f"{n}: {c:,}" for n, c in sample_counts] + ([f"Sum: {tot:,}"] if tot else []),
                      note=None if sample_counts else "Per sample counts not found"))
    steps.append(step("Concatenated scRNA anndata (unfiltered)", concat_count, rem(tot, concat_count), "Merged across samples"))
    bf = qc.get("barcode_filter")
    steps.append(step(f"Barcode filter ({bf})" if bf in ("knee", "knee2") else f"Min genes per cell (>= {qc.get('min_genes')})", filter_count, rem(concat_count, filter_count)))
    steps.append(step(f"Mito filter (pct_counts_mt < {qc.get('pct_mito')})", filtered_count, rem(filter_count, filtered_count)))
    if hashing:
        steps.append(step("Intersection with hashing barcodes", hashing.get("rna_hashing"), rem(filtered_count, hashing.get("rna_hashing"))))
        steps.append(step("Demultiplex filter (HTO)", hashing.get("hashing_demux"), rem(hashing.get("rna_hashing"), hashing.get("hashing_demux")), "Remove negative and multiplet HTOs"))
        steps.append(step("Intersection with guide barcodes", final_count, rem(hashing.get("hashing_demux"), final_count)))
    else:
        steps.append(step("Intersection with guide barcodes", guide_intersection, rem(filtered_count, guide_intersection), "Filtered RNA intersect guide"))
        steps.append(step("Doublet removal (Scrublet)" if qc.get("enable_scrublet") else "Final cells in MuData", final_count, rem(guide_intersection, final_count)))
    return "<div class='flow'>" + "<div class='flow-arrow'>&darr;</div>".join(steps) + "</div>"


def build_dashboard(m: MData, out: Path, blocks_extra: "list[tuple[str, str, list]]", kb_json: "list[dict]", flow: Optional[str],
                    image_dirs: "list[Path]", default_mode: bool = True, title: str = "CRISPR Pipeline Dashboard") -> Path:
    pd, np = _pd(), _np()
    tabs: "dict[str, list]" = defaultdict(list)
    if flow:
        tabs["Filtering Summary"].append(("Barcode filtering flow", flow))
    g = m["gene"]
    tabs["Filtering Summary"].append(("Gene statistics", f"<p>Number of genes detected after filtering: {human_format(g.n_vars)}, mean UMI counts per cell after filtering: "
                                                         f"{human_format(float(row_sums(g.X).mean()))}; cells in MuData: {human_format(m.n_obs)}</p>"))
    for j in sorted(kb_json, key=lambda x: x["prefix"]):
        t = j["table"]
        disp = (f"Total Reads for {j['modality']}: {human_format(t.get('n_processed'))}, Paired Reads Mapped: {human_format(t.get('n_pseudoaligned'))}, "
                f"Alignment Percentage: {_fmt(t.get('p_pseudoaligned'), 1)}%, Total Detected {j['modality']} Barcodes (Unfiltered): {human_format(t.get('numBarcodes'))}")
        tabs[j["modality"]].append((f"Mapping {j['modality']} {j['sample']}", f"<p>{disp}</p>" + table_html(pd.DataFrame({"parameter": list(t), "value": [str(v) for v in t.values()]}), 60)))
    A = as_csr(m["guide"].layers["guide_assignment"])
    freq = pd.DataFrame({"sgRNA": list(map(str, m["guide"].var_names)), "Assignment sum": col_sums(A)})
    tabs["Guide"].append(("Guide assignment", f"<p>Total sgRNA assignment values across all guides: {human_format(freq['Assignment sum'].sum())}, "
                                              f"median per sgRNA: {human_format(freq['Assignment sum'].median())}</p>" + table_html(freq, 40)))
    tabs["Inference"].append(("Guides / cells", f"<p>Mean guides per cell: {human_format(float(row_sums(A).mean()))}, mean cells per guide: {human_format(float(col_sums(A).mean()))}</p>"))
    keys = [("cis_per_guide_results", "Cis"), ("trans_per_guide_results", "Trans")] if default_mode else [("per_guide_results", ""), ("cis_per_guide_results", "")]
    for k, lab in keys:
        r = uns_df(m, k)
        if r is None:
            continue
        pc = next((c for c in ("sceptre_p_value", "p_value", "perturbo_p_value") if c in r.columns), None)
        top = r.sort_values(pc).head(10000) if pc else r
        tabs["Inference"].append((f"{lab} Analysis".strip(), f"<p>Top lowest p-values {len(top)} tested sgRNA-gene pairs ({k})</p>" + table_html(top, 200)))
    for d in image_dirs:
        if not d or not Path(d).is_dir():
            continue
        for png in sorted(Path(d).rglob("*.png")):
            rel = str(png.relative_to(d))
            tab = ("QC" if "additional_qc" in str(d) or any(x in rel for x in ("gene/", "guide/", "intended_target", "trans/")) else
                   "Evaluation" if "evaluation" in str(d) else "Benchmark" if "benchmark" in str(d) else
                   "Guide" if "guide" in png.name.lower() or "sgRNA" in png.name else "Hashing" if "hto" in png.name.lower() else "scRNA")
            tabs[tab].append((png.stem, img_tag(png, rel)))
    for tab, desc, items in blocks_extra:
        tabs[tab].append((desc, "".join(items)))
    order = [t for t in ["Filtering Summary", "scRNA", "Guide", "Hashing", "Inference", "QC", "Evaluation", "Benchmark"] if t in tabs] + [t for t in tabs if t not in
                                                                                                                                    ["Filtering Summary", "scRNA", "Guide", "Hashing", "Inference", "QC", "Evaluation", "Benchmark"]]
    nav = "<div class='tabs'>" + "".join(f"<a href='#{safe_label(t)}'>{t}</a>" for t in order) + "</div>"
    body = "".join(f"<h2 id='{safe_label(t)}'>{t}</h2>" + "".join(f"<h3>{h}</h3>{c}" for h, c in tabs[t]) for t in order)
    html = (f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>CRISPR Dashboard</title>"
            f"{QC_REPORT_CSS}</head><body><h1>{title}</h1><div class='desc'>igvfagent crispr-pipeline (port of {UPSTREAM_REPO}) · engine: "
            f"{', '.join(map(str, m.uns.get('inference_engine', []))) or 'n/a'} · {time.strftime('%Y-%m-%d %H:%M')}</div>{nav}{body}</body></html>")
    p = out / "dashboard.html"
    p.write_text(html)
    print(f"Wrote: {p}")
    return p


def cmd_dashboard(args) -> int:
    out = run_dir(args.label)
    m = read_h5mu(Path(args.mudata))
    hu = read_h5ad(Path(args.hashing_unfiltered_demux)) if args.hashing_unfiltered_demux else None
    figs = dashboard_plots(m, out / "figures", hu) if not args.no_plots else []
    kbj = collect_kb_json([Path(x) for x in (args.ks_dirs or [])])
    flow = None
    if args.gene_ann and args.gene_ann_filtered and args.guide_ann:
        ga, gf, gd = read_h5ad(Path(args.gene_ann)), read_h5ad(Path(args.gene_ann_filtered)), read_h5ad(Path(args.guide_ann))
        thr, _ = barcode_filter_threshold(row_sums(ga.X), args.barcode_filter) if args.barcode_filter in ("knee", "knee2") else (None, None)
        fc = int((row_sums(ga.X) >= thr).sum()) if thr is not None else int((nnz_rows(ga.X) >= args.min_genes).sum())
        sc = [(j["sample"], int(j["table"].get("numBarcodes", 0))) for j in kbj if j["prefix"].startswith("trans-")]
        flow = flow_html(sc, int(ga.n_obs), fc, int(gf.n_obs), len(set(gf.obs_names) & set(gd.obs_names)), m.n_obs,
                         {"barcode_filter": args.barcode_filter, "min_genes": args.min_genes, "pct_mito": args.pct_mito})
    p = build_dashboard(m, out, [], kbj, flow, [out / "figures"] + [Path(x) for x in (args.image_dirs or [])], not args.not_default)
    finish(out, "Dashboard", {"html": str(p), "figures": figs}, [f"Dashboard: `{p}`"])
    return 0


# ---------------------------------------------------------------------------
# run: count matrices -> preprocess -> MuData -> assignment -> pairs -> inference -> QC -> evaluation -> dashboard
# ---------------------------------------------------------------------------

def _load_counts(paths: "list[str]", covariates=None):
    if len(paths) == 1 and Path(paths[0]).is_file() and paths[0].endswith(".h5ad"):
        return read_h5ad(Path(paths[0])), None
    a = concat_ks_dirs(paths, covariates)
    first = sorted(paths, key=lambda x: os.path.basename(os.path.normpath(x)))[0]
    gn = Path(first) / "counts_unfiltered" / "cells_x_genes.genes.names.txt"
    names = [l.strip() for l in open(gn) if l.strip()] if gn.is_file() else None
    return a, names


def run_pipeline(args, out: Path) -> dict:
    pd = _pd()
    plots = not args.no_plots
    stages: "dict[str, Any]" = {}
    t0 = time.time()
    rna, names = _load_counts(args.rna)
    guide, _ = _load_counts(args.guide)
    write_h5ad(rna, out / "concatenated_adata_rna.h5ad", quiet=True)
    stages["concat"] = {"rna_barcodes": int(rna.n_obs), "guide_barcodes": int(guide.n_obs)}
    filt, _, pinfo = preprocess_rna(rna, names, args.min_genes, args.pct_mito, args.reference, args.barcode_filter, out / "figures",
                                    args.upstream_compat, plots)
    write_h5ad(filt, out / "filtered_anndata.h5ad", quiet=True)
    stages["preprocess"] = pinfo
    gtf = read_gtf(Path(args.gtf)) if args.gtf else None
    hashing_f = hashing_u = None
    if args.hashing:
        hraw, _ = _load_counts(args.hashing)
        hashing_f, hashing_u, hinfo = run_hashing(hraw, filt, out / "hashing")
        write_h5ad(hashing_f, out / "concatenated_hashing_demux.h5ad", quiet=True)
        write_h5ad(hashing_u, out / "concatenated_unfiltered_hashing_demux.h5ad", quiet=True)
        stages["hashing"] = hinfo
    meta = read_guide_metadata(Path(args.guide_metadata))
    m = create_mudata(filt, guide, meta, gtf, args.moi, args.capture_method, hashing_f, out / "figures", plots)
    write_h5mu(m, out / "mudata.h5mu", quiet=True)
    stages["create_mudata"] = {"n_cells": m.n_obs, "n_genes": int(m["gene"].n_vars), "n_guides": int(m["guide"].n_vars)}
    if args.scrublet and not args.hashing:
        m, dinfo = remove_doublets(m, out / "figures", plots)
        stages["doublets"] = dinfo
    m, ainfo = assign_guides(m, args.assignment_method, args.assignment_threshold, args.n_em_rep, args.umi_threshold,
                             capture_method=args.capture_method, engine=args.engine, upstream_bin=args.upstream_bin, workdir=out / "assign_work")
    m = filter_genes_by_cells(m, args.min_cells_fraction)
    ainfo["genes_after_filter"] = int(m["gene"].n_vars)
    if args.dual_guide:
        m, ainfo["collapse"] = collapse_guides(m)
    write_h5mu(m, out / "concat_mudata.h5mu", quiet=True)
    stages["guide_assignment"] = ainfo
    strat = args.pairing_strategy
    if strat == "predefined_pairs":
        pairs = pd.read_csv(args.pairs)
    elif strat in ("default", "by_distance"):
        if gtf is None:
            raise SystemExit("--gtf is required for the default / by_distance pairing strategy")
        pairs = create_pairs_to_test(m["guide"].var, m["gene"].var_names, gtf, args.max_target_distance)
    else:
        pairs = None
    if pairs is not None:
        write_tsv(pairs, out / "pairs_to_test.csv", sep=",", quiet=True)
        m_cis, pi = prepare_inference(m, pairs, subset_for_cis=(strat == "default"))
        stages["pairs"] = dict(pi, strategy=strat)
    else:
        m_cis = m
        stages["pairs"] = {"strategy": strat}
    write_h5mu(m_cis, out / "mudata_inference_input.h5mu", quiet=True)
    im, tables, iinfo = run_inference(m_cis, args.inference_method, m if args.inference_method == "default" else None, args.engine,
                                      args.n_resamples, args.side, out / "inference_work", args.upstream_bin, args.seed)
    ip = write_h5mu(im, out / "inference_mudata.h5mu")
    for k, v in tables.items():
        write_tsv(v, out / f"{k}.tsv.gz", quiet=True)
    pgk = "cis_per_guide_results" if "cis_per_guide_results" in tables else "per_guide_results"
    pg, pe = export_outputs(im, tables[pgk])
    write_tsv(pg, out / "per_guide_output.tsv", quiet=True)
    write_tsv(pe, out / "per_element_output.tsv", quiet=True)
    stages["inference"] = iinfo
    stages["qc"] = run_qc(im, out / "additional_qc", plots=plots, sample=args.label)
    stages["evaluate"] = run_evaluate(im, out, gtf, plots, args.num_nodes)
    if args.encode_bed_dir or args.peaks:
        peaks = dict(s.partition("=")[::2] for s in args.peaks) if args.peaks else {
            k: str(Path(args.encode_bed_dir) / v) for k, v in PEAK_META_DICT.items() if (Path(args.encode_bed_dir) / v).is_file()}
        if gtf is not None and peaks:
            try:
                stages["tf_benchmark"] = run_tf_benchmark(im, gtf, peaks, out / "benchmark_output", compat=args.upstream_compat, plots=plots)
            except SystemExit as e:
                stages["tf_benchmark"] = {"error": str(e)}
    figs = dashboard_plots(im, out / "figures", hashing_u) if plots else []
    kbj = collect_kb_json([Path(x) for x in args.rna + args.guide + (args.hashing or []) if Path(x).is_dir()])
    sc = [(j["sample"], int(j["table"].get("numBarcodes", 0))) for j in kbj if j["prefix"].startswith("trans-")]
    fc = pinfo["n_after_barcode_filter"] if args.barcode_filter != "none" else pinfo["n_after_min_genes"]
    hc = {"rna_hashing": stages["hashing"]["n_intersecting"], "hashing_demux": stages["hashing"]["n_singlets"]} if args.hashing else None
    flow = flow_html(sc, pinfo["n_barcodes_in"], fc, pinfo["n_after_mito"], len(set(filt.obs_names) & set(guide.obs_names)), im.n_obs,
                     {"barcode_filter": args.barcode_filter, "min_genes": args.min_genes, "pct_mito": args.pct_mito, "enable_scrublet": args.scrublet}, hc)
    dash = build_dashboard(im, out, [], kbj, flow, [out / "figures", out / "additional_qc", out / "evaluation_output", out / "benchmark_output"],
                           args.inference_method == "default")
    stages["dashboard"] = {"html": str(dash), "figures": len(figs)}
    stages["seconds"] = round(time.time() - t0, 1)
    stages["inference_mudata"] = str(ip)
    return stages


def cmd_run(args) -> int:
    out = run_dir(args.label)
    s = run_pipeline(args, out)
    it = s["qc"].get("intended_target", {})
    rows = [["barcodes in", s["preprocess"]["n_barcodes_in"]], ["after barcode filter", s["preprocess"]["n_after_barcode_filter"]],
            ["after mito filter", s["preprocess"]["n_after_mito"]], ["cells in MuData", s["create_mudata"]["n_cells"]],
            ["cells with a guide", s["guide_assignment"]["cells_with_guide"]], ["genes tested", s["guide_assignment"]["genes_after_filter"]]]
    secs = [md_table(["stage", "cells / genes"], rows),
            md_table(["uns table", "rows"], [[k, v] for k, v in s["inference"]["tables"].items()]),
            f"Engines: {', '.join(s['inference']['engines'])}. Intended-target AUROC {_fmt(it.get('auroc'))}, AUPRC {_fmt(it.get('auprc'))}; "
            f"evaluate_controls AUROC {_fmt(s['evaluate']['controls']['auroc'])}.",
            f"Dashboard: `{s['dashboard']['html']}`; QC report: `{s['qc']['html']}`; MuData: `{s['inference_mudata']}`"]
    finish(out, f"CRISPR pipeline run: {args.label}", s, secs)
    return 0


def cmd_nextflow(args) -> int:
    out = run_dir(args.label)
    cmd = ["nextflow", "run", args.pipeline, "-r", args.revision, "-profile", args.profile, "--input", args.input, "--outdir", args.outdir]
    for kv in args.param or []:
        k, _, v = kv.partition("=")
        cmd += [f"--{k}", v]
    if args.resume:
        cmd.append("-resume")
    ran, text = external(cmd, "nextflow", run=args.execute)
    (out / "nextflow_command.sh").write_text("#!/bin/bash\n" + text + "\n")
    print(f"Wrote: {out / 'nextflow_command.sh'}")
    finish(out, "Upstream Nextflow command", {"command": text, "ran": ran, "defaults": DEFAULTS},
           [f"```bash\n{text}\n```", "Defaults (nextflow.config):\n\n" + md_table(["param", "value"], sorted(DEFAULTS.items()))])
    return 0


# ---------------------------------------------------------------------------
# selftest: synthetic Perturb-seq screen with planted effects
# ---------------------------------------------------------------------------

GUIDE_SEQSPEC_YAML = """!Assay
seqspec_version: 0.2.0
assay_id: MULTISEQ_10XV3_Guide
modalities:
- guide
sequence_spec:
- !Read
  read_id: guide_R1.fq
  name: Read 1
  modality: guide
  primer_id: r1_primer
  min_len: 28
  max_len: 28
  strand: pos
- !Read
  read_id: guide_R2.fq
  modality: guide
  primer_id: r2_primer
  min_len: 83
  max_len: 83
  strand: neg
library_spec:
- !Region
  parent_id: null
  region_id: guide
  region_type: null
  min_len: 111
  max_len: 111
  onlist: null
  regions:
  - !Region
    region_id: r1_primer
    region_type: r1_primer
    sequence_type: fixed
    min_len: 0
    max_len: 0
    regions: null
  - !Region
    region_id: barcode
    region_type: barcode
    sequence: NNNNNNNNNNNNNNNN
    min_len: 16
    max_len: 16
    onlist: !Onlist
      location: remote
      filename: 737K-august-2016.txt
      md5: f62a276e262fdd85262a889d0f48556b
    regions: null
  - !Region
    region_id: umi
    region_type: umi
    min_len: 12
    max_len: 12
    regions: null
  - !Region
    region_id: cdna
    region_type: cdna
    min_len: 20
    max_len: 20
    regions: null
  - !Region
    region_id: common
    region_type: common
    sequence: 'CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC'
    min_len: 63
    max_len: 63
    regions: null
  - !Region
    region_id: r2_primer
    region_type: r2_primer
    min_len: 0
    max_len: 0
    regions: null
"""


def synthetic_screen(d: Path, seed: int = 3) -> dict:
    """Two batches of cells + an ambient barcode tail; 6 targets x 2 guides + 8 NT guides; knockdown x0.2 on the target gene,
    T1 represses genes 100-114 (x0.4); 5% high-mito cells; kb-style output dirs, GTF and guide metadata."""
    np, pd, sp = _np(), _pd(), _sp()
    import anndata as ad  # type: ignore
    rng = np.random.default_rng(seed)
    n_genes = 123
    gids = [f"ENSG{100000 + i:011d}" for i in range(n_genes)]
    syms = [f"GENE{i}" for i in range(n_genes - 3)] + ["MT-CO1", "MT-ND1", "MT-ATP6"]
    starts = [1_000_000 + i * 200_000 for i in range(n_genes)]
    strands = ["+" if i % 2 == 0 else "-" for i in range(n_genes)]
    target_idx = [10, 25, 40, 55, 70, 85]
    trans_idx = list(range(100, 115))
    guides, gmeta = [], []
    letters = "ACGT"
    for t, gi in enumerate(target_idx):
        for k in range(2):
            gid = f"{syms[gi]}_sg{k + 1}"
            guides.append(gid)
            gmeta.append({"guide_id": gid, "spacer": "".join(rng.choice(list(letters), 20)), "targeting": True, "type": "targeting", "guide_chr": "chr1",
                          "guide_start": float(starts[gi] - 300 + 40 * k), "guide_end": float(starts[gi] - 280 + 40 * k), "strand": "+", "pam": "NGG",
                          "intended_target_name": gids[gi], "intended_target_chr": "chr1", "intended_target_start": float(starts[gi]),
                          "intended_target_end": float(starts[gi] + 20_000), "gene_name": syms[gi]})
    for k in range(8):
        gid = f"NTC_{k + 1}"
        guides.append(gid)
        gmeta.append({"guide_id": gid, "spacer": "".join(rng.choice(list(letters), 20)), "targeting": False, "type": "non-targeting", "guide_chr": np.nan,
                      "guide_start": np.nan, "guide_end": np.nan, "strand": np.nan, "pam": "NGG", "intended_target_name": "non-targeting",
                      "intended_target_chr": np.nan, "intended_target_start": np.nan, "intended_target_end": np.nan, "gene_name": np.nan})
    meta = pd.DataFrame(gmeta)
    mp = d / "guide_metadata.tsv"
    meta.to_csv(mp, sep="\t", index=False)
    gtf = d / "genes.gtf"
    with open(gtf, "w") as fh:
        fh.write("##description: synthetic\n")
        for i in range(n_genes):
            attrs = f'gene_id "{gids[i]}.3"; gene_type "protein_coding"; gene_name "{syms[i]}";'
            fh.write(f"chr1\tSYN\tgene\t{starts[i]}\t{starts[i] + 20000}\t.\t{strands[i]}\t.\t{attrs}\n")
            fh.write(f"chr1\tSYN\ttranscript\t{starts[i]}\t{starts[i] + 20000}\t.\t{strands[i]}\t.\t{attrs} transcript_id \"T{i}\";\n")
    base = np.exp(rng.normal(0.3, 0.9, n_genes))
    base[-3:] = [3.0, 2.0, 2.0]
    base[115:120] = 0.001  # rare genes, removed by the > int(0.05 n) cells filter
    truth = {}
    ks_rna, ks_guide = [], []
    for bi, b in enumerate(["B1", "B2"]):
        n_cells, n_empty = 700, 1500
        cells = ["".join(rng.choice(list(letters), 16)) for _ in range(n_cells + n_empty)]
        n_g = np.where(rng.random(n_cells) < 0.8, 1, 2)
        A = np.zeros((n_cells + n_empty, len(guides)), bool)
        for c in range(n_cells):
            A[c, rng.choice(len(guides), n_g[c], replace=False)] = True
        sf = np.r_[np.exp(rng.normal(0, 0.3, n_cells)), np.exp(rng.normal(np.log(0.004), 1.0, n_empty))]
        mu = sf[:, None] * base[None, :] * 18
        for t, gi in enumerate(target_idx):
            hit = A[:, 2 * t] | A[:, 2 * t + 1]
            mu[hit, gi] *= 0.2
            if t == 0:
                mu[np.ix_(hit, trans_idx)] *= 0.4
        himito = np.zeros(n_cells + n_empty, bool)
        himito[rng.choice(n_cells, int(0.05 * n_cells), replace=False)] = True
        mu[himito, -3:] *= 60
        theta = 6.0
        Y = rng.poisson(rng.gamma(theta, mu / theta))
        G = rng.poisson(0.08 * np.r_[np.ones(n_cells), 0.3 * np.ones(n_empty)][:, None] * np.ones(len(guides))[None, :])
        G = G + np.where(A, rng.negative_binomial(4, 4 / (4 + 40), A.shape), 0)
        rows = list(range(n_cells + n_empty))
        rng.shuffle(rows)
        Y, G, A, himito = Y[rows], G[rows], A[rows], himito[rows]
        cells = [cells[i] for i in rows]
        is_cell = np.array([i < n_cells for i in rows])
        rd = d / f"{b}_ks_transcripts_out"
        a = ad.AnnData(sp.csr_matrix(Y.astype(np.float32)), obs=pd.DataFrame(index=pd.Index(cells, name="barcode")),
                       var=pd.DataFrame(index=pd.Index([f"{g}.3" for g in gids], name="gene_id")))
        write_kb_like(a, {"numRecords": int(Y.sum() * 1.3), "numReads": int(Y.sum() * 1.3), "numBarcodes": len(cells) + 500, "numUMIs": int(Y.sum())},
                      {"n_processed": int(Y.sum() * 1.4), "n_pseudoaligned": int(Y.sum() * 1.3), "p_pseudoaligned": 92.9, "n_unique": int(Y.sum())}, rd)
        (rd / "counts_unfiltered" / "cells_x_genes.genes.names.txt").write_text("\n".join(syms) + "\n")
        gd = d / f"{b}_ks_guide_out"
        keep = G.sum(axis=1) > 0
        ag = ad.AnnData(sp.csr_matrix(G[keep].astype(np.float32)), obs=pd.DataFrame(index=pd.Index([c for c, k in zip(cells, keep) if k], name="barcode")),
                        var=pd.DataFrame(index=pd.Index(guides, name="gene_id")))
        write_kb_like(ag, {"numRecords": int(G.sum() * 2), "numReads": int(G.sum() * 2), "numBarcodes": int(keep.sum()), "numUMIs": int(G.sum())},
                      {"n_processed": int(G.sum() * 2), "n_pseudoaligned": int(G.sum() * 1.5), "p_pseudoaligned": 75.0, "n_unique": int(G.sum())}, gd)
        ks_rna.append(str(rd))
        ks_guide.append(str(gd))
        for c, cellflag, arow, hm in zip(cells, is_cell, A, himito):
            truth[f"{c}_{bi}"] = (cellflag, set(np.array(guides)[arow]), hm)
    return {"rna": ks_rna, "guide": ks_guide, "gtf": gtf, "meta": mp, "gids": gids, "syms": syms, "target_idx": target_idx, "trans_idx": trans_idx,
            "guides": guides, "truth": truth, "starts": starts}


def synthetic_fastqs(d: Path, seed: int = 5) -> dict:
    np, pd = _np(), _pd()
    rng = np.random.default_rng(seed)
    L = "ACGT"
    onlist = ["".join(rng.choice(list(L), 16)) for _ in range(200)]
    (d / "onlist.txt").write_text("\n".join(onlist) + "\n")
    feats = {f"g{k}": "".join(rng.choice(list(L), 20)) for k in range(6)}
    pd.DataFrame({"guide_id": list(feats), "spacer": list(feats.values())}).to_csv(d / "fq_guides.tsv", sep="\t", index=False)
    flank = "TTGTGGAAAGGACGAAACACCG"
    r1, r2 = [], []
    truth: Counter = Counter()
    for c in range(30):
        bc = onlist[c]
        for g in rng.choice(list(feats), 2, replace=False):
            for u in range(int(rng.integers(3, 8))):
                umi = "".join(rng.choice(list(L), 12))
                truth[(bc, g)] += 1
                for rep in range(int(rng.integers(1, 4))):
                    sp_ = feats[g]
                    if rep == 1:
                        i = int(rng.integers(0, 20))
                        sp_ = sp_[:i] + ("A" if sp_[i] != "A" else "C") + sp_[i + 1:]
                    bc_ = bc
                    if rep == 2:
                        bc_ = ("A" if bc[0] != "A" else "C") + bc[1:]
                    r1.append(bc_ + umi)
                    r2.append("".join(rng.choice(list(L), 5)) + flank + sp_ + "GTTTTAGAGCTAGAAATAGC")
    for _ in range(40):
        r1.append("".join(rng.choice(list(L), 28)))
        r2.append("".join(rng.choice(list(L), 67)))

    def wfq(path, seqs):
        with gzip.open(path, "wt") as fh:
            for i, s in enumerate(seqs):
                fh.write(f"@r{i}\n{s}\n+\n{'F' * len(s)}\n")
    wfq(d / "g_R1.fastq.gz", r1)
    wfq(d / "g_R2.fastq.gz", r2)
    return {"r1": d / "g_R1.fastq.gz", "r2": d / "g_R2.fastq.gz", "onlist": d / "onlist.txt", "features": feats, "truth": truth, "flank": flank,
            "meta": d / "fq_guides.tsv"}


def _ns(**kw):
    return argparse.Namespace(**kw)


def cmd_selftest(args) -> int:
    checks: "list[tuple[bool, str]]" = []
    made: "list[Path]" = []
    t0 = time.time()

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def latest(label):
        pat = re.compile(rf"^\d{{8}}_\d{{6}}_{re.escape(label)}(_\d+)?$")
        p = sorted((x for x in OUT_ROOT.iterdir() if pat.match(x.name)), key=lambda x: x.stat().st_mtime)[-1]
        made.append(p)
        return p

    before = set(OUT_ROOT.iterdir()) if OUT_ROOT.is_dir() else set()
    try:
        _selftest_body(args, check, latest, t0)
    except Exception as e:  # report, clean up, fail
        import traceback
        traceback.print_exc()
        check(False, f"selftest raised {type(e).__name__}: {e}")
    finally:
        new = (set(OUT_ROOT.iterdir()) - before) if OUT_ROOT.is_dir() else set()
        for p in set(made) | {x for x in new if "_st_" in x.name}:
            shutil.rmtree(p, ignore_errors=True)
        if OUT_ROOT.is_dir() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print(f"\n{len(checks)} checks in {time.time() - t0:.0f} s")
    print("selftest: all checks pass" if ok else f"selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


def _selftest_body(args, check, latest, t0) -> None:
    import tempfile
    np, pd, sp = _np(), _pd(), _sp()
    import anndata as ad  # type: ignore
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        print("\nseqspec + samplesheet + portal")
        yp = d / "guide_seqspec.yml"
        yp.write_text(GUIDE_SEQSPEC_YAML)
        spec = mini_yaml(GUIDE_SEQSPEC_YAML)
        check(spec["sequence_spec"][1]["strand"] == "neg" and spec["library_spec"][0]["regions"][1]["onlist"]["md5"].startswith("f62a")
              and spec["modalities"] == ["guide"], "mini_yaml: !Tags, nested lists of mappings, inline-mapping list items, quoted scalars")
        rep, _ = seqspec_index_kb(load_seqspec(yp), "guide")
        check(rep == "0,0,16:0,16,28:1,63,83", f"seqspec index -t kb: R1 barcode/UMI, R2 (neg strand) feature after the 63-bp common region -> {rep}")
        check(spacer_chemistry(rep) == "0,0,16:0,16,28:1,0,0", "spacer tag: feature triplet becomes '1,0,0' (mappingGuide awk rule)")
        rc = cmd_seqspec(_ns(yaml=str(yp), modalities=["guide"], whitelist="737K.txt", label="st_seqspec"))
        o = latest("st_seqspec")
        ps = pd.read_csv(o / "guide_parsed_seqSpec.txt", sep="\t")
        check(rc == 0 and list(ps.columns) == ["modality", "representation", "barcode_whitelist", "seqspec_file"] and "".join(ps["representation"]) == rep,
              "seqspec: <modality>_parsed_seqSpec.txt; extract_parsed_seqspec = ''.join(representation)")
        ss = pd.DataFrame({"R1_path": ["r1a", "r1b", "r1c"], "R2_path": ["r2a", "r2b", "r2c"], "file_modality": ["scRNA", "gRNA", "gRNA"],
                           "measurement_sets": ["M1", "M1", "M1"], "sequencing_run": ["1", "1", "2"], "lane": ["1", "1", "1"], "seqspec": ["rna.yml", "g.yml", "g.yml"],
                           "barcode_onlist": ["737K.txt", "", ""], "guide_design": ["", "gd.tsv", "gd.tsv"], "barcode_hashtag_map": ["", "", ""]})
        info = parse_samplesheet(ss)
        check(info["modalities"]["grna"] == ["sample_gRNA_M1"] and info["groups"]["sample_gRNA_M1"]["fastqs"] == ["r1b", "r2b", "r1c", "r2c"]
              and info["formula"] == "" and info["guide_design"] == "gd.tsv" and not info["hashing"],
              "samplesheet: runs grouped per modality x measurement set, covariate formula empty for one batch")
        pt = pd.DataFrame({"R1_path": ["IGVFFI0001AAAA"], "R1_md5sum": ["x"], "R2_path": ["IGVFFI0002AAAA"], "file_modality": ["gRNA sequencing"], "measurement_sets": ["M"],
                           "sequencing_run": [1], "lane": [1], "flowcell_id": ["F"], "seqspec": ["IGVFFI0003AAAA"], "barcode_onlist": ["IGVFFI0004AAAA"],
                           "guide_design": ["IGVFFI0005AAAA"], "barcode_hashtag_map": [""]})
        up = portal_to_samplesheet(pt)
        plan = portal_download_plan(pt, "all", d / "dl")
        check(up["file_modality"].tolist() == ["gRNA"] and {p["url"].split("/")[3] for p in plan} == {"sequence-files", "configuration-files", "tabular-files"}
              and len(plan) == 5 and plan[0]["md5"] == "x", "update_samplesheet modality map; download plan routes + md5 for FASTQs")
        fake = {
            "/analysis-sets/IGVFDS1/@@object?format=json": {"construct_library_sets": ["/construct-library-sets/C1/"],
                                                            "input_file_sets": ["/measurement-sets/MS1/", "/auxiliary-sets/AX1/"]},
            "/construct-library-sets/C1/@@embedded?format=json": {"integrated_content_files": [{"accession": "IGVFFIGUIDE", "content_type": "guide RNA sequences", "status": "released"}]},
            "/measurement-sets/MS1/@@object?format=json": {"@id": "/measurement-sets/MS1/", "accession": "MS1", "file_set_type": "experimental data",
                                                           "files": ["/sequence-files/F1/", "/sequence-files/F2/"], "onlist_files": ["/tabular-files/ONL/"],
                                                           "onlist_method": "no combination", "strand_specificity": "5 prime to 3 prime"},
            "/auxiliary-sets/AX1/@@object?format=json": {"@id": "/auxiliary-sets/AX1/", "accession": "AX1", "file_set_type": "gRNA sequencing",
                                                         "measurement_sets": ["/measurement-sets/MS1/"], "files": ["/sequence-files/F3/", "/sequence-files/F4/"]},
        }
        for f, rt, ss_ in (("F1", "R1", "S1"), ("F2", "R2", None), ("F3", "R1", "S2"), ("F4", "R2", None)):
            fake[f"/sequence-files/{f}/@@object?format=json"] = {"@id": f"/sequence-files/{f}/", "content_type": "reads", "illumina_read_type": rt, "status": "released",
                                                                 "sequencing_run": 1, "lane": 1, "flowcell_id": "FC", "index": "i", "content_md5sum": f"md5{f}",
                                                                 "seqspecs": [f"/configuration-files/{ss_}/"] if ss_ else []}
        fake["/configuration-files/S1/@@object?format=json"] = {"upload_status": "validated", "status": "released"}
        fake["/configuration-files/S2/@@object?format=json"] = {"upload_status": "invalidated", "status": "released"}
        rows = generate_per_sample("IGVFDS1", lambda p: fake[p], sgrna_seqspec="fallback_sgrna.yaml")
        check(len(rows) == 2 and {r["file_modality"] for r in rows} == {"scRNA sequencing", "gRNA sequencing"} and rows[1]["seqspec"] == "fallback_sgrna.yaml"
              and rows[0]["seqspec"] == "S1" and rows[1]["measurement_sets"] == "MS1" and rows[0]["barcode_onlist"] == "ONL" and rows[0]["guide_design"] == "IGVFFIGUIDE",
              "portal-samplesheet: validated seqspec used, invalid one replaced by the fallback, aux set inherits the onlist")
        cfg = interface_config(pd.DataFrame({"tab_name": ["scRNA", "Guides", "params"], "batch_name": ["b1, donorA", "b1, donorA", np.nan],
                                             "read1": ["r1.fq", "g1.fq", np.nan], "read2": ["r2.fq", "g2.fq", np.nan],
                                             "variable": ["sequence", "sequence", "QC_pct_mito"], "variable_value": ["", "", "15.0"]}))
        check("QC_pct_mito = 15.0" in cfg and "distance_from_center = 1000000" in cfg and 'fastq_files_guide = [\n    "g1.fq g2.fq"' in cfg and "cov1: ['donorA']" in cfg,
              "interface-config: numeric params, default distance_from_center, FASTQ groups, covariate list")

        print("\nseqSpecCheck + feature reference + Python feature counting")
        fq = synthetic_fastqs(d)
        stats, best = seqspec_check([str(fq["r1"])], [str(fq["r2"])], list(fq["features"].values()), 100000)
        check(best[0]["Best_Config"] == "R2_Fwd" and best[0]["ModePosition"] == 27 and best[0]["TopFlank"] == fq["flank"][-12:],
              f"seqSpecCheck: guides found in R2 forward at position {best[0]['ModePosition']} with flank {best[0]['TopFlank']}")
        g_rev = [s for s in stats if s["Config"] == "R2_Rev"][0]
        check(g_rev["TotalHits"] == 0 and 0 <= [s for s in stats if s["IsWinner"]][0]["Gini"] < 0.6, "seqSpecCheck: reverse orientation has no hits; winner Gini < 0.6")
        mm = kite_mismatch_map({"a": "AAAA", "b": "AAAT"})
        check(mm["AAAA"] == "a" and mm["AAAT"] == "b" and "AAAC" not in mm and mm["CAAA"] == "a", "kite mismatch map: exact kept, colliding 1-mismatch variants dropped")
        a1, insp, ri = count_features([str(fq["r1"]), str(fq["r2"])], "0,0,16:0,16,28:1,27,47", fq["features"], open(fq["onlist"]).read().split())
        got = {(b, f): int(a1[b, f].X.sum()) for (b, f) in fq["truth"]}
        check(got == dict(fq["truth"]) and int(a1.X.sum()) == sum(fq["truth"].values()),
              f"map --engine python: {sum(fq['truth'].values())} planted UMIs recovered exactly (1-mismatch spacer and barcode reads collapse onto the same UMI)")
        a2, _, _ = count_features([str(fq["r1"]), str(fq["r2"])], spacer_chemistry("0,0,16:0,16,28:1,27,47"), fq["features"], open(fq["onlist"]).read().split(),
                                  spacer=fq["flank"][-12:])
        check(int(a2.X.sum()) == int(a1.X.sum()) and ri["p_pseudoaligned"] > 70, "map --engine python with a spacer tag (whole-read search) gives the same counts")
        rc = cmd_map(_ns(modality="guide", fastqs=[str(fq["r1"]), str(fq["r2"])], technology="0,0,16:0,16,28:1,27,47", parsed_seqspec=None, is_10x3v3=False,
                         spacer="", batch="FQ", outdir=None, label="st_map", engine="python", index=None, t2g=None, barcodes=str(fq["onlist"]), threads=1,
                         features=str(fq["meta"]), rev_comp=False, max_reads=0))
        o = latest("st_map")
        check(rc == 0 and (o / "FQ_ks_guide_out" / "counts_unfiltered" / "adata.h5ad").is_file() and (o / "FQ_ks_guide_out" / "run_info.json").is_file(),
              "map: kb-like <batch>_ks_guide_out (counts_unfiltered/adata.h5ad, inspect.json, run_info.json)")
        rc = cmd_feature_ref(_ns(kind="guide", table=str(fq["meta"]), rev_comp=True, spacer="", no_kb=True, label="st_fref"))
        o = latest("st_fref")
        gf = pd.read_csv(o / "guide_features.txt", sep="\t", header=None)
        check(rc == 0 and gf.iloc[0, 0] == revcomp(fq["features"]["g0"]) and (o / "guide_mismatch.fa").is_file(), "feature-ref: reverse-complemented guide_features.txt + kite FASTA")

        print("\nknee filter, MuData, guide assignment")
        W = synthetic_screen(d)
        rna = concat_ks_dirs(W["rna"])
        check(rna.obs_names[0].endswith("_0") and set(rna.obs["batch"]) == {"B1", "B2"} and rna.var_names.name == "gene_id",
              "concat: '<batch>_ks_' -> batch, obs names '<barcode>_<i>' (index_unique '_'), var index name kept")
        names = open(Path(W["rna"][0]) / "counts_unfiltered" / "cells_x_genes.genes.names.txt").read().split()
        filt, _, pinfo = preprocess_rna(rna, names, 500, 15.0, "human", "knee2", None, plots=False)
        tr = W["truth"]
        kept = set(filt.obs_names)
        real = [c for c, v in tr.items() if v[0]]
        empt = [c for c, v in tr.items() if not v[0]]
        hm = [c for c, v in tr.items() if v[0] and v[2]]
        check(sum(c in kept for c in empt) == 0 and sum(c in kept for c in real) > 0.85 * len(real),
              f"knee2 barcode filter: {sum(c in kept for c in empt)}/{len(empt)} empty droplets kept, {sum(c in kept for c in real)}/{len(real)} cells kept")
        check(sum(c in kept for c in hm) == 0 and filt.var_names[0] == W["gids"][0] and filt.var["symbol"].iloc[0] == "GENE0",
              "preprocess: pct_counts_mt < 15 removes the high-mito cells; Ensembl versions stripped, symbols from cells_x_genes.genes.names.txt")
        c1, c2 = elbow_knee_finder(np.arange(1.0, 6), np.array([5, 4.9, 4.8, 1, 0.9]), "basic"), barcode_filter_threshold(np.r_[np.full(50, 1000), np.full(50, 5)], "knee")
        check(c1 is not None and c2[0] == 1000, "elbow_knee_finder: max distance to the chord; knee threshold = count at the knee rank")
        guide = concat_ks_dirs(W["guide"])
        m = create_mudata(filt, guide, read_guide_metadata(W["meta"]), read_gtf(W["gtf"]), "high", "CROP-seq", None, None, plots=False)
        gv = m["guide"].var
        nt_names = sorted(gv.loc[gv["type"] == "non-targeting", "intended_target_name"].unique())
        check(nt_names == [f"non-targeting|{k}" for k in range(1, 5)] and (gv.loc[gv["type"] == "non-targeting", "intended_target_chr"] == "non-targeting").all(),
              "create-mudata: 8 NT guides bucketed into non-targeting|1..4 (median 2 guides per targeting element)")
        k0 = gv.loc[gv["guide_id"] == f"{W['syms'][10]}_sg1", "intended_target_key"].iloc[0]
        check(k0 == f"coords::{W['gids'][10]}::chr1::{W['starts'][10]}::{W['starts'][10] + 20000}" and {"guide_id", "targeting", "type", "intended_target_name", "spacer"} <= set(gv.columns),
              "create-mudata: intended_target_key coords::name::chr::start::end; guide.var carries the design columns")
        check(m["gene"].var["gene_start"].notna().all() and set(m["gene"].obs.columns) >= {"total_gene_umis", "percent_mito", "num_expressed_genes"}
              and "total_guide_umis" in m["guide"].obs.columns and list(m["guide"].uns["moi"]) == ["high"] and m.n_obs == len(set(filt.obs_names) & set(guide.obs_names)),
              "create-mudata: GTF gene coordinates, renamed obs QC columns, guide UMI totals, moi, barcode intersection")
        mp = write_h5mu(m, d / "mudata.h5mu", quiet=True)
        mr = read_h5mu(mp)
        check(set(mr.mod) == {"gene", "guide"} and mr.n_obs == m.n_obs and mr["guide"].var["intended_target_name"].astype(str).tolist() == gv["intended_target_name"].astype(str).tolist(),
              "h5mu round trip without the mudata package (mod/gene, mod/guide, obs, obsmap/varmap, uns)")
        truthA = np.zeros(m["guide"].shape, bool)
        gl = list(map(str, m["guide"].var_names))
        for i, c in enumerate(m["gene"].obs_names):
            for gname in tr[c][1]:
                truthA[i, gl.index(gname)] = True

        def f1(mm_):
            P = as_csr(mm_["guide"].layers["guide_assignment"]).toarray() > 0
            tp = (P & truthA).sum()
            return 2 * tp / (P.sum() + truthA.sum())
        m_s, ai = assign_guides(m, "sceptre", 0.8, 5)
        m_c, _ = assign_guides(m, "cleanser", 1.0)
        m_t, _ = assign_guides(m, "threshold", umi_threshold=5)
        m_x, axi = assign_guides(m, "cleanser", 1.0, engine="external", workdir=d / "ax")
        check((which("cleanser") is not None) or (axi["engine"] == ["python-cleanser"] and f1(m_x) == f1(m_c)),
              f"assign-guides --engine external without cleanser on PATH prints the command and falls back ({axi['engine']})")
        check(f1(m_s) > 0.95 and ai["n_batches"] == 2, f"assign-guides sceptre mixture (per batch, 5 EM restarts, posterior >= 0.8): F1 {f1(m_s):.3f} vs planted")
        check(f1(m_c) > 0.9 and f1(m_t) > 0.9, f"assign-guides cleanser fallback F1 {f1(m_c):.3f}; thresholding (UMI >= 5) F1 {f1(m_t):.3f}")
        mf = filter_genes_by_cells(m_s, 0.05)
        check(mf["gene"].n_vars < m_s["gene"].n_vars and all(W["gids"][i] in set(mf["gene"].var_names) for i in W["target_idx"]),
              f"mudata_concat gene filter: > int(0.05 n) cells ({m_s['gene'].n_vars} -> {mf['gene'].n_vars} genes, targets kept)")

        print("\npairs, chunking, inference")
        pairs = create_pairs_to_test(mf["guide"].var, mf["gene"].var_names, read_gtf(W["gtf"]), 1_000_000)
        own = all(((pairs["guide_id"] == f"{W['syms'][gi]}_sg1") & (pairs["gene_name"] == W["gids"][gi])).any() for gi in W["target_idx"])
        check(own and not pairs["guide_id"].str.startswith("NTC").any() and (pairs.groupby("guide_id").size() <= 11).all(),
              f"create_pairs_to_test: {len(pairs)} guide-gene pairs within 1 Mb; every guide paired with its target; NT guides have no pairs")
        mcis, pinf = prepare_inference(mf, pairs, subset_for_cis=True)
        pt_ = uns_df(mcis, "pairs_to_test")
        check("gene_id" in pt_.columns and "intended_target_key" in pt_.columns and mcis["gene"].n_vars == pinf["n_genes_tested"]
              and mcis["guide"].n_vars == 20 and mcis.n_obs <= mf.n_obs, "prepare_inference: gene_name -> gene_id, target keys, cis subset (tested genes, tested + control guides, cells with guides)")
        check(plan_chunk_ranges(10, 4) == [(0, 4), (4, 7), (7, 10)] and sceptre_chunk_plan(100, 10, 4, "auto", 10 ** 9) == (False, [(0, 10)])
              and sceptre_chunk_plan(100, 10, 4, "force")[1] == [(0, 4), (4, 8), (8, 10)], "chunk plans: balanced PerTurbo ranges; SCEPTRE auto/force")
        z = np.random.default_rng(1).standard_normal(4000)
        check(skewnorm_moment_pvalue(5.0, z) < 1e-4 and 0.3 < skewnorm_moment_pvalue(0.0, z) <= 1.0 and skewnorm_moment_pvalue(-3.0, z, "right") > 0.99,
              "skew-normal resampling approximation: tail, centre and one-sided p-values")
        im, tables, iinfo = run_inference(mcis, "default", mf, "python", 300, "both", d / "work")
        ce, te, tg = tables["cis_per_element_results"], tables["trans_per_element_results"], tables["trans_per_guide_results"]
        own_rows = pd.concat([ce[(ce["gene_id"] == W["gids"][gi]) & (ce["intended_target_name"] == W["gids"][gi])] for gi in W["target_idx"]])
        check(len(own_rows) == 6 and (own_rows["sceptre_p_value"] < 1e-3).all() and (own_rows["perturbo_log2_fc"] < -1.5).all() and (own_rows["sceptre_log2_fc"] < -1.5).all(),
              f"inference default/cis: all 6 targets knocked down (max SCEPTRE p {own_rows['sceptre_p_value'].max():.1e}, max log2FC {own_rows['perturbo_log2_fc'].max():.2f})")
        other = ce[ce["gene_id"] != ce["intended_target_name"]]
        check((other["sceptre_p_value"] < 0.01).mean() < 0.1, f"inference cis: neighbouring genes null (fraction SCEPTRE p < 0.01 = {(other['sceptre_p_value'] < 0.01).mean():.3f})")
        ntt = te[te["intended_target_name"].astype(str).str.startswith("non-targeting|")]
        check(len(ntt) and (ntt["p_value"] < 0.05).mean() < 0.12, f"inference trans: non-targeting elements calibrated (fraction p < 0.05 = {(ntt['p_value'] < 0.05).mean():.3f})")
        t1 = te[te["intended_target_name"] == W["gids"][W["target_idx"][0]]].copy()
        t1["q"] = bh_adjust(t1["p_value"].values)
        tset = {W["gids"][i] for i in W["trans_idx"]}
        hit = t1[t1["gene_id"].isin(tset)]
        miss = t1[~t1["gene_id"].isin(tset) & (t1["gene_id"] != W["gids"][W["target_idx"][0]])]
        down = (miss["q"] < 0.05) & (miss["log2_fc"] < 0)
        check((hit["q"] < 0.05).mean() >= 0.8 and (hit["log2_fc"] < -0.8).all() and down.sum() == 0 and (miss["log2_fc"].abs() < 0.4).all(),
              f"inference trans: T1 represses {int((hit['q'] < 0.05).sum())}/{len(hit)} planted genes (BH < 0.05); no other gene down; the others shift "
              f"up by the library-size composition only (median log2FC {miss['log2_fc'].median():+.2f})")
        check(list(ce.columns) == ["gene_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end", "sceptre_log2_fc",
                                   "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value"] and list(tg.columns) == ["gene_id", "guide_id", "log2_fc", "log2_fc_std", "p_value"]
              and set(im.uns) >= {"cis_per_guide_results", "cis_per_element_results", "trans_per_guide_results", "trans_per_element_results", "pairs_to_test"}
              and "per_guide_results" not in im.uns, "inference outputs: upstream column orders and uns keys (merge_method_results + merge_cis_trans_results)")
        check(len(tg) == 20 * mf["gene"].n_vars and len(te) == 10 * mf["gene"].n_vars, "trans = PerTurbo on all pairs (20 guides / 10 elements x all genes)")
        pgo, peo = export_outputs(im, tables["cis_per_guide_results"])
        check(list(pgo.columns) == ["intended_target_name", "guide_id(s)", "intended_target_chr", "intended_target_start", "intended_target_end", "gene_id", "sceptre_log2_fc",
                                    "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value", "cell_number", "avg_gene_expression"] and peo["guide_id(s)"].str.contains(",").any(),
              "export_output_single: per_guide_output / per_element_output schema; element rows join their guides")
        mg, me = merge_method_results(pd.DataFrame({"gene_id": ["a", "b"], "guide_id": ["g", "g"], "log2_fc": [1, 2], "p_value": [0.1, 0.2]}),
                                      pd.DataFrame({"gene_id": ["a"], "intended_target_name": ["t"], "intended_target_chr": ["chr1"], "intended_target_start": [1],
                                                    "intended_target_end": [2], "log2_fc": [1.0], "p_value": [0.1]}),
                                      pd.DataFrame({"gene_id": ["b", "c"], "guide_id": ["g", "g"], "log2_fc": [3, 4], "p_value": [0.3, 0.4]}),
                                      pd.DataFrame({"gene_id": ["a"], "intended_target_name": ["t"], "intended_target_chr": ["chr1"], "intended_target_start": [1.0],
                                                    "intended_target_end": [2.0], "log2_fc": [2.0], "p_value": [0.2]}))
        check(len(mg) == 3 and mg["sceptre_p_value"].isna().sum() == 1 and len(me) == 1 and me["perturbo_p_value"].iloc[0] == 0.2,
              "merge_method_results: outer guide merge; coordinate-aware element keys (1 == 1.0)")
        gsub = mf.subset_cells(list(mf["gene"].obs_names[:200]))
        one, _, _ = run_inference(gsub, "perturbo", None, "python", 50)
        check(set(one.uns) >= {"per_guide_results", "per_element_results"} and len(uns_df(one, "per_guide_results")) == 20 * gsub["gene"].n_vars,
              "inference perturbo without pairs: all-by-all per_guide_results / per_element_results")

        print("\nQC, evaluation, TF benchmark, dashboard")
        qd = d / "qc"
        qs = run_qc(im, qd, plots=not args.no_plots)
        it = qs["intended_target"]
        check(it["n_guides_tested"] == 12 and it["n_strong_knockdowns"] >= 10 and it["auroc"] > 0.95 and it["n_eval_negatives"] == it["n_eval_positives"],
              f"qc intended target on {qs['results_key']}: {it['n_strong_knockdowns']}/12 strong knockdowns, AUROC {it['auroc']:.3f} (balanced NT negatives)")
        trm = qs["trans"]
        check(trm["n_targeting_guides"] == 12 and trm["n_non_targeting_guides"] == 8 and trm["median_significant_per_guide_nt"] <= 1 and trm["total_significant_tests"] >= 10,
              f"qc trans: BH over all tests; NT guides median {trm['median_significant_per_guide_nt']} significant, total {trm['total_significant_tests']}")
        gmet = pd.read_csv(qd / "gene" / "gene_metrics.tsv", sep="\t")
        check(list(gmet["batch"]) == ["all", "B1", "B2"] and {"umi_median", "genes_q75", "mito_mean"} <= set(gmet.columns) and (qd / "qc_report.html").is_file()
              and (qd / "guide" / "guide_metrics.tsv").is_file(), "qc gene/guide metric tables (overall + per batch) and qc_report.html")
        ev = evaluate_controls(im, d / "ev", plots=False)
        check(ev["auroc"] > 0.95 and ev["n_positives"] == 12 and ev["n_negatives"] == 12 and not ev["warning"],
              f"evaluate_controls: AUROC {ev['auroc']:.3f}, AUPRC {ev['auprc']:.3f} (12 direct-target tests vs 12 NT)")
        bp, bg = igv_tracks(im, read_gtf(W["gtf"]), "cis_per_element_results", "sceptre")
        check(len(bg) == 6 and len(bp) == len(ce) - 6 and list(bp.columns[:3]) == ["chr1", "start1", "end1"], "igv: 6 promoter bedgraph rows, element->gene bedpe for the rest")
        nodes = select_central_nodes(te, 1, "intended_target_name", "gene_id", "log2_fc", 0.1)
        check(nodes == [W["gids"][W["target_idx"][0]]], f"network: central node by weighted degree = T1 ({nodes})")
        ex = perturbation_expression(im, f"{W['syms'][10]}_sg1", W["gids"][10])
        gm_ = ex.groupby("group")["expression"].mean()
        check(set(gm_.index) == {f"{W['syms'][10]}_sg1", NT_LABEL} and gm_[f"{W['syms'][10]}_sg1"] < 0.5 * gm_[NT_LABEL],
              "expression (utlis.get_perturbation_expression): guide cells vs NT-only cells")
        pk = d / "peaks"
        pk.mkdir()
        gtfd = read_gtf(W["gtf"])
        tfid = W["gids"][W["target_idx"][0]]
        rows_ = []
        for i in W["trans_idx"] + [3, 50]:
            tss = W["starts"][i] if i % 2 == 0 else W["starts"][i] + 20000
            rows_.append(("chr1", tss - 200, tss + 200, ".", 800))
        pd.DataFrame(rows_).to_csv(pk / "t1.bed.gz", sep="\t", header=False, index=False)
        tfb = run_tf_benchmark(im, gtfd, {tfid: str(pk / "t1.bed.gz")}, d / "bench", [500, 5000], 0.001, plots=False, name_map={tfid: "T1"})
        b0 = tfb["best_per_tf"][0]
        check(b0["odds_ratio"] > 5 and b0["pvalue"] < 1e-4 and b0["num_targets"] == 17, f"tf-benchmark: planted ChIP targets enriched among T1 trans hits (OR {b0['odds_ratio']:.1f}, p {b0['pvalue']:.1e})")
        allb = pd.read_csv(d / "bench" / "benchmark_tables" / "enrichment_all.tsv", sep="\t")
        pw = promoter_windows(gtfd, 500)
        pw_c = promoter_windows(gtfd, 500, compat=True)
        check(set(allb["promoter_window_width"]) == {500, 5000} and (pw_c["TSS"] != pw["TSS"]).sum() == (np.array(W["starts"]) >= 0).sum() // 2 + 1,
              "tf-benchmark: strand-aware TSS; --upstream-compat reproduces END-as-TSS for every gene")
        kbj = collect_kb_json([Path(x) for x in W["rna"] + W["guide"]])
        dash = build_dashboard(im, d, [], kbj, flow_html([("B1", 100)], 10, 8, 7, 6, 6, {"barcode_filter": "knee2"}), [qd], True)
        txt = dash.read_text()
        check(len(kbj) == 4 and "Filtering Summary" in txt and "Mapping scRNA B1" in txt and "Cis Analysis" in txt, "dashboard: kb JSON mapping blocks, filtering flow, cis/trans inference tables")

        print("\nhashing, doublets, dual guides, control targets")
        rng = np.random.default_rng(9)
        n = 400
        lab = rng.choice(3, n)
        state = rng.random(n)
        H = rng.poisson(3, (n, 3)).astype(float)
        for i in range(n):
            if state[i] < 0.8:
                H[i, lab[i]] += rng.poisson(150)
            elif state[i] < 0.9:
                H[i, lab[i]] += rng.poisson(150)
                H[i, (lab[i] + 1) % 3] += rng.poisson(150)
        ha = ad.AnnData(sp.csr_matrix(H), obs=pd.DataFrame({"batch": ["B1"] * n}, index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=["HTO1", "HTO2", "HTO3"]))
        hf, hu, hinfo = run_hashing(ha, ad.AnnData(np.zeros((n, 1)), obs=pd.DataFrame(index=[f"c{i}" for i in range(n)])), d / "hash")
        exp_type = np.array([f"HTO{l + 1}" if s < 0.8 else ("multiplets" if s < 0.9 else "negative") for l, s in zip(lab, state)])
        acc = (hu.obs["hto_type_split"].astype(str).values == exp_type).mean()
        check(acc > 0.95 and hf.n_obs == int(((hu.obs["hto_type"] != "negative") & (hu.obs["hto_type_split"] != "multiplets")).sum()),
              f"hashing: GMM demultiplex fallback {acc:.3f} accurate; singlets kept, negatives and multiplets dropped")
        rng = np.random.default_rng(4)
        prog = rng.dirichlet(np.ones(60) * 0.3, 4)
        ct = rng.choice(4, 500)
        Xs = rng.poisson(prog[ct] * 800).astype(float)
        dbl = rng.choice(500, 40, replace=False)
        partner = rng.choice(500, 40)
        Xs[dbl] = Xs[dbl] + Xs[partner]
        sc_, pred, thr, _ = scrublet_fallback(sp.csr_matrix(Xs))
        auc_d, _, _ = roc_pr(np.isin(np.arange(500), dbl).astype(int), sc_)
        check(auc_d > 0.8, f"doublets: scrublet fallback scores planted doublets (AUROC {auc_d:.3f}, threshold {thr:.3f})")
        gvd = pd.DataFrame({"guide_id": ["a1", "a2", "b1", "b2"], "intended_target_name": ["A", "A", "B", "B"], "spacer": ["AA", "CC", "GG", "TT"],
                            "intended_target_chr": ["chr1"] * 4, "intended_target_start": [1.0] * 4, "intended_target_end": [2.0] * 4, "targeting": [True] * 4,
                            "type": ["targeting"] * 4}, index=["a1", "a2", "b1", "b2"])
        Ad = np.array([[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0], [1, 0, 0, 0]], float)
        gd_ = ad.AnnData(sp.csr_matrix(Ad * 10), obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]), var=gvd)
        gd_.layers["guide_assignment"] = sp.csr_matrix(Ad)
        md_ = MData({"gene": ad.AnnData(np.ones((4, 2)), obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"])), "guide": gd_})
        col_, cinfo = collapse_guides(md_)
        check(list(col_["guide"].var_names) == ["a1|a2", "b1|b2"] and col_.n_obs == 2 and col_["guide"].var.loc["a1|a2", "spacer"] == "AA|CC",
              "collapse_guides (DUAL_GUIDE): 2 guides on one element -> 'g1|g2'; cells with 2 elements or 1 guide dropped")
        ntv = nt_intended_targets(pd.DataFrame({"guide_id": ["t1", "t2", "n1", "n2", "n3"], "targeting": [True, True, False, False, False],
                                                "type": ["targeting", "targeting", "non-targeting", "non-targeting", "non-targeting"],
                                                "intended_target_name": ["X", "X", "nt", "nt", "nt"]}), "median")
        check(list(ntv["intended_target_name"].astype(str)) == ["X", "X", "non-targeting|1", "non-targeting|1", "non-targeting|2"],
              "create_nt_intended_targets median strategy")

        print("\nCLI wrappers (main([...]) for each subcommand on the synthetic screen)")
        clis = []

        def cli(argv, label, expect):
            try:
                rc_ = main(argv + ["--label", label])
            except SystemExit as e:
                rc_ = e.code
            o_ = latest(label)
            ok_ = rc_ == 0 and (o_ / "report.md").is_file() and all((o_ / x).exists() for x in expect)
            clis.append((argv[0], ok_))
            return o_
        o1 = cli(["concat", "--inputs", *W["rna"]], "st_c_rna", ["concatenated_adata.h5ad"])
        o2 = cli(["concat", "--inputs", *W["guide"]], "st_c_guide", ["concatenated_adata.h5ad"])
        o3 = cli(["preprocess", "--adata", str(o1 / "concatenated_adata.h5ad"), "--gene-names", W["rna"][0]] + (["--no-plots"] if args.no_plots else []),
                 "st_c_pre", ["filtered_anndata.h5ad"])
        o4 = cli(["create-mudata", "--rna", str(o3 / "filtered_anndata.h5ad"), "--guide", str(o2 / "concatenated_adata.h5ad"), "--guide-metadata", str(W["meta"]),
                  "--gtf", str(W["gtf"]), "--no-plots"], "st_c_cm", ["mudata.h5mu"])
        cli(["doublets", "--mudata", str(o4 / "mudata.h5mu"), "--no-plots"], "st_c_dbl", ["mdata_doublets.h5mu"])
        o5 = cli(["assign-guides", "--mudata", str(o4 / "mudata.h5mu"), "--n-em-rep", "2"], "st_c_ag", ["concat_mudata.h5mu"])
        cli(["nt-targets", "--mudata", str(o5 / "concat_mudata.h5mu"), "--strategy", "median"], "st_c_nt", ["mudata_nt_targets.h5mu"])
        o6 = cli(["pairs", "--mudata", str(o5 / "concat_mudata.h5mu"), "--gtf", str(W["gtf"])], "st_c_pairs", ["pairs_to_test.csv", "mudata_inference_input.h5mu"])
        cli(["chunk", "--mudata", str(o6 / "mudata_inference_input.h5mu"), "--chunk-size", "20", "--chunk-mode", "force"], "st_c_chunk", ["chunk_manifest.tsv", "chunks/chunk.000.h5mu"])
        o7 = cli(["inference", "--mudata", str(o6 / "mudata_inference_input.h5mu"), "--trans-mudata", str(o5 / "concat_mudata.h5mu"), "--n-resamples", "100"],
                 "st_c_inf", ["inference_mudata.h5mu", "per_guide_output.tsv", "per_element_output.tsv", "trans_per_element_results.tsv.gz"])
        imp = str(o7 / "inference_mudata.h5mu")
        cli(["merge-results", "--mode", "cis-trans", "--base-mudata", str(o5 / "concat_mudata.h5mu")] + sum(
            [[f"--{k.replace('_', '-')}", str(o7 / f"{k}_results.tsv.gz")] for k in ("cis_per_guide", "cis_per_element", "trans_per_guide", "trans_per_element")], []),
            "st_c_mr", ["inference_mudata.h5mu"])
        cli(["merge-results", "--mode", "export", "--base-mudata", imp], "st_c_mre", ["per_guide_output.tsv"])
        cli(["qc", "--mudata", imp, "--no-plots"], "st_c_qc", ["qc_report.html", "intended_target/intended_target_metrics.tsv", "trans/trans_metrics.tsv"])
        cli(["evaluate", "--mudata", imp, "--gtf", str(W["gtf"]), "--no-plots"], "st_c_ev", ["evaluation_output/cis_sceptre.bedpe", "evaluation_output/controls_evaluation_table.tsv"])
        cli(["expression", "--mudata", imp, "--guide-id", f"{W['syms'][10]}_sg1", "--gene-id", W["gids"][10], "--no-plots"], "st_c_expr", ["perturbation_expression.tsv"])
        cli(["tf-benchmark", "--mudata", imp, "--gtf", str(W["gtf"]), "--peaks", f"{tfid}={pk / 't1.bed.gz'}", "--no-plots"], "st_c_tf", ["benchmark_tables/enrichment_all.tsv"])
        cli(["dashboard", "--mudata", imp, "--ks-dirs", *W["rna"], *W["guide"], "--gene-ann", str(o1 / "concatenated_adata.h5ad"), "--gene-ann-filtered",
             str(o3 / "filtered_anndata.h5ad"), "--guide-ann", str(o2 / "concatenated_adata.h5ad"), "--no-plots"], "st_c_dash", ["dashboard.html"])
        cli(["seqspec-check", "--read1", str(fq["r1"]), "--read2", str(fq["r2"]), "--metadata", str(fq["meta"]), "--no-plots"], "st_c_ssc", ["position_table.csv", "best_config.csv"])
        sp_ = d / "ss.csv"
        ss.to_csv(sp_, index=False)
        cli(["samplesheet", "--samplesheet", str(sp_)], "st_c_ss", ["parse_covariate.csv", "cov_string.txt", "fastq_batches.txt"])
        pt.to_csv(d / "ps.tsv", sep="\t", index=False)
        cli(["portal-download", "--sample", str(d / "ps.tsv")], "st_c_pdl", ["per_sample_local.tsv"])
        pd.DataFrame({"tab_name": ["scRNA"], "batch_name": ["b1"], "read1": ["r1"], "read2": ["r2"], "variable": ["sequence"], "variable_value": [""]}).to_csv(d / "cfg.csv", index=False)
        cli(["interface-config", "--config-table", str(d / "cfg.csv")], "st_c_ic", ["pipeline_input.config"])
        ha.write_h5ad(d / "hash.h5ad")
        ad.AnnData(np.zeros((n, 1)), obs=pd.DataFrame(index=[f"c{i}" for i in range(n)])).write_h5ad(d / "hrna.h5ad")
        cli(["hashing", "--hashing", str(d / "hash.h5ad"), "--rna", str(d / "hrna.h5ad")], "st_c_hash", ["concatenated_hashing_demux.h5ad", "concatenated_unfiltered_hashing_demux.h5ad"])
        bad = [c_ for c_, ok_ in clis if not ok_]
        check(not bad, f"CLI: {len(clis)} subcommand invocations exit 0 with report.md + expected outputs" + (f" (failed: {bad})" if bad else ""))

        print("\nend-to-end run subcommand")
        rc = cmd_run(_ns(label="st_run", rna=W["rna"], guide=W["guide"], hashing=None, guide_metadata=str(W["meta"]), gtf=str(W["gtf"]), min_genes=500,
                         pct_mito=15.0, reference="human", barcode_filter="knee2", moi="high", capture_method="CROP-seq", scrublet=False,
                         assignment_method="sceptre", assignment_threshold=0.8, n_em_rep=3, umi_threshold=5, min_cells_fraction=0.05, dual_guide=False,
                         pairing_strategy="default", pairs=None, max_target_distance=1_000_000, inference_method="default", engine="python", n_resamples=200,
                         side="both", upstream_bin=None, seed=0, num_nodes=1, encode_bed_dir=None, peaks=None, upstream_compat=False, no_plots=args.no_plots))
        o = latest("st_run")
        summ = json.loads((o / "summary.json").read_text())
        check(rc == 0 and (o / "inference_mudata.h5mu").is_file() and (o / "dashboard.html").is_file() and (o / "report.md").is_file()
              and summ["qc"]["intended_target"]["auroc"] > 0.95 and summ["evaluate"]["controls"]["auroc"] > 0.95,
              f"run: count matrices -> dashboard in {summ['seconds']} s; AUROC {summ['qc']['intended_target']['auroc']:.3f}")
        rm = read_h5mu(o / "inference_mudata.h5mu")
        check(list(rm.uns["inference_engine"]) == ["python-perturbo", "python-sceptre"] and "guide_assignment" in rm["guide"].layers
              and (o / "evaluation_output" / "cis_sceptre.bedgraph").is_file(), "run: inference_mudata.h5mu records the engines; evaluation_output tracks written")
        rc = cmd_nextflow(_ns(label="st_nf", pipeline=UPSTREAM_REPO, revision=UPSTREAM_COMMIT, profile="local", input="samplesheet.csv", outdir="out",
                              param=["INFERENCE_method=default"], resume=False, execute=False))
        o = latest("st_nf")
        check(rc == 0 and "--INFERENCE_method default" in (o / "nextflow_command.sh").read_text(), "nextflow: exact upstream command written (not executed)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent crispr-pipeline", description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch ENCODE BED files + example seqspecs from the pinned upstream commit")
    p.add_argument("--dest")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("samplesheet", help="validate a samplesheet; covariates, formula, FASTQ groups")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--from-portal", action="store_true", help="input is a portal per-sample TSV (update_samplesheet.py mapping)")
    p.add_argument("--label", default="samplesheet")
    p.set_defaults(func=cmd_samplesheet)

    p = sub.add_parser("portal-samplesheet", help="IGVF analysis set -> per-sample TSV + samplesheet.csv")
    p.add_argument("--accession", required=True)
    p.add_argument("--rna-seqspec")
    p.add_argument("--sgrna-seqspec")
    p.add_argument("--hash-seqspec")
    p.add_argument("--label")
    p.set_defaults(func=cmd_portal_samplesheet)

    p = sub.add_parser("portal-download", help="download + md5-verify portal files of a per-sample TSV (dry run by default)")
    p.add_argument("--sample", required=True)
    p.add_argument("--file-types", choices=["fastq", "other", "all"], default="all")
    p.add_argument("--output-dir")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--label", default="portal_download")
    p.set_defaults(func=cmd_portal_download)

    p = sub.add_parser("interface-config", help="config table -> pipeline_input.config")
    p.add_argument("--config-table", required=True)
    p.add_argument("--mapping-tabs", nargs="+", default=["scRNA", "Guides", "Hash"])
    p.add_argument("--label", default="interface_config")
    p.set_defaults(func=cmd_interface_config)

    p = sub.add_parser("seqspec", help="seqspec YAML -> kb technology string")
    p.add_argument("--yaml", required=True)
    p.add_argument("--modalities", nargs="+", default=["rna"])
    p.add_argument("--whitelist")
    p.add_argument("--label", default="seqspec")
    p.set_defaults(func=cmd_seqspec)

    p = sub.add_parser("seqspec-check", help="scan R1/R2 for guide / hashtag positions and orientation")
    p.add_argument("--read1", nargs="+", required=True)
    p.add_argument("--read2", nargs="+", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--max-reads", type=int, default=100000)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="seqspec_check")
    p.set_defaults(func=cmd_seqspec_check)

    p = sub.add_parser("feature-ref", help="guide / hashing feature table + kite mismatch FASTA + kb ref")
    p.add_argument("--kind", choices=["guide", "hash"], default="guide")
    p.add_argument("--table", required=True, help="guide_metadata.tsv (guide_id, spacer) or hash metadata (hash_id, sequence)")
    p.add_argument("--rev-comp", action="store_true")
    p.add_argument("--spacer", default="")
    p.add_argument("--no-kb", action="store_true")
    p.add_argument("--label", default="feature_ref")
    p.set_defaults(func=cmd_feature_ref)

    p = sub.add_parser("map", help="kb count wrapper or Python feature counting (guide / hash)")
    p.add_argument("--modality", choices=["rna", "guide", "hash"], required=True)
    p.add_argument("--fastqs", nargs="+", required=True)
    p.add_argument("--technology", help="kb -x string, e.g. 0,0,16:0,16,28:1,0,0")
    p.add_argument("--parsed-seqspec")
    p.add_argument("--is-10x3v3", action="store_true")
    p.add_argument("--spacer", default="")
    p.add_argument("--batch", required=True)
    p.add_argument("--features", help="guide metadata / hash table (python engine)")
    p.add_argument("--rev-comp", action="store_true")
    p.add_argument("--barcodes")
    p.add_argument("--index")
    p.add_argument("--t2g")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--engine", choices=["auto", "kb", "python"], default="auto")
    p.add_argument("--max-reads", type=int, default=0)
    p.add_argument("--outdir")
    p.add_argument("--label")
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("concat", help="<batch>_ks_*_out dirs -> one AnnData")
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--covariates", help="parse_covariate.csv")
    p.add_argument("--output")
    p.add_argument("--label", default="concat")
    p.set_defaults(func=cmd_concat)

    p = sub.add_parser("preprocess", help="knee barcode filter + QC + filters")
    p.add_argument("--adata", required=True)
    p.add_argument("--gene-names", help="cells_x_genes.genes.names.txt or a ks_transcripts_out dir")
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--pct-mito", type=float, default=15.0)
    p.add_argument("--reference", default="human")
    p.add_argument("--barcode-filter", choices=["none", "knee", "knee2"], default="knee2")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="preprocess")
    p.set_defaults(func=cmd_preprocess)

    p = sub.add_parser("create-mudata", help="RNA + guide (+ hashing) -> mudata.h5mu")
    p.add_argument("--rna", required=True)
    p.add_argument("--guide", required=True)
    p.add_argument("--guide-metadata", required=True)
    p.add_argument("--gtf")
    p.add_argument("--hashing")
    p.add_argument("--moi", default="high")
    p.add_argument("--capture-method", default="CROP-seq")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="create_mudata")
    p.set_defaults(func=cmd_create_mudata)

    p = sub.add_parser("hashing", help="hashing filter + demultiplex")
    p.add_argument("--hashing", required=True)
    p.add_argument("--rna", required=True, help="filtered RNA h5ad")
    p.add_argument("--label", default="hashing")
    p.set_defaults(func=cmd_hashing)

    p = sub.add_parser("doublets", help="doublet removal on a MuData")
    p.add_argument("--mudata", required=True)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="doublets")
    p.set_defaults(func=cmd_doublets)

    p = sub.add_parser("assign-guides", help="guide assignment + gene filter (+ dual-guide collapse)")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["sceptre", "cleanser", "threshold"], default="sceptre")
    p.add_argument("--threshold", type=float)
    p.add_argument("--n-em-rep", type=int, default=5)
    p.add_argument("--umi-threshold", type=int, default=5)
    p.add_argument("--batch-key", default="batch")
    p.add_argument("--capture-method", default="CROP-seq")
    p.add_argument("--min-cells-fraction", type=float, default=0.05)
    p.add_argument("--dual-guide", action="store_true")
    p.add_argument("--engine", choices=["python", "external"], default="python", help="external: cleanser binary / upstream assign_grnas_sceptre.R")
    p.add_argument("--upstream-bin", help="upstream bin/ (assign_grnas_sceptre.R) for --engine external")
    p.add_argument("--label", default="assign_guides")
    p.set_defaults(func=cmd_assign_guides)

    p = sub.add_parser("nt-targets", help="control-guide intended targets (same / median)")
    p.add_argument("--mudata", required=True)
    p.add_argument("--strategy", choices=["same", "median"], default="same")
    p.add_argument("--label", default="nt_targets")
    p.set_defaults(func=cmd_nt_targets)

    p = sub.add_parser("pairs", help="pairs_to_test + prepare_inference")
    p.add_argument("--mudata", required=True)
    p.add_argument("--strategy", choices=["default", "by_distance", "predefined_pairs", "all"], default="default")
    p.add_argument("--gtf")
    p.add_argument("--pairs")
    p.add_argument("--limit", type=int, default=1_000_000)
    p.add_argument("--label", default="pairs")
    p.set_defaults(func=cmd_pairs)

    p = sub.add_parser("chunk", help="gene-chunked MuData + manifest")
    p.add_argument("--mudata", required=True)
    p.add_argument("--mode", choices=["sceptre", "perturbo"], default="sceptre")
    p.add_argument("--chunk-size", type=int, default=1000)
    p.add_argument("--chunk-mode", choices=["auto", "off", "force"], default="auto")
    p.add_argument("--auto-threshold-entries", type=int, default=2147483647)
    p.add_argument("--keep-all-guides", action="store_true")
    p.add_argument("--label", default="chunk")
    p.set_defaults(func=cmd_chunk)

    p = sub.add_parser("inference", help="SCEPTRE / PerTurbo cis + trans testing")
    p.add_argument("--mudata", required=True, help="inference input (pairs_to_test in uns)")
    p.add_argument("--trans-mudata", help="full MuData for trans (default method); defaults to --mudata")
    p.add_argument("--method", choices=["sceptre", "perturbo", "sceptre,perturbo", "default"], default="default")
    p.add_argument("--engine", choices=["python", "external"], default="python")
    p.add_argument("--upstream-bin", help="upstream bin/ with inference_sceptre.R and perturbo_inference.py (engine external)")
    p.add_argument("--n-resamples", type=int, default=500)
    p.add_argument("--side", choices=["both", "left", "right"], default="both")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="inference")
    p.set_defaults(func=cmd_inference)

    p = sub.add_parser("merge-results", help="merge method / cis-trans / chunk results; export outputs")
    p.add_argument("--mode", choices=["methods", "cis-trans", "chunks", "legacy", "add-test-results", "export"], required=True)
    p.add_argument("--base-mudata")
    for a in ("sceptre-per-guide", "sceptre-per-element", "perturbo-per-guide", "perturbo-per-element", "cis-per-guide", "cis-per-element",
              "trans-per-guide", "trans-per-element", "chunk-manifest", "cis-mudata", "trans-mudata", "test-results-csv"):
        p.add_argument(f"--{a}")
    p.add_argument("--per-guide-files", nargs="+")
    p.add_argument("--per-element-files", nargs="+")
    p.add_argument("--results-key", default="auto")
    p.add_argument("--label", default="merge_results")
    p.set_defaults(func=cmd_merge_results)

    p = sub.add_parser("qc", help="gene / guide / intended-target / trans QC + HTML report")
    p.add_argument("--mudata", required=True)
    p.add_argument("--results-key", default="auto")
    p.add_argument("--log2fc-col", default="auto")
    p.add_argument("--pvalue-col", default="auto")
    p.add_argument("--fc-threshold", type=float, default=0.4)
    p.add_argument("--pval-threshold", type=float, default=0.05)
    p.add_argument("--fdr-method", choices=["fdr_bh", "bonferroni", "none"], default="fdr_bh")
    p.add_argument("--validated-links", help="TSV guide_target, gene (symbols)")
    p.add_argument("--batch-col", default="batch")
    p.add_argument("--match-symbol", default="auto", help="match intended_target_name to gene symbols: auto / true / false")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="qc")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("evaluate", help="evaluate_controls + IGV tracks + volcano / network plots")
    p.add_argument("--mudata", required=True)
    p.add_argument("--gtf")
    p.add_argument("--num-nodes", type=int, default=1)
    p.add_argument("--min-weight", type=float, default=0.1)
    p.add_argument("--central-nodes", nargs="+")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="evaluate")
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("expression", help="gene expression: guide cells vs NT-only cells")
    p.add_argument("--mudata", required=True)
    p.add_argument("--guide-id", required=True)
    p.add_argument("--gene-id", required=True)
    p.add_argument("--layer")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="expression")
    p.set_defaults(func=cmd_expression)

    p = sub.add_parser("tf-benchmark", help="ENCODE ChIP-seq TF target enrichment of trans results")
    p.add_argument("--mudata", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--encode-bed-dir")
    p.add_argument("--peaks", nargs="+", help="TF_ID=peaks.bed(.gz) (overrides the 5 bundled ENCODE files)")
    p.add_argument("--windows", type=int, nargs="+", default=PROMOTER_WINDOWS)
    p.add_argument("--cutoff", type=float, default=0.001)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="tf_benchmark")
    p.set_defaults(func=cmd_tf_benchmark)

    p = sub.add_parser("dashboard", help="dashboard plots + HTML")
    p.add_argument("--mudata", required=True)
    p.add_argument("--ks-dirs", nargs="+", help="<batch>_ks_*_out dirs (mapping JSON)")
    p.add_argument("--gene-ann")
    p.add_argument("--gene-ann-filtered")
    p.add_argument("--guide-ann")
    p.add_argument("--hashing-unfiltered-demux")
    p.add_argument("--image-dirs", nargs="+")
    p.add_argument("--barcode-filter", default="knee2")
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--pct-mito", type=float, default=15.0)
    p.add_argument("--not-default", action="store_true", help="single test_results layout instead of cis/trans")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="dashboard")
    p.set_defaults(func=cmd_dashboard)

    p = sub.add_parser("run", help="count matrices -> dashboard, all stages")
    p.add_argument("--rna", nargs="+", required=True, help="<batch>_ks_transcripts_out dirs or one .h5ad")
    p.add_argument("--guide", nargs="+", required=True, help="<batch>_ks_guide_out dirs or one .h5ad")
    p.add_argument("--hashing", nargs="+")
    p.add_argument("--guide-metadata", required=True)
    p.add_argument("--gtf")
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--pct-mito", type=float, default=15.0)
    p.add_argument("--reference", default="human")
    p.add_argument("--barcode-filter", choices=["none", "knee", "knee2"], default="knee2")
    p.add_argument("--moi", default="high")
    p.add_argument("--capture-method", default="CROP-seq")
    p.add_argument("--scrublet", action="store_true")
    p.add_argument("--assignment-method", choices=["sceptre", "cleanser", "threshold"], default="sceptre")
    p.add_argument("--assignment-threshold", type=float)
    p.add_argument("--n-em-rep", type=int, default=5)
    p.add_argument("--umi-threshold", type=int, default=5)
    p.add_argument("--min-cells-fraction", type=float, default=0.05)
    p.add_argument("--dual-guide", action="store_true")
    p.add_argument("--pairing-strategy", choices=["default", "by_distance", "predefined_pairs", "all_by_all"], default="default")
    p.add_argument("--pairs")
    p.add_argument("--max-target-distance", type=int, default=1_000_000)
    p.add_argument("--inference-method", choices=["sceptre", "perturbo", "sceptre,perturbo", "default"], default="default")
    p.add_argument("--engine", choices=["python", "external"], default="python")
    p.add_argument("--upstream-bin")
    p.add_argument("--n-resamples", type=int, default=500)
    p.add_argument("--side", choices=["both", "left", "right"], default="both")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--num-nodes", type=int, default=1)
    p.add_argument("--encode-bed-dir")
    p.add_argument("--peaks", nargs="+")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="crispr_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("nextflow", help="print (or --execute) the upstream Nextflow command")
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--profile", default="local")
    p.add_argument("--pipeline", default=UPSTREAM_REPO)
    p.add_argument("--revision", default=UPSTREAM_COMMIT)
    p.add_argument("--param", nargs="+", help="NAME=VALUE nextflow params")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--label", default="nextflow")
    p.set_defaults(func=cmd_nextflow)

    p = sub.add_parser("selftest", help="synthetic screen with planted effects; every subcommand")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    if getattr(args, "match_symbol", None) is not None and args.cmd == "qc":
        args.match_symbol = "auto" if args.match_symbol == "auto" else _truthy(args.match_symbol)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
