#!/usr/bin/env python3
"""Perturb-seq guide/expression pipeline, QC and SCEPTRE-style differential perturbation (port of IGVF-CRISPR/bulk_crispr_pipeline).

Port of https://github.com/IGVF-CRISPR/bulk_crispr_pipeline (branch master, commit
6f4ae09fa7f7dffc1d98b3dc4dc33e73c8ff48ac, 2023-10-25).  No LICENSE file: the
METHODS (definitions, thresholds, column names, output schemas) were read from
the Nextflow workflow, the 13 bin/ scripts, the three .config files and the
nine notebooks, and re-derived in Python; no code was copied.  Relationship:
port.

Despite its name, the repository is NOT a bulk (sorting / proliferation) screen
pipeline: its README calls it "Pipeline single Cell Perturb-seq like".  It is
Lucas Silva Ferreira's pipeline_perturbseq_like (the predecessor of
IGVF_CRISPR_Pipeline) demonstrated on the Gasperini 2019 pilot: kallisto|bustools
(kite) guide and cDNA quantification, per-lane cell QC (knee, mito, Scrublet
doublets), guide binarisation, optional MULTI-seq demultiplexing, cis gene
selection around each targeted element, SCEPTRE high-MOI tests per guide x gene,
BH / Fisher aggregation into a MuData.  There is no MAGeCK / RRA / LFC-of-bins
step in it, so none is invented here.

Definitions, exactly as upstream computes them
  read composition   first 10,000 FASTQ lines; per position the counts of
                     A/C/T/G; bias = sd (ddof 1) of the four counts.
  guide table        rows grouped by Target_name (sorted); Target_name ->
                     "<target>|<n>" (n = 1.. within group);
                     pipeline_id = "<target>|<n>_sgrna_<chr>:<start>:<end>";
                     guide_features.txt = sgRNA_sequences \\t pipeline_id.
  guide counting     kite: every guide plus all its Hamming-1 variants that do
                     not collide with another guide; barcode/UMI/sequence
                     positions from the kallisto technology string
                     (10XV2 0,0,16:0,16,26:1,0,0; 10XV3 0,0,16:0,16,28:1,0,0);
                     barcodes corrected to a whitelist at Hamming 1 when unique;
                     a UMI is counted once per (barcode, guide) and dropped when
                     it maps to more than one guide.
  lane manifest      file_path = <dir>/counts_unfiltered/adata.h5ad; sample =
                     "Guide" if "guide" in path else "scRNA"; lane =
                     <dirname>.split("_")[1] without "L".
  per-lane QC        knee = total UMIs sorted descending; cells with n_genes >=
                     TRANSCRIPTS_UMI_TRHESHOLD (default 100 -- upstream applies
                     the UMI threshold as a minimum GENE count), then n_counts >=
                     knee[EXPECTED_CELL_NUMBER] (default 8000); Ensembl version
                     stripped; percent_mito = MT- counts / total, keep < 0.2;
                     Scrublet doublets removed; guide number_of_nonzero_guides;
                     barcodes kept only if present in both modalities;
                     batch_number = lane.
  gene filter        after concatenating lanes, genes detected in >=
                     int(n_cells * 0.01) cells.
  guide annotation   guide_chr / guide_start / guide_end / guide_number /
                     target_elements parsed from pipeline_id; gene TSS = start
                     (+) or end (-) from the GTF gene records; missing genes
                     -> NOT_FOUND, 0, 0.
  binarisation       guides summed by name (or by target with MERGE), binary =
                     UMIs > GUIDE_UMI_LIMIT (strict, default 5).
  covariates         bath_number (sic), percent_mito,
                     log_number_of_detected_genes = log(n_genes),
                     log_total_gene_count, log_total_guide_count =
                     log(sum of binarised guides + 1).
  MULTI-seq          deMULTIplex: cell BC R1[1..16], UMI R1[17..28], tag
                     R2[1..8] (1-based inclusive); tag -> reference barcode at
                     Hamming <= 1; UMIs counted per cell x barcode; classifyCells
                     (log2, -Inf -> 0, mean-centre, bkde normal kernel, 100-point
                     grid between the 0.1% and 99.9% quantiles, threshold =
                     quantile(c(x[high max], x[low max]), q)); quantile sweep
                     q = 0.01..0.99 by 0.02, findThresh / findQ (singlet
                     proportion local maxima), two rounds removing negatives.
  cis selection      element = target before "|"; "_TSS" stripped; genes on the
                     same chromosome with |TSS - element TSS| < DISTANCE_NEIGHBORS
                     (1,000,000); IN_TRANS "TRUE" -> all genes; elements whose
                     name is not a GTF gene get 10 random genes; ADDGENENAMES
                     always added.  GUIDE_TYPE: "_TSS" POSITIVE_CONTROL,
                     random/scrambled NEGATIVE_CONTROL, "chr" PUTATIVE_ENHANCER.
  testing            per element, guides present in > 30 cells; batch
                     covariate dropped when single-valued; side = DIRECTION
                     (left / right / both); SCEPTRE conditional resampling.
                     Python engine: NB GLM of the gene on the covariates
                     (Poisson IRLS + ML theta), logistic propensity of the
                     guide, B (500) Bernoulli resamples of the guide
                     indicator, efficient-score z; p = resampling p, with a
                     fitted skew-normal tail (skew-t in the R package) below
                     10 / (B + 1); log_fold_change = ln(sum y / sum mu) in
                     guide-positive cells.  Engine sceptre-r writes and runs
                     the upstream run_sceptre_high_moi script instead.
  aggregation        adj_pvalue = BH over every test; significant = adj < 0.01;
                     element p = Fisher combination of its guides' p per gene;
                     sig_not_adj = element p < 0.05.

Upstream bugs corrected by default (reproduce with --upstream-compat)
  * fq_composition keeps every FASTQ line without "@", "+" or "F", so quality
    strings leak into the composition; the port uses sequence lines only.
  * concact_and_pre_filtering never passes --percentage_of_cells_to_include_
    transcript to the gene filter (always 0.01); the port honours it.
  * log_total_gene_count = log(n_genes + 1) (a copy of the detected-genes
    covariate); the port uses log(n_counts + 1).
  * muon_creation writes the TSS+1 end into a column named "end " (trailing
    space), so transcript_end is the gene END; the port writes TSS + 1.
  * PerturbLoader computes an enhancer's coordinate, then discards it and gives
    every non-gene element 10 random genes; the port anchors putative
    enhancers on their guide coordinate (controls still get random genes).
  * multiseq.py drops only nUMI (ncol(dt) is NULL in R), so nUMI_total is
    classified as a barcode; the port drops both total columns.
  * multiseq.py filters calls "Double" instead of "Doublet", so doublets are
    kept; the port removes them.

Subcommands
  assay-spec       Task 1: chemistry string + whitelist URL for 10XV2 / 10XV3 /
                   5' PE, or parse a custom kallisto technology string.
  composition      compositionREADS*: per-position nucleotide bias of R1/R2.
  guide-table      guidePreprocessing: guide xlsx/tsv -> guide_features.txt.
  count-guides     creatingGuideRef + mappingGuide: kite-style guide counting
                   from FASTQ (exact + Hamming-1), barcode whitelist correction,
                   UMI collapsing -> <name>_ks_guide_out/counts_unfiltered/adata.h5ad.
  map-rna          downloadReference + mappingscRNA: runs `kb ref` / `kb count`
                   when kb is on PATH, otherwise prints the exact commands.
  manifest         preprocessing: count directories -> lane manifest.
  prefilter        concact_prefiltering: per-lane QC, doublets, intersection,
                   lane concatenation, gene filter.
  annotate         moun_raw_creation: guide coordinates + gene TSS annotation.
  binarize         merge_bin_and_muon: guide merge, binarisation, covariates.
  multiseq         preprocess_bar_multiseq + MultiSeq: deMULTIplex port.
  perturb-loader   PerturbLoaderGeneration: elements, guide membership, cis
                   gene sets, GUIDE_TYPE, test pairs.
  de               runSceptre: SCEPTRE-style NB score test with conditional
                   resampling (python), Mann-Whitney (Task 3 alternative), or
                   the upstream R script (engine sceptre-r, needs Rscript).
  results          create_anndata_from_sceptre: BH, guide x gene layers,
                   Fisher element aggregation.
  tracks           Task 4: BED + pyGenomeTracks links of guide -> gene tests.
  assign-guides    Task guide-assign: depth-aware Poisson-mixture guide calls
                   next to the fixed UMI threshold.
  cellranger-inputs  Task cellranger: library.csv + feature_ref.csv (+ command).
  rename-fastq     rename_and_symlink: bcl2fastq-style names + symlinks.
  config           parse a Nextflow perturb.config into JSON and validate it.
  inspect          Task 2: summarise a modality bundle / .h5mu (targets, guides).
  run              the chain prefilter -> annotate -> binarize -> (multiseq) ->
                   perturb-loader -> de -> results -> tracks from count data
                   (h5ad) or guide FASTQs.
  selftest         synthetic two-lane screen with planted knock-downs, doublets,
                   high-mito cells, MULTI-seq hashes and FASTQ reads.

Output: Docs/BulkCRISPR/<timestamp>_<label>/.  numpy, pandas, scipy, anndata
required; scrublet, mudata, matplotlib optional (a built-in Scrublet-equivalent
and h5ad modality bundles are used without them).

Usage:
    igvfagent bulk-crispr assay-spec --assay 10XV3
    igvfagent bulk-crispr composition --r1 guide_R1.fastq.gz --r2 guide_R2.fastq.gz --label guide_S1_L1
    igvfagent bulk-crispr guide-table --guides df_from_gasperini_tss.xlsx
    igvfagent bulk-crispr count-guides --guides guide_features.txt --r1 R1.fastq.gz --r2 R2.fastq.gz --chemistry 10XV2 --whitelist 737K-august-2016.txt --name S1_L1
    igvfagent bulk-crispr map-rna --fastq R1.fastq.gz R2.fastq.gz --chemistry 10XV2 --whitelist 737K-august-2016.txt --name S1_L1
    igvfagent bulk-crispr run --guide-table guides.xlsx --gtf genes.gtf --rna-dirs S1_L1_ks_transcripts_out --guide-dirs S1_L1_ks_guide_out --label gasperini
    igvfagent bulk-crispr selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "BulkCRISPR"

UPSTREAM_REPO = "IGVF-CRISPR/bulk_crispr_pipeline"
UPSTREAM_COMMIT = "6f4ae09fa7f7dffc1d98b3dc4dc33e73c8ff48ac"
WL_BASE = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/"
ASSAYS = {
    # name: (kallisto technology string, whitelist URL)
    "10XV2": ("0,0,16:0,16,26:1,0,0", WL_BASE + "737K-august-2016.txt"),
    "10XV3": ("0,0,16:0,16,28:1,0,0", WL_BASE + "3M-february-2018.txt.gz"),
    "5PE": ("0,0,16:0,16,26:1,0,0", WL_BASE + "737K-august-2016.txt"),
}
ASSAY_ALIASES = {"10XV2": "10XV2", "10VX2": "10XV2", "10X_V2": "10XV2", "10XV3": "10XV3", "10X_V3": "10XV3",
                 "5'PE": "5PE", "5PE": "5PE", "5P": "5PE", "10X5": "5PE", "10X5'": "5PE"}
GUIDE_TABLE_COLS = ["sgRNA_ID", "Target_name", "sgRNA_sequences", "chr", "start", "end"]
COVARIATES = ["bath_number", "percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count"]
RESULT_COLS = ["gene_id", "gRNA_id", "pair_type", "p_value", "z_value", "log_fold_change"]
NUMBER_FOR_RANDOM_CONTROLS = 10
MIN_CELLS_PER_GUIDE = 30
NUCS = ["A", "C", "T", "G"]
REQUIRED_CONFIG = ["GUIDE_FEATURES", "CHEMISTRY", "FASTQ_FILES_TRANSCRIPTS", "FASTQ_NAMES_TRANSCRIPTS",
                   "FASTQ_FILES_GUIDES", "FASTQ_NAMES_GUIDES"]
CONFIG_DEFAULTS = {"DISTANCE_NEIGHBORS": 1000000, "IN_TRANS": "FALSE", "EXPECTED_CELL_NUMBER": 8000, "MITO_SPECIE": "hsapiens",
                   "MITO_EXPECTED_PERCENTAGE": 0.2, "PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS": 0.01,
                   "TRANSCRIPTS_UMI_TRHESHOLD": 100, "GUIDE_UMI_LIMIT": 5, "MERGE": False, "DIRECTION": "both",
                   "RUN_MULTISEQ": False, "ADDGENENAMES": "", "BAR_MULTI": [1, 16], "UMI_MULTI": [17, 28], "R2_MULTI_TAG": [1, 8]}

BLUE, ORANGE, GREEN, INK, INK2, AXIS = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#d0cfca"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

log = logging.getLogger("bulk_crispr")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"bulk_crispr_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    base = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d, i = base, 1
    while d.exists():
        i += 1
        d = Path(f"{base}_{i}")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy + scipy + anndata: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


def _ad():
    try:
        import anndata  # type: ignore
        return anndata
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs anndata: pip install anndata") from e


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
    ax.set_xlabel(xlabel, fontsize=9, color=INK2)
    ax.set_ylabel(ylabel, fontsize=9, color=INK2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(labelsize=8, colors=INK2)


def _save(fig, path: Path, figs: Optional[list] = None) -> Path:
    plt = _plt()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"Figure: {path}")
    if figs is not None:
        figs.append(str(path))
    return path


def md_table(headers: "List[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 40) -> str:
    rows = list(rows)
    out = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "---|" * len(headers)]
    for r in rows[:max_rows]:
        out.append("| " + " | ".join(_fmt(x) for x in r) + " |")
    if len(rows) > max_rows:
        out.append(f"| ... {len(rows) - max_rows} more rows | " + " | " * (len(headers) - 1))
    return "\n".join(out)


def _fmt(x, nd: int = 3) -> str:
    if isinstance(x, float):
        if x != x:
            return "NA"
        if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
            return f"{x:.2e}"
        return f"{x:.{nd}f}"
    return str(x)


def _truthy(v) -> bool:
    return str(v).strip().strip("'\"").lower() in {"1", "true", "yes", "t", "y"}


def write_tsv(df, path: Path, index: bool = False, header: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, header=header)
    print(f"TSV: {path}")
    return path


def write_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))
    print(f"JSON: {path}")
    return path


def _json_default(o):
    np = _np()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def finish(d: Path, title: str, sections: "List[str]", summary: dict) -> None:
    rep = d / "report.md"
    head = [f"# {title}", "", f"Port of https://github.com/{UPSTREAM_REPO} (commit {UPSTREAM_COMMIT[:10]}).", ""]
    rep.write_text("\n".join(head + sections) + "\n")
    print(f"Report: {rep}")
    write_json(summary, d / "summary.json")


def open_text(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def which(binary: str) -> Optional[str]:
    return shutil.which(binary)


def run_or_print(cmd: "List[str]", cwd: Optional[Path] = None, dry: bool = False) -> int:
    """Run an external binary if it is on PATH; otherwise print the exact command and return 0."""
    shown = " ".join(cmd)
    if dry or not which(cmd[0]):
        why = "dry run" if dry else f"`{cmd[0]}` not on PATH"
        print(f"Command ({why}; not executed): {shown}")
        return 0
    print(f"Running: {shown}")
    return subprocess.call(cmd, cwd=str(cwd) if cwd else None)


def read_table_any(path: Path):
    """xlsx / csv / tsv reader (upstream reads the guide sheet with pandas.read_excel)."""
    pd = _pd()
    path = Path(path)
    s = path.name.lower()
    if s.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    sep = "," if s.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=sep)


# ---------------------------------------------------------------------------
# Modality bundles (MuData stand-in when mudata is not installed)
# ---------------------------------------------------------------------------

def save_bundle(mods: "Dict[str, Any]", path: Path) -> Path:
    """Write {modality: AnnData} as <path>/ <mod>.h5ad (+ <path>.h5mu when mudata is importable)."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    for k, a in mods.items():
        a = _sanitize(a)
        a.write_h5ad(path / f"{k}.h5ad")
    (path / "modalities.json").write_text(json.dumps(list(mods.keys())))
    print(f"Wrote: {path}")
    try:
        import mudata  # type: ignore
        mudata.MuData({k: _sanitize(v) for k, v in mods.items()}).write(str(path) + ".h5mu")
        print(f"Wrote: {path}.h5mu")
    except Exception:
        pass
    return path


def _sanitize(a):
    """h5ad cannot hold mixed-type object columns; stringify them."""
    pd = _pd()
    for df in (a.obs, a.var):
        for c in list(df.columns):
            if df[c].dtype == object:
                df[c] = df[c].astype(str)
        df.columns = [str(c) for c in df.columns]
    a.obs.index = a.obs.index.astype(str)
    a.var.index = a.var.index.astype(str)
    return a


def load_bundle(path: Path) -> "Dict[str, Any]":
    ad = _ad()
    path = Path(path)
    if path.is_dir():
        names = json.loads((path / "modalities.json").read_text()) if (path / "modalities.json").exists() \
            else [p.stem for p in sorted(path.glob("*.h5ad"))]
        return {n: ad.read_h5ad(path / f"{n}.h5ad") for n in names}
    if path.suffix == ".h5mu":
        try:
            import mudata  # type: ignore
            m = mudata.read_h5mu(str(path))
            return {k: m.mod[k] for k in m.mod}
        except ImportError:
            import h5py  # type: ignore
            try:
                from anndata.io import read_elem  # type: ignore
            except ImportError:
                from anndata.experimental import read_elem  # type: ignore
            out = {}
            with h5py.File(path, "r") as f:
                for k in f["mod"].keys():
                    out[k] = read_elem(f["mod"][k])
            return out
    raise SystemExit(f"not a modality bundle (dir of h5ad) or .h5mu: {path}")


def dense(X):
    np = _np()
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def row_sums(X):
    np = _np()
    return np.asarray(X.sum(axis=1)).ravel()


def col_sums(X):
    np = _np()
    return np.asarray(X.sum(axis=0)).ravel()


def nnz_rows(X):
    np = _np()
    return np.asarray((X > 0).sum(axis=1)).ravel()


def nnz_cols(X):
    np = _np()
    return np.asarray((X > 0).sum(axis=0)).ravel()


# ---------------------------------------------------------------------------
# Task 1: assay specification / kallisto technology strings
# ---------------------------------------------------------------------------

def parse_technology(tech: str) -> "Dict[str, List[Tuple[int, int, int]]]":
    """kallisto bc:umi:seq technology string -> {'bc': [(file, start, stop)], 'umi': [...], 'seq': [...]}.

    stop == 0 means "to the end of the read" (kallisto convention)."""
    parts = tech.strip().split(":")
    if len(parts) != 3:
        raise ValueError(f"technology string must be bc:umi:seq, got {tech!r}")
    out = {}
    for name, p in zip(("bc", "umi", "seq"), parts):
        nums = [int(x) for x in p.split(",") if x.strip() != ""]
        if len(nums) % 3:
            raise ValueError(f"{name} part {p!r} is not a list of file,start,stop triplets")
        out[name] = [tuple(nums[i:i + 3]) for i in range(0, len(nums), 3)]
    return out


def resolve_assay(assay: str, chemistry: str = "", whitelist: str = "") -> "Tuple[str, str]":
    """check_assay_spects.py: 'custom' -> (chemistry, whitelist) as given; known assay -> table entry."""
    key = assay.strip().upper().replace(" ", "")
    if key == "CUSTOM":
        if not chemistry:
            raise ValueError("assay 'custom' needs --chemistry (kallisto bc:umi:seq string)")
        parse_technology(chemistry)
        return chemistry, whitelist
    if key in ASSAY_ALIASES:
        return ASSAYS[ASSAY_ALIASES[key]]
    parse_technology(assay)  # a raw technology string is accepted too
    return assay, whitelist


def cmd_assay_spec(args: argparse.Namespace) -> int:
    try:
        chem, wl = resolve_assay(args.assay, args.chemistry or "", args.whitelist or "")
    except ValueError as e:
        print(f"error: {e}")
        return 2
    t = parse_technology(chem)
    print(f"{chem},{wl}")
    print(json.dumps({"assay": args.assay, "chemistry": chem, "whitelist": wl, "layout": t}))
    return 0


def extract(reads: "List[str]", triplets: "List[Tuple[int, int, int]]") -> str:
    return "".join(reads[f][s:(e if e > 0 else None)] for f, s, e in triplets)


# ---------------------------------------------------------------------------
# compositionREADS*: fq_composition.py
# ---------------------------------------------------------------------------

def read_fastq_head_lines(path: Path, n_lines: int = 10000) -> "List[str]":
    out = []
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i >= n_lines:
                break
            out.append(line.rstrip("\n"))
    return out


def composition_sequences(lines: "List[str]", upstream_compat: bool = False) -> "List[str]":
    if upstream_compat:  # upstream filter: drop lines containing '@', '+' or 'F' (quality strings without them leak)
        return [x for x in lines if "@" not in x and "+" not in x and "F" not in x]
    return [lines[i] for i in range(1, len(lines), 4)]


def compositional_bias(seqs: "List[str]"):
    """Per position: counts of A/C/T/G and their sample sd (upstream subtracts n first, which does not change the sd)."""
    pd, np = _pd(), _np()
    seqs = [s for s in seqs if s]
    L = max((len(s) for s in seqs), default=0)
    rows = []
    for pos in range(L):
        col = [s[pos] for s in seqs if len(s) > pos]
        cnt = {n: sum(1 for c in col if c == n) for n in NUCS}
        vals = np.array([cnt[n] - len(seqs) for n in NUCS], float)
        rows.append({"position": pos, **cnt, "n_reads": len(col), "sd_bias": float(np.std(vals, ddof=1))})
    return pd.DataFrame(rows, columns=["position"] + NUCS + ["n_reads", "sd_bias"])


def cmd_composition(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    figs: list = []
    summ = {"reads": {}}
    secs = ["Per-position compositional bias (sd between the four nucleotide counts) of the first "
            f"{args.n_lines} FASTQ lines, as fq_composition.py.", ""]
    plt = None if args.no_plots else _plt()
    for tag, fq in (("R1", args.r1), ("R2", args.r2)):
        if not fq:
            continue
        seqs = composition_sequences(read_fastq_head_lines(Path(fq), args.n_lines), args.upstream_compat)
        df = compositional_bias(seqs)
        name = f"{safe_label(args.label)}_{tag}"
        write_tsv(df, d / f"{name}_composition.tsv")
        summ["reads"][tag] = {"path": fq, "n_sequences": len(seqs), "length": int(len(df)),
                              "mean_sd_bias": float(df["sd_bias"].mean()) if len(df) else None,
                              "max_sd_bias_position": int(df["sd_bias"].idxmax()) if len(df) else None}
        if plt is not None and len(df):
            fig, ax = plt.subplots(figsize=(12, 2.4))
            ax.plot(df["position"], df["sd_bias"], color=BLUE, lw=1.6)
            ax.set_xlim(0, len(df))
            _style(ax, f"{name}: compositional bias (sd between the 4 nucleotides)", "read position", "sd")
            _save(fig, d / f"{name}_composition.png", figs)
        secs.append(f"- {tag}: {len(seqs)} sequences, length {len(df)}, mean sd {summ['reads'][tag]['mean_sd_bias']:.1f}")
    summ["figures"] = figs
    summ["upstream_compat"] = bool(args.upstream_compat)
    finish(d, "Read composition", secs, summ)
    return 0


# ---------------------------------------------------------------------------
# guidePreprocessing: guide_table_processing.py
# ---------------------------------------------------------------------------

def process_guide_table(df):
    """Target_name -> '<t>|<n>', pipeline_id = '<t>|<n>_sgrna_<chr>:<start>:<end>' (groups sorted by target)."""
    pd = _pd()
    miss = [c for c in ["Target_name", "sgRNA_sequences", "chr", "start", "end"] if c not in df.columns]
    if miss:
        raise SystemExit(f"guide table is missing columns {miss} (expected {GUIDE_TABLE_COLS})")
    parts = []
    for k, v in df.groupby("Target_name", sort=True):
        v = v.copy()
        v["Target_name"] = [f"{k}|{n + 1}" for n in range(len(v))]
        v["pipeline_id"] = [f"{t}_sgrna_{c}:{s}:{e}" for t, c, s, e in zip(v["Target_name"], v["chr"], v["start"], v["end"])]
        parts.append(v)
    return pd.concat(parts)


def cmd_guide_table(args: argparse.Namespace) -> int:
    df = read_table_any(Path(args.guides))
    out = process_guide_table(df)
    d = Path(args.out_dir) if args.out_dir else run_dir(args.label)
    d.mkdir(parents=True, exist_ok=True)
    feat = d / "guide_features.txt"
    out[["sgRNA_sequences", "pipeline_id"]].to_csv(feat, sep="\t", header=False, index=False)
    print(f"Wrote: {feat}")
    write_tsv(out, d / "guide_table_processed.tsv")
    targets = out["Target_name"].str.split("|").str[0]
    summ = {"n_guides": int(len(out)), "n_targets": int(targets.nunique()),
            "guides_per_target": targets.value_counts().to_dict(), "guide_features": str(feat)}
    finish(d, "Guide table", ["Guides targeting the same element share a Target_name; the pipeline numbers them "
                              "<target>|1, <target>|2, ...", "",
                              md_table(["pipeline_id", "sequence"], out[["pipeline_id", "sgRNA_sequences"]].values, 30)], summ)
    return 0


def read_guide_features(path: Path) -> "List[Tuple[str, str]]":
    """guide_features.txt (seq \\t id, no header) or a raw guide table (processed on the fly)."""
    path = Path(path)
    if path.name.lower().endswith((".xlsx", ".xls")) or "Target_name" in path.open("rt", errors="ignore").readline():
        t = process_guide_table(read_table_any(path))
        return [(str(s).upper(), str(i)) for s, i in zip(t["sgRNA_sequences"], t["pipeline_id"])]
    out = []
    with open_text(path) as fh:
        for line in fh:
            if line.strip():
                s, i = line.rstrip("\n").split("\t")[:2]
                out.append((s.strip().upper(), i.strip()))
    return out


# ---------------------------------------------------------------------------
# creatingGuideRef + mappingGuide: kite-style counting
# ---------------------------------------------------------------------------

def hamming1_variants(seq: str) -> "Iterable[str]":
    for i, c in enumerate(seq):
        for b in "ACGT":
            if b != c:
                yield seq[:i] + b + seq[i + 1:]


def build_kite_map(features: "List[Tuple[str, str]]") -> "Tuple[Dict[str, int], Dict[str, int], List[str], dict]":
    """Exact map + Hamming-1 map; variants shared by two guides (or equal to another guide) are removed (kite)."""
    exact: Dict[str, int] = {}
    ids: List[str] = []
    for s, i in features:
        if s in exact:
            raise SystemExit(f"duplicate guide sequence {s} ({ids[exact[s]]} and {i})")
        exact[s] = len(ids)
        ids.append(i)
    mm: Dict[str, int] = {}
    collide = set()
    for s, j in exact.items():
        for v in hamming1_variants(s):
            if v in exact:
                collide.add(v)
                continue
            if v in mm and mm[v] != j:
                collide.add(v)
            else:
                mm[v] = j
    for v in collide:
        mm.pop(v, None)
    return exact, mm, ids, {"n_features": len(ids), "n_mismatch_variants": len(mm), "n_collisions_removed": len(collide)}


def read_whitelist(path: Optional[str]) -> "Optional[set]":
    if not path:
        return None
    with open_text(Path(path)) as fh:
        return {l.strip().split()[0] for l in fh if l.strip()}


def correct_barcode(bc: str, wl: "Optional[set]", cache: dict) -> "Tuple[Optional[str], int]":
    """bustools correct: exact whitelist hit, else the unique whitelist barcode at Hamming 1. Returns (bc, 0 exact / 1 corrected)."""
    if wl is None:
        return bc, 0
    if bc in wl:
        return bc, 0
    if bc in cache:
        return cache[bc]
    hits = [v for v in hamming1_variants(bc) if v in wl]
    res = (hits[0], 1) if len(hits) == 1 else (None, 0)
    cache[bc] = res
    return res


def iter_fastq_pairs(paths: "List[Path]"):
    fhs = [open_text(p) for p in paths]
    try:
        while True:
            recs = []
            for fh in fhs:
                h = fh.readline()
                if not h:
                    return
                s = fh.readline().rstrip("\n")
                fh.readline()
                fh.readline()
                recs.append(s)
            yield recs
    finally:
        for fh in fhs:
            fh.close()


def match_guide(seq: str, exact: dict, mm: dict, lengths: "List[int]") -> "Tuple[Optional[int], int]":
    """Scan the sequence for a guide: exact windows first, then Hamming-1 windows. Ambiguous -> None."""
    for L in lengths:
        hits = {exact[seq[o:o + L]] for o in range(0, len(seq) - L + 1) if seq[o:o + L] in exact}
        if len(hits) == 1:
            return hits.pop(), 0
        if len(hits) > 1:
            return None, -1
    for L in lengths:
        hits = {mm[seq[o:o + L]] for o in range(0, len(seq) - L + 1) if seq[o:o + L] in mm}
        if len(hits) == 1:
            return hits.pop(), 1
        if len(hits) > 1:
            return None, -1
    return None, -2


def count_guides(features, fastqs: "List[Path]", tech: str, whitelist: "Optional[set]" = None, max_reads: int = 0):
    """-> (AnnData cells x guides of UMI counts, stats dict)."""
    pd, np, ad = _pd(), _np(), _ad()
    from scipy import sparse  # type: ignore
    t = parse_technology(tech)
    exact, mm, ids, kstats = build_kite_map(features)
    lengths = sorted({len(s) for s in exact}, reverse=True)
    nfiles = 1 + max(f for part in t.values() for f, _, _ in part)
    if len(fastqs) < nfiles:
        raise SystemExit(f"technology {tech} needs {nfiles} FASTQ files, got {len(fastqs)}")
    st = {"reads": 0, "mapped_exact": 0, "mapped_mismatch": 0, "ambiguous": 0, "unmapped": 0,
          "bc_exact": 0, "bc_corrected": 0, "bc_not_in_whitelist": 0, "umis_multi_guide_dropped": 0}
    umi_feat: Dict[Tuple[str, str], set] = {}
    cache: dict = {}
    for recs in iter_fastq_pairs(fastqs[:nfiles]):
        st["reads"] += 1
        if max_reads and st["reads"] > max_reads:
            st["reads"] -= 1
            break
        g, how = match_guide(extract(recs, t["seq"]).upper(), exact, mm, lengths)
        if g is None:
            st["ambiguous" if how == -1 else "unmapped"] += 1
            continue
        st["mapped_exact" if how == 0 else "mapped_mismatch"] += 1
        bc, corr = correct_barcode(extract(recs, t["bc"]), whitelist, cache)
        if bc is None:
            st["bc_not_in_whitelist"] += 1
            continue
        st["bc_corrected" if corr else "bc_exact"] += 1
        umi_feat.setdefault((bc, extract(recs, t["umi"])), set()).add(g)
    counts: Dict[Tuple[str, int], int] = {}
    for (bc, _u), gs in umi_feat.items():
        if len(gs) > 1:
            st["umis_multi_guide_dropped"] += 1
            continue
        k = (bc, next(iter(gs)))
        counts[k] = counts.get(k, 0) + 1
    bcs = sorted({k[0] for k in counts})
    bi = {b: i for i, b in enumerate(bcs)}
    rows = [bi[b] for b, _ in counts]
    cols = [g for _, g in counts]
    X = sparse.csr_matrix((np.array(list(counts.values()), dtype=np.float32), (rows, cols)), shape=(len(bcs), len(ids)))
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=bcs), var=pd.DataFrame({"feature_name": ids}, index=ids))
    st.update(kstats)
    st["cells_with_guides"] = len(bcs)
    st["umis"] = int(sum(counts.values()))
    st["pct_mapped"] = 100.0 * (st["mapped_exact"] + st["mapped_mismatch"]) / max(st["reads"], 1)
    return a, st


def cmd_count_guides(args: argparse.Namespace) -> int:
    chem, wl_url = resolve_assay(args.chemistry, args.chemistry, args.whitelist or "")
    wl_path = args.whitelist or ""
    feats = read_guide_features(Path(args.guides))
    base = Path(args.out_dir) if args.out_dir else run_dir(args.label or f"count_guides_{args.name}")
    out = base / f"{args.name}_ks_guide_out"
    (out / "counts_unfiltered").mkdir(parents=True, exist_ok=True)
    fastqs = [Path(p) for p in (args.fastq or [])] or [Path(args.r1), Path(args.r2)]
    a, st = count_guides(feats, fastqs, chem, read_whitelist(wl_path), args.max_reads)
    a.write_h5ad(out / "counts_unfiltered" / "adata.h5ad")
    print(f"Wrote: {out / 'counts_unfiltered' / 'adata.h5ad'}")
    write_json(st, out / "kite_stats.json")
    kb_cmds = [f"kb ref -i guide_index.idx -f1 genome.fa.gz -g t2guide.txt --workflow kite guide_features.txt",
               f"kb count -i guide_index.idx -g t2guide.txt --verbose --report --workflow kite -w {wl_path or '<whitelist>'} "
               f"--h5ad -x {chem} -o {args.name}_ks_guide_out -t <cpus> {' '.join(str(p) for p in fastqs)} --overwrite -m 48G"]
    secs = [f"Guide counts for {args.name}: {st['reads']} reads, {st['pct_mapped']:.1f}% mapped "
            f"({st['mapped_exact']} exact, {st['mapped_mismatch']} at Hamming 1), {st['cells_with_guides']} barcodes, {st['umis']} UMIs.",
            "", "Upstream equivalent (kallisto|bustools kite):", "", "```bash", *kb_cmds, "```"]
    st["counts"] = str(out / "counts_unfiltered" / "adata.h5ad")
    st["whitelist_url_for_assay"] = wl_url
    finish(base, f"Guide counting {args.name}", secs, st)
    return 0


# ---------------------------------------------------------------------------
# downloadReference + mappingscRNA (binary wrappers)
# ---------------------------------------------------------------------------

def cmd_map_rna(args: argparse.Namespace) -> int:
    chem, wl_url = resolve_assay(args.chemistry, args.chemistry, args.whitelist or "")
    base = Path(args.out_dir) if args.out_dir else run_dir(args.label or f"map_rna_{args.name}")
    base.mkdir(parents=True, exist_ok=True)
    idx, t2g = args.index or str(base / "transcriptome_index.idx"), args.t2g or str(base / "transcriptome_t2g.txt")
    rc = 0
    cmds = []
    if not args.index:
        cmds.append(["kb", "ref", "-d", args.reference, "-i", idx, "-g", t2g])
    out = base / f"{args.name}_ks_transcripts_out"
    cmds.append(["kb", "count", "-i", idx, "-g", t2g, "--verbose", "--workflow", "kite", "-w", args.whitelist or wl_url,
                 "--h5ad", "-x", chem, "-o", str(out), "-t", str(args.threads), *args.fastq, "--overwrite", "-m", args.memory])
    # (upstream passes --workflow kite to the cDNA count too; kept for faithfulness, override with --workflow)
    if args.workflow:
        cmds[-1][cmds[-1].index("kite")] = args.workflow
    for c in cmds:
        rc = rc or run_or_print(c, dry=args.dry_run)
    ran = bool(which("kb")) and not args.dry_run
    summ = {"commands": [" ".join(c) for c in cmds], "executed": ran, "output": str(out)}
    finish(base, f"scRNA mapping {args.name}", ["kallisto|bustools commands" + (" (executed)" if ran else " (kb not on PATH: printed, not run)"),
                                                 "", "```bash", *summ["commands"], "```"], summ)
    return rc


# ---------------------------------------------------------------------------
# preprocessing.py: lane manifest
# ---------------------------------------------------------------------------

def lane_manifest(dirs: "List[str]"):
    pd = _pd()
    paths = [str(Path(p)) + "/counts_unfiltered/adata.h5ad" for p in dirs]
    df = pd.DataFrame({"file_path": paths})
    df["sample"] = ["Guide" if "guide" in p else "scRNA" for p in paths]
    df["lane"] = [p.split("/")[-3].split("_")[1].replace("L", "") for p in paths]
    return df


def cmd_manifest(args: argparse.Namespace) -> int:
    df = lane_manifest(args.dirs)
    out = Path(args.out) if args.out else run_dir(args.label) / "initial_preprocessing_file_names.txt"
    write_tsv(df, out)
    return 0


# ---------------------------------------------------------------------------
# Doublets: Scrublet (library if installed, else a faithful re-implementation of its defaults)
# ---------------------------------------------------------------------------

def threshold_minimum(values, nbins: int = 256, max_iter: int = 10000) -> float:
    """skimage.filters.threshold_minimum: smooth the histogram until it has two maxima; the minimum between them."""
    np = _np()
    hist, edges = np.histogram(values, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    smooth = hist.astype(float)

    def maxima(h):
        d = np.diff(h)
        out = []
        direction = 0
        for i, x in enumerate(d):
            if x > 0:
                direction = 1
            elif x < 0:
                if direction == 1:
                    out.append(i)
                direction = -1
        return out
    for _ in range(max_iter):
        smooth = np.convolve(np.concatenate([[smooth[0]], smooth, [smooth[-1]]]), np.ones(3) / 3, mode="valid")
        mx = maxima(smooth)
        if len(mx) == 2:
            lo, hi = mx
            return float(centers[lo + int(np.argmin(smooth[lo:hi + 1]))])
        if len(mx) < 2:
            break
    return float("nan")


def scrublet_scores(X, expected_doublet_rate: float = 0.1, sim_doublet_ratio: float = 2.0, n_prin_comps: int = 30,
                    min_counts: int = 3, min_cells: int = 3, min_gene_variability_pctl: float = 85.0, seed: int = 0):
    """Scrublet.scrub_doublets() with its defaults -> (scores, predicted, threshold, simulated scores)."""
    np = _np()
    try:
        import scrublet as scr  # type: ignore
        s = scr.Scrublet(X, expected_doublet_rate=expected_doublet_rate, sim_doublet_ratio=sim_doublet_ratio)
        sc_, pred = s.scrub_doublets(verbose=False)
        return np.asarray(sc_), np.asarray(pred, bool), float(s.threshold_), np.asarray(s.doublet_scores_sim_), "scrublet"
    except Exception:
        pass
    from sklearn.decomposition import PCA  # type: ignore
    from sklearn.neighbors import NearestNeighbors  # type: ignore
    rng = np.random.default_rng(seed)
    E = dense(X).astype(float)
    n_obs = E.shape[0]
    tot = E.sum(1)
    tot[tot == 0] = 1
    target = tot.mean()
    En = E / tot[:, None] * target
    keep = ((E >= min_counts).sum(0) >= min_cells)
    mu = En[:, keep].mean(0)
    var = En[:, keep].var(0)
    vs = np.where(mu > 0, var / np.maximum(mu, 1e-12), 0)  # Fano factor stands in for Scrublet's fitted V-score
    cut = np.percentile(vs, min_gene_variability_pctl) if len(vs) else 0
    genes = np.where(keep)[0][vs >= cut]
    if len(genes) < 2:
        genes = np.where(keep)[0]
    n_sim = int(n_obs * sim_doublet_ratio)
    p1, p2 = rng.integers(0, n_obs, n_sim), rng.integers(0, n_obs, n_sim)
    Es = E[p1] + E[p2]
    ts = Es.sum(1)
    ts[ts == 0] = 1
    Esn = Es / ts[:, None] * target
    A, B = En[:, genes], Esn[:, genes]
    m, sd = A.mean(0), A.std(0)
    sd[sd == 0] = 1
    A, B = (A - m) / sd, (B - m) / sd
    k_pc = max(1, min(n_prin_comps, A.shape[1] - 1, n_obs - 1))
    pca = PCA(n_components=k_pc, random_state=seed, svd_solver="full").fit(A)
    Pa, Pb = pca.transform(A), pca.transform(B)
    P = np.vstack([Pa, Pb])
    is_sim = np.r_[np.zeros(n_obs, bool), np.ones(n_sim, bool)]
    r = n_sim / n_obs
    k = int(round(0.5 * math.sqrt(n_obs)))
    k_adj = max(1, int(round(k * (1 + r))))
    nn = NearestNeighbors(n_neighbors=min(k_adj + 1, len(P))).fit(P)
    idx = nn.kneighbors(P, return_distance=False)[:, 1:]
    kd = is_sim[idx].sum(1)
    rho = expected_doublet_rate
    q = (kd + 1) / (k_adj + 2)
    Ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r))
    obs_s, sim_s = Ld[:n_obs], Ld[n_obs:]
    thr = threshold_minimum(sim_s)
    if thr != thr:
        thr = float(np.percentile(sim_s, 100 * (1 - rho)))
    return obs_s, obs_s > thr, thr, sim_s, "builtin"


# ---------------------------------------------------------------------------
# concact_prefiltering: per-lane QC
# ---------------------------------------------------------------------------

def analyze_batch(a_rna, a_guide, batch, expected_cells: int, mito_pct: float, umi_threshold: int, mito_prefix: str,
                  fig_dir: Optional[Path], figs: list, plots: bool, seed: int = 0):
    """analyze_batch() of concact_and_pre_filtering.py -> (rna, guide, lane stats)."""
    np, pd = _np(), _pd()
    from scipy import sparse  # type: ignore
    a = a_rna.copy()
    if not sparse.issparse(a.X):
        a.X = sparse.csr_matrix(a.X)
    st = {"lane": str(batch), "barcodes_in": int(a.n_obs)}
    knee = np.sort(row_sums(a.X))[::-1]
    above = np.arange(len(knee))[knee > umi_threshold]
    st["knee_cells_above_umi_threshold"] = int(above[-1]) if len(above) else 0
    plt = _plt() if plots else None
    if plt is not None and fig_dir is not None:
        fig_dir.mkdir(parents=True, exist_ok=True)
        from sklearn.decomposition import TruncatedSVD  # type: ignore
        if a.n_obs > 2 and a.n_vars > 2:
            Xs = TruncatedSVD(n_components=2, random_state=seed).fit_transform(a.X)
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.scatter(Xs[:, 0], Xs[:, 1], s=4, alpha=0.5, color=GREEN)
            _style(ax, f"lane {batch}: truncated SVD", "SVD1", "SVD2")
            _save(fig, fig_dir / "svd_batch.png", figs)
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(row_sums(a.X) + 0.5, nnz_rows(a.X) + 0.5, s=4, alpha=0.3, color=GREEN)
        ax.set_xscale("log"); ax.set_yscale("log")
        _style(ax, f"lane {batch}: genes detected vs UMIs", "UMI counts", "genes detected")
        _save(fig, fig_dir / "saturation_batch.png", figs)
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.loglog(np.maximum(knee, 0.5), np.arange(len(knee)) + 1, lw=3, color=GREEN)
        ax.axvline(umi_threshold, color=INK, lw=1.5)
        ax.axhline(max(st["knee_cells_above_umi_threshold"], 1), color=INK, lw=1.5)
        _style(ax, f"lane {batch}: knee plot", "UMI counts", "set of barcodes")
        _save(fig, fig_dir / "knee_plot_batch.png", figs)
    # sc.pp.filter_cells(min_genes=cutoff) -- upstream uses the UMI threshold as the minimum GENE count
    ng = nnz_rows(a.X)
    a = a[ng >= umi_threshold].copy()
    a.obs["n_genes"] = ng[ng >= umi_threshold]
    st["after_min_genes"] = int(a.n_obs)
    ki = min(expected_cells, len(knee) - 1)
    if expected_cells >= len(knee):
        log.warning("EXPECTED_CELL_NUMBER %s >= barcodes %s: using the last knee value (upstream raises IndexError)", expected_cells, len(knee))
    min_counts = float(knee[ki]) if len(knee) else 0.0
    nc = row_sums(a.X)
    a = a[nc >= min_counts].copy()
    a.obs["n_counts"] = nc[nc >= min_counts]
    st["knee_min_counts"] = min_counts
    st["after_knee"] = int(a.n_obs)
    a.var.index = [str(x).split(".")[0] for x in a.var.index]
    fname = a.var["feature_name"].astype(str) if "feature_name" in a.var else pd.Series(a.var.index, index=a.var.index)
    mito = np.asarray(fname.str.startswith(mito_prefix))
    tot = row_sums(a.X)
    a.obs["percent_mito"] = row_sums(a.X[:, mito]) / np.maximum(tot, 1e-12)
    a.obs["n_counts"] = tot
    st["n_mito_genes"] = int(mito.sum())
    if plt is not None and fig_dir is not None:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(a.obs["n_counts"], a.obs["percent_mito"], s=4, alpha=0.4, color=BLUE)
        ax.axhline(mito_pct, color=ORANGE, lw=1)
        _style(ax, f"lane {batch}: mitochondrial fraction", "n_counts", "percent_mito")
        _save(fig, fig_dir / "mito_scatter_batch.png", figs)
    a = a[a.obs["percent_mito"].values < mito_pct].copy()
    st["after_mito"] = int(a.n_obs)
    if plt is not None and fig_dir is not None and a.n_obs:
        fig, axs = plt.subplots(1, 3, figsize=(9, 3))
        for ax, c in zip(axs, ["n_genes", "n_counts", "percent_mito"]):
            ax.violinplot(a.obs[c].astype(float).values, showmedians=True)
            _style(ax, c)
        _save(fig, fig_dir / "box_plot_batch.png", figs)
        frac = np.asarray(a.X.sum(0)).ravel() / max(float(a.X.sum()), 1)
        top = np.argsort(frac)[::-1][:20]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.barh(range(len(top))[::-1], frac[top] * 100, color=BLUE)
        ax.set_yticks(range(len(top))[::-1])
        ax.set_yticklabels(np.asarray(a.var["feature_name"] if "feature_name" in a.var else a.var.index)[top], fontsize=7)
        _style(ax, f"lane {batch}: top 20 genes", "% of counts")
        _save(fig, fig_dir / "top_genes_batch.png", figs)
    if a.n_obs >= 10:
        scores, pred, thr, sim, engine = scrublet_scores(a.X, seed=seed)
    else:
        scores, pred, thr, sim, engine = np.zeros(a.n_obs), np.zeros(a.n_obs, bool), float("nan"), np.zeros(0), "skipped"
    a.obs["doublet_scores"] = scores
    a.obs["predicted_doublets"] = pred
    a.obs["doublet_info"] = a.obs["predicted_doublets"].astype(str)
    st.update({"doublet_engine": engine, "doublet_threshold": thr, "predicted_doublets": int(pred.sum())})
    if plt is not None and fig_dir is not None and len(sim):
        fig, axs = plt.subplots(1, 2, figsize=(8, 3))
        axs[0].hist(scores, bins=40, color=BLUE); axs[0].axvline(thr, color=ORANGE); _style(axs[0], "observed transcriptomes", "doublet score")
        axs[1].hist(sim, bins=40, color=INK2); axs[1].axvline(thr, color=ORANGE); _style(axs[1], "simulated doublets", "doublet score")
        _save(fig, fig_dir / "doublets_batch.png", figs)
    a = a[~a.obs["predicted_doublets"].values].copy()
    st["after_doublets"] = int(a.n_obs)
    g = a_guide.copy()
    g.obs["number_of_nonzero_guides"] = nnz_rows(g.X)
    if plt is not None and fig_dir is not None and g.n_obs:
        fig, ax = plt.subplots(figsize=(3.5, 3))
        ax.violinplot(g.obs["number_of_nonzero_guides"].astype(float).values, showmedians=True)
        _style(ax, f"lane {batch}: guides per barcode", "", "number_of_nonzero_guides")
        _save(fig, fig_dir / "guide_non_zero_batch.png", figs)
    shared = set(g.obs.index).intersection(a.obs.index)
    a = a[[c in shared for c in a.obs.index]].copy()
    g = g[list(a.obs.index)].copy()  # same barcode order in both modalities (MuData aligns them by obs name)
    a.obs["batch_number"] = str(batch)
    g.obs["batch_number"] = str(batch)
    st["cells_shared_with_guides"] = int(a.n_obs)
    return a, g, st


def read_count_h5ad(path: str):
    return _ad().read_h5ad(path)


def cmd_prefilter(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd, np, ad = _pd(), _np(), _ad()
    own = d is None
    d = d or run_dir(args.label)
    if args.manifest:
        man = pd.read_csv(args.manifest, sep="\t", dtype=str)
    else:
        man = lane_manifest(list(args.rna_dirs) + list(args.guide_dirs))
    figs: list = []
    rnas, guides, lanes = [], [], []
    for lane, v in man.groupby("lane", sort=True):
        fg = v[v["sample"] == "Guide"]["file_path"].values
        fr = v[v["sample"] == "scRNA"]["file_path"].values
        if not len(fg) or not len(fr):
            raise SystemExit(f"lane {lane}: need one Guide and one scRNA count file, got {list(v['file_path'])}")
        r, g, st = analyze_batch(read_count_h5ad(fr[0]), read_count_h5ad(fg[0]), lane, args.expected_cell_number,
                                 args.mito_expected_percentage, args.transcripts_umi_threshold, args.mito_prefix,
                                 d / "results_per_lane" / f"lane_{lane}", figs, not args.no_plots, args.seed)
        rnas.append(r); guides.append(g); lanes.append(st)
    rna = ad.concat(rnas, merge="same")
    rna.obs_names_make_unique()
    gd = ad.concat(guides, merge="same")
    gd.obs_names_make_unique()
    pct = 0.01 if args.upstream_compat else args.percentage_of_cells_to_include_transcript
    min_cells = int(rna.n_obs * pct)
    keep = nnz_cols(rna.X) >= min_cells
    rna = rna[:, keep].copy()
    rna.var["n_cells"] = nnz_cols(rna.X)
    out = d / "results_per_lane"
    out.mkdir(parents=True, exist_ok=True)
    _sanitize(rna).write_h5ad(out / "full_raw_scrna_ann_data.h5ad")
    _sanitize(gd).write_h5ad(out / "full_raw_guide_ann_data.h5ad")
    print(f"Wrote: {out / 'full_raw_scrna_ann_data.h5ad'}")
    print(f"Wrote: {out / 'full_raw_guide_ann_data.h5ad'}")
    lt = pd.DataFrame(lanes)
    write_tsv(lt, d / "lane_qc.tsv")
    summ = {"lanes": lanes, "cells": int(rna.n_obs), "genes": int(rna.n_vars), "genes_min_cells": min_cells,
            "gene_filter_fraction": pct, "scrna": str(out / "full_raw_scrna_ann_data.h5ad"),
            "guides": str(out / "full_raw_guide_ann_data.h5ad"), "figures": figs}
    if own:
        finish(d, "Per-lane QC and pre-filtering",
               ["Knee / min-genes / mito / Scrublet filtering per lane, barcode intersection with the guide library, "
                "lane concatenation and the low-expression gene filter.", "",
                md_table(list(lt.columns), lt.values)], summ)
    args._prefilter = summ
    return 0


# ---------------------------------------------------------------------------
# moun_raw_creation: muon_creation.py
# ---------------------------------------------------------------------------

def read_gtf_genes(path: Path):
    """GTF 'gene' records -> chr, start, end, strand, gene_id (no version), gene_name, tss (start if + else end)."""
    pd = _pd()
    rows = []
    rx_id, rx_nm = re.compile(r'gene_id "([^"]+)"'), re.compile(r'gene_name "([^"]+)"')
    with open_text(Path(path)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = rx_id.search(f[8])
            gnm = rx_nm.search(f[8])
            if not gid:
                continue
            rows.append((f[0], int(f[3]), int(f[4]), f[6], gid.group(1).split(".")[0], gnm.group(1) if gnm else gid.group(1)))
    df = pd.DataFrame(rows, columns=["chr", "start", "end", "strand", "gene_id", "gene_name"])
    df["gene_real_start"] = [s if st == "+" else e for s, e, st in zip(df["start"], df["end"], df["strand"])]
    return df


def parse_guide_id(k: str) -> dict:
    after = k.split("|")[1]
    return {"guide_chr": after.split(":")[0].split("_")[-1], "guide_end": after.split(":")[-1],
            "guide_start": after.split(":")[-2], "guide_number": after.split(":")[0].split("_")[0],
            "target_elements": k.split("|")[0]}


def annotate_modalities(guide, rna, genes, upstream_compat: bool = False):
    pd = _pd()
    g = guide.copy()
    parsed = pd.DataFrame([parse_guide_id(k) for k in g.var.index], index=g.var.index)
    for c in parsed.columns:
        g.var[c] = parsed[c].values
    r = rna.copy()
    gi = genes.drop_duplicates("gene_id").set_index("gene_id")
    chrs, starts, ends = [], [], []
    for x in r.var.index:
        if x in gi.index:
            row = gi.loc[x]
            chrs.append(str(row["chr"]))
            starts.append(int(row["gene_real_start"]))
            ends.append(int(row["end"]) if upstream_compat else int(row["gene_real_start"]) + 1)
        else:
            chrs.append("NOT_FOUND"); starts.append(0); ends.append(0)
    r.var["transcript_chr"], r.var["transcript_start"], r.var["transcript_end"] = chrs, starts, ends
    return g, r


def cmd_annotate(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    ad = _ad()
    own = d is None
    d = d or run_dir(args.label)
    genes = read_gtf_genes(Path(args.gtf))
    g, r = annotate_modalities(ad.read_h5ad(args.ann_guide), ad.read_h5ad(args.ann_exp), genes, args.upstream_compat)
    out = save_bundle({"guides": g, "scRNA": r}, d / "raw_mudata_guide_and_transcripts")
    nf = int((r.var["transcript_chr"] == "NOT_FOUND").sum())
    summ = {"bundle": str(out), "guides": int(g.n_vars), "genes": int(r.n_vars), "genes_not_found_in_gtf": nf,
            "targets": sorted(set(g.var["target_elements"]))}
    if own:
        finish(d, "Guide and gene annotation", [f"{g.n_vars} guides parsed from pipeline_id; {r.n_vars - nf} of {r.n_vars} genes "
                                                 "placed at their TSS from the GTF gene records."], summ)
    args._annotate = summ
    return 0


# ---------------------------------------------------------------------------
# merge_bin_and_muon
# ---------------------------------------------------------------------------

def merge_guides(guide, merge: bool):
    """Sum guide UMI columns by feature_name (or by target before '|' with merge). Groups sorted."""
    pd, np = _pd(), _np()
    names = [str(f) for f in (guide.var["feature_name"] if "feature_name" in guide.var else guide.var.index)]
    keys = [n.split("|")[0] for n in names] if merge else names
    X = dense(guide.X)
    uk = sorted(set(keys))
    M = np.zeros((X.shape[0], len(uk)))
    pos = {k: i for i, k in enumerate(uk)}
    for j, k in enumerate(keys):
        M[:, pos[k]] += X[:, j]
    return pd.DataFrame(M, index=guide.obs.index, columns=uk)


def make_covariates(rna, guide_binary, upstream_compat: bool = False):
    pd, np = _pd(), _np()
    cov = pd.DataFrame({"bath_number": rna.obs["batch_number"].astype(str).values}, index=rna.obs.index)
    cov["percent_mito"] = rna.obs["percent_mito"].astype(float).values
    ng = rna.obs["n_genes"].astype(float).values
    cov["log_number_of_detected_genes"] = np.log(ng)
    nc = rna.obs["n_counts"].astype(float).values if "n_counts" in rna.obs else row_sums(rna.X)
    cov["log_total_gene_count"] = np.log(ng + 1) if upstream_compat else np.log(nc + 1)
    cov["log_total_guide_count"] = np.log(guide_binary.sum(1).values + 1)
    return cov


def binarize_modalities(mods, guide_umi_limit: int = 5, merge: bool = False, upstream_compat: bool = False):
    pd, np, ad = _pd(), _np(), _ad()
    from scipy import sparse  # type: ignore
    g0, r0 = mods["guides"], mods["scRNA"]
    merged = merge_guides(g0, merge)
    if set(merged.index) >= set(r0.obs.index):
        merged = merged.loc[list(r0.obs.index)]  # align guide rows to the scRNA cells
    binary = (merged > guide_umi_limit).astype(int)
    cov = make_covariates(r0, binary, upstream_compat)
    names = r0.var["feature_name"].astype(str).values if "feature_name" in r0.var else r0.var.index.values
    first = ~pd.Index(names).duplicated()
    rv = r0.var.copy()
    rv["ensg"] = rv.index.values
    rv.index = names
    rv = rv[first]
    X = r0.X[:, np.where(first)[0]]
    r = ad.AnnData(X=sparse.csr_matrix(X) if not sparse.issparse(X) else X.tocsr(), obs=cov.copy(), var=rv)
    gv = g0.var.copy()
    gv.index = [str(x) for x in (g0.var["feature_name"] if "feature_name" in g0.var else g0.var.index)]
    if merge:
        last = {x.split("|")[0]: x for x in gv.index}
        cols = [last[c] for c in binary.columns]
    else:
        cols = list(binary.columns)
    gvar = gv.loc[cols]
    gvar = gvar[~gvar.index.duplicated()]
    g = ad.AnnData(X=binary.values.astype(np.float32), obs=cov.copy(), var=gvar)
    g.layers["counts"] = merged.values.astype(np.float32)
    return {"guides": g, "scRNA": r}, cov


def cmd_binarize(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd = _pd()
    own = d is None
    d = d or run_dir(args.label)
    mods = load_bundle(Path(args.muon_data))
    out_mods, cov = binarize_modalities(mods, args.guide_umi_limit, _truthy(args.merge), args.upstream_compat)
    figs: list = []
    plt = None if args.no_plots else _plt()
    if plt is not None:
        fig, ax = plt.subplots(figsize=(5, 3.5))
        lv = sorted(cov["bath_number"].unique())
        ax.boxplot([cov.loc[cov["bath_number"] == b, "log_number_of_detected_genes"].values for b in lv], labels=lv)
        _style(ax, "Variation of log(number detected genes) among batches", "bath_number", "log_number_of_detected_genes")
        _save(fig, d / "batch_effect_visualization.png", figs)
    out = save_bundle(out_mods, d / "processed_mudata_guide_and_transcripts")
    g = out_mods["guides"]
    per_guide = pd.Series(g.X.sum(0), index=g.var.index)
    write_tsv(cov, d / "covariates.tsv", index=True)
    summ = {"bundle": str(out), "cells": int(g.n_obs), "guides": int(g.n_vars), "guide_umi_limit": args.guide_umi_limit,
            "merge": _truthy(args.merge), "cells_per_guide": per_guide.astype(int).to_dict(),
            "mean_guides_per_cell": float(g.X.sum(1).mean()) if g.n_obs else 0.0, "figures": figs}
    if own:
        finish(d, "Guide binarisation and covariates",
               [f"Guides binarised at UMI > {args.guide_umi_limit}; mean {summ['mean_guides_per_cell']:.2f} guides per cell.", "",
                md_table(["guide", "cells"], per_guide.astype(int).items())], summ)
    args._binarize = summ
    return 0


# ---------------------------------------------------------------------------
# preprocess_bar_multiseq + MultiSeq: deMULTIplex port
# ---------------------------------------------------------------------------

def multiseq_preprocess(r1: Path, r2: Path, cell_ids: "List[str]", cell=(1, 16), umi=(17, 28), tag=(1, 8)):
    """MULTIseq.preProcess: 1-based inclusive substrings; keep reads whose cell barcode is in cellIDs."""
    keep = set(cell_ids)
    rows = []
    for s1, s2 in iter_fastq_pairs([Path(r1), Path(r2)]):
        c = s1[cell[0] - 1:cell[1]]
        if c in keep:
            rows.append((c, s1[umi[0] - 1:umi[1]], s2[tag[0] - 1:tag[1]]))
    return rows


def multiseq_align(read_table, cell_ids: "List[str]", ref: "List[str]"):
    """MULTIseq.align: tag -> reference barcode at Hamming <= 1 (unique), unique UMIs per cell x barcode.

    Returns a DataFrame cells x (ref barcodes, nUMI, nUMI_total)."""
    pd, np = _pd(), _np()
    exact = {b: i for i, b in enumerate(ref)}
    mm: Dict[str, int] = {}
    bad = set()
    for b, i in exact.items():
        for v in hamming1_variants(b):
            if v in exact or (v in mm and mm[v] != i):
                bad.add(v)
            else:
                mm[v] = i
    for v in bad:
        mm.pop(v, None)
    seen: Dict[Tuple[str, int], set] = {}
    total: Dict[str, set] = {}
    for c, u, t in read_table:
        total.setdefault(c, set()).add(u)
        j = exact.get(t, mm.get(t))
        if j is not None:
            seen.setdefault((c, j), set()).add(u)
    ci = {c: i for i, c in enumerate(cell_ids)}
    M = np.zeros((len(cell_ids), len(ref)))
    for (c, j), us in seen.items():
        M[ci[c], j] = len(us)
    df = pd.DataFrame(M, index=cell_ids, columns=ref)
    df["nUMI"] = M.sum(1)
    df["nUMI_total"] = [len(total.get(c, ())) for c in cell_ids]
    return df


def bkde_normal(x, gridsize: int = 401, tau: float = 4.0):
    """KernSmooth::bkde(x, kernel='normal') with its default oversmoothed bandwidth -> (grid, density)."""
    np = _np()
    x = np.asarray(x, float)
    n = len(x)
    sd = x.std(ddof=1) if n > 1 else 0.0
    del0 = (1.0 / (4.0 * math.pi)) ** (1.0 / 10.0)
    h = del0 * (243.0 / (35.0 * n)) ** 0.2 * sd if sd > 0 else 1e-3
    lo, hi = x.min() - tau * h, x.max() + tau * h
    grid = np.linspace(lo, hi, gridsize)
    # linear binning onto the grid, then a Gaussian kernel sum (as bkde)
    delta = (hi - lo) / (gridsize - 1)
    pos = (x - lo) / delta
    li = np.clip(np.floor(pos).astype(int), 0, gridsize - 2)
    rem = pos - li
    w = np.zeros(gridsize)
    np.add.at(w, li, 1 - rem)
    np.add.at(w, li + 1, rem)
    diff = (grid[:, None] - grid[None, :]) / h
    dens = (np.exp(-0.5 * diff ** 2) / math.sqrt(2 * math.pi)) @ w / (n * h)
    return grid, dens


def local_maxima(x) -> "List[int]":
    """deMULTIplex localMaxima (0-based indices)."""
    np = _np()
    x = np.asarray(x, float)
    if len(x) < 2:
        return []
    y = np.diff(np.r_[-float(2 ** 31 - 1), x]) > 0
    ends, cur = [], 0
    runs = [len(list(g)) for _, g in itertools.groupby(y.tolist())]
    for L in runs:
        cur += L
        ends.append(cur)
    idx = [e - 1 for e in ends[0::2]]
    if x[0] == x[1] and idx:
        idx = idx[1:]
    return idx


def _classify_prep(bar):
    """log2, non-finite -> 0, mean-centre; per barcode the x of the lowest / highest KDE maximum (q-independent)."""
    np = _np()
    with np.errstate(divide="ignore"):
        L = np.log2(bar.values.astype(float))
    L[~np.isfinite(L)] = 0
    L = L - L.mean(0)
    ext = []
    for i in range(L.shape[1]):
        col = L[:, i]
        grid, dens = bkde_normal(col)
        xs = np.linspace(np.quantile(col, 0.001), np.quantile(col, 0.999), 100)
        mx = local_maxima(np.interp(xs, grid, dens))
        ext.append(None if len(mx) <= 1 else (xs[min(mx)], xs[max(mx)]))
    return L, ext


def classify_cells(bar, q: float, prep=None):
    """deMULTIplex classifyCells: threshold = quantile(c(x[high max], x[low max]), q) per barcode."""
    np, pd = _np(), _pd()
    L, ext = prep if prep is not None else _classify_prep(bar)
    calls = np.array(["Negative"] * L.shape[0], dtype=object)
    for i in range(L.shape[1]):
        col = L[:, i]
        if ext[i] is None:
            continue
        a, b = sorted(ext[i])
        thresh = a + q * (b - a)
        hit = np.where(col >= thresh)[0]
        for h in hit:
            calls[h] = bar.columns[i] if calls[h] == "Negative" else "Doublet"
    return pd.Series(calls, index=bar.index)


def find_thresh(call_list: "Dict[float, Any]"):
    pd = _pd()
    rows = []
    for q, calls in call_list.items():
        vc = calls.value_counts()
        n = len(calls)
        pdbl = int(vc.get("Doublet", 0)); pneg = int(vc.get("Negative", 0))
        rows.append({"q": q, "pDoublet": pdbl / n, "pNegative": pneg / n, "pSinglet": (n - pdbl - pneg) / n})
    res = pd.DataFrame(rows)
    ext = [res["q"].iloc[i] for i in local_maxima(res["pSinglet"].values)]
    return res, ext


def find_q(res, extrema) -> float:
    use = res[res["q"].isin(extrema)]
    if not len(use):
        return float(res.loc[res["pSinglet"].idxmax(), "q"])
    return float(use.loc[use["pSinglet"].idxmax(), "q"])


def multiseq_class_columns(bar_table, upstream_compat: bool = False) -> "List[str]":
    """Columns classified: upstream drops only nUMI (ncol(dt) is NULL in R), so nUMI_total is classified too."""
    drop = {"nUMI"} if upstream_compat else {"nUMI", "nUMI_total"}
    return [c for c in bar_table.columns if c not in drop]


def multiseq_drop_classes(upstream_compat: bool = False) -> set:
    """Calls removed from the MuData: upstream tests for 'Double' (never produced), so doublets survive."""
    return {"Double", "Negative"} if upstream_compat else {"Doublet", "Negative"}


def multiseq_classify(bar_table, upstream_compat: bool = False):
    """Two quantile-sweep rounds (q = 0.01..0.99 by 0.02), negatives removed after round 1."""
    np = _np()
    bt = bar_table[multiseq_class_columns(bar_table, upstream_compat)]
    qs = [round(x, 2) for x in np.arange(0.01, 0.99 + 1e-9, 0.02)]
    prep = _classify_prep(bt)
    sweep = {q: classify_cells(bt, q, prep) for q in qs}
    res1, ext1 = find_thresh(sweep)
    q1 = find_q(res1, ext1)
    r1 = classify_cells(bt, q1, prep)
    neg = list(r1.index[r1 == "Negative"])
    bt2 = bt.drop(index=neg)
    q2 = None
    if len(bt2) > 2:
        prep2 = _classify_prep(bt2)
        sweep2 = {q: classify_cells(bt2, q, prep2) for q in qs}
        res2, ext2 = find_thresh(sweep2)
        q2 = find_q(res2, ext2)
        r2 = classify_cells(bt2, q2, prep2)
        neg += list(r2.index[r2 == "Negative"])
        calls = r2[r2 != "Negative"]
    else:
        calls = r1[r1 != "Negative"]
        res2 = None
    final = calls.copy()
    for c in neg:
        final[c] = "Negative"
    return final, {"q_round1": q1, "q_round2": q2, "sweep_round1": res1}


def cmd_multiseq(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd, np, ad = _pd(), _np(), _ad()
    own = d is None
    d = d or run_dir(args.label)
    mods = load_bundle(Path(args.muon_data))
    cells = list(mods["scRNA"].obs.index)
    pd.DataFrame(cells).to_csv(d / "cell_barcode_capturing.csv", index=False, header=False)
    ref = [l.split(",")[0].strip() for l in Path(args.barcodes).read_text().splitlines() if l.strip()]
    rt = multiseq_preprocess(Path(args.r1), Path(args.r2), cells, tuple(args.bar), tuple(args.umi), tuple(args.tag))
    bar = multiseq_align(rt, cells, ref)
    bar.to_csv(d / "bar_table.csv")
    print(f"CSV: {d / 'bar_table.csv'}")
    final, info = multiseq_classify(bar, args.upstream_compat)
    final = final.reindex(bar.index)
    pd.DataFrame({"x": final.values}).to_csv(d / "final_class.csv")
    pd.DataFrame({"x": final.index}).to_csv(d / "final_class_cell_barcode.csv")
    print(f"CSV: {d / 'final_class.csv'}")
    write_tsv(info["sweep_round1"], d / "quantile_sweep_round1.tsv")
    ms = ad.AnnData(X=bar[ref].values.astype(np.float32), obs=pd.DataFrame(index=bar.index), var=pd.DataFrame(index=ref))
    ms.obs["multiseq_class"] = final.values
    drop = multiseq_drop_classes(args.upstream_compat)
    keep = ~ms.obs["multiseq_class"].isin(drop).values
    out = {"scRNA": mods["scRNA"][keep].copy(), "guides": mods["guides"][keep].copy(), "multiseq": ms[keep].copy()}
    for k in ("scRNA", "guides"):
        out[k].obs["multiseq_class"] = ms.obs["multiseq_class"].values[keep]
    bundle = save_bundle(out, d / "processed_mudata_guide_and_transcripts_multiseq_filtered")
    vc = final.value_counts().to_dict()
    figs: list = []
    plt = None if args.no_plots else _plt()
    if plt is not None:
        r = info["sweep_round1"]
        fig, ax = plt.subplots(figsize=(5, 3.5))
        for col, c in (("pSinglet", BLUE), ("pDoublet", ORANGE), ("pNegative", INK2)):
            ax.plot(r["q"], r[col], color=c, label=col)
        ax.axvline(info["q_round1"], ls="--", color=INK)
        ax.legend(frameon=False, fontsize=8)
        _style(ax, "MULTI-seq quantile sweep (round 1)", "q", "proportion")
        _save(fig, d / "class_threshold_check.png", figs)
    summ = {"reads_matched_to_cells": len(rt), "cells": len(cells), "calls": vc, "kept_cells": int(keep.sum()),
            "q_round1": info["q_round1"], "q_round2": info["q_round2"], "bundle": str(bundle), "figures": figs,
            "upstream_compat": bool(args.upstream_compat)}
    if own:
        finish(d, "MULTI-seq demultiplexing", ["deMULTIplex quantile-sweep classification (two rounds).", "",
                                                md_table(["class", "cells"], sorted(vc.items()))], summ)
    args._multiseq = summ
    return 0


# ---------------------------------------------------------------------------
# PerturbLoaderGeneration
# ---------------------------------------------------------------------------

def guide_type(e: str) -> str:
    if "_TSS" in e:
        return "POSITIVE_CONTROL"
    if "random" in e or "scrambled" in e:
        return "NEGATIVE_CONTROL"
    if "chr" in e:
        return "PUTATIVE_ENHANCER"
    return "None"


def _nochr(c) -> str:
    return str(c).replace("chr", "")


def perturb_loader(mods, genes, distance: int = 1_000_000, in_trans: bool = False, add_genes: "Optional[List[str]]" = None,
                   upstream_compat: bool = False, seed: int = 0):
    """-> (element x guide membership, element x tested gene binary, element info)."""
    pd, np = _pd(), _np()
    rng = np.random.RandomState(seed)
    g, r = mods["guides"], mods["scRNA"]
    gids = [str(x) for x in g.var.index]
    elements = list(dict.fromkeys(x.split("|")[0] for x in gids))
    memb = pd.DataFrame([[1 if x.split("|")[0] == e else 0 for x in gids] for e in elements], index=elements, columns=gids)
    coords = {k.split("|")[0]: k.split("_")[-1] for k in gids}
    expr = [str(x) for x in r.var.index]
    ge = genes[genes["gene_name"].isin(set(expr))].drop_duplicates("gene_name")
    add = [x for x in (add_genes or []) if x]
    rows, info = [], []
    for e in elements:
        q = e.replace("_TSS", "")
        typ = guide_type(e)
        anchor = ""
        if in_trans:
            s = set(ge["gene_name"])
            anchor = "trans"
        else:
            hit = ge[ge["gene_name"] == q]
            if len(hit):
                chrom, pos = hit.iloc[0]["chr"], int(hit.iloc[0]["gene_real_start"])
                same = ge[ge["chr"].map(_nochr) == _nochr(chrom)]
                s = set(same.loc[(same["gene_real_start"] - pos).abs() < distance, "gene_name"])
                anchor = f"gene_tss:{chrom}:{pos}"
            elif not upstream_compat and typ != "NEGATIVE_CONTROL" and ":" in coords.get(e, ""):
                c = coords[e].split(":")
                chrom, pos = c[0], int(c[1])
                same = ge[ge["chr"].map(_nochr) == _nochr(chrom)]
                s = set(same.loc[(same["gene_real_start"] - pos).abs() < distance, "gene_name"])
                anchor = f"guide:{chrom}:{pos}"
            else:
                k = min(NUMBER_FOR_RANDOM_CONTROLS, len(ge))
                s = set(ge.sample(k, random_state=rng)["gene_name"]) if k else set()
                anchor = f"random_{k}"
        s = s.union(add)
        rows.append([1 if x in s else 0 for x in expr])
        info.append({"element": e, "GUIDE_TYPE": typ, "anchor": anchor, "n_guides": int(memb.loc[e].sum()), "n_genes_tested": int(sum(rows[-1]))})
    tested = pd.DataFrame(rows, index=elements, columns=expr)
    return memb, tested, pd.DataFrame(info)


def cmd_perturb_loader(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    own = d is None
    d = d or run_dir(args.label)
    mods = load_bundle(Path(args.muon_data))
    genes = read_gtf_genes(Path(args.gtf))
    memb, tested, info = perturb_loader(mods, genes, args.distance_from_guide, _truthy(args.in_trans),
                                        (args.add_gene_names or "").split(","), args.upstream_compat, args.seed)
    pl = d / "perturbdata"
    pl.mkdir(parents=True, exist_ok=True)
    write_tsv(memb, pl / "elements_x_guides.tsv", index=True)
    write_tsv(tested, pl / "elements_x_tested_genes.tsv", index=True)
    write_tsv(info, pl / "elements.tsv")
    (pl / "perturbdata.json").write_text(json.dumps({"muon_data": str(Path(args.muon_data).resolve()), "in_trans": _truthy(args.in_trans),
                                                     "distance_from_guide": args.distance_from_guide}))
    print(f"JSON: {pl / 'perturbdata.json'}")
    summ = {"perturbdata": str(pl), "elements": info.to_dict("records"), "pairs": int(sum(
        info["n_guides"][i] * info["n_genes_tested"][i] for i in range(len(info))))}
    if own:
        finish(d, "Perturbation loader", ["Elements, their guides and the cis genes tested for each.", "",
                                           md_table(list(info.columns), info.values)], summ)
    args._loader = summ
    return 0


# ---------------------------------------------------------------------------
# runSceptre: differential perturbation
# ---------------------------------------------------------------------------

def design_matrix(cov, drop_single_batch: bool = True):
    """Intercept + standardised numeric covariates + batch dummies (batch dropped when single-valued, as runSceptre)."""
    pd, np = _pd(), _np()
    cols = [np.ones(len(cov))]
    names = ["intercept"]
    for c in cov.columns:
        if c == "bath_number":
            lv = sorted(cov[c].astype(str).unique())
            if len(lv) == 1 and drop_single_batch:
                continue
            for l in lv[1:]:
                cols.append((cov[c].astype(str) == l).astype(float).values)
                names.append(f"batch_{l}")
            continue
        v = pd.to_numeric(cov[c], errors="coerce").fillna(0).values.astype(float)
        v[~np.isfinite(v)] = 0
        sd = v.std()
        cols.append((v - v.mean()) / sd if sd > 0 else v * 0)
        names.append(c)
    Z = np.column_stack(cols)
    return Z, names


def _irls(Z, y, family: str, max_iter: int = 50, ridge: float = 1e-6, offset=None):
    np = _np()
    n, p = Z.shape
    off = np.zeros(n) if offset is None else offset
    beta = np.zeros(p)
    if family == "poisson":
        beta[0] = math.log(max(y.mean(), 1e-8))
    for _ in range(max_iter):
        eta = np.clip(Z @ beta + off, -30, 30)
        if family == "poisson":
            mu = np.exp(eta); w = mu; z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        else:
            mu = 1 / (1 + np.exp(-eta)); w = np.maximum(mu * (1 - mu), 1e-9); z = eta - off + (y - mu) / w
        A = Z.T @ (Z * w[:, None]) + ridge * np.eye(p)
        new = np.linalg.solve(A, Z.T @ (w * z))
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    eta = np.clip(Z @ beta + off, -30, 30)
    return beta, (np.exp(eta) if family == "poisson" else 1 / (1 + np.exp(-eta)))


def fit_nb_theta(y, mu) -> float:
    from scipy import optimize, special  # type: ignore
    np = _np()

    def nll(lt):
        th = math.exp(lt)
        return -np.sum(special.gammaln(y + th) - special.gammaln(th) - special.gammaln(y + 1)
                       + th * np.log(th / (th + mu)) + y * np.log(np.maximum(mu, 1e-12) / (th + mu)))
    r = optimize.minimize_scalar(nll, bounds=(math.log(1e-2), math.log(1e4)), method="bounded")
    return float(math.exp(r.x))


def skew_normal_pvalue(z_obs: float, z_null, side: str) -> "Tuple[float, float]":
    """SCEPTRE calibrates the resampled null with a fitted skew distribution (skew-t upstream; skew-normal here).

    The fitted tail is used only beyond the resolution of the resamples (empirical p < 10 / (B + 1)); inside it the
    empirical resampling p is returned, which keeps the bulk of the null exactly calibrated."""
    np = _np()
    from scipy import stats  # type: ignore
    B = len(z_null)
    zl = float((np.sum(z_null <= z_obs) + 1) / (B + 1))
    zr = float((np.sum(z_null >= z_obs) + 1) / (B + 1))
    emp = zl if side == "left" else zr if side == "right" else min(1.0, 2 * min(zl, zr))
    try:
        a, loc, sc = stats.skewnorm.fit(z_null)
        left, right = float(stats.skewnorm.cdf(z_obs, a, loc, sc)), float(stats.skewnorm.sf(z_obs, a, loc, sc))
        p = left if side == "left" else right if side == "right" else min(1.0, 2 * min(left, right))
        if not np.isfinite(p):
            raise ValueError
        if emp >= 10.0 / (B + 1):
            p = emp
    except Exception:
        p = emp
    return max(p, 1e-300), emp


def sceptre_tests(Y, Xg, Z, side: str = "both", B: int = 500, seed: int = 0, gene_names=None, guide_names=None):
    """Conditional-resampling NB score test for every (gene column of Y) x (guide column of Xg).

    Y: cells x genes counts (dense); Xg: cells x guides binary.  Returns list of dicts."""
    np = _np()
    rng = np.random.default_rng(seed)
    n = Y.shape[0]
    nulls = []
    for j in range(Xg.shape[1]):
        x = Xg[:, j].astype(float)
        _, pi = _irls(Z, x, "logit")
        nulls.append((x, (rng.random((B, n)) < pi[None, :]).astype(float)))
    out = []
    for gi in range(Y.shape[1]):
        y = Y[:, gi].astype(float)
        if y.sum() == 0:
            for j in range(Xg.shape[1]):
                out.append({"gene_id": gene_names[gi], "gRNA_id": guide_names[j], "p_value": float("nan"), "z_value": float("nan"),
                            "log_fold_change": float("nan"), "p_empirical": float("nan")})
            continue
        _, mu = _irls(Z, y, "poisson")
        th = fit_nb_theta(y, mu)
        denom = 1 + mu / th
        r, w = (y - mu) / denom, mu / denom
        ZW = Z * w[:, None]
        Inv = np.linalg.pinv(Z.T @ ZW)
        for j, (x, Xn) in enumerate(nulls):
            def zstat(X):
                U = X @ r
                A = X @ ZW
                V = (X * X) @ w - np.einsum("ij,jk,ik->i", A, Inv, A)
                return U / np.sqrt(np.maximum(V, 1e-12))
            zo = float(zstat(x[None, :])[0])
            zn = zstat(Xn)
            zn = zn[np.isfinite(zn)]
            p, pe = skew_normal_pvalue(zo, zn, side)
            t = x > 0
            lfc = float(math.log((y[t].sum() + 0.5) / (mu[t].sum() + 0.5))) if t.any() else float("nan")
            out.append({"gene_id": gene_names[gi], "gRNA_id": guide_names[j], "p_value": p, "z_value": zo,
                        "log_fold_change": lfc, "p_empirical": pe})
    return out


def mannwhitney_tests(Y, Xg, side: str = "both", gene_names=None, guide_names=None):
    """Task 3 alternative: Mann-Whitney U on log1p(CP10k), treated vs untreated cells."""
    np = _np()
    from scipy import stats  # type: ignore
    tot = Y.sum(1)
    N = np.log1p(Y / np.maximum(tot[:, None], 1) * 1e4)
    alt = {"left": "less", "right": "greater", "both": "two-sided"}[side]
    out = []
    for j in range(Xg.shape[1]):
        t = Xg[:, j] > 0
        n1, n2 = int(t.sum()), int((~t).sum())
        for gi in range(Y.shape[1]):
            a, b = N[t, gi], N[~t, gi]
            if n1 == 0 or n2 == 0:
                p, z = float("nan"), float("nan")
            else:
                U, p = stats.mannwhitneyu(a, b, alternative=alt)
                z = (U - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
            ma, mb = Y[t, gi].mean() if n1 else 0, Y[~t, gi].mean() if n2 else 0
            out.append({"gene_id": gene_names[gi], "gRNA_id": guide_names[j], "p_value": float(p), "z_value": float(z),
                        "log_fold_change": float(math.log((ma + 1e-3) / (mb + 1e-3))), "p_empirical": float("nan")})
    return out


SCEPTRE_R = '''setwd("{path}")
library(sceptre)
library(tibble)
library(dplyr)
exp = as.matrix(read.table('gene_exp_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
g_I = as.matrix(read.table('guides_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
cov = read.table('covariates.txt', sep=',', header=TRUE, row.names=1)
if ('bath_number' %in% names(cov)) {{ cov$bath_number = as.factor(cov$bath_number) }}
pairs_test = as_tibble(read.table('pairs.txt', sep=',', header=TRUE, check.names=FALSE))
pairs_test$pair_type = as.factor(pairs_test$pair_type)
result <- run_sceptre_high_moi(gene_matrix = exp, combined_perturbation_matrix = g_I, covariate_matrix = cov,
                               gene_gRNA_group_pairs = pairs_test, side = '{side}')
write.table(result, "results.txt")
'''


def cmd_de(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd, np = _pd(), _np()
    own = d is None
    d = d or run_dir(args.label)
    pl = Path(args.perturbdata)
    meta = json.loads((pl / "perturbdata.json").read_text())
    mods = load_bundle(Path(args.muon_data or meta["muon_data"]))
    memb = pd.read_csv(pl / "elements_x_guides.tsv", sep="\t", index_col=0)
    tested = pd.read_csv(pl / "elements_x_tested_genes.tsv", sep="\t", index_col=0)
    g, r = mods["guides"], mods["scRNA"]
    cov = r.obs[[c for c in COVARIATES if c in r.obs.columns]].copy()
    Z, znames = design_matrix(cov)
    G = dense(g.X)
    gcol = {str(x): i for i, x in enumerate(g.var.index)}
    rcol = {str(x): i for i, x in enumerate(r.var.index)}
    out_dir = d / "sceptre_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    allres, skipped, n_el = [], [], 0
    for e in memb.index:
        guides = [x for x in memb.columns if memb.loc[e, x] == 1 and x in gcol]
        genes = [x for x in tested.columns if tested.loc[e, x] == 1 and x in rcol]
        cells_per = {x: int(G[:, gcol[x]].sum()) for x in guides}
        use = [x for x in guides if cells_per[x] > args.min_cells_per_guide]
        skipped += [{"element": e, "gRNA_id": x, "cells": cells_per[x]} for x in guides if x not in use]
        if not use or not genes:
            continue
        n_el += 1
        ed = out_dir / safe_label(e)
        ed.mkdir(parents=True, exist_ok=True)
        Y = dense(r.X[:, [rcol[x] for x in genes]])
        Xg = G[:, [gcol[x] for x in use]]
        if args.engine == "sceptre-r":
            pd.DataFrame(Y.T, index=genes, columns=r.obs.index).to_csv(ed / "gene_exp_one_gene_guide.txt")
            pd.DataFrame(Xg.T, index=use, columns=g.obs.index).to_csv(ed / "guides_one_gene_guide.txt")
            cm = cov.copy()
            if cm["bath_number"].nunique() == 1:
                del cm["bath_number"]
            cm.to_csv(ed / "covariates.txt")
            pd.DataFrame([[gn, x, "all_elements_test"] for x in use for gn in genes],
                         columns=["gene_id", "gRNA_group", "pair_type"]).to_csv(ed / "pairs.txt", index=False)
            (ed / "r_script.r").write_text(SCEPTRE_R.format(path=str(ed), side=args.direction))
            run_or_print(["Rscript", str(ed / "r_script.r")])
            if (ed / "results.txt").exists():
                rr = pd.read_csv(ed / "results.txt", sep=" ")
                rr["element"] = e
                allres.append(rr)
            continue
        if args.engine == "mannwhitney":
            res = mannwhitney_tests(Y, Xg, args.direction, genes, use)
        else:
            res = sceptre_tests(Y, Xg, Z, args.direction, args.B, args.seed, genes, use)
        rr = pd.DataFrame(res)
        rr["pair_type"] = "all_elements_test"
        rr = rr[RESULT_COLS + ["p_empirical"]]
        rr.to_csv(ed / "results.txt", sep=" ", index=False)
        rr["element"] = e
        allres.append(rr)
    res = pd.concat(allres, ignore_index=True) if allres else pd.DataFrame(columns=RESULT_COLS + ["element"])
    write_tsv(res, d / "de_results.tsv")
    if skipped:
        write_tsv(pd.DataFrame(skipped), d / "guides_below_min_cells.tsv")
    summ = {"engine": args.engine, "side": args.direction, "B": args.B, "elements_tested": n_el, "tests": int(len(res)),
            "guides_skipped_min_cells": len(skipped), "design": znames, "sceptre_out": str(out_dir), "results": str(d / "de_results.tsv")}
    if own:
        finish(d, "Differential perturbation", [f"{len(res)} guide x gene tests over {n_el} elements ({args.engine}, side {args.direction}).", "",
                                                 md_table(RESULT_COLS, res.sort_values("p_value")[RESULT_COLS].values, 25) if len(res) else "no tests"], summ)
    args._de = summ
    return 0


# ---------------------------------------------------------------------------
# create_anndata_from_sceptre
# ---------------------------------------------------------------------------

def bh(p):
    np = _np()
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    if not len(q):
        return out
    o = np.argsort(q)
    ranked = q[o] * len(q) / (np.arange(len(q)) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(len(q))
    res[o] = np.minimum(adj, 1)
    out[ok] = res
    return out


def aggregate_results(res, mods, alpha_guide: float = 0.01, alpha_element: float = 0.05):
    pd, np, ad = _pd(), _np(), _ad()
    from scipy.stats import combine_pvalues  # type: ignore
    res = res.copy()
    res["adj_pvalue"] = bh(res["p_value"].values)
    gv = mods["guides"].var
    for c in ("target_elements", "guide_chr", "guide_start", "guide_end"):
        if c in gv.columns:
            res[c] = res["gRNA_id"].map(gv[c].astype(str).to_dict())
    if "target_elements" not in res.columns:
        res["target_elements"] = res["gRNA_id"].astype(str).str.split("|").str[0]
    res["significant"] = res["adj_pvalue"] < alpha_guide
    guides = list(dict.fromkeys(res["gRNA_id"]))
    genes = list(dict.fromkeys(res["gene_id"]))
    layers = {}
    for col in ("p_value", "adj_pvalue", "z_value", "log_fold_change"):
        layers[col] = res.pivot_table(index="gRNA_id", columns="gene_id", values=col, aggfunc="first").reindex(index=guides, columns=genes)
    rg = ad.AnnData(X=layers["p_value"].values, obs=pd.DataFrame(index=guides), var=pd.DataFrame(index=genes))
    for k in ("adj_pvalue", "z_value", "log_fold_change"):
        rg.layers[k] = layers[k].values
    rg.layers["significant"] = (np.nan_to_num(layers["adj_pvalue"].values, nan=1.0) < alpha_guide).astype(np.int8)
    obs = gv.reindex(guides)
    rg.obs = obs.copy() if len(obs.columns) else rg.obs
    rows = []
    for (el, gene), sub in res.groupby(["target_elements", "gene_id"], sort=False):
        ps = [p for p in sub["p_value"] if np.isfinite(p)]
        pc = float(combine_pvalues(ps).pvalue) if ps else float("nan")
        first = sub.iloc[0]
        rows.append({"element": el, "gene_id": gene, "p_combined": pc, "n_guides": len(ps),
                     "element_chr": first.get("guide_chr", ""), "element_start": first.get("guide_start", ""),
                     "element_end": first.get("guide_end", ""), "mean_log_fold_change": float(sub["log_fold_change"].mean())})
    el = pd.DataFrame(rows)
    el["sig_not_adj"] = el["p_combined"] < alpha_element
    em = el.pivot_table(index="element", columns="gene_id", values="p_combined", aggfunc="first")
    re_ = ad.AnnData(X=em.values, obs=pd.DataFrame(index=[str(x) for x in em.index]), var=pd.DataFrame(index=[str(x) for x in em.columns]))
    info = el.drop_duplicates("element").set_index("element")[["element_chr", "element_start", "element_end"]]
    re_.obs["element"] = re_.obs.index
    for c in info.columns:
        re_.obs[c] = info.reindex(re_.obs.index)[c].astype(str).values
    re_.layers["sig_not_adj"] = (np.nan_to_num(em.values, nan=1.0) < alpha_element).astype(np.int8)
    return res, el, rg, re_


def cmd_results(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd = _pd()
    own = d is None
    d = d or run_dir(args.label)
    src = Path(args.sceptre_dir)
    files = sorted(src.rglob("results.txt"))
    res = pd.concat([pd.read_csv(f, sep=" ") for f in files], ignore_index=True) if files else pd.DataFrame(columns=RESULT_COLS)
    mods = load_bundle(Path(args.muon_data))
    res, el, rg, re_ = aggregate_results(res, mods, args.alpha_guide, args.alpha_element)
    write_tsv(res, d / "result_guides.tsv")
    write_tsv(el, d / "result_elements.tsv")
    bundle = save_bundle({"guides": mods["guides"], "scRNA": mods["scRNA"], "result_guides": rg, "result_elements": re_}, d / "mudata_results")
    sig = res[res["significant"]]
    summ = {"tests": int(len(res)), "significant_guide_gene": int(len(sig)), "elements_sig_not_adj": int(el["sig_not_adj"].sum()),
            "bundle": str(bundle), "result_guides": str(d / "result_guides.tsv"), "result_elements": str(d / "result_elements.tsv"),
            "top": res.sort_values("p_value").head(10)[["gRNA_id", "gene_id", "p_value", "adj_pvalue", "log_fold_change"]].to_dict("records")}
    if own:
        finish(d, "Perturbation results", [f"{len(res)} tests from {len(files)} results.txt files; {len(sig)} guide-gene pairs at BH adj < {args.alpha_guide}.", "",
                                            md_table(["gRNA_id", "gene_id", "p_value", "adj_pvalue", "log_fold_change"],
                                                     res.sort_values("p_value")[["gRNA_id", "gene_id", "p_value", "adj_pvalue", "log_fold_change"]].values, 25)], summ)
    args._results = summ
    return 0


# ---------------------------------------------------------------------------
# Task 4: genome-browser tracks
# ---------------------------------------------------------------------------

def make_tracks(res, rna_var, out: Path, alpha: float = 0.05, all_pairs: bool = False, plots: bool = True, figs: Optional[list] = None):
    pd, np = _pd(), _np()
    out.mkdir(parents=True, exist_ok=True)
    tss = {}
    if "transcript_chr" in rna_var.columns:
        for gname, row in rna_var.iterrows():
            if str(row["transcript_chr"]) != "NOT_FOUND":
                tss[str(gname)] = (str(row["transcript_chr"]), int(row["transcript_start"]))
    r = res.copy()
    r = r[np.isfinite(r["p_value"].astype(float))]
    if not all_pairs:
        r = r[r["adj_pvalue"] < alpha]
    bed, links = [], []
    for _, x in r.iterrows():
        try:
            gc, gs, ge = str(x["guide_chr"]), int(float(x["guide_start"])), int(float(x["guide_end"]))
        except Exception:
            continue
        nl = -math.log10(max(float(x["p_value"]), 1e-300))
        bed.append((gc, gs, ge, f"{x['gRNA_id']}|{x['gene_id']}", min(1000, int(round(nl * 100))), "."))
        if str(x["gene_id"]) in tss:
            tc, tp = tss[str(x["gene_id"])]
            tc = tc if str(tc).startswith("chr") or not gc.startswith("chr") else "chr" + tc
            links.append((gc, gs, ge, tc, tp, tp + 1, round(-math.log10(max(float(x["adj_pvalue"]), 1e-300)), 4)))
    bdf = pd.DataFrame(bed, columns=["chr", "start", "end", "name", "score", "strand"]).sort_values(["chr", "start"])
    ldf = pd.DataFrame(links, columns=["chr1", "start1", "end1", "chr2", "start2", "end2", "score"])
    bdf.to_csv(out / "guide_tests.bed", sep="\t", header=False, index=False)
    ldf.to_csv(out / "guide_gene.links", sep="\t", header=False, index=False)
    print(f"Wrote: {out / 'guide_tests.bed'}")
    print(f"Wrote: {out / 'guide_gene.links'}")
    ini = ("[x-axis]\n\n[guide gene links]\nfile = guide_gene.links\ntitle = guide -> gene (-log10 adj p)\nlinks_type = arcs\n"
           "color = viridis\nline_width = 2\nheight = 4\n\n[spacer]\n\n[guide tests]\nfile = guide_tests.bed\ntitle = tested guides\n"
           "height = 3\ndisplay = collapsed\n")
    (out / "tracks.ini").write_text(ini)
    print(f"Wrote: {out / 'tracks.ini'}")
    plt = _plt() if plots else None
    if plt is not None and len(ldf):
        chrom = ldf["chr1"].value_counts().index[0]
        sub = ldf[ldf["chr1"] == chrom]
        fig, ax = plt.subplots(figsize=(10, 3))
        from matplotlib.patches import Arc  # type: ignore
        for _, x in sub.iterrows():
            a, b = (x["start1"] + x["end1"]) / 2, x["start2"]
            c, w = (a + b) / 2, abs(b - a)
            ax.add_patch(Arc((c, 0), max(w, 1), max(w, 1) * 0.5, theta1=0, theta2=180, color=BLUE, lw=0.5 + min(x["score"], 10) / 3))
        lo = min(sub["start1"].min(), sub["start2"].min()); hi = max(sub["end1"].max(), sub["end2"].max())
        ax.set_xlim(lo - 0.05 * (hi - lo + 1), hi + 0.05 * (hi - lo + 1)); ax.set_ylim(0, max((hi - lo) * 0.3, 1))
        ax.set_yticks([])
        _style(ax, f"{chrom}: significant guide -> gene links", "position (bp)")
        _save(fig, out / f"links_{chrom}.png", figs)
    return bdf, ldf


def cmd_tracks(args: argparse.Namespace, d: Optional[Path] = None) -> int:
    pd = _pd()
    own = d is None
    d = d or run_dir(args.label)
    res = pd.read_csv(args.results, sep="\t")
    mods = load_bundle(Path(args.muon_data))
    figs: list = []
    bdf, ldf = make_tracks(res, mods["scRNA"].var, d / "tracks_dir", args.alpha, args.all_pairs, not args.no_plots, figs)
    region = args.region
    if region:
        run_or_print(["pyGenomeTracks", "--tracks", str(d / "tracks_dir" / "tracks.ini"), "--region", region,
                      "--outFileName", str(d / "tracks_dir" / f"{safe_label(region)}.png")])
    summ = {"bed_rows": int(len(bdf)), "links": int(len(ldf)), "tracks_dir": str(d / "tracks_dir"), "figures": figs}
    if own:
        finish(d, "Genome-browser tracks", [f"{len(bdf)} guide intervals and {len(ldf)} guide -> TSS links "
                                             f"({'all tests' if args.all_pairs else f'adj p < {args.alpha}'}); pyGenomeTracks tracks.ini written."], summ)
    args._tracks = summ
    return 0


# ---------------------------------------------------------------------------
# Task guide-assign: depth-aware guide calls
# ---------------------------------------------------------------------------

def poisson_mixture_calls(C, max_iter: int = 200, init_limit: float = 5.0):
    """Per guide: counts ~ (1-pi) Poisson(s_i l0) + pi Poisson(l1); s_i = the cell's ambient guide depth
    (UMIs summed over its guides at <= init_limit, +1, scaled to mean 1), so background expectations follow depth.

    Returns (binary calls, per-guide parameter table)."""
    np, pd = _np(), _pd()
    from scipy.special import gammaln  # type: ignore
    C = np.asarray(C, float)
    n, m = C.shape
    calls = np.zeros_like(C, dtype=np.int8)
    rows = []
    amb = np.where(C <= init_limit, C, 0).sum(1) + 1.0
    s_all = amb / amb.mean()
    for j in range(m):
        c = C[:, j]
        s = s_all
        pos = c > init_limit
        pi = min(max(pos.mean(), 1e-3), 0.999)
        l0 = max((c[~pos] / s[~pos]).mean() if (~pos).any() else 0.1, 1e-3)
        l1 = max(c[pos].mean() if pos.any() else 10 * l0 + 1, l0 * 2)
        for _ in range(max_iter):
            ll0 = c * np.log(s * l0) - s * l0 - gammaln(c + 1) + math.log(1 - pi)
            ll1 = c * math.log(l1) - l1 - gammaln(c + 1) + math.log(pi)
            mx = np.maximum(ll0, ll1)
            r = np.exp(ll1 - mx) / (np.exp(ll0 - mx) + np.exp(ll1 - mx))
            new_pi = min(max(r.mean(), 1e-4), 1 - 1e-4)
            new_l0 = max(((1 - r) * c).sum() / max(((1 - r) * s).sum(), 1e-12), 1e-4)
            new_l1 = max((r * c).sum() / max(r.sum(), 1e-12), new_l0 * 1.01)
            done = abs(new_pi - pi) < 1e-7 and abs(new_l1 - l1) < 1e-6 and abs(new_l0 - l0) < 1e-6
            pi, l0, l1 = new_pi, new_l0, new_l1
            if done:
                break
        calls[:, j] = ((r > 0.5) & (c >= 1)).astype(np.int8)
        rows.append({"lambda_background": l0, "lambda_signal": l1, "pi": pi, "cells_mixture": int(calls[:, j].sum())})
    return calls, pd.DataFrame(rows)


def cmd_assign_guides(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    d = run_dir(args.label)
    mods = load_bundle(Path(args.muon_data))
    g = mods["guides"]
    C = dense(g.layers["counts"]) if "counts" in g.layers else dense(g.X)
    calls, tab = poisson_mixture_calls(C, init_limit=args.guide_umi_limit)
    thr = (C > args.guide_umi_limit).astype(np.int8)
    tab.insert(0, "guide", [str(x) for x in g.var.index])
    tab["cells_threshold"] = thr.sum(0)
    tab["agreement"] = (calls == thr).mean(0)
    g = g.copy()
    g.layers["binarized"] = calls
    g.layers["threshold_binary"] = thr
    mods = dict(mods)
    mods["guides"] = g
    bundle = save_bundle(mods, d / "mu_with_binary")
    write_tsv(tab, d / "guide_assignment.tsv")
    summ = {"bundle": str(bundle), "guides": int(len(tab)), "mean_agreement": float(tab["agreement"].mean()),
            "cells_mixture_total": int(calls.sum()), "cells_threshold_total": int(thr.sum())}
    finish(d, "Guide assignment", ["Depth-aware two-component Poisson mixture per guide (layer 'binarized') next to the "
                                   f"fixed UMI > {args.guide_umi_limit} rule (layer 'threshold_binary').", "",
                                   md_table(list(tab.columns), tab.values)], summ)
    return 0


# ---------------------------------------------------------------------------
# Task cellranger, rename_and_symlink, config, inspect
# ---------------------------------------------------------------------------

def cmd_cellranger_inputs(args: argparse.Namespace) -> int:
    pd = _pd()
    d = run_dir(args.label)
    t = process_guide_table(read_table_any(Path(args.guides)))
    fr = pd.DataFrame({"id": [x.replace(",", "_") for x in t["pipeline_id"]],
                       "name": t["sgRNA_ID"] if "sgRNA_ID" in t else t["pipeline_id"],
                       "read": args.read, "pattern": args.pattern, "sequence": t["sgRNA_sequences"],
                       "feature_type": "CRISPR Guide Capture",
                       "target_gene_id": t["Target_name"].str.split("|").str[0],
                       "target_gene_name": t["Target_name"].str.split("|").str[0]})
    fr.to_csv(d / "feature_ref.csv", index=False)
    print(f"CSV: {d / 'feature_ref.csv'}")
    lib = pd.DataFrame([{"fastqs": args.rna_fastq_dir, "sample": args.rna_sample, "library_type": "Gene Expression"},
                        {"fastqs": args.guide_fastq_dir, "sample": args.guide_sample, "library_type": "CRISPR Guide Capture"}])
    lib.to_csv(d / "library.csv", index=False)
    print(f"CSV: {d / 'library.csv'}")
    cmd = ["cellranger", "count", f"--id={safe_label(args.label)}", f"--libraries={d / 'library.csv'}",
           f"--transcriptome={args.transcriptome}", f"--feature-ref={d / 'feature_ref.csv'}"]
    rc = run_or_print(cmd, cwd=d, dry=not args.execute)
    finish(d, "Cell Ranger inputs", ["```bash", " ".join(cmd), "```"], {"feature_ref_rows": int(len(fr)), "command": " ".join(cmd)})
    return rc


def cmd_rename_fastq(args: argparse.Namespace) -> int:
    pd = _pd()
    df = pd.read_csv(args.csv)
    d = Path(args.out_dir) if args.out_dir else run_dir(args.label)
    sd = d / "symlink_dir"
    sd.mkdir(parents=True, exist_ok=True)
    names = []
    for _, r in df.iterrows():
        nm = f"{r['Sample']}_S{int(r['Sample_Number']):03d}_L{int(r['Lane_Number']):03d}_{r['Read_Type']}_{int(r['File_Fragment']):03d}.fastq.gz"
        names.append(nm)
        link = sd / nm
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(os.path.abspath(str(r["file_dir"])), link)
    df["New_Name"] = names
    df.to_csv(d / "output.csv", index=False)
    print(f"CSV: {d / 'output.csv'}")
    print(f"Wrote: {sd}")
    return 0


def _parse_nf_value(v: str):
    v = v.strip().rstrip(";").strip()
    if v.startswith("[") and v.endswith("]"):
        inner = v[1:-1].strip()
        if not inner:
            return []
        items = re.findall(r"'[^']*'|\"[^\"]*\"|[^,\s][^,]*", inner)
        return [_parse_nf_value(x) for x in items]
    if (v.startswith("'") and v.endswith("'")) or (v.startswith('"') and v.endswith('"')):
        return v[1:-1]
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def _strip_nf_comment(v: str) -> str:
    """Drop a trailing // comment that is outside quotes (URLs inside quotes survive)."""
    q = None
    for i, ch in enumerate(v):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif v.startswith("//", i):
            return v[:i].strip()
    return v.strip()


def parse_nextflow_config(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("//"):
            continue
        m = re.match(r"params\.([A-Za-z0-9_]+)\s*=\s*(.+)$", s)
        if not m:
            continue
        out[m.group(1)] = _parse_nf_value(_strip_nf_comment(m.group(2)))
    return out


def cmd_config(args: argparse.Namespace) -> int:
    cfg = parse_nextflow_config(Path(args.config).read_text())
    resolved = dict(CONFIG_DEFAULTS)
    resolved.update(cfg)
    missing = [k for k in REQUIRED_CONFIG if k not in cfg]
    issues = []
    ft, fn = resolved.get("FASTQ_FILES_TRANSCRIPTS", []), resolved.get("FASTQ_NAMES_TRANSCRIPTS", [])
    if isinstance(ft, list) and isinstance(fn, list) and len(ft) != len(fn):
        issues.append(f"FASTQ_FILES_TRANSCRIPTS ({len(ft)}) and FASTQ_NAMES_TRANSCRIPTS ({len(fn)}) differ in length")
    gt, gn = resolved.get("FASTQ_FILES_GUIDES", []), resolved.get("FASTQ_NAMES_GUIDES", [])
    if isinstance(gt, list) and isinstance(gn, list) and len(gt) != len(gn):
        issues.append(f"FASTQ_FILES_GUIDES ({len(gt)}) and FASTQ_NAMES_GUIDES ({len(gn)}) differ in length")
    for n in (fn if isinstance(fn, list) else []):
        if len(str(n).split("_")) < 2 or not str(n).split("_")[1].startswith("L"):
            issues.append(f"sample name {n!r} is not <sample>_L<lane> (the lane is parsed from it)")
    try:
        chem, _ = resolve_assay(str(resolved.get("CHEMISTRY", "")))
        resolved["_chemistry_resolved"] = chem
    except Exception as e:
        issues.append(f"CHEMISTRY: {e}")
    if resolved.get("RUN_MULTISEQ") and not all(k in cfg for k in ("R1_MULTI", "R2_MULTI", "BARCODES_MULTIBAR_LIST_MULTI")):
        issues.append("RUN_MULTISEQ is true but R1_MULTI / R2_MULTI / BARCODES_MULTIBAR_LIST_MULTI are not all set")
    if "WHITELIST" not in cfg:
        issues.append("WHITELIST is not set (main.nf passes params.WHITELIST to both kb count calls)")
    out = {"params": cfg, "resolved": resolved, "missing_required": missing, "issues": issues, "valid": not missing}
    print(json.dumps(out, indent=2, default=str))
    if args.out:
        write_json(out, Path(args.out))
    return 0 if not missing else 1


def cmd_inspect(args: argparse.Namespace) -> int:
    pd = _pd()
    mods = load_bundle(Path(args.muon_data))
    info = {k: {"n_obs": int(v.n_obs), "n_vars": int(v.n_vars), "obs": list(v.obs.columns), "var": list(v.var.columns),
                "layers": list(v.layers.keys())} for k, v in mods.items()}
    if "guides" in mods:
        g = mods["guides"]
        te = g.var["target_elements"] if "target_elements" in g.var else pd.Series([str(x).split("|")[0] for x in g.var.index])
        info["target_elements"] = te.value_counts().to_dict()
        info["guide_types"] = pd.Series([guide_type(str(x)) for x in te]).value_counts().to_dict()
        if args.target:
            sel = [str(x) for x, t in zip(g.var.index, te) if t == args.target]
            X = dense(g.X)
            idx = [list(map(str, g.var.index)).index(s) for s in sel]
            cells = (X[:, idx] > 0).any(1) if idx else []
            info["target"] = {"target": args.target, "guides": sel, "cells_with_any_guide": int(sum(cells))}
            if "scRNA" in mods and len(idx):
                r = mods["scRNA"]
                genes = args.genes or []
                if genes:
                    pos = [list(map(str, r.var.index)).index(x) for x in genes if x in set(map(str, r.var.index))]
                    Y = dense(r.X[:, pos])
                    info["target"]["mean_expression_with_vs_without"] = {
                        str(r.var.index[p]): [float(Y[cells, i].mean()) if cells.any() else None, float(Y[~cells, i].mean())]
                        for i, p in enumerate(pos)}
    print(json.dumps(info, indent=2, default=str))
    return 0


# ---------------------------------------------------------------------------
# run: the chain
# ---------------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    pd = _pd()
    t0 = time.time()
    d = run_dir(args.label)
    NS = argparse.Namespace
    steps = {}
    guide_dirs = list(args.guide_dirs or [])
    if args.guide_fastqs:
        if not args.guide_table:
            raise SystemExit("--guide-fastqs needs --guide-table")
        names = args.guide_names or [f"S1_L{i + 1}" for i in range(len(args.guide_fastqs))]
        feats = read_guide_features(Path(args.guide_table))
        chem, _ = resolve_assay(args.chemistry)
        for nm, pair in zip(names, args.guide_fastqs):
            fq = [Path(x) for x in re.split(r"[,\s]+", pair.strip()) if x]
            a, st = count_guides(feats, fq, chem, read_whitelist(args.whitelist))
            out = d / "counts" / f"{nm}_ks_guide_out" / "counts_unfiltered"
            out.mkdir(parents=True, exist_ok=True)
            a.write_h5ad(out / "adata.h5ad")
            print(f"Wrote: {out / 'adata.h5ad'}")
            guide_dirs.append(str(out.parent))
            steps.setdefault("count_guides", {})[nm] = st
    man = lane_manifest(list(args.rna_dirs) + guide_dirs)
    write_tsv(man, d / "initial_preprocessing_file_names.txt")
    pf = NS(manifest=str(d / "initial_preprocessing_file_names.txt"), rna_dirs=[], guide_dirs=[],
            expected_cell_number=args.expected_cell_number, mito_expected_percentage=args.mito_expected_percentage,
            transcripts_umi_threshold=args.transcripts_umi_threshold, mito_prefix=args.mito_prefix,
            percentage_of_cells_to_include_transcript=args.percentage_of_cells_to_include_transcript,
            upstream_compat=args.upstream_compat, no_plots=args.no_plots, seed=args.seed, label="")
    cmd_prefilter(pf, d)
    steps["prefilter"] = pf._prefilter
    an = NS(ann_guide=pf._prefilter["guides"], ann_exp=pf._prefilter["scrna"], gtf=args.gtf, upstream_compat=args.upstream_compat)
    cmd_annotate(an, d)
    steps["annotate"] = an._annotate
    bz = NS(muon_data=an._annotate["bundle"], guide_umi_limit=args.guide_umi_limit, merge=str(bool(args.merge)),
            upstream_compat=args.upstream_compat, no_plots=args.no_plots)
    cmd_binarize(bz, d)
    steps["binarize"] = bz._binarize
    final = bz._binarize["bundle"]
    if args.multiseq_r1:
        ms = NS(muon_data=final, r1=args.multiseq_r1, r2=args.multiseq_r2, barcodes=args.multiseq_barcodes,
                bar=args.bar, umi=args.umi, tag=args.tag, upstream_compat=args.upstream_compat, no_plots=args.no_plots)
        cmd_multiseq(ms, d)
        steps["multiseq"] = ms._multiseq
        final = ms._multiseq["bundle"]
    pl = NS(muon_data=final, gtf=args.gtf, distance_from_guide=args.distance_neighbors, in_trans=args.in_trans,
            add_gene_names=args.add_gene_names, upstream_compat=args.upstream_compat, seed=args.seed)
    cmd_perturb_loader(pl, d)
    steps["perturb_loader"] = pl._loader
    de = NS(perturbdata=pl._loader["perturbdata"], muon_data=final, engine=args.engine, direction=args.direction, B=args.B,
            seed=args.seed, min_cells_per_guide=args.min_cells_per_guide)
    cmd_de(de, d)
    steps["de"] = de._de
    rs = NS(sceptre_dir=de._de["sceptre_out"], muon_data=final, alpha_guide=args.alpha_guide, alpha_element=args.alpha_element)
    cmd_results(rs, d)
    steps["results"] = rs._results
    tr = NS(results=rs._results["result_guides"], muon_data=final, alpha=args.alpha_tracks, all_pairs=False,
            no_plots=args.no_plots, region=None)
    cmd_tracks(tr, d)
    steps["tracks"] = tr._tracks
    lanes = pd.DataFrame(steps["prefilter"]["lanes"])
    res = pd.read_csv(rs._results["result_guides"], sep="\t")
    el = pd.read_csv(rs._results["result_elements"], sep="\t")
    figs = [f for s in steps.values() if isinstance(s, dict) for f in s.get("figures", [])]
    secs = ["## Cell QC per lane", "", md_table(list(lanes.columns), lanes.values), "",
            f"Cells after QC: {steps['prefilter']['cells']}; genes kept: {steps['prefilter']['genes']} "
            f"(detected in >= {steps['prefilter']['genes_min_cells']} cells).", "",
            "## Guides", "", f"{steps['binarize']['guides']} guides, binarised at UMI > {args.guide_umi_limit}; "
            f"mean {steps['binarize']['mean_guides_per_cell']:.2f} guides per cell.", ""]
    if "multiseq" in steps:
        secs += ["## MULTI-seq", "", md_table(["class", "cells"], sorted(steps["multiseq"]["calls"].items())), ""]
    secs += ["## Elements", "", md_table(["element", "GUIDE_TYPE", "anchor", "n_guides", "n_genes_tested"],
                                         [[e["element"], e["GUIDE_TYPE"], e["anchor"], e["n_guides"], e["n_genes_tested"]] for e in steps["perturb_loader"]["elements"]]), "",
             f"## Differential perturbation ({args.engine}, side {args.direction})", "",
             f"{len(res)} guide x gene tests; {int(res['significant'].sum())} at BH adj < {args.alpha_guide}.", "",
             md_table(["gRNA_id", "gene_id", "p_value", "adj_pvalue", "log_fold_change"],
                      res.sort_values("p_value")[["gRNA_id", "gene_id", "p_value", "adj_pvalue", "log_fold_change"]].values, 20), "",
             "## Element-level (Fisher)", "", md_table(["element", "gene_id", "p_combined", "n_guides"],
                                                       el.sort_values("p_combined")[["element", "gene_id", "p_combined", "n_guides"]].values, 15)]
    summ = {"steps": steps, "figures": figs, "seconds": round(time.time() - t0, 1), "upstream_compat": bool(args.upstream_compat),
            "final_bundle": final, "results_bundle": rs._results["bundle"]}
    finish(d, f"Perturb-seq pipeline: {args.label}", secs, summ)
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def _rand_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def _mut(rng, s, k=1):
    s = list(s)
    for i in rng.choice(len(s), k, replace=False):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def synthetic_world(d: Path, seed: int = 11) -> dict:
    """Two lanes, 256 genes (4 cell types), 10 guides (3 TSS targets, 1 enhancer, 1 scrambled control)."""
    pd, np, ad = _pd(), _np(), _ad()
    from scipy import sparse  # type: ignore
    rng = np.random.default_rng(seed)
    genes = []
    for i in range(150):
        genes.append(("1", 100_000 * i + 50_000, "+" if i % 2 == 0 else "-", f"GENE{i:03d}"))
    for i in range(100):
        genes.append(("2", 100_000 * i + 50_000, "+" if i % 2 == 0 else "-", f"GENE{150 + i:03d}"))
    for i in range(6):
        genes.append(("MT", 1000 * i + 500, "+", f"MT-G{i}"))
    gtf = d / "genes.gtf"
    with open(gtf, "w") as fh:
        fh.write("#!genome-build synthetic\n")
        for k, (c, tss, st, nm) in enumerate(genes):
            s, e = (tss, tss + 5000) if st == "+" else (tss - 5000, tss)
            gid = f"ENSG{k:011d}"
            fh.write(f'{c}\tsynth\tgene\t{s}\t{e}\t.\t{st}\t.\tgene_id "{gid}"; gene_version "1"; gene_name "{nm}"; gene_biotype "protein_coding";\n')
            fh.write(f'{c}\tsynth\ttranscript\t{s}\t{e}\t.\t{st}\t.\tgene_id "{gid}"; transcript_id "ENST{k:011d}"; gene_name "{nm}";\n')
    G = len(genes)
    tss_of = {nm: (c, t) for c, t, _, nm in genes}
    # guides: 3 TSS targets, 1 enhancer 30 kb from GENE060, a scrambled control (2 guides each)
    enh_pos = tss_of["GENE060"][1] + 30_000
    targets = [("GENE010_TSS", "chr1", tss_of["GENE010"][1]), ("GENE030_TSS", "chr1", tss_of["GENE030"][1]),
               ("GENE170_TSS", "chr2", tss_of["GENE170"][1]), ("chr1_enh60", "chr1", enh_pos), ("scrambled", "chr3", 91_021_957)]
    seqs = []
    while len(seqs) < 10:
        s = _rand_seq(rng, 20)
        if all(sum(a != b for a, b in zip(s, t)) >= 4 for t in seqs):
            seqs.append(s)
    rows = []
    k = 0
    for t, c, p in targets:
        for j in range(2):
            rows.append({"sgRNA_ID": f"{t}_sg{j + 1}", "Target_name": t, "sgRNA_sequences": seqs[k], "chr": c,
                         "start": p + 20 * j, "end": p + 20 * j + 20})
            k += 1
    gt = pd.DataFrame(rows)
    gt.to_csv(d / "guides.tsv", sep="\t", index=False)
    proc = process_guide_table(gt)
    gid_of_row = dict(zip(proc["sgRNA_ID"], proc["pipeline_id"]))
    guide_ids = [gid_of_row[r["sgRNA_ID"]] for r in rows]
    effects = {}  # guide index -> (gene index, multiplier)
    gidx = {nm: i for i, (_, _, _, nm) in enumerate(genes)}
    for gi, r in enumerate(rows):
        t = r["Target_name"]
        if t.endswith("_TSS"):
            effects[gi] = (gidx[t.replace("_TSS", "")], 0.25)
        elif t == "chr1_enh60":
            effects[gi] = (gidx["GENE060"], 0.5)
    base = rng.lognormal(0, 1.0, G)
    base[[gidx[x] for x in ("GENE010", "GENE030", "GENE170", "GENE060")]] = 6.0
    mt = np.array([nm.startswith("MT-") for _, _, _, nm in genes])
    base[mt] = 3.0
    base = base / base.sum()
    ctype_prog = np.ones((4, G))
    marker_pool = [i for i in range(150, 250) if genes[i][3] not in ("GENE170",)]
    for ct in range(4):
        ctype_prog[ct, marker_pool[ct * 20:(ct + 1) * 20]] = 8.0
    whitelist = set()
    lanes = []
    truth = {"doublet_barcodes": set(), "mito_barcodes": set(), "real": set(), "guides": {}}
    for lane in (1, 2):
        n_real, n_empty = 700, 300
        bcs = []
        while len(bcs) < n_real + n_empty:
            b = _rand_seq(rng, 16)
            if b not in whitelist:
                whitelist.add(b); bcs.append(b)
        ct = rng.integers(0, 4, n_real)
        guides = (rng.random((n_real, 10)) < 0.12)
        depth = rng.lognormal(math.log(3000), 0.3, n_real)
        is_dbl = rng.random(n_real) < 0.06
        partner = rng.integers(0, n_real, n_real)
        is_mito = (rng.random(n_real) < 0.05) & ~is_dbl
        X = np.zeros((n_real + n_empty, G))
        GC = np.zeros((n_real + n_empty, 10))
        for i in range(n_real):
            def profile(ii):
                m = base * ctype_prog[ct[ii]]
                for gj in np.where(guides[ii])[0]:
                    if gj in effects:
                        m[effects[gj][0]] *= effects[gj][1]
                m = m / m.sum()
                return m
            m = profile(i) * depth[i]
            g_on = guides[i].copy()
            if is_dbl[i]:
                m = m + profile(partner[i]) * depth[partner[i]]
                g_on = g_on | guides[partner[i]]
            if is_mito[i]:
                m = m.copy(); m[mt] += m.sum() * 0.5 / mt.sum()
            lam = rng.gamma(5.0, m / 5.0)
            X[i] = rng.poisson(lam)
            GC[i] = rng.poisson(np.where(g_on, 25.0, 0.4))
        for i in range(n_real, n_real + n_empty):
            X[i] = rng.poisson(base * 40)
            GC[i] = rng.poisson(0.3, 10)
        for i in range(n_real):
            truth["real"].add(bcs[i])
            if is_dbl[i]:
                truth["doublet_barcodes"].add(bcs[i])
            if is_mito[i]:
                truth["mito_barcodes"].add(bcs[i])
            truth["guides"][bcs[i]] = guides[i] | (guides[partner[i]] if is_dbl[i] else False)
        var = pd.DataFrame({"feature_name": [nm for _, _, _, nm in genes]}, index=[f"ENSG{k:011d}.1" for k in range(G)])
        rna = ad.AnnData(X=sparse.csr_matrix(X.astype(np.float32)), obs=pd.DataFrame(index=bcs), var=var)
        rd = d / f"S1_L{lane}_ks_transcripts_out" / "counts_unfiltered"
        rd.mkdir(parents=True)
        rna.write_h5ad(rd / "adata.h5ad")
        lanes.append({"lane": lane, "bcs": bcs, "GC": GC, "rna_dir": str(rd.parent)})
    # lane 1 guides as FASTQ (10XV3: 16 bp barcode + 12 bp UMI; R2 = 6 bp + guide + 14 bp)
    L1 = lanes[0]
    r1l, r2l = [], []
    n_mm = n_bcerr = 0
    for i, b in enumerate(L1["bcs"]):
        for gj in range(10):
            for _u in range(int(L1["GC"][i, gj])):
                umi = _rand_seq(rng, 12)
                for _rep in range(int(rng.integers(1, 3))):
                    g = seqs[gj]
                    if rng.random() < 0.1:
                        g = _mut(rng, g); n_mm += 1
                    bb = b
                    if rng.random() < 0.04:
                        bb = _mut(rng, b); n_bcerr += 1
                    r1l.append(bb + umi + "TTTTTTTTTT")
                    r2l.append(_rand_seq(rng, 6) + g + _rand_seq(rng, 14))
    for _ in range(300):  # junk reads
        r1l.append(_rand_seq(rng, 38)); r2l.append(_rand_seq(rng, 40))
    order = rng.permutation(len(r1l))
    for path, seqs_ in ((d / "guide_R1.fastq.gz", r1l), (d / "guide_R2.fastq.gz", r2l)):
        with gzip.open(path, "wt") as fh:
            for n, ix in enumerate(order):
                s = seqs_[ix]
                fh.write(f"@r{n}\n{s}\n+\n{'F' * len(s)}\n")
    decoys = set()
    while len(decoys) < 500:
        decoys.add(_rand_seq(rng, 16))
    (d / "whitelist.txt").write_text("\n".join(sorted(whitelist | decoys)) + "\n")
    # lane 2 guides as an h5ad count matrix
    L2 = lanes[1]
    ga = ad.AnnData(X=sparse.csr_matrix(L2["GC"].astype(np.float32)), obs=pd.DataFrame(index=L2["bcs"]),
                    var=pd.DataFrame({"feature_name": guide_ids}, index=guide_ids))
    g2 = d / "S1_L2_ks_guide_out" / "counts_unfiltered"
    g2.mkdir(parents=True)
    ga.write_h5ad(g2 / "adata.h5ad")
    # MULTI-seq: 4 hashes, cell BC 1-16, UMI 17-28, tag R2 1-8
    hashes = []
    while len(hashes) < 4:
        h = _rand_seq(rng, 8)
        if all(sum(a != b for a, b in zip(h, t)) >= 3 for t in hashes):
            hashes.append(h)
    (d / "multiseq_barcodes.csv").write_text("\n".join(hashes) + "\n")
    m1, m2 = [], []
    ms_truth = {}
    for L in lanes:
        for b in L["bcs"][:700]:
            u = rng.random()
            if u < 0.08:
                cls, hot = "Doublet", list(rng.choice(4, 2, replace=False))
            elif u < 0.14:
                cls, hot = "Negative", []
            else:
                h = int(rng.integers(0, 4)); cls, hot = hashes[h], [h]
            ms_truth[b] = cls
            for h in range(4):
                n = int(rng.poisson(60 if h in hot else 4))
                for _ in range(n):
                    t = hashes[h] if rng.random() > 0.05 else _mut(rng, hashes[h])
                    m1.append(b + _rand_seq(rng, 12)); m2.append(t + _rand_seq(rng, 22))
    for path, s_ in ((d / "ms_R1.fastq.gz", m1), (d / "ms_R2.fastq.gz", m2)):
        with gzip.open(path, "wt") as fh:
            for n, s in enumerate(s_):
                fh.write(f"@m{n}\n{s}\n+\n{'F' * len(s)}\n")
    return {"gtf": gtf, "guides": d / "guides.tsv", "guide_ids": guide_ids, "seqs": seqs, "rows": rows, "lanes": lanes,
            "truth": truth, "ms_truth": ms_truth, "hashes": hashes, "n_mm_reads": n_mm, "n_bcerr_reads": n_bcerr,
            "whitelist": d / "whitelist.txt", "r1": d / "guide_R1.fastq.gz", "r2": d / "guide_R2.fastq.gz",
            "guide_dir_l2": str(g2.parent), "ms_r1": d / "ms_R1.fastq.gz", "ms_r2": d / "ms_R2.fastq.gz",
            "ms_barcodes": d / "multiseq_barcodes.csv", "genes": genes}


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    pd, np, ad = _pd(), _np(), _ad()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    existed = OUT_ROOT.exists()
    before = set(OUT_ROOT.glob("*")) if existed else set()
    NS = argparse.Namespace
    t0 = time.time()
    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            W = synthetic_world(d)
            print(f"\nsynthetic world built ({time.time() - t0:.1f}s)")

            print("\nassay specification")
            t = parse_technology("0,0,16:0,16,26:0,26,0,1,0,0")
            check(t["bc"] == [(0, 0, 16)] and t["umi"] == [(0, 16, 26)] and t["seq"] == [(0, 26, 0), (1, 0, 0)],
                  "parse_technology: bc/umi/seq triplets, multi-part sequence (5' PE config string)")
            check(resolve_assay("10XV3")[0] == "0,0,16:0,16,28:1,0,0" and resolve_assay("10xv2")[1].endswith("737K-august-2016.txt")
                  and resolve_assay("custom", "0,0,16:0,16,26:1,0,0", "wl.txt") == ("0,0,16:0,16,26:1,0,0", "wl.txt"),
                  "resolve_assay: 10XV3 / 10XV2 table, 'custom' passes chemistry + whitelist through")

            print("\nread composition")
            fq = d / "comp.fastq"
            fq.write_text("@a\nACGT\n+\nIIII\n@b\nACGA\n+\nIIII\n")
            c1 = composition_sequences(read_fastq_head_lines(fq), False)
            c2 = composition_sequences(read_fastq_head_lines(fq), True)
            b = compositional_bias(c1)
            check(c1 == ["ACGT", "ACGA"] and len(c2) == 4 and b.loc[0, "A"] == 2 and abs(b.loc[3, "sd_bias"] - np.std([1 - 2, 0 - 2, 1 - 2, 0 - 2], ddof=1)) < 1e-12,
                  "composition: sequence lines only (compat keeps quality lines without @/+/F); sd of the 4 counts")

            print("\nguide table")
            gdir = d / "gt"
            rc = cmd_guide_table(NS(guides=str(W["guides"]), out_dir=str(gdir), label="st"))
            feats = read_guide_features(gdir / "guide_features.txt")
            ids = [i for _, i in feats]
            check(rc == 0 and len(feats) == 10 and ids[0].startswith("GENE010_TSS|1_sgrna_chr1:") and ids[1].startswith("GENE010_TSS|2_")
                  and ids[-1].startswith("scrambled|2_sgrna_chr3:91021977") and ids == sorted(ids, key=lambda x: x.split("|")[0]),
                  "guide-table: <target>|<n>_sgrna_<chr>:<start>:<end>, groups sorted by target")
            ex, mm, _, ks = build_kite_map([("AAAA", "a"), ("AATT", "b")])
            check("AATA" not in mm and "CAAA" in mm and mm["CAAA"] == 0 and ks["n_collisions_removed"] >= 2,
                  "kite map: Hamming-1 variants shared by two guides are removed")

            print("\ncount-guides (lane 1 FASTQ)")
            cg = d / "cg"
            rc = cmd_count_guides(NS(guides=str(gdir / "guide_features.txt"), r1=str(W["r1"]), r2=str(W["r2"]), fastq=None,
                                     chemistry="10XV3", whitelist=str(W["whitelist"]), name="S1_L1", out_dir=str(cg), label=None, max_reads=0))
            st = json.loads((cg / "S1_L1_ks_guide_out" / "kite_stats.json").read_text())
            a1 = ad.read_h5ad(cg / "S1_L1_ks_guide_out" / "counts_unfiltered" / "adata.h5ad")
            L1 = W["lanes"][0]
            truth = pd.DataFrame(L1["GC"], index=L1["bcs"], columns=W["guide_ids"])
            got = pd.DataFrame(dense(a1.X), index=a1.obs.index, columns=a1.var.index).reindex(index=truth.index, columns=truth.columns).fillna(0)
            agree = float((got.values == truth.values).mean())
            check(rc == 0 and agree > 0.995, f"count-guides: UMI counts equal the planted counts in {agree:.4f} of cell x guide entries")
            check(st["mapped_mismatch"] > 0.8 * W["n_mm_reads"] * 0.9 and st["bc_corrected"] > 0.8 * W["n_bcerr_reads"] and st["unmapped"] >= 290,
                  f"count-guides: {st['mapped_mismatch']} Hamming-1 guide reads rescued, {st['bc_corrected']} barcodes whitelist-corrected, {st['unmapped']} junk reads unmapped")

            print("\nmap-rna (kb not on PATH -> command printed)")
            rc = cmd_map_rna(NS(fastq=["R1.fq.gz", "R2.fq.gz"], chemistry="10XV2", whitelist="wl.txt", name="S1_L1", index=None, t2g=None,
                                reference="human", threads=4, memory="48G", workflow=None, out_dir=str(d / "mr"), label=None, dry_run=which("kb") is None))
            s = json.loads((d / "mr" / "summary.json").read_text())
            check(rc == 0 and "kb ref -d human" in s["commands"][0] and "-x 0,0,16:0,16,26:1,0,0" in s["commands"][1],
                  "map-rna: kb ref + kb count commands with the resolved chemistry")

            man = lane_manifest([L["rna_dir"] for L in W["lanes"]] + [str(cg / "S1_L1_ks_guide_out"), W["guide_dir_l2"]])
            check(list(man["sample"]) == ["scRNA", "scRNA", "Guide", "Guide"] and list(man["lane"]) == ["1", "2", "1", "2"],
                  "manifest: sample from 'guide' in the path, lane from <sample>_L<lane>")

            print("\nrun: prefilter -> annotate -> binarize -> multiseq -> perturb-loader -> de -> results -> tracks")
            ra = NS(label="st_bulkcrispr", guide_table=None, guide_fastqs=None, guide_names=None, chemistry="10XV3", whitelist=None,
                    rna_dirs=[L["rna_dir"] for L in W["lanes"]], guide_dirs=[str(cg / "S1_L1_ks_guide_out"), W["guide_dir_l2"]],
                    gtf=str(W["gtf"]), expected_cell_number=650, mito_expected_percentage=0.2, transcripts_umi_threshold=100,
                    mito_prefix="MT-", percentage_of_cells_to_include_transcript=0.01, guide_umi_limit=5, merge=False,
                    multiseq_r1=str(W["ms_r1"]), multiseq_r2=str(W["ms_r2"]), multiseq_barcodes=str(W["ms_barcodes"]),
                    bar=[1, 16], umi=[17, 28], tag=[1, 8], distance_neighbors=1_000_000, in_trans="FALSE", add_gene_names="",
                    engine="sceptre", direction="both", B=199, min_cells_per_guide=30, alpha_guide=0.01, alpha_element=0.05,
                    alpha_tracks=0.05, upstream_compat=False, no_plots=args.no_plots, seed=0)
            rc = cmd_run(ra)
            rd = sorted(OUT_ROOT.glob("*_st_bulkcrispr*"))[-1]
            summ = json.loads((rd / "summary.json").read_text())
            S = summ["steps"]
            check(rc == 0 and (rd / "report.md").is_file(), f"run: report + summary ({summ['seconds']}s)")
            pre = ad.read_h5ad(S["prefilter"]["scrna"])
            cells = set(pre.obs.index)
            tr = W["truth"]
            check(not (cells - tr["real"]), f"prefilter: every empty droplet removed ({len(cells)} cells kept)")
            mito_left = len(cells & tr["mito_barcodes"]) / max(len(tr["mito_barcodes"]), 1)
            check(mito_left < 0.1, f"prefilter: planted high-mito cells removed ({mito_left:.2%} left)")
            dbl_rm = 1 - len(cells & tr["doublet_barcodes"]) / len(tr["doublet_barcodes"])
            singlets = tr["real"] - tr["doublet_barcodes"] - tr["mito_barcodes"]
            lanes_df = pd.DataFrame(S["prefilter"]["lanes"])
            knee_kept = int(lanes_df["after_knee"].sum())
            fp = 1 - len(cells & singlets) / len(singlets)
            check(dbl_rm >= 0.4 and fp < 0.2, f"prefilter: Scrublet ({lanes_df['doublet_engine'].iloc[0]}) removes {dbl_rm:.0%} of planted doublets; "
                  f"{fp:.1%} of singlets lost overall (knee keeps {knee_kept})")
            check((lanes_df["after_knee"] <= 660).all() and (lanes_df["after_min_genes"] <= 700).all(),
                  "prefilter: min_genes = UMI threshold drops the empties; knee[EXPECTED_CELL_NUMBER] caps each lane")
            check(S["prefilter"]["genes_min_cells"] == int(len(cells) * 0.01), "prefilter: genes detected in >= int(n_cells * 0.01) cells")
            mods = load_bundle(Path(S["binarize"]["bundle"]))
            gb, rb = mods["guides"], mods["scRNA"]
            tg = np.array([tr["guides"][b] for b in gb.obs.index]).astype(int)
            order = [W["guide_ids"].index(x) for x in gb.var.index]
            acc = float((dense(gb.X) == tg[:, order]).mean())
            check(acc > 0.98, f"binarize: guide calls (UMI > 5) match the planted guides in {acc:.3f} of entries")
            check(list(rb.obs.columns) == COVARIATES and np.allclose(rb.obs["log_total_gene_count"], np.log(pre.obs.loc[rb.obs.index, "n_counts"].astype(float) + 1)),
                  "binarize: covariates bath_number / percent_mito / log_* ; log_total_gene_count = log(n_counts + 1)")
            ann = load_bundle(Path(S["annotate"]["bundle"]))["scRNA"].var
            gneg = [x for x in W["genes"] if x[2] == "-" and not x[3].startswith("MT")][0]
            row = ann[ann["feature_name"] == gneg[3]].iloc[0]
            gv = load_bundle(Path(S["annotate"]["bundle"]))["guides"].var
            check(int(row["transcript_start"]) == gneg[1] and int(row["transcript_end"]) == gneg[1] + 1
                  and gv.loc[W["guide_ids"][0], "guide_chr"] == "chr1" and gv.loc[W["guide_ids"][0], "target_elements"] == "GENE010_TSS"
                  and gv.loc[W["guide_ids"][1], "guide_number"] == "2",
                  "annotate: minus-strand TSS = gene end, end = TSS + 1; guide_chr / guide_number / target_elements parsed")
            ms = S["multiseq"]
            fin = pd.read_csv(sorted(rd.glob("final_class.csv"))[0])
            cb = pd.read_csv(sorted(rd.glob("final_class_cell_barcode.csv"))[0])
            calls = pd.Series(fin["x"].values, index=cb["x"].values)
            mt = pd.Series({b: W["ms_truth"][b] for b in calls.index})
            sing = mt.isin(W["hashes"])
            acc_s = float((calls[sing] == mt[sing]).mean())
            dbl = mt == "Doublet"
            rec_d = float((calls[dbl] == "Doublet").mean()) if dbl.any() else 1.0
            negs = mt == "Negative"
            rec_n = float((calls[negs] == "Negative").mean()) if negs.any() else 1.0
            check(acc_s > 0.9 and rec_d > 0.6 and rec_n > 0.6, f"multiseq: singlet hash accuracy {acc_s:.3f}, doublet recall {rec_d:.2f}, negative recall {rec_n:.2f}")
            final = load_bundle(Path(summ["final_bundle"]))
            check(not set(final["multiseq"].obs["multiseq_class"]) & {"Doublet", "Negative"} and final["guides"].n_obs == ms["kept_cells"],
                  f"multiseq: Doublet and Negative cells removed from all modalities ({ms['kept_cells']} kept)")
            els = {e["element"]: e for e in S["perturb_loader"]["elements"]}
            check(els["GENE010_TSS"]["GUIDE_TYPE"] == "POSITIVE_CONTROL" and els["GENE010_TSS"]["n_genes_tested"] == 19
                  and els["chr1_enh60"]["GUIDE_TYPE"] == "PUTATIVE_ENHANCER" and els["chr1_enh60"]["anchor"].startswith("guide:chr1")
                  and els["scrambled"]["GUIDE_TYPE"] == "NEGATIVE_CONTROL" and els["scrambled"]["anchor"] == "random_10",
                  "perturb-loader: TSS +/- 1 Mb gene set (19 genes), enhancer anchored on its guide, control gets 10 random genes")
            res = pd.read_csv(S["results"]["result_guides"], sep="\t")
            tgt = res[res["target_elements"].str.endswith("_TSS") & (res["gene_id"] == res["target_elements"].str.replace("_TSS", "", regex=False))]
            check(len(tgt) == 6 and (tgt["p_value"] < 1e-3).all() and (tgt["log_fold_change"] < -0.7).all(),
                  f"de: all 6 TSS guides knock down their gene (max p {tgt['p_value'].max():.1e}, LFC {tgt['log_fold_change'].max():.2f} .. {tgt['log_fold_change'].min():.2f})")
            enh = res[(res["target_elements"] == "chr1_enh60") & (res["gene_id"] == "GENE060")]
            check(len(enh) == 2 and (enh["p_value"] < 0.01).all() and (enh["log_fold_change"] < -0.3).all(),
                  f"de: enhancer guides reduce GENE060 (p {', '.join(f'{p:.1e}' for p in enh['p_value'])})")
            nul = res[(res["target_elements"] == "scrambled")]
            off = res[res["target_elements"].str.endswith("_TSS") & ~res.index.isin(tgt.index)]
            f_nul, f_off = float((nul["p_value"] < 0.05).mean()), float((off["p_value"] < 0.01).mean())
            check(len(nul) == 20 and f_nul <= 0.2 and f_off <= 0.05,
                  f"de: controls calibrated (scrambled p < 0.05 in {f_nul:.0%}; off-target TSS pairs p < 0.01 in {f_off:.1%})")
            try:
                from statsmodels.stats.multitest import multipletests  # type: ignore
                ok_bh = np.allclose(multipletests(res["p_value"], method="fdr_bh")[1], res["adj_pvalue"])
            except ImportError:
                ok_bh = True
            el = pd.read_csv(S["results"]["result_elements"], sep="\t")
            e10 = el[(el["element"] == "GENE010_TSS") & (el["gene_id"] == "GENE010")].iloc[0]
            rb_ = load_bundle(Path(summ["results_bundle"]))
            check(ok_bh and e10["p_combined"] < 1e-6 and e10["n_guides"] == 2 and bool(e10["sig_not_adj"])
                  and set(rb_) == {"guides", "scRNA", "result_guides", "result_elements"} and "adj_pvalue" in rb_["result_guides"].layers,
                  "results: BH adj_pvalue (= statsmodels fdr_bh), Fisher element p, 4-modality results bundle with layers")
            check(S["tracks"]["links"] >= 7 and (Path(S["tracks"]["tracks_dir"]) / "tracks.ini").is_file(),
                  f"tracks: {S['tracks']['links']} guide -> TSS links + BED + pyGenomeTracks ini")
            if not args.no_plots:
                check(len(summ["figures"]) >= 10, f"run: {len(summ['figures'])} figures")

            print("\nupstream-compat switches")
            cov_c = make_covariates(pre, pd.DataFrame(np.zeros((pre.n_obs, 1))), upstream_compat=True)
            check(np.allclose(cov_c["log_total_gene_count"], np.log(pre.obs["n_genes"].astype(float) + 1)),
                  "compat: log_total_gene_count = log(n_genes + 1) as upstream")
            genes = read_gtf_genes(W["gtf"])
            _, rc_ = annotate_modalities(ad.AnnData(np.zeros((1, 1)), var=pd.DataFrame(index=[W["guide_ids"][0]])),
                                         ad.AnnData(np.zeros((1, 1)), var=pd.DataFrame(index=[f"ENSG{W['genes'].index(gneg):011d}"])), genes, True)
            gtf_end = int(genes.loc[genes["gene_name"] == gneg[3], "end"].iloc[0])
            check(int(rc_.var["transcript_end"].iloc[0]) == gtf_end, "compat: transcript_end = gene end (the 'end ' column typo)")
            _, _, info_c = perturb_loader(final, genes, upstream_compat=True)
            check(info_c.set_index("element").loc["chr1_enh60", "anchor"] == "random_10", "compat: enhancer elements get 10 random genes as upstream")
            bar = pd.read_csv(sorted(rd.glob("bar_table.csv"))[0], index_col=0)
            cc, _ = multiseq_classify(bar, upstream_compat=True)
            check(multiseq_class_columns(bar, True)[-1] == "nUMI_total" and "nUMI_total" not in multiseq_class_columns(bar, False)
                  and "Doublet" not in multiseq_drop_classes(True) and len(cc) == len(bar),
                  "compat: nUMI_total classified as a barcode and 'Double' filter (doublets kept) as upstream")

            print("\nstatistics helpers")
            check(local_maxima([1, 3, 2, 5, 4]) == [1, 3] and local_maxima([2, 2, 1, 3, 1]) == [3], "local_maxima: deMULTIplex rle semantics")
            gx, dx = bkde_normal(np.random.default_rng(1).normal(size=500))
            check(abs(np.trapz(dx, gx) - 1) < 0.02, "bkde_normal: density integrates to 1")
            bim = np.r_[np.random.default_rng(2).normal(0, 1, 2000), np.random.default_rng(3).normal(8, 1, 2000)]
            th = threshold_minimum(bim)
            check(2 < th < 6, f"threshold_minimum: {th:.2f} between the modes")

            print("\nassign-guides")
            rc = cmd_assign_guides(NS(muon_data=S["binarize"]["bundle"], guide_umi_limit=5, label="st_bulkcrispr_assign"))
            ag = sorted(OUT_ROOT.glob("*_st_bulkcrispr_assign*"))[-1]
            gm = load_bundle(ag / "mu_with_binary")["guides"]
            acc_m = float((np.asarray(gm.layers["binarized"]) == tg[:, order]).mean())
            check(rc == 0 and acc_m > 0.98, f"assign-guides: Poisson-mixture calls match planted guides in {acc_m:.3f}")

            print("\nhelpers")
            rc = cmd_cellranger_inputs(NS(guides=str(W["guides"]), read="R2", pattern="(BC)", rna_fastq_dir="rna", guide_fastq_dir="gd",
                                          rna_sample="rna", guide_sample="gd", transcriptome="/ref", execute=False, label="st_bulkcrispr_cr"))
            cr = sorted(OUT_ROOT.glob("*_st_bulkcrispr_cr*"))[-1]
            fr = pd.read_csv(cr / "feature_ref.csv")
            check(rc == 0 and len(fr) == 10 and set(fr["feature_type"]) == {"CRISPR Guide Capture"} and (cr / "library.csv").is_file(),
                  "cellranger-inputs: feature_ref.csv (10 guides) + library.csv")
            csvp = d / "fq.csv"
            (d / "x.fastq.gz").write_bytes(b"")
            pd.DataFrame([{"Sample": "S", "Sample_Number": 1, "Lane_Number": 2, "Read_Type": "R1", "File_Fragment": 1, "file_dir": str(d / "x.fastq.gz")}]).to_csv(csvp, index=False)
            rc = cmd_rename_fastq(NS(csv=str(csvp), out_dir=str(d / "ren"), label="x"))
            check(rc == 0 and (d / "ren" / "symlink_dir" / "S_S001_L002_R1_001.fastq.gz").is_symlink(), "rename-fastq: <Sample>_S001_L002_R1_001.fastq.gz symlink")
            cfgp = d / "perturb.config"
            cfgp.write_text("params.GTF_GZ_LINK = 'http://ftp.ensembl.org/x.gtf.gz'\n// comment\nparams.CHEMISTRY = '10XV3'\n"
                            "params.THREADS = 15 // cpus\nparams.IN_TRANS = \"FALSE\"\nparams.CREATE_REF = false\n"
                            "params.FASTQ_FILES_TRANSCRIPTS = ['a_R1.fq.gz a_R2.fq.gz', 'b_R1.fq.gz b_R2.fq.gz']\n"
                            "params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1', 'S1_L2']\nparams.BAR_MULTI= [1,16]\n")
            cfg = parse_nextflow_config(cfgp.read_text())
            check(cfg["GTF_GZ_LINK"] == "http://ftp.ensembl.org/x.gtf.gz" and cfg["THREADS"] == 15 and cfg["CREATE_REF"] is False
                  and cfg["FASTQ_FILES_TRANSCRIPTS"] == ["a_R1.fq.gz a_R2.fq.gz", "b_R1.fq.gz b_R2.fq.gz"] and cfg["BAR_MULTI"] == [1, 16],
                  "config: Nextflow params (URLs, // comments, lists, booleans, ints)")
            rc = cmd_config(NS(config=str(cfgp), out=None))
            check(rc == 1, "config: missing required params -> exit 1")
            rc = cmd_inspect(NS(muon_data=str(summ["final_bundle"]), target="GENE010_TSS", genes=["GENE010"]))
            check(rc == 0, "inspect: modality summary + target subset")
    finally:
        after = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
        for p in after - before:
            shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink()
        if not existed and OUT_ROOT.exists() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.1f}s)")
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent bulk-crispr", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("assay-spec", help="chemistry string + whitelist for an assay (Task 1)")
    p.add_argument("--assay", required=True, help="10XV2, 10XV3, 5PE, custom, or a kallisto bc:umi:seq string")
    p.add_argument("--chemistry", help="with --assay custom: kallisto technology string")
    p.add_argument("--whitelist", help="with --assay custom: whitelist path")
    p.set_defaults(func=cmd_assay_spec)

    p = sub.add_parser("composition", help="per-position nucleotide bias of R1/R2 (fq_composition.py)")
    p.add_argument("--r1", required=True)
    p.add_argument("--r2")
    p.add_argument("--n-lines", type=int, default=10000)
    p.add_argument("--upstream-compat", action="store_true", help="upstream line filter (quality lines can leak in)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="composition")
    p.set_defaults(func=cmd_composition)

    p = sub.add_parser("guide-table", help="guide xlsx/tsv -> guide_features.txt (guide_table_processing.py)")
    p.add_argument("--guides", required=True, help="sgRNA_ID, Target_name, sgRNA_sequences, chr, start, end")
    p.add_argument("--out-dir")
    p.add_argument("--label", default="guide_table")
    p.set_defaults(func=cmd_guide_table)

    p = sub.add_parser("count-guides", help="kite-style guide counting from FASTQ (exact + Hamming-1)")
    p.add_argument("--guides", required=True, help="guide_features.txt or the guide table")
    p.add_argument("--r1")
    p.add_argument("--r2")
    p.add_argument("--fastq", nargs="+", help="FASTQ files in technology-string order (instead of --r1/--r2)")
    p.add_argument("--chemistry", default="10XV3", help="assay name or kallisto bc:umi:seq string")
    p.add_argument("--whitelist", help="barcode whitelist (plain or .gz); Hamming-1 correction when given")
    p.add_argument("--name", default="S1_L1", help="<sample>_L<lane>")
    p.add_argument("--max-reads", type=int, default=0)
    p.add_argument("--out-dir")
    p.add_argument("--label")
    p.set_defaults(func=cmd_count_guides)

    p = sub.add_parser("map-rna", help="kb ref + kb count for the cDNA library (runs if kb is on PATH)")
    p.add_argument("--fastq", nargs="+", required=True)
    p.add_argument("--chemistry", default="10XV3")
    p.add_argument("--whitelist")
    p.add_argument("--name", default="S1_L1")
    p.add_argument("--reference", default="human", help="kb ref -d <reference>")
    p.add_argument("--index", help="existing kallisto index (CUSTOM_REFERENCE_IDX)")
    p.add_argument("--t2g", help="existing t2g (CUSTOM_REFERENCE_T2T)")
    p.add_argument("--threads", type=int, default=10)
    p.add_argument("--memory", default="48G")
    p.add_argument("--workflow", help="kb --workflow (upstream passes kite)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--out-dir")
    p.add_argument("--label")
    p.set_defaults(func=cmd_map_rna)

    p = sub.add_parser("manifest", help="count directories -> lane manifest (preprocessing.py)")
    p.add_argument("--dirs", nargs="+", required=True)
    p.add_argument("--out")
    p.add_argument("--label", default="manifest")
    p.set_defaults(func=cmd_manifest)

    def qc_args(p):
        p.add_argument("--expected-cell-number", type=int, default=8000)
        p.add_argument("--mito-expected-percentage", type=float, default=0.2)
        p.add_argument("--transcripts-umi-threshold", type=int, default=100, help="applied as the minimum genes per cell, as upstream")
        p.add_argument("--mito-prefix", default="MT-", help="mitochondrial gene-name prefix (mt- for mouse)")
        p.add_argument("--percentage-of-cells-to-include-transcript", type=float, default=0.01)
        p.add_argument("--seed", type=int, default=0)

    p = sub.add_parser("prefilter", help="per-lane QC, doublets, lane concatenation, gene filter")
    p.add_argument("--manifest", help="initial_preprocessing_file_names.txt")
    p.add_argument("--rna-dirs", nargs="+", default=[])
    p.add_argument("--guide-dirs", nargs="+", default=[])
    qc_args(p)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="prefilter")
    p.set_defaults(func=cmd_prefilter)

    p = sub.add_parser("annotate", help="guide coordinates + gene TSS (muon_creation.py)")
    p.add_argument("--ann-guide", required=True)
    p.add_argument("--ann-exp", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="annotate")
    p.set_defaults(func=cmd_annotate)

    p = sub.add_parser("binarize", help="guide merge + binarisation + covariates (merge_bin_and_muon.py)")
    p.add_argument("--muon-data", required=True, help="annotated bundle dir or .h5mu")
    p.add_argument("--guide-umi-limit", type=int, default=5)
    p.add_argument("--merge", default="false")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="binarize")
    p.set_defaults(func=cmd_binarize)

    def ms_args(p):
        p.add_argument("--bar", type=int, nargs=2, default=[1, 16], help="cell barcode positions in R1 (1-based, inclusive)")
        p.add_argument("--umi", type=int, nargs=2, default=[17, 28])
        p.add_argument("--tag", type=int, nargs=2, default=[1, 8], help="MULTI-seq tag positions in R2")

    p = sub.add_parser("multiseq", help="MULTI-seq demultiplexing (deMULTIplex port)")
    p.add_argument("--muon-data", required=True)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--barcodes", required=True, help="MULTI-seq barcode list (first CSV column)")
    ms_args(p)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="multiseq")
    p.set_defaults(func=cmd_multiseq)

    p = sub.add_parser("perturb-loader", help="elements, guide membership, cis gene sets (PerturbLoader_generation.py)")
    p.add_argument("--muon-data", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--distance-from-guide", type=int, default=1_000_000)
    p.add_argument("--in-trans", default="FALSE")
    p.add_argument("--add-gene-names", default="")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="perturb_loader")
    p.set_defaults(func=cmd_perturb_loader)

    def de_args(p):
        p.add_argument("--engine", choices=["sceptre", "mannwhitney", "sceptre-r"], default="sceptre")
        p.add_argument("--direction", choices=["left", "right", "both"], default="both")
        p.add_argument("--B", type=int, default=500, help="conditional resamples per guide")
        p.add_argument("--min-cells-per-guide", type=int, default=MIN_CELLS_PER_GUIDE, help="guides need > this many cells")

    p = sub.add_parser("de", help="differential perturbation per element (runSceptre.py)")
    p.add_argument("--perturbdata", required=True, help="perturbdata/ dir from perturb-loader")
    p.add_argument("--muon-data")
    de_args(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="de")
    p.set_defaults(func=cmd_de)

    p = sub.add_parser("results", help="BH + Fisher aggregation into a results bundle (sceptre_anndata_creation.py)")
    p.add_argument("--sceptre-dir", required=True)
    p.add_argument("--muon-data", required=True)
    p.add_argument("--alpha-guide", type=float, default=0.01)
    p.add_argument("--alpha-element", type=float, default=0.05)
    p.add_argument("--label", default="results")
    p.set_defaults(func=cmd_results)

    p = sub.add_parser("tracks", help="BED + pyGenomeTracks links of guide -> gene tests (Task 4)")
    p.add_argument("--results", required=True, help="result_guides.tsv")
    p.add_argument("--muon-data", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--all-pairs", action="store_true")
    p.add_argument("--region", help="chr:start-end to render with pyGenomeTracks when installed")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="tracks")
    p.set_defaults(func=cmd_tracks)

    p = sub.add_parser("assign-guides", help="depth-aware Poisson-mixture guide calls (Task guide-assign)")
    p.add_argument("--muon-data", required=True)
    p.add_argument("--guide-umi-limit", type=int, default=5)
    p.add_argument("--label", default="assign_guides")
    p.set_defaults(func=cmd_assign_guides)

    p = sub.add_parser("cellranger-inputs", help="library.csv + feature_ref.csv for cellranger count (Task cellranger)")
    p.add_argument("--guides", required=True)
    p.add_argument("--read", default="R2")
    p.add_argument("--pattern", default="(BC)")
    p.add_argument("--rna-fastq-dir", required=True)
    p.add_argument("--guide-fastq-dir", required=True)
    p.add_argument("--rna-sample", required=True)
    p.add_argument("--guide-sample", required=True)
    p.add_argument("--transcriptome", required=True)
    p.add_argument("--execute", action="store_true", help="run cellranger when it is on PATH")
    p.add_argument("--label", default="cellranger_inputs")
    p.set_defaults(func=cmd_cellranger_inputs)

    p = sub.add_parser("rename-fastq", help="bcl2fastq-style names + symlinks (rename_and_symlink.py)")
    p.add_argument("--csv", required=True, help="Sample, Sample_Number, Lane_Number, Read_Type, File_Fragment, file_dir")
    p.add_argument("--out-dir")
    p.add_argument("--label", default="rename_fastq")
    p.set_defaults(func=cmd_rename_fastq)

    p = sub.add_parser("config", help="parse + validate a Nextflow perturb.config")
    p.add_argument("--config", required=True)
    p.add_argument("--out")
    p.set_defaults(func=cmd_config)

    p = sub.add_parser("inspect", help="summarise a modality bundle / .h5mu (Task 2)")
    p.add_argument("--muon-data", required=True)
    p.add_argument("--target")
    p.add_argument("--genes", nargs="+")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("run", help="prefilter -> annotate -> binarize -> (multiseq) -> perturb-loader -> de -> results -> tracks")
    p.add_argument("--gtf", required=True)
    p.add_argument("--rna-dirs", nargs="+", required=True, help="<sample>_L<lane>_ks_transcripts_out dirs (kb count output)")
    p.add_argument("--guide-dirs", nargs="+", default=[], help="<sample>_L<lane>_ks_guide_out dirs")
    p.add_argument("--guide-fastqs", nargs="+", help="per lane 'R1,R2' (counted here instead of --guide-dirs)")
    p.add_argument("--guide-names", nargs="+", help="lane names for --guide-fastqs (S1_L1 ...)")
    p.add_argument("--guide-table")
    p.add_argument("--chemistry", default="10XV3")
    p.add_argument("--whitelist")
    qc_args(p)
    p.add_argument("--guide-umi-limit", type=int, default=5)
    p.add_argument("--merge", action="store_true")
    p.add_argument("--multiseq-r1")
    p.add_argument("--multiseq-r2")
    p.add_argument("--multiseq-barcodes")
    ms_args(p)
    p.add_argument("--distance-neighbors", type=int, default=1_000_000)
    p.add_argument("--in-trans", default="FALSE")
    p.add_argument("--add-gene-names", default="")
    de_args(p)
    p.add_argument("--alpha-guide", type=float, default=0.01)
    p.add_argument("--alpha-element", type=float, default=0.05)
    p.add_argument("--alpha-tracks", type=float, default=0.05)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="perturbseq")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic two-lane screen through every subcommand")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
