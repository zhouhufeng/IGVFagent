#!/usr/bin/env python3
"""IGVF single-cell-like CRISPR pipeline on the Gasperini 2019 screen (port of IGVF-CRISPR/Pipeline_Gasperini_2019).

Port of https://github.com/IGVF-CRISPR/Pipeline_Gasperini_2019 (No LICENSE file;
no licence declared), the IGVF CRISPR focus group's processing of the
Gasperini et al. 2019 (Cell 176:377) high-MOI CRISPRi pilot screen with the
"PERTURB-SEQ single-cell like" Nextflow pipeline
(LucasSilvaFerreira/pipeline_perturbseq_like@60cd461, which the upstream
README clones and configures; Perturb_Loader supplies the element/guide/gene
tables).  The upstream repository holds the README (download script, conda
environment, pipeline config), df_from_gasperini_tss.xlsx (98 promoter guides
for 49 genes) and the resulting sample MuData
(mudata/Gasperini_2019_sample_pilot.h5mu, 7,314 cells x 98 guides x 2,127
genes).  Every step the README runs was read -- the Nextflow processes and
their bin/ scripts -- and re-derived in Python (numpy / pandas / scipy /
anndata); R (SCEPTRE, deMULTIplex), kallisto/kb and scrublet become optional
engines.  Relationship: port (methods and thresholds reproduced, no code
copied).

Definitions, as the pipeline computes them (Gasperini config values)
  guide table     group by Target_name; the n-th guide of a target becomes
                  "<Target>|<n>" and pipeline_id "<Target>|<n>_sgrna_<chr>:<start>:<end>";
                  guide_features.txt = sgRNA_sequences <tab> pipeline_id.
  composition     first 10,000 FASTQ lines; per read position the counts of
                  A, C, T, G and their standard deviation (ddof 1) -- the
                  "compositional bias" track.
  cell QC (per lane)  knee = total UMIs sorted descending; cells need
                  >= TRANSCRIPTS_UMI_TRHESHOLD (200) and >= knee[EXPECTED_CELL_NUMBER]
                  (10,000) total UMIs; percent_mito = counts in 'MT-' genes /
                  total, keep < MITO_EXPECTED_PERCENTAGE (0.2); scrublet
                  doublets (expected rate 0.1, simulated ratio 2, k =
                  round(0.5 sqrt(n))) removed; guide number_of_nonzero_guides
                  = guides with > 0 UMIs; barcodes kept only if present in
                  both modalities; obs batch_number = lane.
  lane merge      anndata.concat(merge="same"), obs_names_make_unique; genes
                  kept if detected in >= int(n_cells * 0.01) cells.
  MuData          guides.var: guide_chr = name.split('|')[1].split(':')[0].split('_')[-1],
                  guide_start / guide_end = 2nd-last / last ':' field,
                  guide_number = name.split('|')[1].split('_')[0],
                  target_elements = name.split('|')[0]; scRNA.var:
                  transcript_chr / transcript_start (TSS = gene start on +,
                  end on -) from the GTF gene record (Ensembl id without
                  version), 'NOT_FOUND', 0, 0 otherwise; IGVF additions
                  guides.var['sequences'], uns['elements'],
                  varm['guide_by_element'] as in the sample MuData.
  guide assignment  (optional merge: sum guides of a target) binary =
                  UMI > GUIDE_UMI_LIMIT (3; strict).
  covariates      bath_number (sic), percent_mito, log_number_of_detected_genes
                  = log(n_genes), log_total_gene_count = log(n_counts + 1),
                  log_total_guide_count = log(assigned guides + 1); batch
                  dropped when it has a single level.
  tested pairs    elements = guide-name prefix before '|'; GUIDE_TYPE
                  POSITIVE_CONTROL ('_TSS'), NEGATIVE_CONTROL ('random' /
                  'scrambled'), PUTATIVE_ENHANCER ('chr'); an element's genes
                  are the expressed genes on the same chromosome whose TSS is
                  within DISTANCE_NEIGHBORS (1 Mb) of the TSS of the gene the
                  element is named after ('_TSS' stripped); IN_TRANS = all
                  genes; elements that are not a gene get 10 random genes
                  (NUMBER_FOR_RANDOM_CONTROLS); ADDGENENAMES are appended.
  SCEPTRE         per element, guides with > 30 assigned cells are tested
                  against every tested gene (pair_type 'all_elements_test'),
                  side = DIRECTION ('both').  Python engine: conditional
                  randomisation test -- NB GLM of the gene on the covariates
                  (Poisson start, theta by ML, fixed-theta IRLS), logistic
                  propensity of the guide on the covariates, NB score
                  statistic z (efficient score, nuisance-adjusted), B = 500
                  guide vectors resampled from the propensities, skew-normal
                  fitted to the null z (method of moments; empirical p if the
                  fit fails), p = left / right / 2 x min tail;
                  log_fold_change = 1-D NB MLE of the guide effect with the
                  null fit as offset.  R engine: writes upstream's per-element
                  inputs + run_sceptre_high_moi script and runs Rscript when
                  available.
  results         adj_pvalue = Benjamini-Hochberg over all guide-gene tests;
                  significant = adj_pvalue < 0.01; element-gene p = Fisher's
                  combination of its guides' p-values; sig_not_adj = p < 0.05;
                  MuData modalities result_guides / result_elements.
  Gasperini test  the original paper's model: NB GLM per gRNA group x gene,
                  y ~ indicator(cell has any guide of the group) + guide_count
                  + percent_mito + batch, offset log size factor
                  (total / geometric mean of totals), per-gene theta by ML
                  under the null, likelihood-ratio test (chi-square 1 df);
                  columns as GEO GSE120861 all_deg_results (pairs4merge, beta,
                  intercept, fold_change.transcript_remaining = exp(beta),
                  pvalue.raw, pvalue.empirical vs NTC groups, ...).
  comparison      join to GSE120861_all_deg_results.pilot on gene symbol and
                  gRNA_group = element + '_TSS': Spearman of -log10 p, sign
                  agreement of effects, BH hits at FDR 0.1 in each (2 x 2 +
                  Fisher exact), self-TSS (positive control) recovery.
  MULTI-seq       (deMULTIplex, RUN_MULTISEQ) R1 cell barcode 1-16, UMI 17-28,
                  R2 tag 1-8; tags matched to the reference with Hamming <= 1;
                  bar table = unique UMIs per cell x tag (+ nUMI, nUMI_total);
                  classifyCells: log2, mean-centre, Gaussian KDE (bw.nrd0) on
                  100 points between the 0.1% and 99.9% quantiles, local maxima
                  -> threshold = q-quantile between the lowest-peak maximum and
                  the highest maximum; cells over one threshold = that tag,
                  more = Doublet, none = Negative; findThresh sweeps q =
                  0.01..0.99 by 0.02 and takes the q maximising singlets; two
                  rounds, negatives of round 1 removed; Doublet/Negative cells
                  dropped from the MuData.

Upstream quirks corrected by default (reproduced with --upstream-compat)
  * concact_and_pre_filtering.py filters sc.pp.filter_cells(min_genes=cutoff)
    with the UMI threshold, i.e. >= 200 detected GENES; the port applies it as
    >= 200 UMIs, as the parameter name and the log message say.
  * the lane-merge gene filter ignores --percentage_of_cells_to_include_transcript
    (always 1%); the port honours it (same value by default).
  * merge_bin_and_muon.py computes log_total_gene_count = log(n_genes + 1);
    the port uses log(n_counts + 1).
  * muon_creation.py writes the gene END into transcript_end (it sets a column
    named 'end ' to TSS + 1 and then reads 'end'); the port stores TSS + 1.
  * PerturbLoader parses guide coordinates but then looks the element up as a
    gene, so enhancer elements get 10 random genes; the port uses the guide
    coordinates for elements that are not genes (still random genes for
    non-targeting controls).
  * multiseq.py drops `ncol(dt)` -- base R's t density, so only nUMI is removed
    and nUMI_total is classified as if it were a barcode; the port removes both.

Subcommands
  setup          fetch the upstream guide table + sample MuData (28 MB) and the
                 GEO GSE120861 pilot results (19 MB) into Data/GasperiniPipeline;
                 --whitelist adds 737K-august-2016.txt; --print-raw prints the
                 SRA BAM download + bamtofastq commands (several GB, never run).
  guide-table    xlsx/tsv guide table -> guide_features.txt + pipeline ids.
  composition    per-position nucleotide composition of a FASTQ (+ figure).
  kb             kb ref / kb count (kite) commands for transcripts and guides.
  config         the Nextflow gasperini_sample.config + nextflow run command.
  qc-filter      per-lane cell QC + lane merge from kb h5ad outputs.
  mudata         guide + scRNA AnnData -> raw MuData with coordinates.
  assign         guide binarisation, covariates -> processed MuData.
  pairs          element x guide and element x tested-gene tables.
  sceptre        per-guide cis tests (python CRT or R SCEPTRE).
  results        BH / Fisher aggregation -> mudata_results.h5mu.
  gasperini-test the original Gasperini 2019 NB likelihood-ratio test.
  compare        concordance with the original Gasperini results table.
  multiseq       MULTI-seq sample demultiplexing (deMULTIplex port).
  inspect        check a MuData against the sample-pilot schema.
  run            assign -> pairs -> sceptre -> results [-> gasperini-test -> compare]
                 on a MuData (default: the fetched sample pilot).
  selftest       synthetic two-lane screen with planted knockdowns, doublets,
                 MULTI-seq tags; every subcommand asserted.

Output: Docs/GasperiniPipeline/<timestamp>_<label>/ (report.md, summary.json,
tables, figures).  Needs numpy, pandas, scipy, anndata, h5py; scikit-learn for
doublets; matplotlib optional; mudata / muon / scrublet / R used when present.

Usage:
    igvfagent gasperini-pipeline setup
    igvfagent gasperini-pipeline run --label pilot
    igvfagent gasperini-pipeline run --h5mu Data/GasperiniPipeline/Gasperini_2019_sample_pilot.h5mu --gasperini-test --compare Data/GasperiniPipeline/GSE120861_all_deg_results.pilot.txt.gz
    igvfagent gasperini-pipeline guide-table --table df_from_gasperini_tss.xlsx
    igvfagent gasperini-pipeline qc-filter --rna S1_L1_ks_transcripts_out --guide S1_L1_ks_guide_out --lanes 1
    igvfagent gasperini-pipeline mudata --rna-h5ad full_raw_scrna_ann_data.h5ad --guide-h5ad full_raw_guide_ann_data.h5ad --gtf Homo_sapiens.GRCh38.106.gtf.gz
    igvfagent gasperini-pipeline inspect --h5mu mudata.h5mu
    igvfagent gasperini-pipeline selftest --no-plots
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
from typing import Any, Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "GasperiniPipeline"
DATA_ROOT = ROOT / "Data" / "GasperiniPipeline"

UPSTREAM_REPO = "IGVF-CRISPR/Pipeline_Gasperini_2019"
UPSTREAM_COMMIT = "51b4ced480d312278c8531f4b4710d98448171c4"
PIPELINE_REPO = "LucasSilvaFerreira/pipeline_perturbseq_like"
PIPELINE_COMMIT = "60cd461ae8d1788dfda4a04f686a42b796a79559"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/"
SETUP_FILES = {  # name -> (url, approx bytes, what)
    "df_from_gasperini_tss.xlsx": (RAW + "df_from_gasperini_tss.xlsx", 13983, "98 promoter guides for 49 genes (upstream)"),
    "Gasperini_2019_sample_pilot.h5mu": (RAW + "mudata/Gasperini_2019_sample_pilot.h5mu", 28459747, "sample pilot MuData (upstream)"),
    "GSE120861_all_deg_results.pilot.txt.gz": (GEO + "GSE120861_all_deg_results.pilot.txt.gz", 18819923, "original Gasperini pilot DE results"),
    "GSE120861_grna_groups.pilot.txt.gz": (GEO + "GSE120861_grna_groups.pilot.txt.gz", 38546, "pilot gRNA groups"),
}
WHITELIST = ("737K-august-2016.txt", "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt", 12_533_760)
RAW_DATA = [  # the README's download script (several GB: printed, never run)
    ("pilot_highmoi_screen.1_SI_GA_G1.bam", "https://sra-pub-src-1.s3.amazonaws.com/SRR7967482/pilot_highmoi_screen.1_SI_GA_G1.bam.1", "bam_pilot_scrna_1"),
    ("pilot_highmoi_screen.1_CGTTACCG.grna.bam", "https://sra-pub-src-1.s3.amazonaws.com/SRR7967488/pilot_highmoi_screen.1_CGTTACCG.grna.bam.1", "bam_pilot_guide_1"),
]
BAMTOFASTQ = "https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux"
GTF_URL = "http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
GENOME_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
CONFIG = {  # the README's gasperini_sample.config
    "TRANSCRIPTOME_REFERENCE": "human", "CHEMISTRY": "10XV2", "THREADS": 15, "DISTANCE_NEIGHBORS": 1_000_000,
    "IN_TRANS": "FALSE", "EXPECTED_CELL_NUMBER": 10000, "MITO_SPECIE": "hsapiens", "MITO_EXPECTED_PERCENTAGE": 0.2,
    "PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS": 0.01, "TRANSCRIPTS_UMI_TRHESHOLD": 200, "GUIDE_UMI_LIMIT": 3,
    "MERGE": False, "DIRECTION": "both", "RUN_MULTISEQ": False, "BAR_MULTI": [1, 16], "UMI_MULTI": [17, 28], "R2_MULTI_TAG": [1, 8],
}
NUMBER_FOR_RANDOM_CONTROLS = 10
SCEPTRE_MIN_CELLS = 30
PILOT_SCHEMA = {
    "guides": {"obs": ["number_of_nonzero_guides", "batch_number"],
               "var": ["feature_name", "guide_chr", "guide_end", "guide_start", "guide_number", "target_elements", "sequences"],
               "uns": ["elements"], "varm": ["guide_by_element"]},
    "scRNA": {"obs": ["n_genes", "n_counts", "percent_mito", "doublet_scores", "predicted_doublets", "doublet_info", "batch_number"],
              "var": ["feature_name", "n_cells", "transcript_chr", "transcript_start", "transcript_end"], "uns": [], "varm": []},
}
RESULT_COLS = ["gene_id", "gRNA_id", "pair_type", "p_value", "z_value", "log_fold_change", "p_empirical", "skew_fit_success",
               "n_cells_guide", "element", "guide_type", "engine"]
GASPERINI_COLS = ["pairs4merge", "beta", "intercept", "fold_change.transcript_remaining", "pvalue.raw", "pvalue.empirical",
                  "pvalue.empirical.adjusted", "gene_short_name", "gRNA_group", "site_type", "n_cells_group"]

BLUE, ORANGE, INK, AXIS, GREY = "#2a78d6", "#eb6834", "#0b0b0b", "#d0cfca", "#96a0b3"

log = logging.getLogger("gasperini_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"gasperini_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    if d.exists():
        d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_{os.getpid()}_{int(time.time() * 1e3) % 1000}"
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


def _sp():
    import scipy.sparse as sp  # type: ignore
    return sp


def _ad():
    try:
        import anndata  # type: ignore
        return anndata
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs anndata + h5py") from e


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _elem_io():
    try:
        from anndata.io import read_elem, write_elem  # type: ignore  # anndata >= 0.11
    except ImportError:
        from anndata.experimental import read_elem, write_elem  # type: ignore
    return read_elem, write_elem


def write_tsv(df, path: Path, index: bool = False) -> Path:
    df.to_csv(path, sep="\t", index=index)
    print(f"TSV: {path}")
    return path


def _finite_mean(v) -> float:
    np = _np()
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else float("nan")


def _clean(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def write_json(path: Path, obj: Any) -> Path:
    path.write_text(json.dumps(_clean(obj), indent=2, default=_json_default))
    print(f"JSON: {path}")
    return path


def _json_default(o):
    np = _np()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def write_report(path: Path, lines: List[str]) -> Path:
    path.write_text("\n".join(lines) + "\n")
    print(f"Report: {path}")
    return path


def save_fig(fig, path: Path) -> None:
    fig.savefig(path, dpi=130, bbox_inches="tight")
    print(f"Figure: {path}")


def _dense(X):
    return X.toarray() if hasattr(X, "toarray") else _np().asarray(X)


def _fmt(x, nd=3) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not math.isfinite(x):
        return "NA"
    return f"{x:.{nd}g}" if abs(x) < 1e-3 or abs(x) >= 1e4 else f"{x:.{nd}f}"


def which(binary: str) -> Optional[str]:
    return shutil.which(binary)


# ---------------------------------------------------------------------------
# MuData I/O (mudata when installed, else the h5mu layout written directly)
# ---------------------------------------------------------------------------

def read_h5mu(path) -> Dict[str, Any]:
    try:
        import mudata  # type: ignore
        m = mudata.read_h5mu(str(path))
        return {k: m.mod[k] for k in m.mod}
    except ImportError:
        pass
    import h5py  # type: ignore
    read_elem, _ = _elem_io()
    out = {}
    with h5py.File(str(path), "r") as f:
        order = list(f["mod"].attrs.get("mod-order", list(f["mod"].keys())))
        for k in order:
            k = k.decode() if isinstance(k, bytes) else str(k)
            out[k] = read_elem(f["mod"][k])
    return out


def write_h5mu(path, mods: Dict[str, Any]) -> Path:
    path = Path(path)
    try:
        import mudata  # type: ignore
        mudata.MuData(mods).write(str(path))
        print(f"Wrote: {path}")
        return path
    except ImportError:
        pass
    import h5py  # type: ignore
    np, pd = _np(), _pd()
    _, write_elem = _elem_io()
    names = list(mods)
    obs_names: List[str] = []
    seen = set()
    for a in mods.values():
        for o in a.obs_names:
            if o not in seen:
                seen.add(o)
                obs_names.append(o)
    var_names = [v for a in mods.values() for v in a.var_names]
    with h5py.File(str(path), "w", userblock_size=512) as f:
        for k, a in mods.items():
            write_elem(f, f"mod/{k}", a)
        f["mod"].attrs["mod-order"] = np.array(names, dtype=h5py.string_dtype())
        write_elem(f, "obs", pd.DataFrame(index=pd.Index(obs_names)))
        write_elem(f, "var", pd.DataFrame(index=pd.Index(var_names)))
        oidx = {o: i for i, o in enumerate(obs_names)}
        obsm, obsmap, varm, varmap = {}, {}, {}, {}
        off = 0
        for k, a in mods.items():
            pos = np.array([oidx[o] for o in a.obs_names])
            m = np.zeros(len(obs_names), bool)
            m[pos] = True
            mp = np.zeros(len(obs_names), np.uint32)
            mp[pos] = np.arange(1, a.n_obs + 1, dtype=np.uint32)
            obsm[k], obsmap[k] = m, mp
            vm = np.zeros(len(var_names), bool)
            vm[off:off + a.n_vars] = True
            vp = np.zeros(len(var_names), np.uint32)
            vp[off:off + a.n_vars] = np.arange(1, a.n_vars + 1, dtype=np.uint32)
            varm[k], varmap[k] = vm, vp
            off += a.n_vars
        for key, d in (("obsm", obsm), ("obsmap", obsmap), ("varm", varm), ("varmap", varmap)):
            write_elem(f, key, d)
        write_elem(f, "uns", {})
        f.attrs.update({"encoding-type": "MuData", "encoding-version": "0.1.0", "encoder": "igvfagent", "encoder-version": "1", "axis": 0})
    with open(path, "br+") as fh:
        fh.write(b"MuData (format-version=0.1.0;creator=igvfagent;creator-version=1)")
    print(f"Wrote: {path}")
    return path


def read_counts(path) -> Any:
    """kb count output dir (counts_unfiltered/adata.h5ad) or an .h5ad file."""
    ad = _ad()
    p = Path(path)
    if p.is_dir():
        cand = p / "counts_unfiltered" / "adata.h5ad"
        if not cand.exists():
            raise SystemExit(f"{p} has no counts_unfiltered/adata.h5ad (kb count --h5ad output)")
        p = cand
    a = ad.read_h5ad(str(p))
    if "feature_name" not in a.var.columns:
        a.var["feature_name"] = a.var["gene_name"].astype(str).values if "gene_name" in a.var.columns else a.var_names.astype(str)
    return a


# ---------------------------------------------------------------------------
# Guide table / composition / references
# ---------------------------------------------------------------------------

def guide_features(df):
    """guide_table_processing.py: '<Target>|<n>' names and '<Target>|<n>_sgrna_<chr>:<start>:<end>' pipeline ids."""
    pd = _pd()
    need = {"Target_name", "sgRNA_sequences", "chr", "start", "end"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"guide table lacks columns {sorted(miss)} (needs {sorted(need)})")
    parts = []
    for k, v in df.groupby("Target_name", sort=True):
        v = v.copy()
        v["Target_name"] = [f"{k}|{n + 1}" for n in range(len(v))]
        v["pipeline_id"] = [f"{t}_sgrna_{c}:{s}:{e}" for t, c, s, e in zip(v["Target_name"], v["chr"], v["start"], v["end"])]
        v["element"] = k
        parts.append(v)
    return pd.concat(parts, ignore_index=True)


def read_table(path):
    pd = _pd()
    p = str(path)
    if p.endswith((".xlsx", ".xls")):
        return pd.read_excel(p)
    return pd.read_csv(p, sep="\t" if p.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ",")


def fastq_composition(path, n_lines: int = 10000, upstream_compat: bool = False):
    """fq_composition.py: per-position A/C/T/G counts and their SD across the four nucleotides."""
    pd, np = _pd(), _np()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [fh.readline().rstrip("\n") for _ in range(n_lines)]
    lines = [x for x in lines if x]
    if upstream_compat:  # upstream keeps lines without '@', '+' or 'F' (quality lines without F slip through)
        seqs = [x for x in lines if "@" not in x and "+" not in x and "F" not in x]
    else:
        seqs = [lines[i] for i in range(1, len(lines), 4)]
    L = max((len(s) for s in seqs), default=0)
    arr = np.array([list(s.ljust(L)) for s in seqs]) if seqs else np.zeros((0, 0), dtype="<U1")
    counts = {b: (arr == b).sum(0) for b in "ACTG"}
    df = pd.DataFrame({"position": np.arange(L), **counts})
    df["n_reads"] = len(seqs)
    df["sd"] = df[list("ACTG")].std(axis=1, ddof=1)
    for b in "ACTG":
        df[f"frac_{b}"] = df[b] / max(len(seqs), 1)
    return df


def read_gtf_genes(path):
    """GTF 'gene' records -> gene_id (version stripped), gene_name, chr, start, end, strand, tss."""
    pd = _pd()
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            gid = attrs.get("gene_id", "").split(".")[0]
            start, end = int(f[3]), int(f[4])
            rows.append({"gene_id": gid, "gene_name": attrs.get("gene_name", gid), "chr": f[0], "start": start, "end": end,
                         "strand": f[6], "tss": start if f[6] == "+" else end})
    return pd.DataFrame(rows)


def _chrom(c) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


# ---------------------------------------------------------------------------
# Cell QC (concact_and_pre_filtering.py)
# ---------------------------------------------------------------------------

def _hist_min_threshold(values, nbins: int = 256) -> float:
    """skimage.filters.threshold_minimum: smooth the histogram until bimodal, take the minimum between the peaks."""
    np = _np()
    v = np.asarray(values, float)
    hist, edges = np.histogram(v, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    h = hist.astype(float)
    for _ in range(10000):
        peaks = [i for i in range(1, len(h) - 1) if h[i - 1] < h[i] and h[i + 1] <= h[i]]
        if len(peaks) == 2:
            a, b = peaks
            return float(centers[a + int(np.argmin(h[a:b + 1]))])
        if len(peaks) < 2:
            break
        h = np.convolve(np.r_[h[0], h, h[-1]], np.ones(3) / 3, mode="valid")
    return float(np.quantile(v, 0.9))


def doublet_scores(X, expected_rate: float = 0.1, sim_ratio: float = 2.0, n_pcs: int = 30, seed: int = 0):
    """Scrublet (scrublet package when importable, else a compact port of its kNN doublet score)."""
    np = _np()
    try:
        import scrublet as scr  # type: ignore
        s = scr.Scrublet(X, expected_doublet_rate=expected_rate, sim_doublet_ratio=sim_ratio, random_state=seed)
        sc_, pred = s.scrub_doublets(verbose=False)
        return np.asarray(sc_), np.asarray(pred, bool), "scrublet"
    except ImportError:
        pass
    from sklearn.decomposition import PCA  # type: ignore
    from sklearn.neighbors import NearestNeighbors  # type: ignore
    sp = _sp()
    rng = np.random.default_rng(seed)
    E = sp.csr_matrix(X, dtype=np.float64)
    n, g = E.shape
    tot = np.asarray(E.sum(1)).ravel()
    tot[tot == 0] = 1
    keep = np.asarray(((E >= 3).sum(0))).ravel() >= 3
    En = sp.diags(tot.mean() / tot) @ E
    mu = np.asarray(En.mean(0)).ravel()
    var = np.asarray(En.multiply(En).mean(0)).ravel() - mu ** 2
    fano = np.where(mu > 0, var / np.maximum(mu, 1e-12), 0)
    cand = np.where(keep & (mu > 0))[0]
    if len(cand) < 10:
        cand = np.where(mu > 0)[0]
    genes = cand[fano[cand] >= np.percentile(fano[cand], 85)] if len(cand) > 50 else cand
    n_sim = int(round(n * sim_ratio))
    pa, pb = rng.integers(0, n, n_sim), rng.integers(0, n, n_sim)
    Es = E[pa] + E[pb]
    tsim = np.asarray(Es.sum(1)).ravel()
    tsim[tsim == 0] = 1
    Esn = sp.diags(tot.mean() / tsim) @ Es
    A = _dense(En[:, genes])
    S = _dense(Esn[:, genes])
    m, sd = A.mean(0), A.std(0)
    sd[sd == 0] = 1
    A, S = (A - m) / sd, (S - m) / sd
    k_pc = int(min(n_pcs, A.shape[1] - 1, n - 1))
    pca = PCA(n_components=max(k_pc, 2), random_state=seed).fit(A)
    PA, PS = pca.transform(A), pca.transform(S)
    k = int(round(0.5 * math.sqrt(n)))
    r = n_sim / n
    k_adj = max(int(round(k * (1 + r))), 1)
    nn = NearestNeighbors(n_neighbors=min(k_adj, n + n_sim - 1)).fit(np.vstack([PA, PS]))
    is_sim = np.r_[np.zeros(n, bool), np.ones(n_sim, bool)]

    def score(P, self_in):
        idx = nn.kneighbors(P, n_neighbors=min(k_adj + (1 if self_in else 0), n + n_sim), return_distance=False)
        if self_in:
            idx = idx[:, 1:]
        q = (is_sim[idx].sum(1) + 1) / (idx.shape[1] + 2)
        rho = expected_rate
        return q * rho / r / (1 - rho - q * (1 - rho - rho / r))
    s_obs = score(PA, True)
    s_sim = score(PS, False)
    thr = _hist_min_threshold(s_sim)
    return s_obs, s_obs > thr, "builtin"


def analyze_batch(rna, guide, lane, expected_cells: int, mito_max: float, umi_cutoff: int, upstream_compat: bool = False,
                  run_doublets: bool = True, seed: int = 0) -> Tuple[Any, Any, Dict[str, Any]]:
    np = _np()
    stats: Dict[str, Any] = {"lane": lane, "n_barcodes": int(rna.n_obs)}
    tot = np.asarray(rna.X.sum(1)).ravel()
    knee = np.sort(tot)[::-1]
    stats["n_above_umi_cutoff"] = int((knee > umi_cutoff).sum())
    if upstream_compat:
        ng = np.asarray((rna.X > 0).sum(1)).ravel()
        rna = rna[ng >= umi_cutoff].copy()
    else:
        rna = rna[tot >= umi_cutoff].copy()
    if expected_cells < len(knee):
        kmin = knee[expected_cells]
        t2 = np.asarray(rna.X.sum(1)).ravel()
        rna = rna[t2 >= kmin].copy()
        stats["knee_min_counts"] = float(kmin)
    else:
        stats["knee_min_counts"] = None
        log.warning("lane %s: EXPECTED_CELL_NUMBER %d >= barcodes %d (upstream would raise IndexError); knee filter skipped", lane, expected_cells, len(knee))
    stats["n_after_umi_knee"] = int(rna.n_obs)
    rna.obs["n_genes"] = np.asarray((rna.X > 0).sum(1)).ravel().astype(int)
    rna.var_names = [str(x).split(".")[0] for x in rna.var_names]
    mito = np.array([str(x).startswith("MT-") for x in rna.var["feature_name"]])
    tot = np.asarray(rna.X.sum(1)).ravel().astype(float)
    mt = np.asarray(rna.X[:, np.where(mito)[0]].sum(1)).ravel() if mito.any() else np.zeros(rna.n_obs)
    rna.obs["percent_mito"] = mt / np.maximum(tot, 1)
    rna.obs["n_counts"] = tot
    rna = rna[rna.obs["percent_mito"].values < mito_max].copy()
    stats["n_after_mito"] = int(rna.n_obs)
    if run_doublets and rna.n_obs > 20:
        s, pred, eng = doublet_scores(rna.X, seed=seed)
        stats["doublet_engine"] = eng
    else:
        s, pred = np.zeros(rna.n_obs), np.zeros(rna.n_obs, bool)
        stats["doublet_engine"] = "skipped"
    rna.obs["doublet_scores"] = s
    rna.obs["predicted_doublets"] = pred
    rna.obs["doublet_info"] = pred.astype(str)
    rna = rna[~pred].copy()
    stats["n_doublets"] = int(pred.sum())
    guide = guide.copy()
    guide.obs["number_of_nonzero_guides"] = np.asarray((guide.X > 0).sum(1)).ravel().astype(int)
    shared = set(guide.obs_names) & set(rna.obs_names)
    rna = rna[[c in shared for c in rna.obs_names]].copy()
    guide = guide[[c in shared for c in guide.obs_names]].copy()
    rna.obs["batch_number"] = lane
    guide.obs["batch_number"] = lane
    stats["n_final"] = int(rna.n_obs)
    return rna, guide, stats


def concat_lanes(rnas, guides, pct_cells: float = 0.01, upstream_compat: bool = False):
    ad, np = _ad(), _np()
    r = ad.concat(rnas, merge="same")
    r.obs_names_make_unique()
    g = ad.concat(guides, merge="same")
    g.obs_names_make_unique()
    pct = 0.01 if upstream_compat else pct_cells
    min_cells = int(r.n_obs * pct)
    nc = np.asarray((r.X > 0).sum(0)).ravel()
    r = r[:, nc >= min_cells].copy()
    r.var["n_cells"] = nc[nc >= min_cells].astype(int)
    return r, g


# ---------------------------------------------------------------------------
# MuData creation (muon_creation.py) and guide assignment (merge_bin_and_muon.py)
# ---------------------------------------------------------------------------

def parse_guide_name(name: str) -> Dict[str, str]:
    tail = name.split("|")[1] if "|" in name else name
    head = tail.split(":")[0]
    f = tail.split(":")
    return {"guide_chr": head.split("_")[-1], "guide_end": f[-1], "guide_start": f[-2] if len(f) > 1 else "",
            "guide_number": head.split("_")[0], "target_elements": name.split("|")[0]}


def annotate_guides(guide, sequences: Optional[Dict[str, str]] = None):
    pd, np, sp = _pd(), _np(), _sp()
    parsed = pd.DataFrame([parse_guide_name(str(n)) for n in guide.var_names], index=guide.var_names)
    guide.var["feature_name"] = guide.var_names.astype(str)
    for c in ["guide_chr", "guide_end", "guide_start", "guide_number", "target_elements"]:
        guide.var[c] = parsed[c].values
    guide.var["sequences"] = [(sequences or {}).get(str(n), "") for n in guide.var_names]
    elements = sorted(set(guide.var["target_elements"]))
    col = {e: i for i, e in enumerate(elements)}
    rows = np.arange(guide.n_vars)
    cols = np.array([col[e] for e in guide.var["target_elements"]])
    guide.uns["elements"] = np.array(elements, dtype=object)
    guide.varm["guide_by_element"] = sp.csr_matrix((np.ones(guide.n_vars, np.uint8), (rows, cols)), shape=(guide.n_vars, len(elements)))
    return guide


def annotate_transcripts(rna, genes=None, upstream_compat: bool = False):
    np = _np()
    ids = [str(x).split(".")[0] for x in rna.var_names]
    if genes is None or len(genes) == 0:
        rna.var["transcript_chr"], rna.var["transcript_start"], rna.var["transcript_end"] = "NOT_FOUND", 0, 0
        return rna
    g = genes.drop_duplicates("gene_id").set_index("gene_id")
    chr_, st, en = [], [], []
    for i in ids:
        if i in g.index:
            row = g.loc[i]
            chr_.append(str(row["chr"]))
            st.append(int(row["tss"]))
            en.append(int(row["end"]) if upstream_compat else int(row["tss"]) + 1)
        else:
            chr_.append("NOT_FOUND")
            st.append(0)
            en.append(0)
    rna.var["transcript_chr"] = chr_
    rna.var["transcript_start"] = np.array(st, dtype=np.int64)
    rna.var["transcript_end"] = np.array(en, dtype=np.int64)
    return rna


def binarize_guides(guide, limit: int, merge: bool = False):
    """merge_bin_and_muon.py: optional sum of guides per target, then UMI > GUIDE_UMI_LIMIT."""
    np, sp = _np(), _sp()
    X = sp.csr_matrix(guide.X)
    names = [str(x) for x in guide.var["feature_name"]] if "feature_name" in guide.var else [str(x) for x in guide.var_names]
    if merge:
        keys = [n.split("|")[0] for n in names]
        uniq = sorted(set(keys))
        M = sp.csr_matrix((np.ones(len(keys)), ([uniq.index(k) for k in keys], np.arange(len(keys)))), shape=(len(uniq), len(keys)))
        X = (X @ M.T).tocsr()
        last = {n.split("|")[0]: n for n in names}  # upstream keeps the last guide's coordinates for the merged target
        names = [last[u] for u in uniq]
    B = (X > limit).astype(np.int8)
    return sp.csr_matrix(B), X, names


def make_covariates(rna_obs, n_guides_assigned, upstream_compat: bool = False):
    pd, np = _pd(), _np()
    cov = pd.DataFrame(index=rna_obs.index)
    cov["bath_number"] = rna_obs["batch_number"].values if "batch_number" in rna_obs else 1
    cov["percent_mito"] = rna_obs["percent_mito"].values.astype(float)
    ng = rna_obs["n_genes"].values.astype(float)
    cov["log_number_of_detected_genes"] = np.log(ng)
    base = ng if upstream_compat else rna_obs["n_counts"].values.astype(float)
    cov["log_total_gene_count"] = np.log(base + 1)
    cov["log_total_guide_count"] = np.log(np.asarray(n_guides_assigned, float) + 1)
    return cov


def process_mudata(mods: Dict[str, Any], limit: int, merge: bool = False, upstream_compat: bool = False):
    """Binary guides + covariates in obs of both modalities; scRNA indexed by gene symbol (duplicates dropped)."""
    ad, np = _ad(), _np()
    guide, rna = mods["guides"], mods["scRNA"]
    common = [c for c in rna.obs_names if c in set(guide.obs_names)]
    guide, rna = guide[common].copy(), rna[common].copy()
    B, raw, names = binarize_guides(guide, limit, merge)
    cov = make_covariates(rna.obs, np.asarray(B.sum(1)).ravel(), upstream_compat)
    obs = cov.copy()
    for c in ("n_genes", "n_counts"):
        if c in rna.obs:
            obs[c] = rna.obs[c].values
    rna.var["ensg"] = rna.var_names.astype(str)
    fn = rna.var["feature_name"].astype(str).values
    keep = ~_pd().Series(fn).duplicated().values
    rv = rna.var.iloc[np.where(keep)[0]].copy()
    rv.index = fn[keep]
    rna_p = ad.AnnData(X=rna.X[:, np.where(keep)[0]], obs=obs.copy(), var=rv)
    if merge:
        gvar = guide.var.set_index(guide.var["feature_name"].astype(str)).loc[names].copy()
    else:
        gvar = guide.var.copy()
        gvar.index = names
    g_p = ad.AnnData(X=B, obs=obs.copy(), var=gvar)
    g_p.layers["guide_umi"] = raw
    for k in ("elements",):
        if k in guide.uns:
            g_p.uns[k] = guide.uns[k]
    if not merge and "guide_by_element" in guide.varm:
        g_p.varm["guide_by_element"] = guide.varm["guide_by_element"]
    return {"guides": g_p, "scRNA": rna_p}, cov


# ---------------------------------------------------------------------------
# Tested pairs (PerturbLoader_generation.py)
# ---------------------------------------------------------------------------

def guide_type(e: str) -> Optional[str]:
    if "_TSS" in e:
        return "POSITIVE_CONTROL"
    if "random" in e or "scrambled" in e:
        return "NEGATIVE_CONTROL"
    if "chr" in e:
        return "PUTATIVE_ENHANCER"
    return None


def gene_table_from_rna(rna_var, gtf_genes=None):
    """Expressed genes with chr / TSS: GTF (upstream) or the MuData's transcript_chr / transcript_start."""
    pd = _pd()
    names = list(rna_var.index.astype(str))
    if gtf_genes is not None and len(gtf_genes):
        g = gtf_genes[gtf_genes["gene_name"].isin(set(names))].drop_duplicates("gene_name")
        return pd.DataFrame({"gene_name": g["gene_name"].values, "chr": [_chrom(c) for c in g["chr"]], "tss": g["tss"].astype(int).values})
    v = rna_var[rna_var["transcript_chr"].astype(str) != "NOT_FOUND"]
    return pd.DataFrame({"gene_name": v.index.astype(str), "chr": [_chrom(c) for c in v["transcript_chr"]],
                         "tss": v["transcript_start"].astype(int).values})


def tested_pairs(guide_names: Sequence[str], rna_genes: Sequence[str], genes, distance: int = 1_000_000, in_trans: bool = False,
                 add_genes: Sequence[str] = (), upstream_compat: bool = False, seed: int = 0):
    """element x tested-gene pairs; returns (pairs DataFrame, element->guides dict, element->how dict)."""
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    el_guides: Dict[str, List[str]] = {}
    for g in guide_names:
        el_guides.setdefault(str(g).split("|")[0], []).append(str(g))
    coords = {e: gs[-1].split("_")[-1] for e, gs in el_guides.items()}  # upstream: last guide's "chr:start:end"
    expressed = set(rna_genes)
    gt = genes[genes["gene_name"].isin(expressed)].reset_index(drop=True)
    rows, how = [], {}
    for e in el_guides:
        q = e.replace("_TSS", "")
        if in_trans:
            s, h = set(gt["gene_name"]), "in_trans"
        else:
            hit = gt[gt["gene_name"] == q]
            chrom, pos = None, None
            if len(hit):
                chrom, pos, h = hit["chr"].iloc[0], int(hit["tss"].iloc[0]), "gene_tss"
            elif not upstream_compat and ":" in coords[e]:
                c = coords[e].split(":")
                try:
                    chrom, pos, h = _chrom(c[0]), int(c[1]), "guide_coordinates"
                except (ValueError, IndexError):
                    chrom = None
            if chrom is not None and chrom not in set(gt["chr"]):
                chrom = None  # e.g. non-targeting guides with placeholder coordinates
            if chrom is None:
                k = min(NUMBER_FOR_RANDOM_CONTROLS, len(gt))
                s, h = set(gt["gene_name"].iloc[rng.choice(len(gt), k, replace=False)]) if k else set(), "random_genes"
            else:
                sub = gt[(gt["chr"] == chrom) & (np.abs(gt["tss"] - pos) < distance)]
                s = set(sub["gene_name"])
        s = s | (set(add_genes) & expressed)
        how[e] = h
        for gname in sorted(s):
            rows.append({"element": e, "gene_id": gname, "guide_type": guide_type(e), "selection": h})
    return pd.DataFrame(rows, columns=["element", "gene_id", "guide_type", "selection"]), el_guides, how


# ---------------------------------------------------------------------------
# SCEPTRE-style conditional randomisation test
# ---------------------------------------------------------------------------

def design_matrix(cov, drop_single_batch: bool = True):
    """Intercept + standardised numeric covariates + batch dummies (batch dropped when single-level, as runSceptre.py)."""
    pd, np = _pd(), _np()
    cols = [np.ones(len(cov))]
    names = ["intercept"]
    for c in cov.columns:
        v = cov[c]
        if c in ("bath_number", "batch_number", "batch", "prep_batch"):
            lev = sorted(pd.unique(v.astype(str)))
            if len(lev) <= 1 and drop_single_batch:
                continue
            for lv in lev[1:]:
                cols.append((v.astype(str).values == lv).astype(float))
                names.append(f"{c}={lv}")
            continue
        x = v.values.astype(float)
        x = np.where(np.isfinite(x), x, np.nanmean(x[np.isfinite(x)]) if np.isfinite(x).any() else 0)
        sd = x.std()
        if sd == 0:
            continue
        cols.append((x - x.mean()) / sd)
        names.append(c)
    return np.column_stack(cols), names


def nb_irls(y, X, offset=None, theta: float = float("inf"), beta=None, maxit: int = 30, tol: float = 1e-8):
    np = _np()
    n, p = X.shape
    off = np.zeros(n) if offset is None else offset
    if beta is None:
        beta = np.zeros(p)
        beta[0] = math.log(max(y.mean(), 1e-8)) - (off.mean() if offset is not None else 0)
    for _ in range(maxit):
        eta = np.clip(X @ beta + off, -30, 30)
        mu = np.exp(eta)
        w = mu / (1 + mu / theta) if math.isfinite(theta) else mu
        z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        XtW = X.T * w
        try:
            new = np.linalg.solve(XtW @ X + 1e-8 * np.eye(p), XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(new - beta)) < tol * (1 + np.max(np.abs(beta))):
            beta = new
            break
        beta = new
    mu = np.exp(np.clip(X @ beta + off, -30, 30))
    return beta, mu


def nb_loglik(y, mu, theta: float):
    from scipy.special import gammaln  # type: ignore
    np = _np()
    if not math.isfinite(theta):
        return float(np.sum(y * np.log(np.maximum(mu, 1e-300)) - mu - gammaln(y + 1)))
    return float(np.sum(gammaln(y + theta) - gammaln(theta) - gammaln(y + 1) + theta * np.log(theta / (theta + mu))
                        + y * np.log(np.maximum(mu, 1e-300) / (theta + mu))))


def theta_ml(y, mu, lo: float = 1e-3, hi: float = 1e4) -> float:
    from scipy.optimize import minimize_scalar  # type: ignore
    np = _np()
    if np.var(y) <= np.mean(y):
        return hi
    r = minimize_scalar(lambda lt: -nb_loglik(y, mu, math.exp(lt)), bounds=(math.log(lo), math.log(hi)), method="bounded")
    return float(math.exp(r.x))


def nb_null_fit(y, Z):
    beta, mu = nb_irls(y, Z)
    th = theta_ml(y, mu)
    beta, mu = nb_irls(y, Z, theta=th, beta=beta)
    th = theta_ml(y, mu)
    return mu, th


def logistic_propensity(x, Z, maxit: int = 50):
    np = _np()
    p_ = Z.shape[1]
    beta = np.zeros(p_)
    m = x.mean()
    beta[0] = math.log(max(m, 1e-6) / max(1 - m, 1e-6))
    for _ in range(maxit):
        eta = np.clip(Z @ beta, -30, 30)
        pr = 1 / (1 + np.exp(-eta))
        w = np.maximum(pr * (1 - pr), 1e-10)
        H = (Z.T * w) @ Z + 1e-6 * np.eye(p_)
        step = np.linalg.solve(H, Z.T @ (x - pr))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            break
    return 1 / (1 + np.exp(-np.clip(Z @ beta, -30, 30)))


def skewnorm_mom(null):
    """Method-of-moments skew-normal fit (loc, scale, shape); None when the moments are degenerate."""
    np = _np()
    m, s = float(np.mean(null)), float(np.std(null, ddof=1))
    if not np.isfinite(s) or s <= 0:
        return None
    g = float(np.mean(((null - m) / s) ** 3))
    g = max(min(g, 0.99), -0.99)
    a = abs(g) ** (2 / 3)
    d = math.copysign(math.sqrt((math.pi / 2) * a / (a + ((4 - math.pi) / 2) ** (2 / 3))), g)
    d = max(min(d, 0.995), -0.995)
    alpha = d / math.sqrt(1 - d * d)
    omega = s / math.sqrt(1 - 2 * d * d / math.pi)
    xi = m - omega * d * math.sqrt(2 / math.pi)
    return xi, omega, alpha


def _tail_p(z, null, side: str):
    from scipy.stats import skewnorm  # type: ignore
    np = _np()
    B = len(null)
    emp_l, emp_r = (1 + np.sum(null <= z)) / (B + 1), (1 + np.sum(null >= z)) / (B + 1)
    fit = skewnorm_mom(null)
    ok = fit is not None
    if ok:
        xi, om, al = fit
        pl, pr = float(skewnorm.cdf(z, al, loc=xi, scale=om)), float(skewnorm.sf(z, al, loc=xi, scale=om))
        ok = math.isfinite(pl) and math.isfinite(pr)
    if not ok:
        pl, pr = emp_l, emp_r
    pick = lambda l, r: l if side == "left" else r if side == "right" else min(1.0, 2 * min(l, r))
    return pick(pl, pr), pick(emp_l, emp_r), ok


def nb_effect_mle(y, mu, x, theta: float, maxit: int = 50) -> float:
    """1-D NB MLE of beta in log E[y] = log(mu) + beta * x (x binary)."""
    np = _np()
    idx = x > 0
    if not idx.any():
        return float("nan")
    yy, mm = y[idx], mu[idx]
    if yy.sum() == 0:
        return float("-inf")
    b = math.log(max(yy.sum(), 1e-8) / max(mm.sum(), 1e-12))
    for _ in range(maxit):
        m = mm * math.exp(b)
        den = 1 + m / theta if math.isfinite(theta) else np.ones_like(m)
        score = float(np.sum((yy - m) / den))
        info = float(np.sum(m * (1 + yy / theta) / den ** 2)) if math.isfinite(theta) else float(np.sum(m))
        if info <= 0:
            break
        step = score / info
        b += step
        if abs(step) < 1e-10:
            break
    return float(b)


def crt_guide_genes(Y: Dict[str, Any], x, Z, gene_fits: Dict[str, Tuple[Any, float]], B: int = 500, side: str = "both", seed: int = 0):
    """Test one guide against several genes: observed + B resampled NB score statistics."""
    np = _np()
    rng = np.random.default_rng(seed)
    pr = logistic_propensity(x.astype(float), Z)
    Xt = (rng.random((B, len(x))) < pr[None, :]).astype(np.float64)
    out = []
    for gname, y in Y.items():
        mu, th = gene_fits[gname]
        den = 1 + mu / th if math.isfinite(th) else np.ones_like(mu)
        r = (y - mu) / den
        w = mu / den
        ZtW = Z.T * w
        Hinv = np.linalg.pinv(ZtW @ Z)

        def zstat(xm):
            U = xm @ r
            a = xm @ ZtW.T  # (B, p)
            V = xm @ w - np.einsum("bp,pq,bq->b", a, Hinv, a)
            return U / np.sqrt(np.maximum(V, 1e-12))
        z_obs = float(zstat(x[None, :].astype(np.float64))[0])
        null = zstat(Xt)
        null = null[np.isfinite(null)]
        p, pe, ok = _tail_p(z_obs, null, side)
        out.append({"gene_id": gname, "z_value": z_obs, "p_value": p, "p_empirical": pe, "skew_fit_success": ok,
                    "log_fold_change": nb_effect_mle(y, mu, x, th)})
    return out


def run_sceptre_python(rna, guides, cov, pairs, el_guides, B: int = 500, side: str = "both", min_cells: int = SCEPTRE_MIN_CELLS,
                       seed: int = 0, max_elements: Optional[int] = None):
    pd, np, sp = _pd(), _np(), _sp()
    Z, znames = design_matrix(cov)
    Xg = sp.csc_matrix(guides.X)
    gidx = {str(g): i for i, g in enumerate(guides.var_names)}
    Xe = sp.csc_matrix(rna.X)
    eidx = {str(g): i for i, g in enumerate(rna.var_names)}
    fits: Dict[str, Tuple[Any, float]] = {}
    rows, skipped = [], []
    elements = list(dict.fromkeys(pairs["element"]))
    if max_elements:
        elements = elements[:max_elements]
    for k, e in enumerate(elements):
        genes = [g for g in pairs.loc[pairs["element"] == e, "gene_id"] if g in eidx]
        gtype = guide_type(e)
        for gname in el_guides.get(e, []):
            if gname not in gidx:
                continue
            x = _dense(Xg[:, gidx[gname]]).ravel() > 0
            n1 = int(x.sum())
            if n1 <= min_cells:
                skipped.append({"gRNA_id": gname, "element": e, "n_cells_guide": n1, "reason": f"<= {min_cells} cells"})
                continue
            Y = {}
            for g in genes:
                y = _dense(Xe[:, eidx[g]]).ravel().astype(float)
                if g not in fits:
                    fits[g] = nb_null_fit(y, Z)
                Y[g] = y
            if not Y:
                continue
            for r in crt_guide_genes(Y, x, Z, fits, B=B, side=side, seed=seed + k):
                rows.append({**r, "gRNA_id": gname, "pair_type": "all_elements_test", "n_cells_guide": n1, "element": e,
                             "guide_type": gtype, "engine": "python-crt"})
    res = pd.DataFrame(rows, columns=RESULT_COLS)
    return res, pd.DataFrame(skipped, columns=["gRNA_id", "element", "n_cells_guide", "reason"]), znames


SCEPTRE_R = r'''
setwd("{path}")
library(sceptre)
library(tibble)
library(dplyr)
exp <- as.matrix(read.table('gene_exp_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
g_I <- as.matrix(read.table('guides_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
cov <- read.table('covariates.txt', sep=',', header=TRUE, row.names=1)
if ('bath_number' %in% names(cov)) {{ cov$bath_number <- as.factor(cov$bath_number) }}
pairs_test <- as_tibble(read.table('pairs.txt', sep=',', header=TRUE, check.names=FALSE))
pairs_test$pair_type <- as.factor(pairs_test$pair_type)
result <- run_sceptre_high_moi(gene_matrix = exp, combined_perturbation_matrix = g_I, covariate_matrix = cov,
                               gene_gRNA_group_pairs = pairs_test, side = '{side}')
write.table(result, "results.txt")
'''


def write_sceptre_r_inputs(outdir: Path, rna, guides, cov, pairs, el_guides, side: str = "both", min_cells: int = SCEPTRE_MIN_CELLS,
                           execute: bool = False) -> List[Dict[str, Any]]:
    """runSceptre.py: per-element gene x cell, guide x cell, covariates, pairs, R script; Rscript if present."""
    pd, np = _pd(), _np()
    cov_r = cov.copy()
    if cov_r["bath_number"].nunique() == 1:
        del cov_r["bath_number"]
    cov_r = cov_r[[c for c in ("bath_number", "percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count") if c in cov_r]]
    runs = []
    rscript = which("Rscript")
    for e in dict.fromkeys(pairs["element"]):
        genes = [g for g in pairs.loc[pairs["element"] == e, "gene_id"] if g in set(rna.var_names)]
        gl = [g for g in el_guides.get(e, []) if g in set(guides.var_names)]
        if not genes or not gl:
            continue
        d = outdir / "sceptre_out" / safe_label(e)
        d.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(_dense(rna[:, genes].X).T, index=genes, columns=rna.obs_names).to_csv(d / "gene_exp_one_gene_guide.txt")
        G = pd.DataFrame(_dense(guides[:, gl].X).T, index=gl, columns=guides.obs_names)
        G.to_csv(d / "guides_one_gene_guide.txt")
        cov_r.to_csv(d / "covariates.txt")
        ok = [g for g in gl if G.loc[g].sum() > min_cells]
        pd.DataFrame([[g, gg, "all_elements_test"] for gg in ok for g in genes], columns=["gene_id", "gRNA_group", "pair_type"]).to_csv(d / "pairs.txt", index=False)
        (d / "r_script.r").write_text(SCEPTRE_R.format(path=str(d), side=side))
        cmd = ["Rscript", str(d / "r_script.r")]
        status = "written"
        if execute and rscript:
            rc = subprocess.run(cmd, capture_output=True, text=True)
            status = "ran" if rc.returncode == 0 else f"failed ({rc.stderr.strip()[-200:]})"
        runs.append({"element": e, "dir": str(d), "command": " ".join(cmd), "status": status, "n_guides_tested": len(ok)})
    if not rscript:
        print("Rscript not on PATH: SCEPTRE inputs written; run each element with `Rscript <dir>/r_script.r` (needs katsevich-lab/sceptre v0.1 API).")
    return runs


def read_sceptre_r_results(outdir: Path):
    pd = _pd()
    files = sorted((outdir / "sceptre_out").glob("**/results.txt"))
    if not files:
        return pd.DataFrame(columns=RESULT_COLS)
    df = pd.concat([pd.read_csv(f, sep=" ") for f in files], ignore_index=True)
    if "gRNA_group" in df and "gRNA_id" not in df:
        df = df.rename(columns={"gRNA_group": "gRNA_id"})
    df["engine"] = "R-sceptre"
    return df


# ---------------------------------------------------------------------------
# Results aggregation (sceptre_anndata_creation.py)
# ---------------------------------------------------------------------------

def bh(p):
    np = _np()
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    o = np.argsort(q)
    ranked = q[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.minimum(ranked, 1)
    out[ok] = res
    return out


def aggregate_results(res, guide_var=None):
    from scipy.stats import combine_pvalues  # type: ignore
    pd, np = _pd(), _np()
    res = res.copy()
    res["adj_pvalue"] = bh(res["p_value"].values)
    res["significant"] = res["adj_pvalue"] < 0.01
    if "element" not in res or res["element"].isna().any():
        res["element"] = [str(g).split("|")[0] for g in res["gRNA_id"]]
    rows = []
    for (e, g), v in res.groupby(["element", "gene_id"], sort=False):
        ps = [x for x in v["p_value"].values if np.isfinite(x)]
        p = float(combine_pvalues(ps).pvalue) if ps else float("nan")
        rows.append({"element": e, "gene_id": g, "p_value_fisher": p, "n_guides": len(ps),
                     "mean_log_fold_change": _finite_mean(v["log_fold_change"].values),
                     "guide_type": v["guide_type"].iloc[0] if "guide_type" in v else None})
    el = pd.DataFrame(rows, columns=["element", "gene_id", "p_value_fisher", "n_guides", "mean_log_fold_change", "guide_type"])
    el["sig_not_adj"] = el["p_value_fisher"] < 0.05
    el["adj_pvalue_fisher"] = bh(el["p_value_fisher"].values)
    if guide_var is not None and len(el):
        gv = guide_var.drop_duplicates("target_elements").set_index("target_elements")
        for src, dst in (("guide_chr", "element_chr"), ("guide_start", "element_start"), ("guide_end", "element_end")):
            if src in gv:
                el[dst] = [gv[src].get(e, "") for e in el["element"]]
    return res, el


def results_mudata(res, el, mods):
    ad, pd, np = _ad(), _pd(), _np()
    guides_ = sorted(set(res["gRNA_id"]))
    genes = sorted(set(res["gene_id"]))
    gi, ci = {g: i for i, g in enumerate(guides_)}, {g: i for i, g in enumerate(genes)}

    def mat(df, rkey, col, rix):
        M = np.full((len(rix), len(genes)), np.nan)
        for a, b, v in df[[rkey, "gene_id", col]].itertuples(index=False):
            M[rix[a], ci[b]] = v
        return M
    rvar = mods["scRNA"].var
    var = rvar.loc[[g for g in genes]] if set(genes) <= set(rvar.index) else pd.DataFrame(index=genes)
    gvar = mods["guides"].var
    obs = gvar.loc[guides_] if set(guides_) <= set(gvar.index) else pd.DataFrame(index=guides_)
    rg = ad.AnnData(X=mat(res, "gRNA_id", "p_value", gi), obs=obs.copy(), var=var.copy())
    for c in ("adj_pvalue", "z_value", "log_fold_change"):
        rg.layers[c] = mat(res, "gRNA_id", c, gi)
    rg.layers["significant"] = (np.nan_to_num(rg.layers["adj_pvalue"], nan=1) < 0.01).astype(np.int8)
    els = sorted(set(el["element"]))
    ei = {e: i for i, e in enumerate(els)}
    eobs = pd.DataFrame(index=els)
    for c in ("element_chr", "element_start", "element_end"):
        if c in el:
            eobs[c] = el.drop_duplicates("element").set_index("element")[c].reindex(els).astype(str).values
    re_ = ad.AnnData(X=mat(el, "element", "p_value_fisher", ei), obs=eobs, var=var.copy())
    re_.layers["sig_not_adj"] = (np.nan_to_num(re_.X, nan=1) < 0.05).astype(np.int8)
    return {"guides": mods["guides"], "scRNA": mods["scRNA"], "result_guides": rg, "result_elements": re_}


# ---------------------------------------------------------------------------
# Original Gasperini 2019 test and the comparison
# ---------------------------------------------------------------------------

def size_factors(tot):
    np = _np()
    tot = np.maximum(np.asarray(tot, float), 1)
    return tot / math.exp(np.mean(np.log(tot)))


def gasperini_test(rna, guides, pairs, el_guides, batch=None, percent_mito=None, min_cells: int = 1):
    """NB GLM LRT per gRNA group x gene with size-factor offset and guide_count / percent_mito / batch covariates."""
    from scipy.stats import chi2  # type: ignore
    pd, np, sp = _pd(), _np(), _sp()
    Xg = sp.csc_matrix(guides.X)
    gidx = {str(g): i for i, g in enumerate(guides.var_names)}
    guide_count = np.asarray((Xg > 0).sum(1)).ravel().astype(float)
    Xe = sp.csc_matrix(rna.X)
    eidx = {str(g): i for i, g in enumerate(rna.var_names)}
    tot = np.asarray(Xe.sum(1)).ravel()
    off = np.log(size_factors(tot))
    pm = np.asarray(percent_mito if percent_mito is not None else np.zeros(rna.n_obs), float)
    cols = [np.ones(rna.n_obs), (guide_count - guide_count.mean()) / (guide_count.std() or 1)]
    if pm.std() > 0:
        cols.append((pm - pm.mean()) / pm.std())
    if batch is not None:
        lev = sorted(set(map(str, batch)))
        for lv in lev[1:]:
            cols.append((np.asarray(list(map(str, batch))) == lv).astype(float))
    Z0 = np.column_stack(cols)
    null_fits: Dict[str, Tuple[Any, Any, float, float]] = {}
    rows = []
    for e in dict.fromkeys(pairs["element"]):
        gl = [g for g in el_guides.get(e, []) if g in gidx]
        if not gl:
            continue
        x = (_dense(Xg[:, [gidx[g] for g in gl]]).sum(1) > 0).astype(float).ravel()
        n1 = int(x.sum())
        if n1 < min_cells:
            continue
        for gname in pairs.loc[pairs["element"] == e, "gene_id"]:
            if gname not in eidx:
                continue
            y = _dense(Xe[:, eidx[gname]]).ravel().astype(float)
            if gname not in null_fits:
                b0, mu0 = nb_irls(y, Z0, off)
                th = theta_ml(y, mu0)
                b0, mu0 = nb_irls(y, Z0, off, theta=th, beta=b0)
                null_fits[gname] = (b0, mu0, th, nb_loglik(y, mu0, th))
            b0, mu0, th, ll0 = null_fits[gname]
            Z1 = np.column_stack([Z0, x])
            b1, mu1 = nb_irls(y, Z1, off, theta=th, beta=np.r_[b0, 0.0])
            ll1 = nb_loglik(y, mu1, th)
            stat = max(2 * (ll1 - ll0), 0.0)
            rows.append({"pairs4merge": f"{e}:{gname}", "beta": float(b1[-1]), "intercept": float(b1[0]),
                         "fold_change.transcript_remaining": float(math.exp(b1[-1])), "pvalue.raw": float(chi2.sf(stat, 1)),
                         "gene_short_name": gname, "gRNA_group": e, "site_type": guide_type(e) or "", "n_cells_group": n1, "theta": th})
    df = pd.DataFrame(rows)
    if len(df):
        ntc = df[df["gRNA_group"].map(lambda e: guide_type(e) == "NEGATIVE_CONTROL")]["pvalue.raw"].values
        if len(ntc):
            s = np.sort(ntc)
            df["pvalue.empirical"] = (1 + np.searchsorted(s, df["pvalue.raw"].values, side="right")) / (1 + len(s))
            df["pvalue.empirical.adjusted"] = bh(df["pvalue.empirical"].values)
        else:
            df["pvalue.empirical"] = "not_applicable"
            df["pvalue.empirical.adjusted"] = "not_applicable"
    return df


def compare_to_original(ours, orig, our_p: str, our_effect: str, our_element: str = "element", our_gene: str = "gene_id",
                        group_suffix: str = "_TSS", fdr: float = 0.1):
    from scipy.stats import fisher_exact, spearmanr  # type: ignore
    pd, np = _pd(), _np()
    o = orig.copy()
    o["pvalue.raw"] = pd.to_numeric(o["pvalue.raw"], errors="coerce")
    o["beta"] = pd.to_numeric(o["beta"], errors="coerce")
    a = ours.copy()
    a["_group"] = [e if str(e).endswith(group_suffix) or guide_type(str(e)) == "NEGATIVE_CONTROL" else f"{e}{group_suffix}" for e in a[our_element]]
    m = a.merge(o, left_on=["_group", our_gene], right_on=["gRNA_group", "gene_short_name"], how="inner", suffixes=("", "_orig"))
    m = m[np.isfinite(pd.to_numeric(m[our_p], errors="coerce")) & np.isfinite(m["pvalue.raw"])]
    out: Dict[str, Any] = {"n_ours": int(len(a)), "n_original": int(len(o)), "n_matched": int(len(m))}
    if len(m) < 3:
        return m, out
    p1 = np.maximum(m[our_p].astype(float).values, 1e-300)
    p2 = np.maximum(m["pvalue.raw"].values, 1e-300)
    rho, pr = spearmanr(-np.log10(p1), -np.log10(p2))
    out.update({"spearman_rho_log10p": float(rho), "spearman_p": float(pr)})
    q1, q2 = bh(p1), bh(p2)
    h1, h2 = q1 < fdr, q2 < fdr
    tab = [[int((h1 & h2).sum()), int((h1 & ~h2).sum())], [int((~h1 & h2).sum()), int((~h1 & ~h2).sum())]]
    orr, fp = fisher_exact(tab) if min(len(h1), 1) else (float("nan"), float("nan"))
    out.update({"fdr": fdr, "hits_ours": int(h1.sum()), "hits_original": int(h2.sum()), "hits_both": tab[0][0],
                "contingency": tab, "fisher_odds_ratio": float(orr), "fisher_p": float(fp)})
    eff = m[our_effect].astype(float).values
    both = h1 & h2 & np.isfinite(eff)
    out["sign_agreement_hits"] = float(np.mean(np.sign(eff[both]) == np.sign(m["beta"].values[both]))) if both.any() else None
    fin = np.isfinite(eff)
    out["sign_agreement_all"] = float(np.mean(np.sign(eff[fin]) == np.sign(m["beta"].values[fin]))) if fin.any() else None
    self_tss = m[m["site_type"].astype(str).isin(["selfTSS"])] if "site_type" in m else m.iloc[:0]
    if len(self_tss):
        out["self_tss_pairs"] = int(len(self_tss))
        out["self_tss_hits_ours"] = int((bh(np.maximum(self_tss[our_p].astype(float).values, 1e-300)) < fdr).sum())
    m = m.assign(q_ours=q1, q_original=q2, hit_ours=h1, hit_original=h2)
    return m, out


# ---------------------------------------------------------------------------
# MULTI-seq (deMULTIplex port)
# ---------------------------------------------------------------------------

def multiseq_preprocess(r1: str, r2: str, cell_ids: Sequence[str], cell=(1, 16), umi=(17, 28), tag=(1, 8)):
    """MULTIseq.preProcess: 1-based inclusive ranges on R1 (cell, UMI) and R2 (tag); reads of listed cells only."""
    pd = _pd()
    keep = set(cell_ids)
    opener = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)
    rows = []
    with opener(r1) as f1, opener(r2) as f2:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            s1, s2 = f1.readline().strip(), f2.readline().strip()
            f1.readline(); f1.readline(); f2.readline(); f2.readline()
            c = s1[cell[0] - 1:cell[1]]
            if c in keep:
                rows.append((c, s1[umi[0] - 1:umi[1]], s2[tag[0] - 1:tag[1]]))
    return pd.DataFrame(rows, columns=["Cell", "UMI", "Sample"])


def multiseq_align(reads, cell_ids: Sequence[str], bar_ref: Sequence[str], max_mismatch: int = 1):
    """MULTIseq.align: tag -> reference barcode with Hamming <= 1; unique UMIs per cell x barcode + nUMI, nUMI_total."""
    pd, np = _pd(), _np()
    ref = list(bar_ref)

    def match(t):
        best, bd = None, max_mismatch + 1
        for i, r in enumerate(ref):
            d = sum(1 for a, b in zip(t, r) if a != b) + abs(len(t) - len(r))
            if d < bd:
                best, bd = i, d
        return best if bd <= max_mismatch else None
    cache: Dict[str, Optional[int]] = {}
    tbl = pd.DataFrame(0, index=list(cell_ids), columns=[f"Bar{i + 1}" for i in range(len(ref))])
    total = reads.groupby("Cell")["UMI"].nunique() if len(reads) else pd.Series(dtype=float)
    m = []
    for s in reads["Sample"].values:
        if s not in cache:
            cache[s] = match(s)
        m.append(cache[s])
    r = reads.assign(bar=m).dropna(subset=["bar"])
    if len(r):
        cnt = r.groupby(["Cell", "bar"])["UMI"].nunique()
        for (c, b), v in cnt.items():
            tbl.loc[c, f"Bar{int(b) + 1}"] = v
    tbl["nUMI"] = tbl[[f"Bar{i + 1}" for i in range(len(ref))]].sum(1)
    tbl["nUMI_total"] = total.reindex(tbl.index).fillna(0).astype(int).values
    return tbl


def _bw_nrd0(x) -> float:
    np = _np()
    x = np.asarray(x, float)
    hi = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(hi, iqr / 1.34) if iqr > 0 else hi
    if not lo:
        lo = hi or abs(x[0]) or 1.0
    return 0.9 * lo * len(x) ** (-0.2)


def _local_maxima(v) -> List[int]:
    """deMULTIplex localMaxima: plateau-aware indices of peaks."""
    np = _np()
    y = np.diff(np.r_[-np.inf, v]) > 0
    runs, i = [], 0
    while i < len(y):
        j = i
        while j < len(y) and y[j] == y[i]:
            j += 1
        runs.append(j - i)
        i = j
    cs = np.cumsum(runs)
    idx = list(cs[::2] - 1) if y[0] else list(cs[1::2] - 1)
    if len(v) > 1 and v[0] == v[1] and idx and idx[0] == 0:
        idx = idx[1:]
    return [int(k) for k in idx if 0 <= k < len(v)]


def classify_cells(bar_table, q: float):
    """deMULTIplex classifyCells: per-barcode KDE thresholds at quantile q between the low and high modes."""
    pd, np = _pd(), _np()
    with np.errstate(divide="ignore"):
        X = np.log2(bar_table.values.astype(float))
    X[~np.isfinite(X)] = 0
    X = X - X.mean(0)
    calls = np.array(["Negative"] * X.shape[0], dtype=object)
    for i, name in enumerate(bar_table.columns):
        v = X[:, i]
        if np.std(v) == 0:
            continue
        bw = _bw_nrd0(v)
        grid = np.linspace(np.quantile(v, 0.001), np.quantile(v, 0.999), 100)
        dens = np.exp(-0.5 * ((grid[:, None] - v[None, :]) / bw) ** 2).sum(1)
        ext = _local_maxima(dens)
        if len(ext) <= 1:
            continue
        low = ext[int(np.argmax(dens[ext]))]
        high = max(ext)
        if low == high:
            continue
        lo_s, hi_s = grid[low], grid[high]
        a_, b_ = min(lo_s, hi_s), max(lo_s, hi_s)
        thr = a_ + q * (b_ - a_)  # R quantile(type 7) of two points
        for c in np.where(v >= thr)[0]:
            calls[c] = name if calls[c] == "Negative" else "Doublet"
    return pd.Series(calls, index=bar_table.index)


def find_thresh(bar_table) -> Tuple[float, Any]:
    pd, np = _pd(), _np()
    res = []
    for q in np.round(np.arange(0.01, 0.991, 0.02), 2):
        c = classify_cells(bar_table, q)
        n = len(c)
        res.append({"q": q, "Negative": (c == "Negative").sum() / n, "Doublet": (c == "Doublet").sum() / n,
                    "Singlet": (~c.isin(["Negative", "Doublet"])).sum() / n})
    res = pd.DataFrame(res)
    return float(res.loc[res["Singlet"].idxmax(), "q"]), res


def multiseq_classify(bar_table_full, upstream_compat: bool = False):
    """Two rounds of quantile sweeps; round-1 negatives removed; final calls for every cell."""
    pd = _pd()
    drop = ["nUMI"] if upstream_compat else ["nUMI", "nUMI_total"]
    bt = bar_table_full.drop(columns=[c for c in drop if c in bar_table_full.columns])
    q1, sweep1 = find_thresh(bt)
    r1 = classify_cells(bt, q1)
    neg = list(r1.index[r1 == "Negative"])
    bt2 = bt.drop(index=neg)
    q2, sweep2 = find_thresh(bt2)
    r2 = classify_cells(bt2, q2)
    neg += list(r2.index[r2 == "Negative"])
    final = pd.concat([r2[r2 != "Negative"], pd.Series("Negative", index=neg)])
    return final.reindex(bar_table_full.index).fillna("Negative"), {"q_round1": q1, "q_round2": q2, "sweep1": sweep1, "sweep2": sweep2}


# ---------------------------------------------------------------------------
# inspect
# ---------------------------------------------------------------------------

def inspect_mods(mods) -> List[Dict[str, Any]]:
    np = _np()
    rows = []

    def add(ok, what, detail=""):
        rows.append({"status": "PASS" if ok else "FAIL", "check": what, "detail": detail})
    for mod, spec in PILOT_SCHEMA.items():
        if mod not in mods:
            add(False, f"modality {mod}", f"present: {list(mods)}")
            continue
        a = mods[mod]
        add(True, f"modality {mod}", f"{a.n_obs} x {a.n_vars}")
        for part in ("obs", "var"):
            cols = set(getattr(a, part).columns)
            miss = [c for c in spec[part] if c not in cols]
            add(not miss, f"{mod}.{part} columns", f"missing {miss}" if miss else "all present")
        for k in spec["uns"]:
            add(k in a.uns, f"{mod}.uns['{k}']")
        for k in spec["varm"]:
            add(k in a.varm, f"{mod}.varm['{k}']")
    if "guides" in mods and "scRNA" in mods:
        g, r = mods["guides"], mods["scRNA"]
        add(list(g.obs_names) == list(r.obs_names), "guides and scRNA share barcodes in order", f"{len(set(g.obs_names) & set(r.obs_names))} shared")
        X = _dense(g.X[: min(g.n_obs, 2000)])
        add(bool(np.all(X >= 0) and np.allclose(X, np.round(X))), "guide X holds non-negative integer counts")
        if "target_elements" in g.var and "guide_by_element" in g.varm and "elements" in g.uns:
            M = _dense(g.varm["guide_by_element"])
            els = list(map(str, g.uns["elements"]))
            ok = M.shape == (g.n_vars, len(els)) and all(els[int(np.argmax(M[i]))] == str(t) for i, t in enumerate(g.var["target_elements"]))
            add(ok, "guide_by_element one-hot matches target_elements", f"{len(els)} elements")
        if "feature_name" in g.var:
            bad = [n for n in g.var["feature_name"].astype(str) if "|" not in n or "_sgrna_" not in n]
            add(not bad, "guide names follow '<element>|<n>_sgrna_<chr>:<start>:<end>'", f"{len(bad)} nonconforming")
        if "percent_mito" in r.obs:
            pm = r.obs["percent_mito"].astype(float)
            add(bool((pm >= 0).all() and (pm <= 1).all()), "percent_mito is a fraction in [0, 1]")
    return rows


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:  # noqa: S310 (pinned public URLs)
        shutil.copyfileobj(resp, fh, 1 << 20)
    tmp.replace(dest)
    return dest


def cmd_setup(args: argparse.Namespace) -> int:
    dest = Path(args.dest) if args.dest else DATA_ROOT
    todo = dict(SETUP_FILES)
    if args.no_geo:
        todo = {k: v for k, v in todo.items() if not k.startswith("GSE")}
    if args.whitelist:
        todo[WHITELIST[0]] = (WHITELIST[1], WHITELIST[2], "10x v2 barcode whitelist")
    total = sum(v[1] for v in todo.values())
    print(f"setup: {len(todo)} file(s), ~{total / 1e6:.0f} MB into {dest}")
    bad = 0
    for name, (url, size, what) in todo.items():
        p = dest / name
        if p.exists() and not args.force:
            print(f"  present {p} ({what})")
            continue
        try:
            _download(url, p)
            print(f"Wrote: {p}")
        except Exception as e:
            bad += 1
            print(f"  failed {name}: {e}")
    if args.print_raw:
        print("\n# Raw pilot data (README download script; several GB each -- not downloaded):")
        print(f"wget {BAMTOFASTQ}; chmod +x bamtofastq_linux")
        for name, url, outd in RAW_DATA:
            print(f"wget -O {name} {url}\n./bamtofastq_linux --nthreads=20 {name} {outd}")
        print(f"# references: {GTF_URL}\n#             {GENOME_URL}")
    return 0 if bad == 0 else 1


def cmd_guide_table(args: argparse.Namespace) -> int:
    df = guide_features(read_table(args.table))
    out = run_dir(args.label)
    feat = out / "guide_features.txt"
    df[["sgRNA_sequences", "pipeline_id"]].to_csv(feat, sep="\t", header=False, index=False)
    print(f"Wrote: {feat}")
    write_tsv(df, out / "guide_table_processed.tsv")
    n_el = df["element"].nunique()
    summ = {"n_guides": int(len(df)), "n_elements": int(n_el), "guides_per_element": df.groupby("element").size().value_counts().to_dict(),
            "guide_types": df["element"].map(lambda e: guide_type(e) or "unlabelled").value_counts().to_dict()}
    write_json(out / "summary.json", summ)
    write_report(out / "report.md", [f"# Guide table: {Path(args.table).name}", "", f"{len(df)} guides, {n_el} elements.", "",
                                     "Pipeline ids (`<Target>|<n>_sgrna_<chr>:<start>:<end>`) -- first rows:", ""] +
                 [f"- `{p}` {s}" for p, s in df[["pipeline_id", "sgRNA_sequences"]].head(10).itertuples(index=False)])
    print(f"{len(df)} guides for {n_el} elements -> {feat}")
    return 0


def cmd_composition(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    summ = {}
    for fq in args.fastq:
        df = fastq_composition(fq, args.n_lines, args.upstream_compat)
        tag = safe_label(Path(fq).name)
        write_tsv(df, out / f"{tag}_composition.tsv")
        summ[Path(fq).name] = {"n_reads": int(df["n_reads"].iloc[0]) if len(df) else 0, "length": int(len(df)),
                               "max_sd_position": int(df["sd"].idxmax()) if len(df) else None,
                               "fixed_positions_sd_gt_0.4n": int((df["sd"] > 0.4 * df["n_reads"]).sum()) if len(df) else 0}
        plt = None if args.no_plots else _plt()
        if plt is not None and len(df):
            fig, ax = plt.subplots(figsize=(10, 2.2))
            ax.plot(df["position"], df["sd"], color=BLUE, lw=1.2)
            ax.set_xlim(0, len(df))
            ax.set_title(f"{Path(fq).name}: compositional bias (sd between the 4 nucleotides)", loc="left", fontsize=9)
            ax.set_xlabel("read position")
            ax.set_ylabel("sd of counts")
            save_fig(fig, out / f"{tag}_composition.png")
            plt.close(fig)
    write_json(out / "summary.json", summ)
    write_report(out / "report.md", ["# Read composition", ""] + [f"- {k}: {v['n_reads']} reads x {v['length']} bp; "
                                                                  f"{v['fixed_positions_sd_gt_0.4n']} near-fixed positions" for k, v in summ.items()])
    return 0


def kb_commands(fq_rna: Sequence[str], fq_guide: Sequence[str], features: str, whitelist: str, chemistry: str, threads: int,
                name: str = "S1_L1") -> List[List[str]]:
    return [
        ["kb", "ref", "-d", "human", "-i", "transcriptome_index.idx", "-g", "transcriptome_t2g.txt"],
        ["kb", "ref", "-i", "guide_index.idx", "-f1", "genome.fa.gz", "-g", "t2guide.txt", "--workflow", "kite", features],
        ["kb", "count", "-i", "transcriptome_index.idx", "-g", "transcriptome_t2g.txt", "--verbose", "--workflow", "kite", "-w", whitelist,
         "--h5ad", "-x", chemistry, "-o", f"{name}_ks_transcripts_out", "-t", str(threads), *fq_rna, "--overwrite", "-m", "48G"],
        ["kb", "count", "-i", "guide_index.idx", "-g", "t2guide.txt", "--verbose", "--report", "--workflow", "kite", "-w", whitelist,
         "--h5ad", "-x", chemistry, "-o", f"{name}_ks_guide_out", "-t", str(threads), *fq_guide, "--overwrite", "-m", "48G"],
    ]


def cmd_kb(args: argparse.Namespace) -> int:
    cmds = kb_commands(args.rna_fastqs or ["R1.fastq.gz", "R2.fastq.gz"], args.guide_fastqs or ["gR1.fastq.gz", "gR2.fastq.gz"],
                       args.guide_features, args.whitelist, args.chemistry, args.threads, args.name)
    have = which("kb")
    for c in cmds:
        line = " ".join(c)
        if args.execute and have:
            print(f"$ {line}")
            rc = subprocess.run(c).returncode
            if rc:
                return rc
        else:
            print(line)
    if not have:
        print("# kb (kb-python + kallisto) not on PATH: commands printed only")
    return 0


def cmd_config(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    p = dict(CONFIG)
    lines = [f"params.GTF_GZ_LINK = '{GTF_URL}'", f"params.TRANSCRIPTOME_REFERENCE = \"{p['TRANSCRIPTOME_REFERENCE']}\"",
             f"params.KALLISTO_BIN = '{args.kallisto_bin or which('kallisto') or 'ADD YOUR KALLISTO BIN'}'", f"params.GENOME = '{GENOME_URL}'",
             f"params.GUIDE_FEATURES = '{args.guide_features}'", f"params.CHEMISTRY = \"{args.chemistry}\"", f"params.THREADS = {p['THREADS']}",
             f"params.DISTANCE_NEIGHBORS = {p['DISTANCE_NEIGHBORS']}", f"params.IN_TRANS = \"{p['IN_TRANS']}\"",
             f"params.FASTQ_FILES_TRANSCRIPTS = ['{' '.join(args.rna_fastqs or [])}']", "params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']",
             "params.CUSTOM_REFERENCE = false", "params.CUSTOM_REFERENCE_IDX = ''", "params.CUSTOM_REFERENCE_T2T = ''", "params.CUSTOM_GTF_PATH = ''",
             f"params.FASTQ_FILES_GUIDES = ['{' '.join(args.guide_fastqs or [])}']", "params.FASTQ_NAMES_GUIDES = ['S1_L1']",
             "params.CREATE_REF = false", "params.ADDGENENAMES = ''", f"params.DIRECTION = '{p['DIRECTION']}'", f"params.WHITELIST= '{args.whitelist}'",
             f"params.EXPECTED_CELL_NUMBER = {p['EXPECTED_CELL_NUMBER']}", f"params.MITO_SPECIE = '{p['MITO_SPECIE']}'",
             f"params.MITO_EXPECTED_PERCENTAGE = {p['MITO_EXPECTED_PERCENTAGE']}",
             f"params.PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS={p['PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS']}",
             f"params.TRANSCRIPTS_UMI_TRHESHOLD = {p['TRANSCRIPTS_UMI_TRHESHOLD']}", f"params.GUIDE_UMI_LIMIT = {p['GUIDE_UMI_LIMIT']}",
             "params.MERGE = false", "params.RUN_MULTISEQ = false", "params.R1_MULTI = 'NOT APPLICABLE'", "params.R2_MULTI = 'NOT APPLICABLE'",
             "params.BARCODES_MULTIBAR_LIST_MULTI = 'NOT APPLICABLE'", "params.BAR_MULTI= [1,16]", "params.UMI_MULTI= [17,28]", "params.R2_MULTI_TAG = [1,8]"]
    cfg = out / "gasperini_sample.config"
    cfg.write_text("\n".join(lines) + "\n")
    print(f"Wrote: {cfg}")
    cmd = f"export NXF_VER=22.10.6; nextflow run {args.main_nf} -c {cfg} -w {args.workdir}"
    print(cmd if which("nextflow") else f"# nextflow not on PATH; to run: {cmd}")
    write_json(out / "summary.json", {"config": str(cfg), "command": cmd, "params": p})
    write_report(out / "report.md", ["# Nextflow config (Gasperini sample)", "", "```", *lines, "```", "", f"`{cmd}`"])
    return 0


def cmd_qc_filter(args: argparse.Namespace) -> int:
    ad = _ad()
    if len(args.rna) != len(args.guide):
        raise SystemExit("--rna and --guide need the same number of lanes")
    lanes = args.lanes or [str(i + 1) for i in range(len(args.rna))]
    out = run_dir(args.label)
    rnas, gds, stats = [], [], []
    for i, (r, g) in enumerate(zip(args.rna, args.guide)):
        a, b, st = analyze_batch(read_counts(r), read_counts(g), int(lanes[i]) if str(lanes[i]).isdigit() else lanes[i],
                                 args.expected_cell_number, args.mito_expected_percentage, args.transcripts_umi_threshold,
                                 args.upstream_compat, not args.no_doublets, seed=args.seed)
        rnas.append(a)
        gds.append(b)
        stats.append(st)
        print(f"lane {lanes[i]}: {st['n_barcodes']} barcodes -> {st['n_after_umi_knee']} (UMI/knee) -> {st['n_after_mito']} (mito) "
              f"-> {st['n_final']} (doublets {st['n_doublets']}, shared barcodes)")
    r, g = concat_lanes(rnas, gds, args.percentage_of_cells_to_include_transcript, args.upstream_compat)
    d = out / "results_per_lane"
    d.mkdir(exist_ok=True)
    r.write_h5ad(d / "full_raw_scrna_ann_data.h5ad")
    g.write_h5ad(d / "full_raw_guide_ann_data.h5ad")
    print(f"Wrote: {d / 'full_raw_scrna_ann_data.h5ad'}\nWrote: {d / 'full_raw_guide_ann_data.h5ad'}")
    pd = _pd()
    write_tsv(pd.DataFrame(stats), out / "lane_qc.tsv")
    summ = {"lanes": stats, "n_cells": int(r.n_obs), "n_genes": int(r.n_vars), "n_guides": int(g.n_vars), "upstream_compat": args.upstream_compat}
    write_json(out / "summary.json", summ)
    write_report(out / "report.md", ["# Cell QC and lane merge", "", f"{r.n_obs} cells x {r.n_vars} genes (genes in >= {args.percentage_of_cells_to_include_transcript:.0%} of cells), {g.n_vars} guides.", "",
                                     "| lane | barcodes | after UMI + knee | after mito | doublets | final |", "|---|---|---|---|---|---|"] +
                 [f"| {s['lane']} | {s['n_barcodes']} | {s['n_after_umi_knee']} | {s['n_after_mito']} | {s['n_doublets']} | {s['n_final']} |" for s in stats])
    return 0


def cmd_mudata(args: argparse.Namespace) -> int:
    ad = _ad()
    g = read_counts(args.guide_h5ad)
    r = read_counts(args.rna_h5ad)
    genes = read_gtf_genes(args.gtf) if args.gtf else None
    seqs = None
    if args.guide_features:
        ft = _pd().read_csv(args.guide_features, sep="\t", header=None, names=["seq", "id"])
        seqs = dict(zip(ft["id"], ft["seq"]))
    g = annotate_guides(g, seqs)
    r = annotate_transcripts(r, genes, args.upstream_compat)
    out = run_dir(args.label)
    p = write_h5mu(out / "raw_mudata_guide_and_transcripts.h5mu", {"guides": g, "scRNA": r})
    nf = int((r.var["transcript_chr"] == "NOT_FOUND").sum())
    write_json(out / "summary.json", {"h5mu": str(p), "n_cells": int(r.n_obs), "n_guides": int(g.n_vars), "n_genes": int(r.n_vars),
                                      "genes_without_coordinates": nf, "n_elements": int(len(g.uns["elements"]))})
    write_report(out / "report.md", ["# Raw MuData", "", f"guides {g.n_obs} x {g.n_vars}, scRNA {r.n_obs} x {r.n_vars}; {nf} genes without GTF coordinates.", "", f"`{p}`"])
    return 0


def _load_mods(path):
    mods = read_h5mu(path)
    if "guides" not in mods or "scRNA" not in mods:
        raise SystemExit(f"{path}: needs 'guides' and 'scRNA' modalities (has {list(mods)})")
    return mods


def cmd_assign(args: argparse.Namespace) -> int:
    np = _np()
    mods = _load_mods(args.h5mu)
    proc, cov = process_mudata(mods, args.guide_umi_limit, args.merge, args.upstream_compat)
    out = run_dir(args.label)
    p = write_h5mu(out / "processed_mudata_guide_and_transcripts.h5mu", proc)
    write_tsv(cov, out / "covariates.tsv", index=True)
    B = proc["guides"].X
    per_guide = np.asarray(B.sum(0)).ravel()
    per_cell = np.asarray(B.sum(1)).ravel()
    summ = {"h5mu": str(p), "guide_umi_limit": args.guide_umi_limit, "merge": args.merge, "n_cells": int(B.shape[0]), "n_guides": int(B.shape[1]),
            "mean_guides_per_cell": float(per_cell.mean()), "median_cells_per_guide": float(np.median(per_guide)),
            "guides_over_30_cells": int((per_guide > SCEPTRE_MIN_CELLS).sum())}
    write_json(out / "summary.json", summ)
    plt = None if args.no_plots else _plt()
    if plt is not None:
        fig, ax = plt.subplots(1, 2, figsize=(8, 2.8))
        ax[0].hist(per_cell, bins=range(int(per_cell.max()) + 2), color=BLUE)
        ax[0].set_title("guides per cell", loc="left", fontsize=9)
        ax[1].hist(per_guide, bins=30, color=ORANGE)
        ax[1].set_title("cells per guide", loc="left", fontsize=9)
        save_fig(fig, out / "guide_assignment.png")
        plt.close(fig)
    write_report(out / "report.md", ["# Guide assignment", "", f"UMI > {args.guide_umi_limit}: {summ['mean_guides_per_cell']:.2f} guides per cell, "
                                     f"median {summ['median_cells_per_guide']:.0f} cells per guide; {summ['guides_over_30_cells']} guides with > 30 cells.", "", f"`{p}`"])
    return 0


def cmd_pairs(args: argparse.Namespace) -> int:
    mods = _load_mods(args.h5mu)
    rna = mods["scRNA"]
    rv = rna.var
    if "feature_name" in rv and not set(rv.index.astype(str)) & set(rv["feature_name"].astype(str)):
        rv = rv.set_index(rv["feature_name"].astype(str))
        rv = rv[~rv.index.duplicated()]
    genes = gene_table_from_rna(rv, read_gtf_genes(args.gtf) if args.gtf else None)
    pairs, el_g, how = tested_pairs(list(mods["guides"].var_names), list(rv.index.astype(str)), genes, args.distance, args.in_trans,
                                    args.add_genes or [], args.upstream_compat, args.seed)
    out = run_dir(args.label)
    write_tsv(pairs, out / "tested_pairs.tsv")
    pd = _pd()
    eg = pd.DataFrame([{"element": e, "gRNA_id": g, "guide_type": guide_type(e), "selection": how[e]} for e, gl in el_g.items() for g in gl])
    write_tsv(eg, out / "element_guides.tsv")
    write_json(out / "summary.json", {"n_elements": len(el_g), "n_pairs": int(len(pairs)), "selection": pd.Series(how).value_counts().to_dict(),
                                      "distance": args.distance, "in_trans": args.in_trans})
    write_report(out / "report.md", ["# Tested pairs", "", f"{len(el_g)} elements, {len(pairs)} element-gene pairs (window {args.distance:,} bp)."])
    return 0


def _run_tests(mods, pairs, el_g, cov, args, out: Path):
    if args.engine == "r":
        runs = write_sceptre_r_inputs(out, mods["scRNA"], mods["guides"], cov, pairs, el_g, args.direction, execute=True)
        _pd().DataFrame(runs).to_csv(out / "sceptre_r_runs.tsv", sep="\t", index=False)
        print(f"TSV: {out / 'sceptre_r_runs.tsv'}")
        res = read_sceptre_r_results(out)
        return res, _pd().DataFrame(), []
    return run_sceptre_python(mods["scRNA"], mods["guides"], cov, pairs, el_g, B=args.resamples, side=args.direction, seed=args.seed,
                              max_elements=getattr(args, "max_elements", None))


def _cov_from(mods):
    g = mods["guides"]
    cols = [c for c in ("bath_number", "percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count") if c in g.obs]
    if len(cols) >= 4:
        return g.obs[cols].copy()
    np = _np()
    return make_covariates(mods["scRNA"].obs, np.asarray((mods["guides"].X > 0).sum(1)).ravel())


def cmd_sceptre(args: argparse.Namespace) -> int:
    pd = _pd()
    mods = _load_mods(args.h5mu)
    pairs = pd.read_csv(args.pairs, sep="\t")
    el_g: Dict[str, List[str]] = {}
    for g in mods["guides"].var_names:
        el_g.setdefault(str(g).split("|")[0], []).append(str(g))
    out = run_dir(args.label)
    res, skipped, znames = _run_tests(mods, pairs, el_g, _cov_from(mods), args, out)
    write_tsv(res, out / "results.tsv")
    if len(skipped):
        write_tsv(skipped, out / "skipped_guides.tsv")
    write_json(out / "summary.json", {"engine": args.engine, "n_tests": int(len(res)), "n_guides_skipped": int(len(skipped)),
                                      "covariates": znames, "direction": args.direction, "resamples": args.resamples})
    write_report(out / "report.md", ["# SCEPTRE cis tests", "", f"{len(res)} guide-gene tests ({args.engine}); {len(skipped)} guides with <= 30 cells skipped."])
    return 0


def cmd_results(args: argparse.Namespace) -> int:
    pd = _pd()
    mods = _load_mods(args.h5mu)
    res = pd.read_csv(args.results, sep="\t")
    res, el = aggregate_results(res, mods["guides"].var)
    out = run_dir(args.label)
    write_tsv(res, out / "guide_gene_results.tsv")
    write_tsv(el, out / "element_gene_results.tsv")
    write_h5mu(out / "mudata_results.h5mu", results_mudata(res, el, mods))
    write_json(out / "summary.json", {"n_tests": int(len(res)), "n_significant_fdr01": int(res["significant"].sum()),
                                      "n_element_gene": int(len(el)), "n_element_sig_not_adj": int(el["sig_not_adj"].sum())})
    write_report(out / "report.md", ["# Results", "", f"{int(res['significant'].sum())} of {len(res)} guide-gene tests at BH < 0.01; "
                                     f"{int(el['sig_not_adj'].sum())} of {len(el)} element-gene pairs with Fisher p < 0.05."])
    return 0


def cmd_gasperini_test(args: argparse.Namespace) -> int:
    pd = _pd()
    mods = _load_mods(args.h5mu)
    pairs = pd.read_csv(args.pairs, sep="\t")
    el_g: Dict[str, List[str]] = {}
    for g in mods["guides"].var_names:
        el_g.setdefault(str(g).split("|")[0], []).append(str(g))
    r = mods["scRNA"]
    batch = r.obs["bath_number"] if "bath_number" in r.obs else r.obs.get("batch_number")
    df = gasperini_test(r, mods["guides"], pairs, el_g, batch, r.obs.get("percent_mito"))
    out = run_dir(args.label)
    write_tsv(df, out / "gasperini_nb_results.tsv")
    write_json(out / "summary.json", {"n_tests": int(len(df)), "n_p_lt_1e-3": int((df["pvalue.raw"] < 1e-3).sum()) if len(df) else 0})
    write_report(out / "report.md", ["# Gasperini 2019 NB likelihood-ratio test", "", f"{len(df)} gRNA-group x gene tests."])
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    pd = _pd()
    ours = pd.read_csv(args.results, sep="\t")
    orig = pd.read_csv(args.original, sep="\t")
    if "p_value_fisher" in ours:
        p, eff, el, gene = "p_value_fisher", "mean_log_fold_change", "element", "gene_id"
    elif "pvalue.raw" in ours:
        ours = ours.rename(columns={"pvalue.raw": "p_ours", "beta": "beta_ours"})
        p, eff, el, gene = "p_ours", "beta_ours", "gRNA_group", "gene_short_name"
    else:
        p, eff, el, gene = "p_value", "log_fold_change", "element", "gene_id"
    m, summ = compare_to_original(ours, orig, p, eff, el, gene, args.group_suffix, args.fdr)
    out = run_dir(args.label)
    write_tsv(m, out / "matched_pairs.tsv")
    write_json(out / "summary.json", summ)
    write_report(out / "report.md", _compare_lines(summ))
    return 0


def _compare_lines(s: Dict[str, Any]) -> List[str]:
    L = ["# Comparison with the original Gasperini 2019 results", "", f"Matched pairs: {s.get('n_matched', 0)} (ours {s.get('n_ours')}, original {s.get('n_original')})."]
    if "spearman_rho_log10p" in s:
        L += [f"Spearman rho of -log10 p: {_fmt(s['spearman_rho_log10p'])}; hits at BH < {s['fdr']}: ours {s['hits_ours']}, original {s['hits_original']}, "
              f"both {s['hits_both']} (Fisher OR {_fmt(s['fisher_odds_ratio'])}, p {_fmt(s['fisher_p'])}); effect sign agreement on shared hits "
              f"{_fmt(s['sign_agreement_hits'])}, on all pairs {_fmt(s['sign_agreement_all'])}."]
    return L


def cmd_multiseq(args: argparse.Namespace) -> int:
    pd = _pd()
    cells = [x.strip() for x in open(args.cell_barcodes) if x.strip()]
    ref = pd.read_csv(args.bar_ref, header=None).iloc[:, 0].astype(str).tolist()
    reads = multiseq_preprocess(args.r1, args.r2, cells, tuple(args.bar_multi), tuple(args.umi_multi), tuple(args.r2_multi_tag))
    tbl = multiseq_align(reads, cells, ref)
    calls, info = multiseq_classify(tbl, args.upstream_compat)
    out = run_dir(args.label)
    write_tsv(tbl, out / "bar_table.tsv", index=True)
    write_tsv(pd.DataFrame({"cell_barcode": calls.index, "multiseq_class": calls.values}), out / "final_class.tsv")
    write_tsv(info["sweep1"], out / "threshold_sweep_round1.tsv")
    vc = calls.value_counts().to_dict()
    if args.h5mu:
        mods = _load_mods(args.h5mu)
        keep = [c for c in mods["scRNA"].obs_names if calls.get(c, "Negative") not in ("Doublet", "Negative")]
        tb = tbl.loc[[c for c in keep if c in tbl.index]]
        ad = _ad()
        mm = ad.AnnData(X=tb[[c for c in tb.columns if c.startswith("Bar")]].values.astype(float), obs=pd.DataFrame({"multiseq_class": calls.reindex(tb.index).values}, index=tb.index))
        newmods = {"scRNA": mods["scRNA"][tb.index].copy(), "guides": mods["guides"][tb.index].copy(), "multiseq": mm}
        for k in ("scRNA", "guides"):
            newmods[k].obs["multiseq_class"] = calls.reindex(tb.index).values
        write_h5mu(out / "processed_mudata_guide_and_transcripts_multiseq_filtered.h5mu", newmods)
    write_json(out / "summary.json", {"n_reads": int(len(reads)), "n_cells": len(cells), "calls": vc, "q_round1": info["q_round1"], "q_round2": info["q_round2"]})
    write_report(out / "report.md", ["# MULTI-seq demultiplexing", "", f"{len(reads):,} reads from listed cells; q round 1 = {info['q_round1']}, round 2 = {info['q_round2']}.", ""] +
                 [f"- {k}: {v}" for k, v in sorted(vc.items(), key=lambda kv: -kv[1])])
    return 0


def cmd_inspect(args: argparse.Namespace) -> int:
    pd = _pd()
    path = args.h5mu or str(DATA_ROOT / "Gasperini_2019_sample_pilot.h5mu")
    mods = read_h5mu(path)
    rows = inspect_mods(mods)
    out = run_dir(args.label)
    write_tsv(pd.DataFrame(rows), out / "schema_checks.tsv")
    nf = sum(1 for r in rows if r["status"] == "FAIL")
    for r in rows:
        print(f"  {r['status']}  {r['check']} {r['detail']}")
    write_json(out / "summary.json", {"h5mu": path, "modalities": {k: list(v.shape) for k, v in mods.items()}, "n_checks": len(rows), "n_fail": nf})
    write_report(out / "report.md", [f"# MuData schema check: {Path(path).name}", "", f"{len(rows)} checks, {nf} failed.", "", "| status | check | detail |", "|---|---|---|"] +
                 [f"| {r['status']} | {r['check']} | {r['detail']} |" for r in rows])
    return 0 if nf == 0 else 1


def cmd_run(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    t0 = time.time()
    path = Path(args.h5mu) if args.h5mu else DATA_ROOT / "Gasperini_2019_sample_pilot.h5mu"
    if not path.exists():
        print(f"{path} not found: run `igvfagent gasperini-pipeline setup` or pass --h5mu")
        return 1
    mods = _load_mods(path)
    out = run_dir(args.label)
    print(f"[1/5] guide assignment (UMI > {args.guide_umi_limit}) and covariates")
    proc, cov = process_mudata(mods, args.guide_umi_limit, args.merge, args.upstream_compat)
    write_h5mu(out / "processed_mudata_guide_and_transcripts.h5mu", proc)
    write_tsv(cov, out / "covariates.tsv", index=True)
    print(f"[2/5] tested pairs (cis window {args.distance:,} bp{', in trans' if args.in_trans else ''})")
    rv = proc["scRNA"].var
    genes = gene_table_from_rna(rv, read_gtf_genes(args.gtf) if args.gtf else None)
    pairs, el_g, how = tested_pairs(list(proc["guides"].var_names), list(rv.index.astype(str)), genes, args.distance, args.in_trans,
                                    args.add_genes or [], args.upstream_compat, args.seed)
    write_tsv(pairs, out / "tested_pairs.tsv")
    print(f"      {len(el_g)} elements, {len(pairs)} element-gene pairs ({pd.Series(how).value_counts().to_dict()})")
    print(f"[3/5] SCEPTRE ({args.engine}, B = {args.resamples}, side = {args.direction})")
    res, skipped, znames = _run_tests(proc, pairs, el_g, cov, args, out)
    write_tsv(res, out / "results.tsv")
    if len(skipped):
        write_tsv(skipped, out / "skipped_guides.tsv")
    summ: Dict[str, Any] = {"h5mu": str(path), "n_cells": int(proc["scRNA"].n_obs), "n_genes": int(proc["scRNA"].n_vars),
                            "n_guides": int(proc["guides"].n_vars), "n_elements": len(el_g), "n_pairs": int(len(pairs)),
                            "n_tests": int(len(res)), "n_guides_skipped": int(len(skipped)), "covariates": znames,
                            "params": {k: getattr(args, k) for k in ("guide_umi_limit", "distance", "in_trans", "direction", "resamples", "engine", "upstream_compat", "merge")}}
    el = pd.DataFrame()
    if len(res):
        print("[4/5] BH / Fisher aggregation -> mudata_results.h5mu")
        res, el = aggregate_results(res, proc["guides"].var)
        write_tsv(res, out / "guide_gene_results.tsv")
        write_tsv(el, out / "element_gene_results.tsv")
        write_h5mu(out / "mudata_results.h5mu", results_mudata(res, el, proc))
        summ.update({"n_significant_fdr01": int(res["significant"].sum()), "n_element_gene": int(len(el)),
                     "n_element_sig_not_adj": int(el["sig_not_adj"].sum())})
        for t in ("POSITIVE_CONTROL", "NEGATIVE_CONTROL", "PUTATIVE_ENHANCER"):
            sub = res[res["guide_type"] == t]
            if len(sub):
                summ[f"{t.lower()}_frac_p_lt_0.05"] = float((sub["p_value"] < 0.05).mean())
    gas = None
    if args.gasperini_test:
        print("[5/5] original Gasperini NB likelihood-ratio test")
        r = proc["scRNA"]
        gas = gasperini_test(r, proc["guides"], pairs, el_g, r.obs["bath_number"] if "bath_number" in r.obs else None, r.obs.get("percent_mito"))
        write_tsv(gas, out / "gasperini_nb_results.tsv")
        summ["gasperini_nb_tests"] = int(len(gas))
        if len(gas) and len(el):
            j = el.merge(gas, left_on=["element", "gene_id"], right_on=["gRNA_group", "gene_short_name"])
            if len(j) > 2:
                from scipy.stats import spearmanr  # type: ignore
                summ["sceptre_vs_nb_spearman_log10p"] = float(spearmanr(-np.log10(np.maximum(j["p_value_fisher"], 1e-300)), -np.log10(np.maximum(j["pvalue.raw"], 1e-300)))[0])
    if args.compare:
        orig = pd.read_csv(args.compare, sep="\t")
        cmp_rows = {}
        if len(el):
            m, s = compare_to_original(el, orig, "p_value_fisher", "mean_log_fold_change", group_suffix=args.group_suffix, fdr=args.fdr)
            write_tsv(m, out / "compare_sceptre_vs_original.tsv")
            cmp_rows["sceptre"] = s
        if gas is not None and len(gas):
            m2, s2 = compare_to_original(gas.rename(columns={"pvalue.raw": "p_ours", "beta": "beta_ours"}), orig, "p_ours", "beta_ours",
                                         "gRNA_group", "gene_short_name", args.group_suffix, args.fdr)
            write_tsv(m2, out / "compare_nb_vs_original.tsv")
            cmp_rows["gasperini_nb"] = s2
        summ["comparison"] = cmp_rows
    plt = None if args.no_plots else _plt()
    if plt is not None and len(res):
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
        for t, col in (("NEGATIVE_CONTROL", GREY), (None, BLUE), ("POSITIVE_CONTROL", ORANGE), ("PUTATIVE_ENHANCER", "#1baf7a")):
            sub = res[res["guide_type"].isna()] if t is None else res[res["guide_type"] == t]
            if not len(sub):
                continue
            p = np.sort(np.maximum(sub["p_value"].values, 1e-300))
            e = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
            ax[0].scatter(e, -np.log10(p), s=6, color=col, label=t or "unlabelled")
        lim = ax[0].get_xlim()
        ax[0].plot(lim, lim, color=AXIS, lw=1)
        ax[0].set_xlabel("expected -log10 p")
        ax[0].set_ylabel("observed -log10 p")
        ax[0].legend(fontsize=7, frameon=False)
        ax[0].set_title("QQ by guide type", loc="left", fontsize=9)
        ax[1].scatter(res["log_fold_change"], -np.log10(np.maximum(res["p_value"], 1e-300)), s=5, color=BLUE)
        ax[1].set_xlabel("log fold change (NB MLE)")
        ax[1].set_ylabel("-log10 p")
        ax[1].set_title("volcano", loc="left", fontsize=9)
        save_fig(fig, out / "sceptre_qq_volcano.png")
        plt.close(fig)
    summ["seconds"] = round(time.time() - t0, 1)
    write_json(out / "summary.json", summ)
    L = [f"# Gasperini 2019 pipeline run: {path.name}", "",
         f"Upstream {UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]} (pipeline {PIPELINE_REPO}@{PIPELINE_COMMIT[:7]}).", "",
         f"- cells {summ['n_cells']}, genes {summ['n_genes']}, guides {summ['n_guides']}, elements {summ['n_elements']}",
         f"- guides assigned at UMI > {args.guide_umi_limit}; covariates: {', '.join(znames) if znames else 'n/a'}",
         f"- {summ['n_pairs']} element-gene pairs within {args.distance:,} bp; {summ['n_tests']} guide-gene tests; {summ['n_guides_skipped']} guides skipped (<= 30 cells)"]
    if "n_significant_fdr01" in summ:
        L.append(f"- {summ['n_significant_fdr01']} guide-gene tests at BH < 0.01; {summ['n_element_sig_not_adj']} of {summ['n_element_gene']} element-gene pairs at Fisher p < 0.05")
    for k in ("positive_control_frac_p_lt_0.05", "negative_control_frac_p_lt_0.05", "putative_enhancer_frac_p_lt_0.05"):
        if k in summ:
            L.append(f"- {k.replace('_', ' ')}: {_fmt(summ[k])}")
    if len(el):
        top = el.sort_values("p_value_fisher").head(15)
        L += ["", "## Top element-gene pairs", "", "| element | gene | Fisher p | guides | mean log FC | type |", "|---|---|---|---|---|---|"]
        L += [f"| {r.element} | {r.gene_id} | {_fmt(r.p_value_fisher)} | {r.n_guides} | {_fmt(r.mean_log_fold_change)} | {r.guide_type or ''} |" for r in top.itertuples()]
    for k, s in (summ.get("comparison") or {}).items():
        L += ["", f"## vs original Gasperini ({k})", ""] + _compare_lines(s)[2:]
    write_report(out / "report.md", L)
    print(f"Done in {summ['seconds']}s -> {out}")
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def _synthetic_world(d: Path, seed: int = 11):
    """Two lanes, 300 genes (3 MT-), 6 TSS elements (70% knockdown), 1 enhancer (50%), 2 NTC elements; doublets; empty drops."""
    ad, pd, np, sp = _ad(), _pd(), _np(), _sp()
    rng = np.random.default_rng(seed)
    n_genes = 300
    names = [f"MT-G{i}" if i < 3 else f"GENE{i}" for i in range(n_genes)]
    ids = [f"ENSG{i:011d}" for i in range(n_genes)]
    chrs = ["1" if i < 200 else "2" for i in range(n_genes)]
    tss = np.array([100_000 + (i % 200) * 150_000 for i in range(n_genes)])
    strand = ["+" if i % 2 else "-" for i in range(n_genes)]
    with open(d / "genes.gtf", "w") as fh:
        for i in range(n_genes):
            s, e = (tss[i], tss[i] + 5000) if strand[i] == "+" else (tss[i] - 5000, tss[i])
            fh.write(f"{chrs[i]}\tsyn\tgene\t{s}\t{e}\t.\t{strand[i]}\t.\tgene_id \"{ids[i]}.1\"; gene_name \"{names[i]}\";\n")
    tss_targets = [10, 20, 30, 40, 50, 60]
    enh_target = 120
    guides = []
    for t in tss_targets:
        for k in range(2):
            guides.append({"Target_name": f"GENE{t}_TSS", "sgRNA_sequences": "".join(rng.choice(list("ACGT"), 20)), "chr": "chr1", "start": int(tss[t]) + 50 * k, "end": int(tss[t]) + 50 * k + 20})
    for k in range(2):
        guides.append({"Target_name": "chr1_enh1", "sgRNA_sequences": "".join(rng.choice(list("ACGT"), 20)), "chr": "chr1", "start": int(tss[enh_target]) + 30_000 + 40 * k, "end": int(tss[enh_target]) + 30_020 + 40 * k})
    for nm in ("random_1", "scrambled_2"):
        for k in range(2):
            guides.append({"Target_name": nm, "sgRNA_sequences": "".join(rng.choice(list("ACGT"), 20)), "chr": "NTC", "start": 0, "end": 20})
    gt = pd.DataFrame(guides)
    gt.to_csv(d / "guides.tsv", sep="\t", index=False)
    feats = guide_features(gt)
    gnames = list(feats["pipeline_id"])
    gel = [n.split("|")[0] for n in gnames]
    base = rng.lognormal(0.2, 1.0, n_genes)
    base[:3] = 1.5
    base[tss_targets] = 3.0
    base[enh_target] = 3.0
    lanes = {}
    truth_doublets = []
    for lane in (1, 2):
        n_cells, n_empty = 700, 150
        bcs = [f"L{lane}C{i:05d}" for i in range(n_cells + n_empty)]
        sf = np.r_[rng.lognormal(math.log(4.0), 0.3, n_cells), rng.lognormal(math.log(0.05), 0.3, n_empty)]
        G = np.zeros((n_cells + n_empty, len(gnames)))
        for c in range(n_cells):
            k = rng.poisson(2.5)
            pick = rng.choice(len(gnames), size=min(k, len(gnames)), replace=False)
            G[c, pick] = rng.integers(5, 40, len(pick))
        G += rng.poisson(0.3, G.shape)
        G[0, 0] = 3  # exactly GUIDE_UMI_LIMIT UMIs: not assigned (strict >)
        mu = sf[:, None] * base[None, :]
        ctype = rng.random(n_cells + n_empty) < 0.5  # two cell types with distinct programs (doublets become heterotypic)
        mu[np.ix_(ctype, np.arange(130, 200))] *= 6
        mu[np.ix_(~ctype, np.arange(200, 300))] *= 6
        has = (G > 3)
        for j, t in enumerate(tss_targets):
            idx = [i for i, e in enumerate(gel) if e == f"GENE{t}_TSS"]
            mu[has[:, idx].any(1), t] *= 0.3
        idx = [i for i, e in enumerate(gel) if e == "chr1_enh1"]
        mu[has[:, idx].any(1), enh_target] *= 0.5
        theta = 5.0
        Y = rng.negative_binomial(theta, theta / (theta + mu)).astype(float)
        hi_mito = rng.choice(n_cells, 15, replace=False)
        Y[hi_mito, :3] += rng.poisson(1500, (15, 3))
        dbl = [c for c in rng.choice(n_cells, 25, replace=False) if c not in set(hi_mito)]
        for c in dbl:
            others = np.where((ctype[:n_cells] != ctype[c]) & ~np.isin(np.arange(n_cells), hi_mito))[0]
            o = rng.choice(others)
            Y[c] = Y[c] + Y[o]
        truth_doublets += [bcs[c] for c in dbl]
        rna = ad.AnnData(X=sp.csr_matrix(Y), obs=pd.DataFrame(index=bcs), var=pd.DataFrame({"feature_name": names}, index=[f"{x}.1" for x in ids]))
        gd = ad.AnnData(X=sp.csr_matrix(G), obs=pd.DataFrame(index=bcs), var=pd.DataFrame({"feature_name": gnames}, index=gnames))
        for kind, a in (("transcripts", rna), ("guide", gd)):
            p = d / f"S1_L{lane}_ks_{kind}_out" / "counts_unfiltered"
            p.mkdir(parents=True)
            a.write_h5ad(p / "adata.h5ad")
        lanes[lane] = {"n_cells": n_cells, "hi_mito": [bcs[c] for c in hi_mito]}
    return {"tss_targets": tss_targets, "enh_target": enh_target, "names": names, "feats": feats, "doublets": truth_doublets, "lanes": lanes}


def _synthetic_multiseq(d: Path, rng, n_cells: int = 240, n_tags: int = 4):
    ref = ["".join(rng.choice(list("ACGT"), 8)) for _ in range(n_tags)]
    cells = ["".join(rng.choice(list("ACGT"), 16)) for _ in range(n_cells)]
    truth = {}
    with gzip.open(d / "ms_R1.fq.gz", "wt") as f1, gzip.open(d / "ms_R2.fq.gz", "wt") as f2:
        k = 0
        for i, c in enumerate(cells):
            if i < 20:
                truth[c], tags = "Negative", []
            elif i < 40:
                truth[c], tags = "Doublet", [i % n_tags, (i + 1) % n_tags]
            else:
                truth[c], tags = f"Bar{i % n_tags + 1}", [i % n_tags]
            n_amb = int(rng.integers(3, 8))
            reads = [(t, int(rng.integers(60, 120))) for t in tags] + [(int(rng.integers(0, n_tags)), 1) for _ in range(n_amb)]
            for t, n in reads:
                for _ in range(n):
                    umi = "".join(rng.choice(list("ACGT"), 12))
                    tag = list(ref[t])
                    if rng.random() < 0.2:
                        tag[int(rng.integers(0, 8))] = "N"
                    s1 = c + umi
                    s2 = "".join(tag) + "".join(rng.choice(list("ACGT"), 30))
                    f1.write(f"@m{k}\n{s1}\n+\n{'F' * len(s1)}\n")
                    f2.write(f"@m{k}\n{s2}\n+\n{'F' * len(s2)}\n")
                    k += 1
    (d / "ms_cells.txt").write_text("\n".join(cells) + "\n")
    (d / "ms_ref.csv").write_text("\n".join(ref) + "\n")
    return truth


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    pd, np = _pd(), _np()
    checks: List[Tuple[bool, str]] = []
    t0 = time.time()

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made: List[Path] = []

    def newest(label):
        runs = sorted(OUT_ROOT.glob(f"*_{label}*"))
        made.extend(r for r in runs if r not in made)
        return runs[-1] if runs else None

    ns = argparse.Namespace
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = _synthetic_world(d)

        print("\nguide table + composition")
        rc = cmd_guide_table(ns(table=str(d / "guides.tsv"), label="st_gp_guides"))
        run = newest("st_gp_guides")
        feat = pd.read_csv(run / "guide_features.txt", sep="\t", header=None)
        check(rc == 0 and len(feat) == 18 and "GENE10_TSS|2_sgrna_chr1:1600050:1600070" in set(feat[1]),
              "guide-table: 18 guides, '<Target>|<n>_sgrna_<chr>:<start>:<end>' ids numbered within target")
        rng = np.random.default_rng(3)
        with gzip.open(d / "comp.fq.gz", "wt") as fh:
            for i in range(500):
                s = "ACGTACGTAC" + "".join(rng.choice(list("ACGT"), 20))
                fh.write(f"@r{i}\n{s}\n+\n{'F' * 30}\n")
        rc = cmd_composition(ns(fastq=[str(d / "comp.fq.gz")], n_lines=10000, upstream_compat=False, no_plots=True, label="st_gp_comp"))
        run = newest("st_gp_comp")
        comp = pd.read_csv(run / "comp.fq.gz_composition.tsv", sep="\t")
        check(rc == 0 and comp["sd"].iloc[:10].min() > 200 and comp["sd"].iloc[10:].max() < 60 and comp["n_reads"].iloc[0] == 500,
              "composition: 10 fixed positions have high SD across A/C/T/G, random positions low (2,500 lines -> 500 reads)")

        print("\ncell QC + lane merge")
        rna_dirs = [str(d / f"S1_L{i}_ks_transcripts_out") for i in (1, 2)]
        g_dirs = [str(d / f"S1_L{i}_ks_guide_out") for i in (1, 2)]
        qa = ns(rna=rna_dirs, guide=g_dirs, lanes=["1", "2"], expected_cell_number=699, mito_expected_percentage=0.2,
                transcripts_umi_threshold=200, percentage_of_cells_to_include_transcript=0.01, upstream_compat=False, no_doublets=False, seed=0, label="st_gp_qc")
        rc = cmd_qc_filter(qa)
        run = newest("st_gp_qc")
        ad = _ad()
        rq = ad.read_h5ad(run / "results_per_lane" / "full_raw_scrna_ann_data.h5ad")
        gq = ad.read_h5ad(run / "results_per_lane" / "full_raw_guide_ann_data.h5ad")
        lq = pd.read_csv(run / "lane_qc.tsv", sep="\t")
        empties = sum(1 for o in rq.obs_names if int(o.split("C")[1][:5]) >= 700)
        hm = set(W["lanes"][1]["hi_mito"]) | set(W["lanes"][2]["hi_mito"])
        check(rc == 0 and empties == 0 and not (hm & set(rq.obs_names)), "qc-filter: empty droplets (>= 200 UMIs and >= knee[E], which keeps the E + 1 largest barcodes) and >20% mito cells removed")
        check(set(rq.obs["batch_number"]) == {1, 2} and list(rq.obs_names) == list(gq.obs_names) and "number_of_nonzero_guides" in gq.obs,
              "qc-filter: lanes concatenated with batch_number, barcodes shared by both modalities, number_of_nonzero_guides")
        dbl = set(W["doublets"])
        gone = len(dbl) - len(dbl & set(rq.obs_names))
        called = int(lq["n_doublets"].sum())
        check(gone >= 0.6 * len(dbl) and called < 0.1 * 1400, f"qc-filter: scrublet-style doublets ({called} called, {gone}/{len(dbl)} planted heterotypic doublets removed)")
        check(all(rq.var_names.str.startswith("ENSG")) and not any("." in v for v in rq.var_names) and (rq.var["n_cells"] >= int(rq.n_obs * 0.01)).all(),
              "qc-filter: Ensembl versions stripped, genes kept in >= 1% of cells")
        qa2 = ns(**{**vars(qa), "upstream_compat": True, "no_doublets": True, "transcripts_umi_threshold": 280, "label": "st_gp_qccompat"})
        cmd_qc_filter(qa2)
        run2 = newest("st_gp_qccompat")
        lq2 = pd.read_csv(run2 / "lane_qc.tsv", sep="\t")
        _, _, st280 = analyze_batch(read_counts(rna_dirs[0]), read_counts(g_dirs[0]), 1, 699, 0.2, 280, False, False)
        check(st280["n_after_umi_knee"] == 700 and (lq2["n_after_umi_knee"] < 700).all(),
              f"qc-filter --upstream-compat: the UMI cutoff applied as min_genes removes cells ({lq2['n_after_umi_knee'].tolist()} vs 700 kept at >= 280 UMIs)")

        print("\nMuData + guide assignment")
        rc = cmd_mudata(ns(rna_h5ad=str(run / "results_per_lane" / "full_raw_scrna_ann_data.h5ad"), guide_h5ad=str(run / "results_per_lane" / "full_raw_guide_ann_data.h5ad"),
                           gtf=str(d / "genes.gtf"), guide_features=str(newest("st_gp_guides") / "guide_features.txt"), upstream_compat=False, label="st_gp_mudata"))
        mrun = newest("st_gp_mudata")
        raw = mrun / "raw_mudata_guide_and_transcripts.h5mu"
        mods = read_h5mu(raw)
        gv, rv = mods["guides"].var, mods["scRNA"].var
        n0 = "GENE10_TSS|1_sgrna_chr1:1600000:1600020"
        check(rc == 0 and gv.loc[n0, "guide_chr"] == "chr1" and gv.loc[n0, "guide_start"] == "1600000" and gv.loc[n0, "guide_end"] == "1600020"
              and gv.loc[n0, "guide_number"] == "1" and gv.loc[n0, "target_elements"] == "GENE10_TSS" and len(gv.loc[n0, "sequences"]) == 20,
              "mudata: guide var parsed from the name (chr, start, end, number, target) + sequences")
        e10 = rv.index[rv["feature_name"] == "GENE10"][0]
        e11 = rv.index[rv["feature_name"] == "GENE11"][0]
        check(rv.loc[e10, "transcript_start"] == 1_600_000 and rv.loc[e10, "transcript_end"] == 1_600_001 and rv.loc[e11, "transcript_start"] == 1_750_000,
              "mudata: transcript_start = TSS (gene end on '-', start on '+'), transcript_end = TSS + 1")
        mods_c = {"scRNA": annotate_transcripts(mods["scRNA"].copy(), read_gtf_genes(d / "genes.gtf"), upstream_compat=True)}
        check(mods_c["scRNA"].var.loc[e11, "transcript_end"] == 1_755_000, "mudata --upstream-compat: transcript_end = GTF gene end")
        ins = inspect_mods(mods)
        check(all(r["status"] == "PASS" for r in ins), f"inspect: raw MuData matches the sample-pilot schema ({len(ins)} checks)")
        rc = cmd_assign(ns(h5mu=str(raw), guide_umi_limit=3, merge=False, upstream_compat=False, no_plots=True, label="st_gp_assign"))
        arun = newest("st_gp_assign")
        pm = read_h5mu(arun / "processed_mudata_guide_and_transcripts.h5mu")
        G, U = _dense(pm["guides"].X), _dense(pm["guides"].layers["guide_umi"])
        check(rc == 0 and np.array_equal(G, (U > 3).astype(G.dtype)) and set(np.unique(G)) <= {0, 1},
              "assign: binary = UMI > GUIDE_UMI_LIMIT (3, strict); raw UMIs kept in layers['guide_umi']")
        cov = pd.read_csv(arun / "covariates.tsv", sep="\t", index_col=0)
        o = pm["scRNA"].obs
        check(np.allclose(cov["log_total_gene_count"], np.log(o["n_counts"] + 1)) and np.allclose(cov["log_number_of_detected_genes"], np.log(o["n_genes"]))
              and np.allclose(cov["log_total_guide_count"], np.log(G.sum(1) + 1)), "covariates: log n_genes, log(n_counts + 1), log(assigned guides + 1)")
        pc, covc = process_mudata(mods, 3, False, upstream_compat=True)
        check(np.allclose(covc["log_total_gene_count"], np.log(o["n_genes"].values + 1)), "covariates --upstream-compat: log_total_gene_count = log(n_genes + 1)")
        pmg, _ = process_mudata(mods, 3, merge=True)
        check(pmg["guides"].n_vars == 9 and all("|2_" in v or "|2" in v for v in pmg["guides"].var_names),
              "assign --merge: guides summed per target (9 targets), last guide's coordinates kept")

        print("\ntested pairs")
        rc = cmd_pairs(ns(h5mu=str(arun / "processed_mudata_guide_and_transcripts.h5mu"), gtf=None, distance=1_000_000, in_trans=False, add_genes=None,
                          upstream_compat=False, seed=0, label="st_gp_pairs"))
        prun = newest("st_gp_pairs")
        pr = pd.read_csv(prun / "tested_pairs.tsv", sep="\t")
        g10 = set(pr.loc[pr["element"] == "GENE10_TSS", "gene_id"])
        check(rc == 0 and "GENE10" in g10 and all(abs(int(n[4:]) - 10) <= 6 for n in g10 if n.startswith("GENE")),
              f"pairs: TSS element tests genes within 1 Mb of its gene's TSS ({len(g10)} genes)")
        enh = set(pr.loc[pr["element"] == "chr1_enh1", "gene_id"])
        check("GENE120" in enh and set(pr.loc[pr["element"] == "chr1_enh1", "selection"]) == {"guide_coordinates"}, "pairs: enhancer element uses its guide coordinates (corrected)")
        ntc = pr[pr["element"] == "random_1"]
        check(len(ntc) == 10 and set(ntc["guide_type"]) == {"NEGATIVE_CONTROL"} and set(pr.loc[pr["element"] == "GENE10_TSS", "guide_type"]) == {"POSITIVE_CONTROL"},
              "pairs: NTC elements get 10 random genes; GUIDE_TYPE labels")
        rv2 = pm["scRNA"].var
        pc2, _, how2 = tested_pairs(list(pm["guides"].var_names), list(rv2.index), gene_table_from_rna(rv2), upstream_compat=True)
        check(how2["chr1_enh1"] == "random_genes" and (pc2["element"] == "chr1_enh1").sum() == 10, "pairs --upstream-compat: enhancer element gets 10 random genes")

        print("\nSCEPTRE CRT + aggregation")
        rc = cmd_sceptre(ns(h5mu=str(arun / "processed_mudata_guide_and_transcripts.h5mu"), pairs=str(prun / "tested_pairs.tsv"), engine="python",
                            resamples=300, direction="both", seed=0, label="st_gp_sceptre"))
        srun = newest("st_gp_sceptre")
        res = pd.read_csv(srun / "results.tsv", sep="\t")
        tgt = res[[e.endswith("_TSS") and g == e.replace("_TSS", "") for e, g in zip(res["element"], res["gene_id"])]]
        check(rc == 0 and len(tgt) == 12 and (tgt["p_value"] < 1e-3).all() and (tgt["z_value"] < 0).all(),
              f"sceptre: all 12 positive-control guide->own-gene tests p < 1e-3 with z < 0 (max p {_fmt(tgt['p_value'].max())})")
        check(abs(np.median(tgt["log_fold_change"]) - math.log(0.3)) < 0.4, f"sceptre: log fold change ~ log(0.3) = -1.20 (median {_fmt(np.median(tgt['log_fold_change']))})")
        et = res[(res["element"] == "chr1_enh1") & (res["gene_id"] == "GENE120")]
        check(len(et) == 2 and (et["p_value"] < 0.01).all(), "sceptre: planted enhancer effect (50%) detected for both guides")
        nul = res[(res["guide_type"] == "NEGATIVE_CONTROL") | ((res["guide_type"] == "POSITIVE_CONTROL") & (res["gene_id"] != res["element"].str.replace("_TSS", "")))]
        check(len(nul) > 50 and (nul["p_value"] < 0.05).mean() < 0.12, f"sceptre: null pairs calibrated ({(nul['p_value'] < 0.05).mean():.3f} at p < 0.05, n = {len(nul)})")
        check(set(res.columns) >= {"gene_id", "gRNA_id", "pair_type", "p_value", "z_value", "log_fold_change"} and set(res["pair_type"]) == {"all_elements_test"},
              "sceptre: upstream result columns, pair_type 'all_elements_test'")
        from scipy.stats import combine_pvalues  # type: ignore
        agg, el = aggregate_results(res, pm["guides"].var)
        r10 = el[(el["element"] == "GENE10_TSS") & (el["gene_id"] == "GENE11")].iloc[0]
        ps = res[(res["element"] == "GENE10_TSS") & (res["gene_id"] == "GENE11")]["p_value"].tolist()
        check(abs(r10["p_value_fisher"] - combine_pvalues(ps).pvalue) < 1e-12 and np.allclose(agg["adj_pvalue"], bh(res["p_value"])),
              "results: element p = Fisher combination of its guides; BH over all tests")
        rc = cmd_results(ns(h5mu=str(arun / "processed_mudata_guide_and_transcripts.h5mu"), results=str(srun / "results.tsv"), label="st_gp_results"))
        rrun = newest("st_gp_results")
        rm = read_h5mu(rrun / "mudata_results.h5mu")
        check(rc == 0 and {"result_guides", "result_elements"} <= set(rm) and {"adj_pvalue", "z_value", "log_fold_change", "significant"} <= set(rm["result_guides"].layers)
              and "sig_not_adj" in rm["result_elements"].layers, "results: mudata_results.h5mu with result_guides (4 layers) and result_elements")
        cov_all = _cov_from(pm)
        runs = write_sceptre_r_inputs(d / "rtest", pm["scRNA"], pm["guides"], cov_all, pr, {e: [g for g in pm["guides"].var_names if g.startswith(e + "|")] for e in set(pr["element"])})
        rd = Path(runs[0]["dir"]) if runs else d
        check(runs and all((rd / f).exists() for f in ("gene_exp_one_gene_guide.txt", "guides_one_gene_guide.txt", "covariates.txt", "pairs.txt", "r_script.r"))
              and "run_sceptre_high_moi" in (rd / "r_script.r").read_text(), "sceptre --engine r: per-element SCEPTRE inputs + R script written")

        print("\nGasperini NB test + comparison + run")
        rc = cmd_gasperini_test(ns(h5mu=str(arun / "processed_mudata_guide_and_transcripts.h5mu"), pairs=str(prun / "tested_pairs.tsv"), label="st_gp_nb"))
        nrun = newest("st_gp_nb")
        nb = pd.read_csv(nrun / "gasperini_nb_results.tsv", sep="\t")
        own = nb[[e.endswith("_TSS") and g == e.replace("_TSS", "") for e, g in zip(nb["gRNA_group"], nb["gene_short_name"])]]
        check(rc == 0 and len(own) == 6 and (own["pvalue.raw"] < 1e-8).all() and np.allclose(own["fold_change.transcript_remaining"], 0.3, atol=0.1),
              f"gasperini-test: self-TSS LRT p < 1e-8, fold_change.transcript_remaining ~ 0.3 (median {_fmt(own['fold_change.transcript_remaining'].median())})")
        check(pd.to_numeric(nb["pvalue.empirical"], errors="coerce").notna().all(), "gasperini-test: empirical p-values against the NTC groups")
        orig = nb.copy()
        orig["pvalue.raw"] = np.clip(orig["pvalue.raw"] * np.exp(np.random.default_rng(1).normal(0, 0.5, len(orig))), 0, 1)
        orig["site_type"] = ["selfTSS" if e.replace("_TSS", "") == g else "TSS" for e, g in zip(orig["gRNA_group"], orig["gene_short_name"])]
        orig.to_csv(d / "orig.tsv", sep="\t", index=False)
        pd.DataFrame(el).to_csv(d / "el.tsv", sep="\t", index=False)
        rc = cmd_compare(ns(results=str(d / "el.tsv"), original=str(d / "orig.tsv"), group_suffix="_TSS", fdr=0.1, label="st_gp_compare"))
        crun = newest("st_gp_compare")
        cs = json.loads((crun / "summary.json").read_text())
        check(rc == 0 and cs["n_matched"] > 50 and cs["spearman_rho_log10p"] > 0.3 and cs.get("self_tss_hits_ours", 0) == 6,
              f"compare: matched on element+'_TSS' and gene; rho = {_fmt(cs.get('spearman_rho_log10p'))}, 6/6 self-TSS hits recovered")
        rc = cmd_run(ns(h5mu=str(raw), gtf=None, guide_umi_limit=3, merge=False, distance=1_000_000, in_trans=False, add_genes=None, engine="python",
                        resamples=200, direction="both", seed=0, upstream_compat=False, gasperini_test=True, compare=str(d / "orig.tsv"),
                        group_suffix="_TSS", fdr=0.1, no_plots=args.no_plots, label="st_gp_run"))
        rrun2 = newest("st_gp_run")
        sj = json.loads((rrun2 / "summary.json").read_text())
        check(rc == 0 and sj["positive_control_frac_p_lt_0.05"] > 0.03 and sj["sceptre_vs_nb_spearman_log10p"] > 0.5 and "sceptre" in sj["comparison"]
              and (rrun2 / "mudata_results.h5mu").exists() and (rrun2 / "report.md").exists(),
              f"run: end-to-end on the raw MuData; SCEPTRE vs NB rho {_fmt(sj.get('sceptre_vs_nb_spearman_log10p'))}")

        print("\nMULTI-seq")
        truth = _synthetic_multiseq(d, np.random.default_rng(5))
        rc = cmd_multiseq(ns(r1=str(d / "ms_R1.fq.gz"), r2=str(d / "ms_R2.fq.gz"), cell_barcodes=str(d / "ms_cells.txt"), bar_ref=str(d / "ms_ref.csv"),
                             bar_multi=[1, 16], umi_multi=[17, 28], r2_multi_tag=[1, 8], upstream_compat=False, h5mu=None, label="st_gp_multiseq"))
        mrun2 = newest("st_gp_multiseq")
        fc = pd.read_csv(mrun2 / "final_class.tsv", sep="\t").set_index("cell_barcode")["multiseq_class"]
        acc = np.mean([fc[c] == t for c, t in truth.items() if t.startswith("Bar")])
        dacc = np.mean([fc[c] == "Doublet" for c, t in truth.items() if t == "Doublet"])
        nacc = np.mean([fc[c] == "Negative" for c, t in truth.items() if t == "Negative"])
        check(rc == 0 and acc > 0.9 and dacc > 0.7 and nacc > 0.7, f"multiseq: singlets {acc:.2f}, doublets {dacc:.2f}, negatives {nacc:.2f} recovered (Hamming <= 1 tags)")
        bt = pd.read_csv(mrun2 / "bar_table.tsv", sep="\t", index_col=0)
        check(list(bt.columns[-2:]) == ["nUMI", "nUMI_total"] and (bt["nUMI"] <= bt["nUMI_total"]).all(), "multiseq: bar table with nUMI and nUMI_total columns")
        calls_c, _ = multiseq_classify(bt, upstream_compat=True)
        check(len(calls_c) == len(bt), "multiseq --upstream-compat: nUMI_total kept as a pseudo-barcode (runs)")

        print("\nconfig / kb / setup")
        rc = cmd_config(ns(kallisto_bin=None, guide_features="guide_features.txt", chemistry="10XV2", rna_fastqs=["a_R1.fq.gz", "a_R2.fq.gz"], guide_fastqs=None,
                           whitelist="737K-august-2016.txt", main_nf="pipeline_perturbseq_like/main.nf", workdir="work", label="st_gp_config"))
        cfg = (newest("st_gp_config") / "gasperini_sample.config").read_text()
        check(rc == 0 and "params.GUIDE_UMI_LIMIT = 3" in cfg and "params.TRANSCRIPTS_UMI_TRHESHOLD = 200" in cfg and "params.EXPECTED_CELL_NUMBER = 10000" in cfg,
              "config: README parameters (GUIDE_UMI_LIMIT 3, UMI threshold 200, 10,000 cells)")
        kc = kb_commands(["r1", "r2"], ["g1", "g2"], "guide_features.txt", "wl.txt", "10XV2", 4)
        check(kc[1][-1] == "guide_features.txt" and "kite" in kc[3] and "-x" in kc[2], "kb: kite guide index + kb count commands")
        check(sum(v[1] for v in SETUP_FILES.values()) < 200e6, "setup: default downloads < 200 MB")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.0f}s)")
    print("selftest: all checks pass" if ok else f"selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent gasperini-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch the guide table, sample MuData and GEO pilot results (~48 MB)")
    p.add_argument("--dest", help=f"destination (default {DATA_ROOT})")
    p.add_argument("--no-geo", action="store_true", help="skip the GEO GSE120861 files")
    p.add_argument("--whitelist", action="store_true", help="also fetch 737K-august-2016.txt (~13 MB)")
    p.add_argument("--print-raw", action="store_true", help="print the SRA BAM + bamtofastq commands (several GB; not run)")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("guide-table", help="guide table (xlsx/tsv) -> guide_features.txt")
    p.add_argument("--table", required=True, help="columns Target_name, sgRNA_sequences, chr, start, end")
    p.add_argument("--label", default="guide_table")
    p.set_defaults(func=cmd_guide_table)

    p = sub.add_parser("composition", help="per-position nucleotide composition of FASTQs")
    p.add_argument("--fastq", nargs="+", required=True)
    p.add_argument("--n-lines", type=int, default=10000)
    p.add_argument("--upstream-compat", action="store_true", help="upstream's line filter instead of FASTQ record parsing")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="composition")
    p.set_defaults(func=cmd_composition)

    p = sub.add_parser("kb", help="kb ref / kb count commands (run with --execute when kb is on PATH)")
    p.add_argument("--rna-fastqs", nargs="+")
    p.add_argument("--guide-fastqs", nargs="+")
    p.add_argument("--guide-features", default="guide_features.txt")
    p.add_argument("--whitelist", default="737K-august-2016.txt")
    p.add_argument("--chemistry", default="10XV2")
    p.add_argument("--threads", type=int, default=10)
    p.add_argument("--name", default="S1_L1")
    p.add_argument("--execute", action="store_true")
    p.set_defaults(func=cmd_kb)

    p = sub.add_parser("config", help="write the Nextflow gasperini_sample.config")
    p.add_argument("--guide-features", default="df_from_gasperini_tss.xlsx")
    p.add_argument("--chemistry", default="10XV2")
    p.add_argument("--rna-fastqs", nargs="+")
    p.add_argument("--guide-fastqs", nargs="+")
    p.add_argument("--whitelist", default="737K-august-2016.txt")
    p.add_argument("--kallisto-bin")
    p.add_argument("--main-nf", default="pipeline_perturbseq_like/main.nf")
    p.add_argument("--workdir", default="gasperini_test_01")
    p.add_argument("--label", default="config")
    p.set_defaults(func=cmd_config)

    p = sub.add_parser("qc-filter", help="per-lane cell QC + lane merge")
    p.add_argument("--rna", nargs="+", required=True, help="kb transcript output dirs or .h5ad, one per lane")
    p.add_argument("--guide", nargs="+", required=True, help="kb guide output dirs or .h5ad, one per lane")
    p.add_argument("--lanes", nargs="+", help="lane numbers (default 1..n)")
    p.add_argument("--expected-cell-number", type=int, default=CONFIG["EXPECTED_CELL_NUMBER"])
    p.add_argument("--mito-expected-percentage", type=float, default=CONFIG["MITO_EXPECTED_PERCENTAGE"])
    p.add_argument("--transcripts-umi-threshold", type=int, default=CONFIG["TRANSCRIPTS_UMI_TRHESHOLD"])
    p.add_argument("--percentage-of-cells-to-include-transcript", type=float, default=CONFIG["PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS"])
    p.add_argument("--no-doublets", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="qc_filter")
    p.set_defaults(func=cmd_qc_filter)

    p = sub.add_parser("mudata", help="guide + scRNA AnnData -> raw MuData")
    p.add_argument("--rna-h5ad", required=True)
    p.add_argument("--guide-h5ad", required=True)
    p.add_argument("--gtf", help="GTF(.gz) for gene TSS coordinates")
    p.add_argument("--guide-features", help="guide_features.txt (adds guides.var['sequences'])")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="mudata")
    p.set_defaults(func=cmd_mudata)

    p = sub.add_parser("assign", help="binarise guides, covariates -> processed MuData")
    p.add_argument("--h5mu", required=True)
    p.add_argument("--guide-umi-limit", type=int, default=CONFIG["GUIDE_UMI_LIMIT"])
    p.add_argument("--merge", action="store_true")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="assign")
    p.set_defaults(func=cmd_assign)

    def pair_args(p):
        p.add_argument("--gtf", help="GTF for gene TSS (default: the MuData's transcript_chr / transcript_start)")
        p.add_argument("--distance", type=int, default=CONFIG["DISTANCE_NEIGHBORS"])
        p.add_argument("--in-trans", action="store_true")
        p.add_argument("--add-genes", nargs="+")
        p.add_argument("--seed", type=int, default=0)
        p.add_argument("--upstream-compat", action="store_true")

    p = sub.add_parser("pairs", help="element x guide and element x tested-gene tables")
    p.add_argument("--h5mu", required=True)
    pair_args(p)
    p.add_argument("--label", default="pairs")
    p.set_defaults(func=cmd_pairs)

    def test_args(p):
        p.add_argument("--engine", default="python", choices=["python", "r"])
        p.add_argument("--resamples", type=int, default=500)
        p.add_argument("--direction", default=CONFIG["DIRECTION"], choices=["both", "left", "right"])

    p = sub.add_parser("sceptre", help="per-guide cis tests")
    p.add_argument("--h5mu", required=True, help="processed MuData (assign output)")
    p.add_argument("--pairs", required=True, help="tested_pairs.tsv")
    test_args(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="sceptre")
    p.set_defaults(func=cmd_sceptre)

    p = sub.add_parser("results", help="BH / Fisher aggregation -> mudata_results.h5mu")
    p.add_argument("--h5mu", required=True)
    p.add_argument("--results", required=True)
    p.add_argument("--label", default="results")
    p.set_defaults(func=cmd_results)

    p = sub.add_parser("gasperini-test", help="original Gasperini 2019 NB likelihood-ratio test")
    p.add_argument("--h5mu", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--label", default="gasperini_nb")
    p.set_defaults(func=cmd_gasperini_test)

    def cmp_args(p):
        p.add_argument("--group-suffix", default="_TSS", help="element -> gRNA_group suffix")
        p.add_argument("--fdr", type=float, default=0.1)

    p = sub.add_parser("compare", help="concordance with the original Gasperini results")
    p.add_argument("--results", required=True, help="element_gene_results.tsv, results.tsv or gasperini_nb_results.tsv")
    p.add_argument("--original", required=True, help="GSE120861_all_deg_results.pilot.txt.gz (or at_scale)")
    cmp_args(p)
    p.add_argument("--label", default="compare")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("multiseq", help="MULTI-seq sample demultiplexing")
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--cell-barcodes", required=True, help="one barcode per line (cells of the MuData)")
    p.add_argument("--bar-ref", required=True, help="CSV, first column = MULTI-seq barcodes")
    p.add_argument("--bar-multi", type=int, nargs=2, default=[1, 16])
    p.add_argument("--umi-multi", type=int, nargs=2, default=[17, 28])
    p.add_argument("--r2-multi-tag", type=int, nargs=2, default=[1, 8])
    p.add_argument("--h5mu", help="processed MuData to filter (Doublet / Negative removed)")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="multiseq")
    p.set_defaults(func=cmd_multiseq)

    p = sub.add_parser("inspect", help="check a MuData against the sample-pilot schema")
    p.add_argument("--h5mu", help="default: the fetched sample pilot")
    p.add_argument("--label", default="inspect")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("run", help="assign -> pairs -> sceptre -> results [-> gasperini-test -> compare]")
    p.add_argument("--h5mu", help="MuData with raw guide UMIs (default: the fetched sample pilot)")
    p.add_argument("--guide-umi-limit", type=int, default=CONFIG["GUIDE_UMI_LIMIT"])
    p.add_argument("--merge", action="store_true")
    pair_args(p)
    test_args(p)
    p.add_argument("--max-elements", type=int, help="test only the first N elements (quick look)")
    p.add_argument("--gasperini-test", action="store_true", help="also run the original NB LRT")
    p.add_argument("--compare", help="original results table (GSE120861_all_deg_results.pilot.txt.gz)")
    cmp_args(p)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="gasperini_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic screen, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
