#!/usr/bin/env python3
"""CRISPR-FG scPerturb-seq jamboree tasks: assay specs, Cell Ranger inputs, MuData QC, guide calling, differential perturbation and genome tracks (port of IGVF-CRISPR/CRISPR_FG_JAMBOREE).

Port of https://github.com/IGVF-CRISPR/CRISPR_FG_JAMBOREE (No LICENSE file;
Jupyter notebooks + Nextflow + Python stubs), the April-2023 IGVF CRISPR-FG
jamboree that asked teams to extend the scPERTURB-seq pipeline
LucasSilvaFerreira/pipeline_perturbseq_like.  Every notebook was read:
run_jamboree_v2.ipynb (pipeline config + launch), Task1_assay_specification
(assay -> kallisto read format + whitelist), task_cellranger_v2 (Cell Ranger
library.csv / feature_ref.csv), Task2_mudata_extension (MuData QC, subsetting,
controls, guide-gene distances via IGVF-CRISPR/sccrispr-tools), the guide
assignment task (binarised guide layer), Task3_differential_perturbation_v3
(an alternative differential-perturbation module over the processed MuData)
and Task4_visualization (BED + pyGenomeTracks link tracks from the results
MuData).  The tasks are stubs whose inputs and outputs are defined by the
pipeline, so the pipeline steps that produce and consume them were read too
(filtering_and_lane_merging.py, muon_creation.py, guide_table_processing.py,
PerturbLoader_generation.py, runSceptre.py, sceptre_anndata_creation.py @
60cd461a) and their definitions reproduced.  No LICENSE file upstream, so the
METHODS were ported, not the code.  Relationship: port.

Definitions, exactly as upstream computes them
  assay spec       kallisto bus -x "bc:umi:seq" triplets (file, start, stop);
                   10xv2 0,0,16:0,16,26:1,0,0 + 737K-august-2016, 10xv3
                   0,0,16:0,16,28:1,0,0 + 3M-february-2018; 'custom' uses
                   CHEMISTRY + WHITELIST; an unknown assay stops the pipeline.
                   Output line "chemistry,whitelist" (check_assay_spects.py).
  guide table      rows grouped by Target_name get Target|1, Target|2, ...;
                   pipeline_id = "<Target|n>_sgrna_<chr>:<start>:<end>";
                   guide_features.txt = sgRNA_sequences <tab> pipeline_id.
  guide var        guide_chr / guide_start / guide_end / guide_number /
                   target_elements parsed from the pipeline_id.
  cell QC          per lane: n_genes >= 100; n_counts >= knee[EXPECTED_CELL_
                   NUMBER] (total UMIs sorted decreasing); percent_mito =
                   mito UMIs / UMIs < MITO_EXPECTED_PERCENTAGE (0.2); Scrublet
                   doublets removed (when scrublet is installed); barcodes
                   shared by guide and RNA; batch_number = lane; genes kept
                   in >= int(n_cells * 0.01) cells.
  gene coordinates transcript_start = TSS (start on +, end on -),
                   'NOT_FOUND', 0, 0 when the gene is missing from the GTF.
  binarisation     guide UMI > GUIDE_UMI_LIMIT (5), optionally after summing
                   guides of the same element (--merge); layer 'binarized'.
  covariates       bath_number, percent_mito, log_number_of_detected_genes =
                   log(n_genes), log_total_gene_count, log_total_guide_count.
  guide types      '_TSS' -> POSITIVE_CONTROL, 'random'/'scrambled' ->
                   NEGATIVE_CONTROL, 'chr' -> PUTATIVE_ENHANCER.
  distance         sccrispr-tools calculate_distance: |guide start - gene
                   start| for pairs on the same chromosome (varp 'distance',
                   'same_chrom'); cis = distance < DISTANCE_NEIGHBORS (1 Mb).
  tested pairs     per element, genes on the element's chromosome whose TSS
                   lies within DISTANCE_NEIGHBORS of the element gene's TSS
                   (all genes when IN_TRANS); guides with > 30 assigned cells.
  results          BH over all tests -> adj_pvalue; significant = adj < 0.01;
                   result_guides (guides x genes, NaN = untested; layers
                   adj_pvalue, log_fold_change, significant, z_value);
                   result_elements = Fisher-combined guide p per element,
                   layer sig_not_adj = p < 0.05; mudata_results.h5mu with
                   guides / scRNA / result_guides / result_elements.

Upstream quirks corrected by default (--upstream-compat reproduces them)
  * log_total_gene_count is log(n_genes + 1) upstream (the gene count, not the
    UMI count); the port uses log(n_counts + 1).
  * log_total_guide_count is computed on the binarised matrix upstream; the
    port uses the guide UMIs.
  * Task 2 filters `mdata[n_genes < 6500]` then overwrites the result with
    `mdata[percent_mito < 0.1]`; the port applies both.
  * PerturbLoader tests 10 RANDOM genes when an element name is not a gene
    in the GTF; the port uses the guide coordinates instead.
  * The differential test is SCEPTRE-like but in Python: NB GLM on the
    covariates, score statistic for the guide indicator, conditional
    resampling of the indicator from a logistic propensity (B = 500) and a
    skew-normal fit to the null (sceptre v0.1 fits a skew-t);
    `differential --engine r` writes the upstream per-element inputs and
    run_sceptre_high_moi script and runs them when Rscript is on PATH.

Subcommands
  assay-spec          Task 1: assay name -> kallisto chemistry + whitelist.
  cellranger-inputs   Cell Ranger task: guide table -> feature_ref.csv,
                      library.csv, guide_features.txt + the count command.
  pipeline-config     run_jamboree_v2: perturb.config + nextflow command.
  preprocess          pipeline QC: RNA + guide counts -> raw MuData
                      (guides / scRNA) with coordinates and covariates.
  mudata-qc           Task 2: QC tables/figures, MOI, guide and element
                      coverage, controls, filter summary.
  guide-gene-distance Task 2: guide x gene distance, same-chrom, cis pairs.
  subset              Task 2/3: subset by element / guide category / guide,
                      or keep only some modalities.
  assign-guides       guide task: binarised layer (UMI threshold or
                      covariate-aware Poisson mixture).
  differential        Task 3: differential perturbation -> results MuData.
  tracks              Task 4: BED, bedGraph, links (BEDPE), tracks.ini,
                      arc figures and the pyGenomeTracks commands.
  run                 assign-guides -> differential -> tracks on a MuData.
  selftest            synthetic screen with planted knock-downs.

Output: Docs/CRISPRFGJamboree/<timestamp>_<label>/.  numpy, pandas, scipy,
anndata, h5py required; mudata optional (h5mu read/written through h5py);
matplotlib optional; scanpy optional (HVG table).

Usage:
    igvfagent crispr-fg-jamboree assay-spec --assay 10xv3
    igvfagent crispr-fg-jamboree cellranger-inputs --guide-table df_from_gasperini_tss.xlsx --rna-fastqs R1.fq.gz R2.fq.gz --guide-fastqs G1.fq.gz G2.fq.gz --transcriptome refdata-gex-GRCh38-2020-A
    igvfagent crispr-fg-jamboree pipeline-config --guide-features df_from_gasperini_tss.xlsx --set EXPECTED_CELL_NUMBER=10000
    igvfagent crispr-fg-jamboree preprocess --rna rna.h5ad --guides guides.h5ad --gene-table genes.gtf.gz --expected-cells 5000
    igvfagent crispr-fg-jamboree mudata-qc --mudata raw_mudata_guide_and_transcripts.h5mu
    igvfagent crispr-fg-jamboree assign-guides --mudata raw_mudata_guide_and_transcripts.h5mu --method umi --guide-umi-limit 5
    igvfagent crispr-fg-jamboree differential --mudata mu_with_binary.h5mu --distance 1000000
    igvfagent crispr-fg-jamboree tracks --results mudata_results.h5mu
    igvfagent crispr-fg-jamboree selftest --no-plots
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
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRFGJamboree"

UPSTREAM_REPO = "IGVF-CRISPR/CRISPR_FG_JAMBOREE"
UPSTREAM_COMMIT = "af2af8da297464e5aa2883bb81857342b4e3c369"
PIPELINE_REPO = "LucasSilvaFerreira/pipeline_perturbseq_like"
PIPELINE_COMMIT = "60cd461ae8d1788dfda4a04f686a42b796a79559"
WL_BASE = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/"

# assay -> (kallisto -x string, whitelist file, note)
ASSAYS = {
    "10xv2": ("0,0,16:0,16,26:1,0,0", "737K-august-2016.txt", "10x 3' v2: R1 = 16 bp barcode + 10 bp UMI, R2 = cDNA"),
    "10xv3": ("0,0,16:0,16,28:1,0,0", "3M-february-2018.txt.gz", "10x 3' v3/v3.1: R1 = 16 bp barcode + 12 bp UMI, R2 = cDNA"),
    "10x5p": ("0,0,16:0,16,26:1,0,0", "737K-august-2016.txt", "10x 5' v1/v2, R2-only quantification: R1 = barcode + UMI, R2 = cDNA"),
    "10x5p-pe": ("0,0,16:0,16,26:0,26,0,1,0,0", "737K-august-2016.txt", "10x 5' paired-end: cDNA from R1 after the UMI and all of R2"),
    "dropseq": ("0,0,12:0,12,20:1,0,0", "", "Drop-seq: 12 bp barcode + 8 bp UMI, no whitelist"),
    "celseq2": ("0,6,14:0,0,6:1,0,0", "", "CEL-Seq2: 6 bp UMI then 8 bp barcode"),
}
ASSAY_ALIASES = {"10XV2": "10xv2", "10XV3": "10xv3", "10X5'": "10x5p", "10X5P": "10x5p", "10XV2_5P": "10x5p"}

PIPELINE_PARAMS = [  # run_jamboree_v2.ipynb gasperini_perturb.config, in order
    ("GTF_GZ_LINK", "http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"),
    ("TRANSCRIPTOME_REFERENCE", "human"),
    ("KALLISTO_BIN", "/opt/conda/bin/kallisto"),
    ("GENOME", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"),
    ("GUIDE_FEATURES", "df_from_gasperini_tss.xlsx"),
    ("CHEMISTRY", "10xv2"),
    ("WHITELIST", "737K-august-2016.txt"),
    ("THREADS", 15),
    ("DISTANCE_NEIGHBORS", 1000000),
    ("IN_TRANS", "FALSE"),
    ("FASTQ_FILES_TRANSCRIPTS", []),
    ("FASTQ_NAMES_TRANSCRIPTS", ["S1_L1"]),
    ("FASTQ_FILES_GUIDES", []),
    ("FASTQ_NAMES_GUIDES", ["S1_L1"]),
    ("EXPECTED_CELL_NUMBER", 10000),
    ("MITO_SPECIE", "hsapiens"),
    ("MITO_EXPECTED_PERCENTAGE", 0.2),
    ("PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS", 0.01),
    ("TRANSCRIPTS_UMI_TRHESHOLD", 200),
    ("GUIDE_UMI_LIMIT", 5),
    ("CREATE_REF", False),
]
RESULT_LAYERS = ["adj_pvalue", "log_fold_change", "significant", "z_value"]

INK, INK2, AXIS, SURFACE = "#1f2328", "#57606a", "#8c959f", "#ffffff"
C_SIG, C_NS, C_ELEM = "#c0392b", "#8c959f", "#3b6fb6"


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_fg_jamboree_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
        raise SystemExit("this subcommand needs pandas + numpy + scipy + anndata: pip install 'igvfagent[analysis]'") from e


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
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 40) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "NA"
        v = float(x)
        return f"{v:.{nd}g}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if math.isfinite(f) else None
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(_jsonable(obj), indent=2))
    print(f"JSON: {path}")
    return path


def write_tsv(df, path: Path, header: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, header=header, compression="gzip" if str(path).endswith(".gz") else None)
    print(f"TSV: {path}")
    return path


def write_csv(df, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"CSV: {path}")
    return path


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    print(f"Wrote: {path}")
    return path


def write_report(path: Path, text: str) -> Path:
    path.write_text(text)
    print(f"Report: {path}")
    return path


def dense(X):
    np = _np()
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


# ---------------------------------------------------------------------------
# MuData I/O (mudata optional; h5py + anndata otherwise)
# ---------------------------------------------------------------------------

def read_h5mu(path: Path) -> "dict[str, Any]":
    try:
        import mudata as md  # type: ignore
        m = md.read_h5mu(str(path))
        return {k: m.mod[k].copy() for k in m.mod}
    except ImportError:
        pass
    import h5py  # type: ignore
    try:
        from anndata.io import read_elem  # type: ignore
    except ImportError:
        from anndata.experimental import read_elem  # type: ignore
    with h5py.File(str(path), "r") as h:
        order = list(h["mod"].attrs.get("mod-order", list(h["mod"].keys())))
        return {str(k): read_elem(h[f"mod/{k}"]) for k in order}


def write_h5mu(path: Path, mods: "dict[str, Any]") -> Path:
    """Write a MuData file: mudata when installed, else the same layout via h5py + anndata."""
    pd = _pd()
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import mudata as md  # type: ignore
        md.MuData(mods).write_h5mu(str(path))
        print(f"Wrote: {path}")
        return path
    except ImportError:
        pass
    import h5py  # type: ignore
    try:
        from anndata.io import write_elem  # type: ignore
    except ImportError:
        from anndata.experimental import write_elem  # type: ignore
    obs_names, var_names = [], []
    for a in mods.values():
        obs_names.extend(a.obs_names)
        var_names.extend(a.var_names)
    with h5py.File(str(path), "w") as h:
        h.attrs["encoding-type"] = "MuData"
        h.attrs["encoding-version"] = "0.1.0"
        h.attrs["encoder"] = "igvfagent"
        for k, a in mods.items():
            write_elem(h, f"mod/{k}", a)
        h["mod"].attrs["mod-order"] = list(mods)
        write_elem(h, "obs", pd.DataFrame(index=pd.Index(list(dict.fromkeys(obs_names)), dtype=str)))
        write_elem(h, "var", pd.DataFrame(index=pd.Index(var_names, dtype=str)))
    print(f"Wrote: {path}")
    return path


def pick_mod(mods: dict, candidates: "list[str]") -> str:
    for c in candidates:
        if c in mods:
            return c
    raise SystemExit(f"none of {candidates} in the MuData modalities {list(mods)}")


def guide_key(mods: dict) -> str:
    return pick_mod(mods, ["guides", "guide", "gRNA", "grna"])


def rna_key(mods: dict) -> str:
    return pick_mod(mods, ["scRNA", "rna", "gene", "RNA"])


# ---------------------------------------------------------------------------
# guide naming / annotation
# ---------------------------------------------------------------------------

def read_guide_table(path: Path):
    pd = _pd()
    s = str(path).lower()
    if s.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep="\t" if s.endswith((".tsv", ".txt", ".tsv.gz")) else ",")
    need = {"Target_name", "chr", "start", "end", "sgRNA_sequences"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"guide table lacks columns {sorted(miss)} (needs Target_name, chr, start, end, sgRNA_sequences)")
    return df


def name_guides(df):
    """guide_table_processing.py: Target|n numbering within each target and the pipeline_id."""
    pd = _pd()
    out = []
    for k, v in df.groupby("Target_name", sort=True):
        v = v.copy()
        v["Target_name"] = [f"{k}|{n + 1}" for n in range(v.shape[0])]
        v["pipeline_id"] = [f"{t}_sgrna_{c}:{s}:{e}" for t, c, s, e in zip(v["Target_name"], v["chr"], v["start"], v["end"])]
        out.append(v)
    return pd.concat(out) if out else df.assign(pipeline_id=[])


def parse_guide_id(k: str) -> dict:
    """muon_creation.py parsing of '<Target>|<n>_sgrna_<chr>:<start>:<end>'."""
    try:
        right = k.split("|")[1]
        head = right.split(":")[0]
        return {"guide_chr": head.split("_")[-1], "guide_start": right.split(":")[-2], "guide_end": right.split(":")[-1],
                "guide_number": head.split("_")[0], "target_elements": k.split("|")[0]}
    except IndexError:
        return {"guide_chr": "NOT_FOUND", "guide_start": "0", "guide_end": "0", "guide_number": "1", "target_elements": k.split("|")[0]}


def guide_type(element: str) -> "Optional[str]":
    """PerturbLoader add_test_type."""
    if "_TSS" in element:
        return "POSITIVE_CONTROL"
    if "random" in element or "scrambled" in element:
        return "NEGATIVE_CONTROL"
    if "chr" in element:
        return "PUTATIVE_ENHANCER"
    return None


def annotate_guides(var):
    pd = _pd()
    v = var.copy()
    parsed = pd.DataFrame([parse_guide_id(str(k)) for k in v.index], index=v.index)
    for c in parsed.columns:
        if c not in v.columns:
            v[c] = parsed[c]
    if "feature_name" not in v.columns:
        v["feature_name"] = v.index.astype(str)
    v["guide_type"] = [guide_type(str(e)) or "UNCLASSIFIED" for e in v["target_elements"]]
    return v


# ---------------------------------------------------------------------------
# gene table (GTF or TSV) -> TSS coordinates
# ---------------------------------------------------------------------------

def read_gene_table(path: Path):
    """Genes with chr, start, end, strand, gene_id (no version), gene_name, tss."""
    pd = _pd()
    s = str(path)
    op = gzip.open if s.endswith(".gz") else open
    if ".gtf" in s:
        rows = []
        with op(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
                rows.append({"chr": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                             "gene_id": attrs.get("gene_id", "").split(".")[0], "gene_name": attrs.get("gene_name", attrs.get("gene_id", ""))})
        g = pd.DataFrame(rows)
    else:
        g = pd.read_csv(path, sep="\t")
        g = g.rename(columns={"chrom": "chr", "chromosome": "chr", "gene": "gene_name", "symbol": "gene_name"})
        if "gene_id" not in g.columns:
            g["gene_id"] = g["gene_name"]
        if "strand" not in g.columns:
            g["strand"] = "+"
    g["chr"] = g["chr"].astype(str)
    g["tss"] = [s0 if st == "+" else e0 for s0, e0, st in zip(g["start"], g["end"], g["strand"])]
    return g


def add_gene_coords(var, genes):
    """muon_creation.py: transcript_chr / transcript_start (TSS) / transcript_end, NOT_FOUND 0 0 if absent."""
    v = var.copy()
    by_id = {str(i): r for i, r in zip(genes["gene_id"].astype(str), genes.itertuples())}
    by_name = {str(n): r for n, r in zip(genes["gene_name"].astype(str), genes.itertuples())}
    names = v["feature_name"].astype(str) if "feature_name" in v.columns else v.index.astype(str)
    chrs, starts, ends = [], [], []
    for idx, nm in zip(v.index.astype(str), names):
        r = by_id.get(idx.split(".")[0]) or by_name.get(nm) or by_name.get(idx)
        if r is None:
            chrs.append("NOT_FOUND"); starts.append(0); ends.append(0)
        else:
            chrs.append(str(r.chr).replace("chr", "")); starts.append(int(r.tss)); ends.append(int(r.end))
    v["transcript_chr"], v["transcript_start"], v["transcript_end"] = chrs, starts, ends
    return v


# ---------------------------------------------------------------------------
# Task 1: assay spec
# ---------------------------------------------------------------------------

def resolve_assay(assay: str, chemistry: str = "", whitelist: str = "") -> dict:
    if assay == "custom":
        return {"assay": "custom", "chemistry": chemistry or "custom", "whitelist": whitelist or "custom", "whitelist_url": None,
                "note": "custom: CHEMISTRY and WHITELIST parameters are used as given"}
    key = ASSAY_ALIASES.get(assay.upper(), assay.lower())
    if key not in ASSAYS:
        raise KeyError(assay)
    chem, wl, note = ASSAYS[key]
    return {"assay": key, "chemistry": chem, "whitelist": wl, "whitelist_url": (WL_BASE + wl) if wl else None, "note": note}


def cmd_assay_spec(args: argparse.Namespace) -> int:
    try:
        spec = resolve_assay(args.assay, args.chemistry or "", args.whitelist or "")
    except KeyError:
        print(f"ERROR: assay '{args.assay}' is not in the table ({', '.join(sorted(ASSAYS))}, custom); the pipeline stops here.")
        return 2
    rd = run_dir(args.label)
    if args.download and spec["whitelist_url"]:
        dest = rd / spec["whitelist"]
        try:
            urllib.request.urlretrieve(spec["whitelist_url"], str(dest))
            spec["whitelist_path"] = str(dest)
            print(f"Wrote: {dest}")
        except Exception as e:  # network optional
            spec["whitelist_download_error"] = str(e)
            print(f"whitelist download failed ({e}); URL: {spec['whitelist_url']}")
    wl = spec.get("whitelist_path") or spec["whitelist"] or "none"
    line = f"{spec['chemistry']},{wl}"
    print(line)
    table = [(k, v[0], v[1] or "-", v[2]) for k, v in ASSAYS.items()]
    write_text(rd / "assay_spec.csv", line + "\n")
    write_json({**spec, "csv_line": line, "table": {k: {"chemistry": v[0], "whitelist": v[1], "note": v[2]} for k, v in ASSAYS.items()}},
               rd / "summary.json")
    write_report(rd / "report.md", f"# Assay specification — {spec['assay']}\n\n`kallisto bus -x {spec['chemistry']}` with whitelist "
                 f"`{wl}`.\n\n{spec['note']}\n\nOutput line (check_assay_spects.py contract): `{line}`\n\n"
                 + md_table(["assay", "kallisto -x", "whitelist", "layout"], table) + "\n")
    return 0


# ---------------------------------------------------------------------------
# Cell Ranger task
# ---------------------------------------------------------------------------

FASTQ_RE = re.compile(r"_S\d+_L\d{3}_[RI]\d_\d{3}\.f(ast)?q(\.gz)?$")


def fastq_libraries(paths: "list[str]", library_type: str) -> "list[dict]":
    seen, out = set(), []
    for p in paths:
        for q in str(p).split():
            path = Path(q)
            sample = FASTQ_RE.sub("", path.name)
            key = (str(path.parent), sample)
            if key not in seen:
                seen.add(key)
                out.append({"fastqs": str(path.parent), "sample": sample, "library_type": library_type})
    return out


def cmd_cellranger_inputs(args: argparse.Namespace) -> int:
    pd = _pd()
    rd = run_dir(args.label)
    gt = name_guides(read_guide_table(Path(args.guide_table)))
    write_tsv(gt[["sgRNA_sequences", "pipeline_id"]], rd / "guide_features.txt", header=False)
    fr = pd.DataFrame({"id": gt["pipeline_id"], "name": gt["pipeline_id"], "read": args.read, "pattern": args.pattern,
                       "sequence": gt["sgRNA_sequences"], "feature_type": "CRISPR Guide Capture"})
    if "target_gene_id" in gt.columns:
        fr["target_gene_id"] = gt["target_gene_id"]
        fr["target_gene_name"] = gt.get("target_gene_name", gt["target_gene_id"])
    write_csv(fr, rd / "feature_ref.csv")
    libs = fastq_libraries(args.rna_fastqs or [], "Gene Expression") + fastq_libraries(args.guide_fastqs or [], "CRISPR Guide Capture")
    if libs:
        write_csv(pd.DataFrame(libs, columns=["fastqs", "sample", "library_type"]), rd / "library.csv")
    cmd = ["cellranger", "count", f"--id={args.id}", f"--libraries={rd / 'library.csv'}", f"--transcriptome={args.transcriptome}",
           f"--feature-ref={rd / 'feature_ref.csv'}", f"--localmem={args.localmem}", f"--localcores={args.threads}"]
    if args.chemistry:
        cmd.append(f"--chemistry={args.chemistry}")
    ran = False
    if args.run and libs:
        if shutil.which("cellranger"):
            rc = subprocess.call(cmd, cwd=str(rd))
            ran = rc == 0
        else:
            print("cellranger not on PATH; the command would be:")
    print("  " + " ".join(cmd))
    print(f"Outputs: {rd / args.id / 'outs' / 'filtered_feature_bc_matrix'} (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)")
    write_json({"n_guides": len(gt), "n_targets": int(gt["Target_name"].str.split("|").str[0].nunique()), "libraries": libs,
                "command": cmd, "ran": ran}, rd / "summary.json")
    write_report(rd / "report.md", f"# Cell Ranger inputs\n\n{len(gt)} guides; feature_ref.csv (read {args.read}, pattern "
                 f"`{args.pattern}`), library.csv ({len(libs)} libraries), guide_features.txt.\n\n```bash\n{' '.join(cmd)}\n```\n\n"
                 + md_table(["id", "sequence"], zip(fr["id"], fr["sequence"]), 20) + "\n")
    return 0


# ---------------------------------------------------------------------------
# run_jamboree_v2: pipeline config
# ---------------------------------------------------------------------------

def _nf_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, list):
        return "[" + ", ".join(f"'{x}'" for x in v) + "]"
    return f"'{v}'" if not str(v).startswith('"') else str(v)


def _parse_set(s: str, defaults: "Optional[dict]" = None):
    k, _, v = s.partition("=")
    v = v.strip()
    if isinstance((defaults or {}).get(k.strip()), bool) and v.lower() in ("true", "false"):
        return k.strip(), v.lower() == "true"
    if v.startswith("["):
        return k.strip(), [x.strip().strip("'\"") for x in v.strip("[]").split(",") if x.strip()]
    try:
        return k.strip(), int(v)
    except ValueError:
        try:
            return k.strip(), float(v)
        except ValueError:
            return k.strip(), v


def cmd_pipeline_config(args: argparse.Namespace) -> int:
    params = dict(PIPELINE_PARAMS)
    if args.guide_features:
        params["GUIDE_FEATURES"] = args.guide_features
    if args.rna_fastqs:
        params["FASTQ_FILES_TRANSCRIPTS"] = [" ".join(args.rna_fastqs)]
    if args.guide_fastqs:
        params["FASTQ_FILES_GUIDES"] = [" ".join(args.guide_fastqs)]
    for s in args.set or []:
        k, v = _parse_set(s, dict(PIPELINE_PARAMS))
        params[k] = v
    rd = run_dir(args.label)
    cfg = rd / "perturb.config"
    write_text(cfg, "".join(f"params.{k} = {_nf_value(v)}\n" for k, v in params.items()))
    out_dir = args.work_name
    env = [f"export NXF_VER={args.nxf_ver}"]
    tower = []
    if args.with_tower:
        if os.environ.get("TOWER_ACCESS_TOKEN"):
            tower = ["-with-tower"]
        else:
            print("--with-tower: TOWER_ACCESS_TOKEN is not set in the environment; launching without Tower (the token is never written to disk)")
    cmd = ["nextflow", "run", args.main_nf, "-c", str(cfg)] + tower + ["-with-timeline", str(rd / f"{out_dir}.html")] + \
          (["-resume"] if args.resume else []) + ["-w", out_dir]
    ran = False
    if args.run:
        if shutil.which("nextflow"):
            e = dict(os.environ, NXF_VER=args.nxf_ver)
            ran = subprocess.call(cmd, cwd=str(rd), env=e) == 0
        else:
            print("nextflow not on PATH; the launch would be:")
    print("  " + "; ".join(env) + "; " + " ".join(cmd))
    write_json({"params": params, "command": cmd, "env": env, "ran": ran, "pipeline": f"{PIPELINE_REPO}@{PIPELINE_COMMIT}"}, rd / "summary.json")
    write_report(rd / "report.md", f"# scPerturb-seq pipeline launch\n\nPipeline `{PIPELINE_REPO}` (main.nf).  Config `{cfg.name}`:\n\n"
                 + md_table(["param", "value"], ((k, _nf_value(v)) for k, v in params.items()), 60)
                 + f"\n\n```bash\n{'; '.join(env)}\n{' '.join(cmd)}\n```\n\nTower tokens are read from TOWER_ACCESS_TOKEN only.\n")
    return 0


# ---------------------------------------------------------------------------
# preprocess (filtering_and_lane_merging.py + muon_creation.py)
# ---------------------------------------------------------------------------

def read_counts(path: str):
    """AnnData from .h5ad or a 10x matrix directory / .mtx(.gz) (features x cells)."""
    import anndata as ad  # type: ignore
    pd = _pd()
    p = Path(path)
    if p.suffix == ".h5ad":
        return ad.read_h5ad(str(p))
    import scipy.io  # type: ignore
    d = p if p.is_dir() else p.parent
    mtx = next((d / n for n in ("matrix.mtx.gz", "matrix.mtx") if (d / n).is_file()), p)
    op = gzip.open if str(mtx).endswith(".gz") else open
    with op(mtx, "rb") as fh:
        M = scipy.io.mmread(fh).tocsr()
    feat = next((d / n for n in ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv") if (d / n).is_file()), None)
    bc = next((d / n for n in ("barcodes.tsv.gz", "barcodes.tsv") if (d / n).is_file()), None)
    fdf = pd.read_csv(feat, sep="\t", header=None) if feat else pd.DataFrame({0: [f"f{i}" for i in range(M.shape[0])]})
    bcs = pd.read_csv(bc, sep="\t", header=None)[0].astype(str).tolist() if bc else [f"c{i}" for i in range(M.shape[1])]
    a = ad.AnnData(M.T.tocsr().astype("float32"))
    a.obs_names = bcs
    a.var_names = fdf[0].astype(str).tolist()
    a.var["feature_name"] = fdf[1].astype(str).tolist() if fdf.shape[1] > 1 else a.var_names
    a.var_names_make_unique()
    return a


def knee_threshold(total_counts, expected_cells: int) -> float:
    np = _np()
    knee = np.sort(np.asarray(total_counts, dtype=float))[::-1]
    return float(knee[min(expected_cells, len(knee) - 1)])


def qc_lane(rna, guide, lane, expected_cells: int, mito_max: float, min_genes: int, mito_genes=None, doublets: bool = True):
    """analyze_batch(): per-lane RNA QC, guide nonzero counts, shared barcodes."""
    np = _np()
    stats = {"lane": lane, "barcodes_in": rna.n_obs}
    X = rna.X
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = n_genes >= min_genes
    rna = rna[keep].copy()
    stats["after_min_genes"] = rna.n_obs
    tot = np.asarray(rna.X.sum(axis=1)).ravel()
    thr = knee_threshold(tot, expected_cells)
    rna = rna[tot >= thr].copy()
    stats["knee_min_counts"] = thr
    stats["after_knee"] = rna.n_obs
    names = rna.var["feature_name"].astype(str) if "feature_name" in rna.var.columns else rna.var_names.astype(str)
    if mito_genes:
        is_mt = np.array([(g in mito_genes) or (n in mito_genes) for g, n in zip(rna.var_names.astype(str).str.split(".").str[0], names)])
    else:
        is_mt = np.array([n.upper().startswith("MT-") for n in names])
    tot = np.asarray(rna.X.sum(axis=1)).ravel()
    rna.obs["n_genes"] = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    rna.obs["n_counts"] = tot
    rna.obs["percent_mito"] = np.asarray(rna.X[:, is_mt].sum(axis=1)).ravel() / np.maximum(tot, 1) if is_mt.any() else 0.0
    rna = rna[rna.obs["percent_mito"] < mito_max].copy()
    stats["after_mito"] = rna.n_obs
    stats["n_mito_genes"] = int(is_mt.sum())
    stats["doublets"] = "skipped"
    if doublets:
        try:
            import scrublet as scr  # type: ignore
            s = scr.Scrublet(rna.X)
            rna.obs["doublet_scores"], rna.obs["predicted_doublets"] = s.scrub_doublets()
            rna.obs["doublet_info"] = rna.obs["predicted_doublets"].astype(str)
            rna = rna[~rna.obs["predicted_doublets"].astype(bool)].copy()
            stats["doublets"] = "scrublet"
        except ImportError:
            stats["doublets"] = "skipped (scrublet not installed)"
    stats["after_doublets"] = rna.n_obs
    guide = guide.copy()
    guide.obs["number_of_nonzero_guides"] = np.asarray((guide.X > 0).sum(axis=1)).ravel()
    shared = [c for c in rna.obs_names if c in set(guide.obs_names)]
    rna, guide = rna[shared].copy(), guide[shared].copy()
    rna.obs["batch_number"] = lane
    guide.obs["batch_number"] = lane
    stats["after_barcode_intersection"] = len(shared)
    return rna, guide, stats


def make_covariates(rna, guide_counts, binarized, compat: bool = False):
    """generating_covariates(); corrected log_total_gene_count / log_total_guide_count unless compat."""
    np, pd = _np(), _pd()
    c = pd.DataFrame(index=rna.obs_names)
    c["bath_number"] = rna.obs["batch_number"].values
    c["percent_mito"] = rna.obs["percent_mito"].values
    c["log_number_of_detected_genes"] = np.log(rna.obs["n_genes"].values.astype(float))
    if compat:
        c["log_total_gene_count"] = np.log(rna.obs["n_genes"].values.astype(float) + 1)
        c["log_total_guide_count"] = np.log(np.asarray(binarized.sum(axis=1)).ravel() + 1)
    else:
        c["log_total_gene_count"] = np.log(rna.obs["n_counts"].values.astype(float) + 1)
        c["log_total_guide_count"] = np.log(np.asarray(guide_counts.sum(axis=1)).ravel() + 1)
    return c


def cmd_preprocess(args: argparse.Namespace) -> int:
    import anndata as ad  # type: ignore
    np, pd = _np(), _pd()
    rnas, guides = args.rna, args.guides
    if len(rnas) != len(guides):
        raise SystemExit("--rna and --guides need one entry per lane, in the same order")
    lanes = args.lanes or [str(i + 1) for i in range(len(rnas))]
    mito = set(Path(args.mito_genes).read_text().split()) if args.mito_genes else None
    rd = run_dir(args.label)
    lane_rna, lane_guide, lane_stats = [], [], []
    for r, g, lane in zip(rnas, guides, lanes):
        R, G = read_counts(r), read_counts(g)
        R2, G2, st = qc_lane(R, G, lane, args.expected_cells, args.mito_max, args.min_genes, mito, not args.no_doublets)
        lane_rna.append(R2); lane_guide.append(G2); lane_stats.append(st)
        print(f"[lane {lane}] {st['barcodes_in']} barcodes -> {st['after_barcode_intersection']} cells")
    rna = ad.concat(lane_rna, merge="same") if len(lane_rna) > 1 else lane_rna[0]
    guide = ad.concat(lane_guide, merge="same") if len(lane_guide) > 1 else lane_guide[0]
    rna.obs_names_make_unique(); guide.obs_names_make_unique()
    min_cells = int(rna.n_obs * args.pct_cells_gene)
    ncell = np.asarray((rna.X > 0).sum(axis=0)).ravel()
    rna = rna[:, ncell >= min_cells].copy()
    rna.var["n_cells"] = np.asarray((rna.X > 0).sum(axis=0)).ravel()
    if "feature_name" not in rna.var.columns:
        rna.var["feature_name"] = rna.var_names.astype(str)
    if args.gene_table:
        rna.var = add_gene_coords(rna.var, read_gene_table(Path(args.gene_table)))
    guide.var = annotate_guides(guide.var)
    binar = (dense(guide.X) > args.guide_umi_limit).astype(float)
    cov = make_covariates(rna, guide.X, binar, compat=args.upstream_compat)
    for c in cov.columns:
        rna.obs[c] = cov[c].values
        guide.obs[c] = cov[c].values
    out = write_h5mu(rd / "raw_mudata_guide_and_transcripts.h5mu", {"guides": guide, "scRNA": rna})
    write_tsv(pd.DataFrame(lane_stats), rd / "qc_per_lane.tsv")
    figs = []
    plt = _plt()
    if plt is not None and not args.no_plots:
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
        tot = np.sort(np.asarray(read_counts(rnas[0]).X.sum(axis=1)).ravel())[::-1]
        axes[0].loglog(tot, np.arange(1, len(tot) + 1), color=C_ELEM, lw=2)
        axes[0].axvline(lane_stats[0]["knee_min_counts"], color=INK, lw=1)
        _style(axes[0], "Knee plot (lane 1)", "UMI counts", "barcodes")
        axes[1].scatter(rna.obs["n_counts"], rna.obs["n_genes"], s=4, alpha=0.4, color=C_ELEM)
        _style(axes[1], "Saturation", "UMI counts", "genes detected")
        axes[2].hist(guide.obs["number_of_nonzero_guides"], bins=30, color=C_ELEM)
        _style(axes[2], "Guides per cell (MOI)", "nonzero guides", "cells")
        figs.append(str(_save(fig, rd / "qc_preprocess.png")))
    summ = {"lanes": lane_stats, "n_cells": rna.n_obs, "n_genes": rna.n_vars, "n_guides": guide.n_vars,
            "min_cells_per_gene": min_cells, "upstream_compat": args.upstream_compat, "mudata": str(out), "figures": figs}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", "# Preprocess (pipeline QC)\n\n" + md_table(list(lane_stats[0]), [list(s.values()) for s in lane_stats])
                 + f"\n\n{rna.n_obs} cells x {rna.n_vars} genes (kept in >= {min_cells} cells) and {guide.n_vars} guides -> `{out.name}`.\n")
    return 0


# ---------------------------------------------------------------------------
# Task 2: MuData QC, distances, subsetting
# ---------------------------------------------------------------------------

def guide_gene_distance(gvar, rvar, how: str = "start"):
    """calculate_same_chrom + calculate_distance (sccrispr-tools): |start - start| on the same chromosome."""
    np, pd = _np(), _pd()
    gchr = gvar["guide_chr"].astype(str).str.replace("chr", "", regex=False).to_numpy()
    rchr = rvar["transcript_chr"].astype(str).str.replace("chr", "", regex=False).to_numpy()
    gs = pd.to_numeric(gvar["guide_start"], errors="coerce").to_numpy(dtype=float)
    rs = pd.to_numeric(rvar["transcript_start"], errors="coerce").to_numpy(dtype=float)
    if how == "midpoint":
        gs = (gs + pd.to_numeric(gvar["guide_end"], errors="coerce").to_numpy(dtype=float)) / 2
        rs = (rs + pd.to_numeric(rvar["transcript_end"], errors="coerce").to_numpy(dtype=float)) / 2
    same = (gchr[:, None] == rchr[None, :]) & (rchr[None, :] != "NOT_FOUND")
    dist = np.abs(gs[:, None] - rs[None, :])
    return same, np.where(same, dist, np.nan)


def cmd_guide_gene_distance(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    mods = read_h5mu(Path(args.mudata))
    g, r = mods[guide_key(mods)], mods[rna_key(mods)]
    gv = annotate_guides(g.var)
    same, dist = guide_gene_distance(gv, r.var, args.how)
    ii, jj = np.where(same)
    df = pd.DataFrame({"guide": g.var_names[ii], "target_elements": gv["target_elements"].to_numpy()[ii],
                       "gene": r.var_names[jj], "distance": dist[ii, jj]})
    df["cis"] = df["distance"] < args.distance
    rd = run_dir(args.label)
    write_tsv(df.sort_values(["guide", "distance"]), rd / "guide_gene_distance.tsv.gz")
    cis = df[df["cis"]]
    write_tsv(cis, rd / "cis_pairs.tsv")
    summ = {"n_guides": g.n_vars, "n_genes": r.n_vars, "same_chrom_pairs": len(df), "cis_pairs": len(cis), "distance": args.distance, "how": args.how}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", f"# Guide-gene distances\n\n{len(df)} same-chromosome guide-gene pairs; {len(cis)} within "
                 f"{args.distance:,} bp (cis).\n\n" + md_table(["guide", "gene", "distance"], cis[["guide", "gene", "distance"]].head(30).values) + "\n")
    return 0


def cmd_mudata_qc(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    mods = read_h5mu(Path(args.mudata))
    gk, rk = guide_key(mods), rna_key(mods)
    g, r = mods[gk], mods[rk]
    rd = run_dir(args.label)
    X = r.X
    tot = np.asarray(X.sum(axis=1)).ravel()
    ngen = np.asarray((X > 0).sum(axis=1)).ravel()
    pm = r.obs["percent_mito"].to_numpy(dtype=float) if "percent_mito" in r.obs.columns else np.zeros(r.n_obs)
    cells = pd.DataFrame({"cell": r.obs_names, "n_genes": ngen, "n_counts": tot, "percent_mito": pm})
    write_tsv(cells, rd / "cell_qc.tsv.gz")
    frac = np.asarray((X.multiply(1 / np.maximum(tot, 1)[:, None]) if hasattr(X, "multiply") else X / np.maximum(tot, 1)[:, None]).mean(axis=0)).ravel()
    names = r.var["feature_name"].astype(str).to_numpy() if "feature_name" in r.var.columns else r.var_names.to_numpy()
    top = pd.DataFrame({"gene": names, "mean_fraction_of_counts": frac}).sort_values("mean_fraction_of_counts", ascending=False).head(args.n_top)
    write_tsv(top, rd / "highest_expressed_genes.tsv")
    if args.upstream_compat:
        keep = pm < args.max_mito
    else:
        keep = (ngen < args.max_genes) & (pm < args.max_mito)
    gv = annotate_guides(g.var)
    G = g.X
    nz_guides = np.asarray((G > 0).sum(axis=1)).ravel()
    moi = pd.Series(nz_guides).value_counts().sort_index()
    write_tsv(pd.DataFrame({"number_of_nonzero_guides": moi.index, "cells": moi.values}), rd / "moi_histogram.tsv")
    gcov = pd.DataFrame({"guide": g.var_names, "target_elements": gv["target_elements"].values, "guide_type": gv["guide_type"].values,
                         "n_cells": np.asarray((G > 0).sum(axis=0)).ravel()})
    write_tsv(gcov, rd / "guide_coverage.tsv")
    ecov = gcov.groupby("target_elements", as_index=False).agg(n_guides=("guide", "size"), n_cells=("n_cells", "sum"),
                                                               guide_type=("guide_type", "first"))
    write_tsv(ecov, rd / "element_coverage.tsv")
    ctrl = gcov["guide_type"].fillna("UNCLASSIFIED").value_counts().to_dict()
    hvg = None
    try:
        import scanpy as sc  # type: ignore
        rr = r[keep].copy()
        sc.pp.normalize_total(rr, target_sum=1e4)
        sc.pp.log1p(rr)
        sc.pp.highly_variable_genes(rr)
        hv = rr.var[["highly_variable", "means", "dispersions_norm"]].copy()
        hv.insert(0, "gene", rr.var_names)
        write_tsv(hv, rd / "highly_variable_genes.tsv.gz")
        hvg = int(hv["highly_variable"].sum())
    except Exception as e:  # scanpy optional
        logging.info("HVG skipped: %s", e)
    figs = []
    plt = _plt()
    if plt is not None and not args.no_plots:
        fig, axes = plt.subplots(2, 3, figsize=(12, 6.4))
        for ax, (col, lab) in zip(axes[0], [("n_genes", "genes"), ("n_counts", "UMIs"), ("percent_mito", "fraction mito")]):
            ax.hist(cells[col], bins=40, color=C_ELEM)
            _style(ax, lab, lab, "cells")
        axes[1, 0].scatter(cells["n_counts"], cells["percent_mito"], s=4, alpha=0.4, color=C_ELEM)
        _style(axes[1, 0], "UMIs vs mito", "UMIs", "fraction mito")
        axes[1, 1].bar(moi.index, moi.values, color=C_ELEM)
        _style(axes[1, 1], "MOI", "nonzero guides per cell", "cells")
        axes[1, 2].hist(ecov["n_cells"], bins=20, color=C_ELEM)
        _style(axes[1, 2], "Element coverage", "cells per element", "elements")
        fig.tight_layout()
        figs.append(str(_save(fig, rd / "mudata_qc.png")))
    summ = {"modalities": {k: list(v.shape) for k, v in mods.items()}, "cells_passing_filter": int(keep.sum()),
            "filter": {"max_genes": None if args.upstream_compat else args.max_genes, "max_mito": args.max_mito},
            "median_nonzero_guides": float(np.median(nz_guides)), "guide_types": ctrl, "n_elements": len(ecov),
            "n_highly_variable": hvg, "figures": figs}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", f"# MuData QC\n\nModalities: {', '.join(f'{k} {v.n_obs} x {v.n_vars}' for k, v in mods.items())}.\n\n"
                 f"{int(keep.sum())} of {r.n_obs} cells pass n_genes < {args.max_genes} and percent_mito < {args.max_mito}"
                 f"{' (upstream-compat: mito only)' if args.upstream_compat else ''}.  Median nonzero guides per cell "
                 f"{np.median(nz_guides):.0f}.  Guide categories: {ctrl}.\n\nTop genes:\n\n"
                 + md_table(["gene", "mean fraction"], [(a, _fmt(b, 4)) for a, b in top.values], 20)
                 + "\n\nElement coverage:\n\n" + md_table(["element", "guides", "cells", "type"], ecov[["target_elements", "n_guides", "n_cells", "guide_type"]].values, 30) + "\n")
    return 0


def cmd_subset(args: argparse.Namespace) -> int:
    np = _np()
    mods = read_h5mu(Path(args.mudata))
    gk = guide_key(mods)
    keep_mods = args.keep_modalities or list(mods)
    g = mods[gk]
    gv = annotate_guides(g.var)
    mask = np.ones(g.n_vars, dtype=bool)
    if args.element:
        mask &= gv["target_elements"].isin(args.element).to_numpy()
    if args.category:
        mask &= gv["guide_type"].isin(args.category).to_numpy()
    if args.guide:
        mask &= np.isin(g.var_names, args.guide)
    out_mods = {}
    for k in keep_mods:
        a = mods[k]
        out_mods[k] = a[:, mask].copy() if k == gk else a
    rd = run_dir(args.label)
    out = write_h5mu(rd / (args.out_name or "subset.h5mu"), out_mods)
    info = {"guides_kept": int(mask.sum()), "guides_total": g.n_vars, "modalities": keep_mods, "mudata": str(out)}
    if args.guide and len(args.guide) == 1 and mask.any():
        # guide presence across cells + the counts of its tested (cis) genes
        rk = rna_key(mods)
        r = mods[rk]
        j = int(np.where(mask)[0][0])
        pres = dense(g.X[:, j]).ravel() > 0
        same, dist = guide_gene_distance(gv.iloc[[j]], r.var)
        genes = r.var_names[(dist[0] < args.distance) & same[0]]
        pd = _pd()
        tab = pd.DataFrame(dense(r[:, genes].X), columns=genes, index=r.obs_names)
        tab.insert(0, "guide_present", pres)
        tab.insert(0, "cell", r.obs_names)
        write_tsv(tab, rd / "guide_cells_and_cis_gene_counts.tsv.gz")
        info.update({"cells_with_guide": int(pres.sum()), "cis_genes": list(genes)})
    write_json(info, rd / "summary.json")
    write_report(rd / "report.md", f"# Subset\n\n{info['guides_kept']} of {info['guides_total']} guides; modalities {keep_mods} -> `{out.name}`.\n")
    return 0


# ---------------------------------------------------------------------------
# guide assignment
# ---------------------------------------------------------------------------

def poisson_glm(y, X, offset=None, max_iter: int = 50):
    np = _np()
    n, p = X.shape
    off = np.zeros(n) if offset is None else offset
    beta = np.zeros(p)
    beta[0] = math.log(max(float(np.mean(y)), 1e-8))
    for _ in range(max_iter):
        eta = X @ beta + off
        mu = np.exp(np.clip(eta, -30, 30))
        z = eta - off + (y - mu) / mu
        XtW = X.T * mu
        new = np.linalg.solve(XtW @ X + 1e-8 * np.eye(p), XtW @ z)
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    return beta, X @ beta + off


def mixture_em(y, eta, n_rep: int = 5, rng=None):
    """Two-component Poisson mixture with fixed offsets (reduced EM), best of n_rep random starts."""
    np = _np()
    from scipy.special import gammaln  # type: ignore
    rng = rng or np.random.default_rng(0)
    base = np.exp(np.clip(eta, -30, 30))
    lg = gammaln(y + 1)
    best = None
    for _ in range(n_rep):
        pi, g = rng.uniform(1e-5, 0.1), rng.uniform(math.log(10), math.log(5000))
        prev = -math.inf
        for _ in range(300):
            l0 = math.log1p(-pi) + y * eta - base - lg
            l1 = math.log(pi) + y * (eta + g) - base * math.exp(g) - lg
            lz = np.logaddexp(l0, l1)
            T = np.exp(l1 - lz)
            ll = float(lz.sum())
            pi = min(max(float(T.mean()), 1e-10), 1 - 1e-10)
            num, den = float((T * y).sum()), float((T * base).sum())
            if num <= 0 or den <= 0:
                break
            g = math.log(num / den)
            if abs(ll - prev) < 1e-7 * (1 + abs(ll)):
                break
            prev = ll
        l0 = math.log1p(-pi) + y * eta - base - lg
        l1 = math.log(pi) + y * (eta + g) - base * math.exp(g) - lg
        lz = np.logaddexp(l0, l1)
        res = {"pi": pi, "g": g, "ll": float(lz.sum()), "post": np.exp(l1 - lz)}
        if best is None or res["ll"] > best["ll"]:
            best = res
    return best


def assign_guides(G, obs, method: str = "umi", limit: int = 5, threshold: float = 0.8, seed: int = 0):
    """cells x guides counts -> binary matrix (+ per-guide table)."""
    np, pd = _np(), _pd()
    C = dense(G).astype(float)
    if method == "umi":
        B = (C > limit).astype(float)
        return B, pd.DataFrame({"method": "umi", "n_assigned": B.sum(axis=0)})
    rng = np.random.default_rng(seed)
    cols = [np.ones(C.shape[0])]
    for c in ("log_total_guide_count", "log_number_of_detected_genes", "log_total_gene_count"):
        if c in obs.columns:
            v = obs[c].to_numpy(dtype=float)
            if np.std(v) > 0:
                cols.append((v - v.mean()) / v.std())
    X = np.column_stack(cols)
    B = np.zeros_like(C)
    rows = []
    for j in range(C.shape[1]):
        y = C[:, j]
        if (y > 0).sum() < 10:
            B[:, j] = y > limit
            rows.append({"method": "backup_umi", "pi": np.nan, "g": np.nan})
            continue
        _, eta = poisson_glm(y, X)
        fit = mixture_em(y, eta, rng=rng)
        if fit["g"] <= 0:
            B[:, j] = y > limit
            rows.append({"method": "backup_umi", "pi": fit["pi"], "g": fit["g"]})
        else:
            B[:, j] = fit["post"] >= threshold
            rows.append({"method": "poisson_mixture", "pi": fit["pi"], "g": fit["g"]})
    t = pd.DataFrame(rows)
    t["n_assigned"] = B.sum(axis=0)
    return B, t


def merge_guides_by_element(G, var):
    """merge_data(MERGE=True): sum guides targeting the same element (Target|n -> Target)."""
    np, pd = _np(), _pd()
    el = [str(k).split("|")[0] for k in var.index]
    uniq = list(dict.fromkeys(el))
    M = np.zeros((len(el), len(uniq)))
    for i, e in enumerate(el):
        M[i, uniq.index(e)] = 1
    return dense(G) @ M, uniq


def cmd_assign_guides(args: argparse.Namespace) -> int:
    import anndata as ad  # type: ignore
    np, pd = _np(), _pd()
    sp = _sp()
    mods = read_h5mu(Path(args.mudata))
    gk = guide_key(mods)
    g = mods[gk]
    rd = run_dir(args.label)
    if args.merge:
        C, els = merge_guides_by_element(g.X, g.var)
        var = pd.DataFrame(index=pd.Index(els, dtype=str))
        var["feature_name"] = els
        var["target_elements"] = els
        g = ad.AnnData(sp.csr_matrix(C), obs=g.obs.copy(), var=var)
    B, tab = assign_guides(g.X, g.obs, args.method, args.guide_umi_limit, args.probability_threshold, args.seed)
    g.layers[args.layer] = sp.csr_matrix(B)
    g.obs["number_of_assigned_guides"] = B.sum(axis=1)
    mods[gk] = g
    out = write_h5mu(rd / "mu_with_binary.h5mu", mods)
    tab.insert(0, "guide", g.var_names)
    write_tsv(tab, rd / "guide_assignment_summary.tsv")
    moi = pd.Series(B.sum(axis=1)).value_counts().sort_index()
    summ = {"method": args.method, "layer": args.layer, "guide_umi_limit": args.guide_umi_limit, "merge": args.merge,
            "n_cells": g.n_obs, "n_guides": g.n_vars, "n_assignments": int(B.sum()),
            "cells_with_any_guide": int((B.sum(axis=1) > 0).sum()), "moi_distribution": {int(k): int(v) for k, v in moi.items()},
            "mudata": str(out)}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", f"# Guide assignment ({args.method})\n\nLayer `{args.layer}` added to `{gk}`: "
                 f"{summ['n_assignments']} assignments, {summ['cells_with_any_guide']} of {g.n_obs} cells with a guide.\n\n"
                 + md_table(["guide", "method", "cells assigned"], tab[["guide", "method", "n_assigned"]].values, 30) + "\n")
    return 0


# ---------------------------------------------------------------------------
# Task 3: differential perturbation
# ---------------------------------------------------------------------------

def element_tss(element: str, gvar_rows, rvar, genes_by_name) -> "tuple[str, float]":
    """The TSS used for cis selection: the element's gene (strip _TSS) if known, else the guide coordinate."""
    name = element.replace("_TSS", "")
    if name in genes_by_name:
        return genes_by_name[name]
    r0 = gvar_rows.iloc[0]
    return str(r0["guide_chr"]).replace("chr", ""), float(r0["guide_start"])


def tested_pairs(gv, rvar, distance: int, in_trans: bool, add_genes: "list[str]", min_cells: int, B, compat: bool = False, seed: int = 0):
    """PerturbLoader + runSceptre: element -> genes to test; guides with > min_cells assigned cells."""
    np = _np()
    rng = np.random.default_rng(seed)
    names = rvar["feature_name"].astype(str).to_numpy() if "feature_name" in rvar.columns else rvar.index.astype(str).to_numpy()
    rchr = rvar["transcript_chr"].astype(str).str.replace("chr", "", regex=False).to_numpy()
    rtss = rvar["transcript_start"].astype(float).to_numpy()
    by_name = {n: (c, t) for n, c, t in zip(names, rchr, rtss) if c != "NOT_FOUND"}
    ncells = B.sum(axis=0)
    pairs = []
    for el, rows in gv.groupby("target_elements", sort=False):
        if in_trans:
            genes = set(range(len(names)))
        else:
            nm = el.replace("_TSS", "")
            if compat and nm not in by_name:
                genes = set(rng.choice(len(names), size=min(10, len(names)), replace=False).tolist())
            else:
                c, t = element_tss(el, rows, rvar, by_name)
                genes = set(np.where((rchr == c) & (np.abs(rtss - t) < distance))[0].tolist())
        genes |= {i for i, n in enumerate(names) if n in set(add_genes)}
        for gi in np.where(gv["target_elements"].to_numpy() == el)[0]:
            if ncells[gi] > min_cells:
                for j in sorted(genes):
                    pairs.append((gi, j, el))
    return pairs


def covariate_matrix(obs, rna):
    np = _np()
    cols = [np.ones(rna.n_obs)]
    cands = ["percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count"]
    for c in cands:
        if c in obs.columns:
            v = obs[c].to_numpy(dtype=float)
            if np.isfinite(v).all() and np.std(v) > 0:
                cols.append((v - v.mean()) / v.std())
    if "bath_number" in obs.columns and obs["bath_number"].nunique() > 1:  # dropped when single batch (runSceptre)
        for lv in sorted(obs["bath_number"].unique())[1:]:
            cols.append((obs["bath_number"] == lv).to_numpy(dtype=float))
    if len(cols) == 1:
        tot = np.asarray(rna.X.sum(axis=1)).ravel()
        v = np.log(tot + 1)
        cols.append((v - v.mean()) / max(v.std(), 1e-9))
    return np.column_stack(cols)


def nb_fit(y, Z, max_iter: int = 30):
    """NB2 GLM y ~ Z: Poisson start, method-of-moments theta, IRLS at fixed theta. Returns mu, theta."""
    np = _np()
    _, eta = poisson_glm(y, Z)
    mu = np.exp(np.clip(eta, -30, 30))
    num = float(np.sum(mu ** 2))
    den = float(np.sum((y - mu) ** 2 - mu))
    theta = min(max(num / den, 0.01), 1000.0) if den > 0 else 1000.0
    beta = np.linalg.lstsq(Z, eta, rcond=None)[0]
    for _ in range(max_iter):
        eta = Z @ beta
        mu = np.exp(np.clip(eta, -30, 30))
        w = mu / (1 + mu / theta)
        z = eta + (y - mu) / mu
        ZtW = Z.T * w
        new = np.linalg.solve(ZtW @ Z + 1e-8 * np.eye(Z.shape[1]), ZtW @ z)
        if np.max(np.abs(new - beta)) < 1e-7:
            beta = new
            break
        beta = new
    mu = np.exp(np.clip(Z @ beta, -30, 30))
    return mu, theta


def logistic_propensity(x, Z, max_iter: int = 30):
    np = _np()
    beta = np.zeros(Z.shape[1])
    p0 = min(max(float(x.mean()), 1e-6), 1 - 1e-6)
    beta[0] = math.log(p0 / (1 - p0))
    for _ in range(max_iter):
        eta = Z @ beta
        p = 1 / (1 + np.exp(-eta))
        w = np.maximum(p * (1 - p), 1e-9)
        z = eta + (x - p) / w
        ZtW = Z.T * w
        new = np.linalg.solve(ZtW @ Z + 1e-6 * np.eye(Z.shape[1]), ZtW @ z)
        if np.max(np.abs(new - beta)) < 1e-7:
            beta = new
            break
        beta = new
    return 1 / (1 + np.exp(-(Z @ beta)))


def nb_score_z(Xs, y, mu, theta, Z):
    """Score z-statistics for rows of Xs (k x n) as the treatment indicator, nuisance-adjusted."""
    np = _np()
    w = mu / (1 + mu / theta)
    r = (y - mu) / (1 + mu / theta)
    ZtWZ_inv = np.linalg.inv((Z.T * w) @ Z + 1e-8 * np.eye(Z.shape[1]))
    U = Xs @ r
    A = Xs @ (Z * w[:, None])
    V = (Xs * Xs) @ w - np.einsum("ij,jk,ik->i", A, ZtWZ_inv, A)
    return U / np.sqrt(np.maximum(V, 1e-12))


def resampling_pvalue(z, znull, side: str) -> "tuple[float, str]":
    """Skew-normal fit to the resampled null (sceptre-style), empirical p as fallback/floor."""
    np = _np()
    from scipy import stats  # type: ignore
    B = len(znull)
    emp_r = (1 + np.sum(znull >= z)) / (B + 1)
    emp_l = (1 + np.sum(znull <= z)) / (B + 1)
    try:
        a, loc, scale = stats.skewnorm.fit(znull)
        pr, pl = float(stats.skewnorm.sf(z, a, loc, scale)), float(stats.skewnorm.cdf(z, a, loc, scale))
        how = "skewnorm"
        if not (math.isfinite(pr) and math.isfinite(pl)):
            raise ValueError
    except Exception:
        pr, pl, how = emp_r, emp_l, "empirical"
    if side == "left":
        return pl, how
    if side == "right":
        return pr, how
    return min(1.0, 2 * min(pl, pr)), how


def bh(p):
    np = _np()
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    ranked = p[o] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(adj, 1.0)
    return out


def run_differential(g, r, B, pairs, test: str = "sceptre-nb", side: str = "both", n_resamples: int = 500, seed: int = 0):
    np, pd = _np(), _pd()
    from scipy import stats  # type: ignore
    rng = np.random.default_rng(seed)
    Y = r.X
    Z = covariate_matrix(g.obs if "log_total_guide_count" in g.obs.columns else r.obs, r)
    by_gene: "dict[int, list[tuple[int, str]]]" = {}
    for gi, j, el in pairs:
        by_gene.setdefault(j, []).append((gi, el))
    props: "dict[int, Any]" = {}
    rows = []
    tot = np.asarray(Y.sum(axis=1)).ravel()
    for j, lst in by_gene.items():
        y = dense(Y[:, j]).ravel().astype(float)
        if test == "sceptre-nb":
            mu, theta = nb_fit(y, Z)
        else:
            ln = np.log1p(y / np.maximum(tot, 1) * 1e4)
        for gi, el in lst:
            x = B[:, gi].astype(float)
            n1 = int(x.sum())
            lfc = math.log((y[x > 0].sum() + 0.5) / (mu[x > 0].sum() + 0.5)) if test == "sceptre-nb" else \
                math.log((np.mean((y / np.maximum(tot, 1))[x > 0]) + 1e-9) / (np.mean((y / np.maximum(tot, 1))[x == 0]) + 1e-9))
            if test == "sceptre-nb":
                zval = float(nb_score_z(x[None, :], y, mu, theta, Z)[0])
                if gi not in props:
                    props[gi] = logistic_propensity(x, Z)
                Xs = (rng.uniform(size=(n_resamples, len(x))) < props[gi][None, :]).astype(float)
                znull = nb_score_z(Xs, y, mu, theta, Z)
                p, how = resampling_pvalue(zval, znull, side)
            elif test == "mannwhitney":
                alt = {"both": "two-sided", "left": "less", "right": "greater"}[side]
                res = stats.mannwhitneyu(ln[x > 0], ln[x == 0], alternative=alt)
                p, how = float(res.pvalue), "mannwhitney"
                zval = float(stats.norm.isf(p / (2 if side == "both" else 1)) * (1 if lfc >= 0 else -1))
            else:
                res = stats.ttest_ind(ln[x > 0], ln[x == 0], equal_var=False)
                p = float(res.pvalue)
                if side != "both":
                    p = p / 2 if (res.statistic < 0) == (side == "left") else 1 - p / 2
                zval, how = float(res.statistic), "welch_t"
            rows.append({"gene_id": r.var_names[j], "gRNA_id": g.var_names[gi], "target_elements": el, "pair_type": "all_elements_test",
                         "n_cells_guide": n1, "p_value": p, "z_value": zval, "log_fold_change": lfc, "p_method": how})
    df = pd.DataFrame(rows)
    if len(df):
        df["adj_pvalue"] = bh(df["p_value"].to_numpy())
        df["significant"] = df["adj_pvalue"] < 0.01
    return df


def results_modalities(df, g, r, gv):
    """sceptre_anndata_creation.py: result_guides + result_elements AnnData."""
    import anndata as ad  # type: ignore
    np, pd = _np(), _pd()
    from scipy.stats import combine_pvalues  # type: ignore
    guides = list(dict.fromkeys(df["gRNA_id"]))
    genes = list(dict.fromkeys(df["gene_id"]))
    gi = {k: i for i, k in enumerate(guides)}
    gj = {k: i for i, k in enumerate(genes)}
    mats = {k: np.full((len(guides), len(genes)), np.nan) for k in ["p_value"] + RESULT_LAYERS}
    for row in df.itertuples():
        a, b = gi[row.gRNA_id], gj[row.gene_id]
        for k in mats:
            mats[k][a, b] = float(getattr(row, k))
    obs = gv.loc[guides].copy()
    var = r.var.loc[genes].copy()
    rg = ad.AnnData(mats["p_value"], obs=obs, var=var)
    for k in RESULT_LAYERS:
        rg.layers[k] = mats[k]
    rg.layers["significant"] = np.where(np.isnan(mats["adj_pvalue"]), np.nan, (mats["adj_pvalue"] < 0.01).astype(float))
    els = list(dict.fromkeys(obs["target_elements"]))
    E = np.full((len(els), len(genes)), np.nan)
    for i, e in enumerate(els):
        sub = mats["p_value"][(obs["target_elements"] == e).to_numpy()]
        for j in range(len(genes)):
            v = [x for x in sub[:, j] if not np.isnan(x)]
            if v:
                E[i, j] = combine_pvalues(v).pvalue if len(v) > 1 else v[0]
    first = obs.drop_duplicates("target_elements").set_index("target_elements")
    eobs = pd.DataFrame({"element": els, "element_chr": first.loc[els, "guide_chr"].values,
                         "element_start": first.loc[els, "guide_start"].values, "element_end": first.loc[els, "guide_end"].values},
                        index=pd.Index(els, dtype=str))
    re_ = ad.AnnData(E, obs=eobs, var=var.copy())
    re_.layers["sig_not_adj"] = np.where(np.isnan(E), np.nan, (E < 0.05).astype(float))
    return rg, re_


R_SCEPTRE = """setwd("{d}")
library(sceptre)
library(tibble)
library(dplyr)
exp <- as.matrix(read.table('gene_exp_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
g_I <- as.matrix(read.table('guides_one_gene_guide.txt', sep=',', header=TRUE, row.names=1, check.names=FALSE))
cov <- read.table('covariates.txt', sep=',', header=TRUE, row.names=1)
if ('bath_number' %in% names(cov)) cov$bath_number <- as.factor(cov$bath_number)
pairs_test <- as_tibble(read.table('pairs.txt', sep=',', header=TRUE, check.names=FALSE))
pairs_test$pair_type <- as.factor(pairs_test$pair_type)
result <- run_sceptre_high_moi(gene_matrix = exp, combined_perturbation_matrix = g_I, covariate_matrix = cov,
                               gene_gRNA_group_pairs = pairs_test, side = '{side}')
write.table(result, "results.txt")
"""


def write_r_inputs(rd: Path, g, r, B, pairs, cov, side: str) -> "list[Path]":
    """runSceptre.py layout: one directory per element with the four inputs and r_script.r."""
    pd = _pd()
    out = []
    by_el: "dict[str, list[tuple[int, int]]]" = {}
    for gi, j, el in pairs:
        by_el.setdefault(el, []).append((gi, j))
    for el, lst in by_el.items():
        d = rd / "sceptre_r" / safe_label(el)
        d.mkdir(parents=True, exist_ok=True)
        genes = sorted({j for _, j in lst})
        guides = sorted({gi for gi, _ in lst})
        pd.DataFrame(dense(r.X[:, genes]).T, index=r.var_names[genes], columns=r.obs_names).to_csv(d / "gene_exp_one_gene_guide.txt")
        pd.DataFrame(B[:, guides].T.astype(int), index=g.var_names[guides], columns=r.obs_names).to_csv(d / "guides_one_gene_guide.txt")
        c = cov.copy()
        if "bath_number" in c.columns and c["bath_number"].nunique() == 1:
            del c["bath_number"]
        c.to_csv(d / "covariates.txt")
        pd.DataFrame([(r.var_names[j], g.var_names[gi], "all_elements_test") for gi, j in lst],
                     columns=["gene_id", "gRNA_group", "pair_type"]).to_csv(d / "pairs.txt", index=False)
        (d / "r_script.r").write_text(R_SCEPTRE.format(d=str(d), side=side))
        out.append(d)
    print(f"Wrote: {rd / 'sceptre_r'} ({len(out)} element directories)")
    return out


def cmd_differential(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    mods = read_h5mu(Path(args.mudata))
    gk, rk = guide_key(mods), rna_key(mods)
    g, r = mods[gk], mods[rk]
    gv = annotate_guides(g.var)
    if args.gene_table:
        r.var = add_gene_coords(r.var, read_gene_table(Path(args.gene_table)))
    if "transcript_chr" not in r.var.columns and not args.in_trans:
        raise SystemExit("scRNA.var lacks transcript_chr/transcript_start: pass --gene-table or --in-trans")
    if args.layer in g.layers:
        B = (dense(g.layers[args.layer]) > 0).astype(float)
        bsrc = f"layer {args.layer}"
    else:
        B = (dense(g.X) > args.guide_umi_limit).astype(float)
        bsrc = f"X > {args.guide_umi_limit}"
    pairs = tested_pairs(gv, r.var, args.distance, args.in_trans, args.add_genes or [], args.min_cells_per_guide, B,
                         compat=args.upstream_compat, seed=args.seed)
    rd = run_dir(args.label)
    if not pairs:
        write_json({"n_pairs": 0}, rd / "summary.json")
        write_report(rd / "report.md", "# Differential perturbation\n\nNo testable guide-gene pairs (check --distance, --min-cells-per-guide).\n")
        return 1
    if args.engine == "r":
        cov = g.obs[[c for c in ["bath_number", "percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count"]
                     if c in g.obs.columns]]
        dirs = write_r_inputs(rd, g, r, B, pairs, cov, args.side)
        if shutil.which("Rscript") is None:
            print("Rscript not on PATH; run per element:\n  " + "\n  ".join(f"Rscript {d / 'r_script.r'}" for d in dirs[:5]))
            write_json({"engine": "r", "element_dirs": [str(d) for d in dirs], "ran": False}, rd / "summary.json")
            write_report(rd / "report.md", f"# Differential perturbation (SCEPTRE, R)\n\nInputs for {len(dirs)} elements written; Rscript missing.\n")
            return 0
        for d in dirs:
            subprocess.call(["Rscript", str(d / "r_script.r")])
        res = [pd.read_csv(d / "results.txt", sep=" ") for d in dirs if (d / "results.txt").is_file()]
        df = pd.concat(res, ignore_index=True) if res else pd.DataFrame()
        if "gRNA_group" in df.columns and "gRNA_id" not in df.columns:
            df = df.rename(columns={"gRNA_group": "gRNA_id"})
        if len(df):
            if "log_fold_change" not in df.columns:
                df["log_fold_change"] = np.nan
            df["target_elements"] = [str(x).split("|")[0] for x in df["gRNA_id"]]
            df["adj_pvalue"] = bh(df["p_value"].to_numpy())
            df["significant"] = df["adj_pvalue"] < 0.01
    else:
        df = run_differential(g, r, B, pairs, args.test, args.side, args.n_resamples, args.seed)
    write_tsv(df, rd / "differential_results.tsv.gz")
    rg, re_ = results_modalities(df, g, r, gv)
    out = write_h5mu(rd / "mudata_results.h5mu", {"guides": g, "scRNA": r, "result_guides": rg, "result_elements": re_})
    top = df.sort_values("p_value").head(20)
    summ = {"engine": args.engine, "test": args.test, "side": args.side, "assignment": bsrc, "in_trans": args.in_trans,
            "distance": args.distance, "n_pairs": len(df), "n_guides_tested": int(df["gRNA_id"].nunique()),
            "n_genes_tested": int(df["gene_id"].nunique()), "n_significant": int(df["significant"].sum()),
            "result_guides": list(rg.shape), "result_elements": list(re_.shape), "mudata": str(out),
            "upstream_compat": args.upstream_compat}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", f"# Differential perturbation ({args.test}, side {args.side})\n\n{len(df)} guide-gene tests "
                 f"({summ['n_guides_tested']} guides with > {args.min_cells_per_guide} cells; genes within {args.distance:,} bp of the element"
                 f"{' — in trans: all genes' if args.in_trans else ''}); {summ['n_significant']} significant at BH < 0.01.\n\n"
                 + md_table(["guide", "gene", "cells", "z", "log FC", "p", "adj p"],
                            [(a, b, c, _fmt(d), _fmt(e), _fmt(f), _fmt(h)) for a, b, c, d, e, f, h in
                             top[["gRNA_id", "gene_id", "n_cells_guide", "z_value", "log_fold_change", "p_value", "adj_pvalue"]].values])
                 + f"\n\nResults MuData `{out.name}`: result_guides {rg.n_obs} x {rg.n_vars} (layers {', '.join(RESULT_LAYERS)}), "
                 f"result_elements {re_.n_obs} x {re_.n_vars} (Fisher-combined, sig_not_adj).\n")
    return 0


# ---------------------------------------------------------------------------
# Task 4: tracks
# ---------------------------------------------------------------------------

def results_long(mods: dict):
    """result_guides -> long table of tested pairs with coordinates."""
    np, pd = _np(), _pd()
    rg = mods["result_guides"]
    P = dense(rg.X)
    lay = {k: dense(rg.layers[k]) for k in RESULT_LAYERS if k in rg.layers}
    ii, jj = np.where(~np.isnan(P))
    obs, var = rg.obs, rg.var
    df = pd.DataFrame({"guide": rg.obs_names[ii], "gene": rg.var_names[jj], "p_value": P[ii, jj],
                       "guide_chr": obs["guide_chr"].astype(str).to_numpy()[ii],
                       "guide_start": pd.to_numeric(obs["guide_start"]).to_numpy()[ii],
                       "guide_end": pd.to_numeric(obs["guide_end"]).to_numpy()[ii],
                       "target_elements": obs["target_elements"].astype(str).to_numpy()[ii],
                       "gene_chr": var["transcript_chr"].astype(str).to_numpy()[jj],
                       "gene_start": pd.to_numeric(var["transcript_start"]).to_numpy()[jj],
                       "gene_end": pd.to_numeric(var["transcript_end"]).to_numpy()[jj],
                       "gene_name": (var["feature_name"].astype(str) if "feature_name" in var.columns else var.index.to_series().astype(str)).to_numpy()[jj]})
    for k, M in lay.items():
        df[k] = M[ii, jj]
    return df


def _chr(c: str) -> str:
    c = str(c)
    return c if c.startswith("chr") else f"chr{c}"


def cmd_tracks(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    mods = read_h5mu(Path(args.results))
    df = results_long(mods)
    rd = run_dir(args.label)
    td = rd / "tracks_dir"
    td.mkdir(parents=True, exist_ok=True)
    df["neglog10p"] = -np.log10(np.maximum(df["p_value"].astype(float), 1e-300))
    bed = pd.DataFrame({"chrom": df["guide_chr"].map(_chr), "start": df["guide_start"].astype(int), "end": df["guide_end"].astype(int),
                        "name": df["guide"].str.split("_sgrna").str[0] + "|" + df["gene"].astype(str),
                        "score": np.minimum(1000, np.round(df["neglog10p"] * 100)).astype(int), "strand": "."})
    write_tsv(bed.sort_values(["chrom", "start"]), td / "guide_gene_significance.bed", header=False)
    best = df.sort_values("p_value").drop_duplicates("guide")
    bg = pd.DataFrame({"chrom": best["guide_chr"].map(_chr), "start": best["guide_start"].astype(int), "end": best["guide_end"].astype(int),
                       "value": best["neglog10p"].round(4)})
    write_tsv(bg.sort_values(["chrom", "start"]), td / "guide_min_pvalue.bedgraph", header=False)
    sel = df if args.all_links else df[df["adj_pvalue"] < args.fdr] if "adj_pvalue" in df.columns else df
    tss_e = sel["gene_start"].astype(int) + 1
    links = pd.DataFrame({"chrom1": sel["guide_chr"].map(_chr), "start1": sel["guide_start"].astype(int), "end1": sel["guide_end"].astype(int),
                          "chrom2": sel["gene_chr"].map(_chr), "start2": sel["gene_start"].astype(int), "end2": tss_e,
                          "score": sel["neglog10p"].round(4)})
    write_tsv(links, td / "guide_gene_links.links", header=False)
    ini = (f"[links]\nfile = guide_gene_links.links\ntitle = guide-gene links (-log10 p)\nheight = 3\nlinks_type = arcs\n"
           f"line_width = 1\ncolor = RdYlBu_r\nuse_middle = true\n\n[spacer]\n\n"
           f"[significance]\nfile = guide_gene_significance.bed\ntitle = guides (score = 100 x -log10 p)\nheight = 2\n"
           f"display = collapsed\nlabels = false\n\n[pvalue]\nfile = guide_min_pvalue.bedgraph\ntitle = min -log10 p per guide\n"
           f"height = 2\ncolor = #3b6fb6\n\n[x-axis]\n")
    write_text(td / "tracks.ini", ini)
    cmds, figs = [], []
    plt = _plt()
    for el, sub in df.groupby("target_elements"):
        chrom = _chr(sub["guide_chr"].iloc[0])
        lo = int(min(sub["guide_start"].min(), sub["gene_start"].min())) - args.pad
        hi = int(max(sub["guide_end"].max(), sub["gene_start"].max())) + args.pad
        png = td / f"{safe_label(el)}.png"
        cmds.append(["pyGenomeTracks", "--tracks", str(td / "tracks.ini"), "--region", f"{chrom}:{max(lo, 0)}-{hi}", "-o", str(td / f"{safe_label(el)}_pgt.png")])
        if plt is not None and not args.no_plots:
            fig, ax = plt.subplots(figsize=(7, 2.6))
            g0 = float(sub["guide_start"].mean()) / 1e6
            labelled = set(sub.sort_values("p_value").drop_duplicates("gene").head(3)["gene"])
            done = set()
            for _, rw in sub.iterrows():
                x1, x2 = g0, float(rw["gene_start"]) / 1e6
                h = rw["neglog10p"]
                t = np.linspace(0, np.pi, 60)
                sig = ("adj_pvalue" in rw and rw["adj_pvalue"] < args.fdr)
                ax.plot((x1 + x2) / 2 + (x2 - x1) / 2 * np.cos(t), h * np.sin(t), color=C_SIG if sig else C_NS, lw=1.6 if sig else 0.8)
                ax.plot([x2], [0], marker="|", color=C_SIG if sig else C_NS, ms=6)
                if (sig or rw["gene"] in labelled) and rw["gene"] not in done:
                    done.add(rw["gene"])
                    ax.annotate(str(rw["gene_name"]), (x2, max(h, 0.2)), xytext=(3, 3), textcoords="offset points", fontsize=7,
                                color=C_SIG if sig else INK2)
            ax.axvline(g0, color=C_ELEM, lw=1.5, ls="--")
            ax.set_xlim(lo / 1e6, hi / 1e6)
            ax.set_ylim(bottom=0)
            _style(ax, f"{el}: guide-gene tests (arc height = -log10 p; red = BH < {args.fdr})", f"{chrom} (Mb); dashed = guides", "-log10 p")
            figs.append(str(_save(fig, png)))
    ran = 0
    if shutil.which("pyGenomeTracks") and not args.no_plots:
        for c in cmds:
            ran += subprocess.call(c) == 0
    else:
        print("pyGenomeTracks not on PATH; per-element commands:\n  " + "\n  ".join(" ".join(c) for c in cmds[:5]))
    write_text(td / "pygenometracks_commands.sh", "#!/usr/bin/env bash\n" + "".join(" ".join(c) + "\n" for c in cmds))
    summ = {"n_tests": len(df), "n_links": len(links), "n_elements": int(df["target_elements"].nunique()), "tracks_dir": str(td),
            "figures": figs, "pygenometracks_ran": ran}
    write_json(summ, rd / "summary.json")
    write_report(rd / "report.md", f"# Genome-browser tracks\n\n`tracks_dir/` holds guide_gene_significance.bed ({len(bed)} rows; "
                 f"name = guide|gene, score = 100 x -log10 p), guide_min_pvalue.bedgraph, guide_gene_links.links ({len(links)} "
                 f"pyGenomeTracks links, {'all tests' if args.all_links else f'BH < {args.fdr}'}), tracks.ini and one arc figure per element.\n")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    rc = cmd_assign_guides(argparse.Namespace(mudata=args.mudata, method=args.method, guide_umi_limit=args.guide_umi_limit,
                                              probability_threshold=0.8, merge=False, layer="binarized", seed=args.seed,
                                              label=f"{args.label}_assign"))
    if rc:
        return rc
    a = sorted(OUT_ROOT.glob(f"*_{safe_label(args.label)}_assign"))[-1]
    rc = cmd_differential(argparse.Namespace(mudata=str(a / "mu_with_binary.h5mu"), gene_table=args.gene_table, layer="binarized",
                                             guide_umi_limit=args.guide_umi_limit, distance=args.distance, in_trans=args.in_trans,
                                             add_genes=None, min_cells_per_guide=args.min_cells_per_guide, engine="python",
                                             test=args.test, side=args.side, n_resamples=args.n_resamples, seed=args.seed,
                                             upstream_compat=False, label=f"{args.label}_differential"))
    if rc:
        return rc
    dres = sorted(OUT_ROOT.glob(f"*_{safe_label(args.label)}_differential"))[-1]
    return cmd_tracks(argparse.Namespace(results=str(dres / "mudata_results.h5mu"), all_links=False, fdr=0.01, pad=20000,
                                         no_plots=args.no_plots, label=f"{args.label}_tracks"))


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def synthetic_screen(d: Path, seed: int = 11) -> dict:
    """Two lanes of RNA + guide counts (10x mtx dirs), a GTF, a guide table; planted knock-downs of the element genes."""
    import scipy.io  # type: ignore
    np, pd = _np(), _pd()
    sp = _sp()
    rng = np.random.default_rng(seed)
    n_genes = 320
    chroms = ["1", "2", "3"]
    genes = []
    for i in range(n_genes - 13):
        c = chroms[i % 3]
        genes.append((f"ENSG{i:06d}", f"GENE{i}", c, 1_000_000 + (i // 3) * 150_000, "+" if i % 2 else "-"))
    genes += [(f"ENSGMT{i:02d}", f"MT-CO{i}", "MT", 100 * i + 1, "+") for i in range(13)]
    elements = {"GENE30_TSS": "GENE30", "GENE61_TSS": "GENE61", "GENE92_TSS": "GENE92"}
    effect = {"GENE30": 0.25, "GENE61": 0.3, "GENE92": 1.0}  # GENE92: guides with no effect
    gidx = {g[1]: i for i, g in enumerate(genes)}
    guide_rows = []
    for el, gene in elements.items():
        e = genes[gidx[gene]]
        for k in range(2):
            s = e[3] - 200 + 40 * k
            guide_rows.append({"Target_name": el, "chr": f"chr{e[2]}", "start": s, "end": s + 20, "sgRNA_sequences": "ACGT" * 5})
    for k in range(2):
        guide_rows.append({"Target_name": "random_1", "chr": "chr2", "start": 5_000_000 + 30 * k, "end": 5_000_020 + 30 * k, "sgRNA_sequences": "TTGCA" * 4})
    gt = pd.DataFrame(guide_rows)
    gt_path = d / "guides.csv"
    gt.to_csv(gt_path, index=False)
    named = name_guides(gt)
    gnames = named["pipeline_id"].tolist()
    n_guides = len(gnames)
    base = rng.gamma(1.2, 1.0, size=n_genes)
    base[gidx["GENE30"]] = base[gidx["GENE61"]] = base[gidx["GENE92"]] = 4.0
    base[-13:] *= 2.0
    lanes = []
    truth = {}
    for lane in (1, 2):
        n_cells, n_empty = 450, 120
        cells = [f"L{lane}C{c:04d}" for c in range(n_cells + n_empty)]
        depth = np.concatenate([rng.lognormal(math.log(6000), 0.3, n_cells), rng.lognormal(math.log(60), 0.3, n_empty)])
        # guides: MOI ~ 2 (independent Bernoulli), counts NB high; ambient Poisson low
        assign = rng.uniform(size=(n_cells + n_empty, n_guides)) < 0.12
        assign[n_cells:] = False
        gcount = np.where(assign, rng.negative_binomial(3, 3 / (3 + 40), size=assign.shape) + 8, rng.poisson(0.3, size=assign.shape))
        mult = np.ones((n_cells + n_empty, n_genes))
        for gi, gname in enumerate(gnames):
            el = gname.split("|")[0]
            if el in elements:
                mult[assign[:, gi], gidx[elements[el]]] *= effect[elements[el]]
        w = base[None, :] * mult
        lam = depth[:, None] * w / w.sum(axis=1, keepdims=True)
        X = rng.poisson(lam)
        # 30 high-mito cells
        X[:30, -13:] = X[:30, -13:] * 25 + 20
        ld = d / f"lane{lane}"
        for kind, M, feats in (("rna", X, [(g[0], g[1]) for g in genes]), ("guide", gcount, [(n, n) for n in gnames])):
            dd = ld / kind
            dd.mkdir(parents=True, exist_ok=True)
            with gzip.open(dd / "matrix.mtx.gz", "wb") as fh:
                scipy.io.mmwrite(fh, sp.csr_matrix(M.T).astype(np.int64))
            with gzip.open(dd / "features.tsv.gz", "wt") as fh:
                fh.write("".join(f"{a}\t{b}\t{'Gene Expression' if kind == 'rna' else 'CRISPR Guide Capture'}\n" for a, b in feats))
            with gzip.open(dd / "barcodes.tsv.gz", "wt") as fh:
                fh.write("".join(c + "\n" for c in cells))
        lanes.append(ld)
        truth[lane] = pd.DataFrame(assign, index=cells, columns=gnames)
    gtf = d / "genes.gtf"
    with open(gtf, "w") as fh:
        for gid, gname, c, s, st in genes:
            fh.write(f"{c}\tsynthetic\tgene\t{s}\t{s + 5000}\t.\t{st}\t.\tgene_id \"{gid}.1\"; gene_name \"{gname}\";\n")
    return {"lanes": lanes, "gtf": gtf, "guide_table": gt_path, "gnames": gnames, "truth": truth, "elements": elements, "gidx": gidx}


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    np, pd = _np(), _pd()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def last(label):
        p = sorted(OUT_ROOT.glob(f"*_{label}"))[-1]
        made.append(p)
        return p

    made = []
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_screen(d)

        print("\nassay-spec (Task 1)")
        rc = cmd_assay_spec(argparse.Namespace(assay="10XV2", chemistry=None, whitelist=None, download=False, label="st_jam_assay"))
        a = json.loads((last("st_jam_assay") / "summary.json").read_text())
        check(rc == 0 and a["chemistry"] == "0,0,16:0,16,26:1,0,0" and a["whitelist"] == "737K-august-2016.txt", "assay-spec: 10XV2 -> 0,0,16:0,16,26:1,0,0 + 737K-august-2016")
        rc = cmd_assay_spec(argparse.Namespace(assay="custom", chemistry="0,0,16:0,16,28:1,0,0", whitelist="wl.txt", download=False, label="st_jam_assay2"))
        a2 = json.loads((last("st_jam_assay2") / "summary.json").read_text())
        check(rc == 0 and a2["csv_line"] == "0,0,16:0,16,28:1,0,0,wl.txt", "assay-spec: custom passes CHEMISTRY and WHITELIST through")
        check(cmd_assay_spec(argparse.Namespace(assay="nope", chemistry=None, whitelist=None, download=False, label="x")) == 2,
              "assay-spec: an unknown assay stops the pipeline (exit 2)")

        print("\ncellranger-inputs")
        fq = [str(d / "fq" / "K1000_S1_L001_R1_001.fastq.gz"), str(d / "fq" / "K1000_S1_L001_R2_001.fastq.gz")]
        gq = [str(d / "gq" / "gRNA_S1_L001_R1_001.fastq.gz"), str(d / "gq" / "gRNA_S1_L001_R2_001.fastq.gz")]
        rc = cmd_cellranger_inputs(argparse.Namespace(guide_table=str(W["guide_table"]), rna_fastqs=fq, guide_fastqs=gq, read="R2", pattern="(BC)",
                                                      id="gasperini_01", transcriptome="/ref", localmem=10, threads=5, chemistry=None,
                                                      run=True, label="st_jam_cr"))
        cr = last("st_jam_cr")
        fr = pd.read_csv(cr / "feature_ref.csv")
        lib = pd.read_csv(cr / "library.csv")
        gf = pd.read_csv(cr / "guide_features.txt", sep="\t", header=None)
        check(rc == 0 and list(fr.columns[:6]) == ["id", "name", "read", "pattern", "sequence", "feature_type"] and len(fr) == 8
              and fr["id"].iloc[0] == "GENE30_TSS|1_sgrna_chr1:2499800:2499820", "cellranger-inputs: feature_ref.csv with Target|n_sgrna_chr:start:end ids")
        check(list(lib["library_type"]) == ["Gene Expression", "CRISPR Guide Capture"] and list(lib["sample"]) == ["K1000", "gRNA"],
              "cellranger-inputs: library.csv, one row per FASTQ directory + sample prefix")
        check(gf.shape == (8, 2), "cellranger-inputs: guide_features.txt = sequence <tab> pipeline_id")

        print("\npipeline-config (run_jamboree_v2)")
        os.environ.pop("TOWER_ACCESS_TOKEN", None)
        rc = cmd_pipeline_config(argparse.Namespace(guide_features=str(W["guide_table"]), rna_fastqs=fq, guide_fastqs=gq,
                                                    set=["EXPECTED_CELL_NUMBER=5000", "IN_TRANS=TRUE"], main_nf="pipeline_perturbseq_like/main.nf",
                                                    work_name="gasperini_test_01", resume=True, with_tower=True, nxf_ver="22.10.6",
                                                    run=False, label="st_jam_cfg"))
        cfg = (last("st_jam_cfg") / "perturb.config").read_text()
        check(rc == 0 and "params.EXPECTED_CELL_NUMBER = 5000" in cfg and "params.GUIDE_UMI_LIMIT = 5" in cfg and "params.IN_TRANS = 'TRUE'" in cfg
              and "TOWER" not in cfg, "pipeline-config: notebook defaults + overrides; no Tower token on disk")

        print("\npreprocess")
        lanes = W["lanes"]
        rc = cmd_preprocess(argparse.Namespace(rna=[str(lanes[0] / "rna"), str(lanes[1] / "rna")], guides=[str(lanes[0] / "guide"), str(lanes[1] / "guide")],
                                               lanes=None, expected_cells=440, mito_max=0.2, min_genes=100, mito_genes=None, no_doublets=True,
                                               pct_cells_gene=0.01, gene_table=str(W["gtf"]), guide_umi_limit=5, upstream_compat=False,
                                               no_plots=args.no_plots, label="st_jam_pre"))
        pre = last("st_jam_pre")
        ps = json.loads((pre / "summary.json").read_text())
        mods = read_h5mu(pre / "raw_mudata_guide_and_transcripts.h5mu")
        lane1 = ps["lanes"][0]
        check(rc == 0 and set(mods) == {"guides", "scRNA"} and lane1["after_min_genes"] <= 450 and lane1["after_mito"] <= lane1["after_knee"] - 25,
              f"preprocess: empty droplets (min_genes/knee) and high-mito cells removed ({lane1['barcodes_in']} -> {lane1['after_barcode_intersection']} in lane 1)")
        gm, rm = mods["guides"], mods["scRNA"]
        check(list(gm.var["guide_chr"][:1]) == ["chr1"] and gm.var["target_elements"].iloc[0] == "GENE30_TSS" and
              {"transcript_chr", "transcript_start", "transcript_end"} <= set(rm.var.columns) and (rm.var["transcript_chr"] == "MT").sum() >= 1,
              "preprocess: guide coordinates parsed from pipeline ids; gene TSS from the GTF (version-stripped ids)")
        check({"bath_number", "percent_mito", "log_number_of_detected_genes", "log_total_gene_count", "log_total_guide_count"} <= set(gm.obs.columns)
              and sorted(gm.obs["batch_number"].astype(str).unique()) == ["1", "2"],
              "preprocess: covariates + batch_number per lane in both modalities")
        check(np.allclose(gm.obs["log_total_gene_count"], np.log(rm.obs["n_counts"].astype(float) + 1)),
              "preprocess: log_total_gene_count = log(UMIs + 1) (upstream used n_genes; --upstream-compat)")

        print("\nmudata-qc / distance / subset (Task 2)")
        h5 = str(pre / "raw_mudata_guide_and_transcripts.h5mu")
        rc = cmd_mudata_qc(argparse.Namespace(mudata=h5, n_top=20, max_genes=6500, max_mito=0.1, upstream_compat=False, no_plots=args.no_plots, label="st_jam_qc"))
        q = json.loads((last("st_jam_qc") / "summary.json").read_text())
        check(rc == 0 and q["guide_types"].get("POSITIVE_CONTROL") == 6 and q["guide_types"].get("NEGATIVE_CONTROL") == 2 and q["n_elements"] == 4,
              "mudata-qc: _TSS -> positive, random -> negative controls; element coverage over 4 elements")
        rc = cmd_guide_gene_distance(argparse.Namespace(mudata=h5, how="start", distance=1_000_000, label="st_jam_dist"))
        dd_ = last("st_jam_dist")
        dist = pd.read_csv(dd_ / "guide_gene_distance.tsv.gz", sep="\t")
        g1 = dist[(dist["guide"] == "GENE30_TSS|1_sgrna_chr1:2499800:2499820") & (dist["gene"] == "ENSG000030")]
        check(rc == 0 and len(g1) == 1 and g1["distance"].iloc[0] == abs(2499800 - 2505000) and not (dist["gene"].str.startswith("ENSGMT")).any(),
              "guide-gene-distance: |guide start - gene TSS| on the same chromosome only")
        rc = cmd_subset(argparse.Namespace(mudata=h5, element=None, category=["NEGATIVE_CONTROL"], guide=None, keep_modalities=None, out_name=None,
                                           distance=1_000_000, label="st_jam_sub"))
        sub = read_h5mu(last("st_jam_sub") / "subset.h5mu")
        check(rc == 0 and sub["guides"].n_vars == 2 and sub["scRNA"].n_vars == rm.n_vars, "subset: negative-control guides only")
        rc = cmd_subset(argparse.Namespace(mudata=h5, element=None, category=None, guide=["GENE30_TSS|1_sgrna_chr1:2499800:2499820"],
                                           keep_modalities=None, out_name=None, distance=1_000_000, label="st_jam_sub2"))
        s2 = json.loads((last("st_jam_sub2") / "summary.json").read_text())
        check(rc == 0 and s2["cells_with_guide"] > 50 and "ENSG000030" in s2["cis_genes"], "subset --guide: guide presence per cell + cis gene counts")

        print("\nassign-guides")
        rc = cmd_assign_guides(argparse.Namespace(mudata=h5, method="umi", guide_umi_limit=5, probability_threshold=0.8, merge=False,
                                                  layer="binarized", seed=0, label="st_jam_asg"))
        asg = last("st_jam_asg")
        mb = read_h5mu(asg / "mu_with_binary.h5mu")
        B = dense(mb["guides"].layers["binarized"])
        tr = pd.concat([W["truth"][1], W["truth"][2]]).loc[mb["guides"].obs_names, mb["guides"].var_names].to_numpy()
        acc = float((B == tr).mean())
        check(rc == 0 and acc > 0.99, f"assign-guides umi (> 5 UMIs): {acc:.4f} agreement with the planted assignment")
        rc = cmd_assign_guides(argparse.Namespace(mudata=h5, method="poisson-mixture", guide_umi_limit=5, probability_threshold=0.8, merge=False,
                                                  layer="binarized", seed=0, label="st_jam_asg2"))
        mb2 = read_h5mu(last("st_jam_asg2") / "mu_with_binary.h5mu")
        acc2 = float((dense(mb2["guides"].layers["binarized"]) == tr).mean())
        check(rc == 0 and acc2 > 0.98, f"assign-guides poisson-mixture (covariate offsets): {acc2:.4f} agreement")
        rc = cmd_assign_guides(argparse.Namespace(mudata=h5, method="umi", guide_umi_limit=5, probability_threshold=0.8, merge=True,
                                                  layer="binarized", seed=0, label="st_jam_asg3"))
        mb3 = read_h5mu(last("st_jam_asg3") / "mu_with_binary.h5mu")
        check(rc == 0 and mb3["guides"].n_vars == 4, "assign-guides --merge: guides summed per element before binarising")

        print("\ndifferential (Task 3)")
        rc = cmd_differential(argparse.Namespace(mudata=str(asg / "mu_with_binary.h5mu"), gene_table=None, layer="binarized", guide_umi_limit=5,
                                                 distance=1_000_000, in_trans=False, add_genes=None, min_cells_per_guide=30, engine="python",
                                                 test="sceptre-nb", side="both", n_resamples=300, seed=1, upstream_compat=False, label="st_jam_de"))
        de = last("st_jam_de")
        res = pd.read_csv(de / "differential_results.tsv.gz", sep="\t")
        own = res[res.apply(lambda r: r["target_elements"].replace("_TSS", "") in ("GENE30", "GENE61") and
                            r["gene_id"] == f"ENSG{int(r['target_elements'][4:].split('_')[0]):06d}", axis=1)]
        check(rc == 0 and len(own) == 4 and own["significant"].all() and (own["log_fold_change"] < -0.8).all(),
              f"differential: planted knock-downs significant for all 4 guides (log FC {', '.join(_fmt(x, 2) for x in own['log_fold_change'])})")
        null = res[~res.index.isin(own.index)]
        check(null["significant"].mean() < 0.02 and 0.3 < null["p_value"].median() < 0.7,
              f"differential: other pairs null ({int(null['significant'].sum())} of {len(null)} significant; median p {null['p_value'].median():.2f})")
        g92 = res[res["target_elements"] == "GENE92_TSS"]
        check(len(g92) > 0 and not g92["significant"].any(), "differential: the no-effect element (GENE92) finds nothing")
        rnd = res[res["target_elements"] == "random_1"]
        check(len(rnd) > 0 and not rnd["significant"].any(),
              f"differential: negative controls tested around their guide coordinate ({len(rnd)} pairs, none significant)")
        mr = read_h5mu(de / "mudata_results.h5mu")
        rg = mr["result_guides"]
        check(set(mr) == {"guides", "scRNA", "result_guides", "result_elements"} and set(RESULT_LAYERS) <= set(rg.layers)
              and np.isnan(dense(rg.X)).any(), "differential: mudata_results.h5mu with result_guides (NaN = untested, 4 layers) + result_elements")
        re_ = mr["result_elements"]
        check("sig_not_adj" in re_.layers and list(re_.obs.columns) == ["element", "element_chr", "element_start", "element_end"],
              "differential: result_elements Fisher-combined per element with sig_not_adj")
        rc = cmd_differential(argparse.Namespace(mudata=str(asg / "mu_with_binary.h5mu"), gene_table=None, layer="binarized", guide_umi_limit=5,
                                                 distance=1_000_000, in_trans=False, add_genes=None, min_cells_per_guide=30, engine="python",
                                                 test="mannwhitney", side="left", n_resamples=100, seed=1, upstream_compat=False, label="st_jam_de2"))
        r2 = pd.read_csv(last("st_jam_de2") / "differential_results.tsv.gz", sep="\t")
        check(rc == 0 and r2.set_index(["gRNA_id", "gene_id"]).loc[list(zip(own["gRNA_id"], own["gene_id"])), "significant"].all(),
              "differential --test mannwhitney: the alternative module recovers the same knock-downs")
        rc = cmd_differential(argparse.Namespace(mudata=str(asg / "mu_with_binary.h5mu"), gene_table=None, layer="binarized", guide_umi_limit=5,
                                                 distance=1_000_000, in_trans=False, add_genes=None, min_cells_per_guide=30, engine="r",
                                                 test="sceptre-nb", side="both", n_resamples=100, seed=1, upstream_compat=False, label="st_jam_der"))
        der = last("st_jam_der")
        el_dirs = sorted((der / "sceptre_r").iterdir())
        check(rc == 0 and len(el_dirs) == 4 and all((x / "pairs.txt").is_file() and (x / "r_script.r").is_file() for x in el_dirs),
              "differential --engine r: runSceptre per-element inputs + run_sceptre_high_moi script")

        print("\ntracks (Task 4)")
        rc = cmd_tracks(argparse.Namespace(results=str(de / "mudata_results.h5mu"), all_links=False, fdr=0.01, pad=20000, no_plots=args.no_plots, label="st_jam_tr"))
        trd = last("st_jam_tr") / "tracks_dir"
        bed = pd.read_csv(trd / "guide_gene_significance.bed", sep="\t", header=None)
        lk = pd.read_csv(trd / "guide_gene_links.links", sep="\t", header=None)
        check(rc == 0 and bed.shape[1] == 6 and len(bed) == len(res) and bed[0].str.startswith("chr").all(), "tracks: BED of every test (name = guide|gene)")
        check(len(lk) == int(res["significant"].sum()) and lk.shape[1] == 7 and (trd / "tracks.ini").is_file(),
              f"tracks: {len(lk)} pyGenomeTracks links (BH < 0.01) + tracks.ini")
        if not args.no_plots:
            check(len(list(trd.glob("*.png"))) == 4, "tracks: one arc figure per element")

        print("\nrun")
        rc = cmd_run(argparse.Namespace(mudata=h5, method="umi", guide_umi_limit=5, gene_table=None, distance=1_000_000, in_trans=False,
                                        min_cells_per_guide=30, test="sceptre-nb", side="both", n_resamples=200, seed=2, no_plots=True, label="st_jam_run"))
        for suf in ("assign", "differential", "tracks"):
            last(f"st_jam_run_{suf}")
        check(rc == 0, "run: assign-guides -> differential -> tracks")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    try:
        OUT_ROOT.rmdir()
    except OSError:
        pass
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.0f} s)")
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent crispr-fg-jamboree", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("assay-spec", help="Task 1: assay -> kallisto chemistry string + whitelist")
    p.add_argument("--assay", required=True, help=f"one of {', '.join(ASSAYS)}, or custom")
    p.add_argument("--chemistry", help="kallisto -x string for --assay custom")
    p.add_argument("--whitelist", help="whitelist path for --assay custom")
    p.add_argument("--download", action="store_true", help="download the whitelist into the run directory")
    p.add_argument("--label", default="assay_spec")
    p.set_defaults(func=cmd_assay_spec)

    p = sub.add_parser("cellranger-inputs", help="guide table -> feature_ref.csv, library.csv, guide_features.txt + cellranger count")
    p.add_argument("--guide-table", required=True, help="xlsx/csv/tsv with Target_name, chr, start, end, sgRNA_sequences")
    p.add_argument("--rna-fastqs", nargs="+")
    p.add_argument("--guide-fastqs", nargs="+")
    p.add_argument("--read", default="R2", help="feature_ref read (R1/R2)")
    p.add_argument("--pattern", default="(BC)", help="feature_ref pattern, e.g. (BC) or 5PNNNNNNNNNN(BC)")
    p.add_argument("--id", default="gasperini_01")
    p.add_argument("--transcriptome", default="refdata-gex-GRCh38-2020-A")
    p.add_argument("--localmem", type=int, default=10)
    p.add_argument("--threads", type=int, default=10)
    p.add_argument("--chemistry")
    p.add_argument("--run", action="store_true", help="run cellranger count if it is on PATH")
    p.add_argument("--label", default="cellranger_inputs")
    p.set_defaults(func=cmd_cellranger_inputs)

    p = sub.add_parser("pipeline-config", help="perturb.config + nextflow launch of pipeline_perturbseq_like")
    p.add_argument("--guide-features")
    p.add_argument("--rna-fastqs", nargs="+")
    p.add_argument("--guide-fastqs", nargs="+")
    p.add_argument("--set", nargs="+", help="KEY=VALUE overrides of the config params")
    p.add_argument("--main-nf", default="pipeline_perturbseq_like/main.nf")
    p.add_argument("--work-name", default="gasperini_test_01")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--with-tower", action="store_true", help="add -with-tower when TOWER_ACCESS_TOKEN is set")
    p.add_argument("--nxf-ver", default="22.10.6")
    p.add_argument("--run", action="store_true", help="launch nextflow if it is on PATH")
    p.add_argument("--label", default="pipeline_config")
    p.set_defaults(func=cmd_pipeline_config)

    p = sub.add_parser("preprocess", help="RNA + guide counts -> raw MuData (pipeline QC)")
    p.add_argument("--rna", nargs="+", required=True, help="per-lane RNA counts (.h5ad or 10x matrix dir)")
    p.add_argument("--guides", nargs="+", required=True, help="per-lane guide counts (.h5ad or 10x matrix dir)")
    p.add_argument("--lanes", nargs="+")
    p.add_argument("--expected-cells", type=int, default=10000, help="EXPECTED_CELL_NUMBER (knee index)")
    p.add_argument("--mito-max", type=float, default=0.2, help="MITO_EXPECTED_PERCENTAGE")
    p.add_argument("--min-genes", type=int, default=100)
    p.add_argument("--mito-genes", help="file of mitochondrial gene ids/names (default: symbols starting MT-)")
    p.add_argument("--no-doublets", action="store_true")
    p.add_argument("--pct-cells-gene", type=float, default=0.01, help="PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS")
    p.add_argument("--gene-table", help="GTF(.gz) or TSV with gene_name/gene_id, chr, start, end, strand")
    p.add_argument("--guide-umi-limit", type=int, default=5)
    p.add_argument("--upstream-compat", action="store_true", help="upstream covariate definitions")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="preprocess")
    p.set_defaults(func=cmd_preprocess)

    p = sub.add_parser("mudata-qc", help="Task 2: QC tables and figures of a guides/scRNA MuData")
    p.add_argument("--mudata", required=True)
    p.add_argument("--n-top", type=int, default=20)
    p.add_argument("--max-genes", type=int, default=6500)
    p.add_argument("--max-mito", type=float, default=0.1)
    p.add_argument("--upstream-compat", action="store_true", help="apply only the mito filter (the notebook overwrote the gene filter)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="mudata_qc")
    p.set_defaults(func=cmd_mudata_qc)

    p = sub.add_parser("guide-gene-distance", help="Task 2: guide x gene distances and cis pairs")
    p.add_argument("--mudata", required=True)
    p.add_argument("--how", choices=["start", "midpoint"], default="start")
    p.add_argument("--distance", type=int, default=1_000_000)
    p.add_argument("--label", default="guide_gene_distance")
    p.set_defaults(func=cmd_guide_gene_distance)

    p = sub.add_parser("subset", help="Task 2/3: subset guides by element/category/guide or keep some modalities")
    p.add_argument("--mudata", required=True)
    p.add_argument("--element", nargs="+")
    p.add_argument("--category", nargs="+", help="POSITIVE_CONTROL NEGATIVE_CONTROL PUTATIVE_ENHANCER")
    p.add_argument("--guide", nargs="+")
    p.add_argument("--keep-modalities", nargs="+", help="e.g. guides scRNA (drop result_* modalities)")
    p.add_argument("--out-name")
    p.add_argument("--distance", type=int, default=1_000_000)
    p.add_argument("--label", default="subset")
    p.set_defaults(func=cmd_subset)

    p = sub.add_parser("assign-guides", help="guide calling -> 'binarized' layer")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["umi", "poisson-mixture"], default="umi")
    p.add_argument("--guide-umi-limit", type=int, default=5, help="GUIDE_UMI_LIMIT: assigned when UMI > limit")
    p.add_argument("--probability-threshold", type=float, default=0.8)
    p.add_argument("--merge", action="store_true", help="sum guides of the same element first")
    p.add_argument("--layer", default="binarized")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="assign_guides")
    p.set_defaults(func=cmd_assign_guides)

    def add_de(p):
        p.add_argument("--distance", type=int, default=1_000_000, help="DISTANCE_NEIGHBORS")
        p.add_argument("--in-trans", action="store_true", help="test every gene (IN_TRANS=TRUE)")
        p.add_argument("--min-cells-per-guide", type=int, default=30)
        p.add_argument("--test", choices=["sceptre-nb", "mannwhitney", "welch-t"], default="sceptre-nb")
        p.add_argument("--side", choices=["both", "left", "right"], default="both")
        p.add_argument("--n-resamples", type=int, default=500)
        p.add_argument("--seed", type=int, default=0)
        p.add_argument("--gene-table")
        p.add_argument("--guide-umi-limit", type=int, default=5)

    p = sub.add_parser("differential", help="Task 3: differential perturbation -> mudata_results.h5mu")
    p.add_argument("--mudata", required=True)
    p.add_argument("--layer", default="binarized", help="guide assignment layer (else X > --guide-umi-limit)")
    p.add_argument("--add-genes", nargs="+", help="genes tested for every element (ADDGENENAMES)")
    p.add_argument("--engine", choices=["python", "r"], default="python")
    p.add_argument("--upstream-compat", action="store_true", help="10 random genes for elements not in the gene table")
    add_de(p)
    p.add_argument("--label", default="differential")
    p.set_defaults(func=cmd_differential)

    p = sub.add_parser("tracks", help="Task 4: BED / bedGraph / links / tracks.ini / figures")
    p.add_argument("--results", required=True, help="mudata_results.h5mu")
    p.add_argument("--all-links", action="store_true", help="links for every test, not only BH < --fdr")
    p.add_argument("--fdr", type=float, default=0.01)
    p.add_argument("--pad", type=int, default=20000)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="tracks")
    p.set_defaults(func=cmd_tracks)

    p = sub.add_parser("run", help="assign-guides -> differential -> tracks on a raw MuData")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["umi", "poisson-mixture"], default="umi")
    add_de(p)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="jamboree")
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
