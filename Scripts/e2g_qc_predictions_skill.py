#!/usr/bin/env python3
"""E2G QC, prediction packaging and IGVF Portal submission workflow (port of kaybrand/QC-and-Predictions).

Port of https://github.com/kaybrand/QC-and-Predictions (MIT, Kayla Brand 2026), default branch
`igvf-portal-submission` at c647bdbb9b6843f07034683fb2c10034f9c4dfc0 (2026-08-29T03:49:06Z), plus the
parts of the `CATlas-predictions` (e5e2799c) and `synapse-submission` (92968bb2) branches that differ.
Relationship: port -- every Snakemake rule, Python/R script, config and the igvf_cell_annotation_report
folder was read and the definitions re-derived in Python; no code was copied.  The upstream pipeline takes
QC-filtered clusters from pseudobulk to shareable scE2G data products: quality gating, ATAC/RNA filtering,
scE2G configs, portal-format reformatting, candidate/feature tables, QC aggregation, the Cell Annotation
cache built from the IGVF Portal, and a plan-then-execute IGVF Portal submission of eleven metadata tables
(Synapse manifests on the side).  scE2G itself is `igvfagent sce2g-pipeline` (reused here for its QC
statistics, QC figures, bgzip/tabix and model thresholds); this module is everything around it.

Definitions reproduced from upstream
  quality gate      per (dataset, cluster): cell_count / fragments_total / umi_count summed over the rows of
                    plots/<ds>/<cl>/filtered_cell_subsample_metrics.tsv (sibling of the QC guide), or, for
                    `prefiltered: true`, the QC-guide barcodes joined to <datatables>/<ds>_data/<cl>_per_cell_qc.tsv
                    (num_frags / rna_read_count).  Reasons, first hit wins: user_specified, missing_qc_guide |
                    missing_metrics | missing_per_cell_qc_table, below_min_cell_count (default 100),
                    below_min_fragments_total (2e6), below_min_umi_count (1e6; skipped for ATAC-only clusters),
                    pass.  predictions_on_everything_but_do_not_upload (legacy key process_excluded_no_upload)
                    re-includes quality/user exclusions but never missing-input ones; upload-eligible = passed;
                    manifest-eligible = upload-eligible minus igvf_manifest_excluded.  --legacy-implicit-source
                    reproduces the synapse-submission/CATlas rule (metrics file iff default guide name and file
                    exists, else datatable join; "no stats" = not excluded).  CATlas `atac_frag_file` clusters
                    read only the metrics file written by prefiltered-metrics.
  merged metrics    component metrics grouped by subsample: n_cells / total_fragments / total_RNA_reads summed,
                    mean_frag_per_cell and mean_RNA_per_cell recomputed as ratios, mean_frip / mean_tss
                    n_cells-weighted; values written with 15 significant digits.
  prefiltered       CATlas: one metrics row per cluster from an already-filtered fragments file: n_cells =
                    distinct barcodes, total_fragments = lines, RNA 0, frip / tss NA.
  ATAC filter       full-barcode match (never the 16-bp prefix) of fragment column 4 against the guide's first
                    column, directories annotation-<ct>-IGVF* (comma list = merged cluster; fallback
                    annotation-<ct>), exactly 5 fields per record, every guide barcode must be seen (else exit 1;
                    --allow-missing-barcodes = legacy warning), chromosomes absent from chrom sizes dropped,
                    ordered by chrom-sizes order then start, end; bgzip + tabix -p bed.  --clean keeps bc[:16].
  RNA filter        rna_counts_mtx.h5ad per directory, obs filtered by full name, inner concat; Ensembl ->
                    symbol by summing columns (var gene_symbol, or exact versioned GTF gene_id; any unmatched id
                    is a hard failure; --standard-chromosomes-only keeps chr1-22,X,Y,M); four sanity checks
                    (no duplicate cells/genes, cells == guide barcodes); .mtx directory (genes x cells, gz),
                    .h5ad or dense .csv.gz; package: flat tar.gz of decompressed matrix.mtx/barcodes.tsv/features.tsv.
  reformat          update_scE2G_pred_formats.R header (# Source/Version/GenomeReference IGVFDS0280IQAI/URL/
                    Assays 10x Multiome/SampleAgnostic False/SampleTermName/SampleTermID/CellAnnotation
                    [/ScoreThreshold 'Score >= t'][/ScoreType positive_score unless a _list file]
                    [/Metadata https://data.igvf.org/tabular-files/<alias>]) and column sets: full (feature.-prefixed
                    model features, Score = E2G.Score.qnorm, raw ElementName), thresholded (narrow, ElementName
                    chr:start-end), gene list (TSS +/- 250), element list, element BED from
                    Neighborhoods/EnhancerList.bed (header block + '#ElementChr..' line, bgzip + tabix), bedpe
                    (sort -k1,1 -k2,2n, bgzip, tabix).  --format synapse-legacy = synapse-submission branch
                    (SampleSummaryShort, PRELIMINARY banner, Score.ignoreTPM); --catlas drops normalizedATAC_enh
                    from scATAC full files (CATlas branch).
  candidates        ElementChr..ElementClass, GeneTSS, GeneSymbol, GeneEnsemblID, SampleSummaryShort,
                    E2G_Distance, isSelfPromoter; no header.  features: '# Source: scE2G <model>' header, fixed
                    leading columns then every remaining column (dplyr everything(), so `chr` survives).
  QC aggregation    newest scE2G_predictions_threshold*_stats.tsv per (cluster, model_name) over every cluster
                    directory on disk; cell_count 0 patched from the metrics file; plots as plot_all_qc_stats.R
                    (dashed lines at 2e6 fragments, 100 cells, 1e6 UMIs).
  Cell Annotation   one PseudobulkSet multireport; primary = every input_file_sets @id under /analysis-sets/,
                    principal = all /pseudobulk-sets/; primaries with exactly one sample and an alias are cached;
                    a cluster resolves when every QC-guide subsample has alias suffix <ds>-<cl>-<subsample>
                    (merged clusters: any constituent name); (cl_id, term_id, term_name) must agree; annotation /
                    qualifier from the most-contributing subsample (merged: bare term_name, no qualifier); locked
                    when all primaries are released or a principal carries the annotation.  Snapshot TSV carries
                    portal_fetched_at + a 16-hex sha256 of the cluster set, both enforced on read (24 h).
  stale reformats   delete reformatted files whose '# CellAnnotation/SampleTermName/SampleTermID' header
                    disagrees with the cache (and their .tbi).
  IGVF manifests    eleven tables / fifteen variants (QC_documents, principal_pseudobulk_set,
                    filtered_barcode_list, filtered_atac_fragment_file, filtered_rna_count_matrix,
                    atac_index_file, prediction_set, prediction_tabular_files {full, thresholded, bedpe,
                    elements_bed, genes}, signal_files, elements_bed_index_file, bedpe_index_file) with the exact
                    aliases ("jesse-engreitz:" + '_'.join(parts)), constant fields, controlled vocabulary,
                    reference files, analysis-step versions and derived_from / file_set / input_file_sets
                    links; payload = {aliases, award, lab} + constants + scope fields + row; required columns
                    must be non-empty (False is a value); round = 1 + max(dependency rounds); only Multiome
                    (enabled_families) cluster_model rows; outcomes planned-post / planned-patch / unchanged /
                    deferred / skipped-family-gated / skipped-missing-file / invalid / enabled-check-failed;
                    round<N>_<table>[_<variant>]_{post,patch}.tsv in iu_register.py format (no csv quoting,
                    arrays comma-joined, objects JSON, booleans lower-case, record_id first on patches); a
                    PATCH retires a stale sibling POST, never the reverse; <object_type>.tsv accumulator; sha256
                    payload hash in a SQLite ledger.
  submission        DRY-RUN by default: payloads, iu_register TSVs, upload_plan.json and lineage checks are
                    written and validated; only `manifest --execute` POSTs / PATCHes (REST JSON, arrays split,
                    attachment {"path"} -> base64 data URI, md5sum computed), in round order, to the IGVF
                    sandbox (the Portal redirects sandbox -> staging) unless --production, with
                    IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY from the environment, re-reading every alias
                    afterwards; rows whose dependencies are not uploaded are deferred to the next pass.
  Synapse           manifest diff (owned vs should-exist, longest-cluster-name match, preserved foreign rows,
                    confirm-delete / confirm-overwrite gates) and orphan listings; dry-run unless --execute.

Subcommands
  resolve-exclusions  quality gate -> cluster_stats/<ds>_cluster_stats.tsv + <ds>_pipeline_plan.tsv
  merge-metrics       merge_cluster_metrics.py            prefiltered-metrics  compute_prefiltered_cell_metrics.py
  build-qc-datatables build_qc_datatables.py              filter-atac / filter-rna / package-rna
  sce2g-config        cell_clusters + cluster_metadata tables, scE2G config overlay, model thresholds
  reformat            one file (full | thresholded | gene-list | element-list | element-bed | bedpe)
  candidates / features   candidate E2G pairs / scE2G feature table
  aggregate-qc        all_qc_stats.tsv + QC figures        stale-reformats   header vs Cell Annotation cache
  cell-metadata       multireport -> cache -> per-cluster annotations, snapshot, status (also --catlas-seed)
  cell-annotation-report  igvf_cell_annotation_report (report-all or --qc-guide-dir mode)
  dataset-accessions  dataset -> principal analysis set accession(s)
  washu-report / verify-fragments   CATlas WashU pseudobulk report and fragments md5 verification
  portal-files / compare-archive    primary-pseudobulk file discovery / portal-vs-archive comparison
  manifest            IGVF payloads, rounds, validation, lineage check, upload plan; --execute submits
  patch-submitter-comment  "Version 1" backfill plan         report   per-cluster coverage report.tsv
  synapse-manifest / synapse-orphans   Synapse manifest diff / orphan listing
  distance-depth      CATlas distance-to-TSS vs depth figures
  run                 driver: preflight, warm, local packaging, manifest preview, audit (exit 0/1/2)
  selftest            synthetic datasets with planted QC failures, annotations and files; every subcommand

Usage:
    igvfagent e2g-qc-predictions resolve-exclusions --config igvf0_pipeline_config.yaml
    igvfagent e2g-qc-predictions filter-atac --qc-guide guide.tsv.gz --pseudobulks igvf0/pseudobulks --cell-type k562 --chrom-sizes hg38.sizes --out atac.tsv.gz
    igvfagent e2g-qc-predictions reformat --kind full --input scE2G_predictions.tsv.gz --out ds_cl_scE2G_multiome_powerlaw_v3.e2g.tsv.gz --model multiome_powerlaw_v3 --cell-type K562 --term-id EFO:0002067 --summary K562
    igvfagent e2g-qc-predictions cell-metadata --config cfg.yaml --multireport-json pseudobulk_sets.json
    igvfagent e2g-qc-predictions manifest --config cfg.yaml --cluster-keys igvf0 --label igvf0
    igvfagent e2g-qc-predictions manifest --config cfg.yaml --cluster-keys igvf0 --execute      # sandbox
    igvfagent e2g-qc-predictions run --config cfg.yaml --label igvf0
    igvfagent e2g-qc-predictions selftest --no-plots
"""
from __future__ import annotations

import argparse
import base64
import csv
import fcntl
import glob
import gzip
import hashlib
import io
import json
import logging
import math
import mimetypes
import os
import re
import shutil
import sqlite3
import sys
import tarfile
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "E2GQCPredictions"
DEFAULT_STATE_DB = OUT_ROOT / "igvf_metadata_state.db"

UPSTREAM_REPO = "kaybrand/QC-and-Predictions"
UPSTREAM_COMMIT = "c647bdbb9b6843f07034683fb2c10034f9c4dfc0"
UPSTREAM_BRANCH = "igvf-portal-submission"
UPSTREAM_DATE = "2026-08-29T03:49:06Z"
BRANCH_COMMITS = {"CATlas-predictions": "e5e2799c8b99df0a089dc304e3ef2bd182cd709f",
                  "synapse-submission": "92968bb2487d6729388e3e3216ec1f9479b8691d"}

DEFAULT_QC_GUIDE_NAME = "filtered_barcodes_with_subsamples.tsv.gz"
METRICS_NAME = "filtered_cell_subsample_metrics.tsv"
UNPROCESSABLE_REASONS = frozenset({"missing_qc_guide", "missing_metrics", "missing_per_cell_qc_table"})
PREDICT_EVERYTHING_KEY = "predictions_on_everything_but_do_not_upload"
PREDICT_EVERYTHING_LEGACY_KEY = "process_excluded_no_upload"
CLUSTER_STATS_HEADER = ["dataset", "cluster", "cell_count", "fragments_total", "umi_count", "reason"]
PLAN_HEADER = ["dataset", "cluster", "included", "upload_eligible", "manifest_eligible", "exclusion_reason",
               "cell_count", "fragments_total", "umi_count"]
ATAC_ONLY = ["scATAC_powerlaw_v3"]
MULTIOME_MODEL = "multiome_powerlaw_v3"
SCATAC_MODEL = "scATAC_powerlaw_v3"
CELL_CLUSTERS_HEADER = ["cluster", "rna_matrix_file", "atac_frag_file", "HiC_file", "HiC_type", "HiC_resolution",
                        "alt_TSS", "alt_genes", "model_dir"]
CLUSTER_METADATA_HEADER = ["cluster", "ontology_id", "cell_type", "summary"]
MAX_MEM_MB = 250 * 1000

METRIC_ADDITIVE = ("n_cells", "total_fragments", "total_RNA_reads")
METRIC_RATIOS = {"mean_frag_per_cell": ("total_fragments", "n_cells"), "mean_RNA_per_cell": ("total_RNA_reads", "n_cells")}
METRIC_WEIGHTED = ("mean_frip", "mean_tss")
METRIC_COLUMNS = ["subsample", *METRIC_ADDITIVE, "mean_frag_per_cell", "mean_RNA_per_cell", *METRIC_WEIGHTED]
PER_CELL_QC_COLUMNS = ["analysis_accession", "barcode", "subsample", "rna_read_count", "gene_count", "pct_mito",
                       "pct_ribo", "num_frags", "pct_duplicated_reads", "nucleosomal_signal", "tss_enrichment", "frip"]
DIRNAME_RE = re.compile(r"^annotation-(?P<annotation>.+)-(?P<subsample>IGVFSM\w+)$")
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

GENOME_REFERENCE = "IGVFDS0280IQAI"
SCE2G_URL = "https://github.com/EngreitzLab/scE2G/tree/main"
PORTAL_TABULAR_URL = "https://data.igvf.org/tabular-files/"
REFORMAT_SUFFIXES = (".e2g.tsv.gz", "_element_list.bed.gz", "_gene_list.tsv.gz")
STALE_HEADER_FIELDS = {"# CellAnnotation:": "cell_annotation", "# SampleTermName:": "term_name",
                       "# SampleTermID:": "term_id"}

LOG = logging.getLogger("e2g_qc_predictions")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"e2g_qc_predictions_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    print(f"Log: {path}")
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label))[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir_for(args, default_label: str) -> Path:
    if getattr(args, "out_dir", None):
        d = Path(args.out_dir)
        d.mkdir(parents=True, exist_ok=True)
        return d
    return run_dir(getattr(args, "label", None) or default_label)


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy: pip install 'igvfagent[analysis]'") from e


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


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_jsonable(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return o.item()
        except Exception:
            return str(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_json(path: Path, obj, quiet: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(obj), indent=2, sort_keys=False) + "\n")
    if not quiet:
        print(f"JSON: {path}")
    return path


def write_rows_tsv(path: Path, header: "Sequence[str]", rows: "Iterable[dict]", quiet: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(header), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in header})
    os.replace(tmp, path)
    if not quiet:
        print(f"TSV: {path}")
    return path


def read_rows_tsv(path) -> "List[dict]":
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md_table(headers: "Sequence[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def write_report(out: Path, title: str, lines: "List[str]", summary: dict) -> None:
    rp = out / "report.md"
    rp.write_text(f"# {title}\n\n" + "\n".join(lines) + "\n")
    print(f"Report: {rp}")
    write_json(out / "summary.json", summary)


def md5_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def content_md5(path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _open_text(path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode.replace("t", "") or "r")


def bgzip_tabix(plain: Path, preset: str = "bed") -> Path:
    """bgzip + tabix of `plain` -> plain.gz (+ .tbi); sce2g-pipeline's pysam implementation."""
    from sce2g_pipeline_skill import bgzip_tabix as _bt  # type: ignore
    return _bt(Path(plain), preset=preset)


def fmt_value(v) -> str:
    """data.table::fwrite rendering: NA for missing, TRUE/FALSE, 15 significant digits."""
    if v is None:
        return "NA"
    if isinstance(v, bool) or type(v).__name__ == "bool_":
        return "TRUE" if bool(v) else "FALSE"
    if isinstance(v, float) or type(v).__name__.startswith("float"):
        f = float(v)
        if not math.isfinite(f):
            return "NA" if f != f else ("Inf" if f > 0 else "-Inf")
        if f == int(f) and abs(f) < 1e15:
            return str(int(f))
        return f"{f:.15g}"
    s = str(v)
    return "NA" if s == "nan" else s


# ---------------------------------------------------------------------------
# Config (YAML subset, no PyYAML needed)
# ---------------------------------------------------------------------------

def _yaml_scalar(s: str):
    s = s.strip()
    if s == "" or s in ("~", "null", "Null", "NULL"):
        return None
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        return s[1:-1]
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        return [_yaml_scalar(x) for x in inner.split(",")] if inner else []
    if s.startswith("{") and s.endswith("}"):
        inner = s[1:-1].strip()
        out = {}
        for part in inner.split(","):
            if ":" in part:
                k, v = part.split(":", 1)
                out[k.strip()] = _yaml_scalar(v)
        return out
    low = s.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def _strip_comment(line: str) -> str:
    out, q = [], None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch
        elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
            break
        out.append(ch)
    return "".join(out).rstrip()


def parse_yaml_subset(text: str):
    """Block mappings, block lists of scalars, flow lists / maps and scalars -- the shape every upstream
    *_pipeline_config.yaml uses.  PyYAML is used instead when it is installed."""
    lines = []
    for raw in text.splitlines():
        s = _strip_comment(raw)
        if s.strip():
            lines.append((len(s) - len(s.lstrip(" ")), s.strip()))

    def parse_block(i: int, indent: int):
        if i < len(lines) and lines[i][1].startswith("- "):
            out = []
            while i < len(lines) and lines[i][0] == indent and lines[i][1].startswith("- "):
                out.append(_yaml_scalar(lines[i][1][2:]))
                i += 1
            return out, i
        out = {}
        while i < len(lines) and lines[i][0] == indent:
            key, _, rest = lines[i][1].partition(":")
            key = key.strip().strip("'\"")
            rest = rest.strip()
            i += 1
            if rest:
                out[key] = _yaml_scalar(rest)
            elif i < len(lines) and lines[i][0] > indent:
                out[key], i = parse_block(i, lines[i][0])
            elif i < len(lines) and lines[i][0] == indent and lines[i][1].startswith("- "):
                out[key], i = parse_block(i, indent)
            else:
                out[key] = None
        return out, i

    if not lines:
        return {}
    val, _ = parse_block(0, lines[0][0])
    return val


def load_config(path) -> dict:
    text = Path(path).read_text()
    if str(path).endswith(".json"):
        return json.loads(text)
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text) or {}
    except ImportError:
        return parse_yaml_subset(text)


def resolve_repo_relative(path: str, repo_root) -> str:
    return path if os.path.isabs(path) else os.path.abspath(os.path.join(str(repo_root), path))


def config_output_dir(config: dict, config_path: "Optional[str]" = None) -> str:
    base = Path(config_path).resolve().parent if config_path else Path.cwd()
    return resolve_repo_relative(str(config.get("output_dir", "./results")), base)


def config_bool(value, default=True) -> bool:
    """common.smk _config_bool: survives --config key=false; raises on an unrecognised value."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("false", "0", "no", "off"):
        return False
    if text in ("true", "1", "yes", "on"):
        return True
    raise ValueError(f"expected a boolean-ish value (true/false/1/0/yes/no/on/off), got {value!r}")


def iter_clusters(config: dict):
    for dataset, clusters in (config.get("clusters") or {}).items():
        for cluster, cfg in (clusters or {}).items():
            yield dataset, cluster, cfg or {}


def is_atac_only(cluster_cfg: dict) -> bool:
    return list(cluster_cfg.get("models") or []) == ATAC_ONLY


def resolve_primary_model(models) -> str:
    """Exactly one candidates / feature table per cluster: multiome_powerlaw_v3 when it ran, else scATAC."""
    return MULTIOME_MODEL if MULTIOME_MODEL in (models or []) else SCATAC_MODEL


def annotation_lookup_key(dataset: str, cluster: str, cluster_cfg: dict) -> "Tuple[str, str]":
    """cell_annotations.py: ATAC-only variant clusters resolve under their cell_annotation_key."""
    return dataset, (cluster_cfg or {}).get("cell_annotation_key", cluster)


def determine_mem_mb(input_size_mb: float, gz: bool, attempt: int = 1, min_gb: int = 8) -> int:
    """common.smk determine_mem_mb (ABC's calculator)."""
    size = input_size_mb * (8 if gz else 1)
    return int(min((2 ** (attempt - 1)) * max(2 * size, min_gb * 1000), MAX_MEM_MB))


def model_threshold(scE2G_dir: "Optional[str]", model: str, overrides: "Optional[dict]" = None) -> str:
    """get_model_threshold / _score_threshold: models/<model>/score_threshold_<t> (leading zero restored).
    Falls back to the explicit override, then to the sce2g-pipeline embedded threshold (deviation: upstream
    raises when the marker directory is absent)."""
    if overrides and model in overrides:
        t = str(overrides[model])
    else:
        t = None
        if scE2G_dir:
            matches = glob.glob(os.path.join(str(scE2G_dir), "models", model, "score_threshold_*"))
            if len(matches) > 1:
                raise ValueError(f"expected exactly one score_threshold_* for {model}, found {len(matches)}")
            if matches:
                m = re.search(r"score_threshold_([0-9.]+)$", os.path.basename(matches[0]))
                if not m:
                    raise ValueError(f"couldn't parse a threshold decimal out of {matches[0]}")
                t = m.group(1)
        if t is None:
            try:
                from sce2g_pipeline_skill import MODEL_SPECS  # type: ignore
                t = MODEL_SPECS[model]["score_threshold"]
            except Exception:
                raise ValueError(f"no score_threshold_* marker for model {model!r} and no override")
    if t.startswith("."):
        t = "0" + t
    return t


# ---------------------------------------------------------------------------
# Quality gate (resolve_exclusions.py)
# ---------------------------------------------------------------------------

def predictions_on_everything(exclusion_cfg: dict) -> bool:
    new = exclusion_cfg.get(PREDICT_EVERYTHING_KEY)
    legacy = exclusion_cfg.get(PREDICT_EVERYTHING_LEGACY_KEY)
    if new is not None and legacy is not None:
        print(f"WARNING: both exclusion.{PREDICT_EVERYTHING_KEY} and the deprecated "
              f"{PREDICT_EVERYTHING_LEGACY_KEY} are set -- using {new}", file=sys.stderr)
        return bool(new)
    if legacy is not None:
        print(f"NOTE: exclusion.{PREDICT_EVERYTHING_LEGACY_KEY} is deprecated -- rename it to {PREDICT_EVERYTHING_KEY}",
              file=sys.stderr)
        return bool(legacy)
    return bool(new)


def read_guide_rows(qc_guide_path) -> "List[dict]":
    return read_rows_tsv(qc_guide_path)


def guide_barcodes(qc_guide_path) -> "set":
    """First column of the guide, header skipped (filter_*.py load_passing_barcodes)."""
    out = set()
    with _open_text(qc_guide_path) as fh:
        fh.readline()
        for line in fh:
            line = line.strip()
            if line:
                out.add(line.split("\t")[0])
    return out


def stats_from_subsample_metrics(path) -> "Tuple[int, int, int]":
    c = f = u = 0
    for row in read_rows_tsv(path):
        c += int(float(row["n_cells"]))
        f += int(float(row["total_fragments"]))
        u += int(float(row["total_RNA_reads"]))
    return c, f, u


def stats_from_per_cell_qc_join(qc_guide_path, per_cell_qc_path) -> "Tuple[int, int, int]":
    barcodes = {r["barcode"] for r in read_guide_rows(qc_guide_path)}
    c = f = u = 0
    for row in read_rows_tsv(per_cell_qc_path):
        if row["barcode"] in barcodes:
            c += 1
            f += int(float(row["num_frags"]))
            u += int(float(row["rna_read_count"]))
    return c, f, u


def compute_cluster_stats(datatables_dir: str, data_dir: str, dataset: str, cluster: str, cluster_cfg: dict,
                          legacy_implicit: bool = False) -> "Tuple[Optional[dict], str]":
    has_rna = not is_atac_only(cluster_cfg)
    if "atac_frag_file" in cluster_cfg:  # CATlas: already-filtered fragments, metrics from prefiltered-metrics
        mp = os.path.join(data_dir, "plots", dataset, cluster, METRICS_NAME)
        if not os.path.exists(mp):
            return None, "missing_metrics"
        c, f, u = stats_from_subsample_metrics(mp)
        return {"cell_count": c, "fragments_total": f, "umi_count": u if has_rna else None}, "ok"
    qc_guide = cluster_cfg.get("qc_guide") or ""
    metrics_path = os.path.join(os.path.dirname(qc_guide), METRICS_NAME)
    if legacy_implicit:
        if os.path.basename(qc_guide) == DEFAULT_QC_GUIDE_NAME and os.path.exists(metrics_path):
            c, f, u = stats_from_subsample_metrics(metrics_path)
        else:
            pcq = os.path.join(datatables_dir, f"{dataset}_data",
                               f"{cluster_cfg.get('pseudobulk_annotation', cluster)}_per_cell_qc.tsv")
            if not os.path.exists(qc_guide) or not os.path.exists(pcq):
                return None, "missing_metrics" if not os.path.exists(pcq) else "missing_qc_guide"
            c, f, u = stats_from_per_cell_qc_join(qc_guide, pcq)
    elif not cluster_cfg.get("prefiltered", False):
        if not os.path.exists(metrics_path):
            return None, "missing_metrics"
        c, f, u = stats_from_subsample_metrics(metrics_path)
    else:
        if not os.path.exists(qc_guide):
            return None, "missing_qc_guide"
        pcq = os.path.join(datatables_dir, f"{dataset}_data", f"{cluster}_per_cell_qc.tsv")
        if not os.path.exists(pcq):
            return None, "missing_per_cell_qc_table"
        c, f, u = stats_from_per_cell_qc_join(qc_guide, pcq)
    return {"cell_count": c, "fragments_total": f, "umi_count": u if has_rna else None}, "ok"


def resolve_exclusions(config: dict, legacy_implicit: bool = False):
    """(included, upload_eligible, excluded, stats_by_cluster), keyed by (dataset, cluster)."""
    data_dir = str(config.get("data_dir") or "")
    datatables_dir = config.get("qc_datatables_dir") or os.path.join(data_dir, "datatables")
    exc = config.get("exclusion") or {}
    user = set()
    for name in exc.get("user_specified") or []:
        d, _, c = str(name).partition("/")
        user.add((d, c))
    everything = predictions_on_everything(exc)
    thr = exc.get("auto_thresholds") or {}
    min_cells = thr.get("min_cell_count", 0) or 0
    min_frags = thr.get("min_fragments_total", 0) or 0
    min_umi = thr.get("min_umi_count", 0) or 0
    excluded, stats_by, all_clusters = set(), {}, set()
    for dataset, cluster, cfg in iter_clusters(config):
        key = (dataset, cluster)
        all_clusters.add(key)
        stats, why = compute_cluster_stats(datatables_dir, data_dir, dataset, cluster, cfg, legacy_implicit)
        if key in user:
            excluded.add(key)
            reason = "user_specified"
        elif stats is None:
            if legacy_implicit:
                reason = why  # legacy: no stats yet -> not excluded (only user_specified applies)
            else:
                excluded.add(key)
                reason = why
        elif stats["cell_count"] < min_cells:
            excluded.add(key)
            reason = "below_min_cell_count"
        elif stats["fragments_total"] < min_frags:
            excluded.add(key)
            reason = "below_min_fragments_total"
        elif stats["umi_count"] is not None and stats["umi_count"] < min_umi:
            excluded.add(key)
            reason = "below_min_umi_count"
        else:
            reason = "pass"
        stats_by[key] = {"cell_count": stats["cell_count"] if stats else None,
                         "fragments_total": stats["fragments_total"] if stats else None,
                         "umi_count": stats["umi_count"] if stats else None, "reason": reason}
    if everything:
        if legacy_implicit:
            included = set(all_clusters)
        else:
            processable = {k for k in excluded if stats_by[k]["reason"] not in UNPROCESSABLE_REASONS}
            included = (all_clusters - excluded) | processable
    else:
        included = all_clusters - excluded
    upload = included - excluded
    return included, upload, excluded, stats_by


def manifest_eligible_clusters(config: dict, upload_eligible) -> "set":
    clusters = config.get("clusters") or {}
    return {(d, c) for d, c in upload_eligible if not (clusters[d][c] or {}).get("igvf_manifest_excluded", False)}


def write_cluster_stats_table(dataset: str, stats_by: dict, out_dir, quiet: bool = False) -> Path:
    rows = [{"dataset": ds, "cluster": cl, **row} for (ds, cl), row in sorted(stats_by.items()) if ds == dataset]
    return write_rows_tsv(Path(out_dir) / f"{dataset}_cluster_stats.tsv", CLUSTER_STATS_HEADER, rows, quiet=quiet)


def write_plan_tsv(dataset: str, all_keys, included, upload, manifest, stats_by, out_dir, quiet=False) -> Path:
    rows = []
    for ds, cl in sorted(k for k in all_keys if k[0] == dataset):
        s = stats_by[(ds, cl)]
        rows.append({"dataset": ds, "cluster": cl, "included": "y" if (ds, cl) in included else "n",
                     "upload_eligible": "y" if (ds, cl) in upload else "n",
                     "manifest_eligible": "y" if (ds, cl) in manifest else "n", "exclusion_reason": s["reason"],
                     "cell_count": s["cell_count"], "fragments_total": s["fragments_total"], "umi_count": s["umi_count"]})
    return write_rows_tsv(Path(out_dir) / f"{dataset}_pipeline_plan.tsv", PLAN_HEADER, rows, quiet=quiet)


# ---------------------------------------------------------------------------
# Metrics tables: merge, prefiltered, per-cell QC datatables
# ---------------------------------------------------------------------------

def _fmt15(v) -> str:
    return str(v) if isinstance(v, int) else f"{v:.15g}"


def merge_metrics(paths: "Sequence[str]"):
    acc, seen = {}, defaultdict(list)
    for path in paths:
        rows = read_rows_tsv(path)
        if rows:
            missing = [c for c in METRIC_COLUMNS if c not in rows[0]]
            if missing:
                raise ValueError(f"{path}: missing column(s) {missing}")
        for row in rows:
            sub = row["subsample"]
            seen[sub].append(os.path.basename(os.path.dirname(os.path.abspath(path))))
            cur = acc.setdefault(sub, {**{c: 0 for c in METRIC_ADDITIVE}, **{f"_w_{c}": 0.0 for c in METRIC_WEIGHTED}})
            n = int(float(row["n_cells"]))
            for c in METRIC_ADDITIVE:
                cur[c] += int(float(row[c]))
            for c in METRIC_WEIGHTED:
                cur[f"_w_{c}"] += float(row[c]) * n
    out = []
    for sub in sorted(acc):
        cur = acc[sub]
        r = {"subsample": sub, **{c: cur[c] for c in METRIC_ADDITIVE}}
        for col, (num, den) in METRIC_RATIOS.items():
            r[col] = (cur[num] / cur[den]) if cur[den] else 0
        for col in METRIC_WEIGHTED:
            r[col] = (cur[f"_w_{col}"] / cur["n_cells"]) if cur["n_cells"] else 0
        out.append(r)
    return out, {s: v for s, v in seen.items() if len(v) > 1}


def write_metrics(rows, out_path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    with open(tmp, "w") as fh:
        fh.write("\t".join(METRIC_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r["subsample"] if c == "subsample" else _fmt15(r[c]) for c in METRIC_COLUMNS) + "\n")
    os.replace(tmp, out_path)
    print(f"TSV: {out_path}")
    return out_path


def prefiltered_metrics(frag_file, subsample_name: str) -> dict:
    barcodes, total = set(), 0
    with _open_text(frag_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            total += 1
            barcodes.add(line.split("\t", 4)[3].rstrip("\n"))
    n = len(barcodes)
    return {"subsample": subsample_name, "n_cells": n, "total_fragments": total, "total_RNA_reads": 0,
            "mean_frag_per_cell": total / n if n else "NA", "mean_RNA_per_cell": 0, "mean_frip": "NA", "mean_tss": "NA"}


def per_cell_qc_path(directory) -> "Optional[str]":
    for name in ("per_cell_qc.tsv.gz", "per_cell_qc.tsv"):
        p = os.path.join(directory, name)
        if os.path.exists(p):
            return p
    return None


def annotation_dirs(pseudobulks_root, dataset, annotation) -> "List[str]":
    pat = os.path.join(str(pseudobulks_root), dataset, "pseudobulks", f"annotation-{annotation}-IGVFSM*")
    return sorted(d for d in glob.glob(pat) if os.path.isdir(d))


def discover_annotations(pseudobulks_root, dataset) -> "List[str]":
    base = os.path.join(str(pseudobulks_root), dataset, "pseudobulks")
    if not os.path.isdir(base):
        return []
    found = set()
    for e in os.listdir(base):
        m = DIRNAME_RE.match(e)
        if m and os.path.isdir(os.path.join(base, e)):
            found.add(m.group("annotation"))
    return sorted(found)


def build_one_datatable(pseudobulks_root, dataset, cluster, annotations, out_path) -> "Tuple[int, int, List[str]]":
    sources, skipped = [], []
    for ann in annotations:
        dirs = annotation_dirs(pseudobulks_root, dataset, ann)
        if not dirs:
            skipped.append(f"{ann}:no_directories")
            continue
        for d in dirs:
            p = per_cell_qc_path(d)
            if p is None:
                skipped.append(f"{os.path.basename(d)}:no_per_cell_qc")
            else:
                sources.append(p)
    if not sources:
        return 0, 0, skipped
    header, total = None, 0
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    try:
        with open(tmp, "w", newline="") as out:
            for p in sources:
                with _open_text(p) as fh:
                    first = fh.readline()
                    if header is None:
                        header = first
                        cols = header.rstrip("\n").split("\t")
                        if cols != PER_CELL_QC_COLUMNS:
                            raise ValueError(f"{p}: unexpected columns {cols}")
                        out.write(header)
                    elif first != header:
                        raise ValueError(f"{p}: header differs from the first source in this cluster")
                    for line in fh:
                        out.write(line)
                        total += 1
    except ValueError:
        if tmp.exists():
            tmp.unlink()
        raise
    os.replace(tmp, out_path)
    return total, len(sources), skipped


# ---------------------------------------------------------------------------
# ATAC / RNA filtering (filter_atac_fragments.py, filter_rna_counts.py, package_rna_count_matrix)
# ---------------------------------------------------------------------------

def find_pseudobulk_dirs(pseudobulks_dir, cell_type: str, merged: bool = True) -> "List[str]":
    """annotation-<ct>-IGVF* per comma-separated cell type (merged cluster), fallback annotation-<ct>."""
    out = []
    names = [c.strip() for c in cell_type.split(",")] if merged else [cell_type]
    for ct in names:
        dirs = sorted(d for d in glob.glob(os.path.join(str(pseudobulks_dir), f"annotation-{ct}-IGVF*")) if os.path.isdir(d))
        if dirs:
            out.extend(dirs)
            continue
        fb = os.path.join(str(pseudobulks_dir), f"annotation-{ct}")
        if os.path.isdir(fb):
            out.append(fb)
        else:
            print(f"[warning] No directories found matching annotation-{ct}-IGVF*", file=sys.stderr)
    return out


def read_chrom_order(path) -> "List[str]":
    out = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(line.split()[0])
    return out


def filter_atac(qc_guide, pseudobulks, cell_type, chrom_sizes, out, strict_fields: bool = True,
                allow_missing: bool = False, clean: bool = False) -> dict:
    passing = guide_barcodes(qc_guide)
    dirs = find_pseudobulk_dirs(pseudobulks, cell_type)
    if not dirs:
        raise SystemExit("[error] No fragment directories found.")
    rows, found, per_dir, total = [], set(), [], 0
    for d in dirs:
        fp = os.path.join(d, "fragments.tsv.gz")
        if not os.path.isfile(fp):
            print(f"[warning] fragments.tsv.gz not found in {d}, skipping.", file=sys.stderr)
            continue
        n_pass = n_tot = 0
        with _open_text(fp) as fh:
            for ln, line in enumerate(fh, 1):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if strict_fields and len(f) != 5:
                    raise ValueError(f"Malformed fragment record in {fp}, line {ln}: expected 5 tab-separated fields "
                                     f"(chrom, start, end, barcode, duplicate count), got {len(f)}")
                if len(f) < 4:
                    continue
                n_tot += 1
                bc = f[3]
                if bc in passing:
                    if clean:
                        f[3] = bc[:16]
                    rows.append(f)
                    found.add(bc)
                    n_pass += 1
        per_dir.append({"dir": os.path.basename(d), "retained": n_pass, "total": n_tot})
        total += n_tot
    missing = sorted(passing - found)
    if missing and not allow_missing:
        raise SystemExit(f"[error] {len(missing)} barcode(s) from the QC guide were not found in any fragment file: "
                         + ", ".join(missing[:10]))
    order = read_chrom_order(chrom_sizes)
    rank = {c: i for i, c in enumerate(order)}
    dropped = Counter(r[0] for r in rows if r[0] not in rank)
    kept = [r for r in rows if r[0] in rank]
    kept.sort(key=lambda r: (rank[r[0]], int(r[1]), int(r[2])))
    out = str(out) if str(out).endswith(".gz") else str(out) + ".gz"
    plain = Path(out[:-3])
    plain.parent.mkdir(parents=True, exist_ok=True)
    with open(plain, "w") as fh:
        for r in kept:
            fh.write("\t".join(r) + "\n")
    gz = bgzip_tabix(plain, preset="bed")
    return {"out": str(gz), "n_passing_barcodes": len(passing), "rows_total": total, "rows_retained": len(kept),
            "missing_barcodes": missing, "dropped_by_chrom": dict(dropped), "per_directory": per_dir}


def parse_gtf_map(gtf_path) -> "Tuple[Dict[str, Tuple[str, str]], list]":
    gid, gname = re.compile(r'gene_id "([^"]+)"'), re.compile(r'gene_name "([^"]+)"')
    m, conflicts = {}, []
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            a, b = gid.search(parts[8]), gname.search(parts[8])
            if a and b:
                vid, name = a.group(1), b.group(1)
                if vid in m:
                    if m[vid][1] != name:
                        conflicts.append((vid, m[vid][0], m[vid][1], parts[0], name))
                else:
                    m[vid] = (parts[0], name)
    return m, conflicts


def _collapse(adata, symbols):
    import anndata as ad  # type: ignore
    import scipy.sparse as sps  # type: ignore
    np, pd = _np(), _pd()
    uniq, inv = np.unique(np.asarray(symbols), return_inverse=True)
    S = sps.csr_matrix((np.ones(adata.n_vars), (np.arange(adata.n_vars), inv)), shape=(adata.n_vars, len(uniq)))

    def coll(mat):
        mat = mat if sps.issparse(mat) else sps.csr_matrix(mat)
        return sps.csr_matrix(mat @ S)
    var = pd.DataFrame(index=uniq)
    var.index.name = "gene_symbol"
    return ad.AnnData(X=coll(adata.X), obs=adata.obs.copy(), var=var,
                      layers={k: coll(adata.layers[k]) for k in adata.layers}), int((np.bincount(inv) > 1).sum())


def filter_rna(qc_guide, pseudobulks, cell_type, out, ensembl_ids: bool = False, gtf: "Optional[str]" = None,
               standard_only: bool = False, log_path: "Optional[str]" = None) -> dict:
    try:
        import anndata as ad  # type: ignore
        import scipy.io  # type: ignore
        import scipy.sparse as sps  # type: ignore
    except ImportError:
        raise SystemExit("filter-rna needs anndata + scipy (pip install anndata)")
    np, pd = _np(), _pd()
    if standard_only and not gtf:
        raise SystemExit("[error] --standard-chromosomes-only requires --gtf.")
    if gtf and ensembl_ids:
        raise SystemExit("[error] --gtf cannot be used with --ensemblIDs-as-genes.")
    name = str(out).lower()
    fmt = ("csv" if name.endswith((".csv", ".csv.gz")) else "h5ad" if name.endswith((".h5ad", ".h5"))
           else "mtx" if name.endswith((".mtx", ".mtx.gz")) else None)
    if fmt is None:
        raise SystemExit(f"[error] Cannot infer output format from {out!r}")
    passing = guide_barcodes(qc_guide)
    dirs = find_pseudobulk_dirs(pseudobulks, cell_type)
    adatas, per_dir = [], []
    for d in dirs:
        hp = os.path.join(d, "rna_counts_mtx.h5ad")
        if not os.path.isfile(hp):
            continue
        a = ad.read_h5ad(hp)
        n0 = a.n_obs
        a = a[a.obs_names.isin(list(passing))].copy()
        per_dir.append({"dir": os.path.basename(d), "cells_before": n0, "cells_after": a.n_obs})
        if a.n_obs:
            adatas.append(a)
    if not adatas:
        raise SystemExit("[error] No passing cells found across any directory.")
    combined = ad.concat(adatas, axis=0, join="inner", merge="same")
    n_ens = combined.n_vars
    log_lines, overloaded, unmatched, nonstd = [], 0, [], 0
    if not ensembl_ids:
        if gtf:
            gmap, conflicts = parse_gtf_map(gtf)
            ids = list(combined.var_names)
            unmatched = sorted(e for e in ids if e not in gmap)
            log_lines += [f"GTF file: {gtf}", f"Total Ensembl IDs in anndata: {len(ids)}",
                          f"Matched to GTF (exact versioned): {len(ids) - len(unmatched)}", f"Unmatched: {len(unmatched)}"]
            if unmatched:
                if log_path:
                    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
                    Path(log_path).write_text("\n".join(log_lines + ["", *unmatched]) + "\n")
                raise SystemExit(f"[error] {len(unmatched)} Ensembl IDs in anndata not found in GTF: {unmatched[:10]}")
            if standard_only:
                keep = np.array([gmap[e][0] in STANDARD_CHROMS for e in ids])
                nonstd = int((~keep).sum())
                combined = combined[:, keep].copy()
                ids = list(combined.var_names)
            combined, overloaded = _collapse(combined, [gmap[e][1] for e in ids])
            log_lines += [f"Ensembl IDs on nonstandard chromosomes removed: {nonstd}",
                          f"Unique gene symbols after collapse: {combined.n_vars}", f"Symbols with multiple Ensembl IDs (summed): {overloaded}"]
            if log_path:
                Path(log_path).parent.mkdir(parents=True, exist_ok=True)
                Path(log_path).write_text("\n".join(log_lines) + "\n")
        else:
            combined, overloaded = _collapse(combined, combined.var["gene_symbol"].values)
    errors = []
    if pd.Series(list(combined.obs_names)).duplicated().any():
        errors.append("duplicate cell barcodes")
    if pd.Series(list(combined.var_names)).duplicated().any():
        errors.append("duplicate gene identifiers")
    if combined.n_obs != len(passing):
        errors.append(f"cell count mismatch: matrix has {combined.n_obs} cells, QC guide has {len(passing)}")
    if errors:
        raise SystemExit("[error] Sanity check(s) failed: " + "; ".join(errors))
    if fmt == "csv":
        p = str(out) if str(out).endswith(".gz") else str(out) + ".gz"
        X = combined.X.toarray() if sps.issparse(combined.X) else np.asarray(combined.X)
        df = pd.DataFrame(X, index=combined.obs_names, columns=combined.var_names)
        df.index.name = "barcode"
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, compression="gzip")
        outputs = [p]
    elif fmt == "h5ad":
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        combined.write_h5ad(str(out))
        outputs = [str(out)]
    else:
        base = str(out)
        for ext in (".mtx.gz", ".mtx"):
            if base.lower().endswith(ext):
                base = base[: -len(ext)]
                break
        os.makedirs(base, exist_ok=True)
        X = combined.X if sps.issparse(combined.X) else sps.csr_matrix(combined.X)
        with gzip.open(os.path.join(base, "matrix.mtx.gz"), "wb") as fh:
            scipy.io.mmwrite(fh, X.T)
        with gzip.open(os.path.join(base, "barcodes.tsv.gz"), "wt") as fh:
            fh.writelines(bc + "\n" for bc in combined.obs_names)
        with gzip.open(os.path.join(base, "features.tsv.gz"), "wt") as fh:
            fh.writelines(g + "\n" for g in combined.var_names)
        outputs = [os.path.join(base, n) for n in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]
    return {"outputs": outputs, "n_cells": int(combined.n_obs), "n_genes": int(combined.n_vars), "n_ensembl_in": int(n_ens),
            "symbols_summed": overloaded, "nonstandard_removed": nonstd, "per_directory": per_dir, "format": fmt}


def package_rna_matrix(matrix_dir, tarball) -> Path:
    """Flat tar.gz of the DEcompressed matrix.mtx / barcodes.tsv / features.tsv (Filtered Matrix File spec)."""
    tarball = Path(tarball)
    tarball.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        for n in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
            with gzip.open(os.path.join(str(matrix_dir), n + ".gz"), "rb") as fi, open(os.path.join(td, n), "wb") as fo:
                shutil.copyfileobj(fi, fo)
        with tarfile.open(tarball, "w:gz") as tf:
            for n in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
                tf.add(os.path.join(td, n), arcname=n)
    return tarball


# ---------------------------------------------------------------------------
# scE2G configuration (write_scE2G_config.py, common.smk build_scE2G_config)
# ---------------------------------------------------------------------------

def load_lab_annotations(path) -> dict:
    out = {}
    if path and os.path.exists(path):
        for row in read_rows_tsv(path):
            out[(row["dataset"].strip().upper(), row["lab_celltype"].strip().lower())] = row
    return out


def resolve_cell_type_and_ontology(dataset: str, annotation: str, lab: dict) -> "Tuple[str, str]":
    row = lab.get((dataset.upper(), annotation.lower()))
    if row is None:
        return "TODO: cell_type (no lab_annotations_with_cl.tsv row found)", "TODO: ontology_id"
    cl_term, qual, oid = row.get("CL term", "").strip(), row.get("qualifier", "").strip(), row.get("CL_ID", "").strip()
    return (" ".join(p for p in (cl_term, qual) if p) or "TODO: cell_type (blank CL term/qualifier)",
            oid or "TODO: ontology_id (blank CL_ID)")


def _locked_merge_write(path: Path, header, key_col: str, merge_fn) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(str(path) + ".lock", os.O_CREAT | os.O_RDWR)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        existing = {r[key_col]: r for r in read_rows_tsv(path)} if path.exists() else {}
        merged = merge_fn(existing)
        write_rows_tsv(path, header, [merged[k] for k in sorted(merged)], quiet=True)
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)
    return path


def write_cell_clusters_table(dataset, clusters_cfg, included, out_dir, data_dir_ds, scE2G_dir) -> Path:
    def merge(existing):
        for cl in included:
            cfg = clusters_cfg[cl]
            frag = cfg["atac_frag_file"] if "atac_frag_file" in cfg else os.path.join(data_dir_ds, cl, f"atac_fragments_{dataset}_{cl}.tsv.gz")
            existing[cl] = {"cluster": cl, "rna_matrix_file": "" if is_atac_only(cfg) else
                            os.path.join(data_dir_ds, cl, f"rna_count_matrix_{dataset}_{cl}"),
                            "atac_frag_file": frag, "HiC_file": "", "HiC_type": "", "HiC_resolution": "", "alt_TSS": "",
                            "alt_genes": "", "model_dir": ",".join(os.path.join(str(scE2G_dir), "models", m) for m in cfg["models"])}
        return existing
    return _locked_merge_write(Path(out_dir) / f"{dataset}_cell_clusters.tsv", CELL_CLUSTERS_HEADER, "cluster", merge)


def write_cluster_metadata_table(dataset, clusters_cfg, included, out_dir, lab_annotations_path=None) -> Path:
    lab = load_lab_annotations(lab_annotations_path)

    def merge(existing):
        for cl in included:
            cfg, cur = clusters_cfg[cl], existing.get(cl, {})
            if "cell_type" not in cur or str(cur.get("cell_type", "TODO")).startswith("TODO"):
                ct, oid = resolve_cell_type_and_ontology(dataset, cfg.get("pseudobulk_annotation", cl), lab)
            else:
                ct, oid = cur["cell_type"], cur["ontology_id"]
            existing[cl] = {"cluster": cl, "ontology_id": oid, "cell_type": ct, "summary": cl}
        return existing
    return _locked_merge_write(Path(out_dir) / f"{dataset}_cluster_metadata.tsv", CLUSTER_METADATA_HEADER, "cluster", merge)


def make_paths_absolute(obj, base):
    if isinstance(obj, dict):
        return {k: make_paths_absolute(v, base) for k, v in obj.items()}
    if isinstance(obj, list):
        return [make_paths_absolute(v, base) for v in obj]
    if isinstance(obj, str) and os.path.exists(os.path.join(str(base), obj)):
        return os.path.join(str(base), obj)
    return obj


def build_sce2g_config(config: dict, cell_clusters_table: str, results_dir: str) -> dict:
    scE2G_dir = str(config.get("scE2G_dir") or "")
    base_cfg_path = os.path.join(scE2G_dir, "config", "config.yaml")
    base = load_config(base_cfg_path) if os.path.exists(base_cfg_path) else {}
    base = make_paths_absolute(base, scE2G_dir)
    base.update({"cell_clusters": cell_clusters_table, "results_dir": results_dir, "IGV_dir": results_dir})
    base.update(make_paths_absolute(config.get("scE2G_options") or {}, scE2G_dir))
    return base


def sce2g_command(config_path: str, dataset: str, dry_run: bool = True) -> "List[str]":
    """The snakemake invocation run_pipeline.py stage 2 performs (flags it always bakes in)."""
    cmd = ["snakemake", "-s", "workflow/Snakefile", "--configfile", config_path, "--use-conda", "--conda-prefix",
           "<conda-prefix>", "--rerun-triggers", "mtime", "--keep-going", "--jobs", "4", "--config",
           "pipeline_mode=default", "sce2g_modules=true"]
    return cmd + (["-n"] if dry_run else ["--executor", "slurm", "--profile", "slurm.smk9", "-p"])


# ---------------------------------------------------------------------------
# Portal-format reformatting (update_scE2G_pred_formats.R, reformat.smk, qc_stats.smk bgzip_index_bedpe)
# ---------------------------------------------------------------------------

FULL_MULTIOME_COLS = [
    ("ElementChr", "chr"), ("ElementStart", "start"), ("ElementEnd", "end"), ("ElementName", "name"),
    ("ElementClass", "class"), ("GeneSymbol", "TargetGene"), ("GeneEnsemblID", "TargetGeneEnsembl_ID"),
    ("GeneTSS", "TargetGeneTSS"), ("CellAnnotation", "CellAnnotation"), ("Score", "E2G.Score.qnorm"),
    ("distance", "distance"), ("feature.normalizedATAC_prom", "normalizedATAC_prom"),
    ("feature.numTSSEnhGene", "numTSSEnhGene"), ("feature.numNearbyEnhancers", "numNearbyEnhancers"),
    ("feature.ubiqExpressed", "ubiqExpressed"), ("feature.numCandidateEnhGene", "numCandidateEnhGene"),
    ("feature.ARC.E2G.Score", "ARC.E2G.Score"), ("isSelfPromoter", "isSelfPromoter"), ("Kendall", "Kendall"),
    ("normalizedATAC_enh", "normalizedATAC_enh"), ("RNA_pseudobulkTPM", "RNA_pseudobulkTPM"),
    ("RNA_meanLogNorm", "RNA_meanLogNorm"), ("RNA_percentCellsDetected", "RNA_percentCellsDetected"),
    ("ABC.Score", "ABC.Score")]
FULL_SCATAC_COLS = FULL_MULTIOME_COLS[:16] + [("feature.ABC.Score", "ABC.Score"), ("isSelfPromoter", "isSelfPromoter"),
                                               ("normalizedATAC_enh", "normalizedATAC_enh")]
THRESHOLDED_COLS = [("ElementChr", "chr"), ("ElementStart", "start"), ("ElementEnd", "end"), ("ElementName", "_name2"),
                    ("ElementClass", "class"), ("GeneSymbol", "TargetGene"), ("GeneEnsemblID", "TargetGeneEnsembl_ID"),
                    ("GeneTSS", "TargetGeneTSS"), ("CellAnnotation", "CellAnnotation"), ("Score", "E2G.Score.qnorm"),
                    ("isSelfPromoter", "isSelfPromoter")]
CANDIDATE_COLS = [("ElementChr", "chr"), ("ElementStart", "start"), ("ElementEnd", "end"), ("ElementName", "_name2"),
                  ("ElementClass", "class"), ("GeneTSS", "TargetGeneTSS"), ("GeneSymbol", "TargetGene"),
                  ("GeneEnsemblID", "TargetGeneEnsembl_ID"), ("SampleSummaryShort", "SampleSummaryShort"),
                  ("E2G_Distance", "distance"), ("isSelfPromoter", "isSelfPromoter")]
FEATURE_LEAD = [("ElementChr", "chr"), ("ElementStart", "start"), ("ElementEnd", "end"), ("ElementName", "_name2"),
                ("ElementClass", "class"), ("GeneSymbol", "TargetGene"), ("GeneTSS", "TargetGeneTSS"),
                ("GeneEnsemblID", "TargetGeneEnsembl_ID"), ("isSelfPromoter", "isSelfPromoter"), ("CellType", "CellType"),
                ("E2G_Distance", "distance")]


def read_prediction_table(path):
    pd = _pd()
    return pd.read_csv(path, sep="\t", low_memory=False, na_values=["NA"], keep_default_na=True)


def _select(df, spec: "Sequence[Tuple[str, str]]"):
    pd = _pd()
    missing = [src for _, src in spec if src not in df.columns]
    if missing:
        raise ValueError(f"input is missing column(s) {missing}")
    return pd.DataFrame({out: df[src].to_numpy() for out, src in spec}, columns=[o for o, _ in spec])


def _coord_name(df):
    return df["chr"].astype(str) + ":" + df["start"].astype("Int64").astype(str) + "-" + df["end"].astype("Int64").astype(str)


def reformat_header(method: str, version: str, cell_type, term_id, summary, threshold=None, input_name: str = "",
                    portal_link=None, fmt: str = "portal", cell_type_values=None) -> "List[str]":
    def v(x):
        return "" if x is None else str(x)
    if fmt == "synapse-legacy":
        name = " ".join(cell_type_values) if cell_type_values else v(cell_type)
        h = ["# PRELIMINARY — NOT THE FINAL VERSION", "# Pending official cell type description (SampleSummaryShort)",
             f"# Source: {v(method)}", f"# Version: {v(version)}", f"# GenomeReference: {GENOME_REFERENCE}", f"# URL: {SCE2G_URL}",
             "# Assays: 10x Multiome", "# SampleAgnostic: False", f"# SampleTermName: {name}", f"# SampleTermID: {v(term_id)}",
             f"# SampleSummaryShort: {v(summary)}"]
    else:
        h = [f"# Source: {v(method)}", f"# Version: {v(version)}", f"# GenomeReference: {GENOME_REFERENCE}", f"# URL: {SCE2G_URL}",
             "# Assays: 10x Multiome", "# SampleAgnostic: False", f"# SampleTermName: {v(cell_type)}",
             f"# SampleTermID: {v(term_id)}", f"# CellAnnotation: {v(summary)}"]
    if threshold is not None:
        h.append(f"# ScoreThreshold: {threshold}")
    if "_list" not in input_name:
        h.append("# ScoreType: positive_score")
    if portal_link is not None:
        h += (["# Metadata:", str(portal_link)] if fmt == "synapse-legacy"
              else [f"# Metadata: {PORTAL_TABULAR_URL}{portal_link}"])
    return h


def reformat_kind(input_name: str, method: str, threshold) -> str:
    if "element_list" in input_name:
        return "element_list"
    if "gene_list" in input_name:
        return "gene_list"
    atac = "ATAC" in method
    if threshold is None:
        return "scatac_full" if atac else "multiome_full"
    return "scatac_thresholded" if atac else "multiome_thresholded"


def reformat_table(pred, kind: str, summary=None, all_columns: bool = False, fmt: str = "portal", catlas: bool = False):
    """The column transformation of update_scE2G_pred_formats.R for one input table."""
    _pd()
    df = pred.copy()
    ann_col = "SampleSummaryShort" if fmt == "synapse-legacy" else "CellAnnotation"
    if summary is not None:
        df[ann_col] = summary
    if kind == "element_list":
        if all_columns:
            df = df.rename(columns={"chr": "ElementChr", "start": "ElementStart", "end": "ElementEnd", "class": "ElementClass"})
            df["ElementName"] = df["ElementChr"].astype(str) + ":" + df["ElementStart"].astype(str) + "-" + df["ElementEnd"].astype(str)
            lead = ["ElementChr", "ElementStart", "ElementEnd", "ElementName", "ElementClass"]
            return df[lead + [c for c in df.columns if c not in lead]]
        df["_name2"] = _coord_name(df)
        return _select(df, [("ElementChr", "chr"), ("ElementStart", "start"), ("ElementEnd", "end"),
                            ("ElementName", "_name2"), ("ElementClass", "class")])
    if kind == "gene_list":
        df["TSSStart"] = df["tss"] - 250
        df["TSSEnd"] = df["tss"] + 250
        if all_columns:
            df = df.rename(columns={"chr": "TSSChr", "tss": "TSS", "name": "GeneSymbol", "Ensembl_ID": "GeneEnsemblID",
                                    "strand": "GeneStrand"})
            lead = ["TSSChr", "TSSStart", "TSSEnd", "GeneSymbol", "GeneEnsemblID", "GeneStrand"]
            return df[lead + [c for c in df.columns if c not in lead]]
        return _select(df, [("TSSChr", "chr"), ("TSSStart", "TSSStart"), ("TSSEnd", "TSSEnd"), ("GeneSymbol", "name"),
                            ("GeneEnsemblID", "Ensembl_ID"), ("GeneStrand", "strand")])
    df["_name2"] = _coord_name(df)
    if fmt == "synapse-legacy":
        spec = [(o, ann_col if o == "CellAnnotation" else s) for o, s in THRESHOLDED_COLS]
        spec = [("SampleSummaryShort", s) if o == "CellAnnotation" else (o, s) for o, s in spec]
        if not kind.startswith("scatac"):
            spec.insert(10, ("Score.ignoreTPM", "E2G.Score.qnorm.ignoreTPM"))
        return _select(df, spec)
    if kind == "scatac_full":
        spec = [c for c in FULL_SCATAC_COLS if not (catlas and c[0] == "normalizedATAC_enh")]
        return _select(df, spec)
    if kind == "multiome_full":
        return _select(df, FULL_MULTIOME_COLS)
    return _select(df, THRESHOLDED_COLS)


def write_table_with_header(df, header_lines: "Sequence[str]", out, column_names: bool = True) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    for h in header_lines:
        buf.write(h + "\n")
    if column_names:
        buf.write("\t".join(df.columns) + "\n")
    cols = [df[c].tolist() for c in df.columns]
    for row in zip(*cols) if cols else []:
        buf.write("\t".join(fmt_value(None if (isinstance(x, float) and x != x) else x) for x in row) + "\n")
    data = buf.getvalue().encode()
    if str(out).endswith(".gz"):
        with gzip.open(out, "wb") as fh:
            fh.write(data)
    else:
        out.write_bytes(data)
    return out


def reformat_file(input_path, out, method: str, version: str, cell_type=None, term_id=None, summary=None,
                  threshold=None, portal_link=None, all_columns: bool = False, fmt: str = "portal",
                  catlas: bool = False) -> dict:
    pred = read_prediction_table(input_path)
    name = os.path.basename(str(input_path))
    kind = reformat_kind(name, method, threshold)
    ctv = sorted(pred["CellType"].dropna().astype(str).unique()) if (fmt == "synapse-legacy" and "CellType" in pred.columns) else None
    header = reformat_header(method, version, cell_type, term_id, summary, threshold, name, portal_link, fmt, ctv)
    table = reformat_table(pred, kind, summary, all_columns, fmt, catlas)
    write_table_with_header(table, header, out)
    return {"out": str(out), "kind": kind, "rows": int(len(table)), "columns": list(table.columns), "header": header}


def element_bed_header(version, cell_type, term_id, summary, portal_link) -> "List[str]":
    return ["# Source: ABC (originally EnhancerList.bed)", f"# Version: {version}", f"# GenomeReference: {GENOME_REFERENCE}",
            f"# URL: {SCE2G_URL}", "# Assays: 10x Multiome", "# SampleAgnostic: False", f"# SampleTermName: {cell_type}",
            f"# SampleTermID: {term_id}", f"# CellAnnotation: {summary}", f"# Metadata: {PORTAL_TABULAR_URL}{portal_link}",
            "#ElementChr\tElementStart\tElementEnd\tElementName"]


def write_element_bed(enhancer_list, out_gz, version, cell_type, term_id, summary, portal_link) -> Path:
    """reformat_element_list + bgzip_index_element_list: header block + EnhancerList.bed verbatim, bgzip, tabix."""
    out_gz = str(out_gz)
    plain = Path(out_gz[:-3] if out_gz.endswith(".gz") else out_gz)
    plain.parent.mkdir(parents=True, exist_ok=True)
    with open(plain, "w") as fo:
        fo.write("\n".join(element_bed_header(version, cell_type, term_id, summary, portal_link)) + "\n")
        with _open_text(enhancer_list) as fi:
            shutil.copyfileobj(fi, fo)
    return bgzip_tabix(plain, preset="bed")


def bgzip_index_bedpe(bedpe, out_gz) -> Path:
    """sort -k1,1 -k2,2n (C locale, whole-line last resort) | bgzip; tabix -p bed."""
    with _open_text(bedpe) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]

    def key(ln):
        f = ln.split("\t")
        try:
            s = float(f[1])
        except (IndexError, ValueError):
            s = 0.0
        return (f[0].encode(), s, ln.encode())
    lines.sort(key=key)
    out_gz = str(out_gz)
    plain = Path(out_gz[:-3] if out_gz.endswith(".gz") else out_gz)
    plain.parent.mkdir(parents=True, exist_ok=True)
    plain.write_text("".join(ln + "\n" for ln in lines))
    return bgzip_tabix(plain, preset="bed")


def candidate_pairs(pred, summary: str):
    df = pred.copy()
    df["SampleSummaryShort"] = summary
    df["_name2"] = _coord_name(df)
    return _select(df, CANDIDATE_COLS)


def feature_table(pred):
    df = pred.copy()
    df["_name2"] = _coord_name(df)
    lead = _select(df, FEATURE_LEAD)
    # dplyr: mutate(ElementChr = chr) copies chr, so `chr` itself is never selected and survives in everything()
    used = {s for _, s in FEATURE_LEAD if s not in ("_name2", "chr")} | {"name"}
    rest = [c for c in pred.columns if c not in used]
    for c in rest:
        lead[c] = pred[c].to_numpy()
    return lead


def feature_header(model, cell_type, term_id, summary) -> "List[str]":
    def v(x):
        return "" if x is None else str(x)
    return [f"# Source: scE2G {model}", f"# GenomeReference: {GENOME_REFERENCE}", f"# URL: {SCE2G_URL}", "# Assays: 10x Multiome",
            "# SampleAgnostic: False", f"# SampleTermName: {v(cell_type)}", f"# SampleTermID: {v(term_id)}",
            f"# SampleSummaryShort: {v(summary)}"]


# ---------------------------------------------------------------------------
# QC aggregation (aggregate_qc_stats.py + plot_all_qc_stats_from_merged.R)
# ---------------------------------------------------------------------------

QC_EXCLUDED_DIRS = {"qc_plots", "config", "tmp"}


def cell_count_from_metrics(plots_dir, cluster) -> "Optional[int]":
    p = os.path.join(str(plots_dir), cluster, METRICS_NAME)
    if not os.path.exists(p):
        return None
    return sum(int(float(r["n_cells"])) for r in read_rows_tsv(p))


def aggregate_qc_rows(dataset_dir, plots_dir) -> "List[dict]":
    latest = {}
    if not os.path.isdir(str(dataset_dir)):
        return []
    for cl in os.listdir(str(dataset_dir)):
        if cl in QC_EXCLUDED_DIRS or not os.path.isdir(os.path.join(str(dataset_dir), cl)):
            continue
        for path in glob.glob(os.path.join(str(dataset_dir), cl, "*", "scE2G_predictions_threshold*_stats.tsv")):
            mtime = os.path.getmtime(path)
            rows = read_rows_tsv(path)
            if not rows:
                continue
            row = rows[0]
            if int(float(row.get("cell_count", 0) or 0)) == 0:
                patched = cell_count_from_metrics(plots_dir, row["cluster"])
                if patched:
                    row["cell_count"] = str(patched)
            key = (row["cluster"], row["model_name"])
            if key not in latest or mtime > latest[key][0]:
                latest[key] = (mtime, row)
    return sorted((r for _, r in latest.values()), key=lambda r: (r["cluster"], r["model_name"]))


def qc_plots(all_stats_path, out: Path) -> "Tuple[List[Path], List[str]]":
    pd = _pd()
    try:
        from sce2g_pipeline_skill import qc_figures, qc_summary  # type: ignore
    except Exception:
        return [], []
    stats = pd.read_csv(all_stats_path, sep="\t")
    for c in ("fragments_total", "cell_count", "umi_count"):
        stats[c] = pd.to_numeric(stats[c], errors="coerce").fillna(0)
    ds, warns, _ = qc_summary(stats)
    return qc_figures(stats, ds, None, out), warns


# ---------------------------------------------------------------------------
# Stale reformatted files (stale_reformats.py)
# ---------------------------------------------------------------------------

def embedded_header(path) -> "Optional[dict]":
    found = {}
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                for key in STALE_HEADER_FIELDS:
                    if line.startswith(key):
                        found[key] = line[len(key):].strip()
    except OSError:
        return None
    return found


def reformat_files(cluster_dir, dataset, cluster) -> "List[str]":
    if not os.path.isdir(str(cluster_dir)):
        return []
    prefix = f"{dataset}_{cluster}_"
    return [os.path.join(str(cluster_dir), n) for n in sorted(os.listdir(str(cluster_dir)))
            if n.startswith(prefix) and n.endswith(REFORMAT_SUFFIXES)]


def find_stale(results_dir, dataset, annotation_rows) -> "List[Tuple[str, dict]]":
    stale = []
    for row in annotation_rows:
        if row.get("dataset") != dataset or not row.get("cell_annotation"):
            continue
        want = {k: (row.get(col) or "") for k, col in STALE_HEADER_FIELDS.items()}
        for path in reformat_files(os.path.join(str(results_dir), dataset, row["cluster"]), dataset, row["cluster"]):
            got = embedded_header(path)
            if not got:
                continue
            diffs = {k: (v, want[k]) for k, v in got.items() if v != want[k]}
            if diffs:
                stale.append((path, diffs))
    return stale


def remove_stale(stale) -> "List[str]":
    removed = []
    for path, _ in stale:
        for p in (path, path + ".tbi"):
            if os.path.exists(p):
                os.remove(p)
                removed.append(p)
    return removed


# ---------------------------------------------------------------------------
# SQLite ledger + Cell Annotation cache (igvf_metadata/state.py)
# ---------------------------------------------------------------------------

STATE_SCHEMA = """
CREATE TABLE IF NOT EXISTS clusters (dataset TEXT NOT NULL, cluster TEXT NOT NULL, excluded INTEGER NOT NULL DEFAULT 0,
    exclusion_reason TEXT, PRIMARY KEY (dataset, cluster));
CREATE TABLE IF NOT EXISTS uploads (id INTEGER PRIMARY KEY AUTOINCREMENT, dataset TEXT NOT NULL, cluster TEXT NOT NULL,
    model TEXT NOT NULL DEFAULT '', table_name TEXT NOT NULL, variant TEXT NOT NULL DEFAULT '', alias TEXT NOT NULL UNIQUE,
    payload_hash TEXT NOT NULL, portal_id TEXT, status TEXT NOT NULL DEFAULT 'pending', attempt_count INTEGER NOT NULL DEFAULT 0,
    last_attempt_at TEXT, last_error TEXT, created_at TEXT NOT NULL, updated_at TEXT NOT NULL,
    UNIQUE (dataset, cluster, model, table_name, variant));
CREATE TABLE IF NOT EXISTS cell_annotations (dataset TEXT NOT NULL, cluster TEXT NOT NULL, cell_annotation TEXT NOT NULL,
    cl_id TEXT NOT NULL, term_id TEXT NOT NULL, term_name TEXT NOT NULL, cell_qualifier TEXT, portal_samples TEXT NOT NULL,
    all_primary_released INTEGER NOT NULL, principal_uploaded INTEGER NOT NULL, principal_alias TEXT, fetched_at TEXT NOT NULL,
    UNIQUE (dataset, cluster));
CREATE TABLE IF NOT EXISTS cell_annotations_fetch_log (id INTEGER PRIMARY KEY CHECK (id = 1), fetched_at TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS cell_metadata_primary_pseudobulks (alias TEXT PRIMARY KEY, subsample TEXT, cell_annotation TEXT,
    cl_id TEXT, term_id TEXT, term_name TEXT, cell_qualifier TEXT, status TEXT, fetched_at TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS cell_metadata_principal_pseudobulks (alias TEXT PRIMARY KEY, cell_annotation TEXT, fetched_at TEXT NOT NULL);
"""
CELL_ANNOTATION_COLUMNS = ["dataset", "cluster", "cell_annotation", "cl_id", "term_id", "term_name", "cell_qualifier",
                           "portal_samples", "all_primary_released", "principal_uploaded", "principal_alias"]
CACHE_TTL_HOURS = 24


class State:
    def __init__(self, path):
        path = str(path)
        if os.path.dirname(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
        self.conn = sqlite3.connect(path, timeout=30)
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(STATE_SCHEMA)
        self.conn.commit()

    def close(self):
        self.conn.close()

    def q(self, sql, args=()):
        return [dict(r) for r in self.conn.execute(sql, args).fetchall()]

    def x(self, sql, args=()):
        self.conn.execute(sql, args)
        self.conn.commit()

    # uploads ledger
    def get_upload(self, dataset, cluster, model, table, variant):
        r = self.q("SELECT * FROM uploads WHERE dataset=? AND cluster=? AND model=? AND table_name=? AND variant=?",
                   (dataset, cluster, model or "", table, variant or ""))
        return r[0] if r else None

    def claim_pending(self, dataset, cluster, model, table, variant, alias, h):
        now = _now()
        self.x("""INSERT INTO uploads (dataset, cluster, model, table_name, variant, alias, payload_hash, status, attempt_count,
                  created_at, updated_at) VALUES (?,?,?,?,?,?,?, 'pending', 0, ?, ?)
                  ON CONFLICT(dataset, cluster, model, table_name, variant) DO UPDATE SET payload_hash=excluded.payload_hash,
                  updated_at=excluded.updated_at, status = CASE WHEN status='uploaded' AND payload_hash=excluded.payload_hash
                  THEN status ELSE 'pending' END""", (dataset, cluster, model or "", table, variant or "", alias, h, now, now))
        return self.get_upload(dataset, cluster, model, table, variant)

    def record_result(self, row_id, status, portal_id=None, error=None):
        now = _now()
        self.x("""UPDATE uploads SET status=?, portal_id=COALESCE(?, portal_id), last_error=?, attempt_count=attempt_count+1,
                  last_attempt_at=?, updated_at=? WHERE id=?""", (status, portal_id, error, now, now, row_id))

    def mark_excluded(self, dataset, cluster, reason):
        self.x("""INSERT INTO clusters (dataset, cluster, excluded, exclusion_reason) VALUES (?,?,1,?)
                  ON CONFLICT(dataset, cluster) DO UPDATE SET excluded=1, exclusion_reason=?""", (dataset, cluster, reason, reason))

    def excluded_clusters(self, dataset):
        return {r["cluster"] for r in self.q("SELECT cluster FROM clusters WHERE dataset=? AND excluded=1", (dataset,))}

    # cell metadata cache
    def latest_fetch(self):
        r = self.q("SELECT fetched_at FROM cell_annotations_fetch_log WHERE id=1")
        return r[0]["fetched_at"] if r else None

    def record_fetch(self, now):
        self.x("INSERT INTO cell_annotations_fetch_log (id, fetched_at) VALUES (1, ?) ON CONFLICT(id) DO UPDATE SET fetched_at=excluded.fetched_at", (now,))

    def upsert_primary(self, alias, subsample, ann, cl_id, term_id, term_name, qual, status, now):
        self.x("""INSERT INTO cell_metadata_primary_pseudobulks VALUES (?,?,?,?,?,?,?,?,?) ON CONFLICT(alias) DO UPDATE SET
                  subsample=excluded.subsample, cell_annotation=excluded.cell_annotation, cl_id=excluded.cl_id, term_id=excluded.term_id,
                  term_name=excluded.term_name, cell_qualifier=excluded.cell_qualifier, status=excluded.status,
                  fetched_at=excluded.fetched_at""", (alias, subsample, ann, cl_id, term_id, term_name, qual, status, now))

    def upsert_principal(self, alias, ann, now):
        self.x("""INSERT INTO cell_metadata_principal_pseudobulks VALUES (?,?,?) ON CONFLICT(alias) DO UPDATE SET
                  cell_annotation=excluded.cell_annotation, fetched_at=excluded.fetched_at""", (alias, ann, now))

    def upsert_cell_annotation(self, row: dict, now):
        vals = [row.get(c) for c in CELL_ANNOTATION_COLUMNS]
        vals[8], vals[9] = int(bool(vals[8])), int(bool(vals[9]))
        self.x(f"""INSERT INTO cell_annotations ({', '.join(CELL_ANNOTATION_COLUMNS)}, fetched_at) VALUES ({','.join('?' * 12)})
                  ON CONFLICT(dataset, cluster) DO UPDATE SET {', '.join(f'{c}=excluded.{c}' for c in CELL_ANNOTATION_COLUMNS[2:])},
                  fetched_at=excluded.fetched_at""", (*vals, now))

    def get_cell_annotation(self, dataset, cluster):
        r = self.q("SELECT * FROM cell_annotations WHERE dataset=? AND cluster=?", (dataset, cluster))
        return r[0] if r else None

    def all_cell_annotations(self):
        return self.q("SELECT * FROM cell_annotations")


def payload_hash(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


# ---------------------------------------------------------------------------
# IGVF Portal reads: multireport, classification, Cell Annotation derivation (cell_metadata.py, pseudobulk_sets.py)
# ---------------------------------------------------------------------------

MULTIREPORT_QUERY = ("type=PseudobulkSet&status%21=deleted&limit=all&field=%40id&field=cell_annotation&field=aliases"
                     "&field=cell_type&field=cell_type.term_name&field=cell_type.term_id&field=summary&field=cell_qualifier"
                     "&field=input_file_sets&field=lab&field=samples&field=status")
PORTAL_FILES_FIELDS = ("@id", "accession", "aliases", "status", "file_set_type", "lab", "input_file_sets", "samples", "files",
                       "cell_annotation", "cell_qualifier", "cell_type")
WASHU_LAB = "/labs/yang-li/"
WASHU_FIELDS = ["@id", "aliases", "cell_annotation", "cell_qualifier", "cell_type", "cell_type.term_name",
                "cell_type.term_id", "status", "input_file_sets", "files", "files.content_type", "files.aliases",
                "files.status", "files.href", "files.md5sum"]


def portal_get(path: str):
    """Authenticated read through raw_data_pipeline.portal_json (IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY)."""
    try:
        from raw_data_pipeline import portal_json  # type: ignore
    except Exception:  # pragma: no cover
        from igvfagent.raw_data_pipeline import portal_json  # type: ignore
    return portal_json(path)


def fetch_multireport(query: str, json_path: "Optional[str]" = None) -> "List[dict]":
    if json_path:
        data = json.loads(Path(json_path).read_text())
        return data.get("@graph", []) if isinstance(data, dict) else list(data)
    status, data = portal_get(f"/multireport/?{query}")
    if status != 200 or not data:
        raise SystemExit(f"multireport GET failed (HTTP {status}); pass a saved --multireport-json to run offline")
    return data.get("@graph", [])


def input_file_set_ids(row) -> "List[str]":
    return [e.get("@id", "") for e in row.get("input_file_sets") or [] if isinstance(e, dict)]


def classify_pseudobulk(row) -> "Optional[str]":
    ids = input_file_set_ids(row)
    if not ids:
        return None
    if all(i.startswith("/analysis-sets/") for i in ids):
        return "primary"
    if all(i.startswith("/pseudobulk-sets/") for i in ids):
        return "principal"
    return None


def input_collections(row) -> "List[str]":
    return sorted({i.strip("/").split("/")[0] for i in input_file_set_ids(row) if i})


def principal_analysis_set(row) -> "Optional[str]":
    accs = {e.get("accession") for e in row.get("input_file_sets") or []
            if isinstance(e, dict) and e.get("file_set_type") == "principal analysis" and e.get("accession")}
    return next(iter(accs)) if len(accs) == 1 else None


def cl_id_of(ct) -> "Optional[str]":
    if not isinstance(ct, dict):
        return None
    return (ct.get("@id") or "").strip("/").rsplit("/", 1)[-1] or None


def term_id_of(ct) -> "Optional[str]":
    return (ct.get("term_id") or None) if isinstance(ct, dict) else None


def term_name_of(ct) -> "Optional[str]":
    return (ct.get("term_name") or None) if isinstance(ct, dict) else None


def alias_suffix(alias):
    return alias.split(":", 1)[-1] if alias and ":" in alias else alias


def is_stale(last_fetch, ttl_hours: float = CACHE_TTL_HOURS) -> bool:
    if not last_fetch:
        return True
    return (datetime.now(timezone.utc) - datetime.fromisoformat(last_fetch)).total_seconds() > ttl_hours * 3600


def cache_multireport(state: State, rows: "List[dict]") -> dict:
    """fetch_if_stale's saving half: every primary (exactly one sample, with alias) and principal row."""
    now = _now()
    state.record_fetch(now)
    n_pri = n_prin = no_alias = 0
    for row in rows:
        kind = classify_pseudobulk(row)
        alias = (row.get("aliases") or [None])[0]
        if kind == "principal":
            if row.get("cell_annotation") and alias:
                state.upsert_principal(alias, row["cell_annotation"], now)
                n_prin += 1
            continue
        if kind != "primary":
            continue
        samples = [s.get("accession") for s in row.get("samples") or [] if isinstance(s, dict)]
        if len(samples) != 1:
            continue
        if not alias:
            no_alias += 1
            continue
        ct = row.get("cell_type")
        state.upsert_primary(alias, samples[0], row.get("cell_annotation"), cl_id_of(ct), term_id_of(ct), term_name_of(ct),
                             row.get("cell_qualifier"), row.get("status"), now)
        n_pri += 1
    return {"primary_saved": n_pri, "principal_saved": n_prin, "skipped_no_alias": no_alias, "fetched_at": now}


def subsamples_by_frequency(qc_guide) -> "List[str]":
    counts = Counter(r["subsample"] for r in read_guide_rows(qc_guide))
    return [s for s, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]


def unique_subsamples(qc_guide) -> "List[str]":
    return sorted({r["subsample"] for r in read_guide_rows(qc_guide)})


def merge_sources(cluster_cfg: dict, cluster: str) -> "List[str]":
    raw = str(cluster_cfg.get("pseudobulk_annotation") or "")
    if "," not in raw:
        return [cluster]
    names = [p.strip() for p in raw.split(",") if p.strip()]
    return names or [cluster]


def derive_scopes(state: State, cluster_configs: "Dict[Tuple[str, str], dict]") -> "List[dict]":
    now = _now()
    statuses, local = [], {}
    for (d, c), cfg in sorted(cluster_configs.items()):
        try:
            local[(d, c)] = subsamples_by_frequency(cfg["qc_guide"])
        except (OSError, KeyError) as e:
            statuses.append({"dataset": d, "cluster": c, "resolved": False, "cell_annotation": None,
                             "reason": f"unreadable_qc_guide: {type(e).__name__}: {e}"})
    primaries = {alias_suffix(r["alias"]): r for r in state.q("SELECT * FROM cell_metadata_primary_pseudobulks")}
    principals = {r["cell_annotation"]: r["alias"] for r in state.q("SELECT * FROM cell_metadata_principal_pseudobulks")
                  if r["cell_annotation"] and r["alias"]}
    for (d, c), subs in local.items():
        names = merge_sources(cluster_configs[(d, c)], c)
        rows_by_sub, missing, scope_rows = {}, [], []
        for s in subs:
            hits = [f"{d}-{n}-{s}" for n in names if f"{d}-{n}-{s}" in primaries]
            if hits:
                rows_by_sub[s] = [primaries[h] for h in hits]
                scope_rows.extend(primaries[h] for h in hits)
            else:
                missing.append(s)
        if missing:
            statuses.append({"dataset": d, "cluster": c, "resolved": False, "cell_annotation": None,
                             "reason": f"no_matching_primary_alias: {','.join(sorted(missing))}"})
            continue
        triples = {(r["cl_id"], r["term_id"], r["term_name"]) for r in scope_rows}
        if len(triples) != 1:
            statuses.append({"dataset": d, "cluster": c, "resolved": False, "cell_annotation": None,
                             "reason": f"cell_type_disagreement: {len(triples)} distinct (cl_id, term_id, term_name)"})
            continue
        cl_id, term_id, term_name = next(iter(triples))
        if len(names) > 1:
            ann, qual = term_name, None
        else:
            win = rows_by_sub[subs[0]][0]
            ann, qual = win["cell_annotation"], win["cell_qualifier"]
        blank = [n for n, v in (("cell_annotation", ann), ("cl_id", cl_id), ("term_id", term_id), ("term_name", term_name)) if not v]
        if blank:
            statuses.append({"dataset": d, "cluster": c, "resolved": False, "cell_annotation": None,
                             "reason": f"blank_required_field: {','.join(blank)}"})
            continue
        palias = principals.get(ann)
        state.upsert_cell_annotation({"dataset": d, "cluster": c, "cell_annotation": ann, "cl_id": cl_id, "term_id": term_id,
                                      "term_name": term_name, "cell_qualifier": qual, "portal_samples": ",".join(sorted(subs)),
                                      "all_primary_released": all(r["status"] == "released" for r in scope_rows),
                                      "principal_uploaded": palias is not None, "principal_alias": palias}, now)
        statuses.append({"dataset": d, "cluster": c, "resolved": True, "reason": "ok", "cell_annotation": ann})
    return statuses


def seed_catlas(state: State, rows: "List[dict]", dataset: str = "catlas") -> "List[dict]":
    """seed_catlas_cell_annotations.py: WashU clusters keyed by the last '-' token of the set alias."""
    now, seeded = _now(), []
    for row in rows:
        aliases = row.get("aliases") or []
        if not aliases:
            continue
        cluster = aliases[0].rsplit("-", 1)[-1]
        ct = row.get("cell_type")
        rec = {"dataset": dataset, "cluster": cluster, "cell_annotation": row.get("cell_annotation", ""),
               "cl_id": cl_id_of(ct) or "TODO: ontology_id (not returned by portal)",
               "term_id": term_id_of(ct) or "TODO: ontology_id (not returned by portal)",
               "term_name": term_name_of(ct) or row.get("cell_annotation", ""), "cell_qualifier": row.get("cell_qualifier") or None,
               "portal_samples": row.get("@id"), "all_primary_released": row.get("status") == "released",
               "principal_uploaded": False, "principal_alias": None}
        state.upsert_cell_annotation(rec, now)
        seeded.append(rec)
    return seeded


def shareable_rows(state: State) -> "List[dict]":
    return [r for r in state.all_cell_annotations() if r["all_primary_released"] or r["principal_uploaded"]]


# Snapshot (cell_annotation_snapshot.py)

def cluster_set_digest(keys) -> str:
    return hashlib.sha256(";".join(sorted(f"{d}/{c}" for d, c in keys)).encode()).hexdigest()[:16]


def snapshot_path(output_dir, dataset) -> Path:
    return Path(output_dir) / "igvf_metadata" / f"{dataset}_cell_annotations.tsv"


def write_snapshot(path, rows, fetched_at, digest) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write(f"#portal_fetched_at\t{fetched_at}\n#derived_at\t{_now()}\n#cluster_set_digest\t{digest}\n")
        w = csv.DictWriter(fh, fieldnames=CELL_ANNOTATION_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["dataset"], r["cluster"])):
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in CELL_ANNOTATION_COLUMNS})
    return path


class SnapshotError(RuntimeError):
    pass


def read_snapshot(path, expected_digest=None, max_age_hours: "Optional[float]" = 24) -> dict:
    if not os.path.exists(path):
        raise SnapshotError(f"no CellAnnotation snapshot at {path}")
    prov, data = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                k, _, v = line[1:].rstrip("\n").partition("\t")
                prov[k] = v
            else:
                data.append(line)
    if expected_digest is not None and prov.get("cluster_set_digest") != expected_digest:
        raise SnapshotError(f"{path} was built for a different cluster set")
    if max_age_hours is not None:
        fa = prov.get("portal_fetched_at")
        if not fa:
            raise SnapshotError(f"{path} records no portal_fetched_at")
        age = (datetime.now(timezone.utc) - datetime.fromisoformat(fa)).total_seconds() / 3600
        if age > max_age_hours:
            raise SnapshotError(f"{path} is built on a portal fetch {age:.1f}h old (limit {max_age_hours}h)")
    return {(r["dataset"], r["cluster"]): r for r in csv.DictReader(data, delimiter="\t")}


# ---------------------------------------------------------------------------
# Standalone reports (igvf_cell_annotation_report/*)
# ---------------------------------------------------------------------------

REPORT_JOIN = " | "
ANNOTATION_REPORT_COLUMNS = ["Dataset_Cluster", "Subsamples", "N_Subsamples", "CellAnnotation", "SampleTermID",
                             "SampleTermName", "CellQualifier", "Status", "QCGuideFile"]


def collect_primary_groups(rows) -> "Tuple[Dict[str, List[dict]], dict]":
    groups, skipped = defaultdict(list), Counter()
    for row in rows:
        if classify_pseudobulk(row) != "primary":
            continue
        samples = [s.get("accession") for s in row.get("samples") or [] if isinstance(s, dict)]
        if len(samples) != 1:
            skipped["multi_sample"] += 1
            continue
        alias = (row.get("aliases") or [None])[0]
        if not alias:
            skipped["no_alias"] += 1
            continue
        suf, tail = alias_suffix(alias), f"-{samples[0]}"
        if not suf.endswith(tail):
            skipped["unparseable_alias"] += 1
            continue
        ct = row.get("cell_type")
        groups[suf[: -len(tail)]].append({"subsample": samples[0], "cell_annotation": row.get("cell_annotation"),
                                          "term_id": term_id_of(ct), "term_name": term_name_of(ct),
                                          "cell_qualifier": row.get("cell_qualifier"), "status": row.get("status")})
    return dict(groups), dict(skipped)


def _distinct(vals) -> "List[str]":
    return sorted({v for v in vals if v})


def _group_status(members) -> str:
    return "released" if all(m["status"] == "released" for m in members) else "in progress"


def fallback_annotation_row(dc, members) -> dict:
    return {"Dataset_Cluster": dc, "Subsamples": REPORT_JOIN.join(sorted({m["subsample"] for m in members})),
            "N_Subsamples": len({m["subsample"] for m in members}),
            "CellAnnotation": REPORT_JOIN.join(_distinct(m["cell_annotation"] for m in members)),
            "SampleTermID": REPORT_JOIN.join(_distinct(m["term_id"] for m in members)),
            "SampleTermName": REPORT_JOIN.join(_distinct(m["term_name"] for m in members)),
            "CellQualifier": REPORT_JOIN.join(_distinct(m["cell_qualifier"] for m in members)),
            "Status": _group_status(members), "QCGuideFile": ""}


def find_qc_guide(qc_guide_dir, dataset, cluster) -> "Tuple[Optional[str], Optional[str]]":
    if dataset is None or cluster is None:
        return None, None
    d = os.path.join(str(qc_guide_dir), dataset, cluster)
    if not os.path.isdir(d):
        return None, None
    entries = set(os.listdir(d))
    if DEFAULT_QC_GUIDE_NAME in entries:
        return os.path.join(d, DEFAULT_QC_GUIDE_NAME), DEFAULT_QC_GUIDE_NAME
    cands = sorted(f for f in entries if f.endswith((".tsv", ".tsv.gz")))
    if len(cands) == 1:
        return os.path.join(d, cands[0]), cands[0]
    return None, None


def resolve_annotation_group(dc, members, qc_guide_dir) -> dict:
    dataset, cluster = dc.split("-", 1) if "-" in dc else (None, None)
    gp, gname = find_qc_guide(qc_guide_dir, dataset, cluster)
    if gp is None:
        return fallback_annotation_row(dc, members)
    counts = Counter(r["subsample"] for r in read_guide_rows(gp))
    ordered = [s for s, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])) if n >= 1]
    tids, tnames = _distinct(m["term_id"] for m in members), _distinct(m["term_name"] for m in members)
    if not ordered or len(tids) > 1 or len(tnames) > 1:
        return fallback_annotation_row(dc, members)
    win = next((m for m in members if m["subsample"] == ordered[0]), None)
    if win is None:
        return fallback_annotation_row(dc, members)
    return {"Dataset_Cluster": dc, "Subsamples": REPORT_JOIN.join(ordered), "N_Subsamples": len(ordered),
            "CellAnnotation": win["cell_annotation"] or "", "SampleTermID": tids[0] if tids else "",
            "SampleTermName": tnames[0] if tnames else "", "CellQualifier": win["cell_qualifier"] or "",
            "Status": _group_status(members), "QCGuideFile": gname}


def cell_annotation_report_rows(rows, qc_guide_dir=None) -> "Tuple[List[dict], dict]":
    groups, skipped = collect_primary_groups(rows)
    out = [(resolve_annotation_group(dc, m, qc_guide_dir) if qc_guide_dir else fallback_annotation_row(dc, m))
           for dc, m in sorted(groups.items())]
    return out, skipped


def dataset_accession_mapping(rows) -> "Tuple[Dict[str, List[str]], dict]":
    mapping, skipped = defaultdict(set), Counter()
    for row in rows:
        if classify_pseudobulk(row) != "primary":
            skipped["not_primary"] += 1
            continue
        alias = (row.get("aliases") or [None])[0]
        if not alias:
            skipped["no_alias"] += 1
            continue
        accs = {e.get("accession") for e in row.get("input_file_sets") or [] if isinstance(e, dict) and e.get("accession")}
        if len(accs) != 1:
            skipped["multi_input_file_set"] += 1
            continue
        suf = alias_suffix(alias)
        if "-" not in suf:
            skipped["unparseable_dataset"] += 1
            continue
        mapping[suf.split("-", 1)[0]].add(next(iter(accs)))
    return {d: sorted(a) for d, a in sorted(mapping.items())}, dict(skipped)


WASHU_HEADER = ["PseudobulkSet_ID", "Alias", "CellAnnotation", "CellQualifier", "SampleTermName", "SampleTermID",
                "SampleTermCLID", "Status", "FragmentsFile_Aliases", "FragmentsFile_Status", "FragmentsFile_Href",
                "FragmentsFile_MD5", "N_Fragments_Files"]


def washu_row(row) -> dict:
    ct = row.get("cell_type")
    frags = [f for f in row.get("files", []) or [] if f.get("content_type") == "fragments"]
    f = frags[0] if frags else {}
    return {"PseudobulkSet_ID": row.get("@id", ""), "Alias": (row.get("aliases") or [""])[0],
            "CellAnnotation": row.get("cell_annotation", ""), "CellQualifier": row.get("cell_qualifier") or "",
            "SampleTermName": term_name_of(ct) or "", "SampleTermID": term_id_of(ct) or "", "SampleTermCLID": cl_id_of(ct) or "",
            "Status": row.get("status", ""), "FragmentsFile_Aliases": ",".join(f.get("aliases", []) or []),
            "FragmentsFile_Status": f.get("status", ""), "FragmentsFile_Href": f.get("href", ""),
            "FragmentsFile_MD5": f.get("md5sum", ""), "N_Fragments_Files": len(frags)}


def accession_from_href(href: str) -> "Tuple[str, str]":
    fn = href.rstrip("/").rsplit("/", 1)[-1]
    return fn.split(".", 1)[0], fn


# ---------------------------------------------------------------------------
# Primary-pseudobulk file discovery + archive comparison (portal_files.py, downloader.py, compare_portal_vs_archive.py)
# ---------------------------------------------------------------------------

DEFAULT_FILE_LAB = "/labs/anshul-kundaje/"
TARGET_CONTENT_TYPES = {"fragments": "fragments.tsv.gz", "cell by gene matrix": "rna_counts_mtx.h5ad",
                        "per-cell quality report": "per_cell_qc.tsv.gz"}
EXCLUDED_FILE_STATUSES = frozenset({"revoked", "replaced"})
_ABS_SFN = re.compile(r"^.*/igvf/(?P<dataset>[^/]+)/pseudobulks/(?P<dirname>[^/]+)/(?P<filename>[^/]+)$")
_REL_SFN = re.compile(r"^(?:.*/)?pseudobulks/(?P<dirname>[^/]+)/(?P<filename>[^/]+)$")
_DIRNAME_SUBSAMPLE_ONLY = re.compile(r"^.*?-(?P<subsample>IGVFSM\w+)$")


def annotation_from_alias(alias, dataset, subsample) -> "Optional[str]":
    if not alias or not subsample:
        return None
    suf, tail = alias_suffix(alias), f"-{subsample}"
    if not suf.endswith(tail):
        return None
    suf = suf[: -len(tail)]
    if suf.startswith("pseudobulk-"):
        suf = suf[len("pseudobulk-"):]
    if dataset and suf.startswith(f"{dataset}-"):
        suf = suf[len(dataset) + 1:]
    return suf or None


def resolve_scope(row, file_obj) -> "Tuple[dict, List[str]]":
    sfn = file_obj.get("submitted_file_name") or ""
    reasons = []
    scope = {"submitted_file_name": sfn or None, "dataset": None, "annotation": None, "subsample": None,
             "dirname": None, "filename": None, "rel_path": None}
    if not sfn:
        return scope, ["no_submitted_file_name"]
    m = _ABS_SFN.match(sfn)
    if m:
        scope["dataset"] = m.group("dataset")
    else:
        m = _REL_SFN.match(sfn)
        if not m:
            return scope, ["unparseable_submitted_file_name"]
        reasons.append("relative_submitted_file_name")
        pas = principal_analysis_set(row)
        if pas:
            scope["dataset"] = pas
            reasons.append("dataset_from_principal_analysis_set")
        else:
            reasons.append("no_dataset")
    scope["dirname"], scope["filename"] = m.group("dirname"), m.group("filename")
    alias = (row.get("aliases") or [None])[0]
    dm = DIRNAME_RE.match(scope["dirname"])
    if dm:
        scope["annotation"], scope["subsample"] = dm.group("annotation"), dm.group("subsample")
    else:
        reasons.append("nonstandard_directory_name")
        sm = _DIRNAME_SUBSAMPLE_ONLY.match(scope["dirname"])
        if sm:
            scope["subsample"] = sm.group("subsample")
        else:
            samples = [x.get("accession") for x in row.get("samples") or [] if isinstance(x, dict)]
            if len(samples) == 1:
                scope["subsample"] = samples[0]
        ann = annotation_from_alias(alias, scope["dataset"], scope["subsample"])
        if ann:
            scope["annotation"] = ann
            reasons.append("annotation_from_alias")
        else:
            reasons.append("no_annotation")
    exp = TARGET_CONTENT_TYPES.get(file_obj.get("content_type"))
    if exp and scope["filename"] != exp:
        reasons.append(f"unexpected_filename:{scope['filename']}")
    if alias and scope["subsample"] and scope["subsample"] not in alias:
        reasons.append("alias_subsample_mismatch")
    if scope["dataset"] and scope["dirname"] and scope["filename"]:
        scope["rel_path"] = os.path.join(scope["dataset"], "pseudobulks", scope["dirname"], scope["filename"])
    return scope, reasons


def discover_portal_files(rows, lab=DEFAULT_FILE_LAB, content_types=None, datasets=None) -> "Tuple[List[dict], dict]":
    wanted = set(content_types or TARGET_CONTENT_TYPES)
    report = {"sets_total": len(rows), "sets_by_class": Counter(), "files_by_content_type": Counter(),
              "files_excluded_by_status": Counter(), "upload_status": Counter(), "sets_missing_content_type": Counter(),
              "review_reasons": Counter(), "sets_selected": 0, "sets_with_no_target_files": 0}
    records = []
    for row in rows:
        kind = classify_pseudobulk(row)
        report["sets_by_class"][kind or "unclassified"] += 1
        if kind != "primary":
            continue
        row_lab = (row.get("lab") or {}).get("@id") if isinstance(row.get("lab"), dict) else row.get("lab")
        if lab and row_lab != lab:
            continue
        seen, recs, excl = set(), [], []
        for f in row.get("files") or []:
            ct = f.get("content_type")
            if ct not in wanted:
                continue
            if f.get("status") in EXCLUDED_FILE_STATUSES:
                excl.append(f"{ct}:{f.get('status')}")
                continue
            seen.add(ct)
            scope, reasons = resolve_scope(row, f)
            recs.append({"accession": f.get("accession"), "file_set": row.get("@id"), "principal_analysis_set": principal_analysis_set(row),
                         "lab": row_lab, "content_type": ct, "file_format": f.get("file_format"), "href": f.get("href"),
                         "md5sum": f.get("md5sum"), "file_size": f.get("file_size"), "portal_status": f.get("status"),
                         "upload_status": f.get("upload_status"), **{k: scope[k] for k in ("submitted_file_name", "dataset",
                         "annotation", "subsample", "rel_path")}, "review_reasons": reasons, "set_alias": (row.get("aliases") or [None])[0]})
        if datasets is not None:
            recs = [r for r in recs if r["dataset"] in datasets]
            if not recs:
                continue
        if not recs:
            report["sets_with_no_target_files"] += 1
            continue
        missing = wanted - seen
        if missing:
            report["sets_missing_content_type"][",".join(sorted(missing))] += 1
            for r in recs:
                r["review_reasons"] = r["review_reasons"] + ["set_missing:" + ",".join(sorted(missing))]
        for r in recs:
            report["files_by_content_type"][r["content_type"]] += 1
            report["upload_status"][r["upload_status"]] += 1
            for reason in r["review_reasons"]:
                report["review_reasons"][reason.split(":", 1)[0]] += 1
        for t in excl:
            report["files_excluded_by_status"][t] += 1
        report["sets_selected"] += 1
        records.extend(recs)
    return records, {k: (dict(v) if isinstance(v, Counter) else v) for k, v in report.items()}


def needs_download(rec: dict, dest) -> "Tuple[bool, str]":
    """downloader.needs_download: re-fetch when absent, size differs or the portal md5 no longer matches."""
    if not os.path.exists(dest):
        return True, "absent"
    if rec.get("file_size") not in (None, "") and os.path.getsize(dest) != int(rec["file_size"]):
        return True, f"size {os.path.getsize(dest)} != portal {rec['file_size']}"
    if rec.get("md5sum") and md5_file(dest) != rec["md5sum"]:
        return True, "portal md5 changed since download"
    return False, "unchanged"


def redact(text) -> str:
    return re.sub(r"(https?://[^\s?]+)\?[^\s]*", r"\1?<redacted>", str(text)) if text else text


IDENTICAL, COMPRESSION_ONLY, MEANINGFUL, UNCERTAIN = "identical", "compression_only", "meaningful", "uncertain"


def archive_counterpart(archive_root, dataset, dirname, filename) -> "Optional[str]":
    base = os.path.join(str(archive_root), dataset, "pseudobulks", dirname)
    for cand in ([filename, filename[:-3]] if filename.endswith(".gz") else [filename, filename + ".gz"]):
        p = os.path.join(base, cand)
        if os.path.exists(p):
            return p
    return None


def _tsv_detail(a_path, b_path, barcode_col=None) -> dict:
    detail = {}
    with _open_text(a_path, "rb") as a, _open_text(b_path, "rb") as b:
        ha, hb = a.readline().decode("utf-8", "replace").rstrip("\n"), b.readline().decode("utf-8", "replace").rstrip("\n")
        detail["header_changed"] = ha != hb
        idx = None
        if barcode_col and ha == hb and barcode_col in ha.split("\t"):
            idx = ha.split("\t").index(barcode_col)
        na = nb = 0
        first, ba, bb = None, set(), set()
        while True:
            la, lb = a.readline(), b.readline()
            if not la and not lb:
                break
            if la:
                na += 1
                if idx is not None:
                    ba.add(la.decode().rstrip("\n").split("\t")[idx])
            if lb:
                nb += 1
                if idx is not None:
                    bb.add(lb.decode().rstrip("\n").split("\t")[idx])
            if first is None and la != lb:
                first = max(na, nb)
    detail.update({"rows_portal": na, "rows_archive": nb, "first_differing_line": first})
    if idx is not None:
        detail.update({"barcodes_only_portal": len(ba - bb), "barcodes_only_archive": len(bb - ba)})
    return detail


def _h5ad_detail(a_path, b_path) -> "Tuple[str, dict]":
    try:
        import h5py  # type: ignore
        np = _np()
    except ImportError:
        return UNCERTAIN, {"note": "h5py unavailable"}

    def summarise(path):
        out = {}
        with h5py.File(path, "r") as f:
            for key in ("obs", "var"):
                if key in f:
                    idx = f[key].attrs.get("_index")
                    idx = idx.decode() if isinstance(idx, bytes) else idx
                    if idx and idx in f[key]:
                        vals = f[key][idx][:]
                        h = hashlib.md5()
                        for v in vals:
                            h.update(v if isinstance(v, bytes) else str(v).encode())
                        out[f"n_{key}"], out[f"{key}_names_md5"] = len(vals), h.hexdigest()
            if "X" in f:
                x = f["X"]
                if isinstance(x, h5py.Group):
                    for part in ("data", "indices", "indptr"):
                        if part in x:
                            out[f"X_{part}_md5"] = hashlib.md5(np.ascontiguousarray(x[part][:]).tobytes()).hexdigest()
                else:
                    out["X_md5"] = hashlib.md5(np.ascontiguousarray(x[:]).tobytes()).hexdigest()
        return out
    try:
        pa, ar = summarise(a_path), summarise(b_path)
    except Exception as exc:  # noqa: BLE001
        return UNCERTAIN, {"note": f"h5ad introspection failed: {exc}"}
    diffs = {k for k in set(pa) | set(ar) if pa.get(k) != ar.get(k)}
    return (COMPRESSION_ONLY, {"note": "HDF5 bytes differ but obs/var names and X are identical"}) if not diffs else \
        (MEANINGFUL, {"h5ad_differences": ",".join(sorted(diffs))})


def compare_file(portal_path, archive_path, content_type, characterise: bool = True) -> "Tuple[str, str, dict]":
    ra, rb = md5_file(portal_path), md5_file(archive_path)
    if ra == rb:
        return "unchanged", IDENTICAL, {"raw_md5": ra}
    ca, cb = content_md5(portal_path), content_md5(archive_path)
    if ca == cb:
        return "unchanged", COMPRESSION_ONLY, {"content_md5": ca}
    det = {"content_md5_portal": ca, "content_md5_archive": cb}
    if not characterise:
        return "changed", UNCERTAIN, det
    if content_type == "cell by gene matrix":
        v, extra = _h5ad_detail(portal_path, archive_path)
        return ("unchanged" if v == COMPRESSION_ONLY else "changed"), v, {**det, **extra}
    try:
        return "changed", MEANINGFUL, {**det, **_tsv_detail(portal_path, archive_path,
                                                             "barcode" if content_type == "per-cell quality report" else None)}
    except Exception as exc:  # noqa: BLE001
        return "changed", UNCERTAIN, {**det, "note": str(exc)}


def cluster_decision(statuses: "set") -> str:
    if "not_downloaded" in statuses:
        return "INCOMPLETE_DOWNLOAD"
    if statuses == {"novel"}:
        return "NOVEL_needs_QC_filtering"
    if "novel" in statuses:
        return "PARTIALLY_NOVEL_needs_QC_filtering"
    if "changed" in statuses:
        return "CHANGED_needs_QC_filtering"
    if statuses == {"unchanged"}:
        return "SHAREABLE_unchanged"
    return "REVIEW"


# ---------------------------------------------------------------------------
# Prediction Set submitter_comment backfill (patch_prediction_set_submitter_comment.py)
# ---------------------------------------------------------------------------

VERSION_MARKER = "Version 1"
PREDICTION_SET_QUERY = ("type=PredictionSet&status%21=deleted&limit=all&field=%40id&field=aliases&field=submitter_comment"
                        "&field=uuid&field=accession")


def new_submitter_comment(current) -> "Optional[str]":
    if not current:
        return VERSION_MARKER
    if VERSION_MARKER in current:
        return None
    return f"({VERSION_MARKER}) - {current}"


def plan_submitter_comment(rows, datasets) -> "Tuple[List[dict], dict]":
    to_patch = []
    counts = Counter({"skipped-no-alias": 0, "skipped-wrong-dataset": 0, "unchanged": 0, "set-fresh": 0, "prepended": 0})
    for row in rows:
        alias = (row.get("aliases") or [None])[0]
        if alias is None:
            counts["skipped-no-alias"] += 1
            continue
        ds = alias.split(":", 1)[1].split("_", 1)[0] if ":" in alias else None
        if datasets and ds not in datasets:
            counts["skipped-wrong-dataset"] += 1
            continue
        cur = row.get("submitter_comment")
        new = new_submitter_comment(cur)
        if new is None:
            counts["unchanged"] += 1
            continue
        counts["set-fresh" if not cur else "prepended"] += 1
        to_patch.append({"record_id": row.get("uuid") or row.get("accession") or row.get("@id"), "alias": alias, "old": cur, "new": new})
    return to_patch, dict(counts)


# ---------------------------------------------------------------------------
# Coverage report (generate_report.py)
# ---------------------------------------------------------------------------

REPORT_HEADER = ["dataset", "cluster", "cell_count", "fragments_total", "umi_count", "quality_pass", "exclusion_reason",
                 "has_cell_annotation", "igvf_manifest_excluded", "predictions_generated", "reformatted", "status",
                 "manifest_status", "manifest_rows_expected", "manifest_rows_planned", "manifest_gap_reason"]
GAP_PRIORITY = ["invalid", "enabled-check-failed", "skipped-missing-file", "deferred"]
REAL_ROW_OUTCOMES = frozenset({"planned-post", "planned-patch", "unchanged"})
MANIFEST_GAP_OUTCOMES = frozenset({"skipped-missing-file", "invalid", "enabled-check-failed"})


def derive_status(reason, quality_pass, has_ann, has_pred, reformatted_) -> str:
    if reason in UNPROCESSABLE_REASONS:
        return "missing-qc-guide"
    if not quality_pass:
        return "excluded-quality"
    if not has_pred:
        return "quality-pass-not-yet-processed"
    if not has_ann:
        return "predictions-only-missing-annotation"
    if not reformatted_:
        return "predictions-generated-reformat-pending"
    return "fully-processed"


def summarize_manifest(rows, manifest_excluded, quality_pass, has_ann) -> "Tuple[str, Any, Any, str]":
    if manifest_excluded:
        return "excluded", "", "", "igvf_manifest_excluded"
    if not quality_pass:
        return "excluded", "", "", "failed_quality_gate"
    if not has_ann:
        return "blocked", "", "", "no_cell_annotation"
    if not rows:
        return "not-generated", "", "", ""
    considered = [r for r in rows if r["outcome"] != "skipped-family-gated"]
    planned = [r for r in considered if r["outcome"] in REAL_ROW_OUTCOMES]
    gaps = [r for r in considered if r["outcome"] not in REAL_ROW_OUTCOMES]
    if not gaps:
        return "ready", len(considered), len(planned), ""
    outs = {r["outcome"] for r in gaps}
    worst = next((o for o in GAP_PRIORITY if o in outs), sorted(outs)[0])
    ex = next(r for r in gaps if r["outcome"] == worst)
    suffix = f" (+{len(gaps) - 1} more)" if len(gaps) > 1 else ""
    return "blocked", len(considered), len(planned), f"{worst}: {ex['table']}/{ex['variant'] or '(default)'}{suffix}"


def coverage_report_rows(output_dir, config, annotations: "Dict[Tuple[str, str], dict]", manifest_dir) -> "List[dict]":
    stats = {}
    for p in sorted(glob.glob(os.path.join(str(output_dir), "cluster_stats", "*_cluster_stats.tsv"))):
        for r in read_rows_tsv(p):
            stats[(r["dataset"], r["cluster"])] = r
    coverage = defaultdict(list)
    for p in sorted(glob.glob(os.path.join(str(manifest_dir), "*", "manifest_coverage.tsv"))):
        for r in read_rows_tsv(p):
            coverage[(r["dataset"], r["cluster"])].append(r)
    clusters = config.get("clusters") or {}
    out = []
    for (d, c), s in sorted(stats.items()):
        cfg = (clusters.get(d) or {}).get(c) or {}
        qp = s["reason"] == "pass"
        has_ann = annotation_lookup_key(d, c, cfg) in annotations
        mex = bool(cfg.get("igvf_manifest_excluded", False))
        has_pred = os.path.exists(os.path.join(str(output_dir), "uniformly_processed", "candidate_e2g_pairs",
                                               f"{d}_{c}_candidate_e2g_pairs.tsv.gz"))
        ref = bool(glob.glob(os.path.join(str(output_dir), "uniformly_processed", d, c, f"{d}_{c}_scE2G_*.e2g.tsv.gz")))
        ms, me, mp, mg = summarize_manifest(coverage.get((d, c), []), mex, qp, has_ann)
        out.append({"dataset": d, "cluster": c, "cell_count": s["cell_count"], "fragments_total": s["fragments_total"],
                    "umi_count": s["umi_count"], "quality_pass": "y" if qp else "n", "exclusion_reason": s["reason"],
                    "has_cell_annotation": "y" if has_ann else "n", "igvf_manifest_excluded": "y" if mex else "n",
                    "predictions_generated": "y" if has_pred else "n", "reformatted": "y" if ref else "n",
                    "status": derive_status(s["reason"], qp, has_ann, has_pred, ref), "manifest_status": ms,
                    "manifest_rows_expected": me, "manifest_rows_planned": mp, "manifest_gap_reason": mg})
    return out


# ---------------------------------------------------------------------------
# Synapse manifests (manage_synapse_manifest.py, list_synapse_orphans.py, list_synapse_prediction_orphans.py)
# ---------------------------------------------------------------------------

SYNAPSE_PARENTS = {"filtered_data": "syn53469844", "predictions": "syn53469845", "candidates": "syn75065641",
                   "features": "syn75936892"}
SYNAPSE_NESTED = {"filtered_data": True, "predictions": True, "candidates": False, "features": True}
SYNAPSE_MODEL_TOKENS = {"predictions": ["multiome_powerlaw_v3", "scATAC_powerlaw_v3"],
                        "features": ["multiome_powerlaw_v3", "scATAC_powerlaw_v3"]}
PREDICTION_ORPHANS_EXCLUDE_TOP = ["(Archived)Y2_versions", "GM12878_10XMultiome", "K562_ERR9847049_Multiome", "igvf7", "kasowski"]


def find_cluster_key(path, cluster_keys, legacy_first_match: bool = False):
    norm, base = str(path).replace(os.sep, "/"), os.path.basename(str(path))
    cands = [(d, c) for d, c in sorted(cluster_keys) if f"/{d}/{c}/" in norm or base.startswith(f"{d}_{c}_")]
    if not cands:
        return None
    return cands[0] if legacy_first_match else max(cands, key=lambda dc: len(dc[1]))


def synapse_key(local_path, dataset, cluster, nested) -> "Tuple[str, str, str]":
    b = os.path.basename(str(local_path))
    if nested:
        return dataset, cluster, b
    return dataset, cluster, b if b.startswith(f"{dataset}_{cluster}_") else f"{dataset}_{cluster}_{b}"


def relative_folder_parts(local_path, dataset, cluster) -> "List[str]":
    marker, norm = f"{dataset}/{cluster}/", str(local_path).replace("\\", "/")
    idx = norm.find(marker)
    if idx == -1:
        raise ValueError(f"{local_path} does not contain {marker}")
    rel = os.path.dirname(norm[idx + len(marker):])
    return [dataset, cluster] + [p for p in rel.split("/") if p]


def synapse_plan(product, should_exist_files, owned: "Dict[Tuple[str, str, str], dict]", cluster_keys, preserved_rows=(),
                 parent_id=None, legacy_first_match=False) -> dict:
    nested = SYNAPSE_NESTED.get(product, True)
    tokens = SYNAPSE_MODEL_TOKENS.get(product, [])
    owned = {k: v for k, v in owned.items() if (k[0], k[1]) in cluster_keys and (not tokens or any(t in k[2] for t in tokens))}
    should, unmatched = {}, []
    for lp in should_exist_files:
        m = find_cluster_key(lp, cluster_keys, legacy_first_match)
        if m is None:
            unmatched.append(lp)
            continue
        should[synapse_key(lp, m[0], m[1], nested)] = lp
    to_delete = {k: v for k, v in owned.items() if k not in should}
    new = {k: p for k, p in should.items() if k not in owned}
    over = {k: p for k, p in should.items() if k in owned}
    parent = parent_id or SYNAPSE_PARENTS.get(product, "")
    rows = list(preserved_rows)
    for k, lp in {**new, **over}.items():
        rows.append({"path": lp, "parent": parent + ("/" + "/".join(relative_folder_parts(lp, k[0], k[1])) if nested else "")})
    return {"to_delete": to_delete, "to_upload_new": new, "to_overwrite": over, "manifest_rows": rows, "unmatched": unmatched,
            "abort_overwrite": product == "filtered_data" and bool(over)}


def load_preserved_rows(manifest_path, cluster_keys) -> "List[dict]":
    if not manifest_path or not os.path.exists(manifest_path):
        return []
    return [r for r in read_rows_tsv(manifest_path) if find_cluster_key(r["path"], cluster_keys) is None]


def read_inventory(path) -> "Dict[Tuple[str, str, str], dict]":
    return {(r["dataset"], r["cluster"], r["name"]): r for r in read_rows_tsv(path)}


# ---------------------------------------------------------------------------
# CATlas distance-to-TSS vs depth figures (CATlas-predictions index.html)
# ---------------------------------------------------------------------------

DEPTH_THRESHOLD = 2e6


def distance_depth_table(stats_rows, prediction_files: "Dict[str, str]", group_col: "Optional[str]" = None,
                         neuron_regex: str = r"neuron|MSN|interneuron|excitatory|inhibitory") -> "List[dict]":
    np = _np()
    out = []
    for r in stats_rows:
        cl = r["cluster"]
        row = {"cluster": cl, "model_name": r.get("model_name"), "fragments_total": float(r["fragments_total"] or 0),
               "cell_count": float(r.get("cell_count") or 0), "mean_dist_to_tss": float(r.get("mean_dist_to_tss") or "nan"),
               "group": r.get(group_col) if group_col else r.get("group", "")}
        row["frag_per_cell"] = row["fragments_total"] / row["cell_count"] if row["cell_count"] else float("nan")
        row["median_dist_to_tss"] = row["q25_dist_to_tss"] = row["q75_dist_to_tss"] = row["iqr_dist_to_tss"] = float("nan")
        pf = prediction_files.get(cl)
        if pf and os.path.exists(pf):
            pred = read_prediction_table(pf)
            if "class" in pred.columns:
                pred = pred[pred["class"] != "promoter"]
            col = "distanceToTSS" if "distanceToTSS" in pred.columns else ("distance" if "distance" in pred.columns else None)
            if col and len(pred):
                v = np.asarray(pred[col], dtype=float)
                q25, med, q75 = np.nanpercentile(v, [25, 50, 75])
                row.update({"median_dist_to_tss": float(med), "q25_dist_to_tss": float(q25), "q75_dist_to_tss": float(q75),
                            "iqr_dist_to_tss": float(q75 - q25)})
        label = str(row["group"] or "") + " " + cl
        row["neuron"] = bool(re.search(neuron_regex, label, re.IGNORECASE))
        row["above_threshold"] = row["fragments_total"] >= DEPTH_THRESHOLD
        out.append(row)
    return out


def log_trend(rows, ycol) -> "Optional[Tuple[float, float, int]]":
    np = _np()
    pts = [(math.log10(r["fragments_total"]), r[ycol]) for r in rows if r["above_threshold"] and r["fragments_total"] > 0
           and math.isfinite(r[ycol])]
    if len(pts) < 2:
        return None
    x, y = np.array(pts).T
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept), len(pts)


# ---------------------------------------------------------------------------
# IGVF metadata tables (igvf_metadata/registry.py, context.py, refs.py, tables/*.py)
# ---------------------------------------------------------------------------

IGVF_DEFAULTS = {"lab": "/labs/jesse-engreitz/", "award": "/awards/HG011972/", "alias_prefix": "jesse-engreitz",
                 "enabled_families": ("Multiome",)}
FAMILY_DISPLAY = {"multiome_powerlaw_v3": "Multiome", "multiome_megamap_v3": "Multiome",
                  "scATAC_powerlaw_v3": "scATAC", "scATAC_megamap_v3": "scATAC"}
KUNDAJE_PREFIX = "anshul-kundaje"
PRINCIPAL_PSEUDOBULK_ASV = "/analysis-step-versions/28d9aece-0621-458c-949f-0314c8301228/"
PLOTTING_ASV = "/analysis-step-versions/0077c8e1-f3f7-4e4e-b79e-6e6560820c9b/"
RUN_SCE2G_ASV = "jesse-engreitz:analysis_step_v1_run_scE2G"
REF_ELEMENT_TO_GENE = ["IGVFFI7969JLFC", "IGVFFI0653VCGH", "IGVFFI9573KOZR"]
REF_GENE_QUANT = ["IGVFFI9573KOZR"]
REF_ATAC_FRAGMENTS = ["IGVFFI0653VCGH", "IGVFFI6788CPPS"]
ATAC_FRAGMENT_SPEC_DOC = "/documents/db2a6dd0-cc1d-439e-a610-f9f1d04cfd82/"
RNA_MATRIX_SPEC = "jesse-engreitz:rna_matrix_market_tar_archive_file_format"
BARCODE_LIST_SPEC = "jesse-engreitz:filtered_barcode_membership_file_format"
ARRAY_FIELDS = {"aliases", "derived_from", "reference_files", "file_format_specifications", "documents", "samples",
                "input_file_sets"}
COLLECTIONS = {"tabular_file": "tabular-files", "matrix_file": "matrix-files", "index_file": "index-files",
               "signal_file": "signal-files", "document": "documents", "prediction_set": "prediction-sets",
               "pseudobulk_set": "pseudobulk-sets"}
FILE_OBJECT_TYPES = {"tabular_file", "matrix_file", "index_file", "signal_file"}
MANIFEST_COVERAGE_HEADER = ["dataset", "cluster", "model", "table", "variant", "outcome", "reason"]


def family(model: str) -> str:
    try:
        return FAMILY_DISPLAY[model]
    except KeyError:
        raise ValueError(f"no display-family mapping for model {model!r} -- add it to FAMILY_DISPLAY")


def igvf_config(config: dict) -> dict:
    d = dict(IGVF_DEFAULTS)
    for k, v in ((config or {}).get("igvf") or {}).items():
        if k in d and v is not None:
            d[k] = tuple(v) if k == "enabled_families" else v
    return d


class Ctx:
    def __init__(self, dataset, cluster, model, cfg, igvf, scE2G_dir, data_dir, output_dir, state=None,
                 threshold_overrides=None, results_dir_base=None):
        self.dataset, self.cluster, self.model, self.cfg, self.igvf = dataset, cluster, model, cfg, igvf
        self.scE2G_dir, self.data_dir, self.output_dir, self.state = scE2G_dir, data_dir, output_dir, state
        self.thr = threshold_overrides or {}
        self.results_base = results_dir_base

    def with_model(self, model):
        return Ctx(self.dataset, self.cluster, model, self.cfg, self.igvf, self.scE2G_dir, self.data_dir, self.output_dir,
                   self.state, self.thr, self.results_base)

    @property
    def results_dir(self):
        if self.results_base:
            return str(self.results_base)
        if self.output_dir:
            return os.path.join(self.output_dir, "uniformly_processed")
        return os.path.join(str(self.scE2G_dir), "results", "uniformly_processed")

    @property
    def cluster_dir(self):
        return os.path.join(self.results_dir, self.dataset, self.cluster)

    @property
    def multiome_data_cluster_dir(self):
        base = self.output_dir if self.output_dir else self.data_dir
        return os.path.join(str(base), "multiome_data", self.dataset, self.cluster)

    def alias(self, *parts) -> str:
        return f"{self.igvf['alias_prefix']}:" + "_".join(str(p) for p in parts)

    def threshold(self) -> str:
        return model_threshold(self.scE2G_dir, self.model, self.thr)

    def metadata(self) -> dict:
        row = self.state.get_cell_annotation(self.dataset, self.cluster) if self.state else None
        if row is None:
            raise ValueError(f"{self.dataset}/{self.cluster}: no cached Cell Annotation metadata -- run cell-metadata first")
        return row

    def subsamples_unique(self):
        return unique_subsamples(self.cfg["qc_guide"])

    def subsamples_freq(self):
        return subsamples_by_frequency(self.cfg["qc_guide"])


class Variant:
    def __init__(self, name, row, enabled=None, paths=None, deps=None):
        self.name, self.row = name, row
        self.enabled = enabled or (lambda ctx: True)
        self.paths = paths or (lambda ctx: [])
        self.deps = deps or (lambda ctx: [])


class TableSpec:
    def __init__(self, name, object_type, scope, alias, required, variants, constant=None, scope_fields=None):
        self.name, self.object_type, self.scope, self.alias = name, object_type, scope, alias
        self.required, self.variants = required, variants
        self.constant = constant or {}
        self.scope_fields = scope_fields or (lambda ctx: {})


# path builders (prediction_tabular_files.py, filtered_*.py)
def p_full(c):
    return os.path.join(c.cluster_dir, f"{c.dataset}_{c.cluster}_scE2G_{c.model}.e2g.tsv.gz")


def p_thresholded(c):
    return os.path.join(c.cluster_dir, f"{c.dataset}_{c.cluster}_scE2G_{c.model}_threshold{c.threshold()}.e2g.tsv.gz")


def p_bedpe(c):
    return os.path.join(c.cluster_dir, f"{c.dataset}_{c.cluster}_scE2G_{c.model}_threshold{c.threshold()}.bedpe.gz")


def p_elements_bed(c):
    return os.path.join(c.cluster_dir, f"{c.dataset}_{c.cluster}_element_list.bed.gz")


def p_genes(c):
    return os.path.join(c.cluster_dir, f"{c.dataset}_{c.cluster}_scE2G_multiome_v3_gene_list.tsv.gz")


def p_atac(c):
    return os.path.join(c.multiome_data_cluster_dir, f"atac_fragments_{c.dataset}_{c.cluster}.tsv.gz")


def p_rna_tar(c):
    return os.path.join(c.multiome_data_cluster_dir, f"rna_count_matrix_{c.dataset}_{c.cluster}.tar.gz")


def p_qc_thresholds(c):
    return os.path.join(os.path.dirname(c.cfg["qc_guide"]), "qc_thresholds.tsv")


def build_tables() -> "Dict[str, TableSpec]":
    """The eleven upstream tables, registered in upstream's tables/__init__.py order."""
    T = {}

    def a_qc(c, v=""):
        return c.alias(c.dataset, c.cluster, "QC_thresholds")

    def a_pps(c, v=""):
        return c.alias(c.dataset, c.cluster, "filtered_pseudobulk_set")

    def a_fbl(c, v=""):
        return c.alias(c.dataset, c.cluster, "filtered_barcode_list")

    def a_atac(c, v=""):
        return c.alias(c.dataset, c.cluster, "filtered_ATAC_fragment_file")

    def a_rna(c, v=""):
        return c.alias(c.dataset, c.cluster, "filtered_RNA_count_matrix")

    def a_ps(c, v=""):
        return c.alias(c.dataset, c.cluster, "scE2G", family(c.model), "predictions")

    def a_ptf(c, v):
        return c.alias(c.dataset, c.cluster, "scE2G", family(c.model), "predictions", v)

    def fam(c):
        return family(c.model)

    e2g_specs = lambda c: [c.alias("element_to_gene_interaction_predictions_file_format_pdf"),  # noqa: E731
                           c.alias("element_to_gene_interaction_predictions_file_format_md")]

    def full_row(c):
        parts = [c.alias("scE2G", fam(c), "trained_model"), a_atac(c), a_ptf(c, "elements_bed")]
        if fam(c) == "Multiome":
            parts += [a_rna(c), a_ptf(c, "genes")]
        return {"content_type": "element to gene interactions", "file_format": "tsv",
                "description": f"Full scE2G ({fam(c)}) predictions for {c.dataset} {c.cluster} cells",
                "derived_from": ",".join(parts), "file_format_specifications": e2g_specs(c),
                "reference_files": list(REF_ELEMENT_TO_GENE), "submitted_file_name": p_full(c)}

    def thr_row(c):
        return {"content_type": "element to gene interactions", "file_format": "tsv",
                "description": f"Thresholded scE2G ({fam(c)}) predictions for {c.dataset} {c.cluster} cells",
                "derived_from": a_ptf(c, "full"), "file_format_specifications": e2g_specs(c),
                "reference_files": list(REF_ELEMENT_TO_GENE), "submitted_file_name": p_thresholded(c)}

    def bedpe_row(c):
        return {"content_type": "element to gene interactions", "file_format": "bedpe",
                "description": (f"Bedpe file for genome browser visualization of thresholded scE2G ({fam(c)}) predictions "
                                f"for {c.dataset} {c.cluster} cells"),
                "derived_from": a_ptf(c, "thresholded"), "file_format_specifications": [c.alias("E2G_bedpe_file_format")],
                "reference_files": list(REF_ELEMENT_TO_GENE), "submitted_file_name": p_bedpe(c)}

    def eb_row(c):
        return {"content_type": "elements reference", "file_format": "bed", "file_format_type": "bed3+",
                "description": (f"Annotated elements in scE2G ({fam(c)}) predictions for {c.dataset} {c.cluster} cells "
                                "for IGV visualization"),
                "derived_from": a_atac(c),
                "file_format_specifications": [c.alias("elements_reference_bed_file_format_specification_pdf"),
                                               c.alias("elements_reference_bed_file_format_specification_md")],
                "reference_files": list(REF_ELEMENT_TO_GENE), "submitted_file_name": p_elements_bed(c)}

    def genes_row(c):
        return {"content_type": "gene quantifications", "file_format": "tsv",
                "description": f"Annotated genes in scE2G ({fam(c)}) predictions for {c.dataset} {c.cluster} cells",
                "derived_from": a_rna(c),
                "file_format_specifications": [c.alias("gene_quantifications_file_format_specification_pdf"),
                                               c.alias("gene_quantifications_file_format_specification_md")],
                "reference_files": list(REF_GENE_QUANT), "submitted_file_name": p_genes(c)}

    def full_deps(c):
        d = [("prediction_set", ""), ("filtered_atac_fragment_file", ""), ("prediction_tabular_files", "elements_bed")]
        if fam(c) == "Multiome":
            d += [("prediction_tabular_files", "genes"), ("filtered_rna_count_matrix", "")]
        return d

    ps_dep = [("prediction_set", "")]
    ptf_scope = lambda c: {"md5sum": None, "file_set": a_ps(c)}  # noqa: E731
    T["prediction_tabular_files"] = TableSpec(
        "prediction_tabular_files", "tabular_file", "cluster_model", a_ptf,
        ["aliases", "award", "lab", "content_type", "controlled_access", "file_format", "file_set"],
        [Variant("full", full_row, paths=lambda c: [p_full(c)], deps=full_deps),
         Variant("thresholded", thr_row, paths=lambda c: [p_thresholded(c)],
                 deps=lambda c: ps_dep + [("prediction_tabular_files", "full")]),
         Variant("bedpe", bedpe_row, paths=lambda c: [p_bedpe(c)],
                 deps=lambda c: ps_dep + [("prediction_tabular_files", "thresholded")]),
         Variant("elements_bed", eb_row, paths=lambda c: [p_elements_bed(c)],
                 deps=lambda c: ps_dep + [("filtered_atac_fragment_file", "")]),
         Variant("genes", genes_row, enabled=lambda c: fam(c) == "Multiome", paths=lambda c: [p_genes(c)],
                 deps=lambda c: ps_dep + [("filtered_rna_count_matrix", "")])],
        constant={"controlled_access": False, "filtered": True, "derived_manually": False,
                  "analysis_step_version": RUN_SCE2G_ASV},
        scope_fields=ptf_scope)

    T["signal_files"] = TableSpec(
        "signal_files", "signal_file", "cluster_model", lambda c, v: c.alias(c.dataset, c.cluster, "scE2G", fam(c), "ATAC_bw"),
        ["aliases", "award", "lab", "file_format", "file_set", "content_type"],
        [Variant("atac_bw", lambda c: {"file_format": "bigWig", "content_type": "read-depth signal",
                                       "reference_files": ["IGVFFI7969JLFC", "IGVFFI0653VCGH"],
                                       "strand_specificity": "unstranded", "derived_from": a_atac(c),
                                       "submitted_file_name": os.path.join(c.cluster_dir, "ATAC_norm.bw")},
                 paths=lambda c: [os.path.join(c.cluster_dir, "ATAC_norm.bw")],
                 deps=lambda c: ps_dep + [("filtered_atac_fragment_file", "")])],
        constant={"normalized": True, "analysis_step_version": RUN_SCE2G_ASV}, scope_fields=ptf_scope)

    idx_req = ["aliases", "award", "lab", "file_format", "file_set", "content_type", "controlled_access"]
    idx_const = {"controlled_access": False, "derived_manually": False, "analysis_step_version": RUN_SCE2G_ASV}
    T["bedpe_index_file"] = TableSpec(
        "bedpe_index_file", "index_file", "cluster_model",
        lambda c, v: c.alias(c.dataset, c.cluster, "scE2G", fam(c), "predictions", "bedpe", "index"), idx_req,
        [Variant("bedpe_index", lambda c: {"file_format": "tbi", "content_type": "index", "derived_from": a_ptf(c, "bedpe"),
                                           "description": f"index file for thresholded scE2G {fam(c)} predictions BEDPE file",
                                           "submitted_file_name": p_bedpe(c) + ".tbi"},
                 paths=lambda c: [p_bedpe(c) + ".tbi"], deps=lambda c: ps_dep + [("prediction_tabular_files", "bedpe")])],
        constant=dict(idx_const), scope_fields=ptf_scope)
    T["elements_bed_index_file"] = TableSpec(
        "elements_bed_index_file", "index_file", "cluster_model",
        lambda c, v: c.alias(c.dataset, c.cluster, "scE2G", fam(c), "predictions", "elements_bed", "index"), idx_req,
        [Variant("elements_bed_index", lambda c: {
            "file_format": "tbi", "content_type": "index", "derived_from": a_ptf(c, "elements_bed"),
            "description": (f"index file for annotated candidate elements in scE2G ({fam(c)}) predictions for "
                            f"{c.dataset} {c.cluster} cells"), "submitted_file_name": p_elements_bed(c) + ".tbi"},
            paths=lambda c: [p_elements_bed(c) + ".tbi"], deps=lambda c: ps_dep + [("prediction_tabular_files", "elements_bed")])],
        constant=dict(idx_const), scope_fields=ptf_scope)

    def ps_row(c):
        md = c.metadata()
        return {"input_file_sets": [c.alias("scE2G", fam(c), "model"), a_pps(c)],
                "documents": [c.alias("E2G_prediction_set_file_format_specs_pdf"), c.alias("E2G_prediction_set_file_format_specs_md")],
                "description": f"scE2G {fam(c)} predictions for {c.dataset} {c.cluster} cells",
                "samples": c.subsamples_freq(), "cell_type": f"/sample-terms/{md['cl_id']}/", "cell_qualifier": md["cell_qualifier"]}
    T["prediction_set"] = TableSpec(
        "prediction_set", "prediction_set", "cluster_model", a_ps, ["aliases", "award", "lab"],
        [Variant("", ps_row, deps=lambda c: [("principal_pseudobulk_set", "")])],
        constant={"file_set_type": "element-gene links", "scope": "genome-wide", "submitter_comment": "Version 1"})

    def pps_row(c):
        md = c.metadata()
        return {"cell_type": f"/sample-terms/{md['cl_id']}/", "cell_qualifier": md["cell_qualifier"], "documents": [a_qc(c)],
                "samples": c.subsamples_freq(),
                "input_file_sets": [f"{KUNDAJE_PREFIX}:{c.dataset}-{c.cluster}-{s}" for s in c.subsamples_unique()],
                "description": (f"Filtered datafiles describing a single annotated cell cluster ({c.cluster}); these datafiles "
                                "are inputs to E2G predictive models")}
    T["principal_pseudobulk_set"] = TableSpec(
        "principal_pseudobulk_set", "pseudobulk_set", "cluster", a_pps, ["aliases", "award", "lab", "file_set_type"],
        [Variant("", pps_row, deps=lambda c: [("QC_documents", "")])],
        constant={"file_set_type": "pseudobulk analysis", "merged": True,
                  "submitter_comment": ("This principal pseudobulk set is based on status: 'in progress' data and will be "
                                        "superseded if those data are updated")})

    T["QC_documents"] = TableSpec(
        "QC_documents", "document", "cluster", a_qc, ["aliases", "award", "lab", "document_type", "attachment"],
        [Variant("", lambda c: {"description": f"Quality Control thresholds applied to each cell in {c.cluster}",
                                "document_type": "pipeline parameters", "attachment": {"path": p_qc_thresholds(c)}},
                 paths=lambda c: [p_qc_thresholds(c)])])

    pps_scope = lambda c: {"md5sum": None, "file_set": a_pps(c)}  # noqa: E731
    file_req = ["aliases", "award", "lab", "file_format", "file_set", "content_type", "controlled_access"]
    T["filtered_barcode_list"] = TableSpec(
        "filtered_barcode_list", "tabular_file", "cluster", a_fbl, file_req,
        [Variant("", lambda c: {
            "file_format": "tsv", "content_type": "pseudobulk barcode list", "documents": [a_qc(c)],
            "file_format_specifications": [BARCODE_LIST_SPEC],
            "derived_from": ",".join(f"{KUNDAJE_PREFIX}:{c.dataset}-{c.cluster}-{s}-per_cell_qc_tsv" for s in c.subsamples_unique()),
            "description": (f"Filtered list of 16-bp cell barcodes defining membership in the {c.cluster}; see QC thresholds "
                            "for a record of the filters applied"),
            "submitter_comment": c.cfg.get("submitter_comment"), "analysis_step_version": PLOTTING_ASV,
            "submitted_file_name": c.cfg["qc_guide"]},
            deps=lambda c: [("principal_pseudobulk_set", ""), ("QC_documents", "")])],
        constant={"controlled_access": False, "filtered": True, "derived_manually": False}, scope_fields=pps_scope)

    T["atac_index_file"] = TableSpec(
        "atac_index_file", "index_file", "cluster", lambda c, v: c.alias(c.dataset, c.cluster, "filtered_ATAC_fragment_file_index"),
        file_req,
        [Variant("", lambda c: {"file_format": "tbi", "content_type": "index", "analysis_step_version": PRINCIPAL_PSEUDOBULK_ASV,
                                "derived_from": a_atac(c),
                                "description": (f"ATAC fragment file index used as input for E2G predictions for {c.dataset} "
                                                f"{c.cluster} cells"), "submitted_file_name": p_atac(c) + ".tbi"},
                 paths=lambda c: [p_atac(c) + ".tbi"],
                 deps=lambda c: [("principal_pseudobulk_set", ""), ("filtered_atac_fragment_file", "")])],
        constant={"controlled_access": False, "derived_manually": False}, scope_fields=pps_scope)

    T["filtered_atac_fragment_file"] = TableSpec(
        "filtered_atac_fragment_file", "tabular_file", "cluster", a_atac, file_req,
        [Variant("", lambda c: {
            "file_format": "bed", "file_format_type": "bed3+", "content_type": "fragments",
            "file_format_specifications": ATAC_FRAGMENT_SPEC_DOC, "reference_files": list(REF_ATAC_FRAGMENTS),
            "analysis_step_version": PRINCIPAL_PSEUDOBULK_ASV,
            "derived_from": ",".join([f"{KUNDAJE_PREFIX}:{c.dataset}-{c.cluster}-{s}-fragments_tsv_gz" for s in c.subsamples_unique()]
                                     + [a_fbl(c)]),
            "description": f"Filtered ATAC fragment file containing reads from cells annotated as {c.cluster}",
            "submitted_file_name": p_atac(c)},
            paths=lambda c: [p_atac(c)], deps=lambda c: [("principal_pseudobulk_set", ""), ("filtered_barcode_list", "")])],
        constant={"filtered": True, "controlled_access": False, "derived_manually": False}, scope_fields=pps_scope)

    T["filtered_rna_count_matrix"] = TableSpec(
        "filtered_rna_count_matrix", "matrix_file", "cluster", a_rna, ["aliases", "award", "lab", "file_format", "file_set", "content_type"],
        [Variant("", lambda c: {
            "file_format": "tar", "content_type": "cell by gene matrix", "reference_files": ["/reference-files/IGVFFI9573KOZR/"],
            "analysis_step_version": PRINCIPAL_PSEUDOBULK_ASV,
            "derived_from": ",".join([f"{KUNDAJE_PREFIX}:{c.dataset}-{c.cluster}-{s}-rna_counts_mtx_h5ad" for s in c.subsamples_unique()]
                                     + [a_fbl(c)]),
            "description": f"Filtered Matrix Market file containing RNA transcripts for cells annotated as {c.cluster}",
            "file_format_specifications": [RNA_MATRIX_SPEC], "submitted_file_name": p_rna_tar(c)},
            paths=lambda c: [p_rna_tar(c)], deps=lambda c: [("principal_pseudobulk_set", ""), ("filtered_barcode_list", "")])],
        constant={"filtered": True, "derived_manually": False}, scope_fields=pps_scope)
    return T


TABLES = build_tables()


def build_payload(table: TableSpec, variant: Variant, ctx: Ctx, alias: str) -> dict:
    p = {"aliases": [alias], "award": ctx.igvf["award"], "lab": ctx.igvf["lab"]}
    p.update(table.constant)
    p.update(table.scope_fields(ctx))
    p.update(variant.row(ctx))
    return p


def is_missing(v) -> bool:
    return v is None or (isinstance(v, (str, list, tuple, dict)) and len(v) == 0)


def validate_row(table: TableSpec, variant_name: str, payload: dict) -> None:
    missing = [c for c in table.required if is_missing(payload.get(c))]
    if missing:
        raise ValueError(f"{table.name}/{variant_name}: missing required column(s) {missing}")


def tsv_cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return json.dumps(list(v)) if v and isinstance(v[0], dict) else ",".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v)
    if isinstance(v, bool):
        return str(v).lower()
    return str(v)


def write_iu_tsv(path, rows, record_ids=None) -> "Optional[Path]":
    """iu_register.py create_payloads_from_tsv format: plain tab join, never csv-quoted."""
    if not rows:
        return None
    cols = sorted({k for r in rows for k in r})
    if record_ids is not None:
        cols = ["record_id"] + cols
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows):
            vals = dict(r)
            if record_ids is not None:
                vals["record_id"] = record_ids[i]
            fh.write("\t".join(tsv_cell(vals.get(c)) for c in cols) + "\n")
    return path


def read_iu_tsv(path) -> "Dict[str, dict]":
    if not os.path.exists(path):
        return {}
    lines = Path(path).read_text().splitlines()
    if not lines:
        return {}
    cols = lines[0].split("\t")
    out = {}
    for ln in lines[1:]:
        cells = dict(zip(cols, ln.split("\t")))
        if cells.get("record_id"):
            out[cells["record_id"]] = cells
    return out


def merge_write_iu_tsv(path, rows, record_ids) -> "Optional[Path]":
    existing = read_iu_tsv(path)
    for r, rid in zip(rows, record_ids):
        existing[str(rid)] = {"record_id": str(rid), **{k: tsv_cell(v) for k, v in r.items()}}
    if not existing:
        return None
    cols = ["record_id"] + sorted({k for r in existing.values() for k in r if k != "record_id"})
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in existing.values():
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")
    return Path(path)


def iter_scopes(table: TableSpec, cluster_keys, cluster_configs, igvf):
    for d, c in sorted(cluster_keys):
        cfg = cluster_configs[(d, c)]
        models = [m for m in cfg.get("models") or [] if family(m) in igvf["enabled_families"]] if table.scope == "cluster_model" else [None]
        for m in models:
            yield d, c, m, cfg


def dep_model_key(dep_table: str, model_key: str) -> str:
    return model_key if TABLES[dep_table].scope == "cluster_model" else ""


def compute_round(cache: dict, make_ctx, d, c, model, table_name, variant_name) -> int:
    key = (d, c, model or "", table_name, variant_name)
    if key in cache:
        return cache[key]
    t = TABLES[table_name]
    v = next(x for x in t.variants if x.name == variant_name)
    deps = v.deps(make_ctx(d, c, model))
    r = 1 if not deps else 1 + max(compute_round(cache, make_ctx, d, c, model if TABLES[dt].scope == "cluster_model" else None, dt, dv)
                                   for dt, dv in deps)
    cache[key] = r
    return r


class OfflineReader:
    """Alias lookups answered from a {alias: record} map (default: nothing is live)."""

    def __init__(self, live: "Optional[dict]" = None):
        self.live = live or {}

    def get_by_alias(self, alias, database: bool = False):
        return self.live.get(alias)


class PortalAliasReader:
    """Read-only alias existence check against the submission target (igvf_submission_skill.Portal)."""

    def __init__(self, portal):
        self.portal = portal

    def get_by_alias(self, alias, database: bool = False):
        from igvf_submission_skill import PortalError  # type: ignore
        params = {"frame": "object"}
        if database:
            params["datastore"] = "database"
        try:
            rec = self.portal._request("GET", f"/{alias}/", params=params)
        except PortalError as e:
            if "HTTP 404" in str(e):
                return None
            raise
        return rec or None


def plan_manifest(cluster_keys, cluster_configs, make_ctx, state: State, reader, igvf: dict,
                  table_names: "Optional[Sequence[str]]" = None) -> dict:
    """plan_table over every table: coverage, post/patch entries with rounds, accumulator rows."""
    round_cache, coverage, posts, patches, records = {}, [], [], [], []
    for t in TABLES.values():
        if table_names and t.name not in table_names:
            continue
        for d, c, m, cfg in iter_scopes(t, cluster_keys, cluster_configs, igvf):
            ctx = make_ctx(d, c, m)
            mk = m or ""
            for v in t.variants:
                def cov(outcome, reason=""):
                    coverage.append({"dataset": d, "cluster": c, "model": mk, "table": t.name, "variant": v.name,
                                     "outcome": outcome, "reason": reason})
                try:
                    en = v.enabled(ctx)
                except Exception as e:  # noqa: BLE001
                    cov("enabled-check-failed", str(e))
                    continue
                if not en:
                    cov("skipped-family-gated")
                    continue
                try:
                    missing = [p for p in v.paths(ctx) if not os.path.exists(p)]
                except Exception as e:  # noqa: BLE001
                    cov("enabled-check-failed", str(e))
                    continue
                if missing:
                    cov("skipped-missing-file", ";".join(missing))
                    continue
                alias = t.alias(ctx, v.name)
                deferred = None
                for dt, dv in v.deps(ctx):
                    up = state.get_upload(d, c, dep_model_key(dt, mk), dt, dv)
                    if not up or up["status"] != "uploaded":
                        deferred = f"{dt}/{dv or '(default)'}"
                        break
                rnd = compute_round(round_cache, make_ctx, d, c, m, t.name, v.name)
                try:
                    payload = build_payload(t, v, ctx, alias)
                    validate_row(t, v.name, payload)
                except Exception as e:  # noqa: BLE001
                    cov("invalid", str(e))
                    continue
                h = payload_hash(payload)
                ex = state.get_upload(d, c, mk, t.name, v.name)
                if ex and ex["status"] == "uploaded" and ex["payload_hash"] == h:
                    cov("unchanged", alias)
                    records.append({"dataset": d, "object_type": t.object_type, "record_id": ex["portal_id"], "payload": payload})
                    continue
                row = state.claim_pending(d, c, mk, t.name, v.name, alias, h)
                live = reader.get_by_alias(alias)
                entry = {"row_id": row["id"], "alias": alias, "payload": payload, "table": t.name, "variant": v.name,
                         "object_type": t.object_type, "round": rnd, "dataset": d, "cluster": c, "model": mk,
                         "deferred_on": deferred, "deps": [(dt, dv, dep_model_key(dt, mk)) for dt, dv in v.deps(ctx)]}
                if live is None:
                    posts.append(entry)
                    cov("planned-post", alias)
                else:
                    entry["record_id"] = live.get("uuid") or live.get("accession") or live.get("@id")
                    patches.append(entry)
                    cov("planned-patch", alias)
    return {"coverage": coverage, "posts": posts, "patches": patches, "records": records}


def write_manifest_files(plan: dict, manifest_dir: Path) -> "List[dict]":
    """round<N>_<table>[_<variant>]_{post,patch}.tsv per dataset; PATCH retires a stale sibling POST."""
    groups, this_run = defaultdict(list), set()
    for kind, entries in (("post", plan["posts"]), ("patch", plan["patches"])):
        for e in entries:
            groups[(e["dataset"], e["table"], e["variant"], e["round"], kind)].append(e)
    written = []
    for (d, t, v, r, kind), entries in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][3], kv[0][1], kv[0][2])):
        suffix = f"_{v}" if v else ""
        path = manifest_dir / d / f"round{r}_{t}{suffix}_{kind}.tsv"
        write_iu_tsv(path, [e["payload"] for e in entries], [e["record_id"] for e in entries] if kind == "patch" else None)
        written.append({"path": path, "object_type": entries[0]["object_type"], "kind": kind, "round": r, "rows": len(entries)})
        this_run.add(path)
        sib = manifest_dir / d / f"round{r}_{t}{suffix}_{'post' if kind == 'patch' else 'patch'}.tsv"
        if sib in this_run:
            continue  # deviation: upstream would delete a POST file written moments ago for another cluster
        if sib.exists() and kind == "patch":
            sib.unlink()
        elif sib.exists():
            print(f"WARNING {d}/{sib.name} exists but this round planned a POST -- left both files in place")
    acc = defaultdict(list)
    for rec in plan["records"]:
        acc[(rec["dataset"], rec["object_type"])].append(rec)
    for (d, ot), recs in acc.items():
        merge_write_iu_tsv(manifest_dir / d / f"{ot}.tsv", [x["payload"] for x in recs], [x["record_id"] for x in recs])
    by_ds = defaultdict(list)
    for r in plan["coverage"]:
        by_ds[r["dataset"]].append(r)
    for d, rows in by_ds.items():
        write_rows_tsv(manifest_dir / d / "manifest_coverage.tsv", MANIFEST_COVERAGE_HEADER,
                       sorted(rows, key=lambda r: (r["cluster"], r["table"], r["variant"], r["model"])), quiet=True)
    return written


def attachment_object(path: str) -> dict:
    """igvf_utils' {"path": ...} expansion: {download, type, href: data URI (base64)}."""
    mime = mimetypes.guess_type(path)[0] or ("text/tab-separated-values" if path.endswith(".tsv") else "text/plain")
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return {"download": os.path.basename(path), "type": mime, "href": f"data:{mime};base64,{data}"}


def rest_payload(payload: dict, object_type: str, compute_md5: bool = True) -> dict:
    """iu_register's TSV -> JSON conversion for a direct REST call: arrays split, blanks dropped, attachment inlined,
    md5sum computed from submitted_file_name for file objects."""
    out = {}
    for k, v in payload.items():
        if is_missing(v) and v is not False:
            continue
        if k in ARRAY_FIELDS and isinstance(v, str):
            v = [x for x in v.split(",") if x]
        if k == "attachment" and isinstance(v, dict) and set(v) == {"path"}:
            v = attachment_object(v["path"])
        out[k] = v
    if compute_md5 and object_type in FILE_OBJECT_TYPES and "md5sum" not in out and out.get("submitted_file_name") \
            and os.path.exists(out["submitted_file_name"]):
        out["md5sum"] = md5_file(out["submitted_file_name"])
    return out


def s3_upload(creds: dict, path: str) -> str:
    """Upload a posted file's bytes with the upload_credentials the Portal returned (boto3), else print the command."""
    url = creds.get("upload_url") or creds.get("url")
    try:
        import boto3  # type: ignore
        bucket, _, key = url.replace("s3://", "", 1).partition("/")
        boto3.client("s3", aws_access_key_id=creds["access_key"], aws_secret_access_key=creds["secret_key"],
                     aws_session_token=creds.get("session_token")).upload_file(path, bucket, key)
        return "uploaded"
    except ImportError:
        print(f"boto3 not installed: upload with AWS_ACCESS_KEY_ID=... AWS_SECRET_ACCESS_KEY=... AWS_SESSION_TOKEN=... "
              f"aws s3 cp {path} {url}")
        return "needs_s3_upload"


def execute_plan(plan: dict, client, state: State, reader, log_rows: "List[dict]") -> dict:
    """POST / PATCH in round order; a row whose dependency is not uploaded in the ledger is deferred.
    `client` needs _request(method, path, payload=...) -- igvf_submission_skill.Portal."""
    entries = sorted([("post", e) for e in plan["posts"]] + [("patch", e) for e in plan["patches"]],
                     key=lambda ke: (ke[1]["round"], ke[1]["dataset"], ke[1]["table"], ke[1]["variant"], ke[1]["cluster"]))
    uploaded, failed, deferred = [], [], []
    for kind, e in entries:
        waiting = [f"{dt}/{dv or '(default)'}" for dt, dv, dm in e["deps"]
                   if (state.get_upload(e["dataset"], e["cluster"], dm, dt, dv) or {}).get("status") != "uploaded"]
        if waiting:
            deferred.append({**{k: e[k] for k in ("alias", "table", "variant", "round")}, "waiting_on": waiting[0]})
            continue
        body = rest_payload(e["payload"], e["object_type"])
        try:
            if kind == "post":
                res = client._request("POST", f"/{COLLECTIONS[e['object_type']]}/", payload=body)
            else:
                res = client._request("PATCH", f"/{str(e['record_id']).strip('/')}/", payload=body)
            item = (res.get("@graph") or [res])[0] if isinstance(res, dict) else {}
            if kind == "post" and e["object_type"] in FILE_OBJECT_TYPES and item.get("upload_credentials"):
                s3_upload(item["upload_credentials"], body["submitted_file_name"])
        except Exception as exc:  # noqa: BLE001
            state.record_result(e["row_id"], "failed", error=str(exc)[:500])
            failed.append({"alias": e["alias"], "kind": kind, "error": str(exc)[:300]})
            log_rows.append({"alias": e["alias"], "kind": kind, "result": "error", "detail": str(exc)[:300]})
            continue
        rec = reader.get_by_alias(e["alias"], database=True)
        if rec:
            pid = rec.get("uuid") or rec.get("accession") or rec.get("@id")
            state.record_result(e["row_id"], "uploaded", portal_id=pid)
            uploaded.append({"alias": e["alias"], "kind": kind, "record_id": pid, "table": e["table"], "variant": e["variant"]})
            log_rows.append({"alias": e["alias"], "kind": kind, "result": "uploaded", "detail": pid})
        else:
            state.record_result(e["row_id"], "failed", error="not found on portal after upload attempt")
            failed.append({"alias": e["alias"], "kind": kind, "error": "not readable back"})
            log_rows.append({"alias": e["alias"], "kind": kind, "result": "unverified", "detail": ""})
    return {"uploaded": uploaded, "failed": failed, "deferred": deferred}


EXTERNAL_REF_FIELDS = ("derived_from", "reference_files", "file_format_specifications", "analysis_step_version",
                       "input_file_sets", "documents", "samples", "cell_type")


def external_references(payloads: "Iterable[dict]", own_aliases: "set") -> "List[str]":
    refs = set()
    for p in payloads:
        for f in EXTERNAL_REF_FIELDS:
            v = p.get(f)
            vals = v if isinstance(v, list) else ([x for x in v.split(",") if x] if isinstance(v, str) else [])
            for x in vals:
                if x and x not in own_aliases:
                    refs.add(x)
    return sorted(refs)


def check_lineage(refs: "Sequence[str]", fetch, walk_from: "Sequence[str]" = ()) -> dict:
    """Resolve every external reference (derived_from, input_file_sets, reference_files, ...) with an authenticated
    GET, and walk the Portal graph (portal_lineage.walk) from each dataset's principal analysis set."""
    resolved = {}
    for r in refs:
        status, obj = fetch(f"/{r.strip('/')}/?format=json&frame=object")
        obj = obj or {}
        resolved[r] = {"http_status": status, "@id": obj.get("@id"), "accession": obj.get("accession"),
                       "status": obj.get("status"), "ok": status == 200 and obj.get("status") not in ("deleted", "revoked"),
                       "blocked": status == 403}
    walks = {}
    if walk_from:
        from portal_lineage import walk  # type: ignore
        for acc in walk_from:
            g = walk(acc, fetch=fetch, max_depth=3, max_nodes=200, fetch_qc=False, reverse_search=False, workers=1)
            types = Counter(n.get("type") for n in g.nodes.values())
            walks[acc] = {"nodes": len(g.nodes), "types": dict(types), "blocked": dict(g.blocked),
                          "pseudobulk_sets": sorted(n.get("accession") for n in g.nodes.values() if n.get("type") == "PseudobulkSet")}
    missing = [r for r, v in resolved.items() if not v["ok"] and not v["blocked"]]
    return {"resolved": resolved, "missing": missing, "blocked": [r for r, v in resolved.items() if v["blocked"]], "walks": walks}


def portal_credentials_from_env() -> "Tuple[Optional[str], Optional[str]]":
    key = os.environ.get("IGVF_ACCESS_KEY") or os.environ.get("IGVF_API_KEY")
    secret = os.environ.get("IGVF_SECRET_ACCESS_KEY") or os.environ.get("IGVF_SECRET_KEY")
    return (key, secret) if key and secret else (None, None)


def parse_cluster_keys(raw: str, clusters: dict, excluded_by_dataset=None, strict: bool = True) -> "set":
    """manage_igvf_metadata.parse_cluster_keys: 'dataset/cluster' or a bare 'dataset' (manifest-eligible expansion)."""
    keys = set()
    for tok in [t.strip() for t in (raw or "").split(",") if t.strip()]:
        d, sep, c = tok.partition("/")
        if not sep or not c:
            known = clusters.get(d)
            if known:
                flagged = {x for x, cfg in known.items() if (cfg or {}).get("igvf_manifest_excluded", False)}
                skip = flagged | set((excluded_by_dataset or {}).get(d, set()))
                keys.update((d, x) for x in known if x not in skip)
                continue
            if strict:
                raise SystemExit(f"--cluster-keys: {tok!r} names no dataset in this config (available: {', '.join(sorted(clusters))})")
            keys.add((d, c))
            continue
        if strict and (d not in clusters or c not in (clusters[d] or {})):
            raise SystemExit(f"--cluster-keys: unknown cluster {tok!r}")
        keys.add((d, c))
    return keys


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def _cfg(args) -> "Tuple[dict, Optional[str]]":
    if not getattr(args, "config", None):
        raise SystemExit("--config is required")
    return load_config(args.config), args.config


def _output_dir(args, config) -> str:
    return getattr(args, "output_dir", None) or config_output_dir(config, getattr(args, "config", None))


def cmd_resolve_exclusions(args) -> int:
    config, path = _cfg(args)
    out = out_dir_for(args, "exclusions")
    inc, up, exc, stats = resolve_exclusions(config, args.legacy_implicit_source)
    man = manifest_eligible_clusters(config, up)
    all_keys = {(d, c) for d, c, _ in iter_clusters(config)}
    target = Path(_output_dir(args, config)) if args.write_to_output_dir else out
    for d in sorted({k[0] for k in all_keys}):
        write_cluster_stats_table(d, stats, target / "cluster_stats")
        write_plan_tsv(d, all_keys, inc, up, man, stats, target / "igvf_metadata")
    reasons = Counter(s["reason"] for s in stats.values())
    lines = [f"Configured clusters: {len(all_keys)}; included {len(inc)}, upload-eligible {len(up)}, manifest-eligible "
             f"{len(man)}, excluded {len(exc)}.", "", md_table(["dataset", "cluster", "reason", "cells", "fragments", "UMIs"],
             [(d, c, s["reason"], s["cell_count"], s["fragments_total"], s["umi_count"]) for (d, c), s in sorted(stats.items())])]
    write_report(out, "Quality gate (resolve_exclusions)", lines,
                 {"included": sorted(map(list, inc)), "upload_eligible": sorted(map(list, up)), "manifest_eligible": sorted(map(list, man)),
                  "excluded": sorted(map(list, exc)), "reasons": dict(reasons), "legacy_implicit_source": args.legacy_implicit_source})
    return 0


def cmd_merge_metrics(args) -> int:
    for p in args.metrics:
        if not os.path.exists(p):
            print(f"ERROR: no such file: {p}")
            return 2
    if os.path.exists(args.out) and not args.force:
        print(f"ERROR: {args.out} exists (use --force)")
        return 2
    rows, overlaps = merge_metrics(args.metrics)
    write_metrics(rows, args.out)
    totals = {c: sum(r[c] for r in rows) for c in METRIC_ADDITIVE}
    print(f"merged {len(args.metrics)} component(s) -> {len(rows)} subsample row(s); totals {totals}; "
          f"{len(overlaps)} overlapping subsample(s) summed")
    return 0


def cmd_prefiltered_metrics(args) -> int:
    row = prefiltered_metrics(args.frag_file, args.subsample_name)
    write_rows_tsv(Path(args.out), METRIC_COLUMNS, [row])
    return 0


def cmd_build_qc_datatables(args) -> int:
    cmap = {}
    if args.config:
        for d, c, cfg in iter_clusters(load_config(args.config)):
            cmap.setdefault(d, {})[c] = [a.strip() for a in str(cfg.get("pseudobulk_annotation", c)).split(",") if a.strip()]
    datasets = args.datasets or sorted(d for d in os.listdir(args.pseudobulks_root)
                                       if os.path.isdir(os.path.join(args.pseudobulks_root, d, "pseudobulks")))
    built = failed = 0
    for d in datasets:
        avail = set(discover_annotations(args.pseudobulks_root, d))
        if not avail:
            continue
        clusters = {c: a for c, a in (cmap.get(d) or {}).items() if set(a) & avail}
        covered = {x for a in clusters.values() for x in a}
        for a in sorted(avail - covered):
            clusters[a] = [a]
        for c, anns in sorted(clusters.items()):
            op = Path(args.dest) / "QC_datatables" / f"{d}_data" / f"{c}_per_cell_qc.tsv"
            if op.exists() and not args.force:
                continue
            try:
                n, ns, skipped = build_one_datatable(args.pseudobulks_root, d, c, anns, op)
            except ValueError as e:
                print(f"{d}/{c}: FAILED -- {e}")
                failed += 1
                continue
            if ns:
                built += 1
                print(f"TSV: {op}  ({n} rows from {ns} subsample file(s))")
    print(f"built {built} datatable(s), {failed} failure(s)")
    return 1 if failed else 0


def cmd_filter_atac(args) -> int:
    res = filter_atac(args.qc_guide, args.pseudobulks, args.cell_type, args.chrom_sizes, args.out,
                      strict_fields=not args.lenient_fields, allow_missing=args.allow_missing_barcodes, clean=args.clean)
    print(f"Wrote: {res['out']}")
    print(f"Wrote: {res['out']}.tbi")
    print(f"retained {res['rows_retained']}/{res['rows_total']} fragments; dropped by chrom {res['dropped_by_chrom']}")
    return 0


def cmd_filter_rna(args) -> int:
    res = filter_rna(args.qc_guide, args.pseudobulks, args.cell_type, args.out, args.ensembl_ids_as_genes, args.gtf,
                     args.standard_chromosomes_only, args.log)
    for p in res["outputs"]:
        print(f"Wrote: {p}")
    print(f"{res['n_cells']} cells x {res['n_genes']} genes ({res['n_ensembl_in']} input features, {res['symbols_summed']} symbols summed)")
    if args.package_tar and res["format"] == "mtx":
        print(f"Wrote: {package_rna_matrix(os.path.dirname(res['outputs'][0]), args.package_tar)}")
    return 0


def cmd_package_rna(args) -> int:
    print(f"Wrote: {package_rna_matrix(args.matrix_dir, args.out)}")
    return 0


def cmd_sce2g_config(args) -> int:
    config, path = _cfg(args)
    out = out_dir_for(args, "sce2g_config")
    inc, up, exc, stats = resolve_exclusions(config)
    odir = _output_dir(args, config)
    thresholds, written = [], []
    for d in sorted({k[0] for k in inc}):
        names = {c for ds, c in inc if ds == d}
        cc = write_cell_clusters_table(d, config["clusters"][d], names, out, os.path.join(odir, "multiome_data", d), config.get("scE2G_dir", ""))
        cm = write_cluster_metadata_table(d, config["clusters"][d], names, out, args.lab_annotations)
        scfg = build_sce2g_config(config, str(cc), os.path.join(odir, "uniformly_processed", d))
        written += [cc, cm, write_json(out / f"{d}_scE2G_config.json", scfg, quiet=True)]
        for c in sorted(names):
            for m in config["clusters"][d][c].get("models") or []:
                thresholds.append({"dataset": d, "cluster": c, "model": m,
                                   "score_threshold": model_threshold(config.get("scE2G_dir"), m, config.get("score_thresholds")),
                                   "primary_model": resolve_primary_model(config["clusters"][d][c].get("models"))})
    for p in written:
        print(f"TSV: {p}" if str(p).endswith(".tsv") else f"JSON: {p}")
    write_rows_tsv(out / "model_thresholds.tsv", ["dataset", "cluster", "model", "score_threshold", "primary_model"], thresholds)
    cmd = sce2g_command(str(path), "", dry_run=True)
    write_report(out, "scE2G configuration", ["scE2G is run by `igvfagent sce2g-pipeline run` or, upstream, by:", "",
                                              "```", " ".join(cmd), "```"], {"thresholds": thresholds, "snakemake": cmd})
    return 0


def cmd_reformat(args) -> int:
    if args.kind == "element-bed":
        out = write_element_bed(args.input, args.out, args.version, args.cell_type, args.term_id, args.summary, args.portal_link)
        print(f"Wrote: {out}")
        print(f"Wrote: {out}.tbi")
        return 0
    if args.kind == "bedpe":
        out = bgzip_index_bedpe(args.input, args.out)
        print(f"Wrote: {out}")
        print(f"Wrote: {out}.tbi")
        return 0
    method = args.method or f"scE2G_{args.model}"
    thr = f"Score >= {args.threshold}" if args.threshold is not None else None
    res = reformat_file(args.input, args.out, method, args.version, args.cell_type, args.term_id, args.summary, thr,
                        args.portal_link, args.all_columns, args.format, args.catlas)
    print(f"Wrote: {res['out']}  ({res['kind']}, {res['rows']} rows)")
    return 0


def cmd_candidates(args) -> int:
    t = candidate_pairs(read_prediction_table(args.input), args.summary)
    print(f"Wrote: {write_table_with_header(t, [], args.out, column_names=True)}")
    return 0


def cmd_features(args) -> int:
    t = feature_table(read_prediction_table(args.input))
    print(f"Wrote: {write_table_with_header(t, feature_header(args.model, args.cell_type, args.term_id, args.summary), args.out)}")
    return 0


def cmd_aggregate_qc(args) -> int:
    out = out_dir_for(args, "aggregate_qc")
    rows = aggregate_qc_rows(args.dataset_dir, args.plots_dir)
    header = list(rows[0].keys()) if rows else []
    dest = Path(args.out) if args.out else out / "all_qc_stats.tsv"
    write_rows_tsv(dest, header, rows)
    figs, warns = ([], [])
    if rows and not args.no_plots:
        figs, warns = qc_plots(dest, out)
    write_report(out, "QC statistics (aggregate_qc_stats)", [f"{len(rows)} (cluster, model) rows from {args.dataset_dir}", ""] +
                 [f"- {w}" for w in warns], {"rows": len(rows), "figures": [str(f) for f in figs], "warnings": warns})
    return 0


def _annotation_rows(args) -> "List[dict]":
    if getattr(args, "annotations_tsv", None):
        return read_rows_tsv(args.annotations_tsv)
    st = State(args.state_db)
    try:
        return st.all_cell_annotations()
    finally:
        st.close()


def cmd_stale_reformats(args) -> int:
    rows = _annotation_rows(args)
    total = []
    for d in args.datasets:
        stale = find_stale(args.results_dir, d, rows)
        for p, diffs in stale:
            print(f"STALE {p}: " + ", ".join(f"{k.strip('# :')}: {a!r} -> {b!r}" for k, (a, b) in diffs.items()))
        total += stale
    print(f"TOTAL STALE: {len(total)} file(s)")
    if total and args.delete:
        for p in remove_stale(total):
            print(f"removed {p}")
    return 0


def cmd_cell_metadata(args) -> int:
    out = out_dir_for(args, "cell_metadata")
    state = State(args.state_db)
    try:
        if args.catlas_seed:
            rows = lab_filter(fetch_multireport(washu_query(), args.multireport_json), WASHU_LAB)
            seeded = seed_catlas(state, rows, args.dataset or "catlas")
            write_rows_tsv(out / "seeded_cell_annotations.tsv", CELL_ANNOTATION_COLUMNS, seeded)
            write_report(out, "CATlas Cell Annotation seed", [f"{len(seeded)} WashU cluster(s) seeded"], {"seeded": len(seeded)})
            return 0
        config, path = _cfg(args)
        fetched = None
        if args.multireport_json or (not args.offline and is_stale(state.latest_fetch(), args.ttl_hours)):
            fetched = cache_multireport(state, fetch_multireport(MULTIREPORT_QUERY, args.multireport_json))
        cfgs = {(d, c): cfg for d, c, cfg in iter_clusters(config) if "qc_guide" in cfg}
        statuses = derive_scopes(state, cfgs)
        digest = cluster_set_digest({(d, c) for d, c, _ in iter_clusters(config)})
        snap_root = Path(_output_dir(args, config)) if args.write_snapshot else out
        for d in sorted({k[0] for k in cfgs}):
            rows = [r for r in state.all_cell_annotations() if r["dataset"] == d]
            print(f"TSV: {write_snapshot(snapshot_path(snap_root, d), rows, state.latest_fetch() or '', digest)}")
            write_rows_tsv(snap_root / "igvf_metadata" / f"{d}_cell_annotation_status.tsv",
                           ["dataset", "cluster", "resolved", "reason", "cell_annotation"],
                           [{**s, "resolved": "y" if s["resolved"] else "n"} for s in statuses if s["dataset"] == d])
        write_rows_tsv(out / "cell_annotation_table.tsv", CELL_ANNOTATION_COLUMNS[:2] + CELL_ANNOTATION_COLUMNS[2:7] +
                       ["all_primary_released", "principal_uploaded", "principal_alias"], shareable_rows(state))
        res = sum(1 for s in statuses if s["resolved"])
        write_report(out, "Cell Annotation cache", [f"derived {res}/{len(statuses)} scope(s); fetch: {fetched or 'cache reused'}", "",
                     md_table(["dataset", "cluster", "resolved", "reason"], [(s["dataset"], s["cluster"], s["resolved"], s["reason"]) for s in statuses])],
                     {"statuses": statuses, "fetch": fetched, "digest": digest})
        return 0
    finally:
        state.close()


def lab_filter(rows, lab: str) -> "List[dict]":
    """The multireport's lab.@id filter, re-applied locally (a saved --multireport-json holds every lab)."""
    def lab_of(r):
        v = r.get("lab")
        return v.get("@id") if isinstance(v, dict) else v
    return [r for r in rows if lab_of(r) == lab]


def washu_query(lab: str = WASHU_LAB) -> str:
    from urllib.parse import quote
    return (f"type=PseudobulkSet&status%21=deleted&lab.%40id={quote(lab, safe='')}&limit=all"
            + "".join(f"&field={quote(f)}" for f in WASHU_FIELDS))


def cmd_cell_annotation_report(args) -> int:
    out = out_dir_for(args, "cell_annotation_report")
    rows, skipped = cell_annotation_report_rows(fetch_multireport(MULTIREPORT_QUERY, args.multireport_json), args.qc_guide_dir)
    dest = Path(args.out) if args.out else out / "cell_annotations_by_dataset_cluster.tsv"
    write_rows_tsv(dest, ANNOTATION_REPORT_COLUMNS, rows)
    resolved = sum(1 for r in rows if r["QCGuideFile"])
    multi = sum(1 for r in rows if REPORT_JOIN in r["CellAnnotation"])
    write_report(out, "IGVF cell annotation report", [f"{len(rows)} dataset-cluster rows; {resolved} resolved via a QC guide; "
                 f"{multi} with a multi-value CellAnnotation"], {"rows": len(rows), "resolved": resolved, "multi_value": multi, "skipped": skipped})
    return 0


def cmd_dataset_accessions(args) -> int:
    mapping, skipped = dataset_accession_mapping(fetch_multireport(
        "type=PseudobulkSet&status%21=deleted&limit=all&field=%40id&field=aliases&field=input_file_sets"
        "&field=input_file_sets.%40id&field=input_file_sets.accession&field=samples", args.multireport_json))
    multi = {d: a for d, a in mapping.items() if len(a) > 1}
    for d, a in multi.items():
        print(f"WARNING: {d} maps to {len(a)} accessions {a}")
    write_json(Path(args.json_out) if args.json_out else out_dir_for(args, "dataset_accessions") / "dataset_to_principal_analysis_set_accession.json", mapping)
    return 0


def cmd_washu_report(args) -> int:
    out = out_dir_for(args, "washu_report")
    rows = [washu_row(r) for r in lab_filter(fetch_multireport(washu_query(args.lab), args.multireport_json), args.lab)]
    write_rows_tsv(Path(args.out) if args.out else out / "washu_pseudobulk_report.tsv", WASHU_HEADER, rows)
    return 0


def cmd_verify_fragments(args) -> int:
    rows = read_rows_tsv(args.report)
    os.makedirs(args.outdir, exist_ok=True)
    results = []
    for r in rows:
        href, md5 = r.get("FragmentsFile_Href"), r.get("FragmentsFile_MD5")
        if not href or not md5:
            results.append({**r, "LocalPath": "", "ActualMD5": "", "Match": "SKIPPED_NO_HREF_OR_MD5"})
            continue
        _, fn = accession_from_href(href)
        dest = os.path.join(args.outdir, fn)
        if os.path.exists(dest) and not args.force and md5_file(dest) == md5:
            results.append({**r, "LocalPath": dest, "ActualMD5": md5, "Match": "PASS_CACHED"})
            continue
        if not args.download:
            results.append({**r, "LocalPath": dest, "ActualMD5": md5_file(dest) if os.path.exists(dest) else "",
                            "Match": "FAIL_MD5_MISMATCH" if os.path.exists(dest) else "NOT_DOWNLOADED"})
            continue
        try:
            from raw_data_pipeline import portal_download  # type: ignore
            portal_download(href, Path(dest))
        except Exception as exc:  # noqa: BLE001
            results.append({**r, "LocalPath": dest, "ActualMD5": "", "Match": "DOWNLOAD_FAILED", "error": redact(exc)})
            continue
        got = md5_file(dest)
        results.append({**r, "LocalPath": dest, "ActualMD5": got, "Match": "PASS" if got == md5 else "FAIL_MD5_MISMATCH"})
    cols = (list(rows[0].keys()) if rows else []) + ["LocalPath", "ActualMD5", "Match"]
    write_rows_tsv(Path(args.summary) if args.summary else Path(args.outdir) / "fragments_download_verification.tsv", cols, results)
    bad = [r for r in results if r["Match"] not in ("PASS", "PASS_CACHED")]
    print(f"{len(results) - len(bad)}/{len(results)} verified")
    return 1 if bad else 0


def cmd_portal_files(args) -> int:
    out = out_dir_for(args, "portal_files")
    q = "&".join([("type=PseudobulkSet"), "status%21=deleted", "limit=all"] + [f"field={f.replace('@', '%40')}" for f in PORTAL_FILES_FIELDS])
    recs, rep = discover_portal_files(fetch_multireport(q, args.multireport_json), args.lab or None,
                                      args.content_types, set(args.datasets) if args.datasets else None)
    plan = []
    for r in recs:
        r["local_path"] = os.path.join(args.download_root, r["rel_path"]) if (args.download_root and r["rel_path"]) else ""
        need, why = needs_download(r, r["local_path"]) if r["local_path"] else (False, "needs_review: no local path")
        r["download_reason"], r["download_state"] = why, ("todo" if need else ("done" if why == "unchanged" else "needs_review"))
        r["review_reasons"] = ",".join(r["review_reasons"])
        plan.append(r)
    if args.download:
        from raw_data_pipeline import portal_download  # type: ignore
        for r in plan:
            if r["download_state"] != "todo":
                continue
            try:
                portal_download(r["href"], Path(r["local_path"]))
                ok = not r.get("md5sum") or md5_file(r["local_path"]) == r["md5sum"]
                r["download_state"] = "done" if ok else "md5_mismatch"
            except Exception as exc:  # noqa: BLE001
                r["download_state"], r["download_reason"] = "failed", redact(exc)
    cols = ["accession", "dataset", "annotation", "subsample", "content_type", "file_format", "file_size", "md5sum", "href",
            "portal_status", "upload_status", "rel_path", "local_path", "download_state", "download_reason", "review_reasons",
            "set_alias", "principal_analysis_set", "submitted_file_name"]
    write_rows_tsv(out / "portal_files.tsv", cols, plan)
    write_report(out, "Primary-pseudobulk file discovery", [f"{len(plan)} file(s); {sum(1 for r in plan if r['download_state'] == 'todo')} to download"],
                 {"report": rep, "states": dict(Counter(r["download_state"] for r in plan))})
    return 0


def cmd_compare_archive(args) -> int:
    out = out_dir_for(args, "compare_archive")
    rows = read_rows_tsv(args.discovery)
    results = []
    for r in rows:
        lp = r.get("local_path") or ""
        rec = {k: r.get(k, "") for k in ("dataset", "annotation", "subsample", "content_type", "accession")}
        rec.update({"dirname": os.path.basename(os.path.dirname(lp)) if lp else "", "portal_path": lp, "archive_path": ""})
        if not lp or not os.path.exists(lp):
            rec.update({"status": "not_downloaded", "verdict": "", "detail": f"download_state={r.get('download_state')}"})
        else:
            arch = archive_counterpart(args.archive_root, rec["dataset"], rec["dirname"], os.path.basename(lp))
            if not arch:
                rec.update({"status": "novel", "verdict": "", "detail": "no counterpart in archive"})
            else:
                st, vd, det = compare_file(lp, arch, rec["content_type"], not args.no_characterize)
                rec.update({"archive_path": arch, "status": st, "verdict": vd, "detail": "; ".join(f"{k}={v}" for k, v in det.items())})
        results.append(rec)
    cols = ["dataset", "annotation", "subsample", "dirname", "content_type", "accession", "status", "verdict", "portal_path", "archive_path", "detail"]
    write_rows_tsv(out / "portal_vs_archive.tsv", cols, results)
    by = defaultdict(list)
    for r in results:
        by[(r["dataset"], r["annotation"] or r["dirname"])].append(r)
    roll = [{"dataset": d, "annotation": a, "n_files": len(v), "n_subsamples": len({x["subsample"] for x in v}),
             "statuses": ",".join(sorted({x["status"] for x in v})), "verdicts": ",".join(sorted({x["verdict"] for x in v if x["verdict"]})),
             "decision": cluster_decision({x["status"] for x in v})} for (d, a), v in sorted(by.items())]
    write_rows_tsv(out / "portal_vs_archive_by_cluster.tsv", list(roll[0].keys()) if roll else ["dataset"], roll)
    write_report(out, "Portal vs archive", [md_table(["dataset", "annotation", "decision"], [(r["dataset"], r["annotation"], r["decision"]) for r in roll])],
                 {"decisions": dict(Counter(r["decision"] for r in roll)), "files": dict(Counter(f"{r['status']}/{r['verdict']}" for r in results))})
    return 0


def make_target_portal(production: bool):
    from igvf_submission_skill import Portal  # type: ignore
    key, secret = portal_credentials_from_env()
    if not key:
        raise SystemExit("refusing to contact the Portal for writing: set IGVF_ACCESS_KEY and IGVF_SECRET_ACCESS_KEY "
                         "(or IGVF_API_KEY / IGVF_SECRET_KEY) in the environment")
    return Portal(mode="prod" if production else "sandbox", api_key=key, secret_key=secret)


def ctx_factory(config: dict, state: State, output_dir: str):
    igvf = igvf_config(config)
    clusters = config.get("clusters") or {}

    def make(d, c, model):
        return Ctx(d, c, model, clusters[d][c], igvf, config.get("scE2G_dir"), config.get("data_dir"), output_dir, state,
                   config.get("score_thresholds"), config.get("results_dir_base"))
    return make, igvf


def run_manifest(config: dict, config_path, cluster_keys, excluded_keys, state_db, manifest_dir: Path, out: Path,
                 reader=None, execute_client=None, lineage_fetch=None, lineage_walk: bool = False, table_names=None,
                 output_dir: "Optional[str]" = None, target: str = "sandbox") -> dict:
    """Plan (always), write every artefact, and only when `execute_client` is given submit in round order."""
    state = State(state_db)
    try:
        for d, c in excluded_keys:
            state.mark_excluded(d, c, "resolve_exclusions")
        cluster_keys = set(cluster_keys) - set(excluded_keys)
        clusters = config.get("clusters") or {}
        cfgs = {(d, c): clusters[d][c] for d, c in cluster_keys}
        make, igvf = ctx_factory(config, state, output_dir or config_output_dir(config, config_path))
        reader = reader or OfflineReader()
        plan = plan_manifest(cluster_keys, cfgs, make, state, reader, igvf, table_names)
        files = write_manifest_files(plan, manifest_dir)
        for f in files:
            print(f"TSV: {f['path']}")
        own = {e["alias"] for e in plan["posts"] + plan["patches"]} | {r["payload"]["aliases"][0] for r in plan["records"]}
        allp = [e["payload"] for e in plan["posts"] + plan["patches"]] + [r["payload"] for r in plan["records"]]
        lineage = None
        if lineage_fetch is not None:
            acc = [a for a in ((config.get("igvf") or {}).get("principal_analysis_sets") or {}).values()] if lineage_walk else []
            lineage = check_lineage(external_references(allp, own), lineage_fetch, acc)
            write_json(out / "lineage.json", lineage)
        mode = "prod" if target == "production" else target
        plan_rows = []
        for kind, entries in (("post", plan["posts"]), ("patch", plan["patches"])):
            for e in entries:
                body = rest_payload(e["payload"], e["object_type"], compute_md5=False)
                if isinstance(body.get("attachment"), dict) and "href" in body["attachment"]:
                    body["attachment"] = {**body["attachment"], "href": body["attachment"]["href"][:48] + "..."}
                plan_rows.append({"round": e["round"], "kind": kind, "object_type": e["object_type"], "table": e["table"],
                                  "variant": e["variant"], "dataset": e["dataset"], "cluster": e["cluster"], "model": e["model"],
                                  "alias": e["alias"], "record_id": e.get("record_id"), "waiting_on": e["deferred_on"],
                                  "collection": f"/{COLLECTIONS[e['object_type']]}/", "rest_payload": body})
        plan_rows.sort(key=lambda r: (r["round"], r["dataset"], r["table"], r["variant"], r["cluster"]))
        iu_cmds = [f"iu_register.py -m {mode} -p {f['object_type']} -i {f['path']}" + (" --patch -w" if f["kind"] == "patch" else "")
                   for f in files]
        write_json(out / "upload_plan.json", {"target": target, "executed": execute_client is not None, "rows": plan_rows,
                                               "iu_register_equivalent": iu_cmds})
        with open(out / "payloads.jsonl", "w") as fh:
            for e in plan["posts"] + plan["patches"]:
                fh.write(json.dumps({"alias": e["alias"], "object_type": e["object_type"], "payload": e["payload"]}) + "\n")
        print(f"Wrote: {out / 'payloads.jsonl'}")
        executed = None
        if execute_client is not None:
            log_rows = []
            executed = execute_plan(plan, execute_client, state, reader, log_rows)
            write_rows_tsv(out / "execute_log.tsv", ["alias", "kind", "result", "detail"], log_rows)
        outcomes = Counter(r["outcome"] for r in plan["coverage"])
        return {"plan": plan, "files": [str(f["path"]) for f in files], "outcomes": dict(outcomes), "lineage": lineage,
                "executed": executed, "gaps": sum(outcomes.get(o, 0) for o in MANIFEST_GAP_OUTCOMES),
                "pending": sum(1 for e in plan["posts"] + plan["patches"] if e["deferred_on"])}
    finally:
        state.close()


def _manifest_report(out: Path, res: dict, title: str, executed: bool, target: str) -> None:
    oc = res["outcomes"]
    rows = sorted(Counter((r["table"], r["variant"], r["outcome"]) for r in res["plan"]["coverage"]).items())
    lines = [f"Mode: {'EXECUTE against ' + target if executed else 'DRY-RUN (nothing sent to the Portal)'}.", "",
             f"Outcomes: {dict(oc)}; manifest gaps: {res['gaps']}; rows waiting on a dependency: {res['pending']}.", "",
             md_table(["table", "variant", "outcome", "rows"], [(t, v or "(default)", o, n) for (t, v, o), n in rows], 80), "",
             "Round files (sorting by name is the upload order):", ""] + [f"- `{Path(f).name}`" for f in res["files"]]
    if res.get("lineage"):
        lines += ["", f"Lineage: {len(res['lineage']['resolved'])} external reference(s) resolved, "
                  f"{len(res['lineage']['missing'])} missing, {len(res['lineage']['blocked'])} blocked (403)."]
    if res.get("executed"):
        ex = res["executed"]
        lines += ["", f"Uploaded {len(ex['uploaded'])}, failed {len(ex['failed'])}, deferred {len(ex['deferred'])}."]
    write_report(out, title, lines, {"outcomes": oc, "gaps": res["gaps"], "pending": res["pending"], "files": res["files"],
                                     "executed": res.get("executed"), "lineage_missing": (res.get("lineage") or {}).get("missing"),
                                     "target": target, "dry_run": not executed})


def cmd_manifest(args) -> int:
    config, path = _cfg(args)
    out = out_dir_for(args, "manifest")
    clusters = config.get("clusters") or {}
    state_db = args.state_db or (config.get("igvf") or {}).get("state_db_path") or str(DEFAULT_STATE_DB)
    st = State(state_db)
    excluded_by = {d: st.excluded_clusters(d) for d in clusters}
    st.close()
    keys = parse_cluster_keys(args.cluster_keys, clusters, excluded_by)
    excl = parse_cluster_keys(args.excluded_cluster_keys or "", clusters, strict=False)
    target = "production" if args.production else "sandbox"
    reader, client = OfflineReader(json.loads(Path(args.live_aliases).read_text()) if args.live_aliases else None), None
    if args.execute:
        client = make_target_portal(args.production)
        reader = PortalAliasReader(client)
        print(f"EXECUTE: submitting to {client.base} ({target})")
    elif not args.offline and portal_credentials_from_env()[0]:
        try:
            from igvf_submission_skill import Portal  # type: ignore
            k, s = portal_credentials_from_env()
            reader = PortalAliasReader(Portal(mode="prod" if args.production else "sandbox", api_key=k, secret_key=s))
        except Exception as e:  # noqa: BLE001
            print(f"alias lookups offline ({e})")
    lineage_fetch = None
    if args.lineage_fixture:
        fx = json.loads(Path(args.lineage_fixture).read_text())
        lineage_fetch = lambda p: (200, fx[p.split("?")[0]]) if p.split("?")[0] in fx else (404, None)  # noqa: E731
    elif args.lineage:
        lineage_fetch = portal_get
    manifest_dir = Path(args.manifest_dir) if args.manifest_dir else out / "igvf_manifests"
    res = run_manifest(config, path, keys, excl, state_db, manifest_dir, out, reader, client, lineage_fetch, args.lineage_walk,
                       args.tables or None, target=target)
    _manifest_report(out, res, "IGVF Portal manifests", client is not None, target)
    if res.get("executed") and res["executed"]["failed"]:
        return 1
    return 1 if res["gaps"] else 0


def cmd_patch_submitter_comment(args) -> int:
    out = out_dir_for(args, "submitter_comment")
    rows = fetch_multireport(PREDICTION_SET_QUERY, args.multireport_json)
    to_patch, counts = plan_submitter_comment(rows, set(args.dataset or []))
    f = write_iu_tsv(out / "prediction_set_submitter_comment_patch.tsv", [{"submitter_comment": e["new"]} for e in to_patch],
                     [e["record_id"] for e in to_patch])
    if f:
        print(f"TSV: {f}")
    results = []
    if args.execute and to_patch:
        client = make_target_portal(args.production)
        for e in to_patch:
            try:
                client._request("PATCH", f"/{str(e['record_id']).strip('/')}/", payload={"submitter_comment": e["new"]})
                results.append({"alias": e["alias"], "result": "patched"})
            except Exception as exc:  # noqa: BLE001
                results.append({"alias": e["alias"], "result": f"error: {exc}"})
    write_report(out, "Prediction Set submitter_comment backfill", [f"plan: {counts}; {len(to_patch)} to patch; "
                 f"{'executed' if args.execute else 'DRY-RUN'}"], {"counts": counts, "to_patch": to_patch, "results": results})
    return 0


def cmd_report(args) -> int:
    config, path = _cfg(args)
    out = out_dir_for(args, "coverage_report")
    odir = _output_dir(args, config)
    ann = {(r["dataset"], r["cluster"]): r for r in _annotation_rows(args)}
    rows = coverage_report_rows(odir, config, ann, args.manifest_dir or os.path.join(odir, "igvf_manifests"))
    write_rows_tsv(Path(args.out) if args.out else out / "report.tsv", REPORT_HEADER, rows)
    write_report(out, "Coverage report", [f"{len(rows)} cluster row(s)", "", md_table(["status", "n"], Counter(r["status"] for r in rows).most_common())],
                 {"status": dict(Counter(r["status"] for r in rows)), "manifest_status": dict(Counter(r["manifest_status"] for r in rows))})
    return 0


def cmd_synapse_manifest(args) -> int:
    out = out_dir_for(args, f"synapse_{args.product}")
    keys = parse_cluster_keys(args.cluster_keys, {}, strict=False)
    owned = read_inventory(args.inventory) if args.inventory else {}
    preserved = load_preserved_rows(args.manifest_out, keys)
    plan = synapse_plan(args.product, args.files or [], owned, keys, preserved, args.parent_id, args.legacy_first_match)
    dest = Path(args.manifest_out) if args.manifest_out else out / f"{args.product}_manifest.tsv"
    write_rows_tsv(dest, ["path", "parent"], plan["manifest_rows"])
    for (d, c, n), meta in plan["to_delete"].items():
        print(f"STALE/EXCLUDED (would delete): {d}/{c}/{n} ({meta.get('id')})")
    for (d, c, n) in plan["to_overwrite"]:
        print(f"WOULD OVERWRITE: {d}/{c}/{n}")
    if args.execute:
        if plan["abort_overwrite"] and not args.confirm_overwrite:
            print("confirm-overwrite is False for the shared 'Multiome datasets' folder -- ABORTING real sync")
            return 1
        try:
            import synapseclient  # type: ignore
            import synapseutils  # type: ignore
        except ImportError:
            print("synapseclient not installed: nothing synced (pip install synapseclient)")
            return 1
        syn = synapseclient.login()
        if args.confirm_delete:
            for meta in plan["to_delete"].values():
                syn.delete(meta["id"])
        synapseutils.syncToSynapse(syn, manifestFile=str(dest), dryRun=False)
    write_report(out, f"Synapse manifest ({args.product})", [f"{len(plan['to_upload_new'])} new, {len(plan['to_overwrite'])} overwrite, "
                 f"{len(plan['to_delete'])} stale/excluded; {'EXECUTED' if args.execute else 'dry-run'}"],
                 {"new": len(plan["to_upload_new"]), "overwrite": len(plan["to_overwrite"]), "stale": len(plan["to_delete"]),
                  "unmatched": plan["unmatched"], "preserved": len(preserved)})
    return 0


def cmd_synapse_orphans(args) -> int:
    out = out_dir_for(args, "synapse_orphans")
    inv = read_rows_tsv(args.inventory)
    if args.mode == "cluster":
        in_scope = {(r["dataset"], r["cluster"]) for r in read_rows_tsv(args.in_scope_tsv)}
        rows = [r for r in inv if (r["dataset"], r["cluster"]) not in in_scope]
        for r in rows:
            r["uploaded_by_me"] = "y" if (args.user_id and r.get("created_by") == args.user_id and r.get("modified_by") == args.user_id) else "n"
        cols = ["dataset", "cluster", "name", "id", "created_by", "modified_by", "uploaded_by_me"]
    else:
        covered = set()
        for r in read_rows_tsv(args.manifest):
            parts = r["path"].rstrip("/").split("/")
            covered.add((parts[-3], parts[-2], parts[-1]))
        only = {tuple(t.split("/", 1)) for t in args.only or []}
        excl = set(args.exclude_top or PREDICTION_ORPHANS_EXCLUDE_TOP)
        rows = [r for r in inv if ((r["dataset"], r["cluster"]) in only if only else r["dataset"] not in excl)
                and (only or (r["dataset"], r["cluster"], r["name"]) not in covered)
                and (not args.user_id or r.get("created_by") == args.user_id)]
        cols = ["dataset", "cluster", "name", "id", "created_by", "modified_by"]
    write_rows_tsv(out / "synapse_orphans.tsv", cols, rows)
    write_report(out, "Synapse orphans", [f"{len(rows)} orphan row(s) ({args.mode} mode); read-only, nothing deleted"], {"orphans": len(rows)})
    return 0


def cmd_distance_depth(args) -> int:
    out = out_dir_for(args, "distance_depth")
    stats = [r for p in args.stats for r in read_rows_tsv(p)]
    pmap = {r["cluster"]: r["path"] for r in read_rows_tsv(args.predictions_map)} if args.predictions_map else {}
    gmap = {r["cluster"]: r["group"] for r in read_rows_tsv(args.group_map)} if args.group_map else {}
    for r in stats:
        r["group"] = gmap.get(r["cluster"], r.get("group", ""))
    rows = distance_depth_table(stats, pmap)
    cols = ["cluster", "model_name", "group", "fragments_total", "cell_count", "frag_per_cell", "mean_dist_to_tss", "median_dist_to_tss",
            "q25_dist_to_tss", "q75_dist_to_tss", "iqr_dist_to_tss", "neuron", "above_threshold"]
    write_rows_tsv(out / "distance_vs_depth.tsv", cols, rows)
    trends = {y: log_trend(rows, y) for y in ("mean_dist_to_tss", "median_dist_to_tss", "iqr_dist_to_tss")}
    figs = []
    plt = _plt()
    if plt is not None and not args.no_plots and rows:
        for i, (y, lab) in enumerate((("mean_dist_to_tss", "Mean distance to TSS (bp)"), ("median_dist_to_tss", "Median distance to TSS (bp)"),
                                      ("iqr_dist_to_tss", "IQR of distance to TSS (bp)")), 1):
            fig, ax = plt.subplots(figsize=(6, 4))
            for above in (True, False):
                sub = [r for r in rows if r["above_threshold"] == above and math.isfinite(r[y])]
                ax.scatter([r["fragments_total"] for r in sub], [r[y] for r in sub], s=14, alpha=1.0 if above else 0.25,
                           c=["#792374" if r["neuron"] else "#006479" for r in sub])
            if trends[y]:
                sl, ic, _ = trends[y]
                xs = [DEPTH_THRESHOLD, max(r["fragments_total"] for r in rows)]
                ax.plot(xs, [sl * math.log10(x) + ic for x in xs], color="#52514e", lw=1)
            ax.axvline(DEPTH_THRESHOLD, color="#96a0b3", ls="--", lw=0.8)
            ax.set_xscale("log")
            ax.set_xlabel("Total ATAC fragments")
            ax.set_ylabel(lab)
            p = out / f"fig{i}_{y}_vs_fragments.png"
            fig.savefig(p, bbox_inches="tight", dpi=120)
            plt.close(fig)
            print(f"Figure: {p}")
            figs.append(str(p))
    write_report(out, "Distance to TSS vs sequencing depth", [f"{len(rows)} cluster rows; trends (slope per log10 fragment, clusters "
                 f">= 2e6 fragments): {trends}"], {"trends": trends, "figures": figs,
                 "neuron_vs_non": {k: len([r for r in rows if r['neuron'] == k]) for k in (True, False)}})
    return 0


# ---------------------------------------------------------------------------
# run: the driver (run_pipeline.py): preflight, warm, local packaging, manifest preview, audit
# ---------------------------------------------------------------------------

def package_cluster(config, dataset, cluster, output_dir, ann_row, cm_row, portal_alias, notes: list, reformat: bool) -> "List[str]":
    """The post-scE2G rules for one cluster, run only where the input exists and the output does not."""
    cfg = config["clusters"][dataset][cluster]
    results = os.path.join(output_dir, "uniformly_processed")
    cdir = os.path.join(results, dataset, cluster)
    made, version = [], str(config.get("scE2G_version", ""))
    md = os.path.join(output_dir, "multiome_data", dataset, cluster)
    rna_dir = os.path.join(md, f"rna_count_matrix_{dataset}_{cluster}")
    tar = os.path.join(md, f"rna_count_matrix_{dataset}_{cluster}.tar.gz")
    if not is_atac_only(cfg) and os.path.isdir(rna_dir) and not os.path.exists(tar):
        made.append(str(package_rna_matrix(rna_dir, tar)))
    for model in cfg.get("models") or []:
        t = model_threshold(config.get("scE2G_dir"), model, config.get("score_thresholds"))
        raw = os.path.join(cdir, model, f"scE2G_predictions_threshold{t}.bedpe")
        gz = os.path.join(cdir, f"{dataset}_{cluster}_scE2G_{model}_threshold{t}.bedpe.gz")
        if os.path.exists(raw) and not os.path.exists(gz):
            made.append(str(bgzip_index_bedpe(raw, gz)))
    primary = resolve_primary_model(cfg.get("models"))
    ppred = os.path.join(cdir, primary, "scE2G_predictions.tsv.gz")
    if os.path.exists(ppred):
        cand = os.path.join(results, "candidate_e2g_pairs", f"{dataset}_{cluster}_candidate_e2g_pairs.tsv.gz")
        if not os.path.exists(cand):
            made.append(str(write_table_with_header(candidate_pairs(read_prediction_table(ppred), cm_row.get("summary", cluster)), [], cand)))
        feat = os.path.join(results, "scE2G_features", dataset, cluster, f"{dataset}_{cluster}_scE2G_{primary}_features.tsv.gz")
        if not os.path.exists(feat):
            made.append(str(write_table_with_header(feature_table(read_prediction_table(ppred)),
                                                    feature_header(primary, cm_row.get("cell_type"), cm_row.get("ontology_id"), cm_row.get("summary")), feat)))
    else:
        notes.append(f"{dataset}/{cluster}: no scE2G predictions yet ({ppred}) -- run `igvfagent sce2g-pipeline run` first")
    if not reformat or ann_row is None:
        return made
    name, tid, summ = ann_row.get("term_name"), ann_row.get("term_id"), ann_row.get("cell_annotation")
    for model in cfg.get("models") or []:
        t = model_threshold(config.get("scE2G_dir"), model, config.get("score_thresholds"))
        jobs = [(os.path.join(cdir, model, "scE2G_predictions.tsv.gz"), os.path.join(cdir, f"{dataset}_{cluster}_scE2G_{model}.e2g.tsv.gz"), None, "full"),
                (os.path.join(cdir, model, f"scE2G_predictions_threshold{t}.tsv.gz"),
                 os.path.join(cdir, f"{dataset}_{cluster}_scE2G_{model}_threshold{t}.e2g.tsv.gz"), f"Score >= {t}", "thresholded")]
        if model == MULTIOME_MODEL:
            jobs.append((os.path.join(cdir, model, "scE2G_gene_list.tsv.gz"), os.path.join(cdir, f"{dataset}_{cluster}_scE2G_multiome_v3_gene_list.tsv.gz"), None, "genes"))
        if model == SCATAC_MODEL and config.get("catlas"):
            for meta, var in (("element", "elements"), ("gene", "genes")):
                jobs.append((os.path.join(cdir, model, f"scE2G_{meta}_list.tsv.gz"),
                             os.path.join(cdir, f"{dataset}_{cluster}_scE2G_scATAC_powerlaw_v3_{meta}_list.tsv.gz"), None, var))
        for src, dst, thr, var in jobs:
            if os.path.exists(src) and not os.path.exists(dst):
                reformat_file(src, dst, f"scE2G_{model}", version, name, tid, summ, thr, portal_alias(model, var),
                              catlas=bool(config.get("catlas")))
                made.append(dst)
    el = os.path.join(cdir, "Neighborhoods", "EnhancerList.bed")
    eb = os.path.join(cdir, f"{dataset}_{cluster}_element_list.bed.gz")
    if os.path.exists(el) and not os.path.exists(eb):
        made.append(str(write_element_bed(el, eb, version, name, tid, summ, portal_alias(MULTIOME_MODEL, "elements_bed"))))
    return made


def cmd_run(args) -> int:
    config, path = _cfg(args)
    try:
        mode = args.mode or config.get("pipeline_mode", "default")
        if mode not in ("default", "local_only"):
            raise ValueError(f"pipeline_mode must be 'default' or 'local_only', got {mode!r}")
    except ValueError as e:
        print(f"PREFLIGHT ERROR: {e}")
        return 2
    out = out_dir_for(args, "run")
    odir = _output_dir(args, config)
    lines, degraded, problems, notes = [], [], [], []
    # [0] preflight
    inc, up, exc, stats = resolve_exclusions(config, args.legacy_implicit_source)
    man = manifest_eligible_clusters(config, up)
    all_keys = {(d, c) for d, c, _ in iter_clusters(config)}
    datasets = sorted({k[0] for k in all_keys})
    for d in datasets:
        write_cluster_stats_table(d, stats, Path(odir) / "cluster_stats")
        write_plan_tsv(d, all_keys, inc, up, man, stats, Path(odir) / "igvf_metadata")
    lines += ["## [0] preflight", "", f"mode {mode}; configured {len(all_keys)}, included {len(inc)}, upload-eligible {len(up)}, "
              f"manifest-eligible {len(man)}, excluded {len(exc)}", ""]
    # [1] warm
    state_db = args.state_db or (config.get("igvf") or {}).get("state_db_path") or str(DEFAULT_STATE_DB)
    annotations = {}
    snap_paths = []
    if mode == "default":
        state = State(state_db)
        try:
            cfgs = {(d, c): cfg for d, c, cfg in iter_clusters(config) if "qc_guide" in cfg}
            if args.dry_run:  # upstream: a dry run neither fetches nor derives; it snapshots what is cached
                statuses = [{"dataset": d, "cluster": c, "cell_annotation": None, "reason": "dry-run: not derived",
                             "resolved": state.get_cell_annotation(*annotation_lookup_key(d, c, cfg)) is not None}
                            for (d, c), cfg in sorted(cfgs.items())]
            else:
                try:
                    if args.multireport_json or (not args.offline and is_stale(state.latest_fetch())):
                        cache_multireport(state, fetch_multireport(MULTIREPORT_QUERY, args.multireport_json))
                except SystemExit as e:
                    degraded.append(f"portal fetch failed ({e}); continuing with the existing raw cache")
                statuses = derive_scopes(state, cfgs)
            digest = cluster_set_digest(all_keys)
            for d in datasets:
                rows = [r for r in state.all_cell_annotations() if r["dataset"] == d]
                snap_paths.append(write_snapshot(snapshot_path(odir, d), rows, state.latest_fetch() or "", digest))
                stale = find_stale(os.path.join(odir, "uniformly_processed"), d, rows)
                if stale:
                    notes.append(f"{d}: {len(stale)} reformatted file(s) carry a superseded Cell Annotation" +
                                 (" (dry run: not deleted)" if args.dry_run else " -- deleted for regeneration"))
                    if not args.dry_run:
                        remove_stale(stale)
                write_rows_tsv(Path(odir) / "igvf_metadata" / f"{d}_cell_annotation_status.tsv",
                               ["dataset", "cluster", "resolved", "reason", "cell_annotation"],
                               [{**s, "resolved": "y" if s["resolved"] else "n"} for s in statuses if s["dataset"] == d], quiet=True)
            for d in datasets:
                annotations.update(read_snapshot(snapshot_path(odir, d), digest, None))
            lines += ["## [1] warm", "", f"{sum(1 for s in statuses if s['resolved'])}/{len(statuses)} scope(s) resolved", ""]
        finally:
            state.close()
    reformat_eligible = {(d, c) for d, c in up if annotation_lookup_key(d, c, config["clusters"][d][c]) in annotations}
    manifest_ready = man & reformat_eligible
    # [2] local packaging (scE2G itself: sce2g-pipeline)
    made = []
    try:
        cm_rows = {}
        if not args.dry_run:
            for d in sorted({k[0] for k in inc}):
                names = {c for ds, c in inc if ds == d}
                p = write_cluster_metadata_table(d, config["clusters"][d], names, Path(odir) / "config", args.lab_annotations)
                cm_rows.update({(d, r["cluster"]): r for r in read_rows_tsv(p)})
        make, igvf = ctx_factory(config, None, odir)
        for d, c in sorted(up):
            if args.dry_run:
                continue
            ann = annotations.get(annotation_lookup_key(d, c, config["clusters"][d][c]))
            alias = lambda model, var, d=d, c=c: TABLES["prediction_tabular_files"].alias(make(d, c, model), var)  # noqa: E731
            made += package_cluster(config, d, c, odir, ann, cm_rows.get((d, c), {"summary": c}), alias, notes,
                                    (d, c) in reformat_eligible and mode == "default")
        for d in datasets:
            ddir = os.path.join(odir, "uniformly_processed", d)
            first = next(iter(config["clusters"][d].values()))
            rows = aggregate_qc_rows(ddir, os.path.dirname(os.path.dirname(first.get("qc_guide", ""))) if first.get("qc_guide") else "")
            if rows and not args.dry_run:
                write_rows_tsv(Path(ddir) / "qc_plots" / "all_qc_stats.tsv", list(rows[0].keys()), rows)
    finally:
        for p in snap_paths:
            if os.path.exists(p):
                os.remove(p)
    lines += ["## [2] local packaging", "", f"{len(made)} file(s) written; scE2G command (not run here): `{' '.join(sce2g_command(str(path), ''))}`", ""]
    lines += [f"- {n}" for n in notes]
    # [3] manifest preview + [4] audit
    res = None
    if mode == "default" and not args.dry_run and not any("portal fetch failed" in x for x in degraded):
        if manifest_ready:
            res = run_manifest(config, path, manifest_ready, exc, state_db, Path(odir) / "igvf_manifests", out, OfflineReader(), None,
                               output_dir=odir)
            if res["gaps"]:
                problems.append(f"{res['gaps']} manifest gap row(s)")
        for d, c in sorted(man - manifest_ready):
            problems.append(f"{d}/{c}: manifest-eligible but has no CellAnnotation")
        rows = coverage_report_rows(odir, config, annotations, Path(odir) / "igvf_manifests")
        write_rows_tsv(Path(odir) / "report.tsv", REPORT_HEADER, rows)
        lines += ["## [3] manifest (preview) and [4] audit", "", f"outcomes {res['outcomes'] if res else {}}", ""]
    lines += [f"- DEGRADED: {x}" for x in degraded] + [f"- PROBLEM: {x}" for x in problems]
    rc = 1 if (degraded or problems) else 0
    write_report(out, "E2G QC-and-Predictions run", lines, {"mode": mode, "dry_run": args.dry_run, "exit": rc, "degraded": degraded,
                 "problems": problems, "made": made, "manifest_ready": sorted(map(list, manifest_ready)),
                 "reformat_eligible": sorted(map(list, reformat_eligible))})
    print(f"RESULT: {'INCOMPLETE -- see PROBLEM/DEGRADED lines' if rc else 'complete'} (exit {rc})")
    return rc


# ---------------------------------------------------------------------------
# Self-test: a synthetic world with planted QC failures, annotations and files
# ---------------------------------------------------------------------------

S1, S2, S3 = "IGVFSM0001AAAA", "IGVFSM0002BBBB", "IGVFSM0003CCCC"


def _gz_write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def _metrics_text(rows) -> str:
    return "\t".join(METRIC_COLUMNS) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


def _pred_frame(model: str, n: int = 8, seed: int = 0):
    pd, np = _pd(), _np()
    rng = np.random.RandomState(seed)
    starts = np.arange(n) * 5000 + 1000
    cls = ["promoter" if i == 0 else ("genic" if i % 3 == 0 else "intergenic") for i in range(n)]
    df = pd.DataFrame({"chr": ["chr1"] * (n // 2) + ["chr2"] * (n - n // 2), "start": starts, "end": starts + 500,
                       "class": cls, "TargetGene": [f"GENE{i % 3}" for i in range(n)], "TargetGeneTSS": starts + 20000,
                       "TargetGeneEnsembl_ID": [f"ENSG0000000000{i % 3}" for i in range(n)],
                       "isSelfPromoter": [i == 0 for i in range(n)], "CellType": "k562", "distance": 20000 - np.arange(n) * 10,
                       "E2G.Score.qnorm": np.round(np.linspace(0.05, 0.95, n), 4), "E2G.Score.qnorm.ignoreTPM": np.round(np.linspace(0.06, 0.96, n), 4),
                       "normalizedATAC_prom": rng.rand(n).round(3), "numTSSEnhGene": np.arange(n), "numNearbyEnhancers": np.arange(n) + 1,
                       "ubiqExpressed": [0] * n, "numCandidateEnhGene": np.arange(n) * 2, "ABC.Score": rng.rand(n).round(4),
                       "normalizedATAC_enh": rng.rand(n).round(3)})
    df.insert(3, "name", [f"{c}|{a}:{s}-{s + 500}" for c, a, s in zip(cls, df["chr"], starts)])
    if model.startswith("multiome"):
        df["ARC.E2G.Score"] = rng.rand(n).round(4)
        df["Kendall"] = rng.rand(n).round(4)
        df["RNA_pseudobulkTPM"] = rng.rand(n).round(2) * 10
        df["RNA_meanLogNorm"] = rng.rand(n).round(3)
        df["RNA_percentCellsDetected"] = rng.rand(n).round(3)
    return df


def synthetic_world(d: Path) -> dict:
    pd, np = _pd(), _np()
    W = {"root": d}
    data = d / "data"
    plots = data / "plots" / "ds1"
    bc = {S1: [f"AAAC{i:012d}_{S1}" for i in range(3)], S2: [f"CCCG{i:012d}_{S2}" for i in range(2)]}
    guide_rows = [(b, s, "IGVFDS1105KTIQ") for s in (S1, S2) for b in bc[s]]

    def guide(cl, rows):
        return _gz_write(plots / cl / DEFAULT_QC_GUIDE_NAME, "barcode\tsubsample\tanalysis_accession\n" +
                         "".join("\t".join(r) + "\n" for r in rows))
    W["guide_k562"] = guide("k562", guide_rows)
    (plots / "k562" / METRICS_NAME).write_text(_metrics_text([[S1, 300, 1800000, 1200000, 6000, 4000, 0.4, 8.0],
                                                             [S2, 200, 1200000, 800000, 6000, 4000, 0.2, 6.0]]))
    (plots / "k562" / "qc_thresholds.tsv").write_text("metric\tmin\tmax\nnum_frags\t1000\t100000\n")
    guide("low", [(f"GGGT{i:012d}_{S1}", S1, "IGVFDS1105KTIQ") for i in range(2)])
    (plots / "low" / METRICS_NAME).write_text(_metrics_text([[S1, 50, 5000000, 3000000, 1e5, 6e4, 0.3, 7.0]]))
    (plots / "low" / "qc_thresholds.tsv").write_text("metric\tmin\nnum_frags\t1000\n")
    guide("nometrics", [(f"TTTA{i:012d}_{S1}", S1, "IGVFDS1105KTIQ") for i in range(2)])
    mcf = [(f"ACGT{i:012d}_{S3}", S3, "IGVFDS1105KTIQ") for i in range(4)]
    W["guide_mcf7"] = guide("mcf7", mcf)
    dt = data / "datatables" / "ds1_data"
    dt.mkdir(parents=True, exist_ok=True)
    (dt / "mcf7_per_cell_qc.tsv").write_text("\t".join(PER_CELL_QC_COLUMNS) + "\n" + "".join(
        f"IGVFDS1\t{b}\t{S3}\t100\t50\t1\t1\t{1000 + i}\t0.1\t0.5\t7\t0.3\n" for i, (b, _, _) in enumerate(mcf + [("ZZZZ_x", S3, "")])))
    # per-component metrics for merge-metrics
    (plots / "mcf7_1").mkdir(parents=True, exist_ok=True)
    (plots / "mcf7_2").mkdir(parents=True, exist_ok=True)
    (plots / "mcf7_1" / METRICS_NAME).write_text(_metrics_text([[S3, 30, 300, 60, 10, 2, 0.2, 5.0], [S1, 10, 100, 10, 10, 1, 0.5, 9.0]]))
    (plots / "mcf7_2" / METRICS_NAME).write_text(_metrics_text([[S3, 10, 200, 40, 20, 4, 0.6, 7.0]]))
    # pseudobulks
    pb = d / "pb" / "ds1" / "pseudobulks"
    sizes = d / "chrom.sizes"
    sizes.write_text("chr1\t248956422\nchr2\t242193529\n")
    W["sizes"] = sizes
    for s in (S1, S2):
        rows = []
        for i, b in enumerate(bc[s]):
            rows += [f"chr2\t{900 + i}\t{1000 + i}\t{b}\t1", f"chr1\t{500 - i}\t{600}\t{b}\t2"]
        rows += [f"chr1\t50\t80\tNOISE_{s}\t1", f"chr3\t10\t20\t{bc[s][0]}\t1"]
        _gz_write(pb / f"annotation-k562-{s}" / "fragments.tsv.gz", "\n".join(rows) + "\n")
    try:
        import anndata as ad  # type: ignore
        import scipy.sparse as sps  # type: ignore
        genes = ["ENSG01.1", "ENSG02.1", "ENSG03.1", "ENSG04.1"]
        for s, k in ((S1, 1), (S2, 2)):
            obs = bc[s] + [f"EXTRA{k}"]
            X = sps.csr_matrix(np.array([[k, 1, 5, 2]] * len(obs), dtype=float))
            a = ad.AnnData(X=X, obs=pd.DataFrame(index=obs), var=pd.DataFrame({"gene_symbol": ["GENEA", "GENEA", "GENEC", "GENEB"]}, index=genes))
            (pb / f"annotation-k562-{s}").mkdir(parents=True, exist_ok=True)
            a.write_h5ad(str(pb / f"annotation-k562-{s}" / "rna_counts_mtx.h5ad"))
        W["has_anndata"] = True
    except ImportError:
        W["has_anndata"] = False
    gtf = d / "genes.gtf"
    gtf.write_text("".join(f'{ch}\tx\tgene\t1\t100\t.\t+\t.\tgene_id "{g}"; gene_name "{n}";\n'
                           for ch, g, n in (("chr1", "ENSG01.1", "GENEA"), ("chr2", "ENSG02.1", "GENEA"), ("chrUn_x", "ENSG03.1", "GENEC"),
                                            ("chr1", "ENSG04.1", "GENEB"))))
    W["gtf"] = gtf
    for s in (S1, S2):
        _gz_write(pb / f"annotation-k562-{s}" / "per_cell_qc.tsv.gz", "\t".join(PER_CELL_QC_COLUMNS) + "\n" +
                  "".join(f"IGVFDS1\t{b}\t{s}\t10\t5\t1\t1\t100\t0.1\t0.5\t7\t0.3\n" for b in bc[s]))
    # scE2G models + results tree
    sce2g = d / "scE2G"
    for m, t in ((MULTIOME_MODEL, ".177"), (SCATAC_MODEL, ".174")):
        (sce2g / "models" / m).mkdir(parents=True, exist_ok=True)
        (sce2g / "models" / m / f"score_threshold_{t}").write_text("")
    res = d / "results"
    cdir = res / "uniformly_processed" / "ds1" / "k562"
    for m, t, seed in ((MULTIOME_MODEL, "0.177", 1), (SCATAC_MODEL, "0.174", 2)):
        pred = _pred_frame(m, seed=seed)
        (cdir / m).mkdir(parents=True, exist_ok=True)
        pred.to_csv(cdir / m / "scE2G_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
        thr = pred[pred["E2G.Score.qnorm"] >= float(t)]
        thr.to_csv(cdir / m / f"scE2G_predictions_threshold{t}.tsv.gz", sep="\t", index=False, compression="gzip")
        bp = "".join(f"{r['chr']}\t{r['start']}\t{r['end']}\t{r['chr']}\t{r['TargetGeneTSS']}\t{r['TargetGeneTSS']}\t{r['TargetGene']}"
                     f"\t{r['E2G.Score.qnorm']}\t.\t.\n" for _, r in thr.iloc[::-1].iterrows())
        (cdir / m / f"scE2G_predictions_threshold{t}.bedpe").write_text(bp)
        st = {"n_enh_elements": 5, "n_genes_with_enh": 3, "n_enh_gene_links": 6, "mean_genes_per_enh": 1.2, "mean_enh_per_gene": 2.0,
              "mean_dist_to_tss": 18000 + seed, "mean_enh_width": 500, "n_genes_active_promoter": 3, "n_genes_not_expressed": 0,
              "cluster": "k562", "model_name": m, "fragments_total": 3000000, "cell_count": 500 if m == MULTIOME_MODEL else 0,
              "umi_count": 2000000 if m == MULTIOME_MODEL else 0}
        old = cdir / m / "scE2G_predictions_threshold0.5_stats.tsv"
        write_rows_tsv(old, list(st), [{**st, "n_enh_elements": 999}], quiet=True)
        os.utime(old, (time.time() - 1000, time.time() - 1000))
        write_rows_tsv(cdir / m / f"scE2G_predictions_threshold{t}_stats.tsv", list(st), [st], quiet=True)
    _gz_write(cdir / MULTIOME_MODEL / "scE2G_gene_list.tsv.gz", "chr\ttss\tname\tEnsembl_ID\tstrand\tExpression\n"
              "chr1\t21000\tGENE0\tENSG00000000000\t+\t3.2\nchr1\t26000\tGENE1\tENSG00000000001\t-\t1.1\n")
    _gz_write(cdir / MULTIOME_MODEL / "scE2G_element_list.tsv.gz", "chr\tstart\tend\tclass\tname\tRPM\n"
              "chr1\t1000\t1500\tpromoter\tpromoter|chr1:1000-1500\t3\n")
    (cdir / "Neighborhoods").mkdir(parents=True, exist_ok=True)
    (cdir / "Neighborhoods" / "EnhancerList.bed").write_text("chr1\t1000\t1500\tpromoter|chr1:1000-1500\nchr1\t6000\t6500\tintergenic|chr1:6000-6500\n"
                                                          "chr2\t100\t600\tgenic|chr2:100-600\n")
    (cdir / "ATAC_norm.bw").write_bytes(b"bigwig-placeholder")
    catfrag = _gz_write(d / "catlas" / "IGVFFI3117BGUY.tsv.gz", "chr1\t1\t50\tB1\t1\nchr1\t5\t60\tB2\t1\nchr1\t9\t70\tB1\t1\n")
    lab_ann = d / "lab_annotations.tsv"
    lab_ann.write_text("dataset\tlab_celltype\tCL term\tqualifier\tCL_ID\nDS1\tk562\tK562\t\tEFO:0002067\n")
    W.update({"data": data, "plots": plots, "pb": d / "pb", "sce2g": sce2g, "results": res, "cdir": cdir, "catfrag": catfrag,
              "lab_ann": lab_ann, "bc": bc})
    config = {
        "scE2G_dir": str(sce2g), "data_dir": str(data), "output_dir": str(res), "pseudobulks_root": str(d / "pb"), "scE2G_version": "1.2",
        "pipeline_mode": "default", "igvf": {"principal_analysis_sets": {"ds1": "IGVFDS1105KTIQ"}},
        "exclusion": {"user_specified": [], PREDICT_EVERYTHING_KEY: False,
                      "auto_thresholds": {"min_cell_count": 100, "min_fragments_total": 2000000, "min_umi_count": 1000000}},
        "clusters": {"ds1": {
            "k562": {"pseudobulk_annotation": "k562", "qc_guide": str(W["guide_k562"]), "models": [MULTIOME_MODEL, SCATAC_MODEL]},
            "low": {"pseudobulk_annotation": "low", "qc_guide": str(plots / "low" / DEFAULT_QC_GUIDE_NAME), "models": [MULTIOME_MODEL]},
            "nometrics": {"pseudobulk_annotation": "nometrics", "qc_guide": str(plots / "nometrics" / DEFAULT_QC_GUIDE_NAME), "models": [MULTIOME_MODEL]},
            "k562_ATAC_only": {"pseudobulk_annotation": "k562", "qc_guide": str(W["guide_k562"]), "models": [SCATAC_MODEL],
                               "cell_annotation_key": "k562", "igvf_manifest_excluded": True},
            "mcf7": {"pseudobulk_annotation": "mcf7_1,mcf7_2", "qc_guide": str(W["guide_mcf7"]), "models": [MULTIOME_MODEL], "prefiltered": True}},
            "catlas": {"EMSN_234": {"atac_frag_file": str(catfrag), "models": [SCATAC_MODEL]}}}}
    W["config_path"] = d / "ds1_pipeline_config.json"
    W["config_path"].write_text(json.dumps(config, indent=1))
    W["config"] = config

    def ct(term_id, name):
        return {"@id": f"/sample-terms/{term_id.replace(':', '_')}/", "term_id": term_id, "term_name": name}
    ana = [{"@id": "/analysis-sets/IGVFDS1105KTIQ/", "accession": "IGVFDS1105KTIQ", "file_set_type": "principal analysis"}]
    k562_ct = ct("EFO:0002067", "K562")
    mr = [
        {"@id": "/pseudobulk-sets/IGVFDS0001PSBK/", "aliases": [f"anshul-kundaje:ds1-k562-{S1}"], "cell_annotation": "K562", "cell_qualifier": None,
         "cell_type": k562_ct, "input_file_sets": ana, "samples": [{"accession": S1}], "status": "released", "lab": {"@id": DEFAULT_FILE_LAB},
         "files": [{"accession": "IGVFFI0001FRAG", "content_type": "fragments", "file_format": "tsv", "status": "in progress",
                    "submitted_file_name": f"/oak/x/igvf/ds1/pseudobulks/annotation-k562-{S1}/fragments.tsv.gz", "href": "/f/1", "md5sum": "x",
                    "upload_status": "validated"},
                   {"accession": "IGVFFI0001PCQC", "content_type": "per-cell quality report", "file_format": "tsv", "status": "in progress",
                    "submitted_file_name": f"/oak/x/igvf/ds1/pseudobulks/annotation-k562-{S1}/per_cell_qc.tsv.gz", "href": "/f/2",
                    "upload_status": "validated"}]},
        {"@id": "/pseudobulk-sets/IGVFDS0002PSBK/", "aliases": [f"anshul-kundaje:ds1-k562-{S2}"], "cell_annotation": "K562 derived from HEK293",
         "cell_qualifier": "contaminated", "cell_type": k562_ct, "input_file_sets": ana, "samples": [{"accession": S2}], "status": "in progress"},
        {"@id": "/pseudobulk-sets/IGVFDS0003PSBK/", "aliases": [f"anshul-kundaje:ds1-mcf7_1-{S3}"], "cell_annotation": "MCF-7 MED13L+",
         "cell_qualifier": "MED13L+", "cell_type": ct("EFO:0001203", "MCF-7"), "input_file_sets": ana, "samples": [{"accession": S3}], "status": "released"},
        {"@id": "/pseudobulk-sets/IGVFDS0004PSBK/", "aliases": [f"anshul-kundaje:ds1-mcf7_2-{S3}"], "cell_annotation": "MCF-7 VMP1+",
         "cell_qualifier": "VMP1+", "cell_type": ct("EFO:0001203", "MCF-7"), "input_file_sets": ana, "samples": [{"accession": S3}], "status": "released"},
        {"@id": "/pseudobulk-sets/IGVFDS0005PSBK/", "aliases": [f"anshul-kundaje:ds1-other-{S1}"], "cell_annotation": "HepG2", "cell_type": ct("EFO:0001187", "HepG2"),
         "input_file_sets": ana, "samples": [{"accession": S1}, {"accession": S2}], "status": "released"},
        {"@id": "/pseudobulk-sets/IGVFDS0006PRIN/", "aliases": ["jesse-engreitz:ds1_k562_filtered_pseudobulk_set"], "cell_annotation": "K562",
         "input_file_sets": [{"@id": "/pseudobulk-sets/IGVFDS0001PSBK/"}], "samples": [{"accession": S1}], "status": "in progress"},
        {"@id": "/pseudobulk-sets/IGVFDS3920QDLI/", "aliases": ["yang-li:catlas-human-pseudobulk-EMSN_234"], "lab": {"@id": WASHU_LAB},
         "cell_annotation": "brain EMSN_234 medium spiny neuron (Hsap)", "cell_type": ct("PCL:0051235", "EMSN_234 medium spiny neuron (Hsap)"),
         "input_file_sets": [{"@id": "/curated-sets/IGVFDS9999CURA/"}], "status": "in progress",
         "files": [{"content_type": "fragments", "aliases": ["yang-li:catlas-human-EMSN_234_fragments"], "status": "in progress",
                    "href": "/tabular-files/IGVFFI3117BGUY/@@download/IGVFFI3117BGUY.tsv.gz", "md5sum": md5_file(catfrag)}]},
        {"@id": "/pseudobulk-sets/IGVFDS1981DGID/", "aliases": ["anshul-kundaje:pseudobulk-IGVFDS8428QVAO-endothelial cell-IGVFSM0220TLGT"],
         "lab": {"@id": DEFAULT_FILE_LAB}, "samples": [{"accession": "IGVFSM0220TLGT"}], "status": "in progress",
         "input_file_sets": [{"@id": "/analysis-sets/IGVFDS8428QVAO/", "accession": "IGVFDS8428QVAO", "file_set_type": "principal analysis"}],
         "cell_type": ct("CL:0000115", "endothelial cell"),
         "files": [{"accession": "IGVFFI0694XKXE", "content_type": "cell by gene matrix", "file_format": "h5ad", "status": "in progress",
                    "submitted_file_name": "pseudobulks/annotation_0-IGVFSM0220TLGT/rna_counts_mtx.h5ad", "href": "/f/3"}]},
    ]
    W["multireport"] = d / "multireport.json"
    W["multireport"].write_text(json.dumps({"@graph": mr}))
    return W


class _FakePortal:
    """In-memory stand-in for igvf_submission_skill.Portal used ONLY by the selftest (no network)."""
    base = "memory://portal"

    def __init__(self):
        self.objs, self.calls = {}, []

    def _request(self, method, path, params=None, payload=None, anon=False):
        from igvf_submission_skill import PortalError  # type: ignore
        self.calls.append((method, path, payload))
        key = path.strip("/")
        if method == "POST":
            obj = dict(payload, uuid=f"uuid-{len(self.calls):04d}", **{"@id": f"{path}{len(self.calls)}/"})
            for a in payload.get("aliases", []):
                self.objs[a] = obj
            self.objs[obj["uuid"]] = obj
            return {"@graph": [obj]}
        if method == "PATCH":
            if key not in self.objs:
                raise PortalError(f"HTTP 404 on PATCH {path}")
            self.objs[key].update(payload)
            return {"@graph": [self.objs[key]]}
        if key in self.objs:
            return self.objs[key]
        raise PortalError(f"HTTP 404 on GET {path}")


_IN_SELFTEST = [False]


def cmd_selftest(args) -> int:
    _IN_SELFTEST[0] = True
    try:
        return _selftest(args)
    finally:
        _IN_SELFTEST[0] = False


def _selftest(args) -> int:
    _pd()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    run = main

    t0 = time.time()
    before = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_world(d)
        cfgp, config = str(W["config_path"]), W["config"]
        o = lambda name: str(d / "out" / name)  # noqa: E731
        sdb = str(d / "state.db")

        print("\nconfig parsing")
        y = parse_yaml_subset("a: 1\nb:\n  c: [x, y]\n  d: true\n  e:\n    - p\n    - q\n# c\nf: 'str # not comment'\ng: 2e6\n")
        check(y == {"a": 1, "b": {"c": ["x", "y"], "d": True, "e": ["p", "q"]}, "f": "str # not comment", "g": 2e6},
              "YAML subset: nested maps, flow and block lists, booleans, quoted '#', floats")
        check(config_bool("false") is False and config_bool(0) is False and config_bool("on") is True, "config_bool: 'false' string is False")

        print("\nquality gate")
        inc, up, exc, st = resolve_exclusions(config)
        rs = {k[1]: v["reason"] for k, v in st.items()}
        check(rs == {"k562": "pass", "low": "below_min_cell_count", "nometrics": "missing_metrics", "k562_ATAC_only": "pass",
                     "mcf7": "below_min_cell_count", "EMSN_234": "missing_metrics"}, f"reasons {rs}")
        check(st[("ds1", "mcf7")]["cell_count"] == 4 and st[("ds1", "mcf7")]["fragments_total"] == 1000 + 1001 + 1002 + 1003,
              "prefiltered: QC-guide barcodes joined to the per-cell QC table (extra table row ignored)")
        check(st[("ds1", "k562_ATAC_only")]["umi_count"] is None, "ATAC-only cluster: umi_count None, UMI threshold skipped")
        check(manifest_eligible_clusters(config, up) == {("ds1", "k562")}, "igvf_manifest_excluded keeps the ATAC-only variant out of manifests")
        c2 = json.loads(json.dumps(config))
        c2["exclusion"][PREDICT_EVERYTHING_KEY] = True
        inc2, up2, _, _ = resolve_exclusions(c2)
        check(("ds1", "low") in inc2 and ("ds1", "nometrics") not in inc2 and ("ds1", "low") not in up2,
              "predictions_on_everything: quality failures re-included, missing inputs never, upload still gated")
        _, _, exc3, st3 = resolve_exclusions(config, legacy_implicit=True)
        check(("ds1", "nometrics") not in exc3 and st3[("ds1", "k562")]["reason"] == "pass",
              "--legacy-implicit-source: 'no stats yet' is not an exclusion (synapse-submission / CATlas rule)")
        rc = run(["resolve-exclusions", "--config", cfgp, "--out-dir", o("excl")])
        plan = read_rows_tsv(Path(o("excl")) / "igvf_metadata" / "ds1_pipeline_plan.tsv")
        check(rc == 0 and [r["cluster"] for r in plan if r["manifest_eligible"] == "y"] == ["k562"], "resolve-exclusions: plan TSV")

        print("\nmetrics tables")
        rc = run(["merge-metrics", "--metrics", str(W["plots"] / "mcf7_1" / METRICS_NAME), str(W["plots"] / "mcf7_2" / METRICS_NAME),
                  "--out", o("merged_metrics.tsv")])
        mm = {r["subsample"]: r for r in read_rows_tsv(o("merged_metrics.tsv"))}
        check(rc == 0 and mm[S3]["n_cells"] == "40" and mm[S3]["total_fragments"] == "500" and float(mm[S3]["mean_frag_per_cell"]) == 12.5
              and abs(float(mm[S3]["mean_frip"]) - (0.2 * 30 + 0.6 * 10) / 40) < 1e-12, "merge-metrics: sums, ratio recomputed, n_cells-weighted frip")
        rc = run(["prefiltered-metrics", "--frag-file", str(W["catfrag"]), "--subsample-name", "EMSN_234",
                  "--out", str(W["data"] / "plots" / "catlas" / "EMSN_234" / METRICS_NAME)])
        pm = read_rows_tsv(W["data"] / "plots" / "catlas" / "EMSN_234" / METRICS_NAME)[0]
        check(rc == 0 and pm["n_cells"] == "2" and pm["total_fragments"] == "3" and pm["mean_frip"] == "NA", "prefiltered-metrics (CATlas)")
        _, _, _, st4 = resolve_exclusions(config)
        check(st4[("catlas", "EMSN_234")]["cell_count"] == 2 and st4[("catlas", "EMSN_234")]["umi_count"] is None,
              "atac_frag_file cluster reads the prefiltered metrics file")
        rc = run(["build-qc-datatables", "--pseudobulks-root", str(W["pb"]), "--dest", o("dt"), "--config", cfgp])
        dtab = read_rows_tsv(Path(o("dt")) / "QC_datatables" / "ds1_data" / "k562_per_cell_qc.tsv")
        check(rc == 0 and len(dtab) == 5 and [r["subsample"] for r in dtab] == [S1] * 3 + [S2] * 2, "build-qc-datatables: sorted concat, one header")

        print("\nfiltering")
        md = W["results"] / "multiome_data" / "ds1" / "k562"
        frag_out = md / "atac_fragments_ds1_k562.tsv.gz"
        rc = run(["filter-atac", "--qc-guide", str(W["guide_k562"]), "--pseudobulks", str(W["pb"] / "ds1" / "pseudobulks"),
                  "--cell-type", "k562", "--chrom-sizes", str(W["sizes"]), "--out", str(frag_out)])
        lines = [ln.split("\t") for ln in gzip.open(frag_out, "rt").read().splitlines()]
        check(rc == 0 and len(lines) == 10 and all(r[3] in {b for v in W["bc"].values() for b in v} for r in lines),
              "filter-atac: only guide barcodes (full strings), noise dropped")
        check([r[0] for r in lines] == ["chr1"] * 5 + ["chr2"] * 5 and [int(r[1]) for r in lines[:5]] == sorted(int(r[1]) for r in lines[:5])
              and os.path.exists(str(frag_out) + ".tbi"), "filter-atac: chrom-sizes order, start-sorted, chr3 dropped, tabix index")
        try:
            filter_atac(W["guide_mcf7"], W["pb"] / "ds1" / "pseudobulks", "k562", W["sizes"], d / "x.tsv.gz")
            check(False, "filter-atac: guide barcodes absent from fragments must fail")
        except SystemExit:
            check(True, "filter-atac: guide barcodes absent from every fragment file -> error (legacy branches: warning)")
        if W["has_anndata"]:
            rc = run(["filter-rna", "--qc-guide", str(W["guide_k562"]), "--pseudobulks", str(W["pb"] / "ds1" / "pseudobulks"),
                      "--cell-type", "k562", "--out", str(md / "rna_count_matrix_ds1_k562.mtx"), "--gtf", str(W["gtf"]),
                      "--standard-chromosomes-only", "--log", o("gtf_log.txt")])
            import scipy.io  # type: ignore
            mdir = md / "rna_count_matrix_ds1_k562"
            M = scipy.io.mmread(gzip.open(mdir / "matrix.mtx.gz")).toarray()
            feats = gzip.open(mdir / "features.tsv.gz", "rt").read().split()
            cells = gzip.open(mdir / "barcodes.tsv.gz", "rt").read().split()
            ga = feats.index("GENEA")
            check(rc == 0 and feats == ["GENEA", "GENEB"] and len(cells) == 5 and M.shape == (2, 5)
                  and list(M[ga]) == [2, 2, 2, 3, 3], "filter-rna: GTF collapse sums Ensembl IDs per symbol, non-standard chrom dropped, genes x cells")
            rc = run(["package-rna", "--matrix-dir", str(mdir), "--out", str(md / "rna_count_matrix_ds1_k562.tar.gz")])
            names = tarfile.open(md / "rna_count_matrix_ds1_k562.tar.gz").getnames()
            check(rc == 0 and sorted(names) == ["barcodes.tsv", "features.tsv", "matrix.mtx"], "package-rna: flat tar of decompressed files")
        else:
            check(True, "filter-rna skipped (anndata not installed)")

        print("\nscE2G configuration")
        rc = run(["sce2g-config", "--config", cfgp, "--out-dir", o("sce2g_cfg"), "--lab-annotations", str(W["lab_ann"])])
        cc = {r["cluster"]: r for r in read_rows_tsv(Path(o("sce2g_cfg")) / "ds1_cell_clusters.tsv")}
        cm = {r["cluster"]: r for r in read_rows_tsv(Path(o("sce2g_cfg")) / "ds1_cluster_metadata.tsv")}
        check(rc == 0 and cc["k562_ATAC_only"]["rna_matrix_file"] == "" and cc["k562"]["model_dir"].count("models/") == 2
              and cm["k562"]["cell_type"] == "K562" and cm["k562"]["ontology_id"] == "EFO:0002067" and cm["k562"]["summary"] == "k562",
              "sce2g-config: cell_clusters (ATAC-only has no RNA matrix) and cluster metadata joined from lab annotations")
        check(model_threshold(str(W["sce2g"]), MULTIOME_MODEL) == "0.177" and resolve_primary_model([SCATAC_MODEL]) == SCATAC_MODEL,
              "score threshold from models/<m>/score_threshold_.177 (leading zero restored); primary model rule")
        check(determine_mem_mb(1000, True) == 16000 and determine_mem_mb(1e6, False, 3) == MAX_MEM_MB, "determine_mem_mb (ABC calculator)")

        print("\nreformatting")
        mp = W["cdir"] / MULTIOME_MODEL
        rc = run(["reformat", "--input", str(mp / "scE2G_predictions.tsv.gz"), "--out", o("full.e2g.tsv.gz"), "--model", MULTIOME_MODEL,
                  "--version", "1.2", "--cell-type", "K562", "--term-id", "EFO:0002067", "--summary", "K562",
                  "--portal-link", "jesse-engreitz:ds1_k562_scE2G_Multiome_predictions_full"])
        txt = gzip.open(o("full.e2g.tsv.gz"), "rt").read().splitlines()
        hdr = [ln for ln in txt if ln.startswith("#")]
        cols = txt[len(hdr)].split("\t")
        check(rc == 0 and hdr[0] == "# Source: scE2G_multiome_powerlaw_v3" and "# ScoreType: positive_score" in hdr
              and hdr[-1] == "# Metadata: https://data.igvf.org/tabular-files/jesse-engreitz:ds1_k562_scE2G_Multiome_predictions_full"
              and not any("ScoreThreshold" in h for h in hdr), "reformat full: header block")
        check(cols == [c for c, _ in FULL_MULTIOME_COLS] and txt[len(hdr) + 1].split("\t")[3].startswith("promoter|"),
              "reformat full (Multiome): consortium column order, feature. prefix, raw ElementName")
        check(txt[len(hdr) + 1].split("\t")[cols.index("isSelfPromoter")] == "TRUE" and txt[len(hdr) + 1].split("\t")[8] == "K562",
              "reformat: logicals written TRUE/FALSE, CellAnnotation column set")
        run(["reformat", "--input", str(mp / "scE2G_predictions_threshold0.177.tsv.gz"), "--out", o("thr.e2g.tsv.gz"), "--model", MULTIOME_MODEL,
             "--version", "1.2", "--cell-type", "K562", "--term-id", "EFO:0002067", "--summary", "K562", "--threshold", "0.177"])
        t2 = gzip.open(o("thr.e2g.tsv.gz"), "rt").read().splitlines()
        h2 = [ln for ln in t2 if ln.startswith("#")]
        check("# ScoreThreshold: Score >= 0.177" in h2 and t2[len(h2)].split("\t") == [c for c, _ in THRESHOLDED_COLS]
              and re.match(r"^chr\d:\d+-\d+$", t2[len(h2) + 1].split("\t")[3]), "reformat thresholded: narrow columns, ElementName chr:start-end")
        run(["reformat", "--input", str(W["cdir"] / SCATAC_MODEL / "scE2G_predictions.tsv.gz"), "--out", o("atac_full.e2g.tsv.gz"),
             "--model", SCATAC_MODEL, "--version", "1.2", "--summary", "K562", "--catlas"])
        t3 = gzip.open(o("atac_full.e2g.tsv.gz"), "rt").read().splitlines()
        c3 = [ln for ln in t3 if not ln.startswith("#")][0].split("\t")
        check(c3[16] == "feature.ABC.Score" and "normalizedATAC_enh" not in c3 and "Kendall" not in c3,
              "reformat scATAC full: feature.ABC.Score closes the feature block; --catlas drops normalizedATAC_enh")
        run(["reformat", "--input", str(mp / "scE2G_gene_list.tsv.gz"), "--out", o("genes.tsv.gz"), "--model", MULTIOME_MODEL, "--version", "1.2"])
        g = [ln for ln in gzip.open(o("genes.tsv.gz"), "rt").read().splitlines()]
        check(not any("ScoreType" in ln for ln in g) and g[-2].split("\t") == ["chr1", "20750", "21250", "GENE0", "ENSG00000000000", "+"],
              "reformat gene list: TSS +/- 250, no ScoreType for _list files")
        run(["reformat", "--input", str(mp / "scE2G_predictions_threshold0.177.tsv.gz"), "--out", o("legacy.tsv.gz"), "--model", MULTIOME_MODEL,
             "--version", "1.2", "--summary", "k562", "--threshold", "0.177", "--format", "synapse-legacy"])
        lg = gzip.open(o("legacy.tsv.gz"), "rt").read().splitlines()
        lcols = [ln for ln in lg if not ln.startswith("#")][0].split("\t")
        check(lg[0].startswith("# PRELIMINARY") and "SampleSummaryShort" in lcols and lcols[10] == "Score.ignoreTPM",
              "--format synapse-legacy: PRELIMINARY banner, SampleSummaryShort, Score.ignoreTPM")
        rc = run(["reformat", "--kind", "element-bed", "--input", str(W["cdir"] / "Neighborhoods" / "EnhancerList.bed"), "--out", o("el.bed.gz"),
                  "--version", "1.2", "--cell-type", "K562", "--term-id", "EFO:0002067", "--summary", "K562", "--portal-link", "x:y"])
        eb = gzip.open(o("el.bed.gz"), "rt").read().splitlines()
        check(rc == 0 and eb[0] == "# Source: ABC (originally EnhancerList.bed)" and eb[10] == "#ElementChr\tElementStart\tElementEnd\tElementName"
              and len(eb) == 14 and os.path.exists(o("el.bed.gz") + ".tbi"), "reformat element-bed: header + EnhancerList verbatim, tabix")
        run(["reformat", "--kind", "bedpe", "--input", str(mp / "scE2G_predictions_threshold0.177.bedpe"), "--out", o("p.bedpe.gz")])
        bp = [ln.split("\t") for ln in gzip.open(o("p.bedpe.gz"), "rt").read().splitlines()]
        check([(r[0], int(r[1])) for r in bp] == sorted((r[0], int(r[1])) for r in bp) and os.path.exists(o("p.bedpe.gz") + ".tbi"),
              "bedpe: sort -k1,1 -k2,2n then bgzip + tabix")
        run(["candidates", "--input", str(mp / "scE2G_predictions.tsv.gz"), "--out", o("cand.tsv.gz"), "--summary", "k562"])
        cd = gzip.open(o("cand.tsv.gz"), "rt").read().splitlines()
        check(cd[0].split("\t") == [c for c, _ in CANDIDATE_COLS] and cd[1].split("\t")[8] == "k562", "candidates: columns + SampleSummaryShort")
        run(["features", "--input", str(mp / "scE2G_predictions.tsv.gz"), "--out", o("feat.tsv.gz"), "--model", MULTIOME_MODEL,
             "--cell-type", "K562", "--term-id", "EFO:0002067", "--summary", "k562"])
        ft = gzip.open(o("feat.tsv.gz"), "rt").read().splitlines()
        fcol = ft[8].split("\t")
        check(ft[0] == f"# Source: scE2G {MULTIOME_MODEL}" and fcol[:11] == [c for c, _ in FEATURE_LEAD] and "chr" in fcol[11:]
              and "E2G.Score.qnorm" in fcol, "features: fixed lead columns then everything() (chr survives)")

        print("\nQC aggregation")
        rc = run(["aggregate-qc", "--dataset-dir", str(W["results"] / "uniformly_processed" / "ds1"), "--plots-dir", str(W["plots"]),
                  "--out-dir", o("agg"), "--no-plots"])
        agg = read_rows_tsv(Path(o("agg")) / "all_qc_stats.tsv")
        check(rc == 0 and len(agg) == 2 and all(r["n_enh_elements"] == "5" for r in agg), "aggregate-qc: newest stats.tsv per (cluster, model)")
        check({r["model_name"]: r["cell_count"] for r in agg}[SCATAC_MODEL] == "500", "aggregate-qc: cell_count 0 patched from metrics")
        if not args.no_plots:
            figs, _ = qc_plots(Path(o("agg")) / "all_qc_stats.tsv", Path(o("agg")))
            check(len(figs) >= 5, f"QC figures via sce2g-pipeline ({len(figs)})")

        print("\nCell Annotation cache")
        rc = run(["cell-metadata", "--config", cfgp, "--multireport-json", str(W["multireport"]), "--state-db", sdb, "--out-dir", o("cm")])
        stt = {r["cluster"]: r for r in read_rows_tsv(Path(o("cm")) / "igvf_metadata" / "ds1_cell_annotation_status.tsv")}
        S = State(sdb)
        k = S.get_cell_annotation("ds1", "k562")
        m7 = S.get_cell_annotation("ds1", "mcf7")
        prim = S.q("SELECT alias FROM cell_metadata_primary_pseudobulks")
        S.close()
        check(rc == 0 and k["cell_annotation"] == "K562" and k["cell_qualifier"] is None and k["cl_id"] == "EFO_0002067"
              and k["term_id"] == "EFO:0002067", "derive: most-contributing subsample wins the annotation disagreement")
        check(k["principal_uploaded"] == 1 and k["all_primary_released"] == 0, "locked by an uploaded principal (not all primaries released)")
        check(m7["cell_annotation"] == "MCF-7" and m7["cell_qualifier"] is None, "merged cluster: constituent aliases, bare term_name, no qualifier")
        check(stt["low"]["reason"].startswith("no_matching_primary_alias") and len(prim) == 5,
              "unmatched subsample reported; multi-sample primary and non-primary rows not cached")
        snap = snapshot_path(o("cm"), "ds1")
        digest = cluster_set_digest({(d_, c_) for d_, c_, _ in iter_clusters(config)})
        ok_snap = len(read_snapshot(snap, digest)) >= 2
        try:
            read_snapshot(snap, "0" * 16)
            bad = False
        except SnapshotError:
            bad = True
        check(ok_snap and bad, "snapshot: digest + age enforced on read, foreign cluster set rejected")
        rc = run(["cell-metadata", "--catlas-seed", "--multireport-json", str(W["multireport"]), "--state-db", sdb, "--out-dir", o("seed")])
        S = State(sdb)
        cat = S.get_cell_annotation("catlas", "EMSN_234")
        S.close()
        check(rc == 0 and cat and cat["term_id"] == "PCL:0051235", "CATlas seed: cluster = last alias token, real term id")
        rc = run(["cell-annotation-report", "--multireport-json", str(W["multireport"]), "--out-dir", o("car1")])
        car = {r["Dataset_Cluster"]: r for r in read_rows_tsv(Path(o("car1")) / "cell_annotations_by_dataset_cluster.tsv")}
        check(rc == 0 and car["ds1-k562"]["CellAnnotation"] == "K562 | K562 derived from HEK293" and car["ds1-k562"]["QCGuideFile"] == ""
              and car["ds1-k562"]["Status"] == "in progress", "cell-annotation-report mode 1: distinct values ' | '-joined")
        run(["cell-annotation-report", "--multireport-json", str(W["multireport"]), "--qc-guide-dir", str(W["plots"].parent), "--out-dir", o("car2")])
        car2 = {r["Dataset_Cluster"]: r for r in read_rows_tsv(Path(o("car2")) / "cell_annotations_by_dataset_cluster.tsv")}
        check(car2["ds1-k562"]["CellAnnotation"] == "K562" and car2["ds1-k562"]["Subsamples"] == f"{S1} | {S2}"
              and car2["ds1-k562"]["QCGuideFile"] == DEFAULT_QC_GUIDE_NAME, "cell-annotation-report mode 2: resolved via QC guide")
        run(["dataset-accessions", "--multireport-json", str(W["multireport"]), "--json-out", o("acc.json")])
        acc = json.loads(Path(o("acc.json")).read_text())
        check(acc.get("ds1") == ["IGVFDS1105KTIQ"], f"dataset-accessions: {acc}")
        run(["washu-report", "--multireport-json", str(W["multireport"]), "--out-dir", o("washu")])
        wr = read_rows_tsv(Path(o("washu")) / "washu_pseudobulk_report.tsv")
        wr = [r for r in wr if r["Alias"].startswith("yang-li")]
        (d / "catlas_dl").mkdir()
        shutil.copy(W["catfrag"], d / "catlas_dl" / "IGVFFI3117BGUY.tsv.gz")
        write_rows_tsv(d / "washu.tsv", WASHU_HEADER, wr, quiet=True)
        rc = run(["verify-fragments", "--report", str(d / "washu.tsv"), "--outdir", str(d / "catlas_dl")])
        check(rc == 0 and wr[0]["SampleTermID"] == "PCL:0051235", "washu-report + verify-fragments (md5 PASS_CACHED, no download)")

        print("\nportal file discovery")
        run(["portal-files", "--multireport-json", str(W["multireport"]), "--download-root", str(d / "dl"), "--out-dir", o("pf")])
        pf = read_rows_tsv(Path(o("pf")) / "portal_files.tsv")
        irr = next(r for r in pf if r["accession"] == "IGVFFI0694XKXE")
        check(len(pf) == 3 and irr["annotation"] == "endothelial cell" and irr["dataset"] == "IGVFDS8428QVAO"
              and "annotation_0-IGVFSM0220TLGT" in irr["rel_path"] and "set_missing" in irr["review_reasons"],
              "portal-files: annotation_0 placeholder recovered from the alias, dataset from principal analysis set")
        loc = next(r for r in pf if r["accession"] == "IGVFFI0001PCQC")
        Path(loc["local_path"]).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(W["pb"] / "ds1" / "pseudobulks" / f"annotation-k562-{S1}" / "per_cell_qc.tsv.gz", loc["local_path"])
        arch = d / "archive" / "ds1" / "pseudobulks" / f"annotation-k562-{S1}"
        arch.mkdir(parents=True)
        (arch / "per_cell_qc.tsv").write_bytes(gzip.open(loc["local_path"]).read())
        rc = run(["compare-archive", "--discovery", str(Path(o("pf")) / "portal_files.tsv"), "--archive-root", str(d / "archive"), "--out-dir", o("cmp")])
        cmpr = {r["accession"]: r for r in read_rows_tsv(Path(o("cmp")) / "portal_vs_archive.tsv")}
        check(rc == 0 and cmpr["IGVFFI0001PCQC"]["verdict"] == COMPRESSION_ONLY and cmpr["IGVFFI0001FRAG"]["status"] == "not_downloaded",
              "compare-archive: .tsv.gz vs .tsv of identical content = compression_only")

        print("\ndriver (run)")
        rc = run(["run", "--config", cfgp, "--multireport-json", str(W["multireport"]), "--state-db", sdb, "--out-dir", o("run"),
                  "--lab-annotations", str(W["lab_ann"])])
        cdir = W["cdir"]
        made = [f"ds1_k562_scE2G_{MULTIOME_MODEL}.e2g.tsv.gz", f"ds1_k562_scE2G_{MULTIOME_MODEL}_threshold0.177.e2g.tsv.gz",
                f"ds1_k562_scE2G_{SCATAC_MODEL}_threshold0.174.e2g.tsv.gz", "ds1_k562_scE2G_multiome_v3_gene_list.tsv.gz",
                "ds1_k562_element_list.bed.gz", "ds1_k562_element_list.bed.gz.tbi", f"ds1_k562_scE2G_{MULTIOME_MODEL}_threshold0.177.bedpe.gz.tbi"]
        check(all((cdir / m).exists() for m in made), "run: reformat, element BED, gene list, bedpe.gz + tabix for the annotated cluster")
        check((W["results"] / "uniformly_processed" / "candidate_e2g_pairs" / "ds1_k562_candidate_e2g_pairs.tsv.gz").exists()
              and (W["results"] / "uniformly_processed" / "scE2G_features" / "ds1" / "k562" / f"ds1_k562_scE2G_{MULTIOME_MODEL}_features.tsv.gz").exists()
              and (W["results"] / "uniformly_processed" / "ds1" / "qc_plots" / "all_qc_stats.tsv").exists(), "run: candidates, features, all_qc_stats")
        rpt = {r["cluster"]: r for r in read_rows_tsv(W["results"] / "report.tsv")}
        cov = read_rows_tsv(W["results"] / "igvf_manifests" / "ds1" / "manifest_coverage.tsv")
        check(rc == 0 and rpt["k562"]["status"] == "fully-processed" and rpt["k562"]["manifest_status"] == "ready"
              and rpt["low"]["status"] == "excluded-quality" and rpt["nometrics"]["status"] == "missing-qc-guide"
              and rpt["k562_ATAC_only"]["manifest_status"] == "excluded", "run: exit 0, report.tsv statuses and manifest roll-up")
        check(not snapshot_path(W["results"], "ds1").exists() and len(cov) == 15 and all(r["outcome"] == "planned-post" for r in cov),
              "run: temp snapshot removed; 15 manifest rows planned for the one manifest-ready cluster")
        hdr_k = embedded_header(str(cdir / f"ds1_k562_scE2G_{MULTIOME_MODEL}.e2g.tsv.gz"))
        check(hdr_k["# SampleTermID:"] == "EFO:0002067" and hdr_k["# CellAnnotation:"] == "K562", "reformatted headers carry the portal annotation")

        print("\nIGVF manifests")
        mdir = d / "manifests1"
        live = {"jesse-engreitz:ds1_k562_QC_thresholds": {"uuid": "live-qc-uuid"}}
        (d / "live.json").write_text(json.dumps(live))
        fx = {f"/{r}/": {"@id": f"/x/{r}/", "accession": r, "status": "released"} for r in
              ["IGVFFI7969JLFC", "IGVFFI0653VCGH", "IGVFFI9573KOZR", "IGVFFI6788CPPS", f"anshul-kundaje:ds1-k562-{S1}"]}
        fx["/IGVFDS1105KTIQ/"] = {"@id": "/analysis-sets/IGVFDS1105KTIQ/", "@type": ["AnalysisSet", "FileSet", "Item"], "accession": "IGVFDS1105KTIQ",
                                  "status": "released", "file_set_type": "principal analysis", "input_for": [{"@id": "/pseudobulk-sets/IGVFDS0001PSBK/"}]}
        fx["/pseudobulk-sets/IGVFDS0001PSBK/"] = {"@id": "/pseudobulk-sets/IGVFDS0001PSBK/", "@type": ["PseudobulkSet", "FileSet", "Item"],
                                                  "accession": "IGVFDS0001PSBK", "status": "released"}
        (d / "lineage.json").write_text(json.dumps(fx))
        rc = run(["manifest", "--config", cfgp, "--cluster-keys", "ds1/k562,ds1/low", "--state-db", sdb, "--manifest-dir", str(mdir),
                  "--out-dir", o("man1"), "--live-aliases", str(d / "live.json"), "--lineage-fixture", str(d / "lineage.json"), "--lineage-walk", "--offline"])
        files = sorted(p.name for p in (mdir / "ds1").glob("round*.tsv"))
        expect = ["round1_QC_documents_patch.tsv", "round1_QC_documents_post.tsv", "round2_principal_pseudobulk_set_post.tsv", "round3_filtered_barcode_list_post.tsv",
                  "round3_prediction_set_post.tsv", "round4_filtered_atac_fragment_file_post.tsv", "round4_filtered_rna_count_matrix_post.tsv",
                  "round5_atac_index_file_post.tsv", "round5_prediction_tabular_files_elements_bed_post.tsv",
                  "round5_prediction_tabular_files_genes_post.tsv", "round5_signal_files_atac_bw_post.tsv",
                  "round6_elements_bed_index_file_elements_bed_index_post.tsv", "round6_prediction_tabular_files_full_post.tsv",
                  "round7_prediction_tabular_files_thresholded_post.tsv", "round8_prediction_tabular_files_bedpe_post.tsv",
                  "round9_bedpe_index_file_bedpe_index_post.tsv"]
        check(files == expect, f"fifteen round files, dependency-ordered rounds, live alias -> PATCH ({len(files)})")
        cov1 = read_rows_tsv(mdir / "ds1" / "manifest_coverage.tsv")
        low = Counter(r["outcome"] for r in cov1 if r["cluster"] == "low")
        check(rc == 1 and low["invalid"] >= 2 and low["skipped-missing-file"] >= 5 and low["planned-post"] >= 2,
              f"cluster without annotation/files: invalid + skipped-missing-file gaps reported, exit 1 ({dict(low)})")
        qc = read_iu_tsv(mdir / "ds1" / "round1_QC_documents_patch.tsv")["live-qc-uuid"]
        check(qc["attachment"] == json.dumps({"path": str(W["plots"] / "k562" / "qc_thresholds.tsv")}) and qc["document_type"] == "pipeline parameters",
              "iu_register TSV: record_id first, object JSON-dumped without csv quoting")
        fr = read_rows_tsv(mdir / "ds1" / "round6_prediction_tabular_files_full_post.tsv")[0]
        check(fr["derived_from"] == ("jesse-engreitz:scE2G_Multiome_trained_model,jesse-engreitz:ds1_k562_filtered_ATAC_fragment_file,"
                                     "jesse-engreitz:ds1_k562_scE2G_Multiome_predictions_elements_bed,jesse-engreitz:ds1_k562_filtered_RNA_count_matrix,"
                                     "jesse-engreitz:ds1_k562_scE2G_Multiome_predictions_genes")
              and fr["file_set"] == "jesse-engreitz:ds1_k562_scE2G_Multiome_predictions" and fr["controlled_access"] == "false"
              and fr["reference_files"] == "IGVFFI7969JLFC,IGVFFI0653VCGH,IGVFFI9573KOZR" and fr["analysis_step_version"] == RUN_SCE2G_ASV,
              "full predictions: derived_from chain, file_set, controlled vocab, reference files")
        ps = read_rows_tsv(mdir / "ds1" / "round3_prediction_set_post.tsv")[0]
        pp = read_rows_tsv(mdir / "ds1" / "round2_principal_pseudobulk_set_post.tsv")
        pp = next(r for r in pp if "k562" in r["aliases"])
        check(ps["file_set_type"] == "element-gene links" and ps["samples"] == f"{S1},{S2}" and ps["cell_type"] == "/sample-terms/EFO_0002067/"
              and ps["input_file_sets"] == "jesse-engreitz:scE2G_Multiome_model,jesse-engreitz:ds1_k562_filtered_pseudobulk_set"
              and pp["input_file_sets"] == f"anshul-kundaje:ds1-k562-{S1},anshul-kundaje:ds1-k562-{S2}" and pp["merged"] == "true",
              "prediction / principal pseudobulk sets: samples by frequency, Kundaje primary aliases, sample term")
        at = read_rows_tsv(mdir / "ds1" / "round4_filtered_atac_fragment_file_post.tsv")[0]
        check(at["file_format_type"] == "bed3+" and at["derived_from"].endswith("jesse-engreitz:ds1_k562_filtered_barcode_list")
              and at["reference_files"] == "IGVFFI0653VCGH,IGVFFI6788CPPS", "filtered ATAC fragment file: bed3+, derived_from, reference files")
        lin = json.loads((Path(o("man1")) / "lineage.json").read_text())
        check(lin["resolved"][f"anshul-kundaje:ds1-k562-{S1}"]["ok"] and f"anshul-kundaje:ds1-k562-{S2}" in lin["missing"]
              and "IGVFDS0001PSBK" in lin["walks"]["IGVFDS1105KTIQ"]["pseudobulk_sets"],
              "lineage: derived_from / input_file_sets resolved via the Portal graph (portal_lineage.walk), missing refs flagged")
        up_plan = json.loads((Path(o("man1")) / "upload_plan.json").read_text())
        rounds = [r["round"] for r in up_plan["rows"]]
        check(not up_plan["executed"] and rounds == sorted(rounds) and up_plan["rows"][0]["rest_payload"]["aliases"] == ["jesse-engreitz:ds1_k562_QC_thresholds"]
              and isinstance(next(r for r in up_plan["rows"] if r["table"] == "prediction_tabular_files")["rest_payload"]["derived_from"], list),
              "upload_plan.json: dry-run, round order, REST payloads with arrays split")

        print("\nexecute path (in-memory portal, no network)")
        saved = {k: os.environ.pop(k, None) for k in ("IGVF_ACCESS_KEY", "IGVF_SECRET_ACCESS_KEY", "IGVF_API_KEY", "IGVF_SECRET_KEY")}
        try:
            make_target_portal(False)
            refused = False
        except SystemExit:
            refused = True
        finally:
            for kk, vv in saved.items():
                if vv is not None:
                    os.environ[kk] = vv
        check(refused, "--execute refuses without IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY")
        fake = _FakePortal()
        sdb2 = str(d / "state2.db")
        shutil.copy(sdb, sdb2)
        shutil.copy(sdb, d / "state3.db")
        res = run_manifest(config, cfgp, {("ds1", "k562")}, set(), sdb2, d / "manifests2", Path(o("man1")), PortalAliasReader(fake), fake)
        posts = [c for c in fake.calls if c[0] == "POST"]
        order = [next(e["round"] for e in res["plan"]["posts"] if e["alias"] == c[2]["aliases"][0]) for c in posts]
        body_qc = next(c[2] for c in posts if "QC_thresholds" in c[2]["aliases"][0])
        body_full = next(c[2] for c in posts if c[2]["aliases"][0].endswith("predictions_full"))
        check(len(res["executed"]["uploaded"]) == 15 and not res["executed"]["deferred"] and order == sorted(order),
              "execute: 15 objects POSTed in dependency (round) order in one pass, each re-read by alias")
        check(body_qc["attachment"]["href"].startswith("data:text/tab-separated-values;base64,") and len(body_full["md5sum"]) == 32
              and isinstance(body_full["derived_from"], list) and "file_format_type" not in body_full,
              "REST payload: attachment inlined as data URI, md5sum computed, arrays split")
        res2 = run_manifest(config, cfgp, {("ds1", "k562")}, set(), sdb2, d / "manifests3", Path(o("man1")), PortalAliasReader(fake))
        check(res2["outcomes"] == {"unchanged": 15} and not list((d / "manifests3" / "ds1").glob("round*"))
              and len(read_iu_tsv(d / "manifests3" / "ds1" / "tabular_file.tsv")) == 7, "re-plan: all unchanged, accumulator per object type")
        os.remove(cdir / f"ds1_k562_scE2G_{MULTIOME_MODEL}_threshold0.177.bedpe.gz.tbi")
        c3 = json.loads(json.dumps(config))
        c3["igvf"]["enabled_families"] = ["Multiome", "scATAC"]
        res3 = run_manifest(c3, cfgp, {("ds1", "k562")}, set(), str(d / "state3.db"), d / "manifests4", Path(o("man1")), OfflineReader())
        oc3 = Counter((r["table"], r["variant"], r["model"], r["outcome"]) for r in res3["plan"]["coverage"])
        check(oc3[("prediction_tabular_files", "genes", SCATAC_MODEL, "skipped-family-gated")] == 1
              and oc3[("bedpe_index_file", "bedpe_index", MULTIOME_MODEL, "skipped-missing-file")] == 1
              and oc3[("prediction_tabular_files", "full", SCATAC_MODEL, "planned-post")] == 1,
              "enabled_families scATAC: genes family-gated; missing .tbi -> skipped-missing-file")
        sc_full = next(e for e in res3["plan"]["posts"] if e["model"] == SCATAC_MODEL and e["variant"] == "full")
        check("genes" not in sc_full["payload"]["derived_from"] and ("filtered_rna_count_matrix", "", "") not in sc_full["deps"],
              "scATAC full: no RNA matrix / genes in derived_from or dependencies")

        print("\nstale reformats, report, submitter comment, Synapse, distance-depth")
        rows = [{"dataset": "ds1", "cluster": "k562", "cell_annotation": "K-562 corrected", "term_name": "K562", "term_id": "EFO:0002067"}]
        stale = find_stale(W["results"] / "uniformly_processed", "ds1", rows)
        n_el = len(stale)
        removed = remove_stale([s for s in stale if s[0].endswith("_element_list.bed.gz")])
        check(n_el == 6 and any(p.endswith(".tbi") for p in removed), "stale-reformats: CellAnnotation change flags every reformatted file; .tbi removed too")
        to_patch, cnt = plan_submitter_comment([{"aliases": ["jesse-engreitz:igvf4_a_scE2G_Multiome_predictions"], "submitter_comment": None, "uuid": "u1"},
                                                {"aliases": ["jesse-engreitz:igvf4_b_scE2G_Multiome_predictions"], "submitter_comment": "note", "uuid": "u2"},
                                                {"aliases": ["jesse-engreitz:igvf4_c_scE2G_Multiome_predictions"], "submitter_comment": "(Version 1) - x"},
                                                {"aliases": ["jesse-engreitz:igvf9_d_scE2G_Multiome_predictions"], "submitter_comment": None}], {"igvf4"})
        check([e["new"] for e in to_patch] == ["Version 1", "(Version 1) - note"] and cnt["unchanged"] == 1 and cnt["skipped-wrong-dataset"] == 1,
              "patch-submitter-comment: never stacks the marker, prefixes free text once")
        keys = {("ds1", "h7"), ("ds1", "h7_mesoderm")}
        check(find_cluster_key("/r/ds1_h7_mesoderm_x.tsv.gz", keys) == ("ds1", "h7_mesoderm")
              and find_cluster_key("/r/ds1_h7_mesoderm_x.tsv.gz", keys, legacy_first_match=True) == ("ds1", "h7"),
              "synapse: longest cluster-name match (synapse-submission branch: first match)")
        inv = d / "inv.tsv"
        write_rows_tsv(inv, ["dataset", "cluster", "name", "id", "modifiedBy", "modifiedOn", "created_by", "modified_by"],
                       [{"dataset": "ds1", "cluster": "k562", "name": "old_multiome_powerlaw_v3.tsv.gz", "id": "syn1", "created_by": "u"},
                        {"dataset": "ds1", "cluster": "k562", "name": f"ds1_k562_scE2G_{MULTIOME_MODEL}.e2g.tsv.gz", "id": "syn2", "created_by": "u"},
                        {"dataset": "ds1", "cluster": "k562", "name": "other_model.tsv", "id": "syn3", "created_by": "v"},
                        {"dataset": "igvf7", "cluster": "x", "name": "a.tsv", "id": "syn4", "created_by": "u"}], quiet=True)
        rc = run(["synapse-manifest", "--product", "predictions", "--cluster-keys", "ds1/k562", "--inventory", str(inv), "--out-dir", o("syn"),
                  "--files", str(cdir / f"ds1_k562_scE2G_{MULTIOME_MODEL}.e2g.tsv.gz"), str(cdir / f"ds1_k562_scE2G_{SCATAC_MODEL}.e2g.tsv.gz")])
        sy = json.loads((Path(o("syn")) / "summary.json").read_text())
        check(rc == 0 and sy["new"] == 1 and sy["overwrite"] == 1 and sy["stale"] == 1, "synapse-manifest: new / overwrite / stale, other models invisible")
        write_rows_tsv(d / "inscope.tsv", ["dataset", "cluster"], [{"dataset": "ds1", "cluster": "k562"}], quiet=True)
        run(["synapse-orphans", "--mode", "cluster", "--inventory", str(inv), "--in-scope-tsv", str(d / "inscope.tsv"), "--out-dir", o("orph")])
        orph = read_rows_tsv(Path(o("orph")) / "synapse_orphans.tsv")
        check([r["dataset"] for r in orph] == ["igvf7"], "synapse-orphans: cluster folders with no in-scope replacement")
        stats_p = d / "dd_stats.tsv"
        srows = [{"cluster": f"c{i}", "model_name": SCATAC_MODEL, "fragments_total": f, "cell_count": 100, "mean_dist_to_tss": 20000 + 3000 * math.log10(f)}
                 for i, f in enumerate([1e6, 3e6, 1e7, 3e7, 1e8])]
        write_rows_tsv(stats_p, list(srows[0]), srows, quiet=True)
        write_rows_tsv(d / "pmap.tsv", ["cluster", "path"], [{"cluster": "c2", "path": str(cdir / SCATAC_MODEL / "scE2G_predictions.tsv.gz")}], quiet=True)
        rc = run(["distance-depth", "--stats", str(stats_p), "--predictions-map", str(d / "pmap.tsv"), "--out-dir", o("dd")] + (["--no-plots"] if args.no_plots else []))
        dd = json.loads((Path(o("dd")) / "summary.json").read_text())
        ddt = {r["cluster"]: r for r in read_rows_tsv(Path(o("dd")) / "distance_vs_depth.tsv")}
        check(rc == 0 and abs(dd["trends"]["mean_dist_to_tss"][0] - 3000) < 1e-6 and dd["trends"]["mean_dist_to_tss"][2] == 4
              and ddt["c2"]["median_dist_to_tss"] != "nan", "distance-depth: log10 trend over clusters >= 2e6 fragments, median/IQR from predictions")
        rc = run(["report", "--config", cfgp, "--state-db", sdb, "--manifest-dir", str(W["results"] / "igvf_manifests"), "--out-dir", o("rep")])
        check(rc == 0 and len(read_rows_tsv(Path(o("rep")) / "report.tsv")) == 6, "report: one row per cluster with stats")

    after = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
    for p in after - before:
        shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink()
    n_fail = sum(1 for ok, _ in checks if not ok)
    print(f"\n{len(checks) - n_fail}/{len(checks)} checks pass in {time.time() - t0:.1f}s")
    if n_fail:
        print("selftest: FAILED")
        return 1
    print("selftest: all checks pass")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent e2g-qc-predictions", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, label, config=False, required_config=False):
        p.add_argument("--label", default=label, help="run label (output dir suffix)")
        p.add_argument("--out-dir", help="write here instead of Docs/E2GQCPredictions/<ts>_<label>")
        if config:
            p.add_argument("--config", required=required_config, help="*_pipeline_config.yaml (or .json)")
            p.add_argument("--output-dir", help="override the config's output_dir")
        return p

    p = common(sub.add_parser("resolve-exclusions", help="quality gate -> cluster_stats + pipeline plan"), "exclusions", True, True)
    p.add_argument("--legacy-implicit-source", action="store_true", help="synapse-submission / CATlas stats-source rule")
    p.add_argument("--write-to-output-dir", action="store_true", help="write cluster_stats/ and igvf_metadata/ under output_dir")
    p.set_defaults(func=cmd_resolve_exclusions)

    p = sub.add_parser("merge-metrics", help="merge a merged cluster's component filtered_cell_subsample_metrics.tsv")
    p.add_argument("--metrics", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_merge_metrics)

    p = sub.add_parser("prefiltered-metrics", help="CATlas: metrics row from an already-filtered fragments file")
    p.add_argument("--frag-file", required=True)
    p.add_argument("--subsample-name", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_prefiltered_metrics)

    p = sub.add_parser("build-qc-datatables", help="per-cluster per_cell_qc tables from downloaded pseudobulks")
    p.add_argument("--pseudobulks-root", required=True)
    p.add_argument("--dest", required=True, help="QC_datatables/ is created under here")
    p.add_argument("--datasets", nargs="+")
    p.add_argument("--config", help="pipeline config supplying merged-cluster annotations")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_build_qc_datatables)

    p = sub.add_parser("filter-atac", help="QC-filter + sort + bgzip/tabix ATAC fragments for one cluster")
    for a in ("--qc-guide", "--pseudobulks", "--cell-type", "--chrom-sizes", "--out"):
        p.add_argument(a, required=True)
    p.add_argument("--allow-missing-barcodes", action="store_true", help="legacy branches: warn instead of failing")
    p.add_argument("--lenient-fields", action="store_true", help="legacy branches: skip records with < 4 fields")
    p.add_argument("--clean", action="store_true", help="write the 16-bp barcode only (legacy option)")
    p.set_defaults(func=cmd_filter_atac)

    p = sub.add_parser("filter-rna", help="QC-filter + concatenate RNA h5ad matrices, Ensembl -> symbol collapse")
    for a in ("--qc-guide", "--pseudobulks", "--cell-type", "--out"):
        p.add_argument(a, required=True)
    p.add_argument("--ensembl-ids-as-genes", action="store_true")
    p.add_argument("--gtf")
    p.add_argument("--standard-chromosomes-only", action="store_true")
    p.add_argument("--log")
    p.add_argument("--package-tar", help="also write the Filtered Matrix File tar.gz here (mtx output only)")
    p.set_defaults(func=cmd_filter_rna)

    p = sub.add_parser("package-rna", help="flat tar.gz of decompressed matrix.mtx / barcodes.tsv / features.tsv")
    p.add_argument("--matrix-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_package_rna)

    p = common(sub.add_parser("sce2g-config", help="scE2G cell_clusters / cluster_metadata tables + config overlay"), "sce2g_config", True, True)
    p.add_argument("--lab-annotations", help="lab_annotations_with_cl.tsv (dataset, lab_celltype, CL term, qualifier, CL_ID)")
    p.set_defaults(func=cmd_sce2g_config)

    p = sub.add_parser("reformat", help="portal-format prediction / list / element BED / bedpe file")
    p.add_argument("--kind", choices=["auto", "element-bed", "bedpe"], default="auto",
                   help="auto: full / thresholded / gene list / element list from the input name and --threshold")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--model", default=MULTIOME_MODEL)
    p.add_argument("--method", help="default scE2G_<model>")
    p.add_argument("--version", default="1.2")
    p.add_argument("--cell-type", help="SampleTermName (portal cell_type.term_name)")
    p.add_argument("--term-id", help="SampleTermID (CURIE)")
    p.add_argument("--summary", help="CellAnnotation (SampleSummaryShort in --format synapse-legacy)")
    p.add_argument("--threshold", help="score threshold; header 'ScoreThreshold: Score >= t'")
    p.add_argument("--portal-link", help="the file's IGVF alias")
    p.add_argument("--all-columns", action="store_true")
    p.add_argument("--format", choices=["portal", "synapse-legacy"], default="portal")
    p.add_argument("--catlas", action="store_true", help="CATlas branch: scATAC full files without normalizedATAC_enh")
    p.set_defaults(func=cmd_reformat)

    p = sub.add_parser("candidates", help="candidate E2G pairs table")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    p.set_defaults(func=cmd_candidates)

    p = sub.add_parser("features", help="scE2G feature table with header")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--cell-type")
    p.add_argument("--term-id")
    p.add_argument("--summary")
    p.set_defaults(func=cmd_features)

    p = common(sub.add_parser("aggregate-qc", help="dataset-wide all_qc_stats.tsv + QC figures"), "aggregate_qc")
    p.add_argument("--dataset-dir", required=True)
    p.add_argument("--plots-dir", required=True, help="plots/<dataset> (metrics files for cell_count patching)")
    p.add_argument("--out")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_aggregate_qc)

    p = sub.add_parser("stale-reformats", help="reformatted files whose embedded annotation disagrees with the cache")
    p.add_argument("--results-dir", required=True, help="<output_dir>/uniformly_processed")
    p.add_argument("--datasets", nargs="+", required=True)
    p.add_argument("--state-db", default=str(DEFAULT_STATE_DB))
    p.add_argument("--annotations-tsv")
    p.add_argument("--delete", action="store_true")
    p.set_defaults(func=cmd_stale_reformats)

    p = common(sub.add_parser("cell-metadata", help="IGVF Portal Cell Annotation cache + per-cluster derivation + snapshot"), "cell_metadata", True)
    p.add_argument("--multireport-json", help="saved PseudobulkSet multireport (offline)")
    p.add_argument("--state-db", default=str(DEFAULT_STATE_DB))
    p.add_argument("--offline", action="store_true", help="never GET; derive from the cache")
    p.add_argument("--ttl-hours", type=float, default=CACHE_TTL_HOURS)
    p.add_argument("--write-snapshot", action="store_true", help="write the snapshot under output_dir/igvf_metadata")
    p.add_argument("--catlas-seed", action="store_true", help="seed CATlas (WashU) clusters from the lab-filtered multireport")
    p.add_argument("--dataset", help="dataset label for --catlas-seed (default catlas)")
    p.set_defaults(func=cmd_cell_metadata)

    p = common(sub.add_parser("cell-annotation-report", help="igvf_cell_annotation_report TSV"), "cell_annotation_report")
    p.add_argument("--multireport-json")
    p.add_argument("--qc-guide-dir")
    p.add_argument("--out")
    p.set_defaults(func=cmd_cell_annotation_report)

    p = common(sub.add_parser("dataset-accessions", help="dataset -> principal analysis set accession(s)"), "dataset_accessions")
    p.add_argument("--multireport-json")
    p.add_argument("--json-out")
    p.set_defaults(func=cmd_dataset_accessions)

    p = common(sub.add_parser("washu-report", help="CATlas WashU PseudobulkSet report"), "washu_report")
    p.add_argument("--multireport-json")
    p.add_argument("--lab", default=WASHU_LAB)
    p.add_argument("--out")
    p.set_defaults(func=cmd_washu_report)

    p = sub.add_parser("verify-fragments", help="md5-verify (and with --download fetch) fragments listed in a WashU report")
    p.add_argument("--report", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--download", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--summary")
    p.set_defaults(func=cmd_verify_fragments)

    p = common(sub.add_parser("portal-files", help="primary-pseudobulk file discovery and download plan"), "portal_files")
    p.add_argument("--multireport-json")
    p.add_argument("--lab", default=DEFAULT_FILE_LAB)
    p.add_argument("--content-types", nargs="+")
    p.add_argument("--datasets", nargs="+")
    p.add_argument("--download-root")
    p.add_argument("--download", action="store_true", help="fetch files marked todo (md5-verified)")
    p.set_defaults(func=cmd_portal_files)

    p = common(sub.add_parser("compare-archive", help="portal-downloaded vs archive pseudobulk files"), "compare_archive")
    p.add_argument("--discovery", required=True, help="portal_files.tsv from portal-files")
    p.add_argument("--archive-root", required=True)
    p.add_argument("--no-characterize", action="store_true")
    p.set_defaults(func=cmd_compare_archive)

    p = common(sub.add_parser("manifest", help="IGVF Portal payloads / round TSVs / upload plan; --execute submits"), "manifest", True, True)
    p.add_argument("--cluster-keys", required=True, help="'dataset/cluster,...' or a bare dataset (manifest-eligible expansion)")
    p.add_argument("--excluded-cluster-keys", default="")
    p.add_argument("--state-db")
    p.add_argument("--manifest-dir")
    p.add_argument("--tables", nargs="+", choices=sorted(TABLES))
    p.add_argument("--offline", action="store_true", help="no alias lookups (every row not in the ledger is a POST)")
    p.add_argument("--live-aliases", help="JSON {alias: record} treated as already live (offline PATCH planning)")
    p.add_argument("--lineage", action="store_true", help="resolve every external reference with an authenticated GET")
    p.add_argument("--lineage-walk", action="store_true", help="also walk the Portal graph from igvf.principal_analysis_sets")
    p.add_argument("--lineage-fixture", help="JSON {path: object} answering lineage GETs offline")
    p.add_argument("--execute", action="store_true", help="POST/PATCH for real (credentials from IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY)")
    p.add_argument("--production", action="store_true", help="target api.data.igvf.org instead of the sandbox")
    p.set_defaults(func=cmd_manifest)

    p = common(sub.add_parser("patch-submitter-comment", help="'Version 1' submitter_comment backfill for live Prediction Sets"), "submitter_comment")
    p.add_argument("--multireport-json")
    p.add_argument("--dataset", nargs="+")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--production", action="store_true")
    p.set_defaults(func=cmd_patch_submitter_comment)

    p = common(sub.add_parser("report", help="per-cluster coverage report.tsv"), "coverage_report", True, True)
    p.add_argument("--state-db", default=str(DEFAULT_STATE_DB))
    p.add_argument("--annotations-tsv")
    p.add_argument("--manifest-dir")
    p.add_argument("--out")
    p.set_defaults(func=cmd_report)

    p = common(sub.add_parser("synapse-manifest", help="Synapse path/parent manifest diff for one product"), "synapse")
    p.add_argument("--product", choices=sorted(SYNAPSE_PARENTS), required=True)
    p.add_argument("--cluster-keys", required=True)
    p.add_argument("--files", nargs="+")
    p.add_argument("--inventory", help="TSV of owned Synapse entities (dataset, cluster, name, id, modifiedBy, modifiedOn)")
    p.add_argument("--manifest-out")
    p.add_argument("--parent-id")
    p.add_argument("--legacy-first-match", action="store_true")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--confirm-delete", action="store_true")
    p.add_argument("--confirm-overwrite", action="store_true")
    p.set_defaults(func=cmd_synapse_manifest)

    p = common(sub.add_parser("synapse-orphans", help="read-only Synapse orphan listing (cluster or file level)"), "synapse_orphans")
    p.add_argument("--mode", choices=["cluster", "file"], default="cluster")
    p.add_argument("--inventory", required=True)
    p.add_argument("--in-scope-tsv", help="report.tsv (cluster mode)")
    p.add_argument("--manifest", help="predictions manifest TSV (file mode)")
    p.add_argument("--exclude-top", nargs="+")
    p.add_argument("--only", nargs="+")
    p.add_argument("--user-id")
    p.set_defaults(func=cmd_synapse_orphans)

    p = common(sub.add_parser("distance-depth", help="CATlas distance-to-TSS vs sequencing depth"), "distance_depth")
    p.add_argument("--stats", nargs="+", required=True, help="all_qc_stats.tsv file(s)")
    p.add_argument("--predictions-map", help="TSV cluster, path (thresholded predictions) for median / IQR")
    p.add_argument("--group-map", help="TSV cluster, group")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_distance_depth)

    p = common(sub.add_parser("run", help="driver: preflight, warm, local packaging, manifest preview, audit"), "run", True, True)
    p.add_argument("--mode", choices=["default", "local_only"])
    p.add_argument("--multireport-json")
    p.add_argument("--state-db")
    p.add_argument("--offline", action="store_true")
    p.add_argument("-n", "--dry-run", action="store_true")
    p.add_argument("--legacy-implicit-source", action="store_true")
    p.add_argument("--lab-annotations")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic end-to-end self-test")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    if args.cmd != "selftest" and not _IN_SELFTEST[0]:
        try:
            setup_logging()
        except OSError:
            pass
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
