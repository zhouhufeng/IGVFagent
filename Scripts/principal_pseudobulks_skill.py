#!/usr/bin/env python3
"""Principal pseudobulks from IGVF single-cell multiome data (port of EngreitzLab/generate-principal-pseudobulks).

Port of https://github.com/EngreitzLab/generate-principal-pseudobulks (MIT, Kayla Brand 2026; Snakemake +
Python + R) at commit 80222ccfc026af86fe346f2d2a5992efb1e68cdd (main, 2026-08-18T02:42:22Z).  The
upstream repository QC-filters IGVF multiome *primary* pseudobulks (per cell type x subsample
directories holding fragments.tsv.gz, rna_counts_mtx.h5ad and per_cell_qc.tsv) into released
*principal* pseudobulks: one QC-filtered barcode list, one sorted + tabix-indexed ATAC fragment file
and one gene-symbol RNA count matrix per cluster, plus the per-dataset config table scE2G reads.
Every file was read -- Snakefile, workflow/scripts/filter_atac_fragments.py and filter_rna_counts.py,
plotting_scripts/build_per_cell_qc_datatable.py, explore_qc_thresholds.R and plot_per_cell_qc.R, the
config template, the chrom-sizes reference (embedded here) and the two file specs -- and the
definitions re-derived in Python (numpy / pandas / scipy / anndata / pysam / matplotlib).  No code was
copied.  Relationship: port.  The step upstream takes as given -- building the primary pseudobulks
from IGVF Portal files -- is added here (`fetch`, `build-pseudobulks`), processed-first through
portal_lineage (principal analysis set -> cell annotations; its uniform-pipeline intermediate
analysis sets -> fragments + h5ad), so the chain runs from an IGVF accession.

Definitions, as upstream computes them
  per_cell_qc.tsv   analysis_accession, barcode, subsample, rna_read_count, gene_count, pct_mito,
                    pct_ribo, num_frags, pct_duplicated_reads, nucleosomal_signal, tss_enrichment, frip.
  qc-datatable      (build_per_cell_qc_datatable.py) concatenate per_cell_qc.tsv of every
                    annotation-{cell_type}-IGVF* directory (sorted glob), falling back to
                    annotation-{cell_type}; one header, 0-byte files skipped, header mismatch = error.
  explore           (explore_qc_thresholds.R) defaults rna_min 1000, rna_max Inf, gene_min 1000,
                    gene_max Inf, pct_mt_max 30, pct_ribo_max 100, atac_min 1000, atac_max Inf,
                    tss_enr_min 3, nuc_signal_max 1.5, pct_dup_max 100, frip_min 0.  A cell passes
                    with rna >= min, rna <= max, genes >= / <=, pct_mito < , pct_ribo < , frags >= / <=,
                    nucleosomal_signal < , tss_enrichment >= , pct_dup < , frip >= .  Per threshold:
                    Total = cells failing it, Alone = cells failing only it; always shown: rna_min,
                    atac_min, tss_enr_min, nuc_signal_max, pct_mt_max, plus any overridden or
                    dropping cells; "most stringent" = max Total; --sets "a; b; c"; per-subsample
                    table (kept, dropped, % drop, most stringent threshold).
  qc-filter         (plot_per_cell_qc.R) same flags, defaults as above except gene_min 0, and STRICT
                    inequalities: rna > min, rna < max, genes > / <, pct_mito < , pct_ribo < ,
                    frags > / <, nuc < , tss > , dup < , frip > .  Writes
                    filtered_barcodes_with_subsamples.tsv.gz (barcode, subsample, analysis_accession,
                    rows in subsample order), qc_thresholds.tsv (12 named thresholds, R number format,
                    Inf), filtered_cell_subsample_metrics.tsv (subsample, n_cells, total_fragments,
                    total_RNA_reads, mean_frag_per_cell, mean_RNA_per_cell, mean_frip, mean_tss) and the
                    RNA / ATAC QC figures and cells-per-subsample bars after RNA, ATAC and all QC
                    (PNG here; nUMI_non_MT = rna_read_count x (100 - pct_mito) / 100).
  filter-atac       (filter_atac_fragments.py) QC-guide barcodes (column 1, full string, never
                    truncated to 16 bp) matched against column 4 of every fragments.tsv.gz of the
                    cell type's directories (comma-separated cell types merge clusters); records must
                    have exactly 5 fields; every guide barcode must be seen (else error); chromosomes
                    absent from the chrom-sizes file dropped with a warning; output in chrom-sizes
                    order, sort-bed order within a chromosome (start, end, rest), bgzip + tabix -p bed.
  filter-rna        (filter_rna_counts.py) rna_counts_mtx.h5ad per directory filtered on the full
                    obs_name, concatenated (inner join, merge="same"); Ensembl -> symbol through the
                    GTF (exact versioned gene_id, first gene_name wins, conflicts logged; any unmatched
                    ID is a hard failure; --standard-chromosomes-only keeps chr1-22, X, Y, M) or through
                    var["gene_symbol"]; counts of IDs sharing a symbol summed (symbols sorted, np.unique);
                    sanity checks: no duplicate barcodes or genes, cell count == QC guide; outputs .mtx
                    directory (matrix.mtx.gz genes x cells, barcodes.tsv.gz, features.tsv.gz), .csv.gz
                    (cells x genes) or .h5ad; <out>_gtf_mapping.txt log.
  config-table      (rule make_config_table) cluster, rna_matrix_file, atac_frag_file, HiC_file,
                    HiC_type, HiC_resolution, alt_TSS, alt_genes, model_dir per dataset.
  run               (Snakefile) every (dataset, cell_type) of the config: filter-atac + filter-rna
                    (--gtf transcriptome --standard-chromosomes-only) into
                    {out}/{dataset}/{cell_type}/atac_fragments_{dataset}_{cell_type}.tsv.gz(.tbi) and
                    rna_count_matrix_{dataset}_{cell_type}/, then {out}/config/tables/{dataset}_config.tsv;
                    memory estimate max(4 x input MB (x8 if gz), 8 GB) x 2^(attempt-1), cap 250 GB.

Added here (not in upstream, which starts from primary pseudobulks)
  fetch             portal_lineage.walk + build_plan from any IGVF accession: principal analysis
                    set(s) and their cell-annotation table, and per uniform-pipeline intermediate
                    analysis set (one 10x lane) the best fragments file and RNA h5ad; writes
                    portal_manifest.tsv; downloads with raw_data_pipeline.portal_download within
                    --max-gb (--dry-run lists only).
  build-pseudobulks annotation table + per-lane fragments / h5ad -> primary pseudobulk layout
                    {out}/{dataset}/pseudobulks/annotation-{cell_type}-{subsample}/.  Barcodes become
                    the spec's {16-bp barcode}_{IGVF 10x lane accession}: a 10x GEM suffix (-1) is
                    stripped, ATAC barcodes are translated to RNA (gene-expression) barcodes with
                    --atac-rna-map (10x multiome ATAC<->GEX whitelist pairs), then the lane accession
                    is appended unless a suffix is already present.  per_cell_qc metrics (clean-room,
                    Seurat / Signac definitions): rna_read_count = UMIs, gene_count = genes > 0,
                    pct_mito / pct_ribo = % UMIs in ^MT- / ^RP[SL] genes; num_frags = fragment rows,
                    pct_duplicated_reads = 100 x (sum(col5) - n) / sum(col5); nucleosomal_signal =
                    fragments 147-293 bp / fragments < 147 bp; tss_enrichment (with --tss) = mean
                    Tn5-insertion count per bp in TSS -500..+499 / mean in the +/-(901-1000) flanks,
                    strand-aware (Signac TSSEnrichment; zero flanks -> mean flank); frip (with
                    --peaks) = fraction of fragments overlapping a peak.  A metric that cannot be
                    computed is NA and qc-filter / explore skip an all-NA column (upstream would drop
                    every cell; --upstream-compat keeps the upstream behaviour).
  package-rna       the spec's .tar.gz (matrix.mtx, barcodes.tsv, features.tsv, decompressed).
  validate          every guarantee of FILE_SPEC_QC_FILTERED_BARCODE_LIST.md and
                    FILE_SPEC_RNA_COUNT_MATRIX.md (+ fragments: guide barcodes present, 5 columns,
                    chrom-sizes order, tabix index).
  sce2g-prep        config table -> sce2g-pipeline inputs: tagAlign + fragment_count + cell barcodes
                    (sce2g_pipeline_skill.frag_to_tagalign), RNA pseudobulk TPM / detection / mean
                    log-norm (sce2g_pipeline_skill.rna_features), and config_cell_clusters.tsv for
                    `sce2g-pipeline run --cluster-config`.

Deviations (documented)
  * matrix.mtx is written with the `integer` field the RNA spec states (upstream's scipy mmwrite of a
    float matrix writes `real`); --upstream-compat writes `real`.
  * config-table model_dir lists the embedded sce2g-pipeline model names
    (multiome_powerlaw_v3,scATAC_powerlaw_v3) so the table runs as is; --upstream-compat writes
    upstream's models/... paths.
  * README documents filter_atac_fragments.py --clean (16-bp barcode only) but the script lacks it;
    implemented as documented, off by default.
  * The R scripts' NA comparisons (dplyr drops NA rows) are reproduced in the filters; in explore's
    Total / Alone counts an NA metric counts as failing (R prints NA).
  * Directory discovery also accepts annotation-{cell_type}-* when neither upstream pattern exists
    (subsample labels that are not IGVF accessions).

Subcommands
  fetch, build-pseudobulks, qc-datatable, explore, qc-filter, filter-atac, filter-rna, package-rna,
  config-table, run, validate, sce2g-prep, selftest.

Output: Docs/PrincipalPseudobulks/<timestamp>_<label>/.

Usage:
    igvfagent principal-pseudobulks fetch --accession IGVFDS1244UUGQ --dry-run --max-gb 20
    igvfagent principal-pseudobulks build-pseudobulks --annotations cells.tsv --manifest portal_manifest.tsv --dataset igvf1 --tss tss.bed --peaks peaks.bed
    igvfagent principal-pseudobulks qc-datatable --pseudobulks igvf1/pseudobulks --cell-type k562
    igvfagent principal-pseudobulks explore --meta k562_per_cell_qc.tsv --sets "--tss-min 3" "--tss-min 5 --pct-mt-max 20" --show-subsamples
    igvfagent principal-pseudobulks qc-filter --meta k562_per_cell_qc.tsv --rna-min 1000 --gene-min 1000 --tss-min 3
    igvfagent principal-pseudobulks filter-atac --qc-guide filtered_barcodes_with_subsamples.tsv.gz --pseudobulks igvf1/pseudobulks --cell-type k562
    igvfagent principal-pseudobulks filter-rna --qc-guide filtered_barcodes_with_subsamples.tsv.gz --pseudobulks igvf1/pseudobulks --cell-type k562 --gtf IGVFFI9573KOZR.gtf.gz --standard-chromosomes-only
    igvfagent principal-pseudobulks run --config config_QC_pseudobulks.yaml
    igvfagent principal-pseudobulks validate --qc-guide g.tsv.gz --rna-matrix rna_count_matrix_igvf1_k562 --fragments atac_fragments_igvf1_k562.tsv.gz
    igvfagent principal-pseudobulks sce2g-prep --config-table multiome_data/config/tables/igvf1_config.tsv
    igvfagent principal-pseudobulks selftest --no-plots
"""
from __future__ import annotations

import argparse
import glob
import gzip
import io
import json
import logging
import math
import os
import re
import shutil
import sys
import tarfile
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "PrincipalPseudobulks"
DATA_ROOT = ROOT / "Data" / "PrincipalPseudobulks"

UPSTREAM_REPO = "EngreitzLab/generate-principal-pseudobulks"
UPSTREAM_COMMIT = "80222ccfc026af86fe346f2d2a5992efb1e68cdd"
UPSTREAM_DATE = "2026-08-18T02:42:22Z"

PER_CELL_QC_COLS = ["analysis_accession", "barcode", "subsample", "rna_read_count", "gene_count", "pct_mito",
                    "pct_ribo", "num_frags", "pct_duplicated_reads", "nucleosomal_signal", "tss_enrichment", "frip"]
GUIDE_COLS = ["barcode", "subsample", "analysis_accession"]
CONFIG_TABLE_COLS = ["cluster", "rna_matrix_file", "atac_frag_file", "HiC_file", "HiC_type", "HiC_resolution",
                     "alt_TSS", "alt_genes", "model_dir"]
UPSTREAM_MODEL_DIR = "models/multiome_powerlaw_v3,models/scATAC_powerlaw_v3"
SCE2G_MODEL_DIR = "multiome_powerlaw_v3,scATAC_powerlaw_v3"
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
IGVF_TRANSCRIPTOME = "IGVFFI9573KOZR"   # GENCODE 43 / GRCh38 GTF the Snakefile passes to filter_rna_counts.py

# threshold key -> (per_cell_qc column, explore pass op, qc-filter pass op)
THRESH = {
    "rna_min": ("rna_read_count", ">=", ">"), "rna_max": ("rna_read_count", "<=", "<"),
    "gene_min": ("gene_count", ">=", ">"), "gene_max": ("gene_count", "<=", "<"),
    "pct_mt_max": ("pct_mito", "<", "<"), "pct_ribo_max": ("pct_ribo", "<", "<"),
    "atac_min": ("num_frags", ">=", ">"), "atac_max": ("num_frags", "<=", "<"),
    "tss_enr_min": ("tss_enrichment", ">=", ">"), "nuc_signal_max": ("nucleosomal_signal", "<", "<"),
    "pct_dup_max": ("pct_duplicated_reads", "<", "<"), "frip_min": ("frip", ">=", ">"),
}
THRESH_ORDER = ["rna_min", "rna_max", "gene_min", "gene_max", "pct_mt_max", "pct_ribo_max", "atac_min", "atac_max",
                "tss_enr_min", "nuc_signal_max", "pct_dup_max", "frip_min"]
RNA_KEYS = THRESH_ORDER[:6]
ATAC_KEYS = THRESH_ORDER[6:]
INF = float("inf")
EXPLORE_DEFAULTS = {"rna_min": 1e3, "rna_max": INF, "gene_min": 1000, "gene_max": INF, "pct_mt_max": 30,
                    "pct_ribo_max": 100, "atac_min": 1e3, "atac_max": INF, "tss_enr_min": 3, "nuc_signal_max": 1.5,
                    "pct_dup_max": 100, "frip_min": 0}
FILTER_DEFAULTS = dict(EXPLORE_DEFAULTS, gene_min=0)
ALWAYS_SHOW = ["rna_min", "atac_min", "tss_enr_min", "nuc_signal_max", "pct_mt_max"]
FLAG_MAP = {"--rna-min": "rna_min", "--rna-max": "rna_max", "--gene-min": "gene_min", "--gene-max": "gene_max",
            "--pct-mt-max": "pct_mt_max", "--pct-ribo-max": "pct_ribo_max", "--atac-min": "atac_min",
            "--atac-max": "atac_max", "--tss-min": "tss_enr_min", "--nuc-max": "nuc_signal_max",
            "--pct-dup-max": "pct_dup_max", "--frip-min": "frip_min"}
THRESHOLD_LABELS = {"rna_min": "RNA reads >= ", "rna_max": "RNA reads <= ", "gene_min": "genes >= ",
                    "gene_max": "genes <= ", "pct_mt_max": "% mito < ", "pct_ribo_max": "% ribo < ",
                    "atac_min": "ATAC frags >= ", "atac_max": "ATAC frags <= ", "tss_enr_min": "TSS enrichment >= ",
                    "nuc_signal_max": "nucleosomal signal < ", "pct_dup_max": "% dup reads < ", "frip_min": "FRIP >= "}
THRESHOLD_DOC = {"rna_min": "Minimum RNA reads per cell", "rna_max": "Maximum RNA reads per cell",
                 "gene_min": "Minimum genes per cell", "gene_max": "Maximum genes per cell",
                 "pct_mt_max": "Maximum % mitochondrial", "pct_ribo_max": "Maximum % ribosomal",
                 "atac_min": "Minimum ATAC fragments", "atac_max": "Maximum ATAC fragments",
                 "tss_enr_min": "Minimum TSS enrichment", "nuc_signal_max": "Maximum nucleosomal signal",
                 "pct_dup_max": "Maximum % duplicated reads", "frip_min": "Minimum FRIP"}
# plot_per_cell_qc.R palettes
CP = ["#429130", "#2f9a71", "#159594", "#0096a0", "#0083ab", "#0075b3", "#006eae", "#5b5da3", "#8d4b9b", "#a64791",
      "#b03e67", "#c5373d", "#d8571f", "#e96a00", "#ca9b23"]
GREYS = ["#e5e5e9", "#c5cad7", "#96a0b3", "#6e788d", "#435369", "#1c2a43"]
INK, INK2, AXIS = "#0b0b0b", "#52514e", "#d0cfca"

BARCODE_SPEC_RE = re.compile(r"^[ACGTN]{16}_IGVF[A-Z]{2}\d{4}[A-Z]{4}$")
IGVF_ACC_RE = re.compile(r"^IGVF[A-Z]{2}\d{4}[A-Z]{4}$")
GEM_SUFFIX_RE = re.compile(r"-\d+$")

# reference/IGVF.DACC.GRCh38.chrom.sizes.tsv at the pinned commit (195 contigs, file order = sort order)
_CHROM_SIZES_EMBED = """
chr1:248956422 chr2:242193529 chr3:198295559 chr4:190214555 chr5:181538259 chr6:170805979
chr7:159345973 chr8:145138636 chr9:138394717 chr10:133797422 chr11:135086622 chr12:133275309
chr13:114364328 chr14:107043718 chr15:101991189 chr16:90338345 chr17:83257441 chr18:80373285
chr19:58617616 chr20:64444167 chr21:46709983 chr22:50818468 chrX:156040895 chrY:57227415
chrM:16569 chr1_KI270706v1_random:175055 chr1_KI270707v1_random:32032 chr1_KI270708v1_random:127682 chr1_KI270709v1_random:66860 chr1_KI270710v1_random:40176
chr1_KI270711v1_random:42210 chr1_KI270712v1_random:176043 chr1_KI270713v1_random:40745 chr1_KI270714v1_random:41717 chr2_KI270715v1_random:161471 chr2_KI270716v1_random:153799
chr3_GL000221v1_random:155397 chr4_GL000008v2_random:209709 chr5_GL000208v1_random:92689 chr9_KI270717v1_random:40062 chr9_KI270718v1_random:38054 chr9_KI270719v1_random:176845
chr9_KI270720v1_random:39050 chr11_KI270721v1_random:100316 chr14_GL000009v2_random:201709 chr14_GL000225v1_random:211173 chr14_KI270722v1_random:194050 chr14_GL000194v1_random:191469
chr14_KI270723v1_random:38115 chr14_KI270724v1_random:39555 chr14_KI270725v1_random:172810 chr14_KI270726v1_random:43739 chr15_KI270727v1_random:448248 chr16_KI270728v1_random:1872759
chr17_GL000205v2_random:185591 chr17_KI270729v1_random:280839 chr17_KI270730v1_random:112551 chr22_KI270731v1_random:150754 chr22_KI270732v1_random:41543 chr22_KI270733v1_random:179772
chr22_KI270734v1_random:165050 chr22_KI270735v1_random:42811 chr22_KI270736v1_random:181920 chr22_KI270737v1_random:103838 chr22_KI270738v1_random:99375 chr22_KI270739v1_random:73985
chrY_KI270740v1_random:37240 chrUn_KI270302v1:2274 chrUn_KI270304v1:2165 chrUn_KI270303v1:1942 chrUn_KI270305v1:1472 chrUn_KI270322v1:21476
chrUn_KI270320v1:4416 chrUn_KI270310v1:1201 chrUn_KI270316v1:1444 chrUn_KI270315v1:2276 chrUn_KI270312v1:998 chrUn_KI270311v1:12399
chrUn_KI270317v1:37690 chrUn_KI270412v1:1179 chrUn_KI270411v1:2646 chrUn_KI270414v1:2489 chrUn_KI270419v1:1029 chrUn_KI270418v1:2145
chrUn_KI270420v1:2321 chrUn_KI270424v1:2140 chrUn_KI270417v1:2043 chrUn_KI270422v1:1445 chrUn_KI270423v1:981 chrUn_KI270425v1:1884
chrUn_KI270429v1:1361 chrUn_KI270442v1:392061 chrUn_KI270466v1:1233 chrUn_KI270465v1:1774 chrUn_KI270467v1:3920 chrUn_KI270435v1:92983
chrUn_KI270438v1:112505 chrUn_KI270468v1:4055 chrUn_KI270510v1:2415 chrUn_KI270509v1:2318 chrUn_KI270518v1:2186 chrUn_KI270508v1:1951
chrUn_KI270516v1:1300 chrUn_KI270512v1:22689 chrUn_KI270519v1:138126 chrUn_KI270522v1:5674 chrUn_KI270511v1:8127 chrUn_KI270515v1:6361
chrUn_KI270507v1:5353 chrUn_KI270517v1:3253 chrUn_KI270529v1:1899 chrUn_KI270528v1:2983 chrUn_KI270530v1:2168 chrUn_KI270539v1:993
chrUn_KI270538v1:91309 chrUn_KI270544v1:1202 chrUn_KI270548v1:1599 chrUn_KI270583v1:1400 chrUn_KI270587v1:2969 chrUn_KI270580v1:1553
chrUn_KI270581v1:7046 chrUn_KI270579v1:31033 chrUn_KI270589v1:44474 chrUn_KI270590v1:4685 chrUn_KI270584v1:4513 chrUn_KI270582v1:6504
chrUn_KI270588v1:6158 chrUn_KI270593v1:3041 chrUn_KI270591v1:5796 chrUn_KI270330v1:1652 chrUn_KI270329v1:1040 chrUn_KI270334v1:1368
chrUn_KI270333v1:2699 chrUn_KI270335v1:1048 chrUn_KI270338v1:1428 chrUn_KI270340v1:1428 chrUn_KI270336v1:1026 chrUn_KI270337v1:1121
chrUn_KI270363v1:1803 chrUn_KI270364v1:2855 chrUn_KI270362v1:3530 chrUn_KI270366v1:8320 chrUn_KI270378v1:1048 chrUn_KI270379v1:1045
chrUn_KI270389v1:1298 chrUn_KI270390v1:2387 chrUn_KI270387v1:1537 chrUn_KI270395v1:1143 chrUn_KI270396v1:1880 chrUn_KI270388v1:1216
chrUn_KI270394v1:970 chrUn_KI270386v1:1788 chrUn_KI270391v1:1484 chrUn_KI270383v1:1750 chrUn_KI270393v1:1308 chrUn_KI270384v1:1658
chrUn_KI270392v1:971 chrUn_KI270381v1:1930 chrUn_KI270385v1:990 chrUn_KI270382v1:4215 chrUn_KI270376v1:1136 chrUn_KI270374v1:2656
chrUn_KI270372v1:1650 chrUn_KI270373v1:1451 chrUn_KI270375v1:2378 chrUn_KI270371v1:2805 chrUn_KI270448v1:7992 chrUn_KI270521v1:7642
chrUn_GL000195v1:182896 chrUn_GL000219v1:179198 chrUn_GL000220v1:161802 chrUn_GL000224v1:179693 chrUn_KI270741v1:157432 chrUn_GL000226v1:15008
chrUn_GL000213v1:164239 chrUn_KI270743v1:210658 chrUn_KI270744v1:168472 chrUn_KI270745v1:41891 chrUn_KI270746v1:66486 chrUn_KI270747v1:198735
chrUn_KI270748v1:93321 chrUn_KI270749v1:158759 chrUn_KI270750v1:148850 chrUn_KI270751v1:150742 chrUn_KI270752v1:27745 chrUn_KI270753v1:62944
chrUn_KI270754v1:40191 chrUn_KI270755v1:36723 chrUn_KI270756v1:79590 chrUn_KI270757v1:71251 chrUn_GL000214v1:137718 chrUn_KI270742v1:186739
chrUn_GL000216v2:176608 chrUn_GL000218v1:161147 chrEBV:171823
"""

log = logging.getLogger("principal_pseudobulks")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"principal_pseudobulks_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label))[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    n = 1
    while d.exists():
        n += 1
        d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_{n}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir_for(args: argparse.Namespace, default_label: str) -> Path:
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
        raise SystemExit("principal-pseudobulks needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


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
        raise SystemExit("this step reads/writes .h5ad: pip install anndata") from e


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _sce2g():
    """sce2g_pipeline_skill, imported lazily (tagAlign, RNA pseudobulk features, bgzip/tabix, cluster config)."""
    try:
        import sce2g_pipeline_skill as s  # type: ignore
    except ImportError:
        from igvfagent import sce2g_pipeline_skill as s  # type: ignore
    return s


def _opener(path):
    p = str(path)
    with open(p, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(p, "rt") if gz else open(p, "rt")


def md_table(headers: Sequence[str], rows, max_rows: int = 60) -> str:
    rows = list(rows)
    out = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "---|" * len(headers)]
    for r in rows[:max_rows]:
        out.append("| " + " | ".join(_fmt(x) for x in r) + " |")
    if len(rows) > max_rows:
        out.append(f"| ... {len(rows) - max_rows} more rows |" + " |" * (len(headers) - 1))
    return "\n".join(out)


def _fmt(x, nd: int = 3) -> str:
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        if math.isinf(x):
            return "Inf" if x > 0 else "-Inf"
        if x == int(x) and abs(x) < 1e15:
            return f"{int(x):,}"
        return f"{x:.{nd}g}" if abs(x) < 1e-3 else f"{x:,.{nd}f}"
    if isinstance(x, int) and not isinstance(x, bool):
        return f"{x:,}"
    return str(x)


def rnum(v) -> str:
    """A number as data.table::fwrite writes it: 1000, 1.5, Inf, 15 significant digits."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(v):
        return "NA"
    if math.isinf(v):
        return "Inf" if v > 0 else "-Inf"
    if v == int(v) and abs(v) < 1e15:
        return str(int(v))
    return f"{v:.15g}"


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_jsonable(v) for v in o]
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if math.isnan(f) else ("Inf" if math.isinf(f) and f > 0 else "-Inf" if math.isinf(f) else f)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_summary(out: Path, summary: dict, report_lines: List[str]) -> None:
    (out / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2, default=str))
    (out / "report.md").write_text("\n".join(report_lines) + "\n")
    print(f"Report: {out / 'report.md'}")
    print(f"JSON: {out / 'summary.json'}")


_QUIET = {"record": True, "log": True}


def _record(sub: str, label: str, inputs=(), outputs=()) -> None:
    if not _QUIET["record"]:
        return
    try:
        import _localstore as ls  # type: ignore
        ls.record_analysis("principal_pseudobulks", subcommand=sub, label=label or "",
                           inputs=[str(x) for x in inputs], outputs=[str(x) for x in outputs])
    except Exception:  # noqa: BLE001 - provenance is best effort
        pass


# ---------------------------------------------------------------------------
# Chrom sizes, bgzip / tabix
# ---------------------------------------------------------------------------

def embedded_chrom_sizes() -> List[Tuple[str, int]]:
    out = []
    for tok in _CHROM_SIZES_EMBED.split():
        c, n = tok.rsplit(":", 1)
        out.append((c, int(n)))
    return out


def read_chrom_sizes(path: "Optional[str]") -> List[Tuple[str, int]]:
    """Chrom sizes in file order (the sort order); None -> the embedded IGVF DACC GRCh38 file."""
    if not path:
        return embedded_chrom_sizes()
    out = []
    with _opener(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                f = line.split()
                out.append((f[0], int(f[1]) if len(f) > 1 and f[1].isdigit() else 0))
    return out


def write_chrom_sizes(path: Path) -> Path:
    path.write_text("".join(f"{c}\t{n}\n" for c, n in embedded_chrom_sizes()))
    return path


def bgzip_and_index(plain: Path, out_gz: Path) -> Dict[str, Any]:
    """bgzip + tabix -p bed: pysam if installed, else bgzip/tabix on PATH, else plain gzip (no index)."""
    plain, out_gz = Path(plain), Path(out_gz)
    target_plain = out_gz.with_name(out_gz.name[:-3]) if out_gz.name.endswith(".gz") else out_gz
    if plain != target_plain:
        shutil.move(str(plain), str(target_plain))
    try:
        import pysam  # type: ignore
        res = pysam.tabix_index(str(target_plain), preset="bed", force=True, keep_original=False)
        return {"path": str(res), "index": str(res) + ".tbi", "method": "pysam"}
    except ImportError:
        pass
    if shutil.which("bgzip") and shutil.which("tabix"):
        import subprocess
        subprocess.run(["bgzip", "-f", str(target_plain)], check=True)
        subprocess.run(["tabix", "-f", "-p", "bed", str(out_gz)], check=True)
        return {"path": str(out_gz), "index": str(out_gz) + ".tbi", "method": "htslib"}
    with open(target_plain, "rb") as fi, gzip.open(out_gz, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    target_plain.unlink()
    print(f"pysam / htslib not available: {out_gz} is plain gzip without a tabix index. Would run: "
          f"bgzip {target_plain} && tabix -p bed {out_gz}")
    return {"path": str(out_gz), "index": None, "method": "gzip"}


# ---------------------------------------------------------------------------
# Directory discovery + QC guide (shared by both filter scripts)
# ---------------------------------------------------------------------------

def find_pseudobulk_dirs(pseudobulks_dir, cell_type: str, loose: bool = True) -> List[str]:
    """annotation-{ct}-IGVF* (sorted) per comma-separated ct, else annotation-{ct}; else annotation-{ct}-*."""
    all_dirs: List[str] = []
    for ct in str(cell_type).split(","):
        ct = ct.strip()
        if not ct:
            continue
        pattern = os.path.join(str(pseudobulks_dir), f"annotation-{glob.escape(ct)}-IGVF*")
        dirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))
        if dirs:
            log.info("found %d directories for %s", len(dirs), ct)
            all_dirs.extend(dirs)
            continue
        fallback = os.path.join(str(pseudobulks_dir), f"annotation-{ct}")
        if os.path.isdir(fallback):
            all_dirs.append(fallback)
            continue
        if loose:
            dirs = sorted(d for d in glob.glob(os.path.join(str(pseudobulks_dir), f"annotation-{glob.escape(ct)}-*"))
                          if os.path.isdir(d))
            if dirs:
                print(f"[info] no annotation-{ct}-IGVF* directories; using {len(dirs)} annotation-{ct}-* directories")
                all_dirs.extend(dirs)
                continue
        print(f"[warning] No directories found matching: {pattern}")
    return all_dirs


def load_passing_barcodes(qc_guide) -> set:
    """Column 1 of the QC guide (header skipped), full strings -- never truncated to the 16-bp sequence."""
    bcs = set()
    with _opener(qc_guide) as fh:
        fh.readline()
        for line in fh:
            line = line.strip()
            if line:
                bcs.add(line.split("\t")[0])
    return bcs


# ---------------------------------------------------------------------------
# Step 0: per-cell QC datatable
# ---------------------------------------------------------------------------

def build_qc_datatable(pseudobulks_dir, cell_type: str, out_path: Path) -> dict:
    dirs = find_pseudobulk_dirs(pseudobulks_dir, cell_type)
    if not dirs:
        raise SystemExit(f"[error] No annotation-{cell_type}[-IGVF*] directories found under {pseudobulks_dir}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header, n_rows, used, skipped = None, 0, [], []
    with open(out_path, "w") as out:
        for d in dirs:
            qp = os.path.join(d, "per_cell_qc.tsv")
            if not os.path.exists(qp):
                skipped.append((d, "no per_cell_qc.tsv"))
                continue
            if os.path.getsize(qp) == 0:
                skipped.append((d, "empty per_cell_qc.tsv"))
                continue
            with open(qp) as fh:
                h = fh.readline()
                if header is None:
                    header = h
                    out.write(h)
                elif h != header:
                    out.close()
                    out_path.unlink()
                    raise SystemExit(f"[error] Column mismatch in {qp}:\n  expected: {header.strip()}\n  got:      {h.strip()}")
                for line in fh:
                    out.write(line)
                    n_rows += 1
            used.append(d)
    if header is None:
        out_path.unlink()
        raise SystemExit(f"[error] No usable per_cell_qc.tsv files found for cell type '{cell_type}' under {pseudobulks_dir}")
    for d, why in skipped:
        print(f"[warn] {why} in {d}, skipping")
    return {"cells": n_rows, "directories": len(used), "skipped": skipped, "path": str(out_path)}


# ---------------------------------------------------------------------------
# Step 1: thresholds (explore_qc_thresholds.R + plot_per_cell_qc.R)
# ---------------------------------------------------------------------------

def parse_threshold_string(s: str, defaults: dict) -> dict:
    p = dict(defaults)
    toks = [t for t in str(s or "").strip().split() if t]
    i = 0
    while i < len(toks):
        flag = toks[i]
        if flag not in FLAG_MAP:
            raise SystemExit(f"Unknown flag: {flag}")
        if i + 1 >= len(toks):
            raise SystemExit(f"No value for: {flag}")
        try:
            val = float(toks[i + 1].replace("Inf", "inf"))
        except ValueError:
            raise SystemExit(f"Non-numeric value for {flag}: {toks[i + 1]}")
        p[FLAG_MAP[flag]] = val
        i += 2
    return p


def thresholds_from_args(args, defaults: dict) -> dict:
    p = dict(defaults)
    for k in THRESH_ORDER:
        v = getattr(args, k, None)
        if v is not None:
            p[k] = float(v)
    return p


def _differs(a, b) -> bool:
    a, b = float(a), float(b)
    if math.isinf(a) and math.isinf(b):
        return False
    return not math.isclose(a, b, rel_tol=1.5e-8, abs_tol=0.0)


def read_meta(path):
    """per_cell_qc datatable, subsample as text, stable-sorted by subsample, nUMI_non_MT added."""
    pd = _pd()
    meta = pd.read_csv(path, sep="\t", dtype={"subsample": str, "barcode": str, "analysis_accession": str})
    meta["subsample"] = meta["subsample"].fillna("NA").astype(str)
    for c in PER_CELL_QC_COLS[3:]:
        if c in meta.columns:
            meta[c] = pd.to_numeric(meta[c], errors="coerce")
    if "rna_read_count" in meta.columns and "pct_mito" in meta.columns:
        meta["nUMI_non_MT"] = meta["rna_read_count"] * (100 - meta["pct_mito"]) / 100
    return meta.sort_values("subsample", kind="mergesort").reset_index(drop=True)


def na_columns(meta) -> List[str]:
    return [c for c in PER_CELL_QC_COLS[3:] if c in meta.columns and meta[c].isna().all() and len(meta)]


def pass_mask(meta, p: dict, keys: Sequence[str], mode: str = "final", skip_cols: Sequence[str] = ()):
    """Cells passing thresholds `keys`; mode 'explore' (>= / <=) or 'final' (strict, plot_per_cell_qc.R). NA fails."""
    np = _np()
    ok = np.ones(len(meta), dtype=bool)
    for k in keys:
        col, op_e, op_f = THRESH[k]
        if col in skip_cols:
            continue
        op = op_e if mode == "explore" else op_f
        x = meta[col].to_numpy(dtype=float)
        v = float(p[k])
        with np.errstate(invalid="ignore"):
            r = {">=": x >= v, "<=": x <= v, ">": x > v, "<": x < v}[op]
        ok &= r & ~np.isnan(x)
    return ok


def failure_stats(meta, p: dict, skip_cols: Sequence[str] = ()):
    """build_failure_stats(): per threshold Total (cells failing it) and Alone (cells failing only it)."""
    pd, np = _pd(), _np()
    fail = pd.DataFrame({k: ~pass_mask(meta, p, [k], "explore", skip_cols) for k in THRESH_ORDER})
    n_fail = fail.sum(axis=1).to_numpy()
    total = {k: int(fail[k].sum()) for k in THRESH_ORDER}
    alone = {k: int((fail[k].to_numpy() & (n_fail == 1)).sum()) for k in THRESH_ORDER}
    return fail, total, alone


def label_threshold_set(p: dict, defaults: dict = EXPLORE_DEFAULTS) -> str:
    changed = [k for k in THRESH_ORDER if _differs(defaults[k], p[k])]
    if not changed:
        return "defaults"
    return ", ".join(f"{THRESHOLD_LABELS[k]}{rnum(p[k])}" for k in changed)


def explore_sets(meta, sets: Sequence[dict], show_subsamples: bool = False, skip_cols: Sequence[str] = ()):
    """print_report(): text report + a tidy table (set, threshold, total, alone) + per-subsample rows."""
    n_total = len(meta)
    subs = sorted(meta["subsample"].unique())
    lines = ["=" * 72, "QC THRESHOLD EXPLORATION REPORT", f"Input: {n_total} cells across {len(subs)} subsamples", "=" * 72, ""]
    rows, sub_rows, set_rows = [], [], []
    for si, p in enumerate(sets, 1):
        label = label_threshold_set(p)
        kept = pass_mask(meta, p, THRESH_ORDER, "explore", skip_cols)
        n_kept = int(kept.sum())
        n_drop = n_total - n_kept
        sub_kept = sum(1 for s in subs if kept[(meta["subsample"] == s).to_numpy()].any())
        lines += [f"-- Set {si}: {label}", "-" * 72,
                  f"  Overall:  {n_kept} kept  |  {n_drop} dropped ({100 * n_drop / n_total if n_total else 0:.1f}% of total)",
                  f"  Subsamples: {sub_kept} / {len(subs)} subsamples with >=1 cell passing", ""]
        _, total, alone = failure_stats(meta, p, skip_cols)
        overridden = [k for k in THRESH_ORDER if _differs(EXPLORE_DEFAULTS[k], p[k])]
        show = list(dict.fromkeys(ALWAYS_SHOW + overridden + [k for k in THRESH_ORDER if total[k] > 0]))
        top = max(show, key=lambda k: total[k]) if any(total[k] > 0 for k in show) else None
        # which.max picks the first maximum in `show` order
        if top is not None:
            best = max(total[k] for k in show)
            top = next(k for k in show if total[k] == best)
            lines.append(f"  Most stringent threshold (most cells dropped): {THRESHOLD_LABELS[top]}{rnum(p[top])}  [{total[top]} cells]")
        else:
            lines.append("  No cells dropped.")
        lines += ["", f"  {'Threshold':<26}  {'Total':>8}  {'Alone':>8}", f"  {'-' * 26}  {'-' * 8}  {'-' * 8}"]
        extra = [k for k in show if k not in ALWAYS_SHOW]
        extra = sorted(extra, key=lambda k: -total[k])
        order = [k for k in ALWAYS_SHOW if k in show] + extra
        for k in order:
            mark = " <" if k == top else ""
            lines.append(f"  {THRESHOLD_LABELS[k] + rnum(p[k]):<26}  {total[k]:>8}  {alone[k]:>8}{mark}")
            rows.append({"set": si, "set_label": label, "threshold": k, "label": THRESHOLD_LABELS[k] + rnum(p[k]),
                         "value": p[k], "total": total[k], "alone": alone[k], "most_stringent": k == top})
        set_rows.append({"set": si, "set_label": label, "n_cells": n_total, "kept": n_kept, "dropped": n_drop,
                         "pct_dropped": 100 * n_drop / n_total if n_total else 0.0, "subsamples_with_cells": sub_kept,
                         "n_subsamples": len(subs), "most_stringent": top or ""})
        if show_subsamples:
            lines += ["", "  Per-subsample breakdown:",
                      f"  {'Subsample':<24}  {'Total':>7}  {'Kept':>7}  {'Dropped':>7}  {'% Drop':>6}  Most stringent threshold"]
            for s in subs:
                sm = meta[meta["subsample"] == s].reset_index(drop=True)
                sk = int(pass_mask(sm, p, THRESH_ORDER, "explore", skip_cols).sum())
                _, st, _ = failure_stats(sm, p, skip_cols)
                active = [k for k in THRESH_ORDER if st[k] > 0]
                if active:
                    b = max(st[k] for k in active)
                    t = next(k for k in active if st[k] == b)
                    tl = f"{THRESHOLD_LABELS[t]}{rnum(p[t])} [{st[t]}]"
                else:
                    t, tl = "", "-"
                pct = 100 * (len(sm) - sk) / len(sm) if len(sm) else 0.0
                lines.append(f"  {s:<24}  {len(sm):>7}  {sk:>7}  {len(sm) - sk:>7}  {pct:>5.1f}%  {tl}")
                sub_rows.append({"set": si, "subsample": s, "total": len(sm), "kept": sk, "dropped": len(sm) - sk,
                                 "pct_dropped": pct, "most_stringent": t, "most_stringent_cells": st[t] if t else 0})
        lines.append("")
    lines.append("=" * 72)
    return "\n".join(lines), rows, set_rows, sub_rows


def apply_qc(meta, p: dict, skip_cols: Sequence[str] = ()):
    """plot_per_cell_qc.R filters: (rna mask, atac mask, all mask), strict inequalities."""
    rna = pass_mask(meta, p, RNA_KEYS, "final", skip_cols)
    atac = pass_mask(meta, p, ATAC_KEYS, "final", skip_cols)
    return rna, atac, rna & atac


def subsample_metrics(filt):
    pd = _pd()
    if not len(filt):
        return pd.DataFrame(columns=["subsample", "n_cells", "total_fragments", "total_RNA_reads", "mean_frag_per_cell",
                                     "mean_RNA_per_cell", "mean_frip", "mean_tss"])
    g = filt.groupby("subsample", sort=True)
    m = pd.DataFrame({"n_cells": g.size(), "total_fragments": g["num_frags"].sum(), "total_RNA_reads": g["rna_read_count"].sum(),
                      "mean_frip": g["frip"].mean(), "mean_tss": g["tss_enrichment"].mean()}).reset_index()
    m["mean_frag_per_cell"] = m["total_fragments"] / m["n_cells"]
    m["mean_RNA_per_cell"] = m["total_RNA_reads"] / m["n_cells"]
    return m[["subsample", "n_cells", "total_fragments", "total_RNA_reads", "mean_frag_per_cell", "mean_RNA_per_cell",
              "mean_frip", "mean_tss"]]


def write_r_tsv(df, path: Path, compress: bool = False) -> Path:
    """fwrite(sep='\\t', quote=FALSE): numbers in R format, NA as empty -> here 'NA' for clarity of R reads."""
    lines = ["\t".join(df.columns)]
    for row in df.itertuples(index=False):
        lines.append("\t".join(rnum(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else str(v) for v in row))
    text = "\n".join(lines) + "\n"
    if compress:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)
    return Path(path)


def run_qc_filter(meta_path, out: Path, p: dict, plots: bool = True, upstream_compat: bool = False) -> dict:
    pd = _pd()
    meta = read_meta(meta_path)
    skip = [] if upstream_compat else na_columns(meta)
    for c in skip:
        print(f"[warning] {c} is NA for every cell: its threshold is skipped (--upstream-compat drops every cell)")
    rna, atac, allm = apply_qc(meta, p, skip)
    filt = meta[allm].reset_index(drop=True)
    metrics = subsample_metrics(filt)
    out.mkdir(parents=True, exist_ok=True)
    files = {}
    files["metrics"] = write_r_tsv(metrics, out / "filtered_cell_subsample_metrics.tsv")
    guide = filt[GUIDE_COLS].copy()
    files["guide"] = write_r_tsv(guide, out / "filtered_barcodes_with_subsamples.tsv.gz", compress=True)
    th = pd.DataFrame({"threshold": [THRESHOLD_DOC[k] for k in THRESH_ORDER], "value": [p[k] for k in THRESH_ORDER]})
    files["thresholds"] = write_r_tsv(th, out / "qc_thresholds.tsv")
    for k in ("guide", "thresholds", "metrics"):
        print(f"TSV: {files[k]}")
    n, nk = len(meta), int(allm.sum())
    print(f"\nFiltering summary:\n  Total cells:   {n}\n  Cells kept:    {nk} ({100 * nk / n if n else 0:.1f}%)\n"
          f"  Cells dropped: {n - nk} ({100 * (n - nk) / n if n else 0:.1f}%)\n")
    figs = qc_figures(meta, p, rna, atac, allm, metrics, out) if plots else []
    return {"cells": n, "kept_rna": int(rna.sum()), "kept_atac": int(atac.sum()), "kept": nk, "skipped_thresholds": skip,
            "files": {k: str(v) for k, v in files.items()}, "figures": [str(f) for f in figs], "metrics": metrics,
            "thresholds": p, "guide": guide}


# ---------------------------------------------------------------------------
# QC figures (plot_per_cell_qc.R panels, PNG)
# ---------------------------------------------------------------------------

def _ax_style(ax, title, xlabel="", ylabel=""):
    from matplotlib.ticker import NullFormatter  # type: ignore
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "log":
            axis.set_minor_formatter(NullFormatter())
    ax.set_title(title, fontsize=8, loc="left", color=INK)
    ax.set_xlabel(xlabel, fontsize=7, color=INK2)
    ax.set_ylabel(ylabel, fontsize=7, color=INK2)
    ax.tick_params(labelsize=6, colors=INK2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)


def _vlines(ax, xs, horizontal=False):
    for x in xs:
        if x is None or not math.isfinite(float(x)) or float(x) <= -1e300:
            continue
        (ax.axhline if horizontal else ax.axvline)(x, color=GREYS[4], linestyle="--", linewidth=0.8)


def _density(ax, values, color, log=False, lines=(), lims=None):
    np = _np()
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if log:
        x = x[x > 0]
    if len(x) < 2 or np.ptp(x) == 0:
        ax.text(0.5, 0.5, "n < 2 or constant", transform=ax.transAxes, ha="center", fontsize=6, color=INK2)
        return
    from scipy.stats import gaussian_kde  # type: ignore
    t = np.log10(x) if log else x
    lo, hi = (np.log10(lims[0]), np.log10(lims[1])) if (log and lims) else (t.min(), t.max())
    grid = np.linspace(min(lo, t.min()), max(hi, t.max()), 256)
    try:
        d = gaussian_kde(t)(grid)
    except Exception:  # noqa: BLE001 - singular data
        return
    gx = 10 ** grid if log else grid
    ax.fill_between(gx, d, color=color, alpha=0.5, linewidth=0)
    ax.plot(gx, d, color=color, linewidth=1.2)
    if log:
        ax.set_xscale("log")
        if lims:
            ax.set_xlim(*lims)
    _vlines(ax, [v for v in lines if not (log and float(v) <= 0)])


def _bars(ax, counts, color, title):
    subs = list(counts.index)
    ax.bar(range(len(subs)), counts.values, color=color)
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(int(v)), ha="center", va="bottom", fontsize=6, rotation=45)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs, rotation=45, ha="right", fontsize=6)
    _ax_style(ax, title, "Subsample", "# cells")


def _box(ax, meta, col, color, title, ylabel, lines=(), log=False):
    np = _np()
    subs = list(dict.fromkeys(meta["subsample"]))
    data = []
    for s in subs:
        v = meta.loc[meta["subsample"] == s, col].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if log:
            v = np.log10(v[v > 0])
        data.append(v if len(v) else np.array([np.nan]))
    bp = ax.boxplot(data, patch_artist=True, flierprops={"marker": "o", "markersize": 2, "markerfacecolor": GREYS[2],
                                                         "markeredgecolor": GREYS[2]})
    for b in bp["boxes"]:
        b.set_facecolor(color)
    for m in bp["medians"]:
        m.set_color(INK)
    ax.set_xticks(range(1, len(subs) + 1))
    ax.set_xticklabels(subs, rotation=45, ha="right", fontsize=6)
    lv = [math.log10(x) for x in lines if math.isfinite(float(x)) and float(x) > 0] if log else list(lines)
    _vlines(ax, lv, horizontal=True)
    _ax_style(ax, title, "Subsample", ylabel)


def _scatter(fig, ax, meta, x, y, c, title, xl, yl, clabel, vx=(), hy=(), logy=True):
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
    cmap = LinearSegmentedColormap.from_list("pp", ["#d3a9ce", "#430b4e"])
    m = meta.sort_values(c)
    sc = ax.scatter(m[x], m[y], c=m[c], cmap=cmap, s=3, linewidths=0)
    ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    _vlines(ax, vx)
    _vlines(ax, hy, horizontal=True)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(clabel, fontsize=6)
    _ax_style(ax, title, xl, yl)


def qc_figures(meta, p, rna, atac, allm, metrics, out: Path) -> List[Path]:
    plt = _plt()
    if plt is None:
        print("matplotlib not installed: figures skipped")
        return []
    np = _np()
    figs = []
    counts = meta.groupby("subsample", sort=True).size()
    r, a = CP[9], CP[6]
    # RNA_QC_plots
    fig = plt.figure(figsize=(12, 13))
    gs = fig.add_gridspec(4, 3, hspace=0.75, wspace=0.4)
    _bars(fig.add_subplot(gs[0, 0]), counts, GREYS[2], "Cells per subsample")
    ax = fig.add_subplot(gs[0, 1]); _density(ax, meta["rna_read_count"], r, True, (p["rna_min"], p["rna_max"]), (100, 70000))
    _ax_style(ax, "RNA read count density", "# RNA reads", "Cell density")
    ax = fig.add_subplot(gs[0, 2]); _density(ax, meta["nUMI_non_MT"], r, True, (), (100, 70000))
    _ax_style(ax, "Non-mito read count density", "# non-mito RNA reads", "Cell density")
    ax = fig.add_subplot(gs[1, 0]); _density(ax, meta["gene_count"], r, True, (p["gene_min"], p["gene_max"]), (100, 25000))
    _ax_style(ax, "# genes/cell distribution", "# genes per cell", "Cell density")
    _box(fig.add_subplot(gs[1, 1:]), meta, "gene_count", r, "log10(# genes/cell) per subsample", "log10(# genes per cell)",
         (p["gene_min"], p["gene_max"]), log=True)
    _scatter(fig, fig.add_subplot(gs[2, 0]), meta, "rna_read_count", "gene_count", "pct_mito", "Genes x reads, colored by % mito",
             "# RNA reads per cell", "# genes per cell", "% mito", (p["rna_min"], p["rna_max"]), (p["gene_min"], p["gene_max"]))
    _scatter(fig, fig.add_subplot(gs[2, 1]), meta, "rna_read_count", "gene_count", "pct_ribo", "Genes x reads, colored by % ribo",
             "# RNA reads per cell", "# genes per cell", "% ribo", (p["rna_min"], p["rna_max"]), (p["gene_min"], p["gene_max"]))
    _scatter(fig, fig.add_subplot(gs[2, 2]), meta, "rna_read_count", "pct_mito", "pct_ribo", "Reads x % mito, colored by % ribo",
             "# RNA reads per cell", "% mitochondrial reads", "% ribo", (p["rna_min"], p["rna_max"]), (p["pct_mt_max"],), logy=False)
    ax = fig.add_subplot(gs[3, 0]); _density(ax, meta["pct_mito"], r, True, (p["pct_mt_max"],))
    _ax_style(ax, "% mito density", "% mitochondrial reads", "Cell density")
    ax = fig.add_subplot(gs[3, 1]); _density(ax, meta["pct_ribo"], r, True, (p["pct_ribo_max"],))
    _ax_style(ax, "% ribo density", "% ribosomal reads", "Cell density")
    with np.errstate(divide="ignore", invalid="ignore"):
        gpr = np.log10(meta["gene_count"].to_numpy(float) / meta["rna_read_count"].to_numpy(float))
    ax = fig.add_subplot(gs[3, 2]); _density(ax, gpr, r)
    _ax_style(ax, "log10(genes per read)", "log10(# genes / # reads)", "Cell density")
    figs.append(_savefig(fig, out / "RNA_QC_plots.png"))
    # ATAC_QC_plots
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.7, wspace=0.4)
    _bars(fig.add_subplot(gs[0, 0]), counts, GREYS[2], "Cells per subsample")
    ax = fig.add_subplot(gs[0, 1]); _density(ax, meta["num_frags"], a, False, (p["atac_min"], p["atac_max"]))
    _ax_style(ax, "Fragment density", "# fragments", "Cell density")
    _box(fig.add_subplot(gs[0, 2]), meta, "num_frags", a, "log10(# fragments/cell) per subsample", "log10(# fragments per cell)",
         (p["atac_min"], p["atac_max"]), log=True)
    ax = fig.add_subplot(gs[1, 0]); _density(ax, meta["nucleosomal_signal"], a, False, (p["nuc_signal_max"],))
    _ax_style(ax, "Nucleosomal signal density", "Nucleosomal signal", "Cell density")
    _box(fig.add_subplot(gs[1, 1]), meta, "tss_enrichment", a, "TSS enrichment per subsample", "TSS enrichment", (p["tss_enr_min"],))
    ax = fig.add_subplot(gs[1, 2]); _density2d(fig, ax, meta, p)
    ax = fig.add_subplot(gs[2, 0]); _density(ax, meta["pct_duplicated_reads"], a, False, (p["pct_dup_max"],))
    _ax_style(ax, "% duplicated reads density", "% duplicated ATAC reads", "Cell density")
    ax = fig.add_subplot(gs[2, 1]); _density(ax, meta["frip"], a, False, (p["frip_min"],))
    _ax_style(ax, "FRIP density", "Fraction of reads in peaks (FRIP)", "Cell density")
    _box(fig.add_subplot(gs[2, 2]), meta, "frip", a, "FRIP per subsample", "FRIP", (p["frip_min"],))
    figs.append(_savefig(fig, out / "ATAC_QC_plots.png"))
    for mask, name, title, col in ((rna, "cells_per_subsample_after_RNA_QC.png", "Cells per subsample\nafter RNA QC filtering", GREYS[2]),
                                   (atac, "cells_per_subsample_after_ATAC_QC.png", "Cells per subsample\nafter ATAC QC filtering", GREYS[2]),
                                   (allm, "cells_per_subsample_after_all_QC.png", "Cells per subsample\nafter RNA+ATAC QC filtering", GREYS[3])):
        fig, ax = plt.subplots(figsize=(4, 3.5))
        _bars(ax, meta[mask].groupby("subsample", sort=True).size(), col, title)
        figs.append(_savefig(fig, out / name))
    return figs


def _density2d(fig, ax, meta, p):
    np = _np()
    x = meta["num_frags"].to_numpy(float)
    y = meta["tss_enrichment"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    _ax_style(ax, "Fragments x TSS enrichment", "# fragments", "TSS enrichment")
    if ok.sum() < 3 or np.ptp(y[ok]) == 0:
        return
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
    from scipy.stats import gaussian_kde  # type: ignore
    lx = np.log10(x[ok])
    try:
        k = gaussian_kde(np.vstack([lx, y[ok]]))
    except Exception:  # noqa: BLE001
        return
    gx, gy = np.meshgrid(np.linspace(lx.min(), lx.max(), 80), np.linspace(y[ok].min(), y[ok].max(), 80))
    z = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    cmap = LinearSegmentedColormap.from_list("tss", ["#ffffff", "#c5e5fb", "#002359"], N=8)
    cs = ax.contourf(10 ** gx, gy, z, levels=8, cmap=cmap)
    ax.set_xscale("log")
    _vlines(ax, (p["atac_min"], p["atac_max"]))
    _vlines(ax, (p["tss_enr_min"],), horizontal=True)
    fig.colorbar(cs, ax=ax, fraction=0.046, pad=0.02).ax.tick_params(labelsize=5)


def _savefig(fig, path: Path) -> Path:
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor="white")
    _plt().close(fig)
    print(f"Figure: {path}")
    return path


# ---------------------------------------------------------------------------
# Step 2a: ATAC fragments
# ---------------------------------------------------------------------------

def _stream_filter_fragments(frag_path, passing: set, out_fh, found: set, clean: bool = False) -> Tuple[int, int]:
    n_pass = n_total = 0
    with _opener(frag_path) as fh:
        for ln, line in enumerate(fh, 1):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 5:
                raise ValueError(f"Malformed fragment record in {frag_path}, line {ln}: expected 5 tab-separated fields "
                                 f"(chrom, start, end, barcode, duplicate count), got {len(f)}: {line.rstrip()!r}")
            n_total += 1
            bc = f[3]
            if bc in passing:
                found.add(bc)
                if clean:
                    f[3] = bc.split("_", 1)[0]
                out_fh.write("\t".join(f) + "\n")
                n_pass += 1
    return n_pass, n_total


def sort_fragment_chunks(split_dir: Path, chrom_order: Sequence[str], out_plain: Path) -> Dict[str, int]:
    """Concatenate per-chromosome files in chrom-sizes order, each sorted like sort-bed (start, end, rest)."""
    pd = _pd()
    written = {}
    with open(out_plain, "w") as out:
        for c in chrom_order:
            f = split_dir / f"{c}.bed"
            if not f.exists():
                continue
            df = pd.read_csv(f, sep="\t", header=None, dtype={0: str, 3: str}, keep_default_na=False)
            df = df.sort_values([1, 2, 3, 4], kind="mergesort")
            df.to_csv(out, sep="\t", header=False, index=False)
            written[c] = len(df)
    return written


def filter_atac(qc_guide, pseudobulks, cell_type: str, chrom_sizes: "Optional[str]", out_path: Path,
                clean: bool = False) -> dict:
    passing = load_passing_barcodes(qc_guide)
    dirs = find_pseudobulk_dirs(pseudobulks, cell_type)
    if not dirs:
        raise SystemExit("[error] No fragment directories found.")
    out_path = Path(out_path)
    if not str(out_path).endswith(".gz"):
        out_path = Path(str(out_path) + ".gz")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    found: set = set()
    per_dir = []
    tmpd = Path(tempfile.mkdtemp(prefix="pp_atac_", dir=str(out_path.parent)))
    try:
        raw = tmpd / "passing.tsv"
        with open(raw, "w") as tmp:
            for d in dirs:
                fp = os.path.join(d, "fragments.tsv.gz")
                if not os.path.isfile(fp):
                    print(f"[warning] fragments.tsv.gz not found in {d}, skipping.")
                    continue
                npass, ntot = _stream_filter_fragments(fp, passing, tmp, found, clean)
                per_dir.append({"dir": os.path.basename(d), "retained": npass, "total": ntot})
        missing = passing - found
        if missing:
            ex = "\n".join(f"  {b}" for b in sorted(missing)[:20])
            raise SystemExit(f"[error] {len(missing)} barcode(s) from the QC guide were not found in any fragment file:\n{ex}")
        split = tmpd / "split"
        split.mkdir()
        handles: Dict[str, Any] = {}
        with open(raw) as fh:
            for line in fh:
                c = line.split("\t", 1)[0]
                h = handles.get(c)
                if h is None:
                    h = handles[c] = open(split / f"{c}.bed", "w")
                h.write(line)
        for h in handles.values():
            h.close()
        order = [c for c, _ in read_chrom_sizes(chrom_sizes)]
        dropped = {}
        for c in sorted(set(handles) - set(order)):
            n = sum(1 for _ in open(split / f"{c}.bed"))
            dropped[c] = n
            print(f"[warning] Chrom '{c}' not in chrom sizes file - skipping ({n} rows dropped).")
        plain = tmpd / "sorted.tsv"
        written = sort_fragment_chunks(split, order, plain)
        idx = bgzip_and_index(plain, out_path)
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    tot_pass = sum(x["retained"] for x in per_dir)
    tot = sum(x["total"] for x in per_dir)
    return {"out": idx["path"], "index": idx["index"], "index_method": idx["method"], "directories": per_dir,
            "rows_retained": tot_pass, "rows_total": tot, "rows_written": sum(written.values()),
            "chromosomes_dropped": dropped, "guide_barcodes": len(passing), "barcodes_found": len(found)}


# ---------------------------------------------------------------------------
# Step 2b: RNA count matrix
# ---------------------------------------------------------------------------

_GID = re.compile(r'gene_id "([^"]+)"')
_GNAME = re.compile(r'gene_name "([^"]+)"')


def parse_gtf(gtf_path) -> Tuple[Dict[str, Tuple[str, str]], List[tuple]]:
    """versioned gene_id -> (chrom, gene_name), first occurrence wins; conflicting names recorded."""
    m: Dict[str, Tuple[str, str]] = {}
    conflicts = []
    with _opener(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            gi, gn = _GID.search(parts[8]), _GNAME.search(parts[8])
            if gi and gn:
                vid, name = gi.group(1), gn.group(1)
                if vid in m:
                    if m[vid][1] != name:
                        conflicts.append((vid, m[vid][0], m[vid][1], parts[0], name))
                else:
                    m[vid] = (parts[0], name)
    return m, conflicts


def _collapse(adata, symbols):
    """Sum columns sharing a symbol: X @ S with S (n_ids x n_symbols); symbols sorted (np.unique)."""
    np, sp, pd, ad = _np(), _sp(), _pd(), _ad()
    uniq, inv = np.unique(np.asarray(symbols, dtype=str), return_inverse=True)
    S = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uniq)))

    def col(mat):
        mat = mat if sp.issparse(mat) else sp.csr_matrix(mat)
        return sp.csr_matrix(mat @ S)
    var = pd.DataFrame(index=pd.Index(uniq, name="gene_symbol"))
    return ad.AnnData(X=col(adata.X), obs=adata.obs.copy(), var=var, layers={k: col(adata.layers[k]) for k in adata.layers}), \
        int((np.bincount(inv) > 1).sum())


def collapse_with_gtf(adata, gtf_map, standard_only: bool, log_lines: List[str]):
    np = _np()
    ids = list(map(str, adata.var_names))
    unmatched = [e for e in ids if e not in gtf_map]
    log_lines += [f"Total Ensembl IDs in anndata: {len(ids)}", f"Matched to GTF (exact versioned): {len(ids) - len(unmatched)}",
                  f"Unmatched: {len(unmatched)}"]
    if unmatched:
        log_lines.append(f"\nUnmatched Ensembl IDs ({len(unmatched)}):")
        log_lines += [f"  {e}" for e in sorted(unmatched)]
        return None, {"unmatched": len(unmatched)}
    info: Dict[str, Any] = {"unmatched": 0}
    if standard_only:
        nonstd = [e for e in ids if gtf_map[e][0] not in STANDARD_CHROMS]
        log_lines.append("\n--- Standard chromosomes filter ---")
        std_sym = {gtf_map[e][1] for e in ids if gtf_map[e][0] in STANDARD_CHROMS}
        non_sym = {gtf_map[e][1] for e in nonstd}
        shared, lost = std_sym & non_sym, non_sym - std_sym
        if nonstd:
            log_lines.append(f"Ensembl IDs on nonstandard chromosomes: {len(nonstd)}")
            for c, n in sorted(Counter(gtf_map[e][0] for e in nonstd).items()):
                log_lines.append(f"  {c}: {n} Ensembl IDs")
            log_lines.append(f"\nGene symbols with contributions from BOTH standard and nonstandard chromosomes ({len(shared)}):")
            for s in sorted(shared):
                log_lines.append(f"  {s}: " + ", ".join(f"{e} ({gtf_map[e][0]})" for e in ids if gtf_map[e][1] == s))
            if lost:
                log_lines.append(f"\nGene symbols ONLY on nonstandard chromosomes (removed entirely): {len(lost)}")
                for s in sorted(lost):
                    log_lines.append(f"  {s}: " + ", ".join(e for e in nonstd if gtf_map[e][1] == s))
        else:
            log_lines.append("All Ensembl IDs are on standard chromosomes. No filtering needed.")
        keep = np.array([gtf_map[e][0] in STANDARD_CHROMS for e in ids])
        adata = adata[:, keep].copy()
        ids = list(map(str, adata.var_names))
        log_lines.append(f"Ensembl IDs after chromosome filtering: {len(ids)}")
        info.update({"nonstandard_ids": len(nonstd), "shared_symbols": sorted(shared), "lost_symbols": sorted(lost)})
    syms = [gtf_map[e][1] for e in ids]
    new, n_over = _collapse(adata, syms)
    log_lines += ["\n--- Gene symbol collapse ---", f"Ensembl IDs entering collapse: {len(ids)}",
                  f"Unique gene symbols after collapse: {new.n_vars}", f"Symbols with multiple Ensembl IDs (summed): {n_over}"]
    if n_over:
        by = defaultdict(list)
        for e, s in zip(ids, syms):
            by[s].append(e)
        log_lines.append("\nOverloaded symbols (multiple Ensembl IDs -> same symbol):")
        for s in sorted(k for k, v in by.items() if len(v) > 1):
            log_lines.append(f"  {s} ({len(by[s])} IDs): " + ", ".join(f"{e} ({gtf_map[e][0]})" for e in by[s]))
    info.update({"ids_in": len(ids), "symbols": new.n_vars, "overloaded": n_over})
    return new, info


def rna_sanity_checks(adata, passing: set, ensembl_ids: bool) -> List[str]:
    pd = _pd()
    errors = []
    oc = pd.Series(list(adata.obs_names)).value_counts()
    if (oc > 1).any():
        errors.append(f"Duplicate cell barcodes ({int((oc > 1).sum())} offending barcodes)")
    vc = pd.Series(list(adata.var_names)).value_counts()
    if (vc > 1).any():
        errors.append(f"Duplicate gene identifiers ({int((vc > 1).sum())} offending names)")
    if adata.n_obs != len(passing):
        found = set(adata.obs_names)
        errors.append(f"Cell count mismatch: matrix has {adata.n_obs} cells, QC guide has {len(passing)} "
                      f"(difference: {adata.n_obs - len(passing):+d}); missing {len(passing - found)}, extra {len(found - passing)}")
    return errors


def write_mtx_dir(adata, out_dir: Path, integer: bool = True) -> Dict[str, str]:
    import scipy.io  # type: ignore
    np, sp = _np(), _sp()
    out_dir.mkdir(parents=True, exist_ok=True)
    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    M = sp.coo_matrix(X.T)
    if integer and M.nnz and np.allclose(M.data, np.round(M.data)):
        M = sp.coo_matrix((np.round(M.data).astype(np.int64), (M.row, M.col)), shape=M.shape)
    elif integer and not M.nnz:
        M = M.astype(np.int64)
    paths = {"matrix": out_dir / "matrix.mtx.gz", "barcodes": out_dir / "barcodes.tsv.gz", "features": out_dir / "features.tsv.gz"}
    with gzip.open(paths["matrix"], "wb") as f:
        scipy.io.mmwrite(f, M, field="integer" if M.dtype.kind in "iu" else "real")
    with gzip.open(paths["barcodes"], "wt") as f:
        f.write("".join(f"{b}\n" for b in adata.obs_names))
    with gzip.open(paths["features"], "wt") as f:
        f.write("".join(f"{g}\n" for g in adata.var_names))
    return {k: str(v) for k, v in paths.items()}


def detect_format(out_path: str) -> str:
    n = out_path.lower()
    if n.endswith(".csv.gz") or n.endswith(".csv"):
        return "csv"
    if n.endswith(".h5ad") or n.endswith(".h5"):
        return "h5ad"
    if n.endswith(".mtx.gz") or n.endswith(".mtx"):
        return "mtx"
    raise SystemExit(f"[error] Cannot infer output format from '{out_path}'. Supported: .csv, .csv.gz, .h5ad, .h5, .mtx, .mtx.gz")


def _strip_ext(p: str) -> str:
    for ext in (".mtx.gz", ".mtx", ".csv.gz", ".csv", ".h5ad", ".h5"):
        if p.lower().endswith(ext):
            return p[: -len(ext)]
    return p


def filter_rna(qc_guide, pseudobulks, cell_type: str, out_path: str, gtf: "Optional[str]" = None,
               ensembl_ids: bool = False, standard_only: bool = False, log_path: "Optional[str]" = None,
               integer: bool = True) -> dict:
    ad = _ad()
    if standard_only and not gtf:
        raise SystemExit("[error] --standard-chromosomes-only requires --gtf.")
    if gtf and ensembl_ids:
        raise SystemExit("[error] --gtf cannot be used with --ensemblIDs-as-genes.")
    fmt = detect_format(out_path)
    passing = load_passing_barcodes(qc_guide)
    dirs = find_pseudobulk_dirs(pseudobulks, cell_type)
    if not dirs:
        raise SystemExit("[error] No matching directories found.")
    parts, per_dir = [], []
    for d in dirs:
        h = os.path.join(d, "rna_counts_mtx.h5ad")
        if not os.path.isfile(h):
            print(f"[warning] rna_counts_mtx.h5ad not found in {d}, skipping.")
            continue
        a = ad.read_h5ad(h)
        n0 = a.n_obs
        a = a[a.obs_names.isin(list(passing))].copy()
        per_dir.append({"dir": os.path.basename(d), "cells": n0, "passed": a.n_obs})
        if a.n_obs:
            parts.append(a)
    if not parts:
        raise SystemExit("[error] No passing cells found across any directory.")
    comb = ad.concat(parts, axis=0, join="inner", merge="same")
    info: Dict[str, Any] = {"directories": per_dir, "ensembl_ids_in": int(comb.n_vars)}
    written_log = None
    if not ensembl_ids:
        if gtf:
            lines = ["=== GTF-based gene symbol mapping log ===", f"GTF file: {gtf}", f"Standard chromosomes only: {standard_only}", ""]
            gmap, conflicts = parse_gtf(gtf)
            lines.append(f"Unique versioned Ensembl gene IDs in GTF: {len(gmap)}")
            if conflicts:
                lines.append(f"\nGTF conflicts (same Ensembl ID, different gene names): {len(conflicts)}")
                lines += [f"  {v}: '{pn}' ({pc}) vs '{nn}' ({nc})" for v, pc, pn, nc, nn in conflicts]
            lines.append("")
            new, ginfo = collapse_with_gtf(comb, gmap, standard_only, lines)
            written_log = Path(log_path) if log_path else Path(_strip_ext(out_path) + "_gtf_mapping.txt")
            written_log.parent.mkdir(parents=True, exist_ok=True)
            written_log.write_text("\n".join(lines) + "\n")
            info.update(ginfo)
            info["gtf_conflicts"] = len(conflicts)
            if new is None:
                raise SystemExit(f"[error] {ginfo['unmatched']} Ensembl IDs in anndata not found in GTF (exact versioned match). "
                                 f"See log: {written_log}")
            comb = new
        else:
            if "gene_symbol" not in comb.var.columns:
                raise SystemExit("[error] no --gtf and no var['gene_symbol'] column: pass --gtf or --ensemblIDs-as-genes")
            comb, n_over = _collapse(comb, comb.var["gene_symbol"].astype(str).values)
            info.update({"symbols": int(comb.n_vars), "overloaded": n_over})
    errs = rna_sanity_checks(comb, passing, ensembl_ids)
    if errs:
        raise SystemExit("[error] Sanity check(s) failed:\n  " + "\n  ".join(errs))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    if fmt == "csv":
        pd, np, sp = _pd(), _np(), _sp()
        op = out_path if out_path.endswith(".gz") else out_path + ".gz"
        X = comb.X.toarray() if sp.issparse(comb.X) else np.asarray(comb.X)
        df = pd.DataFrame(X, index=comb.obs_names, columns=comb.var_names)
        df.index.name = "barcode"
        df.to_csv(op, compression="gzip")
        outputs = {"csv": op}
    elif fmt == "h5ad":
        comb.write_h5ad(out_path)
        outputs = {"h5ad": out_path}
    else:
        outputs = write_mtx_dir(comb, Path(_strip_ext(out_path)), integer=integer)
    info.update({"cells": int(comb.n_obs), "genes": int(comb.n_vars), "format": fmt, "outputs": outputs,
                 "log": str(written_log) if written_log else None})
    return info


def package_rna(mtx_dir: Path, out_tar: "Optional[Path]" = None) -> Path:
    """The distributed .tar.gz: matrix.mtx, barcodes.tsv, features.tsv (members decompressed)."""
    mtx_dir = Path(mtx_dir)
    out_tar = Path(out_tar) if out_tar else mtx_dir.with_name(mtx_dir.name + ".tar.gz")
    with tarfile.open(out_tar, "w:gz") as tar:
        for name in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
            src = mtx_dir / (name + ".gz")
            data = gzip.decompress(src.read_bytes()) if src.exists() else (mtx_dir / name).read_bytes()
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            ti.mtime = int(time.time())
            tar.addfile(ti, io.BytesIO(data))
    return out_tar


# ---------------------------------------------------------------------------
# Config table (rule make_config_table) + Snakefile driver
# ---------------------------------------------------------------------------

def config_table(out_root: Path, dataset: str, cell_types: Sequence[str], upstream_compat: bool = False):
    pd = _pd()
    rows = [{"cluster": ct, "rna_matrix_file": str(Path(out_root) / dataset / ct / f"rna_count_matrix_{dataset}_{ct}"),
             "atac_frag_file": str(Path(out_root) / dataset / ct / f"atac_fragments_{dataset}_{ct}.tsv.gz"),
             "HiC_file": "", "HiC_type": "", "HiC_resolution": "", "alt_TSS": "", "alt_genes": "",
             "model_dir": UPSTREAM_MODEL_DIR if upstream_compat else SCE2G_MODEL_DIR} for ct in cell_types]
    return pd.DataFrame(rows, columns=CONFIG_TABLE_COLS)


def determine_mem_mb(input_size_mb: float, gz: bool, attempt: int = 1, min_gb: int = 8, max_mb: int = 250_000) -> int:
    size = input_size_mb * (8 if gz else 1)
    return int(min((2 ** (attempt - 1)) * max(4 * size, min_gb * 1000), max_mb))


def _dir_size_mb(paths) -> float:
    tot = 0
    for p in paths:
        for root, _, files in os.walk(p):
            tot += sum(os.path.getsize(os.path.join(root, f)) for f in files)
    return tot / 1e6


def parse_simple_yaml(text: str) -> dict:
    """The config_QC_pseudobulks.yaml subset: `key: value` scalars and `datasets:` -> {ds: [ct, ...]}."""
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text) or {}
    except ImportError:
        pass
    cfg: Dict[str, Any] = {}
    cur_key = cur_ds = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].rstrip() if not raw.strip().startswith("#") else ""
        if not line.strip():
            continue
        ind = len(line) - len(line.lstrip())
        s = line.strip()
        if ind == 0:
            k, _, v = s.partition(":")
            v = v.strip().strip('"').strip("'")
            cur_key, cur_ds = k.strip(), None
            cfg[cur_key] = v if v else {}
        elif s.startswith("- "):
            item = s[2:].strip().strip('"').strip("'")
            if cur_ds is not None:
                cfg[cur_key][cur_ds].append(item)
            else:
                if not isinstance(cfg.get(cur_key), list):
                    cfg[cur_key] = []
                cfg[cur_key].append(item)
        else:
            k, _, v = s.partition(":")
            v = v.strip()
            if v.startswith("[") and v.endswith("]"):
                cfg[cur_key][k.strip()] = [x.strip().strip('"').strip("'") for x in v[1:-1].split(",") if x.strip()]
                cur_ds = None
            else:
                cur_ds = k.strip()
                cfg[cur_key][cur_ds] = []
    return cfg


def load_run_config(path: str) -> dict:
    text = Path(path).read_text()
    if path.endswith(".json"):
        return json.loads(text)
    return parse_simple_yaml(text)


# ---------------------------------------------------------------------------
# Building primary pseudobulks from Portal files
# ---------------------------------------------------------------------------

def read_table_any(path):
    pd = _pd()
    p = str(path)
    if p.endswith((".h5ad", ".h5")):
        a = _ad().read_h5ad(p, backed="r")
        df = a.obs.copy()
        df.insert(0, "barcode", list(map(str, a.obs_names)))
        return df.reset_index(drop=True)
    sep = "," if p.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)


def _pick_col(df, wanted: "Optional[str]", candidates: Sequence[str], what: str, required: bool = True):
    if wanted:
        if wanted not in df.columns:
            raise SystemExit(f"{what} column {wanted!r} not in {list(df.columns)}")
        return wanted
    low = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    if required:
        raise SystemExit(f"no {what} column found (tried {list(candidates)}); pass it explicitly")
    return None


def safe_cell_type(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(s).strip()).strip("_") or "NA"


def read_barcode_map(path) -> Dict[str, str]:
    """ATAC barcode -> RNA barcode (two columns; 10x multiome ATAC/GEX whitelist pairs, -1 suffix ignored)."""
    m = {}
    with _opener(path) as fh:
        for line in fh:
            f = line.strip().split()
            if len(f) >= 2 and not line.startswith("#"):
                m[GEM_SUFFIX_RE.sub("", f[0])] = GEM_SUFFIX_RE.sub("", f[1])
    return m


def full_barcode(raw: str, lane: str, atac_map: "Optional[Dict[str, str]]" = None) -> str:
    """{16-bp}_{lane}: strip the GEM suffix (-1), translate ATAC -> RNA barcode, append the lane accession
    unless an _suffix is already present (a subpool / lane suffix written by the pipeline is kept)."""
    b = GEM_SUFFIX_RE.sub("", str(raw))
    if "_" in b:
        seq, suf = b.split("_", 1)
        if atac_map:
            seq = atac_map.get(seq, seq)
        return f"{seq}_{suf}"
    if atac_map:
        b = atac_map.get(b, b)
    return f"{b}_{lane}" if lane else b


def gene_symbols_of(adata) -> List[str]:
    for c in ("gene_symbol", "gene_name", "gene_symbols", "feature_name", "gene_names", "symbol"):
        if c in adata.var.columns:
            return list(map(str, adata.var[c].values))
    return list(map(str, adata.var_names))


def rna_cell_qc(adata) -> "Any":
    """rna_read_count (UMIs), gene_count (genes > 0), pct_mito (^MT-), pct_ribo (^RP[SL]) per cell."""
    pd, np, sp = _pd(), _np(), _sp()
    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(np.asarray(adata.X))
    X = sp.csr_matrix(X)
    sym = gene_symbols_of(adata)
    mt = np.array([s.upper().startswith("MT-") for s in sym])
    rb = np.array([bool(re.match(r"^RP[SL]", s.upper())) for s in sym])
    tot = np.asarray(X.sum(axis=1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        pm = np.where(tot > 0, 100 * np.asarray(X[:, mt].sum(axis=1)).ravel() / tot, np.nan)
        pr = np.where(tot > 0, 100 * np.asarray(X[:, rb].sum(axis=1)).ravel() / tot, np.nan)
    X.eliminate_zeros()
    return pd.DataFrame({"rna_read_count": tot, "gene_count": np.diff(X.indptr),
                         "pct_mito": pm, "pct_ribo": pr}, index=list(map(str, adata.obs_names)))


def read_tss(path) -> Dict[str, Tuple[Any, Any]]:
    """TSS bed (chr, start, end[, name, score, strand]) -> per chrom (sorted positions, strand sign).
    The TSS is start for +, end-1 for - when the interval is wider than 1 bp; mid-point otherwise."""
    pd, np = _pd(), _np()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    strand = df[5].astype(str).to_numpy() if df.shape[1] > 5 else np.array(["+"] * len(df))
    s, e = df[1].to_numpy(np.int64), df[2].to_numpy(np.int64)
    pos = np.where(e - s <= 1, s, np.where(strand == "-", e - 1, s))
    if df.shape[1] > 5 and ((e - s) > 2).all() and ((e - s) % 2 == 0).all() and (e - s).min() == (e - s).max():
        pos = (s + e) // 2  # fixed-width windows (e.g. TSS500) centred on the TSS
    out = {}
    for c in pd.unique(df[0]):
        m = (df[0] == c).to_numpy()
        o = np.argsort(pos[m], kind="mergesort")
        out[c] = (pos[m][o], np.where(strand[m][o] == "-", -1, 1))
    return out


def read_peaks_merged(path) -> Dict[str, Tuple[Any, Any]]:
    pd, np = _pd(), _np()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2], dtype={0: str})
    out = {}
    for c, g in df.groupby(0, sort=False):
        s = np.sort(g[1].to_numpy(np.int64))
        e = g[2].to_numpy(np.int64)[np.argsort(g[1].to_numpy(np.int64), kind="mergesort")]
        ms, me = [s[0]], [e[0]]
        for a, b in zip(s[1:], e[1:]):
            if a <= me[-1]:
                me[-1] = max(me[-1], b)
            else:
                ms.append(a)
                me.append(b)
        out[c] = (np.array(ms), np.array(me))
    return out


def atac_cell_counts(frag, cell_index, tss=None, peaks=None):
    """Per-cell sums for num_frags, reads (col5), nucleosome-free / mono, TSS centre / flank, fragments in peaks."""
    np = _np()
    n = len(cell_index)
    ci = frag["cell"].to_numpy()
    s = frag["start"].to_numpy(np.int64)
    e = frag["end"].to_numpy(np.int64)
    ln = e - s
    res = {"num_frags": np.bincount(ci, minlength=n).astype(float),
           "reads": np.bincount(ci, weights=frag["count"].to_numpy(float), minlength=n),
           "nf": np.bincount(ci, weights=(ln < 147).astype(float), minlength=n),
           "mono": np.bincount(ci, weights=((ln >= 147) & (ln < 294)).astype(float), minlength=n),
           "center": np.zeros(n), "flank": np.zeros(n), "in_peaks": np.zeros(n)}
    chroms = frag["chr"].to_numpy()
    for c in np.unique(chroms):
        m = chroms == c
        if tss is not None and c in tss:
            tp, tsg = tss[c]
            for ins in (s[m], e[m] - 1):
                cc = ci[m]
                lo = np.searchsorted(tp, ins - 1000, "left")
                hi = np.searchsorted(tp, ins + 1000, "right")
                k = hi - lo
                if not k.sum():
                    continue
                rep_ins = np.repeat(ins, k)
                rep_cell = np.repeat(cc, k)
                start_of = np.repeat(np.cumsum(k) - k, k)
                j = np.repeat(lo, k) + (np.arange(int(k.sum())) - start_of)
                off = (rep_ins - tp[j]) * tsg[j]
                cen = (off >= -500) & (off <= 499)
                fl = ((off >= -1000) & (off <= -901)) | ((off >= 901) & (off <= 1000))
                res["center"] += np.bincount(rep_cell[cen], minlength=n)
                res["flank"] += np.bincount(rep_cell[fl], minlength=n)
        if peaks is not None and c in peaks:
            ps, pe = peaks[c]
            # a fragment [s, e) overlaps a merged peak [ps, pe) if the last peak starting before e ends after s
            k = np.searchsorted(ps, e[m], "left") - 1
            hit = (k >= 0) & (pe[np.clip(k, 0, None)] > s[m])
            res["in_peaks"] += np.bincount(ci[m][hit], minlength=n)
    return res


def finish_atac_qc(acc, have_tss: bool, have_peaks: bool):
    np = _np()
    nf = acc["num_frags"]
    with np.errstate(divide="ignore", invalid="ignore"):
        dup = np.where(acc["reads"] > 0, 100 * (acc["reads"] - nf) / acc["reads"], np.nan)
        nuc = np.where(acc["nf"] > 0, acc["mono"] / acc["nf"], np.where(acc["mono"] > 0, np.inf, np.nan))
        if have_tss:
            fm = acc["flank"] / 200.0
            pos = fm[(fm > 0) & (nf > 0)]
            fm = np.where(fm > 0, fm, pos.mean() if len(pos) else np.nan)
            tsse = (acc["center"] / 1000.0) / fm
        else:
            tsse = np.full(len(nf), np.nan)
        frip = np.where(nf > 0, acc["in_peaks"] / nf, np.nan) if have_peaks else np.full(len(nf), np.nan)
    return {"num_frags": nf, "pct_duplicated_reads": dup, "nucleosomal_signal": nuc, "tss_enrichment": tsse, "frip": frip}


def read_manifest(path) -> List[dict]:
    pd = _pd()
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = {"lane", "fragments", "rna_matrix"}
    if not need <= set(df.columns):
        raise SystemExit(f"manifest needs columns {sorted(need)} (+ optional analysis_accession); has {list(df.columns)}")
    return df.to_dict(orient="records")


def build_pseudobulks(annotations, lanes: List[dict], out_root: Path, dataset: str, annotation_column=None,
                      barcode_column=None, subsample_column=None, lane_column=None, analysis_column=None,
                      atac_rna_map=None, tss=None, peaks=None, min_cells: int = 1, chunksize: int = 2_000_000) -> dict:
    """Primary pseudobulk layout {out}/{dataset}/pseudobulks/annotation-{ct}[-{subsample}]/ with fragments.tsv.gz,
    rna_counts_mtx.h5ad and per_cell_qc.tsv."""
    pd, np, ad = _pd(), _np(), _ad()
    ann = read_table_any(annotations)
    bcol = _pick_col(ann, barcode_column, ["barcode", "cell_barcode", "barcodes", "cell", "obs_names", "index", "Unnamed: 0"], "barcode")
    ccol = _pick_col(ann, annotation_column, ["cell_type", "celltype", "annotation", "cell_annotation", "cluster",
                                              "cell_type_annotation", "leiden", "seurat_clusters"], "cell-type annotation")
    scol = _pick_col(ann, subsample_column, ["subsample", "sample", "sample_accession", "in_vitro_system", "donor"], "subsample", False)
    lcol = _pick_col(ann, lane_column, ["lane", "library", "file_set", "multiplexed_sample", "gem_well"], "lane", False)
    acol = _pick_col(ann, analysis_column, ["analysis_accession", "analysis_set", "intermediate_analysis"], "analysis", False)
    lane_ids = [str(l["lane"]) for l in lanes]
    amap = read_barcode_map(atac_rna_map) if atac_rna_map else None
    ann = ann.copy()
    raw = ann[bcol].astype(str)
    has_suffix = raw.map(lambda b: "_" in GEM_SUFFIX_RE.sub("", b))
    if lcol:
        lane_of = ann[lcol].astype(str)
    elif has_suffix.all():
        lane_of = raw.map(lambda b: GEM_SUFFIX_RE.sub("", b).split("_", 1)[1])
    elif len(lane_ids) == 1:
        lane_of = pd.Series([lane_ids[0]] * len(ann))
    else:
        raise SystemExit("annotation barcodes carry no lane suffix and there are several lanes: pass --lane-column")
    ann["_full"] = [full_barcode(b, l) for b, l in zip(raw, lane_of)]
    ann["_lane"] = [f.split("_", 1)[1] if "_" in f else l for f, l in zip(ann["_full"], lane_of)]
    ann["_ct"] = ann[ccol].astype(str).map(safe_cell_type)
    ann["_sub"] = ann[scol].astype(str) if scol else ""
    lane_acc = {str(l["lane"]): (l.get("analysis_accession") or str(l["lane"])) for l in lanes}
    ann["_acc"] = ann[acol].astype(str) if acol else ann["_lane"].map(lambda x: lane_acc.get(x, x))
    ann = ann[(ann["_ct"] != "") & (~ann[ccol].astype(str).isin(["", "nan", "NA", "None"]))]
    dups = ann["_full"].duplicated(keep=False)
    if dups.any():
        raise SystemExit(f"{int(dups.sum())} annotation rows share a full barcode ({ann.loc[dups, '_full'].iloc[0]}): "
                         "the lane suffix does not disambiguate them")
    groups = ann.groupby(["_ct", "_sub"], sort=True)
    gkeys = [k for k, g in groups if len(g) >= min_cells]
    small = [(k, len(g)) for k, g in groups if len(g) < min_cells]
    gid = {k: i for i, k in enumerate(gkeys)}
    cells = ann[ann.set_index(["_ct", "_sub"]).index.isin(gkeys)].reset_index(drop=True)
    cells["_gid"] = [gid[(c, s)] for c, s in zip(cells["_ct"], cells["_sub"])]
    cell_idx = {b: i for i, b in enumerate(cells["_full"])}
    pb_root = Path(out_root) / dataset / "pseudobulks"
    gdir = {}
    for (ct, sub), i in gid.items():
        name = f"annotation-{ct}-{safe_cell_type(sub)}" if sub else f"annotation-{ct}"
        if sub and not IGVF_ACC_RE.match(sub):
            log.info("subsample %s is not an IGVF accession: upstream's annotation-%s-IGVF* glob would miss %s", sub, ct, name)
        gdir[i] = pb_root / name
        gdir[i].mkdir(parents=True, exist_ok=True)
    tss_idx = read_tss(tss) if tss else None
    peak_idx = read_peaks_merged(peaks) if peaks else None
    n = len(cells)
    acc = {k: np.zeros(n) for k in ("num_frags", "reads", "nf", "mono", "center", "flank", "in_peaks")}
    tmp_frag = {i: open(gdir[i] / "fragments.unsorted.tsv", "w") for i in gdir}
    rna_parts: Dict[int, list] = defaultdict(list)
    rna_seen = np.zeros(n, dtype=bool)
    lane_stats = []
    try:
        for L in lanes:
            lane = str(L["lane"])
            st = {"lane": lane, "fragments_rows": 0, "fragments_kept": 0, "rna_cells": 0, "rna_kept": 0}
            if L.get("fragments"):
                for chunk in pd.read_csv(L["fragments"], sep="\t", header=None, comment="#", dtype={0: str, 3: str},
                                         chunksize=chunksize, low_memory=False):
                    if chunk.shape[1] < 4:
                        raise SystemExit(f"{L['fragments']}: fragments need >= 4 columns")
                    chunk = chunk.iloc[:, :5].copy()
                    if chunk.shape[1] == 4:
                        chunk[4] = 1
                    chunk.columns = ["chr", "start", "end", "barcode", "count"]
                    st["fragments_rows"] += len(chunk)
                    uniq = pd.unique(chunk["barcode"])
                    tr = {u: full_barcode(u, lane, amap) for u in uniq}
                    fb = chunk["barcode"].map(tr)
                    ix = fb.map(cell_idx)
                    keep = ix.notna().to_numpy()
                    if not keep.any():
                        continue
                    sub = chunk[keep].copy()
                    sub["barcode"] = fb[keep].to_numpy()
                    sub["cell"] = ix[keep].astype(int).to_numpy()
                    st["fragments_kept"] += len(sub)
                    part = atac_cell_counts(sub, cell_idx, tss_idx, peak_idx)
                    for k in acc:
                        acc[k] += part[k]
                    sub["_g"] = cells["_gid"].to_numpy()[sub["cell"].to_numpy()]
                    for g, gg in sub.groupby("_g"):
                        gg[["chr", "start", "end", "barcode", "count"]].to_csv(tmp_frag[int(g)], sep="\t", header=False, index=False)
            if L.get("rna_matrix"):
                a = ad.read_h5ad(L["rna_matrix"])
                st["rna_cells"] = a.n_obs
                names = [full_barcode(b, lane) for b in map(str, a.obs_names)]
                a.obs_names = names
                if a.obs_names.duplicated().any():
                    a = a[~a.obs_names.duplicated()].copy()
                if "gene_symbol" not in a.var.columns:
                    a.var["gene_symbol"] = gene_symbols_of(a)
                m = a.obs_names.isin(list(cell_idx))
                a = a[m].copy()
                st["rna_kept"] = a.n_obs
                ii = np.array([cell_idx[b] for b in a.obs_names], dtype=int)
                rna_seen[ii] = True
                gids = cells["_gid"].to_numpy()[ii]
                for g in np.unique(gids):
                    rna_parts[int(g)].append(a[gids == g].copy())
            lane_stats.append(st)
    finally:
        for h in tmp_frag.values():
            h.close()
    q = finish_atac_qc(acc, tss_idx is not None, peak_idx is not None)
    rna_qc = {}
    summary_rows = []
    order = [c for c, _ in embedded_chrom_sizes()]
    for i, d in gdir.items():
        members = cells[cells["_gid"] == i]
        ii = members.index.to_numpy()
        # fragments: chrom-sizes order (unknown chromosomes after, lexicographic), then start, end
        up = d / "fragments.unsorted.tsv"
        fr = pd.read_csv(up, sep="\t", header=None, dtype={0: str, 3: str}) if up.stat().st_size else \
            pd.DataFrame(columns=[0, 1, 2, 3, 4])
        rank = {c: k for k, c in enumerate(order)}
        if len(fr):
            fr["_r"] = fr[0].map(lambda c: rank.get(c, len(rank)))
            fr = fr.sort_values(["_r", 0, 1, 2, 3], kind="mergesort").drop(columns="_r")
        plain = d / "fragments.tsv"
        fr.to_csv(plain, sep="\t", header=False, index=False)
        up.unlink()
        bgzip_and_index(plain, d / "fragments.tsv.gz")
        parts = rna_parts.get(i, [])
        if parts:
            ga = ad.concat(parts, axis=0, join="inner", merge="same") if len(parts) > 1 else parts[0]
            ga.write_h5ad(d / "rna_counts_mtx.h5ad")
            rq = rna_cell_qc(ga)
            for b, r in rq.iterrows():
                rna_qc[b] = r
        pcq = pd.DataFrame({"analysis_accession": members["_acc"].to_numpy(), "barcode": members["_full"].to_numpy(),
                            "subsample": members["_sub"].replace("", "NA").to_numpy()})
        for k in ("rna_read_count", "gene_count", "pct_mito", "pct_ribo"):
            pcq[k] = [rna_qc[b][k] if b in rna_qc else (0.0 if k in ("rna_read_count", "gene_count") else np.nan)
                      for b in members["_full"]]
        for k in ("num_frags", "pct_duplicated_reads", "nucleosomal_signal", "tss_enrichment", "frip"):
            pcq[k] = q[k][ii]
        pcq = pcq[PER_CELL_QC_COLS]
        pcq.to_csv(d / "per_cell_qc.tsv", sep="\t", index=False, na_rep="NA")
        summary_rows.append({"directory": d.name, "cell_type": members["_ct"].iloc[0], "subsample": members["_sub"].iloc[0],
                             "cells": len(members), "cells_with_rna": int(rna_seen[ii].sum()),
                             "cells_with_fragments": int((q["num_frags"][ii] > 0).sum()), "fragments": int(q["num_frags"][ii].sum()),
                             "umis": float(pcq["rna_read_count"].sum())})
    summ = pd.DataFrame(summary_rows)
    cts = sorted(summ["cell_type"].unique()) if len(summ) else []
    return {"pseudobulks": str(pb_root), "dataset": dataset, "cell_types": cts, "groups": summ, "lanes": lane_stats,
            "annotated_cells": int(len(ann)), "cells_used": n, "groups_below_min_cells": small,
            "tss_enrichment": tss_idx is not None, "frip": peak_idx is not None,
            "atac_rna_map": bool(amap), "columns": {"barcode": bcol, "annotation": ccol, "subsample": scol, "lane": lcol,
                                                    "analysis": acol}}


# ---------------------------------------------------------------------------
# Portal retrieval (processed-first)
# ---------------------------------------------------------------------------

def _lineage():
    try:
        import portal_lineage as pl  # type: ignore
    except ImportError:
        from igvfagent import portal_lineage as pl  # type: ignore
    return pl


def select_portal_inputs(plan: dict, g, include_lab: bool = False) -> dict:
    """Principal analysis set(s) + cell annotations; per intermediate analysis set (10x lane) the best
    fragments file and RNA matrix (h5ad preferred), uniform-pipeline sets first."""
    nodes = g.nodes
    prods = plan.get("products") or {}
    principal = [s for s in plan.get("analysis_sets") or [] if s.get("file_set_type") == "principal analysis"]
    ann = [r for r in prods.get("cell_annotations", []) if r.get("file_set_type") in ("principal analysis", "")] or \
        list(prods.get("cell_annotations", []))
    by_set: Dict[str, dict] = {}
    for p in ("fragments", "rna_matrix"):
        for r in prods.get(p, []):
            fs = r.get("file_set") or ""
            if r.get("file_set_type") == "principal analysis":
                continue
            ent = by_set.setdefault(fs, {"analysis_accession": fs, "uniform": r.get("uniform_pipeline_status") == "completed",
                                         "file_set_type": r.get("file_set_type")})
            if p not in ent:
                if p == "rna_matrix" and str(r.get("file_format")).lower() not in ("h5ad",) and \
                        any(x.get("file_set") == fs and str(x.get("file_format")).lower() == "h5ad" for x in prods.get(p, [])):
                    continue
                ent[p] = r
    lanes = list(by_set.values())
    if not include_lab and any(l["uniform"] for l in lanes):
        lanes = [l for l in lanes if l["uniform"]]
    # lane accession: the (multiplexed) sample of the lane's input measurement sets, else the analysis set
    idx = {v.get("accession"): k for k, v in nodes.items() if v.get("accession")}
    for l in lanes:
        sid = idx.get(l["analysis_accession"])
        inputs = [a for a, b, lab in g.edges if b == sid and lab == "input_for"] + \
                 [b for a, b, lab in g.edges if a == sid and lab == "input_file_sets"]
        samples = sorted({(nodes.get(b) or {}).get("accession") for a, b, lab in g.edges
                          if a in inputs and lab == "samples" and (nodes.get(b) or {}).get("kind") == "sample"} - {None})
        l["lane"] = samples[0] if len(samples) == 1 else l["analysis_accession"]
        l["lane_source"] = "sample" if len(samples) == 1 else "analysis set"
    return {"principal": principal, "annotations": ann, "lanes": lanes}


def fetch_inputs(accession: str, dest: Path, fetch=None, downloader=None, max_gb: float = 20.0, dry_run: bool = False,
                 include_lab: bool = False, have_credentials: bool = False, depth: int = 6) -> dict:
    pd = _pd()
    pl = _lineage()
    if fetch is None:
        fetch = pl.live_fetch()
    g = pl.walk(accession, fetch=fetch, max_depth=depth, fetch_qc=False)
    plan = pl.build_plan(g, have_credentials=have_credentials)
    sel = select_portal_inputs(plan, g, include_lab)
    want = []
    if sel["annotations"]:
        want.append(("cell_annotations", "", sel["annotations"][0]))
    for l in sel["lanes"]:
        for p in ("fragments", "rna_matrix"):
            if l.get(p):
                want.append((p, l["lane"], l[p]))
    budget = float(max_gb)
    got, skipped, rows = {}, [], []
    for prod, lane, r in want:
        name = (r.get("href") or "").rsplit("/", 1)[-1] or f"{r['accession']}.{r.get('file_format') or 'dat'}"
        target = Path(dest) / name
        why = None
        if "controlled" in str(r.get("access")) and not have_credentials:
            why = "controlled access and no IGVF credentials configured"
        elif (r.get("size_gb") or 0) > budget:
            why = f"{r.get('size_gb', 0):.2f} GB exceeds the remaining budget {budget:.2f} GB"
        elif not r.get("href"):
            why = "no download href on the file record"
        if why:
            skipped.append({"product": prod, "lane": lane, "accession": r["accession"], "why": why})
            continue
        budget -= r.get("size_gb") or 0
        if dry_run:
            print(f"  would fetch {prod:17} {r['accession']} {r.get('size_gb', 0):.2f} GB -> {target}")
        else:
            if downloader is None:
                try:
                    from igvfagent.raw_data_pipeline import portal_download as downloader  # type: ignore
                except Exception:
                    from raw_data_pipeline import portal_download as downloader  # type: ignore
            target.parent.mkdir(parents=True, exist_ok=True)
            if not (target.is_file() and target.stat().st_size > 0):
                downloader(r["href"], target)
            print(f"Wrote: {target}")
        got[(prod, lane)] = str(target)
    for l in sel["lanes"]:
        rows.append({"lane": l["lane"], "analysis_accession": l["analysis_accession"],
                     "fragments": got.get(("fragments", l["lane"]), ""), "rna_matrix": got.get(("rna_matrix", l["lane"]), ""),
                     "fragments_accession": (l.get("fragments") or {}).get("accession", ""),
                     "rna_matrix_accession": (l.get("rna_matrix") or {}).get("accession", ""),
                     "uniform_pipeline": l["uniform"], "lane_source": l["lane_source"]})
    manifest = pd.DataFrame(rows, columns=["lane", "analysis_accession", "fragments", "rna_matrix", "fragments_accession",
                                           "rna_matrix_accession", "uniform_pipeline", "lane_source"])
    return {"plan": plan, "selection": sel, "manifest": manifest, "skipped": skipped,
            "annotations_path": got.get(("cell_annotations", "")), "budget_left_gb": budget, "graph": g}


# ---------------------------------------------------------------------------
# Spec validation
# ---------------------------------------------------------------------------

def validate_outputs(qc_guide=None, rna_matrix=None, fragments=None, chrom_sizes=None) -> List[dict]:
    pd = _pd()
    checks: List[dict] = []

    def add(ok, spec, msg):
        checks.append({"ok": bool(ok), "spec": spec, "check": msg})
    guide_set = None
    if qc_guide:
        with open(qc_guide, "rb") as fh:
            add(fh.read(2) == b"\x1f\x8b", "QC guide", "gzip-compressed")
        g = pd.read_csv(qc_guide, sep="\t", dtype=str, keep_default_na=False)
        add(list(g.columns) == GUIDE_COLS, "QC guide", f"header is barcode, subsample, analysis_accession ({list(g.columns)})")
        if "barcode" in g.columns:
            add(not g["barcode"].duplicated().any(), "QC guide", "every barcode unique")
            bad = [b for b in g["barcode"] if not BARCODE_SPEC_RE.match(b)]
            add(not bad, "QC guide", f"barcodes are {{16-bp}}_{{IGVF 10x lane accession}} ({len(bad)} do not match{': ' + bad[0] if bad else ''})")
            guide_set = set(g["barcode"])
        for c in ("subsample", "analysis_accession"):
            if c in g.columns:
                add((g[c] != "").all(), "QC guide", f"{c} filled for every row")
    if rna_matrix:
        import scipy.io  # type: ignore
        d = Path(rna_matrix)

        def pick(n):
            for s in (".gz", ""):
                if (d / (n + s)).exists():
                    return d / (n + s)
            return None
        mp, bp, fp = pick("matrix.mtx"), pick("barcodes.tsv"), pick("features.tsv")
        add(mp and bp and fp, "RNA matrix", "matrix.mtx, barcodes.tsv, features.tsv present")
        if mp and bp and fp:
            with _opener(mp) as fh:
                head = fh.readline().strip()
            add(head == "%%MatrixMarket matrix coordinate integer general", "RNA matrix", f"header {head!r}")
            M = scipy.io.mmread(str(mp))
            bcs = [l.rstrip("\n") for l in _opener(bp)]
            fts = [l.rstrip("\n") for l in _opener(fp)]
            add(all("\t" not in x for x in bcs + fts), "RNA matrix", "barcodes.tsv / features.tsv single-column, no header")
            add(M.shape == (len(fts), len(bcs)), "RNA matrix", f"genes x cells {M.shape} = features x barcodes ({len(fts)}, {len(bcs)})")
            add(len(set(bcs)) == len(bcs), "RNA matrix", "no duplicate cell barcodes")
            add(len(set(fts)) == len(fts), "RNA matrix", "no duplicate gene symbols")
            ens = [f for f in fts if re.match(r"^ENS[A-Z]*G\d+", f)]
            add(not ens, "RNA matrix", f"features are gene symbols, not Ensembl IDs ({len(ens)} Ensembl IDs)")
            if guide_set is not None:
                add(set(bcs) == guide_set, "RNA matrix", f"barcodes.tsv is the QC guide's barcode set ({len(bcs)} vs {len(guide_set)})")
    if fragments:
        order = [c for c, _ in read_chrom_sizes(chrom_sizes)]
        rank = {c: i for i, c in enumerate(order)}
        fr = pd.read_csv(fragments, sep="\t", header=None, comment="#", dtype={0: str, 3: str})
        add(fr.shape[1] == 5, "fragments", f"5 columns ({fr.shape[1]})")
        add(fr[0].isin(list(rank)).all(), "fragments", "every chromosome is in the chrom-sizes file")
        r = fr[0].map(rank).fillna(-1).to_numpy()
        import numpy as np  # type: ignore
        key_ok = bool((np.diff(r) >= 0).all())
        same = np.diff(r) == 0
        pos_ok = bool((np.diff(fr[1].to_numpy())[same] >= 0).all())
        add(key_ok and pos_ok, "fragments", "sorted in chrom-sizes order, then start")
        add(Path(str(fragments) + ".tbi").exists(), "fragments", "tabix index present")
        if guide_set is not None:
            seen = set(fr[3])
            add(guide_set <= seen and seen <= guide_set, "fragments", "barcodes are exactly the QC guide's (every guide barcode has fragments)")
    return checks


# ---------------------------------------------------------------------------
# sce2g-pipeline hand-off
# ---------------------------------------------------------------------------

def sce2g_prep(config_tsv, out: Path, chrom_sizes: "Optional[str]" = None) -> dict:
    pd = _pd()
    s = _sce2g()
    cfg = pd.read_csv(config_tsv, sep="\t", dtype=str, keep_default_na=False)
    out.mkdir(parents=True, exist_ok=True)
    if not chrom_sizes:
        chrom_sizes = str(write_chrom_sizes(out / "IGVF.DACC.GRCh38.chrom.sizes.tsv"))
    res = []
    for r in cfg.to_dict(orient="records"):
        cd = out / r["cluster"]
        ta = s.frag_to_tagalign(r["atac_frag_file"], chrom_sizes, cd)
        M, genes, cells = s.read_rna_matrix(r["rna_matrix_file"])
        _, exp, umis, ncells = s.rna_features(M, genes)
        exp = exp.reset_index().rename(columns={"index": "gene"})
        p = cd / "RNA_pseudobulk.tsv"
        exp.to_csv(p, sep="\t", index=False)
        print(f"TSV: {p}")
        res.append({"cluster": r["cluster"], "fragment_count": ta["fragment_count"], "atac_barcodes": ta["n_barcodes"],
                    "tagAlign": ta["tagAlign"], "rna_cells": ncells, "umis": umis, "genes": len(genes),
                    "rna_pseudobulk": str(p), "tpm_sum": float(exp["RnaPseudobulkTPM"].sum())})
    cc = cfg.copy()
    cc["abc_dir"] = ""
    ccp = out / "config_cell_clusters.tsv"
    cc.to_csv(ccp, sep="\t", index=False)
    print(f"TSV: {ccp}")
    return {"clusters": res, "cell_cluster_config": str(ccp), "chrom_sizes": chrom_sizes}


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def _threshold_args(p, defaults: dict) -> None:
    for flag, key in FLAG_MAP.items():
        p.add_argument(flag, dest=key, type=float, default=None, help=f"default {rnum(defaults[key])}")


def cmd_fetch(args) -> int:
    out = out_dir_for(args, args.accession)
    dest = Path(args.dest) if args.dest else DATA_ROOT / args.accession
    have = bool(os.environ.get("IGVF_API_KEY") or os.environ.get("IGVF_ACCESS_KEY"))
    r = fetch_inputs(args.accession, dest, max_gb=args.max_gb, dry_run=args.dry_run, include_lab=args.include_lab,
                     have_credentials=have, depth=args.depth)
    man = out / "portal_manifest.tsv"
    r["manifest"].to_csv(man, sep="\t", index=False)
    print(f"TSV: {man}")
    sel = r["selection"]
    lines = [f"# Principal pseudobulk inputs for {args.accession}", "",
             f"Processed-first plan via portal_lineage ({r['plan']['n_nodes']} nodes). Verdict: {r['plan']['verdict']}.", "",
             "## Principal analysis sets", ""]
    lines += [f"- {s['accession']} ({s.get('uniform_pipeline_status') or 'no uniform status'}): {s.get('summary') or ''}"
              for s in sel["principal"]] or ["- none visible (unreleased principal sets appear as HTTP 403 in the lineage)"]
    lines += ["", "## Cell annotations", ""]
    lines += [f"- {a['accession']} {a.get('content_type')} ({a.get('size_gb')} GB) in {a.get('file_set')}" for a in sel["annotations"]] or \
        ["- none found: pass your own annotation table to build-pseudobulks"]
    lines += ["", "## Lanes (intermediate analysis sets)", "",
              md_table(["lane", "analysis set", "fragments", "RNA matrix", "uniform"],
                       [[l["lane"], l["analysis_accession"], (l.get("fragments") or {}).get("accession", "-"),
                         (l.get("rna_matrix") or {}).get("accession", "-"), l["uniform"]] for l in sel["lanes"]])]
    if r["skipped"]:
        lines += ["", "## Skipped", ""] + [f"- {s['product']} {s['accession']}: {s['why']}" for s in r["skipped"]]
    lines += ["", "## Next", "", "```",
              f"igvfagent principal-pseudobulks build-pseudobulks --annotations {r['annotations_path'] or '<cells.tsv>'} "
              f"--manifest {man} --dataset {args.accession} --tss <TSS.bed> --peaks <peaks.bed>", "```"]
    write_summary(out, {"accession": args.accession, "dry_run": args.dry_run, "lanes": sel["lanes"], "principal": sel["principal"],
                        "annotations": sel["annotations"], "skipped": r["skipped"], "manifest": str(man),
                        "budget_left_gb": r["budget_left_gb"]}, lines)
    _record("fetch", args.accession, [args.accession], [str(man)])
    return 0


def cmd_build(args) -> int:
    out = out_dir_for(args, args.label or args.dataset)
    if args.manifest:
        lanes = read_manifest(args.manifest)
    else:
        if not args.lane or len(args.lane) != len(args.fragments or []) or len(args.lane) != len(args.rna_matrix or []):
            raise SystemExit("pass --manifest, or equal numbers of --lane, --fragments and --rna-matrix")
        lanes = [{"lane": l, "fragments": f, "rna_matrix": m} for l, f, m in zip(args.lane, args.fragments, args.rna_matrix)]
    root = Path(args.pseudobulk_root) if args.pseudobulk_root else out
    r = build_pseudobulks(args.annotations, lanes, root, args.dataset, args.annotation_column, args.barcode_column,
                          args.subsample_column, args.lane_column, args.analysis_column, args.atac_rna_map, args.tss,
                          args.peaks, args.min_cells)
    gp = out / "pseudobulk_summary.tsv"
    r["groups"].to_csv(gp, sep="\t", index=False)
    print(f"TSV: {gp}")
    cfgp = out / "config_QC_pseudobulks.json"
    cfgp.write_text(json.dumps({"QC_plots_dir": str(out / "plots"), "pseudobulk_dir": str(root), "out_dir": str(out / "multiome_data"),
                                "transcriptome": args.gtf or f"<{IGVF_TRANSCRIPTOME}.gtf.gz>", "chrom_sizes": "",
                                "datasets": {args.dataset: r["cell_types"]}}, indent=2))
    print(f"JSON: {cfgp}")
    lines = [f"# Primary pseudobulks: {args.dataset}", "", f"{r['cells_used']:,} annotated cells in {len(r['groups'])} "
             f"cell type x subsample directories under `{r['pseudobulks']}`.", "",
             md_table(list(r["groups"].columns), r["groups"].values.tolist()), "",
             "## Lanes", "", md_table(["lane", "fragments rows", "kept", "RNA cells", "kept"],
                                      [[l["lane"], l["fragments_rows"], l["fragments_kept"], l["rna_cells"], l["rna_kept"]] for l in r["lanes"]]),
             "", f"TSS enrichment computed: {r['tss_enrichment']}; FRiP computed: {r['frip']}; ATAC->RNA barcode map: {r['atac_rna_map']}.",
             "", "Next: `igvfagent principal-pseudobulks run --config " + str(cfgp) + " --auto-qc`"]
    summ = {k: v for k, v in r.items() if k != "groups"}
    summ.update({"groups": r["groups"].to_dict(orient="records"), "config": str(cfgp)})
    write_summary(out, summ, lines)
    _record("build-pseudobulks", args.dataset, [args.annotations], [r["pseudobulks"]])
    return 0


def cmd_datatable(args) -> int:
    out = out_dir_for(args, args.label or args.cell_type)
    op = Path(args.out) if args.out else out / f"{safe_cell_type(args.cell_type)}_per_cell_qc.tsv"
    r = build_qc_datatable(args.pseudobulks, args.cell_type, op)
    print(f"[info] Wrote {r['cells']} cells from {r['directories']} subsample director{'y' if r['directories'] == 1 else 'ies'} to {op}")
    print(f"TSV: {op}")
    write_summary(out, r, [f"# Per-cell QC datatable: {args.cell_type}", "", f"{r['cells']:,} cells from {r['directories']} "
                           f"directories -> `{op}`."] + [f"- skipped {d}: {w}" for d, w in r["skipped"]])
    return 0


def cmd_explore(args) -> int:
    pd = _pd()
    out = out_dir_for(args, args.label or "explore")
    meta = read_meta(args.meta)
    skip = [] if args.upstream_compat else na_columns(meta)
    if args.sets:
        pieces = [s for part in args.sets for s in part.split(";")]
        sets = [parse_threshold_string(s, EXPLORE_DEFAULTS) for s in pieces if s.strip()] or [dict(EXPLORE_DEFAULTS)]
    else:
        sets = [thresholds_from_args(args, EXPLORE_DEFAULTS)]
    text, rows, set_rows, sub_rows = explore_sets(meta, sets, args.show_subsamples, skip)
    print(text)
    (out / "explore_report.txt").write_text(text + "\n")
    print(f"Wrote: {out / 'explore_report.txt'}")
    pd.DataFrame(rows).to_csv(out / "threshold_breakdown.tsv", sep="\t", index=False)
    pd.DataFrame(set_rows).to_csv(out / "threshold_sets.tsv", sep="\t", index=False)
    print(f"TSV: {out / 'threshold_breakdown.tsv'}")
    print(f"TSV: {out / 'threshold_sets.tsv'}")
    if sub_rows:
        pd.DataFrame(sub_rows).to_csv(out / "per_subsample.tsv", sep="\t", index=False)
        print(f"TSV: {out / 'per_subsample.tsv'}")
    write_summary(out, {"sets": set_rows, "skipped_all_na": skip}, ["# QC threshold exploration", "", "```", text, "```"])
    return 0


def cmd_qc_filter(args) -> int:
    out = out_dir_for(args, args.label or "qc_filter")
    p = thresholds_from_args(args, FILTER_DEFAULTS)
    r = run_qc_filter(args.meta, out, p, plots=not args.no_plots, upstream_compat=args.upstream_compat)
    lines = ["# QC filter (plot_per_cell_qc.R)", "", f"{r['kept']:,} of {r['cells']:,} cells pass (RNA QC {r['kept_rna']:,}, "
             f"ATAC QC {r['kept_atac']:,}).", "", "## Thresholds", "",
             md_table(["threshold", "value"], [[THRESHOLD_DOC[k], rnum(p[k])] for k in THRESH_ORDER]), "",
             "## Cells per subsample after all QC", "", md_table(list(r["metrics"].columns), r["metrics"].values.tolist())]
    if r["skipped_thresholds"]:
        lines += ["", f"All-NA metrics skipped: {', '.join(r['skipped_thresholds'])}"]
    write_summary(out, {k: v for k, v in r.items() if k not in ("metrics", "guide")}, lines)
    _record("qc-filter", args.label or "", [args.meta], [r["files"]["guide"]])
    return 0


def cmd_filter_atac(args) -> int:
    out = out_dir_for(args, args.label or f"atac_{args.cell_type}")
    op = Path(args.out) if args.out else out / f"atac_fragments_{safe_cell_type(args.cell_type)}.tsv.gz"
    r = filter_atac(args.qc_guide, args.pseudobulks, args.cell_type, args.chrom_sizes, op, clean=args.clean)
    print(f"Wrote: {r['out']}")
    if r["index"]:
        print(f"Wrote: {r['index']}")
    lines = [f"# Filtered ATAC fragments: {args.cell_type}", "", f"{r['rows_written']:,} fragments written "
             f"({r['rows_retained']:,} of {r['rows_total']:,} passed the QC guide; {sum(r['chromosomes_dropped'].values()):,} "
             f"on chromosomes absent from the sizes file). All {r['guide_barcodes']:,} guide barcodes found.", "",
             md_table(["directory", "retained", "total"], [[d["dir"], d["retained"], d["total"]] for d in r["directories"]])]
    write_summary(out, r, lines)
    return 0


def cmd_filter_rna(args) -> int:
    out = out_dir_for(args, args.label or f"rna_{args.cell_type}")
    op = args.out or str(out / f"rna_count_matrix_{safe_cell_type(args.cell_type)}.mtx")
    r = filter_rna(args.qc_guide, args.pseudobulks, args.cell_type, op, args.gtf, args.ensembl_ids_as_genes,
                   args.standard_chromosomes_only, args.log, integer=not args.upstream_compat)
    for v in r["outputs"].values():
        print(f"Wrote: {v}")
    if r.get("log"):
        print(f"Wrote: {r['log']}")
    lines = [f"# Filtered RNA count matrix: {args.cell_type}", "", f"{r['cells']:,} cells x {r['genes']:,} "
             f"{'Ensembl IDs' if args.ensembl_ids_as_genes else 'gene symbols'} ({r['format']}); sanity checks passed.", "",
             md_table(["directory", "cells", "passed QC"], [[d["dir"], d["cells"], d["passed"]] for d in r["directories"]])]
    write_summary(out, r, lines)
    return 0


def cmd_package(args) -> int:
    out = out_dir_for(args, args.label or "package")
    t = package_rna(Path(args.rna_matrix), Path(args.out) if args.out else out / (Path(args.rna_matrix).name + ".tar.gz"))
    print(f"Wrote: {t}")
    write_summary(out, {"tar": str(t)}, ["# RNA count matrix archive", "", f"`{t}` (matrix.mtx, barcodes.tsv, features.tsv)"])
    return 0


def cmd_config_table(args) -> int:
    out = out_dir_for(args, args.label or args.dataset)
    df = config_table(Path(args.out_root), args.dataset, args.cell_types, args.upstream_compat)
    p = Path(args.out) if args.out else out / f"{args.dataset}_config.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)
    print(f"TSV: {p}")
    write_summary(out, {"config_table": str(p), "clusters": list(args.cell_types)},
                  ["# Config table", "", md_table(CONFIG_TABLE_COLS[:3], df[CONFIG_TABLE_COLS[:3]].values.tolist())])
    return 0


def run_workflow(cfg: dict, out_default: Path, dry_run: bool = False, auto_qc: bool = False, thresholds=None,
                 plots: bool = False, upstream_compat: bool = False, attempt: int = 1) -> dict:
    qc_dir = Path(cfg.get("QC_plots_dir") or out_default / "plots")
    data_dir = Path(cfg["pseudobulk_dir"])
    out_dir = Path(cfg.get("out_dir") or out_default / "multiome_data")
    if not out_dir.is_absolute() and cfg.get("out_dir"):
        out_dir = out_default / out_dir
    gtf = cfg.get("transcriptome") or None
    if gtf and (gtf.startswith("<") or not Path(gtf).exists()):
        gtf = None
    chrom = cfg.get("chrom_sizes") or None
    if chrom and not Path(chrom).exists():
        chrom = None  # the Snakefile's reference/ path is embedded here
    datasets = cfg["datasets"]
    jobs, results = [], []
    for ds, cts in datasets.items():
        for ct in cts:
            guide = qc_dir / ds / ct / "filtered_barcodes_with_subsamples.tsv.gz"
            pbs = data_dir / ds / "pseudobulks"
            dirs = find_pseudobulk_dirs(pbs, ct) if pbs.exists() else []
            mb = _dir_size_mb(dirs)
            jobs += [{"rule": "qc_filter" if (auto_qc and not guide.exists()) else None, "dataset": ds, "cell_type": ct},
                     {"rule": "atac_fragment_file", "dataset": ds, "cell_type": ct, "mem_mb": determine_mem_mb(mb, True, attempt)},
                     {"rule": "rna_count_matrix", "dataset": ds, "cell_type": ct, "mem_mb": determine_mem_mb(mb, True, attempt)}]
    jobs = [j for j in jobs if j["rule"]] + [{"rule": "make_config_table", "dataset": ds} for ds in datasets]
    if dry_run:
        return {"jobs": jobs, "dry_run": True, "out_dir": str(out_dir)}
    tables = []
    for ds, cts in datasets.items():
        for ct in cts:
            guide = qc_dir / ds / ct / "filtered_barcodes_with_subsamples.tsv.gz"
            pbs = data_dir / ds / "pseudobulks"
            rec: Dict[str, Any] = {"dataset": ds, "cell_type": ct}
            if not guide.exists():
                if not auto_qc:
                    raise SystemExit(f"missing QC guide {guide}: run qc-filter (Step 1) or pass --auto-qc")
                dt = qc_dir / ds / ct / f"{ct}_per_cell_qc.tsv"
                build_qc_datatable(pbs, ct, dt)
                q = run_qc_filter(dt, qc_dir / ds / ct, thresholds or dict(FILTER_DEFAULTS), plots=plots, upstream_compat=upstream_compat)
                rec["qc_kept"] = q["kept"]
                rec["qc_cells"] = q["cells"]
            cdir = out_dir / ds / ct
            a = filter_atac(guide, pbs, ct, chrom, cdir / f"atac_fragments_{ds}_{ct}.tsv.gz")
            rr = filter_rna(str(guide), pbs, ct, str(cdir / f"rna_count_matrix_{ds}_{ct}.mtx"), gtf=gtf,
                            standard_only=bool(gtf), integer=not upstream_compat)
            rec.update({"fragments": a["rows_written"], "atac_file": a["out"], "rna_cells": rr["cells"], "rna_genes": rr["genes"],
                        "rna_dir": str(cdir / f"rna_count_matrix_{ds}_{ct}"), "gene_mapping": "gtf" if gtf else "var gene_symbol"})
            print(f"Wrote: {a['out']}")
            print(f"Wrote: {rec['rna_dir']}")
            results.append(rec)
        tp = out_dir / "config" / "tables" / f"{ds}_config.tsv"
        tp.parent.mkdir(parents=True, exist_ok=True)
        config_table(out_dir, ds, cts, upstream_compat).to_csv(tp, sep="\t", index=False)
        print(f"TSV: {tp}")
        tables.append(str(tp))
    return {"jobs": jobs, "clusters": results, "config_tables": tables, "out_dir": str(out_dir)}


def cmd_run(args) -> int:
    out = out_dir_for(args, args.label or "run")
    if args.config:
        cfg = load_run_config(args.config)
    else:
        if not args.pseudobulk_dir or not args.datasets:
            raise SystemExit("pass --config, or --pseudobulk-dir and --datasets ds:ct1,ct2")
        cfg = {"QC_plots_dir": args.qc_plots_dir, "pseudobulk_dir": args.pseudobulk_dir, "out_dir": args.out_root,
               "transcriptome": args.transcriptome, "chrom_sizes": args.chrom_sizes,
               "datasets": {d.split(":", 1)[0]: [c for c in d.split(":", 1)[1].split(",") if c] for d in args.datasets}}
    for k, v in (("QC_plots_dir", args.qc_plots_dir), ("out_dir", args.out_root), ("transcriptome", args.transcriptome),
                 ("chrom_sizes", args.chrom_sizes)):
        if v and args.config:
            cfg[k] = v
    r = run_workflow(cfg, out, args.dry_run, args.auto_qc, thresholds_from_args(args, FILTER_DEFAULTS), not args.no_plots,
                     args.upstream_compat)
    lines = ["# Principal pseudobulks (Snakefile)", ""]
    lines += ["## Jobs" + (" (dry run)" if args.dry_run else ""), "",
              md_table(["rule", "dataset", "cell type", "mem_mb"], [[j["rule"], j["dataset"], j.get("cell_type", ""), j.get("mem_mb", "")] for j in r["jobs"]])]
    if not args.dry_run:
        lines += ["", "## Clusters", "", md_table(["dataset", "cell type", "fragments", "RNA cells", "genes"],
                                                  [[c["dataset"], c["cell_type"], c["fragments"], c["rna_cells"], c["rna_genes"]] for c in r["clusters"]]),
                  "", "Config tables (feed `sce2g-pipeline run --cluster-config`, or `principal-pseudobulks sce2g-prep`):", ""]
        lines += [f"- `{t}`" for t in r["config_tables"]]
    if args.dry_run:
        for j in r["jobs"]:
            print(f"  job {j['rule']:20} {j['dataset']} {j.get('cell_type', '')} {('mem_mb=' + str(j['mem_mb'])) if 'mem_mb' in j else ''}")
    write_summary(out, r, lines)
    return 0


def cmd_validate(args) -> int:
    out = out_dir_for(args, args.label or "validate")
    checks = validate_outputs(args.qc_guide, args.rna_matrix, args.fragments, args.chrom_sizes)
    for c in checks:
        print(("  ok    " if c["ok"] else "  FAIL  ") + f"[{c['spec']}] {c['check']}")
    bad = [c for c in checks if not c["ok"]]
    _pd().DataFrame(checks).to_csv(out / "spec_checks.tsv", sep="\t", index=False)
    print(f"TSV: {out / 'spec_checks.tsv'}")
    write_summary(out, {"checks": checks, "failed": len(bad)}, ["# File-spec validation", "",
                                                                md_table(["ok", "spec", "check"], [[c["ok"], c["spec"], c["check"]] for c in checks])])
    return 1 if bad else 0


def cmd_sce2g_prep(args) -> int:
    out = out_dir_for(args, args.label or "sce2g_prep")
    r = sce2g_prep(args.config_table, out, args.chrom_sizes)
    lines = ["# sce2g-pipeline inputs", "", md_table(["cluster", "fragments", "ATAC barcodes", "RNA cells", "UMIs"],
                                                     [[c["cluster"], c["fragment_count"], c["atac_barcodes"], c["rna_cells"], c["umis"]] for c in r["clusters"]]),
             "", "Next (after ABC on the tagAlign files):", "", "```",
             f"igvfagent sce2g-pipeline run --cluster-config {r['cell_cluster_config']} --label principal_pseudobulks", "```"]
    write_summary(out, r, lines)
    return 0


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _synthetic(d: Path) -> dict:
    """Two 10x lanes, two subsamples, two cell types; planted low-quality cells, a barcode collision across
    lanes, ATAC barcodes needing the ATAC->RNA map (lane 1), a GEM '-1' suffix, a chromosome absent from the
    chrom sizes, duplicate-symbol / nonstandard-contig genes."""
    pd, np, sp, ad = _pd(), _np(), _sp(), _ad()
    rng = np.random.default_rng(7)
    lanes = ["IGVFSM0001LANA", "IGVFSM0002LANB"]
    acc = {"IGVFSM0001LANA": "IGVFDS0001INTA", "IGVFSM0002LANB": "IGVFDS0002INTB"}
    subs = ["IGVFSM1111SUBA", "IGVFSM2222SUBB"]
    bases = np.array(list("ACGT"))

    def bc():
        return "".join(rng.choice(bases, 16))
    genes = [("ENSG00000000001.1", "chr1", "GENEA"), ("ENSG00000000002.3", "chr1", "GENEB"),
             ("ENSG00000000003.2", "chr2", "DUPSYM"), ("ENSG00000000004.1", "chr3", "DUPSYM"),
             ("ENSG00000000005.1", "chrUn_KI270302v1", "SCAFONLY"), ("ENSG00000000006.1", "chr4", "SHARED"),
             ("ENSG00000000007.1", "chr1_KI270706v1_random", "SHARED"), ("ENSG00000198804.2", "chrM", "MT-CO1"),
             ("ENSG00000142541.17", "chr19", "RPL13A"), ("ENSG00000000010.1", "chr5", "GENEC")]
    gtf = d / "genes.gtf"
    with open(gtf, "w") as fh:
        fh.write("#!genome-build GRCh38\n")
        for gid, ch, name in genes:
            fh.write(f'{ch}\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "{gid}"; gene_name "{name}";\n')
            fh.write(f'{ch}\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "{gid}"; transcript_id "T{gid}"; gene_name "{name}";\n')
    tss = d / "tss.bed"
    tss_pos = [(f"chr{c}", 100_000 * k + 50_000) for c in (1, 2) for k in range(1, 6)]
    tss.write_text("".join(f"{c}\t{p}\t{p + 1}\tg{i}\t0\t+\n" for i, (c, p) in enumerate(tss_pos)))
    peaks = d / "peaks.bed"
    peaks.write_text("".join(f"{c}\t{p - 300}\t{p + 300}\n" for c, p in tss_pos))
    cells, frag_rows = [], {l: [] for l in lanes}
    rna = {l: [] for l in lanes}
    atac_map = {}
    shared_seq = bc()
    k = 0
    for li, lane in enumerate(lanes):
        for ct in ("k562", "hepg2"):
            for si, s in enumerate(subs):
                for j in range(6):
                    seq = shared_seq if (ct == "k562" and si == 0 and j == 0) else bc()
                    bad = {0: None, 1: None, 2: None, 3: None, 4: "lowrna", 5: "lowatac"}[j]
                    cells.append({"barcode": f"{seq}-1", "lane": lane, "cell_type": ct, "subsample": s, "bad": bad or "",
                                  "full": f"{seq}_{lane}", "k": k})
                    aseq = bc() if li == 0 else seq
                    if li == 0:
                        atac_map[aseq] = seq
                    nfr = 50 if bad == "lowatac" else 1500 + 10 * j
                    for f in range(nfr):
                        c, p = tss_pos[(f + k) % len(tss_pos)]
                        near = (f % 3 != 0) and bad != "lowatac"
                        if near:
                            start = p + int(rng.integers(-200, 150))
                        elif f % 60 == 0:   # TSS-flank background (+/- 901-1000 bp)
                            start = p + int(rng.integers(905, 990)) * (1 if f % 120 else -1)
                        else:
                            start = p + int(rng.integers(3000, 40000))
                        ln = 100 if f % 4 else 200
                        frag_rows[lane].append((c, start, start + ln, f"{aseq}-1", 2 if f % 5 == 0 else 1))
                    if j == 0:
                        frag_rows[lane].append(("chrFOO", 10, 110, f"{aseq}-1", 1))
                    umi = 60 if bad == "lowrna" else 3000 + 50 * j
                    rna[lane].append((f"{seq}-1", umi, j, ct))
                    k += 1
    # an unannotated barcode in the fragments + matrix (must be ignored)
    frag_rows[lanes[1]].append(("chr1", 5000, 5100, "TTTTTTTTTTTTTTTT-1", 1))
    rna[lanes[1]].append(("TTTTTTTTTTTTTTTT-1", 999, 0, "k562"))
    man = []
    for lane in lanes:
        fr = pd.DataFrame(frag_rows[lane]).sample(frac=1.0, random_state=1)  # unsorted input
        fp = d / f"{lane}_fragments.tsv.gz"
        fr.to_csv(fp, sep="\t", header=False, index=False, compression="gzip")
        n = len(rna[lane])
        X = np.zeros((n, len(genes)), dtype=np.float32)
        for i, (b, umi, j, ct) in enumerate(rna[lane]):
            w = np.array([30, 20, 10, 10, 3, 5, 2, 8 if j != 4 else 40, 10, 2], dtype=float)
            if ct == "hepg2":
                w[0], w[9] = 2, 30
            X[i] = rng.multinomial(umi, w / w.sum())
        a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[r[0] for r in rna[lane]]),
                       var=pd.DataFrame({"gene_symbol": [g[2] for g in genes]}, index=[g[0] for g in genes]))
        hp = d / f"{lane}_rna.h5ad"
        a.write_h5ad(hp)
        man.append({"lane": lane, "analysis_accession": acc[lane], "fragments": str(fp), "rna_matrix": str(hp)})
    mp = d / "manifest.tsv"
    pd.DataFrame(man).to_csv(mp, sep="\t", index=False)
    amp = d / "atac_to_rna.tsv"
    amp.write_text("".join(f"{a}\t{r}\n" for a, r in atac_map.items()))
    cdf = pd.DataFrame(cells)
    ap = d / "cell_annotations.tsv"
    cdf[["barcode", "lane", "cell_type", "subsample"]].to_csv(ap, sep="\t", index=False)
    return {"cells": cdf, "manifest": mp, "annotations": ap, "atac_map": amp, "gtf": gtf, "tss": tss, "peaks": peaks,
            "lanes": lanes, "subs": subs, "acc": acc, "genes": genes, "shared": shared_seq}


def cmd_selftest(args) -> int:
    global OUT_ROOT
    pd, np, ad = _pd(), _np(), _ad()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def raises(fn, exc=(SystemExit, ValueError)):
        try:
            fn()
        except exc as e:  # noqa: PERF203
            return str(e) or True
        return False
    saved = OUT_ROOT
    _QUIET.update(record=False, log=False)
    plots = not args.no_plots
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        OUT_ROOT = d / "runs"
        try:
            W = _synthetic(d)
            cells = W["cells"]
            print("\nbarcodes")
            check(full_barcode("ACGTACGTACGTACGT-1", "IGVFSM0001LANA") == "ACGTACGTACGTACGT_IGVFSM0001LANA",
                  "full_barcode: GEM -1 stripped, lane accession appended")
            check(full_barcode("AAAA", "L", {"AAAA": "CCCC"}) == "CCCC_L" and full_barcode("AAAA_SUB1", "L", {"AAAA": "CCCC"}) == "CCCC_SUB1",
                  "full_barcode: ATAC->RNA translation; an existing subpool suffix is kept")
            check(len(embedded_chrom_sizes()) == 195 and embedded_chrom_sizes()[24] == ("chrM", 16569) and embedded_chrom_sizes()[-1][0] == "chrEBV",
                  "embedded IGVF DACC GRCh38 chrom sizes: 195 contigs, chr1..chrM first, chrEBV last")

            print("\nbuild-pseudobulks")
            rc = main(["build-pseudobulks", "--annotations", str(W["annotations"]), "--manifest", str(W["manifest"]),
                       "--dataset", "igvf1", "--atac-rna-map", str(W["atac_map"]), "--tss", str(W["tss"]),
                       "--peaks", str(W["peaks"]), "--out-dir", str(d / "build")])
            pbs = d / "build" / "igvf1" / "pseudobulks"
            dirs = sorted(x.name for x in pbs.iterdir())
            check(rc == 0 and dirs == [f"annotation-{ct}-{s}" for ct in ("hepg2", "k562") for s in W["subs"]],
                  f"4 primary pseudobulk directories annotation-{{cell_type}}-{{subsample}} ({len(dirs)})")
            q = pd.read_csv(pbs / f"annotation-k562-{W['subs'][0]}" / "per_cell_qc.tsv", sep="\t")
            check(list(q.columns) == PER_CELL_QC_COLS, "per_cell_qc.tsv has the 12 upstream columns in order")
            check(q["barcode"].map(lambda b: bool(BARCODE_SPEC_RE.match(b))).all() and len(q) == 12,
                  "barcodes are {16-bp}_{IGVF lane accession}; 6 cells x 2 lanes")
            sh = q[q["barcode"].str.startswith(W["shared"])]
            check(len(sh) == 2 and sh["barcode"].nunique() == 2, "the same 16-bp barcode in two lanes stays two cells (lane suffix)")
            check(set(q["analysis_accession"]) == set(W["acc"].values()), "analysis_accession from the manifest lane")
            good = q[q["barcode"].str.contains("_IGVFSM0001LANA")]
            exp_nf = [1501, 1510, 1520, 1530, 1540, 50]
            got_nf = good["num_frags"].tolist()
            check(sorted(got_nf) == sorted(exp_nf), f"lane 1 fragments reach their cells through the ATAC->RNA map ({got_nf}); chrFOO counted")
            check(np.allclose(q["pct_duplicated_reads"], 100 * (np.ceil(q["num_frags"] / 5)) / (q["num_frags"] + np.ceil(q["num_frags"] / 5)), atol=0.2),
                  "pct_duplicated_reads = 100 x (sum(col5) - n) / sum(col5)")
            lo = q[q["num_frags"] == 50]
            check((q.loc[q["num_frags"] > 1000, "tss_enrichment"] > 5).all() and (lo["tss_enrichment"] < 3).all(),
                  "tss_enrichment: planted TSS-proximal cells high, low-ATAC cells low")
            check((q.loc[q["num_frags"] > 1000, "frip"].between(0.5, 0.8)).all() and (lo["frip"] < 0.1).all(),
                  "frip: fraction of fragments in peaks recovers the planted 2/3")
            check(np.allclose(q.loc[q["num_frags"] > 1000, "nucleosomal_signal"], 1 / 3, atol=0.02),
                  "nucleosomal_signal = mono (147-293 bp) / nucleosome-free (< 147 bp) = 1/3 planted")
            check(q["rna_read_count"].min() == 60 and (q["gene_count"] <= 10).all() and q["pct_mito"].max() > 20,
                  "RNA metrics: UMIs, genes > 0, % mito planted high in one cell")
            a = ad.read_h5ad(pbs / f"annotation-k562-{W['subs'][0]}" / "rna_counts_mtx.h5ad")
            check(a.n_obs == 12 and "TTTTTTTTTTTTTTTT_IGVFSM0002LANB" not in a.obs_names, "rna_counts_mtx.h5ad: annotated cells only")
            ff = pd.read_csv(pbs / f"annotation-k562-{W['subs'][0]}" / "fragments.tsv.gz", sep="\t", header=None)
            check(ff.shape[1] == 5 and set(ff[3]) <= set(q["barcode"]), "fragments.tsv.gz: 5 columns, full barcodes")

            print("\nqc-datatable")
            dt = d / "k562_per_cell_qc.tsv"
            r = build_qc_datatable(pbs, "k562", dt)
            check(r["cells"] == 24 and r["directories"] == 2, "concatenates both k562 subsample directories (24 cells)")
            (pbs / "annotation-k562-IGVFSM9999EMPT").mkdir()
            (pbs / "annotation-k562-IGVFSM9999EMPT" / "per_cell_qc.tsv").write_text("")
            r2 = build_qc_datatable(pbs, "k562", d / "k2.tsv")
            check(r2["cells"] == 24 and len(r2["skipped"]) == 1, "0-byte per_cell_qc.tsv skipped")
            (pbs / "annotation-k562-IGVFSM9999EMPT" / "per_cell_qc.tsv").write_text("barcode\tother\nX\t1\n")
            check(bool(raises(lambda: build_qc_datatable(pbs, "k562", d / "k3.tsv"))), "column mismatch is an error")
            shutil.rmtree(pbs / "annotation-k562-IGVFSM9999EMPT")

            print("\nexplore")
            t = pd.DataFrame({"analysis_accession": "A", "barcode": [f"b{i}" for i in range(6)], "subsample": ["s1"] * 3 + ["s2"] * 3,
                              "rna_read_count": [1000, 999, 5000, 5000, 5000, 5000], "gene_count": [2000] * 6,
                              "pct_mito": [1, 1, 40, 1, 1, 1], "pct_ribo": [1] * 6, "num_frags": [5000, 5000, 5000, 10, 5000, 5000],
                              "pct_duplicated_reads": [10] * 6, "nucleosomal_signal": [0.5] * 6,
                              "tss_enrichment": [10, 10, 10, 1, 3, 10], "frip": [0.5] * 6})
            tp = d / "toy.tsv"
            t.to_csv(tp, sep="\t", index=False)
            tm = read_meta(tp)
            _, tot, alo = failure_stats(tm, EXPLORE_DEFAULTS)
            check(tot["rna_min"] == 1 and tot["atac_min"] == 1 and tot["tss_enr_min"] == 1 and tot["pct_mt_max"] == 1,
                  "Total: cells failing each threshold (999 < 1000; 10 frags; TSS 1; 40% mito)")
            check(alo["atac_min"] == 0 and alo["tss_enr_min"] == 0 and alo["rna_min"] == 1,
                  "Alone: the cell failing both ATAC and TSS counts for neither")
            check(int(pass_mask(tm, EXPLORE_DEFAULTS, THRESH_ORDER, "explore").sum()) == 3,
                  "explore keeps rna = 1000 and TSS = 3 (>=): 3 of 6 pass")
            check(label_threshold_set(parse_threshold_string("--tss-min 5 --pct-mt-max 20", EXPLORE_DEFAULTS)) ==
                  "% mito < 20, TSS enrichment >= 5" and label_threshold_set(EXPLORE_DEFAULTS) == "defaults", "set labels as upstream")
            check(bool(raises(lambda: parse_threshold_string("--bogus 1", EXPLORE_DEFAULTS))), "unknown flag rejected")
            rc = main(["explore", "--meta", str(tp), "--sets", "--tss-min 3; --tss-min 5", "--show-subsamples", "--out-dir", str(d / "ex")])
            br = pd.read_csv(d / "ex" / "threshold_breakdown.tsv", sep="\t")
            ps = pd.read_csv(d / "ex" / "per_subsample.tsv", sep="\t")
            txt = (d / "ex" / "explore_report.txt").read_text()
            check(rc == 0 and set(br["set"]) == {1, 2} and int(br[(br.set == 2) & (br.threshold == "tss_enr_min")]["total"].iloc[0]) == 2,
                  "--sets 'a; b': two threshold sets; TSS >= 5 drops 2")
            check(len(ps) == 4 and "Most stringent threshold" in txt and "Per-subsample breakdown" in txt, "per-subsample breakdown written")

            print("\nqc-filter")
            rc = main(["qc-filter", "--meta", str(tp), "--out-dir", str(d / "qf"), "--no-plots"])
            gd = pd.read_csv(d / "qf" / "filtered_barcodes_with_subsamples.tsv.gz", sep="\t")
            check(rc == 0 and list(gd.columns) == GUIDE_COLS and list(gd["barcode"]) == ["b5"],
                  "strict inequalities (plot_per_cell_qc.R): rna = 1000 and TSS = 3 dropped; guide has barcode, subsample, analysis_accession")
            th = (d / "qf" / "qc_thresholds.tsv").read_text().splitlines()
            check(th[0] == "threshold\tvalue" and th[1] == "Minimum RNA reads per cell\t1000" and th[2] == "Maximum RNA reads per cell\tInf"
                  and th[3] == "Minimum genes per cell\t0" and th[10] == "Maximum nucleosomal signal\t1.5", "qc_thresholds.tsv: 12 rows, R number format")
            mt = pd.read_csv(d / "qf" / "filtered_cell_subsample_metrics.tsv", sep="\t")
            check(list(mt.columns) == ["subsample", "n_cells", "total_fragments", "total_RNA_reads", "mean_frag_per_cell",
                                       "mean_RNA_per_cell", "mean_frip", "mean_tss"] and mt["n_cells"].tolist() == [1],
                  "filtered_cell_subsample_metrics.tsv columns")
            t2 = t.copy()
            t2["frip"] = np.nan
            t2.to_csv(d / "toy_na.tsv", sep="\t", index=False, na_rep="NA")
            r_na = run_qc_filter(d / "toy_na.tsv", d / "qf_na", dict(FILTER_DEFAULTS), plots=False)
            r_up = run_qc_filter(d / "toy_na.tsv", d / "qf_up", dict(FILTER_DEFAULTS), plots=False, upstream_compat=True)
            check(r_na["kept"] == 1 and r_na["skipped_thresholds"] == ["frip"] and r_up["kept"] == 0,
                  "all-NA metric skipped (1 kept); --upstream-compat drops every cell as the R filter would")
            gdir = d / "plots" / "igvf1" / "k562"
            r = run_qc_filter(dt, gdir, dict(FILTER_DEFAULTS, gene_min=5, rna_min=500, tss_enr_min=3), plots=plots)
            planted_bad = set(cells[(cells.cell_type == "k562") & (cells.bad != "")]["full"])
            check(r["kept"] == 16 and not (set(r["guide"]["barcode"]) & planted_bad),
                  f"synthetic k562: the 8 planted low-RNA / low-ATAC cells are removed ({r['kept']} of 24 kept)")
            if plots:
                check(all(Path(f).stat().st_size > 5000 for f in r["figures"]) and len(r["figures"]) == 5,
                      "RNA / ATAC QC figures and the three cells-per-subsample bars written")
            guide = gdir / "filtered_barcodes_with_subsamples.tsv.gz"

            print("\nfilter-atac")
            atac_out = d / "multi" / "atac_fragments_igvf1_k562.tsv.gz"
            ra = filter_atac(guide, pbs, "k562", None, atac_out)
            fo = pd.read_csv(ra["out"], sep="\t", header=None, dtype={0: str})
            order = {c: i for i, c in enumerate(c for c, _ in embedded_chrom_sizes())}
            rk = fo[0].map(order).to_numpy()
            srt = (np.diff(rk) >= 0).all() and all((np.diff(g[1].to_numpy()) >= 0).all() for _, g in fo.groupby(0))
            check(srt and set(fo[3]) == set(r["guide"]["barcode"]), "output in chrom-sizes order, start-sorted, exactly the guide's barcodes")
            check(ra["chromosomes_dropped"].get("chrFOO", 0) > 0 and "chrFOO" not in set(fo[0]), "chromosome absent from chrom sizes dropped")
            check(ra["index"] is None or Path(ra["index"]).exists(), f"bgzip + tabix index ({ra['index_method']})")
            pcq = pd.read_csv(dt, sep="\t").set_index("barcode")
            n_exp = int(pcq.loc[list(r["guide"]["barcode"]), "num_frags"].sum()) - int(ra["chromosomes_dropped"].get("chrFOO", 0))
            check(ra["rows_written"] == n_exp, f"every fragment of every passing cell kept ({ra['rows_written']})")
            bad_guide = d / "bad_guide.tsv"
            bad_guide.write_text("barcode\tsubsample\tanalysis_accession\nAAAAAAAAAAAAAAAA_IGVFSM0001LANA\tS\tA\n")
            check("not found in any fragment file" in str(raises(lambda: filter_atac(bad_guide, pbs, "k562", None, d / "x.tsv.gz"))),
                  "a QC-guide barcode missing from every fragment file is an error")
            mal = d / "mal" / "annotation-zz-IGVFSM0000ZZZZ"
            mal.mkdir(parents=True)
            with gzip.open(mal / "fragments.tsv.gz", "wt") as fh:
                fh.write("chr1\t1\t2\tA_B\n")
            check("expected 5 tab-separated fields" in str(raises(lambda: filter_atac(bad_guide, d / "mal", "zz", None, d / "y.tsv.gz"))),
                  "a fragment record without 5 fields is an error")
            rcl = filter_atac(guide, pbs, "k562", None, d / "clean.tsv.gz", clean=True)
            fcl = pd.read_csv(rcl["out"], sep="\t", header=None)
            check(fcl[3].str.len().eq(16).all(), "--clean writes the 16-bp barcode only (README flag)")
            both = filter_atac(guide, pbs, "k562,hepg2", None, d / "merged.tsv.gz")
            check(len(both["directories"]) == 4, "comma-separated cell types merge directory sets")

            print("\nfilter-rna")
            rna_out = d / "multi" / "rna_count_matrix_igvf1_k562.mtx"
            rr = filter_rna(str(guide), pbs, "k562", str(rna_out), gtf=str(W["gtf"]), standard_only=True)
            md = d / "multi" / "rna_count_matrix_igvf1_k562"
            import scipy.io  # type: ignore
            head = gzip.open(md / "matrix.mtx.gz", "rt").readline().strip()
            M = scipy.io.mmread(str(md / "matrix.mtx.gz")).tocsr()
            feats = gzip.open(md / "features.tsv.gz", "rt").read().split()
            bcs = gzip.open(md / "barcodes.tsv.gz", "rt").read().split()
            check(head == "%%MatrixMarket matrix coordinate integer general" and M.shape == (len(feats), len(bcs)),
                  "matrix.mtx.gz: integer MatrixMarket, genes x cells")
            check(feats == sorted(feats) and "SCAFONLY" not in feats and feats.count("DUPSYM") == 1 and "SHARED" in feats,
                  f"features: symbols sorted, scaffold-only gene dropped, duplicate symbol collapsed ({len(feats)} symbols)")
            src = ad.concat([ad.read_h5ad(pbs / f"annotation-k562-{s}" / "rna_counts_mtx.h5ad") for s in W["subs"]])
            src = src[bcs]
            ids = list(src.var_names)
            dup_exp = np.asarray(src.X[:, [ids.index("ENSG00000000003.2"), ids.index("ENSG00000000004.1")]].sum(axis=1)).ravel()
            sh_exp = np.asarray(src.X[:, ids.index("ENSG00000000006.1")].todense()).ravel()
            check(np.allclose(M[feats.index("DUPSYM")].toarray().ravel(), dup_exp) and np.allclose(M[feats.index("SHARED")].toarray().ravel(), sh_exp),
                  "DUPSYM = sum of its two Ensembl IDs; SHARED keeps only its standard-chromosome counts")
            check(set(bcs) == set(r["guide"]["barcode"]) and rr["cells"] == 16, "barcodes.tsv = the QC guide's barcode set")
            lg = Path(rr["log"]).read_text()
            check("Gene symbols ONLY on nonstandard chromosomes (removed entirely): 1" in lg and "SHARED" in lg, "GTF mapping log written")
            bad_gtf = d / "bad.gtf"
            bad_gtf.write_text("\n".join(l for l in W["gtf"].read_text().splitlines() if "ENSG00000000010.1" not in l) + "\n")
            check("not found in GTF" in str(raises(lambda: filter_rna(str(guide), pbs, "k562", str(d / "z.mtx"), gtf=str(bad_gtf)))),
                  "an Ensembl ID missing from the GTF is a hard failure")
            check(bool(raises(lambda: filter_rna(str(guide), pbs, "k562", str(d / "z.mtx"), standard_only=True))),
                  "--standard-chromosomes-only without --gtf rejected")
            re_ = filter_rna(str(guide), pbs, "k562", str(d / "ens.h5ad"), ensembl_ids=True)
            check(re_["genes"] == 10 and ad.read_h5ad(d / "ens.h5ad").var_names[0].startswith("ENSG"), "--ensemblIDs-as-genes keeps IDs (.h5ad)")
            rs = filter_rna(str(guide), pbs, "k562", str(d / "sym.csv.gz"))
            cs = pd.read_csv(d / "sym.csv.gz", index_col=0)
            check(rs["genes"] == 8 and cs.index.name == "barcode" and cs.shape == (16, 8),
                  "no GTF: var gene_symbol collapse (SCAFONLY kept, 8 symbols); .csv.gz cells x genes")
            ru = filter_rna(str(guide), pbs, "k562", str(d / "real.mtx"), gtf=str(W["gtf"]), standard_only=True, integer=False)
            check(gzip.open(d / "real" / "matrix.mtx.gz", "rt").readline().strip().endswith("real general"),
                  "--upstream-compat writes the real field as scipy mmwrite of the float matrix does")

            print("\npackage-rna + validate")
            tar = package_rna(md)
            with tarfile.open(tar) as tf:
                names = sorted(tf.getnames())
                mm = tf.extractfile("matrix.mtx").read(40).decode()
            check(names == ["barcodes.tsv", "features.tsv", "matrix.mtx"] and mm.startswith("%%MatrixMarket"), ".tar.gz holds decompressed members")
            vc = validate_outputs(guide, md, ra["out"])
            check(all(c["ok"] for c in vc) and len(vc) >= 14, f"all {len(vc)} spec checks pass on the outputs")
            dg = d / "dup_guide.tsv.gz"
            with gzip.open(dg, "wt") as fh:
                fh.write("barcode\tsubsample\tanalysis_accession\nAAAAAAAAAAAAAAAA_IGVFSM0001LANA\tS\tA\nAAAAAAAAAAAAAAAA_IGVFSM0001LANA\tS\tA\nACGT\tS\tA\n")
            vb = validate_outputs(dg)
            check(sum(not c["ok"] for c in vb) == 2, "duplicate and malformed barcodes caught by validate")
            rc = main(["validate", "--qc-guide", str(guide), "--rna-matrix", str(md), "--fragments", ra["out"], "--out-dir", str(d / "val")])
            check(rc == 0, "validate subcommand exits 0")

            print("\nconfig-table + run (Snakefile)")
            ctab = config_table(d / "multi", "igvf1", ["k562", "hepg2"])
            check(list(ctab.columns) == CONFIG_TABLE_COLS and ctab["model_dir"].iloc[0] == SCE2G_MODEL_DIR and
                  config_table(d, "x", ["a"], True)["model_dir"].iloc[0] == UPSTREAM_MODEL_DIR, "config table columns; model_dir sce2g names / upstream paths")
            check(determine_mem_mb(100, True) == 8000 and determine_mem_mb(1000, True) == 32000 and determine_mem_mb(1000, True, 2) == 64000
                  and determine_mem_mb(1e5, True) == 250000, "determine_mem_mb as the Snakefile")
            cfgp = d / "config.yaml"
            cfgp.write_text(f'# test\nQC_plots_dir: "{d / "plots"}"\npseudobulk_dir: "{d / "build"}"\nout_dir: "{d / "mdata"}"\n'
                            f'transcriptome: "{W["gtf"]}"\nchrom_sizes: "reference/IGVF.DACC.GRCh38.chrom.sizes.tsv"\n'
                            f'datasets:\n  igvf1:\n    - k562\n    - hepg2\n')
            pc = parse_simple_yaml(cfgp.read_text())
            check(pc["datasets"] == {"igvf1": ["k562", "hepg2"]} and pc["out_dir"] == str(d / "mdata"), "config_QC_pseudobulks.yaml parsed")
            rc = main(["run", "--config", str(cfgp), "--dry-run", "--auto-qc", "--out-dir", str(d / "rdry")])
            sj = json.loads((d / "rdry" / "summary.json").read_text())
            check(rc == 0 and [j["rule"] for j in sj["jobs"]].count("atac_fragment_file") == 2 and
                  [j["rule"] for j in sj["jobs"]].count("qc_filter") == 1 and not (d / "mdata").exists(),
                  "dry run lists jobs (hepg2 needs auto QC) and writes nothing")
            rc = main(["run", "--config", str(cfgp), "--auto-qc", "--gene-min", "5", "--rna-min", "500", "--no-plots", "--out-dir", str(d / "rrun")])
            ctp = d / "mdata" / "config" / "tables" / "igvf1_config.tsv"
            cfg_t = pd.read_csv(ctp, sep="\t", keep_default_na=False)
            check(rc == 0 and list(cfg_t["cluster"]) == ["k562", "hepg2"] and all(Path(x).exists() for x in cfg_t["atac_frag_file"])
                  and all((Path(x) / "matrix.mtx.gz").exists() for x in cfg_t["rna_matrix_file"]),
                  "run: per-cluster fragments + RNA matrix directory + {dataset}_config.tsv")
            check(Path(str(cfg_t["atac_frag_file"].iloc[0]) + ".tbi").exists() or ra["index"] is None, "run: fragments tabix-indexed")
            hg = pd.read_csv(d / "plots" / "igvf1" / "hepg2" / "filtered_barcodes_with_subsamples.tsv.gz", sep="\t")
            check(len(hg) == 16, "--auto-qc wrote the hepg2 QC guide (16 of 24)")

            print("\nsce2g-prep (chain to sce2g-pipeline)")
            try:
                s2 = _sce2g()
                pr = sce2g_prep(ctp, d / "s2")
                rows = s2._cluster_rows(argparse.Namespace(cluster_config=pr["cell_cluster_config"], models=["multiome_powerlaw_v3"]))
                c0 = pr["clusters"][0]
                check(abs(c0["tpm_sum"] - 1e6) < 1e-3 and c0["rna_cells"] == 16 and Path(c0["tagAlign"]).exists(),
                      "tagAlign + RNA pseudobulk TPM (sums to 1e6) via sce2g_pipeline_skill")
                check(c0["fragment_count"] == len(pd.read_csv(cfg_t["atac_frag_file"].iloc[0], sep="\t", header=None)) and c0["atac_barcodes"] == 16,
                      "fragment_count / cell barcodes match the principal pseudobulk")
                check([r_["cluster"] for r_ in rows] == ["k562", "hepg2"] and rows[0]["models"] == ["multiome_powerlaw_v3", "scATAC_powerlaw_v3"]
                      and all(s2.resolve_model(m)["name"] == m for m in rows[0]["models"]),
                      "config_cell_clusters.tsv parses with sce2g-pipeline's --cluster-config reader; model names resolve")
                M2, g2, c2 = s2.read_rna_matrix(cfg_t["rna_matrix_file"].iloc[0])
                check(M2.shape == (len(g2), 16), "sce2g-pipeline reads the RNA matrix directory (genes x cells)")
            except ImportError as e:
                print(f"  skip  sce2g_pipeline_skill not importable: {e}")

            print("\nfetch (processed-first, fixture portal)")
            pl = _lineage()
            base, info = pl.fixture_portal()
            prin = info["ids"]["prin"]
            prin_obj = {"@id": prin, "@type": ["AnalysisSet", "FileSet", "Item"], "accession": "TSTDS0006PRN", "status": "released",
                        "file_set_type": "principal analysis", "input_file_sets": [info["ids"]["lab"], info["ids"]["uni"]],
                        "files": [{"@id": "/tabular-files/TSTFI0701ANN/", "accession": "TSTFI0701ANN", "content_type": "cell annotations",
                                   "file_format": "tsv", "file_size": 2_000_000, "controlled_access": False, "status": "released",
                                   "href": "/tabular-files/TSTFI0701ANN/@@download/TSTFI0701ANN.tsv"}]}

            def fx(url):
                if url.split("?")[0] == prin:
                    return 200, prin_obj
                st_, obj = base(url)
                if st_ == 200 and isinstance(obj, dict) and obj.get("@id") in ("/tabular-files/TSTFI0103FRG/", "/matrix-files/TSTFI0101H5A/"):
                    obj = dict(obj, href=obj["@id"] + "@@download/" + obj["accession"] + (".bed.gz" if "FRG" in obj["accession"] else ".h5ad"))
                if st_ == 200 and isinstance(obj, dict) and obj.get("files"):
                    obj = dict(obj, files=[dict(f, href=f.get("href") or (f["@id"] + "@@download/" + f["accession"]))
                                           if isinstance(f, dict) else f for f in obj["files"]])
                return st_, obj
            got = []

            def fake_dl(href, target):
                got.append(href)
                Path(target).write_text("x")
                return target
            fr_ = fetch_inputs("TSTDS0001RNA", d / "dl", fetch=fx, downloader=fake_dl, max_gb=5.0, dry_run=True)
            man = fr_["manifest"]
            check(list(man["lane"]) == ["TSTSM0001MUX"] and man["fragments_accession"].iloc[0] == "TSTFI0103FRG"
                  and man["rna_matrix_accession"].iloc[0] == "TSTFI0101H5A",
                  "uniform-pipeline lane: fragments + h5ad chosen (lab hdf5 and 13 GB tar not), lane = multiplexed sample")
            check(fr_["selection"]["annotations"] and fr_["selection"]["annotations"][0]["accession"] == "TSTFI0701ANN",
                  "cell annotations taken from the principal analysis set")
            check(not got and any("exceeds the remaining budget" in s["why"] for s in fr_["skipped"]),
                  "--dry-run downloads nothing; files over the --max-gb budget are skipped")
            fr2 = fetch_inputs("TSTDS0001RNA", d / "dl", fetch=fx, downloader=fake_dl, max_gb=10.0)
            check(len(got) == 3 and fr2["annotations_path"] and Path(fr2["manifest"]["fragments"].iloc[0]).exists(),
                  "within budget: annotations, fragments and h5ad downloaded through portal_download")
        finally:
            OUT_ROOT = saved
            _QUIET.update(record=True, log=True)
    n_fail = sum(1 for ok, _ in checks if not ok)
    print(f"\n{len(checks) - n_fail}/{len(checks)} checks pass")
    if n_fail:
        print("selftest: FAILED")
        return 1
    print("selftest: all checks pass")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent principal-pseudobulks", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, label=True):
        if label:
            p.add_argument("--label", default=None, help="run-directory label")
        p.add_argument("--out-dir", default=None, help="output directory (default Docs/PrincipalPseudobulks/<ts>_<label>)")

    p = sub.add_parser("fetch", help="processed-first Portal inputs: principal analysis cell annotations + per-lane fragments / h5ad")
    p.add_argument("--accession", required=True, help="any IGVF accession (principal / intermediate analysis set, measurement set, sample)")
    p.add_argument("--dest", help="download directory (default Data/PrincipalPseudobulks/<accession>)")
    p.add_argument("--max-gb", type=float, default=20.0, help="total download budget in GB")
    p.add_argument("--dry-run", action="store_true", help="plan and list only")
    p.add_argument("--include-lab", action="store_true", help="also use lab (non-uniform) intermediate analysis sets")
    p.add_argument("--depth", type=int, default=6)
    common(p)
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("build-pseudobulks", help="annotations + lane fragments / h5ad -> primary pseudobulk directories")
    p.add_argument("--annotations", required=True, help="cell annotation table (.tsv/.csv/.h5ad obs): barcode + cell type (+ subsample, lane)")
    p.add_argument("--manifest", help="TSV lane, fragments, rna_matrix[, analysis_accession] (from fetch)")
    p.add_argument("--lane", nargs="+", help="lane accessions (with --fragments / --rna-matrix, same order)")
    p.add_argument("--fragments", nargs="+")
    p.add_argument("--rna-matrix", nargs="+")
    p.add_argument("--dataset", required=True)
    p.add_argument("--annotation-column")
    p.add_argument("--barcode-column")
    p.add_argument("--subsample-column")
    p.add_argument("--lane-column")
    p.add_argument("--analysis-column")
    p.add_argument("--atac-rna-map", help="two-column ATAC barcode -> RNA barcode table (10x multiome whitelist pairs)")
    p.add_argument("--tss", help="TSS bed (for tss_enrichment)")
    p.add_argument("--peaks", help="peak bed (for frip)")
    p.add_argument("--gtf", help="recorded in the generated run config as transcriptome")
    p.add_argument("--min-cells", type=int, default=1, help="skip cell type x subsample groups with fewer cells")
    p.add_argument("--pseudobulk-root", help="write {dataset}/pseudobulks here instead of the run directory")
    common(p)
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("qc-datatable", help="Step 0: concatenate per_cell_qc.tsv of a cell type's subsample directories")
    p.add_argument("--pseudobulks", required=True)
    p.add_argument("--cell-type", required=True)
    p.add_argument("--out")
    common(p)
    p.set_defaults(func=cmd_datatable)

    p = sub.add_parser("explore", help="Step 1: Total / Alone cells dropped per threshold (explore_qc_thresholds.R)")
    p.add_argument("--meta", required=True, help="per-cell QC datatable")
    p.add_argument("--sets", nargs="+", help='threshold sets, e.g. "--tss-min 3; --tss-min 5 --frip-min 0.1"')
    p.add_argument("--show-subsamples", action="store_true")
    p.add_argument("--upstream-compat", action="store_true", help="do not skip all-NA metrics")
    _threshold_args(p, EXPLORE_DEFAULTS)
    common(p)
    p.set_defaults(func=cmd_explore)

    p = sub.add_parser("qc-filter", help="Step 1: final thresholds -> QC guide, thresholds, metrics, QC figures (plot_per_cell_qc.R)")
    p.add_argument("--meta", required=True)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--upstream-compat", action="store_true", help="do not skip all-NA metrics")
    _threshold_args(p, FILTER_DEFAULTS)
    common(p)
    p.set_defaults(func=cmd_qc_filter)

    p = sub.add_parser("filter-atac", help="Step 2: QC-guide barcodes -> sorted, bgzipped, tabix-indexed fragments")
    p.add_argument("--qc-guide", required=True)
    p.add_argument("--pseudobulks", required=True)
    p.add_argument("--cell-type", required=True, help="comma-separated to merge clusters")
    p.add_argument("--chrom-sizes", help="sort order (default: embedded IGVF DACC GRCh38 file)")
    p.add_argument("--out")
    p.add_argument("--clean", action="store_true", help="write only the 16-bp barcode")
    common(p)
    p.set_defaults(func=cmd_filter_atac)

    p = sub.add_parser("filter-rna", help="Step 2: QC-guide barcodes -> gene-symbol RNA count matrix")
    p.add_argument("--qc-guide", required=True)
    p.add_argument("--pseudobulks", required=True)
    p.add_argument("--cell-type", required=True)
    p.add_argument("--out", help=".mtx (directory), .csv(.gz) or .h5ad")
    p.add_argument("--gtf", help=f"GTF for Ensembl -> symbol ({IGVF_TRANSCRIPTOME})")
    p.add_argument("--ensembl-ids-as-genes", "--ensemblIDs-as-genes", dest="ensembl_ids_as_genes", action="store_true")
    p.add_argument("--standard-chromosomes-only", action="store_true")
    p.add_argument("--log")
    p.add_argument("--upstream-compat", action="store_true", help="write the MatrixMarket real field")
    common(p)
    p.set_defaults(func=cmd_filter_rna)

    p = sub.add_parser("package-rna", help="RNA matrix directory -> spec .tar.gz (decompressed members)")
    p.add_argument("--rna-matrix", required=True)
    p.add_argument("--out")
    common(p)
    p.set_defaults(func=cmd_package)

    p = sub.add_parser("config-table", help="per-dataset cluster table for scE2G (rule make_config_table)")
    p.add_argument("--out-root", required=True, help="the run's out_dir")
    p.add_argument("--dataset", required=True)
    p.add_argument("--cell-types", nargs="+", required=True)
    p.add_argument("--out")
    p.add_argument("--upstream-compat", action="store_true", help="model_dir as models/... paths")
    common(p)
    p.set_defaults(func=cmd_config_table)

    p = sub.add_parser("run", help="the Snakefile: every (dataset, cell_type) -> fragments + RNA matrix + config tables")
    p.add_argument("--config", help="config_QC_pseudobulks.yaml / .json")
    p.add_argument("--qc-plots-dir")
    p.add_argument("--pseudobulk-dir")
    p.add_argument("--out-root", help="out_dir of the config")
    p.add_argument("--transcriptome", help="GTF")
    p.add_argument("--chrom-sizes")
    p.add_argument("--datasets", nargs="+", help="dataset:ct1,ct2")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--auto-qc", action="store_true", help="build the datatable and run qc-filter where the QC guide is missing")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--upstream-compat", action="store_true")
    _threshold_args(p, FILTER_DEFAULTS)
    common(p)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("validate", help="check outputs against the two file specs")
    p.add_argument("--qc-guide")
    p.add_argument("--rna-matrix")
    p.add_argument("--fragments")
    p.add_argument("--chrom-sizes")
    common(p)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("sce2g-prep", help="config table -> tagAlign, RNA pseudobulk TPM, config_cell_clusters.tsv for sce2g-pipeline")
    p.add_argument("--config-table", required=True)
    p.add_argument("--chrom-sizes")
    common(p)
    p.set_defaults(func=cmd_sce2g_prep)

    p = sub.add_parser("selftest", help="synthetic lanes / cells / genes; every subcommand asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    if args.cmd != "selftest" and _QUIET["log"]:
        print(f"Log: {setup_logging()}")
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
