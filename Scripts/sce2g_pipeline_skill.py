#!/usr/bin/env python3
"""scE2G single-cell enhancer-gene pipeline rewritten in Python (port of EngreitzLab/scE2G).

Port of https://github.com/EngreitzLab/scE2G (MIT, "Engreitz and Andersson Labs"; Snakemake + R +
Python) at commit 7cb2af750fb96006f5d2b7c5475dcff30ea0e6c9 (2026-08-13, tags v1.0/v1.2/v1.3 are
ancestors), together with the pieces of its ENCODE_rE2G submodule (pinned 5930e746, MIT) that scE2G
calls to assemble feature tables, apply and threshold models.  Relationship: port -- every rule and
script was read and the definitions re-derived in numpy/pandas/scipy; no code was copied.  The
existing `sce2g` workbench wraps an upstream checkout and `sce2g-predict` is a clean-room Kendall
toy; this module is the pipeline itself, runnable with no Snakemake, R, Signac, bedtools or
fast_kendall_sc.  The four published v3 models ship inside the module: their feature tables,
score / TPM thresholds and the logistic-regression weights read out of model.pkl and the 23
held-out-chromosome cv_models/*.pkl.  The per-model quantile-normalisation references are data,
fetched by `setup` into Data/scE2G/models/<model>/qnorm_reference.tsv.gz.

Definitions, as upstream computes them
  tagAlign            each fragment (chr, s, e) -> two tags (chr, s, m, N, 1000, +) and
                      (chr, m+1, e, N, 1000, -), m = int((s+e)/2); fragments on chromosomes absent
                      from the chrom-sizes file are dropped; tags sorted in chrom-sizes order, then
                      start, end; bgzip + tabix (pysam).  fragment_count = filtered fragment lines;
                      cell_barcodes = unique barcodes in file order.
  ATAC bigWig         bedtools genomecov -bg of fragments (half-open coverage), optionally scaled by
                      1e6 / fragment_count (ATAC_norm); bigWig when pyBigWig / bedGraphToBigWig exist.
  Kendall pairs       narrowPeak peaks x EnhancerPredictionsAllPutative rows overlapping them
                      (GRanges closed-interval overlap on the numbers as written); PeakName
                      chr-start-end, PairName PeakName_TargetGene, deduplicated.
  ATAC matrix         peak x cell counts of fragments overlapping each unique peak (Signac
                      FeatureMatrix via tabix: fragment [s+1, e] 1-based vs peak [start, end]);
                      cells = RNA-matrix cells (intersected with fragment barcodes when the RNA
                      matrix is not pre-filtered), at most max_cell_count (20,000) sampled.
  RNA features        on RNA counts restricted to those cells: RnaPseudobulkTPM = rowSum / total x
                      1e6; RnaDetectedPercent = fraction of cells with count > 0;
                      mean_log_normalized_rna = rowMean(log1p(count / cell total x 1e4)) (Seurat
                      LogNormalize).  Written as RNA_pseudobulkTPM / RNA_meanLogNorm /
                      RNA_percentCellsDetected.  Gene names are mapped from the GTF gene_name to
                      the TSS500 reference name through a 1:1 Ensembl-ID key (version stripped).
  Kendall             tau-b between a gene's log-normalised expression and a peak's binarised
                      accessibility across cells: (n_c - n_d) / sqrt((n0 - n1)(n0 - n2)), n0 =
                      n(n-1)/2, n1 = sum_t t(t-1)/2 over expression ties, n2 = s(s-1)/2 +
                      (n-s)(n-s-1)/2 with s open cells (Sheth, Qiu et al. 2025 Fig. S1).  Computed
                      here in closed form, n_c - n_d = 2 R_open - s(n+1) with R_open the sum of
                      mid-ranks of the open cells; NA when a margin is constant.
  ARC-E2G             Kendall joined to each ABC pair = max Kendall over overlapping peaks linked
                      to the same gene; over pairs with ABC > 0 and Kendall defined, r = slope of
                      scale(log ABC) ~ scale(Kendall) (= Pearson r), ARC = exp(log ABC + Kendall x
                      sd(log ABC) / sd(Kendall) x r); ARC = ABC elsewhere.  ABC column =
                      powerlaw.Score without Hi-C, ABC.Score with it.
  activity features   (ENCODE_rE2G) numCandidateEnhGene = rank of the element's midpoint counted
                      away from the TSS per gene (0 at the TSS); numTSSEnhGene = number of TSS500
                      windows overlapping [element start, TSS) or [TSS, element end);
                      numNearbyEnhancers / sumNearbyEnhancers = number / summed activity_base of
                      other candidate elements overlapping the element midpoint +/- 5 kb.
  feature table       union of the models' feature tables + the scE2G ARC rows (RNA_meanLogNorm,
                      RNA_pseudobulkTPM, RNA_percentCellsDetected, Kendall [max], ARC.E2G.Score,
                      ABC.Score [sum, dropped if already a feature], normalizedATAC_enh).  External
                      features are merged by closed overlap on (chr, TargetGene) and aggregated per
                      pair (mean / max / sum; NA propagates), renamed input_col -> feature, NA ->
                      fill_value (or the finite mean when fill_value is "mean").
  E2G.Score           1 / (1 + exp(-(b0 + sum_i w_i log(|x_i| + 0.01)))) with inf -> NA -> 0 first;
                      E2G.Score.cv the same with the model held out for that chromosome.
  qnorm               quantile q = (rank(score, zeros as lowest, average ties) - 1) / (N - 1), 0 for
                      zero scores; E2G.Score.qnorm = linear interpolation (extrapolated) of the
                      model's reference quantile -> score table at q, clipped at 0.
  TPM filter          when RNA_pseudobulkTPM exists and the model's tpm_threshold > 0:
                      E2G.Score.qnorm.ignoreTPM = E2G.Score.qnorm, then E2G.Score.qnorm = 0 where
                      TPM < threshold (the IGVF release names these Score.ignoreTPM / Score).
  thresholding        E2G.Score.qnorm >= score_threshold, dropping promoter-class elements unless
                      isSelfPromoter; bedpe chr, start, end, chr, TSS, TSS, gene_name, score, ., .
  QC stats            per cluster x model: elements / genes / links of thresholded non-promoter
                      predictions, mean genes per enhancer, enhancers per gene, distance, width,
                      genes with an active promoter, genes with qnorm score 0; warnings below 2e6
                      fragments, 100 cells, 1e6 UMIs; compared with Sheth, Qiu 2024 reference
                      clusters (fragments_total >= 2e6).
  benchmark           CRISPR pairs merged onto predictions by closed overlap per gene (TSS-universe
                      genes only), NAfilled; AUPRC (sklearn PR curve without its first point),
                      precision at 50 / 70 % recall, precision and recall at the model threshold,
                      each with scipy BCa bootstrap CIs; scATAC_ABC, noTPMfilter and distanceToTSS
                      rows added as upstream does.
  training            (Snakefile_training -> ENCODE_rE2G train_model) unpenalised lbfgs logistic
                      regression on log(|x| + 0.01) of NAfilled CRISPR features, full model plus
                      one model per held-out chromosome; coefficients, per-chromosome log-loss /
                      AUROC / AUPRC; qnorm-ref writes the 10,001-point reference from genome-wide
                      full-model scores (run_e2g_cv_qnormref.py).

Where this differs (all documented, none changes a score)
  * Sampling above max_cell_count uses numpy's RNG (seed 123), not R's sample().
  * Cells absent from the RNA matrix raise an error rather than an R subscript failure; a
    chromosome with no cv model falls back to the full model (upstream crashes).
  * Fragment files that are not coordinate sorted are sorted before tabix indexing.
  * The ABC stage (MACS2 peaks, candidate regions, neighbourhoods, ABC predictions) belongs to the
    ABC repository and is an input here (`run --abc-dir`); ENCODE_rE2G feature selection /
    permutation importance (run_feature_analysis: False by default) and the plotly HTML report
    are not ported -- qc writes the same statistics as TSV + PNG.

Subcommands
  setup            qnorm references (4 models) + TSS500 / gene bounds / QC reference / chrom sizes /
                   gene classes from the pinned commits into Data/scE2G/; --crispr, --gtf and
                   --example add the CRISPR benchmark, the GENCODE v43 GTF and the chr22 fixture.
  models           the embedded v3 models: features, weights, thresholds; --export writes
                   upstream-shaped model directories (feature_table.tsv, threshold files, JSON).
  frag-to-tagalign fragments -> filtered fragments, fragment_count, cell_barcodes, tagAlign(.tbi).
  bigwig           fragments -> ATAC / ATAC_norm bedGraph (+ bigWig when possible).
  kendall-pairs    narrowPeak + ABC AllPutative -> Kendall/Pairs.tsv.gz.
  atac-matrix      pairs + fragments + RNA cells -> peak x cell count matrix (MatrixMarket).
  kendall          pairs + ATAC matrix (or fragments) + RNA -> Pairs.Kendall.tsv.gz, gene
                   expression metrics, umi_count, cell_count.
  arc              ABC AllPutative + Pairs.Kendall -> EnhancerPredictionsAllPutative_ARC.tsv.gz.
  activity-features  ABC AllPutative + EnhancerList (+ gene classes) -> ActivityOnly_features.
  features         activity features + ARC -> feature_table.tsv, external_features_config.tsv,
                   ActivityOnly_plus_external_features, genomewide_features.tsv.gz.
  predict          genome-wide features -> scE2G_predictions(.threshold).tsv.gz (+ .bedpe) per model.
  qnorm-ref        genome-wide features -> a model's qnorm_reference.tsv.gz.
  gene-lists       ABC GeneList / EnhancerList + predictions -> scE2G_gene_list / element_list.
  stats            predictions -> scE2G_predictions_threshold<t>_stats.tsv.
  qc               stats files -> all_qc_stats.tsv, warnings, reference comparison, figures.
  crispr-features  predictions or features + CRISPR benchmark -> ..._features_NAfilled.tsv.gz.
  benchmark        CRISPR feature tables -> crispr_benchmarking_performance_summary.tsv.
  train            CRISPR features -> a new model directory usable by `predict --models`.
  run              the whole per-cluster pipeline from fragments + RNA + ABC outputs.
  validate-example Kendall / RNA features vs the upstream chr22 test fixture.
  validate-release recompute a released IGVF scE2G prediction set (pairs + elements + genes files).
  selftest         synthetic multiome clusters with planted enhancers; every subcommand asserted.

Validation against upstream outputs (2026-09-22)
  * validate-example: tests/expected_output/.../K562_cluster1_chr22p/Kendall/Pairs.Kendall.tsv.gz
    recomputed from the fixture's fragments and RNA matrix: 27,685 / 27,685 pairs, identical NA
    pattern (18,017 defined), Kendall max |diff| 5e-16, RNA TPM / detection / mean log-norm exact.
    (Counting Tn5 insertions instead of fragments gives max |diff| 0.41 -- Signac counts fragments.)
  * validate-release on IGVF K562 scE2G v1.2 (IGVFDS5428HHMB: IGVFFI1706PNVV pairs, IGVFFI3094BAXH
    elements, IGVFFI9905RPTO genes): numCandidateEnhGene / numTSSEnhGene / numNearbyEnhancers
    recomputed, multiome_powerlaw_v3 applied with the embedded weights and the pinned qnorm
    reference: all 11,519,202 released Score.ignoreTPM and Score values reproduced (max |diff|
    6.7e-16).

Output: Docs/scE2GPipeline/<timestamp>_<label>/.  numpy, pandas, scipy required; scikit-learn for
benchmark/train; pysam for bgzip/tabix; matplotlib optional.

Usage:
    igvfagent sce2g-pipeline setup [--crispr] [--gtf] [--example]
    igvfagent sce2g-pipeline models --export Docs/scE2GPipeline/models
    igvfagent sce2g-pipeline run --cluster K562 --fragments atac_fragments.tsv.gz --rna-matrix rna.h5ad --abc-dir ABC/K562 --models multiome_powerlaw_v3 scATAC_powerlaw_v3 --label k562
    igvfagent sce2g-pipeline kendall --pairs Pairs.tsv.gz --fragments atac_fragments.tsv.gz --rna-matrix rna_count_matrix.csv.gz --label chr22
    igvfagent sce2g-pipeline predict --features genomewide_features.tsv.gz --models multiome_powerlaw_v3 --cluster K562
    igvfagent sce2g-pipeline benchmark --crispr-features K562/multiome_powerlaw_v3/EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz
    igvfagent sce2g-pipeline validate-example
    igvfagent sce2g-pipeline validate-release --pairs IGVFFI1706PNVV.tsv.gz --elements IGVFFI3094BAXH.tsv.gz --genes IGVFFI9905RPTO.tsv.gz
    igvfagent sce2g-pipeline selftest --no-plots
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
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "scE2GPipeline"
DATA_ROOT = ROOT / "Data" / "scE2G"
MODELS_DIR = DATA_ROOT / "models"
RES_DIR = DATA_ROOT / "pipeline_resources"

UPSTREAM_REPO = "EngreitzLab/scE2G"
UPSTREAM_COMMIT = "7cb2af750fb96006f5d2b7c5475dcff30ea0e6c9"
UPSTREAM_DATE = "2026-08-13T01:37:47Z"
E2G_REPO = "EngreitzLab/ENCODE_rE2G"
E2G_COMMIT = "5930e746b863e3cab2db0e6d814170027b9a2c97"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
RAW_E2G = f"https://raw.githubusercontent.com/{E2G_REPO}/{E2G_COMMIT}/"

# name -> (url, description); fetched by `setup`
RESOURCES = {
    "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed": (
        RAW + "resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed",
        "gene TSS +/- 250 bp (gene universe, Kendall gene names, numTSSEnhGene)"),
    "CollapsedGeneBounds.hg38.intGENCODEv43.bed": (
        RAW + "resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed", "gene bodies"),
    "reference_qc_metrics_sheth_qiu_2024.tsv": (
        RAW + "resources/reference_qc_metrics_sheth_qiu_2024.tsv", "QC metrics of the Sheth, Qiu 2024 clusters"),
    "GRCh38_EBV.no_alt.chrom.sizes.tsv": (
        RAW_E2G + "reference/GRCh38_EBV.no_alt.chrom.sizes.tsv", "chromosome sizes (ENCODE_rE2G default)"),
    "gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv": (
        RAW_E2G + "resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv",
        "is_ubiquitous_uniform / P2PromoterClass per gene (ENCODE_rE2G gene_classes)"),
}
CRISPR_RESOURCE = ("EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz",
                   RAW + "resources/EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz")
GTF_RESOURCE = ("gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz",
                RAW + "resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz")
EXAMPLE_FILES = {  # chr22 test fixture used by upstream tests/test_sce2g_apply.py
    "atac_fragments.chr22p.tsv.gz": RAW + "resources/example_chr22p_multiome_cluster/atac_fragments.tsv.gz",
    "rna_count_matrix.cluster1.csv.gz": RAW + "resources/example_chr22_multiome_cluster/cluster1/rna_count_matrix.csv.gz",
    "Pairs.Kendall.expected.tsv.gz": RAW + "tests/expected_output/powerlaw_v3_models/K562_cluster1_chr22p/Kendall/Pairs.Kendall.tsv.gz",
    "stats.expected.tsv": RAW + "tests/expected_output/powerlaw_v3_models/K562_cluster1_chr22p/multiome_powerlaw_v3/scE2G_predictions_threshold0.177_stats.tsv",
}

EPSILON = 0.01
MAX_CELL_COUNT = 20000
SCORE_BASE = "E2G.Score"
FINAL_SCORE_COL = "E2G.Score.qnorm"
QC_THRESHOLDS = {"fragments_total": 2e6, "cell_count": 100, "umi_count": 1e6}
ABC_BENCHMARK_THRESHOLD = 0.015
DISTANCE_BENCHMARK_THRESHOLD = -54350
FT_COLS = ["feature", "input_col", "second_input", "aggregate_function", "fill_value", "nice_name"]
EFC_COLS = ["input_col", "source_col", "aggregate_function", "join_by", "source_file"]

# models/<model>/feature_table.tsv (feature, input_col, aggregate_function, fill_value, nice_name)
_FT_COMMON = [
    ("numTSSEnhGene", "numTSSEnhGene", "max", "0", "# TSSs between E and P"),
    ("normalizedATAC_prom", "normalized_atac_prom", "mean", "0", "ATAC signal at P"),
    ("numNearbyEnhancers", "numNearbyEnhancers", "max", "0", "# peaks within 5Kb of E"),
    ("ubiqExpressed", "is_ubiquitous_uniform", "max", "0", "Ubiquitous expression"),
]
_FT_POWERLAW = ("numCandidateEnhGene", "numCandidateEnhGene", "max", "0", "# peaks between E and P")
_FT_MEGAMAP = ("contactFrequency", "hic_contact_pl_scaled_adj", "mean", "0", "Contact frequency")
_FT_ARC = ("ARC.E2G.Score", "ARC.E2G.Score", "mean", "0", "ARC-E2G score")
MODEL_SPECS = {
    "multiome_powerlaw_v3": {"rows": _FT_COMMON + [_FT_POWERLAW, _FT_ARC], "score_threshold": ".177", "tpm_threshold": "1",
                             "pretty": "scE2G (Multiome, powerlaw)", "color": "#792374"},
    "multiome_megamap_v3": {"rows": _FT_COMMON + [_FT_MEGAMAP, _FT_ARC], "score_threshold": ".173", "tpm_threshold": "1",
                            "pretty": "scE2G (Multiome, average Hi-C)", "color": "#9b241c"},
    "scATAC_powerlaw_v3": {"rows": _FT_COMMON + [_FT_POWERLAW, ("ABC.Score", "powerlaw.Score", "sum", "0", "ABC score")],
                           "score_threshold": ".174", "tpm_threshold": "0",
                           "pretty": "scE2G (scATAC-only, powerlaw)", "color": "#006479"},
    "scATAC_megamap_v3": {"rows": _FT_COMMON[:2] + [("numNearbyEnhancers", "numNearbyEnhancers", "max", "0", "# peaks within 5 kb of E")]
                          + _FT_COMMON[3:] + [_FT_MEGAMAP, ("ABC.Score", "ABC.Score", "sum", "0", "ABC score")],
                          "score_threshold": ".187", "tpm_threshold": "0",
                          "pretty": "scE2G (scATAC-only, average Hi-C)", "color": "#00488d"},
}

# ENCODE_rE2G combine_feature_tables_apply.R: rows appended for sc-E2G biosamples
ARC_FEATURE_ROWS = [
    ("RNA_meanLogNorm", "mean_log_normalized_rna", "NA", "mean", "0", "Mean log normalized RNA expression"),
    ("RNA_pseudobulkTPM", "RnaPseudobulkTPM", "NA", "mean", "0", "RNA pseudobulk TPM"),
    ("RNA_percentCellsDetected", "RnaDetectedPercent", "NA", "mean", "0", "RNA percent cells detected"),
    ("Kendall", "Kendall", "NA", "max", "0", "Kendall correlation"),
    ("ARC.E2G.Score", "ARC.E2G.Score", "NA", "mean", "0", "ARC-E2G score"),
    ("ABC.Score", "ABC.Score", "NA", "sum", "0", "ABC score"),
    ("normalizedATAC_enh", "normalized_atac_enh", "NA", "mean", "0", "ATAC signal at E"),
]
# scE2G format_external_features_config_sc.R: how the ARC table is merged
ARC_EXTERNAL_ROWS = [("mean_log_normalized_rna", "mean"), ("RnaPseudobulkTPM", "mean"), ("RnaDetectedPercent", "mean"),
                     ("Kendall", "max"), ("ARC.E2G.Score", "mean")]
# merge_features_with_crispr_data_apply.R: score columns carried onto CRISPR pairs
SCORE_ROWS = ["E2G.Score", "E2G.Score.cv", "E2G.Score.cv.qnorm.ignoreTPM", "E2G.Score.qnorm", "E2G.Score.cv.qnorm",
              "E2G.Score.qnorm.ignoreTPM"]
ACTIVITY_CORE_COLS = ["chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS", "TargetGeneEnsembl_ID",
                      "isSelfPromoter", "isSelfGenic", "CellType", "distance"]
KENDALL_OUT_COLS = ["chr", "start", "end", "TargetGene", "PeakName", "PairName", "mean_log_normalized_rna",
                    "RnaDetectedPercent", "RnaPseudobulkTPM", "Kendall"]
GEX_COLS = ["mean_log_normalized_rna", "RnaDetectedPercent", "RnaPseudobulkTPM"]
STATS_COLS = ["n_enh_elements", "n_genes_with_enh", "n_enh_gene_links", "mean_genes_per_enh", "mean_enh_per_gene",
              "mean_dist_to_tss", "mean_enh_width", "n_genes_active_promoter", "n_genes_not_expressed", "cluster",
              "model_name", "fragments_total", "cell_count", "umi_count"]
QC_LABELS = {
    "fragments_total": "# unique ATAC fragments in cluster", "cell_count": "# cells in cluster",
    "umi_count": "# RNA UMIs in cluster", "frag_per_cell": "Mean unique ATAC fragments per cell",
    "umi_per_cell": "Mean RNA UMIs per cell", "n_enh_gene_links": "# enhancer-gene links",
    "n_enh_elements": "# unique enhancers", "mean_genes_per_enh": "Mean # genes per enhancer",
    "mean_enh_per_gene": "Mean # enhancers per gene", "n_genes_with_enh": "# genes with 1+ enhancer",
    "n_genes_active_promoter": "# genes with accessible promoter", "n_genes_not_expressed": "# genes below TPM threshold",
    "mean_dist_to_tss": "Mean distance to TSS (bp)", "mean_enh_width": "Mean width of enhancer element (bp)",
}
BENCH_COLS = ["cluster", "model", "AUPRC", "AUPRC_95CI_low", "AUPRC_95CI_high", "precision_70_pct_recall",
              "precision_70_pct_recall_95CI_low", "precision_70_pct_recall_95CI_high", "threshold_70_pct_recall",
              "precision_50_pct_recall", "precision_50_pct_recall_95CI_low", "precision_50_pct_recall_95CI_high",
              "threshold_50_pct_recall", "precision_model_threshold", "precision_model_threshold_95CI_low",
              "precision_model_threshold_95CI_high", "recall_model_threshold", "recall_model_threshold_95CI_low",
              "recall_model_threshold_95CI_high", "pct_missing_pos", "pct_missing_neg", "pct_missing_total"]
DEFAULT_TRAIN_PARAMS = {"solver": "lbfgs", "fit_intercept": True, "penalty": None, "max_iter": int(1e8),
                        "class_weight": None, "tol": 1e-4, "warm_start": False, "random_state": 0, "n_jobs": 1}

# Logistic-regression weights extracted from models/<model>/model.pkl (full model) and
# models/<model>/cv_models/model_test_<chr>.pkl (held-out-chromosome models) at the pinned commit:
# [intercept, w_1 .. w_6] in the order of the model's feature list (sklearn LogisticRegression,
# penalty=None, lbfgs; score = 1 / (1 + exp(-(b0 + sum w_i * log(|x_i| + 0.01)))).
MODEL_WEIGHTS = {
    'multiome_megamap_v3': {
        'full': [2.4433587422986625, -1.1744124875247717, -1.008509065729552, 0.20838932583533018, -0.2842633949780195, -0.29086647341883515, 1.7372358636397305],
        'cv': {
            'chr1': [2.501908990214827, -1.149120809126327, -0.8910958564847061, 0.2319251022305382, -0.2849452652604566, -0.28327433186172957, 1.7383367385556137],
            'chr2': [2.5111223412037136, -1.1802715110594733, -1.0168525326395375, 0.2417781112299334, -0.27658564135892383, -0.27082809722051765, 1.7405401932720932],
            'chr3': [2.505318494476338, -1.1926692218565602, -1.1110881580517595, 0.2021102711204486, -0.26372228761428135, -0.32017764785905123, 1.7632782486037084],
            'chr4': [2.612516949791687, -1.1201009931191235, -0.943835345044745, 0.20564112507076446, -0.2910583550680719, -0.2257554162933038, 1.741034331159168],
            'chr5': [2.487481942738878, -1.1951992078685494, -1.0652270452179873, 0.2014027034067702, -0.27687527277177165, -0.29689196620662606, 1.7440266357714316],
            'chr6': [2.618884702424139, -1.1298190430610002, -0.9214363878651195, 0.20004391110348546, -0.27966845887773356, -0.2530832130058837, 1.7463668770505765],
            'chr7': [2.572064961779358, -1.1467707095435997, -1.0615393614693538, 0.20838827968294502, -0.2836354343524355, -0.22975869917912548, 1.7291985750576173],
            'chr8': [2.310167011390833, -1.267452505199429, -0.9913498696254966, 0.19273635763520705, -0.2692987668586971, -0.2958790324791736, 1.6623362767500443],
            'chr9': [2.5002848839451155, -1.1547249271264277, -0.9469360142676797, 0.19262894070674175, -0.28554008175471296, -0.27697597969193005, 1.7399006877055923],
            'chr10': [2.429324267283811, -1.1757750023096296, -1.0239878466935595, 0.20615060543005392, -0.2801310956383961, -0.2961002220758504, 1.734372033285947],
            'chr11': [2.3601154935572435, -1.2073756837008063, -0.6057641264717973, 0.18665936140589182, -0.29554777225461276, -0.33967571004653996, 1.7389686737908563],
            'chr12': [2.591650652721799, -1.147083060173798, -1.2186638915189865, 0.23130109459618775, -0.2780423675863387, -0.24780435585972588, 1.7354582375451],
            'chr13': [2.432824646885325, -1.1716770643330794, -1.0115103232026144, 0.207474289209578, -0.2835555782724288, -0.29105400420425287, 1.7347348879760152],
            'chr14': [2.468901120288694, -1.17455754150228, -0.9834427051711097, 0.21242232314438253, -0.28544172382584576, -0.2892562724347326, 1.7368077405285387],
            'chr15': [2.3463902725329935, -1.1701993683914111, -0.994285642208913, 0.21648409996070628, -0.29106072315875753, -0.3143459730173496, 1.7408891315089141],
            'chr16': [2.39282796714491, -1.1652630393879944, -1.0255524167642016, 0.20761186694959055, -0.2823352359928952, -0.29568368894755354, 1.728031411819456],
            'chr17': [2.3363620241946514, -1.1752497686969487, -1.0601927449861603, 0.19598350707672113, -0.2916583525688803, -0.3209119339483683, 1.7565969083006485],
            'chr18': [2.4209096545386277, -1.165181828813371, -1.0438267998844237, 0.21388495753990466, -0.28459719227710883, -0.2887290783162893, 1.7367489414931654],
            'chr19': [2.446947980742945, -1.1635041038902927, -1.1426590151117992, 0.23930612711502575, -0.31077455517365393, -0.27625698574662466, 1.7468984522363318],
            'chr20': [2.313143995497986, -1.1920551767203942, -1.0280514591866259, 0.20749094919554706, -0.2865788076836745, -0.3404559837842214, 1.7493605366938345],
            'chr21': [2.4380299607857907, -1.170793831290677, -0.9806456734664206, 0.20884253701755542, -0.2841852796031461, -0.2938214477300896, 1.7362902233090314],
            'chr22': [2.476496861526364, -1.2004323862206425, -0.9336557764179433, 0.20950857177681728, -0.2774431215971875, -0.29258608783934764, 1.7327847294680807],
            'chrX': [2.162732748775956, -1.1903550019604194, -1.177163717001577, 0.1791736242394091, -0.2978938306984093, -0.34314706375184506, 1.7494630913056146],
        },
    },
    'multiome_powerlaw_v3': {
        'full': [2.219492420277004, -1.3979194615780126, -0.9776556852269032, 0.17520928853400672, -0.30598394888675984, 0.49466768260841937, 1.7264545609649855],
        'cv': {
            'chr1': [2.1939627273953635, -1.3883447796122232, -0.8648557234900656, 0.20235164662951696, -0.30631057080475, 0.4931871957780802, 1.70629991486508],
            'chr2': [2.2857470570530665, -1.3921999225381785, -0.9735864370024317, 0.20609346437242182, -0.29997284640311034, 0.48524856751549583, 1.7459740625477538],
            'chr3': [2.3108181656556646, -1.4042456129007084, -1.0843843419705599, 0.16691664566244782, -0.2877970984854153, 0.4990188088776576, 1.7411202330762945],
            'chr4': [2.2757519196784606, -1.347637440117777, -0.9202424025621325, 0.17541647622010628, -0.3140511976997994, 0.47238283367169037, 1.7468763190984178],
            'chr5': [2.2686688472246477, -1.4245191753999353, -1.025670116173678, 0.16723905595880356, -0.29964225456128174, 0.5030847868532214, 1.7340657528933423],
            'chr6': [2.42338589526147, -1.3454095160736255, -0.8651066134426846, 0.16727175806457603, -0.2970049876459892, 0.4913576375852409, 1.7723059161034147],
            'chr7': [2.2381748035146942, -1.3715637104808478, -1.010227379794922, 0.17873470809092812, -0.30689244265508303, 0.472361788735751, 1.7289620519640094],
            'chr8': [2.178665028686572, -1.457010826413716, -0.951373323589456, 0.16175030463283324, -0.2907734684904761, 0.4687159849036991, 1.6611322389893373],
            'chr9': [2.240529984389562, -1.390286264742606, -0.9354424170608795, 0.16283057968109368, -0.3065196892165367, 0.5004875437776588, 1.7349192250547616],
            'chr10': [2.1902173590025096, -1.4024841666580707, -0.9969834775176254, 0.17249481950302015, -0.3012279681571069, 0.4941146109072865, 1.7127691834628818],
            'chr11': [2.2245111668261863, -1.3994319170330904, -0.5522362965174883, 0.1483649253813485, -0.3124603928707867, 0.4772358028386606, 1.694513452411317],
            'chr12': [2.322698567607978, -1.3873280168569815, -1.1903262528010277, 0.20481824228387274, -0.30372395451187323, 0.5253867600466093, 1.7791554611810734],
            'chr13': [2.2033112416583442, -1.399572373779503, -0.9811619317260307, 0.17436606039748928, -0.305221515291893, 0.4961949914392068, 1.7222499176944575],
            'chr14': [2.2376496390460434, -1.3793642447375891, -0.9513748412078561, 0.17900265076618477, -0.30771083813365224, 0.4633695532871319, 1.7066354331433113],
            'chr15': [2.1480735235188524, -1.4044507037687146, -0.9747935534546355, 0.18083649631919563, -0.31149621607988753, 0.5155903001336052, 1.727725244422021],
            'chr16': [2.175608700560755, -1.3946933952332263, -0.9970946111786818, 0.17432579402369175, -0.30406396768401817, 0.5041021056690197, 1.7203482153000709],
            'chr17': [2.129307935203679, -1.4260644249638854, -1.0609444409479891, 0.16204592404267268, -0.3129707387705637, 0.5347684891460376, 1.7462682284396454],
            'chr18': [2.1986156074079712, -1.3811509954098555, -1.0116860883491952, 0.18000019411547316, -0.3056511530194176, 0.48317095979262414, 1.720738247010547],
            'chr19': [2.2088642154092737, -1.3793949940728658, -1.0978802422052514, 0.20449900494883982, -0.333841885471096, 0.4575104748663625, 1.7187617336430718],
            'chr20': [2.1723862432114385, -1.4139851151592728, -1.0107916001298083, 0.17167541460726332, -0.3073011725969349, 0.5162761483115608, 1.731160624596585],
            'chr21': [2.2241809989346146, -1.391362575093417, -0.9449510924418998, 0.1755736104891938, -0.3063277233353509, 0.4910719150790711, 1.723141118745967],
            'chr22': [2.2436269636878463, -1.4290712643370207, -0.9065602240594677, 0.17573800460775207, -0.2994740078235776, 0.5029213661131795, 1.723375287996797],
            'chrX': [1.9912154121701195, -1.4316259172612404, -1.1685431616284958, 0.14335488956514253, -0.32032677061667525, 0.5197033182036824, 1.7171867532404865],
        },
    },
    'scATAC_megamap_v3': {
        'full': [3.2984457799101166, -1.378305606448704, -1.281746388622007, 0.19891871258348448, -0.3216364929434469, -0.3641551075479487, 1.93492195907693],
        'cv': {
            'chr1': [3.339964183881814, -1.3467636330336974, -1.1528756682700982, 0.22715894391250024, -0.32079311000541416, -0.3717817617402235, 1.9456481710837323],
            'chr2': [3.2995491910695294, -1.3913711456060365, -1.292201677770545, 0.22844900582210131, -0.31313206761713014, -0.3446933804787357, 1.9183314698493492],
            'chr3': [3.3821750424037753, -1.4021939290871182, -1.3833051973497013, 0.19748385686941075, -0.3005127523228014, -0.4236141500215906, 1.991940105899637],
            'chr4': [3.484202954591855, -1.3177880760111331, -1.212119067991408, 0.19567657841988154, -0.3284763844778752, -0.31349144133394674, 1.9577670194173864],
            'chr5': [3.342422316503871, -1.4004232089480126, -1.3246899539204369, 0.1924665760193672, -0.31252977272766375, -0.36882387878488215, 1.9363688873849272],
            'chr6': [3.4158433580685843, -1.345337885805409, -1.2087081516774734, 0.18654362101068125, -0.31673492976014095, -0.3247520021732308, 1.9244955246944913],
            'chr7': [3.5877451754827563, -1.320957123788604, -1.2808006036539723, 0.20223537264751967, -0.3227174912649355, -0.3367590544035382, 2.0072626989588174],
            'chr8': [3.1679276682225046, -1.4586917561566286, -1.2415738808717838, 0.19209865223108571, -0.3056247605850662, -0.3932327416986415, 1.8889971758403477],
            'chr9': [3.389044982825055, -1.3472421065081235, -1.2209876301140068, 0.1866055068757197, -0.3233235811329021, -0.3693189992540037, 1.9687317215377824],
            'chr10': [3.2970053064326055, -1.3798619300500483, -1.2887481323039354, 0.19843360182578065, -0.31831430071850314, -0.3790189901934776, 1.9445016682560738],
            'chr11': [3.152127744264557, -1.4094330557623334, -1.018829229731441, 0.177984148963964, -0.33621921394095994, -0.3653595874935037, 1.888665439494863],
            'chr12': [3.3931229288799503, -1.3493580800654894, -1.4154752059676783, 0.21807463195200974, -0.3151951805916807, -0.304734076811195, 1.8984800614482493],
            'chr13': [3.289920286104253, -1.3748699539281999, -1.2847133418313657, 0.19764552115731138, -0.32110065763887696, -0.3678084732933296, 1.936564094499893],
            'chr14': [3.3214624173997644, -1.3806807224635238, -1.2733559740395626, 0.20288458113630986, -0.3226966211233365, -0.3582190836278706, 1.9298978259605475],
            'chr15': [3.2318457737877693, -1.3678135162510938, -1.2563852814668555, 0.20475325261119634, -0.3271271885045493, -0.38866974194953474, 1.946797701244067],
            'chr16': [3.230882860089781, -1.3685939983319768, -1.3056451076597144, 0.19829372695216405, -0.3190664485727391, -0.37575377360741175, 1.9282773255803725],
            'chr17': [3.260802324087337, -1.3878577375287566, -1.3518439697563438, 0.18746877591457065, -0.329330813729222, -0.40822424154566295, 1.9828546762903805],
            'chr18': [3.2671534387346846, -1.37137517512261, -1.3151982833843594, 0.20252471915368642, -0.3217876332585654, -0.3634916258298553, 1.9317872549552573],
            'chr19': [3.2199283139982984, -1.405418202943434, -1.438819267482774, 0.21974987251913888, -0.34289928366965944, -0.24621920773595626, 1.8132686219231975],
            'chr20': [3.1751123692953596, -1.3979456883434977, -1.3054648281144479, 0.1957920882626967, -0.3238228477824203, -0.40734757132338806, 1.9418081615272786],
            'chr21': [3.30441550157368, -1.3714531679202406, -1.2368224943968456, 0.198565344880682, -0.3217773967883171, -0.36007398831951964, 1.9299270826281885],
            'chr22': [3.309298684737774, -1.4065762699038329, -1.218063083195366, 0.19942961609972942, -0.3155717076756222, -0.36481797956303263, 1.9244860322068715],
            'chrX': [3.035940206562772, -1.3980921592996727, -1.4679496655319666, 0.1695156803620597, -0.34062028881468426, -0.3943507165186798, 1.9350254175047295],
        },
    },
    'scATAC_powerlaw_v3': {
        'full': [2.993120557537447, -1.5717625402474111, -1.2918132235914868, 0.16387829511164112, -0.3355423988179916, 0.4829454886715295, 1.8416408382036955],
        'cv': {
            'chr1': [2.941347263776975, -1.5595314962348443, -1.1741774836901322, 0.193099308110457, -0.33276063079074103, 0.4887825025269404, 1.8193879072048682],
            'chr2': [3.0506610956769618, -1.571575110315957, -1.2890763156423188, 0.1929340281657668, -0.3288886056192227, 0.4784402103044994, 1.859587792772926],
            'chr3': [3.1147770580874306, -1.5817158308007966, -1.3973625653488742, 0.15869040439526097, -0.3171578334424336, 0.4976093663340716, 1.8704424804359798],
            'chr4': [3.087212760878475, -1.5170689349454483, -1.2308881909083769, 0.16301638906901378, -0.3436775141779383, 0.4678556277894472, 1.877592201807633],
            'chr5': [3.0534353169187627, -1.596182359648272, -1.3298570804512133, 0.1560744757168292, -0.3275882218498133, 0.48884585011168713, 1.8478243864938728],
            'chr6': [3.151351214633753, -1.531576100568508, -1.1953007831589515, 0.1521083504694365, -0.3292855779672764, 0.47741795436374024, 1.8726750540970458],
            'chr7': [3.220907266091516, -1.5265685955069286, -1.2804198783897909, 0.1686233532635038, -0.3382157034109679, 0.4971610266527096, 1.9307391643535743],
            'chr8': [2.9289492152965524, -1.6224471402166767, -1.2438352943934634, 0.15524815144275356, -0.3188636232377743, 0.46379208136664113, 1.7769848447011534],
            'chr9': [3.0610187666290134, -1.5597160589767003, -1.2526253807096808, 0.15253469067765513, -0.33672907386516787, 0.503090400396573, 1.875783044469919],
            'chr10': [2.9636393456049035, -1.5775042297088044, -1.301756006961139, 0.16193865011553443, -0.33130308002387976, 0.48624738923423205, 1.8309480030751193],
            'chr11': [2.8724500698156414, -1.5589558164467423, -0.9721349177559142, 0.14177461555814735, -0.3461993276283949, 0.43016324512935056, 1.7665890294446158],
            'chr12': [3.109309703587873, -1.563290045321738, -1.4284402711923128, 0.19129663702272293, -0.3331661335468229, 0.5114314933507516, 1.88751866144296],
            'chr13': [2.9777907934932055, -1.572982521258498, -1.2959899419242311, 0.16260078858890636, -0.3347700739300987, 0.48546601722994787, 1.8387524072059247],
            'chr14': [2.972951049113419, -1.5548086044996836, -1.2794064329672203, 0.16789376403440595, -0.33663499445400613, 0.4461796838950095, 1.8082142040076787],
            'chr15': [2.946983429496799, -1.57198215000949, -1.2773128104134799, 0.1674463387709713, -0.33961334948571165, 0.5036173198426784, 1.8485527673325692],
            'chr16': [2.950960180907793, -1.573674675044157, -1.323733615638227, 0.16290601260866755, -0.33311176580981816, 0.5048022851372763, 1.844807612853663],
            'chr17': [2.9840911147248224, -1.6081623114584658, -1.3952172358363368, 0.15098925993525458, -0.3426742445242785, 0.5305385286167147, 1.8859829294410702],
            'chr18': [2.961093115952849, -1.5564478251941276, -1.322415855687391, 0.16699207449275003, -0.33464048340356245, 0.47038927574794903, 1.8305552642515586],
            'chr19': [2.814987028686231, -1.5828171318660096, -1.4254088452534266, 0.19110629501217666, -0.35861257845902506, 0.39716493648311063, 1.7433024246496547],
            'chr20': [2.9324566119078352, -1.5837413530555549, -1.3261352702177025, 0.15825910125736764, -0.33674142874036317, 0.49214453786283974, 1.8339837135843478],
            'chr21': [2.9969414126490626, -1.5603039825624088, -1.2430773742713732, 0.16398203266508116, -0.3359388729943667, 0.4737175916102788, 1.8340081854142416],
            'chr22': [2.9985158698887155, -1.6031697344059908, -1.2315203587652443, 0.16405001113123102, -0.3297061465839012, 0.4873546415082759, 1.8320427875849814],
            'chrX': [2.740352136254066, -1.6100754948293474, -1.5003363222952415, 0.13350090217389837, -0.3539326750461904, 0.5032339441930415, 1.8267698305133333],
        },
    },
}

INK, INK2, AXIS, SURFACE, GREY = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb", "#96a0b3"
FALLBACK_COLORS = ["#429130", "#c5373d", "#e96a00", "#ca9b23", "#2a78d6", "#4a3aa7"]

log = logging.getLogger("sce2g_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"sce2g_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label))[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir_for(args: argparse.Namespace, default_label: str) -> Path:
    """--out-dir if given (created), else a fresh timestamped run directory."""
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
        raise SystemExit("this subcommand needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


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


def _style(ax, title: str = "", xlabel: str = "", ylabel: str = "") -> None:
    if title:
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=8)
    ax.set_ylabel(ylabel, color=INK2, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=7)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    try:
        fig.tight_layout()
    except Exception:  # noqa: BLE001
        pass
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "Sequence[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
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
            return "nan"
        if isinstance(x, (int,)) or (isinstance(x, float) and float(x).is_integer() and abs(x) > 1000):
            return f"{int(x):,}"
        return f"{float(x):.{nd}g}" if abs(float(x)) >= 1000 else f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "t", "1", "yes", "y"}


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_tsv(df, path: Path, quiet: bool = False, **kw) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None, **kw)
    if not quiet:
        print(f"TSV: {path}")
    return path


def write_text(path: Path, text: str, quiet: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if not quiet:
        print(f"Wrote: {path}")
    return path


def read_tsv(path, **kw):
    pd = _pd()
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


def write_summary(out: Path, summary: dict, report_lines: "List[str]") -> None:
    (out / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2))
    print(f"JSON: {out / 'summary.json'}")
    (out / "report.md").write_text("\n".join(report_lines) + "\n")
    print(f"Report: {out / 'report.md'}")


def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_chrom_sizes(path) -> "List[Tuple[str, int]]":
    out = []
    with _open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and not f[0].startswith("#"):
                try:
                    out.append((f[0], int(f[1])))
                except ValueError:
                    continue
    return out


def read_bed_like(path, names: "Sequence[str]", ncols: "Optional[int]" = None):
    """Headerless BED-ish file; extra columns are dropped, missing ones filled with NA."""
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    df = df.iloc[:, :len(names)]
    df.columns = list(names)[:df.shape[1]]
    for c in names:
        if c not in df.columns:
            df[c] = pd.NA
    return df


def read_tss500(path):
    return read_bed_like(path, ["chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"])


def model_threshold_str(t: str) -> str:
    """'.177' (threshold file suffix) -> '0.177' (as Python's float() prints it; the filename token)."""
    return str(float(t))


def download(url: str, dest: Path, force: bool = False) -> str:
    dest = Path(dest)
    if dest.is_file() and dest.stat().st_size > 0 and not force:
        return "present"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, 1 << 20)
    tmp.replace(dest)
    print(f"Wrote: {dest}")
    return "fetched"


# ---------------------------------------------------------------------------
# Interval joins
# ---------------------------------------------------------------------------

_KEY_SPAN = 1 << 34  # > any chromosome coordinate; keys are laid end to end on one axis


def keyed_overlap(a, b, a_key: "Sequence[str]", b_key: "Sequence[str]", closed: bool = True,
                  a_cols=("start", "end"), b_cols=("start", "end")):
    """Row-index pairs (ia, ib) of `a` and `b` with equal key columns and overlapping intervals.

    closed=True reproduces GenomicRanges::findOverlaps on IRanges(start, end) built from the numbers
    as written (a.start <= b.end and b.start <= a.end, so book-ended intervals overlap);
    closed=False is bedtools' half-open BED overlap.  Each key (e.g. chr + TargetGene) is shifted
    onto its own stretch of a single axis so the binned join in eqtl_enrichment_skill.overlap_join
    runs once for all keys.
    """
    pd, np = _pd(), _np()
    from eqtl_enrichment_skill import overlap_join  # noqa: WPS433
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    def key(df, cols):
        k = df[cols[0]].astype(str)
        for c in cols[1:]:
            k = k + "\x1f" + df[c].astype(str)
        return k.reset_index(drop=True)

    ka, kb = key(a, list(a_key)), key(b, list(b_key))
    codes, uniq = pd.factorize(pd.concat([ka, kb], ignore_index=True))
    ca, cb = codes[:len(a)].astype(np.int64), codes[len(a):].astype(np.int64)
    add = 1 if closed else 0
    A = pd.DataFrame({"chr": "k", "start": ca * _KEY_SPAN + a[a_cols[0]].to_numpy(np.int64),
                      "end": ca * _KEY_SPAN + a[a_cols[1]].to_numpy(np.int64) + add, "ia": np.arange(len(a))})
    B = pd.DataFrame({"chr": "k", "start": cb * _KEY_SPAN + b[b_cols[0]].to_numpy(np.int64),
                      "end": cb * _KEY_SPAN + b[b_cols[1]].to_numpy(np.int64) + add, "ib": np.arange(len(b))})
    m = overlap_join(A, B)
    if len(m) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    m = m.sort_values(["ia", "ib"], kind="mergesort")
    return m["ia"].to_numpy(np.int64), m["ib"].to_numpy(np.int64)


def group_aggregate(values, groups, n_groups: int, how: str):
    """Aggregate `values` by integer `groups` with R semantics: NA propagates (mean/max/min/sum without na.rm)."""
    pd, np = _pd(), _np()
    s = pd.Series(np.asarray(values, dtype=float))
    g = pd.Series(np.asarray(groups, dtype=np.int64))
    fn = {"mean": "mean", "max": "max", "min": "min", "sum": "sum"}.get(how)
    if fn is None:
        raise SystemExit(f"unsupported aggregate_function {how!r} (mean/max/min/sum)")
    agg = s.groupby(g).agg(fn)
    has_na = s.isna().groupby(g).any()
    agg[has_na.reindex(agg.index).to_numpy()] = np.nan
    out = np.full(n_groups, np.nan)
    out[agg.index.to_numpy()] = agg.to_numpy()
    return out


# ---------------------------------------------------------------------------
# Models: embedded v3 models and user model directories
# ---------------------------------------------------------------------------

def embedded_feature_table(name: str):
    pd = _pd()
    rows = [(f, i, "NA", a, fv, nn) for (f, i, a, fv, nn) in MODEL_SPECS[name]["rows"]]
    return pd.DataFrame(rows, columns=FT_COLS)


def _read_feature_table(path):
    pd = _pd()
    ft = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in FT_COLS:
        if c not in ft.columns:
            ft[c] = "NA" if c == "second_input" else ""
    return ft[FT_COLS]


def _threshold_from_dir(d: Path, prefix: str) -> "Optional[str]":
    hits = sorted(p.name for p in d.glob(prefix + "*"))
    if len(hits) > 1:
        raise SystemExit(f"more than one {prefix}* file in {d}")
    return hits[0][len(prefix):] if hits else None


def _weights_from_pickle(path: Path, features: "Sequence[str]") -> "List[float]":
    import pickle
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with open(path, "rb") as fh:
            m = pickle.load(fh)
    names = list(getattr(m, "feature_names_in_", features))
    if list(names) != list(features):
        raise SystemExit(f"{path}: model features {names} != feature table {list(features)}")
    return [float(m.intercept_[0])] + [float(x) for x in m.coef_[0]]


def resolve_model(spec: str, models_dir: "Optional[Path]" = None) -> dict:
    """An embedded model name (multiome_powerlaw_v3, ...) or a model directory -> model dict.

    Keys: name, features, feature_table (DataFrame), full ([b0, w...]), cv ({chr: [b0, w...]}),
    score_threshold ('.177' style string or None), tpm_threshold (float), qnorm (Path or None),
    polynomial (bool), source.
    """
    models_dir = Path(models_dir) if models_dir else MODELS_DIR
    p = Path(spec)
    if spec in MODEL_SPECS and not p.is_dir():
        ft = embedded_feature_table(spec)
        w = MODEL_WEIGHTS[spec]
        q = models_dir / spec / "qnorm_reference.tsv.gz"
        return {"name": spec, "features": list(ft["feature"]), "feature_table": ft, "full": list(w["full"]),
                "cv": {k: list(v) for k, v in w["cv"].items()}, "score_threshold": MODEL_SPECS[spec]["score_threshold"],
                "tpm_threshold": float(MODEL_SPECS[spec]["tpm_threshold"]), "qnorm": q if q.is_file() else None,
                "qnorm_expected": q, "polynomial": False, "source": "embedded (models/%s at %s)" % (spec, UPSTREAM_COMMIT[:7])}
    if not p.is_dir():
        raise SystemExit(f"unknown model {spec!r}: not one of {sorted(MODEL_SPECS)} and not a directory")
    ft = _read_feature_table(p / "feature_table.tsv")
    features = list(ft["feature"])
    coef_json = p / "model_coefficients.json"
    polynomial = False
    if coef_json.is_file():
        blob = json.loads(coef_json.read_text())
        features = blob.get("model_features", features)
        full, cv = blob["full"], blob.get("cv", {})
        polynomial = bool(blob.get("polynomial", False))
    elif (p / "model.pkl").is_file():
        full = _weights_from_pickle(p / "model.pkl", features)
        cv = {}
        for f in sorted((p / "cv_models").glob("model_test_*.pkl")):
            cv[f.stem[len("model_test_"):]] = _weights_from_pickle(f, features)
    else:
        raise SystemExit(f"{p}: needs model_coefficients.json or model.pkl")
    q = p / "qnorm_reference.tsv.gz"
    if not q.is_file():
        q2 = models_dir / p.name / "qnorm_reference.tsv.gz"
        q = q2 if q2.is_file() else None
    tpm = _threshold_from_dir(p, "tpm_threshold_")
    return {"name": p.name, "features": list(ft["feature"]), "model_features": features, "feature_table": ft,
            "full": [float(x) for x in full], "cv": {k: [float(x) for x in v] for k, v in cv.items()},
            "score_threshold": _threshold_from_dir(p, "score_threshold_"), "tpm_threshold": float(tpm) if tpm else 0.0,
            "qnorm": q, "qnorm_expected": p / "qnorm_reference.tsv.gz", "polynomial": polynomial, "source": str(p)}


def transform_features(df, features: "Sequence[str]", epsilon: float = EPSILON, polynomial: bool = False):
    """log(|x| + epsilon) of the model features (booleans count as 0/1); degree-2 terms if polynomial."""
    np = _np()
    cols = []
    for f in features:
        c = df[f]
        if c.dtype == object:
            low = c.astype(str).str.strip().str.lower().to_numpy()
            num = _pd().to_numeric(c, errors="coerce").to_numpy(float)
            cols.append(np.where(low == "true", 1.0, np.where(low == "false", 0.0, num)))
        else:
            cols.append(pd_numeric(c))
    X = np.column_stack(cols) if cols else np.zeros((len(df), 0))
    M = np.log(np.abs(X) + epsilon)
    if polynomial:
        cols = [M[:, i] for i in range(M.shape[1])]
        extra = []
        for i in range(M.shape[1]):
            for j in range(i, M.shape[1]):
                extra.append(M[:, i] * M[:, j])
        M = np.column_stack([np.ones(len(M))] + cols + extra)
    return M


def logistic_score(M, weights: "Sequence[float]"):
    np = _np()
    w = np.asarray(weights, dtype=float)
    if M.shape[1] + 1 != len(w):
        raise SystemExit(f"model has {len(w) - 1} weights but {M.shape[1]} features")
    from scipy.special import expit  # type: ignore
    return expit(w[0] + M @ w[1:])


# ---------------------------------------------------------------------------
# Fragments -> filtered fragments, counts, barcodes, tagAlign, coverage
# ---------------------------------------------------------------------------

def read_fragments(path, usecols=(0, 1, 2, 3), chunksize: "Optional[int]" = None):
    pd = _pd()
    names = ["chr", "start", "end", "barcode", "count"][:len(usecols)]
    return pd.read_csv(path, sep="\t", header=None, comment="#", usecols=list(usecols), names=names,
                       dtype={"chr": str, "barcode": str}, chunksize=chunksize, low_memory=False)


def _is_coord_sorted(chrom, start) -> bool:
    np = _np()
    if len(chrom) == 0:
        return True
    ch = np.asarray(chrom)
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    blocks = ch[np.concatenate([[0], change])]
    if len(set(blocks)) != len(blocks):
        return False
    st = np.asarray(start)
    d = np.diff(st)
    d[change - 1] = 0
    return bool((d >= 0).all())


def bgzip_tabix(plain: Path, preset: str = "bed") -> Path:
    """Compress with BGZF and write a tabix index (pysam); plain gzip without an index if pysam is missing."""
    plain = Path(plain)
    try:
        import pysam  # type: ignore
        out = pysam.tabix_index(str(plain), preset=preset, force=True, keep_original=False)
        return Path(out)
    except ImportError:
        gz = plain.with_name(plain.name + ".gz")
        with open(plain, "rb") as fi, gzip.open(gz, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        plain.unlink()
        print(f"pysam not installed: {gz} is plain gzip without a tabix index (pip install pysam)")
        return gz


def fragments_to_tags(frag):
    """tagAlign rows (chr, start, end, N, 1000, strand) for a fragments frame, unsorted."""
    pd, np = _pd(), _np()
    s = frag["start"].to_numpy(np.int64)
    e = frag["end"].to_numpy(np.int64)
    mid = (s + e) // 2  # awk int((s+e)/2); coordinates are non-negative
    chrom = frag["chr"].to_numpy()
    tags = pd.DataFrame({"chr": np.concatenate([chrom, chrom]), "start": np.concatenate([s, mid + 1]),
                         "end": np.concatenate([mid, e]), "name": "N", "score": 1000,
                         "strand": np.concatenate([np.full(len(s), "+"), np.full(len(s), "-")])})
    return tags


def sort_by_chrom_order(df, order: "Sequence[str]", cols=("start", "end"), extra: "Sequence[str]" = ()):
    np = _np()
    rank = {c: i for i, c in enumerate(order)}
    r = df["chr"].map(rank).to_numpy()
    keys = [df[c].to_numpy() for c in reversed(list(cols) + list(extra))] + [r]
    idx = np.lexsort(keys)
    return df.iloc[idx].reset_index(drop=True)


def frag_to_tagalign(frag_path, chrom_sizes, cluster_dir: Path, preprocessed: bool = False) -> dict:
    """process_fragment_file + get_cell_barcodes + frag_to_tagAlign for one cluster."""
    pd = _pd()
    cluster_dir = Path(cluster_dir)
    (cluster_dir / "tagAlign").mkdir(parents=True, exist_ok=True)
    (cluster_dir / "Kendall").mkdir(parents=True, exist_ok=True)
    sizes = read_chrom_sizes(chrom_sizes)
    order = [c for c, _ in sizes]
    keep = set(order)
    raw = pd.read_csv(frag_path, sep="\t", header=None, comment="#", dtype={0: str, 3: str}, low_memory=False)
    n_raw = len(raw)
    frag = raw if preprocessed else raw[raw[0].isin(keep)]
    dropped = sorted(set(raw[0].unique()) - keep)
    out = {"fragments_raw": n_raw, "chromosomes_dropped": dropped}
    if not preprocessed:
        if not _is_coord_sorted(frag[0].to_numpy(), frag[1].to_numpy()):
            log.info("fragments not coordinate-sorted; sorting before tabix")
            frag = frag.rename(columns={0: "chr", 1: "start", 2: "end"})
            frag = sort_by_chrom_order(frag, order).rename(columns={"chr": 0, "start": 1, "end": 2})
        plain = cluster_dir / "fragments_filtered.tsv"
        frag.to_csv(plain, sep="\t", header=False, index=False)
        out["fragments_filtered"] = str(bgzip_tabix(plain, "bed"))
        print(f"Wrote: {out['fragments_filtered']}")
    n = len(frag)
    write_text(cluster_dir / "fragment_count.txt", f"{n}\n", quiet=True)
    barcodes = list(dict.fromkeys(frag[3].astype(str)))
    write_text(cluster_dir / "Kendall" / "cell_barcodes.txt", "\n".join(barcodes) + ("\n" if barcodes else ""), quiet=True)
    f4 = frag.iloc[:, :4].copy()
    f4.columns = ["chr", "start", "end", "barcode"]
    tags = fragments_to_tags(f4[f4["chr"].isin(keep)])
    tags = sort_by_chrom_order(tags, order, extra=("strand",))
    plain = cluster_dir / "tagAlign" / "tagAlign.sort"
    tags.to_csv(plain, sep="\t", header=False, index=False)
    out["tagAlign"] = str(bgzip_tabix(plain, "bed"))
    print(f"Wrote: {out['tagAlign']}")
    out.update({"fragment_count": n, "n_barcodes": len(barcodes), "n_tags": len(tags),
                "fragment_count_file": str(cluster_dir / "fragment_count.txt"),
                "cell_barcodes_file": str(cluster_dir / "Kendall" / "cell_barcodes.txt")})
    return out


def fragment_coverage(frag, sizes: "Dict[str, int]", scale: float = 1.0):
    """bedtools genomecov -bg of fragment intervals (half-open), runs of equal non-zero coverage."""
    pd, np = _pd(), _np()
    rows = []
    for chrom, sub in frag.groupby("chr", sort=False):
        s = sub["start"].to_numpy(np.int64)
        e = sub["end"].to_numpy(np.int64)
        if chrom in sizes:
            e = np.minimum(e, sizes[chrom])
        pos = np.concatenate([s, e])
        delta = np.concatenate([np.ones(len(s)), -np.ones(len(e))])
        up, inv = np.unique(pos, return_inverse=True)
        d = np.zeros(len(up))
        np.add.at(d, inv, delta)
        cov = np.cumsum(d)
        st, en, cv = up[:-1], up[1:], cov[:-1]
        keep = cv > 0
        st, en, cv = st[keep], en[keep], cv[keep]
        if len(st) == 0:
            continue
        # merge contiguous runs with identical coverage
        brk = np.ones(len(st), dtype=bool)
        brk[1:] = (st[1:] != en[:-1]) | (cv[1:] != cv[:-1])
        gid = np.cumsum(brk) - 1
        m_st = st[brk]
        m_en = np.zeros(gid[-1] + 1, dtype=np.int64)
        np.maximum.at(m_en, gid, en)
        rows.append(pd.DataFrame({"chr": chrom, "start": m_st, "end": m_en, "value": cv[brk] * scale}))
    if not rows:
        return pd.DataFrame(columns=["chr", "start", "end", "value"])
    return pd.concat(rows, ignore_index=True)


def write_bigwig(bg, sizes: "List[Tuple[str, int]]", bw_path: Path) -> "Optional[Path]":
    try:
        import pyBigWig  # type: ignore
    except ImportError:
        exe = shutil.which("bedGraphToBigWig")
        return None if exe is None else Path("__binary__")
    bw = pyBigWig.open(str(bw_path), "w")
    present = [(c, s) for c, s in sizes if c in set(bg["chr"])]
    bw.addHeader(present)
    for c, _ in present:
        sub = bg[bg["chr"] == c]
        bw.addEntries([c] * len(sub), sub["start"].astype(int).tolist(), ends=sub["end"].astype(int).tolist(),
                      values=sub["value"].astype(float).tolist())
    bw.close()
    return bw_path


# ---------------------------------------------------------------------------
# Kendall pairs, ATAC matrix, RNA matrix, gene-name mapping, Kendall tau-b
# ---------------------------------------------------------------------------

def read_narrowpeak(path):
    return read_bed_like(path, ["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])


def make_kendall_pairs(peaks, allputative):
    """make_kendall_pairs.R: peaks overlapping ABC E-G pairs -> chr,start,end,TargetGene,PeakName,PairName."""
    pd = _pd()
    ia, ib = keyed_overlap(peaks, allputative, ["chr"], ["chr"], closed=True)
    p = peaks.iloc[ia].reset_index(drop=True)
    out = pd.DataFrame({"chr": p["chr"].astype(str), "start": p["start"].astype("int64"), "end": p["end"].astype("int64"),
                        "TargetGene": allputative["TargetGene"].iloc[ib].astype(str).to_numpy()})
    out["PeakName"] = out["chr"] + "-" + out["start"].astype(str) + "-" + out["end"].astype(str)
    out["PairName"] = out["PeakName"] + "_" + out["TargetGene"]
    out = out.sort_values("PairName", kind="mergesort").drop_duplicates("PairName").reset_index(drop=True)
    return out


def choose_cells(rna_cells: "Sequence[str]", fragment_barcodes: "Optional[Sequence[str]]", max_cells: int = MAX_CELL_COUNT,
                 seed: int = 123) -> "List[str]":
    np = _np()
    cells = list(rna_cells)
    if fragment_barcodes is not None:
        fb = set(fragment_barcodes)
        cells = [c for c in dict.fromkeys(cells) if c in fb]
    if len(cells) > max_cells:
        rng = np.random.RandomState(seed)
        cells = [cells[i] for i in rng.choice(len(cells), max_cells, replace=False)]
    return cells


def build_atac_matrix(peaks, frag_path, cells: "Sequence[str]", mode: str = "fragment", chunksize: int = 5_000_000):
    """Signac FeatureMatrix: peaks x cells counts from a fragment file.

    mode 'fragment' (default, Signac via tabix): fragment [s+1, e] overlaps peak [start, end];
    mode 'insertion': count Tn5 cut sites s+1 and e falling in [start, end].
    """
    pd, np, sp = _pd(), _np(), _sp()
    cidx = {c: i for i, c in enumerate(cells)}
    P = pd.DataFrame({"chr": peaks["chr"].astype(str).to_numpy(), "start": peaks["start"].to_numpy(np.int64),
                      "end": peaks["end"].to_numpy(np.int64)})
    rows, cols = [], []
    for ch in read_fragments(frag_path, chunksize=chunksize):
        ch = ch[ch["barcode"].isin(cidx)]
        if len(ch) == 0:
            continue
        ci = ch["barcode"].map(cidx).to_numpy(np.int64)
        s = ch["start"].to_numpy(np.int64)
        e = ch["end"].to_numpy(np.int64)
        chrom = ch["chr"].astype(str).to_numpy()
        if mode == "fragment":
            A = pd.DataFrame({"chr": chrom, "start": s, "end": e + 1})
            Pb = P.assign(end=P["end"])
            cell_of = ci
        elif mode == "insertion":
            A = pd.DataFrame({"chr": np.concatenate([chrom, chrom]), "start": np.concatenate([s + 1, e]),
                              "end": np.concatenate([s + 2, e + 1])})
            Pb = P.assign(end=P["end"] + 1)
            cell_of = np.concatenate([ci, ci])
        else:
            raise SystemExit(f"unknown count mode {mode!r}")
        ia, ib = keyed_overlap(A, Pb, ["chr"], ["chr"], closed=False)
        rows.append(ib)
        cols.append(cell_of[ia])
    r = np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64)
    c = np.concatenate(cols) if cols else np.zeros(0, dtype=np.int64)
    M = sp.coo_matrix((np.ones(len(r), dtype=np.float64), (r, c)), shape=(len(P), len(cells))).tocsr()
    M.sum_duplicates()
    return M


def save_matrix(M, rows: "Sequence[str]", cols: "Sequence[str]", prefix: Path) -> "Dict[str, str]":
    import scipy.io  # type: ignore
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    mtx = Path(str(prefix) + ".mtx.gz")
    with gzip.open(mtx, "wb") as fh:
        scipy.io.mmwrite(fh, M.tocoo())
    rp, cp = Path(str(prefix) + ".peaks.txt"), Path(str(prefix) + ".barcodes.txt")
    rp.write_text("\n".join(rows) + "\n")
    cp.write_text("\n".join(cols) + "\n")
    print(f"Wrote: {mtx}")
    return {"matrix": str(mtx), "rows": str(rp), "cols": str(cp)}


def load_matrix(prefix):
    import scipy.io  # type: ignore
    prefix = str(prefix)
    for suf in (".mtx.gz", ".mtx"):
        if prefix.endswith(suf):
            prefix = prefix[: -len(suf)]
    M = scipy.io.mmread(prefix + ".mtx.gz" if Path(prefix + ".mtx.gz").exists() else prefix + ".mtx").tocsr()
    rows = Path(prefix + ".peaks.txt").read_text().split("\n")[:-1]
    cols = Path(prefix + ".barcodes.txt").read_text().split("\n")[:-1]
    return M, rows, cols


def read_rna_matrix(path):
    """RNA counts as a genes x cells CSC matrix: .h5ad/.h5 (cells x genes, transposed), .csv(.gz) with
    genes as rows, or a 10x directory (matrix.mtx, features/genes.tsv column 1, barcodes.tsv)."""
    pd, np, sp = _pd(), _np(), _sp()
    p = Path(path)
    if p.is_dir():
        import scipy.io  # type: ignore

        def pick(*names):
            for n in names:
                for suf in ("", ".gz"):
                    q = p / (n + suf)
                    if q.exists():
                        return q
            raise SystemExit(f"{p}: missing {names}")
        M = scipy.io.mmread(str(pick("matrix.mtx"))).tocsc()
        genes = pd.read_csv(pick("features.tsv", "genes.tsv"), sep="\t", header=None)[0].astype(str).tolist()
        cells = pd.read_csv(pick("barcodes.tsv"), sep="\t", header=None)[0].astype(str).tolist()
        return M.astype(np.float64), genes, cells
    if p.suffix in (".h5ad", ".h5"):
        import anndata  # type: ignore
        ad = anndata.read_h5ad(str(p))
        X = ad.X
        X = sp.csc_matrix(X.T) if sp.issparse(X) else sp.csc_matrix(np.asarray(X).T)
        return X.astype(np.float64), [str(g) for g in ad.var_names], [str(c) for c in ad.obs_names]
    df = pd.read_csv(p, index_col=0)
    return sp.csc_matrix(df.to_numpy(dtype=np.float64)), [str(g) for g in df.index], [str(c) for c in df.columns]


def rna_features(counts, genes: "Sequence[str]"):
    """compute_kendall.R: umi/cell counts, Seurat LogNormalize, per-gene TPM / detection / mean log-norm."""
    pd, np, sp = _pd(), _np(), _sp()
    X = sp.csc_matrix(counts)
    n_cells = X.shape[1]
    total = float(X.sum())
    colsum = np.asarray(X.sum(axis=0)).ravel()
    # Seurat LogNormalize evaluates log1p(x / colSum * 1e4) in this order; keeping it keeps ties bit-identical
    N = X.copy().astype(np.float64)
    col_of = np.repeat(np.arange(n_cells), np.diff(N.indptr))
    with np.errstate(divide="ignore", invalid="ignore"):
        N.data = np.log1p(N.data / colsum[col_of] * 1e4)
    N.data[~np.isfinite(N.data)] = 0.0
    N = sp.csr_matrix(N)
    rowsum = np.asarray(X.sum(axis=1)).ravel()
    det = np.asarray((X > 0).sum(axis=1)).ravel()
    exp = pd.DataFrame({"mean_log_normalized_rna": np.asarray(N.mean(axis=1)).ravel(),
                        "RnaDetectedPercent": det / n_cells if n_cells else np.nan,
                        "RnaPseudobulkTPM": rowsum / total * 1e6 if total > 0 else np.nan}, index=list(genes))
    return N, exp, total, n_cells


_GTF_NAME = re.compile(r'gene_name "([^"]*)"')
_GTF_ID = re.compile(r'gene_id "([^"]*)"')


def parse_gtf_genes(gtf_path):
    """gene_name / Ensembl gene_id (version suffix stripped) of every 'gene' line, distinct."""
    pd = _pd()
    names, ids = [], []
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "gene":
                continue
            mn, mi = _GTF_NAME.search(f[8]), _GTF_ID.search(f[8])
            names.append(mn.group(1) if mn else None)
            ids.append(re.sub(r"\.\d+$", "", mi.group(1)) if mi else None)
    return pd.DataFrame({"gene_ref_name": names, "Ensembl_ID": ids}).drop_duplicates().reset_index(drop=True)


def build_gene_key(gtf_genes, tss) -> "Dict[str, str]":
    """map_gene_names(): RNA (GTF) gene name -> TSS500 reference name via a 1:1 Ensembl-ID key."""
    t = tss[["name", "Ensembl_ID"]].rename(columns={"name": "abc_name"}).copy()
    t["_o"] = range(len(t))
    m = t.merge(gtf_genes, on="Ensembl_ID", how="left").sort_values("_o", kind="mergesort")
    m = m[m["gene_ref_name"].notna()]
    m = m[m.groupby("Ensembl_ID")["Ensembl_ID"].transform("size") == 1]
    key: Dict[str, str] = {}
    for ref, abc in zip(m["gene_ref_name"], m["abc_name"]):
        key.setdefault(str(ref), str(abc))
    return key


def kendall_tau_b_binary(x, y):
    """Reference tau-b for one expression vector and one 0/1 accessibility vector (Sheth, Qiu Fig. S1)."""
    np = _np()
    from scipy.stats import rankdata  # type: ignore
    x = np.asarray(x, dtype=float)
    y = (np.asarray(y) > 0).astype(float)
    n = len(x)
    r = rankdata(x)
    s = y.sum()
    S = 2.0 * float(r @ y) - s * (n + 1)
    _, t = np.unique(x, return_counts=True)
    n0 = n * (n - 1) / 2.0
    n1 = float((t * (t - 1) / 2.0).sum())
    n2 = s * (s - 1) / 2.0 + (n - s) * (n - s - 1) / 2.0
    d = (n0 - n1) * (n0 - n2)
    return S / math.sqrt(d) if d > 0 else float("nan")


def kendall_for_pairs(rna_norm, gene_idx, atac_bin, peak_idx):
    """Kendall tau-b for each (gene row of rna_norm, peak row of atac_bin) pair; both are row x cell CSR."""
    np, sp = _np(), _sp()
    from scipy.stats import rankdata  # type: ignore
    gene_idx = np.asarray(gene_idx, dtype=np.int64)
    peak_idx = np.asarray(peak_idx, dtype=np.int64)
    R = sp.csr_matrix(rna_norm)
    A = sp.csr_matrix(atac_bin)
    n = R.shape[1]
    out = np.full(len(gene_idx), np.nan)
    if n < 2 or len(gene_idx) == 0:
        return out
    n0 = n * (n - 1) / 2.0
    s_all = np.asarray(A.sum(axis=1)).ravel()
    order = np.argsort(gene_idx, kind="mergesort")
    g_sorted = gene_idx[order]
    bounds = np.flatnonzero(np.diff(g_sorted)) + 1
    for grp in np.split(order, bounds):
        g = gene_idx[grp[0]]
        x = R.getrow(g).toarray().ravel()
        _, t = np.unique(x, return_counts=True)
        n1 = float((t * (t - 1) / 2.0).sum())
        if n0 - n1 <= 0:
            continue
        r = rankdata(x)
        pk = peak_idx[grp]
        R1 = A[pk] @ r
        s = s_all[pk]
        S = 2.0 * R1 - s * (n + 1)
        n2 = s * (s - 1) / 2.0 + (n - s) * (n - s - 1) / 2.0
        den = (n0 - n1) * (n0 - n2)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[grp] = np.where(den > 0, S / np.sqrt(np.where(den > 0, den, 1.0)), np.nan)
    return out


def compute_kendall(pairs, atac_counts, atac_peaks: "Sequence[str]", atac_cells: "Sequence[str]", rna_path,
                    gene_key: "Optional[Dict[str, str]]"):
    """compute_kendall.R end to end. Returns (Pairs.Kendall frame, gene_expression_metrics frame, n_umi, n_cells)."""
    pd, np, sp = _pd(), _np(), _sp()
    counts, genes, rna_cells = read_rna_matrix(rna_path)
    pos = {c: i for i, c in enumerate(rna_cells)}
    missing = [c for c in atac_cells if c not in pos]
    if missing:
        raise SystemExit(f"{len(missing)} ATAC-matrix cells are not in the RNA matrix (e.g. {missing[:3]})")
    counts = counts[:, [pos[c] for c in atac_cells]]
    N, exp, n_umi, n_cells = rna_features(counts, genes)
    if gene_key is None:  # names already match the TSS reference
        row_sub = list(dict.fromkeys(genes))
        new_names = row_sub
    else:
        row_sub = [g for g in dict.fromkeys(genes) if g in gene_key]
        new_names = [gene_key[g] for g in row_sub]
    first_row = {}
    for i, g in enumerate(genes):
        first_row.setdefault(g, i)
    rows = [first_row[g] for g in row_sub]
    N_f = sp.csr_matrix(N)[rows]
    exp_f = exp.iloc[rows].copy()
    exp_f.index = new_names
    gex = pd.DataFrame({"TargetGene": new_names, "RNA_meanLogNorm": exp_f["mean_log_normalized_rna"].to_numpy(),
                        "RNA_pseudobulkTPM": exp_f["RnaPseudobulkTPM"].to_numpy(),
                        "RNA_percentCellsDetected": exp_f["RnaDetectedPercent"].to_numpy()})
    gene_row = {}
    for i, g in enumerate(new_names):
        gene_row.setdefault(g, i)
    peak_row = {}
    for i, p_ in enumerate(atac_peaks):
        peak_row.setdefault(p_, i)
    A = sp.csr_matrix(atac_counts)
    A = (A > 0).astype(np.float64)
    keep = pairs["TargetGene"].astype(str).isin(gene_row) & pairs["PeakName"].astype(str).isin(peak_row)
    pf = pairs[keep].reset_index(drop=True).copy()
    gi = pf["TargetGene"].astype(str).map(gene_row).to_numpy(np.int64)
    pi = pf["PeakName"].astype(str).map(peak_row).to_numpy(np.int64)
    pf["Kendall"] = kendall_for_pairs(N_f, gi, A, pi)
    ex = exp_f.iloc[gi]
    for c in GEX_COLS:
        pf[c] = ex[c].to_numpy()
    pf = pf[KENDALL_OUT_COLS].sort_values(["chr", "start", "end", "TargetGene"], kind="mergesort").reset_index(drop=True)
    return pf, gex, n_umi, n_cells


# ---------------------------------------------------------------------------
# ARC-E2G
# ---------------------------------------------------------------------------

def compute_arc(abc, kendall, abc_score_col: str = "powerlaw.Score"):
    """compute_arc_e2g.R: max overlapping Kendall per ABC pair, ABC x Kendall integration, RNA columns."""
    pd, np = _pd(), _np()
    k = kendall[kendall["Kendall"].notna()].reset_index(drop=True)
    ia, ib = keyed_overlap(abc, k, ["chr", "TargetGene"], ["chr", "TargetGene"], closed=True)
    kend = np.full(len(abc), np.nan)
    if len(ia):
        kend = group_aggregate(k["Kendall"].to_numpy(float)[ib], ia, len(abc), "max")
    out = pd.DataFrame({"chr": abc["chr"].to_numpy(), "start": abc["start"].to_numpy(), "end": abc["end"].to_numpy()})
    out["width"] = abc["end"].to_numpy(np.int64) - abc["start"].to_numpy(np.int64) + 1
    out["strand"] = "*"
    for c in abc.columns:
        if c not in ("chr", "start", "end"):
            out[c] = abc[c].to_numpy()
    out["Kendall"] = kend
    if abc_score_col not in abc.columns:
        raise SystemExit(f"ABC column {abc_score_col!r} missing from the ABC predictions")
    a = abc[abc_score_col].to_numpy(float)
    ok = np.isfinite(a) & (a > 0) & np.isfinite(kend)
    arc = a.copy()
    fit = {"n_used": int(ok.sum()), "r": None, "sd_lnABC": None, "sd_Kendall": None}
    if ok.sum() >= 3:
        la, kk = np.log(a[ok]), kend[ok]
        sd_la, sd_k = float(np.std(la, ddof=1)), float(np.std(kk, ddof=1))
        if sd_la > 0 and sd_k > 0:
            r = float(np.corrcoef(la, kk)[0, 1])
            arc[ok] = np.exp(la + kk * sd_la / sd_k * r)  # beta = 1 / r in the R source
            fit.update({"r": r, "sd_lnABC": sd_la, "sd_Kendall": sd_k})
    out["ARC.E2G.Score"] = arc
    g = k.drop_duplicates("TargetGene").set_index("TargetGene")
    for c in GEX_COLS:
        out[c] = abc["TargetGene"].map(g[c]).to_numpy() if c in g.columns else np.nan
    return out, fit


# ---------------------------------------------------------------------------
# ENCODE_rE2G feature tables as scE2G uses them
# ---------------------------------------------------------------------------

def _mid_int(start, end):
    np = _np()
    return ((np.asarray(start, dtype=float) + np.asarray(end, dtype=float)) / 2).astype(np.int64)


def num_candidate_enh_gene(pred):
    """gen_num_candidate_enh_gene.py: elements counted outward from the TSS per gene (0 at the TSS)."""
    pd, np = _pd(), _np()
    df = pd.DataFrame({"name": pred["name"].to_numpy(), "TargetGene": pred["TargetGene"].to_numpy(),
                       "midpoint": _mid_int(pred["start"], pred["end"]), "tss": pred["TargetGeneTSS"].to_numpy(float),
                       "_o": np.arange(len(pred))})
    df = df.sort_values(["TargetGene", "midpoint"], kind="mergesort")
    val = pd.Series(0, index=df.index, dtype=np.int64)
    up = df["midpoint"] < df["tss"]
    down = df["midpoint"] > df["tss"]
    u = df[up].sort_values("midpoint", ascending=False, kind="mergesort")
    val.loc[u.index] = u.groupby("TargetGene").cumcount().to_numpy() + 1
    d = df[down]
    val.loc[d.index] = d.groupby("TargetGene").cumcount().to_numpy() + 1
    df["NumCandidateEnhGene"] = val
    df = df.sort_values("_o")
    return df[["name", "TargetGene", "NumCandidateEnhGene"]].reset_index(drop=True)


def _count_overlaps_sorted(q_chr, q_start, q_end, t_chr, t_start, t_end, t_val=None):
    """Per query: number (and summed t_val) of targets on the chromosome with t.start < q.end and t.end > q.start."""
    np = _np()
    n = len(q_start)
    cnt = np.zeros(n, dtype=np.int64)
    tot = np.zeros(n) if t_val is not None else None
    q_chr, t_chr = np.asarray(q_chr).astype(str), np.asarray(t_chr).astype(str)
    q_start, q_end = np.asarray(q_start, dtype=np.int64), np.asarray(q_end, dtype=np.int64)
    q_end = np.maximum(q_end, q_start + 1)
    for c in np.unique(q_chr):
        qi = np.flatnonzero(q_chr == c)
        ti = np.flatnonzero(t_chr == c)
        if len(ti) == 0:
            continue
        ts, te = np.asarray(t_start, dtype=np.int64)[ti], np.asarray(t_end, dtype=np.int64)[ti]
        o1, o2 = np.argsort(ts, kind="mergesort"), np.argsort(te, kind="mergesort")
        a = np.searchsorted(ts[o1], q_end[qi], side="left")
        b = np.searchsorted(te[o2], q_start[qi], side="right")
        cnt[qi] = a - b
        if t_val is not None:
            v = np.asarray(t_val, dtype=float)[ti]
            c1 = np.concatenate([[0.0], np.cumsum(v[o1])])
            c2 = np.concatenate([[0.0], np.cumsum(v[o2])])
            tot[qi] = c1[a] - c2[b]
    return cnt, tot


def num_tss_enh_gene(pred, tss):
    """gen_num_tss_enh_gene.py: TSS500 windows overlapping the element-to-TSS span (bedtools intersect -c)."""
    pd, np = _pd(), _np()
    mid = _mid_int(pred["start"], pred["end"])
    tssp = pred["TargetGeneTSS"].to_numpy(float)
    start = pred["start"].to_numpy(np.int64).copy()
    new_end = (mid + pred["distance"].to_numpy(float)).astype(np.int64)
    down = tssp < mid
    new_end[down] = pred["end"].to_numpy(np.int64)[down]
    start[down] = tssp[down].astype(np.int64)
    cnt, _ = _count_overlaps_sorted(pred["chr"], start, new_end, tss["chr"], tss["start"], tss["end"])
    return pd.DataFrame({"name": pred["name"].to_numpy(), "gene": pred["TargetGene"].to_numpy(), "count": cnt})


def nearby_enhancers(pred, enhancer_list, chrom_sizes: "Optional[Dict[str, int]]", distance_bp: int = 5000):
    """gen_num_sum_nearby_enhancers.py: other candidate elements within +/- distance of each element midpoint.

    Returns (num frame name,count; sum frame name,sum) for elements with at least one neighbour."""
    pd, np = _pd(), _np()
    slim = pred[["chr", "start", "end", "name", "activity_base"]].drop_duplicates().reset_index(drop=True)
    mid = _mid_int(enhancer_list["start"], enhancer_list["end"])
    ws = np.maximum(mid - distance_bp, 0)
    we = mid + distance_bp
    if chrom_sizes:
        lim = enhancer_list["chr"].map(chrom_sizes).to_numpy(float)
        we = np.where(np.isfinite(lim), np.minimum(we, np.nan_to_num(lim, nan=0).astype(np.int64)), we)
    cnt, tot = _count_overlaps_sorted(enhancer_list["chr"], ws, we, slim["chr"], slim["start"], slim["end"],
                                      slim["activity_base"].to_numpy(float))
    # remove the element itself (same name) when it overlaps its own window
    q = pd.DataFrame({"qi": np.arange(len(enhancer_list)), "name": enhancer_list["name"].to_numpy(), "ws": ws, "we": we,
                      "qchr": enhancer_list["chr"].astype(str).to_numpy()})
    self_ = q.merge(slim.assign(chr=slim["chr"].astype(str)), on="name")
    hit = (self_["qchr"] == self_["chr"]) & (self_["start"] < self_["we"]) & (self_["end"] > self_["ws"])
    self_ = self_[hit]
    if len(self_):
        sc = self_.groupby("qi").size()
        sa = self_.groupby("qi")["activity_base"].sum()
        cnt[sc.index.to_numpy()] -= sc.to_numpy()
        tot[sa.index.to_numpy()] -= sa.to_numpy()
    names = enhancer_list["name"].to_numpy()
    has = cnt > 0
    num = pd.DataFrame({"name": names[has], "count": cnt[has]})
    sm = pd.DataFrame({"name": names[has], "sum": tot[has]})
    # bedtools groupby output is keyed by enhancer name (one row per name)
    return (num.groupby("name", sort=True, as_index=False)["count"].sum(),
            sm.groupby("name", sort=True, as_index=False)["sum"].sum())


def input_features_of(ft) -> "List[str]":
    vals = list(ft["input_col"]) + list(ft["second_input"])
    out = []
    for v in vals:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            continue
        v = str(v)
        if v in ("", "NA", "nan") or v in out:
            continue
        out.append(v)
    return out


def biosample_feature_table(models: "Sequence[dict]"):
    """combine_feature_tables_apply.R: union of model feature tables + scE2G ARC rows, distinct."""
    pd = _pd()
    df = pd.concat([m["feature_table"] for m in models], ignore_index=True)
    if df["feature"].isin(["ARC.E2G.Score", "Kendall"]).any():
        arc = pd.DataFrame(ARC_FEATURE_ROWS, columns=FT_COLS)
        if (df["feature"] == "ABC.Score").any():
            arc = arc[arc["feature"] != "ABC.Score"]
        df = pd.concat([df, arc], ignore_index=True)
    return df.drop_duplicates().reset_index(drop=True)


def needs_arc(ft) -> bool:
    cols = set(ft["input_col"].astype(str)) | set(ft["second_input"].astype(str))
    return bool(cols & {"Kendall", "ARC.E2G.Score"})


def external_features_config(arc_file: "Optional[str]", user_efc=None):
    """format_external_features_config_sc.R."""
    pd = _pd()
    efc = user_efc if user_efc is not None else pd.DataFrame(columns=EFC_COLS)
    if arc_file:
        if "Pairs.Kendall.tsv.gz" in str(arc_file):
            rows = [("Kendall", "Kendall", "max", "overlap", str(arc_file))]
        else:
            rows = [(c, c, a, "overlap", str(arc_file)) for c, a in ARC_EXTERNAL_ROWS]
        efc = pd.concat([efc, pd.DataFrame(rows, columns=EFC_COLS)], ignore_index=True)
    return efc


def activity_only_features(abc, ft, numcand=None, numtss=None, num5=None, sum5=None, gene_classes=None):
    """activity_only_features.R (pinned ENCODE_rE2G)."""
    pd = _pd()
    abc = abc.copy()
    if "ABC.Score.Numerator" in abc.columns:
        abc["ABC.Numerator"] = abc["ABC.Score.Numerator"]
        abc["ABC.Denominator"] = abc["ABC.Score"] / abc["ABC.Numerator"]
    inputs = input_features_of(ft)
    core = [c for c in ACTIVITY_CORE_COLS if c in abc.columns and c not in inputs]
    out = abc[core].copy()
    for c in inputs:  # select(name, TargetGene, any_of(input_features)) -> input-feature order
        if c in abc.columns and c not in out.columns:
            out[c] = abc[c].to_numpy()
    key = out["name"].astype(str) + "\x1f" + out["TargetGene"].astype(str)
    if numcand is not None and "numCandidateEnhGene" in inputs:
        m = dict(zip(numcand.iloc[:, 0].astype(str) + "\x1f" + numcand.iloc[:, 1].astype(str), numcand.iloc[:, 2]))
        out["numCandidateEnhGene"] = key.map(m).to_numpy()
    if numtss is not None and "numTSSEnhGene" in inputs:
        m = dict(zip(numtss.iloc[:, 0].astype(str) + "\x1f" + numtss.iloc[:, 1].astype(str), numtss.iloc[:, 2]))
        out["numTSSEnhGene"] = key.map(m).to_numpy()
    if num5 is not None and "numNearbyEnhancers" in inputs:
        out["numNearbyEnhancers"] = out["name"].astype(str).map(dict(zip(num5.iloc[:, 0].astype(str), num5.iloc[:, 1]))).to_numpy()
    if sum5 is not None and "sumNearbyEnhancers" in inputs:
        out["sumNearbyEnhancers"] = out["name"].astype(str).map(dict(zip(sum5.iloc[:, 0].astype(str), sum5.iloc[:, 1]))).to_numpy()
    if gene_classes is not None:
        gc = gene_classes.drop_duplicates("TargetGene").set_index("TargetGene")
        for c in inputs:
            if c in gc.columns and c not in out.columns:
                out[c] = out["TargetGene"].map(gc[c]).to_numpy()
    return out


def merge_external_features(df, ft, efc, sources: "Optional[Dict[str, Any]]" = None):
    """merge_external_features.R: add features from external sources by (chr, TargetGene) overlap or gene."""
    pd, np = _pd(), _np()
    needed = [c for c in input_features_of(ft) if c not in df.columns]
    efc = efc[efc["input_col"].isin(needed)]
    for src in list(dict.fromkeys(efc["source_file"].dropna())):
        this = efc[efc["source_file"] == src]
        s = sources[src] if sources and src in sources else read_tsv(src)
        s = s.rename(columns=dict(zip(this["source_col"], this["input_col"])))
        cols = list(this["input_col"])
        if (this["join_by"] == "overlap").any():
            ia, ib = keyed_overlap(df, s, ["chr", "TargetGene"], ["chr", "TargetGene"], closed=True)
            merged_mask = np.zeros(len(df), dtype=bool)
            merged_mask[ia] = True
            new = df.copy()
            for c, how in zip(cols, this["aggregate_function"]):
                new[c] = group_aggregate(pd.to_numeric(s[c], errors="coerce").to_numpy(float)[ib], ia, len(df), how)
            new = pd.concat([new[merged_mask], new[~merged_mask]], ignore_index=True)
            df = new.sort_values(["chr", "start", "end", "TargetGene"], kind="mergesort").reset_index(drop=True)
        elif (this["join_by"] == "TargetGene").all():
            g = s.drop_duplicates("TargetGene").set_index("TargetGene")
            for c in cols:
                df[c] = df["TargetGene"].map(g[c]).to_numpy()
    return df


def get_fill_values(features, ft) -> "Dict[str, float]":
    """get_fill_values.R: fill_value per feature, 'mean' -> mean of the finite values."""
    np = _np()
    out: Dict[str, float] = {}
    for f, v in zip(ft["feature"], ft["fill_value"]):
        if f not in features.columns or f in out:
            continue
        if str(v).strip().lower() == "mean":
            x = pd_numeric(features[f])
            x = x[np.isfinite(x)]
            out[f] = float(x.mean()) if len(x) else float("nan")
        else:
            try:
                out[f] = float(v)
            except (TypeError, ValueError):
                out[f] = float("nan")
    return out


def pd_numeric(s):
    pd, np = _pd(), _np()
    if s.dtype == bool:
        return s.astype(float).to_numpy()
    return pd.to_numeric(s, errors="coerce").to_numpy(float)


def final_features(df, ft):
    """gen_final_features.R: check inputs, interaction terms, rename input_col -> feature, fill NAs."""
    df = df.copy()
    inputs = input_features_of(ft)
    missing = [c for c in inputs if c not in df.columns]
    if missing:
        raise SystemExit(f"required feature columns missing: {missing}")
    sec = ft["second_input"].astype(str)
    intx = ft[~sec.isin(["", "NA", "nan", "None"])]
    for f, a, b in zip(intx["feature"], intx["input_col"], intx["second_input"]):
        df[f] = pd_numeric(df[a]) * pd_numeric(df[b])
    single = ft[~ft["feature"].isin(set(intx["feature"]))]
    for f, a in zip(single["feature"], single["input_col"]):
        if a != f and a in df.columns:
            df = df.rename(columns={a: f})
    for f, v in get_fill_values(df, ft).items():
        if df[f].isna().any():
            df[f] = df[f].where(df[f].notna(), v)
    return df


# ---------------------------------------------------------------------------
# Model application: E2G.Score, cv, qnorm, TPM filter, threshold, bedpe
# ---------------------------------------------------------------------------

def calculate_quantiles(scores):
    """run_e2g_cv.calculate_quantiles: zeros -> quantile 0; others (rank - 1) / (N - 1), zeros ranked lowest."""
    np = _np()
    s = scores.replace(0, np.nan)
    ranks = s.rank(method="average", na_option="top")
    q = (ranks - 1) / (len(scores) - 1) if len(scores) > 1 else ranks * 0
    q[scores == 0] = 0
    return q


def read_qnorm_reference(path):
    ref = read_tsv(path)
    if not {"quantile", "reference_score"} <= set(ref.columns):
        raise SystemExit(f"{path}: expected columns quantile, reference_score")
    return ref


def qnorm_scores(df, ref, crispr_benchmarking: bool = False):
    np = _np()
    from scipy import interpolate  # type: ignore
    f = interpolate.interp1d(ref["quantile"].to_numpy(float), ref["reference_score"].to_numpy(float), kind="linear",
                             fill_value="extrapolate")
    df[SCORE_BASE + ".qnorm"] = np.clip(f(calculate_quantiles(df[SCORE_BASE]).to_numpy(float)), 0, None)
    if crispr_benchmarking:
        df[SCORE_BASE + ".cv.qnorm"] = np.clip(f(calculate_quantiles(df[SCORE_BASE + ".cv"]).to_numpy(float)), 0, None)
    return df


def filter_by_tpm(df, tpm_threshold: float, crispr_benchmarking: bool = False):
    """run_e2g_cv.filter_by_tpm: keep the unfiltered qnorm score as *.ignoreTPM, zero genes below the TPM threshold."""
    np = _np()
    if "RNA_pseudobulkTPM" not in df.columns or float(tpm_threshold) == 0:
        return df
    tpm = pd_numeric(df["RNA_pseudobulkTPM"])
    low = ~(tpm >= float(tpm_threshold))
    df[SCORE_BASE + ".qnorm.ignoreTPM"] = df[SCORE_BASE + ".qnorm"]
    df[SCORE_BASE + ".qnorm"] = np.where(low, 0.0, df[SCORE_BASE + ".qnorm.ignoreTPM"])
    if crispr_benchmarking:
        df[SCORE_BASE + ".cv.qnorm.ignoreTPM"] = df[SCORE_BASE + ".cv.qnorm"]
        df[SCORE_BASE + ".cv.qnorm"] = np.where(low, 0.0, df[SCORE_BASE + ".cv.qnorm.ignoreTPM"])
    return df


def apply_model(features, model: dict, crispr_benchmarking: bool = False, epsilon: float = EPSILON,
                qnorm: bool = True):
    """run_e2g_cv.py: E2G.Score (+ .cv), qnorm against the model reference, TPM filter."""
    np = _np()
    df = features.replace([np.inf, -np.inf], np.nan).fillna(0)
    feats = model["features"]
    missing = [f for f in feats if f not in df.columns]
    if missing:
        raise SystemExit(f"model {model['name']}: features missing from the feature table: {missing}")
    M = transform_features(df, feats, epsilon, model.get("polynomial", False))
    df[SCORE_BASE] = logistic_score(M, model["full"])
    info = {"cv_fallback_chromosomes": []}
    if crispr_benchmarking:
        cv = np.full(len(df), np.nan)
        chrom = df["chr"].astype(str).to_numpy()
        for c in np.unique(chrom):
            idx = np.flatnonzero(chrom == c)
            w = model["cv"].get(c)
            if w is None:
                info["cv_fallback_chromosomes"].append(c)
                w = model["full"]
            cv[idx] = logistic_score(M[idx], w)
        df[SCORE_BASE + ".cv"] = cv
    if qnorm:
        if model.get("qnorm") is None:
            raise SystemExit(f"qnorm reference for {model['name']} not found at {model.get('qnorm_expected')}; "
                             "run `igvfagent sce2g-pipeline setup` (or pass --no-qnorm)")
        df = qnorm_scores(df, read_qnorm_reference(model["qnorm"]), crispr_benchmarking)
    else:
        df[SCORE_BASE + ".qnorm"] = df[SCORE_BASE]
        if crispr_benchmarking:
            df[SCORE_BASE + ".cv.qnorm"] = df[SCORE_BASE + ".cv"]
    df = filter_by_tpm(df, model.get("tpm_threshold", 0.0), crispr_benchmarking)
    return df, info


def threshold_predictions(pred, threshold: float, score_column: str = FINAL_SCORE_COL, include_self_promoter: bool = True):
    """threshold_e2g_predictions.py."""
    f = pred[pd_numeric(pred[score_column]) >= float(threshold)]
    self_prom = f["isSelfPromoter"].map(_truthy) if "isSelfPromoter" in f.columns else False
    if include_self_promoter:
        return f[(f["class"] != "promoter") | self_prom]
    return f[f["class"] != "promoter"]


def predictions_bedpe(pred, score_column: str = FINAL_SCORE_COL):
    """process_model_output.write_connections_bedpe_format (pinned: no TargetGeneIsExpressed filter)."""
    pd = _pd()
    p = pred.drop_duplicates()
    return pd.DataFrame({"chr1": p["chr"], "x1": p["start"], "x2": p["end"], "chr2": p["chr"], "y1": p["TargetGeneTSS"],
                         "y2": p["TargetGeneTSS"], "name": p["TargetGene"].astype(str) + "_" + p["name"].astype(str),
                         "score": p[score_column], "strand1": ".", "strand2": "."})


def qnorm_reference_from_scores(scores, n_scores: int = 10000):
    """run_e2g_cv_qnormref.get_qnorm_quantiles."""
    pd, np = _pd(), _np()
    q = np.linspace(0, 1, int(n_scores), endpoint=True)
    return pd.DataFrame({"quantile": q, "reference_score": np.quantile(np.asarray(scores, dtype=float), q)})


# ---------------------------------------------------------------------------
# Gene / element lists, per-cluster statistics, QC
# ---------------------------------------------------------------------------

def element_gene_lists(abc_gene_list, abc_element_list, pred, cluster: str, tpm_threshold: float, gex=None):
    """generate_element_gene_lists.R."""
    pd, np = _pd(), _np()
    gcols = [c for c in ("normalizedATAC_prom", "ubiqExpressed") if c in pred.columns]
    pg = pred[["TargetGene", "TargetGeneEnsembl_ID"] + gcols].drop_duplicates().rename(
        columns={"TargetGene": "name", "TargetGeneEnsembl_ID": "Ensembl_ID"})
    drop = [c for c in abc_gene_list.columns if c in ("is_ue", "Expression", "cellType") or "RPKM" in c]
    gl = abc_gene_list.drop(columns=drop)
    gl["removed_by_promoter_activity"] = ~gl["name"].isin(set(pg["name"]))
    gl = gl.merge(pg, on=["name", "Ensembl_ID"], how="left")
    gl["CellType"] = cluster
    if "RNA_pseudobulkTPM" in pred.columns and gex is not None:
        rc = [c for c in ("RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected") if c in gex.columns]
        g = gex[["TargetGene"] + rc].rename(columns={"TargetGene": "name"})
        gl = gl.merge(g, on="name", how="left")
        tpm = pd_numeric(gl["RNA_pseudobulkTPM"])
        gl["removed_by_TPM"] = np.where(np.isnan(tpm), False, tpm < float(tpm_threshold))
        gl["in_RNA_matrix"] = ~np.isnan(tpm)
        gl["considered_for_predictions"] = (~gl["removed_by_promoter_activity"]) & (~gl["removed_by_TPM"].astype(bool)) & gl["in_RNA_matrix"]
    else:
        gl["considered_for_predictions"] = ~gl["removed_by_promoter_activity"]
    ecols = [c for c in ("numNearbyEnhancers", "sumNearbyEnhancers", "normalizedATAC_enh") if c in pred.columns]
    pe = pred[["chr", "start", "end"] + ecols].drop_duplicates()
    el = abc_element_list.merge(pe, on=["chr", "start", "end"], how="left")
    return gl, el


def stats_from_predictions(pred_full, pred_thr, score_column: str = FINAL_SCORE_COL) -> dict:
    """get_stats_per_cluster.R (thresholded predictions are filtered to class != promoter first)."""
    np = _np()
    thr = pred_thr[pred_thr["class"] != "promoter"]
    eg = thr[["chr", "start", "end", "TargetGene"]]
    enh = eg[["chr", "start", "end"]].drop_duplicates()
    dist_col = "distanceToTSS" if "distanceToTSS" in thr.columns else "distance"

    def mean(x):
        x = np.asarray(x, dtype=float)
        return float(x.mean()) if len(x) else float("nan")
    return {
        "n_enh_elements": int(len(enh)), "n_genes_with_enh": int(eg["TargetGene"].nunique()), "n_enh_gene_links": int(len(eg)),
        "mean_genes_per_enh": mean(eg.groupby(["chr", "start", "end"]).size()) if len(eg) else float("nan"),
        "mean_enh_per_gene": mean(eg.groupby("TargetGene").size()) if len(eg) else float("nan"),
        "mean_dist_to_tss": mean(thr[dist_col]) if len(thr) else float("nan"),
        "mean_enh_width": mean(enh["end"] - enh["start"]) if len(enh) else float("nan"),
        "n_genes_active_promoter": int(pred_full["TargetGene"].nunique()),
        "n_genes_not_expressed": int(pred_full.loc[pd_numeric(pred_full[score_column]) == 0, "TargetGene"].nunique()),
    }


def model_colors(model_names: "Sequence[str]") -> "Dict[str, str]":
    out, fb = {}, iter(FALLBACK_COLORS)
    for m in model_names:
        out[m] = MODEL_SPECS.get(m, {}).get("color") or next(fb, GREY)
    return out


def qc_summary(stats, reference=None):
    """plot_all_qc_stats.R + qc_report.Rmd: dataset metrics, warnings, comparison with reference clusters."""
    pd, np = _pd(), _np()
    ds = stats[["cluster", "fragments_total", "cell_count", "umi_count"]].drop_duplicates()
    with np.errstate(divide="ignore", invalid="ignore"):
        ds = ds.assign(frag_per_cell=ds["fragments_total"] / ds["cell_count"].where(ds["cell_count"] > 0),
                       umi_per_cell=ds["umi_count"] / ds["cell_count"].where(ds["cell_count"] > 0))
    warnings_ = []
    for col, thr in QC_THRESHOLDS.items():
        sub = ds if col == "fragments_total" else ds[ds[col] > 0]
        bad = sorted(sub.loc[sub[col] < thr, "cluster"].astype(str).unique())
        if bad:
            warnings_.append(f"WARNING: the following identifier(s) have less than {thr:.0e} in {col}: {', '.join(bad)}")
    for col, what in (("cell_count", "Cell counts"), ("umi_count", "RNA UMI counts")):
        empty = sorted(ds.loc[ds[col] == 0, "cluster"].astype(str).unique())
        if empty:
            warnings_.append(f"{what} are not available for: {', '.join(empty)}")
    comp = []
    if reference is not None and len(reference):
        ref = reference[reference["fragments_total"] >= QC_THRESHOLDS["fragments_total"]]
        metrics = [m for m in STATS_COLS[:9]]
        for _, row in stats.iterrows():
            r = ref[ref["model_name"] == row["model_name"]]
            if len(r) == 0:
                continue
            for m in metrics:
                v = float(row[m]) if pd.notna(row[m]) else float("nan")
                rv = pd.to_numeric(r[m], errors="coerce").dropna().to_numpy(float)
                comp.append({"cluster": row["cluster"], "model_name": row["model_name"], "metric": m, "value": v,
                             "reference_n": len(rv), "reference_median": float(np.median(rv)) if len(rv) else float("nan"),
                             "reference_min": float(rv.min()) if len(rv) else float("nan"),
                             "reference_max": float(rv.max()) if len(rv) else float("nan"),
                             "reference_percentile": float((rv <= v).mean() * 100) if len(rv) and math.isfinite(v) else float("nan")})
    return ds, warnings_, pd.DataFrame(comp)


def qc_figures(stats, ds, reference, out: Path) -> "List[Path]":
    plt = _plt()
    if plt is None:
        return []
    np = _np()
    figs = []
    colors = model_colors(list(stats["model_name"].unique()))
    x_vars = ["fragments_total"] + (["umi_count", "cell_count"] if stats["umi_count"].sum() > 0 else [])
    if stats["umi_count"].sum() > 0:
        fig, axes = plt.subplots(1, 2, figsize=(7, 3.4))
        for ax, y in zip(axes, ("frag_per_cell", "umi_per_cell")):
            ax.scatter(ds["cell_count"].clip(lower=1), ds[y], s=18, color="#1c2a43")
            ax.set_xscale("log")
            _style(ax, "", f"{QC_LABELS['cell_count']}\nN = {len(ds)}", QC_LABELS[y])
        figs.append(_save(fig, out / "dataset_metrics.png"))
    else:
        fig, ax = plt.subplots(figsize=(3.5, 3.4))
        ax.scatter(np.random.RandomState(0).uniform(-0.2, 0.2, len(ds)), ds["fragments_total"], s=18, color="#435369")
        ax.set_xticks([])
        _style(ax, "", "", QC_LABELS["fragments_total"])
        figs.append(_save(fig, out / "dataset_metrics.png"))
    sets = {"enhancer_metrics": ["n_enh_gene_links", "n_enh_elements"], "eg_metrics": ["mean_genes_per_enh", "mean_enh_per_gene"],
            "gene_metrics": ["n_genes_with_enh", "n_genes_active_promoter"] + (["n_genes_not_expressed"] if "umi_count" in x_vars else []),
            "distance_and_size_metrics": ["mean_dist_to_tss", "mean_enh_width"]}
    for fname, ys in sets.items():
        fig, axes = plt.subplots(len(ys), len(x_vars), figsize=(3 * len(x_vars) + 0.5, 2.8 * len(ys)), squeeze=False)
        for i, y in enumerate(ys):
            for j, x in enumerate(x_vars):
                ax = axes[i][j]
                ax.axvline(QC_THRESHOLDS[x], color=GREY, ls="--", lw=0.8)
                for m, sub in stats.groupby("model_name"):
                    ax.scatter(sub[x].clip(lower=1), sub[y], s=16, color=colors[m], label=m)
                ax.set_xscale("log")
                _style(ax, "", QC_LABELS[x], QC_LABELS[y])
        axes[0][0].legend(fontsize=6, frameon=False)
        figs.append(_save(fig, out / f"{fname}.png"))
    keys = [k for k in list(QC_LABELS)[5:]]
    ncol = 3
    nrow = int(math.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.4 * nrow), squeeze=False)
    for i, k in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        models = list(stats["model_name"].unique())
        for yi, m in enumerate(models):
            if reference is not None and len(reference) and k in reference.columns:
                rv = reference.loc[(reference["model_name"] == m) & (reference["fragments_total"] >= 2e6), k].dropna()
                if len(rv):
                    ax.scatter(rv, np.full(len(rv), yi) + np.random.RandomState(1).uniform(-0.15, 0.15, len(rv)),
                               s=6, color="#C5CAD7")
            v = stats.loc[stats["model_name"] == m, k].dropna()
            ax.scatter(v, np.full(len(v), yi), s=22, color=colors[m], zorder=3)
        ax.set_yticks(range(len(models)))
        ax.set_yticklabels(models, fontsize=6)
        _style(ax, "", QC_LABELS[k], "")
    for j in range(len(keys), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    figs.append(_save(fig, out / "qc_metric_distributions.png"))
    return figs


# ---------------------------------------------------------------------------
# CRISPR benchmark: feature overlap, training_functions, performance summary
# ---------------------------------------------------------------------------

def crispr_features(pred, crispr, ft, tss, nafill: bool = True, apply_mode: bool = True):
    """merge_features_with_crispr_data_apply.R (apply_mode) / overlap_features_with_crispr_data.R (training).

    Returns (merged table, CRISPR rows without any overlapping prediction)."""
    pd, np = _pd(), _np()
    crispr = crispr.drop(columns=[c for c in ("pair_uid", "merged_uid", "merged_start", "merged_end") if c in crispr.columns])
    cfg = ft[["feature", "aggregate_function", "fill_value"]].copy()
    if apply_mode:
        cfg = pd.concat([cfg, pd.DataFrame({"feature": SCORE_ROWS, "aggregate_function": "mean", "fill_value": "0"})],
                        ignore_index=True)
    cfg = cfg[cfg["feature"].isin(pred.columns)]
    agg = {}
    for f, a in zip(cfg["feature"], cfg["aggregate_function"]):
        agg.setdefault(f, a)
    feats = list(agg)
    crispr = crispr[crispr["measuredGeneSymbol"].isin(set(tss["name"]))].reset_index(drop=True)
    ia, ib = keyed_overlap(crispr, pred, ["chrom", "measuredGeneSymbol"], ["chr", "TargetGene"], closed=True,
                           a_cols=("chromStart", "chromEnd"))
    hit = np.zeros(len(crispr), dtype=bool)
    hit[ia] = True
    merged = crispr[hit].copy()
    uniq = np.flatnonzero(hit)
    remap = np.full(len(crispr), -1)
    remap[uniq] = np.arange(len(uniq))
    for f in feats:
        merged[f] = group_aggregate(pd_numeric(pred[f])[ib], remap[ia], len(uniq), agg[f]) if len(uniq) else []
    missing = crispr[~hit].copy()
    if apply_mode:
        for f in feats:
            missing[f] = np.nan
        out = pd.concat([merged, missing], ignore_index=True)
        out = out.sort_values(["dataset", "chrom", "chromStart", "chromEnd", "measuredGeneSymbol"], kind="mergesort")
    else:
        out = merged.reset_index(drop=True)
    t = tss[["name", "start", "end"]].drop_duplicates("name").set_index("name")
    s_ref, e_ref = out["measuredGeneSymbol"].map(t["start"]), out["measuredGeneSymbol"].map(t["end"])
    center = (out["chromStart"] + out["chromEnd"]) / 2
    by_crispr = (center - (out["startTSS"] + out["endTSS"]) / 2).abs()
    by_ref = (center - (s_ref + e_ref) / 2).abs()
    recomputed = by_crispr.where(s_ref.isna(), by_ref)
    if apply_mode:
        out["distanceToTSS"] = recomputed
        out["distance"] = out["distanceToTSS"]
    elif "distanceToTSS" in out.columns:
        out["distanceToTSS"] = out["distanceToTSS"].where(out["distanceToTSS"].notna(), recomputed)
    if nafill:
        fills = get_fill_values(pred, cfg.assign(nice_name="", input_col="", second_input=""))
        for f, v in fills.items():
            if f in out.columns:
                out[f] = out[f].where(out[f].notna(), v)
    return out.reset_index(drop=True), missing.reset_index(drop=True)


def process_crispr_for_training(df, tss):
    """process_crispr_data.R: rename to chr/start/end/TargetGene, drop Regulated NA, genes in the TSS universe."""
    df = df.rename(columns={"chrom": "chr", "chromStart": "start", "chromEnd": "end", "measuredGeneSymbol": "TargetGene"})
    df = df[df["Regulated"].notna() & (df["Regulated"].astype(str) != "NA")]
    return df[df["TargetGene"].isin(set(tss["name"]))].reset_index(drop=True)


def pr_curve_modified(y_true, y_pred):
    from sklearn.metrics import precision_recall_curve  # type: ignore
    p, r, t = precision_recall_curve(y_true, y_pred)
    return p[1:], r[1:], t


def auc_mod(recall, precision) -> float:
    from sklearn.metrics import auc  # type: ignore
    return 0.0 if len(recall) < 2 else float(auc(recall, precision))


def statistic_aupr(y_true, y_pred) -> float:
    p, r, _ = pr_curve_modified(y_true, y_pred)
    return auc_mod(r, p)


def statistic_delta_aupr(y_true, y_full, y_ablated) -> float:
    return statistic_aupr(y_true, y_ablated) - statistic_aupr(y_true, y_full)


def statistic_precision(y_true, y_pred, target: float = 0.7) -> float:
    np = _np()
    p, r, _ = pr_curve_modified(y_true, y_pred)
    return float(p[np.argsort(np.abs(r - target))[0]])


def threshold_at_target_recall(y_true, y_pred, target: float):
    np = _np()
    p, r, t = pr_curve_modified(y_true, y_pred)
    if np.max(r) > target:
        return float(t[np.argsort(np.abs(r - target))[0] + 1])
    return None


def statistic_precision_at_threshold(y_true, y_pred, threshold) -> float:
    np = _np()
    if threshold is None:
        return 0.0
    p, r, t = pr_curve_modified(y_true, y_pred)
    return float(p[np.argsort(np.abs(t - threshold))[0] - 1])


def statistic_recall_at_threshold(y_true, y_pred, threshold) -> float:
    np = _np()
    p, r, t = pr_curve_modified(y_true, y_pred)
    if threshold < np.max(t):
        return float(r[np.argsort(np.abs(t - threshold))[0] - 1])
    return 0.0


def bootstrap_pvalue(delta: float, boot_distribution) -> float:
    np = _np()
    d = np.asarray(boot_distribution) - np.mean(boot_distribution)
    return float(np.mean(np.abs(d) >= abs(delta)))


def _bootstrap(stat, y, s, n_boot: int, method: str, seed: int):
    np = _np()
    from scipy import stats as st  # type: ignore
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = st.bootstrap((y, s), stat, n_resamples=n_boot, paired=True, confidence_level=0.95, method=method,
                           vectorized=False, random_state=np.random.RandomState(seed))
    return float(np.mean(res.bootstrap_distribution)), float(res.confidence_interval[0]), float(res.confidence_interval[1])


def performance_summary(cluster: str, model_name: str, model_threshold: float, cf, n_boot: int = 1000,
                        method: str = "BCa", seed: int = 0) -> dict:
    """benchmark_performance.performance_summary."""
    np = _np()
    cf = cf.copy()
    reg = cf["Regulated"].map(_truthy)
    pos, neg = cf[reg], cf[~reg]

    def zero(frame, cols):
        m = np.zeros(len(frame), dtype=bool)
        for c in cols:
            if c in frame.columns:
                m |= pd_numeric(frame[c]) == 0
        return int(m.sum())
    if model_name == "distanceToTSS":
        dcol = [c for c in ("distance", "distanceToTSS") if c in cf.columns][0]
        n0p = n0n = 0
        score = -pd_numeric(cf[dcol])
    elif model_name == "scATAC_ABC":
        n0p, n0n = zero(pos, ["ABC.Score"]), zero(neg, ["ABC.Score"])
        score = pd_numeric(cf["ABC.Score"])
    elif model_name.startswith("multiome"):
        n0p, n0n = zero(pos, ["ARC.E2G.Score", "E2G.Score.cv"]), zero(neg, ["ARC.E2G.Score", "E2G.Score.cv"])
        col = "E2G.Score.cv.qnorm.ignoreTPM" if model_name == "multiome_powerlaw_v3_noTPMfilter" else "E2G.Score.cv.qnorm"
        score = pd_numeric(cf[col])
    else:
        n0p, n0n = zero(pos, ["ABC.Score"]), zero(neg, ["ABC.Score"])
        score = pd_numeric(cf["E2G.Score.cv.qnorm"])
    y = reg.to_numpy().astype(np.int64)
    s = np.nan_to_num(score, nan=0.0)
    row = {"cluster": cluster, "model": model_name}
    row["AUPRC"], row["AUPRC_95CI_low"], row["AUPRC_95CI_high"] = _bootstrap(statistic_aupr, y, s, n_boot, method, seed)
    for pct in (70, 50):
        t = threshold_at_target_recall(y, s, pct / 100)
        if t is None:
            vals = (0.0, 0.0, 0.0)
        else:
            vals = _bootstrap(lambda a, b, t=t: statistic_precision_at_threshold(a, b, t), y, s, n_boot, method, seed)
        row[f"precision_{pct}_pct_recall"], row[f"precision_{pct}_pct_recall_95CI_low"], row[f"precision_{pct}_pct_recall_95CI_high"] = vals
        row[f"threshold_{pct}_pct_recall"] = t
    mt = float(model_threshold)
    row["precision_model_threshold"], row["precision_model_threshold_95CI_low"], row["precision_model_threshold_95CI_high"] = \
        _bootstrap(lambda a, b: statistic_precision_at_threshold(a, b, mt), y, s, n_boot, method, seed)
    row["recall_model_threshold"], row["recall_model_threshold_95CI_low"], row["recall_model_threshold_95CI_high"] = \
        _bootstrap(lambda a, b: statistic_recall_at_threshold(a, b, mt), y, s, n_boot, method, seed)
    row["pct_missing_pos"] = n0p / len(pos) if len(pos) else float("nan")
    row["pct_missing_neg"] = n0n / len(neg) if len(neg) else float("nan")
    row["pct_missing_total"] = (n0p + n0n) / len(cf) if len(cf) else float("nan")
    return row


def benchmark_table(entries: "Sequence[dict]", n_boot: int = 1000, method: str = "BCa", seed: int = 0):
    """benchmark_performance.main: entries = [{cluster, model_name, model_threshold, frame}]."""
    pd = _pd()
    entries = list(entries)
    names = [e["model_name"] for e in entries]
    sc = [m for m in dict.fromkeys(names) if m in ("multiome_powerlaw_v3", "scATAC_powerlaw_v3", "multiome_megamap_v3")]
    extra = []
    if sc:
        extra += [dict(e, model_name="scATAC_ABC", model_threshold=ABC_BENCHMARK_THRESHOLD) for e in entries if e["model_name"] == sc[0]]
    if "multiome_powerlaw_v3" in names:
        extra += [dict(e, model_name="multiome_powerlaw_v3_noTPMfilter") for e in entries if e["model_name"] == "multiome_powerlaw_v3"]
    rows = [performance_summary(e["cluster"], e["model_name"], e["model_threshold"], e["frame"], n_boot, method, seed)
            for e in entries + extra]
    if entries:
        rows.append(performance_summary("None", "distanceToTSS", DISTANCE_BENCHMARK_THRESHOLD, (entries + extra)[-1]["frame"],
                                        n_boot, method, seed))
    df = pd.DataFrame(rows, columns=BENCH_COLS)
    return df.sort_values("AUPRC", ascending=False, kind="mergesort").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Training (Snakefile_training -> ENCODE_rE2G train_model.py)
# ---------------------------------------------------------------------------

def get_params(default: dict, override: "Optional[dict]") -> dict:
    """get_params.py: defaults overridden, 'null' -> None, numeric strings converted."""
    out = dict(default)
    if isinstance(override, dict):
        out.update(override)
    out = {k: (None if v == "null" else v) for k, v in out.items()}
    for k, v in list(out.items()):
        if isinstance(v, str):
            try:
                out[k] = int(float(v)) if k in ("max_iter", "random_state", "n_jobs") else float(v)
            except ValueError:
                pass
    return out


def train_models(df, ft, epsilon: float = EPSILON, params: "Optional[dict]" = None, polynomial: bool = False,
                 model_name: str = "E2G"):
    """train_model.train_and_predict: full model + per-held-out-chromosome models, coefficients and metrics."""
    pd, np = _pd(), _np()
    from sklearn.linear_model import LogisticRegression  # type: ignore
    from sklearn.metrics import log_loss, roc_auc_score  # type: ignore
    params = dict(params or DEFAULT_TRAIN_PARAMS)
    params.pop("verbose", None)
    feats = list(ft["feature"])
    df = df.reset_index(drop=True).copy()
    M = transform_features(df, feats, epsilon, False)
    names = list(feats)
    if polynomial:
        from sklearn.preprocessing import PolynomialFeatures  # type: ignore
        poly = PolynomialFeatures(degree=2)
        M = poly.fit_transform(M)
        names = list(poly.get_feature_names_out(feats))
    Y = df["Regulated"].map(_truthy).to_numpy().astype(np.int64)
    full = LogisticRegression(**params).fit(M, Y)
    pf = full.predict_proba(M)[:, 1]
    df[model_name + ".Score_full"] = pf
    coefs = [{"feature": n, "coefficient": float(c), "test_chr": "none"} for n, c in zip(names, full.coef_[0])]
    weights = {"full": [float(full.intercept_[0])] + [float(c) for c in full.coef_[0]], "cv": {}}
    metrics = []
    chroms = np.unique(df["chr"].astype(str))
    if len(chroms) > 1:
        df[model_name + ".Score"] = np.nan
        for c in chroms:
            te = np.flatnonzero(df["chr"].astype(str).to_numpy() == c)
            tr = np.setdiff1d(np.arange(len(Y)), te)
            if len(np.unique(Y[tr])) < 2:
                continue
            m = LogisticRegression(**params).fit(M[tr], Y[tr])
            pt = m.predict_proba(M[te])[:, 1]
            df.loc[te, model_name + ".Score"] = pt
            weights["cv"][c] = [float(m.intercept_[0])] + [float(x) for x in m.coef_[0]]
            coefs += [{"feature": n, "coefficient": float(x), "test_chr": c} for n, x in zip(names, m.coef_[0])]
            ytr, yte = Y[tr], Y[te]
            ptr = m.predict_proba(M[tr])[:, 1]
            has = yte.sum() > 0 and yte.sum() < len(yte)
            metrics.append({"test_chr": c,
                            "log_loss_test_full": log_loss(yte, pf[te], labels=[0, 1]) if yte.sum() > 0 else np.nan,
                            "log_loss_train": log_loss(ytr, ptr, labels=[0, 1]),
                            "log_loss_test": log_loss(yte, pt, labels=[0, 1]) if yte.sum() > 0 else np.nan,
                            "AUROC_test_full": roc_auc_score(yte, pf[te]) if has else np.nan,
                            "AUROC_train": roc_auc_score(ytr, ptr), "AUROC_test": roc_auc_score(yte, pt) if has else np.nan,
                            "AUPRC_test_full": statistic_aupr(yte, pf[te]) if yte.sum() > 0 else np.nan,
                            "AUPRC_train": statistic_aupr(ytr, ptr),
                            "AUPRC_test": statistic_aupr(yte, pt) if yte.sum() > 0 else np.nan,
                            "n_test_pos": int(yte.sum()), "n_test_neg": int(len(yte) - yte.sum()),
                            "n_train_neg": int(len(ytr) - ytr.sum()), "n_train_pos": int(ytr.sum())})
        ok = df[model_name + ".Score"].notna().to_numpy()
        sc = df[model_name + ".Score"].to_numpy(float)
        metrics.append({"test_chr": "all", "log_loss_test_full": log_loss(Y, pf), "log_loss_train": np.nan,
                        "log_loss_test": log_loss(Y[ok], sc[ok], labels=[0, 1]), "AUROC_test_full": roc_auc_score(Y, pf),
                        "AUROC_train": np.nan, "AUROC_test": roc_auc_score(Y[ok], sc[ok]),
                        "AUPRC_test_full": statistic_aupr(Y, pf), "AUPRC_train": np.nan,
                        "AUPRC_test": statistic_aupr(Y[ok], sc[ok]), "n_test_pos": int(Y.sum()),
                        "n_test_neg": int(len(Y) - Y.sum()), "n_train_pos": np.nan, "n_train_neg": np.nan})
    return {"predictions": df, "coefficients": pd.DataFrame(coefs), "metrics": pd.DataFrame(metrics),
            "weights": weights, "model_features": names, "models": full}


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def resource(override: "Optional[str]", name: str, required: bool = True) -> "Optional[Path]":
    if override:
        p = Path(override)
        if not p.exists():
            raise SystemExit(f"not found: {p}")
        return p
    for p in (RES_DIR / name, RES_DIR / "example" / name):
        if p.is_file():
            return p
    if required:
        raise SystemExit(f"{name} not found in {RES_DIR}; run `igvfagent sce2g-pipeline setup` or pass the path")
    return None


def _models_from_args(args) -> "List[dict]":
    md = Path(args.models_dir) if getattr(args, "models_dir", None) else MODELS_DIR
    return [resolve_model(m, md) for m in (args.models or ["multiome_powerlaw_v3", "scATAC_powerlaw_v3"])]


def _count_value(v) -> float:
    """A number, or a file holding one (fragment_count.txt / cell_count.txt / umi_count.txt); 0 if absent."""
    if v is None:
        return 0.0
    p = Path(str(v))
    if p.is_file():
        txt = p.read_text().strip().split()
        return float(txt[0]) if txt else 0.0
    return float(v)


def cmd_setup(args: argparse.Namespace) -> int:
    todo = [(MODELS_DIR / m / "qnorm_reference.tsv.gz", RAW + f"models/{m}/qnorm_reference.tsv.gz") for m in MODEL_SPECS]
    todo += [(RES_DIR / n, u) for n, (u, _) in RESOURCES.items()]
    if args.crispr:
        todo.append((RES_DIR / CRISPR_RESOURCE[0], CRISPR_RESOURCE[1]))
    if args.gtf:
        todo.append((RES_DIR / GTF_RESOURCE[0], GTF_RESOURCE[1]))
    if args.example:
        todo += [(RES_DIR / "example" / n, u) for n, u in EXAMPLE_FILES.items()]
    failed, stat = [], {"fetched": 0, "present": 0}
    for dest, url in todo:
        if args.dry_run:
            print(f"would fetch {url} -> {dest}")
            continue
        try:
            stat[download(url, dest, args.force)] += 1
        except Exception as e:  # noqa: BLE001
            failed.append(f"{dest.name}: {e}")
    if not args.dry_run:
        print(f"scE2G pipeline data in {DATA_ROOT}: fetched {stat['fetched']}, already present {stat['present']}, "
              f"failed {len(failed)} (pinned {UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}, {E2G_REPO}@{E2G_COMMIT[:7]})")
    for f in failed:
        print("  FAILED", f, file=sys.stderr)
    return 1 if failed else 0


def cmd_models(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, "models")
    rows = []
    for name in MODEL_SPECS:
        m = resolve_model(name, Path(args.models_dir) if args.models_dir else None)
        for f, w in zip(["(intercept)"] + m["features"], m["full"]):
            cvw = [v[(["(intercept)"] + m["features"]).index(f)] for v in m["cv"].values()]
            rows.append({"model": name, "feature": f, "weight_full": w, "weight_cv_min": min(cvw), "weight_cv_max": max(cvw),
                         "n_cv_models": len(cvw), "score_threshold": float(m["score_threshold"]),
                         "tpm_threshold": m["tpm_threshold"], "qnorm_reference": str(m["qnorm"] or "missing (run setup)")})
        if args.export:
            d = Path(args.export) / name
            d.mkdir(parents=True, exist_ok=True)
            write_tsv(m["feature_table"], d / "feature_table.tsv", quiet=True)
            (d / f"score_threshold_{m['score_threshold']}").write_text("")
            (d / f"tpm_threshold_{MODEL_SPECS[name]['tpm_threshold']}").write_text("")
            (d / "model_coefficients.json").write_text(json.dumps(
                {"model": name, "model_features": m["features"], "full": m["full"], "cv": m["cv"], "polynomial": False,
                 "epsilon": EPSILON, "source": f"{UPSTREAM_REPO}@{UPSTREAM_COMMIT}"}, indent=1))
            if m["qnorm"]:
                shutil.copy(m["qnorm"], d / "qnorm_reference.tsv.gz")
            print(f"Wrote: {d}")
    tab = pd.DataFrame(rows)
    write_tsv(tab, out / "model_weights.tsv")
    lines = [f"# scE2G v3 models ({UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]})", ""]
    for name, spec in MODEL_SPECS.items():
        sub = tab[tab["model"] == name]
        lines += [f"## {name} -- {spec['pretty']}", "",
                  f"score threshold {model_threshold_str(spec['score_threshold'])} on E2G.Score.qnorm; TPM threshold "
                  f"{spec['tpm_threshold']}", "",
                  md_table(["feature", "weight (full)", "cv min", "cv max"],
                           [[r.feature, f"{r.weight_full:.4f}", f"{r.weight_cv_min:.4f}", f"{r.weight_cv_max:.4f}"]
                            for r in sub.itertuples()]), ""]
    write_summary(out, {"models": {n: {"features": resolve_model(n)["features"], "score_threshold": float(s["score_threshold"]),
                                       "tpm_threshold": float(s["tpm_threshold"])} for n, s in MODEL_SPECS.items()}}, lines)
    return 0


def cmd_frag_to_tagalign(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "tagalign")
    sizes = resource(args.chrom_sizes, "GRCh38_EBV.no_alt.chrom.sizes.tsv")
    res = frag_to_tagalign(args.fragments, sizes, out / args.cluster, preprocessed=args.preprocessed)
    write_summary(out, res, ["# fragments -> tagAlign", "", md_table(["item", "value"], [[k, v] for k, v in res.items()])])
    return 0


def cmd_bigwig(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, "bigwig")
    sizes = read_chrom_sizes(resource(args.chrom_sizes, "GRCh38_EBV.no_alt.chrom.sizes.tsv"))
    sd = dict(sizes)
    frag = read_fragments(args.fragments)
    frag = frag[frag["chr"].isin(sd)]
    n = _count_value(args.fragment_count) if args.fragment_count else float(len(frag))
    d = out / args.cluster
    d.mkdir(parents=True, exist_ok=True)
    made = {}
    for tag, scale in (("ATAC", 1.0), ("ATAC_norm", 1e6 / n if n else 1.0)):
        if tag == "ATAC" and not args.raw:
            continue
        bg = fragment_coverage(frag, sd, scale)
        bg = sort_by_chrom_order(bg, [c for c, _ in sizes])
        bgp = d / f"{tag}.bg"
        bg.to_csv(bgp, sep="\t", header=False, index=False)
        bw = d / f"{tag}.bw"
        r = write_bigwig(bg, sizes, bw)
        if r == Path("__binary__"):
            szp = d / "chrom.sizes"
            pd.DataFrame(sizes).to_csv(szp, sep="\t", header=False, index=False)
            subprocess.check_call([shutil.which("bedGraphToBigWig"), str(bgp), str(szp), str(bw)])
            r = bw
        if r is None:
            print(f"Wrote: {bgp}")
            print(f"bigWig skipped (install pyBigWig or put bedGraphToBigWig on PATH); would run: "
                  f"bedGraphToBigWig {bgp} <chrom.sizes> {bw}")
            made[tag] = str(bgp)
        else:
            print(f"Wrote: {r}")
            made[tag] = str(r)
    write_summary(out, {"fragment_count": n, "tracks": made}, ["# ATAC coverage tracks", "",
                  md_table(["track", "path"], [[k, v] for k, v in made.items()])])
    return 0


def cmd_kendall_pairs(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "kendall_pairs")
    pairs = make_kendall_pairs(read_narrowpeak(args.narrowpeak), read_tsv(args.abc_predictions))
    p = write_tsv(pairs, out / args.cluster / "Kendall" / "Pairs.tsv.gz")
    s = {"pairs": len(pairs), "peaks": int(pairs["PeakName"].nunique()), "genes": int(pairs["TargetGene"].nunique()), "file": p}
    write_summary(out, s, ["# Kendall candidate pairs", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def _cells_for(args, rna_cells, cluster_dir: "Optional[Path]" = None):
    bc = None
    if getattr(args, "rna_unfiltered", False):
        f = getattr(args, "cell_barcodes", None)
        if f:
            bc = [l for l in Path(f).read_text().split("\n") if l]
        else:
            bc = list(dict.fromkeys(read_fragments(args.fragments)["barcode"].astype(str)))
    return choose_cells(rna_cells, bc, args.max_cell_count)


def _rna_cells(path) -> "List[str]":
    p = Path(path)
    if p.suffix in (".h5ad", ".h5"):
        import anndata  # type: ignore
        return [str(c) for c in anndata.read_h5ad(str(p), backed="r").obs_names]
    if p.is_dir():
        return read_rna_matrix(p)[2]
    with _open_text(p) as fh:
        return [c.strip().strip('"') for c in fh.readline().rstrip("\n").split(",")[1:]]


def cmd_atac_matrix(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "atac_matrix")
    pairs = read_tsv(args.pairs)
    peaks = pairs.drop_duplicates("PeakName")[["chr", "start", "end", "PeakName"]].reset_index(drop=True)
    cells = _cells_for(args, _rna_cells(args.rna_matrix))
    M = build_atac_matrix(peaks, args.fragments, cells, args.count_mode)
    files = save_matrix(M, list(peaks["PeakName"]), cells, out / args.cluster / "Kendall" / "atac_matrix")
    s = {"peaks": M.shape[0], "cells": M.shape[1], "nonzero": int(M.nnz), "count_mode": args.count_mode, **files}
    write_summary(out, s, ["# ATAC peak x cell matrix", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def _gene_key_from_args(args) -> "Optional[Dict[str, str]]":
    if getattr(args, "no_gene_mapping", False):
        return None
    gtf = resource(args.gtf, GTF_RESOURCE[0], required=False)
    if gtf is None:
        raise SystemExit("the GENCODE GTF is needed to map RNA gene names (setup --gtf, --gtf PATH), "
                         "or pass --no-gene-mapping if the RNA matrix already uses TSS500 gene names")
    tss = read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"))
    return build_gene_key(parse_gtf_genes(gtf), tss)


def run_kendall_stage(args, pairs, cluster_dir: Path) -> dict:
    """atac matrix (from fragments unless given) + compute_kendall, writing upstream file names."""
    kd = cluster_dir / "Kendall"
    kd.mkdir(parents=True, exist_ok=True)
    if getattr(args, "atac_matrix", None):
        M, peak_names, cells = load_matrix(args.atac_matrix)
    else:
        peaks = pairs.drop_duplicates("PeakName")[["chr", "start", "end", "PeakName"]].reset_index(drop=True)
        cells = _cells_for(args, _rna_cells(args.rna_matrix), cluster_dir)
        M = build_atac_matrix(peaks, args.fragments, cells, getattr(args, "count_mode", "fragment"))
        peak_names = list(peaks["PeakName"])
        save_matrix(M, peak_names, cells, kd / "atac_matrix")
    key = _gene_key_from_args(args)
    kend, gex, n_umi, n_cells = compute_kendall(pairs, M, peak_names, cells, args.rna_matrix, key)
    write_tsv(kend, kd / "Pairs.Kendall.tsv.gz")
    write_tsv(gex, kd / "gene_expression_metrics.tsv.gz")
    write_text(cluster_dir / "umi_count.txt", f"{int(n_umi) if float(n_umi).is_integer() else n_umi}\n", quiet=True)
    write_text(cluster_dir / "cell_count.txt", f"{n_cells}\n", quiet=True)
    k = kend["Kendall"]
    return {"pairs_scored": int(len(kend)), "kendall_defined": int(k.notna().sum()), "kendall_mean": float(k.mean()) if k.notna().any() else None,
            "genes_mapped": int(len(gex)), "umi_count": float(n_umi), "cell_count": int(n_cells), "peaks": len(peak_names),
            "kendall_file": str(kd / "Pairs.Kendall.tsv.gz"), "gene_expression_file": str(kd / "gene_expression_metrics.tsv.gz")}


def cmd_kendall(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "kendall")
    if not args.atac_matrix and not args.fragments:
        raise SystemExit("pass --atac-matrix (from atac-matrix) or --fragments")
    s = run_kendall_stage(args, read_tsv(args.pairs), out / args.cluster)
    write_summary(out, s, ["# Kendall correlation (scE2G compute_kendall)", "",
                           md_table(["item", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in s.items()])])
    return 0


def cmd_arc(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "arc")
    col = "ABC.Score" if args.hic else args.abc_score_col
    arc, fit = compute_arc(read_tsv(args.abc_predictions), read_tsv(args.kendall), col)
    p = write_tsv(arc, out / args.cluster / "ARC" / "EnhancerPredictionsAllPutative_ARC.tsv.gz")
    s = {"pairs": len(arc), "abc_score_col": col, **fit, "file": p}
    write_summary(out, s, ["# ARC-E2G", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def activity_stage(abc, enh_list, ft, tss, sizes: "Optional[Dict[str, int]]", gene_classes, cluster_dir: Path) -> "Tuple[Any, dict]":
    nf = cluster_dir / "new_features"
    nf.mkdir(parents=True, exist_ok=True)
    inputs = input_features_of(ft)
    numcand = numtss = num5 = sum5 = None
    if "numCandidateEnhGene" in inputs:
        numcand = num_candidate_enh_gene(abc)
        write_tsv(numcand, nf / "NumCandidateEnhGene.tsv", quiet=True)
    if "numTSSEnhGene" in inputs:
        numtss = num_tss_enh_gene(abc, tss)
        write_tsv(numtss, nf / "NumTSSEnhGene.tsv", quiet=True)
    if "numNearbyEnhancers" in inputs or "sumNearbyEnhancers" in inputs:
        if enh_list is None:
            raise SystemExit("numNearbyEnhancers needs the ABC EnhancerList (--enhancer-list)")
        num5, sum5 = nearby_enhancers(abc, enh_list, sizes)
        num5.to_csv(nf / "NumEnhancersEG5kb.txt", sep="\t", header=False, index=False)
        sum5.to_csv(nf / "SumEnhancersEG5kb.txt", sep="\t", header=False, index=False)
    act = activity_only_features(abc, ft, numcand, numtss, num5, sum5, gene_classes)
    write_tsv(act, cluster_dir / "ActivityOnly_features.tsv.gz")
    return act, {"pairs": len(act), "columns": list(act.columns)}


def cmd_activity_features(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "activity_features")
    models = _models_from_args(args)
    ft = biosample_feature_table(models)
    gc_path = resource(args.gene_classes, "gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv",
                       required=False)
    sizes_p = resource(args.chrom_sizes, "GRCh38_EBV.no_alt.chrom.sizes.tsv", required=False)
    act, s = activity_stage(read_tsv(args.abc_predictions), read_tsv(args.enhancer_list) if args.enhancer_list else None, ft,
                            read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed")),
                            dict(read_chrom_sizes(sizes_p)) if sizes_p else None, read_tsv(gc_path) if gc_path else None,
                            out / args.cluster)
    write_summary(out, s, ["# ENCODE_rE2G activity-only features", "", f"{s['pairs']:,} pairs", "",
                           "columns: " + ", ".join(s["columns"])])
    return 0


def features_stage(act, ft, arc_path: "Optional[str]", cluster_dir: Path, arc_frame=None, user_efc=None):
    write_tsv(ft, cluster_dir / "feature_table.tsv", quiet=True)
    efc = external_features_config(arc_path if needs_arc(ft) else None, user_efc)
    write_tsv(efc, cluster_dir / "external_features_config.tsv", quiet=True)
    sources = {str(arc_path): arc_frame} if arc_frame is not None and arc_path else None
    plus = merge_external_features(act, ft, efc, sources)
    write_tsv(plus, cluster_dir / "ActivityOnly_plus_external_features.tsv.gz", quiet=True)
    gw = final_features(plus, ft)
    write_tsv(gw, cluster_dir / "genomewide_features.tsv.gz")
    return gw


def cmd_features(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "features")
    ft = biosample_feature_table(_models_from_args(args))
    if needs_arc(ft) and not args.arc:
        raise SystemExit("these models need ARC-E2G features: pass --arc EnhancerPredictionsAllPutative_ARC.tsv.gz")
    gw = features_stage(read_tsv(args.activity_features), ft, args.arc, out / args.cluster,
                        user_efc=read_tsv(args.external_config) if args.external_config else None)
    s = {"pairs": len(gw), "features": list(ft["feature"].unique())}
    write_summary(out, s, ["# genome-wide feature table", "", f"{len(gw):,} pairs; features: " + ", ".join(s["features"])])
    return 0


def predict_stage(gw, models: "Sequence[dict]", cluster_dir: Path, crispr_benchmarking: bool = False, qnorm: bool = True,
                  bedpe_dir: "Optional[Path]" = None, include_self_promoter: bool = True) -> "List[dict]":
    res = []
    for m in models:
        pred, info = apply_model(gw, m, crispr_benchmarking, qnorm=qnorm)
        md = cluster_dir / m["name"]
        write_tsv(pred, md / "scE2G_predictions.tsv.gz")
        thr = m["score_threshold"]
        r = {"model": m["name"], "pairs": len(pred), "predictions": str(md / "scE2G_predictions.tsv.gz"), **info}
        if thr is not None:
            t = threshold_predictions(pred, float(thr), FINAL_SCORE_COL, include_self_promoter)
            tp = md / f"scE2G_predictions_threshold{model_threshold_str(thr)}.tsv.gz"
            write_tsv(t, tp)
            r.update({"threshold": float(thr), "links": len(t), "thresholded": str(tp)})
            if bedpe_dir is not None:
                bp = Path(bedpe_dir) / m["name"] / f"scE2G_predictions_threshold{model_threshold_str(thr)}.bedpe"
                bp.parent.mkdir(parents=True, exist_ok=True)
                predictions_bedpe(t).to_csv(bp, sep="\t", header=False, index=False)
                print(f"Wrote: {bp}")
                r["bedpe"] = str(bp)
        r["mean_score_qnorm"] = float(pred[FINAL_SCORE_COL].mean()) if len(pred) else None
        res.append(r)
    return res


def cmd_predict(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "predict")
    gw = read_tsv(args.features)
    res = predict_stage(gw, _models_from_args(args), out / args.cluster, args.crispr_benchmarking, not args.no_qnorm,
                        (out / args.cluster) if args.bedpe else None, not args.exclude_self_promoter)
    write_summary(out, {"cluster": args.cluster, "models": res},
                  ["# scE2G predictions", "", md_table(["model", "pairs", "threshold", "links"],
                                                       [[r["model"], r["pairs"], r.get("threshold"), r.get("links")] for r in res])])
    return 0


def cmd_qnorm_ref(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "qnorm_ref")
    m = resolve_model(args.model, Path(args.models_dir) if args.models_dir else None)
    gw = read_tsv(args.features).replace([_np().inf, -_np().inf], _np().nan).fillna(0)
    M = transform_features(gw, m["features"], EPSILON, m.get("polynomial", False))
    gw[SCORE_BASE] = logistic_score(M, m["full"])
    ref = qnorm_reference_from_scores(gw[SCORE_BASE], args.n_scores)
    p = write_tsv(ref, out / m["name"] / "qnorm_reference.tsv.gz")
    if args.write_predictions:
        cvm = m["cv"]
        cv = gw[SCORE_BASE].to_numpy().copy()
        for c in gw["chr"].astype(str).unique():
            if c in cvm:
                idx = (gw["chr"].astype(str) == c).to_numpy()
                cv[idx] = logistic_score(M[idx], cvm[c])
        gw[SCORE_BASE + ".cv"] = cv
        gw = filter_by_tpm_raw(gw, float(args.tpm_threshold if args.tpm_threshold is not None else m["tpm_threshold"]))
        write_tsv(gw, out / m["name"] / "genomewide_predictions.tsv.gz")
    s = {"model": m["name"], "n_scores": args.n_scores, "pairs": len(gw), "file": p}
    write_summary(out, s, ["# qnorm reference", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def filter_by_tpm_raw(df, tpm_threshold: float):
    """run_e2g_cv_qnormref.filter_by_tpm: the same filter on the raw (not qnorm) scores, cv included."""
    np = _np()
    if "RNA_pseudobulkTPM" not in df.columns or tpm_threshold == 0:
        return df
    low = ~(pd_numeric(df["RNA_pseudobulkTPM"]) >= tpm_threshold)
    for c in (SCORE_BASE, SCORE_BASE + ".cv"):
        df[c + ".ignoreTPM"] = df[c]
        df[c] = np.where(low, 0.0, df[c + ".ignoreTPM"])
    return df


def cmd_gene_lists(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "gene_lists")
    pred = read_tsv(args.predictions)
    tpm = args.tpm_threshold if args.tpm_threshold is not None else (resolve_model(args.model)["tpm_threshold"] if args.model else 0.0)
    gex = read_tsv(args.gene_expression) if args.gene_expression else None
    gl, el = element_gene_lists(read_tsv(args.gene_list), read_tsv(args.enhancer_list), pred, args.cluster, tpm, gex)
    d = out / args.cluster / (args.model or "model")
    write_tsv(gl, d / "scE2G_gene_list.tsv.gz")
    write_tsv(el, d / "scE2G_element_list.tsv.gz")
    s = {"genes": len(gl), "genes_considered": int(gl["considered_for_predictions"].sum()), "elements": len(el)}
    write_summary(out, s, ["# gene and element lists", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def stats_row(pred_full, pred_thr, cluster: str, model_name: str, fragments, cells, umis) -> dict:
    r = stats_from_predictions(pred_full, pred_thr)
    r.update({"cluster": cluster, "model_name": model_name, "fragments_total": _count_value(fragments),
              "cell_count": _count_value(cells), "umi_count": _count_value(umis)})
    return {k: r[k] for k in STATS_COLS}


def cmd_stats(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, "stats")
    full = read_tsv(args.predictions)
    if args.thresholded:
        thr_df = read_tsv(args.thresholded)
        thr = args.threshold
    else:
        m = resolve_model(args.model_name) if args.threshold is None else None
        thr = float(m["score_threshold"]) if m else args.threshold
        thr_df = threshold_predictions(full, thr)
    row = stats_row(full, thr_df, args.cluster, args.model_name, args.fragment_count, args.cell_count, args.umi_count)
    tag = model_threshold_str(str(thr)) if thr is not None else "NA"
    p = write_tsv(pd.DataFrame([row]), out / args.cluster / args.model_name / f"scE2G_predictions_threshold{tag}_stats.tsv")
    write_summary(out, {"stats": row, "file": p}, ["# prediction statistics", "",
                  md_table(["metric", "value"], [[k, v] for k, v in row.items()])])
    return 0


def qc_stage(stats, reference, out: Path, plots: bool = True) -> dict:
    qd = out / "qc_plots"
    qd.mkdir(parents=True, exist_ok=True)
    write_tsv(stats, qd / "all_qc_stats.tsv")
    ds, warns, comp = qc_summary(stats, reference)
    write_tsv(ds, qd / "dataset_metrics.tsv", quiet=True)
    if len(comp):
        write_tsv(comp, qd / "qc_reference_comparison.tsv")
    write_text(qd / "qc_warnings.txt", "\n".join(warns) + ("\n" if warns else ""), quiet=True)
    figs = qc_figures(stats, ds, reference, qd) if plots else []
    return {"warnings": warns, "comparison": comp, "dataset": ds, "figures": [str(f) for f in figs]}


def cmd_qc(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, "qc")
    stats = pd.concat([read_tsv(p) for p in args.stats], ignore_index=True)
    ref_p = resource(args.reference, "reference_qc_metrics_sheth_qiu_2024.tsv", required=False)
    q = qc_stage(stats, read_tsv(ref_p) if ref_p else None, out, not args.no_plots)
    lines = ["# scE2G prediction QC", ""] + [f"- {w}" for w in q["warnings"]] + ["", md_table(
        STATS_COLS[:11], stats[STATS_COLS[:11]].itertuples(index=False))]
    if len(q["comparison"]):
        c = q["comparison"]
        lines += ["", "## vs Sheth, Qiu 2024 reference clusters", "",
                  md_table(["cluster", "model", "metric", "value", "ref median", "ref percentile"],
                           [[r.cluster, r.model_name, r.metric, _fmt(r.value), _fmt(r.reference_median), _fmt(r.reference_percentile, 1)]
                            for r in c.itertuples()])]
    write_summary(out, {"warnings": q["warnings"], "n_clusters": int(stats["cluster"].nunique()), "figures": q["figures"]}, lines)
    return 0


def cmd_crispr_features(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "crispr_features")
    pred = read_tsv(args.predictions)
    ft = _read_feature_table(args.feature_table) if args.feature_table else biosample_feature_table(_models_from_args(args))
    crispr = read_tsv(resource(args.crispr, CRISPR_RESOURCE[0]))
    tss = read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"))
    tag = "NAfilled" if not args.no_nafill else "withNA"
    merged, missing = crispr_features(pred, crispr, ft, tss, not args.no_nafill, args.mode == "apply")
    d = out / args.cluster / (args.model_name or "features")
    p = write_tsv(merged, d / f"EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{tag}.tsv.gz")
    if args.mode == "training":
        write_tsv(missing, d / f"missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{tag}.tsv.gz")
        write_tsv(process_crispr_for_training(merged, tss),
                  d / f"for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{tag}.tsv.gz")
    s = {"crispr_pairs": len(merged), "without_prediction": len(missing), "regulated": int(merged["Regulated"].map(_truthy).sum()),
         "file": p}
    write_summary(out, s, ["# CRISPR features", "", md_table(["item", "value"], [[k, v] for k, v in s.items()])])
    return 0


def _bench_entries(paths, model_names, thresholds, clusters):
    entries = []
    for i, p in enumerate(paths):
        pp = Path(p)
        model = model_names[i] if model_names else pp.parent.name
        cluster = clusters[i] if clusters else pp.parent.parent.name
        if thresholds:
            t = float(thresholds[i])
        elif model in MODEL_SPECS:
            t = float(MODEL_SPECS[model]["score_threshold"])
        else:
            raise SystemExit(f"no threshold for model {model!r}: pass --model-thresholds")
        entries.append({"cluster": cluster, "model_name": model, "model_threshold": t, "frame": read_tsv(pp)})
    return entries


def cmd_benchmark(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "benchmark")
    df = benchmark_table(_bench_entries(args.crispr_features, args.model_names, args.model_thresholds, args.clusters),
                         args.n_boot, args.ci_method, args.seed)
    write_tsv(df, out / "crispr_benchmarking_performance_summary.tsv")
    plt = _plt()
    if plt is not None and not args.no_plots and len(df):
        fig, ax = plt.subplots(figsize=(6, 0.4 * len(df) + 1))
        lab = (df["cluster"].astype(str) + " / " + df["model"].astype(str)).tolist()[::-1]
        v = df["AUPRC"].to_numpy()[::-1]
        lo, hi = df["AUPRC_95CI_low"].to_numpy()[::-1], df["AUPRC_95CI_high"].to_numpy()[::-1]
        ax.barh(range(len(v)), v, color="#2a78d6")
        ax.errorbar(v, range(len(v)), xerr=[_np().clip(v - lo, 0, None), _np().clip(hi - v, 0, None)], fmt="none", ecolor=INK, lw=0.8)
        ax.set_yticks(range(len(v)))
        ax.set_yticklabels(lab, fontsize=7)
        _style(ax, "AUPRC on CRISPR element-gene pairs (95% bootstrap CI)", "AUPRC", "")
        _save(fig, out / "auprc.png")
    lines = ["# CRISPR benchmark (scE2G benchmark_performance)", "",
             md_table(["cluster", "model", "AUPRC", "95% CI", "precision@70% recall", "precision@threshold", "recall@threshold"],
                      [[r.cluster, r.model, _fmt(r.AUPRC), f"{_fmt(r.AUPRC_95CI_low)}-{_fmt(r.AUPRC_95CI_high)}",
                        _fmt(r.precision_70_pct_recall), _fmt(r.precision_model_threshold), _fmt(r.recall_model_threshold)]
                       for r in df.itertuples()])]
    write_summary(out, {"rows": df.to_dict(orient="records")}, lines)
    return 0


def cmd_train(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, "train")
    if args.feature_table:
        ft = _read_feature_table(args.feature_table)
    else:
        ft = embedded_feature_table(args.model)
    if args.crispr_features:
        df = read_tsv(args.crispr_features)
        if "chrom" in df.columns:
            df = process_crispr_for_training(df, read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed")))
    else:
        if not args.features:
            raise SystemExit("pass --crispr-features (NAfilled CRISPR x features table) or --features genomewide_features.tsv.gz")
        tss = read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"))
        merged, _ = crispr_features(read_tsv(args.features), read_tsv(resource(args.crispr, CRISPR_RESOURCE[0])), ft, tss,
                                    True, apply_mode=False)
        df = process_crispr_for_training(merged, tss)
    params = get_params(DEFAULT_TRAIN_PARAMS, json.loads(args.override_params) if args.override_params else None)
    res = train_models(df, ft, EPSILON, params, args.polynomial)
    d = out / args.name
    (d / "model").mkdir(parents=True, exist_ok=True)
    write_tsv(ft, d / "feature_table.tsv")
    write_tsv(res["predictions"], d / "model" / "training_predictions.tsv")
    write_tsv(res["coefficients"], d / "model" / "model_coefficients.tsv")
    write_tsv(res["metrics"], d / "model" / "performance_metrics.tsv")
    (d / "model" / "training_params.json").write_text(json.dumps(_jsonable(params), indent=1))
    blob = {"model": args.name, "model_features": res["model_features"], "full": res["weights"]["full"], "cv": res["weights"]["cv"],
            "polynomial": bool(args.polynomial), "epsilon": EPSILON, "trained_on": str(args.crispr_features or args.features)}
    (d / "model_coefficients.json").write_text(json.dumps(blob, indent=1))
    print(f"JSON: {d / 'model_coefficients.json'}")
    sc = "E2G.Score" if "E2G.Score" in res["predictions"].columns else "E2G.Score_full"
    y = res["predictions"]["Regulated"].map(_truthy).astype(int).to_numpy()
    s = res["predictions"][sc].fillna(0).to_numpy()
    thr = args.score_threshold if args.score_threshold is not None else threshold_at_target_recall(y, s, 0.7)
    if thr is not None:
        (d / f"score_threshold_{thr:.3f}".replace("_0.", "_.")).write_text("")
    (d / f"tpm_threshold_{args.tpm_threshold:g}").write_text("")
    if args.pickle:
        import pickle
        with open(d / "model.pkl", "wb") as fh:
            pickle.dump(res["models"], fh)
    met = res["metrics"]
    allrow = met[met["test_chr"] == "all"] if len(met) else met
    lines = ["# scE2G model training (ENCODE_rE2G train_model)", "", f"{len(df):,} CRISPR pairs, "
             f"{int(y.sum())} regulated; features: {', '.join(res['model_features'])}", "",
             md_table(["feature", "coefficient (full)"], [[r.feature, f"{r.coefficient:.4f}"]
                                                          for r in res["coefficients"][res["coefficients"]["test_chr"] == "none"].itertuples()])]
    if len(allrow):
        r = allrow.iloc[0]
        lines += ["", f"held-out-chromosome AUPRC {_fmt(r['AUPRC_test'])} (full-model in-sample {_fmt(r['AUPRC_test_full'])})"]
    lines += ["", f"model directory: `{d}` (use with `predict --models {d}`)"]
    write_summary(out, {"model_dir": d, "n_pairs": len(df), "n_regulated": int(y.sum()), "weights": res["weights"],
                        "score_threshold": thr, "metrics_all": allrow.to_dict(orient="records")}, lines)
    return 0


# ---------------------------------------------------------------------------
# run: the per-cluster pipeline (Snakefile rule order)
# ---------------------------------------------------------------------------

ABC_LAYOUT = {"narrowpeak": "Peaks/macs2_peaks.narrowPeak.sorted",
              "abc_predictions": "Predictions/EnhancerPredictionsAllPutative.tsv.gz",
              "enhancer_list": "Neighborhoods/EnhancerList.txt", "gene_list": "Neighborhoods/GeneList.txt"}


def _cluster_rows(args) -> "List[dict]":
    if args.cluster_config:
        cfg = read_tsv(args.cluster_config, dtype=str, keep_default_na=False)
        rows = []
        for r in cfg.to_dict(orient="records"):
            rows.append({"cluster": r["cluster"], "fragments": r["atac_frag_file"], "rna_matrix": r.get("rna_matrix_file") or None,
                         "hic": bool(r.get("HiC_type")), "abc_dir": r.get("abc_dir") or None,
                         "models": [m for m in (r.get("model_dir") or "").split(",") if m] or args.models})
        return rows
    if not args.fragments:
        raise SystemExit("pass --fragments (and --rna-matrix, --abc-dir) or --cluster-config")
    return [{"cluster": args.cluster, "fragments": args.fragments, "rna_matrix": args.rna_matrix, "hic": args.hic,
             "abc_dir": args.abc_dir, "models": args.models}]


def _abc_paths(args, row) -> "Dict[str, Optional[Path]]":
    out = {}
    for k, rel in ABC_LAYOUT.items():
        ov = getattr(args, k, None) if not args.cluster_config else None
        if ov:
            out[k] = Path(ov)
        elif row.get("abc_dir"):
            p = Path(row["abc_dir"]) / rel
            out[k] = p if p.exists() else None
        else:
            out[k] = None
    return out


def abc_command_hint(cluster_rows, out: Path) -> str:
    """make_biosample_config(): the ABC biosample table scE2G generates, and the command that consumes it."""
    pd = _pd()
    rows = [{"biosample": r["cluster"], "DHS": "", "ATAC": str(out / r["cluster"] / "tagAlign" / "tagAlign.sort.gz"),
             "H3K27ac": "", "default_accessibility_feature": "ATAC", "HiC_file": "", "HiC_type": "", "HiC_resolution": ""}
            for r in cluster_rows]
    p = out / "tmp" / "config_abc_biosamples.tsv"
    write_tsv(pd.DataFrame(rows), p, quiet=True)
    return (f"snakemake -s ABC-Enhancer-Gene-Prediction/workflow/Snakefile --config biosamplesTable={p} "
            f"results_dir={out} -j 8   # then rerun with --abc-dir {out}/<cluster>")


def run_cluster(args, row: dict, out: Path, refs: dict, models_dir: Path) -> dict:
    pd = _pd()
    cluster = row["cluster"]
    cd = out / cluster
    models = [resolve_model(m, models_dir) for m in (row["models"] or ["multiome_powerlaw_v3", "scATAC_powerlaw_v3"])]
    ft = biosample_feature_table(models)
    res: Dict[str, Any] = {"cluster": cluster, "models": [m["name"] for m in models]}
    tag = frag_to_tagalign(row["fragments"], refs["chrom_sizes"], cd, preprocessed=args.fragments_preprocessed)
    res["fragments"] = {k: tag[k] for k in ("fragments_raw", "fragment_count", "n_barcodes", "n_tags")}
    if args.igv_tracks:
        frag = read_fragments(row["fragments"])
        sd = refs["sizes"]
        bg = sort_by_chrom_order(fragment_coverage(frag[frag["chr"].isin(sd)], sd, 1e6 / max(tag["fragment_count"], 1)),
                                 list(sd))
        bgp = cd / "ATAC_norm.bg"
        bg.to_csv(bgp, sep="\t", header=False, index=False)
        bw = write_bigwig(bg, list(sd.items()), cd / "ATAC_norm.bw")
        print(f"Wrote: {bw if bw not in (None, Path('__binary__')) else bgp}")
    abc = _abc_paths(args, row)
    if abc["abc_predictions"] is None or abc["narrowpeak"] is None:
        res["status"] = "needs ABC outputs"
        return res
    allp = read_tsv(abc["abc_predictions"])
    arc_path, arc_frame = None, None
    if needs_arc(ft):
        if not row.get("rna_matrix"):
            raise SystemExit(f"{cluster}: models {res['models']} need Kendall/ARC features, which need --rna-matrix")
        pairs = make_kendall_pairs(read_narrowpeak(abc["narrowpeak"]), allp)
        write_tsv(pairs, cd / "Kendall" / "Pairs.tsv.gz")
        kargs = argparse.Namespace(atac_matrix=None, fragments=row["fragments"], rna_matrix=row["rna_matrix"],
                                   rna_unfiltered=args.rna_unfiltered, cell_barcodes=tag["cell_barcodes_file"],
                                   max_cell_count=args.max_cell_count, count_mode=args.count_mode,
                                   no_gene_mapping=args.no_gene_mapping, gtf=args.gtf, tss=str(refs["tss_path"]))
        res["kendall"] = run_kendall_stage(kargs, pairs, cd)
        kend = read_tsv(cd / "Kendall" / "Pairs.Kendall.tsv.gz")
        arc_frame, fit = compute_arc(allp, kend, "ABC.Score" if row.get("hic") else "powerlaw.Score")
        arc_path = str(cd / "ARC" / "EnhancerPredictionsAllPutative_ARC.tsv.gz")
        write_tsv(arc_frame, arc_path)
        res["arc"] = fit
    enh = read_tsv(abc["enhancer_list"]) if abc["enhancer_list"] else None
    act, _ = activity_stage(allp, enh, ft, refs["tss"], refs["sizes"], refs["gene_classes"], cd)
    gw = features_stage(act, ft, arc_path, cd, arc_frame)
    preds = predict_stage(gw, models, cd, bool(refs.get("crispr") is not None), not args.no_qnorm,
                          cd if args.igv_tracks else None)
    res["predictions"] = preds
    gex_path = cd / "Kendall" / "gene_expression_metrics.tsv.gz"
    gex = read_tsv(gex_path) if gex_path.exists() else None
    stats_rows = []
    for m, p in zip(models, preds):
        full = read_tsv(p["predictions"])
        thr = read_tsv(p["thresholded"]) if p.get("thresholded") else full.iloc[0:0]
        if abc["gene_list"] is not None and enh is not None:
            gl, el = element_gene_lists(read_tsv(abc["gene_list"]), enh, full, cluster, m["tpm_threshold"], gex)
            write_tsv(gl, cd / m["name"] / "scE2G_gene_list.tsv.gz", quiet=True)
            write_tsv(el, cd / m["name"] / "scE2G_element_list.tsv.gz", quiet=True)
        has_rna = gex is not None
        srow = stats_row(full, thr, cluster, m["name"], cd / "fragment_count.txt",
                         cd / "cell_count.txt" if has_rna else None, cd / "umi_count.txt" if has_rna else None)
        tagt = model_threshold_str(m["score_threshold"]) if m["score_threshold"] else "NA"
        write_tsv(pd.DataFrame([srow]), cd / m["name"] / f"scE2G_predictions_threshold{tagt}_stats.tsv", quiet=True)
        stats_rows.append(srow)
        if refs.get("crispr") is not None:
            merged, _ = crispr_features(full, refs["crispr"], ft, refs["tss"], True, True)
            write_tsv(merged, cd / m["name"] / "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz", quiet=True)
    res["stats"] = stats_rows
    res["status"] = "done"
    return res


def cmd_run(args: argparse.Namespace) -> int:
    pd = _pd()
    out = out_dir_for(args, args.cluster or "run")
    rows = _cluster_rows(args)
    sizes_p = resource(args.chrom_sizes, "GRCh38_EBV.no_alt.chrom.sizes.tsv")
    tss_p = resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed")
    gc_p = resource(args.gene_classes, "gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv",
                    required=False)
    refs = {"chrom_sizes": sizes_p, "sizes": dict(read_chrom_sizes(sizes_p)), "tss_path": tss_p, "tss": read_tss500(tss_p),
            "gene_classes": read_tsv(gc_p) if gc_p else None,
            "crispr": read_tsv(resource(args.crispr, CRISPR_RESOURCE[0])) if args.crispr else None}
    models_dir = Path(args.models_dir) if args.models_dir else MODELS_DIR
    results = [run_cluster(args, r, out, refs, models_dir) for r in rows]
    pending = [r for r in results if r["status"] != "done"]
    lines = [f"# scE2G pipeline run ({UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}, ported)", ""]
    if pending:
        hint = abc_command_hint(rows, out)
        print("ABC outputs (Peaks/, Predictions/, Neighborhoods/) are required for the Kendall/feature stages and are not "
              "produced by this port. Command that would produce them:\n  " + hint)
        lines += ["ABC outputs missing for: " + ", ".join(r["cluster"] for r in pending), "", "```", hint, "```", ""]
    done = [r for r in results if r["status"] == "done"]
    summary: Dict[str, Any] = {"clusters": results}
    if done:
        stats = pd.DataFrame([s for r in done for s in r["stats"]], columns=STATS_COLS)
        ref_p = resource(args.qc_reference, "reference_qc_metrics_sheth_qiu_2024.tsv", required=False)
        q = qc_stage(stats, read_tsv(ref_p) if ref_p else None, out, not args.no_plots)
        summary["qc_warnings"] = q["warnings"]
        lines += ["## Predictions", "", md_table(["cluster", "model", "pairs", "threshold", "links"],
                  [[r["cluster"], p["model"], p["pairs"], p.get("threshold"), p.get("links")] for r in done for p in r["predictions"]]),
                  "", "## QC", ""] + [f"- {w}" for w in q["warnings"]] + ["", md_table(STATS_COLS[:11],
                                                                                     stats[STATS_COLS[:11]].itertuples(index=False))]
        for r in done:
            if "kendall" in r:
                k = r["kendall"]
                lines += ["", f"{r['cluster']}: Kendall on {k['pairs_scored']:,} peak-gene pairs ({k['kendall_defined']:,} defined), "
                          f"{k['cell_count']} cells, {int(k['umi_count']):,} UMIs; ARC r = {_fmt(r['arc'].get('r'))}"]
        if refs.get("crispr") is not None:
            entries = [{"cluster": r["cluster"], "model_name": p["model"], "model_threshold": p.get("threshold") or 0.0,
                        "frame": read_tsv(out / r["cluster"] / p["model"] /
                                          "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz")}
                       for r in done for p in r["predictions"]]
            bt = benchmark_table(entries, args.n_boot, args.ci_method)
            write_tsv(bt, out / "crispr_benchmarking_performance_summary.tsv")
            summary["benchmark"] = bt.to_dict(orient="records")
            lines += ["", "## CRISPR benchmark", "", md_table(["cluster", "model", "AUPRC"],
                                                              [[r.cluster, r.model, _fmt(r.AUPRC)] for r in bt.itertuples()])]
    write_summary(out, summary, lines)
    return 0


# ---------------------------------------------------------------------------
# validate-example: the upstream chr22 test fixture
# ---------------------------------------------------------------------------

def validate_example(example_dir: Path, gtf: Path, tss_path: Path, count_mode: str = "fragment") -> dict:
    pd, np = _pd(), _np()
    exp = read_tsv(example_dir / "Pairs.Kendall.expected.tsv.gz")
    pairs = exp[["chr", "start", "end", "TargetGene", "PeakName", "PairName"]]
    peaks = pairs.drop_duplicates("PeakName")[["chr", "start", "end", "PeakName"]].reset_index(drop=True)
    rna = example_dir / "rna_count_matrix.cluster1.csv.gz"
    cells = choose_cells(_rna_cells(rna), None)
    M = build_atac_matrix(peaks, example_dir / "atac_fragments.chr22p.tsv.gz", cells, count_mode)
    key = build_gene_key(parse_gtf_genes(gtf), read_tss500(tss_path))
    got, gex, n_umi, n_cells = compute_kendall(pairs, M, list(peaks["PeakName"]), cells, rna, key)
    m = exp.merge(got, on="PairName", suffixes=("_upstream", "_port"), how="outer", indicator=True)
    both = m[m["_merge"] == "both"]
    ku, kp = both["Kendall_upstream"].to_numpy(float), both["Kendall_port"].to_numpy(float)
    na_agree = float((np.isnan(ku) == np.isnan(kp)).mean()) if len(both) else float("nan")
    ok = ~np.isnan(ku) & ~np.isnan(kp)
    res = {"pairs_upstream": int(len(exp)), "pairs_port": int(len(got)), "pairs_matched": int(len(both)),
           "kendall_defined_upstream": int((~np.isnan(ku)).sum()), "kendall_defined_port": int((~np.isnan(kp)).sum()),
           "na_pattern_agreement": na_agree, "kendall_max_abs_diff": float(np.max(np.abs(ku[ok] - kp[ok]))) if ok.any() else None,
           "kendall_pearson": float(np.corrcoef(ku[ok], kp[ok])[0, 1]) if ok.sum() > 2 else None,
           "cells": int(n_cells), "umi_count": float(n_umi), "count_mode": count_mode}
    for c in GEX_COLS:
        a, b = both[c + "_upstream"].to_numpy(float), both[c + "_port"].to_numpy(float)
        res[f"{c}_max_abs_diff"] = float(np.nanmax(np.abs(a - b))) if len(a) else None
    return res


def cmd_validate_example(args: argparse.Namespace) -> int:
    out = out_dir_for(args, "validate_example")
    ex = Path(args.example_dir) if args.example_dir else RES_DIR / "example"
    need = [ex / n for n in ("atac_fragments.chr22p.tsv.gz", "rna_count_matrix.cluster1.csv.gz", "Pairs.Kendall.expected.tsv.gz")]
    gtf = resource(args.gtf, GTF_RESOURCE[0], required=False)
    if any(not p.exists() for p in need) or gtf is None:
        print("upstream chr22 fixture not found; fetch it with `igvfagent sce2g-pipeline setup --example --gtf` "
              f"(looked in {ex})")
        return 0
    tss = resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed")
    res = validate_example(ex, gtf, tss, args.count_mode)
    lines = ["# Kendall / RNA features vs the upstream chr22 fixture", "",
             "tests/expected_output/powerlaw_v3_models/K562_cluster1_chr22p/Kendall/Pairs.Kendall.tsv.gz "
             f"({UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}), recomputed from the same fragments and RNA matrix.", "",
             md_table(["item", "value"], [[k, _fmt(v, 6) if isinstance(v, float) else v] for k, v in res.items()])]
    write_summary(out, res, lines)
    return 0


# ---------------------------------------------------------------------------
# validate-release: recompute a released IGVF scE2G prediction set from its own feature columns
# ---------------------------------------------------------------------------

def release_features(pairs, elements, genes, tss, sizes):
    """Rebuild the model features of an IGVF scE2G release (pairs / elements / genes files) with this port."""
    pd = _pd()
    pr = pairs.rename(columns={"ElementChr": "chr", "ElementStart": "start", "ElementEnd": "end", "GeneSymbol": "TargetGene",
                               "GeneTSS": "TargetGeneTSS", "GeneEnsemblID": "TargetGeneEnsembl_ID", "ElementName": "name"})
    pr = pr.reset_index(drop=True)
    pr["distance"] = ((pr["start"] + pr["end"]) / 2 - pr["TargetGeneTSS"]).abs()
    pr["activity_base"] = 1.0  # sumNearbyEnhancers is not a v3 feature; only the count is used
    enh = elements.rename(columns={"ElementChr": "chr", "ElementStart": "start", "ElementEnd": "end", "ElementName": "name"})
    enh = enh[enh["chr"].isin(set(pr["chr"]))]
    nc = num_candidate_enh_gene(pr)
    nt = num_tss_enh_gene(pr, tss)
    n5, _ = nearby_enhancers(pr, enh, sizes)
    g = genes.drop_duplicates("GeneSymbol").set_index("GeneSymbol")
    f = pd.DataFrame({"chr": pr["chr"].to_numpy(), "TargetGene": pr["TargetGene"].to_numpy()})
    f["numCandidateEnhGene"] = nc["NumCandidateEnhGene"].to_numpy()
    f["numTSSEnhGene"] = nt["count"].to_numpy()
    f["numNearbyEnhancers"] = pr["name"].map(dict(zip(n5["name"], n5["count"]))).fillna(0).to_numpy()
    for c in ("normalizedATAC_prom", "ubiqExpressed", "RNA_pseudobulkTPM"):
        if c in g.columns:
            f[c] = pr["TargetGene"].map(g[c]).to_numpy()
    for c in ("ARC.E2G.Score", "ABC.Score"):
        if c in pr.columns:
            f[c] = pr[c].to_numpy()
    return f


def cmd_validate_release(args: argparse.Namespace) -> int:
    np = _np()
    out = out_dir_for(args, "validate_release")
    pairs = read_tsv(args.pairs, comment="#")
    if args.chrom:
        pairs = pairs[pairs["ElementChr"] == args.chrom]
    elements, genes = read_tsv(args.elements, comment="#"), read_tsv(args.genes, comment="#")
    tss = read_tss500(resource(args.tss, "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed"))
    sizes = dict(read_chrom_sizes(resource(args.chrom_sizes, "GRCh38_EBV.no_alt.chrom.sizes.tsv")))
    f = release_features(pairs, elements, genes, tss, sizes)
    m = resolve_model(args.model, Path(args.models_dir) if args.models_dir else None)
    pred, _ = apply_model(f, m, False, qnorm=m["qnorm"] is not None)
    from scipy.stats import spearmanr  # type: ignore
    res = {"model": m["name"], "pairs": int(len(pred)), "chrom": args.chrom or "all",
           "qnorm_reference": str(m["qnorm"]) if m["qnorm"] else None}
    ign = SCORE_BASE + ".qnorm.ignoreTPM" if SCORE_BASE + ".qnorm.ignoreTPM" in pred.columns else FINAL_SCORE_COL
    for mine, rel in ((SCORE_BASE, "Score.ignoreTPM"), (ign, "Score.ignoreTPM"), (FINAL_SCORE_COL, "Score")):
        if rel not in pairs.columns:
            continue
        a, b = pred[mine].to_numpy(float), pairs[rel].to_numpy(float)
        if mine == SCORE_BASE:
            res["spearman_raw_vs_" + rel] = float(spearmanr(a, b).correlation)
            continue
        d = np.abs(a - b)
        res[f"{mine}_vs_{rel}_max_abs_diff"] = float(d.max()) if len(d) else None
        res[f"{mine}_vs_{rel}_frac_within_1e-6"] = float((d < 1e-6).mean()) if len(d) else None
    if args.write_predictions:
        write_tsv(pred.assign(**{"release_" + c: pairs[c].to_numpy() for c in ("Score", "Score.ignoreTPM") if c in pairs.columns}),
                  out / "recomputed_predictions.tsv.gz")
    lines = ["# Recomputing a released scE2G prediction set", "",
             "Features rebuilt from the release's element, gene and pair files (numCandidateEnhGene, numTSSEnhGene, "
             "numNearbyEnhancers recomputed; ARC.E2G.Score, normalizedATAC_prom, ubiqExpressed and RNA_pseudobulkTPM taken "
             "as released), scored with the embedded weights, qnorm and TPM filter.", "",
             md_table(["item", "value"], [[k, _fmt(v, 6) if isinstance(v, float) else v] for k, v in res.items()])]
    write_summary(out, res, lines)
    return 0


# ---------------------------------------------------------------------------
# selftest: synthetic multiome clusters with planted enhancers
# ---------------------------------------------------------------------------

def _synthetic(tmp: Path, seed: int = 7) -> dict:
    """Two chromosomes, 16 genes, 240 cells in 3 states; per gene 2 planted enhancers whose accessibility follows
    the gene's expression state and 3 null peaks; ABC-shaped outputs, GTF, TSS500, CRISPR truth, QC reference."""
    pd, np = _pd(), _np()
    rng = np.random.RandomState(seed)
    sizes = [("chr1", 3_000_000), ("chr2", 2_000_000)]
    (tmp / "ref").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sizes).to_csv(tmp / "ref" / "chrom.sizes", sep="\t", header=False, index=False)
    genes = []
    for ci, (c, L) in enumerate(sizes):
        for j in range(8):
            tss = 300_000 + j * 300_000 if c == "chr1" else 200_000 + j * 220_000
            gi = len(genes)
            genes.append({"chr": c, "tss": tss, "name": f"GENE{gi}", "ens": f"ENSG{gi:011d}", "state": gi % 3})
    G = pd.DataFrame(genes)
    tss = pd.DataFrame({"chr": G["chr"], "start": G["tss"] - 250, "end": G["tss"] + 250, "name": G["name"], "score": 0,
                        "strand": "+", "Ensembl_ID": G["ens"], "gene_type": "protein_coding"})
    tss.to_csv(tmp / "ref" / "tss500.bed", sep="\t", header=False, index=False)
    # GTF: GENE5 is called OLD5 in the RNA matrix; GENE7's Ensembl ID has two gene names (dropped from the 1:1 key)
    gtf_lines = []
    for g in genes:
        nm = "OLD5" if g["name"] == "GENE5" else g["name"]
        attrs = f'gene_id "{g["ens"]}.3"; gene_type "protein_coding"; gene_name "{nm}"; level 2;'
        gtf_lines.append("\t".join([g["chr"], "HAVANA", "gene", str(g["tss"]), str(g["tss"] + 5000), ".", "+", ".", attrs]))
        gtf_lines.append("\t".join([g["chr"], "HAVANA", "transcript", str(g["tss"]), str(g["tss"] + 5000), ".", "+", ".", attrs]))
        if g["name"] == "GENE7":
            gtf_lines.append("\t".join([g["chr"], "HAVANA", "gene", str(g["tss"]), str(g["tss"] + 10), ".", "+", ".",
                                        f'gene_id "{g["ens"]}.3"; gene_name "GENE7B";']))
    with gzip.open(tmp / "ref" / "genes.gtf.gz", "wt") as fh:
        fh.write("##description: synthetic\n" + "\n".join(gtf_lines) + "\n")
    n_cells, per = 240, 80
    state = np.repeat([0, 1, 2], per)
    cells = [f"c{i:03d}-1" for i in range(n_cells)]
    lib = rng.uniform(0.6, 1.6, n_cells)
    # peaks: promoter, 2 planted enhancers, 3 null peaks per gene
    peaks, truth = [], []
    for g in genes:
        for kind, off in (("promoter", 0), ("enh", -40_000), ("enh", 65_000), ("null", -90_000), ("null", 25_000), ("null", 110_000)):
            ctr = g["tss"] + off + int(rng.randint(-300, 300)) if kind != "promoter" else g["tss"]
            w = 500 if kind == "promoter" else int(rng.randint(300, 700))
            peaks.append({"chr": g["chr"], "start": ctr - w // 2, "end": ctr + w // 2, "gene": g["name"], "kind": kind,
                          "state": g["state"]})
    P = pd.DataFrame(peaks)
    acc = np.zeros((len(P), n_cells), dtype=bool)
    for i, r in P.iterrows():
        if r["kind"] == "enh":
            p = np.where(state == r["state"], 0.8, 0.08)
        elif r["kind"] == "promoter":
            p = np.full(n_cells, 0.5)
        else:
            p = np.full(n_cells, 0.2)
        acc[i] = rng.uniform(size=n_cells) < p
    # RNA counts: state-specific genes; OLD5 naming; MTJUNK outside the gene universe
    rna_rows, rna_names = [], []
    for g in genes:
        lam = np.where(state == g["state"], 6.0, 0.6) * lib
        rna_rows.append(rng.poisson(lam))
        rna_names.append("OLD5" if g["name"] == "GENE5" else g["name"])
    rna_rows.append(rng.poisson(3.0 * lib))
    rna_names.append("MTJUNK")
    R = pd.DataFrame(np.vstack(rna_rows), index=rna_names, columns=cells)
    R.to_csv(tmp / "rna.csv.gz")
    frags = []
    for i, r in P.iterrows():
        for ci in np.flatnonzero(acc[i]):
            for _ in range(1 + int(rng.uniform() < 0.3)):
                s = int(rng.randint(r["start"], r["end"]))
                frags.append((r["chr"], s, s + int(rng.randint(60, 250)), cells[ci], 1))
    for ci in range(n_cells):
        for _ in range(15):
            c, L = sizes[rng.randint(2)]
            s = int(rng.randint(0, L - 400))
            frags.append((c, s, s + int(rng.randint(60, 300)), cells[ci], 1))
    frags += [("chrUn_junk", 100 + k, 300 + k, cells[k], 1) for k in range(5)]
    F = pd.DataFrame(frags, columns=["chr", "start", "end", "barcode", "count"])
    order = {c: i for i, c in enumerate(["chr1", "chr2", "chrUn_junk"])}
    F = F.assign(_o=F["chr"].map(order)).sort_values(["_o", "start", "end"], kind="mergesort").drop(columns="_o")
    with gzip.open(tmp / "fragments.tsv.gz", "wt") as fh:
        fh.write("# id=synthetic\n")
        F.to_csv(fh, sep="\t", header=False, index=False)
    # ABC-shaped outputs: candidate elements = peaks +/- 100 bp; pairs with genes within 1 Mb
    abc_dir = tmp / "abc"
    for sub_ in ("Peaks", "Predictions", "Neighborhoods"):
        (abc_dir / sub_).mkdir(parents=True, exist_ok=True)
    npk = P.assign(name=[f"peak{i}" for i in range(len(P))], score=100, strand=".", signalValue=5.0, pValue=10.0, qValue=8.0,
                   peak=P["end"].sub(P["start"]) // 2)
    npk[["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]].sort_values(
        ["chr", "start"]).to_csv(abc_dir / "Peaks" / "macs2_peaks.narrowPeak.sorted", sep="\t", header=False, index=False)
    E = P.assign(start=P["start"] - 100, end=P["end"] + 100)
    E["class"] = np.where(E["kind"] == "promoter", "promoter", "intergenic")
    E["name"] = E["class"] + "|" + E["chr"] + ":" + E["start"].astype(str) + "-" + E["end"].astype(str)
    E["activity_base"] = acc.sum(axis=1) / n_cells * 10 + 0.1
    rows = []
    for g in genes:
        sub_ = E[(E["chr"] == g["chr"]) & ((E["start"] + E["end"]) / 2 - g["tss"]).abs().le(1_000_000)]
        mid = (sub_["start"] + sub_["end"]) / 2
        dist = (mid - g["tss"]).abs()
        contact = (dist.clip(lower=5000) / 5000) ** -0.87
        num = sub_["activity_base"] * contact
        prom_act = float(E.loc[(E["gene"] == g["name"]) & (E["kind"] == "promoter"), "activity_base"].iloc[0])
        for (idx, e), d, cc, nu in zip(sub_.iterrows(), dist, contact, num):
            rows.append({"chr": e["chr"], "start": e["start"], "end": e["end"], "name": e["name"], "class": e["class"],
                         "activity_base": e["activity_base"], "TargetGene": g["name"], "TargetGeneTSS": g["tss"],
                         "TargetGeneEnsembl_ID": g["ens"], "TargetGeneIsExpressed": True,
                         "isSelfPromoter": bool(e["kind"] == "promoter" and e["gene"] == g["name"]), "isSelfGenic": False,
                         "CellType": "synth", "distance": d, "hic_contact_pl_scaled_adj": cc, "ABC.Score.Numerator": nu,
                         "normalized_atac_enh": e["activity_base"], "normalized_atac_prom": prom_act,
                         "_gene_state": g["state"], "_kind": e["kind"], "_peak_gene": e["gene"]})
    A = pd.DataFrame(rows)
    A["ABC.Score"] = A["ABC.Score.Numerator"] / A.groupby("TargetGene")["ABC.Score.Numerator"].transform("sum")
    A["powerlaw.Score"] = A["ABC.Score"]
    truth = A[["chr", "start", "end", "TargetGene", "_kind", "_peak_gene"]].copy()
    A.drop(columns=["_gene_state", "_kind", "_peak_gene"]).to_csv(abc_dir / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz",
                                                                  sep="\t", index=False)
    E[["chr", "start", "end", "name", "class", "activity_base"]].to_csv(abc_dir / "Neighborhoods" / "EnhancerList.txt", sep="\t",
                                                                          index=False)
    GL = pd.DataFrame({"chr": G["chr"], "start": G["tss"] - 250, "end": G["tss"] + 250, "name": G["name"], "score": 0,
                       "strand": "+", "Ensembl_ID": G["ens"], "gene_type": "protein_coding", "symbol": G["name"], "tss": G["tss"],
                       "Expression": 1.0, "is_ue": False, "cellType": "synth", "ATAC.RPKM.quantile": 0.5,
                       "PromoterActivityQuantile": 0.8})
    GL.to_csv(abc_dir / "Neighborhoods" / "GeneList.txt", sep="\t", index=False)
    pd.DataFrame({"TargetGene": G["name"], "gene_chr": G["chr"], "gene_tss": G["tss"], "gene_tss500bp": G["tss"] + 250,
                  "is_ubiquitous_uniform": [i % 4 == 0 for i in range(len(G))], "P2PromoterClass": False}).to_csv(
        tmp / "ref" / "gene_classes.tsv", sep="\t", index=False)
    # CRISPR truth: planted enhancers regulated, null peaks not; one gene outside the TSS universe
    cr = []
    for i, r in P.iterrows():
        if r["kind"] == "promoter":
            continue
        g = G[G["name"] == r["gene"]].iloc[0]
        cr.append({"dataset": "Synth", "chrom": r["chr"], "chromStart": r["start"], "chromEnd": r["end"],
                   "name": f"{r['gene']}|{r['chr']}:{r['start']}-{r['end']}", "EffectSize": -0.3 if r["kind"] == "enh" else 0.0,
                   "chrTSS": r["chr"], "startTSS": g["tss"], "endTSS": g["tss"] + 1, "measuredGeneSymbol": r["gene"],
                   "Significant": r["kind"] == "enh", "pValueAdjusted": 0.01, "PowerAtEffectSize25": 0.9,
                   "ValidConnection": "TRUE", "CellType": "K562", "Reference": "synthetic", "Regulated": r["kind"] == "enh",
                   "pair_uid": i, "merged_uid": i, "merged_start": r["start"], "merged_end": r["end"]})
    cr.append(dict(cr[0], measuredGeneSymbol="NOTAGENE", name="NOTAGENE|x"))
    cr.append(dict(cr[1], chromStart=2_990_000, chromEnd=2_990_500, name="nowhere", Regulated=False))
    C = pd.DataFrame(cr)
    C.to_csv(tmp / "crispr.tsv.gz", sep="\t", index=False)
    ref = pd.DataFrame([{**{k: v for k, v in zip(STATS_COLS[:9], [36435, 8411, 49924, 1.37, 5.9, 88770, 737, 14389, 3578])},
                         "cluster": f"ref{i}", "model_name": m, "fragments_total": 5e7 * (1 + i), "cell_count": 3000,
                         "umi_count": 5e7} for i in range(3) for m in ("multiome_powerlaw_v3", "scATAC_powerlaw_v3")] +
                       [{**{k: 0 for k in STATS_COLS[:9]}, "cluster": "tiny", "model_name": "multiome_powerlaw_v3",
                         "fragments_total": 1e5, "cell_count": 10, "umi_count": 1e4}])
    ref.to_csv(tmp / "ref" / "qc_reference.tsv", sep="\t", index=False)
    # synthetic qnorm references (the real ones are data fetched by setup)
    md = tmp / "models"
    for m in MODEL_SPECS:
        (md / m).mkdir(parents=True, exist_ok=True)
        q = np.linspace(0, 1, 101)
        pd.DataFrame({"quantile": q, "reference_score": q ** 2}).to_csv(md / m / "qnorm_reference.tsv.gz", sep="\t", index=False)
    return {"sizes": tmp / "ref" / "chrom.sizes", "tss": tmp / "ref" / "tss500.bed", "gtf": tmp / "ref" / "genes.gtf.gz",
            "rna": tmp / "rna.csv.gz", "fragments": tmp / "fragments.tsv.gz", "abc_dir": abc_dir,
            "gene_classes": tmp / "ref" / "gene_classes.tsv", "crispr": tmp / "crispr.tsv.gz", "qc_ref": tmp / "ref" / "qc_reference.tsv",
            "models_dir": md, "peaks": P, "acc": acc, "cells": cells, "frags": F, "abc": A, "truth": truth, "genes": G,
            "rna_df": R, "state": state}


def cmd_selftest(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    import io
    from contextlib import redirect_stdout
    checks: "List[Tuple[bool, str]]" = []

    def check(ok, msg):
        checks.append((bool(ok), msg))
        print(("  ok    " if ok else "  FAIL  ") + msg)

    def run(argv) -> int:
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(argv)
        if rc != 0:
            print(buf.getvalue())
        return rc

    t0 = time.time()
    out_root_existed = OUT_ROOT.exists()
    tmp = Path(tempfile.mkdtemp(prefix="sce2g_pipeline_selftest_"))
    made: List[Path] = []
    try:
        S = _synthetic(tmp)
        noplot = ["--no-plots"] if args.no_plots else []
        common_ref = ["--chrom-sizes", str(S["sizes"]), "--tss", str(S["tss"])]
        md = ["--models-dir", str(S["models_dir"])]

        # --- upstream unit tests (tests/unit/python) -------------------------------------------------
        q = calculate_quantiles(pd.Series([0, 0, 0.2, 0.5, 0.9]))
        check(np.allclose(q.to_numpy(), [0, 0, 0.5, 0.75, 1.0]), "calculate_quantiles == upstream unit test")
        ref = pd.DataFrame({"quantile": [0.0, 0.25, 0.5, 0.75, 1.0], "reference_score": [0.0, 0.1, 0.3, 0.6, 1.0]})
        o = qnorm_scores(pd.DataFrame({"E2G.Score": [0, 0, 0.2, 0.5, 0.9]}), ref, False)
        check(np.allclose(o["E2G.Score.qnorm"], [0, 0, 0.3, 0.6, 1.0]) and "E2G.Score.cv.qnorm" not in o, "qnorm_scores interpolation")
        o = qnorm_scores(pd.DataFrame({"E2G.Score": [0, 0, 0.2, 0.5, 0.9], "E2G.Score.cv": [0, 0.1, 0.2, 0.5, 0.9]}), ref, True)
        check(np.allclose(o["E2G.Score.cv.qnorm"], [0, 0.1, 0.3, 0.6, 1.0]), "qnorm of the cv column")
        d0 = pd.DataFrame({"E2G.Score.qnorm": [0.1, 0.5, 0.9]})
        check(filter_by_tpm(d0.copy(), 1.0).equals(d0), "TPM filter no-op without RNA_pseudobulkTPM")
        d1 = pd.DataFrame({"E2G.Score.qnorm": [0.1, 0.5, 0.9], "RNA_pseudobulkTPM": [0.5, 2.0, 5.0]})
        check(filter_by_tpm(d1.copy(), 0).equals(d1), "TPM filter no-op at threshold 0")
        o = filter_by_tpm(d1.copy(), 1.0)
        check(np.allclose(o["E2G.Score.qnorm"], [0, 0.5, 0.9]) and np.allclose(o["E2G.Score.qnorm.ignoreTPM"], [0.1, 0.5, 0.9]),
              "TPM filter zeroes genes below threshold, keeps .ignoreTPM")
        yt, yp = np.array([0, 1, 0, 1, 1]), np.array([0.2, 0.4, 0.5, 0.6, 0.9])
        check(auc_mod(np.array([1.0]), np.array([0.5])) == 0 and abs(auc_mod(np.array([0.0, 1.0]), np.array([1.0, 1.0])) - 1) < 1e-12,
              "auc_mod single point / two points")
        check(abs(statistic_aupr(np.array([0, 0, 1, 1]), np.array([.1, .2, .8, .9])) - 1) < 1e-12 and
              abs(statistic_aupr(yt, yp) - 0.9027777777777777) < 1e-12, "statistic_aupr == upstream unit tests")
        check(abs(threshold_at_target_recall(yt, yp, 0.5) - 0.5) < 1e-12 and threshold_at_target_recall(yt, yp, 1.5) is None,
              "threshold_at_target_recall")
        check(statistic_precision_at_threshold(yt, yp, None) == 0 and abs(statistic_precision_at_threshold(yt, yp, 0.5) - 2 / 3) < 1e-12
              and statistic_recall_at_threshold(yt, yp, 10.0) == 0, "precision / recall at threshold")

        # --- embedded models --------------------------------------------------------------------------
        m = resolve_model("multiome_powerlaw_v3", S["models_dir"])
        check(m["features"] == ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed",
                                "numCandidateEnhGene", "ARC.E2G.Score"] and abs(m["full"][0] - 2.21949242) < 1e-7
              and abs(m["full"][6] - 1.72645456) < 1e-7 and len(m["cv"]) == 23 and m["tpm_threshold"] == 1.0,
              "embedded multiome_powerlaw_v3: features, pickle intercept/ARC weight, 23 cv models, TPM 1")
        check(all(len(MODEL_WEIGHTS[k]["cv"]) == 23 and len(MODEL_WEIGHTS[k]["full"]) == 7 for k in MODEL_SPECS) and
              [model_threshold_str(MODEL_SPECS[k]["score_threshold"]) for k in MODEL_SPECS] == ["0.177", "0.173", "0.174", "0.187"],
              "four v3 models with thresholds 0.177 / 0.173 / 0.174 / 0.187")
        try:
            from sklearn.linear_model import LogisticRegression  # type: ignore
            X = pd.DataFrame(np.abs(np.random.RandomState(0).normal(1, 1, (30, 6))), columns=m["features"])
            lr = LogisticRegression()
            lr.classes_ = np.array([0, 1])
            lr.coef_ = np.array([m["full"][1:]])
            lr.intercept_ = np.array([m["full"][0]])
            ref_p = lr.predict_proba(np.log(np.abs(X.to_numpy()) + 0.01))[:, 1]
            check(np.allclose(logistic_score(transform_features(X, m["features"]), m["full"]), ref_p), "logistic score == sklearn predict_proba")
        except ImportError:
            check(True, "sklearn missing: predict_proba comparison skipped")
        rc = run(["models", "--out-dir", str(tmp / "o_models"), "--export", str(tmp / "exported")] + md)
        mm = resolve_model(str(tmp / "exported" / "scATAC_megamap_v3"), S["models_dir"])
        check(rc == 0 and mm["score_threshold"] == ".187" and mm["full"] == MODEL_WEIGHTS["scATAC_megamap_v3"]["full"]
              and mm["qnorm"] is not None, "models --export round-trips through a model directory")

        # --- Kendall tau-b ---------------------------------------------------------------------------
        from scipy.stats import kendalltau  # type: ignore
        rs = np.random.RandomState(3)
        okk = True
        for _ in range(20):
            x = rs.poisson(1.5, 60).astype(float)
            y = rs.uniform(size=60) < 0.4
            a, b = kendall_tau_b_binary(x, y), kendalltau(x, y.astype(float)).correlation
            okk &= abs(a - b) < 1e-12
        check(okk, "closed-form Kendall tau-b == scipy.stats.kendalltau (ties in x, binary y)")
        check(math.isnan(kendall_tau_b_binary(np.ones(10), np.arange(10) % 2)) and
              math.isnan(kendall_tau_b_binary(np.arange(10.0), np.zeros(10))), "Kendall NA for a constant margin")

        # --- frag-to-tagalign --------------------------------------------------------------------------
        o = tmp / "o_tag"
        rc = run(["frag-to-tagalign", "--fragments", str(S["fragments"]), "--chrom-sizes", str(S["sizes"]), "--cluster", "c1",
                  "--out-dir", str(o)])
        Fk = S["frags"][S["frags"]["chr"].isin(["chr1", "chr2"])]
        n_frag = int(Path(o / "c1" / "fragment_count.txt").read_text())
        tags = pd.read_csv(o / "c1" / "tagAlign" / "tagAlign.sort.gz", sep="\t", header=None)
        check(rc == 0 and n_frag == len(Fk) and len(tags) == 2 * len(Fk) and "chrUn_junk" not in set(tags[0]),
              "tagAlign: 2 tags per fragment, off-list chromosome dropped, fragment_count")
        f0 = Fk.iloc[0]
        mid = (f0["start"] + f0["end"]) // 2
        t_ok = ((tags[0] == f0["chr"]) & (tags[1] == f0["start"]) & (tags[2] == mid) & (tags[5] == "+")).any() and \
            ((tags[0] == f0["chr"]) & (tags[1] == mid + 1) & (tags[2] == f0["end"]) & (tags[5] == "-")).any()
        srt = all((sub_[1].diff().dropna() >= 0).all() for _, sub_ in tags.groupby(0)) and list(dict.fromkeys(tags[0])) == ["chr1", "chr2"]
        check(t_ok and srt, "tag split at int((s+e)/2) (+: s..m, -: m+1..e), sorted in chrom-sizes order")
        bcs = Path(o / "c1" / "Kendall" / "cell_barcodes.txt").read_text().split()
        check(bcs == list(dict.fromkeys(Fk["barcode"])), "cell_barcodes = unique barcodes in file order")
        check((o / "c1" / "tagAlign" / "tagAlign.sort.gz.tbi").exists() or (o / "c1" / "tagAlign" / "tagAlign.sort.gz").exists(),
              "tagAlign bgzipped (+ tabix index with pysam)")

        # --- bigwig / coverage -------------------------------------------------------------------------
        rc = run(["bigwig", "--fragments", str(S["fragments"]), "--chrom-sizes", str(S["sizes"]), "--raw", "--cluster", "c1",
                  "--out-dir", str(tmp / "o_bw")])
        bg = pd.read_csv(tmp / "o_bw" / "c1" / "ATAC.bg", sep="\t", header=None)
        bgn = pd.read_csv(tmp / "o_bw" / "c1" / "ATAC_norm.bg", sep="\t", header=None)
        tot = float(((bg[2] - bg[1]) * bg[3]).sum())
        check(rc == 0 and abs(tot - float((Fk["end"] - Fk["start"]).sum())) < 1e-6 and
              np.allclose(bgn[3].to_numpy(), bg[3].to_numpy() * 1e6 / len(Fk)), "genomecov -bg: coverage integrates to fragment bp; norm = x 1e6/N")

        # --- run: the whole pipeline -------------------------------------------------------------------
        o_run = tmp / "o_run"
        rc = run(["run", "--cluster", "synth", "--fragments", str(S["fragments"]), "--rna-matrix", str(S["rna"]), "--abc-dir",
                  str(S["abc_dir"]), "--models", "multiome_powerlaw_v3", "scATAC_powerlaw_v3", "--gtf", str(S["gtf"]),
                  "--gene-classes", str(S["gene_classes"]), "--qc-reference", str(S["qc_ref"]), "--crispr", str(S["crispr"]),
                  "--n-boot", "30", "--ci-method", "percentile", "--igv-tracks", "--out-dir", str(o_run)] + common_ref + md + noplot)
        cd = o_run / "synth"
        check(rc == 0 and (o_run / "summary.json").exists() and (o_run / "report.md").exists(), "run completes (report.md, summary.json)")
        kd = read_tsv(cd / "Kendall" / "Pairs.Kendall.tsv.gz")
        check(list(kd.columns) == KENDALL_OUT_COLS, "Pairs.Kendall columns as upstream")
        pairs = read_tsv(cd / "Kendall" / "Pairs.tsv.gz")
        P = S["peaks"]
        exp_n = 0
        for _, r in P.iterrows():
            a = S["abc"]
            hit = a[(a["chr"] == r["chr"]) & (a["start"] <= r["end"]) & (a["end"] >= r["start"])]
            exp_n += hit["TargetGene"].nunique()
        check(len(pairs) == exp_n and pairs["PairName"].is_unique and
              (pairs["PeakName"] == pairs["chr"] + "-" + pairs["start"].astype(str) + "-" + pairs["end"].astype(str)).all(),
              "Kendall pairs = peak x overlapping ABC pairs, deduplicated, PeakName chr-start-end")
        M, pk, cl = load_matrix(cd / "Kendall" / "atac_matrix")
        pr = P.assign(PeakName=P["chr"] + "-" + P["start"].astype(str) + "-" + P["end"].astype(str)).set_index("PeakName")
        i0 = pk.index(pr.index[3])
        r0 = pr.iloc[3]
        brute = S["frags"][(S["frags"]["chr"] == r0["chr"]) & (S["frags"]["start"] < r0["end"]) & (S["frags"]["end"] >= r0["start"])]
        bc = brute["barcode"].value_counts()
        row = M.getrow(i0).toarray().ravel()
        check(all(row[cl.index(b)] == n for b, n in bc.items()) and row.sum() == len(brute),
              "ATAC matrix counts fragments overlapping the peak (tabix semantics), per cell")
        kk = kd.copy()
        kk["_kind"] = "other"
        for _, r in P.iterrows():
            sel = (kk["chr"] == r["chr"]) & (kk["start"] == r["start"]) & (kk["TargetGene"] == r["gene"])
            kk.loc[sel, "_kind"] = r["kind"]
        planted = kk.loc[kk["_kind"] == "enh", "Kendall"].dropna()
        null = kk.loc[kk["_kind"] == "null", "Kendall"].dropna()
        check(len(planted) >= 20 and planted.mean() > 0.25 and abs(null.mean()) < 0.05,
              f"planted enhancers recovered by Kendall (mean {planted.mean():.3f} vs null {null.mean():.3f})")
        check("GENE5" in set(kd["TargetGene"]) and "OLD5" not in set(kd["TargetGene"]) and "GENE7" not in set(kd["TargetGene"]),
              "gene names mapped through the 1:1 Ensembl key (OLD5 -> GENE5; ambiguous GENE7 dropped)")
        gx = read_tsv(cd / "Kendall" / "gene_expression_metrics.tsv.gz")
        Rdf = S["rna_df"]
        g3 = Rdf.loc["GENE3"].to_numpy(float)
        norm = np.log1p(Rdf.to_numpy(float) / Rdf.to_numpy(float).sum(axis=0) * 1e4)
        r3 = gx.set_index("TargetGene").loc["GENE3"]
        check(abs(r3["RNA_pseudobulkTPM"] - g3.sum() / Rdf.to_numpy().sum() * 1e6) < 1e-6 and
              abs(r3["RNA_percentCellsDetected"] - (g3 > 0).mean()) < 1e-12 and
              abs(r3["RNA_meanLogNorm"] - norm[list(Rdf.index).index("GENE3")].mean()) < 1e-12,
              "RNA TPM / percent detected / mean LogNormalize")
        check(int(Path(cd / "umi_count.txt").read_text()) == int(Rdf.to_numpy().sum()) and int(Path(cd / "cell_count.txt").read_text()) == 240,
              "umi_count and cell_count")
        arc = read_tsv(cd / "ARC" / "EnhancerPredictionsAllPutative_ARC.tsv.gz")
        ok_ = arc["Kendall"].notna() & (arc["powerlaw.Score"] > 0)
        la, kv = np.log(arc.loc[ok_, "powerlaw.Score"]), arc.loc[ok_, "Kendall"]
        r_ = np.corrcoef(la, kv)[0, 1]
        exp_arc = np.exp(la + kv * la.std(ddof=1) / kv.std(ddof=1) * r_)
        check(np.allclose(arc.loc[ok_, "ARC.E2G.Score"], exp_arc) and np.allclose(arc.loc[~ok_, "ARC.E2G.Score"], arc.loc[~ok_, "powerlaw.Score"]),
              "ARC-E2G = exp(log ABC + Kendall sd(logABC)/sd(Kendall) r); ABC where Kendall undefined")
        j = arc[arc["Kendall"].notna()].iloc[0]
        cand = kd[(kd["chr"] == j["chr"]) & (kd["TargetGene"] == j["TargetGene"]) & (kd["start"] <= j["end"]) & (kd["end"] >= j["start"])]
        check(abs(cand["Kendall"].max() - j["Kendall"]) < 1e-12 and list(arc.columns[:5]) == ["chr", "start", "end", "width", "strand"],
              "ARC Kendall = max over overlapping peaks of the same gene; GRanges columns")
        abc = S["abc"]
        nc = read_tsv(cd / "new_features" / "NumCandidateEnhGene.tsv")
        g0 = abc[abc["TargetGene"] == "GENE2"].assign(mid=lambda d: ((d["start"] + d["end"]) / 2).astype(int))
        up = g0[g0["mid"] < 900_000].sort_values("mid", ascending=False)
        exp_rank = dict(zip(up["name"], range(1, len(up) + 1)))
        got = dict(zip(nc.loc[nc["TargetGene"] == "GENE2", "name"], nc.loc[nc["TargetGene"] == "GENE2", "NumCandidateEnhGene"]))
        check(all(got[k] == v for k, v in exp_rank.items()), "numCandidateEnhGene counts outward from the TSS")
        nt = read_tsv(cd / "new_features" / "NumTSSEnhGene.tsv")
        e1 = abc.iloc[5]
        m1 = int((e1["start"] + e1["end"]) / 2)
        lo, hi = (e1["TargetGeneTSS"], e1["end"]) if e1["TargetGeneTSS"] < m1 else (e1["start"], int(m1 + e1["distance"]))
        tss_df = read_tss500(S["tss"])
        bn = int(((tss_df["chr"] == e1["chr"]) & (tss_df["start"] < hi) & (tss_df["end"] > lo)).sum())
        check(int(nt.iloc[5]["count"]) == bn and (nt["count"] >= 1).all(), "numTSSEnhGene = TSS500 windows in the E-P span (>= own TSS)")
        n5 = pd.read_csv(cd / "new_features" / "NumEnhancersEG5kb.txt", sep="\t", header=None, names=["name", "n"])
        el = pd.read_csv(S["abc_dir"] / "Neighborhoods" / "EnhancerList.txt", sep="\t")
        e2 = el.iloc[1]
        c2 = int((e2["start"] + e2["end"]) / 2)
        nb = el[(el["chr"] == e2["chr"]) & (el["start"] < c2 + 5000) & (el["end"] > c2 - 5000) & (el["name"] != e2["name"])]
        got_n = dict(zip(n5["name"], n5["n"])).get(e2["name"], 0)
        check(got_n == len(nb), "numNearbyEnhancers = other elements within the midpoint +/- 5 kb")
        gw = read_tsv(cd / "genomewide_features.tsv.gz")
        want = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene", "ARC.E2G.Score",
                "ABC.Score", "RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected", "Kendall", "normalizedATAC_enh"]
        check(all(c in gw.columns for c in want) and "powerlaw.Score" not in gw.columns and gw[want].notna().all().all(),
              "genome-wide features renamed to model names, NA filled")
        ft_run = read_tsv(cd / "feature_table.tsv")
        check(set(ft_run["feature"]) >= {"RNA_pseudobulkTPM", "Kendall", "normalizedATAC_enh"} and (ft_run["feature"] == "ABC.Score").sum() == 1,
              "biosample feature table = model tables + ARC rows (ABC.Score not duplicated)")
        pred = read_tsv(cd / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz")
        ref_s = logistic_score(transform_features(gw, m["features"]), m["full"])
        check(np.allclose(pred["E2G.Score"], ref_s) and "E2G.Score.cv" in pred.columns and "E2G.Score.qnorm.ignoreTPM" in pred.columns,
              "predictions: E2G.Score, cv (CRISPR benchmarking), qnorm, ignoreTPM")
        low = pred["RNA_pseudobulkTPM"] < 1
        check((pred.loc[low, "E2G.Score.qnorm"] == 0).all() and (pred.loc[~low, "E2G.Score.qnorm"] == pred.loc[~low, "E2G.Score.qnorm.ignoreTPM"]).all(),
              "TPM < 1 zeroes E2G.Score.qnorm")
        thr = read_tsv(cd / "multiome_powerlaw_v3" / "scE2G_predictions_threshold0.177.tsv.gz")
        check((thr["E2G.Score.qnorm"] >= 0.177).all() and ((thr["class"] != "promoter") | thr["isSelfPromoter"].map(_truthy)).all()
              and len(thr) == int(((pred["E2G.Score.qnorm"] >= 0.177) & ((pred["class"] != "promoter") | pred["isSelfPromoter"].map(_truthy))).sum()),
              "thresholded links: qnorm >= 0.177, non-self promoters removed")
        bp = pd.read_csv(cd / "multiome_powerlaw_v3" / "scE2G_predictions_threshold0.177.bedpe", sep="\t", header=None)
        check(len(bp) == len(thr.drop_duplicates()) and bp.shape[1] == 10, "bedpe written (10 columns)")
        sat = read_tsv(cd / "scATAC_powerlaw_v3" / "scE2G_predictions.tsv.gz")
        check("E2G.Score.qnorm.ignoreTPM" not in sat.columns, "scATAC model: no TPM filter (tpm_threshold 0)")
        pl = read_tsv(cd / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz")
        gl = read_tsv(cd / "multiome_powerlaw_v3" / "scE2G_gene_list.tsv.gz")
        check("removed_by_TPM" in gl.columns and "considered_for_predictions" in gl.columns and "ATAC.RPKM.quantile" not in gl.columns
              and "is_ue" not in gl.columns, "gene list: TPM / promoter flags, RPKM + is_ue columns dropped")
        st = read_tsv(cd / "multiome_powerlaw_v3" / "scE2G_predictions_threshold0.177_stats.tsv")
        tnp = thr[thr["class"] != "promoter"]
        check(list(st.columns) == STATS_COLS and int(st["n_enh_gene_links"].iloc[0]) == len(tnp) and
              int(st["n_genes_not_expressed"].iloc[0]) == pl.loc[pl["E2G.Score.qnorm"] == 0, "TargetGene"].nunique()
              and int(st["cell_count"].iloc[0]) == 240, "stats: links, genes with score 0, counts")
        qc = read_tsv(o_run / "qc_plots" / "all_qc_stats.tsv")
        warn = (o_run / "qc_plots" / "qc_warnings.txt").read_text()
        comp = read_tsv(o_run / "qc_plots" / "qc_reference_comparison.tsv")
        check(len(qc) == 2 and "fragments_total" in warn and len(comp) == 18 and comp["reference_n"].eq(3).all(),
              "QC: warnings below 2e6 fragments; comparison vs reference clusters with >= 2e6 fragments")
        bt = read_tsv(o_run / "crispr_benchmarking_performance_summary.tsv")
        check(set(bt["model"]) == {"multiome_powerlaw_v3", "scATAC_powerlaw_v3", "scATAC_ABC", "multiome_powerlaw_v3_noTPMfilter",
                                    "distanceToTSS"} and list(bt.columns) == BENCH_COLS, "benchmark rows as upstream (+ABC, noTPM, distance)")
        au = bt.set_index("model")["AUPRC"]
        cf = read_tsv(cd / "multiome_powerlaw_v3" / "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz")
        base = cf["Regulated"].map(_truthy).mean()
        check(au["multiome_powerlaw_v3_noTPMfilter"] > base + 0.2 and au["multiome_powerlaw_v3_noTPMfilter"] > au["distanceToTSS"]
              and au["multiome_powerlaw_v3_noTPMfilter"] > au["scATAC_ABC"],
              f"multiome scE2G recovers planted CRISPR positives (AUPRC {au['multiome_powerlaw_v3_noTPMfilter']:.2f} vs baseline {base:.2f})")
        check("NOTAGENE" not in set(cf["measuredGeneSymbol"]) and cf["E2G.Score.qnorm"].notna().all() and
              (cf.loc[cf["name"] == "nowhere", "E2G.Score.qnorm"] == 0).all(), "CRISPR merge: TSS universe, NAfilled, unmatched pair -> fill 0")

        # --- the same stages one by one -----------------------------------------------------------------
        ab = S["abc_dir"]
        o1 = tmp / "o_steps"
        rc1 = run(["kendall-pairs", "--narrowpeak", str(ab / "Peaks" / "macs2_peaks.narrowPeak.sorted"), "--abc-predictions",
                   str(ab / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz"), "--cluster", "s", "--out-dir", str(o1)])
        pp = o1 / "s" / "Kendall" / "Pairs.tsv.gz"
        rc2 = run(["atac-matrix", "--pairs", str(pp), "--fragments", str(S["fragments"]), "--rna-matrix", str(S["rna"]), "--cluster", "s",
                   "--out-dir", str(o1)])
        rc3 = run(["kendall", "--pairs", str(pp), "--rna-matrix", str(S["rna"]), "--atac-matrix", str(o1 / "s" / "Kendall" / "atac_matrix"),
                   "--gtf", str(S["gtf"]), "--cluster", "s", "--out-dir", str(o1)] + common_ref)
        k1 = read_tsv(o1 / "s" / "Kendall" / "Pairs.Kendall.tsv.gz")
        check(rc1 == rc2 == rc3 == 0 and k1.equals(kd), "kendall-pairs + atac-matrix + kendall == run")
        rc4 = run(["arc", "--abc-predictions", str(ab / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz"), "--kendall",
                   str(o1 / "s" / "Kendall" / "Pairs.Kendall.tsv.gz"), "--cluster", "s", "--out-dir", str(o1)])
        a1 = read_tsv(o1 / "s" / "ARC" / "EnhancerPredictionsAllPutative_ARC.tsv.gz")
        check(rc4 == 0 and np.allclose(a1["ARC.E2G.Score"], arc["ARC.E2G.Score"]), "arc == run")
        mods = ["--models", "multiome_powerlaw_v3", "scATAC_powerlaw_v3"]
        rc5 = run(["activity-features", "--abc-predictions", str(ab / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz"),
                   "--enhancer-list", str(ab / "Neighborhoods" / "EnhancerList.txt"), "--gene-classes", str(S["gene_classes"]),
                   "--cluster", "s", "--out-dir", str(o1)] + mods + common_ref)
        rc6 = run(["features", "--activity-features", str(o1 / "s" / "ActivityOnly_features.tsv.gz"), "--arc",
                   str(o1 / "s" / "ARC" / "EnhancerPredictionsAllPutative_ARC.tsv.gz"), "--cluster", "s", "--out-dir", str(o1)] + mods)
        g1 = read_tsv(o1 / "s" / "genomewide_features.tsv.gz")
        check(rc5 == rc6 == 0 and g1[want].equals(gw[want]), "activity-features + features == run")
        rc7 = run(["predict", "--features", str(o1 / "s" / "genomewide_features.tsv.gz"), "--cluster", "s", "--bedpe", "--out-dir",
                   str(o1)] + mods + md)
        p1 = read_tsv(o1 / "s" / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz")
        check(rc7 == 0 and np.allclose(p1["E2G.Score.qnorm"], pred["E2G.Score.qnorm"]) and "E2G.Score.cv" not in p1.columns,
              "predict == run (no cv without --crispr-benchmarking)")
        rc8 = run(["predict", "--features", str(o1 / "s" / "genomewide_features.tsv.gz"), "--cluster", "nq", "--no-qnorm",
                   "--models", str(tmp / "exported" / "multiome_powerlaw_v3"), "--out-dir", str(o1)])
        p2 = read_tsv(o1 / "nq" / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz")
        check(rc8 == 0 and np.allclose(p2["E2G.Score"], pred["E2G.Score"]), "predict with an exported model directory")
        rc9 = run(["qnorm-ref", "--features", str(o1 / "s" / "genomewide_features.tsv.gz"), "--model", "multiome_powerlaw_v3",
                   "--n-scores", "11", "--write-predictions", "--out-dir", str(o1 / "qn")] + md)
        qr = read_tsv(o1 / "qn" / "multiome_powerlaw_v3" / "qnorm_reference.tsv.gz")
        check(rc9 == 0 and len(qr) == 11 and abs(qr["reference_score"].iloc[-1] - pred["E2G.Score"].max()) < 1e-12 and
              qr["reference_score"].is_monotonic_increasing, "qnorm-ref: quantiles of genome-wide E2G.Score")
        rc10 = run(["gene-lists", "--gene-list", str(ab / "Neighborhoods" / "GeneList.txt"), "--enhancer-list",
                    str(ab / "Neighborhoods" / "EnhancerList.txt"), "--predictions", str(cd / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz"),
                    "--gene-expression", str(cd / "Kendall" / "gene_expression_metrics.tsv.gz"), "--model", "multiome_powerlaw_v3",
                    "--cluster", "s", "--out-dir", str(o1)])
        g2 = read_tsv(o1 / "s" / "multiome_powerlaw_v3" / "scE2G_gene_list.tsv.gz")
        check(rc10 == 0 and g2.drop(columns="CellType").equals(gl.drop(columns="CellType")), "gene-lists == run")
        rc11 = run(["stats", "--predictions", str(cd / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz"), "--model-name",
                    "multiome_powerlaw_v3", "--fragment-count", str(cd / "fragment_count.txt"), "--cell-count", str(cd / "cell_count.txt"),
                    "--umi-count", str(cd / "umi_count.txt"), "--cluster", "s", "--out-dir", str(o1)])
        s2 = read_tsv(o1 / "s" / "multiome_powerlaw_v3" / "scE2G_predictions_threshold0.177_stats.tsv")
        check(rc11 == 0 and s2.drop(columns="cluster").equals(st.drop(columns="cluster")), "stats == run")
        rc12 = run(["qc", "--stats", str(o1 / "s" / "multiome_powerlaw_v3" / "scE2G_predictions_threshold0.177_stats.tsv"),
                    "--reference", str(S["qc_ref"]), "--out-dir", str(o1 / "qc")] + noplot)
        check(rc12 == 0 and (o1 / "qc" / "qc_plots" / "all_qc_stats.tsv").exists(), "qc subcommand")
        rc13 = run(["crispr-features", "--predictions", str(cd / "multiome_powerlaw_v3" / "scE2G_predictions.tsv.gz"), "--crispr",
                    str(S["crispr"]), "--tss", str(S["tss"]), "--cluster", "s", "--model-name", "multiome_powerlaw_v3",
                    "--out-dir", str(o1)] + mods)
        c2_ = read_tsv(o1 / "s" / "multiome_powerlaw_v3" / "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz")
        check(rc13 == 0 and c2_.equals(cf), "crispr-features == run")
        rc14 = run(["benchmark", "--crispr-features", str(o1 / "s" / "multiome_powerlaw_v3" /
                    "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"), "--n-boot", "30", "--ci-method", "BCa",
                    "--out-dir", str(o1 / "bench")] + noplot)
        b2 = read_tsv(o1 / "bench" / "crispr_benchmarking_performance_summary.tsv")
        r2 = b2.set_index("model").loc["multiome_powerlaw_v3"]
        check(rc14 == 0 and r2["AUPRC_95CI_low"] <= r2["AUPRC"] + 1e-9 <= r2["AUPRC_95CI_high"] + 2e-9 and
              abs(r2["AUPRC"] - statistic_aupr(cf["Regulated"].map(_truthy).astype(int), cf["E2G.Score.cv.qnorm"])) < 0.1,
              "benchmark: BCa bootstrap CI brackets the AUPRC")
        # training on the synthetic CRISPR truth (feature overlap in training mode)
        rc15 = run(["train", "--features", str(cd / "genomewide_features.tsv.gz"), "--crispr", str(S["crispr"]), "--tss", str(S["tss"]),
                    "--model", "multiome_powerlaw_v3", "--name", "synth_model", "--override-params", '{"max_iter": 2000}',
                    "--tpm-threshold", "1", "--out-dir", str(o1 / "train")])
        tm = resolve_model(str(o1 / "train" / "synth_model"), S["models_dir"])
        met = read_tsv(o1 / "train" / "synth_model" / "model" / "performance_metrics.tsv")
        coef = read_tsv(o1 / "train" / "synth_model" / "model" / "model_coefficients.tsv")
        w_arc = coef[(coef["test_chr"] == "none") & (coef["feature"] == "ARC.E2G.Score")]["coefficient"].iloc[0]
        check(rc15 == 0 and len(tm["cv"]) == 2 and w_arc > 0 and set(met["test_chr"]) == {"chr1", "chr2", "all"} and
              tm["score_threshold"] is not None and tm["tpm_threshold"] == 1.0,
              "train: full + per-chromosome models, ARC weight > 0, thresholds written")
        rc16 = run(["predict", "--features", str(cd / "genomewide_features.tsv.gz"), "--models", str(o1 / "train" / "synth_model"),
                    "--no-qnorm", "--cluster", "t", "--out-dir", str(o1)])
        check(rc16 == 0 and (o1 / "t" / "synth_model" / "scE2G_predictions.tsv.gz").exists(), "a trained model directory predicts")
        # run without ABC outputs -> prints the ABC command, exit 0
        rc17 = run(["run", "--cluster", "noabc", "--fragments", str(S["fragments"]), "--rna-matrix", str(S["rna"]),
                    "--out-dir", str(tmp / "o_noabc")] + common_ref + md + noplot)
        check(rc17 == 0 and (tmp / "o_noabc" / "tmp" / "config_abc_biosamples.tsv").exists(), "run without ABC outputs: tagAlign + ABC command, exit 0")
        # a default run directory under Docs/, then removed
        before = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
        rc18 = run(["models", "--label", "selftest_models"])
        new = (set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()) - before
        made += list(new)
        check(rc18 == 0 and len(new) == 1 and (list(new)[0] / "report.md").exists(), "default run directory under Docs/scE2GPipeline")
        rc19 = run(["setup", "--dry-run", "--crispr", "--gtf", "--example"])
        check(rc19 == 0, "setup --dry-run")
        check(time.time() - t0 < 110, f"selftest time {time.time() - t0:.1f}s")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        for p in made:
            shutil.rmtree(p, ignore_errors=True)
        if not out_root_existed and OUT_ROOT.exists() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _common(p, label: str, cluster: bool = True) -> None:
    p.add_argument("--label", default=label, help="run label (output directory suffix)")
    p.add_argument("--out-dir", help="write here instead of a new Docs/scE2GPipeline/<timestamp>_<label> directory")
    if cluster:
        p.add_argument("--cluster", default="cluster", help="cluster / biosample name (output subdirectory, CellType)")


def _model_args(p) -> None:
    p.add_argument("--models", nargs="+", help="embedded model names (multiome_powerlaw_v3, multiome_megamap_v3, "
                   "scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories (default: multiome_powerlaw_v3 scATAC_powerlaw_v3)")
    p.add_argument("--models-dir", help=f"where qnorm references live (default {MODELS_DIR})")


def _ref_args(p, gtf: bool = False) -> None:
    p.add_argument("--chrom-sizes", help="chromosome sizes (default: setup's GRCh38_EBV.no_alt.chrom.sizes.tsv)")
    p.add_argument("--tss", help="CollapsedGeneBounds TSS500bp bed (default: setup's intGENCODEv43 file)")
    if gtf:
        p.add_argument("--gtf", help="GENCODE GTF used for the RNA gene names (default: setup --gtf)")
        p.add_argument("--no-gene-mapping", action="store_true", help="RNA row names already are TSS500 gene names")


def _cell_args(p) -> None:
    p.add_argument("--rna-unfiltered", action="store_true", help="RNA matrix has extra cells: keep those seen in the fragments")
    p.add_argument("--cell-barcodes", help="cell_barcodes.txt from frag-to-tagalign (with --rna-unfiltered)")
    p.add_argument("--max-cell-count", type=int, default=MAX_CELL_COUNT)
    p.add_argument("--count-mode", choices=["fragment", "insertion"], default="fragment",
                   help="ATAC matrix: fragments overlapping a peak (Signac, default) or Tn5 insertions in it")


def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent sce2g-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup", help="fetch qnorm references and reference files from the pinned upstream commits")
    p.add_argument("--crispr", action="store_true", help="also the CRISPR benchmark (EPCrisprBenchmark ... intGENCODEv43)")
    p.add_argument("--gtf", action="store_true", help="also the GENCODE v43 GTF (53 MB) for RNA gene-name mapping")
    p.add_argument("--example", action="store_true", help="also the upstream chr22 test fixture (for validate-example)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="list what would be fetched")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("models", help="the embedded v3 models (weights, thresholds); --export writes model directories")
    _common(p, "models", cluster=False)
    p.add_argument("--export", help="write upstream-shaped model directories here")
    p.add_argument("--models-dir")
    p.set_defaults(func=cmd_models)

    p = sub.add_parser("frag-to-tagalign", help="fragments -> filtered fragments, counts, barcodes, tagAlign")
    _common(p, "tagalign")
    p.add_argument("--fragments", required=True)
    p.add_argument("--chrom-sizes")
    p.add_argument("--preprocessed", action="store_true", help="fragments already sorted and restricted to chrom-sizes chromosomes")
    p.set_defaults(func=cmd_frag_to_tagalign)

    p = sub.add_parser("bigwig", help="fragments -> ATAC coverage bedGraph / bigWig (raw with --raw, per-million normalised)")
    _common(p, "bigwig")
    p.add_argument("--fragments", required=True)
    p.add_argument("--chrom-sizes")
    p.add_argument("--fragment-count", help="fragment_count.txt or a number (default: fragments on chrom-sizes chromosomes)")
    p.add_argument("--raw", action="store_true", help="also the unnormalised ATAC track")
    p.set_defaults(func=cmd_bigwig)

    p = sub.add_parser("kendall-pairs", help="narrowPeak + ABC AllPutative -> Kendall/Pairs.tsv.gz")
    _common(p, "kendall_pairs")
    p.add_argument("--narrowpeak", required=True, help="ABC Peaks/macs2_peaks.narrowPeak.sorted")
    p.add_argument("--abc-predictions", required=True, help="ABC Predictions/EnhancerPredictionsAllPutative.tsv.gz")
    p.set_defaults(func=cmd_kendall_pairs)

    p = sub.add_parser("atac-matrix", help="pairs + fragments + RNA cells -> peak x cell count matrix")
    _common(p, "atac_matrix")
    p.add_argument("--pairs", required=True)
    p.add_argument("--fragments", required=True)
    p.add_argument("--rna-matrix", required=True, help="RNA counts (.h5ad, .csv.gz genes x cells, or 10x directory) for the cell list")
    _cell_args(p)
    p.set_defaults(func=cmd_atac_matrix)

    p = sub.add_parser("kendall", help="pairs + ATAC + RNA -> Pairs.Kendall.tsv.gz and gene expression metrics")
    _common(p, "kendall")
    p.add_argument("--pairs", required=True)
    p.add_argument("--rna-matrix", required=True)
    p.add_argument("--atac-matrix", help="prefix/path written by atac-matrix (else built from --fragments)")
    p.add_argument("--fragments")
    _cell_args(p)
    _ref_args(p, gtf=True)
    p.set_defaults(func=cmd_kendall)

    p = sub.add_parser("arc", help="ABC AllPutative + Pairs.Kendall -> ARC-E2G")
    _common(p, "arc")
    p.add_argument("--abc-predictions", required=True)
    p.add_argument("--kendall", required=True, help="Pairs.Kendall.tsv.gz")
    p.add_argument("--abc-score-col", default="powerlaw.Score")
    p.add_argument("--hic", action="store_true", help="cluster has Hi-C: use ABC.Score")
    p.set_defaults(func=cmd_arc)

    p = sub.add_parser("activity-features", help="ENCODE_rE2G numCandidateEnhGene / numTSSEnhGene / nearby enhancers")
    _common(p, "activity_features")
    p.add_argument("--abc-predictions", required=True)
    p.add_argument("--enhancer-list", help="ABC Neighborhoods/EnhancerList.txt")
    p.add_argument("--gene-classes", help="gene_promoter_class ... TSS500bp.tsv (default: setup)")
    _model_args(p)
    _ref_args(p)
    p.set_defaults(func=cmd_activity_features)

    p = sub.add_parser("features", help="activity features + ARC -> genome-wide feature table")
    _common(p, "features")
    p.add_argument("--activity-features", required=True, help="ActivityOnly_features.tsv.gz")
    p.add_argument("--arc", help="EnhancerPredictionsAllPutative_ARC.tsv.gz")
    p.add_argument("--external-config", help="extra external_features_config.tsv rows")
    _model_args(p)
    p.set_defaults(func=cmd_features)

    p = sub.add_parser("predict", help="genome-wide features -> scE2G predictions, thresholded links, bedpe")
    _common(p, "predict")
    p.add_argument("--features", required=True, help="genomewide_features.tsv.gz")
    _model_args(p)
    p.add_argument("--crispr-benchmarking", action="store_true", help="also E2G.Score.cv from the held-out-chromosome models")
    p.add_argument("--no-qnorm", action="store_true", help="E2G.Score.qnorm = E2G.Score (no reference needed)")
    p.add_argument("--bedpe", action="store_true")
    p.add_argument("--exclude-self-promoter", action="store_true", help="drop self-promoter links (upstream default keeps them)")
    p.set_defaults(func=cmd_predict)

    p = sub.add_parser("qnorm-ref", help="genome-wide features -> a model's qnorm_reference.tsv.gz")
    _common(p, "qnorm_ref", cluster=False)
    p.add_argument("--features", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--models-dir")
    p.add_argument("--n-scores", type=int, default=10000)
    p.add_argument("--tpm-threshold", type=float)
    p.add_argument("--write-predictions", action="store_true", help="also genomewide_predictions.tsv.gz (score, cv, TPM filter)")
    p.set_defaults(func=cmd_qnorm_ref)

    p = sub.add_parser("gene-lists", help="ABC gene/element lists + predictions -> scE2G gene and element lists")
    _common(p, "gene_lists")
    p.add_argument("--gene-list", required=True)
    p.add_argument("--enhancer-list", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--gene-expression", help="Kendall/gene_expression_metrics.tsv.gz")
    p.add_argument("--model")
    p.add_argument("--tpm-threshold", type=float)
    p.set_defaults(func=cmd_gene_lists)

    p = sub.add_parser("stats", help="predictions -> per-cluster QC statistics")
    _common(p, "stats")
    p.add_argument("--predictions", required=True)
    p.add_argument("--thresholded")
    p.add_argument("--threshold", type=float)
    p.add_argument("--model-name", required=True)
    p.add_argument("--fragment-count")
    p.add_argument("--cell-count")
    p.add_argument("--umi-count")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("qc", help="stats files -> all_qc_stats, warnings, reference comparison, figures")
    _common(p, "qc", cluster=False)
    p.add_argument("--stats", nargs="+", required=True)
    p.add_argument("--reference", help="reference_qc_metrics_sheth_qiu_2024.tsv (default: setup)")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("crispr-features", help="predictions/features + CRISPR benchmark -> CRISPR feature table")
    _common(p, "crispr_features")
    p.add_argument("--predictions", required=True, help="scE2G_predictions.tsv.gz (apply) or genomewide_features.tsv.gz (training)")
    p.add_argument("--crispr")
    p.add_argument("--tss")
    p.add_argument("--feature-table")
    p.add_argument("--model-name")
    p.add_argument("--mode", choices=["apply", "training"], default="apply")
    p.add_argument("--no-nafill", action="store_true")
    _model_args(p)
    p.set_defaults(func=cmd_crispr_features)

    p = sub.add_parser("benchmark", help="CRISPR feature tables -> AUPRC / precision / recall with bootstrap CIs")
    _common(p, "benchmark", cluster=False)
    p.add_argument("--crispr-features", nargs="+", required=True, help=".../<cluster>/<model>/EPCrispr..._NAfilled.tsv.gz")
    p.add_argument("--model-names", nargs="+")
    p.add_argument("--model-thresholds", nargs="+")
    p.add_argument("--clusters", nargs="+")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--ci-method", default="BCa", choices=["BCa", "percentile", "basic"])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_benchmark)

    p = sub.add_parser("train", help="CRISPR features -> a new scE2G model directory")
    _common(p, "train", cluster=False)
    p.add_argument("--crispr-features", help="CRISPR x features table (training mode, NAfilled)")
    p.add_argument("--features", help="genomewide_features.tsv.gz (merged with --crispr here)")
    p.add_argument("--crispr")
    p.add_argument("--tss")
    p.add_argument("--model", default="multiome_powerlaw_v3", help="embedded model whose feature table to use")
    p.add_argument("--feature-table", help="custom feature_table.tsv")
    p.add_argument("--name", default="custom_model")
    p.add_argument("--polynomial", action="store_true")
    p.add_argument("--override-params", help="JSON of LogisticRegression parameters overriding the defaults")
    p.add_argument("--score-threshold", type=float, help="default: threshold at 70%% recall of the cv score")
    p.add_argument("--tpm-threshold", type=float, default=0.0)
    p.add_argument("--pickle", action="store_true", help="also model.pkl (scikit-learn)")
    p.set_defaults(func=cmd_train)

    p = sub.add_parser("run", help="the per-cluster pipeline from fragments + RNA + ABC outputs")
    _common(p, "run")
    p.add_argument("--cluster-config", help="upstream config_cell_clusters.tsv (+ abc_dir column)")
    p.add_argument("--fragments")
    p.add_argument("--rna-matrix")
    p.add_argument("--abc-dir", help="ABC results for the cluster (Peaks/, Predictions/, Neighborhoods/)")
    for k in ABC_LAYOUT:
        p.add_argument("--" + k.replace("_", "-"), help=f"override {ABC_LAYOUT[k]}")
    p.add_argument("--hic", action="store_true")
    p.add_argument("--fragments-preprocessed", action="store_true")
    p.add_argument("--gene-classes")
    p.add_argument("--qc-reference")
    p.add_argument("--crispr", help="CRISPR benchmark table: also cv scores, CRISPR features and the benchmark")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--ci-method", default="BCa", choices=["BCa", "percentile", "basic"])
    p.add_argument("--igv-tracks", action="store_true", help="ATAC_norm track + bedpe links")
    p.add_argument("--no-qnorm", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    _model_args(p)
    _ref_args(p, gtf=True)
    _cell_args(p)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("validate-example", help="Kendall / RNA features vs the upstream chr22 test fixture")
    _common(p, "validate_example", cluster=False)
    p.add_argument("--example-dir")
    p.add_argument("--gtf")
    p.add_argument("--tss")
    p.add_argument("--count-mode", choices=["fragment", "insertion"], default="fragment")
    p.set_defaults(func=cmd_validate_example)

    p = sub.add_parser("validate-release", help="recompute a released IGVF scE2G prediction set from its own files")
    _common(p, "validate_release", cluster=False)
    p.add_argument("--pairs", required=True, help="all element-gene pairs file (e.g. IGVFFI1706PNVV.tsv.gz)")
    p.add_argument("--elements", required=True, help="element file (e.g. IGVFFI3094BAXH.tsv.gz)")
    p.add_argument("--genes", required=True, help="gene file (e.g. IGVFFI9905RPTO.tsv.gz)")
    p.add_argument("--model", default="multiome_powerlaw_v3")
    p.add_argument("--models-dir")
    p.add_argument("--chrom", help="restrict to one chromosome (qnorm is then not comparable)")
    p.add_argument("--write-predictions", action="store_true", help="also write the recomputed table (large)")
    _ref_args(p)
    p.set_defaults(func=cmd_validate_release)

    p = sub.add_parser("selftest", help="synthetic clusters with planted enhancers; every subcommand asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
