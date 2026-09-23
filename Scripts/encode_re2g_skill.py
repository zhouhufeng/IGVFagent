#!/usr/bin/env python3
"""ENCODE-rE2G enhancer-gene prediction: features, pretrained models, training and feature analysis (port of EngreitzLab/ENCODE_rE2G).

Port of https://github.com/EngreitzLab/ENCODE_rE2G (MIT, Engreitz Lab 2025;
Snakemake + Python + R + bedtools), the logistic-regression enhancer-to-gene
model of Gschwind et al. 2026 ("An encyclopedia of human enhancer-gene
regulatory interactions") built on top of ABC outputs.  Every Snakemake rule,
Python script and R script of both workflows (apply: workflow/Snakefile;
training: workflow/Snakefile_training) was read and re-derived in Python
(pandas + numpy; scikit-learn / scipy when present).  No code was copied.
Relationship: port.  bedtools / csvtk / pigz / GenomicRanges are replaced by
vectorised interval arithmetic; samtools is optional (read counts only).

The nine pretrained models under upstream models/ (feature_table.tsv,
threshold_<t>, model.pkl = a 1 KB sklearn LogisticRegression) are embedded
below as plain coefficients (EMBEDDED_MODELS, extracted from the pickles at
commit d039062b7092de338ed10a77ea208b3d5b06e89c), so scoring needs neither
the pickles nor pickle compatibility.  `verify-upstream` re-scores the
upstream's own expected test output (K562 chr22, dhs_intact_hic) with them.

Definitions, exactly as upstream computes them
  model choice       <accessibility>[_h3k27ac]_<contact>: accessibility =
                     default_accessibility_feature lower-cased (dhs|atac);
                     "_h3k27ac" when the H3K27ac column is filled; contact =
                     powerlaw (no HiC_file) | avg_hic (HiC_type avg) |
                     megamap (HiC_file basename == tissues.hic) | intact_hic
                     (HiC_type hic).  A model_dir column (comma-separated)
                     overrides; a model dir holds feature_table.tsv, model.pkl
                     and exactly one threshold_<t> file.  powerlaw / avg_hic
                     models were never trained upstream -> error.
  midpoint           int((start + end) / 2).
  numCandidateEnhGene  non-promoter pairs of a (TargetGene, TSS) sorted by
                     midpoint; 1 = the candidate closest to the TSS on its
                     side, counting outwards; midpoint == TSS -> 0.
  numTSSEnhGene      number of TSS500bp reference intervals overlapping
                     [start, midpoint + distance) (TSS downstream of E) or
                     [TSS, end) (TSS upstream of E) -- target gene's own TSS
                     included; no overlap -> NA -> fill 0.
  numNearbyEnhancers / sumNearbyEnhancers  window midpoint +/- 5 kb (clipped
                     to [0, chromosome size]) intersected with every row of
                     EnhancerPredictionsAllPutative (E-G rows, promoters
                     included), excluding rows with the same name: row count
                     and sum of activity_base.  Upstream's "| sort -u" is a
                     no-op (it follows the redirect), so an element linked to
                     k genes counts k times; the pretrained models learned
                     that, so it is the default (--dedupe-nearby applies the
                     intended de-duplication).
  ABC.Numerator      ABC.Score.Numerator;  ABC.Denominator = ABC.Score /
                     ABC.Numerator (upstream's definition, i.e. 1/denominator).
  external features  join_by overlap: same chr:TargetGene and closed-interval
                     overlap (GenomicRanges IRanges semantics, touching counts),
                     several hits aggregated by aggregate_function (R NA
                     semantics: any NA -> NA), unmatched -> NA;  join_by
                     TargetGene: left join on the gene.
  final features     feature = input_col * second_input when second_input is
                     set, otherwise input_col renamed to feature; NA filled
                     with fill_value, "mean" = mean over finite values.
  score              X = log(|x| + epsilon), epsilon 0.01, NaN/Inf -> 0 first;
                     ENCODE-rE2G.Score = predict_proba(X)[:, 1] =
                     1 / (1 + exp(-(X.coef + intercept))).
  threshold          score >= t and (class != promoter or isSelfPromoter when
                     include_self_promoter, default True).
  bedpe              TargetGeneIsExpressed rows: chr start end chr TSS TSS
                     TargetGene_name score . .
  stats              num_sequencing_reads (samtools idxstats, chr1-22,X,Y,
                     BAM only), num_enh, num_genes, num_enh_gene_links,
                     num_genes_with_1_enh_min, mean_num_genes_per_enh,
                     mean_num_enh_per_gene(_no_prom), mean_log10_dist_to_tss
                     (log10(0) -> 0), mean_enh_region_size.
  CRISPR overlap     chrom:measuredGeneSymbol + closed-interval overlap with
                     the genome-wide features, aggregate_function per feature;
                     genes outside the TSS universe removed first; unmatched
                     pairs -> missing file; NA distanceToTSS recomputed from the
                     element centre and the TSS reference (or the CRISPR TSS);
                     NA filled with fill values computed on the GENOME-WIDE
                     table ("NAfilled").  process: rename chrom/chromStart/
                     chromEnd/measuredGeneSymbol -> chr/start/end/TargetGene,
                     drop Regulated NA, keep TargetGene in the gene universe.
  training           LogisticRegression(solver lbfgs, penalty None, max_iter
                     1e8, tol 1e-4, random_state 0) on log(|x| + 0.01) of the
                     feature table (optionally degree-2 PolynomialFeatures of the
                     raw features first); full model + leave-one-chromosome-out
                     CV scores; per-chromosome and pooled log loss / AUROC /
                     AUPRC.
  AUPRC              sklearn precision_recall_curve with the first point
                     (recall 1, precision = class balance) dropped, trapezoid
                     auc; < 2 points -> 0.
  precision at 70% recall  threshold = thresholds[i + 1] with i the point whose
                     recall is closest to 0.7 (None if max recall <= 0.7);
                     precision = precision[j - 1] at the threshold closest to it.
  bootstrap          scipy.stats.bootstrap, paired, BCa, 95%, n = 1000;
                     p = mean(|dist - mean(dist)| >= |delta|).
  feature selection  forward: greedily add the feature maximising CV AUPRC
                     (ties -> the later one); backward: greedily drop the
                     feature whose removal keeps AUPRC highest; then refit the
                     sequence with bootstrapped delta AUPRC / delta precision.
  permutation importance  n_repeats (20) shuffles per feature, CV refit,
                     delta AUPRC and delta precision vs the full model.
  all feature sets   every non-empty subset (2^n - 1; upstream runs it when
                     n < 14 and no polynomial), bootstrapped AUPRC, sorted.
  model comparison   CV scores of each trained model; CRISPR pairs missing from
                     the features count as score 0; distance baseline =
                     -|element centre - TSS centre|.

One upstream bug is corrected by default and reproduced with --upstream-compat:
  train_model.py writes log_loss() into AUROC_test_full of the pooled "all"
  row.  Other deliberate differences: merge_external_features.R tests
  row.names(abc) for already-present features (always empty, so everything is
  re-merged); the port tests column names.  Random steps (bootstrap,
  permutation) take --seed; sklearn's verbose flag is forced off.

Subcommands
  setup            fetch reference files (TSS500bp universe, chrom sizes, gene
                   classes, CRISPR benchmark) from the pinned commit; --test-data
                   also fetches the chr22 expected output (~47 MB).
  models           list / export the embedded pretrained models (model dirs
                   with feature_table.tsv, threshold_<t>, model.json).
  select-model     biosample config -> model name(s), dir(s) and threshold.
  features         ABC outputs -> NumCandidateEnhGene / NumTSSEnhGene /
                   NumEnhancersEG5kb / SumEnhancersEG5kb, ActivityOnly features,
                   external features, genomewide_features.tsv.gz.
  apply            genome-wide features x model -> encode_e2g_predictions.tsv.gz,
                   thresholded table, bedpe and stats.
  run              features + apply for every biosample of a config.
  threshold        threshold_e2g_predictions.
  bedpe            process_model_output (IGV bedpe).
  stats            get_stats for a thresholded prediction file.
  qc-plots         generate_plots over many stats files.
  compare-stats    compare_plots: one metric across two sets of runs.
  crispr-features  overlap_features_with_crispr_data + process_crispr_data.
  train            get_params + train_model (full + leave-one-chromosome-out).
  feature-selection  forward / backward sequential selection with bootstrap.
  permutation-importance  permutation_feature_importance.
  all-feature-sets compare_all_feature_sets.
  compare-models   compare_all_models (+ distance baseline) and plots.
  verify-upstream  re-score the upstream expected output with the embedded model.
  selftest         synthetic ABC world with a planted model, every subcommand.

Output: Docs/ENCODErE2G/<timestamp>_<label>/.

Usage:
    igvfagent encode-re2g setup
    igvfagent encode-re2g models --export models_dir
    igvfagent encode-re2g select-model --biosample-config config_biosamples.tsv
    igvfagent encode-re2g features --abc-dir results/K562 --model dhs_intact_hic --tss TSS500bp.bed --chr-sizes sizes.tsv --gene-classes gene_classes.tsv --label K562
    igvfagent encode-re2g apply --features genomewide_features.tsv.gz --model dhs_intact_hic --label K562
    igvfagent encode-re2g run --biosample-config config_biosamples.tsv --abc-results results --label encode
    igvfagent encode-re2g crispr-features --features genomewide_features.tsv.gz --crispr EPCrisprBenchmark.tsv.gz --model dhs_intact_hic --label K562
    igvfagent encode-re2g train --crispr-features for_training.tsv.gz --model dhs_intact_hic --label my_model
    igvfagent encode-re2g feature-selection --crispr-features for_training.tsv.gz --feature-table ft.tsv --direction forward
    igvfagent encode-re2g compare-models --train-dirs runA runB --crispr EPCrisprBenchmark.tsv.gz
    igvfagent encode-re2g verify-upstream --upstream-dir ENCODE_rE2G
    igvfagent encode-re2g selftest --no-plots
"""
from __future__ import annotations

import argparse
import glob
import gzip
import itertools
import json
import logging
import math
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import warnings
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "ENCODErE2G"
DATA_ROOT = ROOT / "Data" / "ENCODErE2G"
RES_DIR = DATA_ROOT / "resources"

UPSTREAM_REPO = "EngreitzLab/ENCODE_rE2G"
UPSTREAM_COMMIT = "d039062b7092de338ed10a77ea208b3d5b06e89c"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"
RESOURCE_FILES = {
    "RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed": "reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed",
    "GRCh38_EBV.no_alt.chrom.sizes.tsv": "reference/GRCh38_EBV.no_alt.chrom.sizes.tsv",
    "gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv":
        "resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv",
    "EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz": "reference/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz",
}
TEST_FILES = {  # upstream CircleCI expected output, K562 chr22, dhs_intact_hic (~47 MB)
    "encode_e2g_predictions.tsv.gz": "tests/expected_output/generic/K562_chr22/dhs_intact_hic/encode_e2g_predictions.tsv.gz",
    "encode_e2g_predictions_threshold0.243.tsv.gz": "tests/expected_output/generic/K562_chr22/dhs_intact_hic/encode_e2g_predictions_threshold0.243.tsv.gz",
    "encode_e2g_predictions_threshold0.243_stats.tsv": "tests/expected_output/generic/K562_chr22/dhs_intact_hic/encode_e2g_predictions_threshold0.243_stats.tsv",
}
MEGAMAP_HIC_FILE = "https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic"
EPSILON = 0.01
MODEL_NAME = "ENCODE-rE2G"
SCORE_COL = MODEL_NAME + ".Score"
NEARBY_WINDOW = 5000
NORMAL_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
FT_COLS = ["feature", "input_col", "second_input", "aggregate_function", "fill_value", "nice_name"]
CORE_COLS = ["chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS", "TargetGeneIsExpressed",
             "TargetGeneEnsembl_ID", "isSelfPromoter", "CellType", "distance"]
REFERENCE_FEATURES = ["numTSSEnhGene", "distance", "activity_base", "TargetGenePromoterActivityQuantile", "numNearbyEnhancers",
                      "sumNearbyEnhancers", "is_ubiquitous_uniform", "P2PromoterClass", "numCandidateEnhGene",
                      "hic_contact_pl_scaled_adj", "ABC.Score", "ABC.Numerator", "ABC.Denominator", "normalized_dhs_prom",
                      "normalized_dhs_enh", "normalized_atac_prom", "normalized_atac_enh", "normalized_h3k27ac_enh",
                      "normalized_h3k27ac_prom"]
# config/config_training.yaml default_params (verbose forced to 0 here; it only prints lbfgs iterations)
DEFAULT_PARAMS = {"solver": "lbfgs", "fit_intercept": True, "penalty": None, "max_iter": 100000000, "class_weight": None,
                  "tol": 1e-4, "warm_start": False, "random_state": 0, "verbose": 0, "n_jobs": 1}
# combine_feature_tables(.R): rows added for sc-E2G feature tables (ARC.E2G.Score or Kendall present)
ARC_ROWS = [("RNA_meanLogNorm", "mean_log_normalized_rna", "NA", "mean", "0", "Mean log normalized RNA expression"),
            ("RNA_pseudobulkTPM", "RnaPseudobulkTPM", "NA", "mean", "0", "RNA pseudobulk TPM"),
            ("RNA_percentCellsDetected", "RnaDetectedPercent", "NA", "mean", "0", "RNA percent cells detected"),
            ("Kendall", "Kendall", "NA", "max", "0", "Kendall correlation"),
            ("ARC.E2G.Score", "ARC.E2G.Score", "NA", "mean", "0", "ARC-E2G score")]
STATS_METRICS = ["num_enh", "num_genes", "num_enh_gene_links", "num_genes_with_1_enh_min", "mean_num_genes_per_enh",
                 "mean_num_enh_per_gene", "mean_num_enh_per_gene_no_prom", "mean_log10_dist_to_tss", "mean_enh_region_size"]
SEQ_DEPTH_METRIC = "num_sequencing_reads"
METRIC_COLS = ["test_chr", "log_loss_test_full", "log_loss_train", "log_loss_test", "AUROC_test_full", "AUROC_train",
               "AUROC_test", "AUPRC_test_full", "AUPRC_train", "AUPRC_test", "n_test_pos", "n_test_neg", "n_train_neg",
               "n_train_pos"]
SEL_COLS = ["aupr", "delta_aupr", "delta_aupr_low", "delta_aupr_high", "pval_aupr", "precision", "delta_precision",
            "delta_precision_low", "delta_precision_high", "pval_precision"]

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
GREY = "#96a0b3"

log = logging.getLogger("encode_re2g")

# Extracted from upstream models/<name>/model.pkl (sklearn LogisticRegression; coef_[0], intercept_[0]),
# models/<name>/feature_table.tsv and the models/<name>/threshold_<t> file name, commit d039062b7092de338ed10a77ea208b3d5b06e89c.
EMBEDDED_MODELS = {
    "atac_h3k27ac_intact_hic": {
        "threshold": 0.259, "threshold_file": "threshold_.259",
        "intercept": 6.148611470040296,
        "coef": [-0.8843764394424506, -0.2108099689178269, -3.1349459623599523, 0.06152711790736186, -0.2572511097045459, 0.7281200037171768, 0.8860255966069794, 1.7155276873119707],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('normalizedH3K27ac_prom', 'normalized_h3k27ac_prom', 'NA', 'mean', '0', 'H3K27ac signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "atac_h3k27ac_megamap": {
        "threshold": 0.203, "threshold_file": "threshold_.203",
        "intercept": 10.64052903271604,
        "coef": [-1.0349326405437151, -1.0667595198853996, -2.8335110651464563, 0.07135358727996773, -0.2760382991178266, 0.856907978770046, -0.8873790467340972, 1.9473838454429782],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('normalizedH3K27ac_prom', 'normalized_h3k27ac_prom', 'NA', 'mean', '0', 'H3K27ac signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "atac_intact_hic": {
        "threshold": 0.232, "threshold_file": "threshold_.232",
        "intercept": 5.656715132901357,
        "coef": [-1.070571341728566, -0.14061769017657152, -1.5167467222518993, 0.03265009559657677, -0.2900958025430249, 0.8274244749870326, 1.092870024539639, 1.494441391090369],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', '', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', '', 'min', '', 'Distance to TSS'),
            ('normalizedATAC_prom', 'normalized_atac_prom', '', 'mean', '0', 'ATAC signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', '', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', '', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', '', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', '', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', '', 'sum', '0', 'ABC score'),
        ],
    },
    "atac_megamap": {
        "threshold": 0.179, "threshold_file": "threshold_.179",
        "intercept": 10.16830611144203,
        "coef": [-1.1993012434930714, -0.9865371813647952, -1.1676428264426517, 0.043383475241728, -0.30365315422700256, 0.9356101107210496, -0.6153906451583226, 1.6825993517427564],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('normalizedATAC_prom', 'normalized_atac_prom', 'NA', 'mean', '0', 'ATAC signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "dhs_h3k27ac_intact_hic": {
        "threshold": 0.298, "threshold_file": "threshold_.298",
        "intercept": 5.951823520244407,
        "coef": [-0.885899895564728, -0.09343694598910249, -2.8476057076769297, 0.03856285780917579, -0.23735343628343067, 0.5039119126387802, 0.5179170699817907, 2.133405111774603],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('activity_prom', 'TargetGenePromoterActivityQuantile', 'NA', 'mean', '0', 'DNase x H3K27ac signal at P (quantile)'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "dhs_h3k27ac_megamap": {
        "threshold": 0.217, "threshold_file": "threshold_.217",
        "intercept": 10.030956536886615,
        "coef": [-1.013872818293418, -0.8298693520653417, -2.453865101173975, 0.03895525291861657, -0.25874065369547977, 0.5921491507767986, -1.158147329262213, 2.4886259369007795],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('activity_prom', 'TargetGenePromoterActivityQuantile', 'NA', 'mean', '0', 'DNase x H3K27ac signal at P (quantile)'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "dhs_intact_hic": {
        "threshold": 0.243, "threshold_file": "threshold_.243",
        "intercept": 5.618330970851473,
        "coef": [-0.956707370025905, -0.07003320677368341, -2.2190334504015223, 0.05063874839875221, -0.24586853191017227, 0.5640778652588982, 0.6248398288161648, 2.0560501356820873],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('normalizedDNase_prom', 'normalized_dhs_prom', 'NA', 'mean', '0', 'DNase signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "dhs_megamap": {
        "threshold": 0.201, "threshold_file": "threshold_.201",
        "intercept": 9.491874512962376,
        "coef": [-1.0904627388710866, -0.735419489396911, -1.809477128534176, 0.05296287699951543, -0.26535482462732857, 0.6501289094492338, -0.9301090484423687, 2.448962859897797],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'min', 'NA', 'Distance to TSS'),
            ('normalizedDNase_prom', 'normalized_dhs_prom', 'NA', 'mean', '0', 'DNase signal at P'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('contactFrequency', 'hic_contact_pl_scaled_adj', 'NA', 'mean', '0', 'Contact frequency'),
            ('ABC.Score', 'ABC.Score', 'NA', 'sum', '0', 'ABC score'),
        ],
    },
    "extended": {
        "threshold": 0.336, "threshold_file": "threshold_0.336",
        "intercept": 5.050679955432233,
        "coef": [-0.6008412196311012, -0.03127453147678739, 4.795323394681517, 2.6243802535622405, -1.9695532990398765, -0.9232468978520423, 6.802962870204077, -0.9096511591334105, 1.2691878472092797, 0.7504220197024996, 0.44178956410753806, -0.18038006423419192, 0.2135280885786818, -0.152594639413453, -0.13435318906773347, -0.25532456036889534, -0.10435776996370495, 0.033277830929179575, 0.058910969583015935, 0.05052219888685203, 0.15030103098675127, -0.152974198874394, -0.17200605288620974, 0.10654796019040284, -0.08525386968077096, 0.44247990935930753, -0.3273894115102607, -0.07892769523019695, -0.058146322545523904, -0.10899971322655855, 0.4764889916781314, 0.17606822987867138, -0.12688845312414998, -0.2545797179478666, -0.05668809959760874, -0.1131362016648394, 0.5005477200439075, 0.17524922962521117, 0.054984923130058566, -0.013161465244964171, 0.8038561827300648, 0.14453664556813728, -3.744184219866797, -0.01370929490366007, -4.519429805322874],
        "feature_table": [
            ('numTSSEnhGene', 'numTSSEnhGene', 'NA', 'max', '0', '# TSSs between E and P'),
            ('distanceToTSS', 'distance', 'NA', 'max', 'NA', 'Distance to TSS'),
            ('activity_enh', 'activity_base_enh', 'NA', 'mean', '0', 'sqrt(DNase x H3K27ac signal at E)'),
            ('3DContact', 'hic_contact_pl_scaled_adj', '', 'mean', '0', 'Contact frequency (ENCODE Hi-C)'),
            ('activity_enh_squared', 'activity_base_enh', 'activity_base_enh', 'mean', '0', 'DNase x H3K27ac signal at E'),
            ('3DContact_squared', 'hic_contact_pl_scaled_adj', 'hic_contact_pl_scaled_adj', 'mean', '0', '(Contact frequency (ENCODE Hi-C))^2'),
            ('activity_prom', 'TargetGenePromoterActivityQuantile', 'NA', 'mean', '0', 'DNase x H3K27ac signal at P (quantile)'),
            ('ABCNumerator', 'ABC.Numerator', 'NA', 'sum', '0', 'ABC_A=DNase-seq x H3K27ac, C=ENCODE Hi-C numerator'),
            ('ABCScore', 'ABC.Score', 'NA', 'sum', '0', 'ABC_A=DNase-seq x H3K27ac, C=ENCODE Hi-C score'),
            ('ABCDenominator', 'ABC.Denominator', 'NA', 'sum', 'mean', 'ABC_A=DNase-seq x H3K27ac, C=ENCODE Hi-C denominator'),
            ('numCandidateEnhGene', 'numCandidateEnhGene', 'NA', 'max', '0', '# peaks between E and P'),
            ('numNearbyEnhancers', 'numNearbyEnhancers', 'NA', 'max', '0', '# peaks within 5Kb of E'),
            ('sumNearbyEnhancers', 'sumNearbyEnhancers', 'NA', 'max', '0', 'Activity of other peaks within 5Kb of E'),
            ('ubiquitousExpressedGene', 'is_ubiquitous_uniform', 'NA', 'max', '0', 'Ubiquitous expression'),
            ('P2PromoterClass', 'P2PromoterClass', 'NA', 'max', '0', 'Promoter class (Bergman et al. 2022)'),
            ('averageCorrWeighted', 'averageCorrWeighted', 'NA', 'mean', 'mean', 'Correlation of ABC scores for gene across cell type'),
            ('H3K4me3_e_max_L_8', 'H3K4me3_e_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at E, signal max'),
            ('H3K4me3_e_grad_max_L_8', 'H3K4me3_e_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at E, gradient max'),
            ('H3K27ac_e_grad_max_L_8', 'H3K27ac_e_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K27ac at E, gradient max'),
            ('DNase_e_grad_max_L_8', 'DNase_e_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: DNase at E, gradient max'),
            ('H3K4me3_e_grad_min_L_8', 'H3K4me3_e_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at E, gradient min'),
            ('H3K27ac_e_grad_min_L_8', 'H3K27ac_e_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: H3K27ac at E, gradient min'),
            ('DNase_e_grad_min_L_8', 'DNase_e_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: DNase at E, gradient min'),
            ('H3K4me3_p_max_L_8', 'H3K4me3_p_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at P, signal max'),
            ('H3K4me3_p_grad_max_L_8', 'H3K4me3_p_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at P, gradient max'),
            ('H3K27ac_p_grad_max_L_8', 'H3K27ac_p_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: H3K27ac at P, gradient max'),
            ('DNase_p_grad_max_L_8', 'DNase_p_grad_max_L_8', 'NA', 'max', '0', 'GraphReg: DNase at P, gradient max'),
            ('H3K4me3_p_grad_min_L_8', 'H3K4me3_p_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: H3K4me3 at P, gradient min'),
            ('H3K27ac_p_grad_min_L_8', 'H3K27ac_p_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: H3K27ac at P, gradient min'),
            ('DNase_p_grad_min_L_8', 'DNase_p_grad_min_L_8', 'NA', 'max', '0', 'GraphReg: DNase at P, gradient min'),
            ('EpiMapScore', 'EpiMapScore', 'NA', 'mean', '0', 'EpiMap score'),
            ('glsCoefficient', 'glsCoefficient', 'NA', 'mean', 'mean', 'DNase E-P correlation (GLS)'),
            ('PEToutsideNormalized', 'PEToutsideNormalized', 'NA', 'max', '0', 'CTCF PET counts spanning E-P pair'),
            ('PETcrossNormalized', 'PETcrossNormalized', 'NA', 'max', '0', 'CTCF PET counts crossing E-P pair'),
            ('promCTCF', 'promCTCF', 'NA', 'max', '0', 'CTCF signal at P'),
            ('enhCTCF', 'enhCTCF', 'NA', 'max', '0', 'CTCF signal at E'),
            ('HiCLoopOutsideNormalized', 'HiCLoopOutsideNormalized', 'NA', 'max', '0', 'Hi-C loops spanning E-P pair'),
            ('HiCLoopCrossNormalized', 'HiCLoopCrossNormalized', 'NA', 'max', '0', 'Hi-C loops crossing E-P pair'),
            ('inTAD', 'inTAD', 'NA', 'max', '0', 'In TAD'),
            ('inCCD', 'inCCD', 'NA', 'max', '0', 'In CCD'),
            ('normalizedEP300_enh', 'normalizedEP300_enh', 'NA', 'mean', '0', 'EP300 signal at E'),
            ('normalizedDNase_enh', 'normalized_dhs_enh', 'NA', 'mean', '0', 'DNase signal at E'),
            ('normalizedDNase_prom', 'normalized_dhs_prom', 'NA', 'mean', '0', 'DNase signal at P'),
            ('normalizedH3K27ac_enh', 'normalized_h3k27ac_enh', 'NA', 'mean', '0', 'H3K27ac signal at E'),
            ('normalizedH3K27ac_prom', 'normalized_h3k27ac_prom', 'NA', 'mean', '0', 'H3K27ac signal at P'),
        ],
    },
}


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"encode_re2g_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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


def _sklearn_ok() -> bool:
    try:
        import sklearn  # type: ignore  # noqa: F401
        return True
    except Exception:
        return False


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
    lines = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
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
        return f"{float(x):.{nd}f}"
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
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if math.isfinite(f) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_tsv(df, path: Path, header: bool = True, quiet: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, header=header, na_rep="NA",
              compression="gzip" if str(path).endswith(".gz") else None)
    if not quiet:
        print(f"TSV: {path}")
    return path


def write_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(_jsonable(obj), indent=2))
    print(f"JSON: {path}")
    return path


def write_report(d: Path, title: str, sections: "list[str]", summary: dict) -> None:
    rp = d / "report.md"
    rp.write_text(f"# {title}\n\n" + "\n\n".join(sections) + "\n")
    print(f"Report: {rp}")
    write_json(summary, d / "summary.json")


def read_tsv(path, **kw):
    pd = _pd()
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "t", "1", "1.0", "yes", "y"}


def _missing(v) -> bool:
    if v is None:
        return True
    if isinstance(v, float) and math.isnan(v):
        return True
    return str(v).strip() in {"", "NA", "nan", "None", "NaN", "<NA>"}


def to_float(s):
    """Series -> float64: bools and TRUE/FALSE strings -> 1/0, unparsable -> NaN."""
    pd, np = _pd(), _np()
    if s.dtype == bool:
        return s.astype(np.float64)
    if s.dtype == object:
        low = s.astype(str).str.strip().str.lower()
        mapped = low.map({"true": 1.0, "false": 0.0})
        num = pd.to_numeric(s, errors="coerce")
        return num.where(mapped.isna(), mapped).astype(np.float64)
    return pd.to_numeric(s, errors="coerce").astype(np.float64)


def to_bool(s):
    np = _np()
    if s.dtype == bool:
        return s
    return s.map(lambda v: (not _missing(v)) and _truthy(v)).astype(bool) if len(s) else s.astype(bool)


def _download(url: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"Wrote: {dest} (cached)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh)
    tmp.replace(dest)
    print(f"Wrote: {dest}")
    return dest


def resource(name: str, override: "Optional[str]" = None, required: bool = True) -> "Optional[Path]":
    if override:
        return Path(override)
    p = RES_DIR / name
    if p.exists():
        return p
    if required:
        raise SystemExit(f"missing {name}: pass it explicitly or run `igvfagent encode-re2g setup`")
    return None


# ---------------------------------------------------------------------------
# Metrics (training_functions.py) -- numpy, verified against sklearn in the selftest
# ---------------------------------------------------------------------------

def precision_recall_curve(y_true, y_score):
    """sklearn.metrics.precision_recall_curve (drop_intermediate=False)."""
    np = _np()
    y = np.asarray(y_true, dtype=np.float64)
    s = np.asarray(y_score, dtype=np.float64)
    order = np.argsort(s, kind="mergesort")[::-1]
    s, y = s[order], y[order]
    distinct = np.where(np.diff(s))[0]
    idx = np.r_[distinct, y.size - 1]
    tps = np.cumsum(y)[idx]
    fps = 1 + idx - tps
    thresholds = s[idx]
    ps = tps + fps
    precision = np.zeros_like(tps)
    np.divide(tps, ps, out=precision, where=(ps != 0))
    recall = np.ones_like(tps) if tps[-1] == 0 else tps / tps[-1]
    return np.hstack((precision[::-1], 1.0)), np.hstack((recall[::-1], 0.0)), thresholds[::-1]


def auc_xy(x, y) -> float:
    """sklearn.metrics.auc: trapezoid, x monotonic in either direction."""
    np = _np()
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    direction = 1.0
    dx = np.diff(x)
    if np.any(dx < 0):
        if np.all(dx <= 0):
            direction = -1.0
        else:
            raise ValueError("x is neither increasing nor decreasing")
    return float(direction * np.sum(dx * (y[1:] + y[:-1]) / 2.0))


def auc_mod(recall, precision) -> float:
    return 0.0 if len(recall) < 2 else auc_xy(recall, precision)


def pr_curve_modified(y_true, y_pred):
    p, r, t = precision_recall_curve(y_true, y_pred)
    return p[1:], r[1:], t


def statistic_aupr(y_true, y_pred) -> float:
    p, r, _ = pr_curve_modified(y_true, y_pred)
    return auc_mod(r, p)


def statistic_delta_aupr(y_true, y_pred_full, y_pred_ablated) -> float:
    return statistic_aupr(y_true, y_pred_ablated) - statistic_aupr(y_true, y_pred_full)


def statistic_precision(y_true, y_pred) -> float:
    """precision at the point whose recall is closest to 70%."""
    np = _np()
    p, r, _ = pr_curve_modified(y_true, y_pred)
    return float(p[int(np.argmin(np.abs(r - 0.7)))])


def threshold_70_pct_recall(y_true, y_pred):
    np = _np()
    p, r, t = pr_curve_modified(y_true, y_pred)
    if np.max(r) > 0.7:
        i = int(np.argmin(np.abs(r - 0.7)))
        return float(t[i + 1])  # upstream offset for the dropped first point
    return None


def statistic_precision_at_threshold(y_true, y_pred, threshold) -> float:
    np = _np()
    if threshold is None:
        return 0.0
    p, r, t = pr_curve_modified(y_true, y_pred)
    i = int(np.argmin(np.abs(t - threshold)))
    return float(p[i - 1])  # upstream index (i == 0 wraps to the last point, as in the original)


def statistic_delta_precision_at_threshold(y_true, y_full, y_ablated, thresh_full, thresh_ablated) -> float:
    return (statistic_precision_at_threshold(y_true, y_ablated, thresh_ablated)
            - statistic_precision_at_threshold(y_true, y_full, thresh_full))


def roc_auc(y_true, y_score) -> float:
    """Mann-Whitney AUROC (average ranks for ties) == sklearn roc_auc_score."""
    np = _np()
    from scipy.stats import rankdata  # type: ignore
    y = np.asarray(y_true).astype(bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(np.asarray(y_score, dtype=np.float64))
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def log_loss(y_true, p) -> float:
    """sklearn log_loss for binary labels (probabilities clipped to [eps, 1 - eps], float64 eps)."""
    np = _np()
    y = np.asarray(y_true, dtype=np.float64)
    q = np.clip(np.asarray(p, dtype=np.float64), np.finfo(np.float64).eps, 1 - np.finfo(np.float64).eps)
    return float(-np.mean(y * np.log(q) + (1 - y) * np.log(1 - q)))


class BootResult:
    def __init__(self, dist, low, high):
        self.bootstrap_distribution = dist
        self.confidence_interval = (low, high)


def bootstrap(data: tuple, stat, n_boot: int, rng) -> BootResult:
    """scipy.stats.bootstrap(paired=True, method='BCa', confidence_level=0.95); percentile fallback without scipy."""
    np = _np()
    data = tuple(np.asarray(d) for d in data)
    try:
        from scipy import stats as st  # type: ignore
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = st.bootstrap(data, stat, n_resamples=n_boot, paired=True, confidence_level=0.95, method="BCa",
                               random_state=rng, vectorized=False)
        return BootResult(np.asarray(res.bootstrap_distribution), float(res.confidence_interval[0]),
                          float(res.confidence_interval[1]))
    except ImportError:
        n = len(data[0])
        dist = np.array([stat(*[d[idx] for d in data]) for idx in (rng.integers(0, n, n) for _ in range(n_boot))])
        return BootResult(dist, float(np.percentile(dist, 2.5)), float(np.percentile(dist, 97.5)))


def bootstrap_pvalue(delta: float, res: BootResult) -> float:
    np = _np()
    centred = res.bootstrap_distribution - res.bootstrap_distribution.mean()
    return float(np.mean(np.abs(centred) >= abs(delta)))


# ---------------------------------------------------------------------------
# Models: embedded coefficients, model dirs, model selection, params, fitting
# ---------------------------------------------------------------------------

def feature_table_df(rows_or_path):
    """feature_table.tsv (feature, input_col, second_input, aggregate_function, fill_value, nice_name) as strings."""
    pd = _pd()
    if isinstance(rows_or_path, (str, Path)):
        df = pd.read_csv(rows_or_path, sep="\t", dtype=str, keep_default_na=False)
    else:
        df = pd.DataFrame(list(rows_or_path), columns=FT_COLS, dtype=str)
    for c in FT_COLS:
        if c not in df.columns:
            df[c] = ""
    return df[FT_COLS + [c for c in df.columns if c not in FT_COLS]].reset_index(drop=True)


def ft_value(v) -> "Optional[str]":
    return None if _missing(v) else str(v).strip()


def combine_feature_tables(tables: "list[Any]"):
    """combine_feature_tables(_apply).R: rbind, add the sc-E2G ARC rows when needed, distinct."""
    pd = _pd()
    df = pd.concat([feature_table_df(t) if not hasattr(t, "columns") else t for t in tables], ignore_index=True)
    if df["feature"].isin(["ARC.E2G.Score", "Kendall"]).any():
        df = pd.concat([df, feature_table_df(ARC_ROWS)], ignore_index=True)
    return df.drop_duplicates().reset_index(drop=True)


def input_features(ft) -> "list[str]":
    out = []
    for v in list(ft["input_col"]) + list(ft["second_input"]):
        v = ft_value(v)
        if v and v not in out:
            out.append(v)
    return out


def poly2(X):
    """sklearn PolynomialFeatures(degree=2) with get_feature_names_out naming ('1', 'a', 'a^2', 'a b')."""
    pd, np = _pd(), _np()
    cols = list(X.columns)
    A = X.to_numpy(dtype=np.float64)
    out = {"1": np.ones(len(X))}
    for i, c in enumerate(cols):
        out[c] = A[:, i]
    for i, j in itertools.combinations_with_replacement(range(len(cols)), 2):
        out[f"{cols[i]}^2" if i == j else f"{cols[i]} {cols[j]}"] = A[:, i] * A[:, j]
    return pd.DataFrame(out, index=X.index)


def design_matrix(df, features: "list[str]", epsilon: float = EPSILON, polynomial: bool = False):
    """X = log(|x| + epsilon) of the feature columns (after degree-2 expansion when polynomial)."""
    pd, np = _pd(), _np()
    missing = [f for f in features if f not in df.columns]
    if missing:
        raise SystemExit(f"features missing from the table: {missing}")
    X = pd.DataFrame({f: to_float(df[f]) for f in features}, index=df.index)
    if polynomial:
        X = poly2(X)
    return np.log(np.abs(X) + epsilon)


def _expit(z):
    np = _np()
    return 1.0 / (1.0 + np.exp(-z))


class Re2gModel:
    """A logistic-regression E2G model: feature table + coefficients (+ optional sklearn object)."""

    def __init__(self, name: str, feature_table, coef, intercept: float, threshold: "Optional[float]",
                 threshold_str: "Optional[str]" = None, source: str = "", sk=None, polynomial: bool = False,
                 design_columns: "Optional[list[str]]" = None):
        np = _np()
        self.name = name
        self.feature_table = feature_table
        self.features = list(feature_table["feature"])
        self.coef = np.asarray(coef, dtype=np.float64)
        self.intercept = float(intercept)
        self.threshold = threshold
        self.threshold_str = threshold_str
        self.source = source
        self.sk = sk
        self.polynomial = polynomial
        self.design_columns = design_columns or self.features
        if len(self.coef) != len(self.design_columns):
            raise SystemExit(f"model {name}: {len(self.coef)} coefficients for {len(self.design_columns)} features")

    def design(self, df, epsilon: float = EPSILON):
        X = design_matrix(df, self.features, epsilon, self.polynomial)
        return X[self.design_columns]

    def predict_proba(self, df, epsilon: float = EPSILON, use_pickle: bool = False):
        np = _np()
        X = self.design(df, epsilon)
        if use_pickle and self.sk is not None:
            return self.sk.predict_proba(X)[:, 1]
        return _expit(X.to_numpy(dtype=np.float64) @ self.coef + self.intercept)

    def to_json(self) -> dict:
        return {"name": self.name, "features": self.features, "design_columns": self.design_columns,
                "coef": [float(c) for c in self.coef], "intercept": self.intercept, "threshold": self.threshold,
                "polynomial": self.polynomial, "epsilon": EPSILON, "source": self.source,
                "transform": "log(|x| + epsilon); p = 1 / (1 + exp(-(X.coef + intercept)))"}

    def write_dir(self, d: Path, write_pickle: bool = False) -> Path:
        d.mkdir(parents=True, exist_ok=True)
        write_tsv(self.feature_table[FT_COLS], d / "feature_table.tsv", quiet=True)
        for old in d.glob("threshold_*"):
            old.unlink()
        if self.threshold is not None:
            (d / f"threshold_{self.threshold_str or format_threshold(self.threshold)}").write_text("")
        (d / "model.json").write_text(json.dumps(self.to_json(), indent=2))
        if write_pickle and self.sk is not None:
            with open(d / "model.pkl", "wb") as fh:
                pickle.dump(self.sk, fh)
        print(f"Wrote: {d}")
        return d


def format_threshold(t: float) -> str:
    return f"{t:.3f}".rstrip("0").rstrip(".") if t is not None else ""


def embedded_model(name: str) -> Re2gModel:
    m = EMBEDDED_MODELS[name]
    return Re2gModel(name, feature_table_df(m["feature_table"]), m["coef"], m["intercept"], m["threshold"],
                     m["threshold_file"].split("_", 1)[1], source=f"embedded:{UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}/models/{name}")


def validate_model_dir(d: Path) -> None:
    """utils.smk::_validate_model_dir (model.json accepted in place of model.pkl)."""
    files = os.listdir(d)
    if "model.pkl" not in files and "model.json" not in files:
        raise SystemExit(f"model.pkl (or model.json) is not provided in the specified model directory {d}")
    if "feature_table.tsv" not in files:
        raise SystemExit(f"feature_table.tsv is not provided in the specified model directory {d}")
    n = sum(s.startswith("threshold_") for s in files)
    if n == 0:
        raise SystemExit(f"A threshold is not provided in the specified model directory {d}")
    if n > 1:
        raise SystemExit(f"More than one threshold is provided in the specified model directory {d}")


def model_threshold(d: Path) -> "tuple[float, str]":
    files = glob.glob(os.path.join(str(d), "threshold_*"))
    assert len(files) == 1, "Should have exactly 1 threshold file in directory"
    s = os.path.basename(files[0]).split("_")[1]
    return float(s), s


def load_model(spec: str, prefer_pickle: bool = False) -> Re2gModel:
    """An embedded model name or a model directory (model.json preferred, else the sklearn model.pkl)."""
    p = Path(spec)
    if spec in EMBEDDED_MODELS and not p.is_dir():
        return embedded_model(spec)
    if not p.is_dir():
        raise SystemExit(f"{spec}: neither an embedded model ({', '.join(EMBEDDED_MODELS)}) nor a model directory")
    validate_model_dir(p)
    ft = feature_table_df(p / "feature_table.tsv")
    thr, thr_s = model_threshold(p)
    sk = None
    if (p / "model.pkl").exists() and (prefer_pickle or not (p / "model.json").exists()):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with open(p / "model.pkl", "rb") as fh:
                    sk = pickle.load(fh)
        except Exception as e:
            if not (p / "model.json").exists():
                raise SystemExit(f"{p}/model.pkl could not be unpickled ({e}); add a model.json with coef/intercept") from e
    if sk is not None:
        names = list(getattr(sk, "feature_names_in_", ft["feature"]))
        poly = names != list(ft["feature"])
        return Re2gModel(p.name, ft, sk.coef_[0], float(sk.intercept_[0]), thr, thr_s, source=str(p / "model.pkl"),
                         sk=sk, polynomial=poly, design_columns=names)
    js = json.loads((p / "model.json").read_text())
    return Re2gModel(p.name, ft, js["coef"], js["intercept"], thr, thr_s, source=str(p / "model.json"),
                     polynomial=bool(js.get("polynomial")), design_columns=js.get("design_columns"))


def model_name_for_row(row: dict, megamap: str = MEGAMAP_HIC_FILE) -> str:
    """utils.smk::_get_biosample_model_dir_from_row -> the models/<name> folder name."""
    params = [str(row["default_accessibility_feature"]).lower()]
    if not _missing(row.get("H3K27ac")):
        params.append("h3k27ac")
    hic = row.get("HiC_file")
    if _missing(hic):
        params.append("powerlaw")
    elif str(row.get("HiC_type")) == "avg":
        params.append("avg_hic")
    elif os.path.basename(str(hic)) == os.path.basename(megamap):
        params.append("megamap")
    elif str(row.get("HiC_type")) == "hic":
        params.append("intact_hic")
    return "_".join(params)


def select_models(biosamples, model_root: "Optional[str]" = None, megamap: str = MEGAMAP_HIC_FILE) -> "list[dict]":
    """utils.smk::expand_biosample_df: one row per (biosample, model dir) with model_dir_base and model_threshold."""
    out = []
    for row in biosamples.to_dict("records"):
        if "model_dir" in row and not _missing(row["model_dir"]):
            specs = [s.strip() for s in str(row["model_dir"]).split(",") if s.strip()]
            for s in specs:
                if Path(s).is_dir():
                    validate_model_dir(Path(s))
        else:
            name = model_name_for_row(row, megamap)
            if model_root:
                folder = os.path.normpath(os.path.join(model_root, name))
                if not os.path.exists(folder):
                    raise SystemExit(f"{folder} not found. Model with input params not supported")
                specs = [folder]
            else:
                if name not in EMBEDDED_MODELS:
                    raise SystemExit(f"models/{name} not found. Model with input params not supported "
                                     f"(trained: {', '.join(sorted(EMBEDDED_MODELS))})")
                specs = [name]
        for s in specs:
            m = load_model(s)
            out.append({"biosample": row["biosample"], "model_dir": s, "model_dir_base": os.path.basename(os.path.normpath(s)),
                        "model_threshold": m.threshold, "model": m, "row": row})
    return out


def get_params(default_params: dict, override_params) -> dict:
    """get_params.py: overrides replace defaults, 'null' -> None, numeric strings -> numbers."""
    final = dict(default_params)
    if isinstance(override_params, str):
        s = override_params.strip()
        if s and s not in ("None", "nan"):
            s = s.replace("'", '"').replace("True", "true").replace("False", "false").replace("None", "null")
            override_params = json.loads(s)
        else:
            override_params = None
    if isinstance(override_params, dict):
        final.update(override_params)
    final = {k: (None if v == "null" else v) for k, v in final.items()}
    for k, v in list(final.items()):
        if isinstance(v, str):
            try:
                final[k] = int(float(v)) if k in ("max_iter", "random_state", "n_jobs") else float(v)
            except ValueError:
                pass
    return final


class FittedLR:
    def __init__(self, coef, intercept, sk=None, columns=None):
        self.coef, self.intercept, self.sk, self.columns = coef, float(intercept), sk, columns

    def predict_proba1(self, X):
        np = _np()
        if self.sk is not None:
            return self.sk.predict_proba(X)[:, 1]
        return _expit(X.to_numpy(dtype=np.float64) @ self.coef + self.intercept)


def _irls(X, y, fit_intercept: bool = True, l2: float = 0.0, max_iter: int = 200, tol: float = 1e-10):
    """Newton-Raphson logistic regression (numpy fallback when scikit-learn is unavailable)."""
    np = _np()
    A = np.column_stack([np.ones(len(X)), X]) if fit_intercept else np.asarray(X)
    w = np.zeros(A.shape[1])
    pen = np.full(A.shape[1], l2)
    if fit_intercept:
        pen[0] = 0.0
    for _ in range(max_iter):
        p = _expit(A @ w)
        g = A.T @ (p - y) + pen * w
        H = (A * (p * (1 - p))[:, None]).T @ A + np.diag(pen) + 1e-10 * np.eye(A.shape[1])
        step = np.linalg.solve(H, g)
        w_new = w - step
        if np.max(np.abs(step)) < tol:
            w = w_new
            break
        w = w_new
    return (w[1:], w[0]) if fit_intercept else (w, 0.0)


def fit_lr(X, y, params: dict) -> FittedLR:
    """LogisticRegression(**params).fit(X, y); numpy Newton fallback (penalty None or l2) without scikit-learn."""
    np = _np()
    y = np.asarray(y).astype(np.int64)
    if _sklearn_ok():
        from sklearn.linear_model import LogisticRegression  # type: ignore
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = LogisticRegression(**params).fit(X, y)
        return FittedLR(m.coef_[0], m.intercept_[0], sk=m, columns=list(X.columns))
    pen = params.get("penalty")
    l2 = 0.0 if pen in (None, "none") else 1.0 / float(params.get("C", 1.0))
    coef, b = _irls(X.to_numpy(dtype=np.float64), y.astype(np.float64), bool(params.get("fit_intercept", True)), l2)
    return FittedLR(coef, b, columns=list(X.columns))


# ---------------------------------------------------------------------------
# Interval helpers (bedtools intersect / GenomicRanges findOverlaps replacements)
# ---------------------------------------------------------------------------

KEY_SPAN = 1 << 32  # > any chromosome length; keyed intervals live on disjoint stretches of one virtual axis


def keyed_overlap(a, b, a_cols: "tuple[str, str, str, str]", b_cols: "tuple[str, str, str, str]", closed: bool = True):
    """Row-index pairs (ia, ib) with equal chr:gene key and overlapping intervals.

    a_cols / b_cols = (chr, gene, start, end).  closed=True reproduces IRanges(start, end) semantics
    (both ends inclusive, so touching intervals overlap); closed=False is BED half-open.
    Uses eqtl_enrichment_skill.overlap_join on a virtual axis: key_id * 2^32 + coordinate.
    """
    pd, np = _pd(), _np()
    from eqtl_enrichment_skill import overlap_join  # type: ignore
    ka = a[a_cols[0]].astype(str).str.cat(a[a_cols[1]].astype(str), sep=":")
    kb = b[b_cols[0]].astype(str).str.cat(b[b_cols[1]].astype(str), sep=":")
    codes, _ = pd.factorize(pd.concat([ka, kb], ignore_index=True))
    ca, cb = codes[:len(a)].astype(np.int64), codes[len(a):].astype(np.int64)
    shared = np.intersect1d(ca, cb)
    ma, mb = np.isin(ca, shared), np.isin(cb, shared)
    if not ma.any():
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    add = 1 if closed else 0
    A = pd.DataFrame({"chr": "k", "start": ca[ma] * KEY_SPAN + to_float(a[a_cols[2]]).to_numpy()[ma].astype(np.int64),
                      "end": ca[ma] * KEY_SPAN + to_float(a[a_cols[3]]).to_numpy()[ma].astype(np.int64) + add,
                      "ia": np.flatnonzero(ma)})
    B = pd.DataFrame({"chr": "k", "start": cb[mb] * KEY_SPAN + to_float(b[b_cols[2]]).to_numpy()[mb].astype(np.int64),
                      "end": cb[mb] * KEY_SPAN + to_float(b[b_cols[3]]).to_numpy()[mb].astype(np.int64) + add,
                      "ib": np.flatnonzero(mb)})
    j = overlap_join(A, B)
    if len(j) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    order = np.lexsort((j["ib"].to_numpy(), j["ia"].to_numpy()))
    return j["ia"].to_numpy().astype(np.int64)[order], j["ib"].to_numpy().astype(np.int64)[order]


AGG_FUNS = {"max", "min", "mean", "sum", "median", "first", "last"}


def aggregate_hits(ia, values, func: str):
    """R-style aggregation of a feature over the hits of each query row: any NA in a group -> NA."""
    pd, np = _pd(), _np()
    func = str(func).strip()
    if func not in AGG_FUNS:
        raise SystemExit(f"unsupported aggregate_function {func!r} (supported: {sorted(AGG_FUNS)})")
    s = to_float(pd.Series(values).reset_index(drop=True))
    g = s.groupby(ia)
    res = getattr(g, func)()
    has_na = s.isna().groupby(ia).any()
    res[has_na.reindex(res.index).to_numpy()] = np.nan
    return res


def merge_by_overlap(base, feats, base_cols, feat_cols, score_cols: "list[str]", agg_funs: "list[str]"):
    """Shared core of overlap_feature_with_abc and merge_feature_to_crispr: (merged, missing)."""
    pd, np = _pd(), _np()
    ia, ib = keyed_overlap(base, feats, base_cols, feat_cols, closed=True)
    base = base.reset_index(drop=True)
    feats = feats.reset_index(drop=True)
    uniq = np.unique(ia)
    merged = base.iloc[uniq].copy()
    for col, fn in zip(score_cols, agg_funs):
        agg = aggregate_hits(ia, feats[col].to_numpy()[ib], fn)
        merged[col] = agg.reindex(uniq).to_numpy()
    missing = base.drop(index=uniq).copy()
    for col in score_cols:
        missing[col] = np.nan
    return merged, missing


# ---------------------------------------------------------------------------
# Genome-wide features (gen_new_features.py, activity_only_features.R, merge_external_features.R,
# gen_final_features.R, get_fill_values.R)
# ---------------------------------------------------------------------------

def add_midpoint(df) -> None:
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype("int64")


def num_candidate_enh_gene(pred):
    """NumCandidateEnhGene: rank of each non-promoter candidate from the TSS, per (TargetGene, TargetGeneTSS) and side."""
    pd, np = _pd(), _np()
    df = pred.sort_values(["chr", "midpoint"], kind="mergesort").reset_index(drop=True)
    out = np.zeros(len(df), dtype=np.int64)
    mid, tss = df["midpoint"].to_numpy(), to_float(df["TargetGeneTSS"]).to_numpy()
    keys = pd.DataFrame({"g": df["TargetGene"].to_numpy(), "t": tss})
    down = np.flatnonzero(mid > tss)
    if down.size:
        out[down] = keys.iloc[down].groupby(["g", "t"], sort=False).cumcount().to_numpy() + 1
    up = np.flatnonzero(mid < tss)[::-1]  # closest to the TSS first
    if up.size:
        out[up] = keys.iloc[up].groupby(["g", "t"], sort=False).cumcount().to_numpy() + 1
    return pd.DataFrame({"name": df["name"].to_numpy(), "TargetGene": df["TargetGene"].to_numpy(), "NumCandidateEnhGene": out})


def _interval_counts(q_start, q_end, t_start, t_end, t_weight=None):
    """For each query [s, e): number (and weight sum) of target intervals with t.start < e and t.end > s."""
    np = _np()
    so = np.argsort(t_start, kind="mergesort")
    eo = np.argsort(t_end, kind="mergesort")
    ss, es = t_start[so], t_end[eo]
    n_before_end = np.searchsorted(ss, q_end, side="left")
    n_end_before_start = np.searchsorted(es, q_start, side="right")
    cnt = n_before_end - n_end_before_start
    if t_weight is None:
        return cnt, None
    w = np.nan_to_num(np.asarray(t_weight, dtype=np.float64))
    cs = np.r_[0.0, np.cumsum(w[so])]
    ce = np.r_[0.0, np.cumsum(w[eo])]
    return cnt, cs[n_before_end] - ce[n_end_before_start]


def read_bed(path, names: "list[str]"):
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    df = df.iloc[:, :len(names)]
    df.columns = names[:df.shape[1]]
    return df


def num_tss_enh_gene(pred, tss_bed):
    """NumTSSEnhGene: TSS500bp intervals between E and P (bedtools intersect -wa -wb | groupby size)."""
    pd, np = _pd(), _np()
    mid = pred["midpoint"].to_numpy().astype(np.int64)
    tss = to_float(pred["TargetGeneTSS"]).to_numpy()
    s = pred["start"].to_numpy().astype(np.int64).copy()
    e = (mid + to_float(pred["distance"]).to_numpy()).astype(np.int64)
    dn = tss < mid  # gene upstream of the enhancer: [TSS, end)
    s[dn] = tss[dn].astype(np.int64)
    e[dn] = pred["end"].to_numpy().astype(np.int64)[dn]
    counts = np.zeros(len(pred), dtype=np.int64)
    chrs = pred["chr"].astype(str).to_numpy()
    for c, sub in tss_bed.groupby(tss_bed["chr"].astype(str)):
        m = chrs == c
        if m.any():
            counts[m], _ = _interval_counts(s[m], e[m], sub["start"].to_numpy().astype(np.int64),
                                            sub["end"].to_numpy().astype(np.int64))
    keep = counts > 0
    return pd.DataFrame({"class": pred["name"].to_numpy()[keep], "gene": pred["TargetGene"].to_numpy()[keep],
                         "0": counts[keep]})


def nearby_enhancers(pred_all, enh_list, chr_sizes: "Optional[dict]" = None, window: int = NEARBY_WINDOW,
                     dedupe: bool = False):
    """NumEnhancersEG5kb / SumEnhancersEG5kb: prediction rows within midpoint +/- window, same name excluded."""
    pd, np = _pd(), _np()
    rows = pred_all[["chr", "start", "end", "name", "activity_base"]]
    if dedupe:
        rows = rows.drop_duplicates()
    el = enh_list[["chr", "start", "end", "name"]].copy()
    add_midpoint(el)
    ws = np.maximum(0, el["midpoint"].to_numpy() - window)
    we = el["midpoint"].to_numpy() + window
    if chr_sizes:
        sizes = el["chr"].map(chr_sizes).to_numpy(dtype=np.float64)
        we = np.where(np.isnan(sizes), we, np.minimum(we, np.nan_to_num(sizes, nan=0).astype(np.int64)))
    cnt = np.zeros(len(el), dtype=np.int64)
    sm = np.zeros(len(el), dtype=np.float64)
    chrs = el["chr"].astype(str).to_numpy()
    act = to_float(rows["activity_base"]).to_numpy()
    for c, idx in rows.groupby(rows["chr"].astype(str)).indices.items():
        m = chrs == c
        if m.any():
            sub = rows.iloc[idx]
            cnt[m], sm[m] = _interval_counts(ws[m], we[m], sub["start"].to_numpy().astype(np.int64),
                                             sub["end"].to_numpy().astype(np.int64), act[idx])
    # remove the element's own rows (they always overlap its own window)
    own = pd.DataFrame({"name": rows["name"].to_numpy(), "a": np.nan_to_num(act),
                        "ov": 1}).groupby("name").agg(n=("ov", "sum"), s=("a", "sum"))
    own = own.reindex(el["name"].to_numpy())
    # only subtract own rows that actually overlap the window (always true unless names are reused elsewhere)
    cnt = cnt - np.nan_to_num(own["n"].to_numpy()).astype(np.int64)
    sm = sm - np.nan_to_num(own["s"].to_numpy())
    keep = cnt > 0
    num = pd.DataFrame({"name": el["name"].to_numpy()[keep], "count": cnt[keep]})
    tot = pd.DataFrame({"name": el["name"].to_numpy()[keep], "sum": sm[keep]})
    # upstream groupby output is sorted by name
    return num.sort_values("name", kind="mergesort").reset_index(drop=True), tot.sort_values("name", kind="mergesort").reset_index(drop=True)


def read_chr_sizes(path) -> dict:
    if not path:
        return {}
    d = read_bed(path, ["chr", "size"])
    return dict(zip(d["chr"].astype(str), d["size"].astype("int64")))


def activity_only_features(abc, ft, cand, tss_cnt, near_n, near_s, gene_classes=None):
    """activity_only_features.R: core ABC columns + requested ABC / new / gene-class features (left joins)."""
    pd = _pd()
    abc = abc.copy()
    if "ABC.Score.Numerator" in abc.columns:
        abc["ABC.Numerator"] = abc["ABC.Score.Numerator"]
        abc["ABC.Denominator"] = to_float(abc["ABC.Score"]) / to_float(abc["ABC.Numerator"])
    inp = input_features(ft)
    core = [c for c in CORE_COLS if c in abc.columns and c not in inp]
    out = abc[core].copy()
    from_abc = [c for c in inp if c in abc.columns and c not in ("name", "TargetGene")]
    if not abc.duplicated(["name", "TargetGene"]).any():
        for c in from_abc:
            out[c] = abc[c].to_numpy()
    else:
        out = out.merge(abc[["name", "TargetGene"] + from_abc], on=["name", "TargetGene"], how="left")
    if "numCandidateEnhGene" in inp:
        t = cand.rename(columns={"NumCandidateEnhGene": "numCandidateEnhGene"})[["name", "TargetGene", "numCandidateEnhGene"]]
        out = out.merge(t, on=["name", "TargetGene"], how="left")
    if "numTSSEnhGene" in inp:
        t = tss_cnt.copy()
        t.columns = ["name", "TargetGene", "numTSSEnhGene"]
        out = out.merge(t, on=["name", "TargetGene"], how="left")
    if "numNearbyEnhancers" in inp:
        out = out.merge(near_n.set_axis(["name", "numNearbyEnhancers"], axis=1), on="name", how="left")
    if "sumNearbyEnhancers" in inp:
        out = out.merge(near_s.set_axis(["name", "sumNearbyEnhancers"], axis=1), on="name", how="left")
    if gene_classes is not None:
        gc_cols = [c for c in inp if c in gene_classes.columns and c != "TargetGene"]
        if gc_cols:
            out = out.merge(gene_classes[["TargetGene"] + gc_cols], on="TargetGene", how="left")
    return out


def read_external_config(path):
    pd = _pd()
    if not path:
        return pd.DataFrame(columns=["input_col", "source_col", "aggregate_function", "join_by", "source_file"])
    cfg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "input_col" not in cfg.columns and "feature" in cfg.columns:
        cfg = cfg.rename(columns={"feature": "input_col"})
    base = Path(path).parent
    # format_external_features_config.R: relative source files resolved against the pipeline dir (here: the config's dir)
    cfg["source_file"] = [str(base / f) if f and not os.path.isabs(f) and (base / f).exists() else f for f in cfg["source_file"]]
    return cfg


def merge_external_features(abc, ft, ext_cfg):
    """merge_external_features.R: join_by overlap (chr:TargetGene, aggregate) or join_by TargetGene."""
    pd = _pd()
    notes = []
    inp = input_features(ft)
    needed = [f for f in inp if f not in abc.columns]
    ext = ext_cfg[ext_cfg["input_col"].isin(needed)]
    for src in [s for s in ext["source_file"].unique() if not _missing(s)]:
        this = ext[ext["source_file"] == src]
        source = read_tsv(src)
        source = source.rename(columns=dict(zip(this["source_col"], this["input_col"])))
        cols = list(this["input_col"])
        if "overlap" in set(this["join_by"]):
            merged, missing = merge_by_overlap(abc, source, ("chr", "TargetGene", "start", "end"),
                                               ("chr", "TargetGene", "start", "end"), cols, list(this["aggregate_function"]))
            notes.append(f"{src}: features overlapping ABC {100.0 * len(merged) / max(len(abc), 1):.2f}%")
            abc = pd.concat([merged, missing], ignore_index=True).sort_values(
                ["chr", "start", "end", "TargetGene"], kind="mergesort").reset_index(drop=True)
        elif (this["join_by"] == "TargetGene").all():
            n_abc = abc["TargetGene"].nunique()
            n_rep = len(set(abc["TargetGene"]) & set(source["TargetGene"]))
            notes.append(f"{src}: genes overlapping ABC {100.0 * n_rep / max(n_abc, 1):.2f}%")
            abc = abc.merge(source[["TargetGene"] + cols], on="TargetGene", how="left")
        else:
            notes.append(f"{src}: features must be added by overlap or TargetGene (skipped)")
    for n in notes:
        log.info(n)
    return abc, notes


def get_fill_values(features, ft) -> dict:
    """get_fill_values.R: fill_value per feature present; 'mean' = mean of finite values; NA -> no fill."""
    np = _np()
    out = {}
    sub = ft[ft["feature"].isin(features.columns)][["feature", "fill_value"]].drop_duplicates()
    for f, v in zip(sub["feature"], sub["fill_value"]):
        v = ft_value(v)
        if v is None:
            out[f] = None
        elif v == "mean":
            x = to_float(features[f]).to_numpy()
            x = x[np.isfinite(x)]
            out[f] = float(x.mean()) if x.size else float("nan")
        else:
            try:
                out[f] = float(v)
            except ValueError:
                out[f] = None
    return out


def fill_na(df, fills: dict):
    for f, v in fills.items():
        if v is not None and f in df.columns and not (isinstance(v, float) and math.isnan(v)):
            if df[f].isna().any():
                df[f] = to_float(df[f]).fillna(v) if df[f].dtype == object else df[f].fillna(v)
    return df


def gen_final_features(df, ft):
    """gen_final_features.R: interaction terms, rename input_col -> feature, fill NAs (full table is written)."""
    inp = input_features(ft)
    miss = [f for f in inp if f not in df.columns]
    if miss:
        raise SystemExit(f"required features are not present: {miss}")
    df = df.copy()
    intx = ft[[ft_value(v) is not None for v in ft["second_input"]]]
    for f, a, b in zip(intx["feature"], intx["input_col"], intx["second_input"]):
        df[f] = to_float(df[a]) * to_float(df[b])
    single = ft[~ft["feature"].isin(set(intx["feature"]))]
    for f, a in zip(single["feature"], single["input_col"]):
        if a in df.columns and a != f:
            df = df.rename(columns={a: f})
    return fill_na(df, get_fill_values(df, ft))


def generate_features(abc_pred, enhancer_list, tss_path, chr_sizes_path, gene_classes_path, ft, ext_cfg, out_dir: Path,
                      dedupe_nearby: bool = False) -> dict:
    """The genomewide_features.smk chain for one biosample; writes upstream-named intermediate files."""
    pd = _pd()
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_all = read_tsv(abc_pred)
    pred = pred_all[pred_all["class"] != "promoter"].copy()
    if len(pred) == 0:
        raise SystemExit("Did not find any enhancers in the Predictions file")
    add_midpoint(pred)
    cand = num_candidate_enh_gene(pred)
    write_tsv(cand, out_dir / "NumCandidateEnhGene.tsv")
    tss = read_bed(tss_path, ["chr", "start", "end", "name", "score", "strand"])
    tss_cnt = num_tss_enh_gene(pred, tss)
    write_tsv(tss_cnt, out_dir / "NumTSSEnhGene.tsv")
    enh = read_tsv(enhancer_list, usecols=["chr", "start", "end", "name"])
    near_n, near_s = nearby_enhancers(pred_all, enh, read_chr_sizes(chr_sizes_path), dedupe=dedupe_nearby)
    write_tsv(near_n, out_dir / "NumEnhancersEG5kb.txt", header=False)
    write_tsv(near_s, out_dir / "SumEnhancersEG5kb.txt", header=False)
    gc = read_tsv(gene_classes_path) if gene_classes_path else None
    act = activity_only_features(pred_all, ft, cand, tss_cnt, near_n, near_s, gc)
    write_tsv(act, out_dir / "ActivityOnly_features.tsv.gz")
    plus, notes = merge_external_features(act, ft, ext_cfg)
    write_tsv(plus, out_dir / "ActivityOnly_plus_external_features.tsv.gz")
    final = gen_final_features(plus, ft)
    write_tsv(final, out_dir / "genomewide_features.tsv.gz")
    write_tsv(ft, out_dir / "feature_table.tsv")
    return {"n_pairs": int(len(final)), "n_enhancer_pairs": int(len(pred)), "n_elements": int(len(enh)),
            "n_numTSSEnhGene_rows": int(len(tss_cnt)), "n_nearby_rows": int(len(near_n)), "external_notes": notes,
            "features_path": str(out_dir / "genomewide_features.tsv.gz"), "features": final}


# ---------------------------------------------------------------------------
# Model application (run_e2g.py, threshold_e2g_predictions.py, process_model_output.py, get_stats.py)
# ---------------------------------------------------------------------------

def make_e2g_predictions(df, model: Re2gModel, epsilon: float = EPSILON, use_pickle: bool = False):
    np = _np()
    df = df.replace([np.inf, -np.inf], np.nan).fillna(0)
    df[SCORE_COL] = model.predict_proba(df, epsilon, use_pickle)
    return df


def threshold_predictions(df, threshold: float, score_col: str = SCORE_COL, include_self_promoter: bool = True):
    f = df[to_float(df[score_col]) >= threshold]
    if include_self_promoter:
        return f[(f["class"] != "promoter") | to_bool(f["isSelfPromoter"])]
    return f[f["class"] != "promoter"]


def bedpe_frame(pred, score_col: str = SCORE_COL):
    pd = _pd()
    pred = pred[to_bool(pred["TargetGeneIsExpressed"])].drop_duplicates()
    return pd.DataFrame({"chr1": pred["chr"], "x1": pred["start"], "x2": pred["end"], "chr2": pred["chr"],
                         "y1": pred["TargetGeneTSS"], "y2": pred["TargetGeneTSS"],
                         "name": pred["TargetGene"].astype(str) + "_" + pred["name"].astype(str),
                         "score": pred[score_col], "strand1": ".", "strand2": "."})


def count_bam_reads(bam: str) -> int:
    """samtools idxstats mapped reads on chr1-22, X, Y (pysam when available)."""
    try:
        import pysam  # type: ignore
        txt = pysam.idxstats(bam)
    except ImportError:
        if not shutil.which("samtools"):
            print(f"samtools not on PATH; would run: samtools idxstats {bam}")
            return 0
        txt = subprocess.check_output(["samtools", "idxstats", bam]).decode()
    total = 0
    for line in txt.strip().splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0] in NORMAL_CHROMOSOMES:
            total += int(f[2])
    return total


def get_num_reads(files: "list[str]") -> int:
    total = 0
    for f in files:
        if not f.endswith(".bam"):
            print("Only support num reads for bam files")
            return 0
        total += count_bam_reads(f)
    return total


def get_stats(df, accessibility: "Optional[list[str]]" = None):
    pd, np = _pd(), _np()
    enh = df.groupby(["chr", "start", "end"]).size()
    genes = df.groupby("TargetGene").size()
    no_prom = df[df["class"] != "promoter"].groupby("TargetGene").size()
    with np.errstate(divide="ignore"):
        ld = np.log10(to_float(df["distanceToTSS"]))
    ld = ld.replace(-np.inf, 0)
    sizes = pd.Series([e - s for (_, s, e) in enh.index]) if len(enh) else pd.Series(dtype=float)
    stats = [(SEQ_DEPTH_METRIC, get_num_reads(accessibility) if accessibility else 0),
             ("num_enh", len(df[["chr", "start", "end"]].drop_duplicates())),
             ("num_genes", len(df["TargetGene"].drop_duplicates())),
             ("num_enh_gene_links", len(df)),
             ("num_genes_with_1_enh_min", int((genes > 0).sum())),
             ("mean_num_genes_per_enh", enh.mean()),
             ("mean_num_enh_per_gene", genes.mean()),
             ("mean_num_enh_per_gene_no_prom", no_prom.mean()),
             ("mean_log10_dist_to_tss", ld.mean()),
             ("mean_enh_region_size", sizes.mean())]
    return pd.DataFrame({"Metric": [s[0] for s in stats], "Value": [float(s[1]) for s in stats]})


# ---------------------------------------------------------------------------
# CRISPR training data (overlap_features_with_crispr_data.R, process_crispr_data.R)
# ---------------------------------------------------------------------------

def overlap_features_with_crispr(crispr, features, ft, tss, nafill: bool = True):
    """merge_feature_to_crispr + distanceToTSS refill + NA fill -> (merged, missing, notes)."""
    pd, np = _pd(), _np()
    crispr = crispr.drop(columns=[c for c in ("pair_uid", "merged_uid", "merged_start", "merged_end") if c in crispr.columns])
    config = ft[ft["feature"].isin(features.columns)]
    agg = {}
    for f, a in zip(config["feature"], config["aggregate_function"]):
        agg.setdefault(f, a)
    missing_genes = set(crispr["measuredGeneSymbol"]) - set(tss["name"])
    notes = [f"Removing CRISPR data for {len(missing_genes)} genes not part of TSS universe"]
    crispr = crispr[~crispr["measuredGeneSymbol"].isin(missing_genes)].reset_index(drop=True)
    cols = list(dict.fromkeys(config["feature"]))
    merged, missing = merge_by_overlap(crispr, features, ("chrom", "measuredGeneSymbol", "chromStart", "chromEnd"),
                                       ("chr", "TargetGene", "start", "end"), cols, [agg[c] for c in cols])
    notes.append(f"Features overlapping crispr: {100.0 * len(merged) / max(len(crispr), 1):.2f}%")
    order = [c for c in ("dataset", "chrom", "chromStart", "chromEnd", "measuredGeneSymbol") if c in merged.columns]
    merged = merged.sort_values(order, kind="mergesort").reset_index(drop=True)
    ref = tss[["name", "start", "end"]].rename(columns={"name": "measuredGeneSymbol", "start": "startTSS_ref", "end": "endTSS_ref"})
    merged = merged.merge(ref, on="measuredGeneSymbol", how="left")
    missing = missing.merge(ref, on="measuredGeneSymbol", how="left")
    if "distanceToTSS" in merged.columns:
        centre = (to_float(merged["chromStart"]) + to_float(merged["chromEnd"])) / 2
        d = to_float(merged["distanceToTSS"])
        from_crispr = (centre - (to_float(merged["startTSS"]) + to_float(merged["endTSS"])) / 2).abs() \
            if "startTSS" in merged.columns else pd.Series(np.nan, index=merged.index)
        from_ref = (centre - (to_float(merged["startTSS_ref"]) + to_float(merged["endTSS_ref"])) / 2).abs()
        na, ref_na = d.isna(), to_float(merged["startTSS_ref"]).isna()
        merged["distanceToTSS"] = np.where(na & ref_na, from_crispr, np.where(na, from_ref, d))
    merged = merged.drop(columns=["startTSS_ref", "endTSS_ref"])
    if nafill:
        merged = fill_na(merged, get_fill_values(features, config))
    return merged, missing, notes


def process_crispr_data(df, genes):
    df = df.rename(columns={"chrom": "chr", "chromStart": "start", "chromEnd": "end", "measuredGeneSymbol": "TargetGene"})
    df = df[~df["Regulated"].map(_missing)]
    return df[df["TargetGene"].isin(set(genes))].reset_index(drop=True)


def regulated_labels(s):
    np = _np()
    return to_bool(s).to_numpy().astype(np.int64)


# ---------------------------------------------------------------------------
# Training and feature analysis (train_model.py, training_functions.py, *_feature_selection.py,
# permutation_feature_importance.py, compare_all_feature_sets.py, compare_all_models.py)
# ---------------------------------------------------------------------------

def training_design(df, feature_list: "list[str]", epsilon: float, polynomial: bool):
    X = design_matrix(df, feature_list, epsilon, polynomial)
    return X.reset_index(drop=True)


def train_and_predict_once(chrs, X, Y, features: "list[str]", params: dict):
    """Leave-one-chromosome-out CV scores (or an in-sample fit when there is one chromosome)."""
    np = _np()
    Xf = X.loc[:, list(features)]
    scores = np.full(len(Y), np.nan)
    chr_list = np.unique(chrs)
    if len(chr_list) > 1:
        for c in chr_list:
            te = chrs == c
            if te.any():
                m = fit_lr(Xf[~te], Y[~te], params)
                scores[te] = m.predict_proba1(Xf[te])
    else:
        scores[:] = fit_lr(Xf, Y, params).predict_proba1(Xf)
    return scores


def train_and_predict(df, ft, params: dict, epsilon: float = EPSILON, polynomial: bool = False,
                      upstream_compat: bool = False) -> dict:
    """train_model.py::train_and_predict: full model, LOCO models, predictions, coefficients, metrics."""
    pd, np = _pd(), _np()
    df = df.reset_index(drop=True).copy()
    X = training_design(df, list(ft["feature"]), epsilon, polynomial)
    Y = regulated_labels(df["Regulated"])
    chrs = df["chr"].astype(str).to_numpy()
    full = fit_lr(X, Y, params)
    probs_full = full.predict_proba1(X)
    df[MODEL_NAME + ".Score_full"] = probs_full
    coef_rows = [pd.DataFrame({"feature": X.columns, "coefficient": full.coef, "test_chr": "none"})]
    metrics, models = [], {"full": full}
    chr_list = np.unique(chrs)
    if len(chr_list) > 1:
        cv = np.full(len(Y), np.nan)
        for c in chr_list:
            te = chrs == c
            m = fit_lr(X[~te], Y[~te], params)
            models[c] = m
            p = m.predict_proba1(X[te])
            cv[te] = p
            ptr = m.predict_proba1(X[~te])
            Yte, Ytr = Y[te], Y[~te]
            n_te_pos, n_tr_pos = int(Yte.sum()), int(Ytr.sum())
            ok = n_te_pos > 0
            both = ok and n_te_pos < len(Yte)
            metrics.append({"test_chr": c,
                            "log_loss_test_full": log_loss(Yte, probs_full[te]) if ok else np.nan,
                            "log_loss_train": log_loss(Ytr, ptr), "log_loss_test": log_loss(Yte, p) if ok else np.nan,
                            "AUROC_test_full": roc_auc(Yte, probs_full[te]) if both else np.nan,
                            "AUROC_train": roc_auc(Ytr, ptr), "AUROC_test": roc_auc(Yte, p) if both else np.nan,
                            "AUPRC_test_full": statistic_aupr(Yte, probs_full[te]) if ok else np.nan,
                            "AUPRC_train": statistic_aupr(Ytr, ptr), "AUPRC_test": statistic_aupr(Yte, p) if ok else np.nan,
                            "n_test_pos": n_te_pos, "n_test_neg": len(Yte) - n_te_pos,
                            "n_train_neg": len(Ytr) - n_tr_pos, "n_train_pos": n_tr_pos})
            coef_rows.append(pd.DataFrame({"feature": X.columns, "coefficient": m.coef, "test_chr": c}))
        df[SCORE_COL] = cv
        metrics.append({"test_chr": "all", "log_loss_test_full": log_loss(Y, probs_full), "log_loss_train": np.nan,
                        "log_loss_test": log_loss(Y, cv),
                        # upstream writes log_loss here; --upstream-compat reproduces it
                        "AUROC_test_full": log_loss(Y, probs_full) if upstream_compat else roc_auc(Y, probs_full),
                        "AUROC_train": np.nan, "AUROC_test": roc_auc(Y, cv),
                        "AUPRC_test_full": statistic_aupr(Y, probs_full), "AUPRC_train": np.nan,
                        "AUPRC_test": statistic_aupr(Y, cv), "n_test_pos": int(Y.sum()), "n_test_neg": int(len(Y) - Y.sum()),
                        "n_train_pos": np.nan, "n_train_neg": np.nan})
    return {"predictions": df, "coefficients": pd.concat(coef_rows, ignore_index=True),
            "metrics": pd.DataFrame(metrics, columns=METRIC_COLS), "models": models, "X": X, "Y": Y}


def _prep_analysis(df, ft, epsilon: float, polynomial: bool):
    pd, np = _pd(), _np()
    df = df.reset_index(drop=True)
    X = training_design(df, list(ft["feature"]), epsilon, polynomial)
    return X, regulated_labels(df["Regulated"]), df["chr"].astype(str).to_numpy(), list(X.columns)


def _eval(Y, s) -> "tuple[float, float]":
    return statistic_aupr(Y, s), statistic_precision_at_threshold(Y, s, threshold_70_pct_recall(Y, s))


def _boot_row(Y, y_last, y_new, n_boot: int, rng) -> dict:
    np = _np()
    t_last, t_new = threshold_70_pct_recall(Y, y_last), threshold_70_pct_recall(Y, y_new)
    rd_a = bootstrap((Y, y_last, y_new), statistic_delta_aupr, n_boot, rng)
    rd_p = bootstrap((Y, y_last, y_new), lambda a, b, c: statistic_delta_precision_at_threshold(a, b, c, t_last, t_new),
                     n_boot, rng)
    da, dp = float(np.mean(rd_a.bootstrap_distribution)), float(np.mean(rd_p.bootstrap_distribution))
    ra = bootstrap((Y, y_new), statistic_aupr, n_boot, rng)
    t = threshold_70_pct_recall(Y, y_new)
    rp = bootstrap((Y, y_new), lambda a, b: statistic_precision_at_threshold(a, b, t), n_boot, rng)
    return {"aupr": float(np.mean(ra.bootstrap_distribution)), "delta_aupr": da, "delta_aupr_low": rd_a.confidence_interval[0],
            "delta_aupr_high": rd_a.confidence_interval[1], "pval_aupr": bootstrap_pvalue(da, rd_a),
            "precision": float(np.mean(rp.bootstrap_distribution)), "delta_precision": dp,
            "delta_precision_low": rd_p.confidence_interval[0], "delta_precision_high": rd_p.confidence_interval[1],
            "pval_precision": bootstrap_pvalue(dp, rd_p)}


def _base_row(Y, y, n_boot: int, rng) -> dict:
    np = _np()
    ra = bootstrap((Y, y), statistic_aupr, n_boot, rng)
    t = threshold_70_pct_recall(Y, y)
    rp = bootstrap((Y, y), lambda a, b: statistic_precision_at_threshold(a, b, t), n_boot, rng)
    return {"aupr": float(np.mean(ra.bootstrap_distribution)), "delta_aupr": 0, "delta_aupr_low": 0, "delta_aupr_high": 0,
            "pval_aupr": 1, "precision": float(np.mean(rp.bootstrap_distribution)), "delta_precision": 0,
            "delta_precision_low": 0, "delta_precision_high": 0, "pval_precision": 1}


def sffs(X, Y, chrs, feature_list: "list[str]", params: dict) -> "tuple[list[str], list[float], list[float]]":
    """forward_sequential_feature_selection.py::SFFS (greedy; ties go to the later feature)."""
    remaining = list(feature_list)
    best, auprs, precs = [], [], []
    for _ in range(len(feature_list)):
        best_k, best_aupr, best_prec = -1, 0.0, 0.0
        for k, f in enumerate(remaining):
            a, p = _eval(Y, train_and_predict_once(chrs, X, Y, best + [f], params))
            if a >= best_aupr:
                best_k, best_aupr, best_prec = k, a, p
        best.append(remaining[best_k])
        auprs.append(best_aupr)
        precs.append(best_prec)
        remaining.pop(best_k)
    return best, auprs, precs


def sffs_significance(X, Y, chrs, feature_list, params, n_boot: int, rng):
    pd, np = _pd(), _np()
    order, _, _ = sffs(X, Y, chrs, feature_list, params)
    y_last = np.ones(len(Y))  # "all_true" baseline
    rows = [dict(feature_added="None", **_base_row(Y, y_last, n_boot, rng))]
    feats = []
    for f in order:
        feats.append(f)
        y_new = train_and_predict_once(chrs, X, Y, feats, params)
        rows.append(dict(feature_added=f, **_boot_row(Y, y_last, y_new, n_boot, rng)))
        y_last = y_new
    return pd.DataFrame(rows, columns=["feature_added"] + SEL_COLS), order


def sbfs(X, Y, chrs, feature_list: "list[str]", params: dict) -> "tuple[list[str], list[float], list[float]]":
    """backward_sequential_feature_selection.py::SBFS -> ['None', removed..., last remaining]."""
    fl = list(feature_list)
    a, p = _eval(Y, train_and_predict_once(chrs, X, Y, fl, params))
    removed, auprs, precs = ["None"], [a], [p]
    for _ in range(len(feature_list) - 1):
        best_k, best_aupr, best_prec = -1, 0.0, 0.0
        for k, f in enumerate(fl):
            a, p = _eval(Y, train_and_predict_once(chrs, X, Y, [g for g in fl if g != f], params))
            if a >= best_aupr:
                best_k, best_aupr, best_prec = k, a, p
        removed.append(fl[best_k])
        auprs.append(best_aupr)
        precs.append(best_prec)
        fl.pop(best_k)
    removed.append(fl[0])
    return removed, auprs, precs


def sbfs_significance(X, Y, chrs, feature_list, params, n_boot: int, rng):
    pd, np = _pd(), _np()
    order, _, _ = sbfs(X, Y, chrs, feature_list, params)
    fl = [f for f in order if f != "None"]
    y_last = train_and_predict_once(chrs, X, Y, fl, params)
    rows = [dict(feature_removed="None", **_base_row(Y, y_last, n_boot, rng))]
    for _ in range(len(fl)):
        to_remove = fl[0]
        if len(fl) == 1:
            y_new = np.ones(len(Y))
        else:
            fl.remove(to_remove)
            y_new = train_and_predict_once(chrs, X, Y, fl, params)
        rows.append(dict(feature_removed=to_remove, **_boot_row(Y, y_last, y_new, n_boot, rng)))
        y_last = y_new
    return pd.DataFrame(rows, columns=["feature_removed"] + SEL_COLS), order


def permutation_feature_importance(X, Y, chrs, feature_list, params, n_repeats: int, rng):
    pd = _pd()
    X = X.copy()
    y_full = train_and_predict_once(chrs, X, Y, feature_list, params)
    t_full = threshold_70_pct_recall(Y, y_full)
    rows = []
    for f in feature_list:
        original = X[f].copy()
        for _ in range(n_repeats):
            X[f] = rng.permutation(X[f].to_numpy())
            y_sh = train_and_predict_once(chrs, X, Y, feature_list, params)
            rows.append({"feature_permuted": f, "delta_aupr": statistic_delta_aupr(Y, y_full, y_sh),
                         "delta_precision": statistic_delta_precision_at_threshold(Y, y_full, y_sh, t_full,
                                                                                   threshold_70_pct_recall(Y, y_sh))})
        X[f] = original
    return pd.DataFrame(rows, columns=["feature_permuted", "delta_aupr", "delta_precision"])


def compare_feature_sets(X, Y, chrs, feature_list, params, n_boot: int, rng):
    pd, np = _pd(), _np()
    n = len(feature_list)
    rows = []
    for i in range(1, 2 ** n):
        bits = [int(b) for b in bin(i)[2:].zfill(n)]
        feats = [f for f, b in zip(feature_list, bits) if b]
        y = train_and_predict_once(chrs, X, Y, feats, params)
        r = bootstrap((Y, y), statistic_aupr, n_boot, rng)
        row = dict(zip(feature_list, bits))
        row.update({"features": str(feats), "n_features": len(feats), "AUPRC": float(np.mean(r.bootstrap_distribution)),
                    "AUPRC_95CI_low": r.confidence_interval[0], "AUPRC_95CI_high": r.confidence_interval[1]})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("AUPRC", ascending=False, kind="mergesort").reset_index(drop=True)


def performance_summary(model_id: str, dataset: str, y_true, y_pred, pct_missing: float, n_boot: int, rng) -> dict:
    np = _np()
    ra = bootstrap((y_true, y_pred), statistic_aupr, n_boot, rng)
    t = threshold_70_pct_recall(y_true, y_pred)
    if t is not None:
        rp = bootstrap((y_true, y_pred), lambda a, b: statistic_precision_at_threshold(a, b, t), n_boot, rng)
    return {"model": model_id, "dataset": dataset, "AUPRC": float(np.mean(ra.bootstrap_distribution)),
            "AUPRC_95CI_low": ra.confidence_interval[0], "AUPRC_95CI_high": ra.confidence_interval[1],
            "precision": 0 if t is None else float(np.mean(rp.bootstrap_distribution)),
            "precision_95CI_low": 0 if t is None else rp.confidence_interval[0],
            "precision_95CI_high": 0 if t is None else rp.confidence_interval[1],
            "threshold_70_pct_recall": t, "pct_missing_elements": pct_missing}


def distance_baseline(crispr):
    """compare_all_models.py 'distance' model: -|element centre - TSS centre| on the raw CRISPR table."""
    np = _np()
    c = crispr.copy()
    c["distance"] = ((to_float(c["chromStart"]) + to_float(c["chromEnd"])) / 2
                     - (to_float(c["startTSS"]) + to_float(c["endTSS"])) / 2).abs()
    c = c[~c["Regulated"].map(_missing) & c["distance"].notna()]
    return regulated_labels(c["Regulated"]), -c["distance"].to_numpy(dtype=np.float64)


# ---------------------------------------------------------------------------
# Figures (generate_plots.py, compare_plots.py, plot_sffs.R, plot_sbfs.R, plot_pfi.R, plot_model_comparison.R)
# ---------------------------------------------------------------------------

def _nice_names(ft, polynomial: bool) -> dict:
    if polynomial or ft is None:
        return {}
    return {f: n for f, n in zip(ft["feature"], ft["nice_name"]) if not _missing(n)}


def _stars(p) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p <= 0.05 else ""


def plot_selection(res, col: str, ft, polynomial: bool, d: Path, stem: str) -> "list[Path]":
    """plot_sffs.R / plot_sbfs.R: delta bars with 95% CI, AUPRC points, significance stars."""
    plt = _plt()
    if plt is None:
        return []
    names = _nice_names(ft, polynomial)
    labels = [names.get(x, x) for x in res[col]]
    out = []
    for metric, xl in (("aupr", "Delta AUPRC"), ("precision", "Delta precision at 70% recall")):
        fig, ax = plt.subplots(figsize=(5.5, max(2.5, 0.32 * len(res) + 1)))
        y = list(range(len(res)))[::-1]
        dv = res[f"delta_{metric}"].astype(float)
        lo, hi = res[f"delta_{metric}_low"].astype(float), res[f"delta_{metric}_high"].astype(float)
        ax.barh(y, dv, color=[ORANGE if v < 0 else BLUE for v in dv], height=0.6)
        ax.errorbar(dv, y, xerr=[(dv - lo).clip(lower=0).fillna(0), (hi - dv).clip(lower=0).fillna(0)], fmt="none", ecolor=INK2, lw=0.8)
        ax.scatter(res[metric].astype(float), y, s=12, color=INK, zorder=3, label=metric.upper() if metric == "aupr" else "precision")
        ref = res[metric].iloc[-1] if col == "feature_added" else res[metric].iloc[0]
        ax.axvline(float(ref), color=GREY, lw=0.8, ls=":")
        for yi, p, h in zip(y, res[f"pval_{metric}"], hi.fillna(0)):
            ax.text(max(float(h), 0) + 0.02, yi, _stars(float(p)), va="center", fontsize=7, color=INK)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(0, color=AXIS, lw=0.8)
        _style(ax, f"{'Forward' if col == 'feature_added' else 'Backward'} selection: {xl}", xl,
               "Feature added" if col == "feature_added" else "Feature removed")
        ax.legend(fontsize=7, frameon=False, loc="lower right")
        out.append(_save(fig, d / f"{stem}_{'auprc' if metric == 'aupr' else 'precision'}.png"))
    return out


def plot_pfi(res, ft, polynomial: bool, n_repeats: int, d: Path) -> "list[Path]":
    """plot_pfi.R: mean delta per permuted feature with t-based 95% CI."""
    plt, np = _plt(), _np()
    if plt is None:
        return []
    from scipy.stats import t as tdist  # type: ignore
    names = _nice_names(ft, polynomial)
    out = []
    for metric, xl in (("delta_aupr", "Delta AUPRC"), ("delta_precision", "Delta precision at 70% recall")):
        g = res.groupby("feature_permuted")[metric]
        s = g.agg(["mean", "std", "count"]).sort_values("mean", ascending=False)
        half = tdist.ppf(0.975, (s["count"] - 1).clip(lower=1)) * s["std"].fillna(0) / np.sqrt(s["count"])
        fig, ax = plt.subplots(figsize=(5.5, max(2.5, 0.32 * len(s) + 1)))
        y = list(range(len(s)))[::-1]
        ax.barh(y, s["mean"], xerr=half, color=BLUE, height=0.6, error_kw={"ecolor": INK2, "lw": 0.8})
        ax.set_yticks(y)
        ax.set_yticklabels([names.get(x, x) for x in s.index], fontsize=7)
        ax.axvline(0, color=AXIS, lw=0.8)
        _style(ax, f"Permutation feature importance (N={n_repeats})", xl, "Feature permuted")
        out.append(_save(fig, d / f"permutation_feature_importance_{'auprc' if metric == 'delta_aupr' else 'precision'}.png"))
    return out


def plot_model_comparison(df, d: Path) -> "list[Path]":
    plt = _plt()
    if plt is None:
        return []
    out = []
    for metric, lab in (("AUPRC", "AUPRC"), ("precision", "Precision at 70% recall")):
        s = df.sort_values(metric, ascending=True)
        lo_col, hi_col = (f"{metric}_95CI_low", f"{metric}_95CI_high")
        fig, ax = plt.subplots(figsize=(6, max(2.5, 0.35 * len(s) + 1)))
        y = range(len(s))
        v = s[metric].astype(float)
        err = [(v - s[lo_col].astype(float)).clip(lower=0).fillna(0), (s[hi_col].astype(float) - v).clip(lower=0).fillna(0)]
        ax.barh(list(y), v, xerr=err, color=[GREY if m == "distance" else BLUE for m in s["model"]], height=0.6,
                error_kw={"ecolor": INK2, "lw": 0.8})
        ax.set_yticks(list(y))
        ax.set_yticklabels([f"{m} ({ds})" for m, ds in zip(s["model"], s["dataset"])], fontsize=7)
        ax.set_xlim(0, 1)
        _style(ax, f"Performance across models: {lab}", lab, "Model")
        out.append(_save(fig, d / f"performance_across_models_{'auprc' if metric == 'AUPRC' else 'precision'}.png"))
    return out


def plot_pr(Y, scores: dict, d: Path, name: str, title: str) -> "list[Path]":
    plt = _plt()
    if plt is None:
        return []
    fig, ax = plt.subplots(figsize=(4.5, 4))
    for (lab, s), c in zip(scores.items(), [BLUE, ORANGE, GREY, INK2]):
        p, r, _ = pr_curve_modified(Y, s)
        ax.plot(r, p, color=c, lw=1.5, label=f"{lab} (AUPRC {statistic_aupr(Y, s):.3f})")
    ax.axhline(float(sum(Y)) / max(len(Y), 1), color=AXIS, lw=0.8, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=7, frameon=False)
    _style(ax, title, "Recall", "Precision")
    return [_save(fig, d / name)]


def plot_stats(stats: dict, d: Path, prefix: str = "") -> "list[Path]":
    """generate_plots.py: sequencing depth density, per-metric distribution and metric vs depth."""
    plt, np = _plt(), _np()
    if plt is None or not stats:
        return []
    out = []
    depth = np.array([s.get(SEQ_DEPTH_METRIC, 0.0) for s in stats.values()], dtype=float)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(depth / 1e6, bins=min(30, max(5, len(depth))), color=BLUE)
    _style(ax, f"{prefix}Num Sequencing Reads (n={len(depth)}, mean {depth.mean() / 1e6:.2f}M, median {np.median(depth) / 1e6:.2f}M)",
           "Millions of reads", "Datasets")
    out.append(_save(fig, d / f"{prefix.strip(': ').lower().replace(' ', '_') or 'all'}_num_sequencing_reads.png"))
    for metric in STATS_METRICS:
        pts = np.array([s.get(metric, np.nan) for s in stats.values()], dtype=float)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 3))
        ok = pts[np.isfinite(pts)]
        if ok.size:
            if ok.size > 1:
                a1.violinplot(ok, showmedians=True)
            a1.boxplot(ok, widths=0.15)
            a1.scatter(np.ones(ok.size) + np.random.default_rng(0).uniform(-0.05, 0.05, ok.size), ok, s=8, color=ORANGE, zorder=3)
        title = metric.replace("_", " ").title()
        _style(a1, f"{prefix}{title} (n={ok.size}, mean {np.nanmean(pts) if ok.size else float('nan'):.2f}, median "
                   f"{np.nanmedian(pts) if ok.size else float('nan'):.2f})", "", title)
        a2.scatter(depth / 1e6, pts, s=10, color=BLUE)
        _style(a2, f"{title} vs Sequencing Depth", "Sequencing Depth (Millions of Reads)", title)
        out.append(_save(fig, d / f"{prefix.strip(': ').lower().replace(' ', '_') or 'all'}_{metric}.png"))
    return out


def outlier_table(stats: dict):
    pd = _pd()
    rows = []
    for metric in STATS_METRICS:
        vals = sorted(((k, v.get(metric)) for k, v in stats.items() if v.get(metric) is not None), key=lambda kv: kv[1])
        for rank, (k, v) in enumerate(reversed(vals[-5:])):
            rows.append({"metric": metric, "position": "top", "rank": rank + 1, "dataset": k, "value": v})
        for rank, (k, v) in enumerate(vals[:5]):
            rows.append({"metric": metric, "position": "bottom", "rank": rank + 1, "dataset": k, "value": v})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_setup(args: argparse.Namespace) -> int:
    RES_DIR.mkdir(parents=True, exist_ok=True)
    for name, rel in RESOURCE_FILES.items():
        _download(RAW + rel, RES_DIR / name, args.force)
    if args.test_data:
        print("fetching the upstream chr22 expected output (~47 MB)")
        for name, rel in TEST_FILES.items():
            _download(RAW + rel, DATA_ROOT / "test_expected" / name, args.force)
    return 0


def _models_from_args(models: "Optional[list[str]]", tables: "Optional[list[str]]" = None):
    ms = [load_model(m) for m in (models or [])]
    fts = [m.feature_table for m in ms] + [feature_table_df(t) for t in (tables or [])]
    if not fts:
        raise SystemExit("give --model (embedded name or model dir) and/or --feature-table")
    return ms, combine_feature_tables(fts)


def cmd_models(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    rows, summ = [], {}
    exp = Path(args.export) if args.export else d / "models"
    for name in sorted(EMBEDDED_MODELS):
        m = embedded_model(name)
        rows.append([name, len(m.features), m.threshold, _fmt(m.intercept), ", ".join(m.features[:10]) + (" ..." if len(m.features) > 10 else "")])
        summ[name] = {"threshold": m.threshold, "n_features": len(m.features), "features": m.features,
                      "coef": dict(zip(m.features, [float(c) for c in m.coef])), "intercept": m.intercept}
        m.write_dir(exp / name, write_pickle=False)
    sec = [f"Pretrained ENCODE-rE2G models embedded from `{UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]}` (models/<name>/model.pkl "
           "coefficients). Each exported directory is a valid model dir for `apply --model DIR` (model.json replaces model.pkl).",
           md_table(["model", "n features", "threshold", "intercept", "features"], rows)]
    write_report(d, "ENCODE-rE2G pretrained models", sec, {"models": summ, "export_dir": str(exp)})
    return 0


def cmd_select_model(args: argparse.Namespace) -> int:
    pd = _pd()
    bs = read_tsv(args.biosample_config)
    sel = select_models(bs, args.model_root, args.megamap_hic)
    d = run_dir(args.label)
    out = pd.DataFrame([{"biosample": r["biosample"], "model_dir": r["model_dir"], "model_dir_base": r["model_dir_base"],
                         "model_threshold": r["model_threshold"]} for r in sel])
    write_tsv(out, d / "config_biosamples_models.tsv")
    write_report(d, "ENCODE-rE2G model selection",
                 ["Model folder = `<dhs|atac>[_h3k27ac]_<powerlaw|avg_hic|megamap|intact_hic>` (utils.smk), or the `model_dir` column.",
                  md_table(list(out.columns), out.itertuples(index=False))],
                 {"selections": out.to_dict("records")})
    return 0


def _abc_paths(abc_dir: "Optional[str]", pred: "Optional[str]", enh: "Optional[str]") -> "tuple[Path, Path]":
    if abc_dir:
        base = Path(abc_dir)
        pred = pred or str(base / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz")
        enh = enh or str(base / "Neighborhoods" / "EnhancerList.txt")
    if not pred or not enh:
        raise SystemExit("give --abc-dir (with Predictions/ and Neighborhoods/) or --abc-predictions and --enhancer-list")
    return Path(pred), Path(enh)


def _refs(args) -> "tuple[Path, Optional[Path], Optional[Path]]":
    tss = resource("RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed", args.tss)
    sizes = resource("GRCh38_EBV.no_alt.chrom.sizes.tsv", args.chr_sizes, required=False)
    gc = resource("gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv", args.gene_classes, required=False)
    return tss, sizes, gc


def cmd_features(args: argparse.Namespace) -> int:
    _, ft = _models_from_args(args.model, args.feature_table)
    pred, enh = _abc_paths(args.abc_dir, args.abc_predictions, args.enhancer_list)
    tss, sizes, gc = _refs(args)
    d = run_dir(args.label)
    bdir = d / safe_label(args.biosample or args.label)
    res = generate_features(pred, enh, tss, sizes, gc, ft, read_external_config(args.external_features_config), bdir,
                            dedupe_nearby=args.dedupe_nearby)
    final = res.pop("features")
    miss = {f: int(final[f].isna().sum()) for f in ft["feature"] if f in final.columns}
    sec = [f"Genome-wide ENCODE-rE2G features for `{pred}`: {res['n_pairs']:,} element-gene pairs "
           f"({res['n_enhancer_pairs']:,} non-promoter), {res['n_elements']:,} candidate elements.",
           md_table(["feature", "input_col", "second_input", "fill_value", "NA left"],
                    [[f, i, s, v, miss.get(f, "absent")] for f, i, s, v in zip(ft["feature"], ft["input_col"], ft["second_input"], ft["fill_value"])]),
           "External features: " + ("; ".join(res["external_notes"]) or "none")]
    write_report(d, "ENCODE-rE2G features", sec, dict(res, na_left=miss, dedupe_nearby=args.dedupe_nearby))
    return 0


def apply_model(features, model: Re2gModel, out_dir: Path, epsilon: float = EPSILON, include_self_promoter: bool = True,
                accessibility: "Optional[list[str]]" = None, use_pickle: bool = False,
                threshold: "Optional[float]" = None) -> dict:
    mdir = out_dir / model.name
    mdir.mkdir(parents=True, exist_ok=True)
    pred = make_e2g_predictions(features.copy(), model, epsilon, use_pickle)
    write_tsv(pred, mdir / "encode_e2g_predictions.tsv.gz")
    t = float(threshold if threshold is not None else model.threshold)
    thr = threshold_predictions(pred, t, SCORE_COL, include_self_promoter)
    ts = str(t)
    write_tsv(thr, mdir / f"encode_e2g_predictions_threshold{ts}.tsv.gz")
    bp = bedpe_frame(thr, SCORE_COL)
    write_tsv(bp, mdir / f"encode_e2g_predictions_threshold{ts}.bedpe", header=False)
    st = get_stats(thr, accessibility)
    write_tsv(st, mdir / f"encode_e2g_predictions_threshold{ts}_stats.tsv")
    return {"model": model.name, "source": model.source, "threshold": t, "n_pairs": int(len(pred)),
            "n_thresholded": int(len(thr)), "n_bedpe": int(len(bp)), "stats": dict(zip(st["Metric"], st["Value"])),
            "score_quantiles": {q: float(pred[SCORE_COL].quantile(q)) for q in (0.5, 0.9, 0.99)},
            "stats_path": str(mdir / f"encode_e2g_predictions_threshold{ts}_stats.tsv")}


def _apply_section(results: "list[dict]") -> str:
    return md_table(["model", "threshold", "pairs scored", "links >= threshold", "genes", "enhancers", "mean enh/gene"],
                    [[r["model"], r["threshold"], f"{r['n_pairs']:,}", f"{r['n_thresholded']:,}", int(r["stats"]["num_genes"]),
                      int(r["stats"]["num_enh"]), _fmt(r["stats"]["mean_num_enh_per_gene"], 2)] for r in results])


def cmd_apply(args: argparse.Namespace) -> int:
    feats = read_tsv(args.features)
    d = run_dir(args.label)
    bdir = d / safe_label(args.biosample or args.label)
    results = [apply_model(feats, load_model(m, args.use_pickle), bdir, args.epsilon, not args.exclude_self_promoter,
                           args.accessibility, args.use_pickle, args.threshold) for m in args.model]
    write_report(d, "ENCODE-rE2G predictions",
                 [f"Scored `{args.features}` with X = log(|x| + {args.epsilon}) and each model's logistic regression; "
                  f"thresholded with include_self_promoter={not args.exclude_self_promoter}.", _apply_section(results)],
                 {"features": args.features, "results": results})
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    bs = read_tsv(args.biosample_config)
    sel = select_models(bs, args.model_root, args.megamap_hic)
    tss, sizes, gc = _refs(args)
    d = run_dir(args.label)
    results, stats = [], {}
    by_bio: dict = {}
    for r in sel:
        by_bio.setdefault(r["biosample"], []).append(r)
    for bio, rs in by_bio.items():
        row = rs[0]["row"]
        abc_dir = row.get("ABC_directory") if not _missing(row.get("ABC_directory")) else (
            str(Path(args.abc_results) / bio) if args.abc_results else None)
        pred, enh = _abc_paths(abc_dir, None, None)
        ft = combine_feature_tables([x["model"].feature_table for x in rs])
        ext = read_external_config(row.get("external_features_config") if not _missing(row.get("external_features_config")) else None)
        bdir = d / safe_label(str(bio))
        res = generate_features(pred, enh, tss, sizes, gc, ft, ext, bdir, dedupe_nearby=args.dedupe_nearby)
        feats = res.pop("features")
        acc_col = str(row.get("default_accessibility_feature", ""))
        acc = [f.strip() for f in str(row.get(acc_col, "")).split(",") if f.strip() and not _missing(f)] if acc_col in row else []
        for x in rs:
            out = apply_model(feats, x["model"], bdir, args.epsilon, not args.exclude_self_promoter, acc)
            out["biosample"] = bio
            results.append(out)
            stats[f"{bio}/{x['model_dir_base']}"] = out["stats"]
    figs = [] if args.no_plots else plot_stats(stats, d)
    sec = [f"{len(by_bio)} biosample(s), {len(results)} model application(s).",
           md_table(["biosample", "model", "threshold", "links", "genes", "enhancers"],
                    [[r["biosample"], r["model"], r["threshold"], f"{r['n_thresholded']:,}", int(r["stats"]["num_genes"]),
                      int(r["stats"]["num_enh"])] for r in results])]
    write_report(d, "ENCODE-rE2G run", sec, {"results": results, "figures": [str(f) for f in figs]})
    return 0


def cmd_threshold(args: argparse.Namespace) -> int:
    df = read_tsv(args.predictions)
    d = run_dir(args.label)
    out = threshold_predictions(df, args.threshold, args.score_column, not args.exclude_self_promoter)
    p = write_tsv(out, d / f"encode_e2g_predictions_threshold{args.threshold}.tsv.gz")
    write_report(d, "ENCODE-rE2G threshold", [f"{len(out):,} of {len(df):,} pairs with {args.score_column} >= {args.threshold}."],
                 {"n_in": int(len(df)), "n_out": int(len(out)), "path": str(p)})
    return 0


def cmd_bedpe(args: argparse.Namespace) -> int:
    df = read_tsv(args.predictions)
    d = run_dir(args.label)
    bp = bedpe_frame(df, args.score_column)
    p = write_tsv(bp, d / (Path(args.predictions).name.replace(".tsv.gz", "").replace(".tsv", "") + ".bedpe"), header=False)
    write_report(d, "ENCODE-rE2G bedpe", [f"{len(bp):,} expressed-gene links written for IGV."], {"n": int(len(bp)), "path": str(p)})
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    df = read_tsv(args.predictions)
    d = run_dir(args.label)
    st = get_stats(df, args.accessibility)
    write_tsv(st, d / (Path(args.predictions).name.replace(".tsv.gz", "").replace(".tsv", "") + "_stats.tsv"))
    write_report(d, "ENCODE-rE2G prediction stats", [md_table(["Metric", "Value"], [[m, _fmt(v, 4)] for m, v in zip(st["Metric"], st["Value"])])],
                 {"stats": dict(zip(st["Metric"], st["Value"]))})
    return 0


def _load_stats(files: "list[str]") -> dict:
    out = {}
    for f in files:
        parts = Path(f).resolve().parts
        key = parts[-3] if len(parts) >= 3 else Path(f).stem
        if key in out:
            key = f"{key}/{parts[-2]}"
        df = read_tsv(f)
        out[key] = dict(zip(df["Metric"], df["Value"].astype(float)))
    return out


def cmd_qc_plots(args: argparse.Namespace) -> int:
    pd = _pd()
    stats = _load_stats(args.stats)
    d = run_dir(args.label)
    groups = {"": stats}
    if args.encode_metadata:
        meta = read_tsv(args.encode_metadata)
        cells, tissues = {}, {}
        for k, v in stats.items():
            acc = next((t for t in re.split(r"[_/]", k) if t.startswith("ENC")), None)
            m = meta[meta["DNase Experiment accession"] == acc]
            (tissues if len(m) and m.iloc[0]["Biosample type"] == "tissue" else cells)[k] = v
        groups = {"All: ": stats, "Cells: ": cells, "Tissues: ": tissues}
    figs = []
    if not args.no_plots:
        for pre, g in groups.items():
            figs += plot_stats(g, d, pre)
        if args.y2ave_metadata:
            figs += plot_y2ave(stats, read_tsv(args.y2ave_metadata), d)
    table = pd.DataFrame([dict(dataset=k, **v) for k, v in stats.items()])
    write_tsv(table, d / "stats_table.tsv")
    write_tsv(outlier_table(stats), d / "outlier_stats.tsv")
    rows = [[m, _fmt(table[m].mean(), 3), _fmt(table[m].median(), 3)] for m in [SEQ_DEPTH_METRIC] + STATS_METRICS if m in table]
    write_report(d, "ENCODE-rE2G QC", [f"{len(stats)} datasets; metrics are E-G pairs after the threshold (get_stats.py).",
                                        md_table(["metric", "mean", "median"], rows)],
                 {"n_datasets": len(stats), "figures": [str(f) for f in figs]})
    return 0


def plot_y2ave(stats: dict, meta, d: Path) -> "list[Path]":
    """generate_plots.py::plot_scatter: metric vs ATAC fragments (nCells x MeanATACFragmentsPerCell)."""
    plt = _plt()
    if plt is None:
        return []
    out = []
    for metric in STATS_METRICS:
        fig, ax = plt.subplots(figsize=(5, 3.5))
        for k, v in stats.items():
            m = meta[meta["CellClusterID"] == k]
            if not len(m):
                continue
            r = m.iloc[0]
            x = float(r["nCells"]) * float(r["MeanATACFragmentsPerCell"])
            ax.scatter([x], [v.get(metric)], s=12, color=ORANGE if "K562" in str(r["ManualAnnotationLabel"]) else BLUE)
        ax.set_xscale("log")
        ax.axvline(2e6, color=ORANGE, lw=0.8)
        _style(ax, metric.replace("_", " ").title(), "Num Fragments", metric)
        out.append(_save(fig, d / f"y2ave_{metric}.png"))
    return out


def cmd_compare_stats(args: argparse.Namespace) -> int:
    pd = _pd()
    rx = re.compile(args.key_regex)

    def keyed(files):
        out = {}
        for f in files:
            m = rx.search(str(f))
            if m:
                df = read_tsv(f)
                out[m.group(0)] = dict(zip(df["Metric"], df["Value"].astype(float)))
        return out
    a, b = keyed(args.a_stats), keyed(args.b_stats)
    keys = sorted(set(a) & set(b))
    d = run_dir(args.label)
    df = pd.DataFrame({"key": keys, args.a_name: [a[k].get(args.metric) for k in keys], args.b_name: [b[k].get(args.metric) for k in keys]})
    write_tsv(df, d / f"compare_{args.metric}.tsv")
    figs = []
    plt = _plt()
    if plt is not None and not args.no_plots and len(df):
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.scatter(df[args.a_name], df[args.b_name], s=12, color=BLUE, label=f"n={len(df)}")
        lim = [0, float(max(df[args.a_name].max(), df[args.b_name].max())) * 1.05]
        ax.plot(lim, lim, color=AXIS, lw=0.8)
        ax.legend(frameon=False, fontsize=7)
        _style(ax, args.metric, args.a_name, args.b_name)
        figs.append(_save(fig, d / f"compare_{args.metric}.png"))
    r = float(df[args.a_name].corr(df[args.b_name])) if len(df) > 2 else float("nan")
    write_report(d, "ENCODE-rE2G run comparison", [f"{len(df)} matched runs on `{args.key_regex}`; Pearson r = {_fmt(r)} for {args.metric}."],
                 {"n_matched": len(df), "pearson_r": r, "figures": [str(f) for f in figs]})
    return 0


def _feature_table_arg(args):
    if getattr(args, "feature_table", None):
        return feature_table_df(args.feature_table)
    if getattr(args, "model", None):
        return load_model(args.model).feature_table
    raise SystemExit("give --feature-table or --model")


def cmd_crispr_features(args: argparse.Namespace) -> int:
    ms, ft = _models_from_args(args.model, args.feature_table)
    feats = read_tsv(args.features)
    crispr = read_tsv(args.crispr)
    tss = read_bed(resource("RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed", args.tss),
                   ["chr", "start", "end", "name", "score", "strand"])
    merged, missing, notes = overlap_features_with_crispr(crispr, feats, ft, tss, nafill=not args.no_na_fill)
    processed = process_crispr_data(merged, tss["name"])
    d = run_dir(args.label)
    stem = Path(args.crispr).name.replace(".tsv.gz", "").replace(".tsv", "") + f".{safe_label(args.dataset or args.label)}_features_" + \
        ("NAfilled" if not args.no_na_fill else "NAnotfilled")
    write_tsv(merged, d / f"{stem}.tsv.gz")
    write_tsv(missing, d / f"missing.{stem}.tsv.gz")
    p = write_tsv(processed, d / f"for_training.{stem}.tsv.gz")
    write_tsv(ft, d / "feature_table.tsv")
    y = regulated_labels(processed["Regulated"]) if len(processed) else []
    summ = {"n_crispr": int(len(crispr)), "n_overlapping": int(len(merged)), "n_missing": int(len(missing)),
            "n_for_training": int(len(processed)), "n_positive": int(sum(y)), "notes": notes, "for_training": str(p),
            "missing": str(d / f"missing.{stem}.tsv.gz")}
    write_report(d, "ENCODE-rE2G CRISPR training data",
                 [f"{len(crispr):,} CRISPR E-G pairs -> {len(merged):,} overlapping feature pairs ({len(missing):,} missing) -> "
                  f"{len(processed):,} for training ({int(sum(y)):,} Regulated).", "\n".join(f"- {n}" for n in notes)], summ)
    return 0


def _params(args) -> dict:
    p = get_params(DEFAULT_PARAMS, args.override_params)
    p["verbose"] = 0
    return p


def cmd_train(args: argparse.Namespace) -> int:
    pd = _pd()
    df = read_tsv(args.crispr_features)
    ft = _feature_table_arg(args)
    params = _params(args)
    res = train_and_predict(df, ft, params, args.epsilon, args.polynomial, args.upstream_compat)
    d = run_dir(args.label)
    (d / "training_params.json").write_text(json.dumps(params, indent=2))
    with open(d / "training_params.pkl", "wb") as fh:
        pickle.dump(params, fh)
    for key, m in res["models"].items():
        if m.sk is not None:
            with open(d / ("model_full.pkl" if key == "full" else f"model_test_{key}.pkl"), "wb") as fh:
                pickle.dump(m.sk, fh)
    write_tsv(res["predictions"], d / "training_predictions.tsv")
    write_tsv(res["coefficients"], d / "model_coefficients.tsv")
    write_tsv(res["metrics"], d / "performance_metrics.tsv")
    Y = res["Y"]
    cv = res["predictions"][SCORE_COL].to_numpy() if SCORE_COL in res["predictions"] else res["predictions"][MODEL_NAME + ".Score_full"].to_numpy()
    thr_raw = threshold_70_pct_recall(Y, cv)
    thr = None if thr_raw is None else round(thr_raw, 3)  # 3 decimals, like the upstream threshold_ files
    full = res["models"]["full"]
    model = Re2gModel(safe_label(args.label), ft, full.coef, full.intercept, thr, None, source=str(d),
                      sk=full.sk, polynomial=args.polynomial, design_columns=list(res["X"].columns))
    model.write_dir(d / "model_dir", write_pickle=True)
    missing = args.missing
    if not missing:
        cand = Path(args.crispr_features).parent / Path(args.crispr_features).name.replace("for_training.", "missing.")
        missing = str(cand) if cand.exists() and cand.name != Path(args.crispr_features).name else None
    figs = [] if args.no_plots else plot_pr(Y, {"LOCO CV": cv, "full model (in-sample)": res["predictions"][MODEL_NAME + ".Score_full"].to_numpy()},
                                            d, "training_pr_curve.png", "ENCODE-rE2G training: CRISPR PR curve")
    allrow = res["metrics"][res["metrics"]["test_chr"] == "all"]
    coefs = res["coefficients"][res["coefficients"]["test_chr"] == "none"]
    summ = {"model": safe_label(args.label), "dataset": args.dataset or "dataset", "crispr_features": args.crispr_features,
            "missing": missing, "n": int(len(Y)), "n_positive": int(Y.sum()), "params": params, "polynomial": args.polynomial,
            "threshold_70_pct_recall_cv": thr_raw, "model_threshold": thr, "coefficients": dict(zip(coefs["feature"], coefs["coefficient"])),
            "intercept": full.intercept, "metrics_all": allrow.iloc[0].to_dict() if len(allrow) else {},
            "backend": "sklearn" if full.sk is not None else "numpy-newton", "figures": [str(f) for f in figs],
            "model_dir": str(d / "model_dir")}
    sec = [f"Logistic regression on {len(Y):,} CRISPR pairs ({int(Y.sum()):,} Regulated), X = log(|x| + {args.epsilon})"
           f"{' after degree-2 polynomial expansion' if args.polynomial else ''}; backend {summ['backend']}.",
           md_table(["feature", "coefficient"], [[f, _fmt(c, 4)] for f, c in zip(coefs["feature"], coefs["coefficient"])] + [["(intercept)", _fmt(full.intercept, 4)]]),
           md_table(["test_chr", "AUPRC_test", "AUROC_test", "log_loss_test", "n_test_pos", "n_test_neg"],
                    [[r.test_chr, _fmt(r.AUPRC_test), _fmt(r.AUROC_test), _fmt(r.log_loss_test), r.n_test_pos, r.n_test_neg]
                     for r in res["metrics"].itertuples()]),
           f"Threshold at 70% recall of the CV scores: {_fmt(thr_raw, 4)} (written, 3 decimals, as `model_dir/threshold_{format_threshold(thr) if thr is not None else ''}`)."]
    write_report(d, "ENCODE-rE2G model training", sec, summ)
    return 0


def _analysis_inputs(args):
    df = read_tsv(args.crispr_features)
    ft = _feature_table_arg(args)
    X, Y, chrs, fl = _prep_analysis(df, ft, args.epsilon, args.polynomial)
    return ft, X, Y, chrs, fl, _params(args), _np().random.default_rng(args.seed)


def cmd_feature_selection(args: argparse.Namespace) -> int:
    ft, X, Y, chrs, fl, params, rng = _analysis_inputs(args)
    if args.direction == "forward":
        res, order = sffs_significance(X, Y, chrs, fl, params, args.n_boot, rng)
        col, stem = "feature_added", "forward_feature_selection"
    else:
        res, order = sbfs_significance(X, Y, chrs, fl, params, args.n_boot, rng)
        col, stem = "feature_removed", "backward_feature_selection"
    d = run_dir(args.label)
    write_tsv(res, d / f"{stem}.tsv")
    figs = [] if args.no_plots else plot_selection(res, col, ft, args.polynomial, d, stem)
    write_report(d, f"ENCODE-rE2G {args.direction} sequential feature selection",
                 [f"{len(fl)} features, {len(Y):,} CRISPR pairs, leave-one-chromosome-out CV, {args.n_boot} BCa bootstrap resamples.",
                  md_table([col, "aupr", "delta_aupr", "95% CI", "p", "precision", "delta_precision"],
                           [[r[col], _fmt(r["aupr"]), _fmt(r["delta_aupr"]), f"[{_fmt(r['delta_aupr_low'])}, {_fmt(r['delta_aupr_high'])}]",
                             _fmt(r["pval_aupr"]), _fmt(r["precision"]), _fmt(r["delta_precision"])] for _, r in res.iterrows()])],
                 {"direction": args.direction, "order": order, "table": res.to_dict("records"), "figures": [str(f) for f in figs]})
    return 0


def cmd_permutation_importance(args: argparse.Namespace) -> int:
    ft, X, Y, chrs, fl, params, rng = _analysis_inputs(args)
    res = permutation_feature_importance(X, Y, chrs, fl, params, args.n_repeats, rng)
    d = run_dir(args.label)
    write_tsv(res, d / "permutation_feature_importance.tsv")
    figs = [] if args.no_plots else plot_pfi(res, ft, args.polynomial, args.n_repeats, d)
    s = res.groupby("feature_permuted")[["delta_aupr", "delta_precision"]].mean().sort_values("delta_aupr")
    write_report(d, "ENCODE-rE2G permutation feature importance",
                 [f"{args.n_repeats} permutations per feature; delta = metric(shuffled) - metric(full), LOCO CV.",
                  md_table(["feature", "mean delta AUPRC", "mean delta precision"], [[i, _fmt(r.delta_aupr, 4), _fmt(r.delta_precision, 4)] for i, r in s.iterrows()])],
                 {"mean_delta": s.to_dict("index"), "figures": [str(f) for f in figs]})
    return 0


def cmd_all_feature_sets(args: argparse.Namespace) -> int:
    args.polynomial = False
    ft, X, Y, chrs, fl, params, rng = _analysis_inputs(args)
    if len(fl) >= 14 and not args.force:
        raise SystemExit(f"{len(fl)} features -> {2 ** len(fl) - 1} models; upstream only runs this for < 14 features (use --force)")
    res = compare_feature_sets(X, Y, chrs, fl, params, args.n_boot, rng)
    d = run_dir(args.label)
    write_tsv(res, d / "all_feature_sets.tsv")
    write_report(d, "ENCODE-rE2G all feature sets",
                 [f"{len(res)} feature subsets, bootstrapped AUPRC (LOCO CV).",
                  md_table(["features", "n", "AUPRC", "95% CI"], [[r["features"], r["n_features"], _fmt(r["AUPRC"]),
                                                                    f"[{_fmt(r['AUPRC_95CI_low'])}, {_fmt(r['AUPRC_95CI_high'])}]"] for _, r in res.head(20).iterrows()])],
                 {"best": res.iloc[0].to_dict(), "n_sets": int(len(res))})
    return 0


def cmd_compare_models(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    rng = np.random.default_rng(args.seed)
    rows = []
    for td in args.train_dirs:
        td = Path(td)
        s = json.loads((td / "summary.json").read_text()) if (td / "summary.json").exists() else {}
        pred = read_tsv(td / "training_predictions.tsv")
        col = SCORE_COL if SCORE_COL in pred.columns else MODEL_NAME + ".Score_full"
        y, p = regulated_labels(pred["Regulated"]), to_float(pred[col]).to_numpy()
        n_miss = 0
        if s.get("missing") and Path(s["missing"]).exists():
            miss = read_tsv(s["missing"])
            miss = miss[~miss["Regulated"].map(_missing)]
            n_miss = len(miss)
            y = np.concatenate([y, regulated_labels(miss["Regulated"])])
            p = np.concatenate([p, np.zeros(n_miss)])
        rows.append(performance_summary(s.get("model", td.name), s.get("dataset", "dataset"), y, p, n_miss / max(len(y), 1), args.n_boot, rng))
    if args.crispr:
        y, p = distance_baseline(read_tsv(args.crispr))
        rows.append(performance_summary("distance", "baseline", y, p, 0.0, args.n_boot, rng))
    df = pd.DataFrame(rows).sort_values("AUPRC", ascending=False, kind="mergesort").reset_index(drop=True)
    d = run_dir(args.label)
    write_tsv(df, d / "performance_across_models.tsv")
    figs = [] if args.no_plots else plot_model_comparison(df, d)
    write_report(d, "ENCODE-rE2G model comparison",
                 ["CV performance on the CRISPR training data; CRISPR pairs without features score 0 (compare_all_models.py).",
                  md_table(["model", "dataset", "AUPRC", "95% CI", "precision@70% recall", "threshold", "% missing"],
                           [[r.model, r.dataset, _fmt(r.AUPRC), f"[{_fmt(r.AUPRC_95CI_low)}, {_fmt(r.AUPRC_95CI_high)}]", _fmt(r.precision),
                             _fmt(r.threshold_70_pct_recall, 4), _fmt(100 * r.pct_missing_elements, 1)] for r in df.itertuples()])],
                 {"table": df.to_dict("records"), "figures": [str(f) for f in figs]})
    return 0


def verify_upstream(expected_dir: Path, upstream_dir: "Optional[Path]" = None) -> dict:
    """Re-score the upstream K562 chr22 expected output with the embedded dhs_intact_hic coefficients."""
    np = _np()
    out: dict = {"expected_dir": str(expected_dir)}
    pred = read_tsv(expected_dir / "encode_e2g_predictions.tsv.gz")
    m = embedded_model("dhs_intact_hic")
    ours = m.predict_proba(pred.replace([np.inf, -np.inf], np.nan).fillna(0), EPSILON)
    diff = np.abs(ours - pred[SCORE_COL].to_numpy())
    out["n_pairs"] = int(len(pred))
    out["max_abs_score_diff"] = float(diff.max())
    thr = threshold_predictions(pred.assign(**{SCORE_COL: ours}), m.threshold, SCORE_COL, True)
    exp_thr = read_tsv(expected_dir / "encode_e2g_predictions_threshold0.243.tsv.gz")
    k1 = set(zip(thr["name"], thr["TargetGene"]))
    k2 = set(zip(exp_thr["name"], exp_thr["TargetGene"]))
    out.update({"n_thresholded": len(thr), "n_thresholded_expected": len(exp_thr), "thresholded_pairs_identical": k1 == k2})
    st = get_stats(thr)
    exp_st = read_tsv(expected_dir / "encode_e2g_predictions_threshold0.243_stats.tsv")
    exp_d = dict(zip(exp_st["Metric"], exp_st["Value"]))
    out["stats"] = {k: {"ours": float(v), "expected": float(exp_d[k])} for k, v in zip(st["Metric"], st["Value"]) if k != SEQ_DEPTH_METRIC}
    out["stats_max_rel_diff"] = max(abs(v["ours"] - v["expected"]) / max(abs(v["expected"]), 1e-12) for v in out["stats"].values())
    if upstream_dir and (upstream_dir / "models").is_dir():
        pk = {}
        for md in sorted((upstream_dir / "models").iterdir()):
            if not (md / "model.pkl").exists() or md.name not in EMBEDDED_MODELS:
                continue
            try:
                um = load_model(str(md), prefer_pickle=True)
            except SystemExit as e:
                pk[md.name] = {"error": str(e)}
                continue
            em = embedded_model(md.name)
            rec = {"coef_identical": bool(np.array_equal(um.coef, em.coef)) and um.intercept == em.intercept,
                   "threshold_identical": um.threshold == em.threshold, "features_identical": um.features == em.features}
            if md.name == "dhs_intact_hic" and um.sk is not None:
                sub = pred.head(200000).replace([np.inf, -np.inf], np.nan).fillna(0)
                rec["max_abs_pickle_vs_numpy"] = float(np.max(np.abs(um.predict_proba(sub, EPSILON, use_pickle=True) - em.predict_proba(sub, EPSILON))))
            pk[md.name] = rec
        out["pickles"] = pk
    return out


def cmd_verify_upstream(args: argparse.Namespace) -> int:
    up = Path(args.upstream_dir) if args.upstream_dir else None
    exp = Path(args.expected_dir) if args.expected_dir else (
        up / "tests/expected_output/generic/K562_chr22/dhs_intact_hic" if up else DATA_ROOT / "test_expected")
    if not (exp / "encode_e2g_predictions.tsv.gz").exists():
        raise SystemExit(f"{exp}: expected output not found (run `igvfagent encode-re2g setup --test-data` or pass --upstream-dir)")
    v = verify_upstream(exp, up)
    d = run_dir(args.label)
    rows = [[k, _fmt(x["ours"], 6), _fmt(x["expected"], 6)] for k, x in v["stats"].items()]
    sec = [f"Upstream CircleCI expected output (K562 chr22, dhs_intact_hic, {v['n_pairs']:,} pairs) re-scored with the embedded "
           f"coefficients: max |score difference| = {v['max_abs_score_diff']:.3g}; thresholded pairs {v['n_thresholded']:,} vs "
           f"{v['n_thresholded_expected']:,} expected, identical = {v['thresholded_pairs_identical']}.",
           md_table(["metric", "ours", "expected"], rows)]
    if "pickles" in v:
        sec.append(md_table(["model", "coef identical", "threshold identical", "features identical", "max |pickle - numpy|"],
                            [[k, r.get("coef_identical"), r.get("threshold_identical"), r.get("features_identical"),
                              r.get("max_abs_pickle_vs_numpy", "")] for k, r in v["pickles"].items()]))
    write_report(d, "ENCODE-rE2G upstream verification", sec, v)
    ok = v["max_abs_score_diff"] < 1e-9 and v["thresholded_pairs_identical"] and v["stats_max_rel_diff"] < 1e-9
    print(f"verify-upstream: {'PASS' if ok else 'MISMATCH'} (max score diff {v['max_abs_score_diff']:.2e}, stats rel diff {v['stats_max_rel_diff']:.2e})")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Self-test: synthetic ABC world with a planted ENCODE-rE2G model
# ---------------------------------------------------------------------------

CUSTOM_FT = [("noiseFeature", "noise_src", "NA", "mean", "0", "Noise (external, overlap)"),
             ("geneNoise", "gene_noise", "NA", "max", "mean", "Gene noise (external, TargetGene)"),
             ("ABCxContact", "ABC.Score", "hic_contact_pl_scaled_adj", "mean", "0", "ABC x contact"),
             ("numNearbyEnhancers", "numNearbyEnhancers", "NA", "max", "0", "# peaks within 5Kb of E"),
             ("P2PromoterClass", "P2PromoterClass", "NA", "max", "0", "Promoter class")]


def synthetic_world(d: Path, seed: int = 5) -> dict:
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    chroms = {"chr1": 3_000_000, "chr2": 2_500_000, "chr3": 2_500_000}
    genes, elements = [], []
    for c, size in chroms.items():
        tss_pos = np.sort(rng.choice(np.arange(200_000, size - 200_000, 1000), 14, replace=False))
        for i, t in enumerate(tss_pos):
            genes.append({"chr": c, "gene": f"G{c[3:]}_{i}", "tss": int(t)})
        pos, k = 20_000, 0
        while pos < size - 20_000:
            w = int(rng.integers(200, 900))
            if all(abs(pos - t) > 2_000 and abs(pos + w - t) > 2_000 for t in tss_pos):
                elements.append({"chr": c, "start": pos, "end": pos + w, "class": "genic" if rng.random() < 0.4 else "intergenic"})
                k += 1
            pos += w + int(rng.integers(600, 45_000))
        for t in tss_pos:
            elements.append({"chr": c, "start": int(t) - 250, "end": int(t) + 250, "class": "promoter"})
    el = pd.DataFrame(elements).sort_values(["chr", "start"]).reset_index(drop=True)
    el["name"] = el["class"] + "|" + el["chr"] + ":" + el["start"].astype(str) + "-" + el["end"].astype(str)
    el["activity_base"] = np.round(rng.lognormal(0.5, 1.0, len(el)), 5)
    g = pd.DataFrame(genes)
    g["prom"] = np.round(rng.lognormal(0, 0.6, len(g)), 5)
    g["q"] = np.round(rng.random(len(g)), 4)
    g["expressed"] = rng.random(len(g)) < 0.85
    g["ubiq"] = rng.random(len(g)) < 0.3
    rows = []
    for gr in g.itertuples(index=False):
        sub = el[(el["chr"] == gr.chr)]
        mid = ((sub["start"] + sub["end"]) / 2).astype(np.int64)
        near = sub[np.abs(mid - gr.tss) <= 1_000_000]
        for e in near.itertuples(index=False):
            m = int((e.start + e.end) / 2)
            dist = abs(m - gr.tss)
            rows.append({"chr": e.chr, "start": e.start, "end": e.end, "name": e.name, "class": e[3], "activity_base": e.activity_base,
                         "TargetGene": gr.gene, "TargetGeneTSS": gr.tss, "TargetGeneIsExpressed": bool(gr.expressed),
                         "TargetGeneEnsembl_ID": "ENSG" + gr.gene, "isSelfPromoter": bool(e[3] == "promoter" and e.start == gr.tss - 250),
                         "CellType": "SYN", "distance": float(dist),
                         "hic_contact_pl_scaled_adj": float(np.round(100.0 / (dist + 5000) * rng.lognormal(0, 0.3), 8)),
                         "normalized_dhs_prom": gr.prom, "TargetGenePromoterActivityQuantile": gr.q})
    pred = pd.DataFrame(rows)
    pred["ABC.Score.Numerator"] = pred["activity_base"] * pred["hic_contact_pl_scaled_adj"]
    pred["ABC.Score"] = pred["ABC.Score.Numerator"] / pred.groupby("TargetGene")["ABC.Score.Numerator"].transform("sum")
    abc = d / "abc" / "SYN"
    (abc / "Predictions").mkdir(parents=True)
    (abc / "Neighborhoods").mkdir(parents=True)
    pred.to_csv(abc / "Predictions" / "EnhancerPredictionsAllPutative.tsv.gz", sep="\t", index=False)
    el[["chr", "start", "end", "name", "class", "activity_base"]].to_csv(abc / "Neighborhoods" / "EnhancerList.txt", sep="\t", index=False)
    tss = pd.DataFrame({"chr": g["chr"], "start": g["tss"] - 250, "end": g["tss"] + 250, "name": g["gene"], "score": 0, "strand": "+"})
    tss.to_csv(d / "tss500.bed", sep="\t", header=False, index=False)
    pd.DataFrame(list(chroms.items())).to_csv(d / "sizes.tsv", sep="\t", header=False, index=False)
    pd.DataFrame({"TargetGene": g["gene"], "gene_chr": g["chr"], "gene_tss": g["tss"], "gene_tss500bp": g["tss"] + 250,
                  "is_ubiquitous_uniform": g["ubiq"], "P2PromoterClass": rng.random(len(g)) < 0.2}).to_csv(d / "gene_classes.tsv", sep="\t", index=False)
    # external features: per-pair noise by overlap (10% missing; two overlapping sources for 20 pairs), per-gene noise (80% of genes)
    npf = pred[["chr", "start", "end", "TargetGene"]].copy()
    npf["noise_val"] = np.round(rng.lognormal(0, 1, len(npf)), 5)
    keep = rng.random(len(npf)) > 0.1
    dup = npf[keep].sample(20, random_state=1).copy()
    dup["noise_val"] = np.round(rng.lognormal(0, 1, len(dup)), 5)
    dup["start"] = dup["start"] + 10
    noise = pd.concat([npf[keep], dup], ignore_index=True)
    noise.to_csv(d / "noise.tsv", sep="\t", index=False)
    gsub = g.sample(frac=0.8, random_state=2)
    pd.DataFrame({"TargetGene": gsub["gene"], "gnoise": np.round(rng.lognormal(0, 1, len(gsub)), 5)}).to_csv(d / "gene_noise.tsv", sep="\t", index=False)
    pd.DataFrame({"input_col": ["noise_src", "gene_noise"], "source_col": ["noise_val", "gnoise"], "aggregate_function": ["mean", "max"],
                  "join_by": ["overlap", "TargetGene"], "source_file": ["noise.tsv", "gene_noise.tsv"]}).to_csv(d / "external.tsv", sep="\t", index=False)
    write_tsv(feature_table_df(CUSTOM_FT), d / "custom_ft.tsv", quiet=True)
    return {"abc": abc, "pred": pred, "el": el, "genes": g, "tss": d / "tss500.bed", "sizes": d / "sizes.tsv", "chroms": chroms,
            "gc": d / "gene_classes.tsv", "ext": d / "external.tsv", "custom_ft": d / "custom_ft.tsv", "noise": noise, "rng": rng}


def plant_crispr(d: Path, W: dict, feats) -> dict:
    """CRISPR-like truth drawn from the embedded dhs_intact_hic model (sharpened, recentred) + edge cases."""
    pd, np = _pd(), _np()
    rng = np.random.default_rng(17)
    m = embedded_model("dhs_intact_hic")
    cand = feats[(feats["class"] != "promoter") & (to_float(feats["distanceToTSS"]) < 600_000)].reset_index(drop=True)
    L = m.design(cand.replace([np.inf, -np.inf], np.nan).fillna(0)).to_numpy() @ m.coef + m.intercept
    lo, hi = -30.0, 30.0
    for _ in range(80):
        shift = (lo + hi) / 2
        lo, hi = (shift, hi) if _expit(1.5 * L + shift).mean() < 0.15 else (lo, shift)
    p = _expit(1.5 * L + shift)
    idx = rng.choice(len(cand), min(1100, len(cand)), replace=False)
    s = cand.iloc[idx].reset_index(drop=True)
    reg = rng.random(len(s)) < p[idx]
    cr = pd.DataFrame({"dataset": "SynthScreen", "chrom": s["chr"], "chromStart": s["start"] + rng.integers(0, 50, len(s)),
                       "chromEnd": s["end"] - rng.integers(0, 50, len(s)), "name": s["TargetGene"] + "|" + s["name"],
                       "EffectSize": np.where(reg, -0.3, 0.0), "chrTSS": s["chr"], "startTSS": s["TargetGeneTSS"], "endTSS": s["TargetGeneTSS"] + 1,
                       "measuredGeneSymbol": s["TargetGene"], "Regulated": reg})
    genes = W["genes"].set_index("gene")
    # element spanning two consecutive non-promoter candidates of one gene -> aggregation
    fp = feats[(feats["class"] != "promoter")].sort_values(["TargetGene", "start"])
    for gname, sub in fp.groupby("TargetGene"):
        if len(sub) >= 4:
            e1, e2 = sub.iloc[1], sub.iloc[2]
            break
    span = {"dataset": "SynthScreen", "chrom": e1["chr"], "chromStart": int(e1["start"]), "chromEnd": int(e2["end"]),
            "name": "span", "EffectSize": -0.5, "chrTSS": e1["chr"], "startTSS": int(genes.loc[gname, "tss"]),
            "endTSS": int(genes.loc[gname, "tss"]) + 1, "measuredGeneSymbol": gname, "Regulated": True}
    e3 = fp[fp["TargetGene"] == gname].iloc[3]
    touch = dict(span, name="touch", chromStart=int(e3["start"]) - 300, chromEnd=int(e3["start"]), Regulated=False)
    extra = [span, touch]
    # CRISPR elements in element-free stretches -> 'missing'
    el = W["el"]
    n_gap = 0
    for c in W["chroms"]:
        ec = el[el["chr"] == c].sort_values("start")
        gaps = ec["start"].to_numpy()[1:] - ec["end"].to_numpy()[:-1]
        for j in np.flatnonzero(gaps > 5000)[:4]:
            s0 = int(ec["end"].to_numpy()[j]) + 2000
            gg = W["genes"][W["genes"]["chr"] == c].iloc[0]
            extra.append({"dataset": "SynthScreen", "chrom": c, "chromStart": s0, "chromEnd": s0 + 500, "name": f"gap{n_gap}",
                          "EffectSize": 0.0, "chrTSS": c, "startTSS": int(gg["tss"]), "endTSS": int(gg["tss"]) + 1,
                          "measuredGeneSymbol": gg["gene"], "Regulated": bool(n_gap % 3 == 0)})
            n_gap += 1
    for k in range(4):
        extra.append(dict(span, name=f"nogene{k}", measuredGeneSymbol="NOT_IN_UNIVERSE"))
    cr = pd.concat([cr, pd.DataFrame(extra)], ignore_index=True)
    cr["pair_uid"] = cr["dataset"] + "|" + cr["name"]
    cr["merged_uid"] = range(len(cr))
    cr["merged_start"], cr["merged_end"] = cr["chromStart"], cr["chromEnd"]
    path = d / "crispr.tsv.gz"
    cr.to_csv(path, sep="\t", index=False)
    return {"path": path, "n_gap": n_gap, "n_sampled": len(s), "span": (e1, e2, gname), "touch": (e3, gname),
            "planted_scale": 1.5, "planted_shift": shift, "base_rate": float(reg.mean())}


def _parse(argv):
    args = build_parser().parse_args(argv)
    return args.func(args)


def _last_run(label: str) -> Path:
    return sorted(OUT_ROOT.glob(f"*_{label}"))[-1]


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    pd, np = _pd(), _np()
    checks = []
    np_ = "--no-plots" if args.no_plots else None

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def flags(*a):
        return [x for x in a if x is not None]

    existed = OUT_ROOT.exists()
    before = set(OUT_ROOT.glob("*")) if existed else set()
    sk = _sklearn_ok()
    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            print("\nmetrics")
            rng = np.random.default_rng(3)
            y = (rng.random(300) < 0.3).astype(int)
            s = np.round(rng.random(300) + 0.4 * y, 2)  # ties on purpose
            if sk:
                from sklearn.metrics import precision_recall_curve as skpr, auc as skauc, roc_auc_score, log_loss as skll  # type: ignore
                from sklearn.preprocessing import PolynomialFeatures  # type: ignore
                p1, r1, t1 = precision_recall_curve(y, s)
                p2, r2, t2 = skpr(y, s)
                check(np.allclose(p1, p2) and np.allclose(r1, r2) and np.allclose(t1, t2), "precision_recall_curve == sklearn (with ties)")
                check(abs(statistic_aupr(y, s) - skauc(r2[1:], p2[1:])) < 1e-12 and abs(roc_auc(y, s) - roc_auc_score(y, s)) < 1e-12
                      and abs(log_loss(y, np.clip(s, 0, 1)) - skll(y, np.clip(s, 0, 1))) < 1e-9,
                      "AUPRC (first point dropped), AUROC and log loss == sklearn")
                Xs = pd.DataFrame(rng.random((5, 3)), columns=["a", "b", "c"])
                pf = PolynomialFeatures(degree=2)
                arr = pf.fit_transform(Xs)
                mine = poly2(Xs)
                check(list(mine.columns) == list(pf.get_feature_names_out(Xs.columns)) and np.allclose(mine.to_numpy(), arr),
                      "poly2 == sklearn PolynomialFeatures(degree=2) values and names")
            check(statistic_aupr([1, 0], [1.0, 1.0]) == 0.0 and threshold_70_pct_recall([1, 1, 1], [0.1, 0.2, 0.3]) is None,
                  "AUPRC of a single-point curve is 0; no 70%-recall threshold when max recall <= 0.7 after the dropped point")
            t = threshold_70_pct_recall(y, s)
            pt = statistic_precision_at_threshold(y, s, t)
            pp, rr, tt = pr_curve_modified(y, s)
            i = int(np.argmin(np.abs(rr - 0.7)))
            check(t == tt[i + 1] and pt == pp[int(np.argmin(np.abs(tt - t))) - 1],
                  f"threshold at 70% recall = thresholds[i+1] ({t}); precision there uses the upstream j-1 index ({pt:.3f})")
            gp = get_params(DEFAULT_PARAMS, "{'C': '0.5', 'penalty': 'l2', 'max_iter': '100.0', 'class_weight': 'null'}")
            check(gp["C"] == 0.5 and gp["max_iter"] == 100 and isinstance(gp["max_iter"], int) and gp["penalty"] == "l2" and gp["class_weight"] is None
                  and gp["solver"] == "lbfgs", "get_params: overrides, numeric strings -> numbers (max_iter int), 'null' -> None")
            cft = combine_feature_tables([embedded_model("dhs_megamap").feature_table, embedded_model("dhs_intact_hic").feature_table,
                                          feature_table_df([("Kendall", "Kendall", "NA", "max", "0", "k")])])
            check(len(cft) == 8 + 1 + 5, f"combine_feature_tables: rbind + distinct + sc-E2G ARC rows ({len(cft)} rows)")

            print("\nmodels")
            check(len(EMBEDDED_MODELS) == 9 and all(len(v["coef"]) == len(v["feature_table"]) for v in EMBEDDED_MODELS.values()),
                  "9 embedded pretrained models, one coefficient per feature_table row")
            bs = pd.DataFrame({"biosample": ["a", "b", "c", "d", "e"], "default_accessibility_feature": ["DHS", "ATAC", "DHS", "ATAC", "DHS"],
                               "H3K27ac": [None, "x.bam", None, None, None],
                               "HiC_file": ["https://x/ENCFF621AIY.hic", MEGAMAP_HIC_FILE, None, "y.hic", "z.hic"],
                               "HiC_type": ["hic", "hic", None, "avg", "hic"], "model_dir": [None, None, None, None, "dhs_megamap"]})
            names = [model_name_for_row(r) for r in bs.to_dict("records")]
            check(names[:4] == ["dhs_intact_hic", "atac_h3k27ac_megamap", "dhs_powerlaw", "atac_avg_hic"],
                  f"model choice from biosample config: {names[:4]}")
            sel = select_models(bs.iloc[[0, 1, 4]])
            check([r["model_dir_base"] for r in sel] == ["dhs_intact_hic", "atac_h3k27ac_megamap", "dhs_megamap"]
                  and [r["model_threshold"] for r in sel] == [0.243, 0.203, 0.201], "select_models: thresholds 0.243 / 0.203; model_dir overrides")
            try:
                select_models(bs.iloc[[2]])
                check(False, "powerlaw model should be rejected")
            except SystemExit as e:
                check("not supported" in str(e), "select_models: untrained powerlaw model rejected (as upstream)")
            em = embedded_model("dhs_intact_hic")
            toy = pd.DataFrame({f: np.abs(rng.normal(1, 2, 50)) for f in em.features})
            manual = 1 / (1 + np.exp(-(np.log(np.abs(toy[em.features].to_numpy()) + 0.01) @ em.coef + em.intercept)))
            check(np.allclose(em.predict_proba(toy), manual, atol=0, rtol=1e-14), "embedded model score = 1/(1+exp(-(log(|x|+0.01).coef + b)))")
            if sk:
                from sklearn.linear_model import LogisticRegression  # type: ignore
                lr = LogisticRegression()
                lr.coef_ = em.coef.reshape(1, -1).copy()
                lr.intercept_ = np.array([em.intercept])
                lr.classes_ = np.array([0, 1])
                lr.feature_names_in_ = np.array(em.features, dtype=object)
                lr.n_features_in_ = len(em.features)
                mdir = d / "pickled_model"
                em.write_dir(mdir)
                (mdir / "model.json").unlink()
                with open(mdir / "model.pkl", "wb") as fh:
                    pickle.dump(lr, fh)
                lm = load_model(str(mdir))
                check(lm.sk is not None and np.max(np.abs(lm.predict_proba(toy, use_pickle=True) - em.predict_proba(toy))) < 1e-14,
                      "model dir with an sklearn model.pkl: predict_proba == embedded numpy math (|diff| < 1e-14)")
            if args.upstream_dir:
                up = Path(args.upstream_dir)
                v = verify_upstream(up / "tests/expected_output/generic/K562_chr22/dhs_intact_hic", up)
                check(v["max_abs_score_diff"] < 1e-9 and v["thresholded_pairs_identical"] and v["stats_max_rel_diff"] < 1e-9
                      and all(r.get("coef_identical") for r in v.get("pickles", {}).values()),
                      f"upstream: pickles == embedded coefficients; chr22 expected scores reproduced (max diff {v['max_abs_score_diff']:.1e}), "
                      f"{v['n_thresholded']} thresholded pairs identical, stats identical")

            print("\nintervals and aggregation")
            A = pd.DataFrame({"chr": ["chr1"] * 3, "g": ["X", "X", "Y"], "s": [100, 500, 100], "e": [200, 600, 200]})
            B = pd.DataFrame({"chr": ["chr1"] * 4, "g": ["X", "X", "X", "Y"], "s": [200, 150, 601, 300], "e": [250, 160, 700, 400], "v": [1.0, 3.0, 5.0, 7.0]})
            ia, ib = keyed_overlap(A, B, ("chr", "g", "s", "e"), ("chr", "g", "s", "e"), closed=True)
            ia2, _ = keyed_overlap(A, B, ("chr", "g", "s", "e"), ("chr", "g", "s", "e"), closed=False)
            check(list(zip(ia, ib)) == [(0, 0), (0, 1)] and list(ia2) == [0],
                  "keyed_overlap: same gene required; closed (IRanges) intervals count touching ends, half-open do not")
            ag = aggregate_hits(np.array([0, 0, 1, 1]), [1.0, np.nan, 2.0, 4.0], "max")
            check(np.isnan(ag[0]) and ag[1] == 4.0, "aggregate: R semantics, any NA in the group -> NA")
            cr_t = pd.DataFrame({"chrom": ["chr1", "chr1"], "chromStart": [100, 900], "chromEnd": [200, 1000], "measuredGeneSymbol": ["X", "X"],
                                 "startTSS": [5000, 5000], "endTSS": [5001, 5001], "dataset": "d", "Regulated": [True, False]})
            ft_t = pd.DataFrame({"chr": ["chr1"], "start": [150], "end": [180], "TargetGene": ["X"], "distanceToTSS": [np.nan], "f": [np.nan]})
            tss_t = pd.DataFrame({"chr": ["chr1"], "start": [3000], "end": [3500], "name": ["X"]})
            mg, ms, _ = overlap_features_with_crispr(cr_t, ft_t, feature_table_df([("distanceToTSS", "distance", "NA", "min", "NA", ""),
                                                                                  ("f", "f", "NA", "max", "7", "")]), tss_t)
            check(len(mg) == 1 and len(ms) == 1 and mg["distanceToTSS"].iloc[0] == abs(150 - 3250) and mg["f"].iloc[0] == 7.0,
                  "crispr overlap: NA distanceToTSS refilled from the TSS reference; NA filled with fill_value; unmatched -> missing")

            print("\nthreshold / bedpe / stats")
            pr = pd.DataFrame({"chr": ["chr1"] * 5, "start": [0, 0, 100, 1000, 2000], "end": [10, 10, 300, 1100, 2100], "name": list("aabcd"),
                               "class": ["intergenic", "intergenic", "promoter", "promoter", "genic"], "TargetGene": ["g1", "g2", "g1", "g2", "g2"],
                               "TargetGeneTSS": [200, 1050, 200, 1050, 1050], "TargetGeneIsExpressed": [True, False, True, True, True],
                               "isSelfPromoter": [False, False, True, False, False], "distanceToTSS": [195.0, 1045.0, 0.0, 1000.0, 1000.0],
                               SCORE_COL: [0.5, 0.3, 0.9, 0.9, 0.1]})
            t1 = threshold_predictions(pr, 0.25)
            t2 = threshold_predictions(pr, 0.25, include_self_promoter=False)
            check(list(t1["name"]) == ["a", "a", "b"] and list(t2["name"]) == ["a", "a"],
                  "threshold: score >= t; promoters kept only when self-promoter (include_self_promoter)")
            bp = bedpe_frame(t1)
            check(list(bp["name"]) == ["g1_a", "g1_b"] and list(bp["y1"]) == [200, 200], "bedpe: expressed genes only, TargetGene_name, TSS as the second anchor")
            st = dict(zip(*get_stats(t1).T.values))
            check(st["num_enh"] == 2 and st["num_genes"] == 2 and st["num_enh_gene_links"] == 3 and abs(st["mean_num_genes_per_enh"] - 1.5) < 1e-12
                  and abs(st["mean_num_enh_per_gene_no_prom"] - 1.0) < 1e-12 and abs(st["mean_enh_region_size"] - 105.0) < 1e-12
                  and abs(st["mean_log10_dist_to_tss"] - (math.log10(195) + math.log10(1045) + 0) / 3) < 1e-12,
                  "get_stats: counts, means, log10(0) -> 0, enhancer size over unique elements")

            print("\nfeatures (synthetic ABC world)")
            W = synthetic_world(d)
            _parse(flags("features", "--abc-dir", str(W["abc"]), "--model", "dhs_intact_hic", "--feature-table", str(W["custom_ft"]),
                         "--tss", str(W["tss"]), "--chr-sizes", str(W["sizes"]), "--gene-classes", str(W["gc"]),
                         "--external-features-config", str(W["ext"]), "--biosample", "SYN", "--label", "st_re2g_features"))
            fr = _last_run("st_re2g_features")
            bdir = fr / "SYN"
            feats = read_tsv(bdir / "genomewide_features.tsv.gz")
            pred = W["pred"]
            check(len(feats) == len(pred) and (fr / "report.md").is_file() and (fr / "summary.json").is_file(),
                  f"features: genomewide table with one row per ABC pair ({len(feats)}), report + summary")
            nonp = pred[pred["class"] != "promoter"].copy()
            nonp["mid"] = ((nonp["start"] + nonp["end"]) / 2).astype(int)
            samp = nonp.sample(200, random_state=4)
            cand = read_tsv(bdir / "NumCandidateEnhGene.tsv").set_index(["name", "TargetGene"])["NumCandidateEnhGene"]
            tssc = read_tsv(bdir / "NumTSSEnhGene.tsv")
            tssc = tssc.set_index([tssc.columns[0], tssc.columns[1]])[tssc.columns[2]]
            tb = read_bed(W["tss"], ["chr", "start", "end", "name", "score", "strand"])
            ok_c = ok_t = True
            for r in samp.itertuples(index=False):
                same = nonp[nonp["TargetGene"] == r.TargetGene]["mid"].to_numpy()
                exp_c = int(((same >= r.mid) & (same < r.TargetGeneTSS)).sum()) if r.mid < r.TargetGeneTSS else int(((same <= r.mid) & (same > r.TargetGeneTSS)).sum())
                ok_c &= int(cand[(r.name, r.TargetGene)]) == exp_c
                a0, a1 = (r.TargetGeneTSS, r.end) if r.TargetGeneTSS < r.mid else (r.start, int(r.mid + r.distance))
                tc = tb[tb["chr"] == r.chr]
                exp_t = int(((tc["start"] < a1) & (tc["end"] > a0)).sum())
                ok_t &= int(tssc.get((r.name, r.TargetGene), 0)) == exp_t
            check(ok_c, "numCandidateEnhGene == brute-force rank from the TSS on each side (200 pairs)")
            check(ok_t, "numTSSEnhGene == brute-force count of TSS500bp intervals in [start, mid+distance) / [TSS, end) (200 pairs)")
            nn = read_tsv(bdir / "NumEnhancersEG5kb.txt", header=None, names=["name", "n"]).set_index("name")["n"]
            ns = read_tsv(bdir / "SumEnhancersEG5kb.txt", header=None, names=["name", "s"]).set_index("name")["s"]
            el = W["el"].copy()
            el["mid"] = ((el["start"] + el["end"]) / 2).astype(int)
            ok_n = True
            for r in el.sample(120, random_state=5).itertuples(index=False):
                ws, we = max(0, r.mid - 5000), min(W["chroms"][r.chr], r.mid + 5000)
                hit = pred[(pred["chr"] == r.chr) & (pred["start"] < we) & (pred["end"] > ws) & (pred["name"] != r.name)]
                ok_n &= int(nn.get(r.name, 0)) == len(hit) and abs(float(ns.get(r.name, 0.0)) - hit["activity_base"].sum()) < 1e-6
            check(ok_n, "num/sumNearbyEnhancers == brute-force count / activity sum of other prediction rows within midpoint +/- 5 kb (120 elements)")
            check(nn.max() > 1 and feats["numNearbyEnhancers"].isna().sum() == 0, "nearby counts include one row per linked gene (upstream sort -u no-op); NA filled 0")
            dd_n, _ = nearby_enhancers(pred, read_tsv(W["abc"] / "Neighborhoods" / "EnhancerList.txt"), W["chroms"], dedupe=True)
            check(dd_n["count"].sum() < nn.sum(), f"--dedupe-nearby counts each neighbouring element once ({dd_n['count'].sum()} < {nn.sum()})")
            act = read_tsv(bdir / "ActivityOnly_features.tsv.gz")
            al = feats.merge(pred[["name", "TargetGene", "hic_contact_pl_scaled_adj"]].rename(columns={"hic_contact_pl_scaled_adj": "hic"}),
                             on=["name", "TargetGene"], how="left")
            f_ok = np.allclose(al["ABCxContact"], al["ABC.Score"] * al["hic"]) and np.allclose(al["contactFrequency"], al["hic"]) \
                and "hic_contact_pl_scaled_adj" not in feats.columns and len(al) == len(feats)
            check(f_ok and "ABC.Score" in act.columns and "numCandidateEnhGene" in act.columns,
                  "final features: interaction term input_col*second_input, input_col renamed to feature")
            gn = pd.read_csv(d / "gene_noise.tsv", sep="\t")
            have = feats["TargetGene"].isin(set(gn["TargetGene"]))
            gm = feats.merge(gn, on="TargetGene", how="left")
            mean_fill = float(gm["gnoise"].mean())
            check(np.allclose(feats.loc[have, "geneNoise"], gm.loc[have, "gnoise"]) and np.allclose(feats.loc[~have, "geneNoise"], mean_fill),
                  f"external TargetGene join; missing genes filled with the 'mean' fill value ({mean_fill:.3f})")
            nz = W["noise"]
            key = feats[["chr", "start", "end", "TargetGene"]].astype(str).agg(":".join, axis=1)
            dup = nz[nz.duplicated(["chr", "end", "TargetGene"], keep=False)]
            dkey = dup["chr"].astype(str) + ":" + (dup["start"] - dup.groupby(["chr", "end", "TargetGene"])["start"].transform(lambda x: x - x.min())).astype(str) \
                + ":" + dup["end"].astype(str) + ":" + dup["TargetGene"]
            exp_mean = dup.groupby(dkey)["noise_val"].mean()
            got = feats.assign(k=key).set_index("k")["noiseFeature"]
            check(len(exp_mean) >= 15 and np.allclose(got.reindex(exp_mean.index).to_numpy(), exp_mean.to_numpy()),
                  f"external overlap join: {len(exp_mean)} pairs with two source intervals aggregated by mean")
            nomatch = ~key.isin(set(nz["chr"].astype(str) + ":" + nz["start"].astype(str) + ":" + nz["end"].astype(str) + ":" + nz["TargetGene"]))
            check(nomatch.sum() > 50 and (feats.loc[nomatch, "noiseFeature"] == 0).all(), f"external overlap join: {int(nomatch.sum())} unmatched pairs -> NA -> fill 0")
            check(np.allclose(act["ABC.Denominator"], pred["ABC.Score"] / pred["ABC.Score.Numerator"]) if "ABC.Denominator" in act else True,
                  "ABC.Denominator = ABC.Score / ABC.Numerator (upstream definition)")

            print("\napply / run / qc")
            _parse(flags("apply", "--features", str(bdir / "genomewide_features.tsv.gz"), "--model", "dhs_intact_hic", "--model", "dhs_megamap",
                         "--biosample", "SYN", "--label", "st_re2g_apply"))
            ar = _last_run("st_re2g_apply")
            sc = read_tsv(ar / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions.tsv.gz")
            exp_sc = em.predict_proba(feats.replace([np.inf, -np.inf], np.nan).fillna(0))
            thr_f = ar / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions_threshold0.243.tsv.gz"
            th = read_tsv(thr_f)
            check(np.allclose(sc[SCORE_COL], exp_sc, rtol=1e-12) and len(th) == int(((sc[SCORE_COL] >= 0.243) & ((sc["class"] != "promoter") | sc["isSelfPromoter"])).sum()),
                  f"apply: ENCODE-rE2G.Score written, {len(th)} pairs >= 0.243 (self-promoters kept)")
            bpf = ar / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions_threshold0.243.bedpe"
            nb = sum(1 for _ in open(bpf))
            check(nb == int(th["TargetGeneIsExpressed"].sum()) and (ar / "SYN" / "dhs_megamap" / "encode_e2g_predictions_threshold0.201_stats.tsv").is_file(),
                  f"apply: bedpe ({nb} expressed links) and stats for both models")
            bsc = pd.DataFrame({"biosample": ["SYN"], "DHS": ["syn.tagAlign.gz"], "ATAC": [None], "H3K27ac": [None], "default_accessibility_feature": ["DHS"],
                                "HiC_file": ["https://x/ENCFF621AIY.hic"], "HiC_type": ["hic"], "HiC_resolution": [5000], "ABC_directory": [str(W["abc"])],
                                "external_features_config": [None]})
            bsc.to_csv(d / "biosamples.tsv", sep="\t", index=False)
            _parse(flags("run", "--biosample-config", str(d / "biosamples.tsv"), "--tss", str(W["tss"]), "--chr-sizes", str(W["sizes"]),
                         "--gene-classes", str(W["gc"]), "--label", "st_re2g_run", np_))
            rr_ = _last_run("st_re2g_run")
            sr = json.loads((rr_ / "summary.json").read_text())
            rsc = read_tsv(rr_ / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions.tsv.gz")
            rj = rsc.merge(sc[["name", "TargetGene", SCORE_COL]], on=["name", "TargetGene"], suffixes=("", "_apply"))
            check(len(sr["results"]) == 1 and sr["results"][0]["model"] == "dhs_intact_hic" and len(rj) == len(sc)
                  and np.allclose(rj[SCORE_COL], rj[SCORE_COL + "_apply"]),
                  "run: model chosen from the biosample config (dhs_intact_hic); scores identical to features + apply")
            if not args.no_plots:
                check(len(sr["figures"]) >= 10, f"run: {len(sr['figures'])} QC figures")
            stats_files = [str(ar / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions_threshold0.243_stats.tsv"),
                           str(ar / "SYN" / "dhs_megamap" / "encode_e2g_predictions_threshold0.201_stats.tsv")]
            _parse(flags("qc-plots", "--stats", *stats_files, "--label", "st_re2g_qc", np_))
            qs = json.loads((_last_run("st_re2g_qc") / "summary.json").read_text())
            check(qs["n_datasets"] == 2 and (_last_run("st_re2g_qc") / "outlier_stats.tsv").is_file(), "qc-plots: stats table + top/bottom-5 outlier table")
            for tag, f in (("ENCSR000AAA", stats_files[0]), ("ENCSR000AAB", stats_files[1])):
                for side in ("a", "b"):
                    (d / side / tag).mkdir(parents=True, exist_ok=True)
                    shutil.copy(f, d / side / tag / "x_stats.tsv")
            _parse(flags("compare-stats", "--a-stats", *[str(p) for p in (d / "a").glob("*/x_stats.tsv")], "--b-stats",
                         *[str(p) for p in (d / "b").glob("*/x_stats.tsv")], "--label", "st_re2g_cmpstats", np_))
            check(json.loads((_last_run("st_re2g_cmpstats") / "summary.json").read_text())["n_matched"] == 2, "compare-stats: runs matched by ENCODE accession")
            _parse(["threshold", "--predictions", str(ar / "SYN" / "dhs_intact_hic" / "encode_e2g_predictions.tsv.gz"), "--threshold", "0.243",
                    "--label", "st_re2g_thr"])
            check(json.loads((_last_run("st_re2g_thr") / "summary.json").read_text())["n_out"] == len(th), "threshold subcommand == apply's thresholded table")

            print("\nCRISPR training data")
            C = plant_crispr(d, W, feats)
            _parse(["crispr-features", "--features", str(bdir / "genomewide_features.tsv.gz"), "--crispr", str(C["path"]), "--tss", str(W["tss"]),
                    "--model", "dhs_intact_hic", "--feature-table", str(W["custom_ft"]), "--dataset", "SYN", "--label", "st_re2g_crispr"])
            cr_run = _last_run("st_re2g_crispr")
            cs = json.loads((cr_run / "summary.json").read_text())
            check(cs["n_missing"] == C["n_gap"] and cs["n_overlapping"] == C["n_sampled"] + 2 and cs["n_crispr"] == C["n_sampled"] + 2 + C["n_gap"] + 4,
                  f"crispr-features: {cs['n_overlapping']} overlapping, {cs['n_missing']} missing (planted {C['n_gap']}), 4 genes outside the TSS universe removed")
            mg = read_tsv(Path(cs["for_training"]))
            e1, e2, gname = C["span"]
            sp = mg[mg["name"] == "span"].iloc[0]
            check(abs(sp["ABC.Score"] - (e1["ABC.Score"] + e2["ABC.Score"])) < 1e-12 and sp["distanceToTSS"] == min(e1["distanceToTSS"], e2["distanceToTSS"])
                  and sp["numTSSEnhGene"] == max(e1["numTSSEnhGene"], e2["numTSSEnhGene"]),
                  "CRISPR element spanning two candidates: ABC.Score summed, distanceToTSS min, numTSSEnhGene max (aggregate_function)")
            tch = mg[mg["name"] == "touch"].iloc[0]
            check(abs(tch["ABC.Score"] - C["touch"][0]["ABC.Score"]) < 1e-12, "CRISPR element touching a candidate's start overlaps it (closed intervals)")
            check(list(mg.columns[:4]) == ["dataset", "chr", "start", "end"] and "TargetGene" in mg.columns and "pair_uid" not in mg.columns,
                  "process_crispr_data: columns renamed to chr/start/end/TargetGene; uid columns dropped")

            print("\ntraining")
            ftr = str(cr_run / f"for_training.crispr.SYN_features_NAfilled.tsv.gz")
            _parse(flags("train", "--crispr-features", ftr, "--model", "dhs_intact_hic", "--dataset", "SYN", "--label", "st_re2g_train", np_))
            tr = _last_run("st_re2g_train")
            ts = json.loads((tr / "summary.json").read_text())
            met = read_tsv(tr / "performance_metrics.tsv")
            tp = read_tsv(tr / "training_predictions.tsv")
            yv = regulated_labels(tp["Regulated"])
            base = yv.mean()
            allr = met[met["test_chr"] == "all"].iloc[0]
            check(list(met.columns) == METRIC_COLS and set(met["test_chr"]) == {"chr1", "chr2", "chr3", "all"} and ts["missing"],
                  "train: per-chromosome + pooled metrics in upstream column order; missing-pairs file found")
            check(allr["AUPRC_test"] > 2.5 * base and allr["AUROC_test"] > 0.8,
                  f"train: LOCO CV AUPRC {allr['AUPRC_test']:.3f} (base rate {base:.3f}), AUROC {allr['AUROC_test']:.3f}")
            Xtr = design_matrix(tp, em.features)
            planted = _expit(C["planted_scale"] * (Xtr.to_numpy() @ em.coef + em.intercept) + C["planted_shift"])
            full_p = tp[MODEL_NAME + ".Score_full"].to_numpy()
            ll = lambda q: float(np.sum(yv * np.log(q) + (1 - yv) * np.log(1 - q)))  # noqa: E731
            check(ll(full_p) >= ll(planted) - 1e-3, f"train: full-model log-likelihood {ll(full_p):.2f} >= planted model {ll(planted):.2f} (MLE)")
            from scipy.stats import spearmanr  # type: ignore
            rho = spearmanr(tp[SCORE_COL], planted).correlation
            check(rho > 0.8, f"train: CV scores rank-correlate with the planted probabilities (Spearman {rho:.3f})")
            cf = np.array([ts["coefficients"][f] for f in em.features])
            wi, bi = _irls(Xtr.to_numpy(), yv.astype(float))
            check(np.max(np.abs(wi - cf)) < 0.05 and abs(bi - ts["intercept"]) < 0.2,
                  f"train: numpy Newton fit agrees with the {ts['backend']} fit (max |coef diff| {np.max(np.abs(wi - cf)):.1e})")
            md = load_model(str(tr / "model_dir"))
            check(np.allclose(md.predict_proba(tp), full_p, rtol=1e-10) and md.threshold == ts["model_threshold"] == round(ts["threshold_70_pct_recall_cv"], 3),
                  "train: model_dir (feature_table + threshold_ + model.json) reproduces the full model; usable by apply")
            _parse(["train", "--crispr-features", ftr, "--model", "dhs_intact_hic", "--label", "st_re2g_train_compat", "--upstream-compat", "--no-plots"])
            mc = read_tsv(_last_run("st_re2g_train_compat") / "performance_metrics.tsv")
            mca = mc[mc["test_chr"] == "all"].iloc[0]
            check(abs(mca["AUROC_test_full"] - mca["log_loss_test_full"]) < 1e-12 and abs(allr["AUROC_test_full"] - allr["log_loss_test_full"]) > 1e-3,
                  "--upstream-compat reproduces log_loss written into AUROC_test_full of the 'all' row")
            _parse(["train", "--crispr-features", ftr, "--model", "dhs_intact_hic", "--label", "st_re2g_train_poly", "--polynomial", "--no-plots",
                    "--override-params", "{'penalty': 'l2', 'C': 1.0}"])
            tpoly = json.loads((_last_run("st_re2g_train_poly") / "summary.json").read_text())
            check(len(tpoly["coefficients"]) == 1 + 8 + 36 and tpoly["params"]["penalty"] == "l2", "train --polynomial: 45 degree-2 terms; override params applied")
            pm = load_model(str(_last_run("st_re2g_train_poly") / "model_dir"))
            check(pm.polynomial and np.allclose(pm.predict_proba(tp), read_tsv(_last_run("st_re2g_train_poly") / "training_predictions.tsv")[MODEL_NAME + ".Score_full"], rtol=1e-10),
                  "polynomial model dir can be applied (model.json keeps the design)")

            print("\nfeature analysis")
            ft4 = d / "ft4.tsv"
            write_tsv(feature_table_df([r for r in EMBEDDED_MODELS["dhs_intact_hic"]["feature_table"] if r[0] in ("ABC.Score", "numTSSEnhGene")]
                                       + [r for r in CUSTOM_FT if r[0] in ("noiseFeature", "geneNoise")]), ft4, quiet=True)
            common = ["--crispr-features", ftr, "--feature-table", str(ft4), "--n-boot", "40", "--seed", "1"]
            _parse(flags("feature-selection", *common, "--direction", "forward", "--label", "st_re2g_ffs", np_))
            fs = json.loads((_last_run("st_re2g_ffs") / "summary.json").read_text())
            tab = pd.DataFrame(fs["table"])
            noise_d = tab[tab["feature_added"].isin(["noiseFeature", "geneNoise"])]["delta_aupr"].abs().max()
            check(fs["order"][0] == "ABC.Score" and len(tab) == 5 and tab.iloc[1]["delta_aupr"] > 0.1 and noise_d < 0.06,
                  f"forward selection: ABC.Score first (delta AUPRC {tab.iloc[1]['delta_aupr']:.3f}), noise features change it by < 0.06 ({noise_d:.3f})")
            _parse(flags("feature-selection", *common, "--direction", "backward", "--label", "st_re2g_bfs", np_))
            bsj = json.loads((_last_run("st_re2g_bfs") / "summary.json").read_text())
            check(bsj["order"][0] == "None" and bsj["order"][-1] == "ABC.Score" and len(bsj["table"]) == 5,
                  f"backward selection: removal order {bsj['order'][1:]} keeps ABC.Score to the end")
            _parse(flags("permutation-importance", "--crispr-features", ftr, "--feature-table", str(ft4), "--n-repeats", "3", "--seed", "1",
                         "--label", "st_re2g_pfi", np_))
            md_ = json.loads((_last_run("st_re2g_pfi") / "summary.json").read_text())["mean_delta"]
            check(md_["ABC.Score"]["delta_aupr"] < -0.1 and abs(md_["noiseFeature"]["delta_aupr"]) < 0.06 and abs(md_["geneNoise"]["delta_aupr"]) < 0.06,
                  f"permutation importance: ABC.Score {md_['ABC.Score']['delta_aupr']:.3f}, noise {md_['noiseFeature']['delta_aupr']:.3f} / {md_['geneNoise']['delta_aupr']:.3f}")
            _parse(["all-feature-sets", *common, "--label", "st_re2g_sets"])
            sets = read_tsv(_last_run("st_re2g_sets") / "all_feature_sets.tsv")
            check(len(sets) == 15 and sets["AUPRC"].is_monotonic_decreasing and sets.iloc[0]["ABC.Score"] == 1 and "ABC.Score" in sets.iloc[0]["features"],
                  "all feature sets: 2^4 - 1 subsets, sorted by AUPRC, best contains ABC.Score")
            ftn = d / "ftn.tsv"
            write_tsv(feature_table_df([r for r in CUSTOM_FT if r[0] in ("noiseFeature", "geneNoise")]), ftn, quiet=True)
            _parse(["train", "--crispr-features", ftr, "--feature-table", str(ftn), "--dataset", "SYN", "--label", "st_re2g_train_noise", "--no-plots"])
            _parse(flags("compare-models", "--train-dirs", str(tr), str(_last_run("st_re2g_train_noise")), "--crispr", str(C["path"]),
                         "--n-boot", "40", "--seed", "2", "--label", "st_re2g_cmp", np_))
            cm = read_tsv(_last_run("st_re2g_cmp") / "performance_across_models.tsv").set_index("model")
            check(len(cm) == 3 and cm.loc["st_re2g_train", "AUPRC"] > cm.loc["st_re2g_train_noise", "AUPRC"] + 0.1
                  and cm.loc["st_re2g_train", "pct_missing_elements"] > 0 and "distance" in cm.index,
                  f"compare-models: trained {cm.loc['st_re2g_train', 'AUPRC']:.3f} > noise-only {cm.loc['st_re2g_train_noise', 'AUPRC']:.3f}; "
                  f"distance baseline {cm.loc['distance', 'AUPRC']:.3f}; missing pairs scored 0")
            _parse(["models", "--label", "st_re2g_models"])
            mm = load_model(str(_last_run("st_re2g_models") / "models" / "extended"))
            check(len(mm.features) == 45 and mm.threshold == 0.336 and np.array_equal(mm.coef, embedded_model("extended").coef),
                  "models: exported model dirs round-trip (extended: 45 features, threshold 0.336)")
    finally:
        after = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
        for p in after - before:
            if "_st_re2g_" in p.name:
                shutil.rmtree(p, ignore_errors=True)
        if not existed and OUT_ROOT.exists() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="igvfagent encode-re2g", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def refs(p):
        p.add_argument("--tss", help="TSS500bp reference bed (default: setup resources)")
        p.add_argument("--chr-sizes", help="chromosome sizes (clips the 5 kb windows; default: setup resources)")
        p.add_argument("--gene-classes", help="gene_promoter_class table (is_ubiquitous_uniform, P2PromoterClass)")

    def training(p, model_required=False):
        p.add_argument("--crispr-features", required=True, help="for_training.*.tsv.gz from crispr-features")
        p.add_argument("--feature-table", help="feature_table.tsv (feature,input_col,second_input,aggregate_function,fill_value,nice_name)")
        p.add_argument("--model", help="take the feature table of this embedded model / model dir")
        p.add_argument("--polynomial", action="store_true", help="degree-2 PolynomialFeatures before the log transform")
        p.add_argument("--epsilon", type=float, default=EPSILON)
        p.add_argument("--override-params", default="", help="dict string overriding the LogisticRegression defaults, e.g. \"{'C': 1, 'penalty': 'l2'}\"")

    p = sub.add_parser("setup", help="fetch reference files from the pinned upstream commit")
    p.add_argument("--test-data", action="store_true", help="also fetch the chr22 expected output (~47 MB) for verify-upstream")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("models", help="list / export the embedded pretrained models")
    p.add_argument("--export", help="directory to write one model dir per model (default: the run dir)")
    p.add_argument("--label", default="models")
    p.set_defaults(func=cmd_models)

    p = sub.add_parser("select-model", help="biosample config -> model(s) and thresholds")
    p.add_argument("--biosample-config", required=True)
    p.add_argument("--model-root", help="directory of model folders (default: the embedded models)")
    p.add_argument("--megamap-hic", default=MEGAMAP_HIC_FILE)
    p.add_argument("--label", default="select_model")
    p.set_defaults(func=cmd_select_model)

    p = sub.add_parser("features", help="ABC outputs -> genome-wide ENCODE-rE2G feature table")
    p.add_argument("--abc-dir", help="ABC biosample dir with Predictions/EnhancerPredictionsAllPutative.tsv.gz and Neighborhoods/EnhancerList.txt")
    p.add_argument("--abc-predictions")
    p.add_argument("--enhancer-list")
    p.add_argument("--model", nargs="+", action="extend", help="embedded model name(s) or model dir(s) whose feature tables to build")
    p.add_argument("--feature-table", nargs="+", action="extend", help="extra feature_table.tsv file(s)")
    p.add_argument("--external-features-config", help="input_col,source_col,aggregate_function,join_by(overlap|TargetGene),source_file")
    p.add_argument("--dedupe-nearby", action="store_true", help="count each neighbouring element once (upstream's intended sort -u)")
    p.add_argument("--biosample")
    refs(p)
    p.add_argument("--label", default="features")
    p.set_defaults(func=cmd_features)

    p = sub.add_parser("apply", help="score genome-wide features with model(s), threshold, bedpe, stats")
    p.add_argument("--features", required=True, help="genomewide_features.tsv.gz")
    p.add_argument("--model", nargs="+", action="extend", required=True)
    p.add_argument("--threshold", type=float, help="override the model threshold")
    p.add_argument("--exclude-self-promoter", action="store_true", help="include_self_promoter = False")
    p.add_argument("--accessibility", nargs="+", action="extend", help="accessibility BAM(s) for num_sequencing_reads")
    p.add_argument("--epsilon", type=float, default=EPSILON)
    p.add_argument("--use-pickle", action="store_true", help="score with the unpickled sklearn model when the model dir has one")
    p.add_argument("--biosample")
    p.add_argument("--label", default="apply")
    p.set_defaults(func=cmd_apply)

    p = sub.add_parser("run", help="features + apply for each biosample of a config")
    p.add_argument("--biosample-config", required=True, help="ABC biosample config (+ optional model_dir, ABC_directory, external_features_config)")
    p.add_argument("--abc-results", help="root with one ABC output dir per biosample (when there is no ABC_directory column)")
    p.add_argument("--model-root")
    p.add_argument("--megamap-hic", default=MEGAMAP_HIC_FILE)
    p.add_argument("--exclude-self-promoter", action="store_true")
    p.add_argument("--dedupe-nearby", action="store_true")
    p.add_argument("--epsilon", type=float, default=EPSILON)
    refs(p)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("threshold", help="threshold a prediction table")
    p.add_argument("--predictions", required=True)
    p.add_argument("--threshold", type=float, required=True)
    p.add_argument("--score-column", default=SCORE_COL)
    p.add_argument("--exclude-self-promoter", action="store_true")
    p.add_argument("--label", default="threshold")
    p.set_defaults(func=cmd_threshold)

    p = sub.add_parser("bedpe", help="IGV bedpe of expressed-gene links")
    p.add_argument("--predictions", required=True)
    p.add_argument("--score-column", default=SCORE_COL)
    p.add_argument("--label", default="bedpe")
    p.set_defaults(func=cmd_bedpe)

    p = sub.add_parser("stats", help="get_stats for a thresholded prediction table")
    p.add_argument("--predictions", required=True)
    p.add_argument("--accessibility", nargs="+", action="extend")
    p.add_argument("--label", default="stats")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("qc-plots", help="QC figures over many *_stats.tsv files")
    p.add_argument("--stats", nargs="+", action="extend", required=True)
    p.add_argument("--encode-metadata", help="ENCODE metadata (DNase Experiment accession, Biosample type) to split cells / tissues")
    p.add_argument("--y2ave-metadata", help="cluster metadata (CellClusterID, nCells, MeanATACFragmentsPerCell, ManualAnnotationLabel)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="qc")
    p.set_defaults(func=cmd_qc_plots)

    p = sub.add_parser("compare-stats", help="one metric across two sets of runs matched by an accession")
    p.add_argument("--a-stats", nargs="+", action="extend", required=True)
    p.add_argument("--b-stats", nargs="+", action="extend", required=True)
    p.add_argument("--a-name", default="V1")
    p.add_argument("--b-name", default="V2")
    p.add_argument("--metric", default="num_enh_gene_links")
    p.add_argument("--key-regex", default=r"ENC[A-Z0-9]+")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="compare_stats")
    p.set_defaults(func=cmd_compare_stats)

    p = sub.add_parser("crispr-features", help="overlap genome-wide features with CRISPR E-G pairs; training table")
    p.add_argument("--features", required=True)
    p.add_argument("--crispr", required=True, help="EPCrisprBenchmark-format TSV")
    p.add_argument("--model", nargs="+", action="extend")
    p.add_argument("--feature-table", nargs="+", action="extend")
    p.add_argument("--tss")
    p.add_argument("--no-na-fill", action="store_true")
    p.add_argument("--dataset", help="dataset name used in output file names")
    p.add_argument("--label", default="crispr_features")
    p.set_defaults(func=cmd_crispr_features)

    p = sub.add_parser("train", help="train a model: full + leave-one-chromosome-out")
    training(p)
    p.add_argument("--missing", help="missing CRISPR pairs file (default: sibling missing.* of --crispr-features)")
    p.add_argument("--dataset")
    p.add_argument("--upstream-compat", action="store_true", help="reproduce log_loss in AUROC_test_full of the 'all' row")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="train")
    p.set_defaults(func=cmd_train)

    p = sub.add_parser("feature-selection", help="forward / backward sequential feature selection with bootstrap")
    training(p)
    p.add_argument("--direction", choices=["forward", "backward"], default="forward")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="feature_selection")
    p.set_defaults(func=cmd_feature_selection)

    p = sub.add_parser("permutation-importance", help="permutation feature importance")
    training(p)
    p.add_argument("--n-repeats", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="permutation_importance")
    p.set_defaults(func=cmd_permutation_importance)

    p = sub.add_parser("all-feature-sets", help="AUPRC of every feature subset")
    training(p)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--force", action="store_true", help="allow >= 14 features")
    p.add_argument("--label", default="all_feature_sets")
    p.set_defaults(func=cmd_all_feature_sets)

    p = sub.add_parser("compare-models", help="compare trained models (+ distance baseline)")
    p.add_argument("--train-dirs", nargs="+", action="extend", required=True, help="run dirs of `train`")
    p.add_argument("--crispr", help="raw CRISPR table for the distance baseline")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="compare_models")
    p.set_defaults(func=cmd_compare_models)

    p = sub.add_parser("verify-upstream", help="re-score the upstream chr22 expected output with the embedded model")
    p.add_argument("--upstream-dir", help="a clone of EngreitzLab/ENCODE_rE2G (also compares the model pickles)")
    p.add_argument("--expected-dir", help="dir with the three expected files (default: setup --test-data location)")
    p.add_argument("--label", default="verify_upstream")
    p.set_defaults(func=cmd_verify_upstream)

    p = sub.add_parser("selftest", help="synthetic ABC world, planted model, every subcommand")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--upstream-dir", help="also verify against an upstream clone")
    p.set_defaults(func=cmd_selftest)
    return ap


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
