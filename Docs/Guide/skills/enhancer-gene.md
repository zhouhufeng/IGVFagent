# Enhancer–gene linkage (E2G)

Enhancer–gene predictions and their benchmarks: ABC, ENCODE-rE2G, scE2G, eQTL and GWAS enrichment, super-enhancer targets, and E2G QC and Portal submission.

**On this page**

- [Enhancer–gene linkage](#enhancergene-linkage)
- [Super-enhancer → target-gene pipeline](#super-enhancer--target-gene-pipeline)
- [scE2G workbench: train, check and benchmark E2G models (`sce2g`)](#sce2g-workbench-train-check-and-benchmark-e2g-models-sce2g)
- [eQTL enrichment benchmark for enhancer-gene predictors (`eqtl-enrich`)](#eqtl-enrichment-benchmark-for-enhancer-gene-predictors-eqtl-enrich)
- [GWAS E2G benchmark (`gwas-e2g`)](#gwas-e2g-benchmark-gwas-e2g)
- [ENCODE-rE2G (`encode-re2g`)](#encode-re2g-encode-re2g)
- [ABC pipeline (`abc-pipeline`)](#abc-pipeline-abc-pipeline)
- [scE2G pipeline (`sce2g-pipeline`)](#sce2g-pipeline-sce2g-pipeline)
- [E2G QC, predictions and IGVF Portal submission (`e2g-qc-predictions`)](#e2g-qc-predictions-and-igvf-portal-submission-e2g-qc-predictions)

## Enhancer–gene linkage

```bash
python3 Scripts/enhancer_gene_linkage_skills.py overview --source catalog --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source encode --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py pull-sets \
  --region chr1:903900-904900 --gene SAMD11 --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets \
  --include-local-catalog --demo-if-empty
python3 Scripts/enhancer_gene_linkage_skills.py write-playbook
```

## Super-enhancer → target-gene pipeline

End-to-end discovery of H3K27ac/BRD4/MED1/P300 ChIP-seq experiments,
ROSE-style super-enhancer calling, ranked-signal "hockey-stick" plots,
multi-track browser views of top SEs, and SE → target-gene linkage via
four parallel streams: Hi-C / ChIA-PET loops, IGVF Catalog rE2G
predictions, TSS proximity windows, and constituent cCRE composition.
Optional motif enrichment (CTCF, AP-1, GATA1, ETS, NFkB, STAT1, FOXA1,
TP53, MYC, SP1) and SE↔cCRE-density correlation. Playbook:
[`Docs/Skills/SE_TARGET_PIPELINE_SKILLS.md`](../../Skills/SE_TARGET_PIPELINE_SKILLS.md).

```bash
# 1) One-shot discovery → SE call → linkage → plots → optional motifs
igvfagent se-targets pipeline --biosample GM12878 --assembly GRCh38 \
  --label gm12878_h3k27ac --max-experiments 1 --top-n 10 \
  --genome /path/to/GRCh38.fa  # motif enrichment is optional

# 2) Just the discovery half (lists candidate SE-driver experiments)
igvfagent se-targets discover --biosample K562 --assembly GRCh38 \
  --target H3K27ac

# 3) Re-run linkage on an existing SE table (skip discovery & calling)
igvfagent se-targets link-targets \
  --se-bed Docs/SETargets/<run>/super_enhancers.bed \
  --label gm12878_relinked --max-ses 100
```

The pipeline is bounded: per-Catalog-call timeouts default to 20 s and
linkage stops after 3 consecutive failures so a slow Catalog endpoint
can never hang the run. Outputs land in `Docs/SETargets/<timestamp>_<label>/`
with `super_enhancers.bed`, ranked SE → target-gene CSVs, the
hockey-stick/Top-SE browser PNGs, and a markdown report.

## scE2G workbench: train, check and benchmark E2G models (`sce2g`)

[scE2G](https://github.com/EngreitzLab/scE2G) (MIT) is a Snakemake + R workflow
with a Singularity container; its scientific asset is the trained model, so this
skill does not reimplement it. It wraps the Engreitz group's "Quick Start on
Building New E2G Models" walkthrough for adding crowdsourced features to the
multiome model, and reimplements the evaluation of
[CRISPR_comparison](https://github.com/EngreitzLab/CRISPR_comparison) (MIT)
clean-room so a new model can be scored against CRISPR element-gene pairs here.

```bash
igvfagent sce2g selftest                                  # fake checkout, synthetic inputs
igvfagent sce2g setup      --repo-dir ~/scE2G             # clone fix/dag-staleness-integration + patches
igvfagent sce2g features   --repo-dir ~/scE2G --name h3k27ac --features K562_features.tsv
igvfagent sce2g configure  --repo-dir ~/scE2G --cluster K562 --name h3k27ac \
    --rna K562_rna_count_matrix.csv.gz --atac-frag K562_atac_fragments.tsv.gz --profile profiles/slurm
igvfagent sce2g check      --repo-dir ~/scE2G --model multiome_h3k27ac
igvfagent sce2g run        --repo-dir ~/scE2G --model multiome_h3k27ac --profile profiles/slurm --execute
igvfagent sce2g predictions --predictions new=results/K562/multiome_h3k27ac/*.e2g.tsv --predictions base=...
igvfagent sce2g benchmark  --predictions new=... --predictions base=... \
    --crispr EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz --pred-config pred_config.txt
```

`setup` applies the walkthrough's four patches idempotently (drop
`conda: "mamba"` from both `Snakefile_training` files, add `SCRIPTS_DIR`, add
`RNA_matrix_filtered` / `max_cell_count`, restore the missing
`multiome_arc_n6.tsv`). `features` renames `ElementChr/Start/End/GeneSymbol` to
`chr/start/end/TargetGene`, replaces spaces in feature names, and writes the
five-column `external_features_config_<name>.tsv` plus a feature table with one
`max / 0 / NA / nice_name` row per feature. `configure` writes the cluster and
model rows in upstream's exact column order and a run script with the Slurm
profile and memory/runtime overrides. `check` catches the mistakes that
otherwise surface hours into training: a feature in the external config that
the feature table never lists, a source file missing a column, a dataset with
no cluster row. `benchmark` scores each CRISPR pair with the pred_config's
aggregate / fill / inverse semantics, reports AUPRC with bootstrap intervals
and precision at 70% recall, and writes `pred_config.txt` + `config.yml` for
the upstream pipeline. With `--all-features` every numeric column of every
table becomes its own predictor, which is how the crowdsourced K562 feature
tables on Synapse (syn73717888, 16 tables, 11 million element-gene rows each)
were scored: tables are streamed and restricted to CRISPR-tested genes, so no
table is loaded whole. `inventory` describes such a folder first, and `merge`
combines batch runs into one ranked table, figure and report. The scorer was
checked against the upstream pipeline itself, run in a container on the same
inputs: per-pair scores identical for all 10,356 K562 CRISPR pairs, and the
AUPRC on the upstream's own definition (`auprc_crispr_comparison`) equal to
four decimals for every predictor tried (benchmark #24). Agent tools:
`sce2g_setup`, `sce2g_features`, `sce2g_configure`, `sce2g_check`,
`sce2g_run`, `sce2g_predictions`, `sce2g_benchmark`, `sce2g_inventory`,
`sce2g_benchmark_merge`, `sce2g_selftest`.

`explain_dataset --download` on IGVF files now works: the portal answers with
a 307 to a pre-signed S3 URL, and forwarding the portal credentials to S3 made
every download fail with HTTP 400; credentials are dropped on cross-host
redirects.

## eQTL enrichment benchmark for enhancer-gene predictors (`eqtl-enrich`)

Port of [EngreitzLab/eQTLEnrichment](https://github.com/EngreitzLab/eQTLEnrichment),
the eQTL benchmark of the ENCODE-rE2G and scE2G papers. It asks whether a
predictor's enhancers are enriched for fine-mapped eQTL variants against 1000G
SNPs, and whether they link those variants to the right eGene. Every rule is
re-derived in pandas/numpy/scipy with no bedtools dependency. Variants and
background SNPs are restricted to distal noncoding regions. Enrichment has
log risk-ratio CIs and hypergeometric p-values. Recall is reported both as
total and as linking to the eGene, across a quantile threshold span and by
eVariant-eGene distance bin. The run also writes aggregated all-matches
curves, enrichment at target recalls with pairwise tests, enhancer set sizes
and heatmaps.

```bash
igvfagent eqtl-enrich setup --synapse            # resources + background SNPs + GTEx SuSiE release
igvfagent eqtl-enrich prepare-gtex --raw GTEx_30tissues_release1.tsv.gz \
    --expression GTEx_median_tpm.gct.gz --out gtex.eqtl.tsv.gz
igvfagent eqtl-enrich run --methods-table methods.tsv --predictions-table predictions.tsv \
    --eqtl gtex.eqtl.tsv.gz --bg-variants all.bg.SNPs.hg38.baseline.v1.1.bed.sorted --label k562_blood
igvfagent eqtl-enrich selftest
```

Two upstream quirks are fixed by default and reproducible with
`--upstream-compat`. The first is the misplaced square root in the SE of the
log enrichment. The second is the unthresholded background counts in the
by-distance tables. Agent tools: `eqtl_enrichment_setup`,
`eqtl_enrichment_prepare_gtex`, `eqtl_enrichment_variants`,
`eqtl_enrichment_run`, `eqtl_enrichment_selftest`.

## GWAS E2G benchmark (`gwas-e2g`)

A Python port of [EngreitzLab/GWAS_E2G_benchmarking](https://github.com/EngreitzLab/GWAS_E2G_benchmarking), the GWAS benchmark from the scE2G and ENCODE-rE2G papers. It tests enhancer-gene predictions against fine-mapped UK Biobank GWAS variants in two ways. **Variant overlap**: for each trait x biosample, what fraction of the distal-noncoding variants with PIP > 0.1 fall in predicted enhancers (recall), and how enriched they are relative to 1000G common SNPs. This is reported at the method's threshold and across a quantile threshold span (enrichment-recall curves). **Gene linking**: how precisely the top-2 predicted genes of each credible set recover the silver-standard causal gene, alone and intersected with PoPS, with PoPS and distance-to-TSS as baselines.

The port keeps the upstream config format (config.yml, methods / predictions / comparisons tables, biosample and trait groups, ALL), output tables and column names, and runs without bedtools or R. Five upstream quirks are corrected by default and reproduced with `--upstream-compat`: the SE formula, quantiles taken over distinct scores, unthresholded single-biosample linking, the gene-universe filter being dropped, and ALL double counting. `setup` fetches the small resources from the pinned commit. The background SNPs come from Synapse.

```bash
igvfagent gwas-e2g setup --traits RBC MCV HbA1c Lym --synapse
igvfagent gwas-e2g validate-config --config config/config.yml
igvfagent gwas-e2g run --config config/config.yml --label sce2g_blood
igvfagent gwas-e2g baseline --label ukbb_baselines
igvfagent gwas-e2g plot --run-dir Docs/GWASE2G/<run> --plot-fixed-scale
igvfagent gwas-e2g selftest --no-plots
```

Agent tools: `gwas_e2g_setup`, `gwas_e2g_validate_config`, `gwas_e2g_variants`, `gwas_e2g_baseline`, `gwas_e2g_run`, `gwas_e2g_plot`.

## ENCODE-rE2G (`encode-re2g`)

A Python port of [EngreitzLab/ENCODE_rE2G](https://github.com/EngreitzLab/ENCODE_rE2G), the logistic-regression enhancer-to-gene model of the ENCODE encyclopedia (Gschwind et al. 2026). Both upstream workflows are covered: the apply workflow (model choice from the ABC biosample config, the new features computed from ABC outputs, external features, final features, scoring, thresholding, IGV bedpe, per-prediction stats and QC plots) and the training workflow (CRISPR overlap, leave-one-chromosome-out training, forward/backward feature selection with bootstrap, permutation importance, all feature subsets, model comparison with a distance baseline). bedtools, csvtk and GenomicRanges are replaced by vectorised interval arithmetic, so only pandas and numpy are needed (scikit-learn and scipy are used when present).

The nine pretrained models are embedded as coefficients taken from the upstream pickles, so scoring does not depend on pickle compatibility. `verify-upstream` re-scores the upstream's own K562 chr22 test output: the scores match to 7e-16, the 7,218 thresholded links are identical, and so is the stats table.

```bash
igvfagent encode-re2g setup                       # TSS universe, chrom sizes, gene classes, CRISPR benchmark
igvfagent encode-re2g select-model --biosample-config config_biosamples.tsv
igvfagent encode-re2g features --abc-dir results/K562 --model dhs_intact_hic --biosample K562
igvfagent encode-re2g apply --features Docs/ENCODErE2G/<run>/K562/genomewide_features.tsv.gz --model dhs_intact_hic
igvfagent encode-re2g crispr-features --features genomewide_features.tsv.gz --crispr Data/ENCODErE2G/resources/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz --model dhs_intact_hic
igvfagent encode-re2g train --crispr-features for_training.*.tsv.gz --model dhs_intact_hic --label my_model
igvfagent encode-re2g feature-selection --crispr-features for_training.*.tsv.gz --model dhs_intact_hic --direction forward
igvfagent encode-re2g selftest --no-plots
```

Agent tools: `encode_re2g_models`, `encode_re2g_select_model`, `encode_re2g_features`, `encode_re2g_apply`, `encode_re2g_run`, `encode_re2g_stats`, `encode_re2g_qc_plots`, `encode_re2g_crispr_features`, `encode_re2g_train`, `encode_re2g_feature_selection`, `encode_re2g_permutation_importance`, `encode_re2g_all_feature_sets`, `encode_re2g_compare_models`, `encode_re2g_verify_upstream`.

## ABC pipeline (`abc-pipeline`)

A full port of [broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) (MIT), pinned at `92ac5036`. It runs the whole Snakemake workflow in Python: MACS2 peaks (the binary is optional; there is also a MACS-like Python fallback), candidate regions (the 150,000 strongest peaks, summit +/-250 bp, blocklist removed, TSS include-list added), neighborhoods (DHS/ATAC/H3K27ac reads from BAM, tagAlign, fragments, bigWig or bedGraph; quantile normalisation to the K562 reference; activity = sqrt(DHS x H3K27ac)), predictions (power law or Hi-C from `.hic`, juicebox, bedpe or average-Hi-C directories), the automatic `abc_thresholds.tsv` threshold, filtering, and QC. It also covers the Hi-C utilities: power-law fit, average Hi-C, juicebox dump, and per-gene Hi-C bedgraphs. Column names and file names match upstream.

It was checked against upstream's own expected test outputs for chr22. On the ATAC tagAlign sample (power law) and the DNase + H3K27ac BAM sample, candidate regions, Counts.bed, GeneList.txt, EnhancerList.txt values, and all 1.9-2.4 M ABC / power-law scores are identical, and so are the thresholded file and GenePredictionStats. The small clean-room model-only scorer `igvfagent abc score` is still available alongside it.

```bash
igvfagent abc-pipeline setup
igvfagent abc-pipeline run --biosample K562 --dhs dnase.bam --h3k27ac h3k27ac.bam --hic-file ENCFF621AIY.hic --hic-type hic --hic-resolution 5000 --label k562
igvfagent abc-pipeline run --biosamples-table config_biosamples.tsv --label batch
igvfagent abc-pipeline predict --enhancers EnhancerList.txt --genes GeneList.txt --chrom-sizes sizes.tsv --accessibility-feature ATAC
igvfagent abc-pipeline powerlaw-fit --hic-dir juicebox_dir --hic-type juicebox
igvfagent abc-pipeline selftest --no-plots
```

Agent tools: `abc_pipeline_setup`, `abc_pipeline_call_peaks`, `abc_pipeline_fragments_to_tagalign`, `abc_pipeline_candidate_regions`, `abc_pipeline_neighborhoods`, `abc_pipeline_predict`, `abc_pipeline_filter`, `abc_pipeline_variant_overlap`, `abc_pipeline_qc`, `abc_pipeline_powerlaw_fit`, `abc_pipeline_average_hic`, `abc_pipeline_split_avg_hic`, `abc_pipeline_juicebox_dump`, `abc_pipeline_hic_bedgraph`, `abc_pipeline_compare`, `abc_pipeline_run`, `abc_pipeline_selftest`.

## scE2G pipeline (`sce2g-pipeline`)

A Python rewrite of [EngreitzLab/scE2G](https://github.com/EngreitzLab/scE2G) (MIT, pinned at 7cb2af7) that runs without Snakemake, R, Signac, bedtools or fast_kendall_sc. It covers every rule: fragments to tagAlign and counts, Kendall peak-gene pairs, the peak x cell ATAC matrix, the Kendall tau-b between accessibility and expression across cells (exact closed form), RNA pseudobulk TPM and detection, ARC-E2G, the ENCODE_rE2G feature tables scE2G assembles, model application with quantile normalisation and the TPM filter, thresholding, gene and element lists, QC statistics against the Sheth, Qiu 2024 reference clusters, the CRISPR benchmark with bootstrap CIs, and training new models. The four published v3 models ship inside the module (feature tables, thresholds and the logistic-regression weights read from model.pkl and the 23 held-out-chromosome models); `setup` fetches their qnorm references.

It matches upstream outputs exactly. On the upstream chr22 test fixture, the port reproduces all 27,685 Kendall pairs (max difference 5e-16). On the IGVF K562 scE2G v1.2 release, it rebuilds the features from the released files and reproduces all 11.5 M `Score` and `Score.ignoreTPM` values (max difference 7e-16). ABC peak calling and predictions come from the ABC pipeline and are an input here (`--abc-dir`).

```bash
igvfagent sce2g-pipeline setup --gtf --crispr --example
igvfagent sce2g-pipeline run --cluster K562 --fragments atac_fragments.tsv.gz --rna-matrix rna.h5ad \
    --abc-dir ABC/K562 --models multiome_powerlaw_v3 scATAC_powerlaw_v3 --crispr Data/scE2G/pipeline_resources/EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz
igvfagent sce2g-pipeline predict --features genomewide_features.tsv.gz --models multiome_powerlaw_v3 --cluster K562 --bedpe
igvfagent sce2g-pipeline validate-example
igvfagent sce2g-pipeline selftest --no-plots
```

Agent tools: `sce2g_pipeline_run`, `sce2g_pipeline_kendall`, `sce2g_pipeline_arc`, `sce2g_pipeline_predict`, `sce2g_pipeline_benchmark`, `sce2g_pipeline_train`, `sce2g_pipeline_qc`, plus setup / models / frag-to-tagalign / kendall-pairs / activity-features / features / crispr-features / validate-example / validate-release.

## E2G QC, predictions and IGVF Portal submission (`e2g-qc-predictions`)

Port of [kaybrand/QC-and-Predictions](https://github.com/kaybrand/QC-and-Predictions) (MIT), the IGVF E2G Pillar pipeline that takes QC-filtered pseudobulk clusters to shareable scE2G enhancer-gene products. Every step is rewritten in Python: the parse-time quality gate (100 cells, 2e6 fragments, 1e6 UMIs; merged, prefiltered and CATlas clusters), ATAC / RNA filtering with full-barcode matching and GTF symbol collapse, scE2G configs, portal-format reformatting (metadata header, consortium column order, element BED and bedpe with bgzip + tabix), candidate and feature tables, dataset-wide QC aggregation and plots, the IGVF Portal Cell Annotation cache and cell annotation report, file discovery and archive comparison, Synapse manifests, and the CATlas distance-vs-depth analysis. scE2G itself runs through `sce2g-pipeline`, whose QC statistics and figures are reused.

The IGVF submission builds all eleven metadata tables (fifteen files per cluster: principal pseudobulk set, prediction set, filtered barcode list, ATAC fragments and index, RNA matrix, five prediction tabular files, ATAC bigWig, two index files, QC document) with the upstream aliases, links and controlled vocabulary, orders them in dependency rounds, validates required fields, resolves external references against the Portal graph, and writes iu_register TSVs, REST payloads and `upload_plan.json`. It is a dry run by default; `--execute` POSTs / PATCHes in round order to the IGVF sandbox (`--production` for the real Portal) with `IGVF_ACCESS_KEY` / `IGVF_SECRET_ACCESS_KEY` from the environment. The `synapse-submission` and `CATlas-predictions` branch variants are available as flags.

```bash
igvfagent e2g-qc-predictions resolve-exclusions --config igvf0_pipeline_config.yaml
igvfagent e2g-qc-predictions cell-metadata --config igvf0_pipeline_config.yaml
igvfagent e2g-qc-predictions run --config igvf0_pipeline_config.yaml --label igvf0
igvfagent e2g-qc-predictions manifest --config igvf0_pipeline_config.yaml --cluster-keys igvf0 --lineage
igvfagent e2g-qc-predictions manifest --config igvf0_pipeline_config.yaml --cluster-keys igvf0 --execute   # sandbox
igvfagent e2g-qc-predictions selftest --no-plots
```

Agent tools: `e2g_qc_resolve_exclusions`, `e2g_qc_merge_metrics`, `e2g_qc_prefiltered_metrics`, `e2g_qc_build_qc_datatables`, `e2g_qc_filter_atac`, `e2g_qc_filter_rna`, `e2g_qc_package_rna`, `e2g_qc_sce2g_config`, `e2g_qc_reformat`, `e2g_qc_candidates`, `e2g_qc_features`, `e2g_qc_aggregate_qc`, `e2g_qc_stale_reformats`, `e2g_qc_cell_metadata`, `e2g_qc_cell_annotation_report`, `e2g_qc_dataset_accessions`, `e2g_qc_washu_report`, `e2g_qc_verify_fragments`, `e2g_qc_portal_files`, `e2g_qc_compare_archive`, `e2g_qc_manifest`, `e2g_qc_patch_submitter_comment`, `e2g_qc_report`, `e2g_qc_synapse_manifest`, `e2g_qc_synapse_orphans`, `e2g_qc_distance_depth`, `e2g_qc_run`, `e2g_qc_selftest`.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
