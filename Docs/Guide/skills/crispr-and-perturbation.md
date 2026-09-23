# CRISPR screens and Perturb-seq

Pooled screens from raw reads, CRISPRi, Flow-FISH, Perturb-seq pipelines (IGVF-CRISPR and pinellolab), SCEPTRE, TF Perturb-seq and the CRISPR jamboree analyses.

**On this page**

- [CRISPR screens from raw reads — three shapes, three tools](#crispr-screens-from-raw-reads--three-shapes-three-tools)
- [CRISPRi / CRISPR-FACS / Perturb-seq](#crispri--crispr-facs--perturb-seq)
- [CRISPRi Flow-FISH screen](#crispri-flow-fish-screen)
- [TF Perturb-seq: calibrated effects to disease and GWAS (`tf-perturb`)](#tf-perturb-seq-calibrated-effects-to-disease-and-gwas-tf-perturb)
- [IGVF CRISPR Pipeline, early Nextflow version (`igvf-crispr-pipeline`)](#igvf-crispr-pipeline-early-nextflow-version-igvf-crispr-pipeline)
- [SCEPTRE for IGVF MuData (`sceptre-igvf`)](#sceptre-for-igvf-mudata-sceptre-igvf)
- [CRISPRi-FlowFISH pipeline (`flowfish-pipeline`)](#crispri-flowfish-pipeline-flowfish-pipeline)
- [Pooled CRISPR screen toolkit (`perturb-tools`)](#pooled-crispr-screen-toolkit-perturb-tools)
- [Perturb-seq pipeline, bulk_crispr_pipeline port (`bulk-crispr`)](#perturb-seq-pipeline-bulk_crispr_pipeline-port-bulk-crispr)
- [CRISPR module dev tools (`crisprdevtools`)](#crispr-module-dev-tools-crisprdevtools)
- [CRISPR-FG scPerturb-seq jamboree (`crispr-fg-jamboree`)](#crispr-fg-scperturb-seq-jamboree-crispr-fg-jamboree)
- [Fishash Table 2 reproduction (`fishash-table2`)](#fishash-table-2-reproduction-fishash-table2)
- [CRISPR SeqSpec (`crispr-seqspec`)](#crispr-seqspec-crispr-seqspec)
- [Gasperini 2019 pipeline (`gasperini-pipeline`)](#gasperini-2019-pipeline-gasperini-pipeline)
- [IGVF CRISPR Perturb-seq pipeline (`crispr-pipeline`)](#igvf-crispr-perturb-seq-pipeline-crispr-pipeline)
- [2nd CRISPR Jamboree analyses (`crispr-jamboree2`)](#2nd-crispr-jamboree-analyses-crispr-jamboree2)
- [3rd CRISPR Jamboree single-cell pipeline (`crispr-jamboree3`)](#3rd-crispr-jamboree-single-cell-pipeline-crispr-jamboree3)
- [4th CRISPR Jamboree inference benchmark (`crispr-jamboree4`)](#4th-crispr-jamboree-inference-benchmark-crispr-jamboree4)
- [Single-cell CRISPR differential expression (`sc-crispr-de`)](#single-cell-crispr-differential-expression-sc-crispr-de)

## CRISPR screens from raw reads — three shapes, three tools

A sorted CRISPR screen is not analysable one bin at a time: the measurement
**is** the comparison between bins, so each tool finds the screen's siblings
from the Portal and analyses them together. Which tool applies depends on the
screen's shape and on its guide **library** — and the wrong tool refuses
rather than returning an undercount.

```bash
# 1. TAIL SORT — bottom20% vs top20%
igvfagent crispr-screen discover IGVFDS5542IBUS      # every bin and replicate
igvfagent crispr-screen analyze  IGVFDS5542IBUS --tail 20 --label ldlr

# 2. LETTERED-BIN GRADIENT — BinA..BinF, no tails (554 Portal MeasurementSets)
igvfagent gradient-screen discover IGVFDS8710ZSOZ
igvfagent gradient-screen analyze  IGVFDS8710ZSOZ --label kitlg

# 3. BASE EDITING — ABE/CBE, following crispr-bean's method
igvfagent bean discover IGVFDS6464SOVZ    # editor + which BEAN inputs IGVF publishes
igvfagent bean count    IGVFDS6464SOVZ    # measure which masking recovers reads
igvfagent bean analyze  IGVFDS6464SOVZ --label ldl_abe
```

**Why three tools rather than one.** A tail sort asks whether a construct is
enriched in the low or the high tail. A six-bin gradient asks what a
construct's mean expression is, computed from all bins at once. And a
base-editing screen cannot be counted by exact guide matching at all,
because the editor edits the guide's own locus: on `IGVFDS6464SOVZ` exact
matching assigns **36.7%** of reads and `bean`'s masked matching assigns
**62.5%**. `crispr-screen` detects a base-editing library and refuses,
naming `bean`, because the base editor is stated only in the library's guide
names — the readout (`gRNA sequencing`) and assay title (`CRISPR FACS
screen`) are identical to an ordinary knockout screen.

**The counting key is chosen by measurement, never assumed.** A
prime-editing library shares one spacer across every variant it installs, so
keying on `spacer` collapses 1,741 pegRNAs onto 52 and credits each one's
reads to an arbitrary sibling. Each tool scores every candidate sequence
column against real reads and reports the table:

| column | separates | found in reads | ambiguous |
|---|---|---|---|
| `rt_template_sequence` | 99.8% | **62.4%** | 0.0% ← chosen |
| `spacer` | 3.0% | 7.3% | 1.7% |
| `peg_sequence` | 99.8% | **0.0%** | 0.0% |

`peg_sequence` separates that library perfectly and appears in nothing — its
entries are longer than the read. Counts are cached per library, so
re-analysing at a different `--tail` returns in seconds instead of minutes.

**`bean` also reports editing activity.** Per-guide self-editing rate is the
quantity BEAN's activity normalisation rests on. Effects come out as
`log2_raw` and `log2_per_edit`, but the p-value stays on `log2_raw`: dividing
by a rate estimated from finite counts is noisiest exactly where it changes
the answer most, so low activity is surfaced as **low power** rather than
divided out. `bean` states in its own output what is *not* available and why, derived from
the library's published columns and the screen's bins rather than from a fixed
list — a hard-coded list of absences becomes a confident description of an
older version of the tool, which is exactly what happened here: it went on
claiming variant-level aggregation and posterior intervals were missing after
`--run-bean` had started producing both.

**`bean analyze --run-bean` hands the counts to the real BEAN.** It writes
BEAN's four input tables, runs `bean create-screen`, then `bean run sorting
variant`, and reports BEAN's per-target posteriors alongside this tool's
per-guide scores. Two constraints are worth knowing before you rely on it,
both measured:

| Constraint | Consequence |
|---|---|
| `bean run sorting` needs an **unsorted bin** as `--control-condition` | Found by shared `construct_library_sets`, not by alias prefix. The depositor does not use one naming scheme per screen: `IGVFDS6464SOVZ`'s tails are `18loci_uptake_Rep*_bottom20`, its unsorted bins are `B*_18loci_Rep*_bulk`. Grouping by the parsed series name split one screen into two and dropped all four — which is why this tool spent several rounds reporting that the deposit had published only tails and BEAN could never model it. It had. `bean discover` reports the control condition up front, before any counting. A screen that genuinely publishes no unsorted bin is still refused rather than given a tail as its control, which would make BEAN succeed while measuring every guide against a baseline that is itself selected. |
| BEAN's activity-normalised model needs `X_bcmatch` | `--scale-by-acc` and `--guide-activity-col` both route through MixtureNormal, which reads `screen.layers['X_bcmatch']` — counts assigned by guide barcode, unpublished by IGVF → `KeyError: 'X_bcmatch'`. The models that avoid it (`--uniform-edit`, `--const-pi`) reject `--guide-activity-col` and write `edit_eff = 1.0` for every guide. So BEAN contributes the posterior; the activity-aware estimate stays this tool's `log2_per_edit`, and the output says which model produced which number. |

## CRISPRi / CRISPR-FACS / Perturb-seq

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
python3 Scripts/crispri_data_skills.py write-playbook
```

## CRISPRi Flow-FISH screen

Clean-room rewrite of
[EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline)
(MIT, Engreitz Lab 2021), following the published methods in
**Fulco 2019** *Nat Genet* and **Nasser 2021** *Nature*. Per-guide
log-normal MLE on bin counts with EM treatment of an "outside" overflow
bin, real-space conversion, per-element Mann-Whitney + Welch t-test +
BH-FDR. See
[`Docs/Skills/FLOWFISH_CRISPR_SKILLS.md`](../../Skills/FLOWFISH_CRISPR_SKILLS.md).

```bash
# Discover IGVF Portal CRISPRi-FlowFISH MeasurementSets
igvfagent flowfish pull-portal --limit 50 --label survey

# Generate a synthetic screen for smoke testing
igvfagent flowfish simulate --out-dir /tmp/ff_smoke \
    --n-elements 60 --guides-per-element 8 --knockdown-frac 0.30 \
    --cells-per-guide 800 --seed 11

# Per-guide log-normal MLE on bin counts
igvfagent flowfish estimate-effects \
    --counts /tmp/ff_smoke/counts.tsv \
    --sortparams /tmp/ff_smoke/sortparams.tsv --label demo

# Real-space + null normalization
igvfagent flowfish real-space \
    --input <ts>_demo_raw_effects.tsv \
    --target-col target --negative-label negative_control \
    --clamp 5 --label demo

# Per-element collapse + significance (MWU + Welch + BH-FDR)
igvfagent flowfish score-elements \
    --effects <ts>_demo_real_space.tsv \
    --target-col target --element-col ElementName \
    --negative-label negative_control \
    --min-guides 5 --fdr 0.05 --min-effect 0.10 --label demo
```

## TF Perturb-seq: calibrated effects to disease and GWAS (`tf-perturb`)

Port of the IGVF TF Perturb-seq consortium's Working Group 3 (disease & GWAS)
jamboree notebook,
[`perturb_seq_analysis_v4.ipynb`](https://github.com/IGVF/tf_perturb_seq/blob/main/docs/jamborees/2026_UTSW/working_groups/wg3_disease_gwas/perturb_seq_analysis_v4.ipynb),
as run on the Huangfu-lab HUES8 definitive-endoderm CRISPRi screen (2,161 TF
targets, 269,491 cells). It starts from the IGVF CRISPR pipeline's *calibrated*
inference tables, not raw counts:

```bash
igvfagent tf-perturb selftest                       # synthetic inputs, planted signals
igvfagent tf-perturb run \
    --calibrated-prefix data/<run>_calibrated_ \
    --mudata data/inference_mudata.h5mu \
    --gwas data/gwas-catalog-associations.tsv \
    --e2g ESC=data/h7_scE2G.e2g.tsv --e2g DE=data/definitive_endoderm_scE2G.e2g.tsv \
    --bigwig GSE213394_DED2-SOX17_hg38.bigwig --label huangfu_de
```

Stages, matching the notebook: guide-library overview; direct-target, cis and
trans significance (adjusted p < 0.05, |log2FC| > 0.2 cis / 1.0 trans,
targeting elements only, the 22M-row trans table streamed and cached); top
trans regulators with their strongest up/down targets; GWAS Catalog SNPs
within 50 kb of TF elements and of trans target genes with Fisher trait
enrichment; scE2G links for perturbed TFs and targets, and GWAS SNPs inside
those E2G elements (flagging assignments that disagree with the catalog's
nearest gene); optional ChIP-seq bigWig support. Output is
`Docs/TFPerturbSeq/<run>/` with `report.md`, every table as CSV, and PNG
figures. Positional overlaps are vectorised, `mudata` and `pyBigWig` are
optional, and the enrichment counts catalog associations exactly as upstream
does, with that caveat stated in the report. Agent tools:
`tf_perturb_seq_analyze`, `tf_perturb_seq_selftest`.

The rest of the [tf_perturb_seq](https://github.com/IGVF/tf_perturb_seq) package is
ported as stage subcommands of the same CLI, so a screen can be taken from the
IGVF CRISPR pipeline's `inference_mudata.h5mu` to calibrated tables without the
consortium's Synapse/Slurm wrappers:

```bash
igvfagent tf-perturb qc-gene   --mudata inference_mudata.h5mu       # Stage-3 gene-mapping QC
igvfagent tf-perturb qc-guide  --mudata inference_mudata.h5mu       # guide capture / guides per cell
igvfagent tf-perturb qc-target --mudata inference_mudata.h5mu       # knockdown, AUROC/AUPRC vs non-targeting
igvfagent tf-perturb calibrate --trans-results perturbo_trans_per_element_output.tsv \
    --mudata inference_mudata.h5mu --prefix run1 --method t-fit     # empirical p, BH, cis/trans tables
igvfagent tf-perturb pathways  --calibrated run1_calibrated_trans_results.tsv --gmt h.all.gmt --prefix run1
igvfagent tf-perturb filter-cells --mudata in.h5mu --out out.h5mu --max-guides 15
igvfagent tf-perturb edist-prep --mudata in.h5mu --out-dir edist/   # energy distance: PCA + guide dict
igvfagent tf-perturb edist-filter --prep-dir edist/                 # DISCO / hypergeometric / K-means outliers
igvfagent tf-perturb edist --prep-dir edist/                        # per-target E-distance vs NT backgrounds
igvfagent tf-perturb filter-guides --mudata in.h5mu --out out.h5mu \
    --targeting-outliers edist/targeting_outlier_table.csv --non-targeting-outliers edist/non_targeting_outlier_table.csv
igvfagent tf-perturb cnmf-export --mudata in.h5mu --out-h5ad perturbnmf.h5ad
igvfagent tf-perturb selftest --no-plots                            # WG3 + every stage on a planted h5mu
```

The energy-distance stages are a CPU rewrite of
[Chikara-Takeuchi/energy_dist_pipeline](https://github.com/Chikara-Takeuchi/energy_dist_pipeline)
that writes the same `pval_edist_full.csv` schema; `edist-validate`,
`edist-summary`, `cnmf-validate` and `pipeline-summary` check and roll up
outputs from either implementation. The self-test builds a 600-cell h5mu with a
planted TF knockdown and a shifted non-targeting guide, and every stage
recovers them. Agent tools: `tf_perturb_qc_gene`, `tf_perturb_qc_guide`,
`tf_perturb_qc_target`, `tf_perturb_calibrate`, `tf_perturb_pathways`,
`tf_perturb_filter_cells`, `tf_perturb_filter_guides`, `tf_perturb_edist_prep`,
`tf_perturb_edist_filter`, `tf_perturb_edist`, `tf_perturb_edist_summary`,
`tf_perturb_edist_validate`, `tf_perturb_cnmf_export`,
`tf_perturb_cnmf_validate`, `tf_perturb_pipeline_summary`.

## IGVF CRISPR Pipeline, early Nextflow version (`igvf-crispr-pipeline`)

A Python port of [IGVF-CRISPR/IGVF_CRISPR_Pipeline](https://github.com/IGVF-CRISPR/IGVF_CRISPR_Pipeline) (commit 6704a0f, 2024-07-18), the first Nextflow version of the IGVF CRISPR Perturb-seq pipeline. It reads a seqspec YAML directly and derives the kallisto-bustools technology string (`seqspec index -t kb`), the reads of each modality and the barcode onlist. It turns guide libraries into kite `guide_features.txt` and builds the `kb ref` / `kb count` commands, running them when kb is installed and printing them otherwise. It then runs the RNA AnnData QC (min_genes 100, min_cells 3, MT-/Mt- and RPS/RPL flags, scanpy QC metrics, knee/violin/scatter plots) and writes the MuData with `transcripts` and `guides` modalities (`ID|sequence` guide names, `number_of_nonzero_guides`, intersected barcodes).

It also adds a pure-Python kite-style guide counter, so the guide arm runs end to end without kallisto. The counter does exact or unique 1-mismatch protospacer matching, onlist barcode correction and distinct-UMI counts. On the upstream example it assigns 1,356 of 1,500 guide reads (1,213 exact, 164 rescued by one mismatch) to 1,147 barcodes. `--upstream-compat` brings back the upstream's positional guide naming, and its genome-as-`-f1` in kite mode.

```bash
igvfagent igvf-crispr-pipeline parse-seqspec --yaml guide.yml --modality guide --directory seqspec_dir
igvfagent igvf-crispr-pipeline guide-features --guide-table gasperini_tss.xlsx
igvfagent igvf-crispr-pipeline count-features --seqspec guide.yml --guide-table guides.xlsx --fastqs guide_R1.fq guide_R2.fq
igvfagent igvf-crispr-pipeline preprocess --adata-rna adata_rna.h5ad --gene-names cells_x_genes.genes.names.txt --reference human
igvfagent igvf-crispr-pipeline create-mdata --adata-rna filtered_anndata.h5ad --adata-guide adata_guide.h5ad --guide-metadata guides.xlsx
igvfagent igvf-crispr-pipeline run --rna-seqspec rna.yml --guide-seqspec guide.yml --guide-table guides.xlsx --guide-fastqs R1.fq R2.fq --adata-rna adata_rna.h5ad --gene-names genes.txt
igvfagent igvf-crispr-pipeline selftest --no-plots
```

Agent tools: `igvf_crispr_parse_seqspec`, `igvf_crispr_guide_features`, `igvf_crispr_kb_ref`, `igvf_crispr_kb_count`, `igvf_crispr_count_features`, `igvf_crispr_preprocess`, `igvf_crispr_create_mdata`, `igvf_crispr_run`.

## SCEPTRE for IGVF MuData (`sceptre-igvf`)

A Python port of [IGVF-CRISPR/sceptreIGVF](https://github.com/IGVF-CRISPR/sceptreIGVF), the R wrapper that runs SCEPTRE on the IGVF CRISPR MuData schema. It converts a MuData (mod/gene counts, mod/guide counts or `guide_assignment` layer, top-level covariates, `pairs_to_test`) into a SCEPTRE analysis, assigns gRNAs (Poisson-mixture EM, thresholding or maximum), runs SCEPTRE's QC, calibration check on non-targeting pairs, power check on positive controls and discovery analysis, and writes results back into MuData `uns`. The column names are the upstream ones: `test_results` for sceptreIGVF, or `per_element_results` / `per_guide_results` (with `cis_` / `trans_` prefixes) for the IGVF CRISPR_Pipeline.

The SCEPTRE statistics are re-implemented in numpy and statsmodels: an NB GLM score test, a conditional-resampling or permutation null reused across genes, and a skew-normal tail fit. This is a documented approximation of sceptre's C++ engine. `r-runner` runs the real R package when Rscript, MuData and sceptreIGVF are installed. `export-mudata` rebuilds the eight benchmark inputs/outputs MuData files.

```bash
igvfagent sceptre-igvf assign-guides --mudata guide_assignment_input.h5mu --label gasperini
igvfagent sceptre-igvf inference --mudata inference_input.h5mu --side left
igvfagent sceptre-igvf inference --mudata inference_input.h5mu --mode pipeline --scope cis
igvfagent sceptre-igvf run --mudata screen.h5mu --label screen
igvfagent sceptre-igvf r-runner --step inference --mudata inference_input.h5mu
igvfagent sceptre-igvf selftest --no-plots
```

Agent tools: `sceptre_igvf_convert`, `sceptre_igvf_assign_guides`, `sceptre_igvf_qc`, `sceptre_igvf_inference`, `sceptre_igvf_calibration_check`, `sceptre_igvf_power_check`, `sceptre_igvf_run`, `sceptre_igvf_export_mudata`, `sceptre_igvf_r_runner`, `sceptre_igvf_compare_results`.

## CRISPRi-FlowFISH pipeline (`flowfish-pipeline`)

A full port of the Engreitz-lab [crispri-flowfish](https://github.com/EngreitzLab/crispri-flowfish) Snakemake workflow, taking a screen from FASTQs to scored enhancers. The steps are: sample-sheet validation and experiment keys; read mapping (bowtie `-v0 --all`, or an exact-match Python mapper when bowtie is missing); guide counts and count tables per PCR and experimental replicate; replicate correlations; and the four sort-parameter loaders. From there it runs `estimate_effect_sizes.R` (weighted average plus a binned log-normal MLE with an EM seventh bin, using R `optim`/`stats4::mle` numerics), real-space normalisation to negative controls, 10-guide windows, TSS qPCR scaling, and collapse to candidate elements. Elements are scored with a Mann-Whitney test and a Student t-test, BH-corrected. The pipeline also does power analysis and writes ScreenData / KnownEnhancers tables. File names and columns match the upstream `results/` tree.

Two upstream quirks are corrected by default and can be reproduced with `--upstream-compat`. First, under R 3.6, `Count[Bin]` indexes by factor level code rather than by bin name. Second, bins and gates are paired by position. `CalculateTilingStatistic.R` depends on a private helper library, so its window test columns are reconstructed. The same applies to the KnownEnhancers formatter.

```bash
igvfagent flowfish-pipeline run --sample-sheet SampleSheet.tsv --design design.txt --sortparams-dir sortParams \
    --fastq-dir fastq --experiment-keycols CellLine --replicate-keycols FlowFISHRep \
    --qpcr qPCR.txt --genelist GeneList.txt --enhancers EnhancerList.bed --label ppif
igvfagent flowfish-pipeline estimate-effects --counts K562-Rep1.bin_counts.txt --sort-params B1_S1.txt
igvfagent flowfish-pipeline score-enhancers --collapsed K562-Rep1.collapse.bed --scaled K562-Rep1.scaled.txt --expt-name K562-Rep1
igvfagent flowfish-pipeline selftest --no-plots
```

Agent tools: `flowfish_pipeline_run`, `flowfish_pipeline_estimate_effects`, `flowfish_pipeline_real_space`, `flowfish_pipeline_windows`, `flowfish_pipeline_tss_kd`, `flowfish_pipeline_normalize_qpcr`, `flowfish_pipeline_collapse`, `flowfish_pipeline_score_enhancers`, `flowfish_pipeline_power`, `flowfish_pipeline_format_screen`, `flowfish_pipeline_count_tables`, `flowfish_pipeline_map_reads`, `flowfish_pipeline_samplesheet`, `flowfish_pipeline_guide_count_plots`.

## Pooled CRISPR screen toolkit (`perturb-tools`)

A port of [IGVF-CRISPR/perturb-tools](https://github.com/IGVF-CRISPR/perturb-tools) (MIT), the IGVF CRISPR working group's AnnData-based package for bulk pooled screens. It covers the whole package: the screen object (count table + guide/sample annotation, TKO sample-name parsing, adding technical replicates), log2(RPM + 1) normalisation, sample-vs-sample and per-replicate log fold changes with mean/median/sd aggregation, the sorting-screen delta-LFC with t-test p-values and the guide-enrichment-along-locus plot, the sample quality report (count distribution, Gini, count and LFC correlation, replicate consistency, outlier jackpot guides), MAGeCK / Excel / CSV export, the PoolQ output reader, guide annotation (protospacer, target, genomic position) and, from the dev branch, sgRNA library design (PAM scan of exons, BsmbI filter, GC / homopolymer flags) and feature pair distances.

Numbers match the upstream tutorials on the TKO HeLa screen to six decimals (e.g. A1BG guide 1: A/B/C.18_8.lfc = -0.341628 / -1.623461 / -0.351306), and log_norm reproduces PoolQ 3.3.2's lognormalized-counts to 4e-15. Upstream bugs (transposed outlier-guide indexing, missed BsmbI motif at position 0, untransposed MAGeCK layer export, flattened CSV matrix) are fixed, with `--upstream-compat` where the upstream behaviour runs.

```bash
igvfagent perturb-tools run --counts readcount-HeLa-lib1 --parse-tko-names --exclude replicate=0 \
    --cond1 18 --cond2 8 --compare-col time --target-col GENE --pos-ctrl core-essential-genes-sym_HGNCID
igvfagent perturb-tools sort-lfc --counts sort_counts.tsv --condit-1 high --condit-2 low --control presort --targets enh1 enh2
igvfagent perturb-tools read-poolq --poolq-dir poolq_out/ --sample-metadata conditions.csv
igvfagent perturb-tools design-library --fasta hg38.fa --gtf gencode.gtf --gene KMT2C --chrom chr7
igvfagent perturb-tools selftest --no-plots
```

Agent tools: `perturb_tools_run`, `perturb_tools_make_screen`, `perturb_tools_lfc_reps`, `perturb_tools_sort_lfc`, `perturb_tools_qc`, `perturb_tools_to_mageck`, `perturb_tools_read_poolq`, `perturb_tools_annotate_guides`, `perturb_tools_design_library`.

## Perturb-seq pipeline, bulk_crispr_pipeline port (`bulk-crispr`)

A port of [IGVF-CRISPR/bulk_crispr_pipeline](https://github.com/IGVF-CRISPR/bulk_crispr_pipeline). Despite the repository name, it is a single-cell Perturb-seq pipeline: Lucas Silva Ferreira's pipeline_perturbseq_like, the predecessor of IGVF_CRISPR_Pipeline, run on the Gasperini 2019 pilot. Every Nextflow process and jamboree task has a subcommand. Guides are counted from FASTQ the way kallisto kite does (exact and Hamming-1 matches, whitelist correction, UMI collapsing). Each lane then goes through cell QC: the knee, the minimum-genes rule, the mitochondrial fraction and Scrublet doublets. After that come guide binarisation (UMI > 5) with the upstream covariates and deMULTIplex MULTI-seq calls. Each element gets its genes within ±1 Mb, and every guide × gene pair gets a SCEPTRE-style conditional-resampling test. Results are aggregated with BH and Fisher, and BED and pyGenomeTracks links are written.

The binaries (kb, cellranger, Rscript/sceptre, pyGenomeTracks) are optional. Each one runs when it is on PATH; otherwise the exact command is printed. `--upstream-compat` brings back seven upstream bugs that are fixed by default: quality lines leaking into read composition, an ignored gene-filter fraction, `log_total_gene_count = log(n_genes+1)`, a transcript_end typo, random genes assigned to enhancers, `nUMI_total` classified as a hash, and a "Double" filter that keeps doublets.

```bash
igvfagent bulk-crispr guide-table --guides df_from_gasperini_tss.xlsx
igvfagent bulk-crispr count-guides --guides guide_features.txt --r1 guide_R1.fastq.gz --r2 guide_R2.fastq.gz --chemistry 10XV2 --whitelist 737K-august-2016.txt --name S1_L1
igvfagent bulk-crispr run --gtf Homo_sapiens.GRCh38.106.gtf --rna-dirs S1_L1_ks_transcripts_out --guide-dirs S1_L1_ks_guide_out --expected-cell-number 8000 --label gasperini_pilot
igvfagent bulk-crispr config --config perturb.config
igvfagent bulk-crispr selftest --no-plots
```

Agent tools: `bulk_crispr_run`, `bulk_crispr_count_guides`, `bulk_crispr_composition`, `bulk_crispr_guide_table`, `bulk_crispr_prefilter`, `bulk_crispr_multiseq`, `bulk_crispr_perturb_loader`, `bulk_crispr_de`, `bulk_crispr_results`, `bulk_crispr_tracks`, `bulk_crispr_assign_guides`, `bulk_crispr_config`, `bulk_crispr_map_rna`, `bulk_crispr_assay_spec`, `bulk_crispr_cellranger_inputs`, `bulk_crispr_inspect`.

## CRISPR module dev tools (`crisprdevtools`)

A port of [IGVF-CRISPR/crisprdevtools](https://github.com/IGVF-CRISPR/crisprdevtools), the helper used to start new modules of IGVF_CRISPR_Pipeline. Upstream has one function, `create_new_module_nextflow(name)`. `new-module` recreates its layout and file contents exactly: `bin/`, `conda_envs/`, `example_data/`, `processes/`, `test/`, `README.md`, `input.config`, a DSL2 `main.nf`, `bin/<name>.py`, `conda_envs/<name>.yaml` and `processes/<name>.nf`. Three upstream bugs are fixed by default and come back with `--upstream-compat`: the conda YAML is indented so it does not parse, every process is named `seqSpecParser`, and the bin script is not executable.

`--template full` also writes a working include plus workflow, an `input.config` params block, `test/test.nf`, example data and an argparse script skeleton. `check-module` validates any module directory: layout, DSL2 header, processes and includes, conda env files, bin shebangs and executable bits, and the params `main.nf` uses.

```bash
igvfagent crisprdevtools new-module --name guide_assignment --dest Modules
igvfagent crisprdevtools new-module --name guide_assignment --template full --dest Modules --force
igvfagent crisprdevtools check-module --path Modules/seqSpecParser
igvfagent crisprdevtools selftest
```

Agent tools: `crisprdevtools_new_module`, `crisprdevtools_check_module`.

## CRISPR-FG scPerturb-seq jamboree (`crispr-fg-jamboree`)

A Python port of the April-2023 IGVF CRISPR-FG jamboree ([IGVF-CRISPR/CRISPR_FG_JAMBOREE](https://github.com/IGVF-CRISPR/CRISPR_FG_JAMBOREE)) together with the pipeline steps its task notebooks plug into ([pipeline_perturbseq_like](https://github.com/LucasSilvaFerreira/pipeline_perturbseq_like)). Every task is a subcommand: assay name to kallisto read format and whitelist (Task 1), Cell Ranger `feature_ref.csv` / `library.csv` from a guide table, the pipeline `perturb.config` and launch command, cell QC into a guides/scRNA MuData, MuData QC with controls, coverage and guide-gene distances (Task 2), guide calling into a `binarized` layer, cis differential perturbation into `mudata_results.h5mu` (Task 3) and BED / links / pyGenomeTracks tracks (Task 4).

The differential module defaults to a SCEPTRE-style NB score test with conditional resampling, and has Mann-Whitney and Welch alternatives. `--engine r` writes the upstream `run_sceptre_high_moi` inputs. External binaries (cellranger, nextflow, Rscript, pyGenomeTracks) are optional: without them the skill prints the exact command. `--upstream-compat` reproduces four upstream quirks (covariate definitions, the overwritten Task 2 filter, random genes for unknown elements).

```bash
igvfagent crispr-fg-jamboree assay-spec --assay 10xv3
igvfagent crispr-fg-jamboree preprocess --rna lane1/rna --guides lane1/guides --gene-table Homo_sapiens.GRCh38.106.gtf.gz --expected-cells 5000
igvfagent crispr-fg-jamboree assign-guides --mudata raw_mudata_guide_and_transcripts.h5mu --method umi
igvfagent crispr-fg-jamboree differential --mudata mu_with_binary.h5mu --distance 1000000
igvfagent crispr-fg-jamboree tracks --results mudata_results.h5mu
igvfagent crispr-fg-jamboree selftest --no-plots
```

Agent tools: `crispr_fg_assay_spec`, `crispr_fg_cellranger_inputs`, `crispr_fg_pipeline_config`, `crispr_fg_preprocess`, `crispr_fg_mudata_qc`, `crispr_fg_guide_gene_distance`, `crispr_fg_subset`, `crispr_fg_assign_guides`, `crispr_fg_differential`, `crispr_fg_tracks`, `crispr_fg_run`.

## Fishash Table 2 reproduction (`fishash-table2`)

A Python port of [IGVF-CRISPR/fishash-table2-reproduction](https://github.com/IGVF-CRISPR/fishash-table2-reproduction). That repository reproduces the SCEPTRE and CLEANSER rows of the Fishash preprint's Table 2, a species-assignment benchmark on four GSE272457 human/mouse barnyard samples (CROP-seq and direct capture, 0 h and 72 h). The skill builds the per-sample inputs and scores them with the authors' cohort and accuracy/stderr definitions against the published values. It also brings its own guide callers. CLEANSER's zero-truncated cs/dc mixtures are re-implemented from the Stan models, with the same priors, bounds, normalisation and median-posterior output, sampled by adaptive Metropolis. SCEPTRE's mixture assignment is approximated with a Poisson GLM plus a reduced EM.

The real `cleanser` binary and the sceptre 0.10.3 R package are optional engines (`--engine cleanser`, `--engine r`). Without them the skill prints the exact upstream command. `compare` checks a scored table against Table 2 and against the upstream repository's own reproduction.

```bash
igvfagent fishash-table2 build-inputs --raw-dir GSE272457 --work-dir work
igvfagent fishash-table2 cleanser --input work/mix0hr_Cropseq_grna_counts.mtx --mode cs --out work/mix0hr_Cropseq_cleanser_cs_posterior.mtx
igvfagent fishash-table2 sceptre-mixture --work-dir work --sample mix0hr_Cropseq --raw-dir GSE272457
igvfagent fishash-table2 score --work-dir work --label table2
igvfagent fishash-table2 run --raw-dir GSE272457 --work-dir work
igvfagent fishash-table2 selftest --no-plots
```

Agent tools: `fishash_build_inputs`, `fishash_cleanser`, `fishash_sceptre_mixture`, `fishash_score`, `fishash_compare`, `fishash_run`.

## CRISPR SeqSpec (`crispr-seqspec`)

Port of [IGVF-CRISPR/CRISPR-SeqSpec](https://github.com/IGVF-CRISPR/CRISPR-SeqSpec), the CRISPR focus group's seqspec descriptions of CROP-seq + 10x v3 + MULTI-seq (one spec per modality) and TAP-seq (Schraivogel 2020). The skill carries a compact index of the four specs at the pinned commit and re-implements, without seqspec or PyYAML, what `seqspec check`, `seqspec index -t kb|chromap`, `seqspec onlist` and `seqspec print` do, for both the 0.2.x and legacy 0.0.x schemas.

It also checks submitted FASTQs against a spec (read lengths, fixed sequences at their coordinates, onlist hit rates). On the upstream files this finds that the guide FASTQs are 150 bp while guide.yml declares 26/43 bp reads, that the 26-bp R1 ends inside the 12-bp UMI, and that 88-95% of cell barcodes are on the whitelist.

```bash
igvfagent crispr-seqspec catalog
igvfagent crispr-seqspec check --spec guide
igvfagent crispr-seqspec index --spec tapseq --modality crispr --tool kb
igvfagent crispr-seqspec fetch && igvfagent crispr-seqspec check-reads --spec guide
```

Agent tools: crispr_seqspec_catalog, crispr_seqspec_check, crispr_seqspec_index, crispr_seqspec_onlist, crispr_seqspec_info, crispr_seqspec_check_reads, crispr_seqspec_fetch.

## Gasperini 2019 pipeline (`gasperini-pipeline`)

Port of [IGVF-CRISPR/Pipeline_Gasperini_2019](https://github.com/IGVF-CRISPR/Pipeline_Gasperini_2019), the IGVF processing of the Gasperini et al. 2019 CRISPRi pilot with the single-cell-like Nextflow pipeline (LucasSilvaFerreira/pipeline_perturbseq_like). Every stage is re-implemented in Python with the config's thresholds: guide table, read composition, kb commands, per-lane cell QC with doublets, MuData creation, guide assignment at UMI > 3, covariates, 1-Mb cis pairs, a SCEPTRE-style conditional randomisation test (or R SCEPTRE when available), BH/Fisher aggregation into mudata_results.h5mu, MULTI-seq demultiplexing, the original Gasperini NB likelihood-ratio test, and a comparison with GEO GSE120861. Six upstream quirks are corrected by default and reproduced with `--upstream-compat`.

On the sample pilot MuData the run takes about 10 s. Each of the 7 self-TSS positive controls that are in the expression matrix comes out a hit, and every hit at BH < 0.1 is also a hit in the original Gasperini results.

```bash
igvfagent gasperini-pipeline setup
igvfagent gasperini-pipeline inspect
igvfagent gasperini-pipeline run --gasperini-test --compare Data/GasperiniPipeline/GSE120861_all_deg_results.pilot.txt.gz
igvfagent gasperini-pipeline selftest
```

Agent tools: gasperini_pipeline_setup, gasperini_pipeline_run, gasperini_pipeline_qc_filter, gasperini_pipeline_mudata, gasperini_pipeline_guide_table, gasperini_pipeline_composition, gasperini_pipeline_multiseq, gasperini_pipeline_compare, gasperini_pipeline_inspect.

## IGVF CRISPR Perturb-seq pipeline (`crispr-pipeline`)

A Python port of [pinellolab/CRISPR_Pipeline](https://github.com/pinellolab/CRISPR_Pipeline), the IGVF consortium's single-cell CRISPR screen pipeline. Every step of the Nextflow workflow is a subcommand: seqspec parsing to kb technology strings, the seqSpecCheck read scan, guide / hashing feature references, kb count wrappers (with a pure-Python guide/hashtag counter), knee barcode filtering and QC, `inference_mudata.h5mu` creation with the consortium guide.var schema (non-targeting buckets, intended_target_key), SCEPTRE-mixture / CLEANSER / threshold guide assignment, cis pairs within 1 Mb, cis + trans inference, hashing demultiplexing, doublet removal, additional QC, evaluate_controls, IGV tracks, the ENCODE TF benchmark and the HTML dashboard.

External tools (kb/kallisto, seqspec, GMM-demux, cleanser, SCEPTRE in R, PerTurbo) run when installed; otherwise the exact upstream command is printed and a Python implementation of the method runs. Outputs keep the upstream column names and uns keys, and `uns['inference_engine']` records which engine produced them.

```bash
igvfagent crispr-pipeline seqspec --yaml guide_seqspec.yml --modalities guide --whitelist 737K-august-2016.txt
igvfagent crispr-pipeline map --modality guide --fastqs g_R1.fastq.gz g_R2.fastq.gz --features guide_metadata.tsv --technology 0,0,16:0,16,28:1,0,0 --barcodes 737K.txt --batch B1 --engine python
igvfagent crispr-pipeline run --rna B1_ks_transcripts_out B2_ks_transcripts_out --guide B1_ks_guide_out B2_ks_guide_out --guide-metadata guide_metadata.tsv --gtf gencode.v46.gtf.gz --label screen
igvfagent crispr-pipeline qc --mudata inference_mudata.h5mu
igvfagent crispr-pipeline tf-benchmark --mudata inference_mudata.h5mu --gtf gencode.v46.gtf.gz --encode-bed-dir Data/CRISPRPipeline/encode_bed_files
igvfagent crispr-pipeline selftest --no-plots
```

Agent tools: `crispr_pipeline_run`, `crispr_pipeline_seqspec`, `crispr_pipeline_seqspec_check`, `crispr_pipeline_feature_ref`, `crispr_pipeline_map`, `crispr_pipeline_concat`, `crispr_pipeline_preprocess`, `crispr_pipeline_create_mudata`, `crispr_pipeline_hashing`, `crispr_pipeline_doublets`, `crispr_pipeline_assign_guides`, `crispr_pipeline_pairs`, `crispr_pipeline_inference`, `crispr_pipeline_merge_results`, `crispr_pipeline_qc`, `crispr_pipeline_evaluate`, `crispr_pipeline_expression`, `crispr_pipeline_tf_benchmark`, `crispr_pipeline_dashboard`, `crispr_pipeline_samplesheet`, `crispr_pipeline_portal_samplesheet`, `crispr_pipeline_portal_download`, `crispr_pipeline_chunk`, `crispr_pipeline_nextflow`, `crispr_pipeline_setup`.

## 2nd CRISPR Jamboree analyses (`crispr-jamboree2`)

A port of [IGVF-CRISPR/CRISPR-JAMBOREE](https://github.com/IGVF-CRISPR/CRISPR-JAMBOREE), the 2024 IGVF CRISPR Jamboree. Each of its single-cell notebooks and scripts is a subcommand here. Guide counting (a STARsolo pseudo-genome, or a pure-Python FASTQ counter with whitelist and UMI 1-mismatch correction), seqspec indexing, guide assignment (UMI threshold, CLEANSER Stan mixtures fitted by MAP, sceptre mixture), the six inference modules on the shared MuData format (R wilcoxon, glm.nb, scanpy wilcoxon / t-test / t-test_overestim_var, a sceptre re-implementation, PerTurbo when installed), the DESeq2-based Perturb-seq simulator, and the evaluation notebooks (AUPRC / AUROC / effect-size correlation, clustergram, volcano, IGV tracks, element-gene network).

On the upstream's own Gasperini output MuDatas, the scanpy, glm.nb and wilcoxon ports reproduce the published p-values and fold changes to 1e-5 or better. Three upstream bugs are fixed, and `--upstream-compat` brings each one back: glm.nb stored exp(b1) as `l2fc`, the simulator wrote the unsimulated object, and the evaluation scored pairs by the raw p-value.

```bash
igvfagent crispr-jamboree2 inspect --mudata gasperini_inference_input.h5mu
igvfagent crispr-jamboree2 infer --mudata gasperini_inference_input.h5mu --methods scanpy-wilcoxon negbinom sceptre --side left
igvfagent crispr-jamboree2 simulate --mudata gasperini_inference_input.h5mu --perturb ENSG00000136856:candidate_enh_3:0.5 --n-null-pairs 50
igvfagent crispr-jamboree2 run --mudata simulation_output.h5mu --label sim
igvfagent crispr-jamboree2 count-guides --guides guides.xlsx --r1 R1.fastq.gz --r2 R2.fastq.gz --whitelist 737K-august-2016.txt
igvfagent crispr-jamboree2 selftest --no-plots
```

Agent tools: `crispr_jamboree2_inspect`, `crispr_jamboree2_guide_reference`, `crispr_jamboree2_count_guides`, `crispr_jamboree2_seqspec_index`, `crispr_jamboree2_assign_guides`, `crispr_jamboree2_infer`, `crispr_jamboree2_simulate`, `crispr_jamboree2_evaluate`, `crispr_jamboree2_volcano`, `crispr_jamboree2_network`, `crispr_jamboree2_run`.

## 3rd CRISPR Jamboree single-cell pipeline (`crispr-jamboree3`)

A stage-by-stage port of [IGVF-CRISPR/CRISPR-jamboree3](https://github.com/IGVF-CRISPR/CRISPR-jamboree3), the September 2024 snapshot of the IGVF single-cell Perturb-seq Nextflow pipeline demonstrated at the third jamboree. Every helper script is a subcommand:

- the configuration-form parser, which writes the shipped `pipeline_input.config` byte for byte;
- seqspec checks and parsing, and the kb mapping commands;
- AnnData concatenation and scRNA QC filtering;
- MuData assembly, Scrublet doublets and GMM hashing demultiplexing;
- CLEANSER / sceptre guide assignment, pairs-to-test preparation and sceptre inference;
- the volcano / network / IGV evaluation and the HTML dashboard.

`run` chains them from count matrices. External binaries (kb, GMM-Demux, cmdstan, the sceptre and PerTurbo packages) are re-implemented, or wrapped and run only when installed. Upstream quirks are corrected by default, and `--upstream-compat` reproduces each one: a pooled CLEANSER fit, positionally assigned guide metadata, `targeting` = TRUE for non-targeting guides, gene symbols in the GTF pairs, and the unused mouse mito prefix.

```bash
igvfagent crispr-jamboree3 configure --config-table configuration.csv
igvfagent crispr-jamboree3 seqspec-parse --yaml multiseq_guide_utsw.yml --modality guide
igvfagent crispr-jamboree3 map --modality guide --seqspec-yaml multiseq_guide_utsw.yml --metadata utsw_guide_metadata_new.xlsx --fastqs "batch_a:R1.fastq.gz R2.fastq.gz"
igvfagent crispr-jamboree3 run --rna rna.h5ad --guide guide.h5ad --hashing hashing.h5ad --guide-metadata utsw_guide_metadata_new.xlsx --gtf gencode.v46.annotation.gtf.gz --pairs user_pairs_to_test.csv
igvfagent crispr-jamboree3 dashboard --mudata inference_mudata.h5mu --gene-ann rna.h5ad --guide-ann guide.h5ad
igvfagent crispr-jamboree3 selftest --no-plots
```

Agent tools: `crispr_jamboree3_configure`, `crispr_jamboree3_seqspec_check`, `crispr_jamboree3_seqspec_parse`, `crispr_jamboree3_map`, `crispr_jamboree3_concat`, `crispr_jamboree3_preprocess`, `crispr_jamboree3_create_mudata`, `crispr_jamboree3_doublets`, `crispr_jamboree3_demultiplex`, `crispr_jamboree3_assign_guides`, `crispr_jamboree3_prepare_inference`, `crispr_jamboree3_infer`, `crispr_jamboree3_evaluate`, `crispr_jamboree3_dashboard`, `crispr_jamboree3_run`.

## 4th CRISPR Jamboree inference benchmark (`crispr-jamboree4`)

A port of [IGVF-CRISPR/CRISPR-Jamboree_2025](https://github.com/IGVF-CRISPR/CRISPR-Jamboree_2025), the 2025 jamboree task. It starts the IGVF CRISPR pipeline from the guide-inference checkpoint, runs sceptre and PerTurbo on a dataset's pairs to test (H9 or WTC11), merges the two methods' results into the pipeline outputs, and scores how well each method separates the control pairs (AUPRC / AUROC, with a precision-recall and ROC grid).

The task bundle on Dropbox has been deleted. The steps were therefore rebuilt from the CRISPR_Pipeline revision the bundle ran (February 2025), and the evaluation script was reconstructed from the notebook's recorded output. A glm.nb comparator is included. PerTurbo runs only when its package is installed.

```bash
igvfagent crispr-jamboree4 prepare --mudata h9_guide_assignment.h5mu --pairs pairs_to_test_h9.csv
igvfagent crispr-jamboree4 infer --mudata mudata_inference_input.h5mu --methods sceptre negbinom
igvfagent crispr-jamboree4 merge --mudata mudata_inference_input.h5mu --sceptre-results test_results.sceptre.tsv --perturbo-results perturbo_test_results.tsv
igvfagent crispr-jamboree4 evaluate --results H9=h9/inference_mudata.h5mu WTC11=wtc11/inference_mudata.h5mu
igvfagent crispr-jamboree4 run --mudata h9_guide_assignment.h5mu --pairs pairs_to_test_h9.csv --label h9
igvfagent crispr-jamboree4 selftest --no-plots
```

Agent tools: `crispr_jamboree4_prepare`, `crispr_jamboree4_infer`, `crispr_jamboree4_merge`, `crispr_jamboree4_evaluate`, `crispr_jamboree4_run`.

## Single-cell CRISPR differential expression (`sc-crispr-de`)

Per-guide differential expression for **Perturb-seq / single-cell CRISPRi
screens**. Every guide is tested independently: the cells carrying it
against the cells with **no detected guide**, one negative-binomial GLM
per gene, then the per-guide p-values are aggregated to a score per
target-gene pair.

Method from
[Gersbachlab-Bioinformatics/sc-crispr-de](https://github.com/Gersbachlab-Bioinformatics/sc-crispr-de)
(MIT). Clean-room reimplementation, no source copied — the upstream is R
(Seurat, `MASS::glm.nb`, brglm2, GenomicRanges) driven by a SLURM array,
and this container has neither R nor a scheduler, so a port was never an
option. What is reproduced is the method.

```bash
pip install 'igvfagent[analysis]'    # + statsmodels, scanpy, anndata

# 1) 10x matrices -> filtered checkpoint with per-cell guide calls.
#    Omit --sgrna when guide features are inside the GEX matrix
#    (CRISPR Guide Capture); CellRanger output is read directly.
igvfagent sc-crispr-de prepare --gex filtered_feature_bc_matrix/ \
    --sgrna sgrna_matrix/ --label screen1

# 2) Per-guide NB GLM. The compute-heavy step: upstream runs one SLURM
#    array task per guide, here guides run across a process pool.
#    --gene-whitelist restricts to a locus when that is the question.
igvfagent sc-crispr-de test --checkpoint <run>.h5ad \
    --latentvar nCount_RNA --workers 8 --label screen1

# 3) Gene level, calibrated against the non-targeting guides.
igvfagent sc-crispr-de aggregate --per-guide <run>_per_guide.tsv \
    --nontargeting NTC_1,NTC_2,... --label screen1

# or all three:
igvfagent sc-crispr-de pipeline --gex ... --sgrna ... --label screen1
```

**Validated on planted signal.** Synthetic 1,300-cell screen, three
targeting guides and twelve NTCs, two genes knocked down by a known
factor:

| planted | expected log2FC | recovered |
|---|---|---|
| GATA1 ×0.25 | −2.00 | −1.89 / −2.22 / −2.03 |
| HBB ×0.40 | −1.32 | −1.18 / −1.25 / −1.38 |

Both reach FDR < 0.05 at gene level; the NTC guides do not.

**Three deliberate differences from upstream**, because a
reimplementation that quietly diverges is worse than none. *Dispersion:*
`MASS::glm.nb` profiles θ by ML with IRLS alternation, statsmodels
estimates α = 1/θ by direct ML — same parameter, different optimiser, so
the last digits differ. *Bias reduction:* `--apply-bias-reduction`
(brglm2, Firth-type) has no equivalent here, so it is **not implemented**
rather than approximated; separated guides are reported as `separated`.
*Aggregation:* upstream calls the external FRACTEL package; this is
α-RRA (Kolde 2012, as in MAGeCK) with an empirical null from the
non-targeting guides. Same family, not the same numbers — the columns are
`rra_*`, never `FRACTEL_*`, so no reader can confuse them. Without enough
NTCs the null is *assumed* rather than measured, and the output says so.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
