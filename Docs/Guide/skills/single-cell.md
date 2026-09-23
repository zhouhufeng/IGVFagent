# Single-cell, multiome and spatial

scRNA-seq, multiome, SPLiT-seq, SHARE-seq, snMCT-seq, cell hashing, Tabula Sapiens, the IGVF single-cell pipeline, principal pseudobulks and Spatial-ATAC-Hi-C.

**On this page**

- [Single-cell, multiome, specialized assays](#single-cell-multiome-specialized-assays)
- [Cross-source multiome survey](#cross-source-multiome-survey)
- [Parse SPLiT-seq pipeline](#parse-split-seq-pipeline)
- [snMCT-seq — RNA and methylation from the same nuclei](#snmct-seq--rna-and-methylation-from-the-same-nuclei)
- [SHARE-seq joint scATAC + scRNA QC](#share-seq-joint-scatac--scrna-qc)
- [Single-cell analysis (UMAP / t-SNE / Leiden / markers)](#single-cell-analysis-umap--t-sne--leiden--markers)
- [MULTI-seq / Cell Hashing demultiplexing](#multi-seq--cell-hashing-demultiplexing)
- [Tabula Sapiens 2.0 (reference human cell atlas)](#tabula-sapiens-20-reference-human-cell-atlas)
- [IGVF single-cell pipeline (`igvf-sc-pipeline`)](#igvf-single-cell-pipeline-igvf-sc-pipeline)
- [Principal pseudobulks (`principal-pseudobulks`)](#principal-pseudobulks-principal-pseudobulks)
- [Spatial-ATAC-Hi-C (spatial 3D genome + chromatin accessibility)](#spatial-atac-hi-c-spatial-3d-genome--chromatin-accessibility)

## Single-cell, multiome, specialized assays

```bash
python3 Scripts/single_cell_data_skills.py smoke --skill all --source encode --limit 5
python3 Scripts/single_cell_data_skills.py manifest --skill scrna --source both --limit 25
python3 Scripts/single_cell_data_skills.py download-examples
python3 Scripts/single_cell_data_skills.py analyze-examples --max-cells 12000
python3 Scripts/single_cell_data_skills.py write-playbook

python3 Scripts/multiome_10x_pipeline.py retrieve --count 20 --fetch-file-details
python3 Scripts/multiome_10x_pipeline.py retrieve \
  --count 20 --fetch-file-details --download-policy all --max-download-gb 30
python3 Scripts/multiome_10x_pipeline.py process-local \
  --file-manifest Data/Manifests/Multiome10x/<files.csv> \
  --download-manifest Data/Manifests/Multiome10x/<download_manifest.csv>
python3 Scripts/multiome_10x_pipeline.py write-playbook
python3 Scripts/multiome_research_demo.py

python3 Scripts/igvf_specialized_data_skills.py smoke --skill all --limit 5
python3 Scripts/igvf_specialized_data_skills.py manifest --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py download-plan --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py write-playbook
```

## Cross-source multiome survey

Surveys six independent public sources for single-cell multiome data
(10x Multiome / SHARE-seq / single-nucleus multiome / SNARE-seq /
Paired-Tag), classifies every file by a training-relevant `kind`
(`matrix_rna`, `matrix_atac`, `fragments`, `peaks`, `annotations`,
`alignments`, `index`, `raw_reads`, `other`), and produces a unified
manifest with a size-capped downloader.

| source        | what it covers                                                       |
|---------------|----------------------------------------------------------------------|
| **IGVF**      | IGVF Portal AnalysisSets/MeasurementSets (10x multiome + SHARE-seq). |
| **ENCODE**    | ENCODE Experiments + Series tagged 10x multiome / SHARE-seq.          |
| **GEO**       | NCBI Gene Expression Omnibus Series via E-utilities + FTP listing.    |
| **CELLxGENE** | CZI CELLxGENE Discover curated H5AD collections.                      |
| **HCA**       | Human Cell Atlas Data Portal projects (Azul `/index/projects`).        |
| **Zenodo**    | Zenodo records with multiome keywords in title / description / files. |

```bash
# Run all six sources in one shot, then build the unified manifest
python3 Scripts/multiome_survey.py survey-all --limit 30 --fetch-files
python3 Scripts/multiome_survey.py manifest --label v1

# Training-relevant slice only (RNA + ATAC + fragments + annotations), 20 GB cap
python3 Scripts/multiome_survey.py download \
    --only matrix_rna,matrix_atac,fragments,annotations \
    --max-download-gb 20

# Per-source surveys
python3 Scripts/multiome_survey.py survey-igvf      --fetch-files
python3 Scripts/multiome_survey.py survey-encode    --fetch-files
python3 Scripts/multiome_survey.py survey-geo       --limit 50 --fetch-files
python3 Scripts/multiome_survey.py survey-cellxgene --fetch-files
python3 Scripts/multiome_survey.py survey-hca       --fetch-files
python3 Scripts/multiome_survey.py survey-zenodo    --per-query 25 --fetch-files

# Inventory on-disk downloads + skill / overview docs
python3 Scripts/multiome_survey.py inventory
python3 Scripts/multiome_survey.py write-playbook
python3 Scripts/multiome_survey.py write-overview
```

A live `survey-all` (May 2026) indexed **401 datasets / 3,837 files /
14.19 TB total** across the six sources; the digest is committed at
`Docs/MultiomeSurvey/LATEST_RUN_SUMMARY.md`. The companion
`Docs/MultiomeSurvey/SOURCES_OVERVIEW.md` covers what each source offers
and points to additional repositories (Allen Brain Cell Atlas, Broad
Single Cell Portal, Synapse, dbGaP, EGA, ArrayExpress, NeMO Archive,
figshare, DDBJ, Terra) where multiome data also lives but isn't
auto-queried by the skill.

## Parse SPLiT-seq pipeline

End-to-end pipeline for Parse-Biosciences combinatorial-barcoding snRNA-seq
(SPLiT-seq) datasets. Handles the SPLiT-seq quirks that generic single-cell
skills don't: multiplexed sub-pools, tens of donors per pool, MatrixMarket
tarball delivery, and per-strain analytical comparisons (the value of the
Mortazavi 8-cube founder atlas in IGVF).

```bash
# Discover SPLiT-seq AnalysisSets in the IGVF Portal
python3 Scripts/splitseq_pipeline.py retrieve --limit 50 --label splitseq_corpus

# Per-pool / per-donor manifest for one or more accessions
python3 Scripts/splitseq_pipeline.py manifest \
  --accessions IGVFDS3222WCZH,IGVFDS6290WNNH --label demo

# Pull files under a size cap, then load into AnnData
python3 Scripts/splitseq_pipeline.py download \
  --manifest Data/Manifests/SPLiTseq/<files_per_pool.csv> --max-download-gb 5
python3 Scripts/splitseq_pipeline.py process-local --label demo

# Full pipeline: QC -> normalize -> Harmony integrate -> UMAP -> Leiden ->
# auto-annotate against bundled mouse marker panels
python3 Scripts/splitseq_pipeline.py analyze \
  --input Data/Cache/SPLiTseq/demo.h5ad --tissue gonad --label demo
python3 Scripts/splitseq_pipeline.py plot \
  --input Data/Cache/SPLiTseq/demo_processed.h5ad --tissue gonad --label demo

# Per-cell-type strain DEG (after donor demultiplexing)
python3 Scripts/splitseq_pipeline.py compare-strains \
  --input Data/Cache/SPLiTseq/demo_processed.h5ad --label 8cube_strain_DEG

# Scaffold a souporcell or vireo demultiplexing run (the agent emits the
# shell script; the demultiplexer itself runs locally)
python3 Scripts/splitseq_pipeline.py demux-script --tool souporcell --n-donors 8
python3 Scripts/splitseq_pipeline.py write-playbook
```

Bundled mouse marker panels: `gonad`, `adrenal`, `brain`, `liver`, `heart`,
`kidney`, `muscle` — covers the full Mortazavi 8-cube founder atlas tissues.

## snMCT-seq — RNA and methylation from the same nuclei

```bash
igvfagent mct discover IGVFDS4826YNLK
igvfagent mct analyze  IGVFDS4826YNLK --label mct_run
```

Both halves are analysed and cross-compared (adjusted Rand index between the
RNA and methylation clusterings, plus a biological-vs-technical covariate
check). Routing this as a transcript assay analyses the RNA and silently
drops the methylation — the half the assay exists for, across 933 Portal
datasets. Methylation ratios are **not** log-normalised, which is the step
that makes a methylome look like an expression matrix.

## SHARE-seq joint scATAC + scRNA QC

Clean-room reimplementation of the QC algorithms in
[broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline)
(MIT, 2021). Consumes IGVF Portal SHARE-seq AnalysisSet deposits
(`fragments.bed[.gz]` + `h5ad`) directly. See
[`Docs/Skills/SHARESEQ_ANALYSIS_SKILLS.md`](../../Skills/SHARESEQ_ANALYSIS_SKILLS.md).

```bash
# Discover IGVF Portal SHARE-seq datasets
igvfagent share pull-portal --limit 50 --label survey

# Round-1/2/3 24-mer barcode demultiplex on a gz FASTQ
igvfagent share demultiplex-bcs --fastq R2.fastq.gz \
    --whitelist whitelist_24mer.txt --out R2.tagged.fastq.gz \
    --shift-correct --label demo

# Per-barcode ATAC QC (TSS enrichment, FRIP)
igvfagent share fragment-qc --fragments fragments.bed.gz \
    --tss-bed tss.bed --peaks-bed peaks.bed --label demo

# Per-barcode RNA QC (UMIs, genes, %MT) from h5ad
igvfagent share rna-qc --h5ad rna.h5ad --label demo

# Joint cell call (Ma 2020 thresholds)
igvfagent share joint-qc --rna-qc <ts>_demo_rna_qc.tsv \
    --atac-qc <ts>_demo_atac_qc.tsv --label demo

# Jaccard multiplet detection
igvfagent share multiplet-detect --fragments fragments.bed.gz --label demo
```

## Single-cell analysis (UMAP / t-SNE / Leiden / markers)

End-to-end Scanpy-driven single-cell workflow: **QC → normalize → HVG
→ PCA → k-NN → UMAP → t-SNE → Leiden clustering → marker-gene DE →
publication figures.** Closes the gap where IGVFagent could discover
and download counts matrices but had to hand off to "use Scanpy or
Seurat" for the actual visualization. Auto-detects input format
(`.h5ad`, 10x `.h5`, 10x `.mtx`, CSV/TSV). Playbook:
[`Docs/Skills/SINGLECELL_ANALYSIS_SKILLS.md`](../../Skills/SINGLECELL_ANALYSIS_SKILLS.md).

```bash
# 1) Full pipeline — one shot, all stages
igvfagent sc-analyze pipeline --input counts.h5ad --label demo \
    --min-genes 200 --min-cells 3 --max-mito 20 \
    --n-hvg 2000 --n-pcs 50 --resolution 1.0 \
    --sample-col sample \
    --highlight-genes CD3D,CD8A,MS4A1,LYZ

# 2) Granular steps — each saves a checkpoint processed.h5ad so the
#    next step can resume from any prior output
igvfagent sc-analyze qc        --input counts.h5ad --label k562 \
    --min-genes 200 --max-mito 20
igvfagent sc-analyze normalize --input processed.h5ad --n-hvg 2000
igvfagent sc-analyze pca       --input processed.h5ad --n-pcs 50
igvfagent sc-analyze umap      --input processed.h5ad --n-neighbors 15
igvfagent sc-analyze tsne      --input processed.h5ad
igvfagent sc-analyze cluster   --input processed.h5ad --resolution 1.0
igvfagent sc-analyze markers   --input processed.h5ad --n-top 25

# 3) Re-render any embedding coloured by a gene or metadata column
igvfagent sc-analyze plot-embedding --input processed.h5ad \
    --embedding umap --color leiden,APOE,TREM2,LDLR
```

Outputs land under `Docs/SingleCell/<timestamp>_<label>/` with the
resumable `processed.h5ad`, `markers.csv`, a markdown report, and
PNGs (QC violins, PCA scree, UMAP and t-SNE coloured by Leiden
cluster / sample / auto-picked top markers, top-3 marker heatmap).
The agent runtime exposes `sc_pipeline`, `sc_qc`, `sc_umap`,
`sc_cluster`, and `sc_plot_embedding` as tools so a single
`igvfagent ask "give me a UMAP of GSE131907 with NK markers"` will
drive the full chain. Live smoke test on the Scanpy PBMC3k tutorial
dataset finishes in ~20 s on 2,700 cells × 32,738 genes and recovers
the canonical T-cell / B-cell / monocyte clusters.

## MULTI-seq / Cell Hashing demultiplexing

Python port of **deMULTIplex2** (Zhu et al., *Nat Methods* 2024 —
[Gartner-Lab/deMULTIplex2](https://github.com/Gartner-Lab/deMULTIplex2)),
the v2 classifier for **MULTI-seq** ([McGinnis et al., *Nat Methods*
2019](https://www.nature.com/articles/s41592-019-0433-8) — the
original method using lipid-tagged 8-nt sample barcodes anchored in
the cell or nuclear membrane via two LMOs;
[Sigma-Aldrich technical
article](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/sequencing/multi-seq-sample-multiplexing-single-cell-analysis-sequencing)
covers the LMO001 reagent kit and protocol). Assigns every cell in a
multiplexed scRNA-seq run to a sample of origin, flags doublets, and
flags negatives — the missing piece between `sc-analyze` (counts →
clusters) and downstream sample-stratified analysis. Playbook:
[`Docs/Skills/MULTISEQ_ANALYSIS_SKILLS.md`](../../Skills/MULTISEQ_ANALYSIS_SKILLS.md).

```bash
# 1) Generate a synthetic 2,000-cell × 6-tag matrix for smoke testing
igvfagent multiseq simulate --n-cells 2000 --n-tags 6 \
    --doublet-rate 0.08 --negative-rate 0.05 --label smoke

# 2) Demultiplex any tag-count matrix (.h5ad / 10x .h5 / .csv / .tsv)
igvfagent multiseq demultiplex --input tag_counts.csv \
    --label demo --prob-cut 0.5 --residual-type rqr

# 3) Standalone diagnostics
igvfagent multiseq histogram --input tag_counts.csv
igvfagent multiseq heatmap   --input tag_counts.csv \
    --calls classifications.csv

# 4) End-to-end with optional accuracy table vs ground truth
igvfagent multiseq pipeline --input tag_counts.csv \
    --ground-truth ground_truth.csv --label end_to_end
```

The algorithm fits, for each sample tag *j* independently, a
two-component negative-binomial mixture via EM:

  * `fit0`  `bc_umi ~ log(tt_umi)`                (off-target background)
  * `fit1`  `(tt_umi - bc_umi) ~ log(tt_umi)`     (background tags in positives)

where `bc_umi` is the focal-tag count and `tt_umi` is the per-cell
total tag UMI. EM is initialized from a cosine-similarity cut, fit
via `statsmodels` NB2 (matching `MASS::glm.nb`), and iterates until
the log-likelihood is stable (≤ 10 iter, tol 1e-3). Cells with
posterior `P(positive) > 0.5` for one tag are **singlets**; ≥2 →
**multiplets**; 0 → **negatives**. Randomized quantile residuals
(Dunn & Smyth 1996) are returned per tag for downstream QC.

The agent runtime exposes `multiseq_demultiplex`, `multiseq_pipeline`,
and `multiseq_simulate` as tools. Outputs go to
`Docs/MultiSeq/<timestamp>_<label>/` with `classifications.csv`,
`posteriors.csv`, `residuals.csv`, `tag_coefficients.csv`, and PNGs
(faceted log-x histograms, mean-tag-count heatmap by call group,
per-tag 4-panel scatter diagnostics). **Live smoke test** on
2,000 simulated cells × 6 tags with 8% doublets + 5% negatives:
**99.85% overall accuracy**, all 159 multiplets correctly flagged,
all 98 negatives correctly flagged. FASTQ → tag-count alignment
(the R package's `readTags` / `alignTags` step) is **not** ported —
hand this skill a count matrix produced by Cell Ranger's Feature
Barcoding workflow or the original `deMULTIplex` aligner.

## Tabula Sapiens 2.0 (reference human cell atlas)

Retrieval and figure-level reproduction of the **Tabula Sapiens 2.0**
human cell atlas (Tabula Sapiens Consortium, *Cell* 2026): 1,136,218
cells across 28 tissues from 24 donors, annotated to 182 fine and 38
broad cell types. Clean-room reimplementation of the analyses in
[czbiohub-sf/tabula-sapiens](https://github.com/czbiohub-sf/tabula-sapiens)
(BSD-3-Clause). Playbook:
[`Docs/Skills/TABULA_SAPIENS_SKILL.md`](../../Skills/TABULA_SAPIENS_SKILL.md).

**Which data routes are open, measured rather than assumed.** The AWS
open-data listing implies the raw reads are public. They are not: the
bucket is *listable* but every `GetObject` returns 403 AccessDenied, on
v1 and v2 alike, with or without `x-amz-request-payer`. That is the data
transfer agreement the paper describes for donor genetic privacy.
`tabula s3-manifest` probes this on each run and reports what it found.
The bucket totals **103.16 TB** — 43.2 TB BAM, 32.1 TB FASTQ, 27.6 TB
STAR intermediates, and just 0.32 TB of count matrices. The open routes
carry the science: figshare (57 GB of processed h5ads), GEO GSE306755
(count matrices + the full 1.1M-cell metadata), and CELLxGENE.

```bash
pip install 'igvfagent[analysis]'

igvfagent tabula pull-geo            # 41 MB metadata — enough for Figure 1
igvfagent tabula status              # local state vs the paper's own numbers
igvfagent tabula overview            # Figure 1: donors, tissues, composition

igvfagent humantfs build-db          # 1,639 curated TFs (Lambert 2018)
igvfagent tabula pull-figshare       # the 57 GB atlas, resumable
igvfagent tabula tf-matrix           # gene x cell-type means
igvfagent tabula tf-specificity --matrix <run>/mean_expression.npz
igvfagent tabula tf-enrichment --tau-table <run>/tf_tau.tsv
igvfagent tabula tf-regulons --tissues Lung,Heart   # SCENIC, clean-room
igvfagent tabula senescence          # Figure 4: CDKN2A+ MKI67- burden
igvfagent tabula sex-de              # Figure 5: pseudobulk sex differences
igvfagent tabula donors --min-age 60 # Figure 6: the ChatTS core, offline
```

**What it reproduces** — `Benchmarks/quake2026_tabula_sapiens` scores
**26/26**:

| Quantity | IGVFagent | Paper |
|---|---:|---:|
| Cells (droplet / FACS) | 1,136,218 (1,093,048 / 43,170) | identical |
| Donors / tissues / fine types / populations | 24 / 28 / 182 / 701 | identical |
| Droplet fine / broad cell types | 175 / 38 | identical |
| Donor ages, sex, age groups | 22–74, 11M/13F, 7/11/6 | identical |
| TFs with zero expression atlas-wide | **2: SHOX, ZBED1** | 2: SHOX, ZBED1 |
| TF specific / non-specific (τ > 0.85, 175 cell types) | 882 / 741 | 890 / 745 |
| GO terms for non-specific TFs (padj < 0.02) | 70 | 69 |

Plus biology that is not merely counting: FOXP3 peaks in regulatory
T cells, FOXI1 in ionocytes, germ-cell TFs in spermatogenic cells, and
all eight of the paper's named ubiquitous TFs fall below τ 0.85.

**A discrepancy worth knowing about.** Figure 3's Methods say the
non-specific TFs were tested "against a background of all 1635
transcription factors". Running both backgrounds on identical input:
the Enrichr default gives 70 terms at padj < 0.02 (paper: 69), the
explicit TF background gives **zero**. gseapy 0.10.5, the cited version,
forwards `background` to the Enrichr API, which ignores it. The
published figure came from the default background. `tf-enrichment`
defaults to that behaviour and offers `--tf-background` for what the
Methods describe.

**Scope.** Starts from the processed deposits, the boundary `share` and
`spatial-hic` also draw. STAR 2.7.11b and CellRanger 7.0.1 are
orchestrated, not reimplemented — they are large C/C++ aligners, and the
published h5ads already carry their output alongside DecontX and scVI
results, which is what lets the clean-room reimplementations here be
*validated* rather than merely run.

License boundary: **Apache-2 throughout**. pySCENIC (GPL-3) and edgeR
(GPL) are reimplemented, not imported.

## IGVF single-cell pipeline (`igvf-sc-pipeline`)

A Python port of the IGVF uniform single-cell multiome pipeline ([IGVF/single-cell-pipeline](https://github.com/IGVF/single-cell-pipeline), with the chromap and kallisto|bustools wrappers from IGVF/atomic-workflows v1.1). Every step of the WDL is a subcommand: SHARE-seq barcode correction (exact / 1-mismatch / +-1 bp shift, poly-G QC), dovetail trimming, barcode-orientation detection, the 10x multiome ATAC<->RNA barcode map, chromap and kb count command construction with subpool handling, chromap log metrics, bulk and per-barcode (ArchR) TSS enrichment, RNA metrics from kb nac layers, joint RNA x ATAC cell calling, and barcode-rank elbow/knee detection with a re-implementation of R's `smooth.spline(spar=1)`.

chromap and kb are optional: when they are not installed, the tool prints the exact command the pipeline would run and exits 0. Every QC step then runs on the h5ad and fragment files those tools produce. `pipeline` takes a Cromwell inputs JSON and runs the whole plan. The Synapse and IGVF-portal helpers read credentials from the environment only, and do a dry run unless asked to execute.

```bash
igvfagent igvf-sc-pipeline correct-fastq --read1 R1.fq.gz --read2 R2.fq.gz --whitelist bc24.txt --sample-type ATAC --prefix lib1
igvfagent igvf-sc-pipeline tss-enrichment --fragments lib1.fragments.tsv.gz --regions tss.bed --prefix lib1
igvfagent igvf-sc-pipeline rna-qc-metrics --h5ad lib1.rna.h5ad --kb-workflow nac --subpool SP1
igvfagent igvf-sc-pipeline joint-qc --rna-metrics rna_barcode_metadata.tsv --atac-metrics lib1.tss_enrichment_barcode_stats.tsv --pkr SP1
igvfagent igvf-sc-pipeline pipeline --inputs-json inputs.json --rna-h5ad lib1.rna.h5ad --fragments lib1.fragments.tsv.gz --tss-bed tss.bed
igvfagent igvf-sc-pipeline selftest --no-plots
```

Agent tools: `igvf_sc_barcode_revcomp_detect`, `igvf_sc_correct_fastq`, `igvf_sc_trim_fastq`, `igvf_sc_tss_enrichment`, `igvf_sc_snapatac2_tsse`, `igvf_sc_rna_qc_metrics`, `igvf_sc_mtx_to_h5ad`, `igvf_sc_modify_barcode_h5`, `igvf_sc_joint_qc`, `igvf_sc_barcode_rank`, `igvf_sc_atac_qc_plots`, `igvf_sc_rna_qc_plots`, `igvf_sc_insert_size_hist`, `igvf_sc_tenx_barcode_map`, `igvf_sc_log_atac`, `igvf_sc_kb_count`, `igvf_sc_kb_index`, `igvf_sc_chromap_align`, `igvf_sc_chromap_index`, `igvf_sc_subpool_fragments`, `igvf_sc_genome_tsv`, `igvf_sc_check_inputs`, `igvf_sc_sample_fastqs`, `igvf_sc_portal_download`, `igvf_sc_synapse_manifest`, `igvf_sc_synapse_upload`, `igvf_sc_synapse_annotations`, `igvf_sc_html_report`, `igvf_sc_pipeline`.

## Principal pseudobulks (`principal-pseudobulks`)

A Python port of [EngreitzLab/generate-principal-pseudobulks](https://github.com/EngreitzLab/generate-principal-pseudobulks) (MIT, pinned at 80222cc), the pipeline that turns IGVF multiome primary pseudobulks into released principal pseudobulks. It covers every step: the per-cell QC datatable, threshold exploration (cells dropped per threshold in total and alone), the final QC filter with its QC guide, threshold record, per-subsample metrics and PNG QC figures, ATAC fragment filtering (sorted, bgzip + tabix), the gene-symbol RNA count matrix (GENCODE 43 GTF, standard chromosomes, duplicate symbols summed, hard failure on unmatched IDs), the per-dataset config table, and the Snakefile driver. It runs without Snakemake, R, bedops or htslib. `validate` checks outputs against both upstream file specs.

It starts from an IGVF accession rather than prebuilt directories. `fetch` walks the Portal lineage, processed-first, to the principal analysis set's cell annotations and each uniform-pipeline lane's fragments and h5ad, within a download budget. `build-pseudobulks` then builds the primary pseudobulk layout, with lane-suffixed barcodes, ATAC-to-RNA multiome barcode translation and per-cell QC. The per-dataset config table feeds `sce2g-pipeline run --cluster-config` directly, and `sce2g-prep` adds the tagAlign and RNA pseudobulk TPM.

```bash
igvfagent principal-pseudobulks fetch --accession IGVFDS1244UUGQ --dry-run --max-gb 20
igvfagent principal-pseudobulks build-pseudobulks --annotations cells.tsv --manifest portal_manifest.tsv --dataset igvf1 --tss tss.bed --peaks peaks.bed
igvfagent principal-pseudobulks explore --meta k562_per_cell_qc.tsv --sets "--tss-min 3" "--tss-min 5 --pct-mt-max 20" --show-subsamples
igvfagent principal-pseudobulks run --config config_QC_pseudobulks.yaml --auto-qc
igvfagent sce2g-pipeline run --cluster-config multiome_data/config/tables/igvf1_config.tsv
```

Agent tools: `principal_pseudobulks_fetch`, `principal_pseudobulks_build`, `principal_pseudobulks_qc_datatable`, `principal_pseudobulks_explore`, `principal_pseudobulks_qc_filter`, `principal_pseudobulks_filter_atac`, `principal_pseudobulks_filter_rna`, `principal_pseudobulks_package_rna`, `principal_pseudobulks_config_table`, `principal_pseudobulks_run`, `principal_pseudobulks_validate`, `principal_pseudobulks_sce2g_prep`, `principal_pseudobulks_selftest`.

## Spatial-ATAC-Hi-C (spatial 3D genome + chromatin accessibility)

Spatially resolved **co-profiling of genome folding and chromatin
accessibility on tissue slides**, after [Wang, Wang, Wang, Youngblood
et al., *Nature Methods* 2026](https://doi.org/10.1038/s41592-026-03217-4)
(GEO **GSE307620**). The assay lays a 50 x 50 microfluidic barcode grid
over a section, giving **2,500 spatial pixels** that each carry a Hi-C
contact set *and* an ATAC fragment set from the same molecules — so
compartments, loops, copy number and gene activity all keep their
tissue coordinates.

Clean-room reimplementation of the computational pipeline in
[wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C)
(MIT), which delegates to scHiCluster, Higashi/cooltools, SnapATAC2 and
NeoLoopFinder. Every one of those algorithms is re-derived here from its
published description; none is imported or vendored. Playbook:
[`Docs/Skills/SPATIAL_ATAC_HIC_SKILLS.md`](../../Skills/SPATIAL_ATAC_HIC_SKILLS.md).

**Scope.** This skill consumes **processed deposits** — contact tables,
`fragments.tsv.gz`, spatial positions — the same boundary `share` draws.
It reads both 4DN `.pairs` and the tabix-indexed contact TSVs GEO
actually ships (`*.hic.fragments.sorted.header.tsv.gz`), normalising
column-name synonyms and picking pixel ids out of compound filenames.
Read alignment is *not* reimplemented: it is the one upstream step with
no tractable clean-room form, and the tools involved
([runHiC](https://github.com/XiaoTaoWang/HiC_pipeline) and
[Trim Galore](https://github.com/FelixKrueger/TrimGalore)) are GPL-3.0,
so they stay external tools rather than dependencies. GSE307620
publishes aligned pairs and fragments directly.

```bash
pip install 'igvfagent[analysis]'   # numpy + scipy + matplotlib

# 0) What does the series deposit?
igvfagent spatial-hic pull-geo --gse GSE307620 \
    --download 'pairs|fragments|positions'

# 1) One barcoded pairs file -> the 50x50 pixel grid.
#    Check `assigned_fraction`: near-zero means --layout is the other way.
igvfagent spatial-hic pixel-demux --pairs sample.pairs.gz \
    --barcode-a barcodes_A.txt --barcode-b barcodes_B.txt --label mouse_R6

# 2) Per-pixel QC. Paper reference band: 25,343-58,403 total contacts,
#    88.1-90.3% cis, 24-33.3% long-range (>=10 kb over cis).
igvfagent spatial-hic qc --pairs-dir <run>/pixels \
    --fragments fragments.tsv.gz --gene-model gencode.gtf.gz

# 3) Gene-level scores in both modalities, on a shared pixel key.
igvfagent spatial-hic gas --fragments fragments.tsv.gz \
    --gene-model gencode.gtf.gz \
    --barcode-a barcodes_A.txt --barcode-b barcodes_B.txt
igvfagent spatial-hic gad --pairs-dir <run>/pixels --gene-model gencode.gtf.gz

# 4) Structure: contact matrix, imputation, compartments, copy number.
#    A single pixel holds only tens of thousands of contacts, so impute
#    before reading structure off one. Paper resolutions: 100 kb for
#    compartments, 25 kb to draw, 10 kb for fine structure.
igvfagent spatial-hic matrix --pairs-dir <run>/pixels --chrom chr2 \
    --resolution 25000 --start 106000000 --end 118000000
igvfagent spatial-hic impute --pairs-dir <run>/pixels --chrom chr2 \
    --resolution 100000 --per-pixel
igvfagent spatial-hic compartment --pairs-dir <run>/pixels \
    --resolution 100000 --gene-model gencode.gtf.gz
igvfagent spatial-hic cnv --pairs-dir <run>/pixels \
    --resolution 5000000 --per-pixel --smooth      # spatial tumour clones

# 5) Loops: quantify at known anchors, ANOVA across clusters, APA.
igvfagent spatial-hic loops --pairs-dir <run>/pixels --bedpe loops.bedpe \
    --clusters clusters.tsv --apa

# 6) Any per-pixel value, back in tissue space.
igvfagent spatial-hic viz --table <run>/gene_activity_score.tsv \
    --column Satb2 --positions tissue_positions.csv
```

Every one of these eleven subcommands is also registered as an agent
tool (`spatial_hic_*`), so the orchestrator can select them from a plain
request — "impute a 25 kb map for chr2", "find the tumour clones in this
section" — rather than needing the CLI spelled out.

**Against the authors' own pipeline**
([wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C)):

| Upstream script | IGVFagent |
|---|---|
| `00.filter-linker.sh` | *not reimplemented* — raw-FASTQ stage |
| `01.bcsplit.sh` / `bcsplit.py` | `pixel-demux` (barcode-B+A offsets follow it) |
| `03.atac_alignment.sh` | *not reimplemented* — runHiC / Trim Galore are GPL-3.0, see above |
| `snapatac2_genecount_matrix.py` | `gas` |
| `Higashi_compartment.sh` | `compartment` |
| `schicluster_impute.sh` | `impute` |
| — | `qc`, `gad`, `matrix`, `cnv`, `loops`, `viz` (from the paper, not in the repo) |

So the boundary is alignment: everything the authors publish downstream
of it has an equivalent here, and the per-pixel QC, GAD scores, CNV,
loop ANOVA and tissue-space rendering that the paper describes but the
repo does not ship are implemented too.

Two things this deliberately does *not* do: it does not call loops de
novo (Peakachu is a trained model — bring its BEDPE), and its CNV bias
correction is a linear fit rather than NeoLoopFinder's Poisson GLM, so
treat copy number as concordant rather than bit-identical.

License boundary: **Apache-2 throughout**. No GPL runtime dependencies.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
