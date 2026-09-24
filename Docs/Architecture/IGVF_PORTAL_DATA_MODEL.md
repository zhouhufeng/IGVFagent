# IGVF Portal data model: how submissions link

How data, results and their context are connected on the IGVF Portal, and how
IGVFagent uses those links. The source is the 22 lab submission diagrams in
`Data/IGVF/DataModels/` (the DACC data-wrangler drafts for each lab's
pipeline), checked on 2026-09-23 against the live Portal schema
(`https://api.data.igvf.org/profiles/<type>.json`). Where a diagram and the
schema disagree, the schema wins; the differences are listed at the end.

IGVFagent implements this in two places:

- `igvfagent processed lineage <accession>` (`Scripts/portal_lineage.py`)
  walks every link below from one accession.
- `igvfagent processed discover` (`Scripts/portal_discover.py`) finds
  accessions by phenotype, tissue, gene or kind of data.

## Object types and what they hold

| Type | Holds | Typical `file_set_type` |
|---|---|---|
| MeasurementSet | raw data from one assay on some samples | experimental data |
| AuxiliarySet | companion sequencing (hashing, barcodes, guides, library quantification) | cell hashing, quantification DNA barcode sequencing, ... |
| AnalysisSet | processed results | intermediate analysis, principal analysis |
| PseudobulkSet | per-cell-type aggregates from a single-cell analysis | |
| ConstructLibrarySet | the designed library a screen or reporter assay used | guide library, reporter library, editing template library, expression vector library, overexpression vector library |
| CuratedSet | reference material: variant lists, elements, guides, external data | variants, elements, guide RNAs, editing templates, barcodes, training data for predictive models, external sequencing data, external data for catalog, QTL, genome, transcriptome, ... |
| ModelSet | a trained model | neural network, variant binding effect, logistic regression, random forest, decision tree, support vector machine |
| PredictionSet | a model's predictions | element-gene links, non-coding variant effects, coding variant effects, coding variant effects on PPI, coding variants pathogenicity, variant TF binding effects, TF binding, disease associations, genetic constraint, chromatin accessibility, gene regulatory networks, spatial gene expression variability, meta analysis |
| Sample (InVitroSystem, Tissue, PrimaryCell, WholeOrganism, MultiplexedSample, TechnicalSample) | biological material, or a *virtual* sample that gives a model its context | |
| File (SequenceFile, AlignmentFile, MatrixFile, TabularFile, SignalFile, ModelFile, ConfigurationFile, ReferenceFile) | one file, in exactly one file set (`file_set`) | |
| AnalysisStepVersion → AnalysisStep → Workflow; SoftwareVersion → Software | how a file was made | |
| Document | protocols, plate maps, seqspec, **file format specifications** | |

## The links

Direction is written as the Portal stores it. `(reverse)` marks a calculated
field the Portal fills from the other side.

### Provenance: data to results to predictions

| Field | From → to | Meaning |
|---|---|---|
| `input_file_sets` / `input_for` (reverse) | AnalysisSet, PseudobulkSet, ModelSet, PredictionSet → any FileSet | what a result was computed from. Chains: MeasurementSet → intermediate AnalysisSet → principal AnalysisSet → PseudobulkSet → PredictionSet |
| `derived_from` / `input_file_for` (reverse) | File → File | file-level provenance (matrix from reads, predictions from a model file and the input data) |
| `analysis_step_version` | File → AnalysisStepVersion | the step (type, input and output content types, parent steps) and software versions that made the file. Model files reach their software through the ModelSet |
| `workflows` | AnalysisSet, File → Workflow | the pipeline, and whether it is the IGVF uniform pipeline |
| `quality_metrics` | File → *QualityMetric | the pipeline's own QC numbers |
| `superseded_by` / `supersedes` | FileSet, File, Sample → newer object | a replacement exists; prefer it (316 file sets were superseded on 2026-09-23) |
| `control_file_sets` / `control_for` (reverse) | FileSet → control FileSet | negative controls (e.g. non-targeting guides) |

### Design: what was built and tested

| Field | From → to | Meaning |
|---|---|---|
| `construct_library_sets` | MeasurementSet, AnalysisSet, Sample → ConstructLibrarySet | the guide, reporter or editing-template library |
| `integrated_content_files` | ConstructLibrarySet → TabularFile / ReferenceFile | the designed content: guide RNA sequences, MPRA sequence designs, editing templates, variants |
| `small_scale_gene_list`, `large_scale_gene_list`, `small_scale_loci_list`, `large_scale_loci_list`, `targeton`, `exon`, `scope`, `guide_type`, `selection_criteria` | ConstructLibrarySet | what the library targets |
| `enrichment_designs`, `barcode_replacement_file` | MeasurementSet → TabularFile | capture / enrichment design, barcode corrections |
| `auxiliary_sets`, `related_measurement_sets` (multiome partner) | MeasurementSet → FileSet | companion data |
| `seqspecs`, `seqspec_of` (reverse), `seqspec_document` | SequenceFile ↔ ConfigurationFile, Document | read structure |
| `onlist_files`, `barcode_map` (MultiplexedSample), `hashtag_barcode_map` (AuxiliarySet) | → TabularFile | barcode onlists and barcode-to-sample maps |

### Samples: the sample tree

| Field | Meaning |
|---|---|
| `sorted_from` + `sorted_from_detail` / `sorted_fractions` (reverse) | FACS / bin sorting: each bin is its own sample, with its own MeasurementSet (HT-Recruit ON/OFF, VAMP-seq bins, Starita low-to-high, base-editing top/bottom 20%) |
| `originated_from` / `origin_of` (reverse) | derivation, e.g. SGE time points from the day-0 population |
| `part_of` / `parts` (reverse) | replicates and sub-samples |
| `pooled_from` / `pooled_in`, `demultiplexed_from` / `demultiplexed_to`, `multiplexed_samples` / `multiplexed_in` | pooling and genetic demultiplexing |
| `time_post_library_delivery` (+ units), `time_post_change` | time points |
| `treatments`, `modifications` (CRISPR modality, e.g. CRISPRi dCas9-KRAB, base editing SpRY-ABE8e, Cas9), `construct_library_sets` | what was done to the cells |
| `virtual` | a stand-in sample for a model or prediction (GTEx tissues, cell types); context, not measured material |

### What a result is about

| Field | On | Meaning |
|---|---|---|
| `associated_phenotypes` | PredictionSet, ConstructLibrarySet | traits and diseases (one prediction set per trait for GWAS-based methods) |
| `assessed_genes` | PredictionSet, ModelSet | genes a prediction or model is about (e.g. the TF of a binding model) |
| `cell_type` | PredictionSet, PseudobulkSet | the cell type of the prediction |
| `targeted_genes`, `functional_assay_mechanisms`, `crispr_screen_biometric` | MeasurementSet, AnalysisSet | the screen's targets and readout (e.g. LDL-C uptake) |
| `external_input_data` | ModelSet → TabularFile | training data held outside IGVF |
| `model_zoo_location`, `model_name`, `model_version`, `prediction_objects` | ModelSet | where the model lives and what it predicts |
| `dbxrefs`, `url` | CuratedSet | external records (e.g. ENCODE `ENCSR` accessions) |
| `file_format_specifications` | File → Document | **the column definitions** of a tabular output |
| `externally_hosted`, `external_host_url` | ModelFile | the model file is hosted elsewhere |
| `publications` | FileSet | citable papers |

## Submission patterns in the diagrams

Each pattern is a diagram in `Data/IGVF/DataModels/`; the examples are live
Portal objects that follow it.

1. **Workflow template** (Workflow example). MeasurementSet reads → `derived_from` MatrixFile in an intermediate AnalysisSet (quantification, filtering) → principal AnalysisSet (merging). Each file carries `analysis_step_version` → AnalysisStep (types quantification / filtering / merging, `parents` for order) → Workflow, with SoftwareVersion → Software (kallisto, scanpy, anndata).
2. **Predictive-model template** (Predictive Model example, Models temp). ModelSet `input_file_sets` → the experimental AnalysisSet it was trained on; model file `derived_from` the training table; PredictionSet `input_file_sets` → ModelSet (and the data it was run on); prediction table `derived_from` the model file; loci lists as `large_scale_loci_list` in a CuratedSet; virtual samples and donors.
3. **TF binding models** (Boyle SEMpl, Allen QBiC-SELEX, Models temp): about 300 ModelSets (`variant binding effect`), each with `assessed_genes` = the TF; training data in a CuratedSet (`training data for predictive models`); PredictionSets of type `variant TF binding effects`. Example: `processed discover --gene GATA1` finds SEMVAR and QBiC-SELEX sets.
4. **Neural-network models on reporter or chromatin data** (Allen BlueSTARR, PRINT / seq2PRINT, Yi eSIG-Net): ModelSets of type `neural network` trained on STARR-seq, SHARE-seq or protein features; virtual samples per cell line and treatment (K562, A549 with DMSO or Dex); chains AnalysisSet (fragments) → AnalysisSet (DNA footprints) → ModelSet → AnalysisSet (sequence attributions, TF binding scores).
5. **Coding-variant predictors and calibration** (Craven ESM-1v, MutPred2, Sunyaev, calibration): PredictionSets of coding variant effects or pathogenicity; ModelFile possibly `externally_hosted`; calibration is an AnalysisStep of type calibration whose output is `derived_from` the raw scores, with PredictionSets taking other PredictionSets and CuratedSets (ClinVar-labelled, SGE / VAMP-seq) as input.
6. **Trait-level predictions** (Boyle T-Land, Dey ColocBoost / cV2F, Liu kidney scorecard, Price, Amariuta and Sunyaev lab models): one PredictionSet per trait with `associated_phenotypes`, input CuratedSets of GWAS or Y2AVE variants (`variants`), virtual samples (GTEx tissues), a file format specification document. Example: IGVFDS0029SZLZ (ColocBoost, coronary artery disease; its tables point to the "Colocboost predictions file format" document).
7. **Sorting screens** (HT-Recruit, Starita, VAMP-seq MultiSTEP; live example IGVFDS6464SOVZ): a ConstructLibrarySet (guide or expression-vector library, `integrated_content_files` = guides / variants) → cells → `sorted_from` fractions → one MeasurementSet per fraction; an AuxiliarySet for barcode-to-variant sequencing feeds an intermediate AnalysisSet; the principal AnalysisSet takes both.
8. **Saturation genome editing** (SGE model, SGE temp, Starita): editing-template library per targeton (`targeton`, `small_scale_gene_list`, `integrated_content_file` in a CuratedSet of editing templates); HAP-1 samples `originated_from` the day-0 population with `time_post_library_delivery` 5 / 13 / 17 days and a Cas9 `modification`; one MeasurementSet per time point plus a negative control; an AnalysisSet over all of them.
9. **MPRA and external sequencing data** (Boyle ENCODE MPRA, Li CATlas, DisCO-VG): a CuratedSet of `external sequencing data` with `dbxrefs` (ENCSR) and URLs, and a CuratedSet of `elements` holding the MPRA sequence designs, both as `input_file_sets` of the principal AnalysisSet; for CATlas, per-cell-type PseudobulkSets (cell annotation, fragments, bigWig, h5ad) built from the external CuratedSet, each with virtual samples.
10. **Perturb-seq** (Yi): guide-library ConstructLibrarySet with a seqspec ConfigurationFile, guide tables as `integrated_content_files`, a CRISPRi `modification`, a plasmid-map Document.
11. **eQTL study** (Weng / Lin / Garber): primary cells and donors (with institutional certificates for controlled access) → MeasurementSets (RNA-seq, WGS) → intermediate AnalysisSets (gene count matrix, genotypes) → principal AnalysisSet (eQTL results).

## Querying well

- **From an accession:** `processed lineage` first. The Portal usually already holds the processed files, their QC, the design and the predictions built on them.
- **From a topic:** `processed discover`. Phenotypes live on PredictionSets (`associated_phenotypes`), not MeasurementSets. Tissues are sample terms (`samples.sample_terms.term_name`); for predictions those samples are virtual. Genes are `assessed_genes` on predictions and models, `targeted_genes` on screens.
- **Exact values:** search filters need the Portal's exact term, and an empty search returns HTTP 404. Read the facet first (`portal facets --field`), or let `discover` match the words. Facets are capped at about 100 terms, so rare values must be counted from the items.
- **Mixed-type searches** (AnalysisSet + PredictionSet) return only the facets the types share; count type-specific fields from the items.
- **Superseded data:** check `superseded_by` and prefer the replacement.
- **Tabular outputs:** read the `file_format_specifications` document before parsing columns.

## Where the drafts differ from the live schema

- PredictionSet `file_set_type` values in the diagrams ("functional effect", "binding effect", "pathogenicity") are drafts; the live enum is the list in the table above (e.g. "non-coding variant effects", "variant TF binding effects", "coding variants pathogenicity").
- ModelSet has no `software_version`; its software is recorded through its files' analysis step versions.
- MeasurementSet has no `control_type` (controls are `control_file_sets` and sample-level); `externally_hosted` is on ModelFile, not TabularFile.
- Field-name misspellings in the drafts (`sotfware_version`, `seqspec_dopcument`, `iinput_file_sets`, `file_format_specs`) are not Portal fields.
- A `workflow` link on AnalysisStep, drawn in several diagrams, is not in the schema; Workflow lists its `analysis_step_versions`.
