# Skills

IGVFagent has 117 skills exposing 599 typed tools. Every skill is an `igvfagent <skill> <subcommand>` command, and the agent calls the same commands as tools.

**On this page**

- [Capabilities](#capabilities)
- [Skill pages](#skill-pages)
- [All skills](#all-skills)

## Capabilities

- IGVF Portal, IGVF Catalog API, IGVF Knowledge Graph (ArangoDB) access.
- ENCODE Portal metadata search and file download.
- Variant annotation against IGVF Catalog evidence (CADD, QTL, phenotypes,
  regulatory elements, predictions).
- **Advanced variant analysis**: integrates Catalog + ENCODE cCRE evidence
  with optional user experimental tables (CRISPRi / MPRA / GWAS), fits
  logistic models, and produces volcano / Miami / evidence-overlap plots
  plus a research-grade markdown report.
- Single-cell RNA-seq, single-cell ATAC-seq, Perturb-seq, and 10x Multiome
  metadata-first workflows.
- Enhancer–gene linkage retrieval and comparison (ABC, rE2G, ENCODE-rE2G,
  catalog-based predictions, eQTL-based linkage).
- MPRA / STARR / BlueSTARR retrieval, summary statistics, and plotting.
- CRISPRi / CRISPR-FACS / Perturb-seq evidence integration with functional
  annotation.
- cCRE (SCREEN) discovery and FAVOR-based variant annotation, plus IGV-like
  browser views.
- Data illustration and interpretation across IGVF and ENCODE search URLs.
- **ChIP-Atlas (Ohta/Oki) reprocessed peak archive** — browse/search/download
  ChIP-seq / ATAC-seq / DNase-seq / Bisulfite peaks across 10 genomes, fetch
  pre-computed Target-Genes tables, and submit TF-enrichment jobs.
- **MaveDB → genomic-coordinate mapping**, including a dedicated **SGE
  (Saturation Genome Editing) cDNA-coordinate path** that handles the full
  HGVS-c grammar (CDS / intronic / 5′UTR / 3′UTR) used by SGE scoresets.
- **Synapse / Sage Bionetworks retrieval** — anonymous metadata walk + search
  and PAT-authenticated download for controlled-access deposits (PsychENCODE,
  AMP-AD/PD, ROSMAP) that IGVF distributes off-Portal.
- **Functional-assay calibration to ACMG/AMP evidence** — bootstrap constrained
  skew-normal mixture fitting + Bayesian (Tavtigian) calibration turns MAVE /
  VAMP-seq / SGE / cell-fitness scores into **PS3 / BS3 evidence strengths**
  (supporting → very strong) with the score window for each, so an assay score
  can be used directly in clinical variant classification.
- **Skill-correctness benchmark suite** — more than twenty-five recent Nature / Cell /
  Science / Nat Genet / Nat Methods / Genome Biol papers reproduced directly
  from public data, each scored against machine-readable ground-truth checks.
  **Scope:** every `run.sh` invokes skills directly as shell commands, with no
  model in the loop, so the suite establishes that the *implementations* are
  correct. It does not evaluate planning, tool selection, or the validity of an
  end-to-end conclusion — see [`Docs/EVALUATION.md`](../../EVALUATION.md). Includes
  **five full single-cell / multiome local reproductions** that download the
  public data and run IGVFagent's real analytical chain end-to-end — Travaglini
  2020 lung (`sc-analyze`), Trevino 2021 cortex multiome (`multiome peak2gene`),
  deMULTIplex2 / Stoeckius cell hashing (`multiseq`), Rosenberg 2018
  SPLiT-seq CNS (`splitseq`), and Ma 2020 SHARE-seq skin (`share`) — with a
  2025 developing-neocortex multiome atlas (Wang, CELLxGENE) in progress. The
  single-cell loaders + QC these exercise are internalized into
  `Scripts/_scload.py` and the skills (see [`Benchmarks/`](../../../Benchmarks/README.md)).
- **Complete MPRA library toolchain** — three clean-room ports of the
  kircherlab suite (Max Schubach / BIH, MIT), covering the whole assay
  lifecycle: **`oligo`** designs the library (region tiling, REF/ALT variant
  oligos, synthesis filters, adapters), **`mpraflow`** turns sequencing into
  per-oligo activity (barcode→oligo assignment, count normalisation, outlier
  removal, replicate aggregation, allelic skew, Lincoln-Petersen complexity),
  and **`mpralib`** does barcode QC (three outlier detectors, replicate
  agreement) plus validation against the **eight IGVF MPRA community file
  standards**. All standard-library-only — no Snakemake, conda, R, bedtools or
  pandas — and verified bit-identical against the upstream scripts.
- **Enhancer / regulatory-REGION annotation** (`enhancer-annot`) — upload an
  enhancer list (xlsx / BED / TSV / CSV) and annotate every **interval**
  against the ENCODE SCREEN cCRE registry: how many cCREs overlap, a ranked
  representative class, bases covered as a union, coverage fraction, and — for
  regions with no overlap — the nearest element and its distance. Distinct from
  the point-level `ccre annotate-variants`. The registry (V4 default,
  2,348,854 elements; V3 selectable) is indexed locally so annotation runs
  offline.
- **Expandable by design — bring your own skills and tools.** A user-extension
  framework auto-discovers your custom tools (one YAML manifest wrapping any
  executable — no code) and custom skills (one Python file → a first-class
  `igvfagent <name>` subcommand) from `~/.igvfagent/` or `UserExtensions/`,
  and absorbs them into the CLI, the `ask` LLM agent, and the UI with zero
  core-code edits. See [Extending IGVFagent](../extending.md).

Direct `python3 Scripts/<module>.py` invocations are listed in [`Scripts/README.md`](../../../Scripts/README.md).

## Skill pages

| Page | Covers |
|---|---|
| [Finding and retrieving data](portal-and-data-access.md) | IGVF Portal, Catalog and ENCODE access, the processed-first Portal lineage, dataset explanation, literature, GEO, Synapse, ChIP-Atlas, the Perturbation Catalogue and local reference databases. |
| [Variants and variant effects](variants.md) | Variant annotation, advanced variant analysis, cCRE and FAVOR, regulatory-region annotation, MaveDB, saturation genome editing and assay calibration to ACMG/AMP evidence. |
| [Single-cell, multiome and spatial](single-cell.md) | scRNA-seq, multiome, SPLiT-seq, SHARE-seq, snMCT-seq, cell hashing, Tabula Sapiens, the IGVF single-cell pipeline, principal pseudobulks and Spatial-ATAC-Hi-C. |
| [CRISPR screens and Perturb-seq](crispr-and-perturbation.md) | Pooled screens from raw reads, CRISPRi, Flow-FISH, Perturb-seq pipelines (IGVF-CRISPR and pinellolab), SCEPTRE, TF Perturb-seq and the CRISPR jamboree analyses. |
| [Enhancer–gene linkage (E2G)](enhancer-gene.md) | Enhancer–gene predictions and their benchmarks: ABC, ENCODE-rE2G, scE2G, eQTL and GWAS enrichment, super-enhancer targets, and E2G QC and Portal submission. |
| [MPRA and STARR-seq](mpra-starr.md) | Massively parallel reporter assays from library design to allelic effects, STARR-seq and scQers. |
| [Bulk genomics, RNA-seq and proteomics](bulk-genomics-proteomics.md) | The ENCODE bulk-genomics pipeline, bulk RNA-seq, and the protein-interaction knowledge graph. |
| [Knowledge graph, warehouse and networks](knowledge-graph.md) | The local knowledge graph that grows with use, the IGVF Catalog mirror, graph traversal, the DuckDB warehouse and network integration. |

## All skills

All 117 `igvfagent` skills, generated from the `SKILLS` registry in [`Scripts/cli.py`](../../../Scripts/cli.py). Run `igvfagent <skill> --help` for its subcommands and options, or `igvfagent tools` for the typed tools the agent calls.

### Finding and retrieving data

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `catalog` | IGVF Catalog (KG) canonical-query layer — universal get-entity, search-region, find-associations (semantic), find-ld, resolve-id, list-sources (clean-room reimpl of IGVF-DACC/igvf-catalog-mcp) | — | [`catalog_query_skill.py`](../../../Scripts/catalog_query_skill.py) |
| `chipatlas` | ChIP-Atlas (Ohta/Oki) reprocessed ChIP/ATAC/DNase/Bisulfite peak archive — browse + per-experiment files + assembled BEDs + Target-Genes + WABI Enrichment (clean-room reimpl of inutano/chip-atlas) | [ChIP-Atlas reprocessed peak archive](portal-and-data-access.md#chip-atlas-reprocessed-peak-archive) | [`chipatlas_skill.py`](../../../Scripts/chipatlas_skill.py) |
| `client` | IGVF Portal / Catalog / Knowledge Graph / ENCODE client | [IGVF / ENCODE / Knowledge Graph client](portal-and-data-access.md#igvf--encode--knowledge-graph-client) | [`igvf_client.py`](../../../Scripts/igvf_client.py) |
| `data` | Catalog / Portal / ENCODE overview + smoke summaries | [Catalog / Portal / ENCODE overviews](portal-and-data-access.md#catalog--portal--encode-overviews) | [`igvf_data_skills.py`](../../../Scripts/igvf_data_skills.py) |
| `document` | Read an uploaded manuscript (PDF/DOCX/text) and derive accessions, assays and a reproduction plan | — | [`document_ingest_skill.py`](../../../Scripts/document_ingest_skill.py) |
| `enrich` | GO + Pathway enrichment validation (ORA via Enrichr; GSEA preranked) over GO_BP/MF/CC + Reactome + KEGG + WikiPathways + MSigDB Hallmark | — | [`enrichment_skill.py`](../../../Scripts/enrichment_skill.py) |
| `explain` | Explain an IGVF / ENCODE accession or search URL | [Data illustration and interpretation](portal-and-data-access.md#data-illustration-and-interpretation) | [`data_illustration_interpretation.py`](../../../Scripts/data_illustration_interpretation.py) |
| `figshare` | figshare retrieval — article / files / download / search. Resolves numeric id, DOI, article URL, or private /s/<token> share link (as printed in paper Data Availability statements). Downloads with md5 verify. Pure urllib + json; the general-purpose counterpart to Zenodo for research-data deposits. | — | [`figshare_skill.py`](../../../Scripts/figshare_skill.py) |
| `frontpage` | Refresh front-page Portal + Knowledge Graph stats | [Catalog / Portal / ENCODE overviews](portal-and-data-access.md#catalog--portal--encode-overviews) | [`igvf_frontpage_summary.py`](../../../Scripts/igvf_frontpage_summary.py) |
| `geo` | NCBI GEO retrieval (search / metadata / download) | [GEO retrieval](portal-and-data-access.md#geo-retrieval) | [`geo_retrieval_skill.py`](../../../Scripts/geo_retrieval_skill.py) |
| `humantfs` | Human Transcription Factors database (Lambert 2018, humantfs.ccbr.utoronto.ca) — local SQLite mirror of all 1,639 curated TFs with DBD family, binding mode, motif status and cross-references; is-tf / lookup / list / motifs / families / export | [Human Transcription Factors database](portal-and-data-access.md#human-transcription-factors-database) | [`humantfs_skill.py`](../../../Scripts/humantfs_skill.py) |
| `intpath` | Integrated pathway data (KEGG + WikiPathways + BioCyc) via IntPath | — | [`intpath_skill.py`](../../../Scripts/intpath_skill.py) |
| `pathway-viz` | Pathway / PPI network figure for a gene list (STRING + Reactome + local KG) | — | [`pathway_viz_skill.py`](../../../Scripts/pathway_viz_skill.py) |
| `pathwaydb` | Pull CURRENT KEGG / Reactome / WikiPathways releases, normalise and integrate them locally | [Pathway databases (KEGG + Reactome + WikiPathways, integrated locally)](portal-and-data-access.md#pathway-databases-kegg--reactome--wikipathways-integrated-locally) | [`pathwaydb_skill.py`](../../../Scripts/pathwaydb_skill.py) |
| `perturb-catalog` | Perturbation Catalogue retrieval: MAVE / CRISPR-screen / Perturb-seq datasets, per-row effects, GSEA, downloads | [Perturbation Catalogue retrieval](portal-and-data-access.md#perturbation-catalogue-retrieval) | [`perturbation_catalog_skill.py`](../../../Scripts/perturbation_catalog_skill.py) |
| `portal` | IGVF Portal canonical-query layer — faceted search, report.tsv, batch-download, endpoint-param introspection (clean-room reimpl of IGVF-DACC/igvf-portal-mcp) | — | [`portal_query_skill.py`](../../../Scripts/portal_query_skill.py) |
| `portal-qc` | Portal-wide QC: the Portal's own audit facets by lab and type, plus provenance integrity it does not surface | — | [`igvf_portal_qc_skill.py`](../../../Scripts/igvf_portal_qc_skill.py) |
| `processed` | What IGVF has ALREADY computed from an accession — uniform-pipeline matrices, fragments, peaks — so a workflow downloads the result instead of recomputing it from raw reads | [Start from what the IGVF Portal already computed (`processed lineage`)](portal-and-data-access.md#start-from-what-the-igvf-portal-already-computed-processed-lineage) | [`processed_first_skill.py`](../../../Scripts/processed_first_skill.py) |
| `ref` | Literature retrieval / validation / study design | [Reference skill (literature retrieval, validation, design)](portal-and-data-access.md#reference-skill-literature-retrieval-validation-design) | [`reference_skill.py`](../../../Scripts/reference_skill.py) |
| `synapse` | Sage Bionetworks Synapse retrieval — entity / children / walk / search / download. Anonymous for public deposits; PAT-authenticated for PsychENCODE, AMP-AD, AMP-PD, and IGVF-controlled Synapse cohorts (set SYNAPSE_AUTH_TOKEN). Pure urllib + json, no synapseclient dep. | [Synapse / Sage Bionetworks retrieval](portal-and-data-access.md#synapse--sage-bionetworks-retrieval) | [`synapse_skill.py`](../../../Scripts/synapse_skill.py) |

### Variants and variant effects

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `alphagenome` | AlphaGenome (google-deepmind/alphagenome): predictions, variant scores, ISM and Atlas scores, with IGVF Catalog rsID/gene resolution | [AlphaGenome predictions and Atlas scores (`alphagenome`)](variants.md#alphagenome-predictions-and-atlas-scores-alphagenome) | [`alphagenome_skill.py`](../../../Scripts/alphagenome_skill.py) |
| `advanced-variant` | Integrated variant scoring + logistic + report | [Advanced variant analysis](variants.md#advanced-variant-analysis) | [`advanced_variant_analysis.py`](../../../Scripts/advanced_variant_analysis.py) |
| `calibrate` | Functional-assay calibration to ACMG/AMP evidence — clean-room reimpl of rosstewart/exCALIBR (Zeiberg et al. bioRxiv 2025). Bootstrap constrained skew-normal mixture EM + Bayesian (Tavtigian) calibration turns MAVE / VAMP-seq / SGE scores into PS3 / BS3 evidence strengths. Subcommands: thresholds, prepare, run, assign, selftest. | [Assay calibration → ACMG/AMP evidence (exCALIBR)](variants.md#assay-calibration--acmgamp-evidence-excalibr) | [`excalibr_skill.py`](../../../Scripts/excalibr_skill.py) |
| `ccre` | cCRE / FAVOR / linkage annotations (variant-level) | [cCRE, FAVOR, IGV-style browser views](variants.md#ccre-favor-igv-style-browser-views) | [`ccre_linkage_annotation_skills.py`](../../../Scripts/ccre_linkage_annotation_skills.py) |
| `enhancer-annot` | Enhancer / regulatory-REGION annotation against the ENCODE SCREEN cCRE registry: upload an enhancer list (xlsx/BED/TSV/CSV), get per-region cCRE overlap, class breakdown and nearest-element distance | [Enhancer / regulatory-region cCRE annotation](variants.md#enhancer--regulatory-region-ccre-annotation) | [`enhancer_annotation_skill.py`](../../../Scripts/enhancer_annotation_skill.py) |
| `mavedb` | MaveDB scoreset → genomic coords (chr/pos/ref/alt) | [MaveDB mapping (incl. SGE cDNA path)](variants.md#mavedb-mapping-incl-sge-cdna-path) | [`mavedb_mapping_skill.py`](../../../Scripts/mavedb_mapping_skill.py) |
| `regulome` | RegulomeDB regulatory rank + per-tissue scores for non-coding variants | — | [`regulome_skill.py`](../../../Scripts/regulome_skill.py) |
| `sge` | Saturation genome editing: variant counts and log2 functional scores from amplicon reads | [Saturation genome editing (SGE)](variants.md#saturation-genome-editing-sge) | [`sge_analysis.py`](../../../Scripts/sge_analysis.py) |
| `variant` | Annotate variants against IGVF Catalog evidence | [Variant annotation](variants.md#variant-annotation) | [`annotate_variant_list.py`](../../../Scripts/annotate_variant_list.py) |
| `variant-list` | Annotate a pasted/file variant list in any notation (chr-pos-ref-alt, rsID, SPDI, HGVS, VCF) via FAVOR + IGVF Catalog, and add the results to the local KG | — | [`variant_list_skill.py`](../../../Scripts/variant_list_skill.py) |
| `variant-verify` | Cross-source verification of variant annotations (FAVOR snapshot vs live ClinVar) | — | [`variant_verify_skill.py`](../../../Scripts/variant_verify_skill.py) |

### Single-cell, multiome and spatial

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `igvf-sc-pipeline` | Port of IGVF/single-cell-pipeline: SHARE-seq/10x multiome barcode correction, trimming, kb/chromap wrappers, TSS enrichment, RNA/ATAC/joint QC, barcode-rank knees, Synapse/portal I/O. | [IGVF single-cell pipeline (`igvf-sc-pipeline`)](single-cell.md#igvf-single-cell-pipeline-igvf-sc-pipeline) | [`igvf_sc_pipeline_skill.py`](../../../Scripts/igvf_sc_pipeline_skill.py) |
| `matrix-qc` | Whole-matrix h5ad summaries that state what they measured (full vs sampled), and threshold aggregation over a table with no warehouse needed | — | [`matrix_qc_skill.py`](../../../Scripts/matrix_qc_skill.py) |
| `mct` | snMCT-seq: RNA and methylation from the same nuclei, analysed and cross-compared | [snMCT-seq — RNA and methylation from the same nuclei](single-cell.md#snmct-seq--rna-and-methylation-from-the-same-nuclei) | [`mct_analysis.py`](../../../Scripts/mct_analysis.py) |
| `multiome` | 10x Multiome retrieval pipeline | [Single-cell, multiome, specialized assays](single-cell.md#single-cell-multiome-specialized-assays) | [`multiome_10x_pipeline.py`](../../../Scripts/multiome_10x_pipeline.py) |
| `multiseq` | MULTI-seq / Cell Hashing demultiplexing (Python port of deMULTIplex2) | [MULTI-seq / Cell Hashing demultiplexing](single-cell.md#multi-seq--cell-hashing-demultiplexing) | [`multiseq_analysis_skill.py`](../../../Scripts/multiseq_analysis_skill.py) |
| `principal-pseudobulks` | Port of EngreitzLab/generate-principal-pseudobulks: IGVF multiome QC guide, filtered fragments, gene-symbol RNA matrix and scE2G config per cluster. | [Principal pseudobulks (`principal-pseudobulks`)](single-cell.md#principal-pseudobulks-principal-pseudobulks) | [`principal_pseudobulks_skill.py`](../../../Scripts/principal_pseudobulks_skill.py) |
| `raw-pipeline` | Process a dataset's raw reads: FASTQ -> count matrix (kallisto\|bustools) -> single-cell analysis | — | [`raw_data_pipeline.py`](../../../Scripts/raw_data_pipeline.py) |
| `sc-analyze` | Single-cell analysis: QC, PCA, UMAP, t-SNE, Leiden, markers, publication figures (Scanpy-driven) | [Single-cell analysis (UMAP / t-SNE / Leiden / markers)](single-cell.md#single-cell-analysis-umap--t-sne--leiden--markers) | [`singlecell_analysis.py`](../../../Scripts/singlecell_analysis.py) |
| `sceps` | scEPS single-cell disease-neighborhood statistics — clean-room reimpl of Genentech/sceps. Integrates GWAS (MAGMA Z) + single-cell atlas; per-neighborhood variance-component d-statistic (GWAS vs matched control genes). Subcommands: estimate. Runs under the scEPS env (scanpy/anndata/statsmodels). | — | [`sceps_skill.py`](../../../Scripts/sceps_skill.py) |
| `scnt-seq` | scNT-seq: split new from old RNA by 4sU T>C metabolic-labelling conversions | — | [`scnt_seq_skill.py`](../../../Scripts/scnt_seq_skill.py) |
| `sctransform` | SCTransform variance-stabilising normalisation (regularised NB regression, Pearson residuals) | — | [`sctransform_skill.py`](../../../Scripts/sctransform_skill.py) |
| `share` | SHARE-seq joint scATAC+scRNA QC (Ma 2020 / Broad pipeline) | [SHARE-seq joint scATAC + scRNA QC](single-cell.md#share-seq-joint-scatac--scrna-qc) | [`share_seq_skill.py`](../../../Scripts/share_seq_skill.py) |
| `singlecell` | Single-cell discovery and example analysis | [Single-cell, multiome, specialized assays](single-cell.md#single-cell-multiome-specialized-assays) | [`single_cell_data_skills.py`](../../../Scripts/single_cell_data_skills.py) |
| `spatial-hic` | Spatial-ATAC-Hi-C — spatially resolved 3D genome + chromatin accessibility on tissue slides (Wang 2026, Nat Methods; GSE307620). 50x50 pixel demux, per-pixel contact QC + TSS enrichment, GAS/GAD gene scores, scHiCluster-style imputation, A/B compartments, per-pixel CNV, loop quantification + APA, and spatial rendering (clean-room reimpl of wangjuan001/Spatial-ATAC-Hi-C) | [Spatial-ATAC-Hi-C (spatial 3D genome + chromatin accessibility)](single-cell.md#spatial-atac-hi-c-spatial-3d-genome--chromatin-accessibility) | [`spatial_atac_hic_skill.py`](../../../Scripts/spatial_atac_hic_skill.py) |
| `specialized` | Specialized IGVF assay catalog (Parse SPLiT-seq, etc.) | [Single-cell, multiome, specialized assays](single-cell.md#single-cell-multiome-specialized-assays) | [`igvf_specialized_data_skills.py`](../../../Scripts/igvf_specialized_data_skills.py) |
| `splitseq` | Parse SPLiT-seq end-to-end pipeline | [Parse SPLiT-seq pipeline](single-cell.md#parse-split-seq-pipeline) | [`splitseq_pipeline.py`](../../../Scripts/splitseq_pipeline.py) |
| `tabula` | Tabula Sapiens 2.0 human cell atlas (Cell 2026) — retrieval from figshare / GEO GSE306755 / CELLxGENE, raw-bucket inventory, and figure reproduction: dataset overview, TF cell-type specificity (tau) | [Tabula Sapiens 2.0 (reference human cell atlas)](single-cell.md#tabula-sapiens-20-reference-human-cell-atlas) | [`tabula_sapiens_skill.py`](../../../Scripts/tabula_sapiens_skill.py) |

### CRISPR screens and Perturb-seq

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `bean` | Base-editing screens: BEAN-style base-edit-aware guide assignment and editing-activity reporting | [CRISPR screens from raw reads — three shapes, three tools](crispr-and-perturbation.md#crispr-screens-from-raw-reads--three-shapes-three-tools) | [`base_editing_screen.py`](../../../Scripts/base_editing_screen.py) |
| `bean-benchmark` | Reproduce the crispr-bean paper (Ryu et al. 2024) on its own deposited data, claim by claim | — | [`bean_benchmark_skill.py`](../../../Scripts/bean_benchmark_skill.py) |
| `bulk-crispr` | Port of IGVF-CRISPR/bulk_crispr_pipeline: kite guide counting, per-lane cell QC + doublets, MULTI-seq, cis gene sets, SCEPTRE-style tests, BH/Fisher results, tracks. | [Perturb-seq pipeline, bulk_crispr_pipeline port (`bulk-crispr`)](crispr-and-perturbation.md#perturb-seq-pipeline-bulk_crispr_pipeline-port-bulk-crispr) | [`bulk_crispr_pipeline_skill.py`](../../../Scripts/bulk_crispr_pipeline_skill.py) |
| `crispr-fg-jamboree` | Port of IGVF-CRISPR/CRISPR_FG_JAMBOREE: assay specs, Cell Ranger inputs, pipeline config, MuData QC, guide calling, cis differential perturbation, genome tracks. | [CRISPR-FG scPerturb-seq jamboree (`crispr-fg-jamboree`)](crispr-and-perturbation.md#crispr-fg-scperturb-seq-jamboree-crispr-fg-jamboree) | [`crispr_fg_jamboree_skill.py`](../../../Scripts/crispr_fg_jamboree_skill.py) |
| `crispr-jamboree2` | Port of IGVF-CRISPR/CRISPR-JAMBOREE (2nd jamboree): guide counting, seqspec, guide assignment, 6 Perturb-seq inference modules, NB simulation, AUPRC/volcano/IGV/network evaluation. | [2nd CRISPR Jamboree analyses (`crispr-jamboree2`)](crispr-and-perturbation.md#2nd-crispr-jamboree-analyses-crispr-jamboree2) | [`crispr_jamboree2_skill.py`](../../../Scripts/crispr_jamboree2_skill.py) |
| `crispr-jamboree3` | Port of IGVF-CRISPR/CRISPR-jamboree3: the single-cell Perturb-seq pipeline stage by stage - config, seqspec, kb, QC, MuData, Scrublet, GMM demux, CLEANSER, sceptre, evaluation, dashboard. | [3rd CRISPR Jamboree single-cell pipeline (`crispr-jamboree3`)](crispr-and-perturbation.md#3rd-crispr-jamboree-single-cell-pipeline-crispr-jamboree3) | [`crispr_jamboree3_skill.py`](../../../Scripts/crispr_jamboree3_skill.py) |
| `crispr-jamboree4` | Port of IGVF-CRISPR/CRISPR-Jamboree_2025 (4th jamboree): checkpointed sceptre / PerTurbo inference, mergedResults outputs and the AUPRC/AUROC control-set evaluation. | [4th CRISPR Jamboree inference benchmark (`crispr-jamboree4`)](crispr-and-perturbation.md#4th-crispr-jamboree-inference-benchmark-crispr-jamboree4) | [`crispr_jamboree4_skill.py`](../../../Scripts/crispr_jamboree4_skill.py) |
| `crispr-pipeline` | Port of pinellolab/CRISPR_Pipeline: IGVF Perturb-seq pipeline - seqspec, mapping, QC, MuData, guide assignment, SCEPTRE/PerTurbo cis+trans inference, evaluation, dashboard. | [IGVF CRISPR Perturb-seq pipeline (`crispr-pipeline`)](crispr-and-perturbation.md#igvf-crispr-perturb-seq-pipeline-crispr-pipeline) | [`crispr_pipeline_skill.py`](../../../Scripts/crispr_pipeline_skill.py) |
| `crispr-screen` | CRISPR FACS screens: reprocess a whole screen (all sorted bins) into per-variant effects | [CRISPR screens from raw reads — three shapes, three tools](crispr-and-perturbation.md#crispr-screens-from-raw-reads--three-shapes-three-tools) | [`crispr_screen_analysis.py`](../../../Scripts/crispr_screen_analysis.py) |
| `crispr-seqspec` | Port of IGVF-CRISPR/CRISPR-SeqSpec: catalogue the CRISPR assay seqspecs, validate a seqspec (seqspec check), print kb/chromap/STARsolo read formats, list onlists, check FASTQs. | [CRISPR SeqSpec (`crispr-seqspec`)](crispr-and-perturbation.md#crispr-seqspec-crispr-seqspec) | [`crispr_seqspec_skill.py`](../../../Scripts/crispr_seqspec_skill.py) |
| `crispr-surf` | Deconvolve a CRISPR tiling screen into regulatory regions (CRISPR-SURF method) | — | [`crispr_surf_skill.py`](../../../Scripts/crispr_surf_skill.py) |
| `crisprdevtools` | Port of IGVF-CRISPR/crisprdevtools: scaffold IGVF CRISPR Nextflow modules (main.nf, processes, bin, conda_envs, test, input.config) and check module layouts. | [CRISPR module dev tools (`crisprdevtools`)](crispr-and-perturbation.md#crispr-module-dev-tools-crisprdevtools) | [`crisprdevtools_skill.py`](../../../Scripts/crisprdevtools_skill.py) |
| `crispresso` | CRISPResso2 — genome-editing outcomes from amplicon reads (indels, substitutions, base-editing) | — | [`crispresso_skill.py`](../../../Scripts/crispresso_skill.py) |
| `crispri` | CRISPRi / CRISPR-FACS / Perturb-seq evidence | [CRISPRi / CRISPR-FACS / Perturb-seq](crispr-and-perturbation.md#crispri--crispr-facs--perturb-seq) | [`crispri_data_skills.py`](../../../Scripts/crispri_data_skills.py) |
| `fishash-table2` | Port of IGVF-CRISPR/fishash-table2-reproduction: CLEANSER and SCEPTRE-mixture guide assignment on human/mouse barnyards, scored against Fishash Table 2. | [Fishash Table 2 reproduction (`fishash-table2`)](crispr-and-perturbation.md#fishash-table-2-reproduction-fishash-table2) | [`fishash_table2_skill.py`](../../../Scripts/fishash_table2_skill.py) |
| `flowfish` | CRISPRi Flow-FISH screen analysis (Fulco/Nasser/Engreitz) | [CRISPRi Flow-FISH screen](crispr-and-perturbation.md#crispri-flow-fish-screen) | [`flowfish_crispr_skill.py`](../../../Scripts/flowfish_crispr_skill.py) |
| `flowfish-pipeline` | Port of EngreitzLab/crispri-flowfish: CRISPRi-FlowFISH from FASTQs to guide counts, binned MLE effect sizes, qPCR scaling and enhancer calls. | [CRISPRi-FlowFISH pipeline (`flowfish-pipeline`)](crispr-and-perturbation.md#crispri-flowfish-pipeline-flowfish-pipeline) | [`flowfish_pipeline_skill.py`](../../../Scripts/flowfish_pipeline_skill.py) |
| `gasperini-pipeline` | Port of IGVF-CRISPR/Pipeline_Gasperini_2019: Gasperini 2019 pilot through the IGVF single-cell-like pipeline - QC, MuData, guide assignment, cis SCEPTRE CRT, original NB test, comparison. | [Gasperini 2019 pipeline (`gasperini-pipeline`)](crispr-and-perturbation.md#gasperini-2019-pipeline-gasperini-pipeline) | [`gasperini_pipeline_skill.py`](../../../Scripts/gasperini_pipeline_skill.py) |
| `gradient-screen` | FACS screens sorted into lettered expression bins (A-F): mean-bin scores per construct | [CRISPR screens from raw reads — three shapes, three tools](crispr-and-perturbation.md#crispr-screens-from-raw-reads--three-shapes-three-tools) | [`gradient_screen_analysis.py`](../../../Scripts/gradient_screen_analysis.py) |
| `guide-map` | Map FASTQ reads to a guide library with imperfect matching (sequencing error and/or base editing) | — | [`guide_map_skill.py`](../../../Scripts/guide_map_skill.py) |
| `igvf-crispr-pipeline` | Port of IGVF-CRISPR/IGVF_CRISPR_Pipeline: seqspec -> kb -x, guide features, kb ref/count, pure-Python guide counting, RNA AnnData QC and MuData assembly. | [IGVF CRISPR Pipeline, early Nextflow version (`igvf-crispr-pipeline`)](crispr-and-perturbation.md#igvf-crispr-pipeline-early-nextflow-version-igvf-crispr-pipeline) | [`igvf_crispr_pipeline_skill.py`](../../../Scripts/igvf_crispr_pipeline_skill.py) |
| `perturb-tools` | Port of IGVF-CRISPR/perturb-tools: pooled CRISPR screen object, RPM normalisation, replicate and sorting-screen LFC, QC, MAGeCK/Excel export, PoolQ reader, guide annotation and design. | [Pooled CRISPR screen toolkit (`perturb-tools`)](crispr-and-perturbation.md#pooled-crispr-screen-toolkit-perturb-tools) | [`perturb_tools_skill.py`](../../../Scripts/perturb_tools_skill.py) |
| `sc-crispr-de` | Single-cell CRISPR differential expression: per-guide negative-binomial GLM against no-guide cells, then alpha-RRA rank aggregation to gene level (method of Gersbach Lab sc-crispr-de) | [Single-cell CRISPR differential expression (`sc-crispr-de`)](crispr-and-perturbation.md#single-cell-crispr-differential-expression-sc-crispr-de) | [`sc_crispr_de_skill.py`](../../../Scripts/sc_crispr_de_skill.py) |
| `sceptre-igvf` | Port of IGVF-CRISPR/sceptreIGVF: SCEPTRE on IGVF CRISPR MuData - guide assignment, calibration, power and discovery tests, results written back to MuData uns. | [SCEPTRE for IGVF MuData (`sceptre-igvf`)](crispr-and-perturbation.md#sceptre-for-igvf-mudata-sceptre-igvf) | [`sceptre_igvf_skill.py`](../../../Scripts/sceptre_igvf_skill.py) |
| `tf-perturb` | IGVF TF Perturb-seq core (port of IGVF/tf_perturb_seq): Stage-3 QC, DEG calibration, pathways, cell/guide filters, energy-distance pipeline, cNMF export, cross-dataset summaries, WG3 disease & GWAS overlay | [TF Perturb-seq: calibrated effects to disease and GWAS (`tf-perturb`)](crispr-and-perturbation.md#tf-perturb-seq-calibrated-effects-to-disease-and-gwas-tf-perturb) | [`tf_perturb_seq_skill.py`](../../../Scripts/tf_perturb_seq_skill.py) |

### Enhancer–gene linkage (E2G)

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `abc` | Activity-by-Contact enhancer→gene prediction from candidate elements + ATAC/H3K27ac signal | [ABC pipeline (`abc-pipeline`)](enhancer-gene.md#abc-pipeline-abc-pipeline) | [`abc_skill.py`](../../../Scripts/abc_skill.py) |
| `abc-pipeline` | Port of broadinstitute/ABC-Enhancer-Gene-Prediction: full ABC pipeline - peaks, candidate regions, neighborhoods, qnorm, power-law/Hi-C predictions, thresholds, QC, Hi-C utilities. | [ABC pipeline (`abc-pipeline`)](enhancer-gene.md#abc-pipeline-abc-pipeline) | [`abc_pipeline_skill.py`](../../../Scripts/abc_pipeline_skill.py) |
| `e2g-qc-predictions` | Port of kaybrand/QC-and-Predictions: E2G cluster QC gate, filtering, scE2G output reformatting, QC stats, Cell Annotation cache and dry-run-first IGVF Portal submission. | [E2G QC, predictions and IGVF Portal submission (`e2g-qc-predictions`)](enhancer-gene.md#e2g-qc-predictions-and-igvf-portal-submission-e2g-qc-predictions) | [`e2g_qc_predictions_skill.py`](../../../Scripts/e2g_qc_predictions_skill.py) |
| `encode-re2g` | Port of EngreitzLab/ENCODE_rE2G: ENCODE-rE2G enhancer-gene features from ABC, 9 embedded pretrained models, thresholds, stats, CRISPR training and feature analysis. | [ENCODE-rE2G (`encode-re2g`)](enhancer-gene.md#encode-re2g-encode-re2g) | [`encode_re2g_skill.py`](../../../Scripts/encode_re2g_skill.py) |
| `enhancer` | Enhancer-gene linkage retrieval and comparison | [Enhancer–gene linkage](enhancer-gene.md#enhancergene-linkage) | [`enhancer_gene_linkage_skills.py`](../../../Scripts/enhancer_gene_linkage_skills.py) |
| `eqtl-enrich` | eQTL enrichment benchmark for E2G predictors (port of EngreitzLab/eQTLEnrichment): distal-noncoding fine-mapped eQTLs vs 1000G background, enrichment / recall (total, linking) across thresholds and distance bins, enrichment-recall curves, heatmaps | [eQTL enrichment benchmark for enhancer-gene predictors (`eqtl-enrich`)](enhancer-gene.md#eqtl-enrichment-benchmark-for-enhancer-gene-predictors-eqtl-enrich) | [`eqtl_enrichment_skill.py`](../../../Scripts/eqtl_enrichment_skill.py) |
| `gwas-e2g` | Port of EngreitzLab/GWAS_E2G_benchmarking: fine-mapped GWAS variant enrichment/recall in predicted enhancers vs 1000G SNPs, curves, gene-linking precision/recall vs silver standard + PoPS | [GWAS E2G benchmark (`gwas-e2g`)](enhancer-gene.md#gwas-e2g-benchmark-gwas-e2g) | [`gwas_e2g_benchmark_skill.py`](../../../Scripts/gwas_e2g_benchmark_skill.py) |
| `open4gene` | Open4Gene peak-to-gene linkage — clean-room reimpl of hbliu/Open4Gene (Liu et al. Science 2025). Hurdle model (logit zero + zero-truncated NB count) linking snATAC peaks to snRNA genes across multiome cells, per cell type, with covariates. Subcommands: link. statsmodels-based; validated vs the R pscl::hurdle reference. | — | [`open4gene_skill.py`](../../../Scripts/open4gene_skill.py) |
| `pgboost` | Combine peak→gene link evidence into one calibrated probability (gradient boosting, leave-one-chromosome-out) | — | [`pgboost_skill.py`](../../../Scripts/pgboost_skill.py) |
| `sce2g` | scE2G workbench: set up / configure / check / run model training with crowdsourced features, describe predictions, CRISPR_comparison-style benchmark | [scE2G workbench: train, check and benchmark E2G models (`sce2g`)](enhancer-gene.md#sce2g-workbench-train-check-and-benchmark-e2g-models-sce2g) | [`sce2g_workbench_skill.py`](../../../Scripts/sce2g_workbench_skill.py) |
| `sce2g-kg` | scE2G element→gene linkages → local KG (bulk, adaptive region tiling, resumable, runs to completion regardless of size) | — | [`sce2g_kg_skill.py`](../../../Scripts/sce2g_kg_skill.py) |
| `sce2g-pipeline` | Port of EngreitzLab/scE2G: the full single-cell E2G pipeline (Kendall, ARC-E2G, v3 models with qnorm/TPM filter, QC, CRISPR benchmark, training) in Python. | [scE2G pipeline (`sce2g-pipeline`)](enhancer-gene.md#sce2g-pipeline-sce2g-pipeline) | [`sce2g_pipeline_skill.py`](../../../Scripts/sce2g_pipeline_skill.py) |
| `sce2g-predict` | scE2G-style enhancer→gene features from paired single-cell ATAC + RNA (Kendall + ABC) | — | [`sce2g_predict_skill.py`](../../../Scripts/sce2g_predict_skill.py) |
| `se-targets` | Super-enhancer → target-gene pipeline (ENCODE) | [Super-enhancer → target-gene pipeline](enhancer-gene.md#super-enhancer--target-gene-pipeline) | [`se_target_pipeline.py`](../../../Scripts/se_target_pipeline.py) |

### MPRA and STARR-seq

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `bcalm` | Barcode-level MPRA activity with empirical-Bayes moderated statistics (BCalm approach) | — | [`bcalm_skill.py`](../../../Scripts/bcalm_skill.py) |
| `eqtl-mpra` | Tissue-filtered fine-mapped eQTLs and MPRA-vs-eQTL correlation / sign-concordance | — | [`eqtl_mpra_concordance_skill.py`](../../../Scripts/eqtl_mpra_concordance_skill.py) |
| `mpra` | MPRA / STARR / BlueSTARR retrieval and analysis | [MPRA / STARR / BlueSTARR](mpra-starr.md#mpra--starr--bluestarr) | [`mpra_data_skills.py`](../../../Scripts/mpra_data_skills.py) |
| `mpraflow` | MPRA counts: barcode->oligo assignment, count normalisation, replicate aggregation, allelic skew (port of kircherlab/MPRAsnakeflow) | [MPRA library design → counts → barcode QC](mpra-starr.md#mpra-library-design--counts--barcode-qc) | [`mpra_snakeflow_skill.py`](../../../Scripts/mpra_snakeflow_skill.py) |
| `mpralib` | MPRA barcode QC: outlier detection (global / oligo-specific / large-expression), replicate agreement, per-barcode activity (port of kircherlab/MPRAlib) | [MPRA library design → counts → barcode QC](mpra-starr.md#mpra-library-design--counts--barcode-qc) | [`mpralib_skill.py`](../../../Scripts/mpralib_skill.py) |
| `oligo` | MPRA oligo library design: tiling, REF/ALT variant oligos, synthesis filters, adapters (port of kircherlab/MPRAOligoDesign) | [MPRA library design → counts → barcode QC](mpra-starr.md#mpra-library-design--counts--barcode-qc) | [`mpra_oligo_design_skill.py`](../../../Scripts/mpra_oligo_design_skill.py) |
| `scqers` | scQers: quantitative single-cell enhancer reporters. Paired oBC (which element) + mBC (how hard it drives) from barcode extraction through bootstrap activity, permutation cell-type specificity, and FDR-controlled CRE calls | [scQers — quantitative single-cell enhancer reporters](mpra-starr.md#scqers--quantitative-single-cell-enhancer-reporters) | [`scqers_skill.py`](../../../Scripts/scqers_skill.py) |
| `starrseq` | STARR-seq allelic test (mpralm clean-room rewrite) | [STARR-seq allelic test](mpra-starr.md#starr-seq-allelic-test) | [`starr_seq_skill.py`](../../../Scripts/starr_seq_skill.py) |

### Bulk genomics, RNA-seq and proteomics

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `biosample-census` | Systematic ENCODE + IGVF census of one biosample (GM12878, K562, ...) with tables and plots | [ENCODE bulk-genomics pipeline](bulk-genomics-proteomics.md#encode-bulk-genomics-pipeline) | [`biosample_census_skill.py`](../../../Scripts/biosample_census_skill.py) |
| `encode` | ENCODE ChIP/ATAC/DNase/Hi-C/ChIA-PET pipeline | [ENCODE bulk-genomics pipeline](bulk-genomics-proteomics.md#encode-bulk-genomics-pipeline) | [`encode_pipeline.py`](../../../Scripts/encode_pipeline.py) |
| `proteomics` | Proteomics & PPI: BioGRID/IntAct/HuRI/Reactome/KEGG/IGVF integration, KG, viz, literature survey | [Proteomics & PPI knowledge graph](bulk-genomics-proteomics.md#proteomics--ppi-knowledge-graph) | [`proteomics_skill.py`](../../../Scripts/proteomics_skill.py) |
| `rnaseq` | Bulk RNA-seq QC / PCA / DEG / DEG→cCRE linkage | [Bulk RNA-seq analysis](bulk-genomics-proteomics.md#bulk-rna-seq-analysis) | [`rnaseq_analysis_skill.py`](../../../Scripts/rnaseq_analysis_skill.py) |

### Knowledge graph, warehouse and networks

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `grn` | Gene regulatory network (dEx) + protein-variant effects from the Catalog | — | [`catalog_grn_skill.py`](../../../Scripts/catalog_grn_skill.py) |
| `kg` | IGVF Knowledge Graph multi-hop traversal | [Knowledge Graph traversal](knowledge-graph.md#knowledge-graph-traversal) | [`kg_traversal_skill.py`](../../../Scripts/kg_traversal_skill.py) |
| `kg-integrate` | Merge the PPI-KG and the ArangoDB mirror into the IGVF integrated KG, incrementally, on one identifier space | — | [`kg_integrate_skill.py`](../../../Scripts/kg_integrate_skill.py) |
| `kg-mirror` | Local IGVF KG mirror (Arango -> Parquet + DuckDB) | [Local IGVF KG mirror (Arango → DuckDB)](knowledge-graph.md#local-igvf-kg-mirror-arango--duckdb) | [`kg_mirror_skill.py`](../../../Scripts/kg_mirror_skill.py) |
| `network` | Network integration — clean-room MILP for context-specific subnetworks (CARNIVAL / Steiner) | [Network integration (clean-room MILP — CARNIVAL + Steiner)](knowledge-graph.md#network-integration-clean-room-milp--carnival--steiner) | [`network_integration_skill.py`](../../../Scripts/network_integration_skill.py) |
| `portal-kg` | Portal → local Knowledge Graph ETL | [Portal → local KG ETL](knowledge-graph.md#portal--local-kg-etl) | [`portal_to_kg_skill.py`](../../../Scripts/portal_to_kg_skill.py) |
| `warehouse` | Central DuckDB warehouse — Silver tier for the IGVF integrated data layer | [Integrated data warehouse (DuckDB Silver tier)](knowledge-graph.md#integrated-data-warehouse-duckdb-silver-tier) | [`warehouse_skill.py`](../../../Scripts/warehouse_skill.py) |

### Agent, provenance and utilities

Skills without a worked example in the guide link to their module; `igvfagent <skill> --help` describes them.

| Skill | What it does | Documented in | Module |
|---|---|---|---|
| `files` | Read and write plain-text files inside the workspace (`write_text_file`, `read_text_file`): never secrets, code or extension directories; content never read from stdin | [Long-running jobs](../jobs.md) | [`workspace_files.py`](../../../Scripts/workspace_files.py) |
| `repro` | Per-paper reproduction records across attempts: outcome, route, agreement with the paper and the authors' outputs, verifier, a self-contained HTML report; publish; post to discussion.genohub.org | [Reproducing a paper](../paper-reproduction.md#records-the-reproductions-page-and-the-forum) | [`reproductions.py`](../../../Scripts/reproductions.py) |
| `paper-code` | Reproduce a paper by running the authors' own code: repository from the Code Availability statement, pinned commit, R/Python environment, unmodified execution, comparison with the authors' rendered output, replay check, per-section report | [Reproducing a paper](../paper-reproduction.md) | [`paper_code_skill.py`](../../../Scripts/paper_code_skill.py) |
| `job` | Durable agent jobs: long analyses planned on disk, gated by harness checks, verified, resumable | [Long-running jobs](../jobs.md) | [`agent_jobs.py`](../../../Scripts/agent_jobs.py) |
| `artifact` | Read back reports/manifests the agent produced (workspace-contained: read / grep / ls) | — | [`artifact_read_skill.py`](../../../Scripts/artifact_read_skill.py) |
| `bench` | Paper → reproduction benchmark. Resolve a publication from a title / URL / DOI / PMID / author+journal+year, harvest its Data Availability statement + accessions from the full text, route them onto an IGVFagent analysis chain, and scaffold a runnable Benchmarks/<paper-id>/ that concordance.py scores. Subcommands: resolve, harvest, route, scaffold, run, score, report, pipeline, selftest, list-routes. | — | [`benchmark_skill.py`](../../../Scripts/benchmark_skill.py) |
| `eval-tiers` | Tier 2 (planning / tool selection) and Tier 3 (conclusion validity) evaluation | — | [`eval_tiers_skill.py`](../../../Scripts/eval_tiers_skill.py) |
| `extauthor` | Author new IGVFagent tools/skills from the agent itself (opt-in: IGVF_ALLOW_AGENT_AUTHORING=1) | — | [`ext_author_skill.py`](../../../Scripts/ext_author_skill.py) |
| `mcp` | Serve the skill registry over the Model Context Protocol so other agents can call IGVF skills | — | [`mcp_server_skill.py`](../../../Scripts/mcp_server_skill.py) |
| `playbook-freeze` | Freeze a recorded session into a deterministic playbook with pinned args and artefact hashes | — | [`playbook_freeze_skill.py`](../../../Scripts/playbook_freeze_skill.py) |
| `project` | Projects + permanent searchable history: file analyses into a named project, recall everything ever produced about an accession, full-text search every past run. Nothing recorded is ever deleted. | — | [`project_skill.py`](../../../Scripts/project_skill.py) |
| `skillcard` | Machine-readable skill specifications (skill cards) with validation status derived from the benchmarks | — | [`skillcard_skill.py`](../../../Scripts/skillcard_skill.py) |
| `submit` | Submit to the IGVF Portal: preflight audit against the recurring DACC findings, revoked-input repair, local file validation, schema templates | — | [`igvf_submission_skill.py`](../../../Scripts/igvf_submission_skill.py) |
| `upstream` | Provenance of absorbed upstream projects: which revision each was built against, its licence, whether it was reimplemented or wrapped, and what has changed there since | — | [`upstream_skill.py`](../../../Scripts/upstream_skill.py) |

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
