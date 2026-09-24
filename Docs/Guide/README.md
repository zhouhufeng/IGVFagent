# IGVFagent documentation

Start with the [project README](../../README.md) for what IGVFagent is and how
to try it. These pages hold the detail.

## Using IGVFagent

| Page | Read it when you want to… |
|---|---|
| [Installation and first run](installation.md) | install locally (pip, pipx or Docker), configure credentials, run the smoke test |
| [The browser UI, accounts and history](web-ui.md) | know what each tab does, how sign-in works, and how past results and projects are kept |
| [LLM backends](llm-backends.md) | choose a model (Claude, OpenAI, Ollama, Claude Code), keep results consistent across models, or drive IGVFagent from another agent |
| [Reproducing a paper](paper-reproduction.md) | run the authors' own code from the paper's repository, compare every printed value and figure with their rendering, repair only the environment |
| [Long-running jobs](jobs.md) | reproduce a paper or run a pipeline as a background job that plans, checks every stage, verifies and resumes |
| [Skills](skills/README.md) | find the skill for a task: the full index of every `igvfagent` command |

## Skills by area

| Page | Covers |
|---|---|
| [Finding and retrieving data](skills/portal-and-data-access.md) | IGVF Portal and Catalog, processed-first Portal lineage, dataset explanation, literature, GEO, Synapse, ChIP-Atlas, Perturbation Catalogue, pathway and TF databases |
| [Variants and variant effects](skills/variants.md) | variant annotation, advanced variant analysis, cCRE and FAVOR, regulatory regions, MaveDB, SGE, exCALIBR calibration |
| [Single-cell, multiome and spatial](skills/single-cell.md) | scRNA-seq, multiome, SPLiT-seq, SHARE-seq, snMCT-seq, cell hashing, Tabula Sapiens, IGVF single-cell pipeline, principal pseudobulks, Spatial-ATAC-Hi-C |
| [CRISPR screens and Perturb-seq](skills/crispr-and-perturbation.md) | pooled screens from raw reads, base editing, CRISPRi, Flow-FISH, IGVF-CRISPR and pinellolab Perturb-seq pipelines, SCEPTRE, TF Perturb-seq, jamborees |
| [Enhancer–gene linkage (E2G)](skills/enhancer-gene.md) | ABC, ENCODE-rE2G, scE2G, eQTL and GWAS benchmarks, super-enhancer targets, E2G QC and Portal submission |
| [MPRA and STARR-seq](skills/mpra-starr.md) | library design, counts, barcode QC, IGVF MPRA standards, STARR-seq allelic tests, scQers |
| [Bulk genomics, RNA-seq and proteomics](skills/bulk-genomics-proteomics.md) | ENCODE bulk pipeline, bulk RNA-seq, the protein-interaction knowledge graph |
| [Knowledge graph, warehouse and networks](skills/knowledge-graph.md) | the local KG that grows with use, the Catalog mirror, traversal, the DuckDB warehouse, network integration |

## Project

| Page | Contents |
|---|---|
| [Architecture](architecture.md) | how the agent is put together, and the repository layout |
| [Reproducibility benchmarks](benchmarks.md) | papers reproduced from public data, and how to run the suite |
| [Extending IGVFagent](extending.md) | your own tools, skills, prompt skills and playbooks |
| [Operating the hosted deployment](deployment.md) | rebuilding the container, credentials, hot copies (operators) |
| [Security](security.md) | credentials and data handling |
| [References and attribution](references.md) | upstream repositories with pinned commits, methods papers, licence policy |
| [What's new](whats-new.md) | the full change log |

## Other documents

- [`Docs/THREAT_MODEL.md`](../THREAT_MODEL.md): what leaves your machine under each backend.
- [`Docs/EVALUATION.md`](../EVALUATION.md): how the agent is evaluated beyond the benchmarks.
- [`Docs/AUTH.md`](../AUTH.md): accounts and sign-in on the hosted deployment.
- [`Docs/PROJECT_SCOPE.md`](../PROJECT_SCOPE.md): what the agent is and is not.
- [`Deploy/README.md`](../../Deploy/README.md): setting up a server.
- [`Benchmarks/README.md`](../../Benchmarks/README.md): the benchmark dashboard.
- [`Docs/Skills/`](../Skills/): the per-skill playbooks the agent reads.
- [`Scripts/README.md`](../../Scripts/README.md): direct `python3 Scripts/…` invocations.

Data overviews and surveys:

- [`Docs/IGVF_PORTAL_DATA_OVERVIEW.md`](../IGVF_PORTAL_DATA_OVERVIEW.md) and [`Docs/IGVF_FRONT_PAGE_DATA_SUMMARY.md`](../IGVF_FRONT_PAGE_DATA_SUMMARY.md): what the IGVF Portal holds.
- [`Docs/IGVF_CATALOG_SMOKE_ANALYSIS.md`](../IGVF_CATALOG_SMOKE_ANALYSIS.md): a smoke analysis of the IGVF Catalog.
- [`Docs/ENCODE_DATA_OVERVIEW.md`](../ENCODE_DATA_OVERVIEW.md) and [`Docs/ENCODE_SMOKE_ANALYSIS.md`](../ENCODE_SMOKE_ANALYSIS.md): the same for ENCODE.
- [`Docs/PATHWAYDB.md`](../PATHWAYDB.md): the local pathway database.
- [`Docs/CRISPRi/CRISPRI_BIOINFORMATICS_APPLICATIONS.md`](../CRISPRi/CRISPRI_BIOINFORMATICS_APPLICATIONS.md): CRISPRi bioinformatics applications.
- [`Docs/SingleCell/Multiome10_survey.md`](../SingleCell/Multiome10_survey.md) and [`Docs/SingleCell/Portal10_survey.md`](../SingleCell/Portal10_survey.md): single-cell and multiome dataset surveys.
- [`Docs/DEPLOYMENT.md`](../DEPLOYMENT.md): legacy workstation deployment notes, superseded by [Installation](installation.md) and [Operating the hosted deployment](deployment.md).

## Where the old README sections went

The README used to hold everything in one 3,700-line page. Every one of its sections is still here; this table maps each old heading to where it lives now.

| Old README section | Now in |
|---|---|
| What IGVF Agent can do | [README](../../README.md#what-igvf-agent-can-do) |
| &nbsp;&nbsp;&nbsp;&nbsp;🎬 Video demos | [README](../../README.md#what-igvf-agent-can-do) |
| 🌐 Try it online — no install required | [README](../../README.md#-try-it-online--no-install-required) |
| What's new | [whats-new.md](whats-new.md) |
| Architecture at a glance | [architecture.md](architecture.md#architecture-at-a-glance) |
| Table of contents | [README.md](README.md) |
| Capabilities | [skills/README.md](skills/README.md#capabilities) |
| Repository layout | [architecture.md](architecture.md#repository-layout) |
| Quick start | [installation.md](installation.md#quick-start) |
| Recreating the hosted container (and why you must) | [deployment.md](deployment.md#recreating-the-hosted-container-and-why-you-must) |
| &nbsp;&nbsp;&nbsp;&nbsp;Check first, then rebuild | [deployment.md](deployment.md#check-first-then-rebuild) |
| &nbsp;&nbsp;&nbsp;&nbsp;It refuses while analyses are running, on purpose | [deployment.md](deployment.md#it-refuses-while-analyses-are-running-on-purpose) |
| &nbsp;&nbsp;&nbsp;&nbsp;Optional: install BEAN for base-editing screens | [deployment.md](deployment.md#optional-install-bean-for-base-editing-screens) |
| &nbsp;&nbsp;&nbsp;&nbsp;Credentials and the shared-deployment settings | [deployment.md](deployment.md#credentials-and-the-shared-deployment-settings) |
| &nbsp;&nbsp;&nbsp;&nbsp;When a hot copy is the right answer instead | [deployment.md](deployment.md#when-a-hot-copy-is-the-right-answer-instead) |
| Choosing an LLM backend | [llm-backends.md](llm-backends.md#choosing-an-llm-backend) |
| &nbsp;&nbsp;&nbsp;&nbsp;Local LLM (free, private, offline) | [llm-backends.md](llm-backends.md#local-llm-free-private-offline) |
| Configuration | [installation.md](installation.md#configuration) |
| Consistent results across LLM backends | [llm-backends.md](llm-backends.md#consistent-results-across-llm-backends) |
| The growing local knowledge graph + database | [skills/knowledge-graph.md](skills/knowledge-graph.md#the-growing-local-knowledge-graph--database) |
| Accounts, approval and sign-in | [web-ui.md](web-ui.md#accounts-approval-and-sign-in) |
| Projects and permanent history | [web-ui.md](web-ui.md#projects-and-permanent-history) |
| Smoke test | [installation.md](installation.md#smoke-test) |
| Skill catalog and usage | [skills/README.md](skills/README.md#all-skills) |
| &nbsp;&nbsp;&nbsp;&nbsp;IGVF / ENCODE / Knowledge Graph client | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#igvf--encode--knowledge-graph-client) |
| &nbsp;&nbsp;&nbsp;&nbsp;Catalog / Portal / ENCODE overviews | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#catalog--portal--encode-overviews) |
| &nbsp;&nbsp;&nbsp;&nbsp;Variant annotation | [skills/variants.md](skills/variants.md#variant-annotation) |
| &nbsp;&nbsp;&nbsp;&nbsp;Advanced variant analysis | [skills/variants.md](skills/variants.md#advanced-variant-analysis) |
| &nbsp;&nbsp;&nbsp;&nbsp;Single-cell, multiome, specialized assays | [skills/single-cell.md](skills/single-cell.md#single-cell-multiome-specialized-assays) |
| &nbsp;&nbsp;&nbsp;&nbsp;Cross-source multiome survey | [skills/single-cell.md](skills/single-cell.md#cross-source-multiome-survey) |
| &nbsp;&nbsp;&nbsp;&nbsp;Parse SPLiT-seq pipeline | [skills/single-cell.md](skills/single-cell.md#parse-split-seq-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;Pathway databases (KEGG + Reactome + WikiPathways, integrated locally) | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#pathway-databases-kegg--reactome--wikipathways-integrated-locally) |
| &nbsp;&nbsp;&nbsp;&nbsp;scQers — quantitative single-cell enhancer reporters | [skills/mpra-starr.md](skills/mpra-starr.md#scqers--quantitative-single-cell-enhancer-reporters) |
| &nbsp;&nbsp;&nbsp;&nbsp;Enhancer–gene linkage | [skills/enhancer-gene.md](skills/enhancer-gene.md#enhancergene-linkage) |
| &nbsp;&nbsp;&nbsp;&nbsp;MPRA / STARR / BlueSTARR | [skills/mpra-starr.md](skills/mpra-starr.md#mpra--starr--bluestarr) |
| &nbsp;&nbsp;&nbsp;&nbsp;MPRA library design → counts → barcode QC | [skills/mpra-starr.md](skills/mpra-starr.md#mpra-library-design--counts--barcode-qc) |
| &nbsp;&nbsp;&nbsp;&nbsp;Enhancer / regulatory-region cCRE annotation | [skills/variants.md](skills/variants.md#enhancer--regulatory-region-ccre-annotation) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPR screens from raw reads — three shapes, three tools | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispr-screens-from-raw-reads--three-shapes-three-tools) |
| &nbsp;&nbsp;&nbsp;&nbsp;Saturation genome editing (SGE) | [skills/variants.md](skills/variants.md#saturation-genome-editing-sge) |
| &nbsp;&nbsp;&nbsp;&nbsp;snMCT-seq — RNA and methylation from the same nuclei | [skills/single-cell.md](skills/single-cell.md#snmct-seq--rna-and-methylation-from-the-same-nuclei) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPRi / CRISPR-FACS / Perturb-seq | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispri--crispr-facs--perturb-seq) |
| &nbsp;&nbsp;&nbsp;&nbsp;SHARE-seq joint scATAC + scRNA QC | [skills/single-cell.md](skills/single-cell.md#share-seq-joint-scatac--scrna-qc) |
| &nbsp;&nbsp;&nbsp;&nbsp;STARR-seq allelic test | [skills/mpra-starr.md](skills/mpra-starr.md#starr-seq-allelic-test) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPRi Flow-FISH screen | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispri-flow-fish-screen) |
| &nbsp;&nbsp;&nbsp;&nbsp;cCRE, FAVOR, IGV-style browser views | [skills/variants.md](skills/variants.md#ccre-favor-igv-style-browser-views) |
| &nbsp;&nbsp;&nbsp;&nbsp;Start from what the IGVF Portal already computed (`processed lineage`) | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#start-from-what-the-igvf-portal-already-computed-processed-lineage) |
| &nbsp;&nbsp;&nbsp;&nbsp;Data illustration and interpretation | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#data-illustration-and-interpretation) |
| &nbsp;&nbsp;&nbsp;&nbsp;Reference skill (literature retrieval, validation, design) | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#reference-skill-literature-retrieval-validation-design) |
| &nbsp;&nbsp;&nbsp;&nbsp;Local IGVF KG mirror (Arango → DuckDB) | [skills/knowledge-graph.md](skills/knowledge-graph.md#local-igvf-kg-mirror-arango--duckdb) |
| &nbsp;&nbsp;&nbsp;&nbsp;Knowledge Graph traversal | [skills/knowledge-graph.md](skills/knowledge-graph.md#knowledge-graph-traversal) |
| &nbsp;&nbsp;&nbsp;&nbsp;Portal → local KG ETL | [skills/knowledge-graph.md](skills/knowledge-graph.md#portal--local-kg-etl) |
| &nbsp;&nbsp;&nbsp;&nbsp;ENCODE bulk-genomics pipeline | [skills/bulk-genomics-proteomics.md](skills/bulk-genomics-proteomics.md#encode-bulk-genomics-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Everything a portal holds for one cell line | [skills/bulk-genomics-proteomics.md](skills/bulk-genomics-proteomics.md#everything-a-portal-holds-for-one-cell-line) |
| &nbsp;&nbsp;&nbsp;&nbsp;Super-enhancer → target-gene pipeline | [skills/enhancer-gene.md](skills/enhancer-gene.md#super-enhancer--target-gene-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;GEO retrieval | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#geo-retrieval) |
| &nbsp;&nbsp;&nbsp;&nbsp;Bulk RNA-seq analysis | [skills/bulk-genomics-proteomics.md](skills/bulk-genomics-proteomics.md#bulk-rna-seq-analysis) |
| &nbsp;&nbsp;&nbsp;&nbsp;Proteomics & PPI knowledge graph | [skills/bulk-genomics-proteomics.md](skills/bulk-genomics-proteomics.md#proteomics--ppi-knowledge-graph) |
| &nbsp;&nbsp;&nbsp;&nbsp;Single-cell analysis (UMAP / t-SNE / Leiden / markers) | [skills/single-cell.md](skills/single-cell.md#single-cell-analysis-umap--t-sne--leiden--markers) |
| &nbsp;&nbsp;&nbsp;&nbsp;Perturbation Catalogue retrieval | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#perturbation-catalogue-retrieval) |
| &nbsp;&nbsp;&nbsp;&nbsp;MULTI-seq / Cell Hashing demultiplexing | [skills/single-cell.md](skills/single-cell.md#multi-seq--cell-hashing-demultiplexing) |
| &nbsp;&nbsp;&nbsp;&nbsp;Integrated data warehouse (DuckDB Silver tier) | [skills/knowledge-graph.md](skills/knowledge-graph.md#integrated-data-warehouse-duckdb-silver-tier) |
| &nbsp;&nbsp;&nbsp;&nbsp;Network integration (clean-room MILP — CARNIVAL + Steiner) | [skills/knowledge-graph.md](skills/knowledge-graph.md#network-integration-clean-room-milp--carnival--steiner) |
| &nbsp;&nbsp;&nbsp;&nbsp;Tabula Sapiens 2.0 (reference human cell atlas) | [skills/single-cell.md](skills/single-cell.md#tabula-sapiens-20-reference-human-cell-atlas) |
| &nbsp;&nbsp;&nbsp;&nbsp;Human Transcription Factors database | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#human-transcription-factors-database) |
| &nbsp;&nbsp;&nbsp;&nbsp;scE2G workbench: train, check and benchmark E2G models (`sce2g`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#sce2g-workbench-train-check-and-benchmark-e2g-models-sce2g) |
| &nbsp;&nbsp;&nbsp;&nbsp;TF Perturb-seq: calibrated effects to disease and GWAS (`tf-perturb`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#tf-perturb-seq-calibrated-effects-to-disease-and-gwas-tf-perturb) |
| &nbsp;&nbsp;&nbsp;&nbsp;eQTL enrichment benchmark for enhancer-gene predictors (`eqtl-enrich`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#eqtl-enrichment-benchmark-for-enhancer-gene-predictors-eqtl-enrich) |
| &nbsp;&nbsp;&nbsp;&nbsp;IGVF CRISPR Pipeline, early Nextflow version (`igvf-crispr-pipeline`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#igvf-crispr-pipeline-early-nextflow-version-igvf-crispr-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;SCEPTRE for IGVF MuData (`sceptre-igvf`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#sceptre-for-igvf-mudata-sceptre-igvf) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPRi-FlowFISH pipeline (`flowfish-pipeline`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispri-flowfish-pipeline-flowfish-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;Pooled CRISPR screen toolkit (`perturb-tools`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#pooled-crispr-screen-toolkit-perturb-tools) |
| &nbsp;&nbsp;&nbsp;&nbsp;IGVF single-cell pipeline (`igvf-sc-pipeline`) | [skills/single-cell.md](skills/single-cell.md#igvf-single-cell-pipeline-igvf-sc-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;GWAS E2G benchmark (`gwas-e2g`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#gwas-e2g-benchmark-gwas-e2g) |
| &nbsp;&nbsp;&nbsp;&nbsp;ENCODE-rE2G (`encode-re2g`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#encode-re2g-encode-re2g) |
| &nbsp;&nbsp;&nbsp;&nbsp;Perturb-seq pipeline, bulk_crispr_pipeline port (`bulk-crispr`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#perturb-seq-pipeline-bulk_crispr_pipeline-port-bulk-crispr) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPR module dev tools (`crisprdevtools`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispr-module-dev-tools-crisprdevtools) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPR-FG scPerturb-seq jamboree (`crispr-fg-jamboree`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispr-fg-scperturb-seq-jamboree-crispr-fg-jamboree) |
| &nbsp;&nbsp;&nbsp;&nbsp;Fishash Table 2 reproduction (`fishash-table2`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#fishash-table-2-reproduction-fishash-table2) |
| &nbsp;&nbsp;&nbsp;&nbsp;CRISPR SeqSpec (`crispr-seqspec`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#crispr-seqspec-crispr-seqspec) |
| &nbsp;&nbsp;&nbsp;&nbsp;Gasperini 2019 pipeline (`gasperini-pipeline`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#gasperini-2019-pipeline-gasperini-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;ABC pipeline (`abc-pipeline`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#abc-pipeline-abc-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;scE2G pipeline (`sce2g-pipeline`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#sce2g-pipeline-sce2g-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;IGVF CRISPR Perturb-seq pipeline (`crispr-pipeline`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#igvf-crispr-perturb-seq-pipeline-crispr-pipeline) |
| &nbsp;&nbsp;&nbsp;&nbsp;2nd CRISPR Jamboree analyses (`crispr-jamboree2`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#2nd-crispr-jamboree-analyses-crispr-jamboree2) |
| &nbsp;&nbsp;&nbsp;&nbsp;3rd CRISPR Jamboree single-cell pipeline (`crispr-jamboree3`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#3rd-crispr-jamboree-single-cell-pipeline-crispr-jamboree3) |
| &nbsp;&nbsp;&nbsp;&nbsp;4th CRISPR Jamboree inference benchmark (`crispr-jamboree4`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#4th-crispr-jamboree-inference-benchmark-crispr-jamboree4) |
| &nbsp;&nbsp;&nbsp;&nbsp;Principal pseudobulks (`principal-pseudobulks`) | [skills/single-cell.md](skills/single-cell.md#principal-pseudobulks-principal-pseudobulks) |
| &nbsp;&nbsp;&nbsp;&nbsp;E2G QC, predictions and IGVF Portal submission (`e2g-qc-predictions`) | [skills/enhancer-gene.md](skills/enhancer-gene.md#e2g-qc-predictions-and-igvf-portal-submission-e2g-qc-predictions) |
| &nbsp;&nbsp;&nbsp;&nbsp;Single-cell CRISPR differential expression (`sc-crispr-de`) | [skills/crispr-and-perturbation.md](skills/crispr-and-perturbation.md#single-cell-crispr-differential-expression-sc-crispr-de) |
| &nbsp;&nbsp;&nbsp;&nbsp;Spatial-ATAC-Hi-C (spatial 3D genome + chromatin accessibility) | [skills/single-cell.md](skills/single-cell.md#spatial-atac-hi-c-spatial-3d-genome--chromatin-accessibility) |
| &nbsp;&nbsp;&nbsp;&nbsp;ChIP-Atlas reprocessed peak archive | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#chip-atlas-reprocessed-peak-archive) |
| &nbsp;&nbsp;&nbsp;&nbsp;MaveDB mapping (incl. SGE cDNA path) | [skills/variants.md](skills/variants.md#mavedb-mapping-incl-sge-cdna-path) |
| &nbsp;&nbsp;&nbsp;&nbsp;Synapse / Sage Bionetworks retrieval | [skills/portal-and-data-access.md](skills/portal-and-data-access.md#synapse--sage-bionetworks-retrieval) |
| &nbsp;&nbsp;&nbsp;&nbsp;Assay calibration → ACMG/AMP evidence (exCALIBR) | [skills/variants.md](skills/variants.md#assay-calibration--acmgamp-evidence-excalibr) |
| Reproducibility benchmark suite | [benchmarks.md](benchmarks.md) |
| Extending IGVFagent | [extending.md](extending.md) |
| &nbsp;&nbsp;&nbsp;&nbsp;Add a custom tool (YAML manifest — no code) | [extending.md](extending.md#add-a-custom-tool-yaml-manifest--no-code) |
| &nbsp;&nbsp;&nbsp;&nbsp;Add a custom skill (Python subcommand) | [extending.md](extending.md#add-a-custom-skill-python-subcommand) |
| &nbsp;&nbsp;&nbsp;&nbsp;Add a prompt skill or playbook | [extending.md](extending.md#add-a-prompt-skill-or-playbook) |
| Deployment with LLM agents | [llm-backends.md](llm-backends.md#driving-igvfagent-from-another-agent) |
| &nbsp;&nbsp;&nbsp;&nbsp;Codex API | [llm-backends.md](llm-backends.md#recipes-by-agent) |
| &nbsp;&nbsp;&nbsp;&nbsp;Claude API | [llm-backends.md](llm-backends.md#recipes-by-agent) |
| &nbsp;&nbsp;&nbsp;&nbsp;Ollama local models | [llm-backends.md](llm-backends.md#recipes-by-agent) |
| Variant lists and your own data | [skills/variants.md](skills/variants.md#variant-lists-and-your-own-data) |
| Workstation notes | [installation.md](installation.md#workstation-notes) |
| Security | [security.md](security.md) |
| References | [references.md](references.md) |
| &nbsp;&nbsp;&nbsp;&nbsp;Reference GitHub repositories | [references.md](references.md#reference-github-repositories) |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Conventions for this table | [references.md](references.md#conventions-for-this-table) |
| &nbsp;&nbsp;&nbsp;&nbsp;Methods papers cited in the skills | [references.md](references.md#methods-papers-cited-in-the-skills) |
| &nbsp;&nbsp;&nbsp;&nbsp;License & attribution policy | [references.md](references.md#license--attribution-policy) |
| License | [README](../../README.md#citing-and-licence) |

---

[Project README](../../README.md)
