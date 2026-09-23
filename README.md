<img src="Docs/Figures/logo.png" alt="IGVF Agent" width="460">

An **auditable, local-execution** AI agent for discovering, retrieving, and
analyzing data from the [IGVF](https://igvf.org/) ecosystem (Portal, Catalog,
Knowledge Graph) and related public resources (ENCODE, FAVOR), with a built-in
**Plan → Action → Results → Evaluation** orchestration loop.

> **"Local" means tool execution, not inference.** Skills run as subprocesses
> on your machine and all files stay in the project folder. With a hosted model
> backend the prompt — including tool output the agent has read — still goes to
> the provider, and federated data access contacts public archives under every
> backend. See [`Docs/THREAT_MODEL.md`](Docs/THREAT_MODEL.md) for the per-mode
> data flow.

![IGVF Agent — system overview](Docs/Figures/IGVFagent_system_overview.png)

_End-to-end view: a knowledge graph and multi-omics data resources feed an orchestration-layer AI agent (Plan → Act → Observe → Refine) backed by short/long-term memory, an execution layer, post-processing, and human-in-the-loop review — producing downstream outputs such as variant prioritization, perturb-seq analysis, enhancer–gene mapping, and fine-mapping._

![IGVF Agent — architecture and skill topology](Docs/Figures/IGVF_agent_archetcture.png)

_Detailed five-layer architecture: user entry points (terminal, NL agent, browser UI) → agent runtime & tool dispatch → 88 skills / 270 typed tools grouped by domain → local persistence (filesystem + DuckDB warehouses) → upstream services. The `network` skill (highlighted) is the apex of the skill DAG — a clean-room MILP reimplementation of CORNETO that reads from the Silver + Bronze warehouses and writes inferred subnetworks back._

## What IGVF Agent can do

**Ask in plain language; it picks the method, runs it locally, and shows its
working.** 88 skills / 270 typed tools, run hosted or installed locally.

![What IGVF Agent can do](Docs/Figures/whatIGVFAgentcando.png)

- **Find and explain IGVF data** — search the Portal and Catalog by assay,
  tissue, gene or accession; say what a dataset actually contains before you
  download it; pull ENCODE, GEO, FAVOR, ChIP-Atlas, Synapse and Figshare too.
- **Analyse raw sequencing end to end** — FASTQ → counts → result, with the
  aligner and reference chosen from the dataset's own metadata rather than
  assumed.
- **CRISPR screens, routed by design not by title** — two-tail FACS screens,
  lettered-bin (A–F) screens, saturation genome editing, CRISPRi, FlowFISH,
  and **base-editing screens** via [crispr-bean](https://github.com/pinellolab/crispr-bean)'s
  masked matching (which recovers **62.5%** of reads where exact matching gets
  36.7%), handing off to the real `bean` binary for its Bayesian model.
- **Variants → function** — annotate lists, reach FAVOR and MaveDB, run
  MPRA/STARR-seq pipelines, link variants to genes through eQTL and enhancer
  evidence, and verify every claim against the record it came from.
- **Single-cell and multiome** — QC, clustering, cell typing, perturb-seq,
  SPLiT-seq, multiome, scE2G enhancer–gene links, spatial ATAC + Hi-C.
- **Regulatory networks** — GRN inference and a clean-room MILP
  reimplementation of CORNETO over the local warehouses.
- **One growing knowledge graph** — the IGVF Catalog mirror, a
  BioGRID/IntAct protein-interaction compendium and everything a session
  learns, merged onto the **same** vertices with per-edge provenance:
  **1.27M edges** and counting.
- **Help you submit to the Portal** — preflight your prediction sets and
  curated sets against the checks the monthly DACC meetings keep raising,
  before the meeting does.
- **Say when it cannot** — an unsupported assay, a revoked input, a model
  that needs a field IGVF does not publish. Refusing with a reason is treated
  as a result, not a failure.

Every run writes its files to your project folder, records which tool produced
what, and reports failed tool calls rather than smoothing them over.

### 🎬 Video demos

| | |
|---|---|
| [![Demo 1](https://img.youtube.com/vi/EQVwIEa-gVg/maxresdefault.jpg)](https://www.youtube.com/watch?v=EQVwIEa-gVg) | [![Demo 2](https://img.youtube.com/vi/c-CyIEArEK8/maxresdefault.jpg)](https://www.youtube.com/watch?v=c-CyIEArEK8) |
| ▶ <https://www.youtube.com/watch?v=EQVwIEa-gVg> | ▶ <https://www.youtube.com/watch?v=c-CyIEArEK8> |
| [![Demo 3](https://img.youtube.com/vi/DXmzSbrZC7E/maxresdefault.jpg)](https://youtu.be/DXmzSbrZC7E) | [![Demo 4](https://img.youtube.com/vi/OnZuOhwh4Qc/maxresdefault.jpg)](https://youtu.be/OnZuOhwh4Qc) |
| ▶ <https://youtu.be/DXmzSbrZC7E> | ▶ <https://youtu.be/OnZuOhwh4Qc> |

## 🌐 Try it online — no install required

**<https://igvfagent.genohub.org>**

A hosted IGVFagent runs on the project's
[Hcloud](https://github.com/zhouhufeng/HCloud) allocation, on the team's own
API key. Nothing to install, no Python environment to build, and **no API key
of your own** — open the URL and start asking questions.

Access is gated by a shared password. Email
[hufengzhou@g.harvard.edu](mailto:hufengzhou@g.harvard.edu) to request it.

**How the hosted instance differs from a local install:**

| | Hosted | Local install |
|---|---|---|
| Setup | none | `pip install -e '.[all]'` |
| LLM cost | paid by the project | your own API key, or free via Ollama |
| Model | Claude Sonnet 5 by default; Haiku 4.5 / Opus 5 / Fable 5 selectable | any backend — Anthropic, OpenAI, Ollama, vLLM, TGI, … |
| Run length | capped (iterations and tokens per turn) | uncapped |
| Workspace | **shared with other users** — see below | private to you |
| Your own data | not for anything sensitive | stays on your machine |

> ⚠️ **The hosted workspace is shared.** All signed-in users read and write one
> `Data/` and `Docs/` tree and one local knowledge graph, so analyses are
> visible to everyone else using the site. Treat it as a public demo sandbox:
> good for exploring IGVF/ENCODE public data, **not** for unpublished or
> sensitive datasets. Install locally for private work.

Heavy or long-running analyses (full multiome pipelines, large downloads) are
better run locally — see [Quick start](#quick-start). Operators: deployment
details are in [`Deploy/README.md`](Deploy/README.md).

## What's new

Most recent work first. Every item below is on `main` and exercised by the
benchmark suite or a worked example in this README.

| Area | Change |
|---|---|
| **Principal pseudobulks** | `principal-pseudobulks` ports EngreitzLab/generate-principal-pseudobulks: IGVF accession to QC-filtered per-cluster fragments, RNA matrices and scE2G config, spec-validated. |
| **Processed-first IGVF Portal lineage** | New `processed lineage` / `processed fetch` (`portal_lineage`, `processed_fetch`): from any accession, a directed walk of every Portal link. That covers analysis, principal, pseudobulk, model and prediction sets, the multiome partner, auxiliary sets, sample barcode maps, seqspecs and the published QC metrics. It gives a start-here table of processed files per product, with access and the Portal's own QC. The agent now calls it first, and `explain --download` fetches processed results instead of raw reads. |
| **4th CRISPR Jamboree benchmark** | The 2025 jamboree task: pairs preparation, sceptre / PerTurbo (plus glm.nb) inference, the pipeline's merged per-guide and per-element outputs, and AUPRC / AUROC against the control pairs across datasets. |
| **3rd CRISPR Jamboree pipeline** | The jamboree-3 IGVF single-cell Perturb-seq pipeline as 14 Python stages (configuration to dashboard) plus an end-to-end `run`, with the upstream quirks fixed and reproducible via `--upstream-compat`. |
| **2nd CRISPR Jamboree analyses** | Port of the 2024 IGVF CRISPR Jamboree: guide counting, CLEANSER/sceptre guide assignment, the six MuData inference modules (validated against the upstream outputs), the DESeq2 Perturb-seq simulator and the evaluation notebooks. |
| **IGVF CRISPR Perturb-seq pipeline** | The consortium CRISPR_Pipeline (seqspec to inference_mudata.h5mu, SCEPTRE/PerTurbo cis+trans, QC, evaluation, TF benchmark, dashboard) runs inside IGVFagent, with Python fallbacks for every external tool. |
| **scE2G pipeline** | scE2G rewritten in Python (`sce2g-pipeline`): Kendall, ARC-E2G, the four embedded v3 models, QC, benchmark and training. It reproduces upstream's chr22 fixture and all 11.5 M released K562 scores exactly. |
| **ABC pipeline** | Runs the whole ABC Snakemake workflow in Python (peaks through thresholded predictions and QC); its chr22 outputs match upstream's expected test outputs exactly. |
| **Gasperini 2019 pipeline** | The IGVF single-cell-like CRISPR pipeline on the Gasperini pilot: QC, MuData, guide assignment, cis SCEPTRE-style tests and the original NB test, benchmarked against GEO GSE120861. |
| **CRISPR SeqSpec** | Catalogue, validate (seqspec check) and index (kb / chromap / STARsolo read formats) the IGVF CRISPR assay seqspecs, and check FASTQs against them. |
| **Fishash Table 2** | `fishash-table2` reruns the CLEANSER and SCEPTRE barnyard species benchmark in Python, with both callers re-implemented, and scores it against the preprint's Table 2. |
| **CRISPR-FG jamboree** | All six jamboree notebooks as `crispr-fg-jamboree` subcommands: assay specs through guide calling, SCEPTRE-style cis tests and pyGenomeTracks links. |
| **CRISPR module dev tools** | Scaffolds IGVF CRISPR Nextflow modules in the pipeline's layout (with the upstream YAML and process-name bugs fixed) and checks existing modules for layout, conda-env and bin-script problems. |
| **Perturb-seq pipeline (bulk_crispr_pipeline port)** | Runs the IGVF-CRISPR Gasperini-pilot pipeline end to end in Python: kite guide counting, lane QC with doublets, MULTI-seq, cis gene sets, SCEPTRE-style tests, BH/Fisher results and browser tracks. |
| **ENCODE-rE2G** | Port of EngreitzLab/ENCODE_rE2G: features from ABC outputs, the nine pretrained models as embedded coefficients (upstream chr22 test output reproduced to 1e-15), CRISPR training and feature analysis. |
| **GWAS E2G benchmark** | Port of EngreitzLab/GWAS_E2G_benchmarking: tests E2G predictions against fine-mapped UK Biobank variants (enrichment, recall, curves) and silver-standard credible-set genes (with PoPS), with no R or bedtools needed. |
| **IGVF single-cell pipeline** | Ports every step of IGVF's multiome pipeline to Python: SHARE-seq barcode correction, TSS enrichment, kb/chromap wrappers, and RNA/ATAC/joint QC with barcode-rank knees. |
| **Pooled CRISPR screen toolkit** | `perturb-tools` ports IGVF-CRISPR/perturb-tools end to end (normalisation, replicate and sorting-screen LFC, QC, MAGeCK export, PoolQ reader, sgRNA library design), matching the upstream tutorials to six decimals. |
| **CRISPRi-FlowFISH pipeline** | Full port of EngreitzLab/crispri-flowfish (mapping, count tables, the R binned-MLE effect sizes, qPCR scaling, windows, peak calling, power), checked by a self-test with planted enhancers. |
| **SCEPTRE for IGVF MuData** | Port of sceptreIGVF: gRNA assignment plus SCEPTRE calibration, power and discovery tests on IGVF CRISPR MuData, with results written to `uns` under the upstream names, and an optional runner for the real R package. |
| **IGVF CRISPR Pipeline (early Nextflow)** | Port of IGVF-CRISPR/IGVF_CRISPR_Pipeline: seqspec to kb technology strings, kite guide references, kb mapping wrappers, a pure-Python guide counter, RNA QC and MuData assembly. |
| **scE2G workbench** | New `sce2g` skill wrapping [EngreitzLab/scE2G](https://github.com/EngreitzLab/scE2G) training with crowdsourced features (setup patches, Synapse feature table → `source_file` + external config + feature table, cluster/model rows, pre-flight `check`, snakemake `run`) and a clean-room [CRISPR_comparison](https://github.com/EngreitzLab/CRISPR_comparison) benchmark (pred_config semantics, AUPRC with bootstrap CI, precision at 70% recall, PR curves; writes the upstream configs too). Eight tools; self-test on a fake checkout plants a good and a random predictor and recovers the ranking. |
| **TF Perturb-seq → disease / GWAS** | New `tf-perturb` skill: port of the IGVF [tf_perturb_seq](https://github.com/IGVF/tf_perturb_seq) WG3 jamboree notebook. Calibrated direct / cis / trans tables → significant effects, top trans regulators, GWAS Catalog overlap with Fisher trait enrichment, scE2G links and GWAS SNPs inside E2G elements, optional ChIP-seq support; report, CSVs, figures. Vectorised overlaps, chunked+cached 22M-row trans filter, `mudata`/`pyBigWig` optional. Self-test plants FOXH1/SOX17 regulators, a T2D SNP block and an enhancer SNP and recovers all of them. |
| **Biosample census + hosted-agent fixes** | New `biosample-census` (`biosample_portal_census`): one command counts everything ENCODE and IGVF hold for a cell line via structured sample-term filters, with tables and plots. Re-authoring an agent extension under its own name is now an update; `explain_dataset` diagnoses a zero-hit search URL instead of writing an empty report; the failure banner shows the exception line, not the traceback header. |
| **Single-cell CRISPR DE** | New `sc-crispr-de` skill (method of [Gersbachlab-Bioinformatics/sc-crispr-de](https://github.com/Gersbachlab-Bioinformatics/sc-crispr-de), MIT). Every guide tested independently against the cells with **no detected guide**, one negative-binomial GLM per gene, then α-RRA aggregation to gene level calibrated on the non-targeting guides. A port was impossible, not merely inconvenient: upstream is R (`MASS::glm.nb`, brglm2, Seurat) on a SLURM array, and the container has neither. Validated on planted knockdowns — GATA1 ×0.25 (expected log2FC −2.00) recovered at −1.89 / −2.22 / −2.03, HBB ×0.40 (−1.32) at −1.18 / −1.25 / −1.38, 630/630 tests converged, NTCs unperturbed. `--apply-bias-reduction` (brglm2) is **not implemented** rather than approximated, and aggregation is α-RRA not FRACTEL — so the columns are `rra_*`, never `FRACTEL_*`. |
| **Tiling-screen deconvolution** | New `crispr-surf` skill. In a tiling screen each guide perturbs a **window**, not a point, so one element makes every nearby guide look active and the signal is smeared — `crispr-screen`, `gradient-screen`, `base-editing-screen` and `flowfish` all *score* guides, none deconvolve. Solves `y = A·beta` with an L1 penalty against a perturbation kernel, significance from the screen's own negative controls. Installing the upstream (AGPL-3.0) was tried first and is impossible — not on PyPI, and its source does not parse under Python 3 (`TabError`). Two elements planted in a 320-guide CRISPRi screen and smeared ±250 bp were each localised to **within ~10 bp** of their true centres, correct sign, exactly two regions, no false positives. |
| **Enhancer–gene prediction, generated** | Two new skills that *compute* links where `enhancer_gene_overview` and `sce2g_kg_pull` only retrieved them. `abc` implements Activity-by-Contact (Fulco 2019 / Nasser 2021, MIT) from candidate elements + bigWigs; Hi-C optional, falling back to the genome-wide power law. `sce2g-predict` adds the single-cell signal ABC cannot see — Kendall correlation between a peak's accessibility and a gene's expression across metacells. Tested on a case that separates them: a true link 50 kb away whose peak co-varies with the gene, against a decoy 2 kb from the TSS that is statistically independent. **ABC alone picks the decoy** (share 0.940 vs 0.051, purely from proximity); the correlation recovers the real link. scE2G's *trained weights* are not reproduced, so the output is labelled features + a transparent score, never "the scE2G score". |
| **IGVF uniform pipeline, completed** | [IGVF/atomic-workflows](https://github.com/IGVF/atomic-workflows) has exactly **two** modules — `igvf-kallisto-bustools` and `igvf-chromap`. kb-python was already wired into `raw_pipeline_run`, so half the official pipeline already ran here; **chromap 0.3.2-r518** is now built into the image (1.5 MB) and completes the tool set: kb for RNA, chromap for ATAC. Its 20 assay directories are *not* a second gap — `spec.md` and `README.md` are 0 bytes upstream. Also fixed: `--workflow nac` was offered in the CLI but could never run, because kb requires `-c1`/`-c2` capture lists that were never passed and `kb ref` was never told to build them. |
| **Metabolic labelling** | New `scnt-seq` skill (method of Qiu *et al.* 2020). Splits new from old RNA by 4sU **T>C conversions** read from a BAM's MD tags — no new binaries, no reference FASTA. Deliberately separate from kb's `nac` workflow, which IGVFagent also has: `nac` infers nascent from **intron content**, this measures **chemical labelling**, and an intronless transcript made an hour ago is new here and mature there. Strand is handled — the conversion presents as A>G on a reverse-strand read, and scoring it as T>C turns real labelling events into background. 400 synthetic reads classified with no errors. |
| **Barcode-level MPRA** | New `bcalm` skill (approach of [kircherlab/BCalm](https://github.com/kircherlab/BCalm)). `mpra_activity` sums an element's barcodes before testing, discarding the spread among them — the best evidence about how noisy that measurement is. Two elements planted with the **same** log2FC of 1.0, differing only in barcode agreement, came out at FDR 7.46e-80 and 1.3e-10: seventy orders of magnitude apart. Per-element variances are moderated by empirical Bayes so a few-barcode element is shrunk toward the trend rather than trusted on its own estimate. |
| **Genome-editing outcomes** | **CRISPResso2 2.3.4** installed (in `/opt/bean-venv` beside BEAN — both need numpy<2) and exposed as `crispresso`. This closes half a gap the base-editing diagnostics reported: bystander edit deconvolution needs CRISPResso2 alignment of the reporter allele. The other half stands — IGVF publishes no reporter column — but the paper's own Zenodo deposit does, so bystander analysis is now reachable there. Validated in a throwaway container before touching the Dockerfile: `bean --help` and `CRISPResso --help` both answer at numpy 1.26.4. |
| **Guide mapping, standalone** | New `guide-map` skill giving CRISPR-Correct's capability. Every primitive already existed — `build_matcher`, `mask_sequence`, the 1-Hamming barcode maps — but embedded inside the screen pipelines with no way to point them at an arbitrary FASTQ and library. Keeps the two kinds of imperfect match separate, because conflating them is a real error: a Hamming budget for **sequencing error**, masking for **base editing**. On 230 reads of known composition a 1-mismatch budget rescued exactly the 100 recoverable reads and left exactly the 30 junk reads unmapped. |
| **Normalisation** | New `sctransform` skill — the one part of Seurat genuinely missing, the rest of its workflow already being ported via scanpy (including Seurat 5 WNN). Regularised NB regression with Pearson residuals (Hafemeister & Satija 2019). The **regularisation** is the contribution, not the NB fit: per-gene θ is noisy, so parameters are smoothed across genes of similar mean expression. Across 300 cells varying 10-fold in depth, mean \|corr(expression, log depth)\| fell **0.175 → 0.063**. Output is residuals, **not counts** — it must not be fed to `sc_crispr_de_test` or `mpra_activity`. |
| **Paper routing** | `"hic"` is a substring of `"which"`. Assay terms were matched with `str.count`, so the Hi-C term fired on essentially every English paper — a VAMP-seq study of secreted proteins routed to the Spatial-ATAC-Hi-C chain, with a scaffolded benchmark, on the strength of one word. Terms now match on word boundaries, verified in both directions. |
| **Transfer ceilings** | The download cap was ten literals across ten skills — 0.3, 1, 2, 2, 2, 5, 5, 20, 20 and 100 GB — so which limit you hit depended on the skill you entered through, and three sat below the size of a single ordinary FASTQ. One shared default now (200 GB, `IGVF_MAX_DOWNLOAD_GB`), plus a free-disk preflight, and a file skipped for size is **logged** rather than recorded only in a CSV column while the run reports success. |
| **One knowledge graph** | New `kg-integrate` skill merges the three graphs that existed side by side without joining: the integrated KG the UI shows (gene nodes keyed by symbol), the proteomics PPI-KG (1.29M BioGRID interactions, also symbols, its `id_map` crosswalk **empty**), and the parquet mirror of the Catalog's ArangoDB (11.5M `proteins_proteins` edges keyed by Ensembl protein). Built on `_localstore`'s own `upsert_node`/`upsert_edge`, so merged edges land on the same vertices as everything the agent records. A 432,448-alias identity index resolves the three identifier spaces over a route that is entirely Ensembl keys — **94.9%** of BioGRID endpoints and **97.8%** of Catalog-PPI endpoints resolve; unresolved ones are kept and flagged, not dropped. Ambiguity is recorded rather than resolved by picking a winner (2,493 symbols map to >1 gene), and synonyms stay opt-in because **7,167** synonym strings are claimed by more than one gene. Incremental per parquet shard, so the still-growing mirror keeps contributing. Graph: 60,627 → **1,274,587** edges. |
| **Data submission** | Two new skills absorbed from the IGVF-Submit prototype, encoding the recurring findings of the monthly IGVF submission meetings (Oct 2025 – Aug 2026) as checks that run in seconds: `submit` (24 rules — missing `description`, `input_file_sets`, `derived_from`, `reference_files`, `file_format_specifications`, `analysis_step_version`, and inputs gone revoked/archived) and `portal-qc` (the Portal's own audit facets by lab and type, plus a measured provenance blind spot). Every rule cites the meeting that produced it, so a failure traces to a conversation rather than reading as an opinion. `crosscheck` answers the question a revoked input actually raises — do I reupload, or just repoint `derived_from`? |
| **Base-editing screens** | New `bean` skill, following [crispr-bean](https://github.com/pinellolab/crispr-bean)'s method: a base editor edits the guide's own locus, so the protospacer read back carries A>G (ABE) or C>T (CBE) changes and exact matching discards those reads — **36.7% of reads assigned vs 62.5%** on `IGVFDS6464SOVZ`. Masks the edited base on both sides as BEAN does rather than allowing free mismatches, detects the editor from the library's guide names *and* by measuring which masking recovers reads, and reports per-guide editing activity so a weakly-editing guide reads as underpowered rather than inactive. `--run-bean` hands those counts to the real BEAN (`create-screen` + `run sorting variant`) and reports its per-target posteriors, — **1,659 target posteriors on `IGVFDS6464SOVZ`**, whose four unsorted bins are found by shared guide library rather than by alias prefix (they are named `B*_18loci_Rep*_bulk`, not `18loci_uptake_*`). BEAN's activity-normalised model remains unavailable: it needs an `X_bcmatch` layer IGVF does not publish, so the posteriors come from BEAN's Normal model and are labelled as not activity-normalised. |
| **Screen routing** | `crispr-screen` and `gradient-screen` now **refuse** a base-editing library and name `bean`. The editor is stated only in the guide library — the readout and assay title are identical to an ordinary knockout screen — so routing metadata alone cannot catch it and the tool that would undercount does. |
| **Lettered-bin screens** | New `gradient-screen` skill for the **554** MeasurementSets sorted into bins A–F rather than two tails, scoring each construct's frequency-weighted mean bin. Its first real run called 27% of the library's 434 non-targeting controls significant; the null is now calibrated against control scatter with a depth-dependent variance model, giving a 0.46% control false-positive rate. |
| **Regulatory evidence** | Gene-centric linkage was querying a region endpoint with one capped request and no target-gene check, so it reported "no kidney records" for WT1 (4,458 target-gene rows, 405 kidney) and attributed other genes' edges to the queried gene. Now paginated to exhaustion (`page=`, since `skip=` is silently ignored), filtered on the row's own `gene` field, with predictions separated from significance and a record-level verifier that rejects cross-record field assembly. |
| **Shared-deployment safety** | The extension panel let **any visitor upload a `.py` that the server then executes**; now off unless the operator opts in. Uploads were pooled in one flat directory where a second `manifest.csv` silently overwrote the first — now per-session. Artefacts are scoped to the run that produced them instead of scanning a shared directory. |
| **Provenance** | The UI shows the build actually serving the page (a content hash of the loaded modules, matching `redeploy.sh`'s), because a stale container was previously indistinguishable from a fixed one. |
| **Enhancer annotation** | New `enhancer-annot` skill: annotate a list of regulatory **regions** (xlsx / BED / TSV / CSV) against the ENCODE SCREEN cCRE registry by interval overlap — element count, ranked class, union coverage, and nearest-element distance for non-overlapping regions. Registry V4 (2,348,854 elements) indexed locally for offline use. |
| **MPRA toolchain** | Three clean-room ports of the kircherlab suite completing the assay lifecycle — `oligo` (library design), `mpraflow` (counts → activity), `mpralib` (barcode QC + IGVF format validation). Standard-library only; verified bit-identical against six upstream scripts. |
| **Reproducibility** | New benchmark `rosen2025_mprasnakeflow` reproduces a **published IGVF artefact byte-for-byte** (210,660 / 210,660 values) and the source paper's stated complexity figures to the digit. Suite now 22 papers. |
| **IGVF file standards** | `mpralib validate` checks any file against the eight IGVF MPRA community schemas (Rosen *et al.* Supplementary Note S1). Every file `mpraflow` writes validates. |
| **Knowledge Graph explorer** | The UI listed graphs in discovery order behind a dropdown, so the IGVF integrated KG — the graph every tool writes into — was hidden behind Proteomics PPI-KG. Now leads, and all graphs are visible at once with what each is for. |
| **Deployment** | `Deploy/redeploy.sh` pulls, rebuilds and recreates the container, then **verifies** that the code actually running inside it matches the checkout. The app image bakes the package in at build time, so `git pull` alone never updated a running site. |
| **Catalog** | Fixed pathway node lookups returning the same record for every ID — `/api/pathways` keys on `id`, not `pathway_id`, and silently drops unknown query parameters. Every node lookup now verifies the response's identity against the request. |

## Architecture at a glance

IGVFagent ships **its own internal orchestrator** — a Plan → Action → Results
loop with an Evaluation Agent that cross-checks discoveries against the
literature, plus an explicit *skill retrieval* layer that selects the right
capability for each task. The agent is not just a bag of CLI tools; it
combines:

- **Inputs** — the IGVF Knowledge Graph (structured: genes, variants,
  regulatory elements, diseases, pathways) and the IGVF Portal (unstructured:
  publications, assay data from ATAC-seq / RNA-seq / CRISPRi / MPRA / 10x
  multiome / Parse SPLiT-seq, metadata, reports).
- **Skill retrieval** — KG queries, database accesses (FAVOR, VEP), file
  parsing (ATAC, RNA), coding tools, and a literature skill that pulls from
  PubMed / bioRxiv / arXiv / Semantic Scholar / OpenAlex.
- **Action execution loop** — KG queries, database calls, coding, and
  literature retrieval feed *data reading*, *analysis*, *tool use*, and
  *error handling* steps.
- **Evaluation Agent** — cross-checks evidence and validates consistency,
  with an optional human-feedback channel.
- **Outputs** — variant scoring & interpretation, multi-omic integration,
  enhancer-gene mapping & GRNs, fine-mapping / GWAS, trajectory inference,
  and cross-tissue / cross-species analyses.
- **Responsible-AI guardrails** — accountability, data provenance, bias &
  fairness review, privacy & consent, AI disclosure, explainability &
  transparency, and system security are first-class citizens of the design.

The agent is also **CLI-first** at the skill layer: every capability exposes
shell-runnable subcommands so the same skills can also be driven by an
external orchestrator (Codex, Claude Code, Ollama-served Qwen, or any
LLM that can invoke `python3 Scripts/...` and read files). Reads, writes,
caches, and logs all stay inside the repository folder for auditability.

In short — **two ways to drive every skill, one shared contract**:

| Mode | Who drives the loop | Best for |
|---|---|---|
| Internal orchestrator | IGVFagent's Plan→Action→Results→Evaluation runner (in-process) | Multi-step integrative analyses with branching and cross-evidence checks |
| External orchestrator | Codex / Claude Code / Ollama / your own harness | Day-to-day shell-level use, scripted pipelines, and CI |

## Table of contents

- [What's new](#whats-new)
- [Try it online (hosted)](#-try-it-online--no-install-required)
- [Capabilities](#capabilities)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Recreating the hosted container (and why you must)](#recreating-the-hosted-container-and-why-you-must)
- [Configuration](#configuration)
- [Smoke test](#smoke-test)
- [Skill catalog and usage](#skill-catalog-and-usage)
  - [IGVF / ENCODE / Knowledge Graph client](#igvf--encode--knowledge-graph-client)
  - [Catalog / Portal / ENCODE overviews](#catalog--portal--encode-overviews)
  - [Variant annotation](#variant-annotation)
  - [Advanced variant analysis](#advanced-variant-analysis)
  - [Single-cell, multiome, specialized assays](#single-cell-multiome-specialized-assays)
  - [Cross-source multiome survey](#cross-source-multiome-survey)
  - [Parse SPLiT-seq pipeline](#parse-split-seq-pipeline)
  - [Enhancer–gene linkage](#enhancergene-linkage)
  - [MPRA / STARR / BlueSTARR](#mpra--starr--bluestarr)
  - [MPRA library design → counts → barcode QC](#mpra-library-design--counts--barcode-qc)
  - [Enhancer / regulatory-region cCRE annotation](#enhancer--regulatory-region-ccre-annotation)
  - [CRISPR screens from raw reads — three shapes, three tools](#crispr-screens-from-raw-reads--three-shapes-three-tools)
  - [Saturation genome editing (SGE)](#saturation-genome-editing-sge)
  - [snMCT-seq — RNA and methylation from the same nuclei](#snmct-seq--rna-and-methylation-from-the-same-nuclei)
  - [CRISPRi / CRISPR-FACS / Perturb-seq](#crispri--crispr-facs--perturb-seq)
  - [cCRE, FAVOR, IGV-style browser views](#ccre-favor-igv-style-browser-views)
  - [Data illustration and interpretation](#data-illustration-and-interpretation)
  - [Reference skill (literature retrieval, validation, design)](#reference-skill-literature-retrieval-validation-design)
  - [Knowledge Graph traversal](#knowledge-graph-traversal)
  - [Portal → local KG ETL](#portal--local-kg-etl)
  - [ENCODE bulk-genomics pipeline](#encode-bulk-genomics-pipeline)
  - [Super-enhancer → target-gene pipeline](#super-enhancer--target-gene-pipeline)
  - [GEO retrieval](#geo-retrieval)
  - [Bulk RNA-seq analysis](#bulk-rna-seq-analysis)
  - [Proteomics & PPI knowledge graph](#proteomics--ppi-knowledge-graph)
  - [Single-cell analysis (UMAP / t-SNE / Leiden / markers)](#single-cell-analysis-umap--t-sne--leiden--markers)
  - [Perturbation Catalogue retrieval](#perturbation-catalogue-retrieval)
  - [MULTI-seq / Cell Hashing demultiplexing](#multi-seq--cell-hashing-demultiplexing)
  - [Integrated data warehouse (DuckDB Silver tier)](#integrated-data-warehouse-duckdb-silver-tier)
  - [Network integration (clean-room MILP — CARNIVAL + Steiner)](#network-integration-clean-room-milp--carnival--steiner)
  - [ChIP-Atlas reprocessed peak archive](#chip-atlas-reprocessed-peak-archive)
  - [MaveDB mapping (incl. SGE cDNA path)](#mavedb-mapping-incl-sge-cdna-path)
  - [Synapse / Sage Bionetworks retrieval](#synapse--sage-bionetworks-retrieval)
  - [Assay calibration → ACMG/AMP evidence (exCALIBR)](#assay-calibration--acmgamp-evidence-excalibr)
- [Reproducibility benchmark suite](#reproducibility-benchmark-suite)
- [Extending IGVFagent](#extending-igvfagent)
  - [Add a custom tool (YAML manifest — no code)](#add-a-custom-tool-yaml-manifest--no-code)
  - [Add a custom skill (Python subcommand)](#add-a-custom-skill-python-subcommand)
  - [Add a prompt skill or playbook](#add-a-prompt-skill-or-playbook)
- [Deployment with LLM agents](#deployment-with-llm-agents)
  - [Codex API](#codex-api)
  - [Claude API](#claude-api)
  - [Ollama local models](#ollama-local-models)
- [Variant lists and your own data](#variant-lists-and-your-own-data)
- [Workstation notes](#workstation-notes)
- [Security](#security)
- [License](#license)

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
- **Skill-correctness benchmark suite** — twenty-two recent Nature / Cell /
  Science / Nat Genet / Nat Methods / Genome Biol papers reproduced directly
  from public data, each scored against machine-readable ground-truth checks.
  **Scope:** every `run.sh` invokes skills directly as shell commands, with no
  model in the loop, so the suite establishes that the *implementations* are
  correct. It does not evaluate planning, tool selection, or the validity of an
  end-to-end conclusion — see [`Docs/EVALUATION.md`](Docs/EVALUATION.md). Includes
  **five full single-cell / multiome local reproductions** that download the
  public data and run IGVFagent's real analytical chain end-to-end — Travaglini
  2020 lung (`sc-analyze`), Trevino 2021 cortex multiome (`multiome peak2gene`),
  deMULTIplex2 / Stoeckius cell hashing (`multiseq`), Rosenberg 2018
  SPLiT-seq CNS (`splitseq`), and Ma 2020 SHARE-seq skin (`share`) — with a
  2025 developing-neocortex multiome atlas (Wang, CELLxGENE) in progress. The
  single-cell loaders + QC these exercise are internalized into
  `Scripts/_scload.py` and the skills (see [`Benchmarks/`](Benchmarks/README.md)).
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
  core-code edits. See [Extending IGVFagent](#extending-igvfagent).

## Repository layout

```
IGVFagent/
├── README.md                            ← you are here
├── LICENSE
├── pyproject.toml                       ← installable package (`pip install -e .[all]`)
├── requirements.txt                     ← minimal pip-only fallback
├── Dockerfile  docker-compose.yml  .dockerignore
├── .env.example                         ← copy to .env and edit locally
├── .gitignore                           ← excludes Data/, Docs/<skill>/2*, .env, caches
│
├── Scripts/                             ← CLI skills + shared runtime
│   ├── cli.py                           ← unified `igvfagent <skill> <subcommand>` dispatcher
│   ├── _agent.py                        ← in-process Plan→Action→Results→Evaluation runtime
│   ├── _llm.py                          ← multi-backend LLM router (Anthropic / OpenAI /
│   │                                       Codex / Ollama / vLLM / TGI / Groq / Together /
│   │                                       DeepInfra / HuggingFace / claude_cli / codex_cli)
│   ├── _tools.py                        ← curated tool registry the runtime exposes
│   ├── _endpoints.py                    ← endpoint resolver (env-overridable)
│   ├── streamlit_app.py                 ← LM-Studio-style browser UI (`igvfagent ui`)
│   │
│   ├── igvf_client.py                   ← Portal / Catalog / KG / ENCODE HTTP client
│   ├── igvf_data_skills.py              ← Catalog / Portal / ENCODE overview + smoke summaries
│   ├── igvf_frontpage_summary.py        ← refresh front-page Portal + KG stats
│   ├── igvf_specialized_data_skills.py  ← specialized assay catalog
│   │
│   ├── annotate_variant_list.py         ← variant-list annotation against Catalog evidence
│   ├── advanced_variant_analysis.py     ← integrated variant scoring + logistic + report
│   │
│   ├── single_cell_data_skills.py       ← scRNA / scATAC discovery + example analysis
│   ├── singlecell_analysis.py           ← Scanpy pipeline: QC, PCA, UMAP, t-SNE,
│   │                                       Leiden, markers, publication figures
│   ├── multiome_10x_pipeline.py         ← 10x Multiome retrieval pipeline
│   ├── multiome_research_demo.py
│   ├── multiome_survey.py               ← cross-source survey: IGVF + ENCODE + GEO +
│   │                                       CELLxGENE + HCA + Zenodo, unified manifest + downloader
│   ├── portal_multiome.py  portal_scrna_10.py
│   ├── putamen_multiome_demo_analysis.py
│   ├── splitseq_pipeline.py             ← Parse SPLiT-seq end-to-end pipeline
│   │
│   ├── ccre_linkage_annotation_skills.py
│   ├── enhancer_gene_linkage_skills.py
│   ├── open4gene_skill.py               ← Open4Gene peak→gene hurdle-model linkage
│   │                                       (clean-room port of hbliu/Open4Gene)
│   ├── sceps_skill.py                   ← scEPS GWAS × single-cell neighborhood
│   │                                       d-statistic (clean-room port of Genentech/sceps)
│   ├── mpra_data_skills.py              ← MPRA / STARR / BlueSTARR
│   ├── mpra_oligo_design_skill.py       ← MPRA oligo LIBRARY DESIGN — tiling,
│   │                                       REF/ALT variant oligos, synthesis
│   │                                       filters, adapters (MPRAOligoDesign port)
│   ├── mpra_snakeflow_skill.py          ← MPRA COUNTS — barcode→oligo assignment,
│   │                                       normalisation, outlier removal, replicate
│   │                                       aggregation, allelic skew, complexity
│   ├── mpralib_skill.py                 ← MPRA BARCODE QC — 3 outlier detectors,
│   │                                       replicate agreement, IGVF format validation
│   ├── _mpra_schemas.py                 ← the 8 IGVF MPRA community file schemas
│   ├── enhancer_annotation_skill.py     ← enhancer / regulatory-REGION annotation
│   │                                       against the SCREEN cCRE registry
│   ├── crispri_data_skills.py           ← CRISPRi / CRISPR-FACS / Perturb-seq
│   ├── encode_pipeline.py               ← ChIP/ATAC/DNase/Hi-C/ChIA-PET pipeline
│   ├── se_target_pipeline.py            ← super-enhancer → target-gene pipeline
│   │
│   ├── kg_traversal_skill.py            ← IGVF Knowledge Graph multi-hop traversal
│   ├── portal_to_kg_skill.py            ← Portal → local KG ETL (SQLite mirror)
│   │
│   ├── geo_retrieval_skill.py           ← NCBI GEO retrieval (search / metadata / download)
│   ├── figshare_skill.py                ← figshare retrieval (article / files / download /
│   │                                       search; id, DOI, URL, or private /s/ token)
│   ├── rnaseq_analysis_skill.py         ← bulk RNA-seq QC / PCA / DEG / DEG→cCRE linkage
│   ├── proteomics_skill.py              ← BioGRID/IntAct/HuRI/Reactome/KEGG + IGVF protein
│   │                                       PPI knowledge graph + per-assay viz + lit survey
│   ├── perturbation_catalog_skill.py    ← Perturbation Catalogue: MAVE / CRISPR-screen /
│   │                                       Perturb-seq retrieval + per-gene pipeline
│   ├── multiseq_analysis_skill.py       ← MULTI-seq / Cell Hashing demultiplexing
│   │                                       (Python port of deMULTIplex2)
│   ├── warehouse_skill.py               ← central DuckDB warehouse (Silver tier
│   │                                       of the integrated data layer)
│   ├── network_integration_skill.py     ← Network integration: clean-room
│   │                                       cvxpy MILP for CARNIVAL +
│   │                                       prize-collecting Steiner
│   │
│   ├── enrichment_skill.py              ← GO / pathway ORA + preranked GSEA (gseapy)
│   ├── chipatlas_skill.py              ← ChIP-Atlas reprocessed peak archive client
│   ├── mavedb_mapping_skill.py         ← MaveDB → genomic coords (+ SGE cDNA path)
│   ├── synapse_skill.py                ← Synapse / Sage Bionetworks retrieval
│   │                                       (anonymous walk/search + PAT download)
│   ├── excalibr_skill.py               ← assay score → ACMG/AMP PS3/BS3 evidence
│   │                                       (clean-room port of rosstewart/exCALIBR)
│   ├── reference_skill.py               ← literature retrieval / validation / study design
│   └── data_illustration_interpretation.py
│
├── Benchmarks/                          ← 22-paper reproducibility suite
│   ├── README.md                        ← suite dashboard + per-paper results
│   ├── OPERATIONS_GUIDE.md  run_all.sh  concordance.py
│   └── <paper-id>/                       ← run.sh + expected.json + figures + README
│
├── Data/                                ← inputs + cached responses (gitignored)
│   ├── Input/VariantList/example_variants.csv
│   ├── Manifests/   Cache/   KG/        ← ETL outputs (gitignored)
│   └── Proteomics/  ← _versions.json + Sources/<src>/  + KG/proteomics.sqlite (gitignored)
│
└── Docs/
    ├── PROJECT_SCOPE.md                 ← what the agent is and isn't
    ├── DEPLOYMENT.md                    ← legacy deploy notes (this README supersedes)
    ├── IGVF_PORTAL_DATA_OVERVIEW.md     IGVF_CATALOG_SMOKE_ANALYSIS.md
    ├── IGVF_FRONT_PAGE_DATA_SUMMARY.md
    ├── ENCODE_DATA_OVERVIEW.md          ENCODE_SMOKE_ANALYSIS.md
    ├── Figures/                         ← architecture diagram, etc.
    ├── CRISPRi/CRISPRI_BIOINFORMATICS_APPLICATIONS.md
    ├── SingleCell/{Multiome10_survey.md, Portal10_survey.md}
    ├── Skills/                          ← per-skill playbooks (one per CLI module)
    │   ├── ADVANCED_VARIANT_ANALYSIS_SKILLS.md
    │   ├── IGVF_PORTAL_DATA_ANALYSIS.md         IGVF_FRONT_PAGE_SUMMARY_SKILLS.md
    │   ├── IGVF_SPECIALIZED_DATA_SKILLS.md      IGVF_KG_TRAVERSAL_SKILLS.md
    │   ├── PORTAL_TO_KG_SKILLS.md
    │   ├── 10X_MULTIOME_SKILLS.md               SINGLE_CELL_ANALYSIS_SKILLS.md
    │   ├── SPLITSEQ_SKILLS.md
    │   ├── ENHANCER_GENE_LINKAGE_SKILLS.md      CCRE_LINKAGE_FAVOR_SKILLS.md
    │   ├── MPRA_ANALYSIS_SKILLS.md              CRISPRI_ANALYSIS_SKILLS.md
    │   ├── ENCODE_PIPELINE_SKILLS.md            SE_TARGET_PIPELINE_SKILLS.md
    │   ├── GEO_RETRIEVAL_SKILLS.md              RNASEQ_ANALYSIS_SKILLS.md
    │   ├── PROTEOMICS_SKILLS.md
    │   ├── PERTURBATION_CATALOG_SKILLS.md
    │   ├── MULTISEQ_ANALYSIS_SKILLS.md
    │   ├── WAREHOUSE_SKILLS.md
    │   ├── NETWORK_INTEGRATION_SKILLS.md
    │   ├── SINGLECELL_ANALYSIS_SKILLS.md
    │   ├── REFERENCE_SKILLS.md
    │   └── DATA_ILLUSTRATION_INTERPRETATION_SKILLS.md
    ├── Logs/                            ← runtime logs (gitignored)
    └── <skill>/<timestamp>_<label>/     ← per-run outputs from every skill (gitignored)
```

Generated outputs (timestamped folders under `Docs/<Skill>/`, manifests under
`Data/Manifests/`, source dumps under `Data/Proteomics/Sources/`, the local
KG mirrors under `Data/KG/` and `Data/Proteomics/KG/`, and caches under
`Data/Cache/`) are gitignored. When a new built-in skill ships, three things
must update together: the script in `Scripts/`, its playbook in
`Docs/Skills/`, and this Repository layout block. (User-supplied skills and
tools are exempt — they are auto-discovered from `~/.igvfagent/` and
`UserExtensions/` with no registry edits; see
[Extending IGVFagent](#extending-igvfagent).)

## Quick start

> Just want to try IGVFagent? Skip all of this and use the hosted instance at
> **<https://igvfagent.genohub.org>** — no install, no API key.
> See [Try it online](#-try-it-online--no-install-required).

**Recommended — native pip install in a virtual env (best for local LLMs):**

```bash
git clone https://github.com/zhouhufeng/IGVFagent.git
cd IGVFagent
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -e '.[ui,llm,analysis]'   # core + UI + LLM SDKs + scanpy stack
# python3 -m pip install -e '.[all]'             # also adds [hic] + [motif] extras

igvfagent --version
igvfagent --help

# Optional: configure credentials / overrides locally (.env is gitignored)
cp .env.example .env
```

Then:

```bash
igvfagent ui                # browser UI at http://127.0.0.1:8501
igvfagent ask "Pull the comprehensive APOE evidence pack."
igvfagent kg gene APOE      # direct CLI to a single skill, no LLM in the loop
```

This is the right path if you already have Ollama running natively
(`ollama serve`) — IGVFagent reaches it at `http://localhost:11434/v1`
and gets full access to your host's RAM and any models you've already
pulled.

**Alternative — `pipx` for a global `igvfagent` (no venv activation):**

```bash
pipx install 'git+https://github.com/zhouhufeng/IGVFagent.git'
pipx inject igvfagent 'igvfagent[analysis,ui,llm]'
igvfagent --help        # works from any directory
```

**Alternative — Docker Compose (self-contained stack, includes Ollama in a container):**

Picks up Ollama-in-a-container by default. Useful for a clean demo
machine with no Python or Ollama set up; **less ideal** if you already
run Ollama natively (the in-container Ollama can't see the host's
models, and Docker Desktop's default 8 GB memory cap is too small for
30B+ models). Requires Docker Desktop running.

```bash
git clone https://github.com/zhouhufeng/IGVFagent.git
cd IGVFagent

# Option A: in-container Ollama (default)
docker compose up -d                          # build agent + start ollama service
docker compose --profile bootstrap up         # one-time: pull qwen3:8b into the container

# Option B: reach back to your host's Ollama (recommended if you already have it)
echo 'OLLAMA_HOST_BASE=http://host.docker.internal:11434/v1' >> .env
echo 'IGVF_LLM_MODEL=qwen3.6:35b-a3b-coding-bf16'           >> .env
docker compose up -d agent                    # only the agent, host Ollama supplies the LLM

open http://127.0.0.1:8501                    # browser UI
docker compose run --rm agent kg gene APOE    # one-shot skill from CLI
docker compose down                           # stop (volumes preserved)
```

`./Data` and `./Docs` are mounted into the container so analyses persist
across restarts. Cloud LLM keys (`ANTHROPIC_API_KEY`, `OPENAI_API_KEY`,
`GROQ_API_KEY`, `TOGETHER_API_KEY`, `DEEPINFRA_API_KEY`, `HF_TOKEN`)
plus the IGVF-specific `IGVF_PORTAL_COOKIE` and
`IGVF_ARANGO_PASSWORD` are forwarded from `.env` automatically.

After installing, the `igvfagent` console command gives you four ways to
drive the same skills:

```bash
# 1) Browser UI — chat input, streaming progress, inline plot rendering
pip install 'igvfagent[ui]'
igvfagent ui                      # opens http://127.0.0.1:8501

# 2) Natural-language CLI — same agent, terminal output
igvfagent ask "Give me the comprehensive APOE evidence pack including\
  literature corroboration, single-cell datasets, and FAVOR annotations."

# 3) Direct skills — every tool addressable as a subcommand
igvfagent kg gene APOE --depth 2 --call-singlecell --call-literature
igvfagent splitseq retrieve --limit 50
igvfagent ref design --data-type parse_split_seq
igvfagent portal-kg pull --tissue macrophage --limit 100

# 4) Introspection
igvfagent backends            # registered LLM providers
igvfagent tools               # the tool catalog the agent runtime sees
igvfagent --help              # full skill list
```

The Streamlit UI exposes the **same** ReAct agent as `igvfagent ask`,
with an interactive sidebar for backend / model / tool-subset selection,
a streaming progress trace as the agent plans and calls tools, and inline
rendering of any artefacts the tools produce.

Supported artefact viewers (rendered directly in the chat, no terminal
round-trip):

- **PNG / JPG / SVG / GIF** — inline gallery (UMAP, dot-plot, volcano, IGV-style snapshots).
- **CSV / TSV** — first 400 rows as a sortable `st.dataframe`, plus a download button.
- **PDF** — embedded base64 iframe (FAVOR / cCRE reports, advanced-variant exports).
- **Markdown reports** — rendered in place; backtick-quoted file paths inside the body are auto-followed so the underlying CSV / JSONL / PDF / PNG that the report references each get their own inline viewer.
- **JSONL** — first 50 records normalized into a DataFrame, falls back to raw JSON.
- **JSON** — pretty-printed up to 200KB.
- **HTML** — sandboxed `st.components.v1.html` embed.
- **TXT / LOG** — first 50KB in a code block.
- Anything else — download button.

When the **Claude Code CLI** backend is selected, the model picker is
restricted to the three current Claude 4.x tiers (`claude-opus-4-7`,
`claude-sonnet-4-6`, `claude-haiku-4-5-20251001`) plus a `(custom...)`
escape hatch — picking a retired model id is no longer possible from
the dropdown.

The legacy `python3 Scripts/<skill>.py …` invocations documented later
in this README continue to work unchanged.

## Recreating the hosted container (and why you must)

**The container does not mount the source tree.** `Scripts/streamlit_app.py`
imports `from igvfagent import ...` — the package `pip install` bakes into
`/opt/venv` at image **build** time. So `git pull` on the host updates the
checkout and changes nothing about the running site, however many times you
pull. This has cost real debugging time: an external tester reported a bug
that was already fixed and pushed, because the container was serving an
image built before the fix.

### Check first, then rebuild

```bash
bash Deploy/redeploy.sh --check     # compare running code against the checkout
bash Deploy/redeploy.sh             # pull, rebuild --no-cache, recreate, verify
```

`--check` hashes the top-level modules inside the container and the same
files in the checkout and prints both. The sidebar shows that same hash, so
you can confirm what a browser is talking to without shell access.

One wording caveat: it prints `STALE` for **any** difference, in either
direction. A container carrying newer code than the checkout is also
reported stale.

### It refuses while analyses are running, on purpose

```
REFUSING to recreate: 2 analysis process(es) are running.
```

`--force-recreate` kills running work. The guard checks for
`sc-analyze`, `raw-pipeline`, `crispr-screen`, `gradient-screen`, `bean`,
`sge`, `mct`, `kg-mirror` and `mirror-kg.sh` **before** rebuilding, so a
long job is not destroyed after a ten-minute image build. Override
deliberately with `FORCE_RECREATE=1`, having decided the running work is
expendable:

```bash
FORCE_RECREATE=1 bash Deploy/redeploy.sh
```

`kg-mirror` keeps per-collection state, so interrupting it loses only the
collection in flight, not the collections already mirrored.

### Optional: install BEAN for base-editing screens

`igvfagent bean` implements BEAN's guide-assignment method independently, and
that is the part that fixes the 36.7% → 62.5% read assignment. It does **not**
reimplement BEAN's Bayesian variant/tiling model, and for effect sizes on a
base-editing screen that model is what you want. To have the real `bean`
available to the agent, build with it:

```bash
IGVF_INSTALL_CRISPR_BEAN=1 bash Deploy/redeploy.sh
```

Off by default: it pulls in torch + pyro and roughly doubles the image, and
most deployments do not analyse base-editing screens.

**If a rebuild is not available yet** (a long mirror or analysis in
flight), `bash Deploy/install-bean.sh` installs BEAN into the running
container. It works, and every step in it is a workaround for a different
layer of the container's hardening, each found the hard way:

| blocker | fix |
|---|---|
| apt's http method drops privileges; `cap_drop: ALL` forbids `setgroups` | `-o APT::Sandbox::User=root` |
| apt cannot write `/var/cache/apt/archives/partial` (`_apt`-owned, `CapEff=0`) | redirect `Dir::Cache::archives` to `/tmp` |
| pip as **root** fails: the venv is `igvf`-owned and root has no `CAP_DAC_OVERRIDE` | run pip as the owner, never root |
| no prebuilt wheel exists, and `CRISPResso2Align.pyx` breaks under Cython 3 | `cython<3` **and** `--no-build-isolation` |
| BEAN uses `np.int_t` / `np.Inf`, removed in NumPy 2 (the app runs 2.x) | its own venv with `numpy<2`, re-pinned last because a dependency overrides it |

BEAN goes in `/workspace/opt/bean-venv`, on the data volume, so it survives
container recreation; the gcc install and the `/usr/local/bin/bean` symlink
do not, so re-run the script after a recreate — or build the image with the
flag above and stop needing it.

**It cannot be added afterwards.** The runtime layer has no compiler, and
BEAN's `bean/mapping/CRISPResso2Align.pyx` needs Cython, so
`pip install crispr-bean` inside a running container fails with
`CompileError: bean/mapping/CRISPResso2Align.pyx`. The builder stage already
has `build-essential` for the `[hic]` extra, so that is where it goes — which
means this is a build-time decision, not a runtime one.

BEAN is AGPL-3.0 and IGVFagent is Apache-2.0. This **installs** it as a
separate program invoked as a subprocess — no linking, no vendoring, no
licence propagation. See [Reference GitHub
repositories](#reference-github-repositories).

### Credentials and the shared-deployment settings

`Deploy/make-live.sh` pushes the IGVF Portal key pair and the ArangoDB
credentials into `Deploy/.env.prod` (mode 600) and then redeploys, verifying
inside the container afterwards. Secrets travel on **stdin**, never in argv,
because argv is visible in `ps` on the host for the life of the call.

```bash
bash Deploy/make-live.sh --check    # report only
bash Deploy/make-live.sh            # upsert credentials, redeploy, verify
```

Environment changes need a **recreate**, not a restart. Two settings matter
on a shared deployment and both default to closed:

| Variable | Default | What it allows |
|---|---|---|
| `IGVF_ALLOW_AGENT_AUTHORING` | `0` | the agent writing tools/skills this host then executes |
| `IGVF_ALLOW_UPLOAD_EXTENSIONS` | `0` | a **visitor** uploading a `.py` this host then executes |
| `IGVF_HEAVY_SLOTS` | `1` | concurrent heavy analyses before queueing |
| `IGVF_PUBLIC_MODE` | `0` | fixed backend + curated model allowlist |

If the deployment is reachable by anyone who has the password, leave the
first two at `0`.

### When a hot copy is the right answer instead

Some fixes reach a running container without a rebuild, which matters when a
multi-day job is in flight:

- **`streamlit_app.py`** — Streamlit re-reads its main script on every
  rerun, so a copied file takes effect on the next interaction.
- **Any CLI skill module** — every tool call runs `igvfagent <skill>` as a
  fresh subprocess, so a copied module is picked up immediately.
- **NOT `_agent.py`, `_llm.py` or `_tools.py`** — these are imported into the
  long-lived Streamlit process and cached in `sys.modules`. A new *tool* in
  particular will not appear in the agent's registry until the process
  restarts, even though its CLI works.

```bash
D=/opt/venv/lib/python3.11/site-packages/igvfagent
docker cp Scripts/<module>.py igvfagent-app:$D/<module>.py
docker exec igvfagent-app sh -c "rm -rf $D/__pycache__"
```

A hot copy is a **stopgap**: it lives in the container, not the image, so a
restart reverts it. Follow it with a real `redeploy.sh` at the next window.

## Choosing an LLM backend

IGVFagent runs the same skills regardless of which model drives the
plan-act loop. Pick the backend that matches your goal:

| Goal | Pick this | How |
|---|---|---|
| Fast Streamlit UI, willing to manage one API key | **Anthropic Claude API** | `export ANTHROPIC_API_KEY=…` then sidebar `Backend = anthropic` |
| Already happy with Claude Code, want analysis-from-chat | **Claude Code as the orchestrator** | `cd IGVFagent && source .venv/bin/activate && claude` — ask in the chat |
| Offline / private / free | **Local Ollama (Qwen / Gemma)** | `ollama pull qwen3:8b && export IGVF_LLM_MODEL=qwen3:8b` |
| Best of both worlds | **Cloud for UI, local for batch** | Anthropic for the UI, Ollama for nightly `igvfagent ask` jobs |

Side-by-side trade-offs:

| | Local Ollama (Qwen / Gemma) | Anthropic Claude API | Claude Code CLI |
|---|---|---|---|
| Cost | free | ~$0.001 – $0.01 / query (Haiku → Sonnet) | covered by your Claude Code plan |
| Latency | 30 – 90s (35B bf16) / 8 – 20s (gemma4:31b) / 2 – 5s (qwen3:8b) | 2 – 10s | 5 – 15s (CLI subprocess + auth) |
| Privacy | prompts stay local; **archive queries still leave** (Catalog/ENCODE/GEO) | prompts + tool output go to Anthropic | prompts + tool output go to Anthropic via Claude Code |
| Tool-call quality | strong on Qwen 3.6 35B / Gemma 4 31B; weak on < 7B | best-in-class (native function calling) | best-in-class (Claude Code drives natively) |
| Setup effort | install Ollama + pull a model | one env var, one sidebar click | already installed if you use Claude Code |
| Multi-turn context | full | full | reset per `claude -p` call |
| Best for | offline / private analyses, batch jobs | day-to-day exploratory queries in the UI | mixed coding + analysis sessions |

Backend resolution rules (so the same env / sidebar settings work
everywhere):

1. Sidebar selection (UI) or `--backend` flag (CLI) wins.
2. Else `IGVF_LLM_BACKEND` env var.
3. Else inferred from the model name (`claude-*` → Anthropic;
   `gpt-*` / `o1-*` / `o3-*` → OpenAI; `qwen*` / `gemma*` / `llama*` /
   `mistral*` → Ollama).
4. Otherwise default: **Ollama with Qwen 3 8B**.

Same logic for the model: explicit `--model` > `IGVF_LLM_MODEL` env >
backend's compiled-in default.

### Local LLM (free, private, offline)

The agent's default backend is **Ollama with Qwen 3 8B** — no API key
required. After installing Ollama:

```bash
ollama serve &              # background daemon
ollama pull qwen3:8b
pip install 'igvfagent[llm]'   # adds the openai SDK (Ollama speaks OpenAI's wire format)
igvfagent ask "Walk the IGVF Knowledge Graph for APOE and return all
   regulatory elements plus the matching IGVF single-cell datasets."
```

`igvfagent models` introspects the local Ollama daemon and lists every
installed model + size, so you can pick one for `--model`:

```bash
igvfagent models           # lists installed Ollama models
igvfagent models --json    # machine-readable
```

Bigger / coding-tuned local models work well too — backend inference is
substring-based, so anything with `qwen` / `gemma` / `llama` / `mistral`
/ `phi` / `deepseek` / `yi` / `command` in the tag auto-routes to Ollama:

```bash
# 35B Qwen 3.6 coding tune (e.g. via Ollama):
igvfagent ask --model qwen3.6:35b-a3b-coding-bf16 \
   "Comprehensive APOE evidence pack with literature corroboration."

# 31B Gemma 4 coding tune:
igvfagent ask --model gemma4:31b-coding-mtp-bf16 \
   "Discover Parse SPLiT-seq AnalysisSets profiling macrophages."

# Pin the model globally so you don't repeat --model:
export IGVF_LLM_BACKEND=ollama
export IGVF_LLM_MODEL=qwen3.6:35b-a3b-coding-bf16
igvfagent ask "..."
```

The same `IGVF_LLM_MODEL` env var is forwarded by the Compose stack —
drop it in a local `.env` and the containerized agent picks it up.

For higher-quality answers, point the agent at Anthropic Claude or OpenAI:

```bash
export ANTHROPIC_API_KEY=...
igvfagent ask --backend anthropic --model claude-sonnet-5 \
   "Compare DRD1 and DRD2 striatal MSN evidence in the local KG."
# Other current Claude ids: claude-opus-5, claude-fable-5, claude-haiku-4-5

export OPENAI_API_KEY=...
igvfagent ask --backend openai --model gpt-4o-mini "..."
```

Every run is persisted under `Docs/Agent/<timestamp>_<label>/` (transcript
JSON + markdown report + artefact paths from each tool call).

**Optional dependency groups** (declared in `pyproject.toml`):

| Extra | Adds | When you need it |
|---|---|---|
| `analysis` | pandas, numpy, scipy, matplotlib, pyarrow, anndata, scanpy, seaborn | Single-cell pipelines, advanced variant analysis, plotting |
| `ui` | streamlit | Browser UI (lands in step 4 of the shipping plan) |
| `llm` | anthropic, openai | Native SDKs for the LLM-driven `igvfagent ask` runner (coming) |
| `all` | analysis + ui + llm | Everything |
| `dev` | pytest, ruff, build | Developer tooling |

## Configuration

The agent reads endpoints and credentials from environment variables. Copy
`.env.example` to `.env` and fill in the values you have access to; `.env`
is gitignored.

Public, read-only endpoints have sensible defaults baked into the scripts —
most reads work without any configuration. Authenticated workflows (for
example unreleased datasets, or knowledge-graph queries that require an
account) need credentials you supply locally. Never commit `.env`, cookies,
tokens, or any other authenticated session material.

Set `IGVF_PROJECT_ROOT` if you want to run the scripts from outside the
repository directory.

## Consistent results across LLM backends

IGVFagent is designed to **reduce cross-backend variance**, not to eliminate
it. Pinned command sequences and playbooks are reproducible because no model
decides anything. Free-form agent planning is not: measured artefact agreement
across repeated free-form runs is well below 1.0, and the mechanisms below
narrow the spread rather than remove it. Treat identical output as a property
of pinned workflows, and auditability — every action recorded as a typed,
re-runnable command — as the property that holds in general.

Four mechanisms reduce the variance:

1. **Identical tool set on every backend.** The runtime exposes one canonical,
   deterministically-ordered tool set to *all* backends, capped at
   `IGVF_LLM_MAX_TOOLS` (default 128, matching OpenAI's function limit;
   Anthropic accepts the full catalogue, and the hosted deployment raises it so
   nothing is hidden). **The cap truncates alphabetically**, so leaving it below
   the catalogue size silently removes whole families of tools — previously
   OpenAI-family models were silently trimmed to 128 while others saw more, so
   the same query could pick a tool that didn't exist elsewhere). Override the
   cap with `IGVF_LLM_MAX_TOOLS`.
2. **Deterministic decoding.** Temperature defaults to 0 and a decoding `seed`
   (`IGVF_LLM_SEED`, default 0) is forwarded to every OpenAI-compatible backend
   (Ollama / vLLM / OpenAI …).
3. **Deterministic router.** Unambiguous query shapes bypass LLM tool-choice
   entirely and run a fixed first tool: a bare gene symbol → `kg_gene`, an rsID
   → `kg_variant`, an IGVF/ENCODE accession or URL → `explain_dataset`, a
   `chrN:start-end` region → `kg_region`.
4. **Templated answer synthesis.** The final answer wraps the model's prose in a
   backend-independent skeleton (tools run + arguments + artefacts produced), so
   the substantive content is identical across models; only the narrative varies.

Every run records a **consistency fingerprint** (system-prompt hash + seed +
tool-set hash) in its transcript. Verify consistency yourself:

```bash
igvfagent consistency                         # offline invariants (no API key)
igvfagent consistency --backends anthropic,ollama   # diff tool-call traces across backends
```

Disable the router or templating with `IGVF_ROUTER=0` / `IGVF_TEMPLATED_ANSWER=0`.

## The growing local knowledge graph + database

Every time IGVFagent **downloads data or analyzes/processes raw data**, a little
more of a *local* knowledge graph and database accumulates — so the agent gets
faster and more self-contained the more you use it, and repeat queries can be
answered from local data. This is the **core default mechanism**, wired into
both the agent loop (every tool call) and the CLI (every direct skill run
auto-harvests its new outputs).

- **Knowledge Graph** — `Data/KG/local_kg.sqlite` (`nodes` + `edges`: genes,
  variants, regions, datasets, studies, analyses, and the relations between
  them — the same graph the `portal-kg` skill builds).
- **Database** — the DuckDB warehouse at `Data/Warehouse/igvf.duckdb`.

```bash
igvfagent localstore stats      # nodes / edges / downloads / analyses logged
igvfagent localstore harvest    # scan Docs/ + downloads and ingest anything new
```

Growth is idempotent (deterministic upserts + a harvest ledger) and safe under
concurrent agent/CLI writers (WAL + busy-timeout). Disable with
`IGVF_LOCALSTORE=0`.

## Accounts, approval and sign-in

The hosted deployment has real accounts. Identity comes from the **Genohub
Discourse community** — signup, email verification and moderation already live
there — and access to the agent requires membership of a Discourse group that
administrators control. Signing up and being approved are two separate steps on
purpose: anyone may join the forum; only approved members may drive an agent
that runs analysis pipelines on a shared machine.

Approving someone is adding them to the group. Revoking is removing them.
There is no second password to distribute and no credential store of our own —
authentication sits in the gateway (`Deploy/auth/gate.py`, stdlib only), never
in the app, because the app is the thing being protected.

Setup, cutover and day-to-day administration: [`Docs/AUTH.md`](Docs/AUTH.md).

A local install has no gateway and therefore no login — one person at the
keyboard, nothing hidden, exactly as before.

## Projects and permanent history

The knowledge graph remembers *what is true*. A separate store remembers *what
was done*: `Data/History/history.sqlite` holds one durable row per agent run —
the question, the answer, the artefact paths, the accessions — plus every skill
invocation and download, all in one FTS5 search index.

Two things follow from that.

**A dataset is analysed once.** Ask about an accession that has been worked on
before and the agent looks it up instead of recomputing it:

```bash
igvfagent project recall IGVFDS5414UFNC   # every result ever produced about it
igvfagent project search "spatial ATAC Hi-C"   # full text over questions AND answers
igvfagent project show Docs/Agent/20260914_150552_…   # replay one session
```

The agent has the same three as tools (`history_recall`, `history_search`,
`history_show`) and its first instruction is to use them before re-running work.

**Results can be grouped into a project that outlives the session.**

```bash
igvfagent project create "Wang 2026 reproduction" --use
igvfagent project add --kind dataset --ref IGVFDS5414UFNC --title "primary set"
igvfagent project items
igvfagent project rename "Wang 2026 reproduction" --to "Spatial ATAC-Hi-C repro"
```

While a project is active, every answer is filed into it automatically — in the
CLI and in the web UI, which has a **🗂️ Project** panel in the sidebar.
Renaming is safe: items reference the project's immutable id, and every former
name stays resolvable, so a reference written down months ago still works.

**Results are shared, organisation is private.** Any past answer is readable
by anyone signed in and is recalled automatically when someone asks about the
same accession — paying twice for the same analysis to hide it from a colleague
helps nobody, and each run still records who produced it. Projects are the
private part: visible to their owner and whoever they share them with
(`igvfagent project share <username>`), where members can add but only the
owner can rename, archive or change membership. `IGVF_HISTORY_SHARED=0` flips
answers back to strict per-user privacy. Without authentication in front,
nothing is filtered at all.

**Nothing here is ever deleted.** That is enforced by `BEFORE DELETE` triggers
on every history table, not by convention — removing an item from a project
marks it removed and keeps the row searchable, and archiving a project hides it
from the listing and is reversible. The store is a separate SQLite file from
the knowledge graph on purpose: the KG is bulk-rebuilt by `kg-integrate` and is
disposable scratch, while history must outlive all of it.

Index work that predates the store (safe to re-run; also rebuilds the index):

```bash
igvfagent project backfill
igvfagent project stats
```

## Smoke test

After `cp .env.example .env` and (optionally) editing it, verify the install:

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/igvf_client.py catalog-api /
python3 Scripts/igvf_data_skills.py overview --limit 5
python3 Scripts/igvf_data_skills.py encode-overview --limit 5
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
```

Expected output locations:

- runtime logs in `Docs/Logs/`
- cached API responses in `Data/`
- manifests in `Data/Manifests/`
- reports in `Docs/<skill>/`

## Skill catalog and usage

`Scripts/README.md` has the complete command list. Quick reference below.

### IGVF / ENCODE / Knowledge Graph client

The thin client behind every skill — direct HTTP calls and AQL.

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/igvf_client.py catalog-api /
python3 Scripts/igvf_client.py catalog-files --limit 25
python3 Scripts/igvf_client.py gene TP53 --limit 10
python3 Scripts/igvf_client.py variant rs58658771 --limit 10
python3 Scripts/igvf_client.py encode-search --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_client.py aql "FOR doc IN genes LIMIT 5 RETURN doc"
```

### Catalog / Portal / ENCODE overviews

High-level inventory and smoke summaries.

```bash
python3 Scripts/igvf_data_skills.py catalog-smoke --limit 10
python3 Scripts/igvf_data_skills.py overview --limit 25
python3 Scripts/igvf_data_skills.py encode-overview --limit 25
python3 Scripts/igvf_data_skills.py encode-smoke --limit 10
python3 Scripts/igvf_data_skills.py encode-export-csv --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_frontpage_summary.py refresh --update-readme
```

### Variant annotation

Annotate any user variant CSV against IGVF Catalog evidence (CADD, QTL,
phenotypes, regulatory elements, predictions). Provide your own variant list;
a tiny illustrative example is shipped at
`Data/Input/VariantList/example_variants.csv`.

```bash
python3 Scripts/annotate_variant_list.py \
  --input Data/Input/VariantList/example_variants.csv \
  --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants \
  --input Data/Input/VariantList/example_variants.csv \
  --max-rows 10
```

### Advanced variant analysis

End-to-end pipeline that combines IGVF Catalog evidence, ENCODE cCRE class,
predicted-functional composites, optional user experimental data, logistic
models, and per-gene Miami / volcano / overlap plots into a research report.
See `Docs/Skills/ADVANCED_VARIANT_ANALYSIS_SKILLS.md` for full details.

```bash
# Annotation + composite + plots only (no experimental data)
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/example_variants.csv \
  --label example_locus

# Joining a user CRISPRi/MPRA/GWAS table and modeling an outcome
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/my_variants.csv \
  --experimental Data/Input/Experimental/my_crispri.csv \
  --outcome BEAN_pval_lt_.05 \
  --gene-list LDLR,PCSK9,APOE \
  --label my_crispri_v1

python3 Scripts/advanced_variant_analysis.py write-playbook
```

Outputs land in `Docs/AdvancedVariantAnalysis/` (annotated CSV, summary stats,
logistic model JSON, markdown report) and `Docs/AdvancedVariantAnalysis/Plots/`
(SVG plots).

### Single-cell, multiome, specialized assays

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

### Cross-source multiome survey

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

### Parse SPLiT-seq pipeline

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

### Pathway databases (KEGG + Reactome + WikiPathways, integrated locally)

Pulls the **current** releases, normalises every gene identifier into one
namespace, merges pathways that different databases describe differently, and
loads the result into the local knowledge graph — so pathway structure is
retrieved once and then queried offline. Each stored fact records which
databases assert it, so agreement between them is visible rather than
flattened into a duplicate. The integration method follows IntPath
(Zhou et al. 2012, cited under [References](#references)); the data is pulled
fresh, not taken from that release.

```bash
igvfagent pathwaydb sources                    # databases, licences, coverage
igvfagent pathwaydb pull                       # download current releases (~180 MB)
igvfagent pathwaydb pull --kgml                # + KEGG maps, for typed relations
igvfagent pathwaydb build --relations          # normalise, merge, load into the KG
igvfagent pathwaydb query --genes CLU,BIN1,PICALM,APOE,TREM2
igvfagent pathwaydb evaluate --agreement       # score the unification criteria
igvfagent pathwaydb status                     # what is cached, how old, which release
```

A pull on 2026-08-29 (KEGG 2026/08/27 · WikiPathways 20260810 · Reactome
current) yielded **211,120 gene-pathway memberships over 4,016 pathways and
14,870 genes**, plus **73,925 typed gene-gene relations** (PPrel 54,002 ·
ECrel 15,996 · GErel 3,197 · PCrel 730) — against 582 pathways in the 2012
IntPath release. Downloads are cached under `Data/References/PathwayDB/`
(gitignored) with a per-file release stamp and SHA-256, so a figure can be
traced back to the exact release behind it.

Pathways are unified across databases by IntPath's longest-common-subsequence
rule on pathway names, corroborated by word overlap and gene membership, and
constrained by Reactome's own parent/child hierarchy so curated fact outranks
string similarity. `igvfagent pathwaydb evaluate` scores each criterion on the
current data — the default reaches zero known-wrong merges where the published
rule alone reaches 3.4%. Full method, CLI reference and measurements:
**[`Docs/PATHWAYDB.md`](Docs/PATHWAYDB.md)**.

BioCyc/HumanCyc is not fetched — verified 2026-08-29, its download requires a
subscription "costing at least $5,000", so it cannot be pulled reproducibly.
If you hold a licence, integrate your own export with
`--extra-gmt HumanCyc=<file.gmt>`; it flows through the same normalisation,
unification and KG ingestion as the fetched sources, and nothing licensed is
redistributed.

### scQers — quantitative single-cell enhancer reporters

Two barcodes per construct are the whole trick. A candidate element drives
mCherry carrying an **mBC**; the same construct constitutively expresses an
**oBC**. In single cells the oBC says *which* element a cell received and the
mBC says *how hard it is driving* — which is what makes the readout
quantitative rather than a sort-based enrichment.

```bash
igvfagent scqers extract-bc --in-r1 mBC_R1.fastq.gz --out-file mBC.txt.gz \
    --start 0 --end 15 --check-seq GCT
igvfagent scqers count-bc --in-file mBC.txt.gz --out-file high_mBC.txt \
    --threshold 300 --plot dist.png
igvfagent scqers subassembly --in-file paired_counts.txt --out-file oBC_mBC.tsv
igvfagent scqers pipeline --counts joined_counts.tsv       # activity → specificity → calls
igvfagent scqers selftest                                  # synthetic data, no inputs
```

Activity is bootstrapped against the library's `minP` and `noP` controls **at
matched sample size**, so a rare element is not beaten by a control that was
merely pooled over more cells. Specificity is a permutation of the
cell-to-cluster assignment, because the best of several clusters has a high
fold-change by construction and a fold-change alone cannot tell you otherwise.
p-values are empirical and BH-corrected within each replicate, and an element
is called only if it clears the FDR in *every* replicate it was measured in.

Method: Shendure lab, *Multiplex profiling of developmental enhancers with
quantitative, single-cell expression reporters*
([shendurelab/scQers](https://github.com/shendurelab/scQers), MIT). The
upstream repository is R and shell shared "for transparency" rather than as a
pipeline; this is the method as a runnable CLI. Barcode extraction and
counting reproduce the upstream example output exactly — all 191 barcodes
above threshold, with identical counts.

### Enhancer–gene linkage

```bash
python3 Scripts/enhancer_gene_linkage_skills.py overview --source catalog --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source encode --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py pull-sets \
  --region chr1:903900-904900 --gene SAMD11 --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets \
  --include-local-catalog --demo-if-empty
python3 Scripts/enhancer_gene_linkage_skills.py write-playbook
```

### MPRA / STARR / BlueSTARR

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
python3 Scripts/mpra_data_skills.py analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
python3 Scripts/mpra_data_skills.py write-playbook
```


### MPRA library design → counts → barcode QC

Three clean-room ports of the kircherlab MPRA suite (Max Schubach / BIH, MIT),
covering the whole assay lifecycle. All standard-library-only — no Snakemake,
conda, R, bedtools, pandas or AnnData — and each verified against the upstream
implementation.

```bash
# ── 1. DESIGN the library ────────────────────────────────────────────────
# Tile regions into oligo windows (3 strategies chosen by region length)
igvfagent oligo tile --regions regions.bed.gz --oligo-length 200 --min-overlap 50

# One REF + one ALT oligo per variant x region — the contrast the assay reads
igvfagent oligo design-variants --regions regions.bed.gz \
    --variants variants.vcf.gz --reference hg38.fa --variant-edge-exclusion 20

# Filters: homopolymers, EcoRI/SbfI sites, simple repeats, TSS, CTCF motifs
igvfagent oligo filter --design design.fa --regions regions.bed.gz

# Everything, one call → synthesis-ready design.fa.gz + per-drop filter log
igvfagent oligo pipeline --regions regions.bed.gz --variants variants.vcf.gz \
    --reference hg38.fa --tile --left-adapter AGGACCGGATCAACT \
    --right-adapter CATTGCGTGAACCGA --label my_library

# ── 2. SEQUENCING → per-oligo ACTIVITY ───────────────────────────────────
# Trust the barcode→oligo assignment (min reads + majority agreement)
igvfagent mpraflow assign-filter --pairs pairs.tsv.gz --minimum 3 --fraction 0.75

# Quantify one replicate, with optional barcode outlier removal
igvfagent mpraflow merge-counts --counts rep1.tsv.gz \
    --assignment assignment.tsv.gz --outlier-detection ratio_mad

# "Would deeper sequencing help?" — Lincoln-Petersen mark-recapture
igvfagent mpraflow complexity --barcode-file reporter_experiment.barcode.tsv.gz

# Everything: per-replicate activity → master table → pooled → allelic skew
igvfagent mpraflow pipeline --assignment assignment.tsv.gz \
    --replicate 1=rep1.tsv.gz --replicate 2=rep2.tsv.gz --replicate 3=rep3.tsv.gz \
    --declaration variants.tsv --threshold 10

# ── 3. BARCODE QC + IGVF FORMAT VALIDATION ───────────────────────────────
# Three outlier detectors from Rosen et al.
igvfagent mpralib outliers --barcode-file barcodes.tsv.gz --method global
igvfagent mpralib outliers --barcode-file barcodes.tsv.gz --method large_expression

# Are those calls reproducible across replicates?
igvfagent mpralib consistency --barcode-file barcodes.tsv.gz --method global

# Validate any file against the 8 IGVF MPRA community standards
igvfagent mpralib schemas
igvfagent mpralib validate --file master_table.tsv.gz --schema reporter_experiment
```

**Verification.** `oligo` reproduces all 35 rows of MPRAOligoDesign's golden
tiling fixture. `mpraflow` is **bit-identical** to six upstream scripts run
side-by-side on the same inputs (`merge_label.py` 420/420 values,
`filterAssignmentTsv.py` 238/238 rows, `generateMasterVariantTable.py` 435/435).
Against the published IGVF artefact for AnalysisSet `IGVFDS1933XFAF` —
where the portal hosts both the input *and* the authors' own pipeline output —
it reproduces **210,660 / 210,660 values exactly**, and reproduces the paper's
stated library-complexity figures to the digit. See
[`Benchmarks/rosen2025_mprasnakeflow/`](Benchmarks/rosen2025_mprasnakeflow/README.md).

Two upstream bugs are fixed and documented (a homopolymer run at a sequence's
end was never checked; `use_most_centered_region` measured distance to the
region *end*), each with a flag to restore the original behaviour for exact
reproduction.

### Enhancer / regulatory-region cCRE annotation

Upload a list of **regions** and annotate every interval against the ENCODE
SCREEN candidate cis-regulatory element registry. This is interval overlap —
distinct from `ccre annotate-variants`, which is point-level.

```bash
# Once: download + index the registry (V4 default, 2,348,854 elements)
igvfagent enhancer-annot build-db --registry V4
igvfagent enhancer-annot db-stats

# Which sheet of the workbook actually holds coordinates?
igvfagent enhancer-annot inspect --input Data/Input/EnhancerList/library.xlsx

# Annotate
igvfagent enhancer-annot annotate \
    --input Data/Input/EnhancerList/library.xlsx \
    --sheet High.sgRNA0820 --label my_library
```

Accepts `.xlsx` / `.xlsm` (any sheet), BED, TSV, CSV, gzipped. Coordinate
columns are auto-detected and a headerless BED is recognised by shape; `.xlsx`
is read with `zipfile` + `xml.etree` so no spreadsheet dependency is needed and
million-row sheets stream. A CRISPR library repeats each enhancer once per
sgRNA, so rows are deduplicated on coordinates by default.

Per region you get: number of overlapping cCREs, a ranked representative class
(PLS > pELS > dELS > CA-H3K4me3 > CA-CTCF > CA-TF > CA > TF), every class seen,
bases covered as a **union** (overlapping elements are not double-counted),
coverage fraction, and — when nothing overlaps — the nearest element and its
distance, which is what separates a cCRE desert from a near miss.

Results are written to `Data/Output/EnhancerAnnotation/<timestamp>_<label>/`
as `annotated_regions.csv`, `regions_without_ccre.csv`, `annotated_regions.bed`
and `summary.json`.

> V3 and V4 use **different class vocabularies** (V3: `DNase-H3K4me3` /
> `CTCF-only` / `CTCF-bound`; V4: `CA-H3K4me3` / `CA-CTCF` / `CA-TF` / `CA` /
> `TF`) and V4 has more than twice as many elements, so annotations are not
> comparable across versions. The registry version is recorded in every
> `summary.json`.

### CRISPR screens from raw reads — three shapes, three tools

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

### Saturation genome editing (SGE)

```bash
igvfagent sge design  IGVFDS4629JYPY
igvfagent sge count   IGVFDS4629JYPY
igvfagent sge score   IGVFDS4629JYPY --label palb2
igvfagent sge analyze IGVFDS4629JYPY
```

SGE reads are a fixed amplicon, so `raw-pipeline` refuses them for a stated
reason: quantifying them against a transcriptome would restate the amplicon
design, not measure anything.

### snMCT-seq — RNA and methylation from the same nuclei

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

### CRISPRi / CRISPR-FACS / Perturb-seq

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
python3 Scripts/crispri_data_skills.py write-playbook
```

### SHARE-seq joint scATAC + scRNA QC

Clean-room reimplementation of the QC algorithms in
[broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline)
(MIT, 2021). Consumes IGVF Portal SHARE-seq AnalysisSet deposits
(`fragments.bed[.gz]` + `h5ad`) directly. See
[`Docs/Skills/SHARESEQ_ANALYSIS_SKILLS.md`](Docs/Skills/SHARESEQ_ANALYSIS_SKILLS.md).

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

### STARR-seq allelic test

Clean-room rewrite of
[gaochengwen/STARR-seq-Data-Analysis](https://github.com/gaochengwen/STARR-seq-Data-Analysis)
(no LICENSE — every line is paraphrased from the published `mpra::mpralm`
methods). Implements TPM counts QC + RLE + Spearman D-stat, per-(SNP,
Allele) aggregation, log activity, and the moderated allelic test with
Smyth-2004 trigamma-inversion eBayes + BH-FDR. See
[`Docs/Skills/STARRSEQ_ANALYSIS_SKILLS.md`](Docs/Skills/STARRSEQ_ANALYSIS_SKILLS.md).

```bash
# Discover IGVF Portal STARR-seq MeasurementSets
igvfagent starrseq pull-portal --limit 50 --label survey

# TPM QC: filter low-expression fragments, RLE matrix, outlier samples
igvfagent starrseq qc --input counts.tsv --label run1

# Collapse barcode-level counts to per-(SNP, Allele)
igvfagent starrseq aggregate --input barcode_counts.tsv --label run1

# Per-fragment log activity + per-SNP allelic skew (descriptive)
igvfagent starrseq activity --input aggregated.tsv --label run1

# mpralm-style allelic test (OLS + eBayes + BH-FDR)
igvfagent starrseq allelic-test --input aggregated.tsv --label run1
```

### CRISPRi Flow-FISH screen

Clean-room rewrite of
[EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline)
(MIT, Engreitz Lab 2021), following the published methods in
**Fulco 2019** *Nat Genet* and **Nasser 2021** *Nature*. Per-guide
log-normal MLE on bin counts with EM treatment of an "outside" overflow
bin, real-space conversion, per-element Mann-Whitney + Welch t-test +
BH-FDR. See
[`Docs/Skills/FLOWFISH_CRISPR_SKILLS.md`](Docs/Skills/FLOWFISH_CRISPR_SKILLS.md).

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

### cCRE, FAVOR, IGV-style browser views

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download \
  --only PLS --download --max-download-gb 1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest \
  --source all --limit 100 --hydrate-limit 50
python3 Scripts/ccre_linkage_annotation_skills.py linkage-download \
  --manifest Data/Manifests/cCRELinkage/<manifest.csv> --only rE2G --download
python3 Scripts/ccre_linkage_annotation_skills.py cosmic-from-favor \
  --region chr19:44851820-44908922
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo \
  --region chr19:44850000-44910000
python3 Scripts/ccre_linkage_annotation_skills.py write-playbook
```

Always run a `*-manifest` command before `*-download`. Full SCREEN cCRE,
rE2G, and single-cell linkage corpora can be many gigabytes.

### Start from what the IGVF Portal already computed (`processed lineage`)

Given any IGVF accession, IGVFagent now walks the links the Portal records
around it before it downloads anything. The links are the ones in the Portal's
data-model diagrams:
- `input_for` and `input_file_sets`, through intermediate, principal,
  pseudobulk, model and prediction sets;
- `related_measurement_sets` (the other modality of a multiome);
- `auxiliary_sets` (MULTI-seq, hashing, guides);
- `samples`, then `barcode_map`, then the curated barcode set;
- `seqspecs`, `onlist_files` and `documents`;
- `derived_from`, for a model's training data;
- `large_scale_loci_list`;
- the `QualityMetric` objects the uniform pipeline published.

The walk is directed. Outputs are followed downward and inputs upward, and
datasets from the same sample pool are listed but not expanded. So the "start
here" table only offers files actually built from your accession. Prediction
sets that apply a model trained on your data to other data are listed
separately from results of your data.

```bash
igvfagent processed lineage IGVFDS9875NBZW          # report.md + plan.json + lineage graph + figure
igvfagent processed fetch IGVFDS9875NBZW --max-gb 20   # the processed files + the Portal's QC, not the reads
igvfagent explain explain IGVFDS9875NBZW --download    # now fetches processed results first (--include-raw for reads)
igvfagent processed selftest                           # offline fixture Portal, 28 checks
```

On IGVFDS9875NBZW, a snRNA-seq 10x multiome with MULTI-seq, the walk finds:
- the uniform-pipeline h5ad, 3.7 GB, with its QC (74.9% pseudoaligned, 93.4%
  of reads on the barcode onlist);
- the public ATAC fragments, 5.6 GB, with 29.3% duplicates;
- cell annotations from the principal analysis;
- the MULTI-seq hashing table and the barcode-to-sample map.

That is about 18 GB of processed files, instead of 155 GB of controlled reads
across the RNA and ATAC measurement sets and the MULTI-seq auxiliary set.
Starting from the raw K562 multiome IGVFDS3910GBDQ, it reaches the scE2G
prediction sets IGVFDS5428HHMB and IGVFDS2032HBUP and their model
IGVFDS6290TXZY. Linked objects the credentials cannot see are reported as
"not visible (403)" rather than dropped. `raw-pipeline plan` and
`processed find/plan` now use the same walk.

Agent tools: `portal_lineage` (the agent's first call for any IGVF accession)
and `processed_fetch`.

### Data illustration and interpretation

```bash
python3 Scripts/data_illustration_interpretation.py explain \
  '<igvf-portal-url>/curated-sets/IGVFDS2544COZH/'
python3 Scripts/data_illustration_interpretation.py explain \
  '<encode-portal-url>/search/?type=Annotation&searchTerm=encode-re2g&status!=archived'
python3 Scripts/data_illustration_interpretation.py explain IGVFDS2544COZH \
  --download --max-download-gb 2
python3 Scripts/data_illustration_interpretation.py write-playbook
```

### Reference skill (literature retrieval, validation, design)

The Reference Skill is the literature arm of the agent: it pulls publications
across PubMed / PMC, bioRxiv, medRxiv, arXiv, Semantic Scholar, and OpenAlex,
weights top-tier journals (Nature / Cell / Science families, NEJM, NAR,
Bioinformatics, Genome Biology / Research, eLife, PNAS), and feeds three
distinct functions used by the internal orchestrator's Evaluation Agent.

```bash
# 1) learn — what does the field do for this topic / assay / biosample?
python3 Scripts/reference_skill.py learn \
  --topic '10x multiome human putamen' --limit 30 --top 15

# 2) validate — has anyone seen these genes / variants / regulatory elements
#               before? cross-check IGVFagent discoveries against literature
python3 Scripts/reference_skill.py validate \
  --input Docs/AdvancedVariantAnalysis/<run>/<label>_summary_stats.csv \
  --context 'putamen Parkinson' --limit-per-item 5

# 3) design — given an IGVF data type, recommend a workflow + cognate studies
#             + matching IGVF Portal AnalysisSets
python3 Scripts/reference_skill.py design \
  --data-type parse_split_seq --assay-title 'Parse SPLiT-seq'

# generic multi-source search and playbook
python3 Scripts/reference_skill.py search --query 'enhancer-gene linkage rE2G ABC' --top 20
python3 Scripts/reference_skill.py write-playbook
```

All API responses are cached under `Data/Cache/References/<source>/` with a
14-day TTL; re-querying is free. Reports under
`Docs/References/<timestamp>_<subcommand>_<label>/`.

### Local IGVF KG mirror (Arango → DuckDB)

The full IGVF Catalog Knowledge Graph in Arango is ~2 TB across 58
collections (25 document + 33 edge). The `kg-mirror` skill streams every
collection except the two planet-scale `variants` tables and the two wide-doc
`genes_coding_variants_scores*` tables (which embed ~80 MB per-variant
score matrices that consistently time out the AQL cursor). Streams via
the read-only AQL cursor API, persists each as zstd-compressed Parquet
shards, then registers a DuckDB warehouse with one view per collection. Lets `igvfagent kg ...` and the
downstream skills run offline against the cached copy. See
[`Docs/Skills/KG_MIRROR_SKILL.md`](Docs/Skills/KG_MIRROR_SKILL.md).

```bash
# Requires Arango read credentials in .env:
#   IGVF_ARANGO_USER=guest
#   IGVF_ARANGO_PASSWORD=guestigvfcatalog

# 1. Inventory the upstream KG
igvfagent kg-mirror inventory

# 2. Mirror everything except the two planet-scale variants tables (default)
igvfagent kg-mirror pull-all

# 3. Or pull one collection at a time (resumable — state on disk)
igvfagent kg-mirror pull --collection genes
igvfagent kg-mirror pull --collection coding_variants --batch-size 10000

# 4. Cap to small/medium collections only (≤ 10 GB each)
igvfagent kg-mirror pull-all --max-collection-bytes 10000000000

# 5. Register Parquet shards as DuckDB views
igvfagent kg-mirror register
igvfagent kg-mirror verify
```

Storage layout:
- `Data/Warehouse/KG/<collection>/<NNNN>.parquet` — ZSTD-compressed shards
- `Data/Warehouse/KG/_state/<collection>.json` — resume cursor
- `Data/Warehouse/igvf_kg_mirror.duckdb` — one `kg_<collection>` view each

### Knowledge Graph traversal

The orchestrator-friendly **comprehensive context** tool. Starts from a
single entity (gene, variant, or genomic region) and iteratively walks the
IGVF Catalog Knowledge Graph (`api.catalogkg.igvf.org` + the underlying
ArangoDB at `db.catalog.igvf.org`), assembling a unified evidence pack that
includes direct neighbors, second-degree relations, and optional cross-skill
enrichment (FAVOR, enhancer-gene linkage, IGVF single-cell datasets, prior
literature). One CLI call → per-relation manifests + JSON evidence pack +
markdown report.

```bash
# Comprehensive APOE traversal — variants, cCREs, transcripts, proteins,
# diseases, pathways, plus enhancer-gene linkage, candidate single-cell
# datasets (matched via the gene's eQTL biological contexts), and prior
# literature
python3 Scripts/kg_traversal_skill.py gene APOE \
  --depth 2 --limit 50 --max-variants 25 --subvariant-limit 10 \
  --call-favor --call-linkage --call-singlecell --call-literature \
  --literature-context Alzheimer cardiovascular --label apoe_full

# Variant-centric: one rsID -> linked genes, cCREs, phenotypes, predictions
python3 Scripts/kg_traversal_skill.py variant rs429358 \
  --call-favor --call-literature --label apoe_e4_variant

# Region-centric: genes + cCREs + linkage in window
python3 Scripts/kg_traversal_skill.py region chr19:44903000-44912000 \
  --call-favor --label apoe_locus

# Direct AQL pass-through
python3 Scripts/kg_traversal_skill.py aql \
  'FOR g IN genes FILTER g.name == "APOE" RETURN g'
python3 Scripts/kg_traversal_skill.py write-playbook
```

Outputs land under `Docs/KGTraversal/<timestamp>_<label>/` — a markdown
report with sectioned evidence per relation type, a complete
`evidence_pack.json`, and one CSV per relation under `Manifests/`. The
manifest CSVs feed directly into the variant-annotation, advanced
variant-analysis, single-cell, SPLiT-seq, and reference skills, which is
exactly the cross-skill composition that the internal Plan → Action →
Results → Evaluation orchestrator chains together.

### Portal → local KG ETL

A persistence and indexing layer that ingests unstructured IGVF Portal
entities (AnalysisSets, MeasurementSets, Samples, Donors, Files,
FileSets, Documents) into a **local SQLite-backed knowledge graph** that
mirrors the IGVF Catalog ArangoDB schema (nodes + edges + properties +
provenance). The local KG can be queried, annotated, enriched from the
remote Catalog, and later exported as `arangoimport`-compatible JSONL for
a future push to the IGVF Catalog. Endpoint URLs and credentials are
resolved through `_endpoints.py` env-var overrides — nothing sensitive is
written to source.

```bash
# 1. Pull a tissue / assay / lab corpus from the Portal (one hop expansion)
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --tissue macrophage --limit 100 --depth 1
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --assay 'Parse SPLiT-seq' --limit 100 --depth 1
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --gene APOE --limit 25 --depth 1

# 2. Text-mine for gene + variant mentions, confirm against the Catalog
python3 Scripts/portal_to_kg_skill.py annotate

# 3. Hydrate confirmed genes from the Catalog (variants, cCREs, diseases…)
python3 Scripts/portal_to_kg_skill.py enrich --symbols APOE,TREM2,LDLR --limit 25

# 4. Read the local KG
python3 Scripts/portal_to_kg_skill.py stats
python3 Scripts/portal_to_kg_skill.py query --gene APOE --limit 50
python3 Scripts/portal_to_kg_skill.py query --node-id analysis_sets/IGVFDS3222WCZH

# 5. Export to ArangoDB-compatible JSONL for an eventual Catalog push
python3 Scripts/portal_to_kg_skill.py export-aql
python3 Scripts/portal_to_kg_skill.py export-cytoscape --limit 500
python3 Scripts/portal_to_kg_skill.py write-playbook
```

Local store: `Data/KG/local_kg.sqlite` (gitignored under `Data/*`).
Exports: `Data/KG/Export/<timestamp>/nodes_<collection>.jsonl` plus
`edges.jsonl`, ready for `arangoimport --type jsonl …`.

### ENCODE bulk-genomics pipeline

End-to-end retrieval, description, peak QC, super-enhancer calling,
SCREEN cCRE integration, and IGV-style browser-SVG visualization for
the major ENCODE bulk assays — ChIP-seq (TF + Histone), ATAC-seq,
DNase-seq, Hi-C, capture Hi-C, ChIA-PET, plus RNA-seq / MNase-seq /
FAIRE-seq / CAGE / RAMPAGE for retrieval & description.

```bash
# 1) Discover experiments by assay × biosample × target
igvfagent encode retrieve --assay 'Histone ChIP-seq' \
  --target H3K27ac --biosample K562 --assembly GRCh38 --limit 50

# 2) Per-file inventory for one or more accessions
igvfagent encode manifest --accessions ENCSR000AKP --label k562_h3k27ac

# 3) Pull files under a size cap, filter by file format
igvfagent encode download \
  --manifest Data/Manifests/ENCODE/<files.csv> \
  --max-download-gb 5 --formats bed bigWig

# 4) Plain-language description of one experiment
igvfagent encode describe --accession ENCSR000AKP

# 5) Peak QC (count, width, score, per-chromosome) from a BED file
igvfagent encode analyze-peaks --bed peaks.bed.gz

# 6) ROSE-style super-enhancer call (H3K27ac / BRD4 / MED1 / P300)
igvfagent encode super-enhancers --bed h3k27ac_peaks.bed \
  --stitching-distance 12500 --tss-bed tss.bed --tss-distance 2000

# 7) Overlay peaks with SCREEN cCRE classes (PLS / pELS / dELS / CTCF)
igvfagent encode integrate-ccre --bed peaks.bed

# 8) IGV-style multi-track SVG browser view
igvfagent encode browser --region chr19:44903000-44912000 \
  --track 'H3K27ac peaks:peaks.bed' \
  --track 'ATAC peaks:atac_peaks.bed' \
  --with-ccre --label apoe_locus
```

#### Everything a portal holds for one cell line

"Summarise the ENCODE and IGVF data for GM12878, with tables and plots" is
one command, not a search URL:

```bash
igvfagent biosample-census --biosample GM12878          # both portals
igvfagent biosample-census --biosample K562 --portal encode --status all
```

It counts every object type on both portals through their **structured
sample-term filters** (`biosample_ontology.term_name` on ENCODE,
`samples.sample_terms.term_name` on IGVF), resolves the term case-insensitively
through the shared ontology id (GM12878 is `EFO:0002784` on both), and writes
`Docs/BiosampleCensus/<run>/` with `report.md`, `totals_by_type.csv`,
`facet_counts.csv`, per-type item tables, and SVG (+PNG) figures for assays,
ChIP targets, labs, annotation types, file formats and release years. Free-text
search is deliberately not used for counting: `searchTerm=GM12878` on the IGVF
Portal matches thousands of unrelated sets. A zero-hit search, which both
portals answer with HTTP 404, is recorded as 0. The agent tool is
`biosample_portal_census`.

`explain_dataset` on a search URL that matches nothing now says *why* --
which filter field the portal does not know (legacy names such as
`biosample_term_name` are rewritten automatically) or which value is
misspelt (`gm12878` vs `GM12878`) -- instead of writing an empty report.

The agent runtime exposes `encode_retrieve`, `encode_describe`,
`encode_super_enhancers`, `encode_integrate_ccre`, and `encode_browser`
as tools, so a single `igvfagent ask` can drive the full pipeline:

```bash
igvfagent ask "Find K562 H3K27ac ChIP-seq experiments on GRCh38, pick \
  the highest-quality one, describe it, and call super-enhancers from \
  its peak BED. Then overlay the super-enhancers against SCREEN cCREs \
  and produce a browser view of the APOE locus."
```

### Super-enhancer → target-gene pipeline

End-to-end discovery of H3K27ac/BRD4/MED1/P300 ChIP-seq experiments,
ROSE-style super-enhancer calling, ranked-signal "hockey-stick" plots,
multi-track browser views of top SEs, and SE → target-gene linkage via
four parallel streams: Hi-C / ChIA-PET loops, IGVF Catalog rE2G
predictions, TSS proximity windows, and constituent cCRE composition.
Optional motif enrichment (CTCF, AP-1, GATA1, ETS, NFkB, STAT1, FOXA1,
TP53, MYC, SP1) and SE↔cCRE-density correlation. Playbook:
[`Docs/Skills/SE_TARGET_PIPELINE_SKILLS.md`](Docs/Skills/SE_TARGET_PIPELINE_SKILLS.md).

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

### GEO retrieval

Search NCBI GEO Series, parse SOFT metadata, list FTP file inventories,
download supplementary files, and produce a tidy sample sheet for any
downstream RNA-seq / ATAC / ChIP analysis. Pure stdlib + `requests` — no
Bioconductor dependency. Playbook:
[`Docs/Skills/GEO_RETRIEVAL_SKILLS.md`](Docs/Skills/GEO_RETRIEVAL_SKILLS.md).

```bash
# 1) Keyword search → top GSEs (title, organism, summary, n samples)
igvfagent geo search --query "GM12878 RNA-seq" --max-results 10

# 2) Pull the SOFT family file and parse series + sample metadata
igvfagent geo series --gse GSE9574 --label nci_breast

# 3) List supplementary files on the GEO FTP site for one Series
igvfagent geo list-files --gse GSE9574

# 4) Download chosen supplementary files (size-capped)
igvfagent geo download --gse GSE9574 --max-download-gb 5 \
  --extensions .txt.gz .tsv.gz .csv.gz

# 5) Write a tidy sample sheet (one row per GSM, condition columns)
igvfagent geo sample-sheet --gse GSE9574 --label nci_breast
```

Outputs land in `Docs/GEO/<timestamp>_<label>/` with the SOFT-derived
metadata JSON, file inventory CSV, sample sheet CSV, and a
human-readable summary. The agent runtime exposes `geo_search`,
`geo_series`, and `geo_download` as tools.

### Bulk RNA-seq analysis

Counts QC, sample PCA + correlation heatmap, differential expression
(pyDESeq2 when available, Welch's t-test on log-CPM with BH FDR
fallback), volcano + MA + top-DEG z-scored heatmap, and DEG → cCRE
linkage via the IGVF Catalog `/api/genes/genomic-elements` endpoint.
Playbook:
[`Docs/Skills/RNASEQ_ANALYSIS_SKILLS.md`](Docs/Skills/RNASEQ_ANALYSIS_SKILLS.md).

```bash
# 1) Counts QC (library size, gene detection, mito %, top genes)
igvfagent rnaseq qc --counts counts.tsv --label k562_vs_gm12878

# 2) PCA + sample correlation heatmap
igvfagent rnaseq pca --counts counts.tsv --sample-sheet samples.csv \
  --label k562_vs_gm12878

# 3) Differential expression (uses pyDESeq2 if installed, else Welch + BH)
igvfagent rnaseq deg --counts counts.tsv --sample-sheet samples.csv \
  --condition-col condition --treated K562 --control GM12878 \
  --label k562_vs_gm12878

# 4) Link significant DEGs to controlling cCREs via IGVF Catalog
igvfagent rnaseq link-cre \
  --deg Docs/RNAseq/<run>/k562_vs_gm12878_deg.csv \
  --label k562_vs_gm12878 --padj-threshold 0.05 --max-genes 200

# 5) End-to-end QC → PCA → DEG → cCRE linkage in one call
igvfagent rnaseq pipeline --counts counts.tsv --sample-sheet samples.csv \
  --condition-col condition --treated K562 --control GM12878 \
  --label k562_vs_gm12878
```

Outputs include `<label>_deg.csv`, `<label>_deg_to_cre.csv`, the
volcano/MA/heatmap/PCA PNGs, and a markdown report under
`Docs/RNAseq/<timestamp>_<label>/`. Chains naturally with `geo` (pull
counts) and `se-targets` (cross-reference up-regulated genes against
SE-driven targets).

### Proteomics & PPI knowledge graph

End-to-end protein/PPI skill: pulls and version-tracks
**BioGRID**, **IntAct**, **HuRI**, **Reactome**, **KEGG**, **UniProt
id-mapping**, and **all 214 IGVF Portal protein-slim MeasurementSets**
(plus the 14 PPI-score files, the DUAL-IPA fluorescence file, and the
UniProt / protein-language-model reference files). Integrates everything
into a local SQLite knowledge graph (`Data/Proteomics/KG/proteomics.sqlite`)
with deduplicated edges keyed on `(id_a, id_b, source, source_id)`.
Provides summary stats, network visualizations, per-IGVF-assay example
figures from real Portal files, and a literature survey for the IGVF
protein assays via the Reference skill. Playbook:
[`Docs/Skills/PROTEOMICS_SKILLS.md`](Docs/Skills/PROTEOMICS_SKILLS.md).

```bash
# 1) Download (only sources that changed upstream)
igvfagent proteomics download --source all
# Or one at a time:
igvfagent proteomics download --source biogrid
igvfagent proteomics download --source intact
igvfagent proteomics download --source huri
igvfagent proteomics download --source reactome
igvfagent proteomics download --source kegg --kegg-max-pathways 400

# 2) Version manifest (local + upstream probe)
igvfagent proteomics versions

# 3) Pull all IGVF Portal protein assays + actual files (semi-qY2H,
#    DUAL-IPA, VAMP-seq, MAVE) from the public S3 bucket
igvfagent proteomics igvf-protein

# 4) Build the integrated PPI-KG (SQLite)
igvfagent proteomics build-kg --sources all

# 5) Summary statistics (per source / evidence type / detection method,
#    top-30 hubs)
igvfagent proteomics kg-stats --label initial

# 6) Network visualizations (degree distribution, top-30 hubs,
#    per-source breakdown, ego graph for a query gene)
igvfagent proteomics kg-visualize --gene TP53 --label tp53

# 7) Literature survey for IGVF protein assays — restricted to the
#    Nature/Cell/Science journal family
igvfagent proteomics assay-survey --label may2026

# 8) Per-assay example histograms generated from real IGVF Portal files
#    (semi-qY2H v1/v2/v3, DUAL-IPA, plus VAMP-seq family / MAVE if
#    those Portal files are pulled)
igvfagent proteomics assay-figures --label demos

# 9) End-to-end orchestrator
igvfagent proteomics pipeline --label may2026 --gene TP53 \
  --sources biogrid,intact,huri,reactome,kegg,igvf

# 10) VAMP-seq deep analysis — pull canonical MaveDB scoresets, run the
#     full Matreyek/Suiter/Clausen/Coyote-Maestas analysis pipeline,
#     and inventory the IGVF Portal raw VAMP-seq experiments
igvfagent proteomics vampseq-pull              # PTEN, TPMT, VKOR, PRKN,
                                               # CYP2C9, NUDT15 from MaveDB
igvfagent proteomics vampseq-analyze --gene PTEN --label pten_deep
igvfagent proteomics vampseq-analyze           # all 6 catalogued targets
igvfagent proteomics vampseq-inventory --label igvf_f9
```

The skill maintains a `Data/Proteomics/_versions.json` manifest with
URL, sha256, and record count per source. `update` only re-fetches
sources where the upstream version differs from the local one. KEGG
calls are throttled (default 0.4 s/req) to respect their TOS; IntAct
defaults to the smaller `intact-micluster.txt` rather than the full
~700 MB `intact.zip`. HuRI uses Ensembl gene IDs; run
`proteomics download --source uniprot` to populate the `id_map` table
for UniProt ↔ Ensembl ↔ Symbol harmonization.

The `vampseq-analyze` subcommand follows the canonical pipeline
distilled from Matreyek *Nat Genet* 2018, Suiter *eLife* 2020, Clausen
*Nat Commun* 2024, and Coyote-Maestas *Nat Commun* 2024 (MultiSTEP) —
producing six publication-grade plots per gene: score distribution,
residue × AA heatmap (the iconic VAMP-seq view), per-residue mean ±
IQR with a domain track, replicate concordance scatter with Pearson
*r*, abundance-class breakdown, and a cumulative ranked-variant curve.
Domain tracks are pre-curated for PTEN (PIP4-bind / Phosphatase / C2 /
C-tail), TPMT, VKOR, PRKN (Ubl / Linker / RING0 / RING1 / IBR / RING2),
CYP2C9, and NUDT15.

The `vampseq-inventory` subcommand decodes the alias scheme on the
IGVF Portal MeasurementSets (`<lab>:<GENE>-DMS-<antibody>-Tile<i>-Replicate<j>-Bin<k>`)
into per-gene coverage matrices: the 144 MultiSTEP sets resolve to
**F9 (Coagulation Factor IX)** across 3 tiles × 4 bins × 4 replicates ×
5 antibody readouts (Light-chain, Heavy-chain, Strep-II-tag, and two
carboxylation-sensitive Gla-domain antibodies); the 36 plain VAMP-seq
sets cover **CYP2C19** and **G6PD**.

### Single-cell analysis (UMAP / t-SNE / Leiden / markers)

End-to-end Scanpy-driven single-cell workflow: **QC → normalize → HVG
→ PCA → k-NN → UMAP → t-SNE → Leiden clustering → marker-gene DE →
publication figures.** Closes the gap where IGVFagent could discover
and download counts matrices but had to hand off to "use Scanpy or
Seurat" for the actual visualization. Auto-detects input format
(`.h5ad`, 10x `.h5`, 10x `.mtx`, CSV/TSV). Playbook:
[`Docs/Skills/SINGLECELL_ANALYSIS_SKILLS.md`](Docs/Skills/SINGLECELL_ANALYSIS_SKILLS.md).

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

### Perturbation Catalogue retrieval

Pulls metadata and per-row perturbation effects from the public
**Perturbation Catalogue** Search API, which indexes ~1,222 datasets
across MAVE (DMS / VAMP-seq family), CRISPR screens (DepMap, Project
Score, Project Achilles, etc.), and Perturb-seq (Replogle 2022,
Nadig 2025, X-Atlas/Orion 2025). Playbook:
[`Docs/Skills/PERTURBATION_CATALOG_SKILLS.md`](Docs/Skills/PERTURBATION_CATALOG_SKILLS.md).

```bash
# 1) Catalogue-wide summary
igvfagent perturb-catalog summary

# 2) Global gene search (one row per perturbed gene, all modalities)
igvfagent perturb-catalog search --query BRCA1 --size 10

# 3) Modality-scoped search with filters
igvfagent perturb-catalog search-modality --modality mave \
    --perturbation-gene-name TP53 \
    --perturbation-position 100_300 \
    --effect-score-name vamp_score \
    --effect-score-value 0.5_1.0 \
    --dataset-limit 25
igvfagent perturb-catalog search-modality --modality crispr-screen \
    --perturbation-gene-name BRCA1 --dataset-limit 25
igvfagent perturb-catalog search-modality --modality perturb-seq \
    --query "lung cancer" --dataset-limit 10

# 4) Full dataset record + paginated per-perturbation rows
igvfagent perturb-catalog dataset --dataset-id replogle_2022_k562_essential_normalized
igvfagent perturb-catalog dataset-rows --modality perturb-seq \
    --dataset-id replogle_2022_k562_essential_normalized \
    --limit 500 --offset 0

# 5) Perturb-seq GSEA hallmark/pathway table
igvfagent perturb-catalog gsea --query BRCA1 --size 50

# 6) Bulk dataset download (auto-detects .csv.gz / .json / .zip)
igvfagent perturb-catalog download --modality perturb-seq \
    --dataset-id replogle_2022_k562_essential_normalized

# 7) End-to-end gene-centric pipeline (the headline command for
#    "what perturbation data exists for gene X?")
igvfagent perturb-catalog pipeline --gene BRCA1 --dataset-limit 10
```

The agent runtime registers `perturb_catalog_summary`,
`perturb_catalog_search`, `perturb_catalog_search_modality`,
`perturb_catalog_dataset`, `perturb_catalog_dataset_rows`,
`perturb_catalog_gsea`, and `perturb_catalog_pipeline` as tools.
Live smoke test on `--gene BRCA1` returns **1,191 CRISPR-screen
datasets** (509 significant), **7 Perturb-seq datasets**, and the
canonical BRCA1 GSEA signature (`HALLMARK_E2F_TARGETS`,
`HALLMARK_G2M_CHECKPOINT`, `HALLMARK_MYC_TARGETS_V1`). Outputs land
under `Docs/Perturbation/<timestamp>_<gene>/` with one JSON per
modality plus a markdown report; bulk downloads go to
`Data/Perturbation/Downloads/<modality>/`.

### MULTI-seq / Cell Hashing demultiplexing

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
[`Docs/Skills/MULTISEQ_ANALYSIS_SKILLS.md`](Docs/Skills/MULTISEQ_ANALYSIS_SKILLS.md).

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

### Integrated data warehouse (DuckDB Silver tier)

The central data warehouse — the "Silver tier" of the integrated
data layer described in
[`Docs/Architecture/INTEGRATED_DATA_LAYER.md`](Docs/Architecture/INTEGRATED_DATA_LAYER.md).
Every other IGVFagent skill becomes a producer that lands QC'd rows
into canonical entity / edge / measurement tables in a single
**DuckDB** database at `Data/Warehouse/igvf.duckdb`, ready to feed
downstream embedding extraction + foundation-model training.
Playbook:
[`Docs/Skills/WAREHOUSE_SKILLS.md`](Docs/Skills/WAREHOUSE_SKILLS.md).

```bash
# 1) Initialise the schema (entities + edges + measurements + provenance)
igvfagent warehouse init

# 2) Pull every available producer into the warehouse
igvfagent warehouse ingest --source all
# or one at a time:
igvfagent warehouse ingest --source proteomics-kg
igvfagent warehouse ingest --source perturb-catalog
igvfagent warehouse ingest --source multiseq
igvfagent warehouse ingest --source mavedb-vampseq

# 3) Stats and cross-skill SQL queries
igvfagent warehouse stats
igvfagent warehouse query \
    "SELECT gene_id, COUNT(*) AS n, AVG(score) FROM vampseq_scores GROUP BY gene_id"
```

DuckDB is chosen over SQLite / Postgres / BigQuery because it is
embedded (no server), columnar (10–100× faster on the analytical
joins this workload runs), reads / writes Parquet natively (so the
Gold-tier embedding tables stay first-class), and ships with vector
functions. Live smoke test (laptop, ~6 s): **32K PPI edges**,
**102K MULTI-seq cells**, **24K VAMP-seq scores**, **28 IGVF
protein-evidence datasets** ingested from five producers. Bulk
`posteriors.csv` etc. are auto-excluded; the tracked outputs are
intentionally lightweight (~14 MB on disk).

### Network integration (clean-room MILP — CARNIVAL + Steiner)

The **integration layer** of the IGVF data warehouse. Implements two
context-specific-subnetwork methods from scratch in cvxpy:

  * **CARNIVAL** — given a signed prior-knowledge graph plus signed
    perturbations + signed measurements, find the minimum-cost
    upstream subnetwork whose vertex signs match the data.
  * **Prize-collecting Steiner** — given per-gene prizes and a PPI,
    find the connected subgraph maximising (prizes − edge costs).

The math is the **CORNETO formulation** ([Rodriguez-Mier et al., *Nat
Mach Intell* 2025](https://www.nature.com/articles/s42256-025-01069-9))
re-implemented in original Apache-2 cvxpy — **no CORNETO source is
imported or vendored.** Algorithms are not copyrightable; source code
is. The full algorithm specification, attribution, and citation list
live in
[`Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`](Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md).
Playbook:
[`Docs/Skills/NETWORK_INTEGRATION_SKILLS.md`](Docs/Skills/NETWORK_INTEGRATION_SKILLS.md).

```bash
pip install 'igvfagent[network]'    # cvxpy + SCIP (free MILP)

# 1) End-to-end self-test: synthetic EGFR → SOS1 → … → MYC cascade
igvfagent network demo --beta 0.05 --solver SCIP
# → CARNIVAL recovers all 6 cascade edges and writes them to the
#   warehouse with upstream='network:demo'.

# 2) Materialise a signed PKN from the proteomics KG
igvfagent network pkn-from-kg --label reactome_pkn

# 3) CARNIVAL on Perturb-seq-style inputs
#    perts.csv: gene,sign   (sign ∈ {-1,+1})
#    degs.csv:  gene,score  (signed log2 fold-change)
igvfagent network carnival \
    --perturbations perts.csv \
    --measurements  degs.csv \
    --pkn-limit 5000 --solver SCIP --label perturb_seq_demo

# 4) Prize-collecting Steiner tree (VAMP-seq abundance prizes on PPI)
igvfagent network steiner --terminals vamp_prizes.csv \
    --pkn-limit 5000 --edge-cost 1.0
```

License boundary: **Apache-2 throughout**. No GPL runtime dependencies.

### Tabula Sapiens 2.0 (reference human cell atlas)

Retrieval and figure-level reproduction of the **Tabula Sapiens 2.0**
human cell atlas (Tabula Sapiens Consortium, *Cell* 2026): 1,136,218
cells across 28 tissues from 24 donors, annotated to 182 fine and 38
broad cell types. Clean-room reimplementation of the analyses in
[czbiohub-sf/tabula-sapiens](https://github.com/czbiohub-sf/tabula-sapiens)
(BSD-3-Clause). Playbook:
[`Docs/Skills/TABULA_SAPIENS_SKILL.md`](Docs/Skills/TABULA_SAPIENS_SKILL.md).

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

### Human Transcription Factors database

Local SQLite mirror of the Lambert/Jolma/Hughes **Human Transcription
Factors** database v1.01 ([humantfs.ccbr.utoronto.ca](https://humantfs.ccbr.utoronto.ca);
Lambert et al., *Cell* 2018): 2,765 assessed proteins of which **1,639
are curated TFs**, with DNA-binding-domain family, binding mode, motif
status and cross-references. Playbook:
[`Docs/Skills/HUMANTFS_SKILL.md`](Docs/Skills/HUMANTFS_SKILL.md).

```bash
igvfagent humantfs build-db                     # ~2 MB, seconds
igvfagent humantfs is-tf --genes "FOXP3,EOMES,ACTB,GAPDH"
igvfagent humantfs lookup --genes SATB2
igvfagent humantfs list --dbd "C2H2 ZF" --with-motif
igvfagent humantfs families
```

`is-tf` distinguishes *unassessed* from *not a TF* — absence from the
database is not evidence against. The data is downloaded at build time,
never vendored; cite Lambert 2018 for derived tables.

### scE2G workbench: train, check and benchmark E2G models (`sce2g`)

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

### TF Perturb-seq: calibrated effects to disease and GWAS (`tf-perturb`)

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

### eQTL enrichment benchmark for enhancer-gene predictors (`eqtl-enrich`)

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

### IGVF CRISPR Pipeline, early Nextflow version (`igvf-crispr-pipeline`)

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

### SCEPTRE for IGVF MuData (`sceptre-igvf`)

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

### CRISPRi-FlowFISH pipeline (`flowfish-pipeline`)

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

### Pooled CRISPR screen toolkit (`perturb-tools`)

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

### IGVF single-cell pipeline (`igvf-sc-pipeline`)

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

### GWAS E2G benchmark (`gwas-e2g`)

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

### ENCODE-rE2G (`encode-re2g`)

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

### Perturb-seq pipeline, bulk_crispr_pipeline port (`bulk-crispr`)

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

### CRISPR module dev tools (`crisprdevtools`)

A port of [IGVF-CRISPR/crisprdevtools](https://github.com/IGVF-CRISPR/crisprdevtools), the helper used to start new modules of IGVF_CRISPR_Pipeline. Upstream has one function, `create_new_module_nextflow(name)`. `new-module` recreates its layout and file contents exactly: `bin/`, `conda_envs/`, `example_data/`, `processes/`, `test/`, `README.md`, `input.config`, a DSL2 `main.nf`, `bin/<name>.py`, `conda_envs/<name>.yaml` and `processes/<name>.nf`. Three upstream bugs are fixed by default and come back with `--upstream-compat`: the conda YAML is indented so it does not parse, every process is named `seqSpecParser`, and the bin script is not executable.

`--template full` also writes a working include plus workflow, an `input.config` params block, `test/test.nf`, example data and an argparse script skeleton. `check-module` validates any module directory: layout, DSL2 header, processes and includes, conda env files, bin shebangs and executable bits, and the params `main.nf` uses.

```bash
igvfagent crisprdevtools new-module --name guide_assignment --dest Modules
igvfagent crisprdevtools new-module --name guide_assignment --template full --dest Modules --force
igvfagent crisprdevtools check-module --path Modules/seqSpecParser
igvfagent crisprdevtools selftest
```

Agent tools: `crisprdevtools_new_module`, `crisprdevtools_check_module`.

### CRISPR-FG scPerturb-seq jamboree (`crispr-fg-jamboree`)

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

### Fishash Table 2 reproduction (`fishash-table2`)

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

### CRISPR SeqSpec (`crispr-seqspec`)

Port of [IGVF-CRISPR/CRISPR-SeqSpec](https://github.com/IGVF-CRISPR/CRISPR-SeqSpec), the CRISPR focus group's seqspec descriptions of CROP-seq + 10x v3 + MULTI-seq (one spec per modality) and TAP-seq (Schraivogel 2020). The skill carries a compact index of the four specs at the pinned commit and re-implements, without seqspec or PyYAML, what `seqspec check`, `seqspec index -t kb|chromap`, `seqspec onlist` and `seqspec print` do, for both the 0.2.x and legacy 0.0.x schemas.

It also checks submitted FASTQs against a spec (read lengths, fixed sequences at their coordinates, onlist hit rates). On the upstream files this finds that the guide FASTQs are 150 bp while guide.yml declares 26/43 bp reads, that the 26-bp R1 ends inside the 12-bp UMI, and that 88-95% of cell barcodes are on the whitelist.

```bash
igvfagent crispr-seqspec catalog
igvfagent crispr-seqspec check --spec guide
igvfagent crispr-seqspec index --spec tapseq --modality crispr --tool kb
igvfagent crispr-seqspec fetch && igvfagent crispr-seqspec check-reads --spec guide
```

Agent tools: crispr_seqspec_catalog, crispr_seqspec_check, crispr_seqspec_index, crispr_seqspec_onlist, crispr_seqspec_info, crispr_seqspec_check_reads, crispr_seqspec_fetch.

### Gasperini 2019 pipeline (`gasperini-pipeline`)

Port of [IGVF-CRISPR/Pipeline_Gasperini_2019](https://github.com/IGVF-CRISPR/Pipeline_Gasperini_2019), the IGVF processing of the Gasperini et al. 2019 CRISPRi pilot with the single-cell-like Nextflow pipeline (LucasSilvaFerreira/pipeline_perturbseq_like). Every stage is re-implemented in Python with the config's thresholds: guide table, read composition, kb commands, per-lane cell QC with doublets, MuData creation, guide assignment at UMI > 3, covariates, 1-Mb cis pairs, a SCEPTRE-style conditional randomisation test (or R SCEPTRE when available), BH/Fisher aggregation into mudata_results.h5mu, MULTI-seq demultiplexing, the original Gasperini NB likelihood-ratio test, and a comparison with GEO GSE120861. Six upstream quirks are corrected by default and reproduced with `--upstream-compat`.

On the sample pilot MuData the run takes about 10 s. Each of the 7 self-TSS positive controls that are in the expression matrix comes out a hit, and every hit at BH < 0.1 is also a hit in the original Gasperini results.

```bash
igvfagent gasperini-pipeline setup
igvfagent gasperini-pipeline inspect
igvfagent gasperini-pipeline run --gasperini-test --compare Data/GasperiniPipeline/GSE120861_all_deg_results.pilot.txt.gz
igvfagent gasperini-pipeline selftest
```

Agent tools: gasperini_pipeline_setup, gasperini_pipeline_run, gasperini_pipeline_qc_filter, gasperini_pipeline_mudata, gasperini_pipeline_guide_table, gasperini_pipeline_composition, gasperini_pipeline_multiseq, gasperini_pipeline_compare, gasperini_pipeline_inspect.

### ABC pipeline (`abc-pipeline`)

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

### scE2G pipeline (`sce2g-pipeline`)

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

### IGVF CRISPR Perturb-seq pipeline (`crispr-pipeline`)

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

### 2nd CRISPR Jamboree analyses (`crispr-jamboree2`)

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

### 3rd CRISPR Jamboree single-cell pipeline (`crispr-jamboree3`)

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

### 4th CRISPR Jamboree inference benchmark (`crispr-jamboree4`)

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

### Principal pseudobulks (`principal-pseudobulks`)

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

### Single-cell CRISPR differential expression (`sc-crispr-de`)

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

### Spatial-ATAC-Hi-C (spatial 3D genome + chromatin accessibility)

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
[`Docs/Skills/SPATIAL_ATAC_HIC_SKILLS.md`](Docs/Skills/SPATIAL_ATAC_HIC_SKILLS.md).

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

### ChIP-Atlas reprocessed peak archive

Browse and pull from the [ChIP-Atlas](https://chip-atlas.org) (Ohta/Oki/DBCLS)
reprocessed archive of public ChIP-seq / ATAC-seq / DNase-seq / Bisulfite-seq
peaks — a clean-room, stdlib-only client over the public HTTP surface. Ten
genomes, four −log10(q) thresholds, per-experiment BigWig/BigBed/BED, assembled
all-peaks BEDs, pre-computed Target-Genes tables, and queued TF-enrichment jobs.
Playbook: [`Docs/Skills/CHIPATLAS_SKILL.md`](Docs/Skills/CHIPATLAS_SKILL.md).

```bash
igvfagent chipatlas list-antigens   --genome hg38 --cell-type Blood
igvfagent chipatlas search          --genome hg38 --antigen GATA1 --cell-type Blood
igvfagent chipatlas target-genes    --antigen H3K4me3 --distance 5000
igvfagent chipatlas submit-enrichment --genes my_genes.txt --genome hg38
```

### MaveDB mapping (incl. SGE cDNA path)

Map [MaveDB](https://www.mavedb.org) multiplexed-assay scoresets to genomic
coordinates via the Ensembl REST API (no UTA / SeqRepo / BLAT dependency). In
addition to the protein-coordinate VAMP-seq path, a dedicated **SGE (Saturation
Genome Editing) path** parses the full HGVS-c grammar used by SGE scoresets
(CDS, intronic, 5′UTR, 3′UTR) and emits VCF-4.2 with score-bearing INFO fields.
This is the path the **Waters 2024 BAP1** and **Buckley 2024 VHL** benchmarks
exercise. Playbook: [`Docs/Skills/MAVEDB_MAPPING_SKILL.md`](Docs/Skills/MAVEDB_MAPPING_SKILL.md).

```bash
igvfagent mavedb map-scoreset --urn urn:mavedb:00000097-0-1   # PTEN VAMP-seq
igvfagent mavedb map-scoreset --urn <BAP1-SGE-urn> --sge      # SGE cDNA path
```

### Synapse / Sage Bionetworks retrieval

Discover and download from [Synapse](https://www.synapse.org) — a clean-room
client over the public REST API (no `synapseclient` dependency). Anonymous read
of entity metadata, child listing, recursive walks, and full-text search; PAT-
authenticated download (`SYNAPSE_AUTH_TOKEN`) for controlled-access deposits
(PsychENCODE, AMP-AD/PD, ROSMAP) that IGVF distributes off-Portal. This skill
powers the **Deng 2024 cortex lentiMPRA** benchmark. Playbook:
[`Docs/Skills/SYNAPSE_RETRIEVAL_SKILLS.md`](Docs/Skills/SYNAPSE_RETRIEVAL_SKILLS.md).

```bash
igvfagent synapse entity   --syn syn21392931                 # metadata (anon)
igvfagent synapse walk     --syn syn21392931 --max-depth 3   # recursive walk
igvfagent synapse search   --query "lentiMPRA cortex" --limit 20
igvfagent synapse download --syn synXXXXXXXX --out-dir Data/Input   # needs PAT
```

### Assay calibration → ACMG/AMP evidence (exCALIBR)

Turn a raw multiplexed-assay score into a **clinically usable evidence
strength**. This is the last mile of the MAVE chain: `mavedb map-scoreset`
gives a variant genomic coordinates, and `calibrate` says how much a given
assay score is actually worth as PS3 / BS3 evidence — supporting, moderate,
strong or very strong — instead of leaving the reader with a bare number.

Clean-room reimplementation of [exCALIBR](https://github.com/rosstewart/exCALIBR)
(MIT), which implements the gene-based calibration method of Zeiberg et al.
(*bioRxiv* 2025.04.29.651326) on top of Tavtigian's Bayesian reading of the
ACMG/AMP guidelines. The chain: label variants into P/LP, B/LB, gnomAD-population
and synonymous samples → fit a **constrained skew-normal mixture** by EM (shared
components, per-sample weights, monotone density-ratio constraint enforced by
binary search inside every M-step) → **bootstrap** it → EM-estimate the prior
P(pathogenic | population) → build LR⁺(score) → solve for **Tavtigian's C** →
emit the score window that earns each evidence strength. Playbook:
[`Docs/Skills/ASSAY_CALIBRATION_SKILL.md`](Docs/Skills/ASSAY_CALIBRATION_SKILL.md).

```bash
# 0) Evidence thresholds alone — instant, no data needed
igvfagent calibrate thresholds --prior 0.1
#    -> C = 348; LR+ 2.08 (supporting) / 4.32 (moderate) / 18.7 (strong) / 348 (very strong)

# 1) Label a scoreset into the four calibration samples
igvfagent calibrate prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
#    or from IGVFagent's own MAVE chain, joined to a ClinVar release:
igvfagent calibrate prepare --mapped Docs/MaveDB/<run>/mapped.tsv \
    --clinvar-tsv variant_summary.txt.gz --gnomad-tsv gnomad_sites.tsv

# 2) Calibrate (long job: progress heartbeat + resumable ledger,
#    2c-vs-3c model selection, calibration JSON + figure)
igvfagent calibrate run --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021 \
    --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
igvfagent calibrate run --table scores.csv --name MSH2 --resume   # continue

# 3) Interpret new variants with the calibration
igvfagent calibrate assign --calibration MSH2_Jia_2021_2c_calibration.json \
    --scores my_variants.csv
#    -> score, evidence_points, acmg_evidence ("PS3 moderate", "BS3 strong", …)

# 4) Validate the port on this machine (~1 min)
igvfagent calibrate selftest
```

Every numeric component of the default path was checked against the upstream
implementation:
Tavtigian's C matches over a 12-prior grid (including the `original` and
`strict` rule variants), the constrained-EM iterates are **bit-identical** for
60 consecutive steps at 2 and 3 components, and the prior EM, LR⁺ → point-range
conversion, and paired model-selection test all agree to machine precision. The
Pillar/IGVF-format loader reproduces upstream's variant labelling exactly
(identical variant-ID sets per sample across ClinVar releases and star
thresholds). On the MSH2 (Jia 2021) example both implementations select the same
model and the same set of evidence strengths, with thresholds inside ~1–2 % of
the score range at equal bootstrap budgets. Unlike upstream the run is fully
seeded, so a rerun reproduces the calibration exactly.

## Reproducibility benchmark suite

IGVFagent ships a **22-paper reproducibility benchmark suite** in
[`Benchmarks/`](Benchmarks/README.md): recent Nature / Cell / Science /
Nat Genet / Nat Methods / Genome Biol papers whose published analyses IGVFagent
reproduces **directly from public data**. Each benchmark is a self-contained
directory — data sources, a deterministic `run.sh`, machine-readable
`expected.json` checks, regenerable figures, and a paper-vs-IGVFagent
`README.md` with a *Concordance / Verdict / Honest caveats* structure. A
stdlib-only scorer (`concordance.py`) turns each run into pass/fail checks.

| # | Paper | Skill exercised | Headline result |
|---|---|---|---|
| ⭐ | **Matreyek 2018** PTEN VAMP-seq *Nat Genet* | `mavedb` | 4/4 concordance checks pass (8,000 variants) |
| 1 | **Waters 2024** BAP1 SGE *Nat Genet* | `mavedb` (SGE path) | LOF +1.1 %, GOF +7.7 % vs paper |
| 2 | **Buckley 2024** VHL SGE *Nat Genet* | `mavedb` (SGE path) | 2,268 / 2,268 variants recovered |
| 3 | **Zou 2024** ChIP-Atlas 3.0 *Nucleic Acids Res* | `chipatlas` | 815 TFs catalogued; GATA1 rank #7 |
| 4 | **Agarwal 2025** lentiMPRA *Nature* | `mpra` | 3/3 discovery artefacts written |
| 5 | **Yao 2024** ENCODE4 CRISPRi *Nat Methods* | `encode` (FCE path) | 368 CRISPR-screen FCEs enumerated |
| 6 | **Mitra 2024** SCARlink multiome *Nat Genet* | `multiome` | 505 multiome AnalysisSets; 4/4 content types |
| 7 | **Weinstock 2024** CD4 CRISPR *Cell Genomics* | `perturb-catalog` + `geo` | 1,197 CRISPR-screen datasets; 99.7 % CRISPRn |
| 8 | **Zheng 2024** in-vivo Perturb-seq *Cell* | `geo` + `sc-analyze` | GSE249416 metadata + 9-file inventory |
| 9 | **Martyn 2025** Variant-FlowFISH *Cell* | `flowfish` | end-to-end chain: 20 elements → 7 Significant |
| 10 | **Joung 2025** TF Perturb-seq *Nat Genet* | `perturb-catalog` | modality scale confirmed (15 datasets) |
| 11 | **Deng 2024** cortex lentiMPRA *Science* | `mpra` + `synapse` | 166-node Synapse walk; 12/12 annotations recovered |
| 12 | **Gschwind 2023 / Sheth 2024** ENCODE-rE2G & scE2G *bioRxiv* | `portal` + `synapse` | ENCODE-rE2G base AUPRC 0.634 reproduces the published 0.634 exactly |
| 13 | **Liu 2025** kidney multiome scorecard *Science* | `open4gene` + `figshare` | port matches R `pscl::hurdle` β to r=1.0; 125,699 unique peaks (exact) |
| 14 | **Travaglini 2020** lung atlas *Nature* | `sc-analyze` (**full local repro**) | 49 clusters vs 46 author types, AMI 0.81; 9/10 markers — 8/8 |
| 15 | **Trevino 2021** cortex multiome *Cell* | `multiome peak2gene` (**full local repro**) | all 26 lineage genes linked; 93 % positive cis links — 4/4 |
| 16 | **deMULTIplex2 / Stoeckius 2018** cell hashing *Genome Biol* | `multiseq` (**full local repro**) | all 8 donor HTOs; 83 % singlets — 7/7 |
| 17 | **Rosenberg 2018** SPLiT-seq CNS *Science* | `splitseq` (**full local repro**) | 156,049-nucleus atlas (exact); all 8 CNS lineages — 4/4 |
| 18 | **Ma 2020** SHARE-seq mouse skin *Cell* | `share` (**full local repro**) | 34,774-cell skin set (exact); all 23 author cell types; Leiden AMI 0.63 — 5/5 |
| 19 | **Wang 2025** developing neocortex multiome *Nature* | `sc-analyze` (**full local repro**) | 232,328-nucleus atlas (exact); 29 author cell types; AMI 0.65 — 7/7 |
| 20 | **Zou/Shi 2026** scEPS GWAS × single-cell *medRxiv* | `sceps` port (**full local repro**) | microglia significantly AD-associated (d=4.76e-5, **Z=4.99, P=6e-7**); d>0 fraction 58 % vs paper 58 % |

| 21 | **Rosen 2025** MPRAsnakeflow *Genome Res / bioRxiv* | `mpraflow` + `mpralib` (**full local repro**) | **210,660 / 210,660 values byte-identical** to the authors' published IGVF artefact; 3/3 stated complexity figures reproduced to the digit — 8/8 |

Benchmarks 14–21 are **full local reproductions** — they download the public
data and run IGVFagent's real single-cell / multiome / GWAS-integration chain,
then score concordance against the authors' own results. The newest, **scEPS**
(Zou/Shi 2026), reimplements the full estimate→cluster→aggregate pipeline
(validated vs upstream, corr 1.0) and reproduces the paper's microglia AD-
association from public SEA-AD + Bellenguez-AD MAGMA scores. Most other
benchmarks run end-to-end as pure online calls; a few (e.g. **Zheng 2024**,
**Deng 2024**) need a user-fetched local file to complete their full chains.

```bash
# Verify the suite works (~60 s)
bash Benchmarks/matreyek2018_pten_vampseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark matreyek2018_pten_vampseq

# Run all online benchmarks and score them
bash Benchmarks/run_all.sh --online-only
.venv/bin/python Benchmarks/concordance.py --all
```
SCIP (PicoMILP backend) is BSD-style, free for any use. Selected
subnetworks flow back into the central DuckDB warehouse as `edges`
rows tagged with `upstream='network:carnival:<label>'` or
`upstream='network:steiner:<label>'` so downstream embedding / 
foundation-model training treats them like any other evidence stream.
Three agent tools registered: `network_demo`, `network_carnival`,
`network_steiner`.

## Extending IGVFagent

IGVFagent is **expandable by design**: the 43 built-in skills and 146+
registered tools are a starting point, not a ceiling. A user-extension
framework absorbs **your own skills and tools** at startup — no core-code
edits, no re-install, no registration step. Anything you drop into an
extension directory becomes a first-class citizen: it shows up in
`igvfagent --help` and `igvfagent tools`, the `igvfagent ask` ReAct agent
can plan with and call it, the Streamlit UI lists it in the tool picker,
and its outputs are harvested into the local knowledge graph exactly like
built-in results.

There are four extension surfaces, from zero-code to full-code:

| Surface | You write | Where it plugs in |
|---|---|---|
| **Custom tool** | one YAML/JSON manifest (no code) | LLM-callable tool in `ask` + UI, wrapping *any* executable or an existing `igvfagent` subcommand |
| **Custom skill** | one Python file with `main()` | first-class `igvfagent <name>` subcommand with the post-run KG harvest |
| **Prompt skill** | one `SKILL.md` file | Claude Code slash-skill orchestrating existing CLI steps |
| **Playbook** | one YAML file | deterministic multi-step tool chain, run via `igvfagent playbook` |

Extensions are discovered from these locations (scanned in order; first
definition of a name wins, and built-in names can never be shadowed):

1. every directory in `$IGVF_USER_EXT_DIR` (`:`-separated) — for testing
   or shared lab locations
2. `~/.igvfagent/` — per-user, works from any checkout
3. `<repo>/UserExtensions/` — per-checkout, committable so a whole lab
   shares one extension set through git

Each location uses the same two subfolders: `tools/` (manifests) and
`skills/` (Python modules). Copy-paste-ready templates live in
[`Docs/Examples/user_extensions/`](Docs/Examples/user_extensions/), and

```bash
igvfagent extensions          # what was discovered, from where, + any
igvfagent extensions --json   # skipped/malformed definitions explained
```

shows exactly what the framework picked up. Discovery is defensive: a
malformed manifest or a skill that fails to import is skipped with a
diagnostic in `igvfagent extensions` and can never break the core CLI or
agent runtime.

### Add a custom tool (YAML manifest — no code)

A *tool* is what the LLM agent calls during `igvfagent ask`. To teach the
agent a new capability, describe your command in a manifest — here a lab
script that computes GC content:

```bash
mkdir -p ~/.igvfagent/tools
cat > ~/.igvfagent/tools/gc_content.yaml <<'EOF'
name: gc_content
description: >
  Compute per-sequence GC content for a FASTA file and write a TSV
  summary. Use when the user asks for GC% or base-composition QC.

# argv of YOUR program — any language, any location
command: ["python3", "/home/me/bin/gc_content.py"]

parameters:               # JSON Schema for what the LLM may pass
  type: object
  properties:
    fasta:  {type: string, description: Path to the input FASTA.}
    window: {type: integer, default: 0}
  required: [fasta]

positional: [fasta]        # emitted as a positional argument
flag_map:
  window: "--window"       # emitted as `--window <value>`
EOF

igvfagent extensions               # confirm it was absorbed
igvfagent tools | grep -A3 gc_content
igvfagent ask "run GC content QC on Data/my_seqs.fa"
```

The only contract your program must honour: read arguments from argv,
print results to stdout, exit 0 on success, and announce output files as
`Report: <path>` / `Manifest: <path>` lines — the same convention every
built-in follows, which is how the agent chains your artefacts into its
next step and how the localstore harvester grows the local KG from them.

Instead of `command:`, use `cli:` to re-surface an existing `igvfagent`
subcommand under your own name, defaults, and description (see the
`kg_gene_quick.yaml` template) — useful for lab-specific shortcuts.
Optional mapping fields mirror the built-in registry: `positional`,
`flag_map` (parameter → flag), `flag_repeat` (list values repeat the
flag), `bool_flags` (bare flags). One manifest may also carry several
definitions under a top-level `tools:` list.

### Add a custom skill (Python subcommand)

A *skill* is a human-facing subcommand. Drop a Python file with a
`main()` into `skills/` and it becomes `igvfagent <name>` (filename stem,
`_` → `-`):

```bash
mkdir -p ~/.igvfagent/skills
cp Docs/Examples/user_extensions/skills/variant_bed_export.py ~/.igvfagent/skills/

igvfagent --help                    # now lists it under "User skills"
igvfagent variant-bed-export my_variants.tsv --out variants.bed
```

Conventions worth copying from the template:

- the module docstring's first line is the description shown in
  `igvfagent --help` and `igvfagent extensions`;
- `main()` owns its own argparse parser and returns an int exit code —
  the dispatcher hands over `sys.argv` exactly as for built-ins;
- print every artefact as `Report: <path>` so downstream steps and the
  local-KG harvest (which runs after user skills too) pick it up.

To expose the same logic to the LLM agent as well, add a small `cli:`- or
`command:`-style manifest next to it — the skill/tool split is the same
one the built-ins use (human CLI surface vs. curated LLM tool surface).

### Add a prompt skill or playbook

The remaining two surfaces need no Python at all:

- **Prompt skills** — drop a `.claude/skills/<name>/SKILL.md` (YAML
  frontmatter + a markdown workflow whose steps are `igvfagent …`
  commands) and Claude Code picks it up with zero registration. The
  eight shipped `igvf-*` skills under `.claude/skills/` are the
  reference implementations; see
  `Docs/Skills/PROMPT_SKILLS_INDEX.md`.
- **Playbooks** — a YAML file listing `steps: [{tool, args}]` over any
  registered tools (including your user tools), with `${param}`
  interpolation and `--param key=value` overrides. Put it in
  `Docs/Playbooks/` and run `igvfagent playbook <name>`; the schema is
  documented at the top of `Scripts/_playbook.py`.

Together the four surfaces mean IGVFagent can absorb a lab's entire
private toolbox — QC scripts, internal APIs, bespoke pipelines — while
keeping the audited, reproducible execution model of the built-ins.

## Deployment with LLM agents

Every skill is a CLI tool, so any orchestration layer that can run shell
commands and read files can drive the agent.

### Codex API

Use Codex as the coding/data agent layer and give it this repository as the
workspace root. Recommended runtime instruction:

```text
Use /path/to/IGVFagent as the project root. Put scripts in Scripts, data in
Data, and logs/reports in Docs/Logs or Docs. Use the local CLI skills
before writing new one-off code.
```

### Claude API

Use Claude with a tool runner that exposes shell commands inside the repo.
Restrict file access to the IGVFagent folder and call the scripts as tools.
Sample tool targets:

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/annotate_variant_list.py --max-rows 10
python3 Scripts/advanced_variant_analysis.py run --input <csv> --label <run-id>
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
```

### Ollama local models

Install Ollama and pull a coding model:

```bash
ollama pull qwen3
ollama pull llama3.1
ollama serve     # starts http://localhost:11434
```

Then connect your preferred local agent runner to the Ollama endpoint. For
workstation use, Qwen-class coding models are useful for command planning and
report drafting; the Python skills do the deterministic data access and
annotation.

## Variant lists and your own data

The variant-annotation skills accept any user-provided CSV via `--input`.
Recognized identifier columns (case-insensitive):

- `rsid` / `rsID` / `dbSNP` — e.g. `rs58658771`
- `chrom`, `pos`, `ref`, `alt` — GRCh38 coordinates and alleles
- `spdi` — NCBI SPDI string
- `hgvs` — e.g. `NC_000019.10:g.44908822C>T`

Additional user columns (locus name, phenotype, prior, notes, experimental
effect, p-value) are preserved through annotation and used by
`advanced_variant_analysis.py` when you pass `--experimental`/`--outcome`.

A minimal `Data/Input/VariantList/example_variants.csv` is included for
smoke testing. Replace it with your own list — the `.gitignore` excludes
everything in this folder except the README and the example CSV.

**Do not commit confidential or pre-publication variant lists.**

## Workstation notes

- Keep large downloads on a disk with enough space. Full cCRE, rE2G, and
  single-cell linkage corpora can be many gigabytes.
- Prefer `*-manifest` commands before `*-download` commands.
- Use `--max-rows` / `--limit` for smoke tests before launching full runs.
- Commit scripts and docs, but do not commit large downloaded data unless
  the project policy explicitly allows it.
- For reproducibility, preserve the generated manifest CSVs and Markdown
  reports under `Data/Manifests/` and `Docs/`.
- Set `IGVF_PROJECT_ROOT` if you want to run the scripts from a directory
  other than the repo, or if the repo lives on an external volume.

## Security

- The repository ships with an aggressive `.gitignore` that excludes runtime
  outputs, caches, logs, and any `.env` / cookie / token files.
- All credentials are read from environment variables; nothing is hardcoded.
- Logs record request URLs and HTTP status codes only — never credential
  headers.
- Never commit `.env`, browser cookie exports, OAuth tokens, API keys, or
  unreleased / pre-publication data.

## References

IGVFagent's analytical skills are clean-room reimplementations that learn
from a number of public reference pipelines and methods papers. Algorithms
are paraphrased from published descriptions and supplementary methods —
no source code is copied verbatim — but every skill below is built on the
shoulders of the work in these repositories. We thank the authors and
maintainers for releasing their code openly.

### Reference GitHub repositories

| IGVFagent skill | Reference repository | License | What we absorb |
|---|---|---|---|
| **SHARE-seq** joint scATAC + scRNA QC (`share`) | [broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline) | MIT (Broad Institute, 2021) | Round-1/2/3 24-mer barcode demultiplex (1-Hamming + ±1 bp shift), `bam_to_fragments` Tn5 +4/−4 shift, TSS enrichment (Ma 2020 formula with 0.2 floor), per-barcode FRIP, joint cell calling thresholds, Jaccard multiplet detection. |
| **STARR-seq** allelic test (`starrseq`) | [gaochengwen/STARR-seq-Data-Analysis](https://github.com/gaochengwen/STARR-seq-Data-Analysis) | No LICENSE (treated as default copyright) — every line is a clean-room rewrite. | TPM-style counts QC, Spearman D-stat outlier flagging, per-(SNP, Allele) aggregation, log activity = log(RNA/DNA), and the `mpra::mpralm`-style allelic test (paraphrased from the underlying limma + voom + eBayes methods). |
| **CRISPRi Flow-FISH** screen (`flowfish`) | [EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline) | MIT (Engreitz Lab, 2021) | Per-guide log-normal MLE on bin-multinomial counts with EM treatment of an "outside" overflow bin, real-space conversion + negative-control rescaling, Mann-Whitney U + Welch t-test per element, BH-FDR, `Significant` / `Regulated` output convention. |
| **Base-editing screens** — base-edit-aware guide assignment + activity (`bean`) | [pinellolab/crispr-bean](https://github.com/pinellolab/crispr-bean) · docs [pinellolab.github.io/crispr-bean](https://pinellolab.github.io/crispr-bean/) | **AGPL-3.0 — no source copied.** IGVFagent is Apache-2.0, so the METHOD was read from `bean/mapping/GuideEditCounter.py` and reimplemented; vendoring would place AGPL over this code. To use BEAN's own model, invoke `bean` as a separate program. | BEAN's `mask_sequence` matching: the edited base is normalised to its product on BOTH sides before comparison (`edited_base=['A>G']` maps every A to G), rather than allowing free mismatches — so sequencing error and similar guides are not absorbed. A base editor self-edits the guide's own locus, so exact matching discards those reads: measured on IGVFDS6464SOVZ (8,192-guide ABE library) **36.7% of reads assigned vs 62.5%**, with C→T recovering only +0.1%, confirming an adenine editor. The editor is detected from the library's own guide names AND by measuring which masking recovers reads. Per-guide self-edit rate is used as the editing-activity estimate — what BEAN's activity normalisation rests on (BEAN = *Base Editing screens' Activity-Normalized variant effect size estimation*). **Handed off, not reimplemented:** `bean analyze --run-bean` writes BEAN's own tables, runs `bean create-screen` and `bean run sorting variant`, and reports its per-target posteriors (`mu`, `mu_sd`, `mu_z`, `n_guides`) — so variant-level aggregation and posterior uncertainty come from BEAN itself, not from a reimplementation. **NOT OBTAINABLE on IGVF-published data** (measured, not assumed): BEAN's accessibility-covariate model (`--scale-by-acc`) *and its own activity normalisation* both route through MixtureNormal, whose data class reads `screen.layers['X_bcmatch']` unconditionally (`bean/preprocessing/data_class.py:309`) — counts assigned by guide barcode, which IGVF does not publish, so both raise `KeyError: 'X_bcmatch'`. The two models that avoid that read (`--uniform-edit`, `--const-pi`) reject `--guide-activity-col` outright, so IGVFagent's own `log2_per_edit` is the only activity-aware estimate available here rather than an approximation of a better one; BEAN contributes the posterior under its Normal model, labelled as not activity-normalised. Reporter-allele and bystander/tiling analysis need CRISPResso2 alignment and a `reporter` column IGVF leaves empty for all 8,192 guides. The bcmatch/semimatch split needs that same barcode, so masked-sequence collisions stay ambiguous instead of being assigned. `bean run sorting` needs an unsorted bin as `--control-condition`, and finding one is a **library** question, not a naming one: this screen's tails are `18loci_uptake_Rep1_bottom20` while its unsorted bins are `B1.1_18loci_Rep1_bulk` (`IGVFDS7513KJPR` and three siblings). Grouping a screen's libraries by the series name parsed from the alias split it in two and lost them, and this table previously stated — wrongly — that `IGVFDS6464SOVZ` had published no unsorted bin at all. Membership is now keyed on `construct_library_sets`, which finds all 20 libraries and still separates the ABE screen (`IGVFDS5242VHWM`) from the CRISPRi screen (`IGVFDS2978DRWC`) sharing both the `18loci` stem and the assay title. `bean run sorting variant --control-condition bulk` then completes on `IGVFDS6464SOVZ`, returning **1,659 target posteriors**. |
| **FASTQ → count matrix** (`raw-pipeline`, `sc-analyze`) | [pachterlab/kb_python](https://github.com/pachterlab/kb_python) (kallisto \| bustools) + [scverse/scanpy](https://github.com/scverse/scanpy) | BSD-2 (kallisto/bustools) · BSD-3 (scanpy) — **invoked, not reimplemented.** | `kb ref` / `kb count` for pseudoalignment and barcode/UMI counting, including the `kite` feature-barcode workflow for guide capture; scanpy for QC, HVG selection, PCA, Leiden and marker tests. This is the only aligner shipped, which is why genomic assays (ATAC-seq, Hi-C) are refused rather than mis-quantified — 40.5% of Portal MeasurementSets are not transcript assays. |
| **MPRA** allelic activity + skew (`mpra`) | [tewhey-lab/MPRASuite](https://github.com/tewhey-lab/MPRASuite) | Apache-2.0 (Tewhey Lab) | DESeq2 NB GLM Wald test (via `pydeseq2`), summit-shift size-factor renormalization (MPRAmodel), allelic-skew paired t-test of per-replicate log2(RNA/DNA) with BH-FDR. |
| **MPRA** QC + counts handling (`mpra`) | [WangLabTHU/esMPRA](https://github.com/WangLabTHU/esMPRA) | Ambiguous (no clean SPDX) — treated as clean-room. | Count-based replicate concordance Pearson r matrix, barcodes-per-oligo and counts-per-oligo histograms. |
| **10x Multiome** analytics (`multiome qc-atac / joint-qc / lsi / wnn / peak2gene / showcase`) | [10XGenomics/analysis_guides](https://github.com/10XGenomics/analysis_guides) (10x Genomics, no LICENSE) + [stuart-lab/signac](https://github.com/stuart-lab/signac) (Stuart Lab, MIT) | Repo has no LICENSE → clean-room; Signac MIT is OK to cite. | Per-barcode TSS enrichment (Signac `TSSEnrichment` ±100 bp center / ±900-1000 bp flank formula), nucleosome signal (mono / NFR length ratio), FRIP from fragments TSV, Signac-convention joint QC thresholds (RNA UMI 1k-25k, ATAC frag 1.8k-100k, TSS>1, nuc<2, FRIP>0.15), TF-IDF + truncated SVD with depth-dim drop (Signac `RunSVD` + `DepthCor`), WNN joint embedding (Seurat 5 / Hao 2021, via `muon`), peak-to-gene correlation. |
| **10x Multiome** Python extensions (`multiome da-peaks / atac-spectral / chromvar / css / multivi`) | [quadbio/scMultiome_analysis_python_vignette](https://github.com/quadbio/scMultiome_analysis_python_vignette) (Treutlein Lab ETH, no LICENSE) | clean-room — algorithms re-implemented from published descriptions only; no source copied | Differential accessibility on TF-IDF peaks (Wilcoxon), snapATAC2-style Jaccard-Laplacian spectral embedding (alternative to LSI), chromVAR-style TF motif deviations with K=50 GC-matched background motif sets, Cluster Similarity Spectrum (He 2020 Genome Biol) batch correction as a BSD-3-friendly alternative to GPL Harmony, optional MultiVI (Ashuach 2023) deep joint VAE via scvi-tools (BSD-3, optional install). SCENIC+ / pycisTopic (academic non-commercial) are documented as external workflows only. |
| **VAMP-seq** abundance scoresets (`proteomics vampseq-analyze`) | [FowlerLab/VAMPseq](https://github.com/FowlerLab/VAMPseq) | MIT (Fowler Lab) | Score-density distribution, residue×AA heatmap, per-residue mean (± IQR) AND per-residue median + 3-residue moving average, N×N replicate concordance matrix, abundance-class bar, nonsense-by-position QC scatter, PyMOL `.pml` overlay. |
| **MaveDB → genomic coords** (`mavedb map-scoreset / showcase`) | [ave-dcd/dcd_mapping](https://github.com/ave-dcd/dcd_mapping) (AVE consortium, MIT) | clean-room reimpl using Ensembl REST API instead of UTA + SeqRepo + BLAT | HGVSp parsing, Ensembl-canonical-transcript lookup, protein-position → genomic-codon mapping via `/map/translation/{ENSP}`, codon validation against the genomic reference, single-nt-change enumeration with ambiguity flags, VCF-4.2 emission with INFO fields (AA_REF/AA_ALT/AA_POS/TRANSCRIPT/SCORE/AMBIG). |
| **MULTI-seq / Cell Hashing** demultiplex (`multiseq`) | [Gartner-Lab/deMULTIplex2](https://github.com/Gartner-Lab/deMULTIplex2) | MIT (Gartner Lab) | NB-GLM tag classifier with randomized-quantile residuals, cosine-based normalization, per-tag histograms + call heatmaps. |
| **MULTI-seq** original method | — (paper: [McGinnis 2019 *Nat Methods*](https://www.nature.com/articles/s41592-019-0433-8)) | n/a | Lipid-modified oligo (LMO) sample barcoding inspiration; used as the upstream context for our demultiplexer. |
| **Network integration** (`network carnival`, `network steiner`) | [saezlab/CORNETO](https://github.com/saezlab/corneto) | GPL (CORNETO) — **runtime dep avoided.** Clean-room MILP reimplementation in pure `cvxpy`. | CARNIVAL signed-perturbation → signed-measurement MILP, Prize-Collecting Steiner-Tree formulation. |
| **Tabula Sapiens 2.0** human cell atlas (`tabula`) | [czbiohub-sf/tabula-sapiens](https://github.com/czbiohub-sf/tabula-sapiens) (paper: Tabula Sapiens Consortium, *Cell* 2026; figshare 27921984, GEO GSE306755) | BSD-3-Clause — clean-room; no source copied. | Streaming per-tissue h5ad aggregation via direct HDF5 CSR row slicing; tspex τ cell-type specificity; GO enrichment of non-specific TFs; SCENIC-style regulon inference (GRNBoost2 co-expression → motif support → AUCell); CDKN2A+MKI67− senescence with donor-balanced strata, per-stratum Wilcoxon DE and cross-donor replication; DecontX variational-EM ambient-RNA model; cNMF consensus modules; pseudobulk negative-binomial sex DE; donor clinical-metadata query (ChatTS core). |
| **Human Transcription Factors** database (`humantfs`) | [humantfs.ccbr.utoronto.ca](https://humantfs.ccbr.utoronto.ca) (Lambert et al., *Cell* 2018) | Upstream terms; data downloaded at build time, never vendored. | Local SQLite mirror of all 1,639 curated TFs with DBD family, binding mode, motif status, CIS-BP motif records and Ensembl/HGNC/Entrez/InterPro/PDB cross-references; is-tf / lookup / list / motifs / families / export. |
| **TF regulon inference** (inside `tabula tf-regulons`) | [aertslab/pySCENIC](https://github.com/aertslab/pySCENIC) | **GPL-3 — runtime dep avoided.** | Clean-room GRNBoost2-style gradient-boosted co-expression, motif-supported regulon pruning, and rank-based AUCell activity. Motif support is weaker than cisTarget (which needs multi-GB ranking databases); the difference is reported in every run's summary. |
| **Pseudobulk differential expression** (inside `tabula sex-de`) | edgeR (Robinson/McCarthy/Smyth) | **GPL — runtime dep avoided.** | Clean-room negative-binomial Wald test on donor pseudobulk with method-of-moments dispersion. Concordant in direction and ranking, not in exact p-values — edgeR shrinks dispersions empirically. |
| **Ambient RNA / modules / specificity** (inside `tabula`) | DecontX (celda), cNMF, tspex | MIT | Clean-room variational EM, consensus NMF with density filtering, and the τ statistic. DecontX and scVI outputs stored in the published h5ads let these be validated against upstream, not merely run. |
| **CRISPR-SURF** tiling-screen deconvolution (`crispr-surf`) | [pinellolab/CRISPR-SURF](https://github.com/pinellolab/CRISPR-SURF) (Hsu et al., *Nat Commun* 2018) | AGPL-3.0 — clean-room reimplementation of the METHOD; no source read or ported. Installing it was attempted first and is impossible: not on PyPI, and its source does not parse under Python 3. | L1-regularised deconvolution of `y = A·beta` against a triangular perturbation kernel sized by nuclease (20 bp Cas9 / 250 bp CRISPRi), empirical null resampled from the screen's own negative-control guides, BH across bins, adjacent significant bins merged into regions. |
| **ABC** Activity-by-Contact enhancer→gene (`abc`) | [broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) (Fulco *et al.* 2019; Nasser *et al.* 2021) | MIT — clean-room; no source copied. | Activity as the geometric mean of accessibility and H3K27ac signal, contact from Hi-C or the genome-wide power law (γ = −0.87), normalised to a per-gene share over a 5 Mb window. Does **not** call peaks or read BAMs. |
| **scE2G** single-cell enhancer→gene features (`sce2g-predict`) | [EngreitzLab/scE2G](https://github.com/EngreitzLab/scE2G) | MIT — clean-room feature computation; the published model's **trained weights are not reproduced**, and the output is labelled features, never "the scE2G score". | Kendall rank correlation of peak accessibility against gene expression across KMeans metacells pooled on RNA profile, combined with an ABC share and distance. |
| **chromap** scATAC alignment | [haowenz/chromap](https://github.com/haowenz/chromap) v0.3.2 | MIT — built from the release tarball and invoked as a separate binary; no linking, no vendoring. | Not reimplemented. Completes [IGVF/atomic-workflows](https://github.com/IGVF/atomic-workflows)' two-module tool set alongside kb-python. |
| **CRISPResso2** genome-editing outcomes (`crispresso`) | [pinellolab/CRISPResso2](https://github.com/pinellolab/CRISPResso2) (Clement *et al.*, *Nat Biotechnol* 2019) v2.3.4 | Installed as a separate program invoked by subprocess — no linking, no vendoring — in `/opt/bean-venv` beside BEAN. | Not reimplemented. Wrapper locates the binary off the app venv's PATH, globs the quantification table whose name carries the run name, and takes a single `A,G` instead of two conversion flags. |
| **scNT-seq** metabolic labelling (`scnt-seq`) | [hongjie7/scNT-seq_pipeline](https://github.com/hongjie7/scNT-seq_pipeline) (Qiu *et al.*, *Nat Methods* 2020) | Clean-room; no source ported. | 4sU T>C conversion counting from BAM MD tags, strand-aware (A>G on reverse-strand reads), with the rate of every other substitution type reported as background and a warning when labelling is not clearly above it. |
| **BCalm** barcode-level MPRA (`bcalm`) | [kircherlab/BCalm](https://github.com/kircherlab/BCalm) | Clean-room; no source ported — BCalm is R on limma. | Per barcode × replicate log2((RNA+1)/(DNA+1)), per-element variances moderated by empirical Bayes with d0/s0² fitted by method of moments on the scaled-F. voom precision weights not applied; concordant ranking, not identical p-values. |
| **SCTransform** normalisation (`sctransform`) | [satijalab/seurat](https://github.com/satijalab/seurat) (Hafemeister & Satija, *Genome Biology* 2019) | MIT — clean-room; no source ported, Seurat being R. | Per-gene NB regression on log10 depth, parameters regularised by a Gaussian kernel over log10(gene mean), Pearson residuals clipped to √(n/30). Seurat's `ksmooth` heuristic and `MASS::theta.ml` are not reproduced. |
| **CRISPR-Correct** imperfect guide mapping (`guide-map`) | [pinellolab/CRISPR-Correct](https://github.com/pinellolab/CRISPR-Correct) | Clean-room; the matching primitives were already in this repo and are now exposed. | Prefix-indexed exact match at every guide length and offset, a Hamming budget for sequencing error, and BEAN-style masking for base editing — kept separate, since free mismatches would absorb sequencing error and cross-map. |
| **sc-CRISPR-DE** single-cell CRISPR differential expression (`sc-crispr-de`) | [Gersbachlab-Bioinformatics/sc-crispr-de](https://github.com/Gersbachlab-Bioinformatics/sc-crispr-de) | MIT (Gersbach Lab Bioinformatics) — clean-room; no source copied. | Per-guide negative-binomial GLM (statsmodels, replacing R `MASS::glm.nb`) of guide-bearing cells against no-guide cells with library size as covariate; process-pool fan-out replacing the SLURM array; α-RRA gene-level aggregation with an empirical null from non-targeting guides, in place of the external FRACTEL package. brglm2 bias reduction is not implemented. |
| **Spatial-ATAC-Hi-C** spatial 3D genome + accessibility (`spatial-hic`) | [wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C) (paper: [Wang 2026 *Nat Methods*](https://doi.org/10.1038/s41592-026-03217-4), GEO GSE307620) | MIT (wangjuan001, 2026) — clean-room; no source copied. | 50×50 microfluidic pixel demultiplex (barcode-B+A concatenation, the upstream `bcsplit.py` offsets), per-pixel cis/trans/long-range contact QC, ArchR-style TSS enrichment, SnapATAC2-style gene activity score (promoter+body Tn5 insertions) and scGAD gene-associated domain score (Hi-C pair ends over gene bodies), scHiCluster convolution+RWR imputation, cooltools-style A/B compartment eigenvector, NeoLoopFinder-style per-pixel CNV + HMM segmentation, MAGIC spatial smoothing, per-pixel loop quantification with one-way-ANOVA cluster-specific loop calling and APA pileup, and 50×50 tissue-space rendering. |
| **Hi-C read processing** (upstream of `spatial-hic`) | [XiaoTaoWang/HiC_pipeline](https://github.com/XiaoTaoWang/HiC_pipeline) (runHiC) | **GPL-3.0 — runtime dep avoided.** | Nothing is imported. runHiC's `.pairs` output is the *input contract* for `spatial-hic`; alignment is the one upstream step with no tractable clean-room form, so it stays an external tool. |
| **Adapter/quality trimming** (upstream of `spatial-hic`) | [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) | **GPL-3.0 — runtime dep avoided.** | Nothing is imported. Documented as the external FASTQ-level step that precedes this skill's inputs. |
| **SPLiT-seq / Parse** pipeline (`splitseq`) | [Chipeyown/SPLiT-seq-Data-Analysis_Toolkit](https://github.com/Chipeyown/SPLiT-seq-Data-Analysis_Toolkit) | MIT (Chipeyown) | Vendored Rd1/Rd2/Rd3 96-well barcode whitelists, knee-plot cell calling, pre/post QC violins, per-Rd1-well summary heatmap, optional Scrublet doublet detection. |
| **Genotype demultiplexing** (referenced in `multiseq`) | [single-cell-genetics/vireo](https://github.com/single-cell-genetics/vireo) · [wheaton5/souporcell](https://github.com/wheaton5/souporcell) | Apache-2.0 / MIT | Documented alternatives when natural genotype variation is available instead of barcoded tags. |
| **Codex / agent runtime** | [openai/codex](https://github.com/openai/codex) | Apache-2.0 | Reference agent runtime; `cli.py` follows its tool-dispatch pattern. |
| **Local IGVF KG mirror** (`kg-mirror`) | [arangodb/arangodb](https://github.com/arangodb/arangodb) (Arango DB hosting the upstream KG) + [duckdb/duckdb](https://github.com/duckdb/duckdb) (local mirror) | Apache-2.0 (Arango community) + MIT (DuckDB) | Read-only AQL cursor streaming of every collection except `variants` + `variants_variants`; persist as ZSTD-Parquet shards under `Data/Warehouse/KG/`; register as DuckDB views in `Data/Warehouse/igvf_kg_mirror.duckdb` for offline querying. |
| **GO + Pathway enrichment validation** (`enrich ora / gsea / go / pathways / showcase`) | [zqfang/GSEApy](https://github.com/zqfang/GSEApy) | BSD-3-Clause (Z. Fang) | Enrichr-proxy ORA over GO_BP/MF/CC + Reactome 2022 + KEGG 2021 Human + WikiPathways 2024 Human + MSigDB Hallmark 2020 (hypergeometric + BH within library); Subramanian-style preranked GSEA (weighted-KS, permutation FDR) for ranked gene-score tables. Used as the validation layer over DEGs, CRISPR hits, and the gene-side of enhancer-gene linkages. |
| **IGVF Portal canonical-query layer** (`portal search / get / schema / list-types / endpoint-params / facets / report / batch-download`) | [IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp) (IGVF DACC, MIT) | clean-room reimpl, stdlib-only — no `igvf-client` PyPI pin | DACC-blessed `/search/` faceted query patterns with the field-filter DSL (dotted embeds, negation, range ops, list values); HTTP Basic auth via `IGVF_ACCESS_KEY`/`IGVF_SECRET_ACCESS_KEY`, or the same key pair in `Docs/Secret/IGVFportalAPI.txt`; `/report.tsv` TSV export; `/batch-download/` manifest of pre-signed S3 URLs (with optional `--fetch`); JSON-schema introspection via `/profiles/<Type>.json`; endpoint-param ↔ search-field map (the DACC UX trick) for any collection. |
| **IGVF Catalog (Knowledge Graph) canonical-query layer** (`catalog get-entity / search-region / find-associations / find-ld / resolve-id / list-sources`) | [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (IGVF DACC, MIT) | clean-room reimpl, stdlib-only — no `httpx` / `mcp` / `pydantic` runtime dep | Universal `get-entity` with 20+ ID auto-detection (rsID / SPDI / HGVS / CA / ENSG / HGNC / Entrez / ENSP / UniProt / MONDO / EFO / GO / HPO / DOID / UBERON / CL / CHEBI / OBA / DB / CHEMBL / CPX / R-HSA / GCST); `search-region` parallel fan-out over genes + variants + genomic-elements with K/M/G suffix region parser; `find-associations` by semantic category (genetic / regulatory / physical / functional / pharmacological / ld / coding / transcription) walking 18 edge endpoints; `find-ld` with r²/D'/ancestry buckets; `resolve-id` cross-reference projection; `list-sources` with per-endpoint source/method introspection; filter DSL with automatic `p_value=lte:5e-8` → `log10pvalue=gte:7.301` translation. |
| **Claude Code prompt-skill suite** (`.claude/skills/igvf-portal-facet-filter`, `igvf-catalog-variant-report`, `igvf-catalog-gene-dossier`, `igvf-catalog-dissect-locus`, `igvf-catalog-regulatory-landscape`, `igvf-catalog-disease-genes`, `igvf-catalog-ld-compare`) | [IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp) + [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (IGVF DACC, MIT) | clean-room paraphrase, retargeted at IGVFagent's `portal` + `catalog` CLI surface | Seven workflow prompt skills (1 portal + 6 catalog) auto-loaded by Claude Code when the user's prompt matches the description. Each is a structured multi-step procedure (resolve identifiers → fan-out per semantic relationship → cross-reference → compile a sectioned report) wired to IGVFagent's CLI commands rather than the upstream MCP tools. Cross-references between skills point downstream to IGVFagent-only follow-ups (`network steiner`, `enrich pathways`, `ccre`, `enhancer`). See `Docs/Skills/PROMPT_SKILLS_INDEX.md` for the suite-level index. |
| **ChIP-Atlas (Ohta/Oki) reprocessed peak archive** (`chipatlas list-genomes / list-qvalues / list-experiment-types / list-antigens / list-cell-types / search / get-experiment / download-experiment / assemble-bed / download-all-peaks / target-genes / submit-enrichment / poll-enrichment / showcase`) | [inutano/chip-atlas](https://github.com/inutano/chip-atlas) (Tazro Inutano Ohta / Shinya Oki / DBCLS, MIT) | clean-room reimpl of the public HTTP surface; stdlib-only (no `httpx`/`mcp`/`pydantic`) | Anonymous polite-1-rps client over three indirected hosts (`chip-atlas.org` JSON browse/search/POST-download, `chip-atlas.dbcls.jp/data` bulk static archive, `dtn1.ddbj.nig.ac.jp/wabi/chipatlas` WABI Enrichment/Diff queue). 10 supported genomes; 4 -log10(q) thresholds (05/10/20/50); browse antigens × cell-class with experiment counts; pull per-experiment BigWig/BigBed/BED; POST a `(genome × ag × cellClass × qval)` tuple to get an assembled all-peaks BED URL; HEAD-probe or stream the bulk `allPeaks_light.{genome}.{qval}.bed.gz` archive; discover and fetch pre-computed Target-Genes tables (e.g. `H3K4me3.5000.tsv` at ±5 kb TSS-proximity); submit + poll WABI Enrichment Analysis jobs for gene-list / BED-region TF over-representation. Cites Zou/Ohta/Oki *Nucleic Acids Res.* 2024 (doi:10.1093/nar/gkae358) and Oki *EMBO Rep.* 2018 (doi:10.15252/embr.201846255). Code MIT-compatible; data NBDC/DBCLS-licensed — we only fetch / link, never redistribute. |
| **Synapse / Sage Bionetworks retrieval** (`synapse entity / children / walk / search / download / write-playbook`) | [Sage-Bionetworks/synapsePythonClient](https://github.com/Sage-Bionetworks/synapsePythonClient) (Sage Bionetworks, Apache-2.0) | clean-room reimpl over the public REST API (`rest-docs.synapse.org`); pure `urllib` + `json`, no `synapseclient` runtime dep | Anonymous-read of entity metadata + annotations + child-listing for projects/folders; depth-capped recursive `walk`; full-text `search`; PAT-authenticated (`SYNAPSE_AUTH_TOKEN`) file download via the `fileHandle` → pre-signed-URL flow for controlled-access deposits (PsychENCODE, AMP-AD/PD, ROSMAP, BrainSpan). Data stays under upstream consortium DUAs — we only fetch with the user's own token, never redistribute. |
| **Open4Gene** peak→gene linkage (`open4gene link`) | [hbliu/Open4Gene](https://github.com/hbliu/Open4Gene) (Liu et al. *Science* 2025, PMID 39913582; **no LICENSE**) | clean-room Python reimpl — no source copied; upstream is R/`pscl::hurdle` | Two-component hurdle model per peak-gene pair: logistic zero component `I(RNA>0) ~ ATAC + covariates` + zero-truncated negative-binomial count component `RNA\|RNA>0 ~ ATAC + covariates`, via statsmodels `Logit` + `TruncatedLFNegativeBinomialP`; per-cell-type / All / Each modes; Spearman; AIC/BIC. **Validated vs the R `pscl::hurdle` reference: zero-component β correlation 1.0, max abs Δ 0.0.** |
| **scEPS** GWAS × single-cell neighborhood d-statistic (`sceps estimate`) | [Genentech/sceps](https://github.com/Genentech/sceps) (Zou/Shi et al. medRxiv 2026; **no LICENSE**) | clean-room Python reimpl — no source copied | Random-walk NAM neighborhood diffusion, per-donor pseudobulk, method-of-moments variance-component model decomposing disease variance into GWAS-gene / mean-expression-matched-control / rest components; per-neighborhood d-statistic (OMEGA_GWAS − OMEGA_CONTROL) with bootstrap disattenuation + delta-method SEs. **Validated vs upstream `test/` fixtures: step size, GWAS-gene count, neighborhood sizes, num-donors, expression variances all match exactly.** |
| **Functional-assay calibration → ACMG/AMP evidence** (`calibrate thresholds / prepare / run / assign / selftest`) | [rosstewart/exCALIBR](https://github.com/rosstewart/exCALIBR) (R. Stewart, Northeastern; MIT) — implements Zeiberg et al. *bioRxiv* 2025.04.29.651326 | clean-room Python reimpl — no source copied; stdlib + numpy/scipy (joblib and the SLURM job-array generator replaced by `concurrent.futures` + a resumable ledger) | Multi-sample skew-normal mixture (components shared across samples, per-sample mixing weights) fitted by EM in Azzalini's (loc, Δ, Γ) parameterisation with truncated-normal moments; monotone density-ratio constraint between adjacent components enforced by binary search on every parameter update; per-sample bootstrap with best-of-N fit selection on held-out likelihood; Saerens-style EM estimate of the population prior; LR⁺ envelope across bootstraps; Tavtigian C = O_PVSt search and C^(points/8) evidence thresholds; LR⁺ → per-strength score ranges with monotonicity repair; paired Wilcoxon / 5th-percentile 2c-vs-3c model selection; ClinVar-star + gnomAD + SpliceAI variant labelling of IGVF/Pillar-format scoresets. **Validated vs upstream: Tavtigian C identical over a 12-prior grid (plus `original` / `strict` variants); constrained-EM iterates bit-identical for 60 steps at K=2 and K=3; prior EM, point-range conversion and model-selection statistics identical to machine precision; scoreset labelling reproduces upstream's variant-ID sets exactly.** |
| **Integrated pathway databases** (`pathwaydb pull / build / query / status / sources`) | [IntPath](https://link.springer.com/article/10.1186/1752-0509-6-S2-S2) (Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. *BMC Syst Biol* 2012; **method, not code**) + [KEGG REST](https://rest.kegg.jp) + [Reactome downloads](https://reactome.org/download-data) + [WikiPathways GMT](https://data.wikipathways.org) | IntPath: method absorbed from the paper only — the published 2012 dataset is not redistributed. KEGG (academic use), Reactome (CC-BY), WikiPathways (CC0) fetched at runtime, never vendored. | IntPath's integration method applied to **current** releases: normalise every gene identifier into one namespace (here NCBI `gene_info`, Entrez → symbol with synonyms), map each database's own relation vocabulary onto one shared set (KEGG's PPrel / ECrel / GPrel / GErel), merge pathways that different databases describe under different names, and record **every** supporting database per fact so cross-database agreement stays queryable instead of collapsing into a duplicate. Reimplements IntPath's LCS name-unification rule (alignment ratio 2xLCS/(\|a\|+\|b\|), the two published acceptance conditions, error-prone word-pair filter, disjoint-set grouping, shortest name as the group name) with a bit-parallel LCS verified against the dynamic program. Extends it where the 2012 rule breaks on modern data: Reactome's own parent/child hierarchy (9,853 pairs) constrains unification so curated fact outranks string similarity, and candidates are corroborated by word overlap and gene membership before merging — measured at zero known-wrong merges vs 3.4% for the published rule alone (`pathwaydb evaluate`). Also adds typed relations parsed from current KEGG KGML maps, Reactome's curated interactor file, per-file release stamps + SHA-256, and a merge audit trail. The 2012 release covered 582 pathways; a current pull yields 4,016 pathways / 211k memberships / 74k typed relations, all loaded into the local knowledge graph. BioCyc (IntPath's third source) is not fetched — its download now requires a subscription — but a licensed local export integrates via `--extra-gmt`. |
| **figshare** data retrieval (`figshare article / files / download / search`) | [figshare API v2](https://docs.figshare.com) (Zenodo-style research-data deposit) | clean-room, urllib + json only | Resolve an article from numeric id, DOI, article URL, or private `/s/<token>` share link; list files (size + md5); md5-verified downloads (single file or whole article); public full-text article search. The general-purpose counterpart to the `synapse` skill for author-deposited supplementary data. |

#### Conventions for this table

**Licence hygiene.** IGVFagent is Apache-2.0. Where an upstream is GPL or
AGPL — BEAN, CORNETO, pySCENIC, edgeR, runHiC, TrimGalore — **no source is
copied**: the method is read and reimplemented, or the tool is invoked as a
separate program, and the *License* column says which. Vendoring copyleft
code would place its licence over this repository.

**Citations only where verified.** Where a formal citation appears it was
checked; otherwise the canonical repository or documentation URL is given,
which is what a reader needs to check the method. No author-year or DOI
string here was reconstructed from memory.

**Standard methods invoked by name** — Benjamini–Hochberg FDR, the adjusted
Rand index, limma's variance-moderation *idea* — are the standard procedure
or the idea at its simplest, not the upstream estimator. The docstrings say
so at each use, because "limma's variance moderation" unqualified would
overstate what `Scripts/_stats.py` does.

**Alternatives deliberately not implemented.** For pooled screens,
count-based negative-binomial models (MAGeCK-style, and BEAN's own `run`)
have materially more power than the per-guide tail-enrichment tests here,
because they model counting noise across guides rather than testing each
guide independently. That is measurable rather than theoretical: on
`IGVFDS6464SOVZ` none of 1,640 positive controls reaches FDR 0.05 under
per-guide testing at four replicates, while the top-ranked controls are the
biologically expected LDLR and HNF4A splice sites with the correct sign and
high editing activity. IGVFagent's contribution on such screens is correct
guide assignment, QC, and an honest statement of what its own test can and
cannot support — not a substitute for the upstream model.

### Methods papers cited in the skills

- **crispr-bean (BEAN)** — "Base Editing screens' Activity-Normalized variant effect size estimation", Pinello Lab. Repository <https://github.com/pinellolab/crispr-bean> (AGPL-3.0), documentation <https://pinellolab.github.io/crispr-bean/>. The base-edit-aware mapping method was read from `bean/mapping/GuideEditCounter.py`; no source is copied. *No paper citation is given here because none was verified from the repository — cite the upstream publication if you use this in a manuscript.*
- **Ma S et al. (2020)** "Chromatin potential identified by shared single-cell profiling of RNA and chromatin." *Cell* 183:1103–1116. doi:[10.1016/j.cell.2020.09.056](https://doi.org/10.1016/j.cell.2020.09.056) — SHARE-seq method.
- **Fulco CP et al. (2019)** "Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations." *Nature Genetics* 51:1664–1669. doi:[10.1038/s41588-019-0538-0](https://doi.org/10.1038/s41588-019-0538-0) — Flow-FISH log-normal bin-MLE method.
- **Nasser J et al. (2021)** "Genome-wide enhancer maps link risk variants to disease genes." *Nature* 593:238–243. doi:[10.1038/s41586-021-03446-x](https://doi.org/10.1038/s41586-021-03446-x) — Flow-FISH at scale + `Significant` / `Regulated` output convention.
- **Arnold CD et al. (2013)** "Genome-wide quantitative enhancer activity maps identified by STARR-seq." *Science* 339:1074–1077. doi:[10.1126/science.1232542](https://doi.org/10.1126/science.1232542) — STARR-seq method.
- **Tewhey R et al. (2016)** "Direct identification of hundreds of expression-modulating variants using a multiplexed reporter assay." *Cell* 165:1519–1529. doi:[10.1016/j.cell.2016.04.027](https://doi.org/10.1016/j.cell.2016.04.027) — Tewhey-lab MPRA method.
- **Smyth GK (2004)** "Linear models and empirical Bayes methods for assessing differential expression in microarray experiments." *Stat Appl Genet Mol Biol* 3:Article 3. doi:[10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027) — eBayes moderation reused in STARR-seq allelic test.
- **Law CW et al. (2014)** "voom: precision weights unlock linear model analysis tools for RNA-seq read counts." *Genome Biol* 15:R29. doi:[10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29) — voom mean-variance methodology.
- **Love MI et al. (2014)** "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biol* 15:550. doi:[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) — DESeq2 NB GLM in MPRA activity + via `pydeseq2`.
- **McGinnis CS et al. (2019)** "MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices." *Nat Methods* 16:619–626. doi:[10.1038/s41592-019-0433-8](https://doi.org/10.1038/s41592-019-0433-8).
- **Zhu Q et al. (2024)** "deMULTIplex2: robust sample demultiplexing for scRNA-seq." *Nat Methods*. — deMULTIplex2.
- **Subramanian A et al. (2005)** "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." *PNAS* 102:15545–15550. doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) — preranked GSEA backing `enrich gsea`.
- **Ashburner M et al. (2000)** "Gene ontology: tool for the unification of biology." *Nat Genet* 25:25–29. doi:[10.1038/75556](https://doi.org/10.1038/75556) — Gene Ontology Consortium; underlying ontology for `enrich go`.
- **Fabregat A et al. (2018)** "The Reactome Pathway Knowledgebase." *Nucleic Acids Res* 46:D649–D655. doi:[10.1093/nar/gkx1132](https://doi.org/10.1093/nar/gkx1132) — Reactome library used in `enrich pathways`.
- **Kanehisa M & Goto S (2000)** "KEGG: Kyoto Encyclopedia of Genes and Genomes." *Nucleic Acids Res* 28:27–30. doi:[10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27) — KEGG library used in `enrich pathways`.
- **Slenter DN et al. (2018)** "WikiPathways: a multifaceted pathway database bridging metabolomics to other omics research." *Nucleic Acids Res* 46:D661–D667. doi:[10.1093/nar/gkx1064](https://doi.org/10.1093/nar/gkx1064) — WikiPathways library used in `enrich pathways`.
- **Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L (2012)** "IntPath — an integrated pathway gene relationship database for model organisms and important pathogens." *BMC Systems Biology* 6(Suppl 2):S2. doi:[10.1186/1752-0509-6-S2-S2](https://doi.org/10.1186/1752-0509-6-S2-S2) — the cross-database pathway integration method (identifier normalisation, shared relation vocabulary, per-fact source provenance) that `pathwaydb` applies to current KEGG / Reactome / WikiPathways releases.
- **Liberzon A et al. (2015)** "The Molecular Signatures Database (MSigDB) hallmark gene set collection." *Cell Syst* 1:417–425. doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) — MSigDB Hallmark library used in `enrich pathways`.
- **Kuleshov MV et al. (2016)** "Enrichr: a comprehensive gene set enrichment analysis web server 2016 update." *Nucleic Acids Res* 44:W90–W97. doi:[10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377) — Enrichr proxy backing `enrich ora`.
- **Fang Z, Liu X, Peltz G (2023)** "GSEApy: a comprehensive package for performing gene set enrichment analysis in Python." *Bioinformatics* 39:btac757. doi:[10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757) — gseapy library powering `enrich`.
- **Zeiberg D, Tejura M, McEwen AE, Fayer S, Pejaver V, Rubin AF, Starita LM, Fowler DM, O'Donnell-Luria A, Radivojac P (2025)** "Gene-based calibration of high-throughput functional assays for clinical variant classification." *bioRxiv* 2025.04.29.651326. doi:[10.1101/2025.04.29.651326](https://doi.org/10.1101/2025.04.29.651326) — the calibration method behind `calibrate` (implemented upstream as exCALIBR).
- **Tavtigian SV et al. (2018)** "Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework." *Genetics in Medicine* 20:1054–1060. doi:[10.1038/gim.2017.210](https://doi.org/10.1038/gim.2017.210) — the C = O_PVSt constant and the C^(points/8) evidence-strength ladder used by `calibrate thresholds`.
- **Brnich SE et al. (2020)** "Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant interpretation framework." *Genome Medicine* 12:3. doi:[10.1186/s13073-019-0690-2](https://doi.org/10.1186/s13073-019-0690-2) — ClinGen SVI framework for converting assay odds-of-pathogenicity into PS3 / BS3 strengths.
- **Richards S et al. (2015)** "Standards and guidelines for the interpretation of sequence variants." *Genetics in Medicine* 17:405–424. doi:[10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30) — the ACMG/AMP rule set whose combining logic `calibrate` reproduces.

### License & attribution policy

- IGVFagent is **Apache-2.0** end-to-end (see [LICENSE](LICENSE)). We accept
  inbound code under MIT, BSD-2, BSD-3, ISC, and Python licenses; we do
  not redistribute GPL or AGPL source at runtime.
- For each absorbed pipeline, the skill source file's docstring names the
  reference repo, the upstream license, and a one-line summary of how the
  algorithm was paraphrased.
- If you spot a method we should attribute differently or a citation we
  missed, please open an issue or PR — we will fix it immediately.

## License

Licensed under the Apache License, Version 2.0 — see [LICENSE](LICENSE).
Copyright 2026 Hufeng Zhou.
