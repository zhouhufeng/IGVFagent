# IGVF Agent

A local, auditable AI agent for discovering, retrieving, and analyzing data
from the [IGVF](https://igvf.org/) ecosystem (Portal, Catalog, Knowledge
Graph) and related public resources (ENCODE, FAVOR), with a built-in
**Plan → Action → Results → Evaluation** orchestration loop.

![IGVF Agent architecture](Docs/Figures/IGVF_agent_archetcture.png)

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

- [Capabilities](#capabilities)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Smoke test](#smoke-test)
- [Skill catalog and usage](#skill-catalog-and-usage)
  - [IGVF / ENCODE / Knowledge Graph client](#igvf--encode--knowledge-graph-client)
  - [Catalog / Portal / ENCODE overviews](#catalog--portal--encode-overviews)
  - [Variant annotation](#variant-annotation)
  - [Advanced variant analysis](#advanced-variant-analysis)
  - [Single-cell, multiome, specialized assays](#single-cell-multiome-specialized-assays)
  - [Parse SPLiT-seq pipeline](#parse-splitseq-pipeline)
  - [Enhancer–gene linkage](#enhancergene-linkage)
  - [MPRA / STARR / BlueSTARR](#mpra--starr--bluestarr)
  - [CRISPRi / CRISPR-FACS / Perturb-seq](#crispri--crispr-facs--perturb-seq)
  - [cCRE, FAVOR, IGV-style browser views](#ccre-favor-igv-style-browser-views)
  - [Data illustration and interpretation](#data-illustration-and-interpretation)
  - [Reference skill (literature retrieval, validation, design)](#reference-skill-literature-retrieval-validation-design)
  - [Knowledge Graph traversal](#knowledge-graph-traversal)
  - [Portal → local KG ETL](#portal--local-kg-etl)
  - [ENCODE bulk-genomics pipeline](#encode-bulk-genomics-pipeline)
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

## Repository layout

```
IGVFagent/
├── README.md                        ← you are here
├── requirements.txt                 ← optional analysis extras
├── .env.example                     ← copy to .env and edit locally
├── .gitignore                       ← excludes Data/, Logs/, .env, caches
├── Scripts/                         ← CLI skills + shared client code
│   ├── README.md                    ← full skill catalog
│   ├── igvf_client.py
│   ├── igvf_data_skills.py
│   ├── igvf_frontpage_summary.py
│   ├── annotate_variant_list.py
│   ├── advanced_variant_analysis.py ← integrated variant analysis pipeline
│   ├── ccre_linkage_annotation_skills.py
│   ├── enhancer_gene_linkage_skills.py
│   ├── mpra_data_skills.py
│   ├── crispri_data_skills.py
│   ├── single_cell_data_skills.py
│   ├── multiome_10x_pipeline.py
│   ├── multiome_research_demo.py
│   ├── igvf_specialized_data_skills.py
│   └── data_illustration_interpretation.py
├── Data/                            ← inputs + cached responses (gitignored)
│   └── Input/VariantList/           ← user-supplied variant CSVs
│       └── example_variants.csv     ← public smoke-test set
└── Docs/
    ├── PROJECT_SCOPE.md             ← what the agent is and isn't
    ├── DEPLOYMENT.md                ← legacy deployment doc (this README supersedes it)
    ├── IGVF_PORTAL_DATA_OVERVIEW.md
    ├── IGVF_CATALOG_SMOKE_ANALYSIS.md
    ├── IGVF_FRONT_PAGE_DATA_SUMMARY.md
    ├── ENCODE_DATA_OVERVIEW.md
    ├── ENCODE_SMOKE_ANALYSIS.md
    ├── CRISPRi/CRISPRI_BIOINFORMATICS_APPLICATIONS.md
    ├── Skills/                      ← per-skill reference docs
    │   ├── ADVANCED_VARIANT_ANALYSIS_SKILLS.md
    │   ├── IGVF_PORTAL_DATA_ANALYSIS.md
    │   ├── IGVF_FRONT_PAGE_SUMMARY_SKILLS.md
    │   ├── IGVF_SPECIALIZED_DATA_SKILLS.md
    │   ├── 10X_MULTIOME_SKILLS.md
    │   ├── SINGLE_CELL_ANALYSIS_SKILLS.md
    │   ├── ENHANCER_GENE_LINKAGE_SKILLS.md
    │   ├── MPRA_ANALYSIS_SKILLS.md
    │   ├── CRISPRI_ANALYSIS_SKILLS.md
    │   ├── CCRE_LINKAGE_FAVOR_SKILLS.md
    │   └── DATA_ILLUSTRATION_INTERPRETATION_SKILLS.md
    └── Logs/                        ← runtime logs (gitignored)
```

Generated outputs (timestamped folders under `Docs/<skill>/`, manifests under
`Data/Manifests/`, caches under `Data/Cache/`) are gitignored by default.

## Quick start

**Easiest (Docker Compose) — full stack including a local Qwen 3 LLM:**

```bash
git clone https://github.com/zhouhufeng/IGVFagent.git
cd IGVFagent
docker compose up -d                          # build + start agent + ollama
docker compose --profile bootstrap up         # one-time: pull qwen3:8b
open http://localhost:8501                    # browser UI
docker compose run --rm agent kg gene APOE    # one-shot skill from CLI
docker compose down                           # stop (volumes preserved)
```

The Compose stack runs IGVFagent and Ollama side-by-side in one network.
The agent reaches Ollama at `http://ollama:11434/v1`. `./Data` and
`./Docs` are mounted into the container so analyses persist across
restarts. Cloud LLM keys are forwarded from a local `.env` file (which
is gitignored) — the agent supports `ANTHROPIC_API_KEY`,
`OPENAI_API_KEY`, `GROQ_API_KEY`, `TOGETHER_API_KEY`, `DEEPINFRA_API_KEY`,
`HF_TOKEN`, plus the IGVF-specific `IGVF_PORTAL_COOKIE` and
`IGVF_ARANGO_PASSWORD` variables.

**Recommended (pipx) — one command, isolated, gives you `igvfagent` globally:**

```bash
pipx install 'git+https://github.com/zhouhufeng/IGVFagent.git'
pipx inject igvfagent 'igvfagent[analysis]'    # optional: scanpy + matplotlib stack
igvfagent --help
```

**Or pip into a virtual env (developer / contributor flow):**

```bash
git clone https://github.com/zhouhufeng/IGVFagent.git
cd IGVFagent
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -e '.[analysis]'        # core + analysis extras
# python3 -m pip install -e '.[all]'           # core + analysis + UI + LLM SDKs
igvfagent --help

# 3) configure endpoints / credentials locally
cp .env.example .env
# edit .env to set IGVF_ARANGO_PASSWORD, optional IGVF_PORTAL_COOKIE, etc.
# .env is gitignored.
```

After installing, the `igvfagent` console command gives you four ways to
drive the same skills:

```bash
# 1) Browser UI — chat input, streaming progress, inline plot rendering
pip install 'igvfagent[ui]'
igvfagent ui                      # opens http://localhost:8501

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
rendering of any artefacts the tools produce (UMAP / dot-plot / volcano
PNG/SVG figures, markdown reports, manifest CSVs).

The legacy `python3 Scripts/<skill>.py …` invocations documented later
in this README continue to work unchanged.

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
| Privacy | data never leaves your machine | prompts go to Anthropic | prompts go to Anthropic via Claude Code |
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
igvfagent ask --backend anthropic --model claude-sonnet-4-5 \
   "Compare DRD1 and DRD2 striatal MSN evidence in the local KG."

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

### CRISPRi / CRISPR-FACS / Perturb-seq

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
python3 Scripts/crispri_data_skills.py write-playbook
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

The agent runtime exposes `encode_retrieve`, `encode_describe`,
`encode_super_enhancers`, `encode_integrate_ccre`, and `encode_browser`
as tools, so a single `igvfagent ask` can drive the full pipeline:

```bash
igvfagent ask "Find K562 H3K27ac ChIP-seq experiments on GRCh38, pick \
  the highest-quality one, describe it, and call super-enhancers from \
  its peak BED. Then overlay the super-enhancers against SCREEN cCREs \
  and produce a browser view of the APOE locus."
```

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

## License

Licensed under the Apache License, Version 2.0 — see [LICENSE](LICENSE).
Copyright 2026 Hufeng Zhou.
