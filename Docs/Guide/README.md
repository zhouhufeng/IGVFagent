# IGVFagent documentation

Start with the [project README](../../README.md) for what IGVFagent is and how
to try it. These pages hold the detail.

## Using IGVFagent

| Page | Read it when you want to… |
|---|---|
| [Installation and first run](installation.md) | install locally (pip, pipx or Docker), configure credentials, run the smoke test |
| [The browser UI, accounts and history](web-ui.md) | know what each tab does, how sign-in works, and how past results and projects are kept |
| [LLM backends](llm-backends.md) | choose a model (Claude, OpenAI, Ollama, Claude Code), keep results consistent across models, or drive IGVFagent from another agent |
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

---

[Project README](../../README.md)
