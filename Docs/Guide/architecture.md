# Architecture

How the agent is put together, and where things live in the repository.

![IGVF Agent — system overview](../Figures/IGVFagent_system_overview.png)

_End-to-end view: a knowledge graph and multi-omics data resources feed an orchestration-layer AI agent (Plan → Act → Observe → Refine) backed by short/long-term memory, an execution layer, post-processing, and human-in-the-loop review — producing downstream outputs such as variant prioritization, perturb-seq analysis, enhancer–gene mapping, and fine-mapping._

![IGVF Agent — architecture and skill topology](../Figures/IGVF_agent_archetcture.png)

_Detailed five-layer architecture: user entry points (terminal, NL agent, browser UI) → agent runtime & tool dispatch → 121 skills / 613 typed tools grouped by domain → local persistence (filesystem + DuckDB warehouses) → upstream services. The `network` skill (highlighted) is the apex of the skill DAG — a clean-room MILP reimplementation of CORNETO that reads from the Silver + Bronze warehouses and writes inferred subnetworks back._

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

## Repository layout

```
IGVFagent/
├── README.md               ← project front page
├── pyproject.toml          ← installable package; Scripts/ is the `igvfagent` package
├── Dockerfile  docker-compose.yml  requirements.txt  LICENSE  .mcp.json
│
├── Scripts/                ← every skill and the runtime
│   ├── cli.py              ← `igvfagent <skill> <subcommand>` dispatcher (the SKILLS registry)
│   ├── _agent.py           ← Plan → Action → Results → Evaluation loop
│   ├── _llm.py             ← LLM backend router (Anthropic, OpenAI, Ollama, vLLM, Claude Code, …)
│   ├── _tools.py           ← the typed tool registry the agent sees
│   ├── _history.py         ← permanent run history and projects
│   ├── _localstore.py      ← the local knowledge graph that grows with use
│   ├── _pathguard.py       ← which files the UI may render
│   ├── streamlit_app.py    ← browser UI (`igvfagent ui`), with kg_sources.py and data_browser.py
│   ├── *_skill.py, *.py    ← one module per skill (see the skill index)
│   ├── ported/             ← paper code rewritten as commands (`igvfagent port`), with provenance
│   ├── promoted/           ← user extensions promoted to built-ins
│   ├── _paper/             ← .docx / .pptx export helpers
│   └── test_*.py           ← offline self-tests, one per skill
│
├── Benchmarks/             ← reproducibility suite (definitions; tracked)
│   ├── README.md           ← suite dashboard
│   ├── concordance.py      ← scorer (`igvfagent bench score`)
│   ├── <paper-id>/         ← expected.json, run.sh, README.md, figures/
│   ├── test_*.py           ← suite-level self-tests
│   ├── _data/              ← public datasets downloaded by run.sh (gitignored)
│   └── results/            ← per-run concordance reports (gitignored)
├── Deploy/                 ← hosted deployment: redeploy.sh, mirror-kg.sh, autosync.sh, gateway (auth/), nginx
├── Docs/
│   ├── Guide/              ← this documentation
│   ├── Skills/             ← per-skill playbooks read by the agent (Markdown)
│   ├── SkillCards/         ← machine-readable skill cards, generated by `igvfagent skillcard` (JSON)
│   ├── Architecture/       ← design references (integration layer, data model)
│   ├── Examples/  Playbooks/  Demos/  Customize/   ← worked examples, YAML playbooks, demo questions
│   ├── Figures/            ← diagrams and logos
│   ├── BugZilla/           ← external test reports
│   ├── Eval/               ← backend reproducibility evaluations
│   ├── upstream.json       ← pinned commit of every absorbed upstream repository
│   ├── AUTH.md  EVALUATION.md  THREAT_MODEL.md  PROJECT_SCOPE.md
│   ├── <Skill>/<timestamp>_<label>/   ← run outputs (gitignored); a few kept as
│   │                                    tracked worked examples: MultiSeq/, BiosampleCensus/,
│   │                                    scE2G/, Network/
│   ├── Benchmark/<timestamp>_<paper-id>/  ← benchmark run output: resolve, harvest, route,
│   │                                    report, port verification, paper-specific runs (gitignored)
│   ├── PaperCode/          ← runs of authors' own code (gitignored)
│   ├── Logs/               ← per-run log files (gitignored)
│   └── Secret/             ← local credentials (gitignored; never commit)
├── Data/                   ← downloads, caches, KG and warehouse files (gitignored)
│   ├── KG/                 ← the local knowledge graph (local_kg.sqlite)
│   ├── Warehouse/          ← Catalog mirror (KG/ Parquet) and DuckDB views
│   ├── Cache/  Manifests/  Jobs/  History/  Uploads/
│   └── <paper>/            ← per-paper downloads used by Benchmarks/<paper-id>/run.sh
└── UserExtensions/         ← optional per-checkout custom tools and skills (absent by default)
```

The authoritative list of skills is the `SKILLS` registry in
[`Scripts/cli.py`](../../Scripts/cli.py), shown with descriptions in the
[skill index](skills/README.md#all-skills). `Scripts/README.md`
lists direct `python3 Scripts/<module>.py` invocations.

Generated outputs are gitignored: timestamped run folders under `Docs/<Skill>/`
(an empty `Docs/<Skill>/` folder only means that skill has not run on this checkout),
benchmark runs under `Docs/Benchmark/`,
manifests under `Data/Manifests/`, the knowledge graph under `Data/KG/`, the
Catalog mirror under `Data/Warehouse/KG/`, and caches under `Data/Cache/`.
A new built-in skill ships as its module in `Scripts/`, a `Scripts/test_*.py`
self-test, a CLI entry and tools, a playbook in `Docs/Skills/`, and, for a
port, a pin in `Docs/upstream.json`. User extensions need none of that; see
[Extending IGVFagent](extending.md).

---

[← Documentation index](README.md) · [Project README](../../README.md)
