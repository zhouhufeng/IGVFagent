# Preparation Skills

How to prepare a fresh IGVFagent checkout for use. Every command in this
doc assumes you are sitting in the project root (`IGVFagent/`) and uses
**relative paths** — never hardcode absolute paths.

## 1. Prerequisites

| What | Why | Check |
|---|---|---|
| Python 3.10–3.12 | Runtime (Python 3.13 breaks numba/scanpy) | `python3 --version` |
| Git | Clone the repo | `git --version` |
| `claude` CLI (optional) | LLM backend = `claude_cli` (reuses your Claude Code login) | `claude --version` |
| ~5 GB free disk | Caches under `Data/Cache/`, manifests, reports | — |

## 2. Clone and create a virtual env

```bash
git clone https://github.com/zhouhufeng/IGVFagent.git
cd IGVFagent
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
```

Everything below uses `.venv/bin/igvfagent …` so you don't have to keep
`source`'ing the venv. If you prefer `source .venv/bin/activate`, the
plain `igvfagent` command works too.

## 3. Choose an install profile

The full `[analysis]` extra pulls in `numba/llvmlite`, which fails to
build on some macOS / Python combinations. Pick the lightest profile
that covers what you'll actually run:

### 3a. Minimum — LLM driver + core data skills (recommended first install)

```bash
.venv/bin/pip install -e '.[llm]'
.venv/bin/pip install matplotlib pandas numpy scipy
```

This is enough to drive every CLI skill that returns CSV/JSON manifests
and reports, plus the variant / MPRA / CRISPRi / enhancer / KG plotting
skills. It is **not** enough for single-cell pipelines (which need
scanpy / anndata).

### 3b. Single-cell stack (only when you need scRNA / scATAC / multiome)

```bash
.venv/bin/pip install -e '.[analysis]'
```

If this fails on `llvmlite`, you have a few options:

- Use Python 3.11 (best macOS-wheel coverage today).
- Skip numba: `pip install scanpy anndata --no-deps` then install the
  remaining deps individually, excluding `numba`.
- Run single-cell skills inside the Compose stack instead
  (`docker compose up -d agent`).

### 3c. Browser UI

```bash
.venv/bin/pip install -e '.[ui,llm]'
```

Adds Streamlit so you can `igvfagent ui`.

### 3d. Everything

```bash
.venv/bin/pip install -e '.[all]'
```

Equivalent to `analysis + ui + llm`. Has the same `llvmlite` caveat as
3b.

## 4. Configure credentials (optional)

Most public endpoints work with no credentials. You only need a `.env`
file if you want any of:

- **Authenticated Portal reads** (unreleased datasets, search of every
  type, KG queries): `IGVF_PORTAL_COOKIE`, `IGVF_ARANGO_PASSWORD`.
- **Cloud LLM backends**: `ANTHROPIC_API_KEY`, `OPENAI_API_KEY`,
  `GROQ_API_KEY`, `TOGETHER_API_KEY`, `DEEPINFRA_API_KEY`, `HF_TOKEN`.

```bash
cp .env.example .env
# edit .env locally — it is gitignored
```

`.env.example` documents every variable.

## 5. Pick an LLM backend

The agent's plan-act loop in `igvfagent ask` needs one of:

| Backend | Setup | When |
|---|---|---|
| `claude_cli` | `claude --version` already works | You use Claude Code anyway. No env var needed; the CLI subprocess reuses your existing login. |
| `anthropic` | `export ANTHROPIC_API_KEY=…` | You want best-in-class function calling and don't mind paying per query. |
| `ollama` | `ollama serve &` + `ollama pull qwen3:8b` | Offline / private / free. |
| `openai` | `export OPENAI_API_KEY=…` + `--backend openai` | You already have an OpenAI key. |

Resolution order: `--backend` flag > `IGVF_LLM_BACKEND` env > inferred
from `--model` name > default (`ollama` + `qwen3:8b`).

## 6. Verify

```bash
.venv/bin/igvfagent --version       # prints 0.1.0
.venv/bin/igvfagent backends        # lists registered providers
.venv/bin/igvfagent tools           # lists the agent's tool registry
.venv/bin/igvfagent client check    # ping Catalog docs + ENCODE
```

`igvfagent client check` is the canonical smoke ping. Expect:

- Catalog docs / API / `llms.txt` → HTTP 200.
- IGVF Portal home → HTTP 403 if `IGVF_PORTAL_COOKIE` is unset (this
  is **expected** and does not affect the public endpoints).

If `client check` returns 200/403 in that pattern, you're ready to run
the test suite — see [`TEST_SKILLS.md`](TEST_SKILLS.md).

## 7. Common pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `Failed to build … llvmlite` | numba wheels missing for your Python | Use Python 3.11, or install without `[analysis]` (§3a). |
| `matplotlib is required for plots` | Installed `[llm]` only | Add `pip install matplotlib pandas numpy scipy`. |
| `Portal returned 403` | No `IGVF_PORTAL_COOKIE` | Expected for public-only use. Set the cookie in `.env` if you need authenticated reads. |
| `Set IGVF_ARANGO_PASSWORD locally before using Knowledge Graph commands.` | KG AQL requires Arango credentials | Either skip `client aql`, or set the password in `.env`. The HTTP KG API (`igvfagent kg gene/variant/region`) works without it. |
| `Claude Code CLI (\`claude\`) is not on PATH.` | `claude_cli` backend chosen without Claude Code installed | Install Claude Code or pick a different backend. |

## 8. Driving the agent well (Karpathy-inspired guidelines)

IGVFagent's `ask` runner is LLM-driven, so the quality of an analysis
depends on how the ask is framed. The four guidelines below are
distilled from
[andrej-karpathy-skills](https://github.com/zhouhufeng/andrej-karpathy-skills)
and apply both to **users prompting the agent** and to **anyone editing
the IGVFagent codebase** (including LLM-driven edits).

These guidelines bias toward caution over speed. For trivial CLI
smokes (e.g. `client check`), use judgment and skip them.

### 8.1 Think before asking

Before you type the `ask`, articulate:

- **What entity?** Gene symbol, rsID, region, accession, or assay?
- **What evidence streams?** Catalog only, or also FAVOR / linkage /
  single-cell / literature?
- **What's the bound?** `--limit`, `--max-rows`, `--depth`.
- **What's success?** A CSV, a markdown report, a specific number?

If you're unsure, ask the agent to *clarify before pulling data*:

> "I want to study APOE in the striatum but I'm not sure whether
> single-cell datasets or MPRA are more relevant. List what's available
> first, then I'll decide."

Anti-pattern: a one-line request that hides three different
interpretations ("study APOE" → variants? expression? regulatory
elements? all three?).

### 8.2 Simplicity first

Send the smallest ask that produces the answer you actually need.

- Don't bundle 5 questions in one prompt — the planner will spawn 5 tool
  calls and you'll wait 5× longer for any signal.
- Don't ask for plots if you only need the CSV — `--label` + the
  generated `*_summary_stats.csv` is often enough.
- Don't pass huge `--limit` values for an exploratory pass; start at 5,
  scale up once you've seen the shape.
- If a single direct CLI skill answers the question, use it instead of
  `igvfagent ask` (which adds an LLM round trip).

Anti-pattern: `--depth 2 --limit 100 --call-favor --call-linkage
--call-singlecell --call-literature` on the first pass of an unfamiliar
gene.

### 8.3 Surgical changes (when editing IGVFagent itself)

If you (or an LLM acting on your behalf) modify the codebase:

- Touch only the lines tied to the stated change.
- Match the existing style (relative paths, Pathlib idioms, the
  `_endpoints.py` resolver pattern, the timestamped output folder
  convention).
- Don't reformat unrelated regions or "improve" adjacent code.
- Flag pre-existing dead code in a comment or message — don't delete
  it unprompted.
- Remove only the imports / helpers that *your* change orphaned.
- Every changed line should trace directly to the requested behavior.

Anti-pattern: a "bug fix" PR that also reformats the surrounding
function, renames variables, and changes the docstring.

### 8.4 Goal-driven execution

Turn a vague ask into a verifiable goal *inside the prompt*. The agent's
plan-act loop runs better when success has a checkable shape:

| Vague | Verifiable |
|---|---|
| "look at APOE" | "return the top-5 APOE coding variants from the IGVF Catalog with CADD > 20 as a CSV" |
| "find regulatory elements" | "list every cCRE in chr19:44,903,000-44,912,000 with class and Catalog support" |
| "make a plot" | "produce a volcano plot of `summary_stats.csv` with `log2fc` vs `-log10(padj)` and threshold lines at p < 0.05" |
| "test the pipeline" | "run `advanced-variant run` on `example_variants.csv` and confirm the report writes ≥1 SVG plot" |

For multi-step work, ask for a brief plan + verification checkpoint
before the run starts:

> "Plan the steps to take an APOE variant list through KG traversal →
> advanced-variant analysis → literature validation. After each step,
> tell me what artefact I should expect under `Docs/` before moving on."

Anti-pattern: "review and improve the splitseq pipeline" — no checkable
success criterion, no scope, no termination condition.

### 8.5 Why this matters here specifically

IGVFagent is a *research* agent: outputs feed downstream interpretation,
sometimes manuscript figures. The cost of a hallucinated dataset
accession or an over-eager `--limit 1000` traversal is real — wasted
analysis time, polluted caches, and audit trails that need to be
unwound. Tight prompts and surgical edits keep `Docs/<Skill>/` clean
and reproducible.
