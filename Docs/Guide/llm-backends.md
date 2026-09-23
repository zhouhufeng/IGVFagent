# LLM backends

Which model drives the agent, how to pick one, how results are kept consistent across models, and how to drive IGVFagent from another agent.

**On this page**

- [Choosing an LLM backend](#choosing-an-llm-backend)
- [Consistent results across LLM backends](#consistent-results-across-llm-backends)
- [Driving IGVFagent from another agent](#driving-igvfagent-from-another-agent)

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
igvfagent ask --backend anthropic --model claude-opus-5-5 \
   "Compare DRD1 and DRD2 striatal MSN evidence in the local KG."
# Other current Claude ids: claude-opus-5-5, claude-sonnet-5, claude-fable-5-1, claude-haiku-4-5

export OPENAI_API_KEY=...
igvfagent ask --backend openai --model gpt-4o-mini "..."
```

Every run is persisted under `Docs/Agent/<timestamp>_<label>/` (transcript
JSON + markdown report + artefact paths from each tool call).

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

## Driving IGVFagent from another agent

Every skill is a shell command that writes its results to files and prints
their paths (`Report:`, `TSV:`, `Figure:` …), so any agent that can run
commands in the repository can use IGVFagent without its own orchestrator.

- **Claude Code.** Open the repository and ask in the chat
  (`cd IGVFagent && source .venv/bin/activate && claude`). It reads the skill
  playbooks in `Docs/Skills/` and calls `igvfagent <skill> …` directly. To have
  the in-process agent use Claude Code as its model instead, pick the
  **Claude Code CLI** backend (see the tables above).
- **Codex or another coding agent.** Give it the repository as the workspace
  and an instruction such as: *"Use the `igvfagent` CLI skills before writing
  new code. Keep data in `Data/` and reports in `Docs/`."*
- **Your own tool runner (Claude API, OpenAI, a local model).** Expose shell
  commands restricted to the repository folder. `igvfagent tools --json` lists
  every tool with its JSON schema and the command it maps to, which you can
  register as functions.

For example:

```bash
igvfagent client check
igvfagent processed lineage IGVFDS9875NBZW       # what the Portal already computed
igvfagent kg gene APOE --depth 2
igvfagent tools --json > tools.json              # schemas for your own runner
```

---

[← Documentation index](README.md) · [Project README](../../README.md)
