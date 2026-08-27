# IGVFagent — data flow, trust boundaries, and threat model

This document exists because "local" and "reproducible" are each doing more
work in casual description than they can carry. It states precisely what
leaves a user's machine under each deployment mode, and separates three
properties that are routinely conflated:

| Property | What it means here | Does IGVFagent provide it? |
|---|---|---|
| **Auditability** | Every action is a typed, recorded, re-runnable command; inputs and outputs are on disk and inspectable | **Yes** — this is the system's strongest and most defensible property |
| **Repeatability** | Re-running the recorded commands produces the same artefacts | **Conditionally** — holds for pinned command sequences; does *not* hold for free-form agent planning, and never holds across upstream data changes |
| **Scientific validity** | The conclusion drawn is biologically correct | **Not established by the current evaluation** — see `Docs/EVALUATION.md` |

Auditability is a property of the *record*. Repeatability is a property of the
*execution*. Validity is a property of the *conclusion*. A system can have the
first without the second, and both without the third.

---

## 1. What "local" does and does not mean

IGVFagent runs its **tool execution** locally: skills are subprocesses on the
user's machine, and every file it reads, writes, caches, or logs stays inside
the project directory. That is real and worth stating.

It does **not** follow that data stays on the machine. Two separate egress
paths exist, and they are independent of each other:

**(a) LLM inference.** Unless the backend is a locally-served model, the prompt
leaves the machine. The prompt is not a bare question — it contains the system
prompt, the full tool catalogue, the conversation history, and **the stdout of
every tool call made so far**, which includes gene symbols, variant IDs,
accessions, file paths, and (since `read_artifact`) the contents of files the
agent has opened.

**(b) Upstream archives.** IGVFagent is a federated data-access system; that is
the point. Even with a fully local model, resolving a query contacts the IGVF
Catalog API, the IGVF Portal, ENCODE, FAVOR, Figshare, GEO, ChIP-Atlas, and
others. The gene you asked about is in those request URLs.

So a locally-served model gives you **local inference**, not **no egress**. The
only configuration with no network egress at all is a local model *and* skills
restricted to already-downloaded data.

---

## 2. Data flow by deployment mode

### Mode A — local install, locally-served model (Ollama / vLLM / TGI)

```
prompt ──▶ localhost model            (stays on machine)
skills ──▶ IGVF / ENCODE / GEO / …    (query terms leave)
files  ──▶ ./Data, ./Docs             (stay on machine)
```

Prompts and analysis artefacts never leave. Query terms do, to public
archives. This is the strongest privacy posture available.

### Mode B — local install, hosted API (Anthropic / OpenAI / Groq / …)

```
prompt + tool stdout + read files ──▶ model provider   (leaves)
skills ──────────────────────────▶ archives           (query terms leave)
files  ──────────────────────────▶ ./Data, ./Docs     (stay)
```

Artefacts stay on disk, but their *contents* can reach the provider whenever
the agent reads them back into context. Do not use this mode for controlled-
access or participant-level data unless your agreement with the provider
permits it.

### Mode C — hosted deployment (`igvfagent.genohub.org`)

Everything in Mode B, plus:

- **The workspace is shared.** `Scripts/_localstore.py` resolves a single
  process-global `IGVF_PROJECT_ROOT`, so every authenticated user reads and
  writes one `Data/`, one `Docs/`, and one local knowledge graph. Your uploads
  and results are visible to other users, and theirs to you.
- **Uploaded data leaves your machine twice** — once to the server, and again
  to the model provider when the agent reads it.
- **Authentication is a single shared password**, so "another user" means
  anyone who has been given it.
- **Extension authoring may be enabled** (`IGVF_ALLOW_AGENT_AUTHORING=1`),
  which lets the agent write Python that the server subsequently executes.

Mode C is a public demonstration sandbox. It is appropriate for exploring
public IGVF/ENCODE data and inappropriate for anything unpublished, embargoed,
participant-level, or otherwise sensitive.

---

## 3. Threat model

**Assets.** The operator's LLM API key; IGVF Portal credentials
(`IGVF_PORTAL_COOKIE`); infrastructure credentials; other users' analysis
outputs; the host itself.

**Adversaries considered.** An unauthenticated internet user; an authenticated
user of the shared deployment acting maliciously; a prompt-injection payload
arriving inside fetched archive content or an uploaded file.

### What is defended

| Control | Where | Against |
|---|---|---|
| No inbound port; outbound-only tunnel | `Deploy/docker-compose.prod.yml` | network reconnaissance |
| Shared-password gate ahead of the app | nginx | anonymous use of the operator's key |
| Artefact-path containment + secrets denylist | `Scripts/_pathguard.py` | reading `Docs/Secret/`, `.env`, keys — via the UI's renderers *and* `read_artifact` |
| Capability drop, `no-new-privileges`, PID and memory ceilings | app container | blast radius of subprocess execution |
| Pinned model, iteration and token caps | public mode | unbounded spend on the operator's key |
| Infrastructure credentials never synced to the host | deploy process | credential theft from the VM |

### What is *not* defended — stated plainly

1. **No per-user isolation.** One workspace, one knowledge graph. An
   authenticated user can read every other user's outputs. Fixing this
   requires making `IGVF_PROJECT_ROOT` per-session; it is currently
   module-level and resolved at import.
2. **Tool execution is subprocess execution.** `Scripts/_tools.py` runs wrapped
   CLIs with model-chosen arguments. Arguments are schema-bound, not free-form
   shell, but the class of risk is real.
3. **Extension authoring is arbitrary code execution by design.** When enabled,
   an authenticated user can have the agent write Python the host then runs.
   This is the intended feature; it is also the largest single risk, and it is
   why the flag defaults to off.
4. **Prompt injection is not mitigated.** Content fetched from archives, and
   files read via `read_artifact`, enter the model's context as trusted text.
   A crafted payload could steer subsequent tool calls.
5. **No rate limiting or per-user quotas.** Caps bound a single run, not the
   number of runs.

---

## 4. Reproducibility, stated precisely

The defensible claim is **auditability**: every agent action is a typed
command, recorded with its arguments, artefacts, and a consistency fingerprint
(system-prompt hash + seed + tool-set hash), and any session can be replayed
as a shell script.

Repeatability is narrower than that record suggests:

- **Pinned command sequences and playbooks are repeatable** — no LLM decides
  anything, so re-running reproduces the artefacts.
- **Free-form agent planning is not.** Measured artefact agreement across runs
  is well below 1.0; the router, seed, canonical tool set, and templated
  synthesis reduce variance without eliminating it.
- **Neither survives upstream change.** IGVF, ENCODE, and GEO records are
  updated and re-released. An identical command can return different data
  next month. Only the *command* is pinned; the *data* is not, unless a
  checksum-verified local copy is used.

Scientific validity is a separate question again, and the current benchmark
suite does not measure it — every `Benchmarks/*/run.sh` invokes skills
directly, with no model in the loop. That suite establishes skill correctness
only. See `Docs/EVALUATION.md`.
