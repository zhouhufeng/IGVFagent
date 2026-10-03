# Installation and first run

Install IGVFagent locally, configure it, and check that it works. To use it without installing anything, see the [hosted instance](../../README.md#-try-it-online--no-install-required).

**On this page**

- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Smoke test](#smoke-test)
- [Workstation notes](#workstation-notes)

## Quick start

> Just want to try IGVFagent? Skip all of this and use the hosted instance at
> **<https://igvfagent.genohub.org>** — no install, no API key.
> See [Try it online](../../README.md#-try-it-online--no-install-required).

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

### Four ways to drive the same skills

After installing, the `igvfagent` console command gives you four ways to
drive the same skills. The browser UI is described in
[The browser UI](web-ui.md).

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

### Direct script calls

The legacy `python3 Scripts/<skill>.py …` invocations shown in the
[skill pages](skills/README.md) continue to work unchanged.

### Optional dependency groups

Declared in `pyproject.toml`:

| Extra | Adds | When you need it |
|---|---|---|
| `analysis` | pandas, numpy, scipy, matplotlib, pyarrow, anndata, scanpy, seaborn | Single-cell pipelines, advanced variant analysis, plotting |
| `ui` | streamlit | Browser UI (`igvfagent ui`) |
| `llm` | anthropic, openai | Native SDKs for the LLM-driven `igvfagent ask` runner |
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

To run every skill's offline self-test, as the
[tests workflow](../../.github/workflows/tests.yml) does on each push (network
blocked, no credentials, a few minutes):

```bash
python3 Scripts/run_selftests.py          # or --list to see what runs
```

Expected output locations:

- runtime logs in `Docs/Logs/`
- cached API responses in `Data/`
- manifests in `Data/Manifests/`
- reports in `Docs/<skill>/`

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

---

[← Documentation index](README.md) · [Project README](../../README.md)
