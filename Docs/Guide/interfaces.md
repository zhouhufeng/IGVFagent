# Ways to use IGVFagent

One agent and one set of skills, reachable six ways: in the browser (hosted or on your own machine), as a chat in the terminal, as a one-line question, as direct commands, from Claude Desktop or an IDE over MCP, and from inside Claude Code or another coding agent.

**On this page**

- [At a glance](#at-a-glance)
- [Which one to use](#which-one-to-use)
- [Hosted: the web UI at igvfagent.genohub.org](#hosted-the-web-ui-at-igvfagentgenohuborg)
- [Local interfaces](#local-interfaces)
  - [Web UI on your machine](#web-ui-on-your-machine)
  - [Terminal chat](#terminal-chat)
  - [One question from the terminal](#one-question-from-the-terminal)
  - [Direct skill commands](#direct-skill-commands)
  - [MCP server: Claude Desktop, IDEs, other agents](#mcp-server-claude-desktop-ides-other-agents)
  - [Claude Code and other coding agents](#claude-code-and-other-coding-agents)
- [What every interface shares](#what-every-interface-shares)

## At a glance

| Interface | Hosted | Local | Who plans the steps | Needs an LLM key | Best for |
|---|:---:|:---:|---|---|---|
| **Web UI** | ✅ [igvfagent.genohub.org](https://igvfagent.genohub.org) | ✅ `igvfagent ui` | IGVFagent's own agent | hosted: no · local: yes, or Ollama | trying it, figures and tables inline, projects and history |
| **Terminal chat** | — | ✅ `bash Scripts/igvfagent_repl.sh` | IGVFagent's own agent | yes, or Ollama | asking question after question without a browser: over SSH, on a cluster |
| **One question** | — | ✅ `igvfagent ask "…"` | IGVFagent's own agent | yes, or Ollama | one answer, in a script or a batch job |
| **Direct skill commands** | — | ✅ `igvfagent <skill> <subcommand>` | you | no | pipelines, CI, exact reproducible runs |
| **MCP server** | — | ✅ `igvfagent mcp serve` | the client's model (Claude Desktop, an IDE, …) | the client's own | IGVF tools inside an assistant you already use |
| **Claude Code / other coding agents** | — | ✅ open the repository | the coding agent | the agent's own | day-to-day work in the repo; the bundled [Claude Code skills](#claude-code-and-other-coding-agents) |

Only the web UI is hosted. Everything else runs where IGVFagent is installed,
which is also where the skills execute and the files are written; see
[Installation and first run](installation.md).

## Which one to use

- **Just looking?** The [hosted web UI](#hosted-the-web-ui-at-igvfagentgenohuborg): no install, no API key.
- **Private or unpublished data, long pipelines, or a cluster?** Install locally. Then pick by how you like to work:
  - in a browser → [web UI on your machine](#web-ui-on-your-machine);
  - in a terminal, question after question → [terminal chat](#terminal-chat);
  - in a terminal, one question at a time or from a script → [`igvfagent ask`](#one-question-from-the-terminal);
  - exact, scriptable steps with no LLM at all → [direct skill commands](#direct-skill-commands).
- **Already live in Claude Desktop, Cursor or another MCP client?** Add the [MCP server](#mcp-server-claude-desktop-ides-other-agents) and ask there.
- **Working in the code with Claude Code or Codex?** [Let it drive the CLI](#claude-code-and-other-coding-agents), with the bundled skills for the common Catalog questions.

## Hosted: the web UI at igvfagent.genohub.org

**<https://igvfagent.genohub.org>** runs IGVFagent on the project's server with
the project's model key. Sign in with a
[Genohub community forum](https://discussion.genohub.org/) account; access
needs approval (email
[hufengzhou@g.harvard.edu](mailto:hufengzhou@g.harvard.edu)).

- Claude Opus 5.5 by default; Sonnet 5, Haiku 4.5 and Fable 5.1 selectable.
- Chat history is private to your account; the knowledge graph is shared and
  grows from everyone's work.
- Runs are capped per turn, and the server is one machine: do not upload
  unpublished or sensitive data.

The tabs, sign-in, projects and history are described in
[The browser UI, accounts and history](web-ui.md); hosted versus local side by
side is in the [README](../../README.md#-try-it-online--no-install-required).

## Local interfaces

All of these need a local install (`pip install -e '.[ui,llm,analysis]'`, see
[Installation](installation.md)). The ones that plan for themselves pick an
LLM automatically: a key in `.env` (`ANTHROPIC_API_KEY`, `OPENAI_API_KEY`, …)
if there is one, otherwise a local Ollama model; see
[LLM backends](llm-backends.md).

### Web UI on your machine

```bash
pip install -e '.[ui]'
igvfagent ui                      # opens http://127.0.0.1:8501
```

The same Streamlit app as the hosted site, with your own model and no run cap:
chat with a live progress trace, figures, tables, PDFs and reports rendered
inline, background jobs, the knowledge-graph explorer and the data viewers.
Tab by tab: [The browser UI](web-ui.md#the-browser-ui).

### Terminal chat

```bash
bash Scripts/igvfagent_repl.sh                                    # default backend and model
IGVF_LLM_BACKEND=ollama IGVF_LLM_MODEL=qwen3:8b bash Scripts/igvfagent_repl.sh   # local, offline
```

A prompt that keeps asking the agent in the terminal, printing each plan, tool
call and result as it goes. Useful where there is no browser: over SSH, inside
`tmux`, on an HPC login node. Each question is a fresh `igvfagent ask` run, so
the agent does not remember the previous question; the web UI keeps a
conversation. Commands inside it:

| Command | Does |
|---|---|
| `:model <name>` | switch model for the next questions |
| `:iter <n>` | cap the Plan → Action steps per question (default 12) |
| `:quiet` | hide or show the step-by-step trace |
| `:history` | list this session's transcripts |
| `:help` | the command list |
| `:q` | leave |

It sources `.env` from the repository and your home directory, so the API key
you set for `ask` works here too. The model is `IGVF_LLM_BACKEND` /
`IGVF_LLM_MODEL` when set; otherwise OpenAI `gpt-5` if there is an OpenAI key,
Claude Opus 5.5 if there is an Anthropic key, and `ask`'s own auto-detection
(a local Ollama model) if there is neither. It uses the repository's `.venv`
when that runs on this machine, then `igvfagent` on your `PATH`; set
`IGVFAGENT_BIN=/path/to/igvfagent` to choose one explicitly.

### One question from the terminal

```bash
igvfagent ask "What has IGVF already computed for IGVFDS9875NBZW?"
igvfagent ask --backend anthropic --model claude-opus-5-5 "Pull the comprehensive APOE evidence pack."
igvfagent ask --tool kg_gene --tool catalog_get_entity "Which enhancers regulate MYC in K562?"
```

The agent behind the web UI, run once: it plans, calls skills, and prints an
answer with the paths of every file it made. Options worth knowing:
`--backend` / `--model` to choose the LLM, `--max-iterations` to cap the loop,
`--tool` (repeatable) to restrict which tools it may use, `--quiet` for the
answer only. Because it is a single command it fits in shell scripts and
SLURM jobs.

### Direct skill commands

```bash
igvfagent --help                                  # every skill
igvfagent kg gene APOE --depth 2                  # one skill, no LLM involved
igvfagent portal-kg pull --tissue macrophage --limit 100
igvfagent tools                                   # the tool catalogue the agent sees
```

Every capability the agent has is also a plain command that writes files and
prints their paths. No model, no API key, nothing non-deterministic: the right
choice for pipelines, CI and runs you want to repeat exactly. Find the command
for a task in the [skill index](skills/README.md).

### MCP server: Claude Desktop, IDEs, other agents

```bash
igvfagent mcp list          # the tools that would be exposed (128 by default)
igvfagent mcp manifest      # a config snippet to paste into your client
igvfagent mcp serve         # run the server on stdio
igvfagent mcp serve --all   # expose every registered tool
```

Serves IGVFagent's tools over the [Model Context Protocol](https://modelcontextprotocol.io),
so an assistant you already use can call them. In Claude Desktop, add to
`claude_desktop_config.json`:

```json
{"mcpServers": {"igvfagent": {"command": "igvfagent", "args": ["mcp", "serve"]}}}
```

Any client that speaks MCP over stdio works the same way (Cursor, VS Code
agents, other agent frameworks). Here the client's model does the planning and
IGVFagent supplies the tools; a tool added to IGVFagent appears over MCP
automatically. The server runs on your machine, so the tools execute and write
their files there.

**IGVFagent installed on another machine** (a workstation or an HPC cluster)
can still serve a client on your laptop: stdio travels over SSH. Give the
client `ssh` as the command, with the full path to that machine's
`igvfagent`:

```json
{"mcpServers": {"igvfagent": {"command": "ssh",
  "args": ["-T", "you@cluster.example.org", "/path/to/envs/igvfagent/bin/igvfagent", "mcp", "serve"]}}}
```

The client cannot answer an SSH password or two-factor prompt, so set up key
login and, where the cluster asks for a one-time code, an SSH `ControlMaster`
connection that you open once in a terminal and the client then reuses. The
tools then run on that machine and write their files there.

To use it from Claude Code running on the same machine as IGVFagent, register
it once: `claude mcp add igvfagent -- /path/to/igvfagent mcp serve`.

### Claude Code and other coding agents

Open the repository in [Claude Code](https://claude.com/claude-code) and ask in
the chat: every skill is a shell command, so it can run them, read the files
they write, and chain them. The repository ships eight Claude Code skills
(in [`.claude/skills/`](../../.claude/skills/)) for questions that come up
often; invoke one by name, e.g. `/igvf-catalog-gene-dossier APOE`:

| Skill | For |
|---|---|
| `igvf-catalog-gene-dossier` | everything about one gene: coordinates, disease links, regulation, interactions, pathways |
| `igvf-catalog-variant-report` | interpreting one variant: associations, eQTL/sQTL effects, coding impact, LD |
| `igvf-catalog-dissect-locus` | from a GWAS hit to candidate causal genes |
| `igvf-catalog-regulatory-landscape` | how a gene is regulated: elements, QTLs, tissue specificity |
| `igvf-catalog-disease-genes` | the genetic architecture of a disease |
| `igvf-catalog-ld-compare` | LD for one variant across the five 1000 Genomes superpopulations |
| `igvf-portal-facet-filter` | narrowing IGVF Portal records facet by facet |
| `igvf-replicate-paper` | reproducing a published paper with IGVFagent and scoring it |

**IGVFagent's tools inside Claude Code, over MCP.** The repository's
[`.mcp.json`](../../.mcp.json) registers the MCP server for everyone who opens
it in Claude Code (Claude Code asks you to approve it once), provided
`igvfagent` is on your `PATH`. If it is not, register the full path for
yourself instead; a personal registration takes precedence:

```bash
claude mcp add igvfagent -- /path/to/envs/igvfagent/bin/igvfagent mcp serve
claude mcp list                   # igvfagent: … ✔ Connected
```

Start a new Claude Code session (servers load at start-up), type `/mcp` to see
the tools, then ask, e.g. *"Using the igvfagent tools, look up APOE in the
IGVF Catalog"*. Claude calls them as `mcp__igvfagent__catalog_get_entity` and
so on. Headless, for scripts: `claude -p "…" --allowedTools mcp__igvfagent`.

**Claude Code as IGVFagent's model.** The other way round: IGVFagent's own
agent (the web UI, `ask`, the terminal chat) plans with Claude Code, billed to
your Claude Code plan instead of an API key:

```bash
igvfagent ask --backend claude_cli "What is the Ensembl ID of APOE?"
IGVF_LLM_BACKEND=claude_cli bash Scripts/igvfagent_repl.sh
```

Codex, Ollama-driven agents and your own harnesses use the same CLI; recipes
per agent are in [Driving IGVFagent from another agent](llm-backends.md#driving-igvfagent-from-another-agent)
and backend trade-offs in [LLM backends](llm-backends.md).

## What every interface shares

- **The same skills.** A capability exists once, as an `igvfagent` command;
  every interface above reaches that same code.
- **The same outputs.** Runs write to `Docs/<Skill>/<timestamp>_<label>/` in
  the project folder, whichever interface started them, and the web UI's data
  viewers show them all.
- **The same growing knowledge graph.** Each run adds what it found to
  `Data/KG/local_kg.sqlite`, so an answer from the terminal today makes a
  later question in the browser faster
  ([how it grows](skills/knowledge-graph.md#the-growing-local-knowledge-graph--database)).
