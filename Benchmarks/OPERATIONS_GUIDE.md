# Benchmark Suite — Common Operations Guide

This guide collects everything that's true for *every* benchmark in
`Benchmarks/`. For per-paper specifics (data download URLs, ground-truth
spot-checks, troubleshooting), see each paper's `OPERATIONS.md`.

## 0. Prerequisites — one-time setup

| Requirement | How to verify |
|---|---|
| You are at the repo root | `pwd` returns `…/IGVFagent` |
| The `.venv/` is built | `.venv/bin/igvfagent --version` returns a version string |
| The console script is wired | `which .venv/bin/igvfagent` resolves |
| Network access to public APIs | `curl -sI https://api.data.igvf.org/profiles/MeasurementSet.json | head -1` returns `HTTP/2 200` |
| 1 GB free disk for outputs | `df -h .` shows ≥ 1 GB free |

For benchmarks that exercise an LLM (UI / `igvfagent ask`), additionally:

| Requirement | How to verify |
|---|---|
| LLM backend configured | `.venv/bin/igvfagent backends` lists at least one with `api_key_env` set |
| For local Ollama | `curl -s http://localhost:11434/api/tags` returns a JSON model list |

## 1. The two execution paths

Every benchmark has two ways to run:

### Path A — Shell script (recommended for the reproducibility claim)

The canonical, deterministic path. Bypasses the LLM entirely. The exact
`igvfagent …` invocations live in `Benchmarks/<paper-id>/run.sh`. The
same commands run every time → byte-equivalent outputs on the same
inputs.

```bash
bash Benchmarks/<paper-id>/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark <paper-id>
```

### Path B — Streamlit UI chat (recommended for exploration)

The LLM-driven path. Boots from `igvfagent ui`, available at
`http://localhost:8501`. Useful for poking at a benchmark interactively
or asking follow-up questions ("OK, now show me the top 10 active
elements by p-value"). Each chat session is captured under
`Docs/Agent/<ts>_<query>/transcript.json` for provenance.

Recommended UI sidebar settings for benchmark questions:

| Setting | Value | Reason |
|---|---|---|
| Backend kind | Local | No cloud API key needed for the bundled local models |
| Model | `qwen3.6:35b-a3b-coding-bf16` | Strongest local model with usable tool-use |
| Max iterations | 12–25 | Bumped for multi-step benchmarks; raise for #14/#15/#28-style synthesis questions |
| Max tokens / turn | 4096 | 8192 if reports get truncated |
| Temperature | 0.0 | Reproducibility + correct tool-arg construction |
| Tool subset | empty | Let the LLM pick |

## 2. Where artefacts land

IGVFagent has two output-layout conventions:

### Convention 1 — Per-run directory (newer skills)

`Docs/<SkillName>/<ts>_<label>/` containing `summary.json`, TSVs, plots.
Used by: `mavedb`, `multiome`, `chipatlas`, `portal`, `catalog`,
`enrich`, `sc-analyze`, `flowfish`.

```
Docs/MaveDB/20260527_175908_matreyek2018_pten_vampseq/
├── summary.json
├── PTEN_mapped.tsv
├── PTEN_mapped.vcf
├── mapping_summary.png
├── mapping_summary.svg
└── showcase_report.md
```

### Convention 2 — Flat-file at the skill root (older skills)

Files named `<ts>_<label>_<purpose>.<ext>` directly under
`Docs/<SkillName>/` (and sometimes `Data/Manifests/<SkillName>/` and
`Data/`). Used by: `mpra pull`, `encode retrieve`, parts of `crispri`.

```
Docs/MPRA/20260527_191820_agarwal2025_lentimpra_portal_many_report.md
Data/Manifests/MPRA/20260527_191820_agarwal2025_lentimpra_portal_many_manifest.csv
Data/20260527_191820_mpra_agarwal2025_lentimpra_portal_many_summary.json
```

`concordance.py` supports both via `extra_search_dirs` in `expected.json`.

### Per-agent-run provenance

When you drive a benchmark through the UI / `igvfagent ask`, an extra
trace lands at:

```
Docs/Agent/<ts>_<query>/
├── transcript.json    Every tool call + every LLM response
└── report.md          User-facing summary
```

This is the auditable trail that distinguishes IGVFagent from a
chat-PDF agent.

## 3. Concordance scoring

```bash
.venv/bin/python Benchmarks/concordance.py --benchmark <paper-id>
.venv/bin/python Benchmarks/concordance.py --all
```

The scorer reads `Benchmarks/<paper-id>/expected.json`, walks the
latest run dir (per the two-convention rule above), and applies four
check types:

| Type | Meaning |
|---|---|
| `range` | A numeric value at a dotted path in the JSON artefact must fall in `[min, max]` |
| `in_set` | A string value must be in the `allowed` set |
| `artefact` | A file (or glob) must exist and be non-empty |
| `row_count_tsv` | A TSV file's row count must fall in `[min, max]` |

### Confirmed vs unconfirmed checks

A check may additionally carry `"confirmed": false`. Benchmarks scaffolded by
`igvfagent bench` use this for every number extracted from a paper's prose —
by regex, or by the optional LLM headline-claim pass. Such a check ships with
`"path": "TODO_SET_JSON_PATH"` and a `provenance` block quoting the sentence
it came from.

**Unconfirmed checks are reported but never scored.** They do not count toward
`n_passed` / `n_total`, and a benchmark whose checks are *all* unconfirmed gets
status `unreviewed` (icon ⊘) rather than `ok` — a run that has verified nothing
must never look green. The suite-level Markdown lists them under each paper with
their source quote.

To promote one into a real reproducibility claim:

1. Run the chain once and open the artefact named by `primary_artefact`.
2. Find the key holding the comparable quantity; put its dotted path in `path`.
3. Check the quoted sentence really states that number for that quantity.
4. Set `"confirmed": true` and tighten `min` / `max` to a defensible tolerance.

Results land in `Benchmarks/results/<ts>_concordance.{json,md}`
(`results/` is gitignored — regenerated on each scoring run).

A clean `ok` status means every declared check passed; `partial` means
some passed; `unreviewed` means only unconfirmed checks exist so nothing
was actually scored; `no_run_found` means the scorer couldn't locate the
run directory (usually: the `run.sh` wasn't executed yet, or the label in
`expected.json` doesn't match what `run.sh` actually wrote).

A few skills (`chipatlas`, `perturb-catalog`, `portal get`) accept no
`--label`, so their run directories are not paper-tagged. Benchmarks on
those routes set `label` to the skill's own default directory token and
carry a `label_note` saying so — scoring will pick up the most recent run
of that skill, whichever paper produced it.

## 4. Running through the LLM backend

For UI-driven benchmark execution, the LLM needs **registered tools**
(the JSON-schema-typed entries in `Scripts/_tools.py`). Not every CLI
subcommand is registered. Coverage at the time of this guide:

| Skill | Registered tool count | Examples |
|---|---:|---|
| `portal` | 13 | search, get, schema, facets, report, batch-download, … |
| `catalog` | 6 | get-entity, search-region, find-associations, find-ld, resolve-id, list-sources |
| `chipatlas` | 8 | list-antigens, search, get-experiment, assemble-bed, target-genes, … |
| `perturb-catalog` | 7 | summary, search, search-modality, dataset, dataset-rows, gsea, pipeline |
| `flowfish` | 5 | estimate-effects, real-space, score-elements, simulate, … |
| `mpra` | 5 | pull, activity, skew, qc, volcano |
| `mavedb` | 2 | map-scoreset, showcase |
| `enrich` | 5 | ora, gsea, go, pathways, showcase |

If your prompt needs a CLI subcommand that isn't registered (e.g.
`mpra portal-manifest`, `multiome retrieve`, `encode retrieve`,
`flowfish pull-portal`), one of:

1. Use the shell-script path (always works).
2. Adapt the prompt to request a *registered* tool that does the same
   discovery work, e.g. ask for `portal_search type=MeasurementSet`
   instead of `mpra_portal_manifest`.
3. Add a tool registration to `Scripts/_tools.py` following the existing
   pattern.

## 5. Common troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `run.sh` exits 77 with "Local input not found" | The local-data step needs a download | Follow the per-paper OPERATIONS.md "Data download" section |
| `concordance.py` reports `no_run_found` | The label in `expected.json` doesn't match what `run.sh` wrote, OR the run hasn't been executed | Re-run `run.sh` first; check `Docs/<skill>/` for the label |
| `concordance.py` reports `partial` with `glob ... matched no non-empty files` | The skill scattered artefacts and the scorer's search dirs are incomplete | Add the missing directory to `extra_search_dirs` in `expected.json` |
| "URN not found" from MaveDB | The scoreset URN in `run.sh` is unverified | See per-paper OPERATIONS.md — most have a `TODO_VERIFY` note + an env-var override |
| Agent UI says "stop max_iterations" | Either real loop, or one of the 10 positional-arg tools mis-dispatched | Already fixed in commit 3155dbc; re-pull / restart streamlit |
| `kg variant <rsid>` returns all-zero relations | Old kg skill — fixed in commit 4e547ab via SPDI pre-resolution | Re-pull / restart |
| Ollama "model not found" | The default model name doesn't match what's installed | `igvfagent models` to list; pass `IGVF_LLM_MODEL=<name>` |

## 6. Reproducibility report

After running benchmarks, the canonical "we reproduced X" artefact is
the concordance Markdown at `Benchmarks/results/<ts>_concordance.md`.
This is what to attach to a paper review, share with a collaborator,
or commit to a public repo as evidence.

```bash
# After running everything:
bash Benchmarks/run_all.sh --online-only
.venv/bin/python Benchmarks/concordance.py --all

# Then share the report:
ls -t Benchmarks/results/*_concordance.md | head -1
```

The report includes:
- Suite summary (total checks passed / total)
- Per-paper status (✓ ok / △ partial / ✗ fail / — no_run_found)
- Per-check detail with the actual measured value vs the expected range

## 7. Conventions used throughout the per-paper guides

* Paths are POSIX-style, relative to the repo root unless otherwise noted.
* `<paper-id>` is the directory name under `Benchmarks/`, e.g. `agarwal2025_lentimpra`.
* `<ts>` is the timestamped prefix `YYYYMMDD_HHMMSS` IGVFagent writes on every artefact.
* `<label>` is what the run.sh passes via `--label`; it's always the `<paper-id>`.
* "Online" = no local input file required; the run.sh pulls everything from public REST APIs.
* "Local" = the run.sh checks for a file under `Data/Benchmarks/<paper-id>/` and exits 77 with download instructions if missing.

That's the shared scaffolding. Each per-paper OPERATIONS.md fills in
the specifics on top of this.
