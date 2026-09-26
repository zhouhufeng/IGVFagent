---
name: igvf-replicate-paper
description: Reproduce a published paper's data and analyses with IGVFagent. Give it anything that identifies the paper — title, URL, DOI, PubMed ID, or author + journal + year — and it pins down the publication, reads its Data Availability statement, routes the deposits onto an IGVFagent analysis chain, and scaffolds a scored benchmark. Use when the user says "replicate this paper", "can IGVFagent reproduce X?", "benchmark this study", or pastes a citation/DOI/URL and asks what IGVFagent can do with it.
argument-hint: <title | URL | DOI | PMID | "author, journal, year">
---

# Replicate a paper with IGVFagent

Drives the `bench` skill end to end: publication → reproduction scaffold → run → concordance report. Everything here is a thin wrapper over `igvfagent bench …` (plus `paper-code` and `extauthor` when a route has to be built), rather than reimplementing their logic.

`$ARGUMENTS` is whatever the user knows about the paper. All of these work:

```
10.1038/s41588-024-01800-z
38969834
https://www.nature.com/articles/s41588-024-01800-z
https://www.biorxiv.org/content/10.1101/2023.11.09.563812v1
Saturation genome editing maps the functional spectrum of pathogenic VHL alleles
```

## Non-stop policy

This skill drives itself to a scored benchmark without pausing to ask the user, with exactly three exceptions — everything else that used to be a stopping point is now handled autonomously and recorded, not silently skipped:

| Situation | Old behaviour | Now |
|---|---|---|
| Paper resolution is `ambiguous` / `low_confidence` | Stop, show candidates, ask the user | **Continue** — auto-pick the top-scored candidate, log the alternates it did not pick (see step 1) |
| No IGVFagent route covers the paper's assay, but the paper has its own repo | Stop at a stub "discovery only" scaffold | **Continue** — reproduce the authors' code directly, and if that alone isn't enough, port it into a permanent IGVFagent tool (step 3b) |
| An unconfirmed (`[UNCONFIRMED]`) numeric check | Stop and hand it to the user to verify by hand | **Continue** — attempt to confirm it programmatically from the run's own artefacts (step 7); if that isn't possible, leave it unconfirmed and move on rather than blocking |
| Paper genuinely `not_found` after broadening the query | Stop, ask for a DOI | **Kept** — no identifier, no analysis; there is nothing to port or auto-pick. Report it and stop *this paper only*. |
| Deposit is controlled-access (dbGaP/EGA) or embargoed | Stop, ask the user | **Kept** — no amount of code-porting substitutes for a missing data-use agreement. Record it plainly in the report and continue with whatever is fetchable; never attempt to bypass access controls. |
| A check is `[UNCONFIRMED]` and cannot be self-verified | (see above) | **Kept in the score** — it stays `unreviewed`, never silently promoted to `ok`. A benchmark run finishing does not mean every claim was validated; the report says which weren't. |

Authoring new code (step 3b) is arbitrary code execution by design (`Docs/THREAT_MODEL.md` calls it "the largest single risk"), so it stays behind `IGVF_ALLOW_AGENT_AUTHORING`. This skill sets that variable **only for the subprocess calls it itself makes** (`export IGVF_ALLOW_AGENT_AUTHORING=1` in front of the specific `extauthor`/`paper-code` invocation, e.g. `IGVF_ALLOW_AGENT_AUTHORING=1 igvfagent extauthor write-skill …`) — it does not change the user's shell or any other session's default, and every write it makes is still the normal, inspectable, revertible extension file under `Scripts/promoted/` or the user extension directory (`igvfagent extauthor list` / `remove` afterward).

## Workflow

### 1. Pin down the paper — auto-resolve, don't wait on it

```bash
igvfagent bench resolve --query "$ARGUMENTS"
```

Add `--author`, `--journal`, `--year` when the user supplied them; they are scored as constraints and break ties.

The command prints a `decision`:

| decision | What to do |
|---|---|
| `resolved` | Continue to step 2. |
| `ambiguous` or `low_confidence` | **Do not stop.** Take `candidates[0]` (highest score). Lock it in deterministically: re-run `bench resolve --doi <candidates[0].doi>` (or `--pmid`). Record the alternates it did not pick — with their scores — in the eventual report/README so the choice is auditable, then continue to step 2. |
| `not_found` | Broaden once before giving up: retry with any `--author`/`--journal`/`--year` the user gave, and if the query was a bare title, try a short web/preprint search for its DOI and re-resolve with `--doi`. If still `not_found`, report that no source has it (or that it may be unpublished) and stop — this decision only, not the rest of a batch. |

Note the `paper_id` it derives — every later command takes `--paper-id`.

### 2. Read the paper's own data statement

```bash
igvfagent bench harvest --paper-id <paper-id>
```

Fetches full text (Europe PMC JATS for PMC open-access; the publisher page for bioRxiv/medRxiv preprints) and extracts the Data/Code Availability statements, repository accessions, assay families, gene symbols, and candidate numeric claims. Add `--no-llm` for deterministic extraction only.

**Check `Full text:` in the output.** If it says `UNAVAILABLE`, the paper is closed-access and everything downstream saw only the title and abstract — say so plainly to the user, because the scaffold will be thin and may not route at all. This is not a stop: continue with whatever the abstract-only harvest yields, and say clearly in the report that it's abstract-only.

**Check `Code:` in the output.** A `github_repo` accession marked `in_data_availability` means the paper names its own analysis repository — carry it into steps 3 and 3b below.

### 3. Route onto an IGVFagent chain

```bash
igvfagent bench route --paper-id <paper-id>
```

Ranks the 15 routes. Analysis routes that matched an assay always outrank the pure-retrieval fallbacks (Synapse / figshare / GEO / Portal). `igvfagent bench list-routes` shows the full table with the benchmark each route was modelled on.

If a route matched, continue to step 4. If nothing matched (`unroutable_reason` explains why — usually closed access, or an assay family outside IGVFagent's covered set), **do not settle for the retrieval-only stub** — go to step 3a.

### 3a. The paper has its own repo: run it directly

If harvest found a `github_repo` accession (marked `in_data_availability`), reproduce the authors' code before falling back to anything weaker:

```bash
igvfagent paper-code pipeline --harvest <run>/harvest.json
```

`paper-code pipeline` does not need the authoring env var on its own — it only builds a throwaway conda/pip env and runs the repo's own script unmodified. Reserve `IGVF_ALLOW_AGENT_AUTHORING` for step 3b, which registers something new.

This clones the repo at a pinned ref, builds its declared R/Python environment, runs its entry point(s) unmodified, and compares against the authors' rendered output where one exists. If it completes with a comparable output, that run *is* your reproduction — go to step 6 with `--route paper-code` semantics (score/report against this run directory; `bench scaffold` isn't needed when `paper-code` already produced a scored run).

If `paper-code pipeline` cannot complete — the repo has no discoverable entry point, the environment fails to build, a workflow step needs a tool IGVFagent doesn't have, or there is no repo named at all but the assay is still outside the 15 routes — go to step 3b. Do not stop and describe the paper as unreproducible; a stub scaffold is the last resort, not the first one reached.

### 3b. Never stop: port the missing piece into a permanent IGVFagent tool

This is the step that used to be "ask the user" or "ship a stub." Instead:

1. **Read the repo's actual analysis code** — the entry script(s) `paper-code inventory <repo>` identified, or the specific step that failed. Use `igvfagent paper-code exec <run_dir> --cmd ...` to poke at the built environment if it's unclear what a step needs.
2. **Write a faithful Python reimplementation** of that logic — not a stub, not a TODO, an actual working port, as a normal IGVFagent-style script (argparse `main()`, writes its artefacts under `Docs/`/`Data/` like every other skill here).
3. **Register it as a real, permanent core tool**, scoped-enabled for this call only:
   ```bash
   IGVF_ALLOW_AGENT_AUTHORING=1 igvfagent extauthor write-skill \
     --name <slug_for_the_method> \
     --description "<what it reproduces, and from which paper/repo>" \
     --source-file <path to the ported .py> \
     --tool-parameters '<JSON Schema for its CLI flags>'
   IGVF_ALLOW_AGENT_AUTHORING=1 igvfagent extauthor validate --name <slug_for_the_method>
   ```
   `write-skill` refuses code that doesn't parse or doesn't define `main()` — that's the safety rail that matters here, not a human sign-off gate. Prefer `extauthor write-tool --cli "paper-code pipeline …"` instead of `write-skill` whenever the repo's own code, run through `paper-code`, is already sufficient — write a *new* skill only for logic that genuinely isn't covered by running the repo directly (e.g. a comparison/aggregation step the paper's code never wrote for you).
4. **Re-route with the new tool** and continue:
   ```bash
   igvfagent bench route    --paper-id <paper-id>      # the new tool now shows up as a candidate route
   igvfagent bench scaffold --paper-id <paper-id> --route <slug_for_the_method>
   ```
5. Fall through to step 6. Note what was authored (name, source paper, files written) in the final report so it's auditable — `igvfagent extauthor list` shows everything currently registered, and it can be retired later with `extauthor remove` if it turns out to be wrong.

A tool authored this way stays in the core going forward: the next paper that needs the same method routes onto it directly, at step 3, with no porting needed.

### 4. Scaffold

```bash
igvfagent bench scaffold --paper-id <paper-id>
```

Writes `Benchmarks/<paper-id>/` with `run.sh`, `expected.json`, `README.md`, `OPERATIONS.md`, `provenance.json`, and registers the id in `Benchmarks/generated.txt`. Override the choice with `--route <name>` (the authored tool from 3b, if that's how you got here); overwrite an existing directory with `--force`.

Steps 1–4 in one go, when nothing needs porting:

```bash
igvfagent bench pipeline --query "$ARGUMENTS"
```

### 5. Read the scaffold, then keep going

Read the generated `run.sh` and `expected.json` and note, in the eventual report, rather than pausing on:

- **`TODO_VERIFY` variables.** The paper's text did not yield them. `run.sh` exits 77 until they are set — set them from the harvested text/provenance yourself where the value is genuinely there; only leave a `TODO_VERIFY` if no value exists anywhere in the harvested material.
- **Local-data requirements.** Controlled-access, embargoed, or R-only formats (`.rds`, `.qs`) cannot be fetched automatically — record this plainly (see the Non-stop policy table) and continue with whatever else the scaffold can run.
- **`[UNCONFIRMED]` checks.** Each has `"path": "TODO_SET_JSON_PATH"` and a `provenance.quote`. Try to confirm them now (step 7) rather than deferring — but a check that can't be confirmed is reported as `unreviewed`, not treated as a blocker.
- **`quote_grounded_in_source: false`** on an LLM claim means the model's quote is not verbatim in the harvested text. Treat that claim as unreliable in the report; don't stop for it.

### 6. Run and score

```bash
igvfagent bench run   --paper-id <paper-id>     # exit 77 = missing local input, not a failure
igvfagent bench score --paper-id <paper-id>
igvfagent bench report --paper-id <paper-id>
```

`report` renders the paper-claim vs IGVFagent-measured table into `Docs/Benchmark/<ts>_<paper-id>/replication_report.md`.

### 7. Confirm what can be confirmed — then finish, don't wait

A green `run.sh` proves the chain executed, not that the paper was reproduced. For each `[UNCONFIRMED]` check, attempt this yourself, right now:

1. Open the artefact named by `primary_artefact` in the run directory.
2. Find the key holding the comparable quantity; put its dotted path in `path`.
3. Verify the `provenance.quote` really states that number for that quantity.
4. If it checks out: set `"confirmed": true` and tighten `min`/`max` to a defensible tolerance, then re-score.
5. If it doesn't check out, or the quantity genuinely isn't in any artefact: leave it `"confirmed": false`. It stays `unreviewed` in the score — that is correct, not a failure to fix. Report the benchmark's true status (`ok` only where every scored check is confirmed; `unreviewed` or `partial` otherwise) and move on. Do not wait for a human pass before calling the run finished.

## What this can and cannot do

**Full analytical reproduction** works when the paper's assay is in a covered family (MAVE/SGE, MPRA/lentiMPRA, CRISPRi Flow-FISH, CRISPR screens, Perturb-seq, scRNA/multiome/SHARE-seq/SPLiT-seq, peak→gene, enhancer→gene, ChIP-Atlas, GWAS × single-cell), **or** when the paper names its own repository — step 3a/3b now reproduces that directly and, where it needed new logic, folds it into the core as a permanent tool rather than stopping.

**Discovery and retrieval only** — with a stub `run.sh` — is now the last resort, reached only when: the data is controlled-access (dbGaP/EGA) or embargoed, no repository is named anywhere in the paper, and no covered assay family matches either. The scaffold says which, in the README's "Honest caveats" section. Do not describe such a benchmark as a reproduction.

**Measured accuracy** on the 21 committed benchmarks (`igvfagent bench selftest --with-router`): resolver 21/22 exact from title alone; router 14/19 exact, 15/19 in the top 3. The misses are papers with no open-access full text — for those, expect to supply the route yourself with `--route`, or reach it through 3a/3b instead.

## Pairs well with

- `igvfagent ref learn --topic <assay>` — what other groups do for the same assay, before deciding what to reproduce.
- `igvfagent calibrate` — after a MaveDB/SGE route, turn the assay scores into ACMG/AMP PS3/BS3 evidence.
- `igvfagent extauthor list` / `igvfagent ext-review` — audit everything a run authored, across every paper, before it accumulates unreviewed.
