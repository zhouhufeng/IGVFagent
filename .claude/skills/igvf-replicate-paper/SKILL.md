---
name: igvf-replicate-paper
description: Reproduce a published paper's data and analyses with IGVFagent. Give it anything that identifies the paper — title, URL, DOI, PubMed ID, or author + journal + year — and it pins down the publication, reads its Data Availability statement, routes the deposits onto an IGVFagent analysis chain, and scaffolds a scored benchmark. Use when the user says "replicate this paper", "can IGVFagent reproduce X?", "benchmark this study", or pastes a citation/DOI/URL and asks what IGVFagent can do with it.
argument-hint: <title | URL | DOI | PMID | "author, journal, year">
---

# Replicate a paper with IGVFagent

Drives the `bench` skill end to end: publication → reproduction scaffold → run → concordance report. Everything here is a thin wrapper over `igvfagent bench …`; run those commands rather than reimplementing their logic.

`$ARGUMENTS` is whatever the user knows about the paper. All of these work:

```
10.1038/s41588-024-01800-z
38969834
https://www.nature.com/articles/s41588-024-01800-z
https://www.biorxiv.org/content/10.1101/2023.11.09.563812v1
Saturation genome editing maps the functional spectrum of pathogenic VHL alleles
```

## Workflow

### 1. Pin down the paper — never skip this

```bash
igvfagent bench resolve --query "$ARGUMENTS"
```

Add `--author`, `--journal`, `--year` when the user supplied them; they are scored as constraints and break ties.

The command prints a `decision`:

| decision | What to do |
|---|---|
| `resolved` | Continue to step 2. |
| `ambiguous` | **Stop and show the user the candidate list.** Ask which one, or ask for a DOI/PMID. Do not guess. |
| `low_confidence` | Same — show candidates, ask for a stronger identifier. |
| `not_found` | Report that no source has it. Ask for a DOI, or whether it is unpublished. |

Note the `paper_id` it derives — every later command takes `--paper-id`.

### 2. Read the paper's own data statement

```bash
igvfagent bench harvest --paper-id <paper-id>
```

Fetches full text (Europe PMC JATS for PMC open-access; the publisher page for bioRxiv/medRxiv preprints) and extracts the Data/Code Availability statements, repository accessions, assay families, gene symbols, and candidate numeric claims. Add `--no-llm` for deterministic extraction only.

**Check `Full text:` in the output.** If it says `UNAVAILABLE`, the paper is closed-access and everything downstream saw only the title and abstract — say so plainly to the user, because the scaffold will be thin and may not route at all.

### 3. Route onto an IGVFagent chain

```bash
igvfagent bench route --paper-id <paper-id>
```

Ranks the 15 routes. Analysis routes that matched an assay always outrank the pure-retrieval fallbacks (Synapse / figshare / GEO / Portal). If nothing matches, the output says why — usually closed access, or a data type outside IGVFagent's covered assay families.

`igvfagent bench list-routes` shows the full table with the benchmark each route was modelled on.

### 4. Scaffold

```bash
igvfagent bench scaffold --paper-id <paper-id>
```

Writes `Benchmarks/<paper-id>/` with `run.sh`, `expected.json`, `README.md`, `OPERATIONS.md`, `provenance.json`, and registers the id in `Benchmarks/generated.txt`. Override the choice with `--route <name>`; overwrite an existing directory with `--force`.

Steps 1–4 in one go:

```bash
igvfagent bench pipeline --query "$ARGUMENTS"
```

### 5. Review before running — this is the human's job, not the tool's

Read the generated `run.sh` and `expected.json` and report to the user:

- **`TODO_VERIFY` variables.** The paper's text did not yield them. `run.sh` exits 77 until they are set via the named env var.
- **Local-data requirements.** Controlled-access, embargoed, or R-only formats (`.rds`, `.qs`) cannot be fetched automatically.
- **Unconfirmed checks.** Every check named `[UNCONFIRMED]` came out of prose via regex or the LLM. Each has `"path": "TODO_SET_JSON_PATH"` and a `provenance.quote`. They are **reported but never scored** — a benchmark with only unconfirmed checks scores `unreviewed`, not `ok`.
- **`quote_grounded_in_source: false`** on an LLM claim means the model's quote is not verbatim in the harvested text. Treat that claim as unreliable.

### 6. Run and score

```bash
igvfagent bench run   --paper-id <paper-id>     # exit 77 = missing local input, not a failure
igvfagent bench score --paper-id <paper-id>
igvfagent bench report --paper-id <paper-id>
```

`report` renders the paper-claim vs IGVFagent-measured table into `Docs/Benchmark/<ts>_<paper-id>/replication_report.md`.

### 7. Promote checks into a real reproducibility claim

A green `run.sh` proves the chain executed. It does **not** prove the paper was reproduced. To make that claim, for each unconfirmed check:

1. Open the artefact named by `primary_artefact` in the run directory.
2. Find the key holding the comparable quantity; put its dotted path in `path`.
3. Verify the `provenance.quote` really states that number for that quantity.
4. Set `"confirmed": true` and tighten `min`/`max` to a defensible tolerance.

Then re-score. Only now does a passing check mean something.

## What this can and cannot do

**Full analytical reproduction** works when the paper's assay is in a covered family: MAVE/SGE, MPRA/lentiMPRA, CRISPRi Flow-FISH, CRISPR screens, Perturb-seq, scRNA/multiome/SHARE-seq/SPLiT-seq, peak→gene, enhancer→gene, ChIP-Atlas, GWAS × single-cell.

**Discovery and retrieval only** — with a stub `run.sh` — when the data is controlled-access (dbGaP, EGA), embargoed, deposited in a format IGVFagent cannot read, or the assay is outside those families. The scaffold says which, in the README's "Honest caveats" section. Do not describe such a benchmark as a reproduction.

**Measured accuracy** on the 21 committed benchmarks (`igvfagent bench selftest --with-router`): resolver 21/22 exact from title alone; router 14/19 exact, 15/19 in the top 3. The misses are papers with no open-access full text — for those, expect to supply the route yourself with `--route`.

## Pairs well with

- `igvfagent ref learn --topic <assay>` — what other groups do for the same assay, before deciding what to reproduce.
- `igvfagent calibrate` — after a MaveDB/SGE route, turn the assay scores into ACMG/AMP PS3/BS3 evidence.
