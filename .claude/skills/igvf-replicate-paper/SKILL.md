---
name: igvf-replicate-paper
description: Reproduce a published paper's data and analyses with IGVFagent. Give it anything that identifies the paper — title, URL, DOI, PubMed ID, or author + journal + year — and it pins down the publication, lists every headline analysis, runs the authors' own code as the reference, rewrites whatever IGVFagent doesn't already cover as verified Python commands absorbed into IGVFagent, and keeps going until every analysis is reproduced. Use when the user says "replicate this paper", "can IGVFagent reproduce X?", "benchmark this study", or pastes a citation/DOI/URL and asks what IGVFagent can do with it.
argument-hint: <title | URL | DOI | PMID | "author, journal, year">
---

# Replicate a paper with IGVFagent

`$ARGUMENTS` is whatever the user knows about the paper (DOI, PMID, URL, bioRxiv link, or title).

## How this differs from Paper2Agent, and what "done" means

Paper2Agent wraps a paper's repository and runs the authors' scripts as they are. IGVFagent does that too, as the **reference**, and then goes further: every analysis IGVFagent does not already have a command for is **rewritten as a Python module, verified against the authors' own output, and registered into IGVFagent** (`igvfagent port register`) so this paper, and every later one, runs it as `igvfagent <name>`. The reproduction then continues on those new commands.

**Done is paper coverage, not a passing check tally.** `igvfagent bench score` prints two things:

- `status` (`ok` / `partial` / `fail`): whether the declared checks passed. One easy count can make this `ok`.
- `reproduction`: `reproduced` only when **every planned analysis** has a passing, confirmed class-A (quantitative vs the paper) or class-B (port vs the authors' code) check tied to it. Counts and file-exists checks never cover an analysis.

Keep working until `reproduction` is `reproduced`, or `reproduced_except_access` (every remaining analysis needs controlled-access, embargoed or undeposited data). `ok` with `incomplete k/N` means **keep going**. Verifying a number on the smallest deposited file is a sanity check, never the reproduction.

## Non-stop policy

Do not pause to ask the user. The only reasons to stop:

| Situation | What to do |
|---|---|
| Paper genuinely `not_found` after broadening the query | Report it; stop this paper only. |
| An analysis needs controlled-access (dbGaP/EGA), embargoed, or undeposited data | `bench plan --block <id> --kind controlled_access\|embargoed\|not_deposited --reason <evidence>`; continue with every other analysis. Never try to bypass access controls. |

Everything else is work, not a stop: an ambiguous paper match (take the top candidate and log the alternates), code in R/Julia/MATLAB/Nextflow, an environment that won't build, a tool IGVFagent lacks, hard-coded paths, multi-GB data (submit with `sbatch`), a port that disagrees with the reference (diagnose and fix, up to 6 attempts per port; then record it as failed and continue).

## The loop

### 1. Resolve, harvest, scaffold

```bash
igvfagent bench resolve --query "$ARGUMENTS"          # ambiguous → re-run with --doi <candidates[0].doi>
igvfagent bench harvest --paper-id <paper-id>
igvfagent bench route   --paper-id <paper-id>
igvfagent bench scaffold --paper-id <paper-id> [--route <name>]
```

If harvest says `Full text: UNAVAILABLE`, find the open preprint (bioRxiv/PMC) and use its text; say so in the report. Record the paper's own repository and data deposits: the Code Availability statement, and if it has none, the authors' GitHub README and the Zenodo/figshare/GEO records. The automated harvester's accessions include third-party tools the methods cite, so check which ones are the paper's own.

### 2. Plan by coverage

```bash
igvfagent bench plan --paper-id <paper-id> --seed
```

`--seed` is a draft from the harvest. Complete it against the paper's **figure list** and the **authors' repository**: one analysis per headline figure panel or claim. Drop Methods-only items such as library cloning or cell counts. For each analysis, record the authors' file(s) that produce it:

```bash
igvfagent bench plan --paper-id <paper-id> --add fig3b_tf_clusters \
    --title "TF perturbations cluster into fibroblast states" --figure "Fig. 3b" \
    --claim "<the paper's sentence>" --value <number> --upstream notebooks/fig3.ipynb
igvfagent bench plan --paper-id <paper-id> --remove <seeded-id-that-is-not-a-result>
```

### 3. Run the authors' code, to get the reference

```bash
igvfagent paper-code pipeline --harvest <run>/harvest.json
```

This pins the repository, builds its environment, and runs it unmodified. Its outputs are the **reference** each port is checked against. Large inputs go through `sbatch` (use the cluster env's `igvfagent`; never run heavy work on a login node). If a step cannot run here after at most 3 environment repairs, don't stop: its reference becomes the authors' deposited intermediate output or published table (Source Data, supplementary tables), and step 4 ports it.

### 4. For each analysis: reuse, or port and absorb

First, reuse. An existing command (`igvfagent --help`, `bench list-routes`) or an earlier port may already do it:

```bash
igvfagent port find --query "<what the analysis computes>"
```

Otherwise port it:

1. Read the authors' code for exactly this analysis (the `--upstream` files) and write a **faithful Python port**. Keep the same algorithm, defaults, filters, normalisation and random seeds. It is an argparse `main(argv=None)` module that writes its artefacts under `Docs/` or `Data/`. No invented science: where their code is ambiguous, follow the paper's Methods and note the choice.
2. Register it into IGVFagent. Set the authoring switch only for this call:
   ```bash
   IGVF_ALLOW_AGENT_AUTHORING=1 igvfagent port register --name <method_slug> \
       --paper-id <paper-id> --repo <repo URL> --commit <pinned SHA> \
       --upstream <repo path> [--upstream ...] --source-file <port.py> \
       --description "<what it computes> (port of <repo>/<file>)" \
       --tool-parameters '<JSON Schema of its flags>' --keywords "<terms>"
   ```
   It is then `igvfagent <method-slug>` on the CLI and a tool in the registry, recorded in `Scripts/ported/registry.json` (paper, repo, commit, upstream files) and marked unreviewed until a human looks at it.
3. Verify it on the **same input** the reference came from:
   ```bash
   igvfagent bench verify-port --paper-id <paper-id> --name <method_slug> --analysis <id> \
       --reference <authors' output> --port-output <port's output> \
       [--key a.b | --column <col> --id-column <id>] [--rtol 1e-6] [--min-corr 0.99]
   ```
   It writes `validation_vs_reference.json` and prints the check to add. Missing rows count as mismatches. On FAIL, read the `worst` rows, fix the port, re-register and re-verify, at most 6 attempts. Then record it as failed and move on.
4. Run the verified command on the paper's **full** data (sbatch when large), and add to `expected.json` checks tied to the analysis with `"analysis": "<id>"`:
   - the class-B `verify-port` check that `verify-port` printed;
   - and, where the paper states a number, a class-A `range` check on the port's output against that number (tolerance from the paper's own precision or seed noise).

### 5. Score, and loop

```bash
igvfagent bench run    --paper-id <paper-id>
igvfagent bench score  --paper-id <paper-id>      # read `reproduction: … k/N`
```

While `reproduction` is `incomplete`, go back to step 4 for the next analysis in state `pending`, `weak` (only counts pass) or `failing`. Do not report the paper as reproduced before this reaches `reproduced` or `reproduced_except_access`.

### 6. Report

```bash
igvfagent bench report --paper-id <paper-id>
igvfagent port list
```

The report's coverage table lists every analysis and its state. In the final message, state:
- the coverage verdict (`k/N`);
- each port absorbed into IGVFagent (name, upstream file, `verify-port` match rate);
- where the reproduction differs from the paper, and why;
- each blocked analysis, with its evidence.

Don't claim more than `bench score` says.

## Pairs well with

- `igvfagent ref learn --topic <assay>`: what other groups do for the same assay.
- `igvfagent calibrate`: after a MaveDB/SGE route, turn assay scores into ACMG/AMP PS3/BS3 evidence.
- `igvfagent ext-review`: audit agent-authored extensions; ported methods are listed by `igvfagent port list` with their provenance and verification.
