# Synapse retrieval skill

Clean-room IGVFagent integration with Sage Bionetworks Synapse
(`https://www.synapse.org`). Pure urllib + json, no `synapseclient`
Python dep required.

## Subcommands

| Command | Purpose |
|---|---|
| `igvfagent synapse entity --syn synXXX` | Fetch entity metadata + annotations. |
| `igvfagent synapse children --syn synXXX` | Paginated child-listing of a project / folder. |
| `igvfagent synapse walk --syn synXXX --max-depth 3` | Recursive tree walk (capped). |
| `igvfagent synapse search --query "<text>"` | Full-text search across public Synapse. |
| `igvfagent synapse download --syn synXXX --out-dir <dir>` | Stream a file entity to disk. |

## Authentication

* For public entities (most IGVF + most upstream-genomics deposits)
  no auth is required — the REST API answers anonymously.
* For **controlled-access** entities (PsychENCODE, AMP-AD, AMP-PD,
  several IGVF-controlled cohorts) you need a **Personal Access
  Token** (PAT). Create one at
  https://www.synapse.org/#!PersonalAccessTokens after accepting the
  relevant Data-Use-Agreement and grant the `view` and `download`
  scopes.
* Export the token in your shell:

      export SYNAPSE_AUTH_TOKEN="eyJ0eXAiOiJKV1Qi..."

  Every `igvfagent synapse ...` call will then automatically include
  `Authorization: Bearer <token>` headers.

## Why this skill

Several IGVF-relevant repositories live on Synapse rather than GEO /
ENCODE Portal:

* **PsychENCODE Consortium** — Deng 2024 *Science* cortex lentiMPRA
  (`syn21392931` — folder name `NeuREs`); cross-references with the
  IGVF cortex multiome lines.
* **AMP-AD + ROSMAP + MSBB** — Alzheimer's cohorts heavily cross-linked
  to IGVF brain assays.
* **AMP-PD** — Parkinson's; relevant to IGVF Corces/Gladstone
  multiome AnalysisSets.
* **BrainSpan v2** — developmental references for the IGVF brain map.
* Several **IGVF-controlled** AnalysisSets (donor-consent restricted)
  are mirrored to Synapse alongside the IGVF Portal.

## Example: Deng 2024 cortex lentiMPRA

    igvfagent synapse entity --syn syn21392931            # NeuREs root folder
    igvfagent synapse children --syn syn21392931          # Data + PEC sub-folders
    igvfagent synapse walk --syn syn21392931 --max-depth 3
        # Enumerate the full PsychENCODE NeuREs project (cap depth=3).

    # After accepting the PsychENCODE Data-Use Agreement:
    export SYNAPSE_AUTH_TOKEN="eyJ0eXAiOi..."
    igvfagent synapse download --syn <file-syn-id> \
        --out-dir Data/Benchmarks/deng2024_cortex_mpra/
