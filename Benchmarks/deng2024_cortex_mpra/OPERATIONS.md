# OPERATIONS — Deng 2024 cortex lentiMPRA + psychiatric variants

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

> **Online discovery works out of the box** (see "Quick run" below).
> The optional analytical step needs per-oligo DNA + RNA fastqs from
> the PsychENCODE Synapse deposit **`syn21392931`** (NeuREs / NOT
> GEO — the scaffold's original GSE236018 was wrong; that's a renal
> fibrosis paper, see `geo series --gse GSE236018` for confirmation).
> The new `igvfagent synapse` skill (Scripts/synapse_skill.py)
> reaches the deposit directly once you've accepted the PsychENCODE
> Data-Use Agreement at https://psychencode.synapse.org/DataAccess
> and exported `SYNAPSE_AUTH_TOKEN`.

## Required inputs

```
Data/Benchmarks/deng2024_cortex_mpra/cortex_mpra_counts.tsv   (required)
Data/Benchmarks/deng2024_cortex_mpra/top_targets.txt          (optional)
```

The optional `top_targets.txt` is a one-gene-per-line list of the most
strongly-activated MPRA target genes; if present, the run will also
fire `enrich ora` on it for the cortex enrichment spot-check.

## How to obtain the inputs (Synapse path, NOT GEO)

```bash
# 1. Accept the PsychENCODE Data-Use Agreement (free, registered account):
#    https://psychencode.synapse.org/DataAccess
# 2. Create a Personal Access Token (PAT) with view + download scopes at:
#    https://www.synapse.org/#!PersonalAccessTokens
# 3. Export the PAT in your shell:
export SYNAPSE_AUTH_TOKEN="eyJ0eXAiOiJKV1Qi..."

# 4. Discover the Deng 2024 deposit structure (anonymous-OK metadata):
.venv/bin/igvfagent synapse entity --syn syn21392931 --annotations
.venv/bin/igvfagent synapse walk --syn syn21392931 --max-depth 3
.venv/bin/igvfagent synapse children --syn syn51090452  # MPRA_CapstoneII

# 5. Identify the paired DNA + RNA fastq pairs you want (see the
#    children manifest). Files are organised as:
#      da-organoid-{dna,rna}-rep1..4_R{1_R3,2}.fastq.gz
#      da-primary-{dna,rna}-rep1..4_R{1_R3,2}.fastq.gz
#    plus per-donor bulk RNA: A2/A3/A5_*/S1_*/S2 _{1,2}.fq.gz

# 6. Download (per-file; PAT auto-applied):
mkdir -p Data/Benchmarks/deng2024_cortex_mpra
.venv/bin/igvfagent synapse download --syn synXXXX \
    --out-dir Data/Benchmarks/deng2024_cortex_mpra

# 7. Process the fastqs through MPRAflow / barcode-counter to get the
#    per-oligo count table cortex_mpra_counts.tsv (paper Methods step).
#    Place at the path below for run.sh's local-step branch:
mv <processed>.tsv Data/Benchmarks/deng2024_cortex_mpra/cortex_mpra_counts.tsv
```

## Quick run

```bash
bash Benchmarks/deng2024_cortex_mpra/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark deng2024_cortex_mpra
```

## What `run.sh` does

```bash
.venv/bin/igvfagent mpra activity --counts <INPUT> --label deng2024_cortex_mpra
.venv/bin/igvfagent mpra volcano  --label deng2024_cortex_mpra

# Optional, if top_targets.txt exists:
.venv/bin/igvfagent enrich ora --genes <TARGETS> --label deng2024_cortex_mpra_pathways
```

## Where artefacts land

```
Docs/MPRA/<ts>_deng2024_cortex_mpra/
├── mpra_activity.tsv     per-oligo log2FC + padj
├── volcano.png           4-panel volcano
├── volcano.svg
└── summary.json

# If enrich step ran:
Docs/Enrichment/<ts>_ora_deng2024_cortex_mpra_pathways/
├── enrichment.tsv
├── enrichment.png
└── summary.json
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `summary.json` exists, non-empty | mpra activity ran to completion |

(Extend after first run by adding a range check on
`n_active_oligos / n_total_oligos` ≈ 0.46.)

## Ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| Total elements tested | ~102,767 |
| Active enhancers called | ~46,802 (≈ 45.5 %) |
| Cell-type-specific patterns | dlEN (deep-layer excitatory) / ulEN (upper-layer) / MGE (medial-ganglionic-eminence interneurons) |
| Disease-variant allelic effects | 164 SCZ + bipolar GWAS variants alter activity |
| Enriched GO terms (target genes) | "neuron differentiation", "synaptic signalling", "axonogenesis" |

After the run:

```bash
# Compute active fraction
.venv/bin/python -c "
import csv
n_active = n_total = 0
for r in csv.DictReader(open('Docs/MPRA/<ts>_*/mpra_activity.tsv'), delimiter='\t'):
    n_total += 1
    try:
        if float(r.get('padj', 1)) < 0.05: n_active += 1
    except: pass
print(f'Active: {n_active:,} / {n_total:,} = {100*n_active/n_total:.1f}%')"
```

Paper target: ≈ 45.5 %.

## Running through the UI

`mpra_activity`, `mpra_volcano`, `enrich_ora` are all registered.
Paste:

```
Run the Deng 2024 cortex lentiMPRA reproducibility benchmark.

The per-oligo counts are at
Data/Benchmarks/deng2024_cortex_mpra/cortex_mpra_counts.tsv.

Step 1: mpra_activity with counts at that path, label="deng2024_cortex_mpra".
Step 2: mpra_volcano with that label.
Step 3: If Data/Benchmarks/deng2024_cortex_mpra/top_targets.txt exists,
        call enrich_ora with genes at that path and label
        "deng2024_cortex_mpra_pathways". Use libraries
        "GO_BP,GO_MF,GO_CC,Reactome_Pathways,KEGG_Human".

Report:
  - fraction of oligos active (paper target ~45%)
  - top 5 cell-type-specific enhancer hits if any cell-type column
    exists in the counts table
  - top 10 enriched pathways (if step 3 ran) — should include
    "neuron differentiation" and "synaptic signalling"
```

UI sidebar: max iterations = 18, temperature = 0.0.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Active fraction reports 0% or 100% | Counts table has wrong column names / negative-control labels | See `mpra activity --help` for expected schema |
| Enrichment step skipped | `top_targets.txt` not present (it's optional) | Create it from the volcano output's top-N hits, or omit the enrichment step |
| `enrich ora` reports 0 enriched terms | Gene symbols use ENSG instead of HGNC | Convert with `catalog get-entity` per gene, or use the `enrich go` subcommand which is more permissive |

## License + provenance

* **Paper data**: Synapse / PsychENCODE Knowledge Portal —
  `syn21392931` (NeuREs project, DOI `10.7303/syn21392931`).
  PsychENCODE Data-Use Agreement required for download; metadata
  + tree-walk are anonymously readable.
* **Code**: IGVFagent Apache-2.0; `synapse_skill.py` is a clean-room
  REST client (urllib only, no `synapseclient` dep).
* **Citation**: Deng C, Whalen S et al. *Science* **384**: eadh0559
  (2024). doi:10.1126/science.adh0559 · PMID:38781390 · PMC12085231
