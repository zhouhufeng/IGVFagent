# OPERATIONS — Deng 2024 cortex lentiMPRA + psychiatric variants

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

> **Local input required.** Download per-oligo DNA + RNA counts from
> GSE236018 first; `run.sh` exits 77 with instructions until the file
> is at the expected path.

## Required inputs

```
Data/Benchmarks/deng2024_cortex_mpra/cortex_mpra_counts.tsv   (required)
Data/Benchmarks/deng2024_cortex_mpra/top_targets.txt          (optional)
```

The optional `top_targets.txt` is a one-gene-per-line list of the most
strongly-activated MPRA target genes; if present, the run will also
fire `enrich ora` on it for the cortex enrichment spot-check.

## How to obtain the inputs

```bash
mkdir -p Data/Benchmarks/deng2024_cortex_mpra
cd Data/Benchmarks/deng2024_cortex_mpra

# GSE236018 supplementary files
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236018
# Look for files named *_oligo_counts.tsv.gz or *_DNA_RNA_counts.tsv.gz

# Once downloaded:
gunzip GSE236018_cortex_oligo_counts.tsv.gz
mv GSE236018_cortex_oligo_counts.tsv cortex_mpra_counts.tsv

# Optional: pull the paper's "active enhancer" target-gene list
# (the paper's Supplementary Table 2 lists target genes for the
# most-active enhancers — copy into top_targets.txt)
cd ../../..
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

* **Paper data**: GEO GSE236018 + PsychENCODE — public.
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Deng C et al. *Science* **384**: eadh0559 (2024).
  doi:10.1126/science.adh0559 · PMID:38781390
