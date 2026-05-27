# Deng 2024 — lentiMPRA in developing human cortex + psychiatric variants

## Citation

Deng C, ..., Pollard KS, Ahituv N. *Science* **384**: eadh0559 (2024).
DOI: [10.1126/science.adh0559](https://doi.org/10.1126/science.adh0559) · PMID: 38781390

## Data

* **GEO**: GSE236018 (primary cortex + organoid MPRA, PsychENCODE)
* **GitHub**: [Whalen-Lab/HumanCortexMPRA](https://github.com/Whalen-Lab/HumanCortexMPRA) — sei / CNN scoring

## Headline workflow (paper)

1. MPRAflow on 102,767 elements across primary cortex + cerebral organoids.
2. 46,802 active enhancers called.
3. 164 psychiatric-disorder GWAS variants alter MPRA activity (SCZ, BP).

## What IGVFagent reproduces

* **Local step**: `mpra activity` / `mpra volcano` on the downloaded DNA + RNA count table.
* `enrich ora` on the top-N target genes → expect "neuron differentiation" + "synaptic signalling" GO terms.

## Ground truth to spot-check

* Active fraction ≈ 46 % of elements (46,802 / 102,767)
* Top GO_BP enrichment should include "neuron differentiation" / "synaptic signalling"
