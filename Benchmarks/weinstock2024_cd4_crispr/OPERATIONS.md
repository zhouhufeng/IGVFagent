# OPERATIONS — Weinstock 2024 CD4+ T cell CRISPR network

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run — fully online, ~30 s

```bash
bash Benchmarks/weinstock2024_cd4_crispr/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark weinstock2024_cd4_crispr
```

## What `run.sh` does

```bash
.venv/bin/igvfagent perturb-catalog search-modality \
    --modality crispr-screen --query KMT2A \
    --label weinstock2024_cd4_crispr

.venv/bin/igvfagent network pkn-from-kg --label weinstock2024_cd4_crispr_pkn
```

Two unrelated tools, chained:

1. `perturb-catalog search-modality` queries for CRISPR-KO screens
   matching KMT2A in the Perturbation Catalogue.
2. `network pkn-from-kg` walks IGVFagent's local proteomics
   knowledge-graph mirror and materializes a signed Prior Knowledge
   Network (PKN) as a SIF file — this is the input substrate for a
   CARNIVAL / Steiner-tree run.

## Where artefacts land

```
Docs/Perturbation/<ts>_search_crispr-screen_KMT2A/
├── search.json
├── search.tsv
└── report.md

Docs/Network/<ts>_pkn_weinstock2024_cd4_crispr_pkn/
├── pkn.sif
├── pkn_summary.json
└── pkn_viz.svg
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `search.json` exists, non-empty | Perturb-Catalogue API returned results |

## Ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| ≥ 1 CRISPR-KO dataset matching KMT2A in CD4+ T-cells | Weinstock 2024 deposited 84 KO screens |
| PKN SIF contains KMT2A and connected to STAT5 / JAK pathway | Paper Fig 4 |
| 211 trans-edges in the LLCB causal network | Paper Table 1 |
| rs45480496 (upstream of KMT2A) is a Th17-enhancer | Paper Fig 5 |

After running:

```bash
# Check KMT2A neighborhood in the PKN
grep -i "kmt2a" Docs/Network/<ts>_*/pkn.sif | head -20
```

## Optional: full causal-network inference

For the complete Weinstock workflow (perturbation footprint → upstream
signalling subnetwork), after `run.sh` finishes:

```bash
# Build a perturbation-effect vector from Weinstock's 84-gene KO data
# (you'd download GSE171737 and compute per-gene log2FC)
# Then:
.venv/bin/igvfagent network carnival \
    --pkn Docs/Network/<ts>_pkn_weinstock2024_cd4_crispr_pkn/pkn.sif \
    --perturbations Data/Benchmarks/weinstock2024_cd4_crispr/perts.tsv \
    --label weinstock2024_carnival

# Or extract a Steiner-tree subnetwork connecting the 84 KO genes:
.venv/bin/igvfagent network steiner \
    --pkn Docs/Network/<ts>_pkn_weinstock2024_cd4_crispr_pkn/pkn.sif \
    --terminals KMT2A,STAT5A,STAT5B,IRF4,BATF,IL2RA \
    --label weinstock2024_steiner
```

Expected Steiner-tree result: KMT2A is connected through STAT5 / JAK
intermediates to IL2RA — recapitulating the paper's Th17-IL2 axis.

## Running through the UI

`perturb_catalog_search_modality` IS registered; `network_*` are not.

For the Perturbation-Catalogue step:

```
Run step 1 of the Weinstock 2024 reproducibility benchmark:
Call perturb_catalog_search_modality with modality="crispr-screen",
query="KMT2A", label="weinstock2024_cd4_crispr".

Then for the top result, call perturb_catalog_dataset with its id.

Report:
  - total CRISPR-KO datasets matching KMT2A
  - top 5 datasets (id, cell type, n_perturbations)
  - whether any are explicitly CD4+ T-cell context
```

For the network step, fall back to the shell (`network pkn-from-kg`
is not yet a registered LLM tool).

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `network pkn-from-kg` errors with "no proteomics KG mirror" | KG mirror hasn't been pulled | Run `igvfagent kg-mirror pull --collection proteins_proteins` first |
| `search.json` returns 0 hits | KMT2A not yet in catalogue | Try `--query IL2RA` (also a paper-relevant gene) |
| PKN SIF has 0 edges | Wrong gene symbols (e.g. mouse vs human) | The proteomics KG is human-only; verify with `igvfagent catalog get-entity KMT2A` |

## License + provenance

* **Paper data**: GEO GSE171737 — public.
* **Code**: IGVFagent Apache-2.0; network skill is a clean-room MILP reimpl over `cvxpy` (no CORNETO GPL runtime dep).
* **Citation**: Weinstock JS et al. *Cell Genomics* **4**: 100693 (2024).
  doi:10.1016/j.xgen.2024.100693
