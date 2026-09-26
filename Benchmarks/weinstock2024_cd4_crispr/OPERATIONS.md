# OPERATIONS — Weinstock 2024 CD4+ T cell CRISPR network

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run — fully online, ~30 s

```bash
bash Benchmarks/weinstock2024_cd4_crispr/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark weinstock2024_cd4_crispr
```

## What `run.sh` does

```bash
igvfagent perturb-catalog pipeline \
    --gene KMT2A --label weinstock2024_cd4_crispr --dataset-limit 50

igvfagent geo series --gse GSE171674
```

`run.sh` prefers `.venv/bin/igvfagent` if it works, else falls back to
`igvfagent` on `$PATH` — the checked-in `.venv` was built on a different
machine, so its shebang can be stale on another host.

1. `perturb-catalog pipeline` queries the Perturbation Catalogue's
   summary + all three modalities (mave, crispr-screen, perturb-seq)
   for KMT2A, and writes everything into one labelled run dir.
2. `geo series` pulls the Weinstock-specific GEO sub-series metadata.

The optional network follow-up (`network pkn-from-kg` / `steiner`,
printed at the end of `run.sh`) is not part of the scored benchmark —
see "Optional: full causal-network inference" below.

## Where artefacts land

```
Docs/Perturbation/<ts>_weinstock2024_cd4_crispr/
├── summary.json
├── global_search.json
├── mave_search.json
├── crispr-screen_search.json      ← primary_artefact in expected.json
├── perturb-seq_search.json
└── report.md

Docs/GEO/<ts>_GSE171674_geo_report.md
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
| `run.sh: .venv/bin/igvfagent: bad interpreter` | `.venv` was built on a different machine (e.g. synced from a laptop) | Already handled — `run.sh` detects this and falls back to `igvfagent` on `$PATH` |
| `network pkn-from-kg` errors with "no proteomics KG mirror" | KG mirror hasn't been pulled | Run `igvfagent kg-mirror pull --collection proteins_proteins` first |
| `crispr-screen_search.json` returns 0 hits | KMT2A not yet in catalogue | Try `--gene IL2RA` (also a paper-relevant gene) |
| PKN SIF has 0 edges | Wrong gene symbols (e.g. mouse vs human) | The proteomics KG is human-only; verify with `igvfagent catalog get-entity KMT2A` |

## License + provenance

* **Paper data**: GEO **GSE171674** (Weinstock's CRISPR sub-series of
  the joint Marson + Pritchard SuperSeries GSE171737) — public.
* **Code**: IGVFagent Apache-2.0; network skill is a clean-room MILP reimpl over `cvxpy` (no CORNETO GPL runtime dep).
* **Citation**: Weinstock JS, Arce MM, Freimer JW, Ota M, Marson A,
  Battle A, Pritchard JK. *Cell Genomics* **4**: 100671 (Nov 2024).
  doi:10.1016/j.xgen.2024.100671 · PMID 39395408 · PMC11605694
  *(Earlier benchmark scaffolds had the DOI mis-spelled as `100693`
  — fixed in this commit.)*
