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

| Signal | Expected (paper) | This port's result |
|---|---|---|
| ≥ 1 CRISPR-KO dataset matching KMT2A in CD4+ T-cells | Weinstock 2024 deposited 84 KO screens | ✓ (discovery-only benchmark above) |
| Edges at \|β\|>0.020 / 0.025 / 0.030 | 350 / 211 / 151 (STAR Methods) | 1,294 / 877 / 618 — ~3.5-4x denser, see `llcb_py/README.md` |
| KMT2A connected to STAT5/JAK/IL2RA pathway | Paper Fig 5 | Not recovered — KMT2A connects to MED12/BCL11B/NFKB2/FOXP1 instead |
| rs45480496 (upstream of KMT2A) is a Th17-enhancer | Paper Fig 5 | Out of scope — needs external GWAS sumstats, not attempted |

## Full causal-network inference (`llcb_py/`)

A from-source Python port of the paper's own LLCB method
(github.com/weinstockj/LLCB, Julia — not R/Stan), fit on the paper's own
raw-count GEO deposit (GSE271788), not a KG analogue:

```bash
cd Benchmarks/weinstock2024_cd4_crispr/llcb_py
pip install pydeseq2   # one-time; everything else is already in this repo's env
python3 01_parse_counts.py
python3 02_normalize.py
python3 03_build_network_input.py
python3 04_fit_llcb.py           # writes Data/Weinstock2024/processed/edges.csv
python3 05_make_figures.py       # updates figures/fig4_*, figures/fig5_*
```

Runtime: well under a minute total on one CPU core. Requires network access
once, to fetch GSE271788's raw counts + series matrix from GEO's FTP and
the 84 gene-symbol→Ensembl mappings from mygene.info (both auto-fetched on
first run into `Data/Weinstock2024/raw/`).

See `llcb_py/README.md` for the full method-by-method comparison against
the paper (`llcb.py`'s module docstring explains every simplification and
one bug it caught and fixed along the way — an unconstrained noise-variance
estimate that silently degenerates on this problem's exactly-square linear
system).

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

For the LLCB causal-network step, there is no registered LLM tool — run
`llcb_py/`'s scripts directly via the shell.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `run.sh: .venv/bin/igvfagent: bad interpreter` | `.venv` was built on a different machine (e.g. synced from a laptop) | Already handled — `run.sh` detects this and falls back to `igvfagent` on `$PATH` |
| `crispr-screen_search.json` returns 0 hits | KMT2A not yet in catalogue | Try `--gene IL2RA` (also a paper-relevant gene) |
| `llcb_py/02_normalize.py` fails with `ModuleNotFoundError: pydeseq2` | not installed | `pip install pydeseq2` (one-time; not otherwise required by this repo) |
| `llcb_py/03_build_network_input.py` warns a gene wasn't found | mygene.info lookup missed a symbol, or GEO's gene ID version differs | Re-run; `01_parse_counts.py` already strips Ensembl version suffixes, so this should be rare |

## License + provenance

* **Paper data used by the discovery-only benchmark**: GEO **GSE171674**
  (public) — the CRISPR sub-series of SuperSeries GSE171737. This is
  *cited by* Weinstock 2024 as related context but is actually Freimer et
  al. 2022's own deposit (PMID 35817986); it is not Weinstock 2024's own
  bulk RNA-seq data.
* **Paper's own data, used by `llcb_py/`**: GEO **GSE271788** — Weinstock
  2024's own raw dedup UMI counts for all 84 KO'd genes across 3-4 donors
  (public, `pubmed_id: 39395408` on the GEO record itself).
* **Code**: IGVFagent Apache-2.0. `llcb_py/` is this session's clean-room
  Python port of the paper's own method.
* **Citation**: Weinstock JS, Arce MM, Freimer JW, Ota M, Marson A,
  Battle A, Pritchard JK. *Cell Genomics* **4**: 100671 (Nov 2024).
  doi:10.1016/j.xgen.2024.100671 · PMID 39395408 · PMC11605694
  *(Earlier benchmark scaffolds had the DOI mis-spelled as `100693`
  — fixed in a previous commit.)*
