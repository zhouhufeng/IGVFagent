# OPERATIONS — Joung 2025 TF Perturb-seq fibroblasts

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run — online catalogue query, ~5 s

```bash
bash Benchmarks/joung2025_tf_perturbseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark joung2025_tf_perturbseq
```

## What `run.sh` does

```bash
.venv/bin/igvfagent perturb-catalog search-modality \
    --modality perturb-seq \
    --query KLF4 \
    --label joung2025_tf_perturbseq
```

Hits the Perturbation Catalogue API (`api.perturbation-catalogue.org`)
and asks: "show me every Perturb-seq dataset where KLF4 was perturbed."
Joung 2025 deposited their TF Perturb-seq screens here, so we expect
hits.

## Where artefacts land

```
Docs/Perturbation/<ts>_search_perturb-seq_KLF4/
├── search.json
├── search.tsv
└── report.md
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `search.json` exists, non-empty | Perturb-Catalogue API returned results |

## Ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| At least 1 Perturb-seq dataset matching KLF4 | Joung 2025 perturbed 1,836 TFs including KLF4 + KLF5 |
| Fibroblast cell-type tag | Paper used dermal fibroblasts |
| KLF4 + KLF5 in same "universal fibroblast state" cluster | Paper Figure 3 |

Open the search.json and verify:

```bash
.venv/bin/python -c "
import json
d = json.load(open('Docs/Perturbation/<ts>_*/search.json'))
print(f'Total hits: {len(d.get(\"results\", []))}')
for r in d.get('results', [])[:5]:
    print(f\"  {r.get('id')}: {r.get('summary', r.get('title',''))[:80]}\")"
```

## Local analytical step (optional)

The full reproducibility test would:

1. Download SCP2169 from the Broad Single-Cell Portal (the paper's
   pre-publication-release distribution — GEO GSE237056 is **under
   embargo until 2027-12-31**, per `geo series --gse GSE237056`).
2. Run `sc-analyze pipeline` on the h5ad to recover the per-perturbation
   pseudobulk DE results.
3. Verify that KLF4 + KLF5 cluster together by trans-effect signature.

```bash
mkdir -p Data/Benchmarks/joung2025_tf_perturbseq
# SCP2169 (public, registered-account download):
#   https://singlecell.broadinstitute.org/single_cell/study/SCP2169
# Embargoed GEO mirror (release scheduled 2027-12-31):
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237056
mv ~/Downloads/joung2025.h5ad Data/Benchmarks/joung2025_tf_perturbseq/

.venv/bin/igvfagent sc-analyze pipeline \
    --input Data/Benchmarks/joung2025_tf_perturbseq/joung2025.h5ad \
    --label joung2025_tf_perturbseq

.venv/bin/igvfagent sc-analyze markers \
    --input Data/Benchmarks/joung2025_tf_perturbseq/joung2025.h5ad \
    --label joung2025_tf_perturbseq
```

Then run `enrich ora` on the KLF4-induced gene set:

```bash
.venv/bin/igvfagent enrich ora \
    --genes <klf4_top_DEGs.txt> \
    --label joung2025_klf4_enrichment
```

Expected top pathways: EMT / Hallmark Hypoxia / Hallmark KRAS Signaling
(reflecting fibroblast-state biology).

## Running through the UI

`perturb_catalog_search_modality` IS registered. Paste:

```
Run the Joung 2025 TF Perturb-seq reproducibility benchmark.

Step 1: Call perturb_catalog_search_modality with modality="perturb-seq",
query="KLF4", label="joung2025_tf_perturbseq".

Step 2: For the top result, call perturb_catalog_dataset with its id
to fetch full metadata.

Step 3: Call perturb_catalog_gsea on the result to see if there's a
precomputed GSEA endpoint.

Report:
  - total hits matching KLF4 in Perturb-seq modality
  - top 5 datasets (id, cell type, n_perturbations)
  - any GSEA pathways already attached to these datasets
```

UI sidebar: max iterations = 15, temperature = 0.0.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `search.json` is empty (no results) | KLF4 isn't yet in the catalogue; or the modality flag is wrong | Try `--query KLF5` or drop `--query` to enumerate all Perturb-seq datasets |
| Perturbation Catalogue API 503 | Service flap | Retry once; the catalogue is hosted on a single small instance |
| `sc-analyze pipeline` chokes on the h5ad | Anndata version mismatch | `pip install -U anndata`; or pass the h5ad through `anndata.read_h5ad` and re-save with current version |

## License + provenance

* **Paper data**: SCP2169 (Broad Single-Cell Portal — public) +
  GEO GSE237056 (embargoed until 2027-12-31 per NCBI's standard
  publication-mandated hold).
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Joung J et al. *Nat Genet* **57**: 828–838 (2025).
  doi:10.1038/s41588-025-02283-2
