# OPERATIONS — Zheng 2024 in-vivo AAV Perturb-seq, mouse cortex

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

> **Local input required.** This benchmark only has a meaningful run
> after you've downloaded GSE249416 (~30K cells, ~500 MB compressed).
> `run.sh` exits 77 cleanly until you place the h5ad at the expected
> path.

## Required input

```
Data/Benchmarks/zheng2024_invivo_perturbseq/GSE249416.h5ad
```

### How to obtain it

GEO doesn't deposit h5ad directly — it ships raw counts + metadata as
separate files. Convert them locally:

```bash
mkdir -p Data/Benchmarks/zheng2024_invivo_perturbseq
cd Data/Benchmarks/zheng2024_invivo_perturbseq

# Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249416
# Download:
#   GSE249416_RAW.tar      (per-sample 10x mtx bundles)
#   GSE249416_sgrna_assignments.tsv
#   GSE249416_metadata.tsv

tar -xf GSE249416_RAW.tar

# Then convert to a single h5ad with anndata:
.venv/bin/python <<EOF
import anndata as ad
import scanpy as sc

# Read each 10x sample, concatenate, add sgRNA assignments
adatas = []
for sample_dir in sorted(__import__('pathlib').Path('.').glob('GSM*')):
    a = sc.read_10x_mtx(sample_dir, var_names='gene_symbols')
    a.obs['sample'] = sample_dir.name
    adatas.append(a)

adata = ad.concat(adatas, join='outer', label='sample')
adata.obs.index = adata.obs.index + '_' + adata.obs['sample']

# Attach sgRNA assignments
import pandas as pd
sgrna = pd.read_csv('GSE249416_sgrna_assignments.tsv', sep='\t', index_col=0)
adata.obs = adata.obs.join(sgrna)

adata.write_h5ad('GSE249416.h5ad', compression='gzip')
print(f'Wrote {adata.shape[0]:,} cells × {adata.shape[1]:,} genes')
EOF

cd ../../..
```

## Quick run

```bash
bash Benchmarks/zheng2024_invivo_perturbseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark zheng2024_invivo_perturbseq
```

## What `run.sh` does

```bash
.venv/bin/igvfagent sc-analyze pipeline --input <INPUT> --label zheng2024_invivo_perturbseq
.venv/bin/igvfagent sc-analyze markers  --input <INPUT> --label zheng2024_invivo_perturbseq
```

`sc-analyze pipeline` chains: QC → normalize → PCA → UMAP → cluster →
markers. Should produce per-cluster marker genes + a UMAP figure that
shows clear neuronal cell types.

## Where artefacts land

```
Docs/SingleCell/<ts>_zheng2024_invivo_perturbseq/
├── qc_summary.json
├── pca_variance.png
├── umap.png
├── leiden_clusters.tsv
├── markers.tsv
└── summary.json
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `summary.json` exists, non-empty | Pipeline ran to completion |

(Extend after first successful run by adding range checks on
`n_cells_pre_qc`, `n_cells_post_qc`, `n_clusters`.)

## Ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| Total cells post-QC | ~25,000–30,000 (paper: ~30K) |
| Number of clusters (Leiden) | 8–15 cortical cell types: Layer-2/3, Layer-4, Layer-5, Layer-6 CT, Layer-6 IT, astrocyte, microglia, OPC, oligodendrocyte |
| Foxg1-KO cells | Should de-repress Tbr1 + Nr2f1 specifically in Layer-6 CT (Fezf2+, Ldb2+) markers |
| Nr2f1-KO cells | Should alter Layer-4 / upper-layer neuron markers |
| Tbr1-KO cells | Should alter Layer-6 IT markers |

After the pipeline runs:

```bash
# Find Layer-6 CT markers and check for Foxg1-perturbed cells in that cluster
.venv/bin/python -c "
import pandas as pd
m = pd.read_csv('Docs/SingleCell/<ts>_*/markers.tsv', sep='\t')
l6_ct = m[m['cluster'].str.contains('Layer-6|L6', case=False, na=False)]
print(l6_ct.head(20))"
```

## Running through the UI

`sc_analyze_pipeline` is registered (along with `sc_analyze_markers`).
Paste:

```
Run the Zheng 2024 in-vivo AAV Perturb-seq reproducibility benchmark.

The h5ad is at Data/Benchmarks/zheng2024_invivo_perturbseq/GSE249416.h5ad.

Step 1: sc_analyze_pipeline with input=that path, label="zheng2024_invivo_perturbseq".

Step 2: sc_analyze_markers with the same input and label.

Step 3: For each of the 4 perturbed TFs (Foxg1, Nr2f1, Tbr1, Tcf4),
report which cluster(s) show the largest expression change in the
perturbed cells (use the cluster markers + perturbation labels in
adata.obs).

Report total cells, n_clusters, top 3 perturbation effects.
```

UI sidebar: max iterations = 30 (this is a multi-step analysis),
temperature = 0.0, max tokens = 8192.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `sc-analyze pipeline` errors with "X dtype mismatch" | The h5ad's X matrix is sparse and a downstream op wants dense | `adata.X = adata.X.toarray()` before saving |
| 0 clusters returned | Resolution too low or too few cells | Pass `--leiden-resolution 0.8 1.0 1.2` to try multiple |
| Markers TSV is empty | No DE detected (likely a normalization issue) | Pre-normalize with `sc-analyze normalize --input ...` first |
| Per-perturbation effect signal is noise | Need pseudobulk DE, not per-cell | Use the `perturb-catalog dataset-rows` tool to get the precomputed pseudobulk |

## License + provenance

* **Paper data**: GEO GSE249416 — public.
* **Code**: IGVFagent Apache-2.0 (sc-analyze is a thin Scanpy wrapper, BSD-3).
* **Citation**: Zheng X et al. *Cell* **187**: 3236–3248 (2024).
  doi:10.1016/j.cell.2024.04.050 · PMID:38772369
