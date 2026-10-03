# Skill: Single-cell analysis (UMAP / t-SNE / Leiden / markers)

End-to-end Scanpy-driven single-cell workflow. Closes the gap where
IGVFagent could discover and download counts matrices but had to hand
off to "run Scanpy or Seurat" for the actual analysis.

## Subcommands

### qc
```
igvfagent sc-analyze qc --input counts.h5ad --label demo \
    --min-genes 200 --min-cells 3 --max-mito 20
```
Loads any of `.h5ad` / 10x `.h5` / `.mtx` / `.csv` / `.tsv`, filters
cells & genes, computes mito %, drops high-mito cells, writes QC
violins, and saves the cleaned `processed.h5ad`.

### normalize
```
igvfagent sc-analyze normalize --input processed.h5ad --n-hvg 2000
```
Total-count normalize, log1p, select top-N highly-variable genes,
scale to unit variance. Saves raw counts in `.raw` for downstream
marker DE.

### pca / umap / tsne / cluster / markers
```
igvfagent sc-analyze pca       --input processed.h5ad --n-pcs 50
igvfagent sc-analyze umap      --input processed.h5ad --n-neighbors 15
igvfagent sc-analyze tsne      --input processed.h5ad --n-pcs 30
igvfagent sc-analyze cluster   --input processed.h5ad --resolution 1.0
igvfagent sc-analyze markers   --input processed.h5ad --n-top 25
```
Each step accepts the `processed.h5ad` from any previous step and
saves a new `processed.h5ad` carrying the additional fields
(`X_pca`, `X_umap`, `X_tsne`, `leiden`, `rank_genes_groups`).

### pipeline (recommended one-shot)
```
igvfagent sc-analyze pipeline --input counts.h5ad --label k562_demo \
    --min-genes 200 --min-cells 3 --max-mito 20 \
    --n-hvg 2000 --n-pcs 50 --resolution 1.0 \
    --sample-col sample \
    --highlight-genes APOE,TREM2,LDLR,GFAP
```
Drives QC → normalize → PCA → neighbors → UMAP → t-SNE → Leiden →
marker DE → figures → markdown report in one go. Saves the
`processed.h5ad` so any later step can resume.

### plot-embedding
```
igvfagent sc-analyze plot-embedding --input processed.h5ad \
    --embedding umap --color leiden,APOE,TREM2 --label apoe_view
```
Re-render UMAP or t-SNE coloured by any combination of obs columns
and genes. Useful for asking the agent "show me APOE expression
overlaid on the UMAP we just made."

## Inputs

| Format | Detection rule |
|---|---|
| AnnData `.h5ad` | suffix `.h5ad` |
| 10x CellRanger | suffix `.h5` |
| 10x sparse | suffix `.mtx` / `.mtx.gz` (looks for `barcodes.tsv` + `features.tsv` in same dir) |
| Text matrix | `.csv` / `.tsv` / `.txt` (and `.gz` variants) — auto-orients to cells×genes |

## Outputs

```
Docs/SingleCell/<ts>_<label>/
    processed.h5ad
    qc_summary.json
    markers.csv
    Plots/
        qc_violins.png
        pca_variance.png
        umap_clusters.png
        umap_sample.png            (if --sample-col given)
        umap_auto_top_markers.png  (auto-picks 4 cluster markers)
        tsne_clusters.png
        heatmap_top_markers.png
    report.md
```

## Cross-skill chaining

- `singlecell` / `multiome` / `splitseq` — discover datasets, then
  feed their downloaded counts straight into `sc-analyze pipeline`.
- `geo download` → `sc-analyze pipeline` — wide net retrieval from
  GEO, then end-to-end analysis.
- `kg` / `proteomics kg-visualize` — after `sc-analyze markers`,
  drop the top marker genes into the KG visualizer to see their
  protein-protein interaction context.

## Dependencies

`scanpy>=1.10`, `anndata>=0.10`, `umap-learn`, `scikit-learn`,
`leidenalg`, `python-igraph`, `matplotlib`. Installed in the
project venv. If you're using an external Python, install with:

```
pip install scanpy 'anndata>=0.10' umap-learn leidenalg python-igraph
```
