# OPERATIONS — Travaglini 2020 Human Lung Cell Atlas (`sc-analyze`)

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run (~8–10 min; first run downloads ~596 MB)

```bash
bash Benchmarks/travaglini2020_lung/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark travaglini2020_lung
```

Expected: `8/8 checks PASSED`.

## Steps `run.sh` performs

1. **Download** the CELLxGENE h5ad (asset `f5568ea3-…h5ad`, collection DOI
   `10.1038/s41586-020-2922-4`) into `Benchmarks/_data/travaglini2020_lung/`.
   Skipped if already present.
2. **`prep_input.py`** — standard CELLxGENE→Scanpy conversion: raw counts
   (`.raw.X`) → `.X`, Ensembl IDs → `feature_name` gene symbols, retain
   `cell_type` / `free_annotation` / `compartment` in `obs`. No filtering.
3. **`sc-analyze pipeline`** — QC → normalize → 2,000 HVG → 50-PC PCA →
   neighbors → UMAP → Leiden (res 1.5) → Wilcoxon markers. Author labels ride
   along via `--sample-col cell_type` for the UMAP overlay.
4. **`make_figures.py`** — scores the Leiden partition against the 46 author
   `cell_type` labels (ARI / AMI / homogeneity / completeness / v-measure),
   checks canonical-marker → cell-type biology, writes `concordance_metrics.json`
   (into both the run dir and `figures/`) plus `fig1_confusion` + `fig2_markers`.

## Where artefacts land

```
Docs/SingleCell/<ts>_travaglini2020_lung/
  qc_summary.json            pipeline summary (n_clusters, n cells, markers_total)
  concordance_metrics.json   ARI/AMI/homogeneity + marker detail  ← scored artefact
  markers.csv                top-25 Wilcoxon markers per Leiden cluster
  processed.h5ad             clustered AnnData (leiden + author cell_type in obs)
  Plots/umap_clusters.png  Plots/umap_sample.png  Plots/heatmap_top_markers.png
```

## Concordance interpretation

| # | Check | Verifies |
|---|---|---|
| 1 | `n_cells` ∈ [65000, 66000] | full atlas loaded |
| 2 | `n_leiden_clusters` ∈ [35, 65] | cluster count ≈ 46 author types |
| 3 | `AMI` ≥ 0.70 | unsupervised partition agrees with author annotation |
| 4 | `homogeneity` ≥ 0.80 | clusters are pure (one author type each) |
| 5 | canonical markers in target top-3 ∈ [8, 10] | marker biology reproduced |
| 6–8 | `processed.h5ad` / `markers.csv` / `qc_summary.json` present | pipeline ran to completion |

## Memory note

`make_figures.py` reads `processed.h5ad` in **backed mode** and pulls only the
marker columns from the prep file, so it runs in < 1 GB RAM. The `sc-analyze`
pipeline itself peaks at ~6–8 GB on 65 k cells; close other large processes if
RAM-constrained.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `prep_input.py` can't find `lung_atlas.h5ad` | download step skipped/failed | re-run `run.sh`, or fetch the asset URL manually |
| pipeline OOM-killed (exit 137) | < 8 GB free | free RAM; or pre-subsample the h5ad |
| `n_clusters` far from 46 | resolution changed | keep `--resolution 1.5` for the scored run |

## License + provenance

* **Data**: CELLxGENE Discover, CC-BY 4.0 — fetch/link only, never redistributed.
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Travaglini KJ et al. *Nature* **587**: 619–625 (2020).
  doi:10.1038/s41586-020-2922-4 · PMID 33208946
