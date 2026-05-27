# Mitra 2024 — SCARlink single-cell multiome regression

## Citation

Mitra S, ..., Leslie CS. *Nat Genet* **56**: 627–636 (2024).
DOI: [10.1038/s41588-024-01689-8](https://doi.org/10.1038/s41588-024-01689-8) · PMC: PMC11018525

## Data

* **GEO**: GSE194122 (10x Multiome PBMC + bone marrow — reanalysis)
* **GitHub**: [snehamitra/SCARlink](https://github.com/snehamitra/SCARlink)

## Headline workflow (paper)

1. 10x Multiome (paired snRNA + snATAC) on PBMC + bone marrow.
2. SCARlink Poisson regression to link each peak to its target gene's expression.
3. Recover chromatin-potential trajectories that distinguish erythroid (GATA1) and monocyte (IRF8) lineages.

## What IGVFagent reproduces

* **Online step** (this benchmark's scope): `multiome retrieve` enumerates IGVF-portal 10x Multiome AnalysisSets and writes a manifest.
* **Optional local step**: if you've downloaded a multiome `.mtx` bundle (or use one of the IGVF AnalysisSets), `multiome process-local` → `joint-qc` → `lsi` → `wnn` → `peak2gene` runs the full Signac-style pipeline. The `peak2gene` output is the direct analogue of SCARlink's per-peak gene-link table.

## Ground truth to spot-check

* GATA1 enhancer should appear as a top correlate in erythroid-lineage cells
* IRF8 enhancer should appear as a top correlate in monocyte-lineage cells
* Manifest should enumerate at least 3 IGVF multiome AnalysisSets
