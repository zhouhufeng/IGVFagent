# 10x Multiome on the IGVF Portal: Inventory, Methods, and Cell/Nature/Science Playbook

This survey complements `Docs/SingleCell/Multiome10_report.md` (per-set
analyses) by (a) summarizing what 10x Multiome data the IGVF Portal hosts as
of May 2026, (b) listing the Cell/Nature/Science papers whose analyses we
mimic, and (c) explaining why we picked the file subset we did and what we
intentionally did not run.

## Inventory: 10x Multiome on IGVF Portal

A live portal query (`GET /search/?type=AnalysisSet&preferred_assay_titles=10x multiome`)
returned the following totals:

- **AnalysisSets** with `preferred_assay_titles=10x multiome`: **1,954**
- **MeasurementSets** with `preferred_assay_titles=10x multiome`: **1,604**
- **Files** with `content_type=annotated sparse peak count matrix`: **506** total
  - 505 stored as `.rds` (R-only)
  - 1 as `h5ad`
- **Files** with `content_type=sparse gene count matrix`: **2,088** total
  - 1,417 as `h5ad`
  - 671 as `.tar` (Cell Ranger ARC triplet: matrix.mtx + features.tsv + barcodes.tsv)
- **Files** with `content_type=cell annotations`: 505 (TSV)
- **Files** with `content_type=fragments`: standard 10x ATAC fragment BED.gz
  per AnalysisSet — typically **~250 MB to >2 GB each**

Each canonical IGVF DACC multiome AnalysisSet ships **four files**:

1. `annotated sparse peak count matrix` (.rds) — joint **cell × peak** ATAC matrix with embedded cell-type annotations (R/Seurat/Signac native).
2. `sparse gene count matrix` (.tar) — **cell × gene** RNA matrix (Cell Ranger ARC output).
3. `fragments` (.bed.gz) — per-fragment ATAC reads with cell barcode.
4. `cell annotations` (.tsv.gz) — barcode-keyed cell-type and QC metadata table.

The dominant matrix format on the portal is the IGVF DACC standard pipeline
(`tar` for RNA + `rds` for ATAC). h5ad-formatted RNA exists and is from a
parallel kb-python pipeline; h5ad-formatted ATAC is essentially absent (n=1
across the entire portal).

## What we downloaded and analyzed

Out of the 30-set portal sample we manifested via
`Scripts/multiome_10x_pipeline.py retrieve --count 30`, we picked the **10
smallest AnalysisSets** by combined RNA-tar + cell-annotations size for
download. For the **3 smallest** we also downloaded the `fragments` BED.gz so
that a real ATAC QC stage is included in the per-set analysis.

We deliberately **skipped** the `.rds` ATAC peak matrices because they require
R (`Matrix::dgCMatrix`); pure-Python readers (`pyreadr` etc.) do not handle the
sparse-matrix dispatch reliably. Joint cell-by-peak analysis with chromatin
embedding (LSI on TF-IDF), peak motif enrichment, and peak-gene linkage is the
piece this survey explicitly does not perform — it is doable but would either
require an R sub-process or the rare h5ad ATAC files.

## The Cell/Nature/Science 10x Multiome analysis playbook

The conventions used in current top-journal multiome papers are remarkably
consistent. The list below maps each pipeline stage to what our script does
(or chooses not to do) on IGVF data.

| Stage                                  | Reference paper                                                              | Our analysis                                                            |
|----------------------------------------|-------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| Cell Ranger ARC pre-processing          | 10x Genomics technical brief                                                  | Inputs to IGVF DACC pipeline; we read its outputs (tar + fragments BED). |
| Per-cell QC: nUMI / nGene / %mito (RNA) | Luecken & Theis *Mol Sys Biol* 2019                                           | Counts/cell, genes/cell, % mito plot per set.                            |
| Per-cell QC: nFrag / TSS / FRiP (ATAC)  | Granja et al. *Nat Genet* 2021 (ArchR); Stuart et al. *Nat Methods* 2021 (Signac) | Fragments BED → per-cell fragment count + barcode-rank knee plot.        |
| Doublet / nucleus filtering             | Wolock et al. (Scrublet)                                                      | Skipped (no scrublet dep); we filter by n_genes ≥ 200.                   |
| Fragment-size distribution              | Buenrostro et al. *Nat Methods* 2013 (ATAC-seq); 10x technical note           | Histogram with mononucleosome (~147 bp) and dinucleosome bands marked.   |
| Library-size normalize + log1p (RNA)    | Stuart et al. *Cell* 2019                                                     | Standard normalize_total(1e4) + log1p.                                   |
| HVG selection                           | scanpy / Seurat                                                                | Top 2000 genes by variance on log-norm.                                  |
| PCA on RNA HVGs                         | Wolf et al. *Genome Biol* 2018 (scanpy); Hao et al. *Cell* 2021 (Seurat WNN)  | Truncated SVD → top 50 PCs.                                              |
| TF-IDF + LSI on ATAC peaks              | Cusanovich et al. *Cell* 2018; Stuart et al. (Signac)                         | Skipped — `.rds` peak matrices not readable in pure Python.              |
| Joint kNN graph (Weighted Nearest Neighbors) | Hao et al. *Cell* 2021 (Seurat v4)                                       | Skipped (no peak matrix); we substitute matched-barcode joint scatter.   |
| Cell-type annotation (marker overlap)   | Tabula Sapiens *Science* 2022; Yao et al. *Nature* 2023                       | We surface **IGVF-provided `cell annotations`** column on the UMAP.      |
| Cluster marker DE (RNA)                 | scanpy `rank_genes_groups`                                                    | Welch's t-test, top 5 markers / cluster.                                 |
| Marker dot plot                         | Standard Seurat/scanpy figure                                                  | Per-set dot plot.                                                        |
| Cell-cycle scoring                      | Tirosh et al. *Science* 2016                                                   | S/G2M scores when symbols available, UMAP overlay, phase counts.         |
| Peak-gene linkage                       | Ma et al. *Cell* 2020 (SHARE-seq DORCs); Mitchel et al. *Nat Genet* 2024 (SCARlink) | Skipped — needs the peak matrix.                                         |
| Motif enrichment (chromVAR)             | Schep et al. *Nat Methods* 2017                                                | Skipped — needs peak matrix + JASPAR/HOMER motif scan.                   |
| Trajectory / pseudotime                 | Wolf et al. PAGA *Genome Biol* 2019; Bergen et al. scVelo *Nat Biotechnol* 2020 | Out of scope; UMAP is qualitative substitute for a small ten-set demo. |

## Per-set figures produced

For each AnalysisSet the script emits, into
`Docs/SingleCell/Plots/Multiome10/<accession>/`:

- `rna_qc_counts_genes_mito.png` — counts/cell vs n_genes coloured by % mito.
- `rna_umap_clusters.png` — RNA UMAP coloured by KMeans cluster.
- `rna_umap_counts.png` — RNA UMAP coloured by log10 counts.
- `rna_umap_celltype.png` — RNA UMAP coloured by IGVF cell-type label, where
  the cell-annotations TSV provides one (top 12 labels).
- `rna_markers_dotplot.png` — top-5 marker genes per cluster, dot plot.
- `rna_cellcycle.png` — Tirosh S vs G2M score scatter + UMAP overlay (only
  when `gene_name` symbols are available in the matrix).
- `atac_fragment_size.png` — fragment-size histogram with mono- and
  dinucleosome lines (only for the 3 sets that include the fragments BED).
- `atac_barcode_rank.png` — barcode-rank fragments-per-barcode plot.
- `joint_rna_vs_atac_depth.png` — per-cell scatter of log10(RNA UMI+1) vs
  log10(ATAC frag+1) on barcodes shared by both modalities.
- `joint_umap_atac_overlay.png` — RNA UMAP coloured by per-cell ATAC depth.

## Landmark Cell/Nature/Science papers using 10x Multiome

| Paper                                                                                   | Year / Venue          | What they did with multiome                                                                     |
|------------------------------------------------------------------------------------------|------------------------|--------------------------------------------------------------------------------------------------|
| Hao et al. — Integrated analysis of multimodal single-cell data                          | *Cell* 2021            | Defined Seurat **WNN** integration on multiome; CITE-seq + RNA + ATAC PBMCs.                      |
| Ma et al. — Chromatin Potential Identified by Shared Single-Cell Profiling                | *Cell* 2020            | SHARE-seq (multiome predecessor); mouse skin DORCs; chromatin priming.                           |
| Stuart et al. — Single-cell chromatin state analysis with **Signac**                      | *Nat Methods* 2021     | Joint scATAC + scRNA pipelines, per-cell weighting, peak gene linkage.                           |
| Granja et al. — **ArchR** for snATAC-seq                                                  | *Nat Genet* 2021       | The de-facto ATAC half of multiome workflows (LSI, gene activity, peak2gene).                    |
| Yao et al. — High-resolution transcriptomic and spatial atlas of whole mouse brain        | *Nature* 2023          | 10x Multiome supplemented MERFISH atlas; 1,687 nuclei across 33 clusters; 2.3M cells snATAC.     |
| Kanemaru et al. — Spatially resolved multiomics of human cardiac niches                   | *Nature* 2023          | snRNA + snATAC (10x Multiome) + Visium across 8 cardiac regions; cellular niches.                 |
| Mitchel et al. — **SCARlink**: gene regulatory model from multi-ome                       | *Nat Genet* 2024       | Per-gene regression linking enhancers to expression at single-cell resolution from 10x Multiome. |
| The IGVF Consortium — Deciphering the impact of genomic variation on function              | *Nature* 2024          | Frames the 10x Multiome data products on the portal as variant-to-function mapping inputs.        |

## Sources

- IGVF consortium framing: [Deciphering the impact of genomic variation on function](https://www.nature.com/articles/s41586-024-07510-0) — *Nature* 2024.
- 10x Genomics product page: [Epi Multiome](https://www.10xgenomics.com/products/epi-multiome).
- WNN / Seurat v4: [Integrated analysis of multimodal single-cell data](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) — Hao et al. *Cell* 2021.
- Signac: [Joint RNA and ATAC analysis with Signac (10x multiomic vignette)](https://stuartlab.org/signac/articles/pbmc_multiomic).
- ArchR: Granja et al. *Nat Genet* 2021.
- SHARE-seq: [Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin](https://www.sciencedirect.com/science/article/pii/S0092867420312538) — Ma et al. *Cell* 2020.
- Whole mouse brain atlas (multiome supplement): [A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain](https://www.nature.com/articles/s41586-023-06812-z) — Yao et al. *Nature* 2023.
- Heart niches: [Spatially resolved multiomics of human cardiac niches](https://www.nature.com/articles/s41586-023-06311-1) — Kanemaru et al. *Nature* 2023.
- Brain vasculature atlas: [Single-cell atlas of the human brain vasculature](https://www.nature.com/articles/s41586-024-07493-y) — *Nature* 2024.
- Multi-ome regression model: [Single-cell multi-ome regression models identify functional and disease-associated enhancers](https://www.nature.com/articles/s41588-024-01689-8) — Mitchel et al. *Nat Genet* 2024 (SCARlink).
- Fragment-size baseline: Buenrostro et al. ATAC-seq *Nat Methods* 2013.
- Cell-cycle markers: Tirosh et al. *Science* 2016.
