# scRNA-seq Survey for the Ten IGVF Portal Files

This survey grounds the analysis in `Portal10_report.md` against (a) the foundational
Cell/Nature/Science papers that defined the assays our files came from, (b) the
IGVF consortium's Nature 2024 framing paper, and (c) the standard scRNA-seq
analysis playbook those papers use.

## What our four h5ad files actually are

We pulled the IGVF Portal record for each accession (`/matrix-files/<acc>/`) and
its parent analysis-set. The portal field `publications` is `None` for all four
files — they have not yet been linked to a peer-reviewed paper as of May 2026.
The assay context, however, is unambiguous:

| accession        | assay         | sample                    | lab / award            | h5ad layers (besides `X`)                                                                                       |
|------------------|---------------|---------------------------|------------------------|------------------------------------------------------------------------------------------------------------------|
| IGVFFI7186HBCN   | scNT-seq2     | iPSC day 16 (in-vitro)    | Hao Wu / U01HG012047   | `spliced`, `unspliced`, `ambiguous`, `total`, `labeled_TC`, `unlabeled_TC`, `sl_TC`/`sn_TC`/`ul_TC`/`un_TC`, etc. |
| IGVFFI5504UBZJ   | scNT-seq2     | iPSC day 4  (in-vitro)    | Hao Wu / U01HG012047   | same as above                                                                                                    |
| IGVFFI6422QRPZ   | scRNA-seq     | mouse                     | IGVF DACC / HG012012   | `mature`, `nascent`, `ambiguous`                                                                                 |
| IGVFFI6218UEMX   | SHARE-seq     | primary human cells       | IGVF DACC / HG012012   | `mature`, `nascent`, `ambiguous`                                                                                 |

Two are **scNT-seq2** runs from the Hao Wu lab pre-processed with the **dynast**
pipeline. Two are **DACC standard pipeline** kallisto|bustools runs (one of them
is the RNA half of a SHARE-seq joint experiment). All four ship with
`spliced/unspliced` (or `mature/nascent`) layers, which is exactly the input
shape used in RNA-velocity-style analyses in *Nature* 2018 / *Nature Methods*
2020/2023.

## The landmark papers behind each method

| Method                | Landmark paper                                                                       | What that paper actually did with the data                                                                                                                                    |
|-----------------------|---------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Standard scRNA-seq    | Wolf, Angerer & Theis. SCANPY. *Genome Biol* 2018 (& Luecken & Theis 2019)            | Defined the QC → normalize → HVG → PCA → kNN → Leiden → UMAP → marker-DE → annotation pipeline that every Cell/Nature/Science paper now follows.                              |
| kallisto\|bustools    | Melsted et al. *Nat Biotechnol* 2021                                                   | Showed that kb-python with `--workflow lamanno` outputs **spliced/unspliced/ambiguous** matrices — the format our DACC files use.                                              |
| RNA velocity          | La Manno et al. *Nature* 2018; Bergen et al. *Nat Biotechnol* 2020 (scVelo)            | Estimated past/future cell states from spliced vs unspliced ratios. Used to draw arrows on UMAPs of differentiating tissues.                                                  |
| Metabolic labeling RNA velocity | Qiu et al. *Nature Methods* 2020 (scNT-seq); Lin et al. 2025 (scNT-seq2)              | Used 4sU labeling to count **newly transcribed** (`labeled_TC`) vs **pre-existing** (`unlabeled_TC`) RNA per cell — addresses limitations of splice-only velocity.            |
| SHARE-seq             | Ma et al. *Cell* 2020                                                                  | Joint scRNA + scATAC in 34,774 mouse-skin cells; defined "domains of regulatory chromatin" (DORCs) and "chromatin potential" — chromatin opens before lineage-specific genes. |
| IGVF consortium       | The IGVF Consortium. *Nature* 2024                                                     | Frames how single-cell maps + perturbation + predictive models combine to score variants. Sets the data-product expectations that make our 10 files reusable assets.          |

## The canonical "publication-style" scRNA-seq analysis used in those papers

Every contemporary Cell/Nature/Science scRNA-seq paper (Tabula Sapiens *Science* 2022,
Human Lung Cell Atlas *Nat Med* 2023, Allen Brain Atlas *Nature* 2023, etc.)
follows essentially the same skeleton. The analyses we have **already done** are
in **plain text**; analyses that the prior versions of `Portal10_report.md`
were missing are **bold**.

1. Per-cell QC: total counts, n_genes_per_cell, % mito → filter cells.
2. Per-gene filter: cells_per_gene, drop ribosomal/mito if needed.
3. **Doublet score** (Scrublet / DoubletFinder).  *(skipped — would need scrublet)*
4. Library-size normalize → log1p (`sc.pp.normalize_total` + `log1p`).
5. **Highly-variable gene plot** (mean vs. variance, with HVGs highlighted).
6. PCA on HVGs.
7. kNN graph → Leiden / KMeans clustering.
8. UMAP and tSNE on PCs, coloured by cluster + by QC covariates.
9. **Cluster marker genes** (Wilcoxon or Welch's t-test, one-vs-rest).
10. **Marker dot plot** + **marker heatmap**: the canonical figures in Cell/Nature scRNA-seq papers.
11. Cell-type annotation by marker overlap with known panels (Tabula Sapiens,
    PanglaoDB, CellTypist) — **out of scope here** (no symbol-aware reference).
12. **Cell-cycle scoring** (Tirosh et al. *Science* 2016 S/G2M gene sets).
13. Velocity / kinetic analysis using `spliced`/`unspliced` (or
    `labeled`/`unlabeled` for scNT-seq) layers — for the IGVF files this is the
    intended downstream use.
14. Trajectory / pseudotime (PAGA, Slingshot) for differentiating systems.

What we **add in this iteration** (matches items 5, 9, 10, 12, 13 above):

- HVG mean-variance scatter with the top-2000 variance-ranked genes coloured.
- One-vs-rest **Welch's t-test** marker detection per KMeans cluster, top-5 markers.
- **Marker dot plot** (size = % cells expressing, colour = mean log-normalized).
- **Marker heatmap** on a per-cell sample (≤40 cells per cluster).
- **Cell-cycle scoring** using Tirosh S and G2/M gene panels — only run on the
  two scNT-seq2 files that carry a `gene_name` column in `var`. Plotted on UMAP.
- **Layer-aware kinetic plot**: per-cell **nascent fraction**
  `nascent / (mature + nascent)` (kb-python files) or **labeled fraction**
  `labeled_TC / total` (scNT-seq2 files), histogram + overlaid on UMAP.

## What we cannot do here without more data

- **DORC / chromatin-priming analysis** (the headline finding of the SHARE-seq
  Cell paper) needs the matched scATAC half. The portal record for
  IGVFFI6218UEMX shows `single-cell ATAC-seq (SHARE-seq), single-cell RNA
  sequencing assay (SHARE-seq)` in its parent analysis-set; the ATAC matrix is
  a sibling file we did not select.
- **RNA velocity vector field** (La Manno *Nature* 2018, scVelo *Nat Biotech*
  2020) needs scVelo, which is heavier and requires careful curation. The
  per-cell nascent-fraction plot we add here captures the same biological signal
  in a single number per cell without committing to a velocity model.
- **Variant-to-cell linkage** (the IGVF consortium's actual goal) needs
  per-donor genotype + cell-type-resolved chromatin/expression QTL pipelines
  layered on top of the cell maps; that is the work the consortium's Mapping
  Centers and Predictive Modeling Projects are building.

## How to interpret each plot in `Portal10_report.md`

- `*_counts_per_cell.png` / `*_genes_per_cell.png`: standard QC. Long left tail
  with a shoulder is a healthy filtered library. A narrow peak near zero
  indicates raw barcode output (this is what `IGVFFI6218UEMX` shows).
- `*_qc_counts_genes_mito.png`: classic "knee" + scatter. Cells with low counts
  AND high % mito are dying.
- `*_pca.png` / `*_tsne_*.png` / `*_umap_*.png`: low-dimensional structure.
  Distinct islands → real cell-state heterogeneity. One big blob with internal
  gradient → continuous variation (typical of differentiating iPSCs).
- `*_top_genes.png`: confirms what mRNA dominates the library. Mt-rRNA tops
  imply mito-biased reads; ribosomal genes dominate normal somatic cells.
- `*_hvg_mean_variance.png` (new): selected HVGs sit on the upper envelope of
  the mean-variance trend.
- `*_markers_dotplot.png` (new): each row = cluster, each column = marker. A
  diagonal pattern with diffuse off-diagonal expression is the publication
  ideal.
- `*_markers_heatmap.png` (new): cell-level expression of cluster markers.
  Block-diagonal = clean clustering; smear = continuous gradients.
- `*_cellcycle_score.png` (new, scNT-seq2 only): S vs G2/M score per cell. iPSCs
  cluster heavily in S/G2M; differentiated cells exit the cycle.
- `*_kinetic_fraction.png` (new): histogram + UMAP overlay of nascent fraction
  (or labeled fraction for scNT-seq2). High nascent/labeled fraction = actively
  transcribing cells; low = quiescent or terminal.

## Sources

- IGVF consortium framing: [Deciphering the impact of genomic variation on function](https://www.nature.com/articles/s41586-024-07510-0) — *Nature* 2024.
- IGVF catalog tooling: [IGVF catalog — from genetic variation to function](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1341/8373948) — *NAR* 2025.
- SHARE-seq: [Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin](https://www.sciencedirect.com/science/article/pii/S0092867420312538) — Ma et al. *Cell* 2020.
- scNT-seq: [Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq](https://www.nature.com/articles/s41592-020-0935-4) — Qiu et al. *Nat Methods* 2020.
- scNT-seq2: [Highly sensitive and scalable time-resolved RNA sequencing in single cells with scNT-seq2](https://pubmed.ncbi.nlm.nih.gov/40502146/) — Lin et al. 2025.
- Metabolic labeling primer: [Time-resolved single-cell RNA-seq using metabolic RNA labelling](https://www.nature.com/articles/s43586-022-00157-z) — *Nat Rev Methods Primers* 2022.
- Standard scRNA-seq pipeline: [Current best practices in single-cell RNA-seq analysis: a tutorial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6582955/) — Luecken & Theis 2019; [Comprehensive QC tools for scRNA-seq data](https://www.nature.com/articles/s41467-022-29212-9) — *Nat Commun* 2022.
- Cell-cycle gene sets: Tirosh et al. *Science* 2016.
- kallisto|bustools / kb-python: [Modular, efficient and constant-memory single-cell RNA-seq preprocessing](https://www.nature.com/articles/s41587-021-00870-2) — Melsted et al. *Nat Biotechnol* 2021.
