# Skill: Single-Cell RNA-seq, Single-Cell ATAC-seq, And Perturb-seq Analysis

Use this skill when the agent needs to find, summarize, download-plan, or analyze single-cell datasets from IGVF Portal or ENCODE.

## Metadata First

Always build a manifest before downloading data. Preserve source database, accession, assay, biosample, file format, output type, assembly, status, and href/download URL.

Run:

```bash
python3 Scripts/single_cell_data_skills.py smoke --skill all --source encode --limit 5
python3 Scripts/single_cell_data_skills.py manifest --skill scrna --source both --limit 25
python3 Scripts/single_cell_data_skills.py download-examples
python3 Scripts/single_cell_data_skills.py analyze-examples --max-cells 12000
```

## Nature/Cell-Style Analysis Pattern

Recent high-impact single-cell RNA/ATAC papers usually follow the same spine: careful cell/sample QC, latent-space construction, cell-state annotation, differential testing, and biological interpretation through regulatory elements, genes, and perturbations. The command above runs a lightweight smoke version with public IGVF examples; full publication workflows should use Scanpy/AnnData, Seurat/Signac, ArchR, or SnapATAC on full matrices.

The built-in example analysis includes a real Matrix Market expression-count pass when local 10x multiome tarballs are available. It scores marker modules and plots cells in a matrix-derived marker-score embedding; this is useful for quick interpretation before installing the heavier tSNE/UMAP stack.

For scRNA-seq:

1. Load raw counts into AnnData/Seurat and preserve raw counts.
2. QC cells by total UMIs, detected genes, mitochondrial/ribosomal fraction, doublet score, donor, batch, and chemistry.
3. Normalize, log-transform, select highly variable genes, regress only justified covariates, run PCA, neighbors, UMAP/tSNE, clustering, marker discovery, and differential expression.
4. Annotate clusters with marker genes and reference mapping; report pseudobulk differential expression for condition contrasts when donors/replicates exist.

For scATAC-seq:

1. Start from fragments, peaks, cell metadata, and genome assembly.
2. QC by fragments per cell, TSS enrichment, FRiP, nucleosome signal, blacklist fraction, and doublets.
3. Build a cell-by-peak matrix, run TF-IDF/LSI, neighbors, UMAP/tSNE, clustering, differential accessibility, motif enrichment, and gene-activity estimates.
4. Link peaks to genes with co-accessibility, multiome paired RNA, eQTL/rE2G evidence, enhancer perturbation evidence, and IGVF Catalog links.

For multiome integration:

1. Use paired cells to jointly inspect RNA expression and ATAC accessibility by cell type.
2. Integrate with weighted nearest neighbor, MultiVI/scVI-style latent models, ArchR gene scores, or Signac bridge integration.
3. Interpret variant loci by asking whether the variant overlaps accessible peaks/cCREs in the relevant cell type, whether the linked gene is expressed in the same cell type, and whether perturbation/eQTL/MPRA/rE2G data support the link.

Reference patterns used to shape this skill:

- Stuart et al., Cell 2019, Comprehensive Integration of Single-Cell Data, DOI: 10.1016/j.cell.2019.05.031.
- Hao et al., Cell 2021, Integrated analysis of multimodal single-cell data, DOI: 10.1016/j.cell.2021.04.048.
- Granja et al., Nature Genetics 2021, ArchR scalable integrative scATAC-seq analysis, DOI: 10.1038/s41588-021-00790-6.
- Ashuach et al., Nature Methods 2023, MultiVI multimodal data integration, DOI: 10.1038/s41592-023-01909-9.

## single-cell RNA-seq

Quantify gene expression across cells/nuclei and compare cell states, tissues, donors, perturbations, or conditions.

Preferred inputs:
- gene-barcode matrix
- h5ad
- matrix hdf5
- fragments are not expected for RNA

Analysis workflow:
1. Build a sample sheet from metadata: accession, biosample, donor, organism, assembly, assay, file accession, file format, and download URL.
2. Load count matrices into AnnData; preserve raw counts before normalization.
3. Run QC for library size, detected genes, mitochondrial/ribosomal fraction, doublets, and batch labels.
4. Normalize/log-transform, select highly variable genes, integrate batches when needed, cluster, annotate cell types, and run differential expression.
5. For IGVF variant/gene interpretation, summarize expression of candidate genes in relevant cell types and conditions.

## single-cell ATAC-seq

Profile chromatin accessibility at cell resolution and connect variants to accessible elements and target genes.

Preferred inputs:
- fragments
- peak calls
- cell-by-peak matrix
- bigWig signal
- h5ad or loom when available

Analysis workflow:
1. Build a manifest for fragments, peak files, cell metadata, genome assembly, and biosample ontology.
2. Run QC for fragments per cell, TSS enrichment, fraction in peaks, blacklist fraction, nucleosome signal, and doublets.
3. Create a cell-by-peak matrix, run TF-IDF/LSI, integrate batches, cluster, and annotate cell types.
4. Intersect variant lists with peaks/cCREs and summarize accessibility around candidate loci.
5. Link peaks to genes using co-accessibility, nearby genes, Catalog element-gene predictions, or matched scRNA data.

## Perturb-seq

Analyze pooled perturbation screens with single-cell expression readouts.

Preferred inputs:
- gene expression matrix
- guide assignment table
- perturbation metadata
- cell metadata
- h5ad when available

Analysis workflow:
1. Build a manifest linking expression matrices, guide assignments, perturbation targets, controls, donors, and conditions.
2. QC cells, guides, perturbation multiplicity, control cells, and batch labels.
3. Assign perturbation status per cell; remove ambiguous or high-multiplicity cells when needed.
4. Run differential expression and pathway/module scoring per perturbation against matched controls.
5. Use IGVF Catalog genes, regulatory elements, and variant links to interpret perturbation targets and downstream effects.

## Reuse Rules

- For scRNA-seq, prefer count matrices or h5ad objects and analyze with Scanpy-compatible AnnData.
- For scATAC-seq, prefer fragments plus peaks/cell metadata and analyze with ArchR, Signac, or SnapATAC-style workflows.
- For Perturb-seq, require both expression and perturbation assignment metadata before modeling effects.
- For IGVF variant interpretation, connect scATAC peaks, scRNA gene expression, Perturb-seq target effects, IGVF Catalog variant-gene evidence, and ENCODE reference context.
- Do not download large files until manifest rows are reviewed.
