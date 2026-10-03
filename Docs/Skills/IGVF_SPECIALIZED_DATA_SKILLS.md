# Skill: IGVF Specialized Data Access, Download, Processing, And Analysis

Use these skills when the agent needs to find and work with specialized IGVF Portal datasets for Parse SPLiT-seq, 10x snATAC-seq with Scale pre-indexing, 10x multiome with MULTI-seq, SHARE-seq, binding effect, and SGE workflows.

## Commands

```bash
python3 Scripts/igvf_specialized_data_skills.py smoke --skill all --limit 5
python3 Scripts/igvf_specialized_data_skills.py manifest --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py download-plan --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py write-playbook
```

## Shared Rules

- Always run a metadata manifest before downloading payloads.
- Preserve accession, `@id`, assay, sample summary, file format, content type, file size, href, download URL, status, and controlled-access fields.
- Prefer released public files for automated smoke runs. Use `IGVF_PORTAL_COOKIE` only for unreleased data that the local logged-in account is authorized to access.
- Keep large payload downloads explicit and size-capped.
- Store raw API JSON under `Data/IGVF/SpecializedIGVF/`, manifests under `Data/Manifests/SpecializedIGVF/`, reports under `Docs/SpecializedIGVF/`, and logs under `Docs/Logs/`.

## 1. Parse SPLiT-seq

Skill key: `parse_split_seq`

Access and analyze Parse Biosciences SPLiT-seq single-cell RNA/nuclei datasets.

Data access:
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': 'Parse SPLiT-seq'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': 'Parse SPLiT-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'Parse SPLiT-seq'}`.
- Query IGVF Portal with `{'type': 'FileSet', 'preferred_assay_titles': 'Parse SPLiT-seq'}`.

Download inputs to prefer:
- gene count matrices
- cell annotations
- h5ad or sparse matrix archives
- barcode/sample demultiplexing metadata

Processing and analysis:
1. Build file/sample manifests from MeasurementSet and AnalysisSet records.
2. Download matrix and cell annotation payloads only after checking file size and access status.
3. Load RNA counts into AnnData, preserve raw counts, and run mitochondrial/ribosomal and library-size QC.
4. Normalize, integrate donors/batches, cluster, annotate cell types, and test gene expression across biological groups.

## 2. 10x snATAC-seq With Scale Pre-indexing

Skill key: `snatac_scale_preindexing`

Retrieve and analyze single-nucleus ATAC-seq libraries that use Scale-style pre-indexing before 10x capture.

Data access:
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': '10x snATAC-seq with Scale pre-indexing'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': '10x snATAC-seq with Scale pre-indexing'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': '10x snATAC-seq with Scale pre-indexing'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': '10x snATAC-seq with Scale pre-indexing', 'content_type': 'fragments'}`.

Download inputs to prefer:
- fragments
- peak calls or peak matrix
- cell annotations
- barcode-to-sample mapping and pre-index metadata

Processing and analysis:
1. Confirm the pre-index/sample mapping before merging fragments.
2. Compute ATAC QC metrics: fragments per nucleus, TSS enrichment, fraction in peaks, blacklist fraction, and doublets.
3. Create or reuse a cell-by-peak matrix; run TF-IDF, LSI, clustering, and cell-type annotation.
4. Intersect peaks/fragments with variants, cCREs, and candidate regulatory elements.

## 3. 10x Multiome With MULTI-seq

Skill key: `multiome_multi_seq`

Retrieve paired RNA/ATAC multiome datasets that use MULTI-seq sample multiplexing.

Data access:
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': '10x multiome with MULTI-seq'}`.
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': '10x multiome with MULTI-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': '10x multiome with MULTI-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'searchTerm': 'MULTI-seq'}`.

Download inputs to prefer:
- RNA sparse gene count matrices
- ATAC peak matrices
- ATAC fragments
- cell annotations with MULTI-seq demultiplexing labels

Processing and analysis:
1. Build paired RNA/ATAC manifests and verify matched cell barcodes and sample labels.
2. Parse MULTI-seq/demultiplexing metadata and remove ambiguous or multiplet assignments.
3. Process RNA with AnnData/Scanpy and ATAC with Signac, ArchR, or SnapATAC-style tools.
4. Integrate modalities by barcode and summarize gene expression, accessibility, cell types, and variant-overlap signals.

## 4. SHARE-seq RNA Modality

Skill key: `share_seq_rna`

Access the RNA/gene-expression side of SHARE-seq paired chromatin-expression datasets.

Data access:
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': 'SHARE-seq', 'assay_titles': 'single-cell RNA sequencing assay'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': 'SHARE-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'SHARE-seq', 'assay_titles': 'single-cell RNA sequencing assay'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'SHARE-seq', 'content_type': 'sparse gene count matrix'}`.

Download inputs to prefer:
- gene count matrix
- cell annotations
- paired ATAC cell identifiers
- sample/donor metadata

Processing and analysis:
1. Confirm SHARE-seq RNA matrices are paired to matching ATAC/chromatin files.
2. Load gene counts into AnnData and run standard RNA QC, normalization, clustering, and differential expression.
3. Keep shared cell identifiers for modality integration.
4. Use expression results to prioritize genes linked to accessible regulatory elements or variants.

## 5. Binding Effect

Skill key: `binding_effect`

Retrieve variant, motif, model, or prediction data that estimate regulatory binding effects.

Data access:
- Query IGVF Portal with `{'type': 'PredictionSet', 'file_set_type': 'binding effect'}`.
- Query IGVF Portal with `{'type': 'FileSet', 'file_set_type': 'binding effect'}`.
- Query IGVF Portal with `{'type': 'ModelSet', 'searchTerm': 'binding effect'}`.
- Query IGVF Portal with `{'type': 'File', 'searchTerm': 'binding effect'}`.

Download inputs to prefer:
- variant-level prediction tables
- model output files
- motif or TF metadata
- reference genome/assembly fields

Processing and analysis:
1. Build variant-level manifests with model, assay, biosample, assembly, and score columns.
2. Normalize score direction and effect allele conventions before joining with variant lists.
3. Annotate variants with predicted TF binding gain/loss, model provenance, and evidence confidence.
4. Integrate binding effect evidence with MPRA, CRISPRi, eQTL, accessibility, and enhancer-gene linkage evidence.

## 6. SHARE-seq ATAC/Chromatin Modality

Skill key: `share_seq_atac`

Access the chromatin-accessibility side of SHARE-seq paired expression-accessibility datasets.

Data access:
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': 'SHARE-seq', 'assay_titles': 'single-cell ATAC-seq'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': 'SHARE-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'SHARE-seq', 'assay_titles': 'single-cell ATAC-seq'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'SHARE-seq', 'content_type': 'fragments'}`.

Download inputs to prefer:
- fragments
- peak calls
- cell-by-peak matrix
- cell annotations with paired RNA identifiers

Processing and analysis:
1. Verify that ATAC barcodes map to the SHARE-seq RNA modality.
2. Run fragments/peak QC and construct cell-by-peak matrices when not already provided.
3. Run TF-IDF/LSI, clustering, motif enrichment, and peak-to-gene or co-accessibility analysis.
4. Intersect peaks with variants and connect accessibility changes to paired RNA gene expression.

## 7. SGE / Saturation Genome Editing

Skill key: `sge`

Retrieve saturation genome editing and other dense variant-function maps from IGVF Portal.

Data access:
- Query IGVF Portal with `{'type': 'MeasurementSet', 'preferred_assay_titles': 'SGE'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': 'SGE'}`.
- Query IGVF Portal with `{'type': 'File', 'preferred_assay_titles': 'SGE'}`.
- Query IGVF Portal with `{'type': 'FileSet', 'preferred_assay_titles': 'SGE'}`.

Download inputs to prefer:
- variant effect tables
- allele/count tables
- guide or edit design tables
- quality-control summaries

Processing and analysis:
1. Build a variant effect manifest with genomic coordinates, alleles, edited sequence, target gene, condition, and score.
2. Normalize variant identifiers to chr-pos-ref-alt and rsID/CAid where available.
3. Summarize score distributions, replicate concordance, controls, and significant/deleterious calls.
4. Join SGE scores to IGVF Catalog variants, genes, enhancer-gene links, MPRA, CRISPRi, and binding-effect evidence.

## 8. SGE Variant Annotation And Integration

Skill key: `sge_variant_annotation`

Analyze SGE variant-function tables in the context of IGVF Catalog/KG variant and gene evidence.

Data access:
- Query IGVF Portal with `{'type': 'File', 'searchTerm': 'variant effect'}`.
- Query IGVF Portal with `{'type': 'File', 'searchTerm': 'functional score'}`.
- Query IGVF Portal with `{'type': 'PredictionSet', 'searchTerm': 'variant effect'}`.
- Query IGVF Portal with `{'type': 'AnalysisSet', 'preferred_assay_titles': 'SGE'}`.

Download inputs to prefer:
- scored variant tables
- functional annotation tables
- gene or target annotations
- links to model/prediction provenance

Processing and analysis:
1. Parse SGE/variant-effect tables and standardize coordinates, alleles, score names, and score directions.
2. Flag coding, splice, promoter, enhancer, and UTR contexts with Catalog/KG and ENCODE annotations.
3. Rank variants by experimental score, predicted binding effect, eQTL/CRISPRi/MPRA evidence, and gene relevance.
4. Write an annotated variant table plus plots for score distribution, locus tracks, and evidence overlap.
