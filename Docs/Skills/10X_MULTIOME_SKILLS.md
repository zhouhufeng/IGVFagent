# Skill: IGVF 10x Multiome Retrieval And Analysis

Use this skill when the agent needs to retrieve, summarize, or prepare paired 10x multiome RNA/ATAC datasets from the IGVF Portal.

## Run

```bash
python3 Scripts/multiome_10x_pipeline.py retrieve --count 20 --fetch-file-details
python3 Scripts/multiome_10x_pipeline.py process-local --file-manifest Data/Manifests/Multiome10x/<files.csv> --download-manifest Data/Manifests/Multiome10x/<download_manifest.csv>
```

The script searches public released IGVF Portal `AnalysisSet` records where `preferred_assay_titles=10x multiome`, preserves full JSON metadata, writes dataset/file/sample manifests, and creates a processing report.

## Inputs To Prefer

- RNA sparse gene count matrices.
- ATAC annotated sparse peak count matrices.
- ATAC fragments.
- Cell annotations for barcode-level modality integration.

## Processing Pattern

1. Build file and sample manifests before downloading payloads.
2. Verify each file is released, public, and aligned to the expected assembly and transcriptome annotation.
3. Load RNA matrices into AnnData/Scanpy and ATAC matrices or fragments into ArchR, Signac, or SnapATAC-style workflows.
4. Join modalities with cell annotations and barcodes.
5. Summarize QC metrics, cell states, accessible peaks, gene expression, and peak-gene or variant-gene context for IGVF interpretation.

## Download Policy

Use `--download-policy processed --max-download-gb N` only after reviewing the file manifest. Raw and processed multiome payloads can be large.
