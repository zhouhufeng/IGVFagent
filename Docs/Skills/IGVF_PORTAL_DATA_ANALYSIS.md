# Skill: IGVF And ENCODE Data Overview And Smoke Analysis

Use this skill when the agent needs to understand what data exists in the IGVF Portal or ENCODE, choose representative input data classes, and perform shallow validation before downloading or deeply analyzing files.

## Data Map

The IGVF Data Portal is the primary ingestion point for raw files, processed files, analysis outputs, predictions, models, and metadata submissions. The IGVF Catalog consumes released portal files and filesets, plus public sources, to build Knowledge Graph tables for genes, variants, genomic elements, phenotypes, and predictions.

ENCODE is a public reference data source for regulatory genomics assays and annotations. The IGVF agent should use ENCODE to contextualize genes, variants, candidate regulatory elements, tissues/cell types, assays, tracks, peaks, and reference annotations.

Major data classes to inspect first:

- `File`: raw and processed files, file formats, file content types, size/status, and download locations.
- `MeasurementSet`: assay-centered experimental inputs and outputs.
- `AnalysisSet`: processed computational outputs and pipeline products.
- `PredictionSet`: prediction outputs such as variant effects and enhancer-gene predictions.
- `ModelSet`: trained or released computational models.
- `ConstructLibrarySet`: MPRA/construct-library linked resources.
- `Sample` and `Donor`: biological context and provenance metadata.
- `Software` and `Document`: reproducibility and protocol context.

Major ENCODE data classes to inspect first:

- `Experiment`: assay metadata for ATAC-seq, ChIP-seq, RNA-seq, reporter assays, and other functional genomics experiments.
- `File`: downloadable signal tracks, peak/region calls, alignments, quantifications, and metadata files.
- `Biosample`: tissue/cell type context, organism, treatment, donor, and ontology metadata.
- `Dataset` and `Annotation`: processed datasets and reference annotations.
- `Reference`: genome/transcriptome references and related resources.
- `Software`: pipeline/software metadata for reproducibility.

## Scripts

Run a public Catalog smoke analysis:

```bash
python3 Scripts/igvf_data_skills.py catalog-smoke --limit 10
```

Run a public ENCODE overview:

```bash
python3 Scripts/igvf_data_skills.py encode-overview --limit 25
```

Run a public ENCODE smoke analysis:

```bash
python3 Scripts/igvf_data_skills.py encode-smoke --limit 10
```

Export an ENCODE metadata table:

```bash
python3 Scripts/igvf_data_skills.py encode-export-csv --type Experiment --limit 100 --param assay_title=ATAC-seq
```

Run an IGVF Portal overview:

```bash
python3 Scripts/igvf_data_skills.py overview --limit 25
```

Run the same command with a local authenticated cookie when you need unreleased datasets:

```bash
export IGVF_PORTAL_COOKIE="<browser session cookie>"
python3 Scripts/igvf_data_skills.py overview --limit 25
```

Run a portal smoke analysis:

```bash
python3 Scripts/igvf_data_skills.py smoke --limit 25
```

Export a small metadata table:

```bash
python3 Scripts/igvf_data_skills.py export-csv --type File --limit 100
```

## Smoke Analysis Strategy

1. Check access and metadata shape with small `limit` values.
2. Summarize totals, statuses, assays, file formats, content types, labs/awards, and examples.
3. For files, inspect metadata first; avoid large downloads until the file class, size, and assay context are clear.
4. For analysis and prediction outputs, preserve links to software, pipelines, input files, and source filesets.
5. Compare portal metadata with public Catalog endpoints when the data is KG-facing.
6. Use ENCODE metadata to identify matching assay, biosample, assembly, file format, and output type before downloading large files.

## ENCODE Smoke Classes

- ATAC-seq experiments: chromatin accessibility context for regulatory elements and variants.
- ChIP-seq experiments: transcription factor and histone-mark evidence.
- RNA-seq experiments: expression context for genes and tissues.
- DNA accessibility experiments: regulatory context for noncoding variant interpretation.
- bigWig files: genome signal tracks for visualization and quantitative summaries.
- BED files: intervals, peaks, and region-level annotations.

## Outputs

- Raw JSON/text responses: `Data/`
- CSV metadata samples: `Data/`
- Reports: `Docs/IGVF_PORTAL_DATA_OVERVIEW.md`, `Docs/IGVF_PORTAL_SMOKE_ANALYSIS.md`, `Docs/IGVF_CATALOG_SMOKE_ANALYSIS.md`, `Docs/ENCODE_DATA_OVERVIEW.md`, and `Docs/ENCODE_SMOKE_ANALYSIS.md`
- Runtime logs: `Docs/Logs/`

## Authentication

Released IGVF Portal data should be accessible without login. Login as your IGVF-authorized user is needed for unreleased datasets. Keep `IGVF_PORTAL_COOKIE` local and never commit it; the script only reads this environment variable and does not write it to tracked files.
