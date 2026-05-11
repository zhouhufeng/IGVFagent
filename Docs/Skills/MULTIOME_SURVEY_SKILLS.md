# Skill: Multiome cross-source survey

Surveys **IGVF Portal**, **ENCODE**, and **GEO** for single-cell multiome
datasets (10x Multiome / SHARE-seq / single-nucleus multiome), produces a
unified file manifest, and downloads the training-relevant files into
``Data/MultiomeSurvey/``.

## Subcommands

### `survey-igvf`

```bash
igvfagent multiome-survey survey-igvf --limit 200 --fetch-files
```

Queries ``preferred_assay_titles=10x multiome`` (and SHARE-seq) on the
IGVF Portal.  Writes:

- `Data/Manifests/MultiomeSurvey/<ts>_<label>_igvf_datasets.csv`
- `Data/Manifests/MultiomeSurvey/<ts>_<label>_igvf_files.csv`
- `Docs/MultiomeSurvey/<ts>_<label>_igvf_report.md`

### `survey-encode`

```bash
igvfagent multiome-survey survey-encode --limit 200 --fetch-files
```

Searches ENCODE Experiments and Series for the multiome / SHARE-seq keywords.

### `survey-geo`

```bash
igvfagent multiome-survey survey-geo --limit 50 --organism 'Homo sapiens' --fetch-files
```

Runs NCBI E-utilities `esearch+esummary` on the `gds` database for the
canonical multiome queries.  Optionally scrapes the GEO FTP listing for
each GSE to expose supplementary files.

### `survey-all`

```bash
igvfagent multiome-survey survey-all --limit 100 --fetch-files
```

Runs all three.

### `manifest`

```bash
igvfagent multiome-survey manifest --label unified_v1
```

Builds a unified file manifest from the **most recent** per-source files
CSVs.  Use this as the input to `download`.

### `download`

```bash
igvfagent multiome-survey download \
    --manifest Data/Manifests/MultiomeSurvey/<ts>_<label>_unified_manifest.csv \
    --only matrix_rna,matrix_atac,fragments,annotations \
    --max-download-gb 20
```

Cap defaults to 5 GB.  ``--only`` accepts a comma-separated list of
kinds: ``matrix_rna``, ``matrix_atac``, ``fragments``, ``peaks``,
``annotations``, ``alignments``, ``index``, ``raw_reads``, ``other``.
``--pattern`` is a case-insensitive regex applied to the filename + URL.

### `inventory`

```bash
igvfagent multiome-survey inventory
```

Walks ``Data/MultiomeSurvey/`` and writes a fresh inventory CSV.

## Output layout

```
Data/
  MultiomeSurvey/
    igvf/<IGVFDS...>/<file>
    encode/<ENCSR...>/<file>
    geo/<GSE...>/<file>
    inventory.csv
  Manifests/MultiomeSurvey/
    <ts>_<label>_<source>_datasets.csv
    <ts>_<label>_<source>_files.csv
    <ts>_<label>_unified_manifest.csv
    <ts>_<label>_unified_manifest_download_log.csv
Docs/
  MultiomeSurvey/
    <ts>_<label>_<source>_report.md
  Skills/MULTIOME_SURVEY_SKILLS.md   (this file)
```

## Privacy

All endpoint URLs are resolved through `Scripts/_endpoints.py`; no URLs,
cookies, or credentials are written to source.  `Data/MultiomeSurvey/`
matches `Data/*` in the repo `.gitignore`, so downloaded payload never
accidentally lands in commits.
