# Skill: GEO retrieval

Wraps NCBI E-utilities + the GEO FTP mirror so the agent can search GEO, pull a Series' metadata + sample sheet, and download supplementary expression matrices from the terminal.

## Subcommands

### `search` — keyword / organism / platform search

```bash
igvfagent geo search --query 'GM12878 RNA-seq' \
    --organism 'Homo sapiens' --limit 25
igvfagent geo search --query 'lymphoblastoid H3K27ac' --study-type 'Expression profiling by high throughput sequencing'
```

### `series` — pull metadata + file inventory for a GSE

```bash
igvfagent geo series --gse GSE9574 --full-samples
```

Writes a markdown report under `Docs/GEO/<ts>_<label>_geo_report.md`, a sample sheet (`<label>_samples.csv`), and a file inventory (`<label>_files.csv`).

### `list-files` — what's available on GEO FTP for a GSE

```bash
igvfagent geo list-files --gse GSE9574
```

### `download` — pull supplementary / matrix files

```bash
igvfagent geo download --gse GSE9574 --only matrix --max-download-gb 1
igvfagent geo download --gse GSE9574 --pattern 'counts.*tsv|tpm.*tsv' --max-download-gb 2
```

By default the ``download`` step pulls everything under ``series/GSEXnnn/GSEX/{matrix,suppl,soft}/``; ``--only`` restricts to one of those categories and ``--pattern`` is a case-insensitive regex applied to the filename.

### `sample-sheet` — clean sample CSV for the rnaseq skill

```bash
igvfagent geo sample-sheet --gse GSE9574 --label gse9574_sheet
```

Reads the SOFT record's full sample block and produces a CSV with one row per GSM, columns expanded from ``characteristics_ch1`` (cell line, treatment, etc.). The output feeds straight into ``igvfagent rnaseq pipeline --sample-sheet <csv>``.

## How this chains with other skills

- After `series` + `download`, hand the sample sheet to `igvfagent rnaseq pipeline` for QC + DEG + plotting.
- The DEG list from rnaseq can then be cross-referenced with `igvfagent kg variant` / `igvfagent encode integrate-ccre` for the controlling regulatory elements.
- For literature corroboration of any GSE record, `igvfagent ref validate --input <degs.csv>`.

## Outputs

- Metadata JSON: `Data/IGVF/GEO/Metadata/<GSE>.json`
- Downloads:     `Data/IGVF/GEO/Downloads/<GSE>/`
- Manifests:     `Data/Manifests/GEO/<ts>_<label>_*.csv`
- Reports:       `Docs/GEO/<ts>_<label>_geo_report.md`