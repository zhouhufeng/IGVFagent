# Multiome survey — run results and download recipes

Source: `Data/Manifests/MultiomeSurvey/20260511_120333_results_v1_unified_manifest.csv`  
Generated: 1778515413.45

## Totals

- **401 multiome datasets** across the six surveyed sources
- **3,837 downloadable files indexed**
- **14.19 TB** total file size if every indexed file were downloaded

## Per-source breakdown

| source | datasets | files | total size | top kinds (count) |
|---|---:|---:|---:|---|
| igvf | 90 | 507 | 1.65 TB | index=131, matrix_rna=99, other=94, fragments=90 |
| encode | 31 | 361 | 2.04 TB | raw_reads=179, matrix_rna=66, alignments=53, index=30 |
| geo | 104 | 774 | 0.00 B | other=659, index=104, fragments=10, annotations=1 |
| cellxgene | 109 | 131 | 619.79 GB | matrix_rna=109, fragments=22 |
| hca | 12 | 1,840 | 9.48 TB | raw_reads=1370, other=362, fragments=101, matrix_atac=7 |
| zenodo | 55 | 224 | 421.55 GB | other=204, fragments=16, annotations=3, raw_reads=1 |

## Size by kind, across all sources

| kind | files | total size |
|---|---:|---:|
| raw_reads | 1,550 | 10.24 TB |
| alignments | 130 | 2.22 TB |
| other | 1,337 | 651.13 GB |
| fragments | 254 | 593.07 GB |
| matrix_rna | 274 | 449.11 GB |
| index | 265 | 70.44 GB |
| matrix_atac | 15 | 1.17 GB |
| annotations | 12 | 52.10 MB |

## Top-10 largest single files

| source | dataset | file | kind | size |
|---|---|---|---|---:|
| igvf | IGVFDS7548PGXM | IGVFFI3288ZMNY | alignments | 54.54 GB |
| cellxgene | 6f7fd0f1-a2ed-4ff1 | dad4819b-4c14-439c-b32a-2c8d68bd22e1 | matrix_rna | 49.54 GB |
| hca | 6601b3d4-ed5a-4e1f | FCA_GND8715519_S1_L001_R2_001.fastq. | raw_reads | 47.81 GB |
| cellxgene | c2876b1b-06d8-4d96 | 92b37feb-aa2c-40d7-bd90-0a9b5ddb3b27 | matrix_rna | 47.29 GB |
| cellxgene | 19f1bec6-eeb7-4133 | 4ec25136-f4bb-4fed-9ed5-133911c1d572 | fragments | 47.04 GB |
| zenodo | 17063631 | pHGGmap_discovery_cohort_ATAC_fragme | fragments | 46.65 GB |
| hca | 6601b3d4-ed5a-4e1f | FCA_GND10287604_S1_L001_R2_001.fastq | raw_reads | 46.29 GB |
| hca | 6601b3d4-ed5a-4e1f | FCA_GND8715519_S1_L001_R1_001.fastq. | raw_reads | 45.85 GB |
| igvf | IGVFDS3722RMQB | IGVFFI3202IZUS | alignments | 44.58 GB |
| encode | ENCSR083HFN | ENCFF863CZA | alignments | 44.14 GB |

## How to download

All commands resolve their endpoints through ``Scripts/_endpoints.py``; no URLs are hard-coded.  Data lands under ``Data/MultiomeSurvey/<source>/<accession>/`` which is gitignored.

### Quickstart — RNA + ATAC + annotations only, 20 GB cap

```bash
igvfagent multiome-survey download \
    --only matrix_rna,matrix_atac,fragments,annotations \
    --max-download-gb 20
```

### Just RNA H5AD (CELLxGENE + IGVF + GEO matrices), 10 GB cap

```bash
igvfagent multiome-survey download \
    --only matrix_rna --max-download-gb 10
```

### Per-source filtered downloads with a regex

```bash
# Just IGVF h5ad RNA + ATAC fragments
igvfagent multiome-survey download --pattern 'IGVFFI.*\.(h5ad|bed\.gz)$' --max-download-gb 5

# CELLxGENE H5AD only
igvfagent multiome-survey download --pattern '\.h5ad$' --max-download-gb 5
```

### Dry-run first to see what you'd pull

```bash
igvfagent multiome-survey download \
    --only matrix_rna,annotations --max-download-gb 5 --dry-run
```

## Per-source download notes

### IGVF

- Datasets indexed: **90**, files indexed: **507**, total: **1.65 TB**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_115543_results_v1_igvf_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_115543_results_v1_igvf_files.csv`
- Example largest files:
    - IGVFFI3288ZMNY (alignments) — 54.54 GB  → `https://api.data.igvf.org/alignment-files/IGVFFI3288ZMNY/@@download/IGVFFI3288ZMNY.bam`
    - IGVFFI3202IZUS (alignments) — 44.58 GB  → `https://api.data.igvf.org/alignment-files/IGVFFI3202IZUS/@@download/IGVFFI3202IZUS.bam`
    - IGVFFI7694PTWA (alignments) — 40.09 GB  → `https://api.data.igvf.org/alignment-files/IGVFFI7694PTWA/@@download/IGVFFI7694PTWA.bam`

### ENCODE

- Datasets indexed: **31**, files indexed: **361**, total: **2.04 TB**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_115826_results_v1_encode_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_115826_results_v1_encode_files.csv`
- Example largest files:
    - ENCFF863CZA (alignments) — 44.14 GB  → `https://www.encodeproject.org/files/ENCFF863CZA/@@download/ENCFF863CZA.bam`
    - ENCFF060RAT (alignments) — 39.38 GB  → `https://www.encodeproject.org/files/ENCFF060RAT/@@download/ENCFF060RAT.bam`
    - ENCFF905UOG (alignments) — 36.06 GB  → `https://www.encodeproject.org/files/ENCFF905UOG/@@download/ENCFF905UOG.bam`

### GEO

- Datasets indexed: **104**, files indexed: **774**, total: **0.00 B**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_115906_results_v1_geo_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_115906_results_v1_geo_files.csv`

### CELLXGENE

- Datasets indexed: **109**, files indexed: **131**, total: **619.79 GB**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_120104_results_v1_cellxgene_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_120104_results_v1_cellxgene_files.csv`
- Example largest files:
    - dad4819b-4c14-439c-b32a-2c8d68bd22e1.h5ad (matrix_rna) — 49.54 GB  → `https://datasets.cellxgene.cziscience.com/dad4819b-4c14-439c-b32a-2c8d68bd22e1.h5ad`
    - 92b37feb-aa2c-40d7-bd90-0a9b5ddb3b27.h5ad (matrix_rna) — 47.29 GB  → `https://datasets.cellxgene.cziscience.com/92b37feb-aa2c-40d7-bd90-0a9b5ddb3b27.h5ad`
    - 4ec25136-f4bb-4fed-9ed5-133911c1d572-fragment.ts (fragments) — 47.04 GB  → `https://datasets.cellxgene.cziscience.com/4ec25136-f4bb-4fed-9ed5-133911c1d572-fragment.ts`

### HCA

- Datasets indexed: **12**, files indexed: **1,840**, total: **9.48 TB**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_120134_results_v1_hca_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_120134_results_v1_hca_files.csv`
- Example largest files:
    - FCA_GND8715519_S1_L001_R2_001.fastq.gz (raw_reads) — 47.81 GB  → ``
    - FCA_GND10287604_S1_L001_R2_001.fastq.gz (raw_reads) — 46.29 GB  → ``
    - FCA_GND8715519_S1_L001_R1_001.fastq.gz (raw_reads) — 45.85 GB  → ``

### ZENODO

- Datasets indexed: **55**, files indexed: **224**, total: **421.55 GB**.
- Datasets manifest: `Data/Manifests/MultiomeSurvey/20260511_120333_results_v1_zenodo_datasets.csv`
- Files manifest:    `Data/Manifests/MultiomeSurvey/20260511_120333_results_v1_zenodo_files.csv`
- Example largest files:
    - pHGGmap_discovery_cohort_ATAC_fragments.tsv.gz (fragments) — 46.65 GB  → `https://zenodo.org/api/records/17063631/files/pHGGmap_discovery_cohort_ATAC_fragments.tsv.`
    - embryo_data_zenodo.tar.gz (other) — 29.62 GB  → `https://zenodo.org/api/records/15281826/files/embryo_data_zenodo.tar.gz/content`
    - hg38.tar.gz (other) — 28.81 GB  → `https://zenodo.org/api/records/20045739/files/hg38.tar.gz/content`
