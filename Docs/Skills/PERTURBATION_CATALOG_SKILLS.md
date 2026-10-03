# Skill: Perturbation Catalogue retrieval

Pulls metadata and per-row perturbation effects from the Perturbation
Catalogue Search API. The catalogue indexes **~1,222 perturbation
datasets** across MAVE (DMS / VAMP-seq family), CRISPR screens, and
Perturb-seq.

## Subcommands

### summary
```
igvfagent perturb-catalog summary
```
Landing-page stats: total datasets, top modalities / tissues / cell
types / diseases / perturbation types. Writes
`Docs/Perturbation/<ts>_summary/{summary.json, report.md}`.

### search (global, gene-level)
```
igvfagent perturb-catalog search --query BRCA1 --size 10
```
Returns one row per perturbed gene with counts in each modality, top
GSEA terms, and faceted tissue/cell-type info. Use this when the
user asks "what does the catalogue have on gene X?".

### search-modality
```
igvfagent perturb-catalog search-modality --modality crispr-screen \
    --perturbation-gene-name BRCA1 --dataset-limit 25
igvfagent perturb-catalog search-modality --modality mave \
    --perturbation-gene-name TP53 \
    --perturbation-position 100_300 \
    --effect-score-name vamp_score
igvfagent perturb-catalog search-modality --modality perturb-seq \
    --query "lung cancer" --dataset-limit 10
```
Modality-scoped search returns datasets + sample rows per dataset.
Supports filters on perturbation gene name, position range (e.g.
`100_300`), score name and value range (e.g. `0.5_1.0`), tissue,
cell line, disease, study year, library generation type, sequencing
platform, and many more (full filter list mirrors the catalogue's
faceted UI).

### dataset
```
igvfagent perturb-catalog dataset --dataset-id <id>
```
Full dataset record (single JSON document).

### dataset-rows
```
igvfagent perturb-catalog dataset-rows --modality mave \
    --dataset-id <id> --limit 500 --offset 0
```
Paginate through the **per-perturbation rows** inside one dataset
(variant- or gRNA-level effect scores).

### gsea
```
igvfagent perturb-catalog gsea --query BRCA1 --size 50
igvfagent perturb-catalog gsea --dataset-id <id>
```
Query the Perturb-seq GSEA endpoint for hallmark/pathway enrichment
tables.

### download
```
igvfagent perturb-catalog download --modality crispr-screen --dataset-id <id>
igvfagent perturb-catalog download --gsea --dataset-id <id>
```
Bulk download of one dataset (or all rows for a modality). Output
extension is auto-sniffed (`.csv.gz`, `.json`, `.zip`, …) and saved
under `Data/Perturbation/Downloads/<modality>/`.

### pipeline (one-shot)
```
igvfagent perturb-catalog pipeline --gene BRCA1
```
Runs summary + global gene search + modality-scoped searches for
each of MAVE / CRISPR-screen / Perturb-seq and writes a markdown
report under `Docs/Perturbation/<ts>_<gene>/`. The headline command
to call when the user asks "what perturbation data exists for
gene X?".

### write-playbook
```
igvfagent perturb-catalog write-playbook
```

## Cross-skill chaining

- `proteomics vampseq-analyze` (canonical MAVE / VAMP-seq scoresets
  from MaveDB) → `perturb-catalog search-modality --modality mave`
  for related Perturbation Catalogue datasets.
- `rnaseq deg` → `perturb-catalog search` for each significant gene
  to pull CRISPR / Perturb-seq evidence.
- `kg gene <SYM>` (IGVF Catalog) → `perturb-catalog pipeline --gene <SYM>`
  to attach perturbation evidence to a gene-centric workflow.
- `proteomics kg-visualize --gene <SYM>` → uses top hubs from the
  Perturb-seq GSEA terms to suggest functional cofactors.

## Output layout

```
Data/Perturbation/
    Searches/<ts>_<scope>.json
    Datasets/<dataset_id>.json
    Downloads/<modality>/<dataset_id>.<auto-sniffed-ext>
    Downloads/perturb-seq-gsea/<id>.csv.gz
Docs/Perturbation/<ts>_<label>/
    summary.json
    global_search.json
    {mave,crispr-screen,perturb-seq}_search.json
    report.md
```
