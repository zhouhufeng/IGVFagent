# OPERATIONS — Zou 2024 ChIP-Atlas 3.0 GATA1 hematopoietic

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run — fully online, ~20 s

```bash
bash Benchmarks/zou2024_chipatlas_gata1/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark zou2024_chipatlas_gata1
```

## What `run.sh` does

Three chained ChIP-Atlas operations:

```bash
# 1. Browse: how many GATA1 experiments for Blood
.venv/bin/igvfagent chipatlas list-antigens \
    --genome hg38 --ag-class "TFs and others" --cell-class Blood --limit 25

# 2. Search experiments by name + cell line
.venv/bin/igvfagent chipatlas search \
    --query "GATA1 K562" --genome hg38 --limit 10

# 3. POST to /download to assemble a unified all-peaks BED URL
.venv/bin/igvfagent chipatlas assemble-bed \
    --genome hg38 --ag-class "TFs and others" --antigen GATA1 \
    --cell-class Blood --qval 05 --label zou2024_chipatlas_gata1
```

## Where artefacts land

```
Docs/ChIPAtlas/<ts>_zou2024_chipatlas_gata1/
├── assemble.json         response from POST /download
└── ...
```

The browser and search calls print to stdout but also save responses
under `Docs/ChIPAtlas/`.

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `assemble.json` exists, non-empty | POST /download succeeded |

## Ground-truth spot-checks

| Signal | Expected (Zou et al. + canonical hematopoietic biology) |
|---|---|
| GATA1 experiments in Blood (`list-antigens`) | Should be among the top 20 by experiment count |
| K562 GATA1 search hits | Should return ≥ 5 reprocessed SRX experiments |
| Assembled BED URL | The `assemble.json` should contain a non-null `url` field pointing at a `.bed.gz` |
| Once the BED is downloaded — peaks overlap | HBB, HBA1, HBA2, KLF1 promoters + the Fulco/Engreitz +24/+58 kb GATA1 enhancers |

To verify the assembled BED actually contains the canonical hematopoietic
loci, download it and intersect:

```bash
# Extract the URL from assemble.json
URL=$(.venv/bin/python -c "import json; print(json.load(open('Docs/ChIPAtlas/<ts>_*/assemble.json'))['response']['url'])")

mkdir -p Data/Benchmarks/zou2024_chipatlas_gata1
curl -s "$URL" -o Data/Benchmarks/zou2024_chipatlas_gata1/gata1_blood.bed.gz

# Check overlap with canonical loci
zcat Data/Benchmarks/zou2024_chipatlas_gata1/gata1_blood.bed.gz | \
    awk '$1=="chr11" && $2 > 5240000 && $3 < 5280000' | head -3
# Should return rows in the HBB cluster

zcat Data/Benchmarks/zou2024_chipatlas_gata1/gata1_blood.bed.gz | \
    awk '$1=="chrX" && $2 > 48780000 && $3 < 48800000' | head -3
# Should return rows around the GATA1 locus
```

## Running through the UI

`chipatlas_list_antigens`, `chipatlas_search`, `chipatlas_assemble_bed`
ARE all registered. Paste:

```
Run the Zou 2024 ChIP-Atlas 3.0 GATA1 hematopoietic case study.

Step 1: chipatlas_list_antigens with genome="hg38",
        ag_class="TFs and others", cell_class="Blood", limit=25.

Step 2: chipatlas_search with query="GATA1 K562", genome="hg38", limit=10.

Step 3: chipatlas_assemble_bed with genome="hg38",
        ag_class="TFs and others", antigen="GATA1", cell_class="Blood",
        qval="05", label="zou2024_chipatlas_gata1".

Report:
  - GATA1's rank in the antigen browser for Blood (how many experiments?)
  - top 5 search hits and their SRX accessions
  - the assembled-BED URL (or null if none)
  - cite Zou/Ohta/Oki Nucleic Acids Res 2024 (doi:10.1093/nar/gkae358)
```

UI sidebar: max iterations = 12, temperature = 0.0.

## Optional: enrichment of the GATA1 hematopoietic target gene list

After downloading the assembled BED and intersecting with TSSs to get
target genes:

```bash
# Build a top-N target gene list (this is a sketch; real workflow needs bedtools + a TSS annotation)
echo -e "HBB\nHBA1\nHBA2\nKLF1\nGATA1\nLMO2\nMYB" > Data/Benchmarks/zou2024_chipatlas_gata1/top_genes.txt

.venv/bin/igvfagent enrich ora \
    --genes Data/Benchmarks/zou2024_chipatlas_gata1/top_genes.txt \
    --label zou2024_gata1_targets_enrichment
```

Expected: dominant GO terms in "hemoglobin complex", "erythrocyte
differentiation", "heme biosynthesis"; dominant Reactome pathway:
"Erythrocytes take up oxygen and release carbon dioxide".

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `list-antigens` returns 0 GATA1 entries for Blood | Cell-class label mismatch | Try `--cell-class "Blood vessel"` or browse with `chipatlas list-cell-types` first |
| `assemble.json` has `url: null` | No matching experiments for the exact (ag × cellClass × qval) tuple | Loosen the q-value: try `--qval 10` (less strict peaks) |
| `chipatlas search "GATA1 K562"` returns experiments without K562 | ChIP-Atlas search is full-text, may match the abstract not the cell type | Use the more specific browse via `list-cell-types --ag-class "TFs and others"` then drill |
| Polite-rate limit (1 req/s) makes things slow | Default is 1 rps | `export CHIPATLAS_MIN_SECONDS=0.3` to go faster (still polite) |

## License + provenance

* **Paper data**: ChIP-Atlas (NBDC/DBCLS) — public, CC-BY 4.0 for data; fetch-only.
* **Code**: IGVFagent Apache-2.0 over upstream MIT (inutano/chip-atlas).
* **Citation**: Zou A, Ohta T, Oki S. *Nucleic Acids Res* **52**: W45–W53 (2024).
  doi:10.1093/nar/gkae358
