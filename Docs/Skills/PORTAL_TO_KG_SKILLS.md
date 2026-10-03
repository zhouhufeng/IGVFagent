# Skill: IGVF Portal → local Knowledge Graph ETL

Pulls unstructured IGVF Portal entities (AnalysisSets, MeasurementSets, Samples, Donors, Files, FileSets) by tissue / gene / assay / lab filters; annotates them against genes and variants; and persists the result to a local SQLite-backed knowledge graph that mirrors the IGVF Catalog ArangoDB schema (nodes + edges + properties + provenance).

## Local store

  Data/KG/local_kg.sqlite     SQLite database (gitignored)
  Data/KG/Export/<timestamp>/ arangoimport-compatible JSONL exports
  Data/Cache/KGLocal/         per-run caches

All endpoint URLs and credentials are resolved through `Scripts/_endpoints.py` env-var overrides. No URLs, passwords, or API keys are written to source.

## Subcommands

### `pull` — fetch Portal entities by filter

```bash
# AnalysisSets profiling macrophages, hydrated + walked one hop
python3 Scripts/portal_to_kg_skill.py pull \
    --type AnalysisSet --tissue macrophage --limit 50 --depth 1

# Multiple types, single tissue
python3 Scripts/portal_to_kg_skill.py pull \
    --types AnalysisSet MeasurementSet --tissue brain --limit 25

# By assay only
python3 Scripts/portal_to_kg_skill.py pull \
    --type AnalysisSet --assay 'Parse SPLiT-seq' --limit 100 --depth 1

# By gene name (free-text searchTerm)
python3 Scripts/portal_to_kg_skill.py pull \
    --type AnalysisSet --gene APOE --limit 25 --depth 1
```

Each Portal entity becomes a node. With `--depth 1` the skill follows `samples`, `donors`, `files`, `input_file_sets`, `workflows`, and `documents` cross-references and ingests them as linked nodes (with appropriate edge labels). Use `--skip-hydrate` to avoid the per-record full hydration call when you only want headline metadata.

### `annotate` — text-mine for genes and variants

```bash
python3 Scripts/portal_to_kg_skill.py annotate
python3 Scripts/portal_to_kg_skill.py annotate \
    --gene-list path/to/HGNC_symbols.txt
python3 Scripts/portal_to_kg_skill.py annotate --skip-catalog-confirm
```

Walks every Portal node in the local KG, scans its description / summary / aliases / title / name fields for gene-symbol tokens, rsID strings, and SPDI/HGVS variant identifiers. Each candidate gene token is confirmed against the IGVF Catalog (`/api/genes?name=…`) and cached in the `gene_cache` table so subsequent runs are free. Confirmed Gene and Variant nodes are added with `mentions_gene` / `mentions_variant` edges back to the Portal node.

Pass `--gene-list` to constrain annotation to a curated set (useful when you only care about a specific GWAS / disease / curated panel). `--skip-catalog-confirm` retains every uppercase-token candidate as an unconfirmed Gene node — fast but noisier.

### `enrich` — hydrate genes from the IGVF Catalog

```bash
# All genes in local KG -> pull variants/cCREs/transcripts/etc.
python3 Scripts/portal_to_kg_skill.py enrich --limit 25

# Just specific symbols
python3 Scripts/portal_to_kg_skill.py enrich --symbols APOE,LDLR,PCSK9
```

For each Gene node, fetches the Catalog's `variants`, `transcripts`, `proteins`, `regulatory_elements`, `diseases`, `pathways`, and `coding_variant_scores` relations and stores them as nodes and edges in the local KG (with `source = igvf_catalog`).

### `query` — read interface

```bash
python3 Scripts/portal_to_kg_skill.py query --gene APOE --limit 50
python3 Scripts/portal_to_kg_skill.py query --tissue macrophage --limit 30
python3 Scripts/portal_to_kg_skill.py query --node-id analysis_sets/IGVFDS3222WCZH
```

### `stats` — counts and recent ingestion summary

```bash
python3 Scripts/portal_to_kg_skill.py stats
```

### `export-aql` — emit ArangoDB-compatible JSONL

```bash
python3 Scripts/portal_to_kg_skill.py export-aql
```

Writes one JSONL file per node collection (`nodes_genes.jsonl`, `nodes_variants.jsonl`, …) plus `edges.jsonl` under `Data/KG/Export/<timestamp>/`. Each record has the ArangoDB convention `_key` (and `_from`/`_to` for edges) so a future push to the actual IGVF Catalog ArangoDB is just `arangoimport --type jsonl …` per file.

### `export-cytoscape` — graph visualization

```bash
python3 Scripts/portal_to_kg_skill.py export-cytoscape --limit 500
```

## Recommended end-to-end flow

```bash
# 1. Pull a tissue-specific corpus from the Portal
python3 Scripts/portal_to_kg_skill.py pull \
    --type AnalysisSet --tissue macrophage --limit 100 --depth 1
python3 Scripts/portal_to_kg_skill.py pull \
    --type MeasurementSet --tissue macrophage --limit 100 --depth 1

# 2. Mine descriptions for gene/variant mentions, confirm against Catalog
python3 Scripts/portal_to_kg_skill.py annotate

# 3. Hydrate every confirmed gene from the Catalog
python3 Scripts/portal_to_kg_skill.py enrich --limit 25

# 4. Inspect the result
python3 Scripts/portal_to_kg_skill.py stats
python3 Scripts/portal_to_kg_skill.py query --gene APOE

# 5. (Optional) Export for a future push to IGVF Catalog ArangoDB
python3 Scripts/portal_to_kg_skill.py export-aql
```

## How this chains with other skills

- After `pull`/`annotate`, the resulting AnalysisSet → Gene edges feed directly into `Scripts/kg_traversal_skill.py` for comprehensive remote-Catalog context on each implicated gene.
- The `single_cell_data_skills.py`, `splitseq_pipeline.py`, and `multiome_10x_pipeline.py` pipelines can ingest the AnalysisSet manifest from `query` to drive downstream analysis.
- `reference_skill.py validate` can be run against the local Gene/Variant nodes to surface prior literature.
- The internal Plan → Action → Results → Evaluation orchestrator uses this skill as the **persistence and indexing layer** for Portal-side evidence so multi-step plans don't re-fetch the same metadata across sessions.