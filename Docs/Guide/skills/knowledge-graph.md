# Knowledge graph, warehouse and networks

The local knowledge graph that grows with use, the IGVF Catalog mirror, graph traversal, the DuckDB warehouse and network integration.

**On this page**

- [The growing local knowledge graph + database](#the-growing-local-knowledge-graph--database)
- [Local IGVF KG mirror (Arango → DuckDB)](#local-igvf-kg-mirror-arango--duckdb)
- [Knowledge Graph traversal](#knowledge-graph-traversal)
- [Portal → local KG ETL](#portal--local-kg-etl)
- [Integrated data warehouse (DuckDB Silver tier)](#integrated-data-warehouse-duckdb-silver-tier)
- [Network integration (clean-room MILP — CARNIVAL + Steiner)](#network-integration-clean-room-milp--carnival--steiner)

## The growing local knowledge graph + database

Every time IGVFagent **downloads data or analyzes/processes raw data**, a little
more of a *local* knowledge graph and database accumulates — so the agent gets
faster and more self-contained the more you use it, and repeat queries can be
answered from local data. This is the **core default mechanism**, wired into
both the agent loop (every tool call) and the CLI (every direct skill run
auto-harvests its new outputs).

- **Knowledge Graph** — `Data/KG/local_kg.sqlite` (`nodes` + `edges`: genes,
  variants, regions, datasets, studies, analyses, and the relations between
  them — the same graph the `portal-kg` skill builds).
- **Database** — the DuckDB warehouse at `Data/Warehouse/igvf.duckdb`.

```bash
igvfagent localstore stats      # nodes / edges / downloads / analyses logged
igvfagent localstore harvest    # scan Docs/ + downloads and ingest anything new
```

Growth is idempotent (deterministic upserts + a harvest ledger) and safe under
concurrent agent/CLI writers (WAL + busy-timeout). Disable with
`IGVF_LOCALSTORE=0`.

## Local IGVF KG mirror (Arango → DuckDB)

The full IGVF Catalog Knowledge Graph in Arango is ~2 TB across 58
collections (25 document + 33 edge). The `kg-mirror` skill streams every
collection except the two planet-scale `variants` tables and the two wide-doc
`genes_coding_variants_scores*` tables (which embed ~80 MB per-variant
score matrices that consistently time out the AQL cursor). Streams via
the read-only AQL cursor API, persists each as zstd-compressed Parquet
shards, then registers a DuckDB warehouse with one view per collection. Lets `igvfagent kg ...` and the
downstream skills run offline against the cached copy. See
[`Docs/Skills/KG_MIRROR_SKILL.md`](../../Skills/KG_MIRROR_SKILL.md).

```bash
# Requires Arango read credentials in .env:
#   IGVF_ARANGO_USER=guest
#   IGVF_ARANGO_PASSWORD=guestigvfcatalog

# 1. Inventory the upstream KG
igvfagent kg-mirror inventory

# 2. Mirror everything except the two planet-scale variants tables (default)
igvfagent kg-mirror pull-all

# 3. Or pull one collection at a time (resumable — state on disk)
igvfagent kg-mirror pull --collection genes
igvfagent kg-mirror pull --collection coding_variants --batch-size 10000

# 4. Cap to small/medium collections only (≤ 10 GB each)
igvfagent kg-mirror pull-all --max-collection-bytes 10000000000

# 5. Register Parquet shards as DuckDB views
igvfagent kg-mirror register
igvfagent kg-mirror verify
```

Storage layout:
- `Data/Warehouse/KG/<collection>/<NNNN>.parquet` — ZSTD-compressed shards
- `Data/Warehouse/KG/_state/<collection>.json` — resume cursor
- `Data/Warehouse/igvf_kg_mirror.duckdb` — one `kg_<collection>` view each

## Knowledge Graph traversal

The orchestrator-friendly **comprehensive context** tool. Starts from a
single entity (gene, variant, or genomic region) and iteratively walks the
IGVF Catalog Knowledge Graph (`api.catalogkg.igvf.org` + the underlying
ArangoDB at `db.catalog.igvf.org`), assembling a unified evidence pack that
includes direct neighbors, second-degree relations, and optional cross-skill
enrichment (FAVOR, enhancer-gene linkage, IGVF single-cell datasets, prior
literature). One CLI call → per-relation manifests + JSON evidence pack +
markdown report.

```bash
# Comprehensive APOE traversal — variants, cCREs, transcripts, proteins,
# diseases, pathways, plus enhancer-gene linkage, candidate single-cell
# datasets (matched via the gene's eQTL biological contexts), and prior
# literature
python3 Scripts/kg_traversal_skill.py gene APOE \
  --depth 2 --limit 50 --max-variants 25 --subvariant-limit 10 \
  --call-favor --call-linkage --call-singlecell --call-literature \
  --literature-context Alzheimer cardiovascular --label apoe_full

# Variant-centric: one rsID -> linked genes, cCREs, phenotypes, predictions
python3 Scripts/kg_traversal_skill.py variant rs429358 \
  --call-favor --call-literature --label apoe_e4_variant

# Region-centric: genes + cCREs + linkage in window
python3 Scripts/kg_traversal_skill.py region chr19:44903000-44912000 \
  --call-favor --label apoe_locus

# Direct AQL pass-through
python3 Scripts/kg_traversal_skill.py aql \
  'FOR g IN genes FILTER g.name == "APOE" RETURN g'
python3 Scripts/kg_traversal_skill.py write-playbook
```

Outputs land under `Docs/KGTraversal/<timestamp>_<label>/` — a markdown
report with sectioned evidence per relation type, a complete
`evidence_pack.json`, and one CSV per relation under `Manifests/`. The
manifest CSVs feed directly into the variant-annotation, advanced
variant-analysis, single-cell, SPLiT-seq, and reference skills, which is
exactly the cross-skill composition that the internal Plan → Action →
Results → Evaluation orchestrator chains together.

## Portal → local KG ETL

A persistence and indexing layer that ingests unstructured IGVF Portal
entities (AnalysisSets, MeasurementSets, Samples, Donors, Files,
FileSets, Documents) into a **local SQLite-backed knowledge graph** that
mirrors the IGVF Catalog ArangoDB schema (nodes + edges + properties +
provenance). The local KG can be queried, annotated, enriched from the
remote Catalog, and later exported as `arangoimport`-compatible JSONL for
a future push to the IGVF Catalog. Endpoint URLs and credentials are
resolved through `_endpoints.py` env-var overrides — nothing sensitive is
written to source.

```bash
# 1. Pull a tissue / assay / lab corpus from the Portal (one hop expansion)
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --tissue macrophage --limit 100 --depth 1
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --assay 'Parse SPLiT-seq' --limit 100 --depth 1
python3 Scripts/portal_to_kg_skill.py pull \
  --type AnalysisSet --gene APOE --limit 25 --depth 1

# 2. Text-mine for gene + variant mentions, confirm against the Catalog
python3 Scripts/portal_to_kg_skill.py annotate

# 3. Hydrate confirmed genes from the Catalog (variants, cCREs, diseases…)
python3 Scripts/portal_to_kg_skill.py enrich --symbols APOE,TREM2,LDLR --limit 25

# 4. Read the local KG
python3 Scripts/portal_to_kg_skill.py stats
python3 Scripts/portal_to_kg_skill.py query --gene APOE --limit 50
python3 Scripts/portal_to_kg_skill.py query --node-id analysis_sets/IGVFDS3222WCZH

# 5. Export to ArangoDB-compatible JSONL for an eventual Catalog push
python3 Scripts/portal_to_kg_skill.py export-aql
python3 Scripts/portal_to_kg_skill.py export-cytoscape --limit 500
python3 Scripts/portal_to_kg_skill.py write-playbook
```

Local store: `Data/KG/local_kg.sqlite` (gitignored under `Data/*`).
Exports: `Data/KG/Export/<timestamp>/nodes_<collection>.jsonl` plus
`edges.jsonl`, ready for `arangoimport --type jsonl …`.

## Integrated data warehouse (DuckDB Silver tier)

The central data warehouse — the "Silver tier" of the integrated
data layer described in
[`Docs/Architecture/INTEGRATED_DATA_LAYER.md`](../../Architecture/INTEGRATED_DATA_LAYER.md).
Every other IGVFagent skill becomes a producer that lands QC'd rows
into canonical entity / edge / measurement tables in a single
**DuckDB** database at `Data/Warehouse/igvf.duckdb`, ready to feed
downstream embedding extraction + foundation-model training.
Playbook:
[`Docs/Skills/WAREHOUSE_SKILLS.md`](../../Skills/WAREHOUSE_SKILLS.md).

```bash
# 1) Initialise the schema (entities + edges + measurements + provenance)
igvfagent warehouse init

# 2) Pull every available producer into the warehouse
igvfagent warehouse ingest --source all
# or one at a time:
igvfagent warehouse ingest --source proteomics-kg
igvfagent warehouse ingest --source perturb-catalog
igvfagent warehouse ingest --source multiseq
igvfagent warehouse ingest --source mavedb-vampseq

# 3) Stats and cross-skill SQL queries
igvfagent warehouse stats
igvfagent warehouse query \
    "SELECT gene_id, COUNT(*) AS n, AVG(score) FROM vampseq_scores GROUP BY gene_id"
```

DuckDB is chosen over SQLite / Postgres / BigQuery because it is
embedded (no server), columnar (10–100× faster on the analytical
joins this workload runs), reads / writes Parquet natively (so the
Gold-tier embedding tables stay first-class), and ships with vector
functions. Live smoke test (laptop, ~6 s): **32K PPI edges**,
**102K MULTI-seq cells**, **24K VAMP-seq scores**, **28 IGVF
protein-evidence datasets** ingested from five producers. Bulk
`posteriors.csv` etc. are auto-excluded; the tracked outputs are
intentionally lightweight (~14 MB on disk).

## Network integration (clean-room MILP — CARNIVAL + Steiner)

The **integration layer** of the IGVF data warehouse. Implements two
context-specific-subnetwork methods from scratch in cvxpy:

  * **CARNIVAL** — given a signed prior-knowledge graph plus signed
    perturbations + signed measurements, find the minimum-cost
    upstream subnetwork whose vertex signs match the data.
  * **Prize-collecting Steiner** — given per-gene prizes and a PPI,
    find the connected subgraph maximising (prizes − edge costs).

The math is the **CORNETO formulation** ([Rodriguez-Mier et al., *Nat
Mach Intell* 2025](https://www.nature.com/articles/s42256-025-01069-9))
re-implemented in original Apache-2 cvxpy — **no CORNETO source is
imported or vendored.** Algorithms are not copyrightable; source code
is. The full algorithm specification, attribution, and citation list
live in
[`Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`](../../Architecture/INTEGRATION_LAYER_REFERENCE.md).
Playbook:
[`Docs/Skills/NETWORK_INTEGRATION_SKILLS.md`](../../Skills/NETWORK_INTEGRATION_SKILLS.md).

```bash
pip install 'igvfagent[network]'    # cvxpy + SCIP (free MILP)

# 1) End-to-end self-test: synthetic EGFR → SOS1 → … → MYC cascade
igvfagent network demo --beta 0.05 --solver SCIP
# → CARNIVAL recovers all 6 cascade edges and writes them to the
#   warehouse with upstream='network:demo'.

# 2) Materialise a signed PKN from the proteomics KG
igvfagent network pkn-from-kg --label reactome_pkn

# 3) CARNIVAL on Perturb-seq-style inputs
#    perts.csv: gene,sign   (sign ∈ {-1,+1})
#    degs.csv:  gene,score  (signed log2 fold-change)
igvfagent network carnival \
    --perturbations perts.csv \
    --measurements  degs.csv \
    --pkn-limit 5000 --solver SCIP --label perturb_seq_demo

# 4) Prize-collecting Steiner tree (VAMP-seq abundance prizes on PPI)
igvfagent network steiner --terminals vamp_prizes.csv \
    --pkn-limit 5000 --edge-cost 1.0
```

License boundary: **Apache-2 throughout**. No GPL runtime dependencies.
SCIP (PicoMILP backend) is BSD-style, free for any use. Selected
subnetworks flow back into the central DuckDB warehouse as `edges`
rows tagged with `upstream='network:carnival:<label>'` or
`upstream='network:steiner:<label>'` so downstream embedding / 
foundation-model training treats them like any other evidence stream.
Three agent tools registered: `network_demo`, `network_carnival`,
`network_steiner`.

> **Current limitations (v1).** The prior-knowledge network is protein-only, built from undirected PPI with every sign set to +1, so variants and regulatory elements enter only as gene-level prizes. CARNIVAL does not yet enforce acyclicity or reachability from the perturbations, so a sign-consistent cycle can explain a measurement on its own. The Steiner model has no connectivity constraint yet, so it selects terminals but no edges. Both follow-ups are specified in [`INTEGRATION_LAYER_REFERENCE.md`](../../Architecture/INTEGRATION_LAYER_REFERENCE.md). Treat results as exploratory until they land.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
