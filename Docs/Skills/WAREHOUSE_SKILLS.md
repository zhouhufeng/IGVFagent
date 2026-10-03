# Skill: Integrated data warehouse

The Silver tier of the IGVF integrated data layer (full architecture in
``Docs/Architecture/INTEGRATED_DATA_LAYER.md``). A DuckDB warehouse at
``Data/Warehouse/igvf.duckdb`` that every other IGVFagent skill lands
QC'd rows into, ready for downstream embedding extraction and
foundation-model training.

## Subcommands

### init
```
igvfagent warehouse init
```
Creates the canonical schema (entities + edges + measurements +
provenance) at ``Data/Warehouse/igvf.duckdb``. Idempotent.

### ingest
```
igvfagent warehouse ingest --source proteomics-kg
igvfagent warehouse ingest --source portal-kg
igvfagent warehouse ingest --source perturb-catalog
igvfagent warehouse ingest --source multiseq
igvfagent warehouse ingest --source mavedb-vampseq
igvfagent warehouse ingest --source all
```
Each source adapter reads an existing skill's local outputs and upserts
into the canonical entity / edge / measurement tables, writing a
``runs`` provenance row.

### stats
```
igvfagent warehouse stats
```
Row counts per table and a tail of the ``runs`` provenance log.

### query
```
igvfagent warehouse query "SELECT COUNT(*) FROM edges GROUP BY relation"
```
Ad-hoc SQL against the warehouse (read-only). Use this to join across
skills — e.g. "VAMP-seq variants in genes that are perturbation-screen
hits, weighted by their PPI neighbourhood size."

### sources
```
igvfagent warehouse sources
```
List the registered producers. Add a new producer by writing a new
function in ``warehouse_skill.py`` and adding it to the ``PRODUCERS``
dict.

## Adding a new producer

```python
def ingest_my_skill(con) -> dict:
    """One-line summary used by `sources`."""
    # Read whatever the upstream skill writes locally
    # Upsert into entity / edge / measurement tables
    # Return {table_name: row_count_added}
    return {"my_table": n}

PRODUCERS["my-skill"] = ingest_my_skill
```

## Why this design

- DuckDB (embedded, columnar, Parquet-native, vector-ready) is 10–100×
  faster than SQLite for the wide analytical joins this workload runs.
- Every skill's outputs stay in their original location (no double
  storage); the warehouse holds canonical references and the small
  structured fields.
- Provenance is enforced at the ingest layer — every run is logged,
  every row carries its `upstream` source.
- Read-only ad-hoc SQL is the universal cross-skill query surface.

See ``Docs/Architecture/INTEGRATED_DATA_LAYER.md`` for the full
three-tier (Bronze / Silver / Gold) architecture, the foundation-model
plan, and the implementation roadmap.
