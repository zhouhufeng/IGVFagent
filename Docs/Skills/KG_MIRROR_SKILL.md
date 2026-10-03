# Skill: Local IGVF Knowledge Graph Mirror

Streams Arango collections via the read-only AQL cursor API and persists
them as Parquet shards on disk, then registers a DuckDB warehouse with
one view per collection. Lets `igvfagent kg ...` and downstream skills
run **offline** against the cached copy.

## Defaults

- Skip list: `variants` (~944 GB doc) and `variants_variants` (~531 GB edge).
  These two account for ~75% of the upstream KG. Override with
  `--skip ""` or `--include-giants`.
- Storage:
  - `Data/Warehouse/KG/<collection>/<NNNN>.parquet` — ZSTD-compressed shards
  - `Data/Warehouse/KG/_state/<collection>.json` — resume cursor
  - `Data/Warehouse/igvf_kg_mirror.duckdb` — DuckDB with `kg_<collection>` views

## Commands

```bash
# 1. Inventory: list all collections with doc count + bytes
igvfagent kg-mirror inventory

# 2. Mirror everything except the two giants
igvfagent kg-mirror pull-all

# 3. Mirror one collection at a time
igvfagent kg-mirror pull --collection genes
igvfagent kg-mirror pull --collection coding_variants --batch-size 10000

# 4. Resume a long-running pull (state is on disk)
igvfagent kg-mirror pull --collection variants_proteins   # resumes automatically

# 5. Cap to small/medium collections only (under 10 GB each)
igvfagent kg-mirror pull-all --max-collection-bytes 10000000000

# 6. Register Parquet shards as DuckDB views
igvfagent kg-mirror register
igvfagent kg-mirror verify
```

## Querying the local mirror

```python
import duckdb
con = duckdb.connect("Data/Warehouse/igvf_kg_mirror.duckdb", read_only=True)
con.sql("SELECT count(*) FROM kg_genes").show()
con.sql("SELECT * FROM kg_proteins LIMIT 5").show()
# Edge collections: same naming convention — kg_genes_pathways, kg_variants_genes, ...
```

## Schema notes

- Every row is a flat dict. Scalar fields (str/int/float/bool) stay as
  scalars. Anything nested (list, dict) is `json.dumps`'d into a string
  column so the Parquet schema stays consistent across batches.
- `_id`, `_key`, `_from`, `_to` are preserved verbatim. Edge collections
  always carry `_from` and `_to` references like `genes/<gene_id>`.

## License

Apache-2.0. Uses only stdlib + pyarrow (Apache-2) + duckdb (MIT). The
Arango HTTP API is hit via basic-auth with the read-only `guest` account.
