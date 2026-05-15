"""Integrated data-warehouse skill — central DuckDB layer.

The "Silver" tier of the IGVF integrated data layer described in
``Docs/Architecture/INTEGRATED_DATA_LAYER.md``. Every other IGVFagent
skill becomes a *producer* that lands QC'd rows in canonical tables;
this skill owns the schema, ingest plumbing, provenance log, and
analytical query surface that downstream embedding / foundation-model
training will consume.

Why DuckDB
----------
DuckDB is embedded (no server), columnar (10–100× faster than SQLite
for analytical joins), and reads / writes Parquet natively — which
means the Gold-tier embedding files are first-class citizens, queryable
without ETL. Ships with vector functions for cosine similarity. The
right shape for "many small tables joined wide for downstream ML".

Subcommands
-----------
    init             Create the schema at Data/Warehouse/igvf.duckdb
    ingest           Pull a named source into the warehouse
                       (proteomics-kg | portal-kg | perturb-catalog |
                        multiseq | mavedb-vampseq | <skill-shortname>)
    stats            Row counts per table + provenance summary
    query            Run an ad-hoc SQL query
    sources          List producers registered in this skill
    write-playbook   Write Docs/Skills/WAREHOUSE_SKILLS.md
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
WAREHOUSE_DIR = ROOT / "Data" / "Warehouse"
WAREHOUSE_DB = WAREHOUSE_DIR / "igvf.duckdb"
DOCS_DIR = ROOT / "Docs" / "Warehouse"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "WAREHOUSE_SKILLS.md"

logger = logging.getLogger("warehouse")


def _duckdb():
    """Soft-import DuckDB so we error helpfully when missing."""
    try:
        import duckdb  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "warehouse requires DuckDB. Install with:\n"
            "  pip install 'duckdb>=0.10'\n"
            f"Original: {exc}"
        )
    return duckdb


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in (s or "run"))


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"warehouse_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def mkdirs() -> None:
    for d in (WAREHOUSE_DIR, DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

# Compact starter schema — full version lives in
# Docs/Architecture/INTEGRATED_DATA_LAYER.md. Adding columns here is a
# non-breaking change because DuckDB supports ALTER TABLE ADD COLUMN.
SCHEMA_SQL = """
-- Entities -----------------------------------------------------------
CREATE TABLE IF NOT EXISTS genes (
    gene_id          TEXT PRIMARY KEY,
    symbol           TEXT,
    chrom            TEXT,
    start_bp         BIGINT,
    end_bp           BIGINT,
    strand           VARCHAR(1),
    biotype          TEXT,
    taxon            INTEGER,
    uniprot_acs      TEXT[],
    source_priority  TEXT,
    ingested_at      TIMESTAMP DEFAULT now()
);

CREATE TABLE IF NOT EXISTS variants (
    variant_id       TEXT PRIMARY KEY,
    chrom            TEXT,
    pos_bp           BIGINT,
    ref              TEXT,
    alt              TEXT,
    rsid             TEXT,
    type             TEXT,
    population_af    DOUBLE,
    cadd             DOUBLE,
    favor_score      DOUBLE,
    consequence      TEXT,
    overlaps_re      TEXT[],
    taxon            INTEGER,
    ingested_at      TIMESTAMP DEFAULT now()
);

CREATE TABLE IF NOT EXISTS proteins (
    uniprot_ac       TEXT PRIMARY KEY,
    symbol           TEXT,
    gene_id          TEXT,
    length_aa        INTEGER,
    organism_taxid   INTEGER,
    sequence_md5     TEXT,
    ingested_at      TIMESTAMP DEFAULT now()
);

CREATE TABLE IF NOT EXISTS regulatory_elements (
    re_id            TEXT PRIMARY KEY,
    chrom            TEXT,
    start_bp         BIGINT,
    end_bp           BIGINT,
    classification   TEXT,
    assembly         TEXT,
    score            DOUBLE,
    source           TEXT,
    ingested_at      TIMESTAMP DEFAULT now()
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id        TEXT PRIMARY KEY,
    biosample        TEXT,
    tissue           TEXT,
    cell_line        TEXT,
    sex              TEXT,
    age              TEXT,
    treatment        TEXT,
    disease          TEXT,
    donor_id         TEXT
);

CREATE TABLE IF NOT EXISTS datasets (
    dataset_id       TEXT PRIMARY KEY,
    upstream_source  TEXT NOT NULL,
    accession        TEXT,
    title            TEXT,
    assay            TEXT,
    modality         TEXT,
    biosample_id     TEXT,
    year             INTEGER,
    n_samples        INTEGER,
    pubmed_id        TEXT,
    license          TEXT,
    download_url     TEXT,
    bronze_path      TEXT,
    ingested_at      TIMESTAMP DEFAULT now()
);

CREATE TABLE IF NOT EXISTS cells (
    cell_id          TEXT PRIMARY KEY,
    dataset_id       TEXT,
    sample_id        TEXT,
    n_genes          INTEGER,
    n_counts         INTEGER,
    pct_counts_mt    DOUBLE,
    leiden           TEXT,
    cell_type_label  TEXT,
    multiseq_call    TEXT,
    ingested_at      TIMESTAMP DEFAULT now()
);

-- Edges (one table for all relations) --------------------------------
CREATE TABLE IF NOT EXISTS edges (
    src_type   TEXT NOT NULL,
    src_id     TEXT NOT NULL,
    dst_type   TEXT NOT NULL,
    dst_id     TEXT NOT NULL,
    relation   TEXT NOT NULL,
    score      DOUBLE,
    pubmed_id  TEXT,
    evidence   TEXT,
    upstream   TEXT NOT NULL,
    UNIQUE(src_type, src_id, dst_type, dst_id, relation, upstream)
);

-- Measurements (one table per modality) ------------------------------
CREATE TABLE IF NOT EXISTS vampseq_scores (
    measurement_id   BIGINT,
    gene_id          TEXT,
    protein_ac       TEXT,
    aa_position      INTEGER,
    aa_wt            VARCHAR(1),
    aa_alt           VARCHAR(1),
    score            DOUBLE,
    se               DOUBLE,
    abundance_class  INTEGER,
    dataset_id       TEXT,
    upstream         TEXT
);

CREATE TABLE IF NOT EXISTS perturbseq_effects (
    measurement_id   BIGINT,
    perturbed_gene   TEXT,
    target_gene      TEXT,
    log2fc           DOUBLE,
    padj             DOUBLE,
    cell_type        TEXT,
    dataset_id       TEXT,
    upstream         TEXT
);

-- Provenance ---------------------------------------------------------
CREATE TABLE IF NOT EXISTS runs (
    run_id        TEXT PRIMARY KEY,
    skill         TEXT,
    subcommand    TEXT,
    args_json     JSON,
    started_at    TIMESTAMP,
    ended_at      TIMESTAMP,
    rows_emitted  BIGINT,
    target_table  TEXT,
    success       BOOLEAN
);

CREATE TABLE IF NOT EXISTS versions (
    source        TEXT PRIMARY KEY,
    version       TEXT,
    upstream_url  TEXT,
    sha256        TEXT,
    fetched_at    TIMESTAMP
);

CREATE TABLE IF NOT EXISTS qc (
    table_name    TEXT,
    metric        TEXT,
    value         DOUBLE,
    passed        BOOLEAN,
    threshold     DOUBLE,
    run_id        TEXT
);
"""


def open_warehouse(read_only: bool = False):
    """Open the DuckDB warehouse, creating directories as needed."""
    duckdb = _duckdb()
    mkdirs()
    return duckdb.connect(str(WAREHOUSE_DB), read_only=read_only)


def init_schema() -> None:
    con = open_warehouse()
    con.execute(SCHEMA_SQL)
    con.close()
    logger.info("Initialised schema at %s", WAREHOUSE_DB)


def _record_run(con, *, skill: str, subcommand: str, args: dict,
                 rows: int, target_table: str, success: bool) -> str:
    # Use microsecond precision so multiple ingests within the same
    # second don't collide on the run_id primary key.
    rid = (f"{time.strftime('%Y%m%d_%H%M%S')}"
            f"_{int((time.time() % 1) * 1_000_000):06d}"
            f"_{safe_label(skill)}_{safe_label(subcommand)}_"
            f"{safe_label(target_table)[:20]}")
    con.execute(
        "INSERT INTO runs VALUES (?, ?, ?, ?, ?, now(), ?, ?, ?)",
        [rid, skill, subcommand, json.dumps(args, default=str),
         time.strftime("%Y-%m-%d %H:%M:%S"), rows, target_table, success],
    )
    return rid


# ---------------------------------------------------------------------------
# Source adapters — each one reads an existing skill's local outputs and
# upserts rows into the warehouse. Add a new function here when you want
# to wire in a new producer.
# ---------------------------------------------------------------------------


def ingest_proteomics_kg(con) -> dict:
    """Pull from Data/Proteomics/KG/proteomics.sqlite → warehouse."""
    src = ROOT / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"
    if not src.exists():
        logger.warning("proteomics KG not found at %s — skip", src)
        return {"proteomics_kg": 0}
    sqlite_con = sqlite3.connect(f"file:{src}?mode=ro", uri=True)
    counts: "dict[str, int]" = {}

    # Proteins from id_map (if populated) — best canonical source
    rows = sqlite_con.execute(
        "SELECT uniprot_ac, entrez, ensembl, symbol FROM id_map "
        "WHERE uniprot_ac IS NOT NULL LIMIT 1000000"
    ).fetchall()
    if rows:
        con.executemany("""
            INSERT OR REPLACE INTO proteins
              (uniprot_ac, symbol, gene_id, organism_taxid)
            VALUES (?, ?, ?, ?)
        """, [(r[0], r[3], r[2], 9606) for r in rows])
        counts["proteins"] = len(rows)

    # Interactions → edges
    inter_rows = sqlite_con.execute("""
        SELECT id_a, id_b, source, id_type, detection_method, evidence_type,
               pubmed_id, confidence_score
        FROM interactions
    """).fetchall()
    if inter_rows:
        edges = [
            ("protein" if r[3] == "uniprot" else r[3],
              r[0],
              "protein" if r[3] == "uniprot" else r[3],
              r[1],
              "interacts_with",
              r[7], r[6], r[5], r[2])
            for r in inter_rows
        ]
        con.executemany("""
            INSERT OR IGNORE INTO edges
              (src_type, src_id, dst_type, dst_id, relation,
               score, pubmed_id, evidence, upstream)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, edges)
        counts["edges_ppi"] = len(edges)

    # IGVF evidence files → datasets table
    ev = sqlite_con.execute("""
        SELECT file_accession, dataset_accession, assay_title, lab,
               content_type, description, url
        FROM igvf_evidence
    """).fetchall()
    if ev:
        ds = [(f"igvf:{r[0]}", "igvf-portal", r[0], r[2], r[2], None,
                None, None, None, None, r[6], None) for r in ev]
        con.executemany("""
            INSERT OR REPLACE INTO datasets
              (dataset_id, upstream_source, accession, title, assay,
               modality, biosample_id, year, n_samples, pubmed_id,
               download_url, bronze_path)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, ds)
        counts["datasets_igvf_protein"] = len(ds)

    sqlite_con.close()
    return counts


def ingest_perturb_catalog(con) -> dict:
    """Pull dataset metadata from Data/Perturbation/Searches/*.json."""
    base = ROOT / "Data" / "Perturbation" / "Searches"
    if not base.exists():
        return {"perturb_catalog_datasets": 0}
    n = 0
    for jf in base.glob("*modality_*.json"):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        for ds in (d.get("datasets") or []):
            meta = (ds.get("dataset") or ds.get("metadata") or {}) \
                if isinstance(ds, dict) else {}
            ds_id = (meta.get("dataset_id") or "")
            if not ds_id:
                continue
            con.execute("""
                INSERT OR REPLACE INTO datasets
                  (dataset_id, upstream_source, accession, title, assay,
                   modality, year, n_samples, pubmed_id, license)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, [
                f"pcat:{ds_id}",
                "perturb-catalogue",
                ds_id,
                meta.get("dataset_study_title")
                  or meta.get("dataset_experiment_title"),
                ", ".join(meta.get("dataset_perturbation_types") or []),
                "perturbation",
                int(meta.get("dataset_study_year") or 0) or None,
                int(meta.get("dataset_number_of_perturbed_samples") or 0)
                  or None,
                None,
                ", ".join(meta.get("dataset_license_labels") or []),
            ])
            n += 1
    return {"perturb_catalog_datasets": n}


def ingest_multiseq(con) -> dict:
    """Pull per-cell classifications + sample assignments from the
    most recent ``Docs/MultiSeq/<run>/classifications.csv``."""
    pd = _pandas_or_none()
    if pd is None:
        return {"multiseq_cells": 0}
    base = ROOT / "Docs" / "MultiSeq"
    runs = sorted(base.glob("*/classifications.csv"),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    n = 0
    for csv in runs[:3]:
        run_id = csv.parent.name
        try:
            df = pd.read_csv(csv, index_col=0)
        except Exception:
            continue
        df = df.reset_index().rename(columns={
            df.index.name or "index": "cell_id",
        })
        if "cell_id" not in df.columns:
            df["cell_id"] = df.iloc[:, 0]
        # Stamp with the run as dataset_id so multiple runs don't collide
        df["cell_id_full"] = run_id + ":" + df["cell_id"].astype(str)
        for _, r in df.iterrows():
            con.execute("""
                INSERT OR REPLACE INTO cells
                  (cell_id, dataset_id, n_counts, multiseq_call)
                VALUES (?, ?, ?, ?)
            """, [r["cell_id_full"], f"multiseq:{run_id}",
                  int(r.get("total_tag_umi") or 0),
                  r.get("droplet_type")])
        n += len(df)
    return {"multiseq_cells": n}


def _pandas_or_none():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError:
        return None


def ingest_portal_kg(con) -> dict:
    """Pull from Data/KG/portal_kg.sqlite → entities + edges."""
    src = ROOT / "Data" / "KG" / "portal_kg.sqlite"
    if not src.exists():
        return {"portal_kg_nodes": 0, "portal_kg_edges": 0}
    sqlite_con = sqlite3.connect(f"file:{src}?mode=ro", uri=True)
    counts: "dict[str, int]" = {"portal_kg_nodes": 0, "portal_kg_edges": 0}
    # Nodes — route by node_type to the relevant entity table
    nodes = sqlite_con.execute(
        "SELECT id, node_type, label, source, properties FROM nodes "
        "LIMIT 1000000"
    ).fetchall()
    for nid, ntype, label, source, _props in nodes:
        if ntype == "gene":
            con.execute(
                "INSERT OR REPLACE INTO genes(gene_id, symbol, source_priority) "
                "VALUES (?, ?, ?)", [nid, label, source])
            counts["portal_kg_nodes"] += 1
        elif ntype == "variant":
            con.execute(
                "INSERT OR REPLACE INTO variants(variant_id, rsid) "
                "VALUES (?, ?)", [nid, label])
            counts["portal_kg_nodes"] += 1
        elif ntype == "protein":
            con.execute(
                "INSERT OR REPLACE INTO proteins(uniprot_ac, symbol) "
                "VALUES (?, ?)", [nid, label])
            counts["portal_kg_nodes"] += 1
    edges = sqlite_con.execute(
        "SELECT from_node, to_node, edge_type, source FROM edges LIMIT 1000000"
    ).fetchall()
    if edges:
        con.executemany("""
            INSERT OR IGNORE INTO edges
              (src_type, src_id, dst_type, dst_id, relation, upstream)
            VALUES (?, ?, ?, ?, ?, ?)
        """, [("node", a, "node", b, rel, src) for a, b, rel, src in edges])
        counts["portal_kg_edges"] = len(edges)
    sqlite_con.close()
    return counts


def ingest_mavedb_vampseq(con) -> dict:
    """VAMP-seq scores from MaveDB CSVs in Data/Proteomics/Sources/MaveDB/."""
    base = ROOT / "Data" / "Proteomics" / "Sources" / "MaveDB"
    if not base.exists():
        return {"vampseq_scores": 0}
    pd = _pandas_or_none()
    if pd is None:
        return {"vampseq_scores": 0}
    n = 0
    for csv in base.glob("urn-mavedb-*.csv"):
        try:
            df = pd.read_csv(csv)
        except Exception:
            continue
        if "hgvs_pro" not in df.columns or "score" not in df.columns:
            continue
        # urn → gene-name map (extend as the proteomics skill catalog grows)
        URN_TO_GENE = {"00000013-a-1": "PTEN", "00000013-b-1": "TPMT",
                         "00000078-a-1": "VKOR", "00001173-a-1": "PRKN",
                         "00000095-a-1": "CYP2C9", "00000054-a-1": "NUDT15"}
        urn_tag = csv.stem.replace("urn-mavedb-", "")
        gene = URN_TO_GENE.get(urn_tag, urn_tag)
        urn_full = "urn:mavedb:" + urn_tag.replace("-", "-")
        for i, row in df.iterrows():
            try:
                # parse "p.Met1Val" → ('M', 1, 'V')
                import re
                m = re.match(
                    r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)?$",
                    str(row.get("hgvs_pro", ""))
                )
                if not m:
                    continue
                AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
                        "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
                        "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
                        "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
                        "Ter":"*"}
                wt = AA3.get(m.group(1)); alt = AA3.get(m.group(3))
                if not wt:
                    continue
                con.execute("""
                    INSERT INTO vampseq_scores
                      (gene_id, aa_position, aa_wt, aa_alt, score, se,
                       abundance_class, dataset_id, upstream)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, [gene, int(m.group(2)), wt, alt or "?",
                      float(row["score"]) if row["score"] == row["score"]
                                            else None,
                      float(row.get("se") or 0) or None,
                      int(row.get("abundance_class") or 0) or None,
                      urn_full, "mavedb"])
                n += 1
            except Exception:
                continue
    return {"vampseq_scores": n}


# Producer registry — wire in new sources here.
PRODUCERS = {
    "proteomics-kg":    ingest_proteomics_kg,
    "portal-kg":        ingest_portal_kg,
    "perturb-catalog":  ingest_perturb_catalog,
    "multiseq":         ingest_multiseq,
    "mavedb-vampseq":   ingest_mavedb_vampseq,
}


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_init(_a) -> int:
    setup_logging("init")
    init_schema()
    print(f"Warehouse initialised: {WAREHOUSE_DB}")
    return 0


def cmd_ingest(args: argparse.Namespace) -> int:
    setup_logging("ingest_" + safe_label(args.source))
    init_schema()
    con = open_warehouse()
    sources = ([args.source]
                if args.source != "all" else list(PRODUCERS.keys()))
    summary: "dict[str, Any]" = {}
    for s in sources:
        fn = PRODUCERS.get(s)
        if fn is None:
            logger.error("Unknown source: %s. Available: %s",
                          s, ", ".join(PRODUCERS))
            continue
        logger.info("Ingesting %s …", s)
        before = _row_counts(con)
        try:
            counts = fn(con)
            success = True
        except Exception as e:
            logger.exception("ingest %s failed: %s", s, e)
            counts = {"error": str(e)}; success = False
        after = _row_counts(con)
        delta = {t: after[t] - before.get(t, 0) for t in after
                  if after[t] - before.get(t, 0) > 0}
        _record_run(con, skill="warehouse", subcommand="ingest",
                     args={"source": s},
                     rows=sum(delta.values()),
                     target_table=",".join(delta.keys()) or "-",
                     success=success)
        summary[s] = {"delta": delta, "result": counts}
    con.close()
    print(json.dumps(summary, indent=2, default=str))
    return 0


def cmd_stats(_a) -> int:
    setup_logging("stats")
    if not WAREHOUSE_DB.exists():
        print("Warehouse not initialised. Run: igvfagent warehouse init")
        return 1
    con = open_warehouse(read_only=True)
    rows = _row_counts(con)
    print("=== Row counts ===")
    for t, n in sorted(rows.items()):
        print(f"  {t:30s} {n:>10,}")
    print("\n=== Last 5 runs ===")
    for r in con.execute(
        "SELECT run_id, skill, subcommand, rows_emitted, success "
        "FROM runs ORDER BY started_at DESC LIMIT 5"
    ).fetchall():
        print(f"  {r[0]:50s} {r[3]:>8,}  ok={r[4]}")
    con.close()
    return 0


def cmd_query(args: argparse.Namespace) -> int:
    setup_logging("query")
    if not WAREHOUSE_DB.exists():
        print("Warehouse not initialised. Run: igvfagent warehouse init")
        return 1
    con = open_warehouse(read_only=True)
    rows = con.execute(args.sql).fetchall()
    cols = [d[0] for d in con.description]
    print("\t".join(cols))
    for r in rows[:args.limit]:
        print("\t".join(str(x) for x in r))
    if len(rows) > args.limit:
        print(f"... ({len(rows) - args.limit:,} more rows)")
    con.close()
    return 0


def cmd_sources(_a) -> int:
    print("Registered producers:")
    for name, fn in PRODUCERS.items():
        doc = (fn.__doc__ or "").strip().splitlines()[0] if fn.__doc__ else ""
        print(f"  {name:20s} {doc}")
    return 0


def _row_counts(con) -> "dict[str, int]":
    """Fetch row counts for every user table in the warehouse."""
    tabs = [r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='main' AND table_type='BASE TABLE'"
    ).fetchall()]
    out: "dict[str, int]" = {}
    for t in tabs:
        try:
            out[t] = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        except Exception:
            out[t] = -1
    return out


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
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
    \"\"\"One-line summary used by `sources`.\"\"\"
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
"""


def cmd_write_playbook(_a) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_TEXT)
    print(f"Wrote: {PLAYBOOK_PATH}")
    return 0


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------


def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(
        prog="warehouse",
        description="Integrated data warehouse (DuckDB) — Silver tier "
                    "for the IGVF integrated data layer.")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("init",
                        help="Create the warehouse schema.")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("ingest",
                        help="Pull a producer's outputs into the warehouse.")
    s.add_argument("--source", required=True,
                    help="Producer name (proteomics-kg | portal-kg | "
                         "perturb-catalog | multiseq | mavedb-vampseq | all).")
    s.set_defaults(func=cmd_ingest)

    s = sub.add_parser("stats",
                        help="Row counts per table + last 5 runs.")
    s.set_defaults(func=cmd_stats)

    s = sub.add_parser("query",
                        help="Run an SQL query (read-only).")
    s.add_argument("sql", help="SQL statement.")
    s.add_argument("--limit", type=int, default=50,
                    help="Max rows printed (default 50).")
    s.set_defaults(func=cmd_query)

    s = sub.add_parser("sources",
                        help="List registered producers.")
    s.set_defaults(func=cmd_sources)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/WAREHOUSE_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
