#!/usr/bin/env python3
"""Build the IGVF integrative knowledge graph from the local Catalog mirror.

The plan is Docs/Plan/datastructure.md. In short: the Catalog's 55 knowledge
collections become a core graph in LadybugDB; the billion-row measurement
collections (LD pairs, variant scores, exon structure) stay in DuckDB over
Parquet. Every snapshot is written as Parquet first, and the graph database is
rebuilt from it, never edited.

Commands:
    kg-build nodes     — node tables for the core collections (except Variant)
    kg-build edges     — edge tables, with endpoint resolution and rejects
    kg-build variants  — Variant nodes: only variants that a core edge touches
    kg-build validate  — counts, duplicate ids, reject rates -> validation.json
    kg-build graph     — load the snapshot into LadybugDB (graph.lbug)
    kg-build bench     — time a fixed set of graph queries
    kg-build all       — nodes, edges, variants, validate, graph, bench

Storage layout (one directory per snapshot):
    Data/Warehouse/KGBuild/<snapshot>/
        nodes/<Table>/<collection>.parquet
        edges/<collection>/<FromTable>__<ToTable>.parquet
        rejects/<collection>.parquet          edges that could not be placed
        manifest.json                         what went where, and why
        validation.json                       written by `validate`
        graph.lbug                            written by `graph`

Node ids are the Catalog's own `_id` (`genes/ENSG…`), so every edge in the
mirror joins without remapping. Each node also carries a `curie`.

License: Apache-2.0. Uses duckdb (MIT) and ladybug (MIT).
"""

from __future__ import annotations

import argparse
import atexit
import csv
import json
import logging
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
LOG_DIR = ROOT / "Docs" / "Logs"
WAREHOUSE_DIR = DATA_DIR / "Warehouse"
KG_DIR = Path(os.environ.get("IGVF_KG_MIRROR_DIR") or WAREHOUSE_DIR / "KG")
BUILD_DIR = Path(os.environ.get("IGVF_KG_BUILD_DIR") or WAREHOUSE_DIR / "KGBuild")
DEFAULT_SNAPSHOT = "kg-2026-10-01"

HUMAN, MOUSE = "NCBITaxon:9606", "NCBITaxon:10090"

# ---------------------------------------------------------------------------
# Schema: node tables
# ---------------------------------------------------------------------------
# Every node table has these columns first, then its typed columns, then
# `props` (the remaining source columns as one JSON string).
COMMON_NODE_COLS: list[tuple[str, str]] = [
    ("id", "STRING"), ("curie", "STRING"), ("name", "STRING"),
    ("taxon", "STRING"), ("source", "STRING"),
]

# Typed columns per table: (graph column, LadybugDB type, source column).
# `end` is reserved in Cypher, so coordinates are pos_start / pos_end.
NODE_TABLES: dict[str, list[tuple[str, str, str]]] = {
    "Gene": [("symbol", "STRING", "symbol"), ("gene_type", "STRING", "gene_type"),
             ("chr", "STRING", "chr"), ("pos_start", "INT64", "start"),
             ("pos_end", "INT64", "end"), ("strand", "STRING", "strand")],
    "Transcript": [("gene_name", "STRING", "gene_name"),
                   ("transcript_type", "STRING", "transcript_type"),
                   ("chr", "STRING", "chr"), ("pos_start", "INT64", "start"),
                   ("pos_end", "INT64", "end"), ("strand", "STRING", "strand")],
    "Protein": [("protein_id", "STRING", "protein_id"),
                ("uniprot_ids", "STRING", "uniprot_ids")],
    "RegulatoryElement": [("element_type", "STRING", "type"), ("chr", "STRING", "chr"),
                          ("pos_start", "INT64", "start"), ("pos_end", "INT64", "end"),
                          ("method", "STRING", "method")],
    "OntologyTerm": [("category", "STRING", ""), ("term_id", "STRING", "term_id"),
                     ("description", "STRING", "description")],
    "Pathway": [], "Complex": [], "Motif": [("tf_name", "STRING", "tf_name")],
    "Drug": [], "Study": [("trait_reported", "STRING", "trait_reported"),
                          ("pmid", "STRING", "pmid")],
    "Donor": [], "Dataset": [],
    "Variant": [("chr", "STRING", "chr"), ("pos", "INT64", "pos"),
                ("ref", "STRING", "ref"), ("alt", "STRING", "alt"),
                ("rsid", "STRING", "rsid"), ("spdi", "STRING", "spdi"),
                ("variation_type", "STRING", "variation_type")],
}

# Catalog node collection -> (node table, default taxon). Variant nodes are
# not listed: they are derived from the edges (see `variants`).
NODE_SOURCES: dict[str, tuple[str, str | None]] = {
    "genes": ("Gene", HUMAN), "mm_genes": ("Gene", MOUSE),
    "transcripts": ("Transcript", HUMAN), "mm_transcripts": ("Transcript", MOUSE),
    "proteins": ("Protein", None),
    "genomic_elements": ("RegulatoryElement", HUMAN),
    "mm_genomic_elements": ("RegulatoryElement", MOUSE),
    "ontology_terms": ("OntologyTerm", None),
    "pathways": ("Pathway", None), "complexes": ("Complex", None),
    "motifs": ("Motif", None), "drugs": ("Drug", None), "studies": ("Study", None),
    "donors": ("Donor", HUMAN), "files_filesets": ("Dataset", None),
}
VARIANT_SOURCES: dict[str, str] = {"variants": HUMAN, "mm_variants": MOUSE}

# Ontology prefix -> Biolink-style category.
ONTOLOGY_CATEGORY: dict[str, str] = {
    "MONDO": "Disease", "Orphanet": "Disease", "DOID": "Disease",
    "HP": "PhenotypicFeature", "EFO": "PhenotypicFeature", "OBA": "PhenotypicFeature",
    "GO": "BiologicalProcessOrActivity", "CHEBI": "ChemicalEntity",
    "CL": "Cell", "PCL": "Cell", "CVCL": "CellLine", "CLO": "CellLine",
    "UBERON": "AnatomicalEntity", "NCIT": "OntologyClass", "PR": "Protein",
}

# ---------------------------------------------------------------------------
# Schema: edge tables
# ---------------------------------------------------------------------------
# Typed edge columns, used whenever the source collection has them.
EDGE_TYPED: list[tuple[str, str]] = [
    ("source", "STRING"), ("method", "STRING"), ("label", "STRING"),
    ("biological_context", "STRING"), ("neg_log10_pvalue", "DOUBLE"),
    ("p_value", "DOUBLE"), ("effect_size", "DOUBLE"), ("log2FC", "DOUBLE"),
    ("z_score", "DOUBLE"), ("score", "DOUBLE"),
    ("posterior_inclusion_probability", "DOUBLE"), ("significant", "BOOL"),
]

# Core edge collection -> Biolink-style predicate (documentation; recorded in
# the manifest). The rel table keeps the collection's name, so every edge is
# traceable to its source collection.
EDGE_PREDICATES: dict[str, str] = {
    "genomic_elements_genes": "regulates",
    "variants_genes": "affects_expression_of",
    "variants_phenotypes": "associated_with",
    "variants_biosamples": "has_regulatory_activity_in",
    "variants_genomic_elements": "located_in_or_affects",
    "variants_biosamples_colocboost": "colocalizes_with",
    "variants_proteins_terms": "affects_binding_in_context",
    "variants_drugs": "associated_with_response_to",
    "variants_drugs_genes": "associated_with_response_to",
    "variants_diseases": "associated_with",
    "variants_diseases_genes": "associated_with",
    "proteins_proteins": "interacts_with", "genes_genes": "interacts_with",
    "mm_genes_mm_genes": "interacts_with", "genes_mm_genes": "orthologous_to",
    "genes_pathways": "participates_in", "pathways_pathways": "subclass_of",
    "gene_products_terms": "has_annotation", "genes_biosamples": "essential_in",
    "genes_transcripts": "transcribed_to", "mm_transcripts_mm_genes": "transcribed_from",
    "transcripts_proteins": "translates_to",
    "ontology_terms_ontology_terms": "subclass_of",
    "genomic_elements_biosamples": "has_regulatory_activity_in",
    "genomic_elements_phenotypes": "associated_with",
    "diseases_genes": "associated_with",
    "complexes_proteins": "has_part", "complexes_terms": "has_annotation",
    "motifs_proteins": "bound_by",
}

# Collections kept out of the graph, with the reason. They stay queryable in
# DuckDB as measurement tables (Docs/Plan/datastructure.md, milestone 3).
MEASUREMENT_ONLY: dict[str, str] = {
    "variants_variants": "LD pairs (5.9 B rows)",
    "variants_coding_variants": "variant -> coding-variant mapping",
    "coding_variants": "coding-variant records",
    "coding_variants_phenotypes": "coding-variant phenotype scores",
    "coding_variants_phenotypes_variant_painting": "endpoints are coding variants",
    "variants_proteins": "variant -> protein scores (296 M); used to resolve "
                         "variants_proteins_terms",
    "genes_structure": "exon structure", "mm_genes_structure": "exon structure",
    "transcripts_genes_structure": "exon structure",
    "mm_transcripts_mm_genes_structure": "exon structure",
    "genes_coding_variants_scores": "per-gene score blobs",
    "genes_coding_variants_scores_grp": "per-gene score blobs",
    "genes_coding_variants_phenotypes_counts": "per-gene counts",
    "genes_amino_acids_IGVF": "per-gene amino-acid tables",
    "variants_IGVF": "per-variant IGVF attributes (17 M); joined later",
    "mm_variants": "mouse variants with no core edge are not graph nodes",
    "variants": "variants with no core edge are not graph nodes",
    "files_filesets_collections": "per-collection file counts",
    "files_filesets_test": "test data",
}

# Source columns left out of a node's `props`: large per-row blobs that stay
# in the measurement tables (FAVOR annotations are ~0.5 KB per variant).
PROPS_EXCLUDE: dict[str, set[str]] = {
    "variants": {"annotations"}, "mm_variants": {"annotations"},
}

_DUCK_TYPE = {"STRING": "VARCHAR", "INT64": "BIGINT", "DOUBLE": "DOUBLE", "BOOL": "BOOLEAN"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"kg_build_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(f"Missing dependency '{name}'. {hint}\nInstall: pip install {name}") from exc


def _q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _lit(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def snap_dir(snapshot: str) -> Path:
    return BUILD_DIR / snapshot


def _duck(threads: int, memory: str, snapshot: str) -> Any:
    duckdb = _require_pkg("duckdb", "Required to build the knowledge graph.")
    tmp = snap_dir(snapshot) / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    # A file-backed staging database, so a 183 M-row edge table fits on disk
    # rather than in memory. It holds only temporary tables.
    stage = tmp / f"stage_{os.getpid()}.duckdb"
    con = duckdb.connect(str(stage))
    atexit.register(lambda: [p.unlink() for p in tmp.glob(f"{stage.name}*") if p.exists()])
    con.execute(f"SET threads={int(threads)}")
    con.execute(f"SET memory_limit={_lit(memory)}")
    con.execute(f"SET temp_directory={_lit(str(tmp))}")
    con.execute("SET preserve_insertion_order=false")
    return con


def _has_shards(collection: str) -> bool:
    d = KG_DIR / collection
    return d.is_dir() and any(d.glob("*.parquet"))


def _src(collection: str) -> str:
    return f"read_parquet({_lit(str(KG_DIR / collection / '*.parquet'))}, union_by_name=true)"


def _columns(con: Any, collection: str) -> list[str]:
    return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {_src(collection)}").fetchall()]


def _inventory() -> dict[str, dict[str, Any]]:
    path = KG_DIR / "_inventory.csv"
    if not path.exists():
        return {}
    return {r["collection"]: r for r in csv.DictReader(open(path))}


def _edge_collections() -> list[str]:
    inv = _inventory()
    cols = [c for c, r in inv.items() if r.get("type") == "edge"] if inv else list(EDGE_PREDICATES)
    return sorted(c for c in cols if c not in MEASUREMENT_ONLY and _has_shards(c))


def _manifest_path(snapshot: str) -> Path:
    return snap_dir(snapshot) / "manifest.json"


def _read_manifest(snapshot: str) -> dict[str, Any]:
    p = _manifest_path(snapshot)
    if p.exists():
        return json.loads(p.read_text())
    return {"snapshot": snapshot, "mirror": str(KG_DIR), "nodes": {}, "edges": {},
            "measurement_only": MEASUREMENT_ONLY}


def _write_manifest(snapshot: str, m: dict[str, Any]) -> None:
    p = _manifest_path(snapshot)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(m, indent=1, sort_keys=True, default=str))
    tmp.replace(p)


def _copy_to(con: Any, select_sql: str, out: Path) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".parquet.tmp")
    con.execute(f"COPY ({select_sql}) TO {_lit(str(tmp))} "
                f"(FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000)")
    tmp.replace(out)
    return con.execute(f"SELECT count(*) FROM read_parquet({_lit(str(out))})").fetchone()[0]


def _typed_expr(src_col: str, typ: str, available: set[str]) -> str:
    if src_col and src_col in available:
        return f"TRY_CAST({_q(src_col)} AS {_DUCK_TYPE[typ]})"
    return f"CAST(NULL AS {_DUCK_TYPE[typ]})"


def _props_expr(available: list[str], used: set[str]) -> str:
    rest = [c for c in available if c not in used and c not in ("_rev", "_key", "_id")]
    if not rest:
        return "CAST(NULL AS VARCHAR)"
    fields = ", ".join(f"{_lit(c)}: {_q(c)}" for c in rest)
    return f"CAST(to_json({{{fields}}}) AS VARCHAR)"


def _curie_expr(table: str) -> str:
    k = "_key"
    if table in ("Gene", "Transcript", "Protein"):
        return f"'ENSEMBL:' || {k}"
    if table == "OntologyTerm":
        return f"regexp_replace({k}, '_', ':')"
    if table == "Pathway":
        return f"'REACT:' || {k}"
    if table == "Complex":
        return f"'ComplexPortal:' || {k}"
    if table == "Drug":
        return f"'PHARMGKB.DRUG:' || {k}"
    if table == "Study":
        return f"CASE WHEN {k} LIKE 'GCST%' THEN 'GWAS:' || {k} ELSE 'igvfkg:' || _id END"
    if table == "Donor":
        return f"CASE WHEN {k} LIKE 'ENCDO%' THEN 'ENCODE:' || {k} ELSE 'igvfkg:' || _id END"
    if table == "Dataset":
        return (f"CASE WHEN {k} LIKE 'ENC%' THEN 'ENCODE:' || {k} "
                f"WHEN {k} LIKE 'IGVF%' THEN 'IGVF:' || {k} ELSE 'igvfkg:' || _id END")
    if table == "Variant":
        return f"CASE WHEN {k} LIKE 'NC\\_%' ESCAPE '\\' THEN 'SPDI:' || {k} ELSE 'igvfkg:' || _id END"
    return "'igvfkg:' || _id"


def _ontology_category_expr() -> str:
    whens = " ".join(f"WHEN {_lit(p)} THEN {_lit(c)}" for p, c in ONTOLOGY_CATEGORY.items())
    return f"CASE regexp_extract(_key, '^([A-Za-z]+)', 1) {whens} ELSE 'OntologyClass' END"


def _taxon_expr(default: str | None, available: set[str]) -> str:
    if default:
        return _lit(default)
    if "organism" in available:
        return (f"CASE WHEN organism ILIKE '%sapiens%' THEN {_lit(HUMAN)} "
                f"WHEN organism ILIKE '%musculus%' THEN {_lit(MOUSE)} "
                f"ELSE CAST(organism AS VARCHAR) END")
    return "CAST(NULL AS VARCHAR)"


def node_select(con: Any, collection: str, table: str, taxon: str | None,
                where: str = "") -> str:
    """SELECT that turns one Catalog node collection into `table`'s columns."""
    cols = _columns(con, collection)
    avail = set(cols)
    used = {"source", "name"} | PROPS_EXCLUDE.get(collection, set())
    parts = [
        "_id AS id",
        f"{_curie_expr(table)} AS curie",
        f"{_typed_expr('name', 'STRING', avail)} AS name",
        f"{_taxon_expr(taxon, avail)} AS taxon",
        f"{_typed_expr('source', 'STRING', avail)} AS source",
    ]
    for col, typ, src in NODE_TABLES[table]:
        if table == "OntologyTerm" and col == "category":
            parts.append(f"{_ontology_category_expr()} AS category")
            continue
        parts.append(f"{_typed_expr(src, typ, avail)} AS {_q(col)}")
        if src:
            used.add(src)
    if table == "Gene" and "symbol" in avail:
        parts[2] = "COALESCE(TRY_CAST(symbol AS VARCHAR), TRY_CAST(name AS VARCHAR)) AS name"
    parts.append(f"{_props_expr(cols, used)} AS props")
    return f"SELECT {', '.join(parts)} FROM {_src(collection)} {where}"


# ---------------------------------------------------------------------------
# nodes
# ---------------------------------------------------------------------------

def cmd_nodes(args: argparse.Namespace) -> int:
    setup_logging()
    con = _duck(args.threads, args.memory, args.snapshot)
    m = _read_manifest(args.snapshot)
    only = set(args.only.split(",")) if args.only else None
    for coll, (table, taxon) in NODE_SOURCES.items():
        if only and coll not in only:
            continue
        if not _has_shards(coll):
            logging.warning("%s: no shards in the mirror, skipped", coll)
            continue
        out = snap_dir(args.snapshot) / "nodes" / table / f"{coll}.parquet"
        t0 = time.time()
        n = _copy_to(con, node_select(con, coll, table, taxon) + " ORDER BY id", out)
        m["nodes"][coll] = {"table": table, "rows": n, "file": str(out.relative_to(snap_dir(args.snapshot))),
                            "seconds": round(time.time() - t0, 1)}
        logging.info("nodes %-24s -> %-18s %12d rows (%.0fs)", coll, table, n, time.time() - t0)
        _write_manifest(args.snapshot, m)
    return 0


# ---------------------------------------------------------------------------
# edges
# ---------------------------------------------------------------------------

def _prefix_table(prefix: str) -> str | None:
    if prefix in NODE_SOURCES:
        return NODE_SOURCES[prefix][0]
    if prefix in VARIANT_SOURCES:
        return "Variant"
    return None


def _node_ids_view(con: Any, snapshot: str, table: str) -> str | None:
    files = sorted((snap_dir(snapshot) / "nodes" / table).glob("*.parquet"))
    if not files:
        return None
    lst = ", ".join(_lit(str(f)) for f in files)
    return f"(SELECT id FROM read_parquet([{lst}]))"


def edge_columns(con: Any, collection: str) -> tuple[list[tuple[str, str]], list[str]]:
    cols = _columns(con, collection)
    typed = [(c, t) for c, t in EDGE_TYPED if c in set(cols)]
    return typed, cols


def cmd_edges(args: argparse.Namespace) -> int:
    setup_logging()
    con = _duck(args.threads, args.memory, args.snapshot)
    m = _read_manifest(args.snapshot)
    only = set(args.only.split(",")) if args.only else None
    edge_colls = _edge_collections()
    inv = _inventory()
    all_edges = set(edge_colls) | {c for c in MEASUREMENT_ONLY
                                   if inv.get(c, {}).get("type") == "edge" and _has_shards(c)}
    for coll in edge_colls:
        if only and coll not in only:
            continue
        t0 = time.time()
        typed, cols = edge_columns(con, coll)
        used = {"_from", "_to"} | {c for c, _ in typed}
        props = _props_expr(cols, used)
        typed_sql = ", ".join(f"TRY_CAST({_q(c)} AS {_DUCK_TYPE[t]}) AS {_q(c)}" for c, t in typed)
        # Stage the collection once with its typed columns, so each endpoint
        # pair below reads a compact table rather than the raw shards.
        con.execute(f"""CREATE OR REPLACE TABLE e AS
            SELECT _from AS raw_from, _to AS raw_to,
                   split_part(_from, '/', 1) AS fp, split_part(_to, '/', 1) AS tp,
                   {typed_sql + ',' if typed_sql else ''} {props} AS props
            FROM {_src(coll)}""")
        n_in = con.execute("SELECT count(*) FROM e").fetchone()[0]
        pairs = con.execute("SELECT fp, tp, count(*) FROM e GROUP BY ALL ORDER BY 3 DESC").fetchall()

        # An endpoint that is itself an edge (`variants_proteins/…`) is a
        # hyperedge: replace it with that edge's own _from, keep its _to in `via`.
        via_cols = ["via"]
        con.execute("ALTER TABLE e ADD COLUMN via VARCHAR")
        for fp, _tp, _n in pairs:
            if fp in all_edges and _has_shards(fp):
                con.execute(f"""UPDATE e SET raw_from = h._from, via = h._to,
                                    fp = split_part(h._from, '/', 1)
                                FROM (SELECT _id, _from, _to FROM {_src(fp)}) h
                                WHERE e.fp = {_lit(fp)} AND e.raw_from = h._id""")
        pairs = con.execute("SELECT fp, tp, count(*) FROM e GROUP BY ALL ORDER BY 3 DESC").fetchall()

        out_dir = snap_dir(args.snapshot) / "edges" / coll
        if out_dir.exists():
            shutil.rmtree(out_dir)
        rej_rows: list[str] = []
        files: dict[str, int] = {}
        typed_names = ", ".join(_q(c) for c, _ in typed)
        payload = (typed_names + ", " if typed_names else "") + ", ".join(via_cols) + ", props"
        # Group prefix pairs by the node-table pair they land in: `genes` and
        # `mm_genes` both map to Gene, and each table pair is one file.
        by_tables: dict[tuple[str, str], list[tuple[str, str]]] = {}
        for fp, tp, _n in pairs:
            ft, tt = _prefix_table(fp), _prefix_table(tp)
            if not ft or not tt:
                reason = ("hyperedge endpoint not in mirror" if fp in all_edges or tp in all_edges
                          else "unknown endpoint type")
                rej_rows.append(f"SELECT raw_from, raw_to, {_lit(reason)} AS reason "
                                f"FROM e WHERE fp = {_lit(fp)} AND tp = {_lit(tp)}")
                continue
            by_tables.setdefault((ft, tt), []).append((fp, tp))
        for (ft, tt), prefix_pairs in sorted(by_tables.items()):
            cond = "(" + " OR ".join(f"(e.fp = {_lit(a)} AND e.tp = {_lit(b)})"
                                     for a, b in prefix_pairs) + ")"
            checks = []
            for side, t in (("raw_from", ft), ("raw_to", tt)):
                if t == "Variant":
                    continue  # Variant nodes are derived from these edges
                ids = _node_ids_view(con, args.snapshot, t)
                if ids is None:
                    raise SystemExit(f"Node table {t} is missing: run `kg-build nodes` first")
                checks.append((side, ids))
            present = " AND ".join(f"e.{s} IN {ids}" for s, ids in checks) or "true"
            fname = f"{ft}__{tt}.parquet"
            files[fname] = _copy_to(
                con,
                f"SELECT raw_from AS from_id, raw_to AS to_id, {payload} FROM e "
                f"WHERE {cond} AND {present} ORDER BY from_id",
                out_dir / fname)
            if checks:
                rej_rows.append(f"SELECT raw_from, raw_to, 'dangling endpoint' AS reason "
                                f"FROM e WHERE {cond} AND NOT ({present})")
        n_rej = 0
        rej_path = snap_dir(args.snapshot) / "rejects" / f"{coll}.parquet"
        if rej_rows:
            n_rej = _copy_to(con, " UNION ALL ".join(rej_rows), rej_path)
        elif rej_path.exists():
            rej_path.unlink()
        n_out = sum(files.values())
        m["edges"][coll] = {
            "predicate": EDGE_PREDICATES.get(coll, "related_to"),
            "rows_in": n_in, "rows_out": n_out, "rejects": n_rej,
            "files": files, "typed": [c for c, _ in typed],
            "endpoint_pairs": [[a, b, c] for a, b, c in pairs],
            "seconds": round(time.time() - t0, 1),
        }
        logging.info("edges %-36s %12d in, %12d out, %10d rejected (%.0fs)",
                     coll, n_in, n_out, n_rej, time.time() - t0)
        _write_manifest(args.snapshot, m)
        con.execute("DROP TABLE e")
    return 0


# ---------------------------------------------------------------------------
# variants
# ---------------------------------------------------------------------------

def cmd_variants(args: argparse.Namespace) -> int:
    setup_logging()
    con = _duck(args.threads, args.memory, args.snapshot)
    m = _read_manifest(args.snapshot)
    files = sorted((snap_dir(args.snapshot) / "edges").glob("*/*.parquet"))
    sides = []
    for f in files:
        ft, tt = f.stem.split("__")
        if ft == "Variant":
            sides.append(f"SELECT from_id AS id FROM read_parquet({_lit(str(f))})")
        if tt == "Variant":
            sides.append(f"SELECT to_id AS id FROM read_parquet({_lit(str(f))})")
    if not sides:
        raise SystemExit("No edge touches a Variant: run `kg-build edges` first")
    t0 = time.time()
    con.execute(f"CREATE OR REPLACE TABLE want AS SELECT DISTINCT id FROM ({' UNION ALL '.join(sides)})")
    n_want = con.execute("SELECT count(*) FROM want").fetchone()[0]
    logging.info("variants referenced by core edges: %d", n_want)
    out_dir = snap_dir(args.snapshot) / "nodes" / "Variant"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    found_total = 0
    for coll, taxon in VARIANT_SOURCES.items():
        if not _has_shards(coll):
            continue
        sel = node_select(con, coll, "Variant", taxon,
                          where="WHERE _id IN (SELECT id FROM want)")
        n = _copy_to(con, sel + " ORDER BY id", out_dir / f"{coll}.parquet")
        found_total += n
        m["nodes"][coll] = {"table": "Variant", "rows": n,
                            "file": f"nodes/Variant/{coll}.parquet"}
        logging.info("variant nodes from %s: %d", coll, n)
    # Referenced variants that no mirrored shard holds (the unfinished
    # ranges): kept as id-only nodes so their edges still load.
    have = ", ".join(_lit(str(p)) for p in sorted(out_dir.glob("*.parquet")))
    cols = COMMON_NODE_COLS + [(c, t) for c, t, _ in NODE_TABLES["Variant"]] + [("props", "STRING")]
    nulls = ", ".join(
        "id" if c == "id" else
        (f"CASE WHEN id LIKE 'mm_variants/%' THEN {_lit(MOUSE)} ELSE {_lit(HUMAN)} END AS taxon"
         if c == "taxon" else
         ("'igvfkg:' || id AS curie" if c == "curie" else
          (f"'{{\"_missing\": \"not in mirror\"}}' AS props" if c == "props" else
           f"CAST(NULL AS {_DUCK_TYPE[t]}) AS {_q(c)}")))
        for c, t in cols)
    missing_sql = (f"SELECT {nulls} FROM want WHERE id NOT IN "
                   f"(SELECT id FROM read_parquet([{have}]))") if have else f"SELECT {nulls} FROM want"
    n_missing = _copy_to(con, missing_sql, out_dir / "_missing.parquet")
    m["variants"] = {"referenced": n_want, "found": found_total, "missing": n_missing,
                     "seconds": round(time.time() - t0, 1)}
    logging.info("variants: %d referenced, %d found, %d missing from the mirror",
                 n_want, found_total, n_missing)
    _write_manifest(args.snapshot, m)
    return 0


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

def cmd_validate(args: argparse.Namespace) -> int:
    setup_logging()
    con = _duck(args.threads, args.memory, args.snapshot)
    m = _read_manifest(args.snapshot)
    sd = snap_dir(args.snapshot)
    checks: list[dict[str, Any]] = []

    def check(name: str, ok: bool, detail: Any) -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        logging.info("%s %s: %s", "PASS" if ok else "FAIL", name, detail)

    inv = _inventory()
    for coll, info in sorted(m.get("nodes", {}).items()):
        if info["table"] == "Variant":
            continue
        want = int(inv[coll]["documents"]) if coll in inv else None
        check(f"nodes:{coll}:count", want is None or info["rows"] >= want,
              {"rows": info["rows"], "inventory": want})
    for table in NODE_TABLES:
        files = sorted((sd / "nodes" / table).glob("*.parquet"))
        if not files:
            continue
        lst = ", ".join(_lit(str(f)) for f in files)
        n, d = con.execute(f"SELECT count(*), count(DISTINCT id) FROM read_parquet([{lst}])").fetchone()
        check(f"nodes:{table}:unique_ids", n == d, {"rows": n, "distinct": d})
    for coll, e in sorted(m.get("edges", {}).items()):
        rate = e["rejects"] / e["rows_in"] if e["rows_in"] else 0.0
        check(f"edges:{coll}:accounted", e["rows_out"] + e["rejects"] == e["rows_in"],
              {"in": e["rows_in"], "out": e["rows_out"], "rejects": e["rejects"]})
        check(f"edges:{coll}:reject_rate", rate <= args.max_reject_rate,
              {"rate": round(rate, 6), "limit": args.max_reject_rate})
    v = m.get("variants")
    if v:
        rate = v["missing"] / v["referenced"] if v["referenced"] else 0.0
        check("variants:missing_rate", rate <= args.max_reject_rate,
              {"missing": v["missing"], "referenced": v["referenced"], "rate": round(rate, 6)})
    n_fail = sum(not c["ok"] for c in checks)
    out = {"snapshot": args.snapshot, "checked": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "passed": len(checks) - n_fail, "failed": n_fail, "checks": checks}
    (sd / "validation.json").write_text(json.dumps(out, indent=1))
    print(f"{len(checks) - n_fail}/{len(checks)} checks pass -> {sd / 'validation.json'}")
    return 0 if n_fail == 0 else 1


# ---------------------------------------------------------------------------
# graph (LadybugDB)
# ---------------------------------------------------------------------------

def _ddl_cols(cols: list[tuple[str, str]]) -> str:
    return ", ".join(f"`{c}` {t}" for c, t in cols)


def cmd_graph(args: argparse.Namespace) -> int:
    setup_logging()
    lb = _require_pkg("ladybug", "Required for the graph database (LadybugDB).")
    m = _read_manifest(args.snapshot)
    sd = snap_dir(args.snapshot)
    path = sd / "graph.lbug"
    if path.exists():
        if not args.rebuild:
            raise SystemExit(f"{path} exists; pass --rebuild to replace it")
        shutil.rmtree(path) if path.is_dir() else path.unlink()
        wal = path.with_name(path.name + ".wal")
        if wal.exists():
            wal.unlink()
    db = lb.Database(str(path), buffer_pool_size=int(args.buffer_gb * 1024 ** 3))
    con = lb.Connection(db, num_threads=args.threads)
    timings: dict[str, float] = {}
    tables_loaded = set()
    for table, typed in NODE_TABLES.items():
        files = sorted((sd / "nodes" / table).glob("*.parquet"))
        if not files:
            continue
        cols = COMMON_NODE_COLS + [(c, t) for c, t, _ in typed] + [("props", "STRING")]
        con.execute(f"CREATE NODE TABLE {table}({_ddl_cols(cols)}, PRIMARY KEY(id))")
        t0 = time.time()
        lst = ", ".join(_lit(str(f)) for f in files)
        con.execute(f"COPY {table} FROM [{lst}]")
        timings[f"node:{table}"] = round(time.time() - t0, 1)
        tables_loaded.add(table)
        logging.info("graph: loaded %s from %d files (%.0fs)", table, len(files), time.time() - t0)
    for coll, e in sorted(m.get("edges", {}).items()):
        files = sorted((sd / "edges" / coll).glob("*.parquet"))
        pairs = [f.stem.split("__") for f in files]
        pairs = [(a, b) for a, b in pairs if a in tables_loaded and b in tables_loaded]
        if not pairs:
            continue
        typed = [(c, t) for c, t in EDGE_TYPED if c in set(e["typed"])]
        cols = typed + [("via", "STRING"), ("props", "STRING")]
        froms = ", ".join(f"FROM {a} TO {b}" for a, b in pairs)
        con.execute(f"CREATE REL TABLE {coll}({froms}, {_ddl_cols(cols)})")
        t0 = time.time()
        for a, b in pairs:
            f = sd / "edges" / coll / f"{a}__{b}.parquet"
            con.execute(f"COPY {coll} FROM {_lit(str(f))} (from={_lit(a)}, to={_lit(b)})")
        timings[f"rel:{coll}"] = round(time.time() - t0, 1)
        logging.info("graph: loaded %s, %d endpoint pairs (%.0fs)", coll, len(pairs), time.time() - t0)
    size = sum(p.stat().st_size for p in ([path] if path.is_file() else path.rglob("*")) if p.is_file())
    m["graph"] = {"path": str(path.relative_to(sd)), "bytes": size, "load_seconds": timings,
                  "built": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _write_manifest(args.snapshot, m)
    print(f"graph.lbug: {size / 1e9:.1f} GB, load {sum(timings.values()):.0f}s")
    return 0


# ---------------------------------------------------------------------------
# bench
# ---------------------------------------------------------------------------

BENCH_QUERIES: list[tuple[str, str]] = [
    ("gene_by_symbol", "MATCH (g:Gene) WHERE g.symbol = 'GATA1' RETURN g.id, g.chr, g.pos_start"),
    ("gene_1hop_elements",
     "MATCH (e:RegulatoryElement)-[r:genomic_elements_genes]->(g:Gene) "
     "WHERE g.symbol = 'MYC' RETURN count(*)"),
    ("variant_gene_pathway_3hop",
     "MATCH (v:Variant)-[:variants_genes]->(g:Gene)-[:genes_pathways]->(p:Pathway) "
     "WHERE g.symbol = 'PCSK9' RETURN count(DISTINCT v), count(DISTINCT p)"),
    ("disease_genes_proteins",
     "MATCH (d:OntologyTerm)-[:diseases_genes]->(g:Gene)-[:genes_transcripts]->(:Transcript)"
     "-[:transcripts_proteins]->(p:Protein) WHERE d.name CONTAINS 'cardiomyopathy' "
     "RETURN count(DISTINCT g), count(DISTINCT p)"),
    ("ppi_2hop",
     "MATCH (a:Protein)-[:proteins_proteins]-(b:Protein)-[:proteins_proteins]-(c:Protein) "
     "WHERE a.protein_id = 'ENSP00000269305' RETURN count(DISTINCT c)"),
]


def cmd_bench(args: argparse.Namespace) -> int:
    setup_logging()
    lb = _require_pkg("ladybug", "Required for the graph database (LadybugDB).")
    sd = snap_dir(args.snapshot)
    db = lb.Database(str(sd / "graph.lbug"), read_only=True,
                     buffer_pool_size=int(args.buffer_gb * 1024 ** 3))
    con = lb.Connection(db, num_threads=args.threads)
    results = []
    for name, q in BENCH_QUERIES:
        runs, rows, err = [], None, None
        for _ in range(args.repeat):
            t0 = time.time()
            try:
                rows = con.execute(q).get_all()
            except Exception as exc:  # report the failure; keep timing the rest
                err = str(exc)[:300]
                break
            runs.append(time.time() - t0)
        results.append({"query": name, "cypher": q, "seconds": [round(r, 4) for r in runs],
                        "result": rows[:5] if rows else rows, "error": err})
        logging.info("bench %-28s %s", name, err or f"{min(runs):.3f}s (first {runs[0]:.3f}s) -> {rows[:3]}")
    m = _read_manifest(args.snapshot)
    m["bench"] = {"run": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}
    _write_manifest(args.snapshot, m)
    return 0 if all(r["error"] is None for r in results) else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_all(args: argparse.Namespace) -> int:
    for step in (cmd_nodes, cmd_edges, cmd_variants):
        rc = step(args)
        if rc:
            return rc
    rc_valid = cmd_validate(args)  # failures are recorded; the build continues
    rc = cmd_graph(args)
    if rc:
        return rc
    return cmd_bench(args) or rc_valid


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the IGVF knowledge graph (Parquet + LadybugDB).")
    parser.add_argument("--snapshot", default=DEFAULT_SNAPSHOT)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    parser.add_argument("--memory", default="48GB", help="DuckDB memory limit")
    sub = parser.add_subparsers(dest="command", required=True)

    for name, help_ in (("nodes", "Node tables for the core collections."),
                        ("edges", "Edge tables, endpoint resolution and rejects.")):
        p = sub.add_parser(name, help=help_)
        p.add_argument("--only", default=None, help="Comma-separated collections.")
    sub.add_parser("variants", help="Variant nodes touched by a core edge.")
    p = sub.add_parser("validate", help="Counts, unique ids, reject rates.")
    p.add_argument("--max-reject-rate", type=float, default=0.01)
    p = sub.add_parser("graph", help="Load the snapshot into LadybugDB.")
    p.add_argument("--rebuild", action="store_true")
    p.add_argument("--buffer-gb", type=float, default=40)
    p = sub.add_parser("bench", help="Time a fixed set of graph queries.")
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--buffer-gb", type=float, default=40)
    p = sub.add_parser("all", help="nodes, edges, variants, validate, graph, bench.")
    p.add_argument("--only", default=None)
    p.add_argument("--max-reject-rate", type=float, default=0.01)
    p.add_argument("--rebuild", action="store_true")
    p.add_argument("--buffer-gb", type=float, default=40)
    p.add_argument("--repeat", type=int, default=3)

    args = parser.parse_args(argv)
    return {
        "nodes": cmd_nodes, "edges": cmd_edges, "variants": cmd_variants,
        "validate": cmd_validate, "graph": cmd_graph, "bench": cmd_bench, "all": cmd_all,
    }[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
