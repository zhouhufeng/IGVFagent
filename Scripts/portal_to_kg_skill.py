#!/usr/bin/env python3
"""IGVF Portal → local Knowledge Graph ETL skill.

Pulls unstructured IGVF Portal entities (AnalysisSets, MeasurementSets,
Samples, Donors, Files, FileSets) by tissue / gene / assay / lab filters,
annotates them against genes and variants, and persists the result to a
local SQLite-backed knowledge graph that mirrors the IGVF Catalog
ArangoDB schema (nodes + edges + properties + provenance).

The local KG is read/write locally and can later be exported as
arangoimport-compatible JSONL for push to the IGVF Catalog. All endpoint
URLs and credentials live exclusively in environment variables resolved
through ``_endpoints.py``; nothing sensitive is written to source.

Subcommands

  pull               Fetch Portal entities by filter and upsert into the
                     local KG, expanding linked samples / donors / files
                     up to ``--depth`` hops.
  annotate           Text-mine description / summary / alias fields of
                     existing Portal nodes for gene-symbol and variant
                     (rsID / SPDI) mentions, confirm against the IGVF
                     Catalog, and create ``mentions_gene`` /
                     ``mentions_variant`` edges.
  enrich             For genes in the local KG, hydrate Catalog evidence
                     (variants, regulatory elements, transcripts,
                     proteins, diseases, pathways) and add as nodes and
                     edges. Caches locally so it does not re-fetch.
  query              Convenience read interface over the local KG (by
                     gene, by tissue, by node id, or by free-text label).
  stats              Per-type node and edge counts plus ingestion summary.
  export-aql         Dump each node-type collection plus edges as
                     arangoimport-compatible JSONL under
                     ``Data/KG/Export/<timestamp>/``.
  export-cytoscape   Dump a cytoscape.js-compatible JSON for visualization.
  write-playbook     Emit ``Docs/Skills/PORTAL_TO_KG_SKILLS.md``.

Local store: ``Data/KG/local_kg.sqlite`` (gitignored under ``Data/*``).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
KG_DIR = DATA_DIR / "KG"
KG_PATH = KG_DIR / "local_kg.sqlite"
KG_CACHE_DIR = DATA_DIR / "Cache" / "KGLocal"
KG_EXPORT_DIR = KG_DIR / "Export"
REPORT_DIR = DOCS_DIR / "PortalToKG"

PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")

USER_AGENT = "IGVFdataAgent-PortalToKG/0.1"

# Heuristic gene-symbol token regex. We later confirm against Catalog so
# false positives (gene-like uppercase words: "PCR", "DNA", "RNA") are
# filtered out before being persisted as Gene nodes.
GENE_TOKEN_RE = re.compile(r"\b[A-Z][A-Z0-9-]{1,9}\b")
RSID_RE = re.compile(r"\brs\d{2,}\b", re.I)
SPDI_RE = re.compile(r"\bNC_\d+\.\d+:\d+:[ACGTN-]+:[ACGTN-]+\b")
HGVS_RE = re.compile(r"\bNC_\d+\.\d+:g\.\d+[ACGT>+-]+\b")

# Common false-positive tokens that match the gene regex but are not genes.
GENE_TOKEN_DENYLIST = {
    "DNA", "RNA", "PCR", "FACS", "WT", "KO", "KI", "ORF", "ATAC", "RNASEQ",
    "HG38", "GRCH38", "MM10", "MM39", "GRCM39", "GENCODE", "CRISPR", "CRISPRI",
    "MPRA", "STARR", "CCRE", "ENCODE", "IGVF", "ABC", "RE2G", "VCF", "BAM",
    "TSV", "CSV", "JSON", "URL", "API", "IGVFDS", "IGVFFI", "IGVFDO", "IGVFSM",
    "ID", "QC", "TSS", "BP", "KB", "MB", "GB", "EFO", "CL", "MONDO", "DOI",
    "HGNC", "ENSG", "ENST", "ENSP", "FFPE", "GEM", "WELL", "POOL", "UMI",
    "USA", "EU", "CRISPR-FACS", "CHIP", "CHIP-SEQ", "STAR-SEQ", "DAY", "CC",
}

# IGVF Portal entity types this skill knows how to ingest.
PORTAL_TYPES = (
    "AnalysisSet", "MeasurementSet", "Sample", "Donor", "File",
    "FileSet", "Document",
)


# ----------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"portal_to_kg_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (KG_DIR, KG_CACHE_DIR, KG_EXPORT_DIR, REPORT_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z") or time.strftime("%Y-%m-%dT%H:%M:%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


# ----------------------------- HTTP helpers ----------------------------------

def request_headers() -> dict[str, str]:
    h = {"User-Agent": USER_AGENT, "Accept": "application/json,*/*"}
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        h["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return h


def http_get(url: str, params: dict | None = None,
             timeout: int = 60) -> tuple[int, Any]:
    if params:
        url = (url + ("&" if "?" in url else "?")
               + urllib.parse.urlencode({k: v for k, v in params.items()
                                            if v is not None}, doseq=True))
    logging.info("GET %s", url)
    req = urllib.request.Request(url, headers=request_headers(), method="GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            content = resp.read()
            try:
                return resp.status, json.loads(content)
            except json.JSONDecodeError:
                return resp.status, content.decode(errors="replace")
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace") if e.fp else ""
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            data = {"http_error_body": body}
        return e.code, data
    except urllib.error.URLError as e:
        return 0, {"network_error": str(e.reason)}


def portal_get(path: str, **params) -> tuple[int, Any]:
    return http_get(PORTAL_API_BASE + path, params=params)


def catalog_get(path: str, **params) -> tuple[int, Any]:
    return http_get(CATALOG_API_BASE + path, params=params)


# ----------------------------- SQLite KG layer -------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS nodes (
    id           TEXT PRIMARY KEY,
    node_type    TEXT NOT NULL,
    label        TEXT,
    source       TEXT NOT NULL,
    source_url   TEXT,
    properties   TEXT,
    ingested_at  TEXT NOT NULL,
    updated_at   TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_nodes_type   ON nodes(node_type);
CREATE INDEX IF NOT EXISTS idx_nodes_label  ON nodes(label);
CREATE INDEX IF NOT EXISTS idx_nodes_source ON nodes(source);

CREATE TABLE IF NOT EXISTS edges (
    id           TEXT PRIMARY KEY,
    from_node    TEXT NOT NULL,
    to_node      TEXT NOT NULL,
    edge_type    TEXT NOT NULL,
    source       TEXT NOT NULL,
    properties   TEXT,
    ingested_at  TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_edges_from ON edges(from_node);
CREATE INDEX IF NOT EXISTS idx_edges_to   ON edges(to_node);
CREATE INDEX IF NOT EXISTS idx_edges_type ON edges(edge_type);

CREATE TABLE IF NOT EXISTS gene_cache (
    symbol       TEXT PRIMARY KEY,
    confirmed    INTEGER NOT NULL,
    payload      TEXT,
    cached_at    TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS run_log (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    cmd          TEXT NOT NULL,
    args         TEXT NOT NULL,
    started_at   TEXT NOT NULL,
    ended_at     TEXT,
    nodes_added  INTEGER DEFAULT 0,
    edges_added  INTEGER DEFAULT 0,
    notes        TEXT
);
"""


def kg_connect() -> sqlite3.Connection:
    KG_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(KG_PATH))
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


def upsert_node(conn: sqlite3.Connection, node_id: str, node_type: str,
                  label: str | None, source: str, source_url: str | None,
                  properties: dict | None) -> bool:
    """Returns True if newly inserted."""
    ts = now()
    cur = conn.execute("SELECT id FROM nodes WHERE id = ?", (node_id,))
    existed = cur.fetchone() is not None
    payload = json.dumps(properties or {}, default=str)
    if existed:
        conn.execute(
            "UPDATE nodes SET node_type=?, label=?, source=?, source_url=?, "
            "properties=?, updated_at=? WHERE id=?",
            (node_type, label, source, source_url, payload, ts, node_id),
        )
        return False
    conn.execute(
        "INSERT INTO nodes(id, node_type, label, source, source_url, "
        "properties, ingested_at, updated_at) VALUES (?,?,?,?,?,?,?,?)",
        (node_id, node_type, label, source, source_url, payload, ts, ts),
    )
    return True


def upsert_edge(conn: sqlite3.Connection, edge_id: str, from_node: str,
                  to_node: str, edge_type: str, source: str,
                  properties: dict | None) -> bool:
    ts = now()
    cur = conn.execute("SELECT id FROM edges WHERE id = ?", (edge_id,))
    if cur.fetchone() is not None:
        conn.execute(
            "UPDATE edges SET from_node=?, to_node=?, edge_type=?, source=?, "
            "properties=? WHERE id=?",
            (from_node, to_node, edge_type, source,
             json.dumps(properties or {}, default=str), edge_id),
        )
        return False
    conn.execute(
        "INSERT INTO edges(id, from_node, to_node, edge_type, source, "
        "properties, ingested_at) VALUES (?,?,?,?,?,?,?)",
        (edge_id, from_node, to_node, edge_type, source,
         json.dumps(properties or {}, default=str), ts),
    )
    return True


def edge_id_for(from_id: str, to_id: str, etype: str) -> str:
    return f"{etype}::{from_id}::{to_id}"


def log_run(conn, cmd: str, args: dict, started_at: str,
             nodes_added: int = 0, edges_added: int = 0,
             notes: str = "") -> None:
    conn.execute(
        "INSERT INTO run_log(cmd, args, started_at, ended_at, nodes_added, "
        "edges_added, notes) VALUES (?,?,?,?,?,?,?)",
        (cmd, json.dumps(args, default=str), started_at, now(),
         nodes_added, edges_added, notes),
    )


# --------------------------- Portal entity ingest ----------------------------

def portal_node_id(entity: dict) -> str | None:
    """Stable ID for a Portal entity, mirroring the Catalog collection layout."""
    coll = portal_collection_for(entity)
    acc = entity.get("accession") or entity.get("uuid") or entity.get("@id")
    if not acc:
        return None
    if isinstance(acc, str) and acc.startswith("/"):
        acc = acc.strip("/").split("/")[-1]
    return f"{coll}/{acc}"


def portal_collection_for(entity: dict) -> str:
    types = entity.get("@type") or []
    if isinstance(types, str):
        types = [types]
    primary = (types[0] if types else "Item").lower()
    return {
        "analysisset": "analysis_sets",
        "measurementset": "measurement_sets",
        "fileset": "file_sets",
        "sample": "samples",
        "donor": "donors",
        "file": "files",
        "matrixfile": "files",
        "tabularfile": "files",
        "alignmentfile": "files",
        "configurationfile": "files",
        "sequencefile": "files",
        "indexfile": "files",
        "document": "documents",
    }.get(primary, "items")


def upsert_portal_entity(conn, entity: dict, source_url: str = "") -> str | None:
    nid = portal_node_id(entity)
    if not nid:
        return None
    ntype = (entity.get("@type") or ["Item"])[0]
    label = (entity.get("accession")
             or entity.get("name")
             or entity.get("title")
             or "")
    upsert_node(conn, nid, ntype, label, "igvf_portal", source_url, entity)
    return nid


def expand_portal_links(conn, entity: dict, parent_id: str,
                          depth: int, visited: set,
                          fetch_remote: bool) -> tuple[int, int]:
    """Follow common Portal cross-references and ingest them."""
    nodes_added = edges_added = 0
    if depth <= 0:
        return 0, 0
    fields_to_follow = [
        ("samples",          "profiles_sample"),
        ("donors",           "from_donor"),
        ("files",            "produces_file"),
        ("input_file_sets",  "derived_from"),
        ("workflows",        "uses_workflow"),
        ("documents",        "documented_by"),
    ]
    for field, edge_label in fields_to_follow:
        children = entity.get(field) or []
        if isinstance(children, dict):
            children = [children]
        for child in children:
            if isinstance(child, str):
                # @id reference like "/samples/IGVFSM…/"
                if not fetch_remote:
                    cid = path_to_node_id(child)
                    if cid and cid not in visited:
                        # Make a stub node — properties only @id
                        upsert_node(conn, cid, infer_type_from_path(child),
                                     cid.split("/")[-1], "igvf_portal_stub",
                                     PORTAL_API_BASE + child,
                                     {"@id": child})
                        nodes_added += 1
                        visited.add(cid)
                    if cid:
                        if upsert_edge(conn, edge_id_for(parent_id, cid, edge_label),
                                        parent_id, cid, edge_label,
                                        "igvf_portal", {}):
                            edges_added += 1
                    continue
                # Hydrate the @id reference
                status, data = http_get(PORTAL_API_BASE + child)
                if status != 200 or not isinstance(data, dict):
                    continue
                cid = upsert_portal_entity(conn, data, PORTAL_API_BASE + child)
                if cid and cid not in visited:
                    nodes_added += 1; visited.add(cid)
                    sub_n, sub_e = expand_portal_links(conn, data, cid,
                                                         depth - 1, visited,
                                                         fetch_remote)
                    nodes_added += sub_n; edges_added += sub_e
                if cid and upsert_edge(conn,
                                          edge_id_for(parent_id, cid, edge_label),
                                          parent_id, cid, edge_label,
                                          "igvf_portal", {}):
                    edges_added += 1
            elif isinstance(child, dict):
                cid = upsert_portal_entity(conn, child)
                if cid:
                    if cid not in visited:
                        nodes_added += 1; visited.add(cid)
                    if upsert_edge(conn,
                                      edge_id_for(parent_id, cid, edge_label),
                                      parent_id, cid, edge_label,
                                      "igvf_portal", {}):
                        edges_added += 1
    return nodes_added, edges_added


def path_to_node_id(at_id: str) -> str | None:
    parts = (at_id or "").strip("/").split("/")
    if len(parts) < 2:
        return None
    coll = parts[0].replace("-", "_")
    return f"{coll}/{parts[1]}"


def infer_type_from_path(at_id: str) -> str:
    parts = (at_id or "").strip("/").split("/")
    if not parts:
        return "Item"
    head = parts[0].rstrip("s")
    return "".join(w.capitalize() for w in head.split("-")) or "Item"


# ----------------------------- pull subcommand --------------------------------

def cmd_pull(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    started = now()
    conn = kg_connect()
    visited: set[str] = set()
    n_added = e_added = 0

    types = args.types or [args.type] if args.type else ["AnalysisSet"]
    queries = []
    for t in types:
        params = {"type": t, "format": "json", "limit": args.limit}
        if args.tissue:
            params["searchTerm"] = args.tissue
        if args.gene and not args.tissue:
            params["searchTerm"] = args.gene
        elif args.gene:
            params["searchTerm"] = f"{args.tissue} {args.gene}"
        if args.assay:
            params["preferred_assay_titles"] = args.assay
        if args.lab:
            params["lab.title"] = args.lab
        if args.status:
            params["status"] = args.status
        queries.append(params)

    seen_ids: set[str] = set()
    for params in queries:
        status, data = portal_get("/search/", **params)
        if status != 200 or not isinstance(data, dict):
            logging.warning("Portal /search/ → %s", status); continue
        graph = data.get("@graph") or []
        logging.info("query %s → %d hits (total=%s)",
                      params, len(graph), data.get("total"))
        for g in graph:
            if not args.skip_hydrate:
                at_id = g.get("@id")
                if at_id:
                    s2, hydrated = http_get(PORTAL_API_BASE + at_id)
                    if s2 == 200 and isinstance(hydrated, dict):
                        g = hydrated
            nid = upsert_portal_entity(conn, g,
                                          source_url=PORTAL_API_BASE + (g.get("@id") or ""))
            if not nid:
                continue
            if nid not in seen_ids:
                seen_ids.add(nid); n_added += 1
            sub_n, sub_e = expand_portal_links(
                conn, g, nid, depth=args.depth, visited=visited,
                fetch_remote=not args.skip_hydrate,
            )
            n_added += sub_n; e_added += sub_e
            conn.commit()

    log_run(conn, "pull", vars(args), started, n_added, e_added)
    conn.commit(); conn.close()

    print(f"Local KG: {KG_PATH}")
    print(f"Nodes added: {n_added}, edges added: {e_added}")
    return KG_PATH


# ----------------------------- annotate subcommand ---------------------------

def is_known_gene(symbol: str, conn) -> bool | None:
    """Returns True/False if cached, None if unknown."""
    cur = conn.execute("SELECT confirmed FROM gene_cache WHERE symbol = ?",
                        (symbol,))
    row = cur.fetchone()
    if row is None:
        return None
    return bool(row[0])


def confirm_gene_via_catalog(symbol: str, conn) -> dict | None:
    cached = is_known_gene(symbol, conn)
    if cached is False:
        return None
    if cached is True:
        cur = conn.execute("SELECT payload FROM gene_cache WHERE symbol = ?",
                            (symbol,))
        row = cur.fetchone()
        return json.loads(row[0]) if row and row[0] else None

    status, data = catalog_get("/api/genes", name=symbol, limit=1)
    rec = None
    if status == 200:
        rows = data if isinstance(data, list) else (data.get("results") or
                                                     data.get("@graph") or [])
        if isinstance(rows, list) and rows:
            rec = rows[0] if isinstance(rows[0], dict) else None
    confirmed = 1 if rec else 0
    conn.execute(
        "INSERT OR REPLACE INTO gene_cache(symbol, confirmed, payload, "
        "cached_at) VALUES (?,?,?,?)",
        (symbol, confirmed, json.dumps(rec or {}, default=str), now()),
    )
    conn.commit()
    return rec


def annotate_node(conn, row: tuple, gene_terms: set[str] | None,
                    confirm: bool) -> tuple[int, int]:
    nid, node_type, label, props_json = row
    if not props_json:
        return 0, 0
    props = json.loads(props_json or "{}")
    text = " ".join(str(props.get(k, "")) for k in
                     ("description", "summary", "title", "name", "aliases",
                       "label"))
    if not text.strip():
        return 0, 0

    n_added = e_added = 0
    # Variant mentions
    for vid_match in (RSID_RE.findall(text) + SPDI_RE.findall(text)
                       + HGVS_RE.findall(text)):
        v_norm = vid_match
        v_node = f"variants/{v_norm}"
        if upsert_node(conn, v_node, "Variant", v_norm, "annotation",
                        None, {"id": v_norm, "raw_match": vid_match}):
            n_added += 1
        if upsert_edge(conn, edge_id_for(nid, v_node, "mentions_variant"),
                        nid, v_node, "mentions_variant", "annotation",
                        {"context": text[:140]}):
            e_added += 1

    # Gene mentions
    candidates = set(GENE_TOKEN_RE.findall(text)) - GENE_TOKEN_DENYLIST
    if gene_terms is not None:
        candidates &= gene_terms
    for sym in sorted(candidates):
        if confirm:
            rec = confirm_gene_via_catalog(sym, conn)
            if not rec:
                continue
            ensg = rec.get("gene_id") or rec.get("ensembl_id") or sym
            g_node = f"genes/{ensg}"
            label_g = rec.get("name") or sym
            if upsert_node(conn, g_node, "Gene", label_g, "igvf_catalog",
                            CATALOG_API_BASE + f"/api/genes?name={sym}", rec):
                n_added += 1
        else:
            g_node = f"genes/{sym}"
            if upsert_node(conn, g_node, "Gene", sym, "annotation_unconfirmed",
                            None, {"symbol": sym}):
                n_added += 1
        if upsert_edge(conn, edge_id_for(nid, g_node, "mentions_gene"),
                        nid, g_node, "mentions_gene", "annotation",
                        {"matched_token": sym, "context": text[:140]}):
            e_added += 1
    return n_added, e_added


def cmd_annotate(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    started = now()
    conn = kg_connect()
    gene_terms: set[str] | None = None
    if args.gene_list:
        gene_terms = {ln.strip() for ln in Path(args.gene_list).read_text().splitlines()
                      if ln.strip()}
    where = "node_type IN ('AnalysisSet','MeasurementSet','FileSet','Sample','Donor','File','Document')"
    if args.where:
        where += f" AND ({args.where})"
    cur = conn.execute(
        f"SELECT id, node_type, label, properties FROM nodes WHERE {where}"
    )
    n_added = e_added = 0
    for row in cur.fetchall():
        sub_n, sub_e = annotate_node(conn, row, gene_terms,
                                       confirm=not args.skip_catalog_confirm)
        n_added += sub_n; e_added += sub_e
        conn.commit()
    log_run(conn, "annotate", vars(args), started, n_added, e_added)
    conn.commit(); conn.close()
    print(f"Annotation added {n_added} nodes, {e_added} edges in {KG_PATH}.")
    return KG_PATH


# ----------------------------- enrich subcommand -----------------------------

GENE_RELATIONS = {
    "variants":              "/api/genes/variants",
    "transcripts":           "/api/genes/transcripts",
    "proteins":              "/api/genes/proteins",
    "regulatory_elements":   "/api/genes/genomic-elements",
    "diseases":              "/api/genes/diseases",
    "pathways":              "/api/genes/pathways",
    "coding_variant_scores": "/api/genes/coding-variants/scores",
}


def hydrate_gene_from_catalog(conn, gene_node_id: str, gene_symbol: str,
                                limit: int) -> tuple[int, int]:
    n_added = e_added = 0
    for rel, path in GENE_RELATIONS.items():
        status, data = catalog_get(path, gene_name=gene_symbol, limit=limit,
                                     verbose="false")
        if status != 200:
            continue
        rows = data if isinstance(data, list) else (data.get("results")
                                                     or data.get("@graph") or [])
        if not isinstance(rows, list):
            continue
        for r in rows:
            if not isinstance(r, dict):
                continue
            target = pick_target_collection(rel, r)
            if not target:
                continue
            t_id, t_type, t_label, t_props = target
            if upsert_node(conn, t_id, t_type, t_label, "igvf_catalog",
                            CATALOG_API_BASE + path, t_props):
                n_added += 1
            etype = catalog_edge_type_for(rel)
            if upsert_edge(conn, edge_id_for(gene_node_id, t_id, etype),
                            gene_node_id, t_id, etype, "igvf_catalog", r):
                e_added += 1
    return n_added, e_added


def pick_target_collection(rel: str, row: dict
                            ) -> tuple[str, str, str, dict] | None:
    if rel == "variants":
        sv = row.get("sequence_variant") or ""
        if isinstance(sv, str) and sv.startswith("variants/"):
            v = sv.split("/", 1)[1]
            return f"variants/{v}", "Variant", v, {"id": v, **row}
    if rel == "transcripts":
        tr = row.get("transcript") or row.get("transcript_id") or ""
        if isinstance(tr, str) and tr:
            tid = tr.split("/")[-1]
            return f"transcripts/{tid}", "Transcript", tid, row
    if rel == "proteins":
        pp = row.get("protein") or ""
        if isinstance(pp, str) and pp.startswith("proteins/"):
            pid = pp.split("/", 1)[1]
            return f"proteins/{pid}", "Protein", pid, row
    if rel == "regulatory_elements":
        ge = row.get("genomic_element") or ""
        if isinstance(ge, str) and ge.startswith("genomic_elements/"):
            geid = ge.split("/", 1)[1]
            return f"genomic_elements/{geid}", "GenomicElement", geid, row
    if rel == "diseases":
        d = row.get("disease") or row.get("term") or row.get("ontology_term")
        if isinstance(d, str):
            did = d.split("/")[-1]
            return f"diseases/{did}", "Disease", did, row
    if rel == "pathways":
        p = row.get("pathway") or row.get("term")
        if isinstance(p, str):
            pid = p.split("/")[-1]
            return f"pathways/{pid}", "Pathway", pid, row
    if rel == "coding_variant_scores":
        cv = (row.get("protein_change") or {}).get("coding_variant_id", "")
        if cv:
            return f"coding_variants/{cv}", "CodingVariant", cv, row
    return None


def catalog_edge_type_for(rel: str) -> str:
    return {
        "variants":              "regulates_or_is_modulated_by",
        "transcripts":           "transcribed_to",
        "proteins":              "translated_to",
        "regulatory_elements":   "regulated_by",
        "diseases":              "associated_with_disease",
        "pathways":              "in_pathway",
        "coding_variant_scores": "has_coding_variant",
    }.get(rel, rel)


def cmd_enrich(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    started = now()
    conn = kg_connect()
    where = "node_type='Gene'"
    if args.symbols:
        sym_set = {s.strip() for s in args.symbols.split(",") if s.strip()}
        in_clause = ",".join("?" * len(sym_set))
        cur = conn.execute(
            f"SELECT id, label FROM nodes WHERE {where} AND label IN ({in_clause})",
            tuple(sym_set),
        )
    else:
        cur = conn.execute(f"SELECT id, label FROM nodes WHERE {where}")
    rows = cur.fetchall()
    logging.info("Enriching %d genes from Catalog", len(rows))
    n_added = e_added = 0
    for nid, sym in rows:
        if not sym:
            continue
        sub_n, sub_e = hydrate_gene_from_catalog(conn, nid, sym, args.limit)
        n_added += sub_n; e_added += sub_e
        conn.commit()
    log_run(conn, "enrich", vars(args), started, n_added, e_added)
    conn.commit(); conn.close()
    print(f"Catalog enrichment added {n_added} nodes, {e_added} edges.")
    return KG_PATH


# ----------------------------- query subcommand ------------------------------

def cmd_query(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    conn = kg_connect()
    out_lines: list[str] = []

    if args.gene:
        sym = args.gene
        # Find gene node by label or canonical id
        cur = conn.execute(
            "SELECT id, label, properties FROM nodes "
            "WHERE node_type='Gene' AND (label=? OR id=?)",
            (sym, f"genes/{sym}"),
        )
        gene_rows = cur.fetchall()
        out_lines.append(f"# Local KG query: gene `{sym}`")
        out_lines.append("")
        if not gene_rows:
            out_lines.append("_no Gene node found locally; "
                              "run `pull` + `annotate` first._")
        for gid, glabel, _ in gene_rows:
            out_lines.append(f"## Node `{gid}`  (label: {glabel})")
            cur2 = conn.execute(
                "SELECT e.edge_type, n.id, n.node_type, n.label "
                "FROM edges e JOIN nodes n ON n.id = e.from_node "
                "WHERE e.to_node = ?", (gid,))
            inbound = cur2.fetchall()
            out_lines.append(f"### Inbound edges ({len(inbound)})")
            for et, src_id, src_type, src_label in inbound[: args.limit]:
                out_lines.append(f"- ({et}) ← `{src_id}` [{src_type}] {src_label}")
            cur2 = conn.execute(
                "SELECT e.edge_type, n.id, n.node_type, n.label "
                "FROM edges e JOIN nodes n ON n.id = e.to_node "
                "WHERE e.from_node = ?", (gid,))
            outbound = cur2.fetchall()
            out_lines.append(f"### Outbound edges ({len(outbound)})")
            for et, tgt_id, tgt_type, tgt_label in outbound[: args.limit]:
                out_lines.append(f"- ({et}) → `{tgt_id}` [{tgt_type}] {tgt_label}")
    elif args.tissue:
        cur = conn.execute(
            "SELECT id, node_type, label, source FROM nodes "
            "WHERE properties LIKE ? OR label LIKE ? LIMIT ?",
            (f"%{args.tissue}%", f"%{args.tissue}%", args.limit),
        )
        out_lines.append(f"# Local KG query: tissue `{args.tissue}`\n")
        for nid, ntype, label, src in cur.fetchall():
            out_lines.append(f"- `{nid}` [{ntype}] {label}  _(source={src})_")
    elif args.node_id:
        cur = conn.execute(
            "SELECT node_type, label, source, source_url, properties "
            "FROM nodes WHERE id = ?", (args.node_id,))
        row = cur.fetchone()
        if not row:
            out_lines.append("_node not found_")
        else:
            ntype, label, source, source_url, properties = row
            out_lines.append(f"# Node `{args.node_id}`")
            out_lines.append(f"- type: `{ntype}`  label: {label}")
            out_lines.append(f"- source: `{source}`  url: `{source_url or ''}`")
            out_lines.append("")
            cur2 = conn.execute(
                "SELECT edge_type, from_node, to_node FROM edges "
                "WHERE from_node = ? OR to_node = ? LIMIT ?",
                (args.node_id, args.node_id, args.limit))
            edges = cur2.fetchall()
            out_lines.append(f"## Edges ({len(edges)})")
            for et, fn, tn in edges:
                arrow = "→" if fn == args.node_id else "←"
                other = tn if fn == args.node_id else fn
                out_lines.append(f"- ({et}) {arrow} `{other}`")
    else:
        raise SystemExit("Provide --gene / --tissue / --node-id.")

    text = "\n".join(out_lines)
    print(text)
    out = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_query.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"\nReport: {out}")
    conn.close()
    return out


# ----------------------------- stats subcommand ------------------------------

def cmd_stats(_args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    conn = kg_connect()
    cur = conn.execute("SELECT node_type, COUNT(*) FROM nodes GROUP BY node_type "
                        "ORDER BY 2 DESC")
    n_by_type = cur.fetchall()
    cur = conn.execute("SELECT edge_type, COUNT(*) FROM edges GROUP BY edge_type "
                        "ORDER BY 2 DESC")
    e_by_type = cur.fetchall()
    cur = conn.execute("SELECT source, COUNT(*) FROM nodes GROUP BY source")
    by_source = cur.fetchall()
    total_nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    total_edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    cached = conn.execute("SELECT COUNT(*) FROM gene_cache "
                            "WHERE confirmed=1").fetchone()[0]
    runs = conn.execute("SELECT cmd, started_at, ended_at, nodes_added, edges_added "
                          "FROM run_log ORDER BY id DESC LIMIT 10").fetchall()

    lines = [f"# Local KG stats — {KG_PATH.name}", "",
             f"Total nodes: **{total_nodes}**, edges: **{total_edges}**",
             f"Catalog-confirmed gene cache entries: **{cached}**",
             "", "## Nodes by type", "", "| Type | n |", "|---|---:|"]
    for t, n in n_by_type:
        lines.append(f"| {t} | {n} |")
    lines += ["", "## Edges by type", "", "| Type | n |", "|---|---:|"]
    for t, n in e_by_type:
        lines.append(f"| {t} | {n} |")
    lines += ["", "## Nodes by source", "", "| Source | n |", "|---|---:|"]
    for s, n in by_source:
        lines.append(f"| {s} | {n} |")
    lines += ["", "## Recent runs", "",
              "| cmd | started | ended | nodes_added | edges_added |",
              "|---|---|---|---:|---:|"]
    for cmd, st, en, na, ea in runs:
        lines.append(f"| {cmd} | {st} | {en} | {na} | {ea} |")
    text = "\n".join(lines)
    print(text)
    out = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_stats.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    conn.close()
    return out


# ----------------------------- export subcommands ----------------------------

def cmd_export_aql(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    conn = kg_connect()
    ts = time.strftime("%Y%m%d_%H%M%S")
    out_dir = KG_EXPORT_DIR / ts
    out_dir.mkdir(parents=True, exist_ok=True)
    # group nodes by collection prefix
    cur = conn.execute("SELECT id, node_type, label, properties FROM nodes")
    by_coll: dict[str, list[dict]] = defaultdict(list)
    for nid, ntype, label, props in cur.fetchall():
        coll = nid.split("/", 1)[0]
        key = nid.split("/", 1)[1]
        rec = {"_key": key, "_type": ntype, "label": label,
               **(json.loads(props or "{}"))}
        by_coll[coll].append(rec)
    for coll, recs in by_coll.items():
        path = out_dir / f"nodes_{coll}.jsonl"
        with path.open("w") as f:
            for r in recs:
                f.write(json.dumps(r, default=str) + "\n")
    cur = conn.execute("SELECT id, from_node, to_node, edge_type, properties "
                        "FROM edges")
    edges_path = out_dir / "edges.jsonl"
    with edges_path.open("w") as f:
        for eid, fn, tn, et, props in cur.fetchall():
            from_coll, from_key = fn.split("/", 1)
            to_coll, to_key = tn.split("/", 1)
            f.write(json.dumps({
                "_key": eid.replace(":", "_"),
                "_from": fn, "_to": tn,
                "edge_type": et,
                **(json.loads(props or "{}")),
            }, default=str) + "\n")
    print(f"Exported {len(by_coll)} node collections + edges to {out_dir}")
    conn.close()
    return out_dir


def cmd_export_cytoscape(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    conn = kg_connect()
    cur = conn.execute("SELECT id, node_type, label FROM nodes LIMIT ?",
                        (args.limit,))
    elements = []
    keep = set()
    for nid, ntype, label in cur.fetchall():
        elements.append({"data": {"id": nid, "label": label or nid,
                                     "type": ntype}})
        keep.add(nid)
    cur = conn.execute("SELECT id, from_node, to_node, edge_type FROM edges")
    for eid, fn, tn, et in cur.fetchall():
        if fn in keep and tn in keep:
            elements.append({"data": {"id": eid, "source": fn, "target": tn,
                                         "label": et}})
    out = KG_EXPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_cytoscape.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"elements": elements}, indent=2))
    print(f"Wrote: {out}  ({len(elements)} elements)")
    conn.close()
    return out


# ----------------------------- playbook subcommand ---------------------------

def cmd_write_playbook(_a) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "PORTAL_TO_KG_SKILLS.md"
    lines = [
        "# Skill: IGVF Portal → local Knowledge Graph ETL",
        "",
        "Pulls unstructured IGVF Portal entities (AnalysisSets, "
        "MeasurementSets, Samples, Donors, Files, FileSets) by tissue / "
        "gene / assay / lab filters; annotates them against genes and "
        "variants; and persists the result to a local SQLite-backed "
        "knowledge graph that mirrors the IGVF Catalog ArangoDB schema "
        "(nodes + edges + properties + provenance).",
        "",
        "## Local store",
        "",
        "  Data/KG/local_kg.sqlite     SQLite database (gitignored)",
        "  Data/KG/Export/<timestamp>/ arangoimport-compatible JSONL exports",
        "  Data/Cache/KGLocal/         per-run caches",
        "",
        "All endpoint URLs and credentials are resolved through "
        "`Scripts/_endpoints.py` env-var overrides. No URLs, passwords, or "
        "API keys are written to source.",
        "",
        "## Subcommands",
        "",
        "### `pull` — fetch Portal entities by filter",
        "",
        "```bash",
        "# AnalysisSets profiling macrophages, hydrated + walked one hop",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --type AnalysisSet --tissue macrophage --limit 50 --depth 1",
        "",
        "# Multiple types, single tissue",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --types AnalysisSet MeasurementSet --tissue brain --limit 25",
        "",
        "# By assay only",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --type AnalysisSet --assay 'Parse SPLiT-seq' --limit 100 --depth 1",
        "",
        "# By gene name (free-text searchTerm)",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --type AnalysisSet --gene APOE --limit 25 --depth 1",
        "```",
        "",
        "Each Portal entity becomes a node. With `--depth 1` the skill "
        "follows `samples`, `donors`, `files`, `input_file_sets`, "
        "`workflows`, and `documents` cross-references and ingests them "
        "as linked nodes (with appropriate edge labels). Use "
        "`--skip-hydrate` to avoid the per-record full hydration call when "
        "you only want headline metadata.",
        "",
        "### `annotate` — text-mine for genes and variants",
        "",
        "```bash",
        "python3 Scripts/portal_to_kg_skill.py annotate",
        "python3 Scripts/portal_to_kg_skill.py annotate \\",
        "    --gene-list path/to/HGNC_symbols.txt",
        "python3 Scripts/portal_to_kg_skill.py annotate --skip-catalog-confirm",
        "```",
        "",
        "Walks every Portal node in the local KG, scans its description / "
        "summary / aliases / title / name fields for gene-symbol tokens, "
        "rsID strings, and SPDI/HGVS variant identifiers. Each candidate "
        "gene token is confirmed against the IGVF Catalog "
        "(`/api/genes?name=…`) and cached in the `gene_cache` table so "
        "subsequent runs are free. Confirmed Gene and Variant nodes are "
        "added with `mentions_gene` / `mentions_variant` edges back to the "
        "Portal node.",
        "",
        "Pass `--gene-list` to constrain annotation to a curated set "
        "(useful when you only care about a specific GWAS / disease / "
        "curated panel). `--skip-catalog-confirm` retains every "
        "uppercase-token candidate as an unconfirmed Gene node — fast "
        "but noisier.",
        "",
        "### `enrich` — hydrate genes from the IGVF Catalog",
        "",
        "```bash",
        "# All genes in local KG -> pull variants/cCREs/transcripts/etc.",
        "python3 Scripts/portal_to_kg_skill.py enrich --limit 25",
        "",
        "# Just specific symbols",
        "python3 Scripts/portal_to_kg_skill.py enrich --symbols APOE,LDLR,PCSK9",
        "```",
        "",
        "For each Gene node, fetches the Catalog's `variants`, `transcripts`, "
        "`proteins`, `regulatory_elements`, `diseases`, `pathways`, and "
        "`coding_variant_scores` relations and stores them as nodes and "
        "edges in the local KG (with `source = igvf_catalog`).",
        "",
        "### `query` — read interface",
        "",
        "```bash",
        "python3 Scripts/portal_to_kg_skill.py query --gene APOE --limit 50",
        "python3 Scripts/portal_to_kg_skill.py query --tissue macrophage --limit 30",
        "python3 Scripts/portal_to_kg_skill.py query --node-id analysis_sets/IGVFDS3222WCZH",
        "```",
        "",
        "### `stats` — counts and recent ingestion summary",
        "",
        "```bash",
        "python3 Scripts/portal_to_kg_skill.py stats",
        "```",
        "",
        "### `export-aql` — emit ArangoDB-compatible JSONL",
        "",
        "```bash",
        "python3 Scripts/portal_to_kg_skill.py export-aql",
        "```",
        "",
        "Writes one JSONL file per node collection (`nodes_genes.jsonl`, "
        "`nodes_variants.jsonl`, …) plus `edges.jsonl` under "
        "`Data/KG/Export/<timestamp>/`. Each record has the ArangoDB "
        "convention `_key` (and `_from`/`_to` for edges) so a future push "
        "to the actual IGVF Catalog ArangoDB is just `arangoimport "
        "--type jsonl …` per file.",
        "",
        "### `export-cytoscape` — graph visualization",
        "",
        "```bash",
        "python3 Scripts/portal_to_kg_skill.py export-cytoscape --limit 500",
        "```",
        "",
        "## Recommended end-to-end flow",
        "",
        "```bash",
        "# 1. Pull a tissue-specific corpus from the Portal",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --type AnalysisSet --tissue macrophage --limit 100 --depth 1",
        "python3 Scripts/portal_to_kg_skill.py pull \\",
        "    --type MeasurementSet --tissue macrophage --limit 100 --depth 1",
        "",
        "# 2. Mine descriptions for gene/variant mentions, confirm against Catalog",
        "python3 Scripts/portal_to_kg_skill.py annotate",
        "",
        "# 3. Hydrate every confirmed gene from the Catalog",
        "python3 Scripts/portal_to_kg_skill.py enrich --limit 25",
        "",
        "# 4. Inspect the result",
        "python3 Scripts/portal_to_kg_skill.py stats",
        "python3 Scripts/portal_to_kg_skill.py query --gene APOE",
        "",
        "# 5. (Optional) Export for a future push to IGVF Catalog ArangoDB",
        "python3 Scripts/portal_to_kg_skill.py export-aql",
        "```",
        "",
        "## How this chains with other skills",
        "",
        "- After `pull`/`annotate`, the resulting AnalysisSet → Gene "
        "edges feed directly into `Scripts/kg_traversal_skill.py` "
        "for comprehensive remote-Catalog context on each implicated gene.",
        "- The `single_cell_data_skills.py`, `splitseq_pipeline.py`, and "
        "`multiome_10x_pipeline.py` pipelines can ingest the AnalysisSet "
        "manifest from `query` to drive downstream analysis.",
        "- `reference_skill.py validate` can be run against the local "
        "Gene/Variant nodes to surface prior literature.",
        "- The internal Plan → Action → Results → Evaluation orchestrator "
        "uses this skill as the **persistence and indexing layer** for "
        "Portal-side evidence so multi-step plans don't re-fetch the same "
        "metadata across sessions.",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI ---------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="IGVF Portal → local Knowledge Graph ETL skill.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("pull", help="Pull Portal entities and ingest into local KG.")
    s.add_argument("--type", default=None, help="Single Portal type, e.g. AnalysisSet.")
    s.add_argument("--types", nargs="*", default=None,
                    help="Multiple Portal types (overrides --type).")
    s.add_argument("--tissue", default=None,
                    help="searchTerm-based tissue / cell-type / biosample.")
    s.add_argument("--gene", default=None,
                    help="searchTerm-based gene mention.")
    s.add_argument("--assay", default=None,
                    help="preferred_assay_titles facet (e.g. 'Parse SPLiT-seq').")
    s.add_argument("--lab", default=None)
    s.add_argument("--status", default="released")
    s.add_argument("--limit", type=int, default=50)
    s.add_argument("--depth", type=int, default=1,
                    help="How many hops of linked entities to follow (0=none).")
    s.add_argument("--skip-hydrate", action="store_true",
                    help="Skip per-entity hydration call (faster, less complete).")
    s.set_defaults(func=cmd_pull)

    s = sub.add_parser("annotate",
                        help="Text-mine genes/variants for existing Portal nodes.")
    s.add_argument("--gene-list", default=None,
                    help="Optional file of allowed gene symbols, one per line.")
    s.add_argument("--skip-catalog-confirm", action="store_true",
                    help="Don't confirm gene tokens against the Catalog.")
    s.add_argument("--where", default=None,
                    help="Extra SQL WHERE clause on the nodes table.")
    s.set_defaults(func=cmd_annotate)

    s = sub.add_parser("enrich",
                        help="Hydrate Catalog evidence for genes in local KG.")
    s.add_argument("--symbols", default=None,
                    help="Comma-separated gene symbols (default = all genes).")
    s.add_argument("--limit", type=int, default=25)
    s.set_defaults(func=cmd_enrich)

    s = sub.add_parser("query", help="Read interface over the local KG.")
    s.add_argument("--gene", default=None)
    s.add_argument("--tissue", default=None)
    s.add_argument("--node-id", default=None)
    s.add_argument("--limit", type=int, default=50)
    s.set_defaults(func=cmd_query)

    s = sub.add_parser("stats", help="Local KG counts and recent ingestion runs.")
    s.set_defaults(func=cmd_stats)

    s = sub.add_parser("export-aql",
                        help="Emit arangoimport-compatible JSONL.")
    s.set_defaults(func=cmd_export_aql)

    s = sub.add_parser("export-cytoscape",
                        help="Emit cytoscape.js-compatible JSON.")
    s.add_argument("--limit", type=int, default=1000)
    s.set_defaults(func=cmd_export_cytoscape)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/PORTAL_TO_KG_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
