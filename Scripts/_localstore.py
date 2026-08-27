"""Growing local Knowledge Graph + database — the integration framework.

Every time IGVFagent downloads data or analyzes/processes raw data, a little
more of a *local* knowledge graph and database accumulates. Over time the agent
gets faster and more self-contained: repeat queries hit the local store first,
and the graph of everything the user has ever touched keeps growing.

Two backing stores, both additive:

* **Knowledge Graph** — SQLite at ``Data/KG/local_kg.sqlite`` (the SAME file the
  ``portal-kg`` skill builds, so its ``nodes`` / ``edges`` tables grow in place).
  We add genes / variants / regions / datasets / tissues as nodes and the
  relations between them (plus analysis-provenance edges) as edges.
* **Database** — a DuckDB warehouse at ``Data/Warehouse/igvf.duckdb`` (best
  effort: structured ``datasets`` + ``runs`` rows). If DuckDB is unavailable the
  KG still grows.

Every write is idempotent (upsert by deterministic id), so harvesting the same
output twice does not double-count. A ``harvest`` scan walks the on-disk
``Data/`` + ``Docs/<skill>/`` + ``Benchmarks/_data`` trees and ingests anything
new, so even direct-CLI and benchmark runs (which bypass the agent loop) feed
the graph.

Pure-stdlib for the KG path (``sqlite3`` + ``json`` + ``hashlib``); DuckDB is
optional.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import time
from pathlib import Path
from typing import Iterable, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
KG_PATH = ROOT / "Data" / "KG" / "local_kg.sqlite"
WAREHOUSE_DB = ROOT / "Data" / "Warehouse" / "igvf.duckdb"

_NOW = lambda: time.strftime("%Y-%m-%dT%H:%M:%S")  # noqa: E731

# ------------------------------- schema -------------------------------------

# Node/edge tables match portal_to_kg_skill so we grow the SAME graph.
_SCHEMA = """
CREATE TABLE IF NOT EXISTS nodes (
    id TEXT PRIMARY KEY, node_type TEXT NOT NULL, label TEXT,
    source TEXT NOT NULL, source_url TEXT, properties TEXT,
    ingested_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_nodes_type  ON nodes(node_type);
CREATE INDEX IF NOT EXISTS idx_nodes_label ON nodes(label);
CREATE TABLE IF NOT EXISTS edges (
    id TEXT PRIMARY KEY, from_node TEXT NOT NULL, to_node TEXT NOT NULL,
    edge_type TEXT NOT NULL, source TEXT NOT NULL, properties TEXT,
    ingested_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_edges_from ON edges(from_node);
CREATE INDEX IF NOT EXISTS idx_edges_to   ON edges(to_node);
-- Growth provenance (new; additive, does not clash with portal_to_kg.run_log)
CREATE TABLE IF NOT EXISTS download_log (
    id TEXT PRIMARY KEY, source TEXT, accession TEXT, url TEXT,
    path TEXT, n_bytes INTEGER, recorded_at TEXT
);
CREATE TABLE IF NOT EXISTS analysis_log (
    id TEXT PRIMARY KEY, skill TEXT, subcommand TEXT, label TEXT,
    inputs TEXT, outputs TEXT, n_nodes INTEGER, n_edges INTEGER, recorded_at TEXT
);
CREATE TABLE IF NOT EXISTS harvest_ledger (
    key TEXT PRIMARY KEY, kind TEXT, harvested_at TEXT
);
"""


def _connect() -> sqlite3.Connection:
    KG_PATH.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(KG_PATH), timeout=30)
    # WAL + a generous busy timeout so concurrent writers (the agent loop and a
    # direct-CLI skill harvesting at the same time) queue instead of erroring.
    con.execute("PRAGMA busy_timeout=30000")
    try:
        con.execute("PRAGMA journal_mode=WAL")
    except sqlite3.OperationalError:
        pass
    con.executescript(_SCHEMA)
    return con


def _nid(node_type: str, name: str) -> str:
    return f"{node_type}:{name.strip().upper()}"


def _digest(*parts: str) -> str:
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]


# ------------------------------- core upserts -------------------------------

def upsert_node(con, node_type: str, name: str, *, source: str,
                label: Optional[str] = None, source_url: Optional[str] = None,
                properties: Optional[dict] = None) -> str:
    """Insert or merge a node; returns its id. Idempotent."""
    nid = _nid(node_type, name)
    now = _NOW()
    row = con.execute("SELECT properties FROM nodes WHERE id=?", (nid,)).fetchone()
    props = properties or {}
    if row:
        prev = {}
        try:
            prev = json.loads(row[0] or "{}")
        except Exception:
            prev = {}
        prev.update(props)
        prev["observations"] = int(prev.get("observations", 1)) + 1
        con.execute("UPDATE nodes SET properties=?, updated_at=?, "
                    "source=COALESCE(source,?) WHERE id=?",
                    (json.dumps(prev), now, source, nid))
    else:
        props["observations"] = 1
        con.execute(
            "INSERT INTO nodes(id,node_type,label,source,source_url,properties,"
            "ingested_at,updated_at) VALUES(?,?,?,?,?,?,?,?)",
            (nid, node_type, label or name, source, source_url,
             json.dumps(props), now, now))
    return nid


def upsert_edge(con, from_id: str, to_id: str, edge_type: str, *,
                source: str, properties: Optional[dict] = None) -> str:
    eid = _digest(from_id, to_id, edge_type, source)
    if not con.execute("SELECT 1 FROM edges WHERE id=?", (eid,)).fetchone():
        con.execute(
            "INSERT INTO edges(id,from_node,to_node,edge_type,source,"
            "properties,ingested_at) VALUES(?,?,?,?,?,?,?)",
            (eid, from_id, to_id, edge_type, source,
             json.dumps(properties or {}), _NOW()))
    return eid


# --------------------------- entity extraction ------------------------------

# NB: no trailing \b — accessions are frequently followed by '_' (a word char),
# e.g. "GSE162170_multiome...", where \b would fail to match.
_ACC_RE = re.compile(r"\b(IGVF[A-Z]{2}[0-9A-Z]{6,}|ENC[A-Z]{2}[0-9A-Z]{6,}|GSE\d{4,})")
_RS_RE = re.compile(r"\brs\d{2,}\b", re.IGNORECASE)
_REGION_RE = re.compile(r"\bchr[0-9XYMT]{1,2}:[0-9,]+-[0-9,]+\b", re.IGNORECASE)


def entities_from(text: str) -> "list[tuple[str, str]]":
    """Extract (node_type, name) pairs from arbitrary text (args, stdout)."""
    out: "list[tuple[str, str]]" = []
    for m in _ACC_RE.findall(text or ""):
        out.append(("dataset", m))
    for m in _RS_RE.findall(text or ""):
        out.append(("variant", m.lower()))
    for m in _REGION_RE.findall(text or ""):
        out.append(("region", m.replace(",", "")))
    return out


# tool_name -> (node_type, arg_key) for the primary entity of a tool call
_TOOL_ENTITY = {
    "kg_gene": ("gene", "symbol"),
    "kg_variant": ("variant", "variant"),
    "kg_region": ("region", "region"),
    "explain_dataset": ("dataset", "accession_or_url"),
    "annotate_variants": ("variant", "input"),
}


# ------------------------------ public API ----------------------------------

def record_download(source: str, *, accession: Optional[str] = None,
                    url: Optional[str] = None, path: Optional[str] = None,
                    n_bytes: Optional[int] = None,
                    title: Optional[str] = None, con=None) -> dict:
    """Record a data download → dataset node in the KG + warehouse row.

    Pass ``con`` to reuse an open connection (avoids self-locking when called
    inside a larger transaction such as :func:`harvest`).
    """
    own = con is None
    con = con or _connect()
    try:
        key = _digest(source, accession or "", url or "", path or "")
        con.execute(
            "INSERT OR REPLACE INTO download_log(id,source,accession,url,path,"
            "n_bytes,recorded_at) VALUES(?,?,?,?,?,?,?)",
            (key, source, accession, url, path, n_bytes, _NOW()))
        n_nodes = 0
        if accession:
            upsert_node(con, "dataset", accession, source=source,
                        source_url=url,
                        properties={"title": title, "path": path,
                                    "n_bytes": n_bytes})
            n_nodes = 1
        if own:
            con.commit()
    finally:
        if own:
            con.close()
    _warehouse_dataset(source, accession, title, url, path)
    return {"recorded": True, "nodes_added": n_nodes}


def record_analysis(skill: str, *, subcommand: str = "", label: str = "",
                    inputs: Optional[Iterable[str]] = None,
                    outputs: Optional[Iterable[str]] = None,
                    entities: "Optional[Iterable[tuple[str, str]]]" = None,
                    text: str = "", con=None) -> dict:
    """Record an analysis/processing step → provenance + entity nodes/edges.

    ``entities`` is an explicit list of (node_type, name); anything found in
    ``text`` (args + stdout) is added too. An ``analysis`` node links to every
    entity it touched, so the graph accumulates what each run was *about*.
    Pass ``con`` to reuse an open connection.
    """
    own = con is None
    con = con or _connect()
    try:
        ents = list(entities or []) + entities_from(text)
        # dedupe
        seen = set()
        uniq = []
        for t, n in ents:
            k = (t, n.upper())
            if k not in seen and n:
                seen.add(k)
                uniq.append((t, n))
        run_id = _digest(skill, subcommand, label, _NOW(), str(len(uniq)))
        analysis_node = upsert_node(
            con, "analysis", run_id, source=f"igvfagent:{skill}",
            label=f"{skill} {subcommand} {label}".strip(),
            properties={"skill": skill, "subcommand": subcommand,
                        "label": label,
                        "inputs": list(inputs or []),
                        "outputs": list(outputs or [])})
        n_edges = 0
        for ntype, name in uniq:
            nid = upsert_node(con, ntype, name, source=f"igvfagent:{skill}")
            upsert_edge(con, analysis_node, nid, "analyzed",
                        source=f"igvfagent:{skill}")
            n_edges += 1
        con.execute(
            "INSERT OR REPLACE INTO analysis_log(id,skill,subcommand,label,"
            "inputs,outputs,n_nodes,n_edges,recorded_at) VALUES(?,?,?,?,?,?,?,?,?)",
            (run_id, skill, subcommand, label,
             json.dumps(list(inputs or [])), json.dumps(list(outputs or [])),
             len(uniq), n_edges, _NOW()))
        if own:
            con.commit()
    finally:
        if own:
            con.close()
    return {"recorded": True, "entities": len(uniq), "edges_added": n_edges}


# --------------------------- artefact ingestion -----------------------------
#
# The framework recorded that a tool ran and regex-scraped its stdout, but
# never opened the tables those tools produced. So an rE2G run could emit 113
# element->gene linkages and the graph would gain an `analyzed` edge and
# nothing else: the biology stayed in a CSV on disk. These recognisers read
# structured artefacts and turn their rows into real entity-to-entity edges,
# which is what makes retrieval actually accumulate into the knowledge graph.
#
# Conservative by design: a table is ingested only when its columns are an
# unambiguous match. An unrecognised table is left alone rather than guessed
# at, because a wrong edge is worse than a missing one.

_MAX_INGEST_BYTES = 12_000_000      # skip very large tables
_MAX_INGEST_ROWS = 20_000


def _norm_header(fieldnames) -> "dict[str, str]":
    """lowercased-stripped -> original, so `Chrom`/`chrom`/` gene ` all match."""
    return {(f or "").strip().lower(): f for f in (fieldnames or [])}


def _region_id(chrom, start, end) -> str:
    return f"{str(chrom).strip()}:{str(start).strip()}-{str(end).strip()}"


def _first(cols: dict, *names):
    """Original column name for the first alias present."""
    for n in names:
        if n in cols:
            return cols[n]
    return None


# ---------------------------------------------------------------------------
# Recognisers. Each is (label, applies(cols) -> bool, emit(row, cols) -> edges)
# where an edge is (src_type, src_name, edge_type, dst_type, dst_name, props).
#
# Ordered most-specific first; the first match wins. A table matching nothing
# is left alone — a wrong edge is worse than a missing one, and guessing at an
# unknown schema is how a knowledge graph fills with noise.
# ---------------------------------------------------------------------------

def _props(row, cols, keys) -> dict:
    out = {}
    for k in keys:
        c = cols.get(k)
        if c and row.get(c) not in (None, ""):
            out[k] = row.get(c)
    return out


def _emit_element_gene(row, cols):
    gene = (row.get(cols["gene"]) or "").strip()
    ch = _first(cols, "chrom", "chr", "element_chr")
    st = _first(cols, "start", "element_start")
    en = _first(cols, "end", "element_end")
    if not (gene and ch and st and en):
        return []
    c, a, b = row.get(ch), row.get(st), row.get(en)
    if not (c and a and b):
        return []
    pr = _props(row, cols, ("score", "source", "tss", "element_type",
                            "cell_type", "model", "effect_size", "method"))
    return [("regulatory_element", _region_id(c, a, b), "regulates",
             "gene", gene, pr)]


def _emit_gene_disease(row, cols):
    g = _first(cols, "gene_symbol", "gene")
    d = _first(cols, "disease")
    if not (g and d):
        return []
    gene, dis = (row.get(g) or "").strip(), (row.get(d) or "").strip()
    if not (gene and dis):
        return []
    pr = _props(row, cols, ("method", "source", "association_status",
                            "orphanet_association_type", "classification",
                            "pmids"))
    return [("gene", gene, "associated_with", "disease",
             dis.replace("_", " "), pr)]


def _emit_gene_pathway(row, cols):
    g, pw = _first(cols, "gene", "gene_symbol"), _first(cols, "pathway")
    if not (g and pw):
        return []
    gene, path = (row.get(g) or "").strip(), (row.get(pw) or "").strip()
    if not (gene and path):
        return []
    pr = _props(row, cols, ("source", "method", "organism"))
    return [("gene", gene, "member_of_pathway", "pathway", path, pr)]


def _emit_ppi(row, cols):
    a, b = _first(cols, "id_a", "protein_a", "symbol_a"), \
           _first(cols, "id_b", "protein_b", "symbol_b")
    if not (a and b):
        return []
    x, y = (row.get(a) or "").strip(), (row.get(b) or "").strip()
    if not (x and y):
        return []
    pr = _props(row, cols, ("source", "detection_method", "evidence_type",
                            "confidence_score", "pubmed_id", "taxon"))
    return [("protein", x, "interacts_with", "protein", y, pr)]


def _emit_qtl(row, cols):
    g = _first(cols, "gene")
    if not g:
        return []
    gene = (row.get(g) or "").strip()
    if not gene:
        return []
    pr = _props(row, cols, ("qtl_type", "effect_size", "neg_log10_pvalue",
                            "biological_context", "biosample_term", "method"))
    v = _first(cols, "variant", "rsid", "variant_id", "_variant", "name")
    vname = (row.get(v) or "").strip() if v else ""
    if not vname:
        return []
    return [("variant", vname, "qtl_affects_gene", "gene", gene, pr)]


def _emit_variant_gene(row, cols):
    g = _first(cols, "gene", "gene_symbol")
    v = _first(cols, "variant", "rsid", "variant_id", "_variant", "raw")
    if not (g and v):
        return []
    gene, var = (row.get(g) or "").strip(), (row.get(v) or "").strip()
    if not (gene and var):
        return []
    return [("variant", var, "variant_in_gene", "gene", gene,
             _props(row, cols, ("source", "method", "consequence")))]


_RECOGNISERS = [
    # (label, required-any groups, emit)
    ("ppi",           lambda c: bool({"id_a", "protein_a", "symbol_a"} & set(c))
                                and bool({"id_b", "protein_b", "symbol_b"} & set(c)),
                      _emit_ppi),
    ("element_gene",  lambda c: "gene" in c and
                                bool({"chrom", "chr", "element_chr"} & set(c)) and
                                bool({"start", "element_start"} & set(c)) and
                                bool({"end", "element_end"} & set(c)),
                      _emit_element_gene),
    ("gene_disease",  lambda c: bool({"gene_symbol", "gene"} & set(c)) and "disease" in c,
                      _emit_gene_disease),
    ("gene_pathway",  lambda c: bool({"gene", "gene_symbol"} & set(c)) and "pathway" in c,
                      _emit_gene_pathway),
    ("qtl",           lambda c: "gene" in c and "qtl_type" in c,
                      _emit_qtl),
    ("variant_gene",  lambda c: bool({"gene", "gene_symbol"} & set(c)) and
                                bool({"variant", "rsid", "variant_id", "_variant", "raw"} & set(c)),
                      _emit_variant_gene),
]


def ingest_table(path, *, source: str, con=None) -> dict:
    """Ingest one delimited artefact into the KG. Returns counts."""
    import csv as _csv

    added = {"nodes": 0, "edges": 0, "rows": 0, "kind": None}
    _edge_ids: "set[str]" = set()
    _node_ids: "set[str]" = set()
    p = Path(path)
    try:
        if not p.is_file() or p.suffix.lower() not in (".csv", ".tsv", ".txt"):
            return added
        if p.stat().st_size > _MAX_INGEST_BYTES:
            return added
    except OSError:
        return added

    own = con is None
    con = con or _connect()
    try:
        delim = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
        with open(p, newline="", encoding="utf-8", errors="replace") as fh:
            reader = _csv.DictReader(fh, delimiter=delim)
            cols = _norm_header(reader.fieldnames)
            if not cols:
                return added
            match = next(((lbl, emit) for lbl, applies, emit in _RECOGNISERS
                          if applies(cols)), None)
            if not match:
                return added                    # unrecognised: leave alone
            label, emit = match
            added["kind"] = label

            for i, row in enumerate(reader):
                if i >= _MAX_INGEST_ROWS:
                    break
                added["rows"] += 1
                try:
                    edges = emit(row, cols)
                except Exception:
                    continue
                for st, sn, et, dt, dn, pr in edges:
                    sid = upsert_node(con, st, sn,
                                      source=pr.get("source") or source,
                                      properties=pr or None)
                    did = upsert_node(con, dt, dn, source=source)
                    _edge_ids.add(upsert_edge(
                        con, sid, did, et,
                        source=pr.get("source") or source, properties=pr))
                    _node_ids.update((sid, did))

        added["nodes"] = len(_node_ids)
        added["edges"] = len(_edge_ids)
        if own:
            con.commit()
    except (OSError, UnicodeDecodeError, _csv.Error):
        return added
    finally:
        if own:
            con.close()
    return added


def ingest_artifacts(paths, *, source: str, con=None) -> dict:
    total = {"nodes": 0, "edges": 0, "rows": 0, "tables": 0, "kinds": []}
    own = con is None
    con = con or _connect()
    try:
        for path in paths or []:
            r = ingest_table(path, source=source, con=con)
            if r.get("edges"):
                total["tables"] += 1
                total["nodes"] += r["nodes"]
                total["edges"] += r["edges"]
                total["rows"] += r["rows"]
                if r.get("kind"):
                    total["kinds"].append(r["kind"])
        if own:
            con.commit()
    finally:
        if own:
            con.close()
    total["kinds"] = sorted(set(total["kinds"]))
    return total


def record_tool_call(tool_name: str, arguments: dict, *,
                     stdout: str = "", artifacts: Optional[dict] = None) -> dict:
    """Convenience hook for the agent runtime: grow the KG from one tool call."""
    entities = []
    spec = _TOOL_ENTITY.get(tool_name)
    if spec and arguments.get(spec[1]):
        entities.append((spec[0], str(arguments[spec[1]])))
    outputs = []
    for paths in (artifacts or {}).values():
        outputs.extend(paths)
    text = " ".join(str(v) for v in (arguments or {}).values()) + "\n" + (stdout or "")
    result = record_analysis(tool_name, subcommand="tool",
                             inputs=list((arguments or {}).values()),
                             outputs=outputs, entities=entities, text=text)
    # Absorb the CONTENT of the tables the tool produced, not just the fact
    # that it ran. Best-effort: a parse failure must never break the run.
    try:
        result["ingested"] = ingest_artifacts(outputs,
                                              source=f"igvfagent:{tool_name}")
    except Exception:
        result["ingested"] = {"error": True}
    return result


# ----------------------------- warehouse (opt) ------------------------------

def _warehouse_dataset(source, accession, title, url, path) -> None:
    """Best-effort append to the DuckDB warehouse `datasets` table."""
    if not accession:
        return
    try:
        import duckdb  # type: ignore
    except Exception:
        return
    try:
        WAREHOUSE_DB.parent.mkdir(parents=True, exist_ok=True)
        con = duckdb.connect(str(WAREHOUSE_DB))
        con.execute("""CREATE TABLE IF NOT EXISTS datasets(
            dataset_id TEXT PRIMARY KEY, upstream_source TEXT, accession TEXT,
            title TEXT, assay TEXT, modality TEXT, biosample_id TEXT,
            year INTEGER, n_samples INTEGER, pubmed_id TEXT, license TEXT,
            download_url TEXT, bronze_path TEXT, ingested_at TIMESTAMP DEFAULT now())""")
        con.execute("DELETE FROM datasets WHERE dataset_id=?", [accession])
        con.execute("INSERT INTO datasets(dataset_id,upstream_source,accession,"
                    "title,download_url,bronze_path) VALUES(?,?,?,?,?,?)",
                    [accession, source, accession, title, url, path])
        con.close()
    except Exception:
        pass  # KG already grew; warehouse is a bonus


# ------------------------------- harvest ------------------------------------

def harvest(roots: "Optional[list[Path]]" = None) -> dict:
    """Scan on-disk outputs and ingest anything new. Idempotent via ledger.

    Captures: download manifests / raw data under Data/ + Benchmarks/_data,
    and analysis run dirs under Docs/<skill>/<ts>_<label>/.
    """
    roots = roots or [ROOT / "Docs", ROOT / "Data" / "Manifests",
                      ROOT / "Benchmarks" / "_data"]
    con = _connect()
    added = {"analyses": 0, "downloads": 0, "nodes_before": 0,
             "nodes_after": 0, "ingested_tables": 0, "ingested_edges": 0}
    added["nodes_before"] = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]

    def _ledger_has(key: str) -> bool:
        return bool(con.execute("SELECT 1 FROM harvest_ledger WHERE key=?",
                                (key,)).fetchone())

    def _ledger_add(key: str, kind: str):
        con.execute("INSERT OR REPLACE INTO harvest_ledger(key,kind,harvested_at)"
                    " VALUES(?,?,?)", (key, kind, _NOW()))

    # analysis run dirs: Docs/<skill>/<ts>_<label>/ with a summary-ish json
    docs = ROOT / "Docs"
    if docs.is_dir():
        for skill_dir in docs.iterdir():
            if not skill_dir.is_dir():
                continue
            for run_dir in skill_dir.glob("2*_*"):
                if not run_dir.is_dir():
                    continue
                key = str(run_dir.relative_to(ROOT))
                if _ledger_has(key):
                    continue
                label = run_dir.name.split("_", 1)[-1]
                outs = [str(p.relative_to(ROOT)) for p in run_dir.rglob("*")
                        if p.is_file()][:50]
                # entities from any json/tsv/csv filenames + label
                text = label + " " + " ".join(p.name for p in run_dir.rglob("*"))
                record_analysis(skill_dir.name, subcommand="harvest",
                                label=label, outputs=outs, text=text, con=con)
                # Absorb the table CONTENT too, not just the provenance. The
                # agent path does this via record_tool_call; without it here a
                # direct `igvfagent kg gene BRCA1` would record that it ran and
                # discard every linkage, disease and pathway row it retrieved.
                try:
                    ing = ingest_artifacts([str(ROOT / o) for o in outs],
                                           source=f"igvfagent:{skill_dir.name}",
                                           con=con)
                    added["ingested_edges"] = (added.get("ingested_edges", 0)
                                               + ing.get("edges", 0))
                    added["ingested_tables"] = (added.get("ingested_tables", 0)
                                                + ing.get("tables", 0))
                except Exception:
                    pass
                _ledger_add(key, "analysis")
                added["analyses"] += 1

    # downloaded datasets under Benchmarks/_data/<paper>/ and Data/Manifests
    droot = ROOT / "Benchmarks" / "_data"
    if droot.is_dir():
        for paper_dir in droot.iterdir():
            if not paper_dir.is_dir():
                continue
            for f in paper_dir.iterdir():
                if not f.is_file():
                    continue
                key = str(f.relative_to(ROOT))
                if _ledger_has(key):
                    continue
                accs = _ACC_RE.findall(f.name)
                acc = accs[0] if accs else None
                record_download("benchmark_data", accession=acc,
                                path=str(f.relative_to(ROOT)),
                                n_bytes=f.stat().st_size, title=paper_dir.name,
                                con=con)
                # link the dataset to its study node so the graph gains edges
                study = upsert_node(con, "study", paper_dir.name,
                                    source="benchmark", label=paper_dir.name)
                if acc:
                    ds = upsert_node(con, "dataset", acc, source="benchmark_data")
                    upsert_edge(con, study, ds, "includes_dataset",
                                source="benchmark")
                _ledger_add(key, "download")
                added["downloads"] += 1

    added["nodes_after"] = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    con.commit(); con.close()
    return added


def backfill(roots=None, *, limit_files: int = 20000) -> dict:
    """Ingest table content from runs already recorded in the harvest ledger.

    harvest() is idempotent by design: a run directory it has seen is skipped
    forever. That is right for provenance, but it means artefacts produced
    before table ingestion existed can never be absorbed — the ledger says
    "seen" and the rows stay on disk. This walks the output tree and ingests
    the tables regardless of ledger state. Upserts make it safe to repeat.
    """
    roots = roots or [ROOT / "Docs", ROOT / "Data" / "Manifests"]
    con = _connect()
    before_n = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    before_e = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    seen = 0
    per_kind: "dict[str, int]" = {}
    try:
        for root in roots:
            if not root.is_dir():
                continue
            for path in root.rglob("*"):
                if seen >= limit_files:
                    break
                if not path.is_file() or path.suffix.lower() not in (".csv", ".tsv"):
                    continue
                seen += 1
                r = ingest_table(path, source=f"igvfagent:backfill", con=con)
                if r.get("edges"):
                    per_kind[r["kind"]] = per_kind.get(r["kind"], 0) + r["edges"]
        con.commit()
        after_n = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
        after_e = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    finally:
        con.close()
    return {"files_scanned": seen, "nodes_added": after_n - before_n,
            "edges_added": after_e - before_e, "by_kind": per_kind}


def stats() -> dict:
    """Current size of the growing local KG + database."""
    con = _connect()
    try:
        def _c(q):
            return con.execute(q).fetchone()[0]
        by_type = dict(con.execute(
            "SELECT node_type, COUNT(*) FROM nodes GROUP BY node_type "
            "ORDER BY 2 DESC").fetchall())
        return {
            "kg_path": str(KG_PATH.relative_to(ROOT)),
            "nodes": _c("SELECT COUNT(*) FROM nodes"),
            "edges": _c("SELECT COUNT(*) FROM edges"),
            "nodes_by_type": by_type,
            "downloads_logged": _c("SELECT COUNT(*) FROM download_log"),
            "analyses_logged": _c("SELECT COUNT(*) FROM analysis_log"),
        }
    finally:
        con.close()


__all__ = ["record_download", "record_analysis", "record_tool_call",
           "ingest_table", "ingest_artifacts", "backfill",
           "upsert_node", "upsert_edge", "harvest", "stats", "entities_from"]
