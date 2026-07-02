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
    return record_analysis(tool_name, subcommand="tool", inputs=list((arguments or {}).values()),
                           outputs=outputs, entities=entities, text=text)


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
    added = {"analyses": 0, "downloads": 0, "nodes_before": 0, "nodes_after": 0}
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
           "upsert_node", "upsert_edge", "harvest", "stats", "entities_from"]
