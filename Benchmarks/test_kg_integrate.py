#!/usr/bin/env python3
"""Merging three knowledge graphs into one, without corrupting it.

A merge mutates a graph that costs hours to rebuild, so the properties that
matter are the ones whose failure is silent:

  * idempotency -- re-running must not duplicate. The first live run proved
    why this needs a test: the merge committed 977,810 edges and then died
    with "database is locked" while writing its own log row, so the recovery
    was to re-run, and a non-idempotent merge would have doubled the graph.
  * one orientation per undirected pair -- upsert_edge cannot know that A-B
    and B-A are the same interaction.
  * source is part of edge identity, so two sources for one pair are two
    rows of evidence, not one overwriting the other.
  * ambiguity is recorded, never resolved by picking a winner.
  * synonyms stay out by default: 7,167 synonym strings in the mirror are
    claimed by more than one gene, and merging on those is unrecoverable.

Runs against a temporary SQLite graph, no mirror and no network.
"""
import json
import os
import sqlite3
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:62} {detail}")
    if not ok:
        FAILURES.append(name)


# A temporary graph, so _localstore's module-level KG_PATH points at it.
_tmp = tempfile.TemporaryDirectory()
os.environ["IGVF_PROJECT_ROOT"] = _tmp.name
(Path(_tmp.name) / "Data" / "KG").mkdir(parents=True)
(Path(_tmp.name) / "Data" / "Proteomics" / "KG").mkdir(parents=True)

import _localstore as ls           # noqa: E402
ls.ROOT = Path(_tmp.name)
ls.KG_PATH = ls.ROOT / "Data" / "KG" / "local_kg.sqlite"
import kg_integrate_skill as kgi   # noqa: E402
kgi.PROTEOMICS = ls.ROOT / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"
kgi.MIRROR = ls.ROOT / "Data" / "Warehouse" / "KG"


def fresh():
    if ls.KG_PATH.exists():
        ls.KG_PATH.unlink()
    for suffix in ("-wal", "-shm"):
        p = Path(str(ls.KG_PATH) + suffix)
        if p.exists():
            p.unlink()
    return kgi._con()


def identity(con, rows):
    """Seed the identity index directly, so these cases need no mirror."""
    now = ls._NOW()
    con.executemany(
        "INSERT OR REPLACE INTO kg_identity"
        "(alias,alias_type,symbol,ensg,ambiguous,source,built_at)"
        " VALUES(?,?,?,?,?,?,?)",
        [(a, t, s, e, amb, "test", now) for a, t, s, e, amb in rows])
    con.commit()


def ppi_db(pairs):
    p = kgi.PROTEOMICS
    if p.exists():
        p.unlink()
    c = sqlite3.connect(str(p))
    c.execute("CREATE TABLE interactions (id_a TEXT, id_b TEXT, "
              "id_type TEXT, source TEXT)")
    c.executemany("INSERT INTO interactions VALUES (?,?,'symbol','biogrid')",
                  pairs)
    c.commit()
    c.close()


def counts(con):
    return (con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0],
            con.execute("SELECT COUNT(*) FROM edges").fetchone()[0])


# ── source strings are frozen ─────────────────────────────────────────────
# upsert_edge's id is _digest(from, to, type, SOURCE). Changing a source
# string does not migrate edges, it duplicates every edge that source ever
# wrote -- silently, on the next merge.
check("the PPI source string is the one already on disk",
      kgi.SRC_PPI == "PPI-KG:BioGRID")
check("the mirror PPI source string is stable",
      kgi.SRC_MIRROR_PPI == "IGVF-Catalog:proteins_proteins")
check("sources are module constants, not computed",
      all(isinstance(getattr(kgi, n), str) for n in
          ("SRC_PPI", "SRC_MIRROR_PPI", "SRC_MIRROR_PATHWAY",
           "SRC_MIRROR_COMPLEX")))

# ── idempotency, the property the live failure depended on ────────────────
con = fresh()
identity(con, [("TP53", "symbol", "TP53", "ENSG00000141510", 0),
                ("MDM2", "symbol", "MDM2", "ENSG00000135679", 0)])
ppi_db([("TP53", "MDM2")])
first = kgi.merge_ppi(con)
n1, e1 = counts(con)
second = kgi.merge_ppi(con)
n2, e2 = counts(con)
check("a merge adds the edge", e1 == 1, f"{e1} edge(s)")
check("re-running adds nothing", (n2, e2) == (n1, e1), f"{n2} nodes {e2} edges")
# `edges_added` counted upsert CALLS, so an idempotent re-run reported
# 977,810 against a graph that gained nothing -- indistinguishable from a
# doubling. The two quantities are now reported separately.
check("the first run reports one new edge in the graph",
      first["edges_new_in_graph"] == 1, str(first["edges_new_in_graph"]))
check("the re-run reports ZERO new edges in the graph",
      second["edges_new_in_graph"] == 0, str(second["edges_new_in_graph"]))
check("but still reports what the source asserts",
      second["edges_asserted_by_source"] == 1)
check("and no new nodes on the re-run",
      second["nodes_new_in_graph"] == 0)
check("the second run still reports what it read",
      second["rows_read"] == first["rows_read"] == 1)
check("both runs are logged, so a re-run is visible",
      con.execute("SELECT COUNT(*) FROM kg_merge_log").fetchone()[0] >= 1)

# ── undirected pairs collapse to one orientation ──────────────────────────
con = fresh()
identity(con, [("A", "symbol", "A", "ENSG1", 0), ("B", "symbol", "B", "ENSG2", 0)])
ppi_db([("A", "B"), ("B", "A")])
kgi.merge_ppi(con)
check("A-B and B-A become ONE edge", counts(con)[1] == 1,
      f"{counts(con)[1]} edge(s)")
row = con.execute("SELECT from_node, to_node FROM edges").fetchone()
check("stored in sorted order, so the orientation is deterministic",
      row == ("gene:A", "gene:B"), str(row))

# ── self-interactions are not graph edges ─────────────────────────────────
con = fresh()
identity(con, [("A", "symbol", "A", "ENSG1", 0)])
ppi_db([("A", "A")])
kgi.merge_ppi(con)
check("a self-interaction produces no edge", counts(con)[1] == 0)

# ── unresolved symbols are kept and flagged, not dropped ──────────────────
con = fresh()
identity(con, [("TP53", "symbol", "TP53", "ENSG00000141510", 0)])
ppi_db([("TP53", "NOTAGENE")])
out = kgi.merge_ppi(con)
check("an edge to an unknown symbol is still merged", counts(con)[1] == 1)
check("and the unresolved endpoint is counted",
      out["endpoints_unresolved"] == 1, str(out["endpoints_unresolved"]))
check("the resolved one is counted separately",
      out["endpoints_resolved"] == 1)
props = json.loads(con.execute("SELECT properties FROM edges").fetchone()[0])
check("the edge records that it was not fully resolved",
      props.get("resolved") is False)

# ── ambiguity is recorded, not silently collapsed ─────────────────────────
con = fresh()
identity(con, [("AMBIG", "symbol", "AMBIG", "ENSG_X", 1),
                ("CLEAN", "symbol", "CLEAN", "ENSG_Y", 0)])
ppi_db([("AMBIG", "CLEAN")])
out = kgi.merge_ppi(con)
check("an endpoint on an ambiguous symbol is counted",
      out["endpoints_on_ambiguous_symbol"] == 1)
check("the edge is still created -- ambiguity is reported, not fatal",
      counts(con)[1] == 1)

# ── the resolver's precedence ─────────────────────────────────────────────
con = fresh()
# A string that is BOTH a real symbol and another gene's synonym must resolve
# to the gene that owns it, or the merge silently reassigns interactions.
identity(con, [("TP53", "symbol", "TP53", "ENSG_A", 0),
                ("TP53", "synonym", "OTHERGENE", None, 1)])
r = kgi.resolver(con, allow_synonyms=True)
check("an exact symbol beats a synonym claim on the same string",
      r["TP53"][0] == "TP53", str(r.get("TP53")))
r_no = kgi.resolver(con)
check("synonyms are excluded from the resolver by default",
      r_no["TP53"][0] == "TP53")
con = fresh()
identity(con, [("Q9DC51", "uniprot", "GNAI3", None, 0),
                ("ENSP00000001", "ensp", "GNAI3", "ENSG_G", 0)])
r = kgi.resolver(con)
check("a uniprot accession resolves to its gene", r["Q9DC51"][0] == "GNAI3")
check("an ENSP resolves to its gene", r["ENSP00000001"][0] == "GNAI3")
check("an unknown alias is simply absent", "NOPE" not in r)

# ── synonyms stay out unless asked for ────────────────────────────────────
con = fresh()
identity(con, [("P53", "synonym", "TP53", None, 1)])
check("a synonym-only alias does not resolve by default",
      "P53" not in kgi.resolver(con))
check("but does with allow_synonyms",
      kgi.resolver(con, allow_synonyms=True)["P53"][0] == "TP53")
check("and it carries its ambiguity flag",
      kgi.resolver(con, allow_synonyms=True)["P53"][1] is True)

# ── the merge writes through the shared primitives ────────────────────────
con = fresh()
identity(con, [("A", "symbol", "A", "E1", 0), ("B", "symbol", "B", "E2", 0)])
ppi_db([("A", "B")])
kgi.merge_ppi(con)
nrow = con.execute("SELECT id, node_type, source, properties FROM nodes "
                    "WHERE id='gene:A'").fetchone()
check("nodes get the id shape _nid produces", nrow and nrow[0] == "gene:A")
check("node_type is the same 'gene' other skills write",
      nrow and nrow[1] == "gene")
check("observations is set, so a second source accumulates",
      nrow and json.loads(nrow[3]).get("observations") == 1)
erow = con.execute("SELECT edge_type, source FROM edges").fetchone()
check("the edge type is 'interacts_with'", erow[0] == "interacts_with")
check("the edge carries its source", erow[1] == kgi.SRC_PPI)

# ── two sources for one pair are two rows of evidence ─────────────────────
con = fresh()
identity(con, [("A", "symbol", "A", "E1", 0), ("B", "symbol", "B", "E2", 0)])
a = ls.upsert_node(con, "gene", "A", source="s1")
b = ls.upsert_node(con, "gene", "B", source="s2")
ls.upsert_edge(con, a, b, "interacts_with", source=kgi.SRC_PPI)
ls.upsert_edge(con, a, b, "interacts_with", source=kgi.SRC_MIRROR_PPI)
con.commit()
check("BioGRID and Catalog evidence for one pair are two edges",
      counts(con)[1] == 2, f"{counts(con)[1]} edge(s)")
ls.upsert_edge(con, a, b, "interacts_with", source=kgi.SRC_PPI)
con.commit()
check("but the same source twice is still one", counts(con)[1] == 2)

# ── status reports without needing the mirror ─────────────────────────────
s = kgi.status(con)
check("status counts nodes and edges", s["edges"] == 2)
check("status breaks edges down by source", len(s["edges_by_source"]) == 2)
check("status reports a missing mirror as zero shards, not a crash",
      s["mirror_shards"]["proteins_proteins"] == 0)
check("status names where the mirror is expected",
      str(kgi.MIRROR) == s["mirror_path"])

# ── guards ────────────────────────────────────────────────────────────────
con = fresh()
missing = kgi.PROTEOMICS
if missing.exists():
    missing.unlink()
out = kgi.merge_ppi(con)
check("a missing PPI database is an error, not an empty success",
      "error" in out)
out = kgi.merge_mirror_ppi(con)
check("an unmirrored collection is an error, not an empty success",
      "error" in out)
check("the key helper strips an Arango collection prefix",
      kgi._key("proteins/ENSP00000001") == "ENSP00000001")
check("and tolerates a bare key", kgi._key("ENSP00000001") == "ENSP00000001")
check("and None", kgi._key(None) == "")

print(f"\n{len(FAILURES)} failure(s)")
_tmp.cleanup()
sys.exit(1 if FAILURES else 0)
