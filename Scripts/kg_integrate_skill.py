#!/usr/bin/env python3
"""Merge the three knowledge graphs IGVFagent maintains into the one it shows.

There are three stores, and until now nothing joined them:

  1. Data/KG/local_kg.sqlite            the "IGVF integrated KG" in the UI.
                                         nodes/edges, gene nodes keyed by
                                         SYMBOL (`gene:TP53`).
  2. Data/Proteomics/KG/proteomics.sqlite  1,292,216 BioGRID interactions,
                                         keyed by symbol, id_map EMPTY.
  3. Data/Warehouse/KG/<collection>/     the ArangoDB mirror as parquet.
                                         11,486,365 proteins_proteins edges
                                         keyed by Ensembl protein (ENSP),
                                         already integrating BioGRID + IntAct.

This module does the merge, and it does it through IGVFagent's existing
graph primitives rather than a parallel path: `_localstore.upsert_node` and
`upsert_edge` write every node and edge, so a merged edge is indistinguishable
in shape from one the agent recorded itself, and `harvest_ledger` makes each
merge incremental.

Three properties of those primitives shape the whole design, and getting any
of them wrong would corrupt a graph that is expensive to rebuild:

* `_nid(type, name)` upper-cases, so `gene:TP53` and `gene:tp53` are one
  node. Gene symbols are case-significant in some nomenclatures but IGVF's
  are not, and the existing 37,744 nodes already live under this rule.

* `upsert_edge`'s id is `_digest(from, to, type, SOURCE)` -- source is part
  of edge identity. So BioGRID evidence and IntAct evidence for the same
  protein pair are two rows, not one overwriting the other. That is the
  behaviour we want for a graph that "keeps merging": re-running a merge is
  idempotent, while adding a NEW source adds evidence instead of replacing
  it. It also means the source string must be stable across runs, or every
  re-run duplicates every edge.

* `upsert_node` merges properties and increments `observations`. An entity
  seen by two sources accumulates rather than flip-flopping, so the
  identifiers a merge learns (ENSG, ENSP, UniProt) persist on the node.

Why the identity layer is not optional
--------------------------------------
The three stores do not share an identifier space: symbol, symbol, and ENSP.
The mirror is the authority for reconciling them, and the route is entirely
on Ensembl keys:

    genes(_key=ENSG, symbol)
      -> genes_transcripts(_from=genes/ENSG, _to=transcripts/ENST)
      -> transcripts_proteins(_from=transcripts/ENST, _to=proteins/ENSP)

188,086 transcript-protein rows, so the index is cheap. What it cannot do is
make the mapping one-to-one: 2,493 of 62,646 symbols map to more than one
ENSG, and 121,455 human proteins sit over 68,881 genes. Ambiguity is
therefore RECORDED on the alias row, not resolved by picking a winner --
`kg_identity.ambiguous` is set, and `status` reports how many edges landed on
an ambiguous symbol so a downstream analysis can exclude them.

Synonyms are deliberately NOT used as aliases by default. The mirror carries
289,002 synonym entries over 238,599 distinct strings, and 7,167 of those
strings are claimed by MORE THAN ONE gene symbol. Admitting them by default
would merge distinct entities on a string match -- the one failure a
knowledge graph cannot recover from, because the evidence that they were
ever separate is gone. `--with-synonyms` exists for callers who know their
input vocabulary; ambiguous synonyms are still flagged.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import sqlite3
import sys
import time
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                        # noqa: E402

logger = logging.getLogger("kg_integrate")

MIRROR = Path(os.environ.get(
    "IGVF_KG_MIRROR",
    str(ls.ROOT / "Data" / "Warehouse" / "KG")))
PROTEOMICS = ls.ROOT / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"

# Source strings are part of edge identity (see module docstring), so they are
# constants: changing one duplicates every edge it ever wrote.
SRC_PPI = "PPI-KG:BioGRID"
SRC_MIRROR_PPI = "IGVF-Catalog:proteins_proteins"
SRC_MIRROR_PATHWAY = "IGVF-Catalog:genes_pathways"
SRC_MIRROR_COMPLEX = "IGVF-Catalog:complexes_proteins"

_IDENTITY_SCHEMA = """
CREATE TABLE IF NOT EXISTS kg_identity (
    alias       TEXT NOT NULL,      -- upper-cased, as _nid would store it
    alias_type  TEXT NOT NULL,      -- ensg | ensp | enst | entrez | hgnc
                                    -- | uniprot | symbol | synonym
    symbol      TEXT NOT NULL,      -- canonical gene symbol (the node name)
    ensg        TEXT,               -- the gene this came from, when known
    ambiguous   INTEGER NOT NULL DEFAULT 0,
    source      TEXT NOT NULL,
    built_at    TEXT NOT NULL,
    PRIMARY KEY (alias, alias_type, symbol)
);
CREATE INDEX IF NOT EXISTS idx_ident_alias  ON kg_identity(alias, alias_type);
CREATE INDEX IF NOT EXISTS idx_ident_symbol ON kg_identity(symbol);
-- Per-merge accounting. The ledger says WHETHER a shard was merged; this says
-- what it produced, so `status` can report coverage without re-reading the
-- mirror and a bad merge can be told from an empty one.
CREATE TABLE IF NOT EXISTS kg_merge_log (
    id          TEXT PRIMARY KEY,
    source      TEXT NOT NULL,
    shard       TEXT,
    rows_read   INTEGER,
    edges_added INTEGER,
    resolved    INTEGER,
    unresolved  INTEGER,
    ambiguous   INTEGER,
    merged_at   TEXT NOT NULL
);
"""


def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _con() -> sqlite3.Connection:
    con = ls._connect()
    con.executescript(_IDENTITY_SCHEMA)
    return con


def _duck():
    """A DuckDB connection with the json extension loaded.

    Needed because every list-valued Arango field -- uniprot_ids, synonyms,
    dbxrefs, pmids -- is mirrored as a JSON STRING, not a parquet LIST. Two
    consequences bite immediately: unnest() raises "can only be applied to
    lists, structs and NULL, not VARCHAR", and len() silently returns the
    CHARACTER count. Measuring synonyms with len() reported 7,298,655 where
    the real entry count is 289,002 -- a wrong number that looks plausible,
    which is the dangerous kind.
    """
    try:
        import duckdb
    except ImportError:
        return None
    con = duckdb.connect()
    try:
        con.execute("INSTALL json; LOAD json;")
    except Exception:
        pass                    # already bundled in most builds
    return con


def _shards(collection: str) -> "list[str]":
    return sorted(glob.glob(str(MIRROR / collection / "*.parquet")))


def _pq(collection: str) -> str:
    """A DuckDB read_parquet over every shard of a mirrored collection.

    Globbed rather than enumerated, so a mirror that is still running -- it
    is, and it adds shards for hours -- is picked up on the next merge
    without changing anything here.
    """
    return f"read_parquet('{MIRROR / collection}/*.parquet')"


def _key(ref: str) -> str:
    """`genes/ENSG00000141510` -> `ENSG00000141510`."""
    return str(ref or "").rsplit("/", 1)[-1]


# ─── identity index ─────────────────────────────────────────────────────────

def build_identity(con, *, with_synonyms: bool = False,
                    organism: str = "Homo sapiens") -> dict:
    """Build the alias -> canonical-symbol index from the mirror.

    The mirror is the authority because it is the only one of the three
    stores that carries several identifier systems for the same entity with
    a provenance for each (GENCODE release, UniProt collection).
    """
    duck = _duck()
    if duck is None:
        return {"error": "duckdb is not installed; cannot read the parquet mirror"}
    if not _shards("genes"):
        return {"error": f"no mirrored genes collection under {MIRROR}"}

    now = ls._NOW()
    rows: "list[tuple]" = []

    # 1. Symbols and the direct cross-references the genes collection carries.
    #    A symbol claimed by more than one gene is marked ambiguous on every
    #    row it produces, so a caller can filter without a second query.
    dup = {r[0] for r in duck.execute(
        f"SELECT symbol FROM {_pq('genes')} WHERE organism = ? "
        f"AND symbol IS NOT NULL GROUP BY symbol HAVING COUNT(*) > 1",
        [organism]).fetchall()}
    genes = duck.execute(
        f"SELECT _key, symbol, entrez, hgnc FROM {_pq('genes')} "
        f"WHERE organism = ? AND symbol IS NOT NULL", [organism]).fetchall()
    for ensg, symbol, entrez, hgnc in genes:
        amb = 1 if symbol in dup else 0
        sym = symbol.strip().upper()
        rows.append((sym, "symbol", sym, ensg, amb, "mirror:genes", now))
        rows.append((str(ensg).upper(), "ensg", sym, ensg, amb, "mirror:genes", now))
        for value, kind in ((entrez, "entrez"), (hgnc, "hgnc")):
            if value:
                # The mirror prefixes these: "ENTREZ:3949", "HGNC:6547".
                bare = str(value).split(":")[-1].strip().upper()
                if bare:
                    rows.append((bare, kind, sym, ensg, amb, "mirror:genes", now))

    # 2. ENSP -> symbol, the hop that lets proteins_proteins be merged at all.
    #    Two joins, both on Ensembl keys, no string matching anywhere.
    prot = duck.execute(f"""
        SELECT DISTINCT
               regexp_replace(tp._to,   '^.*/', '') AS ensp,
               regexp_replace(gt._from, '^.*/', '') AS ensg
        FROM {_pq('transcripts_proteins')} tp
        JOIN {_pq('genes_transcripts')} gt
          ON regexp_replace(gt._to, '^.*/', '')
           = regexp_replace(tp._from, '^.*/', '')
        WHERE tp.organism = ?
    """, [organism]).fetchall()
    by_ensg = {g: (s, (1 if s in dup else 0)) for g, s, *_ in
               ((r[0], r[1].strip().upper()) for r in genes)}
    n_prot_unmapped = 0
    for ensp, ensg in prot:
        hit = by_ensg.get(ensg)
        if not hit:
            n_prot_unmapped += 1
            continue
        sym, amb = hit
        rows.append((str(ensp).upper(), "ensp", sym, ensg, amb,
                      "mirror:transcripts_proteins", now))

    # 3. UniProt accessions, so an external PPI source keyed on UniProt can
    #    be merged later without another index build.
    uni = duck.execute(f"""
        SELECT regexp_replace(_key, '^.*/', '') AS ensp,
               unnest(from_json(uniprot_ids, '["VARCHAR"]')) AS ua
        FROM {_pq('proteins')}
        WHERE organism = ? AND json_array_length(uniprot_ids) > 0
    """, [organism]).fetchall()
    ensp_to_sym = {a: (s, am) for a, t, s, _e, am, *_ in
                   ((r[0], r[1], r[2], r[3], r[4]) for r in rows)
                   if t == "ensp"}
    n_uni = 0
    for ensp, ua in uni:
        hit = ensp_to_sym.get(str(ensp).upper())
        if not hit or not ua:
            continue
        sym, amb = hit
        rows.append((str(ua).strip().upper(), "uniprot", sym, None, amb,
                      "mirror:proteins", now))
        n_uni += 1

    # 4. Synonyms, opt-in. 7.3M entries that collide across genes; admitting
    #    them by default would merge distinct entities on a string match.
    n_syn = 0
    if with_synonyms:
        syn = duck.execute(f"""
            SELECT symbol, unnest(from_json(synonyms, '["VARCHAR"]')) AS s
            FROM {_pq('genes')}
            WHERE organism = ? AND synonyms IS NOT NULL
        """, [organism]).fetchall()
        seen_syn: "dict[str, set]" = {}
        for symbol, s in syn:
            if not s:
                continue
            seen_syn.setdefault(str(s).strip().upper(), set()).add(
                symbol.strip().upper())
        for alias, syms in seen_syn.items():
            amb = 1 if len(syms) > 1 else 0
            for sym in syms:
                rows.append((alias, "synonym", sym, None, amb,
                              "mirror:genes.synonyms", now))
                n_syn += 1

    con.executemany(
        "INSERT OR REPLACE INTO kg_identity"
        "(alias,alias_type,symbol,ensg,ambiguous,source,built_at)"
        " VALUES(?,?,?,?,?,?,?)", rows)
    con.commit()

    counts = dict(con.execute(
        "SELECT alias_type, COUNT(*) FROM kg_identity GROUP BY 1").fetchall())
    return {"rows_written": len(rows), "by_type": counts,
             "ambiguous_symbols": len(dup),
             "ensp_without_a_gene": n_prot_unmapped,
             "uniprot_aliases": n_uni, "synonym_aliases": n_syn,
             "organism": organism, "with_synonyms": with_synonyms}


def resolver(con, *, allow_synonyms: bool = False) -> dict:
    """alias -> (symbol, ambiguous), loaded once for a merge pass.

    Held in memory because a merge does millions of lookups and a per-row
    SQLite query would dominate the runtime. Precedence is explicit: an exact
    symbol beats an ENSG/ENSP mapping, which beats a synonym, so a string
    that is both a real symbol and someone else's synonym resolves to the
    gene that owns it.
    """
    order = ["synonym", "uniprot", "entrez", "hgnc", "enst", "ensp", "ensg",
             "symbol"]
    if not allow_synonyms:
        order.remove("synonym")
    out: "dict[str, tuple]" = {}
    rank = {t: i for i, t in enumerate(order)}
    for alias, atype, symbol, amb in con.execute(
            "SELECT alias, alias_type, symbol, ambiguous FROM kg_identity"):
        if atype not in rank:
            continue
        prev = out.get(alias)
        if prev is None or rank[atype] > prev[2]:
            out[alias] = (symbol, bool(amb), rank[atype])
    return {k: (v[0], v[1]) for k, v in out.items()}


# ─── merges ─────────────────────────────────────────────────────────────────

def _merge_log(con, source, shard, read, added, res, unres, amb):
    con.execute(
        "INSERT OR REPLACE INTO kg_merge_log(id,source,shard,rows_read,"
        "edges_added,resolved,unresolved,ambiguous,merged_at)"
        " VALUES(?,?,?,?,?,?,?,?,?)",
        (ls._digest(source, shard or ""), source, shard, read, added, res,
         unres, amb, ls._NOW()))


def merge_ppi(con, *, limit: int = 0, min_symbol_len: int = 1) -> dict:
    """Merge the proteomics BioGRID interactions as gene-gene edges.

    Both endpoints are already gene symbols, so no resolution is needed to
    place them -- but they are still put through the resolver, because that
    is what records whether a symbol is one the IGVF gene collection knows.
    An edge between two symbols the Catalog has never heard of is kept (it is
    real BioGRID evidence) and flagged, rather than dropped or silently
    treated as equivalent to a resolved one.
    """
    if not PROTEOMICS.exists():
        return {"error": f"no proteomics KG at {PROTEOMICS}"}
    res = resolver(con)
    before = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    nodes_before = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    src = sqlite3.connect(f"file:{PROTEOMICS}?mode=ro", uri=True)
    q = "SELECT id_a, id_b, id_type, source FROM interactions"
    if limit:
        q += f" LIMIT {int(limit)}"
    added = resolved = unresolved = ambiguous = read = 0
    seen_pairs = set()
    for a, b, id_type, source in src.execute(q):
        read += 1
        if not a or not b or len(a) < min_symbol_len or len(b) < min_symbol_len:
            continue
        ua, ub = a.strip().upper(), b.strip().upper()
        if ua == ub:
            continue                       # self-interaction, not a graph edge
        # Undirected: store one orientation so A-B and B-A do not become two
        # edges. upsert_edge cannot know they are the same.
        pair = (ua, ub) if ua <= ub else (ub, ua)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        known = []
        for sym in pair:
            hit = res.get(sym)
            if hit:
                resolved += 1
                if hit[1]:
                    ambiguous += 1
                known.append(hit[0])
            else:
                unresolved += 1
                known.append(sym)
        n1 = ls.upsert_node(con, "gene", known[0], source=SRC_PPI)
        n2 = ls.upsert_node(con, "gene", known[1], source=SRC_PPI)
        ls.upsert_edge(con, n1, n2, "interacts_with", source=SRC_PPI,
                        properties={"id_type": id_type, "evidence": source,
                                     "resolved": all(res.get(s) for s in pair)})
        added += 1
    con.commit()
    gained = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0] - before
    nodes_gained = (con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
                     - nodes_before)
    _merge_log(con, SRC_PPI, None, read, gained, resolved, unresolved, ambiguous)
    con.commit()
    return {"source": SRC_PPI, "rows_read": read,
             "edges_asserted_by_source": added,
             "edges_new_in_graph": gained,
             "nodes_new_in_graph": nodes_gained,
             "endpoints_resolved": resolved, "endpoints_unresolved": unresolved,
             "endpoints_on_ambiguous_symbol": ambiguous,
             "distinct_pairs": len(seen_pairs)}


def merge_mirror_ppi(con, *, shards: int = 0, max_rows_per_shard: int = 0,
                      organism: str = "Homo sapiens") -> dict:
    """Merge the Catalog's own protein-protein edges, resolved to genes.

    Incremental by SHARD: the ledger key is the parquet filename, so a mirror
    that is still growing merges only its new files on the next run. That is
    the mechanism by which this graph "keeps merging and keeps growing"
    without re-reading 11.5 million rows each time.

    Edges are collapsed to gene level because that is the resolution the
    integrated graph's other 60,627 edges use. The protein pair is kept in
    the edge properties, so nothing is lost and a later protein-level view
    can be built from the same rows.
    """
    duck = _duck()
    if duck is None:
        return {"error": "duckdb is not installed; cannot read the parquet mirror"}
    files = _shards("proteins_proteins")
    if not files:
        return {"error": "proteins_proteins is not mirrored yet"}
    res = resolver(con)
    todo = [f for f in files
            if not con.execute("SELECT 1 FROM harvest_ledger WHERE key=?",
                                (f"mirror_ppi:{Path(f).name}",)).fetchone()]
    if shards:
        todo = todo[:shards]
    totals = {"shards_available": len(files), "shards_merged": 0,
               "shards_remaining_after": 0, "rows_read": 0,
               "edges_asserted_by_source": 0, "edges_new_in_graph": 0,
               "endpoints_resolved": 0, "endpoints_unresolved": 0,
               "pairs_same_gene": 0}
    for f in todo:
        q = (f"SELECT _from, _to, source, interaction_type, detection_method, "
             f"confidence_value_biogrid, confidence_value_intact, pmids "
             f"FROM read_parquet('{f}') WHERE organism = ?")
        if max_rows_per_shard:
            q += f" LIMIT {int(max_rows_per_shard)}"
        rows = duck.execute(q, [organism]).fetchall()
        read = added = 0
        before_e = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
        seen = set()
        for (fr, to, source, itype, method, cb, ci, pmids) in rows:
            read += 1
            a, b = _key(fr).upper(), _key(to).upper()
            ha, hb = res.get(a), res.get(b)
            for h in (ha, hb):
                if h:
                    totals["endpoints_resolved"] += 1
                else:
                    totals["endpoints_unresolved"] += 1
            if not (ha and hb):
                continue                  # cannot place it on a gene node
            sa, sb = ha[0], hb[0]
            if sa == sb:
                # Two isoforms of one gene interacting collapses to a
                # self-loop at gene level, which is not information.
                totals["pairs_same_gene"] += 1
                continue
            pair = (sa, sb) if sa <= sb else (sb, sa)
            if pair in seen:
                continue
            seen.add(pair)
            n1 = ls.upsert_node(con, "gene", pair[0], source=SRC_MIRROR_PPI)
            n2 = ls.upsert_node(con, "gene", pair[1], source=SRC_MIRROR_PPI)
            ls.upsert_edge(
                con, n1, n2, "interacts_with", source=SRC_MIRROR_PPI,
                properties={"proteins": [a, b], "evidence": source,
                             "interaction_type": itype,
                             "detection_method": method,
                             "confidence_biogrid": cb,
                             "confidence_intact": ci,
                             "pmids": pmids[:3] if pmids else None})
            added += 1
        gained = (con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
                   - before_e)
        con.execute("INSERT OR REPLACE INTO harvest_ledger(key,kind,harvested_at)"
                    " VALUES(?,?,?)",
                    (f"mirror_ppi:{Path(f).name}", "kg_merge", ls._NOW()))
        _merge_log(con, SRC_MIRROR_PPI, Path(f).name, read, gained, 0, 0, 0)
        con.commit()
        totals["shards_merged"] += 1
        totals["rows_read"] += read
        totals["edges_asserted_by_source"] += added
        totals["edges_new_in_graph"] += gained
    totals["shards_remaining_after"] = len(
        [f for f in files
         if not con.execute("SELECT 1 FROM harvest_ledger WHERE key=?",
                             (f"mirror_ppi:{Path(f).name}",)).fetchone()])
    return totals


def merge_mirror_pathways(con, *, organism: str = "Homo sapiens") -> dict:
    """Gene -> pathway membership from the Catalog, as member_of_pathway.

    The integrated graph already has 937 member_of_pathway edges and 480
    pathway nodes from other skills, and this uses the same edge type and
    node type so the two accumulate on the same vertices instead of forming
    a parallel subgraph.
    """
    duck = _duck()
    if duck is None:
        return {"error": "duckdb is not installed"}
    if not _shards("genes_pathways"):
        return {"error": "genes_pathways is not mirrored yet"}
    res = resolver(con)
    names = {}
    if _shards("pathways"):
        names = {_key(r[0]).upper(): r[1] for r in duck.execute(
            f"SELECT _id, name FROM {_pq('pathways')}").fetchall()}
    rows = duck.execute(
        f"SELECT _from, _to FROM {_pq('genes_pathways')}").fetchall()
    added = unresolved = 0
    before_e = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    for fr, to in rows:
        hit = res.get(_key(fr).upper())
        if not hit:
            unresolved += 1
            continue
        pw = _key(to)
        gid = ls.upsert_node(con, "gene", hit[0], source=SRC_MIRROR_PATHWAY)
        pid = ls.upsert_node(con, "pathway", pw, source=SRC_MIRROR_PATHWAY,
                              label=names.get(pw.upper()) or pw)
        ls.upsert_edge(con, gid, pid, "member_of_pathway",
                        source=SRC_MIRROR_PATHWAY)
        added += 1
    con.commit()
    gained = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0] - before_e
    _merge_log(con, SRC_MIRROR_PATHWAY, None, len(rows), gained, 0,
                unresolved, 0)
    con.commit()
    return {"source": SRC_MIRROR_PATHWAY, "rows_read": len(rows),
             "edges_asserted_by_source": added,
             "edges_new_in_graph": gained,
             "unresolved_genes": unresolved,
             "pathway_names_available": len(names)}


# ─── status ─────────────────────────────────────────────────────────────────

def status(con) -> dict:
    n_nodes = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    n_edges = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    ident = con.execute("SELECT COUNT(*) FROM kg_identity").fetchone()[0]
    by_type = dict(con.execute(
        "SELECT alias_type, COUNT(*) FROM kg_identity GROUP BY 1").fetchall())
    by_src = dict(con.execute(
        "SELECT source, COUNT(*) FROM edges GROUP BY 1 ORDER BY 2 DESC"
    ).fetchall())
    merges = con.execute(
        "SELECT source, COUNT(*), SUM(rows_read), SUM(edges_added) "
        "FROM kg_merge_log GROUP BY 1").fetchall()
    avail = {c: len(_shards(c)) for c in
             ("proteins_proteins", "genes", "genes_pathways",
              "complexes_proteins", "proteins", "transcripts_proteins")}
    done = con.execute("SELECT COUNT(*) FROM harvest_ledger WHERE key LIKE "
                        "'mirror_ppi:%'").fetchone()[0]
    return {"nodes": n_nodes, "edges": n_edges,
             "identity_rows": ident, "identity_by_type": by_type,
             "edges_by_source": by_src,
             "merges": [{"source": s, "runs": n, "rows_read": r,
                          "edges_added": e} for s, n, r, e in merges],
             "mirror_shards": avail,
             "mirror_ppi_shards_merged": done,
             "mirror_ppi_shards_total": avail.get("proteins_proteins", 0),
             "mirror_path": str(MIRROR),
             "proteomics_kg": str(PROTEOMICS) if PROTEOMICS.exists() else None}


# ─── commands ───────────────────────────────────────────────────────────────

def cmd_build_index(args) -> int:
    setup_logging()
    con = _con()
    out = build_identity(con, with_synonyms=args.with_synonyms,
                          organism=args.organism)
    if "error" in out:
        print(f"  {out['error']}")
        return 2
    print(f"  identity index built: {out['rows_written']:,} alias rows")
    for k, v in sorted(out["by_type"].items(), key=lambda kv: -kv[1]):
        print(f"    {k:10} {v:>9,}")
    print(f"  symbols claimed by >1 gene (marked ambiguous): "
          f"{out['ambiguous_symbols']:,}")
    print(f"  ENSP with no gene via transcripts: "
          f"{out['ensp_without_a_gene']:,}")
    if not out["with_synonyms"]:
        print("  synonyms NOT indexed (they collide across genes); "
              "--with-synonyms to include them")
    # record_analysis takes no counts -- it derives provenance from
    # inputs/outputs and extracts entities from `text`.
    ls.record_analysis("kg-integrate", subcommand="build-index",
                        label=args.organism,
                        inputs=[str(MIRROR / "genes"),
                                 str(MIRROR / "transcripts_proteins"),
                                 str(MIRROR / "proteins")],
                        text=f"{out['rows_written']} alias rows, "
                              f"{out['ambiguous_symbols']} ambiguous symbols",
                        con=con)
    con.commit()
    return 0


def cmd_merge(args) -> int:
    setup_logging()
    con = _con()
    if con.execute("SELECT COUNT(*) FROM kg_identity").fetchone()[0] == 0:
        print("  no identity index yet — run: igvfagent kg-integrate "
              "build-index")
        return 2
    outs = []
    if args.what in ("ppi", "all"):
        outs.append(merge_ppi(con, limit=args.limit))
    if args.what in ("mirror-ppi", "all"):
        outs.append(merge_mirror_ppi(con, shards=args.shards,
                                      max_rows_per_shard=args.limit))
    if args.what in ("pathways", "all"):
        outs.append(merge_mirror_pathways(con))
    rc = 0
    for out in outs:
        if "error" in out:
            print(f"  {out['error']}")
            rc = 2
            continue
        for k, v in out.items():
            print(f"    {k:32} {v:,}" if isinstance(v, int)
                  else f"    {k:32} {v}")
        print()
    s = status(con)
    print(f"  graph now: {s['nodes']:,} nodes, {s['edges']:,} edges")
    con.commit()
    ls.record_analysis("kg-integrate", subcommand=f"merge {args.what}",
                        text=f"{s['nodes']} nodes, {s['edges']} edges after "
                              f"merging {args.what}", con=con)
    con.commit()
    return rc


def cmd_status(args) -> int:
    setup_logging()
    s = status(_con())
    print(f"IGVF integrated KG: {s['nodes']:,} nodes, {s['edges']:,} edges")
    print(f"  identity index: {s['identity_rows']:,} aliases "
          f"{s['identity_by_type']}")
    print("  edges by source:")
    for k, v in list(s["edges_by_source"].items())[:12]:
        print(f"    {k:34} {v:>9,}")
    print(f"  mirror at {s['mirror_path']}")
    for k, v in sorted(s["mirror_shards"].items()):
        print(f"    {k:24} {v:>6,} shards")
    print(f"  proteins_proteins merged: "
          f"{s['mirror_ppi_shards_merged']:,} / "
          f"{s['mirror_ppi_shards_total']:,} shards")
    if s["merges"]:
        print("  merge history:")
        for m in s["merges"]:
            print(f"    {m['source']:34} {m['runs']:>4} run(s), "
                  f"{(m['rows_read'] or 0):>10,} read, "
                  f"{(m['edges_added'] or 0):>9,} edges")
    if not s["proteomics_kg"]:
        print("  NOTE: no proteomics PPI-KG on this host")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent kg-integrate",
        description="Merge the PPI-KG and the ArangoDB mirror into the "
                     "IGVF integrated KG, incrementally.")
    sub = p.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build-index", help="Build the alias -> gene-symbol "
                                            "identity index from the mirror.")
    b.add_argument("--with-synonyms", action="store_true",
                   help="Also index gene synonyms. They collide across "
                        "genes; only use it when the input vocabulary is known.")
    b.add_argument("--organism", default="Homo sapiens")

    m = sub.add_parser("merge", help="Merge a source into the graph.")
    m.add_argument("what", choices=["ppi", "mirror-ppi", "pathways", "all"])
    m.add_argument("--limit", type=int, default=0,
                   help="Cap rows read (per shard for mirror-ppi). 0 = all.")
    m.add_argument("--shards", type=int, default=0,
                   help="Cap parquet shards merged this run. 0 = all "
                        "not-yet-merged.")

    sub.add_parser("status", help="Coverage, merge history and what remains.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"build-index": cmd_build_index, "merge": cmd_merge,
            "status": cmd_status}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
