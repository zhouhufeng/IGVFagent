"""What the integrated knowledge graph contains, and what it does not yet.

The Knowledge Graph Explorer reads ``Data/KG/local_kg.sqlite`` only. The
Catalog mirror (``kg-mirror pull``, run by ``Deploy/mirror-kg.sh``) writes
parquet shards under ``Data/Warehouse/KG/<collection>/`` and never touches
that database; the only route from the mirror into the graph is
``kg-integrate merge``. So a collection can be fully mirrored and still
invisible in the Explorer, and nothing on screen said so. This panel puts the
two side by side, per collection:

    mirrored   rows on disk vs documents in the Catalog, and the pull status
               from the resume cursor ``_state/<collection>.json``
    in graph   whether a merger exists for it and whether it has run

and, where the operator allows it, starts a merge in the background.

Read-only against ``local_kg.sqlite`` (``mode=ro``): rendering the panel must
never create tables or take a write lock that a running merge is waiting on.
"""
from __future__ import annotations

import csv
import json
import os
import sqlite3
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

try:
    from igvfagent import _pathguard  # type: ignore
except Exception:
    import _pathguard  # type: ignore

ROOT = _pathguard.project_root()
KG_DB = ROOT / "Data" / "KG" / "local_kg.sqlite"
MIRROR = Path(os.environ.get("IGVF_KG_MIRROR",
                             str(ROOT / "Data" / "Warehouse" / "KG")))
LOG_DIR = ROOT / "Docs" / "Logs"
PID_FILE = LOG_DIR / "kg_integrate_ui.pid"

# Mirror collection -> the kg-integrate step that brings it into the graph.
# The source strings are the ones kg_integrate_skill writes to kg_merge_log.
MERGERS = {
    "proteins_proteins":  ("merge mirror-ppi", "IGVF-Catalog:proteins_proteins"),
    "genes_pathways":     ("merge pathways", "IGVF-Catalog:genes_pathways"),
    "pathways":           ("merge pathways", "IGVF-Catalog:genes_pathways"),
    "complexes":          ("merge complexes", "IGVF-Catalog:complexes_proteins"),
    "complexes_proteins": ("merge complexes", "IGVF-Catalog:complexes_proteins"),
    "genes":              ("build-index", None),
    "transcripts_proteins": ("build-index", None),
    "proteins":           ("build-index", None),
}
ACTIONS = {
    "build-index":      ["kg-integrate", "build-index"],
    "merge ppi":        ["kg-integrate", "merge", "ppi"],
    "merge mirror-ppi": ["kg-integrate", "merge", "mirror-ppi"],
    "merge pathways":   ["kg-integrate", "merge", "pathways"],
    "merge complexes":  ["kg-integrate", "merge", "complexes"],
    "merge all":        ["kg-integrate", "merge", "all"],
}


def _ro(path: Path) -> Optional[sqlite3.Connection]:
    if not path.exists():
        return None
    try:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=5)
        con.execute("PRAGMA busy_timeout=5000")
        return con
    except sqlite3.Error:
        return None


def _q(con, sql: str, default: Any) -> Any:
    try:
        return con.execute(sql).fetchall()
    except sqlite3.Error:
        return default


_CACHE: "dict[str, Any]" = {}
CACHE_SECONDS = 60


def graph_status() -> "dict[str, Any]":
    """Node/edge counts, edges by source, merge history. Never writes.

    Cached for a minute (and until the database file changes): the counts
    take ~1 s on the deployed 2M-edge graph, and Streamlit runs every tab's
    body on every interaction, including a chat message.
    """
    try:
        stamp = KG_DB.stat().st_mtime
    except OSError:
        stamp = None
    hit = _CACHE.get("graph")
    if hit and hit[0] == stamp and time.time() - hit[1] < CACHE_SECONDS:
        return hit[2]
    out = _graph_status()
    _CACHE["graph"] = (stamp, time.time(), out)
    return out


def _graph_status() -> "dict[str, Any]":
    con = _ro(KG_DB)
    if con is None:
        return {"exists": False}
    try:
        nodes = _q(con, "SELECT COUNT(*) FROM nodes", [(0,)])[0][0]
        edges = _q(con, "SELECT COUNT(*) FROM edges", [(0,)])[0][0]
        by_src = _q(con, "SELECT source, COUNT(*) FROM edges GROUP BY 1 "
                         "ORDER BY 2 DESC", [])
        identity = _q(con, "SELECT COUNT(*) FROM kg_identity", [(0,)])[0][0]
        merges = _q(con, "SELECT source, COUNT(*), SUM(rows_read), "
                         "SUM(edges_added), MAX(merged_at) FROM kg_merge_log "
                         "GROUP BY 1", [])
        ppi_done = _q(con, "SELECT COUNT(*) FROM harvest_ledger WHERE key "
                           "LIKE 'mirror_ppi:%'", [(0,)])[0][0]
    finally:
        con.close()
    return {"exists": True, "nodes": nodes, "edges": edges,
            "edges_by_source": by_src, "identity_rows": identity,
            "merges": {m[0]: {"runs": m[1], "rows_read": m[2] or 0,
                              "edges_added": m[3] or 0, "last": m[4]}
                       for m in merges},
            "mirror_ppi_shards_merged": ppi_done}


def mirror_status() -> "list[dict[str, Any]]":
    """One row per Catalog collection: documents, mirrored rows, status."""
    inv: "dict[str, dict]" = {}
    p = MIRROR / "_inventory.csv"
    if p.exists():
        with p.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                inv[r["collection"]] = r
    states: "dict[str, dict]" = {}
    sd = MIRROR / "_state"
    if sd.is_dir():
        for f in sd.glob("*.json"):
            try:
                s = json.loads(f.read_text())
                states[s.get("collection") or f.stem] = s
            except Exception:
                continue
    rows = []
    for name in sorted(set(inv) | set(states)):
        i, s = inv.get(name, {}), states.get(name, {})
        docs = int(i["documents"]) if (i.get("documents") or "").isdigit() else None
        written = int(s.get("rows_written") or 0)
        rows.append({
            "collection": name,
            "kind": i.get("type", ""),
            "catalog_docs": docs,
            "mirrored_rows": written,
            "pct": (100.0 * written / docs) if docs else None,
            "pull": s.get("status") or "not mirrored",
            "gb": (int(i["bytes"]) / 1e9) if (i.get("bytes") or "").isdigit() else None,
        })
    return rows


def coverage(mirror: "list[dict]", graph: "dict[str, Any]") -> "list[dict]":
    """Join the two: where each collection stands on the way into the graph."""
    merges = graph.get("merges", {})
    out = []
    for r in mirror:
        step, source = MERGERS.get(r["collection"], (None, None))
        if step is None:
            state = "no merger yet"
        elif step == "build-index":
            state = ("in identity index" if graph.get("identity_rows")
                     else "needs build-index")
        elif source in merges:
            state = "merged"
        elif r["mirrored_rows"]:
            state = f"ready: {step}"
        else:
            state = "waiting for mirror"
        out.append({**r, "in_graph": state, "step": step or ""})
    return out


def _running_pid() -> Optional[int]:
    try:
        pid = int(PID_FILE.read_text().split()[0])
        os.kill(pid, 0)
        return pid
    except Exception:
        return None


def start_action(action: str, cli_argv: "list[str]") -> "tuple[bool, str]":
    """Start one kg-integrate step detached. One at a time."""
    if action not in ACTIONS:
        return False, f"unknown action {action!r}"
    pid = _running_pid()
    if pid:
        return False, f"a KG merge is already running (pid {pid})"
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"kg_integrate_ui_{time.strftime('%Y%m%d_%H%M%S')}.log"
    fh = log.open("w")
    proc = subprocess.Popen([*cli_argv, *ACTIONS[action]], stdout=fh,
                            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                            cwd=str(ROOT), start_new_session=True)
    PID_FILE.write_text(f"{proc.pid} {log}\n")
    return True, str(log)


def _last_log() -> Optional[Path]:
    logs = sorted(LOG_DIR.glob("kg_integrate_ui_*.log"))
    return logs[-1] if logs else None


def _fmt(n: Optional[float]) -> str:
    return "" if n is None else f"{n:,.0f}"


def render_streamlit_panel(st, *, allow_merge: bool = False,
                           cli_argv: "Optional[list[str]]" = None) -> None:
    st.markdown(
        "### 🗂 Sources\n"
        "The Explorer shows the **integrated graph** only. The Catalog "
        "mirror downloads collections to disk; a collection reaches the "
        "graph only when a merge step exists for it and has run.")
    g = graph_status()
    if not g.get("exists"):
        st.warning(f"No integrated graph at `{KG_DB.relative_to(ROOT)}` yet.")
    else:
        c = st.columns(3)
        c[0].metric("Nodes", f"{g['nodes']:,}")
        c[1].metric("Edges", f"{g['edges']:,}")
        c[2].metric("Identity aliases", f"{g['identity_rows']:,}")
        if g["edges_by_source"]:
            with st.expander("Edges by source", expanded=False):
                st.dataframe([{"source": s, "edges": n}
                              for s, n in g["edges_by_source"]],
                             hide_index=True)

    rows = coverage(mirror_status(), g if g.get("exists") else {})
    if not rows:
        st.info("No Catalog mirror on this host (no `_inventory.csv` or "
                f"resume cursors under `{MIRROR}`).")
    else:
        mirrored = [r for r in rows if r["mirrored_rows"]]
        merged = [r for r in rows if r["in_graph"] in
                  ("merged", "in identity index")]
        done = sum(1 for r in rows if r["pull"] == "done")
        c = st.columns(3)
        c[0].metric("Collections mirrored", f"{done} done",
                    f"{len(mirrored)} with data", delta_color="off")
        c[1].metric("Rows on disk",
                    f"{sum(r['mirrored_rows'] for r in rows) / 1e6:,.0f} M")
        c[2].metric("Reach the graph", f"{len(merged)} of {len(rows)}")
        running = [r["collection"] for r in rows if r["pull"] == "in_progress"]
        if running:
            r0 = next(r for r in rows if r["collection"] == running[0])
            pct = r0["pct"] or 0.0
            st.progress(min(pct / 100.0, 1.0),
                        text=f"Pulling {r0['collection']}: "
                             f"{r0['mirrored_rows']:,} of "
                             f"{_fmt(r0['catalog_docs'])} ({pct:.1f}%)")
        only = st.toggle("Only collections with data on disk", value=True,
                         key="kgsrc_only_mirrored")
        table = [{
            "collection": r["collection"],
            "type": r["kind"],
            "Catalog docs": _fmt(r["catalog_docs"]),
            "mirrored": _fmt(r["mirrored_rows"]),
            "%": "" if r["pct"] is None else f"{r['pct']:.0f}",
            "pull": r["pull"],
            "in graph": r["in_graph"],
        } for r in rows if r["mirrored_rows"] or not only]
        st.dataframe(table, hide_index=True)
        st.caption("**no merger yet**: mirrored, but nothing converts it into "
                   "graph edges, so the Explorer cannot show it. Variant and "
                   "regulatory-element collections are in this state.")

    if g.get("merges"):
        with st.expander("Merge history", expanded=False):
            st.dataframe([{"source": s, **m} for s, m in g["merges"].items()],
                         hide_index=True)

    st.markdown("#### Merge into the graph")
    if not allow_merge or not cli_argv:
        st.caption("Merges are run by the operator on this deployment:\n\n"
                   "```\nigvfagent kg-integrate build-index\n"
                   "igvfagent kg-integrate merge all\n"
                   "igvfagent kg-integrate status\n```")
        return
    pid = _running_pid()
    c1, c2 = st.columns([2, 1])
    with c1:
        action = st.selectbox("Step", list(ACTIONS), key="kgsrc_action",
                              disabled=bool(pid))
    with c2:
        st.write("")
        if st.button("Run", key="kgsrc_run", disabled=bool(pid)):
            ok, msg = start_action(action, cli_argv)
            (st.success if ok else st.error)(
                f"Started — log `{msg}`" if ok else msg)
            pid = _running_pid()
    if pid:
        st.info(f"A merge is running (pid {pid}). Refresh to follow it; the "
                "Explorer shows the new edges once it finishes.")
    log = _last_log()
    if log is not None:
        with st.expander(f"Last merge log · {log.name}", expanded=bool(pid)):
            try:
                st.code("".join(log.read_text().splitlines(True)[-40:])
                        or "(empty)")
            except Exception as e:
                st.caption(f"(could not read {log}: {e})")
