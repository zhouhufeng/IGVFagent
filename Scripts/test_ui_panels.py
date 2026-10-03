"""Offline checks for the run browser and the KG sources panel.

    python3 Scripts/test_ui_panels.py

No Streamlit needed: the panels' data functions are tested directly, on a
temporary workspace (IGVF_PROJECT_ROOT) so nothing real is read or written.
"""
from __future__ import annotations

import importlib
import json
import os
import sqlite3
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

FAIL = []


def check(name: str, cond: bool) -> None:
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
    if not cond:
        FAIL.append(name)


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="uipanels_")).resolve()
    os.environ["IGVF_PROJECT_ROOT"] = str(tmp)
    os.environ["IGVF_KG_MIRROR"] = str(tmp / "Data/Warehouse/KG")
    import _pathguard
    importlib.reload(_pathguard)
    import data_browser as db
    import kg_sources as ks
    importlib.reload(db)
    importlib.reload(ks)

    # ── workspace: two runs, a Portal fetch, static docs, a secret ──
    a = tmp / "Docs/scE2G/20260922_230951_mine"
    b = tmp / "Docs/SingleCell/20260908_231139_theirs"
    for d in (a / "Plots", b, tmp / "Data/Processed/IGVFDS0000AAAA",
              tmp / "Docs/Architecture", tmp / "Docs/Secret/20260101_000000_x"):
        d.mkdir(parents=True, exist_ok=True)
    (a / "report.md").write_text("# r\n")
    (a / "Plots/p.png").write_bytes(b"\x89PNG")
    (a / "api.token").write_text("s3cret")
    (b / "report.md").write_text("# other user\n")
    (tmp / "Data/Processed/IGVFDS0000AAAA/f.tsv").write_text("a\tb\n")
    (tmp / "Docs/Architecture/X.md").write_text("doc")
    (tmp / "Docs/Secret/20260101_000000_x/k.txt").write_text("key")

    rels = {r["rel"] for r in db.all_runs()}
    check("all_runs lists skill runs and Portal fetches",
          {str(a.relative_to(tmp)), str(b.relative_to(tmp)),
           "Data/Processed/IGVFDS0000AAAA"} <= rels)
    check("all_runs skips static docs and Docs/Secret",
          not any(r.startswith(("Docs/Architecture", "Docs/Secret"))
                  for r in rels))
    names = [p.name for p in db.run_files(a)]
    check("run_files puts report.md first", names[:1] == ["report.md"])
    check("run_files never lists secret material", "api.token" not in names)
    check("_run_of maps a nested file to its run",
          db._run_of(a / "Plots/p.png") == a)
    check("_run_of ignores files outside any run",
          db._run_of(tmp / "Docs/Architecture/X.md") is None)

    # Real history store on a temp DB: alice owns run a, bob owns run b.
    os.environ["IGVF_HISTORY_DB"] = str(tmp / "history.sqlite")
    os.environ.pop("IGVF_HISTORY_SHARED", None)
    import _history as H
    importlib.reload(H)
    H.record_session(str(tmp / "Docs/Agent/20260923_000001_alice"),
                     query="alice", meta={"artefacts": [str(a / "report.md"),
                     "/etc/passwd", "Docs/Secret/20260101_000000_x/k.txt"]},
                     owner="alice")
    H.record_session(str(tmp / "Docs/Agent/20260923_000002_bob"),
                     query="bob", meta={"artefacts": [str(b / "report.md")]},
                     owner="bob")

    scoped = {r["rel"] for r in db.runs_for_viewer("alice", H)}
    check("signed-in user sees only runs from their own sessions",
          str(a.relative_to(tmp)) in scoped
          and str(b.relative_to(tmp)) not in scoped)
    check("paths outside the workspace or under Secret never become runs",
          not any("Secret" in r or r.startswith("/") for r in scoped))
    vis = db.visibility_filter("alice", False, H)
    check("the viewer filter admits the user's own run files",
          vis is not None and vis(a / "Plots/p.png"))
    check("...and rejects another user's", not vis(b / "report.md"))
    check("...and files that belong to no run", not vis(tmp / "Docs/Architecture/X.md"))
    check("no filter for admins or a local install",
          db.visibility_filter("alice", True, H) is None
          and db.visibility_filter(None, False, H) is None)

    # ── KG sources ──
    mir = tmp / "Data/Warehouse/KG"
    (mir / "_state").mkdir(parents=True)
    (mir / "_inventory.csv").write_text(
        "type,collection,documents,bytes,skip_default\n"
        "edge,proteins_proteins,100,1000,False\n"
        "edge,genomic_elements_genes,200,2000,False\n"
        "doc,variants,300,3000,False\n"
        "edge,genes_pathways,50,500,False\n")
    for name, rows, st in (("proteins_proteins", 100, "done"),
                           ("genomic_elements_genes", 50, "in_progress"),
                           ("genes_pathways", 50, "done")):
        (mir / "_state" / f"{name}.json").write_text(json.dumps(
            {"collection": name, "rows_written": rows, "status": st}))
    check("graph_status reports a missing graph without creating it",
          ks.graph_status() == {"exists": False} and not ks.KG_DB.exists())

    ks.KG_DB.parent.mkdir(parents=True)
    con = sqlite3.connect(ks.KG_DB)
    con.executescript(
        "CREATE TABLE nodes(id TEXT); CREATE TABLE edges(source TEXT);"
        "CREATE TABLE kg_identity(alias TEXT);"
        "CREATE TABLE kg_merge_log(id TEXT, source TEXT, shard TEXT,"
        " rows_read INT, edges_added INT, resolved INT, unresolved INT,"
        " ambiguous INT, merged_at TEXT);"
        "INSERT INTO nodes VALUES('A'),('B');"
        "INSERT INTO edges VALUES('IGVF-Catalog:proteins_proteins');"
        "INSERT INTO kg_merge_log VALUES('1','IGVF-Catalog:proteins_proteins',"
        "'0.parquet',100,1,0,0,0,'2026-09-23');")
    con.commit()
    con.close()
    g = ks.graph_status()
    check("graph_status reads counts from a graph without harvest_ledger",
          g["nodes"] == 2 and g["edges"] == 1 and g["mirror_ppi_shards_merged"] == 0)
    cov = {r["collection"]: r for r in ks.coverage(ks.mirror_status(), g)}
    check("merged collection is marked merged",
          cov["proteins_proteins"]["in_graph"] == "merged")
    check("mirrored collection whose merger has not run is ready",
          cov["genes_pathways"]["in_graph"] == "ready: merge pathways")
    check("collection with no merger says so",
          cov["genomic_elements_genes"]["in_graph"] == "no merger yet")
    check("unmirrored collection is not mirrored",
          cov["variants"]["pull"] == "not mirrored"
          and cov["variants"]["mirrored_rows"] == 0)
    check("percent mirrored uses Catalog document counts",
          abs(cov["genomic_elements_genes"]["pct"] - 25.0) < 1e-9)

    ok, msg = ks.start_action("rm -rf", ["true"])
    check("start_action refuses anything outside the fixed list", not ok)
    ks.PID_FILE.parent.mkdir(parents=True, exist_ok=True)
    ks.PID_FILE.write_text(f"{os.getpid()} x.log\n")
    ok, msg = ks.start_action("merge all", ["true"])
    check("start_action refuses while another merge is running",
          not ok and "already running" in msg)

    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    print("all checks pass" if not FAIL else f"FAILED: {len(FAIL)}")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
