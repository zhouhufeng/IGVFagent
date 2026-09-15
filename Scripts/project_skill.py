"""Projects and permanent, searchable history.

    igvfagent project create "Spatial ATAC-Hi-C reproduction" --description "..."
    igvfagent project use   "Spatial ATAC-Hi-C reproduction"
    igvfagent project add   --kind dataset --ref IGVFDS5414UFNC --title "primary set"
    igvfagent project items
    igvfagent project rename "Spatial ATAC-Hi-C reproduction" --to "Wang 2026 repro"
    igvfagent project list

    igvfagent project search "spatial ATAC Hi-C"
    igvfagent project recall IGVFDS5414UFNC
    igvfagent project show Docs/Agent/20260914_150552_...

    igvfagent project backfill      # index everything that ran before this existed
    igvfagent project stats

Two things this gives that the filesystem did not.

**Recall.** ``recall <accession>`` answers "what have we already produced about
this dataset" from an index instead of from a re-run: every past agent session,
skill invocation and download that ever mentioned it, newest first, with the
run directory to open. A repeat question costs a lookup rather than an
analysis.

**Projects.** A named container that outlives the session. Anything filed into
a project stays filed; renaming the project does not break references, because
items point at its immutable id and every former name stays resolvable.

Nothing here can delete history: the store's tables carry BEFORE DELETE
triggers, and ``remove`` marks an item removed rather than dropping the row.
Archiving hides a project from the default listing and is reversible.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    from . import _history as H
except ImportError:                                      # direct execution
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _history as H  # type: ignore


def _out(obj, as_json: bool) -> None:
    if as_json:
        print(json.dumps(obj, indent=2, default=str))


def _fail(msg: str) -> int:
    print(f"ERROR: {msg}", file=sys.stderr)
    return 1


# ------------------------------- commands -----------------------------------


def cmd_create(a) -> int:
    res = H.create_project(a.name, a.description or "")
    if res.get("created"):
        print(f"Created project {res['name']!r}  (id {res['id']})")
        if a.use:
            H.set_active(res["id"])
            print(f"Now the active project — new runs are filed into it.")
    else:
        print(f"Project {res['name']!r} already exists (id {res['id']}).")
    _out(res, a.json)
    return 0


def cmd_list(a) -> int:
    rows = H.list_projects(include_archived=a.all)
    if not rows:
        print("No projects yet. Create one:\n"
              "  igvfagent project create \"My study\" --use")
        return 0
    active = H.active_project() or {}
    print(f"{'':2}{'name':40} {'items':>6}  {'updated':19}  id")
    for r in rows:
        mark = "*" if r["id"] == active.get("id") else (
            "~" if r["archived"] else " ")
        print(f"{mark} {r['name'][:40]:40} {r['n_items']:>6}  "
              f"{(r['updated_at'] or '')[:19]:19}  {r['id']}")
    print("\n* = active   ~ = archived")
    _out(rows, a.json)
    return 0


def cmd_use(a) -> int:
    res = H.set_active(a.name if a.name != "-" else None)
    if not res.get("ok"):
        return _fail(res["error"])
    if res.get("active"):
        print(f"Active project: {res['active']}")
    else:
        print("No active project — runs are recorded in history but unfiled.")
    _out(res, a.json)
    return 0


def cmd_rename(a) -> int:
    res = H.rename_project(a.name, a.to)
    if not res.get("ok"):
        return _fail(res["error"])
    print(f"Renamed {res['was']!r} → {res['name']!r}  (id {res['id']} unchanged; "
          f"the old name still resolves)")
    _out(res, a.json)
    return 0


def cmd_describe(a) -> int:
    res = H.set_description(a.name, a.description)
    if not res.get("ok"):
        return _fail(res["error"])
    print("Description updated.")
    _out(res, a.json)
    return 0


def cmd_archive(a) -> int:
    res = H.archive_project(a.name, archived=not a.unarchive)
    if not res.get("ok"):
        return _fail(res["error"])
    print(f"{'Archived' if res['archived'] else 'Restored'} — nothing was "
          f"deleted; contents stay searchable.")
    _out(res, a.json)
    return 0


def _target_project(a) -> "str | None":
    if getattr(a, "project", None):
        return a.project
    act = H.active_project()
    return act["id"] if act else None


def cmd_add(a) -> int:
    proj = _target_project(a)
    if not proj:
        return _fail("no project given and none active — "
                     "`igvfagent project use <name>` or pass --project")
    res = H.add_item(proj, a.kind, a.ref, title=a.title or "", note=a.note or "")
    if not res.get("ok"):
        return _fail(res["error"])
    print(f"Added {a.kind} {a.ref} to {res['project']!r}.")
    _out(res, a.json)
    return 0


def cmd_remove(a) -> int:
    proj = _target_project(a)
    if not proj:
        return _fail("no project given and none active")
    res = H.remove_item(proj, a.kind, a.ref)
    if not res.get("ok"):
        return _fail(res["error"])
    print("Marked removed. The row is kept and stays searchable.")
    _out(res, a.json)
    return 0


def cmd_items(a) -> int:
    proj = _target_project(a)
    if not proj:
        return _fail("no project given and none active")
    rows = H.project_items(proj, include_removed=a.all)
    if not rows:
        print("No items in this project yet.")
        return 0
    for r in rows:
        mark = "~" if r["removed_at"] else " "
        title = (r["title"] or "")[:60]
        print(f"{mark} {(r['added_at'] or '')[:16]:16}  {r['kind']:9} "
              f"{r['ref'][:48]:48}  {title}")
    if a.all:
        print("\n~ = removed from the project (kept in history)")
    _out(rows, a.json)
    return 0


def cmd_search(a) -> int:
    rows = H.search(a.query, limit=a.limit, kind=a.kind or "",
                    project=a.project or "")
    if not rows:
        print("Nothing in history matches that.")
        return 0
    for r in rows:
        print(f"[{r['kind']:8}] {(r['at'] or '')[:16]:16}  "
              f"{(r['title'] or '')[:70]}")
        print(f"           {r['ref']}")
        if r.get("snip"):
            print(f"           {' '.join(r['snip'].split())[:180]}")
    print(f"\n{len(rows)} result(s).")
    _out(rows, a.json)
    return 0


def cmd_recall(a) -> int:
    acc = a.accession.strip().upper()
    rows = H.by_accession(acc, limit=a.limit)
    if not rows:
        print(f"Nothing recorded yet for {acc}. "
              f"If work predates the history store, run "
              f"`igvfagent project backfill` first.")
        return 0
    print(f"{acc} — {len(rows)} recorded result(s), newest first:\n")
    for r in rows:
        print(f"[{r['kind']:8}] {(r['at'] or '')[:16]:16}  "
              f"{(r['title'] or '')[:70]}")
        print(f"           {r['ref']}")
    print("\nOpen an agent session with:  igvfagent project show <ref>")
    _out(rows, a.json)
    return 0


def cmd_show(a) -> int:
    row = H.session(a.ref)
    if not row:
        print(f"No recorded session at {a.ref!r}.")
        return 1
    print(f"# {row['started_at']}   ({row['backend']} / {row['model']})")
    print(f"Query: {row['query']}\n")
    print(f"Iterations {row['iterations']} · tool calls {row['tool_calls']} · "
          f"stop {row['stop_reason']}")
    accs = json.loads(row["accessions"] or "[]")
    if accs:
        print("Accessions: " + ", ".join(accs))
    arts = json.loads(row["artefacts"] or "[]")
    if arts:
        print("Artefacts:")
        for x in arts[:40]:
            print(f"  {x}")
    if not a.brief:
        print("\n--- answer ---\n")
        print(row["answer"] or "(none recorded)")
    _out(dict(row), a.json)
    return 0


def cmd_recent(a) -> int:
    rows = H.recent_sessions(limit=a.limit)
    if not rows:
        print("No sessions recorded yet — try `igvfagent project backfill`.")
        return 0
    for r in rows:
        accs = ", ".join(json.loads(r["accessions"] or "[]")[:3])
        print(f"{(r['started_at'] or '')[:16]:16}  {(r['query'] or '')[:64]}")
        print(f"                  {r['run_dir']}" + (f"   [{accs}]" if accs else ""))
    _out(rows, a.json)
    return 0


def cmd_backfill(a) -> int:
    print("Indexing past agent runs under Docs/Agent/ …")
    s = H.backfill_sessions()
    print(f"  sessions: {s['sessions_added']} added, "
          f"{s['already_indexed']} already indexed")
    print("Indexing skill runs and downloads from the local KG …")
    t = H.backfill_analyses()
    print(f"  analyses/downloads: {t.get('analyses_added', 0)} added"
          + (f"  ({t['note']})" if t.get("note") else ""))
    st = H.stats()
    print(f"\nHistory now holds {st['sessions']} sessions, "
          f"{st['indexed_documents']} searchable documents and "
          f"{st['accessions']} distinct accessions.")
    _out({**s, **t, **st}, a.json)
    return 0


def cmd_reindex(a) -> int:
    res = H.reindex()
    print(f"Rebuilt the search index: {res['reindexed']} history rows "
          f"+ {res.get('analyses_added', 0)} KG rows. No history was altered.")
    _out(res, a.json)
    return 0


def cmd_stats(a) -> int:
    st = H.stats()
    width = max(len(k) for k in st)
    for k, v in st.items():
        print(f"{k:<{width}}  {v}")
    _out(st, a.json)
    return 0


# --------------------------------- CLI --------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent project",
        description="Named projects plus a permanent, searchable history of "
                    "every analysis. Nothing recorded here is ever deleted.")
    p.add_argument("--json", action="store_true",
                   help="Also emit the result as JSON.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("create", help="Create a project.")
    s.add_argument("name")
    s.add_argument("--description", default="")
    s.add_argument("--use", action="store_true",
                   help="Also make it the active project.")
    s.set_defaults(func=cmd_create)

    s = sub.add_parser("list", help="List projects.")
    s.add_argument("--all", action="store_true", help="Include archived.")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("use", help="Set the active project ('-' to clear).")
    s.add_argument("name")
    s.set_defaults(func=cmd_use)

    s = sub.add_parser("rename", help="Rename a project; references survive.")
    s.add_argument("name")
    s.add_argument("--to", required=True)
    s.set_defaults(func=cmd_rename)

    s = sub.add_parser("describe", help="Set a project's description.")
    s.add_argument("name")
    s.add_argument("description")
    s.set_defaults(func=cmd_describe)

    s = sub.add_parser("archive", help="Hide a project (reversible).")
    s.add_argument("name")
    s.add_argument("--unarchive", action="store_true")
    s.set_defaults(func=cmd_archive)

    s = sub.add_parser("add", help="File something into a project.")
    s.add_argument("--project", default="", help="Defaults to the active one.")
    s.add_argument("--kind", default="note", choices=list(H.ITEM_KINDS))
    s.add_argument("--ref", required=True,
                   help="Run directory, accession, artefact path, or free text.")
    s.add_argument("--title", default="")
    s.add_argument("--note", default="")
    s.set_defaults(func=cmd_add)

    s = sub.add_parser("remove", help="Mark an item removed (row is kept).")
    s.add_argument("--project", default="")
    s.add_argument("--kind", required=True, choices=list(H.ITEM_KINDS))
    s.add_argument("--ref", required=True)
    s.set_defaults(func=cmd_remove)

    s = sub.add_parser("items", help="List a project's contents.")
    s.add_argument("--project", default="")
    s.add_argument("--all", action="store_true", help="Include removed items.")
    s.set_defaults(func=cmd_items)

    s = sub.add_parser("search", help="Full-text search across all history.")
    s.add_argument("query")
    s.add_argument("--limit", type=int, default=20)
    s.add_argument("--kind", default="",
                   help="session | analysis | download | project | item")
    s.add_argument("--project", default="", help="Restrict to one project.")
    s.set_defaults(func=cmd_search)

    s = sub.add_parser(
        "recall", help="Everything ever produced about one accession.")
    s.add_argument("accession")
    s.add_argument("--limit", type=int, default=50)
    s.set_defaults(func=cmd_recall)

    s = sub.add_parser("show", help="Replay a recorded agent session.")
    s.add_argument("ref", help="Run directory from search/recall.")
    s.add_argument("--brief", action="store_true", help="Omit the answer text.")
    s.set_defaults(func=cmd_show)

    s = sub.add_parser("recent", help="Most recent sessions.")
    s.add_argument("--limit", type=int, default=20)
    s.set_defaults(func=cmd_recent)

    s = sub.add_parser(
        "backfill",
        help="Index work that predates this store (safe to re-run).")
    s.set_defaults(func=cmd_backfill)

    s = sub.add_parser("reindex", help="Rebuild the search index from the rows.")
    s.set_defaults(func=cmd_reindex)

    s = sub.add_parser("stats", help="What the history store holds.")
    s.set_defaults(func=cmd_stats)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
