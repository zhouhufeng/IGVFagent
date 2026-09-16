"""Permanent history + projects store.

Every agent run already leaves a directory under ``Docs/Agent/<ts>_<slug>/``
with a transcript and a report, and every skill invocation already leaves a row
in the knowledge graph's ``analysis_log``. Neither was *searchable*: the UI's
"past results" list was a filesystem glob that opened up to 300 ``report.md``
files on every sidebar render, and there was no way to ask "what do we already
know about IGVFDS5414UFNC" without re-running the analysis that produced it.

This module adds the missing half:

* **sessions** — one durable row per agent run: the question, the answer, the
  artefacts, the accessions mentioned. Indexed, not globbed.
* **accession_index** — accession → everything ever produced about it, across
  sessions, skill runs and downloads. This is what makes a repeat question
  return in milliseconds instead of minutes.
* **projects** — a named, permanent container. Analyses are added to a project
  and stay there. A project can be renamed without breaking references,
  because references use its immutable id and every former name stays
  resolvable in ``project_names``.
* **search_fts** — one FTS5 index over questions, answers, skill runs and
  project notes, so history is searched by words the way a person remembers it.

**Nothing is ever deleted.** BEFORE DELETE triggers abort on every table that
holds history; "removing" an item from a project sets ``removed_at`` and the
row stays queryable. That is a schema-level guarantee rather than a convention,
so a future skill cannot quietly drop the record.

Deliberately a SEPARATE SQLite file from ``Data/KG/local_kg.sqlite``. The KG is
~800 MB, is bulk-merged and rebuilt by ``kg-integrate``, and is gitignored
scratch that a user may reasonably delete to reclaim disk. History must
outlive all of that, and an FTS index does not belong inside a file that gets
rewritten wholesale.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import threading
import time
from pathlib import Path
from typing import Iterable, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()

HISTORY_DIR = ROOT / "Data" / "History"
HISTORY_PATH = Path(os.environ.get("IGVF_HISTORY_DB")
                    or HISTORY_DIR / "history.sqlite")
AGENT_DIR = ROOT / "Docs" / "Agent"
KG_PATH = ROOT / "Data" / "KG" / "local_kg.sqlite"

# Where the "current project" lives between CLI invocations. A file rather
# than an env var so `igvfagent project use X` in one shell is still in force
# for the Streamlit worker and the next command.
ACTIVE_PATH = HISTORY_DIR / "active_project"

_NOW = lambda: time.strftime("%Y-%m-%dT%H:%M:%S")  # noqa: E731

# Every archive identifier worth keying results on. This started as IGVF and
# ENCODE only, which meant a question about a GEO series was re-analysed every
# time it was asked -- the index simply could not see the accession in it.
#
# Each pattern is anchored on its archive's own prefix rather than something
# generic, because a loose pattern is worse than a narrow one here: a false
# match keys one dataset's results under another dataset's name, and a wrong
# cached answer is far more expensive than a missed one.
_ACCESSION_PATTERNS = (
    r"(?:IGVF|ENC)[A-Z]{2}[0-9A-Z]{6,}",      # IGVF, ENCODE
    r"GS[EM][0-9]{3,}",                        # GEO series / sample
    r"GPL[0-9]{3,}",                           # GEO platform
    r"GDS[0-9]{3,}",                           # GEO dataset
    r"SR[RXPSZ][0-9]{5,}",                     # SRA run/experiment/project/…
    r"ERR[0-9]{5,}|ERX[0-9]{5,}|ERP[0-9]{5,}",  # ENA
    r"DRR[0-9]{5,}|DRX[0-9]{5,}",              # DDBJ
    r"PRJ(?:NA|EB|DB)[0-9]{3,}",               # BioProject
    r"SAM[NED][A-Z]?[0-9]{5,}",                # BioSample
    r"E-[A-Z]{4}-[0-9]+",                      # ArrayExpress / BioStudies
    r"phs[0-9]{6}(?:\.v[0-9]+\.p[0-9]+)?",     # dbGaP
    r"EGA[SDNCF][0-9]{6,}",                    # EGA
    r"syn[0-9]{6,}",                           # Synapse
    r"PXD[0-9]{5,}",                           # PRIDE
)
# No HCA pattern. The Human Cell Atlas keys projects by UUID and short name,
# neither of which has a prefix worth matching -- an attempt at `HCA[0-9A-Z-]+`
# matched HCASMC (human coronary artery smooth muscle cell) in this
# deployment's own history, which would have filed real results under a cell
# type mistaken for an accession. HCA datasets carry GEO/SRA/BioProject
# accessions, and those are matched above.
ACCESSION_RE = re.compile(
    r"\b(?:" + "|".join(_ACCESSION_PATTERNS) + r")\b")

ITEM_KINDS = ("session", "analysis", "artifact", "dataset", "figure",
              "paper", "note")

# ------------------------------- schema -------------------------------------

_SCHEMA = """
CREATE TABLE IF NOT EXISTS sessions (
    id          TEXT PRIMARY KEY,
    run_dir     TEXT UNIQUE,
    started_at  TEXT,
    query       TEXT,
    answer      TEXT,
    backend     TEXT,
    model       TEXT,
    iterations  INTEGER,
    tool_calls  INTEGER,
    stop_reason TEXT,
    artefacts   TEXT,
    accessions  TEXT,
    recorded_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_sessions_started ON sessions(started_at DESC);

CREATE TABLE IF NOT EXISTS projects (
    id          TEXT PRIMARY KEY,
    name        TEXT NOT NULL,
    slug        TEXT NOT NULL,
    description TEXT,
    created_at  TEXT,
    updated_at  TEXT,
    archived    INTEGER NOT NULL DEFAULT 0
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_projects_slug ON projects(slug);

-- Every name a project has ever had, so a reference written down months ago
-- still resolves after a rename.
-- A surrogate key, not (project_id, set_at): two renames inside the same
-- second would collide on that, and an upsert would then overwrite -- i.e.
-- silently lose -- the name it is supposed to preserve. Rows are only ever
-- appended here.
CREATE TABLE IF NOT EXISTS project_names (
    seq        INTEGER PRIMARY KEY AUTOINCREMENT,
    project_id TEXT NOT NULL,
    name       TEXT NOT NULL,
    slug       TEXT NOT NULL,
    set_at     TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_pnames_slug ON project_names(slug);

CREATE TABLE IF NOT EXISTS project_items (
    id         TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    kind       TEXT NOT NULL,
    ref        TEXT NOT NULL,
    title      TEXT,
    note       TEXT,
    added_at   TEXT,
    removed_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_items_project ON project_items(project_id);
CREATE INDEX IF NOT EXISTS idx_items_ref     ON project_items(ref);

CREATE TABLE IF NOT EXISTS accession_index (
    accession TEXT NOT NULL,
    kind      TEXT NOT NULL,
    ref       TEXT NOT NULL,
    title     TEXT,
    at        TEXT,
    PRIMARY KEY (accession, kind, ref)
);
CREATE INDEX IF NOT EXISTS idx_acc_at ON accession_index(accession, at DESC);

-- Append-only audit: creations, renames, adds, soft-removes, archives.
CREATE TABLE IF NOT EXISTS history_events (
    seq        INTEGER PRIMARY KEY AUTOINCREMENT,
    at         TEXT,
    kind       TEXT,
    project_id TEXT,
    detail     TEXT
);

CREATE TABLE IF NOT EXISTS ingest_ledger (
    key         TEXT PRIMARY KEY,
    kind        TEXT,
    ingested_at TEXT
);

CREATE TABLE IF NOT EXISTS project_members (
    project_id TEXT NOT NULL,
    username   TEXT NOT NULL,
    role       TEXT NOT NULL DEFAULT 'member',
    added_at   TEXT,
    PRIMARY KEY (project_id, username)
);
CREATE INDEX IF NOT EXISTS idx_members_user ON project_members(username);
"""

# The search index is DERIVED -- reindex() rebuilds it from the tables -- so
# unlike the history tables it may be dropped and recreated, which is how a
# column gets added to an FTS5 table at all.
_FTS_COLUMNS = ["kind", "ref", "title", "body", "accessions", "at", "owner"]

_FTS_SCHEMA = """
CREATE VIRTUAL TABLE IF NOT EXISTS search_fts USING fts5(
    kind  UNINDEXED,
    ref   UNINDEXED,
    title,
    body,
    accessions,
    at    UNINDEXED,
    owner UNINDEXED,
    tokenize='porter unicode61'
);
"""

# History is append-only at the schema level. A skill that tries to tidy up
# gets an exception, not a silent loss.
_TRIGGERS = """
CREATE TRIGGER IF NOT EXISTS no_delete_sessions BEFORE DELETE ON sessions
BEGIN SELECT RAISE(ABORT, 'history is append-only: sessions cannot be deleted'); END;
CREATE TRIGGER IF NOT EXISTS no_delete_projects BEFORE DELETE ON projects
BEGIN SELECT RAISE(ABORT, 'history is append-only: projects cannot be deleted (archive instead)'); END;
CREATE TRIGGER IF NOT EXISTS no_delete_items BEFORE DELETE ON project_items
BEGIN SELECT RAISE(ABORT, 'history is append-only: project items cannot be deleted (they are marked removed)'); END;
CREATE TRIGGER IF NOT EXISTS no_delete_names BEFORE DELETE ON project_names
BEGIN SELECT RAISE(ABORT, 'history is append-only: former project names are kept so old references resolve'); END;
CREATE TRIGGER IF NOT EXISTS no_delete_events BEFORE DELETE ON history_events
BEGIN SELECT RAISE(ABORT, 'history is append-only: the audit log cannot be deleted'); END;
"""


def connect() -> sqlite3.Connection:
    HISTORY_PATH.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(HISTORY_PATH), timeout=30)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA busy_timeout=30000")
    try:
        con.execute("PRAGMA journal_mode=WAL")
    except sqlite3.OperationalError:
        pass
    con.executescript(_SCHEMA)
    con.executescript(_FTS_SCHEMA)
    con.executescript(_TRIGGERS)
    _migrate(con)
    return con


# Rows written before there were user accounts. They were produced under a
# single shared password by nobody in particular, so they belong to everybody:
# treating them as private would hide the entire existing corpus from every
# user at once and throw away the recall index built over it.
LEGACY_OWNER = ""

# What an unauthenticated deployment calls the person at the keyboard. A local
# install has exactly one user and nothing to hide from them.
LOCAL_OWNER = "local"


def _migrate(con: sqlite3.Connection) -> None:
    """Additive schema upgrades. Never drops a history row.

    Called on every connect and cheap when there is nothing to do: the column
    checks are PRAGMA reads against an already-open handle.
    """
    for table, column in (("sessions", "owner"), ("projects", "owner"),
                          ("project_items", "owner")):
        have = {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
        if column not in have:
            con.execute(f"ALTER TABLE {table} ADD COLUMN {column} "
                        f"TEXT NOT NULL DEFAULT ''")

    # An FTS5 table cannot be ALTERed, so a new column means recreating it.
    # Safe precisely because nothing authoritative lives there.
    #
    # The rows are CARRIED ACROSS rather than regenerated. Regenerating would
    # mean calling reindex(), which re-reads the knowledge graph for the
    # analysis and download documents -- and if that database happened to be
    # missing, moved, or mid-rebuild, those documents would silently vanish
    # from search until someone noticed and ran a backfill. Measured on the
    # live index: 6,478 documents in, 737 out. Copying is lossless, needs no
    # second database, and runs in constant memory.
    have_cols = [r[1] for r in con.execute("PRAGMA table_info(search_fts)")]
    if have_cols and have_cols != _FTS_COLUMNS:
        con.execute("DROP TABLE IF EXISTS _fts_carry")
        con.execute("CREATE TABLE _fts_carry AS "
                    "SELECT kind, ref, title, body, accessions, at "
                    "FROM search_fts")
        con.execute("DROP TABLE search_fts")
        con.executescript(_FTS_SCHEMA)
        # owner '' = the pre-accounts corpus, visible to everyone. Every row
        # that existed before this migration was produced under the shared
        # password, so that is exactly right.
        con.execute("INSERT INTO search_fts"
                    "(kind,ref,title,body,accessions,at,owner) "
                    "SELECT kind, ref, title, body, accessions, at, '' "
                    "FROM _fts_carry")
        con.execute("DROP TABLE _fts_carry")
        con.commit()


def _digest(*parts: str) -> str:
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]


def _slugify(name: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", (name or "").strip().lower()).strip("-")
    return s or "project"


def _rel(path) -> str:
    """Project-relative where possible, so a row survives a re-deploy whose
    absolute prefix differs (/workspace here, /Users/... there)."""
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except Exception:
        return str(path)


# ------------------------------ who is acting -------------------------------
#
# The agent loop and the skills it calls are several layers below the web
# request that carries the signed-in user, and threading a username through
# every one of them would touch code that has nothing to do with identity.
# Streamlit runs each browser session's script in its own thread, so a
# thread-local is both correct here and invisible everywhere else.

_ACTOR = threading.local()


def set_actor(username: "str | None") -> None:
    """Record who the work in THIS thread belongs to. '' means nobody."""
    _ACTOR.name = (username or "").strip()


def actor() -> str:
    name = getattr(_ACTOR, "name", "")
    if name:
        return name
    # For CLI runs, where there are no threads to speak of and an operator may
    # want to attribute a batch.
    return os.environ.get("IGVF_ACTING_USER", "").strip()


def viewer() -> "str | None":
    """Who to filter for, or None when this deployment has no accounts.

    An unauthenticated deployment must behave exactly as it did before: one
    person, everything visible. Only a real username turns filtering on.
    """
    return actor() or None


# ------------------------------ visibility ----------------------------------
#
# "Private per user, shared into projects": you see your own work, plus
# anything filed into a project you own or have been added to, plus the
# pre-accounts corpus that belongs to nobody. Filtering happens in Python
# against a single predicate rather than in six hand-written SQL fragments --
# at this scale (hundreds of sessions, thousands of indexed documents) the
# query cost is irrelevant and one auditable rule beats six chances to get a
# WHERE clause subtly wrong.


def visible_projects(con, viewer: "str | None") -> "set[str] | None":
    """Project ids ``viewer`` may see, or None when nothing is hidden."""
    if viewer is None:
        return None
    rows = con.execute(
        "SELECT id FROM projects WHERE owner = ? OR owner = ? "
        "UNION SELECT project_id FROM project_members WHERE username = ?",
        (viewer, LEGACY_OWNER, viewer))
    return {r[0] for r in rows}


def _shared_refs(con, viewer: "str | None") -> "set[str]":
    """Refs made visible to ``viewer`` by being filed into a shared project."""
    pids = visible_projects(con, viewer)
    if not pids:
        return set()
    marks = ",".join("?" * len(pids))
    return {r[0] for r in con.execute(
        f"SELECT ref FROM project_items WHERE removed_at IS NULL "
        f"AND project_id IN ({marks})", tuple(pids))}


def _visible_to(con, viewer: "str | None"):
    """Build the predicate deciding what ``viewer`` may see.

    ``viewer=None`` -- an unauthenticated or local deployment -- means one
    person at the keyboard and nothing to hide, so everything passes.
    """
    if viewer is None:
        return lambda owner, ref: True
    shared = _shared_refs(con, viewer)
    allowed = {viewer, LEGACY_OWNER}

    def ok(owner: "str | None", ref: "str | None") -> bool:
        return (owner or LEGACY_OWNER) in allowed or (ref or "") in shared

    return ok


# ------------------------------ FTS helpers ---------------------------------


def _index(con, kind: str, ref: str, title: str, body: str,
           accessions: Iterable[str], at: str, owner: str = "") -> None:
    """Replace this document in the search index.

    The FTS row is derived data — it is rebuilt from the tables it summarises
    by :func:`reindex` — so unlike the history tables it is allowed to be
    deleted, and the append-only triggers do not cover it.
    """
    con.execute("DELETE FROM search_fts WHERE kind = ? AND ref = ?", (kind, ref))
    con.execute(
        "INSERT INTO search_fts(kind,ref,title,body,accessions,at,owner) "
        "VALUES(?,?,?,?,?,?,?)",
        (kind, ref, title or "", (body or "")[:200_000],
         " ".join(sorted(set(accessions or ()))), at or _NOW(), owner or ""))


def _note_accessions(con, accessions: Iterable[str], kind: str, ref: str,
                     title: str, at: str) -> int:
    n = 0
    for acc in sorted(set(a.upper() for a in accessions or () if a)):
        con.execute(
            "INSERT OR REPLACE INTO accession_index(accession,kind,ref,title,at) "
            "VALUES(?,?,?,?,?)", (acc, kind, ref, title or "", at or _NOW()))
        n += 1
    return n


def _with_output(title: str, outputs) -> str:
    """Append the first output path to a title, when there is one.

    ``outputs`` arrives as a JSON list (from analysis_log) or already parsed.
    """
    try:
        paths = json.loads(outputs) if isinstance(outputs, str) else list(outputs or [])
    except (ValueError, TypeError):
        paths = []
    for path in paths:
        text = str(path)
        if "/" in text:
            return f"{title} → {text}"
    return title


def _event(con, kind: str, project_id: str = "", detail: str = "") -> None:
    con.execute(
        "INSERT INTO history_events(at,kind,project_id,detail) VALUES(?,?,?,?)",
        (_NOW(), kind, project_id, detail))


# ------------------------------ sessions ------------------------------------


def record_session(run_dir, *, query: str, answer: str = "",
                   meta: Optional[dict] = None, owner: str = "",
                   con=None) -> dict:
    """Record one agent run. Idempotent on ``run_dir``.

    Called from ``_agent._persist_transcript`` so every run lands here as it
    happens; :func:`backfill` covers everything that ran before this existed.
    """
    meta = meta or {}
    own = con is None
    con = con or connect()
    try:
        rel = _rel(run_dir)
        sid = _digest("session", rel)
        name = Path(rel).name
        try:
            started = time.strftime(
                "%Y-%m-%dT%H:%M:%S", time.strptime(name[:15], "%Y%m%d_%H%M%S"))
        except ValueError:
            started = _NOW()
        artefacts = [str(a) for a in (meta.get("artefacts") or [])]
        accs = sorted(set(ACCESSION_RE.findall(
            f"{query}\n{answer}\n{' '.join(artefacts)}")))
        con.execute(
            "INSERT OR REPLACE INTO sessions(id,run_dir,started_at,query,answer,"
            "backend,model,iterations,tool_calls,stop_reason,artefacts,"
            "accessions,recorded_at,owner) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (sid, rel, started, query, answer,
             str(meta.get("backend", "")), str(meta.get("model", "")),
             int(meta.get("iterations", 0) or 0),
             int(meta.get("tool_calls_made", meta.get("tool_calls", 0)) or 0),
             str(meta.get("stop_reason", "")),
             json.dumps(artefacts), json.dumps(accs), _NOW(), owner or ""))
        _index(con, "session", rel, query, f"{query}\n\n{answer}", accs,
               started, owner=owner)
        _note_accessions(con, accs, "session", rel, query, started)

        # Auto-file into whichever project is active, so "everything in this
        # project" needs no bookkeeping from the user.
        pid = active_project_id(con, viewer=owner or None)
        if pid:
            _add_item(con, pid, "session", rel, title=query, note="",
                      owner=owner)
        if own:
            con.commit()
    finally:
        if own:
            con.close()
    return {"recorded": True, "id": sid, "accessions": accs}


# ------------------------------ projects ------------------------------------


def create_project(name: str, description: str = "", owner: str = "",
                   con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        slug = _slugify(name)
        row = con.execute("SELECT * FROM projects WHERE slug = ?",
                          (slug,)).fetchone()
        if row:
            return {"created": False, "id": row["id"], "name": row["name"],
                    "slug": row["slug"],
                    "note": "a project with this name already exists"}
        pid = "prj_" + _digest("project", slug, _NOW())[:10]
        now = _NOW()
        con.execute(
            "INSERT INTO projects(id,name,slug,description,created_at,"
            "updated_at,archived,owner) VALUES(?,?,?,?,?,?,0,?)",
            (pid, name.strip(), slug, description, now, now, owner or ""))
        con.execute("INSERT INTO project_names(project_id,name,slug,set_at) "
                    "VALUES(?,?,?,?)", (pid, name.strip(), slug, now))
        _index(con, "project", pid, name, f"{name}\n{description}", (), now,
               owner=owner)
        _event(con, "project_created", pid, f"{name} (owner: {owner or '-'})")
        if own:
            con.commit()
        return {"created": True, "id": pid, "name": name.strip(), "slug": slug}
    finally:
        if own:
            con.close()


def resolve_project(ref: str, viewer: "str | None" = None,
                    con=None) -> Optional[sqlite3.Row]:
    """Find a project by id, current name/slug, or any name it ever had.

    A viewer only resolves projects they can see, so a name they have no
    access to reads as "no such project" rather than as a permission error --
    which would otherwise confirm the project exists.
    """
    if not ref:
        return None
    own = con is None
    con = con or connect()
    try:
        allowed = visible_projects(con, viewer)
        ref = ref.strip()
        slug = _slugify(ref)
        for sql, arg in (
            ("SELECT * FROM projects WHERE id = ?", ref),
            ("SELECT * FROM projects WHERE slug = ?", slug),
            ("SELECT * FROM projects WHERE lower(name) = lower(?)", ref),
            ("SELECT p.* FROM projects p JOIN project_names n "
             "ON n.project_id = p.id WHERE n.slug = ? "
             "ORDER BY n.seq DESC LIMIT 1", slug),
        ):
            row = con.execute(sql, (arg,)).fetchone()
            if row and (allowed is None or row["id"] in allowed):
                return row
        return None
    finally:
        if own:
            con.close()


def _owner_only(row, viewer: "str | None") -> "dict | None":
    """Guard the operations only a project's owner may perform.

    Sharing grants VISIBILITY, not control. Without this, adding a collaborator
    to a project also handed them the ability to rename it, rewrite its
    description, or archive it out from under the owner. Adding and removing
    items stays open to members -- that is the collaboration the sharing was
    for; the project's identity and its membership list are not.
    """
    if viewer is None or row["owner"] in ("", viewer):
        return None
    return {"ok": False,
            "error": f"only {row['owner']} can do that to {row['name']!r} — "
                     f"you have access to it, but you do not own it"}


def rename_project(ref: str, new_name: str, viewer: "str | None" = None,
                   con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref!r}"}
        denied = _owner_only(row, viewer)
        if denied:
            return denied
        new_slug = _slugify(new_name)
        clash = con.execute(
            "SELECT id FROM projects WHERE slug = ? AND id != ?",
            (new_slug, row["id"])).fetchone()
        if clash:
            return {"ok": False,
                    "error": f"another project already uses the name {new_name!r}"}
        now = _NOW()
        con.execute("UPDATE projects SET name = ?, slug = ?, updated_at = ? "
                    "WHERE id = ?", (new_name.strip(), new_slug, now, row["id"]))
        con.execute("INSERT INTO project_names"
                    "(project_id,name,slug,set_at) VALUES(?,?,?,?)",
                    (row["id"], new_name.strip(), new_slug, now))
        _index(con, "project", row["id"], new_name,
               f"{new_name}\n{row['description'] or ''}", (), now)
        _event(con, "project_renamed", row["id"],
               f"{row['name']} -> {new_name.strip()}")
        if own:
            con.commit()
        return {"ok": True, "id": row["id"], "was": row["name"],
                "name": new_name.strip(), "slug": new_slug}
    finally:
        if own:
            con.close()


def set_description(ref: str, description: str, viewer: "str | None" = None,
                    con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref!r}"}
        denied = _owner_only(row, viewer)
        if denied:
            return denied
        con.execute("UPDATE projects SET description = ?, updated_at = ? "
                    "WHERE id = ?", (description, _NOW(), row["id"]))
        _index(con, "project", row["id"], row["name"],
               f"{row['name']}\n{description}", (), _NOW())
        _event(con, "project_described", row["id"], description[:200])
        if own:
            con.commit()
        return {"ok": True, "id": row["id"]}
    finally:
        if own:
            con.close()


def archive_project(ref: str, archived: bool = True,
                    viewer: "str | None" = None, con=None) -> dict:
    """Hide a project from the default listing. The rows stay — archiving is
    the only "delete" this store has, and it is reversible."""
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref!r}"}
        denied = _owner_only(row, viewer)
        if denied:
            return denied
        con.execute("UPDATE projects SET archived = ?, updated_at = ? WHERE id = ?",
                    (1 if archived else 0, _NOW(), row["id"]))
        _event(con, "project_archived" if archived else "project_unarchived",
               row["id"], row["name"])
        if own:
            con.commit()
        return {"ok": True, "id": row["id"], "archived": bool(archived)}
    finally:
        if own:
            con.close()


def list_projects(include_archived: bool = False,
                  viewer: "str | None" = None, con=None) -> "list[dict]":
    own = con is None
    con = con or connect()
    try:
        sql = ("SELECT p.*, "
               "(SELECT COUNT(*) FROM project_items i "
               " WHERE i.project_id = p.id AND i.removed_at IS NULL) AS n_items "
               "FROM projects p")
        if not include_archived:
            sql += " WHERE p.archived = 0"
        sql += " ORDER BY p.updated_at DESC"
        allowed = visible_projects(con, viewer)
        return [dict(r) for r in con.execute(sql)
                if allowed is None or r["id"] in allowed]
    finally:
        if own:
            con.close()


# ------------------------------- items --------------------------------------


def _add_item(con, project_id: str, kind: str, ref: str, *, title: str = "",
              note: str = "", owner: str = "") -> str:
    iid = _digest("item", project_id, kind, ref)
    now = _NOW()
    existing = con.execute("SELECT id FROM project_items WHERE id = ?",
                           (iid,)).fetchone()
    if existing:
        # Re-adding un-removes, and refreshes the title/note.
        con.execute("UPDATE project_items SET removed_at = NULL, title = ?, "
                    "note = ? WHERE id = ?", (title, note, iid))
    else:
        con.execute("INSERT INTO project_items(id,project_id,kind,ref,title,"
                    "note,added_at,removed_at,owner) "
                    "VALUES(?,?,?,?,?,?,?,NULL,?)",
                    (iid, project_id, kind, ref, title, note, now,
                     owner or ""))
    con.execute("UPDATE projects SET updated_at = ? WHERE id = ?",
                (now, project_id))
    accs = sorted(set(ACCESSION_RE.findall(f"{ref} {title} {note}")))
    _index(con, "item", f"{project_id}/{kind}/{ref}", title or ref,
           f"{ref}\n{title}\n{note}", accs, now, owner=owner)
    _event(con, "item_added", project_id, f"{kind}:{ref}")
    return iid


def add_item(ref_project: str, kind: str, ref: str, *, title: str = "",
             note: str = "", owner: str = "", viewer: "str | None" = None,
             con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref_project!r}"}
        if kind not in ITEM_KINDS:
            return {"ok": False,
                    "error": f"kind must be one of {', '.join(ITEM_KINDS)}"}
        iid = _add_item(con, row["id"], kind, ref, title=title, note=note,
                        owner=owner)
        if own:
            con.commit()
        return {"ok": True, "id": iid, "project": row["name"]}
    finally:
        if own:
            con.close()


def remove_item(ref_project: str, kind: str, ref: str,
                viewer: "str | None" = None, con=None) -> dict:
    """Soft-remove: the row stays and stays searchable, it just stops counting
    as part of the project's current contents."""
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref_project!r}"}
        iid = _digest("item", row["id"], kind, ref)
        cur = con.execute("UPDATE project_items SET removed_at = ? "
                          "WHERE id = ? AND removed_at IS NULL", (_NOW(), iid))
        _event(con, "item_removed", row["id"], f"{kind}:{ref}")
        if own:
            con.commit()
        return {"ok": True, "changed": cur.rowcount}
    finally:
        if own:
            con.close()


def project_items(ref_project: str, include_removed: bool = False,
                  viewer: "str | None" = None, con=None) -> "list[dict]":
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return []
        sql = "SELECT * FROM project_items WHERE project_id = ?"
        if not include_removed:
            sql += " AND removed_at IS NULL"
        sql += " ORDER BY added_at DESC"
        return [dict(r) for r in con.execute(sql, (row["id"],))]
    finally:
        if own:
            con.close()


# ------------------------------ sharing -------------------------------------


def share_project(ref_project: str, username: str, *, role: str = "member",
                  viewer: "str | None" = None, con=None) -> dict:
    """Give another user access to a project, and to everything in it.

    This is the ONLY way work crosses between users: sessions are private to
    whoever ran them, and filing one into a shared project is the deliberate
    act that publishes it to that project's members.
    """
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref_project!r}"}
        denied = _owner_only(row, viewer)
        if denied:
            return denied
        username = username.strip()
        if not username:
            return {"ok": False, "error": "no username given"}
        con.execute("INSERT OR REPLACE INTO project_members"
                    "(project_id,username,role,added_at) VALUES(?,?,?,?)",
                    (row["id"], username, role, _NOW()))
        _event(con, "project_shared", row["id"], f"{username} ({role})")
        if own:
            con.commit()
        return {"ok": True, "project": row["name"], "username": username}
    finally:
        if own:
            con.close()


def unshare_project(ref_project: str, username: str,
                    viewer: "str | None" = None, con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return {"ok": False, "error": f"no project matches {ref_project!r}"}
        denied = _owner_only(row, viewer)
        if denied:
            return denied
        cur = con.execute("DELETE FROM project_members WHERE project_id = ? "
                          "AND username = ?", (row["id"], username.strip()))
        _event(con, "project_unshared", row["id"], username.strip())
        if own:
            con.commit()
        return {"ok": True, "changed": cur.rowcount}
    finally:
        if own:
            con.close()


def project_members(ref_project: str, viewer: "str | None" = None,
                    con=None) -> "list[dict]":
    own = con is None
    con = con or connect()
    try:
        row = resolve_project(ref_project, viewer=viewer, con=con)
        if not row:
            return []
        out = [{"username": row["owner"] or "(pre-accounts)", "role": "owner",
                "added_at": row["created_at"]}]
        out += [dict(r) for r in con.execute(
            "SELECT username, role, added_at FROM project_members "
            "WHERE project_id = ? ORDER BY added_at", (row["id"],))]
        return out
    finally:
        if own:
            con.close()


# --------------------------- active project ---------------------------------


def _active_path(viewer: "str | None") -> Path:
    """Where this user's active-project marker lives.

    Per user, not per deployment: with accounts, one person picking a project
    must not silently redirect everyone else's runs into it.
    """
    if not viewer:
        return ACTIVE_PATH
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", viewer)[:64]
    return HISTORY_DIR / "active" / f"{safe}"


def set_active(ref: Optional[str], viewer: "str | None" = None) -> dict:
    path = _active_path(viewer)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not ref:
        if path.exists():
            path.unlink()
        return {"ok": True, "active": None}
    row = resolve_project(ref, viewer=viewer)
    if not row:
        return {"ok": False, "error": f"no project matches {ref!r}"}
    path.write_text(row["id"])
    return {"ok": True, "active": row["name"], "id": row["id"]}


def active_project_id(con=None, viewer: "str | None" = None) -> Optional[str]:
    """Id of the project new runs are filed into, or None.

    ``IGVF_PROJECT`` wins over the on-disk marker, so a single command or a
    container can be scoped without disturbing the persistent choice.
    """
    env = os.environ.get("IGVF_PROJECT", "").strip()
    if env:
        row = resolve_project(env, viewer=viewer, con=con)
        return row["id"] if row else None
    try:
        pid = _active_path(viewer).read_text().strip()
    except OSError:
        return None
    if not pid:
        return None
    row = resolve_project(pid, viewer=viewer, con=con)
    return row["id"] if row else None


def active_project(con=None, viewer: "str | None" = None) -> Optional[dict]:
    pid = active_project_id(con=con, viewer=viewer)
    if not pid:
        return None
    row = resolve_project(pid, viewer=viewer, con=con)
    return dict(row) if row else None


# ------------------------------- search -------------------------------------


def _fts_query(text: str) -> str:
    """Turn free user text into a safe FTS5 MATCH expression.

    Every token is quoted, so accessions (which contain no operators but do
    contain digits) and stray punctuation such as a trailing ``?`` cannot be
    read as FTS syntax and raise. A trailing ``*`` is preserved as a prefix
    search because that is the one operator worth exposing.
    """
    toks = re.findall(r"[A-Za-z0-9_:\-\.]+\*?", text or "")
    parts = []
    for t in toks:
        star = t.endswith("*")
        core = t.rstrip("*").replace('"', "")
        if core:
            parts.append(f'"{core}"' + ("*" if star else ""))
    return " AND ".join(parts)


def search(query: str, *, limit: int = 20, kind: str = "",
           project: str = "", viewer: "str | None" = None,
           con=None) -> "list[dict]":
    """Full-text search over questions, answers, skill runs and project notes."""
    own = con is None
    con = con or connect()
    try:
        can_see = _visible_to(con, viewer)
        # Overfetch, because rows the viewer may not see are dropped below and
        # a LIMIT applied before that filter would silently return short.
        want = limit if (viewer is None and not project) else limit * 10
        expr = _fts_query(query)
        if not expr:
            sql = ("SELECT kind, ref, title, at, owner, '' AS snip "
                   "FROM search_fts ORDER BY at DESC LIMIT ?")
            rows = con.execute(sql, (want,)).fetchall()
        else:
            sql = ("SELECT kind, ref, title, at, owner, "
                   "snippet(search_fts, 3, '«', '»', ' … ', 18) AS snip "
                   "FROM search_fts WHERE search_fts MATCH ? ")
            args: list = [expr]
            if kind:
                sql += "AND kind = ? "
                args.append(kind)
            sql += ("ORDER BY bm25(search_fts, 0.0, 0.0, 4.0, 1.0, 8.0, 0.0, 0.0) "
                    "LIMIT ?")
            args.append(want)
            try:
                rows = con.execute(sql, args).fetchall()
            except sqlite3.OperationalError:
                return []
        out = [dict(r) for r in rows if can_see(r["owner"], r["ref"])]
        if project:
            prow = resolve_project(project, viewer=viewer, con=con)
            if not prow:
                return []
            refs = {r["ref"] for r in con.execute(
                "SELECT ref FROM project_items WHERE project_id = ? "
                "AND removed_at IS NULL", (prow["id"],))}
            out = [r for r in out if r["ref"] in refs]
        return out[:limit]
    finally:
        if own:
            con.close()


def by_accession(accession: str, limit: int = 50,
                 viewer: "str | None" = None, con=None) -> "list[dict]":
    """Everything ever produced about one accession, newest first.

    The index itself carries no owner -- it points at sessions and skill runs
    -- so visibility is resolved against those: a session row is shown only if
    its session is visible. Skill-run and download rows come from the shared
    knowledge graph and are not per-user work, so they are shown to everyone.
    """
    own = con is None
    con = con or connect()
    try:
        rows = con.execute(
            "SELECT * FROM accession_index WHERE accession = ? "
            "ORDER BY at DESC LIMIT ?",
            (accession.strip().upper(), limit * 10 if viewer else limit))
        out = [dict(r) for r in rows]
        if viewer is not None:
            can_see = _visible_to(con, viewer)
            owners = {r[0]: r[1] for r in con.execute(
                "SELECT run_dir, owner FROM sessions")}
            out = [r for r in out
                   if r["kind"] != "session"
                   or can_see(owners.get(r["ref"], LEGACY_OWNER), r["ref"])]
        return out[:limit]
    finally:
        if own:
            con.close()


def session(run_dir: str, viewer: "str | None" = None,
            con=None) -> Optional[dict]:
    own = con is None
    con = con or connect()
    try:
        rel = _rel(run_dir)
        row = con.execute("SELECT * FROM sessions WHERE run_dir = ? OR id = ?",
                          (rel, rel)).fetchone()
        if not row:
            return None
        if not _visible_to(con, viewer)(row["owner"], row["run_dir"]):
            return None
        return dict(row)
    finally:
        if own:
            con.close()


def recent_sessions(limit: int = 50, viewer: "str | None" = None,
                    con=None) -> "list[dict]":
    own = con is None
    con = con or connect()
    try:
        rows = con.execute(
            "SELECT run_dir, started_at, query, accessions, stop_reason, owner "
            "FROM sessions ORDER BY started_at DESC LIMIT ?",
            (limit * 10 if viewer else limit,))
        can_see = _visible_to(con, viewer)
        return [dict(r) for r in rows
                if can_see(r["owner"], r["run_dir"])][:limit]
    finally:
        if own:
            con.close()


# ------------------------------ backfill ------------------------------------


def _parse_report(path: Path) -> dict:
    """Pull the query and the final answer back out of a run's report.md.

    transcript.json is authoritative when it exists; report.md is the fallback
    for runs whose transcript was never written or has since been pruned.
    """
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return {}
    query, meta = "", {}
    for line in text.splitlines():
        if line.startswith("**Query:**"):
            query = line.split("**Query:**", 1)[1].strip()
        elif line.startswith("- Backend:"):
            m = re.search(r"Backend: `([^`]*)`.*Model: `([^`]*)`", line)
            if m:
                meta["backend"], meta["model"] = m.group(1), m.group(2)
        elif line.startswith("- Iterations:"):
            m = re.search(r"Iterations: (\d+).*Tool calls: (\d+)", line)
            if m:
                meta["iterations"] = int(m.group(1))
                meta["tool_calls_made"] = int(m.group(2))
        elif line.startswith("- Stop reason:"):
            m = re.search(r"`([^`]*)`", line)
            if m:
                meta["stop_reason"] = m.group(1)
    answer = ""
    if "## Final answer" in text:
        answer = text.split("## Final answer", 1)[1]
        if "\n## " in answer:
            answer = answer.split("\n## ", 1)[0]
    return {"query": query, "answer": answer.strip(), "meta": meta}


def backfill_sessions(limit: int = 100_000, con=None) -> dict:
    """Index every past run under Docs/Agent/ that is not indexed yet."""
    own = con is None
    con = con or connect()
    added = skipped = 0
    try:
        try:
            dirs = sorted((d for d in AGENT_DIR.iterdir() if d.is_dir()),
                          key=lambda d: d.name)
        except OSError:
            dirs = []
        known = {r[0] for r in con.execute("SELECT run_dir FROM sessions")}
        for d in dirs[:limit]:
            rel = _rel(d)
            if rel in known:
                skipped += 1
                continue
            query = answer = ""
            meta: dict = {}
            tj = d / "transcript.json"
            if tj.exists():
                try:
                    blob = json.loads(tj.read_text(errors="replace"))
                    query = blob.get("query", "")
                    answer = blob.get("final_answer", "") or ""
                    meta = {k: v for k, v in blob.items()
                            if k not in ("query", "transcript", "final_answer")}
                except (OSError, ValueError):
                    pass
            if not query:
                parsed = _parse_report(d / "report.md")
                query = parsed.get("query", "") or " ".join(
                    d.name.split("_")[2:-1]).replace("-", " ")
                answer = answer or parsed.get("answer", "")
                meta = meta or parsed.get("meta", {})
            record_session(d, query=query, answer=answer, meta=meta, con=con)
            added += 1
        con.commit()
    finally:
        if own:
            con.close()
    return {"sessions_added": added, "already_indexed": skipped}


def backfill_analyses(limit: int = 200_000, con=None) -> dict:
    """Index the knowledge graph's analysis_log and download_log.

    Those tables record every direct-CLI and tool-call run — 6k rows of work
    that never went through the agent loop and so has no Docs/Agent directory.
    Read-only against the KG; nothing is written back to it.
    """
    own = con is None
    con = con or connect()
    added = 0
    try:
        if not KG_PATH.exists():
            return {"analyses_added": 0, "note": "no local KG yet"}
        kg = sqlite3.connect(f"file:{KG_PATH}?mode=ro", uri=True, timeout=30)
        kg.row_factory = sqlite3.Row
        try:
            seen = {r[0] for r in con.execute(
                "SELECT key FROM ingest_ledger WHERE kind = 'analysis'")}
            rows = kg.execute(
                "SELECT id,skill,subcommand,label,inputs,outputs,recorded_at "
                "FROM analysis_log ORDER BY recorded_at DESC LIMIT ?", (limit,))
            for r in rows:
                if r["id"] in seen:
                    continue
                title = " ".join(x for x in (r["skill"], r["subcommand"],
                                             r["label"]) if x).strip()
                # A skill run's id is an opaque digest, so the title carries
                # where its output landed -- that is the part a person can
                # actually open.
                title = _with_output(title, r["outputs"])
                body = f"{title}\n{r['inputs'] or ''}\n{r['outputs'] or ''}"
                accs = sorted(set(ACCESSION_RE.findall(body)))
                at = r["recorded_at"] or _NOW()
                _index(con, "analysis", r["id"], title, body, accs, at)
                _note_accessions(con, accs, "analysis", r["id"], title, at)
                con.execute("INSERT OR REPLACE INTO ingest_ledger"
                            "(key,kind,ingested_at) VALUES(?,?,?)",
                            (r["id"], "analysis", _NOW()))
                added += 1
            for r in kg.execute(
                    "SELECT id,source,accession,url,path,recorded_at "
                    "FROM download_log"):
                key = "dl:" + r["id"]
                if key in seen:
                    continue
                title = f"{r['source']} {r['accession'] or ''}".strip()
                body = f"{title}\n{r['url'] or ''}\n{r['path'] or ''}"
                accs = sorted(set(ACCESSION_RE.findall(body)))
                at = r["recorded_at"] or _NOW()
                _index(con, "download", r["id"], title, body, accs, at)
                _note_accessions(con, accs, "download", r["id"], title, at)
                con.execute("INSERT OR REPLACE INTO ingest_ledger"
                            "(key,kind,ingested_at) VALUES(?,?,?)",
                            (key, "analysis", _NOW()))
                added += 1
        finally:
            kg.close()
        con.commit()
    finally:
        if own:
            con.close()
    return {"analyses_added": added}


def reindex(con=None) -> dict:
    """Rebuild the FTS index from the history tables.

    The index is derived; the tables are the record. If the index is ever
    corrupted or a schema change adds a column, this regenerates it without
    touching a single history row.
    """
    own = con is None
    con = con or connect()
    try:
        con.execute("DELETE FROM search_fts")
        n = 0
        for r in con.execute("SELECT * FROM sessions"):
            accs = json.loads(r["accessions"] or "[]")
            _index(con, "session", r["run_dir"], r["query"],
                   f"{r['query']}\n\n{r['answer']}", accs, r["started_at"],
                   owner=r["owner"])
            n += 1
        for r in con.execute("SELECT * FROM projects"):
            _index(con, "project", r["id"], r["name"],
                   f"{r['name']}\n{r['description'] or ''}", (),
                   r["updated_at"], owner=r["owner"])
            n += 1
        for r in con.execute("SELECT * FROM project_items WHERE removed_at IS NULL"):
            accs = sorted(set(ACCESSION_RE.findall(
                f"{r['ref']} {r['title'] or ''} {r['note'] or ''}")))
            _index(con, "item", f"{r['project_id']}/{r['kind']}/{r['ref']}",
                   r["title"] or r["ref"],
                   f"{r['ref']}\n{r['title'] or ''}\n{r['note'] or ''}",
                   accs, r["added_at"], owner=r["owner"])
            n += 1
        # analysis/download rows live only in the KG, so re-ingest them.
        con.execute("DELETE FROM ingest_ledger WHERE kind = 'analysis'")
        con.commit()
        extra = backfill_analyses(con=con)
        con.commit()
        return {"reindexed": n, **extra}
    finally:
        if own:
            con.close()


def stats(con=None) -> dict:
    own = con is None
    con = con or connect()
    try:
        def one(sql, *a):
            try:
                return con.execute(sql, a).fetchone()[0]
            except sqlite3.Error:
                return 0
        return {
            "db": str(HISTORY_PATH),
            "size_mb": round(HISTORY_PATH.stat().st_size / 1e6, 2)
                       if HISTORY_PATH.exists() else 0.0,
            "sessions": one("SELECT COUNT(*) FROM sessions"),
            "projects": one("SELECT COUNT(*) FROM projects WHERE archived = 0"),
            "projects_archived": one("SELECT COUNT(*) FROM projects WHERE archived = 1"),
            "project_items": one("SELECT COUNT(*) FROM project_items WHERE removed_at IS NULL"),
            "indexed_documents": one("SELECT COUNT(*) FROM search_fts"),
            "accessions": one("SELECT COUNT(DISTINCT accession) FROM accession_index"),
            "accession_rows": one("SELECT COUNT(*) FROM accession_index"),
            "events": one("SELECT COUNT(*) FROM history_events"),
            "active_project": (active_project(con=con) or {}).get("name"),
        }
    finally:
        if own:
            con.close()
