"""The working directory of a local IGVF Agent.

One folder is the workspace: the agent reads the user's files from it and writes its results (``Docs/``,
``Data/``) into it. It is ``IGVF_PROJECT_ROOT`` when that is set and the IGVF Agent folder otherwise.
``igvfagent ask --workdir DIR`` and ``igvfagent ui --workdir DIR`` set it for one run, and the browser UI shows it
and, on a local install, can change it.

Why this module exists (2026-11-19, a collaborator running the UI on an HPC cluster): ``IGVF_PROJECT_ROOT`` was
honoured by the UI, but nothing in the UI showed it, and tools that open a path as given (the submission commands
among them) resolve it against the directory the program was *started* in. Started inside the working directory,
as ``ask`` usually is, everything worked; a UI started from somewhere else did not find the same files. So: one
directory, entered by every front end (``enter``), visible in the UI (``listing``), told to the agent
(``context``), and changeable from the browser only where that is safe (``may_switch``).

Changing it restarts the UI: many modules fix their paths when they are imported, so a new directory takes a new
process. ``igvfagent ui`` supervises the Streamlit server and starts it again when the app leaves a request in the
file named by ``IGVF_UI_RESTART_FILE`` (``request_switch``).

Pure stdlib and no Streamlit import, so it is testable on its own (``Scripts/test_workdir.py``).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

try:
    from . import _pathguard
except ImportError:                                            # run as a script next to its siblings
    import _pathguard  # type: ignore

__all__ = ["root", "is_explicit", "enter", "listing", "context", "validate", "may_switch", "request_switch",
           "read_request", "resolve_input"]

RESTART_ENV = "IGVF_UI_RESTART_FILE"
_LOOPBACK = {"127.0.0.1", "localhost", "::1", ""}
_COUNT_CAP = 1000                 # stop counting a folder's entries here: a folder of 10^6 files is not listed


def _on(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in ("1", "true", "yes", "on")


def root() -> Path:
    """The workspace root (the same answer as the path guard's)."""
    return _pathguard.project_root()


def is_explicit() -> bool:
    """True when the user chose the directory (``IGVF_PROJECT_ROOT`` is set), false for the default, which is the
    IGVF Agent folder itself."""
    return bool(os.environ.get("IGVF_PROJECT_ROOT", "").strip())


def enter(path: "str | os.PathLike") -> Path:
    """Make ``path`` the working directory of this process and of every tool it starts: the workspace root and
    the current directory. Call before any module that fixes its paths at import (``_agent``, ``agent_jobs``).
    Raises ``NotADirectoryError`` when it is not an existing directory."""
    d = Path(path).expanduser()
    if not d.is_dir():
        raise NotADirectoryError(f"not a directory: {path}")
    d = d.resolve()
    os.environ["IGVF_PROJECT_ROOT"] = str(d)
    os.chdir(d)
    return d


def resolve_input(path: "str | os.PathLike") -> str:
    """A user-supplied input path as the caller should open it.

    Absolute paths and paths that exist relative to the current directory are returned unchanged. A relative path
    that does not exist there but does exist under the workspace root is returned as that absolute path: this is
    the case of a tool started from another directory than the workspace. Anything else is returned unchanged, so
    the caller's own "file not found" message names what the user typed."""
    p = Path(str(path)).expanduser()
    if p.is_absolute() or p.exists():
        return str(path)
    if is_explicit():
        cand = root() / p
        if cand.exists():
            return str(cand)
    return str(path)


def _kind(p: Path) -> Optional[str]:
    """``"dir"``, ``"file"``, ``"link"`` (a symlink that leads out of the workspace) or None for what a listing
    leaves out: hidden entries and secret material (the path guard's deny rules)."""
    if p.name.startswith("."):
        return None
    why = _pathguard.why_blocked(p, require_file=False)
    if why == "outside the workspace":
        # Data kept elsewhere and linked in, the usual layout on a cluster. The file tools do not follow it
        # (the workspace is their boundary); tools that take a path as given do.
        try:
            target = p.resolve()
        except (OSError, RuntimeError):
            return None
        if not p.is_symlink() or _pathguard._DENY_SEGMENTS.intersection(seg.lower() for seg in target.parts):
            return None
        return "link"
    if why is not None and not why.startswith("not a "):
        return None
    return "dir" if p.is_dir() else ("file" if p.is_file() else None)


def listing(limit: int = 200) -> "list[dict]":
    """Top-level entries of the working directory, each ``{"name", "kind", "n", "more", "bytes", "target"}``.

    ``kind`` is ``dir``, ``file`` or ``link`` (a symlink to somewhere outside the workspace, with its
    ``target``). Folders first, then files, each by name; the agent's own output folders (``Docs``, ``Data``)
    last. ``n`` is a folder's number of non-hidden entries, counted up to a cap (``more`` is true beyond it)."""
    base = root()
    out: "list[dict]" = []
    try:
        entries = sorted(os.scandir(base), key=lambda e: e.name.lower())
    except OSError:
        return out
    for e in entries:
        p = Path(e.path)
        kind = _kind(p)
        if kind is None:
            continue
        row = {"name": e.name, "kind": kind, "n": 0, "more": False, "bytes": 0, "target": ""}
        try:
            if kind == "link":
                row["target"] = str(p.resolve())
            elif kind == "dir":
                with os.scandir(p) as it:
                    for sub in it:
                        if sub.name.startswith("."):
                            continue
                        row["n"] += 1
                        if row["n"] >= _COUNT_CAP:
                            row["more"] = True
                            break
            else:
                row["bytes"] = e.stat().st_size
        except OSError:
            continue
        out.append(row)
    own = {"Docs", "Data"}
    out.sort(key=lambda r: (r["name"] in own, r["kind"] == "file", r["name"].lower()))
    return out[:limit]


def context(max_names: int = 40) -> Optional[str]:
    """One bracketed line for the agent about the working directory, or None when the user did not choose one
    (the default root is the software's own folder, whose contents say nothing about the user's data)."""
    if not is_explicit():
        return None
    base = root()
    rows = listing(limit=max_names + 1)
    names = []
    for r in rows[:max_names]:
        if r["kind"] == "dir":
            names.append(f"{r['name']}/ ({r['n']}{'+' if r['more'] else ''} entries)")
        elif r["kind"] == "link":
            names.append(f"{r['name']} (a link to {r['target']}, outside the workspace: give tools that full path)")
        else:
            names.append(r["name"])
    holds = ("It holds: " + ", ".join(names) + (", and more" if len(rows) > max_names else "") + ". "
             if names else "It is empty so far. ")
    return (f"[Working directory: `{base}`. This is the workspace: a relative path the user gives is relative to "
            f"it, tools run in it, and results are written under `Docs/` and `Data/` inside it. {holds}"
            f"Look inside a folder with list_artifacts before saying a file is missing.]")


def validate(path: "str | os.PathLike") -> "tuple[Optional[Path], Optional[str]]":
    """``(resolved directory, None)`` if ``path`` may become the working directory from the browser, else
    ``(None, reason)``. Stricter than ``enter``: the UI's request comes over HTTP, so the whole file system and
    the home directory itself are refused (a project folder inside it is fine), as is anything under a secrets
    directory, and the folder must be writable because results are written into it."""
    raw = str(path or "").strip()
    if not raw:
        return None, "enter a folder path"
    d = Path(raw).expanduser()
    if not d.is_absolute():
        return None, "give the full path, starting with / or ~"
    try:
        d = d.resolve()
    except (OSError, RuntimeError) as exc:
        return None, f"cannot resolve it ({exc.__class__.__name__})"
    if not d.is_dir():
        return None, "no such folder"
    if d == Path(d.anchor):
        return None, "that is the whole file system; choose a project folder"
    try:
        if d == Path.home().resolve():
            return None, "that is your home directory; choose a project folder inside it"
    except (OSError, RuntimeError):
        pass
    if _pathguard._DENY_SEGMENTS.intersection(seg.lower() for seg in d.parts):
        return None, "that folder is for credentials or settings"
    if not os.access(d, os.R_OK | os.X_OK):
        return None, "you cannot read that folder"
    if not os.access(d, os.W_OK):
        return None, "that folder is read-only, and results (Docs/, Data/) are written into the working directory"
    return d, None


def may_switch(*, signed_in: bool = False) -> "tuple[bool, str]":
    """May the browser change the working directory? ``(True, "")`` or ``(False, why)``.

    Only on a local install (never on a shared deployment, never for a signed-in account), only when the UI was
    started by ``igvfagent ui`` (which restarts it), and only when the UI listens on this machine alone: a local
    UI has no sign-in, so on ``--host 0.0.0.0`` anyone who can reach the port could point the agent at any folder
    of the account that runs it. ``IGVF_UI_ALLOW_WORKDIR_SWITCH=1`` lifts that last rule."""
    if _on("IGVF_PUBLIC_MODE") or signed_in:
        return False, "shared deployment"
    if not os.environ.get(RESTART_ENV, "").strip():
        return False, "not started by `igvfagent ui`"
    host = os.environ.get("IGVF_UI_HOST", "").strip().lower()
    if host not in _LOOPBACK and not _on("IGVF_UI_ALLOW_WORKDIR_SWITCH"):
        return False, f"the UI is reachable from other machines (--host {host})"
    return True, ""


def request_switch(path: "str | os.PathLike", *, signed_in: bool = False) -> "tuple[bool, str]":
    """Leave a request for ``igvfagent ui`` to restart in ``path``. ``(True, resolved path)`` or
    ``(False, reason)``; the caller then stops the server."""
    ok, why = may_switch(signed_in=signed_in)
    if not ok:
        return False, why
    d, err = validate(path)
    if err or d is None:
        return False, err or "not a folder"
    # Written into the file `igvfagent ui` made (mkstemp: ours, mode 0600), never through a second name next to
    # it: a predictable name in a shared temporary directory is something another account could plant a link at.
    try:
        with open(os.environ[RESTART_ENV], "w") as fh:
            fh.write(json.dumps({"workdir": str(d)}))
    except OSError as exc:
        return False, f"could not write the restart request ({exc.__class__.__name__})"
    return True, str(d)


def read_request(path: "str | os.PathLike | None") -> Optional[str]:
    """The directory a restart request names, or None. Empties the file, so a request is honoured once. Used by
    ``igvfagent ui`` after the server exits; the directory is checked again here, by the process that acts on
    it."""
    if not path:
        return None
    f = Path(path)
    try:
        raw = f.read_text().strip()
        if raw:
            f.write_text("")
    except OSError:
        return None
    if not raw:
        return None
    try:
        want = json.loads(raw).get("workdir")
    except (ValueError, AttributeError):
        return None
    d, err = validate(want or "")
    return str(d) if d is not None and not err else None
