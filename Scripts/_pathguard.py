"""Artefact path containment for the browser UI.

The chat UI renders any file path it finds in tool output or model text:
``_extract_paths_from_text`` scrapes paths out of a report body and the
``_render_*`` helpers then ``read_bytes()`` them into an inline viewer or a
download button. On a single-user laptop that is exactly the feature you
want — a report cites ``Docs/ENCODE/Plots/x.svg`` and the plot appears.

On a **shared deployment** the same feature is an exfiltration path. The
scraper's ``_BARE_ABS_PATH`` matches any absolute path ending in a viewable
extension (``.txt``, ``.json``, ``.log``, ``.md`` …), so a tool that echoes
an unexpected path — or a prompt that talks a model into naming one — turns
into a download button for that file. Two classes of target matter:

* **outside the workspace** — ``/home/ubuntu/.config/…json``, and
* **inside it** — ``Docs/Secret/API-credentials.txt`` is a ``.txt`` living
  *under* the project root, so a root-containment check alone would happily
  serve it.

Hence both rules below: contain to the workspace, then subtract the secret
material inside it. Deny wins over allow.

Pure stdlib, no Streamlit import, so it is unit-testable and usable from any
skill that wants the same check.
"""
from __future__ import annotations

import os
from pathlib import Path

__all__ = ["project_root", "is_safe_artifact", "filter_artifacts", "why_blocked"]


def project_root() -> Path:
    """Workspace root that artefacts must live under."""
    return Path(
        os.environ.get("IGVF_PROJECT_ROOT")
        or Path(__file__).resolve().parents[1]
    ).resolve()


# Path *segments* that never contain renderable artefacts. Matched
# case-insensitively against every component of the resolved path, so
# `Docs/Secret/x.txt`, `docs/secrets/x.txt` and `a/Secretes/b.json` all lose.
_DENY_SEGMENTS = frozenset({
    "secret", "secrets", "secretes",      # incl. the historical typo
    ".ssh", ".aws", ".config", ".git",
    "credentials",
})

# Exact filenames that are never renderable, whatever directory they sit in.
_DENY_NAMES = frozenset({
    ".env", "clouds.yaml", "clouds.yml",
    "api-credentials.txt", "apikeys.txt",
    "id_rsa", "id_ed25519", "authorized_keys",
})

# Suffixes that are never renderable even with an allowed extension upstream.
_DENY_SUFFIXES = (".pem", ".key", ".p12", ".pfx", ".cookie", ".cookies", ".token")

# Substrings in the *filename* that signal credential material.
_DENY_NAME_PARTS = ("apikey", "api_key", "credential", "secret", "password", "token")


def why_blocked(path: "str | Path", *, require_file: bool = True) -> "str | None":
    """Return a short reason the path is not renderable, or ``None`` if it is.

    Resolution follows symlinks, so a symlink planted inside the workspace
    cannot point out of it and slip past the containment check.

    ``require_file=False`` keeps the containment and denylist checks but
    accepts a directory — for callers that list or search a run directory
    rather than reading one file. The security-relevant rules are unchanged;
    only the final "is this a regular file" assertion is relaxed.
    """
    root = project_root()
    try:
        p = Path(path).expanduser().resolve()
    except (OSError, RuntimeError) as e:      # RuntimeError: symlink loop
        return f"unresolvable path ({e.__class__.__name__})"

    # 1) Containment — must sit under the workspace root.
    try:
        p.relative_to(root)
    except ValueError:
        return "outside the workspace"

    # 2) Denylist — secret material *inside* the workspace.
    lowered = [seg.lower() for seg in p.parts]
    if _DENY_SEGMENTS.intersection(lowered):
        return "inside a secrets directory"

    name = p.name.lower()
    if name in _DENY_NAMES:
        return "credential file"
    if name.endswith(_DENY_SUFFIXES):
        return "key material"
    if any(part in name for part in _DENY_NAME_PARTS):
        return "credential-like filename"

    # 3) Must be a real regular file (not a device or FIFO). Directories are
    #    allowed only for callers that opted in via require_file=False.
    if require_file:
        if not p.is_file():
            return "not a regular file"
    elif not (p.is_file() or p.is_dir()):
        return "not a file or directory"

    return None


def is_safe_artifact(path: "str | Path", *, require_file: bool = True) -> bool:
    """True when *path* may be rendered or offered for download."""
    return why_blocked(path, require_file=require_file) is None


def filter_artifacts(paths) -> "list[str]":
    """Keep only the renderable paths, preserving order and dropping dupes."""
    out: "list[str]" = []
    seen: "set[str]" = set()
    for raw in paths or ():
        if is_safe_artifact(raw) and str(raw) not in seen:
            seen.add(str(raw))
            out.append(str(raw))
    return out
