"""User-extension loader — plug your own skills and tools into IGVFagent.

IGVFagent is expandable: beyond the built-in registry, users can drop
their own capabilities into well-known directories and the runtime
absorbs them at startup with **zero core-code edits**. Two extension
kinds are discovered here:

  * **Tools** — declarative YAML / JSON manifests describing a command
    the LLM agent may call (any executable, or an ``igvfagent``
    subcommand). Merged into the ``_tools`` registry, so they appear in
    ``igvfagent tools``, the ``ask`` ReAct loop, and the Streamlit UI
    tool picker exactly like built-ins.
  * **Skills** — plain Python modules exposing ``main()``. Registered
    as first-class ``igvfagent <name>`` subcommands by the CLI
    dispatcher, with the same post-run local-KG harvest as built-ins.

Discovery locations (all optional, scanned in order; first definition
of a name wins):

  1. every directory in ``$IGVF_USER_EXT_DIR`` (``os.pathsep``-separated)
  2. ``~/.igvfagent/``                      (per-user, cross-checkout)
  3. ``<project root>/UserExtensions/``     (per-checkout, committable)

Each location uses the same layout::

    tools/my_tool.yaml      # or .yml / .json — one manifest per file,
                            # or several under a top-level ``tools:`` list
    skills/my_skill.py      # module with main(); subcommand name is the
                            # filename stem with ``_`` -> ``-``

Tool manifest fields (see Docs/Examples/user_extensions/):

  name          required  snake_case identifier the LLM calls
  description   required  one-paragraph natural-language description
  parameters    optional  JSON Schema object (default: no parameters)
  command       one of    argv list (or shell-style string) for any
                          executable on the user's machine
  cli           these     argv tail run as ``igvfagent <cli...>``
  positional    optional  parameter names emitted as positional args
  flag_map      optional  parameter name -> CLI flag ("" = positional)
  flag_repeat   optional  list-valued params repeated as ``--flag v``
  bool_flags    optional  params emitted as bare flags when true

Custom tools should announce artefacts on stdout as ``Report: <path>`` /
``Manifest: <path>`` lines — the same contract every built-in honours —
so the agent can hand paths to the next step.

Everything here is defensive: a malformed manifest or skill file is
skipped with a warning (inspect via ``igvfagent extensions``) and can
never break the core CLI or agent runtime.
"""

from __future__ import annotations

import ast
import json
import logging
import os
import shlex
import sys
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

_TOOL_SUFFIXES = (".yaml", ".yml", ".json")

# Populated by discover_tools() / discover_skills(); surfaced by
# ``igvfagent extensions`` so users can debug a manifest that was skipped.
# De-duplicated on append: discovery re-runs in long-lived processes
# (Streamlit rerenders every interaction), and repeats must not accumulate.
_PROBLEMS: "list[str]" = []


def _note(problem: str) -> None:
    if problem not in _PROBLEMS:
        _PROBLEMS.append(problem)


# --------------------------- Locations --------------------------------------


def _project_root() -> Path:
    env_root = os.environ.get("IGVF_PROJECT_ROOT")
    if env_root:
        return Path(env_root).expanduser()
    return Path(__file__).resolve().parents[1]


def extension_dirs() -> "list[Path]":
    """Ordered, de-duplicated list of user-extension base directories."""
    dirs: "list[Path]" = []
    env = os.environ.get("IGVF_USER_EXT_DIR", "")
    for token in env.split(os.pathsep):
        if token.strip():
            dirs.append(Path(token.strip()).expanduser())
    dirs.append(Path.home() / ".igvfagent")
    dirs.append(_project_root() / "UserExtensions")
    seen: "set[Path]" = set()
    out: "list[Path]" = []
    for d in dirs:
        try:
            key = d.resolve()
        except OSError:
            key = d
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def problems() -> "list[str]":
    """Diagnostics collected during the last discovery pass."""
    return list(_PROBLEMS)


# --------------------------- Tool manifests ----------------------------------


def _load_manifest(path: Path) -> Any:
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".json":
        return json.loads(text)
    try:
        import yaml  # noqa: F401  (optional dep; playbooks already use it)
    except ImportError as exc:
        raise ValueError(
            "PyYAML is not installed — `pip install pyyaml` or use a "
            ".json manifest instead"
        ) from exc
    return yaml.safe_load(text)


def _normalize_tool(raw: Any, source: Path) -> "Optional[dict]":
    """Validate one manifest entry -> plain dict, or None (recorded)."""
    where = str(source)
    if not isinstance(raw, dict):
        _note(f"{where}: entry is not a mapping")
        return None
    name = raw.get("name")
    desc = raw.get("description")
    if not name or not isinstance(name, str):
        _note(f"{where}: missing required `name`")
        return None
    if not desc or not isinstance(desc, str):
        _note(f"{where}: tool `{name}` missing `description`")
        return None
    command = raw.get("command")
    cli = raw.get("cli")
    if bool(command) == bool(cli):
        _note(
            f"{where}: tool `{name}` needs exactly one of `command` "
            f"(any executable) or `cli` (igvfagent subcommand argv)"
        )
        return None
    if isinstance(command, str):
        command = shlex.split(command)
    if isinstance(cli, str):
        cli = shlex.split(cli)
    params = raw.get("parameters") or {"type": "object", "properties": {}}
    if not isinstance(params, dict) or params.get("type", "object") != "object":
        _note(f"{where}: tool `{name}` `parameters` must be a "
              f"JSON Schema object")
        return None
    params.setdefault("type", "object")
    params.setdefault("properties", {})
    return {
        "name":        name,
        "description": desc.strip(),
        "parameters":  params,
        "command":     [str(t) for t in (command or [])],
        "cli":         [str(t) for t in (cli or [])],
        "positional":  [str(p) for p in (raw.get("positional") or [])],
        "flag_map":    {str(k): str(v)
                        for k, v in (raw.get("flag_map") or {}).items()},
        "flag_repeat": [str(p) for p in (raw.get("flag_repeat") or [])],
        "bool_flags":  [str(p) for p in (raw.get("bool_flags") or [])],
        "source":      where,
    }


def discover_tools() -> "list[dict]":
    """Scan every extension dir for tool manifests.

    Returns plain dicts (not ``_tools.Tool`` instances) so this module
    stays import-cycle-free; ``_tools`` converts and merges them.
    """
    found: "list[dict]" = []
    seen: "set[str]" = set()
    for base in extension_dirs():
        tdir = base / "tools"
        if not tdir.is_dir():
            continue
        for path in sorted(tdir.iterdir()):
            if path.suffix.lower() not in _TOOL_SUFFIXES or not path.is_file():
                continue
            try:
                doc = _load_manifest(path)
            except Exception as exc:
                _note(f"{path}: unreadable manifest ({exc})")
                continue
            entries = doc.get("tools") if isinstance(doc, dict) and \
                isinstance(doc.get("tools"), list) else [doc]
            for raw in entries:
                tool = _normalize_tool(raw, path)
                if tool is None:
                    continue
                if tool["name"] in seen:
                    _note(
                        f"{path}: duplicate user tool `{tool['name']}` "
                        f"(first definition wins)")
                    continue
                seen.add(tool["name"])
                found.append(tool)
    return found


# --------------------------- Python skills -----------------------------------


def _first_docstring_line(path: Path) -> str:
    try:
        doc = ast.get_docstring(ast.parse(path.read_text(encoding="utf-8")))
    except Exception:
        return ""
    return (doc or "").strip().splitlines()[0].strip() if doc else ""


def discover_skills() -> "dict[str, dict]":
    """Scan every extension dir for ``skills/*.py`` modules.

    Returns ``{subcommand: {"path": str, "description": str}}``. The
    subcommand is the filename stem with underscores mapped to hyphens
    (``variant_qc.py`` -> ``igvfagent variant-qc``). Files starting
    with ``_`` are treated as private helpers and skipped.
    """
    found: "dict[str, dict]" = {}
    for base in extension_dirs():
        sdir = base / "skills"
        if not sdir.is_dir():
            continue
        for path in sorted(sdir.glob("*.py")):
            if path.name.startswith("_"):
                continue
            name = path.stem.replace("_", "-")
            if name in found:
                _note(
                    f"{path}: duplicate user skill `{name}` "
                    f"(first definition wins: {found[name]['path']})")
                continue
            found[name] = {
                "path":        str(path),
                "description": _first_docstring_line(path) or
                               "(user skill — no docstring)",
            }
    return found


def load_skill(name: str, entry: dict):
    """Import a discovered user skill module from its file path."""
    import importlib.util
    path = Path(entry["path"])
    mod_name = f"igvfagent_user_{path.stem}"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load user skill `{name}` from {path}")
    module = importlib.util.module_from_spec(spec)
    # Registering in sys.modules lets the skill do relative resource
    # lookups and keeps tracebacks readable.
    sys.modules[mod_name] = module
    spec.loader.exec_module(module)
    return module


__all__ = [
    "extension_dirs", "discover_tools", "discover_skills",
    "load_skill", "problems",
]
