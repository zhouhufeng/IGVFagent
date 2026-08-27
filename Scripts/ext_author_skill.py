"""Agent-authored extensions — let IGVFagent write its own tools and skills.

IGVFagent's own agent can call 148 pre-registered tools but cannot *create*
one: every tool call is turned into a fixed argv from a declared schema
(``_tools._build_argv``), and there is no file-write, shell, or code-execution
tool in the registry. That is the entire reason the in-app agent feels less
capable than a coding assistant driving the same model — it is the harness,
not the model.

This skill closes that gap using the extension mechanism that already exists
(``_userext``): the agent writes a tool manifest or a Python skill module into
a user-extension directory, and the runtime absorbs it into the registry with
no restart and no core-code edit. The agent can then call what it just wrote.

**This is deliberately opt-in.** An authored skill is Python that IGVFagent
later executes as a subprocess, so enabling it turns "can chat with the agent"
into "can run code on this host". On a shared deployment behind one shared
password that is a meaningful escalation. Set
``IGVF_ALLOW_AGENT_AUTHORING=1`` only where that is acceptable; every write
path below refuses without it.

Subcommands::

    igvfagent extauthor write-tool  --name t --description d --cli "kg gene"
    igvfagent extauthor write-skill --name s --description d --source "..."
    igvfagent extauthor validate    --name s
    igvfagent extauthor list
    igvfagent extauthor remove      --name s --kind skill

Pure standard library (``json`` + ``re`` + ``compile``); PyYAML is used only
if present, and manifests are written as JSON — which every YAML parser reads —
so the skill never depends on it.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

try:
    from igvfagent import _userext
except Exception:  # pragma: no cover - direct-script execution
    import _userext  # type: ignore

__all__ = ["main"]

# A name must be a safe identifier: it becomes a filename *and* a CLI
# subcommand *and* a tool name the LLM emits. Anchored, so no traversal
# ("../x"), no dotfiles, no extensions, no shell metacharacters.
_NAME_RE = re.compile(r"^[a-z][a-z0-9_]{2,48}$")

_ENABLE_ENV = "IGVF_ALLOW_AGENT_AUTHORING"


class AuthoringDisabled(RuntimeError):
    pass


class InvalidExtension(ValueError):
    pass


def authoring_enabled() -> bool:
    return os.environ.get(_ENABLE_ENV, "").strip().lower() in ("1", "true", "yes", "on")


def _require_enabled() -> None:
    if not authoring_enabled():
        raise AuthoringDisabled(
            f"Extension authoring is disabled. Set {_ENABLE_ENV}=1 to allow "
            "IGVFagent to write new tools/skills. It is off by default because "
            "an authored skill is Python this host will later execute."
        )


def _target_dir() -> Path:
    """First writable extension directory, creating it if needed.

    Uses the same search order as ``_userext`` so anything written here is
    discovered by the normal loader with no extra configuration.
    """
    for d in _userext.extension_dirs():
        try:
            d.mkdir(parents=True, exist_ok=True)
            probe = d / ".write_probe"
            probe.touch()
            probe.unlink()
            return d
        except OSError:
            continue
    raise InvalidExtension(
        "No writable extension directory found; tried: "
        + ", ".join(str(d) for d in _userext.extension_dirs())
    )


def known_subcommands() -> "set[str]":
    """Every `igvfagent <sub>` that actually resolves — built-in and user.

    ``write_tool`` validates against this. Without the check the agent can
    write a manifest pointing at a subcommand it never implemented, which is
    exactly what happened in practice: ``catfile_tool`` -> `igvfagent catfile`,
    ``run_create_variant_csv`` -> `igvfagent create-variant-csv`. Both parsed
    as valid manifests, registered cleanly, and failed only when called.
    """
    names: "set[str]" = set()
    try:
        from igvfagent import cli
    except Exception:  # pragma: no cover
        try:
            import cli  # type: ignore
        except Exception:
            return names
    names |= set(getattr(cli, "SKILLS", {}))
    names |= set(getattr(cli, "TOP_LEVEL", ()))
    try:
        names |= set(_userext.discover_skills())
    except Exception:
        pass
    return names


def _check_name(name: str, *, kind: str) -> str:
    name = (name or "").strip()
    if not _NAME_RE.match(name):
        raise InvalidExtension(
            f"Invalid {kind} name {name!r}. Use lowercase letters, digits and "
            "underscores, 3-49 chars, starting with a letter."
        )
    # Never let an authored extension shadow a built-in tool: _userext gives
    # "first definition wins" to built-ins, so a clash would silently no-op.
    try:
        from igvfagent import _tools
    except Exception:  # pragma: no cover
        import _tools  # type: ignore
    if any(t.name == name for t in _tools.list_tools()):
        raise InvalidExtension(
            f"{name!r} is already a built-in tool name — pick another; a "
            "shadowing extension would be ignored by the loader."
        )
    return name


# --------------------------------------------------------------------------
# write-tool — a declarative manifest wrapping an existing command
# --------------------------------------------------------------------------

def write_tool(*, name: str, description: str,
               cli: "str | None" = None,
               command: "str | None" = None,
               parameters: "str | None" = None,
               positional: "str | None" = None) -> Path:
    """Write a user-tool manifest the agent can immediately call."""
    _require_enabled()
    name = _check_name(name, kind="tool")
    if not (description or "").strip():
        raise InvalidExtension("A tool needs a description — it is the only "
                               "thing the model uses to decide when to call it.")
    if bool(cli) == bool(command):
        raise InvalidExtension(
            "Provide exactly one of --cli (an `igvfagent <...>` subcommand) or "
            "--command (argv for any executable).")

    schema: dict = {"type": "object", "properties": {}}
    if parameters:
        try:
            schema = json.loads(parameters)
        except json.JSONDecodeError as e:
            raise InvalidExtension(f"--parameters is not valid JSON: {e}") from e
        if not isinstance(schema, dict) or schema.get("type") != "object":
            raise InvalidExtension(
                "--parameters must be a JSON Schema object, e.g. "
                '{"type":"object","properties":{"gene":{"type":"string"}}}')

    manifest: dict = {"name": name, "description": description.strip(),
                      "parameters": schema}
    if cli:
        parts = cli.split()
        sub = parts[0] if parts else ""
        known = known_subcommands()
        if known and sub not in known:
            raise InvalidExtension(
                f"`igvfagent {sub}` does not exist, so this tool would "
                f"register successfully and then fail on every call.\n"
                f"If {sub!r} is meant to be NEW functionality, write the "
                f"implementation first with ext_author_skill (which registers "
                f"its own tool), instead of wrapping a command that isn't "
                f"there. To wrap an arbitrary executable rather than an "
                f"igvfagent subcommand, use --command.")
        manifest["cli"] = parts
    else:
        manifest["command"] = command.split() if isinstance(command, str) else command

    # A tool with no parameters can never receive input. That is legitimate for
    # something like `ext_list`, but it was the second half of the observed
    # failure: run_create_variant_csv declared {} and so could not have been
    # given a variant list even if its subcommand had existed.
    if not schema.get("properties"):
        print(f"note: `{name}` declares no parameters, so the model cannot "
              f"pass it any input. If it needs arguments, re-author it with "
              f"--parameters.", file=sys.stderr)
    if positional:
        manifest["positional"] = positional.split()

    d = _target_dir() / "tools"
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.json"          # JSON is valid YAML — no PyYAML needed
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


# --------------------------------------------------------------------------
# write-skill — a Python module registered as `igvfagent <name>`
# --------------------------------------------------------------------------

_SKILL_TEMPLATE = '''"""{description}

Authored by IGVFagent via `extauthor write-skill`.
"""
from __future__ import annotations

import argparse


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="igvfagent {name}",
                                      description={description!r})
    args = parser.parse_args(argv)
    print("{name}: no implementation yet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def write_skill(*, name: str, description: str,
                source: "str | None" = None,
                source_file: "str | None" = None,
                register_tool: bool = True,
                tool_parameters: "str | None" = None,
                flag_map: "str | None" = None) -> Path:
    """Write a Python skill module exposing ``main()``.

    By default this also writes the matching tool manifest, pointing at the
    subcommand the module actually registers. That pairing is the whole point:
    writing the wrapper and the implementation as two independent steps is how
    every previous attempt failed — a manifest naming a subcommand that was
    never implemented registers cleanly and fails only when called.
    """
    _require_enabled()
    name = _check_name(name, kind="skill")
    if not (description or "").strip():
        raise InvalidExtension("A skill needs a description.")

    if source_file:
        src_path = Path(source_file).expanduser()
        if not src_path.is_file():
            raise InvalidExtension(f"--source-file not found: {src_path}")
        source = src_path.read_text()
    if not source:
        source = _SKILL_TEMPLATE.format(name=name, description=description.strip())

    # Refuse to write code that cannot even be parsed — otherwise the loader
    # swallows the ImportError and the skill silently never appears.
    try:
        compile(source, f"<{name}>", "exec")
    except SyntaxError as e:
        raise InvalidExtension(
            f"Refusing to write {name}: source has a syntax error at line "
            f"{e.lineno}: {e.msg}") from e
    if not re.search(r"(?m)^def main\b", source):
        raise InvalidExtension(
            "A skill module must define a top-level `main()` — that is what "
            "the CLI dispatcher calls.")

    d = _target_dir() / "skills"
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.py"
    path.write_text(source if source.endswith("\n") else source + "\n")

    if register_tool:
        # Write the module FIRST so the subcommand exists by the time the
        # manifest referencing it is discovered.
        schema: dict = {"type": "object", "properties": {}}
        if tool_parameters:
            try:
                schema = json.loads(tool_parameters)
            except json.JSONDecodeError as e:
                raise InvalidExtension(
                    f"--tool-parameters is not valid JSON: {e}") from e
        fmap = {}
        if flag_map:
            try:
                fmap = json.loads(flag_map)
            except json.JSONDecodeError as e:
                raise InvalidExtension(f"--flag-map is not valid JSON: {e}") from e
        elif schema.get("properties"):
            # Sensible default: every declared property becomes --property-name,
            # which is what an argparse-based skill expects.
            fmap = {k: "--" + k.replace("_", "-") for k in schema["properties"]}

        manifest = {"name": name, "description": description.strip(),
                    "parameters": schema, "cli": [_subcommand_for(name)]}
        if fmap:
            manifest["flag_map"] = fmap
        td = _target_dir() / "tools"
        td.mkdir(parents=True, exist_ok=True)
        (td / f"{name}.json").write_text(json.dumps(manifest, indent=2) + "\n")

        if not schema.get("properties"):
            print(f"note: `{name}` was registered with no parameters, so the "
                  f"model cannot pass it input. Re-author with "
                  f"--tool-parameters if it takes arguments.", file=sys.stderr)
    return path


# --------------------------------------------------------------------------
# introspection
# --------------------------------------------------------------------------

def list_extensions() -> dict:
    # discover_tools() yields plain dicts (import-cycle-free by design), and
    # discover_skills() is keyed by SUBCOMMAND — the filename stem with
    # underscores mapped to hyphens, so `my_skill.py` -> `my-skill`.
    tools = _userext.discover_tools()
    skills = _userext.discover_skills()
    return {
        "enabled": authoring_enabled(),
        "dirs": [str(d) for d in _userext.extension_dirs()],
        "tools": sorted(t.get("name", "") for t in tools if t.get("name")),
        "skills": sorted(skills.keys()),
        "problems": _userext.problems(),
    }


def _subcommand_for(name: str) -> str:
    """The CLI subcommand a skill file named ``name`` registers as."""
    return name.replace("_", "-")


def validate(name: str) -> dict:
    """Confirm an authored extension is discoverable and (for skills) imports."""
    tools = _userext.discover_tools()
    skills = _userext.discover_skills()
    sub = _subcommand_for(name)

    out = {
        "name": name,
        "registered_as_tool": any(t.get("name") == name for t in tools),
        "registered_as_skill": sub in skills,
        "subcommand": sub if sub in skills else None,
        "importable": None,
        "error": None,
        "problems": _userext.problems(),
    }
    if out["registered_as_skill"]:
        try:
            _userext.load_skill(sub, skills[sub])
            out["importable"] = True
        except Exception as e:
            out["importable"] = False
            out["error"] = f"{type(e).__name__}: {e}"
    return out


def remove(name: str, kind: str) -> "list[str]":
    _require_enabled()
    removed: "list[str]" = []
    for d in _userext.extension_dirs():
        for sub, ext in (("tools", ".json"), ("tools", ".yaml"),
                          ("tools", ".yml"), ("skills", ".py")):
            if kind != "any" and not sub.startswith(kind.rstrip("s")):
                continue
            p = d / sub / f"{name}{ext}"
            # Containment: never delete outside the extension dir, whatever
            # the name resolved to.
            try:
                p.resolve().relative_to(d.resolve())
            except (ValueError, OSError):
                continue
            if p.is_file():
                p.unlink()
                removed.append(str(p))
    return removed


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="igvfagent extauthor",
        description="Author new IGVFagent tools/skills (opt-in; requires "
                    f"{_ENABLE_ENV}=1).")
    sub = parser.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("write-tool", help="Write a tool manifest")
    t.add_argument("--name", required=True)
    t.add_argument("--description", required=True)
    t.add_argument("--cli", help='igvfagent subcommand tail, e.g. "kg gene"')
    t.add_argument("--command", help="argv for any executable")
    t.add_argument("--parameters", help="JSON Schema object")
    t.add_argument("--positional", help="space-separated parameter names")

    s = sub.add_parser("write-skill", help="Write a Python skill module")
    s.add_argument("--name", required=True)
    s.add_argument("--description", required=True)
    s.add_argument("--source", help="Full Python source (must define main())")
    s.add_argument("--source-file", help="Read the source from this path")
    s.add_argument("--tool-parameters",
                   help="JSON Schema object for the auto-registered tool.")
    s.add_argument("--flag-map",
                   help="JSON map of parameter name -> CLI flag. Defaults to "
                        "--param-name for each declared property.")
    s.add_argument("--no-register-tool", action="store_true",
                   help="Write the module only; do not register a tool.")

    v = sub.add_parser("validate", help="Check an extension is registered")
    v.add_argument("--name", required=True)

    sub.add_parser("list", help="List discovered extensions")

    r = sub.add_parser("remove", help="Delete an authored extension")
    r.add_argument("--name", required=True)
    r.add_argument("--kind", default="any", choices=("any", "tool", "skill"))

    args = parser.parse_args(argv)

    try:
        if args.cmd == "write-tool":
            p = write_tool(name=args.name, description=args.description,
                           cli=args.cli, command=args.command,
                           parameters=args.parameters, positional=args.positional)
            print(f"Wrote: {p}")
            print(f"Tool `{args.name}` is now callable — no restart needed.")
        elif args.cmd == "write-skill":
            p = write_skill(name=args.name, description=args.description,
                            source=args.source, source_file=args.source_file,
                            register_tool=not args.no_register_tool,
                            tool_parameters=args.tool_parameters,
                            flag_map=args.flag_map)
            print(f"Wrote: {p}")
            if not args.no_register_tool:
                print(f"Tool `{args.name}` registered and callable now.")
            # Report the SUBCOMMAND, not the filename: the loader maps
            # underscores to hyphens, and the agent calls what this line says.
            print(f"Skill `{args.name}` is now "
                  f"`igvfagent {_subcommand_for(args.name)}`.")
        elif args.cmd == "validate":
            print(json.dumps(validate(args.name), indent=2))
        elif args.cmd == "list":
            print(json.dumps(list_extensions(), indent=2))
        elif args.cmd == "remove":
            gone = remove(args.name, args.kind)
            print(json.dumps({"removed": gone}, indent=2))
    except AuthoringDisabled as e:
        print(f"error: {e}", file=sys.stderr)
        return 3
    except InvalidExtension as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
