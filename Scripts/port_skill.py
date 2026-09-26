"""Absorb a paper's own analysis code into IGVFagent as a built-in method.

Paper2Agent wraps a paper's repository and runs the authors' scripts as they
are. IGVFagent goes one step further when a paper needs an analysis no
existing command covers: the method is rewritten as a Python module, verified
against the authors' own code on the same input, and registered here so every
later paper that needs it calls ``igvfagent <name>`` directly, with no porting.

Subcommands::

    igvfagent port register --name N --paper-id P --repo URL --commit SHA \\
        --upstream path/in/repo.R [--upstream ...] --source-file port.py \\
        --description "..." [--tool-parameters JSON] [--keywords "a,b"]
    igvfagent port find    --query "pseudobulk differential expression"
    igvfagent port list    [--json]
    igvfagent port remove  --name N

``register`` writes ``Scripts/ported/skills/<name>.py`` and
``Scripts/ported/tools/<name>.json``, which ``_userext``/``_tools``/``cli``
load as built-ins, and records provenance (paper, repository, pinned commit,
the upstream files it ports, source hash, static-scan flags) in
``Scripts/ported/registry.json``. A port is usable at once but stays
``reviewed: false`` until a human looks at it; ``list`` shows each port's
verification status from the ``bench verify-port`` runs that name it.

Registering is writing code this host will later run, so it is refused unless
``IGVF_ALLOW_AGENT_AUTHORING=1`` (the same switch ``extauthor`` uses).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

try:
    from igvfagent import ext_author_skill as _ea
except Exception:  # pragma: no cover - direct-script execution
    import ext_author_skill as _ea  # type: ignore

PORTED = HERE / "ported"
REGISTRY = PORTED / "registry.json"
VERIFY_GLOB = "Docs/Benchmark/2*_port_{slug}/validation_vs_reference.json"

__all__ = ["main"]


def _rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def _load_registry() -> dict:
    try:
        return json.loads(REGISTRY.read_text())
    except (OSError, ValueError):
        return {}


def _save_registry(reg: dict) -> None:
    PORTED.mkdir(parents=True, exist_ok=True)
    REGISTRY.write_text(json.dumps(reg, indent=2, sort_keys=True) + "\n")


def _scan(text: str) -> list:
    try:
        import extension_review as _er
        return _er.scan(text)
    except Exception:  # noqa: BLE001 - the scan is advisory
        return []


def _reserved_names() -> set:
    """Built-in commands and tools a port must never shadow."""
    names: set = set()
    try:
        from igvfagent import cli as _cli   # relative imports: only as a package
        names |= {n.replace("-", "_") for n in _cli.SKILLS}
    except Exception:  # noqa: BLE001
        pass
    try:
        import _tools
        ported = set(_load_registry())
        names |= {n for n in _tools._BY_NAME
                  if n not in _tools._USER_TOOL_NAMES and n not in ported}
    except Exception:  # noqa: BLE001
        pass
    return names


def verifications(name: str) -> list:
    """Every ``bench verify-port`` result recorded for this port, newest first."""
    slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")[:40]
    out = []
    for p in sorted(ROOT.glob(VERIFY_GLOB.format(slug=slug)), reverse=True):
        try:
            v = json.loads(p.read_text())["validation_vs_reference"]
        except (OSError, ValueError, KeyError):
            continue
        out.append({"path": str(p.relative_to(ROOT)), "passed": v.get("passed"),
                    "match_rate": v.get("match_rate"), "paper_id": v.get("paper_id"),
                    "analysis": v.get("analysis")})
    return out


def register(*, name: str, paper_id: str, repo: str, commit: str,
             upstream: list, source_file: str, description: str,
             tool_parameters: "str | None" = None, flag_map: "str | None" = None,
             keywords: "str | None" = None) -> dict:
    _ea._require_enabled()
    if not _ea._NAME_RE.match(name or ""):
        raise _ea.InvalidExtension(
            f"Invalid name {name!r}: lowercase letters, digits and underscores, "
            "3-49 chars, starting with a letter.")
    if name in _reserved_names():
        raise _ea.InvalidExtension(f"{name!r} is already a built-in command or tool.")
    if not (description or "").strip():
        raise _ea.InvalidExtension("A port needs a --description.")
    if not upstream:
        raise _ea.InvalidExtension(
            "Name the authors' files this ports with --upstream; a port with no "
            "upstream source is new science, not a port.")
    if not re.fullmatch(r"[0-9a-f]{7,40}", commit or ""):
        raise _ea.InvalidExtension("--commit must be the pinned git SHA of the repository.")
    src = Path(source_file).expanduser()
    if not src.is_file():
        raise _ea.InvalidExtension(f"--source-file not found: {src}")
    source = src.read_text()
    try:
        compile(source, f"<{name}>", "exec")
    except SyntaxError as e:
        raise _ea.InvalidExtension(f"Syntax error at line {e.lineno}: {e.msg}") from e
    if not re.search(r"(?m)^def main\b", source):
        raise _ea.InvalidExtension("A port must define a top-level main(argv=None).")

    schema = json.loads(tool_parameters) if tool_parameters else {"type": "object", "properties": {}}
    if schema.get("type") != "object":
        raise _ea.InvalidExtension("--tool-parameters must be a JSON Schema object.")
    fmap = json.loads(flag_map) if flag_map else {
        k: "--" + k.replace("_", "-") for k in schema.get("properties", {})}

    sha = hashlib.sha256(source.encode()).hexdigest()
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    header = (f"# Ported by `igvfagent port register` on {now[:10]} from {repo} @ {commit[:12]}\n"
              f"# ({', '.join(upstream)}) for {paper_id}. Unreviewed; provenance in "
              f"Scripts/ported/registry.json.\n")
    sdir, tdir = PORTED / "skills", PORTED / "tools"
    sdir.mkdir(parents=True, exist_ok=True)
    tdir.mkdir(parents=True, exist_ok=True)
    (sdir / f"{name}.py").write_text(header + source + ("" if source.endswith("\n") else "\n"))
    manifest = {"name": name, "description": description.strip(), "parameters": schema,
                "cli": [name.replace("_", "-")], "x-ported": {"paper_id": paper_id,
                                                              "repo": repo, "commit": commit}}
    if fmap:
        manifest["flag_map"] = fmap
    (tdir / f"{name}.json").write_text(json.dumps(manifest, indent=2) + "\n")

    reg = _load_registry()
    prev = reg.get(name) or {}
    reg[name] = {
        "paper_id": paper_id, "repo": repo, "commit": commit, "upstream": sorted(upstream),
        "description": description.strip(),
        "keywords": [k.strip() for k in (keywords or "").split(",") if k.strip()],
        "source_sha256": sha, "risk_flags": _scan(source),
        "registered_at": prev.get("registered_at") or now, "updated_at": now,
        "reviewed": False, "reviewer": None,
    }
    _save_registry(reg)
    return {"name": name, "command": f"igvfagent {name.replace('_', '-')}",
            "files": [_rel(sdir / f"{name}.py"), _rel(tdir / f"{name}.json")],
            "risk_flags": reg[name]["risk_flags"]}


def find(query: str, top: int = 5) -> list:
    """Ported methods ranked by word overlap with ``query``: reuse before porting."""
    words = set(re.findall(r"[a-z0-9]{3,}", query.lower()))
    hits = []
    for name, e in _load_registry().items():
        text = " ".join([name.replace("_", " "), e.get("description", ""),
                         " ".join(e.get("keywords") or []), " ".join(e.get("upstream") or [])])
        have = set(re.findall(r"[a-z0-9]{3,}", text.lower()))
        score = len(words & have)
        if score:
            hits.append({"name": name, "score": score, "command": f"igvfagent {name.replace('_', '-')}",
                         "description": e.get("description"), "paper_id": e.get("paper_id")})
    return sorted(hits, key=lambda h: -h["score"])[:top]


def remove(name: str) -> list:
    _ea._require_enabled()
    reg = _load_registry()
    if name not in reg:
        raise _ea.InvalidExtension(f"No port named {name!r}.")
    gone = []
    for p in (PORTED / "skills" / f"{name}.py", PORTED / "tools" / f"{name}.json"):
        if p.is_file():
            p.unlink()
            gone.append(_rel(p))
    del reg[name]
    _save_registry(reg)
    return gone


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="igvfagent port", description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("register", help="Register a Python port as a built-in command + tool.")
    s.add_argument("--name", required=True)
    s.add_argument("--paper-id", required=True)
    s.add_argument("--repo", required=True, help="The authors' repository URL.")
    s.add_argument("--commit", required=True, help="Pinned commit SHA the port follows.")
    s.add_argument("--upstream", action="append", required=True,
                   help="Repo path of the code this ports (repeatable).")
    s.add_argument("--source-file", required=True)
    s.add_argument("--description", required=True)
    s.add_argument("--tool-parameters", help="JSON Schema object for the tool's inputs.")
    s.add_argument("--flag-map", help="JSON {param: --flag}; default --param-name.")
    s.add_argument("--keywords", help="Comma-separated terms `port find` should match.")
    s = sub.add_parser("find", help="Find a ported method to reuse.")
    s.add_argument("--query", required=True)
    s.add_argument("--top", type=int, default=5)
    s = sub.add_parser("list", help="Ported methods with provenance and verification.")
    s.add_argument("--json", action="store_true")
    s = sub.add_parser("remove", help="Retire a port.")
    s.add_argument("--name", required=True)
    a = p.parse_args(argv)

    try:
        if a.cmd == "register":
            out = register(name=a.name, paper_id=a.paper_id, repo=a.repo, commit=a.commit,
                           upstream=a.upstream, source_file=a.source_file,
                           description=a.description, tool_parameters=a.tool_parameters,
                           flag_map=a.flag_map, keywords=a.keywords)
            print(f"Registered {out['command']}  ({', '.join(out['files'])})")
            if out["risk_flags"]:
                print("Static-scan flags (recorded, for review): "
                      + ", ".join(f"{f['flag']}({f['level']})" for f in out["risk_flags"]))
            print("Verify it against the authors' code before relying on it: "
                  "igvfagent bench verify-port --name " + a.name + " …")
            return 0
        if a.cmd == "find":
            hits = find(a.query, a.top)
            for h in hits:
                print(f"{h['score']:3d}  {h['command']:40s} {h['description']}  [{h['paper_id']}]")
            if not hits:
                print("No ported method matches; port the authors' code (igvf-replicate-paper step 4).")
            return 0
        if a.cmd == "list":
            reg = _load_registry()
            rows = [{"name": n, **e, "verifications": verifications(n)} for n, e in sorted(reg.items())]
            if a.json:
                print(json.dumps(rows, indent=2))
                return 0
            for r in rows:
                v = r["verifications"]
                vs = (f"verified {'PASS' if v[0]['passed'] else 'FAIL'} "
                      f"(match_rate {v[0]['match_rate']:.3f})") if v else "not verified"
                print(f"{r['name']:32s} {r['paper_id']:36s} {vs}; "
                      f"{'reviewed' if r.get('reviewed') else 'unreviewed'}")
            if not rows:
                print("No ported methods yet.")
            return 0
        if a.cmd == "remove":
            print("Removed " + ", ".join(remove(a.name)))
            return 0
    except (_ea.AuthoringDisabled, _ea.InvalidExtension, ValueError) as e:
        print(f"port {a.cmd}: {e}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
