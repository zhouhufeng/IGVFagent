"""Skill cards — a machine-readable specification for every IGVFagent skill.

A reviewer of this work asked the question that separates a feature of one
system from reusable scientific infrastructure: *who can create a skill, how is
its scientific validity established, how is it updated when an archive or
method changes, and how can it be reused outside this particular agent?*

A skill card answers those in a form other software can read. Each card
declares the skill's identity and version, its typed interfaces (shell
subcommand and LLM tool schemas), the upstream archives it depends on, the
artefacts it produces, and — the field that matters most — **its validation
status, derived from the benchmark suite rather than asserted**.

Everything is derived from the running system, not hand-maintained:

* identity and description   <- ``cli.SKILLS``
* tool schemas               <- ``_tools`` registry (the exact schemas the
                                model sees)
* upstream resources         <- ``_resolve_endpoint`` calls in the module
* artefact kinds             <- ``Report:``/``Manifest:``/``Wrote:`` emissions
* documentation              <- ``Docs/Skills/*.md``
* **validation**             <- which ``Benchmarks/*/run.sh`` invoke the skill,
                                with the paper DOI/PMID and check count from
                                that benchmark's ``expected.json``

Deriving validation is the point. A hand-written card would say a skill is
validated because its author believed so; this one says a skill is validated
because a named benchmark, reproducing a named paper, executes it and asserts
N machine-checked criteria. A skill no benchmark touches is reported
``unvalidated`` — which is information, not an omission, and is the honest
status of most of the library.

Usage::

    igvfagent skillcard list                    # one line per skill
    igvfagent skillcard show --skill mavedb     # one full card
    igvfagent skillcard export --out cards/     # all cards as JSON
    igvfagent skillcard coverage                # validation coverage report

Pure standard library.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

try:
    from igvfagent import cli as _cli
    from igvfagent import _tools
    from igvfagent import __version__
except Exception:  # pragma: no cover - direct-script execution
    import cli as _cli  # type: ignore
    import _tools  # type: ignore
    try:
        from __init__ import __version__  # type: ignore
    except Exception:
        __version__ = "0.0.0"

__all__ = ["main", "build_card", "build_all_cards", "SCHEMA_VERSION"]

SCHEMA_VERSION = "1.0"

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
SCRIPTS = ROOT / "Scripts"
BENCH = ROOT / "Benchmarks"
SKILL_DOCS = ROOT / "Docs" / "Skills"

_ENDPOINT_RE = re.compile(r'_resolve_endpoint\(\s*"([a-z_]+)"')
_ARTEFACT_RE = re.compile(r'"(Report|Manifest|Wrote|Output):')


# ---------------------------------------------------------------------------
# Derivation helpers
# ---------------------------------------------------------------------------

def _module_path(module: str) -> "Path | None":
    """igvfagent.foo_skill -> Scripts/foo_skill.py"""
    stem = module.rsplit(".", 1)[-1]
    p = SCRIPTS / f"{stem}.py"
    return p if p.is_file() else None


def _module_source(module: str) -> str:
    p = _module_path(module)
    try:
        return p.read_text(encoding="utf-8", errors="replace") if p else ""
    except OSError:
        return ""


def _upstream_resources(src: str) -> "list[str]":
    return sorted(set(_ENDPOINT_RE.findall(src)))


def _artefact_kinds(src: str) -> "list[str]":
    return sorted(set(_ARTEFACT_RE.findall(src)))


def _summary_line(src: str) -> str:
    """First sentence of the module docstring."""
    m = re.match(r'\s*(?:"""|\'\'\')(.+?)(?:\n|"""|\'\'\')', src, re.S)
    return " ".join(m.group(1).split()) if m else ""


def _doc_path(skill: str, module: str) -> "str | None":
    """Best-effort match into Docs/Skills/."""
    if not SKILL_DOCS.is_dir():
        return None
    stem = module.rsplit(".", 1)[-1].replace("_skill", "").replace("_skills", "")
    wanted = {skill.replace("-", ""), stem.replace("_", "")}
    for p in sorted(SKILL_DOCS.glob("*.md")):
        key = p.stem.lower().replace("_", "").replace("skills", "").replace("skill", "")
        if key in wanted or any(w and w in key for w in wanted):
            return str(p.relative_to(ROOT))
    return None


def _tools_for(skill: str) -> "list[dict]":
    """Every LLM tool whose argv begins with this skill's subcommand."""
    out = []
    for t in _tools.list_tools():
        if t.cli and t.cli[0] == skill:
            out.append({
                "name": t.name,
                "subcommand": " ".join(t.cli),
                "description": (t.description or "").lstrip("★ ").strip(),
                "parameters": t.parameters or {"type": "object", "properties": {}},
            })
    return sorted(out, key=lambda d: d["name"])


def _benchmark_index() -> "dict[str, list[dict]]":
    """skill -> benchmarks that execute it.

    Built by reading which igvfagent subcommands (or Scripts/*.py modules) each
    run.sh actually invokes, so the link is evidence of execution rather than a
    curated assertion that can drift from reality.
    """
    idx: "dict[str, list[dict]]" = {}
    if not BENCH.is_dir():
        return idx
    known = set(getattr(_cli, "SKILLS", {}))
    for d in sorted(BENCH.iterdir()):
        run = d / "run.sh"
        if not run.is_file():
            continue
        try:
            text = run.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

        invoked = set(re.findall(r"igvfagent\s+([a-z][a-z0-9-]+)", text))
        for m in re.findall(r"Scripts/([a-z0-9_]+)\.py", text):
            for name, (mod, _desc) in getattr(_cli, "SKILLS", {}).items():
                if mod.rsplit(".", 1)[-1] == m:
                    invoked.add(name)
        invoked &= known
        if not invoked:
            continue

        meta = {}
        exp = d / "expected.json"
        if exp.is_file():
            try:
                meta = json.loads(exp.read_text())
            except (OSError, json.JSONDecodeError):
                meta = {}
        entry = {
            "benchmark": d.name,
            "paper": meta.get("paper", {}),
            "checks": len(meta.get("checks") or []),
            "run": str((d / "run.sh").relative_to(ROOT)),
        }
        for s in invoked:
            idx.setdefault(s, []).append(entry)
    return idx


# ---------------------------------------------------------------------------
# Card construction
# ---------------------------------------------------------------------------

def build_card(skill: str, *, bench_index=None) -> dict:
    skills = getattr(_cli, "SKILLS", {})
    if skill not in skills:
        raise KeyError(f"unknown skill: {skill}")
    module, description = skills[skill]
    src = _module_source(module)
    bench_index = bench_index if bench_index is not None else _benchmark_index()
    evidence = bench_index.get(skill, [])

    return {
        "schema_version": SCHEMA_VERSION,
        "id": skill,
        "version": __version__,
        "title": description,
        "summary": _summary_line(src),
        "module": module,
        "source": str(_module_path(module).relative_to(ROOT)) if _module_path(module) else None,
        "interfaces": {
            "cli": f"igvfagent {skill}",
            "tools": _tools_for(skill),
        },
        "upstream_resources": _upstream_resources(src),
        "artefact_kinds": _artefact_kinds(src),
        "validation": {
            # Derived, never asserted. "unvalidated" is the honest default and
            # is the status of most skills; it means no benchmark executes this
            # skill, not that the skill is known to be wrong.
            "status": "benchmarked" if evidence else "unvalidated",
            "method": ("machine-checked reproduction of published results"
                       if evidence else None),
            "evidence": evidence,
            "n_checks": sum(e["checks"] for e in evidence),
            "scope_note": (
                "Benchmarks execute fixed command sequences with no model in "
                "the loop, so this establishes implementation correctness only "
                "— not that an agent would select this skill unaided. See "
                "Docs/EVALUATION.md."
            ) if evidence else None,
        },
        "documentation": _doc_path(skill, module),
        "license": "Apache-2.0",
    }


def build_all_cards() -> "list[dict]":
    idx = _benchmark_index()
    return [build_card(s, bench_index=idx)
            for s in sorted(getattr(_cli, "SKILLS", {}))]


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_list(args) -> int:
    cards = build_all_cards()
    print(f"{'SKILL':<20} {'TOOLS':>5} {'CHECKS':>6}  VALIDATION   UPSTREAM")
    for c in cards:
        v = c["validation"]
        mark = "benchmarked" if v["status"] == "benchmarked" else "unvalidated"
        up = ",".join(c["upstream_resources"][:3]) or "-"
        print(f"{c['id']:<20} {len(c['interfaces']['tools']):>5} "
              f"{v['n_checks']:>6}  {mark:<12} {up}")
    n_val = sum(1 for c in cards if c["validation"]["status"] == "benchmarked")
    print(f"\n{len(cards)} skills · {n_val} benchmarked · "
          f"{len(cards) - n_val} unvalidated")
    return 0


def cmd_show(args) -> int:
    try:
        print(json.dumps(build_card(args.skill), indent=2))
    except KeyError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0


def cmd_export(args) -> int:
    out = Path(args.out or (ROOT / "Docs" / "SkillCards"))
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    cards = build_all_cards()
    for c in cards:
        (out / f"{c['id']}.json").write_text(json.dumps(c, indent=2) + "\n")
    # A single bundle is what an external agent consumes.
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "generator": f"igvfagent {__version__}",
        "n_skills": len(cards),
        "skills": cards,
    }
    (out / "skillcards.json").write_text(json.dumps(bundle, indent=2) + "\n")
    print(f"Wrote: {out / 'skillcards.json'}")
    print(f"Output: {out}  ({len(cards)} cards)")
    return 0


def cmd_coverage(args) -> int:
    cards = build_all_cards()
    val = [c for c in cards if c["validation"]["status"] == "benchmarked"]
    unval = [c for c in cards if c["validation"]["status"] != "benchmarked"]
    print(f"# Skill validation coverage — igvfagent {__version__}\n")
    print(f"- Skills: **{len(cards)}**")
    print(f"- Benchmarked: **{len(val)}** ({100*len(val)//max(1,len(cards))}%)")
    print(f"- Unvalidated: **{len(unval)}**")
    print(f"- Machine-checked criteria: **{sum(c['validation']['n_checks'] for c in cards)}**\n")
    print("## Benchmarked\n")
    print("| skill | benchmarks | checks |")
    print("|---|---|---|")
    for c in sorted(val, key=lambda c: -c["validation"]["n_checks"]):
        bs = ", ".join(e["benchmark"] for e in c["validation"]["evidence"])
        print(f"| `{c['id']}` | {bs} | {c['validation']['n_checks']} |")
    print("\n## Unvalidated\n")
    print("No benchmark in the suite executes these skills. That is a "
          "statement about coverage, not about correctness.\n")
    print(", ".join(f"`{c['id']}`" for c in unval))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent skillcard",
        description="Machine-readable specification for every skill, with "
                    "validation status derived from the benchmark suite.")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="One line per skill")
    s = sub.add_parser("show", help="Full card for one skill")
    s.add_argument("--skill", required=True)
    e = sub.add_parser("export", help="Write all cards as JSON")
    e.add_argument("--out", help="Output directory (default Docs/SkillCards)")
    sub.add_parser("coverage", help="Validation coverage report (markdown)")

    args = p.parse_args(argv)
    return {"list": cmd_list, "show": cmd_show,
            "export": cmd_export, "coverage": cmd_coverage}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
