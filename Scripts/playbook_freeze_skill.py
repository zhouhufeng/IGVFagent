"""Freeze a recorded session into a deterministic playbook.

Free-form planning and pinned execution sit at opposite ends of a trade-off:
exploration is where a finding comes from, reproducibility is what publishing
it requires. Measured artefact agreement across repeated free-form runs is far
below 1.0, while a pinned playbook reproduces byte-identically.

This makes the gap a workflow rather than a limitation. Explore in natural
language, then freeze the session you liked into a playbook that any user, on
any backend, re-runs deterministically:

    igvfagent playbook-freeze list
    igvfagent playbook-freeze freeze --run <Docs/Agent/…> --out my.yaml
    igvfagent playbook freeze --check my.yaml     # verify a replay

Frozen playbooks carry **expected artefact hashes**. A replay that produces
different bytes is then a detected event rather than a silent divergence —
which matters because the usual cause is not the agent at all, but an upstream
archive that re-released its data. Recording the hash is what separates "this
replayed" from "this replayed and produced the same thing".

Pure standard library.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

__all__ = ["main", "freeze_run", "load_run"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
AGENT_DIR = ROOT / "Docs" / "Agent"
PLAYBOOK_DIR = ROOT / "Docs" / "Playbooks"


def _sha256(path: Path, *, limit: int = 64 * 1024 * 1024) -> "str | None":
    try:
        if not path.is_file() or path.stat().st_size > limit:
            return None
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def load_run(run: str) -> dict:
    """Load a persisted agent run (directory or transcript.json)."""
    p = Path(run).expanduser()
    if not p.is_absolute():
        p = ROOT / p
    if p.is_dir():
        p = p / "transcript.json"
    if not p.is_file():
        raise SystemExit(f"error: no transcript at {p}")
    return json.loads(p.read_text())


def _steps_from(transcript) -> "list[dict]":
    """Every executed tool call, in order, with its arguments.

    Taken from the assistant turns rather than the tool results: the result
    records which tool ran, but only the assistant entry carries the arguments
    that produced it, and a playbook without arguments is not reproducible.
    """
    steps = []
    for entry in transcript or []:
        if entry.get("role") != "assistant":
            continue
        for tc in entry.get("tool_calls") or []:
            name = tc.get("name")
            if not name:
                continue
            steps.append({"tool": name, "args": dict(tc.get("arguments") or {})})
    return steps


def _artefact_hashes(transcript) -> "list[dict]":
    out, seen = [], set()
    for entry in transcript or []:
        arts = entry.get("artifacts") or {}
        paths = ([p for v in arts.values() for p in v]
                 if isinstance(arts, dict) else list(arts))
        for raw in paths:
            p = Path(raw)
            if not p.is_absolute():
                p = ROOT / p
            key = str(p)
            if key in seen:
                continue
            seen.add(key)
            digest = _sha256(p)
            try:
                rel = str(p.relative_to(ROOT))
            except ValueError:
                rel = str(p)
            out.append({"path": rel, "sha256": digest,
                        "bytes": p.stat().st_size if p.is_file() else None})
    return out


def _parameterise(steps, query: str) -> "tuple[list[dict], dict]":
    """Lift repeated literal values into playbook parameters.

    A frozen session is only useful again if it can be pointed at a different
    gene. Values appearing in more than one step, or in the original question,
    become ${parameters} with the observed value as the default — so the
    playbook replays identically by default and generalises on request.
    """
    # Arguments naming a mode, relationship or method are controlled
    # vocabulary: turning them into parameters invites a replay that asks for
    # a relationship the API does not have. Only subject-like arguments —
    # the gene, variant, region or accession a session was *about* — are
    # lifted.
    SUBJECT_ARGS = {"symbol", "gene", "gene_name", "entity_id", "id",
                    "variant", "rsid", "region", "accession",
                    "accession_or_url", "regulator", "response", "urn",
                    "gse", "query", "target"}
    counts: "dict[str, int]" = {}
    for s in steps:
        for k, v in s["args"].items():
            if k not in SUBJECT_ARGS:
                continue
            if isinstance(v, str) and 2 <= len(v) <= 40:
                counts[v] = counts.get(v, 0) + 1

    qlow = (query or "").lower()
    params: "dict[str, dict]" = {}
    name_for: "dict[str, str]" = {}
    for value, n in counts.items():
        in_query = value.lower() in qlow
        if n < 2 and not in_query:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][\w.:\-]*", value):
            continue
        base = re.sub(r"\W+", "_", value.lower())[:24] or "value"
        pname = base
        i = 2
        while pname in params:
            pname = f"{base}_{i}"; i += 1
        params[pname] = {"default": value}
        name_for[value] = pname

    out = []
    for s in steps:
        args = {}
        for k, v in s["args"].items():
            lift = (k in SUBJECT_ARGS and isinstance(v, str) and v in name_for)
            args[k] = f"${{{name_for[v]}}}" if lift else v
        out.append({"tool": s["tool"], "args": args})
    return out, params


def _yaml(obj, indent=0) -> str:
    """Minimal YAML emitter — avoids a PyYAML dependency for writing."""
    pad = "  " * indent
    if isinstance(obj, dict):
        lines = []
        for k, v in obj.items():
            if isinstance(v, (dict, list)) and v:
                lines.append(f"{pad}{k}:")
                lines.append(_yaml(v, indent + 1))
            else:
                lines.append(f"{pad}{k}: {_scalar(v)}")
        return "\n".join(lines)
    if isinstance(obj, list):
        lines = []
        for item in obj:
            if isinstance(item, dict):
                body = _yaml(item, indent + 1)
                first, *rest = body.split("\n")
                lines.append(f"{pad}- {first.strip()}")
                lines.extend(rest)
            else:
                lines.append(f"{pad}- {_scalar(item)}")
        return "\n".join(lines)
    return f"{pad}{_scalar(obj)}"


def _scalar(v) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    s = str(v)
    if s == "" or re.search(r'[:#\[\]{}",\n]|^\s|\s$|^\$', s):
        return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'
    return s


def freeze_run(run: str, *, out: "str | None" = None,
               study: "str | None" = None,
               parameterise: bool = True) -> Path:
    data = load_run(run)
    transcript = data.get("transcript") or []
    steps = _steps_from(transcript)
    if not steps:
        raise SystemExit(
            "error: this run executed no tool calls, so there is nothing to "
            "pin. A playbook freezes a sequence of tool invocations; a session "
            "that only produced prose has none.")

    query = data.get("query") or ""
    params: dict = {}
    if parameterise:
        steps, params = _parameterise(steps, query)

    hashes = _artefact_hashes(transcript)
    doc: dict = {
        "study": study or (query[:90] or "Frozen session"),
        "description": (f"Frozen from a recorded session on "
                        f"{time.strftime('%Y-%m-%d')}. Original question: "
                        f"{query[:200]}"),
    }
    if params:
        doc["parameters"] = params
    doc["steps"] = steps
    doc["synthesis"] = (
        "Summarise the findings from the captured artefacts. State what was "
        "measured and what was not.")
    doc["provenance"] = {
        "frozen_from": str(run),
        "frozen_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "original_backend": data.get("backend"),
        "original_model": data.get("model"),
        "original_stop_reason": data.get("stop_reason"),
        "consistency": data.get("consistency") or {},
    }
    # Expected artefacts, so a replay can detect upstream drift rather than
    # silently producing different data under the same commands.
    doc["expected_artefacts"] = hashes

    path = Path(out) if out else (PLAYBOOK_DIR /
        f"frozen_{re.sub(r'[^a-z0-9]+', '_', (query[:40] or 'session').lower()).strip('_')}.yaml")
    if not path.is_absolute():
        path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_yaml(doc) + "\n")
    return path


def cmd_list(args) -> int:
    if not AGENT_DIR.is_dir():
        print("No recorded sessions yet (Docs/Agent/).")
        return 0
    rows = []
    for d in sorted(AGENT_DIR.iterdir(), reverse=True):
        t = d / "transcript.json"
        if not t.is_file():
            continue
        try:
            data = json.loads(t.read_text())
        except Exception:
            continue
        n = len(_steps_from(data.get("transcript") or []))
        rows.append((d.name, n, (data.get("query") or "")[:52]))
    print(f"{'RUN':<52} {'STEPS':>5}  QUERY")
    for name, n, q in rows[: args.limit]:
        print(f"{name:<52} {n:>5}  {q}")
    print(f"\n{len(rows)} recorded session(s); "
          f"{sum(1 for _, n, _ in rows if n)} have tool calls to freeze")
    return 0


def cmd_freeze(args) -> int:
    path = freeze_run(args.run, out=args.out, study=args.study,
                      parameterise=not args.no_parameters)
    # Count from the structure, not by grepping the YAML we just wrote —
    # the emitter indents list items, so the naive pattern matched nothing.
    data = load_run(args.run)
    n_steps = len(_steps_from(data.get("transcript") or []))
    n_hash = sum(1 for h in _artefact_hashes(data.get("transcript") or [])
                 if h.get("sha256"))
    print(f"Wrote:         {path}")
    print(f"  steps pinned:        {n_steps}")
    print(f"  artefacts hashed:    {n_hash}")
    print(f"\nReplay deterministically with:\n  igvfagent playbook run "
          f"--file {path.relative_to(ROOT) if path.is_relative_to(ROOT) else path}")
    return 0


def cmd_check(args) -> int:
    """Re-hash a playbook's expected artefacts and report drift."""
    p = Path(args.file)
    if not p.is_absolute():
        p = ROOT / p
    text = p.read_text()
    entries = re.findall(r"- path:\s*(\S+)\s*\n\s*sha256:\s*(\S+)", text)
    if not entries:
        print("error: no expected_artefacts in this playbook", file=sys.stderr)
        return 2
    ok = missing = drifted = 0
    for rel, want in entries:
        want = want.strip('"')
        f = ROOT / rel.strip('"')
        got = _sha256(f)
        if got is None:
            print(f"  MISSING  {rel}")
            missing += 1
        elif want in ("null", "None"):
            print(f"  unhashed {rel}")
        elif got != want:
            print(f"  DRIFTED  {rel}")
            print(f"           expected {want[:16]}…  got {got[:16]}…")
            drifted += 1
        else:
            ok += 1
    print(f"\n  {ok} unchanged · {drifted} drifted · {missing} missing")
    if drifted:
        print("  Drift usually means the upstream archive re-released its "
              "data, not that the agent changed. Compare retrieval dates "
              "before assuming a regression.")
    return 1 if (drifted or missing) else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent playbook-freeze",
        description="Convert a recorded session into a deterministic playbook "
                    "with pinned arguments and expected artefact hashes.")
    sub = p.add_subparsers(dest="cmd", required=True)

    l = sub.add_parser("list", help="Recorded sessions available to freeze")
    l.add_argument("--limit", type=int, default=25)

    f = sub.add_parser("freeze", help="Freeze one session into a playbook")
    f.add_argument("--run", required=True,
                   help="Docs/Agent/<run>/ or its transcript.json")
    f.add_argument("--out", help="Output YAML path")
    f.add_argument("--study", help="Playbook title")
    f.add_argument("--no-parameters", action="store_true",
                   help="Pin literal values instead of lifting parameters")

    c = sub.add_parser("check", help="Verify artefacts against a frozen playbook")
    c.add_argument("--file", required=True)

    args = p.parse_args(argv)
    return {"list": cmd_list, "freeze": cmd_freeze, "check": cmd_check}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
