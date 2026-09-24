"""Review queue for agent-authored extensions, and their promotion to built-ins.

Agents on a deployment can author tools and skills (IGVF_ALLOW_AGENT_AUTHORING)
into UserExtensions/. They work there, but they are unreviewed: one of them
read its content from stdin and hung the MCP server of two hosted jobs. This
command is how that code becomes part of IGVF Agent's core, deliberately:

    igvfagent ext-review list [--json]        every extension: kind, what it is,
                                              usage and outcomes in job logs,
                                              safety flags, risk, verdict hint
    igvfagent ext-review show <name>          one extension, with its source
    igvfagent ext-review retire <name> --reason TEXT
                                              move it to _retired_<date>/ (reversible)
    igvfagent ext-review bundle <name>        tar.gz of manifest + module + review
                                              (Docs/Extensions/bundles/), to carry
                                              it from a deployment to a dev checkout
    igvfagent ext-review promote <name|bundle.tar.gz> --reviewer NAME
             [--accept-risk "why"]            copy into Scripts/promoted/ as a
                                              built-in (ships with the code, shadows
                                              the extension), with provenance
    igvfagent ext-review selftest

Promotion is refused when the static scan finds a high-risk pattern (stdin,
eval/exec, a shell, secret access) unless the reviewer accepts it with a
reason, and refused on a public deployment (it must be done in a checkout and
committed). Scripts/test_promoted.py keeps every promoted extension honest:
its manifest loads, its module imports, --help works, and no new high-risk
pattern appeared after review.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tarfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
try:
    from igvfagent import _userext  # type: ignore
except Exception:
    import _userext  # type: ignore

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or HERE.parent).resolve()
JOBS_DIR = Path(os.environ.get("IGVF_JOBS_DIR") or ROOT / "Data" / "Jobs")
REVIEW_DIR = ROOT / "Docs" / "Extensions"
PROMOTED = HERE / "promoted"

# (flag, regex, level). "high" blocks promotion without an accepted reason.
SAFETY = [
    ("reads_stdin", r"sys\.stdin|(?<![\w.])input\s*\(|fileinput", "high"),
    ("eval_exec", r"(?<![\w.])(eval|exec)\s*\(|__import__\s*\(|compile\s*\(", "high"),
    ("shell", r"os\.system|shell\s*=\s*True|os\.popen|pty\.spawn", "high"),
    ("secrets", r"Docs/Secret|/Secret/|API[_-]?KEY|_TOKEN|password|credential", "high"),
    ("subprocess", r"subprocess\.|os\.exec|os\.spawn", "medium"),
    ("deletes", r"shutil\.rmtree|os\.remove|os\.unlink|\.unlink\(|os\.rmdir", "medium"),
    ("writes_files", r"open\([^)]*['\"][wa]b?\+?['\"]|\.write_text\(|\.write_bytes\(|shutil\.(copy|move)", "low"),
    ("network", r"requests\.|urllib\.request|http\.client|socket\.|httpx", "low"),
    ("absolute_paths", r"['\"]/(etc|root|home|usr|var|opt|tmp)/", "medium"),
]
_RANK = {"low": 1, "medium": 2, "high": 3}


def scan(text: str) -> "List[Dict[str, str]]":
    out = []
    for flag, rx, level in SAFETY:
        m = re.search(rx, text or "")
        if m:
            line = text[:m.start()].count("\n") + 1
            out.append({"flag": flag, "level": level, "line": line,
                        "match": text.splitlines()[line - 1].strip()[:140] if text else ""})
    return out


def risk(flags: "List[Dict[str, str]]") -> str:
    return max((f["level"] for f in flags), key=lambda l: _RANK[l], default="none")


def _builtin_names() -> "set[str]":
    try:
        try:
            from igvfagent import cli  # type: ignore
        except Exception:
            import cli  # type: ignore
        skills = set(cli.SKILLS)
    except Exception:  # noqa: BLE001
        skills = set()
    try:
        import _tools
        tools = {t.name for t in _tools._TOOLS if t.name not in _tools._USER_TOOL_NAMES}
    except Exception:  # noqa: BLE001
        tools = set()
    return skills | tools


def _usage() -> "Dict[str, Dict[str, int]]":
    """Calls and outcomes per tool name, from every job's event log."""
    use: Dict[str, Dict[str, int]] = {}
    if not JOBS_DIR.is_dir():
        return use
    for ev in JOBS_DIR.glob("*/events.jsonl"):
        job = ev.parent.name
        try:
            lines = ev.read_text(errors="replace").splitlines()
        except OSError:
            continue
        for ln in lines:
            try:
                e = json.loads(ln)
            except ValueError:
                continue
            name = e.get("text") if e.get("kind") in ("tool", "tool_end") else None
            if not name:
                continue
            u = use.setdefault(name, {"calls": 0, "ok": 0, "failed": 0, "jobs": set()})  # type: ignore[dict-item]
            if e["kind"] == "tool":
                u["calls"] += 1
                u["jobs"].add(job)  # type: ignore[union-attr]
            elif "exit" in e:
                u["ok" if e.get("exit") in (0, "0") else "failed"] += 1
    for u in use.values():
        u["jobs"] = len(u["jobs"])  # type: ignore[arg-type]
    return use


def inventory() -> "List[Dict[str, Any]]":
    builtins = _builtin_names()
    use = _usage()
    tools = {t["name"]: t for t in _userext.discover_tools()}
    skills = _userext.discover_skills()
    rows: List[Dict[str, Any]] = []
    linked_skills = set()
    for name, t in sorted(tools.items()):
        module = None
        if t.get("cli"):
            sk = skills.get(t["cli"][0])
            if sk:
                module = sk["path"]
                linked_skills.add(t["cli"][0])
        src = Path(module).read_text(errors="replace") if module else " ".join(t.get("command") or [])
        flags = scan(src)
        rows.append({"name": name, "kind": "tool", "description": t["description"][:240], "manifest": t.get("source"),
                     "module": module, "command": t.get("command") or None, "cli": t.get("cli") or None,
                     "shadowed_by_builtin": name in builtins, "usage": use.get(name, {"calls": 0, "ok": 0, "failed": 0,
                                                                                      "jobs": 0}),
                     "flags": flags, "risk": risk(flags),
                     "sha256": hashlib.sha256(src.encode()).hexdigest()[:16],
                     "modified": time.strftime("%Y-%m-%d", time.gmtime(Path(t["source"]).stat().st_mtime))
                     if t.get("source") and Path(t["source"]).exists() else None})
    for name, sk in sorted(skills.items()):
        if name in linked_skills:
            continue
        src = Path(sk["path"]).read_text(errors="replace")
        flags = scan(src)
        rows.append({"name": name, "kind": "skill", "description": sk["description"][:240], "manifest": None,
                     "module": sk["path"], "command": None, "cli": [name], "shadowed_by_builtin": name in builtins,
                     "usage": use.get(name.replace("-", "_"), {"calls": 0, "ok": 0, "failed": 0, "jobs": 0}),
                     "flags": flags, "risk": risk(flags), "sha256": hashlib.sha256(src.encode()).hexdigest()[:16],
                     "modified": time.strftime("%Y-%m-%d", time.gmtime(Path(sk["path"]).stat().st_mtime))})
    for r in rows:
        u = r["usage"]
        if r["shadowed_by_builtin"]:
            r["hint"] = "retire (a built-in does this)"
        elif r["risk"] == "high":
            r["hint"] = "fix or retire (high-risk pattern)"
        elif u["calls"] >= 3 and u["failed"] <= u["ok"]:
            r["hint"] = "promotion candidate"
        elif u["calls"] == 0:
            r["hint"] = "unused"
        else:
            r["hint"] = "keep as extension"
    return rows


def _find(name: str) -> Optional[Dict[str, Any]]:
    key = name.replace("-", "_")
    return next((r for r in inventory() if r["name"].replace("-", "_") == key), None)


def write_review() -> "Dict[str, Path]":
    rows = inventory()
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    js = REVIEW_DIR / f"review_{stamp}.json"
    js.write_text(json.dumps(rows, indent=2, default=str))
    L = ["# Extension review", "", f"{len(rows)} agent-authored extension(s) on this deployment.", "",
         "| name | kind | risk | calls (ok/failed, jobs) | hint | flags |", "|---|---|---|---|---|---|"]
    order = {"promotion candidate": 0, "fix or retire (high-risk pattern)": 1, "retire (a built-in does this)": 2,
             "keep as extension": 3, "unused": 4}
    for r in sorted(rows, key=lambda r: (order.get(r["hint"], 9), -r["usage"]["calls"])):
        u = r["usage"]
        L.append(f"| `{r['name']}` | {r['kind']} | {r['risk']} | {u['calls']} ({u['ok']}/{u['failed']}, {u['jobs']}) "
                 f"| {r['hint']} | {', '.join(f['flag'] for f in r['flags'])} |")
    md = REVIEW_DIR / f"review_{stamp}.md"
    md.write_text("\n".join(L) + "\n")
    return {"json": js, "md": md}


def retire(name: str, reason: str) -> Dict[str, Any]:
    r = _find(name)
    if not r:
        return {"ok": False, "error": f"no extension {name}"}
    moved = []
    for p in (r.get("manifest"), r.get("module")):
        if p and Path(p).exists():
            dest = Path(p).parents[1] / f"_retired_{time.strftime('%Y%m%d')}" / Path(p).name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(p, dest)
            moved.append(str(dest))
    if moved:
        note = Path(moved[0]).parent / "RETIRED.md"
        with open(note, "a") as fh:
            fh.write(f"- {time.strftime('%Y-%m-%d %H:%M')} `{r['name']}` ({r['kind']}): {reason} "
                     f"[risk {r['risk']}; flags {', '.join(f['flag'] for f in r['flags']) or 'none'}]\n")
    return {"ok": bool(moved), "moved": moved}


def bundle(name: str) -> Dict[str, Any]:
    r = _find(name)
    if not r:
        return {"ok": False, "error": f"no extension {name}"}
    out = REVIEW_DIR / "bundles" / f"{r['name']}_{r['sha256']}.tar.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(out, "w:gz") as tf:
        for p, arc in ((r.get("manifest"), "tool.json"), (r.get("module"), "skill.py")):
            if p and Path(p).exists():
                tf.add(p, arcname=arc)
        data = json.dumps({k: v for k, v in r.items()}, indent=2, default=str).encode()
        info = tarfile.TarInfo("review.json")
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))
    return {"ok": True, "bundle": str(out)}


def promote(target: str, reviewer: str, accept_risk: str = "") -> Dict[str, Any]:
    if os.environ.get("IGVF_PUBLIC") == "1" or os.environ.get("IGVF_HOSTED") == "1":
        return {"ok": False, "error": "promote in a development checkout and commit it; on a deployment use "
                                      "`ext-review bundle` and copy the bundle over"}
    if not reviewer.strip():
        return {"ok": False, "error": "--reviewer is required: promotion is a reviewed decision"}
    manifest_txt = module_txt = None
    review: Dict[str, Any] = {}
    if target.endswith(".tar.gz") and Path(target).is_file():
        with tarfile.open(target) as tf:
            names = tf.getnames()
            if "tool.json" in names:
                manifest_txt = tf.extractfile("tool.json").read().decode()  # type: ignore[union-attr]
            if "skill.py" in names:
                module_txt = tf.extractfile("skill.py").read().decode()  # type: ignore[union-attr]
            if "review.json" in names:
                review = json.loads(tf.extractfile("review.json").read())  # type: ignore[union-attr]
        name = review.get("name") or Path(target).name.split("_")[0]
    else:
        r = _find(target)
        if not r:
            return {"ok": False, "error": f"no extension or bundle {target}"}
        review, name = r, r["name"]
        if r.get("manifest"):
            manifest_txt = Path(r["manifest"]).read_text()
        if r.get("module"):
            module_txt = Path(r["module"]).read_text()
    if name in _builtin_names() and not (PROMOTED / "tools" / f"{name}.json").exists():
        return {"ok": False, "error": f"`{name}` is already a built-in; retire the extension instead"}
    flags = scan(module_txt or "")
    if risk(flags) == "high" and not accept_risk.strip():
        return {"ok": False, "error": "high-risk pattern(s): " + "; ".join(
            f"{f['flag']} (line {f['line']}: {f['match']})" for f in flags if f["level"] == "high")
                + ". Fix the code, or pass --accept-risk with the reason it is safe."}
    written = []
    prov = {"promoted_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "reviewer": reviewer,
            "source_sha256": hashlib.sha256((module_txt or manifest_txt or "").encode()).hexdigest(),
            "usage_at_review": review.get("usage"), "flags_at_review": [f["flag"] for f in flags],
            "accepted_risk": accept_risk or None}
    if module_txt is not None:
        stem = (json.loads(manifest_txt)["cli"][0] if manifest_txt and json.loads(manifest_txt).get("cli")
                else name).replace("-", "_")
        dst = PROMOTED / "skills" / f"{stem}.py"
        dst.parent.mkdir(parents=True, exist_ok=True)
        header = (f'# Promoted from an agent-authored extension on {prov["promoted_at"][:10]}; reviewed by '
                  f'{reviewer}. Provenance: Scripts/promoted/registry.json.\n')
        dst.write_text(header + module_txt if not module_txt.startswith("# Promoted") else module_txt)
        written.append(str(dst))
    if manifest_txt is not None:
        m = json.loads(manifest_txt)
        m["x-promoted"] = prov
        dst = PROMOTED / "tools" / f"{name}.json"
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(json.dumps(m, indent=2) + "\n")
        written.append(str(dst))
    reg_p = PROMOTED / "registry.json"
    reg = json.loads(reg_p.read_text()) if reg_p.exists() else {}
    reg[name] = {**prov, "files": [os.path.relpath(w, PROMOTED.parent) for w in written]}
    reg_p.write_text(json.dumps(reg, indent=2, sort_keys=True) + "\n")
    return {"ok": True, "written": written,
            "next": "run `python3 Scripts/test_promoted.py`, add a line to Docs/Guide/whats-new.md, commit, deploy; "
                    "then retire the extension on the deployment"}


def selftest() -> int:
    import tempfile
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")

    bad = "import sys\ncontent = sys.stdin.read()\nopen('/etc/x','w')\n"
    fl = {f["flag"] for f in scan(bad)}
    check("stdin, file writes and absolute paths are flagged", {"reads_stdin", "writes_files", "absolute_paths"} <= fl
          and risk(scan(bad)) == "high")
    check("a plain analysis is low risk", risk(scan("import pandas as pd\ndf = pd.read_csv(p)\nprint(df.shape)")) == "none")
    global PROMOTED, JOBS_DIR, REVIEW_DIR
    saved = (PROMOTED, JOBS_DIR, REVIEW_DIR)
    with tempfile.TemporaryDirectory() as td:
        ext = Path(td) / "ext"
        (ext / "tools").mkdir(parents=True)
        (ext / "skills").mkdir()
        (ext / "tools" / "gc_calc.json").write_text(json.dumps({
            "name": "gc_calc", "description": "GC content of a sequence",
            "parameters": {"type": "object", "properties": {"seq": {"type": "string"}}},
            "cli": ["gc-calc"], "flag_map": {"seq": "--seq"}}))
        (ext / "skills" / "gc_calc.py").write_text(
            '"""GC content."""\nimport argparse\n\ndef main(argv=None):\n    a = argparse.ArgumentParser()\n'
            '    a.add_argument("--seq", required=True)\n    s = a.parse_args(argv).seq.upper()\n'
            '    print(round((s.count("G") + s.count("C")) / max(len(s), 1), 3))\n    return 0\n')
        (ext / "skills" / "eater.py").write_text('"""bad"""\nimport sys\ndef main(argv=None):\n    sys.stdin.read()\n')
        jobs = Path(td) / "jobs" / "J1"
        jobs.mkdir(parents=True)
        (jobs / "events.jsonl").write_text("".join(json.dumps(e) + "\n" for e in [
            {"kind": "tool", "text": "gc_calc"}, {"kind": "tool_end", "text": "gc_calc", "exit": 0}] * 3))
        os.environ["IGVF_USER_EXT_DIR"] = str(ext)
        PROMOTED, JOBS_DIR, REVIEW_DIR = Path(td) / "promoted", Path(td) / "jobs", Path(td) / "review"
        try:
            rows = {r["name"]: r for r in inventory()}
            check("inventory pairs a tool with its module and counts its calls",
                  rows["gc_calc"]["module"] and rows["gc_calc"]["usage"]["calls"] == 3
                  and rows["gc_calc"]["hint"] == "promotion candidate")
            check("the stdin reader is high risk", rows["eater"]["risk"] == "high")
            r = promote("eater", "tester")
            check("promotion refuses a high-risk extension without an accepted reason",
                  not r["ok"] and "reads_stdin" in r["error"])
            check("promotion requires a reviewer", not promote("gc_calc", "")["ok"])
            b = bundle("gc_calc")
            r = promote(b["bundle"], "tester")
            reg = json.loads((PROMOTED / "registry.json").read_text())
            check("a bundle promotes into promoted/ with provenance",
                  r["ok"] and (PROMOTED / "tools" / "gc_calc.json").exists()
                  and (PROMOTED / "skills" / "gc_calc.py").exists() and reg["gc_calc"]["reviewer"] == "tester")
            rr = retire("eater", "reads stdin; replaced by built-in files")
            check("retire moves the extension aside (reversible)", rr["ok"] and not (ext / "skills" / "eater.py").exists()
                  and any("_retired_" in m for m in rr["moved"]))
        finally:
            PROMOTED, JOBS_DIR, REVIEW_DIR = saved
            os.environ.pop("IGVF_USER_EXT_DIR", None)
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent ext-review", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list")
    s.add_argument("--json", action="store_true")
    s = sub.add_parser("show")
    s.add_argument("name")
    s = sub.add_parser("retire")
    s.add_argument("name")
    s.add_argument("--reason", required=True)
    s = sub.add_parser("bundle")
    s.add_argument("name")
    s = sub.add_parser("promote")
    s.add_argument("target")
    s.add_argument("--reviewer", required=True)
    s.add_argument("--accept-risk", default="")
    sub.add_parser("selftest")
    a = ap.parse_args(argv)
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "list":
        if a.json:
            print(json.dumps(inventory(), indent=2, default=str))
            return 0
        out = write_review()
        for ln in out["md"].read_text().splitlines()[4:]:
            print(ln)
        print(f"\nReport: {out['md'].relative_to(ROOT) if out['md'].is_relative_to(ROOT) else out['md']}")
        return 0
    if a.cmd == "show":
        r = _find(a.name)
        if not r:
            print(f"no extension {a.name}")
            return 1
        print(json.dumps(r, indent=2, default=str))
        if r.get("module"):
            print("\n--- source ---\n" + Path(r["module"]).read_text()[:6000])
        return 0
    res = {"retire": lambda: retire(a.name, a.reason), "bundle": lambda: bundle(a.name),
           "promote": lambda: promote(a.target, a.reviewer, a.accept_risk)}[a.cmd]()
    print(json.dumps(res, indent=2))
    return 0 if res.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
