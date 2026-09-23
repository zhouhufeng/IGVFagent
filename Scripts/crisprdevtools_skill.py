#!/usr/bin/env python3
"""Nextflow module scaffolding and layout checks for IGVF CRISPR pipeline modules (port of IGVF-CRISPR/crisprdevtools).

Port of https://github.com/IGVF-CRISPR/crisprdevtools (branch main, commit
08dd9ae3356dcaf5f2921b729a72e329eedf9f89, 2024-07-10).  No LICENSE file.  The
upstream is a single 57-line function, create_new_module_nextflow(name), whose
README states the purpose "Tools to help developing modules".  Its behaviour
(the directory layout, file names and file contents it writes) was read and
re-derived here; no code was copied.  Relationship: port.  Because the
upstream is so small, the stated purpose is also implemented faithfully
beyond it: an optional fuller template matching the module layout used by
IGVF-CRISPR/IGVF_CRISPR_Pipeline (Modules/<name>/{main.nf, input.config,
processes/<name>.nf, bin/, conda_envs/<name>.yaml, example_data/, test/}),
and a layout checker for existing modules.

Definitions, exactly as upstream creates them
  directories     <name>/bin, <name>/conda_envs, <name>/example_data,
                  <name>/processes, <name>/test
  README.md       empty
  input.config    empty
  main.nf         "nextflow.enable.dsl=2\\n"
  bin/<name>.py   "#!/usr/bin/env python\\n"
  conda env       conda_envs/<name>.yaml: name <name>; channels conda-forge,
                  bioconda, defaults; dependencies numpy=1.23.5,
                  pandas=2.2.1, pip, pip: openpyxl==3.1.2
  process         processes/<name>.nf: a process with
                  conda "${moduleDir}/conda_envs/<name>.yaml", empty
                  input/output blocks, script calling <name>.py

Upstream bugs corrected by default (reproduce with --upstream-compat)
  * the conda YAML is written from an indented f-string, so every line after
    "name:" carries two extra spaces and the file does not parse as YAML;
    the port writes it flush-left.
  * the process is always named seqSpecParser (copied from the first
    module); the port names it after the module.
  * bin/<name>.py is not made executable, so Nextflow cannot call it from
    bin/; the port sets mode 755.

Subcommands
  new-module    scaffold <dest>/<name> (template minimal = upstream layout
                and contents; template full adds an include + workflow in
                main.nf, a params block in input.config, test/test.nf, an
                argparse skeleton in bin/ and README sections).
  check-module  validate an existing module directory: required dirs,
                DSL2 header, process declarations, conda env references
                resolve and parse, bin scripts have a shebang and are
                executable, includes resolve, input.config params used.
  selftest      scaffolds in a tempdir (both templates, compat on/off) and
                checks them, plus a hand-broken module.

Output: Docs/CRISPRDevTools/<timestamp>_<label>/ (report.md, summary.json;
new-module writes the module there unless --dest is given).  Standard
library only (PyYAML used for YAML checks when installed).

Usage:
    igvfagent crisprdevtools new-module --name guide_assignment --dest Modules
    igvfagent crisprdevtools new-module --name guide_assignment --template full
    igvfagent crisprdevtools check-module --path Modules/seqSpecParser
    igvfagent crisprdevtools selftest
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import stat
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRDevTools"

UPSTREAM_REPO = "IGVF-CRISPR/crisprdevtools"
UPSTREAM_COMMIT = "08dd9ae3356dcaf5f2921b729a72e329eedf9f89"
MODULE_DIRS = ["bin", "conda_envs", "example_data", "processes", "test"]
REQUIRED_DIRS = ["bin", "conda_envs", "processes"]
CONDA_CHANNELS = ["conda-forge", "bioconda", "defaults"]
CONDA_DEPS = ["numpy=1.23.5", "pandas=2.2.1", "pip"]
CONDA_PIP = ["openpyxl==3.1.2"]
NAME_RX = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

log = logging.getLogger("crisprdevtools")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crisprdevtools_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    base = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d, i = base, 1
    while d.exists():
        i += 1
        d = Path(f"{base}_{i}")
    d.mkdir(parents=True, exist_ok=True)
    return d


def md_table(headers: "List[str]", rows: "Iterable[Iterable[Any]]") -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for r in rows:
        out.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |")
    return "\n".join(out)


def finish(d: Path, title: str, sections: "List[str]", summary: dict) -> None:
    rep = d / "report.md"
    rep.write_text("\n".join([f"# {title}", "", f"Port of https://github.com/{UPSTREAM_REPO} (commit {UPSTREAM_COMMIT[:10]}).", ""]
                             + sections) + "\n")
    print(f"Report: {rep}")
    js = d / "summary.json"
    js.write_text(json.dumps(summary, indent=2, default=str))
    print(f"JSON: {js}")


# ---------------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------------

def process_name(module: str, upstream_compat: bool = False) -> str:
    return "seqSpecParser" if upstream_compat else module


def conda_yaml(module: str, upstream_compat: bool = False) -> str:
    lines = [f"name: {module}", "channels:"] + [f"  - {c}" for c in CONDA_CHANNELS] + ["dependencies:"] \
        + [f"  - {x}" for x in CONDA_DEPS] + ["  - pip:"] + [f"    - {x}" for x in CONDA_PIP]
    if upstream_compat:  # the indented f-string: two leading spaces on every line after the first, and a trailing "  "
        return lines[0] + "\n" + "\n".join("  " + l for l in lines[1:]) + "\n  "
    return "\n".join(lines) + "\n"


def process_nf(module: str, upstream_compat: bool = False, template: str = "minimal") -> str:
    name = process_name(module, upstream_compat)
    if template == "full" and not upstream_compat:
        return (f'process {name} {{\n'
                f'    conda "${{moduleDir}}/conda_envs/{module}.yaml"\n\n'
                f'    input:\n        path input_file\n\n'
                f'    output:\n        path "{module}_output.txt", emit: {module}_out\n\n'
                f'    script:\n        """\n        {module}.py --input ${{input_file}} --output {module}_output.txt\n        """\n}}\n')
    body = (f'process {name} {{\n  conda "${{moduleDir}}/conda_envs/{module}.yaml"\n  input:\n  output:\n  script:\n'
            f'        """\n        {module}.py\n        """\n}}\n')
    if upstream_compat:  # the upstream f-string indents every line after the first by two spaces
        ls = body.rstrip("\n").split("\n")
        return ls[0] + "\n" + "\n".join("  " + l for l in ls[1:]) + "\n  "
    return body


def main_nf(module: str, template: str, upstream_compat: bool = False) -> str:
    if template != "full":
        return "nextflow.enable.dsl=2\n"
    name = process_name(module, upstream_compat)
    return (f"nextflow.enable.dsl=2\n\ninclude {{ {name} }} from './processes/{module}'\n\n"
            f"workflow {{\n    {module}_ch = {name}(file(params.{module}_input))\n}}\n")


def input_config(module: str, template: str) -> str:
    if template != "full":
        return ""
    return f'params.{module}_input = "example_data/{module}_example.txt"\n'


def bin_script(module: str, template: str) -> str:
    if template != "full":
        return "#!/usr/bin/env python\n"
    return ("#!/usr/bin/env python\n"
            f'"""{module}: one step of the IGVF CRISPR pipeline."""\n'
            "import argparse\n\n\n"
            "def main():\n"
            "    ap = argparse.ArgumentParser()\n"
            "    ap.add_argument('--input', required=True)\n"
            "    ap.add_argument('--output', required=True)\n"
            "    a = ap.parse_args()\n"
            "    with open(a.input) as fi, open(a.output, 'w') as fo:\n"
            "        fo.write(fi.read())\n\n\n"
            "if __name__ == '__main__':\n"
            "    main()\n")


def readme(module: str, template: str) -> str:
    if template != "full":
        return ""
    return (f"# {module}\n\n## Inputs\n\n- `params.{module}_input`: ...\n\n## Outputs\n\n- `{module}_output.txt`\n\n"
            f"## Run\n\n```bash\nnextflow run main.nf -c input.config\n```\n\n## Test\n\n```bash\nnextflow run test/test.nf -c input.config\n```\n")


def test_nf(module: str, upstream_compat: bool = False) -> str:
    name = process_name(module, upstream_compat)
    return (f"nextflow.enable.dsl=2\n\ninclude {{ {name} }} from '../processes/{module}'\n\n"
            f"workflow {{\n    out = {name}(file(\"${{projectDir}}/../example_data/{module}_example.txt\"))\n"
            f"    out.view {{ f -> assert f.exists(); \"ok: ${{f}}\" }}\n}}\n")


def create_new_module(module: str, dest: Path, template: str = "minimal", upstream_compat: bool = False,
                      force: bool = False) -> "Tuple[Path, List[str]]":
    """create_new_module_nextflow(): the upstream layout (+ the full template when asked)."""
    if not NAME_RX.match(module):
        raise ValueError(f"module name {module!r} must start with a letter and contain only letters, digits and _")
    root = Path(dest) / module
    if root.exists() and any(root.iterdir()) and not force:
        raise FileExistsError(f"{root} exists and is not empty (use --force to overwrite files)")
    written: List[str] = []
    for sd in MODULE_DIRS:
        (root / sd).mkdir(parents=True, exist_ok=True)

    def w(rel: str, text: str) -> None:
        p = root / rel
        p.write_text(text)
        written.append(rel)
    w("README.md", readme(module, template))
    w("input.config", input_config(module, template))
    w("main.nf", main_nf(module, template, upstream_compat))
    w(f"bin/{module}.py", bin_script(module, template))
    if not upstream_compat:
        p = root / "bin" / f"{module}.py"
        p.chmod(p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IRGRP | stat.S_IROTH)
    w(f"conda_envs/{module}.yaml", conda_yaml(module, upstream_compat))
    w(f"processes/{module}.nf", process_nf(module, upstream_compat, template))
    if template == "full":
        w("test/test.nf", test_nf(module, upstream_compat))
        w(f"example_data/{module}_example.txt", "example input\n")
    return root, written


# ---------------------------------------------------------------------------
# YAML (subset) and module checks
# ---------------------------------------------------------------------------

def parse_simple_yaml(text: str) -> Any:
    """PyYAML when installed; else a strict parser for the block-mapping / block-sequence subset conda envs use."""
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text)
    except ImportError:
        pass
    lines = [(len(l) - len(l.lstrip(" ")), l.strip()) for l in text.split("\n")
             if l.strip() and not l.strip().startswith("#")]
    pos = 0

    def scalar(v: str) -> Any:
        v = v.strip()
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "'\"":
            return v[1:-1]
        return v

    def block(indent: int) -> Any:
        nonlocal pos
        if pos >= len(lines):
            return None
        ind, s = lines[pos]
        if ind != indent:
            raise ValueError(f"bad indentation at {s!r} (expected {indent}, got {ind})")
        if s.startswith("- "):
            out = []
            while pos < len(lines) and lines[pos][0] == indent and lines[pos][1].startswith("- "):
                item = lines[pos][1][2:].strip()
                pos += 1
                if item.endswith(":") and ": " not in item:
                    nxt = lines[pos][0] if pos < len(lines) else -1
                    out.append({item[:-1]: block(nxt) if nxt > indent else None})
                elif ": " in item:
                    k, v = item.split(": ", 1)
                    out.append({k: scalar(v)})
                else:
                    out.append(scalar(item))
            if pos < len(lines) and lines[pos][0] > indent:
                raise ValueError(f"bad indentation at {lines[pos][1]!r}")
            return out
        out: Dict[str, Any] = {}
        while pos < len(lines) and lines[pos][0] == indent:
            s = lines[pos][1]
            if s.startswith("- "):
                raise ValueError(f"sequence item {s!r} inside a mapping")
            if s.endswith(":") and ": " not in s:
                pos += 1
                nxt = lines[pos][0] if pos < len(lines) else -1
                out[s[:-1]] = block(nxt) if nxt > indent or (nxt == indent and pos < len(lines) and lines[pos][1].startswith("- ")) else None
            elif ": " in s:
                k, v = s.split(": ", 1)
                if ": " in v:
                    raise ValueError(f"mapping values are not allowed here: {s!r}")
                pos += 1
                out[k] = scalar(v)
            else:
                raise ValueError(f"cannot parse line {s!r}")
        if pos < len(lines) and lines[pos][0] > indent:
            raise ValueError(f"bad indentation at {lines[pos][1]!r} (a continuation of a scalar)")
        return out
    res = block(lines[0][0]) if lines else None
    if pos < len(lines):
        raise ValueError(f"unparsed content from {lines[pos][1]!r}")
    return res


def check_module(path: Path) -> "List[Tuple[str, str, str]]":
    """-> list of (status PASS/WARN/FAIL, check, detail)."""
    path = Path(path)
    out: List[Tuple[str, str, str]] = []

    def add(ok, name, detail="", warn=False):
        out.append(("PASS" if ok else ("WARN" if warn else "FAIL"), name, detail))
    add(path.is_dir(), "module directory exists", str(path))
    if not path.is_dir():
        return out
    for sd in MODULE_DIRS:
        add((path / sd).is_dir(), f"{sd}/ present", "" if sd in REQUIRED_DIRS else "recommended", warn=sd not in REQUIRED_DIRS)
    mn = path / "main.nf"
    txt = mn.read_text() if mn.is_file() else ""
    add(mn.is_file(), "main.nf present")
    add("nextflow.enable.dsl=2" in txt.replace(" ", ""), "main.nf enables DSL2")
    add((path / "input.config").is_file(), "input.config present", warn=True)
    procs = {}
    for nf in sorted((path / "processes").glob("*.nf")) if (path / "processes").is_dir() else []:
        t = nf.read_text()
        names = re.findall(r"^\s*process\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", t, flags=re.M)
        add(bool(names), f"processes/{nf.name} declares a process", ", ".join(names))
        for n in names:
            procs[n] = nf
        if names and nf.stem not in names:
            add(False, f"processes/{nf.name}: process name matches the file", f"declares {names}", warn=True)
        for env in re.findall(r'conda\s+["\']([^"\']+)["\']', t):
            rel = env.replace("${moduleDir}/", "").replace("$moduleDir/", "")
            cand = [path / rel, nf.parent / rel, path / rel.lstrip("./")]
            hit = next((c for c in cand if c.is_file()), None)
            add(hit is not None, f"processes/{nf.name}: conda env resolves", env)
            if hit is not None:
                try:
                    y = parse_simple_yaml(hit.read_text())
                    good = isinstance(y, dict) and {"name", "dependencies"} <= set(y)
                    add(good, f"{hit.relative_to(path)} parses as a conda env", "keys: " + ", ".join(sorted(y)) if isinstance(y, dict) else type(y).__name__)
                except Exception as e:
                    add(False, f"{hit.relative_to(path)} parses as YAML", str(e))
        for scr in re.findall(r"^\s*([A-Za-z0-9_.-]+\.(?:py|R|sh))\b", t, flags=re.M):
            add((path / "bin" / scr).is_file(), f"processes/{nf.name}: script {scr} is in bin/", warn=True)
    add(bool(procs), "at least one process in processes/")
    for inc_names, src in re.findall(r"include\s*\{([^}]*)\}\s*from\s*['\"]([^'\"]+)['\"]", txt):
        target = (path / src)
        f = target if target.suffix == ".nf" else target.with_suffix(".nf")
        add(f.is_file(), f"main.nf include resolves: {src}")
        base = [x.strip().split(" as ")[0].strip() for x in inc_names.split(";")]
        for b in base:
            if b:
                add(b in procs or not f.is_file(), f"included process {b} exists", str(f.name))
    for b in sorted((path / "bin").glob("*")) if (path / "bin").is_dir() else []:
        if b.is_file() and not b.name.startswith("."):
            first = b.read_text(errors="ignore").split("\n", 1)[0]
            add(first.startswith("#!"), f"bin/{b.name} has a shebang", first[:40])
            add(os.access(b, os.X_OK), f"bin/{b.name} is executable", "chmod +x (Nextflow calls bin/ scripts directly)")
    cfg = path / "input.config"
    if cfg.is_file():
        params = set(re.findall(r"params\.([A-Za-z0-9_]+)\s*=", cfg.read_text()))
        used = set(re.findall(r"params\.([A-Za-z0-9_]+)", txt))
        missing = sorted(used - params)
        add(not missing, "params used in main.nf are set in input.config", ", ".join(missing), warn=True)
    return out


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_new_module(args: argparse.Namespace) -> int:
    d = run_dir(args.label or f"new_module_{args.name}")
    dest = Path(args.dest) if args.dest else d
    try:
        root, written = create_new_module(args.name, dest, args.template, args.upstream_compat, args.force)
    except (ValueError, FileExistsError) as e:
        print(f"error: {e}")
        return 2
    print(f"Wrote: {root}")
    res = check_module(root)
    fails = [r for r in res if r[0] == "FAIL"]
    summ = {"module": args.name, "path": str(root), "template": args.template, "upstream_compat": bool(args.upstream_compat),
            "files": written, "dirs": MODULE_DIRS, "check_fail": len(fails), "check": [list(r) for r in res]}
    tree = [f"{args.name}/"] + [f"  {x}/" for x in MODULE_DIRS] + [f"  {x}" for x in written]
    finish(d, f"New Nextflow module: {args.name}", ["```", *tree, "```", "", md_table(["status", "check", "detail"], res)], summ)
    return 0


def cmd_check_module(args: argparse.Namespace) -> int:
    d = run_dir(args.label or f"check_{Path(args.path).name}")
    res = check_module(Path(args.path))
    n = {s: sum(1 for r in res if r[0] == s) for s in ("PASS", "WARN", "FAIL")}
    for r in res:
        print(f"  {r[0]:<5} {r[1]}" + (f"  ({r[2]})" if r[2] else ""))
    finish(d, f"Module check: {Path(args.path).name}", [f"{n['PASS']} pass, {n['WARN']} warn, {n['FAIL']} fail.", "",
                                                         md_table(["status", "check", "detail"], res)],
           {"path": str(args.path), "counts": n, "checks": [list(r) for r in res]})
    return 0 if n["FAIL"] == 0 else 1


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)
    existed = OUT_ROOT.exists()
    before = set(OUT_ROOT.glob("*")) if existed else set()
    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            print("\nupstream layout (--upstream-compat)")
            root, written = create_new_module("my_new_module", d / "compat", "minimal", upstream_compat=True)
            check(sorted(p.name for p in root.iterdir() if p.is_dir()) == sorted(MODULE_DIRS), "compat: bin, conda_envs, example_data, processes, test")
            check((root / "README.md").read_text() == "" and (root / "input.config").read_text() == ""
                  and (root / "main.nf").read_text() == "nextflow.enable.dsl=2\n"
                  and (root / "bin" / "my_new_module.py").read_text() == "#!/usr/bin/env python\n",
                  "compat: empty README/input.config, DSL2 main.nf, shebang-only bin script")
            y = (root / "conda_envs" / "my_new_module.yaml").read_text()
            check(y.startswith("name: my_new_module\n  channels:\n    - conda-forge") and "      - openpyxl==3.1.2" in y,
                  "compat: conda YAML with the upstream f-string indentation")
            try:
                parse_simple_yaml(y)
                bad = False
            except Exception:
                bad = True
            check(bad, "compat: the upstream YAML does not parse (mapping values not allowed / bad indentation)")
            pnf = (root / "processes" / "my_new_module.nf").read_text()
            check(pnf.startswith("process seqSpecParser {") and 'conda "${moduleDir}/conda_envs/my_new_module.yaml"' in pnf and "my_new_module.py" in pnf,
                  "compat: process hard-named seqSpecParser, conda env and script referenced")
            check(not os.access(root / "bin" / "my_new_module.py", os.X_OK), "compat: bin script not executable (as upstream)")
            res = {r[1]: r[0] for r in check_module(root)}
            check(res.get("conda_envs/my_new_module.yaml parses as YAML") == "FAIL" and res.get("bin/my_new_module.py is executable") == "FAIL",
                  "check-module flags the upstream YAML and the non-executable script")

            print("\ncorrected minimal template")
            root, written = create_new_module("guide_assign", d / "m", "minimal")
            env = parse_simple_yaml((root / "conda_envs" / "guide_assign.yaml").read_text())
            check(env == {"name": "guide_assign", "channels": CONDA_CHANNELS, "dependencies": CONDA_DEPS + [{"pip": CONDA_PIP}]},
                  "minimal: conda env parses to name / channels / dependencies (+ pip block)")
            check((root / "processes" / "guide_assign.nf").read_text().startswith("process guide_assign {")
                  and os.access(root / "bin" / "guide_assign.py", os.X_OK), "minimal: process named after the module; bin script executable")
            res = check_module(root)
            check(not [r for r in res if r[0] == "FAIL"], f"minimal: check-module passes ({sum(r[0] == 'PASS' for r in res)} checks)")
            try:
                create_new_module("guide_assign", d / "m", "minimal")
                clash = False
            except FileExistsError:
                clash = True
            check(clash, "new-module refuses to overwrite a non-empty module without --force")
            try:
                create_new_module("3bad-name", d / "m")
                badname = False
            except ValueError:
                badname = True
            check(badname, "new-module rejects names that are not Nextflow identifiers")

            print("\nfull template")
            root, written = create_new_module("de_test", d / "f", "full")
            mn = (root / "main.nf").read_text()
            check("include { de_test } from './processes/de_test'" in mn and "params.de_test_input" in mn
                  and "params.de_test_input" in (root / "input.config").read_text()
                  and (root / "test" / "test.nf").is_file() and (root / "example_data" / "de_test_example.txt").is_file(),
                  "full: include + workflow, input.config params, test/test.nf, example data")
            res = check_module(root)
            check(not [r for r in res if r[0] != "PASS"], f"full: check-module all PASS ({len(res)} checks)")
            import subprocess
            bs = root / "bin" / "de_test.py"
            (d / "in.txt").write_text("hello\n")
            rc = subprocess.call([sys.executable, str(bs), "--input", str(d / "in.txt"), "--output", str(d / "out.txt")])
            check(rc == 0 and (d / "out.txt").read_text() == "hello\n", "full: bin script skeleton runs")

            print("\nbroken module")
            bm = d / "broken"
            (bm / "processes").mkdir(parents=True)
            (bm / "main.nf").write_text("include { foo } from './processes/missing'\n")
            (bm / "processes" / "x.nf").write_text('process y {\n  conda "./conda_envs/none.yaml"\n  script:\n  """\n  x.py\n  """\n}\n')
            res = check_module(bm)
            fails = {r[1] for r in res if r[0] == "FAIL"}
            check({"bin/ present", "conda_envs/ present", "main.nf enables DSL2", "main.nf include resolves: ./processes/missing",
                   "processes/x.nf: conda env resolves"} <= fails, f"check-module: {len(fails)} failures on a broken module")

            print("\nCLI")
            rc = main(["new-module", "--name", "cli_mod", "--dest", str(d / "cli"), "--label", "st_devtools"])
            rc2 = main(["check-module", "--path", str(d / "cli" / "cli_mod"), "--label", "st_devtools_check"])
            check(rc == 0 and rc2 == 0 and (d / "cli" / "cli_mod" / "main.nf").is_file(), "CLI: new-module + check-module write report + summary")
    finally:
        after = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()
        for p in after - before:
            shutil.rmtree(p, ignore_errors=True)
        if not existed and OUT_ROOT.exists() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent crisprdevtools", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("new-module", help="scaffold a Nextflow module directory (create_new_module_nextflow)")
    p.add_argument("--name", required=True, help="module name (a Nextflow identifier)")
    p.add_argument("--dest", help="parent directory (default: the run directory under Docs/CRISPRDevTools)")
    p.add_argument("--template", choices=["minimal", "full"], default="minimal")
    p.add_argument("--upstream-compat", action="store_true", help="reproduce upstream exactly (indented YAML, seqSpecParser name, no chmod)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--label")
    p.set_defaults(func=cmd_new_module)
    p = sub.add_parser("check-module", help="validate a module directory layout")
    p.add_argument("--path", required=True)
    p.add_argument("--label")
    p.set_defaults(func=cmd_check_module)
    p = sub.add_parser("selftest", help="scaffold + check in a tempdir")
    p.add_argument("--no-plots", action="store_true", help="accepted for interface parity (no figures)")
    p.set_defaults(func=cmd_selftest)
    args = ap.parse_args(argv)
    if not logging.getLogger().handlers:
        print(f"Log: {setup_logging()}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
