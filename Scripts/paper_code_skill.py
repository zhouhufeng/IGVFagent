"""Reproduce a paper by running the authors' own code (paper-code skill).

Most papers that matter for IGVF ship their analysis: an R Markdown file, a
Jupyter notebook or a set of scripts in a GitHub repository named in the Code
Availability statement, usually beside the input tables and often beside the
authors' own rendered output (an .html knit, a notebook with stored outputs).
Rebuilding a paper from deposited score tables reproduces a number; running
that code reproduces the paper: every figure, every printed statistic, in the
authors' order, and it can be checked against what the authors rendered.

This follows Paper2Agent's route for code-backed papers (find the repository,
pin a commit, build the declared environment, execute the analysis unmodified,
compare with the authors' rendered outputs, split the result into sections),
as one IGVFagent command the agent and long-running jobs can call.

Subcommands
    find        paper harvest (bench harvest's harvest.json), a DOI/title, or
                --repo owner/name  ->  ranked code repositories
    fetch       clone at a pinned commit into Data/PaperCode/<owner>__<repo>/src
    inventory   entry points (.Rmd/.qmd/.ipynb/.R/.py), inputs they read and
                whether they exist, packages they load, output directories they
                write, the authors' rendered reference outputs, section outline
    env         R: micromamba env (r-base, rmarkdown, pandoc, poppler and
                conda-forge r-<pkg> for each package); Python: uv/venv
    run         execute the entry unmodified in a work copy (Rmd: rmarkdown::
                render with keep_md, error=TRUE so one failing chunk does not
                hide the rest; ipynb: nbconvert --execute --allow-errors;
                scripts: in order). --detach returns at once; job_wait on
                <run>/done.json
    compare     printed values and figure counts against the authors' rendered
                output (identical blocks, numerically close blocks, missing)
    report      report.md + report.html, one section per top-level heading of
                the analysis, with its figures, printed-value agreement and
                errors; summary.json for plan checks
    pipeline    fetch -> inventory -> env -> run -> compare -> report
    status      state of a run directory
    selftest    offline checks (and a real knit when R + rmarkdown exist)

Outputs (stdout lines harvested by the agent): Report:, JSON:, TSV:, Figure:.

Key numbers in summary.json (usable as job plan checks)
    printed.fraction_matched   reference printed blocks reproduced identically
                               or within 1e-6 relative tolerance / reference blocks
    figures.produced           figures the run produced (knitr + saved files)
    figures.reference          images embedded in the authors' rendered output
    chunks.errored             chunks whose output contains an R/Python error
"""
from __future__ import annotations

import argparse
import hashlib
import html as _html
import io
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
CODE_DIR = Path(os.environ.get("IGVF_PAPER_CODE_DIR") or ROOT / "Data" / "PaperCode")
REPORT_DIR = ROOT / "Docs" / "PaperCode"
CRAN = "https://cloud.r-project.org"

R_BASE_PKGS = {"base", "grid", "stats", "utils", "methods", "graphics", "grDevices", "tools", "parallel",
               "splines", "stats4", "tcltk", "compiler", "datasets"}
R_ALWAYS = ["r-base", "r-rmarkdown", "r-knitr", "pandoc", "poppler"]
PY_PIP_NAME = {"sklearn": "scikit-learn", "cv2": "opencv-python", "PIL": "pillow", "yaml": "pyyaml",
               "Bio": "biopython", "skimage": "scikit-image", "bs4": "beautifulsoup4", "umap": "umap-learn",
               "igraph": "python-igraph", "sc": "scanpy", "mpl_toolkits": "matplotlib"}
PY_ALWAYS = ["nbconvert", "nbclient", "ipykernel", "nbformat", "matplotlib"]
_PY_STDLIB = set(getattr(sys, "stdlib_module_names", ())) or {
    "os", "sys", "re", "json", "math", "time", "datetime", "collections", "itertools", "functools", "pathlib",
    "subprocess", "random", "string", "glob", "shutil", "csv", "pickle", "gzip", "io", "typing", "argparse",
    "logging", "copy", "statistics", "warnings", "urllib", "tempfile", "hashlib", "textwrap", "zipfile"}
ENTRY_EXT = (".Rmd", ".rmd", ".qmd", ".ipynb", ".R", ".r", ".py")


def say(kind: str, value: Any) -> None:
    print(f"{kind}: {value}", flush=True)


def _write_json(path: Path, obj: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)
    return path


def _read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return default


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


# ─── find: which repository holds the paper's code ──────────────────────────

_GH = re.compile(r"github\.com/([A-Za-z0-9_.\-]+)/([A-Za-z0-9_.\-]+)")
_CODE_WORDS = re.compile(r"\b(code|script|analys[ie]s|notebook|software|pipeline|R markdown|reproduc)", re.I)


def clean_repo(owner: str, name: str) -> str:
    name = re.sub(r"\.git$", "", name.rstrip(".,;:)]}'\""))
    return f"{owner}/{name}"


def rank_repos(harvest: Dict[str, Any]) -> "List[Dict[str, Any]]":
    """Candidate repositories from a bench harvest, best first.

    The one named in the Code Availability sentence ("Code used for the
    analyses ... is available at github.com/X/Y") wins over tools the paper
    merely used (a subassembly pipeline, a plotting library)."""
    seen: Dict[str, Dict[str, Any]] = {}
    for hit in (harvest.get("accessions") or {}).get("github_repo") or []:
        m = _GH.search("github.com/" + hit.get("value", "")) or _GH.search(hit.get("context", ""))
        if not m:
            continue
        repo = clean_repo(*m.groups())
        ctx = hit.get("context") or ""
        # score the sentence around THIS repo, not the whole window
        i = ctx.lower().find(repo.lower())
        near = ctx[max(0, i - 160): i + 20] if i >= 0 else ctx
        score = 1.0 + (2.0 if hit.get("in_data_availability") else 0.0)
        if re.search(r"(code|analys[ie]s)[^.]{0,80}(used|performed|available)|available[^.]{0,40}code", near, re.I):
            score += 3.0
        if re.search(r"subassembl|pipeline for|tool|package", near, re.I) and not re.search(r"analys", near, re.I):
            score -= 1.0
        score += 0.5 * min(int(hit.get("count") or 1), 4)
        cur = seen.get(repo.lower())
        if not cur or score > cur["score"]:
            seen[repo.lower()] = {"repo": repo, "score": round(score, 2), "context": near.strip()}
    for d in harvest.get("data_availability") or []:
        for m in _GH.finditer(d.get("text", "") if isinstance(d, dict) else str(d)):
            repo = clean_repo(*m.groups())
            seen.setdefault(repo.lower(), {"repo": repo, "score": 2.5, "context": "data availability"})
    return sorted(seen.values(), key=lambda r: -r["score"])


def _gh_api(path: str) -> Any:
    req = urllib.request.Request("https://api.github.com/" + path.lstrip("/"),
                                 headers={"Accept": "application/vnd.github+json", "User-Agent": "igvfagent"})
    tok = os.environ.get("GITHUB_TOKEN")
    if tok:
        req.add_header("Authorization", f"Bearer {tok}")
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def find(harvest_path: Optional[str] = None, repo: Optional[str] = None,
         query: Optional[str] = None) -> Dict[str, Any]:
    cands: List[Dict[str, Any]] = []
    if repo:
        m = _GH.search(repo) or re.fullmatch(r"([\w.\-]+)/([\w.\-]+)", repo.strip())
        if not m:
            return {"ok": False, "error": f"not a GitHub repository: {repo!r}"}
        cands = [{"repo": clean_repo(*m.groups()), "score": 10.0, "context": "given"}]
    elif harvest_path:
        hv = _read_json(Path(harvest_path))
        if not hv:
            return {"ok": False, "error": f"cannot read harvest {harvest_path}"}
        cands = rank_repos(hv)
    elif query:
        try:
            res = _gh_api("search/repositories?q=" + urllib.request.quote(query) + "&per_page=5")
            cands = [{"repo": it["full_name"], "score": 1.0, "context": (it.get("description") or "")[:160]}
                     for it in res.get("items", [])]
        except Exception as e:  # noqa: BLE001
            return {"ok": False, "error": f"GitHub search failed: {e}"}
    return {"ok": bool(cands), "candidates": cands,
            "error": None if cands else "no code repository named in the paper; try --query or --repo"}


# ─── fetch: a pinned clone ──────────────────────────────────────────────────

def repo_dir(repo: str) -> Path:
    return CODE_DIR / repo.replace("/", "__")


def fetch(repo: str, ref: Optional[str] = None) -> Dict[str, Any]:
    d = repo_dir(repo)
    src = d / "src"
    url = f"https://github.com/{repo}.git"
    if not shutil.which("git"):
        return _fetch_tarball(repo, ref, d)
    if not (src / ".git").is_dir():
        src.parent.mkdir(parents=True, exist_ok=True)
        p = subprocess.run(["git", "clone", "--quiet", url, str(src)], capture_output=True, text=True)
        if p.returncode:
            return {"ok": False, "error": f"git clone failed: {p.stderr.strip()[-400:]}"}
    if ref:
        subprocess.run(["git", "-C", str(src), "fetch", "--quiet", "origin", ref], capture_output=True)
        p = subprocess.run(["git", "-C", str(src), "checkout", "--quiet", ref], capture_output=True, text=True)
        if p.returncode:
            return {"ok": False, "error": f"cannot check out {ref}: {p.stderr.strip()[-300:]}"}
    commit = subprocess.run(["git", "-C", str(src), "rev-parse", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    meta = {"repo": repo, "url": f"https://github.com/{repo}", "commit": commit, "src": _rel(src),
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    _write_json(d / "source.json", meta)
    return {"ok": True, **meta}


def _fetch_tarball(repo: str, ref: Optional[str], d: Path) -> Dict[str, Any]:
    """No git on this host: resolve the commit through the GitHub API and
    unpack that exact commit's tarball, so the pin is still recorded."""
    try:
        commit = _gh_api(f"repos/{repo}/commits/{ref or 'HEAD'}")["sha"]
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"cannot resolve {repo}@{ref or 'HEAD'}: {e}"}
    src = d / "src"
    meta = _read_json(d / "source.json", {})
    if not (src.is_dir() and meta.get("commit") == commit):
        req = urllib.request.Request(f"https://codeload.github.com/{repo}/tar.gz/{commit}",
                                     headers={"User-Agent": "igvfagent"})
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read()
        if src.exists():
            shutil.rmtree(src)
        tmp = d / "src.tmp"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            for m in tf.getmembers():  # no absolute paths or parent escapes
                if m.name.startswith("/") or ".." in Path(m.name).parts or m.issym() or m.islnk():
                    continue
                tf.extract(m, tmp)
        (next(tmp.iterdir())).rename(src)
        shutil.rmtree(tmp, ignore_errors=True)
    meta = {"repo": repo, "url": f"https://github.com/{repo}", "commit": commit, "src": _rel(src),
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "via": "tarball"}
    _write_json(d / "source.json", meta)
    return {"ok": True, **meta}


# ─── inventory: what the analysis is, what it needs, what it should produce ─

_FENCE_R = re.compile(r"^\s*```+\s*\{\s*r\b([^}]*)\}\s*$")
_FENCE_ANY = re.compile(r"^\s*```+\s*\{?\s*([A-Za-z]+)?")
_FENCE_END = re.compile(r"^\s*```+\s*$")


def parse_rmd(text: str) -> Dict[str, Any]:
    """Chunks (label, options, code, line) and headings, fence-aware: a `#`
    inside a code chunk is a comment, not a heading."""
    chunks: List[Dict[str, Any]] = []
    heads: List[Dict[str, Any]] = []
    prose: Dict[int, List[str]] = {}
    lines = text.splitlines()
    i, n_unnamed, in_yaml = 0, 0, bool(lines and lines[0].strip() == "---")
    if in_yaml:
        i = 1
        while i < len(lines) and lines[i].strip() != "---":
            i += 1
        i += 1
    while i < len(lines):
        ln = lines[i]
        m = _FENCE_R.match(ln)
        if m:
            hdr = m.group(1).strip().lstrip(",").strip()
            label, opts = "", hdr
            first = hdr.split(",", 1)[0].strip()
            if first and "=" not in first:
                label, opts = first, hdr[len(first):].lstrip(", ")
            if not label:
                n_unnamed += 1
            start = i
            code: List[str] = []
            i += 1
            while i < len(lines) and not _FENCE_END.match(lines[i]):
                code.append(lines[i])
                i += 1
            chunks.append({"index": len(chunks) + 1, "label": label, "auto_label": label or f"unnamed-chunk-{n_unnamed}",
                           "options": opts, "code": "\n".join(code), "line": start + 1,
                           "section": len(heads) - 1 if heads else -1,
                           "include": not re.search(r"include\s*=\s*F(ALSE)?\b", opts),
                           "eval": not re.search(r"eval\s*=\s*F(ALSE)?\b", opts)})
        elif re.match(r"^(#{1,6})\s+\S", ln):
            lvl = len(ln) - len(ln.lstrip("#"))
            heads.append({"level": lvl, "title": ln.lstrip("#").strip(), "line": i + 1})
        elif ln.strip() and heads:
            prose.setdefault(len(heads) - 1, []).append(ln.strip())
        i += 1
    for k, h in enumerate(heads):
        h["prose"] = " ".join(prose.get(k, []))[:600]
    return {"chunks": chunks, "headings": heads}


def sections_of(parsed: Dict[str, Any]) -> "List[Dict[str, Any]]":
    """Top-level sections (the smallest heading level used), each with its chunks."""
    heads = parsed["headings"]
    if not heads:
        return [{"title": "Analysis", "line": 1, "chunks": [c["index"] for c in parsed["chunks"]], "prose": ""}]
    top = min(h["level"] for h in heads)
    secs: List[Dict[str, Any]] = []
    owner = {}
    for k, h in enumerate(heads):
        if h["level"] == top:
            secs.append({"title": h["title"], "line": h["line"], "chunks": [], "prose": h.get("prose", "")})
        owner[k] = max(len(secs) - 1, 0)
    if not secs:
        secs = [{"title": "Analysis", "line": 1, "chunks": [], "prose": ""}]
    for c in parsed["chunks"]:
        secs[owner.get(c["section"], 0)]["chunks"].append(c["index"])
    return secs


_R_LIB = re.compile(r"\b(?:library|require|requireNamespace)\s*\(\s*['\"]?([A-Za-z][A-Za-z0-9.]*)")
_R_NS = re.compile(r"\b([A-Za-z][A-Za-z0-9.]*):::?[A-Za-z_.]")
_R_READ = re.compile(r"\b(?:read\.[a-z]+|read_[a-z]+|fread|readRDS|load|read\.table|readLines|file)\s*\(\s*"
                     r"(?:file\s*=\s*)?['\"]([^'\"]+)['\"]")
_R_WRITE = re.compile(r"\b(?:ggsave|pdf|png|svg|jpeg|tiff|write\.[a-z]+|write_[a-z]+|fwrite|saveRDS|writeLines|"
                      r"sink|cat)\s*\(\s*(?:file(?:name)?\s*=\s*)?['\"]([^'\"]+\.[A-Za-z0-9]{1,5})['\"]")
_PY_IMPORT = re.compile(r"^\s*(?:from\s+([A-Za-z_][\w]*)|import\s+([A-Za-z_][\w]*))", re.M)
_PY_READ = re.compile(r"\b(?:read_csv|read_table|read_excel|read_parquet|read_h5ad|open|loadtxt|load|read)\s*\(\s*"
                      r"['\"]([^'\"]+\.[A-Za-z0-9]{1,6})['\"]")
_PY_WRITE = re.compile(r"\b(?:savefig|to_csv|to_parquet|write|write_h5ad)\s*\(\s*['\"]([^'\"]+\.[A-Za-z0-9]{1,6})['\"]")


def _entry_code(path: Path) -> Tuple[str, Dict[str, Any]]:
    text = path.read_text(errors="replace")
    if path.suffix.lower() in (".rmd", ".qmd"):
        parsed = parse_rmd(text)
        return "\n".join(c["code"] for c in parsed["chunks"]), parsed
    if path.suffix == ".ipynb":
        nb = json.loads(text)
        cells = [("".join(c.get("source") or [])) for c in nb.get("cells", []) if c.get("cell_type") == "code"]
        return "\n".join(cells), {"chunks": [{"index": i + 1, "label": "", "auto_label": f"cell-{i + 1}", "code": s,
                                               "line": 0, "section": -1, "include": True, "eval": True}
                                              for i, s in enumerate(cells)],
                                  "headings": _nb_headings(nb)}
    return text, {"chunks": [{"index": 1, "label": "", "auto_label": path.stem, "code": text, "line": 1,
                              "section": -1, "include": True, "eval": True}], "headings": []}


def _nb_headings(nb: dict) -> "List[Dict[str, Any]]":
    out, code_i = [], 0
    for c in nb.get("cells", []):
        if c.get("cell_type") == "code":
            code_i += 1
            continue
        for ln in "".join(c.get("source") or []).splitlines():
            if re.match(r"^#{1,6}\s+\S", ln):
                out.append({"level": len(ln) - len(ln.lstrip("#")), "title": ln.lstrip("#").strip(),
                            "line": code_i, "prose": ""})
    return out


def language_of(path: Path) -> str:
    return "python" if path.suffix in (".ipynb", ".py") else "r"


def author_package_versions(code: str, ref: Optional[Path]) -> Dict[str, str]:
    """Package versions the authors ran, from their comments (`#'2.2.1'` after
    packageVersion("x")) and from the printed output of their rendering
    (sessionInfo, "x version: 1.2.3"). Recorded, and diffed in the report:
    a difference is the first suspect when a value or figure disagrees."""
    out: Dict[str, str] = {}
    for m in re.finditer(r"packageVersion\(\s*['\"]([\w.]+)['\"]\s*\)\)?\s*#\s*[‘'\"]?([\d][\w.\-]*)", code):
        out[m.group(1)] = m.group(2)
    m = re.search(r"^\s*version\s*#\s*R version ([\d.]+)", code, re.M)
    if m:
        out["R"] = m.group(1)
    if ref and ref.suffix == ".html":
        t = _html.unescape(ref.read_text(errors="replace"))
        for pk, v in re.findall(r"\b([A-Za-z][\w.]*) version: ([\d][\w.\-]*)", t):
            out.setdefault(pk, v)
        m = re.search(r"R version (\d+\.\d+\.\d+)", t)
        if m:
            out.setdefault("R", m.group(1))
        for pk, v in re.findall(r"\b([A-Za-z][\w.]*)_(\d+\.\d+(?:[.\-]\d+)*)\b", t[t.find("other attached"):][:4000]
                                if "other attached" in t else ""):
            out.setdefault(pk, v)
    return out


def reference_output(src: Path, entry: Path) -> Optional[Path]:
    """The authors' own rendered output for an entry, if they committed it."""
    if entry.suffix == ".ipynb":
        nb = _read_json(entry) or {}
        has = any(c.get("outputs") for c in nb.get("cells", []) if c.get("cell_type") == "code")
        return entry if has else None
    for ext in (".html", ".nb.html", ".md", ".pdf"):
        cand = entry.with_suffix(ext)
        if cand.exists() and cand != entry:
            return cand
    return None


def inventory(src: Path, entry: Optional[str] = None) -> Dict[str, Any]:
    files = [p for p in src.rglob("*") if p.is_file() and ".git" not in p.parts]
    entries = [p for p in files if p.suffix in ENTRY_EXT]
    # Notebooks/R Markdown first, larger (the main analysis) first.
    def rank(p: Path) -> Tuple[int, int]:
        k = {".Rmd": 0, ".rmd": 0, ".qmd": 0, ".ipynb": 1}.get(p.suffix, 2)
        return (k, -p.stat().st_size)
    entries.sort(key=rank)
    if entry:
        chosen = [p for p in entries if p.name == entry or str(p.relative_to(src)) == entry]
        if not chosen:
            return {"ok": False, "error": f"entry {entry!r} not found; entries: "
                                          + ", ".join(str(p.relative_to(src)) for p in entries[:20])}
        main = chosen[0]
    elif entries:
        main = entries[0]
    else:
        return {"ok": False, "error": "no .Rmd/.qmd/.ipynb/.R/.py analysis in the repository"}
    code, parsed = _entry_code(main)
    lang = language_of(main)
    if lang == "r":
        pkgs = sorted({m for m in _R_LIB.findall(code)} | {m for m in _R_NS.findall(code)} - R_BASE_PKGS)
        pkgs = [p for p in pkgs if p not in R_BASE_PKGS]
        reads = sorted(set(_R_READ.findall(code)))
        writes = sorted(set(_R_WRITE.findall(code)))
    else:
        mods = {a or b for a, b in _PY_IMPORT.findall(code)}
        local = {p.stem for p in files if p.suffix == ".py"}
        pkgs = sorted(m for m in mods if m and m not in _PY_STDLIB and m not in local)
        reads = sorted(set(_PY_READ.findall(code)))
        writes = sorted(set(_PY_WRITE.findall(code)))
    base = main.parent
    if lang == "r":  # file("x") connections later passed to writeLines/cat/sink are outputs
        for var, path in re.findall(r"([A-Za-z_.][\w.]*)\s*<-\s*file\(\s*['\"]([^'\"]+)['\"]", code):
            if re.search(r"\b(writeLines|cat|sink|write)\s*\([^)]*\b" + re.escape(var) + r"\b", code):
                writes.append(path)
        writes = sorted(set(writes))
    out_dirs = sorted({str(Path(w).parent) for w in writes if str(Path(w).parent) not in (".", "")})
    inputs = [{"path": r, "exists": (base / r).exists() or (src / r).exists(),
               "is_url": bool(re.match(r"https?://", r)),
               "written_by_analysis": r in writes or str(Path(r).parent) in out_dirs} for r in reads]
    inputs = [i for i in inputs if not i["written_by_analysis"]]
    ref = reference_output(src, main)
    ref_stats = reference_stats(ref) if ref else {}
    author_versions = author_package_versions(code, ref)
    reqs = [str(p.relative_to(src)) for p in files
            if p.name in ("requirements.txt", "environment.yml", "environment.yaml", "renv.lock", "DESCRIPTION",
                          "setup.py", "pyproject.toml", "install.R")]
    return {"ok": True, "src": _rel(src), "entry": str(main.relative_to(src)), "language": lang,
            "entries": [str(p.relative_to(src)) for p in entries[:30]], "packages": pkgs, "inputs": inputs,
            "inputs_missing": [i["path"] for i in inputs if not i["exists"] and not i["is_url"]],
            "output_dirs": out_dirs, "outputs_declared": writes, "declared_env_files": reqs,
            "author_versions": author_versions,
            "reference": str(ref.relative_to(src)) if ref else None, "reference_stats": ref_stats,
            "n_chunks": len(parsed["chunks"]),
            "sections": [{"title": s["title"], "chunks": s["chunks"]} for s in sections_of(parsed)]}


# ─── compare: printed values and figures versus the authors' rendering ─────

_NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _norm_block(lines: "List[str]") -> str:
    return "\n".join(re.sub(r"\s+", " ", ln).strip() for ln in lines if ln.strip())


def html_printed_blocks(html_text: str) -> "List[str]":
    """knitr prints results as <pre><code>## ...</code></pre> in the HTML."""
    out = []
    for m in re.finditer(r"<pre(?:\s[^>]*)?><code(?:\s[^>]*)?>(.*?)</code></pre>", html_text, re.S):
        body = _html.unescape(re.sub(r"<[^>]+>", "", m.group(1)))
        lines = body.splitlines()
        if lines and all(ln.startswith("##") or not ln.strip() for ln in lines):
            out.append(_norm_block([ln[2:] for ln in lines]))
    return [b for b in out if b]


def md_printed_blocks(md_text: str) -> "List[str]":
    """Printed output in a knitr keep_md file: fenced blocks of `## ` lines."""
    out, cur, in_fence = [], [], False
    for ln in md_text.splitlines():
        if _FENCE_END.match(ln) or re.match(r"^\s*```+\s*\w*\s*$", ln):
            if in_fence:
                if cur and all(x.startswith("##") or not x.strip() for x in cur):
                    out.append(_norm_block([x[2:] for x in cur]))
                cur, in_fence = [], False
            else:
                in_fence = True
            continue
        if in_fence:
            cur.append(ln)
    return [b for b in out if b]


def nb_printed_blocks(nb: dict) -> "List[str]":
    out = []
    for c in nb.get("cells", []):
        for o in c.get("outputs") or []:
            txt = o.get("text") or (o.get("data") or {}).get("text/plain")
            if txt:
                out.append(_norm_block(("".join(txt) if isinstance(txt, list) else txt).splitlines()))
    return [b for b in out if b]


def _close(a: str, b: str, rtol: float = 1e-6) -> bool:
    if _NUM.sub("#", a) != _NUM.sub("#", b):
        return False
    na, nb_ = [float(x) for x in _NUM.findall(a)], [float(x) for x in _NUM.findall(b)]
    return len(na) == len(nb_) and all(abs(x - y) <= rtol * max(1.0, abs(x), abs(y)) for x, y in zip(na, nb_))


_MESSAGE = re.compile(r"^(Warning|Message|Scale for|Coordinate system|`stat_bin\(\)`|`geom_smooth\(\)`|"
                      r"Picking joint bandwidth|Using .* as id variables|Loading required package|"
                      r"will replace the existing scale|Removed \d+ rows|In addition:|Attaching package|"
                      r"The following objects? (is|are) masked|Registered S3|Saving \d)", re.I)


def is_message(block: str) -> bool:
    """Warnings and notes from libraries (ggplot2 "Removed 8 rows ..."), not
    computed values: their wording changes between package versions, so they
    are reported apart and do not count against the reproduction."""
    lines = [ln for ln in block.splitlines() if ln.strip()]
    return bool(lines) and all(_MESSAGE.search(ln.strip()) or ln.startswith(" ") for ln in lines) \
        and bool(_MESSAGE.search(lines[0]))


def _tokens(b: str) -> "List[str]":
    return sorted(b.replace("\n", " ").split())


def _same_values(a: str, b: str) -> bool:
    """Same printed values in a different layout (a matrix wrapped at another
    console width prints its columns in different line groups)."""
    ta, tb = _tokens(a), _tokens(b)
    if ta == tb:
        return True
    if len(ta) != len(tb):
        return False
    na = sorted(float(x) for x in ta if _NUM.fullmatch(x))
    nb_ = sorted(float(x) for x in tb if _NUM.fullmatch(x))
    wa = sorted(x for x in ta if not _NUM.fullmatch(x))
    wb = sorted(x for x in tb if not _NUM.fullmatch(x))
    return wa == wb and len(na) == len(nb_) and all(abs(x - y) <= 1e-6 * max(1.0, abs(x), abs(y))
                                                    for x, y in zip(na, nb_))


def compare_blocks(ref: "List[str]", got: "List[str]") -> Dict[str, Any]:
    pool = list(got)
    identical, close, layout, missing, msgs, msg_missing, near = 0, 0, 0, [], 0, [], []
    for b in ref:
        if is_message(b):
            msgs += 1
            if b in pool:
                pool.remove(b)
            else:
                msg_missing.append(b[:200])
            continue
        if b in pool:
            pool.remove(b)
            identical += 1
            continue
        k = next((i for i, g in enumerate(pool) if _close(b, g)), None)
        if k is not None:
            pool.pop(k)
            close += 1
            continue
        k = next((i for i, g in enumerate(pool) if _same_values(b, g)), None)
        if k is not None:
            pool.pop(k)
            layout += 1
        else:
            missing.append(b[:300])
            # Same printed structure with other numbers (e.g. a Monte Carlo
            # p-value): keep both side by side with the largest difference.
            skel = _NUM.sub("#", b)
            k = next((i for i, g in enumerate(pool) if _NUM.sub("#", g) == skel), None)
            if k is not None:
                g = pool.pop(k)
                na, ng = [float(x) for x in _NUM.findall(b)], [float(x) for x in _NUM.findall(g)]
                diffs = [(x, y) for x, y in zip(na, ng) if x != y]
                rel = max((abs(x - y) / max(abs(x), abs(y), 1e-300) for x, y in diffs), default=0.0)
                near.append({"reference": b[:300], "produced": g[:300], "values_differing": len(diffs),
                             "values_total": len(na), "max_relative_difference": round(rel, 4),
                             "pairs": [[x, y] for x, y in diffs[:8]]})
    n = len(ref) - msgs
    return {"reference_blocks": n, "produced_blocks": len(got), "identical": identical, "numerically_close": close,
            "same_values_other_layout": layout, "missing": len(missing),
            "fraction_matched": round((identical + close + layout) / n, 4) if n else None,
            "missing_examples": missing[:15], "extra_blocks": len(pool), "same_structure_other_numbers": near[:20],
            "library_messages": {"reference": msgs, "reproduced": msgs - len(msg_missing),
                                 "differing_examples": msg_missing[:6]}}


def reference_stats(ref: Path) -> Dict[str, Any]:
    if ref.suffix == ".ipynb":
        nb = _read_json(ref) or {}
        imgs = sum(1 for c in nb.get("cells", []) for o in c.get("outputs") or []
                   if any(k.startswith("image/") for k in (o.get("data") or {})))
        return {"printed_blocks": len(nb_printed_blocks(nb)), "images": imgs}
    if ref.suffix == ".html":
        t = ref.read_text(errors="replace")
        return {"printed_blocks": len(html_printed_blocks(t)), "images": len(re.findall(r"<img\b", t))}
    return {}


# ─── env: micromamba (R) or uv/venv (Python) ────────────────────────────────

def _mamba_platform() -> str:
    s, m = platform.system(), platform.machine().lower()
    if s == "Darwin":
        return "osx-arm64" if m in ("arm64", "aarch64") else "osx-64"
    return "linux-aarch64" if m in ("arm64", "aarch64") else "linux-64"


def micromamba() -> Path:
    for cand in (os.environ.get("IGVF_MICROMAMBA"), shutil.which("micromamba")):
        if cand and Path(cand).exists():
            return Path(cand)
    dest = CODE_DIR / "bin" / "micromamba"
    if dest.exists():
        return dest
    url = f"https://micro.mamba.pm/api/micromamba/{_mamba_platform()}/latest"
    with urllib.request.urlopen(url, timeout=120) as r:
        data = r.read()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as tf:
        member = tf.getmember("bin/micromamba")
        f = tf.extractfile(member)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(f.read())  # type: ignore[union-attr]
    dest.chmod(0o755)
    return dest


def r_conda_name(pkg: str) -> str:
    return "r-" + pkg.lower()


def env_dir(spec: "List[str]") -> Path:
    h = hashlib.sha256("\n".join(sorted(spec)).encode()).hexdigest()[:12]
    return CODE_DIR / "envs" / h


def build_env(inv: Dict[str, Any], pins: "Sequence[str]" = (), log: Optional[Path] = None) -> Dict[str, Any]:
    """Create (or reuse) the environment the entry needs. Returns its prefix
    and which packages had to come from CRAN/pip because conda had none."""
    logf = open(log, "a") if log else None

    def run(argv: "List[str]", **kw) -> subprocess.CompletedProcess:
        if logf:
            logf.write("$ " + " ".join(argv) + "\n")
            logf.flush()
        if "env" not in kw and "/envs/" in argv[0]:  # the env's compilers/pandoc come first on PATH
            kw["env"] = {**os.environ, "PATH": f"{Path(argv[0]).parent}{os.pathsep}{os.environ.get('PATH', '')}"}
        p = subprocess.run(argv, capture_output=True, text=True, **kw)
        if logf:
            logf.write((p.stdout or "")[-4000:] + (p.stderr or "")[-4000:] + "\n")
            logf.flush()
        return p

    try:
        if inv["language"] == "r":
            cran_pins = [p[5:] for p in pins if p.startswith("cran:")]
            pins = [p for p in pins if not p.startswith("cran:")]
            pinned = {re.split(r"[=<>]", p)[0]: p for p in pins}
            pkgs = list(inv["packages"])
            if "data.table" in pkgs and "reshape2" not in pkgs:
                pkgs.append("reshape2")  # data.table's melt/dcast on a data.frame redirect to reshape2
            wanted = [pinned.get(r_conda_name(p), r_conda_name(p)) for p in pkgs]
            if cran_pins:  # building archived CRAN versions needs compilers
                pins = pins + ["c-compiler", "cxx-compiler", "fortran-compiler", "make"]
            spec = sorted(set(R_ALWAYS + wanted + [p for p in pins if p.split("=")[0] not in
                                                    {w.split("=")[0] for w in wanted}])) + \
                [f"cran:{c}" for c in sorted(cran_pins)]
            prefix = env_dir(spec)
            ok_marker = prefix / ".igvf_env_ok.json"
            if ok_marker.exists():
                return {"ok": True, "prefix": str(prefix), "reused": True, **_read_json(ok_marker, {})}
            mm = micromamba()
            root = CODE_DIR / "mamba-root"
            base = [str(mm), "-r", str(root), "create", "-y", "-p", str(prefix), "-c", "conda-forge", "-c", "bioconda"]
            conda_spec = [x for x in spec if not x.startswith("cran:")]
            p = run(base + conda_spec)
            from_cran: List[str] = []
            if p.returncode:
                # Solve failed on some package: build the base, then add the rest one by one.
                p = run(base + sorted(set(R_ALWAYS + [p for p in pins])))
                if p.returncode:
                    return {"ok": False, "error": "micromamba could not create the R base env; see log",
                            "log": str(log) if log else None}
                inst = [str(mm), "-r", str(root), "install", "-y", "-p", str(prefix), "-c", "conda-forge",
                        "-c", "bioconda"]
                for pkg in pkgs:
                    if run(inst + [pinned.get(r_conda_name(pkg), r_conda_name(pkg))]).returncode and \
                            run(inst + ["bioconductor-" + pkg.lower()]).returncode:
                        from_cran.append(pkg)
            rscript = prefix / "bin" / "Rscript"
            if from_cran:
                expr = (f"options(repos=c(CRAN='{CRAN}')); pk <- c({','.join(repr(x) for x in from_cran)}); "
                        "for (x in pk) if (!requireNamespace(x, quietly=TRUE)) { "
                        "tryCatch(install.packages(x), error=function(e) NULL); "
                        "if (!requireNamespace(x, quietly=TRUE)) { if (!requireNamespace('BiocManager', quietly=TRUE)) "
                        "install.packages('BiocManager'); BiocManager::install(x, ask=FALSE, update=FALSE) } }")
                run([str(rscript), "-e", expr])
            for c in cran_pins:  # e.g. reshape2@1.4.4 -> the CRAN archive tarball
                nm, _, ver = c.partition("@")
                url = f"{CRAN}/src/contrib/Archive/{nm}/{nm}_{ver}.tar.gz"
                run([str(rscript), "-e", f"options(repos=c(CRAN='{CRAN}')); install.packages('{url}', repos=NULL, "
                     "type='source')"])
                have = run([str(rscript), "-e", f"cat(as.character(packageVersion('{nm}')))"]).stdout.strip()
                if have != ver:
                    from_cran.append(f"{nm}@{ver} (FAILED, have {have or 'none'})")
            chk = run([str(rscript), "-e", "pk <- c(" + ",".join(repr(x) for x in pkgs) +
                       "); cat(paste(pk[!vapply(pk, requireNamespace, logical(1), quietly=TRUE)], collapse=','))"])
            missing = [x for x in (chk.stdout or "").strip().split(",") if x]
            missing += [c.split(" ")[0] for c in from_cran if "(FAILED" in c]
            vers = run([str(rscript), "-e", "cat(R.version.string, '\\n'); for (p in c(" +
                        ",".join(repr(x) for x in pkgs + ["rmarkdown", "knitr"]) +
                        ")) if (requireNamespace(p, quietly=TRUE)) cat(p, as.character(packageVersion(p)), '\\n')"])
            info = {"spec": spec, "from_cran": from_cran, "missing": missing, "versions": (vers.stdout or "").strip()}
            if not missing:
                _write_json(ok_marker, info)
            return {"ok": not missing, "prefix": str(prefix), "reused": False, **info,
                    "error": f"packages not installable: {', '.join(missing)}" if missing else None}
        # python
        pip_names = sorted({PY_PIP_NAME.get(m, m) for m in inv["packages"]} | set(PY_ALWAYS) | set(pins))
        prefix = env_dir(["py"] + pip_names)
        ok_marker = prefix / ".igvf_env_ok.json"
        if ok_marker.exists():
            return {"ok": True, "prefix": str(prefix), "reused": True, **_read_json(ok_marker, {})}
        uv = shutil.which("uv")
        if uv:
            p = run([uv, "venv", "--quiet", "--python", os.environ.get("IGVF_PAPER_CODE_PY", "3.11"), str(prefix)])
            if p.returncode:
                p = run([uv, "venv", "--quiet", str(prefix)])
        else:
            p = run([sys.executable, "-m", "venv", str(prefix)])
        if p.returncode:
            return {"ok": False, "error": "could not create a Python venv; see log"}
        py = prefix / "bin" / "python"
        reqs = [ROOT / inv["src"] / r for r in inv.get("declared_env_files", []) if r.endswith("requirements.txt")]
        failed = []
        inst = ([uv, "pip", "install", "--quiet", "--python", str(py)] if uv else [str(py), "-m", "pip", "install", "-q"])
        for r in reqs:
            run(inst + ["-r", str(r)])
        if run(inst + pip_names).returncode:
            for n in pip_names:
                if run(inst + [n]).returncode:
                    failed.append(n)
        vers = run([str(py), "-m", "pip", "freeze"]) if not uv else run([uv, "pip", "freeze", "--python", str(py)])
        info = {"spec": pip_names, "missing": failed, "versions": (vers.stdout or "")[-6000:]}
        if not failed:
            _write_json(ok_marker, info)
        return {"ok": not failed, "prefix": str(prefix), "reused": False, **info,
                "error": f"pip could not install: {', '.join(failed)}" if failed else None}
    finally:
        if logf:
            logf.close()


# ─── run: execute the entry unmodified ──────────────────────────────────────

_KNIT_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
entry <- args[1]; tolerant <- args[2] == "1"
options(repos = c(CRAN = "%s"))
fmt <- rmarkdown::html_document(keep_md = TRUE, self_contained = FALSE)
fmt$knitr$opts_chunk$error <- tolerant
fmt$knitr$opts_chunk$dev <- c("png", "pdf")
fmt$knitr$opts_chunk$dpi <- 110
if (Sys.getenv("IGVF_R_SAMPLE_ROUNDING") == "1") suppressWarnings(RNGkind(sample.kind = "Rounding"))
shims <- strsplit(Sys.getenv("IGVF_R_SHIMS"), ",")[[1]]
if ("melt_reshape" %%in%% shims) {
  # Before R 3.6, reshape2::melt's UseMethod found reshape's exported
  # melt.data.frame on the search path; R 3.6 stopped searching it. Once
  # data.table is attached (it masks melt), put reshape::melt in front, which is
  # what melt() resolved to when the analysis was written. The Rmd is unchanged.
  fmt$knitr$knit_hooks$igvf_compat <- function(before, options, envir) {
    if (before && "package:data.table" %%in%% search() && !("igvf_compat_melt" %%in%% search()))
      attach(list(melt = reshape::melt), name = "igvf_compat_melt", warn.conflicts = FALSE)
    NULL
  }
  fmt$knitr$opts_chunk$igvf_compat <- TRUE
}
t0 <- Sys.time()
res <- tryCatch({ rmarkdown::render(entry, output_format = fmt, envir = new.env(), quiet = FALSE); "ok" },
                error = function(e) paste("render error:", conditionMessage(e)))
writeLines(c(res, format(difftime(Sys.time(), t0, units = "secs"))), "igvf_render_status.txt")
writeLines(capture.output(sessionInfo()), "igvf_sessionInfo.txt")
""" % CRAN


R_SHIMS = {
    "melt_reshape": "melt() resolves to reshape::melt as it did before R 3.6 (the analysis attaches reshape, "
                    "then data.table, and melts data.frames)",
}


def r_shims(inv: Dict[str, Any], code_of_entry: Optional[Path] = None) -> "List[str]":
    """R-compatibility shims this analysis needs, by detection. Each restores a
    behaviour of the R the authors ran without touching their code."""
    if inv.get("language") != "r" or not code_of_entry or code_of_entry.suffix.lower() not in (".rmd", ".qmd"):
        return []
    code = _entry_code(code_of_entry)[0]
    out = []
    lib_order = [m for m in _R_LIB.findall(code)]
    if "reshape" in lib_order and "data.table" in lib_order and re.search(r"(?<![\w.:])melt\s*\(", code) \
            and lib_order.index("reshape") < lib_order.index("data.table"):
        out.append("melt_reshape")
    return out


def _snapshot(d: Path) -> Dict[str, float]:
    return {str(p.relative_to(d)): p.stat().st_mtime for p in d.rglob("*") if p.is_file()}


def run_entry(src: Path, inv: Dict[str, Any], env: Dict[str, Any], run_dir: Path,
              tolerant: bool = True, timeout_min: float = 180, work_name: str = "work",
              inputs: "Optional[Dict[str, str]]" = None) -> Dict[str, Any]:
    work = run_dir / work_name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(src, work, ignore=shutil.ignore_patterns(".git"))
    entry = work / inv["entry"]
    for od in inv.get("output_dirs") or []:
        (entry.parent / od).mkdir(parents=True, exist_ok=True)
    substituted = {}
    for name, path in (inputs or {}).items():
        # New data, same code: the user's file takes the place of the input the
        # analysis reads, under the name the analysis expects.
        dst = entry.parent / name
        shutil.copyfile(path, dst)
        substituted[name] = {"from": str(path), "sha256": hashlib.sha256(dst.read_bytes()).hexdigest()}
    before = _snapshot(work)
    prefix = Path(env["prefix"])
    envvars = dict(os.environ)
    envvars["PATH"] = f"{prefix / 'bin'}{os.pathsep}{envvars.get('PATH', '')}"
    envvars.pop("R_LIBS_USER", None)
    envvars.pop("R_HOME", None)
    envvars["MPLBACKEND"] = "Agg"
    rv = (inv.get("author_versions") or {}).get("R")
    if rv and tuple(int(x) for x in re.findall(r"\d+", rv)[:2]) < (3, 6):
        # R 3.6 changed sample()'s algorithm; the authors' set.seed() draws need the old one.
        envvars["IGVF_R_SAMPLE_ROUNDING"] = "1"
    shims = r_shims(inv, code_of_entry=entry)
    envvars["IGVF_R_SHIMS"] = ",".join(shims)
    t0 = time.time()
    log = run_dir / ("run.log" if work_name == "work" else f"{work_name}.log")
    if inv["language"] == "r" and entry.suffix.lower() in (".rmd", ".qmd"):
        drv = run_dir / "knit_driver.R"
        drv.write_text(_KNIT_DRIVER)
        argv = [str(prefix / "bin" / "Rscript"), "--vanilla", str(drv), entry.name, "1" if tolerant else "0"]
    elif inv["language"] == "r":
        argv = [str(prefix / "bin" / "Rscript"), "--vanilla", entry.name]
    elif entry.suffix == ".ipynb":
        argv = [str(prefix / "bin" / "python"), "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--output", entry.stem + ".executed.ipynb", f"--ExecutePreprocessor.timeout={int(timeout_min * 60)}"]
        argv += (["--allow-errors"] if tolerant else []) + [entry.name]
    else:
        argv = [str(prefix / "bin" / "python"), entry.name]
    with open(log, "w") as lf:
        lf.write("$ " + " ".join(argv) + f"\n(cwd {entry.parent})\n")
        lf.flush()
        try:
            p = subprocess.run(argv, cwd=str(entry.parent), stdout=lf, stderr=subprocess.STDOUT, env=envvars,
                               timeout=timeout_min * 60)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = 124
    after = _snapshot(work)
    new = sorted(k for k, t in after.items() if k not in before or t > before[k] + 1e-6)
    status_txt = (entry.parent / "igvf_render_status.txt")
    render = status_txt.read_text().splitlines()[0] if status_txt.exists() else ("ok" if rc == 0 else f"exit {rc}")
    return {"exit_code": rc, "render": render, "seconds": round(time.time() - t0, 1), "new_files": new,
            "rng_sample_rounding": envvars.get("IGVF_R_SAMPLE_ROUNDING") == "1", "inputs_substituted": substituted,
            "r_shims": shims,
            "log": _rel(log), "work": _rel(work), "argv": argv}


def _stable_bytes(p: Path, root: Optional[Path] = None) -> bytes:
    b = p.read_bytes()
    if root is not None:  # analyses that write getwd() into their outputs
        b = b.replace(str(root.resolve()).encode(), b"<WORKDIR>").replace(str(root).encode(), b"<WORKDIR>")
    if p.suffix.lower() == ".pdf":  # timestamps and ids differ between identical renders
        b = re.sub(rb"/(CreationDate|ModDate)\s*\([^)]*\)", b"", b)
        b = re.sub(rb"/ID\s*\[[^\]]*\]", b"", b)
    return b


def replay_compare(run_dir: Path, inv: Dict[str, Any]) -> Dict[str, Any]:
    """Paper2Agent's replay check: a second execution from the same inputs must
    print the same values and save the same files (PDFs compared without their
    timestamps). Differences mean the analysis is not deterministic here."""
    a, b = run_dir / "work", run_dir / "replay"
    out: Dict[str, Any] = {}
    ma, mb = _produced_md(a, inv), _produced_md(b, inv)
    if ma and mb:
        pa, pb = md_printed_blocks(ma.read_text(errors="replace")), md_printed_blocks(mb.read_text(errors="replace"))
        out["printed_identical"] = sum(1 for x, y in zip(pa, pb) if x == y)
        out["printed_blocks"] = max(len(pa), len(pb))
    ra, rb = _read_json(run_dir / "run.json", {}), _read_json(run_dir / "replay.json", {})
    files = sorted(set(ra.get("new_files", [])) & set(rb.get("new_files", [])))
    files = [f for f in files if Path(f).suffix.lower() in (".pdf", ".png", ".tsv", ".csv", ".txt", ".pml", ".rds")
             and "_files/" not in f and not Path(f).name.startswith("igvf_")]
    same = [f for f in files if _stable_bytes(a / f) == _stable_bytes(b / f)]
    by_path = [f for f in files if f not in same and _stable_bytes(a / f, a) == _stable_bytes(b / f, b)]
    same += by_path
    out.update({"files_compared": len(files), "files_identical": len(same),
                "identical_after_path_normalisation": by_path,
                "files_differing": [f for f in files if f not in same][:30]})
    out["deterministic"] = bool(out.get("printed_identical") == out.get("printed_blocks")
                                and len(same) == len(files))
    return out


# ─── harvest a run: figures, printed values, errors, sections ───────────────

_ERR_LINE = re.compile(r"^##\s*(Error\b|Error in\b)", re.M)
_CASCADE = re.compile(r"object '[^']+' not found|could not find function|undefined columns selected|"
                      r"argument of length 0|cannot open file 'Output/", re.I)
# Known upstream breakages -> the environment repair that fixes them. Only the
# environment is ever changed; the authors' code stays as written.
REPAIR_HINTS = [
    (re.compile(r"melt generic in data\.table has been passed a data\.frame", re.I),
     ["cran:data.table@1.16.4", "cran:reshape2@1.4.4"],
     "data.table >= 1.17 raises an error instead of redirecting melt()/dcast() on a data.frame to reshape2 "
     "(checked: 1.17.8 and 1.18.6 both fail); 1.16.4 still redirects. reshape2 1.4.5 changed melt's "
     "attribute-length check, so 1.4.4 matches the authors."),
    (re.compile(r"dims \[product \d+\] do not match the length of object", re.I), ["cran:reshape2@1.4.4"],
     "reshape2 1.4.5 melt() rejects this table; 1.4.4 is what the analysis was written against."),
    (re.compile(r"there is no package called ['‘]([\w.]+)", re.I), [],
     "a package the analysis loads is missing from the environment: add it with --pin r-<name>."),
]


def _md_chunk_map(md_text: str, chunks: "List[Dict[str, Any]]") -> "List[Tuple[int, int]]":
    """(md line, chunk index) for every echoed code block, by matching its first
    code line to the Rmd chunk it came from (in document order)."""
    firsts = []
    for c in chunks:
        ln = next((x.strip() for x in c["code"].splitlines() if x.strip()), "")
        firsts.append((c["index"], ln))
    out, k, lines = [], 0, md_text.splitlines()
    for i, ln in enumerate(lines):
        if re.match(r"^\s*```+\s*\{?\s*r\b", ln) and i + 1 < len(lines):
            first = lines[i + 1].strip()
            for j in range(k, len(firsts)):
                if firsts[j][1] and firsts[j][1] == first:
                    out.append((i, firsts[j][0]))
                    k = j
                    break
    return out


def md_errors(md_text: str, chunks: "List[Dict[str, Any]]") -> "List[Dict[str, Any]]":
    cmap = _md_chunk_map(md_text, chunks)
    lines = md_text.splitlines()
    by_idx = {c["index"]: c for c in chunks}
    errs = []
    for i, ln in enumerate(lines):
        if not re.match(r"^##\s*Error\b", ln):
            continue
        msg = [ln[2:].strip()]
        for x in lines[i + 1:i + 6]:
            if not x.startswith("##"):
                break
            msg.append(x[2:].strip().lstrip("! ").strip())
        text = " ".join(m for m in msg if m)
        ch = next((c for (l, c) in reversed(cmap) if l < i), None)
        c = by_idx.get(ch) if ch else None
        e = {"message": text[:500], "chunk": ch, "chunk_label": (c or {}).get("auto_label"),
             "rmd_line": (c or {}).get("line"), "kind": "cascade" if _CASCADE.search(text) else "root"}
        for rx, pins, why in REPAIR_HINTS:
            if rx.search(text):
                e["repair"] = {"pins": pins, "why": why}
                break
        errs.append(e)
    return errs


def _produced_md(work: Path, inv: Dict[str, Any]) -> Optional[Path]:
    e = work / inv["entry"]
    for cand in (e.with_suffix(".md"), e.with_suffix(".knit.md")):
        if cand.exists():
            return cand
    return None


def _pdf_to_png(pdf: Path, out: Path, prefix: Optional[Path]) -> Optional[Path]:
    exe = None
    for cand in ([prefix / "bin" / "pdftoppm"] if prefix else []) + [shutil.which("pdftoppm")]:
        if cand and Path(cand).exists():
            exe = str(cand)
            break
    if not exe:
        return None
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    p = subprocess.run([exe, "-png", "-r", "90", "-singlefile", str(pdf), str(stem)], capture_output=True)
    return out if p.returncode == 0 and out.exists() else None


def _norm_label(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def harvest_run(run_dir: Path, inv: Dict[str, Any], env: Dict[str, Any]) -> Dict[str, Any]:
    run = _read_json(run_dir / "run.json", {})
    work = run_dir / "work"
    entry = work / inv["entry"]
    code, parsed = _entry_code(entry)
    chunks = parsed["chunks"]
    new = [work / f for f in run.get("new_files", [])]
    prefix = Path(env["prefix"]) if env.get("prefix") else None
    figs_dir = run_dir / "figures"
    figures: List[Dict[str, Any]] = []
    # knitr's own figures (png; the pdf twin is kept for print quality)
    for f in new:
        if f.suffix.lower() == ".png" and "figure-" in str(f.parent.name) + str(f):
            figures.append({"path": f, "source": "knitr", "png": f})
    # files the analysis saved itself (ggsave/pdf()/savefig)
    for f in new:
        if f.suffix.lower() in (".pdf", ".png", ".svg") and not any(x["path"] == f for x in figures) \
                and "figure-" not in str(f) and f.name != "Rplots.pdf":
            png = f if f.suffix.lower() == ".png" else _pdf_to_png(f, figs_dir / (f.stem + ".png"), prefix)
            figures.append({"path": f, "source": "saved", "png": png})
    # assign each figure to a chunk
    for fg in figures:
        name = fg["path"].name
        fg["chunk"] = None
        if fg["source"] == "knitr":
            base = _norm_label(re.sub(r"-\d+\.png$", "", name))
            for c in chunks:
                if _norm_label(c["auto_label"]) == base:
                    fg["chunk"] = c["index"]
                    break
        else:
            for c in chunks:
                if name in c["code"] or fg["path"].stem in c["code"]:
                    fg["chunk"] = c["index"]
                    break
    printed: List[str] = []
    errors: List[Dict[str, Any]] = []
    if entry.suffix.lower() in (".rmd", ".qmd"):
        md = _produced_md(work, inv)
        if md:
            mdt = md.read_text(errors="replace")
            printed = md_printed_blocks(mdt)
            errors = md_errors(mdt, chunks)
    elif entry.suffix == ".ipynb":
        ex = entry.with_name(entry.stem + ".executed.ipynb")
        nb = _read_json(ex, {}) if ex.exists() else {}
        printed = nb_printed_blocks(nb)
        for i, c in enumerate([c for c in nb.get("cells", []) if c.get("cell_type") == "code"]):
            for o in c.get("outputs") or []:
                if o.get("output_type") == "error":
                    msg = f"{o.get('ename')}: {o.get('evalue')}"[:400]
                    errors.append({"chunk": i + 1, "message": msg,
                                   "kind": "cascade" if o.get("ename") == "NameError" else "root"})
                if o.get("output_type") == "display_data" and "image/png" in (o.get("data") or {}):
                    import base64
                    figs_dir.mkdir(parents=True, exist_ok=True)
                    pth = figs_dir / f"cell{i + 1}_{len(figures) + 1}.png"
                    pth.write_bytes(base64.b64decode(o["data"]["image/png"]))
                    figures.append({"path": pth, "source": "notebook", "png": pth, "chunk": i + 1})
    ref_rel = inv.get("reference")
    cmp: Dict[str, Any] = {}
    if ref_rel:
        ref = ROOT / inv["src"] / ref_rel
        if ref.suffix == ".html":
            cmp = compare_blocks(html_printed_blocks(ref.read_text(errors="replace")), printed)
        elif ref.suffix == ".ipynb":
            cmp = compare_blocks(nb_printed_blocks(_read_json(ref, {})), printed)
    return {"figures": figures, "printed": printed, "errors": errors, "compare": cmp, "parsed": parsed}


# ─── report ─────────────────────────────────────────────────────────────────

def write_report(run_dir: Path, inv: Dict[str, Any], env: Dict[str, Any], src_meta: Dict[str, Any],
                 paper: Optional[str] = None) -> Dict[str, Any]:
    run = _read_json(run_dir / "run.json", {})
    hv = harvest_run(run_dir, inv, env)
    secs = sections_of(hv["parsed"])
    by_chunk: Dict[Any, List[Dict[str, Any]]] = {}
    for fg in hv["figures"]:
        by_chunk.setdefault(fg.get("chunk"), []).append(fg)
    ref_stats = inv.get("reference_stats") or {}
    cmp = hv["compare"]
    summary = {
        "paper": paper, "repo": src_meta.get("repo"), "commit": src_meta.get("commit"), "entry": inv["entry"],
        "language": inv["language"], "render": run.get("render"), "exit_code": run.get("exit_code"),
        "seconds": run.get("seconds"),
        "chunks": {"total": len(hv["parsed"]["chunks"]),
                   "errored": len({e.get("chunk") for e in hv["errors"]}),
                   "root_errors": sum(1 for e in hv["errors"] if e.get("kind") == "root"),
                   "cascade_errors": sum(1 for e in hv["errors"] if e.get("kind") == "cascade")},
        "repair": sorted({p for e in hv["errors"] for p in (e.get("repair") or {}).get("pins", [])}),
        "author_versions_diff": _version_diff(inv.get("author_versions") or {}, env.get("versions") or ""),
        "replay": _read_json(run_dir / "replay_compare.json"),
        "figures": {"produced": len(hv["figures"]), "reference": ref_stats.get("images"),
                    "saved_files": sum(1 for f in hv["figures"] if f["source"] == "saved")},
        "printed": {k: v for k, v in cmp.items() if k != "missing_examples"} if cmp else {},
        "inputs_missing": inv.get("inputs_missing"), "env": {k: env.get(k) for k in ("prefix", "from_cran", "missing")},
        "sections": [], "errors": hv["errors"][:40]}
    L = [f"# Reproduction from the authors' code: {src_meta.get('repo')}", ""]
    if paper:
        L += [f"**Paper:** {paper}", ""]
    L += [f"- Repository: [{src_meta.get('repo')}]({src_meta.get('url')}) at commit `{(src_meta.get('commit') or '')[:12]}`",
          f"- Entry: `{inv['entry']}` ({inv['language'].upper()}), run unmodified"
          + (" with error-tolerant chunks" if run.get("argv") and "1" in run.get("argv", [])[-1:] else ""),
          f"- Render: **{run.get('render')}** in {run.get('seconds')} s; "
          f"{summary['chunks']['errored']} of {summary['chunks']['total']} chunks raised an error "
          f"({summary['chunks']['root_errors']} root cause(s), {summary['chunks']['cascade_errors']} knock-on)",
          f"- Figures: **{summary['figures']['produced']}** produced"
          + (f" (the authors' rendering has {ref_stats.get('images')})" if ref_stats.get("images") is not None else "")
          + f", {summary['figures']['saved_files']} saved by the analysis itself"]
    if cmp:
        L.append(f"- Printed values vs the authors' rendering: **{cmp['identical']} identical + "
                 f"{cmp['numerically_close']} numerically close + {cmp.get('same_values_other_layout', 0)} same "
                 f"values in another layout, of {cmp['reference_blocks']}** "
                 f"({(cmp['fraction_matched'] or 0) * 100:.1f}%); {cmp['missing']} differ")
        lm = cmp.get("library_messages") or {}
        if lm.get("reference"):
            L.append(f"- Library warnings/notes (not values; wording varies by package version): "
                     f"{lm['reproduced']}/{lm['reference']} worded identically")
    else:
        L.append("- The repository has no rendered reference output, so printed values could not be compared.")
    if run.get("inputs_substituted"):
        L.append("- **New data, the authors' code:** " + "; ".join(
            f"`{k}` ← `{v['from']}`" for k, v in run["inputs_substituted"].items())
                 + ". Comparison with the authors' rendering then measures how the results changed, not reproduction.")
    for sh in run.get("r_shims") or []:
        L.append(f"- Compatibility shim `{sh}`: {R_SHIMS.get(sh, '')}. The authors' code is unchanged.")
    if run.get("rng_sample_rounding"):
        L.append(f"- The authors ran R {(inv.get('author_versions') or {}).get('R')}; R 3.6 changed sample(), so "
                 "this run used `RNGkind(sample.kind = \"Rounding\")` to reproduce their seeded random draws")
    rp = summary.get("replay")
    if rp:
        L.append(f"- Replay (second execution): printed values {rp.get('printed_identical')}/{rp.get('printed_blocks')} "
                 f"identical, saved files {rp['files_identical']}/{rp['files_compared']} identical"
                 + (f" ({len(rp.get('identical_after_path_normalisation') or [])} only after replacing the "
                    "working-directory path the analysis writes into them)"
                    if rp.get("identical_after_path_normalisation") else "")
                 + (" — **deterministic**" if rp.get("deterministic") else " — **not deterministic**: "
                    + ", ".join(rp.get("files_differing", [])[:5])))
    if inv.get("inputs_missing"):
        L.append(f"- Inputs referenced but not in the repository: {', '.join(inv['inputs_missing'][:10])}")
    L.append("")
    for si, s in enumerate(secs, 1):
        figs = [f for c in s["chunks"] for f in by_chunk.get(c, [])]
        errs = [e for e in hv["errors"] if e.get("chunk") in s["chunks"]]
        summary["sections"].append({"title": s["title"], "chunks": len(s["chunks"]), "figures": len(figs)})
        L += [f"## {si}. {s['title']}", ""]
        if s.get("prose"):
            L += [f"> {s['prose'][:500]}", ""]
        L.append(f"{len(s['chunks'])} chunk(s), {len(figs)} figure(s).")
        L.append("")
        for f in figs:
            png = f.get("png")
            if png:
                L.append(f"![{Path(f['path']).name}]({_url(png, run_dir)})")
            else:
                L.append(f"- [{Path(f['path']).name}]({_url(f['path'], run_dir)})")
        for e in errs:
            L.append(f"- ⚠️ error: `{e['message']}`")
        L.append("")
    unassigned = by_chunk.get(None, [])
    if unassigned:
        L += ["## Other figures", ""]
        for f in unassigned:
            png = f.get("png")
            L.append(f"![{Path(f['path']).name}]({_url(png, run_dir)})" if png else
                     f"- [{Path(f['path']).name}]({_url(f['path'], run_dir)})")
        L.append("")
    roots = [e for e in hv["errors"] if e.get("kind") != "cascade"]
    if roots:
        L += ["## Root-cause errors", ""]
        for e in roots[:30]:
            where = f"chunk {e.get('chunk')} `{e.get('chunk_label')}` (Rmd line {e.get('rmd_line')})" \
                if e.get("chunk") else "chunk ?"
            L.append(f"- {where}: `{e['message'][:300]}`")
            if e.get("repair"):
                L.append(f"  - environment repair: `{' '.join('--pin ' + p for p in e['repair']['pins'])}` — "
                         f"{e['repair']['why']}")
        casc = len(hv["errors"]) - len(roots)
        if casc:
            L.append(f"- plus {casc} knock-on error(s) (objects the failing chunks would have created).")
        L.append("")
    if summary["author_versions_diff"]:
        L += ["## Package versions: the authors' vs this run", "", "| package | authors | this run |", "|---|---|---|"]
        L += [f"| {k} | {a} | {b} |" for k, (a, b) in sorted(summary["author_versions_diff"].items())] + [""]
    near = cmp.get("same_structure_other_numbers") or []
    if near:
        L += ["## Same output, other numbers", "",
              "These blocks print the same thing as the authors' but some numbers differ. Judge each: a Monte Carlo "
              "p-value or a seeded draw can differ within sampling error when a library consumes random numbers "
              "differently; a deterministic statistic that differs is a real discrepancy.", ""]
        for n in near:
            L.append(f"- {n['values_differing']} of {n['values_total']} numbers differ, largest relative difference "
                     f"{n['max_relative_difference']:.3g}: " + ", ".join(f"{a:g} (authors) vs {b:g}"
                                                                         for a, b in n["pairs"]))
        L.append("")
    other = [b for b in cmp.get("missing_examples") or [] if b not in {n["reference"] for n in near}]
    if other:
        L += ["## Printed values missing from this run (first 15)", ""]
        for b in other:
            L += ["```", b, "```"]
        L.append("")
    L += ["## Environment", "", "```", (env.get("versions") or "")[:3000], "```", "",
          f"Run log: `{run.get('log')}`. Work copy: `{run.get('work')}`."]
    md = run_dir / "report.md"
    md.write_text("\n".join(L))
    _write_json(run_dir / "summary.json", summary)
    _write_json(run_dir / "compare.json", cmp)
    html = _md_to_html("\n".join(L), f"Reproduction: {src_meta.get('repo')}")
    (run_dir / "report.html").write_text(html)
    return {"summary": summary, "report": md, "html": run_dir / "report.html",
            "figures": [f.get("png") or f["path"] for f in hv["figures"]]}


def _url(p: Any, base: Path) -> str:
    return urllib.request.quote(os.path.relpath(p, base))


def _version_diff(author: Dict[str, str], versions_text: str) -> Dict[str, Tuple[str, str]]:
    ours: Dict[str, str] = {}
    for ln in versions_text.splitlines():
        m = re.match(r"^\s*R version ([\d.]+)", ln)
        if m:
            ours["R"] = m.group(1)
            continue
        parts = ln.split()
        if len(parts) == 2:
            ours[parts[0]] = parts[1]
    return {k: (v, ours.get(k, "-")) for k, v in author.items() if ours.get(k) != v and k != "grid"}


def _md_to_html(md: str, title: str) -> str:
    out, in_code = [], False
    for ln in md.splitlines():
        if ln.startswith("```"):
            out.append("</pre>" if in_code else "<pre>")
            in_code = not in_code
            continue
        if in_code:
            out.append(_html.escape(ln))
            continue
        e = _html.escape(ln)
        e = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", r'<figure><img src="\2" alt="\1" loading="lazy"><figcaption>\1</figcaption></figure>', e)
        e = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', e)
        e = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", e)
        e = re.sub(r"`([^`]+)`", r"<code>\1</code>", e)
        m = re.match(r"^(#{1,4})\s+(.*)", e)
        if m:
            e = f"<h{len(m.group(1))}>{m.group(2)}</h{len(m.group(1))}>"
        elif e.startswith("|"):
            cells = [c.strip() for c in e.strip("|").split("|")]
            e = "" if set("".join(cells)) <= set("-: ") else "<div class='row'>" + "".join(
                f"<span>{c}</span>" for c in cells) + "</div>"
        elif e.startswith("- "):
            e = f"<li>{e[2:]}</li>"
        elif e.startswith("&gt; "):
            e = f"<blockquote>{e[5:]}</blockquote>"
        elif e.strip() and not e.startswith("<figure"):
            e = f"<p>{e}</p>"
        out.append(e)
    css = ("body{font-family:system-ui,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;line-height:1.5}"
           "figure{display:inline-block;margin:8px;max-width:48%}figure img{max-width:100%;border:1px solid #ddd}"
           "figcaption{font-size:12px;color:#555}pre{background:#f5f5f5;padding:8px;overflow:auto;font-size:12px}"
           "blockquote{color:#555;border-left:3px solid #ccc;margin:0;padding-left:10px}"
           ".row{display:grid;grid-template-columns:repeat(3,minmax(0,200px));gap:8px;font-size:13px}")
    return f"<!doctype html><html><head><meta charset='utf-8'><title>{_html.escape(title)}</title><style>{css}</style></head><body>" \
           + "\n".join(out) + "</body></html>"


# ─── pipeline and CLI ───────────────────────────────────────────────────────

def _new_run_dir(repo: str) -> Path:
    d = REPORT_DIR / (time.strftime("%Y%m%d_%H%M%S") + "_" + re.sub(r"[^A-Za-z0-9]+", "_", repo).strip("_"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def pipeline(repo: str, *, ref: Optional[str] = None, entry: Optional[str] = None, pins: "Sequence[str]" = (),
             strict: bool = False, paper: Optional[str] = None, run_dir: Optional[Path] = None,
             timeout_min: float = 180, replay: bool = False,
             inputs: "Optional[Dict[str, str]]" = None) -> Dict[str, Any]:
    run_dir = run_dir or _new_run_dir(repo)
    status = run_dir / "status.json"

    def stage(name: str, **kw) -> None:
        _write_json(status, {"state": "running", "stage": name, "at": time.time(), **kw})
        print(f"[paper-code] {name}", flush=True)

    def fail(msg: str, **kw) -> Dict[str, Any]:
        res = {"ok": False, "error": msg, "run_dir": _rel(run_dir), **kw}
        _write_json(status, {"state": "failed", "error": msg, "at": time.time()})
        _write_json(run_dir / "done.json", res)
        return res

    stage("fetch")
    src_meta = fetch(repo, ref)
    if not src_meta.get("ok"):
        return fail(src_meta["error"])
    src = ROOT / src_meta["src"]
    stage("inventory")
    inv = inventory(src, entry)
    _write_json(run_dir / "inventory.json", inv)
    if not inv.get("ok"):
        return fail(inv["error"])
    known = {i["path"] for i in inv.get("inputs") or []}
    bad = [k for k in (inputs or {}) if k not in known]
    if bad:
        return fail(f"--input names must be inputs the analysis reads ({', '.join(sorted(known)[:12])}); got {bad}")
    missing_new = [v for v in (inputs or {}).values() if not Path(v).is_file()]
    if missing_new:
        return fail(f"input files not found: {missing_new}")
    stage("env", packages=inv["packages"])
    env = build_env(inv, pins, log=run_dir / "env.log")
    _write_json(run_dir / "env.json", env)
    if not env.get("ok") and not env.get("prefix"):
        return fail(env.get("error") or "environment failed")
    stage("run", entry=inv["entry"])
    run = run_entry(src, inv, env, run_dir, tolerant=not strict, timeout_min=timeout_min, inputs=inputs)
    _write_json(run_dir / "run.json", run)
    if replay:
        stage("replay")
        _write_json(run_dir / "replay.json", run_entry(src, inv, env, run_dir, tolerant=not strict,
                                                       timeout_min=timeout_min, work_name="replay",
                                                       inputs=inputs))
        _write_json(run_dir / "replay_compare.json", replay_compare(run_dir, inv))
    stage("report")
    rep = write_report(run_dir, inv, env, src_meta, paper)
    res = {"ok": True, "run_dir": _rel(run_dir), "report": _rel(rep["report"]), "html": _rel(rep["html"]),
           "summary": _rel(run_dir / "summary.json"), **{k: rep["summary"][k] for k in ("render", "figures", "printed",
                                                                                          "chunks")},
           "env_warning": env.get("error")}
    _write_json(status, {"state": "done", "at": time.time()})
    _write_json(run_dir / "done.json", res)
    return res


def _announce(res: Dict[str, Any]) -> None:
    rd = ROOT / res["run_dir"]
    if res.get("report"):
        say("Report", res["report"])
        say("Report", res["html"])
        say("JSON", res["summary"])
    for f in sorted((rd / "figures").glob("*.png"))[:80] if (rd / "figures").is_dir() else []:
        say("Figure", _rel(f))
    for f in sorted((rd / "work").rglob("*-1.png"))[:60] if (rd / "work").is_dir() else []:
        say("Figure", _rel(f))


def _detach(argv: "List[str]", run_dir: Path) -> None:
    log = open(run_dir / "pipeline.log", "a")
    cmd = [sys.executable, "-m", "igvfagent.paper_code_skill"] if __package__ else [sys.executable, str(Path(__file__))]
    subprocess.Popen(cmd + argv + ["--run-dir", str(run_dir)], stdout=log, stderr=subprocess.STDOUT,
                     stdin=subprocess.DEVNULL, start_new_session=True, cwd=str(ROOT))


def cmd_selftest() -> int:
    import tempfile
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")

    hv = {"accessions": {"github_repo": [
        {"value": "FowlerLab/VAMPseq.", "count": 2, "in_data_availability": True,
         "context": "Code used for the analyses performed in this work is included as Supplementary Data 5, and also "
                    "available at http://github.com/FowlerLab/VAMPseq. VAMP-seq scores are available at x. Code used "
                    "for subassembly by PacBio is available at http://github.com/shendurelab/AssemblyByPacBio."},
        {"value": "shendurelab/AssemblyByPacBio.", "count": 2, "in_data_availability": True,
         "context": "VAMP-seq scores are available at x. Code used for subassembly by PacBio is available at "
                    "http://github.com/shendurelab/AssemblyByPacBio. General reagents"}]}}
    r = rank_repos(hv)
    check("analysis repository ranks first, trailing '.' stripped", r and r[0]["repo"] == "FowlerLab/VAMPseq")
    rmd = ("---\ntitle: t\n---\n# Setup\nSome prose.\n```{r load data}\n# a comment, not a heading\n"
           "x <- read.table(file = \"in.tsv\", sep='\\t')\nlibrary(ggplot2)\nlibrary(grid)\n```\n"
           "# Results\n```{r, echo=FALSE}\nggsave(\"Output/p1.pdf\")\nprint(1)\n```\n## sub\n```{r plot two}\nplot(1)\n```\n")
    p = parse_rmd(rmd)
    check("chunks parsed, # inside code is not a heading", len(p["chunks"]) == 3 and len(p["headings"]) == 3)
    check("unnamed chunk gets knitr's auto label", p["chunks"][1]["auto_label"] == "unnamed-chunk-1")
    secs = sections_of(p)
    check("top-level sections own their chunks", [s["title"] for s in secs] == ["Setup", "Results"]
          and secs[1]["chunks"] == [2, 3])
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "repo"
        src.mkdir()
        (src / "a.Rmd").write_text(rmd)
        (src / "in.tsv").write_text("a\n1\n")
        (src / "a.html").write_text("<pre><code>## [1] 1</code></pre><pre><code>## [1] 0.3333333</code></pre>"
                                    "<p><img src='x'></p><pre><code class='r'>x &lt;- 1</code></pre>")
        global ROOT
        old = ROOT
        ROOT = Path(td)
        try:
            inv = inventory(src)
        finally:
            ROOT = old
        check("inventory: entry, packages minus base, inputs, output dir", inv["entry"] == "a.Rmd"
              and inv["packages"] == ["ggplot2"] and inv["inputs_missing"] == [] and inv["output_dirs"] == ["Output"])
        check("reference rendering found with its printed blocks and images",
              inv["reference"] == "a.html" and inv["reference_stats"] == {"printed_blocks": 2, "images": 1})
    ref = html_printed_blocks("<pre><code>## [1] 1</code></pre><pre><code>## [1] 0.3333333</code></pre>"
                              "<pre><code>## [1] 7</code></pre>")
    got = md_printed_blocks("```r\nprint(1)\n```\n\n```\n## [1] 1\n```\n\n```\n## [1] 0.33333330001\n```\n")
    c = compare_blocks(ref, got)
    check("compare: identical, numerically close and missing counted",
          (c["identical"], c["numerically_close"], c["missing"]) == (1, 1, 1) and abs(c["fraction_matched"] - 0.6667) < 1e-3)
    check("numbers that differ beyond tolerance do not match", not _close("[1] 0.5", "[1] 0.6"))
    wrapped_a = "a b c\nx 1 2 3\nd\nx 4"
    wrapped_b = "a b\nx 1 2\nc d\nx 3 4"
    c2 = compare_blocks([wrapped_a, "Warning: Removed 8 rows containing missing values (geom_point)."],
                        [wrapped_b, "Warning: Removed 8 rows containing missing values or values outside the scale range"])
    check("a matrix wrapped at another width counts as the same values; library warnings are set apart",
          c2["same_values_other_layout"] == 1 and c2["reference_blocks"] == 1 and c2["fraction_matched"] == 1.0
          and c2["library_messages"] == {"reference": 1, "reproduced": 0, "differing_examples": c2["library_messages"]["differing_examples"]})
    # end-to-end on a Python script with a stub env (system interpreter)
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "repo"
        src.mkdir()
        (src / "analysis.py").write_text("import json\nprint('mean', sum([1,2,3])/3)\nopen('out.txt','w').write('x')\n")
        rd = Path(td) / "run"
        rd.mkdir()
        fake_prefix = Path(td) / "env"
        (fake_prefix / "bin").mkdir(parents=True)
        (fake_prefix / "bin" / "python").symlink_to(sys.executable)
        old = ROOT
        ROOT = Path(td)
        try:
            inv = inventory(src)
            run = run_entry(src, inv, {"prefix": str(fake_prefix)}, rd, timeout_min=2)
            _write_json(rd / "run.json", run)
            rep = write_report(rd, inv, {"prefix": str(fake_prefix), "versions": "stub"},
                               {"repo": "x/y", "url": "https://github.com/x/y", "commit": "abc"})
        finally:
            ROOT = old
        check("python script runs unmodified in a work copy and its outputs are seen",
              run["exit_code"] == 0 and "out.txt" in run["new_files"])
        check("report and summary written", rep["report"].exists() and rep["summary"]["chunks"]["total"] == 1)
    if shutil.which("Rscript") and subprocess.run(["Rscript", "-e", "stopifnot(requireNamespace('rmarkdown'))"],
                                                  capture_output=True).returncode == 0:
        print("  (R + rmarkdown present: a real knit is exercised by `paper-code pipeline`)")
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="igvfagent paper-code", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("find", help="Paper -> ranked code repositories")
    s.add_argument("--harvest", help="harvest.json from `igvfagent bench harvest/pipeline`")
    s.add_argument("--repo", help="owner/name or a GitHub URL")
    s.add_argument("--query", help="GitHub repository search text (title, first author)")
    s = sub.add_parser("fetch", help="Clone at a pinned commit")
    s.add_argument("repo")
    s.add_argument("--ref")
    s = sub.add_parser("inventory", help="Entry points, inputs, packages, reference outputs, sections")
    s.add_argument("repo")
    s.add_argument("--entry")
    s = sub.add_parser("env", help="Build the environment the entry needs")
    s.add_argument("repo")
    s.add_argument("--entry")
    s.add_argument("--pin", action="append", default=[], help="conda/pip pin, e.g. r-reshape2=1.4.4")
    for name in ("pipeline", "run"):
        s = sub.add_parser(name, help="fetch -> inventory -> env -> run -> compare -> report"
                           if name == "pipeline" else "Same as pipeline (the environment is reused)")
        s.add_argument("repo", nargs="?", help="owner/name or GitHub URL (or give --harvest)")
        s.add_argument("--harvest", help="pick the repository from a bench harvest.json")
        s.add_argument("--ref", help="commit/tag to pin (default: current HEAD, recorded)")
        s.add_argument("--entry", help="analysis file inside the repository (default: the main Rmd/notebook)")
        s.add_argument("--pin", action="append", default=[], help="conda/pip pin, e.g. r-reshape2=1.4.4")
        s.add_argument("--strict", action="store_true", help="stop at the first failing chunk")
        s.add_argument("--paper", help="paper title/DOI for the report")
        s.add_argument("--timeout-min", type=float, default=180)
        s.add_argument("--replay", action="store_true", help="execute twice and check the outputs are identical")
        s.add_argument("--input", action="append", default=[], metavar="NAME=PATH",
                       help="run the same code on new data: replace input NAME (as the analysis reads it) with PATH")
        s.add_argument("--detach", action="store_true", help="run in the background; job_wait on <run>/done.json")
        s.add_argument("--run-dir", help=argparse.SUPPRESS)
    s = sub.add_parser("report", help="Re-write the report of a finished run directory")
    s.add_argument("run_dir")
    s = sub.add_parser("status", help="State of a run directory")
    s.add_argument("run_dir")
    sub.add_parser("selftest", help="Offline checks")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "selftest":
        return cmd_selftest()
    if args.cmd == "find":
        res = find(args.harvest, args.repo, args.query)
        for c in res.get("candidates", []):
            print(f"  {c['score']:>5}  {c['repo']}   {c['context'][:120]}")
        if not res["ok"]:
            print(res["error"])
        return 0 if res["ok"] else 1
    if args.cmd == "fetch":
        res = fetch(find(repo=args.repo)["candidates"][0]["repo"], args.ref)
        print(json.dumps(res, indent=2))
        return 0 if res["ok"] else 1
    if args.cmd in ("inventory", "env"):
        repo = find(repo=args.repo)["candidates"][0]["repo"]
        meta = fetch(repo)
        if not meta.get("ok"):
            print(meta["error"])
            return 1
        inv = inventory(ROOT / meta["src"], args.entry)
        if args.cmd == "inventory" or not inv.get("ok"):
            print(json.dumps(inv, indent=2))
            return 0 if inv.get("ok") else 1
        res = build_env(inv, args.pin, log=repo_dir(repo) / "env.log")
        print(json.dumps({k: v for k, v in res.items() if k != "versions"}, indent=2))
        print(res.get("versions", ""))
        return 0 if res.get("ok") else 1
    if args.cmd in ("pipeline", "run"):
        repo = args.repo
        if not repo and args.harvest:
            f = find(harvest_path=args.harvest)
            if not f["ok"]:
                print(f["error"])
                return 1
            repo = f["candidates"][0]["repo"]
            print(f"Code repository (from the paper): {repo}")
        if not repo:
            print("give a repository or --harvest")
            return 2
        repo = find(repo=repo)["candidates"][0]["repo"]
        if args.detach and not args.run_dir:
            rd = _new_run_dir(repo)
            passthru = [a for a in (argv if argv is not None else sys.argv[1:]) if a != "--detach"]
            _detach(passthru, rd)
            print(f"Started in the background: {_rel(rd)}")
            print(f"Wait with job_wait on {_rel(rd / 'done.json')}, or `igvfagent paper-code status {_rel(rd)}`.")
            say("JSON", _rel(rd / "status.json"))
            return 0
        res = pipeline(repo, ref=args.ref, entry=args.entry, pins=args.pin, strict=args.strict, paper=args.paper,
                       run_dir=Path(args.run_dir) if args.run_dir else None, timeout_min=args.timeout_min,
                       replay=args.replay,
                       inputs=dict(x.split("=", 1) for x in args.input) or None)
        print(json.dumps({k: v for k, v in res.items()}, indent=2, default=str))
        if res.get("ok"):
            _announce(res)
        return 0 if res.get("ok") else 1
    if args.cmd == "status":
        rd = Path(args.run_dir)
        rd = rd if rd.is_absolute() else ROOT / rd
        print(json.dumps(_read_json(rd / "done.json") or _read_json(rd / "status.json") or {"state": "unknown"},
                         indent=2))
        return 0
    if args.cmd == "report":
        rd = Path(args.run_dir)
        rd = rd if rd.is_absolute() else ROOT / rd
        inv, env = _read_json(rd / "inventory.json"), _read_json(rd / "env.json", {})
        meta = _read_json(ROOT / inv["src"] / ".." / "source.json", {})
        rep = write_report(rd, inv, env, meta)
        _announce({"run_dir": _rel(rd), "report": _rel(rep["report"]), "html": _rel(rep["html"]),
                   "summary": _rel(rd / "summary.json")})
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
