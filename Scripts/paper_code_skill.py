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
import urllib.parse
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
PY_PIP_NAME = {"pkg_resources": "setuptools", "pyro": "pyro-ppl", "torch": "torch", "Bio": "biopython", "anndata": "anndata",
               "sklearn": "scikit-learn", "cv2": "opencv-python", "PIL": "pillow", "yaml": "pyyaml",
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
    old = _read_json(d / "source.json", {}) or {}
    if old.get("commit") == commit and old.get("src") == _rel(src):
        # Unchanged source: leave the file (and its hash) alone. A job that
        # verified it as evidence must not see it "change" on every fetch.
        return {"ok": True, **old}
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
    if meta.get("commit") == commit and (d / "src").is_dir():
        return {"ok": True, **meta}
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


def notebook_kernel(path: Path) -> Dict[str, str]:
    nb = _read_json(path, {}) or {}
    md = nb.get("metadata") or {}
    ks = md.get("kernelspec") or {}
    li = md.get("language_info") or {}
    return {"name": ks.get("name") or "", "language": (ks.get("language") or li.get("name") or "python").lower(),
            "version": str(li.get("version") or "")}


def language_of(path: Path) -> str:
    if path.suffix == ".ipynb":
        return "r" if notebook_kernel(path)["language"] in ("r", "ir") else "python"
    return "python" if path.suffix == ".py" else "r"


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
            "entries": [str(p.relative_to(src)) for p in entries[:300]], "packages": pkgs, "inputs": inputs,
            "inputs_missing": [i["path"] for i in inputs if not i["exists"] and not i["is_url"]],
            "output_dirs": out_dirs, "outputs_declared": writes, "declared_env_files": reqs,
            "author_versions": author_versions,
            "reference": str(ref.relative_to(src)) if ref else None, "reference_stats": ref_stats,
            "n_chunks": len(parsed["chunks"]),
            "kernel": notebook_kernel(main) if main.suffix == ".ipynb" else None,
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
            if inv.get("r_notebooks") or str(inv.get("entry", "")).endswith(".ipynb"):
                wanted += ["r-irkernel", "jupyter", "nbconvert"]  # R notebooks run on this env's IRkernel
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
        # python: the authors' declared environment first (Paper2Agent: restore
        # a supplied environment before adding anything)
        yml = declared_python_env(ROOT / inv["src"], inv)
        if yml is not None and not os.environ.get("IGVF_PAPER_CODE_IGNORE_DECLARED"):
            return build_env_declared(yml, inv, pins, run)
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


GPU_ONLY = re.compile(r"^(cudatoolkit|cudnn|cuda-[\w-]+|cuda|nccl|pytorch-mutex|pytorch-cuda|[\w-]+-gpu|nvidia-[\w-]+|"
                      r"libcublas[\w-]*|magma[\w-]*|cupy[\w-]*|jaxlib-cuda[\w-]*)$", re.I)


def _dep_name(spec: str) -> "Tuple[str, str]":
    m = re.match(r"^\s*([A-Za-z0-9_.\-]+)\s*(?:[=<>!~]=?\s*([^=\s,;]+))?", spec or "")
    return (m.group(1), m.group(2) or "") if m else (spec, "")


def declared_python_env(src: Path, inv: Dict[str, Any]) -> Optional[Path]:
    """The conda environment file the authors ship for their notebooks, if any:
    one that pins python and brings a Jupyter kernel or runner."""
    try:
        import yaml  # type: ignore
    except Exception:  # noqa: BLE001
        return None
    best, score = None, 0.0
    for rel in inv.get("declared_env_files") or []:
        if not rel.endswith((".yml", ".yaml")):
            continue
        try:
            doc = yaml.safe_load((src / rel).read_text()) or {}
        except Exception:  # noqa: BLE001
            continue
        deps = [d for d in doc.get("dependencies") or [] if isinstance(d, str)]
        names = {_dep_name(d)[0].lower() for d in deps}
        pips = [x for d in doc.get("dependencies") or [] if isinstance(d, dict) for x in d.get("pip") or []]
        names |= {_dep_name(x)[0].lower() for x in pips}
        if "python" not in names:
            continue
        sc = 2.0 + (2.0 if names & {"ipykernel", "jupyter", "nbconvert", "papermill", "jupyterlab", "notebook"} else 0) \
            + len(names) / 100.0 - (1.0 if "/" in rel.strip("/") and rel.count("/") > 3 else 0)
        if sc > score:
            best, score = src / rel, sc
    return best


DECLARED_ENV_BUILDER = "6"  # bump when build_env_declared changes what it installs
BUILD_TOOLS = {"cython", "numpy", "setuptools", "wheel", "setuptools-scm", "setuptools_scm", "pybind11",
               "scikit-build", "cmake", "ninja", "versioneer", "pip"}
PIP_SKIP = {"sklearn"}  # the deprecated "sklearn" dummy package refuses to install; scikit-learn is the real one
# import / pip name -> conda package names to try (conda-forge, bioconda)
CONDA_NAMES = {"moods": ["moods"], "moods-python": ["moods"], "dash_bio": ["dash-bio"], "dash-bio": ["dash-bio"],
               "pybedtools": ["pybedtools"], "pysam": ["pysam"], "pybigwig": ["pybigwig"], "cyvcf2": ["cyvcf2"],
               "htseq": ["htseq"], "pyranges": ["pyranges"], "scanpy": ["scanpy"], "bs4": ["beautifulsoup4"]}


def _conda_install(run, mm: Path, root: Path, prefix: Path, name: str) -> bool:
    """Install one package from conda-forge/bioconda into an existing env
    (for what pip cannot build: C extensions, bioinformatics tools)."""
    for cand in CONDA_NAMES.get(name.lower(), [name.lower(), name.lower().replace("_", "-")]):
        if run([str(mm), "-r", str(root), "install", "-y", "-p", str(prefix), "-c", "conda-forge", "-c", "bioconda",
                cand]).returncode == 0:
            return True
    return False


def _release_archive(name: str, version: str) -> Optional[str]:
    """`name @ <archive URL>` for a version that exists as a tag in the
    project's GitHub repository but was never uploaded to PyPI (e.g. the
    crispr-bean 0.2.9 a paper used). The repository comes from PyPI's
    project URLs."""
    try:
        meta = _http_json(f"https://pypi.org/pypi/{name}/json")
    except Exception:  # noqa: BLE001
        return None
    if version in (meta.get("releases") or {}) and meta["releases"][version]:
        return None  # on PyPI after all: the normal install handles it
    info = meta.get("info") or {}
    urls = list((info.get("project_urls") or {}).values()) + [info.get("home_page") or "", info.get("download_url") or ""]
    repos = [m.group(1) for u in urls for m in [re.search(r"github\.com/([\w.-]+/[\w.-]+?)(?:\.git)?/?$", u or "")] if m]
    for repo in dict.fromkeys(repos):
        for tag in (f"v{version}", version):
            url = f"https://github.com/{repo}/archive/refs/tags/{tag}.tar.gz"
            try:
                req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "igvfagent"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    if r.status < 400:
                        return f"{name} @ {url}"
            except Exception:  # noqa: BLE001
                continue
    return None


def build_env_declared(yml: Path, inv: Dict[str, Any], pins: "Sequence[str]", run) -> Dict[str, Any]:
    """micromamba env from the authors' environment.yml. GPU-only packages are
    dropped on a CPU host; if the exact pins do not solve, pins are relaxed to
    major.minor, then to names (python kept at major.minor). pip dependencies
    install inside the env; a pin that fails is retried unpinned. Every change
    is recorded, so the report can say how this environment differs."""
    import yaml  # type: ignore
    doc = yaml.safe_load(yml.read_text()) or {}
    deps = doc.get("dependencies") or []
    conda = [d for d in deps if isinstance(d, str)]
    pip = [x for d in deps if isinstance(d, dict) for x in d.get("pip") or []]
    # --pin pip:NAME==VER overrides the declared pip dependency (e.g. the paper
    # states it used bean 0.2.9 while the repository's environment says 0.2.5)
    pip_over = {_dep_name(p_[4:])[0].lower(): p_[4:] for p_ in pins if p_.startswith("pip:")}
    explicit = set(pip_over.values())
    pip = [pip_over.pop(_dep_name(x)[0].lower(), x) for x in pip] + list(pip_over.values())
    pins = [p_ for p_ in pins if not p_.startswith("pip:")]
    gpu = bool(shutil.which("nvidia-smi"))
    dropped = [d for d in conda if GPU_ONLY.match(_dep_name(d)[0]) and not gpu]
    conda = [d for d in conda if d not in dropped and _dep_name(d)[0].lower() != "pip"]
    names = {_dep_name(d)[0].lower() for d in conda}
    extra = [x for x in ("pip", "c-compiler", "cxx-compiler", "make", "libxcrypt", "ipykernel", "nbconvert", "nbclient",
                         "nbformat") if x not in names]
    pyver = next((_dep_name(d)[1] for d in conda if _dep_name(d)[0].lower() == "python"), "")
    pymm = ".".join(pyver.split(".")[:2]) if pyver else ""

    def spec_at(level: int) -> "List[str]":
        out = []
        for d in conda:
            n, v = _dep_name(d)
            if level == 0 and v:
                out.append(f"{n}={v.split('=')[0]}")
            elif level == 1 and v:
                mm = ".".join(v.split("=")[0].split(".")[:2])
                out.append(f"{n}={mm}.*" if re.match(r"^\d+\.\d+", mm) else n)
            else:
                out.append(f"python={pymm}.*" if n.lower() == "python" and pymm else n)
        return out + extra + [p_ for p_ in pins if not p_.startswith("cran:")]

    # the builder's version is part of the key: a fix to how envs are built
    # (compilers, fallbacks) must not be hidden by a cached older env
    h = hashlib.sha256((DECLARED_ENV_BUILDER + yml.read_text() + "|".join(pins) + "|".join(pip)).encode()).hexdigest()[:12]
    prefix = CODE_DIR / "envs" / f"declared_{h}"
    ok_marker = prefix / ".igvf_env_ok.json"
    if ok_marker.exists():
        return {"ok": True, "prefix": str(prefix), "reused": True, **_read_json(ok_marker, {})}
    mm = micromamba()
    root = CODE_DIR / "mamba-root"
    level_used = None
    for level in (0, 1, 2):
        if prefix.exists():
            shutil.rmtree(prefix, ignore_errors=True)
        p_ = run([str(mm), "-r", str(root), "create", "-y", "-p", str(prefix), "-c", "conda-forge", "-c", "bioconda",
                  *spec_at(level)])
        if p_.returncode == 0:
            level_used = level
            break
    if level_used is None:
        return {"ok": False, "error": f"the declared environment {yml.name} could not be solved even with pins relaxed",
                "declared": _rel(yml)}
    py = prefix / "bin" / "python"
    pip_failed, pip_relaxed, conda_fallback, pin_failed = [], [], [], []
    # The env is never "activated", so point source builds at its own
    # compilers and headers (the system gcc lacks e.g. crypt.h for Python 3.8).
    benv = {**os.environ, "PATH": f"{prefix / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}",
            "CPATH": str(prefix / "include"), "LIBRARY_PATH": str(prefix / "lib")}
    for var, pats in (("CC", ("*-conda-linux-gnu-cc", "*-apple-darwin*-clang")),
                      ("CXX", ("*-conda-linux-gnu-c++", "*-apple-darwin*-clang++"))):
        hit = [h for pat in pats for h in sorted((prefix / "bin").glob(pat))]
        if hit:
            benv[var] = str(hit[0])
    # distutils links with Python's configured LDSHARED ("gcc -shared ..."),
    # not CC: the system gcc with the env's linker then looks for libc in the
    # wrong place (/lib64 on a Debian host). Link with the env's compiler too.
    if benv.get("CC"):
        rpath = f"-L{prefix / 'lib'} -Wl,-rpath,{prefix / 'lib'}"
        shared = "-dynamiclib -undefined dynamic_lookup" if sys.platform == "darwin" else "-pthread -shared"
        benv["LDSHARED"] = f"{benv['CC']} {shared} {rpath}"
        if benv.get("CXX"):
            benv["LDCXXSHARED"] = f"{benv['CXX']} {shared} {rpath}"
    _run0 = run

    def run(argv: "List[str]", **kw):  # pip gets the build environment
        if argv[:3] == [str(py), "-m", "pip"] and "env" not in kw:
            kw["env"] = benv
        return _run0(argv, **kw)
    if pip:
        if run([str(py), "-m", "pip", "install", "--no-input", *pip]).returncode:
            # Pinned build tools first (a Cython 0.29 package cannot build under
            # Cython 3), so the no-isolation retries below build against them.
            tools = [x for x in pip if _dep_name(x)[0].lower() in BUILD_TOOLS]
            for x in tools:
                run([str(py), "-m", "pip", "install", "--no-input", x])
            for x in pip:
                if x in tools:
                    continue
                if run([str(py), "-m", "pip", "install", "--no-input", x]).returncode == 0:
                    continue
                # Source packages built against the env's own pinned build
                # tools (e.g. Cython 0.29 for a package that breaks on Cython 3)
                if run([str(py), "-m", "pip", "install", "--no-input", "--no-build-isolation", x]).returncode == 0:
                    pip_relaxed.append(f"{x} (built without build isolation)")
                    continue
                n, v = _dep_name(x)
                if n.lower() in PIP_SKIP:
                    continue
                if x in explicit:
                    # A version the user or the paper asked for is never swapped
                    # for another: try the project's own release, else fail loudly.
                    src_spec = _release_archive(n, v) if v and "==" in x else None
                    if src_spec and run([str(py), "-m", "pip", "install", "--no-input", "--no-build-isolation",
                                         src_spec]).returncode == 0:
                        pip_relaxed.append(f"{x} (from the project's release archive: {src_spec.split(' @ ')[-1]})")
                    else:
                        pin_failed.append(x)
                    continue
                if _conda_install(run, mm, root, prefix, n):
                    conda_fallback.append(f"{x} (from conda-forge/bioconda)")
                elif v and run([str(py), "-m", "pip", "install", "--no-input", n]).returncode == 0:
                    pip_relaxed.append(f"{x} -> latest {n}")
                else:
                    pip_failed.append(x)
    # imports the notebooks use that the declared file does not cover
    have = (run([str(py), "-m", "pip", "list", "--format=freeze"]).stdout or "").lower()
    for m in inv.get("packages") or []:
        n = PY_PIP_NAME.get(m, m)
        if n.lower().replace("_", "-") + "==" not in have.replace("_", "-") and \
                run([str(py), "-c", f"import {m}"]).returncode:
            if run([str(py), "-m", "pip", "install", "--no-input", n]).returncode == 0:
                continue
            if _conda_install(run, mm, root, prefix, m):
                conda_fallback.append(f"{m} (imported by the notebooks; from conda-forge/bioconda)")
            else:
                pip_failed.append(f"{n} (imported by the notebooks)")
    vers = run([str(py), "-m", "pip", "freeze"]).stdout or ""
    pv = (run([str(py), "--version"]).stdout or "").strip()
    info = {"declared": _rel(yml), "relaxed_level": ["exact pins", "major.minor pins", "names only"][level_used],
            "dropped_gpu_only": dropped, "pip_failed": pip_failed, "pip_relaxed": pip_relaxed,
            "conda_fallback": conda_fallback, "pin_failed": pin_failed,
            "python": pv, "versions": (pv + "\n" + vers)[-8000:], "missing": pip_failed}
    # Rebuild next time only if a DECLARED dependency failed; modules the
    # notebooks import that no index has (the authors' private helpers) will
    # not appear on a rebuild either.
    if not [x for x in pip_failed if "(imported by the notebooks)" not in x] and not pin_failed:
        _write_json(ok_marker, info)
    if pin_failed:
        return {"ok": False, "prefix": None, "reused": False, **info,
                "error": f"requested version(s) could not be installed: {', '.join(pin_failed)} (not on PyPI and no "
                         "matching release archive was found; pin 'pip:NAME @ URL' to a source archive)"}
    return {"ok": True, "prefix": str(prefix), "reused": False, **info,
            "error": (f"pip could not install: {', '.join(pip_failed)}" if pip_failed else None)}


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
if ("jitter_global_rng" %%in%% shims) {
  # ggplot2 < 3.0 jittered with base::jitter() in the global random stream at
  # draw time; ggplot2 >= 3.0 takes one seed draw and jitters in a private
  # stream, which shifts every later sample()/runif() of a seeded analysis.
  # seed = NULL restores the global-stream behaviour. The Rmd is unchanged.
  prev <- fmt$knitr$knit_hooks$igvf_compat
  fmt$knitr$knit_hooks$igvf_compat <- function(before, options, envir) {
    if (is.function(prev)) prev(before, options, envir)
    if (before && "package:ggplot2" %%in%% search() && !("igvf_compat_jitter" %%in%% search())) {
      pj <- function(width = NULL, height = NULL, seed = NULL) ggplot2::position_jitter(width, height, seed = seed)
      gj <- function(mapping = NULL, data = NULL, stat = "identity", position = "jitter", ..., width = NULL,
                     height = NULL, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
        if (identical(position, "jitter")) position <- ggplot2::position_jitter(width, height, seed = NULL)
        ggplot2::geom_point(mapping = mapping, data = data, stat = stat, position = position, ...,
                            na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes)
      }
      attach(list(position_jitter = pj, geom_jitter = gj), name = "igvf_compat_jitter", warn.conflicts = FALSE)
    }
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
    "jitter_global_rng": "geom_jitter()/position_jitter() draw from the global random stream as in ggplot2 < 3.0 "
                         "(the authors' version), so later seeded draws follow the authors' sequence",
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
    gv = (inv.get("author_versions") or {}).get("ggplot2")
    if gv and re.match(r"^[0-2]\.", gv) and re.search(r"\b(geom_jitter|position_jitter)\s*\(", code) \
            and re.search(r"\bset\.seed\s*\(", code):
        out.append("jitter_global_rng")
    return out


def _runner_python() -> Path:
    """A Python with a current nbconvert + nbclient to drive notebooks:
    IGVF Agent's own interpreter when it has them, else a small cached venv."""
    try:
        import nbclient  # type: ignore  # noqa: F401
        import nbconvert  # type: ignore  # noqa: F401
        return Path(sys.executable)
    except Exception:  # noqa: BLE001
        pass
    venv = CODE_DIR / "envs" / "nb-runner"
    py = venv / "bin" / "python"
    if not py.exists():
        uv = shutil.which("uv")
        if uv:
            subprocess.run([uv, "venv", "--quiet", str(venv)], capture_output=True)
            subprocess.run([uv, "pip", "install", "--quiet", "--python", str(py), "nbconvert", "nbclient",
                            "ipykernel"], capture_output=True)
        else:
            subprocess.run([sys.executable, "-m", "venv", str(venv)], capture_output=True)
            subprocess.run([str(py), "-m", "pip", "install", "-q", "nbconvert", "nbclient", "ipykernel"],
                           capture_output=True)
    return py


_PY_SHIMS = r'''"""Compatibility shims for the authors' notebooks (IGVF Agent paper-code).

Loaded into the kernel at start; the notebooks are not modified. Each shim
restores something the authors' machine had and logs when it fires, so the
report says exactly what was substituted.
"""
import json, os, warnings

def _log(**kw):
    p = os.environ.get("IGVF_SHIM_LOG")
    if p:
        with open(p, "a") as fh:
            fh.write(json.dumps(kw) + "\n")

try:
    import matplotlib.style as _ms
    _orig_use = _ms.use

    def _use(style, *a, **k):
        # A personal style file of the authors (e.g. 'jr') that was never
        # published: fall back to matplotlib's default. Values are unaffected;
        # only the look of the figure differs.
        try:
            return _orig_use(style, *a, **k)
        except OSError:
            names = style if isinstance(style, (list, tuple)) else [style]
            _log(shim="mpl_missing_style", style=[str(n) for n in names])
            warnings.warn(f"IGVF Agent: matplotlib style {style!r} is not available here (the authors' own style "
                          "file); using the default style")
            return _orig_use("default", *a, **k)
    _ms.use = _use
    try:
        import matplotlib.style.core as _msc
        _msc.use = _use
    except Exception:
        pass
except Exception:
    pass


# gpu_to_cpu: code written for the authors' GPUs (--cuda, .cuda(), CUDA
# default tensor types) runs on a CPU-only host. Applied when torch is first
# imported and only if CUDA is unavailable; logged once.
import importlib.abc, importlib.machinery, sys


def _patch_torch(torch):
    try:
        if torch.cuda.is_available():
            return
    except Exception:
        pass
    _orig_sdtt = torch.set_default_tensor_type

    def _sdtt(t):
        name = t if isinstance(t, str) else getattr(t, "__module__", "") + "." + getattr(t, "__name__", "")
        if "cuda" in str(name):
            _log(shim="gpu_to_cpu", call="set_default_tensor_type", requested=str(name))
            cpu = str(name).replace("torch.cuda.", "torch.").replace("cuda.", "")
            return _orig_sdtt(cpu if isinstance(t, str) else getattr(torch, getattr(t, "__name__", "FloatTensor")))
        return _orig_sdtt(t)
    torch.set_default_tensor_type = _sdtt

    def _no_cuda(self, *a, **k):
        if not getattr(_no_cuda, "logged", False):
            _log(shim="gpu_to_cpu", call=".cuda()")
            _no_cuda.logged = True
        return self
    torch.Tensor.cuda = _no_cuda
    try:
        torch.nn.Module.cuda = _no_cuda
    except Exception:
        pass


class _TorchHook(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name != "torch":
            return None
        sys.meta_path.remove(self)
        spec = importlib.machinery.PathFinder.find_spec(name, path)
        if spec is None or spec.loader is None:
            return spec
        orig = spec.loader.exec_module

        def exec_module(module):
            orig(module)
            try:
                _patch_torch(module)
            except Exception:
                pass
        spec.loader.exec_module = exec_module
        return spec


if "torch" in sys.modules:
    _patch_torch(sys.modules["torch"])
elif not any(isinstance(f, _TorchHook) for f in sys.meta_path):
    sys.meta_path.insert(0, _TorchHook())
'''


def write_py_shims(d: Path, log_path: Path) -> Dict[str, str]:
    """The shim module plus a sitecustomize that loads it, for processes that
    are not notebook kernels (workflow rules, scripts). Returns env vars."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "igvf_py_shims.py").write_text(_PY_SHIMS)
    (d / "sitecustomize.py").write_text("try:\n    import igvf_py_shims  # noqa: F401\nexcept Exception:\n    pass\n")
    return {"PYTHONPATH": str(d), "IGVF_SHIM_LOG": str(log_path)}


_NB_DRIVER = r'''"""Execute one notebook with nbclient (IGVF Agent paper-code).

Reads leniently: old notebooks often lack metadata current nbformat requires
(kernelspec.display_name, cell ids); those are filled in memory only, the
file on disk is not touched. The executed copy is always written, also after
a timeout or a kernel death, so whatever ran is kept.
"""
import sys, os, warnings
import nbformat
from nbclient import NotebookClient

src, out, kname, timeout, tolerant = sys.argv[1:6]
warnings.filterwarnings("ignore")
with open(src) as fh:
    nb = nbformat.read(fh, as_version=4)
ks = nb.metadata.setdefault("kernelspec", {})
ks.setdefault("display_name", ks.get("name") or kname)
ks.setdefault("name", kname)
ks.setdefault("language", "python")
try:
    nbformat.validator.normalize(nb)
except Exception:
    pass
client = NotebookClient(nb, kernel_name=kname, timeout=int(timeout), allow_errors=tolerant == "1",
                        resources={"metadata": {"path": os.path.dirname(os.path.abspath(src))}})
rc = 0
try:
    client.execute()
except Exception as e:  # CellExecutionError (strict), CellTimeoutError, DeadKernelError
    print(f"notebook stopped: {type(e).__name__}: {str(e)[:2000]}", file=sys.stderr)
    rc = 1
finally:
    with open(out, "w") as fh:
        nbformat.write(nb, fh)
sys.exit(rc)
'''


def _register_kernel(run_dir: Path, prefix: Path, language: str) -> "Tuple[str, str]":
    """A kernelspec for the authors' environment, private to this run."""
    kdir = run_dir / "_jupyter" / "kernels" / "igvf-authors"
    kdir.mkdir(parents=True, exist_ok=True)
    if language == "r":
        spec = {"argv": [str(prefix / "bin" / "R"), "--slave", "-e", "IRkernel::main()", "--args", "{connection_file}"],
                "display_name": "R (authors' environment)", "language": "R"}
    else:
        shim_dir = run_dir / "_jupyter" / "shims"
        shim_dir.mkdir(parents=True, exist_ok=True)
        (shim_dir / "igvf_py_shims.py").write_text(_PY_SHIMS)
        spec = {"argv": [str(prefix / "bin" / "python"), "-m", "ipykernel_launcher", "-f", "{connection_file}",
                         "--IPKernelApp.exec_lines=['import igvf_py_shims']"],
                "display_name": "Python (authors' environment)", "language": "python",
                "env": {"PYTHONPATH": str(shim_dir), "IGVF_SHIM_LOG": str(run_dir / "_shim_log.jsonl")}}
    (kdir / "kernel.json").write_text(json.dumps(spec, indent=1))
    return "igvf-authors", str(run_dir / "_jupyter")


LINK_OVER_BYTES = int(os.environ.get("IGVF_PAPER_CODE_LINK_BYTES", str(2 * 1024 * 1024)))


def _link_or_copy(s_: str, d: str) -> str:
    """Large files (deposited data) are hard-linked, not copied: a 3.6 GB
    repository would otherwise be duplicated for every run. Analyses write new
    files; one that rewrote an input in place would change the checkout too,
    which the replay check and the pinned commit would then show."""
    try:
        if os.path.getsize(s_) > LINK_OVER_BYTES:
            os.link(s_, d)
            return d
    except OSError:
        pass
    return shutil.copy2(s_, d)


def prepare_work(src: Path, work: Path) -> None:
    shutil.copytree(src, work, ignore=shutil.ignore_patterns(".git"), copy_function=_link_or_copy, symlinks=True)


def _read_shim_log(p: Path) -> "List[Dict[str, Any]]":
    out: List[Dict[str, Any]] = []
    try:
        for ln in p.read_text().splitlines():
            e = json.loads(ln)
            if e not in out:
                out.append(e)
    except (OSError, ValueError):
        pass
    return out


def _snapshot(d: Path) -> Dict[str, float]:
    return {str(p.relative_to(d)): p.stat().st_mtime for p in d.rglob("*") if p.is_file()}


def run_entry(src: Path, inv: Dict[str, Any], env: Dict[str, Any], run_dir: Path,
              tolerant: bool = True, timeout_min: float = 180, work_name: str = "work",
              inputs: "Optional[Dict[str, str]]" = None, reuse_work: bool = False) -> Dict[str, Any]:
    work = run_dir / work_name
    if work.exists() and not reuse_work:
        shutil.rmtree(work) if not work.is_symlink() else work.unlink()
    if not work.exists():
        prepare_work(src, work)
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
        # The notebook is driven by IGVF Agent's own nbconvert/nbclient (the
        # authors' pinned nbconvert may not even import next to newer
        # dependencies), while its code runs in a kernel launched from the
        # authors' environment. Their kernelspec name (e.g. "jy_anbe_py38") is
        # their machine's; ours is registered per run. The notebook is unchanged.
        kname, jpath = _register_kernel(run_dir, prefix, inv["language"])
        envvars["JUPYTER_PATH"] = jpath + os.pathsep + envvars.get("JUPYTER_PATH", "")
        drv = run_dir / "_jupyter" / "igvf_nb_run.py"
        drv.write_text(_NB_DRIVER)
        argv = [str(_runner_python()), str(drv), entry.name, entry.stem + ".executed.ipynb", kname,
                str(int(timeout_min * 60)), "1" if tolerant else "0"]
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
            "r_shims": shims, "py_shims": _read_shim_log(run_dir / "_shim_log.jsonl"),
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
                    # the first error is the root cause; a later NameError is
                    # usually a name the failed cell would have defined
                    errors.append({"chunk": i + 1, "message": msg,
                                   "kind": "cascade" if errors and o.get("ename") == "NameError" else "root"})
                if o.get("output_type") == "display_data" and "image/png" in (o.get("data") or {}):
                    import base64
                    figs_dir.mkdir(parents=True, exist_ok=True)
                    pth = figs_dir / f"cell{i + 1}_{len(figures) + 1}.png"
                    pth.write_bytes(base64.b64decode(o["data"]["image/png"]))
                    figures.append({"path": pth, "source": "notebook", "png": pth, "chunk": i + 1})
    if not errors and run.get("exit_code") not in (0, None):
        executed = (entry.with_name(entry.stem + ".executed.ipynb") if entry.suffix == ".ipynb"
                    else _produced_md(work, inv) if entry.suffix.lower() in (".rmd", ".qmd") else None)
        if entry.suffix not in (".ipynb", ".Rmd", ".rmd", ".qmd") or executed is None or not executed.exists():
            # Nothing ran (or the script died): that is a failure, never "no errors".
            try:
                tail = [ln for ln in (ROOT / run.get("log", "")).read_text(errors="replace").splitlines() if ln.strip()]
            except OSError:
                tail = []
            msg = next((ln for ln in reversed(tail) if re.search(r"Error|error|Exception|No such|not found", ln)),
                       tail[-1] if tail else f"exit {run.get('exit_code')}")
            errors = [{"kind": "root", "chunk": None, "message": f"did not execute (exit {run.get('exit_code')}): "
                                                                 f"{msg.strip()[:300]}"}]
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
    if any(sh.get("shim") == "gpu_to_cpu" for sh in run.get("py_shims") or []):
        L.append("- Compatibility shim `gpu_to_cpu`: code written for the authors' GPUs ran on this CPU-only host "
                 "(CUDA tensor types and .cuda() mapped to CPU). Slower; results can differ in the last digits.")
    for sh in run.get("py_shims") or []:
        if sh.get("shim") == "mpl_missing_style":
            L.append(f"- Compatibility shim `mpl_missing_style`: matplotlib style {sh.get('style')} is the authors' own "
                     "and not published; figures use the default style (values unaffected).")
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


ANALYSIS_EXT = (".ipynb", ".Rmd", ".rmd", ".qmd")


# ─── workflow stage: the authors' pipeline before their notebooks ───────────

def find_workflow(src: Path) -> Optional[Dict[str, str]]:
    for rel in ("workflow/Snakefile", "Snakefile", "workflow/snakefile", "snakefile"):
        if (src / rel).is_file():
            return {"engine": "snakemake", "file": rel}
    for rel in ("main.nf", "workflow/main.nf"):
        if (src / rel).is_file():
            return {"engine": "nextflow", "file": rel}
    return None


def _missing_inputs(text: str) -> "List[str]":
    """Files Snakemake says are missing, from both its layouts (files listed
    directly under "Missing input files for rule X:", or under an indented
    "affected files:")."""
    out, on = [], False
    for ln in text.splitlines():
        t = ln.strip()
        if "Missing input files for rule" in ln:
            on = True
            continue
        if not on:
            continue
        if not t or not ln.startswith((" ", "\t")):
            on = False
            continue
        if re.match(r"^[A-Za-z][\w ]*:\s*\S*", t) and not re.match(r"^[\w./-]+$", t):
            continue  # "output: ...", "affected files:", "wildcards: ..."
        if re.match(r"^[\w./{}*-]+$", t):
            out.append(t)
    return sorted(set(out))[:200]


_MISSING_FILE = [re.compile(r"No such file or directory: '([^']+)'"), re.compile(r"name = '([^']+)'"),
                 re.compile(r"FileNotFoundError: \[Errno 2\] [^:]*: '([^']+)'"), re.compile(r"does not exist: '([^']+)'")]


def missing_files(sub: Path, entry: str, work: Path) -> "List[str]":
    """Workspace-relative files an analysis failed to open (from its root
    errors), resolved against the analysis's own directory."""
    sm = _read_json(sub / "summary.json", {}) or {}
    out = []
    for e in sm.get("errors") or []:
        for rx in _MISSING_FILE:
            for m in rx.finditer(e.get("message") or ""):
                q = Path(m.group(1))
                q = (q if q.is_absolute() else (work / Path(entry).parent / q)).resolve()
                try:
                    rel = str(q.relative_to(work.resolve()))
                except ValueError:
                    continue
                if rel not in out:
                    out.append(rel)
    return out


def _snk(exe: Path, work: Path, wf: Dict[str, str], wd: Path, targets: "Sequence[str]", *extra: str,
         cores: int = 4) -> "List[str]":
    """A snakemake command. Targets go first: in Snakemake 7 `--quiet` takes
    an optional value and would swallow a following target."""
    return [str(exe), *targets, "-s", str((work / wf["file"]).resolve()), "-d", str(wd), "--cores", str(cores),
            "--nolock", *extra]


def _workdirs(work: Path, wf: Dict[str, str]) -> "List[Path]":
    """Where the rules' relative paths live: the repository root, or the
    Snakefile's directory (WorkflowHub layouts keep results/ under workflow/)."""
    out = [work]
    d = (work / wf["file"]).parent
    if d != work:
        out.insert(0 if (d / "results").exists() or (d / "resources").exists() else 1, d)
    return out


def build_targets(work: Path, wf: Dict[str, str], env: Dict[str, Any], targets: "List[str]", log: Path,
                  timeout_min: float, cores: int) -> Dict[str, Any]:
    """Ask the workflow for specific files (what the analyses could not find):
    dry-run each, build every buildable one with --keep-going."""
    prefix = Path(env.get("prefix") or "")
    exe = prefix / "bin" / "snakemake"
    res: Dict[str, Any] = {"requested": targets, "buildable": [], "not_buildable": {}}
    if wf.get("engine") != "snakemake" or not exe.exists() or not targets:
        return res
    envvars = {**os.environ, "PATH": f"{prefix / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}",
               **write_py_shims(log.parent / "_jupyter" / "shims", log.parent / "_shim_log.jsonl")}
    by_wd: Dict[Path, List[str]] = {}
    for t in targets[:120]:
        why: List[str] = []
        for wd in _workdirs(work, wf):
            try:
                rel = os.path.relpath(work / t, wd)
            except ValueError:
                continue
            if rel.startswith(".."):
                continue
            d_ = subprocess.run(_snk(exe, work, wf, wd, [rel], "-n", "--quiet", "rules", cores=cores), cwd=str(wd),
                                env=envvars, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=600)
            txt = (d_.stdout or "") + (d_.stderr or "")
            if d_.returncode == 0:
                by_wd.setdefault(wd, []).append(rel)
                res["buildable"].append(t)
                break
            why = _missing_inputs(txt)[:3] or [next((ln.strip() for ln in txt.splitlines()
                                                    if "MissingRuleException" in ln or "No rule to produce" in ln),
                                                   txt.strip().splitlines()[-1][:160] if txt.strip() else "")]
        else:
            res["not_buildable"][t] = why
    if by_wd:
        t0 = time.time()
        codes = []
        _build_rounds(res, by_wd, exe, work, wf, cores, envvars, log, timeout_min, prefix, codes)
        res["exit"] = max(codes) if codes else None
        res["workdirs"] = [_rel(w) for w in by_wd]
        res["seconds"] = round(time.time() - t0, 1)
        res["built"] = [t for t in res["buildable"] if (work / t).exists()]
    return res


def _pip(py: Path, *args: str) -> "List[str]":
    """pip for an environment's interpreter: its own pip module, else uv."""
    has = subprocess.run([str(py), "-c", "import pip"], capture_output=True, stdin=subprocess.DEVNULL).returncode == 0
    if has:
        return [str(py), "-m", "pip", *args[:1], *(["--no-input"] if args[:1] == ("install",) else []), *args[1:]]
    uv = shutil.which("uv") or "uv"
    return [uv, "pip", *args[:1], "--python", str(py), *args[1:]]


def _build_rounds(res: Dict[str, Any], by_wd: Dict[Path, List[str]], exe: Path, work: Path, wf: Dict[str, str],
                  cores: int, envvars: Dict[str, str], log: Path, timeout_min: float, prefix: Path,
                  codes: "List[int]") -> None:
    """Build the targets; if a rule died on a module the declared environment
    lacks (its authors had it installed without declaring it), install that
    module (pip, then conda-forge/bioconda) and build once more."""
    got_all: List[str] = []
    for attempt in (1, 2, 3, 4):
        start = log.stat().st_size if log.exists() else 0
        round_codes: List[int] = []
        for wd, rels in by_wd.items():
            cmd = _snk(exe, work, wf, wd, rels, "--keep-going", "--rerun-incomplete", cores=cores)
            with open(log, "a") as lf:
                lf.write(f"\n---- targets the analyses asked for (attempt {attempt}) ----\n$ " + " ".join(cmd) + "\n")
                lf.flush()
                try:
                    round_codes.append(subprocess.run(cmd, cwd=str(wd), env=envvars, stdout=lf,
                                                      stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                                      timeout=timeout_min * 60).returncode)
                except subprocess.TimeoutExpired:
                    round_codes.append(124)
        codes[:] = round_codes  # the result is the last round's
        if attempt == 4 or not any(round_codes):
            return
        with open(log, errors="replace") as fh:
            fh.seek(start)
            mods = sorted(set(re.findall(r"ModuleNotFoundError: No module named '([\w.]+)'", fh.read())))
        mods = [m.split(".")[0] for m in mods]
        if not mods:
            return
        py = prefix / "bin" / "python"
        got = []
        # Never upgrade what the authors pinned: everything installed now is a
        # constraint, so pip picks a version compatible with it (or fails).
        # `pip list --format=freeze` names conda-installed packages too (plain
        # `pip freeze` prints them as "torch @ file:///...", which a
        # constraints file cannot use; leaving them out let pip upgrade the
        # authors' conda torch 1.12 to 2.x).
        frz = subprocess.run(_pip(py, "list", "--format=freeze"), capture_output=True, text=True,
                             stdin=subprocess.DEVNULL).stdout or ""
        cons = log.parent / "_constraints.txt"
        cons.write_text("\n".join(ln.split("+")[0] if "==" in ln else ln for ln in frz.splitlines()
                                  if re.match(r"^[A-Za-z0-9_.\-]+==[\w.\-+!]+$", ln.strip())) + "\n")
        for m in sorted(set(mods)):
            pipn = PY_PIP_NAME.get(m, m)
            if subprocess.run(_pip(py, "install", "-c", str(cons), pipn), capture_output=True,
                              stdin=subprocess.DEVNULL).returncode == 0:
                got.append(f"{m} (pip {pipn})")
                continue
            mm = micromamba()
            for cand in CONDA_NAMES.get(m.lower(), [pipn.lower(), m.lower()]):
                if subprocess.run([str(mm), "-r", str(CODE_DIR / "mamba-root"), "install", "-y", "-p", str(prefix),
                                   "-c", "conda-forge", "-c", "bioconda", cand], capture_output=True,
                                  stdin=subprocess.DEVNULL).returncode == 0:
                    got.append(f"{m} (conda {cand})")
                    break
        got_all += got
        res["installed_on_demand"] = got_all
        if not got:
            return


def run_workflow(work: Path, wf: Dict[str, str], env: Dict[str, Any], log: Path, timeout_min: float,
                 cores: int) -> Dict[str, Any]:
    """Dry run (the plan and every missing input), then the real run with
    --keep-going, in the work copy. Rules that declare their own conda envs run
    in the declared environment instead (recorded)."""
    res: Dict[str, Any] = {"engine": wf["engine"], "file": wf["file"]}
    if wf["engine"] != "snakemake":
        res["error"] = f"{wf['engine']} workflows are detected but not run by paper-code yet"
        return res
    prefix = Path(env.get("prefix") or "")
    exe = prefix / "bin" / "snakemake"
    if not exe.exists():
        res["error"] = "snakemake is not in the analysis environment (add it to the declared env or --pin snakemake)"
        return res
    envvars = {**os.environ, "PATH": f"{prefix / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}",
               **write_py_shims(log.parent / "_jupyter" / "shims", log.parent / "_shim_log.jsonl")}
    wd = _workdirs(work, wf)[0]
    res["workdir"] = _rel(wd)
    base = _snk(exe, work, wf, wd, [], cores=cores)
    t0 = time.time()
    dry = subprocess.run(base + ["-n", "--quiet", "rules"], cwd=str(wd), env=envvars, capture_output=True, text=True,
                         stdin=subprocess.DEVNULL, timeout=1800)
    txt = (dry.stdout or "") + (dry.stderr or "")
    jobs = re.findall(r"^\s*(\w[\w.-]*)\s+(\d+)\s*$", txt, re.M)
    res["dry_run"] = {"exit": dry.returncode, "jobs_by_rule": {r: int(n) for r, n in jobs if r not in ("total", "job")},
                      "total_jobs": next((int(n) for r, n in jobs if r == "total"), None),
                      "missing_inputs": _missing_inputs(txt),
                      "conda_rules": len(re.findall(r"^\s*conda\s*:", (work / wf["file"]).read_text(errors="replace"), re.M))}
    targets: List[str] = []
    if dry.returncode != 0 and res["dry_run"]["missing_inputs"]:
        # The full DAG needs files that are not here (raw reads, usually).
        # Run every rule that can be built from what exists instead: often the
        # modelling steps start from the deposited processed objects.
        lr = subprocess.run(base + ["--list-rules"], cwd=str(wd), env=envvars,
                            capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=600)
        rules = [r.strip() for r in (lr.stdout or "").splitlines() if re.match(r"^\s*\w+\s*$", r)]
        blocked: Dict[str, List[str]] = {}
        for r in rules:
            if r == "all":
                continue
            d_ = subprocess.run(_snk(exe, work, wf, wd, [r], "-n", "--quiet", "rules", cores=cores), cwd=str(wd),
                                env=envvars, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=600)
            if d_.returncode == 0:
                targets.append(r)
            else:
                blocked[r] = _missing_inputs((d_.stdout or "") + (d_.stderr or ""))[:5]
        res["targets"] = {"buildable": targets, "blocked": blocked}
    before = _snapshot(work)
    cmd = _snk(exe, work, wf, wd, targets, "--keep-going", "--rerun-incomplete", cores=cores)
    with open(log, "w") as lf:
        lf.write("$ " + " ".join(cmd) + "\n" + txt + "\n---- run ----\n")
        lf.flush()
        if dry.returncode != 0 and not targets:
            rc = dry.returncode  # nothing can be built from what is here
        else:
            try:
                p_ = subprocess.run(cmd, cwd=str(wd), env=envvars, stdout=lf, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, timeout=timeout_min * 60)
                rc = p_.returncode
            except subprocess.TimeoutExpired:
                rc = 124
    out = log.read_text(errors="replace")
    after = _snapshot(work)
    res["run"] = {"exit": rc, "seconds": round(time.time() - t0, 1),
                  "failed_rules": sorted(set(re.findall(r"Error in rule (\S+?):", out)))[:50],
                  "steps_done": len(re.findall(r"\d+ of \d+ steps \(\d+%\) done", out)),
                  "new_files": sum(1 for k, t in after.items() if k not in before or t > before[k] + 1e-6),
                  "log": _rel(log)}
    return res


# ─── raw reads (SRA / ENA) into the paths the workflow reads ─────────────────

ENA_FIELDS = "run_accession,sample_alias,sample_title,experiment_title,library_name,experiment_alias,fastq_bytes,fastq_ftp,fastq_md5"


def _http_json(url: str) -> Any:
    req = urllib.request.Request(url, headers={"User-Agent": "igvfagent"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def resolve_study(acc: str) -> str:
    """GEO series -> its BioProject (ENA knows BioProjects, SRA studies and
    runs directly)."""
    acc = acc.strip()
    if acc.upper().startswith("GSE"):
        es = _http_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term={acc}"
                        "[All Fields]&retmode=json")
        ids = es.get("esearchresult", {}).get("idlist") or []
        if ids:
            sm = _http_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id={ids[0]}"
                            "&retmode=json")
            return sm["result"][ids[0]]["project_acc"]
        raise ValueError(f"no BioProject found for {acc}")
    return acc


def ena_runs(acc: str) -> "List[Dict[str, str]]":
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={urllib.parse.quote(resolve_study(acc))}"
           f"&result=read_run&fields={ENA_FIELDS}&format=tsv&limit=0")
    req = urllib.request.Request(url, headers={"User-Agent": "igvfagent"})
    with urllib.request.urlopen(req, timeout=120) as r:
        lines = r.read().decode().splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def guess_read_layout(src: Path) -> "List[str]":
    """FASTQ path patterns the authors' code builds (to choose --layout)."""
    pats = []
    for f in list(src.rglob("*.py")) + list(src.rglob("*.smk")) + list(src.rglob("Snakefile")) + \
            list(src.rglob("*.R")) + list(src.rglob("*.sh")) + list(src.rglob("*.yaml")) + list(src.rglob("*.yml")):
        if ".git" in f.parts or f.stat().st_size > 2_000_000:
            continue
        for m in re.finditer(r"[\"']([^\"'\n]*\{[^\"'\n]*\}[^\"'\n]*\.(?:fastq|fq)(?:\.gz)?)[\"']",
                             f.read_text(errors="replace")):
            k = f"{m.group(1)}    ({f.relative_to(src)})"
            if k not in pats:
                pats.append(k)
    return pats[:30]


def plan_reads(runs: "List[Dict[str, str]]", layout: str, match: Optional[str] = None, field: str = "library_name",
               include: Optional[str] = None) -> "List[Dict[str, Any]]":
    rx = re.compile(match) if match else None
    inc = re.compile(include) if include else None
    plan = []
    for r in runs:
        key = r.get(field) or ""
        if inc and not inc.search(key):
            continue
        fields = {k: v for k, v in r.items()}
        if rx:
            m = rx.search(key)
            if not m:
                plan.append({"run": r["run_accession"], "key": key, "skip": f"--match did not match {field}"})
                continue
            fields.update({k: v for k, v in m.groupdict().items() if v is not None})
        urls = [u for u in (r.get("fastq_ftp") or "").split(";") if u]
        md5s = (r.get("fastq_md5") or "").split(";")
        sizes = (r.get("fastq_bytes") or "").split(";")
        paired = [i for i, u in enumerate(urls) if re.search(r"_[12]\.f(ast)?q\.gz$", u)]
        use = paired or list(range(len(urls)))
        for n, i in enumerate(use, 1):
            read = re.search(r"_([12])\.f(ast)?q\.gz$", urls[i]).group(1) if paired else str(n)
            try:
                dest = layout.format(**fields, read=read)
            except KeyError as e:
                plan.append({"run": r["run_accession"], "key": key, "skip": f"layout needs {e}"})
                break
            plan.append({"run": r["run_accession"], "key": key, "read": read, "url": "https://" + urls[i],
                         "md5": md5s[i] if i < len(md5s) else "", "bytes": int(sizes[i]) if i < len(sizes) and
                         sizes[i].isdigit() else 0, "dest": dest})
    return plan


def fetch_reads(repo: str, accession: str, layout: str, *, into: str = "", match: Optional[str] = None,
                field: str = "library_name", include: Optional[str] = None, download: bool = False,
                max_gb: float = 200.0, workers: int = 4) -> Dict[str, Any]:
    """Map a study's runs onto the file names the authors' code reads, and
    (with download) fetch them from ENA with md5 checks, resumably. Without
    download it only previews the mapping, so the layout can be checked first."""
    d = repo_dir(repo)
    src = d / "src"
    if not src.is_dir():
        return {"ok": False, "error": f"fetch the repository first (paper-code fetch {repo})"}
    base = (src / into).resolve()
    runs = ena_runs(accession)
    plan = plan_reads(runs, layout, match, field, include)
    todo = [p_ for p_ in plan if not p_.get("skip")]
    for p_ in todo:
        q = (base / p_["dest"]).resolve()
        try:
            q.relative_to(src.resolve())
        except ValueError:
            return {"ok": False, "error": f"layout places {p_['dest']} outside the repository"}
        p_["path"] = q
    dests = [p_["dest"] for p_ in todo]
    dup = sorted({x for x in dests if dests.count(x) > 1})
    total = sum(p_["bytes"] for p_ in todo)
    out: Dict[str, Any] = {"ok": True, "accession": accession, "runs": len(runs), "files": len(todo),
                           "gigabytes": round(total / 1e9, 2), "skipped": [p_ for p_ in plan if p_.get("skip")][:20],
                           "duplicates": dup[:20],
                           "preview": [{"run": p_["run"], "read": p_["read"], "dest": p_["dest"]} for p_ in todo[:12]]}
    if dup:
        out.update(ok=False, error=f"{len(dup)} destination(s) would receive several runs; refine --layout/--match")
        return out
    if not download:
        out["note"] = "preview only: re-run with --yes to download"
        return out
    if total / 1e9 > max_gb:
        return {**out, "ok": False, "error": f"{total / 1e9:.1f} GB exceeds --max-gb {max_gb}"}
    import concurrent.futures

    def one(p_: Dict[str, Any]) -> Dict[str, Any]:
        q: Path = p_["path"]
        if q.exists() and p_["bytes"] and q.stat().st_size == p_["bytes"]:
            return {"dest": p_["dest"], "status": "present"}
        q.parent.mkdir(parents=True, exist_ok=True)
        tmp = q.with_name(q.name + ".part")
        h = hashlib.md5()
        for attempt in range(3):
            try:
                req = urllib.request.Request(p_["url"], headers={"User-Agent": "igvfagent"})
                with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as fh:
                    h = hashlib.md5()
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        h.update(b)
                        fh.write(b)
                break
            except Exception as e:  # noqa: BLE001
                if attempt == 2:
                    return {"dest": p_["dest"], "status": f"failed: {e}"}
        if p_["md5"] and h.hexdigest() != p_["md5"]:
            tmp.unlink(missing_ok=True)
            return {"dest": p_["dest"], "status": "md5 mismatch"}
        os.replace(tmp, q)
        return {"dest": p_["dest"], "status": "downloaded"}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(one, todo))
    ok_n = sum(1 for r in results if r["status"] in ("downloaded", "present"))
    prov = _read_json(d / "data.json", []) or []
    prov.append({"source": f"ENA/SRA {accession}", "layout": layout, "into": into, "match": match,
                 "files": ok_n, "gigabytes": out["gigabytes"],
                 "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                 "file": _rel(base), "checksum": "md5 per file (ENA)"})
    _write_json(d / "data.json", prov)
    manifest = d / f"reads_{re.sub(r'[^A-Za-z0-9]+', '_', accession)}.tsv"
    manifest.write_text("run\tread\tdest\tstatus\n" + "".join(
        f"{p_['run']}\t{p_['read']}\t{p_['dest']}\t{r['status']}\n" for p_, r in zip(todo, results)))
    out.update({"downloaded_or_present": ok_n, "failed": [r for r in results if r["status"] not in
                                                          ("downloaded", "present")][:20],
                "manifest": _rel(manifest), "ok": ok_n == len(todo)})
    return out


# ─── deposited data (Zenodo, URLs) into the checkout ────────────────────────

def _zenodo_id(ref: str) -> Optional[str]:
    m = re.search(r"zenodo\.(\d+)", ref) or re.search(r"records?/(\d+)", ref) or re.fullmatch(r"\s*(\d+)\s*", ref)
    return m.group(1) if m else None


def fetch_data(repo: str, *, zenodo: Optional[str] = None, url: Optional[str] = None, into: str = "",
               files: Optional[str] = None, extract: bool = False) -> Dict[str, Any]:
    """Download a paper's deposited data into its checkout (so every later run
    sees it), verify checksums, optionally extract archives, and record the
    provenance in Data/PaperCode/<repo>/data.json."""
    import fnmatch
    d = repo_dir(repo)
    src = d / "src"
    if not src.is_dir():
        return {"ok": False, "error": f"fetch the repository first (paper-code fetch {repo})"}
    dest = (src / into).resolve()
    try:
        dest.relative_to(src.resolve())
    except ValueError:
        return {"ok": False, "error": "--into must stay inside the repository"}
    dest.mkdir(parents=True, exist_ok=True)
    items: List[Dict[str, Any]] = []
    if zenodo:
        rid = _zenodo_id(zenodo)
        if not rid:
            return {"ok": False, "error": f"not a Zenodo record or DOI: {zenodo}"}
        req = urllib.request.Request(f"https://zenodo.org/api/records/{rid}", headers={"User-Agent": "igvfagent"})
        with urllib.request.urlopen(req, timeout=60) as r:
            rec = json.loads(r.read().decode())
        for f in rec.get("files") or []:
            name = f.get("key") or f.get("filename")
            if files and not fnmatch.fnmatch(name, files):
                continue
            items.append({"name": name, "url": (f.get("links") or {}).get("self") or (f.get("links") or {}).get("download"),
                          "checksum": f.get("checksum"), "size": f.get("size"),
                          "source": f"zenodo:{rec.get('id')} ({rec.get('metadata', {}).get('title', '')[:80]})"})
    elif url:
        items.append({"name": Path(urllib.parse.urlparse(url).path).name or "download", "url": url, "checksum": None,
                      "size": None, "source": url})
    else:
        return {"ok": False, "error": "give --zenodo or --url"}
    done = []
    for it in items:
        out = dest / it["name"]
        if not (out.exists() and it.get("size") and out.stat().st_size == it["size"]):
            req = urllib.request.Request(it["url"], headers={"User-Agent": "igvfagent"})
            h = hashlib.md5()
            tmp = out.with_suffix(out.suffix + ".part")
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as fh:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    h.update(b)
                    fh.write(b)
            want = (it.get("checksum") or "").split(":", 1)[-1] if (it.get("checksum") or "").startswith("md5") else None
            if want and h.hexdigest() != want:
                tmp.unlink()
                return {"ok": False, "error": f"checksum mismatch for {it['name']}"}
            os.replace(tmp, out)
        entry = {"file": _rel(out), "source": it["source"], "checksum": it.get("checksum"),
                 "bytes": out.stat().st_size}
        if extract and re.search(r"\.(zip|tar\.gz|tgz|tar)$", out.name):
            if out.name.endswith(".zip"):
                import zipfile
                with zipfile.ZipFile(out) as z:
                    for m in z.namelist():
                        if m.startswith("/") or ".." in Path(m).parts:
                            continue
                        z.extract(m, dest)
            else:
                with tarfile.open(out) as tf:
                    for m in tf.getmembers():
                        if m.name.startswith("/") or ".." in Path(m.name).parts or m.issym() or m.islnk():
                            continue
                        tf.extract(m, dest)
            entry["extracted_into"] = _rel(dest)
        done.append(entry)
    prov = _read_json(d / "data.json", []) or []
    prov += [{**e, "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())} for e in done]
    _write_json(d / "data.json", prov)
    return {"ok": True, "files": done, "into": _rel(dest)}


def _slug_path(p: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", p).strip("_")[:60]


def _merge_invs(invs: "List[Dict[str, Any]]") -> Dict[str, Any]:
    """One environment per language for all the notebooks that use it."""
    m = dict(invs[0])
    m["packages"] = sorted({p_ for i in invs for p_ in i.get("packages") or []})
    m["r_notebooks"] = any(i["entry"].endswith(".ipynb") for i in invs if i["language"] == "r")
    for i in invs:
        for k, v in (i.get("author_versions") or {}).items():
            m.setdefault("author_versions", {}).setdefault(k, v)
    return m


def _sum_summaries(entries: "List[Dict[str, Any]]") -> Dict[str, Any]:
    tot = {"chunks": {"total": 0, "errored": 0, "root_errors": 0, "cascade_errors": 0},
           "figures": {"produced": 0, "reference": 0, "saved_files": 0},
           "printed": {"reference_blocks": 0, "identical": 0, "numerically_close": 0,
                       "same_values_other_layout": 0, "missing": 0}}
    for e in entries:
        for sect in ("chunks", "figures", "printed"):
            for k in tot[sect]:
                tot[sect][k] += int((e.get(sect) or {}).get(k) or 0)
    pr = tot["printed"]
    got = pr["identical"] + pr["numerically_close"] + pr["same_values_other_layout"]
    pr["fraction_matched"] = round(got / pr["reference_blocks"], 4) if pr["reference_blocks"] else None
    return tot


def _rewrite_links(md: str, prefix: str) -> str:
    return re.sub(r"(!?\[[^\]]*\]\()(?!https?:|data:|#|/)([^)]+)\)", lambda m: f"{m.group(1)}{prefix}{m.group(2)})", md)


def pipeline_multi(repo: str, entries: "List[str]", *, ref: Optional[str], pins: "Sequence[str]", strict: bool,
                   paper: Optional[str], run_dir: Path, entry_timeout_min: float, replay: bool,
                   inputs: "Optional[Dict[str, str]]", workflow: bool = False, workflow_timeout_min: float = 240,
                   cores: int = 4, max_passes: int = 2) -> Dict[str, Any]:
    """Every selected notebook of a repository, in one shared work copy and in
    order (as the authors ran them), each with its own sub-report; one
    combined summary and report for the paper."""
    status = run_dir / "status.json"

    def stage(name: str, **kw) -> None:
        _write_json(status, {"state": "running", "stage": name, "at": time.time(), **kw})
        print(f"[paper-code] {name}", flush=True)

    src_meta = fetch(repo, ref)
    if not src_meta.get("ok"):
        return {"ok": False, "error": src_meta["error"]}
    src = ROOT / src_meta["src"]
    stage("inventory", entries=len(entries))
    invs = []
    for e in entries:
        i = inventory(src, e)
        if i.get("ok"):
            invs.append(i)
    if not invs:
        return {"ok": False, "error": "none of the selected entries could be inventoried"}
    envs: Dict[str, Dict[str, Any]] = {}
    for lang in sorted({i["language"] for i in invs}):
        merged = _merge_invs([i for i in invs if i["language"] == lang])
        stage("env", language=lang, packages=len(merged["packages"]))
        envs[lang] = build_env(merged, pins, log=run_dir / f"env_{lang}.log")
        _write_json(run_dir / f"env_{lang}.json", envs[lang])
    prepare_work(src, run_dir / "work")
    wf_res = None
    wf = find_workflow(src) if workflow else None
    if wf:
        stage("workflow", engine=wf["engine"])
        wenv = envs.get("python") or next(iter(envs.values()))
        wf_res = run_workflow(run_dir / "work", wf, wenv, run_dir / "workflow.log", workflow_timeout_min, cores)
        _write_json(run_dir / "workflow.json", wf_res)
    elif workflow:
        wf_res = {"error": "no Snakefile or main.nf in the repository"}
    if replay:
        prepare_work(src, run_dir / "replay")
    summaries = []
    for k, inv_e in enumerate(invs, 1):
        sub = run_dir / "entries" / f"{k:02d}_{_slug_path(inv_e['entry'])}"
        sub.mkdir(parents=True, exist_ok=True)
        for nm in ("work",) + (("replay",) if replay else ()):
            link = sub / nm
            if not link.exists():
                link.symlink_to(os.path.relpath(run_dir / nm, sub))
        env = envs[inv_e["language"]]
        _write_json(sub / "inventory.json", inv_e)
        _write_json(sub / "env.json", env)
        stage("run", entry=inv_e["entry"], k=k, n=len(invs))
        if not env.get("prefix"):
            summaries.append({"entry": inv_e["entry"], "dir": _rel(sub), "render": "no environment",
                              "error": env.get("error")})
            continue
        mine = {n: v for n, v in (inputs or {}).items() if n in {i["path"] for i in inv_e.get("inputs") or []}}
        run = run_entry(src, inv_e, env, sub, tolerant=not strict, timeout_min=entry_timeout_min, inputs=mine,
                        reuse_work=True)
        _write_json(sub / "run.json", run)
        if replay:
            _write_json(sub / "replay.json", run_entry(src, inv_e, env, sub, tolerant=not strict,
                                                       timeout_min=entry_timeout_min, work_name="replay",
                                                       inputs=mine, reuse_work=True))
            _write_json(sub / "replay_compare.json", replay_compare(sub, inv_e))
        rep_e = write_report(sub, inv_e, env, src_meta, paper)
        sm = rep_e["summary"]
        roots = [e for e in sm.get("errors") or [] if e.get("kind") != "cascade"]
        summaries.append({"entry": inv_e["entry"], "dir": _rel(sub), "render": sm.get("render"),
                          "chunks": sm.get("chunks"), "figures": sm.get("figures"), "printed": sm.get("printed"),
                          "has_reference": bool(inv_e.get("reference")), "inputs_missing": inv_e.get("inputs_missing"),
                          "first_error": roots[0]["message"][:240] if roots else None,
                          "replay": sm.get("replay"), "seconds": sm.get("seconds")})
    # The analyses name what they could not find. The workflow may produce it
    # (e.g. model runs from the deposited processed objects), and another
    # analysis may have written it in this pass (one notebook's table is the
    # next one's input). Build what the workflow can, then re-run every
    # analysis whose missing inputs now exist, until nothing changes.
    passes: List[Dict[str, Any]] = []
    for it in range(1 + max_passes):
        need: Dict[str, List[str]] = {}
        for x in summaries:
            for f_ in missing_files(ROOT / x["dir"], x["entry"], run_dir / "work"):
                need.setdefault(f_, []).append(x["entry"])
        if not need:
            break
        info: Dict[str, Any] = {"pass": it + 2, "missing": len(need)}
        if it == 0 and wf and wf.get("engine") == "snakemake" and not (wf_res or {}).get("error"):
            stage("workflow-targets", files=len(need))
            wenv = envs.get("python") or next(iter(envs.values()))
            bt = build_targets(run_dir / "work", wf, wenv, sorted(need), run_dir / "workflow.log",
                               workflow_timeout_min, cores)
            (wf_res or {}).update({"second_pass": bt})
            _write_json(run_dir / "workflow.json", wf_res)
            info["workflow_built"] = len(bt.get("built") or [])
        now = [f_ for f_ in need if (run_dir / "work" / f_).exists()]
        rerun = sorted({e for f_ in now for e in need[f_]})
        info.update({"now_present": len(now), "rerun": rerun})
        passes.append(info)
        if not rerun:
            break
        for k, inv_e in enumerate(invs, 1):
            if inv_e["entry"] not in rerun:
                continue
            sub = run_dir / "entries" / f"{k:02d}_{_slug_path(inv_e['entry'])}"
            env = envs[inv_e["language"]]
            stage("rerun", entry=inv_e["entry"], pass_=it + 2)
            run = run_entry(src, inv_e, env, sub, tolerant=not strict, timeout_min=entry_timeout_min,
                            reuse_work=True)
            _write_json(sub / "run.json", run)
            sm = write_report(sub, inv_e, env, src_meta, paper)["summary"]
            roots = [e for e in sm.get("errors") or [] if e.get("kind") != "cascade"]
            for x in summaries:
                if x["entry"] == inv_e["entry"]:
                    x.update({"render": sm.get("render"), "chunks": sm.get("chunks"), "figures": sm.get("figures"),
                              "printed": sm.get("printed"), "first_error": roots[0]["message"][:240] if roots else None,
                              "passes": x.get("passes", 1) + 1})
    stage("report")
    tot = _sum_summaries([x for x in summaries if x.get("chunks")])
    clean = sum(1 for x in summaries if x.get("chunks") and not x["chunks"].get("root_errors"))
    summary = {"paper": paper, "repo": src_meta.get("repo"), "commit": src_meta.get("commit"),
               "entry": f"{len(summaries)} analyses", "language": "+".join(sorted(envs)),
               "render": f"{clean}/{len(summaries)} ran without root errors", **tot,
               "entries": summaries, "env": {lang: {k: e.get(k) for k in ("prefix", "declared", "relaxed_level",
                                                                          "dropped_gpu_only", "pip_failed", "missing")}
                                             for lang, e in envs.items()},
               "replay": ({"deterministic": all((x.get("replay") or {}).get("deterministic") for x in summaries
                                                if x.get("replay"))} if replay else None),
               "workflow": wf_res, "data": _read_json(repo_dir(repo) / "data.json"), "passes": passes,
               "still_missing": sorted({f_ for x in summaries for f_ in missing_files(ROOT / x["dir"], x["entry"],
                                                                                      run_dir / "work")})[:300]}
    _write_json(run_dir / "summary.json", summary)
    pr = tot["printed"]
    L = [f"# Reproduction from the authors' code: {src_meta.get('repo')}", ""]
    if paper:
        L += [f"**Paper:** {paper}", ""]
    L += [f"- Repository: [{src_meta.get('repo')}]({src_meta.get('url')}) at commit `{(src_meta.get('commit') or '')[:12]}`",
          f"- **{len(summaries)} analyses** run unmodified in one work copy, in order: **{clean}** without root-cause "
          f"errors; {tot['chunks']['total']} chunks/cells, {tot['chunks']['root_errors']} root errors, "
          f"{tot['chunks']['cascade_errors']} knock-on",
          f"- Printed values vs the authors' stored outputs: **{pr['identical'] + pr['numerically_close'] + pr['same_values_other_layout']}"
          f" of {pr['reference_blocks']}** ({(pr['fraction_matched'] or 0) * 100:.1f}%)",
          f"- Figures: **{tot['figures']['produced']}** produced (the authors' outputs have {tot['figures']['reference']})"]
    for lang, e in envs.items():
        if e.get("declared"):
            L.append(f"- {lang} environment: the authors' `{e['declared']}` ({e.get('relaxed_level')})"
                     + (f"; GPU-only packages dropped on this CPU host: {', '.join(e['dropped_gpu_only'])}"
                        if e.get("dropped_gpu_only") else "")
                     + (f"; pip could not install: {', '.join(e['pip_failed'])}" if e.get("pip_failed") else ""))
    if wf_res:
        if wf_res.get("error"):
            L.append(f"- Workflow: {wf_res['error']}")
        else:
            dr, rn = wf_res.get("dry_run") or {}, wf_res.get("run") or {}
            L.append(f"- Workflow `{wf_res['file']}` ({wf_res['engine']}) ran first: {dr.get('total_jobs')} planned jobs, "
                     f"{rn.get('steps_done')} steps done in {rn.get('seconds')} s, exit {rn.get('exit')}, "
                     f"{rn.get('new_files')} new files"
                     + (f"; failed rules: {', '.join(rn['failed_rules'][:8])}" if rn.get("failed_rules") else "")
                     + (f"; ran the {len(wf_res['targets']['buildable'])} rule(s) buildable from the files present "
                        f"({', '.join(wf_res['targets']['buildable'][:8])}); "
                        f"{len(wf_res['targets']['blocked'])} rule(s) need inputs that are not here"
                        if wf_res.get("targets") else "")
                     + (f"; the analyses asked for {len(wf_res['second_pass']['requested'])} missing file(s), "
                        f"{len(wf_res['second_pass']['buildable'])} buildable by the workflow, "
                        f"{len(wf_res['second_pass'].get('built') or [])} built"
                        + (f" (installed on demand, missing from the declared environment: "
                           f"{', '.join(wf_res['second_pass']['installed_on_demand'])})"
                           if wf_res['second_pass'].get('installed_on_demand') else "")
                        if wf_res.get("second_pass") else "")
                     + (f"; **{len(dr['missing_inputs'])} missing inputs** (deposited data to fetch), e.g. "
                        f"`{dr['missing_inputs'][0]}`" if dr.get("missing_inputs") else ""))
    wf_shims = _read_shim_log(run_dir / "_shim_log.jsonl")
    if any(x.get("shim") == "gpu_to_cpu" for x in wf_shims):
        L.append("- Compatibility shim `gpu_to_cpu` in the workflow: rules written for the authors' GPUs (`--cuda`) ran "
                 "on this CPU-only host (slower; results can differ in the last digits)")
    for ps in summary.get("passes") or []:
        L.append(f"- Pass {ps['pass']}: {ps['missing']} missing input file(s)"
                 + (f", {ps['workflow_built']} built by the workflow" if "workflow_built" in ps else "")
                 + f", {ps['now_present']} now present, {len(ps['rerun'])} analyses re-run")
    if summary.get("still_missing"):
        ext = [f_ for f_ in summary["still_missing"] if f_.startswith(("/", ".."))]
        L.append(f"- **Still missing: {len(summary['still_missing'])} file(s)** the analyses read and neither the "
                 f"repository, its workflow nor fetched deposits provide, e.g. `{summary['still_missing'][0]}`")
    for dd in (summary.get("data") or [])[-5:]:
        L.append(f"- Data: `{dd['file']}` from {dd['source']}")
    L += ["", "| # | analysis | ran | chunks (root errors) | printed matched | figures | first root error |",
          "|---|---|---|---|---|---|---|"]
    for k, x in enumerate(summaries, 1):
        ch, prx, fg = x.get("chunks") or {}, x.get("printed") or {}, x.get("figures") or {}
        mt = ((prx.get("identical") or 0) + (prx.get("numerically_close") or 0) + (prx.get("same_values_other_layout") or 0))
        matched = f"{mt}/{prx['reference_blocks']}" if prx.get("reference_blocks") else "-"
        err = (x.get("first_error") or "").replace("|", "/")[:120]
        L.append(f"| {k} | `{x['entry']}` | {x.get('render')} | {ch.get('total', '-')} ({ch.get('root_errors', '-')}) | "
                 f"{matched} | {fg.get('produced', '-')} | {err} |")
    L.append("")
    for k, x in enumerate(summaries, 1):
        sub = ROOT / x["dir"]
        body = (sub / "report.md").read_text() if (sub / "report.md").exists() else ""
        body = "\n".join(("#" + ln if ln.startswith("#") else ln) for ln in body.splitlines()[1:])
        L += [f"## {k}. `{x['entry']}`", "", _rewrite_links(body, os.path.relpath(sub, run_dir) + "/"), ""]
    md = run_dir / "report.md"
    md.write_text("\n".join(L))
    (run_dir / "report.html").write_text(_md_to_html("\n".join(L), f"Reproduction: {src_meta.get('repo')}"))
    _write_json(run_dir / "compare.json", {"entries": [{"entry": x["entry"], "printed": x.get("printed")}
                                                        for x in summaries]})
    res = {"ok": True, "run_dir": _rel(run_dir), "report": _rel(md), "html": _rel(run_dir / "report.html"),
           "summary": _rel(run_dir / "summary.json"), "render": summary["render"], "figures": tot["figures"],
           "printed": tot["printed"], "chunks": tot["chunks"], "entries": len(summaries)}
    _write_json(status, {"state": "done", "at": time.time()})
    _write_json(run_dir / "done.json", res)
    return res


def pipeline(repo: str, *, ref: Optional[str] = None, entry: Optional[str] = None, pins: "Sequence[str]" = (),
             strict: bool = False, paper: Optional[str] = None, run_dir: Optional[Path] = None,
             timeout_min: float = 180, replay: bool = False,
             inputs: "Optional[Dict[str, str]]" = None, entries: "Optional[List[str]]" = None,
             all_entries: bool = False, max_entries: int = 60, entry_timeout_min: float = 60,
             workflow: bool = False, workflow_timeout_min: float = 240, cores: int = 4) -> Dict[str, Any]:
    run_dir = run_dir or _new_run_dir(repo)
    if workflow and not (entries and len(entries) > 1):
        all_entries = all_entries or not entries  # the workflow stage lives in the multi-analysis path
    if entries and len(entries) > 1 or all_entries or workflow:
        sel = list(entries or [])
        if all_entries and not sel:
            meta = fetch(repo, ref)
            if not meta.get("ok"):
                return {"ok": False, "error": meta["error"]}
            top = inventory(ROOT / meta["src"])
            sel = [e for e in top.get("entries") or [] if e.endswith(ANALYSIS_EXT)][:max_entries]
        if len(sel) > 1 or workflow:
            try:
                res = pipeline_multi(repo, sel, ref=ref, pins=pins, strict=strict, paper=paper, run_dir=run_dir,
                                     entry_timeout_min=entry_timeout_min, replay=replay, inputs=inputs,
                                     workflow=workflow, workflow_timeout_min=workflow_timeout_min, cores=cores)
            except Exception as e:  # noqa: BLE001
                res = {"ok": False, "error": f"{type(e).__name__}: {e}"}
            if not res.get("ok"):
                _write_json(run_dir / "status.json", {"state": "failed", "error": res.get("error"), "at": time.time()})
                _write_json(run_dir / "done.json", {**res, "run_dir": _rel(run_dir)})
            return res
        entry = sel[0] if sel else entry
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
    # several analyses of one repository: shared work copy, per-analysis
    # sub-reports, one combined summary; large files hard-linked
    with tempfile.TemporaryDirectory() as td:
        td = str(Path(td).resolve())
        src = Path(td) / "repo"
        (src / "nb").mkdir(parents=True)
        (src / "big.bin").write_bytes(b"x" * (LINK_OVER_BYTES + 10))
        (src / "nb" / "a.py").write_text("print('a', 1+1)\nopen('a_out.txt','w').write('a')\n")
        (src / "nb" / "b.py").write_text("import missing_module_xyz\n")
        fake = Path(td) / "env"
        (fake / "bin").mkdir(parents=True)
        (fake / "bin" / "python").symlink_to(sys.executable)
        g = globals()
        saved = (g["ROOT"], g["fetch"], g["build_env"])
        g["ROOT"] = Path(td)
        g["fetch"] = lambda repo, ref=None: {"ok": True, "repo": repo, "url": "u", "commit": "c0ffee",
                                             "src": "repo"}
        g["build_env"] = lambda inv, pins=(), log=None: {"ok": True, "prefix": str(fake), "versions": "stub"}
        try:
            rd = Path(td) / "run"
            rd.mkdir()
            res = pipeline_multi("x/y", ["nb/a.py", "nb/b.py"], ref=None, pins=(), strict=False, paper="P",
                                 run_dir=rd, entry_timeout_min=2, replay=False, inputs=None)
            sm = _read_json(rd / "summary.json", {})
            ents = sm.get("entries") or []
            check("several analyses: one summary with each analysis and its own report",
                  res.get("ok") and len(ents) == 2 and all((Path(td) / e["dir"] / "report.md").exists() for e in ents))
            check("files an analysis writes are attributed to it", "nb/a_out.txt" in
                  (_read_json(Path(td) / ents[0]["dir"] / "run.json", {}).get("new_files") or []))
            check("large files are hard-linked into the work copy, not copied",
                  os.stat(rd / "work" / "big.bin").st_ino == os.stat(src / "big.bin").st_ino)
            b_sum = _read_json(Path(td) / ents[1]["dir"] / "summary.json", {})
            check("an analysis that fails to run is a root error, never 'clean'",
                  (b_sum.get("chunks") or {}).get("root_errors") == 1 and "did not execute" in
                  (b_sum.get("errors") or [{}])[0].get("message", "") and "1/2 ran without root errors" in
                  sm.get("render", ""))
            (src / "nb" / "r_reader.py").write_text("print(open('shared.txt').read())\n")
            (src / "nb" / "w_writer.py").write_text("open('shared.txt','w').write('made by writer')\n")
            rd2 = Path(td) / "run2"
            rd2.mkdir()
            pipeline_multi("x/y", ["nb/r_reader.py", "nb/w_writer.py"], ref=None, pins=(), strict=False, paper="P",
                           run_dir=rd2, entry_timeout_min=2, replay=False, inputs=None)
            sm2 = _read_json(rd2 / "summary.json", {})
            reader = next(e for e in sm2["entries"] if e["entry"] == "nb/r_reader.py")
            check("an analysis whose input another analysis writes is re-run once it exists",
                  reader.get("passes") == 2 and not (reader.get("chunks") or {}).get("root_errors")
                  and sm2["passes"][0]["rerun"] == ["nb/r_reader.py"])
            check("the combined report lists every analysis", "`nb/a.py`" in (rd / "report.md").read_text()
                  and "`nb/b.py`" in (rd / "report.md").read_text())
        finally:
            g["ROOT"], g["fetch"], g["build_env"] = saved
    with tempfile.TemporaryDirectory() as td:
        td = str(Path(td).resolve())
        repo_src = Path(td) / "PaperCode" / "o__r" / "src"
        (repo_src / "workflow").mkdir(parents=True)
        (repo_src / "workflow" / "Snakefile").write_text("rule all:\n    input: 'results/x.txt'\n")
        check("a Snakemake workflow is found", find_workflow(repo_src) == {"engine": "snakemake",
                                                                          "file": "workflow/Snakefile"})
        txt = ("Building DAG of jobs...\nMissingInputException in rule map in file workflow/Snakefile, line 3:\n"
               "Missing input files for rule map:\n    output: results/mapped/a.h5ad\n    affected files:\n"
               "        results/raw/A_R1.fastq.gz\n        results/raw/A_R2.fastq.gz\n")
        mi = _missing_inputs(txt)
        check("missing workflow inputs are listed", "results/raw/A_R1.fastq.gz" in mi and "results/raw/A_R2.fastq.gz" in mi)
        sub = Path(td) / "sub"
        sub.mkdir()
        _write_json(sub / "summary.json", {"errors": [
            {"message": "FileNotFoundError: [Errno 2] No such file or directory: '../../results/model_runs/x.csv'"},
            {"message": "OSError: Unable to open file (unable to open file: name = '../../results/mapped/y.h5ad', errno 2)"},
            {"message": "No such file or directory: '/etc/outside'"}]})
        mf = missing_files(sub, "workflow/notebooks/Fig1/a.ipynb", repo_src)
        check("files an analysis could not open are resolved against its own directory",
              mf == ["workflow/results/model_runs/x.csv", "workflow/results/mapped/y.h5ad"])
        payload = Path(td) / "deposit.tsv"
        payload.write_text("a\tb\n1\t2\n")
        g = globals()
        saved = g["CODE_DIR"], g["ROOT"]
        g["CODE_DIR"], g["ROOT"] = Path(td) / "PaperCode", Path(td)
        try:
            res = fetch_data("o/r", url=payload.as_uri(), into="resources/data")
            prov = _read_json(Path(td) / "PaperCode" / "o__r" / "data.json", [])
            check("deposited data is fetched into the checkout with provenance",
                  res["ok"] and (repo_src / "resources" / "data" / "deposit.tsv").read_text().startswith("a\tb")
                  and prov and prov[0]["source"].startswith("file:"))
            check("data cannot be written outside the repository", not fetch_data("o/r", url=payload.as_uri(),
                                                                                    into="../../x")["ok"])
        finally:
            g["CODE_DIR"], g["ROOT"] = saved
    runs = [{"run_accession": "SRR1", "library_name": "LDLvar_rep1_top", "fastq_bytes": "10;12",
             "fastq_ftp": "ftp.x/SRR1_1.fastq.gz;ftp.x/SRR1_2.fastq.gz", "fastq_md5": "a;b"},
            {"run_accession": "SRR2", "library_name": "LDLRCDS_CBE_SpRY_rep2", "fastq_bytes": "5",
             "fastq_ftp": "ftp.x/SRR2_1.fastq.gz", "fastq_md5": "c"},
            {"run_accession": "SRR3", "library_name": "20loci_ATAC_rep1", "fastq_bytes": "1",
             "fastq_ftp": "ftp.x/SRR3_1.fastq.gz", "fastq_md5": "d"}]
    pl = plan_reads(runs, "results/raw/{lib}/{library_name}_R{read}.fastq.gz",
                    r"^(?P<lib>LDLvar|LDLRCDS(?:_CBE_SpRY)?)_(?:rep|plasmid)")
    got = [(x["run"], x.get("dest"), x.get("skip")) for x in pl]
    check("runs map onto the authors' FASTQ names (paired reads, regex groups); non-matching runs are skipped",
          ("SRR1", "results/raw/LDLvar/LDLvar_rep1_top_R2.fastq.gz", None) in got
          and ("SRR2", "results/raw/LDLRCDS_CBE_SpRY/LDLRCDS_CBE_SpRY_rep2_R1.fastq.gz", None) in got
          and any(r == "SRR3" and sk for r, _, sk in got))
    gy = GPU_ONLY
    check("GPU-only conda packages are recognised", all(gy.match(n) for n in ("cudatoolkit", "cudnn", "pytorch-mutex"))
          and not gy.match("pytorch") and not gy.match("numpy"))
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
        s.add_argument("--entry", action="append", default=None,
                       help="analysis file inside the repository (repeatable; default: the main Rmd/notebook)")
        s.add_argument("--all-entries", action="store_true",
                       help="run every notebook / R Markdown file in the repository, in path order")
        s.add_argument("--max-entries", type=int, default=60)
        s.add_argument("--entry-timeout-min", type=float, default=60, help="per analysis, with several")
        s.add_argument("--workflow", action="store_true",
                       help="run the authors' Snakemake workflow first (dry run + real run), then the analyses")
        s.add_argument("--workflow-timeout-min", type=float, default=240)
        s.add_argument("--cores", type=int, default=int(os.environ.get("IGVF_PAPER_CODE_CORES", "4")))
        s.add_argument("--pin", action="append", default=[], help="conda/pip pin, e.g. r-reshape2=1.4.4")
        s.add_argument("--strict", action="store_true", help="stop at the first failing chunk")
        s.add_argument("--paper", help="paper title/DOI for the report")
        s.add_argument("--timeout-min", type=float, default=180)
        s.add_argument("--replay", action="store_true", help="execute twice and check the outputs are identical")
        s.add_argument("--input", action="append", default=[], metavar="NAME=PATH",
                       help="run the same code on new data: replace input NAME (as the analysis reads it) with PATH")
        s.add_argument("--detach", action="store_true", help="run in the background; job_wait on <run>/done.json")
        s.add_argument("--run-dir", help=argparse.SUPPRESS)
    s = sub.add_parser("data", help="Fetch a paper's deposited data (Zenodo record/DOI or URL) into its checkout")
    s.add_argument("repo")
    s.add_argument("--zenodo", help="record id, zenodo.org URL or 10.5281/zenodo.N DOI")
    s.add_argument("--url")
    s.add_argument("--into", default="", help="directory inside the repository (where the analyses read it)")
    s.add_argument("--files", help="glob of record files to take (default: all)")
    s.add_argument("--extract", action="store_true", help="unpack .zip/.tar.gz archives after download")
    s = sub.add_parser("reads", help="Raw reads of a study (SRA/ENA/GEO) into the paths the authors' code reads")
    s.add_argument("repo")
    s.add_argument("--accession", help="BioProject (PRJNA/PRJEB), SRA study (SRP/ERP), GEO series (GSE) or run")
    s.add_argument("--layout", help="destination template relative to --into, with run-table fields and --match "
                                   "groups, e.g. 'results/raw/{lib}/{library_name}_R{read}.fastq.gz'")
    s.add_argument("--into", default="", help="directory inside the repository the layout is relative to")
    s.add_argument("--match", help="regex with named groups applied to --field, e.g. '^(?P<lib>.+?)_(?:rep|plasmid)'")
    s.add_argument("--field", default="library_name")
    s.add_argument("--include", help="only runs whose --field matches this regex")
    s.add_argument("--guess-layout", action="store_true", help="show the FASTQ path patterns in the authors' code")
    s.add_argument("--yes", action="store_true", help="download (default: preview the mapping only)")
    s.add_argument("--max-gb", type=float, default=200)
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
        res = pipeline(repo, ref=args.ref, entry=(args.entry or [None])[0], pins=args.pin, strict=args.strict,
                       paper=args.paper, entries=args.entry, all_entries=args.all_entries,
                       max_entries=args.max_entries, entry_timeout_min=args.entry_timeout_min,
                       workflow=args.workflow, workflow_timeout_min=args.workflow_timeout_min, cores=args.cores,
                       run_dir=Path(args.run_dir) if args.run_dir else None, timeout_min=args.timeout_min,
                       replay=args.replay,
                       inputs=dict(x.split("=", 1) for x in args.input) or None)
        print(json.dumps({k: v for k, v in res.items()}, indent=2, default=str))
        if res.get("ok"):
            _announce(res)
        return 0 if res.get("ok") else 1
    if args.cmd == "reads":
        repo = find(repo=args.repo)["candidates"][0]["repo"]
        if args.guess_layout:
            for ln in guess_read_layout(repo_dir(repo) / "src"):
                print(ln)
            return 0
        if not (args.accession and args.layout):
            print("give --accession and --layout (see --guess-layout)")
            return 2
        res = fetch_reads(repo, args.accession, args.layout, into=args.into, match=args.match, field=args.field,
                          include=args.include, download=args.yes, max_gb=args.max_gb)
        print(json.dumps(res, indent=2, default=str))
        return 0 if res.get("ok") else 1
    if args.cmd == "data":
        res = fetch_data(find(repo=args.repo)["candidates"][0]["repo"], zenodo=args.zenodo, url=args.url,
                         into=args.into, files=args.files, extract=args.extract)
        print(json.dumps(res, indent=2))
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
