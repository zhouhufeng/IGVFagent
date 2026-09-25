"""Reproduction records: one durable, searchable record per paper.

Every finished reproduction or benchmark job leaves its evidence scattered over
the workspace: the job's plan and events (Data/Jobs/<id>/), the paper harvest
(Docs/Benchmark/<run>/), the authors'-code run (Docs/PaperCode/<run>/) and the
concordance scoring (Benchmarks/results/). This module gathers them into one
record per paper, in the persistent data volume, with every attempt kept:

    Data/Reproductions/<paper_id>/record.json          latest attempt + history
    Data/Reproductions/<paper_id>/attempts/<job>.json  each attempt
    Data/Reproductions/<paper_id>/<job>.html           self-contained final report
                                                       (figures embedded), like
                                                       Paper2Agent's delivered HTML

The web UI lists and searches the records (📊 Validation → 📑 Reproductions).
A record is private to its owner until it is published. The owner can post a
summary (headline numbers, key figures, a selected log and a link back) to the
IGVF Agent category on discussion.genohub.org; nothing is posted without an
explicit confirmation, and the post is authored as the owner's forum account.

    igvfagent repro record <job_id>          build/refresh the record of a job
    igvfagent repro backfill                 records for every finished job
    igvfagent repro list [--query TEXT]      search records
    igvfagent repro show <paper_id>
    igvfagent repro publish <paper_id> [--off]
    igvfagent repro post <paper_id> [--yes]  preview (default) or post to the forum
    igvfagent repro selftest
"""
from __future__ import annotations

import argparse
import base64
import html as _html
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
REPRO_DIR = Path(os.environ.get("IGVF_REPRO_DIR") or ROOT / "Data" / "Reproductions")
JOBS_DIR = Path(os.environ.get("IGVF_JOBS_DIR") or ROOT / "Data" / "Jobs")
FORUM_URL = os.environ.get("IGVF_DISCOURSE_URL", "https://discussion.genohub.org").rstrip("/")
FORUM_CATEGORY = int(os.environ.get("IGVF_REPRO_FORUM_CATEGORY", "30"))  # IGVF Agent › IGVF Agent Help
PUBLIC_URL = os.environ.get("IGVF_PUBLIC_URL", "https://igvfagent.genohub.org").rstrip("/")

REPRO_QUERY = re.compile(r"\b(reproduc\w*|replicat\w*|benchmark\w*|re-?run the (?:paper|analysis))\b", re.I)
_LOG_KINDS = {"created", "round", "stage", "verify", "done", "stopped", "failed", "resumed", "budget"}


def _read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return default


def _write_json(p: Path, obj: Any) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, p)


def _abs(rel: str) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else ROOT / p


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", (s or "").lower()).strip("_")[:60] or "paper"


# ─── gather a job's evidence ────────────────────────────────────────────────

def _evidence(job: dict, plan: dict) -> "List[str]":
    ev = [e for st in plan.get("stages") or [] for e in (st.get("evidence") or [])]
    ev += [a for a in job.get("artefacts") or [] if isinstance(a, str)]
    # paths the final report names (runs referenced but not attached as evidence)
    ans = job.get("last_answer") or ""
    ev += re.findall(r"((?:Docs|Benchmarks|Data)/[\w./\-]+)", ans)
    seen, out = set(), []
    for e in ev:
        e = e.rstrip(".,);`'\"")
        if e not in seen:
            seen.add(e)
            out.append(e)
    return out


def _run_dirs(ev: "List[str]", marker: str) -> "List[Path]":
    dirs = []
    for e in ev:
        m = re.match(rf"({marker}/[^/]+)", e)
        if m:
            d = _abs(m.group(1))
            if d.is_dir() and d not in dirs:
                dirs.append(d)
    return dirs


def paper_info(ev: "List[str]", job: dict) -> Dict[str, Any]:
    for d in reversed(_run_dirs(ev, "Docs/Benchmark")):
        res = _read_json(d / "resolution.json") or {}
        hv = _read_json(d / "harvest.json") or {}
        p = res.get("paper") or {}
        if p or hv:
            claims = [{"value": c.get("value"), "unit": c.get("unit"), "section": c.get("section"),
                       "quote": (c.get("quote") or "")[:300], "source": "text"}
                      for c in (hv.get("numeric_claims") or [])[:20]]
            claims += [{"value": c.get("value"), "unit": c.get("unit"), "description": c.get("description"),
                        "quote": (c.get("quote") or "")[:300], "source": "llm",
                        "grounded": c.get("quote_grounded_in_source")}
                       for c in ((hv.get("llm_claims") or {}).get("claims") or [])[:20]]
            code = [h.get("value") for h in (hv.get("accessions") or {}).get("github_repo") or []
                    if h.get("in_data_availability")]
            return {"paper_id": res.get("paper_id") or hv.get("paper_id") or _slug(p.get("title") or ""),
                    "title": p.get("title") or hv.get("title") or job.get("title"),
                    "doi": p.get("doi") or hv.get("doi"), "journal": p.get("journal") or p.get("container"),
                    "year": p.get("year"), "authors": p.get("authors") or p.get("first_author"),
                    "harvest": _rel(d), "claims": claims, "code_repositories": code,
                    "genes": (hv.get("genes") or [])[:12], "assays": [a.get("assay") for a in hv.get("assays") or []][:8]}
    # No harvest: name the paper after a benchmark or the job title.
    bench = next((re.match(r"Benchmarks/([^/]+)/", e).group(1) for e in ev
                  if re.match(r"Benchmarks/([^/]+)/", e) and not e.startswith("Benchmarks/results")), None)
    doi = re.search(r"\b10\.\d{4,9}/[^\s\"'<>)]+", job.get("query") or "")
    return {"paper_id": bench or _slug(job.get("title") or job["id"]), "title": job.get("title"),
            "doi": doi.group(0).rstrip(".") if doi else None, "claims": [], "code_repositories": []}


def _entry_rows(d: Path, s: Dict[str, Any]) -> "List[Dict[str, Any]]":
    if s.get("entries"):
        return [dict(e, run_dir=_rel(d)) for e in s["entries"]]
    return [{"entry": s.get("entry"), "render": s.get("render"), "chunks": s.get("chunks"), "figures": s.get("figures"),
             "printed": s.get("printed"), "replay": s.get("replay"), "run_dir": _rel(d)}]


def code_info(ev: "List[str]") -> Optional[Dict[str, Any]]:
    """The authors'-code result of a job. A job often runs several analyses of
    one repository (one notebook per figure), in one or several paper-code
    runs: every analysis counts, with its latest run."""
    dirs = sorted({d for d in _run_dirs(ev, "Docs/PaperCode") if (d / "summary.json").exists()}, key=lambda d: d.name)
    if not dirs:
        return None
    rows: Dict[str, Dict[str, Any]] = {}
    for d_ in dirs:
        for r in _entry_rows(d_, _read_json(d_ / "summary.json") or {}):
            rows[str(r.get("entry"))] = r
    d = next((x for x in reversed(dirs) if (_read_json(x / "summary.json") or {}).get("entries")), dirs[-1])
    s = _read_json(d / "summary.json") or {}
    run = _read_json(d / "run.json") or {}
    entries = list(rows.values())
    if len(entries) > 1:
        tot = {"chunks": {"total": 0, "errored": 0, "root_errors": 0, "cascade_errors": 0},
               "figures": {"produced": 0, "reference": 0}, "printed": {"reference_blocks": 0, "identical": 0,
                                                                     "numerically_close": 0,
                                                                     "same_values_other_layout": 0, "missing": 0}}
        for e in entries:
            for sect, keys in tot.items():
                for k in keys:
                    keys[k] += int((e.get(sect) or {}).get(k) or 0)
        pr = tot["printed"]
        got = pr["identical"] + pr["numerically_close"] + pr["same_values_other_layout"]
        pr["fraction_matched"] = round(got / pr["reference_blocks"], 4) if pr["reference_blocks"] else None
        s = {**s, **tot, "entry": f"{len(entries)} analyses"}
    claims = []
    for x in reversed(dirs):  # the latest run with a claims ledger
        cs = (_read_json(x / "claims_summary.json") or {}).get("claims")
        if cs:
            claims = cs
            d = x
            break
    return {"claims": claims, "entries": entries, "runs": [_rel(x) for x in dirs],
            "run_dir": _rel(d), "repo": s.get("repo"), "commit": s.get("commit"), "entry": s.get("entry"),
            "language": s.get("language"), "render": s.get("render"), "seconds": s.get("seconds"),
            "chunks": s.get("chunks"), "figures": s.get("figures"), "printed": s.get("printed"),
            "replay": s.get("replay"), "shims": run.get("r_shims") or [],
            "rng_sample_rounding": run.get("rng_sample_rounding"), "sections": s.get("sections"),
            "report_md": _rel(d / "report.md"), "report_html": _rel(d / "report.html")}


def data_info(ev: "List[str]") -> "List[Dict[str, Any]]":
    out = []
    for e in ev:
        if e.startswith("Benchmarks/results/") and e.endswith("_concordance.json"):
            for b in _read_json(_abs(e), []) or []:
                out.append({"benchmark": b.get("paper"), "skill": b.get("skill"), "status": b.get("status"),
                            "n_passed": b.get("n_passed"), "n_total": b.get("n_total"), "run_dir": b.get("run_dir"),
                            "checks": [{"name": c.get("name"), "passed": c.get("passed"), "detail": c.get("detail")}
                                       for c in b.get("checks") or []], "results": e})
    return out


def selected_log(job_id: str, limit: int = 60) -> "List[Dict[str, Any]]":
    """The events worth reading later: rounds, stage verdicts, verifier
    verdicts, the end. Tool-by-tool chatter stays in events.jsonl."""
    out = []
    try:
        lines = (JOBS_DIR / job_id / "events.jsonl").read_text().splitlines()
    except OSError:
        return out
    for ln in lines:
        try:
            e = json.loads(ln)
        except ValueError:
            continue
        if e.get("kind") in _LOG_KINDS:
            out.append({"t": e.get("t"), "kind": e.get("kind"), "text": (e.get("text") or "")[:300],
                        **({"issues": e["issues"][:6]} if e.get("issues") else {})})
    return out[-limit:]


def outcome(rec: Dict[str, Any]) -> str:
    """reproduced | partial | not reproduced | running — from harness facts only."""
    j = rec["job"]
    if j["status"] in ("running", "queued"):
        return "running"
    verdict_ok = (j.get("verdict") or {}).get("verdict") == "pass"
    code, data = rec.get("code"), rec.get("data") or []
    code_ok = bool(code and (code.get("printed") or {}).get("fraction_matched") is not None
                   and code["printed"]["fraction_matched"] >= 0.9 and not (code.get("chunks") or {}).get("root_errors"))
    code_ok = code_ok or bool(code and not (code.get("printed") or {}) and (code.get("figures") or {}).get("produced")
                              and not (code.get("chunks") or {}).get("root_errors"))
    data_ok = bool(data) and all(b.get("status") == "ok" for b in data)
    if code and code.get("claims"):
        vs = [x.get("verdict") for x in code["claims"]]
        if vs and all(v == "reproduced" for v in vs):
            return "reproduced" if (verdict_ok or j.get("standalone")) else "partial"
        return "partial" if any(v in ("reproduced", "partially reproduced") for v in vs) else "not reproduced"
    code_named = bool((rec.get("paper") or {}).get("code_repositories"))
    if code_named and not code:  # the paper's own code was never run: a cross-check, not the reproduction
        return "partial" if (data_ok or j["status"] in ("done", "done_with_blocked")) else "not reproduced"
    if j["status"] == "done" and (verdict_ok or j.get("standalone")) and (code_ok or data_ok):
        return "reproduced"
    if code_ok or data_ok or j["status"] in ("done", "done_with_blocked"):
        return "partial"
    return "not reproduced"


def record_from_job(job_id: str) -> Optional[Dict[str, Any]]:
    job = _read_json(JOBS_DIR / job_id / "job.json")
    if not job:
        return None
    plan = _read_json(JOBS_DIR / job_id / "plan.json", {"stages": []})
    ev = _evidence(job, plan)
    paper = paper_info(ev, job)
    rec: Dict[str, Any] = {
        "paper_id": paper["paper_id"], "paper": {k: v for k, v in paper.items() if k != "paper_id"},
        "job": {"id": job["id"], "title": job.get("title"), "query": (job.get("query") or "")[:2000],
                "owner": job.get("owner") or "", "status": job.get("status"), "orchestrator": job.get("orchestrator"),
                "model": job.get("model"), "rounds": job.get("rounds"), "cost_usd": job.get("cost_usd"),
                "created_at": job.get("created_at"), "finished_at": job.get("finished_at"),
                "verdict": job.get("last_verdict") or {}},
        "stages": [{"id": s["id"], "title": s.get("title"), "status": s["status"],
                    "reason": (s.get("verified") or {}).get("reason"), "notes": (s.get("notes") or [])[-2:],
                    "evidence": (s.get("evidence") or [])[:8]} for s in plan.get("stages") or []],
        "code": code_info(ev), "data": data_info(ev), "log": selected_log(job_id),
        "answer": (JOBS_DIR / job_id / "answer.md").read_text()[:60000]
        if (JOBS_DIR / job_id / "answer.md").exists() else (job.get("last_answer") or "")[:60000],
        "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    rec["route"] = ("authors_code+public_data" if rec["code"] and rec["data"] else
                    "authors_code" if rec["code"] else "public_data" if rec["data"] else "other")
    rec["outcome"] = outcome(rec)
    return _store(rec)


def _dir_for(rec: Dict[str, Any]) -> Path:
    """The same paper keeps one record even when attempts name it differently
    (a benchmark id, a harvest id): match on DOI, then on the code repository."""
    doi = (rec["paper"].get("doi") or "").lower()
    repo = ((rec.get("code") or {}).get("repo") or "").lower()
    if REPRO_DIR.is_dir():
        for d in REPRO_DIR.iterdir():
            r = _read_json(d / "record.json") or {}
            if doi and (r.get("paper") or {}).get("doi", "") and r["paper"]["doi"].lower() == doi:
                return d
            if repo and not doi and ((r.get("code") or {}).get("repo") or "").lower() == repo:
                return d
    return REPRO_DIR / _slug(rec["paper_id"])


def _store(rec: Dict[str, Any]) -> Dict[str, Any]:
    job_id = rec["job"]["id"]
    d = _dir_for(rec)
    rec["paper_id"] = d.name
    html_path = d / f"{job_id}.html"
    rec["html"] = _rel(html_path)
    _write_json(d / "attempts" / f"{job_id}.json", rec)
    html_path.write_text(render_html(rec))
    _refresh_record(d)
    return rec


def _refresh_record(d: Path) -> Dict[str, Any]:
    old = _read_json(d / "record.json", {}) or {}
    attempts = sorted((_read_json(p) for p in (d / "attempts").glob("*.json")),
                      key=lambda r: (r or {}).get("job", {}).get("created_at") or "")
    attempts = [a for a in attempts if a]
    if not attempts:
        return old
    # What any attempt learned about the paper applies to all: once the paper is
    # known to publish its code, a data-only attempt is a cross-check (partial).
    named = sorted({r for a in attempts for r in (a["paper"].get("code_repositories") or []) if r}
                   | {(a.get("code") or {}).get("repo") for a in attempts if a.get("code")} - {None})
    doi = next((a["paper"].get("doi") for a in attempts if a["paper"].get("doi")), None)
    for a in attempts:
        a["paper"]["code_repositories"] = a["paper"].get("code_repositories") or named
        a["paper"]["doi"] = a["paper"].get("doi") or doi
        a["outcome"] = outcome(a)
    rank = {"reproduced": 3, "partial": 2, "running": 1, "not reproduced": 0}
    latest = attempts[-1]
    best = max(attempts, key=lambda a: (rank.get(a["outcome"], 0), a["job"].get("created_at") or ""))
    rec = dict(best)
    rec["latest_job"] = latest["job"]["id"]
    rec["attempts"] = [{"job": a["job"]["id"], "created_at": a["job"].get("created_at"), "outcome": a["outcome"],
                        "route": a["route"], "status": a["job"]["status"], "owner": a["job"].get("owner"),
                        "html": a.get("html")} for a in attempts]
    rec["owners"] = sorted({a["job"].get("owner") or "" for a in attempts})
    rec["published"] = bool(old.get("published"))
    rec["forum"] = old.get("forum") or {}
    spec = tool_spec(rec)
    if spec:
        _write_json(d / "tool.json", spec)
        rec["tool"] = spec["name"]
    elif (d / "tool.json").exists():
        (d / "tool.json").unlink()
    _write_json(d / "record.json", rec)
    return rec


# ─── per-paper tools: the authors' analysis, reusable on new data ───────────
#
# Paper2MCP's end product is the reproduced analysis exposed as tools. Here a
# paper whose authors' code was reproduced becomes one tool, paper_<id>, that
# re-runs that pinned code (same commit, same environment repairs and shims)
# with the user's files in place of the inputs it reads, and reports the
# result section by section. The method stays the authors'; nothing is
# re-implemented.

def _param(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", Path(name).name.lower()).strip("_")[:40] or "input"


def tool_spec(rec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    c = rec.get("code")
    if rec.get("outcome") != "reproduced" or not c or not c.get("run_dir"):
        return None
    rd = _abs(c["run_dir"])
    inv = _read_json(rd / "inventory.json") or {}
    env = _read_json(rd / "env.json") or {}
    if not inv.get("ok"):
        return None
    inputs = {}
    for i in inv.get("inputs") or []:
        if i.get("exists") and not i.get("is_url"):
            k = _param(i["path"])
            while k in inputs or k in ("replay", "label", "detach"):
                k += "_"
            inputs[k] = i["path"]
    pins = [x for x in env.get("spec") or [] if x.startswith("cran:")]
    title = rec["paper"].get("title") or rec["paper_id"]
    secs = [s_["title"] for s_ in (inv.get("sections") or [])][:12]
    name = ("paper_" + re.sub(r"[^a-z0-9_]", "_", rec["paper_id"].lower()))[:64]
    props: Dict[str, Any] = {k: {"type": "string", "description": f"Your file in place of `{v}` (the authors' input "
                                                                   f"of the same format). Omit to keep theirs."}
                             for k, v in inputs.items()}
    props["replay"] = {"type": "boolean", "description": "Execute twice and check the outputs are identical."}
    desc = (f"Run the authors' own analysis from \"{title[:140]}\" ({c.get('repo')} @ {(c.get('commit') or '')[:10]}, "
            f"entry {c.get('entry')}) on NEW DATA. Reproduced by IGVF Agent: {headline(rec)[:160]}. Pass your files in "
            f"place of the inputs it reads ({', '.join(inputs.values())[:200]}); the same code, environment and "
            f"compatibility shims run, and the report gives every section's figures and printed values"
            + (f" (sections: {'; '.join(secs)[:300]})" if secs else "")
            + ". Runs in the background: wait with job_wait on <run_dir>/done.json or paper_code_status.")
    return {"name": name, "description": desc, "paper_id": rec["paper_id"], "repo": c.get("repo"),
            "commit": c.get("commit"), "entry": c.get("entry"), "pins": pins, "inputs": inputs,
            "parameters": {"type": "object", "properties": props},
            "cli": ["repro", "apply", rec["paper_id"], "--detach"],
            "flag_map": {k: "--" + k.replace("_", "-") for k in inputs}, "bool_flags": ["replay"]}


def paper_tools() -> "List[Dict[str, Any]]":
    out = []
    if REPRO_DIR.is_dir():
        for d in sorted(REPRO_DIR.iterdir()):
            t = _read_json(d / "tool.json")
            if t:
                out.append(t)
    return out


def apply(paper_id: str, files: Dict[str, str], *, replay: bool = False, detach: bool = False,
          label: str = "") -> Dict[str, Any]:
    spec = _read_json(REPRO_DIR / _slug(paper_id) / "tool.json")
    if not spec:
        return {"ok": False, "error": f"{paper_id} has no reusable analysis (only papers whose authors' code was "
                                      "reproduced do)"}
    bad = [k for k in files if k not in spec["inputs"]]
    if bad:
        return {"ok": False, "error": f"unknown inputs {bad}; this analysis reads {sorted(spec['inputs'])}"}
    for k, v in files.items():
        if not _abs(v).is_file():
            return {"ok": False, "error": f"{k}: file not found: {v}"}
    argv = ["pipeline", spec["repo"], "--ref", spec["commit"], "--entry", spec["entry"],
            "--paper", f"{spec['paper_id']} — the authors' analysis applied to new data"
                       + (f" ({label})" if label else "")]
    for pin in spec.get("pins") or []:
        argv += ["--pin", pin]
    for k, v in files.items():
        argv += ["--input", f"{spec['inputs'][k]}={_abs(v)}"]
    if replay:
        argv.append("--replay")
    if detach:
        argv.append("--detach")
    try:
        from igvfagent import paper_code_skill as pc  # type: ignore
    except Exception:
        import paper_code_skill as pc  # type: ignore
    rc = pc.main(argv)
    return {"ok": rc == 0, "paper_id": spec["paper_id"], "argv": argv}


def record_from_run(run_dir: str, *, paper_id: Optional[str] = None, doi: Optional[str] = None,
                    title: Optional[str] = None, owner: str = "", note: str = "") -> Optional[Dict[str, Any]]:
    """A record for a standalone `paper-code pipeline` run (no job)."""
    d = _abs(run_dir)
    s = _read_json(d / "summary.json")
    if not s:
        return None
    code = code_info([_rel(d) + "/summary.json"])
    t = title or s.get("paper") or s.get("repo")
    ents = s.get("entries") or []
    clean = sum(1 for e in ents if (e.get("chunks") or {}).get("total") and not (e.get("chunks") or {}).get("root_errors"))
    # a multi-analysis run counts as done when at least one analysis ran clean
    st_ = "done" if s.get("render") == "ok" or clean else "failed"
    rec: Dict[str, Any] = {
        "paper_id": paper_id or _slug(t), "paper": {"title": t, "doi": doi, "claims": [],
                                                    "code_repositories": [s.get("repo")]},
        "job": {"id": "run-" + d.name, "title": f"paper-code {s.get('repo')}", "owner": owner, "status": st_,
                "orchestrator": "paper-code (standalone)", "rounds": None, "standalone": True,
                "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(d.stat().st_mtime)),
                "finished_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime((d / "summary.json").stat().st_mtime)),
                "verdict": {}},
        "stages": [], "code": code, "data": [], "log": [],
        "answer": note, "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "route": "authors_code"}
    rec["outcome"] = outcome(rec)
    return _store(rec)


def backfill() -> "List[str]":
    done = []
    if not JOBS_DIR.is_dir():
        return done
    for jd in sorted(JOBS_DIR.iterdir()):
        j = _read_json(jd / "job.json") or {}
        if j.get("status") in ("running", "queued") or not j:
            continue
        if REPRO_QUERY.search(j.get("query") or "") or any(
                k in json.dumps(_read_json(jd / "plan.json", {})) for k in ("PaperCode", "concordance", "harvest.json")):
            r = record_from_job(j["id"])
            if r:
                done.append(f"{j['id']} -> {r['paper_id']} ({r['outcome']})")
    return done


# ─── list / search / visibility ─────────────────────────────────────────────

def visible(rec: Dict[str, Any], viewer: Optional[str], is_admin: bool) -> bool:
    return viewer is None or is_admin or rec.get("published") or viewer in (rec.get("owners") or [])


def list_records(viewer: Optional[str] = None, is_admin: bool = True, query: str = "") -> "List[Dict[str, Any]]":
    out = []
    if not REPRO_DIR.is_dir():
        return out
    toks = [t for t in re.split(r"\s+", query.lower().strip()) if t]
    for d in REPRO_DIR.iterdir():
        r = _read_json(d / "record.json")
        if not r or not visible(r, viewer, is_admin):
            continue
        hay = " ".join(str(x) for x in (r["paper_id"], r["paper"].get("title"), r["paper"].get("doi"),
                                         r["paper"].get("journal"), r["paper"].get("year"),
                                         (r.get("code") or {}).get("repo"), " ".join(r["paper"].get("genes") or []),
                                         " ".join(a for a in r["paper"].get("assays") or [] if a),
                                         r["outcome"], r["route"], r["job"].get("owner"))).lower()
        if all(t in hay for t in toks):
            out.append(r)
    out.sort(key=lambda r: r["job"].get("finished_at") or r["job"].get("created_at") or "", reverse=True)
    return out


def load(paper_id: str) -> Optional[Dict[str, Any]]:
    return _read_json(REPRO_DIR / _slug(paper_id) / "record.json")


def set_published(paper_id: str, on: bool) -> Dict[str, Any]:
    d = REPRO_DIR / _slug(paper_id)
    r = _read_json(d / "record.json")
    if not r:
        return {"ok": False, "error": f"no record {paper_id}"}
    r["published"] = bool(on)
    _write_json(d / "record.json", r)
    return {"ok": True, "published": r["published"]}


def headline(rec: Dict[str, Any]) -> str:
    c, parts = rec.get("code"), []
    if c and c.get("claims"):
        cl = c["claims"]
        n_ok = sum(1 for x in cl if x.get("verdict") == "reproduced")
        n_part = sum(1 for x in cl if x.get("verdict") == "partially reproduced")
        n_na = sum(1 for x in cl if x.get("verdict") == "not attempted")
        return (f"{n_ok} of {len(cl)} of the paper's results reproduced" + (f", {n_part} partially" if n_part else "")
                + (f", {n_na} not attempted" if n_na else "")
                + ": " + "; ".join(f"{x['title']} ({x['verdict']}, r = {x['primary_pearson']:.3f})"
                                   for x in cl if x.get("primary_pearson") is not None))[:400]
    if c and len(c.get("entries") or []) > 1:
        ents = c["entries"]
        clean = sum(1 for e in ents if (e.get("chunks") or {}).get("total") and not (e.get("chunks") or {}).get("root_errors"))
        parts.append(f"{clean}/{len(ents)} analyses ran without root errors")
    if c:
        pr = c.get("printed") or {}
        if pr.get("reference_blocks"):
            got = (pr.get("identical") or 0) + (pr.get("numerically_close") or 0) + (pr.get("same_values_other_layout") or 0)
            parts.append(f"{got}/{pr['reference_blocks']} printed values match the authors'")
        ch = c.get("chunks") or {}
        if ch.get("total"):
            parts.append(f"{ch['total'] - (ch.get('errored') or 0)}/{ch['total']} chunks ran")
        if (c.get("figures") or {}).get("produced"):
            parts.append(f"{c['figures']['produced']} figures")
        if (c.get("replay") or {}).get("deterministic"):
            parts.append("replay deterministic")
    for b in rec.get("data") or []:
        parts.append(f"{b['benchmark']}: {b['n_passed']}/{b['n_total']} checks")
    return "; ".join(parts) or (rec["job"].get("verdict") or {}).get("verdict") or rec["job"]["status"]


# ─── the self-contained HTML report ─────────────────────────────────────────

_ICON = {"reproduced": "✅", "partial": "🟡", "not reproduced": "❌", "running": "🔄"}


def _md_html(md: str, base: Optional[Path] = None, max_img: int = 1_500_000) -> str:
    """Small Markdown -> HTML with images inlined as data URIs, so the page is
    one portable file."""
    out, in_code, in_list = [], False, False
    for ln in md.splitlines():
        if ln.startswith("```"):
            out.append("</pre>" if in_code else "<pre>")
            in_code = not in_code
            continue
        if in_code:
            out.append(_html.escape(ln))
            continue
        e = _html.escape(ln)

        def img(m: "re.Match") -> str:
            alt, src = m.group(1), urllib.parse.unquote(_html.unescape(m.group(2)))
            p = (base / src) if base and not src.startswith(("http", "data:")) else None
            if p and p.is_file() and p.stat().st_size <= max_img and p.suffix.lower() in (".png", ".jpg", ".jpeg", ".svg"):
                mime = "image/svg+xml" if p.suffix.lower() == ".svg" else f"image/{p.suffix.lower().lstrip('.').replace('jpg', 'jpeg')}"
                src = f"data:{mime};base64," + base64.b64encode(p.read_bytes()).decode()
            return f'<figure><img src="{src}" alt="{alt}" loading="lazy"><figcaption>{alt}</figcaption></figure>'
        e = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", img, e)
        e = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<a href="\2">\1</a>', e)
        e = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", e)
        e = re.sub(r"`([^`]+)`", r"<code>\1</code>", e)
        m = re.match(r"^(#{1,4})\s+(.*)", e)
        if e.startswith(("- ", "* ")):
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{e[2:]}</li>")
            continue
        if in_list:
            out.append("</ul>")
            in_list = False
        if m:
            e = f"<h{len(m.group(1)) + 1}>{m.group(2)}</h{len(m.group(1)) + 1}>"
        elif e.startswith("|"):
            cells = [c.strip() for c in e.strip("|").split("|")]
            e = "" if set("".join(cells)) <= set("-: ") else "<div class='trow'>" + "".join(
                f"<span>{c}</span>" for c in cells) + "</div>"
        elif e.startswith("&gt; "):
            e = f"<blockquote>{e[5:]}</blockquote>"
        elif e.strip() and not e.startswith("<figure"):
            e = f"<p>{e}</p>"
        out.append(e)
    if in_list:
        out.append("</ul>")
    return "\n".join(out)


def _strip_env(md: str) -> str:
    """The per-run environment dump (hundreds of pip lines) and raw stdout
    fragments belong in the run directory, not in the paper's record."""
    md = re.sub(r"\n## Environment\n.*?(?=\n## |\Z)", "\n", md, flags=re.S)
    return re.sub(r"\n## Printed values missing from this run.*?(?=\n## |\Z)", "\n", md, flags=re.S)


def _notebook_digest(c: Dict[str, Any]) -> str:
    """Many analyses: a table (one line each) and the failures grouped by
    cause, not every notebook's report inline."""
    ents = c.get("entries") or []
    groups: Dict[str, List[str]] = {}
    causes = [("a missing input file", r"FileNotFoundError|No such file|Unable to open file|does not exist"),
              ("a missing module", r"ModuleNotFoundError|ImportError|no package called"),
              ("a missing column or key", r"KeyError|not in index|undefined columns"),
              ("did not execute", r"did not execute")]
    for e in ents:
        m = e.get("first_error") or ""
        if m:
            cause = next((k for k, rx in causes if re.search(rx, m)), "another error")
            groups.setdefault(cause, []).append(e.get("entry") or "")
    clean = [e for e in ents if (e.get("chunks") or {}).get("total") and not (e.get("chunks") or {}).get("root_errors")]
    L = [f"{len(ents)} analyses were run unmodified; **{len(clean)}** ran without errors"
         + (f" ({', '.join('`' + Path(e['entry']).name + '`' for e in clean[:6])})" if clean else "") + ".", "",
         "Failures by cause:", ""]
    L += [f"- **{k}**: {len(v)} analyses (e.g. `{Path(v[0]).name}`)" for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))]
    L += ["", "| analysis | ran clean | printed values matched | first error |", "|---|---|---|---|"]
    for e in ents:
        ch, pr = e.get("chunks") or {}, e.get("printed") or {}
        got = (pr.get("identical") or 0) + (pr.get("numerically_close") or 0) + (pr.get("same_values_other_layout") or 0)
        matched = f"{got}/{pr['reference_blocks']}" if pr.get("reference_blocks") else "-"
        ok = "yes" if ch.get("total") and not ch.get("root_errors") else "no"
        err = (e.get("first_error") or "").replace("|", "/")[:90]
        L.append(f"| `{Path(e.get('entry') or '').name}` | {ok} | {matched} | {err} |")
    L.append(f"\nFull per-analysis reports: `{c.get('report_md')}`.")
    return "\n".join(L)


def render_html(rec: Dict[str, Any]) -> str:
    p, j, c = rec["paper"], rec["job"], rec.get("code")
    esc = _html.escape
    h: List[str] = []
    h.append(f"<h1>{esc(p.get('title') or rec['paper_id'])}</h1>")
    meta = " · ".join(esc(str(x)) for x in (p.get("journal"), p.get("year")) if x)
    if p.get("doi"):
        meta += f" · <a href='https://doi.org/{esc(p['doi'])}'>doi:{esc(p['doi'])}</a>"
    h.append(f"<p class='meta'>{meta}</p>")
    h.append(f"<div class='badge {esc(rec['outcome'].replace(' ', '_'))}'>{_ICON.get(rec['outcome'], '')} "
             f"{esc(rec['outcome'])}</div><p><b>{esc(headline(rec))}</b></p>")
    v = j.get("verdict") or {}
    has_claims = bool(c and c.get("claims"))  # then the claims report is the result; notebook tallies are an appendix
    rows = [("Route", rec["route"].replace("_", " ").replace("+", " + ")),
            ("Job", f"{j['id']} · {j['status']}" + (f" · {j['rounds']} round(s)" if j.get("rounds") else "")
                    + f" · {j.get('orchestrator')}"
                    + (f" · ${float(j['cost_usd']):.2f}" if j.get("cost_usd") else "")),
            ("Owner", j.get("owner") or "local"), ("Finished", j.get("finished_at") or "-"),
            ("Independent verifier", (v.get("verdict") or ("- (standalone run)" if j.get("standalone") else "-"))
             + ("; " + "; ".join(v.get("issues") or []) if v.get("issues") else ""))]
    if c:
        rows.insert(1, ("Authors' code", f"{c.get('repo')} @ {(c.get('commit') or '')[:12]} · {c.get('entry')}"
                                         + (f" · shims: {', '.join(c['shims'])}" if c.get("shims") else "")))
    h.append("<table class='kv'>" + "".join(f"<tr><th>{esc(k)}</th><td>{esc(str(v_))}</td></tr>" for k, v_ in rows)
             + "</table>")
    if c and c.get("printed") and not has_claims:
        pr = c["printed"]
        h.append("<h2>Printed values against the authors' rendering</h2><table class='grid'><tr>"
                 "<th>reference blocks</th><th>identical</th><th>numerically close</th><th>same values, other layout</th>"
                 "<th>differ</th><th>matched</th></tr><tr>" + "".join(
                     f"<td>{esc(str(pr.get(k)))}</td>" for k in ("reference_blocks", "identical", "numerically_close",
                                                                 "same_values_other_layout", "missing"))
                 + f"<td><b>{(pr.get('fraction_matched') or 0) * 100:.1f}%</b></td></tr></table>")
    for b in rec.get("data") or []:
        h.append(f"<h2>Benchmark checks: {esc(str(b['benchmark']))} ({b['n_passed']}/{b['n_total']})</h2>"
                 "<table class='grid'><tr><th></th><th>check</th><th>detail</th></tr>" + "".join(
                     f"<tr><td>{'✅' if ch.get('passed') else '❌'}</td><td>{esc(str(ch.get('name')))}</td>"
                     f"<td>{esc(str(ch.get('detail')))}</td></tr>" for ch in b["checks"]) + "</table>")
    if p.get("claims"):
        h.append("<h2>Claims harvested from the paper</h2><p class='meta'>Candidates extracted from the text; the "
                 "report below states which were reproduced.</p><table class='grid'><tr><th>value</th><th>unit</th>"
                 "<th>where</th><th>quote</th></tr>" + "".join(
                     f"<tr><td>{esc(str(cl.get('value')))}</td><td>{esc(str(cl.get('unit') or ''))}</td>"
                     f"<td>{esc(str(cl.get('section') or cl.get('source')))}</td><td>{esc(cl.get('quote') or '')}</td></tr>"
                     for cl in p["claims"][:25]) + "</table>")
    if rec["stages"]:
        h.append("<h2>Plan (harness-verified)</h2><table class='grid'><tr><th></th><th>stage</th><th>harness</th></tr>"
                 + "".join(
            f"<tr><td>{ {'done': '✅', 'blocked': '⛔', 'failed': '❌'}.get(s['status'], '▫️') }</td>"
            f"<td><b>{esc(s['id'])}</b> {esc(s.get('title') or '')}</td><td>{esc(str(s.get('reason') or ''))}"
            + (" — " + esc("; ".join(s["notes"])) if s.get("notes") else "") + "</td></tr>" for s in rec["stages"])
                 + "</table>")
    if rec.get("answer"):
        h.append("<h2>Final report</h2><div class='answer'>" + _md_html(rec["answer"]) + "</div>")
    if c:
        rd = _abs(c["run_dir"])
        claims_md = rd / "REPRODUCTION_REPORT.md"
        if claims_md.exists():
            # the paper's results against ours: verdicts, metrics, deviations
            h.append("<div class='code'>" + _md_html(claims_md.read_text(), base=rd) + "</div>")
        elif len(c.get("entries") or []) > 1:
            h.append("<h2>The repository's analyses</h2>" + _md_html(_notebook_digest(c), base=rd))
        elif _abs(c["report_md"]).exists():
            h.append("<h2>Authors' code, section by section</h2><div class='code'>"
                     + _md_html(_strip_env(_abs(c["report_md"]).read_text()), base=rd) + "</div>")
    if rec.get("log"):
        h.append("<details><summary>Selected run log</summary><pre>" + "\n".join(
            esc(f"{e['t']}  {e['kind']:<8} {e['text']}" + (f"  issues: {'; '.join(e['issues'])}" if e.get("issues") else ""))
            for e in rec["log"]) + "</pre></details>")
    h.append(f"<p class='meta'>Recorded {esc(rec['recorded_at'])} by IGVF Agent · "
             f"<a href='{PUBLIC_URL}/?repro={urllib.parse.quote(rec['paper_id'])}'>open on igvfagent.genohub.org</a></p>")
    css = """
:root{--bg:#fff;--fg:#1d1d1f;--mut:#666;--line:#ddd;--card:#f6f7f9}
@media (prefers-color-scheme:dark){:root{--bg:#141517;--fg:#e8e8ea;--mut:#9a9aa0;--line:#333;--card:#1e1f22}}
body{background:var(--bg);color:var(--fg);font-family:system-ui,-apple-system,sans-serif;max-width:1150px;margin:24px auto;padding:0 16px;line-height:1.5}
a{color:#3b82f6}.meta{color:var(--mut);font-size:14px}h1{font-size:26px;margin-bottom:4px}
.badge{display:inline-block;padding:3px 10px;border-radius:12px;background:var(--card);font-weight:600;margin:8px 0}
.badge.reproduced{background:#dcfce7;color:#14532d}.badge.partial{background:#fef9c3;color:#713f12}.badge.not_reproduced{background:#fee2e2;color:#7f1d1d}
table{border-collapse:collapse;margin:8px 0 16px;font-size:14px;width:100%}th,td{border:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}
th{background:var(--card)}.kv th{width:190px}figure{display:inline-block;margin:6px;max-width:31%;vertical-align:top}
figure img{max-width:100%;border:1px solid var(--line);background:#fff}figcaption{font-size:11px;color:var(--mut);word-break:break-all}
pre{background:var(--card);padding:10px;overflow:auto;font-size:12px}blockquote{color:var(--mut);border-left:3px solid var(--line);margin:0;padding-left:10px}
.trow{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;font-size:13px}
@media (max-width:700px){figure{max-width:100%}.kv th{width:auto}}"""
    return (f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{esc(rec['paper_id'])} reproduction</title><style>{css}</style></head><body>" + "\n".join(h)
            + "</body></html>")


# ─── forum (discussion.genohub.org) ─────────────────────────────────────────

def _forum_key() -> str:
    k = os.environ.get("IGVF_DISCOURSE_API_KEY", "").strip()
    if not k:
        f = ROOT / "Docs" / "Secret" / "discourse-API.txt"
        if f.exists():
            k = f.read_text().strip()
    return k


def key_figures(rec: Dict[str, Any], n: int = 6) -> "List[Path]":
    """One figure per section of the authors'-code report, up to n."""
    c = rec.get("code")
    if not c or not _abs(c["report_md"]).exists():
        return []
    rd = _abs(c["report_md"]).parent
    out: List[Path] = []
    for sec in re.split(r"\n## ", _abs(c["report_md"]).read_text()):
        m = re.search(r"!\[[^\]]*\]\(([^)]+)\)", sec)
        if m:
            p = rd / urllib.parse.unquote(m.group(1))
            if p.is_file() and p.suffix.lower() == ".png" and p.stat().st_size < 3_000_000:
                out.append(p)
        if len(out) >= n:
            break
    return out


def forum_post_text(rec: Dict[str, Any], uploads: "Optional[Dict[str, str]]" = None) -> Dict[str, str]:
    p, j, c = rec["paper"], rec["job"], rec.get("code")
    link = f"{PUBLIC_URL}/?repro={urllib.parse.quote(rec['paper_id'])}"
    title = f"Reproduction: {(p.get('title') or rec['paper_id'])[:180]}"
    L = [f"**{_ICON.get(rec['outcome'], '')} {rec['outcome'].capitalize()}** — {headline(rec)}", ""]
    if p.get("doi"):
        L.append(f"Paper: https://doi.org/{p['doi']}" + (f" ({p.get('journal')}, {p.get('year')})" if p.get("journal") else ""))
    if c:
        L.append(f"Authors' code: [{c['repo']}](https://github.com/{c['repo']}/tree/{c.get('commit')}) · `{c.get('entry')}`, "
                 "run unmodified" + (f" (compatibility shims: {', '.join(c['shims'])})" if c.get("shims") else ""))
    L += [f"Route: {rec['route'].replace('_', ' ').replace('+', ' + ')} · job `{j['id']}` · {j.get('rounds')} round(s)"
          f" · independent verifier: **{(j.get('verdict') or {}).get('verdict', '-')}**",
          f"Full record with every figure: {link}", ""]
    if c and c.get("printed"):
        pr = c["printed"]
        L += ["| printed-value blocks | identical | numerically close | same values, other layout | differ |",
              "|---|---|---|---|---|",
              f"| {pr.get('reference_blocks')} | {pr.get('identical')} | {pr.get('numerically_close')} | "
              f"{pr.get('same_values_other_layout')} | {pr.get('missing')} |", ""]
    for b in rec.get("data") or []:
        L.append(f"**{b['benchmark']}**: {b['n_passed']}/{b['n_total']} benchmark checks")
        L += [f"- {'✅' if ch['passed'] else '❌'} {ch['name']}: {ch['detail']}" for ch in b["checks"][:12]]
        L.append("")
    blocked = [s for s in rec["stages"] if s["status"] in ("blocked", "failed")]
    if blocked:
        L += ["**Blocked or failed stages**"] + [f"- {s['id']}: {'; '.join(s.get('notes') or []) or s.get('reason')}"
                                                 for s in blocked] + [""]
    for name, url in (uploads or {}).items():
        L.append(f"![{name}]({url})")
    if rec.get("log"):
        L += ["", "[details=\"Selected run log\"]", "```"] + [
            f"{e['t']}  {e['kind']:<8} {e['text'][:200]}" for e in rec["log"][-40:]] + ["```", "[/details]"]
    L += ["", "<small>Posted from IGVF Agent (igvfagent.genohub.org).</small>"]
    return {"title": title, "raw": "\n".join(L)}


def _forum_request(path: str, *, user: str, data: Optional[bytes] = None, ctype: str = "application/json") -> Any:
    req = urllib.request.Request(FORUM_URL + path, data=data, method="POST" if data is not None else "GET",
                                 headers={"Api-Key": _forum_key(), "Api-Username": user, "Content-Type": ctype,
                                          "Accept": "application/json", "User-Agent": "igvfagent"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode() or "{}")


def _upload(path: Path, user: str) -> str:
    b = uuid.uuid4().hex
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"type\"\r\n\r\ncomposer\r\n"
            f"--{b}\r\nContent-Disposition: form-data; name=\"synchronous\"\r\n\r\ntrue\r\n"
            f"--{b}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"{path.name}\"\r\n"
            f"Content-Type: image/png\r\n\r\n").encode() + path.read_bytes() + f"\r\n--{b}--\r\n".encode()
    r = _forum_request("/uploads.json", user=user, data=body, ctype=f"multipart/form-data; boundary={b}")
    return r.get("short_url") or r.get("url")


def post_to_forum(paper_id: str, *, as_user: Optional[str] = None, yes: bool = False) -> Dict[str, Any]:
    """Preview (default) or post the record to the forum category. A second
    post for the same paper is a reply in the same topic."""
    d = REPRO_DIR / _slug(paper_id)
    rec = _read_json(d / "record.json")
    if not rec:
        return {"ok": False, "error": f"no record {paper_id}"}
    user = (as_user or rec["job"].get("owner") or os.environ.get("IGVF_DISCOURSE_API_USER") or "system").strip()
    figs = key_figures(rec)
    topic = (rec.get("forum") or {}).get("topic_id")
    if not yes:
        post = forum_post_text(rec, {f.name: f"(upload) {f.name}" for f in figs})
        return {"ok": True, "preview": True, "as_user": user, "category": FORUM_CATEGORY,
                "reply_to_topic": topic, "figures": [_rel(f) for f in figs], **post}
    if not _forum_key():
        return {"ok": False, "error": "no forum API key (IGVF_DISCOURSE_API_KEY)"}
    try:
        uploads = {f.name: _upload(f, user) for f in figs}
        post = forum_post_text(rec, uploads)
        payload = {"raw": post["raw"]}
        if topic:
            payload["topic_id"] = topic
        else:
            payload.update({"title": post["title"], "category": FORUM_CATEGORY})
        r = _forum_request("/posts.json", user=user, data=json.dumps(payload).encode())
    except urllib.error.HTTPError as e:  # type: ignore[attr-defined]
        return {"ok": False, "error": f"forum HTTP {e.code}: {e.read().decode(errors='replace')[:300]}"}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    tid, slug = r.get("topic_id"), r.get("topic_slug")
    url = f"{FORUM_URL}/t/{slug}/{tid}/{r.get('post_number', 1)}"
    rec["forum"] = {"topic_id": tid, "topic_url": f"{FORUM_URL}/t/{slug}/{tid}",
                    "posts": (rec.get("forum") or {}).get("posts", []) + [
                        {"url": url, "job": rec["job"]["id"], "by": user,
                         "at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}]}
    _write_json(d / "record.json", rec)
    return {"ok": True, "url": url, "topic_url": rec["forum"]["topic_url"], "as_user": user}


# ─── Streamlit panel ────────────────────────────────────────────────────────

def render_streamlit_panel(st, viewer: Optional[str], is_admin: bool, focus: Optional[str] = None) -> None:
    st.markdown("#### 📑 Reproductions")
    st.caption("One record per paper, kept across attempts: how it was reproduced, agreement with the paper and the "
               "authors' outputs, the verifier's verdict, every figure. Private to you until published.")
    q = st.text_input("Search", key="repro_q", placeholder="title, DOI, repository, gene, assay, outcome",
                      label_visibility="collapsed")
    recs = list_records(viewer, is_admin, q)
    if not recs:
        st.info("No reproduction records yet. They are written when a reproduction or benchmark job finishes.")
        return
    counts = {k: sum(1 for r in recs if r["outcome"] == k) for k in ("reproduced", "partial", "not reproduced")}
    st.caption(f"{len(recs)} paper(s): ✅ {counts['reproduced']} reproduced · 🟡 {counts['partial']} partial · "
               f"❌ {counts['not reproduced']} not reproduced")
    rows = [{"": _ICON.get(r["outcome"], ""), "paper": (r["paper"].get("title") or r["paper_id"])[:90],
             "year": r["paper"].get("year"), "route": r["route"].replace("_", " "), "result": headline(r)[:120],
             "verifier": (r["job"].get("verdict") or {}).get("verdict", ""), "attempts": len(r.get("attempts") or []),
             "owner": r["job"].get("owner"), "published": "🌐" if r.get("published") else "",
             "forum": "💬" if (r.get("forum") or {}).get("topic_url") else "", "id": r["paper_id"]} for r in recs]
    st.dataframe(rows, hide_index=True, width="stretch")
    ids = [r["paper_id"] for r in recs]
    idx = ids.index(focus) if focus in ids else 0
    pick = st.selectbox("Open", ids, index=idx, key="repro_pick",
                        format_func=lambda i: next((f"{_ICON.get(r['outcome'], '')} {(r['paper'].get('title') or i)[:100]}"
                                                    for r in recs if r["paper_id"] == i), i))
    rec = next(r for r in recs if r["paper_id"] == pick)
    render_record(st, rec, viewer, is_admin)


def render_record(st, rec: Dict[str, Any], viewer: Optional[str], is_admin: bool, key: str = "tab") -> None:
    """One record. `key` separates copies: the same record can be drawn twice on
    a page (the ?repro= banner and the Reproductions tab)."""
    mine = is_admin or viewer is None or viewer in (rec.get("owners") or [])
    attempts = rec.get("attempts") or []
    choice = rec["job"]["id"]
    if len(attempts) > 1:
        choice = st.selectbox("Attempt", [a["job"] for a in reversed(attempts)], key=f"repro_{key}_att_{rec['paper_id']}",
                              format_func=lambda jid: next(f"{a['created_at']} · {a['outcome']} · {a['route']} · {jid}"
                                                           for a in attempts if a["job"] == jid))
    att = _read_json(REPRO_DIR / rec["paper_id"] / "attempts" / f"{choice}.json") or rec
    html_p = _abs(att.get("html") or "")
    c1, c2, c3 = st.columns([1, 1, 2])
    if html_p.is_file():
        c1.download_button("⬇ Report (HTML)", html_p.read_bytes(), file_name=f"{rec['paper_id']}_{choice}.html",
                           mime="text/html", key=f"repro_{key}_dl_{rec['paper_id']}_{choice}")
    if mine:
        if c2.button("🌐 Unpublish" if rec.get("published") else "🌐 Publish", key=f"repro_{key}_pub_{rec['paper_id']}",
                     help="Published records are visible to every signed-in user"):
            set_published(rec["paper_id"], not rec.get("published"))
            st.rerun()
    spec = _read_json(REPRO_DIR / rec["paper_id"] / "tool.json")
    if spec:
        st.info(f"♻️ **Reusable analysis: `{spec['name']}`.** The authors' code for this paper can run on your data. "
                f"Ask the agent, e.g. *\"run {spec['name']} with my file in place of {next(iter(spec['inputs'].values()), 'the input')}\"*, "
                f"or: `igvfagent repro apply {rec['paper_id']} --{next(iter(spec['inputs']), 'input').replace('_', '-')} FILE`. "
                f"Inputs: {', '.join('`' + v + '`' for v in spec['inputs'].values())}.")
    fr = rec.get("forum") or {}
    if fr.get("topic_url"):
        c3.markdown(f"💬 [Forum topic]({fr['topic_url']})")
    if mine:
        with st.expander("💬 Post to discussion.genohub.org (IGVF Agent Help)", expanded=False):
            prev = post_to_forum(rec["paper_id"], as_user=viewer, yes=False)
            st.caption(f"Posts as **{prev.get('as_user')}** in category {prev.get('category')}"
                       + (f", as a reply in the existing topic" if prev.get("reply_to_topic") else ", as a new topic")
                       + f", with {len(prev.get('figures') or [])} key figure(s). The forum is public: check the text.")
            st.markdown(f"**{prev.get('title')}**")
            st.code(prev.get("raw", "")[:6000], language="markdown")
            ok = st.checkbox("I have read this and want to post it", key=f"repro_{key}_ok_{rec['paper_id']}")
            if st.button("Post", key=f"repro_{key}_post_{rec['paper_id']}", disabled=not ok):
                res = post_to_forum(rec["paper_id"], as_user=viewer, yes=True)
                if res.get("ok"):
                    st.success(f"Posted: {res['url']}")
                else:
                    st.error(res.get("error"))
    if html_p.is_file():
        try:
            import streamlit.components.v1 as components
            components.html(html_p.read_text(), height=1100, scrolling=True)
        except Exception:  # noqa: BLE001
            st.markdown(att.get("answer") or "")


# ─── CLI ────────────────────────────────────────────────────────────────────

def selftest() -> int:
    import tempfile
    global ROOT, REPRO_DIR, JOBS_DIR
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")

    saved = (ROOT, REPRO_DIR, JOBS_DIR)
    with tempfile.TemporaryDirectory() as td:
        td = str(Path(td).resolve())
        ROOT, REPRO_DIR, JOBS_DIR = Path(td), Path(td) / "Data" / "Reproductions", Path(td) / "Data" / "Jobs"
        try:
            bd = ROOT / "Docs" / "Benchmark" / "r1"
            bd.mkdir(parents=True)
            _write_json(bd / "resolution.json", {"paper_id": "matreyek2018_multiplex", "paper": {
                "title": "Multiplex assessment of protein variant abundance", "doi": "10.1038/s41588-018-0122-z",
                "journal": "Nat Genet", "year": 2018}})
            _write_json(bd / "harvest.json", {"numeric_claims": [{"value": 8000, "unit": "variants", "quote": "q"}],
                                             "accessions": {"github_repo": [{"value": "FowlerLab/VAMPseq",
                                                                             "in_data_availability": True}]}})
            pc = ROOT / "Docs" / "PaperCode" / "p1"
            (pc / "figures").mkdir(parents=True)
            png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8DwHwAFBQIAX8jx0gAAAABJRU5ErkJggg==")
            (pc / "figures" / "f1.png").write_bytes(png)
            (pc / "report.md").write_text("# R\n\n## 1. Setup\n\n![f1](figures/f1.png)\n")
            _write_json(pc / "summary.json", {"repo": "FowlerLab/VAMPseq", "commit": "a960aa9b", "entry": "a.Rmd",
                                              "chunks": {"total": 67, "errored": 0, "root_errors": 0},
                                              "figures": {"produced": 119}, "replay": {"deterministic": True},
                                              "printed": {"reference_blocks": 110, "identical": 106,
                                                          "numerically_close": 0, "same_values_other_layout": 4,
                                                          "missing": 0, "fraction_matched": 1.0}})
            _write_json(pc / "run.json", {"r_shims": ["melt_reshape"]})
            _write_json(pc / "inventory.json", {"ok": True, "inputs": [
                {"path": "PTEN_variant_data.tsv", "exists": True, "is_url": False},
                {"path": "https://x/y.tsv", "exists": False, "is_url": True}],
                "sections": [{"title": "Setup"}]})
            _write_json(pc / "env.json", {"spec": ["r-base", "cran:reshape2@1.4.4"]})
            jd = JOBS_DIR / "J20260924000000abcd"
            _write_json(jd / "job.json", {"id": jd.name, "title": "Reproduce Matreyek 2018", "query": "Reproduce it",
                                          "owner": "alice", "status": "done", "rounds": 2,
                                          "created_at": "2026-09-24T00:00:00Z", "last_verdict": {"verdict": "pass"}})
            _write_json(jd / "plan.json", {"stages": [
                {"id": "harvest", "title": "h", "status": "done", "evidence": ["Docs/Benchmark/r1/harvest.json"]},
                {"id": "run_code", "title": "c", "status": "done", "evidence": ["Docs/PaperCode/p1/summary.json"]}]})
            (jd / "events.jsonl").write_text(json.dumps({"t": "x", "kind": "stage", "text": "stage run_code verified"})
                                             + "\n" + json.dumps({"t": "y", "kind": "tool", "text": "noise"}) + "\n")
            (jd / "answer.md").write_text("## Done\n\n- all good")
            rec = record_from_job(jd.name)
            check("record built from the job's evidence", rec and rec["paper_id"] == "matreyek2018_multiplex"
                  and rec["route"] == "authors_code" and rec["outcome"] == "reproduced")
            check("headline states the agreement", "110/110 printed values" in headline(rec) and "67/67 chunks" in headline(rec))
            html = (ROOT / rec["html"]).read_text()
            check("HTML is self-contained (figure inlined)", "data:image/png;base64," in html and "Authors' code, section" in html
                  and not Path(rec["html"]).is_absolute())
            check("selected log keeps stage events, drops tool chatter", [e["kind"] for e in rec["log"]] == ["stage"])
            check("private to its owner", list_records("bob", False) == [] and len(list_records("alice", False)) == 1)
            set_published(rec["paper_id"], True)
            check("published records are visible to everyone", len(list_records("bob", False)) == 1)
            check("search matches DOI and repository", len(list_records(None, True, "s41588 vampseq")) == 1
                  and not list_records(None, True, "zzz"))
            jd2 = JOBS_DIR / "J20260925000000abcd"
            _write_json(jd2 / "job.json", {**_read_json(jd / "job.json"), "id": jd2.name, "status": "failed",
                                           "created_at": "2026-09-25T00:00:00Z", "last_verdict": {}})
            _write_json(jd2 / "plan.json", {"stages": [{"id": "harvest", "title": "h", "status": "done",
                                                        "evidence": ["Docs/Benchmark/r1/harvest.json"]}]})
            record_from_job(jd2.name)
            r = load("matreyek2018_multiplex")
            check("attempts accumulate; the best one headlines; published flag kept",
                  len(r["attempts"]) == 2 and r["job"]["id"] == jd.name and r["latest_job"] == jd2.name and r["published"])
            t = paper_tools()
            check("a reproduced paper becomes one reusable tool with its inputs and pins",
                  len(t) == 1 and t[0]["name"] == "paper_matreyek2018_multiplex"
                  and t[0]["inputs"] == {"pten_variant_data_tsv": "PTEN_variant_data.tsv"}
                  and t[0]["pins"] == ["cran:reshape2@1.4.4"] and t[0]["cli"][:3] == ["repro", "apply", "matreyek2018_multiplex"])
            check("apply refuses unknown inputs and missing files",
                  not apply("matreyek2018_multiplex", {"nope": "x"})["ok"]
                  and not apply("matreyek2018_multiplex", {"pten_variant_data_tsv": "missing.tsv"})["ok"])
            prev = post_to_forum("matreyek2018_multiplex", yes=False)
            check("forum post is a preview unless confirmed, authored as the owner",
                  prev["preview"] and prev["as_user"] == "alice" and prev["category"] == FORUM_CATEGORY
                  and "110" in prev["raw"] and "?repro=matreyek2018_multiplex" in prev["raw"])
        finally:
            ROOT, REPRO_DIR, JOBS_DIR = saved
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent repro", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("record")
    s.add_argument("job_id")
    sub.add_parser("backfill")
    s = sub.add_parser("record-run", help="record a standalone paper-code run directory")
    s.add_argument("run_dir")
    s.add_argument("--paper-id")
    s.add_argument("--doi")
    s.add_argument("--title")
    s.add_argument("--owner")
    s.add_argument("--note", help="findings to attach (Markdown), shown as the attempt's report text")
    s.add_argument("--note-file")
    s = sub.add_parser("list")
    s.add_argument("--query", default="")
    s = sub.add_parser("show")
    s.add_argument("paper_id")
    s = sub.add_parser("publish")
    s.add_argument("paper_id")
    s.add_argument("--off", action="store_true")
    s = sub.add_parser("post")
    s.add_argument("paper_id")
    s.add_argument("--as-user")
    s.add_argument("--yes", action="store_true", help="actually post (default: preview only)")
    sub.add_parser("tools", help="reusable per-paper analyses (paper_<id> tools)")
    s = sub.add_parser("apply", help="run a reproduced paper's analysis on new data: --<input> PATH ...")
    s.add_argument("paper_id")
    s.add_argument("--replay", action="store_true")
    s.add_argument("--detach", action="store_true")
    s.add_argument("--label", default="")
    sub.add_parser("selftest")
    a, extra = ap.parse_known_args(argv)
    if extra and a.cmd != "apply":
        ap.error(f"unrecognized arguments: {' '.join(extra)}")
    viewer = os.environ.get("IGVF_ACTING_USER") or None
    admin = os.environ.get("IGVF_ACTING_ADMIN", "1" if not viewer else "0") == "1"
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "record":
        r = record_from_job(a.job_id)
        if not r:
            print(f"no job {a.job_id}")
            return 1
        print(f"{r['paper_id']}: {r['outcome']} — {headline(r)}")
        print(f"Report: {r['html']}")
        return 0
    if a.cmd == "tools":
        for t in paper_tools():
            print(f"{t['name']:<44} {t['repo']}@{(t['commit'] or '')[:10]}  inputs: {', '.join(t['inputs'])}")
        return 0
    if a.cmd == "apply":
        files, it = {}, iter(extra)
        for tok in it:
            if not tok.startswith("--"):
                print(f"unexpected argument {tok!r}")
                return 2
            k, _, v = tok[2:].partition("=")
            files[k.replace("-", "_")] = v or next(it, "")
        r = apply(a.paper_id, files, replay=a.replay, detach=a.detach, label=a.label)
        if not r.get("ok") and r.get("error"):
            print(r["error"])
        return 0 if r.get("ok") else 1
    if a.cmd == "record-run":
        note = Path(a.note_file).read_text() if a.note_file else (a.note or "")
        r = record_from_run(a.run_dir, paper_id=a.paper_id, doi=a.doi, title=a.title, owner=a.owner or "",
                            note=note)
        if not r:
            print(f"no paper-code summary in {a.run_dir}")
            return 1
        print(f"{r['paper_id']}: {r['outcome']} — {headline(r)}")
        print(f"Report: {r['html']}")
        return 0
    if a.cmd == "backfill":
        for ln in backfill():
            print(ln)
        return 0
    if a.cmd == "list":
        for r in list_records(viewer, admin, a.query):
            print(f"{_ICON.get(r['outcome'], ' ')} {r['paper_id']:<40} {r['route']:<24} {headline(r)[:90]}")
        return 0
    if a.cmd == "show":
        r = load(a.paper_id)
        if not r or not visible(r, viewer, admin):
            print(f"no record {a.paper_id}")
            return 1
        print(json.dumps({k: v for k, v in r.items() if k not in ("answer", "log")}, indent=2, default=str)[:12000])
        print(f"Report: {r.get('html')}")
        return 0
    if a.cmd == "publish":
        print(json.dumps(set_published(a.paper_id, not a.off)))
        return 0
    if a.cmd == "post":
        r = post_to_forum(a.paper_id, as_user=a.as_user, yes=a.yes)
        print(json.dumps(r, indent=2)[:8000])
        return 0 if r.get("ok") else 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
