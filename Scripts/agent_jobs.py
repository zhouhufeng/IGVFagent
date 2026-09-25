"""Durable agent jobs: long analyses that keep working until they are verifiably done.

A chat turn is the wrong unit for "reproduce this paper". Before this module,
every message was one bounded agent run, with a 25-iteration cap on the hosted
site. The run lived inside the web request, and its plan lived only in the
reply text. So a reproduction ended "complete" after one pass ("here's where
the reproduction stands"). A follow-up "continue" started cold, and in one
recorded case the agent replied it had no context at all. Closing the tab
ended the run.

A job is the Paper2Agent-style unit instead:

- **Durable.** A job runs in its own process (``igvfagent job worker``),
  not in the web request. Closing the browser does not stop it. Its whole
  state is on disk under ``Data/Jobs/<id>/``, so any process can resume it
  after a crash or a redeploy.
- **Planned on disk.** The agent writes a staged plan with ``plan_set``: id,
  goal, success criteria, evidence and an optional check. It marks progress
  with ``plan_update``. The plan, not the chat text, is what "continue"
  resumes from.
- **Gated.** A stage marked done is verified by the harness, not by the model.
  Its evidence files must exist in the workspace and its check (files,
  json value, table rows, benchmark concordance) must pass. Otherwise the
  stage reverts to failed, with the reason fed back.
- **Rounds, not a turn.** Each round is a full agent run. After it, the harness
  checks the plan. If stages are pending and the budget allows, another round
  starts with a compact state message: the plan, the harness findings and the
  last answer. This compacts context the way a long session does. Budgets
  are wall-clock minutes and rounds, not a small iteration count.
- **Repair loop.** A failed stage is retried up to three times, then the agent
  must mark it blocked with a reason. A blocked stage is reported, never
  passed off as done.
- **Independent verification.** When every stage is done, a fresh model call
  with no shared context reviews the task, the plan, the evidence and the
  final answer, and returns pass or fail with issues. On fail, the issues
  start a repair round, up to two verification cycles.
- **Delegation.** ``delegate_tasks`` runs up to four focused sub-agents in
  parallel, each in a fresh context. Each returns a summary and its
  artefacts. This is the Task-tool pattern.
- **Waiting.** ``job_wait`` blocks on a detached raw-pipeline job or an output
  file (up to hours). A job can therefore start long background work and use
  its result, instead of ending its turn and forgetting it.

Two orchestrators drive the rounds:

- ``internal``: IGVFagent's own loop (``_agent.run``), with the job protocol
  appended to its system prompt and a higher per-round iteration budget.
- ``claude_code``: Claude Code headless (``claude -p``) as the orchestrator,
  authenticated by ``ANTHROPIC_API_KEY``. Every IGVFagent tool is connected
  over MCP (``igvfagent mcp serve``). The session is resumed between rounds
  with ``--resume``, and the same plan gates and verifier apply. On a shared
  deployment it gets IGVFagent's tools plus read-only file tools. It gets no
  shell unless the operator allows agent authoring.

CLI::

    igvfagent job start "reproduce Matreyek 2018 PTEN VAMP-seq" [--orchestrator internal|claude_code]
    igvfagent job list | status JOB | logs JOB | stop JOB | resume JOB
    igvfagent job plan-set --stages-json '[...]'      (inside a job; the agent calls these)
    igvfagent job plan-update --stage s1 --status done --evidence Docs/.../report.md
    igvfagent job plan-show
    igvfagent job delegate --task "..." [--task "..."]
    igvfagent job wait --path Docs/.../done.flag | --raw-job RJ123 [--timeout-min 240]
    igvfagent job selftest
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
JOBS_DIR = Path(os.environ.get("IGVF_JOBS_DIR") or ROOT / "Data" / "Jobs")
log = logging.getLogger("agent_jobs")

TERMINAL = ("done", "done_with_blocked", "failed", "stopped", "budget_exhausted")
WAIT_MAX_MIN = float(os.environ.get("IGVF_JOB_WAIT_MAX_MIN", "240"))
RESUMABLE = ("stopped", "interrupted", "budget_exhausted", "failed", "done_with_blocked")
STAGE_STATES = ("pending", "running", "done", "failed", "blocked")
MAX_STAGE_ATTEMPTS = int(os.environ.get("IGVF_JOB_STAGE_ATTEMPTS", "5"))
MAX_VERIFY_CYCLES = 2
DEFAULT_BUDGET_MIN = float(os.environ.get("IGVF_JOB_BUDGET_MIN", "240"))
DEFAULT_MAX_ROUNDS = int(os.environ.get("IGVF_JOB_MAX_ROUNDS", "12"))
ROUND_ITERATIONS = int(os.environ.get("IGVF_JOB_ROUND_ITERATIONS", "40"))
HEARTBEAT_S = 20
STALE_S = 180
MAX_RUNNING = int(os.environ.get("IGVF_MAX_AGENT_JOBS", "4"))
MAX_DELEGATES = 4

JOB_PROTOCOL = """\
You are running as a LONG-LIVED IGVFagent JOB, not a single chat reply. You may
take many rounds and hours. There is no need to finish in one reply, and ending
early is a failure, not a courtesy.

1. For any task with more than one step, FIRST call plan_set with 3-10 stages.
   Each stage: id, title, goal, success (concrete criteria), evidence (the
   files the stage will produce) and, where possible, a check the harness can
   run: {"kind":"files","paths":[...]} | {"kind":"json","path":P,"key":"a.b",
   "op":">=","value":X} | {"kind":"rows","path":P,"min":N} |
   {"kind":"concordance","benchmark":NAME}.
2. Work ONE stage at a time. When a stage is finished, call plan_update with
   status "done" and the evidence paths you produced. The harness verifies
   the evidence and runs the check. If they fail, the stage reverts to failed
   and you will be told why.
3. When something fails, diagnose it and try a different approach. After 5
   failed attempts, mark the stage "blocked" with the concrete reason
   (missing data, needs credentials, unsupported method). Never mark a stage
   done that was skipped or failed.
4. Long downloads and pipelines are fine: start them, then call job_wait on
   the detached job or output file instead of ending the round.
5. Independent sub-questions can run in parallel with delegate_tasks. Each
   sub-agent starts fresh, so give it everything it needs in the task text.
5b. To save notes, section verdicts or a report as stage evidence, use
   write_text_file (Docs/, Data/ or Benchmarks/ paths) and read_text_file.
   Do not author new tools or extensions to write files, and do not use
   sed edits or sub-agents as a file writer.
6. When every stage is done or blocked, write the final report: what was
   reproduced (with numbers and file paths), what differs from the paper and
   why, and what is blocked.

REPRODUCING A PAPER: CLAIMS FIRST (as Paper2Agent does). A reproduction is
judged on the paper's RESULTS, not on how many notebooks run. Figure notebooks
are plotting layers; the results come from the authors' method.
  1. List 2-5 key quantitative results and, for each, the published table
     that holds it: paper_code_fetch_supplementary gives the Supplementary
     Tables and Source Data (one sheet per result, e.g. "6. LDLvar GWAS BEAN
     result").
  2. For each, run the authors' METHOD the way their workflow does
     (paper_code_exec with the command from their Snakefile/README, on the
     deposited processed data; paper_code_fetch_data for deposits). A needed
     deviation (e.g. masking a sample that failed QC) is passed as `note`.
  3. paper_code_claim compares our table with the published one on the column
     the claim rests on; run the same command again with seed=202 and pass it
     as `noise` so "reproduced" means "within seed noise". A result the paper
     states as a number (a count, r, a percentage) is a claim too: `published`
     + the paper's `quote`, with our value read by `pattern` from what the
     authors' code printed (never typed in).
  4. When a claim falls short, diagnose it scientifically with the authors'
     own tools (rebuild an input the deposit lacks, recover parameters by
     reproducing a deposited intermediate exactly) and document it. A key
     result you could not check is still recorded (paper_code_claim with
     not_attempted = why), so the verdict never overstates coverage.
  5. paper_code_report_claims writes the report: a verdict table per claim,
     metrics, deviations. paper_code_reproduce with all_entries (running every
     notebook) is a supplementary check, not the reproduction.

REPRODUCING A PAPER (the paper's own code is the scientific source of truth)
a. Resolve and harvest the paper first (paper_benchmark / bench harvest), then
   paper_code_find on its harvest.json. If the Code Availability statement
   names a repository, the reproduction IS running that code: call
   paper_code_reproduce with replay=true, and all_entries=true when the
   repository has one notebook per figure (it runs in the background); job_wait on
   <run_dir>/done.json. Never reimplement the authors' analysis yourself, and
   do not present a re-derivation from deposited tables (MaveDB, Portal) as
   the reproduction; that is a separate cross-check stage.
b. Plan stages like: harvest -> find_code -> run_code (check: json
   <run_dir>/summary.json key printed.fraction_matched >= 0.9, or figures.produced
   >= 1 when the authors committed no rendering) -> review_sections -> repair
   (only if needed) -> cross_check (optional: IGVF/MaveDB data route) -> report.
c0. When analyses fail on missing inputs: if the repository has a Snakemake
   workflow, re-run with workflow=true (its dry run lists every missing input);
   fetch the paper's deposited data (Zenodo/figshare in the Data Availability
   statement) with paper_code_fetch_data into the paths the code reads; then
   re-run. When the workflow starts from raw reads, fetch them with
   paper_code_fetch_reads (guess_layout, then preview, then download=true)
   into the paths the authors' code reads. If the paper states software
   versions that differ from the repository's environment (e.g. "version 0.2.9
   of bean was used"), re-run with pin pip:NAME==VERSION; the paper wins.
   Missing data that is not deposited is a blocker, named in the report.
c. Read summary.json and report.md. For chunks that raised errors, diagnose
   from run.log and report.md (a package API change is the usual cause;
   inventory.json lists the versions the authors ran). Repair only the
   environment, with pins (r-<pkg>=<version> or cran:<pkg>@<version>) and
   re-run; at most 3 environment repairs. Never edit the authors' code to make
   it pass; if it cannot run, mark the stage blocked naming the chunk, line and
   error.
d. Verify section by section with delegate_tasks: one fresh sub-agent per
   group of analysis sections, each given the section titles, the figure paths
   and the printed-value agreement for those sections, asked to state which of
   the paper's claims for those sections the outputs support, contradict or
   leave untested. Their verdicts go in the final report.
e. The final report leads with the authors' repository and commit, the
   printed-value agreement (identical / numerically close / missing), the
   figures per section (paths), the errors and blocked chunks, and the
   environment differences from the authors'. Numbers come only from
   summary.json, compare.json and the run's files.
"""


# ─── job store ──────────────────────────────────────────────────────────────

def _now() -> float:
    return time.time()


def _iso(t: Optional[float] = None) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t or _now()))


def job_dir(job_id: str) -> Path:
    if not re.fullmatch(r"J[0-9A-Za-z_-]{6,40}", job_id or ""):
        raise ValueError(f"bad job id {job_id!r}")
    return JOBS_DIR / job_id


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(data, indent=2, default=str))
    os.replace(tmp, path)


def _read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def load_job(job_id: str) -> dict:
    j = _read_json(job_dir(job_id) / "job.json")
    if not j:
        raise FileNotFoundError(f"no job {job_id}")
    return j


def save_job(job: dict) -> None:
    job["updated_at"] = _iso()
    _write_json(job_dir(job["id"]) / "job.json", job)


_JOB_LOCK = threading.Lock()


def update_job(job_id: str, **fields) -> dict:
    # the heartbeat thread and the round loop both write job.json: serialise
    # read-modify-write so neither drops the other's fields
    with _JOB_LOCK:
        j = load_job(job_id)
        j.update(fields)
        save_job(j)
        return j


def event(job_id: str, kind: str, text: str = "", **payload) -> None:
    rec = {"t": _iso(), "kind": kind, "text": text[:2000], **payload}
    with (job_dir(job_id) / "events.jsonl").open("a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def read_events(job_id: str, last: int = 50) -> "List[dict]":
    p = job_dir(job_id) / "events.jsonl"
    try:
        lines = p.read_text().splitlines()[-last:]
    except OSError:
        return []
    out = []
    for l in lines:
        try:
            out.append(json.loads(l))
        except ValueError:
            continue
    return out


def _pid_alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def _container_started() -> float:
    """When this container (its PID 1) started; 0 where there is no /proc.
    A heartbeat older than that came from a worker the restart killed."""
    try:
        return os.stat("/proc/1").st_ctime if os.path.exists("/.dockerenv") else 0.0
    except OSError:
        return 0.0


def effective_status(job: dict) -> str:
    """'interrupted' when a running job's worker is gone (container restart, kill)."""
    st = job.get("status")
    if st in ("running", "queued"):
        # A worker recorded in another container (the one a redeploy replaced)
        # is gone, however fresh its last heartbeat; its pid means nothing here.
        if job.get("host") and job["host"] != socket.gethostname() and not os.environ.get("IGVF_JOBS_SHARED_HOSTS"):
            return "interrupted"
        hb = float(job.get("heartbeat") or 0)
        if not _pid_alive(job.get("pid")) and ((_now() - hb) > STALE_S or hb < _container_started()):
            return "interrupted"
    return st


def list_jobs(viewer: Optional[str] = None, is_admin: bool = True, limit: int = 50) -> "List[dict]":
    if not JOBS_DIR.is_dir():
        return []
    jobs = []
    for d in JOBS_DIR.iterdir():
        j = _read_json(d / "job.json")
        if not j:
            continue
        if viewer and not is_admin and (j.get("owner") or "") != viewer:
            continue
        j["effective_status"] = effective_status(j)
        jobs.append(j)
    jobs.sort(key=lambda j: j.get("created_at") or "", reverse=True)
    return jobs[:limit]


# ─── plan: the job's state machine, verified by the harness ─────────────────

def load_plan(job_id: str) -> dict:
    return _read_json(job_dir(job_id) / "plan.json", {"stages": [], "history": []})


def save_plan(job_id: str, plan: dict) -> None:
    _write_json(job_dir(job_id) / "plan.json", plan)


def _safe_path(p: str) -> Optional[Path]:
    try:
        import _pathguard  # type: ignore
    except Exception:  # pragma: no cover
        _pathguard = None
    q = Path(p)
    if not q.is_absolute():
        q = ROOT / q
    q = q.resolve()
    try:
        q.relative_to(ROOT)
    except ValueError:
        return None
    if _pathguard is not None and _pathguard.why_blocked(q, require_file=False) is not None:
        return None
    return q


def _dig(obj: Any, key: str) -> Any:
    for part in key.split("."):
        if isinstance(obj, dict):
            obj = obj.get(part)
        elif isinstance(obj, list) and part.isdigit() and int(part) < len(obj):
            obj = obj[int(part)]
        else:
            return None
    return obj


def run_check(check: Optional[dict], evidence: "List[str]") -> "Tuple[bool, str]":
    """Harness-side verification of a stage. Returns (passed, reason)."""
    missing = []
    for e in evidence or []:
        q = _safe_path(e)
        if q is None or not q.exists() or (q.is_file() and q.stat().st_size == 0):
            missing.append(e)
    if missing:
        return False, "evidence missing or empty: " + ", ".join(missing[:5])
    if not check:
        return (True, "evidence present") if evidence else (False, "no evidence listed for a done stage")
    kind = check.get("kind")
    if kind == "files":
        bad = [p for p in check.get("paths") or [] if not (_safe_path(p) and _safe_path(p).exists())]
        return (not bad, "all files present" if not bad else "missing: " + ", ".join(bad[:5]))
    if kind == "json":
        q = _safe_path(check.get("path", ""))
        data = _read_json(q) if q else None
        if data is None:
            return False, f"cannot read JSON {check.get('path')}"
        v = _dig(data, check.get("key", ""))
        op, want = check.get("op", "=="), check.get("value")
        try:
            ok = {"==": v == want, "!=": v != want, ">=": float(v) >= float(want), "<=": float(v) <= float(want),
                  ">": float(v) > float(want), "<": float(v) < float(want),
                  "exists": v is not None}.get(op, False)
        except (TypeError, ValueError):
            ok = False
        return ok, f"{check.get('key')} = {v!r} ({op} {want!r})"
    if kind == "rows":
        q = _safe_path(check.get("path", ""))
        if not q or not q.exists():
            return False, f"missing table {check.get('path')}"
        opener = __import__("gzip").open if q.suffix == ".gz" else open
        with opener(q, "rt") as fh:
            n = sum(1 for _ in fh) - (1 if check.get("header", True) else 0)
        return n >= int(check.get("min", 1)), f"{n} rows (need >= {check.get('min', 1)})"
    if kind == "concordance":
        name = str(check.get("benchmark", ""))
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            return False, "bad benchmark name"
        script = ROOT / "Benchmarks" / "concordance.py"
        if not script.exists():
            return False, "Benchmarks/concordance.py not found"
        p = subprocess.run([sys.executable, str(script), "--benchmark", name], cwd=str(ROOT),
                           capture_output=True, text=True, timeout=600)
        tail = (p.stdout + p.stderr).strip().splitlines()[-3:]
        return p.returncode == 0, "concordance: " + " | ".join(tail)
    return False, f"unknown check kind {kind!r}"


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def evidence_hashes(evidence: "List[str]") -> "Dict[str, str]":
    out = {}
    for e in evidence or []:
        q = _safe_path(e)
        if q is not None and q.is_file():
            out[e] = _sha256(q)
    return out


def recheck_done_stages(job_id: str) -> "List[str]":
    """A done stage whose evidence changed or vanished since acceptance is no longer done.

    Paper2Agent's rule: a hash proves agreement with the recorded bytes, so
    changed inputs require re-verification rather than refreshed hashes.
    """
    plan = load_plan(job_id)
    stale = []
    for s in plan.get("stages") or []:
        rec = (s.get("verified") or {}).get("sha256") or {}
        if s["status"] != "done" or not rec:
            continue
        now = evidence_hashes(list(rec))
        changed = [p for p, h in rec.items() if now.get(p) != h]
        if changed:
            s["status"] = "failed"
            s["verified"] = {"passed": False, "reason": "evidence changed since verification: "
                             + ", ".join(changed[:4]), "at": _iso()}
            stale.append(s["id"])
    if stale:
        save_plan(job_id, plan)
        event(job_id, "stale", "re-verification needed: " + ", ".join(stale))
    return stale


def plan_set(job_id: str, stages: "List[dict]") -> dict:
    plan = load_plan(job_id)
    old = {s["id"]: s for s in plan.get("stages") or []}
    out = []
    for i, s in enumerate(stages):
        sid = str(s.get("id") or f"s{i + 1}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,40}", sid):
            raise ValueError(f"bad stage id {sid!r}")
        prev = old.get(sid, {})
        out.append({"id": sid, "title": str(s.get("title") or sid)[:200], "goal": str(s.get("goal") or "")[:1000],
                    "success": str(s.get("success") or "")[:1000],
                    "evidence_expected": [str(x) for x in (s.get("evidence") or [])][:20],
                    "check": s.get("check") if isinstance(s.get("check"), dict) else None,
                    "status": prev.get("status", "pending"), "attempts": prev.get("attempts", 0),
                    "evidence": prev.get("evidence", []), "notes": prev.get("notes", []),
                    "verified": prev.get("verified")})
    plan["stages"] = out
    plan.setdefault("history", []).append({"t": _iso(), "action": "plan_set", "n": len(out)})
    save_plan(job_id, plan)
    event(job_id, "plan", f"plan set: {len(out)} stages", stages=[s["id"] for s in out])
    return plan


def plan_update(job_id: str, stage_id: str, status: str, evidence: "Optional[List[str]]" = None,
                note: str = "") -> "Tuple[dict, str]":
    if status not in STAGE_STATES:
        raise ValueError(f"status must be one of {STAGE_STATES}")
    plan = load_plan(job_id)
    st = next((s for s in plan.get("stages") or [] if s["id"] == stage_id), None)
    if st is None:
        raise ValueError(f"no stage {stage_id!r}; stages: {[s['id'] for s in plan.get('stages') or []]}")
    if evidence:
        st["evidence"] = sorted(set((st.get("evidence") or []) + list(evidence)))
    if note:
        st.setdefault("notes", []).append(f"{_iso()} {note[:500]}")
    msg = ""
    if status == "done":
        ok, reason = run_check(st.get("check"), st.get("evidence") or [])
        st["verified"] = {"passed": ok, "reason": reason, "at": _iso(),
                          "sha256": evidence_hashes(st.get("evidence") or []) if ok else {}}
        if ok:
            st["status"] = "done"
            msg = f"stage {stage_id} verified: {reason}"
        else:
            st["status"] = "failed"
            st["attempts"] = int(st.get("attempts") or 0) + 1
            msg = (f"stage {stage_id} NOT accepted: {reason}. Attempt {st['attempts']} of "
                   f"{MAX_STAGE_ATTEMPTS}; fix it, or mark it blocked with a reason.")
    elif status == "failed":
        st["status"] = "failed"
        st["attempts"] = int(st.get("attempts") or 0) + 1
        msg = f"stage {stage_id} failed (attempt {st['attempts']})"
    else:
        if status == "blocked" and not note:
            raise ValueError("a blocked stage needs a note giving the reason")
        st["status"] = status
        msg = f"stage {stage_id} -> {status}"
    plan.setdefault("history", []).append({"t": _iso(), "action": "update", "stage": stage_id, "status": st["status"]})
    save_plan(job_id, plan)
    event(job_id, "stage", msg, stage=stage_id, status=st["status"])
    return plan, msg


def plan_summary(plan: dict) -> str:
    rows = []
    for s in plan.get("stages") or []:
        mark = {"done": "[x]", "blocked": "[!]", "failed": "[-]", "running": "[>]"}.get(s["status"], "[ ]")
        extra = ""
        if s.get("verified") and not s["verified"].get("passed"):
            extra = f"  <- {s['verified'].get('reason')}"
        if s["status"] == "blocked" and s.get("notes"):
            extra = f"  <- blocked: {s['notes'][-1][21:]}"
        rows.append(f"{mark} {s['id']}: {s['title']} (attempts {s.get('attempts', 0)}){extra}")
    return "\n".join(rows) or "(no plan yet)"


def pending_stages(plan: dict) -> "List[dict]":
    return [s for s in plan.get("stages") or [] if s["status"] not in ("done", "blocked")]


# ─── round runners ──────────────────────────────────────────────────────────

class RoundResult:
    def __init__(self, answer: str = "", stop_reason: str = "", artefacts: "Optional[List[str]]" = None,
                 failed: "Optional[List[dict]]" = None, error: str = "", cost_usd: float = 0.0,
                 session_id: str = ""):
        self.answer, self.stop_reason = answer, stop_reason
        self.artefacts, self.failed = artefacts or [], failed or []
        self.error, self.cost_usd, self.session_id = error, cost_usd, session_id


def _agent_mod():
    try:
        from igvfagent import _agent  # type: ignore
    except Exception:
        import _agent  # type: ignore
    return _agent


class InternalRunner:
    """IGVFagent's own loop, one full run per round, with the job protocol."""

    def __init__(self, backend: Optional[str], model: Optional[str], iterations: int = ROUND_ITERATIONS):
        self.backend, self.model, self.iterations = backend, model, iterations

    def run_round(self, job: dict, prompt: str, history: "List[dict]") -> RoundResult:
        agent = _agent_mod()
        jid = job["id"]

        def cb(ev):
            k, p = ev.kind, ev.payload or {}
            if k == "tool_call_start":
                event(jid, "tool", p.get("name", ""), args=str(p.get("arguments"))[:300])
            elif k == "tool_call_end":
                event(jid, "tool_end", p.get("name", ""), exit=p.get("exit_code"))
            elif k == "llm_call_start":
                event(jid, "think", f"iteration {p.get('iteration')}/{p.get('max_iterations')}")

        try:
            res = agent.run(prompt, backend=self.backend, model=self.model, max_iterations=self.iterations,
                            max_tokens=int(os.environ.get("IGVF_JOB_MAX_TOKENS", "8192")),
                            system_prompt=agent.DEFAULT_SYSTEM_PROMPT + "\n\n" + JOB_PROTOCOL,
                            history=history, callback=cb, enable_router=False)
        except Exception as e:  # noqa: BLE001
            return RoundResult(error=f"{type(e).__name__}: {e}")
        return RoundResult(answer=res.final_answer or "", stop_reason=res.stop_reason,
                           artefacts=list(res.artefacts or []), failed=list(res.failed_calls or []))


def _mcp_config(job: dict) -> Path:
    cfg = job_dir(job["id"]) / "mcp.json"
    igv = shutil.which("igvfagent")
    cmd, args = (igv, ["mcp", "serve"]) if igv else (sys.executable, ["-m", "igvfagent.cli", "mcp", "serve"])
    env = {"IGVF_JOB_ID": job["id"], "IGVF_PROJECT_ROOT": str(ROOT), "IGVF_LLM_MAX_TOOLS": "2000",
           "IGVF_ACTING_USER": job.get("owner") or "", "IGVF_ACTING_ADMIN": "1" if job.get("owner_admin") else "0"}
    _write_json(cfg, {"mcpServers": {"igvfagent": {"command": cmd, "args": args, "env": env}}})
    return cfg


def claude_allowed_tools(allow_shell: bool) -> "List[str]":
    tools = ["mcp__igvfagent", "Read", "Grep", "Glob", "TodoWrite", "Task"]
    if allow_shell:
        tools += ["Bash", "Write", "Edit"]
    return tools


class ClaudeCodeRunner:
    """Claude Code headless as the orchestrator; IGVFagent tools over MCP."""

    def __init__(self, model: Optional[str], max_turns: int = 80, allow_shell: bool = False,
                 binary: Optional[str] = None):
        self.model, self.max_turns, self.allow_shell = model, max_turns, allow_shell
        self.binary = binary or shutil.which("claude") or "claude"

    def argv(self, job: dict, prompt: str) -> "List[str]":
        cfg = _mcp_config(job)
        a = [self.binary, "-p", prompt, "--output-format", "stream-json", "--verbose",
             "--mcp-config", str(cfg), "--strict-mcp-config",
             "--allowedTools", ",".join(claude_allowed_tools(self.allow_shell)),
             "--max-turns", str(self.max_turns), "--append-system-prompt", JOB_PROTOCOL,
             "--permission-mode", "dontAsk" if not self.allow_shell else "acceptEdits"]
        if not self.allow_shell:
            a += ["--disallowedTools", "Bash,Write,Edit,NotebookEdit,WebFetch"]
        if self.model:
            a += ["--model", self.model]
        if job.get("session_id"):
            a += ["--resume", job["session_id"]]
        return a

    @staticmethod
    def parse_stream(lines, on_event: Callable[[str, str, dict], None]) -> RoundResult:
        res = RoundResult()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            try:
                m = json.loads(line)
            except ValueError:
                continue
            t = m.get("type")
            if t == "system" and m.get("subtype") == "init":
                res.session_id = m.get("session_id") or res.session_id
            elif t == "assistant":
                for c in (m.get("message") or {}).get("content") or []:
                    if c.get("type") == "tool_use":
                        on_event("tool", str(c.get("name", "")).replace("mcp__igvfagent__", ""),
                                 {"args": json.dumps(c.get("input"))[:300]})
                    elif c.get("type") == "text" and c.get("text", "").strip():
                        on_event("say", c["text"][:600], {})
            elif t == "result":
                res.answer = m.get("result") or ""
                res.stop_reason = m.get("subtype") or ""
                res.cost_usd = float(m.get("total_cost_usd") or 0.0)
                res.session_id = m.get("session_id") or res.session_id
                # A round that used all its turns is a round boundary, not an
                # LLM failure: the plan carries the work into the next round.
                # Counted as an error, three long rounds in a row failed a job
                # that was making progress (J202609250420500a7a6e).
                if m.get("is_error") and res.stop_reason != "error_max_turns":
                    res.error = res.answer or res.stop_reason
        return res

    def run_round(self, job: dict, prompt: str, history: "List[dict]") -> RoundResult:
        res = self._run_once(job, prompt)
        if job.get("session_id") and res.error and not res.answer.strip() and (
                res.stop_reason == "error_during_execution" or "conversation" in (res.error or "").lower()):
            # The Claude Code session could not be resumed (it lived in a
            # container a redeploy replaced, or expired). The round prompt
            # carries the task, the plan and the open issues, so a fresh
            # session continues from the plan instead of failing the job.
            event(job["id"], "session", f"session {job['session_id'][:8]} could not be resumed; starting a fresh "
                                        "session from the plan")
            update_job(job["id"], session_id="")
            res = self._run_once({**job, "session_id": ""}, prompt)
        return res

    def _run_once(self, job: dict, prompt: str) -> RoundResult:
        jid = job["id"]
        env = dict(os.environ)
        # Keep Claude Code's sessions with the job on the persistent volume, so
        # --resume still works after the container is recreated.
        cdir = job_dir(jid) / "claude"
        cdir.mkdir(parents=True, exist_ok=True)
        env.setdefault("CLAUDE_CONFIG_DIR", str(cdir))
        env.setdefault("IGVF_PROJECT_ROOT", str(ROOT))
        # job_wait may legitimately wait hours; the IGVFagent MCP server bounds
        # every other call itself (IGVF_MCP_TOOL_TIMEOUT), so the client's
        # own per-call timeout (ms) only needs to cover the longest wait.
        env.setdefault("MCP_TOOL_TIMEOUT", str(int((WAIT_MAX_MIN + 10) * 60 * 1000)))
        if not env.get("ANTHROPIC_API_KEY"):
            return RoundResult(error="ANTHROPIC_API_KEY is not set; the Claude Code orchestrator runs on the API key")
        try:
            proc = subprocess.Popen(self.argv(job, prompt), cwd=str(ROOT), env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True, bufsize=1)
        except OSError as e:
            return RoundResult(error=f"cannot start claude: {e}")
        update_job(jid, child_pid=proc.pid)
        res = self.parse_stream(iter(proc.stdout.readline, ""),
                                lambda k, text, p: event(jid, k, text, **p))
        proc.wait()
        if proc.returncode and not res.answer and not res.error:
            res.error = f"claude exited {proc.returncode}"
        return res


# ─── independent verifier ───────────────────────────────────────────────────

VERIFIER_PROMPT = """\
You are an independent verifier. You did not do this work and share none of
its context. Judge only from the evidence below whether the job's final
report is supported: every claimed result must correspond to a done stage
whose evidence exists, numbers must match the evidence excerpts, and nothing
failed or blocked may be presented as done. Missing work that the plan
marked blocked with a real reason is acceptable if the report says so.
Reply with ONLY a JSON object: {"verdict": "pass" | "fail", "issues": ["..."]}.
"""


def _evidence_digest(plan: dict, limit_chars: int = int(os.environ.get("IGVF_VERIFIER_EVIDENCE_CHARS", "60000"))) -> str:
    """What the verifier sees: EVERY evidence path of every stage (with its
    size, or "missing"), then excerpts within a budget, most informative
    first (small JSON summaries, then reports and logs), a fair share per
    file. Listing only the first few files per stage made later evidence
    invisible and got a correct report rejected as unsupported."""
    out, excerpts = [], []
    for s in plan.get("stages") or []:
        out.append(f"## stage {s['id']} [{s['status']}] {s['title']}\nsuccess: {s.get('success')}\n"
                   f"harness: {(s.get('verified') or {}).get('reason', '-')}")
        for e in s.get("evidence") or []:
            q = _safe_path(e)
            if not q or not q.is_file():
                out.append(f"- {e}: (missing)")
                continue
            out.append(f"- {e} ({q.stat().st_size:,} bytes)")
            if q.suffix.lower() in (".md", ".json", ".tsv", ".csv", ".txt", ".log"):
                rank = (0 if q.suffix.lower() == ".json" and q.stat().st_size < 200_000 else
                        1 if q.suffix.lower() == ".md" else 2)
                excerpts.append((rank, e, q))
    listing = "\n".join(out)
    budget = max(limit_chars - len(listing), 4000)
    excerpts.sort(key=lambda x: x[0])
    seen: set = set()
    uniq = [x for x in excerpts if not (x[1] in seen or seen.add(x[1]))]
    per = max(budget // max(len(uniq), 1), 600)
    parts, used = [], 0
    for _, e, q in uniq:
        if used >= budget:
            parts.append(f"(further excerpts omitted: budget; {len(uniq) - len(parts)} file(s) listed above)")
            break
        try:
            text = q.read_text(errors="replace")
        except OSError:
            continue
        take = min(per, budget - used, len(text))
        head = text[:take] + ("\n…" if len(text) > take else "")
        parts.append(f"### {e}\n```\n{head}\n```")
        used += len(head) + len(e) + 20
    return listing + "\n\n# Evidence excerpts\n" + "\n".join(parts)


_REPRO = re.compile(r"\b(reproduc\w*|replicat\w*|re-?run the (?:paper|analysis))\b", re.I)


def code_route_gap(job: dict, plan: dict) -> Optional[str]:
    """A paper reproduction whose harvested Code Availability names a GitHub
    repository must have run that code (a paper-code summary.json in the
    evidence) or blocked the stage with a reason. Reproducing from deposited
    tables alone is a cross-check, not the reproduction."""
    if not _REPRO.search(job.get("query") or ""):
        return None
    ev = [e for st in plan.get("stages") or [] for e in (st.get("evidence") or [])]
    if any("PaperCode" in e and e.endswith(("summary.json", "done.json", "report.md")) for e in ev):
        return None
    if any(st["status"] == "blocked" and re.search(r"code|repositor|paper.code", " ".join(
            [st.get("title") or "", st.get("id") or ""] + list(st.get("notes") or [])), re.I)
           for st in plan.get("stages") or []):
        return None
    repos: "List[str]" = []
    for e in ev:
        if e.endswith("harvest.json"):
            q = _safe_path(e)
            hv = _read_json(q) if q else None
            for hit in ((hv or {}).get("accessions") or {}).get("github_repo") or []:
                if hit.get("in_data_availability"):
                    repos.append(str(hit.get("value", "")).rstrip(".,;"))
    if not repos:
        return None
    return (f"The paper's Code Availability names {', '.join(sorted(set(repos))[:3])}, but no stage ran the authors' "
            "code. Run paper_code_reproduce on it (job_wait on <run_dir>/done.json) and report its summary.json, or "
            "mark a code stage blocked with the reason.")


VERIFIER_FALLBACK = [m.strip() for m in os.environ.get("IGVF_VERIFIER_FALLBACK",
                                                        "claude-opus-5-5,claude-sonnet-5").split(",") if m.strip()]


def _parse_verdict(text: str) -> Optional[dict]:
    """The verifier's JSON, or None when there is no usable verdict (empty
    reply, refusal, prose, JSON cut off at the token limit)."""
    for m in re.finditer(r"\{", text or ""):
        depth = 0
        for j in range(m.start(), len(text)):
            depth += {"{": 1, "}": -1}.get(text[j], 0)
            if depth == 0:
                try:
                    v = json.loads(text[m.start():j + 1])
                except ValueError:
                    break
                if isinstance(v, dict) and str(v.get("verdict", "")).lower() in ("pass", "fail"):
                    return v
                break
    return None


def verify(job: dict, plan: dict, answer: str, llm_call: "Optional[Callable]" = None) -> dict:
    """An independent fresh model call judges the report against the evidence.
    'pass' / 'fail' only come from a real verdict. A model that refuses (a
    safety classifier can trip on genetics evidence), answers nothing or is
    cut off does not count as 'fail': the next model in IGVF_VERIFIER_FALLBACK
    is asked, and if none gives a verdict the result is 'error', never a fail
    that sends the job into repair rounds nobody can satisfy."""
    user = (f"# Task\n{job['query']}\n\n# Plan and evidence\n{_evidence_digest(plan)}\n\n"
            f"# Final report under review\n{answer[:8000]}")
    if llm_call is not None:
        v = _parse_verdict(llm_call(VERIFIER_PROMPT, user) or "")
        if v is None:
            return {"verdict": "error", "issues": ["verifier gave no usable verdict"]}
        return {"verdict": str(v["verdict"]).lower(), "issues": [str(i) for i in (v.get("issues") or [])][:12]}
    try:
        from igvfagent import _llm  # type: ignore
    except Exception:
        import _llm  # type: ignore
    backend = job.get("backend") if job.get("orchestrator") == "internal" else "anthropic"
    first = job.get("verifier_model") or job.get("model") or None
    models = [first] + [m for m in VERIFIER_FALLBACK if m != first] if backend in (None, "anthropic") else [first]
    tried = []
    for model in models:
        try:
            msg = _llm.chat([{"role": "system", "content": VERIFIER_PROMPT}, {"role": "user", "content": user}],
                            backend=backend or None, model=model, max_tokens=4000, temperature=0.0)
            text = msg.content if isinstance(msg.content, str) else json.dumps(msg.content)
            why = getattr(msg, "stop_reason", None) or getattr(msg, "finish_reason", None) or ""
        except Exception as e:  # noqa: BLE001
            tried.append(f"{model or 'default'}: {type(e).__name__}: {str(e)[:120]}")
            continue
        v = _parse_verdict(text)
        if v is not None:
            out = {"verdict": str(v["verdict"]).lower(), "issues": [str(i) for i in (v.get("issues") or [])][:12],
                   "model": model}
            if tried:
                out["fallback_from"] = tried
            return out
        tried.append(f"{model or 'default'}: " + ("refused" if "refus" in str(why).lower() else
                                                  "empty reply" if not (text or "").strip() else "no JSON verdict"))
    # an unavailable verifier must not pass work silently, nor fail it
    return {"verdict": "error", "issues": ["independent verifier unavailable: " + "; ".join(tried)]}


# ─── the worker: rounds until verifiably done ───────────────────────────────

def _continue_prompt(job: dict, plan: dict, issues: "List[str]", last_answer: str, round_no: int) -> str:
    pend = pending_stages(plan)
    parts = [f"JOB CONTINUES — round {round_no}. Original task:\n{job['query']}\n",
             "Plan status (harness-verified):\n" + plan_summary(plan)]
    if issues:
        parts.append("Must address before finishing:\n" + "\n".join(f"- {i}" for i in issues[:12]))
    over = [s["id"] for s in pend if int(s.get("attempts") or 0) >= MAX_STAGE_ATTEMPTS]
    if over:
        parts.append(f"Stages {', '.join(over)} have used {MAX_STAGE_ATTEMPTS} attempts: fix them now with a "
                     "different approach or mark them blocked with the concrete reason.")
    if not plan.get("stages"):
        parts.append("There is no plan yet: call plan_set first.")
    elif pend:
        parts.append(f"Continue with the next pending stage: {pend[0]['id']} — {pend[0]['title']}. "
                     "Do not redo stages marked [x].")
    else:
        parts.append("All stages are done or blocked: write the final report now.")
    if last_answer:
        parts.append("Your previous round ended with:\n" + last_answer[-2500:])
    return "\n\n".join(parts)


def _heartbeat(job_id: str, stop: threading.Event) -> None:
    while not stop.wait(HEARTBEAT_S):
        try:
            update_job(job_id, heartbeat=_now())
        except Exception:  # noqa: BLE001
            pass


def make_runner(job: dict):
    if job.get("orchestrator") == "claude_code":
        return ClaudeCodeRunner(job.get("model"), allow_shell=bool(job.get("allow_shell")))
    return InternalRunner(job.get("backend"), job.get("model"), int(job.get("round_iterations") or ROUND_ITERATIONS))


def run_worker(job_id: str, runner=None, verifier: "Optional[Callable]" = None, sleep=time.sleep) -> dict:
    job = load_job(job_id)
    os.environ["IGVF_JOB_ID"] = job_id
    if job.get("owner"):
        os.environ["IGVF_ACTING_USER"] = job["owner"]
    os.environ["IGVF_ACTING_ADMIN"] = "1" if job.get("owner_admin") else "0"
    try:
        from igvfagent import _history as _h  # type: ignore
    except Exception:
        try:
            import _history as _h  # type: ignore
        except Exception:
            _h = None
    if _h is not None and job.get("owner"):
        try:
            _h.set_actor(job["owner"])
        except Exception:  # noqa: BLE001
            pass
    # wait for a slot
    while sum(1 for j in list_jobs() if j["effective_status"] == "running" and j["id"] != job_id) >= MAX_RUNNING:
        update_job(job_id, status="queued", heartbeat=_now(), pid=os.getpid(), host=socket.gethostname())
        sleep(30)
    if job.get("orchestrator") == "internal" and (job.get("backend") or os.environ.get("IGVF_LLM_BACKEND") or
                                                   "anthropic") == "anthropic":
        # Anthropic has no 128-function limit: a job should see the whole
        # registry, not the parity subset chat uses for OpenAI-family models
        os.environ["IGVF_LLM_MAX_TOOLS"] = os.environ.get("IGVF_JOB_MAX_TOOLS", "700")
    runner = runner or make_runner(job)
    job = update_job(job_id, status="running", pid=os.getpid(), host=socket.gethostname(), heartbeat=_now(),
                     started_at=job.get("started_at") or _iso())
    stop_hb = threading.Event()
    hb = threading.Thread(target=_heartbeat, args=(job_id, stop_hb), daemon=True)
    hb.start()
    deadline = _now() + float(job.get("budget_minutes") or DEFAULT_BUDGET_MIN) * 60
    event(job_id, "start", f"worker {os.getpid()} ({job.get('orchestrator')}), budget "
                           f"{job.get('budget_minutes')} min / {job.get('max_rounds')} rounds")
    history: "List[dict]" = job.get("history") or []
    issues: "List[str]" = job.get("pending_issues") or []
    stale = recheck_done_stages(job_id)
    if stale:
        issues = issues + [f"stage {s} must be re-verified: its evidence changed since it was accepted" for s in stale]
    last_answer = job.get("last_answer") or ""
    llm_errors = 0
    try:
        while True:
            job = load_job(job_id)
            if (job_dir(job_id) / "STOP").exists():
                job = update_job(job_id, status="stopped")
                event(job_id, "stopped", "stop requested")
                break
            rounds = int(job.get("rounds") or 0)
            if rounds >= int(job.get("max_rounds") or DEFAULT_MAX_ROUNDS) or _now() > deadline:
                job = update_job(job_id, status="budget_exhausted", pending_issues=issues, last_answer=last_answer)
                event(job_id, "budget", f"budget exhausted after {rounds} rounds; resume to continue")
                break
            plan = load_plan(job_id)
            prompt = job["query"] if rounds == 0 and not plan.get("stages") else \
                _continue_prompt(job, plan, issues, last_answer, rounds + 1)
            event(job_id, "round", f"round {rounds + 1} starts")
            res = runner.run_round(job, prompt, history[-8:])
            rounds += 1
            rdir = job_dir(job_id) / "rounds"
            _write_json(rdir / f"{rounds:02d}.json", {"prompt": prompt, "answer": res.answer, "stop": res.stop_reason,
                                                      "failed": res.failed, "error": res.error,
                                                      "artefacts": res.artefacts, "cost_usd": res.cost_usd})
            upd: Dict[str, Any] = {"rounds": rounds, "cost_usd": round(float(job.get("cost_usd") or 0) + res.cost_usd, 4)}
            if res.session_id:
                upd["session_id"] = res.session_id
            upd["artefacts"] = sorted(set((job.get("artefacts") or []) + res.artefacts))[-200:]
            job = update_job(job_id, **upd)
            if res.error:
                llm_errors += 1
                event(job_id, "error", res.error)
                if llm_errors > 2:
                    job = update_job(job_id, status="failed", error=res.error)
                    break
                sleep(min(60, 10 * llm_errors))
                continue
            llm_errors = 0
            max_usd = float(job.get("max_usd") or 0)
            if max_usd and float(job.get("cost_usd") or 0) >= max_usd:
                job = update_job(job_id, status="budget_exhausted", last_answer=res.answer)
                event(job_id, "budget", f"cost budget ${max_usd:.2f} reached")
                break
            history = (history + [{"role": "user", "content": prompt[:4000]},
                                  {"role": "assistant", "content": res.answer[-6000:]}])[-12:]
            last_answer = res.answer
            plan = load_plan(job_id)
            issues = [f"tool {f.get('name')} failed: {str(f.get('detail'))[:200]}" for f in res.failed[:6]]
            if pending_stages(plan) or (not plan.get("stages") and rounds == 1 and len(job["query"]) > 200):
                job = update_job(job_id, history=history, pending_issues=issues, last_answer=last_answer)
                continue
            vc = int(job.get("verify_cycles") or 0)
            verdict = verify(job, plan, res.answer, llm_call=verifier)
            gap = code_route_gap(job, plan)
            if gap:  # deterministic: the paper names its code and nobody ran it
                verdict = {"verdict": "fail", "issues": [gap] + list(verdict.get("issues") or [])}
            event(job_id, "verify", f"verifier: {verdict['verdict']}", issues=verdict["issues"])
            if verdict["verdict"] == "fail" and verdict["issues"] and vc < MAX_VERIFY_CYCLES:
                # only a real verdict with reasons starts a repair round; an
                # unavailable verifier ("error") is recorded, not "fixed"
                issues = [f"verifier: {i}" for i in verdict["issues"]]
                job = update_job(job_id, verify_cycles=vc + 1, history=history, pending_issues=issues,
                                 last_answer=last_answer, last_verdict=verdict)
                continue
            blocked = any(s["status"] == "blocked" for s in plan.get("stages") or [])
            status = "done_with_blocked" if blocked else "done"
            if verdict["verdict"] == "fail":
                status = "done_with_blocked"
            # "error": the stages are harness-verified but no independent model
            # would judge the report; the job is done, marked unverified.
            (job_dir(job_id) / "answer.md").write_text(
                f"# {job['query'][:200]}\n\n{res.answer}\n\n## Plan\n\n```\n{plan_summary(plan)}\n```\n\n"
                f"Verifier: **{verdict['verdict']}**" + (": " + "; ".join(verdict["issues"]) if verdict["issues"] else "")
                + "\n")
            job = update_job(job_id, status=status, finished_at=_iso(), last_answer=res.answer, last_verdict=verdict,
                             history=history)
            event(job_id, "done", f"{status} after {rounds} rounds")
            break
    finally:
        stop_hb.set()
        _record_reproduction(job_id)
    return load_job(job_id)


def _record_reproduction(job_id: str) -> None:
    """Keep a durable per-paper record of every finished reproduction or
    benchmark attempt (Data/Reproductions/, 📑 Reproductions in the UI)."""
    try:
        j = load_job(job_id)
        if effective_status(j) in ("running", "queued"):
            return
        try:
            from igvfagent import reproductions as rp  # type: ignore
        except Exception:
            import reproductions as rp  # type: ignore
        plan_text = json.dumps(load_plan(job_id))
        if rp.REPRO_QUERY.search(j.get("query") or "") or any(
                k in plan_text for k in ("PaperCode", "concordance", "harvest.json")):
            rec = rp.record_from_job(job_id)
            if rec:
                event(job_id, "record", f"reproduction record {rec['paper_id']}: {rec['outcome']}")
    except Exception as e:  # noqa: BLE001  (a record must never fail the job)
        try:
            event(job_id, "record", f"reproduction record not written: {type(e).__name__}: {e}")
        except Exception:  # noqa: BLE001
            pass


# ─── starting, stopping, resuming ───────────────────────────────────────────

def _worker_argv(job_id: str) -> "List[str]":
    igv = shutil.which("igvfagent")
    return [igv, "job", "worker", job_id] if igv else [sys.executable, "-m", "igvfagent.cli", "job", "worker", job_id]


def spawn_worker(job_id: str) -> int:
    d = job_dir(job_id)
    try:
        (d / "STOP").unlink()
    except OSError:
        pass
    env = dict(os.environ)
    env["IGVF_JOB_ID"] = job_id
    env.setdefault("IGVF_PROJECT_ROOT", str(ROOT))
    logf = (d / "worker.log").open("a")
    p = subprocess.Popen(_worker_argv(job_id), cwd=str(ROOT), env=env, stdout=logf, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True)
    update_job(job_id, pid=p.pid, heartbeat=_now())
    return p.pid


def create_job(query: str, *, owner: str = "", owner_admin: bool = False, orchestrator: str = "internal",
               backend: Optional[str] = None, model: Optional[str] = None,
               budget_minutes: float = DEFAULT_BUDGET_MIN, max_rounds: int = DEFAULT_MAX_ROUNDS,
               max_usd: float = 0.0, allow_shell: bool = False, title: str = "") -> dict:
    if orchestrator not in ("internal", "claude_code"):
        raise ValueError("orchestrator must be internal or claude_code")
    jid = "J" + time.strftime("%Y%m%d%H%M%S") + uuid.uuid4().hex[:6]
    job = {"id": jid, "query": query, "title": title or query.strip().splitlines()[0][:100], "owner": owner,
           "owner_admin": owner_admin, "orchestrator": orchestrator, "backend": backend, "model": model,
           "budget_minutes": budget_minutes, "max_rounds": max_rounds, "max_usd": max_usd,
           "allow_shell": allow_shell, "round_iterations": ROUND_ITERATIONS, "status": "queued",
           "created_at": _iso(), "rounds": 0, "cost_usd": 0.0}
    job_dir(jid).mkdir(parents=True, exist_ok=True)
    proj = _file_into_active_project(job)
    if proj:
        job["project"], job["project_name"] = proj
    save_job(job)
    save_plan(jid, {"stages": [], "history": []})
    event(jid, "created", f"job created by {owner or 'local'}: {job['title']}"
          + (f" (project {job.get('project_name')})" if proj else ""))
    return job


def _file_into_active_project(job: dict) -> "Optional[Tuple[str, str]]":
    """File a new job into its owner's active project, like a chat answer."""
    try:
        try:
            from igvfagent import _history as H  # type: ignore
        except Exception:
            import _history as H  # type: ignore
        viewer = job.get("owner") or None
        act = H.active_project(viewer=viewer)
        if not act:
            return None
        res = H.add_item(act["id"], "job", job["id"], title=job.get("title") or "", owner=job.get("owner") or "",
                         viewer=viewer)
        return (act["id"], act["name"]) if res.get("ok") else None
    except Exception:  # noqa: BLE001  (no history store: jobs still work)
        return None


def start_job(query: str, **kw) -> dict:
    job = create_job(query, **kw)
    spawn_worker(job["id"])
    return load_job(job["id"])


def stop_job(job_id: str, now: bool = False) -> dict:
    (job_dir(job_id) / "STOP").write_text(_iso())
    job = load_job(job_id)
    if now:
        for pid in (job.get("child_pid"), job.get("pid")):
            if _pid_alive(pid):
                try:
                    os.killpg(os.getpgid(int(pid)), signal.SIGTERM)
                except (OSError, ValueError):
                    pass
        job = update_job(job_id, status="stopped")
    event(job_id, "stop", "stop requested" + (" (immediate)" if now else " (after this round)"))
    return job


def resume_job(job_id: str, extra_minutes: float = 0.0, extra_rounds: int = 0, note: str = "") -> dict:
    job = load_job(job_id)
    st = effective_status(job)
    if st in ("running", "queued"):
        return job
    upd: Dict[str, Any] = {"status": "queued"}
    if extra_rounds or st == "budget_exhausted":
        upd["max_rounds"] = int(job.get("rounds") or 0) + (extra_rounds or DEFAULT_MAX_ROUNDS // 2)
    upd["budget_minutes"] = float(job.get("budget_minutes") or DEFAULT_BUDGET_MIN) + (extra_minutes or DEFAULT_BUDGET_MIN)
    if note:
        upd["pending_issues"] = (job.get("pending_issues") or []) + [f"user: {note}"]
    update_job(job_id, **upd)
    try:
        (job_dir(job_id) / "STOP").unlink()
    except OSError:
        pass
    event(job_id, "resume", "resumed" + (f": {note}" if note else ""))
    spawn_worker(job_id)
    return load_job(job_id)


# ─── tools the agent calls inside a job ─────────────────────────────────────

def _current_job() -> str:
    jid = os.environ.get("IGVF_JOB_ID", "")
    if not jid:
        raise SystemExit("plan tools work inside a long-running job (start one with `igvfagent job start` "
                         "or the job_start tool); this call is not part of a job.")
    return jid


def delegate(tasks: "List[str]", backend: Optional[str] = None, model: Optional[str] = None,
             iterations: int = 25, runner: Optional[Callable[[str], dict]] = None) -> "List[dict]":
    """Run focused sub-agents in parallel, each in a fresh context."""
    if int(os.environ.get("IGVF_JOB_DEPTH", "0")) >= 1:
        raise SystemExit("delegated sub-agents cannot delegate further")
    tasks = [t for t in tasks if t.strip()][:MAX_DELEGATES]

    def one(task: str) -> dict:
        if runner:
            return runner(task)
        agent = _agent_mod()
        os.environ["IGVF_JOB_DEPTH"] = "1"
        r = agent.run(task, backend=backend, model=model, max_iterations=iterations, callback=None,
                      system_prompt=agent.DEFAULT_SYSTEM_PROMPT + "\n\nYou are a delegated sub-agent: do exactly "
                      "the task, then reply with a concise result, the numbers, and the file paths you produced.")
        return {"task": task, "answer": r.final_answer, "artefacts": r.artefacts, "stop": r.stop_reason,
                "failed": len(r.failed_calls or [])}

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(tasks) or 1) as ex:
        return list(ex.map(one, tasks))


def wait_for(path: Optional[str] = None, raw_job: Optional[str] = None, timeout_min: float = 240,
             poll_s: float = 30, sleep=time.sleep) -> "Tuple[bool, str]":
    end = _now() + timeout_min * 60
    last_size = -1
    while _now() < end:
        if path:
            q = _safe_path(path)
            if q is None:
                return False, "path is outside the workspace"
            if q.exists():
                size = q.stat().st_size if q.is_file() else sum(1 for _ in q.iterdir())
                if size == last_size:
                    return True, f"{path} is present and stable"
                last_size = size
        if raw_job:
            if not re.fullmatch(r"[A-Za-z0-9_-]+", raw_job):
                return False, "bad job id"
            p = subprocess.run([*_worker_argv("Jxxxxxxx")[:-3], "raw-pipeline", "status", raw_job, "--tail", "5"],
                               capture_output=True, text=True, timeout=120)
            out = p.stdout + p.stderr
            if "No such job" in out:
                return False, out.strip()
            if not re.search(r"state:\s+running", out):
                return True, "raw-pipeline job finished:\n" + out.strip()[-1500:]
        sleep(poll_s)
    return False, f"still not finished after {timeout_min} min"


# ─── CLI ────────────────────────────────────────────────────────────────────

def _print_job(j: dict) -> None:
    print(f"{j['id']}  {j.get('effective_status') or effective_status(j):18s} {j.get('orchestrator'):11s} "
          f"rounds {j.get('rounds', 0)}  ${float(j.get('cost_usd') or 0):.2f}  {j.get('title', '')[:70]}")


def cmd_start(a) -> int:
    if os.environ.get("IGVF_JOB_ID"):
        print("already inside job " + os.environ["IGVF_JOB_ID"] + ": work through its plan, or use "
              "delegate_tasks for parallel sub-tasks, instead of starting a nested job")
        return 2
    j = start_job(a.query, owner=os.environ.get("IGVF_ACTING_USER", ""),
                  owner_admin=os.environ.get("IGVF_ACTING_ADMIN") == "1", orchestrator=a.orchestrator,
                  backend=a.backend, model=a.model, budget_minutes=a.budget_minutes, max_rounds=a.max_rounds,
                  max_usd=a.max_usd, allow_shell=a.allow_shell)
    print(f"Started job {j['id']} ({j['orchestrator']}); it keeps running in the background.")
    print(f"  follow:  igvfagent job status {j['id']}   |   igvfagent job logs {j['id']}")
    print(f"JSON: {job_dir(j['id']) / 'job.json'}")
    return 0


def cmd_list(a) -> int:
    viewer = os.environ.get("IGVF_ACTING_USER") or None
    for j in list_jobs(viewer, os.environ.get("IGVF_ACTING_ADMIN", "1") == "1" or not viewer, a.limit):
        _print_job(j)
    return 0


def cmd_status(a) -> int:
    j = load_job(a.job)
    j["effective_status"] = effective_status(j)
    _print_job(j)
    print("\nPlan:\n" + plan_summary(load_plan(a.job)))
    if j.get("last_verdict"):
        print(f"\nVerifier: {j['last_verdict'].get('verdict')}  {'; '.join(j['last_verdict'].get('issues') or [])}")
    ans = job_dir(a.job) / "answer.md"
    if ans.exists():
        print(f"\nReport: {ans}")
    elif j.get("last_answer"):
        print("\nLatest:\n" + j["last_answer"][-1500:])
    return 0


def cmd_logs(a) -> int:
    for e in read_events(a.job, a.last):
        print(f"{e['t']}  {e['kind']:9s} {e.get('text', '')[:200]}")
    return 0


def cmd_stop(a) -> int:
    stop_job(a.job, now=a.now)
    print(f"stop requested for {a.job}")
    return 0


def cmd_resume(a) -> int:
    j = resume_job(a.job, a.extra_minutes, a.extra_rounds, a.note or "")
    print(f"{a.job}: {j['status']}")
    return 0


def reverify(job_id: str) -> dict:
    """Run the independent verifier again on a finished job's report (e.g.
    after a verifier refusal) and record the verdict; the job is not re-run."""
    job = load_job(job_id)
    plan = load_plan(job_id)
    answer = job.get("last_answer") or ""
    v = verify(job, plan, answer)
    st = job.get("status")
    finished = bool(plan.get("stages")) and all(x["status"] in ("done", "blocked") for x in plan["stages"])
    # a job that ran out of budget after every stage was done is finished too
    if st in ("done", "done_with_blocked") or (st in ("budget_exhausted", "stopped") and finished):
        blocked = any(x["status"] == "blocked" for x in plan.get("stages") or [])
        st = "done_with_blocked" if (blocked or v["verdict"] == "fail") else "done"
    update_job(job_id, last_verdict=v, status=st)
    event(job_id, "verify", f"re-verified: {v['verdict']}", issues=v["issues"])
    ans = job_dir(job_id) / "answer.md"
    if ans.exists():
        t = re.sub(r"\nVerifier: \*\*[^*]+\*\*.*$", "", ans.read_text(), flags=re.S)
        ans.write_text(t.rstrip() + f"\n\nVerifier: **{v['verdict']}**"
                       + (" (" + v["model"] + ")" if v.get("model") else "")
                       + (": " + "; ".join(v["issues"]) if v["issues"] else "") + "\n")
    _record_reproduction(job_id)
    return {**v, "status": st}


def cmd_reverify(a) -> int:
    r = reverify(a.job)
    print(json.dumps(r, indent=2))
    return 0 if r["verdict"] in ("pass", "fail") else 1


def cmd_resume_interrupted(a) -> int:
    """After a container restart: restart every job whose worker vanished mid-run."""
    n = 0
    for j in list_jobs(limit=200):
        if j["effective_status"] == "interrupted":
            update_job(j["id"], status="stopped")
            resume_job(j["id"], extra_minutes=0.0001, note="resumed after the service restarted")
            print(f"resumed {j['id']}  {j.get('title', '')[:60]}")
            n += 1
    print(f"{n} interrupted job(s) resumed")
    return 0


def cmd_worker(a) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    j = run_worker(a.job)
    print(f"{a.job}: {j.get('status')}")
    return 0


def _loads_lenient(text: str) -> Any:
    """JSON, or the Python-literal form a tool call may render a list as."""
    try:
        return json.loads(text)
    except ValueError:
        import ast
        return ast.literal_eval(text)


def cmd_plan_set(a) -> int:
    stages = _loads_lenient(a.stages_json)
    if isinstance(stages, dict):
        stages = stages.get("stages") or []
    plan = plan_set(_current_job(), stages)
    print(f"plan recorded: {len(plan['stages'])} stages\n" + plan_summary(plan))
    return 0


def cmd_plan_update(a) -> int:
    plan, msg = plan_update(_current_job(), a.stage, a.status, a.evidence or [], a.note or "")
    print(msg + "\n\n" + plan_summary(plan))
    return 0 if not msg.startswith(f"stage {a.stage} NOT") else 1


def cmd_plan_show(a) -> int:
    print(plan_summary(load_plan(_current_job())))
    return 0


def cmd_delegate(a) -> int:
    res = delegate(a.task, backend=a.backend, model=a.model)
    for r in res:
        print(f"### {r['task'][:100]}\n{r['answer']}\n")
        for p in r.get("artefacts") or []:
            print(f"Wrote: {p}")
    return 0


def cmd_wait(a) -> int:
    ok, msg = wait_for(a.path, a.raw_job, a.timeout_min)
    print(msg)
    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="igvfagent job", description="Durable, planned, verified agent jobs.")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("start", help="Start a long-running agent job in the background.")
    s.add_argument("query")
    s.add_argument("--orchestrator", choices=["internal", "claude_code"], default="internal")
    s.add_argument("--backend")
    s.add_argument("--model")
    s.add_argument("--budget-minutes", type=float, default=DEFAULT_BUDGET_MIN)
    s.add_argument("--max-rounds", type=int, default=DEFAULT_MAX_ROUNDS)
    s.add_argument("--max-usd", type=float, default=0.0, help="Stop when this API cost is reached (claude_code).")
    s.add_argument("--allow-shell", action="store_true", help="claude_code only: allow Bash/Write/Edit.")
    s.set_defaults(func=cmd_start)
    s = sub.add_parser("list")
    s.add_argument("--limit", type=int, default=20)
    s.set_defaults(func=cmd_list)
    for name, fn in (("status", cmd_status), ("worker", cmd_worker)):
        s = sub.add_parser(name)
        s.add_argument("job")
        s.set_defaults(func=fn)
    s = sub.add_parser("logs")
    s.add_argument("job")
    s.add_argument("--last", type=int, default=40)
    s.set_defaults(func=cmd_logs)
    s = sub.add_parser("stop")
    s.add_argument("job")
    s.add_argument("--now", action="store_true")
    s.set_defaults(func=cmd_stop)
    s = sub.add_parser("resume")
    s.add_argument("job")
    s.add_argument("--extra-minutes", type=float, default=0.0)
    s.add_argument("--extra-rounds", type=int, default=0)
    s.add_argument("--note", help="Guidance added to the job before it continues.")
    s.set_defaults(func=cmd_resume)
    s = sub.add_parser("reverify", help="Run the independent verifier again on a finished job's report.")
    s.add_argument("job")
    s.set_defaults(func=cmd_reverify)
    s = sub.add_parser("resume-interrupted", help="Restart jobs whose worker died (run after a redeploy).")
    s.set_defaults(func=cmd_resume_interrupted)
    s = sub.add_parser("plan-set", help="(inside a job) record the staged plan")
    s.add_argument("--stages-json", required=True)
    s.set_defaults(func=cmd_plan_set)
    s = sub.add_parser("plan-update", help="(inside a job) update a stage; done is verified by the harness")
    s.add_argument("--stage", required=True)
    s.add_argument("--status", required=True, choices=list(STAGE_STATES))
    s.add_argument("--evidence", nargs="+", action="extend", default=[])
    s.add_argument("--note")
    s.set_defaults(func=cmd_plan_update)
    s = sub.add_parser("plan-show")
    s.set_defaults(func=cmd_plan_show)
    s = sub.add_parser("delegate", help="(inside a job) run focused sub-agents in parallel")
    s.add_argument("--task", nargs="+", action="extend", required=True)
    s.add_argument("--backend")
    s.add_argument("--model")
    s.set_defaults(func=cmd_delegate)
    s = sub.add_parser("wait", help="(inside a job) wait for a detached job or an output file")
    s.add_argument("--path")
    s.add_argument("--raw-job")
    s.add_argument("--timeout-min", type=float, default=240)
    s.set_defaults(func=cmd_wait)
    s = sub.add_parser("selftest", help="Offline test with scripted rounds (no LLM, no network).")
    s.set_defaults(func=lambda a: selftest())
    return p


def _load_env_file() -> None:
    """KEY=VALUE lines from <root>/.env for variables not already set (local installs;
    the hosted container sets them in its environment)."""
    f = ROOT / ".env"
    try:
        lines = f.read_text().splitlines()
    except OSError:
        return
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip().removeprefix("export ").strip() if hasattr(str, "removeprefix") else k.strip(), v.strip().strip('"').strip("'")
        if k and k not in os.environ:
            os.environ[k] = v


def main(argv=None) -> int:
    _load_env_file()
    a = build_parser().parse_args(argv)
    return int(a.func(a) or 0)


# ─── selftest ───────────────────────────────────────────────────────────────

def selftest() -> int:
    import tempfile
    global JOBS_DIR, ROOT
    fails = []

    def check(name, cond):
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
        if not cond:
            fails.append(name)

    tmp = Path(tempfile.mkdtemp(prefix="agentjobs_")).resolve()
    old_root, old_jobs = ROOT, JOBS_DIR
    old_env_root = os.environ.get("IGVF_PROJECT_ROOT")
    ROOT, JOBS_DIR = tmp, tmp / "Data" / "Jobs"
    os.environ["IGVF_PROJECT_ROOT"] = str(tmp)      # the path guard contains evidence to this root
    os.environ.pop("IGVF_JOB_ID", None)
    (tmp / "Docs" / "Run").mkdir(parents=True)
    try:
        # 1. a scripted internal-style runner: plans, does stage 1, fakes stage 2, then fixes it
        class Scripted:
            def __init__(self):
                self.calls = 0

            def run_round(self, job, prompt, history):
                self.calls += 1
                jid = job["id"]
                os.environ["IGVF_JOB_ID"] = jid
                if self.calls == 1:
                    plan_set(jid, [{"id": "s1", "title": "download", "evidence": ["Docs/Run/data.tsv"],
                                    "check": {"kind": "rows", "path": "Docs/Run/data.tsv", "min": 2}},
                                   {"id": "s2", "title": "analyse", "evidence": ["Docs/Run/result.json"],
                                    "check": {"kind": "json", "path": "Docs/Run/result.json", "key": "r",
                                              "op": ">=", "value": 0.8}}])
                    (tmp / "Docs/Run/data.tsv").write_text("a\tb\n1\t2\n3\t4\n")
                    plan_update(jid, "s1", "done", ["Docs/Run/data.tsv"])
                    (tmp / "Docs/Run/result.json").write_text('{"r": 0.5}')
                    _p, msg = plan_update(jid, "s2", "done", ["Docs/Run/result.json"])
                    self.rejected = "NOT accepted" in msg
                    return RoundResult(answer="stage 2 done (claimed)")
                if self.calls == 2:
                    self.saw_issue = "r = 0.5" in prompt and "s2" in prompt
                    (tmp / "Docs/Run/result.json").write_text('{"r": 0.91}')
                    plan_update(jid, "s2", "done", ["Docs/Run/result.json"])
                    return RoundResult(answer="r = 0.91, reproduced")
                return RoundResult(answer="final report: r = 0.91 (Docs/Run/result.json)")

        verdicts = iter(['{"verdict":"fail","issues":["report does not cite data.tsv"]}',
                         '{"verdict":"pass","issues":[]}'])
        j = create_job("reproduce the analysis", owner="alice")
        sc = Scripted()
        out = run_worker(j["id"], runner=sc, verifier=lambda s, u: next(verdicts), sleep=lambda s: None)
        check("a done claim with failing evidence is rejected by the harness", getattr(sc, "rejected", False))
        check("the next round is told exactly why (value and stage)", getattr(sc, "saw_issue", False))
        check("verifier fail triggers a repair round; pass then finishes",
              out["status"] == "done" and out["rounds"] == 3 and out["verify_cycles"] == 1)
        check("answer.md written with plan and verdict",
              "Verifier: **pass**" in (job_dir(j["id"]) / "answer.md").read_text())
        ev = [e["kind"] for e in read_events(j["id"], 200)]
        check("events record plan, stages, rounds and verification",
              all(k in ev for k in ("plan", "stage", "round", "verify", "done")))

        (tmp / "Docs/Run/result.json").write_text('{"r": 0.2}')
        check("changed evidence re-opens a verified stage on resume", recheck_done_stages(j["id"]) == ["s2"]
              and next(s for s in load_plan(j["id"])["stages"] if s["id"] == "s2")["status"] == "failed")
        # 2. budget: a runner that never finishes stops at max_rounds and can be resumed
        class Lazy:
            def run_round(self, job, prompt, history):
                os.environ["IGVF_JOB_ID"] = job["id"]
                if not load_plan(job["id"])["stages"]:
                    plan_set(job["id"], [{"id": "s1", "title": "never done"}])
                return RoundResult(answer="still working")
        j2 = create_job("long task", owner="bob", max_rounds=2)
        o2 = run_worker(j2["id"], runner=Lazy(), verifier=lambda s, u: '{"verdict":"pass"}', sleep=lambda s: None)
        check("round budget stops a job as budget_exhausted (resumable)", o2["status"] == "budget_exhausted"
              and o2["rounds"] == 2)
        # 3. blocked stage with reason finishes as done_with_blocked, not done
        class Blocker:
            def run_round(self, job, prompt, history):
                os.environ["IGVF_JOB_ID"] = job["id"]
                plan_set(job["id"], [{"id": "s1", "title": "needs controlled data"}])
                try:
                    plan_update(job["id"], "s1", "blocked")
                    self.allowed_without_reason = True
                except ValueError:
                    self.allowed_without_reason = False
                plan_update(job["id"], "s1", "blocked", note="dbGaP access required")
                return RoundResult(answer="blocked: needs dbGaP")
        b = Blocker()
        j3 = create_job("controlled task", owner="alice")
        o3 = run_worker(j3["id"], runner=b, verifier=lambda s, u: '{"verdict":"pass"}', sleep=lambda s: None)
        check("blocking needs a reason", b.allowed_without_reason is False)
        check("a blocked stage ends as done_with_blocked", o3["status"] == "done_with_blocked")
        # 3b. a verifier that refuses or returns nothing is "unverified", never
        # a "fail" that sends the job into repair rounds nobody can satisfy
        j3b = create_job("verifier refuses", owner="alice")
        class Done1:
            def run_round(self, job, prompt, history):
                os.environ["IGVF_JOB_ID"] = job["id"]
                plan_set(job["id"], [{"id": "s1", "title": "one"}])
                (tmp / "Docs/Run/x.txt").write_text("ok")
                plan_update(job["id"], "s1", "done", ["Docs/Run/x.txt"])
                return RoundResult(answer="final report")
        o3b = run_worker(j3b["id"], runner=Done1(), verifier=lambda s, u: "", sleep=lambda s: None)
        check("a refusing/empty verifier ends the job done (unverified) without repair rounds",
              o3b["status"] == "done" and o3b.get("verify_cycles", 0) == 0 and o3b["rounds"] == 1
              and (o3b.get("last_verdict") or {}).get("verdict") == "error")
        check("a verdict inside prose is read; a cut-off one is not",
              (_parse_verdict('ok: {"verdict": "pass", "issues": []}') or {}).get("verdict") == "pass"
              and _parse_verdict('{"verdict": "fail", "issues": ["cut') is None)
        many = []
        for i in range(10):
            f = tmp / "Docs/Run" / f"e{i}.json"
            f.write_text(json.dumps({"k": i}))
            many.append(f"Docs/Run/e{i}.json")
        dg = _evidence_digest({"stages": [{"id": "r", "status": "done", "title": "t", "evidence": many}]})
        check("the verifier sees every evidence file of a stage, not just the first few",
              all(e in dg for e in many) and '"k": 9' in dg)
        # 4. stop file stops between rounds
        j4 = create_job("stoppable", owner="alice")
        (job_dir(j4["id"]) / "STOP").write_text("x")
        o4 = run_worker(j4["id"], runner=Lazy(), sleep=lambda s: None)
        check("STOP ends the job before the next round", o4["status"] == "stopped" and o4["rounds"] == 0)
        # 5. evidence outside the workspace or secret paths never counts
        ok, why = run_check(None, ["/etc/passwd"])
        check("evidence outside the workspace is rejected", not ok and "missing" in why)
        # 6. owner scoping
        check("a user lists only their own jobs",
              {x["owner"] for x in list_jobs("alice", is_admin=False)} == {"alice"})
        # 7. interrupted detection
        update_job(j2["id"], status="running", pid=999999, heartbeat=_now() - 600)
        check("a running job whose worker vanished shows as interrupted",
              effective_status(load_job(j2["id"])) == "interrupted")
        # 8. Claude Code stream parsing
        lines = [json.dumps({"type": "system", "subtype": "init", "session_id": "S1"}),
                 json.dumps({"type": "assistant", "message": {"content": [
                     {"type": "tool_use", "name": "mcp__igvfagent__plan_set", "input": {"x": 1}},
                     {"type": "text", "text": "planning"}]}}),
                 json.dumps({"type": "result", "subtype": "success", "result": "all done", "session_id": "S1",
                             "total_cost_usd": 0.42, "is_error": False})]
        seen = []
        r = ClaudeCodeRunner.parse_stream(lines, lambda k, t, p: seen.append((k, t)))
        check("claude stream: session id, tool names, cost and result parsed",
              r.session_id == "S1" and r.cost_usd == 0.42 and r.answer == "all done" and ("tool", "plan_set") in seen)
        r = ClaudeCodeRunner.parse_stream([json.dumps({"type": "result", "subtype": "error_max_turns", "result": "",
                                                       "session_id": "S1", "total_cost_usd": 9.9, "is_error": True})],
                                          lambda k, t, p: None)
        check("claude stream: a round that ran out of turns is a round boundary, not an error",
              r.error == "" and r.stop_reason == "error_max_turns" and r.cost_usd == 9.9)
        jcc = create_job("x", orchestrator="claude_code")
        argv = ClaudeCodeRunner(None).argv({**jcc, "session_id": "S1"}, "go")
        check("claude argv: MCP-only tools, no shell by default, session resumed",
              "--strict-mcp-config" in argv and "Bash" not in argv[argv.index("--allowedTools") + 1]
              and argv[argv.index("--resume") + 1] == "S1")
        # 9. delegation runs tasks in parallel with a fresh runner each
        t0 = time.time()
        res = delegate(["a", "b", "c"], runner=lambda t: (time.sleep(0.3), {"task": t, "answer": t.upper()})[1])
        check("delegate_tasks runs sub-agents in parallel", [x["answer"] for x in res] == ["A", "B", "C"]
              and time.time() - t0 < 0.8)
        # 10. wait for a file
        (tmp / "Docs/Run/flag.txt").write_text("ok")
        ok, msg = wait_for(path="Docs/Run/flag.txt", timeout_min=1, poll_s=0, sleep=lambda s: None)
        check("job_wait returns once an output is present and stable", ok)
    finally:
        ROOT, JOBS_DIR = old_root, old_jobs
        os.environ.pop("IGVF_JOB_ID", None)
        if old_env_root is None:
            os.environ.pop("IGVF_PROJECT_ROOT", None)
        else:
            os.environ["IGVF_PROJECT_ROOT"] = old_env_root
        shutil.rmtree(tmp, ignore_errors=True)
    print("selftest: all checks pass" if not fails else f"selftest: FAILED ({len(fails)})")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
