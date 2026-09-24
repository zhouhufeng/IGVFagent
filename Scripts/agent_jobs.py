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
6. When every stage is done or blocked, write the final report: what was
   reproduced (with numbers and file paths), what differs from the paper and
   why, and what is blocked.
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


def effective_status(job: dict) -> str:
    """'interrupted' when a running job's worker is gone (container restart, kill)."""
    st = job.get("status")
    if st in ("running", "queued"):
        hb = job.get("heartbeat") or 0
        if not _pid_alive(job.get("pid")) and (_now() - float(hb or 0)) > STALE_S:
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
                if m.get("is_error"):
                    res.error = res.answer or res.stop_reason
        return res

    def run_round(self, job: dict, prompt: str, history: "List[dict]") -> RoundResult:
        jid = job["id"]
        env = dict(os.environ)
        env.setdefault("IGVF_PROJECT_ROOT", str(ROOT))
        if not env.get("ANTHROPIC_API_KEY"):
            return RoundResult(error="ANTHROPIC_API_KEY is not set; the Claude Code orchestrator runs on the API key")
        try:
            proc = subprocess.Popen(self.argv(job, prompt), cwd=str(ROOT), env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, bufsize=1)
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


def _evidence_digest(plan: dict, limit_chars: int = 12000) -> str:
    out, used = [], 0
    for s in plan.get("stages") or []:
        out.append(f"## stage {s['id']} [{s['status']}] {s['title']}\nsuccess: {s.get('success')}\n"
                   f"harness: {(s.get('verified') or {}).get('reason', '-')}")
        for e in (s.get("evidence") or [])[:6]:
            q = _safe_path(e)
            if not q or not q.is_file():
                out.append(f"- {e}: (missing)")
                continue
            head = ""
            if q.suffix.lower() in (".md", ".json", ".tsv", ".csv", ".txt", ".log"):
                try:
                    head = q.read_text(errors="replace")[:1200]
                except OSError:
                    head = ""
            chunk = f"- {e} ({q.stat().st_size:,} bytes)\n" + ("```\n" + head + "\n```" if head else "")
            if used + len(chunk) > limit_chars:
                out.append(f"- {e} ({q.stat().st_size:,} bytes; excerpt omitted, budget)")
                continue
            used += len(chunk)
            out.append(chunk)
    return "\n".join(out)


def verify(job: dict, plan: dict, answer: str, llm_call: "Optional[Callable]" = None) -> dict:
    user = (f"# Task\n{job['query']}\n\n# Plan and evidence\n{_evidence_digest(plan)}\n\n"
            f"# Final report under review\n{answer[:8000]}")
    try:
        if llm_call is None:
            try:
                from igvfagent import _llm  # type: ignore
            except Exception:
                import _llm  # type: ignore
            backend = job.get("backend") if job.get("orchestrator") == "internal" else "anthropic"
            msg = _llm.chat([{"role": "system", "content": VERIFIER_PROMPT}, {"role": "user", "content": user}],
                            backend=backend or None, model=job.get("verifier_model") or job.get("model") or None,
                            max_tokens=1500, temperature=0.0)
            text = msg.content if isinstance(msg.content, str) else json.dumps(msg.content)
        else:
            text = llm_call(VERIFIER_PROMPT, user)
        m = re.search(r"\{.*\}", text or "", re.S)
        v = json.loads(m.group(0)) if m else {}
        verdict = "pass" if str(v.get("verdict", "")).lower() == "pass" else "fail"
        return {"verdict": verdict, "issues": [str(i) for i in (v.get("issues") or [])][:12]}
    except Exception as e:  # noqa: BLE001
        # an unavailable verifier must not pass work silently
        return {"verdict": "error", "issues": [f"verifier unavailable: {type(e).__name__}: {e}"]}


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
        update_job(job_id, status="queued", heartbeat=_now(), pid=os.getpid())
        sleep(30)
    if job.get("orchestrator") == "internal" and (job.get("backend") or os.environ.get("IGVF_LLM_BACKEND") or
                                                   "anthropic") == "anthropic":
        # Anthropic has no 128-function limit: a job should see the whole
        # registry, not the parity subset chat uses for OpenAI-family models
        os.environ["IGVF_LLM_MAX_TOOLS"] = os.environ.get("IGVF_JOB_MAX_TOOLS", "700")
    runner = runner or make_runner(job)
    job = update_job(job_id, status="running", pid=os.getpid(), heartbeat=_now(),
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
            event(job_id, "verify", f"verifier: {verdict['verdict']}", issues=verdict["issues"])
            if verdict["verdict"] != "pass" and vc < MAX_VERIFY_CYCLES:
                issues = [f"verifier: {i}" for i in verdict["issues"]] or ["verifier could not confirm the report"]
                job = update_job(job_id, verify_cycles=vc + 1, history=history, pending_issues=issues,
                                 last_answer=last_answer, last_verdict=verdict)
                continue
            blocked = any(s["status"] == "blocked" for s in plan.get("stages") or [])
            status = "done_with_blocked" if blocked else "done"
            if verdict["verdict"] != "pass":
                status = "done_with_blocked"
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
    return load_job(job_id)


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
