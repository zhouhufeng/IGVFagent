"""Web UI for durable agent jobs (agent_jobs.py).

Three pieces for the chat page:

- ``route(query, mode)``: should this message start a background job? In
  ``auto`` mode, tasks that do not fit one reply (reproduce, from raw reads,
  full pipeline, download and analyse, end to end) become jobs; ``always``
  and ``never`` override.
- ``is_continue(query)``: "continue", "keep going", "resume the
  reproduction" resume the user's latest unfinished job from its plan,
  instead of starting cold.
- ``render_panel(st, viewer, is_admin, render_file)``: the Jobs panel. For
  each of the user's recent jobs it shows the state, the harness-verified
  plan checklist, the latest events, the verifier's verdict and, when the
  job is done, its report and files. It has Stop and Resume buttons and
  refreshes itself while a job runs.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, List, Optional

try:
    from igvfagent import agent_jobs as aj  # type: ignore
except Exception:  # checkout
    import agent_jobs as aj  # type: ignore

LONG_TASK = re.compile(
    r"\b(reproduc\w*|replicat\w*|re-?run the (?:paper|analysis|pipeline)|from (?:the )?raw (?:reads|data|fastq)|"
    r"end[- ]to[- ]end|full (?:pipeline|analysis|workflow)|download .{0,60}\b(?:and|then) (?:analy[sz]e|process|run)|"
    r"benchmark (?:it|this|the paper)|agentify|step[- ]by[- ]step (?:analysis|reproduction))\b", re.I)
CONTINUE = re.compile(r"^\s*(?:ok[,.! ]*|yes[,.! ]*|please[,.! ]*)*(?:continue|keep going|go on|carry on|resume|"
                      r"proceed|finish (?:it|the (?:job|reproduction|analysis)))\b", re.I)

STATUS_ICON = {"running": "🔄", "queued": "⏳", "done": "✅", "done_with_blocked": "🟡", "failed": "❌",
               "stopped": "⏹", "budget_exhausted": "⌛", "interrupted": "⚠️"}
STAGE_ICON = {"done": "✅", "blocked": "⛔", "failed": "❌", "running": "🔄", "pending": "▫️"}


def route(query: str, mode: str = "auto") -> bool:
    if mode == "always":
        return True
    if mode == "never":
        return False
    return bool(LONG_TASK.search(query or ""))


def is_continue(query: str) -> bool:
    return bool(CONTINUE.match(query or "")) and len(query) < 400


def latest_unfinished(viewer: Optional[str], is_admin: bool) -> Optional[dict]:
    for j in aj.list_jobs(viewer, is_admin, limit=10):
        if j["effective_status"] not in ("done",):
            return j
    return None


def job_query(query: str, prior: "List[dict]") -> str:
    """The job gets the task plus the recent conversation it refers to."""
    ctx = []
    for m in prior[-6:]:
        c = (m.get("content") or "").strip()
        if c:
            ctx.append(f"[{m.get('role')}] {c[:1500]}")
    return query if not ctx else f"{query}\n\n--- Conversation so far (for context) ---\n" + "\n\n".join(ctx)


def start_from_chat(query: str, prior: "List[dict]", *, owner: str, owner_admin: bool, orchestrator: str,
                    backend: Optional[str], model: Optional[str]) -> dict:
    return aj.start_job(job_query(query, prior), owner=owner, owner_admin=owner_admin,
                        orchestrator=orchestrator, backend=backend, model=model,
                        title=query.strip().splitlines()[0][:100])


def _stage_lines(plan: dict) -> "List[str]":
    out = []
    for s in plan.get("stages") or []:
        line = f"{STAGE_ICON.get(s['status'], '▫️')} **{s['id']}** {s['title']}"
        v = s.get("verified") or {}
        if s["status"] == "failed" and v.get("reason"):
            line += f" — _{v['reason'][:140]}_"
        if s["status"] == "blocked" and s.get("notes"):
            line += f" — _blocked: {s['notes'][-1][21:][:140]}_"
        if s.get("attempts"):
            line += f" (attempt {s['attempts']})"
        out.append(line)
    return out


def render_job(st, j: dict, render_file: Optional[Callable[[str], None]] = None, key: str = "") -> None:
    st_ = j.get("effective_status") or aj.effective_status(j)
    plan = aj.load_plan(j["id"])
    stages = plan.get("stages") or []
    done = sum(1 for s in stages if s["status"] in ("done", "blocked"))
    head = (f"{STATUS_ICON.get(st_, '•')} **{j.get('title', '')[:90]}** · `{j['id']}` · {st_.replace('_', ' ')} · "
            f"round {j.get('rounds', 0)} · {j.get('orchestrator')}")
    if float(j.get("cost_usd") or 0):
        head += f" · ${float(j['cost_usd']):.2f}"
    with st.expander(head, expanded=st_ in ("running", "queued")):
        if stages:
            st.progress(done / len(stages), text=f"{done}/{len(stages)} stages done or blocked")
            st.markdown("\n\n".join(_stage_lines(plan)))
        else:
            st.caption("Planning…" if st_ in ("running", "queued") else "No plan was recorded.")
        v = j.get("last_verdict") or {}
        if v:
            st.caption(f"Independent verifier: **{v.get('verdict')}**"
                       + (" — " + "; ".join(v.get("issues") or [])[:400] if v.get("issues") else ""))
        ev = aj.read_events(j["id"], 8)
        if ev and st_ in ("running", "queued", "interrupted"):
            st.caption("Latest: " + " · ".join(f"{e['kind']} {e.get('text', '')[:60]}" for e in ev[-5:]))
        ans = aj.job_dir(j["id"]) / "answer.md"
        if ans.exists() and st_ not in ("running", "queued"):
            st.markdown(ans.read_text()[:20000])
        elif j.get("last_answer") and st_ not in ("running", "queued"):
            st.markdown(j["last_answer"][-4000:])
        arts = [a for a in (j.get("artefacts") or []) if Path(a).suffix.lower() in (".png", ".svg", ".md", ".tsv", ".csv")]
        if render_file and arts and st_ not in ("running", "queued"):
            with st.expander(f"📁 Files ({len(arts)})", expanded=False):
                for a in arts[-12:]:
                    try:
                        render_file(a)
                    except Exception:  # noqa: BLE001
                        st.caption(a)
        c1, c2, _ = st.columns([1, 1, 4])
        if st_ in ("running", "queued"):
            if c1.button("⏹ Stop", key=f"stop_{j['id']}{key}"):
                aj.stop_job(j["id"])
                st.rerun()
        elif st_ in aj.RESUMABLE:
            if c1.button("▶ Resume", key=f"resume_{j['id']}{key}"):
                aj.resume_job(j["id"])
                st.rerun()


def render_panel(st, viewer: Optional[str], is_admin: bool,
                 render_file: Optional[Callable[[str], None]] = None, limit: int = 5) -> None:
    try:
        jobs = aj.list_jobs(viewer, is_admin, limit=limit)
    except Exception:  # noqa: BLE001
        return
    if not jobs:
        return
    st.markdown("#### 🕒 Your long-running jobs")
    st.caption("Jobs keep running when you close this page. Say “continue” to resume the latest one.")
    for j in jobs:
        render_job(st, j, render_file)
