"""ReAct agent runtime for IGVFagent.

Wires the backend-neutral LLM router (``Scripts/_llm.py``) to the curated
tool registry (``Scripts/_tools.py``) into a single Plan → Action →
Results → Evaluation loop:

  user query  ──>  LLM(plan, tools)
                       │
                       ├── tool_calls?  ──>  for each: _tools.execute(...)
                       │                         │
                       │                         └─> append tool_result, loop
                       │
                       └── final text  ──>  return answer + transcript

The runtime emits a structured event stream through an optional callback
so the Streamlit UI (step 4) and any other front-end can render
incremental progress without re-implementing the loop. By default it
prints a compact human-readable trace to stdout.

Every run also persists a full transcript (JSON) plus a human-readable
markdown report under ``Docs/Agent/<timestamp>_<label>/`` so the audit
trail mirrors what the human-driven shell-based flow produces.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Callable, Optional

logger = logging.getLogger(__name__)


# --------------------------- Dual-mode imports ------------------------------

# Mirror the same dual-mode import pattern the other Scripts/* modules use
# so this works both as `from igvfagent._agent import run` (after
# `pip install -e .`) and as `import _agent` (when Scripts/ is on
# sys.path).
sys.path.insert(0, str(Path(__file__).resolve().parent))


def _load_helpers():
    try:
        from igvfagent import _llm as llm
        from igvfagent import _tools as tools
        from igvfagent import _router as router
        from igvfagent import _localstore as localstore
    except Exception:
        import _llm as llm  # type: ignore[no-redef]
        import _tools as tools  # type: ignore[no-redef]
        import _router as router  # type: ignore[no-redef]
        import _localstore as localstore  # type: ignore[no-redef]
    return llm, tools, router, localstore


def _grow_local_store(localstore, name: str, arguments: dict,
                      result: dict) -> None:
    """Feed one tool call into the growing local KG/DB (core default).

    Never let a store failure break the run; disable with IGVF_LOCALSTORE=0.
    """
    if os.environ.get("IGVF_LOCALSTORE", "1") == "0":
        return
    try:
        localstore.record_tool_call(
            name, arguments or {},
            stdout=result.get("stdout", "") or "",
            artifacts=result.get("artifacts") or {})
    except Exception as e:  # pragma: no cover
        logger.warning("localstore growth skipped for %s: %s", name, e)


# --------------------------- Project paths ----------------------------------

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
AGENT_DIR = ROOT / "Docs" / "Agent"
LOG_DIR = ROOT / "Docs" / "Logs"


# --------------------------- Event + result types --------------------------


@dataclasses.dataclass
class AgentEvent:
    """One step in the agent's run.

    ``kind`` is one of:
      * ``run_start``        — before the first LLM call
      * ``llm_call_start``   — about to call the LLM
      * ``llm_call_end``     — LLM responded (text + any tool_calls)
      * ``tool_call_start``  — about to invoke a tool
      * ``tool_call_end``    — tool returned (exit code + artefact paths)
      * ``final_answer``     — agent decided to stop and wrote text
      * ``error``             — exception during a tool call or LLM call
      * ``run_end``          — run finished (success or otherwise)
    """

    kind: str
    payload: dict


@dataclasses.dataclass
class AgentResult:
    final_answer: str
    iterations: int
    tool_calls_made: int
    stop_reason: str
    transcript: "list[dict]"
    artefacts: "list[str]"
    transcript_path: str
    report_path: str
    backend: str
    model: str


# --------------------------- System prompt ---------------------------------


DEFAULT_SYSTEM_PROMPT = """\
You are IGVFagent, a research assistant for IGVF and ENCODE data.

Workflow:
1. Pick the fewest tools that answer the user. Starting points:
   gene symbol -> kg_gene; IGVF/ENCODE accession or URL -> explain_dataset;
   discover datasets -> portal_kg_pull / splitseq_retrieve / encode_retrieve;
   prior literature check -> ref_validate; study-design -> ref_design.
2. Make 1-4 tool calls. Don't fabricate data; every claim must trace to a
   tool result.
3. Final answer: plain English, Markdown bullets, list the report / manifest
   file paths the tools produced, note caveats, suggest one next CLI call.
"""


# --------------------------- Helpers ---------------------------------------


def _print_callback(event: AgentEvent) -> None:
    """Default human-readable trace to stdout."""
    k = event.kind
    p = event.payload
    if k == "run_start":
        print(f"▶ run_start  backend={p.get('backend')} model={p.get('model')}  "
              f"tools={p.get('n_tools')}")
        c = p.get("consistency") or {}
        if c:
            print(f"  · consistency  seed={c.get('seed')} "
                  f"prompt={c.get('system_prompt_sha1')} "
                  f"toolset={c.get('tool_set_sha1')} ({c.get('n_tools_exposed')} tools)")
    elif k == "llm_call_start":
        print(f"  · llm  iter={p['iteration']}/{p['max_iterations']}  "
              f"messages={p['n_messages']}")
    elif k == "llm_call_end":
        text_preview = (p.get("content") or "").strip().splitlines()[:1]
        text_preview_str = text_preview[0][:120] if text_preview else ""
        tc = p.get("tool_calls") or []
        if tc:
            names = ", ".join(f"{t['name']}({_short_args(t.get('arguments'))})"
                                for t in tc)
            print(f"    -> tool_calls: {names}")
        if text_preview_str:
            print(f"    -> text: {text_preview_str}")
        if p.get("usage"):
            u = p["usage"]
            print(f"    [usage] in={u.get('input_tokens',0)} "
                  f"out={u.get('output_tokens',0)}")
    elif k == "route":
        print(f"  ⇒ route[{p.get('shape')}] -> {p['tool']}"
              f"({_short_args(p.get('arguments'))})  [deterministic]")
    elif k == "tool_call_start":
        tag = " [routed]" if p.get("routed") else ""
        print(f"  · tool {p['name']}({_short_args(p.get('arguments'))}){tag}")
    elif k == "tool_call_end":
        flag = "ok" if p.get("exit_code") == 0 else f"exit={p.get('exit_code')}"
        artefacts = p.get("artifacts") or {}
        n_paths = sum(len(v) for v in artefacts.values())
        print(f"    -> {flag}  artefacts={n_paths}")
    elif k == "final_answer":
        print()
        print("=" * 70)
        print(p["content"])
        print("=" * 70)
    elif k == "error":
        print(f"  ! error in {p.get('where','?')}: {p.get('error')}")
    elif k == "run_end":
        print()
        print(f"✓ run_end  iters={p['iterations']}  tools_called={p['tool_calls_made']}  "
              f"stop={p['stop_reason']}")
        if p.get("transcript_path"):
            print(f"  transcript: {p['transcript_path']}")
        if p.get("report_path"):
            print(f"  report:     {p['report_path']}")


def _short_args(args: Optional[dict]) -> str:
    if not args:
        return ""
    parts = []
    for k, v in list(args.items())[:4]:
        sv = str(v)
        if len(sv) > 30:
            sv = sv[:27] + "..."
        parts.append(f"{k}={sv}")
    if len(args) > 4:
        parts.append(f"+{len(args) - 4}")
    return ", ".join(parts)


def _emit(cb: Optional[Callable[[AgentEvent], None]], kind: str,
            payload: dict) -> None:
    if cb is None:
        return
    try:
        cb(AgentEvent(kind=kind, payload=payload))
    except Exception as e:  # never let a misbehaving callback kill the run
        logger.warning("callback failure on %s: %s", kind, e)


def _format_refusal(model: str, refused: bool = True) -> str:
    """Actionable message when the model returns a safety refusal or an
    empty non-answer on the first step (no tools run, nothing produced)."""
    what = ("declined the request (safety-classifier refusal)"
            if refused else "returned an empty response")
    return (
        f"**The model `{model or '(default)'}` {what} before any tools ran.**\n\n"
        "This is a known failure mode with **Claude Fable 5**: its additional "
        "dual-use safety layer false-positives on IGVFagent's genomics-agent "
        "system prompt and refuses even benign requests (the refusal is "
        "probabilistic, so it is not fixable by rewording the prompt).\n\n"
        "**Use a model without that extra layer** — recommended for all "
        "IGVFagent work:\n"
        "- `claude-sonnet-4-5` (default fallback) or `claude-opus-4-8`\n"
        "- set `IGVF_LLM_MODEL=claude-sonnet-4-5`, pick Sonnet/Opus in the UI "
        "sidebar, or set `IGVF_LLM_FALLBACK_MODEL` to auto-retry on refusal.\n"
    )


def _format_runtime_error(err_msg: str, backend: str, model: str) -> str:
    """Turn a raw exception string into an actionable Markdown report.

    Pattern-matches the most common Ollama / Anthropic / OpenAI failures
    so end users in the Streamlit UI see a useful next step rather than
    an empty 'final answer' panel.
    """
    em = (err_msg or "").lower()
    msg = (
        f"**The agent could not complete the run.**\n\n"
        f"- Backend: `{backend or '(unknown)'}`\n"
        f"- Model: `{model or '(default)'}`\n"
        f"- Underlying error: `{err_msg}`\n\n"
    )
    if "more system memory" in em or "requires more" in em \
            and "memory" in em:
        msg += (
            "### Likely cause: not enough RAM for the chosen model\n\n"
            "Pick a smaller Ollama model in the sidebar (or via "
            "`IGVF_LLM_MODEL` in the Docker `.env`). Roughly memory-safe "
            "options:\n\n"
            "| Model | Approx RAM | Notes |\n"
            "|---|---|---|\n"
            "| `qwen3:0.6b`  | ~1 GB | smallest, weak tool calls |\n"
            "| `qwen2.5:1.5b` | ~2 GB | reasonable for simple tools |\n"
            "| `qwen3:4b` | ~3 GB | good balance |\n"
            "| `qwen3:8b` | ~6 GB | the project default |\n"
            "| `qwen3.6:35b-*` | ~65 GB | needs a real workstation |\n\n"
            "If you're on Docker, also bump container memory: "
            "Docker Desktop → Settings → Resources → Memory ≥ 8 GB.\n"
        )
    elif "401" in em or "unauthor" in em or "invalid api key" in em \
            or "authentication" in em:
        msg += (
            "### Likely cause: missing / wrong API key\n\n"
            "- Anthropic: set `ANTHROPIC_API_KEY` in the environment.\n"
            "- OpenAI:    set `OPENAI_API_KEY`.\n"
            "- For Docker Compose: put the key in a host-local `.env` and "
            "`docker compose up` will forward it.\n"
        )
    elif "404" in em or "not found" in em or "no such model" in em:
        msg += (
            "### Likely cause: model name not recognised by the backend\n\n"
            f"- For Ollama, ensure the model is pulled: "
            f"`ollama pull {model or 'qwen3:8b'}`.\n"
            "- For Anthropic, use a Claude model name (e.g. "
            "`claude-sonnet-4-5`). Putting a Qwen / Gemma name in the "
            "Anthropic backend will always fail.\n"
            "- Run `igvfagent models` to see what your local Ollama "
            "daemon actually has installed.\n"
        )
    elif "429" in em or "rate limit" in em or "rate_limit" in em:
        msg += (
            "### Likely cause: provider rate limit\n\n"
            "Wait ~30 seconds and re-run. For higher throughput either "
            "use a paid tier or switch to a local Ollama model.\n"
        )
    elif "connection" in em or "refused" in em or "timed out" in em \
            or "name resolution" in em:
        msg += (
            "### Likely cause: backend not reachable\n\n"
            "- Ollama: `ollama serve` running? "
            "Check `OLLAMA_HOST_BASE` (default `http://localhost:11434/v1`; "
            "in Docker Compose: `http://ollama:11434/v1`).\n"
            "- Cloud providers: verify network egress is allowed.\n"
        )
    else:
        msg += (
            "### Diagnosis tips\n\n"
            "- Run `igvfagent backends` to confirm your provider is "
            "configured.\n"
            "- Run `igvfagent models` to list Ollama models on the "
            "local daemon.\n"
            "- Check `Docs/Logs/agent_*.log` for the full traceback.\n"
        )
    return msg


def _format_tool_result_for_llm(result: dict, max_chars: int = 3500) -> str:
    """Compact, model-friendly summary of a tool execution."""
    parts: "list[str]" = []
    parts.append(f"exit_code={result.get('exit_code')}")
    artefacts = result.get("artifacts") or {}
    if artefacts:
        for k, vs in artefacts.items():
            for v in vs[:3]:
                parts.append(f"{k}: {v}")
    stdout = (result.get("stdout") or "").strip()
    if stdout:
        # Drop log-prefix lines that are noise to the model.
        useful = "\n".join(
            ln for ln in stdout.splitlines()
            if not ln.startswith("20") or "INFO" not in ln
        ).strip()
        if useful:
            if len(useful) > max_chars:
                useful = useful[: max_chars - 60] + "\n...[truncated]"
            parts.append("stdout:")
            parts.append(useful)
    stderr = (result.get("stderr") or "").strip()
    if stderr and result.get("exit_code") != 0:
        parts.append(f"stderr: {stderr[:600]}")
    return "\n".join(parts)


def _compose_templated_answer(transcript: "list[dict]", artefacts: "list[str]",
                               llm_text: str) -> str:
    """Deterministic answer skeleton from the (backend-independent) tool trace.

    The tools run, their arguments, exit codes, and the artefacts they produce
    are identical across backends for a consistent plan — so templating them
    makes the *substantive* answer identical regardless of which LLM drove the
    loop. The model's free-form prose is kept, clearly labelled, as a summary.
    """
    calls: "list[tuple[str, dict]]" = []
    results: "dict[str, dict]" = {}
    for entry in transcript:
        if entry.get("role") == "assistant":
            for tc in entry.get("tool_calls", []) or []:
                calls.append((tc["name"], tc.get("arguments") or {}))
        elif entry.get("role") == "tool":
            results[entry.get("name", "")] = entry

    lines: "list[str]" = []
    lines.append("## What was run")
    if calls:
        for name, args in calls:
            arg_str = ", ".join(f"{k}={v}" for k, v in list(args.items())[:4])
            lines.append(f"- `{name}({arg_str})`")
    else:
        lines.append("- _(no tools were called)_")

    # Artefacts: deterministic, de-duplicated, sorted.
    uniq = sorted(set(artefacts))
    if uniq:
        lines += ["", "## Artefacts produced"]
        for a in uniq[:40]:
            lines.append(f"- `{a}`")
        if len(uniq) > 40:
            lines.append(f"- …and {len(uniq) - 40} more")

    lines += ["", "## Summary", "",
              (llm_text or "_(no narrative generated)_").strip()]
    return "\n".join(lines)


# --------------------------- Run loop --------------------------------------


def _label_from_query(query: str) -> str:
    digest = hashlib.sha1(query.encode("utf-8")).hexdigest()[:8]
    head = "".join(c if c.isalnum() else "_" for c in query)[:40].strip("_")
    return f"{head or 'query'}_{digest}"


def _persist_transcript(query: str, transcript: "list[dict]",
                          final_answer: str, meta: dict) -> "tuple[str, str]":
    AGENT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    label = _label_from_query(query)
    out_dir = AGENT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    transcript_path = out_dir / "transcript.json"
    transcript_path.write_text(json.dumps({
        "query":         query,
        "transcript":    transcript,
        "final_answer":  final_answer,
        **meta,
    }, indent=2, default=str))
    report_path = out_dir / "report.md"
    lines = [
        f"# IGVFagent run — {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"**Query:** {query}",
        "",
        f"- Backend: `{meta.get('backend','')}`  ·  Model: `{meta.get('model','')}`",
        f"- Iterations: {meta.get('iterations',0)}  ·  Tool calls: {meta.get('tool_calls_made',0)}",
        f"- Stop reason: `{meta.get('stop_reason','')}`",
        f"- Transcript JSON: `{transcript_path.relative_to(ROOT)}`",
        "",
        "## Final answer",
        "",
        final_answer or "_no final text was generated_",
        "",
    ]
    if meta.get("artefacts"):
        lines += ["## Artefacts referenced", ""]
        for a in meta["artefacts"]:
            lines.append(f"- `{a}`")
        lines.append("")
    report_path.write_text("\n".join(lines))
    return str(transcript_path), str(report_path)


def run(
    query: str,
    *,
    backend: Optional[str] = None,
    model: Optional[str] = None,
    max_iterations: int = 8,
    max_tokens: int = 4096,
    temperature: float = 0.0,
    seed: Optional[int] = None,
    enable_router: bool = True,
    templated_answer: bool = True,
    system_prompt: Optional[str] = None,
    extra_context: Optional[str] = None,
    tools_subset: Optional["list[str]"] = None,
    callback: Optional[Callable[[AgentEvent], None]] = _print_callback,
    persist: bool = True,
) -> AgentResult:
    """Run a Plan → Action → Results → Evaluation loop on a user query.

    Returns an :class:`AgentResult` with the final text, full transcript,
    and any artefact paths the wrapped skills announced.
    """
    llm, tools_mod, router, localstore = _load_helpers()

    all_tools = tools_mod.list_tools()
    if tools_subset:
        names = set(tools_subset)
        all_tools = [t for t in all_tools if t.name in names]
    tool_dicts = [t.to_dict() for t in all_tools]

    sys_prompt = system_prompt or DEFAULT_SYSTEM_PROMPT

    # Resolve the decoding seed once (default IGVF_LLM_SEED, else 0) so the
    # same value is used for every LLM call in this run and recorded below.
    if seed is None:
        _env_seed = os.environ.get("IGVF_LLM_SEED", "0")
        try:
            seed = int(_env_seed) if _env_seed != "" else 0
        except ValueError:
            seed = 0

    # Consistency fingerprint: the three inputs that must match for two runs
    # (on any backends) to be comparable — system prompt, decoding seed, and
    # the *canonical* tool set the LLM actually sees (identical across backends).
    canonical = llm.canonical_tools(tool_dicts) or tool_dicts
    canon_names = [t.get("name", "") for t in canonical]
    canon_name_set = set(canon_names)

    # Deterministic router: fixed first tool call(s) for unambiguous query
    # shapes, identical across backends. Only keep routed calls whose tool is
    # actually exposed (survived the canonical cap).
    routed_plan: "list[dict]" = []
    if enable_router and os.environ.get("IGVF_ROUTER", "1") != "0":
        routed_plan = [r for r in router.route(query)
                       if r.get("tool") in canon_name_set]

    consistency = {
        "system_prompt_sha1": hashlib.sha1(sys_prompt.encode()).hexdigest()[:12],
        "seed": seed,
        "n_tools_exposed": len(canon_names),
        "tool_set_sha1": hashlib.sha1(
            ",".join(sorted(canon_names)).encode()).hexdigest()[:12],
        "routed_shapes": [r["shape"] for r in routed_plan],
    }

    user_text = query if not extra_context else f"{query}\n\n{extra_context}"
    messages: "list[dict]" = [{"role": "user", "content": user_text}]

    transcript: "list[dict]" = [{"role": "user", "content": user_text}]
    artefacts: "list[str]" = []
    final_answer = ""
    stop_reason = "max_iterations"
    iters = 0
    tool_calls_made = 0
    chosen_backend = backend or os.environ.get("IGVF_LLM_BACKEND") or ""
    chosen_model = model or ""
    refusal_fallback_used = False
    # No-progress detection: bail out (instead of burning the whole iteration
    # budget) when the model repeats the identical tool-call set, or every tool
    # call fails, for this many consecutive iterations.
    stuck_limit = max(2, int(os.environ.get("IGVF_AGENT_STUCK_LIMIT", "3")))
    last_sig = None
    repeat_count = 0
    consecutive_allfail = 0

    _emit(callback, "run_start", {
        "backend": chosen_backend or "(auto)", "model": chosen_model or "(auto)",
        "n_tools": len(canon_names), "max_iterations": max_iterations,
        "consistency": consistency,
    })

    # --- Deterministic pre-plan: execute routed tool calls first -----------
    # These run identically on every backend, seeding the conversation before
    # the LLM takes over. We synthesize an assistant tool_use turn + tool
    # result so the message history is valid for both Anthropic and OpenAI.
    if routed_plan:
        synth_calls = []
        for i, r in enumerate(routed_plan):
            call_id = f"route_{i}"
            _emit(callback, "route", {"shape": r["shape"], "tool": r["tool"],
                                       "arguments": r["arguments"]})
            _emit(callback, "tool_call_start", {
                "id": call_id, "name": r["tool"], "arguments": r["arguments"],
                "routed": True,
            })
            try:
                result = tools_mod.execute(r["tool"], r["arguments"])
                tool_calls_made += 1
            except Exception as e:  # noqa
                result = {"name": r["tool"], "exit_code": 1, "stdout": "",
                          "stderr": str(e), "artifacts": {}}
            for paths in (result.get("artifacts") or {}).values():
                artefacts.extend(paths)
            _emit(callback, "tool_call_end", {
                "id": call_id, "name": r["tool"],
                "exit_code": result.get("exit_code"),
                "artifacts": result.get("artifacts") or {},
                "routed": True,
            })
            _grow_local_store(localstore, r["tool"], r["arguments"], result)
            synth_calls.append((call_id, r, result))

        # assistant turn announcing the routed calls
        assistant_tc = [{"id": cid, "name": r["tool"], "arguments": r["arguments"]}
                        for cid, r, _ in synth_calls]
        routed_note = ("(router) deterministic tool selection for query shape(s): "
                       + ", ".join(r["shape"] for _, r, _ in synth_calls))
        transcript.append({"role": "assistant", "content": routed_note,
                            "tool_calls": assistant_tc, "routed": True})
        messages.append({"role": "assistant", "content": routed_note,
                         "tool_calls": assistant_tc})
        for cid, r, result in synth_calls:
            summary = _format_tool_result_for_llm(result)
            transcript.append({"role": "tool", "tool_call_id": cid,
                               "name": r["tool"], "content": summary,
                               "exit_code": result.get("exit_code"),
                               "artifacts": result.get("artifacts") or {},
                               "routed": True})
            messages.append({"role": "tool", "tool_call_id": cid,
                             "content": summary})

    for iters in range(1, max_iterations + 1):
        _emit(callback, "llm_call_start", {
            "iteration": iters, "max_iterations": max_iterations,
            "n_messages": len(messages),
        })
        try:
            msg = llm.chat(
                messages=[{"role": "system", "content": sys_prompt}, *messages],
                backend=backend, model=model, tools=tool_dicts,
                max_tokens=max_tokens, temperature=temperature, seed=seed,
            )
        except Exception as e:
            err_msg = str(e)
            _emit(callback, "error",
                  {"where": "llm.chat", "error": err_msg})
            stop_reason = "error"
            final_answer = _format_runtime_error(
                err_msg, chosen_backend or backend or "",
                chosen_model or model or "",
            )
            break
        chosen_backend = msg.backend
        chosen_model = msg.model

        # Record assistant turn (with any tool calls) into the transcript.
        assistant_entry: dict = {
            "role": "assistant",
            "content": msg.content,
            "stop_reason": msg.stop_reason,
            "backend": msg.backend, "model": msg.model,
        }
        if msg.tool_calls:
            assistant_entry["tool_calls"] = [
                {"id": tc.id, "name": tc.name, "arguments": tc.arguments}
                for tc in msg.tool_calls
            ]
        transcript.append(assistant_entry)
        _emit(callback, "llm_call_end", {
            "iteration": iters,
            "content": msg.content,
            "tool_calls": assistant_entry.get("tool_calls", []),
            "stop_reason": msg.stop_reason,
            "usage": msg.usage,
        })

        # A non-answer from the model: either an explicit safety refusal
        # (Fable 5's dual-use classifier false-positives on the genomics
        # system prompt → stop_reason "refusal", empty content) or an empty
        # first response with no tools called. Previously this fell straight
        # into the "no tool_calls → complete" branch below and was reported
        # as a successful run with an empty answer, hiding the real cause.
        # Fall back once to a reliable model, then surface a clear message.
        refused = msg.stop_reason == "refusal"
        empty_start = (iters == 1 and tool_calls_made == 0
                       and not msg.tool_calls
                       and not (msg.content or "").strip())
        if refused or empty_start:
            fallback = os.environ.get("IGVF_LLM_FALLBACK_MODEL",
                                      "claude-sonnet-4-5")
            can_fallback = (
                not refusal_fallback_used
                and (chosen_backend or backend) == "anthropic"
                and (model or chosen_model or "") != fallback
            )
            _emit(callback, "refusal", {
                "iteration": iters,
                "model": chosen_model or model,
                "stop_reason": msg.stop_reason,
                "falling_back_to": fallback if can_fallback else None,
            })
            if can_fallback:
                refusal_fallback_used = True
                model = fallback          # retry this step on the reliable model
                continue
            final_answer = _format_refusal(chosen_model or model or "", refused)
            stop_reason = "refusal"
            break

        if not msg.tool_calls:
            final_answer = msg.content or ""
            stop_reason = "complete"
            break

        # Mirror assistant message into the conversation for the next LLM call.
        messages.append({
            "role": "assistant", "content": msg.content,
            "tool_calls": [{"id": tc.id, "name": tc.name,
                              "arguments": tc.arguments}
                            for tc in msg.tool_calls],
        })

        # Signature of this iteration's requested work, for no-progress
        # detection (identical repeated tool-call sets = spinning).
        iter_sig = tuple(sorted(
            (tc.name, repr(sorted((tc.arguments or {}).items())))
            for tc in msg.tool_calls))
        iter_n = 0
        iter_fail = 0

        # Execute every tool the model asked for.
        for tc in msg.tool_calls:
            _emit(callback, "tool_call_start", {
                "id": tc.id, "name": tc.name, "arguments": tc.arguments,
            })
            try:
                result = tools_mod.execute(tc.name, tc.arguments or {})
                tool_calls_made += 1
            except KeyError:
                result = {
                    "name": tc.name, "exit_code": 127,
                    "stdout": "", "stderr": f"Unknown tool: {tc.name}",
                    "artifacts": {},
                }
            except Exception as e:
                _emit(callback, "error",
                      {"where": f"tool:{tc.name}", "error": str(e)})
                result = {
                    "name": tc.name, "exit_code": 1,
                    "stdout": "", "stderr": str(e),
                    "artifacts": {},
                }
            iter_n += 1
            if int(result.get("exit_code") or 0) != 0:
                iter_fail += 1
            for paths in (result.get("artifacts") or {}).values():
                artefacts.extend(paths)
            _emit(callback, "tool_call_end", {
                "id": tc.id, "name": tc.name,
                "exit_code": result.get("exit_code"),
                "artifacts": result.get("artifacts") or {},
                "stdout_len": len(result.get("stdout") or ""),
            })
            _grow_local_store(localstore, tc.name, tc.arguments or {}, result)

            tool_summary = _format_tool_result_for_llm(result)
            transcript.append({
                "role": "tool", "tool_call_id": tc.id,
                "name": tc.name, "content": tool_summary,
                "exit_code": result.get("exit_code"),
                "artifacts": result.get("artifacts") or {},
            })
            messages.append({
                "role": "tool", "tool_call_id": tc.id,
                "content": tool_summary,
            })

        # No-progress detection: stop early (rather than burning the whole
        # iteration budget) if the model keeps issuing the identical tool-call
        # set, or every tool call fails, for `stuck_limit` iterations in a row.
        if iter_sig == last_sig:
            repeat_count += 1
        else:
            repeat_count = 1
            last_sig = iter_sig
        consecutive_allfail = (consecutive_allfail + 1
                               if iter_n and iter_fail == iter_n else 0)
        if repeat_count >= stuck_limit or consecutive_allfail >= stuck_limit:
            why = ("the model repeated the same tool call(s) without making "
                   "progress" if repeat_count >= stuck_limit
                   else "every tool call failed for several iterations")
            _emit(callback, "stuck", {
                "iteration": iters, "reason": why,
                "repeat_count": repeat_count,
                "consecutive_allfail": consecutive_allfail,
            })
            stop_reason = "stuck"
            break

    # Final-iteration wrap-up. The loop exhausted its tool-call budget while
    # the model was still calling tools, so nothing was ever synthesized. Make
    # one more LLM call with tools DISABLED and an explicit instruction to
    # answer from the evidence already gathered — this rescues the findings
    # (counts, IDs, file paths) instead of discarding them behind a boilerplate
    # "ran out of iterations" message.
    if not final_answer and stop_reason in ("max_iterations", "stuck"):
        wrap_prompt = (
            "You have reached the tool-call budget for this task and cannot "
            "call any more tools. Do NOT request any tools. Using ONLY the "
            "information already gathered above, write the best final answer "
            "you can to the user's original request: summarize what you "
            "found, state concrete results (counts, IDs, file paths, key "
            "facts), say clearly what remains incomplete, and give the user "
            "concrete next steps to finish the task."
        )
        _emit(callback, "wrap_up_start", {"iteration": iters})
        try:
            wrap_msg = llm.chat(
                messages=[{"role": "system", "content": sys_prompt}, *messages,
                          {"role": "user", "content": wrap_prompt}],
                backend=backend, model=model, tools=None,
                max_tokens=max_tokens, temperature=temperature, seed=seed,
            )
        except Exception as e:  # noqa: BLE001 — best-effort rescue
            _emit(callback, "error",
                  {"where": "llm.chat.wrap_up", "error": str(e)})
            wrap_msg = None
        if wrap_msg is not None and (wrap_msg.content or "").strip():
            final_answer = wrap_msg.content
            stop_reason = "max_iterations_wrapped"
            chosen_backend = wrap_msg.backend
            chosen_model = wrap_msg.model
            transcript.append({"role": "user", "content": wrap_prompt,
                               "wrap_up": True})
            transcript.append({
                "role": "assistant", "content": wrap_msg.content,
                "stop_reason": wrap_msg.stop_reason,
                "backend": wrap_msg.backend, "model": wrap_msg.model,
                "wrap_up": True,
            })
            _emit(callback, "wrap_up_end", {
                "content": wrap_msg.content,
                "stop_reason": wrap_msg.stop_reason,
                "usage": wrap_msg.usage,
            })

    # Truncated final answer if the wrap-up call also produced nothing.
    if not final_answer and stop_reason in ("max_iterations",
                                            "max_iterations_wrapped", "stuck"):
        final_answer = (
            "_The agent stopped without producing a final answer "
            "(ran out of iterations, or was not making progress). See the "
            "transcript and any artefacts written so far._"
        )
    elif not final_answer and stop_reason == "error":
        final_answer = (
            "_The agent failed before producing a final answer. "
            "See the transcript / log for the underlying exception._"
        )

    # Deterministic templated synthesis: wrap the model's prose in a
    # backend-independent skeleton (tools run + artefacts) so the substantive
    # answer is identical across LLMs. Only when tools actually ran and the run
    # did not error. Disable with IGVF_TEMPLATED_ANSWER=0.
    if (templated_answer and tool_calls_made > 0 and stop_reason == "complete"
            and os.environ.get("IGVF_TEMPLATED_ANSWER", "1") != "0"):
        final_answer = _compose_templated_answer(transcript, artefacts,
                                                  final_answer)

    transcript_path = ""
    report_path = ""
    if persist:
        transcript_path, report_path = _persist_transcript(
            query, transcript, final_answer,
            meta={
                "backend": chosen_backend, "model": chosen_model,
                "iterations": iters, "tool_calls_made": tool_calls_made,
                "stop_reason": stop_reason, "artefacts": artefacts,
                "consistency": consistency,
            },
        )

    _emit(callback, "final_answer", {"content": final_answer})
    _emit(callback, "run_end", {
        "iterations": iters, "tool_calls_made": tool_calls_made,
        "stop_reason": stop_reason, "artefacts": artefacts,
        "transcript_path": transcript_path, "report_path": report_path,
    })

    return AgentResult(
        final_answer=final_answer,
        iterations=iters,
        tool_calls_made=tool_calls_made,
        stop_reason=stop_reason,
        transcript=transcript,
        artefacts=artefacts,
        transcript_path=transcript_path,
        report_path=report_path,
        backend=chosen_backend,
        model=chosen_model,
    )


__all__ = ["run", "AgentResult", "AgentEvent",
           "DEFAULT_SYSTEM_PROMPT"]
