"""Streamlit UI for IGVFagent.

Browser-based chat front-end that drives the natural-language ReAct
agent (``Scripts/_agent.py``), streams progress events as the agent
plans and calls tools, and renders artefact outputs (UMAP / dot-plot /
volcano PNG/SVG figures, markdown reports, manifest CSVs) inline.

Run via either:

    igvfagent ui                       # preferred (handled by cli.py)
    streamlit run Scripts/streamlit_app.py

The UI is intentionally thin: every action it takes is also reachable
from the terminal. It does not duplicate any logic; the LLM router,
tool registry, and agent loop are all owned by ``_llm.py``,
``_tools.py``, and ``_agent.py`` respectively.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import streamlit as st

# Dual-mode import (installed package OR running from a checkout).
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from igvfagent import _agent, _llm, _tools, __version__
except Exception:
    import _agent  # type: ignore
    import _llm    # type: ignore
    import _tools  # type: ignore
    try:
        from __init__ import __version__  # type: ignore
    except Exception:
        __version__ = "0.1.0"


# --------------------------- Page config -----------------------------------

st.set_page_config(
    page_title="IGVFagent",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# --------------------------- Sidebar config --------------------------------

def _resolve_effective_config(backend_choice: str, model_input: str
                                  ) -> "dict":
    """Compute exactly what _agent.run() will hit — backend + model +
    credential presence + base URL — so the sidebar can display it."""
    import os
    backend_arg = None if backend_choice == "(auto)" else backend_choice
    eff_backend = _llm._resolve_backend(backend_arg, model_input or None)
    default_model = _llm._DEFAULT_MODELS.get(eff_backend, "qwen3:8b")
    eff_model = (model_input or os.environ.get("IGVF_LLM_MODEL")
                 or default_model)
    desc = _llm.describe_backend(eff_backend)
    key_env = desc.get("api_key_env") if isinstance(desc, dict) else None
    key_set = bool(os.environ.get(key_env)) if key_env else None
    base_url = desc.get("base_url") if isinstance(desc, dict) else None
    return {
        "backend":      eff_backend,
        "model":        eff_model,
        "key_env":      key_env,
        "key_set":      key_set,
        "base_url":     base_url,
        "model_source": ("user input" if model_input else
                          ("IGVF_LLM_MODEL env" if os.environ.get("IGVF_LLM_MODEL")
                           else "backend default")),
    }


def _sidebar() -> dict:
    with st.sidebar:
        st.markdown(f"## 🧬 IGVFagent\n_v{__version__}_")
        st.caption(
            "Natural-language interface to the IGVF / ENCODE single-cell, "
            "variant, regulatory-element, and literature stack."
        )
        st.divider()

        st.subheader("Backend")
        backends = ["(auto)"] + _llm.list_backends()
        # Reset the model field whenever the backend changes so a Qwen
        # model name doesn't accidentally carry over into the Anthropic
        # backend (or vice versa).
        prev_backend = st.session_state.get("_prev_backend")
        backend = st.selectbox(
            "LLM provider", backends, index=0,
            help="`(auto)` infers from the model name; defaults to Ollama "
                 "(Qwen 3 8B) when nothing else is configured.",
        )
        if prev_backend is not None and prev_backend != backend:
            st.session_state.pop("_model_input", None)
        st.session_state["_prev_backend"] = backend

        # Suggest a sensible default model placeholder per backend.
        placeholders = {
            "(auto)":   "default for backend (e.g. qwen3:8b)",
            "anthropic": "claude-sonnet-4-5  ← needs ANTHROPIC_API_KEY",
            "openai":    "gpt-4o-mini  ← needs OPENAI_API_KEY",
            "ollama":    "qwen3:8b  (run `igvfagent models` to list yours)",
            "groq":      "llama-3.1-70b-versatile",
            "together":  "meta-llama/Meta-Llama-3.1-70B-Instruct-Turbo",
            "vllm":      "qwen3:8b  (or whatever your vLLM serves)",
            "tgi":       "(your HuggingFace TGI model id)",
        }
        model = st.text_input(
            "Model",
            value="",
            key="_model_input",
            placeholder=placeholders.get(backend, "default for backend"),
        )

        # Resolve and display the effective configuration so the user
        # always sees exactly what the agent run will hit.
        eff = _resolve_effective_config(backend, model)
        st.markdown("**Resolved configuration**")
        st.markdown(
            f"- Backend: `{eff['backend']}`\n"
            f"- Model: `{eff['model']}`  _({eff['model_source']})_"
        )
        if eff["key_env"]:
            icon = ("✅" if eff["key_set"] else
                    ("➖" if eff["backend"] == "ollama" else "❌"))
            st.markdown(
                f"- Credential: `{eff['key_env']}` {icon}"
                + ("" if eff["key_set"] else
                    "  _(unset — set in your shell or .env)_")
            )
        if eff["base_url"]:
            st.markdown(f"- Endpoint: `{eff['base_url']}`")

        # Health check
        if st.button("🔍 Test backend", use_container_width=True):
            with st.spinner("Pinging backend…"):
                try:
                    msg = _llm.chat(
                        messages=[{"role": "user", "content": "ping"}],
                        backend=None if backend == "(auto)" else backend,
                        model=model or None,
                        max_tokens=16, temperature=0.0,
                    )
                    st.session_state["_health_status"] = {
                        "ok": True,
                        "msg": (f"backend `{msg.backend}` · "
                                f"model `{msg.model}` · responded"),
                        "ts": time.strftime("%H:%M:%S"),
                    }
                except Exception as e:
                    st.session_state["_health_status"] = {
                        "ok": False, "msg": str(e),
                        "ts": time.strftime("%H:%M:%S"),
                    }
        hs = st.session_state.get("_health_status")
        if hs:
            if hs["ok"]:
                st.success(f"✅ {hs['msg']}  _(at {hs['ts']})_")
            else:
                st.error(f"❌ {hs['msg']}  _(at {hs['ts']})_")

        with st.expander("Backend env vars", expanded=False):
            for name in _llm.list_backends():
                d = _llm.describe_backend(name)
                bullet = (f"- **{name}** · key=`{d.get('api_key_env','')}`"
                            f"{'  · base=`'+d['base_url']+'`' if d.get('base_url') else ''}")
                st.markdown(bullet)
            if backend == "ollama" or backend == "(auto)":
                with st.spinner("Querying Ollama for installed models…"):
                    try:
                        models = _llm.list_ollama_models(timeout=2)
                    except Exception:
                        models = []
                if models:
                    st.markdown("**Local Ollama models** (paste a name above):")
                    for m in sorted(models, key=lambda r: r.get("name") or ""):
                        sz = (f" — {m['size_gb']:.1f} GB"
                              if m.get("size_gb") else "")
                        st.markdown(f"- `{m.get('name','')}`{sz}")
                else:
                    st.caption(
                        "_(no Ollama daemon reachable — start `ollama serve` "
                        "or set `OLLAMA_HOST_BASE`)_"
                    )

        st.divider()
        st.subheader("Run parameters")
        max_iter = st.slider("Max iterations", 1, 20, 8)
        max_tokens = st.slider("Max tokens / turn", 256, 16384, 4096, 256)
        temperature = st.slider("Temperature", 0.0, 1.5, 0.0, 0.05)

        st.divider()
        st.subheader("Tool subset")
        all_tools = [t.name for t in _tools.list_tools()]
        tool_filter = st.multiselect(
            "Restrict to specific tools",
            options=all_tools, default=[],
            help="Empty = all tools allowed.",
        )

        st.divider()
        if st.button("🗑 Clear conversation", use_container_width=True):
            st.session_state.messages = []
            st.rerun()
    return {
        "backend":     None if backend == "(auto)" else backend,
        "model":       model or None,
        "max_iter":    max_iter,
        "max_tokens":  max_tokens,
        "temperature": temperature,
        "tool_filter": tool_filter or None,
        "effective":   eff,
    }


# --------------------------- Artefact rendering ----------------------------

def _render_artefacts(paths: "list[str]") -> None:
    if not paths:
        return
    seen: "set[str]" = set()
    paths = [p for p in paths if not (p in seen or seen.add(p))]
    images, markdowns, csvs, jsons, others = [], [], [], [], []
    for p in paths:
        low = p.lower()
        if low.endswith((".png", ".svg", ".jpg", ".jpeg")):
            images.append(p)
        elif low.endswith(".md"):
            markdowns.append(p)
        elif low.endswith((".csv", ".tsv")):
            csvs.append(p)
        elif low.endswith(".json"):
            jsons.append(p)
        else:
            others.append(p)

    if images:
        cols = st.columns(min(2, len(images)))
        for idx, img in enumerate(images[:6]):
            try:
                cols[idx % len(cols)].image(
                    img, caption=Path(img).name, use_container_width=True,
                )
            except Exception as e:
                cols[idx % len(cols)].text(f"(image unavailable: {Path(img).name}) {e}")
        if len(images) > 6:
            st.caption(f"… and {len(images) - 6} more image(s) not shown")

    for md in markdowns[:6]:
        with st.expander(f"📄 {Path(md).name}", expanded=False):
            try:
                st.markdown(Path(md).read_text())
            except Exception as e:
                st.text(f"(could not read {md}: {e})")

    for csv in csvs[:6]:
        with st.expander(f"📊 {Path(csv).name}", expanded=False):
            try:
                import pandas as pd
                df = pd.read_csv(csv, nrows=400, sep=None, engine="python")
                st.dataframe(df, use_container_width=True, height=300)
                if len(df) >= 400:
                    st.caption("Preview limited to first 400 rows.")
            except Exception as e:
                st.text(f"(could not read {csv}: {e})")

    for j in jsons[:4]:
        with st.expander(f"📦 {Path(j).name}", expanded=False):
            st.code(j, language=None)
            try:
                txt = Path(j).read_text()
                if len(txt) < 30000:
                    st.json(txt)
                else:
                    st.caption(f"JSON >30KB; not inlined. Path: `{j}`")
            except Exception:
                pass

    for o in others[:8]:
        st.code(o)


# --------------------------- Event rendering -------------------------------

def _short_args(args: dict) -> str:
    if not args:
        return ""
    parts = []
    for k, v in list(args.items())[:5]:
        sv = str(v)
        if len(sv) > 36:
            sv = sv[:33] + "…"
        parts.append(f"`{k}`={sv}")
    return ", ".join(parts)


def _format_event_md(event: _agent.AgentEvent) -> str:
    k = event.kind
    p = event.payload
    if k == "run_start":
        return (f"▶ **run_start** — backend `{p.get('backend')}`, "
                f"model `{p.get('model')}`, {p.get('n_tools')} tools, "
                f"max {p.get('max_iterations')} iters")
    if k == "llm_call_start":
        return (f"🧠 **plan** — iter {p['iteration']}/"
                f"{p['max_iterations']}, {p['n_messages']} messages")
    if k == "llm_call_end":
        bits: "list[str]" = []
        if p.get("content"):
            preview = p['content'].strip().splitlines()[0][:160]
            bits.append(f"💭 {preview}")
        for tc in (p.get("tool_calls") or []):
            bits.append(f"→ will call **`{tc['name']}`** "
                        f"({_short_args(tc.get('arguments') or {})})")
        if p.get("usage"):
            u = p['usage']
            bits.append(f"_tokens in={u.get('input_tokens',0)} "
                        f"out={u.get('output_tokens',0)}_")
        return "\n\n".join(bits) or "_(empty turn)_"
    if k == "tool_call_start":
        return (f"🔧 **tool** `{p['name']}` "
                f"({_short_args(p.get('arguments') or {})})")
    if k == "tool_call_end":
        ec = p.get("exit_code")
        icon = "✅" if ec == 0 else "❌"
        artefacts = p.get("artifacts") or {}
        n = sum(len(v) for v in artefacts.values())
        a_summary = ""
        if artefacts:
            a_summary = "  ·  artefacts: " + ", ".join(
                f"{k}({len(v)})" for k, v in artefacts.items()
            )
        return f"   {icon} exit `{ec}`{a_summary}"
    if k == "error":
        return f"❌ **error** in `{p.get('where','?')}` — {p.get('error')}"
    if k == "final_answer":
        return ""        # rendered separately, full markdown
    if k == "run_end":
        return (f"✓ **run_end** — {p['iterations']} iters, "
                f"{p['tool_calls_made']} tool calls, stop "
                f"`{p['stop_reason']}`")
    return f"_{k}_  `{p}`"


# --------------------------- Main render -----------------------------------

def main() -> None:
    cfg = _sidebar()

    st.title("IGVFagent")
    st.caption(
        "Plan → Action → Results → Evaluation. Every tool call below also "
        "exists as an `igvfagent <skill>` shell command — the UI is just a "
        "different driver of the same skills."
    )

    # Persistent active-config banner so the user always sees what the
    # next chat run will hit. Shows red/amber/green based on whether the
    # configuration is plausibly runnable.
    eff = cfg.get("effective", {})
    banner_cols = st.columns([3, 3, 2, 2])
    banner_cols[0].markdown(f"**Backend** · `{eff.get('backend','?')}`")
    banner_cols[1].markdown(f"**Model** · `{eff.get('model','?')}`")
    banner_cols[2].markdown(
        "**Source** · " + str(eff.get("model_source", "?"))
    )
    if eff.get("key_env"):
        if eff["backend"] == "ollama":
            banner_cols[3].markdown("**Auth** · _none required_")
        elif eff.get("key_set"):
            banner_cols[3].markdown(
                f"**Auth** · ✅ `{eff['key_env']}` set"
            )
        else:
            banner_cols[3].error(f"⚠ `{eff['key_env']}` not set")
    st.divider()

    if "messages" not in st.session_state:
        st.session_state.messages = []

    # Replay prior conversation
    for entry in st.session_state.messages:
        with st.chat_message(entry["role"]):
            st.markdown(entry.get("content", ""))
            if entry.get("artefacts"):
                _render_artefacts(entry["artefacts"])
            if entry.get("meta"):
                st.caption(entry["meta"])

    # Suggestions on first load
    if not st.session_state.messages:
        st.info(
            "💡 Try one of these:\n"
            "- *Give me the comprehensive APOE evidence pack including "
            "literature corroboration and matching IGVF single-cell "
            "datasets.*\n"
            "- *Discover Parse SPLiT-seq datasets profiling macrophages "
            "in mouse and write the per-pool manifest.*\n"
            "- *Explain IGVFDS7013XXYV.*\n"
            "- *Validate APOE, TREM2, and LDLR against published "
            "Alzheimer / cardiovascular literature.*"
        )

    query = st.chat_input(
        "Ask IGVFagent — natural-language query (e.g. 'Walk the KG for APOE')"
    )
    if not query:
        return

    # Pre-flight validation: catch obvious mismatches before we even
    # call the agent. This is what bit you on query 2 (Anthropic
    # backend with a Qwen model name).
    pre_errors: "list[str]" = []
    pre_warnings: "list[str]" = []
    eff_b = (cfg.get("effective") or {}).get("backend") or ""
    eff_m = ((cfg.get("effective") or {}).get("model") or "").lower()
    if eff_b == "anthropic" and "claude" not in eff_m:
        pre_errors.append(
            f"Anthropic backend with non-Claude model `{eff_m}`. "
            f"Either change Backend to `(auto)` / `ollama`, or set Model "
            f"to a Claude name (e.g. `claude-sonnet-4-5`)."
        )
    elif eff_b == "openai" and not any(t in eff_m for t in
                                        ("gpt-", "o1-", "o3-", "o4-",
                                          "codex")):
        pre_warnings.append(
            f"OpenAI backend with model `{eff_m}` — that name doesn't "
            f"look like an OpenAI model. Common picks: `gpt-4o-mini`, "
            f"`gpt-4o`, `o1-preview`."
        )
    if eff_b == "anthropic" and not (cfg.get("effective") or {}).get("key_set"):
        pre_errors.append("`ANTHROPIC_API_KEY` is not set. Export it in "
                            "your shell or in the Compose `.env`.")
    if eff_b == "openai" and not (cfg.get("effective") or {}).get("key_set"):
        pre_errors.append("`OPENAI_API_KEY` is not set.")

    st.session_state.messages.append({"role": "user", "content": query})
    with st.chat_message("user"):
        st.markdown(query)

    if pre_errors:
        with st.chat_message("assistant"):
            st.error("Pre-flight check blocked the run:")
            for e in pre_errors:
                st.markdown(f"- {e}")
            st.markdown(
                "_Fix the configuration in the sidebar (it shows "
                "**Resolved configuration** under the inputs) and try "
                "again._"
            )
        st.session_state.messages.append({
            "role": "assistant",
            "content": ("**Pre-flight check failed:**\n\n"
                          + "\n".join("- " + e for e in pre_errors)),
        })
        return

    with st.chat_message("assistant"):
        for w in pre_warnings:
            st.warning(w)
        # Live event stream container
        status = st.status("Planning…", expanded=True)

        def cb(event: _agent.AgentEvent) -> None:
            line = _format_event_md(event)
            if line:
                status.write(line)
            if event.kind == "tool_call_start":
                status.update(label=f"🔧 {event.payload.get('name','')}",
                              state="running")
            elif event.kind == "llm_call_start":
                status.update(label=f"🧠 plan iter "
                                    f"{event.payload.get('iteration','?')}",
                              state="running")
            elif event.kind == "run_end":
                stop_r = event.payload.get("stop_reason")
                state = "complete" if stop_r == "complete" else "error"
                status.update(label=f"Done · {event.payload.get('iterations')} "
                                    f"iters, "
                                    f"{event.payload.get('tool_calls_made')} "
                                    f"tool calls · stop `{stop_r}`",
                              state=state,
                              expanded=state != "complete")

        try:
            result = _agent.run(
                query,
                backend=cfg["backend"],
                model=cfg["model"],
                max_iterations=cfg["max_iter"],
                max_tokens=cfg["max_tokens"],
                temperature=cfg["temperature"],
                tools_subset=cfg["tool_filter"],
                callback=cb,
            )
        except RuntimeError as exc:
            st.error(str(exc))
            st.markdown(
                "**Hint** — install an LLM SDK and credentials:\n\n"
                "```\n"
                "pip install 'igvfagent[llm]'\n"
                "ollama serve   # in another terminal\n"
                "ollama pull qwen3:8b\n"
                "```\n"
                "Or set `ANTHROPIC_API_KEY` / `OPENAI_API_KEY` and pick "
                "the matching backend in the sidebar."
            )
            return
        except Exception as exc:  # pylint: disable=broad-except
            import traceback
            st.error(f"Agent run failed: {exc}")
            with st.expander("Show full traceback", expanded=False):
                st.code(traceback.format_exc())
            return

        # Final answer (or error). On error stops, render as st.error so
        # the message is impossible to miss; otherwise render markdown.
        if result.stop_reason != "complete" and result.final_answer:
            st.error("Agent run ended before completion. See details below.")
            st.markdown(result.final_answer)
        elif result.final_answer:
            st.markdown(result.final_answer)
        else:
            st.warning("_The agent finished without producing a final answer._")

        # Inline artefact rendering
        if result.artefacts:
            with st.expander(f"📁 Artefacts ({len(result.artefacts)})",
                             expanded=True):
                _render_artefacts(result.artefacts)

        meta_caption = (
            f"backend `{result.backend}`  ·  model `{result.model}`  ·  "
            f"{result.iterations} iter · {result.tool_calls_made} tool calls "
            f"· stop `{result.stop_reason}`"
        )
        if result.report_path:
            meta_caption += f"  ·  report `{result.report_path}`"
        st.caption(meta_caption)

    st.session_state.messages.append({
        "role": "assistant",
        "content": result.final_answer,
        "artefacts": result.artefacts,
        "meta": meta_caption,
    })


def _diagnostics_panel() -> None:
    """Static debug info shown at the bottom of the page so users can
    paste it verbatim when reporting issues."""
    import os
    import platform
    with st.expander("🛠 Diagnostics (paste this verbatim if you hit a bug)",
                       expanded=False):
        rows = [
            ("igvfagent", __version__),
            ("streamlit", getattr(st, "__version__", "?")),
            ("python", f"{platform.python_version()} ({platform.machine()})"),
            ("platform", platform.platform()),
            ("IGVF_PROJECT_ROOT", os.environ.get("IGVF_PROJECT_ROOT", "(unset)")),
            ("IGVF_LLM_BACKEND",  os.environ.get("IGVF_LLM_BACKEND", "(unset)")),
            ("IGVF_LLM_MODEL",    os.environ.get("IGVF_LLM_MODEL", "(unset)")),
            ("OLLAMA_HOST_BASE",  os.environ.get("OLLAMA_HOST_BASE", "(unset)")),
        ]
        for env in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GROQ_API_KEY",
                    "TOGETHER_API_KEY", "DEEPINFRA_API_KEY", "HF_TOKEN"):
            rows.append((env, "✅ set" if os.environ.get(env) else "(unset)"))
        for sdk in ("anthropic", "openai", "pyfaidx", "pyBigWig", "cooler",
                    "hicstraw"):
            try:
                mod = __import__(sdk)
                ver = getattr(mod, "__version__", "?")
                rows.append((f"sdk: {sdk}", ver))
            except Exception:
                rows.append((f"sdk: {sdk}", "(not installed)"))
        rows.append(("registered tools", str(len(_tools.list_tools()))))
        try:
            ollama_models = _llm.list_ollama_models(timeout=2)
            if ollama_models:
                rows.append(("ollama models",
                              ", ".join(m.get("name", "")
                                          for m in ollama_models)[:200]))
            else:
                rows.append(("ollama models", "(daemon unreachable)"))
        except Exception as e:
            rows.append(("ollama models", f"(error: {e})"))
        st.code("\n".join(f"{k:24} {v}" for k, v in rows))


# Streamlit runs `streamlit run <file>` with __name__ == "__main__", so a
# single guard is sufficient. The module is otherwise safe to import (no
# side effects beyond the `st.set_page_config` call at module top, which
# is idempotent on re-imports).
if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # pylint: disable=broad-except
        import traceback
        st.error(f"Streamlit page raised an exception: {e}")
        with st.expander("Full traceback", expanded=True):
            st.code(traceback.format_exc())
    finally:
        try:
            _diagnostics_panel()
        except Exception:
            pass
