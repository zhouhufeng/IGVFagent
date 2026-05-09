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

import os
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


_BACKEND_KIND_LABELS = {
    "local":     "🖥  Local LLM (Ollama)",
    "anthropic": "🤖  Anthropic Claude API",
    "openai":    "⚡  OpenAI / Codex API",
    "advanced":  "🔧  Other (advanced)",
}

_BACKEND_KIND_TO_NAME = {
    "local":     "ollama",
    "anthropic": "anthropic",
    "openai":    "openai",
    # advanced -> resolved from sub-selectbox below
}


def _sidebar_backend_kind() -> str:
    """Render the top-level backend-type radio. Resets dependent state
    when the kind changes."""
    options = list(_BACKEND_KIND_LABELS.keys())
    labels = [_BACKEND_KIND_LABELS[k] for k in options]
    prev = st.session_state.get("_backend_kind", "local")
    idx = options.index(prev) if prev in options else 0
    chosen_label = st.radio("Backend type", labels, index=idx,
                              key="_backend_kind_radio",
                              help="Pick the model provider. Local LLM uses "
                                   "your Ollama daemon; the cloud options need "
                                   "their respective API keys.")
    chosen = options[labels.index(chosen_label)]
    if chosen != prev:
        # Reset model selection when the user switches kind
        for k in ("_local_model_choice", "_anthropic_model_choice",
                   "_openai_model_choice", "_advanced_backend",
                   "_advanced_model"):
            st.session_state.pop(k, None)
        st.session_state["_backend_kind"] = chosen
    return chosen


def _sidebar_local_model_picker() -> "tuple[str, str]":
    """Render the installed-models dropdown + downloadable picker for
    Local LLM. Returns (backend_name, model_name)."""
    with st.spinner("Querying local Ollama daemon…"):
        try:
            installed = _llm.list_ollama_models(timeout=3)
        except Exception:
            installed = []

    if not installed:
        st.warning(
            "Local Ollama daemon not reachable on "
            f"`{_llm._BACKENDS['ollama']['base_url']}`.\n\n"
            "Start it with `ollama serve` in another terminal, or set "
            "`OLLAMA_HOST_BASE` to the right URL."
        )
        st.caption("After Ollama is up, click **Refresh** to repopulate.")
        if st.button("🔄 Refresh installed models", use_container_width=True):
            st.rerun()
        installed_names: "list[str]" = []
    else:
        installed_names = sorted(m.get("name", "") for m in installed
                                    if m.get("name"))

    # Installed models dropdown
    if installed_names:
        labels = []
        for m in sorted(installed, key=lambda r: r.get("name") or ""):
            sz = (f" · {m['size_gb']:.1f} GB"
                  if m.get("size_gb") else "")
            fam = (f" · {m['family']}" if m.get("family") else "")
            labels.append(f"{m['name']}{sz}{fam}")
        prev = st.session_state.get("_local_model_choice")
        idx = next((i for i, n in enumerate(installed_names) if n == prev),
                   0)
        chosen_label = st.selectbox(
            "Installed Ollama models",
            labels, index=idx, key="_local_installed_select",
            help="Models the local Ollama daemon already has on disk.",
        )
        chosen_idx = labels.index(chosen_label)
        chosen_model = installed_names[chosen_idx]
        st.session_state["_local_model_choice"] = chosen_model
    else:
        chosen_model = ""

    # Refresh button to re-query the daemon
    cols = st.columns([2, 1])
    with cols[1]:
        if st.button("🔄", help="Re-query the Ollama daemon",
                       use_container_width=True):
            st.rerun()
    with cols[0]:
        st.caption(f"Endpoint: `{_llm._BACKENDS['ollama']['base_url']}`")

    # Downloadable model picker
    with st.expander("⬇ Download more models", expanded=False):
        st.caption(
            "Pick a model from Ollama's library and pull it into your "
            "local daemon. Big models (≥ 20 GB) can take many minutes — "
            "for those you may prefer running `ollama pull <name>` in a "
            "terminal directly."
        )
        installed_set = set(installed_names)
        available = [(name, gb, note)
                     for name, gb, note in _llm.OLLAMA_LIBRARY
                     if name not in installed_set]
        if not available:
            st.caption("_All curated models already installed._")
        else:
            avail_labels = [f"{n} — ~{gb:.1f} GB ({note})"
                             for n, gb, note in available]
            pick_label = st.selectbox(
                "Available to download", avail_labels,
                key="_download_pick_select",
            )
            pick_name = available[avail_labels.index(pick_label)][0]
            if st.button(f"⬇ Pull `{pick_name}`",
                            use_container_width=True,
                            key="_download_pull_btn"):
                _do_ollama_pull(pick_name)

    return "ollama", chosen_model


def _do_ollama_pull(model_name: str) -> None:
    """Stream-pull an Ollama model with a Streamlit progress bar."""
    progress = st.progress(0.0, text=f"Pulling `{model_name}` …")
    status_box = st.empty()
    last_status = ""
    try:
        for status, pct, total, completed, errored in \
                _llm.pull_ollama_model(model_name):
            if errored:
                st.error(f"Pull failed: {status}")
                break
            last_status = status or last_status
            label = (f"`{model_name}` · {last_status} · "
                     f"{completed/1e9:.2f}/{total/1e9:.2f} GB"
                     if total else f"`{model_name}` · {last_status}")
            try:
                progress.progress(min(max(pct / 100.0, 0.0), 1.0),
                                    text=label)
            except Exception:
                pass
            status_box.caption(label)
        else:
            progress.progress(1.0, text=f"`{model_name}` pulled.")
            st.success(f"✅ `{model_name}` is now available locally.")
            time.sleep(0.5)
            st.rerun()
    except Exception as e:
        st.error(f"Pull failed: {e}")


def _sidebar_anthropic_model_picker() -> "tuple[str, str]":
    """Curated Claude model dropdown + key check."""
    options = list(_llm.ANTHROPIC_MODELS) + ["(custom...)"]
    prev = st.session_state.get("_anthropic_model_choice", options[1])
    idx = options.index(prev) if prev in options else 1
    chosen = st.selectbox("Claude model", options, index=idx,
                            key="_anthropic_model_select")
    if chosen == "(custom...)":
        chosen = st.text_input("Custom Claude model id",
                                 value="", key="_anthropic_model_custom")
    st.session_state["_anthropic_model_choice"] = chosen

    if not os.environ.get("ANTHROPIC_API_KEY"):
        st.error("`ANTHROPIC_API_KEY` is not set. Export it in your shell, "
                  "then restart the UI.")
    else:
        st.caption("✅ `ANTHROPIC_API_KEY` is set.")
    return "anthropic", chosen


def _sidebar_openai_model_picker() -> "tuple[str, str]":
    """Curated OpenAI / Codex model dropdown + key check."""
    options = list(_llm.OPENAI_MODELS) + ["(custom...)"]
    prev = st.session_state.get("_openai_model_choice", options[0])
    idx = options.index(prev) if prev in options else 0
    chosen = st.selectbox("OpenAI / Codex model", options, index=idx,
                            key="_openai_model_select")
    if chosen == "(custom...)":
        chosen = st.text_input("Custom OpenAI model id",
                                 value="", key="_openai_model_custom")
    st.session_state["_openai_model_choice"] = chosen

    if not os.environ.get("OPENAI_API_KEY"):
        st.error("`OPENAI_API_KEY` is not set. Export it in your shell, "
                  "then restart the UI.")
    else:
        st.caption("✅ `OPENAI_API_KEY` is set.")
    return "openai", chosen


def _sidebar_advanced_picker() -> "tuple[str, str]":
    """Free-form backend + model picker for vLLM / TGI / Groq / etc."""
    backends = [b for b in _llm.list_backends()
                 if b not in ("anthropic", "openai", "ollama")]
    backend = st.selectbox("Backend", backends,
                              key="_advanced_backend_select")
    desc = _llm.describe_backend(backend)
    if desc.get("base_url"):
        st.caption(f"Endpoint: `{desc['base_url']}`")
    if desc.get("api_key_env"):
        if os.environ.get(desc["api_key_env"]):
            st.caption(f"✅ `{desc['api_key_env']}` is set.")
        else:
            st.error(f"`{desc['api_key_env']}` is not set.")
    model = st.text_input("Model", value="", key="_advanced_model_input")
    return backend, model


def _sidebar_load_button(backend: str, model: str) -> None:
    """Renders the Load Model button + persistent status."""
    cols = st.columns([2, 1])
    with cols[0]:
        clicked = st.button("⚙  Load Model", type="primary",
                              use_container_width=True,
                              disabled=not (backend and model),
                              help="Ping the backend to confirm credentials, "
                                   "model presence, and (for Ollama) preload "
                                   "the weights into memory so the first chat "
                                   "doesn't pay the cold-load cost.")
    with cols[1]:
        if st.button("🔄", help="Refresh state",
                       use_container_width=True, key="_load_refresh_btn"):
            st.session_state.pop("_loaded_status", None)
            st.rerun()
    if clicked:
        with st.spinner(f"Loading `{model}` via `{backend}` …"):
            t0 = time.time()
            try:
                msg = _llm.chat(
                    messages=[{"role": "user", "content": "READY"}],
                    backend=backend, model=model,
                    max_tokens=8, temperature=0.0,
                )
                dt = time.time() - t0
                st.session_state["_loaded_status"] = {
                    "ok":      True,
                    "backend": msg.backend,
                    "model":   msg.model,
                    "ts":      time.strftime("%H:%M:%S"),
                    "secs":    dt,
                    "preview": (msg.content or "").strip()[:60],
                }
            except Exception as e:
                st.session_state["_loaded_status"] = {
                    "ok":      False,
                    "backend": backend, "model": model,
                    "ts":      time.strftime("%H:%M:%S"),
                    "error":   str(e),
                }
    status = st.session_state.get("_loaded_status")
    if status:
        if status["ok"]:
            st.success(
                f"✅ Loaded · `{status['model']}` via `{status['backend']}`  "
                f"_(in {status['secs']:.1f}s, at {status['ts']})_"
            )
        else:
            st.error(
                f"❌ Load failed for `{status['model']}` via "
                f"`{status['backend']}`  _(at {status['ts']})_\n\n"
                f"`{status['error']}`"
            )


def _sidebar() -> dict:
    with st.sidebar:
        st.markdown(f"## 🧬 IGVFagent\n_v{__version__}_")
        st.caption(
            "Natural-language interface to the IGVF / ENCODE single-cell, "
            "variant, regulatory-element, and literature stack."
        )
        st.divider()

        st.subheader("Model")
        kind = _sidebar_backend_kind()
        if kind == "local":
            backend, model = _sidebar_local_model_picker()
        elif kind == "anthropic":
            backend, model = _sidebar_anthropic_model_picker()
        elif kind == "openai":
            backend, model = _sidebar_openai_model_picker()
        else:
            backend, model = _sidebar_advanced_picker()

        _sidebar_load_button(backend, model)

        # Resolved configuration block — kept for transparency.
        eff = _resolve_effective_config(backend, model)
        with st.expander("Resolved configuration", expanded=False):
            st.markdown(
                f"- Backend: `{eff['backend']}`\n"
                f"- Model: `{eff['model']}`  _({eff['model_source']})_"
            )
            if eff["key_env"]:
                icon = ("✅" if eff["key_set"] else
                        ("➖" if eff["backend"] == "ollama" else "❌"))
                st.markdown(f"- Credential: `{eff['key_env']}` {icon}")
            if eff["base_url"]:
                st.markdown(f"- Endpoint: `{eff['base_url']}`")

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
            help="Empty = all tools allowed. Smaller subsets cut prompt "
                 "size and noticeably speed up local LLMs.",
        )

        st.divider()
        if st.button("🗑 Clear conversation", use_container_width=True):
            st.session_state.messages = []
            st.rerun()

    return {
        "backend":      backend or None,
        "model":        model or None,
        "max_iter":     max_iter,
        "max_tokens":   max_tokens,
        "temperature":  temperature,
        "tool_filter":  tool_filter or None,
        "effective":    eff,
        "kind":         kind,
        "loaded_ok":    bool(st.session_state.get("_loaded_status", {}).get("ok")),
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

    # Active-model banner — LM Studio-style. Big, unmissable, colored
    # by whether a model is loaded.
    status = st.session_state.get("_loaded_status")
    eff = cfg.get("effective", {})
    kind = cfg.get("kind", "local")
    kind_label = _BACKEND_KIND_LABELS.get(kind, kind)
    if status and status.get("ok"):
        st.success(
            f"### {kind_label}  ·  `{status['model']}`  "
            f"·  ✅ Loaded  "
            f"_(in {status['secs']:.1f}s)_\n\n"
            f"All chat queries below will use this model."
        )
    elif status and not status.get("ok"):
        st.error(
            f"### {kind_label}  ·  `{status.get('model','?')}`  "
            f"·  ❌ Load failed\n\n"
            f"Fix the configuration in the sidebar and click "
            f"**⚙  Load Model** again. Error: `{status.get('error','?')}`"
        )
    else:
        st.warning(
            f"### {kind_label}  ·  No model loaded yet\n\n"
            f"Pick a model in the sidebar and click **⚙  Load Model** "
            f"before sending a chat. Until then queries will fail "
            f"or be slow on first call."
        )
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
