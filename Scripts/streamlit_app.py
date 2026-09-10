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

import base64
import os
import re
import sys
import time
import json
from pathlib import Path
from typing import Any

import streamlit as st

# Dual-mode import (installed package OR running from a checkout).
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from igvfagent import _agent, _llm, _pathguard, _tools, _userext, __version__
    from igvfagent import load_dotenv as _load_dotenv
    from igvfagent._stcompat import fit
except Exception:
    import _agent  # type: ignore
    import _llm    # type: ignore
    import _pathguard  # type: ignore
    import _tools  # type: ignore
    import _userext  # type: ignore
    from _stcompat import fit  # type: ignore
    try:
        from __init__ import __version__, load_dotenv as _load_dotenv  # type: ignore
    except Exception:
        __version__ = "0.2.9"
        _load_dotenv = None  # type: ignore

# Belt-and-suspenders: ensure the repo-root .env is loaded into THIS
# Streamlit process even when launched via `streamlit run` directly
# (which bypasses cli.py). Real env vars always win. Without this, the
# Anthropic/OpenAI key silently appears "not set" in the model picker.
if _load_dotenv is not None:
    _load_dotenv()

# KG visualizer — optional, soft-fail if matplotlib/networkx missing.
try:
    from igvfagent import kg_visualizer as _kgviz  # type: ignore
except Exception:
    try:
        import kg_visualizer as _kgviz  # type: ignore
    except Exception:
        _kgviz = None

# Single-cell visualizer — optional, needs scanpy + anndata.
try:
    from igvfagent import sc_visualizer as _scviz  # type: ignore
except Exception:
    try:
        import sc_visualizer as _scviz  # type: ignore
    except Exception:
        _scviz = None

# Network visualizer — optional, needs networkx + matplotlib + pyvis.
try:
    from igvfagent import network_visualizer as _nwviz  # type: ignore
except Exception:
    try:
        import network_visualizer as _nwviz  # type: ignore
    except Exception:
        _nwviz = None

# Spatial-ATAC-Hi-C browser — optional, needs numpy + matplotlib.
try:
    from igvfagent import spatial_hic_visualizer as _sphic  # type: ignore
except Exception:
    try:
        import spatial_hic_visualizer as _sphic  # type: ignore
    except Exception:
        _sphic = None

# Benchmark figure gallery — stdlib only (reads committed PNG/SVG off disk),
# so unlike the other visualizers it has no optional scientific deps.
try:
    from igvfagent import benchmark_visualizer as _bmviz  # type: ignore
except Exception:
    try:
        import benchmark_visualizer as _bmviz  # type: ignore
    except Exception:
        _bmviz = None


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
    "local":      "🖥  Local LLM (Ollama)",
    "anthropic":  "🤖  Anthropic Claude API",
    "openai":     "⚡  OpenAI / Codex API",
    "claude_cli": "🧠  Claude Code CLI (subprocess)",
    "codex_cli":  "💻  Codex CLI (subprocess)",
    "advanced":   "🔧  Other (advanced)",
}

_BACKEND_KIND_TO_NAME = {
    "local":      "ollama",
    "anthropic":  "anthropic",
    "openai":     "openai",
    "claude_cli": "claude_cli",
    "codex_cli":  "codex_cli",
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
        if st.button("🔄 Refresh installed models", **fit(st.button)):
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
                       **fit(st.button)):
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
                            **fit(st.button),
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
    # Preselect the backend's own default rather than a hardcoded index, so
    # reordering ANTHROPIC_MODELS can't silently open the picker on a model
    # other than the one _DEFAULT_MODELS says the backend will actually use.
    default_model = _llm._DEFAULT_MODELS.get("anthropic", "")
    fallback_idx = options.index(default_model) if default_model in options else 0
    prev = st.session_state.get("_anthropic_model_choice", options[fallback_idx])
    idx = options.index(prev) if prev in options else fallback_idx
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


def _sidebar_claude_cli_picker() -> "tuple[str, str]":
    """Picker for the Claude Code CLI backend."""
    ok, info = _llm.claude_cli_available()
    if ok:
        st.caption(f"✅ Claude Code CLI detected — {info}")
    else:
        st.error(
            f"Claude Code CLI not available: {info}\n\n"
            "Install with `npm i -g @anthropic-ai/claude-code`, then "
            "log in via `claude login`. Restart this UI afterward."
        )
        return "claude_cli", ""

    st.caption(
        "Subprocess shell-out to `claude --print` for each LLM turn. "
        "Reuses your Claude Code login (no separate ANTHROPIC_API_KEY "
        "needed). Trade-offs: 5–15s per turn of CLI overhead and "
        "tool-calls are XML-parsed (less robust than native function "
        "calling). For best speed, prefer the **Anthropic Claude API** "
        "backend instead with a real API key."
    )

    # Claude Code CLI runs against whichever models the local `claude`
    # binary supports. Limit the picker to the current Claude 5 tiers so
    # users do not pick a retired/superseded model id. Older ids still
    # work via "(custom...)" if you need to pin one.
    _CLAUDE_CLI_MODELS = (
        "claude-opus-5",
        "claude-sonnet-5",
        "claude-fable-5-1",
    )
    options = ["(use Claude Code's configured default)"] + \
              list(_CLAUDE_CLI_MODELS) + ["(custom...)"]
    prev = st.session_state.get("_claude_cli_model_choice", options[0])
    idx = options.index(prev) if prev in options else 0
    chosen = st.selectbox("Claude model override", options, index=idx,
                            key="_claude_cli_model_select",
                            help="Leave on the default to use whatever "
                                 "Claude Code is configured for, or pick "
                                 "an explicit Claude model name.")
    if chosen == "(use Claude Code's configured default)":
        chosen = ""
    elif chosen == "(custom...)":
        chosen = st.text_input("Custom Claude model id",
                                 value="", key="_claude_cli_model_custom")
    st.session_state["_claude_cli_model_choice"] = chosen or options[0]
    return "claude_cli", chosen or "(default)"


def _sidebar_codex_cli_picker() -> "tuple[str, str]":
    """Picker for the Codex CLI backend."""
    ok, info = _llm.codex_cli_available()
    if ok:
        st.caption(f"✅ Codex CLI detected — {info}")
    else:
        st.error(
            f"Codex CLI not available: {info}\n\n"
            "Install with `npm i -g @openai/codex`, then sign in via "
            "`codex login`. Restart this UI afterward."
        )
        return "codex_cli", ""

    st.caption(
        "Subprocess shell-out to `codex exec` for each LLM turn. "
        "Reuses your Codex CLI login (no separate OPENAI_API_KEY "
        "needed). Same trade-offs as the Claude Code CLI backend: "
        "5–15s per turn of CLI overhead and tool-calls are XML-parsed "
        "(less robust than native function calling). Prefer the "
        "**OpenAI / Codex API** backend with a real key for best speed."
    )

    options = ["(use Codex CLI's configured default)"] + \
              list(_llm.OPENAI_MODELS) + ["(custom...)"]
    prev = st.session_state.get("_codex_cli_model_choice", options[0])
    idx = options.index(prev) if prev in options else 0
    chosen = st.selectbox("Model override", options, index=idx,
                            key="_codex_cli_model_select",
                            help="Leave on the default to use whatever "
                                 "Codex CLI is configured for, or pick "
                                 "an explicit model name.")
    if chosen == "(use Codex CLI's configured default)":
        chosen = ""
    elif chosen == "(custom...)":
        chosen = st.text_input("Custom model id",
                                 value="", key="_codex_cli_model_custom")
    st.session_state["_codex_cli_model_choice"] = chosen or options[0]
    return "codex_cli", chosen or "(default)"


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
                              **fit(st.button),
                              disabled=not (backend and model),
                              help="Ping the backend to confirm credentials, "
                                   "model presence, and (for Ollama) preload "
                                   "the weights into memory so the first chat "
                                   "doesn't pay the cold-load cost.")
    with cols[1]:
        if st.button("🔄", help="Refresh state",
                       **fit(st.button), key="_load_refresh_btn"):
            st.session_state.pop("_loaded_status", None)
            st.rerun()
    if clicked:
        with st.spinner(f"Loading `{model}` via `{backend}` …"):
            t0 = time.time()
            try:
                # 256, not 8: Claude 5 models (Opus 5, Sonnet 5, Fable 5)
                # have thinking on by default, and max_tokens caps thinking
                # + response text together. A tiny budget is spent entirely
                # on thinking, so the probe returned stop_reason=max_tokens
                # with empty content and the status preview came up blank.
                msg = _llm.chat(
                    messages=[{"role": "user", "content": "READY"}],
                    backend=backend, model=model,
                    max_tokens=256, temperature=0.0,
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


def _sidebar_document_upload() -> None:
    """Upload a manuscript or figure and make it available to the agent.

    Uploaded files land in the workspace so every skill can reach them by
    path, which is what lets a user say "reproduce the analysis in this
    paper": `document plan` turns the PDF into accessions, assays and a bench
    chain, and the agent takes it from there.
    """
    st.subheader("📄 Documents")
    st.caption(
        "Upload a manuscript (PDF/DOCX/text) and ask IGVFagent to reproduce "
        "its analysis, or upload figures to discuss. Files are stored in the "
        "workspace and referenced by path."
    )
    ups = st.file_uploader(
        "Upload documents or images",
        type=["pdf", "docx", "txt", "md", "csv", "tsv",
              "png", "jpg", "jpeg", "gif", "svg"],
        accept_multiple_files=True,
        key="_doc_uploader",
        help="PDFs and DOCX are text-extracted for accession/assay detection. "
             "Images are stored and viewable but not OCR'd.",
    )
    if ups and st.button("📥 Add to workspace", **fit(st.button),
                          key="_doc_add"):
        dest_dir = _upload_dir()
        saved = []
        for up in ups:
            # Basename only: an uploaded filename is untrusted input and must
            # never be able to escape the upload directory.
            dest = dest_dir / Path(up.name).name
            dest.write_bytes(up.getbuffer())
            saved.append(dest)
        st.session_state["_uploaded_docs"] = [
            str(p) for p in sorted(dest_dir.iterdir()) if p.is_file()]
        st.success(f"Added {len(saved)} file(s) to this session's "
                    f"workspace (`{dest_dir.relative_to(_PROJECT_ROOT)}`).")

    docs = st.session_state.get("_uploaded_docs") or []
    if not docs:
        # This session's own directory ONLY. Reading the shared parent is
        # what showed every visitor everyone else's uploads before they had
        # asked anything.
        d = _upload_dir()
        if d.is_dir():
            docs = [str(p) for p in sorted(d.iterdir()) if p.is_file()]
            st.session_state["_uploaded_docs"] = docs
    if docs:
        with st.expander(f"In workspace ({len(docs)})", expanded=False):
            for f in docs[:25]:
                st.markdown(f"- `{Path(f).name}`")
            st.caption("Ask, for example: _\"Read the uploaded paper and "
                       "reproduce its analysis\"_ — the agent calls "
                       "`document plan` on it.")


def _session_token() -> str:
    """Stable id for this browser session, used to scope uploads.

    Uploads all landed in one flat Data/Uploads/, with two consequences on a
    shared deployment. Every session listed every file, so a visitor saw the
    filenames of everyone else's uploads and the agent could read them. And
    the destination was the bare basename, so a second visitor uploading
    "manifest.csv" silently OVERWROTE the first -- after which the first
    visitor's next question analysed somebody else's table with nothing to
    indicate it. The second is a correctness failure, not only a privacy one.
    """
    tok = st.session_state.get("_session_token")
    if not tok:
        import secrets
        tok = f"{time.strftime('%Y%m%d-%H%M%S')}-{secrets.token_hex(4)}"
        st.session_state["_session_token"] = tok
    return tok


def _upload_dir() -> Path:
    d = _PROJECT_ROOT / "Data" / "Uploads" / _session_token()
    d.mkdir(parents=True, exist_ok=True)
    return d


def _extension_upload_allowed() -> "tuple[bool, str]":
    """May this deployment accept uploaded tool manifests and skill modules?

    An uploaded .py lands in ~/.igvfagent/skills/ and refresh_user_tools()
    registers it immediately, after which the agent runs it as a subprocess.
    That is code execution inside the container -- which also holds the IGVF
    Portal key -- granted to anyone who can reach the page.

    The panel was rendered unconditionally, so on the hosted deployment
    knowing the shared password was enough. IGVF_ALLOW_AGENT_AUTHORING was
    not consulted here at all; it gates the agent authoring its own tools,
    not a visitor uploading one.

    On a shared deployment this now requires a SEPARATE, deliberate opt-in
    (IGVF_ALLOW_UPLOAD_EXTENSIONS=1). Defaulting it off closes the hole
    without a container rebuild, and a single-user local run is unaffected.
    """
    if not _public_mode():
        return True, ""
    if os.environ.get("IGVF_ALLOW_UPLOAD_EXTENSIONS", "0") == "1":
        return True, ""
    return False, (
        "Uploading tool manifests or skill modules is disabled on this "
        "shared deployment. An uploaded `.py` becomes code this server "
        "executes, so it is enabled only where the operator has set "
        "`IGVF_ALLOW_UPLOAD_EXTENSIONS=1`. Extensions can still be added "
        "locally, or installed on the host by the operator.")


def _sidebar_user_extensions() -> None:
    """User-extension panel: show discovered custom tools/skills and let
    the user install new ones from the browser (saved to ~/.igvfagent/,
    absorbed into the registry immediately — no restart)."""
    st.subheader("🧩 User extensions")
    st.caption(
        "Bring your own tools (YAML/JSON manifest wrapping any script) and "
        "skills (Python modules) — no core-code edits. Templates in "
        "`Docs/Examples/user_extensions/`; tutorial: README → "
        "*Extending IGVFagent*."
    )
    notice = st.session_state.pop("_ext_notice", None)
    if notice:
        st.success(notice)

    user_tools = _userext.discover_tools()
    user_skills = _userext.discover_skills()
    label = (f"Installed: {len(user_tools)} tool(s) · "
             f"{len(user_skills)} skill(s)")
    with st.expander(label, expanded=False):
        st.markdown("**Search locations** (first definition wins):")
        for d in _userext.extension_dirs():
            mark = "✅" if d.is_dir() else "➖"
            st.markdown(f"- {mark} `{d}`")
        if user_tools:
            st.markdown("**Custom tools** — callable by the agent and "
                        "listed in the tool picker above:")
            for t in user_tools:
                st.markdown(f"- `{t['name']}` — {t['description']}")
        if user_skills:
            st.markdown("**Custom skills** — run as `igvfagent <name>` "
                        "in a terminal:")
            for name, entry in user_skills.items():
                st.markdown(f"- `{name}` — {entry['description']}")
        problems = _userext.problems()
        if problems:
            st.markdown("**Skipped definitions:**")
            for p in problems:
                st.warning(p, icon="⚠️")

    allowed, why = _extension_upload_allowed()
    if not allowed:
        st.caption(why)
        return
    uploads = st.file_uploader(
        "Add extension files",
        type=["yaml", "yml", "json", "py"],
        accept_multiple_files=True,
        help="Tool manifests (.yaml/.json) and skill modules (.py). "
             "Installed into ~/.igvfagent/ and picked up immediately — "
             "no restart needed.",
    )
    if uploads and st.button("📦 Install extensions",
                              **fit(st.button)):
        base = Path.home() / ".igvfagent"
        for up in uploads:
            sub = "skills" if up.name.endswith(".py") else "tools"
            dest = base / sub / Path(up.name).name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(up.getbuffer())
        added = _tools.refresh_user_tools()
        st.session_state["_ext_notice"] = (
            f"Installed {len(uploads)} file(s) into `{base}` — "
            f"{added} new tool(s) registered with the agent."
        )
        st.rerun()


def _public_mode() -> bool:
    """True on a shared/hosted deployment (``IGVF_PUBLIC_MODE=1``).

    Hosted at igvfagent.genohub.org the LLM runs on the *operator's* API
    key, so the free-form backend picker stops being a convenience: it lets
    any visitor point the app at Ollama, which isn't running there, and get
    a confusing connection error. Public mode fixes the backend and offers
    a curated **allowlist** of models the operator is willing to pay for —
    visitors still choose, but only from that list. Local single-user runs
    are unaffected.
    """
    return os.environ.get("IGVF_PUBLIC_MODE", "").strip().lower() in ("1", "true", "yes", "on")


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, "").strip() or default)
    except ValueError:
        return default


# Models offered on a shared deployment, in menu order. Override with
# IGVF_PUBLIC_MODELS as a comma-separated list of ids (unknown ids are shown
# with the id as their label, so a new model can be enabled without a code
# change). Prices are per million tokens, input/output, for the operator's
# own cost awareness — the visitor is spending someone else's money.
# Ordered cheapest-first so the default sits at the top and the premium
# models are a deliberate scroll, not an accidental click — on a shared key
# the difference between Sonnet 5 and Fable 5 is ~3.3x per token.
_PUBLIC_MODEL_CATALOG = {
    "claude-sonnet-5":  ("Sonnet 5 — default, fast",      "$3 / $15 per Mtok · 1M context"),
    "claude-haiku-4-5": ("Haiku 4.5 — fastest, cheapest", "$1 / $5 per Mtok · 200K context"),
    "claude-opus-5":    ("Opus 5 — more capable",         "$5 / $25 per Mtok · 1M context · premium"),
    # Fable 5.1, id confirmed against GET /v1/models. The per-token price is
    # deliberately NOT copied over from Fable 5's $10/$50: the models endpoint
    # does not publish pricing and a wrong number here is worse than none,
    # because this line exists to make an expensive click deliberate. Fill it
    # in from the console when known.
    "claude-fable-5-1": ("Fable 5.1 — deepest reasoning",
                          "1M context · most expensive tier"),
}
_PUBLIC_MODELS_DEFAULT = list(_PUBLIC_MODEL_CATALOG)


def _public_model_choices() -> "list[str]":
    raw = os.environ.get("IGVF_PUBLIC_MODELS", "").strip()
    if not raw:
        return list(_PUBLIC_MODELS_DEFAULT)
    picked = [m.strip() for m in raw.split(",") if m.strip()]
    return picked or list(_PUBLIC_MODELS_DEFAULT)


_ORCHESTRATORS = {
    "internal":   ("⚙  Internal (IGVFagent ReAct)", None),
    "claude_cli": ("🧠  External — Claude Code CLI", "claude"),
    "codex_cli":  ("💻  External — Codex CLI",       "codex"),
}

# Which orchestrators a deployment offers. The hosted instance is an
# Anthropic-only test bed, so it ships "internal,claude_cli" and never shows
# Codex — an option that can only ever fail is worse than no option. Local
# installs keep the full list.
_DEFAULT_PUBLIC_ORCHESTRATORS = "internal,claude_cli"


def _allowed_orchestrators() -> "list[str]":
    raw = os.environ.get("IGVF_PUBLIC_ORCHESTRATORS", "").strip()
    if not raw:
        raw = (_DEFAULT_PUBLIC_ORCHESTRATORS if _public_mode()
               else ",".join(_ORCHESTRATORS))
    keys = [k.strip() for k in raw.split(",") if k.strip() in _ORCHESTRATORS]
    return keys or ["internal"]


def _sidebar_orchestrator() -> str:
    """Choose who drives the Plan→Act loop.

    *Internal* is IGVFagent's own loop in ``_agent.py`` using native function
    calling. *External* shells out to a coding CLI — but note what that does
    and doesn't buy you: ``_llm._chat_claude_cli`` runs ``claude --print`` as a
    **text-generation backend**, parsing tool calls back out of the reply. The
    orchestration is still IGVFagent's, and the CLI's own file/shell tools are
    not used. It is a different transport, not a more capable harness.

    An external option whose binary is absent is shown disabled rather than
    hidden, so the choice is explicable instead of mysteriously missing.
    """
    import shutil

    keys = _allowed_orchestrators()
    if len(keys) < 2:
        # Nothing to choose between — don't render a one-option radio.
        st.session_state["_orchestrator"] = keys[0]
        return keys[0]
    avail = {k: (_ORCHESTRATORS[k][1] is None
                 or shutil.which(_ORCHESTRATORS[k][1]) is not None)
             for k in keys}
    labels = [_ORCHESTRATORS[k][0] + ("" if avail[k] else "  — not installed")
              for k in keys]
    prev = st.session_state.get("_orchestrator", "internal")
    idx = keys.index(prev) if prev in keys else 0

    chosen_label = st.radio(
        "Orchestrator", labels, index=idx, key="_orchestrator_radio",
        help="Internal runs IGVFagent's own agent loop with native function "
             "calling. External shells out to a coding CLI as the text "
             "backend — same tools, different transport.")
    chosen = keys[labels.index(chosen_label)]

    if not avail[chosen]:
        st.warning(
            f"`{_ORCHESTRATORS[chosen][1]}` is not on PATH in this "
            "deployment, so this orchestrator cannot run. Falling back to "
            "the internal loop.")
        chosen = "internal"
    st.session_state["_orchestrator"] = chosen
    return chosen


def _sidebar() -> dict:
    public = _public_mode()
    with st.sidebar:
        st.markdown(f"## 🧬 IGVFagent\n_v{__version__}_")
        # Which build is actually serving this page. Without it, a
        # stale container is indistinguishable from a fixed one.
        st.caption(deployed_build_id())
        st.caption(
            "Natural-language interface to the IGVF / ENCODE single-cell, "
            "variant, regulatory-element, and literature stack."
        )
        st.divider()

        st.subheader("Model")
        if public:
            # Backend is fixed by the operator; the model is chosen from a
            # curated allowlist.
            backend = os.environ.get("IGVF_LLM_BACKEND", "anthropic")
            orch = _sidebar_orchestrator()
            if orch != "internal":
                backend = orch
            choices = _public_model_choices()
            default_model = os.environ.get("IGVF_LLM_MODEL", "").strip()
            try:
                idx = choices.index(default_model)
            except ValueError:
                idx = 0
            model = st.selectbox(
                "Model", choices, index=idx,
                format_func=lambda m: _PUBLIC_MODEL_CATALOG.get(m, (m, ""))[0],
                key="_public_model_select",
                help="Models the IGVF team has enabled on this shared "
                     "deployment. All run on the team's API key.",
            )
            blurb = _PUBLIC_MODEL_CATALOG.get(model, ("", ""))[1]
            if blurb:
                st.caption(blurb)
            st.caption(
                "Running on the IGVF team's API key — no key needed. Install "
                "IGVFagent locally to use another backend, including free "
                "local models via Ollama."
            )
            # Downstream (the active-model banner in main()) branches on this;
            # the public path never runs the radio that would otherwise set it.
            kind = "anthropic"
            # Keep the Load button here too. There are no weights to preload
            # for a cloud API, but it is also the *apply* affordance: it pings
            # the backend to confirm the model id and credential, and it sets
            # the `_loaded_status` the active-model banner reads. Without it a
            # visitor picks a model and gets no confirmation that it took.
            _sidebar_load_button(backend, model)
        else:
            kind = _sidebar_backend_kind()
            if kind == "local":
                backend, model = _sidebar_local_model_picker()
            elif kind == "anthropic":
                backend, model = _sidebar_anthropic_model_picker()
            elif kind == "openai":
                backend, model = _sidebar_openai_model_picker()
            elif kind == "claude_cli":
                backend, model = _sidebar_claude_cli_picker()
            elif kind == "codex_cli":
                backend, model = _sidebar_codex_cli_picker()
            else:
                backend, model = _sidebar_advanced_picker()

            _sidebar_load_button(backend, model)

        # Long-running detached work, surfaced where the user can see it.
        _sidebar_jobs()
        _sidebar_history()

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
        if public:
            # Each iteration is a full LLM call, so the iteration cap is the
            # single biggest lever on what one visitor can spend. Ceilings are
            # operator-set; the visitor may lower them but never raise them.
            #
            # The ceiling has to clear the *whole task*, not a typical one: a
            # bulk job ("download all 76 scE2G sets and integrate them") needs
            # roughly one iteration per object plus discovery and synthesis. A
            # ceiling below that doesn't slow the run down, it ends it early
            # with "Agent run ended before completion".
            iter_cap = _env_int("IGVF_PUBLIC_MAX_ITER", 100)
            token_cap = _env_int("IGVF_PUBLIC_MAX_TOKENS", 16384)
            max_iter = st.slider("Max iterations", 1, iter_cap,
                                   min(_env_int("IGVF_PUBLIC_DEFAULT_ITER", 25),
                                       iter_cap))
            max_tokens = st.slider("Max tokens / turn", 256, token_cap,
                                     min(8192, token_cap), 256)
            temperature = 0.0
            st.caption(
                f"Shared deployment: up to {iter_cap} iterations and "
                f"{token_cap} tokens/turn. Raise **Max iterations** for bulk "
                "jobs — each one is a single LLM call, and a long download or "
                "integration run needs roughly one per object."
            )
        else:
            max_iter = st.slider("Max iterations", 1, 300, 12)
            max_tokens = st.slider("Max tokens / turn", 256, 32000, 4096, 256)
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
        _sidebar_document_upload()

        st.divider()
        _sidebar_user_extensions()

        st.divider()
        st.subheader("⬇ Export session")
        if st.session_state.get("messages"):
            st.caption("Download this session's transcript and artefacts. "
                       "File paths shown in answers are on the server; these "
                       "buttons give you the files.")
            _export_controls(key="sidebar")
        else:
            st.caption("Ask something first — exports appear here.")

        st.divider()
        if st.button("🗑 Clear conversation", **fit(st.button)):
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

# Anything ending in one of these is treated as a viewable artefact when
# referenced from a markdown report body.
_VIEWABLE_EXTS = (
    ".png", ".jpg", ".jpeg", ".svg", ".gif",
    ".csv", ".tsv",
    ".md",
    ".json", ".jsonl",
    ".pdf",
    ".html", ".htm",
    ".txt", ".log",
)


# Pull file paths out of arbitrary text (markdown report bodies, etc.).
# Recognizes backtick-quoted paths and bare absolute paths ending in one
# of the known extensions.
_PATH_IN_TICKS = re.compile(r"`([^`\n]+?\.(?:" +
                            "|".join(ext.lstrip(".") for ext in _VIEWABLE_EXTS) +
                            r"))`")
_BARE_ABS_PATH = re.compile(r"(?<![\w/`])(/[^\s`\n]+?\.(?:" +
                            "|".join(ext.lstrip(".") for ext in _VIEWABLE_EXTS) +
                            r"))(?![\w])")
# Relative paths announced inline (e.g. "Browser SVG: Docs/ENCODE/Plots/x.svg").
# Limited to top-level dirs we actually use so we don't accidentally catch
# inline code snippets or URL paths.
_BARE_REL_PATH = re.compile(
    r"(?<![\w/`])((?:Docs|Data|Scripts|tmp)/[^\s`\n]+?\.(?:" +
    "|".join(ext.lstrip(".") for ext in _VIEWABLE_EXTS) +
    r"))(?![\w])"
)

# Project root used to resolve relative artefact paths.
_PROJECT_ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()


# Output DIRECTORIES as an answer names them -- "Visualizations (in
# Docs/SingleCell/<run>/Plots/)" -- carry no file extension, so the
# extension-keyed patterns below never match them. An answer that then lists
# bare filenames underneath gives the renderer nothing at all: the filenames
# have no directory and the directory has no extension. That is how six
# existing plots were described in prose and none of them displayed.
_OUTPUT_DIR_IN_TEXT = re.compile(r"(?:/workspace/)?((?:Docs|Data)/[A-Za-z0-9_.\-]+(?:/[A-Za-z0-9_.\-]+)*)")


def _extract_paths_from_text(text: str) -> "list[str]":
    found: "list[str]" = []
    seen: "set[str]" = set()
    # Directories are collected separately, because the final filter below
    # is filter_artifacts(), whose default demands a regular FILE -- a
    # directory passed through it is dropped as "not a regular file", which
    # is what silently discarded the run folder on the first attempt. They
    # are vetted with require_file=False, the guard's own provision for
    # callers handling a run directory, then appended after that filter.
    dirs: "list[str]" = []
    for m in _OUTPUT_DIR_IN_TEXT.finditer(text or ""):
        rel = m.group(1).rstrip("/")
        if rel in seen:
            continue
        cand = Path(rel)
        if not cand.is_absolute():
            cand = _PROJECT_ROOT / rel
        try:
            if not cand.is_dir():
                continue
            if _pathguard.why_blocked(cand, require_file=False) is None:
                seen.add(rel)
                dirs.append(str(cand))
        except OSError:
            continue
    for rx in (_PATH_IN_TICKS, _BARE_ABS_PATH, _BARE_REL_PATH):
        for m in rx.finditer(text or ""):
            p = m.group(1).strip()
            if p in seen:
                continue
            # Resolve relatives against the project root so a string like
            # "Docs/ENCODE/Plots/foo.svg" is found regardless of CWD.
            candidates = [p]
            if not p.startswith(("/", "~")):
                candidates.append(str(_PROJECT_ROOT / p))
            for cand in candidates:
                try:
                    if Path(cand).is_file():
                        seen.add(p)
                        found.append(cand)
                        break
                except OSError:
                    continue
    # Containment: the regexes above match *any* path with a viewable
    # extension, including ones outside the workspace and credential files
    # inside it. Everything the UI renders funnels through here, so this is
    # the one place the guard has to hold. See Scripts/_pathguard.py.
    return _pathguard.filter_artifacts(found) + dirs


def _download_button(path: str, key_hint: str = "") -> None:
    try:
        data = Path(path).read_bytes()
        st.download_button(
            label=f"⬇ Download {Path(path).name}",
            data=data,
            file_name=Path(path).name,
            key=f"dl_{key_hint}_{path}",
            **fit(st.download_button, stretch=False),
        )
    except Exception as e:
        st.caption(f"(download unavailable: {e})")


def _render_pdf(path: str) -> None:
    """Embed a PDF inline via base64 + iframe; provide a download button."""
    try:
        data = Path(path).read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        st.markdown(
            f'<iframe src="data:application/pdf;base64,{b64}" '
            f'width="100%" height="640" style="border:1px solid #ddd;">'
            f'</iframe>',
            unsafe_allow_html=True,
        )
        _download_button(path, key_hint="pdf")
    except Exception as e:
        st.text(f"(could not embed PDF: {e})")
        _download_button(path, key_hint="pdf_fallback")


def _render_svg(path: str) -> None:
    """Render an SVG file inline.

    Streamlit ≥1.50 dropped SVG support from ``st.image`` (which now
    requires PIL-loadable raster formats), so SVG paths previously
    returned a broken-image placeholder. We embed the raw SVG inside an
    iframe via ``st.components.v1.html``, which honours the SVG's own
    ``<style>`` / ``<defs>`` / gradients and keeps the layout
    self-contained. Height is parsed from the SVG header so tall plots
    (e.g. the rE2G browser view) don't get clipped.
    """
    try:
        svg = Path(path).read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        st.text(f"(SVG unavailable: {Path(path).name}) {e}")
        _download_button(path, key_hint="svg_fallback")
        return
    # Parse height from the outer <svg ... height="N"> attribute.
    h = 600
    m = re.search(r'<svg\b[^>]*\bheight="([0-9]+(?:\.[0-9]+)?)"', svg)
    if m:
        try:
            h = int(float(m.group(1))) + 24
        except Exception:
            h = 600
    # Make the SVG responsive within its iframe by stamping a
    # ``max-width:100%`` style; preserve the intrinsic aspect ratio.
    svg = re.sub(
        r'<svg\b([^>]*)>',
        r'<svg\1 style="max-width:100%;height:auto;display:block">',
        svg, count=1,
    )
    import streamlit.components.v1 as _stc  # type: ignore
    _stc.html(
        '<div style="background:#fff;padding:6px">' + svg + '</div>',
        height=h, scrolling=True,
    )
    st.caption(Path(path).name)
    _download_button(path, key_hint="svg")


def _render_csv_like(path: str, sep: "str|None" = None) -> None:
    try:
        import pandas as pd
        df = pd.read_csv(path, nrows=400, sep=sep, engine="python")
        st.dataframe(df, **fit(st.dataframe), height=300)
        if len(df) >= 400:
            st.caption("Preview limited to first 400 rows.")
        _download_button(path, key_hint="csv")
    except Exception as e:
        st.text(f"(could not read {path}: {e})")
        _download_button(path, key_hint="csv_fallback")


def _render_jsonl(path: str, max_rows: int = 50) -> None:
    try:
        import json
        rows: "list[Any]" = []
        with open(path, "r", encoding="utf-8") as fh:
            for i, line in enumerate(fh):
                if i >= max_rows:
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except Exception:
                    rows.append({"_raw": line[:500]})
        # Try a tabular view if rows are flat dicts.
        try:
            import pandas as pd
            df = pd.json_normalize(rows)
            st.dataframe(df, **fit(st.dataframe), height=300)
        except Exception:
            st.json(rows)
        st.caption(f"First {len(rows)} record(s).")
        _download_button(path, key_hint="jsonl")
    except Exception as e:
        st.text(f"(could not read {path}: {e})")
        _download_button(path, key_hint="jsonl_fallback")




# Match standard markdown image syntax: ![alt text](path)
# Captures: group(1) = alt, group(2) = path
_MD_IMG_RE = re.compile(r"!\[([^\]]*)\]\(([^)\s]+)\)")


def _resolve_md_image_path(raw: str, base_dir: Path) -> "str|None":
    """Resolve a path from a markdown ![alt](raw) reference.

    Tries: (a) absolute path verbatim, (b) relative to the markdown
    file's directory (so a report under Docs/MPRA can reference
    `Plots/foo.png` and we find Docs/MPRA/Plots/foo.png), (c) relative
    to the project root (so `Docs/MPRA/Plots/foo.png` works too).
    """
    if raw.startswith(("http://", "https://", "data:")):
        return raw  # leave network/data URIs alone
    for cand in (Path(raw),
                  base_dir / raw,
                  _PROJECT_ROOT / raw):
        try:
            if cand.is_file() and _pathguard.is_safe_artifact(cand):
                return str(cand)
        except OSError:
            continue
    return None


def _render_markdown_with_images(body: str, *, base_dir: Path) -> "set[str]":
    """Render a markdown body, but replace ``![alt](path)`` refs with
    real Streamlit image widgets so they actually show up in the browser.

    Returns the **set of resolved file paths the helper rendered inline**
    so callers can dedup against their own "linked artefacts" list
    without accidentally stripping images we did NOT render (e.g. SVGs
    that the report mentions only as backtick-wrapped bullet items, with
    no ``![](...)`` markdown image syntax).

    Falls through to a plain ``st.markdown`` render if no image refs are
    present (the common case).
    """
    rendered_inline: "set[str]" = set()
    matches = list(_MD_IMG_RE.finditer(body))
    if not matches:
        st.markdown(body)
        return rendered_inline
    cursor = 0
    for m in matches:
        # Emit any text before this image
        if m.start() > cursor:
            text = body[cursor:m.start()]
            if text.strip():
                st.markdown(text)
        alt = m.group(1) or Path(m.group(2)).name
        resolved = _resolve_md_image_path(m.group(2), base_dir)
        if resolved is None:
            st.caption(f"_(missing image: {m.group(2)})_")
        elif resolved.startswith(("http://", "https://", "data:")):
            st.markdown(f"![{alt}]({resolved})")
        elif resolved.lower().endswith(".svg"):
            _render_svg(resolved)
            rendered_inline.add(resolved)
        else:
            try:
                st.image(resolved, caption=alt, **fit(st.image))
                rendered_inline.add(resolved)
            except Exception as e:
                st.text(f"(image unavailable: {alt}) {e}")
                _download_button(resolved, key_hint=f"mdimg_{m.start()}")
        cursor = m.end()
    # Tail
    if cursor < len(body):
        tail = body[cursor:]
        if tail.strip():
            st.markdown(tail)
    return rendered_inline


def _render_one(path: str, *, depth: int = 0) -> None:
    """Render a single artefact path with the right widget."""
    low = path.lower()
    name = Path(path).name
    if low.endswith(".svg"):
        _render_svg(path)
        return
    if low.endswith((".png", ".jpg", ".jpeg", ".gif")):
        try:
            st.image(path, caption=name, **fit(st.image))
        except Exception as e:
            st.text(f"(image unavailable: {name}) {e}")
        return
    if low.endswith(".pdf"):
        _render_pdf(path)
        return
    if low.endswith(".csv"):
        _render_csv_like(path, sep=",")
        return
    if low.endswith(".tsv"):
        _render_csv_like(path, sep="\t")
        return
    if low.endswith(".jsonl"):
        _render_jsonl(path)
        return
    if low.endswith(".json"):
        try:
            txt = Path(path).read_text()
            if len(txt) < 200_000:
                st.json(txt)
            else:
                st.caption(f"JSON >200KB; download to inspect.")
            _download_button(path, key_hint="json")
        except Exception as e:
            st.text(f"(could not read {path}: {e})")
        return
    if low.endswith(".md"):
        try:
            body = Path(path).read_text()
        except Exception as e:
            st.text(f"(could not read {path}: {e})")
            return
        # Render the body with `![alt](path)` image references rewritten
        # to real Streamlit image widgets (Streamlit's markdown cannot
        # HTTP-fetch local file paths, so the browser shows a broken icon
        # for every `![alt](Plots/foo.png)` style ref in a report). The
        # helper returns the set of file paths it actually rendered
        # inline; we use that for dedup so images that the report
        # mentions only as bullet items (backtick-wrapped paths, NO
        # `![](...)` syntax) still appear in the linked-artefacts panel.
        rendered_inline = _render_markdown_with_images(body, base_dir=Path(path).parent)
        # Render every other referenced artefact (CSV/JSONL/PDF/etc. AND
        # any image we did NOT already render inline).
        if depth == 0:
            referenced = _extract_paths_from_text(body)
            referenced = [r for r in referenced if r not in rendered_inline]
            if referenced:
                # Images get their own dedicated rendering so they're not
                # buried inside an expander.
                imgs_to_render = [
                    r for r in referenced
                    if r.lower().endswith((".png", ".jpg", ".jpeg",
                                            ".gif", ".svg"))
                ]
                non_imgs = [r for r in referenced if r not in imgs_to_render]
                if imgs_to_render:
                    st.caption(
                        f"Linked images in this report ({len(imgs_to_render)}):"
                    )
                    for img in imgs_to_render:
                        _render_one(img, depth=depth + 1)
                if non_imgs:
                    st.caption(
                        f"Other linked artefacts in this report ({len(non_imgs)}):"
                    )
                    for ref in non_imgs[:12]:
                        with st.expander(f"📎 {Path(ref).name}", expanded=False):
                            _render_one(ref, depth=depth + 1)
        return
    if low.endswith((".html", ".htm")):
        try:
            html = Path(path).read_text()
            if len(html) < 1_500_000:
                # Heuristic: pyvis / vis.js / cytoscape / d3 network HTML
                # files look like network graphs and want a taller iframe.
                lower_html = html[:5000].lower()
                is_network_html = any(
                    marker in lower_html
                    for marker in ("vis-network", "vis.js", "cytoscape",
                                   'class="vis-network"', "force-directed",
                                   "d3.forcesimulation")
                )
                height = 820 if is_network_html else 480
                st.components.v1.html(html, height=height, scrolling=True)
            else:
                st.caption("HTML >1.5MB; download to view.")
        except Exception as e:
            st.text(f"(could not read {path}: {e})")
        _download_button(path, key_hint="html")
        return
    if low.endswith((".txt", ".log")):
        try:
            txt = Path(path).read_text()
            if len(txt) < 50_000:
                st.code(txt)
            else:
                st.code(txt[:50_000])
                st.caption(f"Truncated to 50KB of {len(txt):,} bytes.")
        except Exception as e:
            st.text(f"(could not read {path}: {e})")
        _download_button(path, key_hint="txt")
        return
    # Anything else — surface a download button when it exists on disk.
    st.code(path)
    if Path(path).is_file():
        _download_button(path, key_hint="misc")


def deployed_build_id() -> str:
    """A fingerprint of the Python actually loaded, plus the build SHA if set.

    An external tester reported the hosted site still returning a fixed bug
    and asked us to show the deployed commit, because the UI showed only
    "v0.2.9" whether or not the fix was in the image. A git SHA alone would
    not have helped: this container serves an installed package and has no
    repository, and modules are sometimes copied in without a rebuild. So
    the primary identifier is a content hash of the top-level modules --
    the same construction Deploy/redeploy.sh compares against the checkout,
    so the two can be matched by eye.
    """
    try:
        import hashlib
        pkg = Path(_pathguard.__file__).resolve().parent
        # Reproduce Deploy/redeploy.sh's hash EXACTLY -- per-file sha256 hex,
        # the hex strings sorted as text, joined by newlines with a trailing
        # one, then hashed. A near-miss construction would print a value that
        # looks comparable to redeploy.sh's output and never matches it.
        digests = sorted(hashlib.sha256(f.read_bytes()).hexdigest()
                          for f in pkg.glob("*.py"))
        blob = "".join(d + "\n" for d in digests).encode()
        code = hashlib.sha256(blob).hexdigest()[:12]
    except Exception:
        code = "unknown"
    sha = (os.environ.get("IGVF_GIT_SHA") or "").strip()[:7]
    return f"code {code}" + (f" · build {sha}" if sha else "")


# Run directories are stamped <YYYYmmdd>_<HHMMSS>_<label>. Recognising them
# is what stops a parent directory from being expanded into every past run.
_RUN_DIR_RE = re.compile(r"^\d{8}[_-]\d{6}[_-]")


def _looks_like_run_dir(name: str) -> bool:
    return bool(_RUN_DIR_RE.match(name))


def _norm_path(p: str) -> "Optional[Path]":
    """Absolute, symlink-resolved path, or None if it cannot be resolved.

    Paths reach the UI in several spellings of the same file --
    "Docs/x/y.png", "/workspace/Docs/x/y.png", and the absolute host path --
    so deduplicating the raw strings left the same figure listed three times.
    """
    try:
        q = Path(p)
        if not q.is_absolute():
            q = _PROJECT_ROOT / p
        return q.resolve()
    except (OSError, ValueError, RuntimeError):
        return None


def _run_scope(reported: "list[str]") -> "set[Path]":
    """The directories THIS run wrote into, from what its tools reported.

    A run owns the files its own tools announced, and whatever sits inside
    the directories those files are in. It does not own the parent of those
    directories: Docs/KGTraversal holds 15 run directories from previous
    questions, so expanding the parent listed APOE, TP53, BRCA1 and rs7412
    reports as artefacts of a query about six kidney genes.
    """
    scope: "set[Path]" = set()
    for r in reported:
        q = _norm_path(r)
        if q is None:
            continue
        scope.add(q if q.is_dir() else q.parent)
    return scope


def _in_run_scope(path: "Optional[Path]", scope: "set[Path]") -> bool:
    if path is None or not scope:
        return False
    for d in scope:
        try:
            path.relative_to(d)
            return True
        except ValueError:
            continue
    return False


def _collect_run_artefacts(reported: "list[str]", answer: str) -> "list[str]":
    """Artefacts of THIS run: what its tools produced, plus what the answer
    cites from inside that run's own output directories.

    Paths named in the answer are still honoured -- an answer that cites a
    figure should render it -- but only when they fall inside a directory
    this run wrote to. Anything else is another run's file that the model
    happened to mention, and attributing it to this run is wrong twice over:
    it misrepresents provenance, and on a shared deployment it advertises
    files this questioner never asked for.
    """
    scope = _run_scope(reported)
    out: "list[str]" = []
    seen: "set[Path]" = set()

    def add(raw: str) -> None:
        q = _norm_path(raw)
        if q is None or q in seen:
            return
        seen.add(q)
        try:
            out.append(str(q.relative_to(_PROJECT_ROOT)))
        except ValueError:
            out.append(str(q))

    for r in reported:
        add(r)
    for m in _extract_paths_from_text(answer or ""):
        if _in_run_scope(_norm_path(m), scope):
            add(m)
    return out


def _expand_artefact_dirs(paths: "list[str]") -> "list[str]":
    """Replace directory artefacts with the renderable files inside them.

    A run announces the DIRECTORY its outputs landed in ("Output:
    Docs/SingleCell/<run>"), not each file. On the hosted site a directory is
    useless -- there is no filesystem to browse -- so it renders as a folder
    icon and the six plots inside stay invisible.

    One level deep, because the announced path is the run root while the
    figures sit in <run>/Plots. Every candidate still passes _pathguard, so
    expanding a directory cannot surface credential material that a direct
    path would have been refused.
    """
    renderable = (".png", ".jpg", ".jpeg", ".gif", ".svg",
                  ".md", ".csv", ".tsv", ".json")
    out: "list[str]" = []
    seen: "set[str]" = set()

    def _keep(q: str) -> None:
        if q not in seen:
            seen.add(q)
            out.append(q)

    for p in paths:
        try:
            cand = Path(p)
            if not cand.is_absolute():
                cand = _PROJECT_ROOT / p
            if cand.is_dir():
                found: "list[str]" = []
                # One sublevel deep, and ONLY into a subdirectory that looks
                # like part of this run's output (Plots/, Manifests/, ...) --
                # never into a sibling that is itself a run directory. A run
                # directory is named <timestamp>_<label>, so descending into
                # those turned "the folder my results are in" into "every
                # result anyone has ever produced here".
                for c in sorted(cand.iterdir()):
                    if len(found) >= 40:
                        break
                    if (c.is_file() and c.suffix.lower() in renderable
                            and _pathguard.is_safe_artifact(c)):
                        found.append(str(c))
                    elif c.is_dir() and not _looks_like_run_dir(c.name):
                        for g in sorted(c.iterdir()):
                            if len(found) >= 40:
                                break
                            if (g.is_file() and g.suffix.lower() in renderable
                                    and _pathguard.is_safe_artifact(g)):
                                found.append(str(g))
                if found:
                    for f in found:
                        _keep(f)
                    continue
        except OSError:
            pass
        _keep(p)
    return out


def _render_artefacts(paths: "list[str]") -> None:
    if not paths:
        return

    # Turn announced directories into the files inside them, so a run's
    # figures render instead of a folder icon the user cannot open.
    paths = _expand_artefact_dirs(paths)

    # Dedupe AFTER expansion, on the resolved path, keeping first-seen order.
    #
    # Both halves of that matter. Deduping BEFORE expansion cannot see that a
    # directory and a file inside it are the same artefact: a run announcing
    # "Output: <run dir>" and "Report: <run dir>/report.md" -- which is the
    # normal shape -- expands the directory to report.md and then lists the
    # file again. And comparing raw STRINGS misses it even then, because
    # _collect_run_artefacts stores relative paths while _expand_artefact_dirs
    # emits absolute ones, so "Docs/x/report.md" and
    # "/workspace/Docs/x/report.md" are two entries for one file. That is the
    # duplication reported in the Artefacts panel, still present after the
    # earlier normalisation fix because that fix ran on the wrong side of the
    # expansion.
    seen: "set[Path]" = set()
    deduped: "list[str]" = []
    for p in paths:
        key = _norm_path(p)
        if key is None:
            key = Path(p)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(p)
    paths = deduped

    images, svgs, markdowns, tabular, jsonl_files, jsons, pdfs, others = (
        [], [], [], [], [], [], [], []
    )
    for p in paths:
        low = p.lower()
        if low.endswith(".svg"):
            svgs.append(p)
        elif low.endswith((".png", ".jpg", ".jpeg", ".gif")):
            images.append(p)
        elif low.endswith(".md"):
            markdowns.append(p)
        elif low.endswith((".csv", ".tsv")):
            tabular.append(p)
        elif low.endswith(".jsonl"):
            jsonl_files.append(p)
        elif low.endswith(".json"):
            jsons.append(p)
        elif low.endswith(".pdf"):
            pdfs.append(p)
        else:
            others.append(p)

    # Images: render every PNG/JPG/JPEG/GIF inline.
    # Previously we capped at 6 per category, which silently hid plot-heavy
    # runs (encode_pipeline / multiome / mpra literature-demo all emit
    # 6-12+ figures). Show them all in a 2-column gallery.
    if images:
        cols = st.columns(2)
        for idx, img in enumerate(images):
            try:
                # Column methods share st.image's signature, so probe the
                # module-level function for the supported width spelling.
                cols[idx % 2].image(
                    img, caption=Path(img).name, **fit(st.image),
                )
            except Exception as e:
                cols[idx % 2].text(
                    f"(image unavailable: {Path(img).name}) {e}"
                )
                _download_button(img, key_hint=f"img_fb_{idx}")

    # SVGs: render every one via the dedicated SVG renderer.
    for svg_path in svgs:
        _render_svg(svg_path)

    # Markdown reports — render body inline AND chase any file paths the
    # report references so the user sees the underlying CSV / JSONL / PDF
    # / PNG without having to copy a path into a terminal.
    for md in markdowns:
        with st.expander(f"📄 {Path(md).name}", expanded=False):
            _render_one(md, depth=0)

    for path in pdfs:
        with st.expander(f"📕 {Path(path).name}", expanded=False):
            _render_one(path)

    for path in tabular:
        with st.expander(f"📊 {Path(path).name}", expanded=False):
            _render_one(path)

    for path in jsonl_files:
        with st.expander(f"🧾 {Path(path).name}", expanded=False):
            _render_one(path)

    for path in jsons:
        with st.expander(f"📦 {Path(path).name}", expanded=False):
            _render_one(path)

    for o in others:
        with st.expander(f"📁 {Path(o).name or o}", expanded=False):
            _render_one(o)


_JOBS_REFRESH_S = int(os.environ.get("IGVF_UI_JOBS_REFRESH", "5"))


def _job_figures(rec: dict) -> "list[str]":
    """PNGs a job has written so far, newest last.

    Plots do not land in the job's own directory -- the single-cell step
    announces its own ``Output: Docs/SingleCell/<run>`` and writes into
    ``<run>/Plots`` -- so the log is the only link between a job and its
    figures. They are read live, which is what lets them appear one by one
    while the job is still running.
    """
    out: "list[str]" = []
    log = Path(rec.get("log", ""))
    roots: "list[Path]" = []
    try:
        text = log.read_text(errors="replace") if log.exists() else ""
    except OSError:
        text = ""
    for line in text.splitlines():
        for key in ("Output:", "Report:", "Run dir:"):
            if line.startswith(key):
                cand = Path(line.split(":", 1)[1].strip())
                if not cand.is_absolute():
                    cand = _PROJECT_ROOT / cand
                roots.append(cand)
    work = rec.get("work")
    if work:
        roots.append(Path(work))
    seen: "set[str]" = set()
    for root in roots:
        try:
            if not root.is_dir():
                continue
            for sub in (root, root / "Plots"):
                if not sub.is_dir():
                    continue
                for f in sorted(sub.iterdir(), key=lambda q: q.name):
                    if (f.is_file() and f.suffix.lower() == ".png"
                            and str(f) not in seen
                            and _pathguard.is_safe_artifact(f)):
                        seen.add(str(f))
                        out.append(str(f))
        except OSError:
            continue
    return out

# A detached job's progress is invisible without this. Streamlit renders the
# script once per interaction and never polls, so the panel froze at whatever
# the state was when the page loaded: a user watching a 45-minute alignment
# saw a static "running" line and no moving bar, and reasonably concluded
# nothing was happening. st.fragment(run_every=...) re-runs just this panel
# on a timer, leaving the rest of the page (and any in-flight agent run)
# untouched.
if hasattr(st, "fragment"):
    @st.fragment(run_every=f"{_JOBS_REFRESH_S}s")
    def _jobs_fragment() -> None:
        _jobs_body()
else:                                    # Streamlit < 1.37: no fragments
    def _jobs_fragment() -> None:
        _jobs_body()


# ─── Shared results history ────────────────────────────────────────────────
#
# History is deliberately SHARED: one deployment, one shared password, no
# per-user identity, so pretending otherwise would be a fiction. What matters
# is that it is opt-in -- nothing from it renders until someone opens it and
# picks a run. Pushing every past result at every visitor is the bug this
# replaces.
#
# The index comes from Docs/Agent/<ts>_<query-slug>_<hash>/report.md, which
# already records the original question ("**Query:** ...") alongside the run.
# That makes the conversation text the search key, which is how a person
# actually remembers a past analysis -- not by timestamp.

_ACCESSION_IN_TEXT = re.compile(r"\b(?:IGVF|ENC)[A-Z]{2}[0-9A-Z]{6,}\b")


def _history_dir() -> Path:
    return _PROJECT_ROOT / "Docs" / "Agent"


@st.cache_data(ttl=60, show_spinner=False)
def _history_index(_stamp: float) -> "list[dict]":
    """Every past agent run: when, what was asked, which accessions.

    Cached for a minute and keyed on a coarse timestamp, because 100+ runs
    means 100+ small file reads and the sidebar re-renders constantly.
    """
    root = _history_dir()
    out: "list[dict]" = []
    try:
        dirs = sorted((d for d in root.iterdir() if d.is_dir()),
                      key=lambda d: d.name, reverse=True)
    except OSError:
        return out
    for d in dirs[:300]:
        report = d / "report.md"
        query, when = "", d.name[:15]
        try:
            head = report.read_text(errors="replace")[:4000] if report.exists() else ""
        except OSError:
            head = ""
        for line in head.splitlines():
            if line.startswith("**Query:**"):
                query = line.split("**Query:**", 1)[1].strip()
                break
        if not query:
            # Fall back to the slug in the directory name.
            parts = d.name.split("_")
            query = " ".join(parts[2:-1]).replace("-", " ") or d.name
        accs = sorted(set(_ACCESSION_IN_TEXT.findall(query + " " + head)))
        try:
            ts = time.strftime("%Y-%m-%d %H:%M",
                               time.strptime(d.name[:15], "%Y%m%d_%H%M%S"))
        except ValueError:
            ts = when
        out.append({"dir": str(d), "when": ts, "query": query,
                    "accessions": accs,
                    "haystack": (query + " " + " ".join(accs)).lower()})
    return out


# Output DIRECTORIES as they appear in a run report -- no file extension, so
# _extract_paths_from_text (which keys on known extensions) skips them, and
# they are exactly where the figures live.
_DOCS_DIR_IN_TEXT = re.compile(r"(?:/workspace/)?(Docs/[A-Za-z0-9_.\-/]+)")


def _history_figures(entry: dict) -> "list[str]":
    """Figures belonging to one past run, resolved on demand.

    Resolved lazily rather than indexed: 100+ runs would mean walking every
    output tree on disk to build a sidebar list, and only the run someone
    actually opens needs its plots located.
    """
    report = Path(entry["dir"]) / "report.md"
    try:
        text = report.read_text(errors="replace")
    except OSError:
        return []
    cands: "list[str]" = list(_extract_paths_from_text(text))
    for m in _DOCS_DIR_IN_TEXT.finditer(text):
        rel = m.group(1).rstrip("/")
        cand = _PROJECT_ROOT / rel
        try:
            if cand.is_dir():
                cands.append(str(cand))
        except OSError:
            continue
    expanded = _expand_artefact_dirs(cands)
    return [q for q in expanded if q.lower().endswith(
        (".png", ".jpg", ".jpeg", ".gif", ".svg"))]


def _sidebar_history() -> None:
    """Searchable browser over past runs. Collapsed and inert until opened."""
    with st.expander("📚 Past results", expanded=False):
        entries = _history_index(round(time.time() / 60))
        if not entries:
            st.caption("No past runs recorded yet.")
            return
        needle = st.text_input(
            "Search", key="history_query", placeholder="accession or keywords",
            label_visibility="collapsed",
        ).strip().lower()
        shown = [e for e in entries
                 if not needle or all(t in e["haystack"] for t in needle.split())]
        st.caption(f"{len(shown)} of {len(entries)} runs")
        for e in shown[:25]:
            label = e["query"][:58] + ("…" if len(e["query"]) > 58 else "")
            tag = f"  ·  {', '.join(e['accessions'][:2])}" if e["accessions"] else ""
            if st.button(f"{e['when']}{tag}\n{label}",
                         key=f"hist_{e['dir']}", width="stretch"):
                st.session_state["history_open"] = e["dir"]
        if len(shown) > 25:
            st.caption(f"…{len(shown) - 25} more — narrow the search.")


def history_panel() -> None:
    """Render the past run the user selected, if any."""
    opened = st.session_state.get("history_open")
    if not opened:
        return
    # Match on the resolved path, not the raw string: the selection can
    # arrive relative (a restored session, a link) while the index stores
    # absolute paths, and an exact string compare silently finds nothing.
    def _same(a: str, b: str) -> bool:
        try:
            pa, pb = Path(a), Path(b)
            if not pa.is_absolute():
                pa = _PROJECT_ROOT / pa
            if not pb.is_absolute():
                pb = _PROJECT_ROOT / pb
            return pa.resolve() == pb.resolve()
        except OSError:
            return a == b

    entry = next((e for e in _history_index(round(time.time() / 60))
                  if _same(e["dir"], opened)), None)
    if entry is None:
        st.session_state.pop("history_open", None)
        return
    with st.container(border=True):
        head, close = st.columns([8, 1])
        head.markdown(f"**📚 {entry['when']}** — {entry['query']}")
        if close.button("✕", key="hist_close"):
            st.session_state.pop("history_open", None)
            st.rerun()
        if entry["accessions"]:
            st.caption("Accessions: " + ", ".join(entry["accessions"]))
        figs = _history_figures(entry)
        if figs:
            st.caption(f"{len(figs)} figure(s)")
            cols = st.columns(2)
            for i, f in enumerate(figs):
                try:
                    cols[i % 2].image(f, caption=Path(f).name, **fit(st.image))
                except Exception as e:                       # noqa: BLE001
                    cols[i % 2].caption(f"(cannot show {Path(f).name}: {e})")
        else:
            st.caption("No figures were recorded for this run.")
        rp = Path(entry["dir"]) / "report.md"
        if rp.exists():
            with st.expander("Full run report", expanded=False):
                _render_one(str(rp))


def _session_started_at() -> float:
    """Epoch seconds when this browser session first rendered."""
    if "_session_started_at" not in st.session_state:
        st.session_state["_session_started_at"] = time.time()
    return float(st.session_state["_session_started_at"])


def _own_jobs(jobs: "list[dict]") -> "list[dict]":
    """Only jobs this session started.

    The job registry is shared -- one directory on one VM holding every job
    any visitor has ever launched. Rendering all of them meant that merely
    opening the site produced a wall of somebody else's finished figures
    before a single question had been asked: noise for the reader, and on a
    shared deployment, other people's results.

    A job counts as this session's when it started after the session did.
    Crude, but it needs no plumbing through the agent loop and it fails in
    the safe direction -- an unrecognised job is hidden, never shown.
    """
    started = _session_started_at()
    mine = []
    for j in jobs:
        try:
            t = time.mktime(time.strptime(str(j.get("started")),
                                           "%Y-%m-%d %H:%M:%S"))
        except (ValueError, TypeError):
            continue
        if t >= started - 5:            # slack for clock granularity
            mine.append(j)
    return mine


def _live_jobs_body() -> None:
    """Main-area panel: running jobs and the figures they have produced.

    The sidebar list answers "is it running"; this answers "what has it made
    so far". A long analysis writes its plots one at a time over many
    minutes, and on the hosted site there is no filesystem to watch, so
    without this the figures are invisible until someone thinks to ask again.
    """
    try:
        jobs = _own_jobs(_running_jobs())
    except Exception:                                        # noqa: BLE001
        return
    # "stopped" is included deliberately: a job that produced figures and
    # then died still has results worth showing, and hiding them is how a
    # partial run looks like no run at all.
    recent = jobs[:3]
    if not recent:
        return
    for rec in recent:
        figs = _job_figures(rec)
        running = rec.get("state") == "running"
        if not figs and not running:
            continue
        icon = "🔄" if running else "✅"
        st.markdown(f"{icon} **{rec.get('job_id')}** — `{rec.get('accession')}`"
                    f"  ·  {rec.get('state')}")
        prog = rec.get("progress") or {}
        pct = prog.get("percent")
        if isinstance(pct, (int, float)):
            st.progress(min(1.0, max(0.0, pct / 100.0)),
                        text=f"{prog.get('phase','')} · {prog.get('detail','')}")
        elif prog.get("phase"):
            st.caption(f"⏳ {prog.get('phase')} · {prog.get('detail','')}")
        if figs:
            st.caption(f"{len(figs)} figure(s) so far"
                       + (" — more appear as the run continues" if running else ""))
            cols = st.columns(2)
            for i, f in enumerate(figs):
                try:
                    cols[i % 2].image(f, caption=Path(f).name, **fit(st.image))
                except Exception as e:                       # noqa: BLE001
                    cols[i % 2].caption(f"(cannot show {Path(f).name}: {e})")
        elif running:
            st.caption("No figures written yet.")
        if running:
            st.caption(f"auto-refreshing every {_JOBS_REFRESH_S}s · "
                       f"{time.strftime('%H:%M:%S')}")


if hasattr(st, "fragment"):
    @st.fragment(run_every=f"{_JOBS_REFRESH_S}s")
    def _live_jobs_fragment() -> None:
        _live_jobs_body()
else:
    def _live_jobs_fragment() -> None:
        _live_jobs_body()


def live_jobs_panel() -> None:
    """Live view of work THIS session started. Silent otherwise."""
    try:
        jobs = _own_jobs(_running_jobs())
    except Exception:                                        # noqa: BLE001
        return
    if not jobs:
        return
    if any(j.get("state") == "running" for j in jobs):
        _live_jobs_fragment()
    else:
        _live_jobs_body()


def _sidebar_jobs() -> None:
    """Job panel. Auto-refreshing while anything is running."""
    try:
        jobs = _running_jobs()
    except Exception:                                        # noqa: BLE001
        return
    if not jobs:
        return
    if any(j["state"] == "running" for j in jobs):
        _jobs_fragment()
    else:
        _jobs_body()


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
        # Derive the orchestrator from the EVENT PAYLOAD, not session state.
        # The payload's `backend` is what the run actually used, so the label
        # cannot disagree with the run it describes — whereas session state is
        # a UI value that may be unset, stale, or absent from the context this
        # renderer executes in.
        which = ("External" if str(p.get("backend", "")).endswith("_cli")
                 else "Internal")
        return (f"▶ **{which} orchestrator engaged** — Plan → Action → "
                f"Results → Evaluation loop · backend `{p.get('backend')}`, "
                f"model `{p.get('model')}`, {p.get('n_tools')} tools, "
                f"max {p.get('max_iterations')} iters")
    if k == "llm_call_start":
        return (f"🧠 **Plan step** (orchestrator) — iter {p['iteration']}/"
                f"{p['max_iterations']}, {p['n_messages']} messages sent "
                f"to LLM brain")
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
        return (f"✓ **Orchestrator finished** — {p['iterations']} iters, "
                f"{p['tool_calls_made']} tool calls, stop "
                f"`{p['stop_reason']}`")
    return f"_{k}_  `{p}`"


# --------------------------- Main render -----------------------------------


def _running_jobs() -> "list[dict]":
    """Detached raw-pipeline jobs, newest first.

    A long alignment runs outside the request that started it, so without
    this the UI shows nothing at all while 45 GB is being quantified and the
    run looks like it never happened.
    """
    try:
        from igvfagent import raw_data_pipeline as rp
    except Exception:                                        # noqa: BLE001
        try:
            import raw_data_pipeline as rp                   # type: ignore
        except Exception:                                    # noqa: BLE001
            return []
    out = []
    try:
        paths = sorted(rp.JOB_DIR.glob("*.json"),
                       key=lambda q: q.stat().st_mtime, reverse=True)
    except OSError:
        return []
    for jp in paths[:8]:
        try:
            rec = json.loads(jp.read_text())
        except (ValueError, OSError):
            continue
        log = Path(rec.get("log", ""))
        alive = rp._pid_alive(rec.get("pid"), marker="raw-pipeline")
        tail = rp._tail(log, 3) if log.exists() else []
        finished = any("Run dir:" in ln for ln in rp._tail(log, 40)) if log.exists() else False
        rec["state"] = "running" if alive else ("finished" if finished else "stopped")
        rec["tail"] = tail
        rec["progress"] = rp.read_progress(Path(rec.get("work", "")))
        out.append(rec)
    return out


def _heavy_slot_line() -> str:
    """One line on analysis-slot usage, or empty when unavailable."""
    try:
        from igvfagent import _joblock
    except Exception:                                        # noqa: BLE001
        try:
            import _joblock                                  # type: ignore
        except Exception:                                    # noqa: BLE001
            return ""
    try:
        return _joblock.describe()
    except Exception:                                        # noqa: BLE001
        return ""


def _jobs_body() -> None:
    """Render the job list. Re-reads state on every call, so a fragment
    wrapper showing it repeatedly reports current progress rather than a
    snapshot from page load."""
    # Scoped to this session for the same reason the main panel is: the
    # registry is shared across every visitor to the deployment, and one
    # person's accessions are not another's business.
    jobs = _own_jobs(_running_jobs())
    if not jobs:
        st.caption("No analysis jobs in this session.")
        return
    active = [j for j in jobs if j["state"] == "running"]
    header = (f"⚙ Analysis jobs ({len(active)} running)" if active
              else "⚙ Analysis jobs")
    with st.expander(header, expanded=bool(active)):
        for j in jobs:
            icon = {"running": "🔄", "finished": "✅"}.get(j["state"], "⏹")
            st.markdown(f"{icon} **{j.get('job_id')}** — `{j.get('accession')}`  \n"
                        f"<span style='opacity:.7'>{j['state']} · started "
                        f"{j.get('started','?')}</span>",
                        unsafe_allow_html=True)
            prog = j.get("progress") or {}
            pct = prog.get("percent")
            if isinstance(pct, (int, float)):
                st.progress(min(1.0, max(0.0, pct / 100.0)),
                            text=f"{prog.get('phase','')} · {prog.get('detail','')}")
            elif prog.get("phase"):
                st.caption(f"⏳ {prog['phase']} · {prog.get('detail','')}")
            for ln in j.get("tail", []):
                st.caption(ln[:110])
        slots = _heavy_slot_line()
        if slots:
            st.caption(f"⚙ {slots}")
        if active:
            st.caption(f"Refreshing every {_JOBS_REFRESH_S}s · "
                       f"{time.strftime('%H:%M:%S')}")
        st.caption("Jobs keep running after the page closes. Ask "
                   "\u201cstatus of <job id>\u201d for detail.")


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

    # The banner must reflect the ORCHESTRATOR actually selected, not the
    # backend family. In public mode `kind` is pinned to "anthropic", so
    # reading it alone reported "Anthropic Claude API · internal orchestrator"
    # even when the external orchestrator was chosen.
    orch = st.session_state.get("_orchestrator", "internal")
    if orch != "internal":
        kind_label = _ORCHESTRATORS.get(orch, (kind_label, None))[0]

    # Counts were hardcoded ("your 34 IGVFagent skills") and had drifted.
    try:
        n_tools = len(_tools.list_tools())
        n_exposed = len(_llm.canonical_tools(
            [{"name": t.name, "description": t.description or ""}
             for t in _tools.list_tools()]) or [])
    except Exception:
        n_tools = n_exposed = 0
    scale = (f"{n_exposed} of {n_tools} tools exposed"
             if n_tools and n_exposed < n_tools else f"{n_tools} tools")

    # One line, no explanation: which orchestrator is active, and the scale
    # of the tool set. The mechanics of internal vs external belong in the
    # docs, not in a status banner the user reads on every query.
    orch_name = "Internal orchestrator" if orch == "internal" else "External orchestrator"
    how = f"**{orch_name}**  ·  {scale}"

    if status and status.get("ok"):
        st.success(
            f"### {kind_label}  ·  `{status['model']}`  "
            f"·  ✅ Loaded  "
            f"_(in {status['secs']:.1f}s)_\n\n" + how
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

    # ------------------------------------------------------------------
    # Five tabs:
    #   Chat              — the existing LLM-driven flow
    #   Knowledge Graph   — interactive view of the local SQLite KGs
    #   Single-cell       — UMAP / t-SNE / cluster / marker viewer over
    #                       any .h5ad produced by sc-analyze
    #   Network           — context-specific subnetwork views
    #   Spatial           — Spatial-ATAC-Hi-C: per-pixel QC, tissue maps,
    #                       compartments, copy-number clones, loops
    #   Benchmarks        — every committed reproducibility figure, the
    #                       suite dashboard, and the latest concordance run
    # Every tab except Chat is independent of the loaded LLM.
    # ------------------------------------------------------------------
    chat_tab, kg_tab, sc_tab, nw_tab, sp_tab, bm_tab = st.tabs(
        ["💬 Chat", "🕸  Knowledge Graph", "🔬 Single-cell",
         "🔗 Network", "🧬 Spatial", "📊 Benchmarks"]
    )

    with kg_tab:
        if _kgviz is None:
            st.warning(
                "Knowledge Graph visualizer not available — needs "
                "`matplotlib` + `networkx` in this venv. Install with:\n\n"
                "```\npip install matplotlib networkx\n```"
            )
        else:
            try:
                _kgviz.render_streamlit_panel(st)
            except Exception as exc:  # pylint: disable=broad-except
                import traceback
                st.error(f"KG visualizer error: {exc}")
                with st.expander("Traceback", expanded=False):
                    st.code(traceback.format_exc())

    with sc_tab:
        if _scviz is None:
            st.warning(
                "Single-cell visualizer not available — needs `scanpy` + "
                "`anndata` + `matplotlib` in this venv. Install with:\n\n"
                "```\n"
                "pip install scanpy 'anndata>=0.10' umap-learn leidenalg "
                "python-igraph matplotlib\n"
                "```"
            )
        else:
            try:
                _scviz.render_streamlit_panel(st)
            except Exception as exc:  # pylint: disable=broad-except
                import traceback
                st.error(f"Single-cell visualizer error: {exc}")
                with st.expander("Traceback", expanded=False):
                    st.code(traceback.format_exc())

    with nw_tab:
        if _nwviz is None:
            st.warning(
                "Network visualizer not available — needs `networkx`, "
                "`matplotlib`, `pyvis`, and `pandas` in this venv. Install with:\n\n"
                "```\npip install networkx matplotlib pyvis pandas\n```"
            )
        else:
            try:
                _nwviz.render_streamlit_panel(st)
            except Exception as exc:  # pylint: disable=broad-except
                import traceback
                st.error(f"Network visualizer error: {exc}")
                with st.expander("Traceback", expanded=False):
                    st.code(traceback.format_exc())

    with sp_tab:
        if _sphic is None:
            st.warning(
                "Spatial-ATAC-Hi-C browser not available — needs `numpy` and "
                "`matplotlib` in this venv. Install with:\n\n"
                "```\npip install 'igvfagent[analysis]'\n```"
            )
        else:
            try:
                _sphic.render_streamlit_panel(st)
            except Exception as exc:  # pylint: disable=broad-except
                import traceback
                st.error(f"Spatial-ATAC-Hi-C browser error: {exc}")
                with st.expander("Traceback", expanded=False):
                    st.code(traceback.format_exc())

    with bm_tab:
        if _bmviz is None:
            st.warning(
                "Benchmark visualizer not available — "
                "`Scripts/benchmark_visualizer.py` could not be imported."
            )
        else:
            try:
                _bmviz.render_streamlit_panel(st)
            except Exception as exc:  # pylint: disable=broad-except
                import traceback
                st.error(f"Benchmark visualizer error: {exc}")
                with st.expander("Traceback", expanded=False):
                    st.code(traceback.format_exc())

    # ------------------------------------------------------------------
    # Chat tab — history replay + suggestions. Rendered inside a proper
    # `with chat_tab:` context so all child widgets land in the tab.
    # `st.chat_input` is intentionally placed OUTSIDE the tabs at page
    # level — Streamlit only docks it to the viewport bottom when it
    # lives at the top of the script body, not inside a layout
    # container like `st.tabs(...)`.
    # ------------------------------------------------------------------
    if "messages" not in st.session_state:
        st.session_state.messages = []

    with chat_tab:
        # Live view of detached work, above the transcript: a job started in
        # an earlier turn keeps producing figures, and this is where they
        # show up without the user having to ask again.
        history_panel()
        live_jobs_panel()

        # Replay prior conversation
        for entry in st.session_state.messages:
            with st.chat_message(entry["role"]):
                _render_markdown_with_images(entry.get("content", ""),
                                              base_dir=_PROJECT_ROOT)
                if entry.get("artefacts"):
                    _render_artefacts(entry["artefacts"])
                if entry.get("meta"):
                    st.caption(entry["meta"])

        # The turn about to run renders HERE, inside the tab and directly
        # after the replayed transcript. Without this it was emitted at page
        # level -- after the whole tabs widget -- so a question and its
        # answer appeared detached from the conversation, and the question
        # looked like it had never been asked.
        st.session_state["_turn_area"] = st.container()

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

    # Page-level chat input — Streamlit auto-pins this to the bottom of
    # the viewport because it is not inside a tab / column / container.
    submitted = st.chat_input(
        "Ask IGVFagent — natural-language query, or attach a paper (📎)",
        accept_file="multiple",
        file_type=["pdf", "docx", "txt", "md", "csv", "tsv",
                   "png", "jpg", "jpeg", "svg"],
    )
    if not submitted:
        return

    # With accept_file set, chat_input returns an object carrying .text and
    # .files rather than a bare string. Handle both so the widget can be
    # reverted without breaking this path.
    if isinstance(submitted, str):
        query, attachments = submitted, []
    else:
        query = (getattr(submitted, "text", "") or "").strip()
        attachments = list(getattr(submitted, "files", []) or [])

    saved_paths: "list[str]" = []
    if attachments:
        dest_dir = _upload_dir()
        for up in attachments:
            # Basename only: an uploaded filename is untrusted input and must
            # not be able to escape the upload directory.
            dest = dest_dir / Path(up.name).name
            dest.write_bytes(up.getbuffer())
            saved_paths.append(str(dest.relative_to(_PROJECT_ROOT)))
        st.session_state["_uploaded_docs"] = [
            str(q) for q in sorted(dest_dir.iterdir()) if q.is_file()]

    if not query and not saved_paths:
        return

    if saved_paths:
        # Name the files in the query itself. The agent selects tools from the
        # text it is given, so an attachment it is never told about may as well
        # not exist — and document_plan needs the path, not the bytes.
        listed = ", ".join(f"`{p}`" for p in saved_paths)
        if query:
            query = (f"{query}\n\n[Attached file(s), saved in the workspace: "
                     f"{listed}]")
        else:
            query = (f"Read the attached file(s) — {listed} — and summarise "
                     f"what they contain. If a file is a manuscript, use "
                     f"document_plan to extract its accessions, assays and a "
                     f"reproduction plan.")

    # Everything below renders the new turn — keep it inside the chat
    # tab so the new exchange lands above the docked input.
    chat_tab.__enter__()

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
    turn = st.session_state.get("_turn_area") or st.container()
    with turn:
        with st.chat_message("user"):
            st.markdown(query)

    if pre_errors:
        with turn.chat_message("assistant"):
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

    with turn.chat_message("assistant"):
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
            # Prior turns, excluding the one just appended for this query.
            # Artefact paths are folded into the assistant text because a
            # follow-up ("reproduce that analysis") needs the paths, and those
            # live in message metadata rather than in the prose.
            prior = []
            for m in st.session_state.messages[:-1]:
                content = m.get("content") or ""
                if m.get("role") == "assistant" and m.get("artefacts"):
                    paths = "\n".join(f"- {a}" for a in m["artefacts"][:12])
                    content = f"{content}\n\nArtefacts produced:\n{paths}"
                prior.append({"role": m.get("role"), "content": content})

            result = _agent.run(
                query,
                backend=cfg["backend"],
                model=cfg["model"],
                max_iterations=cfg["max_iter"],
                max_tokens=cfg["max_tokens"],
                temperature=cfg["temperature"],
                tools_subset=cfg["tool_filter"],
                history=prior,
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
        # the message is impossible to miss; otherwise render markdown
        # with image-aware chunking so any `![alt](path)` refs the LLM
        # included end up as real widgets, not broken-image icons.
        if result.stop_reason == "complete_with_failures":
            n = getattr(result, "tool_calls_failed", 0)
            st.error(f"{n} tool call(s) did not succeed — this answer is "
                      f"incomplete. The failures are listed at the top of it.")
            _render_markdown_with_images(result.final_answer,
                                          base_dir=_PROJECT_ROOT)
        elif result.stop_reason != "complete" and result.final_answer:
            st.error("Agent run ended before completion. See details below.")
            _render_markdown_with_images(result.final_answer,
                                          base_dir=_PROJECT_ROOT)
        elif result.final_answer:
            _render_markdown_with_images(result.final_answer,
                                          base_dir=_PROJECT_ROOT)
        else:
            st.warning("_The agent finished without producing a final answer._")

        # Inline artefact rendering — combine the agent's declared
        # artefacts with any file paths mentioned in the final answer
        # itself, so the user does not have to copy paths into a terminal
        # to view a referenced CSV / JSONL / PDF / PNG.
        # Only this run's own output. Previously every path mentioned in the
        # answer was added, so an answer that cited an earlier report listed
        # it as an artefact of this run; an external evaluation saw ~70
        # entries for a six-gene query, including APOE, TP53 and BRCA1.
        artefacts = _collect_run_artefacts(list(result.artefacts or []),
                                            result.final_answer or "")
        if artefacts:
            with st.expander(f"📁 Artefacts ({len(artefacts)})",
                             expanded=True):
                _render_artefacts(artefacts)

        meta_caption = (
            (f"**{getattr(result, 'tool_calls_failed', 0)} failed**  ·  "
             if getattr(result, "tool_calls_failed", 0) else "")
            + f"build `{deployed_build_id()}`  ·  "
            f"backend `{result.backend}`  ·  model `{result.model}`  ·  "
            f"{result.iterations} iter · {result.tool_calls_made} tool calls "
            f"· stop `{result.stop_reason}`"
        )
        if result.report_path:
            meta_caption += f"  ·  report `{result.report_path}`"
        st.caption(meta_caption)

        # Hand over the bytes, not a path the user cannot reach.
        answer_md = (f"# {query}\n\n{result.final_answer or ''}\n\n"
                     + ("\n".join(f"- `{a}`" for a in artefacts) if artefacts else "")
                     + f"\n\n_{meta_caption}_\n")
        dl = st.columns([1, 1, 2])
        with dl[0]:
            st.download_button(
                "⬇ Answer (.md)", data=answer_md,
                file_name=f"igvfagent_answer_{time.strftime('%Y%m%d_%H%M%S')}.md",
                mime="text/markdown",
                key=f"dl_ans_{len(st.session_state.get('messages', []))}",
                **fit(st.download_button, stretch=False))
        with dl[1]:
            if result.report_path and _pathguard.is_safe_artifact(result.report_path):
                _download_button(result.report_path, key_hint="agent_report")

    st.session_state.messages.append({
        "role": "assistant",
        "content": result.final_answer,
        "artefacts": artefacts,
        "meta": meta_caption,
    })


def _session_markdown() -> str:
    """The whole conversation as one self-contained markdown document.

    Printing a server-side path — `report /workspace/Docs/Agent/…/report.md` —
    tells a hosted user where a file they cannot reach lives. Users were
    copy-pasting answers out of the browser instead. These exports hand over
    the actual bytes.
    """
    lines = [f"# IGVFagent session", "",
             f"Exported: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
             f"Version: {__version__}", ""]
    for m in st.session_state.get("messages", []):
        who = "You" if m.get("role") == "user" else "IGVFagent"
        lines += [f"## {who}", "", (m.get("content") or "").strip(), ""]
        arts = m.get("artefacts") or []
        if arts:
            lines += ["**Artefacts produced**", ""]
            lines += [f"- `{a}`" for a in arts] + [""]
        if m.get("meta"):
            lines += [f"_{m['meta']}_", ""]
    return "\n".join(lines)


def _session_zip() -> "bytes | None":
    """Session markdown plus every artefact the session produced.

    Artefacts are filtered through _pathguard, so an export cannot become a
    route to files outside the workspace or to secret material inside it.
    """
    import io
    import zipfile
    msgs = st.session_state.get("messages", [])
    if not msgs:
        return None
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("session.md", _session_markdown())
        seen = set()
        for m in msgs:
            for a in (m.get("artefacts") or []):
                if a in seen or not _pathguard.is_safe_artifact(a):
                    continue
                seen.add(a)
                pa = Path(a)
                try:
                    z.write(pa, arcname=f"artefacts/{pa.name}")
                except OSError:
                    continue
    return buf.getvalue()


def _export_controls(*, key: str) -> None:
    """Download buttons for the current session."""
    msgs = st.session_state.get("messages", [])
    if not msgs:
        return
    md = _session_markdown()
    stamp = time.strftime("%Y%m%d_%H%M%S")
    cols = st.columns(2)
    with cols[0]:
        st.download_button("⬇ Session (.md)", data=md,
                           file_name=f"igvfagent_session_{stamp}.md",
                           mime="text/markdown", key=f"dl_md_{key}",
                           **fit(st.download_button, stretch=True))
    with cols[1]:
        blob = _session_zip()
        if blob:
            st.download_button("⬇ Session + artefacts (.zip)", data=blob,
                               file_name=f"igvfagent_session_{stamp}.zip",
                               mime="application/zip", key=f"dl_zip_{key}",
                               **fit(st.download_button, stretch=True))


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
