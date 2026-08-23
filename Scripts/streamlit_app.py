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
from pathlib import Path
from typing import Any

import streamlit as st

# Dual-mode import (installed package OR running from a checkout).
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from igvfagent import _agent, _llm, _tools, _userext, __version__
    from igvfagent import load_dotenv as _load_dotenv
except Exception:
    import _agent  # type: ignore
    import _llm    # type: ignore
    import _tools  # type: ignore
    import _userext  # type: ignore
    try:
        from __init__ import __version__, load_dotenv as _load_dotenv  # type: ignore
    except Exception:
        __version__ = "0.1.0"
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
        "claude-fable-5",
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

    uploads = st.file_uploader(
        "Add extension files",
        type=["yaml", "yml", "json", "py"],
        accept_multiple_files=True,
        help="Tool manifests (.yaml/.json) and skill modules (.py). "
             "Installed into ~/.igvfagent/ and picked up immediately — "
             "no restart needed.",
    )
    if uploads and st.button("📦 Install extensions",
                              use_container_width=True):
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
        elif kind == "claude_cli":
            backend, model = _sidebar_claude_cli_picker()
        elif kind == "codex_cli":
            backend, model = _sidebar_codex_cli_picker()
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
        max_iter = st.slider("Max iterations", 1, 60, 12)
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
        _sidebar_user_extensions()

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


def _extract_paths_from_text(text: str) -> "list[str]":
    found: "list[str]" = []
    seen: "set[str]" = set()
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
    return found


def _download_button(path: str, key_hint: str = "") -> None:
    try:
        data = Path(path).read_bytes()
        st.download_button(
            label=f"⬇ Download {Path(path).name}",
            data=data,
            file_name=Path(path).name,
            key=f"dl_{key_hint}_{path}",
            use_container_width=False,
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
        st.dataframe(df, use_container_width=True, height=300)
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
            st.dataframe(df, use_container_width=True, height=300)
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
            if cand.is_file():
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
                st.image(resolved, caption=alt, use_container_width=True)
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
            st.image(path, caption=name, use_container_width=True)
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


def _render_artefacts(paths: "list[str]") -> None:
    if not paths:
        return
    # Dedupe while preserving order.
    seen: "set[str]" = set()
    paths = [p for p in paths if not (p in seen or seen.add(p))]

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
                cols[idx % 2].image(
                    img, caption=Path(img).name, use_container_width=True,
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
        return (f"▶ **Internal orchestrator engaged** — Plan → Action → "
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
            f"Every chat query runs through IGVFagent's **internal "
            f"orchestrator** (Plan → Action → Results → Evaluation). "
            f"This LLM is the brain at the planning step; your 34 "
            f"IGVFagent skills are the tools the orchestrator calls."
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
    # Three tabs:
    #   Chat              — the existing LLM-driven flow
    #   Knowledge Graph   — interactive view of the local SQLite KGs
    #   Single-cell       — UMAP / t-SNE / cluster / marker viewer over
    #                       any .h5ad produced by sc-analyze
    # The KG and Single-cell tabs are independent of the loaded LLM.
    # ------------------------------------------------------------------
    chat_tab, kg_tab, sc_tab, nw_tab = st.tabs(
        ["💬 Chat", "🕸  Knowledge Graph", "🔬 Single-cell",
         "🔗 Network"]
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
        # Replay prior conversation
        for entry in st.session_state.messages:
            with st.chat_message(entry["role"]):
                _render_markdown_with_images(entry.get("content", ""),
                                              base_dir=_PROJECT_ROOT)
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

    # Page-level chat input — Streamlit auto-pins this to the bottom of
    # the viewport because it is not inside a tab / column / container.
    query = st.chat_input(
        "Ask IGVFagent — natural-language query (e.g. 'Walk the KG for APOE')"
    )
    if not query:
        return

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
        # the message is impossible to miss; otherwise render markdown
        # with image-aware chunking so any `![alt](path)` refs the LLM
        # included end up as real widgets, not broken-image icons.
        if result.stop_reason != "complete" and result.final_answer:
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
        artefacts = list(result.artefacts or [])
        extra = _extract_paths_from_text(result.final_answer or "")
        for p in extra:
            if p not in artefacts:
                artefacts.append(p)
        if artefacts:
            with st.expander(f"📁 Artefacts ({len(artefacts)})",
                             expanded=True):
                _render_artefacts(artefacts)

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
        "artefacts": artefacts,
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
