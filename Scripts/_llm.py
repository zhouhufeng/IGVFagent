"""LLM provider router for IGVFagent.

A single ``chat(messages, ...)`` API that talks to:

  * Anthropic Claude   (``anthropic`` SDK)
  * OpenAI             (``openai`` SDK)
  * Codex API          (``openai`` SDK; optional ``IGVF_CODEX_BASE_URL``)
  * Ollama (local)     (``openai`` SDK pointed at ``http://localhost:11434/v1``)
  * vLLM / HF TGI / DeepInfra / Together / Groq
                       (``openai`` SDK with custom ``IGVF_LLM_BASE_URL``)

The router exposes a backend-neutral message + tool-call shape so the
ReAct agent runner (step 3) and the Streamlit UI (step 4) don't need to
care which provider answered.

Backend selection precedence:

  1. ``backend=`` argument to ``chat()``
  2. ``IGVF_LLM_BACKEND`` environment variable
  3. Inferred from model name (``claude-*`` -> Anthropic, ``gpt-*`` /
     ``o1-*`` / ``o3-*`` -> OpenAI, ``qwen*`` / ``llama*`` / ``mistral*``
     -> Ollama)
  4. Default: Ollama with Qwen 3 8B (``qwen3:8b``)

Credentials are read **only** from environment variables — never from
source. No URLs leak into the wheel; OpenAI-compatible base URLs default
to the public endpoints but every step is overridable.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import os
from typing import Any, Callable, Iterable, Optional

logger = logging.getLogger(__name__)


# --------------------------- Public types -----------------------------------


@dataclasses.dataclass
class ToolCall:
    """A model-issued tool invocation, normalized across backends."""

    id: str
    name: str
    arguments: dict


@dataclasses.dataclass
class Message:
    """The result of a single ``chat()`` call.

    ``content`` is the assistant's text. ``tool_calls`` is non-empty when
    the model decided to invoke one or more registered tools instead of
    (or in addition to) replying with text. ``stop_reason`` is normalized
    across backends to one of: ``end_turn``, ``tool_use``,
    ``max_tokens``, ``stop_sequence``, ``error``, ``other``.
    """

    role: str = "assistant"
    content: str = ""
    tool_calls: "list[ToolCall]" = dataclasses.field(default_factory=list)
    stop_reason: str = "end_turn"
    backend: str = ""
    model: str = ""
    raw: Optional[dict] = None
    usage: Optional[dict] = None


# --------------------------- Backend selection ------------------------------


_DEFAULT_MODELS = {
    "anthropic":  "claude-opus-5",
    "openai":     "gpt-4o-mini",
    "codex":      "gpt-5-codex",
    "ollama":     "qwen3:8b",
    "vllm":       "qwen3:8b",
    "tgi":        "qwen3:8b",
    "groq":       "llama-3.1-70b-versatile",
    "together":   "meta-llama/Meta-Llama-3.1-70B-Instruct-Turbo",
    "deepinfra":  "meta-llama/Meta-Llama-3.1-70B-Instruct",
    # claude_cli: empty default lets Claude Code use whatever it's
    # configured for in `claude config` (no --model override).
    "claude_cli": "",
    # codex_cli: empty default lets the user's `codex` config decide.
    "codex_cli":  "",
}


def _infer_backend(model: Optional[str]) -> str:
    if not model:
        return os.environ.get("IGVF_LLM_BACKEND", "ollama")
    m = model.lower()
    if "claude_cli" in m or "claude-cli" in m or m == "cc":
        return "claude_cli"
    if "codex_cli" in m or "codex-cli" in m or m == "cx":
        return "codex_cli"
    if "claude" in m:
        return "anthropic"
    if any(p in m for p in ("gpt-", "o1-", "o3-", "o4-")):
        return "openai"
    if "codex" in m:
        return "codex"
    if any(p in m for p in ("qwen", "llama", "mistral", "phi", "deepseek",
                              "gemma", "yi", "command")):
        return "ollama"
    return os.environ.get("IGVF_LLM_BACKEND", "ollama")


# Backends whose inference runs on hardware the user controls. Everything
# else transmits the prompt to a third party, whatever it is called.
_LOCAL_BACKENDS = frozenset({"ollama", "vllm", "tgi"})


def _resolve_backend(backend: Optional[str], model: Optional[str]) -> str:
    if backend:
        return backend.lower()
    env = os.environ.get("IGVF_LLM_BACKEND")
    if env:
        return env.lower()
    return _infer_backend(model)


# --------------------------- Tool format adapters ---------------------------

def to_anthropic_tools(tools: Iterable[dict]) -> "list[dict]":
    """Convert internal tool dicts to Anthropic ``tools=[{name, description,
    input_schema}]`` shape."""
    out = []
    for t in tools:
        out.append({
            "name":         t["name"],
            "description":  t.get("description", ""),
            "input_schema": t.get("parameters", {"type": "object",
                                                    "properties": {}}),
        })
    return out


def to_openai_tools(tools: Iterable[dict]) -> "list[dict]":
    """Convert internal tool dicts to OpenAI Chat Completions
    ``tools=[{type:'function', function:{...}}]`` shape."""
    out = []
    for t in tools:
        out.append({
            "type": "function",
            "function": {
                "name":        t["name"],
                "description": t.get("description", ""),
                "parameters":  t.get("parameters", {"type": "object",
                                                      "properties": {}}),
            },
        })
    return out


# --------------------------- Message format adapters ------------------------

def _to_anthropic_messages(messages: Iterable[dict]
                            ) -> "tuple[Optional[str], list[dict]]":
    """Anthropic uses a separate ``system=`` argument; everything else is a
    user/assistant turn with content blocks. Tool results are
    ``content=[{type:'tool_result', tool_use_id, content}]`` user
    messages."""
    system = None
    out = []
    for m in messages:
        role = m.get("role", "user")
        if role == "system":
            system = m.get("content", "") if system is None else \
                     system + "\n\n" + m.get("content", "")
            continue
        if role == "tool":
            out.append({
                "role": "user",
                "content": [{
                    "type":         "tool_result",
                    "tool_use_id":  m.get("tool_call_id", ""),
                    "content":      m.get("content", ""),
                }],
            })
            continue
        if role == "assistant" and m.get("tool_calls"):
            blocks = []
            if m.get("content"):
                blocks.append({"type": "text", "text": m["content"]})
            for tc in m["tool_calls"]:
                blocks.append({
                    "type":  "tool_use",
                    "id":    tc.get("id") if isinstance(tc, dict) else tc.id,
                    "name":  tc.get("name") if isinstance(tc, dict) else tc.name,
                    "input": tc.get("arguments") if isinstance(tc, dict)
                              else tc.arguments,
                })
            out.append({"role": "assistant", "content": blocks})
            continue
        out.append({"role": role, "content": m.get("content", "")})
    return system, out


def _to_openai_messages(messages: Iterable[dict]) -> "list[dict]":
    """OpenAI Chat Completions: system / user / assistant / tool messages
    in a flat list. Tool calls live on the assistant message; tool results
    are ``role='tool'`` with ``tool_call_id``."""
    out = []
    for m in messages:
        role = m.get("role", "user")
        if role == "tool":
            out.append({
                "role":          "tool",
                "tool_call_id":  m.get("tool_call_id", ""),
                "content":       m.get("content", ""),
            })
            continue
        if role == "assistant" and m.get("tool_calls"):
            tool_calls = []
            for tc in m["tool_calls"]:
                tc_id = tc.get("id") if isinstance(tc, dict) else tc.id
                tc_name = tc.get("name") if isinstance(tc, dict) else tc.name
                tc_args = tc.get("arguments") if isinstance(tc, dict) \
                          else tc.arguments
                tool_calls.append({
                    "id":   tc_id,
                    "type": "function",
                    "function": {
                        "name":      tc_name,
                        "arguments": json.dumps(tc_args, default=str),
                    },
                })
            out.append({
                "role":       "assistant",
                "content":    m.get("content", "") or "",
                "tool_calls": tool_calls,
            })
            continue
        out.append({"role": role, "content": m.get("content", "")})
    return out


# --------------------------- Backend implementations ------------------------

# Anthropic models that no longer accept temperature/top_p/top_k. Sending any
# sampling param to these returns a 400 ("temperature is deprecated for this
# model."). Covers Opus 4.7/4.8, the Claude 5 family (Opus 5 / Sonnet 5 /
# Fable 5 / Mythos 5), and any later release. Sonnet 5 rejects only
# *non-default* values, but we send temperature 0.0 (non-default), so it
# belongs here too. Models newer than this list are self-healed at call time:
# _chat_anthropic catches the deprecation 400, strips the sampling params,
# and retries once.
_NO_SAMPLING_PARAM_MODELS = (
    "claude-opus-4-7",
    "claude-opus-4-8",
    "claude-opus-5",
    "claude-sonnet-5",
    "claude-fable-5",
    "claude-mythos-5",
    "claude-mythos-preview",
)


def _accepts_sampling_params(model: str) -> bool:
    """False for models that reject temperature/top_p/top_k (Opus 4.7+, Claude 5 family)."""
    m = (model or "").lower()
    return not any(tag in m for tag in _NO_SAMPLING_PARAM_MODELS)


_SDK_SAMPLING_KWARGS: "Optional[bool]" = None


def _sdk_has_sampling_kwargs() -> bool:
    """True if Messages.create() still accepts temperature/top_p/top_k.

    The anthropic SDK removed all three in 1.0 — passing them raises a
    client-side TypeError, not an API 400. Older models that still accept
    sampling must then get the params via extra_body instead.
    """
    global _SDK_SAMPLING_KWARGS
    if _SDK_SAMPLING_KWARGS is None:
        try:
            import inspect
            from anthropic.resources.messages import Messages
            _SDK_SAMPLING_KWARGS = (
                "temperature" in inspect.signature(Messages.create).parameters)
        except Exception:
            _SDK_SAMPLING_KWARGS = True
    return _SDK_SAMPLING_KWARGS


# One Anthropic client per API key, not one per call. Building it inside
# _chat_anthropic meant a fresh httpx pool — and so a fresh TLS handshake to
# api.anthropic.com — on every iteration of the agent loop.
_ANTHROPIC_CLIENTS: dict = {}


def _anthropic_client(api_key: str):
    import anthropic
    client = _ANTHROPIC_CLIENTS.get(api_key)
    if client is None:
        client = anthropic.Anthropic(api_key=api_key)
        _ANTHROPIC_CLIENTS[api_key] = client
    return client


# --------------------------- Anthropic prompt caching ------------------------
#
# The agent re-sends an identical prefix on every iteration of every turn: the
# full tool catalogue (239 tools, ~45k tokens serialized) plus the system
# prompt. Uncached, that prefix is re-read from cold on each of the 2-4 LLM
# calls a single question costs, which is most of the wait before any text
# appears and most of the bill.
#
# `canonical_tools` already guarantees the tool array is byte-identical across
# calls, users and backends, which is exactly the stability a cache prefix
# needs. Breakpoints, in the order Anthropic assembles the prompt
# (tools -> system -> messages):
#
#   1. last tool      — the catalogue alone stays cached even if the system
#                       prompt is overridden per run.
#   2. system         — covers tools + system together.
#   3. last message   — rolling, so the conversation built up over the loop's
#                       iterations is reused too instead of re-read each time.
#
# Set IGVF_LLM_PROMPT_CACHE=0 to send the prefix uncached.
_CACHE_CONTROL = {"type": "ephemeral"}


def _prompt_cache_enabled() -> bool:
    return os.environ.get("IGVF_LLM_PROMPT_CACHE",
                           "1").strip().lower() not in ("0", "false", "no", "off")


def _mark_cache_breakpoint(blocks: "list[dict]") -> None:
    """Put a cache breakpoint on the last block of ``blocks``, in place."""
    if blocks and isinstance(blocks[-1], dict):
        blocks[-1]["cache_control"] = dict(_CACHE_CONTROL)


def _cacheable_system(system: Optional[str]) -> Any:
    """System prompt as a cache-marked block list (it must be blocks, not a
    bare string, to carry ``cache_control``)."""
    if not system:
        return system
    return [{"type": "text", "text": system, "cache_control": dict(_CACHE_CONTROL)}]


def _mark_last_message_cacheable(msgs: "list[dict]") -> None:
    """Rolling breakpoint at the end of the conversation, in place.

    Only string content is promoted to a block list here; content that is
    already a block list gets the marker on its final block. A message whose
    content is empty is left alone — an empty text block is a 400.
    """
    if not msgs:
        return
    last = msgs[-1]
    content = last.get("content")
    if isinstance(content, str):
        if not content.strip():
            return
        last["content"] = [{"type": "text", "text": content,
                             "cache_control": dict(_CACHE_CONTROL)}]
    elif isinstance(content, list) and content:
        _mark_cache_breakpoint(content)


def _chat_anthropic(messages, *, model, tools, max_tokens, temperature,
                     stop, on_text=None, **kwargs) -> Message:
    try:
        import anthropic
    except ImportError as e:
        raise RuntimeError(
            "Anthropic backend requires the `anthropic` package. "
            "Install with: pip install 'igvfagent[llm]'"
        ) from e
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError(
            "ANTHROPIC_API_KEY not set in the environment."
        )
    client = _anthropic_client(api_key)
    # Anthropic's Messages API has no `seed` param — strip it so the generic
    # kwargs passthrough below doesn't 400. (Determinism on Anthropic comes
    # from temperature 0 on models that still accept it.)
    kwargs.pop("seed", None)
    system, msgs = _to_anthropic_messages(messages)
    cacheable = _prompt_cache_enabled()
    if cacheable:
        _mark_last_message_cacheable(msgs)
    payload: dict = {
        "model":       model,
        "messages":    msgs,
        "max_tokens":  max_tokens,
    }
    # Sampling params (temperature/top_p/top_k) were removed on Opus 4.7+,
    # Fable 5, and Mythos 5 — sending them returns a 400. Only forward them
    # to models that still accept them.
    if _accepts_sampling_params(model):
        payload["temperature"] = temperature
    else:
        for _p in ("top_p", "top_k"):
            kwargs.pop(_p, None)
    if system:
        payload["system"] = _cacheable_system(system) if cacheable else system
    if tools:
        anth_tools = to_anthropic_tools(tools)
        if cacheable:
            _mark_cache_breakpoint(anth_tools)
        payload["tools"] = anth_tools
    if stop:
        payload["stop_sequences"] = stop
    payload.update({k: v for k, v in kwargs.items() if v is not None})

    # anthropic SDK >= 1.0 removed temperature/top_p/top_k from
    # Messages.create() itself; models that still accept them must get them
    # through extra_body, which merges into the request JSON.
    if not _sdk_has_sampling_kwargs():
        moved = {p: payload.pop(p)
                 for p in ("temperature", "top_p", "top_k") if p in payload}
        if moved:
            payload["extra_body"] = {**(payload.get("extra_body") or {}),
                                     **moved}

    def _sampling_in_payload() -> "list[str]":
        eb = payload.get("extra_body") or {}
        return [p for p in ("temperature", "top_p", "top_k")
                if p in payload or p in eb]

    def _create(pl: dict):
        """Issue the request, streaming when the caller wants deltas.

        Streaming is what puts text on the page while it is still being
        generated. With max_tokens at 16384 the blocking form left the browser
        on a spinner for the whole generation, which is most of what "the demo
        is slow" actually was. `.stream()` accumulates into the same final
        Message object, so nothing downstream changes shape.
        """
        if on_text is None:
            return client.messages.create(**pl)
        with client.messages.stream(**pl) as stream:
            for chunk in stream.text_stream:
                if chunk:
                    try:
                        on_text(chunk)
                    except Exception:                       # noqa: BLE001
                        # A rendering failure in the caller must not lose an
                        # answer that is already half-generated.
                        pass
            return stream.get_final_message()

    try:
        resp = _create(payload)
    except anthropic.BadRequestError as e:
        # Self-heal for models newer than _NO_SAMPLING_PARAM_MODELS: the API
        # rejects removed sampling params with 400 ("<param> is deprecated
        # for this model."). Strip them all and retry once.
        err = str(e).lower()
        named = any(p in err for p in ("temperature", "top_p", "top_k"))
        present = _sampling_in_payload()
        if not (named and present and "deprecated" in err):
            raise
        for p in present:
            payload.pop(p, None)
            (payload.get("extra_body") or {}).pop(p, None)
        resp = _create(payload)
    text_parts: "list[str]" = []
    tool_calls: "list[ToolCall]" = []
    for block in (resp.content or []):
        btype = getattr(block, "type", None)
        if btype == "text":
            text_parts.append(getattr(block, "text", ""))
        elif btype == "tool_use":
            tool_calls.append(ToolCall(
                id=getattr(block, "id", ""),
                name=getattr(block, "name", ""),
                arguments=getattr(block, "input", {}) or {},
            ))
    stop_map = {"end_turn":         "end_turn",
                "tool_use":         "tool_use",
                "max_tokens":       "max_tokens",
                "stop_sequence":    "stop_sequence",
                # Surface a safety-classifier refusal as itself instead of
                # collapsing it into "other" — the agent loop needs to tell
                # a refusal apart from a normal empty end_turn to react.
                "refusal":          "refusal",
                "pause_turn":       "pause_turn"}
    return Message(
        content="".join(text_parts),
        tool_calls=tool_calls,
        stop_reason=stop_map.get(getattr(resp, "stop_reason", ""), "other"),
        backend="anthropic",
        model=model,
        raw={"id": getattr(resp, "id", "")},
        usage=getattr(resp, "usage", None) and {
            "input_tokens":  getattr(resp.usage, "input_tokens", 0),
            "output_tokens": getattr(resp.usage, "output_tokens", 0),
            # Without these two a cache that silently stopped working looks
            # exactly like one that works: same answer, same input_tokens
            # field, 10x the bill.
            "cache_write_tokens":
                getattr(resp.usage, "cache_creation_input_tokens", 0) or 0,
            "cache_read_tokens":
                getattr(resp.usage, "cache_read_input_tokens", 0) or 0,
        },
    )


def _chat_openai_compat(messages, *, model, tools, max_tokens, temperature,
                          stop, base_url, api_key_env,
                          api_key_default=None, backend_label="openai",
                          **kwargs) -> Message:
    try:
        import openai
    except ImportError as e:
        raise RuntimeError(
            "OpenAI-compatible backend requires the `openai` package. "
            "Install with: pip install 'igvfagent[llm]'"
        ) from e
    api_key = os.environ.get(api_key_env) or api_key_default
    if not api_key:
        raise RuntimeError(
            f"{api_key_env} not set in the environment."
        )
    # Ollama / vLLM / local servers can take 60-180s to load a multi-GB
    # model on the first request. Honor IGVF_LLM_TIMEOUT (seconds);
    # default 600s for local backends, 120s for cloud.
    default_timeout = (600.0 if backend_label in ("ollama", "vllm", "tgi")
                       else 120.0)
    timeout = float(os.environ.get("IGVF_LLM_TIMEOUT", default_timeout))
    client = openai.OpenAI(api_key=api_key, base_url=base_url,
                            timeout=timeout)
    # OpenAI renamed ``max_tokens`` to ``max_completion_tokens`` for the
    # GPT-5 family and the o1/o3/o4 reasoning models. Classic chat-
    # completions models (gpt-4o*, gpt-4*, gpt-3.5*) still take
    # ``max_tokens``. Other OpenAI-compatible servers (Groq, Together,
    # DeepInfra, Ollama, vLLM, TGI, …) universally accept ``max_tokens``
    # regardless of model name. So we only emit ``max_completion_tokens``
    # for the new OpenAI generations.
    model_lc = (model or "").lower()
    uses_new_param = (
        backend_label in ("openai", "codex", "custom")
        and (model_lc.startswith("gpt-5")
             or model_lc.startswith("o1")
             or model_lc.startswith("o3")
             or model_lc.startswith("o4"))
    )
    payload: dict = {
        "model":       model,
        "messages":    _to_openai_messages(messages),
        "temperature": temperature,
    }
    if uses_new_param:
        payload["max_completion_tokens"] = max_tokens
        # The whole reasoning-class family (gpt-5*, o1*, o3*, o4*) rejects
        # any temperature other than the default 1.0. Coerce silently
        # rather than 400-erroring the user; the dispatcher's tool-call
        # argument coercion is what matters for reproducibility, not the
        # sampling-distribution flatness of the LLM's own free-form text.
        if (model_lc.startswith(("gpt-5", "o1", "o3", "o4"))
                and temperature != 1.0):
            payload["temperature"] = 1.0
    else:
        payload["max_tokens"] = max_tokens
    if tools:
        serialized = to_openai_tools(tools)
        # OpenAI's Chat Completions API caps the ``tools`` array at 128
        # entries (post-GPT-5 generation). IGVFagent has > 128 registered
        # tools today, so we trim — preferring star-prefixed (★) tools
        # first because those are the hand-curated LLM-tool-selection
        # core; the rest fill the remaining slots deterministically.
        OPENAI_TOOLS_MAX = 128
        if (backend_label in ("openai", "codex", "custom")
                and len(serialized) > OPENAI_TOOLS_MAX):
            def _is_starred(entry: dict) -> bool:
                desc = (entry.get("function", {}) or {}).get("description", "") or ""
                return desc.startswith("★") or desc.startswith("★ ")
            starred = [e for e in serialized if _is_starred(e)]
            unstarred = [e for e in serialized if not _is_starred(e)]
            trimmed = (starred + unstarred)[:OPENAI_TOOLS_MAX]
            logging.warning(
                "OpenAI 128-tool cap: trimmed %d -> %d (%d starred kept, "
                "%d unstarred kept, %d dropped)",
                len(serialized), len(trimmed),
                sum(1 for e in trimmed if _is_starred(e)),
                sum(1 for e in trimmed if not _is_starred(e)),
                len(serialized) - len(trimmed))
            serialized = trimmed
        payload["tools"] = serialized
        payload["tool_choice"] = "auto"
    if stop:
        payload["stop"] = stop
    payload.update({k: v for k, v in kwargs.items() if v is not None})

    resp = client.chat.completions.create(**payload)
    choice = resp.choices[0] if resp.choices else None
    if choice is None:
        return Message(stop_reason="error", backend=backend_label,
                       model=model, raw={"empty_response": True})
    msg = choice.message
    tool_calls: "list[ToolCall]" = []
    for tc in (getattr(msg, "tool_calls", None) or []):
        try:
            args = json.loads(tc.function.arguments or "{}")
        except json.JSONDecodeError:
            args = {"_raw_arguments": tc.function.arguments}
        tool_calls.append(ToolCall(
            id=tc.id, name=tc.function.name, arguments=args,
        ))
    finish = getattr(choice, "finish_reason", "") or ""
    stop_map = {"stop":          "end_turn",
                "tool_calls":    "tool_use",
                "function_call": "tool_use",
                "length":        "max_tokens"}
    usage = getattr(resp, "usage", None)
    return Message(
        content=getattr(msg, "content", "") or "",
        tool_calls=tool_calls,
        stop_reason=stop_map.get(finish, "other"),
        backend=backend_label,
        model=model,
        raw={"id": getattr(resp, "id", "")},
        usage=usage and {
            "input_tokens":  getattr(usage, "prompt_tokens", 0),
            "output_tokens": getattr(usage, "completion_tokens", 0),
        },
    )


# --------------------------- Backend table ---------------------------------

# Each entry: backend label -> default base_url + default api_key env var.
# Override any of these via ``IGVF_LLM_BASE_URL`` / ``IGVF_LLM_API_KEY_ENV``.
try:
    from _endpoints import resolve as _resolve_endpoint
except ImportError:  # pragma: no cover
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent))
    from _endpoints import resolve as _resolve_endpoint

_BACKENDS = {
    "anthropic": None,                                     # special-cased
    "openai":    {"base_url": _resolve_endpoint("llm_openai", "OPENAI_BASE_URL"),
                   "api_key_env": "OPENAI_API_KEY"},
    "codex":     {"base_url": _resolve_endpoint("llm_openai", "IGVF_CODEX_BASE_URL"),
                   "api_key_env": "OPENAI_API_KEY"},
    "ollama":    {"base_url": os.environ.get("OLLAMA_HOST_BASE",
                                                "http://localhost:11434/v1"),
                   "api_key_env": "OLLAMA_API_KEY",
                   "api_key_default": "ollama"},
    "vllm":      {"base_url": os.environ.get("IGVF_LLM_BASE_URL",
                                                "http://localhost:8000/v1"),
                   "api_key_env": "IGVF_LLM_API_KEY",
                   "api_key_default": "EMPTY"},
    "tgi":       {"base_url": os.environ.get("IGVF_LLM_BASE_URL",
                                                "http://localhost:3000/v1"),
                   "api_key_env": "IGVF_LLM_API_KEY",
                   "api_key_default": "EMPTY"},
    "groq":      {"base_url": _resolve_endpoint("llm_groq", "GROQ_BASE_URL"),
                   "api_key_env": "GROQ_API_KEY"},
    "together":  {"base_url": _resolve_endpoint("llm_together", "TOGETHER_BASE_URL"),
                   "api_key_env": "TOGETHER_API_KEY"},
    "deepinfra": {"base_url": _resolve_endpoint("llm_deepinfra", "DEEPINFRA_BASE_URL"),
                   "api_key_env": "DEEPINFRA_API_KEY"},
    "huggingface": {
        "base_url": (os.environ.get("IGVF_LLM_BASE_URL")
                       or _resolve_endpoint("llm_hf_inference", "HF_BASE_URL")),
        "api_key_env": "HF_TOKEN",
    },
    # User-provided OpenAI-compatible endpoint
    "custom":    {"base_url": os.environ.get("IGVF_LLM_BASE_URL", ""),
                   "api_key_env": os.environ.get("IGVF_LLM_API_KEY_ENV",
                                                   "IGVF_LLM_API_KEY")},
}


# --------------------------- Cross-backend consistency ----------------------

# The strictest backend (OpenAI Chat Completions) caps the tools array at 128.
# To guarantee every backend sees the *same* tool set — so a query resolves to
# the same tool no matter which LLM drives the loop — we apply one canonical,
# deterministically-ordered selection to ALL backends, not just OpenAI.
# Override with IGVF_LLM_MAX_TOOLS (e.g. raise it if you only ever use
# Anthropic/Ollama and want all 141 exposed — at the cost of OpenAI parity).
_DEFAULT_MAX_TOOLS = 128


def _is_starred_tool(entry: dict) -> bool:
    """A ★-prefixed description marks a hand-curated core tool."""
    desc = entry.get("description", "") or ""
    # tolerate already-serialized OpenAI/Anthropic shapes too
    if not desc and isinstance(entry.get("function"), dict):
        desc = entry["function"].get("description", "") or ""
    return desc.lstrip().startswith("★")


def canonical_tools(tools: "Optional[list[dict]]",
                     max_tools: Optional[int] = None) -> "Optional[list[dict]]":
    """Return one backend-independent, deterministically-ordered tool subset.

    Ordering: starred (core) tools first, then the rest; each group sorted by
    ``name``. Capped to ``max_tools`` (default 128 / IGVF_LLM_MAX_TOOLS). The
    output is identical regardless of the calling backend, which is what makes
    tool selection reproducible across Claude / Codex / Qwen / etc.
    """
    if not tools:
        return tools
    if max_tools is None:
        try:
            max_tools = int(os.environ.get("IGVF_LLM_MAX_TOOLS",
                                            _DEFAULT_MAX_TOOLS))
        except ValueError:
            max_tools = _DEFAULT_MAX_TOOLS

    def _key(e: dict) -> str:
        return (e.get("name")
                or (e.get("function", {}) or {}).get("name", "")
                or "")

    def _family(name: str) -> str:
        """Tool family: the prefix before the first underscore."""
        return name.split("_", 1)[0] if "_" in name else name

    def _round_robin(entries: "list[dict]") -> "list[dict]":
        """Interleave by family, so truncation costs breadth, not a subsystem.

        Sorting purely by name and cutting the tail deletes whole families
        alphabetically. Measured on the deployment at a cap of 200/234, the
        casualties were share_* (5 of 6), spatial_hic_* (7 of 9), tabula_*
        (7 of 11) and starr_* (4 of 5) -- so a question about SHARE-seq or
        Spatial-Hi-C met a subsystem with one tool left, while families
        early in the alphabet kept every one of theirs.

        Interleaving is still fully deterministic (families sorted, members
        sorted within each), which is what parity across backends requires.
        """
        buckets: "dict[str, list[dict]]" = {}
        for e in sorted(entries, key=_key):
            buckets.setdefault(_family(_key(e)), []).append(e)
        out: "list[dict]" = []
        fams = sorted(buckets)
        while any(buckets[f] for f in fams):
            for f in fams:
                if buckets[f]:
                    out.append(buckets[f].pop(0))
        return out

    starred = sorted((e for e in tools if _is_starred_tool(e)), key=_key)
    unstarred = [e for e in tools if not _is_starred_tool(e)]
    # Starred tools stay in name order: they are the curated core and all of
    # them are expected to survive any sane cap. Only the remainder, which is
    # what a cap actually eats into, is interleaved.
    ordered = starred + _round_robin(unstarred)
    if len(ordered) > max_tools:
        dropped = sorted(_key(e) for e in ordered[max_tools:])
        logger.warning(
            "canonical_tools: exposing %d/%d tools identically across all "
            "backends (%d starred kept); dropped for parity: %s",
            max_tools, len(ordered), min(len(starred), max_tools),
            ", ".join(dropped[:20]) + (" …" if len(dropped) > 20 else ""))
        ordered = ordered[:max_tools]
    return ordered


# --------------------------- Public API ------------------------------------

def chat(
    messages: "list[dict]",
    *,
    backend: Optional[str] = None,
    model: Optional[str] = None,
    tools: Optional["list[dict]"] = None,
    max_tokens: int = 4096,
    temperature: float = 0.0,
    stop: Optional["list[str]"] = None,
    seed: Optional[int] = None,
    on_text: Optional[Callable[[str], None]] = None,
    **kwargs,
) -> Message:
    """Backend-neutral chat completion.

    ``on_text`` receives text deltas as they are generated (Anthropic only for
    now; other backends ignore it and return the finished message as before).

    Returns a :class:`Message` with normalized ``content``, ``tool_calls``,
    and ``stop_reason``. Raises ``RuntimeError`` if the resolved backend's
    SDK is missing or the credentials env var is unset.

    For cross-backend consistency the tool list is reduced to one canonical
    deterministic subset (see :func:`canonical_tools`) before dispatch, and a
    decoding ``seed`` (default ``IGVF_LLM_SEED``) is forwarded to every backend
    that supports it (all OpenAI-compatible servers incl. Ollama/vLLM).
    """
    # One identical tool set for every backend.
    tools = canonical_tools(tools)
    # Deterministic decoding seed (OpenAI-compatible backends honor `seed`).
    if seed is None:
        _env_seed = os.environ.get("IGVF_LLM_SEED")
        if _env_seed not in (None, ""):
            try:
                seed = int(_env_seed)
            except ValueError:
                seed = None
    if seed is not None:
        kwargs.setdefault("seed", seed)
    bk = _resolve_backend(backend, model)

    # Local-only mode. "Local tool execution" and "local inference" are
    # different guarantees, and conflating them is how prompts, metadata and
    # findings reach a vendor from a deployment believed to be self-contained.
    # With IGVF_LOCAL_ONLY=1 a hosted backend is refused outright rather than
    # silently used — a warning would be read past.
    if os.environ.get("IGVF_LOCAL_ONLY", "").strip().lower() in ("1", "true", "yes", "on"):
        if bk not in _LOCAL_BACKENDS:
            raise RuntimeError(
                f"IGVF_LOCAL_ONLY=1 refuses the hosted backend {bk!r}: the "
                f"prompt — which accumulates tool output over a session — "
                f"would be transmitted to a third-party provider.\n"
                f"Locally served backends: {', '.join(sorted(_LOCAL_BACKENDS))}.\n"
                f"Note this bounds INFERENCE only; federated data access still "
                f"contacts public archives under every backend, so query terms "
                f"leave the machine regardless (see Docs/THREAT_MODEL.md).")
    chosen_model = (
        model
        or os.environ.get("IGVF_LLM_MODEL")
        or _DEFAULT_MODELS.get(bk, "qwen3:8b")
    )
    logger.info("LLM chat: backend=%s model=%s tools=%d msgs=%d",
                 bk, chosen_model, len(tools or []), len(messages))

    if bk == "anthropic":
        return _chat_anthropic(messages, model=chosen_model, tools=tools,
                                max_tokens=max_tokens,
                                temperature=temperature, stop=stop,
                                on_text=on_text, **kwargs)
    if bk == "claude_cli":
        return _chat_claude_cli(messages, model=chosen_model, tools=tools,
                                  max_tokens=max_tokens,
                                  temperature=temperature, stop=stop,
                                  **kwargs)
    if bk == "codex_cli":
        return _chat_codex_cli(messages, model=chosen_model, tools=tools,
                                 max_tokens=max_tokens,
                                 temperature=temperature, stop=stop,
                                 **kwargs)
    cfg = _BACKENDS.get(bk)
    if not cfg or not cfg.get("base_url"):
        raise RuntimeError(
            f"Unknown or unconfigured backend: {bk}. "
            f"Set IGVF_LLM_BASE_URL / IGVF_LLM_API_KEY_ENV for `custom`."
        )
    return _chat_openai_compat(
        messages, model=chosen_model, tools=tools,
        max_tokens=max_tokens, temperature=temperature, stop=stop,
        base_url=cfg["base_url"],
        api_key_env=cfg["api_key_env"],
        api_key_default=cfg.get("api_key_default"),
        backend_label=bk,
        **kwargs,
    )


def list_backends() -> "list[str]":
    extras = ["claude_cli", "codex_cli"]
    return ["anthropic"] + sorted(set(_BACKENDS) - {"anthropic"}) + extras


def claude_cli_available() -> "tuple[bool, str]":
    """Whether the `claude` (Claude Code) CLI is reachable on PATH.
    Returns (ok, version_or_error_msg)."""
    import shutil
    import subprocess
    if not shutil.which("claude"):
        return False, "`claude` not found on PATH"
    try:
        out = subprocess.run(["claude", "--version"], capture_output=True,
                                text=True, timeout=10, check=False)
        if out.returncode == 0:
            return True, (out.stdout or out.stderr).strip()
        return False, (out.stderr or out.stdout).strip()
    except Exception as e:
        return False, str(e)


def describe_backend(name: str) -> dict:
    name = name.lower()
    if name == "anthropic":
        return {"name": name, "sdk": "anthropic",
                "api_key_env": "ANTHROPIC_API_KEY"}
    cfg = _BACKENDS.get(name)
    if not cfg:
        return {"name": name, "_unknown": True}
    return {"name": name, "sdk": "openai", **cfg}


# --------------------------- Claude Code CLI backend ----------------------

_CLAUDE_CLI_TOOL_PROMPT = """\
You are IGVFagent, an autonomous research assistant. The user is running
this conversation through the Claude Code CLI; we use the CLI as a
backend, not as an interactive coding agent. Do NOT use Claude Code's
own file/edit/bash tools — they will not be executed in this context.

You have access ONLY to the tools listed below. To call one, respond
with EXACTLY this XML block (you may emit multiple in one turn):

<tool_call>
  <name>tool_name_here</name>
  <arguments>{"argname": "value", "another": 42}</arguments>
</tool_call>

When you have enough information to answer the user, respond with:

<final_answer>
Your concise Markdown response. Cite report / manifest paths from the
tool results, flag caveats, and suggest one concrete next CLI call.
</final_answer>

You may emit either tool_calls OR a final_answer per turn, not both.
The system will execute the tool calls and re-prompt you with the
results.

# Tools available

{tools_block}
"""

_PROTOCOL_VIOLATION_NOTE = (
    "**This run did not complete.** The model did not use the tool protocol: "
    "it returned no usable tool call twice, or claimed that tools which are "
    "registered and callable were unavailable. Nothing below "
    "is the result of analysis, and any claim in it about tools being "
    "unavailable is the model's own inference, not a real backend fault. "
    "Re-run the query, or select a different model in the sidebar.\n\n"
    "---\n\n"
)

# Phrases a model uses when it has decided, wrongly, that the backend is
# broken. The harness's own wording for an unregistered tool is "No such tool
# available", so a model reproducing that phrase is quoting an error it never
# received -- there is no code path that hands it one.
_FABRICATED_OUTAGE_RE = __import__("re").compile(
    r"no such tool available"
    r"|no such tool\b"
    r"|tool[- ]availability (problem|issue)"
    r"|tools?\s+(?:are|is|were|was)?\s*(?:not|n't)\s+(?:available|responding)"
    r"|every subsequent tool call .{0,40}failed"
    r"|(?:session|registry) issue on my end",
    __import__("re").IGNORECASE)


def _fabricated_outage(text: str, tools) -> "list[str]":
    """Tool names a final answer calls unavailable that are in fact registered.

    The protocol retry below only fired when a response had NEITHER a
    <tool_call> NOR a <final_answer>. That misses the failure that actually
    reaches a user: the model emits a perfectly well-formed <final_answer>
    whose content is "every tool call returned No such tool available", and
    the loop delivers it as a finished result. On IGVFDS6464SOVZ that answer
    named base_editing_screen_discover, base_editing_screen_analyze,
    crispr_screen_discover and portal_get -- all four present in the running
    registry, all four callable from the CLI at that moment.

    Returning the offending names rather than a bool keeps the retry note
    specific: telling a model that `portal_get` exists is far more corrective
    than telling it its protocol was wrong.
    """
    if not text or not _FABRICATED_OUTAGE_RE.search(text):
        return []
    known = set()
    for t in tools or []:
        name = getattr(t, "name", None) or (
            t.get("name") if isinstance(t, dict) else None)
        if name:
            known.add(name)
    # Only names the text actually mentions, so a generic grumble about some
    # other system does not trigger a retry.
    return sorted(n for n in known
                  if __import__("re").search(r"\b" + n + r"\b", text))


_CLAUDE_CLI_TOOL_CALL_RE = __import__("re").compile(
    r"<tool_call>\s*"
    r"<name>\s*(?P<name>[A-Za-z0-9_\-]+)\s*</name>\s*"
    r"<arguments>\s*(?P<args>.*?)\s*</arguments>\s*"
    r"</tool_call>",
    __import__("re").DOTALL,
)
_CLAUDE_CLI_FINAL_RE = __import__("re").compile(
    r"<final_answer>\s*(.*?)\s*</final_answer>",
    __import__("re").DOTALL,
)


def _claude_cli_render_tools(tools) -> str:
    out = []
    for t in tools or []:
        params = t.get("parameters") or {}
        props = params.get("properties") or {}
        required = set(params.get("required") or [])
        param_lines = []
        for k, v in props.items():
            tag = "" if k in required else "  (optional)"
            ty = v.get("type", "?") if isinstance(v, dict) else "?"
            desc = v.get("description", "") if isinstance(v, dict) else ""
            param_lines.append(f"    - {k}: {ty}{tag}  {desc[:80]}")
        out.append(
            f"## {t['name']}\n{t.get('description','').strip()}\n"
            + ("\nParameters:\n" + "\n".join(param_lines) if param_lines
                else "\n(no parameters)\n")
        )
    return "\n\n".join(out)


def _claude_cli_serialize_messages(messages) -> str:
    """Compact textual rendering of the conversation for the CLI prompt."""
    parts = []
    for m in messages:
        role = m.get("role", "user")
        if role == "system":
            continue   # system goes into the prompt template separately
        if role == "tool":
            parts.append(
                f"[Tool result for call {m.get('tool_call_id','?')}]\n"
                f"{m.get('content','')}\n"
            )
            continue
        if role == "assistant":
            tcs = m.get("tool_calls") or []
            if tcs:
                for tc in tcs:
                    name = tc["name"] if isinstance(tc, dict) else tc.name
                    args = (tc["arguments"] if isinstance(tc, dict)
                            else tc.arguments)
                    parts.append(
                        f"[Assistant called tool]\n<tool_call>\n"
                        f"  <name>{name}</name>\n"
                        f"  <arguments>{json.dumps(args, default=str)}"
                        f"</arguments>\n</tool_call>\n"
                    )
            if m.get("content"):
                parts.append(f"[Assistant said]\n{m['content']}\n")
            continue
        parts.append(f"[User]\n{m.get('content','')}\n")
    return "\n".join(parts)


def _xml_cli_build_prompt(messages, tools) -> str:
    """Shared prompt builder used by both the claude_cli and codex_cli
    subprocess backends."""
    system_chunks = [m.get("content", "") for m in messages
                     if m.get("role") == "system"]
    user_system = "\n\n".join(c for c in system_chunks if c).strip()
    tools_block = _claude_cli_render_tools(tools) if tools else \
                  "_(no tools — produce <final_answer> directly)_"
    framework = _CLAUDE_CLI_TOOL_PROMPT.replace("{tools_block}", tools_block)
    convo = _claude_cli_serialize_messages(messages).strip()
    return (
        framework
        + ("\n\n# Project-specific guidance\n\n" + user_system
            if user_system else "")
        + "\n\n# Conversation so far\n\n"
        + (convo or "_(no prior turns)_")
        + "\n\nRespond now."
    )


def _xml_cli_parse_response(text: str, prefix: str
                              ) -> "tuple[str, list[ToolCall]]":
    """Parse the XML-tool-call framework response into (content,
    tool_calls)."""
    tool_calls: "list[ToolCall]" = []
    for i, m in enumerate(_CLAUDE_CLI_TOOL_CALL_RE.finditer(text)):
        name = m.group("name").strip()
        raw_args = m.group("args").strip()
        try:
            args = json.loads(raw_args) if raw_args else {}
        except json.JSONDecodeError:
            args = {"_raw_arguments": raw_args}
        tool_calls.append(ToolCall(
            id=f"{prefix}_{i+1}", name=name,
            arguments=args if isinstance(args, dict)
                      else {"_raw": str(args)},
        ))
    final_match = _CLAUDE_CLI_FINAL_RE.search(text)
    if final_match:
        content = final_match.group(1).strip()
    elif tool_calls:
        content = ""
    else:
        content = text.strip()
    return content, tool_calls


def _chat_claude_cli(messages, *, model, tools, max_tokens, temperature,
                      stop, **kwargs) -> Message:
    """Subprocess-out to `claude --print` and parse a ReAct-style
    response. Trade-offs vs the native Anthropic API:

      + Re-uses the user's existing Claude Code login (no separate
        ANTHROPIC_API_KEY required).
      + Whatever model Claude Code is configured for is what's used.
      - 5–15s per turn of CLI subprocess + auth overhead.
      - Each call ships the full conversation history (no session
        reuse on the CLI side).
      - Tool calls are text-parsed (XML), not native function calling;
        a malformed model response can drop a turn.
    """
    import shutil
    import subprocess

    if not shutil.which("claude"):
        raise RuntimeError(
            "Claude Code CLI (`claude`) is not on PATH. Install from "
            "https://docs.claude.com/en/docs/claude-code or via "
            "`npm i -g @anthropic-ai/claude-code`."
        )

    prompt = _xml_cli_build_prompt(messages, tools)
    # The prompt goes on STDIN, never in argv. Passing it as an argument
    # works until a conversation grows -- then the whole transcript exceeds
    # the kernel's argv limit and execve fails with
    #     [Errno 7] Argument list too long: 'claude'
    # which killed agent runs mid-flight around iteration 11, after the
    # expensive tool calls had already been made. `claude --print` reads the
    # prompt from stdin when given no positional argument, and stdin has no
    # such limit.
    cmd = ["claude", "--print", "--output-format", "text"]
    # Disable Claude Code's OWN tools. We use the CLI as a text completion
    # engine: the tools it must call are IGVFagent's, described in the
    # prompt and invoked by emitting <tool_call> XML that we parse back. The
    # prompt already tells the model not to reach for Bash/Read/Edit, so
    # this only enforces what it asks for.
    #
    # Defence in depth, NOT a proven bug fix. It was tried as a fix for the
    # IGVFDS5997IVEM run, where the model emitted no XML at all and claimed
    # every tool returned "No such tool available". That hypothesis did not
    # survive testing: with the real 249-tool, 138 KB prompt in this
    # container, claude-sonnet-5 emitted valid <tool_call> XML both with
    # native tools enabled and disabled. The protocol lapse is intermittent
    # model behaviour, and what actually contains it is the retry below.
    #
    # The reason to keep this anyway is exposure: the nested CLI was being
    # handed Bash and Write inside the container that holds the Portal API
    # key. Set IGVF_CLI_NATIVE_TOOLS=1 to restore the old behaviour.
    if os.environ.get("IGVF_CLI_NATIVE_TOOLS", "0") != "1":
        cmd.extend(["--tools", ""])
    if model and model.strip():
        cmd.extend(["--model", model.strip()])

    timeout = float(os.environ.get("IGVF_LLM_TIMEOUT", "600"))
    try:
        result = subprocess.run(
            cmd, input=prompt, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"Claude Code CLI timed out after {timeout:.0f}s "
            f"(set IGVF_LLM_TIMEOUT to override)."
        ) from e

    if result.returncode != 0:
        err = (result.stderr or result.stdout or "").strip()[:600]
        raise RuntimeError(f"`claude` exited with status "
                            f"{result.returncode}: {err}")

    text = result.stdout or ""
    content, tool_calls = _xml_cli_parse_response(text, prefix="cc")

    # Neither a <tool_call> nor a <final_answer> means the model ignored the
    # protocol. The parser falls back to returning the raw prose, which the
    # agent loop then accepts as a finished answer -- so a run that did no
    # work at all was reported as "stop complete" with the model's own
    # invented explanation as the result. Retry once, saying plainly what
    # was missing, before letting that through.
    fabricated = _fabricated_outage(content, tools) if not tool_calls else []
    if fabricated:
        logger.warning("claude_cli: final answer claims %d registered tool(s) "
                       "are unavailable (%s); retrying once",
                       len(fabricated), ", ".join(fabricated[:4]))
    if (not tool_calls
            and (fabricated or not _CLAUDE_CLI_FINAL_RE.search(text))):
        if not fabricated:
            logger.warning("claude_cli: response had no <tool_call> and no "
                           "<final_answer>; retrying once with a corrective note")
        retry_prompt = (
            prompt
            + ("\n\n# These tools exist\n\n"
                "Your previous response stated that tools were unavailable "
                "and named: " + ", ".join(fabricated) + ". Every one of those "
                "is registered and callable in this session. No tool result "
                "saying 'No such tool available' was ever sent to you -- "
                "there is no code path that produces one, so that text was "
                "invented. Call the tool by emitting the <tool_call> XML "
                "block. Do not report a backend outage.\n"
                if fabricated else "")
            + "\n\n# Protocol reminder\n\n"
              "Your previous response contained neither a <tool_call> nor a "
              "<final_answer> block, so it could not be used. You do not "
              "have Claude Code's own tools here and must not try to invoke "
              "them; the only tools that exist are the ones listed above, "
              "and the only way to call one is to emit the <tool_call> XML "
              "block exactly as specified. If a tool call appears to fail, "
              "that is reported to you as a tool result -- it is not a "
              "backend outage. Respond now with either one or more "
              "<tool_call> blocks or a single <final_answer> block."
        )
        try:
            retry = subprocess.run(
                cmd, input=retry_prompt, capture_output=True, text=True,
                timeout=timeout, check=False,
            )
            if retry.returncode == 0:
                rtext = retry.stdout or ""
                rcontent, rcalls = _xml_cli_parse_response(rtext, prefix="cc")
                if rcalls or (_CLAUDE_CLI_FINAL_RE.search(rtext)
                               and not _fabricated_outage(rcontent, tools)):
                    text, content, tool_calls = rtext, rcontent, rcalls
                else:
                    # Say so rather than passing prose off as an answer.
                    return Message(
                        content=_PROTOCOL_VIOLATION_NOTE + rcontent,
                        tool_calls=[], stop_reason="protocol_violation",
                        backend="claude_cli",
                        model=model or "(claude-code default)",
                        raw={"stdout_len": len(rtext), "retried": True},
                    )
        except subprocess.TimeoutExpired:
            logger.warning("claude_cli: protocol retry timed out")

    return Message(
        content=content,
        tool_calls=tool_calls,
        stop_reason="tool_use" if tool_calls else "end_turn",
        backend="claude_cli",
        model=model or "(claude-code default)",
        raw={"stdout_len": len(text)},
    )


# --------------------------- Codex CLI backend ----------------------------

def codex_cli_available() -> "tuple[bool, str]":
    """Whether the `codex` (OpenAI Codex CLI) binary is on PATH.
    Returns (ok, version_or_error_msg)."""
    import shutil
    import subprocess
    if not shutil.which("codex"):
        return False, "`codex` not found on PATH"
    try:
        out = subprocess.run(
            ["codex", "--version"], capture_output=True,
            text=True, timeout=10, check=False,
        )
        if out.returncode == 0:
            return True, (out.stdout or out.stderr).strip()
        return False, (out.stderr or out.stdout).strip()
    except Exception as e:
        return False, str(e)


def _chat_codex_cli(messages, *, model, tools, max_tokens, temperature,
                     stop, **kwargs) -> Message:
    """Subprocess-out to `codex exec` with the same XML tool-call
    framework as claude_cli. Mirrors the trade-off profile: reuses the
    user's Codex CLI login (no separate OPENAI_API_KEY needed) at the
    cost of subprocess + XML-parsing overhead vs the native OpenAI API.
    """
    import shutil
    import subprocess

    if not shutil.which("codex"):
        raise RuntimeError(
            "Codex CLI (`codex`) is not on PATH. Install with "
            "`npm i -g @openai/codex` (or follow "
            "https://github.com/openai/codex)."
        )

    prompt = _xml_cli_build_prompt(messages, tools)

    # `codex exec` runs in non-interactive mode already, so it never
    # blocks on a TTY approval prompt. Older CLIs took
    # `--ask-for-approval never`, but codex >= 0.11 removed that flag
    # from the `exec` subcommand ("unexpected argument"), so we no
    # longer pass it. Our prompt tells Codex to emit XML rather than use
    # its own tools, so no approval is needed.
    cmd = ["codex", "exec"]
    if model and model.strip():
        cmd.extend(["--model", model.strip()])
    # Read the prompt from stdin to avoid argv-length limits.
    cmd.append("-")  # convention: dash means "read prompt from stdin"

    timeout = float(os.environ.get("IGVF_LLM_TIMEOUT", "600"))
    try:
        result = subprocess.run(
            cmd, input=prompt, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"Codex CLI timed out after {timeout:.0f}s "
            f"(set IGVF_LLM_TIMEOUT to override)."
        ) from e

    # Some Codex CLI versions don't accept `-` for stdin or differ on
    # the approval flag. Fall back to passing the prompt as a positional
    # argument (truncated diagnostic on persistent failure).
    if result.returncode != 0 and ("unrecognized" in (result.stderr or "")
                                    or "invalid value" in (result.stderr or "")
                                    or "unexpected argument" in (result.stderr or "")
                                    or "expected one of" in (result.stderr or "")):
        cmd_fallback = ["codex", "exec"]
        if model and model.strip():
            cmd_fallback.extend(["--model", model.strip()])
        cmd_fallback.append(prompt)
        result = subprocess.run(
            cmd_fallback, capture_output=True, text=True,
            timeout=timeout, check=False,
        )

    if result.returncode != 0:
        err = (result.stderr or result.stdout or "").strip()[:600]
        raise RuntimeError(f"`codex exec` exited with status "
                            f"{result.returncode}: {err}")

    text = result.stdout or ""
    content, tool_calls = _xml_cli_parse_response(text, prefix="cx")
    return Message(
        content=content,
        tool_calls=tool_calls,
        stop_reason="tool_use" if tool_calls else "end_turn",
        backend="codex_cli",
        model=model or "(codex-cli default)",
        raw={"stdout_len": len(text)},
    )


def list_ollama_models(base_url: Optional[str] = None,
                         timeout: float = 5.0) -> "list[dict]":
    """Query the configured Ollama daemon for installed models.

    Returns a list of dicts with at least ``name`` (the tag a user would
    pass via ``--model``). Returns an empty list (not raising) on any
    network error so the introspection CLI degrades gracefully.
    """
    import json as _json
    import urllib.error as _urlerr
    import urllib.request as _urlreq

    base = (base_url
            or os.environ.get("OLLAMA_HOST_BASE")
            or _BACKENDS["ollama"]["base_url"])
    # /api/tags lives on the native Ollama port, alongside /v1/* OpenAI compat
    tags_url = base.rstrip("/").replace("/v1", "") + "/api/tags"
    try:
        req = _urlreq.Request(tags_url,
                                headers={"Accept": "application/json"})
        with _urlreq.urlopen(req, timeout=timeout) as resp:
            data = _json.loads(resp.read())
    except (_urlerr.URLError, _urlerr.HTTPError, OSError):
        return []
    out = []
    for m in (data.get("models") or []):
        if not isinstance(m, dict):
            continue
        size = m.get("size") or 0
        out.append({
            "name":     m.get("name") or m.get("model") or "",
            "size":     size,
            "size_gb":  round(size / (1024 ** 3), 2) if size else None,
            "modified": m.get("modified_at") or "",
            "family":   ((m.get("details") or {}).get("family") or ""),
        })
    return out


# Curated Ollama-pullable models with approximate sizes. Used by the UI
# "Download more models" picker so users can grab a model without leaving
# the page. Names match Ollama's library tags.
OLLAMA_LIBRARY = [
    # name                          approx GB   notes
    ("qwen3:0.6b",                       0.5,  "smallest, weak tool calls"),
    ("qwen3:4b",                         2.7,  "fast, decent tool calls"),
    ("qwen3:8b",                         5.2,  "project default; balanced"),
    ("qwen3:14b",                        9.0,  "stronger reasoning"),
    ("qwen3:32b",                       20.0,  "needs ≥ 32 GB free"),
    ("qwen2.5:1.5b",                     1.0,  "tiny baseline"),
    ("qwen2.5:7b",                       4.7,  "Qwen 2.5 standard"),
    ("qwen2.5-coder:7b",                 4.7,  "code-tuned"),
    ("qwen2.5-coder:14b",                9.0,  "code-tuned, stronger"),
    ("gemma3:4b",                        3.3,  "Google Gemma 3"),
    ("gemma3:12b",                       8.1,  ""),
    ("gemma3:27b",                      17.0,  ""),
    ("llama3.1:8b",                      4.7,  "Meta Llama 3.1"),
    ("llama3.1:70b",                    40.0,  "needs ≥ 64 GB free"),
    ("llama3.2:3b",                      2.0,  "tiny Llama 3.2"),
    ("mistral:7b",                       4.1,  "Mistral baseline"),
    ("mistral-small:22b",               13.0,  "Mistral small instruct"),
    ("deepseek-r1:7b",                   4.7,  "reasoning model"),
    ("deepseek-r1:14b",                  9.0,  ""),
    ("phi4:14b",                         9.1,  "Microsoft Phi 4"),
    ("codellama:13b",                    7.4,  "code-tuned Llama"),
]


# Per-backend curated model lists used by the LM Studio-style sidebar
# dropdowns. Real model availability is checked separately (e.g.
# Anthropic API access requires the key + tier).
# The four current Claude tiers. Haiku 4.5 is the only one that still
# accepts temperature/top_p/top_k (see _NO_SAMPLING_PARAM_MODELS), which
# makes it the cheap tier for high-volume runs. Earlier generations are
# reachable through the picker's "(custom...)" field for comparison runs.
ANTHROPIC_MODELS = [
    "claude-opus-5",
    "claude-sonnet-5",
    # Fable 5.1 replaces Fable 5 in the picker. The id was read from
    # GET /v1/models rather than guessed -- it returns
    # "claude-fable-5-1" / "Claude Fable 5.1", and inventing an id is how a
    # picker offers a model that 404s. claude-fable-5 still works if pinned
    # through "(custom...)"; it is simply not offered.
    "claude-fable-5-1",
    "claude-haiku-4-5",
]

OPENAI_MODELS = [
    "gpt-5",
    "gpt-5-codex",
    "gpt-5-mini",
    "gpt-4o",
    "gpt-4o-mini",
    "o3",
    "o3-mini",
    "o1-preview",
    "o1-mini",
]


def pull_ollama_model(model_name: str,
                        base_url: Optional[str] = None,
                        timeout: float = 1800.0):
    """Stream-pull an Ollama model. Yields ``(status, percent, total_bytes,
    completed_bytes)`` tuples as progress arrives so the UI can render a
    progress bar without blocking.

    Returns the final status (``success`` or ``error``) at the end.
    """
    import json as _json
    import urllib.request as _urlreq

    base = (base_url
            or os.environ.get("OLLAMA_HOST_BASE")
            or _BACKENDS["ollama"]["base_url"])
    pull_url = base.rstrip("/").replace("/v1", "") + "/api/pull"
    payload = _json.dumps({"name": model_name, "stream": True}).encode("utf-8")
    req = _urlreq.Request(
        pull_url, data=payload,
        headers={"Content-Type": "application/json",
                 "Accept": "application/x-ndjson"},
        method="POST",
    )
    with _urlreq.urlopen(req, timeout=timeout) as resp:
        for line in resp:
            if not line.strip():
                continue
            try:
                evt = _json.loads(line)
            except _json.JSONDecodeError:
                continue
            yield (
                evt.get("status", ""),
                (float(evt.get("completed", 0)) /
                 float(evt.get("total", 1)) * 100.0
                 if evt.get("total") else 0.0),
                int(evt.get("total") or 0),
                int(evt.get("completed") or 0),
                bool(evt.get("error")),
            )


__all__ = [
    "chat", "Message", "ToolCall",
    "list_backends", "describe_backend", "list_ollama_models",
    "to_anthropic_tools", "to_openai_tools",
]
