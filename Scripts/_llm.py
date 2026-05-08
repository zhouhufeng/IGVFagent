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
from typing import Any, Iterable, Optional

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
    "anthropic": "claude-sonnet-4-5",
    "openai":    "gpt-4o-mini",
    "codex":     "gpt-5-codex",
    "ollama":    "qwen3:8b",
    "vllm":      "qwen3:8b",
    "tgi":       "qwen3:8b",
    "groq":      "llama-3.1-70b-versatile",
    "together":  "meta-llama/Meta-Llama-3.1-70B-Instruct-Turbo",
    "deepinfra": "meta-llama/Meta-Llama-3.1-70B-Instruct",
}


def _infer_backend(model: Optional[str]) -> str:
    if not model:
        return os.environ.get("IGVF_LLM_BACKEND", "ollama")
    m = model.lower()
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

def _chat_anthropic(messages, *, model, tools, max_tokens, temperature,
                     stop, **kwargs) -> Message:
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
    client = anthropic.Anthropic(api_key=api_key)
    system, msgs = _to_anthropic_messages(messages)
    payload: dict = {
        "model":       model,
        "messages":    msgs,
        "max_tokens":  max_tokens,
        "temperature": temperature,
    }
    if system:
        payload["system"] = system
    if tools:
        payload["tools"] = to_anthropic_tools(tools)
    if stop:
        payload["stop_sequences"] = stop
    payload.update({k: v for k, v in kwargs.items() if v is not None})

    resp = client.messages.create(**payload)
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
                "stop_sequence":    "stop_sequence"}
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
    payload: dict = {
        "model":       model,
        "messages":    _to_openai_messages(messages),
        "max_tokens":  max_tokens,
        "temperature": temperature,
    }
    if tools:
        payload["tools"] = to_openai_tools(tools)
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
_BACKENDS = {
    "anthropic": None,                                     # special-cased
    "openai":    {"base_url": "https://api.openai.com/v1",
                   "api_key_env": "OPENAI_API_KEY"},
    "codex":     {"base_url": os.environ.get("IGVF_CODEX_BASE_URL",
                                                "https://api.openai.com/v1"),
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
    "groq":      {"base_url": "https://api.groq.com/openai/v1",
                   "api_key_env": "GROQ_API_KEY"},
    "together":  {"base_url": "https://api.together.xyz/v1",
                   "api_key_env": "TOGETHER_API_KEY"},
    "deepinfra": {"base_url": "https://api.deepinfra.com/v1/openai",
                   "api_key_env": "DEEPINFRA_API_KEY"},
    "huggingface": {
        "base_url": os.environ.get(
            "IGVF_LLM_BASE_URL",
            "https://api-inference.huggingface.co/v1"),
        "api_key_env": "HF_TOKEN",
    },
    # User-provided OpenAI-compatible endpoint
    "custom":    {"base_url": os.environ.get("IGVF_LLM_BASE_URL", ""),
                   "api_key_env": os.environ.get("IGVF_LLM_API_KEY_ENV",
                                                   "IGVF_LLM_API_KEY")},
}


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
    **kwargs,
) -> Message:
    """Backend-neutral chat completion.

    Returns a :class:`Message` with normalized ``content``, ``tool_calls``,
    and ``stop_reason``. Raises ``RuntimeError`` if the resolved backend's
    SDK is missing or the credentials env var is unset.
    """
    bk = _resolve_backend(backend, model)
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
                                temperature=temperature, stop=stop, **kwargs)
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
    return ["anthropic"] + sorted(set(_BACKENDS) - {"anthropic"})


def describe_backend(name: str) -> dict:
    name = name.lower()
    if name == "anthropic":
        return {"name": name, "sdk": "anthropic",
                "api_key_env": "ANTHROPIC_API_KEY"}
    cfg = _BACKENDS.get(name)
    if not cfg:
        return {"name": name, "_unknown": True}
    return {"name": name, "sdk": "openai", **cfg}


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


__all__ = [
    "chat", "Message", "ToolCall",
    "list_backends", "describe_backend", "list_ollama_models",
    "to_anthropic_tools", "to_openai_tools",
]
