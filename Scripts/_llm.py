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
    "anthropic":  "claude-sonnet-4-5",
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
    framework = _CLAUDE_CLI_TOOL_PROMPT.format(tools_block=tools_block)
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
    cmd = ["claude", "--print", "--output-format", "text"]
    if model and model.strip():
        cmd.extend(["--model", model.strip()])
    cmd.append(prompt)

    timeout = float(os.environ.get("IGVF_LLM_TIMEOUT", "600"))
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
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

    # `codex exec` runs in non-interactive mode. Approval mode `never`
    # auto-approves any tool calls Codex's own runtime might want to
    # make — we don't expect any since our prompt explicitly tells it
    # to emit XML rather than use its built-in tools, but the flag
    # keeps the run from blocking on a TTY prompt.
    cmd = ["codex", "exec", "--ask-for-approval", "never"]
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
ANTHROPIC_MODELS = [
    "claude-opus-4-5",
    "claude-sonnet-4-5",
    "claude-sonnet-4",
    "claude-haiku-4-5",
    "claude-3-7-sonnet-latest",
    "claude-3-5-sonnet-latest",
    "claude-3-5-haiku-latest",
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
