# Supplementary Table X — Backend roster and default models

**Model-routing backends registered in IGVFagent (v0.1.0).** A single
tool-calling interface (`_llm.chat`) dispatches to the selected backend.
Column *Transport* gives the tool-calling scheme: the Anthropic Messages API
(tool-use content blocks); the OpenAI-compatible Chat Completions
function-calling interface (structured `tool_calls`); or command-line
passthrough (an external coding-agent client invoked as a subprocess, whose
text output is parsed for `<tool_call>` XML blocks). *Default model* is the
model used when none is specified (`_llm._DEFAULT_MODELS`; backends with no
entry fall back to `qwen3:8b`, and the two CLI backends defer to the external
client's own configured model). Sampling temperature defaults to 0 on
backends that accept it.

| Backend key | Category | Transport | Endpoint (default) | Auth | Default model |
|---|---|---|---|---|---|
| `anthropic` | Hosted API | Anthropic Messages API | api.anthropic.com (SDK) | `ANTHROPIC_API_KEY` | `claude-opus-4-8` |
| `openai` | Hosted API | OpenAI Chat Completions | `https://api.openai.com/v1` | `OPENAI_API_KEY` | `gpt-4o-mini` |
| `codex` | Hosted API | OpenAI Chat Completions | `https://api.openai.com/v1` | `OPENAI_API_KEY` | `gpt-5-codex` |
| `groq` | Hosted API | OpenAI-compatible | `https://api.groq.com/openai/v1` | `GROQ_API_KEY` | `llama-3.1-70b-versatile` |
| `together` | Hosted API | OpenAI-compatible | `https://api.together.xyz/v1` | `TOGETHER_API_KEY` | `meta-llama/Meta-Llama-3.1-70B-Instruct-Turbo` |
| `deepinfra` | Hosted API | OpenAI-compatible | `https://api.deepinfra.com/v1/openai` | `DEEPINFRA_API_KEY` | `meta-llama/Meta-Llama-3.1-70B-Instruct` |
| `huggingface` | Hosted API | OpenAI-compatible | `https://api-inference.huggingface.co/v1` | `HF_TOKEN` | `qwen3:8b` (fallback) |
| `ollama` | Self-hosted | OpenAI-compatible | `http://localhost:11434/v1` | `OLLAMA_API_KEY` | `qwen3:8b` |
| `vllm` | Self-hosted | OpenAI-compatible | `http://localhost:8000/v1` | `IGVF_LLM_API_KEY` | `qwen3:8b` |
| `tgi` | Self-hosted | OpenAI-compatible | `http://localhost:3000/v1` | `IGVF_LLM_API_KEY` | `qwen3:8b` |
| `custom` | Self-hosted | OpenAI-compatible | *(user-set `IGVF_LLM_BASE_URL`)* | `IGVF_LLM_API_KEY` | `qwen3:8b` (fallback) |
| `claude_cli` | CLI passthrough | `claude` subprocess (XML parse) | — (local client) | Claude Code login | *(client default)* |
| `codex_cli` | CLI passthrough | `codex exec` subprocess (XML parse) | — (local client) | Codex/ChatGPT login | *(client default)* |

**Notes.**
- 13 backends in three families: 1 Anthropic-native, 10 OpenAI-compatible
  (6 hosted + 4 self-hosted/local), and 2 command-line passthrough clients.
- Endpoints and auth-env names are overridable via `IGVF_LLM_BASE_URL` /
  `IGVF_LLM_API_KEY` (and per-backend env vars). Values shown are defaults.
- `huggingface` and `custom` have no entry in `_DEFAULT_MODELS`, so they take
  the global fallback (`qwen3:8b`); set an explicit `--model` for these.
- Automatic retries (exponential backoff on HTTP 429 / ≥500) are provided by
  the Anthropic and OpenAI SDK defaults (`max_retries=2`) for the API and
  OpenAI-compatible backends; the CLI passthrough backends have no retry layer
  (the `codex_cli` backend does retry once with an alternate argv on a CLI-flag
  parse error).
- CLI availability is detected at runtime via `shutil.which` + a `--version`
  probe (`claude_cli_available`, `codex_cli_available`).

_Source: `Scripts/_llm.py` (`_DEFAULT_MODELS`, `BACKENDS`) at v0.1.0; roster
reproduced by `igvfagent backends`._
