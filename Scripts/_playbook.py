"""Playbook executor — deterministic multi-step tool chains.

Use case: stop the LLM from being the source of variance for canonical
demo questions. A playbook is a YAML file that pins the exact sequence
of tool calls (and their arguments) the agent should execute. The LLM
only runs at the END to synthesize a final answer, given the captured
artefacts from each scripted step.

Output across LLM backends: **identical artefacts** (the scripted steps
are byte-for-byte reproducible) + **slightly different synthesis prose**
(the LLM still writes the prose, but with all the evidence already in
hand). For the demo, this is what you want — Plot artefacts pin
themselves; the wording can vary safely.

Playbook YAML shape
-------------------
::

    study: "Comprehensive APOE evidence pack"
    description: "Pull KG + MPRA + literature for the user's gene of interest."

    parameters:           # optional; defaults shown
      gene:    { default: APOE }
      disease: { default: "Alzheimer disease" }

    steps:
      - tool: kg_gene
        args: { symbol: "${gene}", depth: 2 }

      - tool: mpra_pull
        args: { source: catalog, limit: 25 }

      - tool: ref_validate
        args: { gene: "${gene}", disease: "${disease}" }

    synthesis: |
      Summarize the findings as a publication-ready evidence pack.
      For each evidence source, briefly state what was found.

Parameter values are substituted via ``${name}`` interpolation in arg
strings before each step runs. Override on the command line with
``--param key=value`` (repeatable).

License: Apache-2.0. Uses stdlib + PyYAML (MIT).
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)


# ─── Helpers ────────────────────────────────────────────────────────────────

def _interpolate(value: Any, params: dict[str, Any]) -> Any:
    """Recursively replace ``${name}`` in strings with parameter values."""
    if isinstance(value, str):
        out = value
        for key, val in params.items():
            out = out.replace(f"${{{key}}}", str(val))
        return out
    if isinstance(value, dict):
        return {k: _interpolate(v, params) for k, v in value.items()}
    if isinstance(value, list):
        return [_interpolate(v, params) for v in value]
    return value


def _coerce_param_value(raw: str) -> Any:
    """Cast common CLI --param strings to int/float/bool/JSON when possible."""
    if raw.lower() in ("true", "false"):
        return raw.lower() == "true"
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        pass
    if raw.startswith(("[", "{", "\"")):
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            pass
    return raw


# ─── Public API ─────────────────────────────────────────────────────────────

def load_playbook(path: Path) -> dict:
    """Parse a YAML playbook file."""
    text = Path(path).read_text(encoding="utf-8")
    pb = yaml.safe_load(text) or {}
    if "steps" not in pb or not isinstance(pb["steps"], list):
        raise ValueError(f"playbook {path} missing top-level `steps:` list")
    return pb


def run_playbook(
    playbook_path: str | Path,
    *,
    params: "dict[str, Any] | None" = None,
    backend: str | None = None,
    model: str | None = None,
    synthesize: bool = True,
    max_tokens: int = 4096,
    temperature: float = 0.0,
    on_step: Any | None = None,
) -> dict:
    """Execute the playbook end-to-end.

    Each step's tool is invoked via ``_tools.execute(name, args)`` — a
    subprocess shim that captures stdout/stderr/exit_code and parses out
    artefact paths from the skill's stdout. The LLM is invoked once at
    the end (if ``synthesize=True``) with a summary of all step results.

    Returns a dict with the same shape as ``AgentResult`` so callers
    treat playbook runs and free-form ask runs interchangeably.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _tools import execute as tools_execute  # noqa: E402

    pb = load_playbook(Path(playbook_path))
    # Merge defaults + overrides
    merged_params: dict[str, Any] = {}
    for pname, pmeta in (pb.get("parameters") or {}).items():
        merged_params[pname] = (pmeta or {}).get("default")
    if params:
        merged_params.update(params)

    started = time.time()
    artefacts: "list[str]" = []
    step_records: "list[dict]" = []

    for idx, raw_step in enumerate(pb.get("steps") or []):
        tool_name = raw_step.get("tool")
        if not tool_name:
            raise ValueError(f"playbook step {idx} missing `tool:`")
        raw_args = raw_step.get("args") or {}
        args = _interpolate(raw_args, merged_params)
        label = raw_step.get("label") or tool_name

        t0 = time.time()
        try:
            result = tools_execute(tool_name, args)
            err = None
        except KeyError as exc:
            result = {"exit_code": -1, "stdout": "",
                      "stderr": f"unknown tool: {exc}",
                      "artifacts": {}, "timed_out": False}
            err = str(exc)
        except Exception as exc:  # pylint: disable=broad-except
            result = {"exit_code": -1, "stdout": "",
                      "stderr": f"execute failed: {exc}",
                      "artifacts": {}, "timed_out": False}
            err = str(exc)
        elapsed = time.time() - t0

        step_artefacts: "list[str]" = []
        for paths in (result.get("artifacts") or {}).values():
            step_artefacts.extend(paths)
        artefacts.extend(step_artefacts)

        rec = {
            "idx": idx, "tool": tool_name, "label": label, "args": args,
            "exit_code": result.get("exit_code"),
            "stdout_tail": (result.get("stdout") or "")[-1200:],
            "stderr_tail": (result.get("stderr") or "")[-400:],
            "artifacts": result.get("artifacts") or {},
            "n_artefacts": len(step_artefacts),
            "elapsed_sec": round(elapsed, 2),
            "error": err,
        }
        step_records.append(rec)
        if on_step:
            try:
                on_step(rec)
            except Exception:
                pass
        logger.info("playbook step %d/%d  %s  exit=%s  artefacts=%d  %.2fs",
                    idx + 1, len(pb["steps"]), tool_name,
                    rec["exit_code"], rec["n_artefacts"], elapsed)

    # ── Optional synthesis (one LLM call, with all step results in context)
    final_answer: str = ""
    synth_used = False
    if synthesize and (pb.get("synthesis") or "").strip():
        try:
            from _agent import DEFAULT_SYSTEM_PROMPT
            from _llm import chat
            sections: "list[str]" = [
                f"# Playbook executed: {pb.get('study', 'untitled')}",
                "",
                pb.get("description", ""),
                "",
                "## Step results",
                "",
            ]
            for sr in step_records:
                sections.append(f"### Step {sr['idx']}: `{sr['tool']}` (exit={sr['exit_code']})")
                sections.append(f"args: `{json.dumps(sr['args'])}`")
                if sr["artifacts"]:
                    sections.append("artefacts:")
                    for k, vs in sr["artifacts"].items():
                        for v in vs:
                            sections.append(f"  - {k}: `{v}`")
                sections.append("stdout tail:")
                sections.append("```")
                sections.append(sr["stdout_tail"] or "(empty)")
                sections.append("```")
                if sr["stderr_tail"]:
                    sections.append("stderr tail:")
                    sections.append("```")
                    sections.append(sr["stderr_tail"])
                    sections.append("```")
                sections.append("")
            sections.append("## Synthesis instruction")
            sections.append(pb["synthesis"])
            user_msg = "\n".join(sections)

            messages = [
                {"role": "system", "content": DEFAULT_SYSTEM_PROMPT},
                {"role": "user", "content": user_msg},
            ]
            resp = chat(
                messages, backend=backend, model=model,
                tools=None, max_tokens=max_tokens, temperature=temperature,
            )
            final_answer = (resp.get("content") or "").strip()
            synth_used = True
        except Exception as exc:  # pylint: disable=broad-except
            logger.warning("synthesis skipped (%s)", exc)
            final_answer = (
                f"_(Synthesis skipped: {exc.__class__.__name__}: {exc})_\n\n"
                f"All scripted steps completed; see artefact list below."
            )

    # Dedupe artefacts while preserving order
    seen: "set[str]" = set()
    dedup_artefacts = []
    for a in artefacts:
        if a not in seen:
            seen.add(a)
            dedup_artefacts.append(a)

    return {
        "study": pb.get("study"),
        "playbook_path": str(playbook_path),
        "params": merged_params,
        "steps": step_records,
        "artefacts": dedup_artefacts,
        "final_answer": final_answer,
        "synthesized": synth_used,
        "elapsed_sec": round(time.time() - started, 2),
        "backend": backend,
        "model": model,
    }


def parse_param_args(values: "list[str] | None") -> dict[str, Any]:
    """Parse repeated --param flags of the form ``key=value``."""
    out: dict[str, Any] = {}
    for v in values or []:
        if "=" not in v:
            raise ValueError(f"--param must be key=value (got {v!r})")
        k, _, raw = v.partition("=")
        out[k.strip()] = _coerce_param_value(raw.strip())
    return out


__all__ = ["load_playbook", "run_playbook", "parse_param_args"]
