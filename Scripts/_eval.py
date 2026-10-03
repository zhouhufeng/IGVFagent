"""Backend-comparison eval harness.

Runs the same question through multiple (backend, model) pairs and
measures **how much they diverge** — both in the artefacts produced and
in the wording of the final answer.

Why this exists: at temperature=0 the same model + same query gives
the same answer, but different models (Claude vs Qwen vs Gemma)
produce wildly different chains of tool calls, different artefact
sets, and different prose. This harness makes that divergence
visible and quantitative so users can decide which questions are
"any-backend safe" and which require Claude.

Metrics
-------
- ``jaccard(artefacts)``  per pair, in [0, 1]. 1.0 = identical artefact
  set; 0.0 = no overlap. The single most important number, because the
  artefacts are what users actually consume.
- ``cos_final_answer``    TF-IDF cosine over the synthesized prose, in
  [0, 1]. High = agents wrote the same answer; low = they disagreed.
  Uses only the stdlib (no scikit-learn dep).
- ``iter_diff``           |Δ iterations|. Captures planning depth diff.
- ``tool_diff``           |Δ tool_calls_made|.

License: Apache-2.0. Uses only stdlib + the existing _agent.run loop.
"""

from __future__ import annotations

import math
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any


def _tokenize(text: str) -> "list[str]":
    """Lower-case ASCII-word tokenizer."""
    return re.findall(r"[a-z0-9]{2,}", (text or "").lower())


def tfidf_cosine(text_a: str, text_b: str) -> float:
    """Two-doc TF-IDF cosine similarity (stdlib only).

    Returns 0.0 when either side is empty.
    """
    a = Counter(_tokenize(text_a))
    b = Counter(_tokenize(text_b))
    if not a or not b:
        return 0.0
    vocab = set(a) | set(b)
    # Naive IDF: log(2/df) where df ∈ {1, 2}
    idf = {w: math.log(2.0 / ((w in a) + (w in b))) + 1.0 for w in vocab}
    va = {w: a.get(w, 0) * idf[w] for w in vocab}
    vb = {w: b.get(w, 0) * idf[w] for w in vocab}
    num = sum(va[w] * vb[w] for w in vocab)
    da = math.sqrt(sum(v * v for v in va.values()))
    db = math.sqrt(sum(v * v for v in vb.values()))
    if da == 0 or db == 0:
        return 0.0
    return float(num / (da * db))


def jaccard(a: "list[str] | set[str]", b: "list[str] | set[str]") -> float:
    sa = set(a or [])
    sb = set(b or [])
    if not sa and not sb:
        return 1.0
    return float(len(sa & sb) / len(sa | sb))


def run_eval(
    question: str,
    *,
    pairs: "list[tuple[str, str | None]]",
    runs: int = 1,
    max_iterations: int = 8,
    max_tokens: int = 4096,
    temperature: float = 0.0,
    quiet: bool = True,
) -> "list[dict]":
    """Run ``runs`` replicates of ``question`` for each (backend, model)
    pair. Returns a flat list of run records."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _agent import run as agent_run  # noqa: E402

    records: "list[dict]" = []
    for backend, model in pairs:
        for run_idx in range(runs):
            label = f"{backend}/{model or 'auto'}#{run_idx}"
            t0 = time.time()
            error: str | None = None
            try:
                result = agent_run(
                    question,
                    backend=backend,
                    model=model,
                    max_iterations=max_iterations,
                    max_tokens=max_tokens,
                    temperature=temperature,
                    callback=None if quiet else __import__("_agent")._print_callback,
                    persist=False,
                )
                rec = {
                    "label": label, "backend": backend, "model": model,
                    "run_idx": run_idx,
                    "final_answer": result.final_answer or "",
                    "artefacts": list(result.artefacts or []),
                    "iterations": result.iterations,
                    "tool_calls_made": result.tool_calls_made,
                    "stop_reason": result.stop_reason,
                    "elapsed_sec": round(time.time() - t0, 2),
                    "error": None,
                }
            except Exception as exc:  # pylint: disable=broad-except
                error = f"{exc.__class__.__name__}: {exc}"
                rec = {
                    "label": label, "backend": backend, "model": model,
                    "run_idx": run_idx, "final_answer": "", "artefacts": [],
                    "iterations": 0, "tool_calls_made": 0,
                    "stop_reason": "error",
                    "elapsed_sec": round(time.time() - t0, 2),
                    "error": error,
                }
            records.append(rec)
    return records


def diff_matrix(records: "list[dict]") -> "list[dict]":
    diffs: "list[dict]" = []
    n = len(records)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = records[i], records[j]
            diffs.append({
                "a": a["label"], "b": b["label"],
                "jaccard_artefacts": jaccard(a["artefacts"], b["artefacts"]),
                "cos_final_answer": tfidf_cosine(a["final_answer"], b["final_answer"]),
                "iter_diff": abs(a["iterations"] - b["iterations"]),
                "tool_diff": abs(a["tool_calls_made"] - b["tool_calls_made"]),
            })
    return diffs


def write_eval_report(
    question: str,
    records: "list[dict]",
    diffs: "list[dict]",
    out_path: Path,
) -> Path:
    """Write a Markdown eval report."""
    lines: "list[str]" = [
        "# Backend-comparison eval report",
        "",
        f"**Question:** {question}",
        "",
        f"**Runs:** {len(records)}  ({len({r['backend'] for r in records})} backends × "
        f"{max(1, len(records) // max(1, len({r['backend'] for r in records})))} replicates each)",
        "",
        "## Per-run summary",
        "",
        "| label | iters | tools | stop | elapsed (s) | artefacts | error |",
        "|---|---:|---:|---|---:|---:|---|",
    ]
    for r in records:
        lines.append(
            f"| `{r['label']}` | {r['iterations']} | {r['tool_calls_made']} | "
            f"`{r['stop_reason']}` | {r['elapsed_sec']} | "
            f"{len(r['artefacts'])} | {r['error'] or ''} |"
        )
    lines += [
        "",
        "## Pairwise divergence",
        "",
        "Higher Jaccard / cosine = more agreement. "
        "Aim for `jaccard_artefacts ≥ 0.7` across backends for "
        "reproducibility-critical questions.",
        "",
        "| A | B | Jaccard(artefacts) | cos(final_answer) | \\|Δiters\\| | \\|Δtools\\| |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for d in diffs:
        lines.append(
            f"| `{d['a']}` | `{d['b']}` | {d['jaccard_artefacts']:.3f} | "
            f"{d['cos_final_answer']:.3f} | {d['iter_diff']} | {d['tool_diff']} |"
        )

    # Per-run artefact lists
    lines += ["", "## Artefacts by run", ""]
    for r in records:
        lines.append(f"### `{r['label']}`")
        if r["artefacts"]:
            for a in r["artefacts"]:
                lines.append(f"- `{a}`")
        else:
            lines.append("_(no artefacts)_")
        lines.append("")

    # Final answers side-by-side (truncated)
    lines += ["## Final answers (first 800 chars)", ""]
    for r in records:
        lines.append(f"### `{r['label']}`")
        ans = (r["final_answer"] or "").strip()
        lines.append("```")
        lines.append(ans[:800] + ("…" if len(ans) > 800 else ""))
        lines.append("```")
        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


__all__ = ["run_eval", "diff_matrix", "write_eval_report",
            "jaccard", "tfidf_cosine"]
