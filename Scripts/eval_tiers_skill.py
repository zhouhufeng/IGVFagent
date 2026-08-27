"""Three-tier evaluation — skill correctness, planning, conclusion validity.

An agent can fail in three independent ways, and one aggregate score conceals
which occurred. This skill implements the two tiers the existing benchmark
suite does not cover.

**Tier 1 — skill correctness.** Already covered by ``Benchmarks/`` and scored
by ``Benchmarks/concordance.py``. Those cases execute fixed command sequences
with no model in the loop, so they establish that implementations are right.
``Benchmarks/taxonomy.py`` reports what each one actually asserts.

**Tier 2 — planning and tool selection** *(here)*. Holds the skills fixed and
measures the decision layer: given a question, does the agent select the
necessary tools, order them so dependencies are respected, bind the arguments
that matter, and recover when a tool fails? Scored against a *gold plan* — the
tool set a domain expert considers necessary and sufficient, with which
arguments are load-bearing and which are free.

Recall is weighted above precision: a missing essential tool invalidates the
conclusion, whereas an extra call mostly costs tokens.

**Tier 3 — conclusion validity** *(here)*. Asks the only question a biologist
cares about: is the answer right? An independent verifier model receives the
question, the final answer, and excerpts of the produced artefacts — but
**not** the reasoning trace, so it cannot be talked into agreeing by the same
chain that produced the error. It judges each claim as supported,
unsupported, or contradicted by the artefacts.

Why the existing metrics cannot do this: ``_eval.py`` computes Jaccard overlap
of artefact sets and TF-IDF cosine of prose. Both measure *agreement*. Two
runs can agree perfectly and both be wrong, and two runs can disagree
textually while reaching the same correct conclusion — "BRCA1 is associated
with Fanconi anemia" versus "BRCA1 variants cause FA complementation group S"
scores poorly on cosine and is biologically concordant.

Usage::

    igvfagent eval-tiers tier2 --case Benchmarks/tier2/kg_gene_basic.json
    igvfagent eval-tiers tier2 --all
    igvfagent eval-tiers tier2 --all --inject-failures
    igvfagent eval-tiers tier3 --case <case.json>
    igvfagent eval-tiers list

Cases are JSON; see ``Benchmarks/tier2/README.md`` for the schema.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

try:
    from igvfagent import _agent, _llm, _pathguard
except Exception:  # pragma: no cover
    import _agent  # type: ignore
    import _llm  # type: ignore
    import _pathguard  # type: ignore

__all__ = ["main", "score_plan", "run_tier2_case", "verify_conclusion"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
# Prefer the workspace copy; fall back to the checkout that ships alongside
# Scripts/. The runtime container installs only the package, not Benchmarks/,
# so on a hosted deployment neither exists — eval-tiers is a development and CI
# tool, and the error below says that rather than reporting a bare path.
_CASE_CANDIDATES = [
    ROOT / "Benchmarks" / "tier2",
    Path(__file__).resolve().parents[1] / "Benchmarks" / "tier2",
]
CASE_DIR = next((d for d in _CASE_CANDIDATES if d.is_dir()),
                _CASE_CANDIDATES[0])
OUT_DIR = ROOT / "Docs" / "EvalTiers"


# ---------------------------------------------------------------------------
# Tier 2 — planning
# ---------------------------------------------------------------------------

class _Recorder:
    """Capture the tool calls the agent chose, in order, with arguments.

    tool_call_start carries both name and arguments; the transcript records
    only the name, so the callback is the only place argument binding can be
    observed.
    """

    def __init__(self) -> None:
        self.calls: "list[dict]" = []
        self.events: "list[str]" = []

    def __call__(self, event) -> None:
        kind = getattr(event, "kind", None) or getattr(event, "type", None)
        payload = getattr(event, "payload", None) or getattr(event, "data", {}) or {}
        if kind:
            self.events.append(str(kind))
        if kind == "tool_call_start":
            self.calls.append({"name": payload.get("name"),
                               "arguments": payload.get("arguments") or {},
                               "routed": bool(payload.get("routed"))})


def score_plan(case: dict, calls: "list[dict]") -> dict:
    """Score observed tool calls against the gold plan."""
    names = [c["name"] for c in calls]
    name_set = set(names)

    required = list(case.get("required_tools") or [])
    # "any_of" groups: satisfied when at least one member was called. Used
    # where several tools are legitimate routes to the same evidence, so a
    # defensible alternative is not scored as a miss.
    any_of = [list(g) for g in (case.get("any_of") or [])]
    forbidden = set(case.get("forbidden_tools") or [])
    req_args = case.get("required_args") or {}
    ordering = [tuple(pair) for pair in (case.get("ordering") or [])]

    hit = [t for t in required if t in name_set]
    missing = [t for t in required if t not in name_set]
    groups_ok = [g for g in any_of if name_set & set(g)]
    groups_missing = [g for g in any_of if not (name_set & set(g))]

    n_required = len(required) + len(any_of)
    n_satisfied = len(hit) + len(groups_ok)
    recall = n_satisfied / n_required if n_required else 1.0

    useful = set(required) | {t for g in any_of for t in g}
    extra = [t for t in names if t not in useful]
    precision = (len(names) - len(extra)) / len(names) if names else 0.0

    # Argument binding: only the arguments the gold plan marks load-bearing,
    # and only for tools that were actually CALLED. Scoring arguments on an
    # uncalled tool double-counts: not calling it is already a recall miss,
    # and where an any_of group offers alternatives the unused alternative is
    # not a defect at all. (Observed: an agent satisfied
    # {catalog_get_entity|kg_gene} via kg_gene and was scored 0.0 on
    # catalog_get_entity's arguments — a scorer bug, not an agent failure.)
    arg_results = []
    for tool, needed in req_args.items():
        made = [c for c in calls if c["name"] == tool]
        if not made:
            continue
        for call in made:
            for arg in needed:
                got = call["arguments"].get(arg)
                arg_results.append({"tool": tool, "arg": arg,
                                    "bound": got not in (None, "", [])})
    arg_ok = sum(1 for a in arg_results if a["bound"])
    arg_score = arg_ok / len(arg_results) if arg_results else 1.0

    # Ordering: dependency respected if `before` first appears earlier than
    # `after`. A pair where either tool was never called is not violated —
    # that is a recall miss, already counted, and double-penalising it would
    # conflate two different failures.
    violations = []
    for before, after in ordering:
        if before in names and after in names:
            if names.index(before) > names.index(after):
                violations.append([before, after])

    forbidden_used = sorted(name_set & forbidden)

    return {
        "n_calls": len(names),
        "calls": names,
        "recall": round(recall, 3),
        "precision": round(precision, 3),
        "required_hit": hit,
        "required_missing": missing,
        "any_of_missing": groups_missing,
        "argument_score": round(arg_score, 3),
        "argument_detail": arg_results,
        "ordering_violations": violations,
        "forbidden_used": forbidden_used,
        "redundant_calls": extra,
        # Recall gates the verdict: a missing essential tool means the
        # conclusion cannot be supported, whatever else went right.
        "passed": (not missing and not groups_missing and not forbidden_used
                   and not violations and arg_score >= 0.999),
    }


def run_tier2_case(case: dict, *, backend=None, model=None,
                   max_iterations=12, inject_failure: "str | None" = None,
                   verbose=False) -> dict:
    """Execute one Tier-2 case and score the plan."""
    rec = _Recorder()
    started = time.time()
    env_backup = None

    if inject_failure:
        # Force a named tool to fail so recovery becomes observable. A run
        # that cannot proceed when an archive errors has a planning defect,
        # not a skill defect, and nothing in the existing suite tests it.
        env_backup = os.environ.get("IGVF_EVAL_FAIL_TOOL")
        os.environ["IGVF_EVAL_FAIL_TOOL"] = inject_failure

    try:
        result = _agent.run(
            case["question"],
            backend=backend, model=model,
            max_iterations=max_iterations,
            callback=rec if not verbose else None,
            persist=False,
        )
    except Exception as e:
        return {"case": case.get("id"), "error": f"{type(e).__name__}: {e}",
                "plan": score_plan(case, rec.calls), "elapsed_s": round(time.time()-started, 1)}
    finally:
        if inject_failure:
            if env_backup is None:
                os.environ.pop("IGVF_EVAL_FAIL_TOOL", None)
            else:
                os.environ["IGVF_EVAL_FAIL_TOOL"] = env_backup

    plan = score_plan(case, rec.calls)
    out = {
        "case": case.get("id"),
        "question": case["question"],
        "injected_failure": inject_failure,
        "stop_reason": result.stop_reason,
        "iterations": result.iterations,
        "elapsed_s": round(time.time() - started, 1),
        "plan": plan,
        "final_answer": result.final_answer,
        "artefacts": result.artefacts,
        "backend": result.backend,
        "model": result.model,
    }
    if inject_failure:
        # Recovery = still reached a conclusion despite the induced failure.
        out["recovered"] = bool(result.final_answer) and result.stop_reason == "complete"
    return out


# ---------------------------------------------------------------------------
# Tier 3 — conclusion validity
# ---------------------------------------------------------------------------

_VERIFIER_PROMPT = """\
You are verifying a scientific claim against evidence. You did NOT produce the
answer and you cannot see how it was produced — judge only what the evidence
supports.

QUESTION ASKED:
{question}

ANSWER GIVEN:
{answer}

EVIDENCE (excerpts from the files the analysis produced):
{evidence}

For each substantive factual claim in the answer, decide:
  SUPPORTED    — the evidence contains it
  UNSUPPORTED  — plausible but absent from the evidence
  CONTRADICTED — the evidence says otherwise

Judge only claims about data. Ignore hedging, suggested next steps, and file
paths. If the evidence is insufficient to judge a claim, mark it UNSUPPORTED
rather than guessing.

Reply as JSON only:
{{"claims": [{{"claim": "...", "verdict": "SUPPORTED|UNSUPPORTED|CONTRADICTED",
 "why": "..."}}],
 "n_supported": 0, "n_unsupported": 0, "n_contradicted": 0,
 "overall": "VALID|PARTIALLY_SUPPORTED|INVALID",
 "reason": "one sentence"}}
"""


def _evidence_excerpts(artefacts, *, per_file=1800, max_files=6) -> str:
    """Read the produced artefacts, contained by _pathguard."""
    out = []
    for path in (artefacts or [])[:max_files]:
        if not _pathguard.is_safe_artifact(path):
            continue
        p = Path(path)
        try:
            with open(p, "rb") as fh:
                head = fh.read(8192)
                if b"\x00" in head:
                    out.append(f"--- {p.name} --- (binary, {p.stat().st_size} bytes)")
                    continue
            text = p.read_text(encoding="utf-8", errors="replace")[:per_file]
        except OSError:
            continue
        out.append(f"--- {p.name} ---\n{text}")
    return "\n\n".join(out) or "(no readable artefacts were produced)"


def verify_conclusion(question: str, answer: str, artefacts,
                      *, backend=None, model=None) -> dict:
    """Independent verification of an answer against its own artefacts."""
    evidence = _evidence_excerpts(artefacts)
    prompt = _VERIFIER_PROMPT.format(question=question, answer=answer or "(none)",
                                      evidence=evidence)
    try:
        msg = _llm.chat([{"role": "user", "content": prompt}],
                        backend=backend, model=model,
                        max_tokens=4000, temperature=0.0)
    except Exception as e:
        return {"error": f"verifier call failed: {type(e).__name__}: {e}"}

    text = ""
    for attr in ("content", "text"):
        v = getattr(msg, attr, None)
        if isinstance(v, str) and v.strip():
            text = v
            break
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        return {"error": "verifier did not return JSON", "raw": text[:400]}
    try:
        verdict = json.loads(m.group(0))
    except json.JSONDecodeError as e:
        return {"error": f"verifier JSON invalid: {e}", "raw": text[:400]}
    verdict["evidence_files"] = len(artefacts or [])
    return verdict


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _load_cases(args) -> "list[dict]":
    global CASE_DIR
    if getattr(args, "case_dir", None):
        d = Path(args.case_dir)
        CASE_DIR = d if d.is_absolute() else ROOT / d
    if args.case:
        p = Path(args.case)
        if not p.is_absolute():
            p = ROOT / p
        return [json.loads(p.read_text())]
    if not CASE_DIR.is_dir():
        raise SystemExit(
            "error: Tier-2 cases not found. eval-tiers is a development and CI "
            "tool and needs the repository checkout — the runtime container "
            "ships only the installed package, not Benchmarks/.\n"
            "  Run it from a clone, or pass --case-dir <path>.\n"
            f"  Looked in: {', '.join(str(d) for d in _CASE_CANDIDATES)}")
    return [json.loads(p.read_text())
            for p in sorted(CASE_DIR.glob("*.json"))]


def cmd_list(args) -> int:
    cases = _load_cases(argparse.Namespace(
        case=None, case_dir=getattr(args, "case_dir", None)))
    print(f"{'CASE':<26} {'REQUIRED TOOLS':<46} QUESTION")
    for c in cases:
        req = ",".join((c.get("required_tools") or [])
                       + ["|".join(g) for g in (c.get("any_of") or [])])
        print(f"{c.get('id',''):<26} {req[:46]:<46} {c['question'][:52]}")
    print(f"\n{len(cases)} Tier-2 case(s) in {CASE_DIR.relative_to(ROOT)}")
    return 0


def cmd_tier2(args) -> int:
    cases = _load_cases(args)
    results = []
    for case in cases:
        print(f"\n=== {case.get('id')} ===")
        r = run_tier2_case(case, backend=args.backend, model=args.model,
                           max_iterations=args.max_iterations)
        results.append(r)
        _print_plan(r)
        if args.inject_failures:
            for tool in (case.get("inject_failures") or [])[:2]:
                print(f"  -- with induced failure of {tool} --")
                rf = run_tier2_case(case, backend=args.backend, model=args.model,
                                    max_iterations=args.max_iterations,
                                    inject_failure=tool)
                rf["case"] = f"{case.get('id')}+fail[{tool}]"
                results.append(rf)
                print(f"     recovered: {rf.get('recovered')}  "
                      f"stop={rf.get('stop_reason')}")
    return _write(results, "tier2", args)


def _print_plan(r) -> None:
    p = r.get("plan") or {}
    print(f"  calls        : {', '.join(p.get('calls') or []) or '(none)'}")
    print(f"  recall       : {p.get('recall')}   precision: {p.get('precision')}")
    print(f"  arguments    : {p.get('argument_score')}")
    if p.get("required_missing"):
        print(f"  MISSING      : {', '.join(p['required_missing'])}")
    if p.get("any_of_missing"):
        print(f"  MISSING any_of: {p['any_of_missing']}")
    if p.get("ordering_violations"):
        print(f"  ORDER        : {p['ordering_violations']}")
    if p.get("forbidden_used"):
        print(f"  FORBIDDEN    : {', '.join(p['forbidden_used'])}")
    print(f"  verdict      : {'PASS' if p.get('passed') else 'FAIL'}"
          f"   stop={r.get('stop_reason')}  {r.get('elapsed_s')}s")


def cmd_tier3(args) -> int:
    cases = _load_cases(args)
    results = []
    for case in cases:
        print(f"\n=== {case.get('id')} (tier 3) ===")
        run = run_tier2_case(case, backend=args.backend, model=args.model,
                             max_iterations=args.max_iterations)
        v = verify_conclusion(case["question"], run.get("final_answer"),
                              run.get("artefacts"),
                              backend=args.verifier_backend or args.backend,
                              model=args.verifier_model)
        results.append({"case": case.get("id"), "plan": run.get("plan"),
                        "verification": v})
        if v.get("error"):
            print(f"  verifier error: {v['error']}")
            continue
        print(f"  supported    : {v.get('n_supported')}")
        print(f"  unsupported  : {v.get('n_unsupported')}")
        print(f"  contradicted : {v.get('n_contradicted')}")
        print(f"  overall      : {v.get('overall')} — {v.get('reason','')[:90]}")
        for c in (v.get("claims") or []):
            if c.get("verdict") != "SUPPORTED":
                print(f"    [{c['verdict']}] {c['claim'][:88]}")
    return _write(results, "tier3", args)


def _write(results, tier, args) -> int:
    ts = time.strftime("%Y%m%d_%H%M%S")
    run = OUT_DIR / f"{ts}_{tier}"
    run.mkdir(parents=True, exist_ok=True)
    out = run / f"{tier}_results.json"
    out.write_text(json.dumps({"tier": tier, "n": len(results),
                               "results": results}, indent=2))
    print(f"\nReport:        {out}")
    if tier == "tier2":
        failed = [r for r in results if not (r.get("plan") or {}).get("passed")]
        print(f"Wrote:         {out}")
        print(f"{len(results) - len(failed)}/{len(results)} case(s) passed")
        return 1 if failed else 0
    bad = [r for r in results
           if (r.get("verification") or {}).get("overall") not in ("VALID", None)]
    print(f"{len(results) - len(bad)}/{len(results)} conclusion(s) fully valid")
    return 1 if bad else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent eval-tiers",
        description="Tier 2 (planning) and Tier 3 (conclusion validity) "
                    "evaluation. Tier 1 lives in Benchmarks/.")
    sub = p.add_subparsers(dest="cmd", required=True)
    lst = sub.add_parser("list", help="List Tier-2 cases")
    lst.add_argument("--case-dir", help="Directory of Tier-2 case JSON files")
    for name in ("tier2", "tier3"):
        s = sub.add_parser(name)
        s.add_argument("--case", help="One case JSON (default: all)")
        s.add_argument("--all", action="store_true", help="All cases (default)")
        s.add_argument("--backend")
        s.add_argument("--model")
        s.add_argument("--max-iterations", type=int, default=12)
        s.add_argument("--case-dir", help="Directory of Tier-2 case JSON files")
        if name == "tier2":
            s.add_argument("--inject-failures", action="store_true",
                           help="Also run each case with induced tool failures")
        else:
            s.add_argument("--verifier-backend",
                           help="Backend for the verifier (default: same)")
            s.add_argument("--verifier-model",
                           help="Model for the verifier. Use a DIFFERENT model "
                                "from the one under test where possible.")
    args = p.parse_args(argv)
    return {"list": cmd_list, "tier2": cmd_tier2, "tier3": cmd_tier3}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
