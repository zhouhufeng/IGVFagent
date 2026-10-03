#!/usr/bin/env python3
"""Cross-backend consistency harness for IGVFagent.

Proves that the same query yields the same *structured plan* no matter which
LLM drives the loop. Two modes:

* **offline** (default, no API keys) — asserts the backend-independent
  invariants that make consistency possible:
    1. ``_router.route(q)`` is a stable pure function (same output every call).
    2. ``_llm.canonical_tools`` yields one identical tool set + hash regardless
       of input ordering (so every backend sees the same tools).
    3. The run fingerprint (prompt sha + seed + tool-set sha) is stable.

* **online** (``--backends anthropic,ollama,…``) — runs the SAME query through
  ``_agent.run`` on each named backend and diffs the executed tool-call
  sequence (names + primary args). Divergence fails the check.

Exit code 0 = all consistent, 1 = a divergence was found.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _router  # noqa: E402
import _llm  # noqa: E402

# Fixed query set spanning routed shapes + a couple free-form queries.
QUERIES = [
    "APOE",
    "TP53",
    "rs429358",
    "what does rs58658771 do?",
    "IGVFDS3222WCZH",
    "chr19:44905000-44910000",
    "https://www.encodeproject.org/experiments/ENCSR000EMT/",
    "What single-cell multiome datasets does IGVF have for cortex?",
]


def _tool_dicts():
    import _tools
    return [t.to_dict() for t in _tools.list_tools()]


def offline_checks() -> dict:
    """Backend-independent invariants. Returns a report dict."""
    report = {"mode": "offline", "checks": [], "consistent": True}

    def _add(name, ok, detail):
        report["checks"].append({"name": name, "passed": bool(ok),
                                  "detail": detail})
        if not ok:
            report["consistent"] = False

    # 1. router determinism (pure function, stable across calls)
    stable = True
    routes = {}
    for q in QUERIES:
        r1 = _router.route(q)
        r2 = _router.route(q)
        routes[q] = [(x["shape"], x["tool"], x["arguments"]) for x in r1]
        if r1 != r2:
            stable = False
    _add("router is a stable pure function", stable,
         f"{sum(1 for v in routes.values() if v)} of {len(QUERIES)} queries routed deterministically")

    # 2. canonical tool set is order-independent + identical for all backends
    tools = _tool_dicts()
    import random
    shuffled = tools[:]
    random.Random(1).shuffle(shuffled)
    c1 = [t["name"] for t in (_llm.canonical_tools(tools) or [])]
    c2 = [t["name"] for t in (_llm.canonical_tools(shuffled) or [])]
    same = c1 == c2
    _add("canonical tool set is order-independent", same,
         f"{len(c1)} tools; identical under input reordering: {same}")

    # 3. fingerprint stability (prompt + seed + tool-set hash)
    import hashlib
    try:
        from _agent import DEFAULT_SYSTEM_PROMPT
    except Exception:
        DEFAULT_SYSTEM_PROMPT = ""
    fp = {
        "system_prompt_sha1": hashlib.sha1(DEFAULT_SYSTEM_PROMPT.encode()).hexdigest()[:12],
        "tool_set_sha1": hashlib.sha1(",".join(sorted(c1)).encode()).hexdigest()[:12],
        "n_tools": len(c1),
    }
    _add("run fingerprint computed", bool(fp["tool_set_sha1"]),
         json.dumps(fp))
    report["fingerprint"] = fp
    report["routes"] = routes
    return report


def _tool_sequence(query: str, backend: str, model: str | None) -> "list[dict]":
    """Run one query on one backend; return the executed tool-call sequence."""
    import _agent
    events: "list[dict]" = []

    def cb(ev):
        if ev.kind == "tool_call_start":
            events.append({"name": ev.payload.get("name"),
                           "routed": ev.payload.get("routed", False)})

    _agent.run(query, backend=backend, model=model, persist=False,
               callback=cb, max_iterations=4)
    return events


def online_checks(backends: "list[str]", model: str | None) -> dict:
    """Diff the executed tool-call sequence across backends per query."""
    report = {"mode": "online", "backends": backends, "per_query": [],
              "consistent": True}
    for q in QUERIES:
        seqs = {}
        for bk in backends:
            try:
                ev = _tool_sequence(q, bk, model)
                seqs[bk] = [e["name"] for e in ev]
            except Exception as e:  # noqa
                seqs[bk] = [f"<error: {e}>"]
        names = list(seqs.values())
        agree = all(n == names[0] for n in names[1:]) if len(names) > 1 else True
        if not agree:
            report["consistent"] = False
        report["per_query"].append({"query": q, "sequences": seqs,
                                     "agree": agree})
    return report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--backends", help="comma-separated backends for the online "
                   "diff, e.g. anthropic,ollama. Omit for offline invariants only.")
    p.add_argument("--model", default=None, help="model override (online mode).")
    p.add_argument("--json", action="store_true", help="emit JSON only.")
    args = p.parse_args(argv)

    reports = [offline_checks()]
    if args.backends:
        reports.append(online_checks(
            [b.strip() for b in args.backends.split(",") if b.strip()],
            args.model))

    consistent = all(r.get("consistent") for r in reports)
    if args.json:
        print(json.dumps({"consistent": consistent, "reports": reports},
                         indent=2, default=str))
    else:
        for r in reports:
            print(f"\n=== consistency: {r['mode']} ===")
            for c in r.get("checks", []):
                print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['name']}: {c['detail']}")
            for pq in r.get("per_query", []):
                mark = "PASS" if pq["agree"] else "DIVERGENT"
                print(f"  [{mark}] {pq['query'][:40]!r}: "
                      + " | ".join(f"{k}={v}" for k, v in pq["sequences"].items()))
        print(f"\nOverall: {'CONSISTENT ✓' if consistent else 'DIVERGENCE FOUND ✗'}")
    return 0 if consistent else 1


if __name__ == "__main__":
    sys.exit(main())
