#!/usr/bin/env python3
"""Two failures observed on the hosted deployment, both about honesty.

A user asked the agent to analyse a base-editing screen. The run made NINE
tool calls, TWO of which failed -- crispr-bean was not installed (exit 1),
and a tool was invoked without its required arguments (exit 2, argparse
usage) -- and it reported `stop_reason=complete`. The user was shown prose
about authoring skills where an analysis had been requested, with nothing
saying that two steps had not run.

The same run also shows why: unable to see `base_editing_screen_analyze`
(built-in tools are frozen at process start, so a newly added one is
invisible until the container is recreated), the agent authored
`crispr_bean_sorting` and `crispr_bean_prepare_and_run` into the SHARED
workspace instead of using the core tool.

These cases are offline.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import _agent  # noqa: E402
import ext_author_skill as ext  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:58} {detail}")
    if not ok:
        FAILURES.append(name)


# ── a run whose tool calls failed must not report success ─────────────────

pre = _agent._failure_preamble([
    {"name": "crispr_bean_prepare_and_run", "exit_code": 1,
     "detail": "ERROR: crispr-bean not found in PATH"},
    {"name": "crispr_bean_sorting", "exit_code": 2,
     "detail": "usage: igvfagent crispr-bean-sorting [-h] --guide-counts ..."},
])
check("the preamble states the answer is incomplete",
      "incomplete" in pre.lower())
check("it counts the failures", "2 tool calls" in pre, pre.splitlines()[0][:60])
check("it names each failed tool", "crispr_bean_sorting" in pre
      and "crispr_bean_prepare_and_run" in pre)
check("it gives each exit code", "exited 1" in pre and "exited 2" in pre)
check("it carries the reason", "crispr-bean not found" in pre)
check("it tells the reader how to read what follows",
      "did run" in pre.lower())

one = _agent._failure_preamble([{"name": "t", "exit_code": 3, "detail": ""}])
check("singular is not written as plural", "1 tool call did" in one,
      one.splitlines()[0][:50])

many = _agent._failure_preamble([{"name": f"t{i}", "exit_code": 1,
                                   "detail": ""} for i in range(9)])
check("a long list is truncated with a count",
      "and 3 more" in many, [l for l in many.splitlines() if "more" in l])

# A deliberate refusal is a nonzero exit and MUST surface: "this is a
# base-editing library, use bean" is exit 3, and burying it is the whole
# problem this is here to prevent.
ref = _agent._failure_preamble([{"name": "crispr_screen_analyze",
                                  "exit_code": 3,
                                  "detail": "REFUSING: base editor (ABE)"}])
check("a deliberate refusal is surfaced, not hidden",
      "REFUSING" in ref and "exited 3" in ref)

check("AgentResult carries the tally",
      "tool_calls_failed" in _agent.AgentResult.__dataclass_fields__)
check("AgentResult carries the failure detail",
      "failed_calls" in _agent.AgentResult.__dataclass_fields__)
check("the tally defaults to zero for callers that do not set it",
      _agent.AgentResult(
          final_answer="", iterations=1, tool_calls_made=0,
          stop_reason="complete", transcript=[], artefacts=[],
          transcript_path="", report_path="", backend="b",
          model="m").tool_calls_failed == 0)


# ── authoring must not shadow a core tool ─────────────────────────────────

# The two skills the agent actually wrote.
for nm, desc in (
    ("crispr_bean_sorting",
     "Run crispr-bean (BEAN) Bayesian variant-effect analysis on FACS "
     "sorting screens"),
    ("crispr_bean_prepare_and_run",
     "Prepare inputs and run bean for a base editing screen"),
):
    r = ext.duplication_refusal(nm, desc)
    check(f"authoring {nm} is refused", r is not None)
    if r:
        check("  the refusal names a core tool to use instead",
              "base_editing_screen" in r, r.splitlines()[2][:60])
        check("  it explains a stale registry, not a missing capability",
              "stale" in r and "recreate the container" in r)
        check("  it offers --force", "--force" in r)

check("a lettered-bin duplicate is refused",
      ext.duplication_refusal(
          "gradient_bin_helper",
          "Score constructs across lettered expression bins of a sorted "
          "screen") is not None)

# The false-positive direction matters as much: a guard that blocks
# legitimate utilities gets forced past and stops being read. Each of these
# matched a core tool on exactly ONE incidental token before the rule
# required two.
for nm, desc in (
    ("tar_extract", "Extract a .tar.gz archive to a directory"),
    ("gz_head", "Print the first lines of a gzipped text file"),
    ("write_text_file", "Write a string to a path"),
    ("fetch_ukbb_gwas", "Download UK Biobank GWAS summary statistics"),
    ("run_ppi_predict", "Predict protein-protein interactions"),
):
    check(f"authoring {nm} is allowed",
          ext.duplication_refusal(nm, desc) is None)

check("--force authors despite an overlap",
      ext.duplication_refusal("crispr_bean_sorting", "run bean",
                               force=True) is None)
check("two shared name tokens are required",
      len(ext.find_similar_core_tools("gz_head", "print lines")) == 0)
check("a genuine duplicate scores above the floor",
      (ext.find_similar_core_tools(
          "crispr_bean_sorting", "run crispr-bean on a sorting screen")
       or [{"score": 0}])[0]["score"] >= 4.0)
check("user extensions are not matched against, only core tools",
      all("crispr_bean" not in h["name"]
          for h in ext.find_similar_core_tools(
              "crispr_bean_sorting", "run bean")))

print(f"\n{29} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
