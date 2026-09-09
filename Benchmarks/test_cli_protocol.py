#!/usr/bin/env python3
"""Does a model that ignores the tool protocol get caught, or reported done?

IGVFDS5997IVEM: the router ran two tools, then the model emitted neither a
<tool_call> nor a <final_answer>. It wrote prose claiming every tool now
returned "No such tool available" -- an inference of its own, not a real
fault. The parser's fallback returned that prose as the answer, and the run
was reported "Done - 1 iters, 2 tool calls - stop complete" having performed
none of the requested analysis.

The lapse is intermittent and could not be reproduced on demand (the same
model, prompt and container emitted valid XML on retest, with native tools
both enabled and disabled). So the fix is not to prevent it but to refuse to
pass it off as a finished answer. These cases stub the subprocess, so no
model is called.
"""
import subprocess as _real_subprocess
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import _llm  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:54} {detail}")
    if not ok:
        FAILURES.append(name)


class FakeProc:
    def __init__(self, out): self.returncode, self.stdout, self.stderr = 0, out, ""


def run_with(outputs):
    """Drive _chat_claude_cli with scripted subprocess stdout."""
    seen = []

    def fake_run(cmd, **kw):
        seen.append(cmd)
        return FakeProc(outputs[min(len(seen) - 1, len(outputs) - 1)])

    real_sub = _real_subprocess
    stub = types.ModuleType("subprocess")
    stub.run = fake_run
    stub.TimeoutExpired = real_sub.TimeoutExpired
    stub.PIPE = real_sub.PIPE
    sys.modules["subprocess"] = stub
    real_which = _llm.shutil.which if hasattr(_llm, "shutil") else None
    import shutil as _sh
    orig_which = _sh.which
    _sh.which = lambda n: "/usr/bin/claude" if n == "claude" else orig_which(n)
    try:
        msg = _llm._chat_claude_cli(
            [{"role": "user", "content": "analyse IGVFDS5997IVEM"}],
            model="claude-sonnet-5", tools=[{"name": "t", "description": "d",
                                            "parameters": {}}],
            max_tokens=1000, temperature=0, stop=None)
        return msg, seen
    finally:
        sys.modules["subprocess"] = real_sub
        _sh.which = orig_which


TOOLCALL = ("<tool_call>\n  <name>crispr_screen_analyze</name>\n"
            "  <arguments>{\"accession\": \"IGVFDS5997IVEM\"}</arguments>\n"
            "</tool_call>")
FINAL = "<final_answer>\nAll done.\n</final_answer>"
# The actual prose from the failed run.
PROSE = ("I hit a tool-availability problem partway through this task and "
         "can't currently proceed further via tool calls - every tool I try "
         "now returns \"No such tool available\", so this looks like a "
         "transient backend issue.")

# 1. A well-behaved turn is untouched, and costs exactly one subprocess call.
msg, seen = run_with([TOOLCALL])
check("tool_call turn parsed", len(msg.tool_calls) == 1, msg.stop_reason)
check("tool_call turn does not retry", len(seen) == 1, f"{len(seen)} call(s)")
msg, seen = run_with([FINAL])
check("final_answer turn parsed", msg.content.strip() == "All done.")
check("final_answer turn does not retry", len(seen) == 1, f"{len(seen)} call(s)")

# 2. Native tools are off by default, on when explicitly asked for.
check("subprocess gets --tools ''", "--tools" in seen[0], " ".join(seen[0][:8]))
import os
os.environ["IGVF_CLI_NATIVE_TOOLS"] = "1"
_msg, seen_native = run_with([FINAL])
check("IGVF_CLI_NATIVE_TOOLS=1 restores native tools",
      "--tools" not in seen_native[0])
os.environ.pop("IGVF_CLI_NATIVE_TOOLS")

# 3. Prose, then a correct retry -> the retry's result is used.
msg, seen = run_with([PROSE, TOOLCALL])
check("protocol lapse triggers exactly one retry", len(seen) == 2,
      f"{len(seen)} call(s)")
check("recovered retry yields the tool call", len(msg.tool_calls) == 1,
      msg.stop_reason)

# 4. Prose twice -> must NOT be reported as a finished answer.
msg, seen = run_with([PROSE, PROSE])
check("two lapses -> stop_reason protocol_violation",
      msg.stop_reason == "protocol_violation", msg.stop_reason)
check("two lapses -> no tool calls invented", msg.tool_calls == [])
check("answer is labelled as not complete",
      "did not complete" in msg.content.lower())
check("the model's outage claim is contradicted",
      "not a real backend fault" in msg.content.lower())
check("the prose is still shown for context", "No such tool" in msg.content)

print(f"\n13 cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
