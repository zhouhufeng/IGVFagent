"""The agent's failure banner must name the error, not the traceback header.

On 2026-09-21 the hosted agent reported `biosample_portal_census exited 1 —
Traceback (most recent call last):` and `ext_author_skill exited 3 — user tool
base_editing_screen_analyze (...) shadows an existing tool; skipped`. Neither
line was the failure: the first was the traceback's header, the second an
import-time log record; the real messages were the exception on the last
stderr line and a REFUSING line on stdout.

    python3 Scripts/test_failure_detail.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _agent  # noqa: E402

FAILED: "list[str]" = []


def check(cond: bool, msg: str) -> None:
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        FAILED.append(msg)


def main() -> int:
    tb = ("user tool `base_editing_screen_analyze` (/workspace/x.json) shadows an existing tool; skipped\n"
          "Traceback (most recent call last):\n"
          '  File "/workspace/Data/UserExtensions/skills/biosample_portal_census.py", line 123, in facets_to_df\n'
          "    for t in f[\"terms\"]:\n"
          "TypeError: string indices must be integers\n")
    d = _agent._failure_detail({"stderr": tb, "stdout": "Report: x"})
    check(d == "TypeError: string indices must be integers", f"traceback -> exception line ({d!r})")

    d = _agent._failure_detail({
        "stderr": "user tool `base_editing_screen_analyze` (/workspace/x.json) shadows an existing tool; skipped\n",
        "stdout": "REFUSING to author 'biosample_portal_census': IGVFagent already has core tools that appear to cover this.\n\n  explain_dataset\n"})
    check(d.startswith("REFUSING to author 'biosample_portal_census'"), f"noise-only stderr -> stdout refusal ({d!r})")

    d = _agent._failure_detail({"stderr": "error: --bed is required\n", "stdout": ""})
    check(d == "error: --bed is required", "plain stderr error kept")

    d = _agent._failure_detail({"stderr": "WARNING:igvfagent._tools:something\nusage: igvfagent x\nx: error: unrecognized arguments\n",
                                "stdout": ""})
    check(d == "usage: igvfagent x", "log-record prefix skipped, first real line kept")

    d = _agent._failure_detail({"stderr": "", "stdout": "Source: encode\nHTTP status: 404\nRESOLVED: no — matched nothing\n"})
    check(d.startswith("RESOLVED: no"), "stdout diagnosis line preferred over its first line")

    d = _agent._failure_detail({"stderr": "", "stdout": ""})
    check(d == "", "empty output -> empty detail")
    check(len(_agent._failure_detail({"stderr": "x" * 500, "stdout": ""})) == 200, "capped at 200 chars")

    if FAILED:
        print(f"\n{len(FAILED)} check(s) failed")
        return 1
    print("\nall checks pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
