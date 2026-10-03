"""Offline tests for durable agent jobs and their chat routing.

    python3 Scripts/test_agent_jobs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_jobs as aj  # noqa: E402
import jobs_ui as ju  # noqa: E402


def main() -> int:
    rc = aj.main(["selftest"])
    fails = []

    def check(name, cond):
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
        if not cond:
            fails.append(name)

    check("reproduce requests route to a job", ju.route("can you reproduce this paper https://doi.org/x"))
    check("pipelines from raw reads route to a job", ju.route("process IGVFDS1 from raw reads end-to-end"))
    check("a lookup stays a single reply", not ju.route("what is APOE?"))
    check("mode always / never override", ju.route("what is APOE?", "always") and
          not ju.route("reproduce the paper", "never"))
    check("'yes please continue with your planned steps' resumes", ju.is_continue(
        "yes please continue with your planed reproduce steps"))
    check("'continue' inside a new question does not", not ju.is_continue("how do I continue a stopped job?"))
    q = ju.job_query("reproduce it", [{"role": "user", "content": "paper at Data/Uploads/x.pdf"},
                                      {"role": "assistant", "content": "It is Wang 2026"}])
    check("the job carries the conversation it refers to", "Data/Uploads/x.pdf" in q and "Wang 2026" in q)
    ok = rc == 0 and not fails
    print("all checks pass" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
