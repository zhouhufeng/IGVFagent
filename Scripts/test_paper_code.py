"""Offline checks for paper reproduction from the authors' code.

- paper_code_skill selftest: repository ranking from a harvest, Rmd parsing,
  sections, inventory, printed-value comparison, a real script run, report.
- agent_jobs.code_route_gap: a reproduction job whose harvested Code
  Availability names a repository cannot pass verification until it has run
  that code (or blocked the stage with a reason).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main() -> int:
    ok = subprocess.run([sys.executable, str(HERE / "paper_code_skill.py"), "selftest"]).returncode == 0
    with tempfile.TemporaryDirectory() as td:
        os.environ["IGVF_PROJECT_ROOT"] = td
        import importlib
        import agent_jobs as aj
        importlib.reload(aj)
        hv = Path(td) / "Docs" / "Benchmark" / "r" / "harvest.json"
        hv.parent.mkdir(parents=True)
        hv.write_text(json.dumps({"accessions": {"github_repo": [
            {"value": "FowlerLab/VAMPseq.", "in_data_availability": True}]}}))
        rel = str(hv.relative_to(td))
        job = {"query": "Reproduce Matreyek 2018 PTEN VAMP-seq"}
        plan = {"stages": [{"id": "harvest", "title": "harvest", "status": "done", "evidence": [rel]},
                           {"id": "data", "title": "MaveDB scores", "status": "done", "evidence": []}]}

        def check(name, cond):
            nonlocal ok
            ok &= bool(cond)
            print(f"  {'ok  ' if cond else 'FAIL'}  {name}")

        gap = aj.code_route_gap(job, plan)
        check("reproduction without running the named code is flagged", gap and "FowlerLab/VAMPseq" in gap)
        plan2 = json.loads(json.dumps(plan))
        plan2["stages"].append({"id": "run_code", "title": "run code", "status": "done",
                                "evidence": ["Docs/PaperCode/x/summary.json"]})
        check("running the code clears it", aj.code_route_gap(job, plan2) is None)
        plan3 = json.loads(json.dumps(plan))
        plan3["stages"].append({"id": "run_code", "title": "run the authors' code", "status": "blocked",
                                "notes": ["2026 repository deleted"]})
        check("a blocked code stage with a reason clears it", aj.code_route_gap(job, plan3) is None)
        check("non-reproduction jobs are not affected",
              aj.code_route_gap({"query": "list K562 datasets"}, plan) is None)
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
