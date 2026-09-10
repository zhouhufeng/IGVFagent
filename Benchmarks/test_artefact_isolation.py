#!/usr/bin/env python3
"""Are the artefacts shown the ones this run produced?

An external evaluation (Docs/BugZilla/
IGVFagent_MultiFlow_Public_Data_Test_Report_2026-09-09.md, section 7.6) saw
about 70 entries in the Artefacts panel during a six-gene kidney query,
including reports for APOE, TP53, BRCA1 and BRCA2, with several paths listed
more than once.

Two causes, both reproduced. Docs/KGTraversal holds one directory per past
question, so expanding that parent yields 30 files from unrelated runs --
exactly the APOE/TP53/BRCA1 set observed. And every path mentioned anywhere
in the answer was added verbatim, so the same file spelled "Docs/x/y.png",
"/workspace/Docs/x/y.png" and as an absolute path appeared three times.

The rule these cases pin: a run owns what its own tools reported, plus what
its answer cites from inside those directories. Nothing else, however
plausibly named.
"""
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))

# streamlit_app imports streamlit at module scope; stub what these functions
# touch so the logic can be tested without a Streamlit runtime.
import importlib.util  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:58} {detail}")
    if not ok:
        FAILURES.append(name)


import streamlit_app as sa  # noqa: E402

ROOT = sa._PROJECT_ROOT


def rel(p):
    return str(Path(p))


# ── run-directory recognition ─────────────────────────────────────────────

check("timestamped run dir recognised",
      sa._looks_like_run_dir("20260909_122712_gene_WT1"))
check("dash-separated stamp recognised",
      sa._looks_like_run_dir("20260909-122712-gene_WT1"))
check("Plots/ is NOT a run dir", not sa._looks_like_run_dir("Plots"))
check("Manifests/ is NOT a run dir", not sa._looks_like_run_dir("Manifests"))
check("a plain name is NOT a run dir", not sa._looks_like_run_dir("KGTraversal"))


# ── path normalisation kills the duplicates ───────────────────────────────

import tempfile  # noqa: E402

tmp = Path(tempfile.mkdtemp(prefix="artefact_test_"))
run = tmp / "20260909_120000_gene_WT1"
(run / "Plots").mkdir(parents=True)
(run / "report.md").write_text("r")
(run / "Plots" / "fig.png").write_bytes(b"\x89PNG")
sibling = tmp / "20260820_091629_gene_APOE"
sibling.mkdir()
(sibling / "report.md").write_text("other run")

# Point the module at this fixture tree.
sa._PROJECT_ROOT = tmp
# The fixture lives outside the real project root, which _pathguard
# correctly refuses. Stub it so these cases test artefact SCOPING and not
# the guard, which Benchmarks/test_guide_key_selection.py covers separately.
sa._pathguard = types.SimpleNamespace(
    is_safe_artifact=lambda p, **k: True,
    filter_artifacts=lambda ps: ps,
    why_blocked=lambda p, **k: None,
)


reported = [str(run / "report.md")]
answer = (f"See {run/'report.md'} and {run/'Plots'/'fig.png'}.\n"
          f"Also /workspace/{(run/'report.md').relative_to(tmp)} again.\n")
got = sa._collect_run_artefacts(reported, answer)
check("the same file in three spellings appears once",
      sum(1 for g in got if g.endswith("report.md")) == 1, str(got))
check("a figure cited from inside the run is kept",
      any(g.endswith("fig.png") for g in got), str(got))

# The heart of it: another run's report, cited in the answer, must not appear.
answer2 = f"Compare with the earlier {sibling/'report.md'} findings."
got2 = sa._collect_run_artefacts(reported, answer2)
check("another run's file cited in the answer is EXCLUDED",
      not any("APOE" in g for g in got2), str(got2))
check("this run's own file still present",
      any(g.endswith("report.md") for g in got2), str(got2))

# No reported artefacts -> no scope -> nothing adopted from the answer.
got3 = sa._collect_run_artefacts([], answer)
check("with nothing reported, nothing is adopted from the answer",
      got3 == [], str(got3))


# ── directory expansion must not cross into sibling runs ──────────────────

# Expanding the RUN directory: its own files and its Plots/ subdir.
exp = sa._expand_artefact_dirs([str(run)])
check("expanding a run dir finds its own files",
      any(e.endswith("report.md") for e in exp), str(exp))
check("expanding a run dir descends into Plots/",
      any(e.endswith("fig.png") for e in exp), str(exp))

# Expanding the PARENT: must not descend into the run directories inside it.
exp2 = sa._expand_artefact_dirs([str(tmp)])
check("expanding a parent does NOT reach into sibling run dirs",
      not any("gene_APOE" in e for e in exp2), str(exp2)[:120])
check("expanding a parent does not reach this run's files either",
      not any(e.endswith("fig.png") for e in exp2), str(exp2)[:120])

# A non-run subdirectory (Plots) is still traversed, which is the case that
# made expansion worth having: answers name the run folder, not each figure.
exp3 = sa._expand_artefact_dirs([str(run)])
check("Plots/ is traversed because it is not a run dir",
      any("fig.png" in e for e in exp3))

import shutil  # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

print(f"\n{16} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
