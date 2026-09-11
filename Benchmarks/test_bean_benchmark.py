#!/usr/bin/env python3
"""The replication harness for Ryu et al. 2024 (crispr-bean).

A benchmark is only worth running if its own arithmetic is right and its
targets cannot drift. Two things are tested here:

  * the metrics, against cases whose answers are known by construction --
    AUPRC especially, because the paper's comparison is a rare-positive
    classification where a wrong tie rule silently inflates the area;
  * the claim table, because a replication that quotes its target from
    memory is not checkable. Every claim carries the sentence it came from.

No network and no data file: `fetch` and `describe` are exercised elsewhere.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import bean_benchmark_skill as bb  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:64} {detail}")
    if not ok:
        FAILURES.append(name)


# ── Spearman ──────────────────────────────────────────────────────────────
# Exact equality is wrong here: the sums of products land at
# 0.9999999999999998, which is correct to 1e-16 and not 1.0. Asserting == 1.0
# tests the floating-point representation, not the statistic.
def close(a, b, tol=1e-9):
    return abs(a - b) <= tol

check("perfect agreement is 1", close(bb.spearman([1, 2, 3, 4], [1, 2, 3, 4]), 1.0))
check("perfect inversion is -1", close(bb.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0))
check("it is rank-based, so a monotone nonlinearity is still 1",
      close(bb.spearman([1, 2, 3, 4], [1, 4, 9, 16]), 1.0))
check("ties are averaged, not broken arbitrarily",
      abs(bb.spearman([1, 1, 2, 2], [1, 1, 2, 2]) - 1.0) < 1e-12)
check("too few points gives nan, not a spurious 1",
      bb.spearman([1, 2], [1, 2]) != bb.spearman([1, 2], [1, 2]))
check("mismatched lengths give nan",
      bb.spearman([1, 2, 3], [1, 2]) != bb.spearman([1, 2, 3], [1, 2]))

# ── AUPRC ─────────────────────────────────────────────────────────────────
check("a perfect ranking is 1.0",
      close(bb.auprc([9, 8, 7, 1, 2, 3], [1, 1, 1, 0, 0, 0]), 1.0))
check("an all-tied score returns the prevalence, not 1.0",
      close(bb.auprc([5, 5, 5, 5], [1, 1, 0, 0]), 0.5))
# The tie rule is the trap: scoring each tied element separately walks a
# staircase that credits precision it has not earned.
check("a tied block is consumed whole, so ties cannot inflate the area",
      close(bb.auprc([5, 5, 1], [1, 0, 0]), 0.5),
      str(bb.auprc([5, 5, 1], [1, 0, 0])))
check("an inverted ranking scores well below prevalence",
      bb.auprc([1, 2, 3, 9, 8, 7], [1, 1, 1, 0, 0, 0]) < 0.5)
check("no positives gives nan rather than 0",
      bb.auprc([1, 2, 3], [0, 0, 0]) != bb.auprc([1, 2, 3], [0, 0, 0]))
check("all positives gives nan rather than 1",
      bb.auprc([1, 2, 3], [1, 1, 1]) != bb.auprc([1, 2, 3], [1, 1, 1]))
check("empty input gives nan", bb.auprc([], []) != bb.auprc([], []))
# Rare positives are the paper's case: 5 positives among 100.
rare_scores = list(range(100, 0, -1))
rare_labels = [1] * 5 + [0] * 95
check("a rare positive class ranked first still scores 1.0",
      close(bb.auprc(rare_scores, rare_labels), 1.0))

# ── the claim table ───────────────────────────────────────────────────────
ids = [c["id"] for c in bb.PAPER_CLAIMS]
check("claims are encoded as data", len(bb.PAPER_CLAIMS) >= 15,
      f"{len(bb.PAPER_CLAIMS)} claims")
check("claim ids are unique", len(ids) == len(set(ids)))
check("every claim cites the sentence it came from",
      all(c.get("source") for c in bb.PAPER_CLAIMS))
check("every citation names a page or figure",
      all(re.search(r"p\.\d+|Fig\.", c["source"]) for c in bb.PAPER_CLAIMS))
check("every claim names the screen it belongs to",
      all(c["screen"] in bb.SCREENS for c in bb.PAPER_CLAIMS))
check("every claim has a tolerance",
      all("tol" in c and c["tol"] >= 0 for c in bb.PAPER_CLAIMS))
# Counts stated exactly in the text must be reproduced exactly; a library
# that has 3,400 guides is not the library the paper describes.
for cid in ("ldlvar.n_guides", "ldlvar.n_variants", "ldlvar.n_negctrl",
            "ldlrcds.n_guides", "ldlrcds.n_negctrl"):
    c = next(x for x in bb.PAPER_CLAIMS if x["id"] == cid)
    check(f"{cid} is exact (tol 0)", c["tol"] == 0)
# Derived counts depend on a threshold over a stochastic fit, so they cannot
# be exact -- but a tolerance wide enough to pass anything tests nothing.
for cid in ("ldlvar.n_significant", "ldlrcds.n_significant"):
    c = next(x for x in bb.PAPER_CLAIMS if x["id"] == cid)
    check(f"{cid} has a tolerance that is loose but not vacuous",
          0 < c["tol"] <= 0.2 * c["value"], f"tol={c['tol']} of {c['value']}")

# ── comparison and reporting ──────────────────────────────────────────────
rows = bb.compare({})
check("nothing measured yet reports 'not measured', never 'agrees'",
      all(r["status"] == "not measured" for r in rows))
rows = bb.compare({"ldlvar.n_guides": 3455, "ldlvar.replicate_rho": 0.86,
                    "ldlvar.auprc_bean": 0.61})
byid = {r["id"]: r for r in rows}
check("an exact hit agrees", byid["ldlvar.n_guides"]["status"] == "agrees")
check("a value inside tolerance agrees",
      byid["ldlvar.replicate_rho"]["status"] == "agrees")
check("a value outside tolerance DIFFERS",
      byid["ldlvar.auprc_bean"]["status"] == "DIFFERS")
check("the delta is reported so the direction is visible",
      abs(byid["ldlvar.auprc_bean"]["delta"] + 0.29) < 1e-9)
text = bb.render(rows, {})
check("the report names the paper", "Ryu" in text and "2024" in text)
check("the report names the data DOI", bb.ZENODO_DOI in text)
check("the report lists what it cannot test",
      "cannot test" in text and "UK Biobank" in text)
check("untestable claims say WHY",
      all(len(why) > 10 for _, _, why in bb.NOT_TESTABLE))

# ── the models, and the agent's ability to reach them ────────────────────
# A hosted run reported that the reporter, the guide barcode, X_bcmatch and
# the accessibility model were all untestable, and recommended pulling raw
# SRA data. Every one of those IS in the paper's deposit. The agent said
# otherwise because no tool for this skill was declared, so it reasoned from
# the IGVF screens instead -- the same class of gap as an unexposed flag.
check("all three of the paper's model variants are defined",
      set(bb.MODELS) == {"bean", "reporter", "uniform"})
check("the full model asks for accessibility scaling",
      "--scale-by-acc" in bb.MODELS["bean"]["flags"])
check("the uniform fallback is the one that drops the reporter",
      set(bb.MODELS["uniform"]["flags"]) == {"--uniform-edit", "--ignore-bcmatch"})
check("BEAN-Reporter sits between them, with neither flag",
      bb.MODELS["reporter"]["flags"] == [])
check("each model maps to the published claim it should be compared with",
      {m["claim"] for m in bb.MODELS.values()} ==
      {"auprc_bean", "auprc_bean_reporter", "auprc_bean_uniform"})
for m in bb.MODELS.values():
    cid = "ldlvar." + m["claim"]
    check(f"{m['claim']} is a claim that actually exists",
          any(c["id"] == cid for c in bb.PAPER_CLAIMS))

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import _tools as _T  # noqa: E402
_bt = [t for t in _T._TOOLS if t.name == "bean_paper_benchmark"]
check("the benchmark is reachable as an agent tool at all", bool(_bt))
if _bt:
    d = _bt[0].description
    check("its description says the IGVF screens are NOT the paper's",
          "NOT the paper's screens" in d)
    check("and that the deposit HAS the barcode/reporter layers",
          "X_bcmatch" in d and "DOES carry" in d)
    check("and pre-empts the wrong 'need raw SRA' conclusion",
          "PRJNA1042659" in d and "NOT needed" in d)
    argv = _T._build_argv(_bt[0], {"subcommand": "run", "screen": "ldlvar",
                                    "model": "all"})
    check("it builds a runnable command",
          argv[-3:] == ["ldlvar", "--model", "all"], " ".join(argv))

print(f"\n{len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
