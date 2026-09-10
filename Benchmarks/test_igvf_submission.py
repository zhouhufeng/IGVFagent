#!/usr/bin/env python3
"""The IGVF Portal submission tooling, absorbed from IGVF-Submit.

The value of these checks is that they encode findings from the monthly IGVF
submission meetings (Oct 2025 - Aug 2026), where the same properties were
missing every month. A rule that silently stops firing is worse than no rule,
because the submitter believes they were checked. So the tests assert the
rules exist, cite their source, and -- for the two pieces of real reasoning,
the allele crosscheck and the array repoint -- that they behave correctly on
constructed inputs.

No network: every case runs against fixtures or pure functions.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import igvf_submission_skill as sub   # noqa: E402
import igvf_portal_qc_skill as qc     # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:60} {detail}")
    if not ok:
        FAILURES.append(name)


# ── credentials: the absorption's one real behaviour change ───────────────
# The standalone tool read IGVF_API_KEY / IGVF_SECRET_KEY. IGVFagent stores
# the same pair as IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY -- what
# make-live.sh writes into .env.prod. Reading only the first pair reports
# "no credentials" on a deployment that holds them.
import os  # noqa: E402


def _with_env(**kw):
    saved = {k: os.environ.get(k) for k in kw}
    for k, v in kw.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    return saved


def _restore(saved):
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


ALL_KEYS = dict(IGVF_API_KEY=None, IGVF_SECRET_KEY=None,
                 IGVF_ACCESS_KEY=None, IGVF_SECRET_ACCESS_KEY=None)

saved = _with_env(**{**ALL_KEYS, "IGVF_API_KEY": "A", "IGVF_SECRET_KEY": "B"})
check("igvf_utils names are honoured", sub._portal_key_pair() == ("A", "B"))
_with_env(**{**ALL_KEYS, "IGVF_ACCESS_KEY": "C",
              "IGVF_SECRET_ACCESS_KEY": "D"})
check("IGVFagent's own names are honoured too",
      sub._portal_key_pair() == ("C", "D"))
_with_env(**{**ALL_KEYS, "IGVF_API_KEY": "A", "IGVF_SECRET_KEY": "B",
              "IGVF_ACCESS_KEY": "C", "IGVF_SECRET_ACCESS_KEY": "D"})
check("igvf_utils names win when both are set, so a submitter following the "
      "IGVF-Submit README is unaffected",
      sub._portal_key_pair() == ("A", "B"))
_with_env(**{**ALL_KEYS, "IGVF_API_KEY": "A"})
check("half a pair does not authenticate",
      not sub.Portal(mode="prod", api_key="A", secret_key="").authenticated)
check("a lone key does not authenticate a Portal",
      not sub.Portal(mode="prod", api_key="A", secret_key=None).authenticated)
check("an explicitly passed pair is not overwritten by the environment",
      (sub.Portal(mode="prod", api_key="X", secret_key="Y").api_key,
       sub.Portal(mode="prod", api_key="X", secret_key="Y").secret_key)
      == ("X", "Y"))
_s2 = _with_env(**{**ALL_KEYS, "IGVF_ACCESS_KEY": "envK",
                    "IGVF_SECRET_ACCESS_KEY": "envS"})
check("passing neither picks up the ambient pair",
      (sub.Portal(mode="prod").api_key,
       sub.Portal(mode="prod").secret_key) == ("envK", "envS"))
check("and an explicit pair still beats the ambient one",
      sub.Portal(mode="prod", api_key="X", secret_key="Y").api_key == "X")
_restore(_s2)
_restore(saved)

# ── sandbox is deprecated ─────────────────────────────────────────────────
# api.sandbox.igvf.org answers HTTP 410. Submission docs and the meeting notes
# still say sandbox, so the intent is redirected rather than left to fail.
check("sandbox redirects to staging", sub.DEPRECATED_MODES["sandbox"] == "staging")
p = sub.Portal(mode="sandbox")
check("a sandbox Portal ends up on staging", p.mode == "staging")
check("and its base host is the staging API",
      p.base == "https://api.staging.igvf.org")
check("prod is untouched", sub.Portal(mode="prod").base ==
      "https://api.data.igvf.org")

# ── dead statuses block release ───────────────────────────────────────────
for st in ("revoked", "deleted", "replaced", "archived"):
    check(f"{st!r} counts as an unusable input", st in sub.DEAD_STATUSES)
check("'released' does not", "released" not in sub.DEAD_STATUSES)
check("'in progress' does not", "in progress" not in sub.DEAD_STATUSES)

# ── every rule cites the meeting that produced it ─────────────────────────
src = Path(sub.__file__).read_text()
rules = set(re.findall(r'Finding\(\s*\n?\s*"([a-z0-9-]+)"', src))
check("the checklist has the rules the meetings produced",
      len(rules) >= 15, f"{len(rules)} distinct rules")
for want in ("description", "input-file-sets", "derived-from", "dead-input",
              "reference-files", "file-format-spec" if
              "file-format-spec" in rules else "format-spec",
              "analysis-step-version", "gzip", "bed3"):
    check(f"rule {want!r} is present", want in rules)
# A fail-level finding a submitter cannot act on is a dead end.
fails = re.findall(r'Finding\(\s*\n?\s*"([a-z0-9-]+)",\s*"fail",(.{0,400}?)\)\)',
                    src, re.S)
missing_fix = [r for r, body in fails if '"fix' not in body and
                body.count('",') < 2]
check("every failure carries a fix or a source",
      len(re.findall(r'meeting', src)) >= 10,
      f"{len(re.findall(r'meeting', src))} meeting citations")

# ── the allele crosscheck: the one piece of real reasoning ────────────────
# A reference-allele mismatch is corrected IN PLACE. The coordinate does not
# move while ref/alt change, so comparing positions alone reports exactly the
# case that matters as "no change".
old = {("chr1", 100): {("A", "G")}, ("chr1", 200): {("C", "T")}}
new = {("chr1", 100): {("G", "A")}, ("chr1", 200): {("C", "T")}}
dropped, changed = sub.revised_positions(old, new)
check("an in-place allele swap is detected as a revision",
      ("chr1", 100) in changed)
check("it is reported as changed, not as dropped", not dropped)
check("an unchanged position is not reported as revised",
      ("chr1", 200) not in changed)

# The Aug 2026 meeting note is precisely this case: "10 alleles mappings are
# incorrect (swapped ref and alt seqs)". A swap must count as a revision by
# default, or the submitter is told their file is unaffected when the input's
# alleles were rewritten under it.
check("a ref/alt swap counts as a revision by default",
      ("chr1", 100) in sub.revised_positions(old, new)[1])
check("--ignore-swaps is opt-in, for an analysis a swap cannot affect",
      ("chr1", 100) not in
      sub.revised_positions(old, new, ignore_swaps=True)[1])
check("normalise_pair makes a swap compare equal",
      sub.normalise_pair(("A", "G")) == sub.normalise_pair(("G", "A")))
check("but it does not collapse genuinely different alleles",
      sub.normalise_pair(("A", "G")) != sub.normalise_pair(("A", "T")))

gone = {("chr1", 100): {("A", "G")}}
d2, c2 = sub.revised_positions(old, gone)
check("a position removed entirely is reported as dropped",
      ("chr1", 200) in d2 and not c2)
check("comparing a map to itself finds nothing",
      sub.revised_positions(old, dict(old)) == (set(), set()))
# A position with no allele information cannot be compared, and guessing
# would manufacture a revision that the data does not support.
noinfo = {("chr1", 100): set(), ("chr1", 200): {("C", "T")}}
check("a position with no allele info is not called revised",
      ("chr1", 100) not in sub.revised_positions(noinfo, old)[1])

# ── schema requirements hide inside disjunctions ──────────────────────────
# prediction_set has NO top-level `required`. Its requirements are in a
# oneOf: lab/award/file_set_type always, then samples OR donors. Reading only
# `required` reported "0 required", emitted an empty template for the profile
# submitters use most, and dropped `samples` -- the exact property the Aug
# 2026 meeting listed as a to-do.
always, groups = sub.schema_required({
    "oneOf": [{"required": ["lab", "award", "file_set_type", "samples"]},
               {"required": ["lab", "award", "file_set_type", "donors"]}]})
check("shared requirements across a oneOf are always required",
      sorted(always) == ["award", "file_set_type", "lab"])
check("the branch-specific ones become an either/or group",
      groups and sorted(sum(groups[0], [])) == ["donors", "samples"])
check("a top-level required list still works",
      sub.schema_required({"required": ["a", "b"]}) == (["a", "b"], []))
check("allOf contributes unconditionally, not as a choice",
      sub.schema_required({"allOf": [{"required": ["x"]},
                                       {"required": ["y"]}]}) == (["x", "y"], []))
check("a single-branch anyOf is not an either/or",
      sub.schema_required({"anyOf": [{"required": ["z"]}]}) == (["z"], []))
check("a profile with no requirements yields nothing, not a crash",
      sub.schema_required({}) == ([], []))
check("branches with no required key are ignored",
      sub.schema_required({"oneOf": [{"properties": {}}]}) == ([], []))
check("duplicates are not repeated",
      sub.schema_required({"required": ["a"],
                            "allOf": [{"required": ["a"]}]}) == (["a"], []))

# ── the tool must not recommend the host it documents as dead ─────────────
skill_src = Path(sub.__file__).read_text()
bad = [l.strip() for l in skill_src.splitlines()
       if "-m sandbox" in l and "deprecated" not in l]
check("nothing recommends rehearsing on sandbox", not bad, str(bad[:1]))
check("staging is what gets recommended", "-m staging" in skill_src)

print(f"\n{len(rules)} rules encoded; {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
