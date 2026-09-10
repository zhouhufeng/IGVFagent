#!/usr/bin/env python3
"""Base-editing-aware guide assignment, following crispr-bean's method.

A collaborator from the lab that produced IGVFDS6464SOVZ flagged our
analysis: "for our base editing screens, it is important to allow for
self-editing in gRNA assignment (allow for A2G edits in the gRNA when
aligning) ... the 30% alignment accuracy seems low". They were right, and
the measurement on their data is:

    exact match only     36.7% of reads assigned
    A>G aware            62.5%          (+25.8 points)
    C>T aware            36.8%          (+0.1, so it is not a CBE)

The method here is crispr-bean's, not an invention: GuideEditCounter
compares mask_sequence(read) against the masked library, normalising the
edited base to its product on BOTH sides, rather than allowing free
mismatches. That distinction is the point of several cases below -- free
mismatches would also absorb sequencing error and cross-map similar guides.

These cases are offline and synthetic.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import raw_data_pipeline as rp  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:56} {detail}")
    if not ok:
        FAILURES.append(name)


# ── masking is BEAN's mask_sequence ───────────────────────────────────────

check("A>G masking maps every A", rp.mask_sequence("AACGT", ("A", "G")) == "GGCGT")
check("C>T masking maps every C", rp.mask_sequence("AACCT", ("C", "T")) == "AATTT")
check("masking leaves other bases alone",
      rp.mask_sequence("TTTT", ("A", "G")) == "TTTT")
check("ABE and CBE are the known editors",
      set(rp.BASE_EDITS) == {"ABE", "CBE"}, str(sorted(rp.BASE_EDITS)))
check("ABE is A>G", rp.BASE_EDITS["ABE"] == ("A", "G"))
check("CBE is C>T", rp.BASE_EDITS["CBE"] == ("C", "T"))


# ── the editor can be read off the library's own guide names ──────────────

check("library naming its guides ABE is detected",
      rp.detect_base_editor(["gRNA__ABCA1_1__ABE_neg_p4_h1",
                              "gRNA__LDLR_2__ABE_pos"]) == "ABE")
check("library naming its guides CBE is detected",
      rp.detect_base_editor(["g1__CBE_x", "g2__CBE_y"]) == "CBE")
check("silent naming returns None rather than guessing",
      rp.detect_base_editor(["guide_1", "guide_2"]) is None)
check("a lone mention does not decide a mixed library",
      rp.detect_base_editor(["ABE_1"] + ["plain"] * 20) is None)


# ── matching: only the editor's own substitution is forgiven ──────────────

SP = "AACAGGTATTTGTGAACTTT"          # a real spacer from IGVFFI4591THXG
LIB = {SP: "g_target", "TTGGCCAATTGGCCAATTGG": "g_other"}
M = rp.build_masked_matcher(LIB, "ABE")

check("an unedited read matches with 0 edits",
      rp.match_read_masked("CCCC" + SP + "CCCC", M) == ("g_target", 0, False))

# One A>G self-edit. It must be at a position that actually holds an A:
# my first attempt edited SP[2], which is a C, so the matcher correctly
# refused it and the test was wrong rather than the code.
_a = SP.index("A")
one = SP[:_a] + "G" + SP[_a + 1:]
check("a single A>G self-edit matches and is counted",
      rp.match_read_masked("CC" + one + "CC", M) == ("g_target", 1, False),
      str(rp.match_read_masked("CC" + one + "CC", M)))

# Three A>G edits.
edited = SP.replace("A", "G", 3)
n_expected = 3
check("multiple A>G edits are counted",
      rp.match_read_masked("CC" + edited + "CC", M)[:2]
      == ("g_target", n_expected),
      str(rp.match_read_masked("CC" + edited + "CC", M)))

# The crucial negative: a substitution the editor CANNOT make is not forgiven.
notedit = SP[:5] + ("T" if SP[5] != "T" else "C") + SP[6:]
g, _n, _a = rp.match_read_masked("CC" + notedit + "CC", M)
check("a non-editor substitution is NOT forgiven",
      g is None, f"got {g}")

# A G>A change is the reverse of the edit and must not match either.
rev = SP.replace("G", "A", 1)
g2, _n2, _a2 = rp.match_read_masked("CC" + rev + "CC", M)
check("the reverse substitution (G>A) is not forgiven", g2 is None, f"got {g2}")

# A CBE matcher must not rescue an ABE read.
Mc = rp.build_masked_matcher(LIB, "CBE")
check("a CBE matcher does not accept an A>G read",
      rp.match_read_masked("CC" + one + "CC", Mc)[0] is None)

try:
    rp.build_masked_matcher(LIB, "XYZ")
    check("an unknown editor is rejected", False, "no error raised")
except ValueError:
    check("an unknown editor is rejected", True)


# ── ambiguity must be reported, never assigned ────────────────────────────

# Two guides differing ONLY at an A/G position collapse under A>G masking.
A1 = "AAAACCCCTTTTGGGGAAAA"
A2 = "GAAACCCCTTTTGGGGAAAA"          # differs from A1 only by A->G at pos 0
AMB = {A1: "g1", A2: "g2"}
Ma = rp.build_masked_matcher(AMB, "ABE")
check("the masked index counts the collision",
      Ma["masked_collisions"] == 1, f"{Ma['masked_collisions']}")
g, _n, amb = rp.match_read_masked("CC" + A2 + "CC", Ma)
check("a read compatible with two guides is flagged ambiguous", amb is True,
      f"guide={g} ambiguous={amb}")
# And an unambiguous read in the same library is still assigned.
g3, n3, amb3 = rp.match_read_masked("CC" + A1 + "CC", Ma)
check("an unambiguous read is still assigned",
      g3 == "g1" and amb3 is False, f"{g3} {amb3}")

check("no match returns (None, 0, False)",
      rp.match_read_masked("TTTTTTTTTTTTTTTTTTTTTTTT", M) == (None, 0, False))


# ── the wrong tool must refuse, not undercount ────────────────────────────
#
# The base-editor signal is in the GUIDE LIBRARY, not the MeasurementSet
# metadata: IGVFDS6464SOVZ's readout is "gRNA sequencing" and its assay is
# "CRISPR FACS screen", identical to an ordinary knockout screen. Routing
# alone therefore cannot catch it, so the tool that would undercount has to.

msg = rp.base_editing_refusal(["gRNA__ABCA1_1__ABE_neg", "gRNA__LDLR__ABE_pos"],
                               "IGVFFI4591THXG")
check("an ABE library produces a refusal", msg is not None)
check("the refusal names the editor and substitution",
      "ABE" in msg and "A>G" in msg)
check("it names the tool to use instead", "bean analyze" in msg)
check("it quantifies the cost", "36.7%" in msg and "62.5%" in msg)
check("it names the library", "IGVFFI4591THXG" in msg)
check("it offers the override", "--force-exact" in msg)

check("a CBE library also refuses",
      "C>T" in (rp.base_editing_refusal(["g__CBE_1", "g__CBE_2"]) or ""))
# The false-positive direction matters as much: prime-editing and ordinary
# libraries must not be blocked.
check("a non-base-editing library does NOT refuse",
      rp.base_editing_refusal(["G137K", "G137N", "G140A"]) is None)
check("a plain knockout library does NOT refuse",
      rp.base_editing_refusal(["sg_TP53_1", "sg_TP53_2"]) is None)
check("an empty library does NOT refuse",
      rp.base_editing_refusal([]) is None)


# ─── BEAN's control condition ──────────────────────────────────────────────
# `bean run sorting` normalises each sorted bin against an unsorted sample and
# defaults --control-condition to "bulk". Handed a screen without one it dies
# with ValueError: No sample has control label `bulk`, which reads like a
# misconfiguration but is a real property of the screen: 18loci_uptake
# (IGVFDS6464SOVZ) sorts bottom20/bottom40/top20/top40 and has no unsorted
# bin at all, while 0426_LDLR137-219 (IGVFDS5542IBUS) does have one. So the
# bulk bins must be COUNTED and included even when --tail selects only the
# 20% pair, and their absence must be reported as a property of the data.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import base_editing_screen as bes  # noqa: E402

# The decision is a pure function so it is testable on a host without BEAN
# installed -- run_bean reaches it only after create-screen has succeeded,
# because the BEAN screen object is worth writing even when `run` must be
# declined.
ctl, why = bes.control_condition(["bottom20", "top20"])
check("no bulk condition -> no control condition", ctl is None)
check("the reason names the missing bin", "bulk" in why and "unsorted" in why)
check("the reason lists what the screen does have",
      "bottom20" in why and "top20" in why)
check("it does not blame configuration",
      "invalid" not in why.lower() and "error" not in why.lower())
check("it says the artefacts were still written",
      "were still written" in why)
check("it explains why a tail is not a substitute",
      "itself selected" in why)
# The refusal must not fire on a screen that HAS a bulk bin, or `bean run`
# never gets a chance on the one design it can actually model.
check("a bulk bin is chosen as the control",
      bes.control_condition(["bottom20", "bulk", "top20"])[0] == "bulk")
check("the portal's 'bulk (unsorted)' spelling is recognised",
      bes.control_condition(["top20", "bulk (unsorted)"])[0]
      == "bulk (unsorted)")
check("case is not load-bearing",
      bes.control_condition(["Bulk", "top20"])[0] == "Bulk")
check("an empty condition list has no control",
      bes.control_condition([])[0] is None)
check("None instead of a list does not crash",
      bes.control_condition(None)[0] is None)
# A bin merely CONTAINING "bulk" late in the name is not an unsorted sample;
# only a leading "bulk" is, which is how the portal spells it.
check("'nonbulk20' is not mistaken for a control",
      bes.control_condition(["nonbulk20", "top20"])[0] is None)

print(f"\n{21 + 10 + 12} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
