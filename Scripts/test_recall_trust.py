"""Whether a recalled answer deserves to be trusted — asserted.

Recall serves past answers back automatically, so a wrong one does not stay a
single bad run: it becomes the premise every later question about that dataset
inherits. On this deployment 34% of recorded sessions did not finish cleanly,
which is the scale of the problem being defended against here.

    python3 Scripts/test_recall_trust.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

_TMP = Path(tempfile.mkdtemp(prefix="igvf-trust-"))
os.environ["IGVF_HISTORY_DB"] = str(_TMP / "history.sqlite")
os.environ["IGVF_PROJECT_ROOT"] = str(_TMP)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _history as H      # noqa: E402
import _agent as A        # noqa: E402

FAILED: "list[str]" = []
ACC = "IGVFDS6290WNNH"


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def main() -> int:
    print("\ntiering, from what the agent already records")
    cases = [
        ({"stop_reason": "complete", "tool_calls_failed": 0}, "", H.TIER_TRUSTED),
        ({"stop_reason": "complete", "tool_calls_failed": 2}, "", H.TIER_PARTIAL),
        ({"stop_reason": "complete_with_failures", "tool_calls_failed": 3}, "", H.TIER_PARTIAL),
        ({"stop_reason": "max_iterations_wrapped", "tool_calls_failed": 0}, "", H.TIER_PARTIAL),
        ({"stop_reason": "error", "tool_calls_failed": 1}, "", H.TIER_FAILED),
        ({"stop_reason": "protocol_violation", "tool_calls_failed": 0}, "", H.TIER_FAILED),
        ({"stop_reason": "unrecognised", "tool_calls_failed": 0}, "", H.TIER_PARTIAL),
        ({"stop_reason": "complete", "tool_calls_failed": 0}, "wrong", H.TIER_WRONG),
        ({"stop_reason": "error", "tool_calls_failed": 9}, "correct", H.TIER_TRUSTED),
    ]
    for row, verdict, want in cases:
        label = f"{row['stop_reason']}/{row['tool_calls_failed']}fail"
        if verdict:
            label += f" + human says {verdict}"
        check(label, H.session_tier(row, verdict), want)

    H.record_session("Docs/Agent/20260101_000000_clean_aaaa1111",
                     query=f"QC {ACC}", answer=f"CLEAN-ANSWER for {ACC}.",
                     meta={"stop_reason": "complete", "tool_calls_failed": 0})
    H.record_session("Docs/Agent/20260102_000000_partial_bbbb2222",
                     query=f"QC {ACC} again", answer=f"PARTIAL-ANSWER for {ACC}.",
                     meta={"stop_reason": "complete_with_failures",
                           "tool_calls_failed": 3})
    H.record_session("Docs/Agent/20260103_000000_broke_cccc3333",
                     query=f"QC {ACC} once more", answer=f"BROKEN-ANSWER for {ACC}.",
                     meta={"stop_reason": "error", "tool_calls_failed": 5})

    print("\nwhat recall actually serves")
    _, block = A.prior_results_for(f"what about {ACC}?")
    check("a clean run is recalled", "CLEAN-ANSWER" in block, True)
    check("a partial run is recalled too", "PARTIAL-ANSWER" in block, True)
    check("...but marked INCOMPLETE", "INCOMPLETE" in block, True)
    check("a run that errored is never recalled", "BROKEN-ANSWER" in block, False)
    check("withholding is disclosed, not silent", "withheld" in block, True)
    check("clean outranks partial regardless of recency",
          block.index("CLEAN-ANSWER") < block.index("PARTIAL-ANSWER"), True)
    check("the block calls itself unverified", "unverified" in block, True)

    print("\nmarking an answer wrong")
    res = H.flag_session("Docs/Agent/20260101_000000_clean_aaaa1111", "wrong",
                         reason="misread the assay", author="alice")
    check("the flag is accepted", res["ok"], True)
    _, block = A.prior_results_for(f"what about {ACC}?")
    check("it stops being served to the next person",
          "CLEAN-ANSWER" in block, False)
    check("but it is NOT deleted — still searchable",
          len(H.search("CLEAN-ANSWER")) >= 1, True)
    check("...and still retrievable in full",
          H.session("Docs/Agent/20260101_000000_clean_aaaa1111") is not None, True)
    check("the reason is on record for the next reader",
          H.latest_verdict("Docs/Agent/20260101_000000_clean_aaaa1111")["reason"],
          "misread the assay")

    print("\nconfirming a run promotes it above unconfirmed ones")
    H.flag_session("Docs/Agent/20260102_000000_partial_bbbb2222", "correct",
                   reason="checked by hand", author="bob")
    _, block = A.prior_results_for(f"what about {ACC}?")
    check("the confirmation is shown to the reader",
          "confirmed correct by bob" in block, True)

    print("\nverdicts are append-only like the rest of history")
    con = H.connect()
    try:
        con.execute("DELETE FROM session_verdicts")
        check("deleting a verdict is blocked", "allowed", "blocked")
    except Exception:
        check("deleting a verdict is blocked", "blocked", "blocked")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
