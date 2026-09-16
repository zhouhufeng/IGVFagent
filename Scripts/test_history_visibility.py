"""Who can see whose work — the rules, asserted.

With accounts, "private per user, shared into projects" stops being a
description and becomes a thing that can be got wrong silently: a missing
filter leaks one lab's unpublished analysis into another's search results, and
an over-eager one hides the 356 sessions this deployment already holds. Both
failures are quiet. Hence this file.

    python3 Scripts/test_history_visibility.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

_TMP = Path(tempfile.mkdtemp(prefix="igvf-vis-"))
os.environ["IGVF_HISTORY_DB"] = str(_TMP / "history.sqlite")
os.environ["IGVF_PROJECT_ROOT"] = str(_TMP)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _history as H  # noqa: E402

FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def refs(rows) -> set:
    return {r.get("ref") or r.get("run_dir") for r in rows}


def main() -> int:
    # --- a session from before accounts existed, plus one each for two users
    H.record_session("Docs/Agent/20260101_000000_legacy_aaaaaaaa",
                     query="legacy shared-password run about IGVFDS0000AAAA",
                     answer="produced a legacy figure")
    H.record_session("Docs/Agent/20260201_000000_alice_bbbbbbbb",
                     query="alice studies IGVFDS1111BBBB",
                     answer="alice's private result", owner="alice")
    H.record_session("Docs/Agent/20260301_000000_bob_cccccccc",
                     query="bob studies IGVFDS2222CCCC",
                     answer="bob's private result", owner="bob")

    # Default policy: RESULTS are shared across the deployment, so the same
    # dataset is not re-analysed once per person, while ORGANISATION (who
    # filed what into which project) stays private.
    print("\nanswers are shared, so nobody pays twice for the same question")
    check("alice sees every recorded run, including bob's",
          refs(H.recent_sessions(viewer="alice")),
          {"Docs/Agent/20260101_000000_legacy_aaaaaaaa",
           "Docs/Agent/20260201_000000_alice_bbbbbbbb",
           "Docs/Agent/20260301_000000_bob_cccccccc"})
    check("an unauthenticated/local deployment sees everything too",
          len(H.recent_sessions(viewer=None)), 3)
    check("alice can search bob's answer text",
          len(H.search("bob's private result", viewer="alice")) >= 1, True)
    check("alice can recall bob's accession — this is the whole point",
          len(H.by_accession("IGVFDS2222CCCC", viewer="alice")), 1)
    check("...and open the run itself",
          H.session("Docs/Agent/20260301_000000_bob_cccccccc",
                    viewer="alice") is not None, True)
    check("the pre-accounts corpus stays readable",
          len(H.search("legacy shared-password", viewer="alice")) >= 1, True)

    print("\nIGVF_HISTORY_SHARED=0 restores strict per-user privacy")
    os.environ["IGVF_HISTORY_SHARED"] = "0"
    check("alice no longer sees bob's run",
          refs(H.recent_sessions(viewer="alice")),
          {"Docs/Agent/20260101_000000_legacy_aaaaaaaa",
           "Docs/Agent/20260201_000000_alice_bbbbbbbb"})
    check("...nor recalls his accession",
          H.by_accession("IGVFDS2222CCCC", viewer="alice"), [])
    check("...nor finds his answer text",
          H.search("bob's private result", viewer="alice"), [])
    check("bob still sees his own", len(H.by_accession("IGVFDS2222CCCC",
                                                        viewer="bob")), 1)
    os.environ["IGVF_HISTORY_SHARED"] = "1"
    check("shared mode restored",
          len(H.by_accession("IGVFDS2222CCCC", viewer="alice")), 1)

    print("\nprojects stay private even though answers are shared")
    H.create_project("Bob's study", "private", owner="bob")
    check("bob sees his project", [p["name"] for p in
          H.list_projects(viewer="bob")], ["Bob's study"])
    check("alice does not", H.list_projects(viewer="alice"), [])
    check("alice cannot resolve it by name",
          H.resolve_project("Bob's study", viewer="alice"), None)
    check("...so renaming it is refused",
          H.rename_project("Bob's study", "Hijacked", viewer="alice")["ok"],
          False)

    print("\nfiling into a project, then sharing it")
    H.add_item("Bob's study", "session",
               "Docs/Agent/20260301_000000_bob_cccccccc",
               title="bob's run", owner="bob", viewer="bob")
    check("the project itself is still invisible to alice",
          H.list_projects(viewer="alice"), [])
    check("...and its contents are not listable by her",
          H.project_items("Bob's study", viewer="alice"), [])
    check("alice cannot share a project she cannot see",
          H.share_project("Bob's study", "alice", viewer="alice")["ok"], False)
    check("bob shares it", H.share_project("Bob's study", "alice",
                                           viewer="bob")["ok"], True)
    check("alice now sees the project", [p["name"] for p in
          H.list_projects(viewer="alice")], ["Bob's study"])
    check("...and its contents",
          len(H.project_items("Bob's study", viewer="alice")), 1)

    print("\nsharing grants visibility, not control")
    H.share_project("Bob's study", "alice", viewer="bob")
    check("a member cannot rename the project",
          H.rename_project("Bob's study", "Hijacked", viewer="alice")["ok"],
          False)
    check("...and the name is untouched",
          H.resolve_project("Bob's study", viewer="bob")["name"],
          "Bob's study")
    check("a member cannot archive it",
          H.archive_project("Bob's study", viewer="alice")["ok"], False)
    check("a member cannot rewrite the description",
          H.set_description("Bob's study", "mine now", viewer="alice")["ok"],
          False)
    check("a member cannot add further members",
          H.share_project("Bob's study", "carol", viewer="alice")["ok"], False)
    check("but a member CAN file work into it",
          H.add_item("Bob's study", "note", "alice-contribution",
                     title="from alice", owner="alice", viewer="alice")["ok"],
          True)
    check("the owner can still rename",
          H.rename_project("Bob's study", "Bob's renamed study",
                           viewer="bob")["ok"], True)
    H.rename_project("Bob's renamed study", "Bob's study", viewer="bob")

    print("\nunsharing takes it back")
    check("only the owner may unshare",
          H.unshare_project("Bob's study", "alice", viewer="alice")["ok"],
          False)
    check("bob unshares", H.unshare_project("Bob's study", "alice",
                                            viewer="bob")["ok"], True)
    check("alice loses sight of the project", H.list_projects(viewer="alice"), [])
    check("...and of its contents",
          H.project_items("Bob's study", viewer="alice"), [])

    print("\nper-user active project")
    H.create_project("Alice's study", owner="alice")
    H.set_active("Alice's study", viewer="alice")
    H.set_active("Bob's study", viewer="bob")
    check("alice's active project is hers",
          (H.active_project(viewer="alice") or {}).get("name"), "Alice's study")
    check("bob's is his, unaffected",
          (H.active_project(viewer="bob") or {}).get("name"), "Bob's study")

    print("\nnew runs file into the runner's own active project")
    H.record_session("Docs/Agent/20260401_000000_alice2_dddddddd",
                     query="another alice run", answer="x", owner="alice")
    check("alice's run landed in alice's project",
          refs(H.project_items("Alice's study", viewer="alice")),
          {"Docs/Agent/20260401_000000_alice2_dddddddd"})
    check("...and not in bob's",
          "Docs/Agent/20260401_000000_alice2_dddddddd"
          in refs(H.project_items("Bob's study", viewer="bob")), False)

    print("\nhistory is still append-only")
    con = H.connect()
    try:
        con.execute("DELETE FROM sessions")
        check("deleting sessions is blocked", "allowed", "blocked")
    except Exception:
        check("deleting sessions is blocked", "blocked", "blocked")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
