"""Upstream provenance stays complete and honest.

Skills here reimplement or wrap published methods whose reference
implementations live in other people's repositories. That provenance used to
be prose in a docstring: which project, sometimes which licence, never which
revision. This asserts the manifest keeps up with the code.

    python3 Scripts/test_upstream.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import upstream_skill as U  # noqa: E402

FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def main() -> int:
    data = U.load()
    entries = {e["repo"]: e for e in data.get("upstreams", [])}
    declared = U.discover()

    print(f"\nmanifest: {len(entries)} upstream(s); "
          f"docstrings declare {len(declared)}")

    print("\nevery upstream named in a docstring is in the manifest")
    # Renames count as the same upstream. CRISPRi-FlowFISH-pipeline became
    # CRISPRi-FlowFISH and the old name 404s; a manifest that could not say so
    # would force a choice between a stale entry and losing the history.
    known = set(entries)
    for e in entries.values():
        known |= set(e.get("former_names") or [])
    missing = sorted(set(declared) - known)
    check("none missing", missing, [])
    if missing:
        print("      run: igvfagent upstream scan")

    print("\nevery entry records what it is and how we relate to it")
    for repo, e in sorted(entries.items()):
        if e.get("relationship") not in U.RELATIONSHIPS:
            check(f"{repo}: relationship is one of {U.RELATIONSHIPS}",
                  e.get("relationship"), "one of the known values")
        elif e.get("relationship") == "unclassified":
            check(f"{repo}: classified", "unclassified", "classified")

    print("\nevery entry is pinned to a revision")
    for repo, e in sorted(entries.items()):
        pin = e.get("pinned") or {}
        if not pin.get("ref"):
            check(f"{repo}: pinned", "no revision", "a commit or tag")

    print("\ncopyleft and unlicensed upstreams are not vendored")
    # An upstream we copied code from binds this repo to its licence. GPL,
    # AGPL and a missing licence are the cases where that matters most, so
    # they are asserted rather than trusted to review.
    for repo, e in sorted(entries.items()):
        lic = (e.get("license") or "").upper()
        risky = (not lic) or lic.startswith(("GPL", "AGPL", "CC-"))
        if risky and e.get("relationship") == "vendored":
            check(f"{repo} ({lic or 'no licence'}) is not vendored",
                  "vendored", "clean-room / reference / wraps-binary")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
