#!/usr/bin/env python3
"""A listing of a path that does not exist must not be a dead end.

On a real run the model asked for `Docs/BaseEditingScreen/18loci_uptake` when
the analysis had written `Docs/BaseEditingScreen/<timestamp>_18loci_uptake`.
`artifact ls` exited 2, the orchestrator recorded "1 tool call did not
succeed", and the whole answer was reported as incomplete -- for a directory
that was one listing of the parent away.

An agent cannot learn the right name from a refusal. It can from a listing.
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Scripts"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:62} {detail}")
    if not ok:
        FAILURES.append(name)


_tmp = tempfile.TemporaryDirectory()
import os  # noqa: E402
# .resolve() matters: on macOS a tempdir under /var resolves to
# /private/var, and _pathguard resolves its root, so an unresolved root
# compares as "outside the workspace".
root = Path(_tmp.name).resolve()
os.environ["IGVF_PROJECT_ROOT"] = str(root)
(root / "Docs" / "BaseEditingScreen" / "20260911_033012_18loci_uptake").mkdir(parents=True)
(root / "Docs" / "BaseEditingScreen" / "ready_check").mkdir()
(root / "Docs" / "BaseEditingScreen" / "20260911_033012_18loci_uptake"
 / "summary.json").write_text("{}")

import _pathguard  # noqa: E402
_pathguard.project_root = lambda: root
import artifact_read_skill as ar  # noqa: E402
ar._root = lambda: root

# ── the exact miss from the live run ─────────────────────────────────────
out = ar.list_artifacts("Docs/BaseEditingScreen/18loci_uptake")
check("a missing path returns a result, not an exception", isinstance(out, dict))
check("it says the path does not exist", out.get("exists") is False)
check("it names the nearest existing directory",
      out.get("nearest_existing") == "Docs/BaseEditingScreen")
check("it lists what is actually there",
      "20260911_033012_18loci_uptake" in out.get("available", []))
check("it suggests the near-miss, which differs only by a timestamp prefix",
      "20260911_033012_18loci_uptake" in out.get("did_you_mean", []),
      str(out.get("did_you_mean")))
check("the note tells the caller what to do next",
      "does not exist" in out.get("note", "")
      and "available" in out.get("note", ""))
check("entries is empty rather than absent, so a caller can iterate safely",
      out.get("entries") == [])

# ── a path that exists still works exactly as before ─────────────────────
ok = ar.list_artifacts("Docs/BaseEditingScreen/20260911_033012_18loci_uptake")
check("an existing directory reports exists=True", ok.get("exists") is True)
check("and lists its files",
      any(e["name"] == "summary.json" for e in ok["entries"]))
check("the parent lists both run directories",
      len(ar.list_artifacts("Docs/BaseEditingScreen")["entries"]) == 2)

# ── the suggestion must not fire on an unrelated name ────────────────────
miss = ar.list_artifacts("Docs/BaseEditingScreen/zzz_unrelated")
check("an unrelated name gets no false suggestion",
      miss.get("did_you_mean") == [], str(miss.get("did_you_mean")))
check("but it still lists what is available",
      len(miss.get("available", [])) == 2)

# ── walking up more than one level ───────────────────────────────────────
deep = ar.list_artifacts("Docs/BaseEditingScreen/nope/deeper/still")
check("it walks up to the first directory that exists",
      deep.get("nearest_existing") == "Docs/BaseEditingScreen",
      str(deep.get("nearest_existing")))

# ── containment is still enforced ────────────────────────────────────────
try:
    ar.list_artifacts("/etc")
    outside = False
except PermissionError:
    outside = True
check("a path outside the workspace is still refused", outside)

print(f"\n{len(FAILURES)} failure(s)")
_tmp.cleanup()
sys.exit(1 if FAILURES else 0)
