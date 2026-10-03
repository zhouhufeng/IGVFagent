#!/usr/bin/env python3
"""Regenerate the counts shown on the signed-out landing page.

The page used to carry "85 skills / 258 typed tools" as literal text in
Deploy/auth/gate.py. Those numbers were right when they were typed and wrong
within weeks — by the time anyone noticed, the real figures were 88 and 270.
A number a human has to remember to update is a number that will be wrong.

The gate runs in its own stdlib-only container and cannot import the agent
package to count for itself, so the counts are written here into a small JSON
file that ships with the gate image. ``test_landing_stats.py`` asserts the
committed file still matches the live catalogue, which turns drift into a
failing test instead of a stale claim on the front page.

    python3 Scripts/gen_landing_stats.py          # rewrite the file
    python3 Scripts/gen_landing_stats.py --check  # report drift, change nothing
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "Deploy" / "auth" / "stats.json"

sys.path.insert(0, str(HERE))


def current() -> dict:
    """Counts of the BUILT-IN catalogue, read from source.

    Deliberately not ``len(_tools._TOOLS)``. That list has agent-authored
    extensions merged into it at import time, so its length depends on which
    machine you ask: 270 in a fresh checkout, 338 on the deployment, and the
    landing page would advertise whichever host last regenerated the file.

    Counting ``_T(`` declarations in the source instead gives the same answer
    everywhere. It is also the more honest number to put on a public page:
    extensions are written by the agent at runtime and are not reviewed — one
    of them is what reported a matrix maximum of 14 for a matrix whose maximum
    was 3,724.
    """
    tools_src = (HERE / "_tools.py").read_text()
    tools = re.findall(r'^    _T\(\n\s*"([a-z0-9_]+)"', tools_src, re.M)

    cli_src = (HERE / "cli.py").read_text()
    skills = re.findall(r'^\s{4}"([a-z0-9-]+)":\s*\(\s*"igvfagent\.',
                        cli_src, re.M)
    return {"skills": len(set(skills)), "tools": len(set(tools))}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    now = current()
    try:
        on_disk = json.loads(OUT.read_text())
    except (OSError, ValueError):
        on_disk = {}

    if "--check" in argv:
        drift = {k: (on_disk.get(k), v) for k, v in now.items()
                 if on_disk.get(k) != v}
        if drift:
            for k, (was, is_) in drift.items():
                print(f"  {k}: landing page says {was}, catalogue has {is_}")
            print("\nRun: python3 Scripts/gen_landing_stats.py")
            return 1
        print(f"landing page is current: {now['skills']} skills, "
              f"{now['tools']} tools")
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(now, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {now['skills']} skills, "
          f"{now['tools']} typed tools")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
