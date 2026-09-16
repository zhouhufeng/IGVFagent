"""The landing page's numbers still match the catalogue.

The signed-out page carried "85 skills / 258 typed tools" as literal text.
They were right when typed and wrong within weeks — the real figures had
reached 88 and 270 before anyone noticed. A number a human must remember to
update is a number that will be wrong, so this turns drift into a failing test
rather than a stale claim on the front page.

    python3 Scripts/test_landing_stats.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_landing_stats as G  # noqa: E402

if __name__ == "__main__":
    rc = G.main(["--check"])
    if rc:
        print("\nThe landing page would show numbers that are no longer true.")
    raise SystemExit(rc)
