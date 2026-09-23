"""Runs the CRISPR-FG jamboree self-test (synthetic two-lane screen, all subcommands asserted).

    python3 Scripts/test_crispr_fg_jamboree.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crispr_fg_jamboree_skill as cj  # noqa: E402


def main() -> int:
    rc = cj.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
