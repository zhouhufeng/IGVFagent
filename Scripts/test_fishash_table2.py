"""Runs the Fishash Table 2 reproduction self-test (synthetic barnyards, all subcommands asserted).

    python3 Scripts/test_fishash_table2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fishash_table2_skill as ft  # noqa: E402


def main() -> int:
    rc = ft.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
