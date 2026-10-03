"""Runs the CRISPR Jamboree 2 self-test (synthetic MuData + FASTQs, all subcommands asserted).

    python3 Scripts/test_crispr_jamboree2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crispr_jamboree2_skill as j2  # noqa: E402


def main() -> int:
    rc = j2.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
