"""Runs the CRISPR Jamboree 4 benchmark self-test (synthetic inference checkpoint + control pairs, all subcommands asserted).

    python3 Scripts/test_crispr_jamboree4.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crispr_jamboree4_skill as j4  # noqa: E402


def main() -> int:
    rc = j4.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
