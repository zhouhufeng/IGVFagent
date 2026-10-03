"""Runs the perturb-tools self-test (synthetic dropout, sorting, PoolQ and guide-design data; all subcommands asserted).

    python3 Scripts/test_perturb_tools.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import perturb_tools_skill as pt  # noqa: E402


def main() -> int:
    rc = pt.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
