"""Runs the scE2G workbench self-test: patches, configs, checks and benchmark asserted.

    python3 Scripts/test_sce2g_workbench.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sce2g_workbench_skill as wb  # noqa: E402


def main() -> int:
    rc = wb.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
