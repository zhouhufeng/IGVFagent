"""Offline self-test of the Portal lineage walk and processed-first plan.

    python3 Scripts/test_processed_first.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import processed_first_skill as pf  # noqa: E402


def main() -> int:
    rc = pf.main(["selftest"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
