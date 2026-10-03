"""Runs the CRISPR Jamboree 3 pipeline self-test (synthetic FASTQs / count matrices / hashing, all stages asserted).

    python3 Scripts/test_crispr_jamboree3.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crispr_jamboree3_skill as j3  # noqa: E402


def main() -> int:
    rc = j3.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
