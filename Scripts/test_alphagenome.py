"""Offline self-test of the AlphaGenome skill (fake backend; no key, no network).

    python3 Scripts/test_alphagenome.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import alphagenome_skill as ag  # noqa: E402


def main() -> int:
    rc = ag.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
