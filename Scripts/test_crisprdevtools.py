"""Runs the crisprdevtools self-test (module scaffolding + layout checks in a tempdir).

    python3 Scripts/test_crisprdevtools.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crisprdevtools_skill as cdt  # noqa: E402


def main() -> int:
    rc = cdt.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
