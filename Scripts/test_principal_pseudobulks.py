"""Runs the principal-pseudobulks self-test (synthetic lanes, cells and genes; every subcommand asserted).

    python3 Scripts/test_principal_pseudobulks.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import principal_pseudobulks_skill as pp  # noqa: E402


def main() -> int:
    rc = pp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
