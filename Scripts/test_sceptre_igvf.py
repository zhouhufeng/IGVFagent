"""Runs the sceptreIGVF port self-test (synthetic high- and low-MOI screens, all subcommands asserted).

    python3 Scripts/test_sceptre_igvf.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sceptre_igvf_skill as sk  # noqa: E402


def main() -> int:
    rc = sk.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
