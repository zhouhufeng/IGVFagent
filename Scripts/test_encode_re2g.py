"""Runs the ENCODE-rE2G port self-test (synthetic ABC world, planted model, all subcommands asserted).

    python3 Scripts/test_encode_re2g.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import encode_re2g_skill as re2g  # noqa: E402


def main() -> int:
    rc = re2g.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
