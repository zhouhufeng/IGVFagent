"""Runs the Gasperini 2019 pipeline self-test (synthetic two-lane screen with planted knockdowns, all subcommands asserted).

    python3 Scripts/test_gasperini_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasperini_pipeline_skill as sk  # noqa: E402


def main() -> int:
    rc = sk.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
