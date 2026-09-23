"""Runs the ABC pipeline self-test (synthetic chromosome, planted enhancers and Hi-C loop, all subcommands asserted).

    python3 Scripts/test_abc_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abc_pipeline_skill as abcp  # noqa: E402


def main() -> int:
    rc = abcp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
