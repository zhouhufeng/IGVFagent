"""Runs the scE2G pipeline self-test (synthetic multiome clusters, every subcommand asserted).

    python3 Scripts/test_sce2g_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sce2g_pipeline_skill as sp  # noqa: E402


def main() -> int:
    rc = sp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
