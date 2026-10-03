"""Runs the CRISPRi-FlowFISH pipeline self-test (synthetic screen with planted enhancers, all subcommands asserted).

    python3 Scripts/test_flowfish_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import flowfish_pipeline_skill as ffp  # noqa: E402


def main() -> int:
    rc = ffp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
