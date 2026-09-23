"""Runs the CRISPR_Pipeline port self-test (synthetic Perturb-seq screen, every subcommand asserted).

    python3 Scripts/test_crispr_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crispr_pipeline_skill as cp  # noqa: E402


def main() -> int:
    rc = cp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
