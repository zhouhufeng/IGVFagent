"""Runs the bulk_crispr_pipeline (Perturb-seq) self-test (synthetic two-lane screen, all subcommands asserted).

    python3 Scripts/test_bulk_crispr_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bulk_crispr_pipeline_skill as bc  # noqa: E402


def main() -> int:
    rc = bc.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
