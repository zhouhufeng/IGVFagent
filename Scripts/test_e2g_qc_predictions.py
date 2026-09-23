"""Runs the E2G QC-and-Predictions self-test (synthetic datasets, every subcommand asserted).

    python3 Scripts/test_e2g_qc_predictions.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import e2g_qc_predictions_skill as eqp  # noqa: E402


def main() -> int:
    rc = eqp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
