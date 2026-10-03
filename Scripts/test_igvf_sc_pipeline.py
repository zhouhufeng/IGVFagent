"""Runs the IGVF single-cell pipeline port self-test (synthetic reads, fragments, matrices; all subcommands asserted).

    python3 Scripts/test_igvf_sc_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import igvf_sc_pipeline_skill as sc  # noqa: E402


def main() -> int:
    rc = sc.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
