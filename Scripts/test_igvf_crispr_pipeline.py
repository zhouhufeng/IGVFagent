"""Runs the IGVF CRISPR pipeline self-test (upstream example seqspecs + synthetic reads / AnnData, all subcommands asserted).

    python3 Scripts/test_igvf_crispr_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import igvf_crispr_pipeline_skill as icp  # noqa: E402


def main() -> int:
    rc = icp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
