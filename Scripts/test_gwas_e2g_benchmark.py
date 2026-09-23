"""Runs the GWAS E2G benchmark self-test (synthetic genome, all subcommands asserted).

    python3 Scripts/test_gwas_e2g_benchmark.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gwas_e2g_benchmark_skill as gw  # noqa: E402


def main() -> int:
    rc = gw.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
