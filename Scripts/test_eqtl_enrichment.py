"""Runs the eQTL enrichment benchmark self-test (synthetic genome, all subcommands asserted).

    python3 Scripts/test_eqtl_enrichment.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eqtl_enrichment_skill as eq  # noqa: E402


def main() -> int:
    rc = eq.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
