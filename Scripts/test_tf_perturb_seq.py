"""Runs the tf-perturb self-test: planted regulators, trait and enhancer SNP recovered.

    python3 Scripts/test_tf_perturb_seq.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tf_perturb_seq_skill as tfp  # noqa: E402


def main() -> int:
    rc = tfp.main(["selftest", "--no-plots"])
    print("all checks pass" if rc == 0 else "FAILED")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
