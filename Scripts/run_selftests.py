"""Run every offline self-test (``Scripts/test_*.py``) and summarise.

    python3 Scripts/run_selftests.py              # all, network blocked
    python3 Scripts/run_selftests.py --list       # what runs and what is skipped
    python3 Scripts/run_selftests.py test_abc_pipeline.py test_ui_panels.py

This is what the GitHub Actions workflow (.github/workflows/tests.yml) runs,
so the README's "tests" badge means exactly: every self-test below passed on
a clean install, with no network access and no credentials.

Network is blocked by pointing every proxy variable at a closed port, so a
test that quietly depends on a live service fails here instead of passing on
a developer's laptop. Tests that need something CI does not have are listed
in SKIP with the reason; nothing is skipped silently.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# test file -> why it does not run in CI. Keep this short and honest.
SKIP: "dict[str, str]" = {
}

DEAD_PROXY = "http://127.0.0.1:9"


def discover() -> "list[str]":
    return sorted(p.name for p in HERE.glob("test_*.py"))


def run_one(name: str, timeout: int) -> "tuple[int, float, str]":
    env = dict(os.environ)
    for k in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy",
              "ALL_PROXY", "all_proxy"):
        env[k] = DEAD_PROXY
    env["NO_PROXY"] = env["no_proxy"] = ""
    env.pop("IGVF_PROJECT_ROOT", None)
    env["PYTHONPATH"] = str(HERE) + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MPLBACKEND", "Agg")
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, str(HERE / name)], cwd=str(ROOT),
                           env=env, capture_output=True, text=True,
                           timeout=timeout)
        rc, out = p.returncode, (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired as e:
        rc = 124
        out = ((e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes)
               else (e.stdout or "")) + f"\nTIMEOUT after {timeout}s"
    return rc, time.time() - t0, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("tests", nargs="*", help="Run only these test files.")
    ap.add_argument("--list", action="store_true",
                    help="Show what would run and what is skipped, and why.")
    ap.add_argument("--timeout", type=int, default=600,
                    help="Per-test timeout in seconds (default 600).")
    args = ap.parse_args(argv)

    names = args.tests or discover()
    run = [n for n in names if n not in SKIP]
    skipped = [n for n in names if n in SKIP]
    if args.list:
        for n in run:
            print(f"  run   {n}")
        for n in skipped:
            print(f"  skip  {n}  — {SKIP[n]}")
        return 0

    failed = []
    for n in run:
        rc, secs, out = run_one(n, args.timeout)
        status = "ok  " if rc == 0 else "FAIL"
        print(f"  {status}  {n}  ({secs:.0f}s)", flush=True)
        if rc != 0:
            failed.append(n)
            tail = "\n".join(out.strip().splitlines()[-25:])
            print("\n".join("        " + l for l in tail.splitlines()), flush=True)
    for n in skipped:
        print(f"  skip  {n}  — {SKIP[n]}")
    print(f"\n{len(run) - len(failed)}/{len(run)} self-tests passed"
          f"{f', {len(skipped)} skipped' if skipped else ''}.")
    if failed:
        print("failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
