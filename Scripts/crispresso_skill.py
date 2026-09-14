#!/usr/bin/env python3
"""CRISPResso2 wrapper: genome-editing outcomes from amplicon reads.

CRISPResso2 (Clement et al., Nat Biotechnol 2019, pinellolab/CRISPResso2) is
run as a SEPARATE PROGRAM via subprocess -- not imported, not vendored. It
lives in /opt/bean-venv beside BEAN because both need numpy<2, and because
BEAN reads the reporter allele through CRISPResso2's aligner.

This wrapper exists for three reasons a bare `CRISPResso` call does not cover:

  * It finds the binary. CRISPResso is deliberately NOT in the app venv, so
    `shutil.which` inside the app process does not see it unless
    /opt/bean-venv/bin is on PATH. The lookup here checks both, and says
    which one it used.
  * It reports WHERE the answer landed. CRISPResso writes a directory of
    tables and figures; the artefact lines this prints are what the agent
    turns into readable output.
  * It surfaces the quantification table rather than leaving the caller to
    guess the filename, which changes with the run name.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "CRISPResso"

# The bean venv is where the Dockerfile installs it; PATH is the fallback for
# a local install that put it somewhere else.
_BEAN_VENV_BIN = Path("/opt/bean-venv/bin")


def crispresso_exe() -> "tuple[str | None, str]":
    cand = _BEAN_VENV_BIN / "CRISPResso"
    if cand.is_file():
        return str(cand), f"{cand} (bean venv)"
    found = shutil.which("CRISPResso")
    if found:
        return found, f"{found} (PATH)"
    return None, ("CRISPResso2 is not installed. It ships with the "
                   "crispr-bean layer: rebuild with "
                   "`IGVF_INSTALL_CRISPR_BEAN=1 bash Deploy/redeploy.sh`.")


def cmd_analyze(args) -> int:
    exe, detail = crispresso_exe()
    if not exe:
        print(detail)
        return 3
    print(f"CRISPResso2: {detail}")

    name = args.name or f"{time.strftime('%Y%m%d_%H%M%S')}_crispresso"
    out = Path(args.output_dir) if args.output_dir else (OUT_DIR / name)
    out.mkdir(parents=True, exist_ok=True)

    cmd = [exe, "--fastq_r1", args.fastq_r1,
           "--amplicon_seq", args.amplicon,
           "--name", name, "--output_folder", str(out)]
    if args.fastq_r2:
        cmd += ["--fastq_r2", args.fastq_r2]
    if args.guide:
        cmd += ["--guide_seq", args.guide]
    if args.base_editor:
        cmd += ["--base_editor_output"]
        # CRISPResso names these separately; a single "A,G" is friendlier for
        # a caller than remembering which flag is which.
        if args.conversion:
            frm, _, to = args.conversion.partition(",")
            if frm and to:
                cmd += ["--conversion_nuc_from", frm.strip(),
                        "--conversion_nuc_to", to.strip()]
    if args.extra_args:
        cmd += args.extra_args.split()

    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
    if r.returncode != 0:
        tail = (r.stderr or r.stdout or "").strip().splitlines()[-6:]
        print(f"CRISPResso exited {r.returncode}:")
        for line in tail:
            print("  " + line)
        return r.returncode

    # Announce the artefacts. The quantification table is the file a reader
    # actually wants, and its name carries the run name, so glob for it
    # rather than making the caller reconstruct it.
    run_dir = next((d for d in sorted(out.glob(f"CRISPResso_on_{name}*"))
                    if d.is_dir()), out)
    print(f"Report: {run_dir}")
    for pattern in ("CRISPResso_quantification_of_editing_frequency.txt",
                    "Quantification_window_nucleotide_percentage_table.txt",
                    "Nucleotide_frequency_table.txt",
                    "CRISPResso_mapping_statistics.txt"):
        for f in run_dir.glob(pattern):
            print(f"Wrote: {f}")
    return 0


def cmd_check(_args) -> int:
    exe, detail = crispresso_exe()
    print(f"CRISPResso2: {detail}")
    if not exe:
        return 3
    r = subprocess.run([exe, "--version"], capture_output=True, text=True,
                        timeout=120)
    ver = (r.stdout + r.stderr).strip().splitlines()
    print("  " + (ver[-1] if ver else "(no version line)"))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent crispresso",
        description="CRISPResso2: quantify genome-editing outcomes from "
                     "amplicon sequencing reads.")
    sub = p.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="Quantify editing from amplicon reads.")
    a.add_argument("--fastq-r1", required=True)
    a.add_argument("--fastq-r2")
    a.add_argument("--amplicon", required=True,
                   help="Reference amplicon SEQUENCE (not a path).")
    a.add_argument("--guide", help="sgRNA spacer, no PAM.")
    a.add_argument("--base-editor", action="store_true")
    a.add_argument("--conversion", help="e.g. 'A,G' for ABE, 'C,T' for CBE.")
    a.add_argument("--name")
    a.add_argument("--output-dir")
    a.add_argument("--extra-args", help="Further CRISPResso flags, verbatim.")
    a.add_argument("--timeout", type=int, default=3600)

    sub.add_parser("check", help="Is CRISPResso2 runnable, and which build?")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"analyze": cmd_analyze, "check": cmd_check}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
