#!/usr/bin/env python3
"""Map sequencing reads to a guide library, tolerating imperfect matches.

The capability CRISPR-Correct (pinellolab) provides: raw FASTQ against a
guide-library table, with matching that survives sequencing error. IGVFagent
already had every primitive -- raw_data_pipeline.build_matcher (prefix index
over construct sequences), mask_sequence (BEAN's edit-aware normalisation),
share_seq_skill's 1-Hamming barcode maps -- but they were EMBEDDED inside the
screen pipelines. There was no way to point them at an arbitrary FASTQ and
library, which is exactly what CRISPR-Correct is for. This exposes them.

TWO KINDS OF IMPERFECT MATCH, and they are not interchangeable:

  * Sequencing error. A read differs from the library by a base that was
    misread. Tolerated with a Hamming distance budget.
  * Base editing. A read differs because the EDITOR changed a base (A>G for
    ABE, C>T for CBE). Handled by masking: normalise the edited base to its
    product in BOTH read and library, so an edited read still matches its own
    guide -- without spending mismatch budget on it.

Masking is not a substitute for a mismatch budget, and the codebase already
reasons about this at raw_data_pipeline.py:368: allowing free mismatches
"would also absorb sequencing error and cross-map". So both are offered
separately, and `--edit` is what you use for a base-editing library.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402
import raw_data_pipeline as rp                              # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "GuideMap"

EDITS = {"ABE": ("A", "G"), "CBE": ("C", "T")}


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_library(path: str, seq_col: Optional[str], id_col: Optional[str]):
    """Guide library table -> {sequence: guide_id}."""
    delim = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    rows = list(csv.DictReader(_open(path), delimiter=delim))
    if not rows:
        raise SystemExit(f"{path} is empty.")
    cols = list(rows[0].keys())

    def pick(given, *names):
        if given:
            if given not in cols:
                raise SystemExit(f"column {given!r} not in {cols}")
            return given
        for n in names:
            for c in cols:
                if c.lower() == n:
                    return c
        return None

    c_seq = pick(seq_col, "spacer", "protospacer", "sequence", "guide_seq", "sgrna")
    c_id = pick(id_col, "guide_id", "guide", "name", "id", "sgrna_id")
    if not c_seq:
        raise SystemExit(f"no guide-sequence column found in {cols}")
    out = {}
    for i, r in enumerate(rows):
        sq = (r.get(c_seq) or "").strip().upper()
        if not sq or set(sq) - set("ACGTN"):
            continue
        out[sq] = (r.get(c_id) or f"guide_{i}").strip()
    if not out:
        raise SystemExit(f"no usable guide sequences in {path}")
    return out, c_seq, c_id


def _hamming_ok(a: str, b: str, budget: int) -> bool:
    if len(a) != len(b):
        return False
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > budget:
                return False
    return True


def cmd_map(args) -> int:
    setup_logging()
    lib, c_seq, c_id = read_library(args.library, args.seq_col, args.id_col)
    logging.info("library: %d guides from column %r", len(lib), c_seq)

    edit = EDITS.get((args.edit or "").upper()) if args.edit else None
    if edit:
        logging.info("base-editing mode %s: masking %s>%s in read and library",
                     args.edit.upper(), *edit)
        masked = {}
        for sq, gid in lib.items():
            masked.setdefault(rp.mask_sequence(sq, edit), gid)
        search = masked
    else:
        search = dict(lib)

    lengths = sorted({len(s) for s in search}, reverse=True)
    counts: Counter = Counter()
    stats = Counter()
    n_reads = 0
    t0 = time.time()

    with _open(args.fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n_reads += 1
            if args.max_reads and n_reads > args.max_reads:
                n_reads -= 1
                break
            seq = line.strip().upper()
            probe = rp.mask_sequence(seq, edit) if edit else seq
            hit = None
            # Exact first, at every guide length, anchored at the offset the
            # library implies. Cheap, and it resolves the large majority.
            for L in lengths:
                if len(probe) < L:
                    continue
                for off in range(0, min(args.max_offset, len(probe) - L) + 1):
                    cand = probe[off:off + L]
                    if cand in search:
                        hit = search[cand]
                        stats["exact"] += 1
                        break
                if hit:
                    break
            if hit is None and args.mismatches > 0:
                for L in lengths:
                    if len(probe) < L:
                        continue
                    cand = probe[:L]
                    for sq, gid in search.items():
                        if len(sq) == L and _hamming_ok(cand, sq, args.mismatches):
                            hit = gid
                            stats["mismatch"] += 1
                            break
                    if hit:
                        break
            if hit is None:
                stats["unmapped"] += 1
            else:
                counts[hit] += 1

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_guidemap"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "guide_counts.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["guide_id", "count"])
        for gid, n in counts.most_common():
            w.writerow([gid, n])

    mapped = stats["exact"] + stats["mismatch"]
    summary = {"reads": n_reads, "mapped": mapped,
               "mapped_fraction": round(mapped / n_reads, 4) if n_reads else 0.0,
               "exact": stats["exact"], "mismatch_rescued": stats["mismatch"],
               "unmapped": stats["unmapped"],
               "guides_detected": len(counts), "guides_in_library": len(lib),
               "edit_mode": args.edit or "none",
               "mismatch_budget": args.mismatches,
               "seconds": round(time.time() - t0, 1)}
    js = out / "guide_map_summary.json"
    js.write_text(json.dumps(summary, indent=2))
    print(f"Report:  {tsv}")
    print(f"Summary: {js}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent guide-map",
        description="Map FASTQ reads to a guide library with imperfect "
                     "matching (sequencing error and/or base editing).")
    sub = p.add_subparsers(dest="command", required=True)
    m = sub.add_parser("map", help="FASTQ + library -> per-guide counts.")
    m.add_argument("--fastq", required=True)
    m.add_argument("--library", required=True, help="CSV/TSV guide library.")
    m.add_argument("--seq-col", help="Guide-sequence column; auto-detected.")
    m.add_argument("--id-col", help="Guide-id column; auto-detected.")
    m.add_argument("--mismatches", type=int, default=1,
                   help="Hamming budget for sequencing error (0 disables).")
    m.add_argument("--edit", choices=sorted(EDITS),
                   help="Base editor: mask its product so edited reads still "
                        "match without spending mismatch budget.")
    m.add_argument("--max-offset", type=int, default=4,
                   help="Search this many bases in from the read start.")
    m.add_argument("--max-reads", type=int, default=0)
    m.add_argument("--label")
    m.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"map": cmd_map}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
