#!/usr/bin/env python3
"""Regression test: every IGVF assay routes to the analysis it needs.

Routing used to be decided from file types, so anything with FASTQ reads
and no matrix went to transcriptome quantification -- correct for an RNA
readout, wrong for the other 40% of the Portal, and silent either way. This
pins the classification so a future change cannot quietly restore that.

Run: python3 Benchmarks/test_assay_routing.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import _assays as A                                            # noqa: E402

# Real datasets whose correct route is known, including the three that were
# each analysed wrongly before being fixed.
CASES = [
    ("IGVFDS3532MONX", {"preferred_assay_titles": ["TAP-seq"],
                         "crispr_screen_readout": "scRNA-seq"}, A.TRANSCRIPT),
    ("IGVFDS9875NBZW", {"preferred_assay_titles":
                         ["10x multiome with MULTI-seq"]}, A.TRANSCRIPT),
    ("IGVFDS6639ECQN", {"preferred_assay_titles": ["RNA-seq"]}, A.TRANSCRIPT),
    ("IGVFDS4629JYPY", {"preferred_assay_titles": ["SGE"]}, A.VARIANT),
    ("IGVFDS6464SOVZ", {"preferred_assay_titles": ["CRISPR FACS screen"],
                         "crispr_screen_readout": "gRNA sequencing"}, A.GUIDE),
    # readout must beat the title: same screen, opposite analysis
    ("facs-with-rna-readout", {"preferred_assay_titles": ["CRISPR FACS screen"],
                                "crispr_screen_readout": "scRNA-seq"},
     A.TRANSCRIPT),
    ("mpra", {"preferred_assay_titles": ["lentiMPRA"]}, A.ELEMENT),
    ("atac", {"preferred_assay_titles": ["ATAC-seq"]}, A.CHROMATIN),
    ("wgs", {"preferred_assay_titles": ["WGS"]}, A.GENOME),
    ("y2h", {"preferred_assay_titles": ["Pooled Y2H"]}, A.PROTEIN),
    # an assay with no rule must be UNKNOWN, never TRANSCRIPT
    ("invented", {"preferred_assay_titles": ["Totally-New-Assay-2027"]},
     A.UNKNOWN),
    ("nothing stated", {}, A.UNKNOWN),
]


def main() -> int:
    failures = 0
    for name, fs, expect in CASES:
        got = A.classify(fs)["route"]
        ok = got == expect
        failures += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name:24} expect={expect:11} got={got}")
    # The property that matters most: nothing unrecognised is ever treated
    # as a transcript library.
    for bogus in ({"preferred_assay_titles": ["No-Such-Assay"]}, {},
                  {"assay_titles": [""]}):
        if A.classify(bogus)["route"] == A.TRANSCRIPT:
            print(f"FAIL  unrecognised assay routed to TRANSCRIPT: {bogus}")
            failures += 1
    print(f"\n{len(CASES)} cases, {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
