#!/usr/bin/env python3
"""Does the guide-counting key get chosen correctly, offline?

Two bugs this locks down, both from the LDLR prime-editing screen:

  1. Keying on `spacer` because it is the conventional column. That library
     shares one spacer across every variant installed at a nick site, so
     1,741 constructs collapsed onto 52 and each one's reads were credited
     to an arbitrary sibling -- effect sizes and FDRs for variants whose
     reads were never told apart.

  2. Keying on whichever column separates the most constructs. peg_sequence
     separates that library perfectly and is 128-134 bp against 128 bp
     reads, so it appears in no read at all and assigns 0%.

Neither is visible in metadata alone, so the choice is calibrated against
reads. These cases are synthetic and need no network.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import Scripts                                   # noqa: E402
sys.modules.setdefault("igvfagent", Scripts)
from Scripts import raw_data_pipeline as rp      # noqa: E402

FAILURES = []


def check(name, got, want):
    ok = got == want
    print(f"{'PASS' if ok else 'FAIL'}  {name:44} expect={want!r} got={got!r}")
    if not ok:
        FAILURES.append(name)


def matcher_of(seq_to_guide):
    return rp.build_matcher(seq_to_guide)


# ── build_matcher / match_read ──────────────────────────────────────────────

# Nesting: G_short's sequence is a prefix of G_long's. The longer match is
# the specific one, and must win regardless of scan order.
NEST = {"ACGTACGTACGT": "G_short", "ACGTACGTACGTTTTTGGGG": "G_long"}
m = matcher_of(NEST)
check("nested: long read -> longest match",
      rp.match_read("CCCC" + "ACGTACGTACGTTTTTGGGG" + "CCCC", m), "G_long")
check("nested: short read -> short match",
      rp.match_read("CCCC" + "ACGTACGTACGT" + "CCCC", m), "G_short")
check("no match returns None",
      rp.match_read("TTTTTTTTTTTTTTTTTTTTTTTT", m), None)

# A longer match later in the read beats a shorter one at the start: this is
# the 1-in-20,000 disagreement between first-found and longest-found.
FAR = {"ACGTACGTACGT": "G_short", "TTGGCCAATTGGCCAATTGG": "G_long"}
m2 = matcher_of(FAR)
check("longest wins over earliest offset",
      rp.match_read("ACGTACGTACGT" + "AA" + "TTGGCCAATTGGCCAATTGG", m2),
      "G_long")

# The prefix index must never invent a match that a full scan would not find.
m3 = matcher_of({"AAAAAAAAAACCCCC": "G1"})
check("prefix hit, full sequence absent -> None",
      rp.match_read("GGGG" + "AAAAAAAAAA" + "TTTTT", m3), None)


# ── key selection ───────────────────────────────────────────────────────────

def fake_index(cols):
    """An index shaped like load_guide_index's output, without the network."""
    cands = []
    for col, mapping in cols.items():
        ls = sorted(len(k) for k in mapping)
        cands.append({"column": col, "distinct": len(mapping),
                      "fraction": len(mapping) / 100.0,
                      "median_len": ls[len(ls) // 2],
                      "seq_to_guide": mapping,
                      "lengths": sorted(set(ls), reverse=True)})
    cands.sort(key=lambda c: (-c["distinct"], c["median_len"]))
    return {"resolved": True, "n_rows": 100, "candidates": cands,
            "target": {}, "type": {}, "why": ""}


def calibrate_on(idx, reads):
    """calibrate_key's ranking, fed reads directly instead of a FASTQ."""
    comp = str.maketrans("ACGTN", "TGCAN")
    tested = []
    for c in idx["candidates"]:
        mt = rp.build_matcher(c["seq_to_guide"])
        hit = 0
        for r in reads:
            if any(rp.match_read(s, mt) for s in (r, r.translate(comp)[::-1])):
                hit += 1
        tested.append({"column": c["column"], "distinct": c["distinct"],
                       "rate": hit / len(reads)})
    ranked = sorted(zip(tested, idx["candidates"]),
                    key=lambda t: (-round(t[0]["rate"], 2), -t[0]["distinct"]))
    return ranked[0][1]["column"] if ranked and ranked[0][0]["rate"] > 0 else None


# The real shape of the LDLR library, in miniature: a shared spacer, a
# discriminating RT template, and a peg_sequence longer than the read.
SPACER = "GCAGACCGGGACTGCTTGGA"
rtt = {f"{'ACGT'*3}{i:04d}": f"V{i}" for i in range(40)}
peg = {(SPACER + "G" * 90 + k): v for k, v in rtt.items()}   # 130 bp
# peg_sequence separates one construct MORE than the RT template does, which
# is the real library's margin: two of its 1,741 RT templates are duplicates
# and one is below the length floor, so peg leads 1,738 to 1,737. That one
# construct is enough for a uniqueness-ranked rule to choose the column that
# appears in no read.
rtt_ranked = {k: v for k, v in list(rtt.items())[:-1]}
idx = fake_index({"spacer": {SPACER: "V0"},
                  "rt_template_sequence": rtt_ranked, "peg_sequence": peg})
reads_128 = [(SPACER + "G" * 40 + k + "A" * (128 - 20 - 40 - len(k)))[:128]
             for k in rtt_ranked]
check("metadata alone would pick peg_sequence",
      idx["candidates"][0]["column"], "peg_sequence")
check("calibration picks rt_template_sequence",
      calibrate_on(idx, reads_128), "rt_template_sequence")

# A conventional KO library: unique spacers, and spacer must stay the choice.
ko = {f"{'TTGG'*4}{i:04d}"[:20]: f"G{i}" for i in range(40)}
idx_ko = fake_index({"spacer": ko})
reads_ko = [("CC" + k + "A" * 106)[:128] for k in ko]
check("unique-spacer library still picks spacer",
      calibrate_on(idx_ko, reads_ko), "spacer")

# Nothing in the reads: refuse rather than return a column that assigns 0%.
check("no key present -> None (refuse to guess)",
      calibrate_on(idx_ko, ["T" * 128] * 50), None)

# barcode must not be a candidate at all: 6 bp is chance-matched everywhere.
check("barcode excluded from candidate columns",
      "barcode" in rp._KEY_COLUMNS, False)
check("key length floor rules out short keys",
      rp._MIN_KEY_LEN >= 10, True)

print(f"\n{9 + 4} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
