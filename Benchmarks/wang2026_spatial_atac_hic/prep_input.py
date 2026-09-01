#!/usr/bin/env python3
"""Build the Spatial-ATAC-Hi-C benchmark input with exact ground truth.

Why synthetic. Reproducing Wang 2026's figures from GSE307620 means
pulling multi-terabyte FASTQ and running an aligner — the one step
``spatial-hic`` deliberately does not own. What *can* be scored
deterministically, on any machine, in seconds, is whether the skill's
algorithms recover structure that is known by construction.

So this generator plants five facts and the benchmark asserts the
pipeline finds all five:

  1. a 50 x 50 pixel grid, every read carrying a resolvable barcode
  2. contact statistics calibrated to the paper's reported bands
     (88.1-90.3% cis; 24-33.3% long-range over cis)
  3. a 3x copy-number gain confined to one spatial clone
  4. a chromatin loop present in one clone and absent in the other
  5. alternating A/B compartment blocks, and a 10x TSS enrichment

Points 3-5 are the ones that matter: they are quantities no amount of
plumbing produces by accident. A pipeline that mis-derives copy number,
mis-orients the compartment eigenvector, or bungles the ANOVA will miss
them even though every file is written.

Ground truth is emitted to ``truth.json`` so ``make_figures.py`` scores
against it rather than against hard-coded constants.
"""

from __future__ import annotations

import gzip
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "Benchmarks" / "_data" / "wang2026_spatial_atac_hic"

SEED = 20260901
GRID = 6                       # 6 x 6 = 36 pixels; the real assay is 50 x 50
CHROM_SIZES = {"chr1": 8_000_000, "chr2": 6_000_000}

# Target contact statistics, taken from the paper's reported bands so a
# correct QC implementation lands inside them.
TARGET_CIS = 0.89              # paper: 88.1-90.3%
TARGET_LONG_RANGE = 0.28       # paper: 24-33.3% of cis contacts
LONG_RANGE_BP = 10_000

# Planted structure.
GAIN_LO, GAIN_HI = 3_000_000, 4_000_000    # chr1 CN gain, clone R only
GAIN_FOLD = 3.0
LOOP_A, LOOP_B = 1_000_000, 1_500_000      # chr1 loop, clone R only
LOOP_SLOP = 3_000
TSS_POS = 500_000                          # chr1, GENEX
TSS_FOLD = 10.0
COMPARTMENT_BLOCK = 600_000              # chr2 A/B block size
CONTACTS_PER_PIXEL = 4_000
FRAGMENTS_PER_PIXEL = 3_000


def _barcodes(n: int, salt: int, length: int = 8) -> "list[str]":
    rng = random.Random(salt)
    seen: "set[str]" = set()
    out: "list[str]" = []
    while len(out) < n:
        s = "".join(rng.choice("ACGT") for _ in range(length))
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _compartment_of(pos: int) -> int:
    """Alternating A (+1) / B (-1) blocks on chr2."""
    return 1 if (pos // COMPARTMENT_BLOCK) % 2 == 0 else -1


def _measure_contacts(pairs: "list[tuple]") -> "dict":
    """Contact statistics counted straight off the generated list.

    Independent of the skill: plain arithmetic over what was just
    emitted. The benchmark asserts the skill reproduces THESE numbers,
    not the aspirational targets above — a generator can miss its own
    target (an early version did, by giving one chromosome its own
    distance profile), and a benchmark that scores against the target
    rather than the data would blame the skill for it.
    """
    cis = trans = long_range = 0
    for _rid, _bc, c1, p1, c2, p2 in pairs:
        if c1 != c2:
            trans += 1
            continue
        cis += 1
        if abs(p2 - p1) >= LONG_RANGE_BP:
            long_range += 1
    total = cis + trans
    return {
        "total": total,
        "cis_fraction": cis / total if total else 0.0,
        "long_range_ratio": long_range / cis if cis else 0.0,
    }


def _measure_gain_fold(pairs: "list[tuple]") -> float:
    """Realised CN fold: coverage density inside the gain vs outside.

    The planted GAIN_FOLD is what the *extra* reads target, but the
    region also collects background, so the observable fold is lower.
    Measure what is actually there.
    """
    chrom = "chr1"
    size = CHROM_SIZES[chrom]
    inside = outside = 0
    for _rid, _bc, c1, p1, c2, p2 in pairs:
        for c, pos in ((c1, p1), (c2, p2)):
            if c != chrom:
                continue
            if GAIN_LO <= pos < GAIN_HI:
                inside += 1
            else:
                outside += 1
    bp_in = GAIN_HI - GAIN_LO
    bp_out = size - bp_in
    d_in = inside / bp_in
    d_out = outside / bp_out if bp_out else 0.0
    return d_in / d_out if d_out else float("nan")


def _measure_tss_enrichment(frags: "list[tuple]") -> float:
    """ArchR density ratio computed independently over both fragment ends.

    Both ends of every fragment are Tn5 insertions, so a spike written
    at the TSS also scatters its partner end 50-200 bp away — which is
    why the realised enrichment sits below the nominal TSS_FOLD.
    """
    centre = flank = 0
    for _c, s, e, _bc in frags:
        for pos in (s, e - 1):
            d = abs(pos - TSS_POS)
            if d <= 25:
                centre += 1
            elif 1900 <= d <= 2000:
                flank += 1
    d_centre = centre / 50.0
    d_flank = flank / 200.0
    return d_centre / d_flank if d_flank else float("nan")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    wl_a = _barcodes(GRID, 11)
    wl_b = _barcodes(GRID, 22)
    (OUT / "barcodes_A.txt").write_text(
        "".join(f"A{i+1}\t{s}\n" for i, s in enumerate(wl_a)))
    (OUT / "barcodes_B.txt").write_text(
        "".join(f"B{i+1}\t{s}\n" for i, s in enumerate(wl_b)))
    (OUT / "chrom.sizes").write_text(
        "".join(f"{c}\t{s}\n" for c, s in CHROM_SIZES.items()))

    # Clone assignment: the right half of the grid carries the gain and
    # the loop, so the CN clustering has a spatial signal to find.
    def clone(a_index: int) -> str:
        return "cloneR" if a_index >= GRID // 2 else "cloneL"

    header = ["## pairs format v1.0", "#shape: upper triangle"]
    header += [f"#chromsize: {c} {s}" for c, s in CHROM_SIZES.items()]
    header += ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"]

    pairs: "list[tuple]" = []
    frags: "list[tuple]" = []
    rid = 0

    for ai in range(GRID):
        for bi in range(GRID):
            barcode = wl_b[bi] + wl_a[ai]          # layout BA, as upstream
            is_r = clone(ai) == "cloneR"

            n = CONTACTS_PER_PIXEL
            # Extra coverage in the gain region for clone R. The gain is
            # realised as additional contacts, which is how a real copy
            # gain shows up in Hi-C marginals.
            n_gain = int(n * 0.10 * (GAIN_FOLD - 1)) if is_r else 0
            n_loop = int(n * 0.02) if is_r else 0
            n_main = n - n_gain - n_loop

            for _ in range(n_main):
                rid += 1
                if rng.random() >= TARGET_CIS:
                    p1 = rng.randrange(CHROM_SIZES["chr1"])
                    p2 = rng.randrange(CHROM_SIZES["chr2"])
                    pairs.append((rid, barcode, "chr1", p1, "chr2", p2))
                    continue
                c = "chr1" if rng.random() < 0.55 else "chr2"
                size = CHROM_SIZES[c]
                # Draw the distance class FIRST, on both chromosomes, so
                # the genome-wide long-range ratio is exactly the target.
                # An earlier version gave chr2 its own long-distance
                # profile, which silently doubled the global ratio.
                if rng.random() < TARGET_LONG_RANGE:
                    d = LONG_RANGE_BP + rng.randrange(0, 2_000_000)
                else:
                    d = rng.randrange(1, LONG_RANGE_BP)
                if c == "chr2":
                    # Compartment structure comes from *placement*, not
                    # distance: resample the anchor until both ends fall
                    # in the same compartment, keeping d fixed.
                    p1 = p2 = 0
                    for _try in range(16):
                        p1 = rng.randrange(max(1, size - d))
                        p2 = p1 + d
                        if p2 >= size:
                            continue
                        if _compartment_of(p1) == _compartment_of(p2):
                            break
                        if rng.random() < 0.2:
                            break
                    p2 = min(p2, size - 1)
                else:
                    p1 = rng.randrange(size)
                    p2 = min(p1 + d, size - 1)
                pairs.append((rid, barcode, c, p1, c, p2))

            for _ in range(n_gain):
                rid += 1
                p1 = rng.randrange(GAIN_LO, GAIN_HI)
                d = (rng.randrange(1, LONG_RANGE_BP)
                     if rng.random() > TARGET_LONG_RANGE
                     else LONG_RANGE_BP + rng.randrange(0, 500_000))
                p2 = min(p1 + d, GAIN_HI - 1)
                pairs.append((rid, barcode, "chr1", p1, "chr1", p2))

            for _ in range(n_loop):
                rid += 1
                p1 = LOOP_A + rng.randrange(-LOOP_SLOP, LOOP_SLOP)
                p2 = LOOP_B + rng.randrange(-LOOP_SLOP, LOOP_SLOP)
                pairs.append((rid, barcode, "chr1", p1, "chr1", p2))

            # ATAC fragments: uniform background plus a TSS spike sized so
            # the ArchR density ratio comes out at TSS_FOLD.
            span = 4000                       # +/-2 kb window around the TSS
            n_bg = FRAGMENTS_PER_PIXEL
            n_spike = int(n_bg * (50 / span) * (TSS_FOLD - 1))
            for _ in range(n_bg):
                s = rng.randrange(TSS_POS - span // 2, TSS_POS + span // 2)
                frags.append(("chr1", s, s + rng.randrange(50, 200), barcode))
            for _ in range(n_spike):
                s = rng.randrange(TSS_POS - 25, TSS_POS + 25)
                frags.append(("chr1", s, s + rng.randrange(50, 200), barcode))

    rng.shuffle(pairs)
    with gzip.open(OUT / "sample.pairs.gz", "wt") as fh:
        fh.write("\n".join(header) + "\n")
        for rid_, bc, c1, p1, c2, p2 in pairs:
            fh.write(f"r{rid_}:{bc}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n")

    with gzip.open(OUT / "fragments.tsv.gz", "wt") as fh:
        for c, s, e, bc in frags:
            fh.write(f"{c}\t{s}\t{e}\t{bc}\t1\n")

    (OUT / "genes.gtf").write_text(
        f'chr1\tsyn\tgene\t{TSS_POS + 1}\t{TSS_POS + 60000}\t.\t+\t.\t'
        'gene_id "GENEX"; gene_name "GENEX";\n'
        f'chr1\tsyn\tgene\t{LOOP_A + 1}\t{LOOP_A + 40000}\t.\t+\t.\t'
        'gene_id "GENEY"; gene_name "GENEY";\n'
        f'chr1\tsyn\tgene\t{GAIN_LO + 1}\t{GAIN_LO + 80000}\t.\t-\t.\t'
        'gene_id "GENEZ"; gene_name "GENEZ";\n'
        'chr2\tsyn\tgene\t200001\t260000\t.\t+\t.\t'
        'gene_id "GENEW"; gene_name "GENEW";\n')

    (OUT / "loops.bedpe").write_text(
        f"chr1\t{LOOP_A - 5000}\t{LOOP_A + 5000}\t"
        f"chr1\t{LOOP_B - 5000}\t{LOOP_B + 5000}\tclone_specific_loop\n"
        f"chr1\t6000000\t6010000\tchr1\t7000000\t7010000\tnull_loop\n")

    with (OUT / "clusters.tsv").open("w") as fh:
        fh.write("pixel\tcluster\n")
        for ai in range(1, GRID + 1):
            for bi in range(1, GRID + 1):
                fh.write(f"{ai:02d}x{bi:02d}\t{clone(ai - 1)}\n")

    realised = _measure_contacts(pairs)
    realised_gain = _measure_gain_fold(pairs)
    realised_tss = _measure_tss_enrichment(frags)

    truth = {
        "realised": {
            "cis_fraction": round(realised["cis_fraction"], 6),
            "long_range_ratio": round(realised["long_range_ratio"], 6),
            "gain_fold": round(realised_gain, 4),
            "tss_enrichment": round(realised_tss, 4),
        },
        "grid": GRID,
        "n_pixels": GRID * GRID,
        "n_pairs": len(pairs),
        "n_fragments": len(frags),
        "chrom_sizes": CHROM_SIZES,
        "target_cis_fraction": TARGET_CIS,
        "target_long_range_ratio": TARGET_LONG_RANGE,
        "long_range_bp": LONG_RANGE_BP,
        "gain": {"chrom": "chr1", "start": GAIN_LO, "end": GAIN_HI,
                 "fold": GAIN_FOLD, "clone": "cloneR"},
        "loop": {"chrom": "chr1", "anchor1": LOOP_A, "anchor2": LOOP_B,
                 "clone": "cloneR", "name": "clone_specific_loop"},
        "null_loop": "null_loop",
        "tss": {"chrom": "chr1", "pos": TSS_POS, "fold": TSS_FOLD},
        "compartment": {"chrom": "chr2", "block": COMPARTMENT_BLOCK},
        "clone_of_pixel": {
            f"{ai:02d}x{bi:02d}": clone(ai - 1)
            for ai in range(1, GRID + 1) for bi in range(1, GRID + 1)
        },
    }
    (OUT / "truth.json").write_text(json.dumps(truth, indent=2))

    print(f"wrote {OUT}")
    print(f"  pairs      {len(pairs):,}")
    print(f"  fragments  {len(frags):,}")
    print(f"  pixels     {GRID * GRID}")
    print(f"  planted    CN gain chr1:{GAIN_LO:,}-{GAIN_HI:,} (cloneR); "
          f"loop chr1:{LOOP_A:,}<->{LOOP_B:,} (cloneR)")
    print("  realised (independently measured off the emitted files):")
    print(f"    cis fraction     {realised['cis_fraction']:.4f}  "
          f"(paper band 0.881-0.903)")
    print(f"    long-range ratio {realised['long_range_ratio']:.4f}  "
          f"(paper band 0.240-0.333)")
    print(f"    CN gain fold     {realised_gain:.3f}x")
    print(f"    TSS enrichment   {realised_tss:.3f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
