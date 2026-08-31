"""MPRA oligo design skill — regions / variants / tiling → synthesis-ready oligos.

Python port of **MPRAOligoDesign** (Max Schubach, Berlin Institute of
Health, Computational Genome Biology / kircherlab —
https://github.com/kircherlab/MPRAOligoDesign, MIT licence,
DOI 10.5281/zenodo.18173304). The original is a Snakemake workflow; this
is a single dependency-free module exposing the same design logic as an
IGVFagent subcommand tree.

If you publish a design produced with this skill, cite the upstream work:

    Rosen JD, Vasanthakumari AD, Salomon K, de Lange N, Dash PM,
    Keukeleire P, Hassan A, Barrera A, Krupkin B, Oualline G, Kircher M,
    Love MI, Schubach M. Uniform processing and analysis of IGVF
    massively parallel reporter assay data with MPRAsnakeflow.
    Genome Research (2025). doi:10.1101/gr.281462.125

What the design does
--------------------
An MPRA library is a pool of short oligos, each carrying a candidate
regulatory sequence plus fixed cloning adapters. Getting from "here are
my regions / variants" to "here is the order form" means four decisions:

1. **Where to put the oligo window.** A region shorter than the oligo is
   centred and padded; a mid-length region gets two tiles pinned to its
   ends; anything longer is tiled from the centre outwards with a
   guaranteed minimum overlap so no position falls in a gap.
2. **What sequence to write.** For a variant design, every region that
   contains the variant yields a REF oligo and an ALT oligo differing
   only at that position — the paired contrast the assay actually reads.
3. **What to throw away.** Sequences that break synthesis or confound
   the readout: long homopolymers, the EcoRI / SbfI sites used for
   cloning, simple repeats, promoters (TSS), and CTCF motifs.
4. **What to bolt on.** The constant left/right adapters.

Upstream layout → this module
-----------------------------
  ==========================================  ==========================
  MPRAOligoDesign                             here
  workflow/rules/tiling.smk (awk strategies)  ``tile``
  workflow/scripts/tiling/centerTiling.py     ``tile`` (centred strategy)
  scripts/oligo_design/getSequencesInclVariants.py  ``design-variants``
  rule oligo_design_regions_getSequences       ``design-regions``
  scripts/oligo_design/filter.py + filterOligos.py  ``filter``
  rule oligo_design_add_adapters               ``adapters``
  workflow/scripts/kMerFilter.py               ``kmer-filter``
  rules/final_design.smk                       ``pipeline`` (final stage)

Deliberate deviations from upstream
-----------------------------------
* **No Snakemake, conda, bedtools, pysam, pyranges, vcfpy, pyfaidx or
  BioPython.** Everything is reimplemented on the standard library, so a
  design runs anywhere IGVFagent runs. bgzip output is gzip-compatible,
  so ``.gz`` inputs are read directly.
* **Homopolymer bug fixed.** Upstream ``nucleotideruns`` never compares
  the final run, so a sequence *ending* in its longest run is scored
  short — ``"ACCCC"`` returns 1 and a pure ``"AAAA"`` returns 0, letting
  homopolymers through the filter. :func:`longest_nucleotide_run` closes
  the run at end-of-sequence.
* **Most-centred region fixed.** Upstream computes
  ``abs(pos - 1 - start - length)``, which is the distance to the region
  *end*, not its centre, so ``use_most_centered_region`` picks the wrong
  region whenever a variant matches several. :func:`_centre_distance`
  uses ``abs((pos - 1) - (start + length / 2))``. Pass
  ``--legacy-centre-metric`` to reproduce upstream output exactly.

Subcommands
-----------
    tile             Region BED → tiled oligo windows (3 strategies)
    design-regions   Region BED + reference → oligo FASTA
    design-variants  VCF + region BED + reference → REF/ALT oligo FASTA
    filter           Apply the five design filters, write a filter log
    adapters         Add constant left/right adapters
    kmer-filter      Pick the least k-mer-redundant sequence per group
    pipeline         tile → design → filter → adapters → final design
    write-playbook   Write Docs/Skills/MPRA_OLIGO_DESIGN_SKILLS.md
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Iterator

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "MPRAOligoDesign"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
REFERENCE_DIR = DATA_DIR / "Reference" / "OligoDesign"

# Upstream defaults (config/example_config.yml + the click defaults).
DEFAULT_OLIGO_LENGTH = 200
DEFAULT_MIN_OVERLAP = 50
DEFAULT_MAX_HOMOPOLYMER = 10
DEFAULT_MAX_REPEAT_FRACTION = 0.25
DEFAULT_VARIANT_EDGE_EXCLUSION = 20

# The two sites used for cloning; an insert containing either cannot be
# cut back out cleanly, so it is dropped.
RESTRICTION_SITES = {"SbfI": "CCTGCA^GG", "EcoRI": "G^AATTC"}

# IUPAC → regex, so a degenerate site spec still matches.
IUPAC = {"B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "K": "[GT]",
         "M": "[AC]", "N": ".", "R": "[AG]", "S": "[CG]",
         "V": "[ACG]", "W": "[AT]", "Y": "[CT]",
         "A": "A", "C": "C", "G": "G", "T": "T"}

STRAND_NAME = {"+": "fwd", "-": "rev", ".": "none"}

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"oligo_design_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


def _run_dir(label: str) -> Path:
    d = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ─── I/O primitives (bgzip is gzip-compatible, so gzip covers both) ─────────

def opener(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def is_gzip(path: str | Path) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(id, sequence)``. Handles multi-line records and .gz."""
    name, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(path: str | Path, records: Iterable[tuple[str, str]]) -> int:
    n = 0
    with opener(path, "wt") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
            n += 1
    return n


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


class BedRecord:
    """A BED line. Extra columns past strand are preserved verbatim."""

    __slots__ = ("chrom", "start", "end", "name", "score", "strand", "extra")

    def __init__(self, chrom, start, end, name=".", score=".", strand=".",
                 extra=()):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = score
        self.strand = strand
        self.extra = list(extra)

    @property
    def length(self) -> int:
        return self.end - self.start

    def fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), self.name,
                str(self.score), self.strand] + [str(x) for x in self.extra]

    def copy(self) -> "BedRecord":
        return BedRecord(self.chrom, self.start, self.end, self.name,
                         self.score, self.strand, list(self.extra))

    def __repr__(self) -> str:
        return f"<BedRecord {self.chrom}:{self.start}-{self.end} {self.name}>"


def read_bed(path: str | Path) -> list[BedRecord]:
    out: list[BedRecord] = []
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                continue
            out.append(BedRecord(
                f[0], f[1], f[2],
                f[3] if len(f) > 3 else ".",
                f[4] if len(f) > 4 else ".",
                f[5] if len(f) > 5 else ".",
                f[6:] if len(f) > 6 else (),
            ))
    return out


def write_bed(path: str | Path, records: Iterable[BedRecord]) -> int:
    n = 0
    with opener(path, "wt") as fh:
        for r in records:
            fh.write("\t".join(r.fields()) + "\n")
            n += 1
    return n


class IntervalIndex:
    """Per-chromosome sorted interval index — the tabix stand-in.

    Upstream fetches simple-repeat / TSS / CTCF annotations through
    ``pysam.TabixFile``, which needs the compiled htslib stack and a .tbi
    beside every BED. The annotation files here are small enough to hold
    in memory, so we sort by start and binary-search instead: same
    answers, no compiled dependency, no index files to keep in sync.
    """

    def __init__(self, path: str | Path | None = None):
        self.by_chrom: dict[str, tuple[list[int], list[tuple[int, int]]]] = {}
        self.path = str(path) if path else None
        self.n = 0
        if path:
            self._load(path)

    def _load(self, path: str | Path) -> None:
        raw: dict[str, list[tuple[int, int]]] = {}
        with opener(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\r\n").split("\t")
                if len(f) < 3:
                    continue
                try:
                    raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
                except ValueError:
                    continue
        for chrom, ivs in raw.items():
            ivs.sort()
            self.by_chrom[chrom] = ([s for s, _ in ivs], ivs)
            self.n += len(ivs)

    def fetch(self, chrom: str, start: int, end: int) -> list[tuple[int, int]]:
        """Every indexed interval overlapping ``[start, end)``."""
        entry = self.by_chrom.get(chrom)
        if entry is None:
            # Tolerate chr-prefix mismatch between design and annotation.
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            entry = self.by_chrom.get(alt)
            if entry is None:
                return []
        starts, ivs = entry
        # Any overlap must start before `end`; walk back from there while
        # intervals can still reach `start`.
        i = bisect.bisect_left(starts, end)
        hits = []
        for j in range(i - 1, -1, -1):
            s, e = ivs[j]
            if e > start:
                hits.append((s, e))
            # Annotation intervals are short; stop once clearly past.
            if s < start - 100_000:
                break
        return hits

    def __bool__(self) -> bool:
        return bool(self.by_chrom)


class Reference:
    """Random-access FASTA reader using a ``.fai`` when one exists.

    Replaces ``pyfaidx.Fasta``. With an index we seek; without one we
    build an equivalent index in memory on first use.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            raise SystemExit(f"Reference FASTA not found: {self.path}")
        if str(self.path).endswith(".gz"):
            raise SystemExit(
                f"{self.path} is gzip-compressed. Random access needs a "
                "plain (or bgzip+faidx) FASTA — decompress it first.")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        fai = Path(str(self.path) + ".fai")
        if fai.exists():
            self._read_fai(fai)
        else:
            self._build_index()
        self._fh = open(self.path, "rb")

    def _read_fai(self, fai: Path) -> None:
        for line in fai.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 5:
                self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))

    def _build_index(self) -> None:
        name = None
        length = offset = line_bases = line_width = 0
        with open(self.path, "rb") as fh:
            pos = 0
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        self.index[name] = (length, offset, line_bases, line_width)
                    name = raw[1:].split()[0].decode()
                    length, offset = 0, pos + len(raw)
                    line_bases = line_width = 0
                else:
                    if line_bases == 0:
                        line_width = len(raw)
                        line_bases = len(raw.rstrip(b"\r\n"))
                    length += len(raw.rstrip(b"\r\n"))
                pos += len(raw)
        if name is not None:
            self.index[name] = (length, offset, line_bases, line_width)

    def chrom_length(self, chrom: str) -> int | None:
        entry = self.index.get(chrom) or self.index.get(self._alt(chrom))
        return entry[0] if entry else None

    @staticmethod
    def _alt(chrom: str) -> str:
        return chrom[3:] if chrom.startswith("chr") else "chr" + chrom

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based, half-open. Out-of-range requests are clipped."""
        entry = self.index.get(chrom)
        if entry is None:
            entry = self.index.get(self._alt(chrom))
        if entry is None:
            raise SystemExit(
                f"Contig {chrom!r} is absent from {self.path.name}. "
                f"Known contigs: {', '.join(list(self.index)[:5])}...")
        length, offset, line_bases, line_width = entry
        start = max(0, start)
        end = min(end, length)
        if end <= start:
            return ""
        newlines = line_width - line_bases
        begin = offset + start // line_bases * line_width + start % line_bases
        stop = offset + end // line_bases * line_width + end % line_bases
        self._fh.seek(begin)
        raw = self._fh.read(stop - begin)
        return raw.decode().replace("\n", "").replace("\r", "").upper()

    def close(self) -> None:
        try:
            self._fh.close()
        except Exception:
            pass


# ─── Tiling ─────────────────────────────────────────────────────────────────
#
# Three strategies, chosen by region length (upstream
# rules/tiling.smk::tiling_splitUpForStrategies):
#
#   len <= centering_max   → one oligo, region centred and padded
#   len <= two_tiles_max   → two oligos pinned to the region's ends
#   otherwise              → centred tiling outwards with min_overlap
#
# The centred generators below are a faithful port of
# workflow/scripts/tiling/centerTiling.py: the arithmetic decides where
# every tile boundary lands, so it is reproduced exactly rather than
# rewritten. Changing it would silently change existing designs.

def tile_region(start: int, end: int, oligo_length: int, min_overlap: int,
                is_left: bool = True) -> Iterator[tuple[int, int]]:
    """One flank of a centred tiling. Port of ``centerTiling.tileRegion``."""
    span = end - start
    div = (span - min_overlap) // (oligo_length - min_overlap)
    if (span - min_overlap) % (oligo_length - min_overlap) != 0:
        div += 1

    if div == 2:                       # just centre + this one flank
        if is_left:
            yield (start, start + oligo_length)
        else:
            yield (end - oligo_length, end)
        return

    if div < 2:
        return

    step = oligo_length - min_overlap
    extension = ((step - (span - min_overlap) % step) % step) // (div - 1)
    first = ((step - (span - min_overlap) % step) % step) % (div - 1)
    for i in range(div):
        tile_left = start
        tile_right = start + oligo_length
        # The centre tile is emitted by the caller; skip it here.
        if not ((i == div - 1 and is_left) or (i == 0 and not is_left)):
            yield (tile_left, tile_right)
        if (i == div - 2 and is_left) or (i == 0 and not is_left):
            start = tile_right - (extension + min_overlap + first)
        else:
            start = tile_right - (extension + min_overlap)


def tile_center_region(start: int, end: int, oligo_length: int,
                       min_overlap: int) -> Iterator[tuple[int, int]]:
    """Centred tile first, then both flanks. Port of ``tileCenterRegion``."""
    center = start + (end - start) / 2
    center_left = int(center - (oligo_length / 2))
    center_right = int(center + (oligo_length / 2))
    yield (center_left, center_right)
    yield from tile_region(start, center_right, oligo_length, min_overlap,
                           is_left=True)
    yield from tile_region(center_left, end, oligo_length, min_overlap,
                           is_left=False)


def _tile_name(base: str, strand: str, i: int, n: int) -> str:
    return f"{base}_{STRAND_NAME.get(strand, 'none')}_tile{i}-{n}"


def tiles_for_centered(rec: BedRecord, oligo_length: int,
                       min_overlap: int) -> list[BedRecord]:
    tiles = list(tile_center_region(rec.start, rec.end, oligo_length,
                                    min_overlap))
    n = len(tiles)
    half = n // 2
    out = []
    for i, (s, e) in enumerate(tiles):
        t = rec.copy()
        t.start, t.end = s, e
        # Upstream numbering: the centre tile (emitted first) takes the
        # middle index so tile numbers read left→right across the region.
        if i == 0:
            idx = half + 1
        elif i > half:
            idx = i + 1
        else:
            idx = i
        t.name = _tile_name(rec.name, rec.strand, idx, n)
        out.append(t)
    return out


def tiles_for_two(rec: BedRecord, oligo_length: int,
                  extension: int = 0) -> list[BedRecord]:
    """Two oligos pinned to the region's ends (optionally slopped)."""
    start = rec.start - extension
    end = rec.end + extension
    left = rec.copy()
    left.start, left.end = start, start + oligo_length
    left.name = _tile_name(rec.name, rec.strand, 1, 2)
    right = rec.copy()
    right.start, right.end = end - oligo_length, end
    right.name = _tile_name(rec.name, rec.strand, 2, 2)
    return [left, right]


def tiles_for_none(rec: BedRecord, oligo_length: int) -> list[BedRecord]:
    """Region shorter than an oligo: centre it and pad to oligo_length."""
    missing = (oligo_length - rec.length) / 2
    t = rec.copy()
    # ceil on the left, floor on the right — matches the upstream awk.
    t.start = rec.start - int(-(-missing // 1))
    t.end = rec.end + int(missing // 1)
    t.name = _tile_name(rec.name, rec.strand, 1, 1)
    return [t]


def two_tiles_max(oligo_length: int, min_overlap: int,
                  configured_max: int, variant_edge_exclusion: int,
                  include_variant_edge: bool) -> int:
    """Port of ``tiling_common.smk::tiling_getMinTwoTiles``."""
    edge = variant_edge_exclusion * 2 if include_variant_edge else 0
    return min(configured_max,
               2 * oligo_length - 2 * min_overlap - edge)


def tile_regions(regions: list[BedRecord], *, oligo_length: int,
                 min_overlap: int, centering_max: int,
                 two_tiles_max_len: int,
                 two_tiles_extension: int = 0) -> tuple[list[BedRecord], dict]:
    out: list[BedRecord] = []
    counts = {"no_tiles": 0, "two_tiles": 0, "centered": 0}
    for rec in regions:
        if rec.length <= centering_max:
            out.extend(tiles_for_none(rec, oligo_length))
            counts["no_tiles"] += 1
        elif rec.length <= two_tiles_max_len:
            out.extend(tiles_for_two(rec, oligo_length, two_tiles_extension))
            counts["two_tiles"] += 1
        else:
            out.extend(tiles_for_centered(rec, oligo_length, min_overlap))
            counts["centered"] += 1
    out.sort(key=lambda r: (r.chrom, r.start))
    return out, counts


# ─── Minimal VCF ────────────────────────────────────────────────────────────

class VcfRecord:
    __slots__ = ("chrom", "pos", "id", "ref", "alt", "qual", "filter",
                 "info", "rest")

    def __init__(self, fields: list[str]):
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4].split(",")
        self.qual = fields[5] if len(fields) > 5 else "."
        self.filter = fields[6] if len(fields) > 6 else "."
        self.info = fields[7] if len(fields) > 7 else "."
        self.rest = fields[8:]

    @property
    def variant_id(self) -> str:
        if self.id and self.id != ".":
            return self.id
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt[0]}"

    def alt_type(self, alt: str) -> str:
        if len(alt) == len(self.ref):
            return "SNV" if len(alt) == 1 else "MNV"
        return "INDEL"

    def line(self, extra_info: dict[str, list[str]] | None = None) -> str:
        info = self.info
        if extra_info:
            add = ";".join(f"{k}={','.join(v)}" for k, v in extra_info.items() if v)
            if add:
                info = add if info in (".", "") else f"{info};{add}"
        return "\t".join([self.chrom, str(self.pos), self.id, self.ref,
                          ",".join(self.alt), self.qual, self.filter, info]
                         + self.rest)


def read_vcf(path: str | Path) -> tuple[list[str], list[VcfRecord]]:
    header, records = [], []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                header.append(line)
                continue
            if not line.strip():
                continue
            records.append(VcfRecord(line.split("\t")))
    return header, records


INFO_LINES = [
    '##INFO=<ID=Region,Number=.,Type=String,Description="Matched regions">',
    '##INFO=<ID=ALT_ID,Number=.,Type=String,Description="Corresponding alt ID(s)">',
    '##INFO=<ID=REF_ID,Number=.,Type=String,Description="Corresponding REF ID(s)">',
]


def _augmented_header(header: list[str]) -> list[str]:
    if not header:
        return ["##fileformat=VCFv4.2"] + INFO_LINES + [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    out = list(header)
    col = out.pop() if out[-1].startswith("#CHROM") else None
    for line in INFO_LINES:
        if line not in out:
            out.append(line)
    if col:
        out.append(col)
    return out


# ─── Design: variants ───────────────────────────────────────────────────────

def _centre_distance(pos0: int, rec: BedRecord, legacy: bool = False) -> float:
    """Distance from a 0-based variant position to a region's centre.

    ``legacy=True`` reproduces upstream's ``abs(pos - start - length)``,
    which measures distance to the region *end*. See the module
    docstring: the corrected metric is the default.
    """
    if legacy:
        return abs(pos0 - rec.start - rec.length)
    return abs(pos0 - (rec.start + rec.length / 2))


def build_alt_sequence(reference: Reference, rec: BedRecord,
                       var: VcfRecord, alt: str) -> tuple[str, str]:
    """Return ``(ref_seq, alt_seq)`` for one region × one ALT allele.

    Port of ``getSequencesInclVariants.getSequences``. For a deletion the
    REF window is *extended* by the length difference so both oligos come
    out the same length — that is what makes REF and ALT comparable on
    the array.
    """
    ref_seq = reference.fetch(rec.chrom, rec.start, rec.end)
    variant_position = var.pos - 1 - rec.start
    if variant_position < 0 or variant_position >= len(ref_seq):
        raise ValueError(
            f"variant {var.variant_id} at {var.chrom}:{var.pos} falls outside "
            f"region {rec.name} ({rec.chrom}:{rec.start}-{rec.end})")

    kind = var.alt_type(alt)
    end = -(len(alt) + 1) if len(alt) > 1 else None

    if kind in ("SNV", "MNV"):
        alt_seq = ref_seq[:variant_position] + alt + ref_seq[variant_position + 1:end]
    else:  # INDEL / INS / DEL
        if len(alt) > len(var.ref):
            alt_seq = ref_seq[:variant_position] + alt + ref_seq[variant_position + 1:end]
        else:
            alt_seq = ref_seq
            extended = reference.fetch(
                rec.chrom, rec.start, rec.end + len(var.ref) - len(alt))
            ref_seq = (extended[:variant_position] + alt
                       + extended[variant_position + len(var.ref):])

    if rec.strand == "-":
        alt_seq = reverse_complement(alt_seq)
        ref_seq = reverse_complement(ref_seq)
    return ref_seq, alt_seq


def design_variants(regions: list[BedRecord], variants: list[VcfRecord],
                    reference: Reference, *,
                    variant_edge_exclusion: int = 0,
                    use_most_centered_region: bool = True,
                    remove_regions_without_variants: bool = False,
                    legacy_centre_metric: bool = False) -> dict[str, Any]:
    """Pair every variant with its region(s) and emit REF/ALT oligos."""
    by_chrom: dict[str, list[BedRecord]] = {}
    for r in regions:
        by_chrom.setdefault(r.chrom, []).append(r)
    for v in by_chrom.values():
        v.sort(key=lambda r: r.start)

    sequences: list[tuple[str, str]] = []      # regions with no variant
    ref_sequences: dict[str, str] = {}
    alt_sequences: dict[str, str] = {}
    design_map: list[dict[str, str]] = []
    kept_regions: dict[str, BedRecord] = {}
    kept_variant_lines: list[str] = []
    removed_variants: list[VcfRecord] = []

    for var in variants:
        pos0 = var.pos - 1
        # The exclusion zone keeps a variant off the very edge of an
        # oligo, where flanking context is truncated.
        candidates = [
            r for r in by_chrom.get(var.chrom, [])
            if (r.start + variant_edge_exclusion) <= pos0
            < (r.end - variant_edge_exclusion)
        ]
        if not candidates:
            removed_variants.append(var)
            continue

        if use_most_centered_region and len(candidates) > 1:
            best = min(_centre_distance(pos0, r, legacy_centre_metric)
                       for r in candidates)
            candidates = [r for r in candidates
                          if _centre_distance(pos0, r, legacy_centre_metric) == best][:1]

        if len(var.alt) > 1:
            raise SystemExit(
                f"Only one ALT allele is supported per record; "
                f"{var.variant_id} has {len(var.alt)}. Split the VCF first.")
        alt = var.alt[0]

        region_names, ref_ids, alt_ids = [], [], []
        for rec in candidates:
            try:
                ref_seq, alt_seq = build_alt_sequence(reference, rec, var, alt)
            except ValueError as e:
                logging.warning("skipped: %s", e)
                continue
            ref_id = f"REF_{rec.name}"
            alt_id = f"ALT_{rec.name}_{var.variant_id}"
            ref_sequences[ref_id] = ref_seq
            alt_sequences[alt_id] = alt_seq
            design_map.append({"Variant": var.variant_id, "Region": rec.name,
                               "REF_ID": ref_id, "ALT_ID": alt_id})
            kept_regions[rec.name] = rec
            region_names.append(rec.name)
            ref_ids.append(ref_id)
            alt_ids.append(alt_id)

        if not region_names:
            removed_variants.append(var)
            continue
        kept_variant_lines.append(var.line(
            {"Region": region_names, "REF_ID": ref_ids, "ALT_ID": alt_ids}))

    region_map = [{"Region": d["Region"], "ID": d[k]}
                  for d in design_map for k in ("REF_ID", "ALT_ID")]

    if not remove_regions_without_variants:
        # Regions no variant landed in still go on the array, as plain
        # reference sequences — they are the design's internal baseline.
        used = {d["Region"] for d in design_map}
        for rec in regions:
            if rec.name in used:
                continue
            seq = reference.fetch(rec.chrom, rec.start, rec.end)
            if rec.strand == "-":
                seq = reverse_complement(seq)
            sequences.append((rec.name, seq))
            region_map.append({"Region": rec.name, "ID": rec.name})

    seen = set()
    region_map = [m for m in region_map
                  if not ((m["Region"], m["ID"]) in seen
                          or seen.add((m["Region"], m["ID"])))]
    seen = set()
    design_map = [d for d in design_map
                  if not (tuple(d.items()) in seen or seen.add(tuple(d.items())))]

    output_regions = (list(kept_regions.values())
                      if remove_regions_without_variants else list(regions))
    removed_regions = [r for r in regions if r.name not in kept_regions]

    return {
        "sequences": sequences,
        "ref_sequences": ref_sequences,
        "alt_sequences": alt_sequences,
        "design_map": design_map,
        "region_map": region_map,
        "regions": output_regions,
        "removed_regions": removed_regions,
        "variant_lines": kept_variant_lines,
        "removed_variants": removed_variants,
    }


def design_regions(regions: list[BedRecord], reference: Reference, *,
                   oligo_length: int | None = None) -> list[tuple[str, str]]:
    """Region BED → oligo FASTA (port of ``oligo_design_regions_getSequences``)."""
    bad = [r for r in regions if oligo_length and r.length != oligo_length]
    if bad:
        detail = ", ".join(f"{r.name} ({r.length})" for r in bad[:5])
        raise SystemExit(
            f"{len(bad)} region(s) are not exactly {oligo_length} bp: {detail}"
            + ("..." if len(bad) > 5 else "")
            + "\nRun `oligo tile` first, or pass --no-length-check.")
    out = []
    for rec in regions:
        seq = reference.fetch(rec.chrom, rec.start, rec.end)
        if rec.strand == "-":
            seq = reverse_complement(seq)
        out.append((rec.name, seq))
    return out


# ─── Filters ────────────────────────────────────────────────────────────────

def longest_nucleotide_run(seq: str) -> int:
    """Longest homopolymer run.

    Upstream's ``nucleotideruns`` never closes the final run, so a
    sequence ending in its longest run is under-reported ("AAAA" → 0).
    Fixed here; see the module docstring.
    """
    longest = 0
    last, run = None, 0
    for base in seq:
        if base == last:
            run += 1
        else:
            longest = max(longest, run)
            last, run = base, 1
    return max(longest, run)


def site_to_regex(site: str) -> str:
    return "".join(IUPAC[c] for c in site.upper() if c in IUPAC)


_RESTRICTION_RE = {name: re.compile(site_to_regex(site), re.IGNORECASE)
                   for name, site in RESTRICTION_SITES.items()}


def filter_sequences(records: Iterable[tuple[str, str]],
                     max_homopolymer: int | None) -> tuple[list[str], dict, list]:
    """Homopolymer + restriction-site filters. Sequence-level, so it also
    applies to a pure-FASTA design with no coordinates."""
    failed, detail = [], []
    reasons = {"hompol": 0, "restrictions": 0}
    for cid, seq in records:
        run = longest_nucleotide_run(seq)
        if max_homopolymer is not None and run > max_homopolymer:
            reasons["hompol"] += 1
            failed.append(cid)
            detail.append({"id": cid, "reason": "homopolymer",
                           "detail": f"run of {run} > {max_homopolymer}"})
            continue
        hit = next((name for name, rx in _RESTRICTION_RE.items()
                    if rx.search(seq)), None)
        if hit:
            reasons["restrictions"] += 1
            failed.append(cid)
            detail.append({"id": cid, "reason": "restriction_site",
                           "detail": f"{hit} site ({RESTRICTION_SITES[hit]})"})
    return failed, reasons, detail


def filter_regions(regions: list[BedRecord], *, repeats: IntervalIndex,
                   tss: IntervalIndex, ctcf: IntervalIndex,
                   max_repeat_fraction: float) -> tuple[list[str], dict, list]:
    """Simple-repeat / TSS / CTCF filters. Coordinate-level, so it only
    applies to designs that came from regions."""
    failed, detail = [], []
    reasons = {"repeats": 0, "TSS": 0, "CTCF": 0}
    for r in regions:
        why = []
        span = float(r.length) or 1.0
        for tstart, tend in repeats.fetch(r.chrom, r.start, r.end):
            overlap = min(tend, r.end) - max(tstart, r.start)
            if overlap / span > max_repeat_fraction:
                reasons["repeats"] += 1
                why.append(f"simple repeat covers {overlap / span:.0%}")
                break
        if tss.fetch(r.chrom, r.start, r.end):
            reasons["TSS"] += 1
            why.append("overlaps an annotated TSS")
        if ctcf.fetch(r.chrom, r.start, r.end):
            reasons["CTCF"] += 1
            why.append("overlaps a CTCF motif")
        if why:
            failed.append(r.name)
            detail.append({"id": r.name, "reason": "region",
                           "detail": "; ".join(why)})
    return failed, reasons, detail


def apply_filters(design_map: list[dict], region_map: list[dict], *,
                  failed_regions: list[str], failed_seqs: list[str],
                  remove_regions_without_variants: bool) -> tuple[list[dict], list[dict], int, int]:
    """Propagate failures through the maps (port of ``filterOligos.write_output``).

    The subtle rule: if a REF oligo fails, its ALT partners must go too —
    an ALT with no REF to compare against is not interpretable.
    """
    total = len({d[k] for d in design_map for k in ("REF_ID", "ALT_ID")}
                | {m["ID"] for m in region_map})

    fr, fs = set(failed_regions), set(failed_seqs)
    if fr:
        region_map = [m for m in region_map if m["Region"] not in fr]
        design_map = [d for d in design_map if d["Region"] not in fr]
    if fs:
        region_map = [m for m in region_map if m["ID"] not in fs]
        ref_to_alt: dict[str, list[str]] = {}
        for d in design_map:
            ref_to_alt.setdefault(d["REF_ID"], []).append(d["ALT_ID"])
        orphaned = {a for ref in fs for a in ref_to_alt.get(ref, [])}
        design_map = [d for d in design_map
                      if d["REF_ID"] not in fs and d["ALT_ID"] not in fs]
        region_map = [m for m in region_map if m["ID"] not in orphaned]
        if remove_regions_without_variants:
            keep = {d[k] for d in design_map for k in ("REF_ID", "ALT_ID")}
            region_map = [m for m in region_map if m["ID"] in keep]

    remaining = len({d[k] for d in design_map for k in ("REF_ID", "ALT_ID")}
                    | {m["ID"] for m in region_map})
    return design_map, region_map, total, total - remaining


def add_adapters(records: Iterable[tuple[str, str]], left: str,
                 right: str) -> Iterator[tuple[str, str]]:
    for cid, seq in records:
        yield cid, f"{left}{seq}{right}"


def kmer_filter(candidates: Iterable[tuple[str, str]],
                reference_seqs: dict[str, str], k: int = 6
                ) -> list[tuple[str, str, int]]:
    """Per group, keep the candidate sharing fewest k-mers with its source.

    Port of ``kMerFilter.py``. Used to pick a scrambled negative control
    that stays as far as possible from the real sequence it was derived
    from. Candidate ids are ``<source_id>_<n>``.
    """
    best: dict[str, tuple[str, int]] = {}
    for cid, seq in candidates:
        source = "_".join(cid.split("_")[:-1])
        if source not in reference_seqs:
            raise SystemExit(f"Sequence {source!r} missing from the input FASTA.")
        ref = reference_seqs[source]
        ref_kmers = {ref[i:i + k] for i in range(len(ref) - k + 1)}
        cand_kmers = {seq[i:i + k] for i in range(len(seq) - k + 1)}
        shared = len(cand_kmers & ref_kmers)
        if source not in best or shared < best[source][1]:
            best[source] = (seq, shared)
    return [(cid, seq, n) for cid, (seq, n) in best.items()]


# ─── Reference annotations ──────────────────────────────────────────────────
#
# The three coordinate filters need annotation BEDs. Upstream ships them
# in reference/; here they are looked up by convention and every one is
# optional — a missing file disables just its filter, loudly, instead of
# failing the run.

ANNOTATION_FILES = {
    "repeats": ("simpleRepeat.bed.gz", "simple repeats (UCSC simpleRepeat)"),
    "tss": ("TSS_pos.bed.gz", "transcription start sites (GENCODE)"),
    "ctcf": ("CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz",
             "CTCF motifs (JASPAR MA0139.1 footprints)"),
}


def load_annotation(kind: str, explicit: str | None) -> IntervalIndex:
    filename, what = ANNOTATION_FILES[kind]
    path = Path(explicit) if explicit else REFERENCE_DIR / filename
    if not path.exists():
        logging.warning(
            "%s filter DISABLED — %s not found at %s. Pass --%s to enable it.",
            kind, what, path, kind.replace("_", "-"))
        return IntervalIndex()
    idx = IntervalIndex(path)
    logging.info("%s: %d intervals from %s", kind, idx.n, path)
    return idx


def _write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    with opener(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_tile(args) -> int:
    setup_logging()
    regions = read_bed(args.regions)
    tt_max = two_tiles_max(args.oligo_length, args.min_overlap,
                           args.two_tiles_max, args.variant_edge_exclusion,
                           args.include_variant_edge)
    tiles, counts = tile_regions(
        regions, oligo_length=args.oligo_length, min_overlap=args.min_overlap,
        centering_max=args.centering_max, two_tiles_max_len=tt_max,
        two_tiles_extension=(args.variant_edge_exclusion
                             if args.include_variant_edge else 0))
    out_dir = _run_dir(args.label or "tile")
    out = out_dir / "regions.tiles.bed.gz"
    write_bed(out, tiles)
    print(f"Regions in:    {len(regions):,}")
    print(f"  no-tiling    {counts['no_tiles']:,}  (<= {args.centering_max} bp, centred + padded)")
    print(f"  two-tiles    {counts['two_tiles']:,}  (<= {tt_max} bp, pinned to ends)")
    print(f"  centred      {counts['centered']:,}  (> {tt_max} bp, tiled with >= {args.min_overlap} bp overlap)")
    print(f"Oligo windows: {len(tiles):,}")
    print(f"Saved:         {out}")
    return 0


def cmd_design_regions(args) -> int:
    setup_logging()
    regions = read_bed(args.regions)
    reference = Reference(args.reference)
    seqs = design_regions(
        regions, reference,
        oligo_length=None if args.no_length_check else args.oligo_length)
    out_dir = _run_dir(args.label or "design_regions")
    fa = out_dir / "design.fa"
    write_fasta(fa, seqs)
    _write_tsv(out_dir / "region_map.tsv.gz",
               [{"Region": r.name, "ID": r.name} for r in regions],
               ["Region", "ID"])
    reference.close()
    print(f"Regions:  {len(regions):,}")
    print(f"Oligos:   {len(seqs):,}")
    print(f"Saved:    {fa}")
    return 0


def cmd_design_variants(args) -> int:
    setup_logging()
    regions = read_bed(args.regions)
    header, variants = read_vcf(args.variants)
    reference = Reference(args.reference)
    res = design_variants(
        regions, variants, reference,
        variant_edge_exclusion=args.variant_edge_exclusion,
        use_most_centered_region=not args.use_all_regions,
        remove_regions_without_variants=args.remove_regions_without_variants,
        legacy_centre_metric=args.legacy_centre_metric)
    out_dir = _run_dir(args.label or "design_variants")

    fa = out_dir / "design.fa"
    write_fasta(fa, list(res["sequences"])
                + sorted(res["ref_sequences"].items())
                + sorted(res["alt_sequences"].items()))
    _write_tsv(out_dir / "variant_region_map.tsv.gz", res["design_map"],
               ["Variant", "Region", "REF_ID", "ALT_ID"])
    _write_tsv(out_dir / "region_map.tsv.gz", res["region_map"],
               ["Region", "ID"])
    write_bed(out_dir / "regions.bed.gz", res["regions"])
    write_bed(out_dir / "regions.removed.bed.gz", res["removed_regions"])
    with opener(out_dir / "variants.vcf.gz", "wt") as fh:
        fh.write("\n".join(_augmented_header(header)) + "\n")
        fh.write("\n".join(res["variant_lines"]) + "\n")
    with opener(out_dir / "variants.removed.vcf.gz", "wt") as fh:
        fh.write("\n".join(header) + "\n")
        for v in res["removed_variants"]:
            fh.write(v.line() + "\n")
    reference.close()

    n_pairs = len(res["design_map"])
    print(f"Variants in:      {len(variants):,}")
    print(f"  paired:         {len(variants) - len(res['removed_variants']):,}")
    print(f"  unplaced:       {len(res['removed_variants']):,}  (no region, or too close to an oligo edge)")
    print(f"REF/ALT pairs:    {n_pairs:,}")
    print(f"Oligos:           {len(res['ref_sequences']) + len(res['alt_sequences']) + len(res['sequences']):,}")
    print(f"Saved:            {out_dir}")
    return 0


def cmd_filter(args) -> int:
    setup_logging()
    seqs = list(read_fasta(args.design))
    failed_seqs, seq_reasons, detail = filter_sequences(
        seqs, args.max_homopolymer_length)

    failed_regions: list[str] = []
    region_reasons = {"repeats": 0, "TSS": 0, "CTCF": 0}
    if args.regions:
        regions = read_bed(args.regions)
        failed_regions, region_reasons, rdetail = filter_regions(
            regions,
            repeats=load_annotation("repeats", args.simple_repeats),
            tss=load_annotation("tss", args.tss_positions),
            ctcf=load_annotation("ctcf", args.ctcf_motifs),
            max_repeat_fraction=args.max_simple_repeat_fraction)
        detail += rdetail

    design_map = _read_tsv(args.variant_map) if args.variant_map else []
    region_map = _read_tsv(args.map) if args.map else [
        {"Region": cid, "ID": cid} for cid, _ in seqs]
    design_map, region_map, total, removed = apply_filters(
        design_map, region_map,
        failed_regions=failed_regions, failed_seqs=failed_seqs,
        remove_regions_without_variants=args.remove_regions_without_variants)

    out_dir = _run_dir(args.label or "filter")
    keep = {m["ID"] for m in region_map} | {
        d[k] for d in design_map for k in ("REF_ID", "ALT_ID")}
    kept = [(cid, s) for cid, s in seqs if cid in keep]
    write_fasta(out_dir / "design_filtered.fa", kept)
    _write_tsv(out_dir / "region_map.tsv.gz", region_map, ["Region", "ID"])
    if design_map:
        _write_tsv(out_dir / "variant_region_map.tsv.gz", design_map,
                   ["Variant", "Region", "REF_ID", "ALT_ID"])
    _write_tsv(out_dir / "filter.log.tsv", detail, ["id", "reason", "detail"])

    print("Failed sequences:")
    print(f"  {seq_reasons['hompol']:>7,} homopolymer runs > {args.max_homopolymer_length}")
    print(f"  {region_reasons['repeats']:>7,} simple repeats > {args.max_simple_repeat_fraction:.0%} of the window")
    print(f"  {region_reasons['TSS']:>7,} TSS overlap")
    print(f"  {region_reasons['CTCF']:>7,} CTCF motif overlap")
    print(f"  {seq_reasons['restrictions']:>7,} EcoRI / SbfI restriction site")
    print(f"Total failed:  {removed:,}")
    print(f"Total passed:  {total - removed:,}")
    print(f"Saved:         {out_dir}")
    return 0


def _read_tsv(path: str | Path) -> list[dict]:
    rows = []
    with opener(path) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\r\n").split("\t"))))
    return rows


def cmd_adapters(args) -> int:
    setup_logging()
    seqs = list(read_fasta(args.design))
    out_dir = _run_dir(args.label or "adapters")
    out = out_dir / "design.adapters.fa"
    n = write_fasta(out, add_adapters(seqs, args.left, args.right))
    total = len(args.left) + len(args.right)
    print(f"Oligos:        {n:,}")
    print(f"Adapters:      left {len(args.left)} bp + right {len(args.right)} bp = {total} bp")
    if seqs:
        print(f"Final length:  {len(seqs[0][1]) + total} bp (insert {len(seqs[0][1])} bp)")
    print(f"Saved:         {out}")
    return 0


def cmd_kmer_filter(args) -> int:
    setup_logging()
    reference_seqs = dict(read_fasta(args.reference_fasta))
    candidates = read_fasta(args.candidates)
    picked = kmer_filter(candidates, reference_seqs, k=args.kmers)
    out_dir = _run_dir(args.label or "kmer_filter")
    out = out_dir / "selected.fa"
    with opener(out, "wt") as fh:
        for cid, seq, n in sorted(picked):
            fh.write(f">{cid} {n}\n{seq}\n" if args.include_overlap
                     else f">{cid}\n{seq}\n")
    print(f"Groups:   {len(picked):,}")
    if picked:
        shared = [n for _, _, n in picked]
        print(f"Shared {args.kmers}-mers with source: "
              f"min {min(shared)}, median {sorted(shared)[len(shared) // 2]}, max {max(shared)}")
    print(f"Saved:    {out}")
    return 0


def cmd_pipeline(args) -> int:
    """tile → design → filter → adapters → final design, in one run dir."""
    setup_logging()
    out_dir = _run_dir(args.label or "pipeline")
    reference = Reference(args.reference)
    regions = read_bed(args.regions)
    summary: dict[str, Any] = {"sample": args.label or "design",
                               "oligo_length": args.oligo_length,
                               "input_regions": len(regions)}

    # 1. tiling
    if args.tile:
        tt_max = two_tiles_max(args.oligo_length, args.min_overlap,
                               args.two_tiles_max, args.variant_edge_exclusion,
                               args.include_variant_edge)
        regions, counts = tile_regions(
            regions, oligo_length=args.oligo_length,
            min_overlap=args.min_overlap, centering_max=args.centering_max,
            two_tiles_max_len=tt_max,
            two_tiles_extension=(args.variant_edge_exclusion
                                 if args.include_variant_edge else 0))
        summary["tiling"] = counts
        write_bed(out_dir / "regions.tiles.bed.gz", regions)
    summary["oligo_windows"] = len(regions)

    # 2. design
    if args.variants:
        header, variants = read_vcf(args.variants)
        res = design_variants(
            regions, variants, reference,
            variant_edge_exclusion=args.variant_edge_exclusion,
            use_most_centered_region=not args.use_all_regions,
            remove_regions_without_variants=args.remove_regions_without_variants,
            legacy_centre_metric=args.legacy_centre_metric)
        seqs = (list(res["sequences"]) + sorted(res["ref_sequences"].items())
                + sorted(res["alt_sequences"].items()))
        design_map, region_map = res["design_map"], res["region_map"]
        design_regions_list = res["regions"]
        summary["variants_in"] = len(variants)
        summary["variants_unplaced"] = len(res["removed_variants"])
        summary["ref_alt_pairs"] = len(design_map)
    else:
        seqs = design_regions(regions, reference,
                              oligo_length=None if args.no_length_check
                              else args.oligo_length)
        design_map = []
        region_map = [{"Region": r.name, "ID": r.name} for r in regions]
        design_regions_list = regions
    summary["oligos_designed"] = len(seqs)

    # 3. filter
    failed_seqs, seq_reasons, detail = filter_sequences(
        seqs, args.max_homopolymer_length)
    failed_regions, region_reasons, rdetail = filter_regions(
        design_regions_list,
        repeats=load_annotation("repeats", args.simple_repeats),
        tss=load_annotation("tss", args.tss_positions),
        ctcf=load_annotation("ctcf", args.ctcf_motifs),
        max_repeat_fraction=args.max_simple_repeat_fraction)
    detail += rdetail
    design_map, region_map, total, removed = apply_filters(
        design_map, region_map, failed_regions=failed_regions,
        failed_seqs=failed_seqs,
        remove_regions_without_variants=args.remove_regions_without_variants)
    keep = {m["ID"] for m in region_map} | {
        d[k] for d in design_map for k in ("REF_ID", "ALT_ID")}
    seqs = [(cid, s) for cid, s in seqs if cid in keep]
    summary["filtered_out"] = removed
    summary["filter_reasons"] = {**seq_reasons, **region_reasons}
    _write_tsv(out_dir / "filter.log.tsv", detail, ["id", "reason", "detail"])

    # 4. adapters + final design, ids namespaced by sample as upstream does
    sample = safe_label(args.label or "design")
    final = [(f"{sample}:{cid}", s) for cid, s in
             add_adapters(seqs, args.left_adapter, args.right_adapter)]
    write_fasta(out_dir / "design.fa.gz", final)
    _write_tsv(out_dir / "region_map.tsv.gz", region_map, ["Region", "ID"])
    if design_map:
        _write_tsv(out_dir / "variant_region_map.tsv.gz", design_map,
                   ["Variant", "Region", "REF_ID", "ALT_ID"])
    summary["oligos_final"] = len(final)
    summary["adapter_bp"] = len(args.left_adapter) + len(args.right_adapter)
    reference.close()

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Input regions:   {summary['input_regions']:,}")
    print(f"Oligo windows:   {summary['oligo_windows']:,}")
    if args.variants:
        print(f"Variants:        {summary['variants_in']:,} "
              f"({summary['variants_unplaced']:,} unplaced)")
        print(f"REF/ALT pairs:   {summary['ref_alt_pairs']:,}")
    print(f"Designed:        {summary['oligos_designed']:,}")
    print(f"Filtered out:    {removed:,}  "
          + ", ".join(f"{k}={v}" for k, v in summary["filter_reasons"].items() if v))
    print(f"FINAL OLIGOS:    {len(final):,}  "
          f"({len(final[0][1]) if final else 0} bp incl. adapters)")
    print(f"Saved:           {out_dir}")
    return 0


PLAYBOOK = """# MPRA Oligo Design — skill playbook

Python port of [MPRAOligoDesign](https://github.com/kircherlab/MPRAOligoDesign)
(Max Schubach, Berlin Institute of Health / kircherlab; MIT).
Cite Rosen *et al.*, *Genome Research* (2025), doi:10.1101/gr.281462.125.

The upstream Snakemake workflow needs conda, bedtools, pysam, pyranges,
vcfpy, pyfaidx and BioPython. This port is standard-library only, so a
design runs wherever IGVFagent runs.

## The four decisions in a design

| Stage | Question | Subcommand |
|---|---|---|
| Tiling | Where does the oligo window sit? | `tile` |
| Sequence | What do we write on the array? | `design-regions` / `design-variants` |
| Filtering | What breaks synthesis or the readout? | `filter` |
| Adapters | What constant flanks get bolted on? | `adapters` |

## Tiling strategies

Chosen automatically by region length:

| Region length | Strategy | Result |
|---|---|---|
| `<= --centering-max` | no tiling | one oligo, region centred and padded |
| `<= two_tiles_max` | two tiles | two oligos pinned to the region's ends |
| longer | centred tiling | centre oligo, then outwards with `>= --min-overlap` |

`two_tiles_max = min(--two-tiles-max, 2*oligo_length - 2*min_overlap - edge)`,
where `edge = 2*--variant-edge-exclusion` when `--include-variant-edge` is set.

## Filters

| Filter | Applies to | Default |
|---|---|---|
| Homopolymer run | every sequence | `> 10` fails |
| EcoRI / SbfI site | every sequence | any occurrence fails |
| Simple repeat | region-derived only | `> 25%` of the window fails |
| TSS overlap | region-derived only | any overlap fails |
| CTCF motif | region-derived only | any overlap fails |

The last three need annotation BEDs under `Data/Reference/OligoDesign/`
(`simpleRepeat.bed.gz`, `TSS_pos.bed.gz`,
`CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz`) or explicit `--simple-repeats`
/ `--tss-positions` / `--ctcf-motifs` paths. A missing file disables that
filter with a warning rather than failing the run.

If a REF oligo fails, its ALT partners are dropped with it — an ALT with
no REF to compare against is not interpretable.

## Usage

```bash
# Tile regions into oligo windows
igvfagent oligo tile --regions regions.bed.gz --oligo-length 200 \\
    --min-overlap 50 --centering-max 200

# Region-only design
igvfagent oligo design-regions --regions regions.tiles.bed.gz \\
    --reference hg38.fa --oligo-length 200

# Variant design: one REF + one ALT oligo per variant x region
igvfagent oligo design-variants --regions regions.bed.gz \\
    --variants variants.vcf.gz --reference hg38.fa \\
    --variant-edge-exclusion 20

# Filter, then add adapters
igvfagent oligo filter --design design.fa --regions regions.bed.gz \\
    --map region_map.tsv.gz --max-homopolymer-length 10
igvfagent oligo adapters --design design_filtered.fa \\
    --left AGGACCGGATCAACT --right CATTGCGTGAACCGA

# Everything in one run directory
igvfagent oligo pipeline --regions regions.bed.gz --variants variants.vcf.gz \\
    --reference hg38.fa --oligo-length 200 --tile --label my_library
```

## Deviations from upstream

* **Homopolymer bug fixed.** Upstream's `nucleotideruns` never closes the
  final run, so `"ACCCC"` scores 1 and `"AAAA"` scores 0 — homopolymers at
  a sequence's end slip through. Fixed here.
* **Most-centred region fixed.** Upstream computes
  `abs(pos - 1 - start - length)`, the distance to the region *end*. The
  corrected centre distance is the default; `--legacy-centre-metric`
  restores upstream behaviour for bit-identical reproduction.
* No Snakemake / conda / bedtools / pysam / pyranges / vcfpy / pyfaidx /
  BioPython. Tabix lookups are replaced by an in-memory interval index.

## Outputs

Everything lands in `Docs/MPRAOligoDesign/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `design.fa.gz` | final oligos, ids namespaced `<sample>:<id>` |
| `regions.tiles.bed.gz` | oligo windows after tiling |
| `variant_region_map.tsv.gz` | Variant → Region → REF_ID → ALT_ID |
| `region_map.tsv.gz` | Region → sequence ID |
| `variants.vcf.gz` | placed variants, INFO annotated with Region/REF_ID/ALT_ID |
| `variants.removed.vcf.gz` | variants no oligo could carry |
| `filter.log.tsv` | one row per dropped sequence, with the reason |
| `summary.json` | counts for every stage |
"""


def cmd_write_playbook(args) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = SKILL_DOC_DIR / "MPRA_OLIGO_DESIGN_SKILLS.md"
    out.write_text(PLAYBOOK)
    print(f"Wrote {out}")
    return 0


def _add_tiling_args(p) -> None:
    p.add_argument("--oligo-length", type=int, default=DEFAULT_OLIGO_LENGTH)
    p.add_argument("--min-overlap", type=int, default=DEFAULT_MIN_OVERLAP,
                   help="Minimum overlap between adjacent tiles.")
    p.add_argument("--centering-max", type=int, default=DEFAULT_OLIGO_LENGTH,
                   help="Regions <= this length get a single centred oligo.")
    p.add_argument("--two-tiles-max", type=int, default=10_000,
                   help="Upper bound for the two-tile strategy (capped by "
                        "2*oligo_length - 2*min_overlap - edge).")
    p.add_argument("--variant-edge-exclusion", type=int,
                   default=DEFAULT_VARIANT_EDGE_EXCLUSION,
                   help="Keep variants this far from an oligo edge.")
    p.add_argument("--include-variant-edge", action="store_true",
                   help="Slop two-tile regions by the edge exclusion.")


def _add_filter_args(p) -> None:
    p.add_argument("--max-homopolymer-length", type=int,
                   default=DEFAULT_MAX_HOMOPOLYMER)
    p.add_argument("--max-simple-repeat-fraction", type=float,
                   default=DEFAULT_MAX_REPEAT_FRACTION)
    p.add_argument("--simple-repeats", default=None)
    p.add_argument("--tss-positions", default=None)
    p.add_argument("--ctcf-motifs", default=None)


def _add_variant_args(p) -> None:
    p.add_argument("--use-all-regions", action="store_true",
                   help="Emit an oligo for every region a variant hits "
                        "(default: only the most centred one).")
    p.add_argument("--remove-regions-without-variants", action="store_true",
                   help="Drop regions no variant landed in (default: keep "
                        "them as plain reference oligos).")
    p.add_argument("--legacy-centre-metric", action="store_true",
                   help="Reproduce upstream's distance-to-end metric when "
                        "picking the most centred region.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="igvfagent oligo",
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("tile", help="Region BED -> tiled oligo windows.")
    p.add_argument("--regions", required=True)
    p.add_argument("--label", default=None)
    _add_tiling_args(p)
    p.set_defaults(func=cmd_tile)

    p = sub.add_parser("design-regions",
                       help="Region BED + reference -> oligo FASTA.")
    p.add_argument("--regions", required=True)
    p.add_argument("--reference", required=True, help="Genome FASTA (uncompressed).")
    p.add_argument("--oligo-length", type=int, default=DEFAULT_OLIGO_LENGTH)
    p.add_argument("--no-length-check", action="store_true",
                   help="Allow regions that are not exactly oligo_length.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_design_regions)

    p = sub.add_parser("design-variants",
                       help="VCF + regions + reference -> REF/ALT oligos.")
    p.add_argument("--regions", required=True)
    p.add_argument("--variants", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--variant-edge-exclusion", type=int,
                   default=DEFAULT_VARIANT_EDGE_EXCLUSION)
    p.add_argument("--label", default=None)
    _add_variant_args(p)
    p.set_defaults(func=cmd_design_variants)

    p = sub.add_parser("filter", help="Apply the five design filters.")
    p.add_argument("--design", required=True, help="Designed oligo FASTA.")
    p.add_argument("--regions", default=None,
                   help="Region BED — enables the repeat/TSS/CTCF filters.")
    p.add_argument("--map", default=None, help="Region -> ID map TSV.")
    p.add_argument("--variant-map", default=None,
                   help="Variant -> Region -> REF/ALT map TSV.")
    p.add_argument("--remove-regions-without-variants", action="store_true")
    p.add_argument("--label", default=None)
    _add_filter_args(p)
    p.set_defaults(func=cmd_filter)

    p = sub.add_parser("adapters", help="Add constant left/right adapters.")
    p.add_argument("--design", required=True)
    p.add_argument("--left", default="", help="Left adapter sequence.")
    p.add_argument("--right", default="", help="Right adapter sequence.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_adapters)

    p = sub.add_parser("kmer-filter",
                       help="Pick the least k-mer-redundant candidate per group.")
    p.add_argument("--candidates", required=True,
                   help="FASTA of candidates, ids '<source_id>_<n>'.")
    p.add_argument("--reference-fasta", required=True,
                   help="FASTA of the source sequences.")
    p.add_argument("--kmers", type=int, default=6)
    p.add_argument("--include-overlap", action="store_true",
                   help="Write the shared k-mer count into the FASTA header.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_kmer_filter)

    p = sub.add_parser("pipeline",
                       help="tile -> design -> filter -> adapters, end to end.")
    p.add_argument("--regions", required=True)
    p.add_argument("--variants", default=None,
                   help="Optional VCF; without it this is a region-only design.")
    p.add_argument("--reference", required=True)
    p.add_argument("--tile", action="store_true",
                   help="Tile the regions first (skip if already oligo-sized).")
    p.add_argument("--no-length-check", action="store_true")
    p.add_argument("--left-adapter", default="")
    p.add_argument("--right-adapter", default="")
    p.add_argument("--label", default=None)
    _add_tiling_args(p)
    _add_filter_args(p)
    _add_variant_args(p)
    p.set_defaults(func=cmd_pipeline)

    p = sub.add_parser("write-playbook",
                       help="Write Docs/Skills/MPRA_OLIGO_DESIGN_SKILLS.md.")
    p.set_defaults(func=cmd_write_playbook)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
