"""MPRA count & assignment skill — sequencing reads → per-oligo activity.

Python port of the analysis core of **MPRAsnakeflow** (Max Schubach,
Berlin Institute of Health at Charité; https://github.com/kircherlab/
MPRAsnakeflow, MIT licence). The upstream project is a Snakemake
workflow of ~12.7k lines (Python + R + rules); this module ports the
algorithmic core — barcode→oligo assignment, count normalisation,
outlier removal, replicate aggregation and the allelic variant table —
into one dependency-free IGVFagent subcommand tree.

If you publish results produced with this skill, cite the upstream work:

    Rosen JD, Vasanthakumari AD, Salomon K, de Lange N, Dash PM,
    Keukeleire P, Hassan A, Barrera A, Krupkin B, Oualline G, Kircher M,
    Love MI, Schubach M. Uniform processing and analysis of IGVF
    massively parallel reporter assay data with MPRAsnakeflow.
    Genome Research (2025). doi:10.1101/gr.281462.125

Where this sits
---------------
``igvfagent oligo`` designs the library; this skill processes what comes
back off the sequencer.

    design  →  ORDER  →  transfect  →  SEQUENCE  →  this skill  →  activity

An MPRA reads out activity indirectly. Each oligo carries a random
barcode; you sequence DNA (how much of each oligo went in) and RNA (how
much transcript came out). Activity is the RNA/DNA ratio. Two things
must happen before that ratio means anything:

1. **Assignment** — learn which barcode belongs to which oligo, from a
   separate barcode-to-oligo sequencing run. A barcode is only trusted
   when enough reads agree (``--minimum``) and they agree strongly
   enough (``--fraction``); barcodes whose reads disagree are
   *ambiguous* and are dropped, not guessed.
2. **Normalisation** — an oligo's counts are summed over its barcodes,
   divided by how many barcodes it had, and scaled to counts per
   million, so an oligo with 40 barcodes is comparable to one with 11.

Upstream layout → this module
-----------------------------
  =============================================  ======================
  MPRAsnakeflow                                  here
  scripts/assignment/filterAssignmentTsv.py      ``assign-filter``
  scripts/count/merge_BC_and_assignment.py       ``assign-stats``
  scripts/count/merge_label.py                   ``merge-counts``
  scripts/count/make_master_tables.R             ``master-table``
  scripts/count/combine_replicates.py            ``combine-replicates``
  scripts/count/merge_replicates_barcode_counts.py  ``barcode-matrix``
  scripts/variants/generateVariantTable.py       ``variant-table``
  rules/experiment/*.smk                         ``pipeline``

Deliberate deviations from upstream
-----------------------------------
* **No Snakemake, conda, R, pandas, numpy or pysam.** The R master-table
  step and the pandas group-by/quantile logic are reimplemented on the
  standard library, matching how ``mpra_data_skills`` already does its
  statistics. Same numbers, no toolchain.
* **Alignment is not reimplemented.** Upstream offers five aligner
  backends (bwa, bbmap, pbmm2, exact, hybrid) to turn assignment reads
  into barcode↔oligo pairs. Reimplementing short-read alignment would be
  absurd, so ``assign-filter`` starts from the barcode-sorted
  ``barcode <tab> oligo <tab> quality`` table those aligners produce.
  ``assign-exact`` covers the alignment-free case directly.
* **Outlier quirks reproduced by default, correctable by flag.** Two
  choices in upstream's ``ratio_mad`` filter define the numbers in
  published MPRAsnakeflow output, so reproducibility wins by default and
  each has an opt-in fix:

  - the test is ``ratio_diff <= times * mad``, one-sided. This is
    deliberate, not a bug: Rosen et al. (2025) describe the rule as
    barcodes "with expression ratios exceeding 5 log2-units *above*
    their oligo-specific median". ``--mad-two-sided`` compares
    ``abs(ratio_diff)`` if you want both tails;
  - the quantile bin edges are ``arange(0, n_bins) / n_bins``, which
    stops at the 95th percentile, so every barcode above it lands in no
    bin, gets a NaN MAD, and is dropped — the top 5% of barcodes by RNA
    count are always removed (``--mad-include-top-bin`` closes the range).

Subcommands
-----------
    count-bc            FASTQ → per-barcode counts
    assign-exact        Assignment reads → barcode↔oligo pairs (exact match)
    assign-filter       Barcode↔oligo pairs → trusted assignment
    assign-stats        Assignment coverage statistics
    merge-counts        Barcode counts + assignment → per-oligo activity
    complexity          Lincoln-Petersen library-complexity estimate
    master-table        Combine replicate tables, filter on barcode count
    combine-replicates  Aggregate replicates into one per-oligo table
    barcode-matrix      Per-barcode DNA/RNA matrix across replicates
    variant-table       REF/ALT oligo pairs → allelic log2 skew
    pipeline            merge-counts → master-table → combine → variants
    write-playbook      Write Docs/Skills/MPRA_SNAKEFLOW_SKILLS.md
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "MPRASnakeflow"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

DEFAULT_MIN_SUPPORT = 3
DEFAULT_FRACTION = 0.75
DEFAULT_BC_THRESHOLD = 10
DEFAULT_SCALING = 10 ** 6
NO_BC = "no_BC"


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"mpra_snakeflow_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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


def opener(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_tsv(path: str | Path, *, header: bool = True,
             names: list[str] | None = None) -> Iterator[dict]:
    with opener(path) as fh:
        cols = names
        if header:
            first = fh.readline().rstrip("\r\n")
            if not first:
                return
            cols = first.lstrip("#").split("\t") if names is None else names
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if cols is None:
                cols = [str(i) for i in range(len(f))]
            yield dict(zip(cols, f))


def write_tsv(path: str | Path, rows: Iterable[dict], cols: list[str]) -> int:
    n = 0
    with opener(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c)) for c in cols) + "\n")
            n += 1
    return n


def _fmt(v: Any) -> str:
    """Full float precision on purpose.

    Rounding happens only where upstream rounds it — the master-table
    step (``PRECISION``). Truncating here would quietly degrade every
    downstream ratio.
    """
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        if math.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return repr(v)
    return str(v)


def read_fastq(path: str | Path) -> Iterator[tuple[str, str, str | None]]:
    """Yield ``(id, seq, qual)`` from FASTQ or FASTA, plain or gzipped."""
    with opener(path) as fh:
        first = fh.readline()
        if not first:
            return
        fh.seek(0) if not str(path).endswith(".gz") else None
        if first.startswith(">"):
            with opener(path) as fh2:
                name, chunks = None, []
                for line in fh2:
                    line = line.rstrip("\r\n")
                    if line.startswith(">"):
                        if name is not None:
                            yield name, "".join(chunks), None
                        name, chunks = line[1:].strip(), []
                    elif line:
                        chunks.append(line)
                if name is not None:
                    yield name, "".join(chunks), None
            return
    with opener(path) as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\r\n")
            fh.readline()
            qual = fh.readline().rstrip("\r\n")
            yield head[1:].strip(), seq, qual


# ─── Statistics helpers (stdlib stand-ins for the pandas/numpy calls) ───────

def quantile(sorted_values: list[float], q: float) -> float:
    """numpy's default 'linear' interpolation, on an already-sorted list."""
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * q
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac


def sample_std(values: list[float]) -> float:
    """Sample standard deviation (ddof=1), matching pandas' default."""
    if len(values) < 2:
        return float("nan")
    return statistics.stdev(values)


def median(values: list[float]) -> float:
    return statistics.median(values) if values else float("nan")


def nanmedian(values: list[float]) -> float:
    """Median skipping NaN — pandas' groupby.transform('median') default.

    Infinities are NOT skipped: upstream produces them from log2(0/x) and
    lets them propagate, so they must stay in the ordering.
    """
    clean = [v for v in values if not math.isnan(v)]
    return statistics.median(clean) if clean else float("nan")


def log2(x: float) -> float:
    if x is None or x <= 0 or math.isnan(x):
        return float("nan")
    return math.log2(x)


# ─── Barcode counting ───────────────────────────────────────────────────────

def count_barcodes(fastq: str | Path, *, start: int = 0,
                   length: int | None = None, min_quality: int = 0
                   ) -> tuple[dict[str, int], dict[str, int]]:
    """Count reads per barcode.

    ``min_quality`` drops a read when any base of its barcode falls below
    that Phred score — a miscalled base invents a barcode that will never
    match the assignment, so it is cheaper to discard the read here.
    """
    counts: dict[str, int] = defaultdict(int)
    stats = {"reads": 0, "low_quality": 0, "ambiguous_base": 0, "too_short": 0}
    for _, seq, qual in read_fastq(fastq):
        stats["reads"] += 1
        end = start + length if length else len(seq)
        bc = seq[start:end]
        if len(bc) < (length or 1):
            stats["too_short"] += 1
            continue
        if "N" in bc.upper():
            stats["ambiguous_base"] += 1
            continue
        if min_quality and qual:
            q = qual[start:end]
            if any(ord(c) - 33 < min_quality for c in q):
                stats["low_quality"] += 1
                continue
        counts[bc.upper()] += 1
    return dict(counts), stats


# ─── Assignment ─────────────────────────────────────────────────────────────

def assign_exact(design_fasta: str | Path, reads: str | Path,
                 barcodes: str | Path | None = None, *,
                 barcode_start: int = 0, barcode_length: int | None = None
                 ) -> tuple[list[tuple[str, str, str]], dict[str, int]]:
    """Alignment-free assignment: exact insert match against the design.

    Upstream reaches the same table via bwa/bbmap/pbmm2; this covers the
    ``mapping_exact`` case without an aligner. Reads whose insert is not
    in the design are reported as ``other``, exactly as the aligner path
    does, so downstream filtering behaves identically.
    """
    lookup: dict[str, str] = {}
    for name, seq, _ in read_fastq(design_fasta):
        lookup[seq.upper()] = name.split()[0]
    logging.info("design: %d unique oligo sequences", len(lookup))

    bc_by_read: dict[str, str] = {}
    if barcodes:
        for name, seq, _ in read_fastq(barcodes):
            end = barcode_start + barcode_length if barcode_length else len(seq)
            bc_by_read[name.split()[0]] = seq[barcode_start:end].upper()

    pairs: list[tuple[str, str, str]] = []
    stats = {"reads": 0, "matched": 0, "other": 0, "no_barcode": 0}
    for name, seq, _ in read_fastq(reads):
        stats["reads"] += 1
        rid = name.split()[0]
        if bc_by_read:
            bc = bc_by_read.get(rid)
            if bc is None:
                stats["no_barcode"] += 1
                continue
        else:
            # No index read: the barcode is the read-name suffix after ':'.
            bc = rid.rsplit(":", 1)[-1].upper()
        oligo = lookup.get(seq.upper())
        if oligo is None:
            stats["other"] += 1
            oligo = "other"
        else:
            stats["matched"] += 1
        pairs.append((bc, oligo, "exact"))
    pairs.sort()
    return pairs, stats


def _most_frequent(counter: dict[str, int]) -> str:
    """Port of ``filterAssignmentTsv.mfreq`` — modal key, ties by order."""
    if not counter:
        return "NA"
    best = max(counter.values())
    for key, value in counter.items():
        if value == best:
            return key
    return "NA"


def filter_assignment(pairs: Iterable[tuple[str, str, str]], *,
                      minimum: int = DEFAULT_MIN_SUPPORT,
                      fraction: float = DEFAULT_FRACTION,
                      report_other: bool = False,
                      report_ambiguous: bool = False
                      ) -> tuple[list[dict], dict[str, int]]:
    """Keep a barcode only when its reads agree on one oligo.

    Port of ``filterAssignmentTsv.py``. A barcode is reported when it has
    at least ``minimum`` reads and at least ``fraction`` of them name the
    same oligo. Because ``fraction`` must exceed 0.5 at most one oligo can
    ever qualify, so taking the modal oligo and testing it is equivalent
    to upstream's running-maximum loop, and clearer.

    Deviation: upstream's ambiguous branch reports
    ``mfreq(cquality[insert])`` where ``insert`` is whatever the loop
    happened to end on — an uninitialised-variable leak. This reports the
    quality of the *modal* oligo, which is the evident intent.
    """
    if fraction <= 0.5:
        raise SystemExit("--fraction must be above 0.5 (a barcode cannot have "
                         "two majority oligos).")

    out: list[dict] = []
    stats = {"barcodes": 0, "assigned": 0, "ambiguous": 0,
             "below_minimum": 0, "other": 0}

    current: str | None = None
    assignments: dict[str, int] = defaultdict(int)
    quality: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    def flush(barcode: str | None) -> None:
        if barcode is None:
            return
        total = sum(assignments.values())
        if total == 0:
            return
        stats["barcodes"] += 1
        oligo = max(assignments, key=lambda k: assignments[k])
        count = assignments[oligo]
        if total >= minimum and count >= minimum and count / total >= fraction:
            if oligo == "other" and not report_other:
                stats["other"] += 1
                return
            stats["assigned"] += 1
            out.append({"barcode": barcode, "oligo": oligo,
                        "quality": _most_frequent(quality[oligo]),
                        "support": f"{count}/{total}"})
        else:
            if total < minimum or count < minimum:
                stats["below_minimum"] += 1
            else:
                stats["ambiguous"] += 1
            if report_ambiguous:
                out.append({"barcode": barcode, "oligo": "ambiguous",
                            "quality": _most_frequent(quality[oligo]),
                            "support": f"{count}/{total}"})

    for barcode, oligo, qual in pairs:
        if barcode != current:
            flush(current)
            if current is not None and barcode < current:
                raise SystemExit(
                    "Input is not sorted by barcode. Sort it first "
                    "(`sort -k1,1`) — the filter streams one barcode at a time.")
            current = barcode
            assignments = defaultdict(int)
            quality = defaultdict(lambda: defaultdict(int))
        assignments[oligo] += 1
        quality[oligo][qual] += 1
    flush(current)
    return out, stats


def assignment_statistics(assignment: list[dict]) -> dict[str, Any]:
    """Barcodes-per-oligo distribution — the headline assignment QC."""
    per_oligo: dict[str, int] = defaultdict(int)
    for row in assignment:
        per_oligo[row["oligo"]] += 1
    counts = sorted(per_oligo.values())
    if not counts:
        return {"oligos": 0, "barcodes": 0}
    return {
        "oligos": len(per_oligo),
        "barcodes": len(assignment),
        "median_barcodes_per_oligo": median([float(c) for c in counts]),
        "mean_barcodes_per_oligo": sum(counts) / len(counts),
        "min_barcodes_per_oligo": counts[0],
        "max_barcodes_per_oligo": counts[-1],
        "q25": quantile([float(c) for c in counts], 0.25),
        "q75": quantile([float(c) for c in counts], 0.75),
        "oligos_below_10_barcodes": sum(1 for c in counts if c < 10),
    }


# ─── Count merging + normalisation (port of scripts/count/merge_label.py) ───

class BarcodeRow:
    __slots__ = ("barcode", "oligo", "dna", "rna", "ratio", "bin")

    def __init__(self, barcode: str, oligo: str, dna: int, rna: int):
        self.barcode = barcode
        self.oligo = oligo
        self.dna = dna
        self.rna = rna
        self.ratio = 0.0
        self.bin = -1


def outlier_removal_by_rna_zscore(rows: list[BarcodeRow], times: float = 3.0
                                  ) -> tuple[list[BarcodeRow], list[BarcodeRow]]:
    """Drop barcodes whose RNA count is an outlier within their oligo.

    Port of ``outlier_removal_by_rna_zscore``. An oligo with a single
    barcode has undefined (NaN) standard deviation, so its z-score
    comparison is False and the barcode is dropped — upstream behaviour,
    reproduced.
    """
    by_oligo: dict[str, list[BarcodeRow]] = defaultdict(list)
    for r in rows:
        by_oligo[r.oligo].append(r)
    kept, removed = [], []
    for group in by_oligo.values():
        vals = [float(r.rna) for r in group]
        mean = sum(vals) / len(vals)
        sd = sample_std(vals)
        for r in group:
            if math.isnan(sd) or sd == 0:
                removed.append(r)
                continue
            if abs((r.rna - mean) / sd) <= times:
                kept.append(r)
            else:
                removed.append(r)
    return kept, removed


def outlier_removal_by_mad(rows: list[BarcodeRow], n_bins: int = 20,
                           times: float = 5.0, *, two_sided: bool = False,
                           include_top_bin: bool = False
                           ) -> tuple[list[BarcodeRow], list[BarcodeRow]]:
    """Drop barcodes whose DNA/RNA ratio deviates from their oligo's median.

    Port of ``outlier_removal_by_mad``. Barcodes are binned by RNA depth
    (deep barcodes are intrinsically less noisy) and, within each bin, the
    median absolute deviation sets the tolerance.

    Faithfulness notes: upstream computes ``log2(dna/rna)`` with numpy, so
    a zero-DNA barcode yields ``-inf`` rather than being excluded, and
    those infinities take part in the per-oligo medians. ``-inf - -inf``
    is NaN, and pandas' median skips NaN while a NaN comparison is False,
    so such barcodes are dropped. All of that is reproduced here — it
    changes which barcodes survive. See the module docstring for the two
    upstream quirks that are reproduced by default.
    """
    if not rows:
        return rows, []

    def safe_log2(num: float, den: float) -> float:
        if den == 0:
            return float("nan")
        if num == 0:
            return float("-inf")
        return math.log2(num / den)

    ratio = {id(r): safe_log2(r.dna, r.rna) for r in rows}

    by_oligo: dict[str, list[BarcodeRow]] = defaultdict(list)
    for r in rows:
        by_oligo[r.oligo].append(r)
    ratio_diff: dict[int, float] = {}
    for group in by_oligo.values():
        med = nanmedian([ratio[id(r)] for r in group])
        for r in group:
            a, b = ratio[id(r)], med
            if math.isnan(a) or math.isnan(b) or (math.isinf(a) and a == b):
                ratio_diff[id(r)] = float("nan")   # -inf - -inf -> NaN
            else:
                ratio_diff[id(r)] = a - b

    # log10(rna); rna == 0 gives -inf, which falls outside every bin.
    log_rna = {id(r): (math.log10(r.rna) if r.rna > 0 else float("-inf"))
               for r in rows}
    finite = sorted(v for v in log_rna.values() if math.isfinite(v))
    if not finite:
        return rows, []
    fractions = [i / n_bins for i in range(n_bins)]
    if include_top_bin:
        fractions.append(1.0)
    edges = sorted({quantile(finite, q) for q in fractions})
    if len(edges) < 2:
        return rows, []

    def find_bin(x: float) -> int:
        # pandas cut(include_lowest=True): [e0, e1], then (e1, e2], ...
        if not math.isfinite(x) or x < edges[0] or x > edges[-1]:
            return -1
        for i in range(1, len(edges)):
            if x <= edges[i]:
                return i - 1
        return -1

    by_bin: dict[int, list[BarcodeRow]] = defaultdict(list)
    for r in rows:
        r.bin = find_bin(log_rna[id(r)])
        if r.bin >= 0:
            by_bin[r.bin].append(r)

    mad_by_bin: dict[int, float] = {}
    for b, group in by_bin.items():
        diffs = [ratio_diff[id(r)] for r in group]
        med = nanmedian(diffs)
        dist = [abs(d - med) if not (math.isnan(d) or math.isnan(med))
                else float("nan") for d in diffs]
        mad_by_bin[b] = nanmedian(dist)

    kept, removed = [], []
    for r in rows:
        mad = mad_by_bin.get(r.bin, float("nan"))
        d = ratio_diff[id(r)]
        if math.isnan(mad) or math.isnan(d):
            removed.append(r)          # NaN comparison is False upstream
            continue
        value = abs(d) if two_sided else d
        (kept if value <= times * mad else removed).append(r)
    return kept, removed


def merge_counts(count_rows: Iterable[dict], assignment_rows: Iterable[dict], *,
                 min_dna: int = 0, min_rna: int = 1,
                 normalize_with_unassigned: bool = False,
                 scaling: float = DEFAULT_SCALING,
                 outlier_method: str | None = None,
                 mad_bins: int = 20, mad_times: float = 5.0,
                 zscore_times: float = 3.0,
                 mad_two_sided: bool = False,
                 mad_include_top_bin: bool = False
                 ) -> dict[str, Any]:
    """Barcode counts + assignment → per-oligo normalised activity.

    Port of ``merge_label.py``. The normalisation is the important part:
    an oligo's counts are summed across its barcodes, divided by how many
    barcodes it had, then scaled per million — so oligos with different
    barcode counts stay comparable.
    """
    pseudo_dna = 1 if min_dna == 0 else 0
    pseudo_rna = 1 if min_rna == 0 else 0

    # A barcode seen against two oligos is evidence of nothing; upstream
    # drops every copy (drop_duplicates keep=False), not just the extras.
    seen: dict[str, str] = {}
    duplicated: set[str] = set()
    for row in assignment_rows:
        bc, oligo = row["barcode"], row["oligo"]
        if bc in seen and seen[bc] != oligo:
            duplicated.add(bc)
        elif bc in seen:
            duplicated.add(bc)
        seen[bc] = oligo
    assoc = {bc: o for bc, o in seen.items() if bc not in duplicated}

    stats: dict[str, Any] = {
        "oligos design": len({o for o in assoc.values()}),
        "barcodes design": len(assoc),
    }

    def _count(v: Any) -> int:
        # Real IGVF reporter-experiment files leave a cell empty when a
        # barcode was not observed in that replicate; treat it as zero
        # rather than crashing (min_rna then drops it, as upstream would).
        if v is None:
            return 0
        v = str(v).strip()
        if v in ("", "NA", "NaN", "nan", "."):
            return 0
        return int(float(v))

    rows: list[BarcodeRow] = []
    for row in count_rows:
        dna = _count(row.get("dna_count"))
        rna = _count(row.get("rna_count"))
        if dna < min_dna or rna < min_rna:
            continue
        bc = row["barcode"]
        rows.append(BarcodeRow(bc, assoc.get(bc, NO_BC), dna, rna))

    stats["barcodes dna/rna"] = len(rows)
    unknown = [r for r in rows if r.oligo == NO_BC]
    stats["unknown barcodes dna/rna"] = len(unknown)
    stats["matched barcodes"] = len(rows) - len(unknown)
    stats["% matched barcodes"] = (
        100.0 * stats["matched barcodes"] / len(rows) if rows else 0.0)

    total_dna = sum(r.dna for r in rows)
    total_rna = sum(r.rna for r in rows)
    stats["total dna counts"] = total_dna
    stats["total rna counts"] = total_rna
    stats["avg dna counts per bc"] = total_dna / len(rows) if rows else 0.0
    stats["avg rna counts per bc"] = total_rna / len(rows) if rows else 0.0

    if not normalize_with_unassigned:
        total_dna -= sum(r.dna for r in unknown)
        total_rna -= sum(r.rna for r in unknown)
        rows = [r for r in rows if r.oligo != NO_BC]

    removed: list[BarcodeRow] = []
    if outlier_method == "ratio_mad":
        rows, removed = outlier_removal_by_mad(
            rows, mad_bins, mad_times,
            two_sided=mad_two_sided, include_top_bin=mad_include_top_bin)
    elif outlier_method == "rna_counts_zscore":
        rows, removed = outlier_removal_by_rna_zscore(rows, zscore_times)
    elif outlier_method:
        raise SystemExit(f"Unknown outlier method {outlier_method!r}.")
    stats["barcode outlier removed"] = len(removed)

    barcode_rows = [{"barcode": r.barcode, "oligo_name": r.oligo,
                     "dna_count": r.dna, "rna_count": r.rna} for r in rows]

    for r in rows:
        r.dna += pseudo_dna
        r.rna += pseudo_rna
    if pseudo_dna:
        total_dna = sum(r.dna for r in rows)
    if pseudo_rna:
        total_rna = sum(r.rna for r in rows)

    grouped: dict[str, list[BarcodeRow]] = defaultdict(list)
    for r in rows:
        grouped[r.oligo].append(r)

    output: list[dict] = []
    for oligo, group in grouped.items():
        n_bc = len(group)
        dna_sum = sum(r.dna for r in group)
        rna_sum = sum(r.rna for r in group)
        dna_norm = (dna_sum / n_bc) / total_dna * scaling if total_dna else float("nan")
        rna_norm = (rna_sum / n_bc) / total_rna * scaling if total_rna else float("nan")
        ratio = rna_norm / dna_norm if dna_norm else float("nan")
        output.append({
            "oligo_name": oligo,
            "dna_counts": dna_sum,
            "rna_counts": rna_sum,
            "dna_normalized": dna_norm,
            "rna_normalized": rna_norm,
            "ratio": ratio,
            "log2FoldChange": log2(ratio),
            "n_bc": n_bc,
        })
    output.sort(key=lambda d: d["oligo_name"])

    stats["oligos dna/rna"] = len(grouped)
    assigned_bcs = sum(len(g) for o, g in grouped.items() if o != NO_BC)
    stats["avg dna/rna barcodes per oligo"] = (
        assigned_bcs / len(grouped) if grouped else 0.0)

    return {"table": output, "statistic": stats,
            "barcodes": barcode_rows,
            "removed_barcodes": [r.barcode for r in removed]}


# ─── Library complexity (Lincoln-Petersen) ─────────────────────────────────
#
# Rosen et al. (2025), "Complexity and sequencing depth estimation": pairs
# of replicates estimate the true library size by mark-recapture. Each
# replicate is a "capture" of the barcode pool; barcodes seen in both are
# the "recaptures". The gap between what you observed and the estimate is
# how much of the library your sequencing missed.

def lincoln_index(n1: int, n2: int, shared: int) -> float:
    """Lincoln-Petersen estimate of total population size."""
    if not shared:
        return float("nan")
    return n1 * n2 / shared


def read_wide_barcode_file(path: str | Path
                           ) -> tuple[dict[str, int], dict[tuple[str, str], int],
                                      list[str], int]:
    """Read an IGVF 'reporter experiment barcode' file.

    The IGVF standard packs all replicates into one wide table
    (``barcode, oligo_name, dna_count_<r>, rna_count_<r>, ...``) and leaves
    both count cells empty where a barcode was not observed.

    Because that layout gives exactly one row per barcode, the counts and
    pairwise overlaps can be accumulated in a single streaming pass — no
    barcode sets are held in memory. That matters: the 240K libraries have
    ~13 M barcodes, where three string sets would cost several GB.

    Returns ``(observed_per_replicate, shared_per_pair, replicates, rows)``.
    """
    with opener(path) as fh:
        cols = fh.readline().rstrip("\r\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        reps = sorted({c[len("dna_count_"):] for c in cols
                       if c.startswith("dna_count_")})
        if not reps:
            raise SystemExit(
                f"No dna_count_<replicate> columns in {path}. Expected an "
                "IGVF 'reporter experiment barcode' file.")
        di = [idx[f"dna_count_{r}"] for r in reps]
        ri = [idx[f"rna_count_{r}"] for r in reps]
        counts = [0] * len(reps)
        pair_idx = [(a, b) for a in range(len(reps)) for b in range(a + 1, len(reps))]
        shared = [0] * len(pair_idx)
        n = 0
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) < len(cols):
                continue
            n += 1
            # "requiring barcodes to have at least one count of DNA and RNA"
            # (Rosen et al.) — both modalities must be present.
            seen = 0
            for k in range(len(reps)):
                d = f[di[k]]
                v = f[ri[k]]
                if d and v and d not in ("0", "NA") and v not in ("0", "NA"):
                    counts[k] += 1
                    seen |= 1 << k
            if seen:
                for j, (a, b) in enumerate(pair_idx):
                    if seen >> a & 1 and seen >> b & 1:
                        shared[j] += 1
    observed = {r: counts[k] for k, r in enumerate(reps)}
    pairs = {(reps[a], reps[b]): shared[j]
             for j, (a, b) in enumerate(pair_idx)}
    return observed, pairs, reps, n


def complexity_metrics(observed: dict[str, int],
                       shared: dict[tuple[str, str], int],
                       reps: list[str]) -> dict[str, Any]:
    """Per-replicate barcode counts, pairwise Lincoln indices, and the
    fraction of the library each replicate is missing."""
    pairs = []
    for i, a in enumerate(reps):
        for b in reps[i + 1:]:
            m = shared.get((a, b), 0)
            pairs.append({"replicates": f"{a}+{b}",
                          "n_a": observed[a], "n_b": observed[b],
                          "shared": m,
                          "lincoln_index": lincoln_index(
                              observed[a], observed[b], m)})
    med_obs = median([float(v) for v in observed.values()])
    lincolns = [p["lincoln_index"] for p in pairs
                if not math.isnan(p["lincoln_index"])]
    med_lin = median(lincolns) if lincolns else float("nan")
    missing = ((1 - med_obs / med_lin) * 100
               if med_lin and not math.isnan(med_lin) else float("nan"))
    return {
        "replicates": reps,
        "observed_barcodes": observed,
        "median_observed_barcodes": med_obs,
        "pairwise": pairs,
        "median_lincoln_index": med_lin,
        "percent_missing": missing,
    }


# ─── Replicate aggregation (port of make_master_tables.R + friends) ─────────

MASTER_COLS = ["replicate", "oligo_name", "dna_counts", "rna_counts",
               "dna_normalized", "rna_normalized", "log2FoldChange", "n_bc"]
PRECISION = 4


def master_table(replicates: list[tuple[str, str]], *,
                 threshold: int = DEFAULT_BC_THRESHOLD
                 ) -> tuple[list[dict], list[dict], list[dict]]:
    """Stack per-replicate tables and filter on barcode support.

    Port of ``make_master_tables.R``. Returns
    ``(all_rows, filtered_rows, per_replicate_statistics)``. Oligos with
    fewer than ``threshold`` barcodes in a replicate are dropped from the
    filtered table: a ratio built on three barcodes is mostly noise.
    """
    rows: list[dict] = []
    for name, path in replicates:
        for r in read_tsv(path):
            if r.get("oligo_name") == NO_BC:
                continue
            row = {"replicate": name, "oligo_name": r["oligo_name"],
                   "dna_counts": int(float(r["dna_counts"])),
                   "rna_counts": int(float(r["rna_counts"])),
                   "n_bc": int(float(r["n_bc"]))}
            for k in ("dna_normalized", "rna_normalized", "ratio",
                      "log2FoldChange"):
                v = r.get(k)
                row[k] = round(float(v), PRECISION) if v not in (None, "", "NA") \
                    else float("nan")
            rows.append(row)

    filtered = [r for r in rows if r["n_bc"] >= threshold]

    def averages(subset: list[dict], name: str) -> list[dict]:
        by_rep: dict[str, list[dict]] = defaultdict(list)
        for r in subset:
            by_rep[r["replicate"]].append(r)
        out = []
        for rep, group in by_rep.items():
            ratios = [g["ratio"] for g in group if not math.isnan(g["ratio"])]
            lfcs = [g["log2FoldChange"] for g in group
                    if not math.isnan(g["log2FoldChange"])]
            out.append({
                "replicate": rep,
                "mean_ratio": sum(ratios) / len(ratios) if ratios else float("nan"),
                "mean_log2FoldChange": sum(lfcs) / len(lfcs) if lfcs else float("nan"),
                "mean_n_bc": sum(g["n_bc"] for g in group) / len(group),
                "BC_filter": name,
            })
        return out

    stats = averages(rows, "None") + averages(filtered, f"n_bc >= {threshold}")
    return rows, filtered, stats


def combine_replicates(master: list[dict],
                       labels: dict[str, str] | None = None) -> list[dict]:
    """Aggregate replicates into one row per oligo.

    Port of ``combine_replicates.py``. Note the two different summaries:
    ``log2FoldChange`` is recomputed from *pooled* counts, while
    ``mean_log2FoldChange`` averages the per-replicate values. They answer
    different questions and upstream reports both.
    """
    total_dna = sum(r["dna_counts"] for r in master)
    total_rna = sum(r["rna_counts"] for r in master)
    scaling = DEFAULT_SCALING

    by_oligo: dict[str, list[dict]] = defaultdict(list)
    for r in master:
        by_oligo[r["oligo_name"]].append(r)

    def mean_of(group: list[dict], key: str) -> float:
        vals = [g[key] for g in group
                if isinstance(g.get(key), (int, float)) and not math.isnan(g[key])]
        return sum(vals) / len(vals) if vals else float("nan")

    out: list[dict] = []
    for oligo, group in by_oligo.items():
        dna_counts = sum(g["dna_counts"] for g in group)
        rna_counts = sum(g["rna_counts"] for g in group)
        dna_norm = dna_counts / total_dna * scaling if total_dna else float("nan")
        rna_norm = rna_counts / total_rna * scaling if total_rna else float("nan")
        row = {
            "oligo_name": oligo,
            "replicates": len(group),
            "dna_counts": dna_counts,
            "rna_counts": rna_counts,
            "dna_normalized": dna_norm,
            "rna_normalized": rna_norm,
            "log2FoldChange": log2(rna_norm / dna_norm) if dna_norm else float("nan"),
            "mean_dna_counts": mean_of(group, "dna_counts"),
            "mean_rna_counts": mean_of(group, "rna_counts"),
            "mean_dna_normalized": mean_of(group, "dna_normalized"),
            "mean_rna_normalized": mean_of(group, "rna_normalized"),
            "mean_log2FoldChange": mean_of(group, "log2FoldChange"),
            "mean_n_bc": mean_of(group, "n_bc"),
            "n_bc": sum(g["n_bc"] for g in group),
        }
        if labels:
            row["Label"] = labels.get(oligo, "")
        out.append(row)
    out.sort(key=lambda d: d["oligo_name"])
    return out


def barcode_matrix(replicates: list[tuple[str, str]], *,
                   threshold: int = DEFAULT_BC_THRESHOLD
                   ) -> tuple[list[dict], list[dict], list[str]]:
    """Per-barcode DNA/RNA counts pivoted across replicates.

    Port of ``merge_replicates_barcode_counts.py``. This is the table the
    downstream per-barcode models (MPRAlib, MPRAnalyze) consume.
    """
    names = [n for n, _ in replicates]
    cells: dict[tuple[str, str], dict[str, int]] = {}
    per_oligo_rep: dict[tuple[str, str], int] = defaultdict(int)
    for name, path in replicates:
        for r in read_tsv(path):
            oligo = r.get("oligo_name")
            if oligo == NO_BC:
                continue
            key = (r["barcode"], oligo)
            cells.setdefault(key, {})
            # Accept the standard suffixed columns or the older plain
            # dna_count/rna_count, so existing files still load.
            dna = r.get(f"dna_count_{name}", r.get("dna_count"))
            rna = r.get(f"rna_count_{name}", r.get("rna_count"))
            cells[key][f"dna_count_{name}"] = int(float(dna))
            cells[key][f"rna_count_{name}"] = int(float(rna))
            per_oligo_rep[(oligo, name)] += 1

    cols = ["barcode", "oligo_name"] + [
        f"{kind}_count_{n}" for n in names for kind in ("dna", "rna")]

    def build(keep: set[tuple[str, str]] | None) -> list[dict]:
        rows = []
        for (bc, oligo), vals in cells.items():
            if keep is not None and (oligo, bc) not in keep:
                continue
            rows.append({"barcode": bc, "oligo_name": oligo, **vals})
        rows.sort(key=lambda d: (d["oligo_name"], d["barcode"]))
        return rows

    passing = {(o, n) for (o, n), c in per_oligo_rep.items() if c >= threshold}
    keep_keys = {(oligo, bc) for (bc, oligo) in cells
                 if any((oligo, n) in passing for n in names)}
    return build(None), build(keep_keys), cols


def variant_table(counts: list[dict], declaration: list[dict]) -> list[dict]:
    """REF/ALT oligo pairs → allelic log2 skew.

    Port of ``generateVariantTable.py``. The declaration file names, for
    each variant, which oligo is REF and which is ALT; the skew is
    ``log2(ratio_ALT / ratio_REF)`` — how much the alternate allele
    changes reporter activity.
    """
    by_name = {r["oligo_name"]: r for r in counts}
    numeric = ("dna_counts", "rna_counts", "dna_normalized",
               "rna_normalized", "ratio", "log2FoldChange", "n_bc")
    out: list[dict] = []
    for d in declaration:
        row = {"ID": d["ID"], "REF": d["REF"], "ALT": d["ALT"]}
        for side in ("REF", "ALT"):
            rec = by_name.get(d[side])
            for k in numeric:
                v = rec.get(k) if rec else 0
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    v = 0
                row[f"{k}_{side}"] = v
        r_ref, r_alt = row["ratio_REF"], row["ratio_ALT"]
        row["log2FoldChange_expression"] = (
            log2(r_alt / r_ref) if r_ref and r_alt else float("nan"))
        out.append(row)
    return out


def master_variant_table(per_replicate: list[list[dict]], *,
                         min_dna: int = 0, min_rna: int = 1) -> list[dict]:
    """Pool per-replicate variant tables into one table per variant.

    Port of ``generateMasterVariantTable.py``. Counts are summed across
    replicates and renormalised against the pooled totals, so this is a
    pooled re-estimate rather than an average of per-replicate skews.
    """
    pseudo_dna = 1 if min_dna == 0 else 0
    pseudo_rna = 1 if min_rna == 0 else 0

    grouped: dict[tuple[str, str, str], dict[str, int]] = {}
    for table in per_replicate:
        for r in table:
            if r["dna_counts_REF"] < min_dna or r["dna_counts_ALT"] < min_dna:
                continue
            if r["rna_counts_REF"] < min_rna or r["rna_counts_ALT"] < min_rna:
                continue
            key = (r["ID"], r["REF"], r["ALT"])
            acc = grouped.setdefault(key, {k: 0 for k in (
                "dna_counts_REF", "rna_counts_REF", "n_bc_REF",
                "dna_counts_ALT", "rna_counts_ALT", "n_bc_ALT")})
            for k in acc:
                acc[k] += int(r[k])

    total_bc = sum(a["n_bc_REF"] + a["n_bc_ALT"] for a in grouped.values())
    total_dna = sum(a["dna_counts_REF"] + a["dna_counts_ALT"]
                    for a in grouped.values()) + total_bc * pseudo_dna
    total_rna = sum(a["rna_counts_REF"] + a["rna_counts_ALT"]
                    for a in grouped.values()) + total_bc * pseudo_rna
    scaling = DEFAULT_SCALING

    out: list[dict] = []
    for (vid, ref, alt), acc in grouped.items():
        row: dict[str, Any] = {"ID": vid, "REF": ref, "ALT": alt, **acc}
        for side in ("REF", "ALT"):
            n_bc = acc[f"n_bc_{side}"]
            for kind, pseudo, total in (("dna", pseudo_dna, total_dna),
                                        ("rna", pseudo_rna, total_rna)):
                if n_bc and total:
                    row[f"{kind}_normalized_{side}"] = (
                        (acc[f"{kind}_counts_{side}"] + pseudo * n_bc)
                        / n_bc / total * scaling)
                else:
                    row[f"{kind}_normalized_{side}"] = 0.0
            dn = row[f"dna_normalized_{side}"]
            row[f"ratio_{side}"] = (
                row[f"rna_normalized_{side}"] / dn if dn else 0.0)
            lfc = log2(row[f"ratio_{side}"])
            row[f"log2FoldChange_{side}"] = 0.0 if math.isnan(lfc) else lfc
        r_ref, r_alt = row["ratio_REF"], row["ratio_ALT"]
        lfc = log2(r_alt / r_ref) if r_ref and r_alt else float("nan")
        row["log2FoldChange_expression"] = 0.0 if math.isnan(lfc) else lfc
        out.append(row)
    out.sort(key=lambda d: d["ID"])
    return out


MASTER_VARIANT_COLS = ["ID", "REF", "ALT",
                       "dna_counts_REF", "rna_counts_REF",
                       "dna_normalized_REF", "rna_normalized_REF",
                       "ratio_REF", "log2FoldChange_REF", "n_bc_REF",
                       "dna_counts_ALT", "rna_counts_ALT",
                       "dna_normalized_ALT", "rna_normalized_ALT",
                       "ratio_ALT", "log2FoldChange_ALT", "n_bc_ALT",
                       "log2FoldChange_expression"]


VARIANT_COLS = (["ID", "REF", "ALT"]
                + [f"{k}_{s}" for s in ("REF", "ALT")
                   for k in ("dna_counts", "rna_counts", "dna_normalized",
                             "rna_normalized", "ratio", "log2FoldChange", "n_bc")]
                + ["log2FoldChange_expression"])


# ─── Commands ───────────────────────────────────────────────────────────────

def _parse_replicates(specs: list[str]) -> list[tuple[str, str]]:
    out = []
    for spec in specs or []:
        if "=" not in spec:
            raise SystemExit(
                f"--replicate expects NAME=PATH, got {spec!r} "
                "(e.g. --replicate rep1=counts_rep1.tsv.gz)")
        name, path = spec.split("=", 1)
        if not Path(path).exists():
            raise SystemExit(f"Replicate {name}: file not found: {path}")
        out.append((name, path))
    return out


def _coerce_count_row(r: dict) -> dict:
    row = {"oligo_name": r["oligo_name"]}
    for k in ("dna_counts", "rna_counts", "n_bc"):
        row[k] = int(float(r[k]))
    for k in ("dna_normalized", "rna_normalized", "ratio", "log2FoldChange"):
        v = r.get(k)
        row[k] = float(v) if v not in (None, "", "NA") else float("nan")
    return row


def cmd_count_bc(args) -> int:
    setup_logging()
    counts, stats = count_barcodes(
        args.fastq, start=args.barcode_start, length=args.barcode_length,
        min_quality=args.min_quality)
    out_dir = _run_dir(args.label or "count_bc")
    out = out_dir / "barcode_counts.tsv.gz"
    rows = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    with opener(out, "wt") as fh:
        for bc, n in rows:
            fh.write(f"{bc}\t{n}\n")
    print(f"Reads:            {stats['reads']:,}")
    print(f"  dropped, N base {stats['ambiguous_base']:,}")
    print(f"  dropped, low Q  {stats['low_quality']:,}")
    print(f"  dropped, short  {stats['too_short']:,}")
    print(f"Barcodes:         {len(counts):,}")
    if counts:
        vals = sorted(counts.values())
        print(f"Reads per barcode: median {median([float(v) for v in vals]):.0f}, "
              f"max {vals[-1]:,}")
    print(f"Saved:            {out}")
    return 0


def cmd_assign_exact(args) -> int:
    setup_logging()
    pairs, stats = assign_exact(
        args.design, args.reads, args.barcodes,
        barcode_start=args.barcode_start, barcode_length=args.barcode_length)
    out_dir = _run_dir(args.label or "assign_exact")
    out = out_dir / "barcode_oligo_pairs.tsv.gz"
    with opener(out, "wt") as fh:
        for bc, oligo, qual in pairs:
            fh.write(f"{bc}\t{oligo}\t{qual}\n")
    print(f"Reads:      {stats['reads']:,}")
    print(f"  matched:  {stats['matched']:,}")
    print(f"  other:    {stats['other']:,}  (insert not in the design)")
    if stats["no_barcode"]:
        print(f"  no BC:    {stats['no_barcode']:,}")
    print(f"Pairs:      {len(pairs):,}  (barcode-sorted, ready for assign-filter)")
    print(f"Saved:      {out}")
    return 0


def cmd_assign_filter(args) -> int:
    setup_logging()
    pairs = ((r["barcode"], r["oligo"], r.get("quality", ""))
             for r in read_tsv(args.pairs, header=False,
                               names=["barcode", "oligo", "quality"]))
    rows, stats = filter_assignment(
        pairs, minimum=args.minimum, fraction=args.fraction,
        report_other=args.report_other, report_ambiguous=args.report_ambiguous)
    out_dir = _run_dir(args.label or "assign_filter")
    out = out_dir / "assignment.tsv.gz"
    with opener(out, "wt") as fh:
        for r in rows:
            fh.write(f"{r['barcode']}\t{r['oligo']}\t{r['quality']}\t{r['support']}\n")
    qc = assignment_statistics([r for r in rows if r["oligo"] != "ambiguous"])
    (out_dir / "statistic.json").write_text(json.dumps(
        {"filter": stats, "coverage": qc}, indent=2))
    print(f"Barcodes seen:    {stats['barcodes']:,}")
    print(f"  assigned:       {stats['assigned']:,}")
    print(f"  ambiguous:      {stats['ambiguous']:,}  (reads disagree)")
    print(f"  below minimum:  {stats['below_minimum']:,}  (< {args.minimum} reads)")
    print(f"  'other':        {stats['other']:,}")
    if qc.get("oligos"):
        print(f"Oligos covered:   {qc['oligos']:,}  "
              f"(median {qc['median_barcodes_per_oligo']:.0f} barcodes each, "
              f"{qc['oligos_below_10_barcodes']:,} below 10)")
    print(f"Saved:            {out}")
    return 0


def cmd_assign_stats(args) -> int:
    setup_logging()
    rows = [{"barcode": r["barcode"], "oligo": r["oligo"]}
            for r in read_tsv(args.assignment, header=False,
                              names=["barcode", "oligo", "quality", "support"])]
    qc = assignment_statistics(rows)
    out_dir = _run_dir(args.label or "assign_stats")
    (out_dir / "assignment_statistic.json").write_text(json.dumps(qc, indent=2))
    for k, v in qc.items():
        print(f"  {k:32s} {v:,.2f}" if isinstance(v, float) else f"  {k:32s} {v:,}")
    print(f"Saved: {out_dir / 'assignment_statistic.json'}")
    return 0


def cmd_merge_counts(args) -> int:
    setup_logging()
    counts = read_tsv(args.counts, header=False,
                      names=["barcode", "dna_count", "rna_count"])
    assignment = read_tsv(args.assignment, header=False,
                          names=["barcode", "oligo", "quality", "support"])
    res = merge_counts(
        counts, assignment, min_dna=args.min_dna_counts,
        min_rna=args.min_rna_counts,
        normalize_with_unassigned=args.include_unassigned_for_normalization,
        scaling=args.scaling, outlier_method=args.outlier_detection,
        mad_bins=args.outlier_mad_bins, mad_times=args.outlier_mad_times,
        zscore_times=args.outlier_zscore_times,
        mad_two_sided=args.mad_two_sided,
        mad_include_top_bin=args.mad_include_top_bin)
    out_dir = _run_dir(args.label or "merge_counts")
    cols = ["oligo_name", "dna_counts", "rna_counts", "dna_normalized",
            "rna_normalized", "ratio", "log2FoldChange", "n_bc"]
    write_tsv(out_dir / "counts.tsv.gz", res["table"], cols)
    write_tsv(out_dir / "barcode_counts.tsv.gz", res["barcodes"],
              ["barcode", "oligo_name", "dna_count", "rna_count"])
    (out_dir / "statistic.json").write_text(json.dumps(res["statistic"], indent=2))
    if res["removed_barcodes"]:
        (out_dir / "removed_barcodes.txt").write_text(
            "\n".join(res["removed_barcodes"]) + "\n")
    s = res["statistic"]
    print(f"Design:            {s['oligos design']:,} oligos / {s['barcodes design']:,} barcodes")
    print(f"Barcodes in data:  {s['barcodes dna/rna']:,}")
    print(f"  matched:         {s['matched barcodes']:,}  ({s['% matched barcodes']:.1f}%)")
    print(f"  unassigned:      {s['unknown barcodes dna/rna']:,}")
    if s["barcode outlier removed"]:
        print(f"  outliers removed:{s['barcode outlier removed']:,}  ({args.outlier_detection})")
    print(f"Oligos quantified: {s['oligos dna/rna']:,}  "
          f"({s['avg dna/rna barcodes per oligo']:.1f} barcodes each on average)")
    print(f"Saved:             {out_dir}")
    return 0


def cmd_complexity(args) -> int:
    """Library complexity from an IGVF reporter-experiment-barcode file."""
    setup_logging()
    observed, shared, reps, n = read_wide_barcode_file(args.barcode_file)
    m = complexity_metrics(observed, shared, reps)
    out_dir = _run_dir(args.label or "complexity")
    (out_dir / "complexity.json").write_text(json.dumps(m, indent=2))
    write_tsv(out_dir / "pairwise_lincoln.tsv", m["pairwise"],
              ["replicates", "n_a", "n_b", "shared", "lincoln_index"])
    print(f"Barcodes in file:      {n:,}")
    print(f"Replicates:            {len(reps)}  ({', '.join(reps)})")
    for r in reps:
        print(f"  {r:12s} {m['observed_barcodes'][r]:>12,} barcodes with DNA and RNA")
    print(f"Median observed:       {m['median_observed_barcodes']:,.0f}")
    for p_ in m["pairwise"]:
        print(f"  Lincoln {p_['replicates']:9s} {p_['lincoln_index']:>12,.0f}  "
              f"(shared {p_['shared']:,})")
    print(f"Median Lincoln index:  {m['median_lincoln_index']:,.0f}")
    print(f"Library missing:       {m['percent_missing']:.1f}%  "
          f"(a replicate sees this much less than the estimated true library)")
    print(f"Saved:                 {out_dir}")
    return 0


def cmd_master_table(args) -> int:
    setup_logging()
    reps = _parse_replicates(args.replicate)
    rows, filtered, stats = master_table(reps, threshold=args.threshold)
    out_dir = _run_dir(args.label or "master_table")
    write_tsv(out_dir / "master_table.all.tsv.gz", rows, MASTER_COLS)
    write_tsv(out_dir / "master_table.tsv.gz", filtered, MASTER_COLS)
    write_tsv(out_dir / "statistic.tsv.gz", stats,
              ["replicate", "mean_ratio", "mean_log2FoldChange",
               "mean_n_bc", "BC_filter"])
    print(f"Replicates:   {len(reps)}  ({', '.join(n for n, _ in reps)})")
    print(f"Rows:         {len(rows):,}")
    print(f"  n_bc >= {args.threshold}: {len(filtered):,}  "
          f"({len(rows) - len(filtered):,} dropped for thin barcode support)")
    print(f"Saved:        {out_dir}")
    return 0


def cmd_combine_replicates(args) -> int:
    setup_logging()
    master = []
    for r in read_tsv(args.input):
        row = {"replicate": r["replicate"], "oligo_name": r["oligo_name"],
               "dna_counts": int(float(r["dna_counts"])),
               "rna_counts": int(float(r["rna_counts"])),
               "n_bc": int(float(r["n_bc"]))}
        for k in ("dna_normalized", "rna_normalized", "log2FoldChange"):
            v = r.get(k)
            row[k] = float(v) if v not in (None, "", "NA") else float("nan")
        master.append(row)
    labels = None
    if args.labels:
        labels = {r["oligo_name"]: r["Label"] for r in
                  read_tsv(args.labels, header=False,
                           names=["oligo_name", "Label"])}
    out = combine_replicates(master, labels)
    cols = ["oligo_name", "replicates", "dna_counts", "rna_counts",
            "dna_normalized", "rna_normalized", "log2FoldChange",
            "mean_dna_counts", "mean_rna_counts", "mean_dna_normalized",
            "mean_rna_normalized", "mean_log2FoldChange", "mean_n_bc", "n_bc"]
    if labels:
        cols.append("Label")
    out_dir = _run_dir(args.label or "combine_replicates")
    write_tsv(out_dir / "combined.tsv.gz", out, cols)
    print(f"Oligos:  {len(out):,}")
    print(f"Saved:   {out_dir / 'combined.tsv.gz'}")
    return 0


def cmd_barcode_matrix(args) -> int:
    setup_logging()
    reps = _parse_replicates(args.replicate)
    allrows, filtered, cols = barcode_matrix(reps, threshold=args.threshold)
    out_dir = _run_dir(args.label or "barcode_matrix")
    write_tsv(out_dir / "barcode_matrix.all.tsv.gz", allrows, cols)
    write_tsv(out_dir / "barcode_matrix.tsv.gz", filtered, cols)
    print(f"Replicates: {len(reps)}")
    print(f"Barcodes:   {len(allrows):,}  ({len(filtered):,} after the "
          f">= {args.threshold} barcodes/oligo filter)")
    print(f"Saved:      {out_dir}")
    return 0


def cmd_variant_table(args) -> int:
    setup_logging()
    counts = [_coerce_count_row(r) for r in read_tsv(args.counts)]
    declaration = list(read_tsv(args.declaration))
    missing = [c for c in ("ID", "REF", "ALT")
               if declaration and c not in declaration[0]]
    if missing:
        raise SystemExit(
            f"Declaration file must have columns ID, REF, ALT (missing: "
            f"{', '.join(missing)}). Got: {list(declaration[0])}")
    out = variant_table(counts, declaration)
    out_dir = _run_dir(args.label or "variant_table")
    write_tsv(out_dir / "variants.tsv.gz", out, VARIANT_COLS)
    scored = [r for r in out if not math.isnan(r["log2FoldChange_expression"])]
    print(f"Variants:      {len(out):,}")
    print(f"  scored:      {len(scored):,}  (both REF and ALT quantified)")
    if scored:
        skews = sorted(abs(r["log2FoldChange_expression"]) for r in scored)
        print(f"  |skew|:      median {median(skews):.3f}, max {skews[-1]:.3f}")
    print(f"Saved:         {out_dir / 'variants.tsv.gz'}")
    return 0


def cmd_pipeline(args) -> int:
    """merge-counts per replicate → master-table → combine → variants."""
    setup_logging()
    reps = _parse_replicates(args.replicate)
    out_dir = _run_dir(args.label or "pipeline")
    assignment = list(read_tsv(args.assignment, header=False,
                               names=["barcode", "oligo", "quality", "support"]))
    summary: dict[str, Any] = {"replicates": [n for n, _ in reps],
                               "bc_threshold": args.threshold}

    per_rep: list[tuple[str, str]] = []
    per_rep_bc: list[tuple[str, str]] = []
    rep_stats = {}
    cols = ["oligo_name", "dna_counts", "rna_counts", "dna_normalized",
            "rna_normalized", "ratio", "log2FoldChange", "n_bc"]
    for name, path in reps:
        counts = read_tsv(path, header=False,
                          names=["barcode", "dna_count", "rna_count"])
        res = merge_counts(
            counts, iter(assignment), min_dna=args.min_dna_counts,
            min_rna=args.min_rna_counts,
            normalize_with_unassigned=args.include_unassigned_for_normalization,
            scaling=args.scaling, outlier_method=args.outlier_detection,
            mad_bins=args.outlier_mad_bins, mad_times=args.outlier_mad_times,
            zscore_times=args.outlier_zscore_times,
            mad_two_sided=args.mad_two_sided,
            mad_include_top_bin=args.mad_include_top_bin)
        p = out_dir / f"counts.{safe_label(name)}.tsv.gz"
        write_tsv(p, res["table"], cols)
        per_rep.append((name, str(p)))
        # Suffix the count columns with the replicate name so even this
        # per-replicate intermediate satisfies the IGVF
        # `reporter_experiment_barcode` schema (checked by
        # `igvfagent mpralib validate`). Unsuffixed dna_count/rna_count
        # are rejected by that schema's additionalProperties: false.
        pb = out_dir / f"barcodes.{safe_label(name)}.tsv.gz"
        rep_rows = [{"barcode": b["barcode"], "oligo_name": b["oligo_name"],
                     f"dna_count_{name}": b["dna_count"],
                     f"rna_count_{name}": b["rna_count"]}
                    for b in res["barcodes"]]
        write_tsv(pb, rep_rows, ["barcode", "oligo_name",
                                 f"dna_count_{name}", f"rna_count_{name}"])
        per_rep_bc.append((name, str(pb)))
        rep_stats[name] = res["statistic"]
    summary["per_replicate"] = rep_stats

    rows, filtered, stats = master_table(per_rep, threshold=args.threshold)
    write_tsv(out_dir / "master_table.all.tsv.gz", rows, MASTER_COLS)
    write_tsv(out_dir / "master_table.tsv.gz", filtered, MASTER_COLS)
    write_tsv(out_dir / "master_statistic.tsv.gz", stats,
              ["replicate", "mean_ratio", "mean_log2FoldChange",
               "mean_n_bc", "BC_filter"])
    summary["master_rows"] = len(rows)
    summary["master_rows_filtered"] = len(filtered)

    labels = None
    if args.labels:
        labels = {r["oligo_name"]: r["Label"] for r in
                  read_tsv(args.labels, header=False,
                           names=["oligo_name", "Label"])}
    combined = combine_replicates(filtered, labels)
    ccols = ["oligo_name", "replicates", "dna_counts", "rna_counts",
             "dna_normalized", "rna_normalized", "log2FoldChange",
             "mean_dna_counts", "mean_rna_counts", "mean_dna_normalized",
             "mean_rna_normalized", "mean_log2FoldChange", "mean_n_bc", "n_bc"]
    if labels:
        ccols.append("Label")
    write_tsv(out_dir / "combined.tsv.gz", combined, ccols)
    summary["oligos_combined"] = len(combined)

    allbc, filtbc, bcols = barcode_matrix(per_rep_bc, threshold=args.threshold)
    write_tsv(out_dir / "barcode_matrix.tsv.gz", filtbc, bcols)

    if args.declaration:
        # Upstream builds a variant table PER REPLICATE from the
        # per-replicate counts (those carry `ratio`), then pools them.
        declaration = list(read_tsv(args.declaration))
        per_rep_variants = []
        for name, path in per_rep:
            counts = [_coerce_count_row(r) for r in read_tsv(path)]
            vt = variant_table(counts, declaration)
            write_tsv(out_dir / f"variants.{safe_label(name)}.tsv.gz",
                      vt, VARIANT_COLS)
            per_rep_variants.append(vt)
        master_variants = master_variant_table(
            per_rep_variants, min_dna=args.min_dna_counts,
            min_rna=args.min_rna_counts)
        write_tsv(out_dir / "variants.tsv.gz", master_variants,
                  MASTER_VARIANT_COLS)
        summary["variants"] = len(declaration)
        summary["variants_scored"] = sum(
            1 for r in master_variants if r["log2FoldChange_expression"])

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Replicates:        {len(reps)}  ({', '.join(n for n, _ in reps)})")
    for name, s in rep_stats.items():
        print(f"  {name:12s} {s['barcodes dna/rna']:>9,} barcodes  "
              f"{s['% matched barcodes']:>5.1f}% assigned  "
              f"{s['oligos dna/rna']:>7,} oligos")
    print(f"Master table:      {len(rows):,} rows "
          f"({len(filtered):,} at n_bc >= {args.threshold})")
    print(f"Oligos combined:   {len(combined):,}")
    if args.declaration:
        print(f"Variants scored:   {summary['variants_scored']:,} "
              f"of {summary['variants']:,}")
    print(f"Saved:             {out_dir}")
    return 0


PLAYBOOK = """# MPRA count & assignment — skill playbook

Python port of the analysis core of
[MPRAsnakeflow](https://github.com/kircherlab/MPRAsnakeflow)
(Max Schubach, Berlin Institute of Health at Charité; MIT).
Cite Rosen *et al.*, *Genome Research* (2025), doi:10.1101/gr.281462.125.

`igvfagent oligo` designs the library; this skill processes what comes
back off the sequencer.

    design → ORDER → transfect → SEQUENCE → this skill → activity

## Why the two-step structure

An MPRA measures activity indirectly. Every oligo carries random
barcodes; you sequence DNA (input) and RNA (output) and the activity is
the RNA/DNA ratio. Two things have to happen first:

1. **Assignment** — a separate sequencing run tells you which barcode
   belongs to which oligo. A barcode is trusted only when at least
   `--minimum` reads agree and at least `--fraction` of them name the
   same oligo. Barcodes whose reads disagree are *ambiguous* and dropped,
   never guessed.
2. **Normalisation** — counts are summed over an oligo's barcodes,
   divided by the number of barcodes, and scaled per million, so an
   oligo with 40 barcodes is comparable to one with 11.

## Pipeline

| Step | Subcommand | In → out |
|---|---|---|
| Count barcodes | `count-bc` | FASTQ → barcode counts |
| Assign (no aligner) | `assign-exact` | reads + design → barcode↔oligo pairs |
| Trust the assignment | `assign-filter` | pairs → assignment |
| Quantify | `merge-counts` | barcode counts + assignment → per-oligo activity |
| Library complexity | `complexity` | barcode file → Lincoln-Petersen estimate |
| Stack replicates | `master-table` | per-replicate tables → master table |
| Aggregate | `combine-replicates` | master → one row per oligo |
| Per-barcode matrix | `barcode-matrix` | per-replicate barcodes → matrix |
| Allelic skew | `variant-table` | REF/ALT pairs → log2 skew |
| Everything | `pipeline` | assignment + replicates → all of the above |

## Usage

```bash
# Assignment, starting from aligner output (barcode-sorted triples)
igvfagent mpraflow assign-filter --pairs pairs.tsv.gz \\
    --minimum 3 --fraction 0.75

# Or skip the aligner entirely when inserts match the design exactly
igvfagent mpraflow assign-exact --design design.fa \\
    --reads inserts.fastq.gz --barcodes index.fastq.gz --barcode-length 15

# Quantify one replicate
igvfagent mpraflow merge-counts --counts rep1_counts.tsv.gz \\
    --assignment assignment.tsv.gz --min-rna-counts 1 \\
    --outlier-detection ratio_mad

# Everything, three replicates, with allelic skew
igvfagent mpraflow pipeline --assignment assignment.tsv.gz \\
    --replicate rep1=rep1.tsv.gz --replicate rep2=rep2.tsv.gz \\
    --replicate rep3=rep3.tsv.gz --declaration variants.tsv \\
    --threshold 10 --label my_experiment
```

## Outlier removal

`--outlier-detection ratio_mad` bins barcodes by RNA depth and trims
those whose DNA/RNA ratio deviates from their oligo's median by more
than `--outlier-mad-times` MADs. `rna_counts_zscore` instead drops
barcodes more than `--outlier-zscore-times` SDs from their oligo's mean
RNA count.

Two upstream quirks are **reproduced by default**, because they define
the numbers in published MPRAsnakeflow output:

| Quirk | Effect | Opt-in fix |
|---|---|---|
| MAD test is one-sided (`ratio_diff <= times*mad`) | deliberate — Rosen et al. define it as ratios "exceeding 5 log2-units *above* their oligo-specific median" | `--mad-two-sided` |
| Bin edges are `arange(0, n_bins)/n_bins`, stopping at the 95th percentile | barcodes above it get no bin, a NaN MAD, and are always dropped — the top 5% by RNA count | `--mad-include-top-bin` |

## Deviations from upstream

* No Snakemake, conda, R, pandas, numpy or pysam — the R master-table
  step and the pandas group-by/quantile logic are reimplemented on the
  standard library.
* **Alignment is not reimplemented.** Upstream offers five aligner
  backends (bwa, bbmap, pbmm2, exact, hybrid) to produce barcode↔oligo
  pairs. `assign-filter` starts from the barcode-sorted
  `barcode <tab> oligo <tab> quality` table they emit; `assign-exact`
  covers the alignment-free case directly. For bwa/bbmap runs, align
  with your existing tooling and feed the result to `assign-filter`.
* Upstream's ambiguous-barcode branch reports the quality of whichever
  oligo the loop happened to end on (an uninitialised-variable leak);
  this reports the modal oligo's quality.

## Reproducing the paper

`complexity` implements the Lincoln-Petersen estimate from Rosen *et al.*
(2025): each replicate is a capture of the barcode pool, barcodes seen in
two replicates are recaptures, and the gap between observed and estimated
barcodes is what the sequencing missed. Verified against the published
figures:

| Dataset | Metric | Paper | This skill |
|---|---|---|---|
| 8K-neurons | median assigned barcodes | 1,444,480 | 1,444,480 |
| 8K-neurons | median Lincoln index | 1,523,572 | 1,523,572 |
| 80K-neurons | median assigned barcodes | 5,459,247 | 5,459,247 |
| 80K-neurons | median Lincoln index | 6,243,618 | 6,243,618 |

The count-aggregation stage reproduces the published IGVF
`reporter experiment` artefact byte-for-byte (210,660/210,660 values). See
`Benchmarks/rosen2025_mprasnakeflow/`.

Note the portal's `reporter experiment` files were produced with
`--min-dna-counts 1`, not the tool default of 0. On the default a DNA
pseudocount of 1 per barcode is added, inflating `dna_counts` by `n_bc` per
oligo and shifting every `log2FoldChange`.

## Key outputs

Everything lands in `Docs/MPRASnakeflow/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `assignment.tsv.gz` | barcode → oligo, quality, read support |
| `counts.<rep>.tsv.gz` | per-oligo DNA/RNA counts, normalised, log2FC, n_bc |
| `master_table.tsv.gz` | all replicates stacked, filtered on n_bc |
| `combined.tsv.gz` | one row per oligo, pooled and per-replicate means |
| `barcode_matrix.tsv.gz` | per-barcode DNA/RNA across replicates |
| `variants.tsv.gz` | REF/ALT pairs with `log2FoldChange_expression` |
| `summary.json` | per-replicate and overall counts |
"""


def cmd_write_playbook(args) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = SKILL_DOC_DIR / "MPRA_SNAKEFLOW_SKILLS.md"
    out.write_text(PLAYBOOK)
    print(f"Wrote {out}")
    return 0


def _add_merge_args(p) -> None:
    p.add_argument("--min-dna-counts", type=int, default=0,
                   help="Minimum DNA counts per barcode (0 uses a pseudocount).")
    p.add_argument("--min-rna-counts", type=int, default=1,
                   help="Minimum RNA counts per barcode (0 uses a pseudocount).")
    p.add_argument("--include-unassigned-for-normalization",
                   action="store_true",
                   help="Keep unassigned barcodes in the normalisation totals.")
    p.add_argument("--scaling", type=float, default=DEFAULT_SCALING)
    p.add_argument("--outlier-detection", default=None,
                   choices=["ratio_mad", "rna_counts_zscore"])
    p.add_argument("--outlier-mad-bins", type=int, default=20)
    p.add_argument("--outlier-mad-times", type=float, default=5.0)
    p.add_argument("--outlier-zscore-times", type=float, default=3.0)
    p.add_argument("--mad-two-sided", action="store_true",
                   help="Compare |ratio_diff| (upstream trims one tail only).")
    p.add_argument("--mad-include-top-bin", action="store_true",
                   help="Close the top quantile bin (upstream always drops "
                        "the top 5%% of barcodes by RNA count).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="igvfagent mpraflow",
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("count-bc", help="FASTQ -> per-barcode counts.")
    p.add_argument("--fastq", required=True)
    p.add_argument("--barcode-start", type=int, default=0)
    p.add_argument("--barcode-length", type=int, default=None)
    p.add_argument("--min-quality", type=int, default=0,
                   help="Drop a read if any barcode base is below this Phred.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_count_bc)

    p = sub.add_parser("assign-exact",
                       help="Assignment reads -> barcode/oligo pairs (exact).")
    p.add_argument("--design", required=True, help="Design FASTA.")
    p.add_argument("--reads", required=True, help="Insert reads (FASTQ/FASTA).")
    p.add_argument("--barcodes", default=None,
                   help="Index read carrying the barcode; without it the "
                        "barcode is taken from the read name after ':'.")
    p.add_argument("--barcode-start", type=int, default=0)
    p.add_argument("--barcode-length", type=int, default=None)
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_assign_exact)

    p = sub.add_parser("assign-filter",
                       help="Barcode/oligo pairs -> trusted assignment.")
    p.add_argument("--pairs", required=True,
                   help="Barcode-sorted 'barcode<TAB>oligo<TAB>quality'.")
    p.add_argument("--minimum", type=int, default=DEFAULT_MIN_SUPPORT,
                   help="Minimum reads supporting an assignment.")
    p.add_argument("--fraction", type=float, default=DEFAULT_FRACTION,
                   help="Fraction of reads that must agree (> 0.5).")
    p.add_argument("--report-other", action="store_true")
    p.add_argument("--report-ambiguous", action="store_true")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_assign_filter)

    p = sub.add_parser("assign-stats", help="Assignment coverage statistics.")
    p.add_argument("--assignment", required=True)
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_assign_stats)

    p = sub.add_parser("merge-counts",
                       help="Barcode counts + assignment -> per-oligo activity.")
    p.add_argument("--counts", required=True,
                   help="'barcode<TAB>dna_count<TAB>rna_count'.")
    p.add_argument("--assignment", required=True)
    p.add_argument("--label", default=None)
    _add_merge_args(p)
    p.set_defaults(func=cmd_merge_counts)

    p = sub.add_parser("complexity",
                       help="Lincoln-Petersen library complexity from an "
                            "IGVF reporter-experiment-barcode file.")
    p.add_argument("--barcode-file", required=True,
                   help="IGVF 'reporter experiment barcode' TSV (wide format).")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_complexity)

    p = sub.add_parser("master-table", help="Stack per-replicate tables.")
    p.add_argument("--replicate", action="append", required=True,
                   metavar="NAME=PATH")
    p.add_argument("--threshold", type=int, default=DEFAULT_BC_THRESHOLD,
                   help="Minimum barcodes per oligo per replicate.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_master_table)

    p = sub.add_parser("combine-replicates",
                       help="Master table -> one row per oligo.")
    p.add_argument("--input", required=True, help="Master table TSV.")
    p.add_argument("--labels", default=None, help="oligo<TAB>label TSV.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_combine_replicates)

    p = sub.add_parser("barcode-matrix",
                       help="Per-barcode DNA/RNA matrix across replicates.")
    p.add_argument("--replicate", action="append", required=True,
                   metavar="NAME=PATH")
    p.add_argument("--threshold", type=int, default=DEFAULT_BC_THRESHOLD)
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_barcode_matrix)

    p = sub.add_parser("variant-table", help="REF/ALT pairs -> allelic skew.")
    p.add_argument("--counts", required=True, help="Per-oligo count table.")
    p.add_argument("--declaration", required=True,
                   help="TSV with columns ID, REF, ALT.")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_variant_table)

    p = sub.add_parser("pipeline",
                       help="merge-counts -> master -> combine -> variants.")
    p.add_argument("--assignment", required=True)
    p.add_argument("--replicate", action="append", required=True,
                   metavar="NAME=PATH")
    p.add_argument("--declaration", default=None)
    p.add_argument("--labels", default=None)
    p.add_argument("--threshold", type=int, default=DEFAULT_BC_THRESHOLD)
    p.add_argument("--label", default=None)
    _add_merge_args(p)
    p.set_defaults(func=cmd_pipeline)

    p = sub.add_parser("write-playbook",
                       help="Write Docs/Skills/MPRA_SNAKEFLOW_SKILLS.md.")
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
