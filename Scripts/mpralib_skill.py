"""MPRA barcode QC skill — outlier detection, activity, replicate agreement.

Python port of the analysis core of **MPRAlib** (Max Schubach, Berlin
Institute of Health at Charité; https://github.com/kircherlab/MPRAlib,
MIT licence). Upstream is an AnnData-backed library; this is a
dependency-free port of the parts that produce the published analyses.

If you publish results produced with this skill, cite the upstream work:

    Rosen JD, Vasanthakumari AD, Salomon K, de Lange N, Dash PM,
    Keukeleire P, Hassan A, Barrera A, Krupkin B, Oualline G, Kircher M,
    Love MI, Schubach M. Uniform processing and analysis of IGVF
    massively parallel reporter assay data with MPRAsnakeflow.
    Genome Research (2025). doi:10.1101/gr.281462.125

Where this sits
---------------
    oligo  →  ORDER  →  SEQUENCE  →  mpraflow  →  THIS SKILL
    design                          counts       barcode QC

``mpraflow`` turns reads into per-oligo activity. This skill asks whether
the individual barcodes behind those numbers can be trusted.

Why barcode outliers matter
---------------------------
Each oligo's activity is an average over its barcodes. A handful of
barcodes with wild RNA counts — a PCR jackpot, an unlucky integration
site — can drag an oligo's ratio far from where its other twenty barcodes
sit. Rosen et al. (2025) define three detectors, each catching a
different failure:

* **Global** — a barcode whose RNA count is extreme for the *whole
  replicate* (|z| > 3). Catches library-wide artefacts.
* **Oligo-specific** — a barcode whose RNA count is extreme *for its own
  oligo* (|z| > 3 within the oligo). Catches a single bad barcode among
  otherwise consistent siblings.
* **Large expression** — a barcode whose activity sits more than 5
  log2-units *above* its oligo's median. One-sided by design: this looks
  for runaway expression, not silence.

Their finding was that removal changes downstream variant calls very
little (correlations > 0.95), so these are diagnostics rather than a
mandatory cleaning step. Treat a high rate as a signal about the
experiment, not as a reason to filter by default.

Upstream layout → this module
-----------------------------
  =========================================  ==========================
  MPRAlib                                    here
  utils/file_validation.validate_tsv_...     ``validate``
  mpralib/schemas/*.json                     ``_mpra_schemas.SCHEMAS``
  _barcode_filter_global_outliers            ``outliers --method global``
  _barcode_filter_oligo_specific_outliers    ``outliers --method oligo``
  _barcode_filter_large_expression_outliers  ``outliers --method large_expression``
  MPRAData._normalize / _compute_activities  ``activity``
  (Figure 4B, replicate agreement)           ``consistency``

Deliberate deviations from upstream
-----------------------------------
* **No AnnData, numpy or pandas.** Upstream materialises a
  barcodes × replicates matrix; the 240K libraries have 20 M barcodes,
  where that costs several GB. This streams the file instead, matching
  how ``mpra_snakeflow_skill`` already works.
* **No ``jsonschema`` dependency.** The eight IGVF schemas use a small,
  fixed subset of draft-07 (type / properties / patternProperties /
  additionalProperties / required / enum / pattern / minLength /
  maxLength / minimum / items / minItems / minProperties / anyOf), so
  that subset is validated directly. The schemas themselves are copied
  verbatim — they define an interchange standard, and a paraphrase would
  be worse than useless.
* pandas' ``std()`` is ddof=1 and skips NaN; both are reproduced
  explicitly, because using the population SD would shift every z-score.

Subcommands
-----------
    outliers      Flag barcode outliers by one of the three methods
    consistency   How often the same barcode is an outlier in every replicate
    activity      Per-barcode normalised DNA/RNA activity (log2 ratio)
    validate      Check a file against an IGVF MPRA standard schema
    schemas       List the standard formats and their required columns
    write-playbook  Write Docs/Skills/MPRALIB_SKILLS.md
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
from typing import Any, Iterator

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "MPRAlib"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

# MPRAlib class defaults (mpradata.py: _SCALING, _PSEUDOCOUNT).
SCALING = 1e6
PSEUDOCOUNT = 1
DEFAULT_ZSCORE = 3.0
DEFAULT_ACTIVITY_TIMES = 5.0

METHODS = ("global", "oligo", "large_expression")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"mpralib_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def _int(v: str) -> int:
    v = v.strip()
    if v in ("", "NA", "NaN", "nan", "."):
        return 0
    return int(float(v))


class BarcodeTable:
    """Streaming reader for an IGVF 'reporter experiment barcode' file.

    One row per barcode, with ``dna_count_<r>`` / ``rna_count_<r>`` per
    replicate and empty cells where a replicate did not see the barcode.
    """

    def __init__(self, path: str | Path):
        self.path = str(path)
        with opener(self.path) as fh:
            cols = fh.readline().rstrip("\r\n").split("\t")
        self.cols = cols
        self.idx = {c: i for i, c in enumerate(cols)}
        self.replicates = sorted({c[len("dna_count_"):] for c in cols
                                  if c.startswith("dna_count_")})
        if not self.replicates:
            raise SystemExit(
                f"No dna_count_<replicate> columns in {path}. Expected an "
                "IGVF 'reporter experiment barcode' file.")
        if "barcode" not in self.idx or "oligo_name" not in self.idx:
            raise SystemExit(
                f"{path} must have 'barcode' and 'oligo_name' columns.")

    def rows(self) -> Iterator[tuple[str, str, list[int], list[int]]]:
        """Yield ``(barcode, oligo, dna_per_replicate, rna_per_replicate)``."""
        di = [self.idx[f"dna_count_{r}"] for r in self.replicates]
        ri = [self.idx[f"rna_count_{r}"] for r in self.replicates]
        b, o = self.idx["barcode"], self.idx["oligo_name"]
        ncol = len(self.cols)
        with opener(self.path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if len(f) < ncol:
                    continue
                yield (f[b], f[o],
                       [_int(f[i]) for i in di],
                       [_int(f[i]) for i in ri])


# ─── Pass 1: per-replicate totals and per-oligo moments ────────────────────

def _welford_std(n: int, mean: float, m2: float) -> float:
    """Sample standard deviation (ddof=1) from streaming moments.

    pandas' ``std()`` defaults to ddof=1; using the population SD instead
    would shift every z-score, so the sample form is explicit here.
    """
    if n < 2:
        return float("nan")
    return math.sqrt(m2 / (n - 1))


class _Moments:
    """Streaming mean / variance (Welford), NaN-skipping by construction."""

    __slots__ = ("n", "mean", "m2")

    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0

    def add(self, x: float) -> None:
        self.n += 1
        d = x - self.mean
        self.mean += d / self.n
        self.m2 += d * (x - self.mean)

    @property
    def std(self) -> float:
        return _welford_std(self.n, self.mean, self.m2)


def scan_totals(table: BarcodeTable, *, need_oligo_moments: bool
                ) -> dict[str, Any]:
    """One pass: replicate totals, observed counts, per-oligo RNA moments."""
    reps = table.replicates
    k = len(reps)
    dna_total = [0] * k
    rna_total = [0] * k
    observed_n = [0] * k
    global_moments = [_Moments() for _ in range(k)]
    oligo_moments: dict[str, list[_Moments]] = {}
    n_rows = 0

    for _, oligo, dna, rna in table.rows():
        n_rows += 1
        om = None
        if need_oligo_moments:
            om = oligo_moments.get(oligo)
            if om is None:
                om = [_Moments() for _ in range(k)]
                oligo_moments[oligo] = om
        for i in range(k):
            d, r = dna[i], rna[i]
            dna_total[i] += d
            rna_total[i] += r
            # MPRAlib's `observed`: (dna + rna) > 0 for that replicate.
            if d + r > 0:
                observed_n[i] += 1
                global_moments[i].add(float(r))
                if om is not None:
                    om[i].add(float(r))
    return {
        "replicates": reps,
        "dna_total": dna_total,
        "rna_total": rna_total,
        "observed_n": observed_n,
        "global": global_moments,
        "oligo": oligo_moments,
        "rows": n_rows,
    }


# ─── The three detectors (port of MPRAlib._barcode_filter_*) ───────────────

def detect_global_outliers(table: BarcodeTable, totals: dict,
                           times_zscore: float = DEFAULT_ZSCORE
                           ) -> tuple[dict[str, set[str]], dict[str, int]]:
    """A barcode whose RNA count is extreme for the whole replicate.

    Port of ``_barcode_filter_global_outliers``: per replicate, z-score
    the RNA counts of every *observed* barcode and flag ``|z| > 3``.
    Non-observed entries are NaN upstream, and ``NaN > z`` is False, so
    they are never flagged.
    """
    reps = totals["replicates"]
    k = len(reps)
    mean = [m.mean for m in totals["global"]]
    std = [m.std for m in totals["global"]]
    flagged: dict[str, set[str]] = {r: set() for r in reps}
    for bc, _, dna, rna in table.rows():
        for i in range(k):
            if dna[i] + rna[i] <= 0:
                continue
            s = std[i]
            if s is None or math.isnan(s) or s == 0:
                continue
            if abs((rna[i] - mean[i]) / s) > times_zscore:
                flagged[reps[i]].add(bc)
    return flagged, {r: totals["observed_n"][i] for i, r in enumerate(reps)}


def detect_oligo_outliers(table: BarcodeTable, totals: dict,
                          times_zscore: float = DEFAULT_ZSCORE
                          ) -> tuple[dict[str, set[str]], dict[str, int]]:
    """A barcode whose RNA count is extreme *for its own oligo*.

    Port of ``_barcode_filter_oligo_specific_outliers``. Upstream does
    ``transform("std").fillna(0).replace(0, 1)``: an oligo with a single
    observed barcode (std NaN) or no variance (std 0) gets std = 1, which
    makes the z-score the raw deviation rather than dividing by zero.
    """
    reps = totals["replicates"]
    k = len(reps)
    om = totals["oligo"]
    flagged: dict[str, set[str]] = {r: set() for r in reps}
    for bc, oligo, dna, rna in table.rows():
        moments = om.get(oligo)
        if moments is None:
            continue
        for i in range(k):
            if dna[i] + rna[i] <= 0:
                continue
            m = moments[i]
            s = m.std
            if s is None or math.isnan(s) or s == 0:
                s = 1.0            # upstream: fillna(0).replace(0, 1)
            if abs((rna[i] - m.mean) / s) > times_zscore:
                flagged[reps[i]].add(bc)
    return flagged, {r: totals["observed_n"][i] for i, r in enumerate(reps)}


def detect_large_expression_outliers(table: BarcodeTable, totals: dict,
                                     times_activity: float = DEFAULT_ACTIVITY_TIMES
                                     ) -> tuple[set[str], int, dict[str, float]]:
    """A barcode whose activity runs far above its oligo's median.

    Port of ``_barcode_filter_large_expression_outliers``. Activity is
    summed across replicates *after* CPM normalisation, so a barcode is
    judged on its whole-experiment behaviour, and the test is one-sided
    (``diff > times``) — this looks for runaway expression, not silence.
    Unlike the other two detectors this yields one verdict per barcode,
    applied to every replicate.
    """
    reps = totals["replicates"]
    k = len(reps)
    # MPRAlib: total = raw total + pseudocount * (observed barcodes)
    denom_d = [totals["dna_total"][i] + PSEUDOCOUNT * totals["observed_n"][i]
               for i in range(k)]
    denom_r = [totals["rna_total"][i] + PSEUDOCOUNT * totals["observed_n"][i]
               for i in range(k)]
    denom_d = [d if d else 1 for d in denom_d]
    denom_r = [d if d else 1 for d in denom_r]

    # Pass A: per-barcode log2 ratio, and per-oligo lists for the median.
    ratios: dict[str, float] = {}
    by_oligo: dict[str, list[float]] = defaultdict(list)
    for bc, oligo, dna, rna in table.rows():
        nd = nr = 0.0
        seen = False
        for i in range(k):
            if dna[i] + rna[i] <= 0:
                continue          # normalisation multiplies by `observed`
            seen = True
            nd += (dna[i] + PSEUDOCOUNT) / denom_d[i] * SCALING
            nr += (rna[i] + PSEUDOCOUNT) / denom_r[i] * SCALING
        if not seen:
            continue              # never observed -> NaN upstream
        if nd == 0 or nr == 0:
            continue              # log2(0) -> -inf -> NaN upstream
        val = math.log2(nr / nd)
        ratios[bc] = val
        by_oligo[oligo].append(val)

    medians = {o: statistics.median(v) for o, v in by_oligo.items() if v}

    # Pass B: flag barcodes above their oligo's median by more than `times`.
    flagged: set[str] = set()
    for bc, oligo, dna, rna in table.rows():
        val = ratios.get(bc)
        if val is None:
            continue
        med = medians.get(oligo)
        if med is None:
            continue
        if (val - med) > times_activity:
            flagged.add(bc)
    considered = len(ratios)
    return flagged, considered, medians


# ─── Commands ───────────────────────────────────────────────────────────────

def _rates(flagged: dict[str, set[str]], denom: dict[str, int]) -> dict[str, float]:
    return {r: (100.0 * len(flagged[r]) / denom[r] if denom.get(r) else 0.0)
            for r in flagged}


def run_outliers(path: str, method: str, *, times_zscore: float,
                 times_activity: float) -> dict[str, Any]:
    table = BarcodeTable(path)
    totals = scan_totals(table, need_oligo_moments=(method == "oligo"))
    reps = table.replicates
    out: dict[str, Any] = {"file": str(path), "method": method,
                           "replicates": reps, "barcodes": totals["rows"],
                           "observed_per_replicate": {
                               r: totals["observed_n"][i]
                               for i, r in enumerate(reps)}}

    if method == "large_expression":
        flagged, considered, _ = detect_large_expression_outliers(
            table, totals, times_activity)
        out["flagged_barcodes"] = len(flagged)
        out["considered_barcodes"] = considered
        out["percent"] = 100.0 * len(flagged) / considered if considered else 0.0
        out["per_replicate_percent"] = {r: out["percent"] for r in reps}
        out["_flagged"] = flagged
        return out

    if method == "global":
        flagged, denom = detect_global_outliers(table, totals, times_zscore)
    else:
        flagged, denom = detect_oligo_outliers(table, totals, times_zscore)
    rates = _rates(flagged, denom)
    out["flagged_per_replicate"] = {r: len(flagged[r]) for r in reps}
    out["per_replicate_percent"] = rates
    # The paper reports "the average among the replicates" for these two.
    out["percent"] = sum(rates.values()) / len(rates) if rates else 0.0
    out["_flagged"] = flagged
    return out


def cmd_outliers(args) -> int:
    setup_logging()
    res = run_outliers(args.barcode_file, args.method,
                       times_zscore=args.times_zscore,
                       times_activity=args.times_activity)
    flagged = res.pop("_flagged")
    out_dir = _run_dir(args.label or f"outliers_{args.method}")
    (out_dir / "outliers.json").write_text(json.dumps(res, indent=2))
    with opener(out_dir / "outlier_barcodes.tsv.gz", "wt") as fh:
        fh.write("replicate\tbarcode\n")
        if isinstance(flagged, set):
            for bc in sorted(flagged):
                fh.write(f"all\t{bc}\n")
        else:
            for r, s in flagged.items():
                for bc in sorted(s):
                    fh.write(f"{r}\t{bc}\n")

    print(f"Barcodes:          {res['barcodes']:,}")
    print(f"Replicates:        {len(res['replicates'])}  "
          f"({', '.join(res['replicates'])})")
    print(f"Method:            {args.method}")
    if args.method == "large_expression":
        print(f"  considered:      {res['considered_barcodes']:,}")
        print(f"  flagged:         {res['flagged_barcodes']:,}  "
              f"({res['percent']:.2f}%)")
    else:
        for r in res["replicates"]:
            print(f"  {r:10s} {res['flagged_per_replicate'][r]:>9,}  "
                  f"({res['per_replicate_percent'][r]:.2f}%)")
        print(f"  mean across replicates: {res['percent']:.2f}%")
    print(f"Saved:             {out_dir}")
    return 0


def cmd_consistency(args) -> int:
    """How often the same barcode is flagged in EVERY replicate.

    Figure 4B of Rosen et al.: a detector that fires on the same barcodes
    in every replicate is describing the library; one that fires on
    different barcodes each time is describing noise.
    """
    setup_logging()
    res = run_outliers(args.barcode_file, args.method,
                       times_zscore=args.times_zscore,
                       times_activity=args.times_activity)
    flagged = res.pop("_flagged")
    if isinstance(flagged, set):
        raise SystemExit(
            "large_expression yields one verdict per barcode across all "
            "replicates, so replicate consistency is not defined for it. "
            "Use --method global or --method oligo.")
    reps = res["replicates"]
    union: set[str] = set()
    for s in flagged.values():
        union |= s
    inter: set[str] = set(flagged[reps[0]])
    for r in reps[1:]:
        inter &= flagged[r]
    # Rosen et al. report "the percentage of barcodes that were
    # consistently identified as outliers across all three biological
    # replicates" — of the barcodes a replicate typically flags, how many
    # are flagged EVERY time. The denominator is therefore the mean
    # per-replicate count, not the union. (Over the union the same data
    # reads 3.0% where the paper reports 7.8%.)
    mean_per_rep = statistics.fmean([len(flagged[r]) for r in reps])
    pct = 100.0 * len(inter) / mean_per_rep if mean_per_rep else 0.0
    pct_union = 100.0 * len(inter) / len(union) if union else 0.0

    out = {"file": str(args.barcode_file), "method": args.method,
           "replicates": reps,
           "flagged_per_replicate": {r: len(flagged[r]) for r in reps},
           "flagged_in_any": len(union),
           "flagged_in_all": len(inter),
           "mean_flagged_per_replicate": mean_per_rep,
           "percent_consistent": pct,
           "percent_consistent_over_union": pct_union}
    out_dir = _run_dir(args.label or f"consistency_{args.method}")
    (out_dir / "consistency.json").write_text(json.dumps(out, indent=2))
    print(f"Method:              {args.method}")
    for r in reps:
        print(f"  {r:10s} {len(flagged[r]):>9,} outlier barcodes")
    print(f"Flagged in any:      {len(union):,}")
    print(f"Flagged in ALL:      {len(inter):,}")
    print(f"Replicate consistency: {pct:.1f}%   "
          f"(flagged-in-all / mean per replicate, as in Rosen et al.)")
    print(f"  over the union:    {pct_union:.1f}%")
    print(f"Saved:               {out_dir}")
    return 0


def cmd_activity(args) -> int:
    """Per-barcode normalised activity, log2(RNA_cpm / DNA_cpm)."""
    setup_logging()
    table = BarcodeTable(args.barcode_file)
    totals = scan_totals(table, need_oligo_moments=False)
    reps = table.replicates
    k = len(reps)
    denom_d = [totals["dna_total"][i] + PSEUDOCOUNT * totals["observed_n"][i]
               for i in range(k)]
    denom_r = [totals["rna_total"][i] + PSEUDOCOUNT * totals["observed_n"][i]
               for i in range(k)]
    denom_d = [d if d else 1 for d in denom_d]
    denom_r = [d if d else 1 for d in denom_r]

    out_dir = _run_dir(args.label or "activity")
    out = out_dir / "barcode_activity.tsv.gz"
    n = 0
    with opener(out, "wt") as fh:
        fh.write("barcode\toligo_name\t"
                 + "\t".join(f"activity_{r}" for r in reps) + "\n")
        for bc, oligo, dna, rna in table.rows():
            vals = []
            for i in range(k):
                if dna[i] + rna[i] <= 0:
                    vals.append("NA")
                    continue
                nd = (dna[i] + PSEUDOCOUNT) / denom_d[i] * SCALING
                nr = (rna[i] + PSEUDOCOUNT) / denom_r[i] * SCALING
                vals.append(f"{math.log2(nr / nd):.6f}" if nd and nr else "NA")
            fh.write(f"{bc}\t{oligo}\t" + "\t".join(vals) + "\n")
            n += 1
    print(f"Barcodes:  {n:,}")
    print(f"Saved:     {out}")
    return 0


PLAYBOOK = """# MPRA barcode QC (MPRAlib) — skill playbook

Python port of the analysis core of
[MPRAlib](https://github.com/kircherlab/MPRAlib) (Max Schubach, Berlin
Institute of Health at Charité; MIT). Cite Rosen *et al.*,
*Genome Research* (2025), doi:10.1101/gr.281462.125.

    oligo  →  ORDER  →  SEQUENCE  →  mpraflow  →  this skill
    design                          counts       barcode QC

## Why barcode outliers matter

An oligo's activity is an average over its barcodes. A few barcodes with
wild RNA counts — a PCR jackpot, an unlucky integration site — can drag
the oligo's ratio away from where its other twenty barcodes sit. The
paper defines three detectors, each catching a different failure:

| Method | Question | Rule |
|---|---|---|
| `global` | extreme for the whole replicate? | \\|z\\| > 3 on RNA counts across all barcodes |
| `oligo` | extreme for its own oligo? | \\|z\\| > 3 on RNA counts within the oligo |
| `large_expression` | runaway activity? | log2 ratio > 5 above the oligo's median |

`large_expression` is one-sided by design — it looks for runaway
expression, not silence — and yields one verdict per barcode rather than
one per replicate.

## Should you filter on these?

Usually **no**. Rosen *et al.* compared variant effects before and after
removal and found correlations above 0.95, concluding that "explicit
outlier removal is not necessary for variant effect analysis when using
BCalm". Treat a high rate as information about the experiment, not as a
cleaning step to apply by reflex.

What the rates *do* tell you (paper, Figure 4A):

| | lentiviral | episomal (plasmid) |
|---|---|---|
| global | 0.67–1.56% | 0.12–0.48% |
| oligo-specific | 0.85–1.86% | 0.85–1.86% |
| large expression | up to 0.01% | 0.05–0.53% |

## Replicate consistency

`consistency` reports what fraction of outlier barcodes are flagged in
*every* replicate. A detector firing on the same barcodes each time is
describing the library; one firing on different barcodes is describing
noise. The paper found episomal assays far more consistent (63.2–83.9%)
than lentiviral ones (5.6–49.0%), and differentiated cells more
consistent than progenitors.

## Usage

```bash
# Outlier rates, one method at a time
igvfagent mpralib outliers --barcode-file reporter_experiment.barcode.tsv.gz \\
    --method global
igvfagent mpralib outliers --barcode-file reporter_experiment.barcode.tsv.gz \\
    --method large_expression

# How reproducible are those calls across replicates?
igvfagent mpralib consistency \\
    --barcode-file reporter_experiment.barcode.tsv.gz --method global

# Per-barcode normalised activity
igvfagent mpralib activity --barcode-file reporter_experiment.barcode.tsv.gz
```

Input is the IGVF **reporter experiment barcode** file, the same one
`mpraflow` consumes — downloadable straight from the portal.

## Validating IGVF standard formats

The IGVF MPRA focus group defined eight interchange formats
(Supplementary Note S1). `validate` checks a file against one and names
the offending column and reason for every bad row:

```bash
igvfagent mpralib schemas          # list the formats + required columns
igvfagent mpralib validate --file master_table.tsv.gz \
    --schema reporter_experiment
```

Every file `igvfagent mpraflow` writes passes: `master_table*.tsv.gz`
against `reporter_experiment`, and `barcode_matrix.tsv.gz` plus the
per-replicate `barcodes.<rep>.tsv.gz` against
`reporter_experiment_barcode`. Published IGVF portal files pass too, so
the check runs both ways — use it before submitting files to the portal,
and to confirm a file someone sent you is the format it claims to be.

One subtlety worth knowing: the count columns are
`anyOf: [integer, string with maxLength 0]` — "a count, or blank if that
replicate never saw the barcode". Blank cells are legal, and 7-9% of rows
in real portal files use them.

## Deviations from upstream

* No AnnData / numpy / pandas. Upstream materialises a
  barcodes × replicates matrix; the 240K libraries have 20 M barcodes,
  where that costs several GB. This streams the file instead.
* pandas' `std()` is ddof=1 and skips NaN — both reproduced explicitly,
  since the population SD would shift every z-score.
* Upstream's `transform("std").fillna(0).replace(0, 1)` for
  oligo-specific z-scores (single-barcode or zero-variance oligos get
  std = 1) is reproduced exactly.
"""


def cmd_write_playbook(args) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = SKILL_DOC_DIR / "MPRALIB_SKILLS.md"
    out.write_text(PLAYBOOK)
    print(f"Wrote {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="igvfagent mpralib",
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")

    def _common(p):
        p.add_argument("--barcode-file", required=True,
                       help="IGVF 'reporter experiment barcode' TSV.")
        p.add_argument("--times-zscore", type=float, default=DEFAULT_ZSCORE)
        p.add_argument("--times-activity", type=float,
                       default=DEFAULT_ACTIVITY_TIMES)
        p.add_argument("--label", default=None)

    p = sub.add_parser("outliers", help="Flag barcode outliers.")
    p.add_argument("--method", choices=METHODS, default="global")
    _common(p)
    p.set_defaults(func=cmd_outliers)

    p = sub.add_parser("consistency",
                       help="Fraction of outliers flagged in every replicate.")
    p.add_argument("--method", choices=("global", "oligo"), default="global")
    _common(p)
    p.set_defaults(func=cmd_consistency)

    p = sub.add_parser("activity",
                       help="Per-barcode normalised log2(RNA/DNA).")
    p.add_argument("--barcode-file", required=True)
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_activity)

    p = sub.add_parser("validate",
                       help="Check a file against an IGVF MPRA standard schema.")
    p.add_argument("--file", required=True, help="TSV/BED to validate (.gz ok).")
    p.add_argument("--schema", required=True,
                   help="Schema name; see `mpralib schemas`.")
    p.add_argument("--max-errors", type=int, default=20,
                   help="Stop collecting detail after this many bad rows.")
    p.add_argument("--max-rows", type=int, default=0,
                   help="Validate only the first N rows (0 = all).")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("schemas",
                       help="List the IGVF MPRA standard formats.")
    p.set_defaults(func=cmd_schemas)

    p = sub.add_parser("write-playbook",
                       help="Write Docs/Skills/MPRALIB_SKILLS.md.")
    p.set_defaults(func=cmd_write_playbook)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 2
    return int(args.func(args) or 0)



# ─── IGVF standard-format validation ───────────────────────────────────────
#
# Port of MPRAlib's utils/file_validation.py. Upstream leans on the
# `jsonschema` package; the eight IGVF schemas use only a small fixed
# subset of draft-07, so that subset is checked directly here rather than
# taking a dependency for it.
#
# TSV carries no types — every cell arrives as a string — so values are
# coerced to the type the schema declares before checking, exactly as
# upstream's _convert_row_value does. A cell that will not convert is
# left as a string so validation reports it, rather than crashing.

# Formats stored without a header row (BED-derived), so columns are
# positional. Same lists as upstream's _get_header_for_schema.
POSITIONAL_HEADERS = {
    "reporter_barcode_to_element_mapping": ["barcode", "oligoName"],
    "reporter_genomic_element": [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand",
        "log2FoldChange", "inputCount", "outputCount",
        "minusLog10PValue", "minusLog10QValue"],
    "reporter_genomic_variant": [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand",
        "log2FoldChange", "inputCountRef", "outputCountRef",
        "inputCountAlt", "outputCountAlt", "minusLog10PValue",
        "minusLog10QValue", "postProbEffect", "CI_lower_95", "CI_upper_95",
        "variantPos", "refAllele", "altAllele"],
}


def _coerce(value: str, spec: dict) -> Any:
    """TSV string → the type the schema declares (upstream _convert_row_value)."""
    t = spec.get("type")
    try:
        if t == "integer":
            return int(value)
        if t == "number":
            return float(value)
        if t == "array":
            import ast
            return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value        # leave it; validation will report the type error
    return value


def _spec_for(key: str, schema: dict) -> dict | None:
    props = schema.get("properties", {})
    if key in props:
        return props[key]
    import re as _re
    for pat, spec in schema.get("patternProperties", {}).items():
        if _re.match(pat, key):
            return spec
    return None


def _check_value(key: str, raw: Any, spec: dict) -> list[str]:
    """Validate one *raw* TSV value against one (sub)schema.

    Coercion happens per branch, not once up front. An ``anyOf`` of
    ``{integer}`` / ``{string, maxLength: 0}`` — the IGVF way of saying
    "a count, or blank if this replicate never saw the barcode" — only
    passes if the integer branch is allowed to parse the string first.
    """
    if "anyOf" in spec:
        errs_all: list[str] = []
        for sub in spec["anyOf"]:
            errs = _check_value(key, raw, sub)
            if not errs:
                return []
            errs_all += errs
        return [f"{key}: matches none of the allowed forms ({'; '.join(errs_all[:2])})"]

    value = _coerce(raw, spec) if isinstance(raw, str) else raw
    t = spec.get("type")
    if t == "integer" and not isinstance(value, int):
        return [f"{key}: expected integer, got {value!r}"]
    if t == "number" and not isinstance(value, (int, float)):
        return [f"{key}: expected number, got {value!r}"]
    if t == "string" and not isinstance(value, str):
        return [f"{key}: expected string, got {value!r}"]
    if t == "array" and not isinstance(value, (list, tuple)):
        return [f"{key}: expected array, got {value!r}"]
    if t == "object" and not isinstance(value, dict):
        return [f"{key}: expected object, got {value!r}"]

    errs: list[str] = []
    if "enum" in spec and value not in spec["enum"]:
        errs.append(f"{key}: {value!r} not one of {spec['enum']}")
    if "pattern" in spec and isinstance(value, str):
        import re as _re
        if not _re.search(spec["pattern"], value):
            errs.append(f"{key}: {value!r} does not match /{spec['pattern']}/")
    if "minLength" in spec and isinstance(value, str) \
            and len(value) < spec["minLength"]:
        errs.append(f"{key}: shorter than {spec['minLength']}")
    if "maxLength" in spec and isinstance(value, str) \
            and len(value) > spec["maxLength"]:
        errs.append(f"{key}: longer than {spec['maxLength']}")
    if "minimum" in spec and isinstance(value, (int, float)) \
            and value < spec["minimum"]:
        errs.append(f"{key}: {value} below minimum {spec['minimum']}")
    if "minItems" in spec and isinstance(value, (list, tuple)) \
            and len(value) < spec["minItems"]:
        errs.append(f"{key}: fewer than {spec['minItems']} items")
    if "items" in spec and isinstance(value, (list, tuple)):
        for i, item in enumerate(value):
            errs += _check_value(f"{key}[{i}]", item, spec["items"])
    return errs


def validate_row(row: dict, schema: dict) -> list[str]:
    """Validate one parsed TSV row. Returns a list of error strings."""
    errs: list[str] = []
    for req in schema.get("required", []):
        if req not in row or row[req] in (None, ""):
            errs.append(f"missing required column {req!r}")
    if schema.get("additionalProperties") is False:
        known = set(schema.get("properties", {}))
        patterns = list(schema.get("patternProperties", {}))
        import re as _re
        for key in row:
            if key in known:
                continue
            if any(_re.match(p, key) for p in patterns):
                continue
            errs.append(f"unexpected column {key!r}")
    minprops = schema.get("minProperties")
    if minprops is not None and len(row) < minprops:
        errs.append(f"row has {len(row)} columns, schema requires >= {minprops}")
    for key, value in row.items():
        if value in (None, ""):
            continue
        spec = _spec_for(key, schema)
        if spec is None:
            continue
        errs += _check_value(key, value, spec)
    return errs


def validate_file(path: str | Path, schema_name: str, *,
                  max_errors: int = 20, max_rows: int = 0) -> dict[str, Any]:
    """Validate a TSV against one IGVF MPRA standard schema."""
    from _mpra_schemas import SCHEMAS
    if schema_name not in SCHEMAS:
        raise SystemExit(
            f"Unknown schema {schema_name!r}. Choose one of: "
            + ", ".join(SCHEMAS))
    schema = SCHEMAS[schema_name]
    header = POSITIONAL_HEADERS.get(schema_name)

    errors: list[dict] = []
    n_rows = 0
    n_bad = 0
    with opener(path) as fh:
        if header is None:
            first = fh.readline().rstrip("\r\n")
            if not first:
                return {"file": str(path), "schema": schema_name,
                        "rows": 0, "rows_invalid": 0, "valid": False,
                        "errors": [{"row": 0, "problems": ["file is empty"]}]}
            header = first.split("\t")
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            n_rows += 1
            row = dict(zip(header, line.split("\t")))
            problems = validate_row(row, schema)
            if problems:
                n_bad += 1
                if len(errors) < max_errors:
                    errors.append({"row": n_rows, "problems": problems})
            if max_rows and n_rows >= max_rows:
                break
    return {"file": str(path), "schema": schema_name, "columns": header,
            "rows": n_rows, "rows_invalid": n_bad,
            "valid": n_bad == 0 and n_rows > 0, "errors": errors}


def cmd_validate(args) -> int:
    setup_logging()
    res = validate_file(args.file, args.schema,
                        max_errors=args.max_errors, max_rows=args.max_rows)
    out_dir = _run_dir(args.label or f"validate_{args.schema}")
    (out_dir / "validation.json").write_text(json.dumps(res, indent=2))
    print(f"File:     {res['file']}")
    print(f"Schema:   {res['schema']}")
    print(f"Columns:  {', '.join(res['columns'])}")
    print(f"Rows:     {res['rows']:,}")
    if res["valid"]:
        print("Result:   VALID — every row matches the IGVF standard")
    else:
        print(f"Result:   INVALID — {res['rows_invalid']:,} row(s) failed")
        for e in res["errors"]:
            print(f"  row {e['row']}: {'; '.join(e['problems'][:3])}")
        if res["rows_invalid"] > len(res["errors"]):
            print(f"  … and {res['rows_invalid'] - len(res['errors']):,} more")
    print(f"Saved:    {out_dir / 'validation.json'}")
    return 0 if res["valid"] else 1


def cmd_schemas(args) -> int:
    from _mpra_schemas import SCHEMAS
    for name, schema in SCHEMAS.items():
        req = schema.get("required", [])
        print(f"\n{name}")
        print(f"  {schema.get('description', schema.get('title', ''))[:96]}")
        if name in POSITIONAL_HEADERS:
            print(f"  no header row; columns are positional: "
                  f"{', '.join(POSITIONAL_HEADERS[name])}")
        if req:
            print(f"  required: {', '.join(req)}")
        if schema.get("patternProperties"):
            print(f"  pattern columns: "
                  f"{', '.join(schema['patternProperties'])}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
