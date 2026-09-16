#!/usr/bin/env python3
"""scQers — single-cell quantitative expression reporters.

Implementation of the scQer method (Shendure lab; "Multiplex profiling of
developmental enhancers with quantitative, single-cell expression reporters",
shendurelab/scQers, MIT). The upstream repository is R and shell, explicitly
shared "for transparency" rather than as a pipeline -- its own README says so.
This is the method as a runnable Python CLI, which is what IGVFagent needs.

WHAT THE ASSAY MEASURES, because the two barcodes are the whole trick.

A candidate cis-regulatory element (CRE) drives mCherry carrying an **mBC**
(mCherry barcode). The same construct constitutively expresses an **oBC**
(oligo barcode) from a separate promoter. In single cells you therefore read:

    oBC  -> WHICH element this cell received   (constitutive: presence)
    mBC  -> HOW HARD that element is driving   (CRE-dependent: activity)

So a cell with a strong enhancer and a cell with a dead one both report their
oBC; only the mBC count separates them. That is what makes the readout
*quantitative* rather than a sort-based enrichment, and it is why every step
below keeps the two barcodes distinguishable.

ACTIVITY is mean mBC UMIs per integration, normalised for how deeply that cell
was sequenced:

    normalised mBC UMI = UMIs_mBC / gex_UMI * mean(gex_UMI)

Without that normalisation a deeply-sequenced cluster looks like a strong
enhancer. The absolute number is meaningless on its own, which is why the
library carries two controls on every experiment -- `minP` (minimal promoter,
the floor for "a promoter with no enhancer") and `noP` (no promoter, the floor
for "nothing at all"). An element is *active* only relative to those.

THE STATISTICS, and why they are resampling rather than a parametric test.

Integrations are not independent replicates of a clean measurement: a CRE is
present at wildly different copy number in different cells, and the number of
integrations differs per element. So:

  * **Activity** is bootstrapped. For each element, resample its integrations
    WITH replacement, n = that element's own integration count, and resample
    minP and noP to the SAME n. Recompute the per-cluster means each time.
    Matching n is the point -- an element seen in 40 integrations and a control
    pooled over 4000 would otherwise be compared at different precision, and
    the control would win on sample size alone.

  * **Specificity** is permuted. Shuffle the cell-type cluster labels across
    cells and recompute the best cluster's fold-change over all other
    clusters. That null says "how cell-type-specific does this element look
    when cell type carries no information", which is exactly the thing a
    fold-change cannot tell you by itself.

  * **p-values are empirical**, (count + 1) / (n + 1) -- never zero, because a
    resampling test cannot resolve below its own resolution, and BH-corrected
    across elements within each biological replicate.

INPUT for the statistics is one row per integration per cell, the shape of the
joined tables deposited under GSE217686 / GSE217689:

    cell_bc  oBC  mBC  CRE_id  CRE_class  UMIs_mBC  gex_UMI  cluster_id  biol_rep

Subcommands:

    extract-bc    FASTQ -> barcodes, fixed position + constant-region check
    count-bc      count reads per barcode (or per barcode pair), threshold
    subassembly   oBC <-> mBC dictionary, keeping only unique pairings
    activity      bootstrap activity per CRE against minP / noP
    specificity   permute cluster labels -> null for cell-type specificity
    call          empirical p-values, BH correction, active/specific calls
    pipeline      activity -> specificity -> call
    selftest      synthetic data end to end, no inputs needed
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "scQers"

log = logging.getLogger("scqers")

# The library's two floors. An element is judged against these, never in the
# absolute -- see the module docstring.
CONTROL_IDS = ("minP", "noP")

DEFAULT_BOOTSTRAPS = 200
DEFAULT_PERMUTATIONS = 200


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")


def _open(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# --------------------------------------------------------------------------
# Barcode extraction
# --------------------------------------------------------------------------


def _fastq_records(path):
    """Yield (name, seq) from a FASTQ, gzipped or not."""
    with _open(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield name[1:].split()[0], seq


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def extract_barcode(seq: str, start: int, end: int, check_seq: str = "",
                    check_offset: int = 0) -> "str | None":
    """One barcode from one read, or None if the read fails its check.

    The constant-region check is what makes fixed-position extraction safe. A
    read with an indel upstream still yields *something* at positions
    [start:end] -- it is just the wrong something, and nothing downstream can
    tell that apart from a real barcode. Requiring the known constant sequence
    immediately after the barcode throws those away instead of silently
    inventing barcodes that were never in the library.
    """
    if len(seq) < end:
        return None
    bc = seq[start:end]
    if "N" in bc:
        return None
    if check_seq:
        at = end + check_offset
        if seq[at:at + len(check_seq)] != check_seq:
            return None
    return bc


def cmd_extract_bc(args) -> int:
    out = Path(args.out_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    kept = total = 0

    if args.in_r2:
        # Paired: one barcode from each mate, emitted as a pair. This is the
        # oBC-mBC subassembly amplicon.
        r2 = _fastq_records(args.in_r2)
        with _open(out, "wt") as fh:
            fh.write(f"{args.r1_name}\t{args.r2_name}\n")
            for (_, s1), (_, s2) in zip(_fastq_records(args.in_r1), r2):
                total += 1
                b1 = extract_barcode(s1, args.start, args.end,
                                     args.check_seq, args.check_offset)
                b2 = extract_barcode(s2, args.start2, args.end2,
                                     args.check_seq2, args.check_offset)
                if b1 is None or b2 is None:
                    continue
                if args.rc_r1:
                    b1 = revcomp(b1)
                if args.rc_r2:
                    b2 = revcomp(b2)
                fh.write(f"{b1}\t{b2}\n")
                kept += 1
    else:
        with _open(out, "wt") as fh:
            fh.write(f"{args.r1_name}\n")
            for _, seq in _fastq_records(args.in_r1):
                total += 1
                bc = extract_barcode(seq, args.start, args.end,
                                     args.check_seq, args.check_offset)
                if bc is None:
                    continue
                if args.rc_r1:
                    bc = revcomp(bc)
                fh.write(f"{bc}\n")
                kept += 1

    frac = kept / total if total else 0.0
    print(f"reads read    : {total:,}")
    print(f"barcodes kept : {kept:,}  ({frac:.1%})")
    if args.check_seq and frac < 0.5 and total:
        print("\nWARNING: over half the reads failed the constant-region "
              f"check for {args.check_seq!r}. Check --start/--end and the "
              "expected constant sequence before trusting this.")
    print(f"written       : {out}")
    ls.record_analysis("scqers", subcommand="extract-bc", label=out.name,
                       inputs=[args.in_r1] + ([args.in_r2] if args.in_r2 else []),
                       outputs=[str(out)])
    return 0


# --------------------------------------------------------------------------
# Counting
# --------------------------------------------------------------------------


def cmd_count_bc(args) -> int:
    counts: Counter = Counter()
    header = None
    with _open(args.in_file) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            key = tuple(line.rstrip("\n").split("\t"))
            if any(k == "" for k in key):
                continue
            counts[key] += 1

    out = Path(args.out_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    ordered = counts.most_common()
    above = [(k, n) for k, n in ordered if n >= args.threshold]

    with _open(out, "wt") as fh:
        fh.write("\t".join(header) + "\tcount\n")
        for key, n in ordered if args.all else above:
            fh.write("\t".join(key) + f"\t{n}\n")

    total_reads = sum(counts.values())
    kept_reads = sum(n for _, n in above)
    print(f"distinct barcodes      : {len(counts):,}")
    print(f"above threshold ({args.threshold}) : {len(above):,}")
    print(f"reads in those         : {kept_reads:,} of {total_reads:,} "
          f"({kept_reads / total_reads:.1%})" if total_reads else "")
    if ordered:
        print(f"top count              : {ordered[0][1]:,}")
        med = ordered[len(ordered) // 2][1]
        print(f"median count           : {med:,}")
    print(f"written                : {out}")

    if args.plot:
        _plot_count_distribution([n for _, n in ordered], args.threshold,
                                 Path(args.plot))
        print(f"plot                   : {args.plot}")
    ls.record_analysis("scqers", subcommand="count-bc", label=out.name,
                       inputs=[args.in_file], outputs=[str(out)])
    return 0


def _plot_count_distribution(counts, threshold, path: Path) -> None:
    """Rank-vs-count on log axes. A real library shows a plateau of true
    barcodes and then a cliff into sequencing-error singletons; the threshold
    belongs at the cliff, and the only honest way to choose it is to look."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not available; skipping plot")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 3.6))
    ax.plot(range(1, len(counts) + 1), counts, lw=1.2, color="#38707f")
    ax.axhline(threshold, color="#c0392b", ls="--", lw=1,
               label=f"threshold = {threshold}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("barcode rank")
    ax.set_ylabel("reads")
    ax.set_title("Barcode count distribution")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------
# Subassembly: oBC <-> mBC
# --------------------------------------------------------------------------


def cmd_subassembly(args) -> int:
    """Build the barcode dictionary, keeping only unambiguous pairings.

    A barcode that maps to two different partners is a chimera, a recombinant,
    or a collision. Keeping it would attribute one element's reporter counts
    to another element -- the single most damaging error available in this
    assay, because it is invisible downstream. The rule is therefore a
    dominance ratio, not a majority: the best partner must beat the runner-up
    by --min-ratio, or the barcode is dropped entirely.
    """
    pairs: Counter = Counter()
    with _open(args.in_file) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i1, i2 = header.index(args.col1), header.index(args.col2)
        except ValueError:
            print(f"ERROR: columns {args.col1!r}/{args.col2!r} not in "
                  f"{header}", file=sys.stderr)
            return 2
        icount = header.index("count") if "count" in header else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i1, i2):
                continue
            n = int(f[icount]) if icount is not None else 1
            pairs[(f[i1], f[i2])] += n

    by_first: "dict[str, Counter]" = defaultdict(Counter)
    for (a, b), n in pairs.items():
        by_first[a][b] += n

    resolved, ambiguous, low = [], 0, 0
    for a, partners in by_first.items():
        ranked = partners.most_common()
        best, best_n = ranked[0]
        if best_n < args.min_count:
            low += 1
            continue
        runner = ranked[1][1] if len(ranked) > 1 else 0
        if runner and best_n < args.min_ratio * runner:
            ambiguous += 1
            continue
        resolved.append((a, best, best_n, runner))

    out = Path(args.out_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    with _open(out, "wt") as fh:
        fh.write(f"{args.col1}\t{args.col2}\tcount\trunner_up_count\n")
        for a, b, n, r in sorted(resolved, key=lambda x: -x[2]):
            fh.write(f"{a}\t{b}\t{n}\t{r}\n")

    print(f"{args.col1} seen            : {len(by_first):,}")
    print(f"uniquely paired          : {len(resolved):,}")
    print(f"dropped, ambiguous       : {ambiguous:,}  "
          f"(best < {args.min_ratio}x runner-up)")
    print(f"dropped, low count       : {low:,}  (< {args.min_count})")
    print(f"written                  : {out}")
    ls.record_analysis("scqers", subcommand="subassembly", label=out.name,
                       inputs=[args.in_file], outputs=[str(out)])
    return 0


# --------------------------------------------------------------------------
# The count table the statistics run on
# --------------------------------------------------------------------------

REQUIRED_COLS = ("cell_bc", "CRE_id", "UMIs_mBC", "cluster_id")


class CountTable:
    """One row per integration per cell, plus the derived normalisation."""

    def __init__(self, rows: "list[dict]"):
        self.rows = rows
        gex = [r["gex_UMI"] for r in rows if r["gex_UMI"] > 0]
        self.mean_gex = sum(gex) / len(gex) if gex else 0.0
        for r in rows:
            # Depth-normalised reporter count. With no gex_UMI column this
            # degrades to the raw count rather than silently dividing by zero.
            r["norm_mBC"] = (r["UMIs_mBC"] / r["gex_UMI"] * self.mean_gex
                             if r["gex_UMI"] > 0 and self.mean_gex else
                             float(r["UMIs_mBC"]))

    @property
    def clusters(self) -> "list[str]":
        return sorted({r["cluster_id"] for r in self.rows})

    @property
    def replicates(self) -> "list[str]":
        return sorted({r["biol_rep"] for r in self.rows})

    def cres(self) -> "list[str]":
        return sorted({r["CRE_id"] for r in self.rows
                       if r["CRE_id"] not in CONTROL_IDS})


def load_counts(path, cluster_groups: "dict | None" = None) -> CountTable:
    rows = []
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in REQUIRED_COLS if c not in header]
        if missing:
            raise SystemExit(
                f"ERROR: {path} is missing required column(s): "
                f"{', '.join(missing)}\nfound: {', '.join(header)}")
        idx = {c: header.index(c) for c in header}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue

            def get(col, default=""):
                return f[idx[col]] if col in idx else default

            try:
                umis = float(get("UMIs_mBC", "0"))
                gex = float(get("gex_UMI", "0") or 0)
            except ValueError:
                continue
            cluster = get("cluster_id")
            if cluster_groups:
                cluster = cluster_groups.get(cluster, "")
                if not cluster:
                    continue   # cluster not assigned to any group: excluded
            rows.append({
                "cell_bc":   get("cell_bc"),
                "oBC":       get("oBC"),
                "mBC":       get("mBC"),
                "CRE_id":    get("CRE_id"),
                "CRE_class": get("CRE_class", "devCRE"),
                "UMIs_mBC":  umis,
                "gex_UMI":   gex,
                "cluster_id": cluster,
                "biol_rep":  get("biol_rep", "A"),
            })
    if not rows:
        raise SystemExit(f"ERROR: no usable rows in {path}")
    return CountTable(rows)


def _summarise(rows, clusters) -> "dict[str, float]":
    """Per-cluster mean normalised activity, and the mean over every OTHER
    cluster. The ratio of those two is the specificity statistic."""
    by_cluster: "dict[str, list[float]]" = defaultdict(list)
    for r in rows:
        by_cluster[r["cluster_id"]].append(r["norm_mBC"])
    total_sum = sum(r["norm_mBC"] for r in rows)
    total_n = len(rows)
    out = {}
    for c in clusters:
        vals = by_cluster.get(c)
        if not vals:
            continue
        mean_in = sum(vals) / len(vals)
        n_out = total_n - len(vals)
        mean_out = (total_sum - sum(vals)) / n_out if n_out else float("nan")
        out[c] = (mean_in, mean_out)
    return out


def _max_cluster(summary) -> "tuple[str, float, float]":
    """The best-expressing cluster, its activity, and its fold-change."""
    best_c, best_mean, best_fc = "", float("-inf"), float("nan")
    for c, (mean_in, mean_out) in summary.items():
        if mean_in > best_mean:
            best_c, best_mean = c, mean_in
            best_fc = (mean_in / mean_out
                       if mean_out and mean_out > 0 else float("nan"))
    return best_c, best_mean, best_fc


# --------------------------------------------------------------------------
# Activity: bootstrap
# --------------------------------------------------------------------------


def bootstrap_activity(tab: CountTable, cre: str, rep: str, n_boot: int,
                       rng: random.Random):
    """Resample integrations for one CRE and its controls at matched n."""
    rows_rep = [r for r in tab.rows if r["biol_rep"] == rep]
    pools = {cre: [r for r in rows_rep if r["CRE_id"] == cre]}
    for ctrl in CONTROL_IDS:
        pools[ctrl] = [r for r in rows_rep if r["CRE_id"] == ctrl]
    if not pools[cre]:
        return None
    n = len(pools[cre])
    clusters = tab.clusters

    obs = {}
    for name, pool in pools.items():
        if pool:
            obs[name] = _max_cluster(_summarise(pool, clusters))

    draws = {name: [] for name in pools}
    for _ in range(n_boot):
        for name, pool in pools.items():
            if not pool:
                continue
            sample = [pool[rng.randrange(len(pool))] for _ in range(n)]
            draws[name].append(_max_cluster(_summarise(sample, clusters)))
    return {"n_integrations": n, "observed": obs, "draws": draws,
            "n_cells": len({r["cell_bc"] for r in pools[cre]})}


def cmd_activity(args) -> int:
    tab = load_counts(args.counts, _cluster_groups(args))
    rng = random.Random(args.seed)
    out_dir = Path(args.out_dir) if args.out_dir else (OUT_DIR / _stamp("activity"))
    out_dir.mkdir(parents=True, exist_ok=True)

    cres = [args.cre] if args.cre else tab.cres()
    reps = tab.replicates
    log.info("%d CRE(s) x %d replicate(s), %d bootstraps",
             len(cres), len(reps), args.bootstraps)

    results = []
    for cre in cres:
        for rep in reps:
            b = bootstrap_activity(tab, cre, rep, args.bootstraps, rng)
            if b is None:
                continue
            own = [d[1] for d in b["draws"][cre]]
            ctrl = [d[1] for name in CONTROL_IDS
                    for d in b["draws"].get(name, [])]
            best_ids = [d[0] for d in b["draws"][cre]]
            modal = Counter(best_ids).most_common(1)[0] if best_ids else ("", 0)
            results.append({
                "CRE_id": cre, "biol_rep": rep,
                "n_integrations": b["n_integrations"],
                "n_cells": b["n_cells"],
                "observed_max_cluster": b["observed"][cre][0],
                "observed_activity": round(b["observed"][cre][1], 4),
                "observed_fold_change": _r(b["observed"][cre][2]),
                "boot_q25": _r(_quantile(own, 0.25)),
                "boot_q50": _r(_quantile(own, 0.50)),
                "boot_q75": _r(_quantile(own, 0.75)),
                "control_q50": _r(_quantile(ctrl, 0.50)) if ctrl else "",
                "modal_max_cluster": modal[0],
                "modal_max_cluster_frac": _r(modal[1] / len(best_ids))
                                          if best_ids else "",
                "p_activity": _r(_empirical_p_activity(own, ctrl)),
            })

    path = out_dir / "activity.tsv"
    _write_tsv(path, results)
    print(f"{len(results)} CRE x replicate result(s)")
    print(f"written: {path}")
    _print_head(results, ("CRE_id", "biol_rep", "n_integrations",
                          "observed_max_cluster", "boot_q50", "control_q50",
                          "p_activity"))
    ls.record_analysis("scqers", subcommand="activity", label=out_dir.name,
                       inputs=[args.counts], outputs=[str(path)])
    return 0


def _empirical_p_activity(own, ctrl) -> float:
    """Mean over bootstrap draws of P(control >= this draw).

    Averaged across the element's own bootstrap distribution rather than
    computed once at the point estimate: that propagates the element's own
    sampling uncertainty into the p-value instead of pretending the observed
    activity was measured exactly.
    """
    if not own or not ctrl:
        return float("nan")
    ctrl_sorted = sorted(ctrl)
    n = len(ctrl_sorted)
    total = 0.0
    for v in own:
        ge = n - _bisect_left(ctrl_sorted, v)
        total += (ge + 1) / (n + 1)
    return total / len(own)


# --------------------------------------------------------------------------
# Specificity: permutation
# --------------------------------------------------------------------------


def cmd_specificity(args) -> int:
    tab = load_counts(args.counts, _cluster_groups(args))
    rng = random.Random(args.seed)
    out_dir = Path(args.out_dir) if args.out_dir else (OUT_DIR / _stamp("specificity"))
    out_dir.mkdir(parents=True, exist_ok=True)

    cres = [args.cre] if args.cre else tab.cres()
    clusters = tab.clusters
    results = []

    for cre in cres:
        for rep in tab.replicates:
            pool = [r for r in tab.rows
                    if r["CRE_id"] == cre and r["biol_rep"] == rep]
            if len(pool) < args.min_integrations:
                continue
            _, obs_act, obs_fc = _max_cluster(_summarise(pool, clusters))
            if not math.isfinite(obs_fc):
                continue

            # Permute the CELL -> cluster mapping, not the rows: two
            # integrations in the same cell must keep the same cluster or the
            # null quietly gains information the real data never had.
            cell_to_cluster = {r["cell_bc"]: r["cluster_id"] for r in pool}
            cells = list(cell_to_cluster)
            labels = [cell_to_cluster[c] for c in cells]

            null = []
            for _ in range(args.permutations):
                rng.shuffle(labels)
                shuffled = dict(zip(cells, labels))
                permuted = [{**r, "cluster_id": shuffled[r["cell_bc"]]}
                            for r in pool]
                _, _, fc = _max_cluster(_summarise(permuted, clusters))
                if math.isfinite(fc):
                    null.append(fc)

            p = ((sum(1 for v in null if v >= obs_fc) + 1) / (len(null) + 1)
                 if null else float("nan"))
            results.append({
                "CRE_id": cre, "biol_rep": rep,
                "n_integrations": len(pool),
                "observed_fold_change": _r(obs_fc),
                "observed_activity": _r(obs_act),
                "null_q50": _r(_quantile(null, 0.50)) if null else "",
                "null_q95": _r(_quantile(null, 0.95)) if null else "",
                "p_specificity": _r(p),
            })

    path = out_dir / "specificity.tsv"
    _write_tsv(path, results)
    print(f"{len(results)} CRE x replicate result(s)")
    print(f"written: {path}")
    _print_head(results, ("CRE_id", "biol_rep", "observed_fold_change",
                          "null_q95", "p_specificity"))
    ls.record_analysis("scqers", subcommand="specificity", label=out_dir.name,
                       inputs=[args.counts], outputs=[str(path)])
    return 0


# --------------------------------------------------------------------------
# Calling
# --------------------------------------------------------------------------


def benjamini_hochberg(pvals) -> "list[float]":
    """BH step-up. Applied WITHIN a replicate, never pooled across them:
    pooling would treat the same element measured twice as two independent
    tests and inflate the discovery count."""
    n = len(pvals)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    out = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        prev = min(prev, pvals[i] * n / k)
        out[i] = min(1.0, prev)
    return out


def cmd_call(args) -> int:
    act = _read_tsv(args.activity)
    spec = _read_tsv(args.specificity) if args.specificity else []
    spec_by = {(r["CRE_id"], r["biol_rep"]): r for r in spec}

    by_rep: "dict[str, list[dict]]" = defaultdict(list)
    for r in act:
        by_rep[r["biol_rep"]].append(r)

    merged = []
    for rep, rows in by_rep.items():
        pa = [_f(r.get("p_activity")) for r in rows]
        qa = benjamini_hochberg(pa)
        ps = [_f(spec_by.get((r["CRE_id"], rep), {}).get("p_specificity"))
              for r in rows]
        qs = benjamini_hochberg(ps) if any(math.isfinite(v) for v in ps) else \
            [float("nan")] * len(rows)
        for r, a, s in zip(rows, qa, qs):
            sp = spec_by.get((r["CRE_id"], rep), {})
            merged.append({
                **r,
                "q_activity": _r(a),
                "p_specificity": sp.get("p_specificity", ""),
                "q_specificity": _r(s),
                "observed_fold_change": sp.get("observed_fold_change",
                                               r.get("observed_fold_change", "")),
                "active": int(a <= args.fdr),
                "specific": int(math.isfinite(s) and s <= args.fdr
                                and _f(sp.get("observed_fold_change")) >= args.min_fc),
            })

    # Reproducibility: an element counts only if it clears the bar in every
    # replicate it was measured in. One replicate is an observation; agreement
    # across replicates is the claim.
    reps_per_cre = defaultdict(list)
    for m in merged:
        reps_per_cre[m["CRE_id"]].append(m)
    for cre, rows in reps_per_cre.items():
        repro_a = all(r["active"] for r in rows) and len(rows) >= args.min_reps
        repro_s = all(r["specific"] for r in rows) and len(rows) >= args.min_reps
        for r in rows:
            r["reproducibly_active"] = int(repro_a)
            r["reproducibly_specific"] = int(repro_s)

    out = Path(args.out_file) if args.out_file else (
        Path(args.activity).parent / "calls.tsv")
    _write_tsv(out, merged)

    n_cre = len(reps_per_cre)
    n_act = sum(1 for rows in reps_per_cre.values() if rows[0]["reproducibly_active"])
    n_spec = sum(1 for rows in reps_per_cre.values() if rows[0]["reproducibly_specific"])
    print(f"elements tested            : {n_cre}")
    print(f"reproducibly active        : {n_act}   (FDR <= {args.fdr}, "
          f"all of >= {args.min_reps} replicates)")
    print(f"reproducibly cell-type-specific : {n_spec}   "
          f"(FDR <= {args.fdr} and fold-change >= {args.min_fc})")
    print(f"written                    : {out}")
    ls.record_analysis("scqers", subcommand="call", label=out.name,
                       inputs=[args.activity], outputs=[str(out)])
    return 0


# --------------------------------------------------------------------------
# Pipeline / selftest
# --------------------------------------------------------------------------


def cmd_pipeline(args) -> int:
    out_dir = Path(args.out_dir) if args.out_dir else (OUT_DIR / _stamp("pipeline"))
    out_dir.mkdir(parents=True, exist_ok=True)
    ns = argparse.Namespace(**vars(args))
    ns.out_dir = str(out_dir)
    ns.cre = None
    if cmd_activity(ns):
        return 1
    if cmd_specificity(ns):
        return 1
    ns.activity = str(out_dir / "activity.tsv")
    ns.specificity = str(out_dir / "specificity.tsv")
    ns.out_file = str(out_dir / "calls.tsv")
    return cmd_call(ns)


def _synth_counts(path: Path, rng: random.Random) -> None:
    """A library with a known answer: one strong ubiquitous element, one
    cluster-specific element, one dead element, plus both controls."""
    clusters = ["pluri", "ecto", "endo", "meso"]
    spec = {
        "CRE_strong":   {c: 40.0 for c in clusters},
        "CRE_specific": {"pluri": 8.0, "ecto": 8.0, "endo": 60.0, "meso": 8.0},
        "CRE_dead":     {c: 3.0 for c in clusters},
        "minP":         {c: 4.0 for c in clusters},
        "noP":          {c: 2.0 for c in clusters},
    }
    with open(path, "wt") as fh:
        fh.write("cell_bc\toBC\tmBC\tCRE_id\tCRE_class\tUMIs_mBC\tgex_UMI\t"
                 "cluster_id\tbiol_rep\n")
        cell = 0
        for rep in ("A", "B"):
            for cre, per_cluster in spec.items():
                n = 60 if cre in CONTROL_IDS else 40
                for i in range(n):
                    cell += 1
                    cl = clusters[rng.randrange(len(clusters))]
                    gex = rng.randint(2000, 8000)
                    lam = per_cluster[cl] * gex / 5000.0
                    umis = max(0, int(rng.gauss(lam, max(1.0, lam * 0.4))))
                    cls = "control" if cre in CONTROL_IDS else "devCRE"
                    fh.write(f"cell{cell}\toBC{cell}\tmBC{cell}\t{cre}\t{cls}\t"
                             f"{umis}\t{gex}\t{cl}\t{rep}\n")


def cmd_selftest(args) -> int:
    out_dir = Path(args.out_dir) if args.out_dir else (OUT_DIR / _stamp("selftest"))
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(0)
    counts = out_dir / "synthetic_counts.tsv"
    _synth_counts(counts, rng)
    print(f"synthetic count table: {counts}\n")

    ns = argparse.Namespace(
        counts=str(counts), out_dir=str(out_dir), cre=None, seed=0,
        bootstraps=100, permutations=100, min_integrations=5,
        cluster_groups=None, fdr=0.05, min_fc=1.5, min_reps=2,
        activity=str(out_dir / "activity.tsv"),
        specificity=str(out_dir / "specificity.tsv"),
        out_file=str(out_dir / "calls.tsv"))
    cmd_activity(ns)
    print()
    cmd_specificity(ns)
    print()
    cmd_call(ns)

    calls = {r["CRE_id"]: r for r in _read_tsv(ns.out_file)}
    print("\n=== expected vs observed ===")
    ok = True
    checks = [
        ("CRE_strong",   "reproducibly_active",   "1"),
        ("CRE_dead",     "reproducibly_active",   "0"),
        ("CRE_specific", "reproducibly_active",   "1"),
        ("CRE_specific", "reproducibly_specific", "1"),
        ("CRE_strong",   "reproducibly_specific", "0"),
    ]
    for cre, field, want in checks:
        got = calls.get(cre, {}).get(field, "?")
        good = str(got) == want
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {cre:14} {field:24} "
              f"expected {want}, got {got}")
    print("\nall checks passed" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------


def _stamp(tag: str) -> str:
    return time.strftime("%Y%m%d_%H%M%S_") + tag


def _r(v, nd=4):
    try:
        return "" if v is None or not math.isfinite(float(v)) else round(float(v), nd)
    except (TypeError, ValueError):
        return ""


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _bisect_left(a, x) -> int:
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _quantile(vals, q: float) -> float:
    if not vals:
        return float("nan")
    s = sorted(vals)
    if len(s) == 1:
        return s[0]
    pos = q * (len(s) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def _cluster_groups(args) -> "dict | None":
    """Optional cluster -> group mapping, as JSON {"group": [ids...]}.

    Analysing at the level of coarse lineages rather than raw Leiden clusters
    is usually what the biology is about, and doing it here means the grouping
    is recorded alongside the result instead of being applied by hand upstream.
    """
    path = getattr(args, "cluster_groups", None)
    if not path:
        return None
    spec = json.loads(Path(path).read_text())
    return {str(c): group for group, members in spec.items() for c in members}


def _write_tsv(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0])
    with open(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _read_tsv(path) -> "list[dict]":
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def _print_head(rows, cols, n=8) -> None:
    if not rows:
        return
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows[:n]))
              for c in cols}
    print()
    print("  " + "  ".join(c.ljust(widths[c]) for c in cols))
    for r in rows[:n]:
        print("  " + "  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))
    if len(rows) > n:
        print(f"  … {len(rows) - n} more")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent scqers",
        description="scQers: quantitative single-cell enhancer reporters — "
                    "barcode extraction through to reproducibly active and "
                    "cell-type-specific CRE calls.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("extract-bc", help="FASTQ -> barcodes at fixed "
                                          "positions, with a constant-region check.")
    s.add_argument("--in-r1", required=True)
    s.add_argument("--in-r2", help="Second mate, for paired barcode amplicons.")
    s.add_argument("--out-file", required=True)
    s.add_argument("--start", type=int, default=0)
    s.add_argument("--end", type=int, default=15)
    s.add_argument("--start2", type=int, default=0)
    s.add_argument("--end2", type=int, default=15)
    s.add_argument("--check-seq", default="",
                   help="Constant sequence expected just after the barcode "
                        "(e.g. GCT). Reads that do not show it are discarded.")
    s.add_argument("--check-seq2", default="")
    s.add_argument("--check-offset", type=int, default=0)
    s.add_argument("--r1-name", default="BC")
    s.add_argument("--r2-name", default="BC2")
    s.add_argument("--rc-r1", action="store_true",
                   help="Reverse-complement the R1 barcode.")
    s.add_argument("--rc-r2", action="store_true")
    s.set_defaults(func=cmd_extract_bc)

    s = sub.add_parser("count-bc", help="Count reads per barcode or barcode pair.")
    s.add_argument("--in-file", required=True)
    s.add_argument("--out-file", required=True)
    s.add_argument("--threshold", type=int, default=300)
    s.add_argument("--all", action="store_true",
                   help="Write every barcode, not just those above threshold.")
    s.add_argument("--plot", help="Write a rank-vs-count distribution here.")
    s.set_defaults(func=cmd_count_bc)

    s = sub.add_parser("subassembly",
                       help="Resolve a barcode dictionary, dropping ambiguous pairs.")
    s.add_argument("--in-file", required=True)
    s.add_argument("--out-file", required=True)
    s.add_argument("--col1", default="oBC")
    s.add_argument("--col2", default="mBC")
    s.add_argument("--min-count", type=int, default=3)
    s.add_argument("--min-ratio", type=float, default=5.0,
                   help="Best partner must beat the runner-up by this factor.")
    s.set_defaults(func=cmd_subassembly)

    def stats_args(sp):
        sp.add_argument("--counts", required=True,
                        help="Joined count table: one row per integration per cell.")
        sp.add_argument("--out-dir")
        sp.add_argument("--cre", help="Restrict to one CRE.")
        sp.add_argument("--seed", type=int, default=0)
        sp.add_argument("--cluster-groups",
                        help='JSON {"lineage": [cluster ids…]} to analyse at '
                             "coarse lineage level.")

    s = sub.add_parser("activity",
                       help="Bootstrap activity per CRE against minP / noP.")
    stats_args(s)
    s.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    s.set_defaults(func=cmd_activity)

    s = sub.add_parser("specificity",
                       help="Permute cluster labels -> cell-type specificity null.")
    stats_args(s)
    s.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    s.add_argument("--min-integrations", type=int, default=5)
    s.set_defaults(func=cmd_specificity)

    s = sub.add_parser("call", help="Empirical p-values -> BH -> final calls.")
    s.add_argument("--activity", required=True)
    s.add_argument("--specificity")
    s.add_argument("--out-file")
    s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--min-fc", type=float, default=1.5)
    s.add_argument("--min-reps", type=int, default=1)
    s.set_defaults(func=cmd_call)

    s = sub.add_parser("pipeline", help="activity -> specificity -> call.")
    stats_args(s)
    s.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    s.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    s.add_argument("--min-integrations", type=int, default=5)
    s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--min-fc", type=float, default=1.5)
    s.add_argument("--min-reps", type=int, default=1)
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("selftest",
                       help="Synthetic data end to end; needs no inputs.")
    s.add_argument("--out-dir")
    s.set_defaults(func=cmd_selftest)

    return p


def main(argv=None) -> int:
    setup_logging()
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
