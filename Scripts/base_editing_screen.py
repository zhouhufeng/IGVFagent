#!/usr/bin/env python3
"""Base-editing screens, following crispr-bean's method where the data allows.

    igvfagent bean discover IGVFDS6464SOVZ
    igvfagent bean count    IGVFDS6464SOVZ
    igvfagent bean analyze  IGVFDS6464SOVZ --label ldl_abe

WHY A SEPARATE PATH FROM crispr-screen. A base editor edits the guide's own
locus as well as its target, so the protospacer sequenced back from the cell
carries substitutions of its own. Exact matching throws those reads away.
Measured on IGVFDS6464SOVZ (ABE, 8,192 guides):

    exact match only     36.7% of reads assigned
    A>G aware            62.5%          (+25.8 points)
    C>T aware            36.8%          (+0.1, so it is not a CBE)

A collaborator from the lab that produced the data flagged the low rate
before we measured it, and named the cause: "it is important to allow for
self-editing in gRNA assignment (allow for A2G edits in the gRNA when
aligning)".

SOURCE. crispr-bean ("Base Editing screens' Activity-Normalized variant
effect size estimation"), Pinello Lab:
    repository  https://github.com/pinellolab/crispr-bean   (AGPL-3.0)
    docs        https://pinellolab.github.io/crispr-bean/
    method read from bean/mapping/GuideEditCounter.py
No BEAN source is copied here. IGVFagent is Apache-2.0 and BEAN is AGPL-3.0,
so the METHOD was read and reimplemented; to use BEAN's own model, invoke
`bean` as a separate program rather than vendoring it.

WHAT IS TAKEN FROM BEAN. The matching method: crispr-bean's GuideEditCounter
compares mask_sequence(read) against the masked library, normalising the
edited base to its product on both sides rather than allowing free
mismatches. Free mismatches would also absorb sequencing error and cross-map
similar guides; masking forgives only the substitution the editor makes. The
self-edit count per read is then used as a per-guide editing-activity
estimate, which is the quantity BEAN's activity normalisation rests on.

WHAT IS NOT, AND CANNOT BE, TAKEN FROM BEAN. Two of BEAN's inputs are not in
what IGVF publishes for this screen:

  guide barcode -- BEAN reads a barcode from R2 and splits counts into
    bcmatch / semimatch / nomatch, which is how it resolves guides that
    collapse onto one masked sequence. The IGVF library table has no barcode
    column, so 13.8% of reads here are compatible with more than one guide
    and are counted as ambiguous rather than assigned to a guess.

  reporter allele -- the `reporter` column is present but empty for all
    8,192 guides, so reporter-allele counting and bystander/tiling analysis
    are impossible from this table. Activity is therefore estimated from
    self-editing of the protospacer, not from a reporter.

And BEAN's `run` step fits a Bayesian variant/tiling model with accessibility
covariates. That is NOT reimplemented here: this scores enrichment between
sorted bins with an activity adjustment, and says so. For the full model,
run BEAN itself on the FASTQs -- this tool prints the command.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import raw_data_pipeline as rp                                  # noqa: E402
import crispr_screen_analysis as cs                             # noqa: E402
from _stats import benjamini_hochberg, moderated_t, _percentile  # noqa: E402

ROOT = rp.ROOT
OUT_DIR = ROOT / "Docs" / "BaseEditingScreen"
LOG_DIR = ROOT / "Docs" / "Logs"

# A guide whose protospacer is never seen edited has not been shown to be
# inactive -- it may simply be shallow. Below this many assigned reads the
# activity estimate is not reported as a rate at all.
MIN_READS_FOR_ACTIVITY = 50


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    p = LOG_DIR / f"bean_{time.strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)-6s %(message)s",
        handlers=[logging.FileHandler(p), logging.StreamHandler(sys.stdout)])
    return p


# ─── Library ────────────────────────────────────────────────────────────────

def load_library(accession: str) -> dict:
    """Spacers, control classes and the base editor for a screen's library."""
    lib = rp.find_guide_library(accession)
    if not lib["resolved"]:
        return {"resolved": False, "why": lib["why"]}
    gf = lib["guide_files"][0]
    idx = rp.load_guide_index(gf["id"] or gf["accession"])
    if not idx["resolved"]:
        return {"resolved": False, "why": idx["why"]}
    spacer = next((c for c in idx["candidates"]
                   if c["column"] in ("spacer", "protospacer")), None)
    if not spacer:
        return {"resolved": False,
                "why": f"{gf['accession']} has no spacer column; a base-editing "
                       f"screen is counted on the protospacer"}
    editor = rp.detect_base_editor(spacer["seq_to_guide"].values())
    return {"resolved": True, "file": gf["accession"], "index": idx,
            "spacer": spacer, "editor": editor,
            "n_guides": idx["n_rows"],
            "type_of": idx["type"], "target_of": idx["target"]}


def calibrate_editor(spacers: "dict[str, str]", accession: str,
                      sample: int = 40000) -> dict:
    """Measure which editor's masking actually recovers reads.

    The library names usually say (every guide here is named ...__ABE_...),
    but naming is a convention and the reads are evidence. Both are reported
    so a disagreement is visible rather than resolved silently.
    """
    reads = []
    comp = str.maketrans("ACGTN", "TGCAN")
    for f in rp.biological_fastqs(accession):
        dest = rp.FASTQ_CACHE / Path(str(f.get("href") or f.get("accession"))).name
        if not dest.exists():
            rp.FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
            rp.portal_download(f["href"], dest)
        import gzip
        op = gzip.open if dest.read_bytes()[:2] == b"\x1f\x8b" else open
        with op(dest, "rt", errors="replace") as fh:      # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    reads.append(line.strip().upper())
                    if len(reads) >= sample:
                        break
        if reads:
            break
    if not reads:
        return {"reads": 0, "tested": [], "chosen": None}

    plain = rp.build_matcher(spacers)
    exact = sum(1 for r in reads
                if rp.match_read(r, plain)
                or rp.match_read(r.translate(comp)[::-1], plain))
    tested = []
    for name in rp.BASE_EDITS:
        mm = rp.build_masked_matcher(spacers, name)
        rec = amb = 0
        for r in reads:
            if (rp.match_read(r, plain)
                    or rp.match_read(r.translate(comp)[::-1], plain)):
                continue
            for s_ in (r, r.translate(comp)[::-1]):
                g, _n, a = rp.match_read_masked(s_, mm)
                if g and not a:
                    rec += 1
                    break
                if a:
                    amb += 1
                    break
        tested.append({"editor": name, "recovered": rec / len(reads),
                        "ambiguous": amb / len(reads),
                        "total": (exact + rec) / len(reads),
                        "masked_collisions": mm["masked_collisions"]})
    tested.sort(key=lambda t: -t["recovered"])
    best = tested[0]
    return {"reads": len(reads), "exact": exact / len(reads),
             "tested": tested,
             "chosen": best["editor"] if best["recovered"] > 0.01 else None}


# ─── Counting with self-edit tracking ───────────────────────────────────────

def count_library(accession: str, masked: dict, plain: dict,
                   max_reads: "Optional[int]" = None) -> dict:
    """Per-guide counts for one sorted bin, tracking self-editing.

    Returns counts, edited_counts (reads whose protospacer carried at least
    one edit), total self-edits per guide, and the read dispositions.
    Ambiguous reads are counted as ambiguous and assigned to nobody: without
    the guide barcode BEAN uses, choosing between the candidates would
    fabricate the number.
    """
    import gzip
    comp = str.maketrans("ACGTN", "TGCAN")
    counts: "Counter[str]" = Counter()
    edited: "Counter[str]" = Counter()
    edit_sum: "Counter[str]" = Counter()
    disp = Counter()
    for f in rp.biological_fastqs(accession):
        dest = rp.FASTQ_CACHE / Path(str(f.get("href") or f.get("accession"))).name
        if not dest.exists():
            rp.FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
            rp.portal_download(f["href"], dest)
        op = gzip.open if dest.read_bytes()[:2] == b"\x1f\x8b" else open
        with op(dest, "rt", errors="replace") as fh:      # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                if max_reads and disp["reads"] >= max_reads:
                    break
                disp["reads"] += 1
                seq = line.strip().upper()
                gid = (rp.match_read(seq, plain)
                       or rp.match_read(seq.translate(comp)[::-1], plain))
                if gid:
                    counts[gid] += 1
                    disp["exact"] += 1
                    continue
                hit = None
                ambiguous = False
                for s_ in (seq, seq.translate(comp)[::-1]):
                    g, n, a = rp.match_read_masked(s_, masked)
                    if g and not a:
                        hit = (g, n)
                        break
                    if a:
                        ambiguous = True
                        break
                if hit:
                    gid, n = hit
                    counts[gid] += 1
                    edited[gid] += 1
                    edit_sum[gid] += n
                    disp["self_edited"] += 1
                elif ambiguous:
                    disp["ambiguous"] += 1
                else:
                    disp["unassigned"] += 1
    return {"counts": counts, "edited": edited, "edit_sum": edit_sum,
             "disposition": dict(disp)}


def guide_activity(per_bin: "dict[str, dict]") -> "dict[str, dict]":
    """Per-guide editing activity, pooled across a screen's bins.

    activity = reads whose protospacer was edited / reads assigned to that
    guide. It is a proxy for how efficiently the editor acts at that guide's
    locus, and it is the quantity an activity normalisation needs.

    Pooled across bins on purpose: activity is a property of the guide and
    the editor, not of the sorted fraction a cell landed in, and per-bin
    estimates are far noisier.
    """
    tot: "Counter[str]" = Counter()
    ed: "Counter[str]" = Counter()
    es: "Counter[str]" = Counter()
    for d in per_bin.values():
        tot.update(d["counts"])
        ed.update(d["edited"])
        es.update(d["edit_sum"])
    out = {}
    for g, n in tot.items():
        enough = n >= MIN_READS_FOR_ACTIVITY
        out[g] = {
            "reads": n,
            "edited_reads": ed.get(g, 0),
            "activity": (ed.get(g, 0) / n) if enough else None,
            "mean_edits_per_edited_read": (es.get(g, 0) / ed[g]) if ed.get(g) else 0.0,
            "activity_estimable": enough,
        }
    return out


def activity_summary(act: "dict[str, dict]", type_of: "dict[str, str]") -> dict:
    """Distribution of activity, split by control class.

    Positive controls are guides whose target is known to move the phenotype.
    If activity were meaningless, controls and variants would look identical;
    that they do not is what makes the normalisation worth applying.
    """
    def vals(pred):
        return [a["activity"] for g, a in act.items()
                if a["activity"] is not None and pred(type_of.get(g, ""))]
    pos = vals(lambda t: "positive control" in (t or "").lower())
    var = vals(lambda t: "positive control" not in (t or "").lower())
    def stat(xs):
        if not xs:
            return None
        xs = sorted(xs)
        return {"n": len(xs), "median": round(_percentile(xs, 0.5), 4),
                 "p10": round(_percentile(xs, 0.10), 4),
                 "p90": round(_percentile(xs, 0.90), 4)}
    return {"positive_controls": stat(pos), "variants": stat(var),
             "n_estimable": sum(1 for a in act.values() if a["activity_estimable"]),
             "n_total": len(act)}


# Below this per-guide editing rate, log2_raw/activity is dominated by the
# divisor's own sampling error and is withheld rather than reported. 2% is an
# order of magnitude under the median editing activity at positive controls
# in the ABE screen measured here, so it excludes only guides with no usable
# rate -- not weakly-editing guides that still have one.
PER_EDIT_ACTIVITY_FLOOR = 0.02


# ─── Activity-normalised scoring ────────────────────────────────────────────

def score_screen(per_bin: "dict[str, dict]", bins: "list[dict]",
                  act: "dict[str, dict]", target_of: "dict[str, str]",
                  type_of: "dict[str, str]", low_pct: int,
                  min_count: int) -> "tuple[list[dict], dict]":
    """Per-guide enrichment between tails, reported raw AND activity-scaled.

    An editing screen measures the phenotype of cells in which the edit was
    MADE. A guide that edits 10% of the time dilutes its own effect roughly
    tenfold, so a raw log-ratio understates the edited cells' phenotype and a
    weakly-editing guide looks inactive when it is only underpowered. That
    conflation is what activity normalisation exists to prevent.

    Both numbers are reported, deliberately:

      log2_raw       what the sorted bins actually show for this guide
      log2_per_edit  log2_raw / activity, the implied effect in edited cells

    log2_per_edit divides by a rate estimated from finite counts, so it is
    noisy exactly where activity is low -- the place it changes the answer
    most. It is therefore reported with the activity it used, is omitted
    entirely when activity is not estimable, and never sets the p-value. The
    test stays on log2_raw, where the sampling model is honest, and low
    activity is surfaced as low power instead of being divided out.
    """
    reps = sorted({b["rep"] for b in bins})
    usable = [r for r in reps
              if any(b["rep"] == r and b["side"] == "bottom" and b["pct"] == low_pct
                     for b in bins)
              and any(b["rep"] == r and b["side"] == "top" and b["pct"] == low_pct
                      for b in bins)]
    ratios: "dict[str, list[float]]" = defaultdict(list)
    pvars: "dict[str, list[float]]" = defaultdict(list)
    ln2 = math.log(2.0)
    for rep in usable:
        lo = next(b for b in bins if b["rep"] == rep and b["side"] == "bottom"
                  and b["pct"] == low_pct)
        hi = next(b for b in bins if b["rep"] == rep and b["side"] == "top"
                  and b["pct"] == low_pct)
        clo = per_bin.get(lo["accession"], {}).get("counts")
        chi = per_bin.get(hi["accession"], {}).get("counts")
        if not clo or not chi:
            continue
        tlo, thi = max(sum(clo.values()), 1), max(sum(chi.values()), 1)
        for g in set(clo) | set(chi):
            a, b_ = clo.get(g, 0), chi.get(g, 0)
            if a + b_ < min_count:
                continue
            ratios[g].append(math.log2(((a + 0.5) / tlo) / ((b_ + 0.5) / thi)))
            pvars[g].append((1.0 / (a + 0.5) + 1.0 / (b_ + 0.5)) / (ln2 ** 2))

    sds = []
    for g, vals in ratios.items():
        if len(vals) >= 3:
            m = sum(vals) / len(vals)
            sds.append(math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1)))
    sd_floor = _percentile(sorted(x for x in sds if x > 0), 0.10) if sds else 0.0

    rows = []
    for g, vals in ratios.items():
        k = len(vals)
        sd_count = math.sqrt(sum(pvars[g]) / k) / math.sqrt(k)
        st = moderated_t(vals, sd_floor, sd_count)
        a = act.get(g, {})
        activity = a.get("activity")
        # log2_raw / activity is unbounded as activity -> 0, and on the
        # LDLR137-219 screen (activity ~0.001, because it is not a base
        # editing screen at all) it printed -436.54 as an effect size. A
        # number like that is not a weak estimate, it is no estimate: the
        # divisor is a rate measured from a handful of reads. Below the floor
        # the ratio is withheld and the low activity is what gets reported,
        # which is the same reasoning that keeps the p-value on log2_raw.
        per_edit = (st["mean"] / activity
                     if (activity or 0) >= PER_EDIT_ACTIVITY_FLOOR else None)
        rows.append({
            "guide": g, "target": target_of.get(g, g),
            "guide_type": type_of.get(g, ""),
            "n_reps": k, "reads": a.get("reads", 0),
            "activity": None if activity is None else round(activity, 4),
            "activity_estimable": a.get("activity_estimable", False),
            "log2_raw": round(st["mean"], 4),
            "log2_per_edit": None if per_edit is None else round(per_edit, 4),
            "sd": None if st["sd"] is None else round(st["sd"], 4),
            "sd_used": None if st["sd_used"] is None else round(st["sd_used"], 4),
            "t": round(st["t"], 4), "df": st["df"], "p_value": st["p_value"],
        })
    for r, q in zip(rows, benjamini_hochberg([r["p_value"] for r in rows])):
        r["fdr"] = q
    rows.sort(key=lambda r: (r["p_value"], -abs(r["log2_raw"])))
    lowact = [r for r in rows
              if r["activity_estimable"] and (r["activity"] or 0) < 0.05]
    return rows, {"replicates_used": usable, "sd_floor": round(sd_floor, 4),
                   "n_low_activity": len(lowact),
                   "test_basis": "log2_raw (activity is reported, not divided out)"}


# ─── Handing counts to the real BEAN ───────────────────────────────────────

def bean_available() -> "tuple[bool, str]":
    """Is the real `bean` runnable, and which build?

    Probed with `--help`, NOT `--version`: BEAN has no --version flag and
    answers "error: unrecognized arguments: --version" with exit 2. An
    agent-authored wrapper used --version as its installed-check and so
    reported "crispr-bean not installed. Install with: pip install
    crispr-bean" against a working installation -- a false negative that is
    indistinguishable, to the reader, from a real missing dependency.
    """
    exe = shutil.which("bean")
    if not exe:
        return False, "not on PATH"
    try:
        r = subprocess.run([exe, "--help"], capture_output=True, text=True,
                            timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        return False, f"{type(e).__name__}: {e}"
    if r.returncode != 0:
        return False, f"`bean --help` exited {r.returncode}"
    subs = ""
    for line in (r.stdout or "").splitlines():
        if "{" in line and "count" in line:
            subs = line.strip()
            break
    return True, f"{exe} ({subs[:70]})" if subs else exe


def write_bean_tables(out: Path, bins: "list[dict]", per_bin: "dict[str, dict]",
                       act: "dict[str, dict]", lib: dict,
                       editor: str) -> dict:
    """Write the three CSVs (plus edits) that `bean create-screen` consumes.

    This is the intended division of labour. IGVFagent contributes the
    base-edit-aware guide assignment -- the part that takes read assignment
    from 36.7% to 62.5% on this screen -- and BEAN contributes the Bayesian
    variant model that is deliberately not reimplemented here. create-screen
    is the seam: gRNA info, sample info, a guide x sample count matrix, and
    optionally per-guide edit counts, which is exactly what counting already
    produced.
    """
    out.mkdir(parents=True, exist_ok=True)
    samples = []
    for b in bins:
        sid = f"rep{b['rep']}_{b['side']}{b['pct']}"
        # BEAN's sorting model needs each bin's position on the sorted
        # phenotype axis as quantiles, not a label: a bottom-20% bin spans
        # [0.0, 0.2] and a top-20% bin spans [0.8, 1.0]. Without these it
        # cannot place the bins relative to each other.
        frac = (b["pct"] or 0) / 100.0
        if b["side"] == "bottom":
            lo, hi = 0.0, frac
        elif b["side"] == "top":
            lo, hi = 1.0 - frac, 1.0
        else:                       # bulk / unsorted spans everything
            lo, hi = 0.0, 1.0
        samples.append({"sample_id": sid, "replicate": b["rep"],
                         "condition": condition_label(b["side"], b["pct"]),
                         "sorting_bin": b["side"], "bin_pct": b["pct"],
                         "lower_quantile": lo, "upper_quantile": hi,
                         "accession": b["accession"]})
    sample_ids = [r["sample_id"] for r in samples]

    guides = sorted(lib["target_of"])
    info_path = out / "bean_gRNA_info.csv"
    with info_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "target", "type", "editing_activity", "editor"])
        for g in guides:
            a = act.get(g, {})
            w.writerow([g, lib["target_of"].get(g, g),
                         lib["type_of"].get(g, ""),
                         "" if a.get("activity") is None else round(a["activity"], 5),
                         editor])

    samples_path = out / "bean_sample_info.csv"
    with samples_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(samples[0]))
        w.writeheader()
        for r in samples:
            w.writerow(r)

    counts_path = out / "bean_gRNA_counts.csv"
    edits_path = out / "bean_edit_counts.csv"
    with counts_path.open("w", newline="") as cf, \
         edits_path.open("w", newline="") as ef:
        cw, ew = csv.writer(cf), csv.writer(ef)
        cw.writerow(["name"] + sample_ids)
        ew.writerow(["name"] + sample_ids)
        for g in guides:
            cw.writerow([g] + [per_bin.get(b["accession"], {})
                                .get("counts", Counter()).get(g, 0)
                                for b in bins])
            ew.writerow([g] + [per_bin.get(b["accession"], {})
                                .get("edited", Counter()).get(g, 0)
                                for b in bins])
    return {"gRNA_info": info_path, "sample_info": samples_path,
             "gRNA_counts": counts_path, "edit_counts": edits_path,
             "n_guides": len(guides), "n_samples": len(sample_ids),
             "conditions": sorted({r["condition"] for r in samples})}


def condition_label(side: str, pct) -> str:
    """BEAN's condition name for a bin. An unsorted bin has no percentage,
    so f"{side}{pct}" spelled it "bulkNone" -- which still matched the
    control-condition test by luck, and would have appeared verbatim in
    BEAN's own output tables."""
    return side if side == "bulk" or pct is None else f"{side}{pct}"


def control_condition(conditions) -> "tuple[str | None, str]":
    """Pick the sample condition `bean run sorting` should normalise against.

    BEAN defaults --control-condition to "bulk" and raises

        ValueError: No sample has control label `bulk`
        (set by `--control-condition`) in ReporterScreen.samples[condition]

    when none is present. That reads like a misconfiguration, but on a real
    screen it is usually a property of the data: the sorting model needs an
    unsorted sample to anchor each guide's baseline abundance, and a screen
    that only sequenced its tails never measured one. Of the two base-editing
    screens reachable here, 18loci_uptake (IGVFDS6464SOVZ) sorts
    bottom20/bottom40/top20/top40 and has no unsorted bin, while
    0426_LDLR137-219 (IGVFDS5542IBUS) does have one -- so this is a per-screen
    fact to report, not a default to paper over. Naming a bin "bulk" that is
    not one would make BEAN produce numbers whose baseline is a tail.
    """
    conds = sorted(set(conditions or []))
    control = next((c for c in conds if c.lower().startswith("bulk")), None)
    if control:
        return control, ""
    return None, (
        "this screen has no unsorted/bulk bin, and BEAN's sorting model needs "
        "one as --control-condition to anchor baseline guide abundance. "
        f"Conditions present: {conds or 'none'}. The counts and the BEAN "
        "screen object were still written, so `bean run` can be pointed at "
        "them once a reference sample is available; naming a tail bin as the "
        "control would give BEAN a baseline that is itself selected.")


def run_bean(out: Path, tables: dict, mode: str = "variant",
              screen_type: str = "sorting", n_iter: int = 0,
              timeout: int = 7200) -> dict:
    """`bean create-screen` then `bean run`, reporting each step honestly."""
    exe = shutil.which("bean")
    if not exe:
        return {"ran": False, "why": "bean not on PATH"}
    steps = []
    prefix = out / "bean_screen"
    cmds = [
        ("create-screen",
         [exe, "create-screen",
          str(tables["gRNA_info"]), str(tables["sample_info"]),
          str(tables["gRNA_counts"]),
          "-e", str(tables["edit_counts"]),
          "-o", str(prefix)]),
    ]
    h5 = Path(f"{prefix}.h5ad")
    for name, cmd in cmds:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        steps.append({"step": name, "exit_code": r.returncode,
                       "stderr": (r.stderr or "").strip()[-400:],
                       "cmd": " ".join(cmd)})
        if r.returncode != 0:
            return {"ran": False, "why": f"{name} exited {r.returncode}",
                     "steps": steps}
    if not h5.exists():
        return {"ran": False, "why": f"{name} exited 0 but wrote no {h5.name}",
                 "steps": steps}
    # `bean run` takes TWO positionals -- {sorting,survival} then
    # {variant,tiling} -- and the column names must be told to it, since our
    # table is not BEAN's own count-samples output. The first attempt passed
    # only "variant" (following the README's abbreviated example) and got
    # "invalid choice: 'variant' (choose from 'sorting', 'survival')".
    #
    # --guide-activity-col is the point of the whole exercise: it is where
    # BEAN takes a per-guide editing rate for its activity normalisation, and
    # the self-edit rate measured during counting goes straight into it.
    # BEAN defaults --control-condition to "bulk" and raises
    # "No sample has control label `bulk`" when none exists. That is a real
    # property of some screens, not a misconfiguration: the 18loci_uptake
    # screen has bottom20/bottom40/top20/top40 and no unsorted bin at all,
    # so its sorting model has nothing to anchor baseline abundance to. Say
    # that plainly instead of passing the ValueError through.
    control, why = control_condition(tables.get("conditions"))
    if not control:
        return {"ran": False, "steps": steps,
                 "screen_h5ad": str(h5), "why": why}
    # Two invocations, tried in order, because BEAN's models differ in what
    # they require of the input -- and the one that consumes editing activity
    # requires something IGVF does not publish.
    #
    #   MixtureNormal (default)  consumes --guide-activity-col, and its data
    #       class reads screen.layers["X_bcmatch"] unconditionally
    #       (bean/preprocessing/data_class.py:309, via
    #       VariantSortingReporterScreenData). X_bcmatch is the count of reads
    #       assigned by GUIDE BARCODE, which IGVF's guide libraries do not
    #       carry -- so this path dies with KeyError: 'X_bcmatch'.
    #
    #   Normal (--uniform-edit)  guards that read behind use_bcmatch
    #       (data_class.py:1120), so --ignore-bcmatch makes it runnable. But
    #       it raises "Can't use the guide activity column while constraining
    #       uniform edit" if --guide-activity-col is passed, and writes
    #       edit_eff = 1.0 for every guide.
    #
    # So on IGVF-published base-editing data, BEAN's Bayesian model is
    # available WITHOUT activity normalisation, or not at all. Both facts are
    # recorded rather than one being hidden: the fallback runs, and the result
    # is labelled with the model that produced it so nobody reads a Normal
    # posterior as an activity-normalised one. IGVFagent's own log2_per_edit
    # stays the activity-aware view, and BEAN contributes the posterior.
    attempts = [
        ("activity", ["--guide-activity-col", "editing_activity"]),
        ("uniform-edit fallback", ["--uniform-edit", "--ignore-bcmatch"]),
    ]
    base = [exe, "run", screen_type, mode, str(h5),
             "--control-condition", control,
             "--replicate-col", "replicate",
             "--condition-col", "condition",
             "--target-col", "target",
             "--sorting-bin-lower-quantile-col", "lower_quantile",
             "--sorting-bin-upper-quantile-col", "upper_quantile"]
    for label, extra in attempts:
        run_dir = out / ("bean_run" if label == "activity"
                          else "bean_run_uniform")
        cmd = base + extra + ["--outdir", str(run_dir)]
        if n_iter:
            cmd += ["--n-iter", str(n_iter)]
        run_dir.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=timeout)
        err = (r.stderr or "").strip()
        steps.append({"step": f"run {screen_type} {mode} [{label}]",
                       "exit_code": r.returncode,
                       "stderr": err[-600:], "cmd": " ".join(cmd)})
        if r.returncode == 0:
            res = {"ran": True, "screen_h5ad": str(h5), "steps": steps,
                    "model": "MixtureNormal" if label == "activity"
                              else "Normal",
                    "activity_normalised": label == "activity",
                    "outdir": str(run_dir), "why": ""}
            if label != "activity":
                res["caveat"] = (
                    "BEAN's activity-normalised MixtureNormal model needs the "
                    "X_bcmatch layer -- read counts assigned by guide barcode "
                    "-- which IGVF's guide libraries do not publish. These "
                    "posteriors come from BEAN's Normal model under "
                    "--uniform-edit, which assumes every guide edits with "
                    "efficiency 1.0 and refuses --guide-activity-col. For the "
                    "activity-aware view use this tool's own log2_per_edit, "
                    "which is reported alongside.")
            res["results"] = read_bean_results(run_dir)
            return res
        # Only fall through on the bcmatch barrier. Any other failure is a
        # real error and must not be masked by a weaker model quietly
        # succeeding in its place.
        if "X_bcmatch" not in err:
            return {"ran": False, "screen_h5ad": str(h5), "steps": steps,
                     "why": f"bean run exited {r.returncode}"}
    return {"ran": False, "screen_h5ad": str(h5), "steps": steps,
             "why": "bean run failed under both the activity and the "
                     "uniform-edit model"}


def read_bean_results(run_dir: Path) -> dict:
    """BEAN's per-target posteriors, if it wrote any.

    BEAN names the file after the model it fitted and nests it under a
    directory named after the input, so the path is not knowable in advance
    -- glob for it rather than construct it.
    """
    hits = sorted(run_dir.glob("**/bean_element_result.*.csv"))
    if not hits:
        return {"n_targets": 0, "path": None}
    rows = []
    with open(hits[0], newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                rows.append({"target": row.get("target", ""),
                              "n_guides": int(row.get("n_guides") or 0),
                              "mu": float(row["mu"]),
                              "mu_sd": float(row["mu_sd"]),
                              "mu_z": float(row["mu_z"])})
            except (KeyError, TypeError, ValueError):
                continue
    rows.sort(key=lambda r: -abs(r["mu_z"]))
    return {"n_targets": len(rows), "path": str(hits[0]), "top": rows[:15]}


# ─── Commands ───────────────────────────────────────────────────────────────

def _bean_command(accession: str, lib_file: str, editor: "Optional[str]") -> str:
    """The real BEAN invocation for this screen, for anyone who wants it.

    This tool does not reimplement `bean run`'s Bayesian model, so it should
    say how to get it rather than leave the impression it has been applied.
    """
    e = {"ABE": "A,G", "CBE": "C,T"}.get(editor or "", "A,G")
    return (f"bean count-samples --input sample_list.csv "
            f"-b {e.split(',')[0]} -f -r "
            f"--guide-info {lib_file}.csv --output-prefix {accession}\n"
            f"  bean qc {accession}.h5ad -o {accession}.masked.h5ad\n"
            f"  bean run variant {accession}.masked.h5ad --scale-by-acc")


def cmd_discover(args: argparse.Namespace) -> int:
    scr = cs.discover_screen(args.accession)
    if "error" in scr:
        print(scr["error"])
        return 2
    lib = load_library(args.accession)
    print(f"Screen:     {scr['series']}  ({len(scr['bins'])} libraries, "
          f"replicates {scr['replicates']})")
    print(f"Bins:       {', '.join(cs._bin_label(b) for b in scr['bin_kinds'])}")
    if not lib["resolved"]:
        print(f"Library:    UNRESOLVED — {lib['why']}")
        return 2
    print(f"Library:    {lib['file']}  ({lib['n_guides']:,} guides)")
    print(f"Editor:     {lib['editor'] or 'not stated in guide names'}"
          f"  (from the library's own guide names)")
    types = Counter(v for v in lib["type_of"].values() if v)
    print(f"Classes:    {dict(types)}")
    idx = lib["index"]
    cols = {c["column"] for c in idx["candidates"]}
    print(f"Sequence columns present: {', '.join(sorted(cols))}")
    for need, why in (("barcode", "BEAN resolves masked-sequence collisions "
                                   "with a guide barcode read from R2"),
                       ("reporter", "BEAN counts reporter alleles for "
                                     "bystander/tiling analysis")):
        have = need in cols
        print(f"  {need:9} {'present' if have else 'NOT PUBLISHED'} — {why}")
    ok, detail = bean_available()
    print(f"\nReal `bean` on this host: "
          f"{'YES — ' + detail if ok else 'NO — ' + detail}")
    if ok:
        # Whether `bean run` can model THIS screen is a property of its bins,
        # not of the install, and it is knowable here -- before a counting
        # run spends an hour to end in a decline. `bean run sorting` needs an
        # unsorted sample as its --control-condition; a screen that sequenced
        # only its tails never measured one.
        ctl, why = control_condition(
            [condition_label(b["side"], b["pct"]) for b in scr["bins"]])
        if ctl:
            print("  `bean analyze --run-bean` will hand the counts to it "
                  f"(create-screen + run), with --control-condition {ctl}.")
        else:
            print("  `bean analyze --run-bean` will write BEAN's screen "
                  "object but STOP before `bean run`:")
            print(f"    {why}")
    print(f"\nFull BEAN pipeline for this screen:\n  "
          f"{_bean_command(args.accession, lib['file'], lib['editor'])}")
    return 0


def cmd_count(args: argparse.Namespace) -> int:
    setup_logging()
    lib = load_library(args.accession)
    if not lib["resolved"]:
        print(lib["why"])
        return 2
    spacers = lib["spacer"]["seq_to_guide"]
    cal = calibrate_editor(spacers, args.accession, args.calibrate_reads)
    print(f"Library:    {lib['file']}  ({lib['n_guides']:,} guides)")
    print(f"Editor named in library: {lib['editor'] or 'none'}")
    if cal["tested"]:
        print(f"Measured on {cal['reads']:,} reads "
              f"(exact match alone: {cal['exact']:.1%}):")
        for t in cal["tested"]:
            print(f"    {t['editor']}  recovers {t['recovered']:>6.1%} more "
                  f"-> {t['total']:>6.1%} total, {t['ambiguous']:>5.1%} ambiguous, "
                  f"{t['masked_collisions']} masked collisions")
    editor = args.editor or cal["chosen"] or lib["editor"]
    if not editor:
        print("No base editor could be determined from the library names or "
              "the reads. Pass --editor ABE or --editor CBE.")
        return 3
    if lib["editor"] and cal["chosen"] and lib["editor"] != cal["chosen"]:
        print(f"  WARNING: the library names say {lib['editor']} but the reads "
              f"favour {cal['chosen']}. Using {editor}; check the library.")
    print(f"Using editor: {editor}")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    setup_logging()
    scr = cs.discover_screen(args.accession)
    if "error" in scr:
        print(scr["error"])
        return 2
    lib = load_library(args.accession)
    if not lib["resolved"]:
        print(lib["why"])
        return 2
    spacers = lib["spacer"]["seq_to_guide"]
    cal = calibrate_editor(spacers, args.accession, args.calibrate_reads)
    editor = args.editor or cal["chosen"] or lib["editor"]
    if not editor:
        print("No base editor determined; pass --editor ABE or --editor CBE.")
        return 3

    bins = [b for b in scr["bins"]
            if b["side"] in ("bottom", "top") and b["pct"] == args.tail]
    if not bins:
        print(f"No bin matches --tail {args.tail}. Present: "
              f"{', '.join(cs._bin_label(b) for b in scr['bin_kinds'])}")
        return 3
    print(f"Screen:     {scr['series']}  ({len(bins)} of {len(scr['bins'])} "
          f"libraries used for bottom{args.tail}%/top{args.tail}%)")
    print(f"Library:    {lib['file']}  ({lib['n_guides']:,} guides)   "
          f"editor: {editor}")
    if cal["tested"]:
        print(f"Assignment measured on {cal['reads']:,} reads: exact "
              f"{cal['exact']:.1%}, with {editor} masking "
              f"{next(t['total'] for t in cal['tested'] if t['editor']==editor):.1%}")

    plain = rp.build_matcher(spacers)
    masked = rp.build_masked_matcher(spacers, editor)
    per_bin: "dict[str, dict]" = {}
    for b in bins:
        d = count_library(b["accession"], masked, plain, args.max_reads)
        per_bin[b["accession"]] = d
        dp = d["disposition"]
        n = max(dp.get("reads", 1), 1)
        print(f"  Rep{b['rep']} {cs._bin_label(b):<12} {dp.get('reads',0):>8,} reads: "
              f"exact {dp.get('exact',0)/n:>5.1%}, self-edited "
              f"{dp.get('self_edited',0)/n:>5.1%}, ambiguous "
              f"{dp.get('ambiguous',0)/n:>5.1%}, unassigned "
              f"{dp.get('unassigned',0)/n:>5.1%}")

    act = guide_activity(per_bin)
    asum = activity_summary(act, lib["type_of"])
    rows, mod = score_screen(per_bin, bins, act, lib["target_of"],
                              lib["type_of"], args.tail, args.min_count)
    if not rows:
        print("\nNothing scored: no replicate had a complete tail pair with "
              "enough reads.")
        return 3

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_{scr['series'][:40]}"
    out = OUT_DIR / label
    out.mkdir(parents=True, exist_ok=True)
    cols = ["guide", "target", "guide_type", "n_reps", "reads", "activity",
            "activity_estimable", "log2_raw", "log2_per_edit", "sd",
            "sd_used", "t", "df", "p_value", "fdr"]
    with (out / "guide_effects.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with (out / "guide_activity.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["guide", "guide_type", "reads", "edited_reads", "activity",
                     "mean_edits_per_edited_read", "activity_estimable"])
        for g, a in sorted(act.items()):
            w.writerow([g, lib["type_of"].get(g, ""), a["reads"],
                         a["edited_reads"], a["activity"],
                         round(a["mean_edits_per_edited_read"], 3),
                         a["activity_estimable"]])

    ctrl = [r for r in rows
            if "positive control" in (r["guide_type"] or "").lower()]
    summary = {
        "screen": scr["series"], "query_set": args.accession,
        "editor": editor,
        "editor_from_library_names": lib["editor"],
        "editor_from_reads": cal.get("chosen"),
        "library": lib["file"], "guides": lib["n_guides"],
        "libraries_counted": len(bins),
        "tail_compared": f"bottom{args.tail}% vs top{args.tail}%",
        "replicates_used": mod["replicates_used"],
        "assignment_exact_only": round(cal.get("exact", 0.0), 4),
        "assignment_with_masking": round(
            next((t["total"] for t in cal["tested"] if t["editor"] == editor),
                 0.0), 4),
        "activity": asum,
        "guides_scored": len(rows),
        "guides_low_activity_under_5pct": mod["n_low_activity"],
        "positive_controls_scored": len(ctrl),
        "positive_controls_significant": sum(1 for r in ctrl if r["fdr"] < 0.05),
        "significant_fdr_0.05": sum(1 for r in rows if r["fdr"] < 0.05),
        "test": mod["test_basis"],
        "not_implemented": [
            "BEAN `run` Bayesian variant/tiling model with accessibility "
            "covariates -- this scores tail enrichment instead",
            "reporter-allele / bystander analysis -- the library's reporter "
            "column is empty for all guides",
            "bcmatch/semimatch split -- the library publishes no guide barcode, "
            "so masked-sequence collisions stay ambiguous",
        ],
        "full_bean_pipeline": _bean_command(args.accession, lib["file"], editor),
    }
    # Hand the base-edit-aware counts to the real BEAN, when asked for and
    # available. This is the point of the split: our mapping, BEAN's model.
    ok, detail = bean_available()
    summary["real_bean_available"] = ok
    summary["real_bean"] = detail
    bean_result = None
    if args.run_bean:
        if not ok:
            print(f"\n  --run-bean requested but the real bean is not usable: "
                  f"{detail}\n  Install it with: bash Deploy/install-bean.sh "
                  f"(or rebuild with IGVF_INSTALL_CRISPR_BEAN=1)")
            summary["bean_run"] = {"ran": False, "why": detail}
        else:
            # BEAN's sorting model needs an unsorted reference to anchor
            # baseline guide abundance, so the bulk bins go in even though
            # the tail comparison does not use them. They have to be COUNTED
            # for that, which the tail-only pass skipped.
            bulk = [b for b in scr["bins"] if b["side"] == "bulk"]
            for b in bulk:
                if b["accession"] in per_bin:
                    continue
                d = count_library(b["accession"], masked, plain, args.max_reads)
                per_bin[b["accession"]] = d
                dp = d["disposition"]; n = max(dp.get("reads", 1), 1)
                print(f"  Rep{b['rep']} bulk (for BEAN)  "
                      f"{dp.get('reads',0):>8,} reads: exact "
                      f"{dp.get('exact',0)/n:>5.1%}, self-edited "
                      f"{dp.get('self_edited',0)/n:>5.1%}")
            tables = write_bean_tables(out, bins + bulk, per_bin, act, lib,
                                        editor)
            print(f"\n  Wrote BEAN inputs: {tables['n_guides']:,} guides x "
                  f"{tables['n_samples']} samples")
            bean_result = run_bean(out, tables, mode=args.bean_mode,
                                    screen_type=args.bean_screen_type,
                                    n_iter=args.bean_iter)
            summary["bean_run"] = bean_result
            for st in bean_result.get("steps", []):
                flag = "ok" if st["exit_code"] == 0 else f"FAILED({st['exit_code']})"
                print(f"    bean {st['step']:16} {flag}")
                if st["exit_code"] != 0 and st.get("stderr"):
                    print(f"      {st['stderr'].splitlines()[-1][:150]}")
            if bean_result["ran"]:
                print(f"  BEAN model ran: {out / 'bean_run'}")
            else:
                print(f"  BEAN did not complete: {bean_result['why']}")

    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    plots = make_plots(out, rows, act, lib["type_of"], scr["series"])

    print()
    for k, v in summary.items():
        if k in ("activity", "not_implemented", "full_bean_pipeline"):
            continue
        print(f"  {k}: {v}")
    print(f"  activity: controls {asum['positive_controls']}")
    print(f"            variants {asum['variants']}")
    if asum["positive_controls"] and asum["variants"]:
        pc = asum["positive_controls"]["median"]
        vr = asum["variants"]["median"]
        print(f"  -> median editing activity {pc:.1%} at positive controls vs "
              f"{vr:.1%} at variants")
    print("\n  NOT reimplemented here:")
    for n in summary["not_implemented"]:
        print(f"    - {n}")
    print(f"\n  For the full model, run BEAN itself:\n    "
          f"{summary['full_bean_pipeline']}")
    withheld = sum(1 for r in rows if r["log2_per_edit"] is None)
    if withheld:
        print(f"\n  /edit withheld for {withheld:,} of {len(rows):,} guides: "
              f"editing activity below {PER_EDIT_ACTIVITY_FLOOR:.0%}, where "
              f"dividing by it reports the divisor's noise, not an effect.")
    br = summary.get("bean_run") or {}
    if br.get("ran"):
        res = br.get("results") or {}
        print(f"\n  BEAN's own model ran: {br['model']}"
              f"{'' if br.get('activity_normalised') else ' (NOT activity-normalised)'}"
              f", {res.get('n_targets', 0):,} targets")
        if br.get("caveat"):
            print(f"    why: {br['caveat']}")
        for t in (res.get("top") or [])[:5]:
            print(f"    {t['target'][:34]:34} mu {t['mu']:>7.3f} "
                  f"z {t['mu_z']:>7.3f}  ({t['n_guides']} guides)")
        if res.get("path"):
            print(f"    full table: {res['path']}")
    print(f"\nTop 15 by significance:")
    print(f"  {'target':30} {'log2':>7} {'act':>6} {'/edit':>8} {'fdr':>9}")
    for r in rows[:15]:
        pe = "-" if r["log2_per_edit"] is None else f"{r['log2_per_edit']:.2f}"
        ac = "-" if r["activity"] is None else f"{r['activity']:.2f}"
        print(f"  {r['target'][:30]:30} {r['log2_raw']:>7.2f} {ac:>6} "
              f"{pe:>8} {r['fdr']:>9.2e}")
    print(f"\nEffects:  {out / 'guide_effects.tsv'}")
    print(f"Activity: {out / 'guide_activity.tsv'}")
    for pth in plots:
        print(f"Plot: {pth}")
    print(f"Output: {out}")
    return 0


def make_plots(out: Path, rows: "list[dict]", act: "dict[str, dict]",
                type_of: "dict[str, str]", series: str) -> "list[Path]":
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return []
    made = []
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    a_ctrl = [a["activity"] for g, a in act.items() if a["activity"] is not None
              and "positive control" in (type_of.get(g, "") or "").lower()]
    a_var = [a["activity"] for g, a in act.items() if a["activity"] is not None
             and "positive control" not in (type_of.get(g, "") or "").lower()]
    if a_var:
        ax[0].hist(a_var, bins=50, color="#4C72B0", alpha=.85,
                   label=f"variants (n={len(a_var)})")
    if a_ctrl:
        ax[0].hist(a_ctrl, bins=30, color="#C44E52", alpha=.85,
                   label=f"positive controls (n={len(a_ctrl)})")
    ax[0].set_xlabel("editing activity (edited reads / assigned reads)")
    ax[0].set_ylabel("guides")
    ax[0].set_title("Self-editing activity per guide")
    ax[0].legend(fontsize=8)

    xs = [r["activity"] for r in rows if r["activity"] is not None]
    ys = [abs(r["log2_raw"]) for r in rows if r["activity"] is not None]
    ax[1].scatter(xs, ys, s=8, alpha=.35, c="#55A868")
    ax[1].set_xlabel("editing activity")
    ax[1].set_ylabel("|log2 raw effect|")
    ax[1].set_title("Effect against activity\n(low activity = low power, "
                     "not evidence of no effect)")

    eff = np.array([r["log2_raw"] for r in rows])
    fdr = np.array([max(r["fdr"], 1e-12) for r in rows])
    isc = np.array(["positive control" in (r["guide_type"] or "").lower()
                    for r in rows])
    ax[2].scatter(eff[~isc], -np.log10(fdr[~isc]), s=10, c="#4C72B0",
                  alpha=.6, label="variants")
    if isc.any():
        ax[2].scatter(eff[isc], -np.log10(fdr[isc]), s=22, c="#C44E52",
                      alpha=.9, label="positive controls")
    ax[2].axhline(-math.log10(0.05), ls="--", lw=1, c="k")
    ax[2].set_xlabel("log2( bottom / top )")
    ax[2].set_ylabel("-log10 FDR")
    ax[2].set_title("Enrichment")
    ax[2].legend(fontsize=8)
    fig.suptitle(f"Base-editing screen — {series}", fontsize=10)
    fig.tight_layout()
    p = out / "base_editing_qc.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    made.append(p)
    return made


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="base_editing_screen",
        description="Base-editing screens, BEAN-style guide assignment.")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("discover", help="Screen, library, editor and what "
                                        "BEAN inputs IGVF publishes.")
    d.add_argument("accession")
    c = sub.add_parser("count", help="Measure which editor's masking "
                                      "recovers reads.")
    c.add_argument("accession")
    c.add_argument("--editor", choices=sorted(rp.BASE_EDITS))
    c.add_argument("--calibrate-reads", type=int, default=40000)
    a = sub.add_parser("analyze", help="Count all bins, estimate activity, "
                                        "score with activity reported.")
    a.add_argument("accession")
    a.add_argument("--editor", choices=sorted(rp.BASE_EDITS))
    a.add_argument("--tail", type=int, default=20)
    a.add_argument("--min-count", type=int, default=10)
    a.add_argument("--max-reads", type=int, default=None)
    a.add_argument("--calibrate-reads", type=int, default=40000)
    a.add_argument("--run-bean", action="store_true",
                    help="Also hand the counts to the real `bean` "
                         "(create-screen + run) for its Bayesian model.")
    a.add_argument("--bean-mode", default="variant",
                    choices=["variant", "tiling"],
                    help="BEAN library design: variant ignores bystander "
                         "edits, tiling models them.")
    a.add_argument("--bean-screen-type", default="sorting",
                    choices=["sorting", "survival"],
                    help="BEAN selection type. A FACS screen is 'sorting'.")
    a.add_argument("--bean-iter", type=int, default=0,
                    help="Override BEAN's --n-iter (0 = its default).")
    a.add_argument("--label")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    return {"discover": cmd_discover, "count": cmd_count,
            "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
