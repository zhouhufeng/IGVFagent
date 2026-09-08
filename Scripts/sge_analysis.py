#!/usr/bin/env python3
"""Saturation genome editing (SGE): variant counts and functional scores.

SGE is not RNA-seq, and quantifying its reads against a transcriptome
answers a different question: on IGVFDS4629JYPY it reported that 99.86% of
counts fall on PALB2, which is the amplicon design restated rather than a
result. What an SGE experiment measures is the *fate of each programmed
variant* -- variants that damage the protein are depleted from the cell
population over time, and the readout is that depletion.

So the analysis is:

1. read the editing-template design (amplicon coordinates per target),
2. fetch the reference amplicon sequence,
3. anchor each read to it and call the substitutions it carries,
4. count variants per sample,
5. compare a LATE timepoint against an EARLY one -- log2(late/early) --
   which is the functional score.

A single timepoint yields counts and cannot yield scores; the commands say
so rather than presenting counts as if they measured function.

Everything is stdlib + numpy/pandas/matplotlib. No aligner: an SGE amplicon
is a fixed-length window whose reads start at a known anchor, so
position-wise comparison is both sufficient and exact, and it is easier to
reason about than a general aligner's gap placement.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import re
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import raw_data_pipeline as rp                                # noqa: E402

ROOT = rp.ROOT
OUT_DIR = ROOT / "Docs" / "SGE"
WORK_DIR = ROOT / "Data" / "SGE"
REF_CACHE = WORK_DIR / "_reference"
LOG_DIR = ROOT / "Docs" / "Logs"

ENSEMBL = os.environ.get("ENSEMBL_REST", "https://rest.ensembl.org")

# Content type of the design file on a construct library set.
_TEMPLATE_CONTENT = "editing templates"


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"sge_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)])
    return log


# ─── Design ─────────────────────────────────────────────────────────────────

def find_design_files(accession: str) -> "list[dict]":
    """Editing-template TSVs for a MeasurementSet, via its library set."""
    st, ms = rp.portal_json(
        f"/search/?type=FileSet&accession={accession}&format=json")
    rows = (ms or {}).get("@graph") or []
    if not rows:
        return []
    st, full = rp.portal_json(f"{rows[0].get('@id')}?format=json")
    out = []
    for cls in (full or {}).get("construct_library_sets") or []:
        for icf in cls.get("integrated_content_files") or []:
            if _TEMPLATE_CONTENT in str(icf.get("content_type", "")).lower():
                out.append({"accession": icf.get("accession"),
                             "id": icf.get("@id"),
                             "library": cls.get("accession")})
    return out


def load_design(file_id: str) -> "list[dict]":
    """Parse an editing-template TSV into per-target dicts."""
    st, obj = rp.portal_json(
        f"{file_id}?format=json" if file_id.startswith("/")
        else f"/tabular-files/{file_id}/?format=json")
    href = (obj or {}).get("href")
    if not href:
        return []
    dest = WORK_DIR / "_design" / Path(href).name
    if not dest.exists():
        rp.portal_download(href, dest)
    text = rp.read_text_maybe_gzip(dest)
    lines = text.splitlines()
    if not lines:
        return []
    hdr = [h.strip() for h in lines[0].split("\t")]
    out = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) >= 3:
            out.append(dict(zip(hdr, parts)))
    return out


def parse_required_edits(spec: str) -> "dict[int, str]":
    """``23626377C,23626344C`` -> {23626377: 'C', ...}.

    These are the fixed edits every template carries (PAM disruption and
    the like). They are present in ~every read by construction, so counting
    them as variants would put two meaningless entries at the top of every
    result table.
    """
    out: "dict[int, str]" = {}
    for tok in (spec or "").split(","):
        tok = tok.strip()
        m = re.match(r"^(\d+)([ACGT])$", tok, re.I)
        if m:
            out[int(m.group(1))] = m.group(2).upper()
    return out


# ─── Reference ──────────────────────────────────────────────────────────────

def fetch_reference(chrom: str, start: int, stop: int) -> str:
    """Amplicon reference sequence from Ensembl, cached on disk."""
    key = f"{chrom}_{start}_{stop}.txt"
    cached = REF_CACHE / key
    if cached.exists():
        return cached.read_text().strip().upper()
    c = chrom.replace("chr", "")
    url = (f"{ENSEMBL}/sequence/region/human/{c}:{start}..{stop}"
           f"?content-type=application/json")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        seq = json.load(r)["seq"].upper()
    REF_CACHE.mkdir(parents=True, exist_ok=True)
    cached.write_text(seq)
    return seq


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ─── Read parsing ───────────────────────────────────────────────────────────

def iter_reads(path: Path, limit: Optional[int] = None):
    """Yield sequences from a FASTQ(.gz), sniffing gzip by magic bytes."""
    op = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
    n = 0
    with op(path, "rt", errors="replace") as fh:            # type: ignore[operator]
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()
                n += 1
                if limit and n >= limit:
                    return


def call_variants(reads, reference: str, anchor_len: int = 25,
                   min_anchor_frac: float = 0.9) -> dict:
    """Count substitutions per (offset, alt) against the reference amplicon.

    Reads are placed by matching the reference's leading ``anchor_len`` bases,
    trying the reverse complement too, because a read pair puts one mate on
    each strand and silently dropping one of them halves the depth.

    Only substitutions are called. Indels in SGE amplicons are rare and a
    position-wise comparison cannot place them; a read whose identity to the
    reference falls below ``min_anchor_frac`` is counted as unparsed rather
    than being forced into a substitution list that would be wrong.
    """
    fwd_anchor = reference[:anchor_len]
    rev_ref = revcomp(reference)
    rev_anchor = rev_ref[:anchor_len]
    counts: "Counter[tuple[int, str]]" = Counter()
    depth = [0] * len(reference)
    stats = Counter()
    def _place(s: str) -> "Optional[tuple[int, int]]":
        """(offset in reference, offset in read) for a placeable read.

        Two directions are tried, because a read pair covers the two ENDS of
        the amplicon: R2 begins at the amplicon start, so the reference's
        leading anchor is found inside the read, while R1 begins somewhere
        downstream and it is the READ's leading anchor that must be found in
        the reference. Checking only the first placed half the depth -- 52%
        of reads were being discarded as unanchored on this dataset.
        """
        off = s.find(fwd_anchor)
        if off >= 0:
            return 0, off
        head = s[:anchor_len]
        rpos = reference.find(head)
        if rpos >= 0:
            return rpos, 0
        return None

    for seq in reads:
        stats["total"] += 1
        placed = None
        s = seq
        for cand in (seq, revcomp(seq)):
            placed = _place(cand)
            if placed is not None:
                s = cand
                break
        if placed is None:
            stats["no_anchor"] += 1
            continue
        rstart, sstart = placed
        s = s[sstart:]
        n = min(len(s), len(reference) - rstart)
        if n < anchor_len:
            stats["too_short"] += 1
            continue
        mism = [(rstart + i, s[i]) for i in range(n)
                if s[i] != reference[rstart + i] and s[i] in "ACGT"]
        if (n - len(mism)) / float(n) < min_anchor_frac:
            stats["low_identity"] += 1
            continue
        stats["parsed"] += 1
        for i in range(rstart, rstart + n):
            depth[i] += 1
        for i, alt in mism:
            counts[(i, alt)] += 1
    return {"counts": counts, "depth": depth, "stats": stats}


# ─── Commands ───────────────────────────────────────────────────────────────

def _resolve_target(design: "list[dict]", target: Optional[str],
                     accession: str) -> "Optional[dict]":
    """Pick the design row, inferring the target from the set alias if needed."""
    if target:
        for row in design:
            if row.get("target", "").upper() == target.upper():
                return row
        return None
    st, ms = rp.portal_json(
        f"/search/?type=FileSet&accession={accession}&format=json")
    rows = (ms or {}).get("@graph") or []
    alias = " ".join((rows[0].get("aliases") or []) if rows else [])
    # "lea-starita:PALB2-X7A-replicate2-day12" -> PALB2_X7A
    m = re.search(r"([A-Z0-9]+)-(X[0-9]+[A-Z]?)", alias, re.I)
    if m:
        want = f"{m.group(1)}_{m.group(2)}".upper()
        for row in design:
            if row.get("target", "").upper() == want:
                return row
    return None


def _reads_for(accession: str, max_gb: float) -> "list[Path]":
    files = [f for f in rp.list_files(accession)
             if str(f.get("file_format", "")).lower() == "fastq"]
    total = sum(float(f.get("file_size") or 0) for f in files) / 1e9
    if max_gb and total > max_gb:
        raise RuntimeError(f"{accession}: reads total {total:.2f} GB, over "
                           f"--max-download-gb {max_gb}")
    out = []
    rp.FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
    for f in files:
        name = Path(str(f.get("href") or f.get("accession"))).name
        dest = rp.FASTQ_CACHE / name
        if not dest.exists():
            print(f"Downloading {f.get('accession')} "
                  f"({rp.gb(f.get('file_size'))} GB)…")
            rp.portal_download(f["href"], dest)
        else:
            print(f"Cached      {f.get('accession')} "
                  f"({rp.gb(f.get('file_size'))} GB)")
        out.append(dest)
    return out


def count_one(accession: str, target: Optional[str], max_reads: Optional[int],
               max_gb: float) -> dict:
    """Variant counts for one sample."""
    designs = find_design_files(accession)
    if not designs:
        raise RuntimeError(
            f"{accession} has no editing-template design reachable: no "
            f"construct_library_sets with a '{_TEMPLATE_CONTENT}' file. "
            f"Without the design there is no amplicon to align to.")
    rows = load_design(designs[0]["id"] or designs[0]["accession"])
    row = _resolve_target(rows, target, accession)
    if row is None:
        names = ", ".join(sorted(r.get("target", "?") for r in rows)[:12])
        raise RuntimeError(
            f"could not determine which target {accession} covers. Pass "
            f"--target. Available: {names}")
    chrom = row["chrom"]
    amp_start, amp_stop = int(row["ampstart"]), int(row["ampstop"])
    reference = fetch_reference(chrom, amp_start, amp_stop)
    required = parse_required_edits(row.get("required_edits", ""))
    fastqs = _reads_for(accession, max_gb)
    agg = None
    for fq in fastqs:
        r = call_variants(iter_reads(fq, limit=max_reads), reference)
        if agg is None:
            agg = r
        else:
            agg["counts"].update(r["counts"])
            agg["depth"] = [a + b for a, b in zip(agg["depth"], r["depth"])]
            agg["stats"].update(r["stats"])
    return {"accession": accession, "target": row["target"], "chrom": chrom,
            "amp_start": amp_start, "amp_stop": amp_stop,
            "edit_start": int(row["editstart"]), "edit_stop": int(row["editstop"]),
            "reference": reference, "required_edits": required,
            "design_file": designs[0]["accession"], **agg}


def variant_table(res: dict, min_count: int = 5) -> "list[dict]":
    """Per-variant rows with genomic coordinates and frequency."""
    ref = res["reference"]
    out = []
    for (i, alt), n in sorted(res["counts"].items(),
                              key=lambda kv: -kv[1]):
        if n < min_count:
            continue
        pos = res["amp_start"] + i
        d = res["depth"][i] or 1
        out.append({
            "chrom": res["chrom"], "pos": pos, "ref": ref[i], "alt": alt,
            "count": n, "depth": d, "freq": round(n / d, 6),
            "offset": i,
            "in_edit_window": res["edit_start"] <= pos <= res["edit_stop"],
            "required_edit": required_hit(res, pos, alt),
        })
    return out


def required_hit(res: dict, pos: int, alt: str) -> bool:
    return res["required_edits"].get(pos) == alt


def _write_tsv(path: Path, rows: "list[dict]", cols: "list[str]") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def cmd_design(args: argparse.Namespace) -> int:
    designs = find_design_files(args.accession)
    if not designs:
        print(f"No editing-template design reachable from {args.accession}.")
        return 2
    print(f"Design files: "
          + ", ".join(f"{d['accession']} (library {d['library']})"
                      for d in designs))
    rows = load_design(designs[0]["id"] or designs[0]["accession"])
    print(f"{len(rows)} targets:")
    for r in rows:
        amp = int(r.get("ampstop", 0)) - int(r.get("ampstart", 0))
        print(f"  {r.get('target','?'):14} {r.get('chrom','?')}:"
              f"{r.get('ampstart','?')}-{r.get('ampstop','?')}  "
              f"amplicon {amp} bp  required_edits={r.get('required_edits','') or '-'}")
    return 0


def cmd_count(args: argparse.Namespace) -> int:
    res = count_one(args.accession, args.target, args.max_reads,
                     args.max_download_gb)
    rows = variant_table(res, args.min_count)
    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_{res['accession']}"
    out = OUT_DIR / label
    cols = ["chrom", "pos", "ref", "alt", "count", "depth", "freq",
            "offset", "in_edit_window", "required_edit"]
    tsv = _write_tsv(out / "variant_counts.tsv", rows, cols)
    st = res["stats"]
    summary = {
        "accession": res["accession"], "target": res["target"],
        "design_file": res["design_file"],
        "amplicon": f"{res['chrom']}:{res['amp_start']}-{res['amp_stop']}",
        "edit_window": f"{res['chrom']}:{res['edit_start']}-{res['edit_stop']}",
        "reads_total": st["total"], "reads_parsed": st["parsed"],
        "reads_no_anchor": st["no_anchor"], "reads_low_identity": st["low_identity"],
        "parsed_fraction": round(st["parsed"] / max(st["total"], 1), 4),
        "median_depth": int(sorted(res["depth"])[len(res["depth"]) // 2]),
        "variants_reported": len(rows),
        "variants_in_edit_window": sum(1 for r in rows if r["in_edit_window"]),
        "required_edits_seen": sum(1 for r in rows if r["required_edit"]),
        "min_count": args.min_count,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"Variant counts: {tsv}")
    print(f"\nNOTE: counts alone are not functional scores. Pair this with an "
          f"early timepoint and run `igvfagent sge score` for log2 enrichment.")
    return 0


def cmd_score(args: argparse.Namespace) -> int:
    late = count_one(args.late, args.target, args.max_reads, args.max_download_gb)
    early = count_one(args.early, args.target, args.max_reads, args.max_download_gb)
    if late["target"] != early["target"]:
        print(f"REFUSING: {args.late} covers {late['target']} but "
              f"{args.early} covers {early['target']} — different amplicons "
              f"cannot be compared.")
        return 3
    import numpy as np

    lt = {(i, a): n for (i, a), n in late["counts"].items()}
    et = {(i, a): n for (i, a), n in early["counts"].items()}
    keys = sorted(set(lt) | set(et))
    l_tot = max(sum(lt.values()), 1)
    e_tot = max(sum(et.values()), 1)
    rows = []
    for (i, alt) in keys:
        lc, ec = lt.get((i, alt), 0), et.get((i, alt), 0)
        if lc + ec < args.min_count:
            continue
        pos = late["amp_start"] + i
        # Frequencies, not raw counts: the two libraries are sequenced to
        # different depths and a raw ratio would measure that instead.
        lf = (lc + args.pseudocount) / l_tot
        ef = (ec + args.pseudocount) / e_tot
        rows.append({
            "chrom": late["chrom"], "pos": pos, "ref": late["reference"][i],
            "alt": alt, "late_count": lc, "early_count": ec,
            "late_freq": round(lf, 8), "early_freq": round(ef, 8),
            "log2_enrichment": round(float(np.log2(lf / ef)), 4),
            "in_edit_window": late["edit_start"] <= pos <= late["edit_stop"],
            "required_edit": required_hit(late, pos, alt),
        })
    # Required edits sit in every read by construction; they are not variants
    # under selection and would anchor the distribution if left in.
    scored = [r for r in rows if not r["required_edit"] and r["in_edit_window"]]
    label = args.label or (f"{time.strftime('%Y%m%d_%H%M%S')}_"
                            f"{late['target']}_score")
    out = OUT_DIR / label
    cols = ["chrom", "pos", "ref", "alt", "late_count", "early_count",
            "late_freq", "early_freq", "log2_enrichment", "in_edit_window",
            "required_edit"]
    _write_tsv(out / "all_variants.tsv", rows, cols)
    tsv = _write_tsv(out / "functional_scores.tsv", scored, cols)
    plots = _plots(out, late, early, scored)
    vals = [r["log2_enrichment"] for r in scored]
    summary = {
        "target": late["target"], "late": args.late, "early": args.early,
        "amplicon": f"{late['chrom']}:{late['amp_start']}-{late['amp_stop']}",
        "late_reads_parsed": late["stats"]["parsed"],
        "early_reads_parsed": early["stats"]["parsed"],
        "variants_scored": len(scored),
        "median_log2": round(float(np.median(vals)), 4) if vals else None,
        "depleted_lt_-1": sum(1 for v in vals if v < -1),
        "enriched_gt_1": sum(1 for v in vals if v > 1),
        "pseudocount": args.pseudocount, "min_count": args.min_count,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"Functional scores: {tsv}")
    for p in plots:
        print(f"Plot: {p}")
    print(f"Output: {out}")
    return 0


def _plots(out: Path, late: dict, early: dict,
            scored: "list[dict]") -> "list[Path]":
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return []
    out.mkdir(parents=True, exist_ok=True)
    made = []
    if scored:
        vals = np.array([r["log2_enrichment"] for r in scored])
        pos = np.array([r["pos"] for r in scored])
        fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
        ax[0].hist(vals, bins=40, color="#4C72B0")
        ax[0].axvline(0, color="k", lw=1, ls="--")
        ax[0].set_xlabel("log2(late / early)")
        ax[0].set_ylabel("variants")
        ax[0].set_title(f"{late['target']} functional score distribution\\n"
                        f"(negative = depleted = damaging)")
        sc = ax[1].scatter(pos, vals, c=vals, cmap="coolwarm_r", s=18)
        ax[1].axhline(0, color="k", lw=1, ls="--")
        ax[1].set_xlabel(f"{late['chrom']} position")
        ax[1].set_ylabel("log2 enrichment")
        ax[1].set_title("score along the amplicon")
        fig.colorbar(sc, ax=ax[1], label="log2")
        fig.tight_layout()
        p = out / "sge_scores.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        made.append(p)
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(range(len(late["depth"])), late["depth"], label="late", lw=1)
    ax.plot(range(len(early["depth"])), early["depth"], label="early", lw=1)
    ax.axvspan(late["edit_start"] - late["amp_start"],
               late["edit_stop"] - late["amp_start"],
               color="orange", alpha=.15, label="edit window")
    ax.set_xlabel("amplicon offset (bp)")
    ax.set_ylabel("reads")
    ax.set_title(f"{late['target']} coverage")
    ax.legend()
    fig.tight_layout()
    p = out / "sge_coverage.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    made.append(p)
    return made


def find_sibling_timepoints(accession: str) -> "list[dict]":
    """Other SGE sets on the same target, so an early mate can be found."""
    st, ms = rp.portal_json(
        f"/search/?type=FileSet&accession={accession}&format=json")
    rows = (ms or {}).get("@graph") or []
    alias = " ".join((rows[0].get("aliases") or []) if rows else [])
    m = re.search(r"([A-Za-z0-9]+)-(X[0-9]+[A-Za-z]?)-replicate(\d+)", alias)
    if not m:
        return []
    gene, tgt, rep = m.group(1), m.group(2), m.group(3)
    st, d = rp.portal_json(
        f"/search/?type=MeasurementSet&preferred_assay_titles=SGE"
        f"&limit=200&format=json")
    out = []
    for r in (d or {}).get("@graph") or []:
        a = " ".join(r.get("aliases") or [])
        if f"{gene}-{tgt}-replicate{rep}" in a and r.get("accession") != accession:
            day = re.search(r"day(\d+)", a)
            out.append({"accession": r.get("accession"), "alias": a,
                         "day": int(day.group(1)) if day else None})
    return sorted(out, key=lambda x: (x["day"] is None, x["day"]))


def cmd_analyze(args: argparse.Namespace) -> int:
    """One-shot: find the design, find the early mate, score."""
    print(f"Target set: {args.accession}")
    designs = find_design_files(args.accession)
    print(f"Design:     {designs[0]['accession'] if designs else 'NOT FOUND'}")
    if not designs:
        print("  Cannot proceed: no editing-template design is linked, so "
              "there is no amplicon reference to call variants against.")
        return 2
    early = args.early
    if not early:
        sibs = find_sibling_timepoints(args.accession)
        print(f"Siblings:   {', '.join(s['accession'] + ' (day' + str(s['day']) + ')' for s in sibs) or 'none found'}")
        cand = [s for s in sibs if s["day"] is not None]
        if cand:
            early = cand[0]["accession"]
            print(f"Early mate: {early} (day {cand[0]['day']}) — earliest "
                  f"timepoint on the same target and replicate")
    if not early:
        print("\nNo early timepoint found, so functional scores are not "
              "possible: SGE measures depletion BETWEEN timepoints. Counting "
              "this sample alone instead.")
        return cmd_count(args)
    args.late, args.early = args.accession, early
    return cmd_score(args)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sge_analysis",
        description="Saturation genome editing: variant counts and scores.")
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("--target", help="Design target, e.g. PALB2_X7A. "
                                          "Inferred from the set alias if omitted.")
        sp.add_argument("--max-reads", type=int, default=None,
                         help="Cap reads per FASTQ (for a quick look).")
        sp.add_argument("--max-download-gb", type=float, default=20.0)
        sp.add_argument("--min-count", type=int, default=5)
        sp.add_argument("--label")
        return sp

    d = sub.add_parser("design", help="List the editing-template targets.")
    d.add_argument("accession")

    c = common(sub.add_parser("count", help="Variant counts for one sample."))
    c.add_argument("accession")

    s = common(sub.add_parser("score", help="log2 enrichment, late vs early."))
    s.add_argument("--late", required=True)
    s.add_argument("--early", required=True)
    s.add_argument("--pseudocount", type=float, default=1.0)

    a = common(sub.add_parser("analyze",
                               help="Find design + early mate, then score."))
    a.add_argument("accession")
    a.add_argument("--early", help="Early timepoint; auto-discovered if omitted.")
    a.add_argument("--pseudocount", type=float, default=1.0)
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    for dd in (OUT_DIR, WORK_DIR, REF_CACHE):
        dd.mkdir(parents=True, exist_ok=True)
    return {"design": cmd_design, "count": cmd_count,
            "score": cmd_score, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
