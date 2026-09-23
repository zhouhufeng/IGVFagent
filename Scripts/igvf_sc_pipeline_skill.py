#!/usr/bin/env python3
"""IGVF single-cell multiome (SHARE-seq / 10x) preprocessing + QC toolkit (port of IGVF/single-cell-pipeline).

Port of https://github.com/IGVF/single-cell-pipeline (MIT, commit 202a21a, a WDL
workflow) together with the two Python wrappers it calls from
https://github.com/IGVF/atomic-workflows (tag v1.1: run_kallisto.py,
run_chromap.py).  Every Python / R / bash / WDL file was read; the analysis
steps were re-derived in Python (numpy + pandas + scipy, optional pysam /
anndata / h5py / matplotlib).  Relationship: port (algorithms translated, source
consulted, no code copied).  The aligners themselves (chromap, kallisto |
bustools) are wrapped: the argument construction is reproduced exactly and the
binary is run only when it is on PATH; otherwise the exact command is printed
and the subcommand exits 0.

Definitions, exactly as upstream computes them
  barcode orientation  first 100 000 reads of the barcode FASTQ, sequence[offset:]
                       looked up in the onlist and in its reverse complement;
                       proportions are divided by 100 000 (not by the reads
                       seen); valid if either >= 0.45; forward chosen when
                       direct >= revcomp.  Output: chromap read-format
                       "bc:<offset>:-1" or "bc:<offset>:-1:-".
  SHARE-seq correction last 99 bp of read 2 hold the three round barcodes in
                       windows [14:24], [52:62], [90:99]; each window is looked
                       up exact at [1:9], then 1-bp left shift [0:8], then 1-bp
                       right shift [2:] (+"N"/"F" when the window is 9 bp)
                       against a dictionary of every onlist 8-mer and all its
                       1-substitution neighbours (ACGTN).  Read 2 starting with
                       10 G -> poly_G_in_first_10bp (checked first); all three
                       windows containing GGGGGGGG -> poly_G_barcode; else
                       mismatch.  RNA: header "<name>_<r1>,<r2>,<r3>[,<pkr>]_<UMI10>",
                       read 2 := r1r2r3 + UMI; ATAC: read 2 := read2[:-99] and a
                       barcode FASTQ r1r2r3.  QC: <prefix>_barcode_qc.txt with
                       library, match, mismatch, poly_G_barcode, poly_G_in_first_10bp.
  dovetail trim        query = revcomp(read2[0:20]); position = read1.rfind(query),
                       else the first window (from the left) of read1 with
                       Levenshtein distance <= 1; both reads cut to [:pos+20].
                       Stats total_reads, untrimmed_reads, trimmed_reads, %trimmed
                       (each pair counts 2).
  TSS enrichment       promoter = TSS(start column) +/- flank (2000); both
                       fragment ends are counted if 0 <= pos-(TSS-flank) <=
                       2*flank-1 (fragments fetched by tabix overlap), strand '-'
                       mirrored.  Bulk: profile + pseudocount 1, normaliser =
                       mean of positions [0:100] and [size-101:size]; smoothed =
                       convolve(profile, ones(w), 'same')/w/normaliser (w = 20);
                       score = max(smoothed).  Per barcode (ArchR heuristic): 301
                       bins = left 100 bp, TSS +/- 50 bp (101), right 100 bp;
                       TSSe = 2*reads_tss/101 / max(0.2, 2*flank_reads/200).
  RNA metrics          scanpy.calculate_qc_metrics: total_counts, genes
                       (n_genes_by_counts); kb "nac": mature + ambiguous + nascent
                       layers summed column-wise (genes summed too, as upstream).
                       Barcodes suffixed "_<subpool>" unless subpool == "none".
  joint QC             RNA barcodes with UMIs >= remove_low_yielding_cells, ATAC
                       barcodes with reads/2 >= it; outer merge (missing -> 0);
                       QC = both / RNA only / ATAC only / neither from
                       (UMIs >= min_umis & genes >= min_genes) and (TSS >= min_tss
                       & frags >= min_frags); QC_count = "<QC> (<n>)".
  barcode rank knee    points (rank, log10 count) with count >= cutoff; elbow =
                       point farthest from the line through the endpoints; the
                       second (knee) point is searched on ranks 1..elbow: smooth
                       spline (R smooth.spline, spar = 1, re-implemented as the
                       same penalised cubic B-spline: .nknots.smspl knots, lambda
                       = ratio * 256^(3*spar-1)), second derivative; if its sign
                       is constant over the middle 80% use all points, else slice
                       between the last >= 0 before and the first >= 0 after its
                       minimum, then the same farthest-point rule.
  10x barcode mapping  when both onlists have the same length: rows
                       "<rna>\\t<revcomp(atac)>" then "<rna>\\t<atac>" (chromap
                       --barcode-translate orientation).
  chromap log -> json  "Number of ..." lines (not threads) and "#" lines, dots
                       removed; percentage_duplicates = 100*dups/(total -
                       unmapped - lowmapq) from the barcode summary, "%5.1f".
  kb count (nac)       kb count --workflow=nac -i idx/index.idx -g idx/t2g.txt
                       -c1 idx/cdna.txt -c2 idx/nascent.txt --sum=total -x FMT
                       -w onlist [-r replacement] --strand S -o OUT --h5ad -t T
                       <l1r1 l1r2 [l1bc] l2r1 ...>; subpool appended to
                       counts_unfiltered[_modified]/adata.h5ad /obs/barcode and
                       cells_x_genes.barcodes.txt; h5ad moved to OUT.h5ad.
  chromap              chromap -x idx/index --read-format FMT -r fasta
                       --remove-pcr-duplicates --remove-pcr-duplicates-at-cell-level
                       --trim-adapters --low-mem --BED|--SAM -l 2000
                       --bc-error-threshold 1 -t T --bc-probability-threshold 0.90
                       -q 30 --barcode-whitelist onlist [--barcode-translate map]
                       ...; subpool appended to fragment column 4 and summary
                       column 1; bgzip + tabix --zero-based --preset bed.
  sample FASTQs        first 10 000 000 reads (40 M lines), re-gzipped.

Deliberate deviations (the upstream behaviour is reproduced with --upstream-compat)
  * correct_fastq.get_barcodes binds r1/r2/r3 to ONE set (`a = b = c = set()`),
    so each round dictionary holds every round's barcodes; and a later
    barcode's 1-mismatch neighbour can overwrite an exact onlist entry (set
    order).  The port keeps the rounds separate and gives exact matches
    priority.
  * compute_tss_enrichment_bulk never fills the per-barcode counters (the
    ArchR bins are computed but not stored) and never writes its per-barcode
    TSV; the port fills them and writes <prefix>.tss_enrichment_barcode_stats.tsv.
  * update-file-annotations returns a 2-tuple for unknown suffixes and crashes
    unpacking three values; the port returns three Nones.
  * log_atac crashes on non-numeric values and divides by zero on an empty
    summary; the port keeps strings and writes null.
  * The WDL default subpool is the string "none", which run_chromap.py /
    run_kallisto.py test only for truthiness, so upstream appends "_none" to
    every barcode (and the bam awk always runs).  The port treats "none" as
    no subpool everywhere, matching qc_rna_extract_metrics.py.
  * sc.pp.calculate_qc_metrics' default percent_top needs >= 500 genes; the
    port computes the two needed columns directly (identical values).

Subcommands
  barcode-revcomp-detect  onlist orientation of the barcode read -> read-format string
  correct-fastq           SHARE-seq R1/R2/R3 barcode correction + poly-G QC
  trim-fastq              dovetail (R1/R2 overlap) trimming
  tss-enrichment          bulk TSS score + profile plot and per-barcode ArchR TSSe
                          (tabix via pysam, or gzip/plain fragment files)
  snapatac2-tsse          snapatac2 import + TSSe + FRiP (optional dependency)
  rna-qc-metrics          per-barcode total_counts / genes from a kb h5ad
  mtx-to-h5ad             mtx + barcodes + genes -> h5ad
  modify-barcode-h5       append _<subpool> to /obs/barcode in place
  joint-qc                joint RNA x ATAC barcode QC table + scatter + density plots
  barcode-rank            elbow / knee of any barcode-rank curve
  atac-qc-plots           fragment barcode-rank plots with elbow/knee
  rna-qc-plots            UMI / gene barcode-rank plots and genes-vs-UMIs
  insert-size-hist        Picard CollectInsertSizeMetrics histogram -> PNG + TSV
  tenx-barcode-map        10x multiome ATAC <-> RNA onlist conversion dictionary
  log-atac                chromap log + barcode summary -> qc_metrics.json
  kb-count / kb-index     kallisto|bustools wrappers (optional binary)
  chromap-align / chromap-index   chromap wrappers (optional binary)
  subpool-fragments       append _<subpool> to fragment / barcode-summary barcodes
  genome-tsv              parse the genome TSV (fasta, kb_nac_idx_tar, chromap_idx_tar)
  check-inputs            classify / fetch inputs (gs:// kept, syn*, https portal)
  sample-fastqs           first N reads of FASTQs
  portal-download         IGVF portal download with IGVF_API_KEY / IGVF_SECRET_KEY
  synapse-manifest        Synapse upload manifest from a results table (dry run)
  synapse-upload          batch upload grouped by parent (dry run unless --execute)
  synapse-annotations     file_type / status / data_type / description annotations
  html-report             tabbed HTML + CSV summary of plots, metrics and logs
  pipeline                the single_cell_pipeline.wdl plan (inputs, mapping, kb,
                          chromap, log) + downstream QC when outputs exist
  selftest                synthetic reads / fragments / matrices with planted signal

Output: Docs/IGVFscPipeline/<timestamp>_<label>/ (or --out-dir).

Usage:
    igvfagent igvf-sc-pipeline barcode-revcomp-detect --fastq bc.fq.gz --onlist 737K-arc-v1.txt --offset 8
    igvfagent igvf-sc-pipeline correct-fastq --read1 R1.fq.gz --read2 R2.fq.gz --whitelist bc24.txt --sample-type ATAC --prefix lib1
    igvfagent igvf-sc-pipeline trim-fastq --read1 R1.fq.gz --read2 R2.fq.gz
    igvfagent igvf-sc-pipeline tss-enrichment --fragments atac.fragments.tsv.gz --regions tss.bed --prefix lib1
    igvfagent igvf-sc-pipeline rna-qc-metrics --h5ad rna.h5ad --kb-workflow nac --subpool SP1
    igvfagent igvf-sc-pipeline joint-qc --rna-metrics rna.tsv --atac-metrics atac.tsv --pkr SP1
    igvfagent igvf-sc-pipeline rna-qc-plots --metrics rna.tsv
    igvfagent igvf-sc-pipeline kb-count --read1 L1_R1.fq.gz --read2 L1_R2.fq.gz --index-dir idx --read-format 1,0,16:1,16,28:0,0,0 --onlist 737K.txt
    igvfagent igvf-sc-pipeline pipeline --inputs-json inputs.json
    igvfagent igvf-sc-pipeline selftest --no-plots
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import io
import json
import logging
import math
import os
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "IGVFscPipeline"

UPSTREAM_REPO = "IGVF/single-cell-pipeline"
UPSTREAM_COMMIT = "202a21a9a57fcd94a82cc36c7e657e1885751ed7"
ATOMIC_REPO = "IGVF/atomic-workflows"
ATOMIC_TAG = "v1.1"
SCRIPT_BASE = "https://github.com/IGVF/single-cell-pipeline/blob/main/tasks/"

BLUE, ORANGE, GREEN, GREY, INK, INK2, AXIS, SURFACE = ("#2a78d6", "#eb6834", "#1baf7a", "#96a0b3",
                                                       "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb")
QC_COLORS = {"both": BLUE, "RNA only": ORANGE, "ATAC only": GREEN, "neither": GREY}

log = logging.getLogger("igvf_sc_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"igvf_sc_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir(args: argparse.Namespace, default_label: str) -> Path:
    if getattr(args, "out_dir", None):
        d = Path(args.out_dir)
        d.mkdir(parents=True, exist_ok=True)
        return d
    return run_dir(getattr(args, "label", None) or default_label)


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "List[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        o = float(o)
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, Path):
        return str(o)
    return o


def finish(d: Path, title: str, summary: dict, body: str = "") -> int:
    """Write report.md + summary.json into the run dir and print the harvest lines."""
    summary = dict(summary)
    summary.setdefault("upstream", {"repo": UPSTREAM_REPO, "commit": UPSTREAM_COMMIT})
    (d / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2))
    lines = [f"# {title}", "", f"Port of https://github.com/{UPSTREAM_REPO} (commit {UPSTREAM_COMMIT[:7]}).", ""]
    if body:
        lines += [body, ""]
    flat = [(k, v) for k, v in summary.items() if not isinstance(v, (dict, list)) and k != "upstream"]
    if flat:
        lines += ["## Summary", "", md_table(["key", "value"], flat), ""]
    (d / "report.md").write_text("\n".join(lines))
    print(f"Report: {d / 'report.md'}")
    print(f"JSON: {d / 'summary.json'}")
    return 0


def write_tsv(df, path: Path, index: bool = False, **kw) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, **kw)
    print(f"TSV: {path}")
    return path


# ---------------------------------------------------------------------------
# Sequence / file helpers
# ---------------------------------------------------------------------------

REV_COMP = str.maketrans("ATGC", "TACG")
REV_COMP_FULL = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Upstream reverse complement (ATGC only; N and lower case kept)."""
    return seq.translate(REV_COMP)[::-1]


def is_gzipped(path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def open_text(path, mode: str = "rt"):
    """Read: sniff gzip magic.  Write: gzip when the name ends with .gz."""
    path = str(path)
    if "r" in mode:
        return gzip.open(path, "rt") if is_gzipped(path) else open(path, "r")
    if path.endswith(".gz"):
        return gzip.open(path, "wt", compresslevel=3)
    return open(path, "w")


def read_fastq(path) -> "Iterator[Tuple[str, str, str]]":
    """Yield (name line, sequence, quality) records, stripped."""
    with open_text(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            yield name.strip(), seq, qual


def read_lines(path) -> "List[str]":
    with open_text(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def which(binary: str) -> Optional[str]:
    return shutil.which(binary)


def run_or_print(cmd: str, binaries: "List[str]", cwd: Optional[Path] = None) -> "Tuple[bool, str]":
    """Run `cmd` through bash -o pipefail if every binary is on PATH; else print it.  Returns (ran, stdout)."""
    missing = [b for b in binaries if not which(b)]
    print(f"Command: {cmd}")
    if missing:
        print(f"Note: {', '.join(missing)} not on PATH -- command not run (this is the exact command the pipeline would execute).")
        return False, ""
    log.info("run: %s", cmd)
    p = subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd], cwd=str(cwd) if cwd else None,
                       capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f"command failed (rc={p.returncode}): {cmd}\n{p.stderr.strip()[-2000:]}")
    return True, p.stdout


# ---------------------------------------------------------------------------
# barcode_revcomp_detect.py
# ---------------------------------------------------------------------------

def barcode_revcomp_detect(fastq, onlist, offset: int, num_reads: int = 100_000, thresh: float = 0.45) -> dict:
    bc = set(read_lines(onlist))
    bcrc = {reverse_complement(b) for b in bc}
    direct = revc = seen = 0
    for i, (_, seq, _) in enumerate(read_fastq(fastq)):
        if i >= num_reads:
            break
        seen += 1
        s = seq[offset:]
        direct += s in bc
        revc += s in bcrc
    direct_prop = direct / num_reads  # upstream divides by the requested number, not the reads seen
    rc_prop = revc / num_reads
    valid = direct_prop >= thresh or rc_prop >= thresh
    forward = direct_prop >= rc_prop
    fmt = f"bc:{offset}:-1" if forward else f"bc:{offset}:-1:-"
    return {"reads_seen": seen, "direct_match_proportion": direct_prop, "revcomp_match_proportion": rc_prop,
            "revcomp_chosen": not forward, "valid": valid, "read_format": fmt if valid else None}


def cmd_barcode_revcomp_detect(args) -> int:
    d = out_dir(args, "barcode_revcomp")
    r = barcode_revcomp_detect(args.fastq, args.onlist, args.offset, args.num_reads, args.threshold)
    qc = d / "barcode_revcomp_qc.txt"
    qc.write_text(f"Direct match proportion: {r['direct_match_proportion']}\n"
                  f"Reverse-complement match proportion: {r['revcomp_match_proportion']}\n"
                  f"Reverse-complement chosen: {r['revcomp_chosen']}\n")
    print(f"Wrote: {qc}")
    if r["valid"]:
        fmt = d / "barcode_read_format.txt"
        fmt.write_text(r["read_format"] + "\n")
        print(f"Wrote: {fmt}")
        print(f"read_format: {r['read_format']}")
    finish(d, "Barcode orientation detection", r,
           "Barcode read matched against the onlist and its reverse complement (upstream barcode_revcomp_detect.py).")
    if not r["valid"]:
        print(f"Error: insufficient barcode match rate: {r['direct_match_proportion']}, {r['revcomp_match_proportion']}")
        return 1
    return 0


# ---------------------------------------------------------------------------
# correct_fastq.py  (SHARE-seq)
# ---------------------------------------------------------------------------

def get_round_barcodes(whitelist, upstream_compat: bool = False) -> "Tuple[List[str], List[str], List[str]]":
    r1, r2, r3 = set(), set(), set()
    if upstream_compat:
        r1 = r2 = r3 = set()  # upstream: one shared set for all three rounds
    for line in read_lines(whitelist):
        r1.add(line[:8])
        r2.add(line[8:16])
        r3.add(line[16:24])
    return sorted(r1), sorted(r2), sorted(r3)


def create_barcode_dict(barcodes: "Iterable[str]", exact_priority: bool = True) -> "Dict[str, str]":
    """Every barcode plus its 1-substitution (ACGTN) neighbours -> barcode."""
    barcodes = list(barcodes)
    d: Dict[str, str] = {}
    for bc in barcodes:
        d[bc] = bc
        for i, base in enumerate(bc):
            for x in "ACGTN":
                if base != x:
                    d[bc[:i] + x + bc[i + 1:]] = bc
    if exact_priority:
        for bc in barcodes:
            d[bc] = bc
    return d


def check_putative_barcode(window: str, bdict: "Dict[str, str]", qual: str) -> "Tuple[Optional[str], str]":
    """Exact location [1:9], then 1-bp left shift [0:8], then right shift [2:] (padded with N/F if 9 bp)."""
    v = bdict.get(window[1:9])
    q = qual[1:9]
    if v is None:
        v = bdict.get(window[:8])
        q = qual[:8]
        if v is None:
            if len(window) < 10:
                v = bdict.get(window[2:] + "N")
                q = qual[2:] + "F"
            else:
                v = bdict.get(window[2:])
                q = qual[2:]
    return v, q


def correct_fastq(read1, read2, out_r1, out_r2, out_bc, whitelist, sample_type: str, prefix: str,
                  qc_path, pkr: Optional[str] = None, upstream_compat: bool = False) -> dict:
    if sample_type not in ("ATAC", "RNA"):
        raise SystemExit("sample_type must be ATAC or RNA")
    b1, b2, b3 = get_round_barcodes(whitelist, upstream_compat)
    d1, d2, d3 = (create_barcode_dict(b, exact_priority=not upstream_compat) for b in (b1, b2, b3))
    n = Counter()
    w1, w2, w3 = open_text(out_r1, "w"), open_text(out_r2, "w"), open_text(out_bc, "w")
    try:
        for (name1, seq1, q1s), (_, seq2, q2s) in zip(read_fastq(read1), read_fastq(read2)):
            bseq, bq = seq2[-99:], q2s[-99:]
            s1, s2, s3 = bseq[14:24], bseq[52:62], bseq[90:99]
            r1, q1 = check_putative_barcode(s1, d1, bq[14:24])
            r2, q2 = check_putative_barcode(s2, d2, bq[52:62])
            r3, q3 = check_putative_barcode(s3, d3, bq[90:99])
            if seq2[:10] == "G" * 10:
                n["poly_G_in_first_10bp"] += 1
            elif r1 and r2 and r3:
                n["match"] += 1
                tag = ",".join(x for x in (r1, r2, r3, pkr) if x)
                if sample_type == "RNA":
                    head = name1.split(" ")[0] + "_" + tag + "_" + seq2[:10]
                    w1.write(f"{head}\n{seq1}\n+\n{q1s}\n")
                    w2.write(f"{head}\n{r1 + r2 + r3 + seq2[:10]}\n+\n{q1 + q2 + q3 + q2s[:10]}\n")
                else:
                    head = name1.split(" ")[0] + "_" + tag
                    w1.write(f"{head}\n{seq1}\n+\n{q1s}\n")
                    w2.write(f"{head}\n{seq2[:-99]}\n+\n{q2s[:-99]}\n")
                    w3.write(f"{head}\n{r1 + r2 + r3}\n+\n{q1 + q2 + q3}\n")
            elif "G" * 8 in s1 and "G" * 8 in s2 and "G" * 8 in s3:
                n["poly_G_barcode"] += 1
            else:
                n["mismatch"] += 1
    finally:
        for w in (w1, w2, w3):
            w.close()
    fields = ["library", "match", "mismatch", "poly_G_barcode", "poly_G_in_first_10bp"]
    vals = [prefix] + [n[f] for f in fields[1:]]
    Path(qc_path).write_text("\t".join(fields) + "\n" + "\t".join(str(v) for v in vals))
    return dict(zip(fields, vals))


def cmd_correct_fastq(args) -> int:
    d = out_dir(args, f"correct_fastq_{args.prefix}")
    ext = ".fastq.gz"
    o1, o2, ob = d / f"{args.prefix}.R1.corrected{ext}", d / f"{args.prefix}.R2.corrected{ext}", d / f"{args.prefix}.barcode.corrected{ext}"
    qc = d / f"{args.prefix}_barcode_qc.txt"
    r = correct_fastq(args.read1, args.read2, o1, o2, ob, args.whitelist, args.sample_type, args.prefix, qc,
                      args.pkr, args.upstream_compat)
    for p in (o1, o2, ob, qc):
        print(f"Wrote: {p}")
    total = sum(v for k, v in r.items() if k != "library")
    r["total_reads"] = total
    r["match_fraction"] = (r["match"] / total) if total else float("nan")
    return finish(d, f"SHARE-seq barcode correction ({args.prefix})", r,
                  "Round barcodes corrected with exact / 1-mismatch / +-1 bp shift lookups (upstream correct_fastq.py).")


# ---------------------------------------------------------------------------
# trim_fastq.py
# ---------------------------------------------------------------------------

def levenshtein_le1(a: str, b: str) -> int:
    """Levenshtein distance capped at 2 (i.e. the score_cutoff=1 behaviour: >1 reported as 2)."""
    try:
        import Levenshtein  # type: ignore
        return Levenshtein.distance(a, b, score_cutoff=1)
    except Exception:
        pass
    la, lb = len(a), len(b)
    if abs(la - lb) > 1:
        return 2
    if la == lb:
        diff = 0
        for x, y in zip(a, b):
            if x != y:
                diff += 1
                if diff > 1:
                    return 2
        return diff
    if la > lb:
        a, b, la, lb = b, a, lb, la
    i = j = 0
    edits = 0
    while i < la and j < lb:
        if a[i] != b[j]:
            edits += 1
            if edits > 1:
                return 2
            j += 1
        else:
            i += 1
            j += 1
    edits += (lb - j) + (la - i)
    return edits if edits <= 1 else 2


def fuzz_align(query: str, seq: str) -> int:
    for i in range(len(seq)):  # first window from the LEFT (upstream comment says rfind, code scans left to right)
        if levenshtein_le1(seq[i:i + len(query)], query) <= 1:
            return i
    return -1


def dovetail_position(seq1: str, seq2: str) -> int:
    query = reverse_complement(seq2[0:20])
    idx = seq1.rfind(query)
    if idx == -1:
        idx = fuzz_align(query, seq1)
    return idx + 20 if idx > -1 else -1


def trim_fastq(read1, read2, out1, out2, stats_path) -> dict:
    total = trimmed = 0
    with open_text(out1, "w") as w1, open_text(out2, "w") as w2:
        for (n1, s1, q1), (n2, s2, q2) in zip(read_fastq(read1), read_fastq(read2)):
            total += 2
            where = dovetail_position(s1, s2)
            if where > -1:
                trimmed += 2
                s1, q1, s2, q2 = s1[:where], q1[:where], s2[:where], q2[:where]
            w1.write(f"{n1}\n{s1}\n+\n{q1}\n")
            w2.write(f"{n2}\n{s2}\n+\n{q2}\n")
    pct = trimmed / total * 100 if total > 0 else 0
    Path(stats_path).write_text("total_reads\tuntrimmed_reads\ttrimmed_reads\t%trimmed\n"
                                + "%i\t%i\t%i\t%0.1f" % (total, total - trimmed, trimmed, pct))
    return {"total_reads": total, "untrimmed_reads": total - trimmed, "trimmed_reads": trimmed, "pct_trimmed": round(pct, 1)}


def cmd_trim_fastq(args) -> int:
    d = out_dir(args, "trim_fastq")
    pre = args.prefix
    o1, o2, st = d / f"{pre}.R1.trimmed.fastq.gz", d / f"{pre}.R2.trimmed.fastq.gz", d / f"{pre}.trimming_stats.txt"
    r = trim_fastq(args.read1, args.read2, o1, o2, st)
    for p in (o1, o2, st):
        print(f"Wrote: {p}")
    return finish(d, "Dovetail trimming", r, "Read 1 / read 2 overlap removed (upstream trim_fastq.py).")


# ---------------------------------------------------------------------------
# compute_tss_enrichment_bulk.py  (bulk + per-barcode ArchR TSSe)
# ---------------------------------------------------------------------------

def read_tss_regions(path, strand_col: int = 4):
    """Columns chr, start, end, strand (strand column 1-based, default 4), like np.loadtxt(usecols=(0,1,2,s-1))."""
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    out = pd.DataFrame({"chr": df[0].astype(str), "start": df[1].astype(int), "end": df[2].astype(int),
                        "strand": df[strand_col - 1].astype(str)})
    return out


def read_fragments(path):
    """chrom, start, end, barcode from a (b)gzip or plain fragment file."""
    pd = _pd()
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3],
                     names=["chr", "start", "end", "barcode"], dtype={0: str, 3: str},
                     compression="gzip" if is_gzipped(path) else None)
    return df


def _tabix_ok(path) -> bool:
    try:
        import pysam  # type: ignore  # noqa: F401
    except Exception:
        return False
    return Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists()


def _insertions_pure(frags, tss, flank: int):
    """(barcode, tss_row, index-in-promoter) for every fragment end inside a promoter window.

    Reproduces the tabix path exactly: a fragment is seen by TSS t when it overlaps
    [max(0, t-flank), t+flank) half-open, and an end at p is counted when
    0 <= p-(t-flank) <= 2*flank-1.  Start ends: t in [p-flank+1, p+flank]; end
    positions: the overlap rule additionally excludes end == t-flank.
    """
    np = _np()
    out_bc, out_row, out_idx = [], [], []
    for chrom, tsub in tss.groupby("chr", sort=False):
        fsub = frags[frags["chr"] == chrom]
        if fsub.empty:
            continue
        order = np.argsort(tsub["start"].to_numpy(), kind="stable")
        t_sorted = tsub["start"].to_numpy()[order]
        rows = tsub.index.to_numpy()[order]
        bcs = fsub["barcode"].to_numpy()
        for pos, hi_off in ((fsub["start"].to_numpy(), flank), (fsub["end"].to_numpy(), flank - 1)):
            lo = np.searchsorted(t_sorted, pos - flank + 1, side="left")
            hi = np.searchsorted(t_sorted, pos + hi_off, side="right")
            cnt = hi - lo
            tot = int(cnt.sum())
            if tot == 0:
                continue
            rep = np.repeat(np.arange(len(pos)), cnt)
            starts = np.repeat(np.cumsum(cnt) - cnt, cnt)
            k = lo[rep] + (np.arange(tot) - starts)
            out_bc.append(bcs[rep])
            out_row.append(rows[k])
            out_idx.append(pos[rep] - (t_sorted[k] - flank))
    if not out_bc:
        return np.array([], dtype=object), np.array([], dtype=int), np.array([], dtype=int)
    return np.concatenate(out_bc), np.concatenate(out_row), np.concatenate(out_idx)


def _insertions_tabix(path, tss, flank: int):
    """Same triples through pysam.TabixFile.fetch, TSS by TSS (the upstream access pattern)."""
    import pysam  # type: ignore
    np = _np()
    tb = pysam.TabixFile(str(path))
    contigs = set(tb.contigs)
    bcs, rows, idxs = [], [], []
    for row, chrom, t in zip(tss.index, tss["chr"], tss["start"]):
        ps, pe = int(t) - flank, int(t) + flank
        if chrom not in contigs:
            continue
        for rec in tb.fetch(chrom, max(0, ps), pe):
            f = rec.split("\t")
            for p in (int(f[1]), int(f[2])):
                if ps <= p <= pe - 1:
                    bcs.append(f[3])
                    rows.append(row)
                    idxs.append(p - ps)
    return np.array(bcs, dtype=object), np.array(rows, dtype=int), np.array(idxs, dtype=int)


def tss_enrichment(fragments, regions, flank: int = 2000, window: int = 20, strand_col: int = 4,
                   engine: str = "auto") -> dict:
    """Bulk TSS profile / score and per-barcode ArchR TSS enrichment."""
    np, pd = _np(), _pd()
    if flank < 151:
        raise SystemExit("--flank must be >= 151 (the ArchR bins need 100 + 101 + 100 bp)")
    tss = read_tss_regions(regions, strand_col).reset_index(drop=True)
    use_tabix = engine == "tabix" or (engine == "auto" and _tabix_ok(fragments))
    if use_tabix:
        bc, row, idx = _insertions_tabix(fragments, tss, flank)
        frags = None
    else:
        frags = read_fragments(fragments)
        bc, row, idx = _insertions_pure(frags, tss, flank)
    size = 2 * flank
    max_window = size - 1
    tss_pos = int(max_window / 2)
    minus = (tss["strand"].to_numpy()[row] == "-") if len(row) else np.zeros(0, bool)
    # bulk profile (index mirrored for '-' TSSs); upstream merge starts from np.ones -> pseudocount 1
    bulk_idx = np.where(minus, max_window - idx, idx)
    bulk = np.ones(size)
    np.add.at(bulk, bulk_idx, 1)
    norm = float(np.mean(bulk[np.r_[0:100, size - 100 - 1:size]]))
    raw_signal = bulk / max(0.2, norm)
    smoothed = np.convolve(bulk, np.ones(window), "same") / window / norm
    score = float(np.max(smoothed))
    # per barcode: 301 ArchR bins
    red = np.full(len(idx), -1)
    left = idx <= 99
    right = (~left) & (idx >= max_window - 99)
    mid = (~left) & (~right) & (idx >= tss_pos - 50) & (idx <= tss_pos + 50)
    red[left] = idx[left]
    red[right] = idx[right] - (max_window - 300)
    red[mid] = idx[mid] - (tss_pos - 50) + 100
    has = red >= 0
    red = np.where(has & minus, 300 - red, red)
    per = pd.DataFrame({"barcode": bc, "bin": red, "has": has})
    promoter_reads = per.groupby("barcode").size() if len(per) else pd.Series(dtype=int)
    binned = per[per["has"]]
    if len(binned):
        tab = binned.groupby(["barcode", "bin"]).size().unstack(fill_value=0).reindex(columns=range(301), fill_value=0)
    else:
        tab = pd.DataFrame(columns=range(301), dtype=int)
    M = tab.to_numpy(dtype=float)
    reads_tss = M[:, 100:201].sum(axis=1)
    flank_norm = 2 * (M[:, np.r_[0:100, 201:301]].sum(axis=1)) / 200
    tsse = 2 * reads_tss / 101 / np.maximum(0.2, flank_norm)
    stats = pd.DataFrame({"barcode": tab.index.astype(str), "reads_tss": reads_tss.astype(int),
                          "reads_tss_flank": M[:, np.r_[0:100, 201:301]].sum(axis=1).astype(int), "tss_enrichment": tsse})
    if frags is None:
        frags = read_fragments(fragments)
    n_frag = frags.groupby("barcode").size()
    allbc = pd.DataFrame({"barcode": n_frag.index.astype(str), "reads_unique": (2 * n_frag.to_numpy()).astype(int)})
    allbc["reads_promoter"] = allbc["barcode"].map(promoter_reads).fillna(0).astype(int)
    allbc = allbc.merge(stats, on="barcode", how="left")
    allbc["reads_tss"] = allbc["reads_tss"].fillna(0).astype(int)
    allbc["reads_tss_flank"] = allbc["reads_tss_flank"].fillna(0).astype(int)
    allbc["tss_enrichment"] = allbc["tss_enrichment"].fillna(0.0)
    allbc = allbc[["barcode", "reads_unique", "reads_promoter", "reads_tss", "tss_enrichment", "reads_tss_flank"]]
    allbc = allbc.sort_values("reads_unique", ascending=False, kind="stable").reset_index(drop=True)
    return {"engine": "tabix" if use_tabix else "pandas", "bulk_profile": bulk, "raw_signal": raw_signal,
            "smoothed": smoothed, "tss_score_bulk": score, "normalization_factor": norm, "per_barcode": allbc,
            "n_tss": len(tss), "n_insertions_in_promoters": int(len(idx)), "flank": flank}


def cmd_tss_enrichment(args) -> int:
    d = out_dir(args, f"tss_{args.prefix}")
    r = tss_enrichment(args.fragments, args.regions, args.flank, args.window, args.strand_col, args.engine)
    pre = args.prefix
    (d / f"{pre}.tss_score_bulk.txt").write_text(f"{r['tss_score_bulk']}\n")
    print(f"Wrote: {d / f'{pre}.tss_score_bulk.txt'}")
    write_tsv(r["per_barcode"], d / f"{pre}.tss_enrichment_barcode_stats.tsv")
    pd = _pd()
    write_tsv(pd.DataFrame({"position": range(-r["flank"], r["flank"]), "counts": r["bulk_profile"],
                            "raw_signal": r["raw_signal"], "smoothed": r["smoothed"]}), d / f"{pre}.tss_profile_bulk.tsv")
    plt = None if args.no_plots else _plt()
    if plt is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(r["raw_signal"], ".", color=INK, ms=2)
        ax.plot(r["smoothed"], color=ORANGE, lw=1.5)
        ax.set_xticks([0, r["flank"], 2 * r["flank"]])
        ax.set_xticklabels([f"-{r['flank']}", "TSS", f"+{r['flank']}"])
        _style(ax, f"Bulk TSS enrichment ({pre}): {r['tss_score_bulk']:.2f}", "", "TSS enrichment")
        _save(fig, d / f"{pre}.tss_enrichment_bulk.png")
    pb = r["per_barcode"]
    summ = {"prefix": pre, "engine": r["engine"], "tss_score_bulk": r["tss_score_bulk"], "n_tss": r["n_tss"],
            "n_barcodes": len(pb), "median_barcode_tss_enrichment": float(pb["tss_enrichment"].median()) if len(pb) else None,
            "insertions_in_promoters": r["n_insertions_in_promoters"], "flank": r["flank"], "window": args.window}
    return finish(d, f"TSS enrichment ({pre})", summ,
                  "Bulk TSS score (max of the smoothed, flank-normalised profile) and per-barcode ArchR TSS enrichment.")


# ---------------------------------------------------------------------------
# snapatac2-tss-enrichment.py (optional dependency)
# ---------------------------------------------------------------------------

def cmd_snapatac2_tsse(args) -> int:
    d = out_dir(args, f"snapatac2_{args.prefix}")
    try:
        import snapatac2 as snap  # type: ignore
    except Exception:
        print("Note: snapatac2 is not installed -- skipping (pip install snapatac2). "
              "Use `tss-enrichment` for the dependency-free ArchR-style TSSe.")
        return finish(d, "snapatac2 TSSe (skipped)", {"skipped": True, "reason": "snapatac2 not installed"})
    sizes = {}
    for ln in read_lines(args.chrom_sizes):
        f = ln.split("\t")
        sizes[f[0]] = int(f[1])
    data = snap.pp.import_data(fragment_file=args.fragments, chrom_sizes=sizes, sorted_by_barcode=False,
                               min_num_fragments=args.min_frag_cutoff, shift_left=0, shift_right=1)
    snap.metrics.tsse(data, args.gtf)
    if not args.no_plots:
        for cut in (args.min_frag_cutoff, 500):
            png = d / f"{args.prefix}.cutoff{cut}.tsse.png"
            snap.pl.tsse(data, min_fragment=cut, width=800, height=1000, show=False, out_file=str(png))
            print(f"Figure: {png}")
    snap.metrics.frip(data, {"tss_frac": args.tss_bed, "promoter_frac": args.promoter_bed})
    tsv = d / f"{args.prefix}.barcode_stats.tsv"
    data.obs.to_csv(tsv, index_label="barcode", sep="\t")
    print(f"TSV: {tsv}")
    data.obs["sample"] = args.prefix
    h5 = d / f"{args.prefix}.snapatac.h5ad"
    data.write(str(h5))
    print(f"Wrote: {h5}")
    return finish(d, "snapatac2 TSSe", {"n_barcodes": int(data.n_obs), "min_frag_cutoff": args.min_frag_cutoff})


# ---------------------------------------------------------------------------
# qc_rna_extract_metrics.py / write_h5ad_from_mtx.py / modify_barcode_h5.py
# ---------------------------------------------------------------------------

def _qc_counts(M):
    """scanpy calculate_qc_metrics -> (total_counts, n_genes_by_counts) per row."""
    np = _np()
    try:
        import scipy.sparse as sp  # type: ignore
        if sp.issparse(M):
            M = M.tocsr()
            return np.asarray(M.sum(axis=1)).ravel(), np.asarray((M > 0).sum(axis=1)).ravel().astype(int)
    except ImportError:
        pass
    M = np.asarray(M)
    return M.sum(axis=1), (M > 0).sum(axis=1).astype(int)


def rna_qc_metrics(adata, kb_workflow: str = "standard", subpool: Optional[str] = "none"):
    pd = _pd()
    if kb_workflow == "standard":
        tot, ng = _qc_counts(adata.X)
    else:
        tot = ng = 0
        for layer in ("mature", "ambiguous", "nascent"):
            if layer not in adata.layers:
                raise SystemExit(f"kb nac h5ad is missing layer '{layer}' (have: {list(adata.layers.keys())})")
            t, g = _qc_counts(adata.layers[layer])
            tot, ng = tot + t, ng + g
    df = pd.DataFrame({"total_counts": tot, "genes": ng}, index=adata.obs_names.astype(str))
    df.index.name = adata.obs_names.name
    if subpool not in (None, "", "none"):
        df.index = df.index.astype(str) + "_" + subpool
        df.index.name = adata.obs_names.name
    return df


def cmd_rna_qc_metrics(args) -> int:
    import anndata  # type: ignore
    d = out_dir(args, "rna_qc_metrics")
    ad = anndata.read_h5ad(args.h5ad)
    df = rna_qc_metrics(ad, args.kb_workflow, args.subpool)
    path = d / (args.out_name or "rna_barcode_metadata.tsv")
    df.to_csv(path, sep="\t")
    print(f"TSV: {path}")
    return finish(d, "RNA barcode metrics", {"kb_workflow": args.kb_workflow, "n_barcodes": len(df),
                                             "total_umis": float(df["total_counts"].sum()),
                                             "median_umis": float(df["total_counts"].median()) if len(df) else None,
                                             "subpool": args.subpool})


def write_h5ad_from_mtx(mtx, barcodes, genes, out):
    import anndata  # type: ignore
    import scipy.io  # type: ignore
    import scipy.sparse as sp  # type: ignore
    X = sp.csr_matrix(scipy.io.mmread(str(mtx)).astype("float32"))  # sc.read_mtx: float32, rows = obs
    ad = anndata.AnnData(X=X)
    ad.obs_names = read_lines(barcodes)
    ad.var_names = read_lines(genes)
    ad.write_h5ad(str(out))
    return ad


def cmd_mtx_to_h5ad(args) -> int:
    d = out_dir(args, "mtx_to_h5ad")
    out = Path(args.out) if args.out else d / "output.h5ad"
    ad = write_h5ad_from_mtx(args.mtx, args.barcodes, args.genes, out)
    print(f"Wrote: {out}")
    return finish(d, "mtx -> h5ad", {"n_obs": ad.n_obs, "n_vars": ad.n_vars, "h5ad": str(out)})


def append_suffix_to_h5ad(h5ad_file, suffix: str, dataset: str = "/obs/barcode") -> int:
    """Append "_<suffix>" to every entry of /obs/barcode (falls back to the obs index dataset)."""
    import h5py  # type: ignore
    np = _np()
    with h5py.File(str(h5ad_file), "r+") as h5:
        key = dataset
        if key not in h5:
            idx = h5["obs"].attrs.get("_index", "_index")
            idx = idx.decode() if isinstance(idx, bytes) else idx
            key = f"/obs/{idx}"
        ds = h5[key]
        vals = [(v.decode("utf-8") if isinstance(v, bytes) else str(v)) + "_" + suffix for v in ds[()]]
        if h5py.check_string_dtype(ds.dtype) is not None and h5py.check_string_dtype(ds.dtype).length is None:
            ds[...] = np.array(vals, dtype=object)
        else:  # fixed-length bytes would truncate -> recreate with the same attributes
            attrs = dict(ds.attrs)
            del h5[key]
            new = h5.create_dataset(key, data=np.array([v.encode() for v in vals], dtype="S"))
            for k, v in attrs.items():
                new.attrs[k] = v
        h5.flush()
    return len(vals)


def cmd_modify_barcode_h5(args) -> int:
    d = out_dir(args, "modify_barcode_h5")
    n = append_suffix_to_h5ad(args.h5ad, args.suffix)
    print(f"Wrote: {args.h5ad}")
    return finish(d, "Subpool suffix appended", {"h5ad": args.h5ad, "suffix": args.suffix, "n_barcodes": n})


# ---------------------------------------------------------------------------
# joint_cell_plotting.py + joint_cell_plotting_density.R
# ---------------------------------------------------------------------------

def _split_lines(path, skip_header: bool = True):
    with open_text(path) as fh:
        if skip_header:
            next(fh, None)
        for ln in fh:
            if ln.strip():
                yield ln.rstrip("\n").rstrip("\r").split("\t")


def joint_metrics(rna_metrics_file, atac_metrics_file, remove_low_yielding_cells: int = 10):
    """Upstream get_metrics: RNA col1 = UMIs, col2 = genes; ATAC col1 = reads (frags = reads/2), col4 = TSSe."""
    pd = _pd()
    rna = {}
    for f in _split_lines(rna_metrics_file):
        umi = int(float(f[1]))
        if umi >= remove_low_yielding_cells:
            rna[f[0]] = (umi, int(float(f[2])))
    atac = {}
    for f in _split_lines(atac_metrics_file):
        if int(float(f[1])) / 2 >= remove_low_yielding_cells:
            atac[f[0]] = (float(f[4]), int(float(f[1])) / 2)
    keys = sorted(set(rna) | set(atac))
    rows = {k: rna.get(k, (0, 0)) + atac.get(k, (0, 0)) for k in keys}
    return pd.DataFrame.from_dict(rows, orient="index", columns=["umis", "genes", "tss", "frags"])


def qc_cells(df, min_umis: int, min_genes: int, min_tss: float, min_frags: int):
    np = _np()
    pr = (df["umis"] >= min_umis) & (df["genes"] >= min_genes)
    pa = (df["tss"] >= min_tss) & (df["frags"] >= min_frags)
    df = df.copy()
    df["QC"] = np.select([pr & pa, pr, pa, (~pr) & (~pa)], ["both", "RNA only", "ATAC only", "neither"], default="neither")
    counts = df["QC"].value_counts()
    df["QC_count"] = [f"{q} ({counts[q]})" for q in df["QC"]]
    return df


def _round_pow10(x: float) -> float:
    return 10 ** math.ceil(math.log10(max(x, 1)))


def plot_joint(df, pkr, cut: dict, path: Path):
    plt = _plt()
    if plt is None or df.empty:
        return None
    lim = _round_pow10(max(df["frags"].max(), df["umis"].max()))
    fig, ax = plt.subplots(figsize=(8, 6))
    for q in ("neither", "ATAC only", "RNA only", "both"):
        s = df[df["QC"] == q]
        if len(s):
            ax.scatter(s["frags"].clip(lower=10), s["umis"].clip(lower=10), s=3, color=QC_COLORS[q], label=s["QC_count"].iloc[0], lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, lim)
    ax.set_ylim(10, lim)
    ax.legend(title="QC", fontsize=7, title_fontsize=8, frameon=False)
    _style(ax, f"Joint Cell Calling ({pkr})", "ATAC Unique Fragments per Barcode", "RNA UMIs per Barcode")
    fig.text(0.5, -0.02, f"ATAC cutoffs: TSS >= {cut['min_tss']}, frags >= {cut['min_frags']}. "
                         f"RNA cutoffs: UMIs >= {cut['min_umis']}, genes >= {cut['min_genes']}", ha="center", fontsize=8, color=INK2)
    return _save(fig, path)


def plot_joint_density(df, pkr, path: Path, radius: float = 0.05):
    """ggpointdensity analogue: colour = neighbours within `radius` log10 units; only if any 'both'."""
    plt = _plt()
    if plt is None or not (df["QC"] == "both").any():
        return None
    np = _np()
    from scipy.spatial import cKDTree  # type: ignore
    s = df[df["QC"].isin(["RNA only", "ATAC only", "both"])]
    pts = np.log10(np.c_[s["frags"].clip(lower=1), s["umis"].clip(lower=1)])
    dens = np.array([len(v) for v in cKDTree(pts).query_ball_point(pts, r=radius)])
    lim = _round_pow10(max(s["frags"].max(), s["umis"].max()))
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
    samba = LinearSegmentedColormap.from_list("sambaNight", ['#1873CC', '#1798E5', '#00BFFF', '#4AC596', '#00CC00',
                                                             '#A2E700', '#FFFF00', '#FFD200', '#FFA500'])
    fig, ax = plt.subplots(figsize=(8.75, 6))
    o = np.argsort(dens)
    sc = ax.scatter(s["frags"].to_numpy()[o], s["umis"].to_numpy()[o], c=dens[o], cmap=samba, s=4, lw=0)
    fig.colorbar(sc, ax=ax, label="n_neighbors")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, lim)
    ax.set_ylim(10, lim)
    _style(ax, f"Joint Cell Calling ({pkr}): Density Plot", "ATAC Unique Fragments per Barcode", "RNA UMIs per Barcode")
    return _save(fig, path)


def cmd_joint_qc(args) -> int:
    d = out_dir(args, f"joint_qc_{args.pkr or 'sample'}")
    df = joint_metrics(args.rna_metrics, args.atac_metrics, args.remove_low_yielding_cells)
    cut = {"min_umis": args.min_umis, "min_genes": args.min_genes, "min_tss": args.min_tss, "min_frags": args.min_frags}
    df = qc_cells(df, **cut)
    csvp = d / "joint_barcode_metadata.csv"
    df.to_csv(csvp)
    print(f"CSV: {csvp}")
    if not args.no_plots:
        plot_joint(df, args.pkr, cut, d / "joint_cell_calling.png")
        plot_joint_density(df, args.pkr, d / "joint_cell_calling_density.png")
    counts = df["QC"].value_counts().to_dict()
    summ = dict(cut, remove_low_yielding_cells=args.remove_low_yielding_cells, n_barcodes=len(df),
                n_both=int(counts.get("both", 0)), n_rna_only=int(counts.get("RNA only", 0)),
                n_atac_only=int(counts.get("ATAC only", 0)), n_neither=int(counts.get("neither", 0)))
    return finish(d, f"Joint cell calling ({args.pkr})", summ,
                  md_table(["QC", "barcodes"], sorted(counts.items(), key=lambda kv: -kv[1])))


# ---------------------------------------------------------------------------
# barcode_rank_functions.R  (elbow / knee)
# ---------------------------------------------------------------------------

def _r_nknots(n: int) -> int:
    """R .nknots.smspl(n)."""
    if n < 50:
        return n
    a1, a2, a3, a4 = math.log2(50), math.log2(100), math.log2(140), math.log2(200)
    if n < 200:
        v = 2 ** (a1 + (a2 - a1) * (n - 50) / 150)
    elif n < 800:
        v = 2 ** (a2 + (a3 - a2) * (n - 200) / 600)
    elif n < 3200:
        v = 2 ** (a3 + (a4 - a3) * (n - 800) / 2400)
    else:
        v = 200 + (n - 3200) ** 0.2
    return int(math.trunc(v + 1e-9))  # guard 2^log2(100) = 99.999...


def smooth_spline_deriv2(x, y, spar: float = 1.0):
    """Second derivative at x of R smooth.spline(x, y, spar=spar) (unique x, unit weights).

    Same model as R: x scaled to [0,1]; cubic B-spline basis on knots at
    .nknots.smspl(n) of the sorted x (boundary knots x3); penalty Sigma_ij =
    int B_i'' B_j''; lambda = r * 256^(3*spar - 1) with r = sum(diag(X'X)[3:(nk-3)])
    / sum(diag(Sigma)[3:(nk-3)]) (sbart.c); coefficients from (X'X + lambda Sigma) c = X'y.
    """
    np = _np()
    from scipy.interpolate import BSpline  # type: ignore
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    o = np.argsort(x, kind="stable")
    xs, ys = x[o], y[o]
    rng = xs[-1] - xs[0]
    xb = (xs - xs[0]) / rng
    n = len(xb)
    nkn = _r_nknots(n)
    sel = np.floor(np.linspace(1, n, nkn)).astype(int) - 1
    inner = xb[sel]
    t = np.r_[[inner[0]] * 3, inner, [inner[-1]] * 3]
    nk = len(t) - 4
    X = BSpline.design_matrix(np.clip(xb, t[3], t[nk]), t, 3).toarray()
    # penalty via 2-point Gauss-Legendre per knot interval (B'' is linear -> exact)
    gx, gw = np.polynomial.legendre.leggauss(2)
    pts, wts = [], []
    for a, b in zip(t[3:nk], t[4:nk + 1]):
        if b > a:
            pts.extend((b - a) / 2 * gx + (a + b) / 2)
            wts.extend((b - a) / 2 * gw)
    pts, wts = np.array(pts), np.array(wts)
    D2 = BSpline(t, np.eye(nk), 3).derivative(2)(pts)
    Sig = D2.T @ (D2 * wts[:, None])
    XtX = X.T @ X
    ratio = np.trace(XtX[2:nk - 3, 2:nk - 3]) / np.trace(Sig[2:nk - 3, 2:nk - 3])
    lam = ratio * 256.0 ** (3 * spar - 1)
    coef = np.linalg.solve(XtX + lam * Sig, X.T @ ys)
    d2 = BSpline(t, coef, 3).derivative(2)(xb) / rng ** 2
    out = np.empty_like(d2)
    out[o] = d2
    return out


def get_vectors(x, y):
    """R get_vectors: slice to the knee region if the 2nd derivative changes sign in the middle 80%."""
    np = _np()
    x, y = np.asarray(x), np.asarray(y)
    sd = smooth_spline_deriv2(x, y, spar=1.0)
    n = len(sd)
    ten = int(round(n * 0.1))
    mid = sd[ten:n - ten]
    if np.all(mid >= 0) or np.all(mid <= 0):
        return x, y, "original"
    k = int(np.argmin(sd)) + 1  # R: second_deriv$x[which.min(...)] -- x are ranks 1..n
    left = np.where(sd[:k] >= 0)[0]
    right = np.where(sd[k - 1:] >= 0)[0]
    if len(left) == 0 or len(right) == 0:
        return x, y, "original"
    e1 = int(left[-1]) + 1
    e2 = k + int(right[0])
    return x[e1 - 1:e2], y[e1 - 1:e2], "sliced"


def elbow_knee_finder(x, y, mode: str = "basic"):
    np = _np()
    x, y = np.asarray(x, float), np.asarray(y, float)
    if mode == "advanced":
        if len(np.unique(x)) < 4:
            return None
        x, y, _ = get_vectors(x, y)
    if len(x) == 0 or len(y) == 0 or len(x) < 2 or x[-1] == x[0]:
        return None
    b = (y[-1] - y[0]) / (x[-1] - x[0])
    a = y[0] - b * x[0]
    dist = np.abs(b * x - y + a) / math.sqrt(b * b + 1)
    i = int(np.argmax(dist))
    return float(x[i]), float(y[i])


def get_elbow_knee_points(x, y) -> "List[float]":
    """[elbow_x, elbow_y(log10), knee_x, knee_y] (length 0, 2 or 4) as in barcode_rank_functions.R."""
    p1 = elbow_knee_finder(x, y, "basic")
    if p1 is None:
        return []
    k = int(p1[0])
    p2 = elbow_knee_finder(list(x)[:k], list(y)[:k], "advanced")
    return list(p1) + (list(p2) if p2 else [])


def barcode_rank(values, cutoff: float):
    """Filter >= cutoff, sort decreasing, ranks 1..n, elbow/knee on log10."""
    np = _np()
    v = np.asarray(values, float)
    v = np.sort(v[v >= cutoff])[::-1]
    rank = np.arange(1, len(v) + 1)
    pts = get_elbow_knee_points(rank, np.log10(v)) if len(v) else []
    return rank, v, pts


def _rank_panels(rank, vals, pts, title_all: str, title_top: str, ylab: str, ylim: float, path: Path):
    plt = _plt()
    if plt is None or not pts:
        return None
    np = _np()
    fig, axes = plt.subplots(2, 1, figsize=(8, 8))
    e = int(pts[0])
    ax = axes[0]
    ax.scatter(rank, vals, s=4, lw=0, color=np.where(rank <= e, "darkblue", "dimgrey"))
    ax.set_yscale("log")
    ax.set_ylim(1, ylim)
    ax.axvline(pts[0], color=INK, lw=0.8)
    ax.axhline(10 ** pts[1], color=INK, lw=0.8)
    ax.annotate(f"({e}, {10 ** pts[1]:.0f})", (pts[0], 10 ** pts[1]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    _style(ax, title_all, f"Barcode rank ({len(rank) - e} low quality cells)", ylab)
    ax = axes[1]
    if len(pts) > 2:
        tr, tv = rank[:e], vals[:e]
        ax.scatter(tr, tv, s=4, lw=0, color=np.where(tr <= pts[2], "darkblue", "dimgrey"))
        ax.set_yscale("log")
        ax.set_ylim(1, ylim)
        ax.axvline(pts[2], color=INK, lw=0.8)
        ax.axhline(10 ** pts[3], color=INK, lw=0.8)
        ax.annotate(f"({int(pts[2])}, {10 ** pts[3]:.0f})", (pts[2], 10 ** pts[3]), xytext=(5, 5), textcoords="offset points", fontsize=8)
        _style(ax, title_top, "Barcode rank", ylab)
    else:
        ax.axis("off")
    fig.tight_layout()
    return _save(fig, path)


def _points_dict(pts) -> dict:
    out = {}
    if len(pts) >= 2:
        out.update(elbow_rank=int(pts[0]), elbow_value=10 ** pts[1])
    if len(pts) >= 4:
        out.update(knee_rank=int(pts[2]), knee_value=10 ** pts[3])
    return out


def _read_metadata(path):
    pd = _pd()
    return pd.read_csv(path, sep=None, engine="python")


def cmd_barcode_rank(args) -> int:
    d = out_dir(args, "barcode_rank")
    df = _read_metadata(args.metrics)
    col = args.column if args.column in df.columns else df.columns[int(args.column) if str(args.column).isdigit() else 1]
    rank, v, pts = barcode_rank(df[col].to_numpy(), args.cutoff)
    res = _points_dict(pts)
    pd = _pd()
    write_tsv(pd.DataFrame({"rank": rank, col: v}), d / "barcode_rank.tsv")
    if not args.no_plots:
        _rank_panels(rank, v, pts, f"{col} per Barcode", f"{col} per Top-Ranked Barcode", f"{col} (log10 scale)",
                     max(10 ** math.ceil(math.log10(max(v.max(), 10))), 10) if len(v) else 10, d / "barcode_rank.png")
    return finish(d, f"Barcode rank ({col})", dict(res, column=col, cutoff=args.cutoff, n_barcodes=int(len(v))))


def cmd_atac_qc_plots(args) -> int:
    d = out_dir(args, "atac_qc_plots")
    df = _read_metadata(args.metrics)
    col = next((c for c in ("unique", "reads_unique", "fragments", "frags") if c in df.columns), None)
    if args.column:
        col = args.column
    if col is None:
        raise SystemExit(f"no fragment column (unique / reads_unique / fragments) in {list(df.columns)}")
    rank, v, pts = barcode_rank(df[col].to_numpy(), args.fragment_cutoff)
    if not args.no_plots:
        _rank_panels(rank, v, pts, "ATAC Fragments per Barcode", "ATAC Fragments per Top-Ranked Barcode",
                     "Fragments per barcode (log10 scale)", 100000, d / "atac_qc_barcode_rank.png")
    return finish(d, "ATAC barcode-rank QC", dict(_points_dict(pts), column=col, fragment_cutoff=args.fragment_cutoff,
                                                   n_barcodes=int(len(v))))


def cmd_rna_qc_plots(args) -> int:
    d = out_dir(args, "rna_qc_plots")
    df = _read_metadata(args.metrics)
    ru, vu, pu = barcode_rank(df["total_counts"].to_numpy(), args.umi_cutoff)
    rg, vg, pg = barcode_rank(df["genes"].to_numpy(), args.gene_cutoff)
    if not args.no_plots:
        _rank_panels(ru, vu, pu, "RNA UMIs per Barcode", "RNA UMIs per Top-Ranked Barcode", "log10(UMIs)", 100000,
                     d / "rna_qc_umi_barcode_rank.png")
        _rank_panels(rg, vg, pg, "RNA Genes per Barcode", "RNA Genes per Top-Ranked Barcode", "log10(genes)", 10000,
                     d / "rna_qc_gene_barcode_rank.png")
        plt = _plt()
        if plt is not None:
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.scatter(df["total_counts"], df["genes"], s=4, lw=0, color="darkblue")
            _style(ax, "RNA Genes vs UMIs", "UMIs", "Genes")
            _save(fig, d / "rna_qc_genes_vs_umis.png")
    summ = {"umi_cutoff": args.umi_cutoff, "gene_cutoff": args.gene_cutoff, "n_barcodes_umi": int(len(vu)),
            "n_barcodes_gene": int(len(vg))}
    summ.update({f"umi_{k}": v for k, v in _points_dict(pu).items()})
    summ.update({f"gene_{k}": v for k, v in _points_dict(pg).items()})
    return finish(d, "RNA barcode-rank QC", summ)


# ---------------------------------------------------------------------------
# plot_insert_size_hist.py
# ---------------------------------------------------------------------------

def read_insert_size_hist(path):
    pd = _pd()
    begin, ins, cnt = False, [], []
    with open_text(path) as fh:
        for ln in fh:
            v = ln.rstrip("\n").split("\t")
            if begin and len(v) == 2:
                ins.append(int(v[0]))
                cnt.append(int(float(v[1])))
            elif v[0] == "insert_size":
                begin = True
    return pd.DataFrame({"insert_size": ins, "count": cnt})


def cmd_insert_size_hist(args) -> int:
    d = out_dir(args, "insert_size")
    df = read_insert_size_hist(args.histogram)
    write_tsv(df, d / "insert_size_histogram.tsv")
    plt = None if args.no_plots else _plt()
    if plt is not None and len(df):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.fill_between(df["insert_size"], df["count"], color=ORANGE, alpha=0.8, lw=0)
        ax.plot(df["insert_size"], df["count"], color=ORANGE)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        _style(ax, f"Insert Size Histogram ({args.pkr})", "Insert size", "Count")
        _save(fig, d / "insert_size_histogram.png")
    tot = int(df["count"].sum()) if len(df) else 0
    mode = int(df.loc[df["count"].idxmax(), "insert_size"]) if len(df) else None
    return finish(d, f"Insert size histogram ({args.pkr})", {"n_pairs": tot, "mode_insert_size": mode, "n_sizes": len(df)})


# ---------------------------------------------------------------------------
# 10x_create_barcode_mapping.wdl
# ---------------------------------------------------------------------------

def tenx_barcode_mapping(atac_onlist, rna_onlist) -> "Optional[List[Tuple[str, str]]]":
    atac, rna = read_lines(atac_onlist), read_lines(rna_onlist)
    if len(atac) != len(rna):
        return None
    rows = [(r, a.translate(REV_COMP_FULL)[::-1]) for a, r in zip(atac, rna)]
    rows += [(r, a) for a, r in zip(atac, rna)]
    return rows


def cmd_tenx_barcode_map(args) -> int:
    d = out_dir(args, "tenx_barcode_map")
    rows = tenx_barcode_mapping(args.atac_onlist, args.rna_onlist)
    if rows is None:
        print("Error: ATAC and RNA onlists differ in length -- no conversion dictionary (upstream leaves it undefined).")
        finish(d, "10x barcode mapping (failed)", {"created": False})
        return 1
    p = d / "barcode_conversion_dict.tsv"
    p.write_text("".join(f"{a}\t{b}\n" for a, b in rows))
    print(f"TSV: {p}")
    return finish(d, "10x ATAC <-> RNA barcode mapping", {"created": True, "n_rows": len(rows), "n_barcodes": len(rows) // 2})


# ---------------------------------------------------------------------------
# task_log_atac.wdl
# ---------------------------------------------------------------------------

def _num(s: str):
    s = s.strip()
    if s.isdigit():
        return int(s)
    try:
        return float(s)
    except ValueError:
        return s


def log_atac(alignment_log, barcode_summary) -> dict:
    rows: List[Tuple[str, str]] = []
    text = Path(alignment_log).read_text().splitlines()
    for ln in text:  # grep "Number of" | grep -v threads | tr -d '.' | sed 's/ /_/g' | sed 's/:_/,/g'
        if "Number of" in ln and "threads" not in ln:
            s = ln.replace(".", "").replace(" ", "_").replace(":_", ",")
            k, _, v = s.partition(",")
            rows.append((k, v))
    for ln in text:  # grep "#" | sed 's/, /\n/g' | tr -d '# ' | sed 's/:/,/g' | tr -d '.'
        if "#" in ln:
            for part in ln.split(", "):
                s = part.replace("#", "").replace(" ", "").replace(":", ",").replace(".", "")
                k, _, v = s.partition(",")
                rows.append((k, v))
    tot = dups = unm = low = 0.0
    with open_text(barcode_summary) as fh:
        next(fh, None)
        for ln in fh:
            f = ln.strip().split(",")
            if len(f) >= 5:
                tot += float(f[1]); dups += float(f[2]); unm += float(f[3]); low += float(f[4])
    den = tot - unm - low
    pdup = "%5.1f" % (100 * dups / den) if den else None
    data: Dict[str, Any] = {}
    for k, v in rows:
        data[k] = _num(v)
    data["percentage_duplicates"] = _num(pdup) if pdup is not None else None
    return data


def cmd_log_atac(args) -> int:
    d = out_dir(args, f"log_atac_{args.prefix}")
    data = log_atac(args.alignment_log, args.barcode_summary)
    p = d / f"{args.prefix}_qc_metrics.json"
    p.write_text(json.dumps(data, indent=4))
    print(f"JSON: {p}")
    return finish(d, "ATAC alignment QC metrics", {k: v for k, v in data.items()},
                  "chromap log + barcode summary parsed as in task_log_atac.wdl.")


# ---------------------------------------------------------------------------
# genome TSV / inputs / sampling / portal download
# ---------------------------------------------------------------------------

GENOME_KEYS = {"fasta": "genome FASTA", "kb_nac_idx_tar": "kallisto|bustools nac index tarball",
               "chromap_idx_tar": "chromap index tarball"}


def read_genome_tsv(path) -> "Dict[str, str]":
    """WDL read_map: two tab-separated columns key -> value."""
    m = {}
    for ln in read_lines(path):
        f = ln.split("\t")
        if len(f) >= 2:
            m[f[0]] = f[1]
    return m


def resolve_genome(genome_tsv, fasta=None, kb_index=None, chromap_index=None) -> dict:
    m = read_genome_tsv(genome_tsv) if genome_tsv else {}
    return {"fasta": fasta or m.get("fasta"), "kb_index_tar_gz": kb_index or m.get("kb_nac_idx_tar"),
            "chromap_index_tar_gz": chromap_index or m.get("chromap_idx_tar"), "annotations": m}


def cmd_genome_tsv(args) -> int:
    d = out_dir(args, "genome_tsv")
    r = resolve_genome(args.genome_tsv, args.fasta, args.kb_index, args.chromap_index)
    for k, v in r["annotations"].items():
        print(f"{k}\t{v}")
    missing = [k for k in ("fasta", "kb_index_tar_gz", "chromap_index_tar_gz") if not r[k]]
    return finish(d, "Genome TSV", dict(r, missing=missing), md_table(["key", "value"], r["annotations"].items()))


def classify_input(path: str) -> str:
    if path.startswith("gs://"):
        return "gcs"
    if path.startswith("syn"):
        return "synapse"
    if path.startswith("https"):
        return "https"
    return "local"


def download_file_via_https(url: str, output_dir, dry_run: bool = False) -> "Tuple[Path, str]":
    """IGVF portal download; basic auth from IGVF_API_KEY / IGVF_SECRET_KEY (environment only)."""
    output_dir = Path(output_dir)
    dest = output_dir / Path(url.split("?")[0]).name
    user, pw = os.environ.get("IGVF_API_KEY"), os.environ.get("IGVF_SECRET_KEY")
    auth = "credentials from env" if (user and pw) else "anonymous"
    if dry_run:
        print(f"Dry run: GET {url} ({auth}) -> {dest}")
        return dest, "dry-run"
    import requests  # type: ignore
    output_dir.mkdir(parents=True, exist_ok=True)
    with requests.Session() as s:
        if user and pw:
            s.auth = (user, pw)
        r = s.get(url, stream=True, timeout=60)
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    print(f"Wrote: {dest}")
    return dest, "downloaded"


def check_input(path: str, dest_dir: Path, download: bool) -> dict:
    kind = classify_input(path)
    rec = {"path": path, "kind": kind, "action": "", "resolved": path}
    if kind == "gcs":
        rec["action"] = "kept (WDL passes gs:// URIs through unchanged)"
    elif kind == "local":
        rec["action"] = "exists" if Path(path).exists() else "MISSING"
    elif kind == "synapse":
        cmd = f"cd {shlex.quote(str(dest_dir))} && synapse get {shlex.quote(path)}"
        if download:
            dest_dir.mkdir(parents=True, exist_ok=True)
            ran, _ = run_or_print(cmd, ["synapse"])
            rec["action"] = "downloaded" if ran else "planned (synapse CLI missing)"
        else:
            print(f"Command: {cmd}")
            rec["action"] = "planned"
    else:
        dest, st = download_file_via_https(path, dest_dir, dry_run=not download)
        rec["action"], rec["resolved"] = st, str(dest)
    return rec


def cmd_check_inputs(args) -> int:
    d = out_dir(args, "check_inputs")
    recs = [check_input(p, d / "files", args.download) for p in args.paths]
    pd = _pd()
    write_tsv(pd.DataFrame(recs), d / "inputs.tsv")
    missing = [r["path"] for r in recs if r["action"] == "MISSING"]
    finish(d, "Input check", {"n_inputs": len(recs), "missing": missing, "download": bool(args.download),
                              "kinds": dict(Counter(r["kind"] for r in recs))},
           md_table(["path", "kind", "action"], [(r["path"], r["kind"], r["action"]) for r in recs]))
    return 1 if missing else 0


def sample_fastq(path, dest, n_reads: int = 10_000_000) -> int:
    n = 0
    with open_text(path) as fh, gzip.open(str(dest), "wt", compresslevel=3) as w:
        for i, ln in enumerate(fh):
            if i >= 4 * n_reads:
                break
            w.write(ln)
            n = i + 1
    return n // 4


def cmd_sample_fastqs(args) -> int:
    d = out_dir(args, "sample_fastqs")
    rows = []
    for p in args.fastqs:
        name = Path(p).name
        dest = d / (name if name.endswith(".gz") else name + ".gz")
        n = sample_fastq(p, dest, args.n_reads)
        print(f"Wrote: {dest}")
        rows.append((name, n))
    return finish(d, "Sampled FASTQs", {"n_reads_requested": args.n_reads, "files": [{"file": a, "reads": b} for a, b in rows]},
                  md_table(["file", "reads written"], rows))


def cmd_portal_download(args) -> int:
    d = out_dir(args, "portal_download")
    dest = Path(args.dest) if args.dest else d / "files"
    recs = []
    for url in args.urls:
        p, st = download_file_via_https(url, dest, dry_run=args.dry_run)
        recs.append({"url": url, "dest": str(p), "status": st})
    return finish(d, "IGVF portal download", {"n": len(recs), "dry_run": bool(args.dry_run), "files": recs,
                                              "credentials": "env" if os.environ.get("IGVF_API_KEY") else "none"})


# ---------------------------------------------------------------------------
# kallisto | bustools  (task_kb_count.wdl, task_kb_index.wdl -> atomic-workflows run_kallisto.py)
# ---------------------------------------------------------------------------

def _maybe_gunzip(path: Optional[str], dest: Path) -> Optional[str]:
    """WDL: *.gz onlist / replacement list decompressed to a .txt next to the run."""
    if not path:
        return None
    if str(path).endswith(".gz"):
        dest.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "rb") as fi, open(dest, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return str(dest)
    return str(path)


def interleave_fastqs(read1: "List[str]", read2: "List[str]", barcode: "Optional[List[str]]" = None) -> "List[str]":
    """paste r1 r2 [bc] per lane, in lane order (l1r1 l1r2 l1bc l2r1 ...)."""
    barcode = barcode or []
    out = []
    for i in range(max(len(read1), len(read2), len(barcode))):
        for lst in (read1, read2, barcode):
            if i < len(lst) and lst[i]:
                out.append(lst[i])
    return out


def kb_count_command(kb_mode: str, index_dir: str, read_format: str, output_dir: str, strand: str, onlist: str,
                     fastqs: "List[str]", threads: int = 4, replacement_list: Optional[str] = None,
                     temp_dir: Optional[str] = None) -> str:
    q = shlex.quote
    tmp = f"--tmp {q(temp_dir)} " if temp_dir else ""
    rep = f"-r {q(replacement_list)} " if replacement_list else ""
    if kb_mode == "nac":
        idx = f"-i {q(index_dir + '/index.idx')} -g {q(index_dir + '/t2g.txt')} -c1 {q(index_dir + '/cdna.txt')} -c2 {q(index_dir + '/nascent.txt')} --sum=total"
    else:  # atomic-workflows v1.1 only implements nac; standard is the documented kb equivalent
        idx = f"-i {q(index_dir + '/index.idx')} -g {q(index_dir + '/t2g.txt')}"
    return (f"kb count --workflow={kb_mode} {tmp}{idx} -x {q(read_format)} -w {q(onlist)} {rep}"
            f"--strand {q(strand)} -o {q(output_dir)} --h5ad -t {threads} " + " ".join(q(f) for f in fastqs))


def kb_postprocess(output_dir: str, subpool: Optional[str], replacement_list: Optional[str]) -> dict:
    """Append the subpool to counts_unfiltered[_modified]/adata.h5ad + barcodes.txt; move h5ad to <out>.h5ad."""
    sfx = "_modified" if replacement_list else ""
    cu = Path(output_dir) / f"counts_unfiltered{sfx}"
    h5 = cu / "adata.h5ad"
    res = {"h5ad": None, "subpool_applied": False}
    if subpool and subpool != "none" and h5.exists():
        append_suffix_to_h5ad(h5, subpool)
        bt = cu / "cells_x_genes.barcodes.txt"
        if bt.exists():
            bt.write_text("".join(f"{ln}_{subpool}\n" for ln in bt.read_text().splitlines()))
        res["subpool_applied"] = True
    if h5.exists():
        dest = Path(str(output_dir).rstrip("/") + ".h5ad")
        shutil.move(str(h5), str(dest))
        res["h5ad"] = str(dest)
    return res


def cmd_kb_count(args) -> int:
    d = out_dir(args, "kb_count")
    onlist = _maybe_gunzip(args.onlist, d / "barcode_inclusion_list.txt")
    repl = _maybe_gunzip(args.replacement_list, d / "replacement_list.txt")
    outp = str(Path(args.output_dir) if args.output_dir else d / "kb_out")
    fq = interleave_fastqs(args.read1, args.read2, args.read_barcode)
    cmd = kb_count_command(args.kb_mode, args.index_dir, args.read_format, outp, args.strand, onlist, fq,
                           args.threads, repl, args.temp_dir)
    Path(outp).mkdir(parents=True, exist_ok=True)
    ran, _ = run_or_print(cmd, ["kb"])
    post = kb_postprocess(outp, args.subpool, repl) if ran else {}
    return finish(d, "kallisto | bustools count", {"ran": ran, "command": cmd, "kb_mode": args.kb_mode,
                                                   "n_fastqs": len(fq), **post})


def kb_index_command(kb_mode: str, genome_fasta: str, gtf: str, output_dir: str, temp_dir: Optional[str] = None) -> str:
    q = shlex.quote
    tmp = f"--tmp {q(temp_dir)} " if temp_dir else ""
    o = output_dir
    if kb_mode == "nac":
        return (f"kb ref {tmp}--workflow=nac -i {q(o + '/index.idx')} -g {q(o + '/t2g.txt')} -c1 {q(o + '/cdna.txt')} "
                f"-c2 {q(o + '/nascent.txt')} -f1 {q(o + '/cdna.fasta')} -f2 {q(o + '/nascent.fasta')} {q(genome_fasta)} {q(gtf)}")
    return f"kb ref {tmp}--workflow=standard -i {q(o + '/index.idx')} -g {q(o + '/t2g.txt')} -f1 {q(o + '/cdna.fasta')} {q(genome_fasta)} {q(gtf)}"


def cmd_kb_index(args) -> int:
    d = out_dir(args, "kb_index")
    o = str(Path(args.output_dir) if args.output_dir else d / "kb_index")
    Path(o).mkdir(parents=True, exist_ok=True)
    cmd = kb_index_command(args.kb_mode, args.genome_fasta, args.gtf, o, args.temp_dir)
    ran, _ = run_or_print(cmd, ["kb"])
    if ran:
        run_or_print(f"tar --exclude='*.fasta' -zcvf {shlex.quote(o + '.tar.gz')} -C {shlex.quote(o)} .", ["tar"])
    return finish(d, "kallisto | bustools index", {"ran": ran, "command": cmd})


# ---------------------------------------------------------------------------
# chromap  (task_chromap.wdl, task_chromap_bam.wdl, task_chromap_index.wdl -> run_chromap.py)
# ---------------------------------------------------------------------------

def chromap_align_command(output: str, index_dir: str, read_format: str, reference_fasta: str, prefix: str,
                          onlist: str, read1: "List[str]", read2: "List[str]", read_barcode: "List[str]",
                          threads: int = 8, barcode_translate: Optional[str] = None, subpool: Optional[str] = None) -> str:
    q = shlex.quote
    bt = f"--barcode-translate {q(barcode_translate)} " if barcode_translate else ""
    common = (f"chromap -x {q(index_dir + '/index')} --read-format {q(read_format)} -r {q(reference_fasta)} "
              f"--remove-pcr-duplicates --remove-pcr-duplicates-at-cell-level --trim-adapters --low-mem "
              f"{'--BED' if output == 'fragments' else '--SAM'} -l 2000 --bc-error-threshold 1 -t {threads} "
              f"--bc-probability-threshold 0.90 -q 30 --barcode-whitelist {q(onlist)} {bt}")
    reads = f"-1 {q(','.join(read1))} -2 {q(','.join(read2))} -b {q(','.join(read_barcode))}"
    if output == "fragments":
        return common + f"-o {q(prefix + '.fragments.tsv')} --summary {q(prefix + '.barcode.summary.csv')} {reads} > {q(prefix + '.log.txt')} 2>&1"
    awk = ""
    if subpool and subpool != "none":
        awk = ("| awk -v suffix=" + q(str(subpool)) +
               " '{OFS=\"\\t\"; for (i=12; i<=NF; i++) if ($i ~ /^CB:Z:/) $i = $i\"_\"suffix; print $0}' ")
    return common + f"--summary {q(prefix + '.barcode.summary.csv')} {reads} -o /dev/stdout {awk}| samtools view -bS - > {q(prefix + '.bam')}"


def subpool_fragments(fragments, subpool: str) -> int:
    """run_chromap.process_fragments: column 4 += _<subpool> (rows with >= 4 columns)."""
    tmp = Path(str(fragments) + ".tmp")
    n = 0
    with open_text(fragments) as fi, open(tmp, "w") as fo:
        for ln in fi:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 4:
                f[3] = f"{f[3]}_{subpool}"
                n += 1
            fo.write("\t".join(f) + "\n")
    shutil.move(str(tmp), str(fragments))
    return n


def subpool_summary(summary_csv, subpool: str) -> int:
    """run_chromap.process_summary: column 1 += _<subpool> except the header."""
    rows = list(csv.reader(open(summary_csv, newline="")))
    with open(summary_csv, "w", newline="") as fo:
        w = csv.writer(fo)
        for i, r in enumerate(rows):
            if i and r:
                r[0] = f"{r[0]}_{subpool}"
            w.writerow(r)
    return max(len(rows) - 1, 0)


def bgzip_tabix(fragments_tsv) -> Optional[str]:
    """bgzip -c + tabix --zero-based --preset bed (pysam if available, else the binaries)."""
    out = str(fragments_tsv) + ".gz"
    try:
        import pysam  # type: ignore
        pysam.tabix_compress(str(fragments_tsv), out, force=True)
        pysam.tabix_index(out, preset="bed", zerobased=True, force=True)
        return out
    except ImportError:
        ran, _ = run_or_print(f"bgzip -c {shlex.quote(str(fragments_tsv))} > {shlex.quote(out)} && "
                              f"tabix --zero-based --preset bed {shlex.quote(out)}", ["bgzip", "tabix"])
        return out if ran else None


def cmd_chromap_align(args) -> int:
    d = out_dir(args, f"chromap_{args.output}")
    onlist = _maybe_gunzip(args.onlist, d / "barcode_inclusion_list.txt")
    prefix = str(d / args.prefix)
    cmd = chromap_align_command(args.output, args.index_dir, args.read_format, args.reference_fasta, prefix, onlist,
                                args.read1, args.read2, args.read_barcode, args.threads, args.barcode_translate, args.subpool)
    need = ["chromap"] + (["samtools"] if args.output == "bam" else [])
    ran, _ = run_or_print(cmd, need)
    res = {"ran": ran, "command": cmd, "output": args.output}
    if ran:
        sp = args.subpool if args.subpool and args.subpool != "none" else None
        if args.output == "fragments":
            if sp:
                subpool_fragments(prefix + ".fragments.tsv", sp)
            res["fragments"] = bgzip_tabix(prefix + ".fragments.tsv")
        else:
            run_or_print(f"samtools index {shlex.quote(prefix + '.bam')}", ["samtools"])
            res["bam"] = prefix + ".bam"
        if sp and Path(prefix + ".barcode.summary.csv").exists():
            subpool_summary(prefix + ".barcode.summary.csv", sp)
    return finish(d, f"chromap align {args.output}", res)


def cmd_chromap_index(args) -> int:
    d = out_dir(args, "chromap_index")
    o = str(Path(args.output_dir) if args.output_dir else d / "chromap_index")
    Path(o).mkdir(parents=True, exist_ok=True)
    fasta = args.genome_fasta
    if fasta.endswith(".gz") and which("chromap"):
        fasta = _maybe_gunzip(fasta, Path(o).parent / Path(fasta).name[:-3])
    cmd = f"chromap -i -r {shlex.quote(fasta)} -o {shlex.quote(o + '/index')}"
    ran, _ = run_or_print(cmd, ["chromap"])
    if ran:
        run_or_print(f"tar -zcvf {shlex.quote(o + '.tar.gz')} -C {shlex.quote(o)} .", ["tar"])
    return finish(d, "chromap index", {"ran": ran, "command": cmd})


def cmd_subpool_fragments(args) -> int:
    d = out_dir(args, "subpool")
    res = {"subpool": args.subpool}
    if args.fragments:
        res["fragments_rows"] = subpool_fragments(args.fragments, args.subpool)
        print(f"Wrote: {args.fragments}")
        if args.bgzip:
            res["fragments_gz"] = bgzip_tabix(args.fragments)
            if res["fragments_gz"]:
                print(f"Wrote: {res['fragments_gz']}")
    if args.summary:
        res["summary_rows"] = subpool_summary(args.summary, args.subpool)
        print(f"Wrote: {args.summary}")
    return finish(d, "Subpool suffix on ATAC barcodes", res)


# ---------------------------------------------------------------------------
# Synapse (upload-to-synapse-y3-results.py, batch_upload_synapse.py, update-file-annotations.py)
# ---------------------------------------------------------------------------

TO_BE_UPLOADED = ["seqspec_atac_onlist_renamed", "seqspec_rna_onlist_renamed", "atac_bam", "atac_bam_log",
                  "atac_chromap_barcode_metadata", "atac_filter_fragments", "atac_filter_fragments_index",
                  "atac_snapatac2_barcode_metadata", "csv_summary", "html_summary", "joint_barcode_metadata",
                  "rna_barcode_metadata", "rna_aggregated_counts_h5ad", "rna_kb_output", "rna_log", "rna_mtx_tar",
                  "rna_mtxs_h5ad"]
ACTIVITY = {
    "seqspec_atac_onlist_renamed": ("Generating onlist for ATAC", "task_seqspec_extract.wdl"),
    "seqspec_rna_onlist_renamed": ("Generating onlist for RNA", "task_seqspec_extract.wdl"),
    "atac_bam": ("Align ATAC and return bam and log", "task_chromap_bam.wdl"),
    "atac_bam_log": ("Align ATAC and return bam log", "task_chromap_bam.wdl"),
    "atac_chromap_barcode_metadata": ("Barcode metadata statistics from Chromap", "task_chromap.wdl"),
    "atac_filter_fragments": ("Raw fragment file from Chromap", "task_chromap.wdl"),
    "atac_filter_fragments_index": ("Index for raw fragment file from Chromap", "task_chromap.wdl"),
    "atac_snapatac2_barcode_metadata": ("Filtered barcode metadata statistics from SnapATAC2", "task_qc_atac.wdl"),
    "csv_summary": ("Pipeline summary in CSV file format", "task_html_report.wdl"),
    "html_summary": ("Pipeline summary in HTML file format", "task_html_report.wdl"),
    "joint_barcode_metadata": ("Joint barcode metadata", "task_joint_qc.wdl"),
    "reference_fasta": ("Reference fasta file", None),
    "rna_barcode_metadata": ("Barcode metadata statistics from KB", "task_kb_count.wdl"),
    "rna_aggregated_counts_h5ad": ("RNA count matrix in h5ad format", "task_kb_count.wdl"),
    "rna_kb_output": ("Total outptu from KB in TAR format", "task_kb_count.wdl"),
    "rna_log": ("Log file from KB", "task_log_rna.wdl"),
    "rna_mtx_tar": ("RNA count matrices in TAR format", "task_kb_count.wdl"),
    "rna_mtxs_h5ad": ("RNA count matrices in h5ad format", "task_kb_count.wdl"),
}
_ATAC_RAW = ["ATAC_barcode", "ATAC_fastq_R1", "ATAC_fastq_R2"]
_RNA_RAW = ["RNA_fastq_R1", "RNA_fastq_R2"]
USED = {
    "seqspec_atac_onlist_renamed": _ATAC_RAW + ["seqspec*"],
    "seqspec_rna_onlist_renamed": _RNA_RAW + ["seqspec*"],
    **{c: _ATAC_RAW + ["seqspec_atac_onlist_renamed"] for c in ("atac_bam", "atac_bam_log", "atac_chromap_barcode_metadata",
                                                                 "atac_filter_fragments", "atac_filter_fragments_index")},
    "atac_snapatac2_barcode_metadata": ["atac_filter_fragments", "atac_filter_fragments_index", "seqspec_atac_onlist_renamed"],
    "joint_barcode_metadata": ["atac_snapatac2_barcode_metadata", "rna_barcode_metadata"],
    **{c: _RNA_RAW + ["seqspec_rna_onlist_renamed"] for c in ("rna_barcode_metadata", "rna_aggregated_counts_h5ad",
                                                              "rna_kb_output", "rna_log", "rna_mtx_tar", "rna_mtxs_h5ad")},
}


def clean_string(s: str, local_root: str, remote_root: str) -> str:
    s = str(s).replace("]", "").replace("[", "").replace(",", ";").replace('"', "")
    return s.replace(remote_root, local_root) if remote_root else s


def activity_and_script(col: str) -> "Tuple[str, str]":
    if col not in ACTIVITY:
        return "", ""
    act, wdl = ACTIVITY[col]
    return f'"{act}"', f'"{SCRIPT_BASE + wdl}"' if wdl else '"input file"'


def synapse_manifest(table, project: str, local_root: str = "", remote_root: str = "", syn=None) -> "List[dict]":
    """Manifest rows path/parent/name/used/executed/activityName.  Without a Synapse client (dry run) the
    per-subpool folder is written as '<project>/<subpool>' and `used` keeps the paths."""
    manifest = []
    for _, row in table.iterrows():
        base = clean_string(row["Subpool"], local_root, remote_root)
        if syn is not None:
            from synapseclient import Folder  # type: ignore
            folder = syn.findEntityId(name=base, parent=project) or syn.store(Folder(name=base, parent=project))
        else:
            folder = f"{project}/{base}"
        for col in table.columns:
            v = row[col]
            if col not in TO_BE_UPLOADED or (isinstance(v, float) and math.isnan(v)):
                continue
            act, script = activity_and_script(col)
            used_cols = USED.get(col)
            if used_cols is None:
                used = '""'
            else:
                items = []
                for c in used_cols:
                    if c == "seqspec*":
                        items.append(clean_string(row.get("seqspec", ""), "", ""))
                    elif c in row.index:
                        items.append(clean_string(row[c], local_root, remote_root))
                if syn is not None:
                    items = [syn.findEntityId(name=Path(p).name.replace(".filter", ""), parent=folder) or p for p in items]
                used = '"' + ";".join(items) + '"'
            name = os.path.basename(str(v)).replace(".filter", "")
            if syn is not None and syn.findEntityId(name=name, parent=folder):
                continue
            manifest.append({"path": str(v).replace(remote_root, local_root) if remote_root else str(v), "parent": folder,
                             "name": name, "used": used, "executed": script, "activityName": act})
    return manifest


def _synapse_client():
    try:
        import synapseclient  # type: ignore
    except Exception:
        return None, "synapseclient not installed (pip install synapseclient)"
    tok = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not tok:
        return None, "SYNAPSE_AUTH_TOKEN not set"
    syn = synapseclient.Synapse()
    syn.login(authToken=tok, silent=True)
    return syn, ""


def cmd_synapse_manifest(args) -> int:
    d = out_dir(args, "synapse_manifest")
    pd = _pd()
    table = pd.read_csv(args.table, sep="\t")
    syn, why = (_synapse_client() if args.execute else (None, "dry run"))
    if args.execute and syn is None:
        print(f"Note: {why} -- writing a dry-run manifest instead.")
    rows = synapse_manifest(table, args.project, args.local_root or "", args.remote_root or "", syn)
    p = Path(args.output) if args.output else d / "manifest.tsv"
    with open(p, "w", newline="") as f:
        cols = ["path", "parent", "name", "used", "executed", "activityName"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar=None, escapechar="\\")
        w.writeheader()
        w.writerows(rows)
    print(f"TSV: {p}")
    return finish(d, "Synapse upload manifest", {"n_rows": len(rows), "n_samples": len(table), "dry_run": syn is None,
                                                 "project": args.project})


def manifest_groups(manifest) -> "List[Tuple[str, List[str]]]":
    """batch_upload_synapse: consecutive rows with the same parent (column 2) form one upload batch."""
    with open(manifest) as f:
        header = f.readline().rstrip("\n")
        groups: List[Tuple[str, List[str]]] = []
        for ln in f:
            pid = ln.split("\t")[1]
            if not groups or groups[-1][0] != pid:
                groups.append((pid, []))
            groups[-1][1].append(ln)
    return groups


def cmd_synapse_upload(args) -> int:
    d = out_dir(args, "synapse_upload")
    header = open(args.manifest).readline()
    groups = manifest_groups(args.manifest)
    syn, why = (_synapse_client() if args.execute else (None, "dry run (pass --execute to upload)"))
    failed, done = [], []
    for pid, lines in groups:
        tmp = d / f"batch_{safe_label(pid)}.tsv"
        tmp.write_text(header + "".join(lines))
        if syn is None:
            print(f"Dry run: would syncToSynapse {len(lines)} file(s) into {pid} ({tmp.name})")
            continue
        try:
            import synapseutils  # type: ignore
            synapseutils.syncToSynapse(syn=syn, manifestFile=str(tmp), sendMessages=False)
            done.append(pid)
        except Exception as e:  # upstream: collect failures, keep going
            log.warning("upload failed for %s: %s", pid, e)
            failed.append(pid)
    fp = d / "failed_ids.txt"
    fp.write_text("".join(f"{x}\n" for x in failed))
    print(f"Wrote: {fp}")
    if syn is None and args.execute:
        print(f"Note: {why}")
    return finish(d, "Synapse batch upload", {"n_batches": len(groups), "uploaded": len(done), "failed": len(failed),
                                              "dry_run": syn is None})


def file_type_status_main(name: str) -> "Tuple[Optional[str], Optional[str], Optional[str]]":
    rules = [(".bam", ("bam", "raw", "yes")), (".log.txt", ("txt", "raw", "no")), (".tsv.gz", ("tsv.gz", "raw", "yes")),
             (".tsv.gz.tbi", (".tbi", "raw", "yes")), (".onlist.txt.gz", ("txt.gz", "raw", "no")),
             (".chromap.barcode.metadata.tsv", ("tsv", "raw", "no")), (".snapatac2.barcode.metadata.tsv", ("tsv", "filtered", "no")),
             (".csv", ("csv", "filtered", "no")), (".html", ("html", "filtered", "no")),
             (".cells_x_genes.total.h5ad", ("h5ad", "raw", "yes")), (".h5ad", ("h5ad", "raw", "no")),
             (".tar.gz", ("tar.gz", "raw", "no")), (".tsv", ("tsv", "raw", "no"))]
    for suf, val in rules:
        if name.endswith(suf):
            return val
    if name.endswith(".txt") and "rna.log" in name:
        return "txt", "raw", "no"
    if name.endswith("seqspec.txt.gz"):
        return "txt.gz", "processed", "no"
    if name.endswith(".yaml"):
        return "yaml", "raw", "no"
    return None, None, None  # upstream returns a 2-tuple here and crashes on unpacking


def file_data_type(name: str) -> str:
    return "atac" if ".atac" in name else "rna" if ".rna" in name else "joint"


def file_description(name: str) -> Optional[str]:
    n = name
    if n.endswith(".tsv.gz"):
        return "[Raw] Fragment file from Chromap"
    if n.endswith(".tsv.gz.tbi"):
        return "[Raw] Fragment file index from Chromap"
    if n.endswith("atac.onlist.txt.gz"):
        return "Barcode onlist file generated by seqspec"
    if n.endswith(".chromap.barcode.metadata.tsv"):
        return "[Raw] Per barcode alignment statistics file from Chromap"
    if n.endswith(".snapatac2.barcode.metadata.tsv"):
        return "[Filtered] Per barcode statistics file from SnapATAC2"
    if n.endswith(".csv") and "joint" in n:
        return "[Filtered] Joint per barcode statistics file"
    if n.endswith(".html"):
        return "[Filtered] HTML summary report"
    if n.endswith(".csv"):
        return "[Filtered] CSV summary report"
    if n.endswith(".cells_x_genes.total.h5ad"):
        return "[Raw] Aggregated(Ambiguous+Spliced+Unspliced) count matrix in h5ad format"
    if n.endswith(".count_matrix.h5ad"):
        return "[Raw] h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total"
    if n.endswith(".mtx.tar.gz"):
        return "[Raw] Tarball containing four separated count matrices in mtx format: Spliced, Unspliced, Ambiguous, and Total"
    if n.endswith(".tar.gz") and "rna.align.kb" in n:
        return "[Raw] Tarball containing all the logs and bus files generated from kb"
    if n.endswith(".txt") and "rna.log" in n:
        return "[Raw] Log file from kb"
    if n.endswith("rna.onlist.txt.gz"):
        return "Barcode onlist file generated by seqspec"
    if n.endswith(".log.txt") and "atac.align" in n:
        return "[Raw] Log file from aligner"
    if n.endswith(".bam"):
        return "[Raw] Aligned bam file from Chromap"
    if "qc" in n and "rna" in n:
        return "[Raw] Per barcode alignment statistics file"
    if n.endswith(".yaml"):
        return "[Raw] seqspec file"
    if n.endswith("_onlist_seqspec.txt.gz") and "atac" in n:
        return "[Processed] ATAC barcode onlist file generated by seqspec"
    if n.endswith("_onlist_seqspec.txt.gz") and "rna" in n:
        return "[Processed] RNA barcode onlist file generated by seqspec"
    return None


def file_annotations(name: str) -> dict:
    ft, st, main = file_type_status_main(name)
    return {"file_type": ft, "status": st, "data_type": file_data_type(name), "file_description": file_description(name),
            "main_output": main}


def cmd_synapse_annotations(args) -> int:
    d = out_dir(args, "synapse_annotations")
    names = list(args.names or [])
    if args.names_file:
        names += read_lines(args.names_file)
    pd = _pd()
    df = pd.DataFrame([dict(name=n, **file_annotations(n)) for n in names])
    write_tsv(df, d / "annotations.tsv")
    applied = 0
    if args.execute and args.root_folder:
        syn, why = _synapse_client()
        if syn is None:
            print(f"Note: {why} -- annotations not applied.")
        else:
            import synapseclient  # type: ignore
            import synapseutils  # type: ignore
            for _, dirs, _ in synapseutils.walk(syn=syn, synId=args.root_folder, includeTypes=["file"]):
                for _, sub_id in dirs:
                    for _, _, files in synapseutils.walk(syn=syn, synId=sub_id, includeTypes=["file"]):
                        for fname, fid in files:
                            ann = syn.get_annotations(fid)
                            syn.set_annotations(synapseclient.Annotations(ann.id, ann.etag, file_annotations(fname)))
                            applied += 1
                        break
    return finish(d, "Synapse file annotations", {"n_files": len(df), "applied": applied,
                                                  "unrecognised": int(df["file_type"].isna().sum()) if len(df) else 0})


# ---------------------------------------------------------------------------
# write_html.py
# ---------------------------------------------------------------------------

def format_number(txt_num) -> str:
    num = float(txt_num)
    if num > 1_000_000_000:
        return f"{int(num / 1_000_000_000)} B"
    if num > 1_000_000:
        return f"{int(round(num, -6) / 1_000_000)} M"
    if num > 1000:
        return f"{int(num / 1000)} K"
    if num < 1:
        return str(round(num))
    return str(num)


def _read_name_value(path) -> "List[Tuple[str, str]]":
    out = []
    for ln in read_lines(path):
        f = ln.split(",")
        out.append((f[0], f[1] if len(f) > 1 else ""))
    return out


def write_html_report(html_path, images: "List[str]", logs: "List[str]", csv_path,
                      atac_metrics: Optional[str] = None, rna_metrics: Optional[str] = None) -> dict:
    def short(p):
        b = os.path.basename(p)
        return b[b.index(".") + 1:] if "." in b else b
    tabs = {"joint": [], "rna": [], "atac": []}
    for im in images:
        n = short(im)
        tabs["atac" if "atac" in n else "rna" if "rna" in n else "joint"].append(im)
    csv_lines: List[str] = []
    h = ['<!DOCTYPE html><html lang="en"><head><title>Results summary</title><meta charset="UTF-8">'
         '<meta name="viewport" content="width=device-width, initial-scale=1.0"><style>'
         'input{display:none} input+label{display:inline-block;border:1px solid #999;background:#EEE;padding:4px 12px;'
         'border-radius:4px 4px 0 0;position:relative;top:1px} input:checked+label{background:#FFF;border-bottom:1px solid transparent}'
         ' input~.tab{display:none;border-top:1px solid #999;padding:12px}'
         + "".join(f"#tab{i}:checked~.tab.content{i}," for i in range(1, 5)) + "#tab5:checked~.tab.content5{display:block}"
         ' img{max-width:100%} table{border-spacing:24px 0}</style></head><body>']
    for i, lab in enumerate(["Joint Plots", "RNA plots", "ATAC plots", "Statistics", "Log files"], 1):
        h.append(f'<input type="radio" name="tabs" id="tab{i}"{" checked" if i == 1 else ""} /><label for="tab{i}">{lab}</label>')

    def pngs(lst):
        for im in lst:
            b64 = base64.b64encode(Path(im).read_bytes()).decode()
            tag = f'<img id="{short(im)}" width="1000" src="data:image/png;base64,{b64}" alt="{os.path.basename(im)}"><br>'
            h.append(tag)
            csv_lines.append(f"{short(im)},{tag}")
    for i, key in enumerate(("joint", "rna", "atac"), 1):
        h.append(f'<div class="tab content{i}">')
        pngs(tabs[key])
        h.append('<a href="#top">Go to top of page</a></div>')
    stats = {}
    if atac_metrics:
        stats["ATAC"] = _read_name_value(atac_metrics)
    if rna_metrics:
        stats["RNA"] = _read_name_value(rna_metrics)
    h.append('<div class="tab content4"><table><th>Summary Statistics</th>')
    for mod, rows in stats.items():
        h.append(f"<tr><td colspan=2>{mod}</td></tr>")
        for k, v in rows:
            csv_lines.append(f"{k},{v}")
            try:
                shown = format_number(float(v))
            except ValueError:
                shown = v
            h.append(f"<tr><td>{k}</td><td id={k}>{shown}</td></tr>")
    h.append("</table></div><div class=\"tab content5\">")
    for lg in logs:
        n = short(lg.split("/")[-1])
        h.append(f"<div id={n}>{lg}</div><br>")
        csv_lines.append(f"{n},{lg}")
    h.append("</div></body></html>")
    for rows in stats.values():  # upstream writes the metrics a second time at the end of the CSV
        csv_lines.extend(f"{k},{v}" for k, v in rows)
    Path(html_path).write_text("".join(h), encoding="utf8")
    Path(csv_path).write_text("".join(ln + "\n" for ln in csv_lines))
    return {"n_images": len(images), "n_joint": len(tabs["joint"]), "n_rna": len(tabs["rna"]), "n_atac": len(tabs["atac"]),
            "n_logs": len(logs), "n_metrics": sum(len(v) for v in stats.values())}


def cmd_html_report(args) -> int:
    d = out_dir(args, "html_report")
    images = list(args.images or []) + (read_lines(args.image_list) if args.image_list else [])
    logs = list(args.logs or []) + (read_lines(args.log_list) if args.log_list else [])
    html, csvp = d / f"{args.prefix}.html", d / f"{args.prefix}.csv"
    r = write_html_report(html, images, logs, csvp, args.atac_metrics, args.rna_metrics)
    print(f"Wrote: {html}")
    print(f"CSV: {csvp}")
    return finish(d, "HTML summary report", r)


# ---------------------------------------------------------------------------
# single_cell_pipeline.wdl  (plan + optional execution + downstream QC)
# ---------------------------------------------------------------------------

WDL_KEYS = {  # single_cell_pipeline.<key> -> argparse dest
    "prefix": "prefix", "subpool": "subpool", "genome_tsv": "genome_tsv", "create_onlist_mapping": "create_onlist_mapping",
    "atac_read1": "atac_read1", "atac_read2": "atac_read2", "fastq_barcode": "fastq_barcode",
    "atac_barcode_inclusion_list": "atac_onlist", "chromap_genome_index_tar_gz": "chromap_index", "genome_fasta": "genome_fasta",
    "atac_read_format": "atac_read_format", "rna_read1": "rna_read1", "rna_read2": "rna_read2",
    "fastq_barcode_rna": "rna_barcode", "rna_barcode_inclusion_list": "rna_onlist", "kb_mode": "kb_mode",
    "rna_read_format": "rna_read_format", "kb_genome_index_tar_gz": "kb_index", "rna_replacement_list": "replacement_list",
    "kb_strand": "kb_strand",
}


PIPELINE_DEFAULTS = {"prefix": "sample", "subpool": "none", "kb_mode": "nac", "kb_strand": "unstranded"}


def _apply_inputs_json(args) -> None:
    data = json.loads(Path(args.inputs_json).read_text())
    for k, v in data.items():
        key = k.split(".", 1)[1] if "." in k else k
        dest = WDL_KEYS.get(key)
        if dest and (getattr(args, dest, None) in (None, [], False)):
            setattr(args, dest, v)


def _as_list(v) -> "List[str]":
    if v is None:
        return []
    if isinstance(v, str):
        return [v] if v else []
    return [str(x) for x in v]


def cmd_pipeline(args) -> int:
    if args.inputs_json:
        _apply_inputs_json(args)
    for k, v in PIPELINE_DEFAULTS.items():
        if getattr(args, k, None) in (None, ""):
            setattr(args, k, v)
    d = out_dir(args, f"pipeline_{args.prefix}")
    steps: List[dict] = []
    q = shlex.quote
    atac1, atac2, atacb = _as_list(args.atac_read1), _as_list(args.atac_read2), _as_list(args.fastq_barcode)
    rna1, rna2, rnab = _as_list(args.rna_read1), _as_list(args.rna_read2), _as_list(args.rna_barcode)
    process_atac, process_rna = bool(atac1 and atac1[0]), bool(rna1 and rna1[0])
    g = resolve_genome(args.genome_tsv, args.genome_fasta, args.kb_index, args.chromap_index)
    steps.append({"step": "genome_tsv", "fasta": g["fasta"], "kb_index": g["kb_index_tar_gz"], "chromap_index": g["chromap_index_tar_gz"]})
    # check_inputs: every non-gs:// input (WDL checks the first element's scheme, then scatters)
    inputs = atac1 + atac2 + atacb + rna1 + rna2 + rnab + [x for x in (args.replacement_list, g["fasta"],
                                                                     g["kb_index_tar_gz"], g["chromap_index_tar_gz"]) if x]
    checks = [check_input(p, d / "inputs", args.download) for p in inputs]
    steps.append({"step": "check_inputs", "n": len(checks), "missing": [c["path"] for c in checks if c["action"] == "MISSING"]})
    resolved = {c["path"]: c["resolved"] for c in checks}
    R = lambda lst: [resolved.get(p, p) for p in lst]  # noqa: E731
    tenx_map = None
    if args.create_onlist_mapping and process_atac and process_rna:
        rows = tenx_barcode_mapping(args.atac_onlist, args.rna_onlist)
        if rows is not None:
            tenx_map = d / "barcode_conversion_dict.tsv"
            tenx_map.write_text("".join(f"{a}\t{b}\n" for a, b in rows))
            print(f"TSV: {tenx_map}")
        steps.append({"step": "mapping_tenx_barcodes", "created": tenx_map is not None})
    rna_h5ad = args.rna_h5ad
    if process_rna:
        out = str(d / args.prefix / "rna")
        idx_dir = str(d / "rna_index_folder")
        onl = _maybe_gunzip(args.rna_onlist, d / "rna_barcode_inclusion_list.txt") if args.rna_onlist else "barcode_inclusion_list.txt"
        rep = _maybe_gunzip(args.replacement_list, d / "replacement_list.txt")
        cmd_tar = f"mkdir -p {q(idx_dir)} {q(out)} && tar xvzf {q(str(g['kb_index_tar_gz']))} --no-same-owner -C {q(idx_dir)}"
        cmd = kb_count_command(args.kb_mode, idx_dir, args.rna_read_format or "", out, args.kb_strand, onl,
                               interleave_fastqs(R(rna1), R(rna2), R(rnab)), args.threads, rep)
        ran = False
        if which("kb") and g["kb_index_tar_gz"] and Path(str(g["kb_index_tar_gz"])).exists():
            run_or_print(cmd_tar, ["tar"])
            ran, _ = run_or_print(cmd, ["kb"])
            if ran:
                post = kb_postprocess(out, args.subpool, rep)
                rna_h5ad = rna_h5ad or post.get("h5ad")
        else:
            print(f"Command: {cmd_tar}")
            run_or_print(cmd, ["kb"]) if not which("kb") else print(f"Command: {cmd}")
        steps.append({"step": "kb_count", "ran": ran, "commands": [cmd_tar, cmd]})
    fragments = args.fragments
    if process_atac:
        prefix = str(d / f"{args.prefix}.atac")
        idx_dir = str(d / "atac_index_folder")
        onl = _maybe_gunzip(args.atac_onlist, d / "atac_barcode_inclusion_list.txt") if args.atac_onlist else "barcode_inclusion_list.txt"
        cmd_tar = f"mkdir -p {q(idx_dir)} && tar xvzf {q(str(g['chromap_index_tar_gz']))} --no-same-owner -C {q(idx_dir)}"
        cmds = [chromap_align_command(o, idx_dir, args.atac_read_format or "", str(g["fasta"]), prefix, onl, R(atac1), R(atac2),
                                      R(atacb), args.threads, str(tenx_map) if tenx_map else None, args.subpool)
                for o in ("fragments", "bam")]
        ran = False
        if which("chromap") and g["chromap_index_tar_gz"] and Path(str(g["chromap_index_tar_gz"])).exists():
            run_or_print(cmd_tar, ["tar"])
            ran, _ = run_or_print(cmds[0], ["chromap"])
            if ran:
                sp = args.subpool if args.subpool and args.subpool != "none" else None
                if sp:
                    subpool_fragments(prefix + ".fragments.tsv", sp)
                    subpool_summary(prefix + ".barcode.summary.csv", sp)
                fragments = fragments or bgzip_tabix(prefix + ".fragments.tsv")
                data = log_atac(prefix + ".log.txt", prefix + ".barcode.summary.csv")
                (d / f"{args.prefix}_qc_metrics.json").write_text(json.dumps(data, indent=4))
                print(f"JSON: {d / f'{args.prefix}_qc_metrics.json'}")
                run_or_print(cmds[1], ["chromap", "samtools"])
        else:
            print(f"Command: {cmd_tar}")
            for c in cmds:
                run_or_print(c, ["chromap"] + (["samtools"] if "--SAM" in c else [])) if not which("chromap") else print(f"Command: {c}")
        steps.append({"step": "chromap", "ran": ran, "commands": [cmd_tar] + cmds})
    # downstream QC on whatever outputs exist
    qc = d / "qc"
    base = ["--out-dir"]
    np_flag = ["--no-plots"] if args.no_plots else []
    rna_tsv = atac_tsv = None
    if rna_h5ad and Path(rna_h5ad).exists():
        _run(["rna-qc-metrics", "--h5ad", str(rna_h5ad), "--kb-workflow", args.kb_mode, "--subpool",
              "none", *base, str(qc / "rna")])
        rna_tsv = qc / "rna" / "rna_barcode_metadata.tsv"
        _run(["rna-qc-plots", "--metrics", str(rna_tsv), *np_flag, *base, str(qc / "rna_plots")])
        steps.append({"step": "qc_rna", "metrics": str(rna_tsv)})
    if fragments and Path(fragments).exists() and args.tss_bed:
        _run(["tss-enrichment", "--fragments", str(fragments), "--regions", args.tss_bed, "--prefix", args.prefix, *np_flag,
              *base, str(qc / "atac")])
        atac_tsv = qc / "atac" / f"{args.prefix}.tss_enrichment_barcode_stats.tsv"
        _run(["atac-qc-plots", "--metrics", str(atac_tsv), *np_flag, *base, str(qc / "atac_plots")])
        steps.append({"step": "qc_atac", "metrics": str(atac_tsv)})
    if rna_tsv and atac_tsv:
        _run(["joint-qc", "--rna-metrics", str(rna_tsv), "--atac-metrics", str(atac_tsv), "--pkr", args.prefix,
              "--min-umis", str(args.min_umis), "--min-genes", str(args.min_genes), "--min-tss", str(args.min_tss),
              "--min-frags", str(args.min_frags), *np_flag, *base, str(qc / "joint")])
        steps.append({"step": "joint_qc", "table": str(qc / "joint" / "joint_barcode_metadata.csv")})
    imgs = sorted(str(p) for p in qc.rglob("*.png")) if qc.exists() else []
    if imgs or rna_tsv or atac_tsv:
        mfile = d / "qc_metrics.csv"
        mrows = []
        for s in (qc / "atac" / "summary.json", qc / "rna" / "summary.json", qc / "joint" / "summary.json"):
            if s.exists():
                for k, v in json.loads(s.read_text()).items():
                    if isinstance(v, (int, float)) and not isinstance(v, bool):
                        mrows.append(f"{s.parent.name}_{k},{v}")
        mfile.write_text("".join(r + "\n" for r in mrows))
        write_html_report(d / f"{args.prefix}.html", imgs, [], d / f"{args.prefix}.csv", str(mfile), None)
        print(f"Wrote: {d / f'{args.prefix}.html'}")
    plan = d / "plan.json"
    plan.write_text(json.dumps(_jsonable(steps), indent=2))
    print(f"JSON: {plan}")
    body = "\n".join(f"- **{s['step']}**" + (f" (ran: {s['ran']})" if "ran" in s else "") for s in steps)
    return finish(d, f"IGVF single-cell pipeline ({args.prefix})",
                  {"prefix": args.prefix, "process_atac": process_atac, "process_rna": process_rna,
                   "n_steps": len(steps), "steps": steps}, body)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _write_fastq(path, recs) -> None:
    with open_text(path, "w") as w:
        for name, seq in recs:
            w.write(f"@{name}\n{seq}\n+\n{'F' * len(seq)}\n")


def _rand_seq(rng, n: int) -> str:
    return "".join(rng.choice(list("ACGT"), n))


def _distinct_barcodes(rng, n: int, k: int = 8, min_dist: int = 3) -> "List[str]":
    out: List[str] = []
    while len(out) < n:
        s = _rand_seq(rng, k)
        if "GGGG" in s:
            continue
        if all(sum(a != b for a, b in zip(s, o)) >= min_dist for o in out):
            out.append(s)
    return out


def cmd_selftest(args) -> int:
    import tempfile
    np, pd = _np(), _pd()
    rng = np.random.default_rng(11)
    checks: List[Tuple[bool, str]] = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    npf = ["--no-plots"] if args.no_plots else []
    with tempfile.TemporaryDirectory() as td:
        T = Path(td)

        def od(name):
            return ["--out-dir", str(T / "runs" / name)]

        def summ(name):
            return json.loads((T / "runs" / name / "summary.json").read_text())

        # ---- barcode orientation
        print("\nbarcode-revcomp-detect")
        onl = _distinct_barcodes(rng, 40, 16, 4)
        (T / "onlist.txt").write_text("\n".join(onl) + "\n")
        recs = [(f"r{i}", _rand_seq(rng, 8) + (reverse_complement(onl[i % 40]) if i < 160 else _rand_seq(rng, 16))) for i in range(200)]
        _write_fastq(T / "bc.fq.gz", recs)
        rc = _run(["barcode-revcomp-detect", "--fastq", str(T / "bc.fq.gz"), "--onlist", str(T / "onlist.txt"), "--offset", "8",
                   "--num-reads", "200", *od("rc")])
        s = summ("rc")
        check(rc == 0 and s["read_format"] == "bc:8:-1:-" and abs(s["revcomp_match_proportion"] - 0.8) < 1e-9,
              f"reverse-complement barcodes detected: {s['read_format']}, rc proportion {s['revcomp_match_proportion']}")
        rc = _run(["barcode-revcomp-detect", "--fastq", str(T / "bc.fq.gz"), "--onlist", str(T / "onlist.txt"), "--offset", "8", *od("rc2")])
        check(rc == 1 and summ("rc2")["valid"] is False, "default 100 000-read denominator: 160 matches in a 200-read file is 'insufficient' (upstream behaviour)")
        _write_fastq(T / "bcf.fq.gz", [(f"r{i}", onl[i % 40]) for i in range(100)])
        _run(["barcode-revcomp-detect", "--fastq", str(T / "bcf.fq.gz"), "--onlist", str(T / "onlist.txt"), "--offset", "0",
              "--num-reads", "100", *od("rc3")])
        check(summ("rc3")["read_format"] == "bc:0:-1", "forward barcodes -> bc:0:-1")

        # ---- SHARE-seq correction
        print("\ncorrect-fastq")
        B = _distinct_barcodes(rng, 24, 8, 3)
        R1b, R2b, R3b = B[:8], B[8:16], B[16:24]
        (T / "wl.txt").write_text("".join(f"{a}{b}{c}\n" for a in R1b for b in R2b for c in R3b))

        def bregion(b1, b2, b3, shift=(0, 0, 0)):
            reg = list(_rand_seq(rng, 99))
            for bc, st, sh in ((b1, 15, shift[0]), (b2, 53, shift[1]), (b3, 91, shift[2])):
                for j, ch in enumerate(bc):
                    if 0 <= st + sh + j < 99:
                        reg[st + sh + j] = ch
            return "".join(reg)
        r1recs, r2recs, expect = [], [], Counter()
        for i in range(60):
            kind = ["exact", "mm", "lshift", "rshift", "polyg10", "polygbc", "junk"][i % 7] if i < 56 else "swap"
            b1, b2, b3 = R1b[i % 8], R2b[(i // 8) % 8], R3b[i % 8]
            umi = _rand_seq(rng, 10).replace("GGGGGGGGGG", "A")
            if kind == "exact":
                reg = bregion(b1, b2, b3)
            elif kind == "mm":
                reg = bregion(b1[:3] + ("A" if b1[3] != "A" else "C") + b1[4:], b2, b3)
            elif kind == "lshift":
                reg = bregion(b1, b2, b3, (-1, -1, -1))
            elif kind == "rshift":
                reg = bregion(b1, b2, b3, (1, 1, 1))
            elif kind == "polyg10":
                reg, umi = bregion(b1, b2, b3), "G" * 10
            elif kind == "polygbc":
                reg = "G" * 99
            elif kind == "junk":
                reg = bregion("TTTTTTTT", "TTTTTTTT", "TTTTTTTT")
            else:  # round-2 window holds a round-1 barcode: only the upstream shared set accepts it
                reg = bregion(b1, R1b[(i + 1) % 8], b3)
            if kind in ("exact", "mm", "lshift", "rshift"):
                expect["match"] += 1
            elif kind == "polyg10":
                expect["poly_G_in_first_10bp"] += 1
            elif kind == "polygbc":
                expect["poly_G_barcode"] += 1
            else:
                expect["mismatch"] += 1
            r1recs.append((f"read{i} extra", _rand_seq(rng, 50)))
            r2recs.append((f"read{i} extra", umi + _rand_seq(rng, 40) + reg))
        _write_fastq(T / "s_R1.fq.gz", r1recs)
        _write_fastq(T / "s_R2.fq.gz", r2recs)
        _run(["correct-fastq", "--read1", str(T / "s_R1.fq.gz"), "--read2", str(T / "s_R2.fq.gz"), "--whitelist", str(T / "wl.txt"),
              "--sample-type", "ATAC", "--prefix", "libA", "--pkr", "PKR1", *od("cf")])
        s = summ("cf")
        check(all(s[k] == expect[k] for k in ("match", "mismatch", "poly_G_barcode", "poly_G_in_first_10bp")),
              f"planted counts recovered: match {s['match']}/{expect['match']}, mismatch {s['mismatch']}, polyG bc {s['poly_G_barcode']}, polyG10 {s['poly_G_in_first_10bp']}")
        qc_txt = (T / "runs" / "cf" / "libA_barcode_qc.txt").read_text().splitlines()
        check(qc_txt[0] == "library\tmatch\tmismatch\tpoly_G_barcode\tpoly_G_in_first_10bp" and qc_txt[1].startswith("libA\t"),
              "QC file schema identical to upstream")
        heads = [r for r in read_fastq(T / "runs" / "cf" / "libA.barcode.corrected.fastq.gz")]
        first = heads[0]
        check(first[0] == f"@read0_{R1b[0]},{R2b[0]},{R3b[0]},PKR1" and first[1] == R1b[0] + R2b[0] + R3b[0],
              "ATAC header '<name>_<r1>,<r2>,<r3>,<pkr>' and corrected barcode read")
        mm = heads[1]
        check(mm[1][:8] == R1b[1], "1-mismatch round-1 barcode corrected to the onlist barcode")
        r2o = next(read_fastq(T / "runs" / "cf" / "libA.R2.corrected.fastq.gz"))
        check(len(r2o[1]) == 50, "ATAC read 2 loses its 99-bp barcode region")
        _run(["correct-fastq", "--read1", str(T / "s_R1.fq.gz"), "--read2", str(T / "s_R2.fq.gz"), "--whitelist", str(T / "wl.txt"),
              "--sample-type", "RNA", "--prefix", "libR", *od("cfr")])
        rr = next(read_fastq(T / "runs" / "cfr" / "libR.R2.corrected.fastq.gz"))
        check(rr[0].startswith(f"@read0_{R1b[0]},{R2b[0]},{R3b[0]}_") and len(rr[1]) == 34 and rr[1][:24] == R1b[0] + R2b[0] + R3b[0],
              "RNA read 2 = r1r2r3 + 10-bp UMI; UMI in header")
        _run(["correct-fastq", "--read1", str(T / "s_R1.fq.gz"), "--read2", str(T / "s_R2.fq.gz"), "--whitelist", str(T / "wl.txt"),
              "--sample-type", "ATAC", "--prefix", "libC", "--upstream-compat", *od("cfc")])
        check(summ("cfc")["match"] == expect["match"] + 4, "--upstream-compat: shared r1/r2/r3 set accepts round-1 barcodes in the round-2 window (+4)")

        # ---- dovetail trim
        print("\ntrim-fastq")
        t1, t2 = [], []
        for i in range(30):
            ins = _rand_seq(rng, 40 if i < 20 else 120)
            ad1, ad2 = "CTGTCTCTTATACACATCT", "CTGTCTCTTATACACATCTG"
            s1 = (ins + ad1 + _rand_seq(rng, 40))[:76]
            s2 = (reverse_complement(ins) + ad2 + _rand_seq(rng, 40))[:76]
            if 10 <= i < 20:  # a mismatch in read 1 inside the overlap -> fuzzy path
                p = 30
                s1 = s1[:p] + ("A" if s1[p] != "A" else "C") + s1[p + 1:]
            t1.append((f"t{i}", s1))
            t2.append((f"t{i}", s2))
        _write_fastq(T / "t_R1.fq", t1)
        _write_fastq(T / "t_R2.fq", t2)
        _run(["trim-fastq", "--read1", str(T / "t_R1.fq"), "--read2", str(T / "t_R2.fq"), *od("tr")])
        s = summ("tr")
        tr1 = [r[1] for r in read_fastq(T / "runs" / "tr" / "sample.R1.trimmed.fastq.gz")]
        check(s["total_reads"] == 60 and s["trimmed_reads"] == 40, f"20 overlapping pairs trimmed (10 exact + 10 via Levenshtein<=1), 10 long inserts kept: {s['trimmed_reads']}/{s['total_reads']}")
        check(all(len(x) == 40 for x in tr1[:20]) and all(len(x) == 76 for x in tr1[20:]), "reads cut to the insert length (rfind position + 20)")
        check(levenshtein_le1("ACGT", "ACGA") == 1 and levenshtein_le1("ACGT", "AGT") == 1 and levenshtein_le1("ACGT", "TGCA") == 2
              and levenshtein_le1("ACGTA", "ACG") == 2, "Levenshtein with cutoff 1 (substitution, deletion, >1)")

        # ---- TSS enrichment
        print("\ntss-enrichment")
        tss_rows, frag_rows = [], []
        chroms = {"chr1": 2_000_000, "chr2": 1_000_000}
        for c, L in chroms.items():
            for t in range(50_000, L - 50_000, 40_000):
                tss_rows.append((c, t, t + 1, "+" if (t // 40_000) % 2 else "-"))
        tss_df = pd.DataFrame(tss_rows)
        tss_df.to_csv(T / "tss.bed", sep="\t", header=False, index=False)
        cells = [f"CELL{i:03d}" for i in range(40)]
        noise = [f"NOISE{i:03d}" for i in range(40)]
        for bc in cells:
            for _ in range(400):
                c, t = tss_rows[rng.integers(len(tss_rows))][:2]
                st = int(t + rng.normal(0, 40))
                frag_rows.append((c, st, st + int(rng.integers(40, 200)), bc))
            for _ in range(300):
                c = "chr1" if rng.random() < 0.66 else "chr2"
                st = int(rng.integers(0, chroms[c] - 500))
                frag_rows.append((c, st, st + int(rng.integers(40, 300)), bc))
        for bc in noise:
            for _ in range(700):
                c = "chr1" if rng.random() < 0.66 else "chr2"
                st = int(rng.integers(0, chroms[c] - 500))
                frag_rows.append((c, st, st + int(rng.integers(40, 300)), bc))
        fr = pd.DataFrame(frag_rows).sort_values([0, 1])
        fr.to_csv(T / "frags.tsv", sep="\t", header=False, index=False)
        fr.to_csv(T / "frags_plain.tsv.gz", sep="\t", header=False, index=False, compression="gzip")
        _run(["tss-enrichment", "--fragments", str(T / "frags_plain.tsv.gz"), "--regions", str(T / "tss.bed"), "--prefix", "lib",
              "--engine", "pandas", *npf, *od("tss")])
        pb = pd.read_csv(T / "runs" / "tss" / "lib.tss_enrichment_barcode_stats.tsv", sep="\t")
        s = summ("tss")
        cm = pb[pb["barcode"].str.startswith("CELL")]["tss_enrichment"].median()
        nm = pb[pb["barcode"].str.startswith("NOISE")]["tss_enrichment"].median()
        check(list(pb.columns[:5]) == ["barcode", "reads_unique", "reads_promoter", "reads_tss", "tss_enrichment"],
              "per-barcode table: barcode, reads_unique, reads_promoter, reads_tss, tss_enrichment (joint-qc reads cols 1 and 4)")
        check(cm > 8 and nm < 3, f"planted TSS signal: cell median TSSe {cm:.1f} vs noise {nm:.2f}")
        check(s["tss_score_bulk"] > 3, f"bulk TSS score {s['tss_score_bulk']:.2f} > 3")
        check(int(pb.loc[pb["barcode"] == "CELL000", "reads_unique"].iloc[0]) == 1400, "reads_unique = 2 x fragments")
        try:
            import pysam  # type: ignore
            pysam.tabix_compress(str(T / "frags.tsv"), str(T / "frags.tsv.gz"), force=True)
            pysam.tabix_index(str(T / "frags.tsv.gz"), preset="bed", zerobased=True, force=True)
            _run(["tss-enrichment", "--fragments", str(T / "frags.tsv.gz"), "--regions", str(T / "tss.bed"), "--prefix", "lib",
                  "--engine", "tabix", "--no-plots", *od("tss_tb")])
            pb2 = pd.read_csv(T / "runs" / "tss_tb" / "lib.tss_enrichment_barcode_stats.tsv", sep="\t")
            m = pb.merge(pb2, on="barcode")
            check(len(m) == len(pb) and (m["reads_tss_x"] == m["reads_tss_y"]).all() and (m["reads_promoter_x"] == m["reads_promoter_y"]).all()
                  and abs(summ("tss_tb")["tss_score_bulk"] - s["tss_score_bulk"]) < 1e-9,
                  "tabix (pysam fetch) and pure-pandas engines give identical counts and bulk score")
        except ImportError:
            print("  skip  pysam not installed: tabix engine not compared")
        # hand-computed single fragment
        pd.DataFrame([("chr1", 10_000, 10_001, "+"), ("chr1", 30_000, 30_001, "-")]).to_csv(T / "tss1.bed", sep="\t", header=False, index=False)
        pd.DataFrame([("chr1", 9_990, 10_030, "A"), ("chr1", 29_000, 29_010, "B"), ("chr1", 31_950, 31_960, "B")]).to_csv(
            T / "fr1.tsv", sep="\t", header=False, index=False)
        r1 = tss_enrichment(T / "fr1.tsv", T / "tss1.bed", 2000, 20, 4, "pandas")
        a = r1["per_barcode"].set_index("barcode")
        check(a.loc["A", "reads_tss"] == 2 and abs(a.loc["A", "tss_enrichment"] - 2 * 2 / 101 / 0.2) < 1e-12,
              "ArchR TSSe = 2*reads_tss/101 / max(0.2, flank norm): 2 TSS insertions, no flank -> 0.198")
        check(a.loc["B", "reads_tss"] == 0 and a.loc["B", "reads_tss_flank"] == 2 and a.loc["B", "reads_promoter"] == 4,
              "'-' strand: +1950 bp insertions land in the (mirrored) left flank bin; 1000 bp ones only in the promoter")
        prof = r1["bulk_profile"]
        check(prof.sum() == 4000 + 6 and prof[1990] == 2 and prof[2030] == 2 and prof[3999 - 1000] == 2, "bulk profile = pseudocount 1 + insertions at promoter offsets")

        # ---- RNA metrics / mtx / suffix
        print("\nrna-qc-metrics / mtx-to-h5ad / modify-barcode-h5")
        import anndata  # type: ignore
        import scipy.sparse as sp  # type: ignore
        n_cells, n_bg, n_genes = 300, 2000, 400
        umis = np.r_[rng.lognormal(np.log(4000), 0.3, n_cells), rng.lognormal(np.log(30), 0.5, n_bg)]
        P = rng.dirichlet(np.ones(n_genes))
        Xtot = np.vstack([rng.multinomial(int(u), P) for u in umis])
        mat = rng.binomial(Xtot, 0.6)
        rest = Xtot - mat
        amb = rng.binomial(rest, 0.3)
        nas = rest - amb
        bcs = [f"CELL{i:04d}" for i in range(n_cells)] + [f"BG{i:05d}" for i in range(n_bg)]
        ad = anndata.AnnData(X=sp.csr_matrix(Xtot.astype("float32")),
                             layers={"mature": sp.csr_matrix(mat.astype("float32")), "ambiguous": sp.csr_matrix(amb.astype("float32")),
                                     "nascent": sp.csr_matrix(nas.astype("float32"))})
        ad.obs_names = bcs
        ad.obs_names.name = "barcode"
        ad.var_names = [f"G{j}" for j in range(n_genes)]
        ad.write_h5ad(T / "rna.h5ad")
        _run(["rna-qc-metrics", "--h5ad", str(T / "rna.h5ad"), "--kb-workflow", "nac", "--subpool", "SP1", *od("rq")])
        rq = pd.read_csv(T / "runs" / "rq" / "rna_barcode_metadata.tsv", sep="\t", index_col=0)
        g_exp = (mat > 0).sum(1) + (amb > 0).sum(1) + (nas > 0).sum(1)
        check(rq.index[0] == "CELL0000_SP1" and list(rq.columns) == ["total_counts", "genes"] and rq.index.name == "barcode",
              "columns total_counts, genes; barcodes suffixed _SP1")
        check(np.allclose(rq["total_counts"].to_numpy(), Xtot.sum(1)) and (rq["genes"].to_numpy() == g_exp).all(),
              "nac: layer sums of total_counts and n_genes_by_counts (genes summed across layers as upstream)")
        try:
            import scanpy as sc  # type: ignore
            qm = sc.pp.calculate_qc_metrics(ad, layer="mature", percent_top=None)[0]
            mine = rna_qc_metrics(anndata.AnnData(X=ad.layers["mature"], obs=ad.obs), "standard", "none")
            check(np.allclose(qm["total_counts"], mine["total_counts"]) and (qm["n_genes_by_counts"].to_numpy() == mine["genes"].to_numpy()).all(),
                  "matches scanpy.pp.calculate_qc_metrics exactly")
        except ImportError:
            print("  skip  scanpy not installed")
        import scipy.io  # type: ignore
        scipy.io.mmwrite(str(T / "m.mtx"), sp.csr_matrix(Xtot[:50]))
        (T / "b.txt").write_text("\n".join(bcs[:50]) + "\n")
        (T / "g.txt").write_text("\n".join(ad.var_names) + "\n")
        _run(["mtx-to-h5ad", "--mtx", str(T / "m.mtx"), "--barcodes", str(T / "b.txt"), "--genes", str(T / "g.txt"), *od("mtx")])
        a2 = anndata.read_h5ad(T / "runs" / "mtx" / "output.h5ad")
        check(a2.shape == (50, n_genes) and a2.obs_names[3] == bcs[3] and float(a2.X.sum()) == float(Xtot[:50].sum()), "mtx -> h5ad shape, names, counts")
        shutil.copy(T / "rna.h5ad", T / "rna_sfx.h5ad")
        _run(["modify-barcode-h5", "--h5ad", str(T / "rna_sfx.h5ad"), "--suffix", "SPX", *od("mod")])
        a3 = anndata.read_h5ad(T / "rna_sfx.h5ad")
        check(a3.obs_names[0] == "CELL0000_SPX" and a3.n_obs == n_cells + n_bg, "/obs/barcode suffixed in place; h5ad still readable")

        # ---- barcode rank / RNA + ATAC plots
        print("\nbarcode-rank / rna-qc-plots / atac-qc-plots")
        rq_plain = T / "rna_metrics.tsv"
        rna_qc_metrics(ad, "standard", "none").to_csv(rq_plain, sep="\t")
        _run(["rna-qc-plots", "--metrics", str(rq_plain), *npf, *od("rp")])
        s = summ("rp")
        check(n_cells <= s["umi_elbow_rank"] <= 3 * n_cells, f"UMI elbow (farthest point from the end-to-end line) at rank {s['umi_elbow_rank']}, past the {n_cells} planted cells")
        check("umi_knee_rank" in s and abs(s["umi_knee_rank"] - n_cells) <= 0.1 * n_cells,
              f"knee on ranks 1..elbow recovers the planted cells: rank {s.get('umi_knee_rank')} (planted {n_cells})")
        xx = np.arange(1, 201, dtype=float)
        check(elbow_knee_finder(xx, np.r_[np.full(100, 3.0), np.full(100, 1.0)])[0] in (100.0, 101.0),
              "line-distance elbow on a step curve sits at the step")
        sd1 = smooth_spline_deriv2(xx, 3 + 0.5 * xx)
        sd2 = smooth_spline_deriv2(xx, (xx / 200) ** 2)
        check(np.all(np.abs(sd1) < 1e-8) and np.all(sd2[20:180] > 0), "smooth.spline re-implementation: f''=0 for a line, > 0 for a convex quadratic")
        sd3 = smooth_spline_deriv2(xx, np.sin(xx / 30.0) + rng.normal(0, 0.01, 200))
        check(np.corrcoef(sd3[20:180], -np.sin(xx[20:180] / 30.0))[0, 1] > 0.9, "spar=1 smoothing keeps the sign pattern of a curved signal")
        check(_r_nknots(49) == 49 and _r_nknots(200) == 100 and _r_nknots(800) == 140 and _r_nknots(3200) == 200,
              ".nknots.smspl breakpoints (49, 100, 140, 200)")
        uq = np.r_[rng.lognormal(np.log(8000), 0.3, 250), rng.lognormal(np.log(40), 0.5, 3000)].astype(int)
        pd.DataFrame({"barcode": [f"A{i}" for i in range(len(uq))], "unique": uq}).to_csv(T / "atac_meta.tsv", sep="\t", index=False)
        _run(["atac-qc-plots", "--metrics", str(T / "atac_meta.tsv"), *npf, *od("ap")])
        check(abs(summ("ap").get("knee_rank", 0) - 250) <= 25, f"ATAC fragment barcode-rank knee at {summ('ap').get('knee_rank')} (planted 250 cells)")
        _run(["barcode-rank", "--metrics", str(rq_plain), "--column", "total_counts", "--no-plots", *od("br")])
        check(summ("br")["elbow_rank"] == s["umi_elbow_rank"], "generic barcode-rank matches the RNA plot elbow")

        # ---- joint QC
        print("\njoint-qc")
        pd.DataFrame({"barcode": ["b1", "b2", "b3", "b4", "b5"], "total_counts": [5000, 5000, 20, 30, 9],
                      "genes": [1000, 1000, 10, 15, 5]}).to_csv(T / "jr.tsv", sep="\t", index=False)
        pd.DataFrame({"barcode": ["b1", "b2", "b3", "b6", "b7"], "reads_unique": [4000, 100, 4000, 18, 4000], "reads_promoter": 0,
                      "reads_tss": 0, "tss_enrichment": [10.0, 10.0, 12.0, 1.0, 2.0]}).to_csv(T / "ja.tsv", sep="\t", index=False)
        _run(["joint-qc", "--rna-metrics", str(T / "jr.tsv"), "--atac-metrics", str(T / "ja.tsv"), "--pkr", "SP1", *npf, *od("jq")])
        jq = pd.read_csv(T / "runs" / "jq" / "joint_barcode_metadata.csv", index_col=0)
        exp = {"b1": "both", "b2": "RNA only", "b3": "ATAC only", "b4": "neither", "b7": "neither"}
        check(list(jq.columns) == ["umis", "genes", "tss", "frags", "QC", "QC_count"] and all(jq.loc[k, "QC"] == v for k, v in exp.items()),
              "QC categories both / RNA only / ATAC only / neither (TSS>=4 & frags>=100; UMIs>=100 & genes>=200)")
        check("b5" not in jq.index and "b6" not in jq.index and jq.loc["b4", "frags"] == 0 and jq.loc["b1", "frags"] == 2000,
              "low-yield barcodes (<10 UMIs, <10 fragments = reads/2) removed; missing modality filled with 0")
        check(jq.loc["b4", "QC_count"] == "neither (2)", "QC_count legend label '<QC> (<n>)'")

        # ---- insert size / 10x / log_atac
        print("\ninsert-size-hist / tenx-barcode-map / log-atac")
        (T / "hist.txt").write_text("## htsjdk\n## METRICS\nMEDIAN\n1\n\n## HISTOGRAM\ninsert_size\tAll_Reads.fr_count\n"
                                    + "".join(f"{i}\t{int(1000 * math.exp(-((i - 150) / 40) ** 2))}\n" for i in range(20, 400)) + "\n")
        _run(["insert-size-hist", "--histogram", str(T / "hist.txt"), "--pkr", "SP1", *npf, *od("ih")])
        check(summ("ih")["mode_insert_size"] == 150 and summ("ih")["n_sizes"] == 380, "Picard histogram parsed after the 'insert_size' header")
        (T / "atac_onl.txt").write_text("AAACGGTT\nACGTACGA\n")
        (T / "rna_onl.txt").write_text("TTTTCCCC\nGGGGAAAA\n")
        _run(["tenx-barcode-map", "--atac-onlist", str(T / "atac_onl.txt"), "--rna-onlist", str(T / "rna_onl.txt"), *od("tx")])
        tx = (T / "runs" / "tx" / "barcode_conversion_dict.tsv").read_text().splitlines()
        check(tx == ["TTTTCCCC\tAACCGTTT", "GGGGAAAA\tTCGTACGT", "TTTTCCCC\tAAACGGTT", "GGGGAAAA\tACGTACGA"],
              "rows '<rna>\\t<revcomp(atac)>' then '<rna>\\t<atac>'")
        (T / "rna_onl3.txt").write_text("TTTTCCCC\n")
        check(_run(["tenx-barcode-map", "--atac-onlist", str(T / "atac_onl.txt"), "--rna-onlist", str(T / "rna_onl3.txt"), *od("tx2")]) == 1,
              "unequal onlists -> no dictionary")
        (T / "chromap.log").write_text("Number of threads: 8\nNumber of reads: 2000000.\nNumber of mapped reads: 1900000.\n"
                                       "Number of uniquely mapped reads: 1700000.\n#reads: 20, #mapped: 19.\nMapped all reads in 3.2s.\n")
        (T / "summary.csv").write_text("barcode,total,duplicate,unmapped,lowmapq\nA,100,20,5,15\nB,300,60,15,65\n")
        _run(["log-atac", "--alignment-log", str(T / "chromap.log"), "--barcode-summary", str(T / "summary.csv"), "--prefix", "lib", *od("la")])
        la = json.loads((T / "runs" / "la" / "lib_qc_metrics.json").read_text())
        check(la["Number_of_reads"] == 2000000 and la["reads"] == 20 and "Number_of_threads" not in la
              and la["percentage_duplicates"] == 26.7, f"chromap log -> json; percentage_duplicates 100*80/300 = {la['percentage_duplicates']}")

        # ---- wrappers
        print("\nkb-count / chromap-align / subpool-fragments")
        (T / "onl.txt.gz").write_bytes(gzip.compress(b"AAAA\nCCCC\n"))
        _run(["kb-count", "--read1", "L1_R1.fq.gz", "L2_R1.fq.gz", "--read2", "L1_R2.fq.gz", "L2_R2.fq.gz", "--index-dir", "idx",
              "--read-format", "1,0,16:1,16,28:0,0,0", "--onlist", str(T / "onl.txt.gz"), "--replacement-list", str(T / "onl.txt.gz"),
              "--strand", "forward", "--output-dir", str(T / "kbo"), *od("kb")])
        kc = summ("kb")["command"]
        check(("--workflow=nac" in kc and "-c1 idx/cdna.txt -c2 idx/nascent.txt --sum=total" in kc and " -r " in kc
               and kc.endswith("L1_R1.fq.gz L1_R2.fq.gz L2_R1.fq.gz L2_R2.fq.gz") and (T / "runs" / "kb" / "barcode_inclusion_list.txt").exists())
              or bool(which("kb")), "kb count command: nac index files, replacement list, interleaved lanes, gunzipped onlist")
        cu = T / "kbo" / "counts_unfiltered_modified"
        cu.mkdir(parents=True)
        shutil.copy(T / "rna.h5ad", cu / "adata.h5ad")
        (cu / "cells_x_genes.barcodes.txt").write_text("CELL0000\nCELL0001\n")
        post = kb_postprocess(str(T / "kbo"), "SP9", "repl.txt")
        check(post["h5ad"] == str(T / "kbo") + ".h5ad" and anndata.read_h5ad(post["h5ad"]).obs_names[0] == "CELL0000_SP9"
              and (cu / "cells_x_genes.barcodes.txt").read_text() == "CELL0000_SP9\nCELL0001_SP9\n",
              "kb post-processing: subpool on _modified h5ad + barcodes.txt, h5ad moved to <out>.h5ad")
        _run(["chromap-align", "--output", "bam", "--index-dir", "idx", "--read-format", "bc:0:15", "--reference-fasta", "hg38.fa",
              "--onlist", str(T / "onlist.txt"), "--read1", "a_R1.fq.gz", "b_R1.fq.gz", "--read2", "a_R2.fq.gz", "b_R2.fq.gz",
              "--read-barcode", "a_I2.fq.gz", "b_I2.fq.gz", "--subpool", "SP1", "--barcode-translate", str(T / "onlist.txt"), *od("cb")])
        cb = summ("cb")["command"]
        check("--SAM -l 2000 --bc-error-threshold 1" in cb and "-1 a_R1.fq.gz,b_R1.fq.gz" in cb and "--barcode-translate" in cb
              and "suffix=SP1" in cb and "samtools view -bS -" in cb, "chromap bam command: --SAM, comma-joined lanes, CB:Z subpool awk")
        shutil.copy(T / "frags.tsv", T / "frags_sp.tsv")
        (T / "summary_sp.csv").write_text((T / "summary.csv").read_text())
        _run(["subpool-fragments", "--fragments", str(T / "frags_sp.tsv"), "--summary", str(T / "summary_sp.csv"), "--subpool", "SP1",
              "--bgzip", *od("spf")])
        fl = (T / "frags_sp.tsv").read_text().splitlines()[0].split("\t")
        sl = (T / "summary_sp.csv").read_text().splitlines()
        check(fl[3].endswith("_SP1") and sl[0].startswith("barcode,") and sl[1].startswith("A_SP1,"),
              "subpool appended to fragment column 4 and summary column 1 (header kept)")
        check(summ("spf").get("fragments_gz") is None or Path(str(T / "frags_sp.tsv.gz.tbi")).exists(), "bgzip + tabix (zero-based bed)")

        # ---- inputs / sampling / portal / synapse / html
        print("\ngenome-tsv / check-inputs / sample-fastqs / portal-download / synapse / html-report")
        (T / "genome.tsv").write_text(f"fasta\t{T / 'hg38.fa'}\nkb_nac_idx_tar\tgs://b/kb.tar.gz\nchromap_idx_tar\tsyn123\n")
        _run(["genome-tsv", "--genome-tsv", str(T / "genome.tsv"), *od("gt")])
        check(summ("gt")["kb_index_tar_gz"] == "gs://b/kb.tar.gz" and summ("gt")["fasta"].endswith("hg38.fa"), "genome TSV keys fasta / kb_nac_idx_tar / chromap_idx_tar")
        rc = _run(["check-inputs", "--paths", str(T / "onlist.txt"), "gs://bucket/x.fq.gz", "syn999",
                   "https://api.data.igvf.org/sequence-files/IGVFFI0000AAAA/@@download/IGVFFI0000AAAA.fastq.gz", *od("ci")])
        ci = pd.read_csv(T / "runs" / "ci" / "inputs.tsv", sep="\t")
        check(rc == 0 and list(ci["kind"]) == ["local", "gcs", "synapse", "https"] and ci["action"].iloc[3] == "dry-run",
              "inputs classified local / gs:// kept / synapse get / portal https (dry run)")
        check(_run(["check-inputs", "--paths", str(T / "nope.fq"), *od("ci2")]) == 1, "missing local input -> rc 1")
        _run(["sample-fastqs", "--fastqs", str(T / "s_R1.fq.gz"), "--n-reads", "7", *od("sf")])
        check(sum(1 for _ in read_fastq(T / "runs" / "sf" / "s_R1.fq.gz")) == 7, "first N reads kept, re-gzipped with the same basename")
        _run(["portal-download", "--urls", "https://api.data.igvf.org/x/@@download/f.txt", "--dry-run", *od("pdl")])
        check(summ("pdl")["files"][0]["status"] == "dry-run", "portal download dry run (no network, env credentials only)")
        res_tab = pd.DataFrame({"Subpool": ["SP1", "SP2"], "ATAC_barcode": ["gs://b/a_I2.fq.gz"] * 2, "ATAC_fastq_R1": ["gs://b/a_R1.fq.gz"] * 2,
                                "ATAC_fastq_R2": ["gs://b/a_R2.fq.gz"] * 2, "RNA_fastq_R1": ["gs://b/r_R1.fq.gz"] * 2,
                                "RNA_fastq_R2": ["gs://b/r_R2.fq.gz"] * 2, "seqspec": ["[\"gs://b/s.yaml\"]"] * 2,
                                "seqspec_atac_onlist_renamed": ["gs://b/SP1.atac.onlist.txt.gz", "gs://b/SP2.atac.onlist.txt.gz"],
                                "atac_bam": ["gs://b/SP1.atac.align.k4.hg38.bam", float("nan")],
                                "rna_barcode_metadata": ["gs://b/SP1.qc.rna.hg38.barcode.metadata.tsv"] * 2,
                                "seqspec_rna_onlist_renamed": ["gs://b/SP1.rna.onlist.txt.gz"] * 2, "unrelated": ["x", "y"]})
        res_tab.to_csv(T / "results.tsv", sep="\t", index=False)
        _run(["synapse-manifest", "--table", str(T / "results.tsv"), "--project", "syn000", "--local-root", "/mnt/",
              "--remote-root", "gs://", *od("sm")])
        man = pd.read_csv(T / "runs" / "sm" / "manifest.tsv", sep="\t", quoting=3)
        bam = man[man["name"] == "SP1.atac.align.k4.hg38.bam"].iloc[0]
        check(len(man) == 7 and bam["path"] == "/mnt/b/SP1.atac.align.k4.hg38.bam" and bam["parent"] == "syn000/SP1"
              and bam["used"].startswith('"/mnt/b/a_I2.fq.gz;') and "task_chromap_bam.wdl" in bam["executed"],
              "manifest rows (NaN skipped, unrelated column ignored), roots rewritten, provenance 'used' list")
        _run(["synapse-upload", "--manifest", str(T / "runs" / "sm" / "manifest.tsv"), *od("su")])
        check(summ("su")["n_batches"] == 2 and summ("su")["dry_run"], "batch upload grouped by parent (2 subpools), dry run")
        _run(["synapse-annotations", "--names", "SP1.atac.fragments.hg38.tsv.gz", "SP1.rna.align.kb.hg38.cells_x_genes.total.h5ad",
              "weird.bin", *od("sa")])
        an = pd.read_csv(T / "runs" / "sa" / "annotations.tsv", sep="\t")
        check(an.loc[0, "file_type"] == "tsv.gz" and an.loc[0, "data_type"] == "atac" and an.loc[1, "main_output"] == "yes"
              and an.loc[1, "file_description"].startswith("[Raw] Aggregated") and pd.isna(an.loc[2, "file_type"]),
              "annotations: type/status/main_output, data_type, description; unknown suffix -> None (upstream crashed)")
        plt = _plt()
        imgs = []
        if plt is not None:
            for nm in ("SP1.atac.tss.png", "SP1.rna.umi.png", "SP1.joint.png"):
                fig, ax = plt.subplots(figsize=(1, 1))
                fig.savefig(T / nm)
                plt.close(fig)
                imgs.append(str(T / nm))
        (T / "am.csv").write_text("tss_score,12.5\nfragments,2500000\n")
        _run(["html-report", "--images", *imgs, "--logs", "gs://b/x.log.txt", "--atac-metrics", str(T / "am.csv"), "--prefix", "SP1", *od("hr")])
        html = (T / "runs" / "hr" / "SP1.html").read_text()
        check(summ("hr")["n_atac"] == (1 if imgs else 0) and "2 M" in html and "data:image/png;base64" in html if imgs else "2 M" in html,
              "HTML report: images split by modality, metrics formatted (2 500 000 -> '2 M')")

        # ---- the WDL plan end to end
        print("\npipeline")
        inputs = {"single_cell_pipeline.prefix": "SP1", "single_cell_pipeline.subpool": "SP1",
                  "single_cell_pipeline.genome_tsv": str(T / "genome.tsv"), "single_cell_pipeline.create_onlist_mapping": True,
                  "single_cell_pipeline.atac_read1": ["gs://b/a_R1.fq.gz"], "single_cell_pipeline.atac_read2": ["gs://b/a_R2.fq.gz"],
                  "single_cell_pipeline.fastq_barcode": ["gs://b/a_I2.fq.gz"], "single_cell_pipeline.atac_barcode_inclusion_list": str(T / "atac_onl.txt"),
                  "single_cell_pipeline.atac_read_format": "bc:8:-1:-", "single_cell_pipeline.rna_read1": ["gs://b/r_R1.fq.gz"],
                  "single_cell_pipeline.rna_read2": ["gs://b/r_R2.fq.gz"], "single_cell_pipeline.rna_barcode_inclusion_list": str(T / "rna_onl.txt"),
                  "single_cell_pipeline.rna_read_format": "0,0,16:0,16,28:1,0,0", "single_cell_pipeline.kb_strand": "forward"}
        (T / "inputs.json").write_text(json.dumps(inputs))
        _write_fastq(T / "dummy.fq", [])
        _run(["pipeline", "--inputs-json", str(T / "inputs.json"), "--rna-h5ad", str(T / "rna.h5ad"), "--fragments", str(T / "frags_plain.tsv.gz"),
              "--tss-bed", str(T / "tss.bed"), "--min-frags", "100", "--min-tss", "4", *npf, *od("pl")])
        steps = {s["step"]: s for s in summ("pl")["steps"]}
        check({"genome_tsv", "check_inputs", "mapping_tenx_barcodes", "kb_count", "chromap", "qc_rna", "qc_atac", "joint_qc"} <= set(steps),
              "plan: genome_tsv, check_inputs, 10x mapping, kb count, chromap, RNA / ATAC / joint QC")
        check(any("--barcode-translate" in c for c in steps["chromap"]["commands"]) and "-r " not in steps["kb_count"]["commands"][1]
              and steps["chromap"]["commands"][1].count("chromap -x") == 1, "chromap gets the 10x conversion dict; no replacement list -> no -r")
        jt = pd.read_csv(steps["joint_qc"]["table"], index_col=0)
        check((jt["QC"] == "both").sum() == 0 and ((jt["QC"] == "RNA only").sum() >= n_cells * 0.9),
              "joint QC runs on pipeline outputs (synthetic RNA and ATAC barcodes are disjoint -> no 'both')")
        check((T / "runs" / "pl" / "SP1.html").exists(), "pipeline HTML summary written")

    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="igvfagent igvf-sc-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add(name, func, help_, label):
        p = sub.add_parser(name, help=help_)
        p.add_argument("--label", default=label)
        p.add_argument("--out-dir", help="write here instead of Docs/IGVFscPipeline/<timestamp>_<label>")
        p.set_defaults(func=func)
        return p

    p = add("barcode-revcomp-detect", cmd_barcode_revcomp_detect, "barcode read orientation vs onlist -> chromap read format", "barcode_revcomp")
    p.add_argument("--fastq", required=True)
    p.add_argument("--onlist", required=True)
    p.add_argument("--offset", type=int, default=0, help="barcode start in the read (10x = 0, 10x multiome ATAC = 8)")
    p.add_argument("--num-reads", type=int, default=100_000)
    p.add_argument("--threshold", type=float, default=0.45)

    p = add("correct-fastq", cmd_correct_fastq, "SHARE-seq barcode correction", "correct_fastq")
    p.add_argument("--read1", required=True)
    p.add_argument("--read2", required=True)
    p.add_argument("--whitelist", required=True, help="R1R2R3 24-mer combinations, one per line")
    p.add_argument("--sample-type", choices=["ATAC", "RNA"], required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--pkr")
    p.add_argument("--upstream-compat", action="store_true", help="shared r1/r2/r3 barcode set and mismatch-overwrites-exact dictionary")

    p = add("trim-fastq", cmd_trim_fastq, "dovetail trimming", "trim_fastq")
    p.add_argument("--read1", required=True)
    p.add_argument("--read2", required=True)
    p.add_argument("--prefix", default="sample")

    p = add("tss-enrichment", cmd_tss_enrichment, "bulk + per-barcode TSS enrichment", "tss")
    p.add_argument("--fragments", required=True)
    p.add_argument("--regions", required=True, help="TSS bed: chr, start (TSS), end, ..., strand")
    p.add_argument("--flank", type=int, default=2000)
    p.add_argument("--window", type=int, default=20)
    p.add_argument("--strand-col", type=int, default=4)
    p.add_argument("--prefix", default="sample")
    p.add_argument("--engine", choices=["auto", "tabix", "pandas"], default="auto")
    p.add_argument("--no-plots", action="store_true")

    p = add("snapatac2-tsse", cmd_snapatac2_tsse, "snapatac2 TSSe + FRiP (optional dependency)", "snapatac2")
    p.add_argument("--fragments", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--tss-bed", required=True)
    p.add_argument("--promoter-bed", required=True)
    p.add_argument("--min-frag-cutoff", type=int, default=100)
    p.add_argument("--prefix", default="sample")
    p.add_argument("--no-plots", action="store_true")

    p = add("rna-qc-metrics", cmd_rna_qc_metrics, "per-barcode total_counts / genes from a kb h5ad", "rna_qc_metrics")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--kb-workflow", choices=["standard", "nac"], default="standard")
    p.add_argument("--subpool", default="none")
    p.add_argument("--out-name")

    p = add("mtx-to-h5ad", cmd_mtx_to_h5ad, "mtx + barcodes + genes -> h5ad", "mtx_to_h5ad")
    p.add_argument("--mtx", required=True)
    p.add_argument("--barcodes", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--out")

    p = add("modify-barcode-h5", cmd_modify_barcode_h5, "append _<suffix> to /obs/barcode in place", "modify_barcode_h5")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--suffix", required=True)

    p = add("joint-qc", cmd_joint_qc, "joint RNA x ATAC barcode QC", "joint_qc")
    p.add_argument("--rna-metrics", required=True)
    p.add_argument("--atac-metrics", required=True)
    p.add_argument("--remove-low-yielding-cells", type=int, default=10)
    p.add_argument("--min-umis", type=int, default=100)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-tss", type=float, default=4)
    p.add_argument("--min-frags", type=int, default=100)
    p.add_argument("--pkr", default="")
    p.add_argument("--no-plots", action="store_true")

    p = add("barcode-rank", cmd_barcode_rank, "elbow / knee of a barcode-rank curve", "barcode_rank")
    p.add_argument("--metrics", required=True)
    p.add_argument("--column", default="1", help="column name or 0-based index (default 1)")
    p.add_argument("--cutoff", type=float, default=10)
    p.add_argument("--no-plots", action="store_true")

    p = add("atac-qc-plots", cmd_atac_qc_plots, "ATAC fragment barcode-rank plots", "atac_qc_plots")
    p.add_argument("--metrics", required=True)
    p.add_argument("--fragment-cutoff", type=float, default=10)
    p.add_argument("--column")
    p.add_argument("--no-plots", action="store_true")

    p = add("rna-qc-plots", cmd_rna_qc_plots, "RNA UMI / gene barcode-rank plots", "rna_qc_plots")
    p.add_argument("--metrics", required=True)
    p.add_argument("--umi-cutoff", type=float, default=10)
    p.add_argument("--gene-cutoff", type=float, default=10)
    p.add_argument("--no-plots", action="store_true")

    p = add("insert-size-hist", cmd_insert_size_hist, "Picard insert-size histogram plot", "insert_size")
    p.add_argument("--histogram", required=True)
    p.add_argument("--pkr", default="")
    p.add_argument("--no-plots", action="store_true")

    p = add("tenx-barcode-map", cmd_tenx_barcode_map, "10x multiome ATAC<->RNA conversion dictionary", "tenx_barcode_map")
    p.add_argument("--atac-onlist", required=True)
    p.add_argument("--rna-onlist", required=True)

    p = add("log-atac", cmd_log_atac, "chromap log + barcode summary -> qc json", "log_atac")
    p.add_argument("--alignment-log", required=True)
    p.add_argument("--barcode-summary", required=True)
    p.add_argument("--prefix", default="sample")

    p = add("kb-count", cmd_kb_count, "kallisto|bustools count (optional binary)", "kb_count")
    p.add_argument("--read1", nargs="+", required=True)
    p.add_argument("--read2", nargs="+", required=True)
    p.add_argument("--read-barcode", nargs="+", default=[])
    p.add_argument("--index-dir", required=True)
    p.add_argument("--read-format", required=True)
    p.add_argument("--onlist", required=True)
    p.add_argument("--kb-mode", choices=["nac", "standard"], default="nac")
    p.add_argument("--strand", default="unstranded")
    p.add_argument("--replacement-list")
    p.add_argument("--subpool", default="none")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--temp-dir")
    p.add_argument("--output-dir")

    p = add("kb-index", cmd_kb_index, "kallisto|bustools ref (optional binary)", "kb_index")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--kb-mode", choices=["nac", "standard"], default="nac")
    p.add_argument("--temp-dir")
    p.add_argument("--output-dir")

    p = add("chromap-align", cmd_chromap_align, "chromap fragments / bam (optional binary)", "chromap")
    p.add_argument("--output", choices=["fragments", "bam"], default="fragments")
    p.add_argument("--index-dir", required=True)
    p.add_argument("--read-format", required=True)
    p.add_argument("--reference-fasta", required=True)
    p.add_argument("--onlist", required=True)
    p.add_argument("--read1", nargs="+", required=True)
    p.add_argument("--read2", nargs="+", required=True)
    p.add_argument("--read-barcode", nargs="+", required=True)
    p.add_argument("--barcode-translate")
    p.add_argument("--subpool", default="none")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--prefix", default="sample.atac")

    p = add("chromap-index", cmd_chromap_index, "chromap index (optional binary)", "chromap_index")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--output-dir")

    p = add("subpool-fragments", cmd_subpool_fragments, "append _<subpool> to fragment / summary barcodes", "subpool")
    p.add_argument("--subpool", required=True)
    p.add_argument("--fragments")
    p.add_argument("--summary")
    p.add_argument("--bgzip", action="store_true", help="bgzip + tabix the fragments afterwards")

    p = add("genome-tsv", cmd_genome_tsv, "parse the genome TSV", "genome_tsv")
    p.add_argument("--genome-tsv", required=True)
    p.add_argument("--fasta")
    p.add_argument("--kb-index")
    p.add_argument("--chromap-index")

    p = add("check-inputs", cmd_check_inputs, "classify / fetch pipeline inputs", "check_inputs")
    p.add_argument("--paths", nargs="+", required=True)
    p.add_argument("--download", action="store_true", help="actually fetch syn* / https inputs")

    p = add("sample-fastqs", cmd_sample_fastqs, "first N reads of FASTQs", "sample_fastqs")
    p.add_argument("--fastqs", nargs="+", required=True)
    p.add_argument("--n-reads", type=int, default=10_000_000)

    p = add("portal-download", cmd_portal_download, "IGVF portal download with env credentials", "portal_download")
    p.add_argument("--urls", nargs="+", required=True)
    p.add_argument("--dest")
    p.add_argument("--dry-run", action="store_true")

    p = add("synapse-manifest", cmd_synapse_manifest, "Synapse upload manifest from a results table", "synapse_manifest")
    p.add_argument("--table", required=True)
    p.add_argument("--project", required=True)
    p.add_argument("--local-root", default="")
    p.add_argument("--remote-root", default="")
    p.add_argument("--output")
    p.add_argument("--execute", action="store_true", help="resolve folders / provenance on Synapse (SYNAPSE_AUTH_TOKEN)")

    p = add("synapse-upload", cmd_synapse_upload, "batch Synapse upload of a manifest", "synapse_upload")
    p.add_argument("--manifest", required=True)
    p.add_argument("--execute", action="store_true")

    p = add("synapse-annotations", cmd_synapse_annotations, "file annotations for Synapse", "synapse_annotations")
    p.add_argument("--names", nargs="+")
    p.add_argument("--names-file")
    p.add_argument("--root-folder")
    p.add_argument("--execute", action="store_true")

    p = add("html-report", cmd_html_report, "tabbed HTML + CSV summary", "html_report")
    p.add_argument("--images", nargs="+", default=[])
    p.add_argument("--image-list")
    p.add_argument("--logs", nargs="+", default=[])
    p.add_argument("--log-list")
    p.add_argument("--atac-metrics")
    p.add_argument("--rna-metrics")
    p.add_argument("--prefix", default="summary")

    p = add("pipeline", cmd_pipeline, "single_cell_pipeline.wdl plan + downstream QC", "pipeline")
    p.add_argument("--inputs-json", help="Cromwell inputs JSON (single_cell_pipeline.<key>)")
    p.add_argument("--prefix", help="analysis-set prefix (default sample)")
    p.add_argument("--subpool", help="subpool suffix (default none)")
    p.add_argument("--genome-tsv")
    p.add_argument("--genome-fasta")
    p.add_argument("--kb-index")
    p.add_argument("--chromap-index")
    p.add_argument("--create-onlist-mapping", action="store_true")
    p.add_argument("--atac-read1", nargs="+")
    p.add_argument("--atac-read2", nargs="+")
    p.add_argument("--fastq-barcode", nargs="+")
    p.add_argument("--atac-onlist")
    p.add_argument("--atac-read-format")
    p.add_argument("--rna-read1", nargs="+")
    p.add_argument("--rna-read2", nargs="+")
    p.add_argument("--rna-barcode", nargs="+")
    p.add_argument("--rna-onlist")
    p.add_argument("--rna-read-format")
    p.add_argument("--kb-mode", choices=["nac", "standard"], help="default nac")
    p.add_argument("--kb-strand", help="default unstranded")
    p.add_argument("--replacement-list")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--download", action="store_true")
    p.add_argument("--rna-h5ad", help="existing kb h5ad for the QC stages")
    p.add_argument("--fragments", help="existing fragment file for the QC stages")
    p.add_argument("--tss-bed", help="TSS bed for tss-enrichment")
    p.add_argument("--min-umis", type=int, default=100)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-tss", type=float, default=4)
    p.add_argument("--min-frags", type=int, default=100)
    p.add_argument("--no-plots", action="store_true")

    p = sub.add_parser("selftest", help="synthetic data with planted signal; every subcommand asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)
    return ap


def _run(argv: "List[str]") -> int:
    """Dispatch without re-initialising logging (used by pipeline + selftest)."""
    args = build_parser().parse_args(argv)
    return args.func(args)


def main(argv: "Optional[List[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
