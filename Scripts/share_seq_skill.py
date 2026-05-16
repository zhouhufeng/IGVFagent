#!/usr/bin/env python3
"""IGVF agent SHARE-seq analytics skill.

End-to-end QC and joint-cell analytics for **SHARE-seq** datasets
(simultaneous single-cell ATAC + single-cell RNA profiling, Ma 2020,
Cell). The math is a clean-room reimplementation of the QC algorithms
in [broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline)
(MIT, 2021), following the published descriptions only.

Scope: this skill is designed to consume the **processed deposits** that
IGVF Portal SHARE-seq AnalysisSets actually publish:

  * ``fragments.bed[.gz]``  — ATAC fragments file (4 cols + count)
  * ``h5ad``                 — sparse RNA gene-count matrix per barcode
  * ``bam`` (optional)       — aligned reads (only used for advanced QC)

It does NOT re-run FASTQ → BAM alignment (IGVF Portal already publishes
those outputs). What it does is:

  pull-portal         Discover IGVF Portal SHARE-seq AnalysisSets/files.
  demultiplex-bcs     Clean-room SHARE-seq round-1/2/3 barcode demultiplex
                      (1-mismatch + ±1bp shift correction, R1+R2+R3 24mer).
  fragment-qc         Per-barcode fragment counts, TSS enrichment,
                      reads-in-peaks (FRIP).
  rna-qc              Per-barcode UMI / gene / mitochondrial-fraction
                      from an h5ad file.
  joint-qc            Merge ATAC + RNA per-barcode tables and apply the
                      Ma-2020 joint-cell pass thresholds.
  multiplet-detect    Pairwise Jaccard multiplet detection on fragments.
  write-playbook      Emit the skill's markdown playbook.

License: Apache-2.0. The upstream pipeline is MIT (compatible). No source
is copied verbatim — the algorithms are paraphrased from the publicly
described pipeline + supplementary methods of Ma et al. 2020.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import logging
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "SHAREseq"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

PORTAL_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_BASE")


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"share_seq_skill_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(
            f"Missing dependency '{name}'. {hint}\nInstall with: pip install {name}"
        ) from exc


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Portal discovery (real IGVF SHARE-seq AnalysisSets)
# ---------------------------------------------------------------------------

def _portal_get(path: str, *, params: dict[str, Any] | None = None) -> Any:
    qs = "?" + urllib.parse.urlencode(params) if params else ""
    url = f"{PORTAL_BASE}{path}{qs}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def cmd_pull_portal(args: argparse.Namespace) -> int:
    """Discover IGVF Portal SHARE-seq AnalysisSets (preferred — these carry
    processed fragments + h5ad) and MeasurementSets (raw FASTQ + seqspec)."""
    setup_logging()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")

    rows = []
    for set_type in ("AnalysisSet", "MeasurementSet"):
        data = _portal_get("/search/", params={
            "type": set_type,
            "preferred_assay_titles": "SHARE-seq",
            "format": "json",
            "limit": str(args.limit),
        })
        for r in data.get("@graph", []):
            files = r.get("files") or []
            rows.append({
                "set_type": set_type,
                "accession": r.get("accession", ""),
                "assay_titles": "|".join(r.get("assay_titles") or []),
                "preferred_assay_titles": "|".join(r.get("preferred_assay_titles") or []),
                "n_files": len(files),
                "status": r.get("status", ""),
                "description": (r.get("description") or "")[:120],
            })
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_portal_share.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                            ["set_type", "accession", "assay_titles", "preferred_assay_titles",
                             "n_files", "status", "description"],
                            delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"SHARE-seq portal manifest: {out} ({len(rows)} sets)")
    save_json(f"share_{args.label}_portal", {"n": len(rows), "out": str(out)})
    return 0


# ---------------------------------------------------------------------------
# Barcode demultiplex (clean-room: round1+round2+round3 + 1mm + ±1 bp shift)
# ---------------------------------------------------------------------------

def _hamming1_neighbors(seq: str) -> Iterable[str]:
    bases = "ACGTN"
    for i, b in enumerate(seq):
        for c in bases:
            if c != b:
                yield seq[:i] + c + seq[i + 1:]


def _load_whitelist(path: Path, *, length: int = 24) -> tuple[dict[str, str], dict[str, str]]:
    """Load a SHARE-seq whitelist file.

    Format: one concatenated R1+R2+R3 24-mer per line (alternatively, 3
    columns of 8-mers). Returns (exact_map, hamming1_map) where keys map
    to the **canonical** 24-mer.
    """
    exact: dict[str, str] = {}
    hamming1: dict[str, str] = {}
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.replace(",", "\t").split("\t") if "\t" in line or "," in line else [line]
            canonical = "".join(parts)[:length]
            if len(canonical) != length:
                continue
            exact[canonical] = canonical
            for nb in _hamming1_neighbors(canonical):
                hamming1.setdefault(nb, canonical)
    return exact, hamming1


def _extract_bc_24mer(seq: str, *, offsets: tuple[int, int, int]) -> str | None:
    """Pull round-1/2/3 8-mers from a SHARE-seq R2 read at fixed offsets."""
    r1, r2, r3 = offsets
    if len(seq) < r3 + 8:
        return None
    return seq[r1:r1 + 8] + seq[r2:r2 + 8] + seq[r3:r3 + 8]


def _correct_24mer(observed: str, exact: dict[str, str], hamming1: dict[str, str]) -> tuple[str | None, str]:
    """Returns (canonical_bc | None, status). Status ∈ {exact, hamming1, miss}."""
    if observed in exact:
        return exact[observed], "exact"
    if observed in hamming1:
        return hamming1[observed], "hamming1"
    return None, "miss"


def cmd_demultiplex_bcs(args: argparse.Namespace) -> int:
    """SHARE-seq barcode demultiplex on a single R2 FASTQ.

    Reads a (gz-compressed) FASTQ, extracts the round-1/2/3 24-mer at
    positions ``r1_offset/r2_offset/r3_offset``, attempts 1-Hamming-mismatch
    correction (also tries ±1 bp positional shift), tags the read header
    with the canonical barcode, and emits a corrected FASTQ + a per-read
    barcode log. Pure stdlib + numpy — no pysam/dnaio needed.
    """
    np = _require_pkg("numpy", "")
    setup_logging()
    exact, hamming1 = _load_whitelist(Path(args.whitelist))
    logging.info("Loaded whitelist: %d exact + %d 1-Hamming neighbors",
                  len(exact), len(hamming1))
    in_fastq = Path(args.fastq)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    counter = collections.Counter()

    open_in = gzip.open if str(in_fastq).endswith(".gz") else open
    open_out = gzip.open if str(out_path).endswith(".gz") else open

    offsets_main = (args.r1_offset, args.r2_offset, args.r3_offset)
    shifts = [0, 1, -1] if args.shift_correct else [0]

    with open_in(in_fastq, "rt") as fin, open_out(out_path, "wt") as fout:
        i = 0
        while True:
            hdr = fin.readline()
            if not hdr:
                break
            seq = fin.readline().rstrip("\n")
            plus = fin.readline()
            qual = fin.readline()
            i += 1
            canonical = None
            status = "miss"
            for s in shifts:
                offsets = (offsets_main[0] + s, offsets_main[1] + s, offsets_main[2] + s)
                obs = _extract_bc_24mer(seq, offsets=offsets)
                if obs is None:
                    continue
                canonical, status = _correct_24mer(obs, exact, hamming1)
                if canonical is not None:
                    break
            counter[status] += 1
            tag = canonical if canonical is not None else "Nx24"
            fout.write(hdr.rstrip("\n") + f"_CB:{tag}\n")
            fout.write(seq + "\n")
            fout.write(plus)
            fout.write(qual)
            if args.max_reads and i >= args.max_reads:
                break
    pct = {k: round(100 * v / max(sum(counter.values()), 1), 3) for k, v in counter.items()}
    summary = {
        "fastq": str(in_fastq), "out": str(out_path),
        "n_reads": int(sum(counter.values())), "counts": dict(counter),
        "percent": pct,
    }
    save_json(f"share_{args.label}_demultiplex", summary)
    print(f"Demultiplexed {summary['n_reads']:,} reads -> {out_path}")
    print(f"  exact: {pct.get('exact', 0)}%, hamming1: {pct.get('hamming1', 0)}%, miss: {pct.get('miss', 0)}%")
    return 0


# ---------------------------------------------------------------------------
# Fragments + per-barcode ATAC QC
# ---------------------------------------------------------------------------

def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_fragments_iter(path: Path) -> Iterable[tuple[str, int, int, str, int]]:
    """Yield (chrom, start, end, barcode, count) from a fragments BED."""
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            barcode = parts[3]
            count = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 1
            yield chrom, start, end, barcode, count


def _load_tss_bed(path: Path) -> dict[str, list[int]]:
    """Load TSS positions per chromosome (one position per line)."""
    tss: dict[str, list[int]] = collections.defaultdict(list)
    with _open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                # Treat the start as the TSS unless explicit strand-aware end is given
                tss[chrom].append(start)
            except ValueError:
                continue
    for chrom in tss:
        tss[chrom].sort()
    return tss


def _load_peaks_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    peaks: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    with _open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            peaks[chrom].append((start, end))
    for chrom in peaks:
        peaks[chrom].sort()
    return peaks


def _interval_hit(sorted_intervals: list[tuple[int, int]], pos: int) -> bool:
    """Return True if `pos` lies in any interval. O(log n) bisection."""
    import bisect
    if not sorted_intervals:
        return False
    starts = [s for s, _ in sorted_intervals]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx >= 0 and sorted_intervals[idx][0] <= pos < sorted_intervals[idx][1]:
        return True
    return False


def _nearest_tss_distance(sorted_tss: list[int], pos: int) -> int:
    import bisect
    if not sorted_tss:
        return 10 ** 12
    idx = bisect.bisect_left(sorted_tss, pos)
    candidates = []
    if idx < len(sorted_tss):
        candidates.append(abs(sorted_tss[idx] - pos))
    if idx > 0:
        candidates.append(abs(sorted_tss[idx - 1] - pos))
    return min(candidates)


def cmd_fragment_qc(args: argparse.Namespace) -> int:
    """Compute per-barcode ATAC QC from a SHARE-seq fragments file.

    Outputs columns: ``barcode, n_fragments, reads_tss_2kb,
    reads_flanks_400bp, tss_enrichment, reads_in_peaks, frip``.
    """
    pd = _require_pkg("pandas", "")
    setup_logging()
    frag_path = Path(args.fragments)
    tss = _load_tss_bed(Path(args.tss_bed)) if args.tss_bed else {}
    peaks = _load_peaks_bed(Path(args.peaks_bed)) if args.peaks_bed else {}

    # Pre-compute per-chrom flank intervals for TSS normalization
    # (first 100 bp + last 100 bp of a ±2 kb window — Ma 2020 convention).
    flank_dist = 2000  # ±2kb window
    flank_width = 100
    tss_intervals: dict[str, list[tuple[int, int]]] = {}
    flank_intervals: dict[str, list[tuple[int, int]]] = {}
    for chrom, positions in tss.items():
        tss_intervals[chrom] = sorted([(max(0, p - flank_dist), p + flank_dist) for p in positions])
        flanks = []
        for p in positions:
            flanks.append((max(0, p - flank_dist), max(0, p - flank_dist) + flank_width))
            flanks.append((p + flank_dist - flank_width, p + flank_dist))
        flank_intervals[chrom] = sorted(flanks)

    # Aggregate per barcode
    per_bc = collections.defaultdict(lambda: {
        "n_fragments": 0, "reads_tss_2kb": 0, "reads_flanks_400bp": 0,
        "reads_in_peaks": 0,
    })
    for chrom, start, end, bc, count in _load_fragments_iter(frag_path):
        rec = per_bc[bc]
        rec["n_fragments"] += count
        midpoint = (start + end) // 2
        if chrom in tss_intervals and _interval_hit(tss_intervals[chrom], midpoint):
            rec["reads_tss_2kb"] += count
        if chrom in flank_intervals and _interval_hit(flank_intervals[chrom], midpoint):
            rec["reads_flanks_400bp"] += count
        if chrom in peaks and _interval_hit(peaks[chrom], midpoint):
            rec["reads_in_peaks"] += count

    rows = []
    for bc, rec in per_bc.items():
        # Ma 2020 TSS enrichment normalization with floor 0.2:
        #   tss_enrichment = (reads_tss / window_size_tss) / max(0.2, reads_flank / window_size_flank)
        # Approximated here with raw reads + a window-width correction factor.
        tss_window = 2 * flank_dist  # 4 kb total around each TSS
        flank_window = 2 * flank_width  # 200 bp per TSS (first 100 + last 100)
        tss_rate = rec["reads_tss_2kb"] / max(tss_window, 1)
        flank_rate = max(rec["reads_flanks_400bp"] / max(flank_window, 1), 0.2 / max(tss_window, 1))
        tss_enrichment = tss_rate / flank_rate if flank_rate > 0 else 0.0
        frip = rec["reads_in_peaks"] / rec["n_fragments"] if rec["n_fragments"] > 0 else 0.0
        rows.append({
            "barcode": bc,
            "n_fragments": rec["n_fragments"],
            "reads_tss_2kb": rec["reads_tss_2kb"],
            "reads_flanks_400bp": rec["reads_flanks_400bp"],
            "tss_enrichment": float(tss_enrichment),
            "reads_in_peaks": rec["reads_in_peaks"],
            "frip": float(frip),
        })
    df = pd.DataFrame(rows).sort_values("n_fragments", ascending=False)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_atac_qc.tsv"
    df.to_csv(out, sep="\t", index=False)
    summary = {
        "fragments": str(frag_path),
        "n_barcodes": int(len(df)),
        "total_fragments": int(df["n_fragments"].sum()),
        "median_fragments_per_bc": float(df["n_fragments"].median()) if len(df) else 0.0,
        "median_tss_enrichment": float(df["tss_enrichment"].median()) if len(df) else 0.0,
        "median_frip": float(df["frip"].median()) if len(df) else 0.0,
    }
    save_json(f"share_{args.label}_atac_qc", summary)
    print(f"Per-barcode ATAC QC: {out}")
    print(f"  {summary['n_barcodes']} barcodes, median fragments={summary['median_fragments_per_bc']:.0f}, "
          f"median TSS enrichment={summary['median_tss_enrichment']:.2f}, median FRIP={summary['median_frip']:.3f}")
    return 0


# ---------------------------------------------------------------------------
# RNA modality QC from h5ad
# ---------------------------------------------------------------------------

def cmd_rna_qc(args: argparse.Namespace) -> int:
    """Per-barcode RNA QC from an h5ad sparse gene-count matrix."""
    anndata = _require_pkg("anndata", "Required to read h5ad RNA matrices.")
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    adata = anndata.read_h5ad(args.h5ad)
    counts = adata.X
    if hasattr(counts, "tocsr"):
        counts = counts.tocsr()
    umi = np.asarray(counts.sum(axis=1)).flatten()
    genes = np.asarray((counts > 0).sum(axis=1)).flatten()
    # Mitochondrial fraction if gene names available
    var_names = list(adata.var_names)
    is_mt = np.array([n.startswith("MT-") or n.startswith("mt-") for n in var_names], dtype=bool)
    if is_mt.any():
        mt_counts = np.asarray(counts[:, is_mt].sum(axis=1)).flatten()
        pct_mt = mt_counts / np.maximum(umi, 1)
    else:
        pct_mt = np.zeros(len(umi))
    df = pd.DataFrame({
        "barcode": list(adata.obs_names),
        "umis": umi.astype(int),
        "genes": genes.astype(int),
        "pct_mt": pct_mt,
    })
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_rna_qc.tsv"
    df.to_csv(out, sep="\t", index=False)
    summary = {
        "h5ad": str(args.h5ad),
        "n_barcodes": int(len(df)),
        "median_umis": float(df["umis"].median()),
        "median_genes": float(df["genes"].median()),
        "median_pct_mt": float(df["pct_mt"].median()) if len(df) else 0.0,
    }
    save_json(f"share_{args.label}_rna_qc", summary)
    print(f"Per-barcode RNA QC: {out}")
    print(f"  {summary['n_barcodes']} barcodes, median UMIs={summary['median_umis']:.0f}, "
          f"median genes={summary['median_genes']:.0f}, median %MT={100*summary['median_pct_mt']:.2f}")
    return 0


# ---------------------------------------------------------------------------
# Joint QC (Ma 2020 pass thresholds)
# ---------------------------------------------------------------------------

def cmd_joint_qc(args: argparse.Namespace) -> int:
    """Merge RNA and ATAC per-barcode tables; call cells per Ma 2020 thresholds."""
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    rna = pd.read_csv(args.rna_qc, sep="\t")
    atac = pd.read_csv(args.atac_qc, sep="\t")
    merged = rna.merge(atac, on="barcode", how="outer").fillna({
        "umis": 0, "genes": 0, "pct_mt": np.nan, "n_fragments": 0,
        "tss_enrichment": np.nan, "frip": np.nan,
    })
    rna_pass = (merged["umis"] >= args.min_umis) & (merged["genes"] >= args.min_genes)
    atac_pass = (merged["tss_enrichment"] >= args.min_tss) & (merged["n_fragments"] >= args.min_frags)
    merged["RNA_pass"] = rna_pass
    merged["ATAC_pass"] = atac_pass
    merged["QC"] = "neither"
    merged.loc[rna_pass & ~atac_pass, "QC"] = "RNA only"
    merged.loc[~rna_pass & atac_pass, "QC"] = "ATAC only"
    merged.loc[rna_pass & atac_pass, "QC"] = "both"
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_joint_qc.tsv"
    merged.to_csv(out, sep="\t", index=False)
    summary = {
        "n_barcodes": int(len(merged)),
        "QC_counts": merged["QC"].value_counts().to_dict(),
        "n_RNA_pass": int(rna_pass.sum()),
        "n_ATAC_pass": int(atac_pass.sum()),
        "thresholds": {"min_umis": args.min_umis, "min_genes": args.min_genes,
                        "min_tss": args.min_tss, "min_frags": args.min_frags},
    }
    save_json(f"share_{args.label}_joint_qc", summary)
    print(f"Joint QC table: {out}")
    for k, v in summary["QC_counts"].items():
        print(f"  {k}: {v}")
    return 0


# ---------------------------------------------------------------------------
# Multiplet detection (Jaccard on shared fragment coordinates)
# ---------------------------------------------------------------------------

def cmd_multiplet_detect(args: argparse.Namespace) -> int:
    """Pairwise Jaccard multiplet detection on barcode × fragment-coord set.

    Builds a per-barcode set of unique (chrom, start) tuples (only barcodes
    with ≥ ``min_fragments``), pairs the top-N barcodes by overlap, and
    flags any pair whose Jaccard exceeds the null cutoff as a multiplet.
    The null is the 99-th percentile of Jaccard across 10,000 random pairs.
    """
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    rng = np.random.default_rng(args.seed)
    bc_sets: dict[str, set[tuple[str, int]]] = collections.defaultdict(set)
    for chrom, start, _, bc, _ in _load_fragments_iter(Path(args.fragments)):
        bc_sets[bc].add((chrom, start))
    bc_sets = {bc: s for bc, s in bc_sets.items() if len(s) >= args.min_fragments}
    barcodes = list(bc_sets.keys())
    logging.info("Multiplet scan: %d barcodes pass min_fragments=%d", len(barcodes), args.min_fragments)
    if len(barcodes) < 4:
        print("Too few barcodes for multiplet scan.")
        return 0

    # Build a null distribution of Jaccard for random pairs
    null = []
    for _ in range(min(args.null_samples, len(barcodes) * (len(barcodes) - 1) // 2)):
        i, j = rng.choice(len(barcodes), size=2, replace=False)
        a, b = bc_sets[barcodes[i]], bc_sets[barcodes[j]]
        union = len(a | b)
        null.append(len(a & b) / union if union > 0 else 0.0)
    cutoff = float(np.percentile(null, 99)) if null else 0.0

    # Sample candidate pairs by sharing at least one fragment via inverted index
    inverted: dict[tuple[str, int], list[int]] = collections.defaultdict(list)
    for idx, bc in enumerate(barcodes):
        for coord in bc_sets[bc]:
            inverted[coord].append(idx)
    cand_pairs: set[tuple[int, int]] = set()
    for coord, bc_list in inverted.items():
        if len(bc_list) < 2:
            continue
        for i in range(len(bc_list)):
            for j in range(i + 1, len(bc_list)):
                if bc_list[i] != bc_list[j]:
                    cand_pairs.add(tuple(sorted((bc_list[i], bc_list[j]))))
                if len(cand_pairs) > args.max_pairs:
                    break
            if len(cand_pairs) > args.max_pairs:
                break
        if len(cand_pairs) > args.max_pairs:
            break

    rows = []
    for i, j in cand_pairs:
        a, b = bc_sets[barcodes[i]], bc_sets[barcodes[j]]
        union = len(a | b)
        jac = len(a & b) / union if union > 0 else 0.0
        if jac >= cutoff:
            rows.append({
                "Barcode1": barcodes[i], "Barcode2": barcodes[j],
                "Common": len(a & b), "JaccardIndex": jac,
                "PrimaryBarcode": barcodes[i] if len(a) >= len(b) else barcodes[j],
                "IsMultiplet": True,
            })
    out_df = pd.DataFrame(rows)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_multiplets.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    print(f"Multiplet pairs: {out} ({len(out_df)} flagged | Jaccard cutoff={cutoff:.4f})")
    save_json(f"share_{args.label}_multiplets",
              {"cutoff": cutoff, "n_multiplet_pairs": int(len(out_df))})
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "SHARESEQ_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Skill: SHARE-seq Joint ATAC + RNA QC

Use this skill for end-to-end QC and joint-cell analytics on SHARE-seq
deposits (Ma et al. *Cell* 2020 — simultaneous scATAC + scRNA per cell,
shared round-1/2/3 split-pool barcodes). The math is a clean-room
reimplementation of the QC algorithms in
[broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline)
(MIT, 2021).

This skill is designed to consume the **processed deposits** that IGVF
Portal SHARE-seq AnalysisSets publish — fragments BED + h5ad — and does
not re-run FASTQ → BAM alignment. Use the upstream Broad WDL pipeline
(or the IGVF Portal ingestion service) for raw-FASTQ processing.

## Commands

```bash
# 1. Discover IGVF Portal SHARE-seq AnalysisSets / MeasurementSets
python3 Scripts/share_seq_skill.py pull-portal --limit 50 --label survey

# 2. (Optional) Demultiplex a raw R2 FASTQ with the SHARE-seq round-1/2/3
#    24-mer barcode + 1-Hamming-mismatch + ±1 bp shift correction
python3 Scripts/share_seq_skill.py demultiplex-bcs \\
    --fastq sample_R2.fastq.gz --whitelist whitelist_24mer.txt \\
    --out sample_R2.tagged.fastq.gz --label demo

# 3. Per-barcode ATAC QC from a fragments BED + (optional) TSS BED + peaks BED
python3 Scripts/share_seq_skill.py fragment-qc \\
    --fragments fragments.bed.gz --tss-bed tss.bed --peaks-bed peaks.bed \\
    --label demo

# 4. Per-barcode RNA QC from an h5ad sparse gene-count matrix
python3 Scripts/share_seq_skill.py rna-qc \\
    --h5ad rna.h5ad --label demo

# 5. Joint per-barcode pass calls (Ma 2020 thresholds)
python3 Scripts/share_seq_skill.py joint-qc \\
    --atac-qc Docs/SHAREseq/<ts>_demo_atac_qc.tsv \\
    --rna-qc Docs/SHAREseq/<ts>_demo_rna_qc.tsv --label demo

# 6. Multiplet detection by pairwise Jaccard over fragment-coord sets
python3 Scripts/share_seq_skill.py multiplet-detect \\
    --fragments fragments.bed.gz --label demo
```

## Input file conventions

| Input | Expected format |
|---|---|
| `fragments.bed[.gz]` | 4-col + count BED (chrom, start, end, barcode, count). Matches the SHARE-seq pipeline `bam_to_fragments.py` schema (`start = read.reference_start + 4; end = start + tlen - 4`). |
| `rna.h5ad` | AnnData with barcode × gene; mitochondrial genes detected by `MT-` / `mt-` prefix. |
| `tss.bed` | BED of TSS positions (any annotation; chrom + start). |
| `peaks.bed` | BED of MACS / iterative-LSI peaks. |

## Default cell-calling thresholds (Ma 2020)

| Modality | Metric | Pass |
|---|---|---|
| RNA | UMIs ≥ `min_umis` (default 100) AND genes ≥ `min_genes` (default 200) | both required |
| ATAC | fragments ≥ `min_frags` (default 100) AND TSS enrichment ≥ `min_tss` (default 4) | both required |
| Joint | RNA_pass AND ATAC_pass | "both" |

## Citation

- **Ma S et al. (2020)** *Cell* 183:1103–1116. "Chromatin potential identified
  by shared single-cell profiling of RNA and chromatin." doi:10.1016/j.cell.2020.09.056
- Pipeline reference: [broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline) (MIT)

## License-clean Python stack

`pandas` (BSD-3), `numpy` (BSD-3), `scipy` (BSD-3), `anndata` (BSD-3),
`statsmodels` (BSD-3). All non-GPL. Only the FASTQ demultiplex is a hand-
rolled scan over gz files using the standard library — no pysam, no
GPL bowtie2 chain at runtime.
""",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="SHARE-seq joint ATAC+RNA QC (clean-room).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("pull-portal", help="Discover IGVF Portal SHARE-seq AnalysisSets/MeasurementSets.")
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--label", default="share_portal_survey")

    p = sub.add_parser("demultiplex-bcs", help="Tag R2 reads with SHARE-seq round-1/2/3 24-mer barcodes.")
    p.add_argument("--fastq", required=True)
    p.add_argument("--whitelist", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--r1-offset", type=int, default=14)
    p.add_argument("--r2-offset", type=int, default=52)
    p.add_argument("--r3-offset", type=int, default=90)
    p.add_argument("--shift-correct", action="store_true")
    p.add_argument("--max-reads", type=int, default=0, help="Stop after this many reads (0 = no limit).")
    p.add_argument("--label", default="share_demultiplex")

    p = sub.add_parser("fragment-qc", help="Per-barcode ATAC QC from a fragments BED.")
    p.add_argument("--fragments", required=True)
    p.add_argument("--tss-bed", default=None)
    p.add_argument("--peaks-bed", default=None)
    p.add_argument("--label", default="share_atac_qc")

    p = sub.add_parser("rna-qc", help="Per-barcode RNA QC from an h5ad.")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--label", default="share_rna_qc")

    p = sub.add_parser("joint-qc", help="Join RNA + ATAC per-barcode tables and call cells.")
    p.add_argument("--rna-qc", required=True)
    p.add_argument("--atac-qc", required=True)
    p.add_argument("--min-umis", type=int, default=100)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-tss", type=float, default=4.0)
    p.add_argument("--min-frags", type=int, default=100)
    p.add_argument("--label", default="share_joint")

    p = sub.add_parser("multiplet-detect", help="Jaccard multiplet scan over fragment-coord sets.")
    p.add_argument("--fragments", required=True)
    p.add_argument("--min-fragments", type=int, default=1000)
    p.add_argument("--max-pairs", type=int, default=200000)
    p.add_argument("--null-samples", type=int, default=10000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--label", default="share_multiplet")

    sub.add_parser("write-playbook", help="Emit the skill's markdown playbook.")

    args = parser.parse_args(argv)
    if args.command == "pull-portal":
        return cmd_pull_portal(args)
    if args.command == "demultiplex-bcs":
        return cmd_demultiplex_bcs(args)
    if args.command == "fragment-qc":
        return cmd_fragment_qc(args)
    if args.command == "rna-qc":
        return cmd_rna_qc(args)
    if args.command == "joint-qc":
        return cmd_joint_qc(args)
    if args.command == "multiplet-detect":
        return cmd_multiplet_detect(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
