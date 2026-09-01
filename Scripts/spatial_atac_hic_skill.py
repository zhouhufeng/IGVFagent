#!/usr/bin/env python3
"""IGVF agent Spatial-ATAC-Hi-C analytics skill.

Spatially resolved co-profiling of 3D genome organisation and chromatin
accessibility on tissue slides, after Wang P., Wang J., Wang Q.,
Youngblood M. W. et al., "Spatial chromatin architecture and
accessibility co-profiling of mammalian tissues", *Nature Methods* 2026
(doi:10.1038/s41592-026-03217-4; GEO **GSE307620**).

The assay lays a 50 x 50 microfluidic barcode grid over a tissue
section, giving 2,500 spatial *pixels*, each carrying both a Hi-C
contact set and an ATAC fragment set. This skill takes those per-pixel
deposits and produces the analyses the paper reports.

Scope — this skill consumes **processed deposits**, the same line
``share_seq_skill`` draws. It starts from:

  * contact tables   4DN ``*.pairs[.gz]`` **or** the tabix-indexed
                     contact TSVs GEO deposits. GSE307620 ships
                     ``*.hic.fragments.sorted.header.tsv.gz``, whose
                     first row is a bare column header rather than a
                     ``#columns:`` comment; both shapes, and the common
                     column-name synonyms (``chr1``/``chrom1``, ...),
                     are handled transparently.
  * ``*.fragments.tsv.gz``    ATAC fragments (chrom, start, end, barcode,
                              count)
  * ``*_tissue_positions_list.csv[.gz]`` spatial barcode -> (row, col)
                              grid, as emitted by AtlasXBrowser
  * a gene model (GTF/GFF/BED) for the gene-level scores

It does **not** re-run FASTQ -> alignment. Read alignment is the one
step in the upstream workflow with no tractable clean-room form
(Cell Ranger ATAC; runHiC, which is GPL-3.0 and therefore excluded from
this Apache-2 codebase by policy). Bring aligned pairs/fragments — the
GEO series publishes them — or run your own aligner first.

Subcommands
-----------
  pull-geo        Discover and fetch the GSE307620 (or any GSE) deposits.
  pixel-demux     Split a barcoded pairs file into the 50 x 50 pixel grid.
  qc              Per-pixel cis/trans/long-range contacts + TSS enrichment.
  gas             Gene activity score matrix from ATAC fragments.
  gad             Gene-associated domain score matrix from Hi-C pairs.
  matrix          Build a binned contact matrix from pairs.
  impute          scHiCluster-style convolution + random-walk imputation.
  compartment     A/B compartment PC1 at 100 kb.
  cnv             Per-pixel copy-number ratio (+ optional HMM segmentation).
  loops           Quantify loops per pixel, ANOVA for cluster-specific loops,
                  APA pileup.
  viz             Render per-pixel values back into tissue space.
  write-playbook  Emit the skill's markdown playbook.

License: Apache-2.0. The reference implementations
(wangjuan001/Spatial-ATAC-Hi-C, MIT; scHiCluster, MIT; cooltools, MIT;
NeoLoopFinder, MIT; SnapATAC2, MIT) are followed from their published
descriptions only — no source is copied, imported or vendored. runHiC
and Trim Galore are GPL-3.0 and are referenced as external tools only,
never as runtime dependencies. Numerics live in ``_spatial_hic_math``.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import logging
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "SpatialATACHiC"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

# The two microfluidic barcode sets are 50 channels each, giving the
# 2,500-pixel grid quoted throughout the paper.
GRID_N = 50

# Long-range threshold used for the paper's Fig. 1f / Extended Data 4.
LONG_RANGE_BP = 10_000


# ---------------------------------------------------------------------------
# Setup helpers
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = (LOG_DIR /
                f"spatial_atac_hic_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_"
                   for ch in (label or "run"))


def run_dir(label: str) -> Path:
    d = REPORT_DIR / f"{timestamp()}_{safe_label(label)}"
    (d / "Plots").mkdir(parents=True, exist_ok=True)
    return d


def _require(name: str, extra: str = "analysis") -> Any:
    """Import a scientific dependency or exit with an actionable hint."""
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(
            f"Missing dependency '{name}'. Install with:\n"
            f"    pip install 'igvfagent[{extra}]'\n"
            f"  (or: pip install {name})"
        ) from exc


def _open_text(path: "str | Path") -> Any:
    """Open plain or gzipped text transparently."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")
    return path


def write_tsv(path: Path, rows: "list[dict]", columns: "list[str] | None" = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(columns or (rows[0].keys() if rows else []))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def _median(values: "list[float]") -> float:
    vals = sorted(v for v in values if v == v)
    if not vals:
        return float("nan")
    mid = len(vals) // 2
    if len(vals) % 2:
        return float(vals[mid])
    return float((vals[mid - 1] + vals[mid]) / 2.0)


# ---------------------------------------------------------------------------
# Pairs I/O
# ---------------------------------------------------------------------------

# 4DN pairs default column order when no `#columns:` header is present.
_DEFAULT_PAIRS_COLUMNS = ["readID", "chrom1", "pos1", "chrom2", "pos2",
                          "strand1", "strand2"]


class PairsHeader:
    """Parsed header of a contact file.

    Covers both shapes this skill meets in the wild: a 4DN ``.pairs``
    file with ``#columns:`` / ``#chromsize:`` comment lines, and a plain
    tabix-indexed TSV whose first row is a bare column header — which is
    what GSE307620 deposits as
    ``*.hic.fragments.sorted.header.tsv.gz``.
    """

    def __init__(self, columns: "list[str]", chromsizes: "dict[str, int]",
                 lines: "list[str]", *, has_header_row: bool = False):
        self.columns = columns
        self.chromsizes = chromsizes
        self.lines = lines
        self.has_header_row = has_header_row

    def index(self, name: str) -> int:
        try:
            return self.columns.index(name)
        except ValueError as exc:
            raise SystemExit(
                f"pairs file has no '{name}' column; found {self.columns}"
            ) from exc


# Column-name synonyms. GEO deposits and pairtools variants label the
# same five fields half a dozen ways; normalise so `header.index("chrom1")`
# works regardless of which spelling the file used.
_COLUMN_ALIASES = {
    "readid": "readID", "read_id": "readID", "name": "readID",
    "qname": "readID", "readname": "readID",
    "chr1": "chrom1", "chrom_1": "chrom1", "chrA": "chrom1",
    "chrom1": "chrom1", "chr_1": "chrom1", "chrom": "chrom1",
    "chr2": "chrom2", "chrom_2": "chrom2", "chrB": "chrom2",
    "chrom2": "chrom2", "chr_2": "chrom2",
    "pos1": "pos1", "position1": "pos1", "start1": "pos1",
    "pos_1": "pos1", "posA": "pos1", "start": "pos1",
    "pos2": "pos2", "position2": "pos2", "start2": "pos2",
    "pos_2": "pos2", "posB": "pos2",
    "strand1": "strand1", "strand2": "strand2",
}


def _normalise_columns(names: "list[str]") -> "list[str]":
    out = []
    for n in names:
        key = n.strip().lstrip("#")
        out.append(_COLUMN_ALIASES.get(key.lower(), key))
    return out


def _looks_like_header_row(fields: "list[str]") -> bool:
    """True when a row carries no integer field.

    Every real contact row has integer positions, so a row with none is
    a column-name row. This is what lets the skill read the plain
    tabix-indexed TSVs GEO deposits alongside true 4DN pairs.
    """
    for f in fields:
        try:
            int(f)
            return False
        except ValueError:
            continue
    return True


def read_pairs_header(path: "str | Path",
                      *, columns_override: "Optional[list[str]]" = None
                      ) -> PairsHeader:
    """Read the header of a pairs / contact-TSV file."""
    columns = list(_DEFAULT_PAIRS_COLUMNS)
    chromsizes: "dict[str, int]" = {}
    lines: "list[str]" = []
    has_header_row = False
    from_comment = False

    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                lines.append(line.rstrip("\n"))
                if line.startswith("#columns:"):
                    columns = _normalise_columns(line.split(":", 1)[1].split())
                    from_comment = True
                elif line.startswith("#chromsize:"):
                    parts = line.split(":", 1)[1].split()
                    if len(parts) >= 2:
                        try:
                            chromsizes[parts[0]] = int(parts[1])
                        except ValueError:
                            pass
                continue
            # First non-comment line. If it has no integer field it is a
            # bare column header (the GEO `*.header.tsv.gz` shape).
            fields = line.rstrip("\n").split("\t")
            if not from_comment and len(fields) > 1 and _looks_like_header_row(fields):
                columns = _normalise_columns(fields)
                has_header_row = True
                lines.append("#columns: " + " ".join(columns))
            break

    if columns_override:
        columns = _normalise_columns(columns_override)
    return PairsHeader(columns, chromsizes, lines,
                       has_header_row=has_header_row)


def iter_pairs(path: "str | Path") -> "Iterator[list[str]]":
    """Yield the data rows of a pairs / contact-TSV file.

    Skips ``#`` comments and, if present, the single bare column-header
    row. Detection is per-file and costs one string scan.
    """
    first = True
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if first:
                first = False
                if len(fields) > 1 and _looks_like_header_row(fields):
                    continue
            yield fields


def read_chrom_sizes(path: "str | Path") -> "dict[str, int]":
    sizes: "dict[str, int]" = {}
    with _open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    if not sizes:
        raise SystemExit(f"no chrom sizes parsed from {path}")
    return sizes


# Contact-file extensions found under a --pairs-dir. `.pairs` is what
# `pixel-demux` writes and what runHiC emits; the `.tsv[.gz]` forms cover
# the tabix-indexed contact tables GEO deposits (GSE307620 ships
# `*.hic.fragments.sorted.header.tsv.gz`).
_CONTACT_SUFFIXES = (".pairs", ".pairs.gz", ".tsv", ".tsv.gz",
                     ".contacts", ".contacts.gz")


def discover_pairs(pairs_dir: "str | Path",
                   *, pattern: "Optional[str]" = None) -> "list[Path]":
    """Contact files under a directory, sorted by name.

    A path to a single file is taken as-is whatever its extension, so an
    oddly-named deposit can always be handed over directly.
    """
    d = Path(pairs_dir)
    if d.is_file():
        return [d]
    if not d.is_dir():
        raise SystemExit(f"not a file or directory: {d}")
    if pattern:
        out = sorted(q for q in d.rglob(pattern) if q.is_file())
        if not out:
            raise SystemExit(f"no files matching {pattern!r} under {d}")
        return out
    out = sorted(q for q in d.rglob("*")
                 if q.is_file() and q.name.endswith(_CONTACT_SUFFIXES))
    if not out:
        raise SystemExit(
            f"no contact files under {d} "
            f"(looked for {', '.join(_CONTACT_SUFFIXES)}). "
            f"Pass --pairs-glob to match a different naming.")
    return out


_PIXEL_ID_RE = re.compile(r"(?<![0-9A-Za-z])(\d{1,3}x\d{1,3})(?![0-9A-Za-z])")


def pixel_name(path: Path) -> str:
    """Pixel id from a per-pixel contact filename.

    Real deposits carry compound extensions --- GSE307620 names its
    contact tables ``*.hic.fragments.sorted.header.tsv.gz`` --- so
    stripping only the final suffix leaves
    ``01x01.hic.fragments.sorted.header`` as the "pixel id", which then
    matches no cluster label and lands nowhere on the tissue grid.

    Prefer an explicit ``AAxBB`` token anywhere in the name; otherwise
    fall back to the segment before the first dot.
    """
    name = path.name
    for suffix in sorted(_CONTACT_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    m = _PIXEL_ID_RE.search(name)
    if m:
        return m.group(1)
    return name.split(".", 1)[0]


# ---------------------------------------------------------------------------
# Barcode / spatial-grid I/O
# ---------------------------------------------------------------------------

def read_barcode_whitelist(path: "str | Path") -> "list[str]":
    """Read a barcode whitelist: one sequence per line, or TSV name<TAB>seq.

    Accepts FASTA too, so the supplementary barcode tables from the paper
    can be handed over in whatever shape they were downloaded in.
    """
    seqs: "list[str]" = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(">"):
                continue
            parts = line.replace(",", "\t").split("\t")
            # Take the last field that looks like DNA — handles both
            # "A1<TAB>ACGT..." and a bare sequence column.
            seq = None
            for field in reversed(parts):
                f = field.strip().upper()
                if f and set(f) <= set("ACGTN"):
                    seq = f
                    break
            if seq:
                seqs.append(seq)
    if not seqs:
        raise SystemExit(f"no barcode sequences parsed from {path}")
    return seqs


def barcode_to_pixel_map(
    barcode_a: "str | Path",
    barcode_b: "str | Path",
    *,
    layout: str = "BA",
) -> "dict[str, str]":
    """Concatenated barcode -> ``AAxBB`` pixel id.

    ``pixel-demux`` names per-pixel pairs files by grid index, while an
    ATAC fragments file is keyed by the raw concatenated barcode. Without
    this translation the GAS matrix (from fragments) and the GAD matrix
    (from pairs) are keyed differently and cannot be joined, plotted on
    the same grid, or handed the same cluster labels.
    """
    wl_a = read_barcode_whitelist(barcode_a)
    wl_b = read_barcode_whitelist(barcode_b)
    len_a, len_b = len(wl_a[0]), len(wl_b[0])
    if any(len(s) != len_a for s in wl_a) or any(len(s) != len_b for s in wl_b):
        raise SystemExit("barcode whitelists must each be fixed-length")
    out: "dict[str, str]" = {}
    for ai, a in enumerate(wl_a, 1):
        for bi, b in enumerate(wl_b, 1):
            key = (b + a) if layout == "BA" else (a + b)
            out[key] = f"{ai:02d}x{bi:02d}"
    return out


def read_spatial_positions(path: "str | Path") -> "dict[str, tuple[int, int, int]]":
    """Read an AtlasXBrowser / Visium ``tissue_positions`` CSV.

    Returns ``barcode -> (in_tissue, row, col)``. Both the headered and
    the classic headerless five-column form are accepted.
    """
    out: "dict[str, tuple[int, int, int]]" = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 4:
                continue
            if not parts[1].lstrip("-").isdigit():
                continue  # header row
            try:
                out[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]))
            except ValueError:
                continue
    if not out:
        raise SystemExit(f"no spatial positions parsed from {path}")
    return out


# ---------------------------------------------------------------------------
# Gene model
# ---------------------------------------------------------------------------

class GeneModel:
    """Gene bodies plus TSS, with an interval index for fast overlap."""

    def __init__(self, genes: "list[dict]"):
        self.genes = genes
        self._by_chrom: "dict[str, list[dict]]" = defaultdict(list)
        for g in genes:
            self._by_chrom[g["chrom"]].append(g)
        self._starts: "dict[str, list[int]]" = {}
        for chrom, gs in self._by_chrom.items():
            gs.sort(key=lambda g: g["start"])
            self._starts[chrom] = [g["start"] for g in gs]
        # Longest gene per chromosome bounds how far back a query must
        # scan to be sure it has not missed an overlap.
        self._max_len = {
            chrom: max((g["end"] - g["start"] for g in gs), default=0)
            for chrom, gs in self._by_chrom.items()
        }

    def __len__(self) -> int:
        return len(self.genes)

    @property
    def chroms(self) -> "list[str]":
        return sorted(self._by_chrom)

    def overlapping(self, chrom: str, pos: int) -> "list[dict]":
        """Genes whose (already-extended) span contains ``pos``."""
        gs = self._by_chrom.get(chrom)
        if not gs:
            return []
        starts = self._starts[chrom]
        lo = bisect.bisect_right(starts, pos)
        span = self._max_len.get(chrom, 0)
        out = []
        i = lo - 1
        while i >= 0:
            g = gs[i]
            if g["start"] < pos - span:
                break
            if g["start"] <= pos < g["end"]:
                out.append(g)
            i -= 1
        return out

    def tss_list(self) -> "list[tuple[str, int]]":
        return [(g["chrom"], g["tss"]) for g in self.genes]


_GTF_GENE_NAME = re.compile(r'gene_name\s+"([^"]+)"')
_GTF_GENE_ID = re.compile(r'gene_id\s+"([^"]+)"')


def read_gene_model(
    path: "str | Path",
    *,
    upstream: int = 2000,
    downstream: int = 0,
    feature: str = "gene",
) -> GeneModel:
    """Read genes from a GTF/GFF or a BED6.

    The span is extended ``upstream`` bases past the TSS on the gene's
    own strand, matching SnapATAC2's ``make_gene_matrix`` convention (the
    paper counts Tn5 insertions "within the promoter region of each gene
    (2 kb upstream of the TSS) and gene body region").

    Parameters
    ----------
    upstream, downstream
        Promoter padding, in bases, applied strand-aware.
    feature
        GTF feature column to keep. Ignored for BED.
    """
    p = Path(path)
    genes: "list[dict]" = []
    is_bed = p.name.replace(".gz", "").lower().endswith((".bed", ".bed6"))

    with _open_text(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if is_bed:
                if len(f) < 3:
                    continue
                chrom, start, end = f[0], int(f[1]), int(f[2])
                name = f[3] if len(f) > 3 else f"{chrom}:{start}"
                strand = f[5] if len(f) > 5 else "+"
            else:
                if len(f) < 9 or f[2] != feature:
                    continue
                chrom, start, end = f[0], int(f[3]) - 1, int(f[4])
                strand = f[6]
                attrs = f[8]
                m = _GTF_GENE_NAME.search(attrs) or _GTF_GENE_ID.search(attrs)
                name = m.group(1) if m else f"{chrom}:{start}"
            if end <= start:
                continue
            tss = start if strand != "-" else end - 1
            if strand != "-":
                lo, hi = start - upstream, end + downstream
            else:
                lo, hi = start - downstream, end + upstream
            genes.append({
                "gene": name, "chrom": chrom,
                "start": max(0, lo), "end": hi,
                "tss": tss, "strand": strand,
            })
    if not genes:
        raise SystemExit(
            f"no {'BED intervals' if is_bed else f'{feature!r} features'} "
            f"parsed from {p}")
    return GeneModel(genes)


# ---------------------------------------------------------------------------
# pull-geo
# ---------------------------------------------------------------------------

def cmd_pull_geo(args: argparse.Namespace) -> int:
    """Discover (and optionally download) the paper's GEO deposits."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import geo_retrieval_skill as geo
    except Exception:
        from igvfagent import geo_retrieval_skill as geo  # type: ignore

    out = run_dir(args.label)
    gse = args.gse

    logging.info("listing GEO supplementary files for %s", gse)
    files = geo.list_gse_files(gse)
    rows = [{"name": f.get("name", ""), "size": f.get("size", ""),
             "url": f.get("url", "")} for f in files]
    tsv = write_tsv(out / "geo_files.tsv", rows, ["name", "size", "url"])

    summary: "dict[str, Any]" = {"gse": gse, "n_files": len(rows)}
    try:
        soft = geo.fetch_soft(gse, "quick")
        series = geo.summarize_series(soft, gse)
        summary["series"] = series
    except Exception as exc:  # pragma: no cover - network shape varies
        logging.warning("could not fetch SOFT metadata: %s", exc)

    downloaded: "list[str]" = []
    if args.download:
        pattern = re.compile(args.download)
        dest_dir = Path(args.dest) if args.dest else (DATA_DIR / "SpatialATACHiC" / gse)
        dest_dir.mkdir(parents=True, exist_ok=True)
        for f in files:
            name = f.get("name", "")
            if not pattern.search(name):
                continue
            dest = dest_dir / name
            if dest.exists() and dest.stat().st_size > 0:
                logging.info("already have %s", name)
                downloaded.append(str(dest))
                continue
            logging.info("downloading %s", name)
            geo.download_file(f["url"], dest)
            downloaded.append(str(dest))
        summary["downloaded"] = downloaded
        summary["dest_dir"] = str(dest_dir)

    js = write_json(out / "geo_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"File listing: {tsv}  ({len(rows)} supplementary files)")
    print(f"Summary: {js}")
    if downloaded:
        print(f"Downloaded {len(downloaded)} file(s) to {summary['dest_dir']}")
    else:
        print("No files downloaded. Re-run with --download '<regex>' to fetch, "
              "e.g. --download 'pairs|fragments|positions'")
    return 0


# ---------------------------------------------------------------------------
# pixel-demux
# ---------------------------------------------------------------------------

def _barcode_from_row(row: "list[str]", *, field: Optional[int],
                      readid_idx: int, sep: str, part: int) -> Optional[str]:
    if field is not None:
        return row[field] if field < len(row) else None
    rid = row[readid_idx] if readid_idx < len(row) else ""
    chunks = rid.split(sep)
    if not chunks:
        return None
    try:
        return chunks[part]
    except IndexError:
        return None


def cmd_pixel_demux(args: argparse.Namespace) -> int:
    """Split one barcoded pairs file into the 50 x 50 pixel grid.

    The barcode is the concatenation of the two microfluidic rounds. The
    upstream ``bcsplit.py`` writes it as barcode-B followed by barcode-A
    (offsets 22:30 then 60:68 of Read2), which is the ``BA`` default
    here; pass ``--layout AB`` if your pairs were built the other way
    round.
    """
    header = read_pairs_header(args.pairs)
    readid_idx = header.index("readID")

    wl_a = read_barcode_whitelist(args.barcode_a)
    wl_b = read_barcode_whitelist(args.barcode_b)
    logging.info("whitelists: A=%d B=%d", len(wl_a), len(wl_b))
    len_a, len_b = len(wl_a[0]), len(wl_b[0])
    if any(len(s) != len_a for s in wl_a) or any(len(s) != len_b for s in wl_b):
        raise SystemExit("barcode whitelists must each be fixed-length")

    idx_a = {s: i + 1 for i, s in enumerate(wl_a)}
    idx_b = {s: i + 1 for i, s in enumerate(wl_b)}

    out = run_dir(args.label)
    pix_dir = out / "pixels"
    pix_dir.mkdir(parents=True, exist_ok=True)

    field = args.barcode_field - 1 if args.barcode_field else None
    handles: "dict[str, Any]" = {}
    counts: "Counter[str]" = Counter()
    n_rows = n_assigned = n_unmatched = 0

    header_block = "\n".join(header.lines) + "\n" if header.lines else ""
    try:
        for row in iter_pairs(args.pairs):
            n_rows += 1
            if args.max_pairs and n_rows > args.max_pairs:
                n_rows -= 1
                break
            bc = _barcode_from_row(row, field=field, readid_idx=readid_idx,
                                   sep=args.readid_sep, part=args.readid_part)
            if not bc:
                n_unmatched += 1
                continue
            bc = bc.strip().upper()
            if args.layout == "BA":
                b_seq, a_seq = bc[:len_b], bc[len_b:len_b + len_a]
            else:
                a_seq, b_seq = bc[:len_a], bc[len_a:len_a + len_b]
            ai, bi = idx_a.get(a_seq), idx_b.get(b_seq)
            if ai is None or bi is None:
                n_unmatched += 1
                continue
            pid = f"{ai:02d}x{bi:02d}"
            fh = handles.get(pid)
            if fh is None:
                fh = (pix_dir / f"{pid}.pairs").open("w")
                if header_block:
                    fh.write(header_block)
                handles[pid] = fh
            fh.write("\t".join(row) + "\n")
            counts[pid] += 1
            n_assigned += 1
    finally:
        for fh in handles.values():
            fh.close()

    rows = [{"pixel": pid, "a_index": int(pid.split("x")[0]),
             "b_index": int(pid.split("x")[1]), "contacts": c}
            for pid, c in sorted(counts.items())]
    tsv = write_tsv(out / "pixel_counts.tsv", rows,
                    ["pixel", "a_index", "b_index", "contacts"])
    summary = {
        "pairs_input": str(args.pairs),
        "layout": args.layout,
        "rows_read": n_rows,
        "assigned": n_assigned,
        "unmatched": n_unmatched,
        "assigned_fraction": (n_assigned / n_rows) if n_rows else 0.0,
        "pixels_with_data": len(counts),
        "grid_capacity": len(wl_a) * len(wl_b),
        "median_contacts_per_pixel": _median([c for c in counts.values()]),
    }
    js = write_json(out / "demux_summary.json", summary)

    print(f"Output dir: {out}")
    print(f"Per-pixel pairs: {pix_dir}  ({len(counts)} pixels)")
    print(f"Pixel counts: {tsv}")
    print(f"Summary: {js}")
    print(f"Assigned {n_assigned:,}/{n_rows:,} pairs "
          f"({summary['assigned_fraction']:.1%}); {n_unmatched:,} unmatched")
    return 0


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------

def pairs_qc(path: "str | Path", *, long_range_bp: int = LONG_RANGE_BP) -> "dict[str, Any]":
    """cis/trans/long-range statistics for one pairs file.

    Mirrors what ``pairtools stats`` reports and what the paper quotes:
    total unique contacts, the cis fraction (88.1-90.3% across their
    samples) and the long-range (>= 10 kb) share of cis contacts
    (24-33.3%).
    """
    header = read_pairs_header(path)
    i_c1 = header.index("chrom1")
    i_c2 = header.index("chrom2")
    i_p1 = header.index("pos1")
    i_p2 = header.index("pos2")

    total = cis = trans = long_range = 0
    for row in iter_pairs(path):
        if len(row) <= max(i_c1, i_c2, i_p1, i_p2):
            continue
        total += 1
        if row[i_c1] != row[i_c2]:
            trans += 1
            continue
        cis += 1
        try:
            dist = abs(int(row[i_p2]) - int(row[i_p1]))
        except ValueError:
            continue
        if dist >= long_range_bp:
            long_range += 1

    return {
        "total_contacts": total,
        "cis_contacts": cis,
        "trans_contacts": trans,
        "long_range_contacts": long_range,
        "cis_fraction": (cis / total) if total else 0.0,
        # The paper's "fraction of long-range contacts over the unique
        # intra-chromosomal contacts" (Fig. 1f) — denominator is cis, not
        # total. Both are reported so neither reading is ambiguous.
        "long_range_ratio": (long_range / cis) if cis else 0.0,
        "long_range_over_total": (long_range / total) if total else 0.0,
    }


def tss_enrichment(
    fragments: "str | Path",
    model: GeneModel,
    *,
    window: int = 2000,
    centre: int = 50,
    flank: int = 100,
    per_barcode: bool = True,
    max_rows: int = 0,
) -> "dict[str, Any]":
    """ArchR-style TSS enrichment from an ATAC fragments file.

    Following the paper's Methods: "the average accessibility within the
    50-bp window centered on each TSS ... normalized by the average
    accessibility of the TSS flanking regions (+/-1,900-2,000 bp)".

    Both fragment ends are counted as Tn5 insertions.
    """
    tss_by_chrom: "dict[str, list[int]]" = defaultdict(list)
    for chrom, pos in model.tss_list():
        tss_by_chrom[chrom].append(pos)
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort()

    half_c = centre // 2
    flank_lo = window - flank  # 1,900 with the defaults

    centre_counts: "Counter[str]" = Counter()
    flank_counts: "Counter[str]" = Counter()
    total_centre = total_flank = 0
    n_rows = 0

    with _open_text(fragments) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            n_rows += 1
            if max_rows and n_rows > max_rows:
                break
            chrom = f[0]
            positions = tss_by_chrom.get(chrom)
            if not positions:
                continue
            bc = f[3] if per_barcode and len(f) > 3 else "__all__"
            try:
                ends = (int(f[1]), int(f[2]) - 1)
            except ValueError:
                continue
            for pos in ends:
                lo = bisect.bisect_left(positions, pos - window)
                hi = bisect.bisect_right(positions, pos + window)
                for t in positions[lo:hi]:
                    d = abs(pos - t)
                    if d <= half_c:
                        centre_counts[bc] += 1
                        total_centre += 1
                    elif flank_lo <= d <= window:
                        flank_counts[bc] += 1
                        total_flank += 1

    # Densities, not raw counts: the centre window is 50 bp and the two
    # flanks together are 2 * 100 bp, so a raw ratio would understate
    # enrichment by 4x.
    centre_bp = float(centre)
    flank_bp = float(2 * flank)

    per_bc: "dict[str, float]" = {}
    if per_barcode:
        for bc in set(centre_counts) | set(flank_counts):
            fl = flank_counts.get(bc, 0) / flank_bp
            ce = centre_counts.get(bc, 0) / centre_bp
            per_bc[bc] = (ce / fl) if fl > 0 else float("nan")

    overall = ((total_centre / centre_bp) / (total_flank / flank_bp)
               if total_flank else float("nan"))
    return {
        "tss_enrichment": overall,
        "centre_insertions": total_centre,
        "flank_insertions": total_flank,
        "per_barcode": per_bc,
        "fragments_read": n_rows,
    }


def cmd_qc(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    logging.info("scanning %d pairs file(s)", len(paths))

    rows: "list[dict]" = []
    for i, p in enumerate(paths, 1):
        stats = pairs_qc(p, long_range_bp=args.long_range_bp)
        stats["pixel"] = pixel_name(p)
        stats["path"] = str(p)
        rows.append(stats)
        if i % 250 == 0:
            logging.info("  %d/%d pixels", i, len(paths))

    cols = ["pixel", "total_contacts", "cis_contacts", "trans_contacts",
            "long_range_contacts", "cis_fraction", "long_range_ratio",
            "long_range_over_total", "path"]
    tsv = write_tsv(out / "pixel_qc.tsv", rows, cols)

    kept = [r for r in rows if r["total_contacts"] >= args.min_contacts]
    summary: "dict[str, Any]" = {
        "n_pixels": len(rows),
        "n_pixels_passing": len(kept),
        "min_contacts": args.min_contacts,
        "long_range_bp": args.long_range_bp,
        "median_total_contacts": _median([r["total_contacts"] for r in kept]),
        "median_cis_fraction": _median([r["cis_fraction"] for r in kept]),
        "median_long_range_ratio": _median([r["long_range_ratio"] for r in kept]),
        "total_contacts": sum(r["total_contacts"] for r in rows),
    }

    if args.fragments and args.gene_model:
        model = read_gene_model(args.gene_model, upstream=0)
        logging.info("TSS enrichment over %d TSS", len(model))
        tss = tss_enrichment(args.fragments, model,
                             max_rows=args.max_fragments)
        summary["tss_enrichment"] = tss["tss_enrichment"]
        summary["fragments_read"] = tss["fragments_read"]
        per_bc = tss.get("per_barcode") or {}
        if per_bc:
            bc_rows = [{"barcode": b, "tss_enrichment": v}
                       for b, v in sorted(per_bc.items())]
            write_tsv(out / "tss_enrichment_per_barcode.tsv", bc_rows,
                      ["barcode", "tss_enrichment"])
            summary["median_tss_enrichment"] = _median(list(per_bc.values()))

    js = write_json(out / "qc_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Per-pixel QC: {tsv}")
    print(f"Summary: {js}")
    print(f"{len(kept)}/{len(rows)} pixels >= {args.min_contacts} contacts; "
          f"median total={summary['median_total_contacts']:.0f}, "
          f"cis={summary['median_cis_fraction']:.1%}, "
          f"long-range={summary['median_long_range_ratio']:.1%}")
    return 0


# ---------------------------------------------------------------------------
# Gene-level score matrices (GAS / GAD)
# ---------------------------------------------------------------------------

def _score_matrix(
    positions_by_key: "Iterable[tuple[str, str, int]]",
    model: GeneModel,
) -> "dict[str, Counter]":
    """Count (key, chrom, pos) insertions into overlapping gene spans."""
    counts: "dict[str, Counter]" = defaultdict(Counter)
    for key, chrom, pos in positions_by_key:
        for g in model.overlapping(chrom, pos):
            counts[key][g["gene"]] += 1
    return counts


def _write_score_matrix(out: Path, name: str,
                        counts: "dict[str, Counter]",
                        *, min_cells: int = 0) -> "tuple[Path, dict]":
    """Write a gene x pixel matrix as TSV, dropping rare genes."""
    gene_cells: "Counter[str]" = Counter()
    for c in counts.values():
        for g in c:
            gene_cells[g] += 1
    genes = sorted(g for g, n in gene_cells.items() if n >= min_cells)
    keys = sorted(counts)

    path = out / name
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene"] + keys)
        for g in genes:
            w.writerow([g] + [counts[k].get(g, 0) for k in keys])

    total = sum(sum(c.values()) for c in counts.values())
    return path, {
        "n_genes": len(genes),
        "n_pixels": len(keys),
        "total_counts": total,
        "median_counts_per_pixel": _median([sum(c.values()) for c in counts.values()]),
    }


def cmd_gas(args: argparse.Namespace) -> int:
    """Gene activity score matrix from ATAC fragments (SnapATAC2-style)."""
    out = run_dir(args.label)
    model = read_gene_model(args.gene_model, upstream=args.upstream)
    logging.info("gene model: %d genes, promoter pad %d bp",
                 len(model), args.upstream)

    keep: "Optional[set[str]]" = None
    if args.barcodes:
        keep = {b.strip() for b in _open_text(args.barcodes).read().split()
                if b.strip()}
        logging.info("restricting to %d barcodes", len(keep))

    # Translate raw barcodes to AAxBB pixel ids so the GAS matrix shares
    # a key space with GAD, the per-pixel QC table and the cluster labels.
    to_pixel: "Optional[dict[str, str]]" = None
    if args.barcode_a and args.barcode_b:
        to_pixel = barcode_to_pixel_map(args.barcode_a, args.barcode_b,
                                        layout=args.layout)
        logging.info("mapping %d barcodes onto the pixel grid", len(to_pixel))
    elif args.barcode_a or args.barcode_b:
        raise SystemExit("--barcode-a and --barcode-b must be given together")

    n_unmapped = 0

    def gen() -> "Iterator[tuple[str, str, int]]":
        nonlocal n_unmapped
        n = 0
        with _open_text(args.fragments) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                n += 1
                if args.max_fragments and n > args.max_fragments:
                    return
                bc = f[3]
                if keep is not None and bc not in keep:
                    continue
                if to_pixel is not None:
                    mapped = to_pixel.get(bc.strip().upper())
                    if mapped is None:
                        n_unmapped += 1
                        continue
                    bc = mapped
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    continue
                # Both fragment ends are Tn5 insertion sites.
                yield bc, f[0], start
                yield bc, f[0], end - 1

    counts = _score_matrix(gen(), model)
    path, stats = _write_score_matrix(out, "gene_activity_score.tsv", counts,
                                      min_cells=args.min_cells)
    stats["source"] = str(args.fragments)
    stats["upstream"] = args.upstream
    stats["keyed_by"] = "pixel" if to_pixel is not None else "barcode"
    if to_pixel is not None:
        stats["unmapped_fragment_ends"] = n_unmapped
    js = write_json(out / "gas_summary.json", stats)

    print(f"Output dir: {out}")
    print(f"GAS matrix: {path}")
    print(f"Summary: {js}")
    print(f"{stats['n_genes']} genes x {stats['n_pixels']} "
          f"{'pixels' if to_pixel is not None else 'barcodes'}, "
          f"{stats['total_counts']:,} counts")
    if to_pixel is None:
        print("NOTE: columns are raw barcodes. Pass --barcode-a/--barcode-b "
              "to key by AAxBB pixel id instead, so this matrix joins with "
              "the GAD matrix, the QC table and `viz`.")
    elif n_unmapped:
        print(f"{n_unmapped:,} fragment end(s) had an off-whitelist barcode")
    return 0


def cmd_gad(args: argparse.Namespace) -> int:
    """Gene-associated domain score matrix from Hi-C pairs.

    Per the paper's Methods, GAD is computed by converting read pairs to
    a BED-like table where each line is one end of a pair, then running
    the same gene-body counting used for GAS. Gene bodies only by
    default (no promoter pad) — the score measures contact density
    across the gene body (Shen 2022, scGAD).
    """
    out = run_dir(args.label)
    model = read_gene_model(args.gene_model, upstream=args.upstream)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    logging.info("gene model: %d genes; %d pairs file(s)", len(model), len(paths))

    counts: "dict[str, Counter]" = {}
    for i, p in enumerate(paths, 1):
        pid = pixel_name(p)
        header = read_pairs_header(p)
        i_c1, i_p1 = header.index("chrom1"), header.index("pos1")
        i_c2, i_p2 = header.index("chrom2"), header.index("pos2")

        def ends(path=p, ic1=i_c1, ip1=i_p1, ic2=i_c2, ip2=i_p2, key=pid):
            n = 0
            for row in iter_pairs(path):
                if len(row) <= max(ic1, ip1, ic2, ip2):
                    continue
                n += 1
                if args.max_pairs and n > args.max_pairs:
                    return
                try:
                    yield key, row[ic1], int(row[ip1])
                    yield key, row[ic2], int(row[ip2])
                except ValueError:
                    continue

        for key, c in _score_matrix(ends(), model).items():
            counts.setdefault(key, Counter()).update(c)
        if i % 250 == 0:
            logging.info("  %d/%d pixels", i, len(paths))

    path, stats = _write_score_matrix(out, "gene_associated_domain_score.tsv",
                                      counts, min_cells=args.min_cells)
    stats["source"] = str(args.pairs_dir)
    js = write_json(out / "gad_summary.json", stats)

    print(f"Output dir: {out}")
    print(f"GAD matrix: {path}")
    print(f"Summary: {js}")
    print(f"{stats['n_genes']} genes x {stats['n_pixels']} pixels, "
          f"{stats['total_counts']:,} counts")
    return 0


# ---------------------------------------------------------------------------
# Contact matrices
# ---------------------------------------------------------------------------

def build_matrix(
    paths: "list[Path]",
    chrom: str,
    resolution: int,
    *,
    chrom_size: int,
    start: int = 0,
    end: Optional[int] = None,
) -> Any:
    """Bin cis contacts on one chromosome (or window) into a dense matrix."""
    np = _require("numpy")
    lo = start
    hi = min(end if end is not None else chrom_size, chrom_size)
    if hi <= lo:
        raise SystemExit(f"empty window {chrom}:{lo}-{hi}")
    n_bins = int(math.ceil((hi - lo) / resolution))
    mat = np.zeros((n_bins, n_bins), dtype=float)

    for p in paths:
        header = read_pairs_header(p)
        i_c1, i_p1 = header.index("chrom1"), header.index("pos1")
        i_c2, i_p2 = header.index("chrom2"), header.index("pos2")
        for row in iter_pairs(p):
            if len(row) <= max(i_c1, i_p1, i_c2, i_p2):
                continue
            if row[i_c1] != chrom or row[i_c2] != chrom:
                continue
            try:
                p1, p2 = int(row[i_p1]), int(row[i_p2])
            except ValueError:
                continue
            if not (lo <= p1 < hi and lo <= p2 < hi):
                continue
            b1 = (p1 - lo) // resolution
            b2 = (p2 - lo) // resolution
            mat[b1, b2] += 1
            if b1 != b2:
                mat[b2, b1] += 1
    return mat


def bin_coverage(
    paths: "list[Path]",
    resolution: int,
    chrom_sizes: "dict[str, int]",
    *,
    chroms: "Optional[list[str]]" = None,
) -> "tuple[Any, list[tuple[str, int, int]]]":
    """Genome-wide per-bin contact coverage (both ends counted)."""
    np = _require("numpy")
    use = chroms or sorted(chrom_sizes)
    offsets: "dict[str, int]" = {}
    bins: "list[tuple[str, int, int]]" = []
    for c in use:
        offsets[c] = len(bins)
        size = chrom_sizes[c]
        for s in range(0, size, resolution):
            bins.append((c, s, min(s + resolution, size)))
    cov = np.zeros(len(bins), dtype=float)

    for p in paths:
        header = read_pairs_header(p)
        i_c1, i_p1 = header.index("chrom1"), header.index("pos1")
        i_c2, i_p2 = header.index("chrom2"), header.index("pos2")
        for row in iter_pairs(p):
            if len(row) <= max(i_c1, i_p1, i_c2, i_p2):
                continue
            for ic, ip in ((i_c1, i_p1), (i_c2, i_p2)):
                c = row[ic]
                off = offsets.get(c)
                if off is None:
                    continue
                try:
                    pos = int(row[ip])
                except ValueError:
                    continue
                if pos >= chrom_sizes[c]:
                    continue
                cov[off + pos // resolution] += 1
    return cov, bins


def cmd_matrix(args: argparse.Namespace) -> int:
    np = _require("numpy")
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    sizes = _resolve_chrom_sizes(args, paths[0])
    if args.chrom not in sizes:
        raise SystemExit(f"{args.chrom} not in chrom sizes ({len(sizes)} entries)")

    mat = build_matrix(paths, args.chrom, args.resolution,
                       chrom_size=sizes[args.chrom],
                       start=args.start, end=args.end)
    npz = out / f"matrix_{args.chrom}_{args.resolution}.npz"
    np.savez_compressed(npz, matrix=mat, chrom=args.chrom,
                        resolution=args.resolution, start=args.start)
    summary = {
        "chrom": args.chrom, "resolution": args.resolution,
        "n_bins": int(mat.shape[0]), "total_contacts": float(mat.sum() / 2),
        "n_pairs_files": len(paths),
    }
    js = write_json(out / "matrix_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Matrix: {npz}  ({mat.shape[0]} x {mat.shape[0]} bins)")
    print(f"Summary: {js}")
    return 0


def _resolve_chrom_sizes(args: argparse.Namespace, sample_pairs: Path) -> "dict[str, int]":
    """Chrom sizes from --chrom-sizes, else the pairs header."""
    if getattr(args, "chrom_sizes", None):
        return read_chrom_sizes(args.chrom_sizes)
    sizes = read_pairs_header(sample_pairs).chromsizes
    if not sizes:
        raise SystemExit(
            "no chrom sizes: the pairs header carries no '#chromsize:' lines, "
            "so pass --chrom-sizes <file>")
    return sizes


# ---------------------------------------------------------------------------
# impute / compartment
# ---------------------------------------------------------------------------

def _load_math():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from igvfagent import _spatial_hic_math as M  # type: ignore
    except Exception:
        import _spatial_hic_math as M  # type: ignore
    return M


def cmd_impute(args: argparse.Namespace) -> int:
    np = _require("numpy")
    M = _load_math()
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    sizes = _resolve_chrom_sizes(args, paths[0])
    if args.chrom not in sizes:
        raise SystemExit(f"{args.chrom} not in chrom sizes")

    imputed: "dict[str, Any]" = {}
    per_pixel = args.per_pixel
    targets = paths if per_pixel else [paths]

    for i, target in enumerate(targets, 1):
        group = [target] if per_pixel else target
        key = pixel_name(target) if per_pixel else "pseudobulk"
        mat = build_matrix(group, args.chrom, args.resolution,
                           chrom_size=sizes[args.chrom],
                           start=args.start, end=args.end)
        if mat.sum() < args.min_contacts:
            continue
        imputed[key] = M.impute_matrix(
            mat, pad=args.pad, rp=args.restart, tol=args.tol,
            zero_diagonal=args.zero_diagonal)
        if per_pixel and i % 100 == 0:
            logging.info("  imputed %d/%d pixels", i, len(targets))

    if not imputed:
        raise SystemExit(
            f"no pixel reached --min-contacts {args.min_contacts} on "
            f"{args.chrom}; lower the threshold or widen the window")

    npz = out / f"imputed_{args.chrom}_{args.resolution}.npz"
    np.savez_compressed(npz, **imputed)
    summary = {
        "chrom": args.chrom, "resolution": args.resolution,
        "pad": args.pad, "restart": args.restart,
        "n_matrices": len(imputed),
        "mode": "per-pixel" if per_pixel else "pseudobulk",
        "n_bins": int(next(iter(imputed.values())).shape[0]),
    }
    js = write_json(out / "impute_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Imputed matrices: {npz}  ({len(imputed)} matrices)")
    print(f"Summary: {js}")
    return 0


def _gene_density_track(model: GeneModel, chrom: str, n_bins: int,
                        resolution: int) -> Any:
    """Genes per bin — used to orient the compartment eigenvector."""
    np = _require("numpy")
    track = np.zeros(n_bins)
    for g in model.genes:
        if g["chrom"] != chrom:
            continue
        b = g["tss"] // resolution
        if 0 <= b < n_bins:
            track[b] += 1
    return track


def cmd_compartment(args: argparse.Namespace) -> int:
    np = _require("numpy")
    M = _load_math()
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    sizes = _resolve_chrom_sizes(args, paths[0])

    chroms = ([c.strip() for c in args.chroms.split(",") if c.strip()]
              if args.chroms else sorted(sizes))
    model = read_gene_model(args.gene_model, upstream=0) if args.gene_model else None

    rows: "list[dict]" = []
    tracks: "dict[str, Any]" = {}
    for chrom in chroms:
        if chrom not in sizes:
            logging.warning("skipping %s: not in chrom sizes", chrom)
            continue
        mat = build_matrix(paths, chrom, args.resolution,
                           chrom_size=sizes[chrom])
        if mat.sum() == 0:
            logging.warning("skipping %s: no contacts", chrom)
            continue
        if args.impute:
            mat = M.impute_matrix(mat, pad=args.pad, zero_diagonal=True)
        orient = (_gene_density_track(model, chrom, mat.shape[0], args.resolution)
                  if model else None)
        pc1 = M.compartment_pc1(mat, orient_track=orient)
        tracks[chrom] = pc1
        for b, v in enumerate(pc1):
            rows.append({
                "chrom": chrom,
                "start": b * args.resolution,
                "end": min((b + 1) * args.resolution, sizes[chrom]),
                "pc1": "" if v != v else round(float(v), 6),
                "compartment": "" if v != v else ("A" if v > 0 else "B"),
            })
        logging.info("%s: %d bins, %d A / %d B", chrom, len(pc1),
                     int((pc1 > 0).sum()), int((pc1 < 0).sum()))

    if not rows:
        raise SystemExit("no chromosome produced a compartment track")

    tsv = write_tsv(out / "compartments.tsv", rows,
                    ["chrom", "start", "end", "pc1", "compartment"])
    npz = out / "compartment_pc1.npz"
    np.savez_compressed(npz, **tracks)

    finite = [r["pc1"] for r in rows if r["pc1"] != ""]
    summary = {
        "resolution": args.resolution,
        "n_chroms": len(tracks),
        "n_bins": len(rows),
        "n_bins_called": len(finite),
        "frac_A": (sum(1 for r in rows if r["compartment"] == "A") /
                   len(finite)) if finite else 0.0,
        "oriented_by": "gene density" if model else "unoriented (sign arbitrary)",
    }
    js = write_json(out / "compartment_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Compartment track: {tsv}")
    print(f"PC1 arrays: {npz}")
    print(f"Summary: {js}")
    print(f"{len(tracks)} chromosomes, {len(finite)}/{len(rows)} bins called, "
          f"{summary['frac_A']:.1%} A")
    if not model:
        print("NOTE: no --gene-model given, so the A/B sign is arbitrary "
              "per chromosome. Pass one to orient A toward gene-dense bins.")
    return 0


# ---------------------------------------------------------------------------
# cnv
# ---------------------------------------------------------------------------

def _read_bedgraph_track(path: "str | Path",
                         bins: "list[tuple[str, int, int]]") -> Any:
    """Map a bedGraph onto the bin list by midpoint lookup."""
    np = _require("numpy")
    index = {(c, s): i for i, (c, s, _e) in enumerate(bins)}
    if not bins:
        return np.zeros(0)
    res = bins[0][2] - bins[0][1]
    out = np.full(len(bins), np.nan)
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.split()
            if len(f) < 4:
                continue
            try:
                start, end, val = int(f[1]), int(f[2]), float(f[3])
            except ValueError:
                continue
            mid = (start + end) // 2
            i = index.get((f[0], (mid // res) * res))
            if i is not None:
                out[i] = val
    return out


def cmd_cnv(args: argparse.Namespace) -> int:
    np = _require("numpy")
    M = _load_math()
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    sizes = _resolve_chrom_sizes(args, paths[0])
    chroms = ([c.strip() for c in args.chroms.split(",") if c.strip()]
              if args.chroms else sorted(sizes))
    chroms = [c for c in chroms if c in sizes]
    if not chroms:
        raise SystemExit("no requested chromosome is present in the chrom sizes")

    # ── pseudobulk track ──────────────────────────────────────────────
    cov, bins = bin_coverage(paths, args.resolution, sizes, chroms=chroms)
    gc = _read_bedgraph_track(args.gc, bins) if args.gc else None
    mappability = (_read_bedgraph_track(args.mappability, bins)
                   if args.mappability else None)
    ratio = M.cnv_ratio(cov, gc=gc, mappability=mappability,
                        ploidy=args.ploidy, min_coverage=args.min_coverage)

    rows: "list[dict]" = []
    seg = None
    if args.segment:
        seg = M.segment_cnv(ratio, sigma=args.sigma, stay=args.stay)
    for i, (c, s, e) in enumerate(bins):
        rows.append({
            "chrom": c, "start": s, "end": e,
            "coverage": float(cov[i]),
            "cn_ratio": "" if ratio[i] != ratio[i] else round(float(ratio[i]), 4),
            "cn_segment": "" if seg is None else int(seg[i]),
        })
    tsv = write_tsv(out / "cnv_pseudobulk.tsv", rows,
                    ["chrom", "start", "end", "coverage", "cn_ratio", "cn_segment"])

    summary: "dict[str, Any]" = {
        "resolution": args.resolution,
        "ploidy": args.ploidy,
        "n_bins": len(bins),
        "n_bins_called": int(np.isfinite(ratio).sum()),
        "gc_corrected": bool(args.gc),
        "mappability_corrected": bool(args.mappability),
        "segmented": bool(args.segment),
        "n_pairs_files": len(paths),
    }

    # ── per-pixel matrix (the paper's 5 Mb single-pixel CN clustering) ─
    per_pixel_path = None
    if args.per_pixel and len(paths) > 1:
        logging.info("per-pixel CN at %d bp over %d pixels",
                     args.resolution, len(paths))
        mat = np.full((len(paths), len(bins)), np.nan)
        names: "list[str]" = []
        for i, p in enumerate(paths):
            names.append(pixel_name(p))
            c_i, _ = bin_coverage([p], args.resolution, sizes, chroms=chroms)
            if c_i.sum() < args.min_pixel_contacts:
                continue
            mat[i] = M.cnv_ratio(c_i, gc=gc, mappability=mappability,
                                 ploidy=args.ploidy,
                                 min_coverage=args.min_coverage)
            if (i + 1) % 250 == 0:
                logging.info("  %d/%d pixels", i + 1, len(paths))

        if args.smooth:
            keep = np.isfinite(mat).any(axis=1)
            if keep.sum() >= 2:
                mat[keep] = M.magic_smooth(mat[keep], t=args.smooth_t,
                                           k=args.smooth_k)

        per_pixel_path = out / "cnv_per_pixel.tsv"
        with per_pixel_path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pixel"] + [f"{c}:{s}" for c, s, _e in bins])
            for i, name in enumerate(names):
                w.writerow([name] + ["" if v != v else round(float(v), 4)
                                     for v in mat[i]])
        summary["n_pixels"] = len(names)
        summary["smoothed"] = bool(args.smooth)

    js = write_json(out / "cnv_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Pseudobulk CN: {tsv}")
    if per_pixel_path:
        print(f"Per-pixel CN: {per_pixel_path}")
    print(f"Summary: {js}")
    print(f"{summary['n_bins_called']}/{len(bins)} bins called at "
          f"{args.resolution:,} bp")
    return 0


# ---------------------------------------------------------------------------
# loops
# ---------------------------------------------------------------------------

def read_bedpe(path: "str | Path") -> "list[dict]":
    """Read loop anchors from a BEDPE."""
    loops: "list[dict]" = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                loops.append({
                    "chrom1": f[0], "start1": int(f[1]), "end1": int(f[2]),
                    "chrom2": f[3], "start2": int(f[4]), "end2": int(f[5]),
                    "name": f[6] if len(f) > 6 else "",
                })
            except ValueError:
                continue
    if not loops:
        raise SystemExit(f"no BEDPE records parsed from {path}")
    return loops


def read_cluster_labels(path: "str | Path") -> "dict[str, str]":
    """Two-column pixel -> cluster TSV/CSV (header optional)."""
    labels: "dict[str, str]" = {}
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line)
            if len(parts) < 2:
                continue
            if parts[0].lower() in ("pixel", "barcode", "cell"):
                continue
            labels[parts[0]] = parts[1]
    if not labels:
        raise SystemExit(f"no cluster labels parsed from {path}")
    return labels


def cmd_loops(args: argparse.Namespace) -> int:
    np = _require("numpy")
    M = _load_math()
    out = run_dir(args.label)
    paths = discover_pairs(args.pairs_dir, pattern=getattr(args, 'pairs_glob', None))
    loops = read_bedpe(args.bedpe)
    res = args.resolution
    logging.info("%d loops, %d pixels, %d bp bins", len(loops), len(paths), res)

    # Only cis loops can be quantified from a per-chromosome matrix.
    cis = [l for l in loops if l["chrom1"] == l["chrom2"]]
    if len(cis) < len(loops):
        logging.warning("dropping %d trans loop(s)", len(loops) - len(cis))
    if not cis:
        raise SystemExit("no cis loops in the BEDPE")

    # ── per-pixel loop strength ───────────────────────────────────────
    strength = np.zeros((len(paths), len(cis)))
    pixels: "list[str]" = []
    for pi, p in enumerate(paths):
        pixels.append(pixel_name(p))
        header = read_pairs_header(p)
        i_c1, i_p1 = header.index("chrom1"), header.index("pos1")
        i_c2, i_p2 = header.index("chrom2"), header.index("pos2")
        # Bin every contact once, then look loops up — cheaper than
        # re-scanning the file per loop.
        binned: "Counter[tuple[str, int, int]]" = Counter()
        for row in iter_pairs(p):
            if len(row) <= max(i_c1, i_p1, i_c2, i_p2):
                continue
            if row[i_c1] != row[i_c2]:
                continue
            try:
                b1, b2 = int(row[i_p1]) // res, int(row[i_p2]) // res
            except ValueError:
                continue
            if b1 > b2:
                b1, b2 = b2, b1
            binned[(row[i_c1], b1, b2)] += 1
        for li, l in enumerate(cis):
            b1 = (l["start1"] + l["end1"]) // 2 // res
            b2 = (l["start2"] + l["end2"]) // 2 // res
            if b1 > b2:
                b1, b2 = b2, b1
            total = 0
            for d1 in range(-args.window, args.window + 1):
                for d2 in range(-args.window, args.window + 1):
                    total += binned.get((l["chrom1"], b1 + d1, b2 + d2), 0)
            strength[pi, li] = total
        if (pi + 1) % 100 == 0:
            logging.info("  %d/%d pixels", pi + 1, len(paths))

    mat_path = out / "loop_strength_per_pixel.tsv"
    with mat_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pixel"] + [
            f"{l['chrom1']}:{l['start1']}-{l['end1']}_{l['start2']}-{l['end2']}"
            for l in cis])
        for i, name in enumerate(pixels):
            w.writerow([name] + [int(v) for v in strength[i]])

    summary: "dict[str, Any]" = {
        "n_loops": len(cis),
        "n_pixels": len(pixels),
        "resolution": res,
        "window_bins": args.window,
        "total_loop_contacts": float(strength.sum()),
    }

    # ── cluster-specific loops via one-way ANOVA ──────────────────────
    diff_path = None
    if args.clusters:
        labels = read_cluster_labels(args.clusters)
        groups: "dict[str, list[int]]" = defaultdict(list)
        for i, name in enumerate(pixels):
            lab = labels.get(name)
            if lab is not None:
                groups[lab].append(i)
        missing = len(pixels) - sum(len(v) for v in groups.values())
        if missing:
            logging.warning("%d pixel(s) had no cluster label", missing)
        if len(groups) < 2:
            raise SystemExit(
                f"need >= 2 clusters with pixels, got {len(groups)}")

        rows: "list[dict]" = []
        for li, l in enumerate(cis):
            arrays = [strength[idx, li] for idx in groups.values()]
            f_stat, p_val = M.one_way_anova(arrays)
            means = {k: float(strength[idx, li].mean())
                     for k, idx in groups.items()}
            top = max(means, key=means.get) if means else ""
            rows.append({
                "chrom": l["chrom1"],
                "start1": l["start1"], "end1": l["end1"],
                "start2": l["start2"], "end2": l["end2"],
                "name": l["name"],
                "F": "" if f_stat != f_stat else round(f_stat, 4),
                "p_value": "" if p_val != p_val else f"{p_val:.6g}",
                "specific": ("yes" if (p_val == p_val and p_val < args.alpha)
                             else "no"),
                "top_cluster": top,
                **{f"mean_{k}": round(v, 3) for k, v in sorted(means.items())},
            })
        diff_path = write_tsv(out / "cluster_specific_loops.tsv", rows)
        n_sig = sum(1 for r in rows if r["specific"] == "yes")
        summary["n_clusters"] = len(groups)
        summary["alpha"] = args.alpha
        summary["n_cluster_specific_loops"] = n_sig
        by_cluster = Counter(r["top_cluster"] for r in rows
                             if r["specific"] == "yes")
        summary["specific_by_cluster"] = dict(by_cluster)

    # ── APA pileup ────────────────────────────────────────────────────
    apa_path = None
    if args.apa:
        sizes = _resolve_chrom_sizes(args, paths[0])
        mats: "dict[str, Any]" = {}
        for chrom in sorted({l["chrom1"] for l in cis}):
            if chrom not in sizes:
                continue
            mats[chrom] = build_matrix(paths, chrom, res,
                                       chrom_size=sizes[chrom])
        anchors = [(l["chrom1"],
                    (l["start1"] + l["end1"]) // 2 // res,
                    (l["start2"] + l["end2"]) // 2 // res) for l in cis]
        pile = M.apa_matrix(mats, anchors, flank=args.apa_flank)
        apa_path = out / "apa_matrix.tsv"
        np.savetxt(apa_path, pile, delimiter="\t", fmt="%.6g")
        if np.isfinite(pile).any():
            c = pile.shape[0] // 2
            corners = [pile[0, 0], pile[0, -1], pile[-1, 0], pile[-1, -1]]
            bg = float(np.nanmean(corners))
            summary["apa_centre"] = float(pile[c, c])
            summary["apa_corner_mean"] = bg
            # A zero corner background makes the ratio undefined, not
            # infinite-and-therefore-great: it usually means the pileup
            # is too sparse to have a background at all. Report the
            # centre and say so rather than emitting a NaN ratio.
            summary["apa_enrichment"] = (float(pile[c, c] / bg) if bg > 0
                                         else None)

    js = write_json(out / "loops_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Loop strength: {mat_path}")
    if diff_path:
        print(f"Cluster-specific loops: {diff_path}  "
              f"({summary['n_cluster_specific_loops']}/{len(cis)} at "
              f"p<{args.alpha})")
    if apa_path:
        print(f"APA pileup: {apa_path}")
        if summary.get("apa_enrichment") is not None:
            print(f"APA centre/corner enrichment: "
                  f"{summary['apa_enrichment']:.2f}x")
        elif "apa_centre" in summary:
            print(f"APA centre={summary['apa_centre']:.3g}; corner "
                  f"background is zero, so no enrichment ratio is defined "
                  f"(the pileup is too sparse for a background estimate)")
    print(f"Summary: {js}")
    return 0


# ---------------------------------------------------------------------------
# viz
# ---------------------------------------------------------------------------

def _pixel_grid_coords(name: str,
                       positions: "Optional[dict[str, tuple[int, int, int]]]"
                       ) -> "Optional[tuple[int, int]]":
    """(row, col) for a pixel, from a positions file or an ``AAxBB`` id."""
    if positions and name in positions:
        _in_tissue, row, col = positions[name]
        return row, col
    m = re.match(r"^(\d+)x(\d+)$", name)
    if m:
        return int(m.group(2)), int(m.group(1))
    return None


def cmd_viz(args: argparse.Namespace) -> int:
    np = _require("numpy")
    matplotlib = _require("matplotlib")
    # Headless by default: this runs on servers and inside the container.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = run_dir(args.label)
    positions = (read_spatial_positions(args.positions)
                 if args.positions else None)

    # ── load the value table ──────────────────────────────────────────
    with _open_text(args.table) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"no data rows in {args.table}")

    key_col = 0
    if args.column not in header:
        # Gene-by-pixel matrices are transposed relative to per-pixel
        # tables: the requested name is a row, not a column.
        gene_rows = {r[0]: r for r in rows}
        if args.column in gene_rows:
            pixels = header[1:]
            values = [float(v) if v not in ("", "NA") else float("nan")
                      for v in gene_rows[args.column][1:]]
        else:
            raise SystemExit(
                f"'{args.column}' is neither a column nor a row of "
                f"{args.table}. Columns: {header[:8]}...")
    else:
        vi = header.index(args.column)
        pixels = [r[key_col] for r in rows]
        values = [float(r[vi]) if len(r) > vi and r[vi] not in ("", "NA")
                  else float("nan") for r in rows]

    grid = np.full((GRID_N, GRID_N), np.nan)
    placed = 0
    for name, val in zip(pixels, values):
        rc = _pixel_grid_coords(name, positions)
        if rc is None:
            continue
        r, c = rc
        # Positions files are 0- or 1-based depending on the exporter;
        # accept both by clamping a 1-based index down.
        if 1 <= r <= GRID_N and 1 <= c <= GRID_N:
            r, c = r - 1, c - 1
        if 0 <= r < GRID_N and 0 <= c < GRID_N:
            grid[r, c] = val
            placed += 1
    if placed == 0:
        raise SystemExit(
            "no pixel could be placed on the grid. Pass --positions with an "
            "AtlasXBrowser tissue_positions CSV, or name pixels 'AAxBB'.")

    finite = grid[np.isfinite(grid)]
    vmin = args.vmin if args.vmin is not None else float(np.percentile(finite, 2))
    vmax = args.vmax if args.vmax is not None else float(np.percentile(finite, 98))

    fig, ax = plt.subplots(figsize=(6.4, 5.6), dpi=200)
    im = ax.imshow(grid, cmap=args.cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest", origin="upper")
    ax.set_title(args.title or f"{args.column}", fontsize=11)
    ax.set_xlabel("barcode A channel")
    ax.set_ylabel("barcode B channel")
    ax.set_xticks([0, GRID_N // 2, GRID_N - 1])
    ax.set_xticklabels([1, GRID_N // 2 + 1, GRID_N])
    ax.set_yticks([0, GRID_N // 2, GRID_N - 1])
    ax.set_yticklabels([1, GRID_N // 2 + 1, GRID_N])
    fig.colorbar(im, ax=ax, shrink=0.82, label=args.column)
    fig.tight_layout()

    png = out / "Plots" / f"spatial_{safe_label(args.column)}.png"
    svg = png.with_suffix(".svg")
    fig.savefig(png)
    fig.savefig(svg)
    plt.close(fig)

    summary = {
        "column": args.column,
        "table": str(args.table),
        "pixels_placed": placed,
        "pixels_in_table": len(pixels),
        "vmin": vmin, "vmax": vmax,
        "value_median": float(np.nanmedian(grid)),
    }
    js = write_json(out / "viz_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Spatial map: {png}")
    print(f"Spatial map (SVG): {svg}")
    print(f"Summary: {js}")
    print(f"Placed {placed}/{len(pixels)} pixels on the {GRID_N}x{GRID_N} grid")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

PLAYBOOK = """# Spatial-ATAC-Hi-C skill

Spatially resolved co-profiling of 3D genome organisation and chromatin
accessibility, after Wang P., Wang J., Wang Q., Youngblood M. W. et al.,
*Nature Methods* 2026 (doi:10.1038/s41592-026-03217-4), GEO **GSE307620**.

The assay puts a 50 x 50 microfluidic barcode grid over a tissue
section: 2,500 spatial pixels, each with a Hi-C contact set *and* an
ATAC fragment set from the same molecules.

## What this skill consumes

Processed deposits only — the same boundary `share_seq` draws:

| Input | Flag | Notes |
|---|---|---|
| Contact table | `--pairs-dir` | 4DN `.pairs`, **or** GEO's tabix-indexed `*.hic.fragments.sorted.header.tsv.gz`; one file or a directory of per-pixel ones |
| ATAC fragments | `--fragments` | `chrom start end barcode count` |
| Spatial positions | `--positions` | AtlasXBrowser `*_tissue_positions_list.csv` |
| Gene model | `--gene-model` | GTF/GFF or BED6 |
| Chrom sizes | `--chrom-sizes` | optional if the pairs header has `#chromsize:` |

**Contact-file shapes.** A 4DN header (`#columns:` / `#chromsize:`), a
bare column-header row, and a headerless TSV all parse identically, and
column-name synonyms (`chr1`/`chrom1`, `position1`/`pos1`, ...) are
normalised. Pixel ids are taken from an `AAxBB` token anywhere in the
filename, so `12x47.hic.fragments.sorted.header.tsv.gz` keys as `12x47`
and joins with the cluster labels and the tissue grid. Use
`--pairs-glob` if a directory holds TSVs you do *not* want swept in.

What GSE307620 actually deposits, per sample (inside `GSE307620_RAW.tar`,
9.9 GB):

```
GSM9228169_MBR6merge.fragments.tsv.gz                    <- --fragments
GSM9228169_MBR6merge.hic.fragments.sorted.header.tsv.gz  <- --pairs-dir
GSM9228169_MouseBrainR6_tissue_positions_list.csv.gz     <- --positions
GSM9228169_MouseBrainR6_scalefactors_json.json.gz
GSM9228169_MouseBrainR6_tissue_{hires,lowres}_image.png.gz
```

Read alignment is **not** reimplemented. It is the one upstream step with
no tractable clean-room form, and the tools involved (Cell Ranger ATAC;
runHiC, GPL-3.0) are excluded from this Apache-2 codebase by policy.
Bring aligned pairs/fragments — GSE307620 publishes them.

## Typical run

```bash
# 0. See what the series deposits, then fetch the parts you want.
igvfagent spatial-hic pull-geo --gse GSE307620 \\
    --download 'pairs|fragments|positions'

# 1. If you have one barcoded pairs file, split it into the pixel grid.
igvfagent spatial-hic pixel-demux --pairs sample.pairs.gz \\
    --barcode-a barcodes_A.txt --barcode-b barcodes_B.txt \\
    --layout BA --label mouse_R6

# 2. Per-pixel QC. The paper reports medians of 25,343-58,403 total
#    contacts, 88.1-90.3% cis, 24-33.3% long-range (>=10 kb over cis).
igvfagent spatial-hic qc --pairs-dir <run>/pixels \\
    --fragments fragments.tsv.gz --gene-model gencode.gtf.gz \\
    --label mouse_R6

# 3. Gene-level scores. GAS from ATAC (promoter + body), GAD from Hi-C
#    pair ends over the gene body.
igvfagent spatial-hic gas --fragments fragments.tsv.gz \\
    --gene-model gencode.gtf.gz --label mouse_R6
igvfagent spatial-hic gad --pairs-dir <run>/pixels \\
    --gene-model gencode.gtf.gz --label mouse_R6

# 4. Imputation, compartments, copy number.
igvfagent spatial-hic impute --pairs-dir <run>/pixels --chrom chr2 \\
    --resolution 100000 --pad 1 --per-pixel --label mouse_R6
igvfagent spatial-hic compartment --pairs-dir <run>/pixels \\
    --resolution 100000 --gene-model gencode.gtf.gz --label mouse_R6
igvfagent spatial-hic cnv --pairs-dir <run>/pixels \\
    --resolution 5000000 --per-pixel --smooth --label tumour

# 5. Loops: quantify at known anchors, test across clusters, pile up.
igvfagent spatial-hic loops --pairs-dir <run>/pixels --bedpe loops.bedpe \\
    --clusters clusters.tsv --apa --label mouse_R6

# 6. Put any per-pixel value back into tissue space.
igvfagent spatial-hic viz --table <run>/gene_activity_score.tsv \\
    --column Satb2 --positions tissue_positions.csv --label mouse_R6
```

## Method notes

**Pixel demux.** The barcode is the two microfluidic rounds
concatenated. Upstream `bcsplit.py` writes barcode-B then barcode-A
(Read2 offsets 22:30 and 60:68), which is the `--layout BA` default.
Pass `--layout AB` if yours is the other way round; the `assigned`
fraction in `demux_summary.json` tells you immediately if you guessed
wrong.

**QC.** `long_range_ratio` uses cis as the denominator, matching the
paper's Fig. 1f ("fraction of long-range contacts over the unique
intra-chromosomal contacts"). `long_range_over_total` is reported
alongside so neither reading is ambiguous.

**TSS enrichment.** ArchR's definition, as the paper specifies:
mean insertion density in a 50 bp window on the TSS over the density in
the +/-1,900-2,000 bp flanks. Densities, not raw counts — the windows
differ in width by 4x.

**Imputation.** scHiCluster's convolution + random-walk-with-restart.
Note the restart term puts `rp` of the walk's mass on the diagonal
(over half of it at the default `rp=0.5`); pass `--zero-diagonal` when
whatever reads the matrix next cares about off-diagonal structure.
Imputation earns its keep at low depth: on a decay toy at 0.3x depth,
off-diagonal correlation with truth goes 0.43 -> 0.85.

**Compartments.** Observed/expected, Pearson correlation, leading
eigenvector. The eigenvector's sign is arbitrary, so pass
`--gene-model` to orient A toward gene-dense bins — otherwise A and B
can flip between chromosomes.

**Copy number.** Binned coverage scaled so the genome median is
`--ploidy` (2 by default; the paper treats every sample as diploid and
reports CN/2). GC and mappability corrections are a linear fit, not
NeoLoopFinder's Poisson GLM — concordant, not bit-identical. `--smooth`
applies MAGIC diffusion, which is what the paper uses before plotting
single-pixel CN in tissue space.

**Loops.** This quantifies loops at anchors you supply and tests them
across clusters (one-way ANOVA, `p < 0.05`, as in the paper). It does
**not** call loops de novo — Peakachu is a trained model, so bring its
BEDPE (or any other caller's).

## Outputs

Everything lands under `Docs/SpatialATACHiC/<timestamp>_<label>/`:
per-pixel TSVs, `.npz` matrix bundles, a `*_summary.json` per command,
and PNG+SVG figures under `Plots/`.

## Provenance and licensing

Apache-2.0. Algorithms are reimplemented from published descriptions;
no upstream source is copied, imported or vendored.

| Capability | Reference implementation | License | Approach |
|---|---|---|---|
| Assay + preprocessing | [wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C) | MIT | clean-room; barcode offsets and linker literals follow the published protocol |
| Contact-matrix build | [XiaoTaoWang/HiC_pipeline](https://github.com/XiaoTaoWang/HiC_pipeline) (runHiC) | **GPL-3.0** | external tool only — never imported. This skill starts from its `.pairs` output |
| Adapter trimming | [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) | **GPL-3.0** | external tool only — upstream of this skill's inputs |
| Imputation | scHiCluster (Zhou 2019 PNAS) | MIT | clean-room convolution + RWR |
| Compartments | cooltools `eigs-cis` | MIT | clean-room O/E + correlation eigenvector |
| CNV / segmentation | NeoLoopFinder `calculate-cnv`, `segment-cnv` | MIT | clean-room; linear bias fit in place of a Poisson GLM |
| Gene activity | SnapATAC2 `make_gene_matrix` | MIT | clean-room promoter+body insertion counting |
| Spatial smoothing | MAGIC (van Dijk 2018) | GPL-2.0 | clean-room diffusion — not imported |
| Loop calling | Peakachu (Salameh 2020) | MIT | **not** reimplemented; bring a BEDPE |

Runtime dependencies: `numpy`, `scipy`, `matplotlib` (all in the
`analysis` extra). No GPL runtime dependencies.
"""


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "SPATIAL_ATAC_HIC_SKILLS.md"
    path.write_text(PLAYBOOK, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_pairs_source(p: argparse.ArgumentParser) -> None:
    p.add_argument("--pairs-dir", required=True,
                   help="A contact file, or a directory of per-pixel ones. "
                        "Accepts 4DN .pairs and the tabix-indexed contact "
                        "TSVs GEO deposits.")
    p.add_argument("--chrom-sizes", default=None,
                   help="chrom.sizes file (optional if the pairs header "
                        "carries #chromsize: lines).")
    p.add_argument("--pairs-glob", default=None,
                   help="Override the contact-file match, e.g. "
                        "'*.hic.fragments.sorted.header.tsv.gz'.")


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(
        prog="spatial-hic",
        description="Spatial-ATAC-Hi-C: spatially resolved 3D genome + "
                    "chromatin accessibility (Wang et al. 2026).")
    sub = parser.add_subparsers(dest="command", required=True)

    # pull-geo
    p = sub.add_parser("pull-geo", help="Discover / fetch GEO deposits.")
    p.add_argument("--gse", default="GSE307620")
    p.add_argument("--download", default=None,
                   help="Regex; matching supplementary files are downloaded.")
    p.add_argument("--dest", default=None)
    p.add_argument("--label", default="spatial_hic_geo")

    # pixel-demux
    p = sub.add_parser("pixel-demux",
                       help="Split a barcoded pairs file into the pixel grid.")
    p.add_argument("--pairs", required=True)
    p.add_argument("--barcode-a", required=True)
    p.add_argument("--barcode-b", required=True)
    p.add_argument("--layout", choices=["BA", "AB"], default="BA",
                   help="Concatenation order in the barcode string "
                        "(default BA, matching upstream bcsplit.py).")
    p.add_argument("--barcode-field", type=int, default=None,
                   help="1-based pairs column holding the barcode.")
    p.add_argument("--readid-sep", default=":")
    p.add_argument("--readid-part", type=int, default=-1,
                   help="Index into the split readID (default last).")
    p.add_argument("--max-pairs", type=int, default=0, help="0 = all.")
    p.add_argument("--label", default="spatial_hic_demux")

    # qc
    p = sub.add_parser("qc", help="Per-pixel contact QC + TSS enrichment.")
    _add_pairs_source(p)
    p.add_argument("--fragments", default=None)
    p.add_argument("--gene-model", default=None)
    p.add_argument("--long-range-bp", type=int, default=LONG_RANGE_BP)
    p.add_argument("--min-contacts", type=int, default=1000)
    p.add_argument("--max-fragments", type=int, default=0)
    p.add_argument("--label", default="spatial_hic_qc")

    # gas
    p = sub.add_parser("gas", help="Gene activity score from ATAC fragments.")
    p.add_argument("--fragments", required=True)
    p.add_argument("--gene-model", required=True)
    p.add_argument("--upstream", type=int, default=2000)
    p.add_argument("--barcodes", default=None,
                   help="Optional whitelist of barcodes to keep.")
    p.add_argument("--barcode-a", default=None,
                   help="With --barcode-b, key columns by AAxBB pixel id "
                        "instead of raw barcode, so the matrix joins with "
                        "the GAD / QC tables.")
    p.add_argument("--barcode-b", default=None)
    p.add_argument("--layout", choices=["BA", "AB"], default="BA")
    p.add_argument("--min-cells", type=int, default=5)
    p.add_argument("--max-fragments", type=int, default=0)
    p.add_argument("--label", default="spatial_hic_gas")

    # gad
    p = sub.add_parser("gad", help="Gene-associated domain score from pairs.")
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--gene-model", required=True)
    p.add_argument("--upstream", type=int, default=0)
    p.add_argument("--min-cells", type=int, default=5)
    p.add_argument("--max-pairs", type=int, default=0)
    p.add_argument("--label", default="spatial_hic_gad")

    # matrix
    p = sub.add_parser("matrix", help="Bin pairs into a contact matrix.")
    _add_pairs_source(p)
    p.add_argument("--chrom", required=True)
    p.add_argument("--resolution", type=int, default=100000)
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=None)
    p.add_argument("--label", default="spatial_hic_matrix")

    # impute
    p = sub.add_parser("impute", help="scHiCluster-style imputation.")
    _add_pairs_source(p)
    p.add_argument("--chrom", required=True)
    p.add_argument("--resolution", type=int, default=100000)
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=None)
    p.add_argument("--pad", type=int, default=1)
    p.add_argument("--restart", type=float, default=0.5)
    p.add_argument("--tol", type=float, default=0.01)
    p.add_argument("--per-pixel", action="store_true")
    p.add_argument("--zero-diagonal", action="store_true")
    p.add_argument("--min-contacts", type=int, default=1)
    p.add_argument("--label", default="spatial_hic_impute")

    # compartment
    p = sub.add_parser("compartment", help="A/B compartment PC1.")
    _add_pairs_source(p)
    p.add_argument("--resolution", type=int, default=100000)
    p.add_argument("--chroms", default=None, help="Comma list (default all).")
    p.add_argument("--gene-model", default=None,
                   help="Orients A toward gene-dense bins.")
    p.add_argument("--impute", action="store_true")
    p.add_argument("--pad", type=int, default=1)
    p.add_argument("--label", default="spatial_hic_compartment")

    # cnv
    p = sub.add_parser("cnv", help="Copy-number ratio from contact coverage.")
    _add_pairs_source(p)
    p.add_argument("--resolution", type=int, default=5000000)
    p.add_argument("--chroms", default=None)
    p.add_argument("--gc", default=None, help="bedGraph of GC content.")
    p.add_argument("--mappability", default=None, help="bedGraph.")
    p.add_argument("--ploidy", type=float, default=2.0)
    p.add_argument("--min-coverage", type=float, default=0.0)
    p.add_argument("--segment", action="store_true", help="HMM segmentation.")
    p.add_argument("--sigma", type=float, default=0.5)
    p.add_argument("--stay", type=float, default=0.99)
    p.add_argument("--per-pixel", action="store_true")
    p.add_argument("--min-pixel-contacts", type=int, default=100)
    p.add_argument("--smooth", action="store_true", help="MAGIC diffusion.")
    p.add_argument("--smooth-t", type=int, default=3)
    p.add_argument("--smooth-k", type=int, default=10)
    p.add_argument("--label", default="spatial_hic_cnv")

    # loops
    p = sub.add_parser("loops", help="Quantify loops, test across clusters.")
    _add_pairs_source(p)
    p.add_argument("--bedpe", required=True)
    p.add_argument("--resolution", type=int, default=10000)
    p.add_argument("--window", type=int, default=1,
                   help="Half-width in bins summed around each anchor.")
    p.add_argument("--clusters", default=None,
                   help="pixel<TAB>cluster; enables the ANOVA test.")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--apa", action="store_true")
    p.add_argument("--apa-flank", type=int, default=5)
    p.add_argument("--label", default="spatial_hic_loops")

    # viz
    p = sub.add_parser("viz", help="Render per-pixel values in tissue space.")
    p.add_argument("--table", required=True)
    p.add_argument("--column", required=True,
                   help="A column name, or a gene (row) in a score matrix.")
    p.add_argument("--positions", default=None)
    p.add_argument("--cmap", default="magma")
    p.add_argument("--vmin", type=float, default=None)
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument("--title", default=None)
    p.add_argument("--label", default="spatial_hic_viz")

    sub.add_parser("write-playbook", help="Emit the skill's markdown playbook.")

    args = parser.parse_args(argv)

    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0

    setup_logging()
    handlers = {
        "pull-geo": cmd_pull_geo,
        "pixel-demux": cmd_pixel_demux,
        "qc": cmd_qc,
        "gas": cmd_gas,
        "gad": cmd_gad,
        "matrix": cmd_matrix,
        "impute": cmd_impute,
        "compartment": cmd_compartment,
        "cnv": cmd_cnv,
        "loops": cmd_loops,
        "viz": cmd_viz,
    }
    fn = handlers.get(args.command)
    if fn is None:
        return 2
    return fn(args)


if __name__ == "__main__":
    sys.exit(main())
