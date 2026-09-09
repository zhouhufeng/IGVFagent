#!/usr/bin/env python3
"""Raw-read processing: an IGVF accession straight through to an analysed matrix.

Asking IGVFagent to "analyse this raw data IGVFDS…" used to produce a
*description* of the dataset and a suggestion of what one would normally run,
because the pieces to actually do it were never connected: the Portal query
lived in one skill, the downloader in another, the single-cell pipeline
wanted a matrix nobody had produced, and nothing in the system turned reads
into that matrix at all.

This module is the missing chain. Given an accession it will

1. resolve the FileSet (authenticated, so unreleased data is visible),
2. inventory the files and decide the cheapest honest route,
3. prefer an existing count matrix -- on the set itself, or on an
   AnalysisSet derived from it -- because aligning reads to reproduce a
   matrix somebody already published is hours of compute for no new
   information,
4. otherwise align the reads with kb-python (kallisto | bustools), and
5. hand the resulting matrix to the single-cell pipeline.

Route selection is reported before anything large is downloaded, and
`plan` does only steps 1-2 so the cost is known before it is spent.

Honesty about what "aligned" means here: kallisto|bustools quantifies
against a transcriptome index. That is the right tool for the scRNA-seq
readout of a CRISPR screen, and `kb count --workflow kite` additionally
assigns feature barcodes (sgRNAs) when a guide library is supplied. It is
not a substitute for a genome aligner where one is genuinely required
(structural work, allele-specific pileups); `plan` says so rather than
producing a matrix that quietly answers a different question.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlsplit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint          # noqa: E402
from _credentials import portal_credentials as _portal_credentials  # noqa: E402
import _assays as _A                                          # noqa: E402

IGVF_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
RUN_DIR = ROOT / "Data" / "RawPipeline"
REF_DIR = ROOT / "Data" / "References" / "kb"
REPORT_DIR = ROOT / "Docs" / "RawPipeline"
LOG_DIR = ROOT / "Docs" / "Logs"

# Matrix-ish content the single-cell pipeline can already read.
_MATRIX_FORMATS = {"h5ad", "h5", "mtx", "tar", "csv", "tsv"}
_MATRIX_HINTS = ("matrix", "count")

# R1 length -> 10x chemistry. v2 is a 16 bp barcode + 10 bp UMI; v3 adds two
# UMI bases. v4/GEM-X is also 28, so 28 stays "v3" unless told otherwise --
# the count matrices are equivalent for our purposes and mislabelling the
# on-list is the only real risk.
_TECH_BY_R1_LEN = {26: "10XV2", 28: "10XV3"}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"raw_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


# ─── HTTP ───────────────────────────────────────────────────────────────────

def _headers() -> "dict[str, str]":
    h = {"Accept": "application/json,*/*", "User-Agent": "IGVFdataAgent/0.1 raw-pipeline"}
    creds = _portal_credentials()
    if creds:
        tok = base64.b64encode(f"{creds[0]}:{creds[1]}".encode()).decode()
        h["Authorization"] = f"Basic {tok}"
    elif os.environ.get("IGVF_PORTAL_COOKIE"):
        h["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return h


class _DropAuthOnCrossHostRedirect(urllib.request.HTTPRedirectHandler):
    """Strip credentials when a download redirects to object storage.

    A Portal `@@download` URL answers 302 to a pre-signed S3 URL that
    carries its own signature in the query string. Forwarding our
    `Authorization: Basic` header to S3 makes it reject the request with
    HTTP 400 -- two credentials for one request -- so every authenticated
    file download fails while the same URL works anonymously. Dropping the
    header on a host change is what makes authenticated downloads work at
    all.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is None:
            return None
        if urlsplit(newurl).netloc != urlsplit(req.full_url).netloc:
            for name in ("Authorization", "Cookie"):
                new.headers.pop(name, None)
                new.headers.pop(name.lower(), None)
                new.unredirected_hdrs.pop(name, None)
        return new


_OPENER = urllib.request.build_opener(_DropAuthOnCrossHostRedirect)


def portal_json(path: str) -> "tuple[int, Any]":
    url = path if path.startswith("http") else f"{IGVF_API_BASE}{path}"
    req = urllib.request.Request(url, headers=_headers())
    try:
        with _OPENER.open(req, timeout=90) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, None
    except Exception as e:                                   # noqa: BLE001
        logging.warning("request failed: %s", e)
        return 0, None


def portal_download(href: str, dest: Path, on_bytes=None) -> Path:
    """Download one Portal file, following the S3 redirect safely.

    ``on_bytes(n)`` is called with each chunk size so a caller can report
    progress while a multi-GB transfer is in flight.
    """
    url = href if href.startswith("http") else f"{IGVF_API_BASE}{href}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=_headers())
    chunk = 4 * 1024 * 1024
    with _OPENER.open(req, timeout=600) as r, tmp.open("wb") as fh:
        while True:
            buf = r.read(chunk)
            if not buf:
                break
            fh.write(buf)
            if on_bytes:
                on_bytes(len(buf))
    tmp.replace(dest)
    return dest


def read_text_maybe_gzip(path: Path) -> str:
    """Decode a file that may or may not actually be gzipped.

    The Portal serves some `.yaml.gz` configuration files uncompressed, so
    trusting the extension raises BadGzipFile on a perfectly good seqspec.
    Sniff the magic bytes instead.
    """
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode(errors="replace")
    return raw.decode(errors="replace")


# ─── Portal inventory ───────────────────────────────────────────────────────

def resolve_file_set(accession: str) -> "Optional[dict]":
    status, data = portal_json(
        f"/search/?type=FileSet&accession={urllib.parse.quote(accession)}&format=json")
    rows = (data or {}).get("@graph") or []
    for row in rows:
        if str(row.get("accession", "")).upper() == accession.upper():
            return row
    if status == 403:
        logging.warning("%s returned 403 — credentials do not reach it", accession)
    return None


def list_files(accession: str) -> "list[dict]":
    """Every File on a FileSet, hydrated enough to download and pair."""
    status, data = portal_json(
        f"/search/?type=File&file_set.accession={urllib.parse.quote(accession)}"
        f"&limit=500&format=json")
    rows = (data or {}).get("@graph") or []
    out = []
    for row in rows:
        at = row.get("@id")
        if not at:
            continue
        st, full = portal_json(f"{at}?format=json")
        out.append(full if isinstance(full, dict) and full.get("accession") else row)
    return out


def classify(files: "list[dict]") -> "dict[str, list[dict]]":
    reads, seqspecs, matrices, other = [], [], [], []
    for f in files:
        fmt = str(f.get("file_format") or "").lower()
        ctype = str(f.get("content_type") or "").lower()
        if fmt == "fastq" or ctype == "reads":
            reads.append(f)
        elif ctype == "seqspec" or (fmt in ("yaml", "json") and "seqspec" in ctype):
            seqspecs.append(f)
        elif fmt == "yaml":
            seqspecs.append(f)
        elif fmt in _MATRIX_FORMATS and any(h in ctype for h in _MATRIX_HINTS):
            matrices.append(f)
        else:
            other.append(f)
    return {"reads": reads, "seqspecs": seqspecs, "matrices": matrices, "other": other}


# Preference order for a published matrix. The single-cell pipeline loads
# .h5ad natively, so picking it over a 13 GB tar of kallisto output is both a
# smaller transfer and one less unpacking step. Taking matrices[0] and hoping
# left the choice to whatever order the Portal happened to return.
_MATRIX_RANK = {"h5ad": 0, "h5": 1, "mtx": 2, "tar": 3, "csv": 4, "tsv": 5}


def rank_matrices(matrices: "list[dict]") -> "list[dict]":
    def key(m):
        fmt = str(m.get("file_format") or "").lower()
        try:
            size = float(m.get("file_size") or 0)
        except (TypeError, ValueError):
            size = 0.0
        return (_MATRIX_RANK.get(fmt, 9), size)
    return sorted(matrices, key=key)


def derived_analysis_sets(accession: str) -> "list[dict]":
    status, data = portal_json(
        f"/search/?type=AnalysisSet"
        f"&input_file_sets.accession={urllib.parse.quote(accession)}"
        f"&limit=100&format=json")
    return (data or {}).get("@graph") or []


def find_guide_library(accession: str) -> dict:
    """Locate the sgRNA library and its protospacer table for a CRISPR screen.

    IGVF records the library on the **AnalysisSet**, not on the
    MeasurementSet -- `construct_library_sets` with an
    `integrated_content_files` entry of content_type "guide RNA sequences".
    Looking on the measurement set (the obvious place) finds nothing and
    invites guessing the library from its name, which is how you end up
    assigning guides against the wrong screen.

    So this walks: measurement set -> derived AnalysisSets -> construct
    library -> guide table. A dataset with no AnalysisSet yet -- one still
    "in progress" -- has no link to find, and that is reported as such
    rather than papered over with a best-guess match.
    """
    out = {"accession": accession, "analysis_sets": [], "libraries": [],
           "guide_files": [], "resolved": False, "why": ""}
    # The accession may itself be the AnalysisSet that carries the link, so
    # look at it directly before looking for sets derived from it.
    candidates = []
    st0, own = portal_json(
        f"/search/?type=FileSet&accession={urllib.parse.quote(accession)}&format=json")
    for row in (own or {}).get("@graph") or []:
        if row.get("@id"):
            candidates.append(row)
    derived = derived_analysis_sets(accession)
    out["analysis_sets"] = [d.get("accession") for d in derived]
    for d in candidates + derived:
        st, full = portal_json(f"{d.get('@id')}?format=json")
        for cls in (full or {}).get("construct_library_sets") or []:
            lib = cls.get("accession")
            if lib and lib not in out["libraries"]:
                out["libraries"].append(lib)
            for icf in cls.get("integrated_content_files") or []:
                if ("guide" in str(icf.get("content_type", "")).lower()
                        and icf.get("accession") not in
                        {g["accession"] for g in out["guide_files"]}):
                    out["guide_files"].append(
                        {"accession": icf.get("accession"),
                         "library": lib,
                         "id": icf.get("@id"),
                         "content_type": icf.get("content_type")})
    if out["guide_files"]:
        out["resolved"] = True
        out["why"] = (f"guide table {out['guide_files'][0]['accession']} via "
                      f"library {out['guide_files'][0]['library']}")
    elif not derived and not any(
            (portal_json(f"{c.get('@id')}?format=json")[1] or {}).get(
                "construct_library_sets") for c in candidates):
        out["why"] = (f"{accession} has no AnalysisSet yet, and IGVF carries "
                      f"the construct_library_sets link there rather than on "
                      f"the measurement set. Until the dataset is processed "
                      f"the library is not recorded anywhere machine-readable "
                      f"-- supply it with --guide-library or --guide-table.")
    else:
        out["why"] = ("AnalysisSets exist but none declares a construct "
                      "library with guide sequences.")
    return out


# ─── Guide library index ────────────────────────────────────────────────────

# Sequence columns that can identify a construct, and whether a match against
# one is specific enough to trust. `barcode` is deliberately absent: the
# LDLR pegRNA library's barcodes are 6 bp, so its 1,740 barcodes cover 42% of
# all possible 6-mers -- scanning reads for them matched 100% of reads and
# 100% of those matches were ambiguous. A barcode is readable only at a known
# offset in a known amplicon, which the metadata does not give us.
_KEY_COLUMNS = ("spacer", "protospacer", "rt_template_sequence", "peg_sequence")
_MIN_KEY_LEN = 10        # shorter than this matches by chance far too often


# Matching a read against a construct library, exactly, without scanning
# every (length, offset) pair. The LDLR RT templates come in 37 distinct
# lengths, so the naive triple loop cost ~8,000 dict lookups per read per
# strand -- about ten hours for the eight libraries of one screen. Indexing
# the constructs by a fixed-length prefix brings that to one lookup per
# offset, ~120 per read, and returns identical counts because a full
# startswith() still confirms every candidate.
_PREFIX_K_MAX = 12


def build_matcher(seq_to_guide: "dict[str, str]") -> dict:
    """Prefix index over construct sequences. Exact, not heuristic."""
    k = min(_PREFIX_K_MAX, min((len(s) for s in seq_to_guide), default=1))
    pref: "dict[str, list[tuple[str, str]]]" = defaultdict(list)
    for sq, gid in seq_to_guide.items():
        pref[sq[:k]].append((sq, gid))
    # Longest first inside each bucket. Two constructs matching at the same
    # offset necessarily share their first k bases, so they land in the same
    # bucket -- which makes this ordering globally longest-first per offset,
    # and stops a short RT template nested in a longer one from winning.
    for v in pref.values():
        v.sort(key=lambda t: -len(t[0]))
    return {"k": k, "pref": dict(pref)}


def match_read(seq: str, m: dict) -> "Optional[str]":
    """The longest construct sequence contained in the read, else None.

    Longest, not first-found: the constructs nest, and a longer match is
    stronger evidence than a shorter one that happens to sit nearer the
    start of the read. Returning at the first hit instead disagreed with a
    full length-ordered scan on 1 read in 20,000 -- rare, but a silent and
    needless difference between the fast path and the obvious one.
    """
    k, pref = m["k"], m["pref"]
    best_gid, best_len = None, 0
    for off in range(0, len(seq) - k + 1):
        for sq, gid in pref.get(seq[off:off + k], ()):
            if len(sq) <= best_len:
                break          # bucket is longest-first: the rest are shorter
            if seq.startswith(sq, off):
                best_gid, best_len = gid, len(sq)
                break
    return best_gid


def load_guide_index(file_id: str) -> dict:
    """Build a sequence -> guide_id index, choosing the discriminating column.

    A CRISPR-KO library has one unique spacer per guide, so the spacer is the
    counting key. A prime-editing library does NOT: the LDLR pegRNA library
    has 1,740 constructs sharing just 52 spacers, because the spacer only
    sets the nick site and the variant lives in the RT template. Counting
    that library by spacer collapses 1,740 pegRNAs onto 52 keys and then
    attributes every read to whichever construct happened to be first --
    confident, precise, meaningless numbers. So pick the key by measuring
    uniqueness rather than assuming a column.
    """
    st, obj = portal_json(
        f"{file_id}?format=json" if file_id.startswith("/")
        else f"/tabular-files/{file_id}/?format=json")
    href = (obj or {}).get("href")
    out = {"resolved": False, "why": "", "key": "", "seq_to_guide": {},
           "lengths": [], "target": {}, "type": {}, "n_rows": 0,
           "n_constructs": 0, "candidates": []}
    if not href:
        out["why"] = f"no href for {file_id}"
        return out
    dest = REF_DIR / "guides" / Path(href).name
    if not dest.exists():
        portal_download(href, dest)
    text = read_text_maybe_gzip(dest)
    lines = text.splitlines()
    if not lines:
        out["why"] = f"{dest.name} is empty"
        return out
    delim = max(("\t", ",", ";"), key=lambda d: len(lines[0].split(d)))
    # csv.reader, not str.split: these tables carry quoted fields containing
    # commas (putative_target_genes is '["ENSG00000130164"]'), and splitting
    # naively shifts every column after it.
    rows = list(csv.reader(lines, delimiter=delim))
    hdr = [h.strip().lower() for h in rows[0]]
    idx = {h: i for i, h in enumerate(hdr)}

    def cell(r: "list[str]", *names: str) -> str:
        for n in names:
            i = idx.get(n)
            if i is not None and len(r) > i:
                return r[i].strip()
        return ""

    recs = []
    for r in rows[1:]:
        gid = cell(r, "guide_id", "guide", "name", "id")
        if not gid:
            continue
        recs.append((gid, r))
    out["n_rows"] = len(recs)
    if not recs:
        out["why"] = f"{dest.name} has no guide_id column or no rows"
        return out

    # Score every candidate column by how many constructs it separates.
    for col in _KEY_COLUMNS:
        if col not in idx:
            continue
        seqs = {}
        for gid, r in recs:
            v = cell(r, col).upper()
            if v and len(v) >= _MIN_KEY_LEN and set(v) <= set("ACGTN"):
                seqs.setdefault(v, gid)
        if seqs:
            ls = sorted(len(k) for k in seqs)
            out["candidates"].append(
                {"column": col, "distinct": len(seqs),
                 "fraction": round(len(seqs) / len(recs), 4),
                 "median_len": ls[len(ls) // 2],
                 "seq_to_guide": seqs,
                 # Longest first: RT templates nest, and the first match
                 # found would otherwise hand the read to the least
                 # specific construct.
                 "lengths": sorted(set(ls), reverse=True)})
    if not out["candidates"]:
        out["why"] = (f"{dest.name} has no usable sequence column "
                      f"(looked for {', '.join(_KEY_COLUMNS)})")
        return out
    for gid, r in recs:
        out["target"][gid] = cell(r, "intended_target_name", "target",
                                  "genomic_element") or gid
        out["type"][gid] = cell(r, "type", "guide_type", "targeting")
    # Discrimination first, then shorter sequences: a key longer than the
    # read can never be found in it. peg_sequence separates this library as
    # well as rt_template_sequence does, but its entries are 128-134 bp
    # against 128 bp reads, so choosing it by uniqueness alone would have
    # produced a 0% assignment rate. calibrate_key() then confirms against
    # real reads rather than trusting this ordering.
    out["candidates"].sort(key=lambda c: (-c["distinct"], c["median_len"]))
    out["resolved"] = True
    return out


def calibrate_key(idx: dict, accession: str, sample: int = 20000) -> dict:
    """Choose the counting key by testing candidates against real reads.

    Uniqueness in the metadata says a key *can* tell constructs apart; only
    the reads say whether it is *present* to be found. For the LDLR pegRNA
    library the measured rates are spacer 7% (22% of hits ambiguous),
    peg_sequence 0% (key longer than the read), rt_template_sequence 62%
    (0.8% ambiguous) -- a ranking no amount of metadata would have revealed.
    """
    comp = str.maketrans("ACGTN", "TGCAN")
    reads: "list[str]" = []
    for f in list_files(accession):
        if str(f.get("file_format", "")).lower() != "fastq":
            continue
        name = Path(str(f.get("href") or f.get("accession"))).name
        dest = FASTQ_CACHE / name
        if not dest.exists():
            FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
            portal_download(f["href"], dest)
        op = gzip.open if dest.read_bytes()[:2] == b"\x1f\x8b" else open
        with op(dest, "rt", errors="replace") as fh:       # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    reads.append(line.strip().upper())
                    if len(reads) >= sample:
                        break
        if reads:
            break
    if not reads:
        return {"chosen": idx["candidates"][0] if idx["candidates"] else None,
                "read_len": 0, "tested": [],
                "why": f"no FASTQ reads available from {accession} to calibrate"}

    read_len = Counter(len(r) for r in reads).most_common(1)[0][0]
    tested = []
    for c in idx["candidates"]:
        mt = build_matcher(c["seq_to_guide"])
        hit = amb = 0
        for r0 in reads:
            found = set()
            for s_ in (r0, r0.translate(comp)[::-1]):
                g = match_read(s_, mt)
                if g:
                    found.add(g)
            if found:
                hit += 1
                if len(found) > 1:
                    amb += 1
        tested.append({**{k: v for k, v in c.items()
                          if k not in ("seq_to_guide", "lengths")},
                       "rate": hit / len(reads),
                       "ambiguous": amb / max(hit, 1)})
    # Best assignment rate wins; discrimination breaks near-ties. A key that
    # is present but cannot separate constructs is useless, and so is one
    # that separates them but never appears.
    ranked = sorted(zip(tested, idx["candidates"]),
                    key=lambda t: (-round(t[0]["rate"], 2), -t[0]["distinct"]))
    return {"chosen": ranked[0][1] if ranked and ranked[0][0]["rate"] > 0 else None,
            "read_len": read_len, "tested": tested, "n_reads": len(reads),
            "why": ""}


def load_guide_table(file_id: str) -> "list[tuple[str, str]]":
    """Download a guide table and return [(guide_id, spacer), ...]."""
    st, obj = portal_json(f"{file_id}?format=json" if file_id.startswith("/")
                          else f"/tabular-files/{file_id}/?format=json")
    href = (obj or {}).get("href")
    if not href:
        return []
    dest = REF_DIR / "guides" / Path(href).name
    if not dest.exists():
        portal_download(href, dest)
    text = read_text_maybe_gzip(dest)
    rows = text.splitlines()
    if not rows:
        return []
    # Sniff the delimiter instead of assuming a tab. IGVF publishes guide
    # tables as BOTH .tsv and .csv -- PALB2's editing templates are tab
    # separated, IGVFFI4591THXG is comma separated -- and splitting a CSV on
    # tabs yields one giant field per line, no "spacer" column, and a silent
    # "0 guides" for a file holding 8,192 of them. Pick whichever delimiter
    # actually splits the header.
    delim = max(("\t", ",", ";"), key=lambda d: len(rows[0].split(d)))
    if len(rows[0].split(delim)) < 2:
        return []
    hdr = [h.strip().lower() for h in rows[0].split(delim)]

    def _col(*names: str) -> "Optional[int]":
        for n in names:
            if n in hdr:
                return hdr.index(n)
        return None

    gi = _col("guide_id", "guide", "name", "id")
    si = _col("spacer", "protospacer", "sequence", "guide_sequence", "sgrna")
    if si is None:
        return []
    if gi is None:
        gi = 0 if si != 0 else 1
    out = []
    for r in rows[1:]:
        parts = r.split(delim)
        if len(parts) > max(gi, si):
            gid, sp = parts[gi].strip(), parts[si].strip().upper()
            if sp and set(sp) <= set("ACGTN"):
                out.append((gid, sp))
    return out


def gb(n: Any) -> float:
    try:
        return round(float(n) / 1e9, 3)
    except (TypeError, ValueError):
        return 0.0


# ─── Read pairing ───────────────────────────────────────────────────────────

def normalize_read_id(rid: str) -> str:
    """Fold the many spellings of a read label onto R1 / R2 / I1.

    seqspec documents in the wild use `R1`, `Read1`, `read_1` and
    `index1` interchangeably. Matching the literal string meant a
    perfectly good spec parsed to `{'READ2': 110}`, which then failed
    chemistry detection for a dataset whose structure was fully described.
    """
    t = re.sub(r"[^a-z0-9]", "", (rid or "").lower())
    if t in ("r1", "read1"):
        return "R1"
    if t in ("r2", "read2"):
        return "R2"
    if t in ("i1", "index1", "index"):
        return "I1"
    if t in ("i2", "index2"):
        return "I2"
    return (rid or "").upper()


def seqspec_reads(text: str) -> "list[dict]":
    """Every !Read block as {id, role, length, files}.

    ``read_id`` is not a reliable label. Across real IGVF seqspecs it is
    ``R1``/``R2``, ``Read1``/``Read2``, or -- on the multiome sets -- the
    *file accession itself*. So the role is derived from the read's length,
    which is a property of the chemistry rather than of the submitter's
    naming: a 26 or 28 bp read is the cell barcode + UMI, a short read is a
    sample index, and the long read is the cDNA insert.
    """
    out = []
    for block in re.split(r"\n- !Read\n", text)[1:]:
        rid = re.search(r"read_id:\s*(\S+)", block)
        mx = re.search(r"max_len:\s*(\d+)", block)
        ids = re.findall(r"file_id:\s*(\S+)", block)
        if not rid:
            continue
        length = int(mx.group(1)) if mx else 0
        name = normalize_read_id(rid.group(1))
        if name in ("R1", "R2", "I1", "I2"):
            role = {"R1": "barcode", "R2": "cdna",
                    "I1": "index", "I2": "index"}[name]
        elif length in _TECH_BY_R1_LEN:
            role = "barcode"
        elif length and length <= 12:
            role = "index"
        elif length >= 40:
            role = "cdna"
        else:
            role = "unknown"
        out.append({"id": rid.group(1), "role": role,
                     "length": length, "files": ids})
    return out


def seqspec_read_files(text: str) -> "dict[str, list[str]]":
    """Map each read to the file accessions it contains, per the seqspec.

    This is the authoritative pairing source. `illumina_read_type` is
    absent on some Files (it is None for every read on IGVFDS6639ECQN),
    so metadata-only pairing declares a fully specified dataset unpairable.
    The seqspec always states which file belongs to which read.
    """
    out: "dict[str, list[str]]" = {}
    for read in seqspec_reads(text):
        if read["role"] == "barcode":
            out.setdefault("R1", []).extend(read["files"])
        elif read["role"] == "cdna":
            out.setdefault("R2", []).extend(read["files"])
    return out


def pair_reads_from_seqspec(reads: "list[dict]", mapping: "dict[str, list[str]]"
                             ) -> "tuple[list[tuple[dict, dict]], list[dict]]":
    """Pair using the seqspec's per-read file lists, in submitted order."""
    by_acc = {str(f.get("accession")): f for f in reads}
    r1 = [by_acc[a] for a in mapping.get("R1", []) if a in by_acc]
    r2 = [by_acc[a] for a in mapping.get("R2", []) if a in by_acc]
    pairs = list(zip(r1, r2))
    used = {id(f) for pair in pairs for f in pair}
    unpaired = [f for f in reads if id(f) not in used]
    return pairs, unpaired


def pair_reads(reads: "list[dict]") -> "tuple[list[tuple[dict, dict]], list[dict]]":
    """Group FASTQs into (R1, R2) pairs on flowcell + lane + sequencing run.

    kb count consumes reads strictly positionally -- barcode file then
    biological file, repeated per pair -- so a mispairing does not error,
    it silently quantifies the wrong barcodes against the wrong cDNA. The
    grouping key is therefore the sequencing identity of the lane, and
    anything that fails to pair is returned separately rather than being
    appended and hoped for.
    """
    buckets: "dict[tuple, dict[str, dict]]" = {}
    for f in reads:
        key = (f.get("flowcell_id"), f.get("lane"), f.get("sequencing_run"))
        rt = str(f.get("illumina_read_type") or "").upper()
        if rt not in ("R1", "R2"):
            continue
        buckets.setdefault(key, {})[rt] = f
    pairs, unpaired = [], []
    for key in sorted(buckets, key=lambda k: tuple(str(x) for x in k)):
        slot = buckets[key]
        if "R1" in slot and "R2" in slot:
            pairs.append((slot["R1"], slot["R2"]))
        else:
            unpaired.extend(slot.values())
    known = {id(f) for pair in pairs for f in pair} | {id(f) for f in unpaired}
    unpaired.extend(f for f in reads if id(f) not in known)
    return pairs, unpaired


def technology_from_seqspec(text: str) -> "tuple[Optional[str], str]":
    """Infer the kb technology string from a seqspec document.

    Returns (technology, evidence). The barcode+UMI length on read 1 is the
    discriminator: 26 bases is 10x v2, 28 is v3.
    """
    reads = seqspec_reads(text)
    barcode = next((r for r in reads if r["role"] == "barcode"), None)
    if barcode and barcode["length"] in _TECH_BY_R1_LEN:
        tech = _TECH_BY_R1_LEN[barcode["length"]]
        return tech, (f"barcode read {barcode['id']} is {barcode['length']} bp "
                      f"-> {tech}")
    lib = re.search(r"lib_struct:\s*(\S+)", text)
    if lib and "10xChromium3" in lib.group(1):
        return "10XV3", f"lib_struct names 10x Chromium 3' ({lib.group(1)})"
    # No barcode-bearing read at all: the library is bulk, and kallisto's
    # BULK technology quantifies it directly. Saying "unknown chemistry"
    # for a single long cDNA read sends the user hunting for a 10x version
    # that was never involved.
    cdna = [r for r in reads if r["role"] == "cdna"]
    if cdna and not barcode:
        return "BULK", (f"one biological read ({cdna[0]['id']}, "
                        f"{cdna[0]['length']} bp) and no barcode read "
                        f"-> bulk, not single-cell")
    shape = {r["id"]: f"{r['length']}bp/{r['role']}" for r in reads}
    return None, (f"seqspec read structure {shape or '{}'} matches no known "
                  f"chemistry; pass --technology explicitly")


# ─── kb-python ──────────────────────────────────────────────────────────────

def kb_available() -> "tuple[bool, str]":
    exe = shutil.which("kb")
    if not exe:
        return False, ("kb-python is not installed. Install the aligner extra: "
                        "pip install 'igvfagent[align]'")
    try:
        out = subprocess.run([exe, "--version"], capture_output=True, text=True,
                              timeout=60)
        ver = (out.stdout + out.stderr).strip().splitlines()
        ver = next((l for l in ver if "kb_python" in l), ver[-1] if ver else "?")
    except Exception as e:                                   # noqa: BLE001
        return False, f"kb found at {exe} but did not run: {e}"
    return True, f"{exe} ({ver})"


def ensure_index(reference: str, dry_run: bool) -> "tuple[Path, Path, list[str]]":
    """Return (index, t2g) for a prebuilt kb reference, downloading if absent.

    Two things here are not optional in a containerised deployment.

    ``--tmp`` puts kb's staging directory on the *same filesystem* as the
    destination. By default kb downloads to ``./tmp`` relative to the working
    directory and then ``os.rename()``s the result into place; when the
    destination is a mounted data volume and the working directory is not,
    that rename is a cross-device link and dies with EXDEV:

        OSError: [Errno 18] Invalid cross-device link:
            'tmp/index.idx' -> '/workspace/Data/References/kb/human/index.idx'

    And the result is *verified* rather than trusted, because kb exits **0**
    on exactly that failure (measured, not assumed). Trusting the exit code
    let a run continue to `kb count` against an index that was never written,
    where it failed far from the real cause with "kallisto index file not
    found".
    """
    ref_root = REF_DIR / safe_label(reference)
    index = ref_root / "index.idx"
    t2g = ref_root / "t2g.txt"
    tmp = ref_root / "tmp"
    cmd = ["kb", "ref", "-d", reference, "-i", str(index), "-g", str(t2g),
           "--tmp", str(tmp)]
    if index.exists() and t2g.exists():
        return index, t2g, []
    ref_root.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return index, t2g, cmd
    logging.info("Building/downloading kb reference %s -> %s", reference, ref_root)
    # kb refuses to start if its tmp directory already exists.
    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)
    subprocess.run(cmd, check=True)
    missing = [str(p) for p in (index, t2g) if not p.exists()]
    if missing:
        raise RuntimeError(
            "kb ref reported success but did not write: " + ", ".join(missing)
            + f". Check write permission and free space on {ref_root}.")
    return index, t2g, cmd


def kb_count_cmd(index: Path, t2g: Path, tech: str, out_dir: Path,
                  pairs: "list[tuple[Path, Path]]", threads: int,
                  workflow: str = "standard",
                  singles: "Optional[list[Path]]" = None) -> "list[str]":
    """Assemble a `kb count` invocation.

    Barcoded chemistries take reads in positional pairs; BULK takes the
    biological reads on their own. Feeding a single-end bulk library through
    the paired path would hand kallisto one file where it expects two and
    quantify against a barcode read that does not exist.
    """
    cmd = ["kb", "count", "-i", str(index), "-g", str(t2g), "-x", tech,
           "-o", str(out_dir), "--h5ad", "-t", str(threads)]
    if workflow != "standard":
        cmd += ["--workflow", workflow]
    if singles:
        # kb refuses BULK without --parity: it cannot know whether two files
        # are a pair or two single-end runs. These are single-end by
        # construction (one biological read, no barcode read), so say so.
        cmd += ["--parity", "single"]
        cmd += [str(p) for p in singles]
    else:
        for r1, r2 in pairs:
            cmd += [str(r1), str(r2)]
    return cmd


# ─── Routing ────────────────────────────────────────────────────────────────

def assay_mismatch(fs: dict) -> "Optional[tuple[str, str, str]]":
    """(label, right_analysis, matched_on) when reads are not transcripts.

    Delegates to Scripts/_assays.py, which classifies every one of the 65
    assay titles on the Portal rather than listing the few that had already
    caused a wrong answer. That inversion is the point: only 59.5% of the
    11,070 MeasurementSets are transcript assays, so a default of "quantify
    it as RNA" was wrong for 4,486 of them, silently. An unrecognised assay
    is now reported as unrecognised instead of being quantified.
    """
    cls = _A.classify(fs or {})
    route = cls["route"]
    if route == _A.TRANSCRIPT:
        return None
    if route == _A.UNKNOWN:
        return (f"an unrecognised assay ({cls['assay']})",
                "unknown -- IGVFagent has no rule for this assay, so it "
                "will not guess. Transcriptome quantification is the wrong "
                "answer for 40% of IGVF assays and would look plausible "
                "here",
                cls["matched_on"])
    return (f"{cls['assay']} (route: {route})",
            f"{cls['analysis']} -- {cls['support']}",
            cls["matched_on"])
    # `crispr_screen_readout` is the decisive field for a CRISPR screen and
    # it is not an assay title, so a title-only check missed it entirely.
    # "gRNA sequencing" means the reads ARE the guide library: quantifying
    # them against a transcriptome gave 1.4% pseudoalignment on
    # IGVFDS6464SOVZ, which is the measurement telling us it was the wrong
    # question. A screen whose readout is scRNA-seq (TAP-seq, Perturb-seq)
    # is the opposite case and must still route to alignment.
    readout = str(fs.get("crispr_screen_readout") or "").strip().lower()
    if readout in ("grna sequencing", "sgrna sequencing", "guide sequencing"):
        return ("a CRISPR screen with gRNA-sequencing readout",
                "per-guide counts and enrichment between sorted "
                "populations -- run `igvfagent raw-pipeline guide-count "
                "<accession>` (tool: raw_pipeline_guide_count)",
                f"crispr_screen_readout={fs.get('crispr_screen_readout')}")
    fields = []
    for key in ("preferred_assay_titles", "assay_titles",
                 "preferred_assay_slims", "assay_slims"):
        v = fs.get(key)
        if isinstance(v, list):
            fields.extend((key, str(x)) for x in v)
        elif v:
            fields.append((key, str(v)))
    for key, raw in fields:
        hit = _NON_TRANSCRIPT_ASSAYS.get(raw.strip().lower())
        if hit:
            return hit[0], hit[1], f"{key}={raw}"
    return None


def decide_route(accession: str, buckets: "dict[str, list[dict]]",
                  force_align: bool, file_set: "Optional[dict]" = None) -> "dict":
    if buckets["matrices"] and not force_align:
        return {"route": "matrix_on_set",
                "matrices": rank_matrices(buckets["matrices"]),
                "why": "the set already carries a count matrix"}
    derived = derived_analysis_sets(accession)
    if derived and not force_align:
        for ds in derived:
            dacc = ds.get("accession")
            if not dacc:
                continue
            dfiles = classify(list_files(dacc))
            if dfiles["matrices"]:
                return {"route": "matrix_derived",
                        "matrices": rank_matrices(dfiles["matrices"]),
                        "analysis_set": dacc,
                        "why": f"AnalysisSet {dacc} derived from {accession} "
                               f"already publishes a count matrix"}
    if buckets["reads"]:
        mism = assay_mismatch(file_set or {})
        if mism and not force_align:
            label, right, matched = mism
            return {"route": "assay_mismatch", "assay": label,
                    "matched_on": matched, "right_analysis": right,
                    "why": (f"this is {label}, whose reads are not a "
                            f"transcript library. Quantifying them against a "
                            f"transcriptome would report which gene the "
                            f"amplicon covers -- the assay design restated, "
                            f"not a result. The readout you want is {right}. "
                            f"Pass force_align=true to quantify anyway.")}
        return {"route": "align", "why": "only raw reads are available; "
                                          "they must be quantified first",
                "derived_sets": [d.get("accession") for d in derived]}
    return {"route": "none", "why": "no reads and no matrix were found"}


def write_json(path: Path, obj: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str))
    return path


# ─── Commands ───────────────────────────────────────────────────────────────

def _inventory(accession: str, force_align: bool) -> "Optional[dict]":
    fs = resolve_file_set(accession)
    if fs is None:
        print(f"RESOLVED: no — {accession} was not returned by the Portal.")
        print("  If it exists but is unreleased, supply credentials that reach "
              "it (IGVF_ACCESS_KEY + IGVF_SECRET_ACCESS_KEY, or "
              "Docs/Secret/IGVFportalAPI.txt) and re-run; `igvfagent "
              "auth-check` shows what is in effect.")
        return None
    files = list_files(accession)
    buckets = classify(files)

    # Pair from the seqspec when there is one -- it states which file
    # belongs to which read, and Files may carry no illumina_read_type at
    # all. Fall back to sequencing metadata, then keep whichever pairing
    # actually produced pairs.
    seqspec_text = ""
    pairs, unpaired = pair_reads(buckets["reads"])
    pairing_source = "file metadata (flowcell + lane + run)"
    if buckets["reads"] and buckets["seqspecs"]:
        # Merge every seqspec, not just the first. A multiome set carries one
        # per modality (RNA, ATAC, the MULTI-seq tag library), so reading only
        # seqspecs[0] describes a third of the reads and leaves the rest
        # looking unpairable.
        merged: "dict[str, list[str]]" = {}
        used_specs = []
        for sp in buckets["seqspecs"]:
            try:
                tmp = RUN_DIR / "_seqspec" / f"{sp.get('accession')}.yaml"
                portal_download(sp.get("href"), tmp)
                text = read_text_maybe_gzip(tmp)
                if not seqspec_text:
                    seqspec_text = text
                for k, v in seqspec_read_files(text).items():
                    merged.setdefault(k, []).extend(v)
                used_specs.append(str(sp.get("accession")))
            except Exception as e:                           # noqa: BLE001
                logging.warning("seqspec %s unavailable: %s", sp.get("accession"), e)
        if merged:
            sp_pairs, sp_unpaired = pair_reads_from_seqspec(buckets["reads"], merged)
            if len(sp_pairs) >= len(pairs):
                pairs, unpaired = sp_pairs, sp_unpaired
                pairing_source = "seqspec " + ", ".join(used_specs)

    route = decide_route(accession, buckets, force_align, file_set=fs)
    total_gb = round(sum(gb(f.get("file_size")) for f in buckets["reads"]), 2)
    print(f"Accession:   {accession}")
    print(f"Type:        {(fs.get('@type') or ['?'])[0]}  |  status: {fs.get('status')}")
    print(f"Summary:     {fs.get('summary') or '-'}")
    print(f"Assay:       {', '.join(fs.get('preferred_assay_titles') or []) or '-'}")
    print(f"Files:       {len(files)}  "
          f"(reads {len(buckets['reads'])}, seqspec {len(buckets['seqspecs'])}, "
          f"matrix {len(buckets['matrices'])}, other {len(buckets['other'])})")
    print(f"Read pairs:  {len(pairs)}  [{pairing_source}]"
          + (f"  |  UNPAIRED {len(unpaired)}" if unpaired else ""))
    print(f"Read bytes:  {total_gb} GB")
    ctrl = [f.get("accession") for f in files if f.get("controlled_access")]
    if ctrl:
        print(f"Controlled:  {len(ctrl)} file(s) need an approved DUA: "
              f"{', '.join(str(c) for c in ctrl[:5])}")
    print(f"ROUTE:       {route['route']} — {route['why']}")
    return {"file_set": fs, "files": files, "buckets": buckets,
            "pairs": pairs, "unpaired": unpaired, "route": route,
            "read_gb": total_gb, "seqspec_text": seqspec_text}


def cmd_plan(args: argparse.Namespace) -> int:
    inv = _inventory(args.accession, args.force_align)
    if inv is None:
        return 2
    ok, detail = kb_available()
    print(f"Aligner:     {'available' if ok else 'MISSING'} — {detail}")
    if inv["route"]["route"] == "align":
        tech, evidence = args.technology, "given with --technology"
        if not tech and inv.get("seqspec_text"):
            tech, evidence = technology_from_seqspec(inv["seqspec_text"])
        elif not tech:
            evidence = "no seqspec on this set; pass --technology"
        print(f"Technology:  {tech or 'UNKNOWN'} — {evidence}")
        print(f"Would download {inv['read_gb']} GB of reads, build/reuse the "
              f"'{args.reference}' index, then run kb count and the "
              f"single-cell pipeline.")
        if args.max_download_gb and inv["read_gb"] > args.max_download_gb:
            print(f"NOTE: that exceeds --max-download-gb {args.max_download_gb}; "
                  f"`run` would stop before downloading.")
    out = write_json(
        REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.accession)}_plan.json",
        {"accession": args.accession, "route": inv["route"],
         "n_files": len(inv["files"]), "read_gb": inv["read_gb"],
         "n_pairs": len(inv["pairs"]), "n_unpaired": len(inv["unpaired"]),
         "aligner_available": ok, "aligner": detail})
    print(f"Plan:        {out}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    label = args.label or f"{stamp}_{safe_label(args.accession)}"

    if getattr(args, "detach", False):
        # Re-invoke this same command without --detach, in its own session.
        argv = [_igvfagent(), "raw-pipeline", "run", args.accession,
                "--label", label,
                "--reference", args.reference,
                "--max-download-gb", str(args.max_download_gb),
                "--workflow", args.workflow,
                "--threads", str(args.threads)]
        if args.technology:
            argv += ["--technology", args.technology]
        if args.force_align:
            argv.append("--force-align")
        if args.skip_analysis:
            argv.append("--skip-analysis")
        if getattr(args, "no_reuse", False):
            argv.append("--no-reuse")
        if args.index and args.t2g:
            argv += ["--index", args.index, "--t2g", args.t2g]
        rec = spawn_detached(argv, label, args.accession)
        print(f"STARTED detached job: {label}")
        print(f"  accession: {args.accession}")
        print(f"  pid:       {rec['pid']}")
        print(f"  log:       {rec['log']}")
        print(f"  work dir:  {rec['work']}")
        print("\nThe job continues after this call returns and after this "
              "conversation ends.")
        print(f"Check on it with:  igvfagent raw-pipeline status {label}")
        return 0

    # A completed identical analysis is worth minutes, not a re-run. This is
    # the demo case: the same accession is entered deliberately, and
    # re-downloading 45 GB to recompute a matrix that already exists on this
    # machine wastes the audience's time as well as the disk.
    rkey = result_key(args.accession, args.technology or "auto",
                      args.workflow, args.reference)
    prior = None if getattr(args, "no_reuse", False) else lookup_result(rkey)
    if prior and not args.force_align and not args.dry_run:
        print(f"REUSING a completed run of this exact analysis "
              f"({prior.get('completed')}).")
        print(f"  matrix:   {prior['matrix']}")
        print(f"  from run: {prior['work']}")
        for k in ("n_processed", "p_pseudoaligned", "shape"):
            if prior.get(k) is not None:
                print(f"  {k}: {prior[k]}")
        print("  Pass no_reuse=true (--no-reuse) to recompute from the reads.")
        matrix = Path(prior["matrix"])
        if not args.skip_analysis:
            sc = [_igvfagent(), "sc-analyze", "pipeline", "--input",
                  str(matrix), "--label", label]
            print("Analysis: " + " ".join(sc))
            subprocess.run(sc, check=False)
        print(f"Run dir: {prior['work']}")
        return 0

    inv = _inventory(args.accession, args.force_align)
    if inv is None:
        return 2
    route = inv["route"]["route"]
    work = RUN_DIR / label
    work.mkdir(parents=True, exist_ok=True)

    matrix: "Optional[Path]" = None
    # Bulk and single-cell need different downstream treatment, and the
    # difference is not cosmetic -- see the guard before the analysis step.
    is_bulk = False

    if route in ("matrix_on_set", "matrix_derived"):
        target = inv["route"]["matrices"][0]
        others = inv["route"]["matrices"][1:]
        print(f"Downloading published matrix {target.get('accession')} "
              f"({target.get('file_format')}, {gb(target.get('file_size'))} GB)"
              + (f" — chosen over {len(others)} other matrix file(s): "
                 + ", ".join(f"{m.get('accession')} ({m.get('file_format')}, "
                              f"{gb(m.get('file_size'))} GB)" for m in others)
                 if others else "") + "…")
        if args.max_download_gb and gb(target.get("file_size")) > args.max_download_gb:
            print(f"STOPPING: that exceeds --max-download-gb "
                  f"{args.max_download_gb}. Nothing was downloaded.")
            return 4
        if args.dry_run:
            print(f"DRY RUN: would download {target.get('href')}")
            matrix = work / (Path(str(target.get("href") or "matrix")).name
                             or "matrix")
        else:
            name = Path(str(target.get("href") or "matrix")).name or "matrix"
            matrix = portal_download(target["href"], work / name)
            print(f"Matrix: {matrix}")

    elif route == "align":
        ok, detail = kb_available()
        if not ok:
            print(f"CANNOT ALIGN: {detail}")
            print("  Nothing was downloaded. Install the extra and re-run, or "
                  "point --index at an existing kallisto index.")
            return 3
        tech = args.technology
        evidence = "given with --technology"
        if not tech and inv.get("seqspec_text"):
            (work / "seqspec.yaml").write_text(inv["seqspec_text"])
            tech, evidence = technology_from_seqspec(inv["seqspec_text"])
        if not tech:
            print(f"CANNOT ALIGN: chemistry unknown — {evidence}")
            return 3
        print(f"Technology: {tech} — {evidence}")

        pairs = inv["pairs"]
        single_end = tech.upper() == "BULK"
        is_bulk = single_end
        if not pairs and not single_end:
            print("CANNOT ALIGN: no R1/R2 pairs could be formed from the "
                  "reads on this set (see UNPAIRED above). Pairing comes "
                  "from the seqspec when present, else flowcell_id + lane + "
                  "sequencing_run.")
            return 3

        if args.max_download_gb and inv["read_gb"] > args.max_download_gb:
            print(f"STOPPING: reads total {inv['read_gb']} GB, over "
                  f"--max-download-gb {args.max_download_gb}. Nothing was "
                  f"downloaded. Raise the cap to proceed.")
            return 4

        if args.index and args.t2g:
            index, t2g, refcmd = Path(args.index), Path(args.t2g), []
        else:
            index, t2g, refcmd = ensure_index(args.reference, args.dry_run)
        if refcmd:
            print("Reference: " + " ".join(refcmd))

        local_pairs = []
        local_singles = []
        fq_dir = work / "fastq"
        fq_dir.mkdir(parents=True, exist_ok=True)
        FASTQ_CACHE.mkdir(parents=True, exist_ok=True)

        def _fetch(f: dict) -> Path:
            """Return a local path for one read file, downloading if needed.

            The bytes live once in FASTQ_CACHE, keyed by file accession, and
            each run links to them. Re-running a dataset therefore costs no
            transfer and no second copy on disk -- which matters for a demo,
            where the same accession is entered deliberately.
            """
            name = Path(str(f.get("href") or f.get("accession"))).name
            cached = FASTQ_CACHE / name
            link = fq_dir / name
            size = float(f.get("file_size") or 0)
            if not cached.exists():
                print(f"Downloading {f.get('accession')} ({gb(size)} GB)…")
                portal_download(f["href"], cached, on_bytes=_tick)
            else:
                print(f"Cached      {f.get('accession')} ({gb(size)} GB)")
                _tick(size)
            if not link.exists():
                try:
                    os.link(cached, link)          # same filesystem: free
                except OSError:
                    link.symlink_to(cached)
            return link
        # Byte-accurate progress across the whole transfer: the download is
        # the long pole (45.6 GB here), and it is the one phase whose total
        # is known up front.
        _all_reads = ([f for pair in pairs for f in pair] if pairs
                      else list(inv["buckets"]["reads"]))
        _total_bytes = sum(float(f.get("file_size") or 0) for f in _all_reads) or 1.0
        _done = {"n": 0.0, "last": 0.0}

        def _tick(n: int) -> None:
            _done["n"] += n
            now = time.time()
            if now - _done["last"] < 2.0:      # cap writes at ~1 every 2s
                return
            _done["last"] = now
            write_progress(work, phase="download",
                           bytes_done=int(_done["n"]),
                           bytes_total=int(_total_bytes),
                           percent=round(100.0 * _done["n"] / _total_bytes, 1),
                           detail=f"{_done['n']/1e9:.1f} / {_total_bytes/1e9:.1f} GB")

        write_progress(work, phase="download", bytes_done=0,
                       bytes_total=int(_total_bytes), percent=0.0,
                       detail=f"0 / {_total_bytes/1e9:.1f} GB")
        if single_end:
            for f in inv["buckets"]["reads"]:
                name = Path(str(f.get("href") or f.get("accession"))).name
                if args.dry_run:
                    print(f"DRY RUN: would fetch {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB)")
                    local_singles.append(fq_dir / name)
                else:
                    local_singles.append(_fetch(f))
        for r1, r2 in pairs:
            pq = []
            for f in (r1, r2):
                name = Path(str(f.get("href") or f.get("accession"))).name
                if args.dry_run:
                    print(f"DRY RUN: would fetch {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB)")
                    pq.append(fq_dir / name)
                else:
                    pq.append(_fetch(f))
            local_pairs.append((pq[0], pq[1]))

        kb_out = work / "kb"
        cmd = kb_count_cmd(index, t2g, tech, kb_out, local_pairs,
                            args.threads, args.workflow,
                            singles=local_singles or None)
        write_progress(work, phase="align", percent=None,
                       detail=f"kb count -x {tech} on {len(local_pairs) or len(local_singles)} "
                              f"input(s); no per-read progress is available")
        print("Aligner: " + " ".join(cmd))
        if args.dry_run:
            print("DRY RUN: would run kb count as printed above.")
            matrix = kb_out / "counts_unfiltered" / "adata.h5ad"
        else:
            subprocess.run(cmd, check=True)
            cand = sorted(kb_out.rglob("*.h5ad"))
            if not cand:
                print("kb count produced no .h5ad — inspect " + str(kb_out))
                return 5
            matrix = cand[0]
            print(f"Matrix: {matrix}")
    elif route == "assay_mismatch":
        print(f"\nNOT RUNNING: {inv['route']['why']}")
        print(f"  assay:      {inv['route'].get('assay')}")
        print(f"  matched on: {inv['route'].get('matched_on')}")
        print(f"  right analysis: {inv['route'].get('right_analysis')}")
        print("  Nothing was downloaded and no matrix was written, because a "
              "gene-count matrix for this assay would answer a different "
              "question than the one asked.")
        return 6
    else:
        print("Nothing to process: no reads and no matrix.")
        return 2

    if matrix:
        write_progress(work, phase="analyse", percent=None,
                       detail=f"matrix ready: {matrix.name}")
        extra = {}
        ri = matrix.parent.parent / "run_info.json"
        if ri.exists():
            try:
                info = json.loads(ri.read_text())
                extra = {"n_processed": info.get("n_processed"),
                         "p_pseudoaligned": info.get("p_pseudoaligned")}
            except (OSError, ValueError):
                pass
        record_result(rkey, matrix, work, extra)

    if matrix and not args.skip_analysis and is_bulk:
        # Never hand a bulk matrix to the single-cell pipeline. It is one
        # sample, so "filter genes seen in fewer than N cells" removes every
        # gene, leaving a (1, 0) matrix that log1p rejects with "Found array
        # with 0 sample(s)" -- an error that reads like the quantification
        # failed when it in fact succeeded. Report the quantification and
        # stop.
        print("\nBulk library: skipping the single-cell pipeline (QC, UMAP, "
              "Leiden all assume many cells; this is one sample).")
        if not args.dry_run:
            summarise_bulk_matrix(matrix, work)
        print("For differential expression across samples, quantify each "
              "sample and use `igvfagent rnaseq`.")
    elif matrix and not args.skip_analysis:
        sc = [_igvfagent(), "sc-analyze", "pipeline", "--input", str(matrix),
              "--label", label]
        print("Analysis: " + " ".join(sc))
        if not args.dry_run:
            subprocess.run(sc, check=False)

    write_progress(work, phase="done", percent=100.0,
                   detail=f"complete: {work}")
    print(f"Run dir: {work}")
    return 0


def summarise_bulk_matrix(matrix: Path, work: Path) -> None:
    """Print what the quantification produced, and save it as a TSV.

    A bulk run has no clustering to report, so without this the run ends
    having written a matrix and said nothing about it.
    """
    try:
        import anndata as ad
        import numpy as np
    except ImportError:
        print(f"Matrix written: {matrix} (install the analysis extra for a "
              f"summary)")
        return
    a = ad.read_h5ad(matrix)
    counts = np.asarray(a.X.sum(axis=0)).ravel()
    detected = int((counts > 0).sum())
    print(f"Samples x genes: {a.shape[0]} x {a.shape[1]}")
    print(f"Total counts:    {int(counts.sum()):,}")
    print(f"Genes detected:  {detected:,}")
    run_info = matrix.parent.parent / "run_info.json"
    if run_info.exists():
        try:
            info = json.loads(run_info.read_text())
            print(f"Reads processed: {info.get('n_processed', 0):,}  "
                  f"pseudoaligned: {info.get('p_pseudoaligned', '?')}%")
        except (ValueError, OSError):
            pass
    names = a.var.get("gene_name")
    labels = list(names) if names is not None else list(a.var_names)
    order = np.argsort(counts)[::-1][:15]
    out = work / "gene_counts.tsv"
    with out.open("w") as fh:
        fh.write("gene\tcount\n")
        for i in np.argsort(counts)[::-1]:
            if counts[i] <= 0:
                break
            fh.write(f"{labels[i]}\t{int(counts[i])}\n")
    print(f"Counts TSV:      {out}")
    print("Top genes:       " + ", ".join(
        f"{labels[i]} ({int(counts[i])})" for i in order[:8]))


def _igvfagent() -> str:
    return shutil.which("igvfagent") or sys.executable


# ─── Detached jobs ──────────────────────────────────────────────────────────
#
# Aligning a real dataset takes tens of minutes to hours: IGVFDS3532MONX is
# 45.6 GB of reads. The agent loop runs a tool as a blocking subprocess with
# no timeout, inside a web request, so a synchronous run holds the browser
# open for the whole job and the answer is lost when the socket drops. That
# is why "analyse this raw data" could only ever be planned, never done.
#
# `--detach` therefore starts the work in its own process, records a job
# file, and returns immediately with an id. `status` reads it back. The agent
# gets a fast tool call; the job outlives the conversation.

JOB_DIR = RUN_DIR / "_jobs"
# Reads are cached by FILE accession, not per run: a Portal file's bytes never
# change, so the second run of a dataset should not re-fetch 45 GB -- and,
# before this, a per-run directory also meant a second copy of it on disk.
FASTQ_CACHE = RUN_DIR / "_fastq"
# Completed runs, keyed by what actually determines the output, so a repeat
# of the same analysis can hand back the matrix instead of recomputing it.
RESULT_INDEX = RUN_DIR / "_results.json"


def result_key(accession: str, tech: str, workflow: str, reference: str) -> str:
    return "|".join([accession.upper(), (tech or "auto").upper(),
                     workflow or "standard", reference or "human"])


def record_result(key: str, matrix: Path, work: Path, extra: dict) -> None:
    try:
        idx = json.loads(RESULT_INDEX.read_text()) if RESULT_INDEX.exists() else {}
    except (OSError, ValueError):
        idx = {}
    idx[key] = {"matrix": str(matrix), "work": str(work),
                "completed": time.strftime("%Y-%m-%d %H:%M:%S"), **extra}
    try:
        RESULT_INDEX.parent.mkdir(parents=True, exist_ok=True)
        tmp = RESULT_INDEX.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(idx, indent=2, sort_keys=True, default=str))
        tmp.replace(RESULT_INDEX)
    except OSError:
        pass


def lookup_result(key: str) -> "Optional[dict]":
    """A previous completed run for the same analysis, if its matrix survives.

    The index is only a hint: the matrix it points at may have been deleted
    to reclaim space, so the file is checked before the entry is trusted.
    """
    try:
        idx = json.loads(RESULT_INDEX.read_text())
    except (OSError, ValueError):
        return None
    rec = idx.get(key)
    if not rec:
        return None
    if not Path(rec.get("matrix", "")).exists():
        return None
    return rec


def _job_path(job_id: str) -> Path:
    return JOB_DIR / f"{safe_label(job_id)}.json"


def write_job(job_id: str, **fields: Any) -> Path:
    JOB_DIR.mkdir(parents=True, exist_ok=True)
    path = _job_path(job_id)
    data = {}
    if path.exists():
        try:
            data = json.loads(path.read_text())
        except ValueError:
            data = {}
    data.update(fields)
    data["job_id"] = job_id
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str))
    return path


def _pid_alive(pid: Any, marker: str = "") -> bool:
    """Is this pid still our job?

    ``os.kill(pid, 0)`` alone is not enough: pids are recycled, so a long
    after a job died some unrelated process can inherit its number and the
    job is reported as still running forever. When ``marker`` is given the
    command line is checked too, which is what makes a "finished" state
    trustworthy.
    """
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    if not marker:
        return True
    try:
        cmdline = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\0", b" ").decode()
    except OSError:
        return True          # no procfs (macOS): fall back to the bare check
    return marker in cmdline


def write_progress(work: Path, **fields: Any) -> None:
    """Record machine-readable progress next to the log.

    A long job is otherwise a silent black box: the UI can only show that
    something is running, never how far along. Written atomically so a
    reader never sees a half-serialised file.
    """
    try:
        work.mkdir(parents=True, exist_ok=True)
        path = work / "progress.json"
        tmp = path.with_suffix(".json.tmp")
        fields["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        tmp.write_text(json.dumps(fields, indent=2, default=str))
        tmp.replace(path)
    except OSError:
        pass


def read_progress(work: Path) -> dict:
    try:
        return json.loads((work / "progress.json").read_text())
    except (OSError, ValueError):
        return {}


def render_bar(pct: float, width: int = 28) -> str:
    pct = max(0.0, min(100.0, float(pct)))
    filled = int(round(width * pct / 100.0))
    return "[" + "#" * filled + "-" * (width - filled) + f"] {pct:5.1f}%"


def _tail(path: Path, n: int = 12) -> "list[str]":
    """Last n meaningful lines: progress bars and INFO chatter are noise."""
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return []
    keep = [ln for ln in lines
            if ln.strip()
            and "B/s" not in ln
            and not re.match(r"^\s*\d+%", ln)
            and " INFO " not in ln]
    return keep[-n:]


def spawn_detached(argv: "list[str]", job_id: str, accession: str) -> dict:
    """Start the run in its own process and return its job record."""
    work = RUN_DIR / job_id
    work.mkdir(parents=True, exist_ok=True)
    log = work / "run.log"
    with log.open("w") as fh:
        proc = subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
    rec = {"accession": accession, "pid": proc.pid, "log": str(log),
           "work": str(work), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
           "argv": argv, "state": "running"}
    write_job(job_id, **rec)
    return rec


def cmd_guide_library(args: argparse.Namespace) -> int:
    """Report whether the sgRNA library for a screen is discoverable yet."""
    r = find_guide_library(args.accession)
    print(f"Accession:      {args.accession}")
    print(f"AnalysisSets:   {', '.join(r['analysis_sets']) or 'none'}")
    print(f"Libraries:      {', '.join(r['libraries']) or 'none'}")
    if r["resolved"]:
        g = r["guide_files"][0]
        print(f"RESOLVED:       yes — {r['why']}")
        idx = load_guide_index(g["id"] or g["accession"])
        print(f"Guide table:    {idx['n_rows']:,} constructs")
        # Report every candidate key's discriminating power. "1,741 guides"
        # alone is what made the LDLR pegRNA library look ready to count
        # when only 52 of its constructs could be told apart by spacer.
        for c in idx["candidates"]:
            print(f"  {c['column']:24} {c['distinct']:>7,} distinct "
                  f"({c['fraction']:>6.1%} of constructs), "
                  f"median {c['median_len']} bp")
        if idx["candidates"]:
            best = max(idx["candidates"], key=lambda c: c["distinct"])
            if best["distinct"] < idx["n_rows"] * 0.9:
                print(f"  NOTE: no column separates all constructs; the best "
                      f"({best['column']}) leaves "
                      f"{idx['n_rows'] - best['distinct']:,} indistinguishable.")
            print("  Which key is actually used is decided against real reads "
                  "at count time, since a key can be unique in the metadata "
                  "and absent from the reads.")
        # kite matches feature barcodes, which for a guide library means the
        # spacer. Recommending it for a prime-editing library would hand back
        # counts collapsed onto the handful of shared nick sites, so the
        # advice is conditional on the spacer actually separating constructs.
        spacer = next((c for c in idx["candidates"]
                       if c["column"] in ("spacer", "protospacer")), None)
        if spacer and spacer["distinct"] >= idx["n_rows"] * 0.9:
            print("\nReady for guide assignment: pass workflow='kite' to "
                  "raw_pipeline_run.")
        else:
            have = spacer["distinct"] if spacer else 0
            print(f"\nNOT suitable for kite: it matches feature barcodes, "
                  f"i.e. spacers, and the spacer separates only {have:,} of "
                  f"{idx['n_rows']:,} constructs here. This is a "
                  f"prime-editing-style library where the spacer sets the "
                  f"nick site and the edit lives elsewhere in the construct.")
            print("Use `raw-pipeline guide-count` or, for a sorted screen, "
                  "`crispr-screen analyze` -- both pick the counting key by "
                  "testing candidates against the reads.")
        return 0
    print(f"RESOLVED:       no")
    print(f"Reason:         {r['why']}")
    return 1


def cmd_assay_coverage(args: argparse.Namespace) -> int:
    """Which IGVF assays are classified, and what share of datasets."""
    from collections import Counter
    st, d = portal_json("/search/?type=MeasurementSet&limit=0&format=json")
    facets = {f.get("field"): f for f in (d or {}).get("facets") or []}
    terms = (facets.get("preferred_assay_titles") or {}).get("terms") or []
    if not terms:
        print("Could not read assay facets from the Portal.")
        return 1
    total = sum(t.get("doc_count", 0) for t in terms)
    by_route: "Counter[str]" = Counter()
    unknown = []
    for t in terms:
        key, n = t.get("key", ""), t.get("doc_count", 0)
        route = _A.ASSAY.get(key.strip().lower())
        if route:
            by_route[route] += n
        else:
            by_route["UNCLASSIFIED"] += n
            unknown.append((n, key))
    print(f"{len(terms)} assay titles, {total:,} MeasurementSets\n")
    for route, n in by_route.most_common():
        analysis, support = _A.ROUTE_GUIDANCE.get(route, ("", ""))
        print(f"  {n:>6} ({100.0*n/total:5.1f}%)  {route}")
        if support:
            print(f"           {support}")
    covered = total - by_route["UNCLASSIFIED"]
    print(f"\nClassified: {covered:,}/{total:,} "
          f"({100.0*covered/total:.1f}% of datasets)")
    print(f"Transcript assays: {by_route[_A.TRANSCRIPT]:,} "
          f"({100.0*by_route[_A.TRANSCRIPT]/total:.1f}%) — the only route "
          f"for which transcriptome quantification is correct.")
    if unknown:
        print("\nUNCLASSIFIED titles (these get an honest refusal, not a guess):")
        for n, k in sorted(unknown, reverse=True):
            print(f"  {n:>5}  {k}")
    return 0


def cmd_guide_count(args: argparse.Namespace) -> int:
    """Count guide occurrences in a gRNA-sequencing library.

    This is the analysis a CRISPR screen with gRNA-sequencing readout
    actually calls for. The reads are the guide library, so the measurement
    is how often each designed guide appears -- and, across sorted
    populations, how that frequency shifts.
    """
    import gzip

    lib = find_guide_library(args.accession)
    if not lib["resolved"]:
        print(f"No guide table reachable for {args.accession}: {lib['why']}")
        return 2
    gf = lib["guide_files"][0]
    # The counting key is chosen by measurement, not assumed to be the
    # spacer: a prime-editing library shares one spacer across every variant
    # it installs, so keying on spacer collapsed 1,741 constructs of the
    # LDLR library onto 52 and attributed each one's reads to an arbitrary
    # sibling. See load_guide_index / calibrate_key.
    idx = load_guide_index(gf["id"] or gf["accession"])
    if not idx["resolved"]:
        print(f"Guide table {gf['accession']}: {idx['why']}")
        return 2
    print(f"Guide library: {gf['library']} -> {gf['accession']} "
          f"({idx['n_rows']:,} constructs)")
    cal = calibrate_key(idx, args.accession, args.calibrate_reads)
    if cal["why"]:
        print(f"  {cal['why']}")
    for t in sorted(cal["tested"], key=lambda t: -t["rate"]):
        mark = ("  <- used" if cal["chosen"]
                and t["column"] == cal["chosen"]["column"] else "")
        print(f"  {t['column']:24} separates {t['fraction']:>6.1%} of "
              f"constructs, found in {t['rate']:>6.1%} of reads, "
              f"{t['ambiguous']:>5.1%} ambiguous{mark}")
    if not cal["chosen"]:
        print("\nNo candidate sequence column appears in these reads. Either "
              "they are not the construct amplicon, or the library table "
              "describes a different assay. Not guessing.")
        return 3
    key = cal["chosen"]
    matcher = build_matcher(key["seq_to_guide"])
    lengths = key["lengths"]
    guides = [(gid, sq) for sq, gid in key["seq_to_guide"].items()]
    print(f"Counting by:    {key['column']}  ({key['distinct']:,} distinct, "
          f"{min(lengths)}-{max(lengths)} bp)")

    files = [f for f in list_files(args.accession)
             if str(f.get("file_format", "")).lower() == "fastq"]
    total_gb = sum(gb(f.get("file_size")) for f in files)
    if args.max_download_gb and total_gb > args.max_download_gb:
        print(f"Reads total {total_gb} GB, over --max-download-gb "
              f"{args.max_download_gb}. Nothing downloaded.")
        return 4
    FASTQ_CACHE.mkdir(parents=True, exist_ok=True)
    local = []
    for f in files:
        name = Path(str(f.get("href") or f.get("accession"))).name
        dest = FASTQ_CACHE / name
        if dest.exists():
            print(f"Cached      {f.get('accession')} ({gb(f.get('file_size'))} GB)")
        else:
            print(f"Downloading {f.get('accession')} ({gb(f.get('file_size'))} GB)…")
            portal_download(f["href"], dest)
        local.append((f.get("accession"), dest))

    comp = str.maketrans("ACGTN", "TGCAN")
    counts: "Counter[str]" = Counter()
    per_file: "dict[str, Counter]" = {}
    stats = Counter()
    for acc, path in local:
        c: "Counter[str]" = Counter()
        op = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
        with op(path, "rt", errors="replace") as fh:      # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                if args.max_reads and stats["reads"] >= args.max_reads:
                    break
                stats["reads"] += 1
                seq = line.strip().upper()
                hit = None
                # Both orientations: a guide amplicon is sequenced from
                # either end depending on the primer, and checking one
                # silently halves the assignable reads.
                for s in (seq, seq.translate(comp)[::-1]):
                    hit = match_read(s, matcher)
                    if hit:
                        break
                if hit:
                    c[hit] += 1
                    stats["assigned"] += 1
                else:
                    stats["unassigned"] += 1
        per_file[str(acc)] = c
        counts.update(c)

    label = args.label or (f"{time.strftime('%Y%m%d_%H%M%S')}_"
                            f"{safe_label(args.accession)}_guides")
    out = REPORT_DIR.parent / "GuideCounts" / label
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    total = max(sum(counts.values()), 1)
    for gid, sq in guides:
        n = counts.get(gid, 0)
        rows.append({"guide_id": gid, key["column"]: sq, "count": n,
                     "freq": round(n / total, 8)})
    rows.sort(key=lambda r: -r["count"])
    with (out / "guide_counts.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh,
                            fieldnames=["guide_id", key["column"], "count", "freq"],
                            delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    detected = sum(1 for r in rows if r["count"] > 0)
    summary = {
        "accession": args.accession, "guide_library": gf["library"],
        "guide_table": gf["accession"], "guides_in_library": idx["n_rows"],
        "counting_key": key["column"],
        "constructs_distinguishable": key["distinct"],
        "read_length": cal["read_len"],
        "reads_examined": stats["reads"], "reads_assigned": stats["assigned"],
        "assignment_rate": round(stats["assigned"] / max(stats["reads"], 1), 4),
        "guides_detected": detected,
        "library_coverage": round(detected / max(len(guides), 1), 4),
        "per_file_assigned": {k: sum(v.values()) for k, v in per_file.items()},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"Guide counts: {out / 'guide_counts.tsv'}")
    print(f"Output: {out}")
    if stats["assigned"] == 0:
        print("\nNO reads matched any designed guide. The reads may not be "
              "this library, or they may carry a constant vector prefix that "
              "shifts the spacer -- check a few reads against the table "
              "before trusting a zero.")
        return 5
    print("\nNote: counts from ONE population are library composition, not "
          "enrichment. A FACS screen scores guides by comparing sorted "
          "against unsorted, so run the partner population too and compare "
          "the freq columns.")
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    JOB_DIR.mkdir(parents=True, exist_ok=True)
    jobs = sorted(JOB_DIR.glob("*.json"))
    if args.job:
        jobs = [_job_path(args.job)] if _job_path(args.job).exists() else []
        if not jobs:
            print(f"No such job: {args.job}")
            return 2
    if not jobs:
        print("No raw-pipeline jobs on record.")
        return 0
    for jp in jobs[-args.limit:]:
        try:
            rec = json.loads(jp.read_text())
        except ValueError:
            continue
        log = Path(rec.get("log", ""))
        alive = _pid_alive(rec.get("pid"), marker="raw-pipeline")
        done = log.exists() and any("Run dir:" in ln for ln in _tail(log, 40))
        state = "running" if alive else ("finished" if done else "stopped")
        if state != rec.get("state"):
            write_job(rec["job_id"], state=state)
        prog = read_progress(Path(rec.get("work", "")))
        print(f"\n=== {rec['job_id']} ===")
        print(f"accession: {rec.get('accession')}   started: {rec.get('started')}")
        print(f"state:     {state}   pid {rec.get('pid')}")
        if prog:
            pct = prog.get("percent")
            bar = render_bar(pct) if isinstance(pct, (int, float)) else "[ working ]"
            print(f"phase:     {prog.get('phase','?'):9} {bar}")
            if prog.get("detail"):
                print(f"           {prog['detail']}   (as of {prog.get('updated','?')})")
        print(f"log:       {log}")
        for ln in _tail(log, args.tail):
            print(f"  | {ln}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="raw_data_pipeline",
        description="Process an IGVF dataset's raw reads into an analysed matrix.")
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("accession", help="IGVF FileSet accession (IGVFDS…)")
        sp.add_argument("--technology", help="kb technology string (e.g. 10XV3). "
                                              "Inferred from the seqspec if omitted.")
        sp.add_argument("--reference", default="human",
                         help="Prebuilt kb reference to use (default: human).")
        sp.add_argument("--force-align", action="store_true",
                         help="Align the reads even if a published matrix exists.")
        sp.add_argument("--max-download-gb", type=float, default=100.0,
                         help="Refuse to download more than this (default 100).")
        return sp

    common(sub.add_parser("plan", help="Report the route and the cost; download nothing."))
    r = common(sub.add_parser("run", help="Execute the route end to end."))
    r.add_argument("--index"); r.add_argument("--t2g")
    r.add_argument("--workflow", default="standard",
                    choices=["standard", "nac", "kite", "kite:10xFB"],
                    help="kb workflow. 'kite' assigns feature barcodes (sgRNAs).")
    r.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    r.add_argument("--label")
    r.add_argument("--dry-run", action="store_true",
                    help="Print every command and download nothing.")
    r.add_argument("--no-reuse", action="store_true",
                    help="Recompute from the reads even if an identical "
                         "analysis was already completed on this machine.")
    r.add_argument("--detach", action="store_true",
                    help="Start the run in its own process and return an id "
                         "immediately. Use for anything large: a synchronous "
                         "run holds the caller open for the whole job.")
    r.add_argument("--skip-analysis", action="store_true",
                    help="Stop at the matrix; do not run the single-cell pipeline.")
    sub.add_parser("assay-coverage",
                    help="Which IGVF assays are classified, and what share "
                         "of Portal datasets each route covers.")

    gc = sub.add_parser("guide-count",
                         help="Count designed guides in a gRNA-sequencing "
                              "library.")
    gc.add_argument("accession")
    gc.add_argument("--max-reads", type=int, default=None)
    gc.add_argument("--max-download-gb", type=float, default=20.0)
    gc.add_argument("--calibrate-reads", type=int, default=20000,
                     help="Reads sampled to choose the counting key.")
    gc.add_argument("--label")

    gl = sub.add_parser("guide-library",
                         help="Is the sgRNA library for this screen "
                              "discoverable yet?")
    gl.add_argument("accession")

    st = sub.add_parser("status", help="Report on detached runs.")
    st.add_argument("job", nargs="?", help="Job id (default: all).")
    st.add_argument("--tail", type=int, default=12,
                     help="Log lines to show per job.")
    st.add_argument("--limit", type=int, default=5,
                     help="How many recent jobs to list.")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    for d in (RUN_DIR, REPORT_DIR, REF_DIR):
        d.mkdir(parents=True, exist_ok=True)
    if args.command == "plan":
        return cmd_plan(args)
    if args.command == "status":
        return cmd_status(args)
    if args.command == "guide-library":
        return cmd_guide_library(args)
    if args.command == "guide-count":
        return cmd_guide_count(args)
    if args.command == "assay-coverage":
        return cmd_assay_coverage(args)
    return cmd_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
