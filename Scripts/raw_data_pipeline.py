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
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlsplit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint          # noqa: E402
from _credentials import portal_credentials as _portal_credentials  # noqa: E402

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


def portal_download(href: str, dest: Path) -> Path:
    """Download one Portal file, following the S3 redirect safely."""
    url = href if href.startswith("http") else f"{IGVF_API_BASE}{href}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=_headers())
    with _OPENER.open(req, timeout=600) as r, tmp.open("wb") as fh:
        shutil.copyfileobj(r, fh, length=4 * 1024 * 1024)
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

def decide_route(accession: str, buckets: "dict[str, list[dict]]",
                  force_align: bool) -> "dict":
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

    route = decide_route(accession, buckets, force_align)
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
        if single_end:
            for f in inv["buckets"]["reads"]:
                name = Path(str(f.get("href") or f.get("accession"))).name
                dest = fq_dir / name
                if args.dry_run:
                    print(f"DRY RUN: would download {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB) -> {dest}")
                elif dest.exists():
                    logging.info("already present: %s", dest)
                else:
                    print(f"Downloading {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB)…")
                    portal_download(f["href"], dest)
                local_singles.append(dest)
        for r1, r2 in pairs:
            p = []
            for f in (r1, r2):
                name = Path(str(f.get("href") or f.get("accession"))).name
                dest = fq_dir / name
                if args.dry_run:
                    print(f"DRY RUN: would download {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB) -> {dest}")
                elif dest.exists():
                    logging.info("already present: %s", dest)
                else:
                    print(f"Downloading {f.get('accession')} "
                          f"({gb(f.get('file_size'))} GB)…")
                    portal_download(f["href"], dest)
                p.append(dest)
            local_pairs.append((p[0], p[1]))

        kb_out = work / "kb"
        cmd = kb_count_cmd(index, t2g, tech, kb_out, local_pairs,
                            args.threads, args.workflow,
                            singles=local_singles or None)
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
    else:
        print("Nothing to process: no reads and no matrix.")
        return 2

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


def _pid_alive(pid: Any) -> bool:
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


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
        alive = _pid_alive(rec.get("pid"))
        done = log.exists() and any("Run dir:" in ln for ln in _tail(log, 40))
        state = "running" if alive else ("finished" if done else "stopped")
        if state != rec.get("state"):
            write_job(rec["job_id"], state=state)
        print(f"\n=== {rec['job_id']} ===")
        print(f"accession: {rec.get('accession')}   started: {rec.get('started')}")
        print(f"state:     {state}   pid {rec.get('pid')}")
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
    r.add_argument("--detach", action="store_true",
                    help="Start the run in its own process and return an id "
                         "immediately. Use for anything large: a synchronous "
                         "run holds the caller open for the whole job.")
    r.add_argument("--skip-analysis", action="store_true",
                    help="Stop at the matrix; do not run the single-cell pipeline.")
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
    return cmd_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
