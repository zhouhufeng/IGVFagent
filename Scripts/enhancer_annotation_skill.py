"""Enhancer / regulatory-region annotation against the ENCODE cCRE registry.

Region-based, not variant-based. The existing ``ccre annotate-variants``
does point lookups: one coordinate, what overlaps it. An enhancer library
is a list of *intervals*, and the questions are different — how much of
this enhancer is covered by a cCRE, which classes, and is a region with
no overlap genuinely intergenic or just short of the nearest element.

    enhancer list (xlsx / bed / tsv / csv)
        -> interval overlap against the SCREEN cCRE registry
        -> per-region annotation + class breakdown + nearest-element fallback

The registry (https://screen.encodeproject.org, files served from the
Weng lab) is downloaded once and indexed into a local SQLite database, so
annotation runs offline afterwards and is reproducible against a pinned
registry version.

Why the reader is hand-written
------------------------------
``.xlsx`` is a zip of XML, so it is read here with ``zipfile`` +
``xml.etree`` rather than pulling in openpyxl. That keeps the skill
inside the project's standard-library-only contract, and it streams:
the reference library this was built against has a 1,048,534-row sheet
that a load-everything reader handles badly.

Subcommands
-----------
    build-db      Download the cCRE registry and index it locally
    db-stats      What the local cCRE database contains
    inspect       Show sheets/columns of an input file before annotating
    annotate      Annotate an enhancer/region list against cCREs
    write-playbook  Write Docs/Skills/ENHANCER_ANNOTATION_SKILLS.md
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import logging
import os
import re
import sqlite3
import sys
import time
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "EnhancerAnnotation"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
CCRE_DIR = DATA_DIR / "Reference" / "cCRE"
CCRE_DB = CCRE_DIR / "ccre.sqlite"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

WENGLAB = _resolve_endpoint("wenglab_dl", "WENGLAB_DL_BASE")
USER_AGENT = "IGVFagent-enhancer-annotation/0.1"

# The registry the user pins. V3 is what screen.encodeproject.org served
# as the stable Registry at time of writing; V4 is the newer build. Both
# are offered because an annotation must be reproducible against the
# version a study actually used.
REGISTRIES = {
    "V3": "Registry-V3/GRCh38-cCREs.bed",
    "V4": "Registry-V4/GRCh38-cCREs.bed",
}
DEFAULT_REGISTRY = "V3"

# Registry classification vocabulary, most-promoter-like first. Order is
# used to pick a single representative class for a region that overlaps
# several elements.
CLASS_RANK = ["PLS", "pELS", "dELS", "DNase-H3K4me3", "CA-CTCF", "CA-TF",
              "CA-only", "TF", "CTCF-only", "CTCF-bound", "Low-DNase"]

CHROM_COL = ("enhancer_chr", "chrom", "chr", "chromosome", "seqnames", "#chrom")
START_COL = ("enhancer_start", "start", "chromstart", "chrstart", "begin")
END_COL = ("enhancer_end", "end", "chromend", "chrend", "stop")
NAME_COL = ("enhancer_name", "name", "id", "element", "region", "enhancer")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"enhancer_annot_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)])
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in str(s))


def _run_dir(label: str) -> Path:
    d = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def opener(path, mode="rt"):
    return gzip.open(str(path), mode) if str(path).endswith(".gz") else open(path, mode)


# ─── Minimal streaming .xlsx reader (stdlib only) ───────────────────────────

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col_index(ref: str) -> int:
    """'BC12' -> 54 (0-based column). Cells can be sparse, so position must
    come from the reference, not from arrival order."""
    n = 0
    for ch in ref:
        if ch.isalpha():
            n = n * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return n - 1


def xlsx_sheets(path) -> "list[str]":
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        return [s.get("name") for s in wb.iter(f"{_NS}sheet")]


def _shared_strings(z) -> "list[str]":
    try:
        raw = z.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    out = []
    for si in ET.fromstring(raw).iter(f"{_NS}si"):
        out.append("".join(t.text or "" for t in si.iter(f"{_NS}t")))
    return out


def read_xlsx(path, sheet: str | None = None, max_rows: int = 0
              ) -> Iterator["list[str]"]:
    """Stream one worksheet as lists of cell strings."""
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = [(s.get("name"), s.get(
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"))
            for s in wb.iter(f"{_NS}sheet")]
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = {r.get("Id"): r.get("Target") for r in rels}
        name = sheet or sheets[0][0]
        match = [(n, rid) for n, rid in sheets if n == name]
        if not match:
            raise SystemExit(f"Sheet {name!r} not found. Available: "
                             + ", ".join(n for n, _ in sheets))
        tgt = target[match[0][1]].lstrip("/")
        if not tgt.startswith("xl/"):
            tgt = "xl/" + tgt
        strings = _shared_strings(z)
        with z.open(tgt) as fh:
            row: dict[int, str] = {}
            width = 0
            for ev, el in ET.iterparse(fh, events=("end",)):
                if el.tag == f"{_NS}c":
                    ref = el.get("r") or ""
                    ci = _col_index(ref) if ref else len(row)
                    t = el.get("t")
                    v = el.find(f"{_NS}v")
                    if t == "inlineStr":
                        node = el.find(f"{_NS}is")
                        txt = "".join(x.text or "" for x in node.iter(f"{_NS}t")) if node is not None else ""
                    elif v is None:
                        txt = ""
                    elif t == "s":
                        try:
                            txt = strings[int(v.text)]
                        except (ValueError, IndexError, TypeError):
                            txt = ""
                    else:
                        txt = v.text or ""
                    row[ci] = txt
                    width = max(width, ci + 1)
                    el.clear()
                elif el.tag == f"{_NS}row":
                    yield [row.get(i, "") for i in range(width)]
                    row = {}
                    el.clear()
                    if max_rows:
                        max_rows -= 1
                        if max_rows <= 0:
                            return


# ─── Region input ───────────────────────────────────────────────────────────

class Region:
    __slots__ = ("name", "chrom", "start", "end", "extra")

    def __init__(self, name, chrom, start, end, extra=None):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.extra = extra or {}

    @property
    def length(self) -> int:
        return self.end - self.start

    def key(self):
        return (self.chrom, self.start, self.end)


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _pick(header: "list[str]", candidates) -> int | None:
    norm = [_norm(h) for h in header]
    for cand in candidates:
        c = _norm(cand)
        for i, h in enumerate(norm):
            if h == c:
                return i
    for cand in candidates:          # substring fallback
        c = _norm(cand)
        for i, h in enumerate(norm):
            if c and c in h:
                return i
    return None


def _norm_chrom(c: str) -> str:
    c = str(c).strip()
    if not c:
        return ""
    if not c.lower().startswith("chr"):
        c = "chr" + c
    return c.replace("chrMT", "chrM")


def rows_from(path: str, sheet: str | None, max_rows: int = 0) -> Iterator["list[str]"]:
    p = str(path).lower()
    if p.endswith((".xlsx", ".xlsm")):
        yield from read_xlsx(path, sheet, max_rows)
        return
    delim = "\t"
    if p.endswith((".csv", ".csv.gz")):
        delim = ","
    with opener(path) as fh:
        n = 0
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            yield line.split(delim)
            n += 1
            if max_rows and n >= max_rows:
                return


def read_regions(path: str, *, sheet: str | None = None,
                 chrom_col: str | None = None, start_col: str | None = None,
                 end_col: str | None = None, name_col: str | None = None,
                 dedupe: bool = True) -> "tuple[list[Region], dict]":
    """Read regions from xlsx / bed / tsv / csv, detecting columns by name.

    A headerless BED is detected by its shape (col2/col3 numeric) rather
    than assumed, so both a plain BED and a spreadsheet with named columns
    work without the caller declaring which it is.
    """
    it = rows_from(path, sheet)
    try:
        first = next(it)
    except StopIteration:
        raise SystemExit(f"{path}: no rows")

    headerless = (len(first) >= 3 and str(first[1]).strip().isdigit()
                  and str(first[2]).strip().isdigit())
    if headerless:
        ci, si, ei = 0, 1, 2
        ni = 3 if len(first) > 3 else None
        header = ["chrom", "start", "end", "name"][:max(3, len(first))]
        rows = [first]
    else:
        header = [str(h) for h in first]
        ci = _pick(header, [chrom_col] if chrom_col else CHROM_COL)
        si = _pick(header, [start_col] if start_col else START_COL)
        ei = _pick(header, [end_col] if end_col else END_COL)
        ni = _pick(header, [name_col] if name_col else NAME_COL)
        rows = []
        missing = [n for n, v in (("chrom", ci), ("start", si), ("end", ei)) if v is None]
        if missing:
            raise SystemExit(
                f"Could not find column(s) {', '.join(missing)} in {path}"
                + (f" sheet {sheet!r}" if sheet else "")
                + f".\nHeader was: {', '.join(header[:14])}"
                + "\nName them explicitly with --chrom-col / --start-col / --end-col.")

    stats = {"rows_read": 0, "skipped_bad_coords": 0, "skipped_no_chrom": 0}
    seen: dict[tuple, Region] = {}
    out: "list[Region]" = []

    def consume(r):
        stats["rows_read"] += 1
        try:
            chrom = _norm_chrom(r[ci])
            start = int(float(str(r[si]).replace(",", "")))
            end = int(float(str(r[ei]).replace(",", "")))
        except (ValueError, IndexError, TypeError):
            stats["skipped_bad_coords"] += 1
            return
        if not chrom or end <= start:
            stats["skipped_no_chrom" if not chrom else "skipped_bad_coords"] += 1
            return
        name = (str(r[ni]).strip() if ni is not None and ni < len(r) and str(r[ni]).strip()
                else f"{chrom}:{start}-{end}")
        reg = Region(name, chrom, start, end)
        if dedupe:
            k = reg.key()
            if k in seen:
                return
            seen[k] = reg
        out.append(reg)

    for r in rows:
        consume(r)
    for r in it:
        consume(r)

    stats["regions"] = len(out)
    stats["header"] = header[:14]
    return out, stats


# ─── cCRE registry: download + local database ───────────────────────────────

def registry_url(version: str) -> str:
    if version not in REGISTRIES:
        raise SystemExit(f"Unknown registry {version!r}. Choose: {', '.join(REGISTRIES)}")
    return f"{WENGLAB}/{REGISTRIES[version]}"


def download_registry(version: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = registry_url(version)
    logging.info("downloading cCRE registry %s", version)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        got = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if total and got % (16 << 20) < (1 << 20):
                logging.info("  %.0f%%  (%.1f MB)", 100 * got / total, got / 1e6)
    tmp.replace(dest)
    logging.info("saved %s (%.1f MB)", dest, dest.stat().st_size / 1e6)
    return dest


def build_db(version: str, bed: Path, db: Path) -> dict:
    db.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(db))
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("""CREATE TABLE IF NOT EXISTS ccre(
        registry TEXT NOT NULL, chrom TEXT NOT NULL,
        start INTEGER NOT NULL, end INTEGER NOT NULL,
        accession TEXT, rdhs TEXT, classes TEXT)""")
    con.execute("""CREATE TABLE IF NOT EXISTS ccre_meta(
        registry TEXT PRIMARY KEY, url TEXT, n_rows INTEGER,
        built_utc TEXT, source_bytes INTEGER)""")
    con.execute("DELETE FROM ccre WHERE registry=?", (version,))
    n = 0
    batch = []
    with opener(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            batch.append((version, _norm_chrom(f[0]), start, end,
                          f[4] if len(f) > 4 else (f[3] if len(f) > 3 else ""),
                          f[3] if len(f) > 3 else "",
                          f[5] if len(f) > 5 else ""))
            n += 1
            if len(batch) >= 50000:
                con.executemany("INSERT INTO ccre VALUES (?,?,?,?,?,?,?)", batch)
                batch.clear()
    if batch:
        con.executemany("INSERT INTO ccre VALUES (?,?,?,?,?,?,?)", batch)
    con.execute("CREATE INDEX IF NOT EXISTS idx_ccre_pos ON ccre(registry,chrom,start)")
    con.execute("INSERT OR REPLACE INTO ccre_meta VALUES (?,?,?,?,?)",
                (version, registry_url(version), n,
                 time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                 bed.stat().st_size))
    con.commit()
    con.close()
    return {"registry": version, "rows": n, "db": str(db)}


class CcreIndex:
    """Per-chromosome sorted intervals with a max-end prefix, so an overlap
    query is two binary searches instead of a scan of 2.3 M rows."""

    def __init__(self, db: Path, version: str, chroms=None):
        if not db.exists():
            raise SystemExit(
                f"No local cCRE database at {db}.\n"
                f"Build it first:  igvfagent enhancer build-db --registry {version}")
        con = sqlite3.connect(str(db))
        q = "SELECT chrom,start,end,accession,classes FROM ccre WHERE registry=?"
        params = [version]
        if chroms:
            q += " AND chrom IN (%s)" % ",".join("?" * len(chroms))
            params += list(chroms)
        q += " ORDER BY chrom,start"
        self.by_chrom: dict[str, dict] = {}
        for chrom, start, end, acc, classes in con.execute(q, params):
            d = self.by_chrom.setdefault(chrom, {"s": [], "e": [], "a": [], "c": [], "m": []})
            d["s"].append(start); d["e"].append(end)
            d["a"].append(acc); d["c"].append(classes)
        con.close()
        self.n = 0
        for d in self.by_chrom.values():
            run = 0
            for e in d["e"]:
                run = max(run, e)
                d["m"].append(run)
            self.n += len(d["s"])
        if not self.n:
            raise SystemExit(
                f"cCRE database has no rows for registry {version}. "
                f"Run:  igvfagent enhancer build-db --registry {version}")

    def overlaps(self, chrom, start, end):
        d = self.by_chrom.get(chrom)
        if not d:
            return []
        s = d["s"]
        hi = bisect.bisect_right(s, end)          # candidates start before region end
        out = []
        for i in range(hi - 1, -1, -1):
            if d["m"][i] <= start:                # no earlier interval can reach
                break
            if d["e"][i] > start:
                out.append((s[i], d["e"][i], d["a"][i], d["c"][i]))
        return out[::-1]

    def nearest(self, chrom, start, end):
        d = self.by_chrom.get(chrom)
        if not d:
            return None
        s = d["s"]
        i = bisect.bisect_left(s, start)
        best = None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(s):
                dist = 0 if (s[j] < end and d["e"][j] > start) else (
                    s[j] - end if s[j] >= end else start - d["e"][j])
                if best is None or dist < best[0]:
                    best = (max(0, dist), d["a"][j], d["c"][j])
        return best


# ─── Annotation ─────────────────────────────────────────────────────────────

def split_classes(s: str) -> "list[str]":
    return [c.strip() for c in str(s).split(",") if c.strip()]


def top_class(classes) -> str:
    for c in CLASS_RANK:
        if c in classes:
            return c
    return sorted(classes)[0] if classes else ""


def annotate_regions(regions, idx: CcreIndex, *, nearest: bool = True) -> "list[dict]":
    out = []
    for r in regions:
        hits = idx.overlaps(r.chrom, r.start, r.end)
        classes, accs, cov = [], [], 0
        spans = []
        for s, e, acc, cl in hits:
            ov = min(e, r.end) - max(s, r.start)
            if ov <= 0:
                continue
            spans.append((max(s, r.start), min(e, r.end)))
            accs.append(acc)
            classes.extend(split_classes(cl))
        # union of covered bases, so overlapping elements are not double-counted
        spans.sort()
        merged = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        cov = sum(b - a for a, b in merged)
        uniq = sorted(set(classes), key=lambda c: CLASS_RANK.index(c) if c in CLASS_RANK else 99)
        row = {
            "name": r.name, "chrom": r.chrom, "start": r.start, "end": r.end,
            "length": r.length,
            "n_ccre": len(accs),
            "ccre_accessions": ";".join(accs[:24]),
            "ccre_classes": ";".join(uniq),
            "top_class": top_class(uniq),
            "ccre_overlap_bp": cov,
            "ccre_overlap_frac": round(cov / r.length, 4) if r.length else 0.0,
            "nearest_ccre": "", "nearest_distance_bp": "",
        }
        if not accs and nearest:
            nb = idx.nearest(r.chrom, r.start, r.end)
            if nb:
                row["nearest_ccre"] = nb[1]
                row["nearest_distance_bp"] = nb[0]
                row["ccre_classes"] = ""
        out.append(row)
    return out


def summarize(rows, regions) -> dict:
    n = len(rows)
    withc = [r for r in rows if r["n_ccre"] > 0]
    cls = {}
    for r in withc:
        for c in r["ccre_classes"].split(";"):
            if c:
                cls[c] = cls.get(c, 0) + 1
    top = {}
    for r in withc:
        if r["top_class"]:
            top[r["top_class"]] = top.get(r["top_class"], 0) + 1
    fracs = sorted(r["ccre_overlap_frac"] for r in withc)
    lens = sorted(r["length"] for r in rows)

    def med(v):
        return v[len(v) // 2] if v else 0
    dists = sorted(int(r["nearest_distance_bp"]) for r in rows
                   if r["nearest_distance_bp"] != "")
    return {
        "regions": n,
        "regions_with_ccre": len(withc),
        "pct_with_ccre": round(100 * len(withc) / n, 2) if n else 0.0,
        "regions_without_ccre": n - len(withc),
        "median_region_length_bp": med(lens),
        "median_overlap_fraction": med(fracs),
        "mean_ccres_per_region": round(sum(r["n_ccre"] for r in rows) / n, 2) if n else 0,
        "class_counts_any": dict(sorted(cls.items(), key=lambda kv: -kv[1])),
        "top_class_counts": dict(sorted(top.items(), key=lambda kv: -kv[1])),
        "median_distance_to_nearest_bp_when_no_overlap": med(dists),
    }


def write_rows(path: Path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


COLS = ["name", "chrom", "start", "end", "length", "n_ccre", "top_class",
        "ccre_classes", "ccre_overlap_bp", "ccre_overlap_frac",
        "ccre_accessions", "nearest_ccre", "nearest_distance_bp"]


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_build_db(args) -> int:
    setup_logging()
    bed = CCRE_DIR / f"GRCh38-cCREs.{args.registry}.bed"
    if bed.exists() and not args.force:
        logging.info("using cached %s (%.1f MB) — pass --force to re-download",
                     bed, bed.stat().st_size / 1e6)
    else:
        download_registry(args.registry, bed)
    info = build_db(args.registry, bed, CCRE_DB)
    print(f"Registry:  {info['registry']}  ({registry_url(args.registry)})")
    print(f"Elements:  {info['rows']:,}")
    print(f"Database:  {info['db']}")
    print(f"Cached BED: {bed}")
    return 0


def cmd_db_stats(args) -> int:
    if not CCRE_DB.exists():
        print(f"No cCRE database yet at {CCRE_DB}.")
        print("Build one:  igvfagent enhancer build-db")
        return 1
    con = sqlite3.connect(str(CCRE_DB))
    print(f"Database: {CCRE_DB}  ({CCRE_DB.stat().st_size/1e6:.1f} MB)")
    for reg, url, n, built, sz in con.execute(
            "SELECT registry,url,n_rows,built_utc,source_bytes FROM ccre_meta"):
        print(f"\n  registry {reg}: {n:,} elements   built {built}")
        print(f"    source: {url}  ({sz/1e6:.1f} MB)")
        rows = con.execute(
            "SELECT classes, COUNT(*) c FROM ccre WHERE registry=? "
            "GROUP BY classes ORDER BY c DESC LIMIT 8", (reg,)).fetchall()
        print("    most common class combinations:")
        for cl, c in rows:
            print(f"      {cl or '(none)':32s} {c:>9,}")
    con.close()
    return 0


def cmd_inspect(args) -> int:
    p = str(args.input).lower()
    print(f"File: {args.input}")
    if p.endswith((".xlsx", ".xlsm")):
        sheets = xlsx_sheets(args.input)
        print(f"Sheets ({len(sheets)}):")
        for s in sheets:
            try:
                head = next(read_xlsx(args.input, s, max_rows=1))
            except StopIteration:
                head = []
            ci = _pick(head, CHROM_COL); si = _pick(head, START_COL); ei = _pick(head, END_COL)
            ok = "coordinates detected" if None not in (ci, si, ei) else "no coordinate columns"
            print(f"  {s:24s} {len(head):>3} cols  — {ok}")
            if head:
                print(f"      {', '.join(str(h)[:22] for h in head[:9])}")
    else:
        rows = rows_from(args.input, None, max_rows=3)
        for i, r in enumerate(rows):
            print(("  header: " if i == 0 else "  row:    ") + " | ".join(str(x)[:20] for x in r[:9]))
    return 0


def cmd_annotate(args) -> int:
    setup_logging()
    regions, rstats = read_regions(
        args.input, sheet=args.sheet, chrom_col=args.chrom_col,
        start_col=args.start_col, end_col=args.end_col, name_col=args.name_col,
        dedupe=not args.keep_duplicates)
    if not regions:
        raise SystemExit("No usable regions were read. Try `enhancer inspect` first.")
    chroms = sorted({r.chrom for r in regions})
    idx = CcreIndex(CCRE_DB, args.registry, chroms=chroms)
    rows = annotate_regions(regions, idx, nearest=not args.no_nearest)
    summ = summarize(rows, regions)
    summ["input"] = str(args.input)
    summ["sheet"] = args.sheet
    summ["registry"] = args.registry
    summ["cCRE_elements_loaded"] = idx.n
    summ.update({k: v for k, v in rstats.items() if k != "header"})

    out_dir = _run_dir(args.label or Path(args.input).stem)
    write_rows(out_dir / "annotated_regions.csv", rows, COLS)
    write_rows(out_dir / "regions_without_ccre.csv",
               [r for r in rows if r["n_ccre"] == 0], COLS)
    (out_dir / "summary.json").write_text(json.dumps(summ, indent=2))
    with open(out_dir / "annotated_regions.bed", "w") as fh:
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\t"
                     f"{int(round(1000*r['ccre_overlap_frac']))}\t.\t"
                     f"{r['top_class'] or 'none'}\n")

    print(f"Input:            {args.input}" + (f"  [sheet {args.sheet}]" if args.sheet else ""))
    print(f"Rows read:        {rstats['rows_read']:,}")
    print(f"Unique regions:   {summ['regions']:,}"
          + ("" if args.keep_duplicates else "  (deduplicated on coordinates)"))
    if rstats["skipped_bad_coords"] or rstats["skipped_no_chrom"]:
        print(f"  skipped:        {rstats['skipped_bad_coords']:,} bad coords, "
              f"{rstats['skipped_no_chrom']:,} no chromosome")
    print(f"cCRE registry:    {args.registry}  ({idx.n:,} elements on {len(chroms)} chromosomes)")
    print(f"With a cCRE:      {summ['regions_with_ccre']:,}  ({summ['pct_with_ccre']}%)")
    print(f"Without:          {summ['regions_without_ccre']:,}"
          + (f"  (median {summ['median_distance_to_nearest_bp_when_no_overlap']:,} bp to nearest)"
             if summ['median_distance_to_nearest_bp_when_no_overlap'] else ""))
    print(f"Median overlap:   {summ['median_overlap_fraction']:.0%} of the region")
    if summ["top_class_counts"]:
        print("Representative class per region:")
        for c, n in list(summ["top_class_counts"].items())[:8]:
            print(f"  {c:16s} {n:>7,}  ({100*n/summ['regions']:.1f}%)")
    print(f"Saved:            {out_dir}")
    return 0


PLAYBOOK = """# Enhancer / regulatory-region annotation — skill playbook

Annotates a list of **regions** against the ENCODE SCREEN candidate
cis-regulatory element (cCRE) registry. This is interval overlap, not
point lookup: `ccre annotate-variants` answers "what overlaps this
coordinate", while this answers "how much of this enhancer is covered by
a cCRE, of which classes, and if none — how far is the nearest".

## Data source

The registry is served from the Weng lab for
[SCREEN](https://screen.encodeproject.org):

| Registry | File |
|---|---|
| `V3` (default) | `Registry-V3/GRCh38-cCREs.bed` |
| `V4` | `Registry-V4/GRCh38-cCREs.bed` |

Downloaded once into `Data/Reference/cCRE/` and indexed into
`Data/Reference/cCRE/ccre.sqlite`, so annotation runs offline afterwards
and is reproducible against a pinned registry version. The endpoint is
resolved at runtime like every other archive — no URL is embedded in
source.

## Usage

```bash
# once — download and index the registry (~64 MB, ~2.3 M elements)
igvfagent enhancer build-db --registry V3
igvfagent enhancer db-stats

# look before you leap: which sheet holds coordinates?
igvfagent enhancer inspect --input Data/Input/EnhancerList/library.xlsx

# annotate
igvfagent enhancer annotate \\
    --input Data/Input/EnhancerList/library.xlsx \\
    --sheet High.sgRNA0820 --label high_enhancers
```

## Input formats

`.xlsx` / `.xlsm`, BED, TSV, CSV, optionally gzipped. Coordinate columns
are detected by name (`Enhancer_chr` / `chrom` / `chr`, `Enhancer_start` /
`start` / `chromStart`, and so on) and a headerless BED is recognised by
shape. Override with `--chrom-col` / `--start-col` / `--end-col` /
`--name-col` when detection is wrong.

`.xlsx` is read with `zipfile` + `xml.etree`, streaming, so no third-party
spreadsheet dependency is required and million-row sheets do not have to
be held in memory.

A CRISPR library lists one row per sgRNA, so the same enhancer repeats
several times. Rows are deduplicated on coordinates by default; pass
`--keep-duplicates` to annotate every row.

## Outputs

Written to `Docs/EnhancerAnnotation/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `annotated_regions.csv` | one row per region with the columns below |
| `regions_without_ccre.csv` | the subset with no overlap, for follow-up |
| `annotated_regions.bed` | BED with overlap fraction as score, class as name |
| `summary.json` | counts, class distribution, coverage statistics |

| Column | Meaning |
|---|---|
| `n_ccre` | number of cCREs overlapping the region |
| `top_class` | one representative class, ranked PLS > pELS > dELS > … |
| `ccre_classes` | every class seen across overlapping elements |
| `ccre_overlap_bp` | bases covered, **union** — overlapping elements are not double-counted |
| `ccre_overlap_frac` | that coverage as a fraction of the region |
| `nearest_ccre` / `nearest_distance_bp` | filled only when nothing overlaps |

## Interpreting the result

A region with no overlapping cCRE is not necessarily inert — it may be a
true negative, or the registry may simply lack coverage in that cell
type. `nearest_distance_bp` is what separates "sits in a cCRE desert"
from "just missed one", so both are reported rather than collapsing the
no-overlap set into a single count.

`top_class` is a convenience for summarising; a region overlapping both a
promoter-like and a distal-enhancer-like element is genuinely both, and
`ccre_classes` keeps that.
"""


def cmd_write_playbook(args) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = SKILL_DOC_DIR / "ENHANCER_ANNOTATION_SKILLS.md"
    out.write_text(PLAYBOOK)
    print(f"Wrote {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent enhancer",
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command")

    q = sub.add_parser("build-db", help="Download and index the cCRE registry.")
    q.add_argument("--registry", default=DEFAULT_REGISTRY, choices=list(REGISTRIES))
    q.add_argument("--force", action="store_true", help="Re-download even if cached.")
    q.set_defaults(func=cmd_build_db)

    q = sub.add_parser("db-stats", help="Show what the local cCRE database holds.")
    q.set_defaults(func=cmd_db_stats)

    q = sub.add_parser("inspect", help="Show sheets/columns before annotating.")
    q.add_argument("--input", required=True)
    q.set_defaults(func=cmd_inspect)

    q = sub.add_parser("annotate", help="Annotate a region list against cCREs.")
    q.add_argument("--input", required=True, help="xlsx / bed / tsv / csv (.gz ok).")
    q.add_argument("--sheet", default=None, help="Worksheet name for xlsx input.")
    q.add_argument("--registry", default=DEFAULT_REGISTRY, choices=list(REGISTRIES))
    q.add_argument("--chrom-col", default=None)
    q.add_argument("--start-col", default=None)
    q.add_argument("--end-col", default=None)
    q.add_argument("--name-col", default=None)
    q.add_argument("--keep-duplicates", action="store_true",
                   help="Annotate every row instead of unique coordinates.")
    q.add_argument("--no-nearest", action="store_true",
                   help="Skip the nearest-element fallback for non-overlapping regions.")
    q.add_argument("--label", default=None)
    q.set_defaults(func=cmd_annotate)

    q = sub.add_parser("write-playbook",
                       help="Write Docs/Skills/ENHANCER_ANNOTATION_SKILLS.md.")
    q.set_defaults(func=cmd_write_playbook)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
