"""ChIP-Atlas (Ohta/Oki) canonical-query skill.

A lightweight, anonymous client for the [ChIP-Atlas](https://chip-atlas.org)
public REST + bulk archive: reprocessed ChIP-seq / ATAC-seq / DNase-seq /
Bisulfite-seq peak calls and metadata for hundreds of thousands of public
SRA experiments, plus the WABI Enrichment Analysis / Diff Analysis job
service at NIG.

Clean-room reimplementation under Apache-2.0 of the public ChIP-Atlas
HTTP surface, modelled on the canonical MCP-server reference at
[inutano/chip-atlas](https://github.com/inutano/chip-atlas) (MIT, Tazro
Inutano Ohta / Shinya Oki / DBCLS, 2017–2026). Stdlib-only — no `httpx` /
`mcp` / `pydantic` runtime dep — and uses three indirected hosts via
``Scripts/_endpoints.py``:

  * ``chip-atlas.org`` — JSON browse / search / POST-download
  * ``chip-atlas.dbcls.jp/data`` — bulk static archive (BigBeds, BigWigs,
    assembled BEDs, ``experimentList.tab``)
  * ``dtn1.ddbj.nig.ac.jp/wabi/chipatlas`` — WABI Enrichment / Diff job
    queue

Conventions
-----------
* ``genome`` ∈ {hg38, hg19, mm10, mm9, rn6, dm6, dm3, ce11, ce10, sacCer3}
* ``agClass`` ∈ {Histone, "TFs and others", "RNA polymerase",
                  "Input control", ATAC-Seq, DNase-seq, Bisulfite-Seq,
                  "Annotation tracks", "Unclassified"}
* ``qval`` ∈ {"05", "10", "20", "50"} — peak-call -log10(q) threshold;
   lower index → looser peaks (more rows), higher index → stricter.
* ``clClass`` / ``clSubClass`` form a 2-level cell-type tree exposed via
   ``/data/sample_types`` and ``/data/chip_antigen`` browser endpoints.
* Per-experiment files live at
   ``…/{genome}/eachData/{bb,bw,bb05,bb10,bb20,bed05,bed10,bed20}/{SRX}.{qval}.{ext}``.

Commands
--------
    chipatlas list-genomes        Supported reference assemblies.
    chipatlas list-qvalues        Allowed -log10(q) thresholds.
    chipatlas list-experiment-types  Browse: agClass × cellClass tree.
    chipatlas list-antigens       Browse: antigens by genome / agClass / cellClass.
    chipatlas list-cell-types     Browse: cell types by genome / agClass.
    chipatlas search              Free-text search over experiments.
    chipatlas get-experiment      Full metadata for one SRX accession.
    chipatlas download-experiment Per-experiment BigBed / BigWig + per-qval BED.
    chipatlas assemble-bed        POST /download → assembled all-peaks BED URL
                                    for a (genome × ag × cellClass × qval) tuple.
    chipatlas download-all-peaks  Bulk allPeaks_light.{genome}.{qval}.bed.gz.
    chipatlas target-genes        Discover antigens with pre-computed Target-Genes
                                    tables, or pull a specific one.
    chipatlas submit-enrichment   POST a gene/region list to the WABI job queue.
    chipatlas poll-enrichment     Poll a WABI job for results / status.
    chipatlas showcase            End-to-end demo over the CTCF / H3K4me3 axis.
    chipatlas write-playbook      Write Docs/Skills/CHIPATLAS_SKILL.md.

License posture
---------------
Apache-2.0. Upstream MCP-server code is MIT (Tazro Ohta) — compatible with
verbatim copy + attribution, but we re-derive the wire contract to stay
stdlib-only. The ChIP-Atlas **data** themselves are licensed by NBDC/DBCLS
(see `https://dbarchive.biosciencedbc.jp/en/chip-atlas/lic.html` — typically
CC-BY 4.0) — we only fetch / link, never redistribute. Cite:
Zou/Ohta/Oki, *Nucleic Acids Res.* 2024 (doi:10.1093/nar/gkae358); Oki
et al., *EMBO Rep.* 2018 (doi:10.15252/embr.201846255).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import sys
import time
import urllib.error
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
REPORT_DIR = DOCS_DIR / "ChIPAtlas"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
SOURCES_DIR = DATA_DIR / "Sources" / "ChIPAtlas"
CACHE_DIR = DATA_DIR / "Cache" / "ChIPAtlas"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

CHIPATLAS_API  = _resolve_endpoint("chipatlas_api",  "CHIPATLAS_API")
CHIPATLAS_DATA = _resolve_endpoint("chipatlas_data", "CHIPATLAS_DATA")
CHIPATLAS_WABI = _resolve_endpoint("chipatlas_wabi", "CHIPATLAS_WABI")
USER_AGENT = "IGVFagent-chipatlas/0.1"

# Sensible polite-client defaults so we don't hammer a small group's server.
_MIN_SECONDS_BETWEEN_REQUESTS = 1.0
_last_request_time: float = 0.0


# Canonical lookups (cached at import time on demand). The upstream
# `list_of_genome.json` and `qval_range.json` endpoints are tiny and
# stable, so a small in-memory cache here means a fresh fork only pays
# for them once per process.

_CANONICAL_AG_CLASSES: tuple[str, ...] = (
    "Histone", "TFs and others", "RNA polymerase", "Input control",
    "ATAC-Seq", "DNase-seq", "Bisulfite-Seq", "Annotation tracks",
    "Unclassified",
)


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"chipatlas_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


# ─── Low-level HTTP ─────────────────────────────────────────────────────────

def _pace() -> None:
    """Be polite: 1 req/s default ceiling, env-tunable."""
    global _last_request_time
    gap = float(os.environ.get(
        "CHIPATLAS_MIN_SECONDS", _MIN_SECONDS_BETWEEN_REQUESTS))
    elapsed = time.time() - _last_request_time
    if elapsed < gap:
        time.sleep(gap - elapsed)
    _last_request_time = time.time()


def _request(url: str, *, method: str = "GET", body: bytes | None = None,
              headers: dict[str, str] | None = None,
              accept: str = "application/json,*/*",
              timeout: int = 120) -> tuple[int, bytes, str]:
    _pace()
    hdrs = {"User-Agent": USER_AGENT, "Accept": accept}
    if headers:
        hdrs.update(headers)
    if body is not None and "Content-Type" not in hdrs:
        hdrs["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, headers=hdrs, method=method)
    logging.info("%s %s", method, url)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            content = r.read()
            ct = r.headers.get("Content-Type", "")
            logging.info("  HTTP %d  %d bytes", r.status, len(content))
            return r.status, content, ct
    except urllib.error.HTTPError as e:
        content = e.read()
        logging.warning("  HTTPError %s  %d bytes", e.code, len(content))
        return e.code, content, e.headers.get("Content-Type", "")
    except urllib.error.URLError as e:
        logging.error("  URLError %s", e.reason)
        return 0, str(e.reason).encode(), "text/plain"


def _get_json(url: str) -> Any:
    status, content, ct = _request(url)
    if status < 200 or status >= 300:
        raise SystemExit(f"GET {url}: HTTP {status}")
    if "json" not in ct:
        # Some chip-atlas endpoints return text/plain JSON.
        try:
            return json.loads(content)
        except Exception:
            raise SystemExit(f"GET {url}: expected JSON, got {ct!r}")
    return json.loads(content)


def _http_head(url: str) -> tuple[int, int]:
    """HEAD: returns (status, content_length)."""
    _pace()
    req = urllib.request.Request(url, method="HEAD",
                                    headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, int(r.headers.get("Content-Length", 0) or 0)
    except urllib.error.HTTPError as e:
        return e.code, 0
    except urllib.error.URLError:
        return 0, 0


def _stream_download(url: str, dest: Path, *,
                       max_bytes: int | None = None) -> int:
    """Stream a URL to disk. Returns bytes written. ``max_bytes`` caps
    the download — useful for sanity-only previews of the giant
    ``allPeaks_light.*.bed.gz`` files (default disabled)."""
    _pace()
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    written = 0
    with urllib.request.urlopen(req, timeout=600) as r, dest.open("wb") as fh:
        while True:
            chunk = r.read(1 << 20)  # 1 MiB
            if not chunk:
                break
            fh.write(chunk)
            written += len(chunk)
            if max_bytes is not None and written >= max_bytes:
                logging.info("Stopping at --max-bytes=%d", max_bytes)
                break
    return written


# ─── Browse / list endpoints ────────────────────────────────────────────────

def cmd_list_genomes(args: argparse.Namespace) -> int:
    setup_logging()
    data = _get_json(f"{CHIPATLAS_API}/data/list_of_genome.json")
    print(f"Genomes ({len(data)}):")
    for g in data:
        print(f"  {g}")
    return 0


def cmd_list_qvalues(args: argparse.Namespace) -> int:
    setup_logging()
    data = _get_json(f"{CHIPATLAS_API}/data/qval_range.json")
    print("Allowed -log10(q) threshold suffixes:")
    for q in data:
        print(f"  {q}  → matches eachData/bb{q}/ + bed{q}/")
    print("\nLower index = looser cutoff (more peaks). 05 / 10 / 20 / 50 ≈ "
            "q < 1e-5 / 1e-10 / 1e-20 / 1e-50.")
    return 0


def cmd_list_experiment_types(args: argparse.Namespace) -> int:
    setup_logging()
    url = (f"{CHIPATLAS_API}/data/experiment_types"
            f"?genome={urllib.parse.quote(args.genome)}"
            f"&clClass={urllib.parse.quote(args.cell_class or 'All cell types')}")
    data = _get_json(url)
    print(f"Experiment types for genome={args.genome!r} "
            f"cellClass={args.cell_class or 'All cell types'!r}:")
    for d in data if isinstance(data, list) else []:
        label = d.get("label") or d.get("id")
        count = d.get("count")
        print(f"  {label:30s}  n={count}")
    return 0


def cmd_list_antigens(args: argparse.Namespace) -> int:
    setup_logging()
    url = (f"{CHIPATLAS_API}/data/chip_antigen"
            f"?genome={urllib.parse.quote(args.genome)}"
            f"&agClass={urllib.parse.quote(args.ag_class)}"
            f"&clClass={urllib.parse.quote(args.cell_class or 'All cell types')}")
    data = _get_json(url)
    items = [d for d in (data if isinstance(data, list) else [])
              if d.get("id") not in ("-", "All")]
    print(f"Antigens for genome={args.genome!r} agClass={args.ag_class!r} "
            f"cellClass={args.cell_class or 'All cell types'!r}: {len(items)}")
    items.sort(key=lambda d: -(d.get("count") or 0))
    for d in items[:args.limit]:
        label = d.get("label") or d.get("id")
        count = d.get("count")
        print(f"  {label:30s}  n={count}")
    if len(items) > args.limit:
        print(f"  ... +{len(items) - args.limit} more (use --limit to expand)")
    return 0


def cmd_list_cell_types(args: argparse.Namespace) -> int:
    setup_logging()
    url = (f"{CHIPATLAS_API}/data/sample_types"
            f"?genome={urllib.parse.quote(args.genome)}"
            f"&agClass={urllib.parse.quote(args.ag_class or 'All experiment types')}")
    data = _get_json(url)
    print(f"Cell-type classes for genome={args.genome!r} "
            f"agClass={args.ag_class or 'All experiment types'!r}:")
    for d in data if isinstance(data, list) else []:
        if d.get("id") in ("-", "All"):
            continue
        label = d.get("label") or d.get("id")
        count = d.get("count")
        print(f"  {label:36s}  n={count}")
    return 0


# ─── Search / metadata ──────────────────────────────────────────────────────

def cmd_search(args: argparse.Namespace) -> int:
    setup_logging()
    params: list[tuple[str, str]] = [("q", args.query)]
    if args.genome:
        params.append(("genome", args.genome))
    params.append(("limit", str(args.limit)))
    qs = urllib.parse.urlencode(params)
    data = _get_json(f"{CHIPATLAS_API}/data/search?{qs}")
    if not isinstance(data, list):
        data = data.get("results", []) if isinstance(data, dict) else []
    print(f"Query {args.query!r}  hits: {len(data)}")
    for d in data[:args.limit]:
        srx = d.get("id") or d.get("SRX") or d.get("experiment")
        title = (d.get("title") or d.get("label") or "")[:80]
        genome = d.get("genome", "")
        antigen = d.get("antigen") or d.get("agSubClass") or ""
        cell = d.get("cell") or d.get("clSubClass") or ""
        print(f"  [{genome:6s} {antigen:14s} {cell:30s}] {srx}  {title}")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_search_{safe_label(args.query)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "search.json"
    out_path.write_text(json.dumps(data, indent=2))
    print(f"Saved JSON:  {out_path}")
    return 0


def cmd_get_experiment(args: argparse.Namespace) -> int:
    setup_logging()
    srx = args.experiment_id
    url = f"{CHIPATLAS_API}/data/exp_metadata.json?expid={urllib.parse.quote(srx)}"
    data = _get_json(url)
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_exp_{safe_label(srx)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "metadata.json"
    out_path.write_text(json.dumps(data, indent=2))
    print(f"Experiment:  {srx}")
    rec = data[0] if isinstance(data, list) and data else (data if isinstance(data, dict) else {})
    for k in ("genome", "agClass", "antigen", "agSubClass",
                "clClass", "clSubClass", "title"):
        v = rec.get(k)
        if v is not None:
            print(f"  {k:14s} {v}")
    print(f"Saved:       {out_path}")
    return 0


# ─── Per-experiment file URLs ───────────────────────────────────────────────

def _eachdata_url(genome: str, kind: str, srx: str,
                    qval: str | None = None) -> str:
    """Compose the static archive URL for a per-experiment file.

    ``kind`` ∈ {bb, bw, bb05, bb10, bb20, bed05, bed10, bed20}.
    For ``bb`` / ``bw`` (qval-agnostic), the file is ``{SRX}.{bw|bb}``.
    For ``bb{NN}`` / ``bed{NN}`` the file is ``{SRX}.{NN}.{ext}``.
    """
    if kind in ("bw", "bb"):
        return f"{CHIPATLAS_DATA}/{genome}/eachData/{kind}/{srx}.{kind}"
    # bb05 / bb10 / bb20 / bed05 / bed10 / bed20
    if kind.startswith("bb"):
        return f"{CHIPATLAS_DATA}/{genome}/eachData/{kind}/{srx}.{kind[-2:]}.bb"
    if kind.startswith("bed"):
        return f"{CHIPATLAS_DATA}/{genome}/eachData/{kind}/{srx}.{kind[-2:]}.bed"
    raise SystemExit(f"Unknown eachData kind: {kind}")


def cmd_download_experiment(args: argparse.Namespace) -> int:
    setup_logging()
    kinds = [k.strip() for k in args.kinds.split(",") if k.strip()]
    if not kinds:
        raise SystemExit("--kinds must list at least one of bb,bw,bb05,bb10,bb20,bed05,bed10,bed20")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_dl_{safe_label(args.experiment_id)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for k in kinds:
        url = _eachdata_url(args.genome, k, args.experiment_id)
        if args.urls_only:
            print(f"  {k:8s} {url}")
            rows.append({"kind": k, "url": url, "status": "url-only"})
            continue
        ext = url.rsplit(".", 1)[-1]
        dest = out_dir / f"{args.experiment_id}.{k}.{ext}"
        try:
            n = _stream_download(url, dest)
            rows.append({"kind": k, "url": url, "path": str(dest),
                          "bytes": n, "status": "ok"})
            print(f"  {k:8s} ok  {n:>12,} B  -> {dest}")
        except urllib.error.HTTPError as e:
            rows.append({"kind": k, "url": url, "status": f"HTTP {e.code}"})
            print(f"  {k:8s} HTTP {e.code}  {url}")
    summary = out_dir / "summary.json"
    summary.write_text(json.dumps(rows, indent=2))
    print(f"Summary:     {summary}")
    return 0


# ─── Assembled BED via POST /download ────────────────────────────────────────

def cmd_assemble_bed(args: argparse.Namespace) -> int:
    """POST a (genome × agClass × agSubClass × clClass × clSubClass × qval)
    condition to ``/download`` and return the assembled all-peaks BED
    URL the upstream service builds."""
    setup_logging()
    payload = {"condition": {
        "genome":     args.genome,
        "agClass":    args.ag_class,
        "agSubClass": args.antigen or "",
        "clClass":    args.cell_class or "All cell types",
        "clSubClass": args.cell_subclass or "",
        "qval":       args.qval,
    }}
    body = json.dumps(payload).encode()
    status, content, ct = _request(
        f"{CHIPATLAS_API}/download",
        method="POST", body=body,
        headers={"Content-Type": "application/json"},
        accept="application/json")
    if status < 200 or status >= 300:
        raise SystemExit(f"POST /download failed: HTTP {status}: "
                          f"{content[:300]!r}")
    data = json.loads(content) if "json" in ct else {}
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_assemble"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "assemble.json"
    out_path.write_text(json.dumps({"request": payload,
                                       "response": data}, indent=2))
    url = data.get("url") if isinstance(data, dict) else None
    print(f"Genome:      {args.genome}")
    print(f"Condition:   {args.ag_class}/{args.antigen or '*'} × "
            f"{args.cell_class or 'All'}/{args.cell_subclass or '*'} @ q={args.qval}")
    if url:
        print(f"Assembled URL: {url}")
        if args.fetch:
            dest = out_dir / "assembled.bed"
            n = _stream_download(url, dest,
                                   max_bytes=args.max_bytes or None)
            print(f"Downloaded:   {n:,} bytes -> {dest}")
    else:
        print("No assembled BED for this combination — "
                "(probably no experiments in this slice).")
    print(f"Saved:       {out_path}")
    return 0


# ─── Bulk all-peaks file ────────────────────────────────────────────────────

def _all_peaks_url(genome: str, qval: str) -> str:
    return (f"{CHIPATLAS_DATA}/{genome}/allPeaks_light/"
            f"allPeaks_light.{genome}.{qval}.bed.gz")


def cmd_download_all_peaks(args: argparse.Namespace) -> int:
    """Download (or just probe) the bulk ``allPeaks_light.{genome}.{qval}
    .bed.gz``. These files are *huge* (GB to tens of GB) — use
    ``--head-only`` for size checks, or ``--max-bytes`` to take a
    sanity-sized prefix."""
    setup_logging()
    url = _all_peaks_url(args.genome, args.qval)
    if args.head_only:
        status, size = _http_head(url)
        print(f"HEAD {url}\n  HTTP {status}  Content-Length: {size:,} bytes "
                f"(~{size/1e9:.1f} GB)")
        return 0
    out_dir = SOURCES_DIR / args.genome / "allPeaks_light"
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"allPeaks_light.{args.genome}.{args.qval}.bed.gz"
    n = _stream_download(url, dest, max_bytes=args.max_bytes or None)
    print(f"Wrote {n:,} bytes -> {dest}")
    return 0


# ─── Target-Genes discovery + per-antigen pull ──────────────────────────────

def cmd_target_genes(args: argparse.Namespace) -> int:
    """Discover which antigens have a pre-computed Target-Genes table
    for the requested genome (``--list``), or fetch one and dump a
    preview (``--antigen <X> --distance <bp>``)."""
    setup_logging()
    if args.list:
        data = _get_json(f"{CHIPATLAS_API}/data/target_genes_analysis.json")
        bucket = data.get(args.genome) if isinstance(data, dict) else None
        if not bucket:
            print(f"No Target-Genes catalogue for genome={args.genome!r}.")
            return 1
        print(f"Antigens with Target-Genes tables for {args.genome!r}: {len(bucket)}")
        for a in bucket[:args.limit]:
            print(f"  {a}")
        if len(bucket) > args.limit:
            print(f"  ... +{len(bucket) - args.limit} more (use --limit)")
        return 0
    # Fetch one antigen's table
    if not args.antigen:
        raise SystemExit("Provide --antigen <name> --distance <bp>, "
                          "or use --list to discover available antigens.")
    distance = str(args.distance)
    # Path convention from the upstream archive layout
    url = (f"{CHIPATLAS_DATA}/{args.genome}/target/"
            f"{urllib.parse.quote(args.antigen)}.{distance}.tsv")
    status, content, ct = _request(url, accept="text/tab-separated-values,*/*")
    if status == 404:
        # Alternate naming with cellClass prefix
        alt = (f"{CHIPATLAS_DATA}/{args.genome}/target/"
                f"{urllib.parse.quote(args.antigen)}.{distance}.tab")
        status, content, ct = _request(alt, accept="text/tab-separated-values,*/*")
        url = alt
    if status < 200 or status >= 300:
        raise SystemExit(f"Target-Genes file not found: {url}  (HTTP {status})")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_targetgenes_{safe_label(args.antigen)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{safe_label(args.antigen)}.{distance}.tsv"
    out_path.write_bytes(content)
    lines = content.decode(errors="replace").splitlines()
    n_rows = max(0, len(lines) - 1)
    header = lines[0] if lines else ""
    cols = header.split("\t")
    print(f"Target-Genes table for antigen={args.antigen!r} at "
            f"±{distance}bp (genome={args.genome}):")
    print(f"  Rows:        {n_rows:,}")
    print(f"  Columns:     {len(cols)}")
    if cols:
        print(f"  First cols:  {', '.join(cols[:6])}"
                + ("  ..." if len(cols) > 6 else ""))
    print(f"  Saved TSV:   {out_path}")
    return 0


# ─── WABI Enrichment / Diff job submission ──────────────────────────────────

def cmd_submit_enrichment(args: argparse.Namespace) -> int:
    """Submit a job to the ChIP-Atlas Enrichment / Diff analysis service
    (``POST /wabi_chipatlas`` on the API host, which proxies to the
    NIG/DDBJ WABI queue). Returns the job id and the polling URL.

    The skill supports the two most common modes:

    * ``--mode genes`` — over-representation of TFs at the regulatory
      regions of a query gene list, against a background gene list.
    * ``--mode regions`` — over-representation of TFs at the query BED
      regions, against a control BED.
    """
    setup_logging()
    payload = {"command": args.mode,
                "title": args.label or f"igvfagent_{int(time.time())}",
                "genome": args.genome,
                "antigenClass": args.ag_class or "TFs and others",
                "cellClass": args.cell_class or "All cell types",
                "threshold": args.qval,
                "userDataKey": "0", "dryRun": "false",
                "permTime": "1",
                "distance": str(args.distance) if args.mode == "genes" else "5000",
                # Inputs: file paths uploaded out-of-band, OR newline-joined symbols/BED
                "userInput": Path(args.query).read_text() if Path(args.query).is_file()
                                else args.query,
                "userBg":   (Path(args.background).read_text()
                              if args.background and Path(args.background).is_file()
                              else (args.background or ""))}
    body = json.dumps(payload).encode()
    status, content, ct = _request(
        f"{CHIPATLAS_API}/wabi_chipatlas",
        method="POST", body=body,
        headers={"Content-Type": "application/json"},
        accept="application/json")
    if status < 200 or status >= 300:
        raise SystemExit(f"WABI submit failed: HTTP {status}: "
                          f"{content[:300]!r}")
    try:
        data = json.loads(content)
    except Exception:
        data = {"raw": content.decode(errors="replace")}
    job_id = (data.get("id") or data.get("requestId")
                or data.get("jobId") or "")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_wabi_submit_{safe_label(args.label or 'job')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "submit.json"
    out_path.write_text(json.dumps({"request": payload,
                                       "response": data}, indent=2))
    print(f"WABI mode:   {args.mode}")
    print(f"Job id:      {job_id or '(see submit.json)'}")
    if job_id:
        print(f"Poll URL:    {CHIPATLAS_WABI}/{job_id}?info=result&format=log")
    print(f"Saved:       {out_path}")
    return 0


def cmd_poll_enrichment(args: argparse.Namespace) -> int:
    """Poll a WABI job for its result. Repeats every --interval seconds
    until either the job is done or --timeout is hit."""
    setup_logging()
    url = (f"{CHIPATLAS_WABI}/{urllib.parse.quote(args.job_id)}"
            f"?info=result&format=log")
    deadline = time.time() + args.timeout
    last_status = ""
    while time.time() < deadline:
        status, content, ct = _request(url, accept="text/plain,application/json,*/*")
        body = content.decode(errors="replace")
        first_line = body.splitlines()[0][:160] if body else ""
        if first_line != last_status:
            logging.info("WABI job %s: HTTP %d  '%s'",
                          args.job_id, status, first_line)
            last_status = first_line
        if status == 200 and ("Status: Finished" in body or
                                 "Finish" in first_line or
                                 "Done" in first_line):
            break
        if status >= 400:
            break
        time.sleep(args.interval)
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_wabi_poll_{safe_label(args.job_id)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "result.log"
    out_path.write_bytes(content)
    print(f"Final HTTP:  {status}")
    print(f"Saved log:   {out_path}")
    return 0


# ─── Showcase ───────────────────────────────────────────────────────────────

def cmd_showcase(args: argparse.Namespace) -> int:
    """End-to-end demo over the CTCF / H3K4me3 axis: list the histone
    antigens for hESC, get the per-experiment SRX bigBed URL for one
    canonical CTCF experiment, and (head-only) probe the bulk
    allPeaks_light file.

    No giant downloads — everything is HEAD or sub-MB JSON / TSV."""
    setup_logging()
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_showcase"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {"genome": args.genome}
    # 1. Genomes
    genomes = _get_json(f"{CHIPATLAS_API}/data/list_of_genome.json")
    summary["genomes_available"] = genomes
    print(f"Genomes available: {genomes}")
    # 2. Histone antigens for hESC / pluripotent stem cells
    url = (f"{CHIPATLAS_API}/data/chip_antigen?genome={urllib.parse.quote(args.genome)}"
            f"&agClass=Histone"
            f"&clClass={urllib.parse.quote(args.cell_class)}")
    antigens = _get_json(url)
    histones = [d for d in antigens if d.get("id") not in ("-", "All")]
    histones.sort(key=lambda d: -(d.get("count") or 0))
    summary["top_histones"] = [{"id": d["id"], "count": d["count"]}
                                  for d in histones[:8]]
    print(f"\nTop histones in {args.genome!r} × {args.cell_class!r}:")
    for d in histones[:8]:
        print(f"  {d['id']:14s} n={d['count']}")
    # 3. allPeaks_light HEAD probe for hg38 q05
    url = _all_peaks_url(args.genome, "05")
    status, size = _http_head(url)
    summary["allPeaks_light"] = {"url": url, "status": status, "bytes": size}
    print(f"\nBulk allPeaks_light HEAD: HTTP {status}  {size:,} bytes "
            f"({size/1e9:.1f} GB)")
    # 4. Single SRX per-experiment file head probe (canonical CTCF SRX)
    canonical_srx = args.canonical_srx
    bb_url = _eachdata_url(args.genome, "bb05", canonical_srx)
    bb_status, bb_size = _http_head(bb_url)
    summary["canonical_srx"] = {"srx": canonical_srx, "bb_url": bb_url,
                                  "status": bb_status, "bytes": bb_size}
    print(f"\nPer-experiment BigBed (SRX={canonical_srx}, q05): HTTP "
            f"{bb_status}  {bb_size:,} bytes")
    # 5. List antigens with Target-Genes tables for this genome
    tg = _get_json(f"{CHIPATLAS_API}/data/target_genes_analysis.json")
    n_tg = len(tg.get(args.genome) or [])
    summary["target_genes_available_antigens"] = n_tg
    print(f"\nAntigens with pre-computed Target-Genes for {args.genome!r}: {n_tg}")
    # Write summary + narrative
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    rep = out_dir / "showcase_report.md"
    rep.write_text(
        f"# ChIP-Atlas skill showcase\n\nGenerated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- Genome: **{args.genome}**\n"
        f"- Cell-type class probed: **{args.cell_class!r}**\n"
        f"- Top histone marks (by experiment count): "
        + ", ".join(f"`{d['id']}` (n={d['count']})" for d in histones[:5])
        + "\n"
        f"- Bulk all-peaks file size: **{size/1e9:.1f} GB** (`{url}`)\n"
        f"- Canonical per-experiment BigBed file size: {bb_size:,} bytes\n"
        f"- Antigens with pre-computed Target-Genes tables: **{n_tg}**\n\n"
        f"## License posture\n\n"
        f"Apache-2.0 IGVFagent ⊃ MIT upstream MCP reference. Data are\n"
        f"NBDC/DBCLS licensed (typically CC-BY 4.0); we only fetch / link.\n"
        f"Cite: Zou/Ohta/Oki *Nucleic Acids Res.* 2024.\n"
    )
    print(f"\nReport: {rep}")
    print(f"Summary: {summary_path}")
    return 0


# ─── Playbook writer ────────────────────────────────────────────────────────

def cmd_write_playbook(args: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "CHIPATLAS_SKILL.md"
    path.write_text("""# Skill: ChIP-Atlas canonical-query layer

A polite anonymous client for the [ChIP-Atlas](https://chip-atlas.org)
public REST + bulk archive — reprocessed ChIP-seq / ATAC-seq /
DNase-seq / Bisulfite-seq peak calls and metadata for hundreds of
thousands of public SRA experiments — plus the WABI Enrichment / Diff
Analysis job queue at NIG/DDBJ.

Clean-room reimplementation under Apache-2.0; stdlib-only. Modelled on
[inutano/chip-atlas](https://github.com/inutano/chip-atlas) (MIT,
Ohta/Oki/DBCLS 2017–2026) — we re-derive the wire contract; no source
copied. The skill addresses three hosts indirected through
``_endpoints.py``:

```
chipatlas_api   → https://chip-atlas.org             JSON browse / search / POST-download
chipatlas_data  → https://chip-atlas.dbcls.jp/data   bulk static archive
chipatlas_wabi  → https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas  Enrichment / Diff job queue
```

## Commands

```bash
# Browse — anonymous, JSON
igvfagent chipatlas list-genomes
igvfagent chipatlas list-qvalues
igvfagent chipatlas list-experiment-types --genome hg38
igvfagent chipatlas list-antigens --genome hg38 --ag-class Histone \\
                                    --cell-class 'Pluripotent stem cell'
igvfagent chipatlas list-cell-types --genome hg38 --ag-class Histone

# Search + fetch metadata
igvfagent chipatlas search --query 'CTCF K562' --genome hg38 --limit 25
igvfagent chipatlas get-experiment SRX150531

# Per-experiment files (BigBed, BigWig, BED at q-thresholds 05/10/20)
igvfagent chipatlas download-experiment SRX150531 --genome hg38 \\
                                                     --kinds bw,bb05,bb10
igvfagent chipatlas download-experiment SRX150531 --genome hg38 \\
                                                     --kinds bw,bb05 --urls-only

# Assemble an all-peaks BED for a (ag × cellClass × qval) tuple
igvfagent chipatlas assemble-bed --genome hg38 \\
                                   --ag-class 'TFs and others' --antigen CTCF \\
                                   --cell-class 'Blood' --qval 05
# add --fetch to also stream the assembled BED to disk

# Bulk all-peaks (HEAD only by default — these are 1–30 GB files!)
igvfagent chipatlas download-all-peaks --genome hg38 --qval 05 --head-only
igvfagent chipatlas download-all-peaks --genome hg38 --qval 05 \\
                                          --max-bytes 100000000   # 100 MB sample

# Pre-computed Target-Genes: discover which antigens have tables, then fetch
igvfagent chipatlas target-genes --genome hg38 --list
igvfagent chipatlas target-genes --genome hg38 --antigen H3K4me3 --distance 5000

# WABI Enrichment Analysis (over-representation of TF peaks at query regions)
igvfagent chipatlas submit-enrichment --mode genes --genome hg38 \\
                                         --query gene_list.txt --background bg_list.txt \\
                                         --label my_gene_enrichment_v1
igvfagent chipatlas poll-enrichment <JOB_ID> --interval 15 --timeout 1800

# One-command showcase
igvfagent chipatlas showcase --genome hg38 --cell-class 'Pluripotent stem cell'
```

## Conventions

| Concept | Allowed values |
|---|---|
| `genome` | `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` |
| `agClass` | `Histone`, `TFs and others`, `RNA polymerase`, `Input control`, `ATAC-Seq`, `DNase-seq`, `Bisulfite-Seq`, `Annotation tracks` |
| `qval` | `05` / `10` / `20` / `50` — peak-call `-log10(q)` threshold; lower digit → looser peaks |
| `kinds` (per-experiment) | `bw` (BigWig), `bb` (BigBed all-peaks), `bb05`/`bb10`/`bb20` (BigBed per-qval), `bed05`/`bed10`/`bed20` (BED per-qval) |
| Cell-type tree | 2-level: `clClass` → `clSubClass`, discoverable via `list-cell-types` |

## What this skill adds over IGVFagent's ENCODE/cCRE/enhancer skills

| Capability | Before | After |
|---|---|---|
| Public SRA/GEO ChIP-seq beyond ENCODE | ❌ | ✓ (hundreds of thousands of SRX) |
| Per-experiment BigBed at fixed q-thresholds | partial (ENCODE only) | ✓ |
| (ag × cellClass × qval) assembled BED | ❌ | ✓ |
| Bulk all-peaks_light archive | ❌ | ✓ (HEAD-probe by default) |
| Pre-computed Target-Genes tables | ❌ | ✓ |
| WABI Enrichment / Diff jobs | ❌ | ✓ |
| Cell-type-aware antigen browser | partial | ✓ (counts per slice) |

## Politeness defaults

The client paces itself at **1 request per second** by default to be
gentle with a small group's servers. Override via `CHIPATLAS_MIN_SECONDS`
(env var, accepts fractional seconds, set `0` only for genuinely-needed
bulk pulls).

## License posture

* **Code**: Apache-2.0 IGVFagent ⊃ MIT upstream (`inutano/chip-atlas`).
* **Data**: NBDC/DBCLS LSDB Archive license (typically CC-BY 4.0). We
  only **fetch** / **link**; never redistribute. Read
  `https://dbarchive.biosciencedbc.jp/en/chip-atlas/lic.html` before
  re-publishing downloaded BED/BigBed bytes.
* **Citation**: Zou A, Ohta T, Oki S. *Nucleic Acids Res.* 2024;
  doi:[10.1093/nar/gkae358](https://doi.org/10.1093/nar/gkae358).
  Oki S et al. *EMBO Rep.* 2018;
  doi:[10.15252/embr.201846255](https://doi.org/10.15252/embr.201846255).

## Pairs well with

- `igvfagent enhancer ...` — once a ChIP-Atlas BED slice is downloaded,
  feed it to IGVFagent's ABC enhancer–gene linkage skill.
- `igvfagent ccre ...` — annotate ChIP-Atlas peaks against ENCODE SCREEN
  cCREs to harmonise with the IGVF-Catalog convention.
- `igvfagent enrich ora` — independent statistical check of any
  Target-Genes table.
- `igvfagent catalog find-associations` — cross-reference TF identity
  → IGVF Catalog edges for downstream regulatory analysis.
""")
    print(f"Wrote: {path}")
    return 0


# ─── argparse plumbing ──────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="chipatlas_skill",
        description="ChIP-Atlas canonical-query skill — clean-room reimpl of "
                     "inutano/chip-atlas (MIT) for IGVFagent.")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("list-genomes", help="Supported reference assemblies.")
    p.set_defaults(func=cmd_list_genomes)

    p = sub.add_parser("list-qvalues",
        help="Allowed -log10(q) peak-call thresholds (05/10/20/50).")
    p.set_defaults(func=cmd_list_qvalues)

    p = sub.add_parser("list-experiment-types",
        help="Browse: agClass counts for a genome × cellClass slice.")
    p.add_argument("--genome", required=True)
    p.add_argument("--cell-class", default="All cell types")
    p.set_defaults(func=cmd_list_experiment_types)

    p = sub.add_parser("list-antigens",
        help="Browse: antigens for a (genome × agClass × cellClass) slice "
              "with per-antigen experiment counts.")
    p.add_argument("--genome", required=True)
    p.add_argument("--ag-class", required=True,
                    help="Histone / 'TFs and others' / ATAC-Seq / DNase-seq / etc.")
    p.add_argument("--cell-class", default="All cell types")
    p.add_argument("--limit", type=int, default=40)
    p.set_defaults(func=cmd_list_antigens)

    p = sub.add_parser("list-cell-types",
        help="Browse: cell-type classes for a (genome × agClass) slice.")
    p.add_argument("--genome", required=True)
    p.add_argument("--ag-class", default="All experiment types")
    p.set_defaults(func=cmd_list_cell_types)

    p = sub.add_parser("search",
        help="Free-text search over ChIP-Atlas experiments.")
    p.add_argument("--query", required=True)
    p.add_argument("--genome", default=None)
    p.add_argument("--limit", type=int, default=25)
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("get-experiment",
        help="Full metadata for one SRX/DRX/ERX accession.")
    p.add_argument("experiment_id")
    p.set_defaults(func=cmd_get_experiment)

    p = sub.add_parser("download-experiment",
        help="Per-experiment file URLs / pulls "
              "(bw / bb / bb05 / bb10 / bb20 / bed05 / bed10 / bed20).")
    p.add_argument("experiment_id")
    p.add_argument("--genome", required=True)
    p.add_argument("--kinds", default="bw,bb05",
                    help="Comma-list of file kinds.")
    p.add_argument("--urls-only", action="store_true",
                    help="Print URLs only — don't download.")
    p.set_defaults(func=cmd_download_experiment)

    p = sub.add_parser("assemble-bed",
        help="POST /download → assembled all-peaks BED for a "
              "(genome × ag × cellClass × qval) tuple.")
    p.add_argument("--genome", required=True)
    p.add_argument("--ag-class", required=True,
                    help="Histone / 'TFs and others' / ATAC-Seq / DNase-seq / ...")
    p.add_argument("--antigen", default=None,
                    help="agSubClass — e.g. CTCF, H3K4me3.")
    p.add_argument("--cell-class", default="All cell types")
    p.add_argument("--cell-subclass", default=None)
    p.add_argument("--qval", default="05",
                    choices=("05", "10", "20", "50"))
    p.add_argument("--fetch", action="store_true",
                    help="Stream the assembled BED to disk too.")
    p.add_argument("--max-bytes", type=int, default=None)
    p.set_defaults(func=cmd_assemble_bed)

    p = sub.add_parser("download-all-peaks",
        help="Bulk allPeaks_light.{genome}.{qval}.bed.gz "
              "(GB-scale — default --head-only).")
    p.add_argument("--genome", required=True)
    p.add_argument("--qval", default="05",
                    choices=("05", "10", "20", "50"))
    p.add_argument("--head-only", action="store_true",
                    help="HEAD-probe only; do not download.")
    p.add_argument("--max-bytes", type=int, default=None,
                    help="Cap the download at this many bytes.")
    p.set_defaults(func=cmd_download_all_peaks)

    p = sub.add_parser("target-genes",
        help="Discover or fetch a pre-computed Target-Genes table.")
    p.add_argument("--genome", required=True)
    p.add_argument("--list", action="store_true",
                    help="List antigens with available tables.")
    p.add_argument("--antigen", default=None)
    p.add_argument("--distance", type=int, default=5000,
                    help="Distance window in bp (default 5000).")
    p.add_argument("--limit", type=int, default=50)
    p.set_defaults(func=cmd_target_genes)

    p = sub.add_parser("submit-enrichment",
        help="Submit an Enrichment Analysis job to the NIG WABI queue.")
    p.add_argument("--mode", required=True,
                    choices=("genes", "regions"))
    p.add_argument("--genome", required=True)
    p.add_argument("--query", required=True,
                    help="Path to a gene-list / BED file, OR the literal text.")
    p.add_argument("--background", default=None,
                    help="Background gene-list / BED file or text.")
    p.add_argument("--ag-class", default="TFs and others")
    p.add_argument("--cell-class", default="All cell types")
    p.add_argument("--qval", default="05",
                    choices=("05", "10", "20", "50"))
    p.add_argument("--distance", type=int, default=5000)
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_submit_enrichment)

    p = sub.add_parser("poll-enrichment",
        help="Poll a WABI Enrichment / Diff job until done or timeout.")
    p.add_argument("job_id")
    p.add_argument("--interval", type=int, default=15)
    p.add_argument("--timeout", type=int, default=1800)
    p.set_defaults(func=cmd_poll_enrichment)

    p = sub.add_parser("showcase",
        help="End-to-end demo — no big downloads, all HEAD/JSON.")
    p.add_argument("--genome", default="hg38")
    p.add_argument("--cell-class", default="Pluripotent stem cell")
    p.add_argument("--canonical-srx", default="SRX150531",
                    help="A small per-experiment SRX used for the BigBed HEAD probe.")
    p.set_defaults(func=cmd_showcase)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/CHIPATLAS_SKILL.md.")
    p.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
