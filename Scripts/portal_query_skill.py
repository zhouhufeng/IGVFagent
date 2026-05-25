"""IGVF Portal canonical-query skill.

A faceted, filter-aware, paginate-safe client for the IGVF Portal that
exposes the same query patterns the IGVF Data Coordinating Center (DACC)
documents as canonical in their MCP server reference implementation
[IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp)
(MIT).

Clean-room reimplementation: we re-derive the HTTP wire contract from the
public portal documentation + the upstream README. No source code is
copied — only the URL paths, parameter names, and filter syntax that
are necessary facts of the portal's REST surface. Where IGVFagent
already had a thinner ``/search?type=...`` helper, this module promotes
it to a full canonical-pattern surface:

* HTTP Basic auth via `IGVF_ACCESS_KEY` + `IGVF_SECRET_ACCESS_KEY`
  (preferred), with the legacy cookie path still honored.
* Faceted search (`limit=0` → `facets[]`).
* Tabular `/report.tsv` export.
* `/batch-download/` for FileSet types.
* JSON schema introspection (`/profiles/<Type>.json`).
* Endpoint-param introspection — return both the Python collection_param
  (snake_case, e.g. `file_set_id`) AND the dotted search_field
  (`file_set.@id`) for every documented filter, so an agent can pivot
  between tools without re-learning fields.
* Field-filter DSL with negation (`field!=value`), range ops
  (`gte:` / `lte:` / `gt:` / `lt:`), and list values.

Commands
--------
    portal search              Faceted search with field_filters DSL.
    portal get                 Get a single item by @id / accession / UUID.
    portal schema              JSON schema for an ItemType.
    portal list-types          Enumerate canonical IGVF ItemTypes.
    portal endpoint-params     collection_param ↔ search_field map for a collection.
    portal facets              Facets-only (`limit=0`) call.
    portal report              `/report.tsv` export.
    portal batch-download      `/batch-download/` for FileSet types.
    portal write-playbook      Write Docs/Skills/PORTAL_QUERY_SKILL.md.

License posture
---------------
Apache-2.0. Clean-room reimpl of IGVF-DACC/igvf-portal-mcp (MIT). Uses
only the stdlib (`urllib`); does NOT pin against the official
`igvf-client` PyPI package — IGVFagent keeps its zero-heavy-dependency
posture.
"""

from __future__ import annotations

import argparse
import base64
import csv
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
REPORT_DIR = DOCS_DIR / "PortalQuery"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
CACHE_DIR = DATA_DIR / "Cache" / "Portal"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

PORTAL_BASE = _resolve_endpoint("portal", "IGVF_PORTAL_BASE")
PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
USER_AGENT = "IGVFagent-portal-query/0.1"

# All HTTP queries target the API subdomain (`api.data.igvf.org` by
# default). The data subdomain is browser-UI traffic and requires a
# session cookie; the API subdomain serves anonymous JSON for
# released items and accepts HTTP Basic auth for restricted items.
_API = PORTAL_API_BASE


# ─── Canonical IGVF item types ──────────────────────────────────────────────
#
# Sourced from the upstream README + the live `/profiles/` index. CamelCase
# is the form expected by ``/search/?type=<ItemType>``. The snake_case
# equivalents live under ``/<snake-form>/`` (e.g. `MeasurementSet` ↔
# `/measurement-sets/`). The mapping is purely typographic — derived
# programmatically by lower-casing + hyphenating + pluralising.

CANONICAL_ITEM_TYPES: list[str] = [
    # File-like
    "AlignmentFile", "ConfigurationFile", "GenomeBrowserAnnotationFile",
    "ImageFile", "IndexFile", "MatrixFile", "ModelFile", "ReferenceFile",
    "SequenceFile", "SignalFile", "TabularFile",
    # FileSet hierarchies
    "AnalysisSet", "AuxiliarySet", "ConstructLibrarySet", "CuratedSet",
    "MeasurementSet", "ModelSet", "PredictionSet", "Workflow",
    # Samples / biospecimens
    "InVitroSystem", "PrimaryCell", "Tissue", "TechnicalSample",
    "WholeOrganism", "Biomarker", "MultiplexedSample",
    # Donor / metadata
    "HumanDonor", "RodentDonor", "Phenotype", "PhenotypeTerm",
    "PhenotypicFeature", "Treatment", "Modification",
    # Constructs / reagents
    "Construct", "CrisprModification", "DegronModification",
    "OpenReadingFrame", "Gene", "Document", "OntologyTerm",
    "AssayTerm", "SampleTerm", "PlatformTerm", "TaxonTerm", "PhenotypeTerm",
    # Lab / org
    "Award", "Lab", "User", "Source", "Software", "SoftwareVersion",
    "Page", "Image", "Publication", "Suggestion",
    # Workflow / analysis
    "Pipeline", "AnalysisStep", "AnalysisStepVersion", "QualityMetric",
    "AnalysisStepRun", "RawSequenceFile",
]


# ─── Auth ───────────────────────────────────────────────────────────────────

def portal_auth() -> tuple[str, str] | None:
    """DACC-blessed HTTP Basic credentials from environment variables.

    Read ``IGVF_ACCESS_KEY`` and ``IGVF_SECRET_ACCESS_KEY``. Returns None
    if either is missing (caller falls back to anonymous / cookie auth).
    """
    user = os.environ.get("IGVF_ACCESS_KEY")
    secret = os.environ.get("IGVF_SECRET_ACCESS_KEY")
    if user and secret:
        return (user, secret)
    return None


def auth_header() -> dict[str, str]:
    """Build the headers dict with either Basic auth or a legacy cookie.

    Order of precedence: HTTP Basic (`IGVF_ACCESS_KEY` +
    `IGVF_SECRET_ACCESS_KEY`) wins; falls back to `IGVF_PORTAL_COOKIE`
    if set; otherwise anonymous.
    """
    headers = {"User-Agent": USER_AGENT,
                "Accept": "application/json,text/tab-separated-values,*/*"}
    auth = portal_auth()
    if auth:
        token = base64.b64encode(f"{auth[0]}:{auth[1]}".encode()).decode()
        headers["Authorization"] = f"Basic {token}"
        return headers
    cookie = os.environ.get("IGVF_PORTAL_COOKIE")
    if cookie:
        headers["Cookie"] = cookie
    return headers


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"portal_query_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


# ─── field_filters DSL ──────────────────────────────────────────────────────

def parse_field_filters(spec: str | None) -> list[tuple[str, str]]:
    """Parse a comma-separated field-filter DSL into ordered (key, value)
    pairs suitable for ``urllib.parse.urlencode``.

    Syntax (mirrors the DACC canonical search syntax):

        ``lab.@id=/labs/j-michael-cherry``           equality on dotted embed
        ``file_set.@id!=/analysis-sets/IGVFDS001``   negation
        ``file_size=gte:1000000``                    range
        ``file_format=bam,bed,bigwig``               list (repeated param)
        ``status=released,archived``                 list

    Multiple clauses are separated by ``;`` (NOT ``,`` because commas
    are reserved for list values inside a clause).

        ``lab.@id=/labs/x;file_format=bam,bed``      → two clauses
    """
    if not spec:
        return []
    out: list[tuple[str, str]] = []
    for clause in spec.split(";"):
        clause = clause.strip()
        if not clause:
            continue
        # Detect negation vs equality
        if "!=" in clause:
            key, raw = clause.split("!=", 1)
            negate = True
        elif "=" in clause:
            key, raw = clause.split("=", 1)
            negate = False
        else:
            raise SystemExit(f"Bad field_filter clause: {clause!r} "
                              "(expected 'field=value' or 'field!=value').")
        key = key.strip()
        # Negation is encoded by appending '!' to the field name
        if negate:
            key = key + "!"
        # List values → repeated params
        values = [v.strip() for v in raw.split(",") if v.strip()]
        if not values:
            continue
        # Range ops pass through unchanged (server-side understands them)
        for v in values:
            out.append((key, v))
    return out


def build_search_url(
    *,
    item_types: Iterable[str] | None = None,
    query: str | None = None,
    field_filters: list[tuple[str, str]] | None = None,
    limit: int | str = 25,
    sort: list[str] | None = None,
    frame: str | None = None,
    format_json: bool = True,
    search_path: str = "/search/",
) -> str:
    """Build a fully-encoded ``/search/?...`` URL with the DACC
    parameter conventions. ``limit`` may be an int or the string
    ``"all"`` (which the portal accepts to disable pagination)."""
    params: list[tuple[str, str]] = []
    for t in (item_types or []):
        params.append(("type", t))
    if query:
        params.append(("searchTerm", query))
    if field_filters:
        params.extend(field_filters)
    if sort:
        for s in sort:
            params.append(("sort", s))
    if frame:
        params.append(("frame", frame))
    params.append(("limit", str(limit)))
    if format_json:
        params.append(("format", "json"))
    qs = urllib.parse.urlencode(params, doseq=True)
    return f"{_API}{search_path}?{qs}"


# ─── Low-level HTTP ─────────────────────────────────────────────────────────

def _request(url: str, *, accept: str | None = None,
              method: str = "GET", body: bytes | None = None,
              extra_headers: dict[str, str] | None = None) -> tuple[int, bytes, str]:
    headers = auth_header()
    if accept:
        headers["Accept"] = accept
    if extra_headers:
        headers.update(extra_headers)
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    logging.info("%s %s", method, url)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            content = r.read()
            ct = r.headers.get("Content-Type", "")
            logging.info("  HTTP %d  %d bytes  ct=%s", r.status, len(content), ct)
            return r.status, content, ct
    except urllib.error.HTTPError as e:
        content = e.read()
        logging.warning("  HTTPError %s  %d bytes", e.code, len(content))
        return e.code, content, e.headers.get("Content-Type", "")
    except urllib.error.URLError as e:
        msg = f"{e.reason}".encode()
        logging.error("  URLError %s", e.reason)
        return 0, msg, "text/plain"


def _ensure_json(status: int, content: bytes, ct: str, *,
                  context: str) -> Any:
    if status < 200 or status >= 300:
        try:
            err = json.loads(content)
            detail = err.get("description") or err.get("detail") or str(err)
        except Exception:
            detail = content[:300].decode(errors="replace")
        raise SystemExit(f"{context}: HTTP {status}: {detail}")
    if "json" not in ct:
        raise SystemExit(f"{context}: expected JSON, got {ct!r}")
    return json.loads(content)


def _write_response(label: str, content: bytes, *,
                     ext: str, out_dir: Path | None = None) -> Path:
    out_dir = out_dir or (REPORT_DIR / time.strftime("%Y%m%d_%H%M%S"))
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{safe_label(label)}.{ext}"
    p.write_bytes(content)
    return p


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_list_types(args: argparse.Namespace) -> int:
    """List the canonical IGVF ItemTypes (CamelCase) and their
    snake-case collection paths."""
    print(f"{'CamelCase':40s} {'Collection':40s}")
    print("-" * 82)
    for t in sorted(set(CANONICAL_ITEM_TYPES)):
        # Pluralise: capital boundaries → '-', then lowercase, then 's'
        import re
        snake = re.sub(r"(?<!^)(?=[A-Z])", "-", t).lower()
        if not snake.endswith("s"):
            snake += "s"
        print(f"{t:40s} /{snake}/")
    return 0


def cmd_get(args: argparse.Namespace) -> int:
    """Get one item by @id, accession, or UUID."""
    setup_logging()
    rid = args.resource_id
    # Normalise: accept '/files/IGVFFI.../', 'IGVFFI...', UUID, or '@id'
    if rid.startswith("/"):
        path = rid if rid.endswith("/") else rid + "/"
    elif rid.startswith("IGVF") or len(rid) == 36:  # accession or UUID
        # Try as a top-level item id (portal will resolve to canonical @id)
        path = f"/{rid}/"
    else:
        path = rid if rid.endswith("/") else rid + "/"
    url = f"{_API}{path}?format=json"
    status, content, ct = _request(url, accept="application/json")
    doc = _ensure_json(status, content, ct, context="portal get")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_get"
    out_path = _write_response(safe_label(rid), content, ext="json",
                                 out_dir=out_dir)
    print(f"@id:        {doc.get('@id', '')}")
    print(f"@type:      {', '.join(doc.get('@type', []))}")
    print(f"uuid:       {doc.get('uuid', '')}")
    print(f"status:     {doc.get('status', '')}")
    print(f"accession:  {doc.get('accession', '')}")
    print(f"summary:    {doc.get('summary', '')[:200]}")
    print(f"Written:    {out_path}")
    return 0


def cmd_schema(args: argparse.Namespace) -> int:
    """Fetch the JSON schema for an ItemType."""
    setup_logging()
    t = args.item_type
    url = f"{_API}/profiles/{t}.json"
    status, content, ct = _request(url, accept="application/json")
    schema = _ensure_json(status, content, ct, context="portal schema")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_schema_{safe_label(t)}"
    out_path = _write_response(t, content, ext="json", out_dir=out_dir)
    print(f"Title:       {schema.get('title', '')}")
    print(f"Description: {(schema.get('description') or '')[:200]}")
    props = list((schema.get("properties") or {}).keys())
    print(f"Properties:  {len(props)}")
    if props:
        # Print up to 20 in 4 columns
        n = min(20, len(props))
        for i in range(0, n, 4):
            print("  " + "  ".join(f"{p:<22}" for p in props[i:i + 4]))
        if len(props) > 20:
            print(f"  ... +{len(props) - 20} more")
    print(f"Schema JSON: {out_path}")
    return 0


def cmd_search(args: argparse.Namespace) -> int:
    """Faceted search with the field_filters DSL."""
    setup_logging()
    item_types = [t.strip() for t in (args.type or "").split(",") if t.strip()]
    field_filters = parse_field_filters(args.field_filters)
    sort = [s.strip() for s in (args.sort or "").split(",") if s.strip()] or None
    url = build_search_url(
        item_types=item_types or None,
        query=args.query,
        field_filters=field_filters,
        limit=("all" if args.limit == "all" else int(args.limit)),
        sort=sort, frame=args.frame,
    )
    status, content, ct = _request(url, accept="application/json")
    data = _ensure_json(status, content, ct, context="portal search")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_search_{safe_label(args.label or 'q')}"
    out_path = _write_response("search", content, ext="json", out_dir=out_dir)
    hits = data.get("@graph", []) or []
    total = data.get("total", len(hits))
    print(f"Total matches: {total:,}")
    print(f"Returned:      {len(hits):,}")
    print(f"Saved JSON:    {out_path}")
    print()
    # Show top 8 hits — accession + type + summary
    for h in hits[:8]:
        acc = h.get("accession") or h.get("@id", "")
        types = ",".join((h.get("@type") or [])[:2])
        summary = (h.get("summary") or h.get("title") or "")[:80]
        print(f"  [{types}] {acc}  {summary}")
    if len(hits) > 8:
        print(f"  ... +{len(hits) - 8} more in JSON.")
    return 0


def cmd_facets(args: argparse.Namespace) -> int:
    """Facets-only call (`limit=0`)."""
    setup_logging()
    item_types = [t.strip() for t in (args.type or "").split(",") if t.strip()]
    field_filters = parse_field_filters(args.field_filters)
    url = build_search_url(
        item_types=item_types or None,
        query=args.query,
        field_filters=field_filters,
        limit=0,
    )
    status, content, ct = _request(url, accept="application/json")
    data = _ensure_json(status, content, ct, context="portal facets")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_facets_{safe_label(args.label or 'q')}"
    out_path = _write_response("facets", content, ext="json", out_dir=out_dir)
    facets = data.get("facets", []) or []
    total = data.get("total", 0)
    print(f"Total matching items: {total:,}")
    print(f"Facet groups:         {len(facets)}")
    print(f"Saved JSON:           {out_path}")
    print()
    for f in facets[:12]:
        title = f.get("title") or f.get("field")
        terms = f.get("terms", []) or []
        n_terms = len(terms)
        top = terms[:5]
        rows = "; ".join(f"{t.get('key', '?')}={t.get('doc_count', 0)}"
                          for t in top)
        print(f"  {title}  ({n_terms} values)")
        if rows:
            print(f"     {rows}")
    if len(facets) > 12:
        print(f"  ... +{len(facets) - 12} more facet groups in JSON.")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    """Fetch the canonical TSV `/report.tsv` export."""
    setup_logging()
    item_types = [t.strip() for t in (args.type or "").split(",") if t.strip()]
    if not item_types:
        raise SystemExit("portal report: --type is required")
    field_filters = parse_field_filters(args.field_filters)
    url = build_search_url(
        item_types=item_types,
        query=args.query,
        field_filters=field_filters,
        limit=("all" if args.limit == "all" else int(args.limit)),
        format_json=False,
        search_path="/report.tsv",
    )
    # Strip the redundant format=json the helper adds (we set format_json=False)
    status, content, ct = _request(url, accept="text/tab-separated-values")
    if status < 200 or status >= 300:
        raise SystemExit(f"portal report: HTTP {status}: {content[:200]!r}")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_report_{safe_label(args.label or item_types[0])}"
    out_path = _write_response(
        "report", content, ext="tsv", out_dir=out_dir)
    # Quick stats — IGVF /report.tsv has a 1-line preamble
    # (timestamp + query URL) followed by the real header row, so we
    # find the first line with multiple tabs as the true column header.
    text = content.decode(errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    header_idx = 0
    for i, ln in enumerate(lines):
        if ln.count("\t") >= 3:  # heuristic: real header has many columns
            header_idx = i
            break
    header = lines[header_idx] if lines else ""
    n_cols = len(header.split("\t")) if header else 0
    n_rows = max(0, len(lines) - header_idx - 1)
    print(f"Rows (excl. header): {n_rows:,}")
    print(f"Columns:             {n_cols}")
    print(f"Saved TSV:           {out_path}")
    if header:
        cols = header.split("\t")
        preview = cols[:8]
        print(f"First {len(preview)} columns: {', '.join(preview)}"
              + ("  ..." if len(cols) > 8 else ""))
    return 0


def cmd_batch_download(args: argparse.Namespace) -> int:
    """Hit the IGVF Portal's `/batch-download/` endpoint for a FileSet
    type. Returns a list of pre-signed S3 download URLs (the portal
    convention), saved as a manifest TXT."""
    setup_logging()
    item_types = [t.strip() for t in (args.type or "").split(",") if t.strip()]
    if not item_types:
        raise SystemExit("portal batch-download: --type is required")
    field_filters = parse_field_filters(args.field_filters)
    # Build a /batch-download/?type=...&... URL
    url = build_search_url(
        item_types=item_types,
        query=args.query,
        field_filters=field_filters,
        limit=("all" if args.limit == "all" else int(args.limit)),
        format_json=False,
        search_path="/batch-download/",
    )
    status, content, ct = _request(url, accept="text/plain")
    if status < 200 or status >= 300:
        raise SystemExit(
            f"portal batch-download: HTTP {status}: {content[:200]!r}")
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_batchdl_{safe_label(args.label or item_types[0])}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir / "batch_download.txt"
    manifest.write_bytes(content)
    n_urls = sum(1 for ln in content.decode(errors="replace").splitlines()
                   if ln.strip() and not ln.startswith("#"))
    print(f"Pre-signed URLs: {n_urls:,}")
    print(f"Manifest:        {manifest}")
    if args.fetch:
        # Optionally pull every file referenced in the manifest
        files_dir = out_dir / "files"
        files_dir.mkdir(exist_ok=True)
        n_done = 0
        for ln in content.decode(errors="replace").splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            fname = ln.rsplit("/", 1)[-1].split("?", 1)[0]
            dest = files_dir / fname
            if dest.exists() and dest.stat().st_size > 0:
                n_done += 1
                continue
            try:
                with urllib.request.urlopen(ln, timeout=120) as r:
                    dest.write_bytes(r.read())
                n_done += 1
            except Exception as e:
                logging.warning("download %s failed: %s", fname, e)
        print(f"Downloaded:      {n_done}/{n_urls} to {files_dir}")
    return 0


def cmd_endpoint_params(args: argparse.Namespace) -> int:
    """Introspect filters for a single collection.

    Returns the `collection_param ↔ search_field` map: the snake-cased
    Python-friendly name an agent would pass to a tool, and the dotted
    `field.@id`-style name the search index understands. This is the
    central UX trick from the upstream MCP server — lets an LLM pivot
    between the typed-tool surface and the raw search surface without
    re-learning fields.
    """
    setup_logging()
    coll = args.collection.strip("/")  # accept 'measurement-sets' or '/measurement-sets/'
    # Hit the collection root with limit=0 — facets reveal documented filters
    url = f"{_API}/{coll}/?limit=0&format=json"
    status, content, ct = _request(url, accept="application/json")
    data = _ensure_json(status, content, ct,
                         context="portal endpoint-params")
    facets = data.get("facets", []) or []
    columns = data.get("columns", {}) or {}
    rows: list[dict[str, str]] = []
    # Field map from facets
    for f in facets:
        field = f.get("field") or ""
        if not field:
            continue
        title = f.get("title") or field
        # Synthesise the snake_case agent-facing name
        # `file_set.@id` -> `file_set_id`; `lab.name` -> `lab_name`;
        # plain `status` -> `status`.
        snake = (field.replace(".@id", "_id")
                       .replace(".", "_"))
        rows.append({
            "agent_param": snake,
            "search_field": field,
            "title":        title,
            "n_terms":      str(len(f.get("terms", []) or [])),
        })
    # Also surface any non-facet columns documented by the collection
    for col, meta in columns.items():
        if any(r["search_field"] == col for r in rows):
            continue
        title = (meta or {}).get("title") if isinstance(meta, dict) else None
        rows.append({"agent_param": col.replace(".", "_"),
                      "search_field": col, "title": title or col,
                      "n_terms": ""})
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_params_{safe_label(coll)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "endpoint_params.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh,
            fieldnames=["agent_param", "search_field", "title", "n_terms"],
            delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{'agent_param':<32s} {'search_field':<32s} {'n_terms':>8s}  title")
    print("-" * 100)
    for r in rows[:30]:
        print(f"{r['agent_param']:<32s} {r['search_field']:<32s} "
              f"{r['n_terms']:>8s}  {r['title']}")
    if len(rows) > 30:
        print(f"... +{len(rows) - 30} more in TSV")
    print(f"\nWrote: {out}")
    return 0


def cmd_write_playbook(args: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "PORTAL_QUERY_SKILL.md"
    path.write_text("""# Skill: IGVF Portal canonical-query layer

A faceted, filter-aware, paginate-safe client for the IGVF Portal that
matches the same query patterns the IGVF Data Coordinating Center (DACC)
documents as canonical in [IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp)
(MIT, Cherry Lab / IGVF DACC, 2024).

Clean-room reimplementation under Apache-2.0 — we re-derive the HTTP wire
contract from the public portal documentation and the upstream README.
No source code is copied. The skill is stdlib-only (`urllib`); it does
**not** depend on the official `igvf-client` PyPI package.

## Commands

```bash
# Search with field_filters
igvfagent portal search --type MeasurementSet \\
                          --field-filters "lab.@id=/labs/j-michael-cherry;file_format=bam,bigwig" \\
                          --limit 25

# Get one item by accession / @id / UUID
igvfagent portal get IGVFDS3909HJKS

# JSON schema for an ItemType
igvfagent portal schema MeasurementSet

# Canonical ItemType enum
igvfagent portal list-types

# collection_param ↔ search_field map (the UX trick)
igvfagent portal endpoint-params measurement-sets

# Facets-only call (limit=0)
igvfagent portal facets --type MeasurementSet --query CRISPR

# Tabular TSV /report.tsv
igvfagent portal report --type SequenceFile \\
                          --field-filters "file_format=fastq" \\
                          --limit 500 --label fastq_v1

# /batch-download/ for FileSet types (manifest of pre-signed S3 URLs)
igvfagent portal batch-download --type AnalysisSet \\
                                  --field-filters "status=released" \\
                                  --label analysis_v1
# add --fetch to actually pull every file referenced in the manifest
```

## Authentication

| Env var | Role |
|---|---|
| `IGVF_ACCESS_KEY` + `IGVF_SECRET_ACCESS_KEY` | DACC-blessed HTTP Basic credentials (preferred) |
| `IGVF_PORTAL_COOKIE` | Legacy cookie auth (still honored) |
| _(neither)_ | Anonymous — public-released items only |

Without auth the portal returns 403 on access-restricted items; the
skill still works for browsing released data.

## The field_filters DSL

A small, predictable filter language matching the DACC search syntax:

| Clause | Wire form |
|---|---|
| `lab.@id=/labs/jb-cherry` | `lab.@id=/labs/jb-cherry` |
| `lab.@id!=/labs/x` | `lab.@id%21=/labs/x`  (`!=` negation) |
| `file_size=gte:1000000` | `file_size=gte%3A1000000` |
| `file_format=bam,bigwig` | `file_format=bam&file_format=bigwig` (repeated param) |
| `status=released,archived` | `status=released&status=archived` |

Multiple clauses joined with `;`:

```
lab.@id=/labs/x;file_format=bam,bed;file_size=gte:1000000
```

## What this skill adds over the legacy `data` / `frontpage` skills

| Capability | Before | After |
|---|---|---|
| HTTP Basic auth | ❌ (cookie only) | ✓ |
| Field-filter DSL with negation + ranges | ❌ | ✓ |
| `/report.tsv` export | ❌ | ✓ |
| `/batch-download/` manifest + fetch | ❌ | ✓ |
| Facets-only call | partial | ✓ |
| Endpoint-param introspection | ❌ | ✓ |
| Canonical ItemType enum | ❌ | ✓ |

## License posture

Apache-2.0. Clean-room reimplementation of IGVF-DACC/igvf-portal-mcp
(MIT, Cherry Lab) — we cite the upstream README + the public portal
HTTP contract as factual sources. The upstream MIT licence would
permit verbatim copy with attribution, but we re-implement to keep the
stdlib-only posture (no `igvf-client==121.0.0` pin).
""")
    print(f"Wrote: {path}")
    return 0


# ─── argparse plumbing ──────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="portal_query_skill",
        description="IGVF Portal canonical-query skill — clean-room reimpl "
                     "of IGVF-DACC/igvf-portal-mcp patterns (MIT).")
    sub = parser.add_subparsers(dest="cmd")

    def _add_filter_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("--type", default=None,
                        help="Comma-separated CamelCase ItemTypes (e.g. "
                              "'MeasurementSet,AnalysisSet').")
        p.add_argument("--query", default=None,
                        help="Free-text search.")
        p.add_argument("--field-filters", default=None,
                        help="DSL: 'lab.@id=/labs/x;file_format=bam,bed'.")
        p.add_argument("--label", default=None)

    p = sub.add_parser("search", help="Faceted search with field_filters DSL.")
    _add_filter_args(p)
    p.add_argument("--limit", default="25",
                    help="Integer or 'all' (default 25).")
    p.add_argument("--sort", default=None,
                    help="Comma list of sort keys (prefix - for descending).")
    p.add_argument("--frame", default=None,
                    help="object / embedded / page (portal default).")
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("get", help="Get one item by @id / accession / UUID.")
    p.add_argument("resource_id",
                    help="@id, accession (IGVFFI...), or UUID.")
    p.set_defaults(func=cmd_get)

    p = sub.add_parser("schema", help="JSON schema for an ItemType.")
    p.add_argument("item_type",
                    help="CamelCase ItemType (e.g. MeasurementSet).")
    p.set_defaults(func=cmd_schema)

    p = sub.add_parser("list-types",
        help="Enumerate canonical IGVF ItemTypes + their collection paths.")
    p.set_defaults(func=cmd_list_types)

    p = sub.add_parser("endpoint-params",
        help="collection_param ↔ search_field map for a collection.")
    p.add_argument("collection",
                    help="Snake-case collection (e.g. measurement-sets).")
    p.set_defaults(func=cmd_endpoint_params)

    p = sub.add_parser("facets", help="Facets-only (limit=0) call.")
    _add_filter_args(p)
    p.set_defaults(func=cmd_facets)

    p = sub.add_parser("report", help="`/report.tsv` export.")
    _add_filter_args(p)
    p.add_argument("--limit", default="all",
                    help="Integer or 'all' (default 'all').")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("batch-download",
        help="`/batch-download/` manifest of pre-signed URLs for a FileSet type.")
    _add_filter_args(p)
    p.add_argument("--limit", default="all")
    p.add_argument("--fetch", action="store_true",
                    help="After fetching the manifest, actually download "
                          "every referenced file.")
    p.set_defaults(func=cmd_batch_download)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/PORTAL_QUERY_SKILL.md.")
    p.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
