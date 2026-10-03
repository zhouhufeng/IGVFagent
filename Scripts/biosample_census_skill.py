#!/usr/bin/env python3
"""Systematic census of everything ENCODE and IGVF hold for one biosample.

Asked "summarise the ENCODE and IGVF data for GM12878, with tables and
plots", the hosted agent had no single tool that answered it. It tried
``explain_dataset`` on an ENCODE search URL (a 404 when the filter matched
nothing), then authored a one-off census skill at runtime, which crashed on
the IGVF Portal's ``file_size`` facet -- whose ``terms`` is a stats *dict*, not
the list every other facet returns -- and could not re-register the fixed
version because the duplication guard matched it against itself. This module
is that census, built in, so the same command runs on the deployment and on a
laptop without anything being authored first.

Two things about the counting are deliberate and worth knowing:

* **Structured filters, not free text.** ``searchTerm=GM12878`` on the IGVF
  Portal matches 7,133 measurement sets -- essentially the whole Portal --
  while the structured ``samples.sample_terms.term_name=GM12878`` matches 5.
  ENCODE's free text is tighter but still counts every record whose
  description mentions the line. Every number here comes from the ontology
  term field the portal itself facets on, so "GM12878" means the sample was
  GM12878, not that the word appears somewhere.
* **A 404 with ``total: 0`` is an answer, not an error.** Both portals return
  HTTP 404 for a search that matches nothing. The census records it as zero.

Subcommand::

    igvfagent biosample-census run --biosample GM12878 [--label x]
        [--portal both|encode|igvf] [--status released|all]
        [--max-items 5000] [--top 20] [--no-plots]

Writes ``Docs/BiosampleCensus/<timestamp>_<label>/`` with ``report.md``,
``totals_by_type.csv``, ``facet_counts.csv``, per-portal item tables, SVG
figures (PNG too when matplotlib is installed) and ``census.json`` recording
every request made. Standard library only; matplotlib is optional.
"""
from __future__ import annotations

import argparse
import base64
import collections
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
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _http  # noqa: F401,E402  (IPv4-preferred resolver, import for side effect)
from _endpoints import resolve as _resolve_endpoint  # noqa: E402
from _credentials import portal_credentials as _portal_credentials  # noqa: E402

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "BiosampleCensus"

ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
IGVF_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
USER_AGENT = "IGVFdataAgent/0.1 biosample-census"
TIMEOUT = 90

# ---------------------------------------------------------------------------
# What to count. Each entry: object type, the structured field that pins the
# biosample, the facets worth tabulating, and whether to pull item rows.
# Fields were verified against both portals on 2026-09-21 -- an unknown field
# is not an error on either portal, it is a silent zero, which is exactly the
# failure this module exists to avoid.
# ---------------------------------------------------------------------------

ENCODE_TYPES: "list[dict]" = [
    {"type": "Experiment",
     "field": "biosample_ontology.term_name",
     "facets": ["assay_title", "target.label", "lab.title", "status",
                "assembly", "award.rfa", "date_released",
                "replicates.library.biosample.treatments.treatment_term_name",
                "files.file_type"],
     "items": True,
     "item_fields": ["accession", "assay_title", "assay_term_name",
                     "target.label", "lab.title", "status", "date_released",
                     "assembly", "award.rfa", "biosample_ontology.term_id",
                     "replicates.library.biosample.treatments.treatment_term_name",
                     "description"]},
    {"type": "FunctionalCharacterizationExperiment",
     "field": "biosample_ontology.term_name",
     "facets": ["assay_title", "lab.title", "status", "assembly"],
     "items": True,
     "item_fields": ["accession", "assay_title", "assay_term_name",
                     "lab.title", "status", "date_released", "assembly",
                     "description"]},
    {"type": "Annotation",
     "field": "biosample_ontology.term_name",
     "facets": ["annotation_type", "lab.title", "status", "assembly",
                "encyclopedia_version", "date_released"],
     "items": True,
     "item_fields": ["accession", "annotation_type", "lab.title", "status",
                     "date_released", "assembly", "encyclopedia_version",
                     "description"]},
    {"type": "ReferenceEpigenome", "field": "biosample_ontology.term_name",
     "facets": ["lab.title", "status"], "items": True,
     "item_fields": ["accession", "lab.title", "status", "description"]},
    {"type": "Series", "field": "biosample_ontology.term_name",
     "facets": ["type", "lab.title", "status"], "items": True,
     "item_fields": ["accession", "@type", "lab.title", "status",
                     "description"]},
    {"type": "Biosample", "field": "biosample_ontology.term_name",
     "facets": ["lab.title", "status", "treatments.treatment_term_name",
                "genetic_modifications.category"],
     "items": False},
    {"type": "File", "field": "biosample_ontology.term_name",
     "facets": ["file_format", "output_type", "output_category", "assembly",
                "status", "lab.title", "file_type"],
     "items": False},
]

IGVF_TYPES: "list[dict]" = [
    {"type": "MeasurementSet", "field": "samples.sample_terms.term_name",
     "facets": ["preferred_assay_titles", "assay_slims", "lab.title", "status",
                "files.file_format", "files.content_type",
                "samples.classifications"], "items": True},
    {"type": "AnalysisSet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "AuxiliarySet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "PredictionSet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "ConstructLibrarySet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "ModelSet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "CuratedSet", "field": "samples.sample_terms.term_name",
     "facets": ["file_set_type", "lab.title", "status"], "items": True},
    {"type": "Sample", "field": "sample_terms.term_name",
     "facets": ["classifications", "lab.title", "status",
                "treatments.treatment_term_name", "modifications.modality"],
     "items": True},
    {"type": "File", "field": "file_set.samples.sample_terms.term_name",
     "facets": ["file_format", "content_type", "assembly", "status",
                "lab.title", "file_set.file_set_type",
                "preferred_assay_titles"], "items": True},
]

IGVF_ITEM_FIELDS = ["accession", "@type", "preferred_assay_titles",
                    "assay_titles", "file_set_type", "lab.title", "status",
                    "release_timestamp", "summary", "content_type",
                    "file_format", "assembly", "file_size", "file_set.@id",
                    "classifications", "sample_terms.term_id",
                    "samples.sample_terms.term_id"]

# Reference palette from the dataviz method: one hue per portal, text in ink
# tokens, recessive axes. ENCODE = slot 1 blue, IGVF = slot 2 orange.
COLOR = {"encode": "#2a78d6", "igvf": "#eb6834"}
INK, INK2, AXIS, SURFACE = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"biosample_census_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def headers_for(portal: str) -> "dict[str, str]":
    h = {"Accept": "application/json,*/*", "User-Agent": USER_AGENT}
    if portal == "igvf":
        creds = _portal_credentials()
        if creds:
            tok = base64.b64encode(f"{creds[0]}:{creds[1]}".encode()).decode()
            h["Authorization"] = f"Basic {tok}"
        elif os.environ.get("IGVF_PORTAL_COOKIE"):
            h["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return h


class Portal:
    """Thin JSON client that also keeps an audit trail of every request."""

    def __init__(self, name: str, base: str) -> None:
        self.name, self.base = name, base
        self.log: "list[dict]" = []

    def url(self, params: "list[tuple[str, str]]") -> str:
        return (f"{self.base}/search/?"
                + urllib.parse.urlencode(params, doseq=True,
                                          quote_via=urllib.parse.quote))

    def get(self, params: "list[tuple[str, str]]") -> "tuple[int, Any]":
        url = self.url(params)
        logging.info("GET %s", url)
        req = urllib.request.Request(url, headers=headers_for(self.name))
        status, data = 0, None
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                status, data = r.status, json.loads(r.read())
        except urllib.error.HTTPError as e:
            status = e.code
            try:
                data = json.loads(e.read())
            except Exception:
                data = {"http_error": status}
        except (urllib.error.URLError, OSError, ValueError) as e:
            status, data = 0, {"network_error": str(e)}
        total = data.get("total") if isinstance(data, dict) else None
        self.log.append({"url": url, "status": status, "total": total})
        time.sleep(0.05)
        return status, data


def is_zero_hit(status: int, data: Any) -> bool:
    """Both portals answer a no-match search with 404 + total 0."""
    return status == 404 and isinstance(data, dict) and data.get("total") == 0


def graph_rows(data: Any) -> "list[dict]":
    if isinstance(data, dict) and isinstance(data.get("@graph"), list):
        return [r for r in data["@graph"] if isinstance(r, dict)]
    return []


# ---------------------------------------------------------------------------
# Facet parsing -- the part that crashed at runtime
# ---------------------------------------------------------------------------

def iter_facet_terms(field: str, terms: Any) -> "Iterable[tuple[str, str, int]]":
    """Yield (field, term, count) from a facet's ``terms`` in any portal shape.

    Seen in the wild: a list of ``{key, doc_count}``; booleans and dates that
    carry ``key_as_string`` beside a numeric ``key``; ``hierarchical`` facets
    whose entries nest a ``subfacet`` with its own terms; and ``stats`` facets
    (IGVF ``file_size``) whose ``terms`` is a dict of min/max/avg. The last one
    is what ``for t in facet["terms"]`` iterated as dict keys.
    """
    if isinstance(terms, dict):
        return  # stats block -- summarised separately by facet_stats()
    if not isinstance(terms, list):
        return
    for t in terms:
        if not isinstance(t, dict):
            continue
        key = t.get("key_as_string", t.get("key"))
        if key is None:
            continue
        count = t.get("doc_count")
        if isinstance(count, bool) or not isinstance(count, int):
            try:
                count = int(count)
            except (TypeError, ValueError):
                continue
        yield field, str(key), count
        sub = t.get("subfacet")
        if isinstance(sub, dict):
            yield from iter_facet_terms(str(sub.get("field") or f"{field}.sub"),
                                        sub.get("terms"))


def facet_table(data: Any, wanted: "Optional[Iterable[str]]" = None
                ) -> "dict[str, list[tuple[str, int]]]":
    """{facet field: [(term, count), ...] sorted desc} for the facets asked for.

    ``wanted`` may name a hierarchical facet's parent or its subfacet field;
    both are collected. With ``wanted=None`` every term facet is returned.
    """
    want = set(wanted) if wanted is not None else None
    out: "dict[str, list[tuple[str, int]]]" = collections.defaultdict(list)
    facets = data.get("facets") if isinstance(data, dict) else None
    for f in facets or []:
        if not isinstance(f, dict) or not f.get("field"):
            continue
        for field, term, count in iter_facet_terms(str(f["field"]), f.get("terms")):
            if want is None or field in want or str(f["field"]) in want:
                out[field].append((term, count))
    for k in out:
        out[k].sort(key=lambda kv: (-kv[1], kv[0]))
    return dict(out)


def facet_stats(data: Any) -> "dict[str, dict]":
    """The ``stats``-type facets (e.g. IGVF file_size) as {field: stats}."""
    out = {}
    facets = data.get("facets") if isinstance(data, dict) else None
    for f in facets or []:
        if isinstance(f, dict) and isinstance(f.get("terms"), dict) and f.get("field"):
            out[str(f["field"])] = {k: v for k, v in f["terms"].items()
                                    if isinstance(v, (int, float))}
    return out


# ---------------------------------------------------------------------------
# Term resolution -- "gm12878" must find GM12878, and we want the ontology id
# ---------------------------------------------------------------------------

def _term_rows(portal: Portal, params: "list[tuple[str, str]]") -> "list[dict]":
    status, data = portal.get(params + [("format", "json"), ("limit", "50"),
                                        ("field", "term_name"), ("field", "term_id"),
                                        ("field", "classification")])
    return graph_rows(data) if status == 200 else []


def _casings(text: str) -> "list[str]":
    t = text.strip()
    out = []
    for v in (t, t.upper(), t.lower(), t.title(), t.capitalize()):
        if v and v not in out:
            out.append(v)
    return out


def resolve_term(portal: Portal, biosample: str,
                 term_id_hint: str = "") -> "dict":
    """Canonical term_name + term_id on this portal, or the input verbatim.

    Exact ``term_name=`` lookups first (both portals are case-sensitive, so
    the obvious casings are tried), then the ontology id when another portal
    already resolved it -- ENCODE and IGVF share EFO/CL/UBERON ids, so
    ``EFO:0002784`` finds GM12878 on both -- then ENCODE's typeahead facet.
    IGVF's ``searchTerm`` on SampleTerm is deliberately NOT used: it returns
    the entire vocabulary (2,060 terms) for any query, so it never resolves
    anything.
    """
    otype = "BiosampleType" if portal.name == "encode" else "SampleTerm"
    want = biosample.strip().lower()
    squash = "".join(ch for ch in want if ch.isalnum())

    def pick(rows: "list[dict]") -> "Optional[dict]":
        exact = [r for r in rows if str(r.get("term_name", "")).lower() == want]
        if not exact:
            exact = [r for r in rows
                     if "".join(ch for ch in str(r.get("term_name", "")).lower()
                                if ch.isalnum()) == squash]
        if not exact:
            return None
        ids = sorted({str(x.get("term_id")) for x in exact if x.get("term_id")})
        return {"term_name": str(exact[0]["term_name"]), "term_id": ", ".join(ids),
                "classification": _scalar(exact[0].get("classification")),
                "resolved": True}

    for casing in _casings(biosample):
        hit = pick(_term_rows(portal, [("type", otype), ("term_name", casing)]))
        if hit:
            return hit
    if term_id_hint:
        for tid in term_id_hint.split(","):
            rows = _term_rows(portal, [("type", otype), ("term_id", tid.strip())])
            if rows and rows[0].get("term_name"):
                r = rows[0]
                return {"term_name": str(r["term_name"]), "term_id": str(r.get("term_id") or tid),
                        "classification": _scalar(r.get("classification")),
                        "resolved": True, "via": f"term_id {tid.strip()}"}
    if portal.name == "encode":
        hit = pick(_term_rows(portal, [("type", otype), ("searchTerm", biosample.strip())]))
        if hit:
            return hit
    else:
        # Case-insensitive match against the sample-term facet of the sets.
        status, data = portal.get([("type", "MeasurementSet"), ("format", "json"), ("limit", "0")])
        for field, term, _count in iter_facet_terms("samples.sample_terms.term_name",
                                                    next((f.get("terms") for f in (data.get("facets") or [])
                                                          if isinstance(f, dict)
                                                          and f.get("field") == "samples.sample_terms.term_name"),
                                                         None) if isinstance(data, dict) else None):
            if term.lower() == want:
                return {"term_name": term, "term_id": "", "classification": "",
                        "resolved": True, "via": "facet"}
    return {"term_name": biosample.strip(), "term_id": "", "classification": "",
            "resolved": False}


# ---------------------------------------------------------------------------
# The census proper
# ---------------------------------------------------------------------------

def census_portal(portal: Portal, spec_list: "list[dict]", term: str, *,
                  status_filter: str, max_items: int) -> "dict":
    totals: "list[dict]" = []
    facet_rows: "list[dict]" = []
    items: "dict[str, list[dict]]" = {}
    stats: "dict[str, dict]" = {}
    for spec in spec_list:
        t = spec["type"]
        base = [("type", t), (spec["field"], term), ("format", "json")]
        # One unfiltered facet request: total + status breakdown + every
        # facet in one round trip, for all statuses.
        st, data = portal.get(base + [("limit", "0")])
        if is_zero_hit(st, data):
            total, ok = 0, True
        elif st == 200 and isinstance(data, dict):
            total, ok = int(data.get("total") or 0), True
        else:
            total, ok = -1, False
        ft = facet_table(data, list(spec["facets"]) + ["status"]) if ok else {}
        by_status = dict(ft.get("status", []))
        totals.append({"portal": portal.name, "type": t, "total": total,
                       "released": by_status.get("released", 0 if ok else -1),
                       "http_status": st, "filter_field": spec["field"],
                       "url": portal.url(base + [("limit", "0")])})
        for field, terms in ft.items():
            for term_v, count in terms:
                facet_rows.append({"portal": portal.name, "type": t,
                                   "facet": field, "term": term_v,
                                   "count": count})
        stats.update({f"{t}.{k}": v for k, v in facet_stats(data).items()})
        if not (spec.get("items") and ok and total > 0):
            continue
        params = list(base)
        if status_filter != "all":
            params.append(("status", status_filter))
        params.append(("limit", str(max_items)))
        fields = spec.get("item_fields") or IGVF_ITEM_FIELDS
        params.extend(("field", f) for f in fields)
        st2, d2 = portal.get(params)
        rows = graph_rows(d2) if not is_zero_hit(st2, d2) else []
        items[t] = [flatten_item(r, fields) for r in rows]
    return {"totals": totals, "facets": facet_rows, "items": items,
            "stats": stats}


def _scalar(v: Any) -> str:
    """Portal fields are inconsistently scalar, list, or embedded object."""
    if v is None:
        return ""
    if isinstance(v, dict):
        for k in ("term_name", "label", "title", "@id", "name", "summary"):
            if v.get(k):
                return _scalar(v[k])
        return ""
    if isinstance(v, list):
        vals = [s for s in (_scalar(x) for x in v) if s]
        return "; ".join(dict.fromkeys(vals))
    return str(v)


def _dig(row: dict, dotted: str) -> Any:
    """Follow a dotted path through dicts and lists, flattening as it goes."""
    cur: Any = row
    for part in dotted.split("."):
        if isinstance(cur, dict):
            cur = cur.get(part)
            continue
        if isinstance(cur, list):
            nxt: "list[Any]" = []
            for c in cur:
                if isinstance(c, dict):
                    v = c.get(part)
                    nxt.extend(v if isinstance(v, list) else [v])
            cur = nxt
            continue
        return None
    return cur


def flatten_item(row: dict, fields: "list[str]") -> "dict[str, str]":
    out = {}
    for f in fields:
        val = _dig(row, f)
        if f == "@type" and isinstance(val, list):
            val = val[0] if val else ""
        out[f] = _scalar(val)
    return out


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def write_csv(path: Path, rows: "list[dict]", cols: "Optional[list[str]]" = None) -> Path:
    cols = cols or (list(rows[0].keys()) if rows else [])
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def top_terms(facet_rows: "list[dict]", portal: str, obj_type: str,
              facet: str, top: int) -> "list[tuple[str, int]]":
    rows = [(r["term"], r["count"]) for r in facet_rows
            if r["portal"] == portal and r["type"] == obj_type and r["facet"] == facet]
    rows.sort(key=lambda kv: (-kv[1], kv[0]))
    if len(rows) > top:
        rest = sum(c for _, c in rows[top:])
        rows = rows[:top] + [(f"Other ({len(rows) - top} more)", rest)]
    return rows


def by_year(items: "list[dict]", key: str = "date_released") -> "list[tuple[str, int]]":
    c: "collections.Counter[str]" = collections.Counter()
    for r in items:
        v = (r.get(key) or "")[:4]
        if v.isdigit():
            c[v] += 1
    if not c:
        return []
    # A year with no releases is still a year: leaving it out of the axis
    # made 2015 and 2024 vanish from the GM12878 timeline.
    lo, hi = int(min(c)), int(max(c))
    return [(str(y), c.get(str(y), 0)) for y in range(lo, hi + 1)]


def count_items(items: "list[dict]", key: str, top: int) -> "list[tuple[str, int]]":
    c: "collections.Counter[str]" = collections.Counter()
    for r in items:
        for part in (r.get(key) or "").split("; "):
            if part:
                c[part] += 1
    rows = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))
    if len(rows) > top:
        rest = sum(v for _, v in rows[top:])
        rows = rows[:top] + [(f"Other ({len(rows) - top} more)", rest)]
    return rows


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]") -> str:
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Figures -- SVG by hand (always), PNG via matplotlib (when present)
# ---------------------------------------------------------------------------

def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _fmt(n: Any) -> str:
    return f"{n:,}" if isinstance(n, int) else str(n)


def hbar_svg_panel(items: "list[tuple[str, int]]", color: str, *, x0: int,
                   y0: int, width: int, title: str) -> "tuple[str, int]":
    """One horizontal-bar panel; returns (svg fragment, height used)."""
    row_h, gap, label_w = 18, 2, 230
    n = len(items)
    height = 30 + n * (row_h + gap) + 10
    vmax = max([c for _, c in items] + [1])
    plot_w = width - label_w - 90
    parts = [f'<text x="{x0}" y="{y0 + 14}" font-size="13" font-weight="600" fill="{INK}">{_esc(title)}</text>']
    y = y0 + 30
    for name, count in items:
        w = max(1, int(plot_w * count / vmax)) if count > 0 else 0
        parts.append(f'<text x="{x0 + label_w - 8}" y="{y + row_h - 5}" font-size="11" '
                     f'text-anchor="end" fill="{INK2}">{_esc(str(name)[:38])}</text>')
        parts.append(f'<rect x="{x0 + label_w}" y="{y + gap}" width="{w}" '
                     f'height="{row_h - gap}" rx="3" fill="{color}"/>')
        parts.append(f'<text x="{x0 + label_w + w + 6}" y="{y + row_h - 5}" '
                     f'font-size="11" fill="{INK}">{_fmt(count)}</text>')
        y += row_h + gap
    parts.append(f'<line x1="{x0 + label_w}" y1="{y0 + 28}" x2="{x0 + label_w}" '
                 f'y2="{y}" stroke="{AXIS}" stroke-width="1"/>')
    return "\n".join(parts), height


def write_hbar_svg(items: "list[tuple[str, int]]", title: str, path: Path,
                   color: str, subtitle: str = "") -> Path:
    width = 760
    body, h = hbar_svg_panel(items, color, x0=20, y0=44, width=width - 20, title=subtitle)
    total_h = 44 + h
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{total_h}" '
           f'viewBox="0 0 {width} {total_h}" font-family="Inter, Helvetica, Arial, sans-serif">',
           f'<rect width="{width}" height="{total_h}" fill="{SURFACE}"/>',
           f'<text x="20" y="26" font-size="16" font-weight="700" fill="{INK}">{_esc(title)}</text>',
           body, "</svg>"]
    path.write_text("\n".join(svg), encoding="utf-8")
    return path


def write_panels_svg(panels: "list[tuple[str, list[tuple[str, int]], str]]",
                     title: str, path: Path) -> Path:
    """Small multiples side by side: totals differ by orders of magnitude
    across portals, so each panel carries its own scale (one axis each,
    never a shared or dual axis)."""
    pw, width = 560, 560 * max(1, len(panels)) + 20
    frags, hmax = [], 0
    for i, (sub, items, color) in enumerate(panels):
        frag, h = hbar_svg_panel(items or [("(nothing found)", 0)], color,
                                 x0=20 + i * pw, y0=44, width=pw - 20, title=sub)
        frags.append(frag)
        hmax = max(hmax, h)
    total_h = 44 + hmax
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{total_h}" '
           f'viewBox="0 0 {width} {total_h}" font-family="Inter, Helvetica, Arial, sans-serif">',
           f'<rect width="{width}" height="{total_h}" fill="{SURFACE}"/>',
           f'<text x="20" y="26" font-size="16" font-weight="700" fill="{INK}">{_esc(title)}</text>',
           *frags, "</svg>"]
    path.write_text("\n".join(svg), encoding="utf-8")
    return path


def write_vbar_svg(items: "list[tuple[str, int]]", title: str, path: Path,
                   color: str, x_label: str) -> Path:
    width, height, left, bottom, top = 760, 320, 60, 60, 44
    n = max(1, len(items))
    plot_w, plot_h = width - left - 20, height - top - bottom
    vmax = max([c for _, c in items] + [1])
    slot = plot_w / n
    bw = max(2, int(slot) - 4)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'viewBox="0 0 {width} {height}" font-family="Inter, Helvetica, Arial, sans-serif">',
             f'<rect width="{width}" height="{height}" fill="{SURFACE}"/>',
             f'<text x="20" y="26" font-size="16" font-weight="700" fill="{INK}">{_esc(title)}</text>',
             f'<line x1="{left}" y1="{top + plot_h}" x2="{width - 20}" y2="{top + plot_h}" stroke="{AXIS}"/>']
    for i, (name, count) in enumerate(items):
        h = int(plot_h * count / vmax)
        x = left + int(i * slot) + 2
        y = top + plot_h - h
        parts.append(f'<rect x="{x}" y="{y}" width="{bw}" height="{h}" rx="3" fill="{color}"/>')
        parts.append(f'<text x="{x + bw / 2:.0f}" y="{y - 4}" font-size="10" text-anchor="middle" fill="{INK}">{_fmt(count)}</text>')
        if n <= 30 or i % 2 == 0:
            parts.append(f'<text x="{x + bw / 2:.0f}" y="{top + plot_h + 14}" font-size="10" '
                         f'text-anchor="middle" fill="{INK2}">{_esc(name)}</text>')
    parts.append(f'<text x="{left + plot_w / 2:.0f}" y="{height - 12}" font-size="11" '
                 f'text-anchor="middle" fill="{INK2}">{_esc(x_label)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def try_png(svg_path: Path, items: "list[tuple[str, int]]", title: str,
            color: str, horizontal: bool = True) -> "Optional[Path]":
    """A PNG twin of the SVG when matplotlib is importable; None otherwise.

    Kept optional on purpose: the metadata skills run on the standard library
    and the SVG is the artefact of record.
    """
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return None
    if not items:
        return None
    names = [str(n)[:40] for n, _ in items]
    vals = [c for _, c in items]
    try:
        if horizontal:
            fig, ax = plt.subplots(figsize=(9, max(2.5, 0.32 * len(items) + 1)))
            ax.barh(names[::-1], vals[::-1], color=color, height=0.72)
            for i, v in enumerate(vals[::-1]):
                ax.text(v, i, f" {v:,}", va="center", fontsize=8, color=INK)
            ax.set_xlabel("count", color=INK2)
        else:
            fig, ax = plt.subplots(figsize=(9, 3.6))
            ax.bar(names, vals, color=color, width=0.72)
            for i, v in enumerate(vals):
                ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8, color=INK)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
        ax.tick_params(colors=INK2, labelsize=8)
        fig.patch.set_facecolor(SURFACE)
        ax.set_facecolor(SURFACE)
        fig.tight_layout()
        out = svg_path.with_suffix(".png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        return out
    except Exception as exc:  # pragma: no cover - plotting backend quirks
        logging.warning("PNG for %s skipped: %s", svg_path.name, exc)
        return None


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_report(*, biosample: str, terms: "dict[str, dict]", status_filter: str,
                 totals: "list[dict]", facets: "list[dict]",
                 items: "dict[str, dict[str, list[dict]]]", stats: "dict[str, dict]",
                 figures: "list[Path]", tables: "list[Path]", top: int,
                 out_dir: Path, portals: "list[str]") -> str:
    L: "list[str]" = []
    L.append(f"# Biosample census: {biosample}")
    L.append("")
    L.append(f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `igvfagent biosample-census`. "
             f"Portals: {', '.join(p.upper() for p in portals)}. "
             f"Item tables use status `{status_filter}`; totals and facet counts cover all statuses.")
    L.append("")
    L.append("## 1. Identity on each portal")
    L.append("")
    rows = []
    for p in portals:
        t = terms.get(p, {})
        rows.append([p.upper(), t.get("term_name", biosample), t.get("term_id") or "-",
                     t.get("classification") or "-",
                     "matched ontology term" if t.get("resolved") else
                     "no ontology match; filtered on the text as given"])
    L.append(md_table(["Portal", "term_name used", "term_id", "classification", "resolution"], rows))
    L.append("")
    L.append("## 2. Totals by object type")
    L.append("")
    L.append("Counts come from the portal's structured sample-term field (column *filter field*), "
             "not free-text search. A dash means the portal could not be reached for that type.")
    L.append("")
    L.append(md_table(["Portal", "Type", "Total (all statuses)", "Released", "Filter field"],
                      [[r["portal"].upper(), r["type"],
                        _fmt(r["total"]) if r["total"] >= 0 else "-",
                        _fmt(r["released"]) if r["released"] >= 0 else "-",
                        r["filter_field"]] for r in totals]))
    L.append("")

    sec = 3

    def facet_section(title: str, portal: str, obj_type: str, facet: str,
                      note: str = "") -> None:
        nonlocal sec
        rows_ = top_terms(facets, portal, obj_type, facet, top)
        if not rows_:
            return
        L.append(f"## {sec}. {title}")
        sec += 1
        L.append("")
        if note:
            L.append(note)
            L.append("")
        L.append(md_table([facet, "count"], [[a, _fmt(b)] for a, b in rows_]))
        L.append("")

    if "encode" in portals:
        facet_section("ENCODE experiments by assay", "encode", "Experiment", "assay_title")
        facet_section("ENCODE ChIP-seq targets", "encode", "Experiment", "target.label",
                      "Targets are shown for experiments that have one (TF and histone ChIP-seq, eCLIP, etc.).")
        facet_section("ENCODE experiments by lab", "encode", "Experiment", "lab.title")
        facet_section("ENCODE experiments by assembly", "encode", "Experiment", "assembly")
        facet_section("ENCODE treatments", "encode", "Experiment",
                      "replicates.library.biosample.treatments.treatment_term_name")
        facet_section("ENCODE functional characterization experiments", "encode",
                      "FunctionalCharacterizationExperiment", "assay_title")
        facet_section("ENCODE annotations by type", "encode", "Annotation", "annotation_type",
                      "Annotations are derived products (cCREs, chromatin states, enhancer-gene links, ...) rather than experiments.")
        facet_section("ENCODE files by format", "encode", "File", "file_format")
        facet_section("ENCODE files by output type", "encode", "File", "output_type")
        exp_items = items.get("encode", {}).get("Experiment", [])
        yrs = by_year(exp_items)
        if yrs:
            L.append(f"## {sec}. ENCODE experiment releases by year")
            sec += 1
            L.append("")
            L.append(md_table(["year", f"experiments ({status_filter})"], [[y, _fmt(c)] for y, c in yrs]))
            L.append("")

    if "igvf" in portals:
        ig_items = items.get("igvf", {})
        set_rows = []
        for t, rows_ in ig_items.items():
            if t in ("File", "Sample"):
                continue
            for r in rows_:
                set_rows.append([r.get("accession", ""), t,
                                 r.get("preferred_assay_titles") or r.get("assay_titles") or r.get("file_set_type") or "",
                                 r.get("lab.title", ""), r.get("status", ""),
                                 (r.get("release_timestamp") or "")[:10],
                                 (r.get("summary") or "")[:110]])
        if set_rows:
            L.append(f"## {sec}. IGVF data sets ({status_filter})")
            sec += 1
            L.append("")
            L.append(md_table(["accession", "type", "assay / set type", "lab", "status", "released", "summary"], set_rows))
            L.append("")
        facet_section("IGVF measurement sets by assay", "igvf", "MeasurementSet", "preferred_assay_titles")
        facet_section("IGVF measurement sets by lab", "igvf", "MeasurementSet", "lab.title")
        facet_section("IGVF files by format", "igvf", "File", "file_format")
        facet_section("IGVF files by content type", "igvf", "File", "content_type")
        samp = ig_items.get("Sample", [])
        if samp:
            L.append(f"## {sec}. IGVF samples ({status_filter})")
            sec += 1
            L.append("")
            L.append(md_table(["accession", "classification", "lab", "status", "summary"],
                              [[r.get("accession", ""), r.get("classifications", ""),
                                r.get("lab.title", ""), r.get("status", ""),
                                (r.get("summary") or "")[:110]] for r in samp]))
            L.append("")
        fs = stats.get("File.file_size")
        if fs and fs.get("sum"):
            L.append(f"IGVF file volume for this sample: {fs['sum'] / 1e9:,.2f} GB across "
                     f"{int(fs.get('count', 0)):,} files (portal `file_size` facet).")
            L.append("")

    L.append(f"## {sec}. Figures")
    sec += 1
    L.append("")
    for f in figures:
        L.append(f"- `{f.name}`")
        if f.suffix == ".svg":
            L.append(f"  ![{f.stem}]({f.name})")
    L.append("")
    L.append(f"## {sec}. Tables")
    sec += 1
    L.append("")
    for t in tables:
        L.append(f"- `{t.name}`")
    L.append("")
    L.append(f"## {sec}. Method and caveats")
    L.append("")
    L.append(f"- ENCODE base `{ENCODE_BASE}`; IGVF Portal API base `{IGVF_API_BASE}`. Every request "
             f"made, with its HTTP status and total, is in `census.json`.")
    L.append("- Both portals answer a search with no matches as HTTP 404 with `total: 0`; the census records that as zero.")
    L.append("- Free-text search is **not** used for counting: on the IGVF Portal `searchTerm=GM12878` matches "
             "thousands of unrelated sets, and on ENCODE it counts any record whose text mentions the line. "
             "Structured sample-term filters are exact.")
    L.append("- ENCODE `File` counts include every processed and raw file attached to the biosample's datasets; "
             "IGVF `File` counts are files whose file set has this sample.")
    L.append("- The IGVF Portal shows only released data without credentials; set `IGVF_ACCESS_KEY` / "
             "`IGVF_SECRET_ACCESS_KEY` to include in-progress sets you are entitled to see.")
    L.append(f"- Output directory: `{out_dir}`")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    portals = ["encode", "igvf"] if args.portal == "both" else [args.portal]
    clients = {"encode": Portal("encode", ENCODE_BASE), "igvf": Portal("igvf", IGVF_API_BASE)}
    ts = time.strftime("%Y%m%d_%H%M%S")
    label = safe_label(args.label or f"{args.biosample}_census")
    out_dir = OUT_ROOT / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    terms: "dict[str, dict]" = {}
    totals: "list[dict]" = []
    facets: "list[dict]" = []
    items: "dict[str, dict[str, list[dict]]]" = {}
    stats: "dict[str, dict]" = {}
    for p in portals:
        c = clients[p]
        hint = terms.get("encode", {}).get("term_id", "") if p == "igvf" else ""
        terms[p] = resolve_term(c, args.biosample, term_id_hint=hint)
        res = census_portal(c, ENCODE_TYPES if p == "encode" else IGVF_TYPES,
                            terms[p]["term_name"], status_filter=args.status,
                            max_items=args.max_items)
        totals += res["totals"]
        facets += res["facets"]
        items[p] = res["items"]
        stats.update(res["stats"])

    # ---- tables
    tables: "list[Path]" = []
    tables.append(write_csv(out_dir / "totals_by_type.csv", totals,
                            ["portal", "type", "total", "released", "filter_field", "http_status", "url"]))
    tables.append(write_csv(out_dir / "facet_counts.csv", facets,
                            ["portal", "type", "facet", "term", "count"]))
    for p, per_type in items.items():
        for t, rows in per_type.items():
            if rows:
                tables.append(write_csv(out_dir / f"{p}_{t.lower()}_items.csv", rows))
    for p in portals:
        for spec in (ENCODE_TYPES if p == "encode" else IGVF_TYPES):
            for facet in spec["facets"]:
                rows = top_terms(facets, p, spec["type"], facet, top=10 ** 6)
                if rows:
                    name = f"table_{p}_{spec['type'].lower()}_{facet.replace('.', '_')}.csv"
                    tables.append(write_csv(out_dir / name,
                                            [{"term": a, "count": b} for a, b in rows]))

    # ---- figures
    figures: "list[Path]" = []
    if not args.no_plots:
        def fig(fname: str, items_: "list[tuple[str, int]]", title: str,
                portal: str, horizontal: bool = True, x_label: str = "") -> None:
            if not items_:
                return
            path = out_dir / fname
            if horizontal:
                write_hbar_svg(items_, title, path, COLOR[portal], subtitle="count")
            else:
                write_vbar_svg(items_, title, path, COLOR[portal], x_label)
            figures.append(path)
            png = try_png(path, items_, title, COLOR[portal], horizontal)
            if png:
                figures.append(png)

        panels = []
        for p in portals:
            panels.append((f"{p.upper()} objects", [(r["type"], r["total"]) for r in totals
                                                    if r["portal"] == p and r["total"] > 0], COLOR[p]))
        if any(items_ for _, items_, _ in panels):
            path = write_panels_svg(panels, f"{args.biosample}: what each portal holds, by object type",
                                    out_dir / "fig1_totals_by_type.svg")
            figures.append(path)
        if "encode" in portals:
            fig("fig2_encode_assays.svg", top_terms(facets, "encode", "Experiment", "assay_title", args.top),
                f"ENCODE {args.biosample} experiments by assay", "encode")
            fig("fig3_encode_targets.svg", top_terms(facets, "encode", "Experiment", "target.label", args.top),
                f"ENCODE {args.biosample} ChIP-seq targets", "encode")
            fig("fig4_encode_labs.svg", top_terms(facets, "encode", "Experiment", "lab.title", args.top),
                f"ENCODE {args.biosample} experiments by lab", "encode")
            fig("fig5_encode_release_timeline.svg", by_year(items.get("encode", {}).get("Experiment", [])),
                f"ENCODE {args.biosample} experiment releases by year ({args.status})", "encode",
                horizontal=False, x_label="year released")
            fig("fig6_encode_annotation_types.svg",
                top_terms(facets, "encode", "Annotation", "annotation_type", args.top),
                f"ENCODE {args.biosample} annotations by type", "encode")
            fig("fig7_encode_file_formats.svg", top_terms(facets, "encode", "File", "file_format", args.top),
                f"ENCODE {args.biosample} files by format", "encode")
        if "igvf" in portals:
            ig_sets = [r for t, rows in items.get("igvf", {}).items() if t not in ("File", "Sample") for r in rows]
            fig("fig8_igvf_assays.svg", count_items(ig_sets, "preferred_assay_titles", args.top)
                or top_terms(facets, "igvf", "MeasurementSet", "preferred_assay_titles", args.top),
                f"IGVF {args.biosample} data sets by assay ({args.status})", "igvf")
            fig("fig9_igvf_labs.svg", count_items(ig_sets, "lab.title", args.top)
                or top_terms(facets, "igvf", "MeasurementSet", "lab.title", args.top),
                f"IGVF {args.biosample} data sets by lab ({args.status})", "igvf")
            fig("fig10_igvf_file_formats.svg", top_terms(facets, "igvf", "File", "file_format", args.top),
                f"IGVF {args.biosample} files by format", "igvf")

    # ---- report + audit trail
    report = build_report(biosample=args.biosample, terms=terms, status_filter=args.status,
                          totals=totals, facets=facets, items=items, stats=stats,
                          figures=figures, tables=tables, top=args.top, out_dir=out_dir,
                          portals=portals)
    report_path = out_dir / "report.md"
    report_path.write_text(report, encoding="utf-8")
    audit = {"biosample": args.biosample, "terms": terms, "status_filter": args.status,
             "totals": totals, "stats": stats,
             "requests": {p: clients[p].log for p in portals},
             "n_items": {p: {t: len(r) for t, r in per.items()} for p, per in items.items()}}
    json_path = out_dir / "census.json"
    json_path.write_text(json.dumps(audit, indent=2, default=str), encoding="utf-8")

    # ---- stdout: identity first, then numbers, then paths the agent can open
    for p in portals:
        t = terms[p]
        how = f"term_id {t['term_id']}" if t.get("resolved") else "no ontology match, text filter"
        print(f"{p.upper()} biosample: {t['term_name']} ({how})")
    unreachable = [r for r in totals if r["total"] < 0]
    for r in totals:
        if r["total"] >= 0:
            print(f"  {r['portal'].upper():6s} {r['type']:38s} {r['total']:>8,}  (released {r['released']:,})")
    for r in unreachable:
        print(f"  {r['portal'].upper():6s} {r['type']:38s} unreachable (HTTP {r['http_status']})")
    for p in portals:
        for t, rows in items.get(p, {}).items():
            if rows:
                print(f"  {p.upper()} {t} items fetched ({args.status}): {len(rows)}")
    print(f"Report: {report_path}")
    for tpath in tables[:4]:
        print(f"CSV: {tpath}")
    if len(tables) > 4:
        print(f"CSV: ... {len(tables) - 4} more table(s) in {out_dir}")
    for f in figures:
        if f.suffix == ".svg":
            print(f"Figure: {f}")
    if any(f.suffix == ".png" for f in figures):
        print(f"PNG twins written beside each SVG ({sum(1 for f in figures if f.suffix == '.png')}).")
    print(f"JSON: {json_path}")
    print(f"Log: {log_path}")
    if unreachable and len(unreachable) == len(totals):
        print("ERROR: no portal answered; nothing was counted.", file=sys.stderr)
        return 1
    return 0


def main(argv: "Optional[list[str]]" = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "run":       # `biosample-census run ...` also works
        argv = argv[1:]
    p = argparse.ArgumentParser(
        prog="igvfagent biosample-census",
        description="Systematic ENCODE + IGVF census for one biosample / cell line, "
                    "with tables and plots.")
    p.add_argument("--biosample", required=True,
                   help="Cell line or tissue term, e.g. GM12878, K562, HepG2, liver.")
    p.add_argument("--label", default="", help="Run-directory label.")
    p.add_argument("--portal", choices=("both", "encode", "igvf"), default="both")
    p.add_argument("--status", default="released",
                   help="Status filter for item tables (default released; 'all' for every status).")
    p.add_argument("--max-items", type=int, default=5000,
                   help="Cap on item rows fetched per object type.")
    p.add_argument("--top", type=int, default=20, help="Top-N terms per table / figure.")
    p.add_argument("--no-plots", action="store_true", help="Tables and report only.")
    return run(p.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
