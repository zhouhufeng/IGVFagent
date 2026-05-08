#!/usr/bin/env python3
"""Reference Skill: literature retrieval, validation, and study design.

Three top-level functions, each exposed as a subcommand:

  1. ``learn``    — retrieve publications/preprints relevant to a topic, assay,
                    or biosample, and summarize what other studies do, how
                    they design experiments, and what plots / tables they
                    typically produce.

  2. ``validate`` — given a list of genes, variants, or regulatory elements
                    discovered by other IGVF-agent skills, search prior
                    literature to check whether the discoveries have been
                    independently reported and how strong the prior evidence
                    is.

  3. ``design``   — given an IGVF data type (assay name or accession), pull
                    cognate published work, recommend a study workflow based
                    on what successful prior studies did, and surface IGVF
                    Catalog/Portal datasets that fit the workflow.

Sources queried (all public, all respected with conservative rate limiting):
  • PubMed / PMC via NCBI E-utilities
  • bioRxiv + medRxiv API
  • arXiv API
  • Semantic Scholar Graph API
  • OpenAlex
  • Crossref (DOI lookup + metadata)

Outputs follow the project pattern: cached responses under
``Data/Cache/References/<source>/``, manifest CSVs under
``Data/Manifests/References/``, and reports under
``Docs/References/<timestamp>/``. Logs in ``Docs/Logs/``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "References"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "References"
CACHE_DIR = DATA_DIR / "Cache" / "References"

USER_AGENT = (
    "IGVFagent-ReferenceSkill/1.0 "
    "(local research tooling; contact: configured locally)"
)

# Conservative pacing per source (seconds between requests)
SOURCE_DELAY = {
    "pubmed":          0.34,   # NCBI: 3 req/s without API key
    "biorxiv":         0.5,
    "arxiv":           3.0,    # arXiv recommends 3-second delay
    "semanticscholar": 1.0,    # 1 req/s on free tier
    "openalex":        0.2,    # OpenAlex is generous; still polite
    "crossref":        0.4,
}

# Curated journals worth highlighting in summaries
TOP_JOURNALS = {
    "Nature", "Nat Methods", "Nat Genet", "Nat Biotechnol", "Nat Commun",
    "Nat Rev Genet", "Nat Cell Biol", "Nat Neurosci", "Nat Med",
    "Cell", "Cell Genom", "Mol Cell", "Cell Syst", "Cell Rep", "Neuron",
    "Cell Stem Cell", "Cancer Cell", "Immunity", "Cell Metab",
    "Science", "Sci Adv", "Sci Transl Med",
    "N Engl J Med", "Lancet",
    "Nucleic Acids Res", "Bioinformatics", "Brief Bioinform",
    "Genome Biol", "Genome Res", "PLoS Comput Biol", "PLoS Genet",
    "eLife", "PNAS",
}

# Canonical visualization vocabulary used to tag papers and to suggest plots
PLOT_VOCAB = {
    # genomics-wide
    "umap", "tsne", "pca", "heatmap", "violin", "boxplot", "scatter",
    # variant + GWAS
    "manhattan", "qq", "miami", "locuszoom", "forest", "volcano",
    # epigenomics
    "track", "metaplot", "tornado", "circos", "chromhmm",
    # single-cell
    "dotplot", "stacked bar", "ridgeplot", "feature plot",
    # networks / graphs
    "network", "venn", "sankey", "alluvial", "upset",
}

# ----------------------------- General utilities ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"reference_skill_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (REPORT_DIR, MANIFEST_DIR, CACHE_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)
    for sub in ("pubmed", "biorxiv", "arxiv", "semanticscholar",
                "openalex", "crossref"):
        (CACHE_DIR / sub).mkdir(parents=True, exist_ok=True)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def cache_key(*parts: str) -> str:
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:24]


def http_get(url: str, params: dict | None = None, headers: dict | None = None,
             accept_json: bool = True, timeout: int = 60,
             source_delay_key: str | None = None) -> tuple[int, str]:
    if params:
        sep = "&" if "?" in url else "?"
        url = f"{url}{sep}{urllib.parse.urlencode(params)}"
    h = {"User-Agent": USER_AGENT}
    if accept_json:
        h["Accept"] = "application/json"
    if headers:
        h.update(headers)
    if source_delay_key and SOURCE_DELAY.get(source_delay_key):
        time.sleep(SOURCE_DELAY[source_delay_key])
    req = urllib.request.Request(url, headers=h)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read().decode("utf-8", errors="replace")
            logging.info("GET %s -> %d (%d bytes)", url, resp.status, len(data))
            return resp.status, data
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace") if e.fp else ""
        logging.warning("GET %s -> HTTP %d", url, e.code)
        return e.code, body
    except Exception as e:
        logging.error("GET %s failed: %s", url, e)
        return 0, ""


def cached_get(source: str, url: str, params: dict | None = None,
               headers: dict | None = None, accept_json: bool = True,
               ttl_days: int = 14) -> tuple[int, str]:
    key_str = url + "|" + json.dumps(params or {}, sort_keys=True)
    key = cache_key(source, key_str)
    cache_path = CACHE_DIR / source / f"{key}.txt"
    if cache_path.exists():
        age_days = (time.time() - cache_path.stat().st_mtime) / 86400
        if age_days <= ttl_days:
            return 200, cache_path.read_text()
    status, body = http_get(url, params=params, headers=headers,
                            accept_json=accept_json,
                            source_delay_key=source)
    if status == 200 and body:
        cache_path.write_text(body)
    return status, body


# ------------------------------- PubMed client --------------------------------

def pubmed_search(query: str, limit: int = 25,
                  date_from: str | None = None,
                  date_to: str | None = None,
                  journals: list[str] | None = None) -> list[dict]:
    """Query PubMed via NCBI E-utilities. Returns a list of records with
    title, authors, journal, year, doi, pmid, pmcid, abstract."""
    base = _resolve_endpoint("pubmed_eutils", "PUBMED_EUTILS_BASE")

    parts = [f"({query})"]
    if journals:
        joined = " OR ".join(f"\"{j}\"[Journal]" for j in journals)
        parts.append(f"({joined})")
    if date_from or date_to:
        df = (date_from or "1900").replace("-", "/")
        dt = (date_to or "3000").replace("-", "/")
        parts.append(f"(\"{df}\"[Date - Publication] : \"{dt}\"[Date - Publication])")
    full_query = " AND ".join(parts)

    status, body = cached_get(
        "pubmed", f"{base}/esearch.fcgi",
        params={"db": "pubmed", "term": full_query, "retmax": limit,
                "retmode": "json", "sort": "relevance"},
    )
    if status != 200 or not body:
        logging.warning("PubMed esearch returned %d", status)
        return []
    try:
        ids = json.loads(body).get("esearchresult", {}).get("idlist", [])
    except Exception:
        return []
    if not ids:
        return []

    status, body = cached_get(
        "pubmed", f"{base}/esummary.fcgi",
        params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
    )
    summaries: dict[str, dict] = {}
    if status == 200 and body:
        try:
            data = json.loads(body).get("result", {})
            for pid in ids:
                if pid in data:
                    summaries[pid] = data[pid]
        except Exception:
            pass

    status, body = cached_get(
        "pubmed", f"{base}/efetch.fcgi",
        params={"db": "pubmed", "id": ",".join(ids), "retmode": "xml",
                "rettype": "abstract"},
        accept_json=False,
    )
    abstracts: dict[str, str] = {}
    if status == 200 and body:
        try:
            root = ET.fromstring(body)
            for art in root.iter("PubmedArticle"):
                pmid_el = art.find(".//PMID")
                if pmid_el is None:
                    continue
                pmid = pmid_el.text or ""
                texts = []
                for ab in art.iter("AbstractText"):
                    label = ab.attrib.get("Label", "")
                    txt = "".join(ab.itertext()).strip()
                    if not txt:
                        continue
                    texts.append(f"[{label}] {txt}" if label else txt)
                if texts:
                    abstracts[pmid] = " ".join(texts)
        except ET.ParseError:
            pass

    results = []
    for pid in ids:
        s = summaries.get(pid, {})
        authors = s.get("authors", []) or []
        author_names = [a.get("name", "") for a in authors[:5]]
        ids_arr = s.get("articleids", []) or []
        doi = next((x.get("value", "") for x in ids_arr if x.get("idtype") == "doi"), "")
        pmcid = next((x.get("value", "") for x in ids_arr if x.get("idtype") == "pmc"), "")
        results.append({
            "source": "pubmed",
            "pmid": pid,
            "doi": doi,
            "pmcid": pmcid,
            "title": s.get("title", ""),
            "journal": s.get("fulljournalname") or s.get("source", ""),
            "year": (s.get("pubdate", "") or "")[:4],
            "authors": "; ".join(author_names) + ("…" if len(authors) > 5 else ""),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
            "abstract": abstracts.get(pid, "")[:2400],
        })
    return results


# ------------------------------ bioRxiv client --------------------------------

def biorxiv_search(query: str, limit: int = 25, server: str = "biorxiv",
                   days_back: int = 730) -> list[dict]:
    """Crossref-backed search of bioRxiv/medRxiv preprints."""
    base = _resolve_endpoint("crossref", "CROSSREF_API_BASE")
    member = "246" if server == "biorxiv" else "31795"  # Cold Spring Harbor / medRxiv
    today = time.strftime("%Y-%m-%d")
    from_date = time.strftime(
        "%Y-%m-%d",
        time.gmtime(time.time() - days_back * 86400),
    )
    status, body = cached_get(
        "crossref", f"{base}/works",
        params={
            "query": query,
            "filter": (f"member:{member},from-pub-date:{from_date},"
                        f"until-pub-date:{today}"),
            "rows": limit,
            "select": "title,author,issued,DOI,abstract,container-title,subtype",
        },
    )
    if status != 200 or not body:
        return []
    try:
        items = json.loads(body).get("message", {}).get("items", [])
    except Exception:
        return []
    out = []
    for it in items:
        title = (it.get("title") or [""])[0]
        authors = it.get("author", []) or []
        author_names = [
            f"{a.get('given','')} {a.get('family','')}".strip()
            for a in authors[:5]
        ]
        issued = (it.get("issued") or {}).get("date-parts", [[None]])
        year = str(issued[0][0]) if issued and issued[0][0] else ""
        doi = it.get("DOI", "")
        abs_html = it.get("abstract", "") or ""
        abstract = re.sub(r"<[^>]+>", " ", abs_html).strip()
        out.append({
            "source": server,
            "doi": doi,
            "title": title,
            "journal": (it.get("container-title") or [server])[0],
            "year": year,
            "authors": "; ".join(author_names) + ("…" if len(authors) > 5 else ""),
            "url": f"https://doi.org/{doi}" if doi else "",
            "abstract": abstract[:2400],
        })
    return out


# ------------------------------- arXiv client ---------------------------------

def arxiv_search(query: str, limit: int = 25,
                 categories: list[str] | None = None) -> list[dict]:
    base = _resolve_endpoint("arxiv_api", "ARXIV_API_BASE")
    q = f'all:"{query}"'
    if categories:
        cat = " OR ".join(f"cat:{c}" for c in categories)
        q = f"({q}) AND ({cat})"
    status, body = cached_get(
        "arxiv", base,
        params={"search_query": q, "start": 0, "max_results": limit,
                "sortBy": "relevance", "sortOrder": "descending"},
        accept_json=False,
    )
    if status != 200 or not body:
        return []
    out = []
    try:
        ns = {"a": "http://www.w3.org/2005/Atom"}
        root = ET.fromstring(body)
        for e in root.findall("a:entry", ns):
            title = (e.findtext("a:title", default="", namespaces=ns) or "").strip()
            summary = (e.findtext("a:summary", default="", namespaces=ns) or "").strip()
            published = e.findtext("a:published", default="", namespaces=ns) or ""
            authors = [a.findtext("a:name", default="", namespaces=ns) or ""
                       for a in e.findall("a:author", ns)]
            link_id = e.findtext("a:id", default="", namespaces=ns) or ""
            arxiv_id = link_id.rsplit("/", 1)[-1]
            out.append({
                "source": "arxiv",
                "arxiv_id": arxiv_id,
                "doi": "",
                "title": re.sub(r"\s+", " ", title),
                "journal": "arXiv",
                "year": published[:4],
                "authors": "; ".join(authors[:5]) + ("…" if len(authors) > 5 else ""),
                "url": link_id,
                "abstract": summary[:2400],
            })
    except ET.ParseError:
        pass
    return out


# --------------------------- Semantic Scholar client --------------------------

def semanticscholar_search(query: str, limit: int = 25) -> list[dict]:
    base = _resolve_endpoint("semanticscholar", "SS_API_BASE")
    fields = ("title,abstract,authors.name,venue,year,externalIds,url,"
              "citationCount,influentialCitationCount,referenceCount,tldr")
    status, body = cached_get(
        "semanticscholar", f"{base}/paper/search",
        params={"query": query, "limit": min(limit, 100), "fields": fields},
    )
    if status != 200 or not body:
        return []
    try:
        items = json.loads(body).get("data", []) or []
    except Exception:
        return []
    out = []
    for it in items:
        ext = it.get("externalIds", {}) or {}
        authors = [a.get("name", "") for a in (it.get("authors") or [])[:5]]
        tldr = (it.get("tldr") or {}).get("text", "")
        out.append({
            "source": "semanticscholar",
            "ss_id": it.get("paperId", ""),
            "doi": ext.get("DOI", ""),
            "pmid": ext.get("PubMed", ""),
            "title": it.get("title", ""),
            "journal": it.get("venue", ""),
            "year": str(it.get("year") or ""),
            "authors": "; ".join(authors) + ("…" if len(it.get("authors") or []) > 5 else ""),
            "url": it.get("url", ""),
            "abstract": (it.get("abstract") or "")[:2400],
            "tldr": tldr,
            "citation_count": it.get("citationCount", 0),
            "influential_count": it.get("influentialCitationCount", 0),
        })
    return out


def semanticscholar_recommendations(paper_id: str, limit: int = 20) -> list[dict]:
    """Return papers similar to ``paper_id`` (DOI, PMID, S2 paper ID)."""
    base = _resolve_endpoint("semanticscholar", "SS_API_BASE")
    fields = "title,abstract,authors.name,venue,year,externalIds,url"
    status, body = cached_get(
        "semanticscholar",
        f"{base.replace('/graph/v1','')}/recommendations/v1/papers/forpaper/{paper_id}",
        params={"limit": limit, "fields": fields},
    )
    if status != 200 or not body:
        return []
    try:
        items = json.loads(body).get("recommendedPapers", []) or []
    except Exception:
        return []
    out = []
    for it in items:
        ext = it.get("externalIds", {}) or {}
        authors = [a.get("name", "") for a in (it.get("authors") or [])[:5]]
        out.append({
            "source": "semanticscholar_rec",
            "ss_id": it.get("paperId", ""),
            "doi": ext.get("DOI", ""),
            "pmid": ext.get("PubMed", ""),
            "title": it.get("title", ""),
            "journal": it.get("venue", ""),
            "year": str(it.get("year") or ""),
            "authors": "; ".join(authors),
            "url": it.get("url", ""),
            "abstract": (it.get("abstract") or "")[:2400],
        })
    return out


# ------------------------------- OpenAlex client ------------------------------

def openalex_search(query: str, limit: int = 25,
                    concept_ids: list[str] | None = None) -> list[dict]:
    base = _resolve_endpoint("openalex", "OPENALEX_BASE")
    flt = []
    if concept_ids:
        flt.append("concepts.id:" + "|".join(concept_ids))
    params = {"search": query, "per-page": min(limit, 50)}
    if flt:
        params["filter"] = ",".join(flt)
    status, body = cached_get("openalex", f"{base}/works", params=params)
    if status != 200 or not body:
        return []
    try:
        items = json.loads(body).get("results", [])
    except Exception:
        return []
    out = []
    for it in items[:limit]:
        ids = it.get("ids", {}) or {}
        authors = [a.get("author", {}).get("display_name", "")
                   for a in (it.get("authorships") or [])[:5]]
        host = (it.get("host_venue") or it.get("primary_location", {}) or {})
        venue = host.get("display_name") or (host.get("source") or {}).get("display_name", "")
        abstract = invert_abstract(it.get("abstract_inverted_index", {}))
        out.append({
            "source": "openalex",
            "openalex_id": ids.get("openalex", ""),
            "doi": ids.get("doi", "").replace("https://doi.org/", ""),
            "pmid": ids.get("pmid", "").rsplit("/", 1)[-1] if ids.get("pmid") else "",
            "title": it.get("title", ""),
            "journal": venue,
            "year": str(it.get("publication_year") or ""),
            "authors": "; ".join(authors),
            "url": it.get("doi") or ids.get("openalex", ""),
            "abstract": abstract[:2400],
            "citation_count": it.get("cited_by_count", 0),
        })
    return out


def invert_abstract(idx: dict) -> str:
    if not idx:
        return ""
    positions: dict[int, str] = {}
    for word, locs in idx.items():
        for p in locs:
            positions[p] = word
    if not positions:
        return ""
    return " ".join(positions[i] for i in sorted(positions))


# ------------------------------- Aggregation ---------------------------------

def dedup_records(records: Iterable[dict]) -> list[dict]:
    """Merge records that share a DOI or PMID, preferring richer entries."""
    seen: "OrderedDict[str, dict]" = OrderedDict()
    for r in records:
        key = (r.get("doi") or "").lower() or f"pmid:{r.get('pmid','')}" \
              or f"ss:{r.get('ss_id','')}" or f"ax:{r.get('arxiv_id','')}" \
              or r.get("title", "").lower()[:120]
        if not key:
            continue
        if key in seen:
            cur = seen[key]
            if not cur.get("abstract") and r.get("abstract"):
                cur["abstract"] = r["abstract"]
            for k in ("citation_count", "influential_count", "tldr"):
                if not cur.get(k) and r.get(k):
                    cur[k] = r[k]
            sources = cur.get("source", "")
            if r.get("source") and r["source"] not in sources:
                cur["source"] = sources + "+" + r["source"]
        else:
            seen[key] = dict(r)
    return list(seen.values())


def score_relevance(rec: dict, query_terms: list[str]) -> float:
    text = " ".join([rec.get("title", ""), rec.get("abstract", ""),
                     rec.get("tldr", "") or ""]).lower()
    if not text.strip():
        return 0.0
    score = 0.0
    for t in query_terms:
        t = t.lower()
        if not t:
            continue
        if t in (rec.get("title", "") or "").lower():
            score += 3.0
        if t in (rec.get("abstract", "") or "").lower():
            score += 1.0
    if (rec.get("journal") or "") in TOP_JOURNALS:
        score += 1.5
    if (rec.get("citation_count") or 0) > 100:
        score += 1.0
    if (rec.get("citation_count") or 0) > 1000:
        score += 1.0
    return score


def detect_plot_types(text: str) -> list[str]:
    txt = (text or "").lower()
    found = []
    for p in PLOT_VOCAB:
        if p in txt:
            found.append(p)
    return sorted(set(found))


def write_manifest(records: list[dict], path: Path) -> None:
    if not records:
        path.write_text("source,title,journal,year,doi,pmid,url,abstract\n")
        return
    cols = ["source", "title", "journal", "year", "authors", "doi", "pmid",
            "pmcid", "ss_id", "arxiv_id", "openalex_id", "citation_count",
            "influential_count", "url", "tldr", "abstract"]
    rows = []
    for r in records:
        rows.append({c: r.get(c, "") for c in cols})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(rows)


def trim_abstract(s: str, max_chars: int = 350) -> str:
    s = (s or "").strip()
    if len(s) <= max_chars:
        return s
    return s[: max_chars - 1].rstrip() + "…"


# ------------------------------- Subcommands ----------------------------------

def cmd_search(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    sources = args.sources or ["pubmed", "biorxiv", "arxiv",
                                "semanticscholar", "openalex"]
    all_recs: list[dict] = []
    if "pubmed" in sources:
        all_recs += pubmed_search(args.query, limit=args.limit,
                                  date_from=args.date_from,
                                  date_to=args.date_to,
                                  journals=args.journals or None)
    if "biorxiv" in sources:
        all_recs += biorxiv_search(args.query, limit=args.limit, server="biorxiv")
    if "medrxiv" in sources:
        all_recs += biorxiv_search(args.query, limit=args.limit, server="medrxiv")
    if "arxiv" in sources:
        all_recs += arxiv_search(args.query, limit=args.limit,
                                 categories=args.arxiv_categories or None)
    if "semanticscholar" in sources:
        all_recs += semanticscholar_search(args.query, limit=args.limit)
    if "openalex" in sources:
        all_recs += openalex_search(args.query, limit=args.limit)

    deduped = dedup_records(all_recs)
    terms = [t for t in re.split(r"[\s,]+", args.query) if len(t) > 2]
    deduped.sort(key=lambda r: score_relevance(r, terms), reverse=True)

    label = (args.label or args.query).replace(" ", "_")[:60]
    out_dir = REPORT_DIR / f"{ts}_search_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = MANIFEST_DIR / f"{ts}_search_{label}_manifest.csv"
    write_manifest(deduped, manifest)

    report = out_dir / "search_report.md"
    lines = [f"# Reference Search: `{args.query}`", "",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "",
             f"Records (deduplicated): **{len(deduped)}**", "",
             f"Manifest: `{manifest.relative_to(ROOT)}`", "",
             "## Top results", ""]
    for i, r in enumerate(deduped[: args.top], 1):
        plots = detect_plot_types(r.get("abstract", ""))
        cite = r.get("citation_count") or ""
        lines += [
            f"### {i}. {r.get('title','(untitled)')}",
            f"- **{r.get('journal','')}** · {r.get('year','')} · "
            f"sources: `{r.get('source','')}` · citations: {cite}",
            f"- Authors: {r.get('authors','')}",
            (f"- DOI: [{r['doi']}](https://doi.org/{r['doi']})"
             if r.get('doi') else ""),
            (f"- URL: <{r['url']}>" if r.get("url") else ""),
            (f"- Plot types mentioned: {', '.join(plots)}" if plots else ""),
            "",
            f"  > {trim_abstract(r.get('abstract',''))}",
            "",
        ]
    report.write_text("\n".join(lines))
    logging.info("Wrote: %s", report)
    print(f"Report:   {report}")
    print(f"Manifest: {manifest}")
    return report


def cmd_learn(args: argparse.Namespace) -> Path:
    """Find similar studies and extract methods + visualization patterns."""
    setup_logging(); mkdirs()
    ts = timestamp()

    journals = args.journals or sorted(TOP_JOURNALS)
    pubs = pubmed_search(args.topic, limit=args.limit, journals=journals)
    bio = biorxiv_search(args.topic, limit=args.limit)
    ax = arxiv_search(args.topic, limit=args.limit,
                      categories=["q-bio.GN", "q-bio.QM"])
    ss = semanticscholar_search(args.topic, limit=args.limit)
    deduped = dedup_records(pubs + bio + ax + ss)
    terms = [t for t in re.split(r"[\s,]+", args.topic) if len(t) > 2]
    deduped.sort(key=lambda r: score_relevance(r, terms), reverse=True)
    deduped = deduped[: args.top]

    plot_counts: dict[str, int] = defaultdict(int)
    method_phrases: list[tuple[str, str]] = []
    method_re = re.compile(
        r"\b(?:we (?:used|applied|performed|analy[sz]ed|integrated|identified)"
        r"|using\s+\w+|based on\s+\w+)[^.]{4,180}\.", re.I,
    )
    for r in deduped:
        for p in detect_plot_types(r.get("abstract", "")):
            plot_counts[p] += 1
        for m in method_re.findall(r.get("abstract", "") or ""):
            method_phrases.append((r.get("title", ""), m.strip()))

    label = (args.label or args.topic).replace(" ", "_")[:60]
    out_dir = REPORT_DIR / f"{ts}_learn_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = MANIFEST_DIR / f"{ts}_learn_{label}_manifest.csv"
    write_manifest(deduped, manifest)

    plot_table = "| Plot type | # papers |\n|---|---:|\n" + "\n".join(
        f"| {p} | {n} |" for p, n in sorted(plot_counts.items(),
                                              key=lambda x: -x[1])
    ) if plot_counts else "_No recognizable plot vocabulary in retrieved abstracts._"

    method_lines = []
    for title, phrase in method_phrases[: args.method_examples]:
        method_lines.append(f"- _{title[:90]}_: {phrase}")
    if not method_lines:
        method_lines = ["_No method-style sentences extracted from abstracts._"]

    report = out_dir / "learn_report.md"
    lines = [
        f"# Learn From The Field: `{args.topic}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Top-ranked papers retained: **{len(deduped)}**  ·  "
        f"Manifest: `{manifest.relative_to(ROOT)}`",
        "",
        "## Most-cited / top-journal hits",
        "",
    ]
    for i, r in enumerate(deduped[:10], 1):
        lines += [
            f"### {i}. {r.get('title','(untitled)')}",
            f"- **{r.get('journal','')}** · {r.get('year','')} · "
            f"citations: {r.get('citation_count','')}",
            (f"- DOI: https://doi.org/{r['doi']}" if r.get("doi") else ""),
            f"- TL;DR: {r.get('tldr') or trim_abstract(r.get('abstract',''))}",
            "",
        ]
    lines += ["## Visualization patterns observed", "", plot_table, "",
              "## Methodology phrasing observed", "", *method_lines, "",
              "## Suggested figure recipe (consensus)", ""]
    consensus = [p for p, n in plot_counts.items() if n >= 2][:6]
    if consensus:
        lines.append("Based on the retrieved literature, a publication-style "
                     "figure pack for this topic should typically include:")
        for p in consensus:
            lines.append(f"- **{p}** plot — appears in ≥2 retrieved papers.")
    else:
        lines.append("Not enough abstracts had explicit plot-type vocabulary "
                     "to derive a consensus. Default panel: UMAP, dot-plot, "
                     "violin, heatmap, volcano.")
    report.write_text("\n".join(lines))
    logging.info("Wrote: %s", report)
    print(f"Report:   {report}")
    print(f"Manifest: {manifest}")
    return report


def cmd_validate(args: argparse.Namespace) -> Path:
    """Cross-check a list of genes/variants/regions against prior literature."""
    setup_logging(); mkdirs()
    ts = timestamp()

    items = read_validate_input(args.input)
    if not items:
        raise SystemExit(f"No items found in {args.input}")
    context = " ".join(args.context or [])

    label = (args.label or Path(args.input).stem).replace(" ", "_")[:60]
    out_dir = REPORT_DIR / f"{ts}_validate_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = MANIFEST_DIR / f"{ts}_validate_{label}_manifest.csv"

    rows: list[dict] = []
    item_summaries: list[dict] = []
    for it in items[: args.max_items]:
        q = it["term"]
        if context:
            q = f"{q} AND ({context})"
        recs = pubmed_search(q, limit=args.limit_per_item)
        recs += semanticscholar_search(q, limit=args.limit_per_item)
        recs = dedup_records(recs)
        terms = [it["term"]] + [c for c in (args.context or [])]
        recs.sort(key=lambda r: score_relevance(r, terms), reverse=True)
        for r in recs[: args.limit_per_item]:
            row = {"query_term": it["term"], "category": it.get("category", ""),
                   **r}
            rows.append(row)
        item_summaries.append({
            "term": it["term"],
            "category": it.get("category", ""),
            "n_hits": len(recs),
            "top_papers": [(r.get("title", "")[:120], r.get("doi", ""),
                            r.get("year", ""))
                           for r in recs[:3]],
        })

    write_manifest(rows, manifest)

    lines = [f"# Validate Discoveries Against Literature",
             "",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
             "",
             f"Input: `{args.input}`  ·  Items checked: {len(item_summaries)}",
             (f"Context: `{context}`" if context else ""),
             f"Manifest: `{manifest.relative_to(ROOT)}`",
             "",
             "## Summary",
             "",
             "| Term | Category | Prior hits | Top paper |",
             "|---|---|---:|---|"]
    for s in item_summaries:
        top = s["top_papers"][0] if s["top_papers"] else ("(none)", "", "")
        lines.append(
            f"| `{s['term']}` | {s['category']} | {s['n_hits']} | "
            f"{top[0]}{' ('+top[2]+')' if top[2] else ''} |"
        )
    lines += ["", "## Per-item details", ""]
    for s in item_summaries:
        lines += [f"### {s['term']}  _({s['category'] or 'unspecified'})_",
                  f"Prior literature hits: **{s['n_hits']}**", ""]
        if not s["top_papers"]:
            lines.append("_No prior literature found at the configured "
                          "search limits._\n")
            continue
        for title, doi, year in s["top_papers"]:
            lines.append(
                f"- {title}{' ('+year+')' if year else ''}"
                f"{' · https://doi.org/'+doi if doi else ''}"
            )
        lines.append("")
    report = out_dir / "validate_report.md"
    report.write_text("\n".join(lines))
    logging.info("Wrote: %s", report)
    print(f"Report:   {report}")
    print(f"Manifest: {manifest}")
    return report


def read_validate_input(path: str) -> list[dict]:
    p = Path(path)
    if not p.exists():
        return []
    if p.suffix.lower() in (".csv", ".tsv"):
        sep = "," if p.suffix.lower() == ".csv" else "\t"
        out = []
        with p.open() as f:
            rdr = csv.DictReader(f, delimiter=sep)
            for row in rdr:
                term = (row.get("term") or row.get("gene") or row.get("variant")
                        or row.get("rsid") or row.get("element")
                        or row.get("name") or "").strip()
                if not term:
                    continue
                cat = (row.get("category") or row.get("type")
                       or guess_category(term))
                out.append({"term": term, "category": cat})
        return out
    return [{"term": t.strip(), "category": guess_category(t.strip())}
            for t in p.read_text().splitlines() if t.strip()]


def guess_category(term: str) -> str:
    if re.match(r"^rs\d+$", term, re.I):
        return "variant"
    if re.match(r"^chr[0-9XYM]+:\d+", term):
        return "region"
    if re.match(r"^[A-Z][A-Z0-9-]{1,12}$", term):
        return "gene"
    return ""


# --------- Curated study-design workflows for the `design` subcommand --------

WORKFLOWS = {
    "10x_multiome": {
        "label": "10x Multiome (snATAC + snRNA, paired)",
        "phases": [
            "Discover IGVF AnalysisSets with preferred_assay_titles=10x multiome.",
            "Per-sample QC: TSS enrichment, fragment count, %mito, gene complexity.",
            "Doublet removal (e.g. AMULET on ATAC, scDblFinder on RNA).",
            "Normalization (CP10k log1p for RNA; TF-IDF for ATAC).",
            "Joint embedding: WNN (Seurat) or MultiVI / GLUE.",
            "Cell-type annotation (reference mapping or marker-based).",
            "Differential accessibility / expression by cell type and condition.",
            "Peak-to-gene linking (Signac LinkPeaks, ArchR addPeak2GeneLinks, ABC).",
            "Variant overlap with cell-type-resolved cCREs.",
        ],
        "plots": ["UMAP by cell type", "ATAC TSS-enrichment violin",
                   "marker gene dot plot", "peak track plot",
                   "peak-to-gene heatmap", "differential MA plot"],
        "queries": ["10x multiome single-cell", "snATAC snRNA paired",
                     "WNN multiomic clustering"],
    },
    "parse_split_seq": {
        "label": "Parse SPLiT-seq (combinatorial barcoding snRNA-seq)",
        "phases": [
            "Demultiplex by donor genotype (souporcell, vireo, demuxlet).",
            "Per-cell QC and filtering (genes detected, %mito, doublet score).",
            "Normalize, log1p, HVG, scale, PCA.",
            "Batch / pool integration (Harmony, scVI) — important for SPLiT-seq runs.",
            "UMAP, Leiden clustering.",
            "Cell-type annotation (Azimuth, CellTypist, marker-based).",
            "Per-strain / per-condition differential expression.",
            "If multiplexed across genotypes/strains: eQTL or strain-effect modelling.",
        ],
        "plots": ["UMAP by cell type / strain / pool", "marker dot plot",
                   "stacked composition bar", "DEG volcano",
                   "strain-effect heatmap"],
        "queries": ["Parse SPLiT-seq", "combinatorial barcoding single-nucleus",
                     "split-pool ligation transcriptomics"],
    },
    "mpra": {
        "label": "MPRA / STARR / BlueSTARR",
        "phases": [
            "Library design and tile mapping to genomic coordinates.",
            "Read-count normalization (DNA vs RNA, replicate-aware).",
            "Activity calling (MPRAnalyze, mpralm, BeanCounter).",
            "QC: replicate correlation, barcode complexity, scrambled controls.",
            "Variant-effect estimation for allelic constructs.",
            "Annotation against cCRE / TF motifs / Catalog evidence.",
        ],
        "plots": ["replicate scatter", "allelic effect volcano",
                   "tile activity heatmap", "TF motif enrichment"],
        "queries": ["massively parallel reporter assay variants",
                     "saturation mutagenesis MPRA enhancer",
                     "BlueSTARR STARR-seq"],
    },
    "crispri": {
        "label": "CRISPRi / CRISPR-FACS / Perturb-seq screens",
        "phases": [
            "Library design and guide-to-element mapping.",
            "Read-count quantification and QC (guide diversity, replicate "
            "correlation).",
            "Hit calling (MAGeCK, BEAN, CASA, scMAGeCK for single-cell).",
            "Effect annotation against cCRE classes and gene-element linkage.",
            "Cell-type specificity if Perturb-seq.",
            "Cross-modality validation with MPRA, ATAC, GWAS.",
        ],
        "plots": ["guide rank plot", "MAGeCK volcano", "Perturb-seq UMAP",
                   "element-by-gene effect heatmap"],
        "queries": ["CRISPRi enhancer screen", "Perturb-seq",
                     "non-coding variant CRISPR FACS"],
    },
    "enhancer_gene": {
        "label": "Enhancer-gene linkage",
        "phases": [
            "Discover available linkage methods: ABC, rE2G, eQTL, "
            "Hi-C/Promoter-Capture, scE2G.",
            "Pull per-cell-type predictions and harmonize coordinates.",
            "Compare methods (ABC vs rE2G vs eQTL agreement).",
            "Annotate variants of interest with linked genes per cell type.",
            "Validate against experimental CRISPRi / MPRA evidence.",
        ],
        "plots": ["upset diagram of method overlap", "linkage score heatmap",
                   "track plot with peak-gene loops", "tissue-by-method bar"],
        "queries": ["enhancer gene linkage activity by contact",
                     "rE2G prediction enhancer", "eQTL enhancer mapping"],
    },
}


def cmd_design(args: argparse.Namespace) -> Path:
    """Recommend a workflow + cognate published studies + IGVF data."""
    setup_logging(); mkdirs()
    ts = timestamp()

    key = args.data_type.lower().replace("-", "_").replace(" ", "_")
    if key not in WORKFLOWS:
        # fuzzy match
        for k in WORKFLOWS:
            if k in key or key in k:
                key = k; break
    if key not in WORKFLOWS:
        raise SystemExit(
            f"Unknown data type: {args.data_type!r}. Known: "
            + ", ".join(WORKFLOWS.keys())
        )
    wf = WORKFLOWS[key]
    label = (args.label or key)[:60]
    out_dir = REPORT_DIR / f"{ts}_design_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_query = " OR ".join(f"\"{q}\"" for q in wf["queries"])
    ref_recs = pubmed_search(ref_query, limit=args.lit_limit,
                              journals=sorted(TOP_JOURNALS))
    ref_recs += semanticscholar_search(wf["queries"][0], limit=args.lit_limit)
    ref_recs = dedup_records(ref_recs)
    ref_recs.sort(key=lambda r: score_relevance(r, wf["queries"]), reverse=True)
    ref_manifest = MANIFEST_DIR / f"{ts}_design_{label}_lit_manifest.csv"
    write_manifest(ref_recs[:50], ref_manifest)

    igvf_hits = igvf_portal_query(args.assay_title or args.data_type,
                                   limit=args.igvf_limit)
    igvf_path = MANIFEST_DIR / f"{ts}_design_{label}_igvf_manifest.csv"
    write_igvf_manifest(igvf_hits, igvf_path)

    lines = [f"# Workflow Design: {wf['label']}",
             "",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
             "",
             "## Recommended workflow phases",
             ""]
    for i, phase in enumerate(wf["phases"], 1):
        lines.append(f"{i}. {phase}")
    lines += ["", "## Plots and tables to produce", ""]
    for p in wf["plots"]:
        lines.append(f"- {p}")
    lines += ["",
               "## Cognate published studies",
               "",
               f"Manifest: `{ref_manifest.relative_to(ROOT)}` (top {min(50, len(ref_recs))} of {len(ref_recs)})",
               ""]
    for i, r in enumerate(ref_recs[: args.top], 1):
        lines += [
            f"### {i}. {r.get('title','(untitled)')}",
            f"- **{r.get('journal','')}** · {r.get('year','')} · "
            f"citations: {r.get('citation_count','')}",
            (f"- https://doi.org/{r['doi']}" if r.get("doi") else ""),
            f"- {r.get('tldr') or trim_abstract(r.get('abstract',''))}",
            "",
        ]
    lines += ["",
               "## Matching IGVF Portal datasets",
               "",
               f"Manifest: `{igvf_path.relative_to(ROOT)}`",
               "",
               "| Accession | Type | Assay | Lab | Status | Description |",
               "|---|---|---|---|---|---|"]
    for h in igvf_hits[:25]:
        desc = (h.get("description") or h.get("summary") or "")[:120]
        lines.append(
            f"| {h.get('accession','')} | {h.get('@type','')} | "
            f"{h.get('preferred_assay_titles','') or h.get('assay_titles','')} | "
            f"{h.get('lab','')} | {h.get('status','')} | {desc} |"
        )
    if not igvf_hits:
        lines.append("| _no IGVF Portal results — try a different assay_title_ "
                      "| | | | | |")

    report = out_dir / "design_report.md"
    report.write_text("\n".join(lines))
    logging.info("Wrote: %s", report)
    print(f"Report:   {report}")
    print(f"Lit manifest:  {ref_manifest}")
    print(f"IGVF manifest: {igvf_path}")
    return report


def igvf_portal_query(assay_title: str, limit: int = 25) -> list[dict]:
    base = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
    params = {"type": "AnalysisSet", "format": "json", "limit": limit,
              "preferred_assay_titles": assay_title}
    status, body = cached_get("crossref",  # generic pool
                              f"{base}/search/", params=params)
    if status == 200 and body:
        try:
            data = json.loads(body)
            graph = data.get("@graph", [])
            simplified = []
            for g in graph[:limit]:
                lab = (g.get("lab") or {})
                simplified.append({
                    "@type": (g.get("@type") or [""])[0],
                    "accession": g.get("accession", ""),
                    "preferred_assay_titles":
                        ", ".join(g.get("preferred_assay_titles", []) or []),
                    "assay_titles":
                        ", ".join(g.get("assay_titles", []) or []),
                    "lab": (lab.get("title") or lab.get("@id", ""))
                              if isinstance(lab, dict) else str(lab),
                    "status": g.get("status", ""),
                    "description": g.get("description", ""),
                    "summary": g.get("summary", ""),
                })
            return simplified
        except Exception as e:
            logging.warning("Portal parse failed: %s", e)
    # Fallback: no results
    return []


def write_igvf_manifest(hits: list[dict], path: Path) -> None:
    cols = ["accession", "@type", "preferred_assay_titles", "assay_titles",
            "lab", "status", "summary", "description"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for h in hits:
            w.writerow({c: h.get(c, "") for c in cols})


# ----------------------------- write-playbook ---------------------------------

def cmd_write_playbook(args: argparse.Namespace) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "REFERENCE_SKILLS.md"
    lines = [
        "# Skill: Reference (literature retrieval, validation, study design)",
        "",
        "Use this skill when an IGVF agent run needs to consult prior "
        "literature: to scope a problem, to validate discoveries, or to "
        "design a study workflow that mirrors successful published practice.",
        "",
        "## Sources",
        "",
        "- PubMed / PMC (NCBI E-utilities)",
        "- bioRxiv + medRxiv (via Crossref)",
        "- arXiv (Atom API)",
        "- Semantic Scholar (paper search + recommendation)",
        "- OpenAlex (open scholarly graph)",
        "",
        "Top-tier journals are weighted in ranking: Nature/Cell/Science "
        "families, NEJM, Nucleic Acids Research, Bioinformatics, Genome "
        "Biology / Genome Research, eLife, PNAS.",
        "",
        "## Subcommands",
        "",
        "### 1. `learn` — what does the field do?",
        "",
        "```bash",
        "python3 Scripts/reference_skill.py learn \\",
        "    --topic '10x multiome human putamen' --limit 30 --top 15",
        "```",
        "",
        "Returns a manifest of top-ranked papers and a report extracting:",
        "- recurring methodology phrasing,",
        "- visualization vocabulary used in abstracts (UMAP / dot-plot / "
        "volcano / track / heatmap …),",
        "- a consensus figure recipe.",
        "",
        "### 2. `validate` — has anyone seen this before?",
        "",
        "```bash",
        "python3 Scripts/reference_skill.py validate \\",
        "    --input Docs/AdvancedVariantAnalysis/<run>/<label>_summary_stats.csv \\",
        "    --context 'putamen Parkinson' --limit-per-item 5",
        "```",
        "",
        "Input is any CSV with `gene` / `variant` / `rsid` / `element` / "
        "`term` columns (case-insensitive). Each row is searched against "
        "PubMed and Semantic Scholar; the manifest links each discovery to "
        "the prior-evidence papers, and the report ranks the strength of "
        "prior support.",
        "",
        "### 3. `design` — what should our study look like?",
        "",
        "```bash",
        "python3 Scripts/reference_skill.py design \\",
        "    --data-type parse_split_seq --assay-title 'Parse SPLiT-seq'",
        "```",
        "",
        "Known data types: " + ", ".join(f"`{k}`" for k in WORKFLOWS) + ".",
        "",
        "Outputs: a recommended pipeline (phases + plots/tables), a literature "
        "manifest of cognate published studies, and a manifest of IGVF Portal "
        "AnalysisSets that match the assay.",
        "",
        "### `search` — generic multi-source search",
        "",
        "```bash",
        "python3 Scripts/reference_skill.py search \\",
        "    --query 'enhancer-gene linkage rE2G ABC' --top 20",
        "```",
        "",
        "## Caching",
        "",
        "All API responses are cached under `Data/Cache/References/<source>/` "
        "with a 14-day TTL. Re-running the same query is free.",
        "",
        "## Outputs",
        "",
        "- Reports: `Docs/References/<timestamp>_<subcommand>_<label>/*.md`",
        "- Manifests: `Data/Manifests/References/`",
        "- Logs: `Docs/Logs/reference_skill_*.log`",
        "",
        "## How this skill chains with other IGVF agent skills",
        "",
        "1. After `data_illustration_interpretation.py` produces a dataset "
        "summary, call `learn` with the assay/biosample to get the typical "
        "analysis recipe.",
        "2. After `advanced_variant_analysis.py` writes a discovery table, "
        "feed the `summary_stats.csv` (or any gene/variant CSV) into "
        "`validate`.",
        "3. Before starting a new study, call `design` with the planned IGVF "
        "data type to seed the workflow and surface matching IGVF datasets.",
    ]
    path.write_text("\n".join(lines))
    logging.info("Wrote: %s", path)
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI ----------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reference skill: retrieve, validate, and design from "
                    "the published literature."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search",
                        help="Multi-source keyword search.")
    s.add_argument("--query", required=True)
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--top", type=int, default=20)
    s.add_argument("--label", default="")
    s.add_argument("--sources", nargs="*",
                    choices=["pubmed", "biorxiv", "medrxiv", "arxiv",
                             "semanticscholar", "openalex"])
    s.add_argument("--journals", nargs="*", default=None)
    s.add_argument("--date-from", default=None)
    s.add_argument("--date-to", default=None)
    s.add_argument("--arxiv-categories", nargs="*", default=None)
    s.set_defaults(func=cmd_search)

    s = sub.add_parser("learn",
                        help="Find what other studies do for a topic / assay.")
    s.add_argument("--topic", required=True)
    s.add_argument("--limit", type=int, default=30)
    s.add_argument("--top", type=int, default=15)
    s.add_argument("--method-examples", type=int, default=20)
    s.add_argument("--journals", nargs="*", default=None)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_learn)

    s = sub.add_parser("validate",
                        help="Check IGVF agent discoveries against the "
                             "published literature.")
    s.add_argument("--input", required=True,
                    help="CSV/TSV/text file of genes/variants/elements.")
    s.add_argument("--context", nargs="*", default=None,
                    help="Free-text context terms (e.g. tissue, disease).")
    s.add_argument("--limit-per-item", type=int, default=5)
    s.add_argument("--max-items", type=int, default=50)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_validate)

    s = sub.add_parser("design",
                        help="Workflow recommendation for an IGVF data type.")
    s.add_argument("--data-type", required=True,
                    help="One of: " + ", ".join(WORKFLOWS.keys()))
    s.add_argument("--assay-title", default=None,
                    help="IGVF Portal preferred_assay_titles to filter by.")
    s.add_argument("--lit-limit", type=int, default=25)
    s.add_argument("--igvf-limit", type=int, default=25)
    s.add_argument("--top", type=int, default=10)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_design)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/REFERENCE_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
