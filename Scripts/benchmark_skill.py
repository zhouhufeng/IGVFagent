#!/usr/bin/env python3
"""Benchmark Skill: publication -> reproduction scaffold -> run -> concordance.

This skill automates what was previously hand-built for each of the papers
under ``Benchmarks/``: take whatever the user knows about a publication
(title, URL, DOI, PubMed ID, authors, journal, year), pin down *which* paper
that is, read the paper's own Data Availability statement, work out which
IGVFagent analytical chain can reproduce it, and write a runnable benchmark
directory that the existing ``Benchmarks/concordance.py`` scorer understands.

Seven subcommands, meant to be run in order (or all at once via ``pipeline``):

  1. ``resolve``     — any identifier or free text -> ONE canonical paper
                       record (Crossref + PubMed + Europe PMC + bioRxiv +
                       OpenAlex + Semantic Scholar), with a confidence score
                       and explicit disambiguation candidates when the match
                       is not clean.

  2. ``harvest``     — pull the paper's full text (Europe PMC full-text XML,
                       covering both PMC open-access and bioRxiv/medRxiv
                       preprints), extract the Data / Code Availability
                       statements, every repository accession it mentions,
                       the assay families it used, and candidate numeric
                       ground-truth claims.

  3. ``route``       — map accessions + assay families onto an IGVFagent
                       skill chain, using a routing table seeded from the
                       benchmarks already in this repository.

  4. ``scaffold``    — write ``Benchmarks/<paper-id>/`` with README.md,
                       OPERATIONS.md, expected.json, run.sh and
                       provenance.json, and register the paper in
                       ``Benchmarks/generated.txt``.

  5. ``run``         — execute the generated ``run.sh``.

  6. ``score``       — delegate to ``Benchmarks/concordance.py``.

  7. ``report``      — render the paper-claim vs IGVFagent-measured
                       comparison.

Plus ``pipeline`` (1-4, optionally 5-7 with ``--execute``), ``selftest``
(re-derive the already-committed benchmarks from their DOIs and score the
resolver + router against them), and ``list-routes``.

Two deliberate honesty constraints are built in:

* **Ground truth extracted from prose is never silently trusted.** Every
  check this skill proposes from the paper's text is written into
  ``expected.json`` with ``"confirmed": false`` and a ``provenance`` block
  quoting the sentence it came from. ``concordance.py`` tallies unconfirmed
  checks separately and never lets them turn a run green. A human promotes
  a check by setting ``"confirmed": true``.

* **Routes that cannot actually reproduce the paper say so.** When the data
  is controlled-access, embargoed, or in a format IGVFagent cannot read, the
  generated ``run.sh`` exits 77 (the suite's "skipped — missing local input"
  convention) with download instructions, rather than pretending.

Outputs follow the project pattern: cached HTTP under
``Data/Cache/References/<source>/``, per-run artefacts under
``Docs/Benchmark/<timestamp>_<paper-id>/``, logs in ``Docs/Logs/``.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

# Reuse the reference skill's cached, rate-limited HTTP layer and its
# multi-source search clients rather than duplicating them.
import reference_skill as _ref  # noqa: E402


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "Benchmark"
BENCH_DIR = ROOT / "Benchmarks"
GENERATED_LIST = BENCH_DIR / "generated.txt"

# Conservative pacing, matching reference_skill's convention.
_ref.SOURCE_DELAY.setdefault("europepmc", 0.34)
_ref.SOURCE_DELAY.setdefault("crossref", 0.2)


# ----------------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / (
        f"benchmark_skill_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def _norm_title(s: str) -> str:
    """Aggressively normalise a title for comparison."""
    s = re.sub(r"<[^>]+>", " ", s or "")
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


_STOPWORDS = {
    "a", "an", "and", "are", "as", "at", "be", "by", "for", "from", "in",
    "into", "is", "its", "of", "on", "or", "the", "to", "with", "using",
    "via", "that", "this", "we", "their", "human", "cell", "cells", "gene",
    "genes", "analysis", "study", "high", "novel", "new", "reveals",
    "identifies", "across", "during", "between", "single",
}


def _title_tokens(s: str) -> set:
    return {t for t in _norm_title(s).split() if t not in _STOPWORDS and len(t) > 2}


def _jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / float(len(a | b))


def _surname(authors: str) -> str:
    """First author surname from a free-text author string."""
    if not authors:
        return ""
    first = re.split(r"[;,]", authors.strip())[0].strip()
    parts = [p for p in first.split() if p]
    if not parts:
        return ""
    # "Buckley M" -> Buckley;  "M Buckley" -> Buckley
    if len(parts) >= 2 and len(parts[-1]) <= 3 and parts[-1].isupper():
        cand = parts[0]
    else:
        cand = parts[-1] if len(parts[-1]) > 3 else parts[0]
    return re.sub(r"[^A-Za-z]", "", cand)


def _write_json(path: Path, obj: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))
    return path


def _read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def _runs_for(paper_id: str) -> List[Path]:
    """All run dirs for a paper, newest first."""
    if not REPORT_DIR.is_dir():
        return []
    return sorted(
        (p for p in REPORT_DIR.glob(f"2*_{paper_id}") if p.is_dir()),
        key=lambda p: p.name, reverse=True,
    )


def _latest_run(paper_id: str) -> Optional[Path]:
    runs = _runs_for(paper_id)
    return runs[0] if runs else None


def _load_stage(paper_id: str, stage: str, explicit: Optional[str] = None) -> Any:
    """Load ``resolution.json`` / ``harvest.json`` / ``routing.json``.

    Looks at an explicit path first, then the benchmark dir's
    ``provenance.json``, then the newest ``Docs/Benchmark/<ts>_<paper-id>/``.
    """
    if explicit:
        return _read_json(Path(explicit))
    # Newest run dir wins. Each subcommand writes its own timestamped dir, so
    # the newest dir is not necessarily the one holding this stage — scan back
    # through them. Run dirs are checked before the benchmark's provenance.json
    # so that re-running `harvest` and re-scaffolding picks up the new data
    # instead of the copy frozen into the previous scaffold.
    for run in _runs_for(paper_id):
        if (run / f"{stage}.json").is_file():
            return _read_json(run / f"{stage}.json")
    prov = BENCH_DIR / paper_id / "provenance.json"
    if prov.is_file():
        blob = _read_json(prov)
        if stage in blob:
            return blob[stage]
    raise SystemExit(
        f"No {stage}.json found for paper-id {paper_id!r}. "
        f"Run `igvfagent bench {'resolve' if stage == 'resolution' else stage}` first, "
        f"or pass --{stage}-json <path>."
    )


# ----------------------------------------------------------------------------
# 1. Identifier classification
# ----------------------------------------------------------------------------

DOI_RE = re.compile(r"\b(10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+)")
PMID_RE = re.compile(r"^(?:pmid[:\s]*)?(\d{7,9})$", re.I)
PMCID_RE = re.compile(r"\b(PMC\d{6,9})\b", re.I)
ARXIV_RE = re.compile(r"\b(\d{4}\.\d{4,5})(v\d+)?\b")

# Publisher URL slug -> DOI prefix. Enough to turn a pasted article URL into
# a DOI without an extra network round-trip.
_URL_DOI_PATTERNS = [
    # https://www.nature.com/articles/s41588-024-01800-z
    (re.compile(r"nature\.com/articles/([a-z0-9\-]+)", re.I), "10.1038/{0}"),
    # https://www.biorxiv.org/content/10.1101/2023.11.09.563812v1
    (re.compile(r"(?:bio|med)rxiv\.org/content/(10\.\d{4,9}/[^\s?#]+?)(?:v\d+)?(?:\.full|\.pdf)?$", re.I), "{0}"),
    # https://doi.org/10.1126/science.adh0559
    (re.compile(r"doi\.org/(10\.\d{4,9}/[^\s?#]+)", re.I), "{0}"),
]

_URL_ID_PATTERNS = [
    (re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d{7,9})", re.I), "pmid"),
    (re.compile(r"ncbi\.nlm\.nih\.gov/pmc/articles/(PMC\d+)", re.I), "pmcid"),
    (re.compile(r"europepmc\.org/(?:article|abstract)/[a-z]+/(\d{7,9})", re.I), "pmid"),
    (re.compile(r"arxiv\.org/abs/(\d{4}\.\d{4,5})", re.I), "arxiv"),
]


def classify_identifier(raw: str) -> Dict[str, Any]:
    """Classify a user-supplied string into a typed identifier.

    Returns ``{"kind": ..., "value": ..., "raw": ...}`` where kind is one of
    ``doi`` / ``pmid`` / ``pmcid`` / ``arxiv`` / ``free_text``.
    """
    s = (raw or "").strip().strip("<>").rstrip(".,;")
    if not s:
        return {"kind": "free_text", "value": "", "raw": raw}

    if s.lower().startswith(("http://", "https://", "www.")):
        for pat, kind in _URL_ID_PATTERNS:
            m = pat.search(s)
            if m:
                return {"kind": kind, "value": m.group(1), "raw": raw,
                        "note": "extracted from URL"}
        for pat, tmpl in _URL_DOI_PATTERNS:
            m = pat.search(s)
            if m:
                doi = tmpl.format(*m.groups())
                doi = re.sub(r"(v\d+)?(\.full|\.pdf)?$", "", doi)
                return {"kind": "doi", "value": doi, "raw": raw,
                        "note": "extracted from URL"}
        m = DOI_RE.search(s)
        if m:
            return {"kind": "doi", "value": m.group(1).rstrip("."), "raw": raw,
                    "note": "DOI found inside URL"}
        return {"kind": "free_text", "value": s, "raw": raw,
                "note": "URL not recognised — falling back to text search"}

    m = PMCID_RE.search(s)
    if m and len(s) <= 15:
        return {"kind": "pmcid", "value": m.group(1).upper(), "raw": raw}
    m = DOI_RE.search(s)
    if m:
        return {"kind": "doi", "value": m.group(1).rstrip("."), "raw": raw}
    m = PMID_RE.match(s)
    if m:
        return {"kind": "pmid", "value": m.group(1), "raw": raw}
    if re.match(r"^arxiv[:\s]", s, re.I):
        m = ARXIV_RE.search(s)
        if m:
            return {"kind": "arxiv", "value": m.group(1), "raw": raw}
    return {"kind": "free_text", "value": s, "raw": raw}


# ----------------------------------------------------------------------------
# 2. Metadata clients
# ----------------------------------------------------------------------------

def crossref_work(doi: str) -> Optional[Dict[str, Any]]:
    """Canonical publisher metadata for a DOI."""
    base = _resolve_endpoint("crossref", "CROSSREF_BASE")
    from urllib.parse import quote
    status, body = _ref.cached_get("crossref", f"{base}/works/{quote(doi, safe='')}")
    if status != 200 or not body:
        return None
    try:
        msg = json.loads(body).get("message") or {}
    except Exception:
        return None
    if not msg:
        return None
    authors = []
    first_family = ""
    for a in msg.get("author") or []:
        nm = " ".join(x for x in [a.get("family"), a.get("given")] if x)
        if nm:
            authors.append(nm)
        if not first_family and a.get("family"):
            first_family = a["family"]
    date = (msg.get("published-print") or msg.get("published-online")
            or msg.get("issued") or {})
    parts = (date.get("date-parts") or [[None]])[0]
    year = parts[0] if parts else None
    titles = msg.get("title") or []
    containers = msg.get("container-title") or []
    return {
        "source": "crossref",
        "doi": (msg.get("DOI") or doi).lower(),
        "title": titles[0] if titles else "",
        "authors": "; ".join(authors),
        # Crossref gives family/given structurally — keep the surname so the
        # paper-id derivation never has to guess which token is the family name.
        "first_author_family": first_family,
        "journal": containers[0] if containers else "",
        "year": year,
        "type": msg.get("type"),
        "publisher": msg.get("publisher"),
        "volume": msg.get("volume"),
        "issue": msg.get("issue"),
        "page": msg.get("page"),
        "article_number": msg.get("article-number"),
        "abstract": re.sub(r"<[^>]+>", " ", msg.get("abstract") or "").strip(),
        "url": msg.get("URL"),
        "is_preprint": (msg.get("type") == "posted-content"),
        "relation": msg.get("relation") or {},
    }


def europepmc_query(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """Europe PMC search — the one index that covers PMC *and* preprints."""
    base = _resolve_endpoint("europepmc", "EUROPEPMC_BASE")
    status, body = _ref.cached_get(
        "europepmc", f"{base}/search",
        params={"query": query, "format": "json", "pageSize": limit,
                "resultType": "core"},
    )
    if status != 200 or not body:
        return []
    try:
        results = json.loads(body).get("resultList", {}).get("result", []) or []
    except Exception:
        return []
    out = []
    for r in results:
        out.append({
            "source": "europepmc",
            "title": r.get("title") or "",
            "authors": r.get("authorString") or "",
            "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title")
                        or r.get("bookOrReportDetails", {}).get("publisher")
                        or r.get("source") or "",
            "year": int(r["pubYear"]) if str(r.get("pubYear", "")).isdigit() else None,
            "doi": (r.get("doi") or "").lower(),
            "pmid": r.get("pmid") or "",
            "pmcid": r.get("pmcid") or "",
            "epmc_id": r.get("id") or "",
            "epmc_source": r.get("source") or "",
            "is_open_access": (r.get("isOpenAccess") == "Y"),
            "has_fulltext": (r.get("hasTextMinedTerms") == "Y"
                              or r.get("inEPMC") == "Y" or r.get("inPMC") == "Y"),
            "in_epmc": (r.get("inEPMC") == "Y"),
            "abstract": r.get("abstractText") or "",
            "citation_count": r.get("citedByCount"),
        })
    return out


def pubmed_by_id(pmid: str) -> Optional[Dict[str, Any]]:
    """esummary for one PMID, plus its PMCID if the article is in PMC."""
    hits = europepmc_query(f"EXT_ID:{pmid} AND SRC:MED", limit=3)
    for h in hits:
        if str(h.get("pmid")) == str(pmid):
            return h
    base = _resolve_endpoint("pubmed_eutils", "PUBMED_EUTILS_BASE")
    status, body = _ref.cached_get(
        "pubmed", f"{base}/esummary.fcgi",
        params={"db": "pubmed", "id": pmid, "retmode": "json"},
    )
    if status != 200 or not body:
        return None
    try:
        rec = json.loads(body).get("result", {}).get(str(pmid))
    except Exception:
        return None
    if not rec:
        return None
    doi = ""
    for aid in rec.get("articleids") or []:
        if aid.get("idtype") == "doi":
            doi = (aid.get("value") or "").lower()
    return {
        "source": "pubmed",
        "title": rec.get("title") or "",
        "authors": "; ".join(a.get("name", "") for a in (rec.get("authors") or [])),
        "journal": rec.get("fulljournalname") or rec.get("source") or "",
        "year": int(rec["pubdate"][:4]) if str(rec.get("pubdate", ""))[:4].isdigit() else None,
        "doi": doi,
        "pmid": str(pmid),
        "abstract": "",
    }


def biorxiv_details(doi: str) -> Optional[Dict[str, Any]]:
    """bioRxiv/medRxiv detail record — also reveals the published DOI."""
    base = _resolve_endpoint("biorxiv_api", "BIORXIV_API_BASE")
    for server in ("biorxiv", "medrxiv"):
        status, body = _ref.cached_get(
            "biorxiv", f"{base}/details/{server}/{doi}", accept_json=True)
        if status != 200 or not body:
            continue
        try:
            coll = json.loads(body).get("collection") or []
        except Exception:
            continue
        if not coll:
            continue
        rec = coll[-1]
        return {
            "source": server,
            "title": rec.get("title") or "",
            "authors": rec.get("authors") or "",
            "journal": server,
            "year": int(str(rec.get("date", ""))[:4]) if str(rec.get("date", ""))[:4].isdigit() else None,
            "doi": (rec.get("doi") or doi).lower(),
            "version": rec.get("version"),
            "abstract": rec.get("abstract") or "",
            "published_doi": (rec.get("published") or "").lower()
                              if rec.get("published") not in (None, "NA") else "",
        }
    return None


def _merge_records(*records: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Merge metadata records, first non-empty value wins per field."""
    out: Dict[str, Any] = {}
    seen_sources = []
    for rec in records:
        if not rec:
            continue
        seen_sources.append(rec.get("source"))
        for k, v in rec.items():
            if k == "source":
                continue
            if v in (None, "", [], {}):
                continue
            if k not in out or out[k] in (None, "", [], {}):
                out[k] = v
    out["sources"] = [s for s in seen_sources if s]
    return out


# ----------------------------------------------------------------------------
# 3. Candidate scoring / disambiguation
# ----------------------------------------------------------------------------

def score_candidate(rec: Dict[str, Any], want: Dict[str, Any]) -> Tuple[float, List[str]]:
    """Score a candidate against the user's stated constraints.

    Returns (score in [0, 1], list of human-readable reasons).
    """
    reasons: List[str] = []
    score = 0.0
    weight = 0.0

    if want.get("title"):
        w = 0.55
        weight += w
        a, b = _norm_title(rec.get("title", "")), _norm_title(want["title"])
        if a and b and (a == b):
            score += w
            reasons.append("title exact")
        else:
            j = _jaccard(_title_tokens(rec.get("title", "")), _title_tokens(want["title"]))
            score += w * j
            reasons.append(f"title overlap {j:.2f}")

    if want.get("author"):
        w = 0.2
        weight += w
        auth_l = (rec.get("authors") or "").lower()
        wants = [a.strip().lower() for a in re.split(r"[;,]", want["author"]) if a.strip()]
        hits = sum(1 for a in wants if a and a in auth_l)
        frac = hits / float(len(wants)) if wants else 0.0
        score += w * frac
        reasons.append(f"author match {hits}/{len(wants)}")

    if want.get("journal"):
        w = 0.15
        weight += w
        jr = _norm_title(rec.get("journal", ""))
        jw = _norm_title(want["journal"])
        if jr and jw and (jw in jr or jr in jw):
            score += w
            reasons.append("journal match")
        else:
            reasons.append(f"journal mismatch ({rec.get('journal','?')})")

    if want.get("year"):
        w = 0.1
        weight += w
        try:
            dy = abs(int(rec.get("year") or 0) - int(want["year"]))
        except (TypeError, ValueError):
            dy = 99
        if dy == 0:
            score += w
            reasons.append("year exact")
        elif dy == 1:
            score += w * 0.6
            reasons.append("year ±1 (preprint/issue lag)")
        else:
            reasons.append(f"year mismatch ({rec.get('year')})")

    if weight <= 0:
        return 0.0, ["no constraints to score against"]
    return score / weight, reasons


def derive_paper_id(rec: Dict[str, Any], assays: Optional[List[str]] = None) -> str:
    """``<surname><year>_<slug>`` — matches the existing directory convention."""
    surname = (re.sub(r"[^A-Za-z]", "", rec.get("first_author_family") or "")
                or _surname(rec.get("authors", "")) or "paper").lower()
    year = rec.get("year") or ""
    slug_src: List[str] = []
    for a in (assays or [])[:1]:
        slug_src.append(re.sub(r"[^a-z0-9]+", "", a.lower())[:14])
    for tok in _norm_title(rec.get("title", "")).split():
        if tok in _STOPWORDS or len(tok) < 4:
            continue
        if tok in slug_src:
            continue
        slug_src.append(tok)
        if len(slug_src) >= 2:
            break
    slug = "_".join(s for s in slug_src if s) or "benchmark"
    return f"{surname}{year}_{slug}"[:60]


# ----------------------------------------------------------------------------
# 4. Full-text harvest
# ----------------------------------------------------------------------------

# Repository accession patterns. Order matters only for readability; every
# pattern is applied. ``controlled`` marks repositories whose data cannot be
# fetched without an approved application — the scaffolder turns those into
# exit-77 guards rather than pretending the chain will run.
ACCESSION_PATTERNS: List[Dict[str, Any]] = [
    {"key": "geo_series",    "re": r"\bGSE\d{4,7}\b",                     "repo": "NCBI GEO"},
    {"key": "geo_sample",    "re": r"\bGSM\d{4,9}\b",                     "repo": "NCBI GEO"},
    {"key": "sra_project",   "re": r"\bSRP\d{5,9}\b",                     "repo": "NCBI SRA"},
    {"key": "bioproject",    "re": r"\bPRJ(?:NA|EB|DB)\d{4,9}\b",         "repo": "BioProject"},
    {"key": "igvf_dataset",  "re": r"\bIGVFDS[A-Z0-9]{6,}\b",             "repo": "IGVF Portal"},
    {"key": "igvf_file",     "re": r"\bIGVFFI[A-Z0-9]{6,}\b",             "repo": "IGVF Portal"},
    {"key": "encode_series", "re": r"\bENCSR[A-Z0-9]{6}\b",               "repo": "ENCODE"},
    {"key": "encode_file",   "re": r"\bENCFF[A-Z0-9]{6}\b",               "repo": "ENCODE"},
    {"key": "synapse",       "re": r"\bsyn\d{7,9}\b",                     "repo": "Synapse"},
    # Papers cite either the scoreset URN (…-a-1) or just the experiment URN
    # (…-a); the trailing scoreset index is optional here and normalised later.
    {"key": "mavedb_urn",    "re": r"\burn:mavedb:\d{8}-[a-z0-9]+(?:-\d+)?\b", "repo": "MaveDB"},
    {"key": "figshare_doi",  "re": r"\b10\.6084/m9\.figshare\.\d+(?:\.v\d+)?\b", "repo": "figshare"},
    {"key": "zenodo_doi",    "re": r"\b10\.5281/zenodo\.\d+\b",           "repo": "Zenodo"},
    {"key": "zenodo_record", "re": r"zenodo\.org/records?/(\d+)",         "repo": "Zenodo"},
    {"key": "arrayexpress",  "re": r"\bE-(?:MTAB|GEOD|PROT)-\d+\b",       "repo": "ArrayExpress/BioStudies"},
    {"key": "pride",         "re": r"\bPXD\d{6,9}\b",                     "repo": "PRIDE"},
    {"key": "dbgap",         "re": r"\bphs\d{6}(?:\.v\d+\.p\d+)?\b",      "repo": "dbGaP", "controlled": True},
    {"key": "ega",           "re": r"\bEGA[SD]\d{11}\b",                  "repo": "EGA", "controlled": True},
    {"key": "scp",           "re": r"\bSCP\d{3,6}\b",                     "repo": "Single Cell Portal"},
    {"key": "cellxgene",     "re": r"cellxgene\.cziscience\.com/collections/([0-9a-f\-]{36})", "repo": "CELLxGENE"},
    {"key": "github_repo",   "re": r"github\.com/([\w.\-]+/[\w.\-]+)",    "repo": "GitHub"},
]

# Assay-family detection. ``spec`` is a specificity multiplier: naming a
# particular protocol ("SHARE-seq", "cell hashing") is far stronger evidence
# for a route than a generic modality ("single-cell RNA"), which nearly every
# paper in this corpus mentions somewhere. Without this, generic terms drown
# out the specific one that actually determines the analysis chain.
ASSAY_TERMS: List[Dict[str, Any]] = [
    {"assay": "lentiMPRA",          "spec": 2.5, "terms": ["lentimpra", "lentiviral mpra"]},
    {"assay": "MPRA",               "spec": 2.0, "terms": ["mpra", "massively parallel reporter assay"]},
    {"assay": "STARR-seq",          "spec": 2.5, "terms": ["starr-seq", "starrseq"]},
    {"assay": "SGE",                "spec": 2.5, "terms": ["saturation genome editing", "sge"]},
    {"assay": "VAMP-seq",           "spec": 2.5, "terms": ["vamp-seq", "vampseq"]},
    {"assay": "DMS",                "spec": 2.0, "terms": ["deep mutational scan", "multiplexed assay of variant effect",
                                                "mave", "variant effect map"]},
    # Martyn 2025 calls the same assay "Variant-EFFECTS" / "flow-sorting
    # experiments with CRISPR targeting screens" and never says "Flow-FISH".
    {"assay": "Flow-FISH",          "spec": 2.5, "terms": ["flow-fish", "flowfish",
                                                "flow-sorting", "variant-effects"]},
    {"assay": "CRISPRi screen",     "spec": 1.5, "terms": ["crispri screen", "crispri", "crispr interference"]},
    {"assay": "CRISPR screen",      "spec": 1.2, "terms": ["crispr screen", "crispr knockout screen", "crisprn"]},
    {"assay": "Perturb-seq",        "spec": 2.5, "terms": ["perturb-seq", "perturbseq", "crop-seq", "crispr droplet"]},
    {"assay": "10x Multiome",       "spec": 2.0, "terms": ["10x multiome", "multiome", "joint rna and atac",
                                                "single-nucleus multiome"]},
    {"assay": "SHARE-seq",          "spec": 2.5, "terms": ["share-seq", "shareseq"]},
    {"assay": "SPLiT-seq",          "spec": 2.5, "terms": ["split-seq", "splitseq", "split-pool"]},
    # "demultiplex" is deliberately absent: every multiplexed scRNA paper uses
    # the word, so it fired on atlases that merely mention demultiplexing.
    {"assay": "cell hashing",       "spec": 2.5, "terms": ["cell hashing", "multi-seq", "multiseq",
                                                "hashtag oligo", "hashsolo", "demuxlet"]},
    {"assay": "scRNA-seq",          "spec": 1.0, "terms": ["single-cell rna", "scrna-seq", "snrna-seq",
                                                "single-nucleus rna", "cell atlas"]},
    {"assay": "scATAC-seq",         "spec": 1.0, "terms": ["scatac", "snatac", "single-cell atac"]},
    {"assay": "ChIP-seq",           "spec": 1.5, "terms": ["chip-seq", "chipseq", "chip-atlas"]},
    {"assay": "enhancer-gene",      "spec": 1.5, "terms": ["enhancer-gene", "enhancer to gene", "e2g",
                                                "activity-by-contact", "abc model", "peak-to-gene",
                                                "peak2gene"]},
    {"assay": "GWAS",               "spec": 1.0, "terms": ["genome-wide association", "gwas", "magma",
                                                "polygenic", "fine-mapping"]},
    {"assay": "eQTL",               "spec": 1.0, "terms": ["eqtl", "sqtl", "caqtl"]},
    {"assay": "proteomics",         "spec": 1.5, "terms": ["mass spectrometry", "proteomic", "interactome"]},
]
_ASSAY_SPEC = {a["assay"]: a.get("spec", 1.0) for a in ASSAY_TERMS}

# Not every accession is equally diagnostic. Almost every paper in this corpus
# deposits *something* in GEO, so a GSE barely narrows the route; a MaveDB URN
# or a Synapse id essentially determines it.
ACCESSION_WEIGHT = {
    "mavedb_urn": 8.0, "synapse": 6.0, "figshare_doi": 6.0, "zenodo_doi": 5.0,
    "zenodo_record": 5.0, "igvf_dataset": 5.0, "igvf_file": 5.0,
    "encode_series": 5.0, "cellxgene": 5.0, "scp": 4.0, "arrayexpress": 3.0,
    "pride": 4.0, "geo_series": 1.5, "geo_sample": 1.0, "bioproject": 1.0,
    "sra_project": 1.0, "dbgap": 1.0, "ega": 1.0, "github_repo": 0.5,
}

# Repositories that essentially every assay deposits into. Their presence says
# "this paper released sequencing data" — it does NOT say which assay produced
# it, so they must never decide *between* analysis routes. They are still
# reported as evidence (and still score as a *primary* match, which is what the
# geo_retrieval / synapse / figshare fallback routes are for); they just carry
# no weight when they appear as a route's `support_accessions`.
#
# Without this, a paper's GEO + Synapse ids were credited only to whichever
# routes happened to list those repos (`mpra` was the sole assay route with
# any `support_accessions`), which sent an open-access CRISPRi enhancer screen
# — CRISPRi mentioned 201x, MPRA 6x — to the MPRA chain.
GENERIC_DEPOSIT_ACCESSIONS = {
    "geo_series", "geo_sample", "sra_project", "bioproject",
    "dbgap", "ega", "arrayexpress", "synapse", "github_repo",
}

# Section titles come through in two flavours: prose ("Data availability",
# "Data and code availability") and JATS sec-type/notes-type slugs
# ("data-availability", "code_availability"), so the separator must allow
# whitespace, hyphens, and underscores.
_SEP = r"[\s\-_]+"
_DATA_SECTION_RE = re.compile(
    r"(data|code|software|materials?)" + _SEP
    + r"(and" + _SEP + r"(code|data|software)" + _SEP + r")?availability"
    r"|availability" + _SEP + r"of" + _SEP + r"(data|code)"
    r"|accession" + _SEP + r"(codes?|numbers?)"
    r"|data" + _SEP + r"access(ibility)?\b",
    re.I,
)

# Deterministic numeric-claim extraction. Only quantities with an explicit,
# unambiguous unit are captured — anything vaguer is left to the LLM draft
# pass, which is always marked unconfirmed anyway.
_COUNT_UNITS = (
    "variants|SNVs|SNV|single.nucleotide variants|nucleotide variants|alleles|"
    "missense variants|substitutions|edits|"
    "cells|nuclei|peaks|elements|enhancers|oligos|guides|sgRNAs|"
    "gRNAs|genes|donors|samples|clusters|cell types|celltypes|pairs|"
    "regulatory elements|candidate enhancers|scoresets|datasets|experiments|"
    "perturbations|transcription factors|TFs|loci|SNPs|reads|links|clones"
)
# Journals routinely interpose an italicised gene symbol or a qualifier between
# the count and its unit ("2,268 *VHL* SNVs", "8,000 unique PTEN variants"), so
# allow up to two short intervening tokens.
_COUNT_RE = re.compile(
    r"([\d]{1,3}(?:,\d{3})+|\b\d{2,9}\b)\s+"
    r"(?:(?:unique|distinct|total|possible|individual|[A-Z][A-Za-z0-9\-]{0,11})\s+){0,2}"
    r"(" + _COUNT_UNITS + r")\b"
)
_METRIC_RE = re.compile(
    r"\b(AUPRC|AUROC|AUC|precision|recall|sensitivity|specificity|"
    r"Pearson(?:'s)?\s*r|Spearman(?:'s)?\s*(?:rho|r)|correlation|R2|R\^2|r2)\b"
    r"[^.\n]{0,30}?(\d\.\d{2,4})", re.I)


def fetch_fulltext(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Fetch full text, best source first; degrade to abstract-only if none.

    Two sources, because no single one covers both halves of this corpus:

    * **Europe PMC** ``/{PMCID}/fullTextXML`` — JATS XML for open-access
      articles in PMC. Verified: this endpoint takes the PMCID directly; the
      ``/{source}/{id}/`` form 404s.
    * **bioRxiv / medRxiv** ``/content/<doi>v<n>.full`` — HTML. Europe PMC
      indexes preprints (``PPR…`` ids) but does **not** serve their full text,
      so preprints need the publisher page.
    """
    base = _resolve_endpoint("europepmc", "EUROPEPMC_BASE")
    out: Dict[str, Any] = {"fulltext_source": None, "xml": "", "text": "",
                            "sections": [], "degraded_reason": None}

    pmcids: List[str] = []
    if rec.get("pmcid"):
        pmcids.append(str(rec["pmcid"]))
    if rec.get("doi") and not pmcids:
        for hit in europepmc_query('DOI:"%s"' % rec["doi"], limit=3):
            if hit.get("pmcid"):
                pmcids.append(hit["pmcid"])

    for pmcid in pmcids:
        status, body = _ref.cached_get(
            "europepmc", f"{base}/{pmcid}/fullTextXML", accept_json=False)
        if status == 200 and body and body.lstrip().startswith("<"):
            out["fulltext_source"] = f"europepmc:{pmcid}/fullTextXML"
            out["xml"] = body
            break

    if out["xml"]:
        try:
            out["sections"], out["text"] = _parse_jats(out["xml"])
        except Exception as e:  # malformed XML — keep the raw text
            logging.warning("JATS parse failed (%s); using raw text", e)
            out["text"] = re.sub(r"<[^>]+>", " ", out["xml"])
        return out

    # Preprint fallback.
    doi = (rec.get("doi") or "").lower()
    if doi.startswith("10.1101/"):
        version = rec.get("version") or ""
        versions = [str(version)] if version else ["1"]
        for server in ("biorxiv", "medrxiv"):
            for v in versions:
                url = f"https://www.{server}.org/content/{doi}v{v}.full"
                status, body = _ref.cached_get(
                    "biorxiv", url, accept_json=False,
                    headers={"Accept": "text/html"})
                if status == 200 and body and len(body) > 5000:
                    out["fulltext_source"] = url
                    out["sections"], out["text"] = _parse_html_article(body)
                    return out

    out["degraded_reason"] = (
        "no open-access full text available (Europe PMC serves JATS only for "
        "PMC open-access articles; this paper is closed-access, not yet "
        "indexed, or a preprint whose publisher page could not be read) — "
        "harvest ran on title + abstract only"
    )
    out["text"] = " ".join(
        x for x in [rec.get("title", ""), rec.get("abstract", "")] if x)
    out["sections"] = [{"title": "Abstract", "text": rec.get("abstract", "")}]
    return out


_TAG_STRIP_RE = re.compile(
    r"<(script|style|nav|footer|header)\b[^>]*>.*?</\1>", re.I | re.S)
_HEADING_RE = re.compile(r"<h([1-4])\b[^>]*>(.*?)</h\1>", re.I | re.S)


def _parse_html_article(html: str) -> Tuple[List[Dict[str, str]], str]:
    """Split a preprint HTML page into (sections, flat text) on its headings."""
    html = _TAG_STRIP_RE.sub(" ", html)

    def _plain(chunk: str) -> str:
        chunk = re.sub(r"<[^>]+>", " ", chunk)
        chunk = (chunk.replace("&nbsp;", " ").replace("&amp;", "&")
                       .replace("&lt;", "<").replace("&gt;", ">")
                       .replace("&#x2019;", "'").replace("&quot;", '"'))
        return " ".join(chunk.split())

    marks = [(m.start(), m.end(), _plain(m.group(2)))
             for m in _HEADING_RE.finditer(html)]
    sections: List[Dict[str, str]] = []
    if not marks:
        return [{"title": "(full text)", "text": _plain(html)}], _plain(html)
    for i, (_, end, title) in enumerate(marks):
        stop = marks[i + 1][0] if i + 1 < len(marks) else len(html)
        text = _plain(html[end:stop])
        if text:
            sections.append({"title": title or "(untitled section)", "text": text})
    return sections, " ".join(s["text"] for s in sections)


def _parse_jats(xml: str) -> Tuple[List[Dict[str, str]], str]:
    """Split JATS full-text XML into (sections, flat text)."""
    root = ET.fromstring(xml)

    def _flat(el) -> str:
        return " ".join("".join(el.itertext()).split())

    sections: List[Dict[str, str]] = []
    for tag, label in (("abstract", "Abstract"),):
        for el in root.iter(tag):
            txt = _flat(el)
            if txt:
                sections.append({"title": label, "text": txt})

    body = root.find(".//body")
    if body is not None:
        for sec in body.iter("sec"):
            t_el = sec.find("title")
            title = _flat(t_el) if t_el is not None else ""
            text = _flat(sec)
            if text:
                sections.append({"title": title or "(untitled section)", "text": text})

    # Back matter carries most journals' Data Availability statements.
    for bm in root.iter("back"):
        for sec in bm.iter("sec"):
            t_el = sec.find("title")
            title = _flat(t_el) if t_el is not None else ""
            text = _flat(sec)
            if text:
                sections.append({"title": title or "(back matter)", "text": text})
        for notes in bm.iter("notes"):
            text = _flat(notes)
            if text:
                sections.append({"title": notes.get("notes-type") or "notes",
                                  "text": text})

    flat = " ".join(s["text"] for s in sections)
    return sections, flat


def extract_accessions(text: str, sections: List[Dict[str, str]]) -> Dict[str, List[Dict[str, Any]]]:
    """Deterministic accession extraction, with the sentence each came from."""
    found: Dict[str, List[Dict[str, Any]]] = {}
    for spec in ACCESSION_PATTERNS:
        pat = re.compile(spec["re"], re.I if spec["key"] != "mavedb_urn" else 0)
        hits: Dict[str, Dict[str, Any]] = {}
        for m in pat.finditer(text):
            val = m.group(1) if m.groups() else m.group(0)
            val = val.strip()
            norm = val.upper() if spec["key"] not in ("synapse", "github_repo",
                                                       "cellxgene", "mavedb_urn",
                                                       "figshare_doi", "zenodo_doi") else val
            if norm in hits:
                hits[norm]["count"] += 1
                continue
            start = max(0, m.start() - 200)
            end = min(len(text), m.end() + 200)
            hits[norm] = {
                "value": norm,
                "repo": spec["repo"],
                "controlled": bool(spec.get("controlled")),
                "count": 1,
                "context": text[start:end].strip(),
                "in_data_availability": False,
            }
        if hits:
            found[spec["key"]] = list(hits.values())

    # Flag accessions that appear inside a Data/Code Availability section —
    # those are the primary deposits, not incidental citations.
    da_text = " ".join(s["text"] for s in sections
                        if _DATA_SECTION_RE.search(s.get("title", "")))
    if da_text:
        for entries in found.values():
            for e in entries:
                if e["value"] in da_text or e["value"].lower() in da_text.lower():
                    e["in_data_availability"] = True
    return found


def extract_assays(text: str) -> List[Dict[str, Any]]:
    low = text.lower()
    out = []
    for spec in ASSAY_TERMS:
        n = sum(low.count(t) for t in spec["terms"])
        if n:
            out.append({"assay": spec["assay"], "mentions": n})
    out.sort(key=lambda d: d["mentions"], reverse=True)
    return out


def extract_data_availability(sections: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [
        {"title": s.get("title", ""), "text": s["text"][:4000]}
        for s in sections if _DATA_SECTION_RE.search(s.get("title", ""))
    ]


def _sentences(text: str) -> List[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text) if s.strip()]


def extract_numeric_claims(sections: List[Dict[str, str]],
                            max_claims: int = 40) -> List[Dict[str, Any]]:
    """Deterministic candidate ground-truth claims, with their source sentence.

    These are *candidates*: the regex knows a number and a unit stood next to
    each other, not that the number is the paper's headline result. Everything
    produced here lands in ``expected.json`` as ``confirmed: false``.
    """
    priority = {"abstract": 0, "result": 1, "discussion": 3, "method": 2}

    def _rank(title: str) -> int:
        t = (title or "").lower()
        for k, v in priority.items():
            if k in t:
                return v
        return 4

    claims: List[Dict[str, Any]] = []
    seen = set()
    for sec in sorted(sections, key=lambda s: _rank(s.get("title", ""))):
        for sent in _sentences(sec["text"]):
            if len(sent) > 600:
                continue
            for m in _COUNT_RE.finditer(sent):
                raw, unit = m.group(1), m.group(2)
                try:
                    val = int(raw.replace(",", ""))
                except ValueError:
                    continue
                if val < 2:
                    continue
                key = (val, unit.lower())
                if key in seen:
                    continue
                seen.add(key)
                claims.append({
                    "kind": "count",
                    "value": val,
                    "unit": unit,
                    "section": sec.get("title", ""),
                    "quote": sent[:400],
                })
            for m in _METRIC_RE.finditer(sent):
                metric, raw = m.group(1), m.group(2)
                try:
                    val = float(raw)
                except ValueError:
                    continue
                key = (val, metric.lower())
                if key in seen:
                    continue
                seen.add(key)
                claims.append({
                    "kind": "metric",
                    "value": val,
                    "unit": metric,
                    "section": sec.get("title", ""),
                    "quote": sent[:400],
                })
            if len(claims) >= max_claims:
                return claims
    return claims


# ----------------------------------------------------------------------------
# 5. LLM draft pass (optional, always marked unconfirmed)
# ----------------------------------------------------------------------------

_LLM_PROMPT = """You are helping build a REPRODUCIBILITY BENCHMARK for a paper.

Below are the paper's abstract, its Data Availability statement, the assay
families detected in its text, and numeric quantities a regex found next to a
unit. Your job is to pick the small set of numbers that are the paper's own
HEADLINE, CHECKABLE RESULTS — the kind another group would try to reproduce.

Rules:
- Only propose a claim if the supplied quote actually states it. Never infer.
- Prefer counts of the primary analysed entity (variants scored, cells in the
  final atlas, elements tested, links called) and headline performance metrics.
- Skip sample-prep trivia, reagent volumes, read counts, and citation numbers.
- At most 6 claims.

Return ONLY a JSON array, no prose, each element:
  {"description": "...", "value": <number>, "unit": "...",
   "tolerance_frac": <0.0-0.5>, "quote": "<verbatim sentence from the input>"}

ABSTRACT:
%(abstract)s

DATA AVAILABILITY:
%(availability)s

DETECTED ASSAYS: %(assays)s

REGEX CANDIDATES:
%(candidates)s
"""


def llm_draft_claims(harvest: Dict[str, Any], *, backend: Optional[str] = None,
                      model: Optional[str] = None) -> Dict[str, Any]:
    """Ask the configured LLM to rank which numeric claims are headline results.

    Failure here is never fatal — the deterministic claims stand on their own.
    """
    result: Dict[str, Any] = {"enabled": True, "claims": [], "error": None,
                               "backend": backend, "model": model}
    try:
        try:
            from igvfagent import _llm  # type: ignore
        except Exception:
            import _llm  # type: ignore
    except Exception as e:
        result["enabled"] = False
        result["error"] = f"LLM layer unavailable: {e}"
        return result

    cands = harvest.get("numeric_claims") or []
    if not cands:
        result["error"] = "no regex candidates to rank"
        return result

    prompt = _LLM_PROMPT % {
        "abstract": (harvest.get("abstract") or "")[:4000],
        "availability": " ".join(
            d["text"] for d in harvest.get("data_availability") or [])[:2500]
            or "(none found)",
        "assays": ", ".join(a["assay"] for a in harvest.get("assays") or []) or "(none)",
        "candidates": json.dumps(
            [{k: c[k] for k in ("kind", "value", "unit", "section", "quote")}
             for c in cands[:30]], indent=1)[:8000],
    }
    try:
        msg = _llm.chat(
            [{"role": "user", "content": prompt}],
            backend=backend, model=model, tools=[],
            max_tokens=2048, temperature=0.0,
        )
    except Exception as e:
        result["enabled"] = False
        result["error"] = f"LLM call failed: {e}"
        return result

    text = (msg.content or "").strip()
    m = re.search(r"\[[\s\S]*\]", text)
    if not m:
        result["error"] = "LLM response contained no JSON array"
        return result
    try:
        parsed = json.loads(m.group(0))
    except Exception as e:
        result["error"] = f"LLM JSON parse failed: {e}"
        return result

    quotes = {c["quote"] for c in cands}
    for item in parsed if isinstance(parsed, list) else []:
        if not isinstance(item, dict) or "value" not in item:
            continue
        q = str(item.get("quote", ""))
        # Guard against fabricated quotes: keep the flag, do not drop silently.
        grounded = any(q[:80] and q[:80] in full for full in quotes)
        result["claims"].append({
            "description": str(item.get("description", ""))[:200],
            "value": item.get("value"),
            "unit": str(item.get("unit", ""))[:40],
            "tolerance_frac": item.get("tolerance_frac", 0.1),
            "quote": q[:400],
            "quote_grounded_in_source": grounded,
        })
    return result


# ----------------------------------------------------------------------------
# 6. Routing table
# ----------------------------------------------------------------------------
#
# Every route below is derived from a benchmark that already exists in
# `Benchmarks/` — the CLI invocations are the ones those run.sh files actually
# use, not invented syntax. `match` scores a route against the harvest;
# `vars` declares the shell variables run.sh needs and where to fill them from.

IGVF = '"$IGVF"'

ROUTES: List[Dict[str, Any]] = [
    {
        "name": "mavedb_scoreset",
        "title": "MaveDB scoreset → genomic coordinates (variant-effect map)",
        "skill_output_dir": "MaveDB",
        "primary_artefact": "summary.json",
        "reference_benchmark": "buckley2024_vhl / waters2024_bap1 / matreyek2018_pten_vampseq",
        "match": {"accessions": ["mavedb_urn"],
                   "assays": ["SGE", "VAMP-seq", "DMS"]},
        "vars": [
            {"name": "URN", "from": "accessions.mavedb_urn",
             "prompt": "MaveDB scoreset URN (find it at mavedb.org)"},
            {"name": "GENE", "from": "genes",
             "prompt": "HGNC symbol of the assayed gene"},
        ],
        "steps": [
            f'{IGVF} mavedb map-scoreset --urn "$URN" --gene "$GENE" --label "$LABEL"',
            f'{IGVF} catalog get-entity "$GENE" || true',
            f'{IGVF} catalog find-associations "$GENE" --relationship genetic --limit 10 || true',
        ],
        "checks": [
            {"name": "MaveDB scoreset downloaded + parsed (rows in)",
             "type": "range", "path": "n_rows_in", "min": 1},
            {"name": "Per-variant TSV written + non-empty",
             "type": "artefact", "filename": "*_mapped.tsv"},
            {"name": "VCF artefact written",
             "type": "artefact", "filename": "*_mapped.vcf"},
        ],
        "followups": [
            "calibrate run  — turn the assay scores into ACMG/AMP PS3/BS3 evidence",
        ],
    },
    {
        "name": "mpra",
        "title": "MPRA / lentiMPRA activity + skew",
        "skill_output_dir": "MPRA",
        "primary_artefact": "summary.json",
        "reference_benchmark": "agarwal2025_lentimpra / deng2024_cortex_mpra",
        "match": {"assays": ["lentiMPRA", "MPRA", "STARR-seq"],
                   "support_accessions": ["geo_series", "synapse", "igvf_dataset"]},
        "vars": [
            {"name": "COUNTS", "from": "local_input",
             "default": "Benchmarks/_data/{paper_id}/oligo_counts.tsv",
             "prompt": "per-oligo DNA + RNA count table"},
        ],
        "local_input": {
            "var": "COUNTS",
            "hint": "Per-oligo DNA + RNA counts from the paper's deposit "
                    "(GEO supplementary file, or the Synapse/PsychENCODE folder).",
        },
        "steps": [
            f'{IGVF} mpra portal-manifest --limit 50 --label "$LABEL" || true',
            {"cmd": f'{IGVF} mpra activity --counts "$COUNTS" --label "$LABEL"',
             "needs_local": True},
            {"cmd": f'{IGVF} mpra qc --counts "$COUNTS" --label "$LABEL"',
             "needs_local": True},
            {"cmd": f'{IGVF} mpra volcano --label "$LABEL"', "needs_local": True},
        ],
        "checks": [
            {"name": "MPRA discovery manifest written",
             "type": "artefact", "filename": "*manifest*.csv"},
        ],
        "extra_search_dirs": ["Data/Manifests/MPRA", "Data"],
    },
    {
        "name": "flowfish",
        "title": "CRISPRi Flow-FISH element scoring",
        "skill_output_dir": "FlowFISH",
        "primary_artefact": "summary.json",
        "reference_benchmark": "martyn2025_variant_flowfish",
        "match": {"assays": ["Flow-FISH"]},
        "vars": [],
        "steps": [
            f'{IGVF} flowfish pull-portal --limit 500 --label "$LABEL"',
            'SIM_DIR="Benchmarks/_data/$LABEL/sim"',
            'mkdir -p "$SIM_DIR"',
            f'{IGVF} flowfish simulate --out-dir "$SIM_DIR" '
            f'--n-elements 20 --guides-per-element 5 --knockdown-frac 0.5 '
            f'--cells-per-guide 200 --seed 42',
            f'{IGVF} flowfish estimate-effects --counts "$SIM_DIR/counts.tsv" '
            f'--sortparams "$SIM_DIR/sortparams.tsv" --label "${{LABEL}}_pipeline"',
        ],
        "checks": [
            {"name": "Flow-FISH portal enumeration written",
             "type": "artefact", "filename": "*.tsv"},
        ],
        "notes": "The chain scaffolds with a SIMULATED screen (the smoke-test "
                  "path from martyn2025_variant_flowfish). To actually reproduce "
                  "the paper, replace the simulate step with the paper's own "
                  "guide × sorted-bin counts, then continue with real-space and "
                  "score-elements.",
    },
    {
        "name": "perturb_catalog",
        "title": "Perturbation Catalogue census (CRISPR screen / Perturb-seq)",
        "skill_output_dir": "Perturbation",
        "primary_artefact": "summary.json",
        "reference_benchmark": "weinstock2024_cd4_crispr / joung2025_tf_perturbseq",
        "match": {"assays": ["CRISPR screen", "CRISPRi screen", "Perturb-seq"]},
        "vars": [
            {"name": "MODALITY", "from": "const",
             "default": "crispr-screen",
             "prompt": "perturb-catalog modality (crispr-screen | perturb-seq | mave)"},
        ],
        "steps": [
            f'{IGVF} perturb-catalog summary',
            f'{IGVF} perturb-catalog search-modality --modality "$MODALITY" '
            f'--dataset-limit 100',
        ],
        "checks": [
            {"name": "Perturbation Catalogue summary written",
             "type": "artefact", "filename": "summary.json"},
        ],
        # Neither subcommand accepts --label, so the run dir is not
        # paper-tagged; score against the skill's own default dir name.
        "label_hint": "summary",
    },
    {
        "name": "geo_retrieval",
        "tier": "retrieval",
        "title": "GEO series metadata + supplementary-file inventory",
        "skill_output_dir": "SingleCell",
        "primary_artefact": "summary.json",
        "reference_benchmark": "zheng2024_invivo_perturbseq",
        "match": {"accessions": ["geo_series"]},
        "vars": [
            {"name": "GSE", "from": "accessions.geo_series",
             "prompt": "GEO series accession"},
        ],
        "steps": [f'{IGVF} geo series --gse "$GSE"'],
        "checks": [
            {"name": "GEO series metadata retrieved",
             "type": "artefact", "filename": "*.json"},
        ],
        "notes": "Retrieval only. Pair with an analytical route once the "
                  "supplementary files are in a format IGVFagent reads "
                  "(h5ad / mtx / TSV — not .rds or .qs).",
    },
    {
        "name": "sc_analyze",
        "title": "Single-cell atlas reproduction (QC → HVG → PCA → UMAP → Leiden → markers)",
        "skill_output_dir": "SingleCell",
        "primary_artefact": "summary.json",
        "reference_benchmark": "travaglini2020_lung / wang2025_neocortex_multiome",
        "match": {"assays": ["scRNA-seq"],
                   "accessions": ["cellxgene"],
                   "support_accessions": ["geo_series"]},
        "vars": [
            {"name": "H5AD", "from": "local_input",
             "default": "Benchmarks/_data/{paper_id}/atlas.h5ad",
             "prompt": "the paper's cell-by-gene matrix as .h5ad"},
            {"name": "CELLTYPE_COL", "from": "const", "default": "cell_type",
             "prompt": "obs column holding the authors' own cell-type labels"},
        ],
        "local_input": {
            "var": "H5AD",
            "hint": "CELLxGENE distributes most published atlases as .h5ad "
                    "(cellxgene.cziscience.com → collection → download). For "
                    "GEO deposits, convert the matrix to .h5ad first.",
        },
        "steps": [
            {"cmd": f'{IGVF} sc-analyze pipeline --input "$H5AD" --label "$LABEL" '
                     f'--resolution 1.0 --n-hvg 2000 --n-pcs 50 --skip-tsne '
                     f'--sample-col "$CELLTYPE_COL"',
             "needs_local": True},
        ],
        "checks": [
            {"name": "Single-cell pipeline summary written",
             "type": "artefact", "filename": "summary.json"},
        ],
        "notes": "Concordance against the authors' labels (ARI / AMI / "
                  "homogeneity) is the reproducibility claim — add those as "
                  "range checks once you have run it once.",
    },
    {
        "name": "multiome_peak2gene",
        "title": "10x Multiome peak→gene cis-regulatory linkage",
        "skill_output_dir": "Multiome10x",
        "primary_artefact": "summary.json",
        "reference_benchmark": "trevino2021_cortex_multiome / mitra2024_scarlink",
        "match": {"assays": ["10x Multiome", "scATAC-seq", "enhancer-gene"],
                   "support_accessions": ["igvf_dataset", "geo_series"]},
        "vars": [],
        "steps": [
            f'{IGVF} multiome retrieve --count 5 --label "$LABEL"',
        ],
        "checks": [
            {"name": "Multiome AnalysisSet enumeration written",
             "type": "artefact", "filename": "summary.json"},
        ],
        "notes": "For a full local reproduction add `multiome process-local "
                  "--input <bundle>` then `multiome peak2gene --label $LABEL`.",
    },
    {
        "name": "share_seq",
        "title": "SHARE-seq joint RNA + ATAC QC",
        "skill_output_dir": "SHAREseq",
        "primary_artefact": "summary.json",
        "reference_benchmark": "ma2020_shareseq",
        "match": {"assays": ["SHARE-seq"]},
        "vars": [
            {"name": "H5AD", "from": "local_input",
             "default": "Benchmarks/_data/{paper_id}/shareseq_rna.h5ad",
             "prompt": "SHARE-seq RNA matrix as .h5ad"},
        ],
        "local_input": {"var": "H5AD",
                         "hint": "SHARE-seq RNA DGE from the paper's GEO series."},
        "steps": [
            {"cmd": f'{IGVF} share rna-qc --h5ad "$H5AD" --label "$LABEL"',
             "needs_local": True},
        ],
        "checks": [{"name": "SHARE-seq QC summary written",
                     "type": "artefact", "filename": "summary.json"}],
    },
    {
        "name": "splitseq",
        "title": "SPLiT-seq combinatorial-barcode atlas",
        "skill_output_dir": "SPLiTseq",
        "primary_artefact": "summary.json",
        "reference_benchmark": "rosenberg2018_splitseq",
        "match": {"assays": ["SPLiT-seq"]},
        "vars": [
            {"name": "DGE", "from": "local_input",
             "default": "Benchmarks/_data/{paper_id}/dge.mtx",
             "prompt": "SPLiT-seq DGE matrix"},
        ],
        "local_input": {"var": "DGE",
                         "hint": "SPLiT-seq DGE from the paper's GEO series."},
        "steps": [
            {"cmd": f'{IGVF} splitseq analyze --input "$DGE" --label "$LABEL"',
             "needs_local": True},
        ],
        "checks": [{"name": "SPLiT-seq analysis summary written",
                     "type": "artefact", "filename": "summary.json"}],
    },
    {
        "name": "multiseq",
        "title": "Cell hashing / MULTI-seq demultiplexing",
        "skill_output_dir": "MultiSeq",
        "primary_artefact": "summary.json",
        "reference_benchmark": "demultiplex2_stoeckius",
        "match": {"assays": ["cell hashing"]},
        "vars": [
            {"name": "TAGS", "from": "local_input",
             "default": "Benchmarks/_data/{paper_id}/tags.csv",
             "prompt": "cell × barcode tag count matrix (CSV)"},
        ],
        "local_input": {"var": "TAGS",
                         "hint": "HTO / MULTI-seq tag count matrix from the paper's deposit."},
        "steps": [
            {"cmd": f'{IGVF} multiseq demultiplex --input "$TAGS" --label "$LABEL"',
             "needs_local": True},
        ],
        "checks": [{"name": "Demultiplexing summary written",
                     "type": "artefact", "filename": "summary.json"}],
    },
    {
        "name": "chipatlas",
        "title": "ChIP-Atlas TF / regulatory-element census",
        "skill_output_dir": "ChIPAtlas",
        "primary_artefact": "summary.json",
        "reference_benchmark": "zou2024_chipatlas_gata1",
        "match": {"assays": ["ChIP-seq"]},
        "vars": [
            {"name": "ANTIGEN", "from": "genes", "prompt": "TF / antigen symbol"},
            {"name": "GENOME", "from": "const", "default": "hg38",
             "prompt": "genome assembly"},
            {"name": "AG_CLASS", "from": "const", "default": "TFs and others",
             "prompt": "antigen class (Histone | TFs and others | ATAC-Seq | DNase-seq)"},
            {"name": "CELL_CLASS", "from": "const", "default": "All cell types",
             "prompt": "ChIP-Atlas cell class"},
        ],
        "steps": [
            f'{IGVF} chipatlas list-antigens --genome "$GENOME" '
            f'--ag-class "$AG_CLASS" --cell-class "$CELL_CLASS" --limit 40',
            f'{IGVF} chipatlas search --query "$ANTIGEN" --genome "$GENOME" --limit 10',
            f'{IGVF} chipatlas assemble-bed --genome "$GENOME" '
            f'--ag-class "$AG_CLASS" --antigen "$ANTIGEN" '
            f'--cell-class "$CELL_CLASS" --qval 05',
        ],
        "checks": [{"name": "assemble-bed response written",
                     "type": "artefact", "filename": "assemble.json"}],
        # No chipatlas subcommand accepts --label (see zou2024_chipatlas_gata1).
        "label_hint": "assemble",
    },
    {
        "name": "encode_fce",
        "title": "ENCODE functional-characterization experiment enumeration",
        "skill_output_dir": "ENCODE",
        "primary_artefact": "summary.json",
        "reference_benchmark": "yao2024_encode4_crispri",
        "match": {"accessions": ["encode_series"],
                   "assays": ["CRISPRi screen", "CRISPR screen"]},
        "vars": [
            {"name": "ASSAY", "from": "const", "default": "CRISPR screen",
             "prompt": "ENCODE assay title"},
        ],
        "steps": [
            f'{IGVF} encode retrieve --assay "$ASSAY" --label "$LABEL"',
        ],
        "checks": [{"name": "ENCODE enumeration written",
                     "type": "artefact", "filename": "summary.json"}],
    },
    {
        "name": "synapse",
        "tier": "retrieval",
        "title": "Synapse / PsychENCODE deposit walk",
        "skill_output_dir": "Synapse",
        "primary_artefact": "summary.json",
        "reference_benchmark": "deng2024_cortex_mpra",
        "match": {"accessions": ["synapse"]},
        "vars": [
            {"name": "SYN", "from": "accessions.synapse",
             "prompt": "Synapse entity id"},
        ],
        "steps": [
            f'{IGVF} synapse walk --syn "$SYN" --max-depth 3 --label "$LABEL"',
        ],
        "checks": [{"name": "Synapse walk written",
                     "type": "artefact", "filename": "summary.json"}],
        "notes": "Public folders read anonymously; controlled cohorts "
                  "(PsychENCODE, AMP-AD, AMP-PD) need SYNAPSE_AUTH_TOKEN and "
                  "an accepted data-use agreement.",
    },
    {
        "name": "figshare",
        "tier": "retrieval",
        "title": "figshare / Zenodo deposit retrieval",
        "skill_output_dir": "Figshare",
        "primary_artefact": "summary.json",
        "reference_benchmark": "liu2025_open4gene",
        "match": {"accessions": ["figshare_doi", "zenodo_doi", "zenodo_record"]},
        "vars": [
            {"name": "ARTICLE", "from": "accessions.figshare_doi",
             "prompt": "figshare article id or DOI"},
        ],
        "steps": [
            f'{IGVF} figshare article --id "$ARTICLE" --label "$LABEL"',
            f'{IGVF} figshare files --id "$ARTICLE" || true',
        ],
        "checks": [{"name": "figshare article metadata written",
                     "type": "artefact", "filename": "*.json"}],
    },
    {
        "name": "portal_discovery",
        "tier": "retrieval",
        "title": "IGVF Portal faceted discovery (fallback route)",
        "skill_output_dir": "Portal",
        "primary_artefact": "summary.json",
        "reference_benchmark": "mitra2024_scarlink",
        "match": {"accessions": ["igvf_dataset", "igvf_file"]},
        "vars": [
            {"name": "ACCESSION", "from": "accessions.igvf_dataset",
             "prompt": "IGVF Portal accession"},
        ],
        # `portal get` takes the id positionally and accepts no --label.
        "steps": [
            f'{IGVF} portal get "$ACCESSION"',
        ],
        "checks": [{"name": "Portal record retrieved",
                     "type": "artefact", "filename": "*.json"}],
        "label_hint": "get",
    },
]


def route(harvest: Dict[str, Any], *, top: int = 3) -> Dict[str, Any]:
    """Rank routes against a harvest. Accession evidence outweighs prose."""
    acc = harvest.get("accessions") or {}
    assays = {a["assay"]: a["mentions"] for a in (harvest.get("assays") or [])}
    total_assay_mentions = float(sum(assays.values())) or 1.0
    ranked: List[Dict[str, Any]] = []

    for spec in ROUTES:
        score = 0.0
        why: List[str] = []

        def _acc_score(key: str, *, support: bool = False) -> float:
            entries = acc.get(key) or []
            if not entries:
                return 0.0
            in_da = any(e.get("in_data_availability") for e in entries)
            generic = support and key in GENERIC_DEPOSIT_ACCESSIONS
            why.append(f"{key}={entries[0]['value']}"
                        + (" (Data Availability)" if in_da else " (in text)")
                        + (" [generic deposit, unscored]" if generic else ""))
            # A repository every assay uses cannot discriminate between assay
            # routes — report it, but do not let it rank them.
            if generic:
                return 0.0
            # Being named in a Data/Code Availability section means it is the
            # paper's own deposit, not a citation of someone else's.
            return ACCESSION_WEIGHT.get(key, 3.0) * (2.0 if in_da else 1.0)

        # Assay evidence first — it is what actually determines the chain.
        assay_score = 0.0
        for a in spec["match"].get("assays", []) or []:
            if a in assays:
                n = assays[a]
                # Three signals, in decreasing order of trust:
                #   spec       how diagnostic the term itself is
                #   log(n)     how much the paper talks about it — log-scaled so
                #              a heavily-used assay separates from an incidental
                #              mention without a single term running away
                #   share      that assay's fraction of ALL assay mentions, which
                #              is what actually distinguishes "this paper IS a
                #              CRISPRi screen" from "it cites one"
                assay_score += (3.0 * _ASSAY_SPEC.get(a, 1.0)
                                 + 1.2 * math.log10(1.0 + n)
                                 + 2.0 * (n / total_assay_mentions))
                why.append(f"assay:{a}×{n}")
        score += assay_score

        # Primary accessions are diagnostic on their own (a MaveDB URN means
        # a MaveDB scoreset, full stop).
        for key in spec["match"].get("accessions", []) or []:
            score += _acc_score(key)

        # Supporting accessions only count once an assay has already matched.
        # Otherwise a paper that merely deposits in GEO or Synapse would be
        # routed to whichever assay route happens to list those repositories —
        # which sent a lung cell atlas to the MPRA chain.
        if assay_score > 0:
            for key in spec["match"].get("support_accessions", []) or []:
                score += _acc_score(key, support=True)

        if score <= 0:
            continue
        ranked.append({
            "route": spec["name"],
            "title": spec["title"],
            "score": round(score, 2),
            "evidence": why,
            "tier": spec.get("tier", "analysis"),
            "assay_matched": assay_score > 0,
            "skill_output_dir": spec["skill_output_dir"],
            "reference_benchmark": spec.get("reference_benchmark"),
            "requires_local_input": bool(spec.get("local_input")),
            "notes": spec.get("notes"),
        })

    # Analysis routes that matched an assay outrank pure-retrieval routes
    # regardless of score. Retrieval (Synapse / figshare / GEO / Portal) is
    # what you fall back to when the assay is unknown — a paper depositing in
    # Synapse should still be routed to its analysis chain when we can tell
    # what the assay was. Within each tier, score decides.
    ranked.sort(key=lambda d: (
        0 if (d["tier"] == "analysis" and d["assay_matched"]) else 1,
        -d["score"],
    ))
    controlled = sorted({
        e["repo"] for entries in acc.values() for e in entries if e.get("controlled")
    })
    return {
        "plans": ranked[:top],
        "all_scored": ranked,
        "selected": ranked[0]["route"] if ranked else None,
        "controlled_access_repos": controlled,
        "unroutable_reason": None if ranked else (
            "No accession pattern and no known assay family matched. Either the "
            "full text was unavailable (check harvest.degraded_reason) or this "
            "paper's data type is outside IGVFagent's covered assay families."
        ),
    }


def get_route(name: str) -> Dict[str, Any]:
    for spec in ROUTES:
        if spec["name"] == name:
            return spec
    raise SystemExit(f"unknown route: {name!r}. See `igvfagent bench list-routes`.")


# ----------------------------------------------------------------------------
# 7. Scaffolding
# ----------------------------------------------------------------------------

def _resolve_vars(spec: Dict[str, Any], harvest: Dict[str, Any],
                   paper_id: str) -> List[Dict[str, Any]]:
    """Fill each route variable from the harvest; flag the ones that need a human."""
    acc = harvest.get("accessions") or {}
    out = []
    for v in spec.get("vars") or []:
        name, src = v["name"], v.get("from")
        value, origin = None, "unresolved"
        if src == "const":
            value, origin = v.get("default"), "route default"
        elif src == "local_input":
            value = (v.get("default") or "").format(paper_id=paper_id)
            origin = "local path (must be downloaded)"
        elif src == "genes":
            genes = harvest.get("genes") or []
            if genes:
                value, origin = genes[0], "gene symbol extracted from title/abstract"
        elif src and src.startswith("accessions."):
            key = src.split(".", 1)[1]
            entries = acc.get(key) or []
            entries = sorted(entries,
                              key=lambda e: (e.get("in_data_availability"), e["count"]),
                              reverse=True)
            if entries:
                value = entries[0]["value"]
                origin = ("Data Availability statement"
                           if entries[0].get("in_data_availability") else "paper text")
                # A MaveDB *experiment* URN (urn:mavedb:00000675-a) is not a
                # scoreset; map-scoreset needs the scoreset URN. Papers cite
                # the experiment form often enough to be worth handling, but
                # "-1" is an assumption, so say so rather than hide it.
                if key == "mavedb_urn" and re.search(r"-[a-z0-9]+$", value) \
                        and not re.search(r"-\d+$", value):
                    origin += (f" — paper cited the experiment URN {value}; "
                                f"assuming scoreset -1 (VERIFY on mavedb.org)")
                    value = value + "-1"
        out.append({"name": name, "value": value, "origin": origin,
                     "prompt": v.get("prompt", ""),
                     "needs_human": value is None})
    return out


def _render_run_sh(paper_id: str, spec: Dict[str, Any],
                    variables: List[Dict[str, Any]],
                    resolution: Dict[str, Any]) -> str:
    paper = resolution.get("paper") or {}
    lines = [
        "#!/usr/bin/env bash",
        f"# {paper_id} — {spec['title']}",
        f"# Paper: {paper.get('title', '?')}",
        f"#        {paper.get('journal', '?')} {paper.get('year', '?')} "
        f"· doi:{paper.get('doi', '?')}",
        "#",
        "# GENERATED by `igvfagent bench scaffold`. Review before trusting the",
        "# numbers: variables marked TODO_VERIFY below were not resolvable from",
        "# the paper's text and need a human to fill in.",
        "set -euo pipefail",
        'HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'ROOT="$(cd "$HERE/../.." && pwd)"',
        'cd "$ROOT"',
        f'LABEL="{paper_id}"',
        'IGVF=".venv/bin/igvfagent"',
        "",
    ]

    for v in variables:
        env_name = f"{paper_id.split('_')[0].upper()}_{v['name']}"
        if v["needs_human"]:
            lines.append(f"# TODO_VERIFY: {v['prompt']}")
            lines.append(f'{v["name"]}="${{{env_name}:-}}"')
            lines.append(f'if [ -z "${v["name"]}" ]; then')
            lines.append(f'    echo "[{paper_id}] {v["name"]} is unset — {v["prompt"]}."')
            lines.append(f'    echo "  Set it with: export {env_name}=<value>"')
            lines.append("    exit 77")
            lines.append("fi")
        else:
            lines.append(f"# {v['name']}: from {v['origin']}")
            lines.append(f'{v["name"]}="${{{env_name}:-{v["value"]}}}"')
        lines.append("")

    local = spec.get("local_input")
    has_local = bool(local)
    if has_local:
        var = local["var"]
        lines += [
            f'# Local-input guard — the suite convention is exit 77 for "skipped,',
            f'# missing local input" so run_all.sh reports it honestly.',
            f'if [ ! -f "${var}" ]; then',
            f'    echo "[{paper_id}] Local input not found: ${var}"',
            f'    echo "  {local["hint"]}"',
            f'    echo "  Online discovery steps still ran; analytical steps skipped."',
            "    exit 77",
            "fi",
            "",
        ]

    for step in spec["steps"]:
        cmd = step["cmd"] if isinstance(step, dict) else step
        lines.append(cmd)
        lines.append("")

    lines += [
        'echo ""',
        f'echo "== {paper_id} benchmark complete =="',
        f'echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark {paper_id}"',
        "",
    ]
    return "\n".join(lines)


def _render_expected(paper_id: str, spec: Dict[str, Any],
                      resolution: Dict[str, Any], harvest: Dict[str, Any],
                      variables: List[Dict[str, Any]]) -> Dict[str, Any]:
    paper = resolution.get("paper") or {}
    checks: List[Dict[str, Any]] = []

    # Structural checks the route guarantees when it runs. These are the only
    # checks written as confirmed — they assert "the chain produced its
    # artefacts", never "the paper's number was reproduced".
    for c in spec.get("checks", []):
        chk = dict(c)
        chk["confirmed"] = True
        chk["provenance"] = {"kind": "structural",
                              "note": "asserted by the route, not by the paper"}
        checks.append(chk)

    # Everything derived from the paper's prose is unconfirmed by construction.
    for claim in (harvest.get("numeric_claims") or [])[:12]:
        val = claim["value"]
        tol = 0.05 if claim["kind"] == "count" else 0.15
        checks.append({
            "name": f"[UNCONFIRMED] paper claims {val} {claim['unit']}",
            "type": "range",
            "path": "TODO_SET_JSON_PATH",
            "min": val * (1 - tol),
            "max": val * (1 + tol),
            "expected": val,
            "confirmed": False,
            "provenance": {
                "kind": "regex_from_fulltext",
                "section": claim.get("section", ""),
                "quote": claim.get("quote", ""),
            },
        })
    for claim in (harvest.get("llm_claims") or {}).get("claims", [])[:6]:
        val = claim.get("value")
        if not isinstance(val, (int, float)):
            continue
        tol = float(claim.get("tolerance_frac") or 0.1)
        checks.append({
            "name": f"[UNCONFIRMED] {claim.get('description') or 'LLM-proposed claim'}",
            "type": "range",
            "path": "TODO_SET_JSON_PATH",
            "min": val * (1 - tol),
            "max": val * (1 + tol),
            "expected": val,
            "confirmed": False,
            "provenance": {
                "kind": "llm_draft",
                "backend": (harvest.get("llm_claims") or {}).get("backend"),
                "model": (harvest.get("llm_claims") or {}).get("model"),
                "quote": claim.get("quote", ""),
                "quote_grounded_in_source": claim.get("quote_grounded_in_source"),
            },
        })

    data_source: Dict[str, Any] = {"type": spec["name"]}
    for v in variables:
        if v["value"] is not None:
            data_source[v["name"].lower()] = v["value"]

    out: Dict[str, Any] = {
        "paper": {
            "title": paper.get("title"),
            "authors": paper.get("authors"),
            "journal": paper.get("journal"),
            "year": paper.get("year"),
            "doi": paper.get("doi"),
            "pmid": paper.get("pmid") or None,
            "pmcid": paper.get("pmcid") or None,
        },
        "generated_by": {
            "tool": "igvfagent bench scaffold",
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
            "route": spec["name"],
            "resolver_confidence": resolution.get("confidence"),
            "fulltext_source": harvest.get("fulltext_source"),
            "review_required": True,
        },
        "data_source": data_source,
        "skill_output_dir": spec["skill_output_dir"],
        # concordance.py locates the run via Docs/<skill>/2*_*<label>*/. Most
        # skills accept --label so the paper-id lands in the dir name; a few
        # (chipatlas, perturb-catalog, portal get) do not, and for those the
        # route declares the default dir token the skill actually writes.
        "label": spec.get("label_hint") or paper_id,
        "primary_artefact": spec.get("primary_artefact", "summary.json"),
        "checks": checks,
    }
    if spec.get("label_hint"):
        out["label_note"] = (
            f"This route's CLI steps do not accept --label, so the run "
            f"directory is not paper-tagged; scoring matches the skill's "
            f"default '{spec['label_hint']}' dir and will pick up the most "
            f"recent run of that skill, whichever paper produced it."
        )
    if spec.get("extra_search_dirs"):
        out["extra_search_dirs"] = spec["extra_search_dirs"]
    return out


def _render_readme(paper_id: str, spec: Dict[str, Any],
                    resolution: Dict[str, Any], harvest: Dict[str, Any],
                    variables: List[Dict[str, Any]],
                    routing: Dict[str, Any]) -> str:
    paper = resolution.get("paper") or {}
    unresolved = [v for v in variables if v["needs_human"]]
    L = [
        f"# {paper_id}",
        "",
        "> **Draft — generated by `igvfagent bench scaffold`.** Structural checks "
        "are live; every check tagged `[UNCONFIRMED]` in `expected.json` came out "
        "of the paper's prose and has *not* been verified by a human. Promote a "
        "check by setting a real `path` and `\"confirmed\": true`.",
        "",
        "## Paper",
        "",
        f"**{paper.get('title', '?')}**  ",
        f"{paper.get('authors', '?')}  ",
        f"*{paper.get('journal', '?')}* {paper.get('year', '?')} · "
        f"doi:[{paper.get('doi', '?')}](https://doi.org/{paper.get('doi', '')})"
        + (f" · PMID {paper.get('pmid')}" if paper.get("pmid") else "")
        + (f" · {paper.get('pmcid')}" if paper.get("pmcid") else ""),
        "",
        f"Resolver confidence: **{resolution.get('confidence', 0):.2f}** "
        f"({resolution.get('decision', '?')}).",
        "",
        "## Data sources found in the paper",
        "",
    ]
    acc = harvest.get("accessions") or {}
    if acc:
        L += ["| Repository | Accession | In Data Availability | Mentions |",
              "|---|---|---|---|"]
        for key, entries in sorted(acc.items()):
            for e in sorted(entries, key=lambda x: x["count"], reverse=True)[:6]:
                L.append(f"| {e['repo']} | `{e['value']}` | "
                          f"{'✓' if e['in_data_availability'] else '—'} | {e['count']} |")
    else:
        L.append("_No repository accessions found "
                  + ("(full text was unavailable — see caveats)."
                     if harvest.get("degraded_reason") else "in the full text.")
                  + "_")
    L += ["", "## What IGVFagent does", "",
          f"**Route:** `{spec['name']}` — {spec['title']}  ",
          f"**Skill output dir:** `Docs/{spec['skill_output_dir']}/`  ",
          f"**Modelled on:** {spec.get('reference_benchmark', '—')}", ""]
    if spec.get("notes"):
        L += [f"> {spec['notes']}", ""]
    L += ["```bash", f"bash Benchmarks/{paper_id}/run.sh",
          f".venv/bin/python Benchmarks/concordance.py --benchmark {paper_id}",
          "```", ""]

    L += ["## Concordance", "",
          "_Not yet run. Execute `run.sh` then `concordance.py`, and paste the "
          "paper-claim vs IGVFagent-measured table here._", "",
          "## Honest caveats", ""]
    if harvest.get("degraded_reason"):
        L.append(f"* **Full text unavailable.** {harvest['degraded_reason']} "
                  "Accession and claim extraction therefore saw only the abstract, "
                  "so this scaffold is almost certainly incomplete.")
    if unresolved:
        L.append("* **Unresolved inputs.** " + ", ".join(
            f"`{v['name']}` ({v['prompt']})" for v in unresolved)
            + " — `run.sh` exits 77 until these are supplied.")
    if spec.get("local_input"):
        L.append(f"* **Local data required.** {spec['local_input']['hint']}")
    if spec.get("label_hint"):
        L.append("* **Run directory is not paper-tagged.** This route's CLI "
                  "steps accept no `--label`, so `concordance.py` matches the "
                  f"skill's default `{spec['label_hint']}` directory and will "
                  "score whichever run of this skill is newest — not "
                  "necessarily this paper's.")
    if routing.get("controlled_access_repos"):
        L.append("* **Controlled access.** The paper deposits in "
                  + ", ".join(routing["controlled_access_repos"])
                  + ", which needs an approved application — IGVFagent cannot "
                    "fetch it for you.")
    n_unconf = sum(1 for c in
                    (harvest.get("numeric_claims") or [])[:12]) + len(
                        (harvest.get("llm_claims") or {}).get("claims", [])[:6])
    if n_unconf:
        L.append(f"* **{n_unconf} unconfirmed checks.** Extracted from prose by "
                  "regex and (where enabled) an LLM. They do not count toward a "
                  "pass until a human sets their JSON path and confirms them.")
    L.append("* This scaffold proves the *chain runs and produces artefacts*. It "
             "does not by itself prove the paper's result was reproduced — that "
             "claim needs confirmed checks with real JSON paths.")
    L += ["", "## Provenance", "",
          "`provenance.json` in this directory holds the full resolve / harvest / "
          "route record, including every source URL consulted.", ""]
    return "\n".join(L)


def _render_operations(paper_id: str, spec: Dict[str, Any],
                        variables: List[Dict[str, Any]]) -> str:
    L = [
        f"# {paper_id} — Operations",
        "",
        "Shared prerequisites live in `Benchmarks/OPERATIONS_GUIDE.md`; this file "
        "only covers what is specific to this paper.",
        "",
        "## 1. Inputs",
        "",
        "| Variable | Value | Origin | Override env var |",
        "|---|---|---|---|",
    ]
    prefix = paper_id.split("_")[0].upper()
    for v in variables:
        val = v["value"] if v["value"] is not None else "**TODO_VERIFY**"
        L.append(f"| `{v['name']}` | `{val}` | {v['origin']} | `{prefix}_{v['name']}` |")
    L += ["", "## 2. Run", "",
          "```bash", f"bash Benchmarks/{paper_id}/run.sh",
          f".venv/bin/python Benchmarks/concordance.py --benchmark {paper_id}",
          "```", "",
          "Exit 77 means a required input is missing — the message names which.",
          "", "## 3. Promoting unconfirmed checks", "",
          "`expected.json` ships the paper's prose-derived numbers as "
          "`\"confirmed\": false` with `\"path\": \"TODO_SET_JSON_PATH\"`. To turn "
          "one into a real reproducibility claim:",
          "",
          "1. Run the chain once and open the artefact named by `primary_artefact`.",
          "2. Find the key that holds the comparable quantity; put its dotted path "
          "in `path`.",
          "3. Check the quoted sentence in `provenance.quote` really states that "
          "number for that quantity.",
          "4. Set `\"confirmed\": true` and tighten `min`/`max` to the tolerance "
          "you are willing to defend.",
          "",
          "## 4. Troubleshooting", "",
          "See §5 of `Benchmarks/OPERATIONS_GUIDE.md`.", ""]
    return "\n".join(L)


def scaffold(paper_id: str, resolution: Dict[str, Any], harvest: Dict[str, Any],
              routing: Dict[str, Any], *, route_name: Optional[str] = None,
              force: bool = False) -> Dict[str, Any]:
    name = route_name or routing.get("selected")
    if not name:
        raise SystemExit(
            "No route selected and none could be inferred. "
            + (routing.get("unroutable_reason") or "")
            + " Pass --route <name> explicitly (`igvfagent bench list-routes`)."
        )
    spec = get_route(name)
    target = BENCH_DIR / paper_id
    if target.exists() and not force:
        raise SystemExit(
            f"{target.relative_to(ROOT)} already exists. Pass --force to overwrite."
        )
    target.mkdir(parents=True, exist_ok=True)

    variables = _resolve_vars(spec, harvest, paper_id)

    run_sh = target / "run.sh"
    run_sh.write_text(_render_run_sh(paper_id, spec, variables, resolution))
    run_sh.chmod(0o755)

    expected = _render_expected(paper_id, spec, resolution, harvest, variables)
    _write_json(target / "expected.json", expected)
    (target / "README.md").write_text(
        _render_readme(paper_id, spec, resolution, harvest, variables, routing))
    (target / "OPERATIONS.md").write_text(
        _render_operations(paper_id, spec, variables))
    _write_json(target / "provenance.json", {
        "paper_id": paper_id,
        "resolution": resolution,
        "harvest": harvest,
        "routing": routing,
        "route_selected": name,
        "variables": variables,
    })

    # Register in the generated-benchmark list that run_all.sh picks up.
    existing = []
    if GENERATED_LIST.is_file():
        existing = [ln.strip() for ln in GENERATED_LIST.read_text().splitlines()
                    if ln.strip() and not ln.startswith("#")]
    if paper_id not in existing:
        existing.append(paper_id)
        GENERATED_LIST.write_text(
            "# Benchmarks scaffolded by `igvfagent bench scaffold`.\n"
            "# run_all.sh appends these to its paper list. One paper-id per line.\n"
            + "\n".join(sorted(existing)) + "\n")

    n_unconf = sum(1 for c in expected["checks"] if not c.get("confirmed"))
    return {
        "paper_id": paper_id,
        "dir": str(target.relative_to(ROOT)),
        "route": name,
        "files": ["run.sh", "expected.json", "README.md", "OPERATIONS.md",
                   "provenance.json"],
        "n_checks": len(expected["checks"]),
        "n_unconfirmed": n_unconf,
        "unresolved_vars": [v["name"] for v in variables if v["needs_human"]],
        "requires_local_input": bool(spec.get("local_input")),
    }


# ----------------------------------------------------------------------------
# Gene-symbol extraction (feeds the mavedb / chipatlas routes)
# ----------------------------------------------------------------------------

_GENE_RE = re.compile(r"\b([A-Z][A-Z0-9]{2,9})\b")
_GENE_STOP = {
    "DNA", "RNA", "PCR", "SNP", "SNV", "UTR", "CRISPR", "MPRA", "ATAC", "CHIP",
    "GWAS", "QTL", "EQTL", "AUC", "AUPRC", "FDR", "TSS", "CDS", "WGS", "WES",
    "IGVF", "ENCODE", "GEO", "NCBI", "EMBL", "EBI", "USA", "NIH", "HGNC",
    "CAS9", "GFP", "FACS", "UMAP", "PCA", "TSNE", "SGRNA", "GRNA", "HTO",
    "ONLY", "THIS", "WITH", "FROM", "THAT", "WERE", "HAVE", "MORE", "SUCH",
    "BOTH", "EACH", "THAN", "THEN", "WHEN", "ALSO", "INTO", "OVER", "ALL",
    "SGE", "DMS", "MAVE", "VAMP", "ABC", "E2G", "TAD", "CTCF",
    # Assay/method names that look like gene symbols in title case.
    "EFFECTS", "SEQ", "FISH", "MPRA", "STARR", "SHARE", "SPLIT", "HTO",
    "NGS", "TSV", "CSV", "JSON", "HTTP", "HTTPS", "FASTQ", "BAM", "VCF",
}


def extract_genes(rec: Dict[str, Any], text: str) -> List[str]:
    """Gene symbols, ranked by prominence in title > abstract > body."""
    scores: Dict[str, float] = {}
    for field, weight in ((rec.get("title", ""), 10.0),
                           (rec.get("abstract", ""), 3.0)):
        for m in _GENE_RE.finditer(field or ""):
            g = m.group(1)
            if g in _GENE_STOP or g.isdigit():
                continue
            scores[g] = scores.get(g, 0.0) + weight
    for m in _GENE_RE.finditer(text[:200000]):
        g = m.group(1)
        if g in _GENE_STOP or g.isdigit():
            continue
        scores[g] = scores.get(g, 0.0) + 0.05
    return [g for g, _ in sorted(scores.items(), key=lambda kv: kv[1], reverse=True)[:12]]


# ----------------------------------------------------------------------------
# Subcommand implementations
# ----------------------------------------------------------------------------

def do_resolve(args: argparse.Namespace) -> Dict[str, Any]:
    want: Dict[str, Any] = {
        "title": args.title, "author": args.author,
        "journal": args.journal, "year": args.year,
    }
    ident = {"kind": "free_text", "value": "", "raw": ""}
    query = args.query or args.doi or args.pmid or args.pmcid or args.url or args.title
    if not query:
        raise SystemExit(
            "Nothing to resolve. Pass --query (title / URL / DOI / PMID) or any "
            "of --doi / --pmid / --pmcid / --url / --title."
        )
    if args.doi:
        ident = {"kind": "doi", "value": args.doi, "raw": args.doi}
    elif args.pmid:
        ident = {"kind": "pmid", "value": str(args.pmid), "raw": str(args.pmid)}
    elif args.pmcid:
        ident = {"kind": "pmcid", "value": args.pmcid, "raw": args.pmcid}
    elif args.url:
        ident = classify_identifier(args.url)
    else:
        ident = classify_identifier(query)
        if ident["kind"] == "free_text" and not want["title"]:
            want["title"] = ident["value"]

    logging.info("Identifier classified as %s: %s", ident["kind"], ident["value"])

    candidates: List[Dict[str, Any]] = []
    consulted: List[str] = []

    if ident["kind"] == "doi":
        doi = ident["value"]
        cr = crossref_work(doi)
        consulted.append("crossref")
        epmc = europepmc_query(f'DOI:"{doi}"', limit=5)
        consulted.append("europepmc")
        brx = biorxiv_details(doi) if (cr or {}).get("is_preprint") or not cr else None
        if brx:
            consulted.append(brx["source"])
        merged = _merge_records(cr, epmc[0] if epmc else None, brx)
        if merged:
            merged.setdefault("doi", doi)
            candidates.append(merged)
        # A preprint that has since been published: surface both.
        pub_doi = (brx or {}).get("published_doi")
        if pub_doi and pub_doi != doi:
            pub = crossref_work(pub_doi)
            if pub:
                pub["note"] = "published version of the preprint you supplied"
                candidates.append(pub)

    elif ident["kind"] in ("pmid", "pmcid"):
        if ident["kind"] == "pmid":
            rec = pubmed_by_id(ident["value"])
        else:
            hits = europepmc_query(f'PMCID:"{ident["value"]}"', limit=3)
            rec = hits[0] if hits else None
        consulted.append("europepmc/pubmed")
        if rec and rec.get("doi"):
            cr = crossref_work(rec["doi"])
            consulted.append("crossref")
            rec = _merge_records(rec, cr)
        if rec:
            candidates.append(rec)

    else:
        q = want["title"] or ident["value"]
        parts = [q]
        if want["author"]:
            parts.append(want["author"])
        joint = " ".join(parts)
        try:
            candidates += europepmc_query(joint, limit=args.limit)
            consulted.append("europepmc")
        except Exception as e:
            logging.warning("europepmc search failed: %s", e)
        for fn, label in ((_ref.pubmed_search, "pubmed"),
                           (_ref.openalex_search, "openalex"),
                           (_ref.semanticscholar_search, "semanticscholar")):
            try:
                candidates += fn(joint, limit=args.limit)
                consulted.append(label)
            except Exception as e:
                logging.warning("%s search failed: %s", label, e)
        try:
            candidates += _ref.biorxiv_search(joint, limit=args.limit)
            consulted.append("biorxiv")
        except Exception as e:
            logging.warning("biorxiv search failed: %s", e)
        candidates = _ref.dedup_records(candidates)

    if not candidates:
        return {"decision": "not_found", "identifier": ident,
                 "sources_consulted": consulted, "confidence": 0.0,
                 "candidates": [],
                 "reason": "No source returned a record for this input."}

    # An exact identifier is its own evidence; free text has to be scored.
    scored: List[Dict[str, Any]] = []
    for c in candidates:
        if ident["kind"] in ("doi", "pmid", "pmcid") and not any(want.values()):
            s, reasons = 1.0, [f"exact {ident['kind']} lookup"]
        else:
            s, reasons = score_candidate(c, want)
            if ident["kind"] in ("doi", "pmid", "pmcid"):
                s = max(s, 0.9)
                reasons.append(f"exact {ident['kind']} lookup")
        scored.append({"record": c, "score": round(s, 3), "reasons": reasons})
    scored.sort(key=lambda d: d["score"], reverse=True)

    best = scored[0]
    margin = best["score"] - (scored[1]["score"] if len(scored) > 1 else 0.0)
    if best["score"] >= 0.85 and (len(scored) == 1 or margin >= 0.15):
        decision = "resolved"
    elif best["score"] >= 0.55:
        decision = "ambiguous"
    else:
        decision = "low_confidence"

    paper = best["record"]
    paper_id = args.paper_id or derive_paper_id(paper)
    out = {
        "decision": decision,
        "confidence": best["score"],
        "margin_over_runner_up": round(margin, 3),
        "identifier": ident,
        "constraints": {k: v for k, v in want.items() if v},
        "sources_consulted": sorted(set(consulted)),
        "paper_id": paper_id,
        "paper": paper,
        "match_reasons": best["reasons"],
        "candidates": [
            {"score": c["score"], "title": c["record"].get("title"),
             "authors": (c["record"].get("authors") or "")[:120],
             "journal": c["record"].get("journal"), "year": c["record"].get("year"),
             "doi": c["record"].get("doi"), "pmid": c["record"].get("pmid"),
             "reasons": c["reasons"]}
            for c in scored[:args.top]
        ],
    }
    return out


def do_harvest(args: argparse.Namespace, resolution: Dict[str, Any]) -> Dict[str, Any]:
    rec = resolution.get("paper") or {}
    ft = fetch_fulltext(rec)
    sections = ft["sections"]
    text = ft["text"]

    harvest: Dict[str, Any] = {
        "paper_id": resolution.get("paper_id"),
        "doi": rec.get("doi"),
        "title": rec.get("title"),
        # Crossref abstracts arrive with JATS indentation intact; collapse it so
        # downstream matching and the LLM prompt see clean prose.
        "abstract": " ".join((rec.get("abstract") or "").split()),
        "fulltext_source": ft["fulltext_source"],
        "degraded_reason": ft["degraded_reason"],
        "n_sections": len(sections),
        "n_chars": len(text),
        "section_titles": [s.get("title", "") for s in sections],
        "data_availability": extract_data_availability(sections),
        "accessions": extract_accessions(text, sections),
        "assays": extract_assays(text),
        "genes": extract_genes(rec, text),
        "numeric_claims": extract_numeric_claims(sections),
        "llm_claims": {"enabled": False, "claims": [],
                        "error": "disabled with --no-llm"},
    }
    if not args.no_llm:
        harvest["llm_claims"] = llm_draft_claims(
            harvest, backend=args.backend, model=args.model)
    return harvest


def do_run(paper_id: str) -> Dict[str, Any]:
    script = BENCH_DIR / paper_id / "run.sh"
    if not script.is_file():
        raise SystemExit(f"No {script.relative_to(ROOT)} — run `bench scaffold` first.")
    print(f"== running {script.relative_to(ROOT)} ==")
    proc = subprocess.run(["bash", str(script)], cwd=str(ROOT))
    status = {0: "ok", 77: "skipped_missing_local_input"}.get(
        proc.returncode, "failed")
    return {"paper_id": paper_id, "returncode": proc.returncode, "status": status}


def do_score(paper_id: str) -> Dict[str, Any]:
    scorer = BENCH_DIR / "concordance.py"
    py = ROOT / ".venv" / "bin" / "python"
    exe = str(py) if py.is_file() else sys.executable
    proc = subprocess.run([exe, str(scorer), "--benchmark", paper_id],
                           cwd=str(ROOT), capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
    results = sorted((BENCH_DIR / "results").glob("*_concordance.json"),
                      key=lambda p: p.name, reverse=True)
    payload = _read_json(results[0]) if results else []
    mine = next((r for r in payload if r.get("paper") == paper_id), None)
    return {"paper_id": paper_id, "returncode": proc.returncode,
             "concordance": mine, "report": str(results[0]) if results else None}


def _render_report(paper_id: str, expected: Dict[str, Any],
                    score: Dict[str, Any]) -> str:
    conc = score.get("concordance") or {}
    checks = {c.get("name"): c for c in conc.get("checks", [])}
    L = [f"# Replication report — {paper_id}", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M:%S %Z')} by "
         "`igvfagent bench report`.", ""]
    p = expected.get("paper", {})
    L += [f"**{p.get('title', '?')}**  ",
          f"*{p.get('journal', '?')}* {p.get('year', '?')} · doi:{p.get('doi', '?')}",
          "", f"Route: `{(expected.get('data_source') or {}).get('type', '?')}` · "
              f"artefacts under `Docs/{expected.get('skill_output_dir')}/`", ""]

    conf = [c for c in expected.get("checks", []) if c.get("confirmed")]
    unconf = [c for c in expected.get("checks", []) if not c.get("confirmed")]

    L += ["## Confirmed checks", ""]
    if conf:
        L += ["| Check | Expected | Measured | Result |", "|---|---|---|---|"]
        for c in conf:
            r = checks.get(c.get("name")) or {}
            exp = c.get("expected", f"[{c.get('min')}, {c.get('max')}]"
                        if c.get("type") == "range" else c.get("type"))
            L.append(f"| {c.get('name')} | {exp} | {r.get('detail', '—')} | "
                      f"{'✓' if r.get('passed') else '✗'} |")
    else:
        L.append("_None._")

    L += ["", "## Unconfirmed checks (paper-derived, not human-verified)", ""]
    if unconf:
        L += ["| Paper claim | Source | Quote |", "|---|---|---|"]
        for c in unconf:
            prov = c.get("provenance") or {}
            quote = (prov.get("quote") or "")[:160].replace("|", r"\|")
            L.append(f"| {c.get('name')} | {prov.get('kind', '?')} | {quote} |")
        L += ["", "These do **not** count toward the pass tally. Promote one by "
                  "setting its `path` to a real key in the run artefact and "
                  "flipping `confirmed` to `true`.", ""]
    else:
        L.append("_None._")

    L += ["", "## Verdict", "",
          f"- Status: **{conc.get('status', 'not scored')}**",
          f"- Confirmed checks passed: **{conc.get('n_passed', 0)} / "
          f"{conc.get('n_total', 0)}**",
          f"- Unconfirmed checks awaiting review: **{len(unconf)}**", ""]
    if conc.get("status") == "ok" and not conf:
        L.append("> No confirmed checks exist, so a green status here means only "
                  "that the chain ran — not that the paper was reproduced.")
    return "\n".join(L)


# ----------------------------------------------------------------------------
# selftest — validate the resolver + router against the committed benchmarks
# ----------------------------------------------------------------------------

def do_selftest(args: argparse.Namespace) -> Dict[str, Any]:
    """Re-derive committed benchmarks from their DOIs and score the machinery.

    Two objective questions, both answerable without a human:

    * **Resolver:** take the DOI from a committed ``expected.json``, get the
      canonical title from Crossref, feed *only that title* back through free-text
      resolution, and check the DOI round-trips.
    * **Router:** harvest the paper and check the top route's
      ``skill_output_dir`` matches the one the committed benchmark uses.
    """
    cases = []
    for exp_path in sorted(BENCH_DIR.glob("*/expected.json")):
        spec = _read_json(exp_path)
        paper = spec.get("paper") or {}
        doi = (paper.get("doi") or "").split(";")[0].strip()
        if not doi:
            continue
        cases.append({
            "paper_id": exp_path.parent.name,
            "doi": doi,
            "skill_output_dir": spec.get("skill_output_dir"),
        })
    if args.limit:
        cases = cases[: args.limit]

    results = []
    for case in cases:
        row: Dict[str, Any] = dict(case)
        cr = crossref_work(case["doi"])
        if not cr or not cr.get("title"):
            row["resolver"] = "no_crossref_record"
            results.append(row)
            continue
        row["title"] = cr["title"]

        ns = argparse.Namespace(
            query=cr["title"], doi=None, pmid=None, pmcid=None, url=None,
            title=cr["title"], author=None, journal=None, year=None,
            limit=10, top=5, paper_id=None,
        )
        try:
            res = do_resolve(ns)
        except SystemExit as e:
            row["resolver"] = f"error: {e}"
            results.append(row)
            continue
        got = (res.get("paper") or {}).get("doi", "") or ""
        row["resolved_doi"] = got
        row["resolver_confidence"] = res.get("confidence")
        row["resolver"] = "hit" if got.lower() == case["doi"].lower() else (
            "hit_in_candidates" if any(
                (c.get("doi") or "").lower() == case["doi"].lower()
                for c in res.get("candidates", [])) else "miss")

        if args.with_router:
            hv = do_harvest(argparse.Namespace(no_llm=True, backend=None, model=None),
                             res)
            rt = route(hv)
            row["routed"] = rt.get("selected")
            top_dirs = [p["skill_output_dir"] for p in rt.get("plans", [])]
            row["routed_skill_dir"] = top_dirs[0] if top_dirs else None
            if case["skill_output_dir"] is None:
                row["router"] = "no_ground_truth"
            elif row["routed_skill_dir"] == case["skill_output_dir"]:
                row["router"] = "hit"
            elif case["skill_output_dir"] in top_dirs:
                row["router"] = "hit_in_top3"
            else:
                row["router"] = "miss"
            row["fulltext"] = hv.get("fulltext_source") or "abstract-only"
        results.append(row)

    def _tally(key):
        vals = [r.get(key) for r in results if r.get(key)]
        return {v: vals.count(v) for v in sorted(set(vals))}

    return {"n_cases": len(results),
             "resolver_tally": _tally("resolver"),
             "router_tally": _tally("router") if args.with_router else None,
             "cases": results}


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def _print_resolution(res: Dict[str, Any]) -> None:
    print()
    print(f"Decision:   {res['decision']}  (confidence {res.get('confidence', 0):.2f})")
    if res["decision"] == "not_found":
        print(f"Reason:     {res.get('reason')}")
        return
    p = res.get("paper") or {}
    print(f"Paper id:   {res.get('paper_id')}")
    print(f"Title:      {p.get('title')}")
    print(f"Authors:    {(p.get('authors') or '')[:110]}")
    print(f"Journal:    {p.get('journal')} {p.get('year')}")
    print(f"DOI:        {p.get('doi')}"
          + (f"   PMID {p.get('pmid')}" if p.get("pmid") else "")
          + (f"   {p.get('pmcid')}" if p.get("pmcid") else ""))
    print(f"Why:        {'; '.join(res.get('match_reasons', []))}")
    if res["decision"] != "resolved" and len(res.get("candidates", [])) > 1:
        print()
        print("Ambiguous — candidates:")
        for i, c in enumerate(res["candidates"], 1):
            print(f"  {i}. [{c['score']:.2f}] {(c.get('title') or '')[:88]}")
            print(f"       {c.get('journal')} {c.get('year')} · {c.get('doi')}")
        print()
        print("Disambiguate with --doi / --pmid, or add --author / --journal / --year.")


def _print_harvest(hv: Dict[str, Any]) -> None:
    print()
    print(f"Full text:  {hv.get('fulltext_source') or 'UNAVAILABLE'}")
    if hv.get("degraded_reason"):
        print(f"  ! {hv['degraded_reason']}")
    print(f"Sections:   {hv.get('n_sections')}  ({hv.get('n_chars', 0):,} chars)")
    da = hv.get("data_availability") or []
    print(f"Data availability statements: {len(da)}")
    for d in da[:2]:
        print(f"  [{d['title']}] {d['text'][:220]}…")
    acc = hv.get("accessions") or {}
    n = sum(len(v) for v in acc.values())
    print(f"Accessions: {n}")
    for key, entries in sorted(acc.items()):
        marks = ", ".join(
            e["value"] + ("*" if e["in_data_availability"] else "")
            for e in entries[:6])
        print(f"  {key:15s} {marks}")
    if n:
        print("  (* = named in a Data/Code Availability section)")
    print("Assays:     " + ", ".join(
        f"{a['assay']}×{a['mentions']}" for a in (hv.get("assays") or [])[:8]))
    print("Genes:      " + ", ".join((hv.get("genes") or [])[:10]))
    print(f"Numeric claim candidates: {len(hv.get('numeric_claims') or [])} "
          f"(all unconfirmed)")
    llm = hv.get("llm_claims") or {}
    if llm.get("claims"):
        print(f"LLM-proposed headline claims: {len(llm['claims'])}")
        for c in llm["claims"]:
            flag = "" if c.get("quote_grounded_in_source") else "  [quote NOT found in source]"
            print(f"  - {c.get('value')} {c.get('unit')}: {c.get('description')}{flag}")
    elif llm.get("error"):
        print(f"LLM draft: skipped ({llm['error']})")


def _print_routing(rt: Dict[str, Any]) -> None:
    print()
    if not rt.get("plans"):
        print("No route matched.")
        print(f"  {rt.get('unroutable_reason')}")
        return
    print("Routes (best first):")
    for p in rt["plans"]:
        flag = "  [needs local data]" if p["requires_local_input"] else ""
        print(f"  [{p['score']:5.2f}] {p['route']:20s} {p['title']}{flag}")
        print(f"           evidence: {', '.join(p['evidence'])}")
        print(f"           modelled on: {p.get('reference_benchmark')}")
    if rt.get("controlled_access_repos"):
        print(f"  ! controlled-access deposits: "
              f"{', '.join(rt['controlled_access_repos'])}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Resolve a publication, then scaffold / run / score a "
                    "reproduction of it with IGVFagent's analysis skills.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def _ident_args(s):
        s.add_argument("--query", help="Title, URL, DOI, PMID, or free text.")
        s.add_argument("--doi")
        s.add_argument("--pmid")
        s.add_argument("--pmcid")
        s.add_argument("--url")
        s.add_argument("--title")
        s.add_argument("--author", help="Surname(s), comma-separated.")
        s.add_argument("--journal")
        s.add_argument("--year", type=int)
        s.add_argument("--limit", type=int, default=15,
                        help="Records to pull per source for free-text search.")
        s.add_argument("--top", type=int, default=5,
                        help="Candidates to show when ambiguous.")
        s.add_argument("--paper-id", help="Override the derived benchmark id.")

    def _llm_args(s):
        s.add_argument("--no-llm", action="store_true",
                        help="Deterministic extraction only; skip the LLM "
                             "headline-claim draft.")
        s.add_argument("--backend", help="LLM backend override.")
        s.add_argument("--model", help="LLM model override.")

    s = sub.add_parser("resolve", help="Identifier or free text -> one paper.")
    _ident_args(s)

    s = sub.add_parser("harvest", help="Full text -> accessions + claims.")
    s.add_argument("--paper-id", required=True)
    s.add_argument("--resolution-json", help="Explicit resolution.json path.")
    _llm_args(s)

    s = sub.add_parser("route", help="Accessions + assays -> IGVFagent chain.")
    s.add_argument("--paper-id", required=True)
    s.add_argument("--harvest-json")
    s.add_argument("--top", type=int, default=3)

    s = sub.add_parser("scaffold", help="Write Benchmarks/<paper-id>/.")
    s.add_argument("--paper-id", required=True)
    s.add_argument("--resolution-json")
    s.add_argument("--harvest-json")
    s.add_argument("--routing-json")
    s.add_argument("--route", help="Force a specific route name.")
    s.add_argument("--force", action="store_true",
                    help="Overwrite an existing benchmark directory.")

    s = sub.add_parser("run", help="Execute Benchmarks/<paper-id>/run.sh.")
    s.add_argument("--paper-id", required=True)

    s = sub.add_parser("score", help="Score with Benchmarks/concordance.py.")
    s.add_argument("--paper-id", required=True)

    s = sub.add_parser("report", help="Render the replication report.")
    s.add_argument("--paper-id", required=True)

    s = sub.add_parser("pipeline",
                        help="resolve -> harvest -> route -> scaffold "
                             "(add --execute for run -> score -> report).")
    _ident_args(s)
    _llm_args(s)
    s.add_argument("--route", help="Force a specific route name.")
    s.add_argument("--force", action="store_true")
    s.add_argument("--execute", action="store_true",
                    help="Also run, score, and report. Off by default: review "
                         "the scaffold before trusting it.")

    s = sub.add_parser("selftest",
                        help="Re-derive the committed benchmarks from their DOIs "
                             "and score the resolver (and optionally the router).")
    s.add_argument("--limit", type=int, default=0,
                    help="Only test the first N benchmarks.")
    s.add_argument("--with-router", action="store_true",
                    help="Also fetch full text and score routing (slower).")

    sub.add_parser("list-routes", help="Show the routing table.")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    log_path = setup_logging()
    _ref.mkdirs()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    logging.info("Log file: %s", log_path)

    ts = timestamp()

    if args.cmd == "list-routes":
        print(f"{len(ROUTES)} routes\n")
        for spec in ROUTES:
            print(f"{spec['name']:22s} {spec['title']}")
            print(f"{'':22s} skill dir: Docs/{spec['skill_output_dir']}/ · "
                  f"modelled on {spec.get('reference_benchmark', '—')}")
            m = spec["match"]
            print(f"{'':22s} matches: "
                  f"accessions={m.get('accessions') or '—'} "
                  f"assays={m.get('assays') or '—'}")
            if spec.get("local_input"):
                print(f"{'':22s} needs local data: {spec['local_input']['hint']}")
            print()
        return 0

    if args.cmd == "selftest":
        out = do_selftest(args)
        print()
        print(f"Cases: {out['n_cases']}")
        print(f"Resolver: {out['resolver_tally']}")
        if out.get("router_tally"):
            print(f"Router:   {out['router_tally']}")
        print()
        for r in out["cases"]:
            mark = {"hit": "✓", "hit_in_candidates": "~",
                     "hit_in_top3": "~"}.get(r.get("resolver"), "✗")
            line = f"  {mark} {r['paper_id']:32s} resolver={r.get('resolver')}"
            if r.get("router"):
                line += f"  router={r.get('router')} ({r.get('routed')})"
            print(line)
        out_dir = REPORT_DIR / f"{ts}_selftest"
        _write_json(out_dir / "selftest.json", out)
        print(f"\nWrote {(out_dir / 'selftest.json').relative_to(ROOT)}")
        return 0

    if args.cmd in ("resolve", "pipeline"):
        res = do_resolve(args)
        _print_resolution(res)
        paper_id = res.get("paper_id") or "unresolved"
        out_dir = REPORT_DIR / f"{ts}_{paper_id}"
        _write_json(out_dir / "resolution.json", res)
        print(f"\nWrote {(out_dir / 'resolution.json').relative_to(ROOT)}")
        if res["decision"] == "not_found":
            return 1
        if args.cmd == "resolve":
            if res["decision"] != "resolved":
                print("\nNot a clean match — confirm before scaffolding.")
                return 3
            print(f"\nNext: igvfagent bench harvest --paper-id {paper_id}")
            return 0
    else:
        paper_id = args.paper_id

    if args.cmd == "harvest" or args.cmd == "pipeline":
        if args.cmd == "harvest":
            res = _load_stage(paper_id, "resolution", getattr(args, "resolution_json", None))
            out_dir = REPORT_DIR / f"{ts}_{paper_id}"
        hv = do_harvest(args, res)
        _print_harvest(hv)
        _write_json(out_dir / "harvest.json", hv)
        print(f"\nWrote {(out_dir / 'harvest.json').relative_to(ROOT)}")
        if args.cmd == "harvest":
            print(f"\nNext: igvfagent bench route --paper-id {paper_id}")
            return 0

    if args.cmd == "route" or args.cmd == "pipeline":
        if args.cmd == "route":
            res = _load_stage(paper_id, "resolution")
            hv = _load_stage(paper_id, "harvest", getattr(args, "harvest_json", None))
            out_dir = REPORT_DIR / f"{ts}_{paper_id}"
        rt = route(hv, top=getattr(args, "top", 3))
        _print_routing(rt)
        _write_json(out_dir / "routing.json", rt)
        print(f"\nWrote {(out_dir / 'routing.json').relative_to(ROOT)}")
        if args.cmd == "route":
            if not rt.get("plans"):
                return 4
            print(f"\nNext: igvfagent bench scaffold --paper-id {paper_id}")
            return 0

    if args.cmd == "scaffold" or args.cmd == "pipeline":
        if args.cmd == "scaffold":
            res = _load_stage(paper_id, "resolution", args.resolution_json)
            hv = _load_stage(paper_id, "harvest", args.harvest_json)
            rt = (_read_json(Path(args.routing_json)) if args.routing_json
                   else _load_stage(paper_id, "routing"))
        sc = scaffold(paper_id, res, hv, rt, route_name=args.route, force=args.force)
        print()
        print(f"Scaffolded {sc['dir']}  (route: {sc['route']})")
        print(f"  files:        {', '.join(sc['files'])}")
        print(f"  checks:       {sc['n_checks']} "
              f"({sc['n_unconfirmed']} unconfirmed, awaiting human review)")
        if sc["unresolved_vars"]:
            print(f"  TODO_VERIFY:  {', '.join(sc['unresolved_vars'])} "
                  f"— run.sh exits 77 until these are set")
        if sc["requires_local_input"]:
            print("  local data:   required (see README caveats)")
        print()
        print("Review the scaffold, then:")
        print(f"  bash Benchmarks/{paper_id}/run.sh")
        print(f"  .venv/bin/python Benchmarks/concordance.py --benchmark {paper_id}")
        if args.cmd == "scaffold" or not args.execute:
            return 0

    if args.cmd == "run" or (args.cmd == "pipeline" and args.execute):
        rr = do_run(paper_id)
        print(f"\nrun.sh -> {rr['status']} (rc={rr['returncode']})")
        if args.cmd == "run":
            return 0 if rr["status"] != "failed" else 1

    if args.cmd == "score" or (args.cmd == "pipeline" and args.execute):
        sr = do_score(paper_id)
        if args.cmd == "score":
            return 0

    if args.cmd == "report" or (args.cmd == "pipeline" and args.execute):
        if args.cmd == "report":
            sr = do_score(paper_id)
        expected = _read_json(BENCH_DIR / paper_id / "expected.json")
        md = _render_report(paper_id, expected, sr)
        out_dir = REPORT_DIR / f"{ts}_{paper_id}"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "replication_report.md").write_text(md)
        print()
        print(md)
        print(f"\nWrote {(out_dir / 'replication_report.md').relative_to(ROOT)}")
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
