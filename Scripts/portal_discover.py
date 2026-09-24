#!/usr/bin/env python3
"""Find IGVF Portal data by topic: phenotype, tissue, gene, prediction type, library.

`processed lineage` answers "what is linked to this accession?". This answers
the question before it: "which accessions are about coronary artery disease /
heart / GATA1 / element-gene links?" It uses the link fields the lab
submission diagrams (Data/IGVF/DataModels, summarised in
Docs/Architecture/IGVF_PORTAL_DATA_MODEL.md) put on each object type:

    PredictionSet   associated_phenotypes, assessed_genes, file_set_type,
                    samples.sample_terms (virtual samples = model context)
    ModelSet        file_set_type, assessed_genes, samples
    MeasurementSet  samples.sample_terms, targeted_genes,
                    construct_library_sets.file_set_type, preferred_assay_titles
    AnalysisSet     samples.sample_terms, file_set_type
    CuratedSet      file_set_type (variants, elements, guide RNAs, training
                    data for predictive models, external sequencing data ...)
    ConstructLibrarySet  file_set_type (guide / reporter / editing template
                    library), small_scale_gene_list, associated_phenotypes

Words are matched to the values the Portal actually uses, not guessed: for each
filter the facet of that field is read (and, when the Portal caps a facet at
~100 terms, counted from the items), the closest values are chosen by shared
words and stems plus a small ontology synonym map (cardiac ~ heart), and the
search then filters on those exact values. A search that matches nothing is a
result, reported with the nearest values, not an error.

    igvfagent processed discover --phenotype "coronary artery disease"
    igvfagent processed discover --tissue heart --types MeasurementSet,AnalysisSet
    igvfagent processed discover --gene GATA1
    igvfagent processed discover --prediction-type "element-gene links" --tissue liver
"""
from __future__ import annotations

import json
import re
import urllib.parse
from typing import Any, Callable, Dict, List, Optional, Tuple

try:
    from igvfagent import portal_query_skill as pq  # type: ignore
except Exception:  # direct execution
    import portal_query_skill as pq  # type: ignore

Fetch = Callable[[str], "Tuple[int, Any]"]

DEFAULT_TYPES = ("PredictionSet", "ModelSet", "MeasurementSet", "AnalysisSet", "CuratedSet",
                 "ConstructLibrarySet")

# criterion -> {type: facet field}. A type without the field is skipped for
# that criterion (so --phenotype does not search MeasurementSets).
FIELDS: "Dict[str, Dict[str, str]]" = {
    "phenotype": {"PredictionSet": "associated_phenotypes.term_name",
                  "ConstructLibrarySet": "associated_phenotypes.term_name"},
    "tissue": {t: "samples.sample_terms.term_name" for t in
               ("PredictionSet", "ModelSet", "MeasurementSet", "AnalysisSet", "CuratedSet")},
    "gene": {"PredictionSet": "assessed_genes.symbol", "ModelSet": "assessed_genes.symbol",
             "MeasurementSet": "targeted_genes.symbol",
             "ConstructLibrarySet": "small_scale_gene_list.symbol"},
    "prediction_type": {"PredictionSet": "file_set_type", "ModelSet": "file_set_type"},
    "library_type": {"MeasurementSet": "construct_library_sets.file_set_type",
                     "ConstructLibrarySet": "file_set_type"},
    "assay": {"MeasurementSet": "preferred_assay_titles", "AnalysisSet": "assay_titles"},
    "curated_type": {"CuratedSet": "file_set_type"},
}
RESULT_FIELDS = ("accession", "@id", "@type", "file_set_type", "summary", "status", "lab.title",
                 "associated_phenotypes.term_name", "samples.sample_terms.term_name",
                 "preferred_assay_titles", "superseded_by", "scope", "model_name")


def _live() -> Fetch:
    try:
        from igvfagent.raw_data_pipeline import portal_json  # type: ignore
    except Exception:
        from raw_data_pipeline import portal_json  # type: ignore
    return portal_json


def _search(fetch: Fetch, params: "List[Tuple[str, str]]") -> dict:
    status, data = fetch("/search/?" + urllib.parse.urlencode(params + [("format", "json")], doseq=True))
    if status == 200 and isinstance(data, dict):
        return data
    # the Portal answers an empty search with 404; with portal_json the body is dropped
    return {"@graph": [], "total": 0, "facets": [], "_status": status}


def facet_values(fetch: Fetch, typ: str, field: str, base: "List[Tuple[str, str]]") -> "Dict[str, int]":
    """Every value of `field` for `typ` under `base` filters, with counts."""
    data = _search(fetch, [("type", typ), *base, ("limit", "0")])
    facet = next((f for f in data.get("facets") or [] if f.get("field") == field), None)
    terms = {str(t.get("key")): int(t.get("doc_count") or 0) for t in (facet or {}).get("terms") or []}
    if facet is None or len(terms) >= 90:
        # capped or absent facet: count from the items themselves
        items = _search(fetch, [("type", typ), *base, ("field", field), ("limit", "10000")]).get("@graph") or []
        counts: "Dict[str, int]" = {}
        for it in items:
            level: "List[Any]" = [it]
            for part in field.split("."):
                nxt: "List[Any]" = []
                for x in level:
                    y = x.get(part) if isinstance(x, dict) else None
                    nxt.extend(y if isinstance(y, list) else [y] if y is not None else [])
                level = nxt
            for v in {str(x) for x in level if not isinstance(x, (dict, list))}:
                counts[v] = counts.get(v, 0) + 1
        if counts:
            terms = counts
    return terms


# Enumerated types: a near miss is a different kind of data, so only an exact
# value or a whole-phrase containment counts ("element-gene links" must not
# fall back to "spatial gene expression variability" on the word "gene").
CATEGORICAL = {"prediction_type", "library_type", "curated_type", "assay"}


def match_terms(text: str, values: "Dict[str, int]", top: int = 6, strict: bool = False) -> "List[str]":
    """Portal values closest to what the user typed (exact first)."""
    low = text.strip().lower()
    exact = [v for v in values if v.lower() == low]
    if exact and strict:
        return exact
    if exact:
        # "heart" should also reach "heart left ventricle"; whole-word containment only
        word = re.compile(r"(?<![a-z0-9])" + re.escape(low) + r"(?![a-z0-9])")
        wider = sorted((v for v in values if v not in exact and word.search(v.lower())),
                       key=lambda v: -values[v])
        return (exact + wider)[:max(top, 12)]
    if strict:
        return [v for v in values if low in v.lower() or v.lower() in low][:top]
    want = pq._query_tokens(text)
    if not want:                                  # short words (e.g. a gene symbol)
        return [v for v in values if low in v.lower()][:top]
    scored = sorted(((pq._closeness(want, v), values[v], v) for v in values), reverse=True)
    best = scored[0][0] if scored else 0
    # keep only terms nearly as close as the best one: "coronary artery disease"
    # should reach "coronary artery disorder", not every term sharing "disease"
    return [v for sc, _n, v in scored if sc > 0 and sc >= 0.75 * best][:top]


def discover(criteria: "Dict[str, str]", types: "Optional[List[str]]" = None, limit: int = 25,
             fetch: "Optional[Fetch]" = None, include_superseded: bool = False) -> dict:
    fetch = fetch or _live()
    types = list(types or DEFAULT_TYPES)
    out: dict = {"criteria": criteria, "types": {}, "skipped": {}}
    for typ in types:
        applicable = {c: FIELDS[c][typ] for c in criteria if typ in FIELDS.get(c, {})}
        if len(applicable) < len(criteria):
            out["skipped"][typ] = sorted(set(criteria) - set(applicable))
            continue
        base: "List[Tuple[str, str]]" = []
        matched: "Dict[str, List[str]]" = {}
        near: "Dict[str, List[str]]" = {}
        for crit, field in applicable.items():
            vals = facet_values(fetch, typ, field, base)
            m = match_terms(criteria[crit], vals, strict=crit in CATEGORICAL)
            if not m:
                near[crit] = [v for v, _ in sorted(vals.items(), key=lambda kv: -kv[1])[:8]]
                break
            matched[crit] = m
            base += [(field, v) for v in m]       # repeated field = OR within, AND across fields
        if near:
            out["types"][typ] = {"total": 0, "matched": matched, "no_match": near, "items": []}
            continue
        params = [("type", typ), *base, *[("field", f) for f in RESULT_FIELDS], ("limit", str(limit)),
                  ("sort", "-release_timestamp")]
        data = _search(fetch, params)
        items = []
        for it in data.get("@graph") or []:
            if it.get("superseded_by") and not include_superseded:
                continue
            items.append({
                "accession": it.get("accession"), "type": typ, "file_set_type": it.get("file_set_type") or "",
                "summary": (it.get("summary") or "")[:200], "status": it.get("status"),
                "lab": (it.get("lab") or {}).get("title") if isinstance(it.get("lab"), dict) else "",
                "phenotypes": sorted({p.get("term_name") for p in it.get("associated_phenotypes") or []
                                      if isinstance(p, dict) and p.get("term_name")}),
                "sample_terms": sorted({t.get("term_name") for s in it.get("samples") or [] if isinstance(s, dict)
                                        for t in s.get("sample_terms") or [] if isinstance(t, dict)
                                        and t.get("term_name")})[:6],
                "assays": it.get("preferred_assay_titles") or [], "scope": it.get("scope") or "",
                "model_name": it.get("model_name") or ""})
        out["types"][typ] = {"total": int(data.get("total") or 0), "matched": matched, "items": items}
    return out


def to_markdown(res: dict) -> str:
    crit = ", ".join(f"{k.replace('_', ' ')} = {v!r}" for k, v in res["criteria"].items())
    lines = [f"# IGVF Portal data for {crit}", ""]
    total = sum(t["total"] for t in res["types"].values())
    lines += [f"{total} matching object(s) across {len(res['types'])} type(s). Each accession can be walked "
              "with `igvfagent processed lineage <accession>` to reach its processed files, models, "
              "training data, samples and QC.", ""]
    for typ, t in res["types"].items():
        lines += [f"## {typ} ({t['total']})", ""]
        if t.get("matched"):
            lines += ["Matched Portal terms: " + "; ".join(f"{k}: {', '.join(v)}" for k, v in t["matched"].items()), ""]
        if t.get("no_match"):
            for k, v in t["no_match"].items():
                lines += [f"No {k} value matched. Most common values here: {', '.join(v)}", ""]
            continue
        if not t["items"]:
            lines += ["(none)", ""]
            continue
        lines += ["| accession | kind | summary | phenotypes / samples | lab |", "|---|---|---|---|---|"]
        for it in t["items"]:
            ctx = "; ".join(x for x in (", ".join(it["phenotypes"][:3]), ", ".join(it["sample_terms"][:3]),
                                        ", ".join(it["assays"][:2])) if x)
            lines.append(f"| {it['accession']} | {it['file_set_type']} | {it['summary'][:110]} | {ctx} | "
                         f"{it['lab'] or ''} |")
        if t["total"] > len(t["items"]):
            lines.append(f"\n{t['total'] - len(t['items'])} more; raise --limit to list them.")
        lines.append("")
    if res.get("skipped"):
        lines += ["Types without a field for every criterion were skipped: "
                  + "; ".join(f"{k} (no {', '.join(v)})" for k, v in res["skipped"].items()), ""]
    return "\n".join(lines)


def fixture_fetch() -> "Tuple[Fetch, list]":
    """Offline Portal with a few prediction and measurement sets."""
    preds = [
        {"accession": "TSTDS1PRDCAD", "@type": ["PredictionSet"], "file_set_type": "disease associations",
         "summary": "ColocBoost CAD colocalisation", "associated_phenotypes": [{"term_name": "coronary artery disease"}],
         "samples": [{"sample_terms": [{"term_name": "heart left ventricle"}]}], "status": "released"},
        {"accession": "TSTDS2PRDOLD", "@type": ["PredictionSet"], "file_set_type": "disease associations",
         "summary": "old CAD", "associated_phenotypes": [{"term_name": "coronary artery disease"}],
         "superseded_by": ["/prediction-sets/TSTDS1PRDCAD/"], "status": "released"},
        {"accession": "TSTDS3PRDLDL", "@type": ["PredictionSet"], "file_set_type": "non-coding variant effects",
         "summary": "LDL", "associated_phenotypes": [{"term_name": "low density lipoprotein cholesterol measurement"}],
         "samples": [{"sample_terms": [{"term_name": "liver"}]}], "status": "released"},
    ]
    ms = [{"accession": "TSTDS4MSHRT", "@type": ["MeasurementSet"], "file_set_type": "experimental data",
           "summary": "snRNA-seq of heart", "samples": [{"sample_terms": [{"term_name": "heart left ventricle"}]}],
           "preferred_assay_titles": ["10x multiome"], "status": "released"}]
    objs = {"PredictionSet": preds, "MeasurementSet": ms}
    calls: list = []

    def value_of(it, field):
        level = [it]
        for part in field.split("."):
            nxt = []
            for x in level:
                y = x.get(part) if isinstance(x, dict) else None
                nxt.extend(y if isinstance(y, list) else [y] if y is not None else [])
            level = nxt
        return {str(x) for x in level if not isinstance(x, (dict, list))}

    def fetch(url: str):
        calls.append(url)
        q = urllib.parse.parse_qs(url.split("?", 1)[1])
        typ = q.get("type", [""])[0]
        items = objs.get(typ, [])
        for k, vals in q.items():
            if k in ("type", "limit", "format", "field", "sort"):
                continue
            items = [it for it in items if value_of(it, k) & set(vals)]
        if q.get("limit") == ["0"]:
            facets = []
            for crit in FIELDS.values():
                f = crit.get(typ)
                if f and f not in [x["field"] for x in facets]:
                    counts: "Dict[str, int]" = {}
                    for it in items:
                        for v in value_of(it, f):
                            counts[v] = counts.get(v, 0) + 1
                    facets.append({"field": f, "terms": [{"key": k, "doc_count": n} for k, n in counts.items()]})
            return 200, {"@graph": [], "total": len(items), "facets": facets}
        if not items:
            return 404, None
        return 200, {"@graph": items, "total": len(items)}
    return fetch, calls


def selftest() -> "List[Tuple[bool, str]]":
    fetch, calls = fixture_fetch()
    checks = []
    r = discover({"phenotype": "CAD coronary artery disease"}, fetch=fetch)
    pr = r["types"].get("PredictionSet", {})
    checks.append((pr.get("matched", {}).get("phenotype") == ["coronary artery disease"],
                   "phenotype words map to the Portal term 'coronary artery disease'"))
    checks.append(([i["accession"] for i in pr.get("items", [])] == ["TSTDS1PRDCAD"],
                   "superseded prediction sets are left out by default"))
    checks.append(("MeasurementSet" not in r["types"] or r["types"]["MeasurementSet"]["total"] == 0,
                   "types without the phenotype link are not claimed as matches"))
    r2 = discover({"tissue": "cardiac"}, types=["MeasurementSet", "PredictionSet"], fetch=fetch)
    checks.append((r2["types"]["MeasurementSet"]["matched"].get("tissue") == ["heart left ventricle"]
                   and r2["types"]["MeasurementSet"]["items"][0]["accession"] == "TSTDS4MSHRT",
                   "'cardiac' reaches 'heart left ventricle' through the synonym map"))
    r3 = discover({"phenotype": "migraine"}, types=["PredictionSet"], fetch=fetch)
    checks.append((r3["types"]["PredictionSet"]["total"] == 0 and bool(r3["types"]["PredictionSet"].get("no_match")),
                   "no matching phenotype: reported with the nearest values, not an error"))
    checks.append((match_terms("heart", {"heart": 5, "heart left ventricle": 3, "hearty": 1, "liver": 9})
                   == ["heart", "heart left ventricle"], "a tissue reaches its sub-terms by whole word"))
    checks.append((match_terms("element-gene links", {"spatial gene expression variability": 3}, strict=True) == [],
                   "categorical types never fall back to a different type sharing a word"))
    md = to_markdown(r)
    checks.append(("TSTDS1PRDCAD" in md and "processed lineage" in md, "markdown lists results and the next step"))
    return checks


if __name__ == "__main__":
    ok = True
    for good, msg in selftest():
        ok &= good
        print(("  ok    " if good else "  FAIL  ") + msg)
    raise SystemExit(0 if ok else 1)
