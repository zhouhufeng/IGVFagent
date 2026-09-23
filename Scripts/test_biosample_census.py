"""Offline checks for the biosample census and the search-URL diagnosis.

These pin the exact shapes that broke at runtime on 2026-09-21, none of
which need the network:

* an IGVF ``stats`` facet whose ``terms`` is a dict (``file_size``) -- the
  runtime-authored census iterated it as a list and crashed;
* ``hierarchical`` facets with a nested ``subfacet``; boolean / date facets
  that carry ``key_as_string``;
* legacy ENCODE / IGVF filter fields, which match nothing and answer 404;
* re-authoring an extension under its own name, which the duplication
  guard used to refuse as a duplicate of itself.

    python3 Scripts/test_biosample_census.py
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import biosample_census_skill as census  # noqa: E402
import data_illustration_interpretation as explain  # noqa: E402

FAILED: "list[str]" = []


def check(cond: bool, msg: str) -> None:
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        FAILED.append(msg)


IGVF_FILE_FACETS = {
    "total": 11,
    "facets": [
        {"field": "file_format", "type": "terms",
         "terms": [{"key": "fastq", "doc_count": 8}, {"key": "tsv", "doc_count": 3}]},
        {"field": "file_size", "type": "stats",
         "terms": {"count": 11, "min": 102.0, "max": 5e9, "avg": 1e8, "sum": 1.1e9}},
        {"field": "controlled_access", "type": "terms",
         "terms": [{"key": 0, "key_as_string": "false", "doc_count": 11}]},
        {"field": "file_set.assay_slims", "type": "hierarchical",
         "terms": [{"key": "gene expression", "doc_count": 5,
                    "subfacet": {"field": "file_set.assay_titles",
                                 "terms": [{"key": "MPRA", "doc_count": 5}]}}]},
        {"field": "status", "type": "terms",
         "terms": [{"key": "released", "doc_count": 10}, {"key": "in progress", "doc_count": 1}]},
        {"field": "broken", "terms": "not a list"},
        "not a dict",
    ],
}


def main() -> int:
    print("\nfacet parsing survives every portal shape")
    ft = census.facet_table(IGVF_FILE_FACETS)
    check(ft["file_format"] == [("fastq", 8), ("tsv", 3)], "plain terms list, sorted desc")
    check("file_size" not in ft, "stats dict is not iterated as terms")
    check(ft["controlled_access"] == [("false", 11)], "key_as_string wins over numeric key")
    check(ft["file_set.assay_slims"] == [("gene expression", 5)], "hierarchical parent kept")
    check(ft["file_set.assay_titles"] == [("MPRA", 5)], "hierarchical subfacet flattened")
    check("broken" not in ft, "non-list terms skipped")
    st = census.facet_stats(IGVF_FILE_FACETS)
    check(st["file_size"]["sum"] == 1.1e9 and st["file_size"]["count"] == 11, "stats facet summarised")
    wanted = census.facet_table(IGVF_FILE_FACETS, ["status", "file_set.assay_slims"])
    check(set(wanted) == {"status", "file_set.assay_slims", "file_set.assay_titles"},
          "wanted filter keeps the subfacet of a wanted parent")

    print("\nzero-hit 404 is a count of zero")
    check(census.is_zero_hit(404, {"total": 0, "notification": "No results found"}), "404 + total 0")
    check(not census.is_zero_hit(404, {"http_error": 404}), "bare 404 is not a zero hit")
    check(not census.is_zero_hit(200, {"total": 0}), "200 with total 0 is not a 404")

    print("\nitem flattening follows nested lists")
    row = {"accession": "ENCSR1", "target": {"label": "CTCF"}, "@type": ["Experiment", "Dataset"],
           "replicates": [{"library": {"biosample": {"treatments": [
               {"treatment_term_name": "DMSO"}, {"treatment_term_name": "IFN"}]}}},
               {"library": {"biosample": {"treatments": [{"treatment_term_name": "DMSO"}]}}}],
           "biosample_ontology": {"term_id": "EFO:0002784"}}
    flat = census.flatten_item(row, ["accession", "target.label", "@type",
                                     "replicates.library.biosample.treatments.treatment_term_name",
                                     "biosample_ontology.term_id", "missing.field"])
    check(flat["target.label"] == "CTCF", "embedded object -> label")
    check(flat["@type"] == "Experiment", "@type -> first entry")
    check(flat["replicates.library.biosample.treatments.treatment_term_name"] == "DMSO; IFN",
          "nested lists flattened and de-duplicated")
    check(flat["biosample_ontology.term_id"] == "EFO:0002784", "dotted path into dict")
    check(flat["missing.field"] == "", "missing path is empty, not an error")

    print("\ntables")
    yrs = census.by_year([{"date_released": "2011-05-01"}, {"date_released": "2013-01-01"},
                          {"date_released": ""}, {"date_released": "n/a"}])
    check(yrs == [("2011", 1), ("2012", 0), ("2013", 1)], "gap years appear as zero")
    rows = [{"portal": "encode", "type": "Experiment", "facet": "assay_title", "term": t, "count": c}
            for t, c in (("a", 5), ("b", 4), ("c", 3), ("d", 2))]
    check(census.top_terms(rows, "encode", "Experiment", "assay_title", 2) ==
          [("a", 5), ("b", 4), ("Other (2 more)", 5)], "top-N folds the remainder into Other")
    check(census.count_items([{"lab.title": "X; Y"}, {"lab.title": "X"}], "lab.title", 10) ==
          [("X", 2), ("Y", 1)], "semicolon-joined values counted separately")

    print("\nfigures are well-formed SVG and never divide by zero")
    import xml.dom.minidom
    with tempfile.TemporaryDirectory() as td:
        p1 = census.write_hbar_svg([("a", 0), ("b", 0)], "t", Path(td) / "h.svg", "#2a78d6")
        p2 = census.write_vbar_svg([], "t", Path(td) / "v.svg", "#2a78d6", "year")
        p3 = census.write_panels_svg([("E", [("x", 3)], "#2a78d6"), ("I", [], "#eb6834")], "t",
                                     Path(td) / "p.svg")
        for p in (p1, p2, p3):
            try:
                xml.dom.minidom.parse(str(p))
                check(True, f"{p.name} parses")
            except Exception as exc:
                check(False, f"{p.name} parses ({exc})")
        svg = p1.read_text()
        check("&lt;" not in svg and 'fill="#2a78d6"' in svg, "series colour applied")

    print("\nsearch-URL normalisation")
    pairs, notes = explain.normalize_search_query(
        "encode", [("type", "Experiment"), ("biosample_term_name", "GM12878"),
                   ("biosample_type!", "tissue"), ("assay_title", "TF ChIP-seq")])
    check(dict(pairs)["biosample_ontology.term_name"] == "GM12878", "legacy field rewritten")
    check(dict(pairs)["biosample_ontology.classification!"] == "tissue", "negated legacy field keeps the !")
    check(dict(pairs)["assay_title"] == "TF ChIP-seq" and len(notes) == 2, "valid field untouched, notes emitted")
    pairs, notes = explain.normalize_search_query("igvf", [("samples.summary", "GM12878")])
    check(dict(pairs) == {"samples.sample_terms.term_name": "GM12878"}, "IGVF legacy field rewritten")
    src, url, label = explain.build_json_url(
        "https://www.encodeproject.org/report/?type=Experiment&biosample_term_name=GM12878&format=json")
    check(src == "encode" and "/search/" in url and url.count("format=json") == 1
          and "biosample_ontology.term_name=GM12878" in url, "report view -> search JSON, one format param")
    check(explain.is_search_url("https://www.encodeproject.org/search/?type=Experiment"), "search URL detected")
    check(not explain.is_search_url("ENCSR000EMT"), "accession is not a search URL")
    check(not explain.is_search_url("https://www.encodeproject.org/experiments/ENCSR000EMT/"), "object URL is not a search URL")

    print("\nempty-search diagnosis (facet probe stubbed)")
    probe = {"facets": [
        {"field": "biosample_ontology.term_name", "terms": [{"key": "GM12878", "doc_count": 1},
                                                            {"key": "K562", "doc_count": 1}]},
        {"field": "assay_title", "terms": [{"key": "TF ChIP-seq", "doc_count": 1}]},
        {"field": "status", "terms": [{"key": "released", "doc_count": 1}]},
    ]}
    real_fetch = explain.fetch_json
    explain.fetch_json = lambda source, url: (200, probe)  # type: ignore
    try:
        body = {"total": 0, "notification": "No results found", "filters": [
            {"field": "type", "term": "Experiment"},
            {"field": "biosample_ontology.term_name", "term": "gm12878"},
            {"field": "biosample_term_name", "term": "GM12878"},
            {"field": "assay_titel", "term": "TF ChIP-seq"},
            {"field": "assay_title", "term": "ChIP-seq TF"},
            {"field": "status", "term": "released"}]}
        lines = "\n".join(explain.diagnose_empty_search("encode", body))
    finally:
        explain.fetch_json = real_fetch  # type: ignore
    check("spells it `GM12878`" in lines, "case mismatch named")
    check("retired field" in lines and "biosample_ontology.term_name" in lines, "legacy field named with replacement")
    check("`assay_titel` is not a facet" in lines and "assay_title" in lines, "typo in field name -> closest facet")
    check("No Experiment has `assay_title=ChIP-seq TF`" in lines and "TF ChIP-seq" in lines, "bad value -> closest value")
    check("biosample_portal_census" in lines, "points at the census tool")
    no_filters = explain.diagnose_empty_search("encode", {"total": 0})
    check(any("no filters" in ln for ln in no_filters), "handles a body without filters")

    print("\ndownload redirects drop credentials only when they leave the host")
    import urllib.request
    h = explain._AuthStrippingRedirectHandler()
    req = urllib.request.Request("https://api.data.igvf.org/tabular-files/X/@@download/X.tsv.gz",
                                 headers={"Authorization": "Basic abc", "Cookie": "c=1", "User-Agent": "t"})
    cross = h.redirect_request(req, None, 307, "Temporary Redirect", {},
                               "https://igvf-public.s3.amazonaws.com/x.tsv.gz?X-Amz-Signature=1")
    same = h.redirect_request(req, None, 307, "Temporary Redirect", {}, "https://api.data.igvf.org/other/")
    check(cross is not None and not cross.has_header("Authorization") and not cross.has_header("Cookie")
          and cross.has_header("User-agent"), "S3 redirect: Authorization and Cookie removed, other headers kept")
    check(same is not None and same.has_header("Authorization"), "same-host redirect keeps Authorization")

    print("\nextension re-authoring is an update, not a duplicate")
    import ext_author_skill as ext  # noqa: E402
    with tempfile.TemporaryDirectory() as td:
        env = dict(os.environ)
        os.environ["IGVF_USER_EXT_DIR"] = td
        os.environ["IGVF_ALLOW_AGENT_AUTHORING"] = "1"
        os.environ["HOME"] = td  # keep ~/.igvfagent out of the picture
        try:
            src1 = "def main(argv=None):\n    print('v1')\n    return 0\n"
            src2 = src1.replace("v1", "v2")
            rc1 = ext.main(["write-skill", "--name", "zebra_tally", "--description",
                            "Tally zebra stripes in a text file", "--source", src1,
                            "--tool-parameters", '{"type":"object","properties":{"label":{"type":"string"}}}'])
            check(rc1 == 0, "first authoring accepted")
            check(ext.is_update("zebra_tally"), "name now counts as an existing extension")
            rc2 = ext.main(["write-skill", "--name", "zebra_tally", "--description",
                            "Tally zebra stripes in a text file", "--source", src2,
                            "--tool-parameters", '{"type":"object","properties":{"top":{"type":"integer"}}}'])
            check(rc2 == 0, "re-authoring under the same name accepted")
            skill = Path(td) / "skills" / "zebra_tally.py"
            check("v2" in skill.read_text(), "new source on disk")
            check((Path(td) / "skills" / "zebra_tally.py.prev").read_text().count("v1") == 1,
                  "previous source kept as .prev")
            manifest = json.loads((Path(td) / "tools" / "zebra_tally.json").read_text())
            check("top" in manifest["parameters"]["properties"], "manifest updated")
            check(not [s for s in ext._userext.discover_skills() if s.endswith(".prev")],
                  ".prev files are invisible to the loader")
            hits = ext.find_similar_core_tools("zebra_tally", "Tally zebra stripes in a text file")
            check(not hits, "guard does not score the proposal against its own earlier version")
            rc3 = ext.main(["write-skill", "--name", "biosample_portal_census", "--description",
                            "x", "--source", src1])
            check(rc3 in (2, 3), "a built-in tool name is still refused")
        finally:
            os.environ.clear()
            os.environ.update(env)

    if FAILED:
        print(f"\n{len(FAILED)} check(s) failed")
        return 1
    print("\nall checks pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
