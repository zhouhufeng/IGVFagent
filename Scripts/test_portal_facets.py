"""Offline checks for `portal facets` edge cases seen in real agent runs.

    python3 Scripts/test_portal_facets.py

1. The IGVF Portal answers a zero-result search with HTTP 404 and an empty
   Search body. That must be read as 0 results (exit 0), with suggestions
   of real values, not as a failed tool call.
2. A search over several types returns only their shared facets, so a
   type-specific field (file_set_type) must be counted from the items.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault("IGVF_PROJECT_ROOT", tempfile.mkdtemp(prefix="pfacets_"))
import portal_query_skill as q  # noqa: E402

FAIL = []


def check(name, cond):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
    if not cond:
        FAIL.append(name)


EMPTY = {"@type": ["Search"], "@graph": [], "facets": [], "total": 0}
ITEMS = {"@type": ["Search"], "total": 3, "@graph": [
    {"@type": ["AnalysisSet"], "file_set_type": "principal analysis",
     "samples": [{"sample_terms": [{"term_name": "heart left ventricle"}]}]},
    {"@type": ["AnalysisSet"], "file_set_type": "intermediate analysis",
     "samples": [{"sample_terms": [{"term_name": "heart"}, {"term_name": "aorta"}]},
                 {"sample_terms": [{"term_name": "heart"}]}]},
    {"@type": ["PredictionSet"], "file_set_type": "principal analysis"},
]}
SHARED = {"@type": ["Search"], "@graph": [], "total": 3,
          "facets": [{"field": "type", "title": "Object Type", "terms": []}]}
TERMS = {"@type": ["Search"], "@graph": [], "total": 9,
         "facets": [{"field": "samples.sample_terms.term_name", "terms": [
             {"key": "heart", "doc_count": 5}, {"key": "K562", "doc_count": 40},
             {"key": "anterior left cardiac atrium", "doc_count": 2}]}]}


def fake(routes):
    def _request(url, **kw):
        for needle, body, status in routes:
            if needle(url):
                return status, json.dumps(body).encode(), "application/json"
        raise AssertionError(f"unexpected URL {url}")
    return _request


def run(argv_ns):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = q.cmd_facets(argv_ns)
    return rc, buf.getvalue()


def ns(**kw):
    base = dict(type=None, query=None, field_filters=None, field=None, label="t")
    base.update(kw)
    return argparse.Namespace(**base)


def main() -> int:
    q.setup_logging = lambda *a, **k: None
    q.REPORT_DIR = Path(os.environ["IGVF_PROJECT_ROOT"]) / "Docs" / "PortalQuery"

    data = q._ensure_json(404, json.dumps(EMPTY).encode(), "application/json",
                          context="t")
    check("a 404 carrying an empty Search is read as zero results",
          data.get("total") == 0)
    try:
        q._ensure_json(404, b'{"@type": ["HTTPNotFound"]}', "application/json",
                       context="t")
        check("a real 404 still fails", False)
    except SystemExit:
        check("a real 404 still fails", True)

    q._request = fake([
        (lambda u: "cardiac" in u, EMPTY, 404),
        (lambda u: "limit=0" in u, TERMS, 200),
    ])
    rc, out = run(ns(type="MeasurementSet",
                     field_filters="samples.sample_terms.term_name=cardiac muscle cell"))
    check("zero-result filtered facets exit 0", rc == 0)
    check("...say so plainly", "No items match these filters" in out)
    check("...and suggest real cardiac/heart terms first",
          out.find("heart=5") != -1 and out.find("K562") == -1
          and out.find("anterior left cardiac atrium=2") != -1)

    q._request = fake([
        (lambda u: "field=file_set_type" in u, ITEMS, 200),
        (lambda u: "limit=0" in u, SHARED, 200),
    ])
    rc, out = run(ns(type="AnalysisSet,PredictionSet", field="file_set_type"))
    check("a facet absent on a mixed-type search is counted from items (exit 0)",
          rc == 0 and "principal analysis=2" in out
          and "intermediate analysis=1" in out)

    q._request = fake([(lambda u: True, ITEMS, 200)])
    counts, seen, n = q._count_field_from_items(
        ["AnalysisSet"], None, [], "samples.sample_terms.term_name")
    check("nested lists are flattened, each item counted once per value",
          counts == {"heart left ventricle": 1, "heart": 1, "aorta": 1}
          and seen == 2 and n == 3)

    print("all checks pass" if not FAIL else f"FAILED: {len(FAIL)}")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
