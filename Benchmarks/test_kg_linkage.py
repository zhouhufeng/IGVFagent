#!/usr/bin/env python3
"""Gene-centric regulatory evidence: no false absences, no borrowed genes.

Written from an external evaluation (Docs/BugZilla/
IGVFagent_MultiFlow_Public_Data_Test_Report_2026-09-09.md) that found the
gene workflow reporting "no kidney-specific regulatory records" for WT1 and
HNF4A, and attributing other genes' element-to-gene edges to GATA3 and SOX9.

Both were reproduced, and the cause was in three places at once:

  fetch_linkage_for_region issued ONE capped request and never paged, so
  every result was truncated and no caller could tell;

  it never compared a row's own `gene` field to the queried gene, and the
  endpoint returns edges for elements OVERLAPPING the region -- of the first
  25 rows over the WT1 locus, 23 distinct target genes appear and none is
  WT1;

  listify() ended in `return [data]`, so a 404 body came back as one data
  row -- and the primary endpoint it called, /api/regulatory-regions/genes,
  404s for every region.

Measured after the fix, exhaustively and gene-filtered: WT1 4,458 target-gene
rows of which 405 are kidney/renal; HNF4A 2,398 and 5; GATA3 7,934 and 600;
SOX9 2,934 and 185. The external audit's own limit=500 lower bounds (390 /
346 / 235 / 181 target-gene rows) were reproduced exactly before the fix.

These cases are offline: the API shape is stubbed, so they pin the logic
rather than the database's current contents.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import kg_traversal_skill as kg  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:58} {detail}")
    if not ok:
        FAILURES.append(name)


# ── an API error must never become a data row ──────────────────────────────

check("404 body is not a data row",
      kg.listify({"message": "Not found", "code": "NOT_FOUND"}) == [])
check("error-shaped body detected", kg.is_error_body({"error": "boom"}))
check("{'detail': ...} detected", kg.is_error_body({"detail": "nope"}))
check("a real single record still survives",
      len(kg.listify({"gene": "genes/ENSG1", "class": "prediction"})) == 1)
check("a wrapped list is unwrapped",
      len(kg.listify({"results": [{"a": 1}, {"b": 2}]})) == 2)
check("non-dict junk is dropped", kg.listify([1, "x", {"a": 1}]) == [{"a": 1}])


# ── target gene extraction ────────────────────────────────────────────────

# The real field is a prefixed string. Splitting on "." matched nothing,
# which is indistinguishable from "this gene has no records".
check("prefixed string id", kg.target_gene_id({"gene": "genes/ENSG00000184937"})
      == "ENSG00000184937")
check("versioned id keeps the base",
      kg.target_gene_id({"gene": "genes/ENSG00000184937.12"}) == "ENSG00000184937")
check("embedded object id",
      kg.target_gene_id({"gene": {"_id": "ENSG00000101076"}}) == "ENSG00000101076")
check("missing gene field -> empty, not a crash",
      kg.target_gene_id({}) == "")


# ── pagination: page, never skip ──────────────────────────────────────────

class Stub:
    """Stands in for the Catalog. Records how it was called."""

    def __init__(self, total, page_cap=500, ignore_skip=True, fail_on=None):
        self.total, self.cap = total, page_cap
        self.ignore_skip, self.fail_on = ignore_skip, fail_on
        self.calls = []

    def __call__(self, path, **params):
        self.calls.append(params)
        if self.fail_on is not None and len(self.calls) > self.fail_on:
            return 503, {"message": "boom", "code": "ERR"}
        lim = min(int(params.get("limit", self.cap)), self.cap)
        pg = int(params.get("page", 0))
        # The real API ignores skip/offset; a pager relying on them loops.
        if self.ignore_skip and ("skip" in params or "offset" in params):
            pg = 0
        start = pg * lim
        return 200, [{"gene": f"genes/ENSG{i % 3}", "i": i}
                     for i in range(start, min(start + lim, self.total))]


def with_stub(stub, **kw):
    orig = kg.catalog_get
    kg.catalog_get = stub
    try:
        return kg.catalog_paged("/api/genomic-elements/genes", **kw)
    finally:
        kg.catalog_get = orig


rows, meta = with_stub(Stub(1250), page_limit=500)
check("pages until exhausted", len(rows) == 1250, f"{len(rows)} rows")
check("reports exhausted, not truncated", meta["truncated"] is False,
      meta["stopped_because"])
check("page count is right", meta["pages"] == 3, f"pages={meta['pages']}")

rows, meta = with_stub(Stub(100_000), page_limit=500, max_pages=3)
check("max_pages stops and MARKS TRUNCATED", meta["truncated"] is True,
      meta["stopped_because"])
check("truncated run still returns what it got", len(rows) == 1500,
      f"{len(rows)} rows")

# An endpoint repeating page 1 must not spin forever.
rows, meta = with_stub(Stub(50, page_cap=10, ignore_skip=True), page_limit=10)
check("repeating endpoint terminates", len(rows) <= 50)

rows, meta = with_stub(Stub(5000, fail_on=2), page_limit=500)
check("an HTTP error mid-page marks truncation",
      meta["truncated"] and "503" in meta["stopped_because"],
      meta["stopped_because"])
st = Stub(1200)
with_stub(st, page_limit=500)
check("pager uses page=, not skip=",
      all("page" in c and "skip" not in c and "offset" not in c
          for c in st.calls), f"{st.calls[:2]}")


# ── exact target-gene filtering, the GATA3 / SOX9 case ────────────────────

def linkage_with(total, gene_id=None, **kw):
    orig = kg.catalog_get
    kg.catalog_get = Stub(total)
    try:
        return kg.fetch_linkage_for_region("chr1:1-2", gene_id=gene_id, **kw)
    finally:
        kg.catalog_get = orig


# Stub rows cycle ENSG0/1/2, so a third of them target the queried gene.
L = linkage_with(900, gene_id="ENSG1")
m = L["meta"]
check("only the queried gene's rows are kept",
      all(kg.target_gene_id(r) == "ENSG1" for r in L["region_predictions"]),
      f"{m['rows_for_target_gene']} kept")
check("other genes' rows are counted as dropped",
      m["rows_dropped_other_genes"] == 900 - m["rows_for_target_gene"],
      f"dropped={m['rows_dropped_other_genes']}")
check("distinct target genes reported", m["distinct_target_genes"] == 3)
check("target_gene_id recorded in meta", m["target_gene_id"] == "ENSG1")

# Unfiltered is allowed for a REGION question, but must be self-declaring.
L2 = linkage_with(30)
check("unfiltered keeps every row", len(L2["region_predictions"]) == 30)
check("unfiltered declares no target gene",
      L2["meta"]["target_gene_id"] is None)
check("unfiltered has no rows_for_target_gene count",
      "rows_for_target_gene" not in L2["meta"])

# The dead endpoint must not be called at all.
calls = []


def spy(path, **params):
    calls.append(path)
    return 200, []


orig = kg.catalog_get
kg.catalog_get = spy
try:
    kg.fetch_linkage_for_region("chr1:1-2", exhaustive=True)
finally:
    kg.catalog_get = orig
check("the 404 endpoint is never requested",
      all("regulatory-regions/genes" not in c for c in calls), str(set(calls)))


# ── absence claims and prediction/observation split ───────────────────────

note = kg.evidence_note({"rows_for_target_gene": 0, "truncated": True,
                          "stopped_because": "hit max_pages=3",
                          "observed_count": 0, "prediction_count": 0})
check("a truncated empty result forbids an absence claim",
      "TRUNCATED" in note and "absence cannot be concluded" in note, note[:70])
note2 = kg.evidence_note({"rows_for_target_gene": 0, "truncated": False,
                           "pages": 4, "observed_count": 0,
                           "prediction_count": 0})
check("an exhausted empty result may state absence",
      "TRUNCATED" not in note2 and "exhausted" in note2, note2[:70])


def classed(rows):
    orig = kg.catalog_get
    kg.catalog_get = lambda p, **kw: (200, rows)
    try:
        return kg.fetch_linkage_for_region("chr1:1-2", exhaustive=False,
                                            limit=len(rows) + 1)["meta"]
    finally:
        kg.catalog_get = orig


m = classed([{"gene": "genes/G", "class": "prediction", "score": 0.9,
              "p_value_adj": None, "significant": None},
             {"gene": "genes/G", "class": "observed data", "significant": False,
              "p_value_adj": 0.93}])
check("predictions and observations are counted apart",
      (m["observed_count"], m["prediction_count"]) == (1, 1),
      f"obs={m['observed_count']} pred={m['prediction_count']}")


# ── record-level answer verification ──────────────────────────────────────
#
# The report's remaining P0: "Each reported gene, element, biosample, score,
# effect, adjusted P value, significance status, and source URL should map to
# the same source record. The verifier should reject cross-record field
# assembly." That failure is invisible to spot-checking, because every value
# in the sentence is genuine -- only the combination is not.

SRC = [
    {"gene": "genes/ENSG00000184937", "genomic_element": "genomic_elements/E1",
     "biological_context": "left kidney from ENCDO633IHH", "method": "ENCODE-rE2G",
     "class": "prediction", "source_url": "https://x/ENCFF899XZX",
     "score": 0.7911, "p_value_adj": None, "significant": None},
    {"gene": "genes/ENSG00000184937", "genomic_element": "genomic_elements/E2",
     "biological_context": "liver", "method": "ENCODE-rE2G",
     "class": "prediction", "source_url": "https://x/ENCFF000AAA",
     "score": 0.2200, "p_value_adj": None, "significant": None},
    {"gene": "genes/ENSG00000184937", "genomic_element": "genomic_elements/E3",
     "biological_context": "CD8-positive T cell", "method": "CRISPR screen",
     "class": "observed data", "source_url": "https://x/IGVFFI1",
     "score": None, "p_value_adj": 0.938, "significant": False},
]
RECS = kg.evidence_records(SRC)

check("every record gets an id", len({r["record_id"] for r in RECS}) == 3)
check("record ids are stable across calls",
      [r["record_id"] for r in kg.evidence_records(SRC)]
      == [r["record_id"] for r in RECS])

v = kg.verify_claim({"gene": "ENSG00000184937",
                     "biological_context": "left kidney from ENCDO633IHH",
                     "score": 0.7911, "source_url": "https://x/ENCFF899XZX"}, RECS)
check("a faithful row verifies", v["ok"] and v["record_id"] == RECS[0]["record_id"])
check("rounding is tolerated",
      kg.verify_claim({"gene": "ENSG00000184937", "score": 0.791}, RECS)["ok"])
check("the genes/ prefix is tolerated",
      kg.verify_claim({"gene": "genes/ENSG00000184937",
                       "score": 0.7911}, RECS)["ok"])

# The headline case: kidney from record 1, score from record 2.
v = kg.verify_claim({"gene": "ENSG00000184937",
                     "biological_context": "left kidney from ENCDO633IHH",
                     "score": 0.2200}, RECS)
check("CROSS-RECORD assembly is rejected", not v["ok"], v["reason"])
check("and is named as such, not 'not found'", v["reason"] == "CROSS_RECORD",
      v["reason"])
check("it reports that each value did occur somewhere",
      set(v["fields_found_elsewhere"]) ==
      {"gene", "biological_context", "score"}, str(v["fields_found_elsewhere"]))

# Source URL borrowed from another record is the same failure.
v = kg.verify_claim({"biological_context": "left kidney from ENCDO633IHH",
                     "source_url": "https://x/ENCFF000AAA"}, RECS)
check("a borrowed source URL is rejected", v["reason"] == "CROSS_RECORD",
      v["reason"])

v = kg.verify_claim({"gene": "ENSG00000184937",
                     "biological_context": "renal cortex", "score": 0.99}, RECS)
check("invented values are NOT_IN_SOURCE", v["reason"] == "NOT_IN_SOURCE",
      v["reason"])
check("and the unmatched fields are named",
      set(v.get("unmatched_fields", [])) == {"biological_context", "score"},
      str(v.get("unmatched_fields")))

# A null adjusted p-value quoted as a number: the significance confusion.
check("a null p_value_adj quoted as a number is rejected",
      not kg.verify_claim({"gene": "ENSG00000184937", "score": 0.7911,
                           "p_value_adj": 0.01}, RECS)["ok"])
check("null matches null",
      kg.verify_claim({"gene": "ENSG00000184937", "score": 0.7911,
                       "p_value_adj": None, "significant": None},
                      RECS)["ok"])
check("significant=False is not treated as absent",
      kg.verify_claim({"genomic_element": "E3", "significant": False},
                      RECS)["ok"])
check("significant=True on a False record is rejected",
      not kg.verify_claim({"genomic_element": "E3", "significant": True},
                          RECS)["ok"])
check("a claim with no citable field is rejected",
      kg.verify_claim({"notes": "looks good"}, RECS)["reason"]
      == "NO_CITABLE_FIELDS")

t = kg.verify_table([
    {"gene": "ENSG00000184937", "score": 0.7911},
    {"gene": "ENSG00000184937", "biological_context": "left kidney from ENCDO633IHH",
     "score": 0.2200},
    {"gene": "ENSG00000184937", "biological_context": "renal cortex"},
], RECS)
check("verify_table counts verified and rejected",
      (t["n_claims"], t["n_verified"], t["n_rejected"]) == (3, 1, 2), str(t)[:80])
check("verify_table separates the two failure kinds",
      (t["cross_record"], t["not_in_source"]) == (1, 1),
      f"cross={t['cross_record']} missing={t['not_in_source']}")
check("verify_table is not ok when any row fails", t["ok"] is False)

print(f"\n{28 + 18} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
