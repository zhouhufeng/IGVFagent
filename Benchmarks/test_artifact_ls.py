#!/usr/bin/env python3
"""A listing of a path that does not exist must not be a dead end.

On a real run the model asked for `Docs/BaseEditingScreen/18loci_uptake` when
the analysis had written `Docs/BaseEditingScreen/<timestamp>_18loci_uptake`.
`artifact ls` exited 2, the orchestrator recorded "1 tool call did not
succeed", and the whole answer was reported as incomplete -- for a directory
that was one listing of the parent away.

An agent cannot learn the right name from a refusal. It can from a listing.
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Scripts"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:62} {detail}")
    if not ok:
        FAILURES.append(name)


_tmp = tempfile.TemporaryDirectory()
import os  # noqa: E402
# .resolve() matters: on macOS a tempdir under /var resolves to
# /private/var, and _pathguard resolves its root, so an unresolved root
# compares as "outside the workspace".
root = Path(_tmp.name).resolve()
os.environ["IGVF_PROJECT_ROOT"] = str(root)
(root / "Docs" / "BaseEditingScreen" / "20260911_033012_18loci_uptake").mkdir(parents=True)
(root / "Docs" / "BaseEditingScreen" / "ready_check").mkdir()
(root / "Docs" / "BaseEditingScreen" / "20260911_033012_18loci_uptake"
 / "summary.json").write_text("{}")

import _pathguard  # noqa: E402
_pathguard.project_root = lambda: root
import artifact_read_skill as ar  # noqa: E402
ar._root = lambda: root

# ── the exact miss from the live run ─────────────────────────────────────
out = ar.list_artifacts("Docs/BaseEditingScreen/18loci_uptake")
check("a missing path returns a result, not an exception", isinstance(out, dict))
check("it says the path does not exist", out.get("exists") is False)
check("it names the nearest existing directory",
      out.get("nearest_existing") == "Docs/BaseEditingScreen")
check("it lists what is actually there",
      "20260911_033012_18loci_uptake" in out.get("available", []))
check("it suggests the near-miss, which differs only by a timestamp prefix",
      "20260911_033012_18loci_uptake" in out.get("did_you_mean", []),
      str(out.get("did_you_mean")))
check("the note tells the caller what to do next",
      "does not exist" in out.get("note", "")
      and "available" in out.get("note", ""))
check("entries is empty rather than absent, so a caller can iterate safely",
      out.get("entries") == [])

# ── a path that exists still works exactly as before ─────────────────────
ok = ar.list_artifacts("Docs/BaseEditingScreen/20260911_033012_18loci_uptake")
check("an existing directory reports exists=True", ok.get("exists") is True)
check("and lists its files",
      any(e["name"] == "summary.json" for e in ok["entries"]))
check("the parent lists both run directories",
      len(ar.list_artifacts("Docs/BaseEditingScreen")["entries"]) == 2)

# ── the suggestion must not fire on an unrelated name ────────────────────
miss = ar.list_artifacts("Docs/BaseEditingScreen/zzz_unrelated")
check("an unrelated name gets no false suggestion",
      miss.get("did_you_mean") == [], str(miss.get("did_you_mean")))
check("but it still lists what is available",
      len(miss.get("available", [])) == 2)

# ── walking up more than one level ───────────────────────────────────────
deep = ar.list_artifacts("Docs/BaseEditingScreen/nope/deeper/still")
check("it walks up to the first directory that exists",
      deep.get("nearest_existing") == "Docs/BaseEditingScreen",
      str(deep.get("nearest_existing")))

# ── containment is still enforced ────────────────────────────────────────
try:
    ar.list_artifacts("/etc")
    outside = False
except PermissionError:
    outside = True
check("a path outside the workspace is still refused", outside)


# ─── ranking must be over the whole file, never a sample ──────────────────
# From the 10 Sep retest: the hosted agent reported a "top 3 by score" for
# GATA3, SOX9 and WT1 that it had obtained from bounded grep hits and file
# excerpts. The rows were genuine and correctly attributed; the RANKING was
# false. For GATA3 it named 0.9909405 as the top score while 0.9999999981 sat
# in the same file. There was no tool that sorted an artefact, so the model
# used what existed and described the result as a ranking.
import csv as _csv  # noqa: E402

_f = root / "Docs" / "links.csv"
with open(_f, "w", newline="") as fh:
    w = _csv.writer(fh)
    w.writerow(["target", "biosample", "score", "source"])
    w.writerow(["GATA3", "kidney glomerular epithelial cell", "0.9999999981", "TOP1"])
    w.writerow(["GATA3", "renal cortical epithelial cell", "0.9999999974", "TOP2"])
    w.writerow(["GATA3", "kidney", "0.9999988818", "TOP3"])
    w.writerow(["GATA3", "adrenal gland", "0.9999999999", "ADRENAL"])
    w.writerow(["GATA3", "kidney", "0.9909405", "SAMPLED1"])
    w.writerow(["GATA3", "kidney", "not_a_number", "BADROW"])
    for i in range(500):
        w.writerow(["GATA3", "liver", "0.5", f"L{i}"])

r = ar.rank_artifact("Docs/links.csv", column="score", n=3,
                      where="kidney,renal", exclude="adrenal",
                      where_column="biosample")
check("ranking finds the true top record, not a sampled one",
      r["rows"][0]["source"] == "TOP1", str(r["values"][:1]))
check("and the true 2nd and 3rd",
      [x["source"] for x in r["rows"]] == ["TOP1", "TOP2", "TOP3"],
      str([x["source"] for x in r["rows"]]))
# 'kidney' alone misses 'renal cortical epithelial cell'; 'renal' alone also
# matches 'adrenal gland'. The audit needed both terms and an exclusion.
one_term = ar.rank_artifact("Docs/links.csv", column="score", n=3,
                             where="kidney", where_column="biosample")
check("a single 'kidney' term misses the renal-cortex record",
      "TOP2" not in [x["source"] for x in one_term["rows"]])
naive = ar.rank_artifact("Docs/links.csv", column="score", n=1,
                          where="renal", where_column="biosample")
check("a bare 'renal' term would wrongly admit adrenal",
      naive["rows"][0]["source"] == "ADRENAL")
check("which is why exclude exists and beats it",
      r["rows"][0]["source"] != "ADRENAL")
# Unparseable scores must be counted, never silently ranked as zero.
check("rows whose score does not parse are excluded and counted",
      r["n_excluded_unparseable"] == 1, str(r["n_excluded_unparseable"]))
check("the ranking reports how much it scanned",
      r["n_scanned"] == 506, str(r["n_scanned"]))
check("and how many it actually ranked", r["n_ranked"] == 4, str(r["n_ranked"]))
check("it asserts the ranking was complete", r["ranking_is_complete"] is True)
# A missing column must be an error, not an empty ranking that reads as
# "there is nothing here".
bad = ar.rank_artifact("Docs/links.csv", column="nosuchcol", n=3)
check("a missing ranking column is an error, not an empty result",
      "error" in bad and "nosuchcol" in bad["error"])
check("and the error lists the columns that do exist",
      "score" in str(bad.get("columns", [])))
asc = ar.rank_artifact("Docs/links.csv", column="score", n=1, ascending=True)
check("ascending ranks lowest-first", asc["values"][0] == 0.5)


# ─── a view limit is not a retrieval limit ────────────────────────────────
# The hosted answer said Catalog retrieval was "TRUNCATED" when the traversal
# had been exhaustive and only its own read_artifact view was cut. The word
# came straight from read_artifact's header, which said TRUNCATED without
# saying what had been truncated.
_big = root / "Docs" / "big.txt"
_big.write_text("x" * 500_000)
out = ar.read_artifact("Docs/big.txt", max_bytes=1000)
check("a cut response is flagged", out["truncated"] is True)
check("and carries an unambiguous alias", out["view_truncated"] is True)
check("whose meaning says it is about the VIEW, not the retrieval",
      "says nothing about whether the data in the file was completely "
      "retrieved" in out["truncation_meaning"])
full = ar.read_artifact("Docs/links.csv")
check("an untruncated read says so", full["truncated"] is False)
check("and its meaning field is unambiguous too",
      full["truncation_meaning"] == "full file returned")
# The rendered header is what the model actually copies into its answer.
src = Path(ar.__file__).read_text()
check("the printed header says VIEW TRUNCATED, not bare TRUNCATED",
      "VIEW TRUNCATED" in src)
check("and says it is a display limit",
      "display limit of" in src and "NOT a limit on how much data was" in src)
check("and points at the tool that reads the whole file",
      "artifact top" in src)

# ── a truthful record of an absence is not a validation failure ───────────
# The GSE213151 audit showed raw_rna_matrix_listed=false and
# atac_peak_matrix_listed=false as FAILED checks. Both values were correct --
# GEO does not supply those files -- so the manifest was valid and the report
# told the user their file was broken.
_man = root / "Docs" / "manifest.csv"
_man.write_text(
    "sample_id,cell_line,day,rna_gsm,atac_gsm,raw_rna,peak_mtx\n"
    "AN1_d7,AN1,d7,G1,A1,false,false\n"
    "AN1_d26,AN1,d26,G2,A2,false,false\n"
    "BJFF_d26,BJFF,d26,G3,A3,false,false\n")
au = ar.audit_manifest("Docs/manifest.csv", unique="sample_id,rna_gsm,atac_gsm",
                        pair="rna_gsm:atac_gsm", group="cell_line,day",
                        absent_ok="raw_rna,peak_mtx")
_st = {c["check"]: c["status"] for c in au["checks"]}
check("an all-false declared-absent column is a LIMITATION, not a failure",
      _st.get("raw_rna") == "limitation", str(_st.get("raw_rna")))
check("and so is the second one", _st.get("peak_mtx") == "limitation")
check("the manifest as a whole does not FAIL", au["counts"]["fail"] == 0,
      str(au["counts"]))
check("the verdict says valid-with-limitations, not failed",
      "PASS WITH LIMITATIONS" in au["verdict"])
check("and explains that a limitation is not a defect in the file",
      "not defects in the file" in au["verdict"]
      or "NOT a malformed manifest" in str(au["checks"]))
# A single-timepoint group is a study-design limitation, not a bad file.
check("a group present at one timepoint is flagged as a limitation",
      _st.get("single-timepoint group: BJFF") == "limitation")
# The distinction must still catch REAL malformation.
_bad = root / "Docs" / "bad.csv"
_bad.write_text("sample_id,rna_gsm,atac_gsm\nS1,G1,A1\nS1,G1,A2\n")
au2 = ar.audit_manifest("Docs/bad.csv", unique="sample_id,rna_gsm",
                         pair="rna_gsm:atac_gsm")
check("duplicate ids are still a FAIL", au2["counts"]["fail"] >= 1)
check("and the verdict says so", au2["verdict"].startswith("FAIL"))
# require_true is the opt-in for columns where false really is a defect.
au3 = ar.audit_manifest("Docs/manifest.csv", require_true="raw_rna")
check("require_true turns the same column into a failure",
      any(c["status"] == "fail" for c in au3["checks"]))
check("so the caller chooses which absences matter",
      au["counts"]["fail"] == 0 and au3["counts"]["fail"] >= 1)

print(f"\n{len(FAILURES)} failure(s)")
_tmp.cleanup()
sys.exit(1 if FAILURES else 0)
