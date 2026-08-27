"""Cross-source verification of variant annotations.

An agent's conclusions inherit its sources' errors, and the most dangerous
source error is a stale snapshot: it returns a confident, well-formed answer
that is simply out of date. This skill re-checks annotations against a second,
independent source and reports where they disagree.

The motivating case is real. Annotating 25 APOB stopgain variants through
FAVOR and then re-checking ClinVar directly showed material divergence — for
rs145143533, FAVOR reported *Benign* with the trait *familial
hypercholesterolemia 2*, while current ClinVar reports **Pathogenic** with
*familial hypobetalipoproteinemia*. That is a reversal in both clinical
significance and the implied direction of effect, and nothing in the
single-source output signals it. An analysis built on the first answer would
be confidently wrong.

Two verification axes:

* **agreement** — does the second source report the same clinical
  significance and the same disease?
* **mechanistic coherence** — for a truncating (nonsense/frameshift) variant,
  is the asserted phenotype consistent with loss of function? APOB stopgains
  are expected to cause hypobetalipoproteinemia; an assertion of
  hypercholesterolemia (which arises from LDLR-binding-defective missense
  alleles) is worth flagging for review rather than silently accepting.

Disagreement is reported, never resolved. ClinVar is treated as current rather
than correct; the point is to surface the divergence for a human, not to pick
a winner automatically.

Usage::

    igvfagent variant-verify clinvar --input Data/Input/VariantList/apob_lof_rsids.txt
    igvfagent variant-verify clinvar --variants "rs145143533 rs121918390"
    igvfagent variant-verify clinvar --annotated Docs/VariantList/<run>/annotated_variants.csv

Pure standard library (urllib + json + csv).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

try:
    from igvfagent import variant_list_skill as vls
except Exception:  # pragma: no cover
    import variant_list_skill as vls  # type: ignore

__all__ = ["main", "clinvar_lookup", "compare"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
OUT_DIR = ROOT / "Docs" / "VariantVerify"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Truncating consequences, for the coherence check.
_TRUNCATING = ("stopgain", "nonsense", "frameshift", "stoploss", "startloss")

# Phenotype direction for APOB-style lipid genes. Kept explicit and narrow:
# a general phenotype-direction ontology is out of scope, and guessing would
# produce confident nonsense on genes this table does not cover.
_DIRECTION = {
    "hypobeta": "loss-of-function",
    "hypercholesterol": "not-loss-of-function",
}


def _get(url: str, retries: int = 2):
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except Exception:
            if attempt == retries:
                return None
            time.sleep(1.0 + attempt)
    return None


def clinvar_lookup(rsid: str) -> dict:
    """Current ClinVar significance and traits for one rsID.

    esearch on a bare rsID sometimes returns a UID whose esummary carries no
    result block, so the accession-scoped query is tried first and the bare
    term kept as a fallback. A lookup that fails returns status "unresolved"
    rather than an empty record, because "no data" and "no assertion" mean
    very different things here.
    """
    ids = []
    for term in (f"{rsid}[Variant ID]", rsid):
        s = _get(f"{EUTILS}/esearch.fcgi?db=clinvar&term="
                 f"{urllib.parse.quote(term)}&retmode=json")
        ids = ((s or {}).get("esearchresult") or {}).get("idlist") or []
        if ids:
            break
    if not ids:
        return {"rsid": rsid, "status": "not_in_clinvar"}

    d = _get(f"{EUTILS}/esummary.fcgi?db=clinvar&id={ids[0]}&retmode=json")
    res = (d or {}).get("result") or {}
    key = next((k for k in res if k != "uids"), None)
    if not key:
        return {"rsid": rsid, "status": "unresolved", "uid": ids[0]}

    r = res[key]
    gc = r.get("germline_classification") or {}
    traits = [t.get("trait_name") for t in (gc.get("trait_set") or [])
              if t.get("trait_name")]
    traits = [t for t in traits if t.strip().lower() != "not provided"]
    mc = r.get("molecular_consequence_list") or []
    return {
        "rsid": rsid,
        "status": "ok",
        "uid": ids[0],
        "title": r.get("title"),
        "significance": gc.get("description"),
        "traits": traits,
        "consequence": "; ".join(mc) if isinstance(mc, list) else str(mc),
        "last_evaluated": gc.get("last_evaluated"),
    }


def _direction(text: str) -> "str | None":
    low = (text or "").lower()
    hits = {v for k, v in _DIRECTION.items() if k in low}
    if len(hits) == 1:
        return hits.pop()
    return "both" if len(hits) > 1 else None


def _norm_sig(s: str) -> str:
    return re.sub(r"[^a-z]", "", (s or "").lower())


def compare(favor_row: dict, cv: dict) -> dict:
    """One variant: FAVOR vs current ClinVar."""
    f_sig = (favor_row.get("favor_clnsig") or "").replace("_", " ").strip()
    f_dis = (favor_row.get("favor_clndn") or "").replace("_", " ").strip()
    c_sig = cv.get("significance") or ""
    c_dis = "; ".join(cv.get("traits") or [])
    conseq = (favor_row.get("favor_genecode_comprehensive_exonic_category")
              or cv.get("consequence") or "")

    flags = []
    if cv.get("status") != "ok":
        flags.append(f"clinvar_{cv.get('status')}")
    else:
        if f_sig and c_sig and _norm_sig(f_sig) != _norm_sig(c_sig):
            flags.append("significance_differs")
        fd, cd = _direction(f_dis), _direction(c_dis)
        if fd and cd and fd != cd:
            flags.append("direction_of_effect_differs")
        if f_dis and c_dis and not (set(f_dis.lower().split())
                                    & set(c_dis.lower().split())):
            flags.append("disease_terms_disjoint")

    # Mechanistic coherence, independent of agreement.
    if any(t in conseq.lower() for t in _TRUNCATING):
        for label, source in (("favor", f_dis), ("clinvar", c_dis)):
            if _direction(source) in ("not-loss-of-function", "both"):
                flags.append(f"truncating_variant_with_non_LoF_phenotype[{label}]")

    return {
        "variant": favor_row.get("raw") or favor_row.get("rsid"),
        "rsid": favor_row.get("favor_rsid") or favor_row.get("rsid"),
        "consequence": conseq,
        "favor": {"significance": f_sig or None, "diseases": f_dis or None},
        "clinvar": {"significance": c_sig or None, "diseases": c_dis or None,
                    "status": cv.get("status"),
                    "last_evaluated": cv.get("last_evaluated")},
        "flags": sorted(set(flags)),
        "agree": not flags,
    }


def _rows_from_annotated(path: Path) -> "list[dict]":
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _rows_from_variants(text: str) -> "list[dict]":
    parsed, _fail = vls.parse_variants(text)
    out = []
    for rec in parsed:
        rs = rec.get("rsid")
        if not rs:
            continue                      # ClinVar lookup needs an rsID
        out.append({"raw": rec["raw"], "rsid": rs, "favor_rsid": rs})
    return out


def cmd_clinvar(args) -> int:
    if args.annotated:
        p = Path(args.annotated)
        if not p.is_absolute():
            p = ROOT / p
        rows = _rows_from_annotated(p)
    elif args.input:
        p = Path(args.input)
        if not p.is_absolute():
            p = ROOT / p
        rows = _rows_from_variants(p.read_text(encoding="utf-8", errors="replace"))
    elif args.variants:
        rows = _rows_from_variants(args.variants)
    else:
        print("error: provide --annotated, --input or --variants", file=sys.stderr)
        return 2

    rows = [r for r in rows if (r.get("favor_rsid") or r.get("rsid"))]
    if not rows:
        print("error: no rsIDs available; ClinVar verification needs rsIDs "
              "(chr-pos-ref-alt input can be mapped first with "
              "`variant-list annotate`)", file=sys.stderr)
        return 2
    if args.max_rows:
        rows = rows[: args.max_rows]

    results = []
    for r in rows:
        rs = r.get("favor_rsid") or r.get("rsid")
        cv = clinvar_lookup(rs)
        results.append(compare(r, cv))
        time.sleep(args.sleep)

    ts = time.strftime("%Y%m%d_%H%M%S")
    run = OUT_DIR / f"{ts}_{re.sub(r'[^A-Za-z0-9_.-]', '_', args.label)}"
    run.mkdir(parents=True, exist_ok=True)
    (run / "verification.json").write_text(json.dumps(
        {"label": args.label, "n": len(results), "results": results}, indent=2))

    disagree = [r for r in results if not r["agree"]]
    report = run / "report.md"
    report.write_text(_render(args.label, results, disagree))

    print(f"Verified:      {len(results)} variant(s) against current ClinVar")
    print(f"Disagreements: {len(disagree)}")
    for r in disagree[:10]:
        print(f"  {str(r['rsid']):14s} {','.join(r['flags'])}")
    print(f"Report:        {report}")
    print(f"Wrote:         {run / 'verification.json'}")
    # Non-zero when anything diverges, so a pipeline can gate on it.
    return 1 if disagree else 0


def _render(label, results, disagree) -> str:
    L = [f"# Variant verification — {label}", "",
         f"- Variants checked: **{len(results)}**",
         f"- Disagreements with current ClinVar: **{len(disagree)}**", "",
         "ClinVar is treated as *current*, not as *correct*. Divergence is "
         "reported for human review, not resolved automatically.", "",
         "## Disagreements", "",
         "| variant | rsid | FAVOR | current ClinVar | flags |",
         "|---|---|---|---|---|"]
    for r in disagree:
        f, c = r["favor"], r["clinvar"]
        L.append(f"| `{r['variant']}` | {r['rsid']} | "
                 f"{f['significance'] or '-'} / {(f['diseases'] or '-')[:40]} | "
                 f"{c['significance'] or '-'} / {(c['diseases'] or '-')[:40]} | "
                 f"{', '.join(r['flags'])} |")
    if not disagree:
        L.append("| — | — | — | — | none |")
    L += ["", "## All variants", "",
          "| variant | consequence | FAVOR sig | ClinVar sig | agree |",
          "|---|---|---|---|---|"]
    for r in results:
        L.append(f"| `{r['variant']}` | {r['consequence'] or '-'} | "
                 f"{r['favor']['significance'] or '-'} | "
                 f"{r['clinvar']['significance'] or '-'} | "
                 f"{'yes' if r['agree'] else '**no**'} |")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent variant-verify",
        description="Re-check variant annotations against an independent "
                    "source and report divergence.")
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("clinvar", help="Verify against current ClinVar")
    c.add_argument("--annotated", help="annotated_variants.csv from variant-list")
    c.add_argument("--input", help="File of variants (rsIDs)")
    c.add_argument("--variants", help="Inline variants (rsIDs)")
    c.add_argument("--label", default="verify")
    c.add_argument("--max-rows", type=int, default=0)
    c.add_argument("--sleep", type=float, default=0.4)
    args = p.parse_args(argv)
    return cmd_clinvar(args)


if __name__ == "__main__":
    raise SystemExit(main())
