"""Annotate a variant list in any common notation, and land it in the KG.

Three gaps made "functionally annotate these variants and add them to the
knowledge graph" impossible before this skill, and the agent's own
user-extension directory is the evidence: it had authored
``run_create_variant_csv`` and ``annotate_kg_variants_tool`` for itself,
because —

1. ``annotate_variant_list`` requires a **CSV file path**. A list pasted into
   chat cannot reach it: the agent has no tool that writes a data file.
2. Its identifier parsing reads **separate CSV columns** (``chrom``, ``pos``,
   ``ref``, ``alt``). The single-token notations people actually paste —
   ``2-21001846-G-A``, ``chr2:21001846:G>A`` — parse as nothing.
3. Annotation results were never written to the local knowledge graph, so
   "add them as vertices" had no implementation at all.

This skill accepts variants **inline or from a file**, in any of the notations
below, normalises them to the Catalog's identifier form, annotates them
through the existing pipeline, and upserts each one into
``Data/KG/local_kg.sqlite`` as a ``variant`` node with edges to whatever the
annotation resolved (genes, diseases, phenotypes).

Accepted notations (auto-detected per token, mixed lists are fine)::

    2-21001846-G-A            chr2-21001846-G-A
    2:21001846:G:A            chr2:21001846:G:A
    2_21001846_G_A            chr2:g.21001846G>A
    chr2:21001846G>A          rs763341676
    NC_000002.12:21001845:G:A (SPDI, already 0-based)
    2 21001846 . G A          (VCF record)

Coordinates are treated as **1-based** (the near-universal convention for
``chr-pos-ref-alt`` text) and converted to the Catalog's 0-based
``variant_id``, matching ``annotate_variant_list.infer_variant_params``. SPDI
input is passed through unchanged because SPDI is already 0-based — silently
re-decrementing it would shift every position by one.

Usage::

    igvfagent variant-list annotate --variants "2-21001846-G-A rs763341676"
    igvfagent variant-list annotate --input my_variants.txt --label mystudy
    igvfagent variant-list parse    --variants "..."      # dry-run the parser
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from pathlib import Path

try:
    from igvfagent import annotate_variant_list as avl
    from igvfagent import _localstore as ls
    from igvfagent import ccre_linkage_annotation_skills as ccre
except Exception:  # pragma: no cover - direct-script execution
    import annotate_variant_list as avl  # type: ignore
    import _localstore as ls  # type: ignore
    import ccre_linkage_annotation_skills as ccre  # type: ignore

__all__ = ["main", "parse_variant_token", "parse_variants"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
OUT_DIR = ROOT / "Docs" / "VariantList"

_RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
# SPDI: RefSeq accession : position : ref : alt  (position already 0-based)
_SPDI_RE = re.compile(r"^(N[CGTWM]_\d+\.\d+):(\d+):([ACGTN]*):([ACGTN]*)$", re.IGNORECASE)
# chr-pos-ref-alt with -, :, or _ separators; optional chr prefix; optional
# g. and > forms.  chr2:g.21001846G>A / chr2:21001846G>A / 2-21001846-G-A
_DELIM_RE = re.compile(
    r"^(?:chr)?([0-9]{1,2}|[XYM]|MT)[-:_]g?\.?(\d+)[-:_]([ACGTN]+)[-:_>]([ACGTN]+)$",
    re.IGNORECASE)
_HGVS_RE = re.compile(
    r"^(?:chr)?([0-9]{1,2}|[XYM]|MT)[-:_]g?\.?(\d+)([ACGTN]+)>([ACGTN]+)$",
    re.IGNORECASE)


class VariantParseError(ValueError):
    pass


class _SkipCatalog(Exception):
    """Internal: --sources excluded the Catalog."""


def _norm_chrom(c: str) -> str:
    c = c.upper().lstrip("CHR")
    return "MT" if c == "M" else c


def parse_variant_token(token: str) -> dict:
    """Parse one variant in any supported notation.

    Returns a dict with ``raw``, ``kind``, and the Catalog query parameters
    (``rsid`` or ``variant_id``). Raises VariantParseError if unrecognised —
    callers report those rather than silently dropping them, because a
    silently skipped variant is indistinguishable from one with no evidence.
    """
    t = (token or "").strip().strip(",;\"'")
    if not t:
        raise VariantParseError("empty token")

    if _RSID_RE.match(t):
        return {"raw": t, "kind": "rsid", "rsid": t.lower()}

    m = _SPDI_RE.match(t)
    if m:
        # SPDI positions are already 0-based — pass through untouched.
        return {"raw": t, "kind": "spdi", "variant_id": t}

    for rx, kind in ((_DELIM_RE, "chr-pos-ref-alt"), (_HGVS_RE, "hgvs-g")):
        m = rx.match(t)
        if m:
            chrom, pos, ref, alt = m.groups()
            try:
                pos0 = int(pos) - 1          # 1-based input -> 0-based Catalog
            except ValueError as e:
                raise VariantParseError(f"bad position in {t!r}") from e
            if pos0 < 0:
                raise VariantParseError(f"position must be >= 1 in {t!r}")
            return {
                "raw": t, "kind": kind,
                "chrom": _norm_chrom(chrom), "pos": int(pos),
                "ref": ref.upper(), "alt": alt.upper(),
                "variant_id": f"chr{_norm_chrom(chrom)}:{pos0}:{ref.upper()}:{alt.upper()}",
            }

    raise VariantParseError(f"unrecognised variant notation: {t!r}")


def _vcf_records(text: str) -> "list[str]":
    """Pull CHROM/POS/REF/ALT out of VCF-style whitespace records."""
    out = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) >= 5 and f[1].isdigit() and re.fullmatch(r"[ACGTN]+", f[3], re.I):
            out.append(f"{f[0]}-{f[1]}-{f[3]}-{f[4]}")
    return out


def _tokenize(text: str) -> "list[str]":
    """Split a blob into variant tokens, line by line.

    Per-line rather than whole-blob, because a file may legitimately mix
    notations: an earlier version detected a VCF-shaped line and then parsed
    *only* VCF records, silently discarding every rsID and chr-pos-ref-alt in
    the same file. Comments and blank lines are dropped here so that a `#`
    header does not surface as a batch of unparseable tokens.
    """
    out: "list[str]" = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()        # strip comments
        if not line:
            continue
        vcf = _vcf_records(line)                     # a VCF record for THIS line
        out.extend(vcf if vcf else
                   [t for t in re.split(r"[\s,]+", line) if t])
    return out


def parse_variants(text: str) -> "tuple[list[dict], list[dict]]":
    """Parse a free-form blob. Returns (parsed, failures)."""
    tokens = _tokenize(text)

    parsed: "list[dict]" = []
    failures: "list[dict]" = []
    seen: "set[str]" = set()
    for tok in tokens:
        if not tok.strip():
            continue
        try:
            rec = parse_variant_token(tok)
        except VariantParseError as e:
            failures.append({"raw": tok, "error": str(e)})
            continue
        key = rec.get("variant_id") or rec.get("rsid")
        if key in seen:
            continue
        seen.add(key)
        parsed.append(rec)
    return parsed, failures


def _favor_query(rec: dict) -> "tuple[str, str] | None":
    """(kind, value) for FAVOR, or None if this record cannot address it.

    FAVOR's /v1/variants takes the **1-based** `chr-pos-ref-alt` form — the
    exact notation people paste — so the raw coordinates are used here, not
    the 0-based Catalog variant_id. Passing the Catalog id would silently
    query the wrong position.
    """
    if rec["kind"] == "rsid":
        return ("rsid", rec["rsid"])
    if rec.get("chrom") and rec.get("pos"):
        return ("variant", f"{rec['chrom']}-{rec['pos']}-{rec['ref']}-{rec['alt']}")
    return None


def _annotate_favor(rec: dict) -> dict:
    q = _favor_query(rec)
    if not q:
        return {"favor_status": "no_favor_identifier"}
    try:
        status, data = ccre.fetch_favor_variant(*q)
    except Exception as e:
        return {"favor_status": f"error: {type(e).__name__}"}
    if status != 200:
        return {"favor_status": f"http_{status}"}
    flat = ccre.flatten_favor(data)
    flat["favor_status"] = "ok"
    return flat


def _to_csv_row(rec: dict) -> dict:
    """Shape a parsed record into the column names annotate_row expects."""
    if rec["kind"] == "rsid":
        return {"rsID": rec["rsid"]}
    if rec["kind"] == "spdi":
        return {"SPDI": rec["variant_id"]}
    return {"chrom": rec["chrom"], "pos": str(rec["pos"]),
            "ref": rec["ref"], "alt": rec["alt"]}


def _kg_ingest(records: "list[dict]", *, label: str) -> dict:
    """Upsert annotated variants into the local KG as nodes + edges."""
    con = ls._connect()
    added = {"variant_nodes": 0, "gene_edges": 0, "disease_edges": 0,
             "regulatory_edges": 0}
    try:
        run_node = ls.upsert_node(
            con, "analysis", f"variant_list_{label}",
            source="igvfagent:variant-list",
            label=f"variant-list annotate {label}",
            properties={"skill": "variant-list", "n_variants": len(records)})

        for rec in records:
            name = rec.get("rsid") or rec.get("variant_id")
            props = {k: v for k, v in rec.items()
                     if k in ("kind", "chrom", "pos", "ref", "alt", "raw")}
            props.update({k: v for k, v in (rec.get("annotation") or {}).items()
                          if v not in ("", None)})
            vid = ls.upsert_node(con, "variant", name,
                                 source="igvfagent:variant-list",
                                 properties=props)
            ls.upsert_edge(con, run_node, vid, "analyzed",
                            source="igvfagent:variant-list")
            added["variant_nodes"] += 1

            ann = rec.get("annotation") or {}
            # Entity-to-entity edges. Before this, every edge in the graph was
            # analysis->entity, so the KG was a provenance index rather than a
            # biological network. Genes/diseases/regulatory elements are drawn
            # from BOTH sources, because a FAVOR-only run must still produce
            # structure — the Catalog field names alone yielded zero edges.
            for gene in _genes_from(ann):
                gid = ls.upsert_node(con, "gene", gene,
                                     source="igvfagent:variant-list")
                ls.upsert_edge(con, vid, gid, "variant_in_gene",
                                source=ann.get("_gene_src", "annotation"))
                added["gene_edges"] += 1
            for dis in _diseases_from(ann):
                did = ls.upsert_node(con, "disease", dis,
                                     source="igvfagent:variant-list")
                ls.upsert_edge(con, vid, did, "associated_with",
                                source="clinvar/igvf_catalog")
                added["disease_edges"] += 1
            for kind, elem in _regulatory_from(ann):
                rid = ls.upsert_node(con, "regulatory_element", elem,
                                     source="favor",
                                     properties={"element_type": kind})
                ls.upsert_edge(con, vid, rid, "overlaps_element", source="favor")
                added["regulatory_edges"] = added.get("regulatory_edges", 0) + 1
        con.commit()
    finally:
        con.close()
    return added


def _genes_from(ann: dict) -> "list[str]":
    """Gene symbols from either source.

    FAVOR's `geneinfo` is NCBI-style `SYMBOL:GeneID` (e.g. `APOB:338`), so the
    id has to be stripped or every variant would create a node named
    "APOB:338" that never matches the gene node anything else creates.
    """
    out: "list[str]" = []
    for raw in _split_multi(ann.get("favor_geneinfo", "")):
        sym = raw.split(":")[0].strip()
        if sym:
            out.append(sym)
    out += _split_multi(ann.get("igvf_kg_qtl_genes", ""))
    seen, uniq = set(), []
    for g in out:
        if g.upper() not in seen:
            seen.add(g.upper()); uniq.append(g)
    return uniq


def _diseases_from(ann: dict) -> "list[str]":
    """ClinVar disease names (FAVOR `clndn`) plus Catalog phenotypes.

    ClinVar joins names with `|` and uses underscores for spaces; both are
    normalised so `Familial_hypercholesterolemia` and the Catalog's
    `Familial hypercholesterolemia` collapse to one node.
    """
    out = []
    for raw in _split_multi(ann.get("favor_clndn", "")):
        name = raw.replace("_", " ").strip()
        if name and name.lower() not in ("not provided", "not specified", "na"):
            out.append(name)
    out += [d.replace("_", " ") for d in _split_multi(ann.get("igvf_kg_phenotypes", ""))]
    seen, uniq = set(), []
    for d in out:
        if d.lower() not in seen:
            seen.add(d.lower()); uniq.append(d)
    return uniq


def _regulatory_from(ann: dict) -> "list[tuple[str, str]]":
    """Regulatory elements FAVOR reports the variant overlapping."""
    out = []
    for field, kind in (("favor_cage_enhancer", "CAGE_enhancer"),
                         ("favor_cage_promoter", "CAGE_promoter"),
                         ("favor_genehancer", "GeneHancer"),
                         ("favor_super_enhancer", "super_enhancer")):
        val = str(ann.get(field, "") or "").strip()
        if val and val.lower() not in ("na", "none", "nan", "0", ""):
            out.append((kind, f"{kind}:{val}"[:120]))
    return out


def _split_multi(value: str) -> "list[str]":
    if not value:
        return []
    return [v.strip() for v in re.split(r"[;,|]", str(value)) if v.strip()][:25]


def cmd_parse(args) -> int:
    text = _read_input(args)
    parsed, failures = parse_variants(text)
    print(json.dumps({"parsed": parsed, "failures": failures,
                      "n_parsed": len(parsed), "n_failed": len(failures)},
                     indent=2))
    return 0 if parsed else 2


def _read_input(args) -> str:
    if getattr(args, "input", None):
        p = Path(args.input).expanduser()
        if not p.is_absolute():
            p = ROOT / p
        if not p.is_file():
            raise SystemExit(f"error: no such file: {p}")
        return p.read_text(encoding="utf-8", errors="replace")
    if getattr(args, "variants", None):
        return args.variants
    raise SystemExit("error: provide --variants or --input")


def cmd_annotate(args) -> int:
    text = _read_input(args)
    parsed, failures = parse_variants(text)
    if not parsed:
        print("error: no variants could be parsed", file=sys.stderr)
        for f in failures[:10]:
            print(f"  {f['raw']}: {f['error']}", file=sys.stderr)
        return 2

    label = args.label or "variants"
    ts = time.strftime("%Y%m%d_%H%M%S")
    run_dir = OUT_DIR / f"{ts}_{avl.safe_label(label)}"
    run_dir.mkdir(parents=True, exist_ok=True)

    sources = {x.strip().lower() for x in (args.sources or "favor,catalog").split(",")}
    limit = args.max_rows if args.max_rows and args.max_rows > 0 else len(parsed)
    # The constant is DEFAULT_ENDPOINTS; an empty dict here would make
    # every annotation silently return nothing, so fail loudly instead.
    endpoints = getattr(avl, "DEFAULT_ENDPOINTS", None)
    if not endpoints:
        raise SystemExit("error: annotate_variant_list.DEFAULT_ENDPOINTS "
                          "not found — cannot annotate")

    annotated = []
    for rec in parsed[:limit]:
        row = _to_csv_row(rec)
        annotation = {}
        try:
            if "catalog" not in sources:
                raise _SkipCatalog()
            annotation, _evidence = avl.annotate_row(
                row, endpoints=endpoints, use_cache=not args.no_cache,
                sleep_seconds=args.sleep)
        except _SkipCatalog:
            annotation = {"igvf_kg_query_status": "skipped"}
        except Exception as e:                      # network/API failure
            annotation = {"igvf_kg_query_status": f"error: {type(e).__name__}"}
        if "favor" in sources:
            annotation.update(_annotate_favor(rec))
        rec = dict(rec)
        rec["annotation"] = annotation
        annotated.append(rec)

    # Artefacts first, so a KG failure never loses the annotation work.
    csv_path = run_dir / "annotated_variants.csv"
    cols = ["raw", "kind", "variant_id", "rsid", "chrom", "pos", "ref", "alt"]
    ann_cols = sorted({k for r in annotated for k in (r.get("annotation") or {})})
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + ann_cols, extrasaction="ignore")
        w.writeheader()
        for r in annotated:
            row = {k: r.get(k, "") for k in cols}
            row.update(r.get("annotation") or {})
            w.writerow(row)

    json_path = run_dir / "annotated_variants.json"
    json_path.write_text(json.dumps(
        {"label": label, "n_input": len(parsed), "n_annotated": len(annotated),
         "failures": failures, "variants": annotated}, indent=2))

    kg = {"skipped": True}
    if not args.no_kg:
        try:
            kg = _kg_ingest(annotated, label=avl.safe_label(label))
        except Exception as e:
            kg = {"error": f"{type(e).__name__}: {e}"}

    report = run_dir / "report.md"
    report.write_text(_render_report(label, annotated, failures, kg))

    print(f"Parsed:   {len(parsed)} variant(s); {len(failures)} unparsed")
    if failures:
        for f in failures[:5]:
            print(f"  unparsed: {f['raw']}  ({f['error']})")
    print(f"KG:       {json.dumps(kg)}")
    print(f"Report:        {report}")
    print(f"Manifest:      {csv_path}")
    print(f"Wrote:         {json_path}")
    return 0


def _render_report(label, annotated, failures, kg) -> str:
    lines = [f"# Variant list annotation — {label}", "",
             f"- Variants annotated: **{len(annotated)}**",
             f"- Unparsed tokens: **{len(failures)}**",
             f"- Knowledge-graph ingest: `{json.dumps(kg)}`", "",
             "## Variants", "",
             "| input | parsed as | catalog id | status |", "|---|---|---|---|"]
    for r in annotated:
        st = (r.get("annotation") or {}).get("igvf_kg_query_status", "")
        lines.append(f"| `{r['raw']}` | {r['kind']} | "
                     f"`{r.get('variant_id') or r.get('rsid')}` | {st} |")
    if failures:
        lines += ["", "## Unparsed", ""]
        lines += [f"- `{f['raw']}` — {f['error']}" for f in failures]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="igvfagent variant-list",
        description="Annotate a variant list in any common notation and add "
                    "the results to the local knowledge graph.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    for name, help_ in (("annotate", "Annotate and ingest into the KG"),
                         ("parse", "Dry-run the parser only")):
        sp = sub.add_parser(name, help=help_)
        sp.add_argument("--variants", help="Inline list: any notation, any "
                                            "separator (space/comma/newline).")
        sp.add_argument("--input", help="Path to a file of variants "
                                         "(txt/tsv/vcf/csv).")
        sp.add_argument("--label", default="variants")
        if name == "annotate":
            sp.add_argument("--max-rows", type=int, default=0,
                            help="0 = all.")
            sp.add_argument("--no-kg", action="store_true",
                            help="Annotate but do not touch the KG.")
            sp.add_argument("--sources", default="favor,catalog",
                            help="Comma list: favor, catalog. Default both.")
            sp.add_argument("--no-cache", action="store_true")
            sp.add_argument("--sleep", type=float, default=0.1)

    args = parser.parse_args(argv)
    return cmd_parse(args) if args.cmd == "parse" else cmd_annotate(args)


if __name__ == "__main__":
    raise SystemExit(main())
