"""Tissue-filtered eQTLs, and MPRA-vs-eQTL concordance.

Two capability gaps a user hit trying to compare HepG2 allelic MPRA against
GTEx-liver fine-mapped eQTLs: there was no way to restrict eQTLs to a tissue,
and nothing computed correlation or concordance between the two assays.

**Tissue filtering is client-side, and that is not a shortcut.** The Catalog's
`/api/genes/variants` endpoint accepts `biological_context=liver` and returns
zero rows, and accepts `study=studies/QTS000015` and returns rows from other
studies — both verified against live data. Filtering therefore happens here,
after paging, and the row counts always report how many were seen before
filtering so a narrow result is never mistaken for a complete one.

**Multi-gene variants.** A fine-mapped variant is frequently credible for
several genes. `--pip-rule max` keeps the single highest-PIP gene per variant,
which is what the user asked for; `--pip-rule all` keeps every pair.

**Concordance.** Once MPRA allelic effects and eQTL effect sizes are joined on
a normalised variant id, the skill reports Pearson and Spearman correlation
plus a sign-concordance rate — the fraction of shared variants whose MPRA and
eQTL effects point the same way. Sign concordance is the more robust statistic
here: the two assays measure different things on different scales, so
agreement in direction is meaningful where agreement in magnitude may not be.

Correlations are reported with n and a bootstrap CI, because the intersection
of an MPRA library with one tissue's credible sets is usually small, and a
correlation from a handful of variants is noise wearing a decimal point.

Usage::

    igvfagent eqtl-mpra tissues --gene APOB
    igvfagent eqtl-mpra eqtl --genes APOB,PCSK9 --tissue liver --pip-rule max
    igvfagent eqtl-mpra concordance --genes APOB,PCSK9,LDLR \\
        --tissue liver --biosample HepG2 --significant-only

Pure standard library: urllib, csv, json, statistics, random.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import statistics
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

__all__ = ["main", "fetch_eqtls", "normalise_variant", "concordance"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
OUT_DIR = ROOT / "Docs" / "eQTLMPRA"
CATALOG = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")

# The MPRA variant-biosample collection runs to at least 20,000 rows and is
# coordinate-ordered, so a shallow page cap silently confines every result to
# chr1. Raised, with the cap reported whenever it binds — a truncated pull
# that looks complete is what made the first runs return zero overlap.
_MAX_PAGES = 400


def _get(path: str, **params):
    url = f"{CATALOG}{path}?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return json.load(r)
    except Exception:
        return []


def _paged(path: str, *, max_rows: int, page_size: int = 100, **params):
    """Page until exhausted or max_rows. Pages genuinely advance here."""
    rows = []
    for page in range(_MAX_PAGES):
        if len(rows) >= max_rows:
            return rows, False              # caller's row budget reached
        batch = _get(path, limit=str(min(page_size, max_rows - len(rows))),
                     page=str(page), verbose="false", **params)
        if not isinstance(batch, list) or not batch:
            return rows, False              # empty page: nothing left
        rows.extend(batch)
        if len(batch) < page_size:
            return rows, False              # short page: collection exhausted
    # Fell out of the loop: the page cap bound before the collection ended, so
    # this pull is truncated and any "n of raw n" is the cap, not the total.
    return rows, len(rows) < max_rows


# ---------------------------------------------------------------------------
# Variant identity
# ---------------------------------------------------------------------------

_CHR_BY_ACC = {f"NC_{i:06d}": (str(i) if i < 23 else {23: "X", 24: "Y"}.get(i, str(i)))
               for i in range(1, 25)}


def normalise_variant(value: str) -> str:
    """A comparable key for a variant across sources.

    MPRA and eQTL records both use RefSeq-accession SPDI
    (`NC_000002.12:21221274:T:C`), sometimes prefixed `variants/`. Reducing to
    `chr:pos:ref:alt` lets the two tables join. Anything unrecognised is
    returned stripped rather than dropped, so a join miss is visible as a
    non-matching key instead of a silently discarded row.
    """
    v = (value or "").strip()
    if v.startswith("variants/"):
        v = v.split("/", 1)[1]
    parts = v.split(":")
    if len(parts) == 4:
        acc, pos, ref, alt = parts
        chrom = _CHR_BY_ACC.get(acc.split(".")[0])
        if chrom:
            return f"chr{chrom}:{pos}:{ref.upper()}:{alt.upper()}"
    return v


# ---------------------------------------------------------------------------
# eQTL retrieval
# ---------------------------------------------------------------------------

def fetch_eqtls(genes, *, tissue: str = "", study: str = "",
                max_rows: int = 2000, pip_rule: str = "max",
                min_pip: float = 0.0) -> dict:
    """Fine-mapped eQTLs for genes, optionally restricted to one tissue."""
    seen_ctx: "dict[str, int]" = {}
    kept, n_raw = [], 0
    want_ctx = (tissue or "").strip().lower()
    want_study = (study or "").strip().lower()

    for gene in genes:
        rows, _trunc = _paged("/api/genes/variants", gene_name=gene,
                              max_rows=max_rows)
        for r in rows:
            if str(r.get("label") or r.get("method")).lower() != "eqtl":
                continue
            n_raw += 1
            ctx = str(r.get("biological_context") or "")
            seen_ctx[ctx] = seen_ctx.get(ctx, 0) + 1
            # Client-side: the endpoint ignores biological_context= and study=.
            if want_ctx and want_ctx not in ctx.lower():
                continue
            if want_study and want_study not in str(r.get("study") or "").lower():
                continue
            pip = _f(r.get("posterior_inclusion_probability"))
            if pip is not None and pip < min_pip:
                continue
            kept.append({
                "variant": normalise_variant(str(r.get("sequence_variant") or "")),
                "raw_variant": str(r.get("sequence_variant") or ""),
                "gene": str(r.get("gene") or "").split("/")[-1],
                "gene_query": gene,
                "tissue": ctx,
                "study": str(r.get("study") or ""),
                "pip": pip,
                "effect_size": _f(r.get("effect_size")),
                "z_score": _f(r.get("z_score")),
                "p_value": _f(r.get("p_value")),
                "source_url": str(r.get("source_url") or ""),
            })

    collapsed = _collapse_by_pip(kept) if pip_rule == "max" else kept
    return {"rows": collapsed, "n_raw_eqtl": n_raw, "n_after_filter": len(kept),
            "n_after_pip_rule": len(collapsed),
            "contexts_seen": dict(sorted(seen_ctx.items(),
                                          key=lambda kv: -kv[1])[:25]),
            "tissue": tissue, "study": study, "pip_rule": pip_rule}


def _collapse_by_pip(rows):
    """Highest-PIP gene per variant.

    A fine-mapped variant is often credible for several genes; keeping every
    pair would let one variant contribute repeatedly to a correlation and
    overstate n.
    """
    best: "dict[str, dict]" = {}
    for r in rows:
        k = r["variant"]
        cur = best.get(k)
        if cur is None or (r["pip"] or -1) > (cur["pip"] or -1):
            best[k] = r
    return list(best.values())


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# MPRA retrieval
# ---------------------------------------------------------------------------

def fetch_mpra(*, biosample: str = "", significant_only: bool = False,
               max_rows: int = 5000) -> dict:
    """Variant-level allelic MPRA effects, optionally one biosample."""
    terms = _resolve_ontology(biosample) if biosample else set()
    rows, truncated = _paged("/api/variants/biosamples", method="MPRA",
                             max_rows=max_rows)
    n_raw = len(rows)
    out = []
    for r in rows:
        bios = str(r.get("biosample") or "").split("/")[-1]
        if terms and bios not in terms:
            continue
        if significant_only and str(r.get("significant", "")).lower() not in ("true", "1"):
            continue
        out.append({
            # The field is `variant`. `name` is a relationship label
            # ("modulates regulatory activity of") and falling back to it
            # silently produced one bogus key for every row, so every join
            # returned zero — which reads as "no shared variants" rather than
            # as a bug. Never fall back to a field that is not an identifier.
            "variant": normalise_variant(str(r.get("variant")
                                              or r.get("sequence_variant") or "")),
            "genomic_element": str(r.get("genomic_element") or "").split("/")[-1],
            "postProbEffect": _f(r.get("postProbEffect")),
            "biosample": bios,
            "log2FC": _f(r.get("log2FC")),
            "neg_log10_p": _f(r.get("neg_log10_pvalue")),
            "significant": str(r.get("significant", "")),
        })
    return {"rows": [r for r in out if r["variant"]], "n_raw_mpra": n_raw,
            "truncated_at_page_cap": truncated,
            "biosample": biosample, "ontology_terms": sorted(terms)}


def _resolve_ontology(spec: str) -> "set[str]":
    tail = (spec or "").split("/")[-1].replace(":", "_")
    if tail and tail[0].isalpha() and "_" in tail and tail.split("_")[-1].isdigit():
        return {tail}
    out = set()
    for row in _get("/api/ontology-terms", name=spec.lower(), limit="25") or []:
        if row.get("term_id"):
            out.add(str(row["term_id"]))
    return out


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def _pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else None


def _rank(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _spearman(xs, ys):
    return _pearson(_rank(xs), _rank(ys)) if len(xs) >= 3 else None


def _boot_ci(xs, ys, fn, *, iters=2000, seed=0):
    """Percentile bootstrap CI. Reported because these intersections are small
    and a bare r from a handful of variants invites over-reading."""
    if len(xs) < 5:
        return None
    rng = random.Random(seed)
    stats = []
    n = len(xs)
    for _ in range(iters):
        idx = [rng.randrange(n) for _ in range(n)]
        v = fn([xs[i] for i in idx], [ys[i] for i in idx])
        if v is not None:
            stats.append(v)
    if len(stats) < iters // 4:
        return None
    stats.sort()
    return [round(stats[int(0.025 * len(stats))], 4),
            round(stats[int(0.975 * len(stats)) - 1], 4)]


def concordance(mpra_rows, eqtl_rows) -> dict:
    """Join on variant and compare MPRA allelic effect with eQTL effect size."""
    eq = {r["variant"]: r for r in eqtl_rows if r.get("variant")}
    pairs = []
    for m in mpra_rows:
        e = eq.get(m["variant"])
        if not e:
            continue
        if m.get("log2FC") is None or e.get("effect_size") is None:
            continue
        pairs.append((m, e))

    # Diagnose a zero/small overlap rather than reporting a bare 0, which
    # reads as "these assays disagree" when the real cause is that the two
    # row sets never covered the same chromosomes.
    def _chrs(rows):
        return {r["variant"].split(":")[0] for r in rows
                if r.get("variant") and ":" in r["variant"]}
    mchr, echr = _chrs(mpra_rows), _chrs(eqtl_rows)
    coverage = {"mpra_chromosomes": sorted(mchr)[:8],
                "eqtl_chromosomes": sorted(echr)[:8],
                "chromosomes_in_common": sorted(mchr & echr)[:8]}

    xs = [m["log2FC"] for m, _ in pairs]
    ys = [e["effect_size"] for _, e in pairs]
    same = sum(1 for a, b in zip(xs, ys) if a * b > 0)
    nonzero = sum(1 for a, b in zip(xs, ys) if a * b != 0)

    return {
        "n_mpra": len(mpra_rows), "n_eqtl": len(eqtl_rows),
        "n_shared_variants": len(pairs),
        "pearson_r": round(_pearson(xs, ys), 4) if _pearson(xs, ys) is not None else None,
        "pearson_ci95": _boot_ci(xs, ys, _pearson),
        "spearman_rho": round(_spearman(xs, ys), 4) if _spearman(xs, ys) is not None else None,
        "spearman_ci95": _boot_ci(xs, ys, _spearman),
        "sign_concordance": round(same / nonzero, 4) if nonzero else None,
        "n_sign_evaluable": nonzero,
        "coverage": coverage,
        "no_overlap_reason": (
            "the two row sets share no chromosome — widen the gene set or "
            "raise --mpra-max-rows; the MPRA collection is coordinate-ordered "
            "and cannot be filtered server-side"
            if not (mchr & echr) and mpra_rows and eqtl_rows else None),
        "pairs": [{"variant": m["variant"], "gene": e["gene"],
                   "tissue": e["tissue"], "mpra_log2FC": m["log2FC"],
                   "eqtl_effect_size": e["effect_size"], "eqtl_pip": e["pip"],
                   "same_direction": (m["log2FC"] * e["effect_size"]) > 0}
                  for m, e in pairs],
    }


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _run_dir(label: str) -> Path:
    d = OUT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cmd_tissues(args) -> int:
    """What tissues exist for these genes — answerable before filtering."""
    res = fetch_eqtls([g.strip() for g in args.genes.split(",") if g.strip()],
                      max_rows=args.max_rows, pip_rule="all")
    print(f"eQTL rows scanned: {res['n_raw_eqtl']}")
    print("Available biological_context values:")
    for ctx, n in res["contexts_seen"].items():
        print(f"  {n:>5}  {ctx}")
    if not res["contexts_seen"]:
        print("  (none — these genes have no eQTL edges in the Catalog)")
    return 0


def cmd_eqtl(args) -> int:
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    res = fetch_eqtls(genes, tissue=args.tissue, study=args.study,
                      max_rows=args.max_rows, pip_rule=args.pip_rule,
                      min_pip=args.min_pip)
    d = _run_dir(args.label)
    _write_csv(d / "eqtls.csv", res["rows"])
    (d / "summary.json").write_text(json.dumps(
        {k: v for k, v in res.items() if k != "rows"}, indent=2))
    print(f"eQTL rows (raw):        {res['n_raw_eqtl']}")
    print(f"After tissue/study:     {res['n_after_filter']}"
          + (f"  [tissue={args.tissue!r}]" if args.tissue else ""))
    print(f"After --pip-rule {args.pip_rule}:  {res['n_after_pip_rule']}")
    if args.tissue and not res["n_after_filter"]:
        print(f"NOTE: no rows matched tissue {args.tissue!r}. Contexts present: "
              + ", ".join(list(res["contexts_seen"])[:8]))
    print(f"Manifest:      {d / 'eqtls.csv'}")
    print(f"Wrote:         {d / 'summary.json'}")
    return 0


def cmd_concordance(args) -> int:
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    eq = fetch_eqtls(genes, tissue=args.tissue, study=args.study,
                     max_rows=args.max_rows, pip_rule=args.pip_rule,
                     min_pip=args.min_pip)
    mp = fetch_mpra(biosample=args.biosample,
                    significant_only=args.significant_only,
                    max_rows=args.mpra_max_rows)
    res = concordance(mp["rows"], eq["rows"])

    d = _run_dir(args.label)
    _write_csv(d / "eqtls.csv", eq["rows"])
    _write_csv(d / "mpra.csv", mp["rows"])
    _write_csv(d / "shared_variants.csv", res["pairs"])
    meta = {"genes": genes, "tissue": args.tissue, "biosample": args.biosample,
            "significant_only": args.significant_only,
            "pip_rule": args.pip_rule,
            "eqtl": {k: v for k, v in eq.items() if k != "rows"},
            "mpra": {k: v for k, v in mp.items() if k != "rows"},
            "stats": {k: v for k, v in res.items() if k != "pairs"}}
    (d / "summary.json").write_text(json.dumps(meta, indent=2))
    (d / "report.md").write_text(_render(meta, res))

    if mp.get("truncated_at_page_cap"):
        print(f"WARNING: the MPRA pull hit the page cap at {mp['n_raw_mpra']} "
              f"rows — the collection is larger. Raise --mpra-max-rows; the "
              f"rows retrieved are the FIRST by coordinate, not a sample.")
    print(f"MPRA rows:              {res['n_mpra']} (raw {mp['n_raw_mpra']}"
          + (f", biosample={mp['ontology_terms']}" if mp['ontology_terms'] else "") + ")")
    print(f"eQTL rows:              {res['n_eqtl']} (raw {eq['n_raw_eqtl']}"
          + (f", tissue={args.tissue!r}" if args.tissue else "") + ")")
    print(f"Shared variants:        {res['n_shared_variants']}")
    if res["n_shared_variants"] < 5:
        cov = res.get("coverage") or {}
        print("NOTE: too few shared variants for a meaningful correlation.")
        print(f"      MPRA chromosomes: {cov.get('mpra_chromosomes')}")
        print(f"      eQTL chromosomes: {cov.get('eqtl_chromosomes')}")
        if res.get("no_overlap_reason"):
            print(f"      Cause: {res['no_overlap_reason']}")
            print("      The Catalog serves MPRA variant-biosample rows in "
                  "coordinate order and honours no filter on this endpoint "
                  "(biosample=, variant= are accepted and ignored), so "
                  "reaching a distant locus requires paging the whole "
                  "collection. Prefer --genes on chromosomes the MPRA pages "
                  "already cover, or supply a local MPRA table.")
    else:
        print(f"Pearson r:              {res['pearson_r']}  CI95={res['pearson_ci95']}")
        print(f"Spearman rho:           {res['spearman_rho']}  CI95={res['spearman_ci95']}")
        print(f"Sign concordance:       {res['sign_concordance']} "
              f"(n={res['n_sign_evaluable']})")
    print(f"Report:        {d / 'report.md'}")
    print(f"Manifest:      {d / 'shared_variants.csv'}")
    return 0


def _write_csv(path: Path, rows) -> None:
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _render(meta, res) -> str:
    s = meta["stats"]
    L = [f"# MPRA vs eQTL concordance", "",
         f"- Genes: `{', '.join(meta['genes'])}`",
         f"- eQTL tissue filter: `{meta['tissue'] or '(none)'}`",
         f"- MPRA biosample: `{meta['biosample'] or '(none)'}`"
         + (f" -> {meta['mpra']['ontology_terms']}" if meta['mpra'].get('ontology_terms') else ""),
         f"- Multi-gene rule: `{meta['pip_rule']}`", "",
         "## Counts", "",
         f"| stage | n |", "|---|---|",
         f"| eQTL rows retrieved | {meta['eqtl']['n_raw_eqtl']} |",
         f"| after tissue/study filter | {meta['eqtl']['n_after_filter']} |",
         f"| after PIP rule | {meta['eqtl']['n_after_pip_rule']} |",
         f"| MPRA rows retrieved | {meta['mpra']['n_raw_mpra']} |",
         f"| MPRA after filters | {s['n_mpra']} |",
         f"| **shared variants** | **{s['n_shared_variants']}** |", "",
         "## Concordance", ""]
    if s["n_shared_variants"] < 5:
        L += ["Too few shared variants for a meaningful statistic. The "
              "intersection of an MPRA library with one tissue's credible sets "
              "is often small; widen the gene set before interpreting.", ""]
    else:
        L += [f"| statistic | value | 95% CI |", "|---|---|---|",
              f"| Pearson r | {s['pearson_r']} | {s['pearson_ci95']} |",
              f"| Spearman rho | {s['spearman_rho']} | {s['spearman_ci95']} |",
              f"| Sign concordance | {s['sign_concordance']} | n={s['n_sign_evaluable']} |",
              "",
              "Sign concordance is the more robust statistic: the two assays "
              "measure different quantities on different scales, so agreement "
              "in direction carries meaning where agreement in magnitude may "
              "not.", ""]
    L += ["## Caveats", "",
          "- Tissue and study filters are applied client-side. The Catalog "
          "endpoint accepts `biological_context=` and `study=` and does not "
          "honour them, so counts before and after filtering are both reported.",
          "- eQTL effect sizes come from EBI eQTL Catalogue SuSiE credible "
          "sets as ingested into the IGVF Catalog; check `source_url` for the "
          "exact release.",
          "- MPRA and eQTL variants are joined on normalised `chr:pos:ref:alt` "
          "(GRCh38). Variants that fail to normalise keep their raw id and "
          "will not join, which is visible rather than silent.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent eqtl-mpra",
        description="Tissue-filtered fine-mapped eQTLs, and MPRA-vs-eQTL "
                    "correlation / concordance.")
    sub = p.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("tissues", help="List eQTL tissues available for genes")
    t.add_argument("--genes", required=True)
    t.add_argument("--max-rows", type=int, default=1000)

    for name, helptext in (("eqtl", "Fetch tissue-filtered eQTLs"),
                            ("concordance", "MPRA vs eQTL concordance")):
        c = sub.add_parser(name, help=helptext)
        c.add_argument("--genes", required=True)
        c.add_argument("--tissue", default="",
                       help="Substring match on biological_context, e.g. liver")
        c.add_argument("--study", default="", help="e.g. studies/QTS000015")
        c.add_argument("--pip-rule", default="max", choices=("max", "all"),
                       help="max = highest-PIP gene per variant")
        c.add_argument("--min-pip", type=float, default=0.0)
        c.add_argument("--max-rows", type=int, default=2000)
        c.add_argument("--label", default=f"{name}_run")
        if name == "concordance":
            c.add_argument("--biosample", default="",
                           help="MPRA biosample name or ontology id (HepG2)")
            c.add_argument("--significant-only", action="store_true")
            c.add_argument("--mpra-max-rows", type=int, default=5000)

    args = p.parse_args(argv)
    return {"tissues": cmd_tissues, "eqtl": cmd_eqtl,
            "concordance": cmd_concordance}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
