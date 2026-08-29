"""IntPath — integrated pathway data from KEGG, WikiPathways and BioCyc.

Live pathway APIs answer "which pathways contain this gene" one database at a
time, and each names the same biology differently. IntPath (Zhou et al.) does
the harder part in advance: it normalises gene identifiers across databases,
maps each source's own relation vocabulary onto a common one, and merges
pathways that are the same pathway under different names — so a query returns
integrated evidence rather than three overlapping answers to compare by hand.

The dataset this skill loads is a two-file release:

    genes  : pathway <TAB> gene <TAB> source          (23,873 rows, 582 pathways)
    pairs  : geneA <TAB> geneB <TAB> relation <TAB> pathway <TAB> source
                                                     (71,116 rows, 5,748 genes)

The `source` field records **every** database supporting a fact, so a row
reading "KEGG WikiPathways" is corroborated evidence rather than a duplicate.
That distinction is the reason to prefer this over querying each API in turn.

Relation types follow the KEGG vocabulary IntPath normalises onto:

    PPrel  protein-protein          ECrel  enzyme-enzyme (sequential catalysis)
    GPrel  gene-product / complex   GErel  gene-expression regulation

Usage::

    igvfagent intpath load --pairs <file> --genes <file>   # ingest into the KG
    igvfagent intpath query --genes CLU,BIN1,PICALM        # pathways + partners
    igvfagent intpath status

The data files are not redistributed here. Point `--pairs` / `--genes` at a
local copy, or set ``IGVF_INTPATH_DIR``; without them the skill says so rather
than returning an empty result that looks like an answer.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

__all__ = ["main", "load_pairs", "load_memberships", "ingest"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = Path(os.environ.get("IGVF_INTPATH_DIR")
                or ROOT / "Data" / "Reference" / "IntPath")
OUT_DIR = ROOT / "Docs" / "IntPath"

# IntPath's normalised relation vocabulary -> the KG's edge types. Kept
# explicit: collapsing all four into "interacts_with" would discard the
# distinction between a physical interaction and transcriptional regulation,
# which is most of what makes the relation worth recording.
RELATION_EDGE = {
    "PPrel": ("interacts_with", "protein-protein interaction"),
    "ECrel": ("enzyme_relation", "sequential catalysis (shared metabolite)"),
    "GPrel": ("in_same_complex", "gene product / complex membership"),
    "GErel": ("regulates_expression", "gene-expression regulation"),
}


def _resolve(spec: "str | None", *names) -> "Path | None":
    if spec:
        p = Path(spec).expanduser()
        if not p.is_absolute():
            p = ROOT / spec
        return p if p.is_file() else None
    for n in names:
        p = DATA_DIR / n
        if p.is_file():
            return p
    return None


def _split_sources(field: str) -> "list[str]":
    """A source field may list several databases: 'KEGG WikiPathways'.

    Split rather than kept whole, so corroboration is queryable — otherwise
    'KEGG WikiPathways' and 'WikiPathways KEGG' look like different sources
    when they record the same agreement.
    """
    return sorted({s for s in (field or "").split() if s})


def load_pairs(path: Path, *, genes_filter=None) -> "list[dict]":
    out = []
    want = {g.upper() for g in genes_filter} if genes_filter else None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            a, b, rel, pathway, src = f[0], f[1], f[2], f[3], f[4]
            if want and a.upper() not in want and b.upper() not in want:
                continue
            out.append({"gene_a": a, "gene_b": b,
                        "relation": rel.split()[0] if rel else "",
                        "relations_all": rel,
                        "pathway": pathway,
                        "sources": _split_sources(src),
                        "n_sources": len(_split_sources(src))})
    return out


def load_memberships(path: Path, *, genes_filter=None) -> "list[dict]":
    out = []
    want = {g.upper() for g in genes_filter} if genes_filter else None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            pathway, gene, src = f[0], f[1], f[2]
            if want and gene.upper() not in want:
                continue
            out.append({"pathway": pathway, "gene": gene,
                        "sources": _split_sources(src),
                        "n_sources": len(_split_sources(src))})
    return out


# ---------------------------------------------------------------------------
# Knowledge-graph ingestion
# ---------------------------------------------------------------------------

def ingest(pairs, memberships, *, label: str = "intpath",
           min_sources: int = 1) -> dict:
    """Write IntPath edges into the local knowledge graph."""
    try:
        from igvfagent import _localstore as ls
    except Exception:
        try:
            import _localstore as ls  # type: ignore
        except Exception:
            return {"error": "local store unavailable"}

    added = collections.Counter()
    con = ls._connect()
    try:
        run = ls.upsert_node(
            con, "analysis", f"intpath_{label}", source="igvfagent:intpath",
            label=f"IntPath ingest {label}",
            properties={"skill": "intpath", "n_pairs": len(pairs),
                        "n_memberships": len(memberships)})

        for m in memberships:
            if m["n_sources"] < min_sources:
                continue
            g = ls.upsert_node(con, "gene", m["gene"], source="IntPath")
            pw = ls.upsert_node(con, "pathway", m["pathway"], source="IntPath")
            ls.upsert_edge(con, g, pw, "member_of_pathway",
                            source="IntPath",
                            properties={"databases": ",".join(m["sources"]),
                                        "n_sources": m["n_sources"]})
            added["member_of_pathway"] += 1

        for p in pairs:
            if p["n_sources"] < min_sources:
                continue
            etype, _desc = RELATION_EDGE.get(p["relation"],
                                              ("interacts_with", ""))
            a = ls.upsert_node(con, "gene", p["gene_a"], source="IntPath")
            b = ls.upsert_node(con, "gene", p["gene_b"], source="IntPath")
            ls.upsert_edge(con, a, b, etype, source="IntPath",
                            properties={"relation": p["relation"],
                                        "pathway": p["pathway"],
                                        "databases": ",".join(p["sources"]),
                                        "n_sources": p["n_sources"]})
            added[etype] += 1

        ls.upsert_edge(con, run, run, "analyzed", source="igvfagent:intpath")
        con.commit()
    finally:
        con.close()
    return dict(added)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _missing(pairs_path, genes_path) -> int:
    print("error: IntPath data files not found.\n"
          f"  Looked in: {DATA_DIR}\n"
          "  Expected: sapiensIntPathGenePairs and sapiensIntPathGenes\n"
          "  Pass --pairs/--genes explicitly, or set IGVF_INTPATH_DIR.\n"
          "  The data is not redistributed with IGVFagent.", file=sys.stderr)
    return 2


def cmd_status(args) -> int:
    pairs_path = _resolve(args.pairs, "sapiensIntPathGenePairs")
    genes_path = _resolve(args.genes_file, "sapiensIntPathGenes")
    print(f"IntPath data dir: {DATA_DIR}")
    print(f"  gene pairs:   {pairs_path or '(not found)'}")
    print(f"  memberships:  {genes_path or '(not found)'}")
    if not (pairs_path and genes_path):
        return 2
    pairs = load_pairs(pairs_path)
    mems = load_memberships(genes_path)
    dbs = collections.Counter(s for p in pairs for s in p["sources"])
    rels = collections.Counter(p["relation"] for p in pairs)
    print(f"\n  {len(pairs):,} gene pairs · "
          f"{len({p['pathway'] for p in pairs}):,} pathways · "
          f"{len({p['gene_a'] for p in pairs} | {p['gene_b'] for p in pairs}):,} genes")
    print(f"  {len(mems):,} memberships · "
          f"{len({m['pathway'] for m in mems}):,} pathways")
    print("\n  source databases:")
    for d, n in dbs.most_common():
        print(f"    {d:16s} {n:>7,}")
    print("\n  relation types:")
    for r, n in rels.most_common():
        _e, desc = RELATION_EDGE.get(r, ("", "unmapped"))
        print(f"    {r:8s} {n:>7,}  {desc}")
    corroborated = sum(1 for p in pairs if p["n_sources"] > 1)
    print(f"\n  {corroborated:,} pair(s) supported by more than one database")
    return 0


def cmd_load(args) -> int:
    pairs_path = _resolve(args.pairs, "sapiensIntPathGenePairs")
    genes_path = _resolve(args.genes_file, "sapiensIntPathGenes")
    if not (pairs_path and genes_path):
        return _missing(pairs_path, genes_path)
    pairs = load_pairs(pairs_path)
    mems = load_memberships(genes_path)
    print(f"Loaded:        {len(pairs):,} pairs, {len(mems):,} memberships")
    if args.dry_run:
        print("(dry run — nothing written)")
        return 0
    added = ingest(pairs, mems, label=args.label, min_sources=args.min_sources)
    print(f"KG:            {json.dumps(added)}")
    return 0 if "error" not in added else 1


def cmd_query(args) -> int:
    pairs_path = _resolve(args.pairs, "sapiensIntPathGenePairs")
    genes_path = _resolve(args.genes_file, "sapiensIntPathGenes")
    if not (pairs_path and genes_path):
        return _missing(pairs_path, genes_path)
    want = [g.strip() for g in args.genes.split(",") if g.strip()]
    pairs = load_pairs(pairs_path, genes_filter=want)
    mems = load_memberships(genes_path, genes_filter=want)

    ts = time.strftime("%Y%m%d_%H%M%S")
    run = OUT_DIR / f"{ts}_{args.label}"
    run.mkdir(parents=True, exist_ok=True)
    import csv as _csv
    with open(run / "edges.csv", "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["source", "target", "kind", "evidence", "score",
                    "relation", "pathway"])
        for p in pairs:
            w.writerow([p["gene_a"], p["gene_b"], "ppi",
                        "IntPath:" + "+".join(p["sources"]),
                        p["n_sources"], p["relation"], p["pathway"]])
        for m in mems:
            w.writerow([m["gene"], m["pathway"], "pathway",
                        "IntPath:" + "+".join(m["sources"]),
                        m["n_sources"], "", m["pathway"]])

    by_pw = collections.Counter(m["pathway"] for m in mems)
    covered = {m["gene"] for m in mems} | {p["gene_a"] for p in pairs} \
        | {p["gene_b"] for p in pairs}
    missing = [g for g in want if g.upper() not in {c.upper() for c in covered}]

    print(f"Genes queried: {len(want)}")
    print(f"Pathways:      {len(by_pw)}")
    print(f"Gene pairs:    {len(pairs)}")
    if missing:
        print(f"Not in IntPath: {', '.join(missing)}  "
              f"(the release covers 5,748 genes; absence here is not evidence "
              f"of absence in biology)")
    print("\nTop pathways:")
    for pw, n in by_pw.most_common(args.top):
        gs = sorted({m["gene"] for m in mems if m["pathway"] == pw})
        print(f"  {n:>2} genes  {pw[:56]:<56} {','.join(gs[:5])}")
    print(f"\nManifest:      {run / 'edges.csv'}")
    print("  (edges.csv uses the network-edge shape, so `localstore harvest` "
          "absorbs it)")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent intpath",
        description="Integrated pathway data (KEGG + WikiPathways + BioCyc).")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, helptext in (("status", "Show what data is available"),
                            ("load", "Ingest the full release into the KG"),
                            ("query", "Pathways and partners for genes")):
        s = sub.add_parser(name, help=helptext)
        s.add_argument("--pairs", help="sapiensIntPathGenePairs path")
        s.add_argument("--genes-file", help="sapiensIntPathGenes path")
        if name == "query":
            s.add_argument("--genes", required=True)
            s.add_argument("--top", type=int, default=12)
            s.add_argument("--label", default="intpath_query")
        if name == "load":
            s.add_argument("--label", default="release")
            s.add_argument("--min-sources", type=int, default=1,
                           help="Keep only facts supported by >= N databases")
            s.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)
    return {"status": cmd_status, "load": cmd_load, "query": cmd_query}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
