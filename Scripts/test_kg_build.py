"""kg-build turns a Catalog mirror into a graph that keeps every edge accounted for.

Against a small stand-in mirror shaped like the real one (Arango `_id`s,
shards whose column sets differ, human and mouse collections, an edge whose
endpoint is another edge), this checks that:
  - each node collection lands in its table with a CURIE, a taxon and its
    other columns in `props`; human and mouse genes share the Gene table;
  - an edge whose endpoint is missing goes to the rejects with a reason,
    and every input edge is either written or rejected, never dropped;
  - an endpoint that is itself an edge (`variants_proteins/…`) is resolved
    to that edge's variant, with its protein kept in `via`;
  - only variants a core edge touches become Variant nodes, and a touched
    variant the mirror lacks still becomes an id-only node;
  - the measurement-only collections never reach the graph;
  - `validate` flags a reject rate over its limit;
  - the LadybugDB build loads every node and edge and answers a 3-hop query.

    python3 Scripts/test_kg_build.py
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

TMP = Path(tempfile.mkdtemp(prefix="kgbuild_test_"))
os.environ["IGVF_PROJECT_ROOT"] = str(TMP)
os.environ["IGVF_KG_MIRROR_DIR"] = str(TMP / "KG")
os.environ["IGVF_KG_BUILD_DIR"] = str(TMP / "KGBuild")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kg_build_skill as B  # noqa: E402
import duckdb  # noqa: E402

FAILS: list[str] = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"{'PASS' if ok else 'FAIL'} {name}: got {got!r}" + ("" if ok else f", want {want!r}"))
    if not ok:
        FAILS.append(name)


def shard(collection: str, idx: int, rows: list[dict]) -> None:
    d = TMP / "KG" / collection
    d.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    cols = sorted({k for r in rows for k in r})
    vals = ", ".join("(" + ", ".join(
        "NULL" if r.get(c) is None else
        (str(r[c]) if isinstance(r[c], (int, float)) and not isinstance(r[c], bool) else
         ("true" if r[c] is True else "false" if r[c] is False else B._lit(str(r[c]))))
        for c in cols) + ")" for r in rows)
    con.execute(f"COPY (SELECT * FROM (VALUES {vals}) t({', '.join(B._q(c) for c in cols)})) "
                f"TO {B._lit(str(d / f'{idx:04d}.parquet'))}")


def node(coll: str, key: str, **kw) -> dict:
    return {"_id": f"{coll}/{key}", "_key": key, "_rev": "r", **kw}


def edge(coll: str, key: str, frm: str, to: str, **kw) -> dict:
    return {"_id": f"{coll}/{key}", "_key": key, "_rev": "r", "_from": frm, "_to": to, **kw}


def build_mirror() -> None:
    shard("genes", 0, [node("genes", "ENSG1", symbol="PCSK9", name="PCSK9", chr="chr1", start=100, end=200,
                            gene_type="protein_coding", source="GENCODE", hgnc="HGNC:20001")])
    # second shard has a column the first lacks
    shard("genes", 1, [node("genes", "ENSG2", symbol="LDLR", name="LDLR", chr="chr19", start=5, end=9,
                            source="GENCODE", entrez="3949")])
    shard("mm_genes", 0, [node("mm_genes", "ENSMUSG1", symbol="Pcsk9", name="Pcsk9", chr="chr4",
                               start=1, end=2, source="GENCODE")])
    shard("pathways", 0, [node("pathways", "R-HSA-1", name="Lipid metabolism", source="Reactome")])
    shard("ontology_terms", 0, [node("ontology_terms", "MONDO_0005148", name="hypercholesterolemia",
                                     term_id="MONDO:0005148"),
                                node("ontology_terms", "EFO_0004611", name="LDL cholesterol")])
    shard("proteins", 0, [node("proteins", "ENSP1", name="PCSK9_HUMAN", organism="Homo sapiens",
                               protein_id="ENSP1")])
    shard("variants", 0, [node("variants", "NC_000001.11:150:A:G", chr="chr1", pos=150, ref="A", alt="G",
                               rsid='["rs1"]', spdi="NC_000001.11:150:A:G", annotations="{big}"),
                          node("variants", "NC_000001.11:999:C:T", chr="chr1", pos=999, ref="C", alt="T")])
    shard("genes_pathways", 0, [edge("genes_pathways", "1", "genes/ENSG1", "pathways/R-HSA-1", source="Reactome"),
                                edge("genes_pathways", "2", "genes/ENSG2", "pathways/R-HSA-1", source="Reactome"),
                                edge("genes_pathways", "3", "genes/ENSG404", "pathways/R-HSA-1", source="Reactome")])
    shard("variants_genes", 0, [
        edge("variants_genes", "a", "variants/NC_000001.11:150:A:G", "genes/ENSG1",
             neg_log10_pvalue=12.5, effect_size=0.3, source="GTEx"),
        # this variant is referenced but not in the mirror
        edge("variants_genes", "b", "variants/NC_000002.12:7:G:C", "genes/ENSG2",
             neg_log10_pvalue="4.0", source="GTEx"),
    ])
    shard("variants_phenotypes", 0, [edge("variants_phenotypes", "p", "variants/NC_000001.11:150:A:G",
                                          "ontology_terms/EFO_0004611", p_value=1e-9)])
    shard("genes_mm_genes", 0, [edge("genes_mm_genes", "o", "genes/ENSG1", "mm_genes/ENSMUSG1")])
    shard("diseases_genes", 0, [edge("diseases_genes", "d", "ontology_terms/MONDO_0005148", "genes/ENSG2"),
                                edge("diseases_genes", "x", "weird/thing", "genes/ENSG2")])
    # measurement-only collection used to resolve the hyperedge below
    shard("variants_proteins", 0, [edge("variants_proteins", "vp1", "variants/NC_000001.11:150:A:G",
                                        "proteins/ENSP1", score=0.9)])
    shard("variants_proteins_terms", 0, [
        edge("variants_proteins_terms", "t1", "variants_proteins/vp1", "ontology_terms/EFO_0004611"),
        edge("variants_proteins_terms", "t2", "variants_proteins/vp_missing", "ontology_terms/EFO_0004611"),
    ])
    shard("variants_variants", 0, [edge("variants_variants", "ld", "variants/NC_000001.11:150:A:G",
                                        "variants/NC_000001.11:999:C:T", r2=0.8)])
    inv = TMP / "KG" / "_inventory.csv"
    rows = [("doc", "genes", 2), ("doc", "mm_genes", 1), ("doc", "pathways", 1),
            ("doc", "ontology_terms", 2), ("doc", "proteins", 1), ("doc", "variants", 2),
            ("edge", "genes_pathways", 3), ("edge", "variants_genes", 2),
            ("edge", "variants_phenotypes", 1), ("edge", "genes_mm_genes", 1),
            ("edge", "diseases_genes", 2), ("edge", "variants_proteins", 1),
            ("edge", "variants_proteins_terms", 2), ("edge", "variants_variants", 1)]
    inv.write_text("type,collection,documents,bytes,skip_default\n" +
                   "".join(f"{t},{c},{n},0,False\n" for t, c, n in rows))


def main() -> int:
    build_mirror()
    common = ["--snapshot", "t", "--threads", "2", "--memory", "2GB"]
    check("nodes exit", B.main(common + ["nodes"]), 0)
    check("edges exit", B.main(common + ["edges"]), 0)
    check("variants exit", B.main(common + ["variants"]), 0)

    sd = TMP / "KGBuild" / "t"
    m = json.loads((sd / "manifest.json").read_text())
    con = duckdb.connect()

    def rows(sql: str):
        return con.execute(sql).fetchall()

    genes = rows(f"SELECT id, curie, name, taxon, pos_start, props FROM read_parquet('{sd}/nodes/Gene/*.parquet') ORDER BY id")
    check("Gene rows (human + mouse)", [g[0] for g in genes], ["genes/ENSG1", "genes/ENSG2", "mm_genes/ENSMUSG1"])
    check("Gene curie", genes[0][1], "ENSEMBL:ENSG1")
    check("Gene taxon mouse", genes[2][3], B.MOUSE)
    check("Gene coordinates typed", genes[0][4], 100)
    check("Gene props keep other columns", json.loads(genes[1][5]).get("entrez"), "3949")
    ont = rows(f"SELECT curie, category FROM read_parquet('{sd}/nodes/OntologyTerm/*.parquet') ORDER BY curie")
    check("ontology curie + category", ont, [("EFO:0004611", "PhenotypicFeature"), ("MONDO:0005148", "Disease")])

    gp = m["edges"]["genes_pathways"]
    check("dangling gene rejected", (gp["rows_out"], gp["rejects"]), (2, 1))
    rej = rows(f"SELECT raw_from, reason FROM read_parquet('{sd}/rejects/genes_pathways.parquet')")
    check("reject reason", rej, [("genes/ENSG404", "dangling endpoint")])
    dg = m["edges"]["diseases_genes"]
    check("unknown endpoint rejected", (dg["rows_out"], dg["rejects"]), (1, 1))
    for coll, e in m["edges"].items():
        check(f"{coll}: in = out + rejects", e["rows_in"], e["rows_out"] + e["rejects"])

    vg = rows(f"SELECT from_id, neg_log10_pvalue FROM read_parquet('{sd}/edges/variants_genes/Variant__Gene.parquet') ORDER BY 1")
    check("string p-value cast to DOUBLE", vg[1][1], 4.0)
    vpt = rows(f"SELECT from_id, to_id, via FROM read_parquet('{sd}/edges/variants_proteins_terms/*.parquet')")
    check("hyperedge resolved to variant, protein in via", vpt,
          [("variants/NC_000001.11:150:A:G", "ontology_terms/EFO_0004611", "proteins/ENSP1")])
    check("unresolvable hyperedge rejected", m["edges"]["variants_proteins_terms"]["rejects"], 1)
    check("measurement-only not built", sorted(p.name for p in (sd / "edges").iterdir()) ,
          sorted(["genes_pathways", "variants_genes", "variants_phenotypes", "genes_mm_genes",
                  "diseases_genes", "variants_proteins_terms"]))

    var = rows(f"SELECT id, curie, chr, pos, props FROM read_parquet('{sd}/nodes/Variant/*.parquet') ORDER BY id")
    check("only touched variants", [v[0] for v in var],
          ["variants/NC_000001.11:150:A:G", "variants/NC_000002.12:7:G:C"])
    check("variant curie", var[0][1], "SPDI:NC_000001.11:150:A:G")
    check("annotations blob left out of props", "annotations" in json.loads(var[0][4] or "{}"), False)
    check("missing variant kept id-only", json.loads(var[1][4]), {"_missing": "not in mirror"})
    check("variant counts", (m["variants"]["referenced"], m["variants"]["found"], m["variants"]["missing"]), (2, 1, 1))

    check("validate fails on 1/3 reject rate", B.main(common + ["validate"]), 1)
    check("validate passes with a loose limit", B.main(common + ["validate", "--max-reject-rate", "0.6"]), 0)

    check("graph exit", B.main(common + ["graph", "--buffer-gb", "0.5"]), 0)
    import ladybug as lb
    db = lb.Database(str(sd / "graph.lbug"), read_only=True, buffer_pool_size=512 * 1024 ** 2)
    g = lb.Connection(db)
    check("graph node count", g.execute("MATCH (n) RETURN count(*)").get_all()[0][0], 9)
    check("graph edge count", g.execute("MATCH ()-[r]->() RETURN count(*)").get_all()[0][0], 8)
    check("3-hop variant -> gene -> pathway",
          g.execute("MATCH (v:Variant)-[:variants_genes]->(g:Gene)-[:genes_pathways]->(p:Pathway) "
                    "RETURN v.id, g.symbol, p.name ORDER BY v.id").get_all(),
          [["variants/NC_000001.11:150:A:G", "PCSK9", "Lipid metabolism"],
           ["variants/NC_000002.12:7:G:C", "LDLR", "Lipid metabolism"]])
    check("ortholog edge human -> mouse",
          g.execute("MATCH (a:Gene)-[:genes_mm_genes]->(b:Gene) RETURN a.symbol, b.taxon").get_all(),
          [["PCSK9", B.MOUSE]])
    check("rebuild refused without --rebuild", _raises(lambda: B.main(common + ["graph"])), True)

    print(f"\n{'ALL PASS' if not FAILS else f'{len(FAILS)} FAILED: ' + ', '.join(FAILS)}")
    return 1 if FAILS else 0


def _raises(fn) -> bool:
    try:
        fn()
    except SystemExit:
        return True
    return False


if __name__ == "__main__":
    sys.exit(main())
