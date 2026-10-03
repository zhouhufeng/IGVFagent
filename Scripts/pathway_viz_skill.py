"""Pathway and interaction networks for a gene list.

Asked "visualise the pathway networks of these genes", the agent had no route:
`enrich_pathways` needs a gene list but depends on `gseapy`, which is not
installed anywhere; nothing drew a network from genes alone; and the agent
worked around the gap by authoring a `write_text_file` tool pointing at an
`igvfagent write-text-file` subcommand that does not exist — so it failed
three times and exhausted its budget without producing a figure.

This skill takes the gene list **inline** and produces the figure directly:

    igvfagent pathway-viz network --genes CLU,BIN1,PICALM,SORL1,ABCA7

Edges come from three sources, each optional and each labelled on the figure
so a viewer can tell evidence apart:

* **STRING** — scored protein-protein interactions (public API, no key)
* **Reactome** — shared pathway membership, drawn as gene↔pathway edges
* **the local knowledge graph** — whatever prior IGVFagent runs accumulated

Rendering uses networkx + matplotlib, which are already required by the
analysis extra, so this adds no new hard dependency. `gseapy` is not needed;
where it would have been used for over-representation statistics, a
hypergeometric test is computed directly.

Nothing is invented. A gene with no edges is drawn as an isolated node rather
than dropped, because "not connected in these sources" and "not in the figure"
mean different things to a reader.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

__all__ = ["main", "fetch_string_edges", "fetch_reactome_pathways",
           "build_graph"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
OUT_DIR = ROOT / "Docs" / "PathwayViz"

STRING_API = "https://string-db.org/api"
REACTOME_API = "https://reactome.org/ContentService"
_UA = "igvfagent (https://github.com/zhouhufeng/IGVFagent)"


def _get(url: str, *, timeout: int = 60) -> "str | None":
    try:
        req = urllib.request.Request(url, headers={"User-Agent": _UA})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace")
    except Exception:
        return None


def parse_genes(spec: str) -> "list[str]":
    """Accept a comma/space/newline list, or a path to a file of genes."""
    if not spec:
        return []
    p = Path(spec).expanduser()
    if not p.is_absolute():
        p = ROOT / spec
    if p.is_file():
        spec = p.read_text(encoding="utf-8", errors="replace")
    out, seen = [], set()
    for tok in re.split(r"[\s,;|]+", spec.strip()):
        tok = tok.strip().strip('"\'')
        if tok and tok.upper() not in seen and not tok.startswith("#"):
            seen.add(tok.upper())
            out.append(tok)
    return out


# ---------------------------------------------------------------------------
# Edge sources
# ---------------------------------------------------------------------------

def fetch_string_edges(genes, *, species: int = 9606,
                       min_score: float = 0.4) -> "list[dict]":
    """Scored protein-protein interactions from STRING.

    `min_score` defaults to STRING's own "medium confidence" cutoff. Drawing
    every edge STRING will return produces a hairball in which a 0.15 text-
    mining coincidence looks identical to an experimentally supported complex.
    """
    if len(genes) < 2:
        return []
    q = urllib.parse.urlencode({
        "identifiers": "\r".join(genes), "species": species,
        "caller_identity": "igvfagent"})
    body = _get(f"{STRING_API}/tsv/network?{q}")
    if not body:
        return []
    rows = list(csv.DictReader(body.splitlines(), delimiter="\t"))
    out = []
    upper = {g.upper() for g in genes}
    for r in rows:
        a = (r.get("preferredName_A") or "").strip()
        b = (r.get("preferredName_B") or "").strip()
        try:
            score = float(r.get("score") or 0)
        except ValueError:
            continue
        # STRING expands to neighbours; keep only edges within the query set,
        # otherwise the figure answers a different question than was asked.
        if not (a and b) or a.upper() not in upper or b.upper() not in upper:
            continue
        if score < min_score:
            continue
        out.append({"source": a, "target": b, "score": round(score, 3),
                    "kind": "ppi", "evidence": "STRING",
                    "experimental": float(r.get("escore") or 0),
                    "database": float(r.get("dscore") or 0),
                    "textmining": float(r.get("tscore") or 0)})
    return out


def fetch_reactome_pathways(genes, *, max_per_gene: int = 6) -> "list[dict]":
    """Reactome pathway membership, as gene→pathway edges."""
    out = []
    for g in genes:
        body = _get(f"{REACTOME_API}/data/mapping/UniProt/"
                    f"{urllib.parse.quote(g)}/pathways?species=9606")
        if not body:
            body = _get(f"{REACTOME_API}/search/query?"
                        + urllib.parse.urlencode(
                            {"query": g, "species": "Homo sapiens",
                             "types": "Pathway", "cluster": "true"}))
            if not body:
                continue
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                continue
            names = []
            for grp in (data.get("results") or []):
                for e in (grp.get("entries") or [])[:max_per_gene]:
                    nm = re.sub(r"<[^>]+>", "", e.get("name") or "").strip()
                    if nm:
                        names.append(nm)
            for nm in names[:max_per_gene]:
                out.append({"source": g, "target": nm, "kind": "pathway",
                            "evidence": "Reactome", "score": 1.0})
            continue
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            continue
        for p in (data if isinstance(data, list) else [])[:max_per_gene]:
            nm = (p.get("displayName") or "").strip()
            if nm:
                out.append({"source": g, "target": nm, "kind": "pathway",
                            "evidence": "Reactome", "score": 1.0,
                            "stId": p.get("stId")})
    return out


KEGG_API = "https://rest.kegg.jp"


def _kegg_gene_ids(genes) -> "dict[str, str]":
    """Symbol -> KEGG hsa: id.

    Resolved from KEGG's own human gene list rather than by guessing an
    Entrez id, and matched on the symbol field only. A substring match on the
    description would map CLU onto "iron-sulfur cluster" — KEGG's find
    endpoint returns exactly that for CLU, so the naive approach silently
    produces the wrong gene.
    """
    body = _get(f"{KEGG_API}/list/hsa", timeout=90)
    if not body:
        return {}
    want = {g.upper(): g for g in genes}
    out = {}
    for line in body.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        kid, names = parts[0], parts[3]
        symbols = [n.strip().upper() for n in names.split(";")[0].split(",")]
        for sym in symbols:
            if sym in want and want[sym] not in out:
                out[want[sym]] = kid
    return out


def fetch_kegg_pathways(genes, *, max_per_gene: int = 6) -> "tuple[list[dict], dict]":
    """KEGG pathway membership as gene->pathway edges, plus pathway ids.

    Returns (edges, {pathway_name: hsa_id}) — the ids let the caller download
    KEGG's own rendered pathway map, which is the diagram most biologists
    picture when they say "KEGG pathway".
    """
    ids = _kegg_gene_ids(genes)
    if not ids:
        return [], {}
    names_body = _get(f"{KEGG_API}/list/pathway/hsa", timeout=60) or ""
    pw_name = {}
    for line in names_body.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pw_name[parts[0].replace("path:", "")] = \
                parts[1].split(" - ")[0].strip()

    edges, pathway_ids = [], {}
    for gene, kid in ids.items():
        body = _get(f"{KEGG_API}/link/pathway/{kid}", timeout=60)
        if not body:
            continue
        hits = []
        for line in body.splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1].startswith("path:"):
                pid = parts[1].replace("path:", "")
                hits.append(pid)
        for pid in hits[:max_per_gene]:
            label = pw_name.get(pid, pid)
            edges.append({"source": gene, "target": label, "kind": "pathway",
                          "evidence": "KEGG", "score": 1.0, "kegg_id": pid})
            pathway_ids[label] = pid
    return edges, pathway_ids


def download_kegg_maps(pathway_ids: dict, out_dir: Path, *, limit: int = 5) -> "list[Path]":
    """KEGG's own rendered pathway diagrams for the top pathways."""
    saved = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for label, pid in list(pathway_ids.items())[:limit]:
        try:
            req = urllib.request.Request(f"{KEGG_API}/get/{pid}/image",
                                         headers={"User-Agent": _UA})
            with urllib.request.urlopen(req, timeout=90) as r:
                data = r.read()
            if not data.startswith(b"\x89PNG"):
                continue
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)[:48]
            path = out_dir / f"kegg_{pid}_{safe}.png"
            path.write_bytes(data)
            saved.append(path)
        except Exception:
            continue
    return saved


def fetch_local_kg_edges(genes) -> "list[dict]":
    """Edges already accumulated in the local knowledge graph."""
    try:
        from igvfagent import _localstore as ls
    except Exception:
        try:
            import _localstore as ls  # type: ignore
        except Exception:
            return []
    out = []
    upper = {g.upper() for g in genes}
    try:
        con = ls._connect()
        rows = con.execute(
            "SELECT from_node, to_node, edge_type, source FROM edges "
            "WHERE edge_type IN ('interacts_with','member_of_pathway',"
            "'regulates','associated_with')").fetchall()
        con.close()
    except Exception:
        return []
    for frm, to, etype, src in rows:
        a = str(frm).split(":", 1)[-1]
        b = str(to).split(":", 1)[-1]
        if a.upper() in upper or b.upper() in upper:
            out.append({"source": a, "target": b, "kind": etype,
                        "evidence": f"local KG ({src})", "score": 0.9})
    return out


# ---------------------------------------------------------------------------
# Enrichment without gseapy
# ---------------------------------------------------------------------------

def _log_choose(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def hypergeometric_p(k, K, n, N) -> float:
    """P(X >= k): k hits in n draws, K successes in a population of N.

    Computed directly rather than via gseapy, which is not installed anywhere
    and was the reason the enrichment path failed. In log space, because the
    binomials overflow at genome scale.
    """
    if k <= 0 or K <= 0 or n <= 0 or N <= 0 or k > min(K, n):
        return 1.0
    total = 0.0
    for i in range(k, min(K, n) + 1):
        lp = _log_choose(K, i) + _log_choose(N - K, n - i) - _log_choose(N, n)
        total += math.exp(lp)
    return min(1.0, max(0.0, total))


def pathway_enrichment(pathway_edges, n_genes: int,
                       *, genome_size: int = 20000,
                       typical_pathway: int = 60) -> "list[dict]":
    """Over-representation per pathway, by hypergeometric test.

    The pathway sizes here come from Reactome's own membership only for the
    queried genes, so the population parameters are approximations and are
    labelled as such — an exact test needs full pathway sizes, which this
    endpoint does not return.
    """
    counts: "dict[str, list]" = {}
    for e in pathway_edges:
        counts.setdefault(e["target"], []).append(e["source"])
    out = []
    for pw, members in counts.items():
        k = len(set(members))
        p = hypergeometric_p(k, typical_pathway, n_genes, genome_size)
        out.append({"pathway": pw, "n_genes": k,
                    "genes": sorted(set(members)), "p_approx": p})
    out.sort(key=lambda r: (r["p_approx"], -r["n_genes"]))
    return out


def ingest_to_kg(genes, edges, *, label: str) -> dict:
    """Write the network back into the local knowledge graph.

    Retrieval that does not accumulate is retrieval done twice. These edges —
    STRING interactions, KEGG and Reactome membership — are exactly the
    biological structure the graph exists to hold, so a pathway search
    performed once should answer the next question locally.

    Sources are preserved on every edge, so a later reader can tell a STRING
    text-mining association from a curated KEGG membership rather than
    inheriting an undifferentiated blob.
    """
    try:
        from igvfagent import _localstore as ls
    except Exception:
        try:
            import _localstore as ls  # type: ignore
        except Exception:
            return {"error": "local store unavailable"}

    added = {"gene_nodes": 0, "pathway_nodes": 0,
             "interacts_with": 0, "member_of_pathway": 0}
    try:
        con = ls._connect()
        run = ls.upsert_node(
            con, "analysis", f"pathway_viz_{label}",
            source="igvfagent:pathway-viz",
            label=f"pathway network {label}",
            properties={"skill": "pathway-viz", "n_genes": len(genes)})

        gene_ids = {}
        for g in genes:
            gid = ls.upsert_node(con, "gene", g, source="igvfagent:pathway-viz")
            gene_ids[g] = gid
            ls.upsert_edge(con, run, gid, "analyzed",
                            source="igvfagent:pathway-viz")
        added["gene_nodes"] = len(gene_ids)

        seen_pw = set()
        for e in edges:
            src, tgt = e.get("source"), e.get("target")
            if not (src and tgt):
                continue
            props = {k: v for k, v in e.items()
                     if k in ("score", "evidence", "kegg_id", "stId",
                              "experimental", "database", "textmining")}
            if e.get("kind") == "ppi":
                a = gene_ids.get(src) or ls.upsert_node(
                    con, "gene", src, source="igvfagent:pathway-viz")
                b = gene_ids.get(tgt) or ls.upsert_node(
                    con, "gene", tgt, source="igvfagent:pathway-viz")
                # Source is the EVIDENCE database, matching what the
                # artefact recogniser records when harvest re-reads
                # edges.csv. Using a skill-specific label instead would
                # produce a second, near-identical edge for the same fact,
                # because the edge id hashes the source.
                ls.upsert_edge(con, a, b, "interacts_with",
                                source=e.get("evidence", "STRING"),
                                properties=props)
                added["interacts_with"] += 1
            elif e.get("kind") == "pathway":
                g = gene_ids.get(src) or ls.upsert_node(
                    con, "gene", src, source="igvfagent:pathway-viz")
                pw = ls.upsert_node(con, "pathway", tgt,
                                    source=e.get("evidence", "pathway"),
                                    properties=props)
                if tgt not in seen_pw:
                    seen_pw.add(tgt)
                ls.upsert_edge(con, g, pw, "member_of_pathway",
                                source=e.get("evidence", "pathway"),
                                properties=props)
                added["member_of_pathway"] += 1
        added["pathway_nodes"] = len(seen_pw)
        con.commit()
        con.close()
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    return added


# ---------------------------------------------------------------------------
# Graph + figure
# ---------------------------------------------------------------------------

def build_graph(genes, edges):
    import networkx as nx
    G = nx.Graph()
    for g in genes:
        G.add_node(g, kind="gene")
    for e in edges:
        for node, default in ((e["source"], "gene"), (e["target"], None)):
            if node not in G:
                kind = default or ("pathway" if e["kind"] == "pathway" else "gene")
                G.add_node(node, kind=kind)
        G.add_edge(e["source"], e["target"], **{
            k: v for k, v in e.items() if k not in ("source", "target")})
    return G


def draw(G, out_png: Path, *, title: str, layout: str = "spring"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import networkx as nx

    if G.number_of_nodes() == 0:
        return None
    pos = {"spring": lambda g: nx.spring_layout(g, seed=0, k=0.9),
           "circular": nx.circular_layout,
           "kamada": nx.kamada_kawai_layout}.get(layout,
            lambda g: nx.spring_layout(g, seed=0, k=0.9))(G)

    genes = [n for n, d in G.nodes(data=True) if d.get("kind") == "gene"]
    pathways = [n for n, d in G.nodes(data=True) if d.get("kind") != "gene"]

    fig, ax = plt.subplots(figsize=(13, 10))
    ppi = [(u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "ppi"]
    pw = [(u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "pathway"]
    other = [(u, v) for u, v, d in G.edges(data=True)
             if d.get("kind") not in ("ppi", "pathway")]

    if ppi:
        widths = [0.8 + 3.2 * G[u][v].get("score", 0.5) for u, v in ppi]
        nx.draw_networkx_edges(G, pos, edgelist=ppi, width=widths,
                               edge_color="#1B5EA8", alpha=0.75, ax=ax)
    if pw:
        nx.draw_networkx_edges(G, pos, edgelist=pw, width=1.0,
                               edge_color="#0F848C", alpha=0.35,
                               style="dashed", ax=ax)
    if other:
        nx.draw_networkx_edges(G, pos, edgelist=other, width=1.0,
                               edge_color="#C26B1F", alpha=0.5, ax=ax)

    nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color="#EFF4F9",
                           edgecolors="#0B2E59", linewidths=2.0,
                           node_size=1500, ax=ax)
    if pathways:
        nx.draw_networkx_nodes(G, pos, nodelist=pathways, node_color="#E6F4F1",
                               edgecolors="#0F848C", linewidths=1.2,
                               node_shape="s", node_size=900, ax=ax)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in genes},
                            font_size=10, font_weight="bold", ax=ax)
    if pathways:
        nx.draw_networkx_labels(
            G, pos, labels={n: (n[:28] + "…" if len(n) > 28 else n)
                            for n in pathways},
            font_size=7, font_color="#0F848C", ax=ax)

    handles = [mpatches.Patch(color="#1B5EA8", label="STRING PPI (width = score)"),
               mpatches.Patch(color="#0F848C", label="Reactome pathway membership"),
               mpatches.Patch(color="#C26B1F", label="local knowledge graph")]
    ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=False)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170)
    fig.savefig(out_png.with_suffix(".svg"))
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

def cmd_network(args) -> int:
    genes = parse_genes(args.genes)
    if not genes:
        print("error: no genes parsed from --genes", file=sys.stderr)
        return 2

    sources = {s.strip().lower() for s in (args.sources or "string,reactome,kegg,local").split(",")}
    edges: "list[dict]" = []
    counts = {}
    if "string" in sources:
        e = fetch_string_edges(genes, min_score=args.min_score)
        counts["STRING"] = len(e); edges += e
    if "reactome" in sources:
        e = fetch_reactome_pathways(genes, max_per_gene=args.max_pathways)
        counts["Reactome"] = len(e); edges += e
    kegg_ids: dict = {}
    if "kegg" in sources:
        e, kegg_ids = fetch_kegg_pathways(genes, max_per_gene=args.max_pathways)
        counts["KEGG"] = len(e); edges += e
    if "local" in sources:
        e = fetch_local_kg_edges(genes)
        counts["local KG"] = len(e); edges += e

    G = build_graph(genes, edges)
    ts = time.strftime("%Y%m%d_%H%M%S")
    label = re.sub(r"[^A-Za-z0-9_.-]", "_", args.label or "pathway_network")
    run = OUT_DIR / f"{ts}_{label}"
    run.mkdir(parents=True, exist_ok=True)

    with open(run / "edges.csv", "w", newline="", encoding="utf-8") as fh:
        cols = ["source", "target", "kind", "evidence", "score"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(edges)

    pw_edges = [e for e in edges if e["kind"] == "pathway"]
    enrich = pathway_enrichment(pw_edges, len(genes)) if pw_edges else []
    (run / "enrichment.json").write_text(json.dumps(enrich, indent=2))

    png = draw(G, run / "network.png",
               title=args.title or f"Pathway / interaction network — {len(genes)} genes",
               layout=args.layout)

    # Absorb into the growing knowledge graph before drawing: the figure is
    # one view of this evidence, the graph is where it persists.
    kg = {"skipped": True}
    if not args.no_kg:
        kg = ingest_to_kg(genes, edges, label=label)

    kegg_maps = []
    if kegg_ids and not args.no_kegg_maps:
        # KEGG's own rendered diagram is what most biologists mean by "the
        # KEGG pathway"; reproducing it from edges would be a worse picture.
        top = {r["pathway"]: kegg_ids[r["pathway"]]
               for r in enrich if r["pathway"] in kegg_ids}
        kegg_maps = download_kegg_maps(top or kegg_ids, run / "kegg_maps",
                                       limit=args.kegg_map_limit)

    isolated = [g for g in genes if G.degree(g) == 0]
    (run / "report.md").write_text(_render(genes, counts, G, enrich, isolated))

    print(f"Genes:         {len(genes)}")
    for k, v in counts.items():
        print(f"  {k:12s} {v} edge(s)")
    print(f"Nodes/edges:   {G.number_of_nodes()} / {G.number_of_edges()}")
    if not kg.get("skipped"):
        print(f"KG:            {json.dumps(kg)}")
    if isolated:
        print(f"Unconnected:   {', '.join(isolated)}  "
              f"(no edge in the selected sources — shown as isolated nodes)")
    if enrich:
        print("Top pathways:")
        for r in enrich[:5]:
            print(f"  {r['n_genes']:>2} genes  p≈{r['p_approx']:.2e}  {r['pathway'][:58]}")
    if png:
        print(f"Report:        {run / 'report.md'}")
        print(f"Wrote:         {png}")
        print(f"Wrote:         {png.with_suffix('.svg')}")
    for m in kegg_maps:
        print(f"Wrote:         {m}")
    print(f"Manifest:      {run / 'edges.csv'}")
    return 0


def _render(genes, counts, G, enrich, isolated) -> str:
    L = [f"# Pathway / interaction network — {len(genes)} genes", "",
         f"Genes: `{', '.join(genes)}`", "",
         "## Edge sources", "", "| source | edges |", "|---|---|"]
    L += [f"| {k} | {v} |" for k, v in counts.items()]
    L += ["", f"Graph: **{G.number_of_nodes()}** nodes, "
              f"**{G.number_of_edges()}** edges.", ""]
    if isolated:
        L += [f"**Unconnected genes:** {', '.join(isolated)} — no edge in the "
              f"selected sources. They are drawn as isolated nodes rather than "
              f"omitted, because absence of an edge and absence from the "
              f"figure are different statements.", ""]
    if enrich:
        L += ["## Pathway over-representation", "",
              "| pathway | genes | p (approx) |", "|---|---|---|"]
        for r in enrich[:15]:
            L.append(f"| {r['pathway'][:60]} | {', '.join(r['genes'][:6])} "
                     f"| {r['p_approx']:.2e} |")
        L += ["", "*p is a hypergeometric approximation: Reactome's mapping "
              "endpoint returns membership for the queried genes but not full "
              "pathway sizes, so the population parameters are estimated. "
              "Treat the ranking as indicative and the values as approximate.*",
              ""]
    L += ["## Caveats", "",
          "- STRING edges are filtered to the queried gene set; neighbours "
          "STRING would add are excluded so the figure answers the question "
          "asked.",
          "- Edge width encodes STRING confidence. A thin edge may rest on "
          "text-mining alone; `edges.csv` carries the per-channel subscores.",
          "- Pathway membership is drawn as gene↔pathway edges, not as "
          "gene↔gene: two genes in one pathway are not thereby known to "
          "interact.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent pathway-viz",
        description="Draw a pathway / interaction network for a gene list.")
    sub = p.add_subparsers(dest="cmd", required=True)
    n = sub.add_parser("network", help="Build and draw the network")
    n.add_argument("--genes", required=True,
                   help="Gene symbols (comma/space separated) or a file path")
    n.add_argument("--sources", default="string,reactome,kegg,local",
                   help="Comma list: string, reactome, kegg, local")
    n.add_argument("--min-score", type=float, default=0.4,
                   help="STRING confidence cutoff (default 0.4 = medium)")
    n.add_argument("--max-pathways", type=int, default=6,
                   help="Reactome pathways per gene")
    n.add_argument("--layout", default="spring",
                   choices=("spring", "circular", "kamada"))
    n.add_argument("--title", default="")
    n.add_argument("--label", default="pathway_network")
    n.add_argument("--no-kegg-maps", action="store_true",
                   help="Skip downloading KEGG's rendered pathway diagrams")
    n.add_argument("--kegg-map-limit", type=int, default=4)
    n.add_argument("--no-kg", action="store_true",
                   help="Draw the figure without writing to the knowledge graph")
    args = p.parse_args(argv)
    return cmd_network(args)


if __name__ == "__main__":
    raise SystemExit(main())
