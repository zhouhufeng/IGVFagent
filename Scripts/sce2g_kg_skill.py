"""scE2G → local Knowledge Graph bulk ingestion.

Pulls single-cell / bulk enhancer(element)→gene regulatory linkages from the
IGVF **Catalog** endpoint ``/api/genomic-elements/genes`` and loads them as
``regulates`` edges into the local KG (``Data/KG/local_kg.sqlite``), so they
are available for downstream non-coding-variant → gene analysis.

Why this is a *deterministic* skill and not an LLM loop
------------------------------------------------------
The bulk pull is a large, mechanical job (potentially hours). It belongs in
one tool call that runs to completion and streams progress, rather than being
decomposed into hundreds of LLM turns that hit an iteration cap. This is the
IGVFagent pattern for "long jobs that must finish no matter how long".

Completeness under a capped API
-------------------------------
The endpoint hard-caps every response at 500 rows and *ignores* ``skip``, so
ordinary pagination cannot retrieve a window with > 500 links. Completeness is
achieved by **adaptive region tiling**: query a window; if it returns the full
page (>= cap, so possibly truncated) and is wider than ``--min-window``, split
it in half and recurse; otherwise it is complete and every row is ingested. A
window that is still at the cap at ``--min-window`` is logged as a possible
truncation (honest reporting) rather than silently dropped.

Resumable + progress
--------------------
Every completed leaf window is recorded in the KG ``harvest_ledger``; a resumed
run skips finished windows. Node/edge upserts are idempotent, so re-running is
always safe. A heartbeat line is printed for every window and a rollup every
``--heartbeat`` windows, so a multi-hour run reports what it is doing.

License: Apache-2.0.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    from igvfagent import catalog_query_skill as _cat
    from igvfagent import _localstore as _ls
except Exception:  # running from a checkout
    import catalog_query_skill as _cat  # type: ignore
    import _localstore as _ls            # type: ignore

E2G_PATH = "/api/genomic-elements/genes"
PAGE_CAP = 500          # server hard-caps responses at 500 rows
EDGE_TYPE = "regulates"
EDGE_SOURCE = "IGVF-Catalog:scE2G"

# GRCh38 chromosome lengths (bp). The Catalog reports GRCh38 coordinates.
GRCH38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def _fetch_window(chrom: str, start: int, end: int) -> "list[dict]":
    """One capped query for a region. Returns [] on any transient error
    (logged) so a single bad window never aborts the whole job."""
    region = f"{chrom}:{start}-{end}"
    try:
        rows = _cat._catalog_get(
            E2G_PATH,
            {"region": region, "limit": str(PAGE_CAP), "verbose": "true"},
        )
    except SystemExit as exc:            # _catalog_get raises SystemExit on HTTP errs
        print(f"  ! {region}: query failed ({exc}); skipping", flush=True)
        return []
    except Exception as exc:             # noqa: BLE001
        print(f"  ! {region}: {type(exc).__name__}: {exc}; skipping", flush=True)
        return []
    return rows if isinstance(rows, list) else []


def _ingest_rows(con, rows: "list[dict]") -> "tuple[int, int, int]":
    """Upsert element/gene nodes and regulates edges. Returns
    (edges, elements, genes) newly touched."""
    n_edge = n_elem = n_gene = 0
    for r in rows:
        ge = r.get("genomic_element") or {}
        gn = r.get("gene") or {}
        elem_name = ge.get("_id") or ge.get("name")
        gene_name = gn.get("name") or gn.get("_id")
        if not elem_name or not gene_name:
            continue
        elem_id = _ls.upsert_node(
            con, "genomic_element", str(elem_name), source="IGVF-Catalog",
            label=ge.get("name"),
            properties={k: ge.get(k) for k in ("chr", "start", "end", "type")
                        if ge.get(k) is not None})
        n_elem += 1
        gene_id = _ls.upsert_node(
            con, "gene", str(gene_name), source="IGVF-Catalog",
            properties={"ensembl_id": gn.get("_id"), "chr": gn.get("chr"),
                        "start": gn.get("start"), "end": gn.get("end")})
        n_gene += 1
        _ls.upsert_edge(
            con, elem_id, gene_id, EDGE_TYPE, source=EDGE_SOURCE,
            properties={k: r.get(k) for k in
                        ("method", "class", "score", "p_value",
                         "biosample_term", "biological_context",
                         "files_filesets", "source", "source_url")
                        if r.get(k) is not None})
        n_edge += 1
    return n_edge, n_elem, n_gene


def _fmt(n: int) -> str:
    return f"{n:,}"


def _hms(sec: float) -> str:
    sec = int(sec)
    return f"{sec // 3600:02d}:{(sec % 3600) // 60:02d}:{sec % 60:02d}"


def cmd_pull(args: argparse.Namespace) -> int:
    # Resolve the target chromosomes / region.
    if args.region and args.region.lower() != "all":
        if ":" in args.region:
            chrom, span = args.region.split(":", 1)
            s, e = span.replace(",", "").split("-")
            seed = [(chrom, int(s), int(e))]
        else:
            chrom = args.region
            seed = [(chrom, 0, GRCH38.get(chrom, 250_000_000))]
    else:
        chroms = ([f"chr{c}" for c in args.chromosomes.split(",")]
                  if args.chromosomes else list(GRCH38))
        seed = [(c, 0, GRCH38[c]) for c in chroms if c in GRCH38]

    con = _ls._connect()
    resume = not args.no_resume

    def ledger_done(key: str) -> bool:
        return bool(con.execute(
            "SELECT 1 FROM harvest_ledger WHERE key=?", (key,)).fetchone())

    def ledger_mark(key: str):
        con.execute(
            "INSERT OR REPLACE INTO harvest_ledger(key,kind,harvested_at) "
            "VALUES(?,?,?)", (key, "sce2g_window", _ls._NOW()))

    t0 = time.time()
    tot_edges = tot_elem = tot_gene = 0
    windows_done = windows_split = windows_truncated = queries = 0
    heartbeat = max(1, int(args.heartbeat))

    print(f"▶ scE2G → local KG bulk pull  ({E2G_PATH})", flush=True)
    print(f"  targets: {', '.join(sorted({c for c, _, _ in seed}))}"
          f"  · min-window={_fmt(args.min_window)}bp · resume={resume}"
          f"  · db={_ls.KG_PATH}", flush=True)

    # DFS over windows using an explicit stack (chrom, start, end).
    stack = list(reversed(seed))
    while stack:
        if args.max_windows and windows_done >= args.max_windows:
            print(f"  (stopped: --max-windows={args.max_windows} reached)",
                  flush=True)
            break
        chrom, start, end = stack.pop()
        if end <= start:
            continue
        key = f"sce2g:{chrom}:{start}-{end}"
        if resume and ledger_done(key):
            continue

        rows = _fetch_window(chrom, start, end)
        queries += 1
        width = end - start

        if len(rows) >= PAGE_CAP and width > args.min_window:
            # Possibly truncated — subdivide and recurse.
            mid = start + width // 2
            stack.append((chrom, mid, end))
            stack.append((chrom, start, mid))
            windows_split += 1
            if args.verbose:
                print(f"  · {chrom}:{_fmt(start)}-{_fmt(end)} "
                      f"→ {len(rows)} (cap) SPLIT", flush=True)
            continue

        # Leaf window — ingest everything.
        ne, nl, ng = _ingest_rows(con, rows)
        con.commit()
        ledger_mark(key)
        con.commit()
        tot_edges += ne
        tot_elem += nl
        tot_gene += ng
        windows_done += 1
        if len(rows) >= PAGE_CAP:
            windows_truncated += 1
            print(f"  ⚠ {chrom}:{_fmt(start)}-{_fmt(end)} still at cap "
                  f"({PAGE_CAP}) at min-window — possible truncation", flush=True)

        if windows_done % heartbeat == 0:
            el = time.time() - t0
            rate = tot_edges / el if el else 0
            print(f"  · {chrom}:{_fmt(start)}-{_fmt(end)}  "
                  f"windows={_fmt(windows_done)} (+{_fmt(windows_split)} split)  "
                  f"edges={_fmt(tot_edges)}  queries={_fmt(queries)}  "
                  f"{_hms(el)}  {rate:.0f} edges/s", flush=True)

    con.commit()
    n_nodes = con.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    n_edges = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    n_reg = con.execute(
        "SELECT COUNT(*) FROM edges WHERE edge_type=? AND source=?",
        (EDGE_TYPE, EDGE_SOURCE)).fetchone()[0]
    con.close()

    el = time.time() - t0
    # Persist a summary artefact.
    out_dir = Path(_ls.KG_PATH).resolve().parent.parent.parent / "Docs" / "KGTraversal"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "endpoint": E2G_PATH, "edge_type": EDGE_TYPE, "source": EDGE_SOURCE,
        "windows_ingested": windows_done, "windows_split": windows_split,
        "windows_truncated": windows_truncated, "api_queries": queries,
        "edges_touched": tot_edges, "elements_touched": tot_elem,
        "genes_touched": tot_gene, "kg_total_nodes": n_nodes,
        "kg_total_edges": n_edges, "kg_regulates_edges": n_reg,
        "elapsed_sec": round(el, 1),
    }
    sfile = out_dir / "sce2g_kg_ingest_summary.json"
    sfile.write_text(json.dumps(summary, indent=2))

    print(f"\n✓ scE2G ingest complete in {_hms(el)}", flush=True)
    print(f"  regulates edges (this source) in KG: {_fmt(n_reg)}", flush=True)
    print(f"  KG totals: {_fmt(n_nodes)} nodes / {_fmt(n_edges)} edges", flush=True)
    print(f"  windows: {_fmt(windows_done)} ingested, {_fmt(windows_split)} split, "
          f"{_fmt(windows_truncated)} truncated · {_fmt(queries)} API queries", flush=True)
    if windows_truncated:
        print(f"  ⚠ {windows_truncated} window(s) hit the {PAGE_CAP}-row cap at "
              f"--min-window; lower --min-window to recover those.", flush=True)
    print(f"  Summary: {sfile}", flush=True)
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent sce2g-kg",
        description="Bulk-ingest scE2G element→gene linkages from the IGVF "
                    "Catalog into the local Knowledge Graph.")
    sub = p.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("pull", help="Adaptive-tiling bulk pull into the KG.")
    pl.add_argument("--region", default="all",
                    help="'all' (whole genome, default), a chromosome "
                         "('chr19'), or a locus ('chr19:44900000-45000000').")
    pl.add_argument("--chromosomes", default=None,
                    help="Comma-separated chroms to restrict 'all' "
                         "(e.g. '19,20,X'). Prefixed automatically.")
    pl.add_argument("--min-window", type=int, default=20000,
                    help="Smallest window (bp) before a capped page is "
                         "accepted as-is (default 20000).")
    pl.add_argument("--heartbeat", type=int, default=25,
                    help="Print a rollup every N ingested windows (default 25).")
    pl.add_argument("--max-windows", type=int, default=0,
                    help="Stop after N ingested windows (0 = unlimited; for tests).")
    pl.add_argument("--no-resume", action="store_true",
                    help="Ignore the harvest ledger and re-scan all windows.")
    pl.add_argument("--verbose", action="store_true",
                    help="Log every split window.")
    pl.set_defaults(func=cmd_pull)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
