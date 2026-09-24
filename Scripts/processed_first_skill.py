#!/usr/bin/env python3
"""Find what IGVF has already computed, before computing it again.

A MeasurementSet carries raw reads. It also carries ``input_for`` — the
AnalysisSets built from it — and each of those carries
``uniform_pipeline_status``. When that says ``completed``, the IGVF uniform
pipeline has already aligned, counted and assembled the matrices, and the
outputs are sitting on the Portal.

Nothing in this codebase read either field. ``raw-pipeline plan`` inspected a
MeasurementSet's own files, saw FASTQs, and planned an alignment — for data
IGVF had already aligned. On IGVFDS9875NBZW that is 16 controlled FASTQs, about
70 GB and hours of compute, to reproduce a 3.7 GB h5ad that already exists and
downloads in minutes:

    IGVFDS9875NBZW  ->  input_for  ->  IGVFDS5497CMDG  (uniform, completed)
        IGVFFI8531JIMM   cell by gene matrix   h5ad   3.7 GB
        IGVFFI7525HRSN   fragments             bed    5.6 GB   (public)
        IGVFFI9390KPDV   alignments            bam     38 GB   (controlled)
        IGVFFI1093ZZEP   kallisto matrix       tar     13 GB

Access matters as much as size. Raw reads on a controlled set need credentials
and an approved data-use agreement; the fragments file above is not controlled
and anyone can fetch it. So "already processed" is sometimes the difference
between an analysis that runs and one that cannot run at all.

    igvfagent processed find IGVFDS9875NBZW
    igvfagent processed find IGVFDS9875NBZW --json
    igvfagent processed plan IGVFDS9875NBZW --want matrix

``plan`` answers the question that actually matters at the start of a
workflow: given what I want to do, should I download a processed file or run
the pipeline?

``lineage`` walks the whole Portal graph around an accession (portal_lineage.py):
analysis sets through any number of ``input_for`` hops (intermediate ->
principal -> predictions), the multiome partner, auxiliary sets, the
multiplexed sample and its barcode map, seqspecs and onlists, documents, the
QC metric objects the uniform pipeline published, model sets and prediction
sets, and every linked object the credentials cannot see. It writes a report
with a "start here" table, the Portal's own QC numbers, a lineage figure and
the graph as JSON.

    igvfagent processed lineage IGVFDS9875NBZW
    igvfagent processed fetch IGVFDS9875NBZW --want rna_matrix,fragments --max-gb 20
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

try:
    from igvfagent.raw_data_pipeline import portal_json      # type: ignore
except Exception:                                            # direct execution
    from raw_data_pipeline import portal_json                # type: ignore

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "Processed"
log = logging.getLogger("processed")

# What a downstream analysis actually wants, best first. An h5ad loads
# natively into the single-cell path; a 13 GB tar of kallisto output is the
# same information behind an unpacking step, and a BAM is upstream of both.
MATRIX_RANK = (
    ("cell by gene matrix", "h5ad"),
    ("sparse gene count matrix", "h5ad"),
    ("cell by gene matrix", "tar"),
    ("kallisto cell by gene matrix", "tar"),
)
WANTS = {
    "matrix":    ("cell by gene matrix", "sparse gene count matrix",
                  "kallisto cell by gene matrix"),
    "fragments": ("fragments",),
    "alignments": ("alignments",),
    "peaks":     ("peaks", "pseudoreplicated peaks"),
    "index":     ("index",),
}


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")


def _get(path: str) -> "dict | None":
    status, data = portal_json(path)
    if status != 200 or not isinstance(data, dict):
        return None
    return data


def analysis_sets_for(accession: str) -> "list[dict]":
    """AnalysisSets built from this accession, via ``input_for``.

    Falls back to a search when the object does not carry ``input_for`` --
    older records and some types do not -- so a caller gets the same answer
    whichever shape the Portal returns.
    """
    obj = _get(f"/{accession}/?format=json") or _get(
        f"/search/?type=FileSet&accession={urllib.parse.quote(accession)}"
        f"&format=json")
    if not obj:
        return []
    if obj.get("@graph"):
        obj = obj["@graph"][0]

    out = list(obj.get("input_for") or [])
    if not out:
        data = _get(f"/search/?type=AnalysisSet"
                    f"&input_file_sets.accession={urllib.parse.quote(accession)}"
                    f"&format=json") or {}
        out = data.get("@graph") or []
    return out


def processed_files(accession: str) -> dict:
    """Every processed output reachable from this accession.

    Backed by the full graph walk, so analysis sets several ``input_for`` hops
    away (principal analyses, prediction sets) are included, not only the
    first hop. Falls back to the one-hop lookup if the walk fails.
    """
    try:
        import portal_lineage as pl  # noqa: E402
        g = pl.walk(accession, fetch=portal_json, max_depth=4, fetch_qc=False)
        sets = []
        for k, n in g.nodes.items():
            if n.get("kind") == "fileset" and n.get("type") in ("AnalysisSet", "PredictionSet", "ModelSet"):
                sets.append({"accession": n["accession"], "uniform_pipeline_status": n.get("uniform_pipeline_status") or "",
                             "status": n.get("status") or "", "aliases": n.get("aliases") or [],
                             "workflows": [w.get("accession") or w.get("name") for w in n.get("workflows") or []],
                             "summary": (n.get("summary") or "")[:160], "type": n.get("type"),
                             "file_set_type": n.get("file_set_type")})
        by_set = {k: n for k, n in g.nodes.items() if n.get("kind") == "fileset"}
        files = []
        for k, f in g.nodes.items():
            if f.get("kind") != "file" or f.get("product") in ("raw_reads", "seqspec"):
                continue
            fs = by_set.get(f.get("file_set") or "")
            if not fs or fs.get("type") not in ("AnalysisSet", "PredictionSet", "ModelSet"):
                continue
            files.append({"accession": f.get("accession"), "content_type": f.get("content_type") or "",
                          "file_format": f.get("file_format") or "", "file_size": f.get("file_size"),
                          "controlled_access": f.get("controlled_access"), "href": f.get("href") or "",
                          "s3_uri": f.get("s3_uri") or "", "from_analysis_set": fs["accession"],
                          "uniform": fs.get("uniform_pipeline_status") == "completed"})
        if sets or files:
            return {"accession": accession, "analysis_sets": sets, "files": files,
                    "blocked": [{"@id": k, "http_status": v} for k, v in g.blocked.items()]}
    except Exception as exc:  # pragma: no cover
        log.debug("graph walk failed, one-hop fallback: %s", exc)
    return _processed_files_one_hop(accession)


def _processed_files_one_hop(accession: str) -> dict:
    """Every processed output reachable from this accession (one input_for hop)."""
    sets, files = [], []
    for a in analysis_sets_for(accession):
        acc = a.get("accession")
        if not acc:
            continue
        full = _get(f"/analysis-sets/{acc}/?format=json") or {}
        info = {
            "accession": acc,
            "uniform_pipeline_status": (full.get("uniform_pipeline_status")
                                        or a.get("uniform_pipeline_status") or ""),
            "status": full.get("status") or a.get("status") or "",
            "aliases": full.get("aliases") or a.get("aliases") or [],
            "workflows": [w.get("accession") or w.get("@id")
                          for w in full.get("workflows") or []],
            "summary": (full.get("summary") or "")[:160],
        }
        sets.append(info)
        for f in full.get("files") or []:
            files.append({
                "accession": f.get("accession"),
                "content_type": f.get("content_type") or "",
                "file_format": f.get("file_format") or "",
                "file_size": f.get("file_size"),
                "controlled_access": f.get("controlled_access"),
                "href": f.get("href") or "",
                "s3_uri": f.get("s3_uri") or "",
                "from_analysis_set": acc,
                "uniform": info["uniform_pipeline_status"] == "completed",
            })
    return {"accession": accession, "analysis_sets": sets, "files": files}


def rank_for(files: "list[dict]", want: str) -> "list[dict]":
    wanted = WANTS.get(want, (want,))
    hits = [f for f in files
            if any(w in (f["content_type"] or "").lower() for w in wanted)]
    # Completed uniform-pipeline output first -- it is the one IGVF stands
    # behind -- then the preferred format, then smaller files.
    def key(f):
        ct, fmt = (f["content_type"] or "").lower(), (f["file_format"] or "").lower()
        try:
            pref = MATRIX_RANK.index((ct, fmt))
        except ValueError:
            pref = len(MATRIX_RANK)
        return (not f["uniform"], pref, f.get("file_size") or 0)
    return sorted(hits, key=key)


def _gb(n) -> str:
    try:
        return f"{float(n) / 1e9:.1f} GB"
    except (TypeError, ValueError):
        return "?"


# --------------------------------------------------------------------------


def cmd_find(args) -> int:
    res = processed_files(args.accession)
    sets, files = res["analysis_sets"], res["files"]

    if args.json:
        print(json.dumps(res, indent=2))
    if not sets:
        print(f"{args.accession}: no analysis sets — nothing has been "
              f"processed from this yet.\n"
              f"Running the pipeline yourself is the only route.")
        return 0

    print(f"{args.accession} — {len(sets)} analysis set(s), "
          f"{len(files)} processed file(s)\n")
    for s in sets:
        flag = ("UNIFORM PIPELINE COMPLETED"
                if s["uniform_pipeline_status"] == "completed"
                else f"pipeline: {s['uniform_pipeline_status'] or 'not stated'}")
        print(f"  {s['accession']}  [{flag}]")
        if s["workflows"]:
            print(f"    workflow: {', '.join(s['workflows'])}")
        if s["summary"]:
            print(f"    {s['summary']}")
    print()
    if files:
        print(f"  {'file':16} {'content type':30} {'fmt':6} {'size':>9}  access")
        for f in sorted(files, key=lambda x: -(x.get("file_size") or 0)):
            acc = ("controlled" if f["controlled_access"] is True
                   else "public" if f["controlled_access"] is False
                   else "not stated")
            print(f"  {f['accession'] or '?':16} {f['content_type'][:30]:30} "
                  f"{f['file_format'][:6]:6} {_gb(f['file_size']):>9}  {acc}")

    best = rank_for(files, "matrix")
    if best:
        b = best[0]
        print(f"\n  A processed matrix already exists: {b['accession']} "
              f"({b['file_format']}, {_gb(b['file_size'])}).")
        print(f"  Fetch it rather than re-running alignment:")
        print(f"    igvfagent client download-file {b['accession']}")
    ls.record_analysis("processed", subcommand="find", label=args.accession,
                       inputs=[args.accession],
                       outputs=[f["accession"] for f in files if f["accession"]])
    return 0


def cmd_plan(args) -> int:
    res = processed_files(args.accession)
    ranked = rank_for(res["files"], args.want)
    uniform = [s for s in res["analysis_sets"]
               if s["uniform_pipeline_status"] == "completed"]

    print(f"{args.accession} — want: {args.want}\n")
    if ranked:
        b = ranked[0]
        print(f"  USE THE EXISTING FILE. {b['accession']} "
              f"({b['content_type']}, {b['file_format']}, "
              f"{_gb(b['file_size'])})")
        print(f"    from analysis set {b['from_analysis_set']}"
              + ("  [uniform pipeline, completed]" if b["uniform"] else ""))
        # Three states, not two. The Portal omits controlled_access on some
        # files, and reading a missing field as "public" would send someone
        # to a download that 403s. Say "not stated" and let them check.
        if b["controlled_access"] is True:
            print("    controlled access — needs credentials and an approved "
                  "data-use agreement")
        elif b["controlled_access"] is False:
            print("    public — no credentials needed")
        else:
            print("    access not stated on the file record — check before "
                  "assuming it is downloadable")
        print(f"\n    igvfagent client download-file {b['accession']}")
        if len(ranked) > 1:
            print(f"\n  {len(ranked) - 1} other candidate(s): "
                  + ", ".join(f"{f['accession']} ({f['file_format']})"
                              for f in ranked[1:4]))
        print("\n  Running the pipeline from raw reads would recompute this.")
    elif uniform:
        print(f"  The uniform pipeline completed ({uniform[0]['accession']}), "
              f"but produced nothing matching '{args.want}'.")
        print(f"  Available: "
              + ", ".join(sorted({f['content_type'] for f in res['files']})))
        print(f"\n  Run the pipeline for this product:")
        print(f"    igvfagent raw-pipeline plan {args.accession}")
    else:
        print("  Nothing processed. Run the pipeline:")
        print(f"    igvfagent raw-pipeline plan {args.accession}")
    if args.json:
        print()
        print(json.dumps({"want": args.want, "candidates": ranked,
                          "uniform_sets": uniform}, indent=2))
    ls.record_analysis("processed", subcommand="plan", label=args.accession,
                       inputs=[args.accession],
                       outputs=[f["accession"] for f in ranked[:1]])
    return 0


def _have_credentials() -> bool:
    try:
        from _credentials import portal_credentials  # type: ignore
        return bool(portal_credentials())
    except Exception:
        return False


def _run_dir(label: str) -> Path:
    import time as _t
    d = OUT_DIR / f"{_t.strftime('%Y%m%d_%H%M%S')}_{label}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cmd_lineage(args) -> int:
    import portal_lineage as pl  # noqa: E402
    wants = [w.strip() for w in (args.want or "").split(",") if w.strip()] or None
    g = pl.walk(args.accession, fetch=portal_json, max_depth=args.depth, max_nodes=args.max_nodes,
                fetch_qc=not args.no_qc, fanout=args.fanout)
    if not g.nodes or (len(g.nodes) == 1 and g.blocked):
        code = next(iter(g.blocked.values()), 0) if g.blocked else 0
        print(f"{args.accession}: the Portal did not return this accession (HTTP {code}).")
        if code == 403:
            print("  It exists but is not visible with the current credentials (unreleased or controlled). "
                  "Run `igvfagent auth-check`.")
        return 2
    plan = pl.build_plan(g, wants=wants, have_credentials=_have_credentials())
    out = _run_dir(f"{args.accession}_lineage")
    (out / "lineage_graph.json").write_text(json.dumps(g.to_json(), indent=2, default=str))
    (out / "plan.json").write_text(json.dumps(plan, indent=2, default=str))
    fig = None if args.no_plots else pl.draw(g, out / "lineage.png")
    rep = pl.write_report(plan, g, out / "report.md", fig)
    r = plan["root"]
    print(f"{r['accession']} ({r.get('type')}): {plan['n_nodes']} linked objects, {plan['n_edges']} links, "
          f"{len(plan['blocked'])} not visible")
    if plan.get("node_cap_reached"):
        print(f"  note: the {args.max_nodes}-object cap was reached; rerun with --max-nodes 600 to see distant links")
    if plan["verdict"] == "processed":
        print(f"PROCESSED RESULTS EXIST: start from these instead of {plan['raw_gb']:.1f} GB of raw reads")
    elif plan["verdict"] == "raw_only":
        print(f"ONLY RAW DATA: {plan['raw_gb']:.1f} GB of reads; no processed outputs are linked yet")
    for sh in plan["start_here"]:
        qc = "; ".join(f"{k} {v}" for k, v in list(sh["qc_headline"].items())[:4])
        print(f"  {sh['product']:22} {sh['accession']:16} {sh['file_format']:6} {sh['size_gb']:8.2f} GB  "
              f"{sh['access']:12} from {sh['file_set']} ({sh['file_set_type'] or '-'}"
              f"{', uniform' if sh['uniform_pipeline_status'] == 'completed' else ''})" + (f"  QC: {qc}" if qc else ""))
    use = [m for m in plan.get("predictions_and_models") or [] if m.get("relation") == "model_use"]
    for m in plan.get("predictions_and_models") or []:
        if m.get("relation") != "model_use":
            print(f"  {m['type']:14} {m['accession']}  {m.get('file_set_type') or ''}  "
                  f"{m.get('scope') or m.get('model_name') or ''}  (inputs: {', '.join(m['inputs'][:3])})")
    if use:
        print(f"  {len(use)} prediction set(s) apply a model trained on this data to other data: "
              + ", ".join(m["accession"] for m in use[:8]) + (" ..." if len(use) > 8 else ""))
    docs = [b for b in plan["blocked"] if b["type"] == "Document"]
    for b in plan["blocked"]:
        if b["type"] != "Document":
            print(f"  NOT VISIBLE    {b['accession']} ({b['type']}, HTTP {b['http_status']}): {b['meaning']}")
    if docs:
        print(f"  NOT VISIBLE    {len(docs)} document(s): HTTP 403 (not released)")
    print(f"Report: {rep}")
    print(f"JSON: {out / 'plan.json'}")
    print(f"JSON: {out / 'lineage_graph.json'}")
    if fig:
        print(f"Figure: {fig}")
    ls.record_analysis("processed", subcommand="lineage", label=args.accession, inputs=[args.accession],
                       outputs=[s["accession"] for s in plan["start_here"]])
    return 0


def cmd_fetch(args) -> int:
    import portal_lineage as pl  # noqa: E402
    try:
        from igvfagent.raw_data_pipeline import portal_download  # type: ignore
    except Exception:
        from raw_data_pipeline import portal_download  # type: ignore
    wants = [w.strip() for w in args.want.split(",") if w.strip()]
    g = pl.walk(args.accession, fetch=portal_json, max_depth=args.depth, fetch_qc=True)
    creds = _have_credentials()
    plan = pl.build_plan(g, wants=wants, have_credentials=creds)
    dest = ROOT / "Data" / "Processed" / args.accession
    dest.mkdir(parents=True, exist_ok=True)
    budget = args.max_gb
    got, skipped = [], []
    for sh in plan["start_here"]:
        name = sh["href"].rsplit("/", 1)[-1] if sh.get("href") else f"{sh['accession']}.{sh['file_format']}"
        target = dest / name
        if "controlled" in sh["access"] and not creds:
            skipped.append((sh, "controlled access and no IGVF credentials configured"))
            continue
        if sh["size_gb"] > budget:
            skipped.append((sh, f"{sh['size_gb']:.1f} GB exceeds the remaining budget {budget:.1f} GB"))
            continue
        if not sh.get("href"):
            skipped.append((sh, "no download href on the file record"))
            continue
        if args.dry_run:
            print(f"  would fetch {sh['product']:18} {sh['accession']} {sh['size_gb']:.2f} GB -> {target}")
            continue
        if target.is_file() and target.stat().st_size > 0:
            print(f"  already here {target}")
        else:
            print(f"  fetching {sh['product']} {sh['accession']} ({sh['size_gb']:.2f} GB) ...", flush=True)
            portal_download(sh["href"], target)
        budget -= sh["size_gb"]
        got.append((sh, target))
        print(f"Wrote: {target}")
    # the Portal's QC attachments (kb_info.json, barcode summaries ...) are small; take them too
    qc_saved = []
    if not args.dry_run and not args.no_qc:
        for q in (n for n in g.nodes.values() if n.get("kind") == "qc"):
            (dest / "qc").mkdir(exist_ok=True)
            base = q["@id"].rstrip("/")
            (dest / "qc" / f"{q['type']}_{q['accession']}.json").write_text(json.dumps(q, indent=2, default=str))
            for field, href in (q.get("attachments") or {}).items():
                t = dest / "qc" / f"{q['accession']}_{href.rsplit('/', 1)[-1]}"
                try:
                    portal_download(f"{base}/{href}", t)
                    qc_saved.append(t)
                except Exception as exc:  # noqa: BLE001
                    print(f"  QC attachment {href}: {exc}")
    for sh, why in skipped:
        print(f"  SKIPPED {sh['product']} {sh['accession']}: {why}")
    (dest / "fetch_manifest.json").write_text(json.dumps({"accession": args.accession, "fetched": [
        {"product": sh["product"], "accession": sh["accession"], "path": str(t), "size_gb": sh["size_gb"]} for sh, t in got],
        "skipped": [{"product": sh["product"], "accession": sh["accession"], "why": why} for sh, why in skipped],
        "qc": [str(x) for x in qc_saved]}, indent=2))
    print(f"JSON: {dest / 'fetch_manifest.json'}")
    return 0 if got or args.dry_run else 1


def cmd_discover(args) -> int:
    """Topic search: which Portal objects are about this phenotype / tissue / gene / type."""
    try:
        from igvfagent import portal_discover as pdisc  # type: ignore
    except Exception:
        import portal_discover as pdisc  # type: ignore
    crit = {k: v for k, v in (("phenotype", args.phenotype), ("tissue", args.tissue), ("gene", args.gene),
                               ("prediction_type", args.prediction_type), ("library_type", args.library_type),
                               ("assay", args.assay), ("curated_type", args.curated_type)) if v}
    if not crit:
        print("give at least one of --phenotype --tissue --gene --prediction-type --library-type --assay "
              "--curated-type")
        return 2
    types = [t.strip() for t in (args.types or "").split(",") if t.strip()] or None
    res = pdisc.discover(crit, types=types, limit=args.limit, include_superseded=args.include_superseded)
    import re as _re
    out = _run_dir("discover_" + _re.sub(r"[^A-Za-z0-9._-]+", "_", "_".join(str(v) for v in crit.values()))[:60])
    (out / "discover.json").write_text(json.dumps(res, indent=2))
    (out / "report.md").write_text(pdisc.to_markdown(res))
    for typ, t in res["types"].items():
        m = "; ".join(f"{k}: {', '.join(v)}" for k, v in (t.get("matched") or {}).items())
        print(f"{typ:20s} {t['total']:5d}  {m}")
        for it in t["items"][:5]:
            print(f"    {it['accession']}  {it['file_set_type'][:28]:28s} {it['summary'][:80]}")
        for k, v in (t.get("no_match") or {}).items():
            print(f"    no {k} match; nearest: {', '.join(v[:5])}")
    print(f"JSON: {out / 'discover.json'}")
    print(f"Report: {out / 'report.md'}")
    return 0


def cmd_selftest(args) -> int:
    import tempfile
    import portal_lineage as pl  # noqa: E402
    import portal_discover as pdisc  # noqa: E402
    checks: "list" = []
    with tempfile.TemporaryDirectory() as td:
        pl.selftest(checks, Path(td))
    for good, msg in pdisc.selftest():
        checks.append((good, msg))
        print(("  ok    " if good else "  FAIL  ") + msg)
    ok = all(c for c, _ in checks)
    print("selftest: all checks pass" if ok else f"selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent processed",
        description="What IGVF has already computed from an accession — so a "
                    "workflow starts from processed output instead of "
                    "recomputing it from raw reads.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("find", help="Analysis sets and processed files for an "
                                    "accession.")
    s.add_argument("accession")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_find)

    s = sub.add_parser("plan", help="Download an existing file, or run the "
                                    "pipeline? Answers for one product.")
    s.add_argument("accession")
    s.add_argument("--want", default="matrix",
                   choices=sorted(WANTS), help="What the workflow needs.")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_plan)

    s = sub.add_parser("lineage", help="Walk every Portal link from an accession: processed outputs, QC, "
                                       "partners, samples, predictions; report + figure + JSON.")
    s.add_argument("accession")
    s.add_argument("--want", default="", help="Comma-separated products to list (default: all), e.g. rna_matrix,fragments")
    s.add_argument("--depth", type=int, default=4)
    s.add_argument("--max-nodes", type=int, default=300)
    s.add_argument("--fanout", type=int, default=25, help="Max links followed out of one shared hub object")
    s.add_argument("--no-qc", action="store_true", help="Skip fetching QualityMetric objects")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_lineage)

    s = sub.add_parser("fetch", help="Download the best processed file per product (and the Portal QC).")
    s.add_argument("accession")
    s.add_argument("--want", default="rna_matrix,atac_matrix,fragments,peaks,cell_annotations,demultiplexing,predictions")
    s.add_argument("--max-gb", type=float, default=20.0, help="Total download budget in GB")
    s.add_argument("--depth", type=int, default=4)
    s.add_argument("--no-qc", action="store_true")
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_fetch)

    s = sub.add_parser("discover", help="Topic search: Portal objects about a phenotype, tissue, gene, "
                                        "prediction type, library type or assay, matched to the Portal's own terms.")
    s.add_argument("--phenotype", help="e.g. 'coronary artery disease', 'LDL cholesterol'")
    s.add_argument("--tissue", help="sample term, e.g. heart, liver, K562")
    s.add_argument("--gene", help="assessed / targeted gene symbol")
    s.add_argument("--prediction-type", help="PredictionSet/ModelSet type, e.g. 'element-gene links'")
    s.add_argument("--library-type", help="construct library type, e.g. 'guide library', 'reporter library'")
    s.add_argument("--assay", help="assay title, e.g. 'MPRA', '10x multiome'")
    s.add_argument("--curated-type", help="CuratedSet type, e.g. variants, elements")
    s.add_argument("--types", default="", help="Comma-separated Portal types (default: prediction, model, "
                                               "measurement, analysis, curated and construct library sets)")
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--include-superseded", action="store_true")
    s.set_defaults(func=cmd_discover)

    s = sub.add_parser("selftest", help="Offline test of the graph walk and plan on a fixture Portal.")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_selftest)
    return p


def main(argv=None) -> int:
    setup_logging()
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
