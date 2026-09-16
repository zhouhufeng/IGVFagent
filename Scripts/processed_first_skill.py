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
    """Every processed output reachable from this accession."""
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
    return p


def main(argv=None) -> int:
    setup_logging()
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
