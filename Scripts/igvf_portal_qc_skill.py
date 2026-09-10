#!/usr/bin/env python3
"""
igvf-qc -- portal-wide quality control for IGVF data.

Where igvf_sub.py checks one submission before you send it, this checks what is
already on the portal: across all labs, all file types, at whatever scope you
ask for.

Two different sources of truth are combined:

  1. The portal's OWN audits, which are already computed and exposed as search
     facets (`audit.ERROR.category`, `audit.NOT_COMPLIANT.category`,
     `audit.INTERNAL_ACTION.category`). Reading facet counts gives a complete
     picture in one request instead of 150,000.

  2. Provenance integrity, re-ranked. The portal already audits most of this as
     `mismatched status`, but at INTERNAL_ACTION severity -- the DCC-facing
     band -- mixed in with ~4,000 other status mismatches, and without naming
     the replacement file.

     It also has one blind spot, measured against prod in Sep 2026: of 55 live
     files derived from an unusable file, 30 carried no audit at all, and all
     30 were the same shape -- `status=released` deriving from
     `status=archived`. Every other combination is audited. That looks
     deliberate (archived provenance may be considered legitimate), so this
     tool reports revoked and archived provenance separately and lets you pick
     the policy with --dead-status.

Zero dependencies beyond igvf_sub.py sitting alongside it.
"""

from __future__ import annotations

import argparse
import base64
import collections
import csv
import http.client
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import igvf_submission_skill as sub                              # noqa: E402
from igvf_submission_skill import (  # noqa: E402
    Portal, PortalError, DEAD_STATUSES, acc_of, die, fail, info, ok, warn, _c,
)

__version__ = "1.0.0"

def _default_cache() -> str:
    """Where the portal-wide graph pull is cached.

    ~/.cache is wrong in the deployment: the container runs as `igvf`, whose
    home may not be writable, and a per-user cache means each caller re-pulls
    a graph that takes minutes and is identical for everyone. IGVF_DATA_ROOT
    is the data volume the rest of IGVFagent writes to, so the cache is shared
    and survives container recreation. Falls back to ~/.cache standalone.
    """
    root = os.environ.get("IGVF_QC_CACHE")
    if root:
        return root
    for base in (os.environ.get("IGVF_DATA_ROOT"), "/workspace", None):
        if base and os.path.isdir(base) and os.access(base, os.W_OK):
            return os.path.join(base, "Data", "Cache", "igvf-qc")
    return os.path.expanduser("~/.cache/igvf-qc")


CACHE_DIR = _default_cache()
DEFAULT_MAX_AGE = 6 * 3600          # re-pull the graph at most every 6 hours

# Severity ordering the portal uses, most to least urgent.
SEVERITIES = ["ERROR", "NOT_COMPLIANT", "WARNING", "INTERNAL_ACTION"]
SEVERITY_LABEL = {
    "ERROR": "ERROR",
    "NOT_COMPLIANT": "NOT COMPLIANT",
    "WARNING": "WARNING",
    "INTERNAL_ACTION": "DCC ACTION",
}

# Abstract types worth sweeping by default. Each covers its concrete subtypes.
SURVEY_TYPES = ["File", "FileSet", "Sample", "Donor", "Document", "Software",
                "SoftwareVersion", "Workflow", "AnalysisStep", "AnalysisStepVersion"]


# ---------------------------------------------------------------------------
# Bulk fetch with an on-disk cache
# ---------------------------------------------------------------------------

def cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def _stream_get(portal: Portal, query: str, dest: str, attempts: int = 3) -> dict:
    """GET a search into a file, then parse it.

    The complete File scan is tens of megabytes over a chunked connection, and
    an interrupted read raises IncompleteRead after minutes of transfer. Stream
    it to disk and retry rather than buffering it all in memory.
    """
    url = portal.base.rstrip("/") + "/search/?" + query
    headers = {"Accept": "application/json", "User-Agent": f"igvf-qc/{__version__}"}
    if portal.authenticated:
        token = base64.b64encode(
            f"{portal.api_key}:{portal.secret_key}".encode()).decode()
        headers["Authorization"] = f"Basic {token}"
    tmp = dest + ".part"
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers=headers, method="GET")
            with urllib.request.urlopen(req, timeout=900) as resp, open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out, 1 << 20)
            with open(tmp) as fh:
                return json.load(fh)
        except (http.client.IncompleteRead, urllib.error.URLError,
                json.JSONDecodeError, TimeoutError) as exc:
            last = exc
            if attempt < attempts:
                wait = 2 ** attempt
                print(warn(f"transfer failed ({type(exc).__name__}); "
                           f"retrying in {wait}s"), file=sys.stderr)
                time.sleep(wait)
        finally:
            if os.path.exists(tmp) and attempt == attempts:
                os.remove(tmp)
    raise PortalError(f"could not fetch {url}: {last}")


def bulk_search(portal: Portal, item_type: str, fields, max_age=DEFAULT_MAX_AGE,
                refresh=False, extra=None, quiet=False):
    """One complete `limit=all` scan, cached on disk.

    `limit=all` rather than paging because the portal caps paging depth at
    99,999 and File alone has ~150,000 objects, so paging cannot enumerate it.
    Keep `fields` minimal: adding a few embedded properties took one scan from
    37 MB to 661 MB.
    """
    key = f"{portal.mode}.{item_type}.{'-'.join(sorted(fields))}"
    if extra:
        key += "." + "-".join(f"{k}={v}" for k, v in sorted(extra.items()))
    path = cache_path(key.replace("/", "_").replace("@", "") + ".json")

    if not refresh and os.path.exists(path):
        age = time.time() - os.path.getmtime(path)
        if age < max_age:
            if not quiet:
                print(info(f"using cached {item_type} ({age/60:.0f} min old; "
                           f"--refresh to re-pull)"), file=sys.stderr)
            with open(path) as fh:
                return json.load(fh)

    if not quiet:
        print(info(f"scanning all {item_type} objects on {portal.mode} "
                   f"(one request, can take a minute)..."), file=sys.stderr)
    query = "&".join([f"type={item_type}", "limit=all", "format=json"]
                     + [f"field={urllib.parse.quote(f, safe='.')}" for f in fields]
                     + [qparam(k, v) for k, v in (extra or {}).items()])
    res = _stream_get(portal, query, path)
    graph = res.get("@graph", [])
    with open(path, "w") as fh:
        json.dump(graph, fh)
    if not quiet:
        print(info(f"{len(graph):,} {item_type} objects cached"), file=sys.stderr)
    return graph


def qparam(key: str, value) -> str:
    """One query parameter, URL-encoded. Lab titles contain spaces and commas."""
    return f"{urllib.parse.quote(str(key), safe='.@')}=" \
           f"{urllib.parse.quote(str(value), safe='.@')}"


def facets_for(portal: Portal, item_type: str, extra=None) -> dict:
    """Facet terms for a type, without pulling any objects (limit=0)."""
    query = f"type={urllib.parse.quote(item_type)}&limit=0&format=json"
    if extra:
        query += "&" + "&".join(qparam(k, v) for k, v in extra.items())
    res = portal._request("GET", "/search/?" + query)
    out = {"total": res.get("total", 0), "facets": {}}
    for f in res.get("facets", []):
        # A facet's terms are usually [{"key":..., "doc_count":...}], but some
        # facet types return bare strings or omit the count.
        terms = {}
        for term in f.get("terms") or []:
            if not isinstance(term, dict):
                continue
            key, count = term.get("key"), term.get("doc_count")
            if key is None or not count:
                continue
            terms[str(key)] = count
        out["facets"][f.get("field")] = {"title": f.get("title"), "terms": terms}
    return out


def acc(obj: dict) -> str:
    """Accession from an object that may only carry @id."""
    return obj.get("accession") or acc_of(obj.get("@id", ""))


def lab_name(obj: dict) -> str:
    lab = obj.get("lab")
    if isinstance(lab, dict):
        return lab.get("title") or acc_of(lab.get("@id", "")) or "?"
    return acc_of(lab) if lab else "?"


def bar(n: int, total: int, width: int = 24) -> str:
    if not total:
        return ""
    filled = int(round(width * n / total))
    return "█" * filled + "─" * (width - filled)


# ---------------------------------------------------------------------------
# qc summary -- portal inventory and audit rollup
# ---------------------------------------------------------------------------

def cmd_summary(portal: Portal, args) -> int:
    types = args.type or SURVEY_TYPES
    extra = {"lab.title": args.lab} if args.lab else None
    scope = f" (lab.title={args.lab})" if args.lab else ""
    print(_c(f"IGVF portal inventory -- {portal.mode}{scope}", "1"))
    print()
    rows, all_audits = [], collections.Counter()
    for t in types:
        try:
            f = facets_for(portal, t, extra)
        except PortalError as exc:
            print(warn(f"{t}: {exc}"))
            continue
        total = f["total"]
        audit_total = 0
        for sev in SEVERITIES:
            terms = f["facets"].get(f"audit.{sev}.category", {}).get("terms", {})
            for cat, n in terms.items():
                all_audits[(sev, t, cat)] += n
                if sev in ("ERROR", "NOT_COMPLIANT"):
                    audit_total += n
        rows.append((t, total, audit_total, f))
    w = max((len(r[0]) for r in rows), default=8)
    print(f"  {'type':<{w}}  {'objects':>9}  {'flagged':>8}")
    for t, total, flagged, _ in rows:
        mark = _c(f"{flagged:>8}", "31;1") if flagged else f"{flagged:>8}"
        print(f"  {t:<{w}}  {total:>9,}  {mark}")
    print()
    print(info("'flagged' counts ERROR + NOT COMPLIANT audits raised by the portal."))

    for t, total, _, f in rows:
        st = f["facets"].get("status", {}).get("terms", {})
        if not st:
            continue
        dead = sum(v for k, v in st.items() if k in DEAD_STATUSES)
        if dead or args.verbose:
            print()
            print(_c(f"  {t} status", "1"))
            for k, v in sorted(st.items(), key=lambda kv: -kv[1]):
                flagd = "  <- unusable as an input" if k in DEAD_STATUSES else ""
                print(f"    {k:<20} {v:>8,}  {bar(v, total)}{flagd}")

    if all_audits:
        print()
        print(_c("portal audits by severity", "1"))
        for sev in SEVERITIES:
            items = {(t, c): n for (s, t, c), n in all_audits.items() if s == sev}
            if not items:
                continue
            head = fail if sev == "ERROR" else warn if sev == "NOT_COMPLIANT" else info
            print()
            print(head(f"{SEVERITY_LABEL[sev]}  ({sum(items.values()):,} total)"))
            for (t, c), n in sorted(items.items(), key=lambda kv: -kv[1]):
                print(f"    {n:>7,}  {t:<12} {c}")
    return 0


# ---------------------------------------------------------------------------
# qc audits -- drill into the portal's own audits
# ---------------------------------------------------------------------------

def cmd_audits(portal: Portal, args) -> int:
    item_type = args.item_type
    extra = {}
    if args.lab:
        extra["lab.title"] = args.lab
    f = facets_for(portal, item_type, extra or None)
    print(_c(f"portal audits on {item_type}"
             + (f" for lab.title={args.lab}" if args.lab else ""), "1"))
    print(f"  {f['total']:,} objects in scope")
    print()
    found = False
    for sev in SEVERITIES:
        if args.severity and sev != args.severity:
            continue
        terms = f["facets"].get(f"audit.{sev}.category", {}).get("terms", {})
        if not terms:
            continue
        found = True
        head = fail if sev == "ERROR" else warn if sev == "NOT_COMPLIANT" else info
        print(head(SEVERITY_LABEL[sev]))
        for cat, n in sorted(terms.items(), key=lambda kv: -kv[1]):
            print(f"    {n:>7,}  {cat}")
            if args.examples:
                q = (f"type={urllib.parse.quote(item_type)}&limit={args.examples}"
                     f"&format=json&"
                     + qparam(f"audit.{sev}.category", cat)
                     + "&field=accession&field=status&field=lab.title")
                if args.lab:
                    q += "&" + qparam("lab.title", args.lab)
                try:
                    res = portal._request("GET", "/search/?" + q)
                    for g in res.get("@graph", []):
                        print(f"             {g.get('accession')}  "
                              f"{g.get('status'):<12} {lab_name(g)}")
                    print(f"             {portal.ui_base}/search/?type={item_type}"
                          f"&audit.{sev}.category={cat.replace(' ', '+')}")
                except PortalError as exc:
                    print(f"             (examples unavailable: {exc})")
        print()
    if not found:
        print(ok("no audits in this scope."))
    return 0


# ---------------------------------------------------------------------------
# qc provenance -- the check the portal does not run
# ---------------------------------------------------------------------------

def build_file_graph(portal: Portal, args):
    # Deliberately minimal: accession is derivable from @id, and each extra
    # embedded property multiplies the scan size.
    fields = ["@id", "status", "lab", "derived_from"]
    graph = bulk_search(portal, "File", fields, refresh=args.refresh,
                        max_age=0 if args.refresh else DEFAULT_MAX_AGE)
    return {x["@id"]: x for x in graph}


def cmd_provenance(portal: Portal, args) -> int:
    by_id = build_file_graph(portal, args)
    dead_statuses = set(args.dead_status or DEAD_STATUSES)
    dead = {i: x.get("status") for i, x in by_id.items()
            if x.get("status") in dead_statuses}
    print()
    print(_c("provenance integrity", "1"))
    print(f"  {len(by_id):,} files, {len(dead):,} in a status counted as unusable "
          f"({', '.join(sorted(dead_statuses))})")

    offenders = []
    for x in by_id.values():
        if x.get("status") in dead_statuses:
            continue
        bad = [r for r in (x.get("derived_from") or []) if r in dead]
        if bad:
            offenders.append((x, bad))
    if args.lab:
        needle = args.lab.lower()
        offenders = [(x, b) for x, b in offenders if needle in lab_name(x).lower()]
    if args.released_only:
        offenders = [(x, b) for x, b in offenders if x.get("status") == "released"]

    if not offenders:
        print()
        print(ok("no live file depends on an unusable input in this scope."))
        return 0

    # Split by the status of the worst input. A revoked input is a retraction;
    # an archived one may be legitimate retired provenance, and the portal
    # deliberately does not audit released-from-archived.
    def worst(bad):
        return "revoked" if any(dead[b] == "revoked" for b in bad) else "archived"

    revoked_prov = [(x, b) for x, b in offenders if worst(b) == "revoked"]
    archived_prov = [(x, b) for x, b in offenders if worst(b) != "revoked"]

    print(f"  {len(offenders):,} live file(s) derive from one")
    print()

    def block(title, items, level, note):
        if not items:
            return
        head = fail if level == "fail" else warn
        rel = [o for o in items if o[0].get("status") == "released"]
        print(head(f"{title}: {len(items)} file(s), {len(rel)} already released"))
        print(f"       {note}")
        print()
        per_lab = collections.Counter(lab_name(x) for x, _ in items)
        for lab, n in per_lab.most_common():
            nrel = sum(1 for x, _ in items
                       if lab_name(x) == lab and x.get("status") == "released")
            print(f"    {n:>5}  ({nrel} released)  {lab}")
        print()
        for x, bad in sorted(items, key=lambda o: (lab_name(o[0]), acc(o[0])))[:args.limit]:
            first = by_id[bad[0]]
            more = f"  +{len(bad)-1} more" if len(bad) > 1 else ""
            # Severity follows the block, not the row: an archived-provenance
            # row is a question even when the file is released.
            if x.get("status") == "released":
                lead = head("")
            else:
                lead = warn("") if level == "fail" else info("")
            print(f"{lead} {acc(x):<16} {x.get('status'):<12} {lab_name(x)[:30]:<30}")
            print(f"        derived_from {acc(first)} [{first.get('status')}]{more}")
            sup = first.get("superseded_by") or []
            if sup:
                print(f"        replacement: {', '.join(acc_of(s) for s in sup)}")
        if len(items) > args.limit:
            print(info(f"    ...{len(items) - args.limit} more (raise --limit, or --csv)"))
        print()

    block("REVOKED provenance", revoked_prov, "fail",
          "A revoked input was retracted, usually for a data error. Repoint "
          "derived_from to the superseding file. The portal audits these as "
          "'mismatched status' under DCC ACTION.")
    block("ARCHIVED provenance", archived_prov, "warn",
          "Archived inputs may be legitimate retired provenance -- the portal "
          "does NOT audit released-from-archived, so treat this as a question, "
          "not a defect. Confirm the archived file is still the right citation.")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["accession", "status", "lab", "provenance_class",
                         "dead_input", "dead_input_status", "n_dead_inputs",
                         "all_dead_inputs"])
            for x, bad in offenders:
                d = by_id[bad[0]]
                wr.writerow([acc(x), x.get("status"), lab_name(x), worst(bad),
                             acc(d), d.get("status"), len(bad),
                             " ".join(acc_of(b) for b in bad)])
        print(ok(f"wrote {args.csv} ({len(offenders)} rows)"))
    return 1 if revoked_prov else 0


# ---------------------------------------------------------------------------
# qc uploads -- files whose bytes never validated
# ---------------------------------------------------------------------------

def cmd_uploads(portal: Portal, args) -> int:
    extra = {"lab.title": args.lab} if args.lab else None
    f = facets_for(portal, "File", extra)
    terms = f["facets"].get("upload_status", {}).get("terms", {})
    print(_c("upload status" + (f" -- lab.title={args.lab}" if args.lab else ""), "1"))
    print(f"  {f['total']:,} files in scope")
    print()
    bad_total = 0
    for k, v in sorted(terms.items(), key=lambda kv: -kv[1]):
        healthy = k in ("validated", "validation exempted")
        if not healthy:
            bad_total += v
        mark = ok("") if healthy else fail("")
        print(f"{mark} {k:<22} {v:>8,}  {bar(v, f['total'])}")
    print()
    if bad_total:
        print(fail(f"{bad_total:,} file(s) have not validated. 'invalidated' means the "
                   f"portal checked the bytes and rejected them; 'pending' means it has "
                   f"not finished; 'file not found' means the upload never landed."))
        for state in ("file not found", "invalidated"):
            if terms.get(state):
                print()
                print(_c(f"  examples: {state}", "1"))
                q = (f"type=File&limit={args.examples}&format=json&"
                     + qparam("upload_status", state)
                     + "&field=accession&field=status&field=lab.title"
                       "&field=submitted_file_name")
                if args.lab:
                    q += "&" + qparam("lab.title", args.lab)
                try:
                    res = portal._request("GET", "/search/?" + q)
                    for g in res.get("@graph", []):
                        print(f"      {g.get('accession')}  {g.get('status'):<12} "
                              f"{lab_name(g)[:28]:<28} {g.get('submitted_file_name') or ''}")
                except PortalError as exc:
                    print(f"      (unavailable: {exc})")
        return 1
    print(ok("every file in scope has validated."))
    return 0


# ---------------------------------------------------------------------------
# qc filesets -- file sets missing the things that block release
# ---------------------------------------------------------------------------

def cmd_filesets(portal: Portal, args) -> int:
    fields = ["@id", "accession", "status", "lab.title", "description",
              "input_file_sets", "files", "file_set_type", "samples", "donors"]
    graph = bulk_search(portal, args.item_type, fields, refresh=args.refresh,
                        max_age=0 if args.refresh else DEFAULT_MAX_AGE)
    if args.lab:
        needle = args.lab.lower()
        graph = [x for x in graph if needle in lab_name(x).lower()]

    checks = collections.OrderedDict([
        ("no description", lambda x: not (x.get("description") or "").strip()),
        ("no files", lambda x: not x.get("files")),
        ("no samples or donors", lambda x: not x.get("samples") and not x.get("donors")),
    ])
    if args.item_type == "PredictionSet":
        checks["no input_file_sets"] = lambda x: not x.get("input_file_sets")

    print()
    print(_c(f"{args.item_type} completeness"
             + (f" -- lab {args.lab}" if args.lab else ""), "1"))
    print(f"  {len(graph):,} object(s) in scope")
    print()
    hits = {}
    for label, test in checks.items():
        bad = [x for x in graph if test(x)]
        hits[label] = bad
        mark = ok("") if not bad else fail("") if len(bad) > len(graph) * 0.2 else warn("")
        pct = (100.0 * len(bad) / len(graph)) if graph else 0
        print(f"{mark} {label:<24} {len(bad):>6,} / {len(graph):,}  ({pct:.0f}%)")

    for label, bad in hits.items():
        if not bad:
            continue
        print()
        print(_c(f"  {label}", "1"))
        for x in bad[:args.limit]:
            print(f"    {x.get('accession'):<16} {x.get('status'):<12} "
                  f"{lab_name(x)[:32]:<32} {x.get('file_set_type') or ''}")
        if len(bad) > args.limit:
            print(info(f"    ...{len(bad) - args.limit} more"))

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["accession", "status", "lab", "file_set_type", "problem"])
            for label, bad in hits.items():
                for x in bad:
                    wr.writerow([x.get("accession"), x.get("status"), lab_name(x),
                                 x.get("file_set_type"), label])
        print()
        print(ok(f"wrote {args.csv}"))
    return 1 if any(hits.values()) else 0


# ---------------------------------------------------------------------------
# qc lab -- one lab's scorecard
# ---------------------------------------------------------------------------

def cmd_lab(portal: Portal, args) -> int:
    lab = args.lab_title
    print(_c(f"QC scorecard -- {lab}  ({portal.mode})", "1"))
    print()
    for t in ("File", "FileSet"):
        f = facets_for(portal, t, {"lab.title": lab})
        if not f["total"]:
            print(warn(f"no {t} objects found for lab.title={lab!r}. "
                       f"Exact title required, e.g. 'Xihong Lin, HSPH'."))
            continue
        print(_c(f"{t}: {f['total']:,} objects", "1"))
        st = f["facets"].get("status", {}).get("terms", {})
        print("  status: " + ", ".join(f"{k}={v:,}" for k, v in
                                      sorted(st.items(), key=lambda kv: -kv[1])))
        up = f["facets"].get("upload_status", {}).get("terms", {})
        if up:
            print("  upload: " + ", ".join(f"{k}={v:,}" for k, v in
                                          sorted(up.items(), key=lambda kv: -kv[1])))
        any_audit = False
        for sev in SEVERITIES:
            terms = f["facets"].get(f"audit.{sev}.category", {}).get("terms", {})
            for cat, n in sorted(terms.items(), key=lambda kv: -kv[1]):
                any_audit = True
                head = fail if sev == "ERROR" else warn if sev == "NOT_COMPLIANT" else info
                print("  " + head(f"{SEVERITY_LABEL[sev]:<14} {n:>5,}  {cat}"))
        if not any_audit:
            print("  " + ok("no portal audits"))
        print()
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

EPILOG = """\
examples:
  igvf-qc summary                              # portal inventory + audit rollup
  igvf-qc audits File --severity ERROR --examples 5
  igvf-qc provenance --released-only           # released data with dead provenance
  igvf-qc uploads                              # files whose bytes never validated
  igvf-qc filesets PredictionSet --csv gaps.csv
  igvf-qc lab "Xihong Lin, HSPH"

Scope any command to one lab with --lab "<exact lab.title>".
The File graph is cached under ~/.cache/igvf-qc; --refresh re-pulls it.
"""


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="igvfagent portal-qc",
        description="Portal-wide quality control for IGVF data.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--mode", default=None,
                    help="prod | staging | hostname (default: $IGVF_MODE or prod)")
    ap.add_argument("-V", "--version", action="version", version=f"igvf-qc {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("summary", help="portal inventory and audit rollup")
    p.add_argument("--type", nargs="*", help=f"default: {' '.join(SURVEY_TYPES)}")
    p.add_argument("--lab", help="restrict to one lab.title")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="show status breakdown even when nothing is unusable")
    p.set_defaults(fn=cmd_summary)

    p = sub.add_parser("audits", help="drill into the portal's own audits")
    p.add_argument("item_type", nargs="?", default="File")
    p.add_argument("--severity", choices=SEVERITIES)
    p.add_argument("--lab")
    p.add_argument("--examples", type=int, default=0,
                   help="show N example accessions per category")
    p.set_defaults(fn=cmd_audits)

    p = sub.add_parser("provenance",
                       help="live files derived from revoked/archived files")
    p.add_argument("--lab")
    p.add_argument("--released-only", action="store_true",
                   help="only already-released offenders")
    p.add_argument("--limit", type=int, default=40)
    p.add_argument("--csv", help="write the full list to a CSV")
    p.add_argument("--refresh", action="store_true", help="re-pull the file graph")
    p.add_argument("--dead-status", nargs="*",
                   help=f"statuses to treat as unusable (default: "
                        f"{' '.join(sorted(DEAD_STATUSES))})")
    p.set_defaults(fn=cmd_provenance)

    p = sub.add_parser("uploads", help="files whose upload never validated")
    p.add_argument("--lab")
    p.add_argument("--examples", type=int, default=5)
    p.set_defaults(fn=cmd_uploads)

    p = sub.add_parser("filesets", help="file sets missing release-blocking properties")
    p.add_argument("item_type", nargs="?", default="PredictionSet")
    p.add_argument("--lab")
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--csv")
    p.add_argument("--refresh", action="store_true")
    p.set_defaults(fn=cmd_filesets)

    p = sub.add_parser("lab", help="one lab's QC scorecard")
    p.add_argument("lab_title", help='exact lab.title, e.g. "Xihong Lin, HSPH"')
    p.set_defaults(fn=cmd_lab)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    portal = Portal(mode=args.mode)
    try:
        return args.fn(portal, args)
    except PortalError as exc:
        die(str(exc))
    except KeyboardInterrupt:
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
