#!/usr/bin/env python3
"""Score IGVFagent benchmark runs against per-paper ``expected.json``.

For each benchmark directory under ``Benchmarks/<paper-id>/``:

1. Read ``expected.json`` (ground-truth ranges + qualitative checks).
2. Locate the latest output directory under ``Docs/<skill>/2*_<label>/``
   that matches the benchmark's label.
3. Open the canonical artefact (usually ``summary.json`` or a TSV).
4. Apply each declared check; tally pass / fail / skip.
5. Judge paper coverage (``reproduction``): ``status`` only says whether the
   declared checks passed, which a single easy count can satisfy. When
   ``expected.json`` lists the paper's ``analyses`` (``igvfagent bench plan``),
   a paper is ``reproduced`` only when every analysis has a passing,
   confirmed class-A/B check (``Benchmarks/taxonomy.py``) tied to it by
   ``"analysis": <id>``. Counts and file-exists checks never cover one.
6. Emit a per-paper JSON + a suite-level Markdown summary under
   ``Benchmarks/results/<ts>_concordance.{json,md}``.

The scorer is intentionally minimal — pure stdlib, no pandas required —
because the benchmark plan is meant to be reproducible by anyone who
clones IGVFagent and runs ``bash Benchmarks/run_all.sh``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from taxonomy import _check_class  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = Path(__file__).resolve().parent
RESULTS = BENCHMARKS / "results"
DOCS = ROOT / "Docs"


def _match_run_dir(base: Path, label: str,
                   require: str | None = None) -> Path | None:
    """Apply the two run-dir conventions against one base directory.

    ``require`` names the paper's primary artefact. A shared base such as
    ``Docs/Benchmark/`` also holds the benchmark tool's own per-paper dirs
    (resolve, harvest, route, report, port verification) under the same
    label, so the newest dir that holds the artefact wins over newer ones
    that do not. With no such dir the newest match is returned as before.
    """
    if not base.is_dir():
        return None
    # Convention 1: dir match
    dirs = sorted(
        (p for p in base.glob(f"2*_*{label}*") if p.is_dir()),
        key=lambda p: p.name,
        reverse=True,
    )
    if dirs:
        if require:
            for d in dirs:
                if (d / require).is_file():
                    return d
        return dirs[0]
    # Convention 2: flat-file match
    flat = sorted(
        (p for p in base.glob(f"2*_*{label}*") if p.is_file()),
        key=lambda p: p.name,
        reverse=True,
    )
    if flat:
        # Return the base dir itself; checks use the label to glob for files.
        return base
    return None


def latest_run_dir(skill_dir_name: str, label: str,
                    extra_search_dirs: "list[Path]" = (),
                    require: str | None = None) -> Path | None:
    """Find the most recent ``Docs/<skill>/2*_<label>*/`` directory.

    Two conventions are supported:

    1. Per-run directory: ``Docs/<skill>/<ts>_<label>/`` containing
       ``summary.json``, TSVs, plots, etc. (the convention used by
       the newer skills: mavedb, multiome, chipatlas, portal, catalog).
    2. Flat-file convention: ``Docs/<skill>/<ts>_<label>_*`` files
       sitting directly under the skill dir (older skills like the
       legacy ``mpra pull``). For these we return the skill dir itself
       and the artefact checks must match by glob pattern.

    Some skills (e.g. ``encode retrieve``) never write anything under
    ``Docs/<skill>/`` at all — their only output is a manifest under
    ``Data/Manifests/<skill>/``. When ``Docs/<skill>/`` has no match,
    fall back to each of ``extra_search_dirs`` (declared per-paper in
    ``expected.json``) using the same two conventions.
    """
    hit = _match_run_dir(DOCS / skill_dir_name, label, require)
    if hit is not None:
        return hit
    for extra in extra_search_dirs:
        hit = _match_run_dir(Path(extra), label, require)
        if hit is not None:
            return hit
    return None


def read_artefact(d: Path, filename: str) -> Any:
    """Read a JSON or TSV artefact from a run directory.

    A ``Docs/…`` or ``Data/…`` name is a glob from the repo root instead
    (newest match), for artefacts written outside the route's run dir, such
    as ``bench verify-port``'s ``validation_vs_reference.json``.
    """
    if filename.startswith(("Docs/", "Data/")):
        hits = sorted(ROOT.glob(filename), key=lambda q: q.stat().st_mtime)
        if not hits:
            return None
        d, filename = hits[-1].parent, hits[-1].name
    p = d / filename
    if not p.is_file():
        return None
    if filename.endswith(".json"):
        try:
            return json.loads(p.read_text())
        except Exception:
            return None
    if filename.endswith((".tsv", ".txt")):
        try:
            with p.open() as fh:
                rdr = csv.DictReader(fh, delimiter="\t")
                return list(rdr)
        except Exception:
            return None
    return p.read_text()


def check_range(value: Any, spec: dict) -> tuple[bool, str]:
    """Hard metric: a numeric value must fall in [min, max] (inclusive)."""
    if value is None:
        return False, "value is None"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return False, f"non-numeric value: {value!r}"
    lo = spec.get("min")
    hi = spec.get("max")
    exp = spec.get("expected")
    msgs: list[str] = []
    ok = True
    if lo is not None and v < lo:
        ok = False
        msgs.append(f"value {v} < min {lo}")
    if hi is not None and v > hi:
        ok = False
        msgs.append(f"value {v} > max {hi}")
    if ok:
        return True, f"{v} in [{lo}, {hi}]" + (f" (expected ≈{exp})" if exp is not None else "")
    return False, "; ".join(msgs)


def check_set_member(value: Any, spec: dict) -> tuple[bool, str]:
    """Qualitative: a string value must be in an allowed set."""
    allowed = spec.get("allowed") or spec.get("top_term_must_be_in") or []
    if not allowed:
        return False, "no allowed set declared"
    if value in allowed:
        return True, f"value {value!r} in allowed set"
    return False, f"value {value!r} not in allowed set {allowed}"


def check_artefact_exists(d: Path, spec: dict,
                           extra_search_dirs: "list[Path]" = ()) -> tuple[bool, str]:
    """Existence check: a named file must be present in the run dir.

    ``filename`` may be a bare name (matched as ``d/<name>``) or a glob
    pattern with ``*`` (matched against the run dir and against any
    extra search dirs — useful for older skills that scatter artefacts
    across ``Docs/<skill>/``, ``Data/Manifests/<skill>/``, and ``Data/``).
    """
    fname = spec.get("filename") or spec.get("artefact") or spec.get("glob")
    if not fname:
        return False, "no filename declared"
    candidates: list[Path] = [d] + [Path(x) for x in (extra_search_dirs or [])]
    if "*" in fname or "?" in fname:
        for c in candidates:
            if not c.is_dir():
                continue
            hits = sorted(c.glob(fname),
                           key=lambda p: p.stat().st_mtime, reverse=True)
            non_empty = [p for p in hits if p.stat().st_size > 0]
            if non_empty:
                p = non_empty[0]
                rel = p.relative_to(ROOT) if p.is_absolute() else p
                return True, f"{fname} matched → {rel} ({p.stat().st_size:,} bytes)"
        return False, (f"glob {fname!r} matched no non-empty files in any of "
                         f"{[str(c.relative_to(ROOT) if c.is_absolute() else c) for c in candidates]}")
    for c in candidates:
        p = c / fname
        if p.is_file() and p.stat().st_size > 0:
            rel = p.relative_to(ROOT) if p.is_absolute() else p
            return True, f"{fname} present ({p.stat().st_size:,} bytes) at {rel}"
    return False, f"{fname} not found in {[str(c.relative_to(ROOT) if c.is_absolute() else c) for c in candidates]}"


def resolve_glob_artefact(d: Path, pattern: str,
                           extra_search_dirs: "list[Path]" = ()) -> Path | None:
    """Return the newest non-empty file matching ``pattern`` across the run
    dir and any extra search dirs (same resolution order as
    ``check_artefact_exists``)."""
    candidates: list[Path] = [d] + [Path(x) for x in (extra_search_dirs or [])]
    for c in candidates:
        if not c.is_dir():
            continue
        hits = sorted(c.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
        non_empty = [p for p in hits if p.stat().st_size > 0]
        if non_empty:
            return non_empty[0]
    return None


def read_csv_rows(p: Path) -> list[dict] | None:
    try:
        with p.open(newline="") as fh:
            return list(csv.DictReader(fh))
    except Exception:
        return None


def get_path(obj: Any, dotted: str) -> Any:
    """Walk a dotted key path through a nested dict / list of dicts."""
    if obj is None:
        return None
    cur: Any = obj
    for part in dotted.split("."):
        if cur is None:
            return None
        if isinstance(cur, list):
            # Allow numeric index for list elements
            try:
                idx = int(part)
                cur = cur[idx] if 0 <= idx < len(cur) else None
                continue
            except ValueError:
                pass
            # Otherwise, treat as a key match on the first dict in the list
            cur = next(
                (e.get(part) for e in cur if isinstance(e, dict) and part in e),
                None,
            )
            continue
        if isinstance(cur, dict):
            cur = cur.get(part)
            continue
        return None
    return cur


# The only reasons an analysis may stay unreproduced without holding the paper
# back from "reproduced_except_access": the data cannot lawfully be had. Any
# other blocker (runtime, environment, language, data size) is work to do.
ACCESS_BLOCKERS = ("controlled_access", "embargoed", "not_deposited")


def judge_coverage(spec: dict, checks: list[dict]) -> dict:
    """Coverage of the paper's planned analyses by passing class-A/B checks.

    ``checks`` are the scored check results, each carrying the spec's
    ``analysis`` id and ``class``. Returns ``reproduction`` (reproduced /
    reproduced_except_access / incomplete / unplanned), ``coverage`` ("k/N")
    and one row per analysis.
    """
    analyses = spec.get("analyses") or []
    if not analyses:
        return {"reproduction": "unplanned", "coverage": None, "analyses": []}
    rows = []
    for a in analyses:
        aid = a.get("id")
        mine = [c for c in checks if c.get("analysis") == aid]
        strong = [c for c in mine if c["passed"] and c.get("class") in ("A", "B")]
        blocker = (a.get("blocker") or {}).get("kind") if isinstance(a.get("blocker"), dict) \
            else a.get("blocker")
        if strong:
            state = "reproduced"
        elif blocker in ACCESS_BLOCKERS:
            state = "blocked_access"
        elif mine and any(c["passed"] for c in mine):
            state = "weak"          # only counts / artefacts pass
        elif mine:
            state = "failing"
        else:
            state = "pending"
        rows.append({"id": aid, "title": a.get("title"), "state": state,
                     "blocker": blocker, "tool": a.get("tool"),
                     "n_checks": len(mine), "n_strong_passed": len(strong)})
    n = len(rows)
    k = sum(r["state"] == "reproduced" for r in rows)
    if k == n:
        verdict = "reproduced"
    elif k and all(r["state"] in ("reproduced", "blocked_access") for r in rows):
        verdict = "reproduced_except_access"
    else:
        verdict = "incomplete"
    return {"reproduction": verdict, "coverage": f"{k}/{n}", "analyses": rows}


def score_benchmark(paper_dir: Path) -> dict:
    """Score one paper's benchmark; return a structured result dict."""
    exp_path = paper_dir / "expected.json"
    if not exp_path.is_file():
        return {"paper": paper_dir.name, "status": "no_expected_json"}
    spec = json.loads(exp_path.read_text())
    skill = spec.get("skill_output_dir")            # e.g. "MaveDB", "ChIPAtlas"
    label = spec.get("label")
    if not skill or not label:
        return {"paper": paper_dir.name, "status": "expected_json_missing_skill_or_label"}
    # Extra search roots for skills that scatter artefacts. expected.json
    # may declare ``extra_search_dirs: ["Data/Manifests/MPRA", "Data"]`` etc.
    extras = [ROOT / x for x in (spec.get("extra_search_dirs") or [])]
    run_dir = latest_run_dir(skill, label, extras,
                             require=spec.get("primary_artefact"))
    if run_dir is None:
        return {"paper": paper_dir.name, "status": "no_run_found",
                "skill": skill, "label": label, **judge_coverage(spec, [])}
    result: dict[str, Any] = {
        "paper":   paper_dir.name,
        "skill":   skill,
        "label":   label,
        "run_dir": str(run_dir.relative_to(ROOT)),
        "checks":  [],
    }
    # Load the primary artefact once if declared
    primary_artefact = spec.get("primary_artefact", "summary.json")
    payload = read_artefact(run_dir, primary_artefact)

    for chk in spec.get("checks", []):
        ctype = chk.get("type")
        name = chk.get("name", "(unnamed)")
        # Checks scaffolded by `igvfagent bench` out of a paper's prose carry
        # ``"confirmed": false`` until a human sets a real JSON path and vouches
        # for them. They are reported, never scored — otherwise the suite would
        # be grading itself against text extraction.
        if chk.get("confirmed") is False:
            result.setdefault("unconfirmed", []).append({
                "name": name, "type": ctype,
                "expected": chk.get("expected"),
                "provenance": chk.get("provenance"),
            })
            continue
        try:
            # A chain that writes several summaries (spatial-hic emits
            # qc_ / demux_ / gas_ / gad_ / compartment_ / cnv_summary.json)
            # cannot be checked from the primary artefact alone. A check may
            # name its own file with ``"artefact": "demux_summary.json"``;
            # without this the path was looked up in the WRONG json and the
            # check failed as "value is None" rather than saying so.
            if chk.get("artefact"):
                src = read_artefact(run_dir, chk["artefact"])
                searched = [run_dir]
                for extra in extras:
                    if src is not None:
                        break
                    src = read_artefact(extra, chk["artefact"])
                    searched.append(extra)
                if src is None:
                    where = ", ".join(str(p) for p in searched)
                    result["checks"].append({
                        "name": name, "type": ctype, "passed": False,
                        "detail": f"artefact {chk['artefact']!r} not found in "
                                   f"{where}",
                        "class": _check_class(chk),
                        "analysis": chk.get("analysis")})
                    continue
            else:
                src = payload

            if ctype == "range":
                v = get_path(src, chk["path"])
                ok, msg = check_range(v, chk)
            elif ctype == "in_set":
                v = get_path(src, chk["path"])
                ok, msg = check_set_member(v, chk)
            elif ctype == "artefact":
                ok, msg = check_artefact_exists(run_dir, chk, extras)
            elif ctype == "row_count_tsv":
                rows = read_artefact(run_dir, chk["filename"])
                v = len(rows) if isinstance(rows, list) else None
                ok, msg = check_range(v, chk)
            elif ctype == "csv_row_count":
                pattern = chk.get("glob") or chk.get("filename")
                p = resolve_glob_artefact(run_dir, pattern, extras)
                if p is None:
                    ok, msg = False, f"glob {pattern!r} matched no file"
                else:
                    rows = read_csv_rows(p)
                    v = len(rows) if isinstance(rows, list) else None
                    ok, msg = check_range(v, chk)
                    msg = f"{msg} (in {p.relative_to(ROOT)})"
            elif ctype == "csv_value_present":
                pattern = chk.get("glob") or chk.get("filename")
                p = resolve_glob_artefact(run_dir, pattern, extras)
                if p is None:
                    ok, msg = False, f"glob {pattern!r} matched no file"
                else:
                    rows = read_csv_rows(p)
                    col, val = chk.get("column"), chk.get("value")
                    n_match = sum(1 for r in (rows or [])
                                   if isinstance(r, dict) and r.get(col) == val)
                    min_rows = chk.get("min_rows", 1)
                    ok = n_match >= min_rows
                    msg = (f"{n_match} rows with {col}={val!r} "
                            f"(need ≥{min_rows}) in {p.relative_to(ROOT)}")
            else:
                ok, msg = False, f"unknown check type {ctype!r}"
        except Exception as e:
            ok, msg = False, f"check failed with exception: {e}"
        result["checks"].append({"name": name, "type": ctype,
                                   "passed": ok, "detail": msg,
                                   "class": _check_class(chk),
                                   "analysis": chk.get("analysis")})
    n_total = len(result["checks"])
    n_pass = sum(1 for c in result["checks"] if c["passed"])
    result["n_passed"] = n_pass
    result["n_total"] = n_total
    result["n_unconfirmed"] = len(result.get("unconfirmed") or [])
    if n_total == 0:
        # Only unconfirmed checks exist — the chain may have run, but nothing
        # has been verified. Never report that as a pass.
        result["status"] = "unreviewed" if result["n_unconfirmed"] else "fail"
    else:
        result["status"] = "ok" if n_pass == n_total \
                             else ("partial" if n_pass > 0 else "fail")
    result.update(judge_coverage(spec, result["checks"]))
    return result


def render_markdown(results: list[dict], ts: str) -> str:
    out: list[str] = []
    out.append(f"# IGVFagent benchmark concordance — {ts}\n")
    n_tot = sum(r.get("n_total", 0) for r in results)
    n_pass = sum(r.get("n_passed", 0) for r in results)
    n_ok = sum(1 for r in results if r.get("status") == "ok")
    n_part = sum(1 for r in results if r.get("status") == "partial")
    n_fail = sum(1 for r in results
                  if r.get("status") in ("fail", "no_run_found",
                                            "no_expected_json"))
    n_unconf = sum(r.get("n_unconfirmed", 0) for r in results)
    out.append(f"**Suite summary:** {n_pass} / {n_tot} checks passed across "
                f"{len(results)} papers — {n_ok} clean ✓ · {n_part} partial · "
                f"{n_fail} unscored or failed.\n")
    if n_unconf:
        out.append(f"**{n_unconf} unconfirmed checks** were reported but not "
                    f"scored. These were extracted from paper text by "
                    f"`igvfagent bench` and need a human to set a real JSON path "
                    f"and flip `\"confirmed\": true` before they count.\n")
    n_repro = sum(1 for r in results if r.get("reproduction") == "reproduced")
    out.append(f"**Paper coverage:** {n_repro} / {len(results)} papers reproduced "
               f"across every planned analysis. `ok` below means only that the "
               f"declared checks passed.\n")
    out.append("| Paper | Status | Checks | Reproduction | Run dir |")
    out.append("|---|---|---|---|---|")
    for r in results:
        st = r.get("status", "?")
        icon = {"ok": "✓", "partial": "△", "fail": "✗", "unreviewed": "⊘",
                  "no_run_found": "—", "no_expected_json": "?"}.get(st, "?")
        cks = f"{r.get('n_passed', 0)}/{r.get('n_total', 0)}"
        if r.get("n_unconfirmed"):
            cks += f" (+{r['n_unconfirmed']} unconfirmed)"
        rd = r.get("run_dir", "—")
        rep = r.get("reproduction", "?")
        if r.get("coverage"):
            rep += f" {r['coverage']}"
        out.append(f"| `{r['paper']}` | {icon} {st} | {cks} | {rep} | `{rd}` |")
    out.append("")
    for r in results:
        out.append(f"## {r['paper']}\n")
        out.append(f"- Status: **{r.get('status', '?')}** · reproduction: "
                   f"**{r.get('reproduction', '?')}**"
                   + (f" ({r['coverage']} analyses)" if r.get("coverage") else ""))
        for a in r.get("analyses", []):
            extra = f", blocker: {a['blocker']}" if a.get("blocker") else ""
            out.append(f"  - analysis `{a['id']}` — {a['state']}{extra}"
                       + (f" — {a['title']}" if a.get("title") else ""))
        if r.get("run_dir"):
            out.append(f"- Run dir: `{r['run_dir']}`")
        for c in r.get("checks", []):
            mark = "✓" if c["passed"] else "✗"
            out.append(f"  - {mark} **{c['name']}** ({c['type']}, class "
                       f"{c.get('class', '?')}): {c['detail']}")
        for c in r.get("unconfirmed", []):
            prov = c.get("provenance") or {}
            # `provenance` is normally {"kind": ..., "quote": ...}, but a
            # plain string is the obvious thing to write and used to crash
            # the whole report here. Accept both.
            if isinstance(prov, str):
                prov = {"kind": prov}
            elif not isinstance(prov, dict):
                prov = {"kind": str(prov)}
            out.append(f"  - ⊘ **{c['name']}** ({c.get('type')}): NOT SCORED — "
                        f"unconfirmed, source: {prov.get('kind', '?')}")
            if prov.get("quote"):
                out.append(f"      > {prov['quote'][:220]}")
        out.append("")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--benchmark", help="One paper-id to score (e.g. waters2024_bap1).")
    g.add_argument("--all", action="store_true", help="Score every benchmark.")
    args = p.parse_args(argv)

    if args.all:
        papers = sorted(d for d in BENCHMARKS.iterdir()
                         if d.is_dir() and (d / "expected.json").is_file())
    else:
        target = BENCHMARKS / args.benchmark
        if not (target / "expected.json").is_file():
            raise SystemExit(f"No Benchmarks/{args.benchmark}/expected.json")
        papers = [target]

    results = [score_benchmark(d) for d in papers]
    RESULTS.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    json_out = RESULTS / f"{ts}_concordance.json"
    md_out   = RESULTS / f"{ts}_concordance.md"
    json_out.write_text(json.dumps(results, indent=2, default=str))
    md_out.write_text(render_markdown(results, ts))

    # Console table
    print(f"\n=== Concordance — {ts} ===")
    for r in results:
        st = r.get("status", "?")
        cks = f"{r.get('n_passed', 0)}/{r.get('n_total', 0)}"
        rep = r.get("reproduction", "?")
        if r.get("coverage"):
            rep += f" {r['coverage']}"
        print(f"  [{st:>12s}] {r['paper']:35s}  {cks:>7s}  reproduction: {rep}")
    print(f"\nSummary: {json_out}")
    print(f"Markdown: {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
