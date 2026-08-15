#!/usr/bin/env python3
"""Score IGVFagent benchmark runs against per-paper ``expected.json``.

For each benchmark directory under ``Benchmarks/<paper-id>/``:

1. Read ``expected.json`` (ground-truth ranges + qualitative checks).
2. Locate the latest output directory under ``Docs/<skill>/2*_<label>/``
   that matches the benchmark's label.
3. Open the canonical artefact (usually ``summary.json`` or a TSV).
4. Apply each declared check; tally pass / fail / skip.
5. Emit a per-paper JSON + a suite-level Markdown summary under
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


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = Path(__file__).resolve().parent
RESULTS = BENCHMARKS / "results"
DOCS = ROOT / "Docs"


def latest_run_dir(skill_dir_name: str, label: str) -> Path | None:
    """Find the most recent ``Docs/<skill>/2*_<label>*/`` directory.

    Two conventions are supported:

    1. Per-run directory: ``Docs/<skill>/<ts>_<label>/`` containing
       ``summary.json``, TSVs, plots, etc. (the convention used by
       the newer skills: mavedb, multiome, chipatlas, portal, catalog).
    2. Flat-file convention: ``Docs/<skill>/<ts>_<label>_*`` files
       sitting directly under the skill dir (older skills like the
       legacy ``mpra pull``). For these we return the skill dir itself
       and the artefact checks must match by glob pattern.
    """
    base = DOCS / skill_dir_name
    if not base.is_dir():
        return None
    # Convention 1: dir match
    dirs = sorted(
        (p for p in base.glob(f"2*_*{label}*") if p.is_dir()),
        key=lambda p: p.name,
        reverse=True,
    )
    if dirs:
        return dirs[0]
    # Convention 2: flat-file match
    flat = sorted(
        (p for p in base.glob(f"2*_*{label}*") if p.is_file()),
        key=lambda p: p.name,
        reverse=True,
    )
    if flat:
        # Return the skill dir; checks use the label to glob for files.
        return base
    return None


def read_artefact(d: Path, filename: str) -> Any:
    """Read a JSON or TSV artefact from a run directory."""
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
    fname = spec.get("filename") or spec.get("artefact")
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
    run_dir = latest_run_dir(skill, label)
    if run_dir is None:
        return {"paper": paper_dir.name, "status": "no_run_found",
                "skill": skill, "label": label}
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

    # Extra search roots for skills that scatter artefacts. expected.json
    # may declare ``extra_search_dirs: ["Data/Manifests/MPRA", "Data"]`` etc.
    extras = [ROOT / x for x in (spec.get("extra_search_dirs") or [])]

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
            if ctype == "range":
                v = get_path(payload, chk["path"])
                ok, msg = check_range(v, chk)
            elif ctype == "in_set":
                v = get_path(payload, chk["path"])
                ok, msg = check_set_member(v, chk)
            elif ctype == "artefact":
                ok, msg = check_artefact_exists(run_dir, chk, extras)
            elif ctype == "row_count_tsv":
                rows = read_artefact(run_dir, chk["filename"])
                v = len(rows) if isinstance(rows, list) else None
                ok, msg = check_range(v, chk)
            else:
                ok, msg = False, f"unknown check type {ctype!r}"
        except Exception as e:
            ok, msg = False, f"check failed with exception: {e}"
        result["checks"].append({"name": name, "type": ctype,
                                   "passed": ok, "detail": msg})
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
    out.append("| Paper | Status | Checks | Run dir |")
    out.append("|---|---|---|---|")
    for r in results:
        st = r.get("status", "?")
        icon = {"ok": "✓", "partial": "△", "fail": "✗", "unreviewed": "⊘",
                  "no_run_found": "—", "no_expected_json": "?"}.get(st, "?")
        cks = f"{r.get('n_passed', 0)}/{r.get('n_total', 0)}"
        if r.get("n_unconfirmed"):
            cks += f" (+{r['n_unconfirmed']} unconfirmed)"
        rd = r.get("run_dir", "—")
        out.append(f"| `{r['paper']}` | {icon} {st} | {cks} | `{rd}` |")
    out.append("")
    for r in results:
        out.append(f"## {r['paper']}\n")
        out.append(f"- Status: **{r.get('status', '?')}**")
        if r.get("run_dir"):
            out.append(f"- Run dir: `{r['run_dir']}`")
        for c in r.get("checks", []):
            mark = "✓" if c["passed"] else "✗"
            out.append(f"  - {mark} **{c['name']}** ({c['type']}): {c['detail']}")
        for c in r.get("unconfirmed", []):
            prov = c.get("provenance") or {}
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
        print(f"  [{st:>12s}] {r['paper']:35s}  {cks}")
    print(f"\nSummary: {json_out}")
    print(f"Markdown: {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
