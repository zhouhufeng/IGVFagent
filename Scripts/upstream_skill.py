#!/usr/bin/env python3
"""Track the upstream projects this repo absorbed, and notice when they move.

Dozens of skills here reimplement or wrap a published method whose reference
implementation lives in somebody else's repository. Until now that provenance
existed only as prose in a module docstring: which project, sometimes which
licence, never which *version*. So there was no way to answer either of the
two questions that actually matter.

**What did we build against?** "Clean-room reimplementation of the ABC model"
names a method, not a revision. If upstream changed its scoring a year ago,
nothing here records whether we match the old behaviour or the new one.

**Has it changed since?** Nobody was watching. A method paper gets a v2, a
pipeline fixes a bug in its normalisation, a licence changes from MIT to
AGPL — all invisible.

``Docs/upstream.json`` is the manifest: one entry per upstream, with the
commit or tag we built against, the licence, and the RELATIONSHIP — whether we
reimplemented it clean-room, vendored it, or shell out to its binary. That
last field is not bookkeeping: an AGPL upstream we reimplemented without
reading is a different legal position from one we vendored, and the manifest
is where that distinction is written down rather than remembered.

    igvfagent upstream scan      # discover declarations, refresh the manifest
    igvfagent upstream list      # what we absorbed, and from which revision
    igvfagent upstream pin ABC   # record the revision we are current with
    igvfagent upstream check     # what moved upstream since we pinned it

``check`` is the point of the exercise. It asks GitHub for each project's
current head and latest release and reports the drift, so "upstream changed"
becomes something you find out on purpose rather than by surprise.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _localstore as ls                                    # noqa: E402

ROOT = ls.ROOT
MANIFEST = ROOT / "Docs" / "upstream.json"
SCRIPTS = ROOT / "Scripts"
log = logging.getLogger("upstream")

GITHUB_API = "https://api.github.com"

# How this repo relates to the upstream. The distinction is legal as well as
# technical, so it is recorded rather than inferred.
RELATIONSHIPS = (
    "clean-room",      # method reimplemented without reading their source
    "port",            # their algorithm translated, source consulted
    "vendored",        # their code copied in
    "wraps-binary",    # we shell out to their installed program
    "wraps-sdk",       # we call their published library / hosted API as a
                       # dependency; none of their code is in this repo
    "reference",       # consulted for data formats or conventions only
    "unclassified",
)

_REPO_RE = re.compile(r"github\.com/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)")
# No licence is inferred from docstring text. A first attempt matched the
# first SPDX-looking token in each module header and called it the upstream's
# licence -- but these headers state IGVFagent's OWN licence first, under a
# "License posture" heading, with the upstream's named separately. That read
# chipatlas, corneto and tabula-sapiens as mislabelled when all three were
# meticulously correct, and would have had me "fix" three accurate docstrings.
# GitHub is authoritative for what the upstream is licensed under; ask it.


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")


def load() -> dict:
    try:
        return json.loads(MANIFEST.read_text())
    except (OSError, ValueError):
        return {"upstreams": []}


def save(data: dict) -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    data["upstreams"].sort(key=lambda e: e["repo"].lower())
    MANIFEST.write_text(json.dumps(data, indent=2) + "\n")


def _cli_names() -> "dict[str, str]":
    """module filename -> CLI subcommand, from the dispatcher."""
    src = (SCRIPTS / "cli.py").read_text()
    out = {}
    for m in re.finditer(r'"([a-z0-9-]+)":\s*\(\s*"igvfagent\.([a-z0-9_]+)"', src):
        out[m.group(2) + ".py"] = m.group(1)
    return out


def discover() -> "dict[str, dict]":
    """Upstream declarations found in module docstrings.

    Reads only the first 4 kB of each file: a repository named that far down is
    a passing mention, not a statement of provenance, and treating every URL in
    a long module as an upstream would fill the manifest with noise.
    """
    cli = _cli_names()
    found: "dict[str, dict]" = {}
    for path in sorted(SCRIPTS.glob("*.py")):
        head = path.read_text(errors="replace")[:4000]
        for m in _REPO_RE.finditer(head):
            repo = m.group(1).rstrip(".,);:")
            if repo.lower().startswith("zhouhufeng/"):
                continue
            e = found.setdefault(repo, {"repo": repo, "modules": [],
                                        "skills": []})
            if path.name not in e["modules"]:
                e["modules"].append(path.name)
            skill = cli.get(path.name)
            if skill and skill not in e["skills"]:
                e["skills"].append(skill)
    return found


def _gh(path: str) -> "dict | list | None":
    req = urllib.request.Request(
        GITHUB_API + path,
        headers={"Accept": "application/vnd.github+json",
                 "User-Agent": "igvfagent-upstream"})
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 403:
            log.warning("GitHub rate limit or access denied for %s — set "
                        "GITHUB_TOKEN to raise the limit", path)
        elif e.code == 404:
            log.warning("not found upstream: %s", path)
        return None
    except Exception as exc:
        log.warning("%s: %s", path, exc)
        return None


# --------------------------------------------------------------------------


def cmd_scan(args) -> int:
    data = load()
    by_repo = {e["repo"]: e for e in data["upstreams"]}
    found = discover()
    added = updated = 0
    for repo, info in found.items():
        e = by_repo.get(repo)
        if e is None:
            # A new entry is deliberately UNPINNED and unclassified. Guessing
            # either would put a claim in the manifest that nobody checked,
            # which is the failure this file exists to prevent.
            e = {"repo": repo, "skills": [], "modules": [],
                 "license": "",
                 "relationship": "unclassified",
                 "pinned": None, "notes": ""}
            data["upstreams"].append(e)
            by_repo[repo] = e
            added += 1
        before = json.dumps(e, sort_keys=True)
        e["modules"] = sorted(set(e.get("modules", [])) | set(info["modules"]))
        e["skills"] = sorted(set(e.get("skills", [])) | set(info["skills"]))

        if json.dumps(e, sort_keys=True) != before:
            updated += 1
    save(data)
    print(f"{len(found)} upstream(s) declared in docstrings")
    print(f"  {added} new, {updated} updated -> {MANIFEST.relative_to(ROOT)}")
    unpinned = [e for e in data["upstreams"] if not e.get("pinned")]
    uncls = [e for e in data["upstreams"] if e.get("relationship") == "unclassified"]
    if unpinned:
        print(f"\n  {len(unpinned)} not pinned to a revision — "
              f"`igvfagent upstream pin <repo>`")
    if uncls:
        print(f"  {len(uncls)} unclassified — set relationship "
              f"({', '.join(RELATIONSHIPS[:-1])}) in the manifest")
    return 0


def cmd_list(args) -> int:
    data = load()
    rows = data["upstreams"]
    if args.unpinned:
        rows = [e for e in rows if not e.get("pinned")]
    if not rows:
        print("Nothing recorded. Run: igvfagent upstream scan")
        return 0
    print(f"{len(rows)} upstream project(s)\n")
    for e in rows:
        pin = e.get("pinned") or {}
        ref = (pin.get("ref") or "")[:12] or "NOT PINNED"
        print(f"  {e['repo']:46} {e.get('license','?'):10} "
              f"{e.get('relationship','?'):14} {ref}")
        if e.get("skills"):
            print(f"    skills: {', '.join(e['skills'])}")
        if args.verbose and e.get("modules"):
            print(f"    modules: {', '.join(e['modules'])}")
        if args.verbose and e.get("notes"):
            print(f"    {e['notes']}")
    return 0


def cmd_pin(args) -> int:
    data = load()
    matches = [e for e in data["upstreams"]
               if args.repo.lower() in e["repo"].lower()
               or args.repo.lower() in [s.lower() for s in e.get("skills", [])]]
    if not matches:
        print(f"No upstream matches {args.repo!r}. `igvfagent upstream list`",
              file=sys.stderr)
        return 1
    if len(matches) > 1 and not args.all:
        print(f"{args.repo!r} matches {len(matches)}: "
              + ", ".join(e["repo"] for e in matches)
              + "\nBe more specific, or pass --all.", file=sys.stderr)
        return 1
    for e in matches:
        info = _gh(f"/repos/{e['repo']}")
        if not info:
            print(f"  {e['repo']}: could not reach GitHub")
            continue
        branch = info.get("default_branch", "main")
        commit = _gh(f"/repos/{e['repo']}/commits/{branch}") or {}
        sha = commit.get("sha", "")
        when = ((commit.get("commit") or {}).get("committer") or {}).get("date", "")
        rel = _gh(f"/repos/{e['repo']}/releases/latest") or {}
        e["pinned"] = {
            "ref": sha, "kind": "commit", "branch": branch,
            "upstream_date": when,
            "release": rel.get("tag_name") or "",
            "pinned_at": time.strftime("%Y-%m-%d"),
        }
        if info.get("license") and not e.get("license"):
            e["license"] = (info["license"] or {}).get("spdx_id", "")
        print(f"  pinned {e['repo']} @ {sha[:12]} ({branch}, {when[:10]})"
              + (f", release {rel.get('tag_name')}" if rel.get("tag_name") else ""))
    save(data)
    return 0


def cmd_check(args) -> int:
    data = load()
    rows = [e for e in data["upstreams"] if e.get("pinned")]
    skipped = [e for e in data["upstreams"] if not e.get("pinned")]
    if not rows:
        print("Nothing is pinned yet, so there is nothing to compare against.")
        print("  igvfagent upstream pin <repo>")
        return 0

    moved, still, gone = [], [], []
    print(f"checking {len(rows)} pinned upstream(s)…\n")
    for e in rows:
        pin = e["pinned"]
        branch = pin.get("branch") or "main"
        commit = _gh(f"/repos/{e['repo']}/commits/{branch}")
        if commit is None:
            gone.append(e)
            continue
        head = commit.get("sha", "")
        if head == pin.get("ref"):
            still.append(e)
            continue
        cmp_ = _gh(f"/repos/{e['repo']}/compare/{pin['ref']}...{head}") or {}
        rel = _gh(f"/repos/{e['repo']}/releases/latest") or {}
        moved.append({
            "entry": e, "head": head,
            "ahead": cmp_.get("ahead_by"),
            "date": ((commit.get("commit") or {}).get("committer") or {}).get("date", ""),
            "release": rel.get("tag_name") or "",
            "new_release": bool(rel.get("tag_name")
                                and rel.get("tag_name") != pin.get("release")),
        })

    if moved:
        print(f"MOVED SINCE WE PINNED ({len(moved)}):\n")
        for m in moved:
            e = m["entry"]
            ahead = (f"{m['ahead']} commit(s) ahead" if m["ahead"] is not None
                     else "ahead by an unknown amount")
            print(f"  {e['repo']}")
            print(f"    pinned {e['pinned']['ref'][:12]} "
                  f"({e['pinned'].get('upstream_date','')[:10]}) -> "
                  f"head {m['head'][:12]} ({m['date'][:10]})  {ahead}")
            if m["new_release"]:
                print(f"    NEW RELEASE: {m['release']} "
                      f"(we pinned at {e['pinned'].get('release') or 'no release'})")
            print(f"    used by: {', '.join(e.get('skills') or e.get('modules', []))}"
                  f"   [{e.get('relationship','?')}]")
            print(f"    https://github.com/{e['repo']}/compare/"
                  f"{e['pinned']['ref'][:12]}...{m['head'][:12]}")
            print()
    if still:
        print(f"UNCHANGED ({len(still)}): "
              + ", ".join(e["repo"] for e in still))
    if gone:
        print(f"\nUNREACHABLE ({len(gone)}): "
              + ", ".join(e["repo"] for e in gone)
              + "\n  renamed, made private, or deleted — worth knowing either way")
    if skipped:
        print(f"\nNOT PINNED, so not checked ({len(skipped)}): "
              + ", ".join(e["repo"] for e in skipped))

    if args.json:
        print()
        print(json.dumps({"moved": [{"repo": m["entry"]["repo"],
                                      "ahead": m["ahead"],
                                      "head": m["head"],
                                      "new_release": m["new_release"]}
                                     for m in moved]}, indent=2))
    ls.record_analysis("upstream", subcommand="check",
                       label=f"{len(moved)} moved",
                       outputs=[str(MANIFEST)])
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent upstream",
        description="Which upstream projects this repo absorbed, which "
                    "revision each was built against, and what has changed "
                    "there since.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scan", help="Find upstream declarations in docstrings "
                                    "and refresh the manifest.")
    s.set_defaults(func=cmd_scan)

    s = sub.add_parser("list", help="What we absorbed, and from which revision.")
    s.add_argument("--unpinned", action="store_true",
                   help="Only those with no recorded revision.")
    s.add_argument("--verbose", action="store_true")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("pin", help="Record the revision we are current with.")
    s.add_argument("repo", help="Repo substring or skill name.")
    s.add_argument("--all", action="store_true",
                   help="Pin every match, not just a unique one.")
    s.set_defaults(func=cmd_pin)

    s = sub.add_parser("check", help="What moved upstream since we pinned it.")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_check)
    return p


def main(argv=None) -> int:
    setup_logging()
    args = build_parser().parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
