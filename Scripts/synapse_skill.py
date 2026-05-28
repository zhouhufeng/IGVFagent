"""Synapse (Sage Bionetworks) data-discovery + download skill.

A clean-room IGVFagent skill that hits the public Synapse REST API
(https://rest-docs.synapse.org) over pure urllib + json — no dependency
on the upstream ``synapseclient`` Python library.

Covers:
* anonymous-read of public entity metadata, annotations, and
  child-listing for projects / folders (no auth required);
* authenticated download of file payloads via the standard Personal
  Access Token (PAT) carried in ``SYNAPSE_AUTH_TOKEN`` env var —
  required for controlled-access deposits like PsychENCODE +
  IGVF-controlled Synapse folders.

Why a Synapse skill matters for IGVFagent:
* IGVF deposits several controlled-access cohorts to Synapse rather
  than GEO / ENCODE Portal.
* PsychENCODE Consortium (Deng 2024 cortex lentiMPRA, AMP-AD,
  AMP-PD, BrainSpan v2) all distribute via Synapse.
* AMP-PD, ROSMAP, MSBB and similar AMP-AD studies live on Synapse and
  are frequent downstream targets for IGVF cross-reference.

CLI surface::

    igvfagent synapse entity     --syn synXXXXX [--auth-token / env]
    igvfagent synapse children   --syn synXXXXX [--types folder,file,table]
    igvfagent synapse search     --query <q> [--limit N]
    igvfagent synapse download   --syn synXXXXX [--out-dir DIR]
    igvfagent synapse walk       --syn synXXXXX [--max-depth N] [--max-children N]
    igvfagent synapse write-playbook
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Paths + logging
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DOC_DIR = ROOT / "Docs/Synapse"
MANIFEST_DIR = ROOT / "Data/Manifests/Synapse"
DOWNLOAD_DIR = ROOT / "Data/Synapse"
LOG_DIR = ROOT / "Docs/Logs"
PLAYBOOK_PATH = ROOT / "Docs/Skills/SYNAPSE_RETRIEVAL_SKILLS.md"

# Synapse REST endpoints (public docs at https://rest-docs.synapse.org).
SYNAPSE_REPO = "https://repo-prod.prod.sagebase.org/repo/v1"
SYNAPSE_FILE = "https://repo-prod.prod.sagebase.org/file/v1"
USER_AGENT = "IGVFagent/0.1 synapse-skill (+https://github.com/zhouhufeng/IGVFagent)"


# ---------------------------------------------------------------------------
# Logging + filesystem helpers
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = LOG_DIR / f"synapse_skill_{ts}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(logfile), logging.StreamHandler()],
    )
    logging.info("Log file: %s", logfile)
    return logfile


def mkdirs() -> None:
    for d in (DOC_DIR, MANIFEST_DIR, DOWNLOAD_DIR):
        d.mkdir(parents=True, exist_ok=True)


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("._-")
    return s or "synapse"


# ---------------------------------------------------------------------------
# Synapse REST primitives
# ---------------------------------------------------------------------------

def _auth_headers(token: Optional[str] = None) -> "dict[str, str]":
    """Build request headers, opportunistically including PAT from env."""
    h = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    tok = token or os.environ.get("SYNAPSE_AUTH_TOKEN") \
        or os.environ.get("SYNAPSE_PAT")
    if tok:
        h["Authorization"] = f"Bearer {tok.strip()}"
    return h


def synapse_get(path: str, *, base: str = SYNAPSE_REPO,
                  token: Optional[str] = None,
                  params: Optional[dict] = None,
                  timeout: int = 60) -> tuple[int, Any]:
    """GET <base>/<path> with optional params. Returns (status, json|text)."""
    url = base + (path if path.startswith("/") else "/" + path)
    if params:
        qs = urllib.parse.urlencode(
            {k: v for k, v in params.items() if v is not None}, doseq=True
        )
        url = f"{url}?{qs}"
    req = urllib.request.Request(url, headers=_auth_headers(token))
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", "replace")
            try:
                return resp.status, json.loads(raw)
            except json.JSONDecodeError:
                return resp.status, raw
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            payload = body
        return e.code, payload


def synapse_post(path: str, body: dict, *, base: str = SYNAPSE_REPO,
                   token: Optional[str] = None,
                   timeout: int = 60) -> tuple[int, Any]:
    url = base + (path if path.startswith("/") else "/" + path)
    data = json.dumps(body).encode("utf-8")
    headers = _auth_headers(token)
    headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", "replace")
            try:
                return resp.status, json.loads(raw)
            except json.JSONDecodeError:
                return resp.status, raw
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            payload = body
        return e.code, payload


def _is_synapse_id(s: str) -> bool:
    return bool(re.match(r"^syn\d+$", s.strip()))


def _short_type(t: str) -> str:
    """Convert Synapse concreteType -> short tag (folder, file, project, table)."""
    if not t:
        return ""
    return t.split(".")[-1].lower().replace("entity", "")


# ---------------------------------------------------------------------------
# Higher-level operations
# ---------------------------------------------------------------------------

def get_entity(syn_id: str, *, token: Optional[str] = None) -> dict:
    """Fetch the core entity metadata for synXXXXX."""
    if not _is_synapse_id(syn_id):
        raise SystemExit(f"Invalid Synapse ID {syn_id!r} (expect 'synNNNN').")
    sc, data = synapse_get(f"/entity/{syn_id}", token=token)
    if sc == 401 or sc == 403:
        raise SystemExit(
            f"Synapse entity {syn_id} is access-controlled (HTTP {sc}). "
            "Set SYNAPSE_AUTH_TOKEN to your Personal Access Token "
            "after accepting the relevant terms-of-use on synapse.org."
        )
    if sc >= 400 or not isinstance(data, dict):
        raise SystemExit(f"Failed to fetch {syn_id}: HTTP {sc} — {data!r}")
    return data


def get_annotations(syn_id: str, *, token: Optional[str] = None) -> dict:
    """Fetch the user-defined annotations for an entity."""
    sc, data = synapse_get(f"/entity/{syn_id}/annotations2", token=token)
    if sc >= 400:
        return {}
    return data if isinstance(data, dict) else {}


def list_children(syn_id: str, *,
                    types: Optional[list[str]] = None,
                    sort_by: str = "NAME",
                    sort_dir: str = "ASC",
                    max_pages: int = 50,
                    token: Optional[str] = None) -> tuple[int, list[dict]]:
    """Paginated child-listing for a project/folder.

    Returns (totalChildCount, list_of_child_summaries).
    """
    if not _is_synapse_id(syn_id):
        raise SystemExit(f"Invalid Synapse ID {syn_id!r}.")
    types = types or ["folder", "file", "table", "link"]
    body = {
        "parentId":              syn_id,
        "includeTypes":          types,
        "sortBy":                sort_by,
        "sortDirection":         sort_dir,
        "includeTotalChildCount": True,
    }
    children: list[dict] = []
    total = 0
    next_token = None
    pages = 0
    while pages < max_pages:
        if next_token:
            body["nextPageToken"] = next_token
        sc, data = synapse_post("/entity/children", body, token=token)
        if sc == 401 or sc == 403:
            raise SystemExit(
                f"Children listing for {syn_id} is access-controlled "
                f"(HTTP {sc}). Provide SYNAPSE_AUTH_TOKEN."
            )
        if sc >= 400 or not isinstance(data, dict):
            raise SystemExit(
                f"Children listing for {syn_id} failed: HTTP {sc} — {data!r}"
            )
        page = data.get("page") or []
        children.extend(page)
        total = data.get("totalChildCount", total or len(children))
        next_token = data.get("nextPageToken")
        pages += 1
        if not next_token:
            break
        # Be polite — Synapse occasionally rate-limits anonymous bursts.
        time.sleep(0.10)
    return total, children


def walk_entity(syn_id: str, *, max_depth: int = 5,
                  max_children_per_node: int = 200,
                  token: Optional[str] = None) -> list[dict]:
    """Depth-first walk producing one flat list of `{depth, id, name, type, ...}`.

    Hard-caps both depth and per-node children to avoid runaway traversals
    on enormous PsychENCODE-scale projects.
    """
    out: list[dict] = []

    def _walk(node_id: str, depth: int):
        if depth > max_depth:
            return
        try:
            total, children = list_children(
                node_id, max_pages=1 + max_children_per_node // 100,
                token=token,
            )
        except SystemExit as e:
            logging.warning("walk: %s", e)
            return
        for c in children[:max_children_per_node]:
            short_t = _short_type(c.get("type") or "")
            out.append({
                "depth":              depth,
                "id":                 c.get("id"),
                "name":                c.get("name"),
                "type":                short_t,
                "version":             c.get("versionNumber"),
                "modifiedOn":          c.get("modifiedOn"),
                "benefactorId":        c.get("benefactorId"),
                "parent":              node_id,
            })
            if short_t in ("folder", "project"):
                _walk(c.get("id"), depth + 1)

    root = get_entity(syn_id, token=token)
    out.append({
        "depth":             0,
        "id":                root.get("id"),
        "name":              root.get("name"),
        "type":              _short_type(root.get("concreteType") or ""),
        "version":           root.get("versionNumber"),
        "modifiedOn":        root.get("modifiedOn"),
        "benefactorId":      root.get("benefactorId"),
        "parent":            root.get("parentId"),
    })
    _walk(syn_id, 1)
    return out


def search_synapse(query: str, *, limit: int = 25,
                     token: Optional[str] = None) -> dict:
    """Full-text search across the public Synapse catalogue.

    POST /search expects a `SearchQuery` envelope (see rest-docs.synapse.org).
    """
    body = {
        "queryTerm":     [query],
        "size":          int(limit),
        "returnFields":  [
            "Id", "name", "nodeType", "owner", "version_label",
            "consortium", "concreteType", "modified_on",
        ],
    }
    sc, data = synapse_post("/search", body, token=token, timeout=90)
    if sc >= 400:
        raise SystemExit(f"Synapse search failed: HTTP {sc} — {data!r}")
    return data if isinstance(data, dict) else {}


def file_download_url(syn_id: str, *, token: Optional[str] = None,
                         redirect: bool = True) -> str:
    """Return a fully-resolved download URL for a file entity.

    For controlled-access entities the underlying redirect URL is signed
    against the caller's PAT (so you can stream it without leaking creds).
    """
    sc, data = synapse_get(
        f"/entity/{syn_id}/file",
        token=token,
        params={"redirect": "false"},   # we want the URL, not the body
    )
    # The endpoint returns the redirect URL as a *string* (sometimes wrapped
    # in JSON, depending on Accept handling). Handle both.
    if isinstance(data, str):
        url = data.strip().strip('"')
        if url.startswith("http"):
            return url
    if isinstance(data, dict) and "downloadUrl" in data:
        return data["downloadUrl"]
    raise SystemExit(
        f"Could not derive download URL for {syn_id}: HTTP {sc} — {data!r}"
    )


def stream_to_disk(url: str, dest: Path, *, token: Optional[str] = None,
                     timeout: int = 600) -> int:
    """Stream a file to `dest`. Returns bytes written.

    For Synapse-signed URLs (returned by file_download_url) you do NOT need
    to re-send the Bearer token — the URL is already pre-signed.
    """
    headers = {"User-Agent": USER_AGENT}
    # Only re-send PAT if the URL is on the Synapse domain (rare; the file
    # endpoint usually pre-signs to S3).
    if "synapse.org" in url or "sagebase.org" in url:
        h = _auth_headers(token)
        h["User-Agent"] = USER_AGENT
        headers = h
    req = urllib.request.Request(url, headers=headers)
    written = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(req, timeout=timeout) as resp, dest.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            written += len(chunk)
    return written


# ---------------------------------------------------------------------------
# Output rendering
# ---------------------------------------------------------------------------

def write_csv(path: Path, rows: list[dict], cols: Optional[list[str]] = None):
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data: Any):
    path.write_text(json.dumps(data, indent=2, default=str))


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_entity(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    syn = args.syn
    data = get_entity(syn, token=args.auth_token)
    ann = get_annotations(syn, token=args.auth_token) if args.annotations else {}

    out_json = DOC_DIR / f"{timestamp()}_{safe_label(syn)}_entity.json"
    write_json(out_json, {"entity": data, "annotations": ann})

    # Short on-screen summary
    short_t = _short_type(data.get("concreteType") or "")
    print(f"Synapse entity:  {syn}")
    print(f"  name:          {data.get('name')}")
    print(f"  type:          {short_t}")
    print(f"  parentId:      {data.get('parentId')}")
    print(f"  benefactorId:  {data.get('benefactorId')}  (access-control root)")
    print(f"  modifiedOn:    {data.get('modifiedOn')}")
    if ann:
        ann_kv = (ann or {}).get("annotations") or {}
        if ann_kv:
            print("  annotations:")
            for k, v in list(ann_kv.items())[:10]:
                val = v.get("value") if isinstance(v, dict) else v
                print(f"    {k} = {val}")
    print(f"Metadata JSON:   {out_json}")
    return out_json


def cmd_children(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    syn = args.syn
    types = [t.strip() for t in (args.types or "folder,file,table,link").split(",")]
    total, children = list_children(syn, types=types, token=args.auth_token)
    rows = []
    for c in children:
        rows.append({
            "id":               c.get("id"),
            "name":             c.get("name"),
            "type":             _short_type(c.get("type") or ""),
            "version":          c.get("versionNumber"),
            "modifiedOn":       c.get("modifiedOn"),
            "benefactorId":     c.get("benefactorId"),
        })
    label = safe_label(args.label or syn)
    out_csv = MANIFEST_DIR / f"{timestamp()}_{label}_children.csv"
    write_csv(out_csv, rows,
                cols=["id", "name", "type", "version",
                      "modifiedOn", "benefactorId"])
    print(f"Parent: {syn}")
    print(f"Total children: {total}")
    print(f"Returned: {len(rows)}")
    print(f"Manifest CSV: {out_csv}")
    return out_csv


def cmd_walk(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    rows = walk_entity(
        args.syn, max_depth=args.max_depth,
        max_children_per_node=args.max_children,
        token=args.auth_token,
    )
    label = safe_label(args.label or args.syn)
    out_csv = MANIFEST_DIR / f"{timestamp()}_{label}_walk.csv"
    write_csv(out_csv, rows,
                cols=["depth", "id", "name", "type",
                      "version", "parent", "benefactorId", "modifiedOn"])
    print(f"Root: {args.syn}")
    print(f"Nodes (incl. root): {len(rows)}")
    # Per-type counts
    from collections import Counter
    cnts = Counter(r["type"] for r in rows)
    for t, n in cnts.most_common():
        print(f"  {n:5d}  {t}")
    print(f"Manifest CSV: {out_csv}")
    return out_csv


def cmd_search(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    data = search_synapse(args.query, limit=args.limit, token=args.auth_token)
    hits = data.get("hits") or []
    rows = []
    for h in hits:
        rows.append({
            "id":           h.get("id"),
            "name":         h.get("name"),
            "type":         h.get("node_type") or h.get("nodeType"),
            "owner":        h.get("owner_id") or h.get("owner"),
            "modified_on":  h.get("modified_on"),
        })
    label = safe_label(args.label or args.query)
    out_csv = MANIFEST_DIR / f"{timestamp()}_{label}_search.csv"
    write_csv(out_csv, rows,
                cols=["id", "name", "type", "owner", "modified_on"])
    print(f"Query: {args.query!r}")
    print(f"Total found (reported): {data.get('found')}")
    print(f"Hits returned: {len(rows)}")
    for r in rows[:10]:
        print(f"  {r['id']}  {r['name']!s:60.60s}  ({r['type']})")
    print(f"Manifest CSV: {out_csv}")
    return out_csv


def cmd_download(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    syn = args.syn
    meta = get_entity(syn, token=args.auth_token)
    if _short_type(meta.get("concreteType") or "") != "file":
        raise SystemExit(
            f"{syn} is a {_short_type(meta.get('concreteType') or '')!r}, "
            "not a file. Use `synapse children` or `walk` to enumerate "
            "downloadable file IDs first."
        )
    name = meta.get("name") or f"{syn}.bin"
    out_dir = Path(args.out_dir or DOWNLOAD_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / safe_label(name)
    url = file_download_url(syn, token=args.auth_token)
    logging.info("Downloading %s -> %s", syn, dest)
    n = stream_to_disk(url, dest, token=args.auth_token)
    print(f"Downloaded {n:,} bytes  →  {dest}")
    return dest


def cmd_write_playbook(_a: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_MARKDOWN)
    print(f"Wrote playbook: {PLAYBOOK_PATH}")
    return PLAYBOOK_PATH


PLAYBOOK_MARKDOWN = """\
# Synapse retrieval skill

Clean-room IGVFagent integration with Sage Bionetworks Synapse
(`https://www.synapse.org`). Pure urllib + json, no `synapseclient`
Python dep required.

## Subcommands

| Command | Purpose |
|---|---|
| `igvfagent synapse entity --syn synXXX` | Fetch entity metadata + annotations. |
| `igvfagent synapse children --syn synXXX` | Paginated child-listing of a project / folder. |
| `igvfagent synapse walk --syn synXXX --max-depth 3` | Recursive tree walk (capped). |
| `igvfagent synapse search --query "<text>"` | Full-text search across public Synapse. |
| `igvfagent synapse download --syn synXXX --out-dir <dir>` | Stream a file entity to disk. |

## Authentication

* For public entities (most IGVF + most upstream-genomics deposits)
  no auth is required — the REST API answers anonymously.
* For **controlled-access** entities (PsychENCODE, AMP-AD, AMP-PD,
  several IGVF-controlled cohorts) you need a **Personal Access
  Token** (PAT). Create one at
  https://www.synapse.org/#!PersonalAccessTokens after accepting the
  relevant Data-Use-Agreement and grant the `view` and `download`
  scopes.
* Export the token in your shell:

      export SYNAPSE_AUTH_TOKEN="eyJ0eXAiOiJKV1Qi..."

  Every `igvfagent synapse ...` call will then automatically include
  `Authorization: Bearer <token>` headers.

## Why this skill

Several IGVF-relevant repositories live on Synapse rather than GEO /
ENCODE Portal:

* **PsychENCODE Consortium** — Deng 2024 *Science* cortex lentiMPRA
  (`syn21392931` — folder name `NeuREs`); cross-references with the
  IGVF cortex multiome lines.
* **AMP-AD + ROSMAP + MSBB** — Alzheimer's cohorts heavily cross-linked
  to IGVF brain assays.
* **AMP-PD** — Parkinson's; relevant to IGVF Corces/Gladstone
  multiome AnalysisSets.
* **BrainSpan v2** — developmental references for the IGVF brain map.
* Several **IGVF-controlled** AnalysisSets (donor-consent restricted)
  are mirrored to Synapse alongside the IGVF Portal.

## Example: Deng 2024 cortex lentiMPRA

    igvfagent synapse entity --syn syn21392931            # NeuREs root folder
    igvfagent synapse children --syn syn21392931          # Data + PEC sub-folders
    igvfagent synapse walk --syn syn21392931 --max-depth 3
        # Enumerate the full PsychENCODE NeuREs project (cap depth=3).

    # After accepting the PsychENCODE Data-Use Agreement:
    export SYNAPSE_AUTH_TOKEN="eyJ0eXAiOi..."
    igvfagent synapse download --syn <file-syn-id> \\
        --out-dir Data/Benchmarks/deng2024_cortex_mpra/
"""


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------

def _add_common(p: argparse.ArgumentParser):
    p.add_argument(
        "--auth-token", default=None,
        help="Personal Access Token (PAT). Defaults to $SYNAPSE_AUTH_TOKEN "
              "or $SYNAPSE_PAT.",
    )


def main() -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent synapse",
        description="Sage Bionetworks Synapse retrieval + walk + download "
                      "(clean-room, urllib-only).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("entity", help="Fetch entity metadata.")
    sp.add_argument("--syn", required=True,
                      help="Synapse ID (e.g. syn21392931).")
    sp.add_argument("--annotations", action="store_true",
                      help="Also fetch user-defined annotations.")
    _add_common(sp)
    sp.set_defaults(func=cmd_entity)

    sp = sub.add_parser("children", help="List children of a project/folder.")
    sp.add_argument("--syn", required=True)
    sp.add_argument("--types", default="folder,file,table,link",
                      help="Comma-separated entity types to include.")
    sp.add_argument("--label", default=None)
    _add_common(sp)
    sp.set_defaults(func=cmd_children)

    sp = sub.add_parser("walk", help="Recursive walk (depth-capped).")
    sp.add_argument("--syn", required=True)
    sp.add_argument("--max-depth", type=int, default=3)
    sp.add_argument("--max-children", type=int, default=200,
                      help="Max children fetched per node.")
    sp.add_argument("--label", default=None)
    _add_common(sp)
    sp.set_defaults(func=cmd_walk)

    sp = sub.add_parser("search", help="Full-text search across Synapse.")
    sp.add_argument("--query", required=True)
    sp.add_argument("--limit", type=int, default=25)
    sp.add_argument("--label", default=None)
    _add_common(sp)
    sp.set_defaults(func=cmd_search)

    sp = sub.add_parser("download", help="Download a file entity to disk.")
    sp.add_argument("--syn", required=True)
    sp.add_argument("--out-dir", default=None)
    _add_common(sp)
    sp.set_defaults(func=cmd_download)

    sp = sub.add_parser("write-playbook",
                          help="Emit Docs/Skills/SYNAPSE_RETRIEVAL_SKILLS.md.")
    sp.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    return int(args.func(args) is None)


if __name__ == "__main__":
    sys.exit(main())
