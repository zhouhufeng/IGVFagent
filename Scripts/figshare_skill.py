"""figshare data-discovery + download skill.

A clean-room IGVFagent skill that hits the public figshare REST API
(https://docs.figshare.com — API v2) over pure urllib + json, with no
dependency on any third-party figshare client.

figshare is, alongside Zenodo, one of the two dominant general-purpose
repositories for the raw / derived data deposited alongside research
publications (supplementary tables, per-record statistics, model outputs,
processed matrices). Reproducibility benchmarks in IGVFagent repeatedly need
to pull an author-deposited figshare item given only the DOI, the numeric
article id, or a *private share link* (``figshare.com/s/<token>``) printed in
a paper's Data Availability statement.

Covers:
* resolve an article from a numeric id, a full DOI, an article URL, or a
  private share token / ``/s/<token>`` link;
* list every file in an article (name, size, md5) without downloading;
* download one file, or every file in an article, to disk with md5 verify;
* full-text search across public figshare articles.

CLI surface::

    igvfagent figshare article  --id <id|token|url|doi> [--label L]
    igvfagent figshare files    --id <id|token|url|doi> [--label L]
    igvfagent figshare download --id <id|token|url|doi> [--file-id N] [--out-dir DIR]
    igvfagent figshare search   --query <q> [--limit N]
    igvfagent figshare write-playbook

Auth: public figshare items need no token. A personal token (for private or
embargoed items you own) may be supplied via ``--token`` or the
``FIGSHARE_TOKEN`` env var.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Paths + logging
# ---------------------------------------------------------------------------

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOC_DIR = ROOT / "Docs/Figshare"
MANIFEST_DIR = ROOT / "Data/Manifests/Figshare"
DOWNLOAD_DIR = ROOT / "Data/Figshare"
LOG_DIR = ROOT / "Docs/Logs"
PLAYBOOK_PATH = ROOT / "Docs/Skills/FIGSHARE_RETRIEVAL_SKILLS.md"

API = "https://api.figshare.com/v2"
WEB = "https://figshare.com"
USER_AGENT = "IGVFagent/0.1 figshare-skill (+https://github.com/zhouhufeng/IGVFagent)"


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = LOG_DIR / f"figshare_skill_{ts}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(logfile), logging.StreamHandler()],
    )
    logging.info("Log file: %s", logfile)
    return logfile


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def _headers(token: Optional[str]) -> "dict[str, str]":
    h = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if token:
        h["Authorization"] = f"token {token}"
    return h


def _get_json(url: str, token: Optional[str] = None, timeout: int = 60) -> Any:
    req = urllib.request.Request(url, headers=_headers(token))
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post_json(url: str, payload: dict, token: Optional[str] = None,
               timeout: int = 60) -> Any:
    data = json.dumps(payload).encode("utf-8")
    h = _headers(token)
    h["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=h, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


# ---------------------------------------------------------------------------
# Article resolution — from id / URL / DOI / private share token
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(r"^[0-9a-f]{16,}$", re.I)


def resolve_article(spec: str, token: Optional[str] = None,
                    timeout: int = 60) -> "tuple[int, Optional[str]]":
    """Return (article_id, private_link_token) from a flexible spec:
    numeric id, article URL, DOI, ``/s/<token>`` share URL, or bare token."""
    spec = spec.strip()

    # bare numeric id
    if spec.isdigit():
        return int(spec), None

    # /s/<token> share link (URL or bare token)
    m = re.search(r"/s/([0-9a-f]{16,})", spec, re.I)
    share = m.group(1) if m else (spec if _TOKEN_RE.match(spec) else None)
    if share:
        aid = _resolve_share_token(share, timeout)
        return aid, share

    # title-only URL with no id and no token
    if spec.startswith("http"):
        raise SystemExit(f"No article id found in URL: {spec!r}")

    # /articles/<...>/<id> URL
    m = re.search(r"/articles/(?:[^/]+/)*(\d+)", spec)
    if m:
        return int(m.group(1)), None

    # DOI (…/m9.figshare.<id>[.v<n>])
    m = re.search(r"figshare\.(\d+)", spec)
    if m:
        return int(m.group(1)), None

    raise SystemExit(f"Could not resolve a figshare article from: {spec!r}")


def _resolve_share_token(share: str, timeout: int) -> int:
    """Resolve a private share token to an article id via the web redirect.

    Note: figshare's ``/s/<token>`` viewer is served through Cloudflare and
    typically returns a JS challenge (HTTP 202, no id) to non-browser clients,
    so opaque tokens frequently cannot be resolved programmatically. In that
    case supply the numeric article id (or DOI) explicitly together with
    ``--private-link <token>``; ``api.figshare.com`` honours the token on the
    per-article endpoints and is not Cloudflare-walled."""
    try:
        req = urllib.request.Request(f"{WEB}/s/{share}",
                                     headers={"User-Agent": USER_AGENT,
                                              "Accept": "text/html"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            final = resp.geturl()
            body = resp.read(1 << 18).decode("utf-8", "replace")
        for pat in (r"/articles/(?:[^/]+/)*(\d{5,})",
                    r'"(?:article_id|id)"\s*:\s*(\d{5,})',
                    r"figshare\.(\d{5,})"):
            m = re.search(pat, final) or re.search(pat, body)
            if m:
                return int(m.group(1))
    except urllib.error.HTTPError as e:
        logging.debug("share resolve -> HTTP %s", e.code)
    raise SystemExit(
        f"Could not resolve private share token {share!r} to an article id "
        "(figshare's /s/ viewer is Cloudflare-protected). Re-run with the "
        "numeric article id or DOI plus the token, e.g.:\n"
        f"  igvfagent figshare files --id <ARTICLE_ID> --private-link {share}")


def _article_url(aid: int, share: Optional[str]) -> str:
    u = f"{API}/articles/{aid}"
    return u + (f"?private_link={share}" if share else "")


def _files_url(aid: int, share: Optional[str]) -> str:
    u = f"{API}/articles/{aid}/files"
    return u + (f"?private_link={share}" if share else "")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _fetch_article(spec: str, token: Optional[str], timeout: int,
                   private_link: Optional[str] = None):
    if private_link and str(spec).strip().isdigit():
        aid, share = int(spec), private_link
    else:
        aid, share = resolve_article(spec, token, timeout)
        share = share or private_link
    meta = _get_json(_article_url(aid, share), token, timeout)
    # files sometimes only appear on the dedicated endpoint
    files = meta.get("files")
    if not files:
        try:
            files = _get_json(_files_url(aid, share), token, timeout)
            meta["files"] = files
        except urllib.error.HTTPError:
            meta["files"] = []
    return aid, share, meta


def cmd_article(args: argparse.Namespace) -> None:
    setup_logging()
    aid, share, meta = _fetch_article(args.id, args.token, args.timeout, getattr(args,'private_link',None))
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = DOC_DIR / f"{_ts()}_article_{aid}{('_'+args.label) if args.label else ''}.json"
    out.write_text(json.dumps(meta, indent=2))
    print(f"figshare article:  {aid}")
    print(f"  title:      {meta.get('title','')}")
    print(f"  doi:        {meta.get('doi','')}")
    authors = ", ".join(a.get("full_name", "") for a in meta.get("authors", [])[:6])
    print(f"  authors:    {authors}")
    print(f"  published:  {meta.get('published_date','') or meta.get('created_date','')}")
    print(f"  files:      {len(meta.get('files') or [])}")
    if share:
        print(f"  private_link: {share}")
    print(f"Metadata JSON: {out}")


def cmd_files(args: argparse.Namespace) -> None:
    setup_logging()
    aid, share, meta = _fetch_article(args.id, args.token, args.timeout, getattr(args,'private_link',None))
    files = meta.get("files") or []
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    man = MANIFEST_DIR / f"{_ts()}_article_{aid}{('_'+args.label) if args.label else ''}_files.csv"
    with man.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file_id", "name", "size_bytes", "size_mb",
                    "computed_md5", "download_url"])
        for f in files:
            w.writerow([f.get("id"), f.get("name"), f.get("size"),
                        round((f.get("size") or 0) / 1e6, 3),
                        f.get("computed_md5") or f.get("supplied_md5") or "",
                        f.get("download_url") or ""])
    total = sum((f.get("size") or 0) for f in files)
    print(f"figshare article {aid}: {len(files)} files, {total/1e6:.1f} MB total")
    for f in files:
        print(f"  [{f.get('id')}] {(f.get('size') or 0)/1e6:8.2f} MB  {f.get('name')}")
    print(f"Manifest CSV: {man}")


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_one(f: dict, out_dir: Path, token: Optional[str],
                  share: Optional[str], timeout: int) -> Path:
    url = f.get("download_url")
    if share and "private_link" not in url:
        url = url + (("&" if "?" in url else "?") + f"private_link={share}")
    dest = out_dir / f.get("name")
    req = urllib.request.Request(url, headers=_headers(token))
    with urllib.request.urlopen(req, timeout=timeout) as resp, dest.open("wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    want = f.get("computed_md5") or f.get("supplied_md5")
    if want:
        got = _md5(dest)
        status = "OK" if got == want else f"MISMATCH (got {got})"
        logging.info("md5 %s: %s", dest.name, status)
    logging.info("downloaded %s (%.1f MB)", dest, dest.stat().st_size / 1e6)
    return dest


def cmd_download(args: argparse.Namespace) -> None:
    setup_logging()
    aid, share, meta = _fetch_article(args.id, args.token, args.timeout, getattr(args,'private_link',None))
    files = meta.get("files") or []
    if args.file_id is not None:
        files = [f for f in files if str(f.get("id")) == str(args.file_id)]
        if not files:
            raise SystemExit(f"file-id {args.file_id} not in article {aid}")
    if args.name_contains:
        files = [f for f in files if args.name_contains.lower() in (f.get("name") or "").lower()]
    out_dir = Path(args.out_dir) if args.out_dir else DOWNLOAD_DIR / f"article_{aid}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(files)} file(s) from figshare article {aid} -> {out_dir}")
    for f in files:
        p = _download_one(f, out_dir, args.token, share, args.timeout)
        print(f"  {p.stat().st_size/1e6:8.2f} MB  {p}")


def cmd_search(args: argparse.Namespace) -> None:
    setup_logging()
    res = _post_json(f"{API}/articles/search",
                     {"search_for": args.query, "page_size": args.limit},
                     args.token, args.timeout)
    print(f"figshare search {args.query!r}: {len(res)} results")
    for r in res:
        print(f"  [{r.get('id')}] {r.get('title','')[:80]}  ({r.get('doi','')})")


def cmd_write_playbook(args: argparse.Namespace) -> None:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(__doc__)
    print(f"Wrote {PLAYBOOK_PATH}")


def _add_common(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--token", default=os.environ.get("FIGSHARE_TOKEN"),
                    help="figshare personal token (env FIGSHARE_TOKEN); "
                         "only needed for private/embargoed items.")
    sp.add_argument("--timeout", type=int, default=120)


def main() -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent figshare",
        description="figshare retrieval skill (clean-room, urllib-only).")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("article", help="Fetch article metadata.")
    sp.add_argument("--id", required=True, help="id | /s/<token> | URL | DOI")
    sp.add_argument("--private-link", default=None, help="figshare private share token (use with numeric --id).")
    sp.add_argument("--label", default=None)
    _add_common(sp); sp.set_defaults(func=cmd_article)

    sp = sub.add_parser("files", help="List files in an article.")
    sp.add_argument("--id", required=True)
    sp.add_argument("--private-link", default=None, help="figshare private share token (use with numeric --id).")
    sp.add_argument("--label", default=None)
    _add_common(sp); sp.set_defaults(func=cmd_files)

    sp = sub.add_parser("download", help="Download file(s) from an article.")
    sp.add_argument("--id", required=True)
    sp.add_argument("--private-link", default=None, help="figshare private share token (use with numeric --id).")
    sp.add_argument("--file-id", default=None, help="Download only this file id.")
    sp.add_argument("--name-contains", default=None,
                    help="Download only files whose name contains this string.")
    sp.add_argument("--out-dir", default=None)
    _add_common(sp); sp.set_defaults(func=cmd_download)

    sp = sub.add_parser("search", help="Full-text search public articles.")
    sp.add_argument("--query", required=True)
    sp.add_argument("--limit", type=int, default=25)
    _add_common(sp); sp.set_defaults(func=cmd_search)

    sp = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/FIGSHARE_RETRIEVAL_SKILLS.md.")
    sp.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    return int(args.func(args) is None)


if __name__ == "__main__":
    sys.exit(main())
