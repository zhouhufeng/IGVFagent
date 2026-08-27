"""Read back the artefacts the agent just produced.

Every skill announces its outputs on stdout as ``Report: <path>`` /
``Manifest: <path>``, and the agent faithfully reports those paths — but until
now it could not *open* them. There is no file-reading tool anywhere in the
registry, so the model only ever knew what a skill's stdout summary happened
to print. That is why answers stop at "18 diseases" instead of naming them:
the list exists on disk, in a file the agent wrote, that it cannot read.

This skill closes that loop. Reads are contained by ``_pathguard`` — the same
guard the browser UI uses for inline rendering — so a path can only resolve
inside the workspace, and never to secret material within it
(``Docs/Secret/``, ``.env``, ``*.pem``, credential-like names). Binary files
are refused rather than dumped as mojibake, and output is byte-capped so a
200 MB manifest cannot blow the model's context window in one call.

Subcommands::

    igvfagent artifact read  --path Docs/KGTraversal/<run>/report.md
    igvfagent artifact read  --path <p> --head 40
    igvfagent artifact grep  --pattern BRCA1 --path Docs/KGTraversal
    igvfagent artifact ls    --path Docs/KGTraversal/<run>

Pure standard library.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

try:
    from igvfagent import _pathguard
except Exception:  # pragma: no cover - direct-script execution
    import _pathguard  # type: ignore

__all__ = ["main"]

# Big enough for a full report, small enough that one call cannot evict the
# conversation. Callers wanting more should page with --head / --tail.
_DEFAULT_MAX_BYTES = 200_000
_HARD_MAX_BYTES = 2_000_000


def _root() -> Path:
    return _pathguard.project_root()


def _resolve(path: str, *, allow_dir: bool = False) -> Path:
    """Resolve a user/model-supplied path against the workspace, then guard it.

    ``allow_dir`` is for ls/grep, which legitimately target a run directory;
    containment and the secrets denylist still apply either way.
    """
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = _root() / p
    reason = _pathguard.why_blocked(p, require_file=not allow_dir)
    if reason:
        raise PermissionError(f"refusing to read {path!r}: {reason}")
    return p.resolve()


def _looks_binary(chunk: bytes) -> bool:
    # A NUL in the first block is the classic, cheap binary test; it is what
    # git uses. Avoids trying to decode an .h5ad or a .parquet as text.
    return b"\x00" in chunk


def read_artifact(path: str, *, max_bytes: int = _DEFAULT_MAX_BYTES,
                  head: "int | None" = None,
                  tail: "int | None" = None) -> dict:
    p = _resolve(path)
    size = p.stat().st_size
    max_bytes = max(1, min(int(max_bytes), _HARD_MAX_BYTES))

    with open(p, "rb") as fh:
        probe = fh.read(8192)
        if _looks_binary(probe):
            return {"path": str(p.relative_to(_root())), "bytes": size,
                    "binary": True, "text": None,
                    "note": "binary file — not decoded. Use a skill that "
                            "understands this format."}
        fh.seek(0)
        raw = fh.read(max_bytes + 1)

    truncated = len(raw) > max_bytes
    text = raw[:max_bytes].decode("utf-8", errors="replace")

    lines = text.splitlines()
    if head is not None:
        lines = lines[:max(0, int(head))]
        text = "\n".join(lines)
    elif tail is not None:
        lines = lines[-max(0, int(tail)):]
        text = "\n".join(lines)

    return {
        "path": str(p.relative_to(_root())),
        "bytes": size,
        "binary": False,
        "truncated": truncated,
        "returned_lines": len(lines),
        "text": text,
    }


def grep_artifacts(pattern: str, *, path: str = "Docs",
                   max_hits: int = 50, ignore_case: bool = True) -> dict:
    import re
    base = _resolve(path, allow_dir=True) if path else _root()
    flags = re.IGNORECASE if ignore_case else 0
    rx = re.compile(pattern, flags)
    hits: "list[dict]" = []
    scanned = 0

    targets = [base] if base.is_file() else sorted(
        q for q in base.rglob("*")
        if q.is_file() and _pathguard.is_safe_artifact(q))
    for q in targets:
        if len(hits) >= max_hits:
            break
        try:
            with open(q, "rb") as fh:
                if _looks_binary(fh.read(8192)):
                    continue
            scanned += 1
            with open(q, "r", encoding="utf-8", errors="replace") as fh:
                for i, line in enumerate(fh, 1):
                    if rx.search(line):
                        hits.append({"file": str(q.relative_to(_root())),
                                      "line": i, "text": line.rstrip()[:300]})
                        if len(hits) >= max_hits:
                            break
        except OSError:
            continue
    return {"pattern": pattern, "files_scanned": scanned,
            "hits": hits, "truncated": len(hits) >= max_hits}


def list_artifacts(path: str = "Docs", *, limit: int = 200) -> dict:
    base = _resolve(path, allow_dir=True) if path else _root()
    if base.is_file():
        return {"path": str(base.relative_to(_root())), "entries": [
            {"name": base.name, "bytes": base.stat().st_size, "dir": False}]}
    entries = []
    for q in sorted(base.iterdir())[:max(1, int(limit))]:
        try:
            entries.append({"name": q.name, "dir": q.is_dir(),
                             "bytes": (q.stat().st_size if q.is_file() else None)})
        except OSError:
            continue
    return {"path": str(base.relative_to(_root())), "entries": entries}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="igvfagent artifact",
        description="Read back files the agent produced, contained to the "
                    "workspace.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("read", help="Read a text artefact")
    r.add_argument("--path", required=True)
    r.add_argument("--max-bytes", type=int, default=_DEFAULT_MAX_BYTES)
    r.add_argument("--head", type=int)
    r.add_argument("--tail", type=int)

    g = sub.add_parser("grep", help="Search inside workspace artefacts")
    g.add_argument("--pattern", required=True)
    g.add_argument("--path", default="Docs")
    g.add_argument("--max-hits", type=int, default=50)

    l = sub.add_parser("ls", help="List a directory of artefacts")
    l.add_argument("--path", default="Docs")
    l.add_argument("--limit", type=int, default=200)

    args = parser.parse_args(argv)
    try:
        if args.cmd == "read":
            out = read_artifact(args.path, max_bytes=args.max_bytes,
                                head=args.head, tail=args.tail)
            # Print the text plainly so it lands in the model's tool result
            # as readable content rather than a JSON-escaped blob.
            if out.get("text") is not None:
                hdr = (f"# {out['path']}  ({out['bytes']:,} bytes"
                       + (", TRUNCATED" if out.get("truncated") else "") + ")")
                print(hdr)
                print(out["text"])
            else:
                print(json.dumps(out, indent=2))
        elif args.cmd == "grep":
            print(json.dumps(grep_artifacts(args.pattern, path=args.path,
                                             max_hits=args.max_hits), indent=2))
        elif args.cmd == "ls":
            print(json.dumps(list_artifacts(args.path, limit=args.limit),
                              indent=2))
    except PermissionError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except FileNotFoundError:
        print(f"error: no such file: {getattr(args, 'path', '')}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
