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
import csv
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


def _nearest_existing(p: Path) -> "tuple[Path, list[str]]":
    """Walk up to the first directory that exists, and list what is in it.

    A listing of a path that does not exist is a dead end for an agent: it
    has no way to learn the right name, so it either gives up or guesses
    again. This happened on a real run -- the model asked for
    `Docs/BaseEditingScreen/18loci_uptake` when the run had written
    `Docs/BaseEditingScreen/20260911_033012_18loci_uptake`, got exit 2, and
    reported the whole answer as incomplete. The directory it wanted was one
    `ls` of the parent away.
    """
    cur = p
    for _ in range(4):
        cur = cur.parent
        if cur.is_dir():
            try:
                return cur, sorted(q.name for q in cur.iterdir())[:60]
            except OSError:
                return cur, []
    return p, []


def _did_you_mean(name: str, candidates: "list[str]") -> "list[str]":
    """Candidates that contain, or are contained by, the requested name.

    Substring rather than edit distance because the miss is nearly always a
    missing timestamp prefix or a truncated suffix, not a typo.
    """
    n = name.lower()
    return [c for c in candidates
            if n and (n in c.lower() or c.lower() in n)][:8]


def list_artifacts(path: str = "Docs", *, limit: int = 200) -> dict:
    # Resolve WITHOUT requiring existence first, so a path that is merely
    # absent can be answered helpfully instead of raising. Containment is
    # still enforced below by _resolve.
    want = Path(path).expanduser()
    if not want.is_absolute():
        want = _root() / want
    if path and not want.exists():
        parent, siblings = _nearest_existing(want)
        guesses = _did_you_mean(want.name, siblings)
        try:
            rel_parent = str(parent.relative_to(_root()))
        except ValueError:
            rel_parent = str(parent)
        return {"path": path, "exists": False,
                 "note": (f"{path!r} does not exist. Nearest existing "
                           f"directory is {rel_parent!r}; its contents are "
                           f"listed under 'available'."
                           + (f" Closest matches: {guesses}." if guesses else "")),
                 "nearest_existing": rel_parent,
                 "did_you_mean": guesses,
                 "available": siblings,
                 "entries": []}
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
    return {"path": str(base.relative_to(_root())), "exists": True,
             "entries": entries}


def _open_maybe_gz(path):
    """Text handle for a plain or gzipped artefact.

    Local rather than imported from igvf_submission_skill: artifact reading is
    a core capability and should not acquire a dependency on a submission
    skill to open a file.
    """
    import gzip
    with open(path, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    return (gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz
            else open(path, "r", encoding="utf-8", errors="replace"))


def _sniff_delim(sample: str) -> str:
    """Tab or comma, whichever is more frequent in the header line.

    csv.Sniffer guesses from the whole sample and gets confused by quoted
    free-text fields, which these manifests have (biosample names contain
    commas). The header is the reliable signal.
    """
    head = (sample.splitlines() or [""])[0]
    return "\t" if head.count("\t") >= head.count(",") else ","


def rank_artifact(path: str, *, column: str, n: int = 10,
                   ascending: bool = False,
                   where: "str | None" = None,
                   where_column: "str | None" = None,
                   exclude: "str | None" = None) -> dict:
    """Top-N rows of a delimited artefact by one column, over the WHOLE file.

    This exists because its absence produced false claims. Asked for "the
    highest-scoring kidney records", the agent had only grep_artifacts
    (bounded hits) and read_artifact (truncated view), so it returned
    whatever those surfaced and described it as a ranking. Independently
    audited against the Catalog, the rows were genuine and correctly
    attributed -- but for GATA3 it reported 0.9909405 as the top score when
    0.9999999981 was present in the same file, and the same happened for SOX9
    and WT1.

    So the ranking is done here, deterministically, rather than inferred from
    a sample:

      * the entire file is read, not a head or a grep window;
      * rows whose ranking column does not parse as a number are counted and
        excluded, never silently treated as zero;
      * `n_scanned` and `n_ranked` are returned so a caller can see the
        ranking covered everything;
      * ties are broken by first appearance, so repeated runs agree.

    `where` filters BEFORE ranking, which is what "highest-scoring kidney
    record" actually asks for -- filter, then rank. It takes a COMMA-SEPARATED
    list and keeps a row matching ANY term, with `exclude` removing rows
    matching any of its own terms, because one substring cannot express the
    question these manifests are actually asked:

        kidney                 misses "renal cortical epithelial cell"
        renal                  also matches "ADRENAL gland"

    The independent audit of this data required "an explicit kidney or renal
    cortex context, excluding adrenal records", which is
    `where="kidney,renal" exclude="adrenal"` -- and the auditors noted that
    their own first pass used bare `renal` and had to be corrected for exactly
    the adrenal collision. Exclusion is applied after inclusion, so a term
    appearing in both loses.
    """
    p = _resolve(path)
    delim = None
    rows: "list[dict]" = []
    with _open_maybe_gz(p) as fh:
        sample = fh.read(64_000)
        fh.seek(0)
        delim = _sniff_delim(sample)
        reader = csv.DictReader(fh, delimiter=delim)
        fieldnames = reader.fieldnames or []
        if column not in fieldnames:
            return {"error": (f"no column {column!r} in {path}; columns are "
                               f"{fieldnames[:25]}"),
                     "columns": fieldnames}
        if where_column and where_column not in fieldnames:
            return {"error": (f"no column {where_column!r} in {path}; columns "
                               f"are {fieldnames[:25]}")}
        n_scanned = n_unparseable = n_filtered_out = 0
        for row in reader:
            n_scanned += 1
            if where or exclude:
                hay = ((row.get(where_column) or "") if where_column
                        else " ".join(str(v) for v in row.values())).lower()
                keep = True
                if where:
                    keep = any(t.strip().lower() in hay
                                for t in where.split(",") if t.strip())
                if keep and exclude:
                    keep = not any(t.strip().lower() in hay
                                    for t in exclude.split(",") if t.strip())
                if not keep:
                    n_filtered_out += 1
                    continue
            raw = row.get(column)
            try:
                val = float(raw)
            except (TypeError, ValueError):
                n_unparseable += 1
                continue
            rows.append({"_value": val, **row})
    # Stable: sort only on the value, so equal values keep file order.
    rows.sort(key=lambda r: r["_value"], reverse=not ascending)
    top = rows[:max(1, int(n))]
    return {"path": str(p.relative_to(_root())), "column": column,
             "ascending": ascending, "where": where,
             "where_column": where_column, "exclude": exclude,
             "n_scanned": n_scanned,
             "n_ranked": len(rows),
             "n_excluded_unparseable": n_unparseable,
             "n_excluded_by_filter": n_filtered_out,
             "ranking_is_complete": True,
             "rows": [{k: v for k, v in r.items() if k != "_value"}
                       for r in top],
             "values": [r["_value"] for r in top]}


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

    t = sub.add_parser("top", help="Top-N rows by a column, ranked over the "
                                    "WHOLE file (use this for any 'highest' "
                                    "or 'top' claim, never grep)")
    t.add_argument("--path", required=True)
    t.add_argument("--column", required=True)
    t.add_argument("--n", type=int, default=10)
    t.add_argument("--ascending", action="store_true")
    t.add_argument("--where", help="Keep only rows containing this substring "
                                    "before ranking.")
    t.add_argument("--where-column", help="Restrict --where to one column.")
    t.add_argument("--exclude", help="Comma-separated terms; drop rows "
                                      "matching any (e.g. adrenal).")

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
        elif args.cmd == "top":
            print(json.dumps(rank_artifact(
                args.path, column=args.column, n=args.n,
                ascending=args.ascending, where=args.where,
                where_column=args.where_column,
                exclude=args.exclude), indent=2))
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
