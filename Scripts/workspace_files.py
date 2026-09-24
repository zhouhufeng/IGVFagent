"""Read and write plain-text files inside the workspace (write-text-file,
read-text-file).

Hosted agents (Claude Code jobs over MCP included) have no shell and no native
file writer, yet a job must save section verdicts, notes and reports as stage
evidence. Before this built-in existed, agents authored their own writer
extension; one read its content from stdin when --content was missing, and
run from the MCP server that child swallowed the server's JSON-RPC stream and
hung every later call.

Rules: the path must resolve inside the workspace and pass the artefact path
guard (never Docs/Secret or credential files); source code, .git and the
extension directories are not writable; content comes only from --content
(never stdin); at most IGVF_WRITE_MAX_BYTES (default 5 MB) per write.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional, Sequence

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
MAX_BYTES = int(os.environ.get("IGVF_WRITE_MAX_BYTES", str(5 * 1024 * 1024)))
READONLY_TOP = {"Scripts", ".git", "Deploy", ".github", "src"}
TEXT_SUFFIXES = {".md", ".txt", ".json", ".tsv", ".csv", ".yaml", ".yml", ".html", ".log", ".bed", ".R", ".r",
                 ".py", ".Rmd", ".sh", ".pml", ".xml", ".tex", ".ini", ".cfg", ".toml", ".fa", ".fasta", ".gtf", ".vcf"}


def resolve(path: str, *, for_write: bool) -> "tuple[Optional[Path], Optional[str]]":
    q = Path(path)
    q = (q if q.is_absolute() else ROOT / q).resolve()
    try:
        rel = q.relative_to(ROOT)
    except ValueError:
        return None, f"{path}: outside the workspace"
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        try:
            from igvfagent import _pathguard  # type: ignore
        except Exception:
            import _pathguard  # type: ignore
        why = _pathguard.why_blocked(q, require_file=False)
    except Exception:  # noqa: BLE001
        why = "path guard unavailable" if "secret" in str(rel).lower() else None
    if why and for_write and "not a file" in why:
        why = None  # a new file: only the secret/credential rules apply
    if why:
        return None, f"{path}: {why}"
    if for_write:
        if rel.parts and (rel.parts[0] in READONLY_TOP or "UserExtensions" in rel.parts):
            return None, f"{path}: {'UserExtensions' if 'UserExtensions' in rel.parts else rel.parts[0]}/ is not writable from a tool"
        if q.suffix not in TEXT_SUFFIXES:
            return None, f"{path}: only text files ({', '.join(sorted(TEXT_SUFFIXES))})"
    return q, None


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent files", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write")
    w.add_argument("--path", required=True)
    w.add_argument("--content", required=True, help="full text (never read from stdin)")
    w.add_argument("--append", action="store_true")
    r = sub.add_parser("read")
    r.add_argument("--path", required=True)
    r.add_argument("--offset", type=int, default=0, help="first line (0-based)")
    r.add_argument("--limit", type=int, default=400, help="lines to return")
    a = ap.parse_args(argv)
    if a.cmd == "write":
        q, err = resolve(a.path, for_write=True)
        if err:
            print(json.dumps({"ok": False, "error": err}))
            return 2
        data = a.content if a.content.endswith("\n") or not a.content else a.content + "\n"
        if len(data.encode()) > MAX_BYTES:
            print(json.dumps({"ok": False, "error": f"content over {MAX_BYTES} bytes"}))
            return 2
        q.parent.mkdir(parents=True, exist_ok=True)
        with open(q, "a" if a.append else "w", encoding="utf-8") as fh:
            fh.write(data)
        rel = str(q.relative_to(ROOT))
        print(json.dumps({"ok": True, "path": rel, "bytes": q.stat().st_size, "append": a.append}))
        print(f"Wrote: {rel}")
        return 0
    q, err = resolve(a.path, for_write=False)
    if err or not q.is_file():
        print(json.dumps({"ok": False, "error": err or f"{a.path}: no such file"}))
        return 2
    lines = q.read_text(errors="replace").splitlines()
    chunk = lines[a.offset:a.offset + a.limit]
    print(f"# {q.relative_to(ROOT)} lines {a.offset + 1}-{a.offset + len(chunk)} of {len(lines)}")
    print("\n".join(chunk))
    return 0


if __name__ == "__main__":
    sys.exit(main())
