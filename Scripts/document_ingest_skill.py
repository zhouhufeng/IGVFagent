"""Read an uploaded manuscript and turn it into a reproduction plan.

A user with a paper PDF wants to ask "reproduce this analysis". The pieces
already existed — `bench harvest` scans text for repository accessions and
assay families, `bench route` maps those onto an IGVFagent chain, `bench
scaffold` writes a runnable benchmark — but nothing could turn a PDF into the
text those steps consume. Papers reached the system only through a DOI or PMID
that resolved to Europe PMC full text, so an unpublished, embargoed, or
paywalled manuscript could not be used at all.

This skill closes that gap:

    PDF / DOCX / text  ->  extracted text  ->  accessions + assays
                       ->  suggested `igvfagent bench` chain

Extraction is deliberately conservative about what it claims. A scanned PDF
with no text layer yields no characters, and that is reported as "no
extractable text — this looks like a scanned image" rather than as a paper
with no data. Silently returning zero accessions for an unreadable file is the
failure mode most likely to be believed.

Images are accepted and stored, but **not** read: there is no OCR here, and
pretending otherwise would produce confident empty results. They are kept
alongside the run so a user can point the agent at a figure and discuss it.

Usage::

    igvfagent document read    --path paper.pdf
    igvfagent document analyze --path paper.pdf        # + accessions/assays
    igvfagent document plan    --path paper.pdf        # + bench chain

Requires `pypdf` for PDFs (declared in the `analysis` extra); DOCX and plain
text need nothing beyond the standard library.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

__all__ = ["main", "extract_text", "analyze_text"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
OUT_DIR = ROOT / "Docs" / "Documents"
UPLOAD_DIR = ROOT / "Data" / "Uploads"

_IMAGE_EXT = {".png", ".jpg", ".jpeg", ".gif", ".tif", ".tiff", ".bmp", ".webp"}
_TEXT_EXT = {".txt", ".md", ".text", ".csv", ".tsv", ".xml", ".html", ".htm"}


class UnreadableDocument(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def _pdf_text(path: Path, max_pages: int = 300) -> "tuple[str, dict]":
    try:
        import pypdf
    except ImportError as e:
        raise UnreadableDocument(
            "reading PDFs needs `pypdf`. Install with: "
            "pip install 'igvfagent[analysis]'  (or pip install pypdf)") from e
    try:
        reader = pypdf.PdfReader(str(path))
    except Exception as e:
        raise UnreadableDocument(f"could not open PDF: {type(e).__name__}: {e}") from e

    if getattr(reader, "is_encrypted", False):
        try:
            reader.decrypt("")          # many papers are "encrypted" with no password
        except Exception:
            raise UnreadableDocument("PDF is password-protected")

    pages = reader.pages[:max_pages]
    chunks = []
    for p in pages:
        try:
            chunks.append(p.extract_text() or "")
        except Exception:
            chunks.append("")
    text = "\n".join(chunks)
    meta = {"pages": len(reader.pages), "pages_read": len(pages),
            "chars": len(text)}
    if len(text.strip()) < 200:
        # A PDF of page images has no text layer. Saying "0 accessions found"
        # here would be a true statement about a document we never read.
        raise UnreadableDocument(
            f"no extractable text ({len(text.strip())} chars from "
            f"{len(pages)} page(s)) — this looks like a scanned image PDF. "
            "OCR is not available; supply a text-layer PDF or paste the "
            "Data Availability section directly.")
    return text, meta


def _docx_text(path: Path) -> "tuple[str, dict]":
    """DOCX is a zip of XML; strip tags rather than add a dependency."""
    try:
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8", errors="replace")
    except Exception as e:
        raise UnreadableDocument(f"could not open DOCX: {type(e).__name__}") from e
    xml = re.sub(r"</w:p>", "\n", xml)
    text = re.sub(r"<[^>]+>", "", xml)
    text = re.sub(r"[ \t]+", " ", text)
    if len(text.strip()) < 100:
        raise UnreadableDocument("DOCX contained no readable text")
    return text, {"chars": len(text)}


def extract_text(path) -> "tuple[str, dict]":
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = ROOT / p
    if not p.is_file():
        raise UnreadableDocument(f"no such file: {p}")
    ext = p.suffix.lower()
    if ext == ".pdf":
        text, meta = _pdf_text(p)
    elif ext == ".docx":
        text, meta = _docx_text(p)
    elif ext in _TEXT_EXT:
        text = p.read_text(encoding="utf-8", errors="replace")
        meta = {"chars": len(text)}
    elif ext in _IMAGE_EXT:
        raise UnreadableDocument(
            f"{ext} is an image. There is no OCR in this skill, so no text can "
            "be extracted. The file is stored and can be viewed, but it cannot "
            "be searched for accessions.")
    else:
        raise UnreadableDocument(f"unsupported file type: {ext or '(none)'}")
    meta.update({"path": str(p), "type": ext.lstrip(".")})
    return text, meta


# ---------------------------------------------------------------------------
# Analysis — reuses the benchmark skill's own patterns
# ---------------------------------------------------------------------------

def _bench():
    try:
        from igvfagent import benchmark_skill as bs
    except Exception:
        import benchmark_skill as bs  # type: ignore
    return bs


def analyze_text(text: str) -> dict:
    """Accessions, assays and the data-availability statement.

    Deliberately reuses ACCESSION_PATTERNS and ASSAY_TERMS from the benchmark
    skill rather than restating them, so a document analysed here and a paper
    harvested from Europe PMC are read by exactly the same rules.
    """
    bs = _bench()
    low = text.lower()

    accessions = {}
    for pat in getattr(bs, "ACCESSION_PATTERNS", []):
        hits = re.findall(pat["re"], text, re.I)
        if not hits:
            continue
        flat = []
        for h in hits:
            flat.append(h if isinstance(h, str) else next((x for x in h if x), ""))
        uniq = sorted({h for h in flat if h})
        accessions[pat["key"]] = {"repo": pat["repo"], "ids": uniq[:50],
                                   "n": len(uniq),
                                   "controlled": bool(pat.get("controlled"))}

    assays = []
    for a in getattr(bs, "ASSAY_TERMS", []):
        hit = [t for t in a["terms"] if t in low]
        if hit:
            assays.append({"assay": a["assay"], "spec": a.get("spec", 1.0),
                            "matched": hit,
                            "mentions": sum(low.count(t) for t in hit)})
    # Specificity first: naming a particular protocol is stronger evidence for
    # a route than a generic modality every paper mentions somewhere.
    assays.sort(key=lambda a: (-a["spec"], -a["mentions"]))

    return {
        "accessions": accessions,
        "n_accessions": sum(v["n"] for v in accessions.values()),
        "assays": assays,
        "data_availability": _data_availability(text),
        "doi": _first(re.findall(r"\b10\.\d{4,9}/[^\s,;)\]]+", text)),
        "pmid": _first(re.findall(r"PMID:?\s*(\d{6,8})", text, re.I)),
    }


def _first(seq):
    return seq[0].rstrip(".") if seq else None


def _data_availability(text: str) -> "list[str]":
    """Pull the data-availability paragraphs — the densest source of accessions."""
    out = []
    pat = re.compile(
        r"(data\s+(?:and\s+code\s+)?availability|availability\s+of\s+data"
        r"|accession\s+codes?|data\s+access)", re.I)
    for m in pat.finditer(text):
        chunk = " ".join(text[m.start(): m.start() + 1200].split())
        out.append(chunk)
        if len(out) >= 3:
            break
    return out


def suggest_chain(analysis: dict) -> dict:
    """Map accessions + assays onto an IGVFagent bench chain."""
    bs = _bench()
    assays = [a["assay"] for a in analysis["assays"][:6]]
    try:
        routed = bs.route({"assays": [{"assay": a} for a in analysis["assays"]],
                           "accessions": analysis["accessions"]}, top=3)
    except Exception:
        routed = None

    ids = []
    for key, v in analysis["accessions"].items():
        for i in v["ids"][:3]:
            ids.append((v["repo"], i, v["controlled"]))

    steps = []
    if analysis.get("doi"):
        steps.append(f"igvfagent bench resolve --query {analysis['doi']}")
    for repo, acc, controlled in ids[:6]:
        if controlled:
            steps.append(f"# {acc} ({repo}) is controlled-access — needs an "
                         f"approved application, cannot be fetched")
        elif repo == "NCBI GEO":
            steps.append(f"igvfagent geo retrieve --accession {acc}")
        elif repo == "IGVF Portal":
            steps.append(f"igvfagent explain --target {acc}")
        elif repo == "ENCODE":
            steps.append(f"igvfagent explain --target {acc}")
        elif repo == "MaveDB":
            steps.append(f"igvfagent mavedb map --urn {acc}")
    if analysis.get("doi"):
        steps.append(f"igvfagent bench pipeline --query {analysis['doi']} "
                     f"--execute   # full scaffold -> run -> score")
    return {"assays": assays, "routed": routed, "suggested_commands": steps}


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _run_dir(label: str) -> Path:
    d = OUT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{re.sub(r'[^A-Za-z0-9_.-]', '_', label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cmd_read(args) -> int:
    try:
        text, meta = extract_text(args.path)
    except UnreadableDocument as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    d = _run_dir(args.label or Path(args.path).stem)
    (d / "text.txt").write_text(text)
    (d / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"Extracted:     {meta['chars']:,} characters"
          + (f" from {meta['pages']} page(s)" if meta.get("pages") else ""))
    print(f"Wrote:         {d / 'text.txt'}")
    if args.head:
        print("\n" + " ".join(text.split())[: args.head])
    return 0


def cmd_analyze(args) -> int:
    try:
        text, meta = extract_text(args.path)
    except UnreadableDocument as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    a = analyze_text(text)
    d = _run_dir(args.label or Path(args.path).stem)
    (d / "text.txt").write_text(text)
    (d / "analysis.json").write_text(json.dumps({"meta": meta, **a}, indent=2))
    (d / "report.md").write_text(_render(meta, a, None))
    _print_analysis(meta, a)
    print(f"Report:        {d / 'report.md'}")
    print(f"Wrote:         {d / 'analysis.json'}")
    return 0


def cmd_plan(args) -> int:
    try:
        text, meta = extract_text(args.path)
    except UnreadableDocument as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    a = analyze_text(text)
    chain = suggest_chain(a)
    d = _run_dir(args.label or Path(args.path).stem)
    (d / "text.txt").write_text(text)
    (d / "analysis.json").write_text(json.dumps({"meta": meta, **a,
                                                  "chain": chain}, indent=2))
    (d / "report.md").write_text(_render(meta, a, chain))
    _print_analysis(meta, a)
    print("\nSuggested reproduction chain:")
    for s in chain["suggested_commands"]:
        print(f"  {s}")
    if not chain["suggested_commands"]:
        print("  (no runnable step — no public accession was found in the text)")
    print(f"\nReport:        {d / 'report.md'}")
    return 0


def _print_analysis(meta, a) -> None:
    print(f"Extracted:     {meta['chars']:,} characters"
          + (f" from {meta['pages']} page(s)" if meta.get("pages") else ""))
    print(f"DOI:           {a.get('doi') or '(none found)'}")
    print(f"Accessions:    {a['n_accessions']}")
    for key, v in a["accessions"].items():
        flag = "  [CONTROLLED ACCESS]" if v["controlled"] else ""
        print(f"  {v['repo']:22s} {', '.join(v['ids'][:5])}{flag}")
    print(f"Assays:        {', '.join(x['assay'] for x in a['assays'][:8]) or '(none detected)'}")
    if a["data_availability"]:
        print(f"Data availability statement found "
              f"({len(a['data_availability'])} passage(s))")


def _render(meta, a, chain) -> str:
    L = [f"# Document analysis — {Path(meta['path']).name}", "",
         f"- Type: `{meta['type']}`",
         f"- Extracted: **{meta['chars']:,}** characters"
         + (f" from {meta['pages']} pages" if meta.get("pages") else ""),
         f"- DOI: `{a.get('doi') or 'not found'}`", "",
         "## Accessions", ""]
    if a["accessions"]:
        L += ["| repository | ids | controlled |", "|---|---|---|"]
        for key, v in a["accessions"].items():
            L.append(f"| {v['repo']} | `{'`, `'.join(v['ids'][:8])}` | "
                     f"{'**yes**' if v['controlled'] else 'no'} |")
    else:
        L.append("None found. If this paper deposits data, the statement may "
                 "sit in a supplement rather than the main PDF.")
    L += ["", "## Assays detected", ""]
    L += ([f"- **{x['assay']}** (specificity {x['spec']}, {x['mentions']} mentions)"
           for x in a["assays"][:10]] or ["None detected."])
    if a["data_availability"]:
        L += ["", "## Data availability (as written)", ""]
        L += [f"> {p[:900]}" for p in a["data_availability"]]
    if chain:
        L += ["", "## Suggested reproduction chain", "", "```bash"]
        L += chain["suggested_commands"] or ["# no public accession found"]
        L += ["```"]
    L += ["", "## Caveats", "",
          "- Accessions and assays are detected by the same patterns "
          "`igvfagent bench harvest` uses, so a document read here and a paper "
          "resolved from Europe PMC are interpreted identically.",
          "- Detection is textual. An accession that appears only in a figure, "
          "a supplementary file, or a scanned page will not be found.",
          "- Controlled-access deposits are flagged but cannot be retrieved "
          "without an approved application.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent document",
        description="Read an uploaded manuscript (PDF/DOCX/text) and derive a "
                    "reproduction plan from it.")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, helptext in (("read", "Extract text only"),
                            ("analyze", "Extract + find accessions and assays"),
                            ("plan", "Analyze + suggest a bench chain")):
        s = sub.add_parser(name, help=helptext)
        s.add_argument("--path", required=True)
        s.add_argument("--label", default="")
        if name == "read":
            s.add_argument("--head", type=int, default=0,
                           help="Print the first N characters.")
    args = p.parse_args(argv)
    return {"read": cmd_read, "analyze": cmd_analyze, "plan": cmd_plan}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
