#!/usr/bin/env python3
"""Gene-regulatory-network (dEx) and protein-variant effect queries against
the IGVF Catalog API.

Two Catalog collections that no other IGVFagent skill reaches:

``/api/gene-regulatory-network``
    Differential-expression gene regulatory network. One row per
    regulator-element -> response-gene edge, carrying ``log2FC``,
    ``neg_log10_pvalue``, ``significant``, ``crispr_modality`` and the
    perturbed ``genomic_element`` coordinates. Queryable from either end
    (``--regulator`` / ``--response``), which makes three distinct
    questions answerable:

      * regulator -> targets    (``--regulator EOMES``)
      * target -> regulators    (``--response ARHGEF3``)
      * self-edge check         (both, same symbol) — a TF perturbation
        that changes the TF's own expression. Not all TFs have one.

``/api/proteins/variants``
    Sequence-variant -> protein effects: allele-specific binding and
    motif-disruption calls from SEMVAR, ADASTRA and friends.

Ported from the IGVF QA notebooks in
``Docs/References/IGVF_Catalog-dev/``. Differences from the notebooks,
all deliberate:

  * Defaults to the production Catalog (``api.catalogkg.igvf.org``) like
    every other skill here. ``--host dev`` selects the dev/demo API the
    notebooks use. Both collections are served by both hosts.
  * REST only. The notebooks also open a direct ArangoDB connection with
    a hardcoded ``guest`` password; that needs ``python-arango`` and adds
    a credential to the repo for no capability gain, so it is dropped.
  * Pagination is the notebooks' ``paginated()`` helper, with the
    page-walk bounded by ``--max-results`` so a broad filter cannot run
    away.

Usage
-----
    igvfagent grn network --regulator EOMES --p-value 0.05
    igvfagent grn network --response ARHGEF3
    igvfagent grn network --regulator KLF2 --response KLF2      # self-edge
    igvfagent grn network --response CCR7 --method 'CRISPR screen'
    igvfagent grn protein-variants --protein ELF2
    igvfagent grn protein-variants --protein TP53 --method ADASTRA
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "CatalogGRN"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

USER_AGENT = "IGVFagent-catalog-grn/0.1"

GRN_PATH = "/api/gene-regulatory-network"
PROTEIN_VARIANTS_PATH = "/api/proteins/variants"

# Observed ``method`` values on these two collections. Not a closed set —
# the Catalog adds methods as datasets land — so these are documentation
# for the caller, never a client-side filter.
GRN_METHODS = ("Perturb-seq", "CRISPR screen")
PROTEIN_VARIANT_METHODS = ("SEMVAR", "ADASTRA")


def resolve_base(host: str) -> str:
    """``prod`` -> the production Catalog API, ``dev`` -> the dev/demo API.

    ``IGVF_CATALOG_API_BASE`` overrides the prod host (consistent with
    catalog_query_skill); ``IGVF_CATALOG_API_DEV_BASE`` overrides dev.
    """
    if host == "dev":
        return _resolve_endpoint("catalog_api_dev", "IGVF_CATALOG_API_DEV_BASE")
    return _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"catalog_grn_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


# ─── HTTP ───────────────────────────────────────────────────────────────────

def _request(url: str, *, timeout: float = 120.0,
             retries: int = 3) -> "tuple[int, bytes, str]":
    """GET with retry on transport errors and 5xx.

    The Catalog routinely takes 15-25s per page on these two collections,
    and an occasional request times out outright. Without a retry a single
    blip aborts a page-walk that is otherwise minutes deep, so transport
    errors and 5xx get a bounded exponential backoff. 4xx is returned as-is
    — it is a filter problem the caller needs to see, not a flake.
    """
    req = urllib.request.Request(url, headers={
        "User-Agent": USER_AGENT,
        "Accept": "application/json",
    })
    last: "tuple[int, bytes, str]" = (0, b"no attempt made", "text/plain")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return (resp.getcode(), resp.read(),
                        resp.headers.get("Content-Type", ""))
        except urllib.error.HTTPError as e:
            body, ct = e.read(), e.headers.get("Content-Type", "")
            if e.code < 500:
                return e.code, body, ct
            last = (e.code, body, ct)
            logging.warning("  HTTP %s (attempt %d/%d)",
                            e.code, attempt + 1, retries)
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            reason = getattr(e, "reason", e)
            last = (0, str(reason).encode(), "text/plain")
            logging.warning("  transport error %s (attempt %d/%d)",
                            reason, attempt + 1, retries)
        if attempt < retries - 1:
            time.sleep(2 ** attempt)
    logging.error("  giving up after %d attempts", retries)
    return last


def _get_json(base: str, path: str,
              params: "list[tuple[str, str]]") -> Any:
    qs = urllib.parse.urlencode(params, doseq=True)
    url = f"{base}{path}?{qs}" if qs else f"{base}{path}"
    logging.info("GET %s", url)
    status, content, ct = _request(url)
    if status < 200 or status >= 300:
        try:
            err = json.loads(content)
            detail = err.get("message") or err.get("detail") or str(err)
        except Exception:
            detail = content[:300].decode(errors="replace")
        # The Catalog rejects an unfiltered collection scan with a 400
        # naming the acceptable filters. Pass that through verbatim — it
        # is the most useful thing we could say.
        raise SystemExit(f"{path}: HTTP {status}: {detail}")
    if "json" not in ct:
        raise SystemExit(f"{path}: expected JSON, got {ct!r}")
    return json.loads(content)


def paginated(base: str, path: str, params: "list[tuple[str, str]]", *,
              page_size: int = 500, max_results: int = 10000) -> "list[dict]":
    """Walk 0-indexed pages until a short page, ``max_results``, or empty.

    Mirrors the notebooks' ``paginated()`` helper. A full page is only a
    *hint* that more exist, so the loop stops on the first short page
    rather than trusting any total count.
    """
    out: "list[dict]" = []
    page = 0
    while len(out) < max_results:
        batch = _get_json(base, path,
                          list(params) + [("page", str(page)),
                                          ("limit", str(page_size))])
        if not isinstance(batch, list):
            raise SystemExit(
                f"{path}: expected a JSON array, got {type(batch).__name__}")
        out.extend(batch)
        logging.info("  page %d: %d rows (running total %d)",
                     page, len(batch), len(out))
        if len(batch) < page_size:
            break
        page += 1
    return out[:max_results]


# ─── Shaping ────────────────────────────────────────────────────────────────

def flatten(rec: dict) -> dict:
    """One level of dict-flattening with a ``parent_child`` key prefix.

    The notebooks do this so nested ``genomic_element`` coordinates land
    as real columns instead of a stringified dict.
    """
    flat: dict = {}
    for k, v in rec.items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                flat[f"{k}_{kk}"] = vv
        else:
            flat[k] = v
    return flat


def write_outputs(rows: "list[dict]", label: str, *,
                  query: dict) -> "tuple[Path, Path, Optional[Path]]":
    """Write JSON + Markdown (+ TSV when non-empty). Returns the paths."""
    stamp = time.strftime("%Y%m%d_%H%M%S")
    outdir = REPORT_DIR / f"{stamp}_{label}"
    outdir.mkdir(parents=True, exist_ok=True)

    flat = [flatten(r) for r in rows]

    raw = outdir / "records.json"
    raw.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    tsv: Optional[Path] = None
    if flat:
        cols: "list[str]" = []
        for r in flat:
            for k in r:
                if k not in cols:
                    cols.append(k)
        tsv = outdir / "records.tsv"
        with tsv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(flat)

    report = outdir / "report.md"
    report.write_text(_render_report(flat, label, query=query, tsv=tsv),
                      encoding="utf-8")
    return report, raw, tsv


def _render_report(flat: "list[dict]", label: str, *,
                   query: dict, tsv: Optional[Path]) -> str:
    lines = [f"# Catalog query — {label}", ""]
    lines.append("## Query")
    lines.append("")
    for k, v in query.items():
        if v not in (None, ""):
            lines.append(f"- **{k}**: `{v}`")
    lines.append("")
    lines.append(f"**Rows returned:** {len(flat)}")
    lines.append("")

    if not flat:
        lines += [
            "No rows matched. This is a real answer, not an error — the "
            "filter combination has no edges in this Catalog build. Worth "
            "trying: drop `--p-value`, widen `--method`, or check the "
            "symbol is the Catalog's preferred name.",
            "",
        ]
        return "\n".join(lines)

    # Significance / effect-size summary, when the GRN columns are present.
    sig = [r for r in flat if r.get("significant") is True]
    if any("significant" in r for r in flat):
        lines.append(f"**Significant edges:** {len(sig)} / {len(flat)}")
        lines.append("")
    l2 = [r["log2FC"] for r in flat
          if isinstance(r.get("log2FC"), (int, float))]
    if l2:
        lines.append(f"**log2FC range:** {min(l2):.3f} … {max(l2):.3f}")
        lines.append("")

    for col, title in (("method", "Methods"),
                       ("source", "Sources"),
                       ("crispr_modality", "CRISPR modality"),
                       ("biological_context", "Biological contexts"),
                       ("label", "Effect labels")):
        counts: "dict[str, int]" = {}
        for r in flat:
            v = r.get(col)
            if v not in (None, ""):
                counts[str(v)] = counts.get(str(v), 0) + 1
        if counts:
            lines.append(f"## {title}")
            lines.append("")
            for v, n in sorted(counts.items(), key=lambda kv: -kv[1])[:12]:
                lines.append(f"- {v} — {n}")
            lines.append("")

    lines.append("## First rows")
    lines.append("")
    preview_cols = [c for c in (
        "response_gene", "genomic_element_regulator_gene",
        "genomic_element_chr", "genomic_element_start", "genomic_element_end",
        "log2FC", "neg_log10_pvalue", "significant", "crispr_modality",
        "sequence_variant", "protein_complex", "label", "name",
        "method", "source", "biological_context",
    ) if any(c in r for r in flat)]
    if preview_cols:
        lines.append("| " + " | ".join(preview_cols) + " |")
        lines.append("|" + "|".join(["---"] * len(preview_cols)) + "|")
        for r in flat[:15]:
            lines.append("| " + " | ".join(
                str(r.get(c, "")) for c in preview_cols) + " |")
        lines.append("")
    if tsv:
        lines.append(f"Full table: `{tsv}`")
        lines.append("")
    return "\n".join(lines)


# ─── Subcommands ────────────────────────────────────────────────────────────

def cmd_network(args: argparse.Namespace) -> int:
    if not (args.regulator or args.response or args.method
            or args.files_fileset):
        raise SystemExit(
            "Give at least one filter: --regulator, --response, --method, "
            "or --files-fileset. The Catalog rejects an unfiltered scan of "
            "this collection.")

    params: "list[tuple[str, str]]" = []
    if args.regulator:
        params.append(("regulator_gene_name", args.regulator))
    if args.response:
        params.append(("response_gene_name", args.response))
    if args.method:
        params.append(("method", args.method))
    if args.files_fileset:
        params.append(("files_fileset", args.files_fileset))
    if args.p_value is not None:
        # The Catalog's comparison-filter syntax, straight from the
        # notebooks: p_value=lte:0.05.
        params.append(("p_value", f"lte:{args.p_value}"))

    base = resolve_base(args.host)
    rows = paginated(base, GRN_PATH, params,
                     page_size=args.page_size, max_results=args.max_results)

    self_edge = bool(args.regulator and args.response
                     and args.regulator == args.response)
    label = args.label or _default_label(
        "grn", args.regulator, args.response, self_edge=self_edge)

    report, raw, tsv = write_outputs(rows, label, query={
        "endpoint": f"{base}{GRN_PATH}",
        "regulator_gene_name": args.regulator,
        "response_gene_name": args.response,
        "method": args.method,
        "files_fileset": args.files_fileset,
        "p_value": f"lte:{args.p_value}" if args.p_value is not None else None,
        "host": args.host,
        "query_type": ("self-edge" if self_edge
                       else "regulator -> targets" if args.regulator
                       else "target -> regulators" if args.response
                       else "method/fileset scan"),
    })

    sig = sum(1 for r in rows if r.get("significant") is True)
    print(f"Rows:          {len(rows)}  (significant: {sig})")
    if self_edge and not rows:
        print("Self-edge:     none found — expected for many TFs.")
    print(f"Report:        {report}")
    print(f"Records:       {raw}")
    if tsv:
        print(f"Manifest:      {tsv}")
    return 0


def cmd_protein_variants(args: argparse.Namespace) -> int:
    if not (args.protein or args.method or args.files_fileset):
        raise SystemExit(
            "Give at least one filter: --protein, --method, or "
            "--files-fileset. The Catalog rejects an unfiltered scan of "
            "this collection.")

    params: "list[tuple[str, str]]" = []
    if args.protein:
        params.append(("protein_name", args.protein))
    if args.method:
        params.append(("method", args.method))
    if args.source:
        params.append(("source", args.source))
    if args.files_fileset:
        params.append(("files_fileset", args.files_fileset))

    base = resolve_base(args.host)
    rows = paginated(base, PROTEIN_VARIANTS_PATH, params,
                     page_size=args.page_size, max_results=args.max_results)

    label = args.label or _default_label(
        "protvar", args.protein, args.method)
    report, raw, tsv = write_outputs(rows, label, query={
        "endpoint": f"{base}{PROTEIN_VARIANTS_PATH}",
        "protein_name": args.protein,
        "method": args.method,
        "source": args.source,
        "files_fileset": args.files_fileset,
        "host": args.host,
    })

    print(f"Rows:          {len(rows)}")
    print(f"Report:        {report}")
    print(f"Records:       {raw}")
    if tsv:
        print(f"Manifest:      {tsv}")
    return 0


def _default_label(prefix: str, *parts: Any, self_edge: bool = False) -> str:
    bits = [prefix] + [str(p) for p in parts if p]
    if self_edge:
        bits.append("selfedge")
    slug = "_".join(bits)
    return "".join(c if (c.isalnum() or c in "_-") else "-" for c in slug)[:60]


# ─── CLI ────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent grn",
        description="IGVF Catalog gene-regulatory-network (dEx) and "
                    "protein-variant effect queries.")
    sub = p.add_subparsers(dest="command", required=True)

    def _common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--host", choices=("prod", "dev"), default="prod",
                        help="Catalog API to query. 'prod' (default) is "
                             "api.catalogkg.igvf.org; 'dev' is the "
                             "dev/demo API used by the IGVF QA notebooks "
                             "(pre-release collections, no uptime "
                             "guarantee).")
        sp.add_argument("--page-size", type=int, default=500,
                        help="Rows per API request (default 500).")
        sp.add_argument("--max-results", type=int, default=10000,
                        help="Stop after this many rows (default 10000).")
        sp.add_argument("--files-fileset",
                        help="Restrict to one Catalog fileset, e.g. "
                             "files_filesets/IGVFFI4081AUXT.")
        sp.add_argument("--label", help="Slug for the output directory.")

    n = sub.add_parser(
        "network",
        help="Differential-expression gene regulatory network edges.")
    n.add_argument("--regulator",
                   help="Regulator (TF) gene symbol — returns its targets.")
    n.add_argument("--response",
                   help="Response gene symbol — returns its regulators.")
    n.add_argument("--method",
                   help="Assay method, e.g. " +
                        " / ".join(f"'{m}'" for m in GRN_METHODS) + ".")
    n.add_argument("--p-value", type=float, default=None,
                   help="Keep edges with p <= this value (sent as "
                        "p_value=lte:<v>). Typical: 0.05.")
    _common(n)
    n.set_defaults(func=cmd_network)

    v = sub.add_parser(
        "protein-variants",
        help="Sequence-variant effects on proteins (allele-specific "
             "binding, motif disruption).")
    v.add_argument("--protein", help="Protein / gene symbol, e.g. ELF2.")
    v.add_argument("--method",
                   help="Scoring method, e.g. " +
                        " / ".join(f"'{m}'" for m in PROTEIN_VARIANT_METHODS)
                        + ".")
    v.add_argument("--source", help="Originating source, e.g. IGVF, ADASTRA.")
    _common(v)
    v.set_defaults(func=cmd_protein_variants)

    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
