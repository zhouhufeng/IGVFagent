#!/usr/bin/env python3
"""Super-enhancer → target-gene pipeline.

End-to-end ENCODE-driven workflow that takes a cell line / tissue
(GM12878, K562, HepG2, liver, hippocampus, …), discovers the relevant
H3K27ac ChIP-seq + 3D-chromatin (Hi-C / ChIA-PET / capture Hi-C)
experiments, downloads the peak files, calls super-enhancers ROSE-style,
and links each super-enhancer to candidate target genes by combining

  * 3D loops          — anchor-overlap with Hi-C / ChIA-PET bedpe
  * SCREEN cCREs      — class assignment (PLS / ELS) within each SE
  * IGVF Catalog rE2G — region → gene predictions per biosample
  * Proximity         — nearest-TSS fallback when nothing else fires

Each (SE, gene) pair is scored by how many of those four methods support
it; the result is a ranked target-gene table, a per-SE drill-down, an
SE↔gene network plot, and a browser SVG for any gene of interest.

Subcommands

  pipeline      One-command run: discover → download → call SEs → link
                targets → write reports.
  discover      Just discovery: list candidate ChIP-seq + 3D-chromatin
                experiments for a biosample.
  link-targets  Given an existing SE BED + optional loops bedpe, run
                the linkage step in isolation.
  write-playbook  Emit Docs/Skills/SE_TARGET_PIPELINE_SKILLS.md.

Outputs land under ``Docs/SETargets/<timestamp>_<label>/`` with manifest
CSVs, the SE BED, the linked-genes CSV, plots, and a markdown report.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

# Reuse battle-tested helpers from encode_pipeline; lazy import so the
# package layout is the only thing that has to be right.
try:
    from igvfagent import encode_pipeline as ep
except Exception:
    import encode_pipeline as ep  # type: ignore


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "SETargets"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "SETargets"
DOWNLOAD_DIR = DATA_DIR / "IGVF" / "SETargets" / "Downloads"

ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
WENGLAB_DL = _resolve_endpoint("wenglab_dl")

USER_AGENT = "IGVFdataAgent-SE2Target/0.1"

# Default enhancer-mark target for SE calling. H3K27ac is canonical;
# BRD4 / MED1 / MED12 / P300 also work and follow the same code path.
DEFAULT_SE_TARGET = "H3K27ac"

# 3D-chromatin assay families that produce per-experiment .bedpe loops
THREE_D_ASSAYS = (
    "Hi-C", "in situ Hi-C", "intact Hi-C",
    "capture Hi-C",
    "ChIA-PET",
)

# Known method names we'll look for in IGVF Catalog regulatory-regions
# /genes responses to spot enhancer-gene predictions per biosample.
LINKAGE_METHODS = ("ENCODE-rE2G", "rE2G", "ABC", "scE2G")


# --------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"se_target_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (DATA_DIR, REPORT_DIR, PLOT_DIR, MANIFEST_DIR,
              DOWNLOAD_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# --------------------------- HTTP helpers ----------------------------------

def fetch_json(url: str, timeout: int = 60) -> "tuple[int, Any]":
    logging.info("GET %s", url)
    req = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT,
                        "Accept": "application/json"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace") if e.fp else ""
        try:
            return e.code, json.loads(body)
        except json.JSONDecodeError:
            return e.code, {"http_error_body": body[:600]}
    except urllib.error.URLError as e:
        return 0, {"network_error": str(e.reason)}
    except Exception as e:  # socket.timeout and friends
        return 0, {"network_error": str(e)}


def encode_get(path: str, **params) -> "tuple[int, Any]":
    url = ENCODE_BASE + path
    if params:
        url = (url + ("&" if "?" in url else "?")
               + urllib.parse.urlencode({k: v for k, v in params.items()
                                            if v is not None}, doseq=True))
    return fetch_json(url)


def catalog_get(path: str, timeout: int = 30, **params
                  ) -> "tuple[int, Any]":
    url = CATALOG_API_BASE + path
    if params:
        url = (url + ("&" if "?" in url else "?")
               + urllib.parse.urlencode({k: v for k, v in params.items()
                                            if v is not None}))
    return fetch_json(url, timeout=timeout)


def write_csv(path: Path, rows: "list[dict]",
                cols: "Optional[list[str]]" = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --------------------------- Discovery -------------------------------------

def discover_chipseq(biosample: str, target: str = DEFAULT_SE_TARGET,
                       assembly: str = "GRCh38",
                       limit: int = 25) -> "list[dict]":
    """Find Histone ChIP-seq experiments for a biosample × target ×
    assembly. Returns simplified records ranked by released-status +
    file-count."""
    rows: "list[dict]" = []
    for facet in ("Histone ChIP-seq", "TF ChIP-seq", "ChIP-seq"):
        sc, data = encode_get(
            "/search/", **{
                "type":          "Experiment",
                "format":        "json",
                "limit":         limit,
                "status":        "released",
                "assay_title":   facet,
                "biosample_ontology.term_name": biosample,
                "target.label":  target,
                "assembly":      assembly,
            },
        )
        if sc == 200 and isinstance(data, dict):
            for g in (data.get("@graph") or []):
                if isinstance(g, dict):
                    g["_assay_title"] = facet
                    rows.append(g)
        time.sleep(0.1)
    # Deduplicate by accession
    by_acc: "dict[str, dict]" = {}
    for r in rows:
        acc = r.get("accession", "")
        if acc and acc not in by_acc:
            by_acc[acc] = r
    return list(by_acc.values())


def discover_chromatin_accessibility(biosample: str,
                                          assembly: str = "GRCh38",
                                          limit: int = 25) -> "list[dict]":
    """DNase-seq / ATAC-seq experiments for the same biosample."""
    out: "list[dict]" = []
    for facet in ("DNase-seq", "ATAC-seq"):
        sc, data = encode_get(
            "/search/", **{
                "type":          "Experiment",
                "format":        "json",
                "limit":         limit,
                "status":        "released",
                "assay_title":   facet,
                "biosample_ontology.term_name": biosample,
                "assembly":      assembly,
            },
        )
        if sc == 200 and isinstance(data, dict):
            for g in (data.get("@graph") or []):
                if isinstance(g, dict):
                    g["_assay_title"] = facet
                    out.append(g)
        time.sleep(0.1)
    return out


def discover_3d(biosample: str, assembly: str = "GRCh38",
                 limit: int = 25) -> "list[dict]":
    """Hi-C / capture Hi-C / ChIA-PET experiments for the biosample."""
    rows: "list[dict]" = []
    for facet in THREE_D_ASSAYS:
        sc, data = encode_get(
            "/search/", **{
                "type":          "Experiment",
                "format":        "json",
                "limit":         limit,
                "status":        "released",
                "assay_title":   facet,
                "biosample_ontology.term_name": biosample,
                "assembly":      assembly,
            },
        )
        if sc == 200 and isinstance(data, dict):
            for g in (data.get("@graph") or []):
                if isinstance(g, dict):
                    g["_assay_title"] = facet
                    rows.append(g)
        time.sleep(0.1)
    return rows


def hydrate(accession: str) -> dict:
    sc, data = encode_get(f"/experiments/{accession}/", format="json")
    return data if sc == 200 and isinstance(data, dict) else {
        "accession": accession, "_http_status": sc,
    }


def pick_peak_file(meta: dict,
                    preferred_outputs=("replicated peaks", "stable peaks",
                                          "peaks"),
                    preferred_format="bed narrowPeak") -> "Optional[dict]":
    """Choose the best peaks file from an experiment for SE calling.
    Prefers replicated > stable > peaks; narrowPeak > broadPeak BED."""
    files = meta.get("files") or []
    candidates = []
    for f in files:
        if not isinstance(f, dict) or f.get("status") != "released":
            continue
        ftype = (f.get("file_type") or "").lower()
        otype = (f.get("output_type") or "").lower()
        if "narrowpeak" not in ftype and "broadpeak" not in ftype:
            continue
        score = 0
        for i, want in enumerate(preferred_outputs):
            if want in otype:
                score += (10 - i)
        if "narrowpeak" in ftype:
            score += 3
        f_simplified = {**f, "_score": score}
        candidates.append(f_simplified)
    if not candidates:
        return None
    candidates.sort(key=lambda r: -r["_score"])
    return candidates[0]


def pick_loops_file(meta: dict) -> "Optional[dict]":
    """Pick a .bedpe loops file from a Hi-C / ChIA-PET experiment."""
    files = meta.get("files") or []
    for f in files:
        if not isinstance(f, dict) or f.get("status") != "released":
            continue
        fmt = (f.get("file_format") or "").lower()
        otype = (f.get("output_type") or "").lower()
        if "bedpe" in fmt or "loops" in otype or "interactions" in otype:
            return f
    return None


def download_file(url: str, dest: Path) -> int:
    logging.info("Download %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    written = 0
    with urllib.request.urlopen(req, timeout=180) as resp, \
         dest.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk); written += len(chunk)
    return written


# --------------------------- Linkage methods -------------------------------

def link_via_loops(se_regions: "list[dict]", loops_bedpe: Path
                    ) -> "list[dict]":
    """Use Hi-C / ChIA-PET loops to link SEs to target gene candidates.

    For each SE, find loops where one anchor overlaps the SE; the other
    anchor's nearest gene is a target candidate. Returns rows with
    method='loops' and a 'distance' field (loop-anchor distance).
    """
    loops = ep._parse_bedpe(loops_bedpe)
    out: "list[dict]" = []
    for se in se_regions:
        for L in loops:
            # Anchor 1 inside SE?
            a1 = (L["chrom1"] == se["chrom"] and
                  L["end1"] >= se["start"] and L["start1"] <= se["end"])
            a2 = (L["chrom2"] == se["chrom"] and
                  L["end2"] >= se["start"] and L["start2"] <= se["end"])
            if not (a1 or a2):
                continue
            target_chrom = L["chrom2"] if a1 else L["chrom1"]
            target_start = L["start2"] if a1 else L["start1"]
            target_end = L["end2"] if a1 else L["end1"]
            target_mid = (target_start + target_end) // 2
            distance = abs(target_mid - ((se["start"] + se["end"]) // 2))
            out.append({
                "se_id":       se["name"],
                "se_chrom":    se["chrom"],
                "se_start":    se["start"],
                "se_end":      se["end"],
                "target_chrom": target_chrom,
                "target_start": target_start,
                "target_end":   target_end,
                "method":      "loops",
                "loop_score":  L.get("score", 0.0),
                "loop_distance": distance,
                "_anchor":     "anchor1" if a1 else "anchor2",
            })
    return out


def link_via_catalog(se_regions: "list[dict]", limit_per_se: int = 25,
                        timeout_per_call: int = 20,
                        max_consecutive_failures: int = 3,
                        max_ses: int = 100) -> "list[dict]":
    """Query the IGVF Catalog regulatory-regions/genes endpoint per SE
    region. Returns enhancer-gene predictions (rE2G / ABC / scE2G).

    Hard-caps:
      - ``timeout_per_call``: per-HTTP-call timeout (sec). Catalog can
        be slow / occasionally hang; we never block longer than this
        per call.
      - ``max_consecutive_failures``: if N catalog calls in a row time
        out or error, give up the catalog step entirely (proximity +
        loops still produce target candidates).
      - ``max_ses``: cap on number of SEs queried — protects against
        hour-long runs when there are thousands of SEs.
    """
    out: "list[dict]" = []
    consecutive_failures = 0
    queried = 0
    aborted = False
    for i, se in enumerate(se_regions[:max_ses]):
        region = f"{se['chrom']}:{se['start']}-{se['end']}"
        got_for_this_se = False
        for path in ("/api/regulatory-regions/genes",
                       "/api/genomic-elements/genes"):
            t0 = time.time()
            sc, data = catalog_get(path, region=region,
                                          limit=limit_per_se,
                                          timeout=timeout_per_call)
            elapsed = time.time() - t0
            queried += 1
            if sc != 200:
                consecutive_failures += 1
                logging.info("catalog SE %d/%d %s -> HTTP %s after %.1fs "
                              "(consecutive failures: %d)",
                              i + 1, min(len(se_regions), max_ses),
                              path, sc, elapsed, consecutive_failures)
                if consecutive_failures >= max_consecutive_failures:
                    logging.warning(
                        "Aborting catalog linkage after %d consecutive "
                        "failures; falling back to loops + proximity only.",
                        consecutive_failures,
                    )
                    aborted = True
                    break
                continue
            consecutive_failures = 0
            rows = (data if isinstance(data, list)
                    else (data.get("results") or data.get("@graph") or []))
            for r in rows:
                if not isinstance(r, dict):
                    continue
                gene_path = r.get("gene") or r.get("target_gene") or ""
                gene_id = (gene_path.split("/", 1)[1]
                            if isinstance(gene_path, str) and "/" in gene_path
                            else gene_path)
                out.append({
                    "se_id":      se["name"],
                    "se_chrom":   se["chrom"],
                    "se_start":   se["start"],
                    "se_end":     se["end"],
                    "gene_id":    gene_id,
                    "method":     r.get("method", path),
                    "score":      r.get("score", ""),
                    "biological_context":
                        r.get("biological_context", ""),
                    "source":     r.get("source", ""),
                })
                got_for_this_se = True
            logging.info("catalog SE %d/%d %s -> %d hits in %.1fs",
                          i + 1, min(len(se_regions), max_ses),
                          path, len(rows), elapsed)
            if got_for_this_se:
                break  # don't try the alternate path if we already got hits
            time.sleep(0.05)
        if aborted:
            break
    logging.info("catalog linkage: queried %d, got %d gene rows%s",
                  queried, len(out), " (aborted early)" if aborted else "")
    return out


def link_via_proximity(se_regions: "list[dict]",
                          window_kb: int = 500,
                          per_se_limit: int = 5) -> "list[dict]":
    """Catalog gene query within ±window_kb of each SE's midpoint."""
    out: "list[dict]" = []
    half = window_kb * 1000
    for se in se_regions:
        mid = (se["start"] + se["end"]) // 2
        region = f"{se['chrom']}:{max(mid - half, 0)}-{mid + half}"
        sc, data = catalog_get("/api/genes", region=region,
                                  limit=per_se_limit * 5)
        if sc != 200:
            continue
        rows = (data if isinstance(data, list)
                else (data.get("results") or data.get("@graph") or []))
        # Score by distance to TSS
        cands = []
        for g in rows[: per_se_limit * 5]:
            if not isinstance(g, dict):
                continue
            tss = g.get("start") or g.get("start_position") or 0
            try:
                dist = abs(int(tss) - mid)
            except (TypeError, ValueError):
                dist = half
            cands.append((dist, g))
        cands.sort(key=lambda x: x[0])
        for dist, g in cands[:per_se_limit]:
            out.append({
                "se_id":   se["name"],
                "se_chrom": se["chrom"],
                "se_start": se["start"],
                "se_end":  se["end"],
                "gene_id": g.get("name") or g.get("symbol", ""),
                "method":  "proximity",
                "distance_to_tss": dist,
                "score":   max(0.0, 1.0 - dist / half),
            })
        time.sleep(0.05)
    return out


def link_via_ccre_classes(se_regions: "list[dict]",
                              ccre_bed: Path) -> "list[dict]":
    """For each SE, find which cCRE classes (PLS / pELS / dELS / CTCF)
    fall inside it. Reports per-SE composition (not per-gene) but is
    useful evidence that a region is enhancer-like (mostly ELS) vs
    promoter-bound (mostly PLS)."""
    annotated, _ = ep.overlap_peaks_with_ccres(
        [{**se, "name": se["name"]} for se in se_regions], ccre_bed,
    )
    by_se: "dict[str, dict]" = {}
    for r in annotated:
        sid = r["name"]
        bucket = by_se.setdefault(sid, {"counts": Counter()})
        bucket["counts"][r.get("ccre_class", "none")] += 1
    out: "list[dict]" = []
    for se in se_regions:
        bucket = by_se.get(se["name"], {"counts": Counter()})
        out.append({
            "se_id":   se["name"],
            "se_chrom": se["chrom"],
            "se_start": se["start"],
            "se_end":  se["end"],
            "ccre_classes":
                ", ".join(f"{k}:{v}" for k, v in bucket["counts"].items()),
            "n_pls":  bucket["counts"].get("PLS", 0),
            "n_pels": bucket["counts"].get("pELS", 0),
            "n_dels": bucket["counts"].get("dELS", 0),
            "n_ctcf": bucket["counts"].get("CTCF-only", 0),
        })
    return out


# --------------------------- Aggregation -----------------------------------

def aggregate_se_targets(loops_links: "list[dict]",
                            catalog_links: "list[dict]",
                            proximity_links: "list[dict]"
                            ) -> "list[dict]":
    """Combine the three gene-level evidence streams into a per-(SE, gene)
    table with method counts and a composite support score."""
    by_pair: "dict[tuple[str, str], dict]" = {}

    for r in catalog_links:
        key = (r["se_id"], r.get("gene_id", ""))
        if not key[1]:
            continue
        rec = by_pair.setdefault(key, {
            "se_id": r["se_id"], "gene_id": key[1],
            "se_chrom": r["se_chrom"], "se_start": r["se_start"],
            "se_end": r["se_end"], "methods": set(),
            "loop_score": 0.0, "catalog_method": "",
            "biological_context": "", "distance_to_tss": None,
        })
        rec["methods"].add(f"catalog:{r.get('method','')}")
        rec["catalog_method"] = r.get("method", "")
        if r.get("biological_context"):
            rec["biological_context"] = r["biological_context"]

    for r in proximity_links:
        key = (r["se_id"], r.get("gene_id", ""))
        if not key[1]:
            continue
        rec = by_pair.setdefault(key, {
            "se_id": r["se_id"], "gene_id": key[1],
            "se_chrom": r["se_chrom"], "se_start": r["se_start"],
            "se_end": r["se_end"], "methods": set(),
            "loop_score": 0.0, "catalog_method": "",
            "biological_context": "", "distance_to_tss": None,
        })
        rec["methods"].add("proximity")
        rec["distance_to_tss"] = r.get("distance_to_tss")

    # Loops link SEs to anchor regions; we don't always know the gene at the
    # other anchor without a TSS BED. Carry them as separate evidence
    # rows and join to (SE, gene) pairs whose proximity / catalog gene is
    # close to the loop's other anchor.
    loop_anchor_by_se: "dict[str, list[dict]]" = defaultdict(list)
    for r in loops_links:
        loop_anchor_by_se[r["se_id"]].append(r)

    rows: "list[dict]" = []
    for (se_id, gene), rec in by_pair.items():
        # Check if any loop anchor is near (within 50kb) this gene's
        # implied position (proxied by the SE region for now).
        loop_support = False
        loop_score = 0.0
        for L in loop_anchor_by_se.get(se_id, []):
            # Without a TSS BED we can't be precise; treat any loop hit
            # on the SE as weak supporting evidence for ALL its target
            # genes from other methods.
            loop_support = True
            loop_score = max(loop_score, float(L.get("loop_score", 0)))
        if loop_support:
            rec["methods"].add("loops")
            rec["loop_score"] = loop_score

        rec["methods"] = sorted(rec["methods"])
        rec["n_methods"] = len(rec["methods"])
        rec["support_score"] = (
            (3.0 if any(m.startswith("catalog:") for m in rec["methods"])
                 else 0.0)
            + (2.0 if "loops" in rec["methods"] else 0.0)
            + (1.0 if "proximity" in rec["methods"] else 0.0)
        )
        rows.append(rec)
    rows.sort(key=lambda r: (-r["support_score"], r["se_id"], r["gene_id"]))
    return rows


# --------------------------- Plots -----------------------------------------

def plot_se_target_network(targets: "list[dict]", out_path: Path,
                              top_n: int = 60) -> None:
    """Bipartite SE↔gene network. SEs left, genes right, edge weight =
    support score."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    data = [r for r in targets if r.get("gene_id")][:top_n]
    if not data:
        return

    ses = sorted({r["se_id"] for r in data})
    genes = sorted({r["gene_id"] for r in data})
    se_y = {s: i for i, s in enumerate(ses)}
    gene_y = {g: i for i, g in enumerate(genes)}

    fig_h = max(4.0, 0.22 * max(len(ses), len(genes)))
    fig, ax = plt.subplots(figsize=(9, fig_h))
    for r in data:
        ax.plot([0, 1], [se_y[r["se_id"]], gene_y[r["gene_id"]]],
                  color={3: "#D55E00", 4: "#D55E00", 5: "#D55E00",
                         6: "#D55E00", 2: "#0072B2",
                         1: "#999999", 0: "#cccccc"}.get(
                              int(r.get("support_score", 0)), "#888"),
                  alpha=min(0.95, 0.3 + 0.15 * r.get("support_score", 0)),
                  lw=0.7)
    for s, y in se_y.items():
        ax.text(-0.02, y, s, ha="right", va="center", fontsize=7)
    for g, y in gene_y.items():
        ax.text(1.02, y, g, ha="left", va="center", fontsize=7)
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylim(-0.5, max(len(ses), len(genes)) - 0.5)
    ax.invert_yaxis()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Super-enhancers", "Target genes"])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"SE → target network (top {len(data)} edges)")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight",
                  facecolor="white")
    plt.close(fig)


# --------------------------- Subcommands -----------------------------------

def cmd_discover(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    label = safe_label(args.label or args.biosample)

    chip = discover_chipseq(args.biosample, target=args.target,
                              assembly=args.assembly, limit=args.limit)
    access = discover_chromatin_accessibility(args.biosample,
                                                  assembly=args.assembly,
                                                  limit=args.limit)
    chrm3d = discover_3d(args.biosample, assembly=args.assembly,
                            limit=args.limit)

    def _summarize(rows, label):
        return [{
            "category":     label,
            "accession":    r.get("accession", ""),
            "assay_title":  r.get("_assay_title") or r.get("assay_title", ""),
            "target":       (r.get("target") or {}).get("label", "")
                              if isinstance(r.get("target"), dict) else "",
            "biosample":    (r.get("biosample_ontology") or
                              {}).get("term_name", "")
                              if isinstance(r.get("biosample_ontology"), dict)
                              else "",
            "lab":          (r.get("lab") or {}).get("title", "")
                              if isinstance(r.get("lab"), dict) else "",
            "url":          ENCODE_BASE + (r.get("@id") or ""),
            "description":  (r.get("description") or "")[:160],
        } for r in rows]

    summary = (_summarize(chip, f"ChIP-seq ({args.target})")
               + _summarize(access, "DNase/ATAC")
               + _summarize(chrm3d, "3D chromatin"))
    out = MANIFEST_DIR / f"{ts}_{label}_discovery.csv"
    write_csv(out, summary,
              cols=["category", "accession", "assay_title", "target",
                     "biosample", "lab", "url", "description"])
    print(f"Discovered {len(chip)} ChIP / {len(access)} DNase-ATAC / "
          f"{len(chrm3d)} 3D experiments for biosample={args.biosample}")
    print(f"Manifest: {out}")
    return out


def cmd_pipeline(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    biosample = args.biosample
    label = safe_label(args.label or f"{biosample}_{args.target}_se2target")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    plot_dir = out_dir / "Plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ----- 1. Discovery -------------------------------------------------
    print(f"▶ Step 1/5  Discovering experiments for {biosample}")
    chip = discover_chipseq(biosample, target=args.target,
                              assembly=args.assembly, limit=args.limit)
    if not chip:
        raise SystemExit(
            f"No {args.target} ChIP-seq experiments found for "
            f"biosample={biosample!r} on {args.assembly}. "
            f"Try `igvfagent encode retrieve --assay 'Histone ChIP-seq' "
            f"--target {args.target} --biosample {biosample}` to "
            f"explore."
        )
    chrm3d = discover_3d(biosample, assembly=args.assembly,
                            limit=args.limit) if args.include_3d else []

    # Pick the best experiment (most files, replicated, released)
    chip_hydrated = []
    for h in chip[: min(args.experiment_limit, len(chip))]:
        meta = hydrate(h.get("accession", ""))
        chip_hydrated.append(meta)
        time.sleep(0.1)
    chip_hydrated.sort(key=lambda m: -len(m.get("files") or []))
    pick_meta = chip_hydrated[0] if chip_hydrated else None
    if not pick_meta:
        raise SystemExit("Could not hydrate any ChIP-seq experiments.")
    pick_acc = pick_meta.get("accession", "")
    print(f"   → ChIP-seq pick: {pick_acc} "
          f"({len(pick_meta.get('files') or [])} files)")

    # ----- 2. Download peaks BED ---------------------------------------
    print(f"▶ Step 2/5  Downloading peak BED for {pick_acc}")
    peak_file = pick_peak_file(pick_meta)
    if not peak_file:
        raise SystemExit(f"No suitable peaks BED file in {pick_acc}.")
    pf_url = ENCODE_BASE + (peak_file.get("href") or "")
    # Preserve the URL's original filename so .gz/.bgz detection works
    # downstream in parse_bed (which keys on the filename suffix).
    url_filename = (peak_file.get("href") or "").rsplit("/", 1)[-1] \
                    or "peaks.bed.gz"
    pf_dest = DOWNLOAD_DIR / f"{pick_acc}_{url_filename}"
    pf_dest.parent.mkdir(parents=True, exist_ok=True)
    if not pf_dest.exists() or pf_dest.stat().st_size == 0:
        download_file(pf_url, pf_dest)
    print(f"   → peaks BED: {pf_dest} "
          f"({pf_dest.stat().st_size / (1024**2):.1f} MB)")

    # ----- 3. Call super-enhancers --------------------------------------
    print(f"▶ Step 3/5  Calling super-enhancers")
    rows, _ = ep.parse_bed(pf_dest)
    if not rows:
        raise SystemExit("No peaks parsed from downloaded BED.")
    stitched = ep.stitch_peaks(
        rows, stitching_distance=args.stitching_distance,
    )
    stitched.sort(key=lambda r: r["score_sum"], reverse=True)
    n_se = ep.find_inflection([r["score_sum"] for r in stitched])
    se_regions: "list[dict]" = []
    for i, r in enumerate(stitched[:n_se]):
        se_regions.append({
            "chrom": r["chrom"], "start": r["start"], "end": r["end"],
            "name":  f"SE_{i+1}", "score": r["score_sum"],
            "rank":  i + 1, "n_peaks": r["n_peaks"],
        })
    se_bed = out_dir / f"{label}_super_enhancers.bed"
    with se_bed.open("w") as f:
        for s in se_regions:
            f.write(f"{s['chrom']}\t{s['start']}\t{s['end']}\t"
                     f"{s['name']}\t{s['score']:.3f}\t.\n")
    print(f"   → {len(se_regions)} super-enhancers above the inflection")

    # ----- 4a. cCRE class composition per SE ---------------------------
    print(f"▶ Step 4/5  Linking SEs to target genes")
    print(f"   (a) cCRE class composition")
    ccre_bed = ep._ensure_screen_ccres()
    ccre_summary = link_via_ccre_classes(se_regions, ccre_bed)
    write_csv(out_dir / f"{label}_se_ccre_composition.csv", ccre_summary)

    # ----- 4b. Catalog enhancer-gene linkage ---------------------------
    print(f"   (b) IGVF Catalog rE2G / ABC predictions")
    catalog_links = link_via_catalog(se_regions,
                                          limit_per_se=args.per_se_genes)
    print(f"       {len(catalog_links)} catalog SE→gene rows")

    # ----- 4c. Loop overlap (if 3D experiments + bedpe available) ------
    loops_links: "list[dict]" = []
    loops_dest: Optional[Path] = None
    if args.include_3d and chrm3d:
        for cm in chrm3d:
            meta3d = hydrate(cm.get("accession", ""))
            lf = pick_loops_file(meta3d)
            if lf:
                lf_url = ENCODE_BASE + (lf.get("href") or "")
                # Preserve the URL filename so .gz / .bedpe / .bedpe.gz
                # is detected correctly when parsing.
                lf_filename = (lf.get("href") or "").rsplit("/", 1)[-1] \
                                or "loops.bedpe"
                loops_dest = DOWNLOAD_DIR / (
                    f"{cm.get('accession','3d')}_{lf_filename}"
                )
                if not loops_dest.exists() or loops_dest.stat().st_size == 0:
                    try:
                        download_file(lf_url, loops_dest)
                    except Exception as e:
                        logging.warning("loop download failed: %s", e)
                        loops_dest = None
                if loops_dest:
                    print(f"   (c) Loops file: {loops_dest.name}")
                    loops_links = link_via_loops(se_regions, loops_dest)
                    print(f"       {len(loops_links)} loop-based SE↔anchor "
                          f"associations")
                    break
            time.sleep(0.05)

    # ----- 4d. Proximity fallback --------------------------------------
    print(f"   (d) Proximity (±{args.proximity_kb} kb) fallback")
    proximity_links = link_via_proximity(se_regions,
                                              window_kb=args.proximity_kb,
                                              per_se_limit=args.per_se_genes)
    print(f"       {len(proximity_links)} proximity SE→gene rows")

    # ----- 4e. Aggregate -----------------------------------------------
    targets = aggregate_se_targets(loops_links, catalog_links, proximity_links)
    write_csv(out_dir / f"{label}_se_targets.csv", targets)
    write_csv(out_dir / f"{label}_se_targets_loops.csv", loops_links)
    write_csv(out_dir / f"{label}_se_targets_catalog.csv", catalog_links)
    write_csv(out_dir / f"{label}_se_targets_proximity.csv", proximity_links)

    # ----- 5. Plots + report -------------------------------------------
    print(f"▶ Step 5/5  Plots + report")
    plot_se_target_network(targets, plot_dir / f"{label}_se_target_network.png",
                              top_n=args.top_n_edges)

    if args.gene:
        # Per-gene browser SVG: pick the best SE for this gene
        gene_rows = [r for r in targets if r["gene_id"].upper() == args.gene.upper()]
        if gene_rows:
            best = gene_rows[0]
            region = (f"{best['se_chrom']}:"
                       f"{max(best['se_start'] - 50000, 0)}-"
                       f"{best['se_end'] + 50000}")
            from io import StringIO  # noqa
            tracks = [f"super_enhancers:{se_bed}"]
            if loops_dest:
                # Render loops as a track too
                tracks.append(f"loops:{loops_dest}")
            try:
                ns = argparse.Namespace(
                    region=region, track=tracks, with_ccre=True,
                    ccre_bed=None, width=1000, label=f"{label}_{args.gene}",
                )
                ep.cmd_browser(ns)
            except Exception as e:
                logging.warning("Browser SVG failed: %s", e)

    report = out_dir / f"{label}_report.md"
    n_pairs = len(targets)
    n_genes = len({r["gene_id"] for r in targets if r.get("gene_id")})
    n_high = sum(1 for r in targets if r.get("support_score", 0) >= 4)
    lines = [
        f"# Super-Enhancer → Target-Gene pipeline",
        f"## Biosample: `{biosample}` · Target: `{args.target}` · Assembly: `{args.assembly}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Inputs",
        "",
        f"- ChIP-seq experiment: `{pick_acc}`  "
            f"_(picked from {len(chip)} candidates)_",
        f"- Peaks BED: `{pf_dest.relative_to(ROOT)}`",
        (f"- 3D-chromatin experiment: `{cm.get('accession','')}` "
         f"(loops: `{loops_dest.name if loops_dest else 'none'}`)"
         if loops_links else "- 3D chromatin: _not used_"),
        f"- SCREEN cCRE BED: `{ccre_bed.relative_to(ROOT)}`",
        "",
        "## Super-enhancer call",
        "",
        f"- Stitching distance: `{args.stitching_distance}` bp",
        f"- Stitched regions ranked by signal: **{len(stitched):,}**",
        f"- Above the inflection (called as super-enhancers): "
        f"**{len(se_regions):,}**",
        f"- BED: `{se_bed.relative_to(ROOT)}`",
        "",
        "## Top 20 super-enhancers",
        "",
        "| rank | region | score_sum | n_peaks |",
        "|---:|---|---:|---:|",
    ]
    for s in se_regions[:20]:
        lines.append(
            f"| {s['rank']} | {s['chrom']}:{s['start']}-{s['end']} | "
            f"{s['score']:.2f} | {s['n_peaks']} |"
        )

    lines += ["",
               "## SE → target-gene linkage",
               "",
               f"- Catalog rE2G/ABC predictions: **{len(catalog_links)}**",
               f"- Loop-based anchor associations: **{len(loops_links)}** "
               f"_(via {loops_dest.name if loops_dest else 'no loops'})_",
               f"- Proximity (±{args.proximity_kb} kb): "
               f"**{len(proximity_links)}**",
               f"- Aggregated unique (SE, gene) pairs: **{n_pairs:,}** "
               f"covering **{n_genes:,}** distinct genes",
               f"- High-confidence pairs (≥ 4 support): **{n_high:,}**",
               ""]

    lines += ["## Top 30 (SE, gene) pairs by support",
              "",
              "| SE | gene | methods | support | catalog method | "
              "biological context |",
              "|---|---|---|---:|---|---|"]
    for r in targets[:30]:
        lines.append(
            f"| {r['se_id']} ({r['se_chrom']}:{r['se_start']}-{r['se_end']}) "
            f"| {r['gene_id']} | {', '.join(r.get('methods',[]))} "
            f"| {r.get('support_score', 0):.0f} "
            f"| {r.get('catalog_method','')} "
            f"| {r.get('biological_context','')} |"
        )

    lines += ["",
               "## Outputs",
               "",
               f"- Aggregated target table: "
                  f"`{(out_dir/f'{label}_se_targets.csv').relative_to(ROOT)}`",
               f"- Per-method tables: "
                  f"`{label}_se_targets_loops.csv`, "
                  f"`{label}_se_targets_catalog.csv`, "
                  f"`{label}_se_targets_proximity.csv`",
               f"- cCRE class composition per SE: "
                  f"`{label}_se_ccre_composition.csv`",
               f"- Network plot: `Plots/{label}_se_target_network.png`",
               "",
               "## Suggested follow-ups",
               "",
               "- For a specific candidate gene, drill in with "
                  "`igvfagent kg gene <SYMBOL> --depth 2 --call-literature`.",
               "- For locus-level visualization, "
                  "`igvfagent encode browser --region "
                  "chrN:start-end --track 'SE:<bed>' --with-ccre`.",
               "- For literature corroboration of the high-support SE→gene "
                  "list, "
                  "`igvfagent ref validate --input <gene_list.csv> "
                  f"--context '{biosample} super-enhancer'`.",
               ]
    report.write_text("\n".join(lines))
    print(f"\n✓ Run dir:    {out_dir}")
    print(f"  Report:     {report}")
    print(f"  SE targets: {(out_dir / f'{label}_se_targets.csv')}")
    return report


def cmd_link_targets(args: argparse.Namespace) -> Path:
    """Standalone linkage step. Inputs: existing SE BED (any source) +
    optional loops bedpe + optional cCRE BED. Useful for re-running just
    the linkage layer on a custom SE call from outside this pipeline."""
    setup_logging(); mkdirs()
    ts = timestamp()
    label = safe_label(args.label or "link_targets")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    se_rows, _ = ep.parse_bed(Path(args.se_bed))
    se_regions = [{
        "chrom": r["chrom"], "start": r["start"], "end": r["end"],
        "name":  r.get("name", f"SE_{i+1}"),
        "score": r.get("score", 0.0), "rank": i + 1,
        "n_peaks": 0,
    } for i, r in enumerate(se_rows)]

    catalog = link_via_catalog(se_regions, limit_per_se=args.per_se_genes)
    proximity = link_via_proximity(se_regions, window_kb=args.proximity_kb,
                                       per_se_limit=args.per_se_genes)
    loops = (link_via_loops(se_regions, Path(args.loops_bedpe))
              if args.loops_bedpe else [])

    aggregated = aggregate_se_targets(loops, catalog, proximity)
    write_csv(out_dir / f"{label}_se_targets.csv", aggregated)
    write_csv(out_dir / f"{label}_loops.csv", loops)
    write_csv(out_dir / f"{label}_catalog.csv", catalog)
    write_csv(out_dir / f"{label}_proximity.csv", proximity)
    plot_se_target_network(aggregated,
                              out_dir / f"{label}_se_target_network.png")

    print(f"Run dir: {out_dir}")
    print(f"Aggregated rows: {len(aggregated)}")
    return out_dir / f"{label}_se_targets.csv"


def cmd_write_playbook(_a) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "SE_TARGET_PIPELINE_SKILLS.md"
    lines = [
        "# Skill: Super-enhancer → target-gene pipeline",
        "",
        "End-to-end workflow that takes any ENCODE-supported cell line / "
        "tissue (GM12878, K562, HepG2, liver, hippocampus, cortex, …), "
        "discovers H3K27ac (or BRD4 / MED1 / P300) ChIP-seq plus optional "
        "3D-chromatin (Hi-C / ChIA-PET / capture Hi-C) experiments, "
        "downloads the peak files, calls super-enhancers ROSE-style, and "
        "links each super-enhancer to candidate target genes.",
        "",
        "Linkage combines four evidence streams per (SE, gene) pair:",
        "",
        "1. **3D loops** — Hi-C / ChIA-PET / capture-Hi-C anchors that "
        "overlap the SE; the other anchor's region implies a target gene.",
        "2. **IGVF Catalog rE2G / ABC predictions** — region-level "
        "regulatory-region → gene maps in the same biosample.",
        "3. **Proximity** — nearest gene TSS within ±N kb of the SE midpoint.",
        "4. **SCREEN cCRE classes** — per-SE composition (PLS / pELS / "
        "dELS / CTCF) as quality / specificity evidence (not direct "
        "linkage but useful for filtering).",
        "",
        "Each (SE, gene) pair is scored by how many of those streams "
        "support it. The top pairs are written to a ranked CSV and "
        "rendered as a bipartite SE↔gene network.",
        "",
        "## Subcommands",
        "",
        "### `pipeline` — one-command workflow",
        "",
        "```bash",
        "# Default: GM12878 H3K27ac, no 3D",
        "igvfagent se-targets pipeline --biosample GM12878",
        "",
        "# Include Hi-C / ChIA-PET 3D chromatin layer + focus on APOE",
        "igvfagent se-targets pipeline \\",
        "    --biosample GM12878 --target H3K27ac --assembly GRCh38 \\",
        "    --include-3d --gene APOE --label gm12878_apoe_run",
        "",
        "# Different cell line, different mark",
        "igvfagent se-targets pipeline --biosample K562 --target BRD4",
        "igvfagent se-targets pipeline --biosample 'liver' --target H3K27ac",
        "```",
        "",
        "Outputs land under `Docs/SETargets/<timestamp>_<label>/`:",
        "",
        "  - `<label>_super_enhancers.bed`",
        "  - `<label>_se_targets.csv`           # aggregated, sorted by support",
        "  - `<label>_se_targets_loops.csv`     # loop-anchor evidence rows",
        "  - `<label>_se_targets_catalog.csv`   # rE2G / ABC predictions",
        "  - `<label>_se_targets_proximity.csv` # nearest-TSS rows",
        "  - `<label>_se_ccre_composition.csv`  # PLS/pELS/dELS counts per SE",
        "  - `Plots/<label>_se_target_network.png`",
        "  - `<label>_report.md`                # human-readable summary",
        "",
        "### `discover` — list candidate experiments",
        "",
        "```bash",
        "igvfagent se-targets discover --biosample GM12878 --target H3K27ac",
        "```",
        "",
        "Writes a manifest of candidate ChIP-seq + DNase/ATAC + 3D-"
        "chromatin experiments for the biosample (no downloads, no SE "
        "calling).",
        "",
        "### `link-targets` — standalone linkage step",
        "",
        "```bash",
        "igvfagent se-targets link-targets \\",
        "    --se-bed path/to/super_enhancers.bed \\",
        "    --loops-bedpe path/to/hic_loops.bedpe",
        "```",
        "",
        "Useful when you've called SEs from a custom upstream pipeline and "
        "just want the IGVFagent linkage layer on top.",
        "",
        "## How this chains with other skills",
        "",
        "- After `pipeline`, drill into any candidate gene with "
        "`igvfagent kg gene <SYMBOL> --depth 2 --call-literature`.",
        "- For locus-level visualization, "
        "`igvfagent encode browser --region chrN:start-end "
        "--track 'SE:<bed>' --with-ccre`.",
        "- For literature corroboration of the high-support gene list, "
        "`igvfagent ref validate --input <gene_list.csv>`.",
        "- The internal Plan → Action → Results → Evaluation orchestrator "
        "exposes `se_targets_pipeline` as a single tool so a natural-"
        "language query like *'For GM12878, call super-enhancers from "
        "H3K27ac and tell me which genes they target via Hi-C'* "
        "fires the whole workflow in one shot.",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI -------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Super-enhancer → target-gene pipeline (ENCODE → SE → "
                    "loops + cCRE + Catalog → ranked target genes).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("pipeline",
                        help="One-command workflow.")
    s.add_argument("--biosample", required=True,
                    help="ENCODE biosample term name (GM12878, K562, "
                         "liver, hippocampus, …).")
    s.add_argument("--target", default=DEFAULT_SE_TARGET,
                    help=f"Enhancer-mark ChIP-seq target. Default "
                         f"{DEFAULT_SE_TARGET}.")
    s.add_argument("--assembly", default="GRCh38")
    s.add_argument("--limit", type=int, default=15,
                    help="Max experiments per assay group to consider.")
    s.add_argument("--experiment-limit", type=int, default=5,
                    help="How many ChIP-seq experiments to hydrate before "
                         "picking the best one.")
    s.add_argument("--include-3d", action="store_true",
                    help="Also pull a Hi-C / ChIA-PET loops file for "
                         "loop-based linkage.")
    s.add_argument("--stitching-distance", type=int, default=12500)
    s.add_argument("--proximity-kb", type=int, default=500,
                    help="Half-window for the proximity-based linkage "
                         "fallback (default ±500 kb).")
    s.add_argument("--per-se-genes", type=int, default=8,
                    help="Cap on per-SE candidate genes per method.")
    s.add_argument("--gene", default=None,
                    help="Optional gene of interest — also renders a "
                         "browser SVG of its top-supported SE.")
    s.add_argument("--top-n-edges", type=int, default=80,
                    help="Network plot edge cap.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("discover",
                        help="List candidate experiments without "
                             "downloading anything.")
    s.add_argument("--biosample", required=True)
    s.add_argument("--target", default=DEFAULT_SE_TARGET)
    s.add_argument("--assembly", default="GRCh38")
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_discover)

    s = sub.add_parser("link-targets",
                        help="Standalone linkage step.")
    s.add_argument("--se-bed", required=True,
                    help="BED of super-enhancers to link.")
    s.add_argument("--loops-bedpe", default=None,
                    help="Optional Hi-C / ChIA-PET / cHi-C loops bedpe.")
    s.add_argument("--proximity-kb", type=int, default=500)
    s.add_argument("--per-se-genes", type=int, default=8)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_link_targets)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/SE_TARGET_PIPELINE_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
