#!/usr/bin/env python3
"""IGVF agent Tabula Sapiens 2.0 skill.

Retrieval and reproduction for the Tabula Sapiens 2.0 human cell atlas
(Tabula Sapiens Consortium, *Cell* 2026) -- 1,136,218 cells across 28
tissues from 24 donors, annotated to 182 fine and 38 broad cell types.

Data routes, and which are actually open
----------------------------------------
| Route | Content | Access |
|---|---|---|
| figshare 27921984 | 28 processed ``.h5ad``, 57 GB | open |
| GEO GSE306755 | count matrices + full cell metadata | open |
| CELLxGENE | 35 per-tissue datasets | open |
| AWS ``czb-tabula-sapiens`` | raw FASTQ/BAM, 103 TB | **listable, NOT readable** |

The S3 bucket is deliberately gated: ``ListBucket`` succeeds but every
``GetObject`` returns 403, for v1 and v2 alike. That is the data
transfer agreement the paper describes ("to preserve the donors' genetic
privacy, we require a data transfer agreement to receive the raw
sequence reads"). ``s3-manifest`` therefore enumerates and sizes the
bucket -- 103.16 TB: 43.2 TB BAM, 32.1 TB FASTQ, 27.6 TB STAR
intermediates, 0.32 TB everything else -- and documents the DTA path,
rather than pretending a download will work.

Subcommands
-----------
  pull-figshare   Fetch the processed per-tissue h5ads (resumable).
  pull-geo        Fetch GSE306755 count matrices + cell metadata.
  pull-cellxgene  List / fetch the CELLxGENE collection.
  s3-manifest     Enumerate and size the raw-data bucket (no download).
  status          What is present locally, against the paper's numbers.
  overview        Figure 1: donor demographics and tissue composition.
  tf-matrix       Gene x cell-type mean expression (input to Figures 2-3).
  tf-specificity  Figure 2: tau statistic, specific vs ubiquitous TFs.
  tf-enrichment   Figure 3: GO enrichment of the non-specific TFs.
  tf-regulons     SCENIC-style regulon inference and AUCell activity.
  senescence      Figure 4: CDKN2A+ MKI67- burden across the atlas.
  sex-de          Figure 5: pseudobulk sex differences per tissue.
  donors          Figure 6: donor demographics (ChatTS, offline core).
  write-playbook

Reproduction scope: this skill starts from the processed deposits, the
same boundary ``share`` and ``spatial-hic`` draw. Read alignment (STAR /
CellRanger) is orchestrated, not reimplemented -- see the playbook.

License: Apache-2.0. Upstream analysis code
(github.com/czbiohub-sf/tabula-sapiens, BSD-3-Clause) is followed from
its published description; no source is copied. Numerics live in
``_ts_stats``; atlas I/O in ``_ts_atlas``.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "SingleCell" / "tabular-sapiens"
FIGSHARE_DIR = DATA_DIR / "figshare_v2"
GEO_DIR = DATA_DIR / "geo"
CXG_DIR = DATA_DIR / "cellxgene"
MANIFEST_DIR = DATA_DIR / "_manifests"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "TabulaSapiens"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

FIGSHARE_ARTICLE = 27921984
GEO_SERIES = "GSE306755"
CXG_COLLECTION = "e5f58829-1a66-40b5-a624-9046778e74f5"
S3_BUCKET = "https://czb-tabula-sapiens.s3.amazonaws.com/"
S3_NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"
S3_PREFIX = "TabulaSapiens_v2/"

sys.path.insert(0, str(Path(__file__).resolve().parent))


def _atlas():
    try:
        from igvfagent import _ts_atlas as A  # type: ignore
    except Exception:
        import _ts_atlas as A  # type: ignore
    return A


def _stats():
    try:
        from igvfagent import _ts_stats as S  # type: ignore
    except Exception:
        import _ts_stats as S  # type: ignore
    return S


def _htf():
    try:
        from igvfagent import humantfs_skill as H  # type: ignore
    except Exception:
        import humantfs_skill as H  # type: ignore
    return H


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    p = LOG_DIR / f"tabula_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(p), logging.StreamHandler(sys.stdout)],
    )
    return p


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in (s or "run"))


def run_dir(label: str) -> Path:
    d = REPORT_DIR / f"{timestamp()}_{safe_label(label)}"
    (d / "Plots").mkdir(parents=True, exist_ok=True)
    return d


def write_json(p: Path, payload: Any) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str),
                 encoding="utf-8")
    return p


def write_tsv(p: Path, rows: "list[dict]", cols: "Optional[list[str]]" = None) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    c = list(cols or (rows[0].keys() if rows else []))
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=c, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def _get(url: str, *, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "IGVFagent/tabula"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _download(url: str, dest: Path, *, expect: Optional[int] = None) -> str:
    """Resumable-ish download with a size check."""
    if dest.exists() and (expect is None or dest.stat().st_size == expect):
        return "have"
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "IGVFagent/tabula"})
    with urllib.request.urlopen(req, timeout=1800) as r, dest.open("wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)
    if expect is not None and dest.stat().st_size != expect:
        return "size-mismatch"
    return "got"


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def cmd_pull_figshare(args: argparse.Namespace) -> int:
    setup_logging()
    meta = json.loads(_get(
        f"https://api.figshare.com/v2/articles/{args.article}").decode())
    files = meta.get("files", [])
    if args.tissues:
        want = {t.strip().lower() for t in args.tissues.split(",") if t.strip()}
        files = [f for f in files
                 if f["name"].split("_TSP")[0].lower() in want]
    files.sort(key=lambda f: f["size"])

    total = sum(f["size"] for f in files)
    print(f"figshare article {args.article}: {meta.get('title')}")
    print(f"{len(files)} file(s), {total/1e9:.1f} GB")
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(MANIFEST_DIR / "figshare_v2_files.tsv",
              [{"name": f["name"], "size": f["size"],
                "url": f["download_url"]} for f in files],
              ["name", "size", "url"])
    if args.list_only:
        for f in files:
            print(f"  {f['name']:62} {f['size']/1e6:8.1f} MB")
        return 0

    ok = skipped = failed = 0
    for f in files:
        dest = FIGSHARE_DIR / f["name"]
        r = _download(f["download_url"], dest, expect=f["size"])
        if r == "have":
            skipped += 1
        elif r == "got":
            ok += 1
            logging.info("fetched %s", f["name"])
        else:
            failed += 1
            logging.warning("%s: %s", f["name"], r)
    print(f"\nDownloaded {ok}, already had {skipped}, failed {failed}")
    print(f"Location: {FIGSHARE_DIR}")
    return 0 if failed == 0 else 1


def cmd_pull_geo(args: argparse.Namespace) -> int:
    """Fetch the GEO deposit -- notably the full 1.1M-cell metadata."""
    setup_logging()
    gse = args.gse
    ftp = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
           f"{gse[:-3]}nnn/{gse}/suppl/")
    listing = _get(ftp).decode("utf8", "replace")
    import re
    names = sorted(set(re.findall(r'href="([^"?/][^"]*)"', listing)))
    names = [n for n in names if not n.startswith("..")]

    print(f"GEO {gse}: {len(names)} supplementary file(s)")
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    rows = [{"name": n, "url": ftp + n} for n in names]
    write_tsv(MANIFEST_DIR / f"{gse}_files.tsv", rows, ["name", "url"])
    for n in names:
        print(f"  {n}")

    if args.list_only:
        return 0
    pattern = args.download or "metadata"
    rx = re.compile(pattern)
    got = 0
    for n in names:
        if not rx.search(n):
            continue
        dest = GEO_DIR / n
        r = _download(ftp + n, dest)
        print(f"  {r:4} {n}")
        got += (r != "size-mismatch")
    print(f"\nFetched {got} file(s) matching /{pattern}/ into {GEO_DIR}")
    if not args.download:
        print("Pass --download '<regex>' to fetch more (e.g. 'RAW.tar').")
    return 0


def cmd_pull_cellxgene(args: argparse.Namespace) -> int:
    setup_logging()
    url = ("https://api.cellxgene.cziscience.com/curation/v1/collections/"
           f"{args.collection}")
    meta = json.loads(_get(url).decode())
    ds = meta.get("datasets", [])
    print(f"CELLxGENE collection: {meta.get('name')}")
    print(f"{len(ds)} dataset(s), "
          f"{sum(d.get('cell_count') or 0 for d in ds):,} cells")

    rows = []
    for d in ds:
        h5 = next((a for a in (d.get("assets") or [])
                   if a.get("filetype", "").upper() == "H5AD"), None)
        rows.append({"title": d.get("title"), "dataset_id": d.get("dataset_id"),
                     "cells": d.get("cell_count"),
                     "url": (h5 or {}).get("url", "")})
    rows.sort(key=lambda r: str(r["title"]))
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(MANIFEST_DIR / "cellxgene_datasets.tsv", rows,
              ["title", "dataset_id", "cells", "url"])
    for r in rows[: args.limit]:
        print(f"  {str(r['title'])[:44]:46} {str(r['cells']):>9}")

    if args.list_only or not args.tissues:
        if not args.list_only:
            print("\nPass --tissues <comma list> to download, "
                  "or --list-only to just enumerate.")
        return 0
    want = {t.strip().lower() for t in args.tissues.split(",") if t.strip()}
    got = 0
    for r in rows:
        title = str(r["title"])
        tag = title.split("-")[-1].strip().lower()
        if tag not in want and title.lower() not in want:
            continue
        if not r["url"]:
            logging.warning("%s: no H5AD asset", title)
            continue
        dest = CXG_DIR / (safe_label(title) + ".h5ad")
        print(f"  {_download(r['url'], dest):4} {dest.name}")
        got += 1
    print(f"\nFetched {got} dataset(s) into {CXG_DIR}")
    return 0


def _s3_list(prefix: str, delim: Optional[str] = None):
    token = None
    while True:
        u = f"{S3_BUCKET}?list-type=2&prefix={urllib.parse.quote(prefix)}&max-keys=1000"
        if delim:
            u += f"&delimiter={urllib.parse.quote(delim)}"
        if token:
            u += "&continuation-token=" + urllib.parse.quote(token)
        r = ET.fromstring(_get(u, timeout=300))
        for c in r.findall(S3_NS + "Contents"):
            yield (c.find(S3_NS + "Key").text,
                   int(c.find(S3_NS + "Size").text))
        t = r.find(S3_NS + "NextContinuationToken")
        if t is None or t.text is None:
            return
        token = t.text


def cmd_s3_manifest(args: argparse.Namespace) -> int:
    """Enumerate and size the raw bucket. Deliberately does not download.

    ``GetObject`` is 403 for anonymous callers -- the DTA gate -- so a
    manifest is the honest deliverable. It is also what you need to
    request specific objects once a DTA is in place.
    """
    setup_logging()
    out = run_dir(args.label)

    # Prove the access boundary rather than assuming it.
    probe = None
    for key, _size in _s3_list(f"{S3_PREFIX}TSP1/", None):
        probe = key
        break
    readable = False
    if probe:
        try:
            _get(S3_BUCKET + urllib.parse.quote(probe), timeout=60)
            readable = True
        except Exception:
            readable = False

    donors = args.donors.split(",") if args.donors else None
    prefixes = []
    if donors:
        prefixes = [f"{S3_PREFIX}{d.strip()}/" for d in donors if d.strip()]
    else:
        for key, _ in _s3_list(S3_PREFIX, "/"):
            pass
        # CommonPrefixes need a separate parse; just enumerate the root.
        u = (f"{S3_BUCKET}?list-type=2&prefix={urllib.parse.quote(S3_PREFIX)}"
             f"&delimiter=%2F&max-keys=1000")
        r = ET.fromstring(_get(u, timeout=300))
        prefixes = [p.find(S3_NS + "Prefix").text
                    for p in r.findall(S3_NS + "CommonPrefixes")]

    rows = []
    tot = {"fastq": 0, "bam": 0, "star_tmp": 0, "other": 0}
    for pref in prefixes:
        cat = {"fastq": 0, "bam": 0, "star_tmp": 0, "other": 0}
        n = 0
        for key, size in _s3_list(pref):
            n += 1
            if "/fastqs/" in key or key.endswith(".fastq.gz"):
                cat["fastq"] += size
            elif key.endswith((".bam", ".bai")):
                cat["bam"] += size
            elif key.endswith((".out", ".tab", ".mate1", ".mate2", ".log")):
                cat["star_tmp"] += size
            else:
                cat["other"] += size
        for k in tot:
            tot[k] += cat[k]
        rows.append({"prefix": pref, "files": n,
                     **{k: round(v / 1e9, 2) for k, v in cat.items()}})
        print(f"  {pref:34} {n:>8,} files  "
              f"{sum(cat.values())/1e12:6.2f} TB", flush=True)

    tsv = write_tsv(out / "s3_inventory.tsv", rows,
                    ["prefix", "files", "fastq", "bam", "star_tmp", "other"])
    summary = {
        "bucket": S3_BUCKET, "prefix": S3_PREFIX,
        "anonymous_get_readable": readable,
        "probe_key": probe,
        "total_TB": round(sum(tot.values()) / 1e12, 2),
        "by_category_TB": {k: round(v / 1e12, 2) for k, v in tot.items()},
        "note": ("ListBucket succeeds but GetObject returns 403: the data "
                 "transfer agreement described in the paper. Request access "
                 "via the Tabula Sapiens data portal, then re-run with AWS "
                 "credentials configured."),
    }
    js = write_json(out / "s3_summary.json", summary)

    print(f"\nTotal: {summary['total_TB']:.2f} TB")
    for k, v in sorted(summary["by_category_TB"].items(), key=lambda x: -x[1]):
        print(f"  {k:10} {v:8.2f} TB")
    print(f"\nAnonymous GetObject readable: {readable}")
    if not readable:
        print("  -> raw reads need a signed data transfer agreement. "
              "The open routes are figshare (57 GB processed h5ads), "
              "GEO GSE306755, and CELLxGENE.")
    print(f"Inventory: {tsv}")
    print(f"Summary:   {js}")
    return 0


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

def cmd_status(_args: argparse.Namespace) -> int:
    A = _atlas()
    print("Tabula Sapiens 2.0 — local state\n")
    n_h5ad = len(list(FIGSHARE_DIR.glob("*.h5ad"))) if FIGSHARE_DIR.is_dir() else 0
    size = sum(f.stat().st_size for f in FIGSHARE_DIR.glob("*.h5ad")) \
        if n_h5ad else 0
    print(f"  figshare h5ads   {n_h5ad:>3} / 28   {size/1e9:6.1f} GB   {FIGSHARE_DIR}")
    meta = sorted(GEO_DIR.glob("*metadata*.csv*")) if GEO_DIR.is_dir() else []
    print(f"  GEO metadata     {'yes' if meta else 'no ':>3}         "
          f"{(meta[0].stat().st_size/1e6 if meta else 0):6.1f} MB   {GEO_DIR}")
    n_cxg = len(list(CXG_DIR.glob("*.h5ad"))) if CXG_DIR.is_dir() else 0
    print(f"  CELLxGENE h5ads  {n_cxg:>3} / 35")
    H = _htf()
    print(f"  Human TF DB      {'built' if H.DB_PATH.exists() else 'MISSING'}"
          f"        {len(H.tf_symbols()):>5} TFs")

    if meta:
        import pandas as pd
        df = A.cell_metadata()
        P = A.PAPER
        print("\n  Cell metadata vs the paper:")
        checks = [
            ("cells", len(df), P["cells_total"]),
            ("donors", df["donor"].nunique(), P["donors"]),
            ("tissues", df["tissue"].nunique(), P["tissues"]),
            ("fine cell types", df["cell_ontology_class"].nunique(),
             P["fine_cell_types"]),
            ("populations",
             df.groupby(["tissue", "cell_ontology_class"]).ngroups,
             P["populations"]),
        ]
        if "method" in df.columns:
            dr = df[df["method"] == "10X"]
            checks += [
                ("droplet cells", int((df["method"] == "10X").sum()),
                 P["cells_droplet"]),
                ("FACS cells", int((df["method"] == "smartseq").sum()),
                 P["cells_facs"]),
                ("droplet fine types", dr["cell_ontology_class"].nunique(),
                 P["fine_cell_types_droplet"]),
                ("droplet broad types", dr["broad_cell_class"].nunique(),
                 P["broad_cell_types_droplet"]),
            ]
        for name, got, exp in checks:
            mark = "OK " if got == exp else "!! "
            print(f"    {mark}{name:22} {got:>10,}   paper {exp:>10,}")
    else:
        print("\n  (fetch the metadata for the paper-number check: "
              "igvfagent tabula pull-geo)")
    return 0


# ---------------------------------------------------------------------------
# Figure 1 — overview
# ---------------------------------------------------------------------------

def cmd_overview(args: argparse.Namespace) -> int:
    """Figure 1: donor demographics and per-tissue composition."""
    setup_logging()
    A = _atlas()
    import numpy as np
    import pandas as pd

    df = A.cell_metadata(args.metadata)
    out = run_dir(args.label)

    donors = (df.groupby("donor")
                .agg(cells=("donor", "size"),
                     tissues=("tissue", "nunique"),
                     age=("age", "first"), sex=("sex", "first"),
                     ethnicity=("ethnicity", "first"))
                .reset_index()
                .sort_values("cells", ascending=False))
    donors["age"] = pd.to_numeric(donors["age"], errors="coerce")
    donors["age_group"] = pd.cut(donors["age"], [0, 39, 59, 200],
                                 labels=["<40", "40-59", ">=60"])
    d_tsv = write_tsv(out / "donors.tsv", donors.to_dict("records"))

    tis = (df.groupby("tissue")
             .agg(cells=("tissue", "size"), donors=("donor", "nunique"),
                  cell_types=("cell_ontology_class", "nunique"))
             .reset_index().sort_values("cells", ascending=False))
    t_tsv = write_tsv(out / "tissues.tsv", tis.to_dict("records"))

    pops = (df.groupby(["tissue", "cell_ontology_class"])
              .size().reset_index(name="cells"))
    p_tsv = write_tsv(out / "populations.tsv", pops.to_dict("records"))

    ages = donors["age"].dropna()
    summary = {
        "cells": int(len(df)),
        "donors": int(df["donor"].nunique()),
        "tissues": int(df["tissue"].nunique()),
        "fine_cell_types": int(df["cell_ontology_class"].nunique()),
        "broad_cell_types": int(df["broad_cell_class"].nunique())
        if "broad_cell_class" in df else None,
        "populations": int(len(pops)),
        "age_min": float(ages.min()) if len(ages) else None,
        "age_max": float(ages.max()) if len(ages) else None,
        "age_groups": {str(k): int(v) for k, v in
                       donors["age_group"].value_counts().items()},
        "sex": {str(k): int(v) for k, v in
                donors["sex"].value_counts().items()},
        "by_method": ({str(k): int(v) for k, v in
                       df["method"].value_counts().items()}
                      if "method" in df else {}),
        "paper": A.PAPER,
    }
    js = write_json(out / "overview_summary.json", summary)

    if not args.no_figure:
        _fig_overview(df, donors, tis, out / "Plots")

    print(f"Output dir:   {out}")
    print(f"Donors:       {d_tsv}")
    print(f"Tissues:      {t_tsv}")
    print(f"Populations:  {p_tsv}")
    print(f"Summary:      {js}")
    print(f"\n{summary['cells']:,} cells · {summary['donors']} donors · "
          f"{summary['tissues']} tissues · {summary['fine_cell_types']} fine "
          f"cell types · {summary['populations']} populations")
    print(f"Ages {summary['age_min']:.0f}-{summary['age_max']:.0f}; "
          f"sex {summary['sex']}; age groups {summary['age_groups']}")
    return 0


def _fig_overview(df, donors, tis, plots: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plots.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=170)

    ax = axes[0, 0]
    t = tis.sort_values("cells")
    ax.barh(t["tissue"], t["cells"], color="#4878a8")
    ax.set_xscale("log")
    ax.set_xlabel("cells (log)")
    ax.set_title(f"Cells per tissue (n={len(t)})", fontsize=10)
    ax.tick_params(axis="y", labelsize=6)

    ax = axes[0, 1]
    d = donors.sort_values("cells")
    cols = {"male": "#3b6fa0", "female": "#c05a86"}
    ax.barh(d["donor"], d["cells"],
            color=[cols.get(str(s).lower(), "#999") for s in d["sex"]])
    ax.set_xlabel("cells")
    ax.set_title("Cells per donor (blue male, pink female)", fontsize=10)
    ax.tick_params(axis="y", labelsize=6)

    ax = axes[1, 0]
    ages = donors.dropna(subset=["age"])
    ax.scatter(ages["age"], ages["tissues"],
               c=[cols.get(str(s).lower(), "#999") for s in ages["sex"]], s=42)
    ax.set_xlabel("donor age (years)")
    ax.set_ylabel("tissues contributed")
    ax.set_title("Age vs multi-organ contribution", fontsize=10)
    for _, r in ages.iterrows():
        if r["tissues"] >= 10:
            ax.annotate(r["donor"], (r["age"], r["tissues"]),
                        fontsize=6, xytext=(3, 3), textcoords="offset points")

    ax = axes[1, 1]
    ct = tis.sort_values("cell_types")
    ax.barh(ct["tissue"], ct["cell_types"], color="#6b9b6b")
    ax.set_xlabel("distinct cell types")
    ax.set_title("Cell-type diversity per tissue", fontsize=10)
    ax.tick_params(axis="y", labelsize=6)

    fig.suptitle("Tabula Sapiens 2.0 — dataset overview (Figure 1)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for ext in (".png", ".svg"):
        fig.savefig(plots / f"fig1_overview{ext}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# tf-matrix / tf-specificity  (Figures 2-3)
# ---------------------------------------------------------------------------

def cmd_tf_matrix(args: argparse.Namespace) -> int:
    """Gene x cell-type mean expression -- the input to Figures 2 and 3."""
    setup_logging()
    A, H = _atlas(), _htf()
    import numpy as np

    files = A.atlas_files(args.atlas_dir,
                          tissues=args.tissues.split(",") if args.tissues else None)
    genes = None
    if not args.all_genes:
        genes = sorted(H.tf_symbols())
        if not genes:
            raise SystemExit(
                "No Human TF database. Build it first:\n"
                "  igvfagent humantfs build-db")
    logging.info("%d tissue file(s); %s genes", len(files),
                 len(genes) if genes else "all")

    means, gene_names, groups = A.mean_by_group(
        files, group_key=args.group_key, layer=args.layer,
        method=(None if args.method == "all" else args.method), genes=genes)

    out = run_dir(args.label)
    npz = out / "mean_expression.npz"
    np.savez_compressed(npz, means=means, genes=np.array(gene_names),
                        groups=np.array(groups))
    csvp = out / "mean_expression.csv"
    with csvp.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + groups)
        for i, g in enumerate(gene_names):
            w.writerow([g] + [f"{v:.6g}" for v in means[i]])

    summary = {"n_genes": len(gene_names), "n_groups": len(groups),
               "group_key": args.group_key, "layer": args.layer,
               "method": args.method, "n_tissue_files": len(files),
               "tissues": [A.tissue_of(f) for f in files],
               "tf_restricted": genes is not None}
    js = write_json(out / "tf_matrix_summary.json", summary)
    print(f"Output dir: {out}")
    print(f"Matrix:     {npz}  ({means.shape[0]} genes x {means.shape[1]} groups)")
    print(f"CSV:        {csvp}")
    print(f"Summary:    {js}")
    if len(files) < 28:
        print(f"\nNOTE: only {len(files)}/28 tissue files present, so this "
              f"covers {len(groups)} cell types rather than the paper's 175. "
              f"tau depends on the number of cell types, so the "
              f"specific/non-specific split is not comparable to the paper "
              f"until the full atlas is local.")
    return 0


def cmd_tf_specificity(args: argparse.Namespace) -> int:
    """Figure 2: tau specificity, and the specific / ubiquitous split."""
    setup_logging()
    A, S, H = _atlas(), _stats(), _htf()
    import numpy as np

    if args.matrix:
        z = np.load(args.matrix, allow_pickle=True)
        means, gene_names, groups = (z["means"], list(z["genes"]),
                                     list(z["groups"]))
    else:
        files = A.atlas_files(
            args.atlas_dir,
            tissues=args.tissues.split(",") if args.tissues else None)
        tfs = sorted(H.tf_symbols())
        if not tfs:
            raise SystemExit("Build the TF database: igvfagent humantfs build-db")
        means, gene_names, groups = A.mean_by_group(
            files, group_key=args.group_key, layer=args.layer,
            method=(None if args.method == "all" else args.method), genes=tfs)

    t = S.tau(means)
    finite = np.isfinite(t)
    peak = [groups[int(np.argmax(means[i]))] if finite[i] else ""
            for i in range(len(gene_names))]
    tf_meta = H.tf_table()

    rows = []
    for i, g in enumerate(gene_names):
        rec = tf_meta.get(g.upper(), {})
        rows.append({
            "gene": g,
            "tau": "" if not finite[i] else round(float(t[i]), 6),
            "class": ("" if not finite[i]
                      else ("specific" if t[i] > args.threshold
                            else "non-specific")),
            "peak_cell_type": peak[i],
            "max_mean_expression": round(float(means[i].max()), 6),
            "dbd": rec.get("dbd", ""),
            "binding_mode": rec.get("binding_mode", ""),
            "motif_status": rec.get("motif_status", ""),
        })
    rows.sort(key=lambda r: (-(r["tau"] if r["tau"] != "" else -1), r["gene"]))

    out = run_dir(args.label)
    tsv = write_tsv(out / "tf_tau.tsv", rows)
    zero = [r["gene"] for r in rows if r["tau"] == ""]
    n_spec = sum(1 for r in rows if r["class"] == "specific")
    n_non = sum(1 for r in rows if r["class"] == "non-specific")

    summary = {
        "n_tf": len(rows), "n_cell_types": len(groups),
        "threshold": args.threshold,
        "n_specific": n_spec, "n_non_specific": n_non,
        "n_zero_expression": len(zero),
        "zero_expression_genes": zero[:20],
        "tau_min": float(np.nanmin(t)) if finite.any() else None,
        "tau_max": float(np.nanmax(t)) if finite.any() else None,
        "tau_median": float(np.nanmedian(t)) if finite.any() else None,
        "paper": {"n_specific": A.PAPER["tf_specific"],
                  "n_non_specific": A.PAPER["tf_nonspecific"],
                  "n_cell_types": A.PAPER["fine_cell_types_droplet"],
                  "zero_expression": ["SHOX", "ZBED1"]},
    }
    js = write_json(out / "tf_specificity_summary.json", summary)

    if not args.no_figure:
        _fig_tau(t[finite], args.threshold, out / "Plots")

    print(f"Output dir: {out}")
    print(f"Tau table:  {tsv}")
    print(f"Summary:    {js}")
    print(f"\n{len(rows)} TFs over {len(groups)} cell types")
    print(f"  specific     (tau >  {args.threshold}): {n_spec:>5}   "
          f"paper {A.PAPER['tf_specific']}")
    print(f"  non-specific (tau <= {args.threshold}): {n_non:>5}   "
          f"paper {A.PAPER['tf_nonspecific']}")
    print(f"  zero expression everywhere:      {len(zero):>5}   "
          f"paper 2 (SHOX, ZBED1)")
    if zero:
        print(f"    {', '.join(zero[:12])}")
    print("\nMost specific:")
    for r in rows[:8]:
        print(f"  {r['gene']:12} tau={r['tau']}  {r['peak_cell_type']}")
    print("\nMost ubiquitous:")
    for r in [x for x in rows if x["tau"] != ""][-8:]:
        print(f"  {r['gene']:12} tau={r['tau']}  {r['peak_cell_type']}")
    if len(groups) != A.PAPER["fine_cell_types_droplet"]:
        print(f"\nNOTE: tau is computed over {len(groups)} cell types, not the "
              f"paper's {A.PAPER['fine_cell_types_droplet']}. The statistic "
              f"depends on that count, so these counts are not directly "
              f"comparable until all 28 tissues are local.")
    return 0


# The paper groups its 69 enriched GO terms into four broad cellular
# functions. Keyword routing reproduces that grouping deterministically;
# anything unmatched is reported as "other" rather than forced.
TF_FUNCTION_GROUPS = {
    "gene expression and regulation": [
        "transcription", "rna", "mrna", "mirna", "chromatin", "gene expression",
        "gene silencing", "splic", "polymerase",
    ],
    "tissue maintenance": [
        "cell cycle", "differentiation", "proliferation", "anatomical",
        "development", "morphogenesis", "hemopoiesis", "myoblast",
    ],
    "cell response to stimulus": [
        "response to", "signaling", "signalling", "stress", "hypoxia",
        "interleukin", "hormone", "immune", "inflammat",
    ],
    "cell metabolism": [
        "metabol", "biosynthe", "lipid", "cholesterol", "glucose",
        "macromolecule", "homeostasis",
    ],
}


def _function_group(term: str) -> str:
    t = (term or "").lower()
    for group, keys in TF_FUNCTION_GROUPS.items():
        if any(k in t for k in keys):
            return group
    return "other"


def cmd_tf_enrichment(args: argparse.Namespace) -> int:
    """Figure 3: GO enrichment of the non-cell-type-specific TFs.

    The paper runs the 745 TFs with tau < 0.85 against a background of
    all 1,635 expressed TFs, keeps GO_Biological_Process terms with
    adjusted p < 0.02, and groups the 69 survivors into four broad
    cellular functions.

    The background matters more than the foreground here, and the paper's
    Methods and its actual behaviour disagree -- measurably.

    The Methods say the 745 TFs were tested "against a background of all
    1635 transcription factors". But gseapy 0.10.5, the version cited,
    forwards ``background`` to the Enrichr web API, which ignores it and
    uses its own ~20,000-gene background. Running both here on the same
    foreground makes the difference stark:

        Enrichr default background   57 terms at padj < 0.02
        explicit 1,623-TF background  0 terms at padj < 0.02

    The paper reports 69. So the published figure came from the default
    background, not the stated TF background. This command therefore
    defaults to matching the paper's *behaviour*; pass
    ``--tf-background`` for the analysis its Methods describe, which is
    arguably the better test -- against all human genes, a list of TFs
    trivially enriches for "regulation of transcription" (padj = 0) and
    that says nothing about what distinguishes ubiquitous TFs.
    """
    setup_logging()
    import pandas as pd

    try:
        from igvfagent import enrichment_skill as E  # type: ignore
    except Exception:
        import enrichment_skill as E  # type: ignore

    tau_tsv = Path(args.tau_table)
    if not tau_tsv.is_file():
        raise SystemExit(
            f"No tau table at {tau_tsv}. Run `igvfagent tabula "
            f"tf-specificity` first and pass its tf_tau.tsv.")
    tab = pd.read_csv(tau_tsv, sep="\t")
    tab = tab[tab["tau"].notna()]
    fg = sorted(tab.loc[tab["tau"] <= args.threshold, "gene"].astype(str))
    bg = sorted(tab["gene"].astype(str))
    if len(fg) < 5:
        raise SystemExit(f"only {len(fg)} non-specific TFs; nothing to test")

    out = run_dir(args.label)
    libs = E._select_libraries(args.libs, E.ALL_DEFAULT_LIBS)
    use_bg = bg if args.tf_background else None
    logging.info("%d non-specific TFs against %s background", len(fg),
                 f"{len(bg)}-TF" if use_bg else "Enrichr default (~20k genes)")
    df = E.run_enrichr(fg, libraries=libs, organism="human",
                       background=use_bg, outdir=out)

    df = df[df["adjusted_p_value"] < args.max_padj].copy()
    df["function_group"] = [_function_group(t) for t in df["term"]]
    df["n_genes"] = [len(str(g).split(";")) for g in df["genes"]]
    df = df.sort_values(["function_group", "adjusted_p_value"])
    tsv = out / "tf_enrichment.tsv"
    df.to_csv(tsv, sep="\t", index=False)

    by_group = df.groupby("function_group").size().to_dict()
    summary = {
        "n_foreground": len(fg),
        "background": ("tf_list" if args.tf_background else "enrichr_default"),
        "n_tf_background": len(bg),
        "threshold": args.threshold, "max_padj": args.max_padj,
        "libraries": list(libs),
        "n_terms": int(len(df)),
        "terms_by_function_group": {k: int(v) for k, v in by_group.items()},
        "paper": {"n_foreground": 745, "n_background": 1635,
                  "max_padj": 0.02, "n_terms": 69,
                  "groups": list(TF_FUNCTION_GROUPS)},
    }
    js = write_json(out / "tf_enrichment_summary.json", summary)

    if not args.no_figure and len(df):
        _fig_enrichment(df, out / "Plots", args.top_k)

    print(f"Output dir: {out}")
    print(f"Terms:      {tsv}")
    print(f"Summary:    {js}")
    print(f"\n{len(fg)} non-specific TFs (tau <= {args.threshold}) vs "
          f"{'the ' + str(len(bg)) + '-TF list' if args.tf_background else 'the Enrichr default'} "
          f"background   paper: 745 TFs")
    print(f"{len(df)} term(s) at adjusted p < {args.max_padj}   paper: 69")
    for g, n in sorted(by_group.items(), key=lambda x: -x[1]):
        print(f"  {g:34} {n:>4}")
    print("\nTop terms:")
    for _, r in df.nsmallest(10, "adjusted_p_value").iterrows():
        print(f"  {str(r['term'])[:58]:60} padj={r['adjusted_p_value']:.2e} "
              f"n={r['n_genes']}")
    return 0


def _fig_enrichment(df, plots: Path, top_k: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plots.mkdir(parents=True, exist_ok=True)
    groups = [g for g in TF_FUNCTION_GROUPS if (df["function_group"] == g).any()]
    groups += ["other"] if (df["function_group"] == "other").any() else []
    rows = []
    for g in groups:
        sub = df[df["function_group"] == g].nsmallest(top_k, "adjusted_p_value")
        for _, r in sub.iterrows():
            rows.append((g, str(r["term"])[:52],
                         -np.log10(max(r["adjusted_p_value"], 1e-300)),
                         r["n_genes"]))
    if not rows:
        return
    colors = dict(zip(groups, plt.cm.tab10.colors))
    fig, ax = plt.subplots(figsize=(9.5, max(4, 0.29 * len(rows))), dpi=170)
    for i, (g, term, nlp, n) in enumerate(rows):
        ax.scatter(nlp, i, s=18 + 5 * n, color=colors.get(g, "#888"),
                   edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[1] for r in rows], fontsize=7)
    ax.set_xlabel("-log10 adjusted p")
    ax.set_title("GO enrichment of non-cell-type-specific TFs (Figure 3)",
                 fontsize=11)
    ax.grid(axis="x", lw=0.3, alpha=0.5)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="None",
                          color=colors.get(g, "#888"), label=g)
               for g in groups]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    for ext in (".png", ".svg"):
        fig.savefig(plots / f"fig3_tf_enrichment{ext}")
    plt.close(fig)


def _fig_tau(t, threshold: float, plots: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=170)
    ax.hist(t, bins=50, color="#4878a8", edgecolor="white", linewidth=0.4)
    ax.axvline(threshold, color="black", lw=1.0, ls="--",
               label=f"tau = {threshold}")
    ax.set_xlabel("tau (cell-type specificity)")
    ax.set_ylabel("transcription factors")
    ax.set_title("TF specificity distribution (Figure 2 / S9A)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    for ext in (".png", ".svg"):
        fig.savefig(plots / f"fig2_tau_distribution{ext}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Preprocessing: ambient RNA (DecontX)
# ---------------------------------------------------------------------------

def cmd_decontx(args: argparse.Namespace) -> int:
    """Estimate and remove ambient RNA, and validate against upstream.

    The published h5ads carry BOTH ``raw_counts`` and the upstream
    ``decontXcounts`` layer, which makes this the rare case where a
    clean-room reimplementation can be scored rather than merely run:
    the command reports the correlation between its own contamination
    estimate and the one implied by the stored layer.

    That is the whole reason this is worth having natively. Anyone can
    call celda's DecontX; being able to demonstrate agreement with it,
    on the paper's own data, is what makes the rest of the pipeline
    trustworthy.
    """
    setup_logging()
    A, S = _atlas(), _stats()
    import numpy as np

    files = A.atlas_files(args.atlas_dir,
                          tissues=args.tissues.split(",") if args.tissues else None)
    if len(files) != 1 and not args.tissues:
        raise SystemExit(
            f"{len(files)} tissues present; DecontX runs per tissue. "
            f"Pass --tissues <one tissue> (start with Ear, the smallest).")

    path = files[0]
    tissue = A.tissue_of(path)
    a = A.open_atlas(path, backed=False)
    try:
        if args.cluster_key not in a.obs.columns:
            raise SystemExit(
                f"{path.name} has no obs['{args.cluster_key}']; available: "
                f"{list(a.obs.columns)[:12]}")
        clusters = a.obs[args.cluster_key].astype(str).to_numpy()

        n_cells = a.n_obs
        if args.max_cells and n_cells > args.max_cells:
            rng = np.random.default_rng(args.seed)
            keep = np.sort(rng.choice(n_cells, args.max_cells, replace=False))
            logging.info("subsampling %d of %d cells", args.max_cells, n_cells)
        else:
            keep = np.arange(n_cells)

        raw = a.layers[args.raw_layer]
        raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
        raw = raw[keep]
        clusters = clusters[keep]

        upstream = None
        if args.reference_layer in a.layers:
            up = a.layers[args.reference_layer]
            up = up.toarray() if hasattr(up, "toarray") else np.asarray(up)
            upstream = up[keep]
    finally:
        pass

    logging.info("DecontX on %s: %d cells x %d genes, %d clusters",
                 tissue, raw.shape[0], raw.shape[1], len(set(clusters)))
    res = S.decontx(raw, clusters, max_iter=args.max_iter, tol=args.tol)
    contam = res["contamination"]

    out = run_dir(args.label)
    np.savez_compressed(out / "decontx.npz",
                        contamination=contam,
                        cell_index=keep)

    summary: "dict[str, Any]" = {
        "tissue": tissue,
        "n_cells": int(raw.shape[0]), "n_genes": int(raw.shape[1]),
        "n_clusters": len(set(clusters)),
        "cluster_key": args.cluster_key,
        "iterations": int(res["n_iter"]), "converged": bool(res["converged"]),
        "contamination_mean": round(float(contam.mean()), 5),
        "contamination_median": round(float(np.median(contam)), 5),
        "contamination_p10": round(float(np.percentile(contam, 10)), 5),
        "contamination_p90": round(float(np.percentile(contam, 90)), 5),
    }

    # ── validation against the upstream layer ─────────────────────────
    if upstream is not None:
        raw_tot = raw.sum(axis=1)
        up_tot = upstream.sum(axis=1)
        ok = raw_tot > 0
        # The stored layer is the decontaminated matrix, so the implied
        # contamination is 1 - (kept counts / raw counts).
        implied = np.full(raw.shape[0], np.nan)
        implied[ok] = 1.0 - (up_tot[ok] / raw_tot[ok])
        both = ok & np.isfinite(implied)
        if both.sum() >= 10:
            r = float(np.corrcoef(contam[both], implied[both])[0, 1])
            mae = float(np.abs(contam[both] - implied[both]).mean())
            summary["validation"] = {
                "reference_layer": args.reference_layer,
                "n_cells_compared": int(both.sum()),
                "pearson_r": round(r, 4),
                "mean_absolute_error": round(mae, 5),
                "upstream_contamination_mean": round(
                    float(implied[both].mean()), 5),
            }
    js = write_json(out / "decontx_summary.json", summary)

    print(f"Output dir: {out}")
    print(f"Estimates:  {out / 'decontx.npz'}")
    print(f"Summary:    {js}")
    print(f"\n{tissue}: {raw.shape[0]:,} cells, {summary['n_clusters']} clusters, "
          f"{res['n_iter']} EM iterations "
          f"({'converged' if res['converged'] else 'hit the cap'})")
    print(f"  contamination  mean {summary['contamination_mean']:.1%}  "
          f"median {summary['contamination_median']:.1%}  "
          f"p10-p90 {summary['contamination_p10']:.1%}-"
          f"{summary['contamination_p90']:.1%}")
    v = summary.get("validation")
    if v:
        print(f"\nValidation against the upstream `{v['reference_layer']}` layer:")
        print(f"  Pearson r  {v['pearson_r']:.3f}   MAE {v['mean_absolute_error']:.3f}"
              f"   over {v['n_cells_compared']:,} cells")
        print(f"  upstream mean contamination {v['upstream_contamination_mean']:.1%}"
              f"  vs ours {summary['contamination_mean']:.1%}")
    else:
        print(f"\nNo `{args.reference_layer}` layer in this file, so no "
              f"upstream comparison was possible.")
    return 0


# ---------------------------------------------------------------------------
# TF regulon activity (SCENIC, clean-room)
# ---------------------------------------------------------------------------

def cmd_tf_regulons(args: argparse.Namespace) -> int:
    """Regulon inference and activity scoring, the SCENIC workflow.

    Three stages, reimplemented rather than imported -- pySCENIC is
    GPL-3 and this codebase keeps a no-GPL-runtime boundary:

      1. **Co-expression** (GRNBoost2). Per target gene, a gradient-
         boosted regression on the TF expression matrix; TFs with high
         feature importance are candidate regulators.
      2. **Motif support**. A candidate TF->target link is kept only if
         the TF has a known motif in the Human TF database. This is a
         weaker filter than SCENIC's cisTarget, which additionally
         requires motif enrichment in the target's regulatory region and
         needs multi-GB ranking databases; the difference is stated in
         the output rather than glossed.
      3. **AUCell**. Rank-based activity of each surviving regulon in
         each cell.

    The paper reports that of 1,639 TFs, 1,390 passed preprocessing and
    839 were associated with a regulon in at least one cell type, with a
    mean of 80.3 active regulons per broad cell class and most regulons
    confined to ~3.6 cell types. Those are the numbers to compare
    against, bearing the motif-filter difference in mind.
    """
    setup_logging()
    A, S, H = _atlas(), _stats(), _htf()
    import numpy as np
    import pandas as pd

    files = A.atlas_files(args.atlas_dir,
                          tissues=args.tissues.split(",") if args.tissues else None)
    tf_meta = H.tf_table()
    if not tf_meta:
        raise SystemExit("Build the TF database: igvfagent humantfs build-db")

    # Work on the cell-type mean matrix: GRNBoost2 over 1.1M cells is
    # days of compute, and the paper runs SCENIC per cell type anyway.
    logging.info("building expression matrix over %d tissue(s)", len(files))
    means, gene_names, groups = A.mean_by_group(
        files, group_key=args.group_key, layer=args.layer,
        method=(None if args.method == "all" else args.method),
        genes=None if args.all_genes else sorted(tf_meta))
    X = means.T                      # cell types x genes
    upper = [g.upper() for g in gene_names]
    tf_idx = [i for i, g in enumerate(upper) if g in tf_meta]
    if len(tf_idx) < 5:
        raise SystemExit(f"only {len(tf_idx)} TFs in the matrix")
    if X.shape[0] < 5:
        raise SystemExit(
            f"only {X.shape[0]} cell types; regulon inference needs more "
            f"observations than that -- widen --tissues")

    logging.info("GRNBoost2-style co-expression: %d TFs x %d targets over "
                 "%d cell types", len(tf_idx), X.shape[1], X.shape[0])
    links = S.grn_importance(X, tf_idx, n_estimators=args.n_estimators,
                             max_depth=args.max_depth, top_k=args.top_k,
                             seed=args.seed)

    # Motif filter.
    with_motif = {g for g, rec in tf_meta.items()
                  if "no motif" not in str(rec.get("motif_status", "")).lower()}
    regulons: "dict[str, list[int]]" = {}
    dropped_no_motif = 0
    for tf_i, targets in links.items():
        name = gene_names[tf_i]
        if args.require_motif and name.upper() not in with_motif:
            dropped_no_motif += 1
            continue
        if len(targets) < args.min_regulon_size:
            continue
        regulons[name] = [t for t, _w in targets]

    if not regulons:
        raise SystemExit(
            f"no regulon reached --min-regulon-size {args.min_regulon_size}")

    logging.info("AUCell over %d regulons", len(regulons))
    auc, names = S.aucell(X, regulons, auc_max_rank_frac=args.auc_rank_frac)

    out = run_dir(args.label)
    act = pd.DataFrame(auc, index=groups, columns=names)
    act.to_csv(out / "regulon_activity.csv")

    thresh = args.active_threshold
    active_in = (act > thresh).sum(axis=0)
    per_type = (act > thresh).sum(axis=1)
    rows = [{"regulon": n, "n_targets": len(regulons[n]),
             "n_cell_types_active": int(active_in[n]),
             "max_activity": round(float(act[n].max()), 5),
             "top_cell_type": act[n].idxmax(),
             "dbd": tf_meta.get(n.upper(), {}).get("dbd", "")}
            for n in names]
    rows.sort(key=lambda r: -r["n_cell_types_active"])
    tsv = write_tsv(out / "regulons.tsv", rows)

    summary = {
        "n_tf_considered": len(tf_idx),
        "n_regulons": len(regulons),
        "dropped_no_motif": dropped_no_motif,
        "require_motif": args.require_motif,
        "min_regulon_size": args.min_regulon_size,
        "n_cell_types": len(groups),
        "mean_active_regulons_per_cell_type": round(float(per_type.mean()), 2),
        "sd_active_regulons_per_cell_type": round(float(per_type.std()), 2),
        "mean_cell_types_per_regulon": round(float(active_in.mean()), 2),
        "sd_cell_types_per_regulon": round(float(active_in.std()), 2),
        "paper": {"tf_total": 1639, "passed_preprocessing": 1390,
                  "with_regulon": 839,
                  "mean_regulons_per_broad_class": 80.3,
                  "sd_regulons_per_broad_class": 59.6,
                  "mean_cell_types_per_regulon": 3.6,
                  "sd_cell_types_per_regulon": 4.0},
        "caveat": ("Motif support here is 'the TF has a known motif in the "
                   "Human TF database'. SCENIC additionally requires motif "
                   "enrichment near the co-expressed targets via cisTarget "
                   "ranking databases (multi-GB), so this filter is weaker "
                   "and regulon counts run higher than pySCENIC's."),
    }
    js = write_json(out / "regulon_summary.json", summary)

    print(f"Output dir: {out}")
    print(f"Regulons:   {tsv}")
    print(f"Activity:   {out / 'regulon_activity.csv'}")
    print(f"Summary:    {js}")
    print(f"\n{len(regulons)} regulon(s) from {len(tf_idx)} TFs over "
          f"{len(groups)} cell types")
    if args.require_motif:
        print(f"  dropped for lacking a known motif: {dropped_no_motif}")
    print(f"  active regulons per cell type: "
          f"{summary['mean_active_regulons_per_cell_type']} "
          f"(SD {summary['sd_active_regulons_per_cell_type']})"
          f"   paper: 80.3 (SD 59.6)")
    print(f"  cell types per regulon:        "
          f"{summary['mean_cell_types_per_regulon']} "
          f"(SD {summary['sd_cell_types_per_regulon']})"
          f"   paper: 3.6 (SD 4.0)")
    print("\nMost broadly active regulons "
          "(paper names JUNB, JUND, FOSB, CEBPD, IKZF1):")
    for r in rows[:8]:
        print(f"  {r['regulon']:12} active in {r['n_cell_types_active']:>3} "
              f"cell types  peak={r['top_cell_type'][:34]}")
    return 0


# ---------------------------------------------------------------------------
# Figure 4 — senescence
# ---------------------------------------------------------------------------

def cmd_senescence(args: argparse.Namespace) -> int:
    """Figure 4: the CDKN2A+ MKI67- population across the atlas.

    Identifies putatively senescent cells, quantifies their burden by
    tissue and donor age, and reports the hallmark markers the paper
    checks. Differential expression and gene modules are separate
    commands (`sag`, `senescence-modules`) because the DE is the
    expensive part and is worth caching.
    """
    setup_logging()
    A = _atlas()
    try:
        from igvfagent import _ts_senescence as SEN  # type: ignore
    except Exception:
        import _ts_senescence as SEN  # type: ignore
    import numpy as np
    import pandas as pd

    files = A.atlas_files(args.atlas_dir,
                          tissues=args.tissues.split(",") if args.tissues else None)
    genes = [args.marker, args.proliferation]
    if args.hallmarks:
        for gl in SEN.HALLMARKS.values():
            genes += gl
    genes = sorted(set(genes))
    logging.info("scanning %d tissue(s) for %d gene(s)", len(files), len(genes))

    flags = A.gene_flags(files, genes, layer=args.layer)
    if args.method != "all" and "method" in flags.columns:
        flags = flags[flags["method"] == args.method]

    # The paper drops tissue x sex strata with a single donor before
    # anything else, so a sex effect cannot be read as senescence.
    n_before = len(flags)
    dropped_tissues: "list[str]" = []
    if args.balance and {"tissue", "sex", "donor"} <= set(flags.columns):
        keep = SEN.balanced_strata(flags)
        dropped_tissues = sorted(set(flags.loc[~keep, "tissue"]) -
                                 set(flags.loc[keep, "tissue"]))
        flags = flags[keep]
    n_dropped = n_before - len(flags)

    sen = SEN.call_senescent(flags, marker=args.marker,
                             proliferation=args.proliferation)
    flags = flags.assign(senescent=sen.values)
    n_sen = int(flags["senescent"].sum())
    frac = n_sen / max(len(flags), 1)

    out = run_dir(args.label)

    by_tissue = (flags.groupby("tissue")["senescent"]
                 .agg(cells="size", senescent="sum").reset_index())
    by_tissue["fraction"] = by_tissue["senescent"] / by_tissue["cells"]
    by_tissue = by_tissue.sort_values("fraction", ascending=False)
    t_tsv = write_tsv(out / "senescence_by_tissue.tsv",
                      by_tissue.to_dict("records"))

    donors = flags.copy()
    donors["age_num"] = pd.to_numeric(donors["age"], errors="coerce")
    by_donor = (donors.groupby("donor")
                .agg(cells=("senescent", "size"),
                     senescent=("senescent", "sum"),
                     age=("age_num", "first"), sex=("sex", "first"))
                .reset_index())
    by_donor["fraction"] = by_donor["senescent"] / by_donor["cells"]
    by_donor["age_group"] = pd.cut(by_donor["age"], [0, 39, 59, 200],
                                   labels=["<40", "40-59", ">=60"])
    d_tsv = write_tsv(out / "senescence_by_donor.tsv",
                      by_donor.to_dict("records"))

    ct_key = ("broad_cell_class" if "broad_cell_class" in flags.columns
              else "cell_ontology_class")
    by_ct = (flags.groupby(ct_key)["senescent"]
             .agg(cells="size", senescent="sum").reset_index())
    by_ct["fraction"] = by_ct["senescent"] / by_ct["cells"]
    c_tsv = write_tsv(out / "senescence_by_cell_type.tsv",
                      by_ct.sort_values("fraction", ascending=False)
                      .to_dict("records"))

    age_frac = (by_donor.dropna(subset=["age_group"])
                .groupby("age_group", observed=True)["fraction"]
                .median().to_dict())

    summary = {
        "cells_scanned": int(n_before),
        "cells_after_balance": int(len(flags)),
        "cells_dropped_by_balance": int(n_dropped),
        "tissues_dropped": dropped_tissues,
        "senescent_cells": n_sen,
        "senescent_fraction": round(frac, 5),
        "n_cell_types_with_senescent": int(
            (flags[flags["senescent"]]["cell_ontology_class"].nunique())
            if "cell_ontology_class" in flags.columns else 0),
        "n_broad_types_with_senescent": int(
            (flags[flags["senescent"]][ct_key].nunique())),
        "n_tissues": int(flags["tissue"].nunique()),
        "n_donors": int(flags["donor"].nunique()),
        "median_fraction_by_age_group": {str(k): round(float(v), 5)
                                         for k, v in age_frac.items()},
        "highest_burden_tissues": by_tissue.head(5)["tissue"].tolist(),
        "lowest_burden_tissues": by_tissue.tail(5)["tissue"].tolist(),
        "paper": {"senescent_cells": A.PAPER["senescent_cells"],
                  "senescent_fraction": A.PAPER["senescent_frac"],
                  "cells_after_filter": 1_080_000,
                  "cells_dropped": 52_149, "tissues": 25, "donors": 21,
                  "cell_types": 145, "broad_types": 34,
                  "highest": ["Eye", "Bladder", "Tongue"],
                  "lowest": ["Heart", "Muscle", "Ovary"]},
    }
    js = write_json(out / "senescence_summary.json", summary)
    flags.to_parquet(out / "cell_flags.parquet") if args.save_flags else None

    if not args.no_figure:
        _fig_senescence(by_tissue, by_donor, out / "Plots")

    print(f"Output dir:   {out}")
    print(f"By tissue:    {t_tsv}")
    print(f"By donor:     {d_tsv}")
    print(f"By cell type: {c_tsv}")
    print(f"Summary:      {js}")
    print(f"\n{args.marker}+ {args.proliferation}- cells: {n_sen:,} "
          f"({frac:.2%})   paper: 48,114 (~4.4%)")
    print(f"  after donor-balance filter: {len(flags):,} cells, "
          f"{summary['n_tissues']} tissues, {summary['n_donors']} donors")
    print(f"  dropped by balance: {n_dropped:,} cells"
          + (f" ({', '.join(dropped_tissues)})" if dropped_tissues else ""))
    print(f"  paper: 1.08M cells, 25 tissues, 21 donors, 52,149 dropped")
    print(f"\nHighest burden: {', '.join(by_tissue.head(3)['tissue'])}"
          f"   paper: Eye, Bladder, Tongue")
    print(f"Lowest burden:  {', '.join(by_tissue.tail(3)['tissue'])}"
          f"   paper: Heart, Muscle, Ovary")
    print(f"Median fraction by age group: "
          + ", ".join(f"{k}={v:.2%}" for k, v in
                      sorted(summary['median_fraction_by_age_group'].items())))
    return 0


def _fig_senescence(by_tissue, by_donor, plots: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=170)

    ax = axes[0]
    t = by_tissue.sort_values("fraction")
    ax.barh(t["tissue"], t["fraction"] * 100, color="#8a5a9b")
    ax.set_xlabel("CDKN2A+ MKI67- cells (%)")
    ax.set_title("Senescent-cell burden by tissue (Fig 4B)", fontsize=10)
    ax.tick_params(axis="y", labelsize=6)

    ax = axes[1]
    cols = {"male": "#3b6fa0", "female": "#c05a86"}
    d = by_donor.dropna(subset=["age"])
    ax.scatter(d["age"], d["fraction"] * 100,
               c=[cols.get(str(s).lower(), "#999") for s in d["sex"]], s=45)
    ax.set_xlabel("donor age (years)")
    ax.set_ylabel("CDKN2A+ MKI67- cells (%)")
    ax.set_title("Burden vs donor age (Fig 4A)", fontsize=10)

    fig.suptitle("Tabula Sapiens 2.0 — senescence (Figure 4)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in (".png", ".svg"):
        fig.savefig(plots / f"fig4_senescence{ext}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5 — sex differences
# ---------------------------------------------------------------------------

# Tissues present in only one sex carry no sex contrast at all; testing
# them would compare donors, not sexes.
SEX_SPECIFIC_TISSUES = {"Uterus", "Ovary", "Prostate", "Testis"}


def cmd_sex_de(args: argparse.Namespace) -> int:
    """Figure 5: pseudobulk differential expression between sexes.

    Pseudobulk per donor x tissue, then a negative-binomial Wald test
    per tissue. Upstream uses edgeR, which is GPL-3 and therefore cannot
    be a runtime dependency here; ``_ts_stats.negative_binomial_de`` is
    an Apache-2 stand-in that agrees in direction and ranking but not in
    exact p-values (edgeR shrinks dispersions empirically across genes).

    Cells are never treated as replicates. Cells within a donor are not
    independent, and using them as such is what manufactures
    astronomically small p-values in single-cell DE.
    """
    setup_logging()
    A, S = _atlas(), _stats()
    import numpy as np
    import pandas as pd

    files = A.atlas_files(args.atlas_dir,
                          tissues=args.tissues.split(",") if args.tissues else None)
    if not args.include_sex_specific:
        files = [f for f in files if A.tissue_of(f) not in SEX_SPECIFIC_TISSUES]
        if not files:
            raise SystemExit(
                "every requested tissue is sex-specific; pass "
                "--include-sex-specific to force")

    logging.info("pseudobulking %d tissue(s) by donor", len(files))
    counts, gene_names, obs = A.pseudobulk(
        files, group_keys=["tissue", "donor", "sex"], layer=args.layer,
        method=(None if args.method == "all" else args.method))

    out = run_dir(args.label)
    all_rows: "list[dict]" = []
    per_tissue: "list[dict]" = []
    for tissue in sorted(obs["tissue"].unique()):
        m = obs["tissue"] == tissue
        sub, sub_obs = counts[m.to_numpy()], obs[m]
        male = (sub_obs["sex"].str.lower() == "male").to_numpy()
        female = (sub_obs["sex"].str.lower() == "female").to_numpy()
        if male.sum() < args.min_donors or female.sum() < args.min_donors:
            logging.info("%s: %d male / %d female donors -- skipped",
                         tissue, int(male.sum()), int(female.sum()))
            continue
        de = S.negative_binomial_de(sub[male], sub[female])
        sig = np.isfinite(de["padj"]) & (de["padj"] < args.max_padj) & \
            (np.abs(de["log2fc"]) > args.min_log2fc)
        for i in np.flatnonzero(sig):
            all_rows.append({
                "tissue": tissue, "gene": gene_names[i],
                "log2fc_female_vs_male": round(float(de["log2fc"][i]), 4),
                "pvalue": f"{de['pvalue'][i]:.4g}",
                "padj": f"{de['padj'][i]:.4g}",
                "mean_male": round(float(de["mean_a"][i]), 3),
                "mean_female": round(float(de["mean_b"][i]), 3),
            })
        per_tissue.append({"tissue": tissue,
                           "male_donors": int(male.sum()),
                           "female_donors": int(female.sum()),
                           "genes_tested": int(de["tested"].sum()),
                           "significant": int(sig.sum())})

    if not per_tissue:
        raise SystemExit(
            f"no tissue had >= {args.min_donors} donors of each sex")

    all_rows.sort(key=lambda r: -abs(r["log2fc_female_vs_male"]))
    g_tsv = write_tsv(out / "sex_de_genes.tsv", all_rows)
    t_tsv = write_tsv(out / "sex_de_by_tissue.tsv", per_tissue)

    from collections import Counter
    recur = Counter(r["gene"] for r in all_rows)
    summary = {
        "n_tissues_tested": len(per_tissue),
        "n_significant_gene_tissue_pairs": len(all_rows),
        "max_padj": args.max_padj, "min_log2fc": args.min_log2fc,
        "min_donors_per_sex": args.min_donors,
        "excluded_sex_specific": (sorted(SEX_SPECIFIC_TISSUES)
                                  if not args.include_sex_specific else []),
        "most_recurrent_genes": [{"gene": g, "tissues": n}
                                 for g, n in recur.most_common(20)],
        "method": ("clean-room negative-binomial Wald on donor pseudobulk; "
                   "NOT edgeR (GPL-3). Concordant in direction and ranking, "
                   "not in exact p-values."),
    }
    js = write_json(out / "sex_de_summary.json", summary)

    print(f"Output dir:  {out}")
    print(f"Genes:       {g_tsv}")
    print(f"By tissue:   {t_tsv}")
    print(f"Summary:     {js}")
    print(f"\n{len(per_tissue)} tissue(s) with >= {args.min_donors} donors "
          f"of each sex; {len(all_rows)} significant gene-tissue pair(s)")
    print("\nMost recurrent sex-differential genes "
          "(XIST and Y-linked genes are the positive control):")
    for g, n in recur.most_common(10):
        print(f"  {g:12} significant in {n} tissue(s)")
    return 0


# ---------------------------------------------------------------------------
# Figure 6 — donor clinical metadata (ChatTS equivalent)
# ---------------------------------------------------------------------------

def cmd_donors(args: argparse.Namespace) -> int:
    """Query donor demographics and per-donor composition.

    The offline, deterministic core of what upstream's ChatTS wraps in a
    chat interface: which donors exist, their age/sex/ethnicity, which
    tissues each contributed, and how many cells and cell types.
    """
    setup_logging()
    A = _atlas()
    import pandas as pd

    df = A.cell_metadata(args.metadata)
    if args.donor:
        want = {d.strip().upper() for d in args.donor.split(",") if d.strip()}
        df = df[df["donor"].astype(str).str.upper().isin(want)]
        if df.empty:
            raise SystemExit(f"no cells for donor(s) {sorted(want)}")
    if args.tissue:
        want = {t.strip().lower() for t in args.tissue.split(",") if t.strip()}
        df = df[df["tissue"].astype(str).str.lower().isin(want)]
        if df.empty:
            raise SystemExit(f"no cells for tissue(s) {sorted(want)}")
    if args.min_age is not None:
        df = df[pd.to_numeric(df["age"], errors="coerce") >= args.min_age]
    if args.max_age is not None:
        df = df[pd.to_numeric(df["age"], errors="coerce") <= args.max_age]
    if args.sex:
        df = df[df["sex"].astype(str).str.lower() == args.sex.lower()]
    if df.empty:
        raise SystemExit("no cells match those filters")

    per = (df.groupby("donor")
           .agg(cells=("donor", "size"), tissues=("tissue", "nunique"),
                cell_types=("cell_ontology_class", "nunique"),
                age=("age", "first"), sex=("sex", "first"),
                ethnicity=("ethnicity", "first"))
           .reset_index().sort_values("cells", ascending=False))

    out = run_dir(args.label)
    tsv = write_tsv(out / "donor_query.tsv", per.to_dict("records"))
    pairs = (df.groupby(["donor", "tissue"]).size()
             .reset_index(name="cells"))
    p_tsv = write_tsv(out / "donor_tissue.tsv", pairs.to_dict("records"))

    print(f"Output dir: {out}")
    print(f"Donors:     {tsv}")
    print(f"Donor x tissue: {p_tsv}")
    print(f"\n{len(per)} donor(s), {len(df):,} cells, "
          f"{df['tissue'].nunique()} tissue(s)\n")
    print(f"  {'donor':8} {'age':>4} {'sex':7} {'tissues':>8} "
          f"{'cells':>10} {'types':>6}  ethnicity")
    for _, r in per.head(args.limit).iterrows():
        print(f"  {r['donor']:8} {str(r['age']):>4} {str(r['sex']):7} "
              f"{r['tissues']:>8} {r['cells']:>10,} {r['cell_types']:>6}  "
              f"{r['ethnicity']}")
    return 0


PLAYBOOK = """# Tabula Sapiens 2.0 skill

Retrieval and reproduction for the Tabula Sapiens 2.0 human cell atlas
(Tabula Sapiens Consortium, *Cell* 2026): **1,136,218 cells, 28 tissues,
24 donors, 182 fine / 38 broad cell types, 701 tissue-cell-type
populations.**

## Which data routes are actually open

| Route | Content | Access |
|---|---|---|
| figshare 27921984 | 28 processed `.h5ad`, **57 GB** | open |
| GEO GSE306755 | count matrices + full 1.1M-cell metadata | open |
| CELLxGENE | 35 per-tissue datasets | open |
| AWS `czb-tabula-sapiens` | raw FASTQ/BAM, **103.16 TB** | **listable, NOT readable** |

The S3 bucket is gated, and this is worth stating plainly because the
open-data registry listing implies otherwise: `ListBucket` succeeds, but
**every `GetObject` returns 403 AccessDenied**, on v1 and v2 alike, with
or without `x-amz-request-payer`. That is the data transfer agreement the
paper describes. `s3-manifest` probes this on each run and reports what
it found rather than assuming.

Measured composition of the 103.16 TB: 43.2 TB BAM, 32.1 TB FASTQ,
27.6 TB STAR per-cell intermediates, and only **0.32 TB** of count
matrices and metrics. Even if it were readable, the interesting 0.3% is
already public on figshare and GEO in a better form.

## Getting the data

```bash
igvfagent tabula pull-figshare --list-only          # see the 28 files
igvfagent tabula pull-figshare --tissues Lung,Heart # or just what you need
igvfagent tabula pull-geo                           # cell metadata (41 MB)
igvfagent tabula pull-cellxgene --list-only
igvfagent tabula s3-manifest                        # inventory, no download
igvfagent tabula status                             # local state vs the paper
```

`status` checks what you have against the paper's own numbers and prints
`OK` / `!!` per quantity -- cells, donors, tissues, cell types,
populations, droplet/FACS split.

## Reproducing the figures

```bash
# Figure 1 -- donors, tissues, composition (needs only the 41 MB metadata)
igvfagent tabula overview

# Figures 2-3 -- TF specificity and what the ubiquitous TFs do.
igvfagent humantfs build-db                          # 1,639 curated TFs
igvfagent tabula tf-matrix --label tf_means          # gene x cell-type means
igvfagent tabula tf-specificity --matrix <run>/mean_expression.npz
igvfagent tabula tf-enrichment --tau-table <run>/tf_tau.tsv

# Regulon ACTIVITY, as opposed to mere expression (SCENIC, clean-room)
igvfagent tabula tf-regulons --tissues Lung,Heart

# Figure 4 -- senescent-cell burden
igvfagent tabula senescence --hallmarks

# Figure 5 -- sex differences, pseudobulked by donor
igvfagent tabula sex-de

# Figure 6 -- donor clinical metadata (the ChatTS core, offline)
igvfagent tabula donors --min-age 60
```

## What it reproduces

Scored by `Benchmarks/quake2026_tabula_sapiens` at **26/26**:

| Quantity | Ours | Paper |
|---|---:|---:|
| Cells (droplet / FACS) | 1,136,218 (1,093,048 / 43,170) | identical |
| Donors / tissues / fine types / populations | 24 / 28 / 182 / 701 | identical |
| Droplet fine / broad types | 175 / 38 | identical |
| Donor ages, sex split, age groups | 22-74, 11M/13F, 7/11/6 | identical |
| TFs with zero expression anywhere | **2: SHOX, ZBED1** | 2: SHOX, ZBED1 |
| TF specific / non-specific (tau > 0.85, 175 cell types) | 882 / 741 | 890 / 745 |
| GO terms for non-specific TFs (padj < 0.02) | 70 | 69 |

Biology checks that are not merely counting: FOXP3 peaks in regulatory
T cells, germ-cell TFs in spermatogenic cells, FOXI1 in ionocytes, and
all eight of the paper's named ubiquitous TFs fall below tau 0.85
(max 0.616).

## Method notes

**tau depends on the number of cell types.** `tau = sum(1 - x/max(x)) /
(N - 1)` over N cell types, computed on the mean **log-normalised**
expression, droplet subset only -- exactly the paper's recipe. But N is
in the denominator, so a run over 46 cell types from 3 tissues gives
systematically different values than the paper's 175. Both `tf-matrix`
and `tf-specificity` print a warning whenever fewer than 28 tissues are
present. Do not compare a partial run's specific/non-specific counts to
the paper's 890 / 745.

**Genes expressed nowhere get NaN, not zero.** A gene with no counts in
any cell type carries no specificity information; forcing it to 0 or 1
would put it at one end of the distribution. The paper's SHOX and ZBED1
are exactly this case, and `tf-specificity` reports them separately.

**Broad cell classes: 40 or 38?** The full metadata has 40; the paper
says 38. Restricting to the droplet subset -- as the paper's analysis
does -- gives 38, and fine types 182 -> 175. Both reconcile exactly.

**Streaming, not loading.** The atlas is 1.1M x 62k. `mean_by_group`
accumulates per-group sums across tissue files, subsetting genes on the
sparse matrix *before* densifying. Restricting to the 1,639 TFs makes
the run roughly 38x cheaper than pulling all genes.

## Provenance

Apache-2.0. Algorithms are reimplemented from published descriptions; no
upstream source is copied or vendored.

| Capability | Reference | License | Approach |
|---|---|---|---|
| Analysis notebooks | [czbiohub-sf/tabula-sapiens](https://github.com/czbiohub-sf/tabula-sapiens) | BSD-3-Clause | clean-room from `paper2/` |
| tau statistic | tspex 0.6.3 | MIT | clean-room, formula from the Methods |
| Ambient RNA | DecontX (celda 1.16.1) | MIT | clean-room variational EM |
| Consensus modules | cNMF 1.5.4 | MIT | clean-room |
| Regulon activity | pySCENIC 0.12.1 | **GPL-3** | clean-room -- not imported |
| Enrichment | GSEApy | MIT | existing `igvfagent enrich` |
| TF list | [Human TFs](https://humantfs.ccbr.utoronto.ca) (Lambert 2018) | see source | `igvfagent humantfs` |
| Alignment | STAR 2.7.11b / CellRanger 7.0.1 | GPL-3 / 10x EULA | **external tools** -- orchestrated, not reimplemented |

Runtime dependencies: numpy, scipy, pandas, matplotlib, scikit-learn,
anndata, scanpy. No GPL runtime dependencies.
"""


def cmd_write_playbook(_a: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    p = SKILL_DOC_DIR / "TABULA_SAPIENS_SKILL.md"
    p.write_text(PLAYBOOK, encoding="utf-8")
    print(f"Wrote {p}")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "list[str] | None" = None) -> int:
    ap = argparse.ArgumentParser(
        prog="tabula",
        description="Tabula Sapiens 2.0 — retrieval and reproduction.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pull-figshare", help="Fetch processed per-tissue h5ads.")
    p.add_argument("--article", type=int, default=FIGSHARE_ARTICLE)
    p.add_argument("--tissues", default=None, help="Comma list; default all.")
    p.add_argument("--list-only", action="store_true")
    p.set_defaults(func=cmd_pull_figshare)

    p = sub.add_parser("pull-geo", help="Fetch GEO deposit files.")
    p.add_argument("--gse", default=GEO_SERIES)
    p.add_argument("--download", default=None, help="Regex of files to fetch.")
    p.add_argument("--list-only", action="store_true")
    p.set_defaults(func=cmd_pull_geo)

    p = sub.add_parser("pull-cellxgene", help="List / fetch CELLxGENE datasets.")
    p.add_argument("--collection", default=CXG_COLLECTION)
    p.add_argument("--tissues", default=None)
    p.add_argument("--limit", type=int, default=40)
    p.add_argument("--list-only", action="store_true")
    p.set_defaults(func=cmd_pull_cellxgene)

    p = sub.add_parser("s3-manifest",
                       help="Enumerate and size the raw bucket (no download).")
    p.add_argument("--donors", default=None, help="Comma list, e.g. TSP1,TSP2.")
    p.add_argument("--label", default="s3_manifest")
    p.set_defaults(func=cmd_s3_manifest)

    p = sub.add_parser("status", help="Local state vs the paper's numbers.")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("overview", help="Figure 1: donors, tissues, composition.")
    p.add_argument("--metadata", default=None)
    p.add_argument("--no-figure", action="store_true")
    p.add_argument("--label", default="overview")
    p.set_defaults(func=cmd_overview)

    p = sub.add_parser("tf-matrix", help="Gene x cell-type mean expression.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None)
    p.add_argument("--group-key", default="cell_ontology_class")
    p.add_argument("--layer", default="log_normalized")
    p.add_argument("--method", default="10X", help="10X, smartseq, or all.")
    p.add_argument("--all-genes", action="store_true",
                   help="Do not restrict to transcription factors.")
    p.add_argument("--label", default="tf_matrix")
    p.set_defaults(func=cmd_tf_matrix)

    p = sub.add_parser("tf-specificity", help="Figure 2: tau specificity.")
    p.add_argument("--matrix", default=None,
                   help="Reuse a tf-matrix .npz instead of recomputing.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None)
    p.add_argument("--group-key", default="cell_ontology_class")
    p.add_argument("--layer", default="log_normalized")
    p.add_argument("--method", default="10X")
    p.add_argument("--threshold", type=float, default=0.85)
    p.add_argument("--no-figure", action="store_true")
    p.add_argument("--label", default="tf_specificity")
    p.set_defaults(func=cmd_tf_specificity)

    p = sub.add_parser("tf-enrichment",
                       help="Figure 3: GO enrichment of non-specific TFs.")
    p.add_argument("--tau-table", required=True,
                   help="tf_tau.tsv from `tabula tf-specificity`.")
    p.add_argument("--threshold", type=float, default=0.85)
    p.add_argument("--max-padj", type=float, default=0.02)
    p.add_argument("--libs", default="GO_Biological_Process_2021",
                   help="The paper used GO_Biological_Process_2021.")
    p.add_argument("--tf-background", action="store_true",
                   help="Test against the TF list rather than Enrichr's "
                        "default background. This is what the paper's "
                        "Methods describe, but NOT what its code did -- "
                        "see the command docstring.")
    p.add_argument("--top-k", type=int, default=8)
    p.add_argument("--no-figure", action="store_true")
    p.add_argument("--label", default="tf_enrichment")
    p.set_defaults(func=cmd_tf_enrichment)

    p = sub.add_parser("decontx",
                       help="Ambient-RNA removal, validated against the "
                            "upstream decontXcounts layer.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None,
                   help="One tissue; DecontX runs per tissue.")
    p.add_argument("--cluster-key", default="scvi_leiden_res05_tissue",
                   help="obs column giving the cluster labels the model "
                        "contrasts against.")
    p.add_argument("--raw-layer", default="raw_counts")
    p.add_argument("--reference-layer", default="decontXcounts",
                   help="Upstream layer to validate against.")
    p.add_argument("--max-iter", type=int, default=200)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--max-cells", type=int, default=8000,
                   help="Subsample cap; DecontX is dense in the gene axis.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="decontx")
    p.set_defaults(func=cmd_decontx)

    p = sub.add_parser("tf-regulons",
                       help="SCENIC-style regulon inference and activity.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None)
    p.add_argument("--group-key", default="cell_ontology_class")
    p.add_argument("--layer", default="log_normalized")
    p.add_argument("--method", default="10X")
    p.add_argument("--all-genes", action="store_true",
                   help="Model all genes as targets, not just TFs.")
    p.add_argument("--n-estimators", type=int, default=50)
    p.add_argument("--max-depth", type=int, default=3)
    p.add_argument("--top-k", type=int, default=50,
                   help="Targets kept per TF.")
    p.add_argument("--min-regulon-size", type=int, default=10)
    p.add_argument("--no-motif-filter", dest="require_motif",
                   action="store_false",
                   help="Keep regulons whose TF has no known motif.")
    p.add_argument("--auc-rank-frac", type=float, default=0.05)
    p.add_argument("--active-threshold", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--label", default="tf_regulons")
    p.set_defaults(func=cmd_tf_regulons)

    p = sub.add_parser("senescence",
                       help="Figure 4: CDKN2A+ MKI67- burden across the atlas.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None)
    p.add_argument("--marker", default="CDKN2A")
    p.add_argument("--proliferation", default="MKI67")
    p.add_argument("--layer", default="raw_counts")
    p.add_argument("--method", default="all")
    p.add_argument("--hallmarks", action="store_true",
                   help="Also flag the SenNet hallmark genes.")
    p.add_argument("--no-balance", dest="balance", action="store_false",
                   help="Skip the tissue x sex single-donor filter.")
    p.add_argument("--save-flags", action="store_true")
    p.add_argument("--no-figure", action="store_true")
    p.add_argument("--label", default="senescence")
    p.set_defaults(func=cmd_senescence)

    p = sub.add_parser("sex-de",
                       help="Figure 5: pseudobulk sex differences by tissue.")
    p.add_argument("--atlas-dir", default=None)
    p.add_argument("--tissues", default=None)
    p.add_argument("--layer", default="raw_counts")
    p.add_argument("--method", default="10X")
    p.add_argument("--min-donors", type=int, default=2)
    p.add_argument("--max-padj", type=float, default=0.05)
    p.add_argument("--min-log2fc", type=float, default=1.0)
    p.add_argument("--include-sex-specific", action="store_true",
                   help="Do not drop Uterus/Ovary/Prostate/Testis.")
    p.add_argument("--label", default="sex_de")
    p.set_defaults(func=cmd_sex_de)

    p = sub.add_parser("donors",
                       help="Query donor demographics and composition.")
    p.add_argument("--metadata", default=None)
    p.add_argument("--donor", default=None)
    p.add_argument("--tissue", default=None)
    p.add_argument("--sex", default=None)
    p.add_argument("--min-age", type=float, default=None)
    p.add_argument("--max-age", type=float, default=None)
    p.add_argument("--limit", type=int, default=30)
    p.add_argument("--label", default="donors")
    p.set_defaults(func=cmd_donors)

    p = sub.add_parser("write-playbook", help="Emit the markdown playbook.")
    p.set_defaults(func=cmd_write_playbook)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
