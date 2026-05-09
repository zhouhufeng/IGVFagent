#!/usr/bin/env python3
"""Proteomics & PPI skill.

Pulls and version-tracks well-maintained protein-protein interaction
datasets (BioGRID, IntAct, HuRI), pathway databases (Reactome, KEGG),
and IGVF Portal protein assays / files. Integrates everything into a
local SQLite proteomics knowledge graph (PPI-KG) that is incrementally
update-aware. Also runs summary statistics, generates network
visualizations, surveys the Nature/Cell/Science literature for IGVF
protein assays (VAMP-seq family, MAVE, semi-qY2H, DUAL-IPA) via the
existing Reference skill, and produces representative example figures
per assay from real IGVF Portal files.

Subcommands
-----------
    download <source>           Fetch latest from one source
                                  (biogrid|intact|huri|reactome|kegg|igvf|uniprot|all)
    versions                     Local-vs-upstream version table
    update                       Pull only sources where upstream is newer
    igvf-protein                 IGVF Portal protein-assay manifest + PPI file pull
    build-kg                     Ingest downloaded sources into local SQLite KG
    kg-stats                     Summary statistics over the integrated KG
    kg-visualize <gene>          Degree distribution, top hubs, ego graph
    assay-survey                 Reference skill literature scan per IGVF protein assay
    assay-figures                Per-assay example plots from real IGVF Portal files
    pipeline                     End-to-end: download → kg → stats → viz → survey
    write-playbook               Write Docs/Skills/PROTEOMICS_SKILLS.md

Outputs land under
  Data/Proteomics/Sources/<src>/...
  Data/Proteomics/KG/proteomics.sqlite
  Docs/Proteomics/<ts>_<label>/...
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable, Optional

# ---------------------------------------------------------------------------
# Project paths and endpoints
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "Proteomics"
SOURCES_DIR = DATA_DIR / "Sources"
KG_DIR = DATA_DIR / "KG"
KG_PATH = KG_DIR / "proteomics.sqlite"
DOCS_DIR = ROOT / "Docs" / "Proteomics"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "PROTEOMICS_SKILLS.md"

BIOGRID_BASE = _resolve_endpoint("biogrid", "BIOGRID_BASE")
INTACT_BASE = _resolve_endpoint("intact_ftp", "INTACT_BASE")
HURI_BASE = _resolve_endpoint("huri", "HURI_BASE")
REACTOME_BASE = _resolve_endpoint("reactome_dl", "REACTOME_BASE")
KEGG_BASE = _resolve_endpoint("kegg_rest", "KEGG_BASE")
UNIPROT_IDMAP_BASE = _resolve_endpoint("uniprot_idmap", "UNIPROT_IDMAP_BASE")
PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")

logger = logging.getLogger("proteomics")

# ---------------------------------------------------------------------------
# Generic utilities
# ---------------------------------------------------------------------------


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"proteomics_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in (s or "run"))


def mkdirs() -> None:
    for d in (DATA_DIR, SOURCES_DIR, KG_DIR, DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def http_request(url: str, *, method: str = "GET", timeout: int = 60,
                 headers: Optional[dict] = None) -> tuple[int, bytes, dict]:
    h = {"User-Agent": "IGVFagent-Proteomics/0.1",
         "Accept": "application/json,text/plain,*/*"}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, headers=h, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read(), dict(resp.headers)
    except urllib.error.HTTPError as e:
        body = b""
        try:
            body = e.read()
        except Exception:
            pass
        return e.code, body, dict(getattr(e, "headers", {}) or {})
    except Exception as e:
        logger.warning("HTTP error %s: %s", url, e)
        return 0, b"", {}


def http_head(url: str, *, timeout: int = 30) -> dict:
    """Return response headers via HEAD; falls back to a 1-byte ranged GET."""
    s, _, h = http_request(url, method="HEAD", timeout=timeout)
    if s in (200, 204) and h:
        return h
    s, _, h = http_request(url, headers={"Range": "bytes=0-0"}, timeout=timeout)
    if s in (200, 206):
        return h
    return {}


def http_download(url: str, dest: Path, *, chunk: int = 1 << 20,
                  timeout: int = 600) -> int:
    """Stream URL → dest (atomic via .part)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    h = {"User-Agent": "IGVFagent-Proteomics/0.1"}
    req = urllib.request.Request(url, headers=h)
    n = 0
    with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as f:
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            f.write(buf)
            n += len(buf)
    tmp.replace(dest)
    logger.info("Downloaded %s (%.1f MB)", dest.name, n / (1 << 20))
    return n


def sha256_of(path: Path, *, max_bytes: Optional[int] = None) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        if max_bytes:
            h.update(f.read(max_bytes))
        else:
            for buf in iter(lambda: f.read(1 << 20), b""):
                h.update(buf)
    return h.hexdigest()


def open_text(path: Path):
    """Open path for text reading, transparently decompressing .gz."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Version tracking
# ---------------------------------------------------------------------------

VERSION_FILE = DATA_DIR / "_versions.json"


def load_versions() -> dict:
    if VERSION_FILE.exists():
        try:
            return json.loads(VERSION_FILE.read_text())
        except Exception:
            return {}
    return {}


def save_versions(v: dict) -> None:
    VERSION_FILE.parent.mkdir(parents=True, exist_ok=True)
    VERSION_FILE.write_text(json.dumps(v, indent=2, sort_keys=True))


def update_version(source: str, *, version: str, url: str,
                    path: Optional[Path] = None,
                    n_records: Optional[int] = None) -> None:
    v = load_versions()
    entry = {
        "version": version,
        "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "url": url,
    }
    if path and path.exists():
        entry["path"] = str(path)
        entry["sha256"] = sha256_of(path, max_bytes=64 << 20)
        entry["size_bytes"] = path.stat().st_size
    if n_records is not None:
        entry["n_records"] = n_records
    v[source] = entry
    save_versions(v)


# ---------------------------------------------------------------------------
# BioGRID
# ---------------------------------------------------------------------------


class _LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: "list[str]" = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for k, v in attrs:
                if k == "href":
                    self.links.append(v)


def biogrid_latest_release() -> Optional[str]:
    """Scrape Release-Archive listing for newest BIOGRID-X.Y.Z dir."""
    s, body, _ = http_request(f"{BIOGRID_BASE}/Release-Archive/")
    if s != 200:
        logger.warning("BioGRID release-archive HTTP %s", s)
        return None
    p = _LinkParser()
    p.feed(body.decode("utf-8", "replace"))
    versions = []
    for href in p.links:
        m = re.search(r"BIOGRID-(\d+\.\d+\.\d+)", href)
        if m:
            versions.append(m.group(1))
    if not versions:
        return None
    versions.sort(key=lambda x: tuple(int(p) for p in x.split(".")))
    return versions[-1]


def biogrid_download(target_version: Optional[str] = None) -> Optional[Path]:
    ver = target_version or biogrid_latest_release()
    if not ver:
        logger.error("BioGRID: could not determine latest release.")
        return None
    out = SOURCES_DIR / "BioGRID"
    out.mkdir(parents=True, exist_ok=True)
    fname = f"BIOGRID-ALL-{ver}.tab3.zip"
    url = f"{BIOGRID_BASE}/Release-Archive/BIOGRID-{ver}/{fname}"
    dest = out / fname
    if dest.exists():
        logger.info("BioGRID v%s already cached at %s", ver, dest)
    else:
        logger.info("Downloading BioGRID v%s ...", ver)
        try:
            http_download(url, dest, timeout=900)
        except Exception as e:
            logger.error("BioGRID download failed: %s", e)
            return None
    update_version("biogrid", version=ver, url=url, path=dest)
    return dest


def biogrid_iter_human(zip_path: Path, *, max_rows: int = 0) -> Iterable[dict]:
    """Stream BioGRID tab3 rows, filtered to human-human."""
    with zipfile.ZipFile(zip_path) as zf:
        members = [n for n in zf.namelist() if n.endswith(".tab3.txt")]
        if not members:
            return
        with zf.open(members[0]) as fp:
            text = io.TextIOWrapper(fp, encoding="utf-8", errors="replace")
            header = text.readline().rstrip("\r\n").lstrip("#").split("\t")
            idx = {c.strip(): i for i, c in enumerate(header)}
            n = 0
            for line in text:
                cells = line.rstrip("\r\n").split("\t")
                if len(cells) < len(header):
                    continue
                try:
                    if cells[idx["Organism Name Interactor A"]] != "Homo sapiens":
                        continue
                    if cells[idx["Organism Name Interactor B"]] != "Homo sapiens":
                        continue
                except KeyError:
                    org_a = cells[idx.get("Organism ID Interactor A", -1)] \
                        if idx.get("Organism ID Interactor A", -1) >= 0 else ""
                    org_b = cells[idx.get("Organism ID Interactor B", -1)] \
                        if idx.get("Organism ID Interactor B", -1) >= 0 else ""
                    if org_a != "9606" or org_b != "9606":
                        continue
                yield {
                    "id_a": cells[idx["Official Symbol Interactor A"]] or "",
                    "id_b": cells[idx["Official Symbol Interactor B"]] or "",
                    "id_type": "symbol",
                    "uniprot_a": (cells[idx.get("SWISS-PROT Accessions Interactor A", -1)]
                                  if idx.get("SWISS-PROT Accessions Interactor A", -1) >= 0 else "").split("|")[0],
                    "uniprot_b": (cells[idx.get("SWISS-PROT Accessions Interactor B", -1)]
                                  if idx.get("SWISS-PROT Accessions Interactor B", -1) >= 0 else "").split("|")[0],
                    "source_id": cells[idx.get("#BioGRID Interaction ID", -1)]
                                 if idx.get("#BioGRID Interaction ID", -1) >= 0
                                 else cells[idx.get("BioGRID Interaction ID", -1)],
                    "detection_method": cells[idx.get("Experimental System", -1)]
                                       if idx.get("Experimental System", -1) >= 0 else "",
                    "evidence_type": cells[idx.get("Experimental System Type", -1)]
                                    if idx.get("Experimental System Type", -1) >= 0 else "",
                    "pubmed_id": cells[idx.get("Publication Source", -1)]
                                if idx.get("Publication Source", -1) >= 0 else "",
                    "score": cells[idx.get("Score", -1)] if idx.get("Score", -1) >= 0 else "",
                    "throughput": cells[idx.get("Throughput", -1)]
                                 if idx.get("Throughput", -1) >= 0 else "",
                }
                n += 1
                if max_rows and n >= max_rows:
                    return


# ---------------------------------------------------------------------------
# IntAct
# ---------------------------------------------------------------------------


def intact_latest_version() -> str:
    h = http_head(f"{INTACT_BASE}/psimitab/intact-micluster.txt")
    return h.get("Last-Modified", "current").strip() or "current"


def intact_download() -> Optional[Path]:
    out = SOURCES_DIR / "IntAct"
    out.mkdir(parents=True, exist_ok=True)
    url = f"{INTACT_BASE}/psimitab/intact-micluster.txt"
    dest = out / "intact-micluster.txt"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        logger.info("IntAct cached at %s (%.1f MB)",
                    dest, dest.stat().st_size / (1 << 20))
    else:
        logger.info("Downloading IntAct micluster (this can take a few minutes) ...")
        try:
            http_download(url, dest, timeout=1800)
        except Exception as e:
            logger.error("IntAct download failed: %s", e)
            return None
    update_version("intact", version=intact_latest_version(), url=url, path=dest)
    return dest


_UNIPROT_RE = re.compile(r"uniprotkb:([A-Z0-9]+)")
_PMID_RE = re.compile(r"pubmed:(\d+)")


def intact_iter_human(path: Path, *, max_rows: int = 0) -> Iterable[dict]:
    """Stream PSI-MITAB rows filtered to human-human."""
    with open_text(path) as fp:
        header = fp.readline().rstrip("\r\n").split("\t")
        n = 0
        for line in fp:
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < 14:
                continue
            taxid_a = cells[9]
            taxid_b = cells[10]
            if "9606" not in taxid_a or "9606" not in taxid_b:
                continue
            up_a_m = _UNIPROT_RE.search(cells[0])
            up_b_m = _UNIPROT_RE.search(cells[1])
            if not (up_a_m and up_b_m):
                continue
            pmid_m = _PMID_RE.search(cells[8])
            try:
                conf = float(cells[14].split(":")[-1]) if len(cells) > 14 and cells[14] else None
            except Exception:
                conf = None
            yield {
                "id_a": up_a_m.group(1),
                "id_b": up_b_m.group(1),
                "id_type": "uniprot",
                "source_id": cells[13] if len(cells) > 13 else "",
                "detection_method": cells[6],
                "evidence_type": "experimental",
                "pubmed_id": pmid_m.group(1) if pmid_m else "",
                "score": str(conf) if conf is not None else "",
            }
            n += 1
            if max_rows and n >= max_rows:
                return


# ---------------------------------------------------------------------------
# HuRI
# ---------------------------------------------------------------------------

HURI_FILES = ["HuRI.tsv", "HI-union.tsv", "Lit-BM.tsv"]


def huri_download() -> "list[Path]":
    out = SOURCES_DIR / "HuRI"
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for f in HURI_FILES:
        url = f"{HURI_BASE}/{f}"
        dest = out / f
        if dest.exists() and dest.stat().st_size > 1000:
            logger.info("HuRI %s cached", f)
        else:
            try:
                http_download(url, dest, timeout=300)
            except Exception as e:
                logger.warning("HuRI %s failed: %s", f, e)
                continue
        paths.append(dest)
    if paths:
        first = paths[0]
        h = http_head(f"{HURI_BASE}/{HURI_FILES[0]}")
        update_version("huri", version=h.get("Last-Modified", "static"),
                       url=f"{HURI_BASE}/HuRI.tsv", path=first)
    return paths


def huri_iter(path: Path) -> Iterable[dict]:
    """HuRI files are 2-column Ensembl-gene-id TSVs."""
    with path.open("r") as fp:
        for line in fp:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0].strip(), parts[1].strip()
            if not a or not b or a == "Symbol":
                continue
            yield {
                "id_a": a, "id_b": b,
                "id_type": "ensembl_gene",
                "evidence_type": "experimental",
                "detection_method": "Y2H",
            }


# ---------------------------------------------------------------------------
# Reactome
# ---------------------------------------------------------------------------

REACTOME_FILES = {
    "pathways":      "ReactomePathways.txt",
    "uniprot2path":  "UniProt2Reactome.txt",
    "ensembl2path":  "Ensembl2Reactome.txt",
    "interactions":  "interactors/reactome.homo_sapiens.interactions.tab-delimited.txt",
}


def reactome_download() -> dict:
    out = SOURCES_DIR / "Reactome"
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for k, f in REACTOME_FILES.items():
        url = f"{REACTOME_BASE}/{f}"
        dest = out / Path(f).name
        if dest.exists() and dest.stat().st_size > 0:
            logger.info("Reactome %s cached", dest.name)
        else:
            try:
                http_download(url, dest, timeout=300)
            except Exception as e:
                logger.warning("Reactome %s failed: %s", dest.name, e)
                continue
        paths[k] = dest
    if paths:
        h = http_head(f"{REACTOME_BASE}/{REACTOME_FILES['pathways']}")
        update_version("reactome", version=h.get("Last-Modified", "current"),
                       url=f"{REACTOME_BASE}/{REACTOME_FILES['pathways']}",
                       path=list(paths.values())[0])
    return paths


def reactome_iter_pathways(path: Path) -> Iterable[dict]:
    """ReactomePathways.txt: pathway_id <tab> name <tab> organism."""
    with path.open("r") as fp:
        for line in fp:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3:
                continue
            pid, name, organism = parts[:3]
            if "Homo sapiens" not in organism:
                continue
            yield {"pathway_id": pid, "name": name, "organism": organism,
                   "source": "reactome"}


def reactome_iter_membership(path: Path) -> Iterable[dict]:
    """UniProt2Reactome.txt: UniProt <tab> ReactomeID <tab> URL <tab> name <tab>
       evidence <tab> organism."""
    with path.open("r") as fp:
        for line in fp:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 6:
                continue
            up, pid, _, name, _, organism = parts[:6]
            if "Homo sapiens" not in organism:
                continue
            yield {"pathway_id": pid, "uniprot_ac": up.strip(), "name": name,
                   "source": "reactome"}


def reactome_iter_interactions(path: Path) -> Iterable[dict]:
    """reactome.homo_sapiens.interactions.tab-delimited.txt — UniProt-keyed
    pathway-derived interactions for human."""
    with path.open("r") as fp:
        header = fp.readline().rstrip("\r\n").split("\t")
        idx = {c.strip(): i for i, c in enumerate(header)}
        a_col = next((idx[k] for k in idx if "Interactor 1 uniprot id" in k), 0)
        b_col = next((idx[k] for k in idx if "Interactor 2 uniprot id" in k), 1)
        for line in fp:
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < max(a_col, b_col) + 1:
                continue
            a, b = cells[a_col].strip(), cells[b_col].strip()
            a = a.split(":")[-1] if a else a
            b = b.split(":")[-1] if b else b
            if not a or not b:
                continue
            yield {"id_a": a, "id_b": b, "id_type": "uniprot",
                   "evidence_type": "curated",
                   "detection_method": "pathway-derived"}


# ---------------------------------------------------------------------------
# KEGG (REST API)
# ---------------------------------------------------------------------------


def kegg_download(*, max_pathways: int = 400, throttle: float = 0.4) -> Path:
    out = SOURCES_DIR / "KEGG"
    out.mkdir(parents=True, exist_ok=True)
    s, body, _ = http_request(f"{KEGG_BASE}/list/pathway/hsa")
    if s != 200:
        logger.error("KEGG list/pathway/hsa HTTP %s", s)
        return out
    pathways = []
    for line in body.decode("utf-8", "replace").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pid = parts[0].strip()
            name = parts[1].strip()
            pathways.append({"pathway_id": pid, "name": name,
                             "source": "kegg", "organism": "Homo sapiens"})
    (out / "pathways_hsa.json").write_text(json.dumps(pathways, indent=2))
    membership = []
    pathways_to_pull = pathways[:max_pathways]
    for i, p in enumerate(pathways_to_pull):
        time.sleep(throttle)
        s, body, _ = http_request(
            f"{KEGG_BASE}/link/hsa/{p['pathway_id']}", timeout=30)
        if s != 200:
            continue
        for line in body.decode("utf-8", "replace").splitlines():
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = parts[1].strip()
            membership.append({"pathway_id": p["pathway_id"],
                                "gene_id": gene, "source": "kegg"})
        if (i + 1) % 25 == 0:
            logger.info("KEGG %d/%d pathways", i + 1, len(pathways_to_pull))
    (out / "pathway_membership_hsa.json").write_text(
        json.dumps(membership, indent=2))
    update_version("kegg", version=time.strftime("%Y-%m-%d"),
                   url=f"{KEGG_BASE}/list/pathway/hsa",
                   path=out / "pathways_hsa.json",
                   n_records=len(membership))
    return out


def kegg_iter_pathways(p: Path) -> Iterable[dict]:
    f = p / "pathways_hsa.json"
    if f.exists():
        for d in json.loads(f.read_text()):
            yield d


def kegg_iter_membership(p: Path) -> Iterable[dict]:
    f = p / "pathway_membership_hsa.json"
    if f.exists():
        for d in json.loads(f.read_text()):
            yield d


# ---------------------------------------------------------------------------
# UniProt id-mapping (optional ID harmonization)
# ---------------------------------------------------------------------------


def uniprot_idmap_download() -> Optional[Path]:
    out = SOURCES_DIR / "UniProt"
    out.mkdir(parents=True, exist_ok=True)
    fname = "HUMAN_9606_idmapping_selected.tab.gz"
    dest = out / fname
    url = f"{UNIPROT_IDMAP_BASE}/{fname}"
    if dest.exists() and dest.stat().st_size > 100_000:
        logger.info("UniProt idmap cached")
    else:
        logger.info("Downloading UniProt idmap (~80MB) ...")
        try:
            http_download(url, dest, timeout=1800)
        except Exception as e:
            logger.error("UniProt idmap download failed: %s", e)
            return None
    h = http_head(url)
    update_version("uniprot_idmap", version=h.get("Last-Modified", "current"),
                   url=url, path=dest)
    return dest


def uniprot_idmap_iter(path: Path, *, max_rows: int = 0) -> Iterable[dict]:
    """idmapping_selected.tab columns: 0:UniProtKB-AC 1:UniProtKB-ID 2:GeneID
    3:RefSeq 4:GI 5:PDB 6:GO 7:UniRef100 ... 18:Ensembl 19:Ensembl_TRS
    20:Ensembl_PRO 21:Additional 22:Alt_taxon 23:Tax_ID."""
    with open_text(path) as fp:
        n = 0
        for line in fp:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 19:
                continue
            uac = parts[0].strip()
            entrez = parts[2].split(";")[0].strip() if parts[2] else ""
            ensembl = parts[18].split(";")[0].strip() if parts[18] else ""
            yield {"uniprot_ac": uac, "entrez": entrez, "ensembl": ensembl}
            n += 1
            if max_rows and n >= max_rows:
                return


# ---------------------------------------------------------------------------
# IGVF Portal — protein assays + PPI files
# ---------------------------------------------------------------------------


def _portal_get(path: str) -> dict:
    url = f"{PORTAL_API_BASE}{path}"
    s, body, _ = http_request(url, timeout=60,
                              headers={"Accept": "application/json"})
    if s != 200:
        logger.warning("Portal HTTP %s on %s", s, path)
        return {}
    try:
        return json.loads(body)
    except Exception:
        return {}


def igvf_protein_download() -> dict:
    """Pull the 214 protein-slim MeasurementSets, the 14 PPI score files,
    1 stability file, and the protein reference files."""
    out = SOURCES_DIR / "IGVF"
    out.mkdir(parents=True, exist_ok=True)
    files_dir = out / "ppi_files"
    files_dir.mkdir(parents=True, exist_ok=True)

    summary: "dict[str, Any]" = {}
    # 1) Protein-slim MeasurementSets, full graph
    ms = _portal_get("/search/?type=MeasurementSet&assay_slims=protein"
                     "&format=json&limit=500")
    sets = ms.get("@graph", []) or []
    summary["measurement_sets"] = sets
    (out / "measurement_sets.json").write_text(json.dumps(sets, indent=2))
    logger.info("IGVF: %d protein-slim MeasurementSets", len(sets))

    # 2) PPI-score files
    ppi = _portal_get("/search/?content_type=protein+to+protein+interaction+score"
                      "&format=json&limit=200")
    files = ppi.get("@graph", []) or []
    summary["ppi_files_meta"] = files
    (out / "ppi_files_meta.json").write_text(json.dumps(files, indent=2))
    logger.info("IGVF: %d PPI-score files", len(files))

    # 3) Protein-stability fluorescence file
    sf = _portal_get("/search/?content_type=protein+stability+fluorescence+score"
                      "&format=json&limit=20")
    sfiles = sf.get("@graph", []) or []
    summary["stability_files_meta"] = sfiles
    (out / "stability_files_meta.json").write_text(json.dumps(sfiles, indent=2))

    # 4) Protein reference files
    pr = _portal_get("/search/?content_type=proteins&format=json&limit=50")
    refs = pr.get("@graph", []) or []
    summary["protein_reference_meta"] = refs
    (out / "protein_reference_meta.json").write_text(json.dumps(refs, indent=2))

    # 5) Protein language model weights
    plm = _portal_get("/search/?content_type=protein+language+model"
                       "&format=json&limit=50")
    plms = plm.get("@graph", []) or []
    summary["protein_language_models"] = plms
    (out / "protein_language_models.json").write_text(json.dumps(plms, indent=2))

    # 6) Pull the actual PPI / stability files (small TSVs/CSVs).
    # The search summary doesn't include `s3_uri` / `href`, so we fetch
    # each file's detail object first and then pull from the public S3
    # bucket (igvf-public is world-readable).
    pulled = []
    for f in files + sfiles:
        access = f.get("accession") or "FILE"
        atid = f.get("@id") or f"/files/{access}/"
        detail = _portal_get(f"{atid}?format=json")
        s3_uri = detail.get("s3_uri", "")
        href = detail.get("href", "")
        # Convert s3://igvf-public/<key> -> public HTTPS URL
        dl = ""
        if s3_uri.startswith("s3://igvf-public/"):
            key = s3_uri[len("s3://igvf-public/"):]
            dl = f"https://igvf-public.s3.amazonaws.com/{key}"
        elif href:
            dl = f"{PORTAL_API_BASE}{href}"
        if not dl:
            logger.warning("IGVF file %s: no download URL", access)
            continue
        # Match the actual extension on the upstream file
        ext = (Path(s3_uri or href).suffix
               or "." + (detail.get("file_format") or "tsv"))
        if (s3_uri or href).endswith(".gz") and not ext.endswith(".gz"):
            ext = ext + ".gz"
        local = files_dir / f"{access}{ext}"
        if local.exists() and local.stat().st_size > 0:
            pulled.append(str(local))
            continue
        try:
            http_download(dl, local, timeout=300)
            pulled.append(str(local))
        except Exception as e:
            logger.warning("IGVF file %s download failed (%s): %s",
                           access, dl, e)
    summary["pulled_files"] = pulled
    update_version("igvf", version=time.strftime("%Y-%m-%d"),
                   url=PORTAL_API_BASE,
                   path=out / "measurement_sets.json",
                   n_records=len(sets))
    return summary


# ---------------------------------------------------------------------------
# SQLite KG
# ---------------------------------------------------------------------------

KG_SCHEMA = """
CREATE TABLE IF NOT EXISTS proteins (
  pid TEXT PRIMARY KEY,
  id_type TEXT NOT NULL,
  symbol TEXT,
  organism_taxid INTEGER DEFAULT 9606,
  source_set TEXT
);
CREATE TABLE IF NOT EXISTS interactions (
  rowid INTEGER PRIMARY KEY AUTOINCREMENT,
  id_a TEXT NOT NULL,
  id_b TEXT NOT NULL,
  id_type TEXT NOT NULL,
  source TEXT NOT NULL,
  source_id TEXT NOT NULL DEFAULT '',
  detection_method TEXT,
  evidence_type TEXT,
  pubmed_id TEXT,
  confidence_score REAL,
  taxon INTEGER DEFAULT 9606,
  UNIQUE(id_a, id_b, source, source_id)
);
CREATE INDEX IF NOT EXISTS idx_int_a ON interactions(id_a);
CREATE INDEX IF NOT EXISTS idx_int_b ON interactions(id_b);
CREATE INDEX IF NOT EXISTS idx_int_src ON interactions(source);
CREATE INDEX IF NOT EXISTS idx_int_type ON interactions(id_type);
CREATE TABLE IF NOT EXISTS pathways (
  pathway_id TEXT PRIMARY KEY,
  source TEXT NOT NULL,
  name TEXT,
  organism_taxid INTEGER DEFAULT 9606
);
CREATE TABLE IF NOT EXISTS pathway_membership (
  pathway_id TEXT NOT NULL,
  member_id TEXT NOT NULL,
  id_type TEXT NOT NULL,
  source TEXT NOT NULL,
  UNIQUE(pathway_id, member_id, source)
);
CREATE INDEX IF NOT EXISTS idx_pm_path ON pathway_membership(pathway_id);
CREATE INDEX IF NOT EXISTS idx_pm_mem ON pathway_membership(member_id);
CREATE TABLE IF NOT EXISTS igvf_evidence (
  file_accession TEXT PRIMARY KEY,
  dataset_accession TEXT,
  assay_title TEXT,
  lab TEXT,
  content_type TEXT,
  description TEXT,
  url TEXT
);
CREATE TABLE IF NOT EXISTS id_map (
  uniprot_ac TEXT,
  entrez TEXT,
  ensembl TEXT,
  symbol TEXT,
  PRIMARY KEY (uniprot_ac)
);
CREATE INDEX IF NOT EXISTS idx_idmap_entrez ON id_map(entrez);
CREATE INDEX IF NOT EXISTS idx_idmap_ensembl ON id_map(ensembl);
CREATE INDEX IF NOT EXISTS idx_idmap_symbol ON id_map(symbol);
CREATE TABLE IF NOT EXISTS versions (
  source TEXT PRIMARY KEY,
  version TEXT,
  downloaded_at TEXT,
  url TEXT,
  sha256 TEXT,
  n_records INTEGER
);
"""


def open_kg() -> sqlite3.Connection:
    KG_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(KG_PATH))
    conn.executescript(KG_SCHEMA)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    return conn


def _ingest_interactions(conn: sqlite3.Connection, rows: Iterable[dict],
                         source: str, *, batch: int = 5000) -> int:
    cur = conn.cursor()
    n = 0
    buf = []
    for r in rows:
        a = (r.get("id_a") or "").strip()
        b = (r.get("id_b") or "").strip()
        if not a or not b:
            continue
        if a > b:
            a, b = b, a
        buf.append((
            a, b, r.get("id_type", "symbol"), source,
            r.get("source_id", "") or "",
            r.get("detection_method", "") or "",
            r.get("evidence_type", "") or "",
            r.get("pubmed_id", "") or "",
            float(r["score"]) if r.get("score") not in (None, "", "-") and
            _isfloat(r.get("score")) else None,
        ))
        if len(buf) >= batch:
            cur.executemany("""
              INSERT OR IGNORE INTO interactions
                (id_a,id_b,id_type,source,source_id,detection_method,
                 evidence_type,pubmed_id,confidence_score)
              VALUES (?,?,?,?,?,?,?,?,?)
            """, buf)
            n += len(buf)
            buf.clear()
            conn.commit()
    if buf:
        cur.executemany("""
          INSERT OR IGNORE INTO interactions
            (id_a,id_b,id_type,source,source_id,detection_method,
             evidence_type,pubmed_id,confidence_score)
          VALUES (?,?,?,?,?,?,?,?,?)
        """, buf)
        n += len(buf)
    conn.commit()
    return n


def _isfloat(s) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False


def _ingest_pathways(conn: sqlite3.Connection, rows: Iterable[dict]) -> int:
    cur = conn.cursor()
    n = 0
    for r in rows:
        cur.execute("""
          INSERT OR REPLACE INTO pathways(pathway_id,source,name)
          VALUES (?,?,?)
        """, (r["pathway_id"], r.get("source", "?"), r.get("name", "")))
        n += 1
    conn.commit()
    return n


def _ingest_pathway_membership(conn: sqlite3.Connection,
                               rows: Iterable[dict],
                               id_field: str, id_type: str) -> int:
    cur = conn.cursor()
    n = 0
    buf = []
    for r in rows:
        mid = r.get(id_field) or r.get("uniprot_ac") or r.get("gene_id")
        if not mid:
            continue
        buf.append((r["pathway_id"], mid, id_type, r.get("source", "?")))
        if len(buf) >= 5000:
            cur.executemany("""
              INSERT OR IGNORE INTO pathway_membership
                (pathway_id,member_id,id_type,source)
              VALUES (?,?,?,?)
            """, buf)
            n += len(buf)
            buf.clear()
            conn.commit()
    if buf:
        cur.executemany("""
          INSERT OR IGNORE INTO pathway_membership
            (pathway_id,member_id,id_type,source)
          VALUES (?,?,?,?)
        """, buf)
        n += len(buf)
    conn.commit()
    return n


def _ingest_idmap(conn: sqlite3.Connection, rows: Iterable[dict]) -> int:
    cur = conn.cursor()
    n = 0
    buf = []
    for r in rows:
        buf.append((r["uniprot_ac"], r.get("entrez", ""),
                    r.get("ensembl", ""), r.get("symbol", "")))
        if len(buf) >= 5000:
            cur.executemany("""
              INSERT OR REPLACE INTO id_map(uniprot_ac,entrez,ensembl,symbol)
              VALUES (?,?,?,?)
            """, buf)
            n += len(buf); buf.clear(); conn.commit()
    if buf:
        cur.executemany("""
          INSERT OR REPLACE INTO id_map(uniprot_ac,entrez,ensembl,symbol)
          VALUES (?,?,?,?)
        """, buf)
        n += len(buf)
    conn.commit()
    return n


def _ingest_igvf_evidence(conn: sqlite3.Connection, igvf_dir: Path) -> int:
    cur = conn.cursor()
    n = 0
    for jf in ("ppi_files_meta.json", "stability_files_meta.json",
               "protein_reference_meta.json", "protein_language_models.json"):
        p = igvf_dir / jf
        if not p.exists():
            continue
        try:
            files = json.loads(p.read_text())
        except Exception:
            continue
        for f in files:
            cur.execute("""
              INSERT OR REPLACE INTO igvf_evidence
                (file_accession,dataset_accession,assay_title,lab,content_type,
                 description,url)
              VALUES (?,?,?,?,?,?,?)
            """, (
                f.get("accession") or "",
                (f.get("file_set", {}) or {}).get("accession", ""),
                ", ".join((f.get("file_set", {}) or {}).get("preferred_assay_titles", []) or []),
                ((f.get("lab") or {}).get("title", "") if isinstance(f.get("lab"), dict) else f.get("lab", "")) or "",
                f.get("content_type", "") or "",
                f.get("description", "") or "",
                f.get("@id", "") or "",
            ))
            n += 1
    conn.commit()
    return n


# ---------------------------------------------------------------------------
# Build-KG orchestrator
# ---------------------------------------------------------------------------


def build_kg(*, sources: "list[str]", max_rows: int = 0) -> dict:
    mkdirs()
    conn = open_kg()
    counts = {}
    try:
        if "biogrid" in sources:
            zp = next((SOURCES_DIR / "BioGRID").glob("BIOGRID-ALL-*.tab3.zip"), None)
            if zp:
                logger.info("Ingest BioGRID %s ...", zp.name)
                counts["biogrid"] = _ingest_interactions(
                    conn, biogrid_iter_human(zp, max_rows=max_rows), "biogrid")
            else:
                logger.warning("BioGRID zip not present; run 'download biogrid' first")

        if "intact" in sources:
            ip = SOURCES_DIR / "IntAct" / "intact-micluster.txt"
            if ip.exists():
                logger.info("Ingest IntAct ...")
                counts["intact"] = _ingest_interactions(
                    conn, intact_iter_human(ip, max_rows=max_rows), "intact")
            else:
                logger.warning("IntAct not present")

        if "huri" in sources:
            for fname in HURI_FILES:
                hp = SOURCES_DIR / "HuRI" / fname
                if hp.exists():
                    logger.info("Ingest HuRI %s ...", fname)
                    counts[f"huri:{fname}"] = _ingest_interactions(
                        conn, huri_iter(hp), f"huri:{fname.replace('.tsv','')}")

        if "reactome" in sources:
            rd = SOURCES_DIR / "Reactome"
            if (rd / "ReactomePathways.txt").exists():
                _ingest_pathways(conn, reactome_iter_pathways(
                    rd / "ReactomePathways.txt"))
            if (rd / "UniProt2Reactome.txt").exists():
                _ingest_pathway_membership(
                    conn, reactome_iter_membership(rd / "UniProt2Reactome.txt"),
                    id_field="uniprot_ac", id_type="uniprot")
            if (rd / "reactome.homo_sapiens.interactions.tab-delimited.txt").exists():
                counts["reactome"] = _ingest_interactions(
                    conn, reactome_iter_interactions(
                        rd / "reactome.homo_sapiens.interactions.tab-delimited.txt"),
                    "reactome")

        if "kegg" in sources:
            kd = SOURCES_DIR / "KEGG"
            if (kd / "pathways_hsa.json").exists():
                _ingest_pathways(conn, kegg_iter_pathways(kd))
                _ingest_pathway_membership(
                    conn, kegg_iter_membership(kd),
                    id_field="gene_id", id_type="kegg_gene")

        if "uniprot" in sources:
            up = SOURCES_DIR / "UniProt" / "HUMAN_9606_idmapping_selected.tab.gz"
            if up.exists():
                logger.info("Ingest UniProt idmap ...")
                counts["uniprot_idmap"] = _ingest_idmap(
                    conn, uniprot_idmap_iter(up, max_rows=max_rows))

        if "igvf" in sources:
            ip = SOURCES_DIR / "IGVF"
            if ip.exists():
                counts["igvf"] = _ingest_igvf_evidence(conn, ip)
    finally:
        conn.commit()
        conn.close()
    return counts


# ---------------------------------------------------------------------------
# Stats and visualization
# ---------------------------------------------------------------------------


def kg_stats() -> dict:
    if not KG_PATH.exists():
        return {"error": "KG not built yet — run proteomics build-kg first."}
    conn = sqlite3.connect(str(KG_PATH))
    try:
        c = conn.cursor()
        out: "dict[str, Any]" = {}
        out["interactions_total"] = c.execute(
            "SELECT COUNT(*) FROM interactions").fetchone()[0]
        out["proteins_distinct"] = c.execute(
            "SELECT COUNT(DISTINCT id) FROM ("
            "SELECT id_a AS id FROM interactions UNION "
            "SELECT id_b AS id FROM interactions)").fetchone()[0]
        out["pathways_total"] = c.execute(
            "SELECT COUNT(*) FROM pathways").fetchone()[0]
        out["pathway_membership"] = c.execute(
            "SELECT COUNT(*) FROM pathway_membership").fetchone()[0]
        out["igvf_evidence"] = c.execute(
            "SELECT COUNT(*) FROM igvf_evidence").fetchone()[0]
        out["per_source"] = dict(c.execute(
            "SELECT source, COUNT(*) FROM interactions GROUP BY source"
        ).fetchall())
        out["per_id_type"] = dict(c.execute(
            "SELECT id_type, COUNT(*) FROM interactions GROUP BY id_type"
        ).fetchall())
        out["per_evidence_type"] = dict(c.execute(
            "SELECT evidence_type, COUNT(*) FROM interactions "
            "GROUP BY evidence_type"
        ).fetchall())
        out["per_detection_method"] = dict(c.execute(
            "SELECT detection_method, COUNT(*) FROM interactions "
            "GROUP BY detection_method "
            "ORDER BY 2 DESC LIMIT 30"
        ).fetchall())
        out["top_hubs"] = c.execute("""
            SELECT id, COUNT(*) AS deg FROM (
              SELECT id_a AS id FROM interactions UNION ALL
              SELECT id_b AS id FROM interactions
            ) GROUP BY id ORDER BY deg DESC LIMIT 30
        """).fetchall()
        return out
    finally:
        conn.close()


def kg_visualize(label: str, *, gene: Optional[str] = None,
                 ego_max_neighbors: int = 60) -> Path:
    """Generate degree distribution, top-hubs bar plot, optional ego graph."""
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        raise RuntimeError(f"matplotlib required for kg-visualize: {e}")
    if not KG_PATH.exists():
        raise RuntimeError("KG not built yet — run build-kg first.")
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}"
    plots = out / "Plots"
    plots.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(KG_PATH))
    c = conn.cursor()

    # Degree distribution -----------------------------------------------------
    rows = c.execute("""
        SELECT id, COUNT(*) AS deg FROM (
          SELECT id_a AS id FROM interactions UNION ALL
          SELECT id_b AS id FROM interactions
        ) GROUP BY id
    """).fetchall()
    if rows:
        degs = [r[1] for r in rows]
        cnt = Counter(degs)
        x = sorted(cnt.keys())
        y = [cnt[k] for k in x]
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.scatter(x, y, s=12, alpha=0.7, c="#2C7FB8")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Degree (log)"); ax.set_ylabel("Number of proteins (log)")
        ax.set_title(f"PPI-KG degree distribution (n={len(rows):,} proteins)")
        ax.grid(alpha=0.3, linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "degree_distribution.png", dpi=150)
        plt.close(fig)

    # Top hubs ----------------------------------------------------------------
    top = c.execute("""
        SELECT id, COUNT(*) AS deg FROM (
          SELECT id_a AS id FROM interactions UNION ALL
          SELECT id_b AS id FROM interactions
        ) GROUP BY id ORDER BY deg DESC LIMIT 30
    """).fetchall()
    if top:
        names = [r[0][:18] for r in top][::-1]
        vals = [r[1] for r in top][::-1]
        fig, ax = plt.subplots(figsize=(7, 8))
        ax.barh(range(len(names)), vals, color="#E6550D")
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("Degree (interaction partners)")
        ax.set_title("Top-30 PPI-KG hub proteins")
        fig.tight_layout()
        fig.savefig(plots / "top_hubs.png", dpi=150)
        plt.close(fig)

    # Per-source breakdown ----------------------------------------------------
    src = c.execute(
        "SELECT source, COUNT(*) FROM interactions GROUP BY source "
        "ORDER BY 2 DESC").fetchall()
    if src:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.bar(range(len(src)), [r[1] for r in src],
                color="#31A354")
        ax.set_xticks(range(len(src)))
        ax.set_xticklabels([r[0] for r in src], rotation=30, ha="right",
                            fontsize=8)
        ax.set_ylabel("Interactions"); ax.set_yscale("log")
        ax.set_title("Interactions per source")
        fig.tight_layout()
        fig.savefig(plots / "per_source.png", dpi=150)
        plt.close(fig)

    # Ego graph ---------------------------------------------------------------
    ego_path = None
    if gene:
        nbrs = c.execute("""
            SELECT id_b AS id FROM interactions WHERE id_a=?
            UNION
            SELECT id_a AS id FROM interactions WHERE id_b=?
        """, (gene, gene)).fetchall()
        nbrs = [r[0] for r in nbrs][:ego_max_neighbors]
        if nbrs:
            try:
                import networkx as nx  # type: ignore
                G = nx.Graph()
                for nb in nbrs:
                    G.add_edge(gene, nb)
                # Add edges among neighbors so we can see clusters
                if len(nbrs) > 1:
                    placeholders = ",".join("?" * len(nbrs))
                    inner = c.execute(f"""
                        SELECT id_a, id_b FROM interactions
                        WHERE id_a IN ({placeholders}) AND id_b IN ({placeholders})
                    """, nbrs + nbrs).fetchall()
                    for a, b in inner:
                        G.add_edge(a, b)
                pos = nx.spring_layout(G, k=0.6, seed=7)
                fig, ax = plt.subplots(figsize=(8, 8))
                node_colors = ["#E6550D" if n == gene else "#3182BD"
                                for n in G.nodes()]
                node_sizes = [600 if n == gene else 120 for n in G.nodes()]
                nx.draw_networkx_edges(G, pos, alpha=0.35, ax=ax)
                nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                                        node_size=node_sizes, ax=ax,
                                        edgecolors="white", linewidths=0.6)
                labels = {n: n for n in G.nodes() if n == gene
                           or G.degree(n) >= 3}
                nx.draw_networkx_labels(G, pos, labels=labels, font_size=7,
                                        ax=ax)
                ax.set_title(f"Ego graph: {gene}  (k={len(nbrs)} neighbors)")
                ax.axis("off")
                fig.tight_layout()
                ego_path = plots / f"ego_{safe_label(gene)}.png"
                fig.savefig(ego_path, dpi=150)
                plt.close(fig)
            except Exception as e:
                logger.warning("Ego graph skipped (networkx?): %s", e)
    conn.close()

    # Markdown summary --------------------------------------------------------
    md = ["# Proteomics PPI-KG visualization", "",
          f"- Output dir: `{out}`",
          f"- KG: `{KG_PATH}`", "",
          "## Plots",
          "- `Plots/degree_distribution.png` — log-log degree distribution",
          "- `Plots/top_hubs.png` — top-30 hub proteins by degree",
          "- `Plots/per_source.png` — interaction count per source database"]
    if ego_path:
        md.append(f"- `Plots/{ego_path.name}` — ego graph for `{gene}`")
    (out / "report.md").write_text("\n".join(md) + "\n")
    return out


# ---------------------------------------------------------------------------
# Per-assay example figures (from real IGVF Portal files)
# ---------------------------------------------------------------------------


_PREFERRED_SCORE_PATTERNS = (
    # Y2H growth / fluorescence readouts
    "per_well_integrated_intensity",
    "integrated_intensity",
    "manual_score_growth",
    "colony_area",
    "fluorescence",
    "intensity",
    "FCAS",
    # MAVE / VAMP-seq
    "abundance", "stability_score", "score", "fitness",
    "log2FoldChange", "logFC", "vamp_score", "abundance_score",
    # General
    "value", "effect", "dms_score",
)


def _read_score_column(path: Path) -> "list[float]":
    """Read the most informative numeric column from a TSV/CSV (gz-aware).

    Strategy:
      1. Sniff the header. For each cell, if its lower-cased name contains
         one of `_PREFERRED_SCORE_PATTERNS`, mark as candidate.
      2. Materialize a sample (up to 50k rows) and pick the candidate
         with the largest coefficient of variation (std / |mean|). This
         skips constant columns like assay_id.
      3. Fall back to the highest-variance numeric column overall.
    """
    name = path.name.lower()
    delim = "," if name.endswith(".csv") or name.endswith(".csv.gz") else "\t"
    try:
        with open_text(path) as fp:
            reader = csv.reader(fp, delimiter=delim)
            try:
                header = next(reader)
            except StopIteration:
                return []
            preferred = []
            for i, h in enumerate(header):
                hl = h.strip().lower()
                if any(pat.lower() in hl for pat in _PREFERRED_SCORE_PATTERNS):
                    preferred.append(i)
            # Materialize a sample for variance sniffing
            cols: "dict[int, list[float]]" = {}
            n_rows = 0
            sample_cap = 50000
            for row in reader:
                for i, v in enumerate(row):
                    try:
                        x = float(v)
                    except Exception:
                        continue
                    cols.setdefault(i, []).append(x)
                n_rows += 1
                if n_rows >= sample_cap:
                    break
            if not cols:
                return []
            def _cv(vals):
                if len(vals) < 2:
                    return 0.0
                m = sum(vals) / len(vals)
                v = sum((x - m) ** 2 for x in vals) / max(1, len(vals) - 1)
                std = v ** 0.5
                if abs(m) < 1e-12:
                    return std
                return std / abs(m)
            # Choose: best preferred-col by CV; else best overall
            scope = preferred if preferred else list(cols.keys())
            best = max(scope, key=lambda i: _cv(cols.get(i, [])))
            sample = cols[best]
            # If we hit the cap, keep streaming the rest of the file for that col
            if n_rows >= sample_cap:
                with open_text(path) as fp2:
                    r2 = csv.reader(fp2, delimiter=delim)
                    next(r2, None)
                    for j, row in enumerate(r2):
                        if j < sample_cap:
                            continue
                        if len(row) > best:
                            try:
                                sample.append(float(row[best]))
                            except Exception:
                                pass
            return sample
    except Exception as e:
        logger.warning("read_score_column %s: %s", path, e)
        return []


ASSAY_FIGURE_SPECS = [
    {"assay": "VAMP-seq (MultiSTEP)", "kind": "stability",
     "color": "#54278F",
     "title": "VAMP-seq (MultiSTEP) abundance score distribution"},
    {"assay": "VAMP-seq", "kind": "stability",
     "color": "#7570B3",
     "title": "VAMP-seq abundance score distribution"},
    {"assay": "MAVE", "kind": "fitness",
     "color": "#1B9E77",
     "title": "MAVE variant-effect score distribution"},
    {"assay": "Arrayed semi-qY2H v1", "kind": "ppi",
     "color": "#D95F02",
     "title": "semi-qY2H v1 interaction-strength distribution"},
    {"assay": "Arrayed semi-qY2H v2", "kind": "ppi",
     "color": "#E7298A",
     "title": "semi-qY2H v2 interaction-strength distribution"},
    {"assay": "Arrayed semi-qY2H v3", "kind": "ppi",
     "color": "#A6761D",
     "title": "semi-qY2H v3 interaction-strength distribution"},
    {"assay": "DUAL-IPA", "kind": "ppi",
     "color": "#666666",
     "title": "DUAL-IPA fluorescence readout distribution"},
]


def _representative_file_for_assay(assay: str) -> Optional[Path]:
    """Pick a representative IGVF-pulled file for the assay term.

    Prefers `preferred_assay_titles` exact match on the file metadata,
    falls back to description keyword match.
    """
    igvf = SOURCES_DIR / "IGVF"
    files_dir = igvf / "ppi_files"
    if not files_dir.exists():
        return None
    meta_files = []
    for jf in ("ppi_files_meta.json", "stability_files_meta.json"):
        p = igvf / jf
        if p.exists():
            try:
                meta_files.extend(json.loads(p.read_text()))
            except Exception:
                pass

    # 1) Exact preferred_assay_titles match
    for f in meta_files:
        titles = f.get("preferred_assay_titles") or []
        if assay in titles:
            access = f.get("accession", "")
            for p in files_dir.iterdir():
                if access and p.name.startswith(access):
                    return p

    # 2) Description-based fallback for VAMP-seq family / MAVE
    asl = assay.lower()
    for f in meta_files:
        desc = (f.get("description") or "").lower()
        if asl in desc:
            access = f.get("accession", "")
            for p in files_dir.iterdir():
                if access and p.name.startswith(access):
                    return p
    return None


def assay_figures(label: str) -> Path:
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        raise RuntimeError(f"matplotlib required: {e}")
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}_assay_figures"
    plots = out / "Plots"
    plots.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for spec in ASSAY_FIGURE_SPECS:
        f = _representative_file_for_assay(spec["assay"])
        if not f:
            logger.info("No local file found for %s", spec["assay"])
            summary_rows.append({"assay": spec["assay"], "file": None,
                                  "n_values": 0})
            continue
        vals = _read_score_column(f)
        if not vals:
            logger.info("No numeric column in %s for %s", f.name, spec["assay"])
            summary_rows.append({"assay": spec["assay"], "file": f.name,
                                  "n_values": 0})
            continue
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
        ax.hist(vals, bins=40, color=spec["color"], edgecolor="white",
                 alpha=0.85)
        ax.set_xlabel("Score")
        ax.set_ylabel("Count")
        ax.set_title(f"{spec['title']}\n(IGVF file {f.stem}; n={len(vals)})",
                      fontsize=10)
        ax.grid(alpha=0.25, linestyle=":")
        fig.tight_layout()
        png = plots / f"{safe_label(spec['assay'])}.png"
        fig.savefig(png, dpi=150)
        plt.close(fig)
        summary_rows.append({"assay": spec["assay"], "file": f.name,
                              "n_values": len(vals),
                              "min": min(vals), "max": max(vals),
                              "mean": sum(vals) / len(vals)})

    md = ["# IGVF protein-assay example figures", "",
          f"- Output dir: `{out}`",
          "- All figures derived from the IGVF Portal protein-assay files "
          "pulled by `proteomics igvf-protein`.", "",
          "## Per-assay summary", "",
          "| Assay | File | n | min | mean | max |",
          "|---|---|---|---|---|---|"]
    for r in summary_rows:
        if r["n_values"] == 0:
            md.append(f"| {r['assay']} | {r['file'] or '-'} | 0 | - | - | - |")
        else:
            md.append(f"| {r['assay']} | `{r['file']}` | {r['n_values']} | "
                       f"{r['min']:.3g} | {r['mean']:.3g} | {r['max']:.3g} |")
    (out / "report.md").write_text("\n".join(md) + "\n")
    (out / "summary.json").write_text(json.dumps(summary_rows, indent=2))
    return out


# ---------------------------------------------------------------------------
# Literature survey via the Reference skill
# ---------------------------------------------------------------------------

ASSAY_TERMS = [
    "VAMP-seq MultiSTEP",
    "VAMP-seq",
    "MAVE multiplexed assay variant effect",
    "semi-quantitative yeast two-hybrid",
    "DUAL-IPA",
]


def assay_survey(label: str, *, max_per_assay: int = 20) -> Path:
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}_assay_survey"
    out.mkdir(parents=True, exist_ok=True)
    # We import the reference skill module and call its building blocks.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import reference_skill as ref  # type: ignore
    except Exception as e:
        raise RuntimeError(f"reference_skill import failed: {e}")
    journals = ["Nature", "Cell", "Science",
                 "Nat Genet", "Nat Methods", "Nat Biotechnol",
                 "Nat Commun", "Mol Cell", "Mol Syst Biol"]
    aggregate: "dict[str, list[dict]]" = {}
    for term in ASSAY_TERMS:
        try:
            pubmed = ref.pubmed_search(term, limit=max_per_assay)
        except Exception as e:
            logger.warning("pubmed_search(%s) failed: %s", term, e)
            pubmed = []
        try:
            sem = ref.semanticscholar_search(term, limit=max_per_assay)
        except Exception as e:
            logger.warning("semanticscholar_search(%s) failed: %s", term, e)
            sem = []
        try:
            oa = ref.openalex_search(term, limit=max_per_assay)
        except Exception as e:
            logger.warning("openalex_search(%s) failed: %s", term, e)
            oa = []
        merged = ref.dedup_records(list(pubmed) + list(sem) + list(oa))
        # Filter to high-impact journals
        filt = []
        for r in merged:
            j = (r.get("journal") or r.get("venue") or "")
            if any(jn.lower() in j.lower() for jn in journals):
                filt.append(r)
        if not filt:
            filt = merged[:max_per_assay]
        aggregate[term] = filt[:max_per_assay]
    (out / "literature_survey.json").write_text(json.dumps(aggregate, indent=2))

    md = ["# IGVF protein-assay literature survey", "",
          f"- Output: `{out}/literature_survey.json`",
          f"- Journals filtered to: {', '.join(journals)}", ""]
    for term, recs in aggregate.items():
        md.append(f"## {term}  (n={len(recs)})")
        md.append("")
        for r in recs:
            title = r.get("title", "?")
            year = r.get("year") or r.get("publication_year") or "?"
            j = r.get("journal") or r.get("venue") or "?"
            doi = r.get("doi") or r.get("DOI") or ""
            md.append(f"- **{title}**  ({j}, {year})"
                       + (f"  doi:{doi}" if doi else ""))
        md.append("")
    (out / "literature_survey.md").write_text("\n".join(md) + "\n")
    return out


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_download(args: argparse.Namespace) -> int:
    mkdirs()
    sources = [s.strip().lower() for s in args.source.split(",") if s.strip()]
    if "all" in sources:
        sources = ["biogrid", "intact", "huri", "reactome", "kegg", "igvf"]
    setup_logging("download_" + "_".join(sources))
    for src in sources:
        if src == "biogrid":
            biogrid_download(target_version=args.biogrid_version)
        elif src == "intact":
            intact_download()
        elif src == "huri":
            huri_download()
        elif src == "reactome":
            reactome_download()
        elif src == "kegg":
            kegg_download(max_pathways=args.kegg_max_pathways)
        elif src == "uniprot":
            uniprot_idmap_download()
        elif src == "igvf":
            igvf_protein_download()
        else:
            logger.warning("Unknown source: %s", src)
    print(json.dumps(load_versions(), indent=2))
    return 0


def cmd_versions(_a) -> int:
    mkdirs()
    setup_logging("versions")
    local = load_versions()
    print(json.dumps(local, indent=2))
    # Probe upstream-latest where cheap
    upstream = {}
    try:
        upstream["biogrid"] = biogrid_latest_release()
    except Exception:
        pass
    upstream["intact"] = intact_latest_version()
    print("\nUpstream probes:")
    print(json.dumps(upstream, indent=2))
    return 0


def cmd_update(_a) -> int:
    mkdirs()
    setup_logging("update")
    local = load_versions()
    # BioGRID
    upstream = biogrid_latest_release()
    if upstream and local.get("biogrid", {}).get("version") != upstream:
        logger.info("BioGRID: %s -> %s", local.get("biogrid", {}).get("version"),
                     upstream)
        biogrid_download(target_version=upstream)
    else:
        logger.info("BioGRID: up to date")
    # IntAct
    cur_lm = intact_latest_version()
    if local.get("intact", {}).get("version") != cur_lm:
        logger.info("IntAct: stale, refreshing")
        intact_download()
    # HuRI is mostly static
    huri_download()
    reactome_download()
    print(json.dumps(load_versions(), indent=2))
    return 0


def cmd_igvf_protein(_a) -> int:
    mkdirs()
    setup_logging("igvf_protein")
    summary = igvf_protein_download()
    print(json.dumps({k: (len(v) if isinstance(v, list) else v)
                      for k, v in summary.items()}, indent=2))
    return 0


def cmd_build_kg(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("build_kg")
    sources = [s.strip().lower() for s in args.sources.split(",") if s.strip()]
    if "all" in sources:
        sources = ["biogrid", "intact", "huri", "reactome", "kegg",
                   "uniprot", "igvf"]
    counts = build_kg(sources=sources, max_rows=args.max_rows)
    print(json.dumps(counts, indent=2))
    return 0


def cmd_kg_stats(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("kg_stats")
    s = kg_stats()
    out = DOCS_DIR / f"{timestamp()}_{safe_label(args.label or 'kg_stats')}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "kg_stats.json").write_text(json.dumps(s, indent=2, default=str))
    md = ["# Proteomics PPI-KG summary", "",
          f"- DB: `{KG_PATH}`",
          f"- Total interactions: **{s.get('interactions_total', 0):,}**",
          f"- Distinct proteins: **{s.get('proteins_distinct', 0):,}**",
          f"- Pathways: **{s.get('pathways_total', 0):,}**  "
          f"(membership rows: {s.get('pathway_membership', 0):,})",
          f"- IGVF Portal protein evidence files: "
          f"**{s.get('igvf_evidence', 0):,}**", "",
          "## Per source",
          "| Source | Interactions |", "|---|---|"]
    for src, n in s.get("per_source", {}).items():
        md.append(f"| {src} | {n:,} |")
    md += ["", "## Per evidence type",
           "| Type | n |", "|---|---|"]
    for k, n in s.get("per_evidence_type", {}).items():
        md.append(f"| {k or '(unspecified)'} | {n:,} |")
    md += ["", "## Top-10 hub proteins (by degree)"]
    for pid, deg in s.get("top_hubs", [])[:10]:
        md.append(f"- `{pid}` — degree {deg}")
    (out / "report.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0


def cmd_kg_visualize(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("kg_visualize_" + (args.gene or "all"))
    out = kg_visualize(args.label or "viz", gene=args.gene,
                        ego_max_neighbors=args.max_neighbors)
    print(f"Output: {out}")
    return 0


def cmd_assay_survey(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("assay_survey")
    out = assay_survey(args.label or "survey", max_per_assay=args.max_per_assay)
    print(f"Output: {out}")
    return 0


def cmd_assay_figures(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("assay_figures")
    out = assay_figures(args.label or "assays")
    print(f"Output: {out}")
    return 0


def cmd_pipeline(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("pipeline_" + (args.label or "run"))
    sources = [s.strip().lower() for s in args.sources.split(",") if s.strip()]
    if "all" in sources:
        sources = ["biogrid", "intact", "huri", "reactome", "kegg", "igvf"]
    # Step 1 — download
    if not args.skip_download:
        for src in sources:
            try:
                if src == "biogrid": biogrid_download()
                elif src == "intact": intact_download()
                elif src == "huri": huri_download()
                elif src == "reactome": reactome_download()
                elif src == "kegg": kegg_download(max_pathways=args.kegg_max_pathways)
                elif src == "igvf": igvf_protein_download()
                elif src == "uniprot": uniprot_idmap_download()
            except Exception as e:
                logger.warning("download(%s) failed: %s", src, e)
    # Step 2 — build KG
    counts = build_kg(sources=sources, max_rows=args.max_rows)
    logger.info("Ingest counts: %s", counts)
    # Step 3 — stats
    cmd_kg_stats(argparse.Namespace(label=(args.label or "pipeline") + "_stats"))
    # Step 4 — visualize
    if args.gene:
        kg_visualize((args.label or "pipeline") + "_viz",
                      gene=args.gene, ego_max_neighbors=args.max_neighbors)
    else:
        kg_visualize((args.label or "pipeline") + "_viz",
                      gene=None, ego_max_neighbors=args.max_neighbors)
    # Step 5 — assay figures
    if not args.skip_assay_figures:
        try:
            assay_figures((args.label or "pipeline") + "_assays")
        except Exception as e:
            logger.warning("assay_figures: %s", e)
    # Step 6 — literature survey
    if not args.skip_literature:
        try:
            assay_survey((args.label or "pipeline") + "_lit",
                          max_per_assay=args.max_per_assay)
        except Exception as e:
            logger.warning("assay_survey: %s", e)
    print("Pipeline done.")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: Proteomics & PPI knowledge graph

End-to-end proteomics skill for IGVFagent. Aggregates **BioGRID**,
**IntAct**, **HuRI**, **Reactome**, **KEGG**, **UniProt id-mapping**,
and **IGVF Portal protein assays** into a local SQLite knowledge graph
that is incrementally update-aware. Provides summary statistics,
network visualizations (degree distribution, hub bar plot, ego graph),
literature surveys for VAMP-seq family / MAVE / semi-qY2H / DUAL-IPA via
the Reference skill, and per-assay example figures from real Portal
files.

## Subcommands

### download
```
igvfagent proteomics download --source all
igvfagent proteomics download --source biogrid
igvfagent proteomics download --source intact
igvfagent proteomics download --source huri
igvfagent proteomics download --source reactome
igvfagent proteomics download --source kegg --kegg-max-pathways 400
igvfagent proteomics download --source igvf
igvfagent proteomics download --source uniprot   # idmap, optional
```
The skill stores `Data/Proteomics/_versions.json` with version, URL, sha256,
and record count per source. BioGRID is auto-resolved against the latest
`Release-Archive` listing. IntAct uses `Last-Modified` on
`psimitab/intact-micluster.txt`. HuRI is mostly static.

### versions / update
```
igvfagent proteomics versions
igvfagent proteomics update
```
`update` only re-fetches sources where the upstream version differs from
the local one (delta refresh).

### igvf-protein
```
igvfagent proteomics igvf-protein
```
Pulls all 214 protein-slim MeasurementSets, the 14 PPI-score files, the
DUAL-IPA stability file, the UniProt protein references, and the protein
language model files from the IGVF Portal. Saves metadata JSON + actual
TSV/CSV file content under `Data/Proteomics/Sources/IGVF/`.

### build-kg
```
igvfagent proteomics build-kg --sources all
igvfagent proteomics build-kg --sources biogrid,reactome,igvf
```
Ingests downloaded sources into `Data/Proteomics/KG/proteomics.sqlite`.
Schema: `interactions`, `proteins`, `pathways`, `pathway_membership`,
`igvf_evidence`, `id_map`, `versions`. Edges deduped on
`(id_a, id_b, source, source_id)`.

### kg-stats
```
igvfagent proteomics kg-stats --label initial
```
Summary report under `Docs/Proteomics/<ts>_initial/`:
total interactions, distinct proteins, pathways, IGVF evidence,
breakdowns by source / evidence type / detection method, top-30 hubs.

### kg-visualize
```
igvfagent proteomics kg-visualize --gene TP53 --label tp53
igvfagent proteomics kg-visualize --label snapshot
```
Produces `Plots/degree_distribution.png`, `Plots/top_hubs.png`,
`Plots/per_source.png`, and (when `--gene` given) an ego graph
`Plots/ego_<GENE>.png` (requires networkx).

### assay-survey
```
igvfagent proteomics assay-survey --label may2026
```
Calls the Reference skill (PubMed + Semantic Scholar + OpenAlex) to
retrieve current studies on:
- VAMP-seq (MultiSTEP)
- VAMP-seq
- MAVE
- semi-qY2H v1 / v2 / v3 (semi-quantitative yeast two-hybrid)
- DUAL-IPA
Filters to Nature / Cell / Science family journals (Nat Genet, Nat Methods,
Nat Biotechnol, Mol Cell, Mol Syst Biol, Nat Commun) and writes
`literature_survey.md` and `.json`.

### assay-figures
```
igvfagent proteomics assay-figures --label demos
```
Generates per-assay example histograms from the actual IGVF Portal files
pulled by `igvf-protein`. One PNG per assay under
`Docs/Proteomics/<ts>_demos_assay_figures/Plots/`.

### pipeline
```
igvfagent proteomics pipeline --label may2026 --gene TP53 \\
  --sources biogrid,intact,huri,reactome,kegg,igvf
```
End-to-end: download → build-kg → kg-stats → kg-visualize → assay-figures
→ assay-survey. Use `--skip-download` to reuse cached files,
`--skip-literature` to skip the network-dependent literature scan,
`--max-rows N` to cap parser ingestion (useful for smoke tests).

### write-playbook
```
igvfagent proteomics write-playbook
```
Writes this document to `Docs/Skills/PROTEOMICS_SKILLS.md`.

## Storage

```
Data/Proteomics/
  _versions.json
  Sources/
    BioGRID/BIOGRID-ALL-X.Y.Z.tab3.zip
    IntAct/intact-micluster.txt
    HuRI/{HuRI.tsv, HI-union.tsv, Lit-BM.tsv}
    Reactome/{ReactomePathways.txt, UniProt2Reactome.txt,
              reactome.homo_sapiens.interactions.tab-delimited.txt}
    KEGG/{pathways_hsa.json, pathway_membership_hsa.json}
    UniProt/HUMAN_9606_idmapping_selected.tab.gz
    IGVF/{measurement_sets.json, ppi_files_meta.json,
          stability_files_meta.json, protein_reference_meta.json,
          protein_language_models.json, ppi_files/IGVFFI*}
  KG/proteomics.sqlite

Docs/Proteomics/<ts>_<label>/
  report.md
  Plots/*.png
  literature_survey.md (if assay-survey)
  per_assay_figures/Plots/*.png (if assay-figures)
```

## Cross-skill chaining

- `proteomics build-kg` → `kg gene <SYM>` to combine the proteomics PPI
  graph with the IGVF Catalog Knowledge Graph (variants, regulatory
  elements, pathways).
- `proteomics assay-survey` → `ref design --query "<assay>"` to draft
  a study design that mirrors the canonical analysis flow used by
  recent Nature/Cell/Science studies.
- `rnaseq deg` significant DEGs → query
  `proteomics kg-visualize --gene <SYM>` to see what proteins the
  upregulated genes interact with.
- `se-targets pipeline` → for each SE-target gene, look up upstream
  PPI partners with `proteomics kg-visualize` to suggest functional
  cofactors.

## Notes

- KEGG REST is throttled (default 0.4s/req) to respect their TOS;
  bulk distribution requires a license.
- IntAct full PSI-MITAB (`intact.zip`) is ~700MB — the skill defaults
  to the smaller `intact-micluster.txt`. Set `INTACT_BASE` env to
  override.
- HuRI uses Ensembl gene IDs. Run `proteomics download --source uniprot`
  to populate `id_map` for cross-resource harmonization.
"""


def cmd_write_playbook(_a) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_TEXT)
    print(f"Wrote: {PLAYBOOK_PATH}")
    return 0


# ---------------------------------------------------------------------------
# CLI wiring
# ---------------------------------------------------------------------------


def main(argv: Optional["list[str]"] = None) -> int:
    p = argparse.ArgumentParser(
        prog="proteomics",
        description="Proteomics & PPI: multi-source ingestion, KG, "
                    "analysis, viz, IGVF assay literature & figures.")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("download", help="Fetch latest from one source.")
    s.add_argument("--source", required=True,
                    help="biogrid|intact|huri|reactome|kegg|igvf|uniprot|all "
                         "(comma-separated)")
    s.add_argument("--biogrid-version", default=None,
                    help="Pin a specific BioGRID release (e.g. 4.4.244).")
    s.add_argument("--kegg-max-pathways", type=int, default=400)
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("versions",
                        help="Show local-vs-upstream versions.")
    s.set_defaults(func=cmd_versions)

    s = sub.add_parser("update",
                        help="Refresh only sources that differ upstream.")
    s.set_defaults(func=cmd_update)

    s = sub.add_parser("igvf-protein",
                        help="Pull IGVF Portal protein-assay manifest + files.")
    s.set_defaults(func=cmd_igvf_protein)

    s = sub.add_parser("build-kg",
                        help="Ingest downloaded sources into local SQLite KG.")
    s.add_argument("--sources", default="all")
    s.add_argument("--max-rows", type=int, default=0,
                    help="Cap rows per source (0 = unlimited; useful for smoke "
                         "tests).")
    s.set_defaults(func=cmd_build_kg)

    s = sub.add_parser("kg-stats", help="Summary stats on the integrated KG.")
    s.add_argument("--label", default=None)
    s.set_defaults(func=cmd_kg_stats)

    s = sub.add_parser("kg-visualize",
                        help="Degree/hub/ego plots for the integrated KG.")
    s.add_argument("--gene", default=None,
                    help="Symbol/UniProt to draw an ego graph around.")
    s.add_argument("--max-neighbors", type=int, default=60)
    s.add_argument("--label", default=None)
    s.set_defaults(func=cmd_kg_visualize)

    s = sub.add_parser("assay-survey",
                        help="Literature survey for IGVF protein assays "
                              "via the Reference skill.")
    s.add_argument("--label", default=None)
    s.add_argument("--max-per-assay", type=int, default=20)
    s.set_defaults(func=cmd_assay_survey)

    s = sub.add_parser("assay-figures",
                        help="Per-assay example figures from real IGVF Portal "
                              "files.")
    s.add_argument("--label", default=None)
    s.set_defaults(func=cmd_assay_figures)

    s = sub.add_parser("pipeline",
                        help="End-to-end: download → kg → stats → viz → "
                              "figures → literature.")
    s.add_argument("--sources", default="all")
    s.add_argument("--label", default="pipeline")
    s.add_argument("--gene", default=None)
    s.add_argument("--max-rows", type=int, default=0)
    s.add_argument("--max-neighbors", type=int, default=60)
    s.add_argument("--max-per-assay", type=int, default=20)
    s.add_argument("--kegg-max-pathways", type=int, default=400)
    s.add_argument("--skip-download", action="store_true")
    s.add_argument("--skip-literature", action="store_true")
    s.add_argument("--skip-assay-figures", action="store_true")
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/PROTEOMICS_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
