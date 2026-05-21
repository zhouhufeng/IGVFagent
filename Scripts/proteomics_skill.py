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
            # Skip rows where UniProt is missing (small molecules, complexes
            # that Reactome flags with `-` in the UniProt column).
            if a in ("-", "") or b in ("-", ""):
                continue
            a = a.split(":")[-1] if a else a
            b = b.split(":")[-1] if b else b
            if not a or not b:
                continue
            # Only accept canonical UniProt accessions (skip the numeric
            # Reactome internal pseudo-IDs that occasionally leak through).
            if not (a[0].isalpha() and b[0].isalpha()):
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
# VAMP-seq deep analysis (MaveDB scoresets + IGVF Portal raw inventory)
# ---------------------------------------------------------------------------

# Curated catalog of canonical published VAMP-seq scoresets on MaveDB.
# Each entry: target gene → MaveDB URN of the scoreset, PubMed ID, and
# domain track (start, end, label) for residue-level annotation in the
# per-position plots.
MAVEDB_VAMPSEQ_CATALOG = {
    "PTEN": {
        "urn":    "urn:mavedb:00000013-a-1",
        "paper":  "Matreyek et al. Nat Genet 2018  PMID 29785012",
        "uniprot": "P60484",
        "pdb_id":  "1d5r",                       # FowlerLab canonical
        "length": 403,
        "domains": [(1, 13, "PIP4-bind"),
                     (14, 185, "Phosphatase"),
                     (186, 353, "C2"),
                     (354, 403, "C-tail")],
    },
    "TPMT": {
        "urn":    "urn:mavedb:00000013-b-1",
        "paper":  "Matreyek et al. Nat Genet 2018  PMID 29785012",
        "uniprot": "P51580",
        "pdb_id":  "2bzg",                       # FowlerLab canonical
        "length": 245,
        "domains": [(1, 245, "Methyltransferase")],
    },
    "VKOR": {
        "urn":    "urn:mavedb:00000078-a-1",   # VKOR abundance, Suiter eLife 2020
        "paper":  "Suiter et al. eLife 2020  PMID 33198913",
        "uniprot": "Q9BQB6",
        "length": 163,
        "domains": [(1, 163, "VKOR fold")],
    },
    "PRKN": {
        "urn":    "urn:mavedb:00001173-a-1",   # Parkin VAMP-seq, Clausen Nat Commun 2024
        "paper":  "Clausen et al. Nat Commun 2024  PMID 38378758",
        "uniprot": "O60260",
        "length": 465,
        "domains": [(1, 76, "Ubl"),
                     (77, 144, "Linker"),
                     (145, 215, "RING0"),
                     (216, 327, "RING1"),
                     (328, 377, "IBR"),
                     (378, 465, "RING2")],
    },
    "CYP2C9": {
        "urn":    "urn:mavedb:00000095-a-1",
        "paper":  "Amorosi et al. Genome Med 2021  PMID 33648532",
        "uniprot": "P11712",
        "length": 490,
        "domains": [(1, 490, "P450")],
    },
    "NUDT15": {
        "urn":    "urn:mavedb:00000054-a-1",
        "paper":  "Suiter et al. PNAS 2020  PMID 32094184",
        "uniprot": "Q9NV35",
        "length": 164,
        "domains": [(1, 164, "Nudix")],
    },
}

# Categorical bin -> human label and color, matching the Matreyek 2018
# `abundance_class` column convention (0=lowest .. 4=hyperabundant).
ABUNDANCE_CLASS_LABELS = {
    0: ("low", "#762A83"),
    1: ("low-int", "#9970AB"),
    2: ("intermediate", "#C2A5CF"),
    3: ("WT-like", "#5AAE61"),
    4: ("hyper-abund", "#1B7837"),
}

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")    # 20 AAs + stop


def download_mavedb_scoreset(urn: str, dest_dir: Optional[Path] = None) -> Path:
    """Pull the score CSV for a MaveDB scoreset URN from the public REST API."""
    dest_dir = dest_dir or (SOURCES_DIR / "MaveDB")
    dest_dir.mkdir(parents=True, exist_ok=True)
    safe = urn.replace(":", "-").replace("/", "_")
    dest = dest_dir / f"{safe}.csv"
    url = f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores"
    if not dest.exists() or dest.stat().st_size < 100:
        logger.info("MaveDB: downloading %s -> %s", urn, dest.name)
        try:
            http_download(url, dest, timeout=180)
        except Exception as e:
            logger.error("MaveDB %s download failed: %s", urn, e)
            return dest
    return dest


_HGVS_PRO_RE = re.compile(
    r"^p\.\(?(?P<wt>[A-Z][a-z]{2}|=|Ter)(?P<pos>\d+)"
    r"(?P<alt>[A-Z][a-z]{2}|=|Ter|del|fs|\*)?\)?$"
)
_AA_3TO1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*",
}


def _parse_hgvs_pro(s: str) -> Optional[dict]:
    """Parse e.g. `p.Met1Val` -> {pos:1, wt:'M', alt:'V', kind:'missense'}."""
    if not s or s == "NA":
        return None
    m = _HGVS_PRO_RE.match(s.strip())
    if not m:
        return None
    wt3 = m.group("wt")
    alt3 = m.group("alt")
    pos = int(m.group("pos"))
    wt = _AA_3TO1.get(wt3, "?")
    if alt3 in (None, "=", ""):
        return {"pos": pos, "wt": wt, "alt": wt, "kind": "synonymous"}
    if alt3 == "Ter" or alt3 == "*":
        return {"pos": pos, "wt": wt, "alt": "*", "kind": "nonsense"}
    if alt3 in ("del", "fs"):
        return {"pos": pos, "wt": wt, "alt": "-", "kind": "indel"}
    alt = _AA_3TO1.get(alt3, "?")
    if alt == "?":
        return None
    if alt == wt:
        kind = "synonymous"
    else:
        kind = "missense"
    return {"pos": pos, "wt": wt, "alt": alt, "kind": kind}


_BIOPHYSICAL_COLS = (
    "rsa", "bfactor", "hydro1", "hydro2", "grantham",
    "AA1_psic", "AA2_psic", "delta_psic", "evolutionary_coupling_avg",
    "volume", "polarity", "tm", "ddg",
)


def _read_mavedb_csv(path: Path) -> "list[dict]":
    """Stream MaveDB scoreset CSV, parse hgvs_pro, attach numeric columns.

    Captures per-replicate scores, normal-approx 95 % CI bounds, and any
    biophysical feature columns (RSA, B-factor, hydrophobicity scales,
    Grantham distance, PSIC conservation, evolutionary coupling, Tm,
    Rosetta ΔΔG) that the upstream scoreset publishes. Used by
    ``analyze_vampseq_scoreset`` for the Fowler-lab-style
    feature-correlation panel.
    """
    rows = []

    def _num(x):
        if x in (None, "", "NA"):
            return None
        try:
            return float(x)
        except Exception:
            return None

    with path.open("r", encoding="utf-8", errors="replace") as fp:
        reader = csv.DictReader(fp)
        for r in reader:
            parsed = _parse_hgvs_pro(r.get("hgvs_pro", ""))
            if not parsed:
                continue
            score_f = _num(r.get("score"))
            if score_f is None:
                continue
            sd_f  = _num(r.get("sd"))
            se_f  = _num(r.get("se"))
            lo_f  = _num(r.get("lower_ci"))
            hi_f  = _num(r.get("upper_ci"))
            # Normal-approx fallback when explicit CI not provided
            if (lo_f is None or hi_f is None) and se_f is not None:
                lo_f = score_f - 1.96 * se_f
                hi_f = score_f + 1.96 * se_f
            try: cls = int(float(r.get("abundance_class", "0") or 0))
            except Exception: cls = -1
            reps = []
            for k in ("score1","score2","score3","score4",
                      "score5","score6","score7","score8"):
                v = _num(r.get(k))
                if v is not None:
                    reps.append(v)
            # Optional biophysical features (any subset may be present)
            biophys = {c: _num(r.get(c)) for c in _BIOPHYSICAL_COLS
                        if c in r and _num(r.get(c)) is not None}
            rows.append({**parsed, "score": score_f, "sd": sd_f, "se": se_f,
                          "lower_ci": lo_f, "upper_ci": hi_f,
                          "abundance_class": cls, "replicates": reps,
                          "biophys": biophys,
                          "expts": r.get("expts","")})
    return rows


def _matplotlib():
    import matplotlib  # type: ignore
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore
    return matplotlib, plt, LinearSegmentedColormap


def analyze_vampseq_scoreset(csv_path: Path, *, gene: str,
                               domains: "list[tuple]" = (),
                               length: Optional[int] = None,
                               label: Optional[str] = None,
                               pdb_id: Optional[str] = None) -> Path:
    """Deep VAMP-seq analysis on a MaveDB-format scoreset.

    Generates:
      1. distribution.png      — overlaid score densities (missense / syn / nonsense)
      2. heatmap.png           — residue × AA matrix (the iconic VAMP-seq plot)
      3. per_position.png      — per-residue mean ± IQR with domain track
      4. replicate_corr.png    — pairwise replicate scatter + Pearson r
      5. abundance_class.png   — categorical breakdown bar
      6. cumulative.png        — sorted-rank cumulative variant plot
      7. summary.json + report.md
    """
    _, plt, LSC = _matplotlib()
    import numpy as np  # used by several blocks below
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label or gene)}_vampseq"
    plots = out / "Plots"
    plots.mkdir(parents=True, exist_ok=True)

    rows = _read_mavedb_csv(csv_path)
    if not rows:
        raise RuntimeError(f"No usable rows parsed from {csv_path}")
    miss = [r["score"] for r in rows if r["kind"] == "missense"]
    syn  = [r["score"] for r in rows if r["kind"] == "synonymous"]
    non  = [r["score"] for r in rows if r["kind"] == "nonsense"]

    L = length or max(r["pos"] for r in rows)

    summary = {
        "gene": gene,
        "csv": str(csv_path),
        "n_rows": len(rows),
        "n_missense": len(miss),
        "n_synonymous": len(syn),
        "n_nonsense": len(non),
        "score_min": min(r["score"] for r in rows),
        "score_max": max(r["score"] for r in rows),
        "score_median": sorted(miss)[len(miss)//2] if miss else None,
        "length": L,
    }

    # --- 1. Distribution ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = 50
    if miss:
        ax.hist(miss, bins=bins, color="#1F77B4", alpha=0.65,
                  label=f"missense (n={len(miss)})", density=True)
    if syn:
        ax.hist(syn, bins=bins, color="#2CA02C", alpha=0.65,
                  label=f"synonymous (n={len(syn)})", density=True)
    if non:
        ax.hist(non, bins=bins, color="#D62728", alpha=0.65,
                  label=f"nonsense (n={len(non)})", density=True)
    ax.axvline(0.0, color="black", linestyle=":", lw=1, alpha=0.7)
    ax.axvline(1.0, color="black", linestyle=":", lw=1, alpha=0.7)
    ax.text(0.02, ax.get_ylim()[1]*0.92, "nonsense=0", fontsize=8)
    ax.text(1.02, ax.get_ylim()[1]*0.92, "WT=1", fontsize=8)
    ax.set_xlabel("Abundance score"); ax.set_ylabel("Density")
    ax.set_title(f"{gene} VAMP-seq score distribution")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(plots / "distribution.png", dpi=150)
    plt.close(fig)

    # --- 2. Heatmap (residue × AA) ------------------------------------------
    import math
    grid = [[float("nan")] * 21 for _ in range(L)]
    for r in rows:
        if r["kind"] not in ("missense", "synonymous", "nonsense"):
            continue
        if r["pos"] < 1 or r["pos"] > L:
            continue
        try:
            j = AA_ORDER.index(r["alt"])
        except ValueError:
            continue
        v = r["score"]
        if grid[r["pos"]-1][j] != grid[r["pos"]-1][j]:
            grid[r["pos"]-1][j] = v
        else:
            grid[r["pos"]-1][j] = (grid[r["pos"]-1][j] + v) / 2
    cmap = LSC.from_list("vamp", ["#67001F", "#D6604D", "#F7F7F7",
                                    "#92C5DE", "#053061"], N=256)
    fig_h = max(8, L / 25)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))
    arr = [[v if v == v else None for v in row] for row in grid]
    import numpy as np
    M = np.array([[v if v is not None else np.nan for v in row]
                   for row in arr], dtype=float)
    vmax = max(1.4, np.nanpercentile(M, 99))
    vmin = min(-0.2, np.nanpercentile(M, 1))
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                    interpolation="nearest", origin="upper")
    ax.set_xticks(range(21)); ax.set_xticklabels(AA_ORDER, fontsize=7)
    ystep = max(1, L // 25)
    ax.set_yticks(list(range(0, L, ystep)))
    ax.set_yticklabels([str(p+1) for p in range(0, L, ystep)], fontsize=7)
    ax.set_xlabel("Substituted amino acid")
    ax.set_ylabel(f"{gene} residue position (1–{L})")
    ax.set_title(f"{gene} VAMP-seq abundance heatmap (n={len(rows):,} variants)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Abundance score (0=null, 1=WT)")
    fig.tight_layout()
    fig.savefig(plots / "heatmap.png", dpi=160)
    plt.close(fig)

    # --- 3. Per-position mean +/- IQR with domain track ---------------------
    pos_means = [None] * L
    pos_q1 = [None] * L
    pos_q3 = [None] * L
    by_pos = defaultdict(list)
    for r in rows:
        if r["kind"] != "missense":
            continue
        if 1 <= r["pos"] <= L:
            by_pos[r["pos"]].append(r["score"])
    for p, vs in by_pos.items():
        vs2 = sorted(vs)
        n = len(vs2)
        pos_means[p-1] = sum(vs2)/n
        pos_q1[p-1] = vs2[n//4] if n >= 4 else min(vs2)
        pos_q3[p-1] = vs2[max(0, (3*n)//4 - 1)] if n >= 4 else max(vs2)
    fig, (ax_d, ax) = plt.subplots(2, 1, figsize=(11, 4),
                                     gridspec_kw={"height_ratios":[0.18,1]},
                                     sharex=True)
    ax_d.set_ylim(0, 1)
    palette = ["#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"]
    for i, (a, b, name) in enumerate(domains or [(1, L, "full length")]):
        ax_d.add_patch(plt.Rectangle((a-1, 0.2), b - a + 1, 0.6,
                                        color=palette[i % len(palette)],
                                        alpha=0.85))
        ax_d.text((a+b)/2, 0.5, name, ha="center", va="center",
                    fontsize=7, color="white", weight="bold")
    ax_d.set_yticks([]); ax_d.set_xlim(0, L)
    xs = [i+1 for i in range(L) if pos_means[i] is not None]
    ms = [pos_means[i] for i in range(L) if pos_means[i] is not None]
    q1 = [pos_q1[i] for i in range(L) if pos_means[i] is not None]
    q3 = [pos_q3[i] for i in range(L) if pos_means[i] is not None]
    ax.fill_between(xs, q1, q3, color="#3182BD", alpha=0.25, label="IQR")
    ax.plot(xs, ms, color="#08519C", lw=0.8, label="mean")
    ax.axhline(1.0, color="green", linestyle=":", lw=1, alpha=0.6)
    ax.axhline(0.0, color="red", linestyle=":", lw=1, alpha=0.6)
    ax.set_xlim(0, L); ax.set_ylim(min(-0.2, min(q1) if q1 else 0),
                                       max(1.4, max(q3) if q3 else 1.2))
    ax.set_xlabel(f"{gene} residue position")
    ax.set_ylabel("Missense abundance (per-residue)")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.25, linestyle=":")
    ax.set_title(f"{gene} VAMP-seq — per-residue abundance with domain track")
    fig.tight_layout()
    fig.savefig(plots / "per_position.png", dpi=160)
    plt.close(fig)

    # --- 4. Replicate concordance (pairwise score1 vs score2) ---------------
    pairs = []
    for r in rows:
        rps = r.get("replicates") or []
        if len(rps) >= 2:
            pairs.append((rps[0], rps[1]))
    if pairs:
        a = [p[0] for p in pairs]; b = [p[1] for p in pairs]
        ma = sum(a)/len(a); mb = sum(b)/len(b)
        cov = sum((x-ma)*(y-mb) for x, y in pairs)
        sa = sum((x-ma)**2 for x in a) ** 0.5
        sb = sum((y-mb)**2 for y in b) ** 0.5
        pearson = cov / (sa*sb) if sa*sb > 0 else 0.0
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(a, b, s=4, alpha=0.4, c="#08519C")
        lo = min(min(a), min(b)); hi = max(max(a), max(b))
        ax.plot([lo, hi], [lo, hi], "k:", lw=0.8)
        ax.set_xlabel("Replicate 1 score"); ax.set_ylabel("Replicate 2 score")
        ax.set_title(f"{gene} replicate concordance  Pearson r = {pearson:.3f}  (n={len(pairs):,})")
        ax.grid(alpha=0.25, linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "replicate_corr.png", dpi=150)
        plt.close(fig)
        summary["replicate_pearson"] = pearson
        summary["replicate_n_pairs"] = len(pairs)

    # --- 5. Abundance class breakdown ---------------------------------------
    cls_counts = Counter()
    for r in rows:
        if r["kind"] != "missense":
            continue
        cls_counts[r.get("abundance_class", -1)] += 1
    if cls_counts:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        keys = sorted(cls_counts.keys())
        labels = [ABUNDANCE_CLASS_LABELS.get(k, (str(k), "#888"))[0] for k in keys]
        colors = [ABUNDANCE_CLASS_LABELS.get(k, (str(k), "#888"))[1] for k in keys]
        vals = [cls_counts[k] for k in keys]
        ax.bar(range(len(keys)), vals, color=colors, edgecolor="white")
        ax.set_xticks(range(len(keys))); ax.set_xticklabels(labels, fontsize=9)
        for i, v in enumerate(vals):
            ax.text(i, v, f" {v:,}", ha="center", va="bottom", fontsize=8)
        ax.set_ylabel("Missense variants")
        ax.set_title(f"{gene} VAMP-seq abundance classes")
        ax.grid(alpha=0.25, axis="y", linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "abundance_class.png", dpi=150)
        plt.close(fig)
        summary["abundance_class_counts"] = {str(k): v for k, v in cls_counts.items()}

    # --- 6. Cumulative ranked variants (with CI band) -----------------------
    # FowlerLab repo emits per-variant CIs alongside the ranked-score plot.
    # Pair the cumulative curve with a shaded 95 % CI band when MaveDB or
    # the SE column let us compute it.
    miss_rows = [r for r in rows if r["kind"] == "missense"]
    miss_rows.sort(key=lambda r: r["score"])
    if miss_rows:
        fig, ax = plt.subplots(figsize=(7, 4))
        xs2 = list(range(len(miss_rows)))
        ys = np.array([r["score"] for r in miss_rows])
        ax.plot(xs2, ys, color="#08519C", lw=1, label="score")
        # Shaded 95 % CI band (sorted-rank order)
        if all(r.get("lower_ci") is not None and r.get("upper_ci") is not None
               for r in miss_rows[:50]):
            los = np.array([r["lower_ci"] or r["score"] for r in miss_rows])
            his = np.array([r["upper_ci"] or r["score"] for r in miss_rows])
            ax.fill_between(xs2, los, his, color="#08519C", alpha=0.18,
                              label="95 % CI")
        ax.axhline(0.5, color="#888", linestyle=":", lw=0.6)
        n_low = int((ys < 0.5).sum())
        ax.fill_between(xs2[:n_low], [0]*n_low, ys[:n_low],
                          color="#FCBBA1", alpha=0.6,
                          label=f"low (<0.5) n={n_low}")
        ax.set_xlabel("Variant rank"); ax.set_ylabel("Abundance score")
        ax.set_title(f"{gene} missense variants ranked by abundance")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.25, linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "cumulative.png", dpi=150)
        plt.close(fig)
        summary["frac_low_abundance"] = n_low / max(len(miss_rows), 1)

    # --- 7. Nonsense scores by position (Fowler-lab QC plot) ----------------
    # The canonical "do truncations crash the score?" check — scatter every
    # nonsense variant's score against its residue index, with WT=1 and
    # nonsense=0 anchors. Truncations clustered near 0 = healthy library.
    non_rows = [r for r in rows if r["kind"] == "nonsense"]
    if non_rows:
        fig, ax = plt.subplots(figsize=(9, 3.2))
        ax.scatter([r["pos"] for r in non_rows],
                     [r["score"] for r in non_rows],
                     s=10, alpha=0.7, c="#D62728", edgecolor="white",
                     linewidth=0.4)
        ax.axhline(1.0, color="green", linestyle=":", lw=0.8, label="WT = 1")
        ax.axhline(0.0, color="red", linestyle=":", lw=0.8, label="nonsense = 0")
        ax.set_xlim(0, L)
        ax.set_xlabel(f"{gene} residue position")
        ax.set_ylabel("Nonsense abundance")
        ax.set_title(f"{gene} truncation (nonsense) score by position "
                       f"— n={len(non_rows)} stop codons")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.25, linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "nonsense_by_position.png", dpi=150)
        plt.close(fig)
        summary["nonsense_median"] = float(np.median([r["score"]
                                                          for r in non_rows]))

    # --- 8. Per-position moving-average median (Fowler-lab style) -----------
    # FowlerLab uses ma(score, n=3); we render a small panel alongside the
    # mean ± IQR view rather than replacing it.
    pos_median = [None] * L
    by_pos2 = defaultdict(list)
    for r in rows:
        if r["kind"] == "missense" and 1 <= r["pos"] <= L:
            by_pos2[r["pos"]].append(r["score"])
    for p, vs in by_pos2.items():
        s2 = sorted(vs)
        pos_median[p-1] = s2[len(s2)//2]
    xs3 = [i+1 for i in range(L) if pos_median[i] is not None]
    ms_med = [pos_median[i] for i in range(L) if pos_median[i] is not None]
    if ms_med:
        # 3-residue centred moving average over the median
        win = 3
        ma_y = []
        for k in range(len(ms_med)):
            lo = max(0, k - win // 2)
            hi = min(len(ms_med), k + win // 2 + 1)
            ma_y.append(np.mean(ms_med[lo:hi]))
        fig, ax = plt.subplots(figsize=(11, 3))
        ax.plot(xs3, ms_med, color="#9ECAE1", lw=0.6, alpha=0.8,
                 label="per-residue median")
        ax.plot(xs3, ma_y, color="#08519C", lw=1.1,
                 label="3-residue moving average")
        ax.axhline(1.0, color="green", ls=":", lw=0.5, alpha=0.6)
        ax.axhline(0.0, color="red", ls=":", lw=0.5, alpha=0.6)
        ax.set_xlim(0, L)
        ax.set_xlabel(f"{gene} residue position")
        ax.set_ylabel("Per-residue missense median")
        ax.set_title(f"{gene} — per-residue median (Fowler-lab style)")
        ax.legend(loc="lower left", fontsize=8)
        ax.grid(alpha=0.25, linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "per_position_median_ma.png", dpi=150)
        plt.close(fig)

    # --- 9. N×N replicate scatter matrix (when ≥3 reps) ---------------------
    # Generalises the rep-1-vs-rep-2 plot to the full pairwise grid for
    # scoresets that publish 3+ replicate columns. PTEN has 8.
    rep_arrays: "list[list[float]]" = []
    max_reps = max((len(r.get("replicates") or []) for r in rows), default=0)
    if max_reps >= 3:
        for k in range(max_reps):
            rep_arrays.append([
                r["replicates"][k]
                for r in rows
                if len(r.get("replicates") or []) > k
                and r["replicates"][k] is not None
            ])
        # Pad: align on the row index so pairwise scatter shares cells
        aligned = []
        for r in rows:
            reps = r.get("replicates") or []
            if len(reps) >= max_reps:
                aligned.append(reps[:max_reps])
        if aligned and max_reps <= 8:
            A = np.array(aligned, dtype=float)
            n_panels = max_reps
            fig, axes = plt.subplots(n_panels, n_panels,
                                       figsize=(2 * n_panels, 2 * n_panels))
            corrs = np.eye(n_panels)
            for i in range(n_panels):
                for j in range(n_panels):
                    ax = axes[i][j]
                    if i == j:
                        ax.hist(A[:, i], bins=40, color="#08519C", alpha=0.85)
                        ax.set_xticks([]); ax.set_yticks([])
                        ax.set_title(f"rep{i+1}", fontsize=8)
                    elif i > j:
                        ax.scatter(A[:, j], A[:, i], s=2, alpha=0.3,
                                     c="#3182BD")
                        ax.set_xticks([]); ax.set_yticks([])
                    else:
                        ma_ = ~np.isnan(A[:, i]) & ~np.isnan(A[:, j])
                        if ma_.sum() >= 3:
                            xa, ya = A[ma_, j], A[ma_, i]
                            mx, my = xa.mean(), ya.mean()
                            cov = ((xa - mx) * (ya - my)).sum()
                            sx = ((xa - mx) ** 2).sum() ** 0.5
                            sy = ((ya - my) ** 2).sum() ** 0.5
                            r = cov / (sx * sy) if sx * sy > 0 else 0.0
                            corrs[i, j] = corrs[j, i] = r
                            ax.text(0.5, 0.5, f"r = {r:.2f}",
                                     ha="center", va="center",
                                     transform=ax.transAxes,
                                     fontsize=10, color="#08519C",
                                     fontweight="bold")
                        ax.set_xticks([]); ax.set_yticks([])
            fig.suptitle(f"{gene} VAMP-seq replicate matrix (n_reps={max_reps})",
                           fontsize=11)
            fig.tight_layout()
            fig.savefig(plots / "replicate_matrix.png", dpi=150)
            plt.close(fig)
            summary["replicate_matrix_min_r"] = float(corrs[np.triu_indices(
                n_panels, k=1)].min())
            summary["replicate_matrix_max_r"] = float(corrs[np.triu_indices(
                n_panels, k=1)].max())

    # --- 10. Biophysical-feature Spearman ρ panel ---------------------------
    # When the scoreset publishes biophysical columns (RSA, B-factor,
    # hydrophobicity, Grantham, PSIC conservation, ΔΔG, Tm, …), correlate
    # each one with the abundance score — Fowler-lab Fig 3 / 4 staple.
    feature_corr: "dict[str, float]" = {}
    feature_n:    "dict[str, int]" = {}
    for col in _BIOPHYSICAL_COLS:
        pairs = [(r["biophys"][col], r["score"])
                 for r in rows
                 if r.get("biophys", {}).get(col) is not None
                 and r["kind"] == "missense"]
        if len(pairs) < 10:
            continue
        xs_, ys_ = np.array([p[0] for p in pairs]), \
                    np.array([p[1] for p in pairs])
        # Spearman ρ via rankdata to avoid scipy dependency
        def _rank(a):
            order = a.argsort()
            ranks = np.empty(len(a))
            ranks[order] = np.arange(len(a))
            return ranks
        rx, ry = _rank(xs_), _rank(ys_)
        mx, my = rx.mean(), ry.mean()
        cov = ((rx - mx) * (ry - my)).sum()
        sx, sy = ((rx - mx) ** 2).sum() ** 0.5, ((ry - my) ** 2).sum() ** 0.5
        rho = cov / (sx * sy) if sx * sy > 0 else 0.0
        feature_corr[col] = float(rho)
        feature_n[col] = len(pairs)
    if feature_corr:
        fig, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(feature_corr))))
        keys = list(feature_corr.keys())
        vals = [feature_corr[k] for k in keys]
        bar_colors = ["#D62728" if v < 0 else "#1F77B4" for v in vals]
        ax.barh(range(len(keys)), vals, color=bar_colors)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels([f"{k} (n={feature_n[k]})" for k in keys],
                             fontsize=8)
        ax.axvline(0, color="black", lw=0.4)
        ax.set_xlabel("Spearman ρ vs abundance score")
        ax.set_title(f"{gene} — biophysical feature correlation")
        ax.set_xlim(-1, 1)
        ax.grid(alpha=0.25, axis="x", linestyle=":")
        fig.tight_layout()
        fig.savefig(plots / "feature_correlations.png", dpi=150)
        plt.close(fig)
        summary["biophysical_spearman"] = feature_corr

    # --- 11. PyMOL colorscale .pml (opt-in via pdb_id) ----------------------
    # FowlerLab emits PDB-mapped abundance overlays; here we generate a
    # ready-to-source .pml that loads the user-specified structure and
    # paints the canonical blue → white → red scale onto chain B-factor
    # set from each residue's median missense abundance score.
    if pdb_id:
        pml_lines = [
            f"# IGVFagent VAMP-seq abundance overlay for {gene}",
            f"# Source: {csv_path.name}",
            f"load https://files.rcsb.org/download/{pdb_id}.cif, struct",
            "hide everything",
            "show cartoon",
            "color grey80",
        ]
        # Reset b-factors using a dict of {resi: score_median}
        pml_lines.append("alter struct, b = 0")
        for p, vs in by_pos2.items():
            med = float(np.median(vs))
            pml_lines.append(
                f"alter struct and resi {p}, b = {med:.4f}")
        pml_lines += [
            "rebuild",
            "spectrum b, blue_white_red, struct, minimum=0.0, maximum=1.4",
            "bg_color white",
            "ray 1200, 900",
            f"png Plots/{safe_label(gene)}_structure.png, dpi=150",
        ]
        (plots / f"{safe_label(gene)}_structure.pml").write_text(
            "\n".join(pml_lines) + "\n")
        summary["pymol_pml"] = str(plots / f"{safe_label(gene)}_structure.pml")
        summary["pdb_id"] = pdb_id

    # --- Report -------------------------------------------------------------
    md = ["# VAMP-seq deep analysis — " + gene, "",
          f"- Source: `{csv_path.name}`  ({summary['n_rows']:,} variants)",
          f"- Length parsed: **{L}** residues",
          f"- Missense / Synonymous / Nonsense: "
          f"**{len(miss):,} / {len(syn):,} / {len(non):,}**",
          f"- Score range: **{summary['score_min']:.2f} … {summary['score_max']:.2f}**"
          f" (anchors: nonsense=0, WT=1)",
          ]
    if "replicate_pearson" in summary:
        md.append(f"- Replicate Pearson r: **{summary['replicate_pearson']:.3f}** "
                  f"(n={summary['replicate_n_pairs']:,} variant pairs)")
    if "frac_low_abundance" in summary:
        md.append(f"- Fraction of missense with low abundance "
                  f"(score < 0.5): **{summary['frac_low_abundance']*100:.1f}%**")
    md += ["",
           "## Plots",
           "- `Plots/distribution.png` — score densities by variant kind",
           "- `Plots/heatmap.png` — residue × AA abundance matrix (the iconic VAMP-seq view)",
           "- `Plots/per_position.png` — per-residue mean ± IQR with domain track",
           "- `Plots/per_position_median_ma.png` — per-residue median + 3-residue moving average",
           "- `Plots/replicate_corr.png` — replicate-1 vs replicate-2 concordance",
           "- `Plots/abundance_class.png` — categorical class breakdown",
           "- `Plots/cumulative.png` — variants sorted by abundance rank (with 95 % CI band)",
           "- `Plots/nonsense_by_position.png` — truncation-score QC (Fowler-lab style)"]
    if "replicate_matrix_min_r" in summary:
        md.append("- `Plots/replicate_matrix.png` — N×N replicate scatter "
                  f"(r ∈ [{summary['replicate_matrix_min_r']:.2f}, "
                  f"{summary['replicate_matrix_max_r']:.2f}])")
    if "biophysical_spearman" in summary:
        md.append("- `Plots/feature_correlations.png` — biophysical-feature Spearman ρ")
    if pdb_id and summary.get("pymol_pml"):
        md.append(f"- `{summary['pymol_pml']}` — PyMOL `.pml` overlay "
                  f"for PDB **{pdb_id}**  "
                  f"(open in PyMOL and source to render)")
    if cls_counts:
        md += ["", "## Abundance class counts (missense)",
               "| Class | Count |", "|---|---|"]
        for k in sorted(cls_counts.keys()):
            lab = ABUNDANCE_CLASS_LABELS.get(k, (str(k), ""))[0]
            md.append(f"| {k} ({lab}) | {cls_counts[k]:,} |")
    (out / "report.md").write_text("\n".join(md) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2,
                                                   default=str))
    return out


# --- IGVF Portal raw VAMP-seq inventory analysis ---------------------------

_ALIAS_RE = re.compile(
    r"(?P<gene>[A-Z0-9]+)-DMS-(?P<antibody>[A-Za-z0-9\-]+?)"
    r"-Tile(?P<tile>\d+)-Replicate(?P<rep>\d+)-Bin(?P<bin>\d+)",
    re.IGNORECASE,
)


def _decode_igvf_alias(alias: str) -> Optional[dict]:
    """Decode aliases like
       'lea-starita:F9-DMS-Light-chain-antibody-Tile1-Replicate1-Bin1-FileSet'."""
    if not alias:
        return None
    bare = alias.split(":", 1)[-1]
    m = _ALIAS_RE.search(bare)
    if not m:
        return None
    return {"gene": m.group("gene").upper(),
            "antibody": m.group("antibody").rstrip("-"),
            "tile": int(m.group("tile")),
            "replicate": int(m.group("rep")),
            "bin": int(m.group("bin"))}


def inventory_igvf_vampseq(label: str) -> Path:
    """Decode the 144 VAMP-seq (MultiSTEP) MeasurementSets into a tile × bin
    × replicate × antibody coverage matrix and render heatmaps."""
    _, plt, _ = _matplotlib()
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}_igvf_vampseq"
    plots = out / "Plots"
    plots.mkdir(parents=True, exist_ok=True)

    ms_path = SOURCES_DIR / "IGVF" / "measurement_sets.json"
    if not ms_path.exists():
        # Pull on the fly
        igvf_protein_download()
    ms = json.loads(ms_path.read_text())
    rows = []
    for m in ms:
        titles = m.get("preferred_assay_titles") or []
        if "VAMP-seq (MultiSTEP)" not in titles and "VAMP-seq" not in titles:
            continue
        for a in (m.get("aliases") or []):
            d = _decode_igvf_alias(a)
            if d:
                rows.append({**d, "accession": m.get("accession"),
                              "lab": (m.get("lab") or {}).get("title", "")
                                if isinstance(m.get("lab"), dict)
                                else m.get("lab", ""),
                              "summary": (m.get("summary") or "")[:200],
                              "preferred_assay_titles": ", ".join(titles)})
                break

    inv_path = out / "inventory.json"
    inv_path.write_text(json.dumps(rows, indent=2))

    # Per-gene tile × bin coverage
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)

    md = ["# IGVF Portal raw VAMP-seq inventory", "",
          f"- {len(rows):,} MeasurementSets decoded across **{len(by_gene)} target genes**",
          ""]
    for gene, recs in sorted(by_gene.items(), key=lambda x: -len(x[1])):
        tiles = sorted({r["tile"] for r in recs})
        bins = sorted({r["bin"] for r in recs})
        reps = sorted({r["replicate"] for r in recs})
        antibodies = sorted({r["antibody"] for r in recs})
        md += [f"## {gene}",
                f"- {len(recs)} MeasurementSets",
                f"- Tiles: {tiles}",
                f"- Bins: {bins}",
                f"- Replicates: {reps}",
                f"- Antibody readouts: {antibodies}",
                ""]
        # Tile × bin presence-count matrix (summed over rep × antibody)
        mat = [[0] * len(bins) for _ in tiles]
        for r in recs:
            ti = tiles.index(r["tile"]); bi = bins.index(r["bin"])
            mat[ti][bi] += 1
        import numpy as np
        M = np.array(mat)
        fig, ax = plt.subplots(figsize=(max(4, len(bins)*0.7),
                                          max(2.4, len(tiles)*0.6)))
        im = ax.imshow(M, aspect="auto", cmap="YlGnBu")
        for i in range(len(tiles)):
            for j in range(len(bins)):
                ax.text(j, i, str(mat[i][j]), ha="center", va="center",
                          fontsize=8,
                          color="white" if mat[i][j] > M.max()/2 else "black")
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels([f"Bin{b}" for b in bins], fontsize=8)
        ax.set_yticks(range(len(tiles)))
        ax.set_yticklabels([f"Tile{t}" for t in tiles], fontsize=8)
        ax.set_title(f"{gene}: VAMP-seq MeasurementSet count by tile × bin")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        fig.tight_layout()
        fig.savefig(plots / f"{safe_label(gene)}_tile_bin_matrix.png", dpi=150)
        plt.close(fig)

        # Antibody × tile heatmap (replicate completeness)
        if len(antibodies) > 1:
            ant_mat = [[0] * len(tiles) for _ in antibodies]
            for r in recs:
                ai = antibodies.index(r["antibody"]); ti = tiles.index(r["tile"])
                ant_mat[ai][ti] += 1
            A = np.array(ant_mat)
            fig, ax = plt.subplots(figsize=(max(4, len(tiles)*0.8),
                                              max(2.4, len(antibodies)*0.6)))
            im = ax.imshow(A, aspect="auto", cmap="Purples")
            for i in range(len(antibodies)):
                for j in range(len(tiles)):
                    ax.text(j, i, str(ant_mat[i][j]), ha="center",
                              va="center", fontsize=8,
                              color="white" if ant_mat[i][j] > A.max()/2
                              else "black")
            ax.set_xticks(range(len(tiles)))
            ax.set_xticklabels([f"Tile{t}" for t in tiles], fontsize=8)
            ax.set_yticks(range(len(antibodies)))
            ax.set_yticklabels(antibodies, fontsize=8)
            ax.set_title(f"{gene}: VAMP-seq antibody × tile coverage")
            fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
            fig.tight_layout()
            fig.savefig(plots / f"{safe_label(gene)}_antibody_tile_matrix.png",
                          dpi=150)
            plt.close(fig)

    md += ["",
           "## Decoded fields per MeasurementSet",
           "Aliases follow the pattern "
           "`<lab>:<GENE>-DMS-<antibody>-Tile<i>-Replicate<j>-Bin<k>-FileSet`",
           "and are decoded into per-gene tile / bin / replicate / antibody slots.",
           "",
           f"Full inventory: `inventory.json` ({len(rows):,} rows)"]
    (out / "report.md").write_text("\n".join(md) + "\n")
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


def cmd_vampseq_pull(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("vampseq_pull_" + (args.gene or "all"))
    targets = ([args.gene.upper()] if args.gene
                else list(MAVEDB_VAMPSEQ_CATALOG.keys()))
    pulled = []
    for g in targets:
        meta = MAVEDB_VAMPSEQ_CATALOG.get(g)
        if not meta:
            logger.warning("No MaveDB URN curated for %s", g)
            continue
        p = download_mavedb_scoreset(meta["urn"])
        if p.exists() and p.stat().st_size > 0:
            pulled.append({"gene": g, "csv": str(p),
                            "size_bytes": p.stat().st_size,
                            "urn": meta["urn"], "paper": meta["paper"]})
    print(json.dumps(pulled, indent=2))
    return 0


def cmd_vampseq_analyze(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("vampseq_analyze_" + (args.gene or "all"))
    targets = ([args.gene.upper()] if args.gene
                else list(MAVEDB_VAMPSEQ_CATALOG.keys()))
    outs = []
    for g in targets:
        meta = MAVEDB_VAMPSEQ_CATALOG.get(g)
        if not meta:
            logger.warning("No catalog entry for %s", g)
            continue
        csv = download_mavedb_scoreset(meta["urn"])
        if not csv.exists() or csv.stat().st_size < 100:
            logger.warning("CSV for %s missing or empty", g)
            continue
        try:
            out = analyze_vampseq_scoreset(
                csv, gene=g,
                domains=meta.get("domains") or (),
                length=meta.get("length"),
                label=args.label or g,
                pdb_id=args.pdb_id or
                  (meta.get("pdb_id") if isinstance(meta, dict) else None))
            outs.append(str(out))
        except Exception as e:
            logger.error("Analyze failed for %s: %s", g, e)
    print(json.dumps(outs, indent=2))
    return 0



def cmd_vampseq_showcase(args) -> int:
    """Single-command showcase: download MaveDB scoreset + full Fowler-lab
    plot suite + composite publication figure + deep narrative report."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _composite_figure import build_composite_figure
    mkdirs()
    gene = (getattr(args, "gene", None) or "PTEN").upper()
    setup_logging("vampseq_showcase_" + gene)
    meta = MAVEDB_VAMPSEQ_CATALOG.get(gene)
    if not meta:
        raise SystemExit(
            f"No MaveDB VAMP-seq catalog entry for {gene}. "
            f"Available: {sorted(MAVEDB_VAMPSEQ_CATALOG)}"
        )
    csv_path = download_mavedb_scoreset(meta["urn"])
    if not csv_path.exists() or csv_path.stat().st_size < 100:
        raise SystemExit(f"Could not download {gene} scoreset from MaveDB.")
    logger.info("Running deep analysis for %s (urn=%s) ...", gene, meta["urn"])
    out_dir = analyze_vampseq_scoreset(
        csv_path, gene=gene,
        domains=meta.get("domains") or (),
        length=meta.get("length"),
        label=getattr(args, "label", None) or (gene + "_showcase"),
        pdb_id=getattr(args, "pdb_id", None) or meta.get("pdb_id"),
    )
    plots_dir = out_dir / "Plots"

    composite_layout = [
        ("Score distribution",           "distribution.png",           (0, 0)),
        ("Residue x amino-acid heatmap", "heatmap.png",                (0, 1)),
        ("Per-residue mean +/- IQR",     "per_position.png",           (0, 2)),
        ("Per-residue median + 3-aa MA", "per_position_median_ma.png", (1, 0)),
        ("Replicate concordance",        "replicate_corr.png",         (1, 1)),
        ("N x N replicate matrix",       "replicate_matrix.png",       (1, 2)),
        ("Abundance class summary",      "abundance_class.png",        (2, 0)),
        ("Nonsense-by-position QC",      "nonsense_by_position.png",   (2, 1)),
        ("Cumulative ranked variants",   "cumulative.png",             (2, 2)),
    ]
    composite = build_composite_figure(
        plots_dir,
        out_path=plots_dir / "composite_publication_figure.png",
        title=f"{gene} VAMP-seq -- comprehensive abundance landscape",
        subtitle=(f"Source: MaveDB {meta['urn']}  -  {meta.get('paper', '')}  -  "
                   f"length={meta.get('length', '?')} aa"),
        layout=composite_layout, n_rows=3, n_cols=3,
        panel_w=5.5, panel_h=4.5,
    )

    report = out_dir / "showcase_report.md"
    lines = [
        f"# {gene} VAMP-seq -- comprehensive showcase report",
        "",
        f"Generated: {timestamp()}",
        f"Source: MaveDB `{meta['urn']}`  -  paper: {meta.get('paper', 'unknown')}",
        f"Protein length: **{meta.get('length', '?')} aa**  -  UniProt: `{meta.get('uniprot', '?')}`",
        "",
        "## Figure suite",
        "",
        ("![composite](Plots/composite_publication_figure.png)"
         if composite else "_(composite figure not generated)_"),
        "",
        "### Individual panels",
        "",
        "- `Plots/distribution.png` -- overlaid missense / synonymous / nonsense score densities",
        "- `Plots/heatmap.png` -- residue x amino-acid abundance matrix",
        "- `Plots/per_position.png` -- per-residue mean abundance with IQR + domain track",
        "- `Plots/per_position_median_ma.png` -- per-residue median + 3-aa moving average",
        "- `Plots/replicate_corr.png` -- rep-1 vs rep-2 concordance scatter",
        "- `Plots/replicate_matrix.png` -- N x N pairwise replicate matrix",
        "- `Plots/abundance_class.png` -- abundant / hypomorph / null classification bars",
        "- `Plots/nonsense_by_position.png` -- nonsense-control QC vs position",
        "- `Plots/cumulative.png` -- cumulative-ranked variant curve with 95% CI",
    ]
    if (plots_dir / "feature_correlations.png").is_file():
        lines.append("- `Plots/feature_correlations.png` -- Spearman rho vs RSA, B-factor, hydrophobicity, Grantham, PSIC")
    lines += [
        "",
        "## How to read the figure suite",
        "",
        "1. **Distribution plot:** sanity check the assay -- synonymous variants cluster at ~1.0, nonsense at ~0, missense spans both.",
        "2. **Residue x AA heatmap:** vertical columns of low abundance across all AAs at a single residue mark structurally critical positions.",
        "3. **Per-residue mean/median + domain track:** runs of buried (intolerant) residues correlate with secondary-structure elements.",
        "4. **Replicate concordance:** Pearson r > 0.85 between any two replicates is the typical QC gate.",
        "5. **Nonsense-by-position:** confirms the NMD calibration; expect uniformly low abundance except near the C-terminus.",
        "6. **Biophysical correlations (if present):** RSA usually wins for stability-driven loss of function.",
        "",
        f"All artefacts under: `{out_dir}`",
    ]
    report.write_text("\n".join(lines))
    print(f"Report: {report}")
    print(f"Output dir: {out_dir}")
    if composite:
        print(f"Composite figure: {composite}")
        print(f"Composite (SVG): {composite.with_suffix('.svg')}")
    for ppng in sorted(plots_dir.glob("*.png")):
        print(f"Wrote plot: {ppng}")
    return 0


def cmd_vampseq_inventory(args: argparse.Namespace) -> int:
    mkdirs()
    setup_logging("vampseq_inventory")
    out = inventory_igvf_vampseq(args.label or "igvf_vampseq")
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

### vampseq-pull / vampseq-analyze / vampseq-inventory
```
# 1) Pull canonical published scoresets from MaveDB (PTEN, TPMT, VKOR,
#    PRKN, CYP2C9, NUDT15) — full per-replicate score CSVs
igvfagent proteomics vampseq-pull
igvfagent proteomics vampseq-pull --gene PTEN

# 2) Deep analysis — produces the Fowler-lab-style plot suite per gene
igvfagent proteomics vampseq-analyze --gene PTEN --label pten_deep \
    --pdb-id 1d5r                          # optional PyMOL .pml export
igvfagent proteomics vampseq-analyze       # all 6 catalogued targets

# 3) Inventory the IGVF Portal raw VAMP-seq experiments by decoding the
#    alias scheme (<lab>:<GENE>-DMS-<antibody>-Tile<i>-Replicate<j>-Bin<k>).
#    Produces tile×bin and antibody×tile coverage matrices per gene.
igvfagent proteomics vampseq-inventory --label igvf_f9
```
The deep analysis follows the canonical VAMP-seq pipeline distilled from
Matreyek et al. *Nat Genet* 2018, Suiter *eLife* 2020, Clausen *Nat Commun*
2024, and Coyote-Maestas *Nat Commun* 2024 (MultiSTEP), augmented to
match the public **FowlerLab/VAMPseq** analysis Rmd plot suite:

  1. **Distribution** — overlay missense / synonymous / nonsense densities
     anchored at WT=1, nonsense=0.
  2. **Residue × AA heatmap** — the iconic VAMP-seq view; one cell per
     (position, substituted AA), color = abundance score.
  3. **Per-position mean ± IQR with domain track** — annotates which
     domain each residue belongs to (e.g. PTEN phosphatase / C2 / C-tail,
     PRKN Ubl / RING0 / RING1 / IBR / RING2).
  4. **Per-position median + 3-residue moving average** — Fowler-lab
     `ma(score, n=3)` style smoothing alongside the per-residue median.
  5. **Replicate concordance** — Pearson r between rep-1 and rep-2 per
     variant.
  6. **N × N replicate scatter matrix** — emitted when ≥ 3 replicates
     are present (e.g. PTEN with 8 reps). Lower triangle = scatter,
     diagonal = per-rep histogram, upper triangle = Pearson r.
  7. **Abundance class breakdown** — categorical bar (low / low-int /
     intermediate / WT-like / hyper-abundant), matching the Matreyek
     2018 `abundance_class` convention.
  8. **Cumulative ranked variants with 95 % CI band** — sorted score
     curve, CI envelope from `lower_ci` / `upper_ci` (or normal-approx
     from `se`), low-abundance fraction (score < 0.5) shaded.
  9. **Nonsense-by-position scatter** — the canonical "do truncations
     crash the score?" QC plot from Matreyek 2018 Fig 1.
  10. **Biophysical-feature Spearman ρ panel** — emitted only when the
      scoreset carries RSA / B-factor / hydrophobicity / Grantham /
      PSIC conservation / ΔΔG / Tm columns (most public MaveDB
      scoresets don't; the FowlerLab supplementary PTEN / TPMT tables
      do).
  11. **PyMOL .pml export** — opt-in via `--pdb-id`. Writes a ready-to-
      source `.pml` that loads the structure, sets per-residue
      B-factor to the median missense abundance score, and applies the
      blue → white → red colorscale. Default PDB IDs are pinned to
      the canonical Fowler-lab structures (PTEN → `1d5r`,
      TPMT → `2bzg`).

Catalogued MaveDB targets (URN, paper, length, domains) are in
`MAVEDB_VAMPSEQ_CATALOG` in `proteomics_skill.py` — extend this dict to
analyze additional published scoresets.

The `vampseq-inventory` command decodes the IGVF Portal `aliases` field
to build a coverage matrix across the 144 MultiSTEP MeasurementSets
(currently all targeting **F9 / Coagulation Factor IX** across 3 tiles ×
4 bins × 4 replicates × 5 antibody readouts) and the 36 plain VAMP-seq
sets (CYP2C19, G6PD).

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

    s = sub.add_parser("vampseq-pull",
                        help="Download published VAMP-seq scoresets from MaveDB "
                              "(PTEN, TPMT, VKOR, PRKN, CYP2C9, NUDT15).")
    s.add_argument("--gene", default=None,
                    help="Gene symbol to pull. Default = all curated entries.")
    s.set_defaults(func=cmd_vampseq_pull)

    s = sub.add_parser("vampseq-analyze",
                        help="Deep VAMP-seq analysis: distribution, residue×AA "
                              "heatmap, per-residue mean/median (+ moving avg), "
                              "replicate concordance + N×N matrix, abundance "
                              "class, ranked variants with 95 % CI, "
                              "nonsense-by-position QC, biophysical-feature "
                              "Spearman ρ, and an optional PyMOL .pml export.")
    s.add_argument("--gene", default=None,
                    help="Gene to analyze (PTEN/TPMT/VKOR/PRKN/CYP2C9/NUDT15) "
                         "or omit to analyze all available.")
    s.add_argument("--label", default=None)
    s.add_argument("--pdb-id", default=None,
                    help="Optional PDB id (e.g. 1d5r for PTEN). When provided, "
                         "writes a PyMOL .pml that loads the structure and "
                         "applies the blue→white→red abundance colorscale.")
    s.set_defaults(func=cmd_vampseq_analyze)

    s = sub.add_parser("vampseq-showcase",
                        help="COMPREHENSIVE one-command VAMP-seq demo. "
                              "Pulls a MaveDB scoreset, runs the full "
                              "Fowler-lab plot suite (10 plots), builds a "
                              "9-panel publication composite figure, and "
                              "writes a deep narrative report. THIS IS THE "
                              "RIGHT TOOL FOR ANY VAMP-SEQ DEMO QUESTION.")
    s.add_argument("--gene", default="PTEN",
                    help="Gene to showcase (default PTEN; also TPMT, VKOR, "
                          "PRKN, CYP2C9, NUDT15).")
    s.add_argument("--label", default=None)
    s.add_argument("--pdb-id", default=None,
                    help="Override PDB id used for PyMOL overlay.")
    s.set_defaults(func=cmd_vampseq_showcase)

    s = sub.add_parser("vampseq-inventory",
                        help="Inventory the IGVF Portal raw VAMP-seq "
                              "MeasurementSets (tile × bin × replicate × "
                              "antibody coverage).")
    s.add_argument("--label", default=None)
    s.set_defaults(func=cmd_vampseq_inventory)

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
