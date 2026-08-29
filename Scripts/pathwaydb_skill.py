"""pathwaydb — pull current pathway databases and integrate them locally.

Live pathway APIs answer "which pathways contain this gene" one database at a
time, and each names the same biology differently. IntPath (Zhou et al.) showed
what to do about that: normalise gene identifiers to one namespace, map each
source's own relation vocabulary onto a common one, and record *every*
supporting database per fact so agreement between databases stays visible
instead of collapsing into a duplicate row.

This skill applies that method to **current releases**. The published IntPath
dataset is a snapshot from over a decade ago — 582 pathways and 71,116 gene
pairs, against 372 KEGG pathways, ~160k Reactome human memberships and a
WikiPathways release that has moved to a different distribution entirely
today. Its original pipeline cannot simply be re-run: the WikiPathways SOAP
webservice it consumed has been retired (``listPathways`` now 404s). So the
method is reimplemented here against the routes that are actually supported
now, and the release identifier of every download is recorded so a result can
be traced to the exact data it came from.

Sources, all open and key-free:

    KEGG          REST (rest.kegg.jp) — membership, plus typed relations
                  parsed from KGML pathway maps
    Reactome      bulk download — membership (NCBI2Reactome_All_Levels) and
                  curated gene-gene interactions
    WikiPathways  current GMT release (data.wikipathways.org)

BioCyc/HumanCyc, IntPath's third source, is deliberately **not** included:
its flat-file download moved behind a paid subscription, so it cannot be
pulled reproducibly by an open agent. ``pathwaydb sources`` says so rather
than quietly returning two databases where a user expects three.

Identifiers are normalised through NCBI ``gene_info`` — the authoritative
Entrez -> symbol map, which also supplies synonyms, so a database that reports
a deprecated symbol still lands on the current one. Using it rather than any
single pathway database's own gene list keeps the namespace independent of the
sources being merged.

Usage::

    igvfagent pathwaydb pull                  # download current releases
    igvfagent pathwaydb build                 # normalise, integrate, ingest
    igvfagent pathwaydb build --relations     # + typed KEGG/Reactome edges
    igvfagent pathwaydb status                # what is cached, and how old
    igvfagent pathwaydb query --genes CLU,BIN1,PICALM
    igvfagent pathwaydb sources               # coverage and provenance

``build`` writes edges into the local knowledge graph, so pathway structure
pulled once is available to every later query without re-fetching.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

__all__ = ["main", "pull", "build", "integrate", "GeneIndex"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
CACHE_DIR = Path(os.environ.get("IGVF_PATHWAYDB_DIR")
                 or ROOT / "Data" / "Reference" / "PathwayDB")
OUT_DIR = ROOT / "Docs" / "PathwayDB"

KEGG_REST = "https://rest.kegg.jp"
REACTOME_DL = "https://reactome.org/download/current"
WIKIPATHWAYS_GMT = "https://data.wikipathways.org/current/gmt/"
NCBI_GENE_INFO = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
                  "Mammalia/Homo_sapiens.gene_info.gz")

USER_AGENT = "igvfagent-pathwaydb (functional genomics agent; contact via repo)"

# One relation vocabulary across databases, following the KEGG terms IntPath
# normalised onto. Kept distinct rather than collapsed into "interacts_with":
# the difference between a physical interaction and transcriptional regulation
# is most of what makes a typed edge worth storing.
RELATION_EDGE = {
    "PPrel": ("interacts_with", "protein-protein interaction"),
    "ECrel": ("enzyme_relation", "sequential catalysis (shared metabolite)"),
    "GPrel": ("in_same_complex", "gene product / complex membership"),
    "GErel": ("regulates_expression", "gene-expression regulation"),
    "PCrel": ("acts_on_compound", "protein-compound relation"),
}


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def _fetch(url: str, dest: Path, *, timeout: int = 600,
           max_age_days: "float | None" = None) -> "Path | None":
    """Download to ``dest`` unless a fresh copy is already cached.

    Re-downloading a 98 MB Reactome file on every invocation would make the
    skill unusable interactively, so freshness is checked first; ``pull
    --force`` is the way to override.
    """
    if dest.is_file() and dest.stat().st_size and max_age_days is not None:
        age = (time.time() - dest.stat().st_mtime) / 86400.0
        if age <= max_age_days:
            return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        if not tmp.stat().st_size:
            tmp.unlink(missing_ok=True)
            return None
        tmp.replace(dest)      # atomic, so an interrupted pull never leaves
        return dest            # a truncated file that later parses as valid
    except Exception:
        tmp.unlink(missing_ok=True)
        return None


def _sha256(path: Path, *, limit: int = 1 << 26) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        read = 0
        while read < limit:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
            read += len(chunk)
    return h.hexdigest()[:16]


def kegg_release(cache: Path) -> str:
    f = _fetch(f"{KEGG_REST}/info/kegg", cache / "kegg_info.txt",
               timeout=60, max_age_days=7)
    if not f:
        return "unknown"
    text = f.read_text(errors="replace")
    for line in text.splitlines():
        if "Release" in line:
            return line.strip()
    # `/info/kegg` no longer carries a Release line; it dates each database
    # instead. The pathway row's date is the stamp that matters here.
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0] == "pathway":
            return f"Release {parts[-1]}"
    return "unknown"


def wikipathways_release() -> "tuple[str, str] | tuple[None, None]":
    """Find the current human GMT: (filename, release date)."""
    try:
        req = urllib.request.Request(WIKIPATHWAYS_GMT,
                                     headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as r:
            html = r.read().decode("utf-8", errors="replace")
    except Exception:
        return None, None
    names = [n for n in re.findall(r'href="([^"]*\.gmt)"', html)
             if "Homo_sapiens" in n]
    if not names:
        return None, None
    name = sorted(names)[-1]
    m = re.search(r"(\d{8})", name)
    return name, (m.group(1) if m else "unknown")


def pull(cache: Path, *, sources, force: bool = False,
         kgml: bool = False, log=print) -> dict:
    """Download current releases. Returns a manifest of what landed."""
    age = None if force else 7.0
    # Seeded from any existing manifest: a partial pull (`--sources kegg`)
    # must not erase the provenance of files it did not touch, or `status`
    # would under-report a cache that is actually complete.
    manifest: "dict[str, dict]" = {}
    prior = cache / "manifest.json"
    if prior.is_file():
        try:
            manifest = json.loads(prior.read_text())
        except Exception:
            manifest = {}

    def _record(key, path, url, extra=None):
        if path and path.is_file():
            manifest[key] = {"file": str(path.relative_to(ROOT)) if
                             str(path).startswith(str(ROOT)) else str(path),
                             "url": url, "bytes": path.stat().st_size,
                             "sha256_16": _sha256(path), **(extra or {})}
            log(f"  {key:<26} {path.stat().st_size/1e6:>8.1f} MB")
        else:
            log(f"  {key:<26} {'FAILED':>8}")

    log("Normaliser")
    _record("ncbi_gene_info",
            _fetch(NCBI_GENE_INFO, cache / "Homo_sapiens.gene_info.gz",
                   max_age_days=age), NCBI_GENE_INFO)

    if "kegg" in sources:
        rel = kegg_release(cache)
        log(f"KEGG  ({rel})")
        _record("kegg_genes", _fetch(f"{KEGG_REST}/list/hsa",
                                     cache / "kegg_hsa.tsv", max_age_days=age),
                f"{KEGG_REST}/list/hsa", {"release": rel})
        _record("kegg_pathways", _fetch(f"{KEGG_REST}/list/pathway/hsa",
                                        cache / "kegg_pathways.tsv",
                                        max_age_days=age),
                f"{KEGG_REST}/list/pathway/hsa", {"release": rel})
        _record("kegg_links", _fetch(f"{KEGG_REST}/link/pathway/hsa",
                                     cache / "kegg_links.tsv", max_age_days=age),
                f"{KEGG_REST}/link/pathway/hsa", {"release": rel})
        if kgml:
            n = pull_kgml(cache, force=force, log=log)
            manifest["kegg_kgml"] = {"maps": n,
                                     "dir": str((cache / "kgml").name)}

    if "reactome" in sources:
        log("Reactome (current)")
        _record("reactome_pathways",
                _fetch(f"{REACTOME_DL}/NCBI2Reactome_All_Levels.txt",
                       cache / "NCBI2Reactome_All_Levels.txt",
                       max_age_days=age),
                f"{REACTOME_DL}/NCBI2Reactome_All_Levels.txt")
        _record("reactome_hierarchy",
                _fetch(f"{REACTOME_DL}/ReactomePathwaysRelation.txt",
                       cache / "ReactomePathwaysRelation.txt",
                       max_age_days=age),
                f"{REACTOME_DL}/ReactomePathwaysRelation.txt")
        _record("reactome_interactions",
                _fetch(f"{REACTOME_DL}/interactors/"
                       "reactome.homo_sapiens.interactions.tab-delimited.txt",
                       cache / "reactome_interactions.tsv", max_age_days=age),
                f"{REACTOME_DL}/interactors/"
                "reactome.homo_sapiens.interactions.tab-delimited.txt")

    if "wikipathways" in sources:
        name, date = wikipathways_release()
        log(f"WikiPathways ({date or 'unknown release'})")
        if name:
            _record("wikipathways_gmt",
                    _fetch(WIKIPATHWAYS_GMT + name, cache / "wikipathways.gmt",
                           max_age_days=age),
                    WIKIPATHWAYS_GMT + name, {"release": date, "name": name})
        else:
            log(f"  {'wikipathways_gmt':<26} {'FAILED':>8}")

    manifest["_pulled_at"] = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                    time.gmtime())}
    (cache / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def pull_kgml(cache: Path, *, force: bool = False, log=print) -> int:
    """Download KGML for every human KEGG map.

    KGML is where KEGG's *typed* relations live — the membership endpoints
    give only gene-in-pathway. One request per map, paced: KEGG asks for
    polite use and will throttle a burst, which would leave a half-empty
    cache that silently under-reports relations.
    """
    paths = _fetch(f"{KEGG_REST}/list/pathway/hsa", cache / "kegg_pathways.tsv",
                   max_age_days=7.0)
    if not paths:
        return 0
    ids = [l.split("\t")[0].replace("path:", "")
           for l in paths.read_text(errors="replace").splitlines() if l.strip()]
    d = cache / "kgml"
    d.mkdir(parents=True, exist_ok=True)
    got = 0
    for i, pid in enumerate(ids, 1):
        dest = d / f"{pid}.xml"
        if dest.is_file() and dest.stat().st_size and not force:
            got += 1
            continue
        if _fetch(f"{KEGG_REST}/get/{pid}/kgml", dest, timeout=90):
            got += 1
        time.sleep(0.34)
        if i % 50 == 0:
            log(f"  KGML {i}/{len(ids)}")
    log(f"  KGML {got}/{len(ids)} maps cached")
    return got


# ---------------------------------------------------------------------------
# Identifier normalisation
# ---------------------------------------------------------------------------

class GeneIndex:
    """Entrez/symbol/synonym -> current HGNC-style symbol.

    Built from NCBI ``gene_info`` so the namespace is independent of the
    databases being merged. Synonyms are included but never allowed to
    override an official symbol: NAT2's synonym list contains "AAC2", and a
    synonym that collides with another gene's real symbol must lose.
    """

    def __init__(self, path: "Path | None" = None):
        self.by_entrez: "dict[str, str]" = {}
        self.official: "set[str]" = set()
        self.by_alias: "dict[str, str]" = {}
        self.by_ensembl: "dict[str, str]" = {}
        if path and path.is_file():
            self._load(path)

    def _load(self, path: Path) -> None:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[0] != "9606":
                    continue
                entrez, symbol, synonyms = f[1], f[2], f[4]
                if not symbol or symbol == "-":
                    continue
                self.by_entrez[entrez] = symbol
                # dbXrefs carries the Ensembl gene id. Reactome's interactor
                # file ships an entirely empty Entrez column in the current
                # release, so ENSG is the only usable key there.
                if len(f) > 5 and f[5] and f[5] != "-":
                    for x in f[5].split("|"):
                        if x.startswith("Ensembl:"):
                            self.by_ensembl.setdefault(x[8:], symbol)
                self.official.add(symbol.upper())
                if synonyms and synonyms != "-":
                    for s in synonyms.split("|"):
                        s = s.strip().upper()
                        if s and s not in self.by_alias:
                            self.by_alias[s] = symbol

    def resolve(self, token: str) -> "str | None":
        """Map one identifier onto a current symbol, or None."""
        t = (token or "").strip()
        if not t:
            return None
        if t.isdigit():
            return self.by_entrez.get(t)
        for prefix in ("entrez:", "ncbigene:", "geneid:"):
            if t.lower().startswith(prefix):
                return self.by_entrez.get(t.split(":", 1)[1])
        if t.upper().startswith("ENSG"):
            return self.by_ensembl.get(t.split(".")[0])
        u = t.upper()
        if u in self.official:
            return t if t in self.official or t.upper() == u else u
        alias = self.by_alias.get(u)
        if alias:
            return alias
        # Unrecognised but symbol-shaped: keep it rather than drop the row.
        # An unmapped symbol is still a usable fact; a dropped one is silent
        # loss, which is worse in a merge whose whole point is coverage.
        return t if re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]{0,20}", t) else None

    def __len__(self) -> int:
        return len(self.by_entrez)


# ---------------------------------------------------------------------------
# Per-database parsing -> a common membership row
# ---------------------------------------------------------------------------

def _kegg_symbol(field: str) -> str:
    """KEGG's name field is 'SYM; description'.

    The symbol is taken from before the semicolon only. Matching against the
    description instead maps CLU onto 'iron-sulfur cluster', which KEGG's own
    find endpoint really does return.
    """
    return field.split(";")[0].split(",")[0].strip()


def parse_kegg(cache: Path, idx: GeneIndex) -> "list[dict]":
    genes = cache / "kegg_hsa.tsv"
    paths = cache / "kegg_pathways.tsv"
    links = cache / "kegg_links.tsv"
    if not (genes.is_file() and paths.is_file() and links.is_file()):
        return []
    sym = {}
    for line in genes.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) >= 4:
            entrez = f[0][4:] if f[0].startswith("hsa:") else f[0]
            sym[f[0]] = idx.by_entrez.get(entrez) or _kegg_symbol(f[3])
    name = {}
    for line in paths.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            name[f[0].replace("path:", "")] = clean_name(f[1].split(" - ")[0])
    out = []
    for line in links.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) < 2 or not f[1].startswith("path:"):
            continue
        g, pid = sym.get(f[0]), f[1].replace("path:", "")
        if g and pid in name:
            out.append({"pathway": name[pid], "gene": g, "source": "KEGG",
                        "pathway_id": pid})
    return out


def parse_reactome(cache: Path, idx: GeneIndex) -> "list[dict]":
    f = cache / "NCBI2Reactome_All_Levels.txt"
    if not f.is_file():
        return []
    out = []
    with open(f, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 or p[5] != "Homo sapiens":
                continue
            g = idx.resolve(p[0])
            if g:
                out.append({"pathway": clean_name(p[3]), "gene": g,
                            "source": "Reactome", "pathway_id": p[1]})
    return out


def parse_wikipathways(cache: Path, idx: GeneIndex) -> "list[dict]":
    f = cache / "wikipathways.gmt"
    if not f.is_file():
        return []
    out = []
    for line in f.read_text(errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        # GMT descriptor is "Name%WikiPathways_YYYYMMDD%WPnnn%Homo sapiens"
        pathway = clean_name(parts[0].split("%")[0])
        pid = parts[1]
        m = re.search(r"(WP\d+)", parts[0] + " " + pid)
        if m:
            pid = m.group(1)
        for ident in parts[2:]:
            g = idx.resolve(ident)
            if g:
                out.append({"pathway": pathway, "gene": g,
                            "source": "WikiPathways", "pathway_id": pid})
    return out


# ---------------------------------------------------------------------------
# Typed gene-gene relations
# ---------------------------------------------------------------------------

def parse_kegg_relations(cache: Path, idx: GeneIndex) -> "list[dict]":
    """Typed relations from cached KGML maps.

    A KGML entry can name several genes (a component family drawn as one
    box), so a relation between two entries expands to the cross-product of
    their members — which is how KEGG itself renders it.
    """
    d = cache / "kgml"
    if not d.is_dir():
        return []
    pathnames = {}
    pf = cache / "kegg_pathways.tsv"
    if pf.is_file():
        for line in pf.read_text(errors="replace").splitlines():
            f = line.split("\t")
            if len(f) >= 2:
                pathnames[f[0].replace("path:", "")] = f[1].split(" - ")[0].strip()

    out = []
    for xml in sorted(d.glob("*.xml")):
        try:
            root = ET.parse(xml).getroot()
        except Exception:
            continue
        pid = xml.stem
        pname = pathnames.get(pid, root.get("title") or pid)
        members: "dict[str, list[str]]" = {}
        for e in root.findall("entry"):
            if e.get("type") != "gene":
                continue
            syms = []
            for token in (e.get("name") or "").split():
                if token.startswith("hsa:"):
                    g = idx.by_entrez.get(token[4:])
                    if g:
                        syms.append(g)
            gfx = e.find("graphics")
            if not syms and gfx is not None and gfx.get("name"):
                s = idx.resolve(gfx.get("name").split(",")[0].strip())
                if s:
                    syms.append(s)
            if syms:
                members[e.get("id")] = syms
        for rel in root.findall("relation"):
            rtype = rel.get("type") or ""
            if rtype not in RELATION_EDGE:
                continue
            a_list = members.get(rel.get("entry1")) or []
            b_list = members.get(rel.get("entry2")) or []
            subtypes = ",".join(s.get("name") or "" for s in rel.findall("subtype"))
            for a in a_list:
                for b in b_list:
                    if a != b:
                        out.append({"gene_a": a, "gene_b": b,
                                    "relation": rtype, "subtype": subtypes,
                                    "pathway": pname, "source": "KEGG"})
    return out


def parse_reactome_interactions(cache: Path, idx: GeneIndex) -> "list[dict]":
    """Curated gene-gene interactions from Reactome's interactor file."""
    f = cache / "reactome_interactions.tsv"
    if not f.is_file():
        return []
    out = []
    with open(f, encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        cols = [c.strip().lower() for c in header.lstrip("#").split("\t")]

        def _col(*names):
            for n in names:
                for i, c in enumerate(cols):
                    if n in c:
                        return i
            return None

        i_a = _col("interactor 1 entrez")
        i_b = _col("interactor 2 entrez")
        e_a = _col("interactor 1 ensembl")
        e_b = _col("interactor 2 ensembl")
        i_ctx = _col("pathway", "interaction context")
        i_type = _col("interaction type")
        if (i_a is None and e_a is None) or (i_b is None and e_b is None):
            return []

        def _gene(parts, entrez_col, ens_col):
            """Entrez first, Ensembl second.

            The Entrez column is empty throughout the current release, so
            falling back is not an edge case here — it is the normal path.
            The Ensembl field lists transcripts and proteins alongside the
            gene, so only ENSG entries are considered.
            """
            if entrez_col is not None and len(parts) > entrez_col:
                v = parts[entrez_col].strip()
                if v and v != "-":
                    g = idx.resolve(v.split(":")[-1])
                    if g:
                        return g
            if ens_col is not None and len(parts) > ens_col:
                for tok in parts[ens_col].split("|"):
                    tok = tok.split(":")[-1].strip()
                    if tok.startswith("ENSG"):
                        g = idx.resolve(tok)
                        if g:
                            return g
            return None

        for line in fh:
            p = line.rstrip("\n").split("\t")
            a = _gene(p, i_a, e_a)
            b = _gene(p, i_b, e_b)
            if not a or not b or a == b:
                continue
            out.append({"gene_a": a, "gene_b": b, "relation": "PPrel",
                        "subtype": (p[i_type] if i_type is not None
                                    and len(p) > i_type else ""),
                        "pathway": (p[i_ctx] if i_ctx is not None
                                    and len(p) > i_ctx else ""),
                        "source": "Reactome"})
    return out


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

_SUFFIX = re.compile(r"\s*[-(]\s*(homo sapiens|human)[^)]*\)?\s*$", re.I)

# WikiPathways' GMT generator flattens the HTML entity &#39; to a bare " 39 ",
# shipping "Alzheimer 39 s disease". Repaired here rather than left to reach a
# figure legend, and rather than unescaping entities generally — the file
# contains no actual entities to unescape, only this already-broken form.
_APOS = re.compile(r"\s+39\s+s\b")


def clean_name(name: str) -> str:
    return _APOS.sub("'s", (name or "").strip())


def canon_pathway(name: str) -> str:
    """Fold names that differ only cosmetically for the same pathway.

    Conservative on purpose. KEGG's 'Apoptosis - Homo sapiens (human)' and
    WikiPathways' 'Apoptosis' are the same pathway and should merge; two
    genuinely different names for related biology are left separate rather
    than guessed at, because a wrong merge silently fabricates agreement
    between databases, which is exactly the signal this tool reports.
    """
    n = _SUFFIX.sub("", (name or "").strip())
    n = re.sub(r"[^A-Za-z0-9]+", " ", n).strip().lower()
    return n


def merge_by_overlap(rows, *, jaccard: float = 0.7, min_genes: int = 5):
    """Merge pathways across databases whose gene sets nearly coincide.

    Name matching alone is far too weak: it merges only ~40 of KEGG's 354
    pathways against WikiPathways, while Reactome describes the same biology
    with entirely different wording ('Adherens junction' vs 'Adherens
    junctions interactions'). Gene membership is the stronger and more
    defensible criterion — two pathways from different databases annotating
    nearly the same genes are the same pathway, whatever they are called.

    Only *cross-database* pairs are merged. Two pathways within one database
    are that database's own deliberate distinction (Reactome's hierarchy
    nests broad parents over specific children) and collapsing them would
    destroy real structure while inventing corroboration that does not exist.

    Returns ``{(source, canon): group_key}`` for every pathway that merged.
    """
    sets: "dict[tuple, set]" = {}
    for r in rows:
        canon = canon_pathway(r.get("pathway", ""))
        gene = (r.get("gene") or "").strip().upper()
        if canon and gene:
            sets.setdefault((r["source"], canon), set()).add(gene)
    keys = [k for k, v in sets.items() if len(v) >= min_genes]

    by_db: "dict[str, list]" = {}
    for k in keys:
        by_db.setdefault(k[0], []).append(k)

    parent: "dict[tuple, tuple]" = {k: k for k in keys}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            # Deterministic representative, so a rebuild groups identically.
            lo, hi = (ra, rb) if ra <= rb else (rb, ra)
            parent[hi] = lo

    dbs = sorted(by_db)
    for i, da in enumerate(dbs):
        for db in dbs[i + 1:]:
            # Index the larger side, scan the smaller: candidate lookup then
            # touches only pathways that actually share a gene, instead of
            # every pair of the two databases.
            small, large = ((by_db[da], by_db[db])
                            if len(by_db[da]) <= len(by_db[db])
                            else (by_db[db], by_db[da]))
            index: "dict[str, list]" = {}
            for k in large:
                for g in sets[k]:
                    index.setdefault(g, []).append(k)
            for k in small:
                a = sets[k]
                counts: "dict[tuple, int]" = {}
                for g in a:
                    for other in index.get(g, ()):
                        counts[other] = counts.get(other, 0) + 1
                for other, shared in counts.items():
                    b = sets[other]
                    # Jaccard >= t forces |shared| >= t * max(|a|,|b|), so a
                    # cheap size check discards most candidates first.
                    if shared < jaccard * max(len(a), len(b)):
                        continue
                    if shared / float(len(a | b)) >= jaccard:
                        union(k, other)
    return {k: find(k) for k in keys if find(k) != k or
            sum(1 for x in keys if find(x) == k) > 1}


def integrate(rows, *, overlap_jaccard: "float | None" = 0.7) -> "list[dict]":
    """Merge memberships across databases, keeping every supporting source."""
    groups = (merge_by_overlap(rows, jaccard=overlap_jaccard)
              if overlap_jaccard else {})
    # A merged group is addressed by its representative's canonical name, so
    # rows from either database land in the same bucket.
    remap = {k: v[1] for k, v in groups.items()}

    merged: "dict[tuple, dict]" = {}
    display: "dict[str, str]" = {}
    for r in rows:
        canon = canon_pathway(r.get("pathway", ""))
        canon = remap.get((r["source"], canon), canon)
        gene = (r.get("gene") or "").strip()
        if not canon or not gene:
            continue
        # Prefer the shortest display name: KEGG's is the same pathway with a
        # species suffix, and the bare name reads better in a figure legend.
        cur = display.get(canon)
        if cur is None or len(r["pathway"].strip()) < len(cur):
            display[canon] = r["pathway"].strip()
        e = merged.setdefault((canon, gene.upper()),
                              {"gene": gene, "sources": set(),
                               "pathway_ids": set()})
        e["sources"].add(r["source"])
        if r.get("pathway_id"):
            e["pathway_ids"].add(f"{r['source']}:{r['pathway_id']}")
    out = []
    for (canon, _g), e in merged.items():
        out.append({"pathway": display[canon], "pathway_canon": canon,
                    "gene": e["gene"], "sources": sorted(e["sources"]),
                    "n_sources": len(e["sources"]),
                    "pathway_ids": sorted(e["pathway_ids"])})
    return out


def integrate_relations(rows) -> "list[dict]":
    merged: "dict[tuple, dict]" = {}
    for r in rows:
        a, b = r["gene_a"], r["gene_b"]
        key = (min(a.upper(), b.upper()), max(a.upper(), b.upper()),
               r["relation"])
        e = merged.setdefault(key, {"gene_a": a, "gene_b": b,
                                    "relation": r["relation"],
                                    "sources": set(), "pathways": set(),
                                    "subtypes": set()})
        e["sources"].add(r["source"])
        if r.get("pathway"):
            e["pathways"].add(r["pathway"])
        if r.get("subtype"):
            e["subtypes"].add(r["subtype"])
    out = []
    for e in merged.values():
        out.append({"gene_a": e["gene_a"], "gene_b": e["gene_b"],
                    "relation": e["relation"], "sources": sorted(e["sources"]),
                    "n_sources": len(e["sources"]),
                    "subtype": ";".join(sorted(s for s in e["subtypes"] if s))[:200],
                    "pathways": sorted(e["pathways"])[:5]})
    return out


# ---------------------------------------------------------------------------
# Knowledge-graph ingestion
# ---------------------------------------------------------------------------

def prune_stale(con, build_id: str, *, ls) -> int:
    """Drop pathway nodes left behind by an earlier pathwaydb build.

    A pathway's node identity is its name, so when a source corrects a name
    the rebuild writes a new node and the old one is orphaned — visible in
    the KG as a duplicate pathway nobody links to. Every pathway node this
    skill writes carries the id of the build that wrote it, so anything
    stamped by a *previous* pathwaydb build and not by this one is stale.

    Scoped to nodes this skill stamped: pathway nodes written by pathway-viz
    or by an artefact harvest carry no stamp and are never touched.
    """
    stale = []
    for nid, props in con.execute(
            "SELECT id, properties FROM nodes WHERE node_type='pathway'"):
        try:
            d = json.loads(props or "{}")
        except Exception:
            continue
        b = d.get("pathwaydb_build")
        if b and b != build_id:
            stale.append(nid)
    for nid in stale:
        con.execute("DELETE FROM edges WHERE from_node=? OR to_node=?",
                    (nid, nid))
        con.execute("DELETE FROM nodes WHERE id=?", (nid,))
    return len(stale)


def ingest(memberships, relations, *, label: str = "current",
           min_sources: int = 1, release: str = "",
           build_id: str = "", prune: bool = False) -> dict:
    try:
        from igvfagent import _localstore as ls
    except Exception:
        try:
            import _localstore as ls  # type: ignore
        except Exception:
            return {"error": "local store unavailable"}

    added = collections.Counter()
    con = ls._connect()
    try:
        run = ls.upsert_node(
            con, "analysis", f"pathwaydb_{label}", source="igvfagent:pathwaydb",
            label=f"pathwaydb ingest {label}",
            properties={"skill": "pathwaydb", "release": release,
                        "n_memberships": len(memberships),
                        "n_relations": len(relations)})

        for m in memberships:
            if m["n_sources"] < min_sources:
                continue
            # The source is the evidence database, matching what pathway-viz
            # and the harvest recognisers write, so a pathway edge ingested
            # from here collapses onto the same edge rather than duplicating
            # it under a second label.
            src = "+".join(m["sources"])
            g = ls.upsert_node(con, "gene", m["gene"], source=src)
            pw = ls.upsert_node(con, "pathway", m["pathway"], source=src,
                                properties={"pathwaydb_build": build_id})
            ls.upsert_edge(con, g, pw, "member_of_pathway", source=src,
                           properties={"databases": ",".join(m["sources"]),
                                       "n_sources": m["n_sources"],
                                       "pathway_ids": ",".join(m["pathway_ids"])[:300]})
            added["member_of_pathway"] += 1

        for r in relations:
            if r["n_sources"] < min_sources:
                continue
            etype, _d = RELATION_EDGE.get(r["relation"], ("interacts_with", ""))
            src = "+".join(r["sources"])
            a = ls.upsert_node(con, "gene", r["gene_a"], source=src)
            b = ls.upsert_node(con, "gene", r["gene_b"], source=src)
            ls.upsert_edge(con, a, b, etype, source=src,
                           properties={"relation": r["relation"],
                                       "subtype": r.get("subtype", ""),
                                       "databases": ",".join(r["sources"]),
                                       "n_sources": r["n_sources"],
                                       "pathways": ",".join(r.get("pathways", []))[:300]})
            added[etype] += 1

        if prune:
            n = prune_stale(con, build_id, ls=ls)
            if n:
                added["_pruned_stale_pathways"] = n
        ls.upsert_edge(con, run, run, "analyzed", source="igvfagent:pathwaydb")
        con.commit()
    finally:
        con.close()
    return dict(added)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _cache_dir(args) -> Path:
    c = Path(args.cache) if getattr(args, "cache", None) else CACHE_DIR
    return c if c.is_absolute() else (ROOT / c)


def _load_index(cache: Path) -> GeneIndex:
    p = cache / "Homo_sapiens.gene_info.gz"
    if not p.is_file():
        p = _fetch(NCBI_GENE_INFO, p, max_age_days=30.0) or p
    return GeneIndex(p if p.is_file() else None)


def cmd_pull(args) -> int:
    cache = _cache_dir(args)
    sources = {s.strip().lower() for s in args.sources.split(",") if s.strip()}
    print(f"Cache: {cache}\n")
    m = pull(cache, sources=sources, force=args.force, kgml=args.kgml)
    ok = [k for k, v in m.items() if not k.startswith("_") and v.get("bytes")]
    print(f"\n{len(ok)} file(s) cached · manifest at {cache / 'manifest.json'}")
    return 0 if ok else 1


def build(cache: Path, *, sources, relations: bool = False,
          overlap_jaccard: "float | None" = 0.7, log=print) -> dict:
    idx = _load_index(cache)
    log(f"Normaliser:    {len(idx):,} human genes from NCBI gene_info")

    rows, per_source = [], {}
    for key, fn in (("kegg", parse_kegg), ("reactome", parse_reactome),
                    ("wikipathways", parse_wikipathways)):
        if key not in sources:
            continue
        r = fn(cache, idx)
        per_source[key] = len(r)
        rows += r
        log(f"  {key:<14}{len(r):>9,} memberships")

    rel_rows = []
    if relations:
        for key, fn in (("kegg", parse_kegg_relations),
                        ("reactome", parse_reactome_interactions)):
            if key not in sources:
                continue
            r = fn(cache, idx)
            rel_rows += r
            log(f"  {key:<14}{len(r):>9,} typed relations")

    merged = integrate(rows, overlap_jaccard=overlap_jaccard)
    merged_rel = integrate_relations(rel_rows)
    return {"memberships": merged, "relations": merged_rel,
            "per_source": per_source, "genes_indexed": len(idx)}


def cmd_build(args) -> int:
    cache = _cache_dir(args)
    srcs = {s.strip().lower() for s in args.sources.split(",") if s.strip()}
    if not args.no_pull:
        if not (cache / "manifest.json").is_file():
            print("No cache found — pulling first.\n")
            pull(cache, sources=srcs, kgml=args.relations)
            print()
        elif args.relations and "kegg" in srcs and not (cache / "kgml").is_dir():
            # Typed KEGG relations live only in KGML. Without this the build
            # would quietly report Reactome relations alone and look complete.
            print("Typed relations requested but no KGML cached — fetching.\n")
            pull_kgml(cache)
            print()
    res = build(cache, sources=srcs, relations=args.relations,
                overlap_jaccard=(None if args.no_overlap_merge
                                 else args.overlap))
    merged, rel = res["memberships"], res["relations"]
    if not merged:
        print("error: nothing parsed — run `igvfagent pathwaydb pull` first",
              file=sys.stderr)
        return 1

    multi = sum(1 for m in merged if m["n_sources"] > 1)
    pathways = len({m["pathway_canon"] for m in merged})
    genes = len({m["gene"].upper() for m in merged})

    ts = time.strftime("%Y%m%d_%H%M%S")
    run = OUT_DIR / f"{ts}_{args.label}"
    run.mkdir(parents=True, exist_ok=True)
    with open(run / "memberships.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "target", "kind", "evidence", "score", "ids"])
        for m in merged:
            w.writerow([m["gene"], m["pathway"], "pathway",
                        "+".join(m["sources"]), m["n_sources"],
                        ";".join(m["pathway_ids"])])
    if rel:
        with open(run / "relations.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["source", "target", "kind", "evidence", "score", "subtype"])
            for r in rel:
                w.writerow([r["gene_a"], r["gene_b"], r["relation"],
                            "+".join(r["sources"]), r["n_sources"],
                            r.get("subtype", "")])

    release = _release_string(cache)
    summary = {"built_at": ts, "release": release,
               "per_source": res["per_source"],
               "genes_indexed": res["genes_indexed"],
               "memberships": len(merged), "pathways": pathways,
               "genes": genes, "multi_source_memberships": multi,
               "relations": len(rel)}
    (run / "summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\nIntegrated:    {len(merged):,} gene-pathway memberships")
    print(f"  pathways:    {pathways:,}")
    print(f"  genes:       {genes:,}")
    print(f"  corroborated by >1 database: {multi:,} "
          f"({100.0*multi/max(1,len(merged)):.1f}%)")
    if rel:
        by = collections.Counter(r["relation"] for r in rel)
        print(f"  relations:   {len(rel):,} "
              f"({', '.join(f'{k} {v:,}' for k, v in by.most_common())})")

    if not args.no_kg:
        added = ingest(merged, rel, label=args.label,
                       min_sources=args.min_sources, release=release,
                       build_id=ts, prune=not args.no_prune)
        print(f"KG:            {json.dumps(added)}")
    print(f"Manifest:      {run}")
    return 0


def _release_string(cache: Path) -> str:
    m = cache / "manifest.json"
    if not m.is_file():
        return "unknown"
    try:
        d = json.loads(m.read_text())
    except Exception:
        return "unknown"
    bits = []
    kg = (d.get("kegg_links") or {}).get("release")
    if kg:
        bits.append(kg.replace("Release ", "KEGG "))
    wp = (d.get("wikipathways_gmt") or {}).get("release")
    if wp:
        bits.append(f"WikiPathways {wp}")
    if d.get("reactome_pathways"):
        bits.append("Reactome " + (d.get("_pulled_at") or {}).get("utc", "")[:10])
    return " · ".join(bits) or "unknown"


def cmd_status(args) -> int:
    cache = _cache_dir(args)
    print(f"Cache: {cache}")
    m = cache / "manifest.json"
    if not m.is_file():
        print("\nNothing cached. Run: igvfagent pathwaydb pull")
        return 1
    d = json.loads(m.read_text())
    print(f"Pulled: {(d.get('_pulled_at') or {}).get('utc', '?')}")
    print(f"Release: {_release_string(cache)}\n")
    for k, v in d.items():
        if k.startswith("_"):
            continue
        if "bytes" in v:
            p = ROOT / v["file"] if not Path(v["file"]).is_absolute() else Path(v["file"])
            age = ((time.time() - p.stat().st_mtime) / 86400.0
                   if p.is_file() else float("nan"))
            print(f"  {k:<24} {v['bytes']/1e6:>8.1f} MB   {age:>4.1f} d old   "
                  f"{v.get('sha256_16','')}")
        else:
            print(f"  {k:<24} {v}")
    runs = sorted(OUT_DIR.glob("*/summary.json")) if OUT_DIR.is_dir() else []
    if runs:
        s = json.loads(runs[-1].read_text())
        print(f"\nLast build: {s['memberships']:,} memberships · "
              f"{s['pathways']:,} pathways · {s['genes']:,} genes · "
              f"{s.get('relations', 0):,} relations")
    else:
        print("\nNo build yet. Run: igvfagent pathwaydb build")
    return 0


def cmd_sources(args) -> int:
    print("Pathway databases integrated by this skill\n")
    rows = [
        ("KEGG", "rest.kegg.jp", "open, polite rate limit",
         "membership + typed KGML relations"),
        ("Reactome", "reactome.org/download/current", "open, CC-BY",
         "membership + curated interactions"),
        ("WikiPathways", "data.wikipathways.org/current/gmt", "open, CC0",
         "membership (GMT release)"),
    ]
    for name, host, lic, gives in rows:
        print(f"  {name:<14}{host:<38}{lic}")
        print(f"  {'':<14}{gives}")
    print("\n  BioCyc/HumanCyc  not included — flat-file download now requires")
    print("                   a paid subscription, so it cannot be pulled")
    print("                   reproducibly. IntPath's 2011 release included it;")
    print("                   `igvfagent intpath` still serves that snapshot.")
    print("\n  Identifiers normalised via NCBI gene_info (Entrez -> symbol,")
    print("  with synonyms), independent of the databases being merged.")
    return 0


def cmd_query(args) -> int:
    cache = _cache_dir(args)
    genes = [g.strip() for g in (args.genes or "").replace(",", " ").split() if g.strip()]
    if not genes:
        print("error: --genes required", file=sys.stderr)
        return 2
    idx = _load_index(cache)
    want = {(idx.resolve(g) or g).upper() for g in genes}

    sources = {s.strip().lower() for s in args.sources.split(",") if s.strip()}
    rows = []
    for key, fn in (("kegg", parse_kegg), ("reactome", parse_reactome),
                    ("wikipathways", parse_wikipathways)):
        if key in sources:
            rows += [r for r in fn(cache, idx) if r["gene"].upper() in want]
    if not rows:
        print("No cached data for those genes. Run: igvfagent pathwaydb pull")
        return 1

    merged = integrate(rows, overlap_jaccard=None)
    by_pathway: "dict[str, dict]" = {}
    for m in merged:
        e = by_pathway.setdefault(m["pathway"], {"genes": set(), "sources": set()})
        e["genes"].add(m["gene"].upper())
        e["sources"].update(m["sources"])

    ordered = sorted(by_pathway.items(),
                     key=lambda kv: (-len(kv[1]["genes"]), -len(kv[1]["sources"]),
                                     kv[0]))
    print(f"{len(want)} gene(s) · {len(ordered)} pathway(s)\n")
    print(f"  {'genes':<6} {'db':<24} pathway")
    for name, e in ordered[:args.top]:
        hits = ",".join(sorted(e["genes"]))
        print(f"  {len(e['genes']):<6} {'+'.join(sorted(e['sources'])):<24} "
              f"{name[:70]}")
        if args.verbose:
            print(f"  {'':<31} {hits}")
    if len(ordered) > args.top:
        print(f"\n  ... {len(ordered) - args.top} more (raise --top)")

    shared = [n for n, e in ordered if len(e["genes"]) >= max(2, len(want) // 2)]
    if shared:
        print(f"\nShared by >= {max(2, len(want)//2)} of the query genes: "
              f"{len(shared)}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent pathwaydb",
        description="Pull current KEGG / Reactome / WikiPathways releases, "
                    "normalise identifiers, integrate them, and store the "
                    "result in the local knowledge graph.")
    sub = p.add_subparsers(dest="cmd", required=True)

    def _common(s, *, with_sources=True):
        if with_sources:
            s.add_argument("--sources", default="kegg,reactome,wikipathways",
                           help="Comma-separated subset to use")
        s.add_argument("--cache", help="Download cache directory "
                                       "(default Data/Reference/PathwayDB)")
        return s

    pu = _common(sub.add_parser("pull", help="Download current releases"))
    pu.add_argument("--force", action="store_true",
                    help="Re-download even if the cache is fresh")
    pu.add_argument("--kgml", action="store_true",
                    help="Also cache KGML maps (typed KEGG relations; "
                         "~370 requests, a few minutes)")

    b = _common(sub.add_parser("build", help="Normalise, integrate, ingest"))
    b.add_argument("--relations", action="store_true",
                   help="Also build typed gene-gene relations")
    b.add_argument("--label", default="current")
    b.add_argument("--min-sources", type=int, default=1,
                   help="Only ingest facts supported by >= N databases")
    b.add_argument("--no-kg", action="store_true",
                   help="Write files but do not touch the knowledge graph")
    b.add_argument("--overlap", type=float, default=0.7,
                   help="Gene-set Jaccard above which pathways from different "
                        "databases are treated as the same pathway")
    b.add_argument("--no-overlap-merge", action="store_true",
                   help="Merge on pathway name only, never on gene overlap")
    b.add_argument("--no-prune", action="store_true",
                   help="Keep pathway nodes from earlier pathwaydb builds "
                        "that this build no longer produces")
    b.add_argument("--no-pull", action="store_true",
                   help="Fail rather than downloading if the cache is empty")

    _common(sub.add_parser("status", help="What is cached and how old"),
            with_sources=False)
    sub.add_parser("sources", help="Databases, licences and coverage")

    q = _common(sub.add_parser("query", help="Pathways for a gene list"))
    q.add_argument("--genes", required=True)
    q.add_argument("--top", type=int, default=25)
    q.add_argument("--verbose", action="store_true",
                   help="Also list which query genes hit each pathway")

    args = p.parse_args(argv)
    return {"pull": cmd_pull, "build": cmd_build, "status": cmd_status,
            "query": cmd_query, "sources": cmd_sources}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
