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


def parse_gmt(path: Path, label: str, idx: GeneIndex) -> "list[dict]":
    """Read any GMT as one more source database.

    This is how a database that cannot be fetched anonymously still gets
    integrated. BioCyc/HumanCyc -- IntPath's third source -- is the case that
    motivated it: its bulk download now requires a paid subscription, so
    IGVFagent cannot pull it, but a user who holds a licence can export the
    pathways they have rights to and point this at the file. The export then
    goes through exactly the same identifier normalisation, unification and
    knowledge-graph ingestion as the fetched sources, and shows up in the
    per-fact source list like any other database.

    Nothing licensed is redistributed: the file stays on the user's machine.
    """
    if not path.is_file():
        return []
    out = []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        pathway = clean_name(parts[0].split("%")[0])
        pid = parts[1].strip()
        for ident in parts[2:]:
            g = idx.resolve(ident.strip())
            if g:
                out.append({"pathway": pathway, "gene": g, "source": label,
                            "pathway_id": pid})
    return out


def _parse_extra(specs, idx) -> "tuple[list, dict]":
    """``--extra-gmt LABEL=path`` entries -> membership rows."""
    rows, counts = [], {}
    for spec in specs or ():
        if "=" not in spec:
            print(f"warning: --extra-gmt expects LABEL=path, got {spec!r}",
                  file=sys.stderr)
            continue
        label, _, raw = spec.partition("=")
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = ROOT / raw
        r = parse_gmt(path, label.strip(), idx)
        if not r:
            print(f"warning: no gene sets read from {path}", file=sys.stderr)
        counts[label.strip()] = len(r)
        rows += r
    return rows, counts

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


# ---------------------------------------------------------------------------
# Pathway-name unification — IntPath's method (Zhou et al. 2012, BMC Syst Biol
# 6(Suppl 2):S2), reimplemented from the published description
# ---------------------------------------------------------------------------
#
# IntPath identifies two pathway names as the same pathway by longest common
# subsequence, and reports LCS as more reliable for this than gene-set
# overlap. Its acceptance rule, verbatim from the paper:
#
#     alignment score = number of aligned characters (the LCS length)
#     alignment ratio = 2 x score / (len(a) + len(b))
#
#     accept if EITHER
#       (1) score > len(shorter) - 1  AND  ratio >= 0.5
#       (2) ratio > 0.91
#
# then filter with an "error-prone words pair list", group with a disjoint-set
# structure, and name each group by its SHORTEST member name.

# Pairs of words that make two different pathways look alike. The paper names
# the first two; "VEGF signaling pathway" vs "EGFR1 Signaling Pathway" scores
# ratio 0.933, which clears rule (2) on its own, so without this filter the
# rule actively merges them. Extend via --error-pairs.
ERROR_PRONE_PAIRS = [
    ("egfr", "vegf"), ("t cell", "b cell"), ("type i", "type ii"),
    ("alpha", "beta"), ("her2", "egfr"), ("insulin", "igf"),
    ("mapk", "erk5"), ("wnt", "shh"), ("tgf beta", "bmp"),
    # Measured false merges: "eukaryotic transcription initiation" scores
    # ratio 0.94 against "eukaryotic translation initiation".
    ("transcription", "translation"), ("replication", "repair"),
    ("anabolism", "catabolism"), ("import", "export"),
    ("il 2", "il 4"), ("smooth muscle", "cardiac muscle"),
    ("male", "female"), ("mitochondrial", "cytosolic"),
]

# Words that invert or restrict a pathway's meaning. The paper's error-prone
# list is pairwise; this is the unary case, and it is needed because modern
# Reactome names the negative and defective forms of a pathway explicitly:
# "Apoptosis" and "Suppression of apoptosis" satisfy the published rule (the
# shorter is a subsequence of the longer, ratio 0.58 >= 0.5) yet mean opposite
# things. If one name carries such a qualifier and the other does not, they
# are different pathways.
QUALIFIERS = (
    "suppression", "defective", "non", "negative", "positive", "anti",
    "resistance", "deficiency", "aberrant", "abnormal", "loss", "gain",
    "inhibition", "dysregulated", "impaired", "reduced", "escape",
    # Measured: "downregulation of tgf beta receptor signaling" merged with
    # "tgf beta receptor signaling", and "diseases of base excision repair"
    # with "base excision repair".
    "downregulation", "upregulation", "diseases", "defects", "novo",
    "checkpoints", "mitotic", "meiotic",
    # Stage words: "dna replication" and "dna replication pre initiation" are
    # a pathway and one phase of it. Only a one-sided occurrence rejects, so
    # two names that both name the same stage still merge.
    "initiation", "preinitiation", "elongation", "termination", "formation",
)

# Tokens carrying a number distinguish members of a family: vitamin B12 is
# not vitamin B6, IL-2 is not IL-4, type I is not type II. Character-level
# similarity cannot see this -- "vitamin b12 metabolism" and "vitamin b6
# metabolism" align almost perfectly -- so it is checked explicitly.
_ROMAN = {"i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x"}


def _numeric_mismatch(a: str, b: str) -> bool:
    def marks(name):
        out = set()
        for t in name.split():
            if any(ch.isdigit() for ch in t) or t in _ROMAN:
                out.add(t)
        return out
    ma, mb = marks(a), marks(b)
    return bool(ma or mb) and ma != mb


def lcs_length(a: str, b: str) -> int:
    """Longest common subsequence length, bit-parallel.

    The plain dynamic program is O(len(a) x len(b)) per pair, and unification
    compares millions of name pairs; this is the Crochemore/Iliopoulos/Pinzon
    bit-vector formulation, which does a machine word of the DP table at a
    time. Verified identical to the DP over random strings.
    """
    m = len(a)
    if not m or not b:
        return 0
    pm: "dict[str, int]" = {}
    for i, ch in enumerate(a):
        pm[ch] = pm.get(ch, 0) | (1 << i)
    full = (1 << m) - 1
    v = full
    for ch in b:
        u = v & pm.get(ch, 0)
        v = ((v + u) | (v - u)) & full
    return m - bin(v).count("1")


def alignment_ratio(score: int, a: str, b: str) -> float:
    n = len(a) + len(b)
    return (2.0 * score / n) if n else 0.0


def _error_prone(a: str, b: str, pairs) -> bool:
    """True if the pair trips the error-prone word list.

    Checked in both directions: one name carrying one partner while the other
    carries the other partner is the signature of a near-name collision
    between genuinely different pathways.
    """
    for w1, w2 in pairs:
        if (w1 in a and w2 in b and w1 not in b and w2 not in a):
            return True
        if (w2 in a and w1 in b and w2 not in b and w1 not in a):
            return True
    return False


def names_related(a: str, b: str, *, pairs=ERROR_PRONE_PAIRS) -> bool:
    """IntPath's two-condition acceptance rule, plus the mismatch filter."""
    if a == b:
        return True
    shorter = min(len(a), len(b))
    longer = max(len(a), len(b))
    # Both rules need score >= (len(a)+len(b))/4 at best, and score can never
    # exceed the shorter length, so wildly different lengths cannot qualify.
    if shorter * 4 < (shorter + longer):
        return False
    score = lcs_length(a, b)
    ratio = alignment_ratio(score, a, b)
    return _accept(a, b, score, ratio, pairs)


def _qualifier_mismatch(a: str, b: str) -> bool:
    ta, tb = set(a.split()), set(b.split())
    for q in QUALIFIERS:
        if (q in ta) != (q in tb):
            return True
    return False


def _accept(a: str, b: str, score: int, ratio: float, pairs) -> bool:
    """The published acceptance rule plus two guards it needs on current data.

    Rule (1) — the shorter name is a subsequence of the longer, lengths within
    about 3x — is permissive when the shorter name is a single generic word.
    Reactome (which IntPath did not include) has top-level pathways literally
    named "Disease" and "Metabolism", and "disease" is a subsequence of
    "alzheimer disease" at ratio 0.58, so the rule as published merges
    Alzheimer, Chagas and prion disease into one group through that shared
    parent. Requiring two tokens before rule (1) applies confines
    single-word names to the much stricter rule (2).
    """
    shorter_name = a if len(a) <= len(b) else b
    shorter = len(shorter_name)
    multiword = len(shorter_name.split()) >= 2
    ok = ((multiword and score > shorter - 1 and ratio >= 0.5)
          or (ratio > 0.91))
    if not ok:
        return False
    if _qualifier_mismatch(a, b) or _numeric_mismatch(a, b):
        return False
    return not _error_prone(a, b, pairs)


def reactome_ancestry(cache: Path) -> "set[tuple]":
    """Canonical-name pairs that Reactome itself declares parent and child.

    Two names in an ancestor/descendant relation are a deliberate distinction
    by the source curators, not two labels for one pathway, so they must never
    be unified however similar the strings look. Reactome publishes the
    hierarchy, so this is read rather than inferred — the heuristics only need
    to cover cases where no database states the answer.
    """
    rel = cache / "ReactomePathwaysRelation.txt"
    names_f = cache / "NCBI2Reactome_All_Levels.txt"
    if not (rel.is_file() and names_f.is_file()):
        return set()

    name: "dict[str, str]" = {}
    with open(names_f, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6 and f[5] == "Homo sapiens" and f[1] not in name:
                name[f[1]] = canon_pathway(clean_name(f[3]))

    children: "dict[str, list]" = {}
    for line in rel.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[0].startswith("R-HSA") and f[1].startswith("R-HSA"):
            children.setdefault(f[0], []).append(f[1])

    # Transitive closure, so a grandparent is blocked from merging with a
    # grandchild too. The human hierarchy is small enough to walk directly.
    pairs = set()
    for root in list(children):
        stack, seen = list(children.get(root, ())), set()
        while stack:
            c = stack.pop()
            if c in seen:
                continue
            seen.add(c)
            a, b = name.get(root), name.get(c)
            if a and b and a != b:
                pairs.add((a, b))
                pairs.add((b, a))
            stack.extend(children.get(c, ()))
    return pairs

def merge_by_name(rows, *, within_db: bool = False,
                  pairs=ERROR_PRONE_PAIRS, forbidden=None, report=None):
    """Group pathway names by IntPath's best-hit + disjoint-set unification.

    For each name in one database's list, the best hit in another database's
    list is the candidate with the highest alignment ratio; the acceptance
    rule then decides whether that best hit is actually the same pathway.

    ``within_db`` also compares names inside one database, which IntPath does
    deliberately. It is off by default here because Reactome — which IntPath
    did not include — is explicitly hierarchical, so 'Signaling by EGFR' and
    'Signaling by EGFR in Cancer' are a parent and child the database means to
    keep apart, not two names for one thing.
    """
    keys: "dict[str, set]" = {}
    for r in rows:
        canon = canon_pathway(r.get("pathway", ""))
        if canon:
            keys.setdefault(r["source"], set()).add(canon)
    by_db = {db: sorted(v) for db, v in keys.items()}

    parent: "dict[tuple, tuple]" = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = (ra, rb) if ra <= rb else (rb, ra)
            parent[hi] = lo

    for db, names in by_db.items():
        for n in names:
            find((db, n))

    dbs = sorted(by_db)
    comparisons = []
    for i, da in enumerate(dbs):
        for dbb in dbs[i + 1:]:
            comparisons.append((da, dbb))
        if within_db:
            comparisons.append((da, da))

    for da, dbb in comparisons:
        xs, ys = by_db[da], by_db[dbb]
        for x in xs:
            best, best_ratio = None, 0.0
            lx = len(x)
            for y in ys:
                if da == dbb and y == x:
                    continue
                ly = len(y)
                # Cheap length gate before the LCS: the acceptance rule can
                # never fire when one name is more than ~3x the other.
                if min(lx, ly) * 4 < lx + ly:
                    continue
                sc = lcs_length(x, y)
                ra = alignment_ratio(sc, x, y)
                if ra > best_ratio:
                    best, best_ratio, best_sc = y, ra, sc
            if best is None:
                continue
            if forbidden and (x, best) in forbidden:
                continue
            if _accept(x, best, best_sc, best_ratio, pairs):
                union((da, x), (dbb, best))
                if report is not None:
                    report.append({"a_db": da, "a": x, "b_db": dbb,
                                   "b": best, "score": best_sc,
                                   "ratio": round(best_ratio, 4)})

    groups: "dict[tuple, tuple]" = {}
    for k in list(parent):
        groups[k] = find(k)
    return groups

# ---------------------------------------------------------------------------
# Alternative unification criteria
# ---------------------------------------------------------------------------
#
# LCS is IntPath's published choice and is the default, but it is one signal
# among several and it is purely lexical: it cannot tell that two differently
# worded names annotate the same genes, and it cannot tell that two similarly
# worded names annotate different ones. These are the other criteria worth
# measuring against it. None is assumed better -- `pathwaydb evaluate` scores
# them on the same data.

def pathway_gene_sets(rows) -> "dict[tuple, set]":
    sets: "dict[tuple, set]" = {}
    for r in rows:
        canon = canon_pathway(r.get("pathway", ""))
        gene = (r.get("gene") or "").strip().upper()
        if canon and gene:
            sets.setdefault((r["source"], canon), set()).add(gene)
    return sets


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def overlap_coefficient(a: set, b: set) -> float:
    """|A n B| / min(|A|,|B|) -- Szymkiewicz-Simpson.

    Tolerates size asymmetry where Jaccard does not. That matters here: a
    Reactome sub-pathway of 12 genes fully contained in a KEGG pathway of 200
    scores 1.0 here and 0.06 by Jaccard. Whether that should count as "the
    same pathway" is exactly the question -- it is a containment signal, not
    an equivalence signal, so it is offered as a distinct criterion rather
    than folded into the Jaccard one.
    """
    m = min(len(a), len(b))
    return (len(a & b) / m) if m else 0.0


def _log_hypergeom_sf(k: int, K: int, n: int, N: int) -> float:
    """log P(X >= k) for the hypergeometric, in log space.

    Computed with lgamma rather than factorials so pathway sizes in the
    thousands do not overflow, and without scipy so the base install stays
    dependency-free -- the same reason the enrichment code works this way.
    """
    from math import lgamma, exp, log
    def lc(n_, r_):
        if r_ < 0 or r_ > n_:
            return float("-inf")
        return lgamma(n_ + 1) - lgamma(r_ + 1) - lgamma(n_ - r_ + 1)
    denom = lc(N, n)
    total = 0.0
    hi = min(K, n)
    for i in range(k, hi + 1):
        t = lc(K, i) + lc(N - K, n - i) - denom
        total += exp(t)
    return log(total) if total > 0 else -745.0


def token_similarity(a: str, b: str) -> float:
    """Word-level Jaccard over content tokens.

    Character-level LCS is blind to word order and happily aligns letters
    across unrelated words, which is how "disease" scores 0.58 against
    "alzheimer disease". Comparing word sets instead removes that failure
    mode, at the cost of missing pure spelling variants.
    """
    stop = {"pathway", "pathways", "signaling", "signalling", "of", "the",
            "in", "and", "a", "by"}
    ta = {t for t in a.split() if t not in stop}
    tb = {t for t in b.split() if t not in stop}
    if not ta or not tb:
        return 0.0
    return len(ta & tb) / len(ta | tb)


def candidate_pairs(sets, *, cross_db_only: bool = True):
    """Pathway pairs sharing at least one gene, keyed for the gene criteria.

    Only pairs with a shared gene can score above zero on any gene-based
    criterion, so an inverted index over genes enumerates every candidate
    without touching the quadratic number of pairs that share nothing.
    """
    index: "dict[str, list]" = {}
    for k, genes in sets.items():
        for g in genes:
            index.setdefault(g, []).append(k)
    seen = set()
    for keys in index.values():
        if len(keys) > 400:      # a gene in hundreds of pathways contributes
            continue             # noise, not signal, and dominates the cost
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                if cross_db_only and a[0] == b[0]:
                    continue
                pair = (a, b) if a <= b else (b, a)
                if pair not in seen:
                    seen.add(pair)
                    yield pair


def merge_by_genes(rows, *, criterion: str = "jaccard", cutoff: float = 0.7,
                   within_db: bool = False, forbidden=None, report=None):
    """Unify pathways whose gene sets agree, by one of three criteria."""
    sets = pathway_gene_sets(rows)
    universe = len({g for v in sets.values() for g in v})
    parent: "dict[tuple, tuple]" = {k: k for k in sets}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = (ra, rb) if ra <= rb else (rb, ra)
            parent[hi] = lo

    for a, b in candidate_pairs(sets, cross_db_only=not within_db):
        ga, gb = sets[a], sets[b]
        if len(ga) < 5 or len(gb) < 5:
            continue
        if forbidden and (a[1], b[1]) in forbidden:
            continue
        if criterion == "jaccard":
            score = jaccard(ga, gb)
            ok = score >= cutoff
        elif criterion == "containment":
            score = overlap_coefficient(ga, gb)
            ok = score >= cutoff
        elif criterion == "hypergeom":
            shared = len(ga & gb)
            if shared < 3:
                continue
            score = -_log_hypergeom_sf(shared, len(ga), len(gb), universe)
            # cutoff is -log(p); 30 is about p < 1e-13, strict enough that
            # chance overlap between large pathways does not qualify.
            ok = score >= cutoff
        else:
            raise ValueError(f"unknown criterion: {criterion}")
        if ok:
            union(a, b)
            if report is not None:
                report.append({"a_db": a[0], "a": a[1], "b_db": b[0],
                               "b": b[1], "score": round(score, 4),
                               "ratio": round(score, 4)})
    return {k: find(k) for k in sets}


def merge_by_tokens(rows, *, cutoff: float = 0.75, within_db: bool = False,
                    forbidden=None, report=None):
    sets = pathway_gene_sets(rows)
    names: "dict[str, list]" = {}
    for db, canon in sets:
        names.setdefault(db, []).append(canon)
    parent: "dict[tuple, tuple]" = {k: k for k in sets}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = (ra, rb) if ra <= rb else (rb, ra)
            parent[hi] = lo

    dbs = sorted(names)
    for i, da in enumerate(dbs):
        others = dbs[i + 1:] + ([da] if within_db else [])
        for dbb in others:
            for x in names[da]:
                for y in names[dbb]:
                    if da == dbb and y <= x:
                        continue
                    if forbidden and (x, y) in forbidden:
                        continue
                    sc = token_similarity(x, y)
                    if sc >= cutoff and not _qualifier_mismatch(x, y) \
                            and not _numeric_mismatch(x, y) \
                            and not _error_prone(x, y, ERROR_PRONE_PAIRS):
                        union((da, x), (dbb, y))
                        if report is not None:
                            report.append({"a_db": da, "a": x, "b_db": dbb,
                                           "b": y, "score": round(sc, 4),
                                           "ratio": round(sc, 4)})
    return {k: find(k) for k in sets}

def merge_consensus(rows, *, within_db: bool = False, forbidden=None,
                    report=None, votes_needed: int = 2,
                    pairs=ERROR_PRONE_PAIRS):
    """Unify pathways only where independent signals agree.

    Measured on the current releases, no single criterion is sufficient:

      LCS name similarity   highest recall, but its false positives are all
                            specialisations ("cell cycle" vs "cell cycle
                            mitotic") that read as near-identical strings
      token overlap         no hierarchy violations, but misses spelling
                            variants that share no whole word
      gene-set Jaccard      does not separate the two classes at all -- true
                            merges run as low as 0.034 while nested pairs
                            reach 0.471, so any threshold trades one error
                            for the other

    So this asks three cheap, genuinely independent questions -- do the names
    align, do they share words, do they annotate the same genes -- and merges
    on agreement rather than on any one of them. The candidate set still comes
    from the LCS rule, which is what keeps recall; the votes are what keep
    precision.

    The hard constraint comes first and is not a vote: when a source database
    declares two pathways to be a parent and a child, that is curated fact and
    outranks every similarity score computed here.
    """
    sets = pathway_gene_sets(rows)
    by_db: "dict[str, list]" = {}
    for db, canon in sets:
        by_db.setdefault(db, []).append(canon)
    for db in by_db:
        by_db[db].sort()

    parent: "dict[tuple, tuple]" = {k: k for k in sets}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = (ra, rb) if ra <= rb else (rb, ra)
            parent[hi] = lo

    dbs = sorted(by_db)
    comparisons = [(dbs[i], dbs[j]) for i in range(len(dbs))
                   for j in range(i + 1, len(dbs))]
    if within_db:
        comparisons += [(d, d) for d in dbs]

    for da, dbb in comparisons:
        for x in by_db[da]:
            lx = len(x)
            best, best_ratio, best_sc = None, 0.0, 0
            best_tok, best_tok_score = None, 0.0
            for y in by_db[dbb]:
                if da == dbb and y == x:
                    continue
                ly = len(y)
                tk = token_similarity(x, y)
                if tk > best_tok_score:
                    best_tok, best_tok_score = y, tk
                if min(lx, ly) * 4 < lx + ly:
                    continue
                sc = lcs_length(x, y)
                ra = alignment_ratio(sc, x, y)
                if ra > best_ratio:
                    best, best_ratio, best_sc = y, ra, sc

            # Two independent nominations: the closest name alignment and the
            # closest word overlap. LCS alone never proposes "signaling by
            # hippo" for "hippo signaling pathway" -- the words are reordered,
            # so the characters do not align -- and word overlap alone never
            # proposes a pure spelling variant. Voting on the union of both
            # candidate sets is what lifts recall without lowering precision.
            for cand in {c for c in (best, best_tok) if c}:
                if forbidden and (x, cand) in forbidden:
                    continue
                if _qualifier_mismatch(x, cand) or _numeric_mismatch(x, cand) \
                        or _error_prone(x, cand, pairs):
                    continue
                sc = lcs_length(x, cand)
                ra = alignment_ratio(sc, x, cand)
                ga, gb = sets[(da, x)], sets[(dbb, cand)]
                tok = token_similarity(x, cand)
                gj = jaccard(ga, gb) if ga and gb else 0.0
                # Identical names need no corroboration -- there is nothing to
                # disambiguate -- so they pass on the name vote alone.
                if x == cand:
                    votes = 3
                else:
                    votes = sum((ra > 0.91, tok >= 0.5, gj >= 0.20))
                if votes >= votes_needed:
                    union((da, x), (dbb, cand))
                    if report is not None:
                        report.append({"a_db": da, "a": x, "b_db": dbb,
                                       "b": cand, "score": sc,
                                       "ratio": round(ra, 4),
                                       "tokens": round(tok, 3),
                                       "gene_jaccard": round(gj, 3),
                                       "votes": votes})
    return {k: find(k) for k in sets}

METHODS = ("consensus", "lcs", "tokens", "jaccard", "containment",
           "hypergeom", "none")

_DEFAULT_CUTOFF = {"jaccard": 0.7, "containment": 0.8, "hypergeom": 30.0,
                   "tokens": 0.75}


def _run_method(m: str, rows, *, within_db=False, forbidden=None, report=None,
                cutoff=None):
    if m == "none":
        return {}
    if m == "consensus":
        return merge_consensus(rows, within_db=within_db, forbidden=forbidden,
                               report=report)
    if m == "lcs":
        return merge_by_name(rows, within_db=within_db, forbidden=forbidden,
                             report=report)
    if m == "tokens":
        return merge_by_tokens(rows, cutoff=cutoff or _DEFAULT_CUTOFF["tokens"],
                               within_db=within_db, forbidden=forbidden,
                               report=report)
    if m in ("jaccard", "containment", "hypergeom"):
        return merge_by_genes(rows, criterion=m,
                              cutoff=cutoff or _DEFAULT_CUTOFF[m],
                              within_db=within_db, forbidden=forbidden,
                              report=report)
    raise ValueError(f"unknown method: {m}")


def integrate(rows, *, method: str = "consensus", within_db: bool = False,
              forbidden=None, report=None) -> "list[dict]":
    """Merge memberships across databases, keeping every supporting source.

    ``method`` names one criterion, or several joined by ``+`` to take their
    union (e.g. ``lcs+jaccard``). ``pathwaydb evaluate`` measures each one on
    the current data rather than assuming which is best.
    """
    groups: "dict[tuple, tuple]" = {}
    for m in [x.strip() for x in method.split("+") if x.strip()]:
        g = _run_method(m, rows, within_db=within_db, forbidden=forbidden,
                        report=report)
        # Each criterion returns its own disjoint-set map over the same keys;
        # composing them takes the union of the criteria, which is what
        # "lcs+jaccard" should mean. Intersection is a different question and
        # is answered by `evaluate`, not by silently changing this.
        for k, v in g.items():
            if v == k:
                continue
            groups[k] = groups.get(v, v)
    # A merged group is addressed by its representative's canonical name, so
    # rows from either database land in the same bucket.
    remap = {k: v[1] for k, v in groups.items() if v != k}

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
          method: str = "consensus", within_db: bool = False,
          use_hierarchy: bool = True, extra_gmt=None, log=print) -> dict:
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

    extra_rows, extra_counts = _parse_extra(extra_gmt, idx)
    for label, n in extra_counts.items():
        log(f"  {label:<14}{n:>9,} memberships (local)")
    rows += extra_rows
    per_source.update(extra_counts)

    rel_rows = []
    if relations:
        for key, fn in (("kegg", parse_kegg_relations),
                        ("reactome", parse_reactome_interactions)):
            if key not in sources:
                continue
            r = fn(cache, idx)
            rel_rows += r
            log(f"  {key:<14}{len(r):>9,} typed relations")

    forbidden = reactome_ancestry(cache) if use_hierarchy else set()
    if forbidden:
        log(f"  hierarchy    {len(forbidden)//2:>9,} parent/child pairs "
            f"protected from merging")
    report: "list[dict]" = []
    merged = integrate(rows, method=method, within_db=within_db,
                       forbidden=forbidden, report=report)
    merged_rel = integrate_relations(rel_rows)
    return {"memberships": merged, "relations": merged_rel,
            "per_source": per_source, "genes_indexed": len(idx),
            "merge_report": report}


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
                method=args.merge, within_db=args.within_db,
                use_hierarchy=not args.no_hierarchy,
                extra_gmt=args.extra_gmt)
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

    rep = res.get("merge_report") or []
    if rep:
        with open(run / "merges.csv", "w", newline="", encoding="utf-8") as fh:
            cols = []
            for r in rep:
                for k in r:
                    if k not in cols:
                        cols.append(k)
            w = csv.DictWriter(fh, fieldnames=cols, restval="")
            w.writeheader()
            w.writerows(rep)

    release = _release_string(cache)
    summary = {"built_at": ts, "release": release,
               "per_source": res["per_source"],
               "genes_indexed": res["genes_indexed"],
               "memberships": len(merged), "pathways": pathways,
               "genes": genes, "multi_source_memberships": multi,
               "relations": len(rel), "name_merges": len(rep)}
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

    merged = integrate(rows, method=args.merge)
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


def cmd_evaluate(args) -> int:
    """Score each unification criterion on the current data.

    There is no gold standard for "these two pathways are the same", so this
    reports two signals that are objective rather than a single invented
    accuracy number:

      hierarchy violations  merges that Reactome itself declares to be a
                            parent and a child. These are known-wrong: the
                            source curators deliberately kept them apart.
                            Lower is better; this is a precision signal.

      identical-name recall of the cross-database pairs whose names are
                            already identical -- unambiguously the same
                            pathway -- how many does the criterion find?
                            This is a recall signal, and it is BIASED IN
                            FAVOUR of the name-based criteria, which get it
                            for free. It is meaningful for comparing the
                            gene-based criteria with each other.

    Both are proxies and are reported as such. The guard is disabled during
    evaluation so violations can be counted rather than silently prevented.
    """
    cache = _cache_dir(args)
    idx = _load_index(cache)
    srcs = {x.strip().lower() for x in args.sources.split(",") if x.strip()}
    rows = []
    for key, fn in (("kegg", parse_kegg), ("reactome", parse_reactome),
                    ("wikipathways", parse_wikipathways)):
        if key in srcs:
            rows += fn(cache, idx)
    if not rows:
        print("error: nothing cached — run `igvfagent pathwaydb pull` first",
              file=sys.stderr)
        return 1

    forbidden = reactome_ancestry(cache)
    sets = pathway_gene_sets(rows)
    by_db: "dict[str, set]" = {}
    for db, canon in sets:
        by_db.setdefault(db, set()).add(canon)
    dbs = sorted(by_db)
    identical = set()
    for i, da in enumerate(dbs):
        for dbb in dbs[i + 1:]:
            for n in by_db[da] & by_db[dbb]:
                identical.add(frozenset({(da, n), (dbb, n)}))

    print(f"Pathways:      {len(sets):,} across {len(dbs)} databases")
    print(f"Reactome hierarchy: {len(forbidden)//2:,} parent/child pairs "
          f"(known-different)")
    print(f"Identical names across databases: {len(identical):,} pairs "
          f"(known-same)\n")

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    results = {}
    print(f"  {'criterion':<14}{'merges':>9}{'hier.viol':>11}"
          f"{'viol.rate':>11}{'ident.found':>13}{'recall':>9}")
    for m in methods:
        report: "list[dict]" = []
        t0 = time.time()
        try:
            _run_method(m, rows, forbidden=None, report=report,
                        cutoff=args.cutoff)
        except ValueError as e:
            print(f"  {m:<14} {e}")
            continue
        pairs = {frozenset({(r["a_db"], r["a"]), (r["b_db"], r["b"])})
                 for r in report}
        viol = sum(1 for r in report if (r["a"], r["b"]) in forbidden)
        found = len(pairs & identical)
        rate = (100.0 * viol / len(report)) if report else 0.0
        rec = (100.0 * found / len(identical)) if identical else 0.0
        results[m] = {"pairs": pairs, "merges": len(report),
                       "violations": viol, "found": found,
                       "seconds": round(time.time() - t0, 1)}
        print(f"  {m:<14}{len(report):>9,}{viol:>11,}{rate:>10.1f}%"
              f"{found:>13,}{rec:>8.1f}%")

    if len(results) > 1 and args.agreement:
        print("\n  pairwise agreement (Jaccard over proposed merges):")
        ms = list(results)
        for i, a in enumerate(ms):
            for b in ms[i + 1:]:
                pa, pb = results[a]["pairs"], results[b]["pairs"]
                u = len(pa | pb)
                print(f"    {a:<13} vs {b:<13} "
                      f"{(len(pa & pb)/u if u else 0):.3f}  "
                      f"({len(pa & pb):,} shared)")

    if args.out:
        out = Path(args.out)
        if not out.is_absolute():
            out = ROOT / args.out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(
            {m: {k: v for k, v in r.items() if k != "pairs"}
             for m, r in results.items()}, indent=2))
        print(f"\n  wrote {out}")
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
    b.add_argument("--merge", default="consensus",
                   help="Pathway unification criterion, or several joined by "
                        "'+' for their union. One of: " + ", ".join(METHODS) +
                        ". Default consensus (IntPath's LCS rule as the "
                        "candidate generator, confirmed by independent "
                        "signals); 'lcs' is the published rule alone.")
    b.add_argument("--extra-gmt", action="append", metavar="LABEL=PATH",
                   help="Integrate a local GMT as another source database, "
                        "e.g. --extra-gmt HumanCyc=~/humancyc.gmt. Repeatable. "
                        "The route for databases that cannot be fetched "
                        "anonymously, such as BioCyc/HumanCyc.")
    b.add_argument("--no-hierarchy", action="store_true",
                   help="Do not read Reactome's own parent/child hierarchy "
                        "when deciding which pathways may be unified")
    b.add_argument("--within-db", action="store_true",
                   help="Also unify names inside one database, as IntPath "
                        "does. Off by default because Reactome's hierarchy "
                        "means parent and child pathways to stay distinct")
    b.add_argument("--no-prune", action="store_true",
                   help="Keep pathway nodes from earlier pathwaydb builds "
                        "that this build no longer produces")
    b.add_argument("--no-pull", action="store_true",
                   help="Fail rather than downloading if the cache is empty")

    _common(sub.add_parser("status", help="What is cached and how old"),
            with_sources=False)
    sub.add_parser("sources", help="Databases, licences and coverage")

    e = _common(sub.add_parser("evaluate",
                                help="Score each unification criterion"))
    e.add_argument("--methods", default="lcs,tokens,jaccard,containment,"
                                        "hypergeom",
                   help="Comma-separated criteria to score")
    e.add_argument("--cutoff", type=float,
                   help="Override the criterion's cutoff")
    e.add_argument("--agreement", action="store_true",
                   help="Also print pairwise agreement between criteria")
    e.add_argument("--out", help="Write the scores to a JSON file")

    q = _common(sub.add_parser("query", help="Pathways for a gene list"))
    q.add_argument("--genes", required=True)
    q.add_argument("--top", type=int, default=25)
    q.add_argument("--merge", default="consensus",
                   help="Pathway unification criterion (default consensus)")
    q.add_argument("--verbose", action="store_true",
                   help="Also list which query genes hit each pathway")

    args = p.parse_args(argv)
    return {"pull": cmd_pull, "build": cmd_build, "status": cmd_status,
            "query": cmd_query, "sources": cmd_sources,
            "evaluate": cmd_evaluate}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
