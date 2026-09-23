#!/usr/bin/env python3
"""Seqspec catalogue, validator and read-format indexer for CRISPR assays (port of IGVF-CRISPR/CRISPR-SeqSpec).

Port of https://github.com/IGVF-CRISPR/CRISPR-SeqSpec (No LICENSE file; the
repository declares no licence), the IGVF CRISPR focus group's collection of
seqspec (https://github.com/pachterlab/seqspec) descriptions of CRISPR
single-cell assays: a CROP-seq + 10x v3 + MULTI-seq experiment split into one
YAML per modality (rna.yml, guide.yml, multiseq.yml, seqspec 0.2.0) and the
TAP-seq screen of Schraivogel et al. 2020 (spec.yaml, seqspec 0.0.3), each
next to 1,000-read FASTQs and the barcode / guide onlists.  The upstream
repository holds data, not code, so the skill re-implements -- in plain Python,
no seqspec, PyYAML or jsonschema needed -- what the upstream README asks of a
submission and what `seqspec check`, `seqspec index`, `seqspec onlist` and
`seqspec print` do with such a file.  Relationship: port (definitions taken
from the seqspec specification the files are written against; no code was
copied).  Because the upstream has no licence, the module embeds only a
compact derived JSON index of the specs (region trees with run-length encoded
sequences) and fetches the original files from the pinned commit on demand.

Definitions reproduced
  YAML          seqspec's tagged YAML (!Assay / !Region / !Read / !Onlist) is
                parsed by a small block-YAML reader (PyYAML is used when
                installed); tags are kept as '__tag__'.
  schemas       seqspec 0.2.x: Assay{seqspec_version, assay_id, name, doi, date,
                description, modalities, lib_struct, sequence_protocol,
                sequence_kit, library_protocol, library_kit, sequence_spec,
                library_spec}; Read{read_id, name, modality, primer_id,
                min_len, max_len, strand}; Region{parent_id, region_id,
                region_type, name, sequence_type, sequence, min_len, max_len,
                onlist, regions}; Onlist{location, filename, md5}.
                seqspec 0.0.x: Assay{assay, sequencer, seqspec_version, name,
                doi, publication_date, description, modalities, lib_struct,
                assay_spec}; reads are the FASTQ-named child regions of each
                modality region (their children are laid out 5'->3' in read
                coordinates); Onlist{filename, md5}.
  enums         modality in {dna, rna, tag, protein, atac, crispr};
                sequence_type in {fixed, random, onlist, joined}; strand in
                {pos, neg}; region_type in the seqspec 0.2 vocabulary
                (barcode, umi, cdna, gdna, sgrna_target, poly_A, linker,
                truseq_read1/2, nextera_read1/2, illumina_p5/p7, index5/7,
                custom_primer, ...); onlist location in {local, remote}.
  check         C01 required keys and types, C02 enums, C03 every modality has
                a library region whose region_id is the modality, C04 unique
                region_ids within a modality, C05 unique read_ids, C06 a read's
                primer_id exists in its modality's library, C07 a read's
                modality is declared, C08 onlist files exist (spec dir; name or
                name +/- '.gz'), C09 md5 is 32 hex chars and matches the file,
                C10 min_len <= len(sequence) <= max_len for leaf regions,
                C11 alphabet: fixed = ACGT only, random = all N or all X,
                onlist = all N, C12 min_len <= max_len, C13 a parent's min/max
                length equals the sum of its children's, C14 read max_len <=
                sequence available after the primer (sum of downstream
                max_len), C15 read FASTQ files exist, C16 parent_id matches the
                enclosing region (0.2.x).  Errors fail the check; warnings do
                not.
  projection    (seqspec `index`) leaves of the modality region in 5'->3'
                order; strand pos reads take the leaves after the primer,
                strand neg reads the leaves before it in reverse; each region
                occupies [start, start + max_len) and the list is cut at the
                read's max_len (the last region is truncated).  0.0.x reads
                lay their children out from 0 the same way.
  kb            "<barcodes>:<umis>:<features>" where each group is a comma
                list of "<fastq index>,<start>,<stop>" (0-based, stop
                exclusive), "-1,-1,-1" when a group is empty; barcodes are
                region types containing 'barcode', UMIs 'umi', features cdna /
                gdna (or, for crispr/guide modalities, sgrna / sgrna_target /
                crispr regions when present, else cdna).
  chromap       "--read-format bc:<start>:<stop-1>,r1:<start>:<stop-1>[,r2:...]"
                (inclusive ends) plus -1/-2/--barcode FASTQ arguments; barcodes
                must come from one FASTQ and features from at most two.
  starsolo      --soloType CB_UMI_Simple --soloCBstart/--soloCBlen/
                --soloUMIstart/--soloUMIlen (1-based starts) from the barcode
                read.
  onlist        every region with sequence_type onlist -> (modality, region,
                filename, location, md5, local path if present).
  check-reads   the README's submission rule (a ~1M-read, or here any, FASTQ
                pair per modality): FASTQs exist, read lengths inside
                [min_len, max_len], fixed regions match at their projected
                coordinates (Hamming <= 1), barcode / guide slices found in
                their onlists.

Subcommands
  catalog      every spec in the upstream repo (embedded index, a local clone
               with --from-dir, or the pinned commit with --live): assay,
               version, modalities, reads, regions, onlists, kb strings, check
               result -> catalog.tsv / reads.tsv / regions.tsv / onlists.tsv /
               catalog.json.
  check        validate a seqspec (path or catalogue name) -> issues.tsv.
  index        kb / chromap / starsolo / tab read-format string for a modality.
  onlist       onlist filenames (and local paths) of a spec.
  info         print the library structure of a spec (seqspec print-like).
  check-reads  validate FASTQs against a spec (lengths, fixed regions, onlists).
  fetch        download the upstream files at the pinned commit into
               Data/CRISPRSeqSpec (about 20 MB; --specs-only for the YAMLs).
  selftest     synthetic 0.2.x and 0.0.x specs + FASTQs with planted defects.

Output: Docs/CRISPRSeqSpec/<timestamp>_<label>/ (report.md, summary.json and
tables).  Pure standard library + pandas; matplotlib is not used.

Usage:
    igvfagent crispr-seqspec catalog
    igvfagent crispr-seqspec check --spec multiseq_10xv3_cropseq/guide.yml
    igvfagent crispr-seqspec index --spec TAPseq_Schraivogel2020/spec.yaml --modality crispr --tool kb
    igvfagent crispr-seqspec onlist --spec rna
    igvfagent crispr-seqspec info --spec tapseq
    igvfagent crispr-seqspec fetch
    igvfagent crispr-seqspec check-reads --spec Data/CRISPRSeqSpec/multiseq_10xv3_cropseq/guide.yml
    igvfagent crispr-seqspec selftest
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRSeqSpec"
DATA_ROOT = ROOT / "Data" / "CRISPRSeqSpec"

UPSTREAM_REPO = "IGVF-CRISPR/CRISPR-SeqSpec"
UPSTREAM_COMMIT = "890435311b0c36859c97f6afbceb3c61910f63b2"
UPSTREAM_DATE = "2024-07-12T16:49:26Z"
RAW = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/"

# Files at the pinned commit: path -> (bytes, md5).
REPO_FILES = {
    "README.md": (237, "86a082a50c8cab81dec5867b062c6a63"),
    "TAPseq_Schraivogel2020/10x_bc_whitelist_737k_201608.txt.gz": (2237444, "33438767d3070dcbdf27c1f34cef3cf6"),
    "TAPseq_Schraivogel2020/R1.fastq.gz": (531046, "1c421e8315fa930dd68911cfb8ce4bf7"),
    "TAPseq_Schraivogel2020/R2.fastq.gz": (795505, "02cbdde6ea9940524e5e0a6c45edf415"),
    "TAPseq_Schraivogel2020/guides_schraivogel2020.txt.gz": (43679, "1fa3e81b3110a599c7fdd549019ac1d4"),
    "TAPseq_Schraivogel2020/spec.yaml": (3621, "f5deb88a9d4580a171abe003ce3d3a7b"),
    "multiseq_10xv3_cropseq/3M-february-2018.txt.gz": (12211647, "417292fed237e038f4d7e20b0ddac53e"),
    "multiseq_10xv3_cropseq/graphical_repr.png": (72031, "305e6eb6dcd63681f0fe6bfd8d11be42"),
    "multiseq_10xv3_cropseq/guide.yml": (2343, "7b638fcfc35807d4088221d643d7ead5"),
    "multiseq_10xv3_cropseq/guide_R1.fq": (551533, "7b6e007a55eb53a597558b0cd707ba27"),
    "multiseq_10xv3_cropseq/guide_R2.fq": (551533, "3813c83654729d251c67f57e1d8e57d3"),
    "multiseq_10xv3_cropseq/multiseq.yml": (2074, "d1e79372c338e641aebeeb39201c364b"),
    "multiseq_10xv3_cropseq/multiseq_R1.fq": (552688, "2f86d5b57e04a8739f68487f41f7bffa"),
    "multiseq_10xv3_cropseq/multiseq_R2.fq": (552688, "a3fd59696ed4a8349f2c983bd7520b56"),
    "multiseq_10xv3_cropseq/readme.md": (1591, "60edceb0aa2f102466001abf7fb8428b"),
    "multiseq_10xv3_cropseq/rna.yml": (2221, "42cff7c17d65cfa776370ce022c08b75"),
    "multiseq_10xv3_cropseq/scRNA_R1.fq": (551559, "44d4c5179eafa4e0b820a758534edc4a"),
    "multiseq_10xv3_cropseq/scRNA_R2.fq": (551559, "c784375089042d26bc13fe15b2eff26b"),
}
SPEC_FILES = [p for p in REPO_FILES if p.endswith((".yml", ".yaml"))]
SPEC_ALIASES = {"rna": "multiseq_10xv3_cropseq/rna.yml", "guide": "multiseq_10xv3_cropseq/guide.yml",
                "multiseq": "multiseq_10xv3_cropseq/multiseq.yml", "tapseq": "TAPseq_Schraivogel2020/spec.yaml",
                "tap-seq": "TAPseq_Schraivogel2020/spec.yaml", "cropseq": "multiseq_10xv3_cropseq/guide.yml"}

# Compact derived index of the upstream specs (generated by `_compact()` from the
# pinned commit).  Region = [region_id, region_type, sequence_type, sequence (RLE),
# min_len, max_len, onlist or null, children]; onlist = [filename, location, md5].
CATALOG_INDEX_JSON = r"""[{"path":"TAPseq_Schraivogel2020/spec.yaml","sha256":"9fc911955c4e1f94b0e9e22460dca64d5ef6e62b5bf732d7590c088298093779","legacy":true,"meta":{"assay":"TAP-seq (CROP-seq vector)","sequencer":"Illumina","seqspec_version":"0.0.3","name":"TAPseq_Schraivogel2020","doi":"https://doi.org/10.1038/s41592-020-0837-5","publication_date":"01 June 2020","description":"TAP-seq enhancer screens from Schraivogel et al., 2020 using the CROP-seq vector for guide delivery","modalities":["rna","crispr"],"lib_struct":"https://doi.org/10.21203/rs.3.pex-864/v1"},"reads":[],"library":[["rna","","","",0,1024,null,[["R1.fastq.gz","","","",0,1024,null,[["barcode","","onlist","N16",16,16,["10x_bc_whitelist_737k_201608.txt.gz",null,null],[]],["umi","","random","N10",10,10,null,[]]]],["R2.fastq.gz","","","",0,1024,null,[["cdna","","random","",58,58,null,[]]]]]],["crispr","","","",0,1024,null,[["R1.fastq.gz","","","",0,1024,null,[["barcode","","onlist","N16",16,16,["10x_bc_whitelist_737k_201608.txt.gz",null,null],[]],["umi","","random","N10",10,10,null,[]]]],["R2.fastq.gz","","","",0,1024,null,[["cdna","","random","",18,18,null,[]],["sgrna","","onlist","N20",20,20,["guides_schraivogel2020.txt.gz",null,null],[]],["tapseq_primer","","fixed","TGTG2A3G2ACGA3CAC2",20,20,null,[]]]]]]]},{"path":"multiseq_10xv3_cropseq/guide.yml","sha256":"c52797e3232b5453cc024bb202f642e5301e6d302e09a6a9a58ce8271d897a1a","legacy":false,"meta":{"seqspec_version":"0.2.0","assay_id":"CROPSEQ_10XV3_MULTISEQ","name":"CROPSEQ_10XV3_MULTISEQ","doi":"XXX","date":"08 July 2022","description":"CROP-seq USING 10XV3 and Multiseq","modalities":["guide"],"lib_struct":"https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html","sequence_protocol":"Not-specified","sequence_kit":"Not-specified","library_protocol":"10xv3 RNA + crop-seq + multiseq","library_kit":"Not-specified"},"reads":[["guide_R1.fq","Read 1","guide","r1_primer",26,26,"pos"],["guide_R2.fq","Read 2","guide","r2_primer",43,43,"neg"]],"library":[["guide",null,null,"A16N28XA16",150,150,null,[["r1_primer","r1_primer","fixed","A16",16,16,null,[]],["barcode","barcode","onlist","N16",16,16,["3M-february-2018.txt","remote","updateit"],[]],["umi","umi","fixed","N12",12,12,null,[]],["cdna","cdna","random","N20",20,20,null,[]],["common","common","fixed","CT2GTG2A3G2ACGA3CAC2G",23,23,null,[]],["r2_primer","r2_primer","fixed","A16",16,16,null,[]]]]]},{"path":"multiseq_10xv3_cropseq/multiseq.yml","sha256":"ab24a469e4061fcbc7c74a06768ac7dd723c53ffa83a70d34b37930dcc7e5148","legacy":false,"meta":{"seqspec_version":"0.2.0","assay_id":"10xv3","name":"10xv3","doi":"","date":"9 July 2024","description":"Multi-seq","modalities":["multiseq"],"lib_struct":"https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html","sequence_protocol":"Not-specified","sequence_kit":"Not-specified","library_protocol":"","library_kit":"Not-specified"},"reads":[["multiseq_R1.fq","Read 1","multiseq","r1_primer",26,26,"pos"],["multiseq_R2.fq","Read 2","multiseq","r2_primer",150,150,"neg"]],"library":[["multiseq",null,null,"A16N28XA16",150,150,null,[["r1_primer","r1_primer","fixed","A16",16,16,null,[]],["cell_barcode","cell_barcode","onlist","N16",16,16,["3M-february-2018.txt","remote","updateit"],[]],["umi","umi","fixed","N12",12,12,null,[]],["cdna","cdna","random","N8",8,8,null,[]],["r2_primer","r2_primer","fixed","A16",16,16,null,[]]]]]},{"path":"multiseq_10xv3_cropseq/rna.yml","sha256":"4e7735bce2c9e0a329d987c6ac6f5797ec0cda70d6d42f2d67f5dd2fbb3652b8","legacy":false,"meta":{"seqspec_version":"0.2.0","assay_id":"10xv3","name":"10xv3","doi":"https://doi.org/10.1126/science.aam8999","date":"15 March 2018","description":"10x Genomics v3 single-cell rnaseq","modalities":["rna"],"lib_struct":"https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html","sequence_protocol":"Not-specified","sequence_kit":"Not-specified","library_protocol":"10xv3 RNA","library_kit":"Not-specified"},"reads":[["scRNA_R1.fq","Read 1","rna","r1_primer",28,28,"pos"],["scRNA_R2.fq","Read 2","rna","r2_primer",150,150,"neg"]],"library":[["rna",null,null,"A16N28XA16",150,150,null,[["r1_primer","r1_primer","fixed","A16",16,16,null,[]],["barcode","barcode","onlist","N16",16,16,["3M-february-2018.txt","remote","updateit!"],[]],["umi","umi","fixed","N12",12,12,null,[]],["cdna","cdna","random","X150",1,150,null,[]],["r2_primer","r2_primer","fixed","A16",16,16,null,[]]]]]}]"""

MODALITIES = {"dna", "rna", "tag", "protein", "atac", "crispr"}
SEQUENCE_TYPES = {"fixed", "random", "onlist", "joined"}
STRANDS = {"pos", "neg"}
ONLIST_LOCATIONS = {"local", "remote"}
REGION_TYPES = {
    "atac", "barcode", "cdna", "crispr", "custom_primer", "dna", "fastq", "fastq_link", "gdna", "hic",
    "illumina_p5", "illumina_p7", "index5", "index7", "linker", "ME1", "ME2", "methyl", "named",
    "nextera_read1", "nextera_read2", "poly_A", "poly_G", "poly_T", "poly_C", "protein", "rna", "s5", "s7",
    "sgrna_target", "tag", "truseq_read1", "truseq_read2", "umi",
}
ASSAY_KEYS_V02 = ["seqspec_version", "assay_id", "name", "doi", "date", "description", "modalities", "lib_struct",
                  "sequence_protocol", "sequence_kit", "library_protocol", "library_kit", "sequence_spec", "library_spec"]
ASSAY_KEYS_V00 = ["assay", "sequencer", "seqspec_version", "name", "doi", "publication_date", "description",
                  "modalities", "lib_struct", "assay_spec"]
READ_KEYS = {"read_id": str, "name": str, "modality": str, "primer_id": str, "min_len": int, "max_len": int, "strand": str}
REGION_KEYS_V02 = ["parent_id", "region_id", "region_type", "name", "sequence_type", "sequence", "min_len", "max_len",
                   "onlist", "regions"]
REGION_KEYS_V00 = ["region_id", "region_type", "name", "sequence_type", "sequence", "min_len", "max_len", "onlist", "regions"]
FEATURE_TYPES = {"cdna", "gdna", "dna", "rna"}
GUIDE_FEATURE_TYPES = {"sgrna", "sgrna_target", "crispr", "guide"}
FASTQ_RE = re.compile(r"\.(fq|fastq)(\.gz)?$", re.I)

log = logging.getLogger("crispr_seqspec")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_seqspec_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_json(path: Path, obj: Any) -> Path:
    path.write_text(json.dumps(obj, indent=2, default=str))
    print(f"JSON: {path}")
    return path


def write_tsv(path: Path, rows: List[Dict[str, Any]], columns: Optional[List[str]] = None) -> Path:
    cols = columns or (list(rows[0].keys()) if rows else [])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"TSV: {path}")
    return path


def write_report(path: Path, lines: List[str]) -> Path:
    path.write_text("\n".join(lines) + "\n")
    print(f"Report: {path}")
    return path


# ---------------------------------------------------------------------------
# Minimal tagged block-YAML reader (the subset seqspec files use)
# ---------------------------------------------------------------------------

_INT_RE = re.compile(r"^[-+]?\d+$")
_FLOAT_RE = re.compile(r"^[-+]?(\d+\.\d*|\.\d+|\d+)([eE][-+]?\d+)?$")


def _scalar(tok: str) -> Any:
    t = tok.strip()
    if t == "" or t in ("null", "Null", "NULL", "~"):
        return None
    if t[0] == "'" and t.endswith("'") and len(t) >= 2:
        return t[1:-1].replace("''", "'")
    if t[0] == '"' and t.endswith('"') and len(t) >= 2:
        return bytes(t[1:-1], "utf-8").decode("unicode_escape")
    if t in ("true", "True", "TRUE"):
        return True
    if t in ("false", "False", "FALSE"):
        return False
    if _INT_RE.match(t):
        return int(t)
    if _FLOAT_RE.match(t) and "." in t:
        return float(t)
    if t.startswith("[") and t.endswith("]"):
        inner = t[1:-1].strip()
        return [] if not inner else [_scalar(x) for x in inner.split(",")]
    return t


def _strip_comment(line: str) -> str:
    out, q = [], None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
            break
        out.append(ch)
    return "".join(out).rstrip()


class _Lines:
    def __init__(self, text: str):
        self.items: List[Tuple[int, str]] = []
        for raw in text.replace("\t", "  ").splitlines():
            s = _strip_comment(raw)
            if not s.strip() or s.strip() in ("---", "..."):
                continue
            self.items.append((len(s) - len(s.lstrip(" ")), s.strip()))


def _split_tag(rest: str) -> Tuple[Optional[str], str]:
    rest = rest.strip()
    if rest.startswith("!"):
        parts = rest.split(None, 1)
        return parts[0][1:], (parts[1] if len(parts) > 1 else "")
    return None, rest


def _split_key(s: str) -> Optional[Tuple[str, str]]:
    m = re.match(r"^('[^']*'|\"[^\"]*\"|[^:'\"]+?)\s*:(\s+|$)(.*)$", s)
    if not m:
        return None
    key = m.group(1)
    if key[0] in "'\"":
        key = key[1:-1]
    return key.strip(), m.group(3)


def _parse_block(L: _Lines, i: int, indent: int) -> Tuple[Any, int]:
    if i >= len(L.items):
        return None, i
    ind, s = L.items[i]
    if s.startswith("- ") or s == "-":
        return _parse_seq(L, i, ind)
    return _parse_map(L, i, ind)


def _parse_value(L: _Lines, i: int, ind: int, rest: str) -> Tuple[Any, int]:
    """Value after 'key:' or '- ' on line i-1 (already consumed); rest is the inline text."""
    tag, rest = _split_tag(rest)
    if rest == "":
        if i < len(L.items):
            nind, ns = L.items[i]
            if nind > ind or (nind == ind and (ns.startswith("- ") or ns == "-") and tag is None):
                val, i = _parse_block(L, i, nind)
                if tag and isinstance(val, dict):
                    val["__tag__"] = tag
                return val, i
        return ({"__tag__": tag} if tag else None), i
    return _scalar(rest), i


def _parse_map(L: _Lines, i: int, ind: int) -> Tuple[Dict[str, Any], int]:
    out: Dict[str, Any] = {}
    while i < len(L.items):
        cind, s = L.items[i]
        if cind != ind or s.startswith("- ") or s == "-":
            if cind < ind or (cind == ind and (s.startswith("- ") or s == "-")):
                break
            raise ValueError(f"bad indentation near: {s!r}")
        kv = _split_key(s)
        if kv is None:
            raise ValueError(f"expected 'key: value', got {s!r}")
        k, rest = kv
        out[k], i = _parse_value(L, i + 1, ind, rest)
    return out, i


def _parse_seq(L: _Lines, i: int, ind: int) -> Tuple[List[Any], int]:
    out: List[Any] = []
    while i < len(L.items):
        cind, s = L.items[i]
        if cind != ind or not (s.startswith("- ") or s == "-"):
            break
        body = s[1:].strip()
        tag, body = _split_tag(body)
        if body and _split_key(body):
            # inline first key of a mapping item: re-read it at the item's inner indent
            inner = ind + 2
            L.items[i] = (inner, body)
            val, i = _parse_map(L, i, inner)
            if tag:
                val["__tag__"] = tag
            out.append(val)
        elif body == "":
            if i + 1 < len(L.items) and L.items[i + 1][0] > ind:
                val, i = _parse_block(L, i + 1, L.items[i + 1][0])
                if tag and isinstance(val, dict):
                    val["__tag__"] = tag
                out.append(val)
            else:
                out.append({"__tag__": tag} if tag else None)
                i += 1
        else:
            out.append(_scalar(body))
            i += 1
    return out, i


def parse_yaml(text: str) -> Any:
    """Parse seqspec YAML; PyYAML (with tag passthrough) when importable, else the built-in reader."""
    try:
        import yaml  # type: ignore

        class _L(yaml.SafeLoader):
            pass

        def _tagged(loader, suffix, node):
            d = loader.construct_mapping(node, deep=True)
            d["__tag__"] = suffix
            return d
        _L.add_multi_constructor("!", _tagged)
        doc = yaml.load(text, Loader=_L)
        # PyYAML resolves unquoted dates such as 2020-06-01; keep them as strings like the fallback
        return doc
    except ImportError:
        pass
    L = _Lines(text)
    if not L.items:
        return None
    ind, s = L.items[0]
    tag = None
    if s.startswith("!") and len(s.split()) == 1:
        tag = s[1:]
        L.items.pop(0)
    doc, _ = _parse_block(L, 0, L.items[0][0] if L.items else 0)
    if tag and isinstance(doc, dict):
        doc["__tag__"] = tag
    return doc


# ---------------------------------------------------------------------------
# Spec model
# ---------------------------------------------------------------------------

def _rle(seq: Optional[str]) -> str:
    if not seq:
        return ""
    out, prev, n = [], seq[0], 0
    for ch in seq:
        if ch == prev:
            n += 1
        else:
            out.append(prev + (str(n) if n > 1 else ""))
            prev, n = ch, 1
    out.append(prev + (str(n) if n > 1 else ""))
    return "".join(out)


def _unrle(s: str) -> str:
    return "".join(ch * (int(n) if n else 1) for ch, n in re.findall(r"([A-Za-z])(\d*)", s or ""))


def spec_version(doc: Dict[str, Any]) -> str:
    return str(doc.get("seqspec_version") or "")


def is_legacy(doc: Dict[str, Any]) -> bool:
    """seqspec 0.0.x (assay_spec with FASTQ-named read regions) vs 0.1+/0.2 (sequence_spec + library_spec)."""
    return "assay_spec" in doc and "library_spec" not in doc


def rtype(region: Dict[str, Any]) -> str:
    t = region.get("region_type")
    return str(t) if t else str(region.get("region_id") or "")


def children(region: Dict[str, Any]) -> List[Dict[str, Any]]:
    r = region.get("regions")
    return [c for c in r if isinstance(c, dict)] if isinstance(r, list) else []


def leaves(region: Dict[str, Any]) -> List[Dict[str, Any]]:
    ch = children(region)
    if not ch:
        return [region]
    out: List[Dict[str, Any]] = []
    for c in ch:
        out.extend(leaves(c))
    return out


def walk(region: Dict[str, Any], depth: int = 0, parent: Optional[Dict[str, Any]] = None):
    yield region, depth, parent
    for c in children(region):
        yield from walk(c, depth + 1, region)


def library_regions(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    key = "assay_spec" if is_legacy(doc) else "library_spec"
    lib = doc.get(key) or []
    return [r for r in lib if isinstance(r, dict)]


def modality_region(doc: Dict[str, Any], modality: str) -> Optional[Dict[str, Any]]:
    for r in library_regions(doc):
        if r.get("region_id") == modality:
            return r
    return None


def modalities(doc: Dict[str, Any]) -> List[str]:
    m = doc.get("modalities") or []
    return [str(x) for x in m] if isinstance(m, list) else [str(m)]


def reads(doc: Dict[str, Any], modality: Optional[str] = None) -> List[Dict[str, Any]]:
    """Reads of a spec; for 0.0.x the FASTQ-named child regions of each modality region become reads."""
    out: List[Dict[str, Any]] = []
    if is_legacy(doc):
        for mreg in library_regions(doc):
            mod = str(mreg.get("region_id"))
            for c in children(mreg):
                rid = str(c.get("region_id") or "")
                if FASTQ_RE.search(rid) or rtype(c) in ("fastq", "fastq_link") or children(c):
                    total = sum(int(x.get("max_len") or 0) for x in children(c))
                    tmin = sum(int(x.get("min_len") or 0) for x in children(c))
                    out.append({"read_id": rid, "name": c.get("name") or rid, "modality": mod, "primer_id": None,
                                "min_len": tmin, "max_len": total, "strand": "pos", "_region": c, "_legacy": True})
    else:
        for r in doc.get("sequence_spec") or []:
            if isinstance(r, dict):
                out.append(dict(r))
    if modality is not None:
        out = [r for r in out if str(r.get("modality")) == modality]
    return out


def project_read(doc: Dict[str, Any], read: Dict[str, Any]) -> List[Dict[str, Any]]:
    """seqspec `index`: regions covered by a read with [start, stop) read coordinates, cut at read max_len."""
    if read.get("_legacy"):
        seq = leaves(read["_region"]) if children(read["_region"]) else []
    else:
        mreg = modality_region(doc, str(read.get("modality")))
        if mreg is None:
            return []
        lv = leaves(mreg)
        ids = [str(r.get("region_id")) for r in lv]
        pid = str(read.get("primer_id"))
        if pid in ids:
            i0 = ids.index(pid)
            lo, hi = i0, i0
        else:
            sub = next((r for r, _, _ in walk(mreg) if str(r.get("region_id")) == pid), None)
            if sub is None:
                return []
            sl = [str(x.get("region_id")) for x in leaves(sub)]
            lo, hi = ids.index(sl[0]), ids.index(sl[-1])
        seq = lv[hi + 1:] if str(read.get("strand", "pos")) == "pos" else list(reversed(lv[:lo]))
    out, pos = [], 0
    rmax = int(read.get("max_len") or 0)
    for r in seq:
        if pos >= rmax:
            break
        ln = int(r.get("max_len") or 0)
        stop = min(pos + ln, rmax)
        out.append({"region_id": r.get("region_id"), "region_type": rtype(r), "sequence_type": r.get("sequence_type"),
                    "sequence": r.get("sequence") or "", "start": pos, "stop": stop, "truncated": stop < pos + ln,
                    "onlist": r.get("onlist"), "min_len": r.get("min_len"), "max_len": ln})
        pos = stop
    return out


def _is_barcode(t: str) -> bool:
    return "barcode" in t.lower()


def _feature_types(modality: str, regions: List[Dict[str, Any]], override: Optional[List[str]] = None) -> set:
    if override:
        return {x.lower() for x in override}
    types = {str(r["region_type"]).lower() for r in regions}
    if modality.lower() in ("crispr", "guide") and types & GUIDE_FEATURE_TYPES:
        return GUIDE_FEATURE_TYPES
    return FEATURE_TYPES


def index_modality(doc: Dict[str, Any], modality: str, dedupe: bool = False) -> List[Tuple[Dict[str, Any], List[Dict[str, Any]]]]:
    out, seen = [], set()
    for rd in reads(doc, modality):
        regs = project_read(doc, rd)
        if dedupe:
            regs = [r for r in regs if r["region_id"] not in seen]
        seen.update(r["region_id"] for r in regs)
        out.append((rd, regs))
    return out


def format_index(doc: Dict[str, Any], modality: str, tool: str = "kb", dedupe: bool = False,
                 feature_types: Optional[List[str]] = None) -> Dict[str, Any]:
    idx = index_modality(doc, modality, dedupe=dedupe)
    allregs = [r for _, regs in idx for r in regs]
    ftypes = _feature_types(modality, allregs, feature_types)
    bcs, umis, feats = [], [], []
    for i, (rd, regs) in enumerate(idx):
        for r in regs:
            t = str(r["region_type"]).lower()
            if _is_barcode(t):
                bcs.append((i, r))
            elif t == "umi":
                umis.append((i, r))
            elif t in ftypes:
                feats.append((i, r))
    files = [str(rd.get("read_id")) for rd, _ in idx]
    dup = [rid for rid in {r["region_id"] for r in allregs} if sum(1 for x in allregs if x["region_id"] == rid) > 1]
    notes = []
    if dup:
        notes.append(f"regions covered by more than one read (read longer than the insert): {', '.join(sorted(map(str, dup)))}"
                     " -- use --dedupe to keep each region in the first read that covers it")
    res: Dict[str, Any] = {"modality": modality, "tool": tool, "fastqs": files, "notes": notes}
    if tool == "kb":
        g = lambda lst: ",".join(f"{i},{r['start']},{r['stop']}" for i, r in lst) or "-1,-1,-1"
        res["string"] = f"{g(bcs)}:{g(umis)}:{g(feats) if feats else '-1,-1,-1'}"
    elif tool == "chromap":
        bc_files = sorted({i for i, _ in bcs})
        f_files = []
        for i, _ in feats:
            if i not in f_files:
                f_files.append(i)
        if len(bc_files) > 1:
            raise ValueError("chromap supports barcodes from one FASTQ only")
        if len(f_files) > 2:
            raise ValueError("chromap supports genomic reads from at most two FASTQs")
        parts = [f"bc:{r['start']}:{r['stop'] - 1}" for _, r in bcs]
        for n, fi in enumerate(f_files, 1):
            rr = [r for i, r in feats if i == fi]
            parts.append(f"r{n}:{rr[0]['start']}:{rr[-1]['stop'] - 1}")
        args = []
        for n, fi in enumerate(f_files, 1):
            args += [f"-{n}", files[fi]]
        if bc_files:
            args += ["--barcode", files[bc_files[0]]]
        res["string"] = ",".join(parts)
        res["command_args"] = " ".join(args + ["--read-format", res["string"]])
    elif tool == "starsolo":
        if not bcs:
            raise ValueError("no barcode region: STARsolo CB_UMI_Simple needs one")
        bi, b = bcs[0]
        s = f"--soloType CB_UMI_Simple --soloCBstart {b['start'] + 1} --soloCBlen {b['stop'] - b['start']}"
        u = [r for i, r in umis if i == bi]
        if u:
            s += f" --soloUMIstart {u[0]['start'] + 1} --soloUMIlen {u[0]['stop'] - u[0]['start']}"
        cdna = [files[i] for i, _ in feats]
        res["string"] = s
        res["command_args"] = f"--readFilesIn {' '.join(dict.fromkeys(cdna))} {files[bi]} {s}"
    elif tool == "tab":
        res["string"] = "\n".join(f"{files[i]}\t{r['region_id']}\t{r['region_type']}\t{r['start']}\t{r['stop']}"
                                  for i, (_, regs) in enumerate(idx) for r in regs)
    else:
        raise ValueError(f"unknown tool {tool!r} (kb, chromap, starsolo, tab)")
    res["regions"] = [{"fastq": files[i], "fastq_index": i, **{k: r[k] for k in ("region_id", "region_type", "start", "stop", "truncated")}}
                      for i, (_, regs) in enumerate(idx) for r in regs]
    return res


def onlists(doc: Dict[str, Any], spec_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    out = []
    for mreg in library_regions(doc):
        for r, depth, parent in walk(mreg):
            ol = r.get("onlist")
            if isinstance(ol, dict) or (r.get("sequence_type") == "onlist" and ol is not None):
                ol = ol if isinstance(ol, dict) else {}
                fn = ol.get("filename")
                local = _resolve_file(spec_dir, fn) if (spec_dir and fn) else None
                out.append({"modality": mreg.get("region_id"), "region_id": r.get("region_id"), "region_type": rtype(r),
                            "filename": fn, "location": ol.get("location"), "md5": ol.get("md5"),
                            "local_path": str(local) if local else None})
    return out


def _resolve_file(spec_dir: Optional[Path], name: Optional[str]) -> Optional[Path]:
    if not spec_dir or not name:
        return None
    cands = [spec_dir / name]
    cands.append(spec_dir / (name[:-3] if name.endswith(".gz") else name + ".gz"))
    for c in cands:
        if c.exists():
            return c
    return None


# ---------------------------------------------------------------------------
# check (seqspec check equivalent)
# ---------------------------------------------------------------------------

def _issue(out: List[Dict[str, Any]], level: str, code: str, where: str, msg: str) -> None:
    out.append({"level": level, "code": code, "where": where, "message": msg})


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_spec(doc: Any, spec_dir: Optional[Path] = None, check_md5: bool = True) -> List[Dict[str, Any]]:
    iss: List[Dict[str, Any]] = []
    if not isinstance(doc, dict):
        _issue(iss, "error", "C01", "assay", "document is not a mapping (!Assay)")
        return iss
    legacy = is_legacy(doc)
    keys = ASSAY_KEYS_V00 if legacy else ASSAY_KEYS_V02
    ver = spec_version(doc)
    if doc.get("__tag__") not in (None, "Assay"):
        _issue(iss, "error", "C01", "assay", f"root tag is !{doc.get('__tag__')}, expected !Assay")
    for k in keys:
        if k not in doc:
            _issue(iss, "error", "C01", "assay", f"missing required key '{k}' (seqspec {ver or '?'})")
    for k in ("name", "description", "doi", "lib_struct"):
        if k in doc and doc[k] is not None and not isinstance(doc[k], str):
            _issue(iss, "error", "C01", f"assay.{k}", f"'{k}' must be a string, got {type(doc[k]).__name__}")
    if not re.match(r"^\d+\.\d+\.\d+$", ver or ""):
        _issue(iss, "error", "C01", "assay.seqspec_version", f"seqspec_version {ver!r} is not X.Y.Z")
    if isinstance(doc.get("doi"), str) and doc["doi"] and not (doc["doi"].startswith("http") or doc["doi"].startswith("10.")):
        _issue(iss, "warning", "C01", "assay.doi", f"doi {doc['doi']!r} is not a DOI or URL")
    mods = modalities(doc)
    if not mods:
        _issue(iss, "error", "C01", "assay.modalities", "no modalities")
    for m in mods:
        if m not in MODALITIES:
            _issue(iss, "error", "C02", "assay.modalities", f"modality {m!r} not in {sorted(MODALITIES)}")
    if len(set(mods)) != len(mods):
        _issue(iss, "error", "C02", "assay.modalities", "duplicated modalities")
    lib = library_regions(doc)
    lib_ids = [str(r.get("region_id")) for r in lib]
    for m in mods:
        if m not in lib_ids:
            _issue(iss, "error", "C03", f"modality {m}", f"no {'assay_spec' if legacy else 'library_spec'} region with region_id '{m}'")
    rkeys = REGION_KEYS_V00 if legacy else REGION_KEYS_V02
    for mreg in lib:
        mid = str(mreg.get("region_id"))
        seen: Dict[str, int] = {}
        for r, depth, parent in walk(mreg):
            rid = str(r.get("region_id"))
            where = f"{mid}/{rid}"
            if r.get("__tag__") not in (None, "Region"):
                _issue(iss, "error", "C01", where, f"tagged !{r.get('__tag__')}, expected !Region")
            for k in rkeys:
                if k not in r:
                    _issue(iss, "error", "C01", where, f"region missing key '{k}'")
            seen[rid] = seen.get(rid, 0) + 1
            st, seq = r.get("sequence_type"), r.get("sequence")
            mn, mx = r.get("min_len"), r.get("max_len")
            for k, v in (("min_len", mn), ("max_len", mx)):
                if not isinstance(v, int) or isinstance(v, bool) or v < 0:
                    _issue(iss, "error", "C01", where, f"{k} must be a non-negative integer, got {v!r}")
            rt = r.get("region_type")
            if not legacy:
                if not rt:
                    _issue(iss, "error", "C02", where, "region_type is empty/null")
                elif rt not in REGION_TYPES and not (depth == 0 and rt in MODALITIES):
                    _issue(iss, "error", "C02", where, f"region_type {rt!r} not in the seqspec vocabulary")
            elif not rt:
                _issue(iss, "warning", "C02", where, "region_type is empty (seqspec 0.0.x allows it; kb/chromap typing falls back to region_id)")
            is_leaf = not children(r)
            if st not in SEQUENCE_TYPES:
                if is_leaf or not legacy and st is not None:
                    _issue(iss, "error" if is_leaf else "warning", "C02", where, f"sequence_type {st!r} not in {sorted(SEQUENCE_TYPES)}")
                elif not legacy:
                    _issue(iss, "warning", "C02", where, f"sequence_type {st!r} on a parent region (seqspec uses 'joined')")
            if isinstance(mn, int) and isinstance(mx, int) and mn > mx:
                _issue(iss, "error", "C12", where, f"min_len {mn} > max_len {mx}")
            if is_leaf and isinstance(seq, str) and isinstance(mn, int) and isinstance(mx, int) and st in ("fixed", "random", "onlist") and seq:
                if not (mn <= len(seq) <= mx):
                    _issue(iss, "error", "C10", where, f"len(sequence) = {len(seq)} outside [min_len {mn}, max_len {mx}]")
            if is_leaf and isinstance(seq, str) and seq:
                s = seq.upper()
                if st == "fixed" and set(s) - set("ACGT"):
                    _issue(iss, "error", "C11", where, f"fixed sequence contains non-ACGT symbols ({''.join(sorted(set(s) - set('ACGT')))}); use sequence_type random/onlist")
                elif st == "random" and not (set(s) <= {"N"} or set(s) <= {"X"}):
                    _issue(iss, "error", "C11", where, "random sequence must be all N or all X")
                elif st == "onlist" and set(s) - {"N"}:
                    _issue(iss, "error", "C11", where, "onlist sequence must be all N")
            elif is_leaf and st in ("fixed",) and not seq:
                _issue(iss, "error", "C11", where, "fixed region without a sequence")
            if st == "onlist" and not isinstance(r.get("onlist"), dict):
                _issue(iss, "error", "C08", where, "sequence_type onlist but no !Onlist")
            ol = r.get("onlist")
            if isinstance(ol, dict):
                fn = ol.get("filename")
                if not fn:
                    _issue(iss, "error", "C01", where, "onlist without filename")
                if not legacy:
                    if ol.get("location") not in ONLIST_LOCATIONS:
                        _issue(iss, "error", "C02", where, f"onlist location {ol.get('location')!r} not in {sorted(ONLIST_LOCATIONS)}")
                md5 = ol.get("md5")
                if md5 is not None and not re.match(r"^[0-9a-fA-F]{32}$", str(md5)):
                    _issue(iss, "error", "C09", where, f"onlist md5 {md5!r} is not a 32-character hex digest")
                elif md5 is None:
                    _issue(iss, "warning", "C09", where, "onlist md5 is null")
                if spec_dir is not None and fn:
                    p = _resolve_file(spec_dir, fn)
                    if p is None:
                        lvl = "warning" if ol.get("location") == "remote" else "error"
                        _issue(iss, lvl, "C08", where, f"onlist file {fn} not found in {spec_dir}")
                    else:
                        if p.name != fn:
                            _issue(iss, "warning", "C08", where, f"onlist {fn} found only as {p.name} (compression mismatch)")
                        if check_md5 and md5 and re.match(r"^[0-9a-fA-F]{32}$", str(md5)) and _md5(p) != str(md5).lower():
                            _issue(iss, "error", "C09", where, f"md5 of {p.name} does not match the spec")
            ch = children(r)
            if ch and isinstance(mn, int) and isinstance(mx, int):
                smin = sum(int(c.get("min_len") or 0) for c in ch)
                smax = sum(int(c.get("max_len") or 0) for c in ch)
                if (smin, smax) != (mn, mx) and not (legacy and depth <= 1):
                    _issue(iss, "error" if not legacy else "warning", "C13", where,
                           f"min/max_len {mn}/{mx} != sum of children {smin}/{smax}")
            if not legacy and depth > 0 and parent is not None and r.get("parent_id") != parent.get("region_id"):
                _issue(iss, "error", "C16", where, f"parent_id {r.get('parent_id')!r} != enclosing region {parent.get('region_id')!r}")
        for rid, n in seen.items():
            if n > 1:
                _issue(iss, "error", "C04", mid, f"region_id {rid!r} used {n} times")
    rds = reads(doc)
    # 0.0.x specs name the read regions after the FASTQ, so the same name may recur across modalities
    rkey = (lambda r: (str(r.get("modality")), str(r.get("read_id")))) if legacy else (lambda r: str(r.get("read_id")))
    rids = [rkey(r) for r in rds]
    for rid in sorted(set(rids), key=str):
        if rids.count(rid) > 1:
            _issue(iss, "error", "C05", "sequence_spec", f"read_id {rid!r} repeated")
    for rd in rds:
        where = f"read {rd.get('read_id')}"
        if not rd.get("_legacy"):
            if rd.get("__tag__") not in (None, "Read"):
                _issue(iss, "error", "C01", where, f"tagged !{rd.get('__tag__')}, expected !Read")
            for k, typ in READ_KEYS.items():
                if k not in rd:
                    _issue(iss, "error", "C01", where, f"read missing key '{k}'")
                elif not isinstance(rd[k], typ) or isinstance(rd[k], bool):
                    _issue(iss, "error", "C01", where, f"'{k}' must be {typ.__name__}, got {rd[k]!r}")
            if rd.get("strand") not in STRANDS:
                _issue(iss, "error", "C02", where, f"strand {rd.get('strand')!r} not in {sorted(STRANDS)}")
            if isinstance(rd.get("min_len"), int) and isinstance(rd.get("max_len"), int) and rd["min_len"] > rd["max_len"]:
                _issue(iss, "error", "C12", where, "min_len > max_len")
        mod = str(rd.get("modality"))
        if mod not in mods:
            _issue(iss, "error", "C07", where, f"modality {mod!r} is not declared in assay.modalities")
        mreg = modality_region(doc, mod)
        if not rd.get("_legacy") and mreg is not None:
            ids = {str(r.get("region_id")) for r, _, _ in walk(mreg)}
            if str(rd.get("primer_id")) not in ids:
                _issue(iss, "error", "C06", where, f"primer_id {rd.get('primer_id')!r} not found in library region '{mod}'")
            else:
                regs = project_read(doc, dict(rd, max_len=10 ** 9))
                avail = sum(r["stop"] - r["start"] for r in regs)
                if isinstance(rd.get("max_len"), int) and rd["max_len"] > avail:
                    _issue(iss, "warning", "C14", where, f"read max_len {rd['max_len']} exceeds the {avail} bp available after primer {rd.get('primer_id')}")
                cut = project_read(doc, rd)
                for r in cut:
                    if r["truncated"] and str(r["region_type"]).lower() in ("umi",) or (r["truncated"] and _is_barcode(str(r["region_type"]))):
                        _issue(iss, "warning", "C14", where, f"read ends inside {r['region_type']} region {r['region_id']} ({r['stop'] - r['start']} of {r['max_len']} bp read)")
        if spec_dir is not None:
            if _resolve_file(spec_dir, str(rd.get("read_id"))) is None:
                _issue(iss, "warning", "C15", where, f"FASTQ {rd.get('read_id')} not found in {spec_dir}")
    return iss


# ---------------------------------------------------------------------------
# Loading specs (path, embedded index, pinned commit)
# ---------------------------------------------------------------------------

def _compact_region(r: Dict[str, Any]) -> list:
    ol = r.get("onlist")
    olc = [ol.get("filename"), ol.get("location"), ol.get("md5")] if isinstance(ol, dict) else None
    return [r.get("region_id"), r.get("region_type"), r.get("sequence_type"), _rle(r.get("sequence") or ""),
            r.get("min_len"), r.get("max_len"), olc, [_compact_region(c) for c in children(r)]]


def _expand_region(c: list, parent: Optional[str] = None) -> Dict[str, Any]:
    rid, rt, st, seq, mn, mx, ol, ch = c
    r = {"__tag__": "Region", "parent_id": parent, "region_id": rid, "region_type": rt, "name": rid, "sequence_type": st,
         "sequence": _unrle(seq), "min_len": mn, "max_len": mx,
         "onlist": ({"__tag__": "Onlist", "filename": ol[0], "location": ol[1], "md5": ol[2]} if ol else None),
         "regions": [_expand_region(x, rid) for x in ch] or None}
    return r


def _compact(doc: Dict[str, Any], path: str, text: str) -> Dict[str, Any]:
    meta_keys = [k for k in doc if k not in ("sequence_spec", "library_spec", "assay_spec", "__tag__")]
    return {"path": path, "sha256": hashlib.sha256(text.encode()).hexdigest(), "legacy": is_legacy(doc),
            "meta": {k: doc[k] for k in meta_keys},
            "reads": [[r.get("read_id"), r.get("name"), r.get("modality"), r.get("primer_id"), r.get("min_len"), r.get("max_len"), r.get("strand")]
                      for r in (doc.get("sequence_spec") or [])],
            "library": [_compact_region(r) for r in library_regions(doc)]}


def _expand(entry: Dict[str, Any]) -> Dict[str, Any]:
    doc: Dict[str, Any] = {"__tag__": "Assay", **entry["meta"]}
    lib = [_expand_region(c) for c in entry["library"]]
    if entry["legacy"]:
        doc["assay_spec"] = lib
    else:
        doc["sequence_spec"] = [{"__tag__": "Read", "read_id": a, "name": b, "modality": c, "primer_id": d,
                                 "min_len": e, "max_len": f, "strand": g} for a, b, c, d, e, f, g in entry["reads"]]
        doc["library_spec"] = lib
    return doc


def embedded_index() -> List[Dict[str, Any]]:
    return json.loads(CATALOG_INDEX_JSON)


def fetch_file(rel: str, dest_root: Path = DATA_ROOT, force: bool = False) -> Path:
    dest = dest_root / rel
    if dest.exists() and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = RAW + rel
    with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 (pinned GitHub raw URL)
        data = resp.read()
    exp = REPO_FILES.get(rel)
    if exp and hashlib.md5(data).hexdigest() != exp[1]:
        raise RuntimeError(f"md5 mismatch for {rel} (expected {exp[1]})")
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)
    return dest


def resolve_spec(spec: str, allow_fetch: bool = True) -> Tuple[Dict[str, Any], Optional[Path], str, str]:
    """Return (doc, spec_dir, display name, source) for a path or a catalogue name/alias."""
    p = Path(spec).expanduser()
    if p.exists() and p.is_file():
        text = p.read_text()
        return parse_yaml(text), p.parent.resolve(), str(p), "file"
    key = SPEC_ALIASES.get(spec.lower(), spec)
    entries = {e["path"]: e for e in embedded_index()}
    if key not in entries:
        matches = [k for k in entries if spec.lower() in k.lower()]
        if len(matches) == 1:
            key = matches[0]
        else:
            raise SystemExit(f"spec {spec!r} is neither a file nor a catalogue entry ({', '.join(sorted(entries))}; aliases {', '.join(sorted(SPEC_ALIASES))})")
    local = DATA_ROOT / key
    if local.exists():
        return parse_yaml(local.read_text()), local.parent, key, "fetched copy"
    if allow_fetch:
        try:
            got = fetch_file(key)
            return parse_yaml(got.read_text()), got.parent, key, "pinned commit"
        except Exception as e:  # offline -> embedded index
            log.info("fetch of %s failed (%s); using the embedded index", key, e)
    return _expand(entries[key]), None, key, "embedded index"


# ---------------------------------------------------------------------------
# FASTQ checks
# ---------------------------------------------------------------------------

def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fastq(path: Path, n: int = 10000) -> List[str]:
    seqs = []
    with _open(path) as fh:
        while len(seqs) < n:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            seqs.append(s)
    return seqs


def load_onlist(path: Path) -> set:
    out = set()
    with _open(path) as fh:
        for line in fh:
            toks = line.strip().replace(",", "\t").split()
            for t in reversed(toks):
                if t and set(t.upper()) <= set("ACGTN"):
                    out.add(t.upper())
                    break
    return out


def _ham(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


def check_reads(doc: Dict[str, Any], spec_dir: Path, n: int = 10000, fastq_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    rows = []
    fdir = fastq_dir or spec_dir
    onl_cache: Dict[str, set] = {}
    for mod in modalities(doc):
        for rd, regs in index_modality(doc, mod):
            fq = _resolve_file(fdir, str(rd.get("read_id")))
            base = {"modality": mod, "read_id": rd.get("read_id")}
            if fq is None:
                rows.append({**base, "check": "fastq_exists", "region_id": "", "value": 0, "n": 0, "status": "FAIL",
                             "detail": f"not found in {fdir}"})
                continue
            seqs = read_fastq(fq, n)
            lens = [len(s) for s in seqs]
            mn, mx = int(rd.get("min_len") or 0), int(rd.get("max_len") or 0)
            inr = sum(1 for x in lens if mn <= x <= mx)
            frac = inr / max(len(lens), 1)
            rows.append({**base, "check": "read_length", "region_id": "", "value": round(frac, 4), "n": len(lens),
                         "status": "PASS" if frac >= 0.95 else "FAIL",
                         "detail": f"lengths {min(lens) if lens else 0}-{max(lens) if lens else 0} vs spec [{mn}, {mx}]"})
            for r in regs:
                sl = [s[r["start"]:r["stop"]] for s in seqs]
                if r["sequence_type"] == "fixed" and r["sequence"] and set(r["sequence"].upper()) <= set("ACGT"):
                    ref = r["sequence"].upper()[: r["stop"] - r["start"]]
                    hit = sum(1 for x in sl if _ham(x.upper(), ref) <= 1) / max(len(sl), 1)
                    rows.append({**base, "check": "fixed_match", "region_id": r["region_id"], "value": round(hit, 4), "n": len(sl),
                                 "status": "PASS" if hit >= 0.5 else "FAIL", "detail": f"{ref} at {r['start']}-{r['stop']} (Hamming <= 1)"})
                ol = r.get("onlist")
                if isinstance(ol, dict) and ol.get("filename"):
                    op = _resolve_file(spec_dir, ol["filename"])
                    if op is None:
                        rows.append({**base, "check": "onlist_match", "region_id": r["region_id"], "value": "", "n": 0,
                                     "status": "SKIP", "detail": f"onlist {ol['filename']} not available"})
                        continue
                    if str(op) not in onl_cache:
                        onl_cache[str(op)] = load_onlist(op)
                    onl = onl_cache[str(op)]
                    L = {len(x) for x in onl}
                    w = r["stop"] - r["start"]
                    k = min(L | {w}) if L else w
                    trimmed = {x[:k] for x in onl} if L and L != {w} else onl
                    hit = sum(1 for x in sl if x.upper()[:k] in trimmed) / max(len(sl), 1)
                    rows.append({**base, "check": "onlist_match", "region_id": r["region_id"], "value": round(hit, 4), "n": len(sl),
                                 "status": "PASS" if hit >= 0.5 else "WARN",
                                 "detail": f"{op.name}: {len(onl):,} entries (length {sorted(L)[:3]}), slice {r['start']}-{r['stop']}" + (f", compared on {k}-bp prefixes" if L != {w} else "")})
    return rows


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

def spec_summary(doc: Dict[str, Any], path: str, spec_dir: Optional[Path]) -> Dict[str, Any]:
    iss = check_spec(doc, spec_dir, check_md5=spec_dir is not None)
    mods = modalities(doc)
    kb = {}
    for m in mods:
        try:
            kb[m] = format_index(doc, m, "kb")["string"]
        except Exception as e:
            kb[m] = f"error: {e}"
    return {
        "path": path, "assay_id": doc.get("assay_id") or doc.get("assay"), "name": doc.get("name"),
        "seqspec_version": spec_version(doc), "schema": "0.0.x (assay_spec)" if is_legacy(doc) else "0.2.x (library_spec)",
        "doi": doc.get("doi"), "date": doc.get("date") or doc.get("publication_date"), "description": doc.get("description"),
        "modalities": mods, "n_reads": len(reads(doc)),
        "n_regions": sum(1 for m in library_regions(doc) for _ in walk(m)),
        "onlists": sorted({o["filename"] for o in onlists(doc) if o["filename"]}),
        "kb": kb, "n_errors": sum(1 for i in iss if i["level"] == "error"),
        "n_warnings": sum(1 for i in iss if i["level"] == "warning"), "issues": iss,
    }


def build_catalog(source: str = "embedded", from_dir: Optional[Path] = None) -> List[Tuple[Dict[str, Any], str, Optional[Path], str]]:
    out = []
    if from_dir:
        for p in sorted(Path(from_dir).rglob("*")):
            if p.suffix in (".yml", ".yaml") and ".git" not in p.parts:
                doc = parse_yaml(p.read_text())
                if isinstance(doc, dict) and ("library_spec" in doc or "assay_spec" in doc):
                    out.append((doc, str(p.relative_to(from_dir)), p.parent, "local dir"))
        return out
    for e in embedded_index():
        if source == "live":
            text = urllib.request.urlopen(RAW + e["path"], timeout=60).read().decode()  # noqa: S310
            if hashlib.sha256(text.encode()).hexdigest() != e["sha256"]:
                log.warning("%s at the pinned commit differs from the embedded index", e["path"])
            out.append((parse_yaml(text), e["path"], (DATA_ROOT / e["path"]).parent if (DATA_ROOT / e["path"]).exists() else None, "pinned commit"))
        else:
            local = DATA_ROOT / e["path"]
            if local.exists():
                out.append((parse_yaml(local.read_text()), e["path"], local.parent, "fetched copy"))
            else:
                out.append((_expand(e), e["path"], None, "embedded index"))
    return out


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_catalog(args: argparse.Namespace) -> int:
    t0 = time.time()
    specs = build_catalog("live" if args.live else "embedded", Path(args.from_dir) if args.from_dir else None)
    if not specs:
        print("no seqspec files found")
        return 1
    out = run_dir(args.label)
    cat, rds, regs, ols, full = [], [], [], [], []
    for doc, path, sdir, src in specs:
        s = spec_summary(doc, path, sdir)
        s["source"] = src
        full.append(s)
        for m in s["modalities"]:
            cat.append({"spec": path, "assay_id": s["assay_id"], "name": s["name"], "seqspec_version": s["seqspec_version"],
                        "schema": s["schema"], "modality": m, "doi": s["doi"], "date": s["date"],
                        "n_reads": len(reads(doc, m)), "kb_string": s["kb"].get(m), "onlists": ";".join(s["onlists"]),
                        "check_errors": s["n_errors"], "check_warnings": s["n_warnings"], "source": src})
        for rd in reads(doc):
            rds.append({"spec": path, "read_id": rd.get("read_id"), "modality": rd.get("modality"), "primer_id": rd.get("primer_id"),
                        "strand": rd.get("strand"), "min_len": rd.get("min_len"), "max_len": rd.get("max_len"),
                        "regions_read": ",".join(f"{r['region_id']}[{r['start']}:{r['stop']}]" for r in project_read(doc, rd))})
        for mreg in library_regions(doc):
            for r, depth, parent in walk(mreg):
                ol = r.get("onlist") if isinstance(r.get("onlist"), dict) else {}
                regs.append({"spec": path, "modality": mreg.get("region_id"), "region_id": r.get("region_id"),
                             "parent": parent.get("region_id") if parent else "", "depth": depth, "region_type": r.get("region_type"),
                             "sequence_type": r.get("sequence_type"), "min_len": r.get("min_len"), "max_len": r.get("max_len"),
                             "sequence": (r.get("sequence") or "")[:60], "onlist": ol.get("filename")})
        for o in onlists(doc, sdir):
            ols.append({"spec": path, **o})
    write_tsv(out / "catalog.tsv", cat)
    write_tsv(out / "reads.tsv", rds)
    write_tsv(out / "regions.tsv", regs)
    write_tsv(out / "onlists.tsv", ols, ["spec", "modality", "region_id", "region_type", "filename", "location", "md5", "local_path"])
    write_json(out / "catalog.json", full)
    summ = {"upstream": UPSTREAM_REPO, "commit": UPSTREAM_COMMIT, "n_specs": len(full), "n_modality_rows": len(cat),
            "specs": [{k: s[k] for k in ("path", "assay_id", "seqspec_version", "modalities", "kb", "onlists", "n_errors", "n_warnings", "source")} for s in full]}
    write_json(out / "summary.json", summ)
    lines = [f"# CRISPR-SeqSpec catalogue ({UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]})", "",
             f"{len(full)} seqspec file(s), {len(cat)} modality row(s).  kb strings are `seqspec index -t kb` equivalents; "
             "check counts are `seqspec check` equivalents (file checks only where the files are present locally).", "",
             "| spec | assay | version | modality | kb read format | onlists | errors / warnings |", "|---|---|---|---|---|---|---|"]
    for c in cat:
        lines.append(f"| {c['spec']} | {c['assay_id']} | {c['seqspec_version']} | {c['modality']} | `{c['kb_string']}` | {c['onlists']} | {c['check_errors']} / {c['check_warnings']} |")
    lines += ["", "## Upstream files at the pinned commit", "", "| file | bytes |", "|---|---|"]
    lines += [f"| {p} | {s:,} |" for p, (s, _) in sorted(REPO_FILES.items())]
    write_report(out / "report.md", lines)
    print(f"{len(full)} specs: " + "; ".join(f"{s['path']} [{','.join(s['modalities'])}] errors={s['n_errors']}" for s in full))
    print(f"Done in {time.time() - t0:.1f}s -> {out}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    doc, sdir, name, src = resolve_spec(args.spec, allow_fetch=not args.offline)
    if args.spec_dir:
        sdir = Path(args.spec_dir)
    iss = check_spec(doc, sdir, check_md5=not args.no_md5)
    out = run_dir(args.label)
    write_tsv(out / "issues.tsv", iss, ["level", "code", "where", "message"])
    ne = sum(1 for i in iss if i["level"] == "error")
    nw = len(iss) - ne
    for i in iss:
        print(f"  [{i['level']}] {i['code']} {i['where']}: {i['message']}")
    ok = ne == 0
    write_json(out / "summary.json", {"spec": name, "source": src, "spec_dir": str(sdir) if sdir else None,
                                      "seqspec_version": spec_version(doc), "valid": ok, "n_errors": ne, "n_warnings": nw})
    lines = [f"# seqspec check: {name}", "", f"Source: {src}; seqspec {spec_version(doc)}; file checks in: {sdir or 'not run (no local copy)'}", "",
             f"**{'VALID' if ok else 'INVALID'}** -- {ne} error(s), {nw} warning(s).", "",
             "| level | code | where | message |", "|---|---|---|---|"]
    lines += [f"| {i['level']} | {i['code']} | {i['where']} | {i['message']} |" for i in iss]
    write_report(out / "report.md", lines)
    print(f"{name}: {'valid' if ok else 'INVALID'} ({ne} errors, {nw} warnings)")
    return 0 if (ok or not args.strict) else 1


def cmd_index(args: argparse.Namespace) -> int:
    doc, sdir, name, src = resolve_spec(args.spec, allow_fetch=not args.offline)
    mods = [args.modality] if args.modality else modalities(doc)
    out = run_dir(args.label)
    results = []
    for m in mods:
        if m not in modalities(doc) and modality_region(doc, m) is None:
            print(f"modality {m!r} not in {modalities(doc)}")
            return 1
        try:
            res = format_index(doc, m, args.tool, dedupe=args.dedupe, feature_types=args.feature_types)
        except ValueError as e:
            print(f"{name} [{m}] {args.tool}: cannot format ({e}){'' if args.dedupe else '; try --dedupe'}")
            return 1
        results.append(res)
        print(f"{name} [{m}] {args.tool}: {res['string']}")
        if res.get("command_args"):
            print(f"  args: {res['command_args']}")
        for n in res["notes"]:
            print(f"  note: {n}")
    rows = [dict(modality=r["modality"], **x) for r in results for x in r["regions"]]
    write_tsv(out / "index_regions.tsv", rows, ["modality", "fastq", "fastq_index", "region_id", "region_type", "start", "stop", "truncated"])
    write_json(out / "summary.json", {"spec": name, "source": src, "tool": args.tool, "dedupe": args.dedupe,
                                      "results": [{k: r[k] for k in r if k != "regions"} for r in results]})
    lines = [f"# seqspec index -t {args.tool}: {name}", "", f"Source: {src}", ""]
    for r in results:
        lines += [f"## {r['modality']}", "", f"FASTQs (index order): {', '.join(r['fastqs'])}", "", "```", r["string"], "```", ""]
        if r.get("command_args"):
            lines += [f"Arguments: `{r['command_args']}`", ""]
        lines += [f"- note: {n}" for n in r["notes"]]
    write_report(out / "report.md", lines)
    return 0


def cmd_onlist(args: argparse.Namespace) -> int:
    doc, sdir, name, src = resolve_spec(args.spec, allow_fetch=not args.offline)
    rows = onlists(doc, sdir)
    if args.modality:
        rows = [r for r in rows if r["modality"] == args.modality]
    out = run_dir(args.label)
    write_tsv(out / "onlists.tsv", rows, ["modality", "region_id", "region_type", "filename", "location", "md5", "local_path"])
    for r in rows:
        print(f"{r['modality']}\t{r['region_id']}\t{r['filename']}\t{r['local_path'] or '(not local)'}")
    write_json(out / "summary.json", {"spec": name, "source": src, "onlists": rows,
                                      "filenames": sorted({r['filename'] for r in rows if r['filename']})})
    write_report(out / "report.md", [f"# seqspec onlist: {name}", "", "| modality | region | file | location | md5 | local |", "|---|---|---|---|---|---|"] +
                 [f"| {r['modality']} | {r['region_id']} | {r['filename']} | {r['location']} | {r['md5']} | {r['local_path'] or ''} |" for r in rows])
    return 0


def render_info(doc: Dict[str, Any]) -> List[str]:
    lines = [f"{doc.get('name')} (assay {doc.get('assay_id') or doc.get('assay')}, seqspec {spec_version(doc)})",
             f"  modalities: {', '.join(modalities(doc))}"]
    for mreg in library_regions(doc):
        lines.append(f"  library [{mreg.get('region_id')}]")
        for r, depth, _ in walk(mreg):
            if depth == 0:
                continue
            seq = r.get("sequence") or ""
            lines.append(f"  {'  ' * depth}- {r.get('region_id')} ({rtype(r)}, {r.get('sequence_type')}, {r.get('min_len')}-{r.get('max_len')} bp)"
                         + (f" {seq[:30]}{'...' if len(seq) > 30 else ''}" if seq else "")
                         + (f" onlist={r['onlist'].get('filename')}" if isinstance(r.get("onlist"), dict) else ""))
        for rd in reads(doc, str(mreg.get("region_id"))):
            regs = project_read(doc, rd)
            bar = "|".join(f"{x['region_id']}:{x['start']}-{x['stop']}" for x in regs)
            lines.append(f"    read {rd.get('read_id')} ({rd.get('strand')}, {rd.get('min_len')}-{rd.get('max_len')} bp, primer {rd.get('primer_id')}): {bar}")
    return lines


def cmd_info(args: argparse.Namespace) -> int:
    doc, sdir, name, src = resolve_spec(args.spec, allow_fetch=not args.offline)
    print(f"# {name} ({src})")
    for line in render_info(doc):
        print(line)
    return 0


def cmd_check_reads(args: argparse.Namespace) -> int:
    doc, sdir, name, src = resolve_spec(args.spec, allow_fetch=not args.offline)
    if args.spec_dir:
        sdir = Path(args.spec_dir)
    if sdir is None:
        print("no local copy of the spec's files: run `igvfagent crispr-seqspec fetch` or pass a spec path / --spec-dir")
        return 1
    rows = check_reads(doc, sdir, n=args.n_reads, fastq_dir=Path(args.fastq_dir) if args.fastq_dir else None)
    out = run_dir(args.label)
    write_tsv(out / "read_checks.tsv", rows, ["modality", "read_id", "check", "region_id", "value", "n", "status", "detail"])
    nfail = sum(1 for r in rows if r["status"] == "FAIL")
    for r in rows:
        print(f"  {r['status']:4s} {r['modality']}/{r['read_id']} {r['check']} {r['region_id']} {r['value']} ({r['detail']})")
    write_json(out / "summary.json", {"spec": name, "spec_dir": str(sdir), "n_checks": len(rows), "n_fail": nfail,
                                      "n_warn": sum(1 for r in rows if r["status"] == "WARN")})
    write_report(out / "report.md", [f"# FASTQ vs seqspec: {name}", "", f"{len(rows)} checks, {nfail} failed (first {args.n_reads:,} reads per FASTQ).", "",
                                     "| status | modality | read | check | region | value | detail |", "|---|---|---|---|---|---|---|"] +
                 [f"| {r['status']} | {r['modality']} | {r['read_id']} | {r['check']} | {r['region_id']} | {r['value']} | {r['detail']} |" for r in rows])
    return 0 if nfail == 0 or not args.strict else 1


def cmd_fetch(args: argparse.Namespace) -> int:
    dest = Path(args.dest) if args.dest else DATA_ROOT
    todo = [p for p in REPO_FILES if (not args.specs_only or p.endswith((".yml", ".yaml", ".md")))]
    total = sum(REPO_FILES[p][0] for p in todo)
    print(f"Fetching {len(todo)} file(s), {total / 1e6:.1f} MB, from {UPSTREAM_REPO}@{UPSTREAM_COMMIT[:7]} into {dest}")
    bad = 0
    for p in todo:
        try:
            got = fetch_file(p, dest, force=args.force)
            print(f"Wrote: {got}")
        except Exception as e:
            bad += 1
            print(f"  failed {p}: {e}")
    return 0 if bad == 0 else 1


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

SYN_V02 = """!Assay
seqspec_version: 0.2.0
assay_id: SYN_CROPSEQ
name: SYN_CROPSEQ
doi: https://doi.org/10.0000/synthetic
date: 01 January 2024
description: synthetic CROP-seq guide library  # comment
modalities:
- crispr
lib_struct: ''
sequence_protocol: Not-specified
sequence_kit: Not-specified
library_protocol: synthetic
library_kit: Not-specified
sequence_spec:
- !Read
  read_id: syn_R1.fq.gz
  name: Read 1
  modality: crispr
  primer_id: truseq_read1
  min_len: 28
  max_len: 28
  strand: pos
- !Read
  read_id: syn_R2.fq.gz
  name: Read 2
  modality: crispr
  primer_id: truseq_read2
  min_len: 43
  max_len: 43
  strand: neg
library_spec:
- !Region
  parent_id: null
  region_id: crispr
  region_type: crispr
  name: crispr
  sequence_type: joined
  sequence: ''
  min_len: 138
  max_len: 138
  onlist: null
  regions:
  - !Region
    parent_id: crispr
    region_id: truseq_read1
    region_type: truseq_read1
    name: truseq_read1
    sequence_type: fixed
    sequence: ACACTCTTTCCCTACACGACGCTCTTCCGATCT
    min_len: 33
    max_len: 33
    onlist: null
    regions: null
  - !Region
    parent_id: crispr
    region_id: barcode
    region_type: barcode
    name: barcode
    sequence_type: onlist
    sequence: NNNNNNNNNNNNNNNN
    min_len: 16
    max_len: 16
    onlist: !Onlist
      location: local
      filename: syn_onlist.txt
      md5: __MD5__
    regions: null
  - !Region
    parent_id: crispr
    region_id: umi
    region_type: umi
    name: umi
    sequence_type: random
    sequence: NNNNNNNNNNNN
    min_len: 12
    max_len: 12
    onlist: null
    regions: null
  - !Region
    parent_id: crispr
    region_id: sgrna_target
    region_type: sgrna_target
    name: sgrna_target
    sequence_type: onlist
    sequence: NNNNNNNNNNNNNNNNNNNN
    min_len: 20
    max_len: 20
    onlist: !Onlist
      location: local
      filename: syn_guides.txt
      md5: null
    regions: null
  - !Region
    parent_id: crispr
    region_id: scaffold
    region_type: linker
    name: scaffold
    sequence_type: fixed
    sequence: CTTGTGGAAAGGACGAAACACCG
    min_len: 23
    max_len: 23
    onlist: null
    regions: null
  - !Region
    parent_id: crispr
    region_id: truseq_read2
    region_type: truseq_read2
    name: truseq_read2
    sequence_type: fixed
    sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    min_len: 34
    max_len: 34
    onlist: null
    regions: null
"""

SYN_V00 = """!Assay
assay: SYN-TAP
sequencer: Illumina
seqspec_version: 0.0.3
name: SYN_TAP
doi: https://doi.org/10.0000/tap
publication_date: 01 June 2020
description: synthetic legacy spec
modalities:
- rna
- crispr
lib_struct: https://example.org
assay_spec:
- !Region
  region_id: rna
  region_type: ''
  name: rna
  sequence_type: ''
  sequence: ''
  min_len: 0
  max_len: 1024
  onlist: null
  regions:
  - !Region
    region_id: R1.fastq.gz
    region_type: ''
    name: R1.fastq.gz
    sequence_type: ''
    sequence: ''
    min_len: 0
    max_len: 1024
    onlist: null
    regions:
    - !Region
      region_id: barcode
      region_type: ''
      name: barcode
      sequence_type: onlist
      sequence: NNNNNNNNNNNNNNNN
      min_len: 16
      max_len: 16
      onlist: !Onlist
        filename: wl.txt.gz
        md5: null
      regions: null
    - !Region
      region_id: umi
      region_type: ''
      name: umi
      sequence_type: random
      sequence: NNNNNNNNNN
      min_len: 10
      max_len: 10
      onlist: null
      regions: null
  - !Region
    region_id: R2.fastq.gz
    region_type: ''
    name: R2.fastq.gz
    sequence_type: ''
    sequence: ''
    min_len: 0
    max_len: 1024
    onlist: null
    regions:
    - !Region
      region_id: cdna
      region_type: ''
      name: cdna
      sequence_type: random
      sequence: ''
      min_len: 58
      max_len: 58
      onlist: null
      regions: null
- !Region
  region_id: crispr
  region_type: ''
  name: crispr
  sequence_type: ''
  sequence: ''
  min_len: 0
  max_len: 1024
  onlist: null
  regions:
  - !Region
    region_id: R1.fastq.gz
    region_type: ''
    name: R1.fastq.gz
    sequence_type: ''
    sequence: ''
    min_len: 0
    max_len: 1024
    onlist: null
    regions:
    - !Region
      region_id: barcode
      region_type: ''
      name: barcode
      sequence_type: onlist
      sequence: NNNNNNNNNNNNNNNN
      min_len: 16
      max_len: 16
      onlist: !Onlist
        filename: wl.txt.gz
        md5: null
      regions: null
    - !Region
      region_id: umi
      region_type: ''
      name: umi
      sequence_type: random
      sequence: NNNNNNNNNN
      min_len: 10
      max_len: 10
      onlist: null
      regions: null
  - !Region
    region_id: R2.fastq.gz
    region_type: ''
    name: R2.fastq.gz
    sequence_type: ''
    sequence: ''
    min_len: 0
    max_len: 1024
    onlist: null
    regions:
    - !Region
      region_id: cdna
      region_type: ''
      name: cdna
      sequence_type: random
      sequence: ''
      min_len: 18
      max_len: 18
      onlist: null
      regions: null
    - !Region
      region_id: sgrna
      region_type: ''
      name: sgrna
      sequence_type: onlist
      sequence: NNNNNNNNNNNNNNNNNNNN
      min_len: 20
      max_len: 20
      onlist: !Onlist
        filename: guides.txt
        md5: null
      regions: null
"""


def _rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def cmd_selftest(args: argparse.Namespace) -> int:
    import random
    import shutil
    import tempfile
    rng = random.Random(7)
    checks: List[Tuple[bool, str]] = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made: List[Path] = []

    def newest(label):
        runs = sorted(OUT_ROOT.glob(f"*_{label}"))
        if runs:
            made.append(runs[-1])
        return runs[-1] if runs else None

    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        print("\nYAML reader")
        y = parse_yaml("!Assay\na: 1\nb: 'x y'\nc:\n- p\n- q\nd: !Onlist\n  filename: f.txt\n  md5: null\ne: [1, 2]\n")
        check(y == {"__tag__": "Assay", "a": 1, "b": "x y", "c": ["p", "q"], "d": {"filename": "f.txt", "md5": None, "__tag__": "Onlist"}, "e": [1, 2]},
              "tagged block YAML: scalars, quoted strings, same-indent sequences, nested tagged maps, flow lists")
        emb = embedded_index()
        check(len(emb) == 4 and {e["path"] for e in emb} == set(SPEC_FILES), "embedded index covers the 4 upstream seqspec files")

        print("\nupstream specs (embedded index)")
        g = _expand(next(e for e in emb if e["path"].endswith("guide.yml")))
        check(format_index(g, "guide", "kb")["string"] == "0,0,16:0,16,26:1,23,43",
              "guide.yml kb: barcode 0-16, UMI cut to 16-26 by the 26-bp R1, protospacer at R2 23-43 after the 23-bp scaffold")
        gi = check_spec(g)
        codes = {i["code"] for i in gi if i["level"] == "error"}
        check({"C02", "C11", "C09"} <= codes, "guide.yml check flags modality 'guide' / region types (C02), fixed UMI of N (C11), md5 'updateit' (C09)")
        check(any(i["code"] == "C14" and "umi" in i["message"] for i in gi), "guide.yml check warns that R1 (26 bp) ends inside the 12-bp UMI")
        r = _expand(next(e for e in emb if e["path"].endswith("rna.yml")))
        check(format_index(r, "rna", "kb")["string"] == "0,0,16:0,16,28:1,0,150", "rna.yml kb = 10x v3 layout 0,0,16:0,16,28:1,0,150")
        m = _expand(next(e for e in emb if e["path"].endswith("multiseq.yml")))
        mk = format_index(m, "multiseq", "kb")
        mkd = format_index(m, "multiseq", "kb", dedupe=True)
        check(mk["notes"] and mkd["string"] == "0,0,16:0,16,26:1,0,8",
              "multiseq.yml: 150-bp R2 re-reads barcode/UMI (noted); --dedupe gives 0,0,16:0,16,26:1,0,8")
        t = _expand(next(e for e in emb if e["path"].endswith("spec.yaml")))
        check(is_legacy(t) and format_index(t, "crispr", "kb")["string"] == "0,0,16:0,16,26:1,18,38"
              and format_index(t, "rna", "kb")["string"] == "0,0,16:0,16,26:1,0,58",
              "TAP-seq 0.0.3 spec: rna 0,0,16:0,16,26:1,0,58; crispr feature = sgrna 1,18,38")
        check(format_index(t, "rna", "chromap")["string"] == "bc:0:15,r1:0:57", "chromap read-format bc:0:15,r1:0:57 (inclusive ends)")
        check(sorted({o["filename"] for o in onlists(t)}) == ["10x_bc_whitelist_737k_201608.txt.gz", "guides_schraivogel2020.txt.gz"],
              "onlist filenames of the TAP-seq spec")

        print("\nsynthetic 0.2.0 spec + FASTQs")
        wl = [rnd(16) for _ in range(200)]
        guides = [rnd(20) for _ in range(30)]
        (d / "syn_onlist.txt").write_text("\n".join(wl) + "\n")
        (d / "syn_guides.txt").write_text("\n".join(f"g{i}\t{s}" for i, s in enumerate(guides)) + "\n")
        md5 = _md5(d / "syn_onlist.txt")
        (d / "spec.yml").write_text(SYN_V02.replace("__MD5__", md5))
        scaffold = "CTTGTGGAAAGGACGAAACACCG"
        r1, r2 = [], []
        for i in range(400):
            bc = wl[i % len(wl)] if i < 360 else rnd(16)
            r1.append(bc + rnd(12))
            r2.append(scaffold + guides[i % len(guides)])
        for name, seqs in (("syn_R1.fq.gz", r1), ("syn_R2.fq.gz", r2)):
            with gzip.open(d / name, "wt") as fh:
                for i, s in enumerate(seqs):
                    fh.write(f"@r{i}\n{s}\n+\n{'F' * len(s)}\n")
        doc, sdir, _, _ = resolve_spec(str(d / "spec.yml"))
        iss = check_spec(doc, sdir)
        check(not [i for i in iss if i["level"] == "error"], f"valid synthetic spec passes check ({len(iss)} warnings)")
        k = format_index(doc, "crispr", "kb")
        check(k["string"] == "0,0,16:0,16,28:1,23,43", f"auto feature = sgrna_target for crispr modality: {k['string']}")
        k2 = format_index(doc, "crispr", "kb", feature_types=["linker"])
        check(k2["string"].endswith(":1,0,23"), "--feature-types overrides the feature region")
        s = format_index(doc, "crispr", "starsolo")["string"]
        check("--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12" in s, "STARsolo CB/UMI arguments")
        rows = check_reads(doc, sdir)
        get = lambda c, rid="": next((x for x in rows if x["check"] == c and x["region_id"] == rid), {})
        check(abs(float(get("onlist_match", "barcode").get("value", 0)) - 0.9) < 1e-9, "check-reads: 90% of planted barcodes found in the onlist")
        check(float(get("fixed_match", "scaffold").get("value", 0)) == 1.0 and float(get("onlist_match", "sgrna_target").get("value", 0)) == 1.0,
              "check-reads: scaffold fixed match 100%, guides 100% in the guide onlist")
        check(all(x["status"] == "PASS" for x in rows if x["check"] == "read_length"), "check-reads: read lengths inside [min_len, max_len]")

        print("\nplanted defects")
        bad = (SYN_V02.replace("__MD5__", "a" * 32).replace("modality: crispr\n  primer_id: truseq_read2", "modality: crispr\n  primer_id: no_primer")
               .replace("sequence: NNNNNNNNNNNN\n    min_len: 12", "sequence: NNNNNNNNNNNNN\n    min_len: 12")
               .replace("region_id: scaffold\n    region_type: linker", "region_id: umi\n    region_type: linker")
               .replace("min_len: 138", "min_len: 137").replace("strand: neg", "strand: minus"))
        (d / "bad.yml").write_text(bad)
        bdoc, bdir, _, _ = resolve_spec(str(d / "bad.yml"))
        bc = {i["code"] for i in check_spec(bdoc, bdir) if i["level"] == "error"}
        for code, what in (("C06", "primer_id missing"), ("C09", "md5 mismatch"), ("C10", "sequence longer than max_len"),
                           ("C04", "duplicate region_id"), ("C13", "parent length != children"), ("C02", "strand enum")):
            check(code in bc, f"check catches {what} ({code})")
        (d / "syn_onlist.txt").unlink()
        c8 = [i for i in check_spec(doc, sdir) if i["code"] == "C08"]
        check(c8 and c8[0]["level"] == "error", "check catches a missing local onlist (C08)")
        (d / "legacy").mkdir()
        (d / "legacy" / "spec.yaml").write_text(SYN_V00)
        ldoc, ldir, _, _ = resolve_spec(str(d / "legacy" / "spec.yaml"))
        check(format_index(ldoc, "crispr", "kb")["string"] == "0,0,16:0,16,26:1,18,38", "legacy 0.0.3 layout parsed from YAML text")

        print("\nsubcommands")
        rc = cmd_check(argparse.Namespace(spec=str(d / "bad.yml"), spec_dir=None, no_md5=False, offline=True, strict=True, label="st_seqspec_check"))
        run = newest("st_seqspec_check")
        check(rc == 1 and run and (run / "issues.tsv").exists() and (run / "report.md").exists(), "check --strict exits 1 on an invalid spec and writes issues.tsv/report.md")
        rc = cmd_index(argparse.Namespace(spec=str(d / "legacy" / "spec.yaml"), modality="crispr", tool="kb", dedupe=False, feature_types=None, offline=True, label="st_seqspec_index"))
        run = newest("st_seqspec_index")
        sj = json.loads((run / "summary.json").read_text()) if run else {}
        check(rc == 0 and sj.get("results", [{}])[0].get("string") == "0,0,16:0,16,26:1,18,38", "index subcommand writes the kb string to summary.json")
        rc = cmd_catalog(argparse.Namespace(live=False, from_dir=None, label="st_seqspec_catalog"))
        run = newest("st_seqspec_catalog")
        ok = False
        if run:
            import pandas as pd
            cat = pd.read_csv(run / "catalog.tsv", sep="\t")
            ok = rc == 0 and len(cat) == 5 and set(cat["modality"]) == {"rna", "guide", "multiseq", "crispr"}
        check(ok, "catalog: 4 specs x modalities = 5 rows (rna, guide, multiseq, TAP-seq rna + crispr)")
        rc = cmd_catalog(argparse.Namespace(live=False, from_dir=str(d / "legacy"), label="st_seqspec_catalog_dir"))
        newest("st_seqspec_catalog_dir")
        check(rc == 0, "catalog --from-dir scans a local directory")
        rc = cmd_onlist(argparse.Namespace(spec=str(d / "spec.yml"), modality=None, offline=True, label="st_seqspec_onlist"))
        run = newest("st_seqspec_onlist")
        check(rc == 0 and run and "syn_guides.txt" in (run / "onlists.tsv").read_text(), "onlist subcommand lists both onlists")
        rc = cmd_check_reads(argparse.Namespace(spec=str(d / "spec.yml"), spec_dir=None, fastq_dir=None, n_reads=1000, offline=True, strict=False, label="st_seqspec_reads"))
        newest("st_seqspec_reads")
        check(rc == 0, "check-reads subcommand runs with the onlist missing (SKIP rows)")
        check(any("read syn_R1.fq.gz" in line for line in render_info(doc)), "info renders library + read layout")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent crispr-seqspec", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def spec_arg(p):
        p.add_argument("--spec", required=True, help="seqspec YAML path, or a catalogue entry / alias (rna, guide, multiseq, tapseq)")
        p.add_argument("--offline", action="store_true", help="never fetch from GitHub; use the embedded index")

    p = sub.add_parser("catalog", help="catalogue every seqspec in the upstream repo")
    p.add_argument("--live", action="store_true", help="re-read the specs from the pinned commit")
    p.add_argument("--from-dir", help="scan a local directory (e.g. a clone) instead")
    p.add_argument("--label", default="catalog")
    p.set_defaults(func=cmd_catalog)

    p = sub.add_parser("check", help="validate a seqspec (seqspec check)")
    spec_arg(p)
    p.add_argument("--spec-dir", help="directory holding the onlists / FASTQs (default: the spec's directory)")
    p.add_argument("--no-md5", action="store_true", help="skip md5 verification of local onlists")
    p.add_argument("--strict", action="store_true", help="exit 1 when there are errors")
    p.add_argument("--label", default="check")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("index", help="read-format string for a tool (seqspec index -t)")
    spec_arg(p)
    p.add_argument("--modality", help="modality (default: all)")
    p.add_argument("--tool", default="kb", choices=["kb", "chromap", "starsolo", "tab"])
    p.add_argument("--dedupe", action="store_true", help="keep each region only in the first read that covers it")
    p.add_argument("--feature-types", nargs="+", help="region types to report as the feature (default: cdna/gdna, or sgRNA regions for crispr)")
    p.add_argument("--label", default="index")
    p.set_defaults(func=cmd_index)

    p = sub.add_parser("onlist", help="onlist filenames of a spec")
    spec_arg(p)
    p.add_argument("--modality")
    p.add_argument("--label", default="onlist")
    p.set_defaults(func=cmd_onlist)

    p = sub.add_parser("info", help="print the library structure")
    spec_arg(p)
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("check-reads", help="validate FASTQs against a spec")
    spec_arg(p)
    p.add_argument("--spec-dir", help="directory with the onlists (default: the spec's directory)")
    p.add_argument("--fastq-dir", help="directory with the FASTQs (default: --spec-dir)")
    p.add_argument("--n-reads", type=int, default=10000)
    p.add_argument("--strict", action="store_true")
    p.add_argument("--label", default="check_reads")
    p.set_defaults(func=cmd_check_reads)

    p = sub.add_parser("fetch", help="download the upstream files at the pinned commit")
    p.add_argument("--dest", help=f"destination (default {DATA_ROOT})")
    p.add_argument("--specs-only", action="store_true")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("selftest", help="synthetic specs + FASTQs, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true", help="accepted for interface parity (no figures are drawn)")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
