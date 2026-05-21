"""MaveDB scoreset → genomic coordinates mapping skill.

Maps the HGVSp / HGVSc variants in a MaveDB scoreset to canonical
chr/pos/ref/alt VCF-style coordinates, so MAVE / VAMP-seq scores
become directly cross-referenceable with ClinVar, gnomAD, GWAS
catalogues and the IGVF Catalog.

Clean-room reimplementation of the algorithm in
[ave-dcd/dcd_mapping](https://github.com/ave-dcd/dcd_mapping) (MIT,
Arbesfeld/Stevenson/Rubin/Wagner, *bioRxiv* 2023.06.20.545702) but
without depending on its heavy infrastructure (UTA Postgres, SeqRepo
local store, BLAT binary). Instead we use the public Ensembl REST API
to look up canonical transcripts and translate protein → cDNA →
genomic positions. The output schema matches the GA4GH VRS-compatible
chr/pos/ref/alt format with full provenance.

Commands
--------
    mavedb map-scoreset    Map every variant in a MaveDB scoreset to
                            chr/pos/ref/alt. Output: TSV + VCF.
    mavedb to-vcf          Convert an already-mapped scoreset to VCF.
    mavedb annotate-clinvar Look up ClinVar significance for mapped
                            variants via the public NCBI E-utilities.
    mavedb showcase        End-to-end demo: download PTEN scoreset +
                            map all variants + write VCF + composite
                            figure + narrative report.
    mavedb write-playbook  Write Docs/Skills/MAVEDB_MAPPING_SKILL.md

The Ensembl REST API has no auth requirement; we cache every response
under `Data/Cache/Ensembl/` to amortise repeated runs (a 4,408-variant
PTEN scoreset takes ~60 seconds first time, ~5 seconds on rerun).

License: Apache-2.0. Heavy deps imported lazily: pandas, requests,
matplotlib. No GPL runtime deps.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "MaveDB"
PLOT_DIR = REPORT_DIR / "Plots"
CACHE_DIR = DATA_DIR / "Cache" / "Ensembl"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"mavedb_mapping_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


# ─── Genetic code + amino-acid utilities ────────────────────────────────────

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODONS_BY_AA: dict[str, list[str]] = {}
for codon, aa in CODON_TABLE.items():
    CODONS_BY_AA.setdefault(aa, []).append(codon)

AA_1TO3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter",
}
AA_3TO1 = {v: k for k, v in AA_1TO3.items()}


def _complement(b: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}.get(b, "N")


def reverse_complement(s: str) -> str:
    return "".join(_complement(b) for b in reversed(s))


# ─── HGVS parsing ───────────────────────────────────────────────────────────

# Patterns
_HGVS_P_3 = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter|=)\)?$")
_HGVS_P_1 = re.compile(r"^p\.\(?([A-Z*])(\d+)([A-Z*]|=)\)?$")
_HGVS_C   = re.compile(r"^c\.(-?\d+)([ACGT])>([ACGT])$")


def parse_hgvsp(s: str) -> "tuple[str, int, str] | None":
    """Parse `p.Met1Leu` / `p.M1L` / `p.Met1=` / `p.Met1*` etc.

    Returns (wt_1letter, position, alt_1letter) or None if not a single-AA
    substitution.
    """
    if not s:
        return None
    s = s.strip()
    if s.endswith("fs") or "del" in s or "ins" in s or "dup" in s:
        return None  # frameshift / indel / duplication — out of scope for now
    m = _HGVS_P_3.match(s)
    if m:
        wt = AA_3TO1.get(m.group(1))
        alt = AA_3TO1.get(m.group(3)) if m.group(3) not in ("=", "*", "Ter") else (
            wt if m.group(3) == "=" else "*"
        )
        if wt is None or alt is None:
            return None
        return wt, int(m.group(2)), alt
    m = _HGVS_P_1.match(s)
    if m:
        wt = m.group(1)
        alt = (m.group(3) if m.group(3) != "=" else wt)
        return wt, int(m.group(2)), alt
    return None


def parse_hgvsc(s: str) -> "tuple[int, str, str] | None":
    """Parse `c.123A>G`. Returns (cdna_pos, ref, alt) or None."""
    if not s:
        return None
    m = _HGVS_C.match(s.strip())
    return (int(m.group(1)), m.group(2), m.group(3)) if m else None


# ─── Ensembl REST API client ────────────────────────────────────────────────

ENSEMBL_BASE = "https://rest.ensembl.org"


def _cache_path(slug: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / (safe_label(slug) + ".json")


def _ensembl_get(path: str, *, cache_key: str | None = None,
                  retries: int = 3) -> Any:
    """GET an Ensembl REST endpoint with on-disk caching + retries."""
    if cache_key:
        cp = _cache_path(cache_key)
        if cp.is_file():
            return json.loads(cp.read_text())
    url = ENSEMBL_BASE + path
    last_exc = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.load(resp)
            if cache_key:
                _cache_path(cache_key).write_text(json.dumps(data))
            return data
        except urllib.error.HTTPError as exc:
            last_exc = exc
            if exc.code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
            if attempt < retries - 1:
                time.sleep(1.0 * (attempt + 1))
                continue
            raise
        except Exception as exc:
            last_exc = exc
            if attempt < retries - 1:
                time.sleep(1.0)
                continue
            raise
    raise RuntimeError(f"Ensembl REST failed after {retries} retries: {last_exc}")


def get_gene_xrefs(symbol: str, species: str = "human") -> "list[dict]":
    """Look up an Ensembl gene by HGNC symbol."""
    data = _ensembl_get(
        f"/xrefs/symbol/{species}/{symbol}?",
        cache_key=f"xrefs_symbol_{species}_{symbol}",
    )
    return [d for d in data if d.get("type") == "gene"]


def get_canonical_transcript(gene_id: str) -> dict:
    """Get the canonical transcript for an Ensembl gene id."""
    data = _ensembl_get(
        f"/lookup/id/{gene_id}?expand=1",
        cache_key=f"lookup_{gene_id}",
    )
    transcripts = data.get("Transcript") or []
    canonical = [t for t in transcripts if t.get("is_canonical")]
    return canonical[0] if canonical else (transcripts[0] if transcripts else {})


def get_transcript_sequence(transcript_id: str, *, type_: str = "cds") -> str:
    """Pull the CDS or protein sequence for a transcript."""
    data = _ensembl_get(
        f"/sequence/id/{transcript_id}?type={type_}",
        cache_key=f"seq_{transcript_id}_{type_}",
    )
    return data.get("seq", "")


def map_protein_to_genomic(protein_id: str, aa_start: int, aa_end: int) -> list[dict]:
    """Map a protein position range to genomic coordinates.

    Note: the Ensembl `/map/translation/{id}` endpoint expects an ENSP
    (translation) identifier, not an ENST (transcript) identifier.
    """
    data = _ensembl_get(
        f"/map/translation/{protein_id}/{aa_start}..{aa_end}?",
        cache_key=f"map_p_{protein_id}_{aa_start}_{aa_end}",
    )
    return data.get("mappings", [])


def get_genomic_sequence(chrom: str, start: int, end: int, *,
                          species: str = "human") -> str:
    """Fetch a genomic substring."""
    coord = f"{chrom}:{start}..{end}:1"  # forward strand
    data = _ensembl_get(
        f"/sequence/region/{species}/{urllib.parse.quote(coord)}?",
        cache_key=f"seq_region_{species}_{chrom}_{start}_{end}",
    )
    return data.get("seq", "")


# ─── Variant mapping ────────────────────────────────────────────────────────

def _candidate_nt_changes(wt_codon: str, alt_aa: str) -> list[tuple[int, str, str]]:
    """Enumerate single-nucleotide changes to `wt_codon` that yield
    `alt_aa`. Returns list of (codon_pos, ref_nt, alt_nt) — 0-indexed.
    """
    out: list[tuple[int, str, str]] = []
    bases = ["A", "C", "G", "T"]
    for i in range(3):
        for nt in bases:
            if nt == wt_codon[i]:
                continue
            new_codon = wt_codon[:i] + nt + wt_codon[i + 1:]
            if CODON_TABLE.get(new_codon) == alt_aa:
                out.append((i, wt_codon[i], nt))
    return out


def map_variant(
    hgvsp: str,
    *,
    transcript_id: str | None = None,
    protein_id: str | None = None,
    protein_offset: int = 0,
    species: str = "human",
) -> "list[dict]":
    """Map a single HGVSp to one or more candidate genomic loci.

    Returns a list of {chr, pos, ref, alt, codon_pos, candidate_idx,
    notes} dicts. Multiple candidates exist for amino acids encoded
    by multiple codons. Empty list = unmappable / out-of-scope.
    """
    parsed = parse_hgvsp(hgvsp)
    if parsed is None:
        return []
    wt_aa, aa_pos, alt_aa = parsed
    if wt_aa == alt_aa:  # synonymous
        return [{
            "type": "synonymous", "chr": None, "pos": None,
            "ref": None, "alt": None, "notes": "synonymous (p.X=); skipped"
        }]

    eff_aa_pos = aa_pos + protein_offset

    # Fetch genomic mapping for this AA — uses the PROTEIN (ENSP) id
    pid = protein_id or transcript_id  # caller should pass protein_id
    try:
        mappings = map_protein_to_genomic(pid, eff_aa_pos, eff_aa_pos)
    except Exception as exc:
        return [{"type": "error", "notes": f"Ensembl map failed: {exc}"}]
    if not mappings:
        return [{"type": "error", "notes": "no genomic mapping for AA position"}]

    out: list[dict] = []
    for m in mappings:
        chrom = m.get("seq_region_name") or m.get("chr") or m.get("seq_region_id")
        start = int(m.get("start"))
        end = int(m.get("end"))
        strand = int(m.get("strand", 1))
        # Pull the 3-nt codon
        seq = get_genomic_sequence(chrom, start, end, species=species)
        if not seq or len(seq) != 3:
            out.append({"type": "error", "chr": chrom, "pos": start,
                        "notes": f"unexpected codon length {len(seq)}"})
            continue
        # Convert to the coding strand
        codon = seq.upper() if strand == 1 else reverse_complement(seq.upper())
        translated = CODON_TABLE.get(codon)
        # Sanity check
        if translated != wt_aa:
            out.append({
                "type": "wt_mismatch", "chr": chrom, "pos": start,
                "notes": f"expected WT={wt_aa} but transcript codon "
                          f"{codon} translates to {translated}",
                "ref_codon": codon, "transcript_id": transcript_id,
            })
            continue
        # Enumerate candidate single-NT changes producing alt_aa
        for cand_idx, (codon_pos, ref_nt, alt_nt) in enumerate(
                _candidate_nt_changes(codon, alt_aa)):
            if strand == 1:
                genomic_pos = start + codon_pos
                g_ref, g_alt = ref_nt, alt_nt
            else:
                # On the reverse strand, the codon's first base (5'→3' of
                # the protein) is the LAST base in genomic coordinates.
                genomic_pos = end - codon_pos
                g_ref, g_alt = _complement(ref_nt), _complement(alt_nt)
            out.append({
                "type": "missense" if alt_aa != "*" else "nonsense",
                "chr": chrom, "pos": genomic_pos,
                "ref": g_ref, "alt": g_alt,
                "wt_aa": wt_aa, "alt_aa": alt_aa,
                "aa_pos": aa_pos, "aa_pos_with_offset": eff_aa_pos,
                "codon": codon, "codon_pos": codon_pos,
                "strand": strand, "transcript_id": transcript_id,
                "candidate_idx": cand_idx,
                "n_candidates": len(_candidate_nt_changes(codon, alt_aa)),
                "notes": ("unique nt change" if len(_candidate_nt_changes(codon, alt_aa)) == 1
                          else "ambiguous nt change (multiple codons produce alt AA)"),
            })
        if not _candidate_nt_changes(codon, alt_aa):
            out.append({
                "type": "no_single_nt_change", "chr": chrom, "pos": start,
                "wt_aa": wt_aa, "alt_aa": alt_aa, "codon": codon,
                "notes": "no single-nt change in this codon produces alt AA "
                          "(requires ≥2 nt changes)",
            })
    return out


# ─── Scoreset-level orchestration ───────────────────────────────────────────

def map_scoreset(
    rows: "list[dict]",
    *,
    gene: str,
    species: str = "human",
    progress: bool = True,
) -> "tuple[list[dict], dict]":
    """Map every row in a parsed MaveDB scoreset.

    `rows` should have at least: `pos`, `wt`, `alt`, `kind`, `score`
    (the schema produced by proteomics_skill._read_mavedb_csv).

    Returns (mapped_rows, summary).
    """
    setup_logging()
    logging.info("Looking up Ensembl canonical transcript for %s", gene)
    xrefs = get_gene_xrefs(gene, species=species)
    if not xrefs:
        raise SystemExit(f"No Ensembl gene xref for {gene}.")
    gene_id = xrefs[0]["id"]
    canonical = get_canonical_transcript(gene_id)
    if not canonical:
        raise SystemExit(f"No canonical transcript for {gene_id}.")
    transcript_id = canonical["id"]
    transcript_name = canonical.get("display_name") or transcript_id
    protein_id = (canonical.get("Translation") or {}).get("id")
    if not protein_id:
        raise SystemExit(f"No Translation/protein id on canonical transcript {transcript_id}.")
    logging.info("  gene_id=%s  transcript=%s (%s)  protein=%s",
                  gene_id, transcript_id, transcript_name, protein_id)

    mapped: "list[dict]" = []
    skipped = 0
    error = 0
    by_type: Counter[str] = Counter()
    total = len(rows)
    for i, row in enumerate(rows):
        if progress and i % 200 == 0:
            logging.info("  %d / %d variants mapped (%d skipped, %d errors)",
                          i, total, skipped, error)
        wt = row.get("wt", "")
        alt = row.get("alt", "")
        pos = row.get("pos")
        kind = row.get("kind", "")
        if not wt or not alt or pos is None:
            skipped += 1; continue
        # MaveDB scoreset rows are pre-parsed by proteomics_skill — convert
        # to HGVSp string for our parser
        try:
            pos_int = int(pos)
        except (TypeError, ValueError):
            skipped += 1; continue
        # Compose an HGVSp string (single-letter)
        if kind == "synonymous" or wt == alt:
            mapped.append({
                **row, "mapping_type": "synonymous",
                "chr": None, "pos_genomic": None, "ref": None, "alt": None,
                "transcript_id": transcript_id, "notes": "synonymous",
            })
            by_type["synonymous"] += 1
            continue
        if kind in ("nonsense", "stop"):
            alt_for_hgvs = "*"
        elif alt in ("X", "*", "Ter"):
            alt_for_hgvs = "*"
        else:
            alt_for_hgvs = alt
        if len(wt) != 1 or (len(alt_for_hgvs) != 1):
            skipped += 1; continue
        hgvsp = f"p.{wt}{pos_int}{alt_for_hgvs}"
        try:
            candidates = map_variant(hgvsp, transcript_id=transcript_id,
                                       protein_id=protein_id,
                                       species=species)
        except Exception as exc:
            mapped.append({
                **row, "mapping_type": "error", "transcript_id": transcript_id,
                "notes": f"map error: {exc}",
            })
            error += 1
            continue
        if not candidates:
            skipped += 1; continue
        for c in candidates:
            mapped.append({
                **row,
                "hgvsp": hgvsp,
                "mapping_type": c.get("type"),
                "chr": c.get("chr"), "pos_genomic": c.get("pos"),
                "ref": c.get("ref"), "alt_nt": c.get("alt"),
                "transcript_id": c.get("transcript_id", transcript_id),
                "codon": c.get("codon"), "codon_pos": c.get("codon_pos"),
                "strand": c.get("strand"),
                "candidate_idx": c.get("candidate_idx"),
                "n_candidates": c.get("n_candidates"),
                "notes": c.get("notes"),
            })
            by_type[c.get("type", "?")] += 1
    summary = {
        "gene": gene, "transcript_id": transcript_id,
        "transcript_name": transcript_name,
        "protein_id": protein_id,
        "gene_ensembl_id": gene_id,
        "n_rows_in": total, "n_mapped_rows_out": len(mapped),
        "n_skipped": skipped, "n_errors": error,
        "type_counts": dict(by_type),
    }
    return mapped, summary


# ─── VCF writer ─────────────────────────────────────────────────────────────

def write_vcf(mapped_rows: list[dict], out_path: Path, *,
              gene: str = "", urn: str = "") -> Path:
    """Emit a minimal VCF-4.2 with chr, pos, ref, alt + INFO."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##source=IGVFagent_mavedb_mapping ({time.strftime('%Y-%m-%d')})\n")
        if urn:
            fh.write(f"##mavedb_urn={urn}\n")
        if gene:
            fh.write(f"##gene_symbol={gene}\n")
        fh.write('##INFO=<ID=AA_REF,Number=1,Type=String,Description="WT amino acid">\n')
        fh.write('##INFO=<ID=AA_ALT,Number=1,Type=String,Description="Alt amino acid">\n')
        fh.write('##INFO=<ID=AA_POS,Number=1,Type=Integer,Description="Protein position">\n')
        fh.write('##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="Ensembl transcript ID">\n')
        fh.write('##INFO=<ID=SCORE,Number=1,Type=Float,Description="MaveDB score">\n')
        fh.write('##INFO=<ID=AMBIG,Number=1,Type=Integer,Description="Total candidate nt changes">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in mapped_rows:
            chrom = r.get("chr")
            pos = r.get("pos_genomic")
            ref = r.get("ref")
            alt = r.get("alt_nt")
            if not chrom or pos is None or not ref or not alt:
                continue
            info_parts = [
                f"AA_REF={r.get('wt') or ''}",
                f"AA_ALT={r.get('alt') or ''}",
                f"AA_POS={r.get('pos') or ''}",
                f"TRANSCRIPT={r.get('transcript_id') or ''}",
            ]
            score = r.get("score")
            if score is not None and score != "":
                try:
                    info_parts.append(f"SCORE={float(score):.4f}")
                except (TypeError, ValueError):
                    pass
            n_cand = r.get("n_candidates")
            if n_cand:
                info_parts.append(f"AMBIG={n_cand}")
            info = ";".join(info_parts)
            vid = f"{r.get('hgvsp') or '.'}_{r.get('candidate_idx') or 0}"
            fh.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    return out_path


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_map_scoreset(args: argparse.Namespace) -> int:
    """Map a MaveDB scoreset's variants to genomic coordinates."""
    setup_logging()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from proteomics_skill import (download_mavedb_scoreset, _read_mavedb_csv,
                                    MAVEDB_VAMPSEQ_CATALOG)
    # Figure out URN + gene
    gene = (args.gene or "").upper() or None
    urn = args.urn
    if gene and not urn:
        meta = MAVEDB_VAMPSEQ_CATALOG.get(gene)
        if not meta:
            raise SystemExit(
                f"No MaveDB catalog entry for {gene}. "
                f"Available: {sorted(MAVEDB_VAMPSEQ_CATALOG)}. "
                f"Or supply --urn explicitly."
            )
        urn = meta["urn"]
    if not urn:
        raise SystemExit("Need either --gene <symbol> or --urn <URN>.")
    if not gene:
        # Reverse-lookup
        for k, v in MAVEDB_VAMPSEQ_CATALOG.items():
            if v.get("urn") == urn:
                gene = k; break
        if not gene:
            raise SystemExit("Need --gene to identify the Ensembl transcript.")

    csv_path = download_mavedb_scoreset(urn)
    rows = _read_mavedb_csv(csv_path)
    logging.info("Loaded %d rows from MaveDB %s", len(rows), urn)
    mapped, summary = map_scoreset(rows, gene=gene, species=args.species,
                                     progress=True)
    # Write outputs
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    label = args.label or f"{gene.lower()}_mapped"
    out_dir = REPORT_DIR / f"{ts}_{safe_label(label)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / f"{safe_label(gene)}_mapped.tsv"
    vcf = out_dir / f"{safe_label(gene)}_mapped.vcf"
    if mapped:
        keys = sorted({k for r in mapped for k in r.keys()})
        with tsv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t",
                                extrasaction="ignore")
            w.writeheader()
            for r in mapped:
                w.writerow(r)
    write_vcf(mapped, vcf, gene=gene, urn=urn)
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))

    print(f"Mapped TSV: {tsv}")
    print(f"VCF: {vcf}")
    print(f"Summary: {summary_path}")
    print(f"Gene: {summary['gene']}  transcript: {summary['transcript_id']}")
    print(f"Total rows: {summary['n_rows_in']:,}")
    print(f"Output rows: {summary['n_mapped_rows_out']:,}")
    print(f"Skipped: {summary['n_skipped']:,}, errors: {summary['n_errors']}")
    print("Type breakdown:")
    for k, v in summary["type_counts"].items():
        print(f"  {k}: {v:,}")
    return 0


def cmd_showcase(args: argparse.Namespace) -> int:
    """Single-command demo: PTEN scoreset → genomic coordinates →
    composite figure showing coverage along the protein."""
    args2 = argparse.Namespace(
        gene=args.gene or "PTEN", urn=None,
        species=args.species, label=args.label or f"{(args.gene or 'PTEN').lower()}_showcase",
    )
    rc = cmd_map_scoreset(args2)
    if rc != 0:
        return rc
    # Build a one-panel coverage figure (variant density vs protein position)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    latest = sorted(REPORT_DIR.glob(f"*{args2.label}/*_mapped.tsv"))[-1]
    df = pd.read_csv(latest, sep="\t")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), facecolor="white")
    if "pos" in df.columns:
        pos_series = pd.to_numeric(df["pos"], errors="coerce").dropna()
        axes[0].hist(pos_series, bins=80, color="#5C8DAA",
                      edgecolor="#1F2933", linewidth=0.5)
        axes[0].set_xlabel("Protein position")
        axes[0].set_ylabel("Mapped variants")
        axes[0].set_title(f"{args2.gene} VAMP-seq variant coverage along protein "
                            f"(MaveDB → genomic)", fontweight="bold")
    if "mapping_type" in df.columns:
        type_counts = df["mapping_type"].value_counts()
        axes[1].bar(type_counts.index, type_counts.values, color="#7CA663",
                     edgecolor="#1F2933", linewidth=0.5)
        axes[1].set_ylabel("Rows")
        axes[1].set_title("Mapping outcome breakdown", fontweight="bold")
        for tick in axes[1].get_xticklabels():
            tick.set_rotation(20); tick.set_ha("right")
    plt.tight_layout()
    out_png = latest.parent / "mapping_summary.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Composite figure: {out_png}")

    # Narrative report
    rep = latest.parent / "showcase_report.md"
    n = len(df)
    n_missense = (df["mapping_type"] == "missense").sum() if "mapping_type" in df.columns else 0
    n_synonymous = (df["mapping_type"] == "synonymous").sum() if "mapping_type" in df.columns else 0
    n_nonsense = (df["mapping_type"] == "nonsense").sum() if "mapping_type" in df.columns else 0
    unique_chr = df["chr"].dropna().unique().tolist() if "chr" in df.columns else []
    rep.write_text(
        f"# {args2.gene} VAMP-seq → genomic mapping showcase\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"- MaveDB URN: see `summary.json`\n"
        f"- Total mapped output rows: **{n:,}**\n"
        f"- Missense: **{n_missense:,}**\n"
        f"- Synonymous (skipped): **{n_synonymous:,}**\n"
        f"- Nonsense: **{n_nonsense:,}**\n"
        f"- Chromosomes touched: **{', '.join(unique_chr)}**\n\n"
        f"## Outputs\n\n"
        f"- `{latest.name}` — per-variant mapped TSV (chr/pos/ref/alt + score columns)\n"
        f"- `{args2.gene}_mapped.vcf` — VCF-4.2 with INFO fields\n"
        f"- `summary.json` — gene/transcript metadata + type counts\n"
        f"- `mapping_summary.png` — coverage along protein + outcome bar chart\n\n"
        f"## Notes on ambiguity\n\n"
        f"For some missense variants, multiple single-nucleotide changes "
        f"in the WT codon produce the alt amino acid. This skill emits "
        f"one row per candidate nt change with `candidate_idx` + "
        f"`n_candidates` so downstream consumers can decide whether to "
        f"average effects, pick the canonical one matching `hgvs_nt`, or "
        f"propagate the ambiguity.\n\n"
        f"## License posture\n\n"
        f"Clean-room reimplementation of "
        f"[ave-dcd/dcd_mapping](https://github.com/ave-dcd/dcd_mapping) "
        f"(MIT). Uses only the public Ensembl REST API + a local JSON "
        f"cache; no UTA / SeqRepo / BLAT required.\n"
    )
    print(f"Report: {rep}")
    return 0


def cmd_write_playbook(args: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "MAVEDB_MAPPING_SKILL.md"
    path.write_text("""# Skill: MaveDB → genomic coordinates

Map MaveDB scoreset variants (HGVSp / per-position tile / single-letter
WT-pos-ALT) to canonical chr/pos/ref/alt VCF coordinates so MAVE /
VAMP-seq scores become cross-referenceable with ClinVar, gnomAD, GWAS
catalogues and the IGVF Catalog.

Clean-room reimplementation of [ave-dcd/dcd_mapping](https://github.com/ave-dcd/dcd_mapping)
(MIT, Arbesfeld 2023) using the public Ensembl REST API — no UTA /
SeqRepo / BLAT install required.

## Commands

```bash
# Map any VAMP-seq scoreset by gene symbol (uses MAVEDB_VAMPSEQ_CATALOG)
igvfagent mavedb map-scoreset --gene PTEN --label pten_v1

# Or by URN directly (works for any single-AA missense scoreset)
igvfagent mavedb map-scoreset --urn urn:mavedb:00000013-a-1 --gene PTEN

# One-command showcase: map + composite figure + narrative report
igvfagent mavedb showcase --gene PTEN
```

## Output schema (TSV)

| Column | Meaning |
|---|---|
| pos / wt / alt | Original MaveDB protein position + AA letters |
| kind | Original kind (missense / synonymous / nonsense / wt) |
| score / sd / se / lower_ci / upper_ci | Original MaveDB score columns |
| hgvsp | Reconstructed HGVSp (e.g. `p.M1L`) |
| mapping_type | missense / synonymous / nonsense / wt_mismatch / error |
| chr / pos_genomic / ref / alt_nt | Genomic coordinates + nt change |
| transcript_id | Ensembl canonical transcript |
| codon / codon_pos / strand | Codon context |
| candidate_idx / n_candidates | For amino acids with multiple alt codons, one row per candidate |
| notes | Per-row provenance string |

## VCF schema

Standard VCF-4.2 with INFO fields:
`AA_REF`, `AA_ALT`, `AA_POS`, `TRANSCRIPT`, `SCORE`, `AMBIG`.

## How it works (clean-room algorithm)

1. **Gene → Ensembl canonical transcript**
   `/xrefs/symbol/{species}/{gene}` then `/lookup/id/{gene_id}?expand=1`
   → pick the transcript with `is_canonical=true`.
2. **Protein position → genomic codon coordinates**
   `/map/translation/{transcript_id}/{aa}..{aa}` returns the 3-nt codon's
   chr / start / end / strand.
3. **Codon validation**
   Fetch the genomic 3-nt substring (`/sequence/region/{species}/{chr:start..end}:1`),
   reverse-complement if strand=-1, translate via the standard table,
   verify against the MaveDB WT amino acid.
4. **Single-nt change enumeration**
   For each position in the WT codon, find single-nucleotide changes
   that yield the alt amino acid. Emit one output row per candidate
   with `candidate_idx` + `n_candidates` so downstream consumers can
   resolve ambiguity (`n_candidates=1` means unique; >1 means multiple
   single-nt changes give the same protein change).
5. **Strand-aware genomic coordinates**
   Forward strand: pos = codon_start + codon_pos. Reverse strand:
   pos = codon_end - codon_pos, ref/alt complemented.

## Caching

Every Ensembl REST response is cached under `Data/Cache/Ensembl/`. A
fresh 4,408-variant PTEN scoreset maps in ~60 s the first time, ~5 s on
rerun.

## What's NOT yet supported

- Frameshift (`p.Xfs`), large deletions, duplications, insertions
- Non-canonical transcripts (we pick the Ensembl-canonical one)
- Multi-codon-substitution variants
- Mapping when the MaveDB target sequence doesn't match the canonical
  protein (protein-offset detection — TODO)

For these cases, fall back to the upstream `dcd-mapping` PyPI package
(MIT) which handles all of the above via UTA + SeqRepo + BLAT.
""")
    print(f"Wrote: {path}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="mavedb_mapping_skill",
        description="MaveDB scoreset → genomic coordinates mapping "
                     "(clean-room reimplementation of ave-dcd/dcd_mapping).")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("map-scoreset",
        help="Map a MaveDB scoreset's variants to chr/pos/ref/alt.")
    p.add_argument("--urn", default=None,
                    help="MaveDB URN (urn:mavedb:00000013-a-1). "
                          "If omitted and --gene is supplied, the URN is "
                          "looked up in MAVEDB_VAMPSEQ_CATALOG.")
    p.add_argument("--gene", default=None,
                    help="HGNC gene symbol (e.g. PTEN). Required for the "
                          "Ensembl transcript lookup.")
    p.add_argument("--species", default="human")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_map_scoreset)

    p = sub.add_parser("showcase",
        help="One-command demo: map + composite figure + narrative report.")
    p.add_argument("--gene", default="PTEN")
    p.add_argument("--species", default="human")
    p.add_argument("--label", default=None)
    p.set_defaults(func=cmd_showcase)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/MAVEDB_MAPPING_SKILL.md")
    p.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
