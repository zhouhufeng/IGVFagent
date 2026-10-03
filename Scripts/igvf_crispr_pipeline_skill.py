#!/usr/bin/env python3
"""IGVF CRISPR Perturb-seq pipeline, early Nextflow version: seqspec parsing, guide references, kb mapping, AnnData QC, MuData assembly (port of IGVF-CRISPR/IGVF_CRISPR_Pipeline).

Port of https://github.com/IGVF-CRISPR/IGVF_CRISPR_Pipeline (No LICENSE file;
Nextflow DSL2 + Python helper scripts, pinned at commit 6704a0f, 2024-07-18).
Every module and workflow was read: Modules/{downloadReference, mappingscRNA,
seqSpecParser} and Workflows/{main_pipeline_guide, main_pipeline_rna,
MergeTwoWorkflows, process_anndata}.  Because the upstream declares no
licence, only the METHODS (definitions, thresholds, column names, output
schemas, command lines) were re-derived here in Python; no code was copied.
Relationship: port.

Definitions, exactly as upstream computes them
  seqspec reads        `seqspec info` -> the read_id of every sequence_spec
                       !Read whose modality equals the requested modality,
                       in file order, comma-joined (parsing_guide_metadata.py
                       process_reads).  Here the YAML is read directly (a
                       small parser for the seqspec subset; the `seqspec`
                       binary is used instead when on PATH, --engine).
  kb technology (-x)   `seqspec index -m <mod> -t kb -r <reads>`: for each
                       read (index = position in the read list), the leaves of
                       the modality's library_spec after the read's primer
                       (strand pos) or before it, reversed (strand neg); the
                       leaves are laid end to end by max_len and clipped to
                       [0, read max_len); barcode / umi / feature
                       (cdna, gdna, protein, tag, sgrna_target) coordinates
                       are written "file,start,stop", grouped "BC:UMI:FEATURE",
                       a missing group as -1,-1,-1.  rna.yml -> 0,0,16:0,16,28:1,0,150.
  onlist               the first line of the YAML starting with "filename"
                       gives the whitelist; barcode_whitelist = join(directory, it).
  parsed seqspec       <modality>_parsed_seqSpec.txt, TSV with columns
                       modality, representation, barcode_whitelist, seqspec_dir.
  extract              the -x string = concatenation (no separator) of the
                       representation column (extract_parsed_seqspec.py).
  guide features       guide table (xlsx; tsv/csv also accepted) ->
                       guide_features.txt, tab-separated, no header, columns
                       sgRNA_sequences, sgRNA_ID (guide_features.py).
  process reads        strip all whitespace from the ';'-separated read-1 and
                       read-2 lists, zip them and emit "dir/r1 dir/r2 ..."
                       (process_reads.py).
  kb ref               transcriptome: kb ref -d <species> -i transcriptome_index.idx
                       -g transcriptome_t2g.txt --kallisto <k_bin>;  guides (kite):
                       kb ref -i guide_index.idx -f1 <fasta> -g t2guide.txt
                       --kallisto <k_bin> --workflow kite guide_features.txt.
  kb count             kb count -i <idx> -g <t2g> --verbose -w <whitelist> --h5ad
                       --kallisto <k_bin> -x <tech> -o ks_{transcripts,guide}_out
                       -t <cpus> <fastqs> --overwrite -m 48G   (cpus 6; the
                       MergeTwoWorkflows variant uses -m 30G).
  preprocess           var_names <- first column of the gene-names file (count
                       must equal n_vars, else error); knee plot of log1p(total
                       UMI) per barcode, sorted descending (knee_plot_rna.png);
                       obs.batch_number = 1; sc.pp.filter_cells(min_genes=100);
                       sc.pp.filter_genes(min_cells=3); var.mt = name startswith
                       "MT-" if reference == "human" else "Mt-"; var.ribo = startswith
                       RPS/RPL; sc.pp.calculate_qc_metrics(qc_vars=[mt, ribo],
                       log1p=True); violin of n_genes_by_counts / total_counts /
                       pct_counts_mt and scatter total_counts vs n_genes_by_counts
                       coloured by pct_counts_mt; -> filtered_anndata.h5ad.
  create MuData        guide var_names <- sgRNA_ID + "|" + sgRNA_sequences;
                       obs.number_of_nonzero_guides = count of guides with X > 0;
                       obs.batch_number = 1; guide knee plot (knee_plot_guide.png);
                       barcodes = intersection of RNA and guide obs_names;
                       MuData {'transcripts': rna[bc], 'guides': guide[bc]} -> mudata.h5mu.

Deliberate deviations (each reversible with --upstream-compat where it matters)
  * create-mdata: upstream assigns the "ID|sequence" names to the guide matrix
    by POSITION; the port matches them by sgRNA_ID when the matrix's var_names
    are guide IDs (kb/kite output) and falls back to position otherwise.  The
    upstream intersects barcodes through a Python set (arbitrary order); the
    port keeps the RNA obs order.  --upstream-compat: positional names, set order.
  * kb-ref --workflow kite: upstream passes the downloaded hg38 genome as -f1,
    but in kite mode -f1 is where kb WRITES the mismatch FASTA, so the genome
    download is unused (and overwritten).  The port writes -f1 guide_mismatch.fa;
    --upstream-compat passes --genome as -f1.
  * The multiseq example uses region_type cell_barcode; seqspec's kb formatter
    only recognises "barcode".  The port treats cell_barcode as a barcode and
    reports it (--strict-region-types reproduces the binary).
  * The multiseq example declares a 150-nt read 2 over a 52-nt insert, so read 2
    runs through the UMI and barcode too (0,0,16,1,20,36:0,16,26,1,8,20:1,0,8);
    the port reproduces the seqspec layout rule rather than second-guessing it.
  * preprocess: scanpy's default percent_top (50,100,200,500) raises on < 500
    genes; the tiers that fit are kept (identical output on real data).
  * The Workflows/*/3M-february-2018.txt.gz copies in the upstream repo are
    truncated gzips (~2.9 M of 6.79 M barcodes); count-features warns and uses
    the readable part (the Modules/seqSpecParser copy is complete).
  * count-features is an addition, not an upstream step: a pure-Python stand-in
    for `kb count --workflow kite` (exact or unique Hamming-1 feature match as
    kite's mismatch index does, whitelist exact or unique Hamming-1 barcode
    correction, distinct-UMI counting) so the guide arm runs without kallisto.

Subcommands
  parse-seqspec     seqspec YAML -> <modality>_parsed_seqSpec.txt (+ reads, -x).
  extract-seqspec   parsed-seqspec TSV -> the kb -x string (stdout + JSON).
  guide-features    guide table -> guide_features.txt.
  process-reads     ';'-separated read-1 / read-2 lists -> the kb fastq argument.
  download-genome   wget of the genome FASTA used by the upstream kite step
                    (prints the command; --execute to run it; hg38 is ~1 GB).
  kb-ref            kb ref for the transcriptome (-d) or guides (kite); runs kb
                    when on PATH, else prints the command and exits 0.
  kb-count          kb count for the rna or guide modality, -x from a parsed
                    seqspec; same binary policy.
  count-features    pure-Python guide counting from paired FASTQs + tech string
                    -> counts_unfiltered/adata.h5ad.
  preprocess        RNA AnnData QC (process_anndata/PreprocessAnnData).
  create-mdata      RNA + guide AnnData + guide metadata -> MuData
                    (process_anndata/CreateMuData).
  run               the chain: seqspec parse -> guide features -> guide counts
                    (kb or count-features) -> preprocess -> create-mdata.
  selftest          upstream example seqspecs + synthetic FASTQs / AnnData with
                    planted signal; every subcommand asserted.

Output: Docs/IGVFCRISPRPipeline/<timestamp>_<label>/.  pandas + numpy required;
anndata + scipy for count/MuData steps; scanpy for preprocess; mudata optional
(an h5mu is written with h5py + anndata otherwise); matplotlib optional.

Usage:
    igvfagent igvf-crispr-pipeline parse-seqspec --yaml guide.yml --modality guide --directory seqspec_dir
    igvfagent igvf-crispr-pipeline extract-seqspec --file guide_parsed_seqSpec.txt
    igvfagent igvf-crispr-pipeline guide-features --guide-table gasperini_tss.xlsx
    igvfagent igvf-crispr-pipeline process-reads --dir fastqs --reads1 "a_R1.fq.gz; b_R1.fq.gz" --reads2 "a_R2.fq.gz; b_R2.fq.gz"
    igvfagent igvf-crispr-pipeline kb-ref --mode kite --guide-features guide_features.txt
    igvfagent igvf-crispr-pipeline kb-count --modality guide --index guide_index.idx --t2g t2guide.txt --parsed-seqspec guide_parsed_seqSpec.txt --fastqs R1.fq R2.fq
    igvfagent igvf-crispr-pipeline count-features --seqspec guide.yml --guide-table guides.xlsx --fastqs guide_R1.fq guide_R2.fq
    igvfagent igvf-crispr-pipeline preprocess --adata-rna adata_rna.h5ad --gene-names cells_x_genes.genes.names.txt --reference human
    igvfagent igvf-crispr-pipeline create-mdata --adata-rna filtered_anndata.h5ad --adata-guide adata_guide.h5ad --guide-metadata guides.xlsx
    igvfagent igvf-crispr-pipeline run --adata-rna adata_rna.h5ad --gene-names genes.txt --guide-seqspec guide.yml --guide-table guides.xlsx --guide-fastqs R1.fq R2.fq
    igvfagent igvf-crispr-pipeline selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "IGVFCRISPRPipeline"

UPSTREAM_REPO = "IGVF-CRISPR/IGVF_CRISPR_Pipeline"
UPSTREAM_COMMIT = "6704a0f16e4de35185db7e358268aabc75fe5256"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
BARCODE_TYPES = {"barcode"}
BARCODE_TYPES_LENIENT = {"barcode", "cell_barcode"}
UMI_TYPES = {"umi"}
FEATURE_TYPES = {"cdna", "gdna", "protein", "tag", "sgrna_target"}
PARSED_COLS = ["modality", "representation", "barcode_whitelist", "seqspec_dir"]
KB_OUT = {"rna": "ks_transcripts_out", "guide": "ks_guide_out"}

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"

log = logging.getLogger("igvf_crispr_pipeline")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"igvf_crispr_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    base = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d, k = base, 1
    while d.exists():
        k += 1
        d = Path(f"{base}_{k}")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


def _ad():
    try:
        import anndata as ad  # type: ignore
        return ad
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs anndata + scipy: pip install anndata scipy") from e


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "List[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def write_json(obj: Any, path: Path) -> Path:
    path.write_text(json.dumps(obj, indent=2, default=str))
    print(f"JSON: {path}")
    return path


def write_report(path: Path, title: str, sections: "List[str]") -> Path:
    head = [f"# {title}", "",
            f"Port of [{UPSTREAM_REPO}](https://github.com/{UPSTREAM_REPO}) @ `{UPSTREAM_COMMIT[:7]}` (no licence declared; methods re-derived).", ""]
    path.write_text("\n".join(head + sections) + "\n")
    print(f"Report: {path}")
    return path


def _open_text(path: "str | Path"):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


# ---------------------------------------------------------------------------
# seqspec YAML (subset parser: mappings, block lists, !Tags, scalars)
# ---------------------------------------------------------------------------

_KEY_RE = re.compile(r"^([A-Za-z0-9_\-. ]+?):(?:\s+(.*))?$")


def _scalar(tok: str) -> Any:
    t = tok.strip()
    if t.startswith("!"):  # a tag on a scalar line, e.g. "!Onlist" with nothing after it
        parts = t.split(None, 1)
        t = parts[1] if len(parts) > 1 else ""
    if t == "" or t in ("null", "~", "Null", "NULL"):
        return None
    if len(t) >= 2 and t[0] == t[-1] and t[0] in "'\"":
        inner = t[1:-1]
        return inner.replace("''", "'") if t[0] == "'" else inner
    if t in ("true", "True"):
        return True
    if t in ("false", "False"):
        return False
    if re.fullmatch(r"[-+]?\d+", t):
        return int(t)
    if re.fullmatch(r"[-+]?(\d+\.\d*|\.\d+|\d+)([eE][-+]?\d+)?", t):
        return float(t)
    return t


def _strip_comment(line: str) -> str:
    out, q = [], None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"" and (i == 0 or line[i - 1] in " :-"):
            q = ch
        elif ch == "#" and (i == 0 or line[i - 1] == " "):
            break
        out.append(ch)
    return "".join(out).rstrip()


def parse_seqspec_yaml(text: str) -> Any:
    """Parse the block-YAML subset seqspec files use (custom !Assay/!Read/!Region/!Onlist tags become plain dicts)."""
    lines: "List[List[Any]]" = []
    for raw in text.splitlines():
        s = _strip_comment(raw.replace("\t", "  "))
        if not s.strip() or s.strip() in ("---", "..."):
            continue
        lines.append([len(s) - len(s.lstrip(" ")), s.strip()])
    if lines and lines[0][1].startswith("!") and " " not in lines[0][1]:
        lines = lines[1:]  # document-level tag (!Assay)

    def is_item(t: str) -> bool:
        return t == "-" or t.startswith("- ")

    def node(i: int, ind: int) -> "Tuple[Any, int]":
        return (seq(i, ind) if is_item(lines[i][1]) else mapping(i, ind))

    def child(i: int, ind: int, allow_same_indent_list: bool) -> "Tuple[Any, int]":
        if i < len(lines):
            ni, nt = lines[i]
            if ni > ind:
                return node(i, ni)
            if allow_same_indent_list and ni == ind and is_item(nt):
                return seq(i, ind)
        return None, i

    def seq(i: int, ind: int) -> "Tuple[list, int]":
        out = []
        while i < len(lines) and lines[i][0] == ind and is_item(lines[i][1]):
            rest = lines[i][1][1:].lstrip(" ")
            if rest.startswith("!"):
                parts = rest.split(None, 1)
                rest = parts[1] if len(parts) > 1 else ""
            if rest == "":
                v, i = child(i + 1, ind, False)
                out.append(v)
            elif _KEY_RE.match(rest):
                lines[i] = [ind + len(lines[i][1]) - len(rest), rest]
                v, i = mapping(i, lines[i][0])
                out.append(v)
            else:
                out.append(_scalar(rest))
                i += 1
        return out, i

    def mapping(i: int, ind: int) -> "Tuple[dict, int]":
        out: "Dict[str, Any]" = {}
        while i < len(lines) and lines[i][0] == ind and not is_item(lines[i][1]):
            m = _KEY_RE.match(lines[i][1])
            if not m:
                raise ValueError(f"seqspec YAML: cannot parse line {lines[i][1]!r}")
            key, val = m.group(1).strip(), (m.group(2) or "").strip()
            if val.startswith("!") and " " not in val:
                val = ""
            if val == "":
                v, i = child(i + 1, ind, True)
                out[key] = v
            else:
                out[key] = _scalar(val)
                i += 1
        return out, i

    if not lines:
        return {}
    val, _ = node(0, lines[0][0])
    return val


def load_seqspec(path: "str | Path") -> dict:
    with _open_text(path) as fh:
        spec = parse_seqspec_yaml(fh.read())
    if not isinstance(spec, dict) or "sequence_spec" not in spec or "library_spec" not in spec:
        raise ValueError(f"{path}: not a seqspec (needs sequence_spec and library_spec)")
    return spec


def seqspec_reads(spec: dict, modality: str) -> "List[str]":
    """process_reads: read_ids of the sequence_spec entries of this modality, file order."""
    return [str(r.get("read_id")) for r in (spec.get("sequence_spec") or []) if str(r.get("modality")) == modality]


def _leaves(region: dict) -> "List[dict]":
    subs = region.get("regions") or []
    if not subs:
        return [region]
    out: "List[dict]" = []
    for r in subs:
        out.extend(_leaves(r))
    return out


def _libspec(spec: dict, modality: str) -> dict:
    for r in spec.get("library_spec") or []:
        if str(r.get("region_id")) == modality:
            return r
    raise ValueError(f"no library_spec region with region_id {modality!r}")


def read_region_coordinates(spec: dict, modality: str, read_id: str) -> "List[dict]":
    """Leaves covered by one read, as {region_id, region_type, start, stop} clipped to the read length."""
    read = next((r for r in spec["sequence_spec"] if str(r.get("read_id")) == read_id), None)
    if read is None:
        raise ValueError(f"read {read_id!r} not in sequence_spec")
    leaves = _leaves(_libspec(spec, modality))
    ids = [str(l.get("region_id")) for l in leaves]
    primer = str(read.get("primer_id"))
    if primer not in ids:
        raise ValueError(f"primer {primer!r} of read {read_id!r} is not a leaf of {modality!r}")
    idx = ids.index(primer)
    rgns = leaves[:idx][::-1] if str(read.get("strand", "pos")) == "neg" else leaves[idx + 1:]
    rlen = int(read.get("max_len") or 0)
    out, prev = [], 0
    for r in rgns:
        nxt = prev + int(r.get("max_len") or 0)
        start, stop = prev, min(nxt, rlen)
        prev = nxt
        if start >= rlen:
            break
        out.append({"region_id": str(r.get("region_id")), "region_type": str(r.get("region_type")), "start": start, "stop": stop})
    return out


def kb_technology(spec: dict, modality: str, reads: "Optional[List[str]]" = None, strict: bool = False) -> "Tuple[str, dict]":
    """seqspec index -t kb equivalent.  Returns (x_string, details)."""
    reads = reads if reads is not None else seqspec_reads(spec, modality)
    bc_types = BARCODE_TYPES if strict else BARCODE_TYPES_LENIENT
    bcs, umis, feats, per_read, lenient_hits = [], [], [], {}, []
    for fi, rid in enumerate(reads):
        coords = read_region_coordinates(spec, modality, rid)
        per_read[rid] = coords
        for c in coords:
            rt = c["region_type"].lower()
            s = f"{fi},{c['start']},{c['stop']}"
            if rt in bc_types:
                bcs.append(s)
                if rt not in BARCODE_TYPES:
                    lenient_hits.append(c["region_id"])
            elif rt in UMI_TYPES:
                umis.append(s)
            elif rt in FEATURE_TYPES:
                feats.append(s)
    x = ":".join(",".join(g) if g else "-1,-1,-1" for g in (bcs, umis, feats))
    return x, {"reads": reads, "coordinates": per_read, "cell_barcode_as_barcode": lenient_hits}


def extract_whitelist_file(yaml_path: "str | Path") -> "Optional[str]":
    """First line starting with 'filename' -> its value (upstream's line scan, not a YAML walk)."""
    with _open_text(yaml_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("filename"):
                parts = s.split(":", 1)
                if len(parts) > 1:
                    return parts[1].strip()
    return None


def _seqspec_binary_index(yaml_path: str, modality: str, reads: "List[str]") -> "Optional[str]":
    exe = shutil.which("seqspec")
    if not exe:
        return None
    cmd = [exe, "index", "-m", modality, "-t", "kb", "-r", ",".join(reads), yaml_path]
    log.info("run: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        log.warning("seqspec index failed: %s", res.stderr)
        return None
    return res.stdout.replace("\n", "")


def parse_seqspec_file(yaml_path: "str | Path", modalities: "List[str]", directory: "str | Path",
                       engine: str = "auto", strict: bool = False) -> "Tuple[Any, dict]":
    """parsing_guide_metadata.py main(): one row per modality."""
    pd = _pd()
    spec = load_seqspec(yaml_path)
    wl = extract_whitelist_file(yaml_path)
    rows, info = [], {}
    for m in modalities:
        reads = seqspec_reads(spec, m)
        x_py, det = kb_technology(spec, m, reads, strict=strict)
        x_bin = _seqspec_binary_index(str(yaml_path), m, reads) if engine in ("auto", "seqspec") else None
        if engine == "seqspec" and x_bin is None:
            raise SystemExit("--engine seqspec: the seqspec binary is not on PATH (pip install seqspec)")
        x = x_bin if x_bin is not None else x_py
        rows.append([m, x, os.path.join(str(directory), wl) if wl else "", str(directory)])
        info[m] = {"reads": reads, "representation": x, "python_representation": x_py, "seqspec_binary": x_bin,
                   "engine": "seqspec" if x_bin is not None else "python", "coordinates": det["coordinates"],
                   "cell_barcode_as_barcode": det["cell_barcode_as_barcode"], "onlist": wl}
    return pd.DataFrame(rows, columns=PARSED_COLS), {"assay_id": spec.get("assay_id"), "modalities": info}


def extract_parsed_seqspec(path: "str | Path") -> "Optional[str]":
    """extract_parsed_seqspec.py: join the representation column; None when the file is missing."""
    pd = _pd()
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path, sep="\t")
    return "".join(df["representation"].astype(str))


def parse_technology(x: str) -> "Dict[str, List[Tuple[int, int, int]]]":
    groups = x.strip().split(":")
    if len(groups) != 3:
        raise ValueError(f"technology string {x!r} is not BC:UMI:FEATURE")
    out = {}
    for name, g in zip(("barcode", "umi", "feature"), groups):
        nums = [int(v) for v in g.split(",") if v != ""]
        if len(nums) % 3:
            raise ValueError(f"technology group {g!r} is not a list of file,start,stop triples")
        trip = [tuple(nums[i:i + 3]) for i in range(0, len(nums), 3)]
        out[name] = [t for t in trip if t[0] >= 0]  # type: ignore[misc]
    return out  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# guide features, read lists
# ---------------------------------------------------------------------------

def read_guide_table(path: "str | Path"):
    pd = _pd()
    p = str(path)
    if p.endswith((".xlsx", ".xls")):
        df = pd.read_excel(p, header=0)
    elif p.endswith((".csv", ".csv.gz")):
        df = pd.read_csv(p)
    else:
        df = pd.read_csv(p, sep="\t")
    missing = {"sgRNA_ID", "sgRNA_sequences"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: guide table needs columns sgRNA_ID and sgRNA_sequences (missing {sorted(missing)})")
    return df


def write_guide_features(table, out: Path) -> Path:
    table[["sgRNA_sequences", "sgRNA_ID"]].to_csv(out, sep="\t", header=False, index=False)
    print(f"Wrote: {out}")
    return out


def read_guide_features(path: "str | Path") -> "List[Tuple[str, str]]":
    out = []
    with _open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                out.append((parts[0].strip().upper(), parts[1].strip()))
    return out


def pair_reads(directory: str, reads1: str, reads2: str) -> "List[str]":
    r1 = re.sub(r"\s", "", reads1).split(";")
    r2 = re.sub(r"\s", "", reads2).split(";")
    out = []
    for a, b in zip(r1, r2):
        out += [os.path.join(directory, a), os.path.join(directory, b)]
    return out


# ---------------------------------------------------------------------------
# binaries: kb ref / kb count / wget
# ---------------------------------------------------------------------------

def _kallisto_arg() -> str:
    return shutil.which("kallisto") or "$(which kallisto)"


def kb_ref_command(mode: str, species: str = "human", guide_features: str = "guide_features.txt", fasta: "Optional[str]" = None,
                   out_dir: "Optional[Path]" = None) -> "List[str]":
    o = (lambda n: str(out_dir / n)) if out_dir else (lambda n: n)
    if mode == "transcriptome":
        return ["kb", "ref", "-d", species, "-i", o("transcriptome_index.idx"), "-g", o("transcriptome_t2g.txt"), "--kallisto", _kallisto_arg()]
    return ["kb", "ref", "-i", o("guide_index.idx"), "-f1", fasta or o("guide_mismatch.fa"), "-g", o("t2guide.txt"),
            "--kallisto", _kallisto_arg(), "--workflow", "kite", guide_features]


def kb_count_command(index: str, t2g: str, whitelist: str, technology: str, out: str, fastqs: "List[str]",
                     threads: int = 6, memory: str = "48G") -> "List[str]":
    return ["kb", "count", "-i", index, "-g", t2g, "--verbose", "-w", whitelist, "--h5ad", "--kallisto", _kallisto_arg(),
            "-x", technology, "-o", out, "-t", str(threads)] + list(fastqs) + ["--overwrite", "-m", memory]


def run_or_print(cmd: "List[str]", binary: str, dry_run: bool) -> "Tuple[str, int]":
    """Run cmd when `binary` is on PATH (and not dry_run); otherwise print it.  Returns (status, returncode)."""
    text = " ".join(shlex.quote(c) if not c.startswith("$(") else c for c in cmd)
    print(f"Command: {text}")
    if dry_run or not shutil.which(binary):
        why = "--dry-run" if dry_run else f"`{binary}` is not on PATH"
        print(f"not executed ({why}); run the command above where {binary} is installed")
        return "printed", 0
    log.info("run: %s", text)
    res = subprocess.run(cmd)
    return ("ran" if res.returncode == 0 else "failed"), res.returncode


# ---------------------------------------------------------------------------
# count-features: pure-Python kite-style guide counting
# ---------------------------------------------------------------------------

_BASES = "ACGTN"


def hamming1(seq: str) -> "Iterable[str]":
    for i, ch in enumerate(seq):
        for b in _BASES:
            if b != ch:
                yield seq[:i] + b + seq[i + 1:]


def feature_lookup(features: "List[Tuple[str, str]]", mismatch: bool = True) -> "Tuple[Dict[str, int], Dict[str, int]]":
    """exact map sequence -> feature index, and the unique Hamming-1 map (collisions dropped, as kite does)."""
    exact = {s: i for i, (s, _) in enumerate(features)}
    one: "Dict[str, int]" = {}
    bad = set()
    if mismatch:
        for i, (s, _) in enumerate(features):
            for v in hamming1(s):
                if v in exact:
                    continue
                if v in one and one[v] != i:
                    bad.add(v)
                else:
                    one[v] = i
        for v in bad:
            one.pop(v, None)
    return exact, one


def _fastq_seqs(path: str) -> "Iterable[str]":
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()


def _cut(reads: "List[str]", spec: "List[Tuple[int, int, int]]") -> str:
    parts = []
    for f, s, e in spec:
        r = reads[f]
        parts.append(r[s:] if e == 0 else r[s:e])
    return "".join(parts)


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def count_features(fastqs: "List[str]", technology: str, features: "List[Tuple[str, str]]", whitelist: "Optional[str]" = None,
                   mismatch: bool = True, both_strands: bool = False):
    """Returns (AnnData barcodes x features of distinct-UMI counts, stats dict)."""
    np, pd, ad = _np(), _pd(), _ad()
    import scipy.sparse as sp  # type: ignore
    tech = parse_technology(technology)
    nfiles = 1 + max(f for g in tech.values() for f, _, _ in g)
    if len(fastqs) % nfiles:
        raise ValueError(f"{len(fastqs)} FASTQs given but the technology string uses {nfiles} files per set")
    exact, one = feature_lookup(features, mismatch)
    stats = {"reads": 0, "feature_exact": 0, "feature_1mm": 0, "feature_revcomp": 0, "no_feature": 0}
    records = []  # (barcode, umi, feature_idx)
    for k in range(0, len(fastqs), nfiles):
        its = [_fastq_seqs(p) for p in fastqs[k:k + nfiles]]
        for reads in zip(*its):
            stats["reads"] += 1
            reads = list(reads)
            bc, umi, feat = _cut(reads, tech["barcode"]), _cut(reads, tech["umi"]), _cut(reads, tech["feature"])
            fi = exact.get(feat)
            if fi is not None:
                stats["feature_exact"] += 1
            else:
                fi = one.get(feat)
                if fi is not None:
                    stats["feature_1mm"] += 1
                elif both_strands:
                    rc = _revcomp(feat)
                    fi = exact.get(rc, one.get(rc))
                    if fi is not None:
                        stats["feature_revcomp"] += 1
            if fi is None:
                stats["no_feature"] += 1
                continue
            records.append((bc, umi, fi))
    stats["barcodes_observed"] = len({r[0] for r in records})
    corr: "Dict[str, Optional[str]]" = {}
    if whitelist:
        observed = {r[0] for r in records}
        neigh: "Dict[str, List[str]]" = {}
        for b in observed:
            for v in hamming1(b):
                neigh.setdefault(v, []).append(b)
        in_wl_exact, near = set(), {}
        n_wl = 0
        try:
            with _open_text(whitelist) as fh:
                for line in fh:
                    w = line.strip().upper()
                    if not w:
                        continue
                    n_wl += 1
                    if w in observed:
                        in_wl_exact.add(w)
                    for b in neigh.get(w, ()):
                        near.setdefault(b, set()).add(w)
        except EOFError:  # truncated .gz (the upstream Workflows/*/3M-february-2018.txt.gz copies are cut short)
            stats["whitelist_truncated"] = True
            log.warning("whitelist %s is a truncated gzip; used the %d barcodes before the cut", whitelist, n_wl)
            print(f"warning: whitelist {whitelist} is truncated; used the first {n_wl} barcodes")
        stats["whitelist_barcodes_read"] = n_wl
        for b in observed:
            if b in in_wl_exact:
                corr[b] = b
            elif len(near.get(b, ())) == 1:
                corr[b] = next(iter(near[b]))
            else:
                corr[b] = None
        stats["barcodes_whitelist_exact"] = len(in_wl_exact)
        stats["barcodes_corrected_1mm"] = sum(1 for b, c in corr.items() if c is not None and c != b)
        stats["barcodes_rejected"] = sum(1 for c in corr.values() if c is None)
    triples = set()
    reads_kept = 0
    for bc, umi, fi in records:
        b = corr.get(bc, bc) if whitelist else bc
        if b is None:
            continue
        reads_kept += 1
        triples.add((b, umi, fi))
    stats["reads_assigned"] = reads_kept
    stats["umis"] = len(triples)
    bcs = sorted({t[0] for t in triples})
    bidx = {b: i for i, b in enumerate(bcs)}
    cnt: "Dict[Tuple[int, int], int]" = {}
    for b, _, fi in triples:
        key = (bidx[b], fi)
        cnt[key] = cnt.get(key, 0) + 1
    rows = np.array([k[0] for k in cnt], dtype=np.int64)
    cols = np.array([k[1] for k in cnt], dtype=np.int64)
    vals = np.array(list(cnt.values()), dtype=np.float32)
    X = sp.csr_matrix((vals, (rows, cols)), shape=(len(bcs), len(features)), dtype=np.float32)
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(bcs, name="barcode")),
                       var=pd.DataFrame({"sequence": [s for s, _ in features]}, index=pd.Index([i for _, i in features], name="gene_id")))
    stats["barcodes"] = len(bcs)
    return adata, stats


# ---------------------------------------------------------------------------
# AnnData QC and MuData
# ---------------------------------------------------------------------------

def knee_frame(X, obs_names):
    np, pd = _np(), _pd()
    tot = np.asarray(X.sum(1)).ravel()
    df = pd.DataFrame({"sum": tot, "barcodes": list(obs_names)}).sort_values("sum", ascending=False, kind="mergesort").reset_index(drop=True)
    df["sum_log"] = np.log1p(df["sum"])
    return df


def plot_knee(knee, path: Path, title: str = "Knee Plot") -> "Optional[Path]":
    plt = _plt()
    if plt is None:
        return None
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(knee.index, knee["sum_log"], marker="o", linestyle="-", markersize=3, color=BLUE, linewidth=1)
    _style(ax, title, "Barcode Index", "Log of UMI Counts")
    return _save(fig, path)


def preprocess_adata(adata, gene_names: "Optional[List[str]]", min_genes: int = 100, min_cells: int = 3, reference: str = "human",
                     mt_prefix: "Optional[str]" = None):
    """PreprocessAnnData; returns (filtered adata, knee frame, info)."""
    try:
        import scanpy as sc  # type: ignore
    except ImportError as e:  # pragma: no cover
        raise SystemExit("preprocess needs scanpy (pip install scanpy)") from e
    np = _np()
    if gene_names is not None:
        if len(gene_names) != adata.shape[1]:
            raise ValueError("The number of gene names does not match the number of variables in adata_rna")
        adata.var_names = gene_names
    knee = knee_frame(adata.X, adata.obs_names)
    n0 = adata.shape
    adata.obs["batch_number"] = 1
    sc.pp.filter_cells(adata, min_genes=min_genes)
    n1 = adata.shape
    sc.pp.filter_genes(adata, min_cells=min_cells)
    prefix = mt_prefix if mt_prefix is not None else ("MT-" if reference == "human" else "Mt-")
    adata.var["mt"] = adata.var_names.str.startswith(prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # scanpy's default percent_top (50, 100, 200, 500) raises on < 500 genes; keep the tiers that fit (identical otherwise)
    top = [t for t in (50, 100, 200, 500) if t <= adata.shape[1]] or None
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=top)
    info = {"cells_in": int(n0[0]), "genes_in": int(n0[1]), "cells_after_min_genes": int(n1[0]), "cells_out": int(adata.shape[0]),
            "genes_out": int(adata.shape[1]), "min_genes": min_genes, "min_cells": min_cells, "reference": reference, "mt_prefix": prefix,
            "n_mt_genes": int(adata.var["mt"].sum()), "n_ribo_genes": int(adata.var["ribo"].sum()),
            "median_total_counts": float(np.median(adata.obs["total_counts"])) if adata.shape[0] else float("nan"),
            "median_n_genes_by_counts": float(np.median(adata.obs["n_genes_by_counts"])) if adata.shape[0] else float("nan"),
            "median_pct_counts_mt": float(np.median(adata.obs["pct_counts_mt"])) if adata.shape[0] else float("nan")}
    return adata, knee, info


def plot_qc(adata, out: Path) -> "List[Path]":
    plt = _plt()
    if plt is None or adata.shape[0] == 0:
        return []
    np = _np()
    figs = []
    keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6))
    rng = np.random.default_rng(0)
    for ax, k in zip(axes, keys):
        v = adata.obs[k].to_numpy(dtype=float)
        parts = ax.violinplot([v], showmedians=True)
        for b in parts["bodies"]:
            b.set_facecolor(BLUE); b.set_alpha(0.35)
        ax.scatter(1 + rng.uniform(-0.2, 0.2, v.size), v, s=2, color=INK2, alpha=0.5)  # jitter 0.4
        ax.set_xticks([])
        _style(ax, k)
    figs.append(_save(fig, out / "violin_plot_scrna.png"))
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sca = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"], c=adata.obs["pct_counts_mt"], s=5, cmap="viridis")
    fig.colorbar(sca, ax=ax, label="pct_counts_mt")
    _style(ax, "QC scatter", "total_counts", "n_genes_by_counts")
    figs.append(_save(fig, out / "scatter_plot_scrna.png"))
    return figs


def guide_var_names(adata_guide, meta, upstream_compat: bool = False) -> "Tuple[List[str], str]":
    names = (meta["sgRNA_ID"].astype(str) + "|" + meta["sgRNA_sequences"].astype(str)).tolist()
    ids = meta["sgRNA_ID"].astype(str).tolist()
    current = [str(v) for v in adata_guide.var_names]
    if not upstream_compat and set(current) <= set(ids) and len(set(ids)) == len(ids):
        m = dict(zip(ids, names))
        return [m[c] for c in current], "by_id"
    if len(names) != adata_guide.shape[1]:
        raise ValueError(f"guide metadata has {len(names)} rows but the guide matrix has {adata_guide.shape[1]} features")
    return names, "positional"


def create_mdata(adata_rna, adata_guide, meta, upstream_compat: bool = False):
    """CreateMuData core; returns (rna_sub, guide_sub, knee, info)."""
    np = _np()
    names, how = guide_var_names(adata_guide, meta, upstream_compat)
    adata_guide.var_names = names
    X = adata_guide.X
    adata_guide.obs["number_of_nonzero_guides"] = np.asarray((X > 0).sum(axis=1)).ravel()
    adata_guide.obs["batch_number"] = 1
    knee = knee_frame(X, adata_guide.obs_names)
    gset = set(adata_guide.obs_names)
    if upstream_compat:
        inter = list(set(adata_rna.obs_names).intersection(adata_guide.obs_names))
    else:
        inter = [b for b in adata_rna.obs_names if b in gset]
    rna_sub = adata_rna[inter, :].copy()
    guide_sub = adata_guide[inter, :].copy()
    info = {"rna_barcodes": int(adata_rna.shape[0]), "guide_barcodes": int(adata_guide.shape[0]), "intersecting_barcodes": len(inter),
            "guides": int(adata_guide.shape[1]), "genes": int(adata_rna.shape[1]), "guide_naming": how,
            "median_nonzero_guides": float(np.median(guide_sub.obs["number_of_nonzero_guides"])) if len(inter) else float("nan"),
            "cells_with_guide": int((guide_sub.obs["number_of_nonzero_guides"] > 0).sum()) if len(inter) else 0}
    return rna_sub, guide_sub, knee, info


def _elem_io():
    try:
        from anndata.io import read_elem, write_elem  # type: ignore  # anndata >= 0.11
    except ImportError:
        from anndata.experimental import read_elem, write_elem  # type: ignore
    return read_elem, write_elem


def write_h5mu(path: Path, mods: "Dict[str, Any]") -> str:
    """MuData file with the given modalities; mudata when installed, else an h5mu written with h5py + anndata."""
    try:
        import mudata  # type: ignore
        mudata.MuData(mods).write(str(path))
        return "mudata"
    except ImportError:
        pass
    np, pd = _np(), _pd()
    import h5py  # type: ignore
    _, write_elem = _elem_io()
    names = list(mods)
    obs_index = list(mods[names[0]].obs_names)
    var_index: "List[str]" = []
    for n in names:
        var_index += list(mods[n].var_names)
    with h5py.File(str(path), "w") as f:
        f.attrs["encoding-type"] = "MuData"
        f.attrs["encoding-version"] = "0.1.0"
        f.attrs["encoder"] = "igvfagent"
        f.attrs["encoder_version"] = "1"
        write_elem(f, "obs", pd.DataFrame(index=pd.Index(obs_index)))
        write_elem(f, "var", pd.DataFrame(index=pd.Index(var_index)))
        obsm, varm, off = {}, {}, 0
        omap, vmap = f.create_group("obsmap"), f.create_group("varmap")
        for n in names:
            a = mods[n]
            obsm[n] = np.isin(np.array(obs_index, dtype=object), np.array(list(a.obs_names), dtype=object))
            vm = np.zeros(len(var_index), dtype=bool); vm[off:off + a.shape[1]] = True
            varm[n] = vm
            pos = {b: i + 1 for i, b in enumerate(a.obs_names)}
            omap.create_dataset(n, data=np.array([pos.get(b, 0) for b in obs_index], dtype=np.uint32))
            vv = np.zeros(len(var_index), dtype=np.uint32); vv[off:off + a.shape[1]] = np.arange(1, a.shape[1] + 1)
            vmap.create_dataset(n, data=vv)
            off += a.shape[1]
        write_elem(f, "obsm", obsm)
        write_elem(f, "varm", varm)
        write_elem(f, "obsp", {})
        write_elem(f, "varp", {})
        write_elem(f, "uns", {})
        g = f.create_group("mod")
        for n in names:
            write_elem(g, n, mods[n])
        g.attrs["mod-order"] = np.array(names, dtype=object).astype("S")
    return "h5py"


def read_h5mu_mod(path: Path, name: str):
    import h5py  # type: ignore
    read_elem, _ = _elem_io()
    with h5py.File(str(path), "r") as f:
        return read_elem(f["mod"][name])


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------

def cmd_parse_seqspec(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    directory = args.directory or str(Path(args.yaml).resolve().parent)
    df, info = parse_seqspec_file(args.yaml, args.modality, directory, args.engine, args.strict_region_types)
    name = f"{'_'.join(args.modality)}_parsed_seqSpec.txt"
    df.to_csv(out / name, sep="\t", index=False)
    print(f"TSV: {out / name}")
    rows = []
    for m, d in info["modalities"].items():
        wl = os.path.join(directory, d["onlist"]) if d["onlist"] else ""
        rows.append([m, ",".join(d["reads"]), f"`{d['representation']}`", d["engine"], wl,
                     "yes" if wl and os.path.exists(wl) else ("gz only" if wl and os.path.exists(wl + ".gz") else "no")])
        print(f"{m}: reads={','.join(d['reads'])} x={d['representation']}")
    notes = []
    for m, d in info["modalities"].items():
        if d["cell_barcode_as_barcode"]:
            notes.append(f"- {m}: region(s) {d['cell_barcode_as_barcode']} have region_type cell_barcode; treated as barcode "
                         "(the seqspec binary's kb formatter only recognises `barcode`; --strict-region-types reproduces it).")
        if d["seqspec_binary"] is not None and d["seqspec_binary"] != d["python_representation"]:
            notes.append(f"- {m}: seqspec binary `{d['seqspec_binary']}` differs from the Python parser `{d['python_representation']}`.")
    coord_rows = [[m, r, c["region_id"], c["region_type"], c["start"], c["stop"]]
                  for m, d in info["modalities"].items() for r, cs in d["coordinates"].items() for c in cs]
    summary = {"yaml": str(args.yaml), "assay_id": info["assay_id"], "directory": directory, "parsed_seqspec": str(out / name),
               "modalities": {m: {k: v for k, v in d.items() if k != "coordinates"} for m, d in info["modalities"].items()}}
    write_json(summary, out / "summary.json")
    write_report(out / "report.md", f"seqspec parse: {Path(args.yaml).name}", [
        f"Assay `{info['assay_id']}`.  Output `{name}` (columns {', '.join(PARSED_COLS)}).", "",
        md_table(["modality", "reads", "kb -x", "engine", "barcode_whitelist", "whitelist present"], rows), "",
        "## Read-to-region coordinates", "", md_table(["modality", "read", "region", "type", "start", "stop"], coord_rows), "",
        *(["## Notes", ""] + notes if notes else []),
        "The Nextflow step expects the onlist uncompressed next to the seqspec (`<dir>/3M-february-2018.txt`)."])
    return 0


def cmd_extract_seqspec(args: argparse.Namespace) -> int:
    x = extract_parsed_seqspec(args.file)
    if x is None:
        print("Cannot find the parsed_seqspec file")
        return 1
    print(x)
    return 0


def cmd_guide_features(args: argparse.Namespace) -> int:
    out = Path(args.out) if args.out else run_dir(args.label) / "guide_features.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    t = read_guide_table(args.guide_table)
    write_guide_features(t, out)
    if not args.out:
        d = out.parent
        dup_seq = int(t["sgRNA_sequences"].duplicated().sum()); dup_id = int(t["sgRNA_ID"].duplicated().sum())
        lens = t["sgRNA_sequences"].astype(str).str.len().value_counts().to_dict()
        write_json({"guide_table": str(args.guide_table), "guides": int(len(t)), "duplicate_sequences": dup_seq, "duplicate_ids": dup_id,
                    "length_counts": {int(k): int(v) for k, v in lens.items()}, "guide_features": str(out)}, d / "summary.json")
        write_report(d / "report.md", "Guide features", [
            f"{len(t)} guides from `{args.guide_table}` -> `guide_features.txt` (sgRNA_sequences, sgRNA_ID; tab, no header).", "",
            f"Duplicate sequences: {dup_seq}; duplicate IDs: {dup_id}; lengths: {lens}.", "",
            md_table(list(t.columns[:6]), t.iloc[:10, :6].values.tolist())])
    return 0


def cmd_process_reads(args: argparse.Namespace) -> int:
    print(" ".join(pair_reads(args.dir, args.reads1, args.reads2)))
    return 0


def cmd_download_genome(args: argparse.Namespace) -> int:
    out = Path(args.out)
    cmd = ["wget", args.url, "-O", str(out)]
    print("note: hg38.fa.gz is ~1 GB; with kite (kb ref --workflow kite) the upstream passes it as -f1, which kb overwrites, "
          "so it is only needed for kb-ref --upstream-compat")
    if not args.execute:
        print(f"Command: {' '.join(cmd)}")
        print("not executed (pass --execute to download)")
        return 0
    if shutil.which("wget"):
        status, rc = run_or_print(cmd, "wget", False)
        return rc
    import urllib.request
    urllib.request.urlretrieve(args.url, str(out))
    print(f"Wrote: {out}")
    return 0


def cmd_kb_ref(args: argparse.Namespace) -> int:
    out = Path(args.out_dir) if args.out_dir else run_dir(args.label)
    out.mkdir(parents=True, exist_ok=True)
    gf = args.guide_features
    if args.mode == "kite":
        if args.guide_table and not gf:
            gf = str(write_guide_features(read_guide_table(args.guide_table), out / "guide_features.txt"))
        if not gf:
            raise SystemExit("kb-ref --mode kite needs --guide-features or --guide-table")
    fasta = args.genome if (args.upstream_compat and args.genome) else None
    cmd = kb_ref_command(args.mode, args.species, gf or "guide_features.txt", fasta, out)
    status, rc = run_or_print(cmd, "kb", args.dry_run)
    write_json({"mode": args.mode, "command": cmd, "status": status, "returncode": rc}, out / "summary.json")
    write_report(out / "report.md", f"kb ref ({args.mode})", ["```bash", " ".join(cmd), "```", "", f"Status: {status}."])
    return rc


def cmd_kb_count(args: argparse.Namespace) -> int:
    out = Path(args.out_dir) if args.out_dir else run_dir(args.label)
    out.mkdir(parents=True, exist_ok=True)
    tech = args.technology or (extract_parsed_seqspec(args.parsed_seqspec) if args.parsed_seqspec else None)
    if not tech:
        raise SystemExit("kb-count needs --technology or a readable --parsed-seqspec")
    wl = args.whitelist
    if not wl and args.parsed_seqspec:
        wl = str(_pd().read_csv(args.parsed_seqspec, sep="\t")["barcode_whitelist"].iloc[0]).strip()
    fastqs = list(args.fastqs or [])
    if args.fastq_dir:
        fastqs = pair_reads(args.fastq_dir, args.reads1 or "", args.reads2 or "")
    if not fastqs:
        raise SystemExit("kb-count needs --fastqs or --fastq-dir with --reads1/--reads2")
    cmd = kb_count_command(args.index, args.t2g, wl or "whitelist.txt", tech, str(out / KB_OUT[args.modality]), fastqs, args.threads, args.memory)
    status, rc = run_or_print(cmd, "kb", args.dry_run)
    h5 = out / KB_OUT[args.modality] / "counts_unfiltered" / "adata.h5ad"
    if h5.is_file():
        print(f"Wrote: {h5}")
    write_json({"modality": args.modality, "technology": tech, "whitelist": wl, "command": cmd, "status": status, "returncode": rc,
                "adata": str(h5) if h5.is_file() else None}, out / "summary.json")
    write_report(out / "report.md", f"kb count ({args.modality})", ["```bash", " ".join(cmd), "```", "", f"Status: {status}."])
    return rc


def _features_from_args(args, out: Path) -> "Tuple[List[Tuple[str, str]], Optional[Any]]":
    if getattr(args, "guide_features", None):
        return read_guide_features(args.guide_features), None
    t = read_guide_table(args.guide_table)
    write_guide_features(t, out / "guide_features.txt")
    return [(str(s).upper(), str(i)) for s, i in zip(t["sgRNA_sequences"], t["sgRNA_ID"])], t


def _tech_from_args(args) -> "Tuple[str, Optional[str], dict]":
    if args.technology:
        return args.technology, args.whitelist, {}
    if args.parsed_seqspec:
        wl = args.whitelist or str(_pd().read_csv(args.parsed_seqspec, sep="\t")["barcode_whitelist"].iloc[0]).strip()
        return extract_parsed_seqspec(args.parsed_seqspec) or "", wl, {}
    if args.seqspec:
        spec = load_seqspec(args.seqspec)
        x, det = kb_technology(spec, args.modality)
        wl = args.whitelist
        if not wl:
            f = extract_whitelist_file(args.seqspec)
            if f:
                for cand in (Path(args.seqspec).parent / f, Path(args.seqspec).parent / (f + ".gz")):
                    if cand.is_file():
                        wl = str(cand)
                        break
        return x, wl, det
    raise SystemExit("count-features needs --technology, --parsed-seqspec or --seqspec")


def _run_count_features(args, out: Path, fastqs: "List[str]") -> "Tuple[Path, dict]":
    feats, _ = _features_from_args(args, out)
    tech, wl, _ = _tech_from_args(args)
    if args.no_whitelist:
        wl = None
    adata, stats = count_features(fastqs, tech, feats, wl, mismatch=not args.exact, both_strands=args.both_strands)
    cu = out / "counts_unfiltered"
    cu.mkdir(parents=True, exist_ok=True)
    h5 = cu / "adata.h5ad"
    adata.write_h5ad(h5)
    print(f"Wrote: {h5}")
    stats.update({"technology": tech, "whitelist": wl, "features": len(feats), "adata": str(h5),
                  "features_detected": int((adata.X.sum(0) > 0).sum()) if adata.shape[0] else 0})
    return h5, stats


def cmd_count_features(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    h5, stats = _run_count_features(args, out, args.fastqs)
    ad = _ad()
    a = ad.read_h5ad(h5)
    tot = _np().asarray(a.X.sum(0)).ravel()
    top = sorted(zip(a.var_names, tot), key=lambda t: -t[1])[:15]
    pd = _pd()
    tab = pd.DataFrame({"feature": a.var_names, "umis": tot, "barcodes": _np().asarray((a.X > 0).sum(0)).ravel()})
    tab.to_csv(out / "feature_totals.tsv", sep="\t", index=False)
    print(f"TSV: {out / 'feature_totals.tsv'}")
    if not args.no_plots and a.shape[0]:
        plot_knee(knee_frame(a.X, a.obs_names), out / "knee_plot_guide.png")
    write_json(stats, out / "summary.json")
    write_report(out / "report.md", "Guide counts (kite-style, pure Python)", [
        f"Technology `{stats['technology']}`; whitelist `{stats['whitelist']}`.", "",
        md_table(["metric", "value"], [[k, v] for k, v in stats.items() if isinstance(v, (int, float))]), "",
        "## Top features by UMIs", "", md_table(["feature", "UMIs"], [[f, int(v)] for f, v in top]), "",
        "Reads are split by the technology string; the feature is matched exactly or by a unique Hamming-1 neighbour (kite's mismatch "
        "index); barcodes are kept if on the whitelist or one unique mismatch away; counts are distinct (barcode, UMI, feature) triples."])
    return 0


def _load_adata(path: str):
    ad = _ad()
    return ad.read_h5ad(path)


def _read_gene_names(path: "Optional[str]") -> "Optional[List[str]]":
    if not path:
        return None
    pd = _pd()
    return pd.read_csv(path, header=None)[0].astype(str).tolist()


def _do_preprocess(args, out: Path) -> "Tuple[Path, dict]":
    a = _load_adata(args.adata_rna)
    a, knee, info = preprocess_adata(a, _read_gene_names(args.gene_names), args.min_genes, args.min_cells, args.reference, args.mt_prefix)
    figs = []
    if not args.no_plots:
        f = plot_knee(knee, out / "knee_plot_rna.png")
        figs += [f] if f else []
        figs += plot_qc(a, out)
    h5 = out / "filtered_anndata.h5ad"
    a.write_h5ad(h5)
    print(f"Wrote: {h5}")
    qc = a.obs.reset_index().rename(columns={"index": "barcode"})
    qc.to_csv(out / "cell_qc.tsv", sep="\t", index=False)
    print(f"TSV: {out / 'cell_qc.tsv'}")
    info.update({"filtered_anndata": str(h5), "figures": [str(p) for p in figs]})
    return h5, info


def cmd_preprocess(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    _, info = _do_preprocess(args, out)
    write_json(info, out / "summary.json")
    write_report(out / "report.md", "RNA AnnData QC", [
        md_table(["metric", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in info.items() if not isinstance(v, list)]), "",
        f"filter_cells(min_genes={info['min_genes']}) then filter_genes(min_cells={info['min_cells']}); mt prefix `{info['mt_prefix']}`, "
        "ribo RPS/RPL; calculate_qc_metrics(qc_vars=['mt','ribo'], log1p=True)."])
    return 0


def _do_create_mdata(adata_rna_path: str, adata_guide_path: str, meta_path: str, out: Path, upstream_compat: bool, no_plots: bool) -> dict:
    rna = _load_adata(adata_rna_path)
    guide = _load_adata(adata_guide_path)
    meta = read_guide_table(meta_path)
    rsub, gsub, knee, info = create_mdata(rna, guide, meta, upstream_compat)
    figs = []
    if not no_plots and len(knee):
        f = plot_knee(knee, out / "knee_plot_guide.png")
        figs += [f] if f else []
    h5mu = out / "mudata.h5mu"
    info["writer"] = write_h5mu(h5mu, {"transcripts": rsub, "guides": gsub})
    print(f"Wrote: {h5mu}")
    gsub.obs[["number_of_nonzero_guides", "batch_number"]].reset_index().rename(columns={"index": "barcode"}).to_csv(
        out / "guide_cell_summary.tsv", sep="\t", index=False)
    print(f"TSV: {out / 'guide_cell_summary.tsv'}")
    info.update({"mudata": str(h5mu), "figures": [str(p) for p in figs]})
    return info


def cmd_create_mdata(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    info = _do_create_mdata(args.adata_rna, args.adata_guide, args.guide_metadata, out, args.upstream_compat, args.no_plots)
    write_json(info, out / "summary.json")
    write_report(out / "report.md", "MuData assembly", [
        md_table(["metric", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in info.items() if not isinstance(v, list)]), "",
        "Modalities `transcripts` and `guides` over the intersecting barcodes; guide features named `sgRNA_ID|sgRNA_sequences`."])
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    out = run_dir(args.label)
    steps, summary = [], {"steps": {}}
    directory = None
    # 1. seqspec
    for mod, y in (("rna", args.rna_seqspec), ("guide", args.guide_seqspec)):
        if not y:
            continue
        directory = str(Path(y).resolve().parent)
        df, info = parse_seqspec_file(y, [mod], directory, args.engine)
        p = out / f"{mod}_parsed_seqSpec.txt"
        df.to_csv(p, sep="\t", index=False)
        print(f"TSV: {p}")
        summary["steps"][f"seqspec_{mod}"] = {"representation": info["modalities"][mod]["representation"], "reads": info["modalities"][mod]["reads"]}
        steps.append([f"seqspec ({mod})", f"`{info['modalities'][mod]['representation']}`"])
    # 2. guide features
    if args.guide_table:
        t = read_guide_table(args.guide_table)
        write_guide_features(t, out / "guide_features.txt")
        steps.append(["guide features", f"{len(t)} guides"])
    # 3. guide counts
    adata_guide = args.adata_guide
    if not adata_guide and args.guide_fastqs:
        if not args.guide_table:
            raise SystemExit("run: counting guides needs --guide-table")
        gdir = out / KB_OUT["guide"]
        gdir.mkdir(exist_ok=True)
        if shutil.which("kb") and shutil.which("kallisto") and not args.python_counts:
            ref = kb_ref_command("kite", guide_features=str(out / "guide_features.txt"), out_dir=out)
            s1, rc1 = run_or_print(ref, "kb", False)
            tech = extract_parsed_seqspec(out / "guide_parsed_seqSpec.txt") if args.guide_seqspec else args.technology
            wl = args.whitelist or (str(_pd().read_csv(out / "guide_parsed_seqSpec.txt", sep="\t")["barcode_whitelist"].iloc[0]) if args.guide_seqspec else "")
            cnt = kb_count_command(str(out / "guide_index.idx"), str(out / "t2guide.txt"), wl, tech or "", str(gdir), args.guide_fastqs)
            s2, rc2 = run_or_print(cnt, "kb", False)
            adata_guide = str(gdir / "counts_unfiltered" / "adata.h5ad")
            steps.append(["guide counts", f"kb ref/count ({s1}/{s2})"])
        else:
            ns = argparse.Namespace(guide_features=None, guide_table=args.guide_table, technology=args.technology, whitelist=args.whitelist,
                                    parsed_seqspec=None, seqspec=args.guide_seqspec, modality="guide", no_whitelist=args.no_whitelist,
                                    exact=False, both_strands=False)
            h5, st = _run_count_features(ns, gdir, args.guide_fastqs)
            adata_guide = str(h5)
            summary["steps"]["count_features"] = st
            steps.append(["guide counts", f"count-features: {st['reads']} reads, {st['umis']} UMIs, {st['barcodes']} barcodes"])
    # 4. preprocess
    filtered = None
    if args.adata_rna:
        ns = argparse.Namespace(adata_rna=args.adata_rna, gene_names=args.gene_names, min_genes=args.min_genes, min_cells=args.min_cells,
                                reference=args.reference, mt_prefix=args.mt_prefix, no_plots=args.no_plots)
        h5, info = _do_preprocess(ns, out)
        filtered = str(h5)
        summary["steps"]["preprocess"] = info
        steps.append(["preprocess", f"{info['cells_in']} -> {info['cells_out']} cells, {info['genes_in']} -> {info['genes_out']} genes"])
    # 5. MuData
    if filtered and adata_guide and args.guide_table and Path(adata_guide).is_file():
        info = _do_create_mdata(filtered, adata_guide, args.guide_table, out, args.upstream_compat, args.no_plots)
        summary["steps"]["create_mdata"] = info
        steps.append(["create MuData", f"{info['intersecting_barcodes']} barcodes x ({info['genes']} genes, {info['guides']} guides)"])
    summary["run_dir"] = str(out)
    write_json(summary, out / "summary.json")
    write_report(out / "report.md", "IGVF CRISPR pipeline run", [md_table(["step", "result"], steps), "",
                 "Steps not given inputs were skipped.  Order follows the upstream workflows: seqSpecParser -> creatingGuideRef -> "
                 "mappingGuide / mappingscRNA -> PreprocessAnnData -> CreateMuData."])
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

# Upstream example seqspecs (Modules/seqSpecParser/example_data), compacted; the rna expected output is the upstream's
# own Modules/mappingscRNA/example_data/rna_parsed_seqSpec.txt.
_REGION = "  - !Region\n    parent_id: {p}\n    region_id: {i}\n    region_type: {t}\n    name: {i}\n    sequence_type: {st}\n    sequence: {s}\n    min_len: {mn}\n    max_len: {mx}\n    onlist: {ol}\n    regions: null\n"
_ONLIST = "!Onlist\n      location: remote\n      filename: {f}\n      md5: updateit!"


def _example_yaml(mod: str, r1: int, r2: int, umi: int, cdna: "Tuple[int, int]", common: "Optional[str]" = None, bc_type: str = "barcode",
                  onlist: str = "3M-february-2018.txt") -> str:
    head = (f"!Assay\nseqspec_version: 0.2.0\nassay_id: test_{mod}\nname: test_{mod}\ndoi: 'XXX'\ndate: 08 July 2022\n"
            f"description: example {mod}  # comment\nmodalities:\n- {mod}\nlib_struct: https://teichlab.github.io/x.html\n"
            f"sequence_spec:\n- !Read\n  read_id: {mod}_R1.fq\n  name: Read 1\n  modality: {mod}\n  primer_id: r1_primer\n  min_len: {r1}\n"
            f"  max_len: {r1}\n  strand: pos\n- !Read\n  read_id: {mod}_R2.fq\n  name: Read 2\n  modality: {mod}\n  primer_id: r2_primer\n"
            f"  min_len: {r2}\n  max_len: {r2}\n  strand: neg\nlibrary_spec:\n- !Region\n  parent_id: null\n  region_id: {mod}\n"
            f"  region_type: null\n  name: null\n  sequence_type: null\n  sequence: AAAANNNN \n  min_len: 150\n  max_len: 150\n  onlist: null\n  regions:\n")
    body = _REGION.format(p=mod, i="r1_primer", t="r1_primer", st="fixed", s="A" * 16, mn=16, mx=16, ol="null")
    body += _REGION.format(p=mod, i="barcode", t=bc_type, st="onlist", s="N" * 16, mn=16, mx=16, ol=_ONLIST.format(f=onlist))
    body += _REGION.format(p=mod, i="umi", t="umi", st="fixed", s="N" * umi, mn=umi, mx=umi, ol="null")
    body += _REGION.format(p=mod, i="cdna", t="cdna", st="random", s="N" * cdna[1], mn=cdna[0], mx=cdna[1], ol="null")
    if common:
        body += _REGION.format(p=mod, i="common", t="common", st="fixed", s=common, mn=len(common), mx=len(common), ol="null")
    body += _REGION.format(p=mod, i="r2_primer", t="r2_primer", st="fixed", s="A" * 16, mn=16, mx=16, ol="null")
    return head + body


COMMON = "CTTGTGGAAAGGACGAAACACCG"


def _synthetic_guides(d: Path, rng) -> dict:
    """10 guides, 60 cells; each cell carries one planted guide (40 cells) or two (10) or none (10); FASTQs in the guide.yml layout."""
    np, pd = _np(), _pd()
    bases = np.array(list("ACGT"))
    guides = ["".join(rng.choice(bases, 20)) for _ in range(10)]
    ids = [f"G{i}_sg{1 + i % 2}" for i in range(10)]
    table = pd.DataFrame({"sgRNA_ID": ids, "Target_name": [i.split("_")[0] for i in ids], "sgRNA_sequences": guides,
                          "chr": "chr1", "start": np.arange(10) * 1000, "end": np.arange(10) * 1000 + 20})
    cells = ["".join(rng.choice(bases, 16)) for _ in range(60)]
    wl = cells + ["".join(rng.choice(bases, 16)) for _ in range(500)]
    wl_path = d / "3M-february-2018.txt"
    wl_path.write_text("\n".join(wl) + "\n")
    truth: "Dict[str, Dict[int, int]]" = {}
    r1, r2 = [], []
    k = 0
    for ci, cb in enumerate(cells):
        assigned = [] if ci >= 50 else ([ci % 10] if ci < 40 else [ci % 10, (ci + 3) % 10])
        truth[cb] = {}
        for g in assigned:
            n_umi = 3 + (ci % 4)
            truth[cb][g] = n_umi
            for u in range(n_umi):
                umi = "".join(rng.choice(bases, 12))
                for rep in range(2):  # each UMI sequenced twice (collapsed)
                    feat = guides[g]
                    if rep == 1 and u == 0:  # one read per molecule with a single guide mismatch -> Hamming-1 rescue
                        pos = int(rng.integers(0, 20))
                        feat = feat[:pos] + ("A" if feat[pos] != "A" else "C") + feat[pos + 1:]
                    bc = cb
                    if rep == 1 and u == 1:  # a barcode sequencing error -> whitelist correction
                        bc = cb[:5] + ("T" if cb[5] != "T" else "G") + cb[6:]
                    r1.append(bc + umi + "".join(rng.choice(bases, 122)))
                    r2.append(COMMON + feat + "".join(rng.choice(bases, 150 - 43)))
                    k += 1
        # background reads with no guide (filtered)
        r1.append(cb + "".join(rng.choice(bases, 134)))
        r2.append("".join(rng.choice(bases, 150)))
    fq1, fq2 = d / "guide_R1.fq", d / "guide_R2.fq"
    for p, seqs in ((fq1, r1), (fq2, r2)):
        with open(p, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f"@r{i}\n{s}\n+\n{'F' * len(s)}\n")
    # gzip copy of read 1 to exercise .gz input
    with gzip.open(str(fq1) + ".gz", "wt") as fh:
        fh.write(fq1.read_text())
    xlsx = d / "guides.xlsx"
    table.to_excel(xlsx, index=False)
    tsv = d / "guides.tsv"
    table.to_csv(tsv, sep="\t", index=False)
    return {"table": table, "xlsx": xlsx, "tsv": tsv, "cells": cells, "truth": truth, "fq1": fq1, "fq2": fq2, "whitelist": wl_path,
            "guides": guides, "ids": ids}


def _synthetic_rna(d: Path, cells: "List[str]", rng):
    """RNA counts for the 60 guide cells + 40 RNA-only cells; 13 cells planted with < 100 genes, 30 genes seen in < 3 cells."""
    np, pd, ad = _np(), _pd(), _ad()
    import scipy.sparse as sp  # type: ignore
    extra = ["".join(rng.choice(np.array(list("ACGT")), 16)) for _ in range(40)]
    obs = cells + extra
    genes = [f"MT-ND{i}" for i in range(1, 6)] + [f"RPL{i}" for i in range(1, 6)] + [f"RPS{i}" for i in range(1, 4)] + \
            ["Mt-co1", "ATP5F1"] + [f"GENE{i:04d}" for i in range(600)] + [f"RARE{i:02d}" for i in range(30)]
    n, g = len(obs), len(genes)
    X = np.zeros((n, g), dtype=np.float32)
    low = set(range(0, 100, 9)) | {5}  # 13 low-complexity barcodes
    for i in range(n):
        ng = 40 if i in low else 250
        cols = rng.choice(np.arange(g - 30), ng, replace=False)
        X[i, cols] = rng.poisson(3, ng) + 1
        X[i, :5] = rng.poisson(10 if i % 2 == 0 else 1, 5)  # mito genes
    for j in range(g - 30, g):  # rare genes: at most 2 cells
        X[rng.choice(n, 2, replace=False), j] = 1
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=obs), var=pd.DataFrame(index=[f"ENSG{i:011d}" for i in range(g)]))
    h5 = d / "adata_rna.h5ad"
    a.write_h5ad(h5)
    gn = d / "cells_x_genes.genes.names.txt"
    gn.write_text("\n".join(genes) + "\n")
    return {"h5ad": h5, "genes": gn, "X": X, "obs": obs, "gene_list": genes, "low": low}


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    np, pd = _np(), _pd()
    checks: "List[Tuple[bool, str]]" = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made: "List[Path]" = []

    def newest(label: str) -> Path:
        p = sorted(OUT_ROOT.glob(f"*_{label}"))[-1]
        made.append(p)
        return p

    rng = np.random.default_rng(11)
    before = set(OUT_ROOT.glob("*_st_icp_*")) if OUT_ROOT.is_dir() else set()
    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            print("\nseqspec YAML parser")
            rna_y = d / "rna.yml"; rna_y.write_text(_example_yaml("rna", 28, 150, 12, (1, 150)))
            guide_y = d / "guide.yml"; guide_y.write_text(_example_yaml("guide", 26, 43, 12, (20, 20), COMMON))
            guide2_y = d / "guide2.yml"; guide2_y.write_text(_example_yaml("guide", 26, 43, 10, (20, 20), "CATCCACCAAACCCATACGTCCGAGT"))
            ms_y = d / "multiseq.yml"; ms_y.write_text(_example_yaml("multiseq", 26, 150, 12, (8, 8), bc_type="cell_barcode"))
            spec = load_seqspec(rna_y)
            check(spec["assay_id"] == "test_rna" and spec["modalities"] == ["rna"] and spec["doi"] == "XXX" and spec["description"] == "example rna",
                  "parser: scalars, quoted strings, same-indent lists, trailing comments")
            leaves = _leaves(spec["library_spec"][0])
            check([l["region_id"] for l in leaves] == ["r1_primer", "barcode", "umi", "cdna", "r2_primer"]
                  and leaves[1]["onlist"]["filename"] == "3M-february-2018.txt" and leaves[0]["onlist"] is None,
                  "parser: nested !Region lists and the !Onlist mapping")
            check(seqspec_reads(spec, "rna") == ["rna_R1.fq", "rna_R2.fq"] and seqspec_reads(spec, "guide") == [], "process_reads: read_ids of the modality")
            x_rna, _ = kb_technology(spec, "rna")
            check(x_rna == "0,0,16:0,16,28:1,0,150", f"kb -x rna = {x_rna} (upstream rna_parsed_seqSpec.txt: 0,0,16:0,16,28:1,0,150)")
            x_g, _ = kb_technology(load_seqspec(guide_y), "guide")
            check(x_g == "0,0,16:0,16,26:1,23,43", f"kb -x guide = {x_g} (UMI clipped at read-1 length 26; protospacer after the 23-nt common on read 2)")
            x_g2, _ = kb_technology(load_seqspec(guide2_y), "guide")
            check(x_g2 == "0,0,16:0,16,26:1,26,43", f"kb -x guide2 (10x v2, 26-nt common) = {x_g2}")
            x_ms, det = kb_technology(load_seqspec(ms_y), "multiseq")
            x_ms_strict, _ = kb_technology(load_seqspec(ms_y), "multiseq", strict=True)
            # the upstream multiseq.yml declares a 150-nt read 2 over a 52-nt insert, so read 2 also spans UMI and barcode (seqspec does the same)
            check(x_ms == "0,0,16,1,20,36:0,16,26,1,8,20:1,0,8" and x_ms_strict == "-1,-1,-1:0,16,26,1,8,20:1,0,8"
                  and det["cell_barcode_as_barcode"] == ["barcode", "barcode"],
                  f"kb -x multiseq = {x_ms} (cell_barcode as barcode; read 2 of 150 nt runs through the 52-nt insert); strict {x_ms_strict}")
            check(extract_whitelist_file(rna_y) == "3M-february-2018.txt", "onlist: first 'filename' line")
            check(parse_technology(x_g) == {"barcode": [(0, 0, 16)], "umi": [(0, 16, 26)], "feature": [(1, 23, 43)]}
                  and parse_technology("-1,-1,-1:0,0,12:1,0,0")["barcode"] == [], "parse_technology: triples, -1 groups dropped")

            print("\nparse-seqspec / extract-seqspec")
            rc = cmd_parse_seqspec(argparse.Namespace(yaml=str(guide_y), modality=["guide"], directory="example_data", engine="python",
                                                      strict_region_types=False, label="st_icp_seqspec"))
            run = newest("st_icp_seqspec")
            ps = pd.read_csv(run / "guide_parsed_seqSpec.txt", sep="\t")
            check(rc == 0 and list(ps.columns) == PARSED_COLS and ps.iloc[0]["barcode_whitelist"] == "example_data/3M-february-2018.txt"
                  and ps.iloc[0]["representation"] == x_g and (run / "report.md").is_file(), "parse-seqspec: 4-column TSV, whitelist joined to directory")
            check(extract_parsed_seqspec(run / "guide_parsed_seqSpec.txt") == x_g and extract_parsed_seqspec(d / "missing.txt") is None,
                  "extract-seqspec: representation column; missing file -> None")
            two = pd.concat([ps, ps.assign(representation="0,0,1:0,1,2:1,0,3")])
            two.to_csv(d / "two.txt", sep="\t", index=False)
            check(extract_parsed_seqspec(d / "two.txt") == x_g + "0,0,1:0,1,2:1,0,3", "extract-seqspec: rows concatenated without a separator (upstream behaviour)")

            print("\nguide-features / process-reads")
            S = _synthetic_guides(d, rng)
            out_gf = d / "gf.txt"
            rc = cmd_guide_features(argparse.Namespace(guide_table=str(S["xlsx"]), out=str(out_gf), label="x"))
            gf = pd.read_csv(out_gf, sep="\t", header=None)
            check(rc == 0 and gf.shape == (10, 2) and gf[0].tolist() == S["guides"] and gf[1].tolist() == S["ids"],
                  "guide-features: xlsx -> sgRNA_sequences<TAB>sgRNA_ID, no header")
            check(read_guide_table(S["tsv"])["sgRNA_ID"].tolist() == S["ids"], "guide-features: TSV guide tables accepted")
            pr = pair_reads("fq", "a_R1.fq.gz; b_R1.fq.gz ;c_R1.fq.gz", " a_R2.fq.gz;b_R2.fq.gz; c_R2.fq.gz")
            check(pr == ["fq/a_R1.fq.gz", "fq/a_R2.fq.gz", "fq/b_R1.fq.gz", "fq/b_R2.fq.gz", "fq/c_R1.fq.gz", "fq/c_R2.fq.gz"],
                  "process-reads: whitespace stripped, ';' split, R1/R2 interleaved under the directory")

            print("\nkb ref / kb count / download-genome (binary optional)")
            has_kb = bool(shutil.which("kb"))
            cmd = kb_ref_command("kite", guide_features="guide_features.txt")
            check(cmd[:3] == ["kb", "ref", "-i"] and cmd[cmd.index("-f1") + 1] == "guide_mismatch.fa" and cmd[-3:] == ["--workflow", "kite", "guide_features.txt"],
                  "kb-ref kite: -i guide_index.idx -f1 guide_mismatch.fa -g t2guide.txt --workflow kite guide_features.txt")
            cmd = kb_ref_command("transcriptome", "mouse")
            check(cmd[:4] == ["kb", "ref", "-d", "mouse"] and "transcriptome_t2g.txt" in cmd, "kb-ref transcriptome: -d <species>")
            cc = kb_count_command("i.idx", "t2g.txt", "wl.txt", x_rna, "ks_transcripts_out", ["R1.fq", "R2.fq"])
            check(cc[cc.index("-x") + 1] == x_rna and cc[-5:] == ["R1.fq", "R2.fq", "--overwrite", "-m", "48G"] and cc[cc.index("-t") + 1] == "6",
                  "kb-count: -x from seqspec, fastqs then --overwrite -m 48G, 6 threads")
            rc = cmd_kb_ref(argparse.Namespace(mode="kite", species="human", guide_features=None, guide_table=str(S["xlsx"]), genome="genome.fa.gz",
                                               upstream_compat=True, out_dir=None, label="st_icp_kbref", dry_run=True))
            run = newest("st_icp_kbref")
            sref = json.loads((run / "summary.json").read_text())
            check(rc == 0 and sref["status"] == "printed" and sref["command"][sref["command"].index("-f1") + 1] == "genome.fa.gz"
                  and (run / "guide_features.txt").is_file(), "kb-ref --dry-run --upstream-compat: prints the command with the genome as -f1, exit 0")
            rc = cmd_kb_count(argparse.Namespace(modality="guide", index="guide_index.idx", t2g="t2guide.txt", whitelist=None, technology=None,
                                                 parsed_seqspec=str(newest("st_icp_seqspec") / "guide_parsed_seqSpec.txt"), fastqs=[str(S["fq1"]), str(S["fq2"])],
                                                 fastq_dir=None, reads1=None, reads2=None, threads=6, memory="30G", out_dir=None, label="st_icp_kbcount",
                                                 dry_run=True))
            run = newest("st_icp_kbcount")
            skc = json.loads((run / "summary.json").read_text())
            check(rc == 0 and skc["status"] == "printed" and skc["technology"] == x_g and skc["whitelist"] == "example_data/3M-february-2018.txt"
                  and skc["command"][skc["command"].index("-o") + 1].endswith("ks_guide_out"),
                  "kb-count: -x and -w read from the parsed seqspec; ks_guide_out; printed when kb is absent / --dry-run")
            rc = cmd_download_genome(argparse.Namespace(url=HG38_URL, out=str(d / "genome.fa.gz"), execute=False))
            check(rc == 0 and not (d / "genome.fa.gz").exists(), "download-genome: prints the wget command, nothing downloaded without --execute")

            print("\ncount-features (pure-Python kite)")
            exact, one = feature_lookup([("AAAA", "a"), ("AAAT", "b"), ("CCCC", "c")])
            check("AAAC" not in one and "AAAG" not in one and one.get("CCCA") == 2 and "AAAA" not in one,
                  "feature_lookup: Hamming-1 neighbours shared by two guides are dropped (kite collision rule)")
            fq_gz = [str(S["fq1"]) + ".gz", str(S["fq2"])]
            adata_g, st = count_features(fq_gz, x_g, [(s, i) for s, i in zip(S["guides"], S["ids"])], str(S["whitelist"]))
            exp = {(cb, S["ids"][g]): n for cb, gs in S["truth"].items() for g, n in gs.items()}
            got = {}
            Xc = adata_g.X.tocoo()
            for r, c, v in zip(Xc.row, Xc.col, Xc.data):
                got[(adata_g.obs_names[r], adata_g.var_names[c])] = int(v)
            check(got == exp, f"count-features: distinct-UMI counts equal the planted truth for all {len(exp)} cell-guide pairs")
            n_mm = sum(len(gs) for gs in S["truth"].values())
            check(st["feature_1mm"] == n_mm and st["barcodes_corrected_1mm"] == 50 and st["no_feature"] == 60,
                  f"count-features: {st['feature_1mm']} mismatched guide reads rescued, {st['barcodes_corrected_1mm']} barcodes corrected, 60 guide-less reads dropped")
            check(adata_g.shape == (50, 10) and list(adata_g.var_names) == S["ids"], "count-features: 50 guide-carrying barcodes x 10 features in guide_features order")
            adata_ex, st_ex = count_features(fq_gz, x_g, [(s, i) for s, i in zip(S["guides"], S["ids"])], None, mismatch=False)
            check(st_ex["feature_1mm"] == 0 and st_ex["feature_exact"] == st["feature_exact"] and st_ex["barcodes"] > 50,
                  "count-features --exact / no whitelist: mismatches not rescued, uncorrected barcodes kept")
            rc = cmd_count_features(argparse.Namespace(fastqs=[str(S["fq1"]), str(S["fq2"])], guide_features=None, guide_table=str(S["xlsx"]),
                                                       technology=None, parsed_seqspec=None, seqspec=str(guide_y), modality="guide", whitelist=str(S["whitelist"]),
                                                       no_whitelist=False, exact=False, both_strands=False, no_plots=args.no_plots, label="st_icp_count"))
            run = newest("st_icp_count")
            scf = json.loads((run / "summary.json").read_text())
            check(rc == 0 and scf["umis"] == sum(exp.values()) and (run / "counts_unfiltered" / "adata.h5ad").is_file() and scf["technology"] == x_g,
                  "count-features CLI: -x from the seqspec, counts_unfiltered/adata.h5ad written")

            print("\npreprocess")
            R = _synthetic_rna(d, S["cells"], rng)
            try:
                import scanpy  # type: ignore  # noqa: F401
                has_sc = True
            except ImportError:
                has_sc = False
            if has_sc:
                rc = cmd_preprocess(argparse.Namespace(adata_rna=str(R["h5ad"]), gene_names=str(R["genes"]), min_genes=100, min_cells=3, reference="human",
                                                       mt_prefix=None, no_plots=args.no_plots, label="st_icp_pre"))
                run = newest("st_icp_pre")
                sp_ = json.loads((run / "summary.json").read_text())
                fa = _load_adata(str(run / "filtered_anndata.h5ad"))
                X = R["X"]
                keep_c = (X > 0).sum(1) >= 100
                keep_g = (X[keep_c] > 0).sum(0) >= 3
                check(rc == 0 and fa.shape == (int(keep_c.sum()), int(keep_g.sum())) and sp_["cells_in"] - sp_["cells_out"] == len(R["low"]),
                      f"preprocess: filter_cells(min_genes=100) drops the {len(R['low'])} planted low-complexity barcodes; filter_genes(min_cells=3) -> {fa.shape[1]} genes")
                check(not any(v.startswith("RARE") for v in fa.var_names) and fa.var_names[0] == "MT-ND1", "preprocess: gene names applied; rare genes removed")
                Xf = X[keep_c][:, keep_g]
                gl = np.array(R["gene_list"])[keep_g]
                mt = np.char.startswith(gl.astype(str), "MT-")
                pct = 100 * Xf[:, mt].sum(1) / Xf.sum(1)
                check(int(fa.var["mt"].sum()) == 5 and not fa.var.loc["Mt-co1", "mt"] and int(fa.var["ribo"].sum()) == 8
                      and np.allclose(fa.obs["pct_counts_mt"].to_numpy(), pct, atol=1e-4),
                      "preprocess: human mt = MT- (5 genes, Mt-co1 excluded), ribo RPS/RPL (8), pct_counts_mt recomputed exactly")
                check({"batch_number", "n_genes", "n_genes_by_counts", "log1p_total_counts", "pct_counts_ribo", "total_counts_mt", "pct_counts_in_top_500_genes"} <= set(fa.obs.columns)
                      and (fa.obs["batch_number"] == 1).all() and {"n_cells", "mt", "ribo", "mean_counts", "pct_dropout_by_counts"} <= set(fa.var.columns),
                      "preprocess: calculate_qc_metrics columns (log1p=True), batch_number = 1")
                a2, _, info2 = preprocess_adata(_load_adata(str(R["h5ad"])), R["gene_list"], reference="mouse")
                check(info2["mt_prefix"] == "Mt-" and info2["n_mt_genes"] == 1, "preprocess: non-human reference uses the Mt- prefix (upstream rule)")
                try:
                    preprocess_adata(_load_adata(str(R["h5ad"])), R["gene_list"][:-1])
                    check(False, "preprocess: gene-name count mismatch raises")
                except ValueError:
                    check(True, "preprocess: gene-name count mismatch raises ValueError")
                if not args.no_plots:
                    check(len(sp_["figures"]) == 3, "preprocess: knee, violin and scatter figures")

                print("\ncreate-mdata")
                rc = cmd_create_mdata(argparse.Namespace(adata_rna=str(run / "filtered_anndata.h5ad"), adata_guide=str(newest("st_icp_count") / "counts_unfiltered" / "adata.h5ad"),
                                                         guide_metadata=str(S["xlsx"]), upstream_compat=False, no_plots=args.no_plots, label="st_icp_mdata"))
                run_m = newest("st_icp_mdata")
                sm = json.loads((run_m / "summary.json").read_text())
                expected_bc = [b for b in fa.obs_names if b in set(adata_g.obs_names)]
                g = read_h5mu_mod(run_m / "mudata.h5mu", "guides")
                t = read_h5mu_mod(run_m / "mudata.h5mu", "transcripts")
                check(rc == 0 and sm["intersecting_barcodes"] == len(expected_bc) and list(g.obs_names) == expected_bc == list(t.obs_names),
                      f"create-mdata: {len(expected_bc)} intersecting barcodes (guide cells passing RNA QC), RNA order kept, both modalities")
                check(list(g.var_names) == [f"{i}|{s}" for i, s in zip(S["ids"], S["guides"])] and sm["guide_naming"] == "by_id",
                      "create-mdata: guide features named sgRNA_ID|sgRNA_sequences")
                nz = {cb: len(gs) for cb, gs in S["truth"].items()}
                check(all(int(v) == nz[b] for b, v in zip(g.obs_names, g.obs["number_of_nonzero_guides"])) and (g.obs["batch_number"] == 1).all(),
                      "create-mdata: number_of_nonzero_guides = planted guides per cell (1 or 2); batch_number = 1")
                import h5py  # type: ignore
                with h5py.File(str(run_m / "mudata.h5mu"), "r") as f:
                    enc = f.attrs.get("encoding-type")
                    enc = enc.decode() if isinstance(enc, bytes) else enc
                    check(enc == "MuData" and set(f["mod"].keys()) == {"transcripts", "guides"}, f"create-mdata: h5mu with mod/transcripts + mod/guides ({sm['writer']})")
                a_rna, a_g = _load_adata(str(run / "filtered_anndata.h5ad")), _load_adata(str(newest("st_icp_count") / "counts_unfiltered" / "adata.h5ad"))
                a_g.var_names = [f"x{i}" for i in range(a_g.shape[1])]
                _, gs2, _, info_c = create_mdata(a_rna, a_g, S["table"], upstream_compat=False)
                check(info_c["guide_naming"] == "positional" and gs2.var_names[3] == f"{S['ids'][3]}|{S['guides'][3]}",
                      "create-mdata: positional naming when the matrix has no guide IDs (upstream behaviour)")
                perm = list(reversed(range(10)))
                a_g2 = _load_adata(str(newest("st_icp_count") / "counts_unfiltered" / "adata.h5ad"))[:, perm].copy()
                _, gs3, _, _ = create_mdata(_load_adata(str(run / "filtered_anndata.h5ad")), a_g2, S["table"], upstream_compat=False)
                a_g3 = _load_adata(str(newest("st_icp_count") / "counts_unfiltered" / "adata.h5ad"))[:, perm].copy()
                _, gs4, _, _ = create_mdata(_load_adata(str(run / "filtered_anndata.h5ad")), a_g3, S["table"], upstream_compat=True)
                check(gs3.var_names[0] == f"{S['ids'][9]}|{S['guides'][9]}" and gs4.var_names[0] == f"{S['ids'][0]}|{S['guides'][0]}",
                      "create-mdata: reordered matrix named by ID; --upstream-compat reproduces the positional (mis)naming")

                print("\nrun (end to end)")
                rc = cmd_run(argparse.Namespace(rna_seqspec=str(rna_y), guide_seqspec=str(guide_y), guide_table=str(S["xlsx"]), adata_guide=None,
                                                guide_fastqs=[str(S["fq1"]), str(S["fq2"])], technology=None, whitelist=None, no_whitelist=False,
                                                python_counts=True, adata_rna=str(R["h5ad"]), gene_names=str(R["genes"]), min_genes=100, min_cells=3,
                                                reference="human", mt_prefix=None, engine="python", upstream_compat=False, no_plots=True, label="st_icp_run"))
                run_r = newest("st_icp_run")
                sr = json.loads((run_r / "summary.json").read_text())
                check(rc == 0 and sr["steps"]["seqspec_guide"]["representation"] == x_g and sr["steps"]["count_features"]["umis"] == sum(exp.values())
                      and sr["steps"]["create_mdata"]["intersecting_barcodes"] == len(expected_bc) and (run_r / "mudata.h5mu").is_file()
                      and (run_r / "rna_parsed_seqSpec.txt").is_file(),
                      "run: seqspec -> guide features -> counts -> preprocess -> MuData, same numbers as the individual steps")
            else:
                print("  skip  preprocess / create-mdata / run: scanpy not installed")
    finally:
        leftovers = (set(OUT_ROOT.glob("*_st_icp_*")) - before) if OUT_ROOT.is_dir() else set()
        for p in set(made) | leftovers:
            shutil.rmtree(p, ignore_errors=True)
        if OUT_ROOT.is_dir() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[List[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent igvf-crispr-pipeline", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("parse-seqspec", help="seqspec YAML -> <modality>_parsed_seqSpec.txt with the kb -x string and whitelist")
    p.add_argument("--yaml", required=True, help="seqspec YAML (one modality per file, upstream convention)")
    p.add_argument("--modality", nargs="+", required=True, help="modalities to parse, e.g. rna guide multiseq")
    p.add_argument("--directory", help="seqspec directory joined to the onlist filename (default: the YAML's directory)")
    p.add_argument("--engine", choices=["auto", "python", "seqspec"], default="auto", help="auto: seqspec binary when on PATH, else the Python parser")
    p.add_argument("--strict-region-types", action="store_true", help="only region_type 'barcode' is a barcode (seqspec binary behaviour)")
    p.add_argument("--label", default="seqspec")
    p.set_defaults(func=cmd_parse_seqspec)

    p = sub.add_parser("extract-seqspec", help="parsed-seqspec TSV -> kb -x string")
    p.add_argument("--file", required=True)
    p.set_defaults(func=cmd_extract_seqspec)

    p = sub.add_parser("guide-features", help="guide table (xlsx/tsv/csv) -> guide_features.txt")
    p.add_argument("--guide-table", required=True, help="table with sgRNA_ID and sgRNA_sequences columns")
    p.add_argument("--out", help="output path (default: a run dir with report)")
    p.add_argument("--label", default="guide_features")
    p.set_defaults(func=cmd_guide_features)

    p = sub.add_parser("process-reads", help="';'-separated read lists -> interleaved kb fastq argument")
    p.add_argument("--dir", required=True)
    p.add_argument("--reads1", required=True)
    p.add_argument("--reads2", required=True)
    p.set_defaults(func=cmd_process_reads)

    p = sub.add_parser("download-genome", help="wget the genome FASTA the upstream kite step takes (~1 GB for hg38)")
    p.add_argument("--url", default=HG38_URL)
    p.add_argument("--out", default="genome.fa.gz")
    p.add_argument("--execute", action="store_true", help="actually download (otherwise print the command)")
    p.set_defaults(func=cmd_download_genome)

    p = sub.add_parser("kb-ref", help="kb ref: transcriptome index (-d) or guide kite index")
    p.add_argument("--mode", choices=["transcriptome", "kite"], default="transcriptome")
    p.add_argument("--species", default="human", help="kb ref -d value (transcriptome mode)")
    p.add_argument("--guide-features", help="guide_features.txt (kite)")
    p.add_argument("--guide-table", help="guide table; guide_features.txt is written from it (kite)")
    p.add_argument("--genome", help="genome FASTA passed as -f1 with --upstream-compat")
    p.add_argument("--upstream-compat", action="store_true", help="pass --genome as -f1 as the upstream does")
    p.add_argument("--out-dir")
    p.add_argument("--dry-run", action="store_true", help="print the command only")
    p.add_argument("--label", default="kb_ref")
    p.set_defaults(func=cmd_kb_ref)

    p = sub.add_parser("kb-count", help="kb count for rna or guide reads, -x from a parsed seqspec")
    p.add_argument("--modality", choices=["rna", "guide"], default="rna")
    p.add_argument("--index", required=True)
    p.add_argument("--t2g", required=True)
    p.add_argument("--parsed-seqspec", help="<modality>_parsed_seqSpec.txt (gives -x and -w)")
    p.add_argument("--technology", help="kb -x string (overrides --parsed-seqspec)")
    p.add_argument("--whitelist")
    p.add_argument("--fastqs", nargs="+", help="FASTQs in kb order (R1 R2 R1 R2 ...)")
    p.add_argument("--fastq-dir", help="with --reads1/--reads2 ';'-separated lists (MergeTwoWorkflows style)")
    p.add_argument("--reads1")
    p.add_argument("--reads2")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--memory", default="48G")
    p.add_argument("--out-dir")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--label", default="kb_count")
    p.set_defaults(func=cmd_kb_count)

    p = sub.add_parser("count-features", help="pure-Python kite-style guide counting -> counts_unfiltered/adata.h5ad")
    p.add_argument("--fastqs", nargs="+", required=True, help="FASTQs in technology order (R1 R2 [R1 R2 ...]); .gz allowed")
    p.add_argument("--guide-table", help="guide table with sgRNA_ID / sgRNA_sequences")
    p.add_argument("--guide-features", help="guide_features.txt (sequence<TAB>ID) instead of --guide-table")
    p.add_argument("--seqspec", help="seqspec YAML (gives -x and the onlist next to it)")
    p.add_argument("--modality", default="guide")
    p.add_argument("--parsed-seqspec")
    p.add_argument("--technology")
    p.add_argument("--whitelist", help="barcode onlist (plain or .gz)")
    p.add_argument("--no-whitelist", action="store_true", help="keep all barcodes uncorrected")
    p.add_argument("--exact", action="store_true", help="exact feature matches only (no Hamming-1 rescue)")
    p.add_argument("--both-strands", action="store_true", help="also try the reverse complement of the feature")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="count_features")
    p.set_defaults(func=cmd_count_features)

    def add_pre(p):
        p.add_argument("--gene-names", help="cells_x_genes.genes.names.txt; one name per var, first column")
        p.add_argument("--min-genes", type=int, default=100)
        p.add_argument("--min-cells", type=int, default=3)
        p.add_argument("--reference", default="human", help="'human' -> MT- prefix, anything else -> Mt-")
        p.add_argument("--mt-prefix", help="override the mitochondrial gene prefix (e.g. mt- for mouse)")

    p = sub.add_parser("preprocess", help="RNA AnnData QC (PreprocessAnnData)")
    p.add_argument("--adata-rna", required=True)
    add_pre(p)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="preprocess")
    p.set_defaults(func=cmd_preprocess)

    p = sub.add_parser("create-mdata", help="RNA + guide AnnData + guide metadata -> MuData (CreateMuData)")
    p.add_argument("--adata-rna", required=True, help="filtered_anndata.h5ad")
    p.add_argument("--adata-guide", required=True, help="guide counts h5ad (kb or count-features)")
    p.add_argument("--guide-metadata", required=True, help="guide table (xlsx/tsv) with sgRNA_ID, sgRNA_sequences")
    p.add_argument("--upstream-compat", action="store_true", help="positional guide naming and set-ordered barcodes")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="mudata")
    p.set_defaults(func=cmd_create_mdata)

    p = sub.add_parser("run", help="seqspec -> guide features -> guide counts -> preprocess -> MuData")
    p.add_argument("--rna-seqspec")
    p.add_argument("--guide-seqspec")
    p.add_argument("--guide-table")
    p.add_argument("--adata-guide", help="existing guide counts h5ad (skips counting)")
    p.add_argument("--guide-fastqs", nargs="+")
    p.add_argument("--technology", help="guide kb -x string when no --guide-seqspec")
    p.add_argument("--whitelist")
    p.add_argument("--no-whitelist", action="store_true")
    p.add_argument("--python-counts", action="store_true", help="use count-features even if kb + kallisto are installed")
    p.add_argument("--adata-rna", help="RNA counts h5ad (e.g. kb count counts_unfiltered/adata.h5ad)")
    add_pre(p)
    p.add_argument("--engine", choices=["auto", "python", "seqspec"], default="auto")
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="igvf_crispr_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="upstream example seqspecs + synthetic reads/AnnData, all subcommands asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
