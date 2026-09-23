"""AlphaGenome predictions and Atlas scores from inside IGVFagent (wraps google-deepmind/alphagenome).

https://github.com/google-deepmind/alphagenome — the official AlphaGenome API
client (Apache-2.0). Relationship: **wraps-sdk**. AlphaGenome is a hosted model
run by Google DeepMind; nothing here reimplements it. This skill turns the
SDK's calls into IGVFagent commands that take the notations IGVF users already
have (rsIDs, SPDI, chr-pos-ref-alt, gene symbols), resolve them through the
IGVF Catalog, and write the usual run folder: report.md, summary.json, TSVs
and figures, with every artefact announced on stdout.

Access and terms (printed in every report)
    An API key is required: https://alphagenome.google/api (free for
    non-commercial use). Outputs are for non-commercial research only, must
    not be used to train other models, and must not be used for clinical
    decisions (https://alphagenome.google/terms). The key is read from
    ALPHAGENOME_API_KEY, or from Docs/Secret/ALPHAGENOME_API_KEY.txt, and is
    never printed or logged.

Python
    The SDK needs Python >= 3.10; IGVFagent supports 3.9. When the SDK is not
    importable here, every API subcommand re-runs itself under a 3.10+
    interpreter that has it: $IGVF_ALPHAGENOME_PYTHON, else
    ~/.igvfagent/alphagenome-venv (created by `igvfagent alphagenome setup`).

Definitions (from the SDK, alphagenome.models.dna_client / variant_scorers)
    Sequence lengths      16KB=16,384  100KB=131,072  500KB=524,288  1MB=1,048,576
                          (a prediction window must be one of these; regions
                          are centred and resized to it)
    Coordinates           Interval: 0-based, half-open (like BED and the IGVF
                          Catalog's gene start/end). Variant: 1-based position.
                          SPDI and Catalog variant positions are 0-based, so
                          AlphaGenome position = SPDI position + 1.
    Output types          ATAC CAGE DNASE RNA_SEQ CHIP_HISTONE CHIP_TF
                          SPLICE_SITES SPLICE_SITE_USAGE SPLICE_JUNCTIONS
                          CONTACT_MAPS PROCAP
    Organisms             HOMO_SAPIENS (hg38), MUS_MUSCULUS (mm10)
    Recommended scorers   variant_scorers.RECOMMENDED_VARIANT_SCORERS, e.g.
                          DNASE/ATAC/CHIP_TF/CAGE/PROCAP = CenterMask 501 bp
                          DIFF_LOG2_SUM, CHIP_HISTONE = 2,001 bp, RNA_SEQ =
                          GeneMaskLFC, SPLICE_* and POLYADENYLATION; *_ACTIVE
                          variants use ACTIVE_SUM. At most 20 per request.
    ISM                   in silico mutagenesis: every alternative base at
                          every position of --ism-interval (3 variants/bp),
                          sent in 10-bp chunks by the SDK.
    Atlas                 pre-computed variant scores (incl. AVI, the
                          AlphaGenome Variant Impact score) for the human
                          genome; queried per variant or per interval.

Subcommands
    setup             check the SDK, the interpreter and the key; --install
                      creates the 3.10+ environment; --ping calls the API
    metadata          every predicted track (output type, assay, ontology
                      CURIE, biosample); --search to find e.g. heart tracks
    predict-interval  predicted tracks over a region or gene
    predict-variant   REF vs ALT tracks around a variant
    score-variants    recommended (or chosen) variant scores, tidy table
    score-interval    interval scores (gene-level) for a region
    ism               in silico mutagenesis over a short interval
    atlas-scorers     the Atlas scorers available
    atlas-variants    pre-computed Atlas scores for variants
    atlas-interval    pre-computed Atlas scores for every variant in a region
    selftest          offline checks with a fake backend (no key, no network)

Usage:
    igvfagent alphagenome setup --install --ping
    igvfagent alphagenome metadata --search heart
    igvfagent alphagenome score-variants --variants rs429358 chr22:36201698:A>C --ontology UBERON:0000948
    igvfagent alphagenome predict-variant --variant rs429358 --outputs RNA_SEQ DNASE --ontology UBERON:0000948
    igvfagent alphagenome predict-interval --gene APOE --outputs RNA_SEQ --ontology UBERON:0002107
    igvfagent alphagenome ism --ism-interval chr19:44908670-44908700 --scorer DNASE --ontology UBERON:0000948
    igvfagent alphagenome atlas-variants --variants rs429358 --genes APOE
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
OUT_ROOT = ROOT / "Docs" / "AlphaGenome"
LOG_DIR = ROOT / "Docs" / "Logs"
CACHE_DIR = ROOT / "Data" / "Cache" / "alphagenome"

KEY_ENV = ("ALPHAGENOME_API_KEY",)
KEY_FILE = ROOT / "Docs" / "Secret" / "ALPHAGENOME_API_KEY.txt"
PY_ENV = "IGVF_ALPHAGENOME_PYTHON"
VENV = Path(os.environ.get("IGVF_ALPHAGENOME_VENV")
            or Path.home() / ".igvfagent" / "alphagenome-venv")
REEXEC_FLAG = "IGVF_ALPHAGENOME_CHILD"

SEQ_LENGTHS = {"16KB": 2 ** 14, "100KB": 2 ** 17, "500KB": 2 ** 19, "1MB": 2 ** 20}
OUTPUT_TYPES = ("ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
                "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
                "CONTACT_MAPS", "PROCAP")
ORGANISMS = ("HOMO_SAPIENS", "MUS_MUSCULUS")
MAX_SCORERS = 20
MAX_VARIANTS = 1000          # per call; the API is meant for 1000s, not 1M+
MAX_ISM_BP = 100             # 3 variants per bp
MAX_ATLAS_BP = 20000

TERMS = ("AlphaGenome outputs are for non-commercial research use only, must not be "
         "used to train other models, and must not be used for clinical "
         "decision-making (https://alphagenome.google/terms). Cite Avsec et al., "
         "Nature 2026 (https://www.nature.com/articles/s41586-025-10014-0).")

logger = logging.getLogger("alphagenome")


# ─── plumbing ───────────────────────────────────────────────────────────────

def safe_label(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "run").strip("_")[:60] or "run"


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def setup_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    p = LOG_DIR / f"alphagenome_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(p), logging.StreamHandler(sys.stderr)])
    print(f"Log: {p}")


def api_key() -> Optional[str]:
    for k in KEY_ENV:
        v = (os.environ.get(k) or "").strip()
        if v:
            return v
    try:
        v = KEY_FILE.read_text().strip()
        return v or None
    except OSError:
        return None


def key_help() -> str:
    return ("No AlphaGenome API key. Get one (free for non-commercial use) at "
            "https://alphagenome.google/api, then either `export "
            "ALPHAGENOME_API_KEY=...` or save it to "
            "Docs/Secret/ALPHAGENOME_API_KEY.txt (gitignored). On the hosted "
            "deployment the operator adds ALPHAGENOME_API_KEY to Deploy/.env.prod.")


def sdk_available() -> bool:
    if sys.version_info < (3, 10):
        return False
    try:
        import alphagenome  # type: ignore  # noqa: F401
        from alphagenome.models import dna_client  # type: ignore  # noqa: F401
        return True
    except Exception:
        return False


def sdk_python() -> Optional[str]:
    """A 3.10+ interpreter that has the SDK, for re-running this module."""
    cand = os.environ.get(PY_ENV)
    if cand and Path(cand).exists():
        return cand
    for name in ("bin/python", "Scripts/python.exe"):
        p = VENV / name
        if p.exists():
            return str(p)
    return None


def reexec(argv: Sequence[str]) -> int:
    """Run this subcommand under the SDK interpreter; stream its output."""
    py = sdk_python()
    if not py:
        print("The AlphaGenome SDK needs Python >= 3.10 and is not available to "
              f"this interpreter (Python {sys.version.split()[0]}). Run "
              "`igvfagent alphagenome setup --install` once, or point "
              f"{PY_ENV} at a Python 3.10+ that has `pip install alphagenome`.")
        return 2
    env = dict(os.environ)
    env[REEXEC_FLAG] = "1"
    env["IGVF_PROJECT_ROOT"] = str(ROOT)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parent) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.call([py, str(Path(__file__).resolve()), *argv], env=env)


# ─── IGVF-side resolution: variants, intervals, genes ───────────────────────

_CATALOG = None


def catalog_base() -> str:
    global _CATALOG
    if _CATALOG is None:
        try:
            import _endpoints  # type: ignore
            _CATALOG = _endpoints.resolve("catalog_api").rstrip("/")
        except Exception:
            _CATALOG = "https://api.catalogkg.igvf.org"
    return _CATALOG


def _get_json(url: str, timeout: int = 60) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "IGVFagent-alphagenome"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


_NC_CHROM = {**{f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)},
             "NC_000023": "chrX", "NC_000024": "chrY", "NC_012920": "chrM"}
_DELIM = re.compile(r"^(?:chr)?([0-9]{1,2}|[XYM]|MT)[-:_]g?\.?(\d+)[-:_]([ACGTN]+)[-:_>]([ACGTN]+)$", re.I)
_HGVS = re.compile(r"^(?:chr)?([0-9]{1,2}|[XYM]|MT)[-:_]g?\.?(\d+)([ACGTN]+)>([ACGTN]+)$", re.I)
_SPDI = re.compile(r"^(NC_\d{6})\.\d+:(\d+):([ACGTN]*):([ACGTN]*)$", re.I)
_RSID = re.compile(r"^rs\d+$", re.I)
_REGION = re.compile(r"^(?:chr)?([0-9]{1,2}|[XYM]|MT):([\d,]+)-([\d,]+)$", re.I)


def _chrom(c: str) -> str:
    c = c.upper()
    c = c[3:] if c.startswith("CHR") else c
    return "chrM" if c in ("M", "MT") else f"chr{c}"


def parse_variant(token: str, resolve_rsid=None) -> Dict[str, Any]:
    """{chrom, position (1-based), ref, alt, raw}. rsIDs go through the Catalog."""
    t = (token or "").strip().strip(",;\"'")
    m = _DELIM.match(t) or _HGVS.match(t)
    if m:
        c, pos, ref, alt = m.groups()
        return {"raw": t, "chrom": _chrom(c), "position": int(pos),
                "ref": ref.upper(), "alt": alt.upper()}
    m = _SPDI.match(t)
    if m:
        acc, pos0, ref, alt = m.groups()
        if acc.upper() not in _NC_CHROM:
            raise ValueError(f"unknown RefSeq chromosome in {t!r}")
        return {"raw": t, "chrom": _NC_CHROM[acc.upper()], "position": int(pos0) + 1,
                "ref": ref.upper(), "alt": alt.upper()}
    if _RSID.match(t):
        return (resolve_rsid or catalog_rsid)(t)
    raise ValueError(f"unrecognised variant notation: {t!r} (use rsID, SPDI, "
                     "chr-pos-ref-alt or chr:pos:ref>alt)")


def catalog_rsid(rsid: str) -> Dict[str, Any]:
    url = f"{catalog_base()}/api/variants?{urllib.parse.urlencode({'rsid': rsid.lower()})}"
    data = _get_json(url)
    rows = data if isinstance(data, list) else data.get("data", [])
    if not rows:
        raise ValueError(f"{rsid} not found in the IGVF Catalog")
    alts = [r for r in rows if r.get("ref") and r.get("alt")]
    if len(alts) > 1:
        logger.info("%s has %d alleles in the Catalog; using %s>%s (list the others "
                    "explicitly to score them)", rsid, len(alts), alts[0]["ref"], alts[0]["alt"])
    r = alts[0] if alts else rows[0]
    return {"raw": rsid, "chrom": r["chr"], "position": int(r["pos"]) + 1,
            "ref": r["ref"].upper(), "alt": r["alt"].upper(), "rsid": rsid.lower()}


def catalog_gene(symbol: str, organism: str = "Homo sapiens") -> Dict[str, Any]:
    url = (f"{catalog_base()}/api/genes?"
           + urllib.parse.urlencode({"name": symbol, "organism": organism}))
    data = _get_json(url)
    rows = data if isinstance(data, list) else data.get("data", [])
    rows = [r for r in rows if (r.get("name") or "").upper() == symbol.upper()] or rows
    if not rows:
        raise ValueError(f"gene {symbol!r} not found in the IGVF Catalog")
    r = rows[0]
    return {"gene_id": r.get("_id"), "name": r.get("name"), "chrom": r["chr"],
            "start": int(r["start"]), "end": int(r["end"]), "strand": r.get("strand")}


def parse_region(text: str) -> Tuple[str, int, int]:
    """chr:start-end, 0-based half-open (BED / AlphaGenome Interval.from_str)."""
    m = _REGION.match((text or "").strip())
    if not m:
        raise ValueError(f"bad region {text!r}; use chr:start-end (0-based, half-open)")
    c, s, e = m.groups()
    s, e = int(s.replace(",", "")), int(e.replace(",", ""))
    if e <= s:
        raise ValueError(f"region end must be > start in {text!r}")
    return _chrom(c), s, e


def window(chrom: str, start: int, end: int, length: int) -> Tuple[str, int, int]:
    """Centre a region in a supported-length prediction window."""
    if end - start > length:
        raise ValueError(f"region is {end - start:,} bp; the largest window is "
                         f"{max(SEQ_LENGTHS.values()):,} bp")
    centre = (start + end) // 2
    s = max(0, centre - length // 2)
    return chrom, s, s + length


def read_variants(tokens: Sequence[str], path: Optional[str],
                  resolve_rsid=None) -> Tuple[List[dict], List[dict]]:
    raw: List[str] = list(tokens or [])
    if path:
        for line in Path(path).read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            f = line.split()
            if len(f) >= 5 and f[1].isdigit() and re.fullmatch(r"[ACGTN]+", f[3], re.I):
                raw.extend(f"{f[0]}-{f[1]}-{f[3]}-{a}" for a in f[4].split(","))  # VCF
            else:
                raw.extend(t for t in re.split(r"[\s,]+", line) if t)
    out, bad, seen = [], [], set()
    for t in raw:
        try:
            v = parse_variant(t, resolve_rsid)
        except Exception as e:  # report, never drop silently
            bad.append({"raw": t, "error": str(e)})
            continue
        key = (v["chrom"], v["position"], v["ref"], v["alt"])
        if key not in seen:
            seen.add(key)
            out.append(v)
    return out, bad


def vid(v: dict) -> str:
    return f"{v['chrom']}:{v['position']}:{v['ref']}>{v['alt']}"


def resolve_region(args) -> Tuple[str, int, int, str]:
    """(chrom, start, end, label) from --interval or --gene."""
    if getattr(args, "gene", None):
        g = catalog_gene(args.gene, "Mus musculus" if args.organism == "MUS_MUSCULUS"
                         else "Homo sapiens")
        return g["chrom"], g["start"], g["end"], f"{g['name']} ({g['gene_id']})"
    if getattr(args, "interval", None):
        c, s, e = parse_region(args.interval)
        return c, s, e, f"{c}:{s}-{e}"
    raise SystemExit("give --interval chr:start-end or --gene SYMBOL")


# ─── backend: the only code that touches the SDK ────────────────────────────

class SdkBackend:
    """Thin adapter over alphagenome; returns pandas/numpy, never SDK objects."""

    def __init__(self, key: str):
        from alphagenome.models import dna_client  # type: ignore
        self._dna_client = dna_client
        self._key = key
        self._model = None
        self._atlas = None

    @property
    def model(self):
        if self._model is None:
            self._model = self._dna_client.create(self._key)
        return self._model

    @property
    def atlas(self):
        if self._atlas is None:
            from alphagenome.atlas import atlas  # type: ignore
            self._atlas = atlas.create(self._key)
        return self._atlas

    # helpers
    def _org(self, organism: str):
        return getattr(self._dna_client.Organism, organism)

    def _outs(self, names: Sequence[str]):
        from alphagenome.models import dna_output  # type: ignore
        return [getattr(dna_output.OutputType, n) for n in names]

    @staticmethod
    def _interval(c, s, e):
        from alphagenome.data import genome  # type: ignore
        return genome.Interval(chromosome=c, start=s, end=e)

    @staticmethod
    def _variant(v):
        from alphagenome.data import genome  # type: ignore
        return genome.Variant(chromosome=v["chrom"], position=v["position"],
                              reference_bases=v["ref"], alternate_bases=v["alt"])

    def _scorers(self, names: Sequence[str], organism: str):
        from alphagenome.models import variant_scorers as vs  # type: ignore
        if not names:
            return vs.get_recommended_scorers(self._org(organism).to_proto())
        bad = [n for n in names if n not in vs.RECOMMENDED_VARIANT_SCORERS]
        if bad:
            raise ValueError(f"unknown scorer(s) {bad}; choose from "
                             f"{sorted(vs.RECOMMENDED_VARIANT_SCORERS)}")
        return [vs.RECOMMENDED_VARIANT_SCORERS[n] for n in names]

    @staticmethod
    def _tracks(output) -> Dict[str, dict]:
        out = {}
        for name in OUTPUT_TYPES:
            td = getattr(output, name.lower(), None)
            if td is None or getattr(td, "values", None) is None:
                continue
            out[name] = {"values": td.values, "metadata": td.metadata.reset_index(drop=True),
                         "resolution": int(getattr(td, "resolution", 1) or 1)}
        return out

    # API
    def metadata(self, organism: str):
        df = self.model.output_metadata(self._org(organism)).concatenate()
        df["output_type"] = df["output_type"].map(lambda x: getattr(x, "name", str(x)))
        return df.reset_index(drop=True)

    def predict_interval(self, c, s, e, outputs, ontology, organism):
        o = self.model.predict_interval(self._interval(c, s, e), organism=self._org(organism),
                                        requested_outputs=self._outs(outputs),
                                        ontology_terms=list(ontology) or None)
        return self._tracks(o)

    def predict_variant(self, c, s, e, v, outputs, ontology, organism):
        o = self.model.predict_variant(self._interval(c, s, e), self._variant(v),
                                       organism=self._org(organism),
                                       requested_outputs=self._outs(outputs),
                                       ontology_terms=list(ontology) or None)
        return self._tracks(o.reference), self._tracks(o.alternate)

    def score_variants(self, windows, variants, scorers, organism):
        from alphagenome.models import variant_scorers as vs  # type: ignore
        res = self.model.score_variants([self._interval(*w) for w in windows],
                                        [self._variant(v) for v in variants],
                                        variant_scorers=self._scorers(scorers, organism),
                                        organism=self._org(organism), progress_bar=False)
        return vs.tidy_scores(res)

    def score_interval(self, c, s, e, organism):
        from alphagenome.models import variant_scorers as vs  # type: ignore
        res = self.model.score_interval(self._interval(c, s, e), organism=self._org(organism))
        return vs.tidy_scores(res)

    def ism(self, win, ism_iv, scorers, organism):
        from alphagenome.models import variant_scorers as vs  # type: ignore
        res = self.model.score_ism_variants(self._interval(*win), self._interval(*ism_iv),
                                            variant_scorers=self._scorers(scorers, organism),
                                            organism=self._org(organism), progress_bar=False)
        return vs.tidy_scores(res)

    def atlas_scorers(self):
        import pandas as pd
        rows = [{"scorer": k, "is_signed": m.is_signed,
                 "n_tracks": int(len(m.track_metadata))}
                for k, m in self.atlas.scorer_metadata().items()]
        return pd.DataFrame(rows)

    def atlas_variants(self, variants, scorers, ontology, genes):
        res = self.atlas.query_variants([self._variant(v) for v in variants],
                                        requested_scorers=list(scorers),
                                        ontology_terms=list(ontology) or None,
                                        gene_names=list(genes) or None, progress_bar=False)
        return flatten_atlas(res)

    def atlas_interval(self, c, s, e, scorers, ontology, genes):
        res = self.atlas.query_interval(self._interval(c, s, e),
                                        requested_scorers=list(scorers),
                                        ontology_terms=list(ontology) or None,
                                        gene_names=list(genes) or None, progress_bar=False)
        return flatten_atlas(res)


def flatten_atlas(res) -> "Any":
    """{scorer: AnnData(variants[x genes] x tracks)} -> one long table."""
    import numpy as np
    import pandas as pd
    frames = []
    keep = [c for c in ("name", "ontology_curie", "biosample_name", "biosample_type",
                        "Assay title", "strand", "gtex_tissue", "transcription_factor",
                        "histone_mark")]
    for scorer, ad in res.items():
        X = np.asarray(ad.X)
        Q = np.asarray(ad.layers["quantiles"]) if "quantiles" in getattr(ad, "layers", {}) else None
        obs = ad.obs.reset_index(drop=True) if ad.obs is not None and len(ad.obs.columns) else \
            pd.DataFrame({"variant": [None] * X.shape[0]})
        var = ad.var.reset_index(drop=True)
        vcols = [c for c in keep if c in var.columns]
        for j in range(X.shape[1]):
            f = obs.copy()
            f["variant"] = f["variant"].map(lambda x: str(x) if x is not None else None)
            f["scorer"] = scorer
            for c in vcols:
                f["track_" + c.replace(" ", "_").lower() if c != "name" else "track_name"] = var.at[j, c]
            f["raw_score"] = X[:, j]
            if Q is not None:
                f["quantile_score"] = Q[:, j]
            frames.append(f)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def get_backend(offline: bool = False):
    if offline or os.environ.get("IGVF_ALPHAGENOME_FAKE") == "1":
        return FakeBackend()
    key = api_key()
    if not key:
        raise SystemExit(key_help())
    return SdkBackend(key)


# ─── outputs ────────────────────────────────────────────────────────────────

def write_report(d: Path, title: str, lines: List[str], summary: dict) -> None:
    (d / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    body = [f"# {title}", "", *lines, "", "---", "", f"_{TERMS}_", ""]
    (d / "report.md").write_text("\n".join(body))
    print(f"JSON: {d / 'summary.json'}")
    print(f"Report: {d / 'report.md'}")


def write_tsv(df, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)
    print(f"TSV: {path}")


def md_table(df, n: int = 15) -> List[str]:
    if df is None or not len(df):
        return ["(none)"]
    df = df.head(n)
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:.4g}" if isinstance(v, float) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return out


def track_label(meta_row) -> str:
    parts = [str(meta_row.get("name", ""))]
    for c in ("biosample_name", "Assay title", "strand"):
        v = meta_row.get(c)
        if v is not None and str(v) not in ("", "nan", "."):
            parts.append(str(v))
    return " · ".join(dict.fromkeys(p for p in parts if p))


def summarise_tracks(tracks: Dict[str, dict], win_start: int, lo: int, hi: int):
    """Per-track mean / max / sum over [lo, hi) (genome coordinates)."""
    import numpy as np
    import pandas as pd
    rows = []
    for ot, t in tracks.items():
        vals, res = np.asarray(t["values"], dtype=float), t["resolution"]
        if vals.ndim != 2:            # contact maps are 3-D; summarise the diagonal band
            continue
        a, b = max(0, (lo - win_start) // res), max(1, -(-(hi - win_start) // res))
        seg = vals[a:b]
        for j in range(vals.shape[1]):
            m = t["metadata"].iloc[j].to_dict() if j < len(t["metadata"]) else {}
            col = seg[:, j]
            rows.append({"output_type": ot, "track": m.get("name", j),
                         "ontology_curie": m.get("ontology_curie"),
                         "biosample_name": m.get("biosample_name"),
                         "assay": m.get("Assay title"), "strand": m.get("strand"),
                         "mean": float(col.mean()) if col.size else float("nan"),
                         "max": float(col.max()) if col.size else float("nan"),
                         "sum": float(col.sum()) if col.size else float("nan")})
    return pd.DataFrame(rows)


def plot_tracks(tracks: Dict[str, dict], win_start: int, lo: int, hi: int, path: Path,
                title: str, alt: Optional[Dict[str, dict]] = None, max_tracks: int = 8,
                mark: Optional[int] = None) -> Optional[Path]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return None
    panels = []
    for ot, t in tracks.items():
        vals = np.asarray(t["values"], dtype=float)
        if vals.ndim != 2:
            continue
        for j in range(min(vals.shape[1], max_tracks - len(panels))):
            panels.append((ot, j))
        if len(panels) >= max_tracks:
            break
    if not panels:
        return None
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 1.6 * len(panels) + 0.6),
                             sharex=True, squeeze=False)
    for ax, (ot, j) in zip(axes[:, 0], panels):
        t = tracks[ot]
        res = t["resolution"]
        a, b = max(0, (lo - win_start) // res), max(1, -(-(hi - win_start) // res))
        x = win_start + np.arange(a, a + len(np.asarray(t["values"])[a:b])) * res
        ax.fill_between(x, np.asarray(t["values"], dtype=float)[a:b, j], color="#4C72B0",
                        alpha=0.75, lw=0, label="REF" if alt else None, step="post")
        if alt and ot in alt:
            ax.plot(x, np.asarray(alt[ot]["values"], dtype=float)[a:b, j], color="#C44E52",
                    lw=0.9, label="ALT", drawstyle="steps-post")
        if mark is not None:
            ax.axvline(mark, color="k", lw=0.6, ls="--")
        meta = t["metadata"].iloc[j].to_dict() if j < len(t["metadata"]) else {}
        ax.set_ylabel(ot, fontsize=7)
        ax.set_title(track_label(meta), fontsize=7, loc="left")
        ax.tick_params(labelsize=6)
    if alt:
        axes[0, 0].legend(fontsize=6, loc="upper right")
    axes[-1, 0].set_xlabel("position (bp)", fontsize=7)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def filter_ontology(df, ontology: Sequence[str], search: Optional[str] = None):
    if df is None or not len(df):
        return df
    if ontology and "ontology_curie" in df.columns:
        df = df[df["ontology_curie"].isin(list(ontology))]
    if search:
        s = search.lower()
        cols = [c for c in ("biosample_name", "ontology_curie", "track_name", "gtex_tissue")
                if c in df.columns]
        mask = False
        for c in cols:
            mask = mask | df[c].astype(str).str.lower().str.contains(s, regex=False)
        df = df[mask]
    return df


def top_scores(df, n: int = 25):
    if df is None or not len(df):
        return df
    key = "quantile_score" if "quantile_score" in df.columns and df["quantile_score"].notna().any() \
        else "raw_score"
    return df.reindex(df[key].abs().sort_values(ascending=False).index).head(n)


# ─── commands ───────────────────────────────────────────────────────────────

def cmd_setup(args) -> int:
    ok = True
    print(f"Python here:        {sys.version.split()[0]}")
    here = sdk_available()
    print(f"SDK importable here: {'yes' if here else 'no'}")
    py = sdk_python()
    if not here:
        print(f"SDK interpreter:    {py or '(none)'}")
    if not here and not py and args.install:
        py = install_venv()
        ok = py is not None
    if not here and not py:
        print("  -> run `igvfagent alphagenome setup --install` to create a Python "
              f"3.10+ environment at {VENV} with the SDK.")
        ok = False
    key = api_key()
    print(f"API key:            {'found' if key else 'missing'}"
          f"{' (' + ('environment' if any(os.environ.get(k) for k in KEY_ENV) else 'Docs/Secret file') + ')' if key else ''}")
    if not key:
        print("  -> " + key_help())
        ok = False
    if args.ping and key and (here or py):
        if here:
            b = SdkBackend(key)
            n = len(b.metadata("HOMO_SAPIENS"))
            print(f"API reachable:      yes ({n:,} human tracks)")
        else:
            rc = reexec(["setup", "--ping"])
            ok = ok and rc == 0
    print("\n" + TERMS)
    return 0 if ok else 2


def install_venv() -> Optional[str]:
    uv = shutil.which("uv")
    VENV.parent.mkdir(parents=True, exist_ok=True)
    try:
        if uv:
            subprocess.check_call([uv, "venv", "-q", "-p", "3.11", str(VENV)])
            subprocess.check_call([uv, "pip", "install", "-q", "--python",
                                   str(VENV / "bin" / "python"), "alphagenome", "matplotlib"])
        else:
            base = next((shutil.which(f"python3.{m}") for m in (13, 12, 11, 10)
                         if shutil.which(f"python3.{m}")), None)
            if not base:
                print("Neither `uv` nor a python3.10+ is on PATH; install one, or set "
                      f"{PY_ENV} to a Python 3.10+ with `pip install alphagenome`.")
                return None
            subprocess.check_call([base, "-m", "venv", str(VENV)])
            subprocess.check_call([str(VENV / "bin" / "python"), "-m", "pip", "install",
                                   "-q", "alphagenome", "matplotlib"])
    except subprocess.CalledProcessError as e:
        print(f"install failed: {e}")
        return None
    print(f"Wrote: {VENV}")
    return sdk_python()


def cmd_metadata(args) -> int:
    b = get_backend(args.offline)
    df = b.metadata(args.organism)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if args.output_type:
        df = df[df["output_type"].astype(str).str.upper().str.contains("|".join(args.output_type))]
    shown = filter_ontology(df, args.ontology or [], args.search)
    d = run_dir(args.label or f"metadata_{args.search or args.organism}")
    write_tsv(shown, d / "tracks.tsv")
    counts = shown.groupby("output_type").size().rename("tracks").reset_index() if len(shown) else shown
    bios = (shown.groupby(["ontology_curie", "biosample_name"]).size().rename("tracks")
            .reset_index().sort_values("tracks", ascending=False)) if len(shown) and \
        "biosample_name" in shown.columns else shown
    write_tsv(bios, d / "biosamples.tsv")
    lines = [f"- Organism: {args.organism}", f"- Tracks matching: {len(shown):,} of {len(df):,}",
             f"- Filter: search={args.search!r} ontology={args.ontology or []}", "",
             "## Tracks per output type", "", *md_table(counts, 20), "",
             "## Biosamples (use the ontology CURIE in --ontology)", "", *md_table(bios, 25)]
    write_report(d, "AlphaGenome track metadata", lines,
                 {"organism": args.organism, "n_tracks": int(len(shown)),
                  "n_total": int(len(df)), "search": args.search})
    return 0


def cmd_predict_interval(args) -> int:
    c, s, e, label = resolve_region(args)
    L = SEQ_LENGTHS[args.length]
    w = window(c, s, e, L)
    if not args.ontology and not args.all_tracks:
        raise SystemExit("pass --ontology CURIE(s) (see `alphagenome metadata --search`) "
                         "or --all-tracks; all tracks for a 1 Mb window is a very large response")
    b = get_backend(args.offline)
    tracks = b.predict_interval(*w, args.outputs, args.ontology or [], args.organism)
    d = run_dir(args.label or f"predict_{label.split()[0]}")
    summ = summarise_tracks(tracks, w[1], s, e)
    write_tsv(summ, d / "track_summary.tsv")
    plot_tracks(tracks, w[1], s, e, d / "tracks.png", f"AlphaGenome · {label}")
    top = summ.sort_values("mean", ascending=False) if len(summ) else summ
    lines = [f"- Region: {label} = {c}:{s:,}-{e:,} ({e - s:,} bp)",
             f"- Prediction window: {w[0]}:{w[1]:,}-{w[2]:,} ({args.length})",
             f"- Outputs: {', '.join(args.outputs)}; ontology: {', '.join(args.ontology or ['all'])}",
             f"- Tracks returned: {len(summ):,}", "",
             "## Highest mean signal over the region", "",
             *md_table(top[["output_type", "track", "biosample_name", "mean", "max"]]
                       if len(top) else top)]
    write_report(d, f"AlphaGenome prediction · {label}", lines,
                 {"region": [c, s, e], "window": list(w), "outputs": args.outputs,
                  "ontology": args.ontology, "n_tracks": int(len(summ))})
    return 0


def cmd_predict_variant(args) -> int:
    vs, bad = read_variants([args.variant], None)
    if not vs:
        raise SystemExit(f"could not parse --variant: {bad[0]['error'] if bad else args.variant}")
    v = vs[0]
    L = SEQ_LENGTHS[args.length]
    w = window(v["chrom"], v["position"] - 1, v["position"], L)
    if not args.ontology and not args.all_tracks:
        raise SystemExit("pass --ontology CURIE(s) or --all-tracks")
    b = get_backend(args.offline)
    ref, alt = b.predict_variant(*w, v, args.outputs, args.ontology or [], args.organism)
    lo, hi = v["position"] - 1 - args.radius, v["position"] + args.radius
    rs, as_ = summarise_tracks(ref, w[1], lo, hi), summarise_tracks(alt, w[1], lo, hi)
    import numpy as np
    m = rs.merge(as_, on=["output_type", "track", "ontology_curie", "biosample_name", "assay",
                          "strand"], suffixes=("_ref", "_alt"))
    if len(m):
        m["log2fc_sum"] = np.log2((m["sum_alt"] + 1) / (m["sum_ref"] + 1))
        m["diff_sum"] = m["sum_alt"] - m["sum_ref"]
        m = m.reindex(m["log2fc_sum"].abs().sort_values(ascending=False).index)
    d = run_dir(args.label or f"variant_{vid(v)}")
    write_tsv(m, d / "ref_alt_summary.tsv")
    plot_tracks(ref, w[1], lo, hi, d / "ref_alt.png", f"AlphaGenome · {vid(v)}",
                alt=alt, mark=v["position"])
    lines = [f"- Variant: {vid(v)} (input {v['raw']})",
             f"- Window: {w[0]}:{w[1]:,}-{w[2]:,} ({args.length}); summarised within "
             f"±{args.radius:,} bp", f"- Outputs: {', '.join(args.outputs)}; ontology: "
             f"{', '.join(args.ontology or ['all'])}", "",
             "## Largest ALT vs REF changes (log2 fold change of the summed signal)", "",
             *md_table(m[["output_type", "track", "biosample_name", "sum_ref", "sum_alt",
                          "log2fc_sum"]] if len(m) else m)]
    write_report(d, f"AlphaGenome variant prediction · {vid(v)}", lines,
                 {"variant": v, "window": list(w), "n_tracks": int(len(m))})
    return 0


def cmd_score_variants(args) -> int:
    vs, bad = read_variants(args.variants or [], args.input)
    if not vs:
        raise SystemExit("no variants parsed" + (f": {bad}" if bad else ""))
    if len(vs) > args.max_variants:
        raise SystemExit(f"{len(vs)} variants > --max-variants {args.max_variants}; the "
                         "API suits thousands of predictions, not millions — split the list "
                         "or use atlas-variants for pre-computed scores")
    if len(args.scorers or []) > MAX_SCORERS:
        raise SystemExit(f"at most {MAX_SCORERS} scorers per request")
    L = SEQ_LENGTHS[args.length]
    wins = [window(v["chrom"], v["position"] - 1, v["position"], L) for v in vs]
    b = get_backend(args.offline)
    df = b.score_variants(wins, vs, args.scorers or [], args.organism)
    full = df
    df = filter_ontology(df, args.ontology or [], args.search)
    d = run_dir(args.label or f"score_{len(vs)}variants")
    write_tsv(df, d / "variant_scores.tsv")
    top = top_scores(df, 30)
    cols = [c for c in ("variant_id", "variant_scorer", "gene_name", "track_name",
                        "biosample_name", "raw_score", "quantile_score") if c in df.columns]
    per_var = (df.assign(_a=df["quantile_score"].abs() if "quantile_score" in df.columns
                         else df["raw_score"].abs())
               .groupby("variant_id")["_a"].max().rename("max_abs_score").reset_index()
               .sort_values("max_abs_score", ascending=False)) if len(df) else df
    write_tsv(per_var, d / "variant_summary.tsv")
    lines = [f"- Variants scored: {len(vs)}" + (f"; unparsed: {len(bad)}" if bad else ""),
             f"- Scorers: {', '.join(args.scorers) if args.scorers else 'recommended set'}",
             f"- Window: {args.length} centred on each variant",
             f"- Score rows: {len(df):,} (of {len(full):,} before the ontology/search filter)",
             "", "## Strongest effects (by |quantile score| where available)", "",
             *md_table(top[cols] if len(top) else top, 30), "",
             "## Per variant", "", *md_table(per_var, 50)]
    if bad:
        lines += ["", "## Not parsed", "", *[f"- `{x['raw']}`: {x['error']}" for x in bad]]
    write_report(d, "AlphaGenome variant scores", lines,
                 {"n_variants": len(vs), "unparsed": bad, "n_rows": int(len(df)),
                  "scorers": args.scorers or "recommended"})
    return 0


def cmd_score_interval(args) -> int:
    c, s, e, label = resolve_region(args)
    w = window(c, s, e, SEQ_LENGTHS[args.length])
    b = get_backend(args.offline)
    df = filter_ontology(b.score_interval(*w, args.organism), args.ontology or [], args.search)
    d = run_dir(args.label or f"interval_{label.split()[0]}")
    write_tsv(df, d / "interval_scores.tsv")
    lines = [f"- Region: {label}; window {w[0]}:{w[1]:,}-{w[2]:,} ({args.length})",
             f"- Score rows: {len(df):,}", "", "## Highest scores", "",
             *md_table(top_scores(df, 25))]
    write_report(d, f"AlphaGenome interval scores · {label}", lines,
                 {"region": [c, s, e], "window": list(w), "n_rows": int(len(df))})
    return 0


def cmd_ism(args) -> int:
    c, s, e = parse_region(args.ism_interval)
    if e - s > args.max_bp:
        raise SystemExit(f"ISM over {e - s} bp is {3 * (e - s)} variants; the cap is "
                         f"--max-bp {args.max_bp}")
    w = window(c, s, e, SEQ_LENGTHS[args.length])
    b = get_backend(args.offline)
    df = filter_ontology(b.ism(w, (c, s, e), args.scorer or [], args.organism),
                         args.ontology or [], args.search)
    d = run_dir(args.label or f"ism_{c}_{s}")
    write_tsv(df, d / "ism_scores.tsv")
    mat = None
    if len(df):
        import numpy as np
        import pandas as pd
        parts = df["variant_id"].astype(str).str.extract(r":(\d+):([ACGTN]+)>([ACGTN]+)")
        df = df.assign(position=parts[0].astype(int), alt=parts[2])
        mat = (df.groupby(["position", "alt"])["raw_score"].mean().unstack("alt")
               .reindex(columns=["A", "C", "G", "T"]))
        mat.index.name = "position"
        mat.reset_index().to_csv(d / "ism_matrix.tsv", sep="\t", index=False)
        print(f"TSV: {d / 'ism_matrix.tsv'}")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(mat)), 2.4))
            vmax = float(np.nanmax(np.abs(mat.values))) or 1.0
            im = ax.imshow(mat.T.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax.set_yticks(range(4))
            ax.set_yticklabels(["A", "C", "G", "T"])
            ax.set_xticks(range(0, len(mat), max(1, len(mat) // 12)))
            ax.set_xticklabels(mat.index[::max(1, len(mat) // 12)], fontsize=6, rotation=45)
            fig.colorbar(im, ax=ax, label="mean raw score")
            ax.set_title(f"ISM {c}:{s:,}-{e:,}", fontsize=8)
            fig.tight_layout()
            fig.savefig(d / "ism_heatmap.png", dpi=140)
            plt.close(fig)
            print(f"Figure: {d / 'ism_heatmap.png'}")
        except Exception:
            pass
    sens = (mat.abs().max(axis=1).sort_values(ascending=False).rename("max_abs_effect")
            .reset_index().head(15)) if mat is not None else None
    lines = [f"- ISM interval: {c}:{s:,}-{e:,} ({e - s} bp, {3 * (e - s)} variants)",
             f"- Scorer(s): {', '.join(args.scorer) if args.scorer else 'recommended set'}; "
             f"ontology: {', '.join(args.ontology or ['all'])}", "",
             "## Most sensitive positions (mean over matched tracks)", "", *md_table(sens)]
    write_report(d, f"AlphaGenome in silico mutagenesis · {c}:{s}-{e}", lines,
                 {"ism_interval": [c, s, e], "n_rows": int(len(df))})
    return 0


def cmd_atlas_scorers(args) -> int:
    b = get_backend(args.offline)
    df = b.atlas_scorers()
    d = run_dir(args.label or "atlas_scorers")
    write_tsv(df, d / "atlas_scorers.tsv")
    write_report(d, "AlphaGenome Atlas scorers", [f"- Scorers: {len(df)}", "", *md_table(df, 60)],
                 {"n_scorers": int(len(df))})
    return 0


def _atlas_scorers(b, requested):
    if requested:
        return list(requested)
    return list(b.atlas_scorers()["scorer"])


def cmd_atlas_variants(args) -> int:
    vs, bad = read_variants(args.variants or [], args.input)
    if not vs:
        raise SystemExit("no variants parsed" + (f": {bad}" if bad else ""))
    b = get_backend(args.offline)
    df = b.atlas_variants(vs, _atlas_scorers(b, args.scorers), args.ontology or [], args.genes or [])
    d = run_dir(args.label or f"atlas_{len(vs)}variants")
    write_tsv(df, d / "atlas_scores.tsv")
    lines = [f"- Variants: {len(vs)}" + (f"; unparsed: {len(bad)}" if bad else ""),
             f"- Rows: {len(df):,}", "", "## Strongest pre-computed effects", "",
             *md_table(top_scores(df, 30))]
    write_report(d, "AlphaGenome Atlas variant scores", lines,
                 {"n_variants": len(vs), "n_rows": int(len(df)), "unparsed": bad})
    return 0


def cmd_atlas_interval(args) -> int:
    c, s, e, label = resolve_region(args)
    if e - s > args.max_bp:
        raise SystemExit(f"{e - s:,} bp > --max-bp {args.max_bp:,} (every position has 3 "
                         "variants; narrow the region)")
    b = get_backend(args.offline)
    df = b.atlas_interval(c, s, e, _atlas_scorers(b, args.scorers), args.ontology or [],
                          args.genes or [])
    d = run_dir(args.label or f"atlas_{label.split()[0]}")
    write_tsv(df, d / "atlas_scores.tsv")
    lines = [f"- Region: {label} ({e - s:,} bp)", f"- Rows: {len(df):,}", "",
             "## Strongest pre-computed effects", "", *md_table(top_scores(df, 30))]
    write_report(d, f"AlphaGenome Atlas · {label}", lines,
                 {"region": [c, s, e], "n_rows": int(len(df))})
    return 0


# ─── offline fake backend + selftest ────────────────────────────────────────

class FakeBackend:
    """Deterministic stand-in with the SdkBackend's return shapes."""

    TRACKS = [("UBERON:0000948", "heart", "total RNA-seq"),
              ("UBERON:0002107", "liver", "total RNA-seq"),
              ("CL:0000746", "cardiac muscle cell", "DNase-seq")]

    def _meta(self, ot):
        import pandas as pd
        return pd.DataFrame([{"name": f"{ot}_{i}", "ontology_curie": o, "biosample_name": b,
                              "Assay title": a, "strand": "."}
                             for i, (o, b, a) in enumerate(self.TRACKS)])

    def metadata(self, organism):
        import pandas as pd
        return pd.concat([self._meta(ot).assign(output_type=ot) for ot in ("RNA_SEQ", "DNASE")],
                         ignore_index=True)

    def _tracks(self, s, e, outputs, bump=None):
        import numpy as np
        L = e - s
        out = {}
        for ot in outputs:
            v = np.tile(np.linspace(0.1, 1.0, len(self.TRACKS)), (L, 1))
            if bump is not None:
                v[max(0, bump - s - 50):bump - s + 50, 0] += 5.0     # heart track only
            out[ot] = {"values": v, "metadata": self._meta(ot), "resolution": 1}
        return out

    def predict_interval(self, c, s, e, outputs, ontology, organism):
        return self._tracks(s, e, outputs)

    def predict_variant(self, c, s, e, v, outputs, ontology, organism):
        return self._tracks(s, e, outputs), self._tracks(s, e, outputs, bump=v["position"] - 1)

    def score_variants(self, windows, variants, scorers, organism):
        import pandas as pd
        rows = []
        for i, v in enumerate(variants):
            for o, b, a in self.TRACKS:
                raw = (2.0 if b == "heart" else 0.1) * (i + 1)
                rows.append({"variant_id": vid(v), "variant_scorer": "RNA_SEQ", "gene_name": "APOE",
                             "track_name": f"{b}_{a}", "ontology_curie": o, "biosample_name": b,
                             "raw_score": raw, "quantile_score": min(0.999, raw / 5)})
        return pd.DataFrame(rows)

    def score_interval(self, c, s, e, organism):
        import pandas as pd
        return pd.DataFrame([{"gene_name": "APOE", "track_name": b, "ontology_curie": o,
                              "biosample_name": b, "raw_score": 1.0 + k}
                             for k, (o, b, a) in enumerate(self.TRACKS)])

    def ism(self, win, iv, scorers, organism):
        import pandas as pd
        c, s, e = iv
        rows = []
        for p in range(s + 1, e + 1):
            for alt in "ACGT":
                rows.append({"variant_id": f"{c}:{p}:N>{alt}", "ontology_curie": "UBERON:0000948",
                             "biosample_name": "heart",
                             "raw_score": 3.0 if (p == s + 5 and alt == "T") else 0.01})
        return pd.DataFrame(rows)

    def atlas_scorers(self):
        import pandas as pd
        return pd.DataFrame([{"scorer": "AVI", "is_signed": False, "n_tracks": 1},
                             {"scorer": "RNA_SEQ", "is_signed": True, "n_tracks": 3}])

    def atlas_variants(self, variants, scorers, ontology, genes):
        import pandas as pd
        return pd.DataFrame([{"variant": vid(v), "scorer": sc, "raw_score": 0.5,
                              "quantile_score": 0.9} for v in variants for sc in scorers])

    def atlas_interval(self, c, s, e, scorers, ontology, genes):
        import pandas as pd
        return pd.DataFrame([{"variant": f"{c}:{p}:A>G", "scorer": sc, "raw_score": 0.1}
                             for p in range(s + 1, min(e, s + 5) + 1) for sc in scorers])


def cmd_selftest(args) -> int:
    fails = []

    def check(name, cond):
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
        if not cond:
            fails.append(name)

    global OUT_ROOT, CACHE_DIR, KEY_FILE
    os.environ["IGVF_ALPHAGENOME_QUIET"] = "1"     # no log files from the selftest
    tmp = Path(tempfile.mkdtemp(prefix="alphagenome_selftest_"))
    OUT_ROOT, CACHE_DIR, KEY_FILE = tmp / "out", tmp / "cache", tmp / "Secret" / "KEY.txt"
    fake_rsid = lambda t: {"raw": t, "chrom": "chr19", "position": 44908684, "ref": "T",
                           "alt": "C", "rsid": t}

    print("notation")
    v = parse_variant("chr22:36201698:A>C")
    check("chr:pos:ref>alt keeps the 1-based position",
          (v["chrom"], v["position"], v["ref"], v["alt"]) == ("chr22", 36201698, "A", "C"))
    check("2-21001846-G-A gains the chr prefix", parse_variant("2-21001846-G-A")["chrom"] == "chr2")
    sp = parse_variant("NC_000019.10:44908683:T:C")
    check("SPDI's 0-based position becomes 1-based (rs429358 = chr19:44908684)",
          (sp["chrom"], sp["position"]) == ("chr19", 44908684))
    check("rsIDs resolve through the Catalog resolver",
          parse_variant("rs429358", fake_rsid)["position"] == 44908684)
    try:
        parse_variant("not-a-variant")
        check("bad notation raises", False)
    except ValueError:
        check("bad notation raises", True)
    vcf = tmp / "v.vcf"
    vcf.write_text("#CHROM POS ID REF ALT\nchr1\t100\t.\tA\tG,T\n# c\nrs1 chr2:5:C>G junk\n")
    got, bad = read_variants([], str(vcf), fake_rsid)
    check("VCF multi-allelic rows split, notations mixed per line, dupes dropped",
          len(got) == 4 and [b["raw"] for b in bad] == ["junk"])
    check("region parses 0-based half-open", parse_region("chr1:1,000-2,000") == ("chr1", 1000, 2000))
    c, s, e = window("chr1", 1000, 2000, SEQ_LENGTHS["16KB"])
    check("window is centred and exactly the supported length",
          e - s == 16384 and s <= 1000 and e >= 2000)
    try:
        window("chr1", 0, 2_000_000, SEQ_LENGTHS["1MB"])
        check("regions over 1 Mb are refused", False)
    except ValueError:
        check("regions over 1 Mb are refused", True)

    print("\ncredentials")
    saved = {k: os.environ.pop(k, None) for k in KEY_ENV}
    check("no key -> None", api_key() is None)
    KEY_FILE.parent.mkdir(parents=True)
    KEY_FILE.write_text("  secret-value-123 \n")
    check("key file is read and stripped", api_key() == "secret-value-123")
    os.environ["ALPHAGENOME_API_KEY"] = "env-wins"
    check("environment beats the file", api_key() == "env-wins")
    buf_ok = "secret-value-123" not in key_help() and "env-wins" not in key_help()
    check("help text never contains the key", buf_ok)
    for k, val in saved.items():
        if val is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = val

    print("\ncommands (fake backend, no network)")
    import contextlib
    import io

    def run(argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(argv + ["--offline"])
        return rc, out.getvalue()

    def arte(out, kind):
        return [l.split(": ", 1)[1] for l in out.splitlines() if l.startswith(kind + ": ")]

    rc, out = run(["metadata", "--search", "heart"])
    tracks = arte(out, "TSV")
    import pandas as pd
    t = pd.read_csv(tracks[0], sep="\t") if tracks else pd.DataFrame()
    check("metadata --search heart keeps only heart tracks",
          rc == 0 and len(t) == 2 and set(t["biosample_name"]) == {"heart"})
    rc, out = run(["predict-interval", "--interval", "chr19:44905790-44909393",
                   "--outputs", "RNA_SEQ", "--ontology", "UBERON:0000948", "--length", "16KB"])
    check("predict-interval writes a summary, a figure and a report",
          rc == 0 and arte(out, "TSV") and arte(out, "Report") and
          (arte(out, "Figure") or args.no_plots))
    try:
        run(["predict-interval", "--interval", "chr1:1-100", "--outputs", "RNA_SEQ"])
        check("predict-interval without --ontology is refused", False)
    except SystemExit:
        check("predict-interval without --ontology is refused", True)
    rc, out = run(["predict-variant", "--variant", "chr19:44908684:T>C", "--outputs", "RNA_SEQ",
                   "--ontology", "UBERON:0000948", "--length", "16KB"])
    ra = pd.read_csv(arte(out, "TSV")[0], sep="\t")
    check("predict-variant ranks the planted ALT gain in heart first",
          rc == 0 and ra.iloc[0]["biosample_name"] == "heart" and ra.iloc[0]["log2fc_sum"] > 1)
    rc, out = run(["score-variants", "--variants", "chr19:44908684:T>C", "chr22:36201698:A>C",
                   "--search", "heart"])
    sv = pd.read_csv(arte(out, "TSV")[0], sep="\t")
    check("score-variants filters to heart and reports both variants",
          rc == 0 and set(sv["biosample_name"]) == {"heart"} and sv["variant_id"].nunique() == 2)
    try:
        run(["score-variants", "--variants", *[f"chr1:{i}:A>G" for i in range(1, 6)],
             "--max-variants", "3"])
        check("score-variants enforces --max-variants", False)
    except SystemExit:
        check("score-variants enforces --max-variants", True)
    rc, out = run(["score-interval", "--interval", "chr19:44905790-44909393", "--length", "16KB"])
    check("score-interval runs", rc == 0 and arte(out, "TSV"))
    rc, out = run(["ism", "--ism-interval", "chr19:44908670-44908690", "--length", "16KB"])
    mat = [p for p in arte(out, "TSV") if p.endswith("ism_matrix.tsv")]
    m = pd.read_csv(mat[0], sep="\t") if mat else pd.DataFrame()
    check("ism builds a position x base matrix with the planted hotspot on top",
          rc == 0 and len(m) == 20 and m.set_index("position")["T"].idxmax() == 44908675)
    try:
        run(["ism", "--ism-interval", "chr1:0-500"])
        check("ism caps the interval width", False)
    except SystemExit:
        check("ism caps the interval width", True)
    rc, out = run(["atlas-scorers"])
    check("atlas-scorers lists scorers", rc == 0 and arte(out, "TSV"))
    rc, out = run(["atlas-variants", "--variants", "chr19:44908684:T>C"])
    av = pd.read_csv(arte(out, "TSV")[0], sep="\t")
    check("atlas-variants defaults to every Atlas scorer", rc == 0 and set(av["scorer"]) == {"AVI", "RNA_SEQ"})
    rc, out = run(["atlas-interval", "--interval", "chr19:44908680-44908690", "--scorers", "AVI"])
    check("atlas-interval runs", rc == 0 and arte(out, "TSV"))
    rep = Path(arte(out, "Report")[0]).read_text() if arte(out, "Report") else ""
    check("every report carries the AlphaGenome terms", "non-commercial" in rep)

    if sdk_available():
        print("\nSDK objects (installed, no network)")
        from alphagenome.data import genome  # type: ignore
        iv = SdkBackend._interval("chr19", 44900000, 44900000 + SEQ_LENGTHS["16KB"])
        var = SdkBackend._variant(parse_variant("chr19:44908684:T>C"))
        check("SDK Interval/Variant build from our records",
              isinstance(iv, genome.Interval) and var.position == 44908684)
        import numpy as np
        from alphagenome.data import track_data  # type: ignore
        from alphagenome.models import dna_output  # type: ignore
        meta = pd.DataFrame({"name": ["heart_rna", "liver_rna"], "strand": [".", "."],
                             "ontology_curie": ["UBERON:0000948", "UBERON:0002107"],
                             "biosample_name": ["heart", "liver"]})
        td = track_data.TrackData(values=np.ones((128, 2), dtype=np.float32) * [[1.0, 3.0]],
                                  metadata=meta, resolution=1, interval=iv.resize(128))
        tr = SdkBackend._tracks(dna_output.Output(rna_seq=td))
        summ = summarise_tracks(tr, iv.resize(128).start, iv.resize(128).start,
                                iv.resize(128).end)
        check("real SDK TrackData converts and summarises (per-track mean)",
              list(tr) == ["RNA_SEQ"] and
              summ.set_index("biosample_name")["mean"].round(3).to_dict() == {"heart": 1.0, "liver": 3.0})
    else:
        print("\n  skip  SDK object checks (SDK not importable here)")

    shutil.rmtree(tmp, ignore_errors=True)
    print("\nselftest: all checks pass" if not fails else f"\nselftest: FAILED ({len(fails)})")
    return 0 if not fails else 1


# ─── CLI ────────────────────────────────────────────────────────────────────

API_COMMANDS = {"metadata", "predict-interval", "predict-variant", "score-variants",
                "score-interval", "ism", "atlas-scorers", "atlas-variants", "atlas-interval"}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="igvfagent alphagenome",
                                 description="AlphaGenome predictions and Atlas scores "
                                             "(google-deepmind/alphagenome).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, org=True):
        if org:
            p.add_argument("--organism", choices=ORGANISMS, default="HOMO_SAPIENS")
        p.add_argument("--label")
        p.add_argument("--offline", action="store_true", help=argparse.SUPPRESS)
        return p

    def onto(p):
        p.add_argument("--ontology", nargs="+", action="extend", default=[],
                       help="Ontology CURIE(s), e.g. UBERON:0000948 (heart), CL:0000746. "
                            "`alphagenome metadata --search <tissue>` lists them.")
        return p

    def length(p, default="1MB"):
        p.add_argument("--length", choices=list(SEQ_LENGTHS), default=default,
                       help="Prediction window (default %(default)s).")
        return p

    def region(p):
        p.add_argument("--interval", help="chr:start-end, 0-based half-open")
        p.add_argument("--gene", help="gene symbol, resolved through the IGVF Catalog")
        return p

    s = sub.add_parser("setup", help="Check SDK, interpreter and key.")
    s.add_argument("--install", action="store_true",
                   help="Create a Python 3.10+ environment with the SDK if needed.")
    s.add_argument("--ping", action="store_true", help="Call the API once (metadata).")
    s.set_defaults(func=cmd_setup)

    p = common(onto(sub.add_parser("metadata", help="Every predicted track.")))
    p.add_argument("--search", help="Substring of biosample, CURIE or track name.")
    p.add_argument("--output-type", nargs="+", action="extend", choices=OUTPUT_TYPES)
    p.set_defaults(func=cmd_metadata)

    p = common(length(onto(region(sub.add_parser("predict-interval", help="Predicted tracks over a region.")))))
    p.add_argument("--outputs", nargs="+", action="extend", choices=OUTPUT_TYPES, required=True)
    p.add_argument("--all-tracks", action="store_true")
    p.set_defaults(func=cmd_predict_interval)

    p = common(length(onto(sub.add_parser("predict-variant", help="REF vs ALT tracks around a variant."))))
    p.add_argument("--variant", required=True)
    p.add_argument("--outputs", nargs="+", action="extend", choices=OUTPUT_TYPES, required=True)
    p.add_argument("--radius", type=int, default=1000,
                   help="Summarise within +/- this many bp of the variant.")
    p.add_argument("--all-tracks", action="store_true")
    p.set_defaults(func=cmd_predict_variant)

    p = common(length(onto(sub.add_parser("score-variants", help="Variant effect scores."))))
    p.add_argument("--variants", nargs="+", action="extend", default=[])
    p.add_argument("--input", help="File of variants (any notation, or VCF).")
    p.add_argument("--scorers", nargs="+", action="extend", default=[],
                   help="Recommended scorer names (DNASE, RNA_SEQ, SPLICE_SITES, ...). "
                        "Default: the full recommended set.")
    p.add_argument("--search", help="Keep rows whose biosample/track/tissue contains this.")
    p.add_argument("--max-variants", type=int, default=MAX_VARIANTS)
    p.set_defaults(func=cmd_score_variants)

    p = common(length(onto(region(sub.add_parser("score-interval", help="Interval (gene-level) scores.")))))
    p.add_argument("--search")
    p.set_defaults(func=cmd_score_interval)

    p = common(length(onto(sub.add_parser("ism", help="In silico mutagenesis."))))
    p.add_argument("--ism-interval", required=True, help="chr:start-end, 0-based half-open")
    p.add_argument("--scorer", nargs="+", action="extend", default=[], help="Scorer name(s); default recommended set.")
    p.add_argument("--search")
    p.add_argument("--max-bp", type=int, default=MAX_ISM_BP)
    p.set_defaults(func=cmd_ism)

    p = common(sub.add_parser("atlas-scorers", help="Atlas scorers available."), org=False)
    p.set_defaults(func=cmd_atlas_scorers)

    for name, fn, helptext in (("atlas-variants", cmd_atlas_variants, "Atlas scores for variants."),
                               ("atlas-interval", cmd_atlas_interval, "Atlas scores in a region.")):
        p = common(onto(sub.add_parser(name, help=helptext)), org=False)
        if name == "atlas-variants":
            p.add_argument("--variants", nargs="+", action="extend", default=[])
            p.add_argument("--input")
        else:
            region(p)
            p.add_argument("--max-bp", type=int, default=MAX_ATLAS_BP)
            p.set_defaults(organism="HOMO_SAPIENS")
        p.add_argument("--scorers", nargs="+", action="extend", default=[], help="Default: every Atlas scorer.")
        p.add_argument("--genes", nargs="+", action="extend", default=[], help="Restrict gene-level scores.")
        p.set_defaults(func=fn)

    s = sub.add_parser("selftest", help="Offline checks (no key, no network).")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_selftest)
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    if (args.cmd in API_COMMANDS and not getattr(args, "offline", False)
            and os.environ.get("IGVF_ALPHAGENOME_FAKE") != "1" and not sdk_available()):
        if os.environ.get(REEXEC_FLAG):
            print("AlphaGenome SDK not importable in the SDK interpreter either; "
                  "re-run `igvfagent alphagenome setup --install`.")
            return 2
        return reexec(argv)
    if args.cmd != "selftest" and not os.environ.get("IGVF_ALPHAGENOME_QUIET"):
        setup_logging()
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
