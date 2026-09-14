#!/usr/bin/env python3
"""RegulomeDB: regulatory rank and per-tissue scores for non-coding variants.

Built as the reachable part of TLand-predict (rnsherpa/TLand-predict, MIT),
which annotates non-coding variants across organs. TLand's full pipeline needs
Snakemake, conda, Sei's GPU model and TLand's own trained weights -- none of
which are here, and the weights are the scientific asset rather than the code.
RegulomeDB is its other input, it is a REST API, and it already carries the
organ-specific signal: `tissue_specific_scores` gives a probability per tissue
alongside the genome-wide rank.

WHAT IT ADDS. IGVFagent could already say whether a variant is conserved and
deleterious (FAVOR: CADD, GERP), what the Catalog knows about it, and what
ClinVar calls it. None of those say whether it sits in something REGULATORY.
RegulomeDB's rank (1a best .. 7 worst) is built from ChIP-seq, chromatin
accessibility, footprints, motifs and QTL evidence, so it answers a different
question: not "is this base important" but "is this base in a regulatory
element, in which tissues".

WHY THE RESPONSE IS NOT SIMPLY RETURNED. One variant returns ~1.3 MB, almost
all of it the `@graph` of every supporting ENCODE experiment. Handing that to
a model would flood the context for one variant and crowd out the answer, so
the score, the evidence flags and the tissue scores are extracted, and the
supporting experiments are COUNTED and summarised by assay rather than
reproduced.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint          # noqa: E402,F401
import _localstore as ls                                     # noqa: E402

ROOT = ls.ROOT
OUT_DIR = ROOT / "Docs" / "RegulomeDB"
BASE = "https://www.regulomedb.org/regulome-search/"
USER_AGENT = "IGVFagent/regulome (+https://github.com/zhouhufeng/IGVFagent)"

# Rank 1a-2c carry direct experimental evidence; 3a-6 are progressively
# weaker; 7 is "no data". Stated so a caller reading a bare "4" knows where
# that sits rather than guessing.
RANK_MEANING = {
    "1a": "eQTL + TF binding + matched motif + matched footprint + open chromatin",
    "1b": "eQTL + TF binding + any motif + footprint + open chromatin",
    "1c": "eQTL + TF binding + matched motif + open chromatin",
    "1d": "eQTL + TF binding + any motif + open chromatin",
    "1f": "eQTL + TF binding / open chromatin",
    "2a": "TF binding + matched motif + matched footprint + open chromatin",
    "2b": "TF binding + any motif + footprint + open chromatin",
    "2c": "TF binding + matched motif + open chromatin",
    "3a": "TF binding + any motif + open chromatin",
    "4":  "TF binding + open chromatin",
    "5":  "TF binding OR open chromatin",
    "6":  "motif hit only",
    "7":  "no regulatory evidence",
}


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])


def _get(region: str, assembly: str, timeout: int) -> dict:
    url = (BASE + "?" + urllib.parse.urlencode(
        {"regions": region, "genome": assembly, "format": "json"}))
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                 "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def parse_variant(tok: str) -> "str | None":
    """'chr1:39492461', 'chr1-39492461-A-G' or a 1-bp region -> a region."""
    t = tok.strip()
    if not t or t.startswith("#"):
        return None
    if ":" in t and "-" in t.split(":", 1)[1]:
        return t                                    # already a region
    for sep in (":", "-", "\t"):
        if sep in t:
            f = t.replace("\t", sep).split(sep)
            if len(f) >= 2 and f[1].isdigit():
                c = f[0] if f[0].startswith("chr") else "chr" + f[0]
                p = int(f[1])
                # RegulomeDB wants a half-open 1-bp window.
                return f"{c}:{p - 1}-{p}"
    return None


def summarise(d: dict) -> dict:
    score = d.get("regulome_score") or {}
    feats = d.get("features") or {}
    tissue = score.get("tissue_specific_scores") or {}
    graph = d.get("@graph") or []
    assays = Counter()
    biosamples = Counter()
    for g in graph:
        a = (g.get("assay_title") or g.get("assay_term_name") or "")
        if a:
            assays[a] += 1
        b = ((g.get("biosample_ontology") or {}).get("term_name") or "")
        if b:
            biosamples[b] += 1
    top_tissue = sorted(((k, float(v)) for k, v in tissue.items()
                         if str(v).replace(".", "", 1).isdigit()),
                        key=lambda kv: -kv[1])
    rank = str(score.get("ranking") or "")
    return {
        "rank": rank,
        "rank_meaning": RANK_MEANING.get(rank, ""),
        "probability": score.get("probability"),
        "chip_seq": feats.get("ChIP"),
        "chromatin_accessibility": feats.get("Chromatin_accessibility"),
        "footprint": feats.get("Footprint"),
        "motif": feats.get("PWM"),
        "qtl": feats.get("QTL"),
        "n_supporting_experiments": len(graph),
        "top_assays": "; ".join(f"{k}({v})" for k, v in assays.most_common(4)),
        "n_biosamples": len(biosamples),
        "top_tissues": "; ".join(f"{k}={v:.3f}" for k, v in top_tissue[:5]),
        "_tissue_scores": {k: v for k, v in top_tissue},
    }


def cmd_annotate(args) -> int:
    setup_logging()
    toks: "list[str]" = []
    if args.variants:
        toks += [t for t in args.variants.replace(",", " ").split() if t]
    if args.input:
        toks += [l.strip() for l in Path(args.input).read_text().splitlines()
                 if l.strip() and not l.startswith("#")]
    regions = []
    for t in toks:
        r = parse_variant(t)
        if r:
            regions.append((t, r))
        else:
            logging.warning("could not parse %r as a variant/region", t)
    if not regions:
        raise SystemExit("no parsable variants. Use chr1:39492461 or "
                          "chr1-39492461-A-G, or --input with one per line.")
    logging.info("%d variant(s)", len(regions))

    rows, tissue_rows = [], []
    for i, (tok, region) in enumerate(regions, 1):
        try:
            d = _get(region, args.assembly, args.timeout)
            s = summarise(d)
        except Exception as e:                              # noqa: BLE001
            logging.warning("%s failed: %s", tok, e)
            rows.append({"query": tok, "region": region, "rank": "",
                         "rank_meaning": f"lookup failed: {type(e).__name__}"})
            continue
        ts = s.pop("_tissue_scores")
        rows.append({"query": tok, "region": region, **s})
        for t, v in ts.items():
            tissue_rows.append({"query": tok, "tissue": t, "score": f"{v:.5f}"})
        if i % 5 == 0:
            logging.info("  %d/%d", i, len(regions))
        time.sleep(args.sleep)              # be a polite API client

    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_regulome"
    out = Path(args.out_dir) if args.out_dir else (OUT_DIR / label)
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "regulome_annotations.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in rows[0]})
    ttsv = out / "regulome_tissue_scores.tsv"
    if tissue_rows:
        with ttsv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t",
                               fieldnames=["query", "tissue", "score"])
            w.writeheader()
            for r in tissue_rows:
                w.writerow(r)
    js = out / "regulome_summary.json"
    js.write_text(json.dumps({"variants": len(rows),
                              "with_rank": sum(1 for r in rows if r.get("rank")),
                              "tissue_rows": len(tissue_rows),
                              "assembly": args.assembly}, indent=2))
    print(f"Report:  {tsv}")
    if tissue_rows:
        print(f"Wrote:   {ttsv}")
    print(f"Summary: {js}")
    for r in rows[:10]:
        print(f"  {r['query']:22} rank={r.get('rank',''):3} "
              f"p={r.get('probability','')}  {str(r.get('rank_meaning',''))[:52]}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvfagent regulome",
        description="RegulomeDB regulatory rank + per-tissue scores for "
                     "non-coding variants.")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("annotate", help="Variants -> rank, evidence, tissues.")
    a.add_argument("--variants", help="chr1:39492461 or chr1-39492461-A-G, "
                                       "comma/space separated.")
    a.add_argument("--input", help="File with one variant per line.")
    a.add_argument("--assembly", default="GRCh38")
    a.add_argument("--sleep", type=float, default=0.4)
    a.add_argument("--timeout", type=int, default=120)
    a.add_argument("--label")
    a.add_argument("--out-dir")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return {"annotate": cmd_annotate}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
