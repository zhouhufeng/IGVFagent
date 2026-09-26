"""Recompute the paper's headline numbers from its own deposited tables.

Liu et al., Science 2025 (doi:10.1126/science.adp4753) deposit their derived
results on Figshare 26299093 (CC BY). This script recomputes, from those
files, the quantities the paper states for Figs. 2B, 3A, 3D, 4B-C, 4F, 5B-D
and 6B-E, and writes them to Data/liu2025/verify/derived/derived_metrics.json
for concordance.py. Where the paper does not say how a quantity was computed,
the choice made here is recorded next to the number ("definition").

Inputs (downloaded by run.sh): Data/liu2025/figshare/*.
"""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
FS = ROOT / "Data/liu2025/figshare"
OUT = ROOT / "Data/liu2025/verify/derived"
SCORECARD_TSV = ROOT / "Data/liu2025/out/scorecard.tsv.gz"


def scorecard() -> pd.DataFrame:
    if not SCORECARD_TSV.is_file():
        import openpyxl
        SCORECARD_TSV.parent.mkdir(parents=True, exist_ok=True)
        ws = openpyxl.load_workbook(FS / "Scorecard.xlsx", read_only=True).worksheets[0]
        with gzip.open(SCORECARD_TSV, "wt", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            for row in ws.iter_rows(values_only=True):
                w.writerow(["" if v is None else v for v in row[:46]])
    d = pd.read_csv(SCORECARD_TSV, sep="\t", low_memory=False)
    return d.loc[:, ~d.columns.str.startswith("Unnamed")]


def in_intervals(chroms, pos, iv: pd.DataFrame, half_open: bool) -> np.ndarray:
    """pos inside any interval of iv (columns c, s, e) on the same chromosome."""
    hit = np.zeros(len(pos), bool)
    chroms = np.asarray(chroms)
    pos = np.asarray(pos)
    for c in np.unique(chroms):
        P = iv[iv.c == c].sort_values("s")
        if P.empty:
            continue
        st, en = P.s.to_numpy(), P.e.to_numpy()
        m = chroms == c
        p = pos[m]
        i = np.searchsorted(st, p - (1 if half_open else 0), side="right") - 1
        h = np.zeros(len(p), bool)
        for k in range(8):                                   # overlapping intervals
            j = np.clip(i - k, 0, None)
            lo_ok = (st[j] < p) if half_open else (st[j] <= p)
            h |= (i - k >= 0) & lo_ok & (p <= en[j])
        hit[m] = h
    return hit


def main() -> int:
    M: dict = {}

    # ---- Fig. 2B: ASE genes (RASQUAL, FDR 1%); the paper counts protein-coding genes
    ase = {}
    for tag in ("Tubule", "Glomeruli"):
        d = pd.read_csv(FS / f"ASE_{tag}.txt.gz", sep="\t", usecols=["GeneID", "GeneType"])
        ase[tag] = set(d.GeneID[d.GeneType == "protein_coding"])
    union = ase["Tubule"] | ase["Glomeruli"]
    M["fig2b"] = {
        "definition": "distinct protein_coding GeneID in the deposited significant-association tables",
        "n_ase_genes_tubule": len(ase["Tubule"]), "n_ase_genes_glomeruli": len(ase["Glomeruli"]),
        "n_ase_genes_union": len(union),
        "ase_shared_fraction_of_union": len(ase["Tubule"] & ase["Glomeruli"]) / len(union),
        "paper_implied_shared_fraction_of_union": (8573 + 8725 - 10398) / 10398,
    }

    # ---- Fig. 3A / S5G: bulk allele-specific accessibility
    b = pd.read_csv(FS / "bASA.txt.gz", sep="\t",
                    usecols=["Feature_ID", "rs_ID", "Effect_size", "Log10_Benjamini_Hochberg_Qvalue",
                             "Distance", "MAF"])
    loc = b.assign(local=b.Distance == 0).groupby("Feature_ID").local.agg(["any", "all"])
    n_ld = int((loc["any"] & ~loc["all"]).sum())
    M["fig3a"] = {
        "n_basa_peaks": int(b.Feature_ID.nunique()), "n_basa_snps": int(b.rs_ID.nunique()),
        "min_maf": float(b.MAF.min()),
        "n_basa_peaks_local_and_distal": n_ld,
        "fraction_basa_peaks_local_and_distal": n_ld / b.Feature_ID.nunique(),
        "paper_implied_fraction_local_and_distal": 8317 / 19083,
    }

    # ---- Fig. 3D: bASA vs ASE effect sizes (join unspecified in the paper)
    bb = b.sort_values("Log10_Benjamini_Hochberg_Qvalue").drop_duplicates("rs_ID")
    a = pd.read_csv(FS / "ASE_Tubule.txt.gz", sep="\t",
                    usecols=["RSID", "Effect_size", "Log10_Benjamini_Hochberg_Qvalue"])
    aa = a.sort_values("Log10_Benjamini_Hochberg_Qvalue").drop_duplicates("RSID")
    mm = bb.merge(aa, left_on="rs_ID", right_on="RSID", suffixes=("_b", "_a"))
    M["fig3d"] = {
        "definition": "SNPs significant in both bASA and tubule ASE; strongest association per SNP; "
                      "Spearman of RASQUAL effect sizes (pi)",
        "n_shared_snps": int(len(mm)),
        "basa_ase_spearman_correlation": float(spearmanr(mm.Effect_size_b, mm.Effect_size_a)[0]),
    }

    # ---- Fig. 4B-D: single-nucleus ASA
    s = pd.read_csv(FS / "snASA.txt.gz", sep="\t")
    n_ct = s.groupby("ID").Celltype.nunique()
    M["fig4b"] = {
        "definition": "SNP identity = rsID",
        "n_snasa_snps": int(len(n_ct)), "n_celltype_specific": int((n_ct == 1).sum()),
        "n_shared": int((n_ct > 1).sum()),
        "snasa_celltype_specific_fraction": float((n_ct == 1).mean()),
    }
    peaks = pd.read_csv(FS / "snATAC_peaks.bed.gz", sep="\t", header=None, names=["c", "s", "e"])
    u = s.drop_duplicates("ID")
    inpk = in_intervals("chr" + u.CHR.astype(str), u.POS.to_numpy(), peaks, half_open=True)
    M["fig4c"] = {
        "definition": "snASA SNP inside the deposited Human.Kidney.OpenChromatin.snATAC peak set (BED)",
        "n_in_peaks": int(inpk.sum()), "snasa_in_peak_fraction": float(inpk.mean()),
    }
    # ---- Fig. 4F: PT snASA vs bASA
    pt = s[s.Celltype == "PT"].assign(alt_fraction=lambda x: x["ALT.Total"] / x["REF.ALT.Total"])
    bl = b[b.Distance == 0].sort_values("Log10_Benjamini_Hochberg_Qvalue").drop_duplicates("rs_ID")
    m4 = pt.merge(bl, left_on="ID", right_on="rs_ID")
    ball = pt.merge(bb, left_on="ID", right_on="rs_ID")
    M["fig4f"] = {
        "definition": "PT snASA ALT read fraction vs bASA effect size (pi) for the SNP inside its own "
                      "bASA peak (local pairs); strongest peak per SNP",
        "n_shared_snps": int(len(m4)),
        "snasa_pt_vs_basa_spearman_correlation": float(spearmanr(m4.alt_fraction, m4.Effect_size)[0]),
        "alt_definition_all_pairs_spearman_correlation": float(spearmanr(ball.alt_fraction, ball.Effect_size)[0]),
    }

    # ---- Fig. 5B-C / S8: Open4Gene links (FDR < 0.01, deposited)
    o = pd.read_csv(FS / "Open4Gene_sig.txt.gz", sep="\t", low_memory=False)
    o["pair"] = o.peak + "|" + o.gene
    allc, ct = o[o.Celltype == "AllCell"], o[o.Celltype != "AllCell"]
    up = o.drop_duplicates("pair")
    peaks_per_gene = up.groupby("gene").peak.nunique()
    genes_per_peak = up.groupby("peak").gene.nunique()
    zsig = o["zinb.res.zero.p.All.fdr"] < 0.01
    M["fig5c"] = {
        "n_links": int(len(up)), "n_peaks": int(o.peak.nunique()), "n_genes": int(o.gene.nunique()),
        "n_celltype_links": int(ct.pair.nunique()),
        "fraction_celltype_links_not_in_allcell": float((~ct.pair.drop_duplicates().isin(allc.pair)).mean()),
        "fraction_peaks_single_gene": float((genes_per_peak == 1).mean()),
        "median_peaks_per_gene": float(peaks_per_gene.median()),
        "top_gene_by_peaks": str(peaks_per_gene.idxmax()), "n_peaks_top_gene": int(peaks_per_gene.max()),
        "positive_beta_fraction": float((o.loc[zsig, "zinb.res.zero.beta"] > 0).mean()),
        "positive_beta_definition": "zero-component beta > 0 among rows with zero-component FDR < 0.01",
    }

    # ---- Fig. 5D / S9A: GWAS variants (P < 5e-8, non-MHC, GRCh37) in Open4Gene peaks
    gw = pd.read_csv(FS / "eGFRcrea_GWAS_Multi.txt.gz", sep="\t", usecols=["MarkerName", "CHR", "POS", "P.value"])
    gw = gw[(gw["P.value"] < 5e-8) & ~((gw.CHR == 6) & (gw.POS > 25e6) & (gw.POS < 35e6))]
    pk = up[["peak", "gene"]].copy()
    pk[["c", "s", "e"]] = pk.peak.str.split("-", expand=True)
    pk["s"], pk["e"] = pk.s.astype(int), pk.e.astype(int)
    iv = pk.drop_duplicates("peak")
    rows = []
    for c, g in gw.groupby("CHR"):
        P = iv[iv.c == f"chr{c}"].sort_values("s")
        st, en, nm = P.s.to_numpy(), P.e.to_numpy(), P.peak.to_numpy()
        pos = g.POS.to_numpy()
        i = np.searchsorted(st, pos, side="right") - 1
        for k in range(8):
            j = np.clip(i - k, 0, None)
            ok = (i - k >= 0) & (st[j] <= pos) & (pos <= en[j])
            rows.append(pd.DataFrame({"MarkerName": g.MarkerName.to_numpy()[ok], "peak": nm[j[ok]]}))
    hits = pd.concat(rows).drop_duplicates().merge(pk[["peak", "gene"]], on="peak")
    per_var = hits.groupby("MarkerName").gene.nunique()
    M["fig5d"] = {
        "definition": "multi-ancestry GWAS P < 5e-8 outside chr6:25-35 Mb, inside an Open4Gene peak (1-based inclusive)",
        "n_gwas_variants_in_open4gene_peaks": int(len(per_var)), "n_target_genes": int(hits.gene.nunique()),
        "gwas_variant_single_gene_fraction": float((per_var == 1).mean()),
    }

    # ---- Fig. 6: Kidney Disease Genetic Scorecard
    d = scorecard()
    v, gcol, lcol = "Variant RSID", "Target gene symbol", "Location relative to gene"
    reg = d[(d["Prioritization Score"] >= 10) & (d[lcol] != "CDS")]
    top = d.sort_values("Prioritization Score", ascending=False, kind="mergesort").iloc[0]
    M["fig6b"] = {
        "definition": "regulatory = variant-gene rows with Prioritization Score >= 10 whose variant is not in the target CDS",
        "n_regulatory_variants": int(reg[v].nunique()), "n_regulatory_target_genes": int(reg[gcol].nunique()),
        "top_variant": str(top[v]), "top_variant_gene": str(top[gcol]),
        "top_prioritization_score": int(top["Prioritization Score"]),
    }
    cds = d[d[lcol] == "CDS"]
    per_gene = cds.groupby(gcol)[v].nunique()
    M["fig6d"] = {
        "definition": "coding variants = scorecard rows located in the target gene's CDS (the scorecard covers "
                      "957 GWAS loci; the paper's 1,363-variant table S24 is not deposited)",
        "n_coding_variants": int(cds[v].nunique()), "n_coding_genes": int(len(per_gene)),
        "fraction_genes_multiple_coding_variants": float((per_gene >= 2).mean()),
        "top_coding_gene": str(per_gene.idxmax()), "n_coding_variants_top_gene": int(per_gene.max()),
    }
    both = set(reg[gcol]) & set(d.loc[d["Gene CDS"] == 1, gcol])
    val = set(reg[reg[gcol].isin(both) & ((reg.eGFRcys == 1) | (reg.BUN == 1))][gcol]) | \
        set(d[d[gcol].isin(both) & ((d["Gene CDS eGFRcys"] == 1) | (d["Gene CDS BUN"] == 1))][gcol])
    M["fig6e"] = {
        "definition": "genes with a regulatory variant (score >= 10) and the Gene CDS flag; validated = a "
                      "regulatory or coding variant also associated with eGFRcys or BUN",
        "n_genes_coding_and_regulatory": len(both), "n_validated": len(val),
        "validated_fraction": len(val) / len(both),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "derived_metrics.json").write_text(json.dumps(M, indent=2))
    print(json.dumps(M, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
