"""Verify the two ports absorbed from this paper against the authors' code.

  liu2025_gwas_loci  vs hbliu/Kidney_Epi_Pri eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh
                        (authors' awk Steps 1-2; authors' Step 3 formatting + Step 4 R loop)
  liu2025_open4gene  vs hbliu/Open4Gene R/Open4Gene.R run on its bundled test data
                        (a 62,278-cell slice of the paper's kidney multiome), and vs the
                        paper's published Open4Gene summary statistics for the same pairs

Each comparison goes through `igvfagent bench verify-port`; the resulting
validation_vs_reference.json is copied to Data/liu2025/verify/<tag>/ so that
expected.json can point at a stable path per check.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "Data/liu2025"
REF_O4G = D / "ref_authors/open4gene"
OUT_O4G = D / "out/open4gene"
IGVF = shutil.which("igvfagent") or "igvfagent"
PAPER = "liu2025_kidney_multiome"


def verify(tag: str, *args: str) -> str:
    cmd = [IGVF, "bench", "verify-port", "--paper-id", PAPER, *args]
    res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    text = res.stdout + res.stderr
    verdict = next((l for l in text.splitlines() if l.startswith(("PASS", "FAIL"))), text[-400:])
    wrote = re.search(r"^Wrote (\S+validation_vs_reference\.json)", text, re.M)
    if not wrote:
        sys.exit(f"verify-port {tag} produced no artefact:\n{text}")
    dest = D / "verify" / tag
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / wrote.group(1), dest / "validation_vs_reference.json")
    print(f"{tag:34s} {verdict}")
    time.sleep(1.1)          # verify-port names its output directory by the second
    return verdict


def gwas_tables() -> None:
    p = pd.read_csv(D / "out/gwas_loci_cli/significant_non_mhc.tsv.gz", sep="\t")
    p.assign(ID=p.RSID + "_" + p.SNP)[["ID", "P"]].to_csv(D / "out/gwas_loci_cli/port_ID_P.tsv",
                                                          sep="\t", index=False)


def open4gene_tables() -> None:
    R = pd.concat([pd.read_csv(REF_O4G / f"R_{c}.res.txt", sep="\t") for c in ("All", "Each")])
    P = pd.concat([pd.read_csv(OUT_O4G / f"port_{c}.res.txt", sep="\t") for c in ("All", "Each")])
    for t in (R, P):
        t["ID"] = t.Peak + "|" + t.Gene + "|" + t.Celltype
    R.to_csv(REF_O4G / "R_all.tsv", sep="\t", index=False)
    P.to_csv(OUT_O4G / "port_all.tsv", sep="\t", index=False)
    # count component is only identifiable where R's own SE is finite and < 10
    ident = set(R.ID[R["hurdle.Res.count.se"] < 10])
    R[R.ID.isin(ident)].to_csv(REF_O4G / "R_identifiable_count.tsv", sep="\t", index=False)
    P[P.ID.isin(ident)].to_csv(OUT_O4G / "port_identifiable_count.tsv", sep="\t", index=False)
    # the paper's own published Open4Gene statistics for the same peak-gene pairs
    cols = ["gene", "peak", "Celltype", "TotalCellNum", "ExpressCellNum", "OpenCellNum",
            "hurdle.Res.zero.beta", "zse", "zz", "hurdle.Res.zero.p", "hurdle.Res.count.beta", "cse", "cz",
            "cp", "aic", "bic", "spearman.rho", "sp"]
    rows = []
    for chunk in pd.read_csv(D / "figshare/Open4Gene_all.txt.gz", sep="\t", header=0, names=cols,
                             chunksize=2_000_000, na_values="NA"):
        chunk["Celltype"] = chunk.Celltype.replace({"AllCell": "All"})
        chunk["ID"] = chunk.peak + "|" + chunk.gene + "|" + chunk.Celltype
        rows.append(chunk[chunk.ID.isin(set(P.ID))])
    pd.concat(rows).to_csv(OUT_O4G / "published_paper_testpairs.tsv", sep="\t", index=False)


def main() -> int:
    gwas_tables()
    open4gene_tables()
    g, a = D / "ref_authors", D / "out/gwas_loci_cli"
    verify("gwas_filter_vs_awk", "--name", "liu2025_gwas_loci", "--analysis", "fig1a_gwas_significant",
           "--reference", str(g / "ref_ID_P.tsv"), "--port-output", str(a / "port_ID_P.tsv"),
           "--column", "P", "--id-column", "ID", "--rtol", "1e-6", "--atol", "1e-300",
           "--min-match-rate", "0.999")
    verify("loci_merged_vs_R", "--name", "liu2025_gwas_loci", "--analysis", "fig1b_independent_loci",
           "--reference", str(g / "step34/clumped_leads_merged.tsv"),
           "--port-output", str(a / "clumped_leads_merged.tsv"),
           "--column", "Merged", "--id-column", "SNP", "--rtol", "0", "--atol", "0", "--min-match-rate", "1.0")
    verify("loci_cM_vs_plink", "--name", "liu2025_gwas_loci", "--analysis", "fig1b_independent_loci",
           "--reference", str(g / "step34/clumped_leads_merged.tsv"),
           "--port-output", str(a / "clumped_leads_merged.tsv"),
           "--column", "cM", "--id-column", "SNP", "--rtol", "1e-4", "--atol", "1e-4", "--min-match-rate", "0.99")
    o = ["--name", "liu2025_open4gene", "--analysis", "fig5_open4gene_method", "--id-column", "ID"]
    verify("o4g_zero_beta_vs_R", *o, "--reference", str(REF_O4G / "R_all.tsv"),
           "--port-output", str(OUT_O4G / "port_all.tsv"), "--column", "hurdle.Res.zero.beta",
           "--rtol", "1e-6", "--atol", "1e-6", "--min-match-rate", "1.0")
    verify("o4g_zero_p_vs_R", *o, "--reference", str(REF_O4G / "R_all.tsv"),
           "--port-output", str(OUT_O4G / "port_all.tsv"), "--column", "hurdle.Res.zero.p",
           "--rtol", "1e-4", "--atol", "1e-9", "--min-match-rate", "1.0")
    verify("o4g_count_z_identifiable_vs_R", *o, "--reference", str(REF_O4G / "R_identifiable_count.tsv"),
           "--port-output", str(OUT_O4G / "port_identifiable_count.tsv"), "--column", "hurdle.Res.count.z",
           "--rtol", "0", "--atol", "0.01", "--min-match-rate", "0.95")
    verify("o4g_count_beta_all_vs_R", *o, "--reference", str(REF_O4G / "R_all.tsv"),
           "--port-output", str(OUT_O4G / "port_all.tsv"), "--column", "hurdle.Res.count.beta",
           "--rtol", "1e-4", "--atol", "1e-5", "--min-match-rate", "0.0")
    verify("o4g_zero_beta_vs_published", *o, "--reference", str(OUT_O4G / "published_paper_testpairs.tsv"),
           "--port-output", str(OUT_O4G / "port_all.tsv"), "--column", "hurdle.Res.zero.beta",
           "--rtol", "1e-6", "--atol", "1e-6", "--min-match-rate", "1.0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
