# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/hbliu/Kidney_Epi_Pri @ 552c461f85b9
# (eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh) for liu2025_kidney_multiome. Unreviewed; provenance in Scripts/ported/registry.json.
"""Independent GWAS loci from summary statistics (Liu et al., Science 2025).

Python port of hbliu/Kidney_Epi_Pri eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh,
the Susztak-lab locus definition the Science 2025 eGFRcrea GWAS reuses:

  1. keep variants with P < 5e-8                          (awk, Step1)
  2. drop the MHC, chr6 25-35 Mb (strict bounds)          (awk, Step2)
  3. plink --clump --clump-p1 5e-8 --clump-r2 0.1
     --clump-kb 10000 against 1000 Genomes EUR, SNP ids CHR:POS  (Step3)
  4. order clumped leads by CHR, POS; start a new locus when the chromosome
     changes or the genetic distance to the previous lead is >= 0.1 cM;
     keep the lowest-P lead per locus                    (R loop, Step4)

Step 3 shells out to plink 1.9 exactly as the authors did. Step 4 reads
genetic positions by linear interpolation in a HapMap-II GRCh37 map, which
is what `plink --cm-map` does for the authors' 1000GP_Phase3 combined_b37
maps (the same HapMap-II map). Steps 1-2 alone run without plink
(--stop-after-filter).
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

P_SIG = 5e-8
MHC_CHR, MHC_LO, MHC_HI = 6, 25_000_000, 35_000_000


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def filter_significant(gwas: str, chr_col: str, pos_col: str, p_col: str,
                       id_col: str) -> tuple[pd.DataFrame, dict]:
    """Steps 1-2. Returns the non-MHC significant variants and counts."""
    n_total = 0
    keep = []
    usecols = [id_col, chr_col, pos_col, p_col]
    for chunk in pd.read_csv(gwas, sep="\t", usecols=usecols, chunksize=2_000_000,
                             dtype={chr_col: str}, low_memory=False):
        n_total += len(chunk)
        # awk compares the field numerically; "5e-300"-style strings and
        # values below double range ("1e-400") parse to floats / 0 here.
        p = pd.to_numeric(chunk[p_col], errors="coerce")
        keep.append(chunk[p < P_SIG].assign(**{p_col: p[p < P_SIG]}))
    sig = pd.concat(keep, ignore_index=True)
    chrom = pd.to_numeric(sig[chr_col], errors="coerce")
    in_mhc = (chrom == MHC_CHR) & (sig[pos_col] > MHC_LO) & (sig[pos_col] < MHC_HI)
    counts = {"n_variants_tested": int(n_total),
              "n_significant_all": int(len(sig)),
              "n_significant_mhc": int(in_mhc.sum()),
              "n_significant_non_mhc": int((~in_mhc).sum())}
    out = sig[~in_mhc].copy()
    out["CHR"] = pd.to_numeric(out[chr_col], errors="coerce").astype("Int64")
    out["POS"] = out[pos_col].astype(int)
    out["SNP"] = out["CHR"].astype(str) + ":" + out["POS"].astype(str)
    out["P"] = out[p_col]
    out["RSID"] = out[id_col]
    return out[["SNP", "CHR", "POS", "P", "RSID"]], counts


def run_clump(sig: pd.DataFrame, bfile: str, plink: str, workdir: Path) -> pd.DataFrame:
    """Step 3: plink 1.9 clumping, flags as in the authors' script."""
    workdir.mkdir(parents=True, exist_ok=True)
    assoc = workdir / "forPlinkClumping.txt"
    sig[["SNP", "P"]].to_csv(assoc, sep="\t", index=False)
    out = workdir / "PlinkClumping.p5e8.r20.1"
    cmd = [plink, "--bfile", bfile, "--clump", str(assoc),
           "--clump-p1", "5e-8", "--clump-r2", "0.1", "--clump-kb", "10000",
           "--out", str(out)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    clumped = pd.read_csv(f"{out}.clumped", sep=r"\s+")
    return clumped


def load_map(map_path: str) -> dict:
    """HapMap-II map: chr position rate cM (Eagle layout) or per-chr files."""
    m = pd.read_csv(map_path, sep=r"\s+", compression="infer")
    m.columns = [c.lower() for c in m.columns]
    cm_col = [c for c in m.columns if "cm" in c and "rate" not in c][-1]
    chr_col = [c for c in m.columns if c.startswith("chr")][0]
    pos_col = [c for c in m.columns if "position" in c or c == "pos"][0]
    maps = {}
    for c, g in m.groupby(chr_col):
        g = g.sort_values(pos_col)
        maps[str(c).replace("chr", "")] = (g[pos_col].to_numpy(), g[cm_col].to_numpy())
    return maps


def cm_of(maps: dict, chrom, pos) -> float:
    pos_arr, cm_arr = maps[str(chrom)]
    return float(np.interp(pos, pos_arr, cm_arr))


def merge_leads(clumped: pd.DataFrame, maps: dict, max_cm: float = 0.1) -> pd.DataFrame:
    """Step 4: the authors' R loop merging consecutive leads within 0.1 cM."""
    leads = clumped[["CHR", "BP", "SNP", "P"]].rename(columns={"BP": "POS"}).copy()
    leads["cM"] = [cm_of(maps, c, p) for c, p in zip(leads["CHR"], leads["POS"])]
    leads = leads.sort_values(["CHR", "POS"], kind="mergesort").reset_index(drop=True)
    merged = np.empty(len(leads), dtype=int)
    locus = 1
    merged[0] = locus
    for i in range(1, len(leads)):
        if (leads.at[i, "CHR"] == leads.at[i - 1, "CHR"]
                and abs(leads.at[i, "cM"] - leads.at[i - 1, "cM"]) < max_cm):
            merged[i] = locus
        else:
            locus += 1
            merged[i] = locus
    leads["Merged"] = merged
    return leads


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gwas", required=True, help="GWAS summary statistics (tsv/tsv.gz)")
    ap.add_argument("--chr-col", default="CHR")
    ap.add_argument("--pos-col", default="POS")
    ap.add_argument("--p-col", default="P.value")
    ap.add_argument("--id-col", default="MarkerName")
    ap.add_argument("--bfile", help="plink bfile prefix of the LD reference (SNP ids CHR:POS)")
    ap.add_argument("--plink", default=shutil.which("plink") or "plink")
    ap.add_argument("--genetic-map", help="HapMap-II GRCh37 genetic map")
    ap.add_argument("--stop-after-filter", action="store_true",
                    help="run only steps 1-2 (no plink / map needed)")
    ap.add_argument("--out", required=True, help="output directory")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sig, counts = filter_significant(a.gwas, a.chr_col, a.pos_col, a.p_col, a.id_col)
    sig.to_csv(out / "significant_non_mhc.tsv.gz", sep="\t", index=False)
    summary = {"method": "port of hbliu/Kidney_Epi_Pri eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh",
               "input": str(a.gwas), **counts}
    if not a.stop_after_filter:
        if not (a.bfile and a.genetic_map):
            ap.error("--bfile and --genetic-map are required unless --stop-after-filter")
        clumped = run_clump(sig, a.bfile, a.plink, out / "plink")
        maps = load_map(a.genetic_map)
        leads = merge_leads(clumped, maps)
        leads.to_csv(out / "clumped_leads_merged.tsv", sep="\t", index=False)
        top = (leads.sort_values("P", kind="mergesort").groupby("Merged", sort=True)
               .head(1).sort_values(["CHR", "POS"]))
        top.to_csv(out / "independent_loci.tsv", sep="\t", index=False)
        summary.update({
            "n_significant_in_reference": int(sig["SNP"].isin(
                pd.read_csv(f"{a.bfile}.bim", sep=r"\s+", header=None, usecols=[1])[1]).sum()),
            "n_clumped_leads": int(len(leads)),
            "n_independent_loci": int(leads["Merged"].nunique()),
        })
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
