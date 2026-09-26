"""Parse the paper's own raw count deposit (GEO GSE271788) into a clean
genes x samples count matrix plus a sample metadata table.

Input:  Data/Weinstock2024/raw/GSE271788_dedup_counts.txt
        (featureCounts output the authors deposited verbatim on GEO;
         311 samples = 84 KO'd genes x 3 donors + 59 AAVS1 control samples)
Output: Data/Weinstock2024/processed/counts.csv.gz   (genes x samples, raw)
        Data/Weinstock2024/processed/sample_meta.csv (sample, donor, ko, is_control)
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "Data/Weinstock2024/raw/GSE271788_dedup_counts.txt"
OUT_DIR = ROOT / "Data/Weinstock2024/processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COL_RE = re.compile(r"output/bam/dedup/Donor_(\d+)_(.+)\.dedup\.bam")


def main() -> None:
    df = pd.read_csv(RAW, sep="\t", skiprows=1)
    meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    sample_cols = [c for c in df.columns if c not in meta_cols]

    records = []
    rename = {}
    for c in sample_cols:
        m = COL_RE.match(c)
        if not m:
            raise ValueError(f"unrecognised sample column: {c}")
        donor, ko = m.group(1), m.group(2)
        is_control = ko.startswith("AAVS1")
        sample_id = f"Donor{donor}_{ko}"
        rename[c] = sample_id
        records.append(
            {
                "sample": sample_id,
                "donor": donor,
                "ko": "AAVS1" if is_control else ko,
                "is_control": is_control,
            }
        )

    meta = pd.DataFrame.from_records(records)
    assert meta["sample"].is_unique, "duplicate sample ids after renaming"

    counts = df[["Geneid"] + sample_cols].rename(columns=rename)
    counts["Geneid"] = counts["Geneid"].str.replace(r"\.\d+$", "", regex=True)
    counts = counts.set_index("Geneid")

    print(f"genes: {counts.shape[0]}  samples: {counts.shape[1]}")
    print(meta["is_control"].value_counts())
    print(meta["donor"].value_counts().sort_index())

    counts.to_csv(OUT_DIR / "counts.csv.gz")
    meta.to_csv(OUT_DIR / "sample_meta.csv", index=False)
    print(f"wrote {OUT_DIR / 'counts.csv.gz'} and {OUT_DIR / 'sample_meta.csv'}")


if __name__ == "__main__":
    main()
