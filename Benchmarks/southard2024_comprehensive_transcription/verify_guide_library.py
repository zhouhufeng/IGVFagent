#!/usr/bin/env python3
"""Verify paper-claimed guide-library counts against the authors' own
deposited data (Zenodo 10.5281/zenodo.15200179), rather than trusting the
paper's prose.

Downloads fibroblast_CRISPRa_aggr_total_guide_umis.h5 (58MB; a small
slice of the ~8-46GB processed/raw deposits under the same Zenodo
community) and reads its guide-identity axis directly.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA_DIR = ROOT / "Data" / "Benchmarks" / "southard2024_comprehensive_transcription"
H5_PATH = DATA_DIR / "fibroblast_CRISPRa_aggr_total_guide_umis.h5"
H5_URL = ("https://zenodo.org/api/records/15200179/files/"
          "fibroblast_CRISPRa_aggr_total_guide_umis.h5/content")
OUT_PATH = DATA_DIR / "guide_library_verification.json"


def main() -> int:
    import h5py

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not H5_PATH.exists():
        print(f"Downloading {H5_URL} -> {H5_PATH}", file=sys.stderr)
        urllib.request.urlretrieve(H5_URL, H5_PATH)

    with h5py.File(H5_PATH, "r") as f:
        guides = [g.decode() for g in f["guide_umis/index_level1"][:]]
        n_droplets = int(f["guide_umis/index_level0"].shape[0])

    n_guides = len(guides)
    n_non_targeting = sum(1 for g in guides if g.startswith("non_targeting_"))
    n_tf_targets = len({g.split("_")[0] for g in guides
                         if not g.startswith("non_targeting_")})

    out = {
        "source": {
            "zenodo_doi": "10.5281/zenodo.15200179",
            "zenodo_record": "Hs27-CRISPRa-TFs",
            "file": "fibroblast_CRISPRa_aggr_total_guide_umis.h5",
        },
        "n_guides_total": n_guides,
        "n_non_targeting_guides": n_non_targeting,
        "n_distinct_tf_targets": n_tf_targets,
        "n_droplets_pre_filter": n_droplets,
        "paper_claims": {
            "guides (Methods, final library)": 10979,
            "non-targeting negative controls (Methods)": 78,
            "transcription factors targeted (Summary)": 1836,
        },
    }
    OUT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
