#!/usr/bin/env python3
"""What analysis does each IGVF assay actually call for.

The failure this exists to end: routing was decided from FILE types, so any
dataset with FASTQ reads and no published matrix was sent to transcriptome
quantification. That is right for an RNA readout and wrong for everything
else, and the Portal has 65 distinct `preferred_assay_titles` across 11,070
MeasurementSets. Patching them one at a time as they failed -- SGE, then
gRNA-sequencing screens -- meant the default stayed wrong for the rest.

So the default is inverted here. An assay is classified explicitly or it is
reported as UNRECOGNISED; nothing falls through to "quantify it as RNA".
Being told "I do not know how to analyse Pooled Y2H" is a usable answer.
Being handed gene counts for it is not.

`READOUT` is consulted before the title, because for a CRISPR screen the
title says how cells were selected ("CRISPR FACS screen") while
`crispr_screen_readout` says what was sequenced -- and only the second
determines the analysis. The same FACS screen is a transcript library when
its readout is scRNA-seq and a guide library when it is gRNA sequencing.

Counts in the comments are datasets on the Portal as of 2026-09, kept so a
future reader can see which gaps actually matter.
"""

from __future__ import annotations

from typing import Optional

# Route names. `transcript` is the only one that may go to kallisto.
TRANSCRIPT = "transcript"      # RNA readout: quantify against a transcriptome
GUIDE = "guide"                # reads are an sgRNA library
VARIANT = "variant"            # amplicon of a variant library (SGE/MAVE)
ELEMENT = "element"            # MPRA/STARR: element activity from barcodes
CHROMATIN = "chromatin"        # ATAC/ChIP/Hi-C: genomic, not transcript
GENOME = "genome"              # WGS/methylation
PROTEIN = "protein"            # Y2H/IPA interaction assays
UNKNOWN = "unknown"

# How to actually analyse each route, and what IGVFagent can do today.
ROUTE_GUIDANCE = {
    TRANSCRIPT: ("transcript quantification then single-cell or bulk analysis",
                 "supported: raw_pipeline_run"),
    GUIDE: ("per-guide counts, and enrichment between selected populations",
            "supported: raw_pipeline_guide_count"),
    VARIANT: ("per-variant functional scores from amplicon variant calling",
              "supported for SGE: sge_analyze"),
    ELEMENT: ("per-element activity from barcode counts",
              "partially supported: the mpra_* and starr_* tools"),
    CHROMATIN: ("peak or contact analysis against the genome, not a "
                "transcriptome",
                "NOT supported end-to-end; no aligner for genomic reads"),
    GENOME: ("variant calling or methylation calling against the genome",
             "NOT supported; no genomic aligner"),
    PROTEIN: ("interaction scoring, not sequence quantification",
              "NOT supported from raw reads"),
}

# crispr_screen_readout -> route. Checked BEFORE the title.
READOUT = {
    "scrna-seq": TRANSCRIPT,                        # 1573 datasets
    "scrna-seq including guide capture": TRANSCRIPT,       # 2
    "snatac-seq": CHROMATIN,                              # 1
    "grna sequencing": GUIDE,                       # 560
    "sgrna sequencing": GUIDE,
    "guide sequencing": GUIDE,
    "endogenous allelic sequencing": VARIANT,       # 535
    "exogenous allelic sequencing": VARIANT,        # 55
}

# preferred_assay_titles / assay_titles -> route.
ASSAY = {
    # --- RNA readout -----------------------------------------------------
    "10x multiome": TRANSCRIPT,                     # 1982
    "snmct-seq": TRANSCRIPT,                        # 933
    "cc-perturb-seq": TRANSCRIPT,                   # 818
    "share-seq": TRANSCRIPT,                        # 682
    "perturb-seq": TRANSCRIPT,                      # 507
    "parse split-seq": TRANSCRIPT,                  # 354
    "sge-rna": TRANSCRIPT,                          # 318
    "tap-seq": TRANSCRIPT,                          # 200
    "rna-seq": TRANSCRIPT,                          # 196
    "10x multiome with multi-seq": TRANSCRIPT,      # 142
    "snm3c-seq": TRANSCRIPT,                        # 90
    "scnt-seq2": TRANSCRIPT,                        # 79
    "mtscmultiome": TRANSCRIPT,                     # 76
    "morf-share-seq": TRANSCRIPT,                   # 73
    "scnt-seq3": TRANSCRIPT,                        # 25
    "parse perturb-seq": TRANSCRIPT,                # 19
    "spatial transcriptomics": TRANSCRIPT,          # 18
    "scrna-seq": TRANSCRIPT,                        # 16
    "scmultiome-nt-seq": TRANSCRIPT,                # 12
    "ont drna": TRANSCRIPT,                         # 8
    "in vivo perturb-seq": TRANSCRIPT,              # 4
    "scnt-seq": TRANSCRIPT,                         # 4
    "crop-seq": TRANSCRIPT,                         # 2
    "multiome perturb-seq": TRANSCRIPT,             # 2
    # --- variant libraries ----------------------------------------------
    "sge": VARIANT,                                 # 1167
    "immune-sge": VARIANT,                          # 841
    "variant-effects": VARIANT,                     # 488
    "label-seq": VARIANT,                           # 196
    "vamp-seq (multistep)": VARIANT,                # 144
    "vamp-seq": VARIANT,                            # 56
    "mave": VARIANT,                                # 21
    "varaccess": VARIANT,                           # 10
    "variant painting via fluorescence": VARIANT,   # 1
    "variant painting via immunostaining": VARIANT, # 1
    # --- reporter assays -------------------------------------------------
    "mpra": ELEMENT,                                # 114
    "starr-seq": ELEMENT,                           # 52
    "in vivo mpra": ELEMENT,                        # 28
    "electroporated mpra": ELEMENT,                 # 27
    "lentimpra": ELEMENT,                           # 27
    "aav-mpra": ELEMENT,                            # 24
    "mpra (scqer)": ELEMENT,                        # 6
    # --- CRISPR screens: title says selection, readout says what was
    #     sequenced. Defaulted to GUIDE only when the readout is absent,
    #     because a screen with no stated readout is far more often a
    #     guide library than a transcript one.
    "crispr flowfish screen": GUIDE,                # 366
    "crispr facs screen": GUIDE,                    # 300
    "proliferation crispr screen": GUIDE,           # 79
    "sccrispr screen": TRANSCRIPT,                  # 24 (single-cell readout)
    "morf screen": GUIDE,                           # 12
    "crispr macs screen": GUIDE,                    # 9
    "migration crispr screen": GUIDE,               # 6
    # --- chromatin / genome ----------------------------------------------
    "10x snatac-seq with scale pre-indexing": CHROMATIN,   # 144
    "atac-seq": CHROMATIN,                          # 115
    "hicar": CHROMATIN,                             # 20
    "cut&run": CHROMATIN,                           # 16
    "tf chip-seq": CHROMATIN,                       # 16
    "hi-c": CHROMATIN,                              # 12
    "histone chip-seq": CHROMATIN,                  # 8
    "ont fiber-seq": CHROMATIN,                     # 8
    "snatac-seq": CHROMATIN,                        # 4
    "wgs": GENOME,                                  # 75
    "bisulfite-seq": GENOME,                        # 26
    "ont direct wgs": GENOME,                       # 8
    # --- protein interaction ---------------------------------------------
    "pooled y2h": PROTEIN,                          # 46
    "arrayed semi-qy2h v1": PROTEIN,
    "arrayed semi-qy2h v2": PROTEIN,
    "arrayed semi-qy2h v3": PROTEIN,
    "dual-ipa": PROTEIN,                            # 1
}


def classify(file_set: dict) -> dict:
    """Route for a FileSet, with what it was decided from.

    Returns {route, matched_on, assay, analysis, support}. ``route`` is
    UNKNOWN when nothing matches -- deliberately, so the caller can say so
    instead of quantifying an unrecognised assay as RNA.
    """
    readout = str((file_set or {}).get("crispr_screen_readout") or "").strip()
    if readout and readout.lower() in READOUT:
        route = READOUT[readout.lower()]
        return _pack(route, f"crispr_screen_readout={readout}", readout)

    titles = []
    for key in ("preferred_assay_titles", "assay_titles"):
        v = (file_set or {}).get(key)
        if isinstance(v, list):
            titles.extend((key, str(x)) for x in v)
        elif v:
            titles.append((key, str(v)))
    for key, raw in titles:
        route = ASSAY.get(raw.strip().lower())
        if route:
            return _pack(route, f"{key}={raw}", raw)
    seen = ", ".join(raw for _, raw in titles) or "none stated"
    return _pack(UNKNOWN, f"no rule for: {seen}", seen)


def _pack(route: str, matched_on: str, assay: str) -> dict:
    analysis, support = ROUTE_GUIDANCE.get(
        route, ("unknown -- this assay has no classification yet",
                "NOT supported; the assay is unrecognised"))
    return {"route": route, "matched_on": matched_on, "assay": assay,
            "analysis": analysis, "support": support}


def coverage() -> dict:
    """How many assay titles and readouts are classified."""
    return {"assay_titles_classified": len(ASSAY),
            "readouts_classified": len(READOUT),
            "routes": sorted(set(ASSAY.values()) | set(READOUT.values()))}


__all__ = ["classify", "coverage", "ASSAY", "READOUT", "ROUTE_GUIDANCE",
           "TRANSCRIPT", "GUIDE", "VARIANT", "ELEMENT", "CHROMATIN",
           "GENOME", "PROTEIN", "UNKNOWN"]
