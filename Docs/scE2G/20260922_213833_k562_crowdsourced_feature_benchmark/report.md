# E2G predictor benchmark on CRISPR element-gene pairs: k562_crowdsourced_feature_benchmark

Generated 2026-09-22 21:38 by `igvfagent sce2g merge` from 6 benchmark run(s): 20260922_212801_k562_features_A, 20260922_212801_k562_features_B, 20260922_212801_k562_features_C, 20260922_212801_k562_features_D_epcot, 20260922_212801_k562_re2g_published_anchor, 20260922_212801_k562_sce2g_reference.

CRISPR pairs: 10,356, of which 471 regulated (random-baseline precision 0.0455). Each predictor's score for a CRISPR pair is the aggregate of the predicted elements overlapping the tested element for the same gene (max unless stated), `fill_value` when none overlaps; p-value / FDR / distance style columns are inverted so that higher means more likely. `auprc_negated` is the AUPRC of the negated score, so a feature that works in the other direction is visible.

## Ranked predictors

| rank | predictor | AUPRC | 95% CI | precision @ 70% recall | AUPRC negated | x random | pairs overlapped |
|---|---|---|---|---|---|---|---|
| 1 | rE2G_ext | 0.7581 | 0.725-0.796 | 0.699 | 0.0234 | 16.7 | 1.000 |
| 2 | rE2G_ext_noEP300 | 0.7303 | 0.692-0.770 | 0.672 | 0.0234 | 16.1 | 1.000 |
| 3 | rE2G_base | 0.6341 | 0.592-0.689 | 0.543 | 0.0235 | 13.9 | 1.000 |
| 4 | scE2G_ignoreTPM | 0.5891 | 0.542-0.635 | 0.473 | 0.0275 | 12.9 | 0.891 |
| 5 | scE2G | 0.5304 | 0.486-0.577 | 0.254 | 0.0321 | 11.7 | 0.891 |
| 6 | ARC_E2G | 0.4952 | 0.454-0.548 | 0.336 | 0.0276 | 10.9 | 0.891 |
| 7 | ABC | 0.4910 | 0.450-0.545 | 0.352 | 0.0276 | 10.8 | 0.891 |
| 8 | baseline:E2G_Distance | 0.4181 | 0.376-0.470 | 0.273 | 0.0278 | 9.2 | 0.861 |
| 9 | Pinloop:Pinloop | 0.3569 | 0.308-0.404 | 0.197 | 0.0274 | 7.8 | 0.861 |
| 10 | Signac:Signac_Score | 0.2267 | 0.187-0.267 | 0.081 | 0.0286 | 5.0 | 0.799 |
| 11 | Signac:Signac_pvalue | 0.2121 | 0.171-0.252 | 0.069 | 0.0312 | 4.7 | 0.799 |
| 12 | SCENT:SCENT_beta | 0.1746 | 0.148-0.205 | 0.047 | 0.0366 | 3.8 | 0.449 |
| 13 | EPCOT:EP300 signal at E | 0.1585 | 0.135-0.190 | 0.096 | 0.0289 | 3.5 | 0.861 |
| 14 | EPCOT_woK562:EP300 signal at E | 0.1383 | 0.118-0.166 | 0.086 | 0.0289 | 3.0 | 0.861 |
| 15 | EPCOT_woK562:H3K4me1 signal at E | 0.1095 | 0.093-0.133 | 0.087 | 0.0299 | 2.4 | 0.861 |
| 16 | EPCOT_woK562:H3K27ac signal at E | 0.1084 | 0.093-0.130 | 0.087 | 0.0288 | 2.4 | 0.861 |
| 17 | EPCOT:H3K4me1 signal at E | 0.0933 | 0.080-0.109 | 0.083 | 0.0304 | 2.0 | 0.861 |
| 18 | EPCOT:H3K27ac signal at E | 0.0910 | 0.080-0.108 | 0.084 | 0.0314 | 2.0 | 0.861 |
| 19 | EPCOT:H3K4me1 signal at P | 0.0795 | 0.064-0.101 | 0.049 | 0.0370 | 1.8 | 0.861 |
| 20 | Cicero:Cicero_Score | 0.0767 | 0.066-0.090 | 0.046 | 0.0358 | 1.7 | 0.529 |
| 21 | EPCOT_woK562:H3K4me1 signal at P | 0.0738 | 0.061-0.090 | 0.056 | 0.0337 | 1.6 | 0.861 |
| 22 | EPCOT_woK562:GRO-seq reverse signal at E | 0.0727 | 0.066-0.083 | 0.076 | 0.0316 | 1.6 | 0.861 |
| 23 | EPCOT_woK562:GRO-seq forward signal at E | 0.0699 | 0.063-0.080 | 0.072 | 0.0327 | 1.5 | 0.861 |
| 24 | EPCOT_woK562:H3K9ac signal at E | 0.0685 | 0.061-0.078 | 0.070 | 0.0373 | 1.5 | 0.861 |
| 25 | ENCODEstats:disp | 0.0683 | 0.054-0.085 | 0.043 | 0.0463 | 1.5 | 0.861 |
| 26 | EPCOT_woK562:GRO-cap reverse signal at E | 0.0674 | 0.061-0.076 | 0.075 | 0.0314 | 1.5 | 0.861 |
| 27 | EPCOT_woK562:H3K4me2 signal at E | 0.0674 | 0.060-0.076 | 0.073 | 0.0364 | 1.5 | 0.861 |
| 28 | EPCOT_woK562:Bru-seq signal at E | 0.0658 | 0.058-0.076 | 0.066 | 0.0342 | 1.4 | 0.861 |
| 29 | EPCOT_woK562:NET-CAGE reverse signal at E | 0.0657 | 0.059-0.074 | 0.071 | 0.0333 | 1.4 | 0.861 |
| 30 | EPCOT:GRO-cap reverse signal at E | 0.0653 | 0.059-0.073 | 0.072 | 0.0342 | 1.4 | 0.861 |
| 31 | EPCOT:H3K4me2 signal at E | 0.0648 | 0.058-0.073 | 0.072 | 0.0325 | 1.4 | 0.861 |
| 32 | EPCOT:GRO-cap forward signal at E | 0.0638 | 0.058-0.072 | 0.070 | 0.0342 | 1.4 | 0.861 |
| 33 | EPCOT_woK562:NET-CAGE forward signal at E | 0.0638 | 0.057-0.072 | 0.071 | 0.0339 | 1.4 | 0.861 |
| 34 | EPCOT:GRO-seq reverse signal at E | 0.0633 | 0.058-0.072 | 0.067 | 0.0342 | 1.4 | 0.861 |
| 35 | EPCOT:H2AFZ signal at E | 0.0630 | 0.055-0.073 | 0.065 | 0.0356 | 1.4 | 0.861 |
| 36 | EPCOT_woK562:GRO-cap forward signal at E | 0.0622 | 0.056-0.071 | 0.066 | 0.0349 | 1.4 | 0.861 |
| 37 | EPCOT:H3K27me3 signal at P | 0.0611 | 0.050-0.075 | 0.048 | 0.0429 | 1.3 | 0.861 |
| 38 | EPCOT:NET-CAGE reverse signal at E | 0.0607 | 0.054-0.069 | 0.061 | 0.0361 | 1.3 | 0.861 |
| 39 | EPCOT_woK562:H3K79me2 signal at E | 0.0604 | 0.054-0.069 | 0.059 | 0.0381 | 1.3 | 0.861 |
| 40 | EPCOT:H3K9me3 signal at E | 0.0601 | 0.051-0.073 | 0.050 | 0.0462 | 1.3 | 0.861 |
| 41 | SCENT:SCENT_Score | 0.0599 | 0.053-0.068 | 0.082 | 0.0339 | 1.3 | 0.449 |
| 42 | EPCOT_woK562:H2AFZ signal at E | 0.0595 | 0.053-0.067 | 0.060 | 0.0376 | 1.3 | 0.861 |
| 43 | ArchR:ArchR_FDR | 0.0582 | 0.049-0.073 | 0.045 | 0.0414 | 1.3 | 0.602 |
| 44 | EPCOT:GRO-seq forward signal at E | 0.0577 | 0.053-0.065 | 0.062 | 0.0353 | 1.3 | 0.861 |
| 45 | EPCOT_woK562:H3K36me3 signal at P | 0.0577 | 0.045-0.075 | 0.047 | 0.0430 | 1.3 | 0.861 |
| 46 | EPCOT:H3K36me3 signal at P | 0.0571 | 0.048-0.072 | 0.049 | 0.0419 | 1.3 | 0.861 |
| 47 | EPCOT:H3K9me3 signal at P | 0.0568 | 0.050-0.066 | 0.055 | 0.0366 | 1.2 | 0.861 |
| 48 | ChromHMM:CrhmmBool-Enh_Strong | 0.0567 | 0.051-0.063 | 0.045 | 0.0401 | 1.2 | 0.861 |
| 49 | EPCOT_woK562:H3K4me2 signal at P | 0.0554 | 0.048-0.064 | 0.046 | 0.0443 | 1.2 | 0.861 |
| 50 | ArchR:ArchR_Score | 0.0550 | 0.046-0.067 | 0.045 | 0.0483 | 1.2 | 0.602 |
| 51 | EPCOT:H3K4me3 signal at E | 0.0547 | 0.049-0.062 | 0.055 | 0.0434 | 1.2 | 0.861 |
| 52 | EPCOT:H3K9ac signal at E | 0.0546 | 0.049-0.062 | 0.051 | 0.0427 | 1.2 | 0.861 |
| 53 | EPCOT:NET-CAGE forward signal at E | 0.0545 | 0.050-0.061 | 0.062 | 0.0363 | 1.2 | 0.861 |
| 54 | EPCOT:H3K79me2 signal at E | 0.0517 | 0.047-0.059 | 0.045 | 0.0475 | 1.1 | 0.861 |
| 55 | EPCOT_woK562:H3K4me3 signal at E | 0.0517 | 0.046-0.058 | 0.052 | 0.0494 | 1.1 | 0.861 |
| 56 | EPCOT:Bru-seq signal at E | 0.0504 | 0.046-0.057 | 0.058 | 0.0377 | 1.1 | 0.861 |
| 57 | phyloP:max_phyloP | 0.0502 | 0.044-0.057 | 0.050 | 0.0403 | 1.1 | 0.861 |
| 58 | EPCOT:NET-CAGE forward signal at P | 0.0486 | 0.042-0.060 | 0.043 | 0.0493 | 1.1 | 0.861 |
| 59 | pLI_LOEUF:pLI | 0.0480 | 0.042-0.055 | 0.047 | 0.0434 | 1.1 | 0.861 |
| 60 | pLI_LOEUF:LOEUF | 0.0479 | 0.043-0.056 | 0.047 | 0.0430 | 1.1 | 0.861 |
| 61 | Motif:UniqueMotifDensityJaspar2026 | 0.0478 | 0.043-0.059 | 0.046 | 0.0441 | 1.1 | 0.861 |
| 62 | EPCOT_woK562:H3K9me3 signal at P | 0.0477 | 0.042-0.054 | 0.051 | 0.0412 | 1.1 | 0.861 |
| 63 | phastCons:mean_phastCons | 0.0472 | 0.043-0.053 | 0.050 | 0.0413 | 1.0 | 0.861 |
| 64 | ChromHMM:CrhmmBool-Acetylated_Regions | 0.0470 | 0.043-0.052 | 0.045 | 0.0442 | 1.0 | 0.861 |
| 65 | ChromHMM:CrhmmBool-Enh_ESC | 0.0464 | 0.042-0.051 | 0.045 | 0.0447 | 1.0 | 0.861 |
| 66 | EPCOT_woK562:H3K27me3 signal at P | 0.0461 | 0.038-0.057 | 0.043 | 0.0514 | 1.0 | 0.861 |
| 67 | EPCOT_woK562:RNA-seq signal at E | 0.0460 | 0.042-0.051 | 0.049 | 0.0418 | 1.0 | 0.861 |
| 68 | TFgene:is_tf | 0.0458 | 0.042-0.051 | 0.045 | 0.0452 | 1.0 | 0.861 |
| 69 | Alu:OverlapsAluS | 0.0455 | 0.041-0.050 | 0.045 | 0.0455 | 1.0 | 0.861 |
| 70 | ChromHMM:CrhmmBool-Heterochromatin_Regions | 0.0453 | 0.041-0.049 | 0.045 | 0.0457 | 1.0 | 0.861 |
| 71 | Alu:OverlapsAluY | 0.0453 | 0.041-0.049 | 0.045 | 0.0458 | 1.0 | 0.861 |
| 72 | Alu:OverlapsAlu | 0.0451 | 0.041-0.049 | 0.045 | 0.0459 | 1.0 | 0.861 |
| 73 | ChromHMM:CrhmmBool-Quiescent_Regions | 0.0450 | 0.041-0.049 | 0.045 | 0.0460 | 1.0 | 0.861 |
| 74 | Alu:OverlapsAluJ | 0.0448 | 0.041-0.049 | 0.045 | 0.0466 | 1.0 | 0.861 |
| 75 | ChromHMM:CrhmmBool-Enh_Mesenchyme | 0.0448 | 0.041-0.049 | 0.045 | 0.0462 | 1.0 | 0.861 |
| 76 | ChromHMM:CrhmmBool-Transcribed_Region | 0.0448 | 0.041-0.049 | 0.045 | 0.0475 | 1.0 | 0.861 |
| 77 | ChromHMM:CrhmmBool-Enh_Transcribed | 0.0448 | 0.041-0.049 | 0.045 | 0.0466 | 1.0 | 0.861 |
| 78 | EPCOT:H3K36me3 signal at E | 0.0447 | 0.040-0.050 | 0.046 | 0.0483 | 1.0 | 0.861 |
| 79 | ChromHMM:CrhmmBool-Promoter_Bivalent | 0.0446 | 0.041-0.049 | 0.045 | 0.0474 | 1.0 | 0.861 |
| 80 | ChromHMM:CrhmmBool-Enh_Weak | 0.0446 | 0.041-0.048 | 0.045 | 0.0466 | 1.0 | 0.861 |
| 81 | ChromHMM:CrhmmBool-RepressivePolycomb_Regions | 0.0444 | 0.040-0.048 | 0.045 | 0.0471 | 1.0 | 0.861 |
| 82 | ChromHMM:CrhmmBool-Promoter_Active | 0.0443 | 0.040-0.049 | 0.045 | 0.0469 | 1.0 | 0.861 |
| 83 | ChromHMM:CrhmmBool-TSS | 0.0443 | 0.040-0.048 | 0.045 | 0.0478 | 1.0 | 0.861 |
| 84 | ChromHMM:CrhmmBool-Other | 0.0441 | 0.040-0.048 | 0.045 | 0.0493 | 1.0 | 0.861 |
| 85 | phyloP:min_phyloP | 0.0439 | 0.039-0.049 | 0.047 | 0.0423 | 1.0 | 0.861 |
| 86 | Motif:MotifDensityJaspar2026 | 0.0438 | 0.040-0.051 | 0.045 | 0.0451 | 1.0 | 0.861 |
| 87 | Motif:CTCFMotifHitJaspar2026 | 0.0437 | 0.040-0.048 | 0.045 | 0.0508 | 1.0 | 0.861 |
| 88 | EPCOT_woK562:H3K9me3 signal at E | 0.0436 | 0.039-0.051 | 0.038 | 0.0651 | 1.0 | 0.861 |
| 89 | EPCOT:RNA-seq signal at E | 0.0433 | 0.040-0.048 | 0.049 | 0.0426 | 0.9 | 0.861 |
| 90 | Motif:UniqueMotifCountsJaspar2026 | 0.0428 | 0.038-0.048 | 0.047 | 0.0443 | 0.9 | 0.861 |
| 91 | EPCOT:NET-CAGE reverse signal at P | 0.0426 | 0.037-0.051 | 0.042 | 0.0500 | 0.9 | 0.861 |
| 92 | EPCOT_woK562:H3K36me3 signal at E | 0.0423 | 0.037-0.048 | 0.041 | 0.0526 | 0.9 | 0.861 |
| 93 | EPCOT_woK562:H4K20me1 signal at P | 0.0423 | 0.038-0.048 | 0.044 | 0.0535 | 0.9 | 0.861 |
| 94 | EPCOT_woK562:NET-CAGE reverse signal at P | 0.0423 | 0.037-0.052 | 0.041 | 0.0530 | 0.9 | 0.861 |
| 95 | Motif:MotifCountsJaspar2026 | 0.0420 | 0.037-0.047 | 0.046 | 0.0450 | 0.9 | 0.861 |
| 96 | sHet:genebayes_shet | 0.0413 | 0.038-0.046 | 0.046 | 0.0447 | 0.9 | 0.861 |
| 97 | EPCOT_woK562:NET-CAGE forward signal at P | 0.0412 | 0.036-0.050 | 0.043 | 0.0502 | 0.9 | 0.861 |
| 98 | EPCOT_woK562:EP300 signal at P | 0.0408 | 0.036-0.049 | 0.040 | 0.0630 | 0.9 | 0.861 |
| 99 | EPCOT:EP300 signal at P | 0.0407 | 0.035-0.051 | 0.040 | 0.0722 | 0.9 | 0.861 |
| 100 | EPCOT_woK562:H2AFZ signal at P | 0.0401 | 0.036-0.045 | 0.045 | 0.0473 | 0.9 | 0.861 |
| 101 | ENCODEstats:std | 0.0397 | 0.036-0.044 | 0.042 | 0.0509 | 0.9 | 0.861 |
| 102 | EPCOT:H3K4me2 signal at P | 0.0392 | 0.034-0.046 | 0.043 | 0.0601 | 0.9 | 0.861 |
| 103 | EPCOT:H2AFZ signal at P | 0.0383 | 0.034-0.043 | 0.044 | 0.0522 | 0.8 | 0.861 |
| 104 | ENCODEstats:mean | 0.0375 | 0.033-0.045 | 0.041 | 0.0560 | 0.8 | 0.861 |
| 105 | EPCOT_woK562:Bru-seq signal at P | 0.0374 | 0.034-0.042 | 0.042 | 0.0513 | 0.8 | 0.861 |
| 106 | EPCOT:RNA-seq signal at P | 0.0374 | 0.034-0.043 | 0.042 | 0.0546 | 0.8 | 0.861 |
| 107 | EPCOT:H4K20me1 signal at P | 0.0373 | 0.033-0.043 | 0.041 | 0.0552 | 0.8 | 0.861 |
| 108 | EPCOT_woK562:H3K27me3 signal at E | 0.0372 | 0.033-0.041 | 0.040 | 0.0605 | 0.8 | 0.861 |
| 109 | EPCOT:H3K4me3 signal at P | 0.0372 | 0.034-0.042 | 0.042 | 0.0553 | 0.8 | 0.861 |
| 110 | EPCOT:GRO-cap reverse signal at P | 0.0371 | 0.034-0.043 | 0.041 | 0.0569 | 0.8 | 0.861 |
| 111 | EPCOT_woK562:RNA-seq signal at P | 0.0364 | 0.033-0.040 | 0.040 | 0.0544 | 0.8 | 0.861 |
| 112 | EPCOT_woK562:H3K4me3 signal at P | 0.0363 | 0.033-0.040 | 0.043 | 0.0531 | 0.8 | 0.861 |
| 113 | EPCOT:H3K79me2 signal at P | 0.0363 | 0.033-0.041 | 0.044 | 0.0562 | 0.8 | 0.861 |
| 114 | EPCOT:GRO-seq reverse signal at P | 0.0360 | 0.033-0.041 | 0.041 | 0.0555 | 0.8 | 0.861 |
| 115 | EPCOT:H3K27me3 signal at E | 0.0360 | 0.032-0.040 | 0.040 | 0.0659 | 0.8 | 0.861 |
| 116 | EPCOT:H3K9ac signal at P | 0.0357 | 0.032-0.040 | 0.041 | 0.0611 | 0.8 | 0.861 |
| 117 | EPCOT:Bru-seq signal at P | 0.0355 | 0.033-0.039 | 0.041 | 0.0553 | 0.8 | 0.861 |
| 118 | EPCOT:GRO-cap forward signal at P | 0.0354 | 0.032-0.039 | 0.042 | 0.0544 | 0.8 | 0.861 |
| 119 | EPCOT_woK562:GRO-seq reverse signal at P | 0.0353 | 0.032-0.040 | 0.041 | 0.0556 | 0.8 | 0.861 |
| 120 | EPCOT:GRO-seq forward signal at P | 0.0352 | 0.032-0.039 | 0.043 | 0.0535 | 0.8 | 0.861 |
| 121 | EPCOT_woK562:GRO-seq forward signal at P | 0.0348 | 0.032-0.038 | 0.041 | 0.0558 | 0.8 | 0.861 |
| 122 | EPCOT_woK562:H3K79me2 signal at P | 0.0348 | 0.032-0.038 | 0.040 | 0.0550 | 0.8 | 0.861 |
| 123 | EPCOT:H4K20me1 signal at E | 0.0346 | 0.031-0.038 | 0.038 | 0.0739 | 0.8 | 0.861 |
| 124 | EPCOT:H3K27ac signal at P | 0.0341 | 0.031-0.038 | 0.040 | 0.0639 | 0.8 | 0.861 |
| 125 | EPCOT_woK562:GRO-cap reverse signal at P | 0.0336 | 0.030-0.038 | 0.041 | 0.0609 | 0.7 | 0.861 |
| 126 | EPCOT_woK562:GRO-cap forward signal at P | 0.0334 | 0.030-0.037 | 0.040 | 0.0601 | 0.7 | 0.861 |
| 127 | EPCOT_woK562:H3K27ac signal at P | 0.0333 | 0.030-0.037 | 0.039 | 0.0852 | 0.7 | 0.861 |
| 128 | EPCOT_woK562:H3K9ac signal at P | 0.0333 | 0.030-0.037 | 0.041 | 0.0612 | 0.7 | 0.861 |
| 129 | EPCOT_woK562:H4K20me1 signal at E | 0.0326 | 0.029-0.036 | 0.037 | 0.0860 | 0.7 | 0.861 |

## Reading the table

- A predictor at the random baseline carries no information about which tested elements regulate their gene.
- `pairs overlapped` below 1.0 means some CRISPR elements had no predicted element for that gene in the table; those pairs scored `fill_value` and count against the predictor, exactly as in CRISPR_comparison.
- Bootstrap intervals resample CRISPR pairs with replacement; overlapping intervals mean the ranking between two predictors is not settled by this dataset.

## Files

- `benchmark_summary.tsv`
- `auprc_ranked.png`

