# SPLiT-seq barcode whitelists

The three rounds of ligation barcodes used in Rosenberg et al.,
*Science* 2018 (Split-pool ligation-based transcriptome sequencing).
Each round has 96 8-mer barcodes addressed by 96-well-plate position
(A1–H12). The combined cell barcode is Rd1+Rd2+Rd3 (24 bp); cells are
identified by the Cartesian product (96³ ≈ 884K theoretical cells per
sublibrary).

Vendored from the Chipeyown SDAT toolkit
(https://github.com/Chipeyown/SPLiT-seq-Data-Analysis_Toolkit), MIT-
licensed. Sequences match the Parse Biosciences Evercode WT v1
chemistry.

| File | What |
|---|---|
| `barcodes_rd1.tsv` | 96 Rd1 ligation barcodes (oligo-dT wells 1–48, random-hexamer wells 49–96; the two halves of each well share the same sample). |
| `barcodes_rd2.tsv` | 96 Rd2 ligation barcodes. |
| `barcodes_rd3.tsv` | 96 Rd3 ligation barcodes. |

Columns: `well_id` (A1–H12), `well_index` (1–96), `primer_type`,
`barcode_seq` (the 8-mer used for cell-barcode reconstruction),
`oligo_sequence` (full oligo as ordered).
