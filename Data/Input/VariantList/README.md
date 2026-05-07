# Variant lists

Drop your variant CSV files here and pass them to the analysis skills with
`--input Data/Input/VariantList/<your_file>.csv`.

A minimal `example_variants.csv` is included for smoke testing. Replace it
with your own list.

## Expected columns

Variant identity may be supplied via any of:

- `rsid` / `rsID` / `dbSNP` — e.g. `rs58658771`
- `chrom`, `pos`, `ref`, `alt` — GRCh38 coordinates and alleles
- `spdi` — NCBI SPDI string
- `hgvs` — e.g. `NC_000019.10:g.44908822C>T`

Additional user columns (locus name, phenotype, prior, notes) are preserved
and passed through to the annotated output.

## Privacy

Do **not** commit confidential or pre-publication variant lists. The default
`.gitignore` excludes everything in this folder except `README.md` and
`example_variants.csv`.
