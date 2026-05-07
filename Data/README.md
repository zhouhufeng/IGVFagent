# Data

Downloaded IGVF Portal, IGVF Catalog, IGVF Knowledge Graph, ENCODE, and FAVOR
outputs are written here at runtime.

Large generated files and cached responses are ignored by Git by default. Keep
only small documentation fixtures or intentional sample data under version
control.

## Input variant lists

Place your own variant CSVs under `Data/Input/VariantList/`. Filenames must
have a `.csv` extension. A tiny illustrative file is shipped at
`Data/Input/VariantList/example_variants.csv`.

**Never commit confidential or pre-publication variant lists** — they are
gitignored by default.
