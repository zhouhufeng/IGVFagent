# Human Transcription Factors (HTF) database skill

Local SQLite mirror of the Lambert/Jolma/Hughes **Human Transcription
Factors** database v1.01 (https://humantfs.ccbr.utoronto.ca; Lambert et
al., *Cell* 2018, doi:10.1016/j.cell.2018.01.029).

2,765 assessed proteins, of which **1,639 are curated TFs**, each with
its DNA-binding-domain family, binding mode, motif status and
cross-references. This is the TF list the Tabula Sapiens 2.0 paper uses,
and the reference any regulon or cell-type-specificity analysis needs.

## Build it once

```bash
igvfagent humantfs build-db              # ~2 MB, seconds
igvfagent humantfs build-db --with-pwms  # also fetch position weight matrices
igvfagent humantfs info
```

Lands at `Data/HumanTFs/human_tfs.sqlite`. The build refuses to
overwrite without `--force`, and warns if the TF count is not 1,639 --
the number every downstream claim is calibrated to.

## Query it

```bash
igvfagent humantfs is-tf --genes "FOXP3, EOMES, ACTB, GAPDH, IKZF1"
igvfagent humantfs lookup --genes SATB2
igvfagent humantfs list --dbd "C2H2 ZF" --limit 20
igvfagent humantfs list --with-motif --out known_motif_tfs.tsv
igvfagent humantfs motifs --gene GATA1
igvfagent humantfs families
igvfagent humantfs export --format symbol --out tf_symbols.txt
igvfagent humantfs export --format ensembl
```

`is-tf` is the workhorse: paste any gene list and it partitions into
TFs, non-TFs and *unassessed* — the third bucket matters, because "not
in the database" is not the same as "not a TF".

## As a library

```python
from igvfagent.humantfs_skill import tf_symbols, tf_table
tfs = tf_symbols()          # {'AATF', 'ABL1', ...}
rec = tf_table()['GATA1']   # dbd, binding_mode, motif_status, ...
```

## Notes

**Zinc fingers dominate.** C2H2 ZF is by far the largest DBD family, and
a large share of those TFs are computationally predicted with no
validated motif — which is exactly why the Tabula Sapiens paper flags
ubiquitous zinc-finger TFs as an understudied group. `list --dbd "C2H2
ZF"` plus `--with-motif` gives you the split.

**Motif status is not binding evidence.** `Known motif` means a motif
exists in CIS-BP for that TF, from direct or inferred evidence; it says
nothing about whether the TF is active in your cells. Pair it with
regulon activity (`igvfagent tabula tf-regulons`) before claiming
function.

**The data is downloaded, not vendored.** This repository stores no HTF
content; `build-db` fetches it at run time. Cite Lambert 2018 and honour
the upstream terms when redistributing derived tables.
