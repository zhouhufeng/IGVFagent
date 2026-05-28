# OPERATIONS — Buckley 2024 VHL SGE

> **URN-pending.** Same situation as `waters2024_bap1` — the MaveDB URN
> in `run.sh` is unverified; the pre-flight curl check makes the run
> exit 77 cleanly until you set `VHL_URN` to a verified scoreset URN.

For shared prerequisites + the two execution paths + concordance scoring,
see `Benchmarks/OPERATIONS_GUIDE.md`. For URN-discovery technique, see
`waters2024_bap1/OPERATIONS.md` § "How to find the real URNs".

## Quick recipe to find the right URN

```bash
curl -s 'https://api.mavedb.org/api/v1/experiments/' | \
    .venv/bin/python -c "
import json, sys
for exp in json.load(sys.stdin):
    title = (exp.get('title') or '')
    if 'VHL' in title and 'SGE' in title.upper():
        print(exp.get('urn'), '-', title[:80])
"
```

For each candidate experiment URN, open `https://www.mavedb.org/#/experiments/<urn>`
and copy the scoreset URN(s). Pre-flight with
`curl -sI "https://api.mavedb.org/api/v1/score-sets/<URN>/scores" | head -1`.

## Quick run

```bash
VHL_URN="urn:mavedb:00001183-a-1" \
    bash Benchmarks/buckley2024_vhl/run.sh

.venv/bin/python Benchmarks/concordance.py --benchmark buckley2024_vhl
```

`run.sh` also fires two `catalog` cross-references (anonymous, online):

```bash
.venv/bin/igvfagent catalog get-entity VHL
.venv/bin/igvfagent catalog find-associations VHL --relationship genetic --limit 10
```

These exercise the IGVF Catalog node + edge endpoints for the gene
(VHL is at ENSG00000134086, chr3:10,141,778-10,153,667 GRCh38) and
return its catalogued disease + phenotype associations.

## Where artefacts will land

```
Docs/MaveDB/<ts>_buckley2024_vhl/
├── VHL_mapped.tsv
├── VHL_mapped.vcf
├── summary.json
└── ...

Docs/CatalogQuery/<ts>_get_VHL/entity.json
Docs/CatalogQuery/<ts>_assoc_VHL_genetic/associations.json
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `n_rows_in ∈ [1500, 6000]` (expect ≈ 2,268) | Full VHL SGE scoreset downloaded |
| 2 | `gene = "VHL"` | Symbol resolved |
| 3 | `VHL_mapped.vcf` exists, non-empty | VCF emission worked |
| 4 | `VHL_mapped.tsv` has 500-20,000 rows | Codon mapping yielded output |

## Ground-truth spot-checks

| Signal | Expected |
|---|---|
| Chromosome | All rows show `chr3` (VHL locus, distinct from BAP1's chr3 location) |
| Position range | Approximately `10,141,000-10,154,000` GRCh38 |
| Total SNVs | 2,268 (paper) |
| Functional buckets | LOF1 / LOF2 / Intermediate / Neutral (paper Fig 2) |
| Type-2A pheo position 188 (e.g. p.Tyr98His) | Intermediate bucket score ∈ [–0.4, –0.22] |
| Type-1 truncating variants | LOF1 score < –1.26 |
| Catalog `find-associations` for `VHL` (genetic) | Should return associations to renal cell carcinoma, von Hippel-Lindau disease, pheochromocytoma |

## Running through the UI

```
Run the Buckley 2024 VHL SGE reproducibility benchmark in three steps:

1) Call mavedb_map_scoreset with urn="urn:mavedb:00001183-a-1" (or my
   verified VHL_URN), gene="VHL", label="buckley2024_vhl".

2) Call catalog_get_entity with id="VHL" to fetch the IGVF Catalog
   record for the gene.

3) Call catalog_find_associations with entity_id="VHL",
   relationship="genetic", limit=10, to cross-reference the SGE
   results against catalogued disease associations.

Report:
  - VHL_mapped.tsv row count
  - chromosome (should be chr3 around position 10.1M)
  - the IGVF Catalog's catalogued ENSG/HGNC/Entrez identifiers
  - the top 5 disease associations and their sources (Orphanet/ClinGen/etc.)
```

UI sidebar: max iterations = 12 (three tool calls + synthesis).

## Troubleshooting

Inherits everything from `waters2024_bap1/OPERATIONS.md`. Additional
VHL-specific notes:

| Symptom | Cause | Fix |
|---|---|---|
| Mapped positions look like chr3:52M (BAP1 range) | Ensembl resolved a different VHL paralog | Set `--species human` explicitly; verify `summary.json["transcript_id"]` is ENST00000256474 |
| `catalog find-associations VHL` returns 0 rows | The IGVF Catalog's mirror is mid-update | Drop `--relationship genetic` to use `--relationship all` and re-check |

## License + provenance

* **Paper data** (MaveDB scoreset): MaveDB CC-BY 4.0.
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Buckley M et al. *Nat Genet* **56**: 1446–1455 (2024).
  doi:10.1038/s41588-024-01800-z · PMID:38969834
