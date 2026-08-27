# Variant lists

Input for `igvfagent variant-list annotate`, which annotates a set of variants
against FAVOR and the IGVF Catalog and adds the results to the local knowledge
graph as variant vertices with edges to genes, diseases and regulatory
elements.

You do **not** need to prepare a CSV. Paste the variants directly:

```bash
igvfagent variant-list annotate --variants "2-21001846-G-A rs763341676"
igvfagent variant-list annotate --input Data/Input/VariantList/apob_lof_hg38.txt
igvfagent variant-list parse    --input <file>       # dry-run the parser only
```

In the browser UI, paste the list into chat and ask for annotation — the agent
calls the same skill.

## Supported notations

Auto-detected per token; a single file may mix them freely, and duplicates
across notations collapse to one variant.

| Notation | Example |
|---|---|
| `chr-pos-ref-alt` | `2-21001846-G-A` |
| with `chr` prefix | `chr2-21001846-G-A` |
| colon-separated | `2:21001846:G:A` |
| underscore-separated | `2_21002893_G_A` |
| HGVS genomic | `chr2:g.21002815G>A` |
| dbSNP rsID | `rs763341676` |
| SPDI | `NC_000002.12:21001845:G:A` |
| VCF record | `2	21001846	.	G	A` |
| CSV with columns | `rsID` / `SPDI` / `chrom,pos,ref,alt` |

Separators may be spaces, commas or newlines. Lines beginning with `#` are
treated as comments.

**Coordinate convention.** `chr-pos-ref-alt` and HGVS inputs are read as
**1-based** — the standard for these notations — and converted internally to
the Catalog's 0-based identifier. SPDI is passed through unchanged because it
is already 0-based. FAVOR is queried with the 1-based form it expects. Mixing
these up shifts every position by one, so the distinction is enforced rather
than inferred.

## Files here

| File | Contents |
|---|---|
| `apob_lof_hg38.txt` | 25 APOB loss-of-function (stopgain) variants, GRCh38, `chr-pos-ref-alt` |
| `apob_lof_rsids.txt` | the **same 25 variants** as rsIDs |
| `apob_lof_mixed_notation.txt` | the same variants written in seven different notations — the parser fixture |
| `gene-variants-APOB-LoF.tsv` | FAVOR annotation export for the 25 stopgain variants |
| `gene-variants-APOB.tsv` | FAVOR annotation export for a broader APOB variant set |
| `example_variants.csv` | minimal CSV for smoke-testing the column-based path |

The APOB set is a worked example: all 25 are stopgain variants in *APOB*,
whose loss of function underlies familial hypobetalipoproteinaemia, and
several carry ClinVar assertions for familial hypercholesterolaemia. It is
small enough to annotate in seconds and biologically coherent enough that the
resulting knowledge-graph subgraph is worth looking at.

The TSVs contain public annotation columns only (allele frequencies, GENCODE
categories, ClinVar assertions, predicted scores) — no participant-level data.

## Adding your own

Drop a file here and pass it with `--input`. Note that `.gitignore` tracks
only the example files listed above; anything else you add stays local, which
is intended — put your own variant lists here without risk of committing them.
