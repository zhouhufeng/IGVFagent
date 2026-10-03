# Skill: MaveDB → genomic coordinates

Map MaveDB scoreset variants (HGVSp / per-position tile / single-letter
WT-pos-ALT) to canonical chr/pos/ref/alt VCF coordinates so MAVE /
VAMP-seq scores become cross-referenceable with ClinVar, gnomAD, GWAS
catalogues and the IGVF Catalog.

Clean-room reimplementation of [ave-dcd/dcd_mapping](https://github.com/ave-dcd/dcd_mapping)
(MIT, Arbesfeld 2023) using the public Ensembl REST API — no UTA /
SeqRepo / BLAT install required.

## Commands

```bash
# Map any VAMP-seq scoreset by gene symbol (uses MAVEDB_VAMPSEQ_CATALOG)
igvfagent mavedb map-scoreset --gene PTEN --label pten_v1

# Or by URN directly (works for any single-AA missense scoreset)
igvfagent mavedb map-scoreset --urn urn:mavedb:00000013-a-1 --gene PTEN

# One-command showcase: map + composite figure + narrative report
igvfagent mavedb showcase --gene PTEN
```

## Output schema (TSV)

| Column | Meaning |
|---|---|
| pos / wt / alt | Original MaveDB protein position + AA letters |
| kind | Original kind (missense / synonymous / nonsense / wt) |
| score / sd / se / lower_ci / upper_ci | Original MaveDB score columns |
| hgvsp | Reconstructed HGVSp (e.g. `p.M1L`) |
| mapping_type | missense / synonymous / nonsense / wt_mismatch / error |
| chr / pos_genomic / ref / alt_nt | Genomic coordinates + nt change |
| transcript_id | Ensembl canonical transcript |
| codon / codon_pos / strand | Codon context |
| candidate_idx / n_candidates | For amino acids with multiple alt codons, one row per candidate |
| notes | Per-row provenance string |

## VCF schema

Standard VCF-4.2 with INFO fields:
`AA_REF`, `AA_ALT`, `AA_POS`, `TRANSCRIPT`, `SCORE`, `AMBIG`.

## How it works (clean-room algorithm)

1. **Gene → Ensembl canonical transcript**
   `/xrefs/symbol/{species}/{gene}` then `/lookup/id/{gene_id}?expand=1`
   → pick the transcript with `is_canonical=true`.
2. **Protein position → genomic codon coordinates**
   `/map/translation/{transcript_id}/{aa}..{aa}` returns the 3-nt codon's
   chr / start / end / strand.
3. **Codon validation**
   Fetch the genomic 3-nt substring (`/sequence/region/{species}/{chr:start..end}:1`),
   reverse-complement if strand=-1, translate via the standard table,
   verify against the MaveDB WT amino acid.
4. **Single-nt change enumeration**
   For each position in the WT codon, find single-nucleotide changes
   that yield the alt amino acid. Emit one output row per candidate
   with `candidate_idx` + `n_candidates` so downstream consumers can
   resolve ambiguity (`n_candidates=1` means unique; >1 means multiple
   single-nt changes give the same protein change).
5. **Strand-aware genomic coordinates**
   Forward strand: pos = codon_start + codon_pos. Reverse strand:
   pos = codon_end - codon_pos, ref/alt complemented.

## Caching

Every Ensembl REST response is cached under `Data/Cache/Ensembl/`. A
fresh 4,408-variant PTEN scoreset maps in ~60 s the first time, ~5 s on
rerun.

## What's NOT yet supported

- Frameshift (`p.Xfs`), large deletions, duplications, insertions
- Non-canonical transcripts (we pick the Ensembl-canonical one)
- Multi-codon-substitution variants
- Mapping when the MaveDB target sequence doesn't match the canonical
  protein (protein-offset detection — TODO)

For these cases, fall back to the upstream `dcd-mapping` PyPI package
(MIT) which handles all of the above via UTA + SeqRepo + BLAT.
