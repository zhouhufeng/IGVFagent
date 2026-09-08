# Saturation genome editing (`igvfagent sge`)

SGE is **not** RNA-seq. Its reads are a fixed amplicon covering one exon
tile, so quantifying them against a transcriptome answers a different
question — on IGVFDS4629JYPY that produced "99.86% of counts map to PALB2",
which is the amplicon design restated, not a result. The same number comes
back however the experiment turned out.

What SGE measures is the **fate of each programmed variant**: variants that
damage the protein are depleted from the cell population between an early
and a late timepoint. The score is that depletion.

```bash
igvfagent sge design  IGVFDS4629JYPY          # which targets does the library cover?
igvfagent sge analyze IGVFDS4629JYPY          # find design + early mate, score
igvfagent sge count   IGVFDS4629JYPY          # one sample's variant composition
igvfagent sge score --late IGVFDS4629JYPY --early IGVFDS5691KKNE
```

## How a score is produced

1. **Design** — the MeasurementSet links a construct library whose
   `integrated_content_files` include an *editing templates* TSV: per target,
   the amplicon coordinates (`ampstart`/`ampstop`), the edited window
   (`editstart`/`editstop`), and the `required_edits` every template carries.
2. **Reference** — the amplicon sequence is fetched from Ensembl for those
   coordinates and cached under `Data/SGE/_reference`.
3. **Variant calling** — each read is anchored to the amplicon and compared
   position-wise. No aligner is needed: the amplicon is a fixed window and
   reads start at known anchors, which is both exact and easier to reason
   about than a general aligner's gap placement.
4. **Early mate** — the sample alias
   (`lea-starita:PALB2-X7A-replicate2-day12`) identifies gene, target and
   replicate, so the matching earliest timepoint is found automatically.
5. **Score** — `log2(late_freq / early_freq)` per variant. **Negative =
   depleted = damaging.**

Frequencies, not raw counts: the two libraries are sequenced to different
depths, and a raw ratio would measure that instead of biology.

## What is excluded from scores, and why

- **`required_edits`** sit in essentially every read by construction (PAM
  disruption and the like). They are not under selection and would anchor
  the distribution; they are reported in `all_variants.tsv` and dropped from
  `functional_scores.tsv`.
- **Variants outside the edit window** — not programmed, so not scoreable.
- **Reads below 90% identity** are counted as `low_identity`, not forced
  into a substitution list. Only substitutions are called; a position-wise
  comparison cannot place an indel, and pretending otherwise would invent
  variants.

## Measured on PALB2 X7A (replicate 2, day 12 vs day 5)

```
355 variants scored     median log2 −0.10     IQR −0.32 … 0.17
19 depleted (< −1)       1 enriched (> 1)      range −2.03 … 1.11
reads parsed: 393,451 late / 279,826 early
```

Both of the design's `required_edits` were recovered, which is the check
that the amplicon reference is aligned correctly — get the coordinates wrong
and they vanish. The depleted variants cluster at chr16:23626291–23626328
rather than scattering, which is what a real functional signal looks like.

Outputs under `Docs/SGE/<label>/`: `functional_scores.tsv`,
`all_variants.tsv`, `summary.json`, `sge_scores.png` (score distribution and
score-along-amplicon) and `sge_coverage.png`.

## Limits

- **One timepoint gives counts, not scores.** `analyze` says so and falls
  back to counting rather than presenting composition as function.
- **Substitutions only** — no indel calling.
- **No amino-acid annotation yet.** Scores are per nucleotide position; a
  missense/nonsense/synonymous breakdown needs the transcript reading frame
  mapped onto the amplicon.
- **No replicate combination.** Each `analyze` scores one replicate pair;
  the published PALB2 analysis combines several.
