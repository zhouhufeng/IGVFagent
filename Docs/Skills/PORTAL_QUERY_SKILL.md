# Skill: IGVF Portal canonical-query layer

A faceted, filter-aware, paginate-safe client for the IGVF Portal that
matches the same query patterns the IGVF Data Coordinating Center (DACC)
documents as canonical in [IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp)
(MIT, Cherry Lab / IGVF DACC, 2024).

Clean-room reimplementation under Apache-2.0 — we re-derive the HTTP wire
contract from the public portal documentation and the upstream README.
No source code is copied. The skill is stdlib-only (`urllib`); it does
**not** depend on the official `igvf-client` PyPI package.

## Commands

```bash
# Search with field_filters
igvfagent portal search --type MeasurementSet \
                          --field-filters "lab.@id=/labs/j-michael-cherry;file_format=bam,bigwig" \
                          --limit 25

# Get one item by accession / @id / UUID
igvfagent portal get IGVFDS3909HJKS

# JSON schema for an ItemType
igvfagent portal schema MeasurementSet

# Canonical ItemType enum
igvfagent portal list-types

# collection_param ↔ search_field map (the UX trick)
igvfagent portal endpoint-params measurement-sets

# Facets-only call (limit=0)
igvfagent portal facets --type MeasurementSet --query CRISPR

# Tabular TSV /report.tsv
igvfagent portal report --type SequenceFile \
                          --field-filters "file_format=fastq" \
                          --limit 500 --label fastq_v1

# /batch-download/ for FileSet types (manifest of pre-signed S3 URLs)
igvfagent portal batch-download --type AnalysisSet \
                                  --field-filters "status=released" \
                                  --label analysis_v1
# add --fetch to actually pull every file referenced in the manifest
```

## Authentication

| Env var | Role |
|---|---|
| `IGVF_ACCESS_KEY` + `IGVF_SECRET_ACCESS_KEY` | DACC-blessed HTTP Basic credentials (preferred) |
| `IGVF_PORTAL_COOKIE` | Legacy cookie auth (still honored) |
| _(neither)_ | Anonymous — public-released items only |

Without auth the portal returns 403 on access-restricted items; the
skill still works for browsing released data.

## The field_filters DSL

A small, predictable filter language matching the DACC search syntax:

| Clause | Wire form |
|---|---|
| `lab.@id=/labs/jb-cherry` | `lab.@id=/labs/jb-cherry` |
| `lab.@id!=/labs/x` | `lab.@id%21=/labs/x`  (`!=` negation) |
| `file_size=gte:1000000` | `file_size=gte%3A1000000` |
| `file_format=bam,bigwig` | `file_format=bam&file_format=bigwig` (repeated param) |
| `status=released,archived` | `status=released&status=archived` |

Multiple clauses joined with `;`:

```
lab.@id=/labs/x;file_format=bam,bed;file_size=gte:1000000
```

## What this skill adds over the legacy `data` / `frontpage` skills

| Capability | Before | After |
|---|---|---|
| HTTP Basic auth | ❌ (cookie only) | ✓ |
| Field-filter DSL with negation + ranges | ❌ | ✓ |
| `/report.tsv` export | ❌ | ✓ |
| `/batch-download/` manifest + fetch | ❌ | ✓ |
| Facets-only call | partial | ✓ |
| Endpoint-param introspection | ❌ | ✓ |
| Canonical ItemType enum | ❌ | ✓ |

## License posture

Apache-2.0. Clean-room reimplementation of IGVF-DACC/igvf-portal-mcp
(MIT, Cherry Lab) — we cite the upstream README + the public portal
HTTP contract as factual sources. The upstream MIT licence would
permit verbatim copy with attribution, but we re-implement to keep the
stdlib-only posture (no `igvf-client==121.0.0` pin).
