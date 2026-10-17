# IGVF Portal data submission in IGVF Agent: progress

2026-10-17 · zhouhufeng

## Summary

Submitters can now check, audit, create and upload IGVF Portal objects with one tool, from the command line or the web app. Writes rehearse on staging first and send only on an explicit `--apply` or the user's own Submit click.

- **CLI:** `igvfagent submit` has 23 commands covering every stage, from key checks and object ordering to creating, uploading, patching and auditing.
- **The Portal's own audits** now run locally, on Portal objects and on unsubmitted drafts. They reproduce every finding a client can check (327 of 327) and report none the Portal does not.
- **Web app (igvfagent.genohub.org):** users paste their own Portal keys, upload data to a private folder, let the agent dry-run the submission and click **Submit to staging**. Prod submission is switched off.
- **Agent tools:** nine read-only `igvf_submit_*` tools, so the agent can audit and prepare but never write.

Status on 2026-10-17: everything below is deployed to igvfagent.genohub.org. This public repository's code is updated separately; the newest pieces (`sheet-*`, `portal-audit`, the web key store) will follow in a later release.

## Timeline

| Date | What landed |
|---|---|
| 2026-09-10 | IGVF-Submit absorbed as the `submit` and `portal-qc` skills: `doctor`, `audit` (24 DACC-meeting rules), `inputs`, `crosscheck`, `validate`, `template`, `enum`, `patch` |
| 2026-09-11 | `plan` (dependency order), `check-revoked`, `status` |
| 2026-09-23 | IGVF data model from the lab submission diagrams; lineage context |
| 2026-10-15 | `portal-qc`: provenance by reverse `derived_from` lookup (6 min instead of a full scan), lab slug match, CSV output |
| 2026-10-16 | `submit create`: POSTs drafts in dependency order, validates against the live schema, reuses existing objects, uploads file contents, per-host Portal keys |
| 2026-10-16 | Web submission: keys pasted in chat are removed before the model sees them and stored encrypted per user; "My data" folder; Submit-to-staging button |
| 2026-10-16 | Fix: audit file checks were off after `create` reused the `FILE_TYPES` name; `patch` waits until the Portal shows the write |
| 2026-10-16 | IGVF-DACC igvf-metadata-submitter (Google Sheets) ported as `sheet-template`, `sheet-get`, `sheet-post`, `sheet-patch`, `sheet-put`, `sheet-upload` for TSV/CSV/XLSX |
| 2026-10-16 | Submitter overview, `Docs/Logs/SubmissionSkills/IGVF-Agent-submission-skills.md` |
| 2026-10-17 | igvfd's own audit checkers ported as `submit portal-audit`, with staging and prod fixtures and an agreement report |

## What exists now

### Command line: `igvfagent submit`

| Stage | Commands |
|---|---|
| Set up | `doctor` (keys, host, account, labs it may submit for) |
| Prepare | `plan`, `template`, `enum`, `sheet-template`, `sheet-get` |
| Check before upload | `validate` (gzip, md5, ragged rows, BED coordinates), `portal-audit --drafts` |
| Submit | `create` (drafts and `iu_register`-style sheets, `--upload`), `sheet-post`, `sheet-upload` |
| Fix | `patch` (incl. `--repoint`), `sheet-patch`, `sheet-put` (admin) |
| Check after submission | `audit` (24 DACC rules), `portal-audit` (the Portal's own checkers, `--compare` against `@@audit-self`), `inputs`, `check-revoked`, `crosscheck`, `status`, `get`, `fileset-files` |
| Portal-wide QC | `igvfagent portal-qc summary / audits / provenance / uploads / filesets / lab` |

Every write is a dry run unless `--apply` is given, and prod also needs `--yes`. The CLI needs Python 3.9 and no packages.

### The Portal's own audits (`portal-audit`)

The igvfd audit checkers (`src/igvfd/audit`) are ported unchanged into `Scripts/igvfd_audit_rules.py`. `Scripts/igvfd_audits.py` gives them what the server would: `request.embed` becomes a GET on the Portal API, or a lookup of a draft by alias.

Agreement with the Portal's computed audits (`Benchmarks/submission_staging/igvfd_audit_agreement.md`):

| | Objects | Portal findings | Matched | Missed | Extra | Recall | Precision |
|---|---|---|---|---|---|---|---|
| All | 186 | 345 | 327 | 18 | 0 | 94.8% | 100% |
| prod | 127 | 273 | 265 | 8 | 0 | 97.1% | 100% |
| staging | 59 | 72 | 62 | 10 | 0 | 86.1% | 100% |

All 18 misses need data a client cannot have:
- 14 incorrect-file-bucket findings read the server's S3 storage.
- 3 ontology findings need the Portal's ontology database.
- 1 status finding concerns a linked object these keys cannot read; it is reported as `UNVERIFIABLE`, not as a pass.

Tests: `Benchmarks/test_igvfd_audits.py`, with fixtures in `Benchmarks/submission_staging/igvfd_fixtures/`.

### Web app

1. The user pastes their Portal Access Key ID and Secret. `Scripts/portal_keystore.py` removes them from the message before the model or any transcript sees them, checks which host they work on, and stores them encrypted under `IGVF_KEYSTORE_SECRET`.
2. The user uploads data and sheets to **My data**, a private per-user folder of up to 100 MB per file.
3. The agent checks, audits and dry-runs `create` with the user's own keys.
4. The user clicks **Submit to staging**. Only that click writes. Prod needs `IGVF_ALLOW_PROD_SUBMIT=1` on the server and is off.

### Agent tools

`igvf_submit_audit`, `igvf_submit_inputs`, `igvf_submit_crosscheck`, `igvf_submit_validate`, `igvf_submit_fileset_files`, `igvf_submit_template`, `igvf_submit_sheet_template`, `igvf_submit_sheet_get`, `igvf_submit_create` (dry run only). None of them can write to the Portal.

### Verified on staging

Lab `xihong-lin`, award `HG012064`:
- Objects posted by `create`: `TSTDS63011942`, `TSTFI69604689`, and one document. They match what the official `iu_register` posted (`TSTDS45289259`, `TSTFI32661857`) on every submitted field.

## Open items

1. **`portal-audit` is not an agent tool yet.** The CLI has it, but the agent still audits with the 24-rule `igvf_submit_audit`. The next step is an `igvf_submit_portal_audit` tool (read-only), and making the web flow run it on drafts before the Submit button is offered.
2. **File format specification document:** writing it for the user is planned as `spec-doc`.
3. **Prod submission from the web app:** off until staging use by real submitters has been reviewed.
4. **Large raw data through the web app:** uploads are capped at 100 MB per file. Larger files go through `create --upload` or `sheet-upload` from the command line.
5. **Server-only audits:** file bucket and ontology checks cannot run client-side and are reported as skipped. `--compare` shows the Portal's own result for submitted objects.

## Where things live

| What | Where |
|---|---|
| Submission CLI | `Scripts/igvf_submission_skill.py` (`igvfagent submit`) |
| Sheet commands | `Scripts/igvf_metadata_sheet.py` |
| Portal audits | `Scripts/igvfd_audits.py`, `Scripts/igvfd_audit_rules.py` |
| Portal-wide QC | `Scripts/igvf_portal_qc_skill.py` (`igvfagent portal-qc`) |
| Key store | `Scripts/portal_keystore.py` |
| Agent tools | `Scripts/_tools.py` (`igvf_submit_*`) |
| Submitter overview | `Docs/Logs/SubmissionSkills/IGVF-Agent-submission-skills.md` |
| Submitter's guide | `Docs/Skills/igvf-submission.md` |
| Tests | `Benchmarks/test_igvf_submission.py`, `Benchmarks/test_igvfd_audits.py` |
