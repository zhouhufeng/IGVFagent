# Skill: IGVF Front-Page Portal And Knowledge Graph Summary

Use this skill when a user asks for an IGVF data overview, a current Portal/KG summary, or a README front-page stats refresh.

## Command

```bash
python3 Scripts/igvf_frontpage_summary.py refresh --update-readme
```

## What It Summarizes

- IGVF Portal object-type totals from the configured portal API.
- IGVF Catalog API smoke evidence classes from the configured catalog API.
- IGVF Knowledge Graph / ArangoDB collection counts when `IGVF_ARANGO_PASSWORD` is set locally.

## Outputs

- Root `README.md` generated block between `IGVF_FRONT_PAGE_STATS_START` and `IGVF_FRONT_PAGE_STATS_END`.
- Stable Markdown report: `Docs/IGVF_FRONT_PAGE_DATA_SUMMARY.md`.
- Stable machine-readable JSON: `Data/Summaries/igvf_frontpage_summary.json`.
- Timestamped JSON snapshots in `Data/Summaries/`.
- Runtime logs in `Docs/Logs/`.

## Update Policy

Run manually before releases, demos, or pushes to GitHub. For periodic refreshes on a workstation, use cron or another scheduler from the repo root:

```bash
cd /path/to/IGVFdataAgent
python3 Scripts/igvf_frontpage_summary.py refresh --update-readme
```

Do not store IGVF Portal cookies, ArangoDB passwords, or other credentials in the repository.
