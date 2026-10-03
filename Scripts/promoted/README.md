# Promoted extensions

Agent-authored extensions that a reviewer promoted to built-ins with
`igvfagent ext-review promote`. `tools/*.json` are tool manifests, and
`skills/*.py` are the modules behind them. They ship with the code, count as
built-ins (a same-named user extension is shadowed), and are checked by
`Scripts/test_promoted.py`. `registry.json` records who reviewed each one,
when, the source hash, its usage at review, and any accepted risk. See
[Docs/Guide/extending.md](../../Docs/Guide/extending.md#from-extension-to-core).
