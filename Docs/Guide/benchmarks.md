# Reproducibility benchmarks

Published analyses that IGVFagent reproduces from public data. The up-to-date dashboard is [`Benchmarks/README.md`](../../Benchmarks/README.md).

IGVFagent ships a **reproducibility benchmark suite** in
[`Benchmarks/`](../../Benchmarks/README.md): recent Nature / Cell / Science /
Nat Genet / Nat Methods / Genome Biol papers whose published analyses IGVFagent
reproduces **directly from public data**. Each benchmark is a self-contained
directory — data sources, a deterministic `run.sh`, machine-readable
`expected.json` checks, regenerable figures, and a paper-vs-IGVFagent
`README.md` with a *Concordance / Verdict / Honest caveats* structure. A
stdlib-only scorer (`concordance.py`) turns each run into pass/fail checks.

The benchmark table lives in one place, the dashboard in
[`Benchmarks/README.md`](../../Benchmarks/README.md): each paper's headline
result, the skills it exercises, whether it ran through the IGVFagent MCP
server, and a link to its page. Each paper has one directory. Parts of a
reproduction that used to be separate benchmarks sit in subfolders of their
paper, for example `liu2025_kidney_multiome/open4gene/`.

Most benchmarks download public data and run IGVFagent's real analysis chain,
then score concordance against the authors' own results. A few need a
user-side credential or a large local download to complete; their READMEs say
which.

```bash
# Verify the suite works (~60 s)
bash Benchmarks/matreyek2018_pten_vampseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark matreyek2018_pten_vampseq

# Run all online benchmarks and score them
bash Benchmarks/run_all.sh --online-only
.venv/bin/python Benchmarks/concordance.py --all
```

---

[← Documentation index](README.md) · [Project README](../../README.md)
