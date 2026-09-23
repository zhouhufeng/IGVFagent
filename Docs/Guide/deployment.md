# Operating the hosted deployment

For operators of a shared server such as igvfagent.genohub.org. Server setup is in [`Deploy/README.md`](../../Deploy/README.md) and sign-in administration in [`Docs/AUTH.md`](../AUTH.md).

## Recreating the hosted container (and why you must)

**The container does not mount the source tree.** `Scripts/streamlit_app.py`
imports `from igvfagent import ...` — the package `pip install` bakes into
`/opt/venv` at image **build** time. So `git pull` on the host updates the
checkout and changes nothing about the running site, however many times you
pull. This has cost real debugging time: an external tester reported a bug
that was already fixed and pushed, because the container was serving an
image built before the fix.

### Check first, then rebuild

```bash
bash Deploy/redeploy.sh --check     # compare running code against the checkout
bash Deploy/redeploy.sh             # pull, rebuild --no-cache, recreate, verify
```

`--check` hashes the top-level modules inside the container and the same
files in the checkout and prints both. The sidebar shows that same hash, so
you can confirm what a browser is talking to without shell access.

One wording caveat: it prints `STALE` for **any** difference, in either
direction. A container carrying newer code than the checkout is also
reported stale.

### It refuses while analyses are running, on purpose

```
REFUSING to recreate: 2 analysis process(es) are running.
```

`--force-recreate` kills running work. The guard checks for
`sc-analyze`, `raw-pipeline`, `crispr-screen`, `gradient-screen`, `bean`,
`sge`, `mct`, `kg-mirror` and `mirror-kg.sh` **before** rebuilding, so a
long job is not destroyed after a ten-minute image build. Override
deliberately with `FORCE_RECREATE=1`, having decided the running work is
expendable:

```bash
FORCE_RECREATE=1 bash Deploy/redeploy.sh
```

`kg-mirror` keeps per-collection state, so interrupting it loses only the
collection in flight, not the collections already mirrored.

### Optional: install BEAN for base-editing screens

`igvfagent bean` implements BEAN's guide-assignment method independently, and
that is the part that fixes the 36.7% → 62.5% read assignment. It does **not**
reimplement BEAN's Bayesian variant/tiling model, and for effect sizes on a
base-editing screen that model is what you want. To have the real `bean`
available to the agent, build with it:

```bash
IGVF_INSTALL_CRISPR_BEAN=1 bash Deploy/redeploy.sh
```

Off by default: it pulls in torch + pyro and roughly doubles the image, and
most deployments do not analyse base-editing screens.

**If a rebuild is not available yet** (a long mirror or analysis in
flight), `bash Deploy/install-bean.sh` installs BEAN into the running
container. It works, and every step in it is a workaround for a different
layer of the container's hardening, each found the hard way:

| blocker | fix |
|---|---|
| apt's http method drops privileges; `cap_drop: ALL` forbids `setgroups` | `-o APT::Sandbox::User=root` |
| apt cannot write `/var/cache/apt/archives/partial` (`_apt`-owned, `CapEff=0`) | redirect `Dir::Cache::archives` to `/tmp` |
| pip as **root** fails: the venv is `igvf`-owned and root has no `CAP_DAC_OVERRIDE` | run pip as the owner, never root |
| no prebuilt wheel exists, and `CRISPResso2Align.pyx` breaks under Cython 3 | `cython<3` **and** `--no-build-isolation` |
| BEAN uses `np.int_t` / `np.Inf`, removed in NumPy 2 (the app runs 2.x) | its own venv with `numpy<2`, re-pinned last because a dependency overrides it |

BEAN goes in `/workspace/opt/bean-venv`, on the data volume, so it survives
container recreation; the gcc install and the `/usr/local/bin/bean` symlink
do not, so re-run the script after a recreate — or build the image with the
flag above and stop needing it.

**It cannot be added afterwards.** The runtime layer has no compiler, and
BEAN's `bean/mapping/CRISPResso2Align.pyx` needs Cython, so
`pip install crispr-bean` inside a running container fails with
`CompileError: bean/mapping/CRISPResso2Align.pyx`. The builder stage already
has `build-essential` for the `[hic]` extra, so that is where it goes — which
means this is a build-time decision, not a runtime one.

BEAN is AGPL-3.0 and IGVFagent is Apache-2.0. This **installs** it as a
separate program invoked as a subprocess — no linking, no vendoring, no
licence propagation. See [Reference GitHub
repositories](references.md#reference-github-repositories).

### Credentials and the shared-deployment settings

`Deploy/make-live.sh` pushes the IGVF Portal key pair and the ArangoDB
credentials into `Deploy/.env.prod` (mode 600) and then redeploys, verifying
inside the container afterwards. Secrets travel on **stdin**, never in argv,
because argv is visible in `ps` on the host for the life of the call.

```bash
bash Deploy/make-live.sh --check    # report only
bash Deploy/make-live.sh            # upsert credentials, redeploy, verify
```

Environment changes need a **recreate**, not a restart. Two settings matter
on a shared deployment and both default to closed:

| Variable | Default | What it allows |
|---|---|---|
| `IGVF_ALLOW_AGENT_AUTHORING` | `0` | the agent writing tools/skills this host then executes |
| `IGVF_ALLOW_UPLOAD_EXTENSIONS` | `0` | a **visitor** uploading a `.py` this host then executes |
| `IGVF_HEAVY_SLOTS` | `1` | concurrent heavy analyses before queueing |
| `IGVF_PUBLIC_MODE` | `0` | fixed backend + curated model allowlist |

If the deployment is reachable by anyone who has the password, leave the
first two at `0`.

### When a hot copy is the right answer instead

Some fixes reach a running container without a rebuild, which matters when a
multi-day job is in flight:

- **`streamlit_app.py`** — Streamlit re-reads its main script on every
  rerun, so a copied file takes effect on the next interaction.
- **Any CLI skill module** — every tool call runs `igvfagent <skill>` as a
  fresh subprocess, so a copied module is picked up immediately.
- **NOT `_agent.py`, `_llm.py` or `_tools.py`** — these are imported into the
  long-lived Streamlit process and cached in `sys.modules`. A new *tool* in
  particular will not appear in the agent's registry until the process
  restarts, even though its CLI works.

```bash
D=/opt/venv/lib/python3.11/site-packages/igvfagent
docker cp Scripts/<module>.py igvfagent-app:$D/<module>.py
docker exec igvfagent-app sh -c "rm -rf $D/__pycache__"
```

A hot copy is a **stopgap**: it lives in the container, not the image, so a
restart reverts it. Follow it with a real `redeploy.sh` at the next window.

---

[← Documentation index](README.md) · [Project README](../../README.md)
