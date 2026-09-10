#!/usr/bin/env bash
# Install crispr-bean (BEAN) into a RUNNING igvfagent container.
#
#   bash Deploy/install-bean.sh            # do it
#   bash Deploy/install-bean.sh --check    # report only
#
# WHEN TO USE THIS. Prefer the image: `IGVF_INSTALL_CRISPR_BEAN=1 bash
# Deploy/redeploy.sh` builds BEAN in the builder stage, which is durable and
# needs none of the workarounds below. This script exists for the case where
# a rebuild is not available yet -- a long mirror or analysis in flight -- and
# BEAN is needed for testing today.
#
# WHY EACH STEP LOOKS ODD. Five attempts, each blocked by a different layer of
# the container's hardening. All five workarounds are load-bearing:
#
#  1. apt's http method drops privileges via setgroups, which `cap_drop: ALL`
#     plus `no-new-privileges` forbids -> "Operation not permitted", rc=100.
#     Fixed by APT::Sandbox::User=root.
#
#  2. apt cannot write /var/cache/apt/archives/partial (owned by _apt, and
#     CapEff=0 so even uid 0 has no CAP_DAC_OVERRIDE) -> every .deb download
#     fails EACCES. Fixed by redirecting Dir::Cache::archives to /tmp.
#
#  3. pip must NOT run as root. /opt/venv and the target tree are owned by
#     igvf:igvf 755, and with CapEff=0 root cannot write a directory it does
#     not own -- while igvf, the owner, can. Running pip as root fails with
#     "Permission denied: .../cython.py".
#
#  4. crispr-bean has NO prebuilt wheel (`pip download --only-binary` returns
#     "from versions: none"), so a C compiler is unavoidable, and its
#     bean/mapping/CRISPResso2Align.pyx needs Cython < 3 ("'uint_t' is not a
#     type identifier" under Cython 3) with build isolation OFF -- otherwise
#     pip builds in an isolated env with its own Cython 3 and ignores the pin.
#
#  5. BEAN needs NumPy 1.x: the .pyx uses np.int_t/np.uint_t and the code uses
#     np.Inf, all removed in NumPy 2.0. The app's own numpy is 2.x and scanpy
#     depends on it, so BEAN gets its OWN venv -- it is invoked as a
#     subprocess and does not need the app's interpreter. A dependency also
#     re-installs numpy 2.x AFTER the pin, so numpy<2 is reinstalled last.
#
# The venv lives under /workspace, a bind mount to the data volume, so it
# SURVIVES container recreation rather than needing a reinstall each rebuild.
set -euo pipefail

CONTAINER="${CONTAINER:-igvfagent-app}"
VENV="${BEAN_VENV:-/workspace/opt/bean-venv}"
CHECK_ONLY=0
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

say() { printf '\n=== %s ===\n' "$1"; }
dex()  { docker exec "$CONTAINER" sh -c "$1"; }
dexr() { docker exec -u root "$CONTAINER" sh -c "$1"; }

say "current state"
if dex "command -v bean >/dev/null 2>&1"; then
    echo "  bean on PATH: $(dex 'command -v bean')"
    dex "bean --help 2>&1 | head -3" | sed 's/^/  /'
    dex "$VENV/bin/python -c 'import numpy;print(\"  bean venv numpy:\", numpy.__version__)'" 2>/dev/null || true
    echo "  already installed — nothing to do."
    exit 0
fi
echo "  bean: not installed"
dex "/opt/venv/bin/python -c 'import numpy;print(\"  app numpy:\", numpy.__version__)'"

if [ "$CHECK_ONLY" = 1 ]; then
    printf '\n--check complete. Re-run without --check to install.\n'
    exit 0
fi

say "1. build toolchain (apt, sandboxed as root, cache redirected)"
dexr 'mkdir -p /tmp/apt/archives/partial /tmp/apt/lists/partial'
dexr 'apt-get -o APT::Sandbox::User=root \
        -o Dir::Cache::archives=/tmp/apt/archives \
        -o Dir::State::lists=/tmp/apt/lists update -qq' || true
dexr 'apt-get -o APT::Sandbox::User=root \
        -o Dir::Cache::archives=/tmp/apt/archives \
        -o Dir::State::lists=/tmp/apt/lists \
        install -y --no-install-recommends build-essential >/dev/null'
dex 'command -v gcc && gcc --version | head -1' | sed 's/^/  /'

say "2. isolated venv with NumPy 1.x (pip as the OWNER, never root)"
dex "mkdir -p \$(dirname $VENV) && python3 -m venv $VENV"
dex "$VENV/bin/pip install -q --upgrade pip 'cython<3' 'numpy<2' setuptools wheel"
dex "$VENV/bin/python -c 'import numpy;print(\"  build numpy:\", numpy.__version__)'"

say "3. crispr-bean, no build isolation so the Cython pin is honoured"
dex "$VENV/bin/pip install --no-build-isolation crispr-bean" | tail -2 | sed 's/^/  /'

say "4. pin NumPy back down (a dependency re-installs 2.x over the pin)"
dex "$VENV/bin/pip install -q 'numpy<2' 2>&1 | tail -2" | sed 's/^/  /' || true
dex "$VENV/bin/python -c 'import numpy;print(\"  runtime numpy:\", numpy.__version__)'"

say "5. expose on PATH and verify"
dexr "ln -sf $VENV/bin/bean /usr/local/bin/bean"
dex 'command -v bean' | sed 's/^/  /'
for sub in count run qc; do
    if dex "bean $sub --help >/dev/null 2>&1"; then
        echo "  bean $sub: OK"
    else
        echo "  bean $sub: FAILED" >&2
        exit 1
    fi
done
dex "$VENV/bin/python -c 'import bean, anndata, pyro, torch; print(\"  bean, anndata, pyro, torch import OK\")'"

cat <<'DONE'

DONE. Two caveats, both real:

  zarr wants numpy>=2 and has 1.x here. bean/anndata/pyro/torch and
  `bean count|run|qc --help` were all verified working, but a zarr code path
  reached during a full `bean run` on real data has NOT been exercised.

  This lives in the container plus a venv on the data volume. The venv
  survives recreation; the gcc install and the /usr/local/bin/bean symlink do
  NOT. Re-run this script after a recreate, or build the image with
  IGVF_INSTALL_CRISPR_BEAN=1 and stop needing it.
DONE
