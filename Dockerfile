# syntax=docker/dockerfile:1.7
#
# IGVFagent — multi-stage container image.
#
#   builder  : installs the package with the [all] extras (analysis +
#              ui + llm) into an isolated /opt/venv. Build deps live
#              here only and never reach the runtime layer.
#   runtime  : lean python:3.11-slim layer that copies /opt/venv and
#              runs as a non-root `igvf` user out of /workspace.
#              IGVF_PROJECT_ROOT pins all skills to /workspace so
#              mounted ./Data and ./Docs persist analyses across runs.
#
# Build:
#   docker build -t igvfagent:latest .
#
# Run the UI:
#   docker run --rm -p 8501:8501 \
#     -v "$PWD/Data:/workspace/Data" -v "$PWD/Docs:/workspace/Docs" \
#     igvfagent:latest
#
# Run a one-shot skill:
#   docker run --rm igvfagent:latest kg gene APOE --depth 1 --limit 5

ARG PYTHON_VERSION=3.11

# ---------------------------- builder stage --------------------------------
FROM python:${PYTHON_VERSION}-slim-bookworm AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Build deps for scientific wheels that occasionally need a compiler.
#
# libcurl4-openssl-dev + zlib1g-dev are required by the [hic] extra:
# hic-straw ships no wheel, so pip compiles src/straw.cpp, which does
# `#include <curl/curl.h>` (remote .hic reading over HTTP) and links zlib
# for block decompression. Without the headers `pip install '.[all]'`
# dies with "fatal error: curl/curl.h: No such file or directory".
# Builder-stage only — the runtime layer needs just the shared libs.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        curl \
        ca-certificates \
        libcurl4-openssl-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY pyproject.toml LICENSE README.md ./
COPY Scripts ./Scripts

# Isolated venv -> easy to copy into the slim runtime layer.
RUN python -m venv /opt/venv \
 && /opt/venv/bin/pip install --upgrade pip \
 && /opt/venv/bin/pip install '.[all]'

# crispr-bean (BEAN) for base-editing screens. OPT-IN, because it pulls in
# torch + pyro for its Bayesian variant/tiling model and roughly doubles the
# image; most deployments do not analyse base-editing screens.
#
#   docker compose build --build-arg INSTALL_CRISPR_BEAN=1
#
# It cannot be installed into the RUNTIME layer at all: that layer has no
# compiler, and BEAN's bean/mapping/CRISPResso2Align.pyx needs Cython and a
# C toolchain, so `pip install crispr-bean` there dies with
# "CompileError: bean/mapping/CRISPResso2Align.pyx". Here in the builder
# stage build-essential is already present for the [hic] extra.
#
# BEAN is AGPL-3.0 and IGVFagent is Apache-2.0. This INSTALLS it as a
# separate program invoked as a subprocess -- no linking, no vendoring, no
# licence propagation. `igvfagent bean` implements BEAN's guide-assignment
# method independently; the installed `bean` binary is what provides the
# Bayesian model that is deliberately not reimplemented.
# INTO ITS OWN VENV, not /opt/venv. BEAN cannot share the app's interpreter
# environment: bean/mapping/CRISPResso2Align.pyx uses np.int_t / np.uint_t and
# the code uses np.Inf, all removed in NumPy 2.0, so BEAN needs numpy<2 --
# while scanpy / numba / anndata in /opt/venv are resolved against numpy 2.x.
# Pinning /opt/venv down to numpy<2 to satisfy BEAN would break the stack the
# rest of the app runs on. Since `bean` is invoked as a SUBPROCESS
# (base_editing_screen.bean_available uses shutil.which), a separate venv on
# PATH is all that is needed and nothing has to agree about numpy.
#
# Each pin below is load-bearing, and all three were established the hard way
# in Deploy/Dockerfile.bean029:
#   numpy<2   the .pyx types above
#   cython<3  CRISPResso2Align.pyx does not compile under Cython 3 --
#             "'uint_t' is not a type identifier", which is exactly how the
#             unpinned version of this step failed
#   --no-build-isolation  or pip builds in a fresh overlay carrying its own
#             Cython 3 and BOTH pins are silently ignored
#   torch from the CPU index  BEFORE crispr-bean, or its resolver pulls the
#             default CUDA wheels: nvidia-cublas, cudnn, nccl and the rest,
#             several GB onto a VM with no GPU. Measured: 1.8 GB with the CPU
#             wheel against a multi-GB tree without it.
#   numpy<2 AGAIN, and zarr<3, AFTER crispr-bean  its dependency tree
#             reinstalls numpy 2.x over the pin (measured: 1.26.4 -> 2.4.6),
#             and `bean --help` then dies on `np.Inf` at import. zarr 3
#             requires numpy>=2, so it is pinned back with it.
#
# Verified end to end before this was committed: `bean run --help` answers and
# bean.read_h5ad loads the paper's own deposit as a ReporterScreen (3451, 40)
# with X_bcmatch present.
#
# Verified with `bean --help`, NOT `bean --version`: BEAN has no --version
# flag and answers "error: unrecognized arguments: --version" with exit 2, so
# the previous line here could never have passed even had the build worked.
# chromap — the scATAC half of the IGVF uniform pipeline.
#
# IGVF/atomic-workflows (MIT) has exactly TWO modules: igvf-kallisto-bustools
# and igvf-chromap. kb-python already provides the first, so this binary is
# what completes IGVFagent's coverage of the official pipeline's tool set:
# kb for RNA, chromap for ATAC.
#
# Built from the release tarball because there is no apt package, no PyPI
# package, and the GitHub release ships source only. It is a small, cheap
# build -- 1.5 MB binary, and the toolchain is already here for the [hic]
# extra -- which is why this is worth doing where STAR (a ~30 GB-RAM index
# build) is not. Verified before committing: chromap --version reports
# 0.3.2-r518.
# Unconditional, deliberately. A build ARG would need the runtime COPY to be
# conditional too, and the usual workaround -- stubbing the path so COPY
# succeeds -- puts a file named `chromap` on PATH that is not chromap. That is
# the exact false-positive that `bean --version` and kb_available already had
# to be hardened against: a tool that looks installed and fails on use is
# worse than one that is honestly absent.
RUN cd /tmp \
 && curl -sL --max-time 300 \
        https://github.com/haowenz/chromap/archive/refs/tags/v0.3.2.tar.gz \
        -o chromap.tgz \
 && tar xzf chromap.tgz \
 && cd chromap-0.3.2 \
 && make -j"$(nproc)" \
 && install -m 0755 chromap /usr/local/bin/chromap \
 && cd /tmp && rm -rf chromap-0.3.2 chromap.tgz \
 && /usr/local/bin/chromap --version

# CRISPResso2 is installed into this SAME venv, not the app venv. It needs the
# same numpy<2 (BEAN already vendors its CRISPResso2Align.pyx), and BEAN reads
# the reporter allele THROUGH CRISPResso2 alignment -- bystander edit
# deconvolution is unavailable without it, which is a gap base_editing_screen
# reports today. It is not on PyPI, so it installs from the repo, and
# --no-build-isolation keeps the cython<3 pin above in force. Verified in a
# throwaway container before this was committed: `bean --help` and
# `CRISPResso --help` both answer, with numpy 1.26.4 and torch 2.4.1+cpu.
ARG INSTALL_CRISPR_BEAN=0
RUN if [ "$INSTALL_CRISPR_BEAN" = "1" ]; then \
        python -m venv /opt/bean-venv \
     && /opt/bean-venv/bin/pip install --no-cache-dir --upgrade pip \
     && /opt/bean-venv/bin/pip install --no-cache-dir "numpy<2" "cython<3" \
     && /opt/bean-venv/bin/pip install --no-cache-dir \
            --index-url https://download.pytorch.org/whl/cpu "torch==2.4.1" \
     && /opt/bean-venv/bin/pip install --no-cache-dir --no-build-isolation \
            crispr-bean \
     && /opt/bean-venv/bin/pip install --no-cache-dir --no-build-isolation \
            git+https://github.com/pinellolab/CRISPResso2.git \
     && /opt/bean-venv/bin/pip install --no-cache-dir "numpy<2" "zarr<3" \
     && /opt/bean-venv/bin/bean --help > /dev/null \
     && /opt/bean-venv/bin/bean run --help > /dev/null \
     && /opt/bean-venv/bin/CRISPResso --help > /dev/null \
     && /opt/bean-venv/bin/python -c "import numpy, torch, bean; assert numpy.__version__.startswith('1.'), numpy.__version__; print('bean venv: numpy', numpy.__version__, '| torch', torch.__version__)"; \
    else \
        echo "crispr-bean NOT installed (INSTALL_CRISPR_BEAN=0)."; \
        mkdir -p /opt/bean-venv; \
    fi


# ---------------------------- runtime stage --------------------------------
FROM python:${PYTHON_VERSION}-slim-bookworm AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:/opt/bean-venv/bin:$PATH" \
    IGVF_PROJECT_ROOT=/workspace \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true

# Minimal runtime libs:
#   libgomp1   numpy/scipy/sklearn-style scientific stack
#   curl       healthcheck against the Streamlit HTTP probe
#   zlib1g     bustools links libz.so.1 (kb-python's bundled binary).
#              The slim base already carries it via CPython's zlib module,
#              so this is belt-and-braces against a base-image change
#              silently breaking `kb count` -- the aligner needs no other
#              system library. Verified by reading DT_NEEDED off the
#              shipped ELFs: kallisto wants only glibc + libstdc++, and
#              notably NOT libhdf5 (the darwin build does, which is why it
#              fails on a mac without homebrew hdf5 but works here).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
        curl \
        ca-certificates \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Optional: the Claude Code CLI, which backs the "External orchestrator"
# choice in the UI (_llm._chat_claude_cli shells out to `claude --print`).
# Without it that option can only ever report "not installed", because the
# backend needs the *program* on PATH — an ANTHROPIC_API_KEY alone does not
# put it there. Headless `--print` authenticates from ANTHROPIC_API_KEY, so
# no extra credential is required. Build with
# --build-arg INSTALL_CLAUDE_CLI=0 to skip it and keep the image lean.
ARG INSTALL_CLAUDE_CLI=1
RUN if [ "$INSTALL_CLAUDE_CLI" = "1" ]; then \
        apt-get update && apt-get install -y --no-install-recommends \
            nodejs npm \
        && npm install -g --no-fund --no-audit @anthropic-ai/claude-code \
        && npm cache clean --force \
        && rm -rf /var/lib/apt/lists/* ; \
    fi

# Non-root user with a writable home + workspace. Claude Code writes its
# config under $HOME, so the home directory is load-bearing, not cosmetic.
RUN useradd --create-home --shell /bin/bash --uid 1000 igvf \
 && mkdir -p /workspace/Data /workspace/Docs \
 && chown -R igvf:igvf /workspace

COPY --from=builder --chown=igvf:igvf /opt/venv /opt/venv
# Empty when INSTALL_CRISPR_BEAN=0, so this COPY is unconditional and
# costs nothing in the default build.
COPY --from=builder --chown=igvf:igvf /opt/bean-venv /opt/bean-venv
# chromap is a single 1.5 MB binary, always built above.
COPY --from=builder /usr/local/bin/chromap /usr/local/bin/chromap

USER igvf
WORKDIR /workspace

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["igvfagent"]
# Default command: launch the browser UI bound to all interfaces. Override
# at runtime to drive any other skill, e.g.
#   docker run --rm igvfagent:latest kg gene APOE
CMD ["ui", "--host", "0.0.0.0", "--no-browser"]
