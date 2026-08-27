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


# ---------------------------- runtime stage --------------------------------
FROM python:${PYTHON_VERSION}-slim-bookworm AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    IGVF_PROJECT_ROOT=/workspace \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true

# Minimal runtime libs:
#   libgomp1   numpy/scipy/sklearn-style scientific stack
#   curl       healthcheck against the Streamlit HTTP probe
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
        curl \
        ca-certificates \
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
