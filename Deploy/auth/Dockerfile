# The authentication gate. Deliberately the smallest container in the stack:
# it decides who may drive an agent that runs analysis CLIs as subprocesses, so
# its attack surface should be a single Python file and nothing else.
#
# No `pip install` anywhere below. gate.py is stdlib-only, which means this
# image has no dependency tree to audit, to patch, or to be compromised
# through -- the property that matters most for the one component that cannot
# be allowed to fail open.
FROM python:3.12-alpine

RUN adduser -D -u 1000 gate
WORKDIR /app
COPY --chown=gate:gate Deploy/auth/gate.py /app/gate.py
# The gate serves the pages an unauthenticated visitor sees, so it needs the
# logo. Built from the repo root (see the compose build context) purely so
# these two files are reachable.
COPY --chown=gate:gate Docs/Figures/logo.png Docs/Figures/logo-mark.png /app/branding/

USER gate
EXPOSE 9000

# No HEALTHCHECK curl: the image has no curl and adding one would mean adding
# a package. Python is already here.
HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
    CMD python3 -c "import urllib.request,sys; \
sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:9000/_auth/health', timeout=3).status == 200 else 1)"

ENTRYPOINT ["python3", "-u", "/app/gate.py"]
