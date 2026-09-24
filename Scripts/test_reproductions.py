"""Offline checks for per-paper reproduction records (reproductions.py):
records from a job's evidence, the self-contained HTML, owner privacy and
publishing, attempts accumulating, and forum posts that stay previews until
confirmed."""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.exit(subprocess.run([sys.executable, str(HERE / "reproductions.py"), "selftest"]).returncode)
