"""Offline checks for the extension review queue (extension_review.py)."""
import subprocess
import sys
from pathlib import Path

sys.exit(subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "extension_review.py"), "selftest"]).returncode)
