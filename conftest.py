"""Ensures the repo root (and therefore the `rna_masshunter` package) is
importable regardless of how pytest is invoked."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
