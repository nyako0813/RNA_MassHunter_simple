"""CLI entrypoint for the simple RNA_MassHunter pipeline.

New file (仕様書 §6, §8.4), originally named `simple_main.py` and renamed to
`main.py` so the entry script is invoked the same way as in the other
projects. Independent of the main RNA_MassHunter repository's much larger
`main.py` — a thin argparse wrapper around `rna_masshunter.simple_pipeline.run()`.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rna_masshunter import simple_pipeline
from rna_masshunter.trna_library import load_trna_library

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = REPO_ROOT / "config.yaml"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="RNA_MassHunter simple pipeline")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH), help="Path to config.yaml (default: config.yaml in the repository root)")
    parser.add_argument(
        "--list-trna", action="store_true",
        help="List available tRNA types from data/trna_library.yaml and exit.",
    )
    return parser.parse_args(argv)


def _list_trna(root: Path = REPO_ROOT) -> int:
    library = load_trna_library(root / "data" / "trna_library.yaml")
    if not library:
        print("data/trna_library.yaml is empty or missing.", file=sys.stderr)
        return 1
    for entry_id in sorted(library):
        entry = library[entry_id]
        print(f"{entry_id}\t{entry.get('amino_acid')}\t{entry.get('anticodon')}\t{entry.get('length')}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.list_trna:
        return _list_trna()

    result = simple_pipeline.run(args.config)

    for warning in result.get("warnings", []):
        print(f"[{warning['Level']}] {warning['Source']}: {warning['Message']}", file=sys.stderr)

    output_path = result.get("output_path")
    if output_path:
        print(f"Wrote {output_path}")
    else:
        print("Excel output was disabled or skipped; no report was written.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
