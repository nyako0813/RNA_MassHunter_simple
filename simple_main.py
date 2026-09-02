"""CLI entrypoint for the simple RNA_MassHunter pipeline.

New file (仕様書 §6, §8.4). Independent of the original repository's
main.py — a thin argparse wrapper around `rna_masshunter.simple_pipeline.run()`.
"""
from __future__ import annotations

import argparse
import sys

from rna_masshunter import simple_pipeline


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="RNA_MassHunter simple pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
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
