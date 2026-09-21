"""Regenerate the `conserved_modifications` of data/trna_library.yaml
(claude_code/conserved_modifications_confidence_tier_spec.md).

Not part of the pipeline — a reproducible one-off generator. Running it is
idempotent: hand-curated single-label entries (archaeosine G+ at 15,
agmatidine C+ at the Ile(CAU) wobble) are kept and given their confidence
tier; every `modification_candidates` entry is dropped and rebuilt from the
position/base rules below. The file header comment is left untouched.

    python tools/add_conserved_modifications.py            # rewrite the yaml
    python tools/add_conserved_modifications.py --check    # exit 1 if it would change

Standard (Sprinzl) positions are estimated from the library's own running
numbering: 34 = `wobble_position`, 37 = wobble + 3, and 55 / 58 are counted
back from the discriminator base (standard 73). That assumes a standard
T-arm; when it doesn't hold the base at the estimated position simply fails
the rule's base check and nothing is added (safe direction).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

LIBRARY_PATH = Path(__file__).resolve().parent.parent / "data" / "trna_library.yaml"

CONFIRMED = "confirmed"
HIGH_PROBABILITY = "high_probability"

ARCHAEOSINE_NOTE = (
    "position 15 is G in nearly all elongator tRNAs; archaeosine is the default hypothesis but not guaranteed "
    "(observed mass loss in intact-tRNA MS suggests occasional absence)"
)

# position/base/candidates/note per rule. Rules whose candidate sets come from the main
# RNA_MassHunter repository's rule_sets/ are named in the note. `note_extra` optionally
# appends a sentence to the note, only for tRNAs whose anticodon 3rd base (standard
# position 36) equals `when_anticodon_3rd_base` — a claim that is only true for those.
# "methylation" (generic class, no defined mass) from archaea_A58_methylation is
# deliberately omitted: m1A / m6A already cover the same +CH2 mass.
RULES: list[dict[str, Any]] = [
    {
        "position": 55, "base": "U", "candidates": ["Y"],
        "note": "standard position 55 (discriminator - 18); base is U; pseudouridine is near-universal "
                "(archaea_general.yaml universal_U55_pseudouridine). Note: pseudouridine is isobaric with U, "
                "so a mass match here cannot distinguish it from unmodified U.",
    },
    {
        "position": 34, "base": "U",
        "candidates": ["cnm5U", "cmnm5U", "mnm5U", "mnm5s2U", "s2U", "s4U", "mcm5U", "mcm5s2U", "ncm5U", "ncm5Um"],
        "note": "standard position 34 (wobble); base is U; candidate set from "
                "rule_sets/methanosarcina_acetivorans.yaml ma_U34_main_target.",
    },
    {
        "position": 37, "base": "A", "candidates": ["t6A", "ms2t6A", "hn6A", "ms2hn6A"],
        "note": "standard position 37 (wobble + 3); base is A; candidate set from "
                "rule_sets/methanosarcina_acetivorans.yaml ma_A37_t6A.",
        "note_extra": {
            "when_anticodon_3rd_base": "U",
            "text": "Consistent with the anticodon 3rd base = U correlation (all 16 such tRNAs have A37).",
        },
    },
    {
        "position": 58, "base": "A", "candidates": ["m1A", "m6A"],
        "note": "standard position 58 (discriminator - 15); base is A; candidate set from "
                "rule_sets/archaea_general.yaml archaea_A58_methylation (generic 'methylation' omitted: "
                "m1A/m6A cover the same +CH2 mass).",
    },
    {
        "position": 37, "base": "G", "candidates": ["imG-14", "imG", "imG2"], "only_for_id_prefix": "tRNA-Phe-GAA-",
        "note": "standard position 37 (wobble + 3); base is G; wyosine-pathway modifications are near-universal at "
                "position 37 of tRNA-Phe(GAA) in eukaryotes and archaea (G37 is their substrate). Additional support: "
                "this mass group (321.11 / 335.12 Da) was detected at multiple peaks with sufficient intensity in "
                "RNA_MassHunter_simple_report.xlsx.",
    },
]


def discriminator_position(seq: str) -> int:
    """Running number (1-based) of standard position 73 (the discriminator base)."""
    if seq.endswith("CCA"):
        return len(seq) - 3
    return len(seq)


def canonical_positions(entry: dict[str, Any]) -> dict[int, int]:
    """Standard position -> this entry's running number. (Standard 15 equals
    running 15 in this library — verified for all 58 entries — so it is not
    listed.)"""
    wobble = entry["wobble_position"]
    d = discriminator_position(entry["sequence"])
    return {34: wobble, 37: wobble + 3, 55: d - 18, 58: d - 15}


def base_at(entry: dict[str, Any], running_position: int) -> str | None:
    seq = entry["sequence"]
    return seq[running_position - 1] if 1 <= running_position <= len(seq) else None


def rule_entries(entry: dict[str, Any]) -> list[dict[str, Any]]:
    """The `modification_candidates` entries the rules assign to one tRNA."""
    positions = canonical_positions(entry)
    found = []
    for rule in RULES:
        prefix = rule.get("only_for_id_prefix")
        if prefix and not entry["id"].startswith(prefix):
            continue
        position = positions[rule["position"]]
        if base_at(entry, position) != rule["base"]:
            continue
        note = rule["note"]
        extra = rule.get("note_extra")
        if extra and entry["anticodon"][2] == extra["when_anticodon_3rd_base"]:
            note = f"{note} {extra['text']}"
        found.append({
            "position": position,
            "modification_candidates": list(rule["candidates"]),
            "confidence": HIGH_PROBABILITY,
            "note": note,
        })
    return found


def _curated_entry(item: dict[str, Any]) -> dict[str, Any]:
    """Tier a hand-curated single-label entry, keeping its own note unless it is the archaeosine default."""
    label = item.get("modification")
    if label == "G+":
        return {"position": item["position"], "modification": "G+", "confidence": HIGH_PROBABILITY, "note": ARCHAEOSINE_NOTE}
    if label == "C+":
        return {"position": item["position"], "modification": "C+", "confidence": CONFIRMED, "note": item["note"]}
    raise ValueError(f"unexpected curated conserved modification {item!r}; add its tier to _curated_entry()")


def rebuild(entry: dict[str, Any]) -> list[dict[str, Any]]:
    curated = [_curated_entry(item) for item in entry.get("conserved_modifications", []) if "modification_candidates" not in item]
    return sorted(curated + rule_entries(entry), key=lambda item: item["position"])  # stable: curated before rule at a tie


def render(text: str) -> str:
    split = text.index("trna_library:")
    header, body = text[:split], text[split:]
    data = yaml.safe_load(body)
    for entry in data["trna_library"]:
        entry["conserved_modifications"] = rebuild(entry)
    return header + yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=10000)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="exit 1 instead of writing if the yaml is out of date")
    args = parser.parse_args(argv)
    current = LIBRARY_PATH.read_text(encoding="utf-8")
    updated = render(current)
    if args.check:
        return 0 if updated == current else 1
    LIBRARY_PATH.write_text(updated, encoding="utf-8")
    print(f"wrote {LIBRARY_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
