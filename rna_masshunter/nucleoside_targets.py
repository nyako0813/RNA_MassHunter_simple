"""Known-nucleoside mass targets for Nuclease P1 complete-digestion mode.

New module (仕様書 §24, Phase 10). Nuclease P1 fully hydrolyzes RNA down to
individual ribonucleoside 5'-monophosphates (or, after alkaline
phosphatase, free nucleosides) rather than leaving oligomer fragments —
so this mode compares observed MS1 peaks directly against a fixed universe
of *known* nucleoside masses (the 4 standard bases plus every modified
nucleoside catalogued in data/modifications.yaml), instead of the
fragment/ΔDa-driven flow the rest of the pipeline uses for oligomer
digestion (RNase A/T1/T2 etc.).

Deliberately scoped to known-nucleoside mass matching only — no
Formula Candidate search or unknown-modification inference is attempted
here (仕様書 §24: unlike mass_comparison.py's ΔDa exploration, P1 mode
only ever asks "does this peak match one of ~100+ known nucleosides",
never "what could explain this residual mass difference").
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import DEFAULT_PHOSPHATE_MASS
from rna_masshunter.models import Modification
from rna_masshunter.warnings_manager import add_warning

# Standard ribonucleoside elemental formulas (free nucleoside — base + ribose,
# no phosphate). Monoisotopic masses computed via ElementalComposition and
# cross-checked against published values in tests/test_nucleoside_targets.py:
#   A 267.0968, C 243.0855, G 283.0917, U 244.0695 Da.
_STANDARD_NUCLEOSIDE_FORMULAS: dict[str, dict[str, int]] = {
    "A": {"C": 10, "H": 13, "N": 5, "O": 4},
    "C": {"C": 9, "H": 13, "N": 3, "O": 5},
    "G": {"C": 10, "H": 13, "N": 5, "O": 5},
    "U": {"C": 9, "H": 12, "N": 2, "O": 6},
}
_STANDARD_NUCLEOSIDE_NAMES: dict[str, str] = {
    "A": "Adenosine", "C": "Cytidine", "G": "Guanosine", "U": "Uridine",
}


@dataclass(frozen=True)
class NucleosideTarget:
    label: str          # e.g. "A", "m1A"
    name: str            # display name, e.g. "Adenosine", "1-methyladenosine"
    category: str        # "standard" | "modified"
    nucleoside_mass: float  # free nucleoside monoisotopic mass (no phosphate)


def build_standard_nucleoside_targets() -> list[NucleosideTarget]:
    """The 4 unmodified ribonucleosides, computed from their elemental
    formulas rather than read from data/base_masses.yaml — that file's
    `rna_residue_masses` are residue masses *within an oligomer chain*
    (already net of the water lost to phosphodiester bond formation), which
    is not the same quantity as a fully-hydrolyzed free nucleoside's mass."""
    return [
        NucleosideTarget(
            label=base, name=_STANDARD_NUCLEOSIDE_NAMES[base], category="standard",
            nucleoside_mass=ElementalComposition(formula).exact_mass,
        )
        for base, formula in _STANDARD_NUCLEOSIDE_FORMULAS.items()
    ]


def _modification_display_name(modification: Modification) -> str:
    raw = modification.raw or {}
    return str(raw.get("name") or modification.symbol or modification.id)


def build_modified_nucleoside_targets(
    modifications: list[Modification],
    warnings: list[dict[str, Any]] | None = None,
) -> list[NucleosideTarget]:
    """One target per data/modifications.yaml entry, using its
    `modified_nucleoside_mass_mono` field directly (already the free
    modified-nucleoside monoisotopic mass — no formula computation needed).
    Entries missing that field are skipped with a warning rather than
    guessed at."""
    targets: list[NucleosideTarget] = []
    for modification in modifications:
        raw = modification.raw or {}
        mass = raw.get("modified_nucleoside_mass_mono")
        if mass is None:
            if warnings is not None:
                add_warning(warnings, "INFO", "nucleoside_targets", "Modification has no modified_nucleoside_mass_mono; skipped for P1 mode.", modification.id or modification.symbol)
            continue
        try:
            mass = float(mass)
        except (TypeError, ValueError):
            if warnings is not None:
                add_warning(warnings, "WARNING", "nucleoside_targets", "modified_nucleoside_mass_mono is not numeric; skipped for P1 mode.", modification.id or modification.symbol)
            continue
        targets.append(NucleosideTarget(
            label=str(modification.symbol or modification.id),
            name=_modification_display_name(modification),
            category="modified",
            nucleoside_mass=mass,
        ))
    return targets


def build_nucleoside_target_universe(
    modifications: list[Modification],
    warnings: list[dict[str, Any]] | None = None,
) -> list[NucleosideTarget]:
    """Standard bases + modified nucleosides, in that order."""
    return build_standard_nucleoside_targets() + build_modified_nucleoside_targets(modifications, warnings)


def theoretical_mass_for_target(
    target: NucleosideTarget,
    dephosphorylated: bool,
    base_masses: dict[str, Any],
) -> float:
    """仕様書 §24: 既存の `config.alkaline_phosphatase.enabled` をそのまま
    流用する（P1モード専用の新規configキーは増やさない）。
    dephosphorylated=True（AP有効）なら遊離ヌクレオシド質量そのもの、
    False（AP無効）なら5'-一リン酸化分（+phosphate）を加える。
    `base_masses.yaml`の`constants.phosphate`を優先し、無ければ
    `masses.DEFAULT_PHOSPHATE_MASS`にフォールバックする（他のモジュールと
    同じ既存の値取得パターンに揃える）。"""
    if dephosphorylated:
        return target.nucleoside_mass
    phosphate = float((base_masses or {}).get("constants", {}).get("phosphate", DEFAULT_PHOSPHATE_MASS))
    return target.nucleoside_mass + phosphate
