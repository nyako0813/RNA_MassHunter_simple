"""Modification-hypothesis theoretical-mass check (hypothesis_mass_check_spec.md).

Given hypotheses of the form "this position carries this modification" (or
a dinucleotide / base+element-adduct variant), compute each hypothesis'
theoretical neutral mass and report whether raw MS1 peaks exist at that
mass. Deliberately a presence check only (spec §1, §7): nothing here infers
which modification an observed mass "is", and no candidates are generated
for masses that don't match.

Two hypothesis sources feed the same matching logic (spec §3):
- `conserved_modifications` on the tRNA selected via `sequence.trna_type`
  (data/trna_library.yaml) -> source "trna_library_default", each entry
  carrying a `confidence` tier and either one `modification` or a
  `modification_candidates` list (expanded into independent hypotheses)
- `config.hypothesis_check.targets` -> source "config_manual"

Three target forms (spec §3.2, §4):
- `label`: one catalog nucleoside, mass taken as-is (spec §4.1) from
  nucleoside_targets.build_nucleoside_target_universe (the 4 standard
  bases + data/modifications.yaml).
- `components` (2 labels) + `linkage`: dinucleotide (spec §4.2).
- `base` + `add_elements`: catalog nucleoside plus simple monoisotopic
  element additions (spec §4.3).

Dinucleotide mass reuses the elemental-composition model of the main
repository's `p1_sap_dinucleotide_candidates.py` (LINKAGE_COMPOSITIONS /
CONDENSATION_ADJUSTMENT, copied here because that repository isn't an
importable dependency of this one) instead of a second, independently
derived formula: N1 + N2 + linkage - H2O, where linkage is HPO3 (normal
phosphodiester) or HPO2S (phosphorothioate). That is numerically identical
to the spec's `N1 + N2 + HPO3 - H2O (+ S - O for phosphorothioate)`.
5'/3' order is not distinguished (mass can't tell them apart, spec §7).

Peaks are matched on neutral mass (spec §5): for each charge 1..max_charge
the observed neutral mass is `neutral_mass_from_mz(mz, charge, polarity)`
and must be within `mass_tolerance_ppm` of the theoretical mass.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass
from typing import Any

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES, mz_from_neutral_mass, neutral_mass_from_mz
from rna_masshunter.mass_shift_ms1_search import build_sorted_peak_index
from rna_masshunter.models import Modification, Peak, RunConfig
from rna_masshunter.nucleoside_targets import build_nucleoside_target_universe
from rna_masshunter.observed_mass import assign_peak_ids
from rna_masshunter.warnings_manager import add_warning

SOURCE_TRNA_LIBRARY_DEFAULT = "trna_library_default"
SOURCE_CONFIG_MANUAL = "config_manual"

# Same compositions as the main repository's p1_sap_dinucleotide_candidates.py.
_LINKAGE_COMPOSITIONS = {
    "phosphodiester": ElementalComposition({"H": 1, "O": 3, "P": 1}),
    "phosphorothioate": ElementalComposition({"H": 1, "O": 2, "P": 1, "S": 1}),
}
_CONDENSATION_ADJUSTMENT = ElementalComposition.delta({"H": -2, "O": -1})


@dataclass
class Hypothesis:
    name: str
    source: str
    description: str
    theoretical_mass: float
    # "confirmed" / "high_probability" for trna_library defaults; "" for
    # config.hypothesis_check.targets (ad-hoc hypotheses carry no tier).
    confidence: str = ""


@dataclass
class HypothesisCheckRow:
    hypothesis_name: str
    source: str
    formula_description: str
    theoretical_mass: float
    match_found: bool
    confidence: str = ""
    peak_id: str | None = None
    charge: int | None = None
    observed_mass: float | None = None
    delta_da: float | None = None
    delta_ppm: float | None = None
    intensity: float | None = None
    rt: float | None = None
    scan_count: int = 1
    rt_range: tuple[float, float] | None = None


# --- theoretical masses --------------------------------------------------------

def dinucleotide_mass(mass_n1: float, mass_n2: float, linkage: str = "phosphodiester") -> float:
    """spec §4.2: N1 + N2 + linkage - H2O, from the two free-nucleoside masses."""
    key = str(linkage).strip().lower()
    if key not in _LINKAGE_COMPOSITIONS:
        raise ValueError(f"linkage must be one of {sorted(_LINKAGE_COMPOSITIONS)}, got {linkage!r}")
    return mass_n1 + mass_n2 + _LINKAGE_COMPOSITIONS[key].exact_mass + _CONDENSATION_ADJUSTMENT.exact_mass


def base_plus_elements_mass(base_mass: float, add_elements: dict[str, int]) -> float:
    """spec §4.3: catalog mass plus the monoisotopic mass of each added element."""
    total = float(base_mass)
    for element, count in add_elements.items():
        if element not in MONOISOTOPIC_ATOMIC_MASSES:
            raise ValueError(f"Unsupported element in add_elements: {element}")
        if isinstance(count, bool) or not isinstance(count, int):
            raise ValueError(f"add_elements[{element}] must be an integer, got {count!r}")
        total += MONOISOTOPIC_ATOMIC_MASSES[element] * count
    return total


def build_label_mass_catalog(
    modifications: list[Modification],
    warnings: list[dict[str, Any]] | None = None,
) -> dict[str, float]:
    """label -> free-nucleoside mass (spec §4.1): reuses nucleoside_targets'
    universe rather than re-reading data/modifications.yaml separately."""
    return {t.label: t.nucleoside_mass for t in build_nucleoside_target_universe(modifications, warnings)}


def _lookup(catalog: dict[str, float], label: Any, context: str) -> float:
    key = str(label)
    if key not in catalog:
        raise ValueError(f"{context}: label {key!r} was not found in the nucleoside catalog (4 standard bases + data/modifications.yaml)")
    return catalog[key]


def _format_elements(add_elements: dict[str, int]) -> str:
    return " ".join(f"{'+' if count >= 0 else '-'}{element}{abs(count) if abs(count) != 1 else ''}" for element, count in add_elements.items())


def hypothesis_from_target(target: dict[str, Any], catalog: dict[str, float], source: str = SOURCE_CONFIG_MANUAL) -> Hypothesis:
    """Build one Hypothesis from a `hypothesis_check.targets` entry. Raises
    ValueError on a malformed entry (unknown label/element/linkage, or not
    exactly one of the three target forms)."""
    if not isinstance(target, dict):
        raise ValueError(f"hypothesis_check target must be a mapping, got {target!r}")
    forms = [key for key in ("label", "components", "base") if target.get(key) not in (None, "", [])]
    if len(forms) != 1:
        raise ValueError(f"hypothesis_check target {target.get('name')!r} must specify exactly one of label / components / base (got {forms or 'none'})")
    form = forms[0]
    name = str(target.get("name") or "").strip()
    context = f"hypothesis_check target {name or target!r}"

    if form == "label":
        label = str(target["label"])
        mass = _lookup(catalog, label, context)
        return Hypothesis(name or label, source, f"{label} ({mass:.4f})", mass)

    if form == "components":
        components = target["components"]
        if not isinstance(components, (list, tuple)) or len(components) != 2:
            raise ValueError(f"{context}: components must be a list of exactly 2 labels")
        linkage = str(target.get("linkage") or "phosphodiester").strip().lower()
        mass = dinucleotide_mass(_lookup(catalog, components[0], context), _lookup(catalog, components[1], context), linkage)
        description = f"{components[0]} + {components[1]} dinucleotide, {linkage} linkage"
        return Hypothesis(name or f"{components[0]}-{components[1]} ({linkage})", source, description, mass)

    base = str(target["base"])
    add_elements = target.get("add_elements")
    if not isinstance(add_elements, dict) or not add_elements:
        raise ValueError(f"{context}: base requires a non-empty add_elements mapping")
    mass = base_plus_elements_mass(_lookup(catalog, base, context), add_elements)
    return Hypothesis(name or f"{base} {_format_elements(add_elements)}", source, f"{base} ({catalog[base]:.4f}) {_format_elements(add_elements)}", mass)


def hypotheses_from_config_targets(config: RunConfig, catalog: dict[str, float]) -> list[Hypothesis]:
    targets = (config.hypothesis_check or {}).get("targets") or []
    return [hypothesis_from_target(target, catalog, SOURCE_CONFIG_MANUAL) for target in targets]


_VALID_CONFIDENCE = {"confirmed", "high_probability"}


def hypotheses_from_conserved_modifications(
    trna_entry: dict[str, Any] | None,
    catalog: dict[str, float],
    modifications: list[Modification],
    warnings: list[dict[str, Any]] | None = None,
) -> list[Hypothesis]:
    """One hypothesis per label of each `conserved_modifications` entry of the
    selected tRNA (spec §3.1). An entry gives its label either as a single
    `modification` or as a `modification_candidates` list; every candidate
    becomes its own independent hypothesis sharing the entry's position and
    `confidence` (they are mutually exclusive in biology, but the check just
    tests each one). A label missing from the catalog is skipped with an
    ERROR warning rather than aborting the run (a hand-edited library typo
    shouldn't take the whole report down). A modification whose target_bases
    don't include the base at `position` gets a WARNING but is still checked
    — the entry is a hypothesis, not a guarantee."""
    if not trna_entry:
        return []
    sequence = str(trna_entry.get("sequence") or "")
    target_bases_by_label = {str(m.symbol or m.id): [str(b).upper() for b in (m.target_bases or [])] for m in modifications}
    hypotheses: list[Hypothesis] = []
    for item in trna_entry.get("conserved_modifications") or []:
        position = item.get("position")
        labels = ([item["modification"]] if item.get("modification") else []) + list(item.get("modification_candidates") or [])
        confidence = str(item.get("confidence") or "")
        if confidence and confidence not in _VALID_CONFIDENCE and warnings is not None:
            add_warning(warnings, "WARNING", "hypothesis_check", f"{trna_entry.get('id')}: position {position} has unknown confidence {confidence!r} (expected one of {sorted(_VALID_CONFIDENCE)}).")
        base_at_position = sequence[position - 1].upper() if isinstance(position, int) and 1 <= position <= len(sequence) else None
        for label in map(str, labels):
            if label not in catalog:
                if warnings is not None:
                    add_warning(warnings, "ERROR", "hypothesis_check", f"{trna_entry.get('id')}: conserved modification {label!r} at position {position} is not in the nucleoside catalog; skipped.")
                continue
            expected = target_bases_by_label.get(label)
            if warnings is not None and base_at_position and expected and base_at_position not in expected:
                add_warning(warnings, "WARNING", "hypothesis_check", f"{trna_entry.get('id')}: position {position} is {base_at_position} but {label} targets {expected}.")
            mass = catalog[label]
            hypotheses.append(Hypothesis(
                name=f"{trna_entry.get('id')} position {position} = {label}",
                source=SOURCE_TRNA_LIBRARY_DEFAULT,
                description=f"{label} ({mass:.4f})" + (f" — {item['note']}" if item.get("note") else ""),
                theoretical_mass=mass,
                confidence=confidence,
            ))
    return hypotheses


# --- matching ------------------------------------------------------------------

def check_hypotheses(hypotheses: list[Hypothesis], peaks: list[Peak], config: RunConfig) -> list[HypothesisCheckRow]:
    """spec §5/§6: one or more rows per hypothesis — one row per matching
    (peak, charge), or a single Match_Found=No row when nothing matches."""
    hc = config.hypothesis_check or {}
    tolerance_ppm = float(hc.get("mass_tolerance_ppm", 15) or 15)
    max_charge = int(hc.get("max_charge", 4) or 4)
    polarity = str((config.instrument or {}).get("polarity", "negative") or "negative").lower()

    peak_ids = assign_peak_ids(peaks)
    id_by_object = {id(peak): peak_ids[index] for index, peak in enumerate(peaks)}
    index = build_sorted_peak_index(peaks)

    rows: list[HypothesisCheckRow] = []
    for hypothesis in hypotheses:
        theoretical = hypothesis.theoretical_mass
        window_da = theoretical * tolerance_ppm / 1e6
        matches: list[HypothesisCheckRow] = []
        for charge in range(1, max_charge + 1):
            # A neutral-mass window of ±window_da is ±window_da/charge in m/z.
            center = mz_from_neutral_mass(theoretical, charge, polarity)
            half_width = window_da / charge
            lo = bisect.bisect_left(index.mzs, center - half_width)
            hi = bisect.bisect_right(index.mzs, center + half_width)
            for peak in index.peaks[lo:hi]:
                observed = neutral_mass_from_mz(peak.mz, charge, polarity)
                delta_da = observed - theoretical
                delta_ppm = delta_da / theoretical * 1e6
                if abs(delta_ppm) > tolerance_ppm:
                    continue
                matches.append(HypothesisCheckRow(
                    hypothesis_name=hypothesis.name, source=hypothesis.source,
                    formula_description=hypothesis.description, theoretical_mass=theoretical,
                    match_found=True, confidence=hypothesis.confidence, peak_id=id_by_object.get(id(peak)), charge=charge,
                    observed_mass=observed, delta_da=delta_da, delta_ppm=delta_ppm,
                    intensity=peak.intensity, rt=peak.rt,
                    scan_count=getattr(peak, "scan_count", 1), rt_range=getattr(peak, "rt_range", None),
                ))
        if matches:
            matches.sort(key=lambda row: (row.charge, abs(row.delta_ppm)))
            rows.extend(matches)
        else:
            rows.append(HypothesisCheckRow(
                hypothesis_name=hypothesis.name, source=hypothesis.source,
                formula_description=hypothesis.description, theoretical_mass=theoretical, match_found=False,
                confidence=hypothesis.confidence,
            ))
    return rows
