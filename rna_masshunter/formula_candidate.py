"""Elemental-composition ("Formula Candidate") enumeration for a ΔDa value.

New module (no vendored counterpart) implementing 仕様書
docs/design/RNA_MassHunter_再設計_実装仕様書.md §8.1, §12.

Deliberately does not import `rna_masshunter.modifications` — Formula
Candidate and Modification Candidate must stay fully independent code paths
(仕様書 §12.4, Claude Code rule #2).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from rna_masshunter.elemental_composition import ElementalComposition, _ELEMENT_ORDER
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES


@dataclass(frozen=True)
class FormulaCandidateResult:
    formula: str
    exact_mass: float
    mass_error_da: float
    mass_error_ppm: float | None
    atom_count: int
    elements: dict[str, int]


def calculate_delta_mass(observed_mass: float, theoretical_mass: float) -> float:
    """ΔDa = observed - theoretical."""
    return observed_mass - theoretical_mass


@lru_cache(maxsize=4096)
def _enumerate_cached(
    mass_difference_rounded: float,
    tolerance_da: float,
    max_total_atoms: int,
    element_limits: tuple[tuple[str, int], ...],
    elements: tuple[str, ...],
) -> tuple[FormulaCandidateResult, ...]:
    limits = dict(element_limits)
    results: list[FormulaCandidateResult] = []
    counts: dict[str, int] = {}

    def dfs(idx: int, remaining_atoms: int, running_mass: float) -> None:
        if idx == len(elements):
            error = running_mass - mass_difference_rounded
            if abs(error) <= tolerance_da:
                comp = ElementalComposition.delta(dict(counts))
                results.append(FormulaCandidateResult(
                    formula=comp.canonical_string(),
                    exact_mass=comp.exact_mass,
                    mass_error_da=error,
                    mass_error_ppm=None,
                    atom_count=sum(abs(v) for v in counts.values()),
                    elements=comp.to_dict(),
                ))
            return
        element = elements[idx]
        limit = min(limits.get(element, max_total_atoms), remaining_atoms)
        atomic_mass = MONOISOTOPIC_ATOMIC_MASSES[element]
        for count in range(-limit, limit + 1):
            if count == 0:
                counts.pop(element, None)
            else:
                counts[element] = count
            dfs(idx + 1, remaining_atoms - abs(count), running_mass + atomic_mass * count)
        counts.pop(element, None)

    dfs(0, max_total_atoms, 0.0)
    results.sort(key=lambda r: (abs(r.mass_error_da), r.formula))
    return tuple(results)


def enumerate_formula_candidates(
    mass_difference: float,
    tolerance_da: float,
    *,
    max_total_atoms: int = 7,
    element_limits: dict[str, int] | None = None,
    elements: tuple[str, ...] = _ELEMENT_ORDER,
    max_candidates: int | None = None,
) -> list[FormulaCandidateResult]:
    """Enumerate elemental compositions (signed atom counts, C/H/N/O/P/S/Se
    by default) that explain `mass_difference` (Da) within `tolerance_da`,
    subject to sum(|count|) <= max_total_atoms.

    Ranking is by abs(mass_error_da) ascending only — atom_count is never
    used as a ranking key (仕様書 §12.3, Claude Code rule #3: candidates with
    3+ atoms are not penalized just for having more atoms).

    ΔDa is rounded to a fixed grid (0.0001 Da) before the DFS search so that
    repeated calls with near-identical mass differences (common across many
    Mass Comparison rows) hit an lru_cache instead of re-running the search
    (仕様書 §12.2).
    """
    resolved_limits = element_limits if element_limits is not None else {e: max_total_atoms for e in elements}
    limits_key = tuple(sorted((e, resolved_limits.get(e, max_total_atoms)) for e in elements))
    rounded = round(mass_difference, 4)
    cached = _enumerate_cached(rounded, tolerance_da, max_total_atoms, limits_key, tuple(elements))
    results = list(cached)
    if max_candidates is not None:
        results = results[:max_candidates]
    return results


def rank_formula_candidates(candidates: list[FormulaCandidateResult]) -> list[FormulaCandidateResult]:
    """Explicit abs(mass_error_da)-ascending sort, callable independently of
    enumerate_formula_candidates (which already returns a sorted list)."""
    return sorted(candidates, key=lambda c: (abs(c.mass_error_da), c.formula))


def get_recommended_formula(candidates: list[FormulaCandidateResult]) -> FormulaCandidateResult | None:
    """Top-ranked candidate, or None if the list is empty."""
    return candidates[0] if candidates else None
