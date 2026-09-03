"""Elemental-composition ("Formula Candidate") enumeration for a ΔDa value.

New module (no vendored counterpart) implementing 仕様書
docs/design/RNA_MassHunter_再設計_実装仕様書.md §8.1, §12.

Deliberately does not import `rna_masshunter.modifications` — Formula
Candidate and Modification Candidate must stay fully independent code paths
(仕様書 §12.4, Claude Code rule #2).
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
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


@dataclass(frozen=True)
class _Composition:
    formula: str
    exact_mass: float
    atom_count: int
    elements: tuple[tuple[str, int], ...]


@lru_cache(maxsize=64)
def _composition_universe(
    max_total_atoms: int,
    element_limits: tuple[tuple[str, int], ...],
    elements: tuple[str, ...],
) -> tuple[tuple[float, ...], tuple[_Composition, ...]]:
    """Enumerate *every* elemental composition within sum(|count|) <=
    max_total_atoms for a given (max_total_atoms, element_limits, elements)
    config, independent of any target ΔDa — the fixed "universe" of
    achievable mass shifts. `enumerate_formula_candidates` below binary-
    searches this (sorted-by-mass, cached) universe instead of re-running
    the DFS for every ΔDa query.

    2026-09-03 rewrite: `enumerate_formula_candidates` previously ran a
    fresh DFS per call (with an lru_cache keyed on ΔDa rounded to 0.0001
    Da — 仕様書 §12.2's original approach). That works fine when only the
    ~50 post-truncation rows per fragment x charge need a search, but once
    Mass Comparison started running Formula/Modification search on *every*
    matched peak *before* truncating (so truncation can prioritize rows
    with a candidate — see mass_comparison.py), the number of DFS calls
    could reach the tens of thousands per pipeline run and each call rarely
    hit the old cache (real ΔDa values are continuous, not repeating to
    0.0001 Da) — a real-data run that used to take ~1m17s never finished in
    over an hour. The DFS itself doesn't depend on the target ΔDa at all,
    only `tolerance_da` gates which of its results count as a "match" — so
    computing the whole reachable-mass universe once per (max_total_atoms,
    element_limits, elements) config and reusing it via bisect for every
    query turns O(queries * DFS) into O(one DFS + queries * log(universe)).
    """
    limits = dict(element_limits)
    results: list[_Composition] = []
    counts: dict[str, int] = {}

    def dfs(idx: int, remaining_atoms: int, running_mass: float) -> None:
        if idx == len(elements):
            comp = ElementalComposition.delta(dict(counts))
            results.append(_Composition(
                formula=comp.canonical_string(),
                exact_mass=comp.exact_mass,
                atom_count=sum(abs(v) for v in counts.values()),
                elements=tuple(sorted(comp.to_dict().items())),
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
    results.sort(key=lambda c: c.exact_mass)
    masses = tuple(c.exact_mass for c in results)
    return masses, tuple(results)


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

    Internally this binary-searches a cached, precomputed universe of every
    achievable composition for the given (max_total_atoms, element_limits,
    elements) — see `_composition_universe`'s docstring for why (仕様書
    §12.2, revised 2026-09-03). `mass_difference` itself is used exactly as
    given (no rounding) — precision is not traded away by this change,
    only the old per-call rounding-for-caching is gone since it is no
    longer needed.
    """
    resolved_limits = element_limits if element_limits is not None else {e: max_total_atoms for e in elements}
    limits_key = tuple(sorted((e, resolved_limits.get(e, max_total_atoms)) for e in elements))
    masses, universe = _composition_universe(max_total_atoms, limits_key, tuple(elements))

    lo = bisect_left(masses, mass_difference - tolerance_da)
    hi = bisect_right(masses, mass_difference + tolerance_da)

    results = [
        FormulaCandidateResult(
            formula=comp.formula,
            exact_mass=comp.exact_mass,
            mass_error_da=comp.exact_mass - mass_difference,
            mass_error_ppm=None,
            atom_count=comp.atom_count,
            elements=dict(comp.elements),
        )
        for comp in universe[lo:hi]
    ]
    results.sort(key=lambda r: (abs(r.mass_error_da), r.formula))
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
