"""Fragment x charge x Peak matching and Mass Comparison row generation.

New module (仕様書 §6, §8.2, §15). Orchestrates existing, independently
reusable pieces:

- `mass_shift_ms1_search.build_sorted_peak_index` / `find_peaks_near_mz` for
  the fragment/charge -> peak binary search.
- `formula_candidate.enumerate_formula_candidates` for Formula Candidates.
- `modifications.find_modifications_by_mass_shift` for Modification
  Candidates.
- `ms2_support.find_ms2_support` (imported lazily, see below) for the
  optional MS2 reference columns.

Formula Candidate and Modification Candidate are computed independently and
written to separate fields (仕様書 §12.4, Claude Code rule #2) —
`formula_candidate.py` is never told about `modifications.py` or vice versa;
only this module combines their outputs into one row.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rna_masshunter.elemental_composition import _ELEMENT_ORDER
from rna_masshunter.formula_candidate import (
    FormulaCandidateResult,
    calculate_delta_mass,
    enumerate_formula_candidates,
    get_recommended_formula,
)
from rna_masshunter.mass_shift_ms1_search import build_sorted_peak_index, find_peaks_near_mz
from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import Fragment, Modification, Peak, RunConfig
from rna_masshunter.modifications import find_modifications_by_mass_shift
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass
from rna_masshunter.observed_mass import assign_peak_ids


@dataclass
class MassComparisonRow:
    peak_id: str
    fragment_id: str
    sequence: str
    charge: int
    intensity: float
    observed_mz: float
    observed_mass: float
    theoretical_mass: float
    delta_da: float
    delta_ppm: float
    recommended_formula: str | None
    formula_candidates: str
    known_modification: str | None
    modification_candidates: str
    ms2_spectrum_id: str | None = None
    ms2_matched_ion_count: int = 0
    ms2_matched_ions: str = ""
    rt: float | None = None
    scan_id: str | None = None
    warnings: list[str] = field(default_factory=list)


def _modification_display_name(modification: Modification) -> str:
    raw = modification.raw or {}
    return str(raw.get("name") or modification.symbol or modification.id)


def _modification_short_label(modification: Modification) -> str:
    return str(modification.symbol or modification.id)


def build_modification_candidates(
    delta_da: float,
    modifications: list[Modification],
    tolerance_da: float,
) -> list[Modification]:
    """仕様書 §13: known-modification search, sorted by |Δ| ascending
    (`find_modifications_by_mass_shift` itself does not sort)."""
    matches = find_modifications_by_mass_shift(modifications, delta_da, tolerance_da)
    return sorted(matches, key=lambda m: abs(m.mass_shift_from_unmodified - delta_da))


def _format_formula_candidates(candidates: list[FormulaCandidateResult]) -> str:
    parts = []
    for candidate in candidates:
        sign = "+" if candidate.mass_error_da >= 0 else "-"
        parts.append(f"{candidate.formula} [{sign}{abs(candidate.mass_error_da):.5f} Da]")
    return "; ".join(parts)


def _format_modification_candidates(delta_da: float, candidates: list[Modification]) -> str:
    parts = []
    for modification in candidates:
        diff = modification.mass_shift_from_unmodified - delta_da
        sign = "+" if diff >= 0 else "-"
        parts.append(f"{_modification_short_label(modification)} [Δ={sign}{abs(diff):.4f} Da]")
    return "; ".join(parts)


def build_mass_comparison_rows(
    fragments: list[Fragment],
    peaks: list[Peak],
    modifications: list[Modification],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None = None,
    ms2_spectra: list[Any] | None = None,
    ms2_ion_index: dict[str, list[Any]] | None = None,
) -> list[MassComparisonRow]:
    """仕様書 §8.2, §7 (data flow) の手順に対応:

    1. build_sorted_peak_index(peaks) でピーク集合を1回だけソート
    2. fragment x charge ごとに理論m/zを計算し、find_peaks_near_mz でマッチ抽出
    3. マッチごとに観測中性質量・ΔDa・Δppmを計算
    4. formula_candidate.enumerate_formula_candidates / build_modification_candidates
       を独立に呼び出す
    5. (ms2_spectra/ms2_ion_index が渡され、かつ config.ms2_annotation.enabled
       の場合のみ) ms2_support.find_ms2_support で参考情報を追加
    """
    mc_config = config.mass_comparison or {}
    if not mc_config.get("enabled", True):
        return []
    if not fragments or not peaks:
        return []

    fm_config = config.fragment_mapping or {}
    min_charge = int(fm_config.get("min_charge", 1) or 1)
    max_charge = int(fm_config.get("max_charge", 8) or 8)
    polarity = str(fm_config.get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str(config.instrument.get("polarity", "negative") or "negative").lower()
    tolerance_ppm = float(mc_config.get("mz_tolerance_ppm", 10) or 10)

    fc_config = config.formula_candidate or {}
    formula_enabled = bool(fc_config.get("enabled", True))
    tolerance_da = float(fc_config.get("mass_tolerance_da", 0.01) or 0.01)
    max_total_atoms = int(fc_config.get("max_total_atoms", 7) or 7)
    element_limits = fc_config.get("element_limits") or None
    elements = tuple(fc_config.get("elements") or _ELEMENT_ORDER)
    max_candidates = fc_config.get("max_candidates_per_match")
    include_zero_delta = bool(fc_config.get("include_zero_delta_as_match", True))

    peak_id_by_object = {id(peak): peak_id for peak, peak_id in zip(peaks, assign_peak_ids(peaks).values())}
    sorted_index = build_sorted_peak_index(peaks)

    # ms2_support is imported lazily so mass_comparison.py has no hard
    # import-time dependency on ms2_support.py (added in a later phase,
    # 仕様書 §22 phase 2.5) — with ms2_spectra/ms2_ion_index left at their
    # defaults, this branch is never taken and the import never happens.
    ms2_enabled = bool((config.ms2_annotation or {}).get("enabled")) and bool(ms2_spectra) and ms2_ion_index is not None
    find_ms2_support = None
    if ms2_enabled:
        from rna_masshunter.ms2_support import find_ms2_support as _find_ms2_support
        find_ms2_support = _find_ms2_support

    rows: list[MassComparisonRow] = []
    for fragment in fragments:
        for charge in range(min_charge, max_charge + 1):
            theoretical_mz = theoretical_mz_from_mass(fragment.unmodified_mass, charge, polarity)
            for match in find_peaks_near_mz(sorted_index, theoretical_mz, tolerance_ppm):
                peak = match.peak
                observed_mass = neutral_mass_from_mz(peak.mz, charge, polarity)
                delta_da = calculate_delta_mass(observed_mass, fragment.unmodified_mass)
                delta_ppm = ppm_error(observed_mass, fragment.unmodified_mass)

                if not include_zero_delta and round(delta_da, 5) == 0.0:
                    continue

                recommended_formula = None
                formula_candidates_str = ""
                if formula_enabled:
                    candidates = enumerate_formula_candidates(
                        delta_da,
                        tolerance_da,
                        max_total_atoms=max_total_atoms,
                        element_limits=element_limits,
                        elements=elements,
                        max_candidates=max_candidates,
                    )
                    top = get_recommended_formula(candidates)
                    recommended_formula = top.formula if top else None
                    formula_candidates_str = _format_formula_candidates(candidates)

                modification_matches = build_modification_candidates(delta_da, modifications, tolerance_da)
                known_modification = _modification_display_name(modification_matches[0]) if modification_matches else None
                modification_candidates_str = _format_modification_candidates(delta_da, modification_matches)

                row = MassComparisonRow(
                    peak_id=peak_id_by_object.get(id(peak), ""),
                    fragment_id=fragment.fragment_id,
                    sequence=fragment.sequence,
                    charge=charge,
                    intensity=peak.intensity,
                    observed_mz=peak.mz,
                    observed_mass=observed_mass,
                    theoretical_mass=fragment.unmodified_mass,
                    delta_da=delta_da,
                    delta_ppm=delta_ppm,
                    recommended_formula=recommended_formula,
                    formula_candidates=formula_candidates_str,
                    known_modification=known_modification,
                    modification_candidates=modification_candidates_str,
                    rt=peak.rt,
                    scan_id=peak.scan_id,
                )

                if find_ms2_support is not None:
                    support = find_ms2_support(
                        fragment.fragment_id, charge, peak.mz, ms2_spectra, ms2_ion_index, config,
                    )
                    if support is not None:
                        row.ms2_spectrum_id = support.spectrum_id
                        row.ms2_matched_ion_count = support.matched_ion_count
                        row.ms2_matched_ions = support.matched_ions

                rows.append(row)

    return rows
