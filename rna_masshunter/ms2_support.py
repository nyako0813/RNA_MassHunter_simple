"""MS2 reference-information support for Mass Comparison rows.

New module (仕様書 §6, §8.5, §14A). Strictly "参考情報の付加" (reference
info only) — never influences Recommended Formula / Known Modification, and
never attempts precursor "rescue", evidence-level scoring, or modification
localization (仕様書 §14A のスコープ外の点、Claude Code rule #14). Only unmodified
d/w/a/z theoretical ions (from `ms2_extraction.generate_theoretical_ms2_ions`)
are ever matched against observed MS2 peaks.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rna_masshunter.models import Fragment, MS2SpectrumInfo, RunConfig, TheoreticalMS2Ion
from rna_masshunter.ms1_mapping import ppm_error
from rna_masshunter.ms2_extraction import generate_theoretical_ms2_ions


@dataclass(frozen=True)
class MS2SupportResult:
    spectrum_id: str
    matched_ion_count: int
    matched_ions: str


def build_ms2_ion_index(
    theoretical_fragments: list[Fragment],
    config: RunConfig,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> dict[str, list[TheoreticalMS2Ion]]:
    """`generate_theoretical_ms2_ions()` の結果を parent_fragment_id ごとに
    グルーピングする。simple_pipeline.py 側で1回だけ計算して使い回す想定
    (Mass Comparisonの行数分は再計算しない、仕様書 §8.5)。"""
    ions = generate_theoretical_ms2_ions(theoretical_fragments, config, base_masses, warnings)
    index: dict[str, list[TheoreticalMS2Ion]] = {}
    for ion in ions:
        index.setdefault(ion.parent_fragment_id, []).append(ion)
    return index


def _select_precursor_spectrum(
    charge: int,
    observed_precursor_mz: float,
    ms2_spectra: list[MS2SpectrumInfo],
    tolerance_ppm: float,
) -> MS2SpectrumInfo | None:
    candidates = []
    for spectrum in ms2_spectra:
        if spectrum.precursor_mz is None:
            continue
        if spectrum.precursor_charge is not None and spectrum.precursor_charge != charge:
            continue
        if abs(ppm_error(observed_precursor_mz, spectrum.precursor_mz)) > tolerance_ppm:
            continue
        candidates.append(spectrum)
    if not candidates:
        return None
    # 複数一致時は base_peak_intensity 最大のものを採用する（タイブレークの
    # 単純化、原則2を踏まえ複雑な優先順位付けはしない。仕様書 §8.5）。
    return max(candidates, key=lambda s: s.base_peak_intensity or 0.0)


def _ion_label(ion: TheoreticalMS2Ion) -> str:
    return f"{ion.ion_type}{len(ion.ion_sequence)}"


def find_ms2_support(
    fragment_id: str,
    charge: int,
    observed_precursor_mz: float,
    ms2_spectra: list[MS2SpectrumInfo] | None,
    ms2_ion_index: dict[str, list[TheoreticalMS2Ion]] | None,
    config: RunConfig,
) -> MS2SupportResult | None:
    """仕様書 §8.5 の手順:

    1. precursor_match_tolerance_ppm 以内・charge一致（またはNone）のMS2ス
       ペクトルを探す。複数あれば base_peak_intensity 最大を採用。
    2. 一致するスペクトルが無ければ None（MS1のみのMass Comparison行は
       そのまま残り、MS2列が空欄になるだけ）。
    3. 採用したスペクトルの観測ピークを、fragment_id の未修飾理論d/w/a/z
       イオン（ms2_ion_index）に対し mz_tolerance_ppm 以内でマッチさせる。
    4. マッチしたイオンをppm誤差の絶対値昇順で返す。
    """
    if not ms2_spectra or not ms2_ion_index:
        return None

    ms2_config = config.ms2_annotation or {}
    precursor_tolerance_ppm = float(ms2_config.get("precursor_match_tolerance_ppm", 20) or 20)
    fragment_tolerance_ppm = float(ms2_config.get("mz_tolerance_ppm", 20) or 20)

    spectrum = _select_precursor_spectrum(charge, observed_precursor_mz, ms2_spectra, precursor_tolerance_ppm)
    if spectrum is None:
        return None

    ions = ms2_ion_index.get(fragment_id, [])
    matched: list[tuple[TheoreticalMS2Ion, float]] = []
    for ion in ions:
        best_err_ppm: float | None = None
        for peak_mz, _peak_intensity in spectrum.peaks:
            err_ppm = ppm_error(peak_mz, ion.theoretical_mz)
            if abs(err_ppm) <= fragment_tolerance_ppm and (best_err_ppm is None or abs(err_ppm) < abs(best_err_ppm)):
                best_err_ppm = err_ppm
        if best_err_ppm is not None:
            matched.append((ion, best_err_ppm))

    matched.sort(key=lambda pair: abs(pair[1]))
    matched_ions_str = "; ".join(
        f"{_ion_label(ion)} [{'+' if err_ppm >= 0 else '-'}{abs(err_ppm):.1f} ppm]"
        for ion, err_ppm in matched
    )

    return MS2SupportResult(
        spectrum_id=spectrum.spectrum_id,
        matched_ion_count=len(matched),
        matched_ions=matched_ions_str,
    )
