"""P1 complete-digestion mode: match observed MS1 peaks directly against the
known-nucleoside mass universe (仕様書 §24, Phase 10).

Unlike mass_comparison.py's wide `max_delta_da` search (built to explore an
*unknown* ΔDa against oligomer fragments), this uses a narrow ppm tolerance
match against a small, fixed set of *known* target masses — there is no
"explain this residual mass difference" step here (§24 scope: known-
nucleoside matching only). `mass_shift_ms1_search.find_peaks_near_mz` (the
same binary-search engine `mass_comparison.py` deliberately avoided for its
narrow-ppm-around-the-*unmodified*-mass limitation) is exactly the right
tool here, since every target IS a specific known mass.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rna_masshunter.mass_shift_ms1_search import build_sorted_peak_index, find_peaks_near_mz
from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import Peak, RunConfig
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass
from rna_masshunter.nucleoside_targets import NucleosideTarget, theoretical_mass_for_target
from rna_masshunter.observed_mass import assign_peak_ids


@dataclass
class NucleosideComparisonRow:
    peak_id: str
    charge: int
    intensity: float
    observed_mz: float
    observed_mass: float
    target_label: str
    target_name: str
    target_category: str  # "standard" | "modified"
    theoretical_mass: float
    delta_da: float
    delta_ppm: float
    rt: float | None = None
    scan_id: str | None = None
    scan_count: int = 1
    rt_range: tuple[float, float] | None = None
    warnings: list[str] = field(default_factory=list)


def build_nucleoside_comparison_rows(
    peaks: list[Peak],
    targets: list[NucleosideTarget],
    config: RunConfig,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[NucleosideComparisonRow]:
    """仕様書 §24 の手順:

    1. config.alkaline_phosphatase.enabled から脱リン酸化の有無を決定
       （P1モード専用のconfigキーは増やさず既存キーを流用）
    2. 各target（標準4塩基 + modifications.yaml収載の修飾ヌクレオシド）
       について、charge範囲(fragment_mappingを流用)ごとに理論m/zを計算
    3. mass_shift_ms1_search.find_peaks_near_mz で狭いppm許容差
       (mass_comparison.nucleoside_mz_tolerance_ppm、既定10ppm)内の
       ピークを検索する（mass_comparison.pyの広いmax_delta_da窓とは
       異なり、既知の個別質量に対する直接照合のため）
    4. マッチごとに観測中性質量・ΔDa・Δppmを計算して行を生成する
    """
    mc_config = config.mass_comparison or {}
    if not mc_config.get("enabled", True):
        return []
    if not peaks or not targets:
        return []

    fm_config = config.fragment_mapping or {}
    min_charge = int(fm_config.get("min_charge", 1) or 1)
    max_charge = int(fm_config.get("max_charge", 8) or 8)
    polarity = str(fm_config.get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str(config.instrument.get("polarity", "negative") or "negative").lower()

    tolerance_ppm = float(mc_config.get("nucleoside_mz_tolerance_ppm", 10) or 10)
    dephosphorylated = bool((config.alkaline_phosphatase or {}).get("enabled", False))

    peak_id_by_object = {id(peak): peak_id for peak, peak_id in zip(peaks, assign_peak_ids(peaks).values())}
    sorted_index = build_sorted_peak_index(peaks)

    rows: list[NucleosideComparisonRow] = []
    for target in targets:
        theoretical_mass = theoretical_mass_for_target(target, dephosphorylated, base_masses)
        for charge in range(min_charge, max_charge + 1):
            theoretical_mz = theoretical_mz_from_mass(theoretical_mass, charge, polarity)
            for match in find_peaks_near_mz(sorted_index, theoretical_mz, tolerance_ppm):
                peak = match.peak
                observed_mass = neutral_mass_from_mz(peak.mz, charge, polarity)
                delta_da = observed_mass - theoretical_mass
                delta_ppm = ppm_error(observed_mass, theoretical_mass)
                rows.append(NucleosideComparisonRow(
                    peak_id=peak_id_by_object.get(id(peak), ""),
                    charge=charge,
                    intensity=peak.intensity,
                    observed_mz=peak.mz,
                    observed_mass=observed_mass,
                    target_label=target.label,
                    target_name=target.name,
                    target_category=target.category,
                    theoretical_mass=theoretical_mass,
                    delta_da=delta_da,
                    delta_ppm=delta_ppm,
                    rt=peak.rt,
                    scan_id=peak.scan_id,
                    scan_count=getattr(peak, "scan_count", 1),
                    rt_range=getattr(peak, "rt_range", None),
                ))

    return rows
