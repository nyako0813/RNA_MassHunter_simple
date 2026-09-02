"""Shared, efficient MS1 mass-shift search.

Peaks are sorted once by m/z; each fragment x shift x charge combination is
looked up via binary search instead of a full linear scan, so this stays
fast even with many candidate mass shifts (known modifications or simple
elemental deltas) against tens of thousands of MS1 peaks.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass
from typing import Any

from rna_masshunter.ms1_mapping import ppm_error, _confidence


@dataclass
class SortedPeakIndex:
    peaks: list[Any]
    mzs: list[float]


def build_sorted_peak_index(peaks: list[Any]) -> SortedPeakIndex:
    ordered = sorted(peaks, key=lambda p: float(getattr(p, "mz")))
    return SortedPeakIndex(peaks=ordered, mzs=[float(getattr(p, "mz")) for p in ordered])


@dataclass
class PeakMatch:
    peak: Any
    observed_mz: float
    error_ppm: float
    error_da: float
    confidence: str


def find_peaks_near_mz(index: SortedPeakIndex, theoretical_mz: float, tolerance_ppm: float) -> list[PeakMatch]:
    """Binary-search `index` for peaks within `tolerance_ppm` of theoretical_mz."""
    if not index.mzs or theoretical_mz <= 0:
        return []
    window = theoretical_mz * tolerance_ppm / 1e6
    lo = bisect.bisect_left(index.mzs, theoretical_mz - window)
    hi = bisect.bisect_right(index.mzs, theoretical_mz + window)
    matches: list[PeakMatch] = []
    for peak in index.peaks[lo:hi]:
        observed_mz = float(getattr(peak, "mz"))
        error_ppm = ppm_error(observed_mz, theoretical_mz)
        if abs(error_ppm) > tolerance_ppm:
            continue
        error_da = observed_mz - theoretical_mz
        tier = getattr(peak, "tier", None)
        confidence = _confidence(error_ppm, tolerance_ppm, tier)
        matches.append(PeakMatch(peak=peak, observed_mz=observed_mz, error_ppm=error_ppm, error_da=error_da, confidence=confidence))
    return matches
