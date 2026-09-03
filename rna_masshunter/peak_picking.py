from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from rna_masshunter.models import Peak
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra


def extract_ms1_peaks(mzml_path: str | Path, reconstruction_config: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> list[Peak]:
    rt_min = reconstruction_config.get("rt_min")
    rt_max = reconstruction_config.get("rt_max")
    mz_min = float(reconstruction_config.get("mz_min", 0))
    mz_max = float(reconstruction_config.get("mz_max", float("inf")))
    intensity_threshold = float(reconstruction_config.get("intensity_threshold", 0))
    peaks: list[Peak] = []

    for spectrum in iter_spectra(mzml_path):
        if int(spectrum.get("ms level", 0)) != 1:
            continue
        rt = _rt_minutes(spectrum)
        if rt is not None and rt_min is not None and rt < float(rt_min):
            continue
        if rt is not None and rt_max is not None and rt > float(rt_max):
            continue
        mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
        intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
        if mz_array.size != intensity_array.size:
            continue
        mask = (mz_array >= mz_min) & (mz_array <= mz_max) & (intensity_array >= intensity_threshold)
        for mz_value, intensity in zip(mz_array[mask], intensity_array[mask], strict=False):
            peaks.append(Peak(mz=float(mz_value), intensity=float(intensity), rt=rt, scan_id=str(spectrum.get("id", ""))))
    return peaks


def merge_adjacent_profile_points(peaks: list[Peak], tolerance_ppm: float) -> list[Peak]:
    """Collapse raw profile-mode sample points that belong to the same true
    ion back down to one representative point per ion.

    New function (project addition, not vendored — this project's real
    mzML input turned out to be profile-mode rather than centroided; see
    docs/design/RNA_MassHunter_再設計_実装仕様書.md §14B). `extract_ms1_peaks`
    above does no centroiding, so a single ion's continuous m/z envelope in
    profile data comes back as several `Peak`s a few micro-Da apart instead
    of one — inflating 04_Observed_Mass/06_Mass_Comparison with near-
    duplicate rows for what is physically one observation.

    Deliberately minimal: this is chain-linked ppm clustering + "take the
    tallest point in the cluster", not real centroiding (no baseline
    subtraction, no isotope handling) — that is out of scope for this
    project (仕様書 §14: "簡易版ではclassify_peak_tiers...必須としない"
    applies in the same spirit here).

    Clustering is scoped to peaks sharing the same `scan_id` (i.e. the same
    RT): profile-point splitting happens *within* one scan's m/z trace.
    Peaks from different scans are never merged with each other here, even
    if close in m/z — a chromatographic peak legitimately spanning several
    scans is a different, meaningful signal (peak shape over RT), not the
    artifact this function targets.

    `peaks` are grouped by `scan_id`, sorted by m/z within each group, and
    walked once: consecutive peaks within `tolerance_ppm` of the
    previously-seen point in the current cluster are chained together: for
    each finished cluster, the single highest-intensity `Peak` is kept as
    its representative and the rest are dropped.
    """
    if not peaks:
        return []

    by_scan: dict[Any, list[Peak]] = {}
    for peak in peaks:
        by_scan.setdefault(peak.scan_id, []).append(peak)

    merged: list[Peak] = []
    for scan_peaks in by_scan.values():
        ordered = sorted(scan_peaks, key=lambda p: p.mz)
        cluster = [ordered[0]]
        for peak in ordered[1:]:
            window = cluster[-1].mz * tolerance_ppm / 1e6
            if peak.mz - cluster[-1].mz <= window:
                cluster.append(peak)
            else:
                merged.append(max(cluster, key=lambda p: p.intensity))
                cluster = [peak]
        merged.append(max(cluster, key=lambda p: p.intensity))
    return merged


def merge_peaks_across_scans(peaks: list[Peak], tolerance_ppm: float, max_scan_gap: int) -> list[Peak]:
    """Collapse the same ion's repeated per-scan detections, across a run of
    nearby scans, into one representative row.

    New function (project addition; 仕様書 §14C). This is the sequel to
    `merge_adjacent_profile_points` above (§14B): after that function has
    already collapsed *within-scan* profile splitting, a real chromatographic
    peak still legitimately produces one detection of the same ion in each
    of several *consecutive* MS1 scans as it elutes (confirmed against this
    project's real data: e.g. the same "AG" fragment detected across
    scan cycles 388 and 402-405 within ~0.08 min). Kept unmerged, that is
    what was still inflating 06_Mass_Comparison with several near-identical
    rows per fragment even after §14B. This function merges *that* — but
    only within a bounded number of skipped scans (`max_scan_gap`), so two
    genuinely separate detection events (e.g. the same compound eluting
    twice, far apart) are not silently collapsed into one.

    Deliberately the last extension in this direction (per user instruction,
    2026-09-03): no smoothing, no baseline estimation, no real EIC/XIC apex
    detection — just chain-linked ppm clustering across a scan-index window,
    the same minimal philosophy as §14B.

    Unlike `merge_adjacent_profile_points`, information IS lost by merging
    here (which scan, which exact RT) unless preserved explicitly — so each
    merged `Peak.scan_count` records how many scans contributed, and
    `Peak.rt_range` records the (min, max) RT span, while the ordinary
    `rt`/`scan_id` fields keep the representative (highest-intensity)
    scan's own values. Callers that surface `Peak`s in a report should show
    `scan_count`/`rt_range` alongside `rt`/`scan_id` so this merge stays
    visible rather than silently discarding reproducibility information.

    Scan order/gap is derived from the RT of each distinct `scan_id` present
    in `peaks` itself (peaks are extracted from MS1 scans in acquisition
    order, so ranking distinct scan_ids by RT reconstructs that order) —
    there is no separate "full scan list" input. A peak joins a cluster only
    if it is within `tolerance_ppm` of the cluster's current chain endpoint
    (by m/z, same rule as §14B) AND its scan rank is within `max_scan_gap`
    of at least one scan rank already in the cluster.

    Caveat (inherent to any ppm+scan-gap clustering, not fixable without
    real XIC peak detection): two distinct co-eluting compounds that happen
    to fall within `tolerance_ppm` of each other could be merged if their
    detections overlap in scan range. `max_scan_gap` is a provisional
    default (2) — see 仕様書 §14C for the recommendation to check its
    distribution against real data and spot-check boundary cases before
    trusting it unreviewed on a new dataset.
    """
    if not peaks:
        return []

    scan_rt: dict[Any, float] = {}
    for peak in peaks:
        if peak.rt is not None and peak.scan_id not in scan_rt:
            scan_rt[peak.scan_id] = peak.rt
    scan_rank = {scan_id: rank for rank, scan_id in enumerate(sorted(scan_rt, key=lambda sid: scan_rt[sid]))}

    def _rank(peak: Peak) -> int | None:
        return scan_rank.get(peak.scan_id)

    ordered = sorted(peaks, key=lambda p: p.mz)
    clusters: list[list[Peak]] = [[ordered[0]]]
    cluster_ranks: list[set[int]] = [{r for r in [_rank(ordered[0])] if r is not None}]

    for peak in ordered[1:]:
        cluster, ranks = clusters[-1], cluster_ranks[-1]
        window = cluster[-1].mz * tolerance_ppm / 1e6
        mz_ok = peak.mz - cluster[-1].mz <= window
        rank = _rank(peak)
        gap_ok = rank is None or not ranks or any(abs(rank - r) <= max_scan_gap for r in ranks)
        if mz_ok and gap_ok:
            cluster.append(peak)
            if rank is not None:
                ranks.add(rank)
        else:
            clusters.append([peak])
            cluster_ranks.append({r for r in [rank] if r is not None})

    merged: list[Peak] = []
    for cluster in clusters:
        representative = max(cluster, key=lambda p: p.intensity)
        rts = [p.rt for p in cluster if p.rt is not None]
        rt_range = (min(rts), max(rts)) if rts else None
        merged.append(replace(representative, scan_count=len(cluster), rt_range=rt_range))
    return merged
