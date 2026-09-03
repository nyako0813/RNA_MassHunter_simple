import pytest

from rna_masshunter.models import Peak
from rna_masshunter.peak_picking import merge_adjacent_profile_points, merge_peaks_across_scans


def test_mountain_shaped_cluster_merges_to_apex_point():
    # A single true ion sampled as a rising/falling profile trace within one
    # scan; the apex (highest intensity) should survive as the sole
    # representative.
    peaks = [
        Peak(mz=500.00000, intensity=100.0, rt=1.0, scan_id="scan_1"),
        Peak(mz=500.00002, intensity=400.0, rt=1.0, scan_id="scan_1"),
        Peak(mz=500.00004, intensity=900.0, rt=1.0, scan_id="scan_1"),  # apex
        Peak(mz=500.00006, intensity=350.0, rt=1.0, scan_id="scan_1"),
        Peak(mz=500.00008, intensity=90.0, rt=1.0, scan_id="scan_1"),
    ]

    merged = merge_adjacent_profile_points(peaks, tolerance_ppm=10)

    assert len(merged) == 1
    assert merged[0].intensity == 900.0
    assert merged[0].mz == pytest.approx(500.00004)


def test_distinct_ions_within_tolerance_window_but_far_in_chain_stay_separate():
    # Two well-separated true peaks in the same scan must not be merged into
    # one just because tolerance_ppm is generous.
    peaks = [
        Peak(mz=500.00000, intensity=800.0, rt=1.0, scan_id="scan_1"),
        Peak(mz=600.00000, intensity=700.0, rt=1.0, scan_id="scan_1"),
    ]

    merged = merge_adjacent_profile_points(peaks, tolerance_ppm=10)

    assert len(merged) == 2
    assert {p.mz for p in merged} == {500.0, 600.0}


def test_peaks_in_different_scans_are_never_merged():
    # Same m/z, same instant almost, but different scans (different RT) —
    # this is chromatographic peak width across scans, not profile-mode
    # splitting within one scan, so it must be left alone.
    peaks = [
        Peak(mz=500.00001, intensity=800.0, rt=1.00, scan_id="scan_1"),
        Peak(mz=500.00001, intensity=850.0, rt=1.05, scan_id="scan_2"),
    ]

    merged = merge_adjacent_profile_points(peaks, tolerance_ppm=10)

    assert len(merged) == 2
    assert {p.scan_id for p in merged} == {"scan_1", "scan_2"}


def test_chain_linkage_merges_a_run_even_if_endpoints_exceed_tolerance():
    # Each consecutive pair is within tolerance of its neighbor, but the
    # first and last points, taken directly, would not be — chain linkage
    # (compare to the last point added to the cluster) should still merge
    # the whole run into one representative.
    tol_ppm = 10
    step = 500.0 * tol_ppm / 1e6 * 0.5  # half the tolerance window per step
    mzs = [500.0 + i * step for i in range(6)]
    peaks = [Peak(mz=mz, intensity=float(10 * (i + 1)), rt=1.0, scan_id="scan_1") for i, mz in enumerate(mzs)]

    merged = merge_adjacent_profile_points(peaks, tolerance_ppm=tol_ppm)

    assert len(merged) == 1
    assert merged[0].intensity == 60.0  # the last (highest-intensity) point


def test_empty_input_returns_empty_list():
    assert merge_adjacent_profile_points([], tolerance_ppm=10) == []


def test_single_peak_passthrough():
    peaks = [Peak(mz=500.0, intensity=123.0, rt=1.0, scan_id="scan_1")]
    merged = merge_adjacent_profile_points(peaks, tolerance_ppm=10)
    assert merged == peaks


# --- merge_peaks_across_scans (§14C) ----------------------------------------

def _scan_peak(scan_index: int, mz: float, intensity: float) -> Peak:
    return Peak(mz=mz, intensity=intensity, rt=1.0 + scan_index * 0.02, scan_id=f"scan_{scan_index}")


def test_same_ion_across_consecutive_scans_merges_with_scan_count_and_rt_range():
    peaks = [_scan_peak(i, 500.0, float(10 * (i + 1))) for i in range(5)]  # ranks 0-4, apex at rank 4

    merged = merge_peaks_across_scans(peaks, tolerance_ppm=10, max_scan_gap=2)

    assert len(merged) == 1
    result = merged[0]
    assert result.intensity == 50.0
    assert result.scan_count == 5
    assert result.rt_range == pytest.approx((1.0, 1.08))
    assert result.scan_id == "scan_4"  # representative = highest-intensity scan


def _peaks_with_scan_rank_gap(gap: int) -> list[Peak]:
    """Build a target ion (m/z 500) detected at scan-rank 0 and scan-rank
    `gap`, with distractor peaks (far-away m/z, never chain-linkable to the
    target) filling every rank in between — otherwise those intermediate
    scan_ids simply would not exist in the list, and the *rank* gap would
    collapse to 1 regardless of how far apart the chosen scan indices look
    (rank is assigned only to scan_ids actually present, see
    merge_peaks_across_scans' docstring)."""
    peaks = [_scan_peak(0, 500.0, 800.0)]
    for i in range(1, gap):
        peaks.append(_scan_peak(i, 900.0, 50.0))
    peaks.append(_scan_peak(gap, 500.0, 850.0))
    return peaks


def test_gap_beyond_max_scan_gap_keeps_detections_separate():
    peaks = _peaks_with_scan_rank_gap(5)

    merged = merge_peaks_across_scans(peaks, tolerance_ppm=10, max_scan_gap=2)

    target_ion_results = [p for p in merged if p.mz == 500.0]
    assert len(target_ion_results) == 2
    assert all(p.scan_count == 1 for p in target_ion_results)


def test_gap_exactly_at_max_scan_gap_boundary_merges():
    peaks = _peaks_with_scan_rank_gap(2)  # rank gap == 2 == max_scan_gap

    merged = merge_peaks_across_scans(peaks, tolerance_ppm=10, max_scan_gap=2)

    target_ion_results = [p for p in merged if p.mz == 500.0]
    assert len(target_ion_results) == 1
    assert target_ion_results[0].scan_count == 2


def test_gap_one_beyond_max_scan_gap_boundary_stays_separate():
    peaks = _peaks_with_scan_rank_gap(3)  # rank gap == 3 > max_scan_gap == 2

    merged = merge_peaks_across_scans(peaks, tolerance_ppm=10, max_scan_gap=2)

    target_ion_results = [p for p in merged if p.mz == 500.0]
    assert len(target_ion_results) == 2


def test_distinct_mz_in_adjacent_scans_never_merges_regardless_of_scan_gap():
    peaks = [_scan_peak(0, 500.0, 800.0), _scan_peak(1, 600.0, 700.0)]

    merged = merge_peaks_across_scans(peaks, tolerance_ppm=10, max_scan_gap=2)

    assert len(merged) == 2
    assert {p.scan_count for p in merged} == {1}


def test_merge_peaks_across_scans_empty_and_single_passthrough():
    assert merge_peaks_across_scans([], tolerance_ppm=10, max_scan_gap=2) == []
    single = [Peak(mz=500.0, intensity=1.0, rt=1.0, scan_id="scan_0")]
    merged = merge_peaks_across_scans(single, tolerance_ppm=10, max_scan_gap=2)
    assert len(merged) == 1
    assert merged[0].scan_count == 1
    assert merged[0].rt_range == (1.0, 1.0)
