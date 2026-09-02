from pathlib import Path

import pytest

from rna_masshunter.masses import load_base_masses
from rna_masshunter.models import Fragment, MS2SpectrumInfo, RunConfig
from rna_masshunter.ms2_support import build_ms2_ion_index, find_ms2_support

REPO_ROOT = Path(__file__).resolve().parent.parent


def _config() -> RunConfig:
    return RunConfig(
        instrument={"polarity": "negative"},
        ms2_annotation={
            "enabled": True,
            "mz_tolerance_ppm": 20,
            "precursor_match_tolerance_ppm": 20,
            "use_theoretical_fragments": True,
            "ion_series": ["d", "w", "a", "z"],
            "charge_states": [1],
            "min_ion_length": 1,
            "max_ion_length": None,
        },
    )


def _fragment(fragment_id: str, sequence: str) -> Fragment:
    return Fragment(
        fragment_id=fragment_id, target_id="TEST", sequence=sequence,
        start=1, end=len(sequence), standard_start=1, standard_end=len(sequence),
        enzyme="RNase_T1", missed_cleavages=0, terminal_form="default",
        unmodified_mass=0.0,
    )


def _spectrum(spectrum_id, precursor_mz, precursor_charge, peaks, base_peak_intensity=1000.0):
    return MS2SpectrumInfo(
        spectrum_id=spectrum_id, scan_index=1, rt=5.0,
        precursor_mz=precursor_mz, precursor_charge=precursor_charge, precursor_intensity=1000.0,
        num_peaks=len(peaks), base_peak_mz=peaks[0][0] if peaks else None,
        base_peak_intensity=base_peak_intensity, total_ion_current=sum(p[1] for p in peaks),
        peaks=peaks,
    )


@pytest.fixture
def base_masses():
    return load_base_masses(REPO_ROOT / "data" / "base_masses.yaml")


@pytest.fixture
def synthetic_fragments():
    return [_fragment("FRAG_A", "GACU"), _fragment("FRAG_B", "GGCU")]


@pytest.fixture
def ion_index(synthetic_fragments, base_masses):
    return build_ms2_ion_index(synthetic_fragments, _config(), base_masses)


def test_build_ms2_ion_index_groups_by_fragment(synthetic_fragments, base_masses):
    index = build_ms2_ion_index(synthetic_fragments, _config(), base_masses)
    assert set(index) == {"FRAG_A", "FRAG_B"}
    assert all(ion.parent_fragment_id == "FRAG_A" for ion in index["FRAG_A"])
    # d/w/a/z series over a 4-mer, 3 cut points -> at least a handful of ions
    assert len(index["FRAG_A"]) > 0


def test_no_matching_precursor_returns_none(ion_index):
    spectra = [_spectrum("scan_1", precursor_mz=9999.0, precursor_charge=2, peaks=[(500.0, 100.0)])]
    result = find_ms2_support("FRAG_A", 2, 100.0, spectra, ion_index, _config())
    assert result is None


def test_unknown_fragment_id_returns_no_matched_ions(ion_index):
    spectra = [_spectrum("scan_1", precursor_mz=100.0, precursor_charge=2, peaks=[(500.0, 100.0)])]
    result = find_ms2_support("FRAG_NOT_PRESENT", 2, 100.0, spectra, ion_index, _config())
    assert result is not None
    assert result.matched_ion_count == 0
    assert result.matched_ions == ""


def test_matched_ions_sorted_by_ppm_error(ion_index):
    target_ion = ion_index["FRAG_A"][0]
    other_ion = ion_index["FRAG_A"][1]
    # exact hit for target_ion, a small offset for other_ion so it ranks second
    offset_mz = other_ion.theoretical_mz * (1 + 5e-6)
    spectra = [_spectrum(
        "scan_1", precursor_mz=250.0, precursor_charge=1,
        peaks=[(offset_mz, 200.0), (target_ion.theoretical_mz, 500.0)],
    )]

    result = find_ms2_support("FRAG_A", 1, 250.0, spectra, ion_index, _config())

    assert result is not None
    assert result.spectrum_id == "scan_1"
    assert result.matched_ion_count >= 2
    first_label = result.matched_ions.split(";")[0]
    assert first_label.startswith(f"{target_ion.ion_type}{len(target_ion.ion_sequence)}")


def test_multiple_candidate_spectra_prefers_highest_base_peak_intensity(ion_index):
    target_ion = ion_index["FRAG_A"][0]
    weak = _spectrum("weak", precursor_mz=250.0, precursor_charge=1, peaks=[(target_ion.theoretical_mz, 50.0)], base_peak_intensity=50.0)
    strong = _spectrum("strong", precursor_mz=250.0, precursor_charge=1, peaks=[(target_ion.theoretical_mz, 900.0)], base_peak_intensity=900.0)

    result = find_ms2_support("FRAG_A", 1, 250.0, [weak, strong], ion_index, _config())

    assert result is not None
    assert result.spectrum_id == "strong"


def test_precursor_charge_none_is_treated_as_wildcard(ion_index):
    target_ion = ion_index["FRAG_A"][0]
    spectra = [_spectrum("scan_1", precursor_mz=250.0, precursor_charge=None, peaks=[(target_ion.theoretical_mz, 500.0)])]

    result = find_ms2_support("FRAG_A", 3, 250.0, spectra, ion_index, _config())

    assert result is not None
    assert result.matched_ion_count >= 1


def test_no_spectra_or_no_index_returns_none(ion_index):
    assert find_ms2_support("FRAG_A", 1, 250.0, [], ion_index, _config()) is None
    assert find_ms2_support("FRAG_A", 1, 250.0, None, ion_index, _config()) is None
    assert find_ms2_support("FRAG_A", 1, 250.0, [_spectrum("s", 250.0, 1, [(1.0, 1.0)])], {}, _config()) is None
