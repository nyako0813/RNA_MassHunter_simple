import pytest

from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Peak, RunConfig
from rna_masshunter.nucleoside_comparison import build_nucleoside_comparison_rows
from rna_masshunter.nucleoside_targets import NucleosideTarget, build_standard_nucleoside_targets


def _config(**overrides) -> RunConfig:
    config = RunConfig(
        instrument={"polarity": "negative"},
        fragment_mapping={"enabled": True, "min_charge": 1, "max_charge": 2, "polarity": "auto"},
        mass_comparison={"enabled": True, "nucleoside_mz_tolerance_ppm": 10},
        alkaline_phosphatase={"enabled": True},  # dephosphorylated by default in these tests
    )
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


@pytest.fixture
def standard_targets() -> list[NucleosideTarget]:
    return build_standard_nucleoside_targets()


def test_exact_match_to_standard_nucleoside(standard_targets):
    adenosine = next(t for t in standard_targets if t.label == "A")
    charge = 1
    mz = mz_from_neutral_mass(adenosine.nucleoside_mass, charge, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0, rt=2.0, scan_id="scan_1")]

    rows = build_nucleoside_comparison_rows(peaks, standard_targets, _config(), base_masses={})

    matches = [r for r in rows if r.target_label == "A" and r.charge == charge]
    assert len(matches) == 1
    assert matches[0].delta_da == pytest.approx(0.0, abs=1e-6)
    assert matches[0].target_category == "standard"
    assert matches[0].peak_id == "PK0001"


def test_no_match_beyond_narrow_tolerance(standard_targets):
    adenosine = next(t for t in standard_targets if t.label == "A")
    charge = 1
    # 14 Da off — the whole point of the narrow-tolerance design is that
    # this must NOT match A; it's simply not in the known universe unless a
    # matching modified-nucleoside target also has this exact mass.
    mz = mz_from_neutral_mass(adenosine.nucleoside_mass + 14.0, charge, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]

    rows = build_nucleoside_comparison_rows(peaks, standard_targets, _config(), base_masses={})

    assert rows == []


def test_phosphorylated_mode_shifts_theoretical_mass(standard_targets):
    adenosine = next(t for t in standard_targets if t.label == "A")
    charge = 1
    base_masses = {"constants": {"phosphate": 79.966331}}
    # Peak matches the phosphorylated (5'-monophosphate) mass, not the free
    # nucleoside mass.
    mz = mz_from_neutral_mass(adenosine.nucleoside_mass + 79.966331, charge, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]

    config_dephospho = _config(alkaline_phosphatase={"enabled": True})
    config_phospho = _config(alkaline_phosphatase={"enabled": False})

    assert build_nucleoside_comparison_rows(peaks, standard_targets, config_dephospho, base_masses) == []
    rows = build_nucleoside_comparison_rows(peaks, standard_targets, config_phospho, base_masses)
    matches = [r for r in rows if r.target_label == "A"]
    assert len(matches) == 1
    assert matches[0].delta_da == pytest.approx(0.0, abs=1e-6)
    assert matches[0].theoretical_mass == pytest.approx(adenosine.nucleoside_mass + 79.966331, abs=1e-6)


def test_modified_target_matches_independently_of_standard(standard_targets):
    modified = NucleosideTarget(label="m1A", name="1-methyladenosine", category="modified", nucleoside_mass=281.1124)
    targets = standard_targets + [modified]
    charge = 1
    mz = mz_from_neutral_mass(281.1124, charge, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]

    rows = build_nucleoside_comparison_rows(peaks, targets, _config(), base_masses={})

    assert len(rows) == 1
    assert rows[0].target_label == "m1A"
    assert rows[0].target_category == "modified"


def test_no_peaks_or_no_targets_returns_empty_list(standard_targets):
    assert build_nucleoside_comparison_rows([], standard_targets, _config(), base_masses={}) == []
    peaks = [Peak(mz=300.0, intensity=1000.0)]
    assert build_nucleoside_comparison_rows(peaks, [], _config(), base_masses={}) == []


def test_disabled_mass_comparison_returns_empty_list(standard_targets):
    peaks = [Peak(mz=300.0, intensity=1000.0)]
    config = _config(mass_comparison={"enabled": False, "nucleoside_mz_tolerance_ppm": 10})
    assert build_nucleoside_comparison_rows(peaks, standard_targets, config, base_masses={}) == []


def test_scan_count_and_rt_range_carried_through_from_merged_peak(standard_targets):
    adenosine = next(t for t in standard_targets if t.label == "A")
    charge = 1
    mz = mz_from_neutral_mass(adenosine.nucleoside_mass, charge, "negative")
    peak = Peak(mz=mz, intensity=1000.0, rt=2.0, scan_id="scan_1", scan_count=5, rt_range=(1.9, 2.1))

    rows = build_nucleoside_comparison_rows([peak], standard_targets, _config(), base_masses={})

    matches = [r for r in rows if r.target_label == "A"]
    assert matches[0].scan_count == 5
    assert matches[0].rt_range == (1.9, 2.1)
