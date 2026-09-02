import pytest

from rna_masshunter.mass_comparison import build_mass_comparison_rows
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Fragment, Modification, Peak, RunConfig


def _fragment(fragment_id: str, sequence: str, mass: float) -> Fragment:
    return Fragment(
        fragment_id=fragment_id,
        target_id="TEST",
        sequence=sequence,
        start=1,
        end=len(sequence),
        standard_start=1,
        standard_end=len(sequence),
        enzyme="RNase_T1",
        missed_cleavages=0,
        terminal_form="default",
        unmodified_mass=mass,
    )


def _config(**overrides) -> RunConfig:
    config = RunConfig(
        instrument={"polarity": "negative"},
        fragment_mapping={"enabled": True, "mz_tolerance_ppm": 10, "min_charge": 1, "max_charge": 3, "polarity": "auto"},
        formula_candidate={
            "enabled": True, "max_total_atoms": 7,
            "elements": ["C", "H", "N", "O", "P", "S", "Se"],
            "element_limits": {"C": 3, "H": 7, "N": 3, "O": 4, "P": 1, "S": 1, "Se": 1},
            "mass_tolerance_da": 0.01, "max_candidates_per_match": 10,
            "include_zero_delta_as_match": True,
        },
        mass_comparison={"enabled": True, "mz_tolerance_ppm": 10},
        ms2_annotation={"enabled": False},
    )
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


@pytest.fixture
def synthetic_fragments():
    return [
        _fragment("FRAG_A", "GACU", 1200.15),
        _fragment("FRAG_B", "GGCU", 1300.25),
    ]


@pytest.fixture
def synthetic_modifications():
    return [
        Modification(
            id="m1A", symbol="m1A", mass_shift_from_unmodified=14.0156,
            category="biological", target_bases=["A"], raw={"name": "1-methyladenosine"},
        ),
        Modification(
            id="ho5C", symbol="ho5C", mass_shift_from_unmodified=15.9949,
            category="biological", target_bases=["C"], raw={"name": "5-hydroxycytidine"},
        ),
    ]


def test_exact_match_row_has_zero_delta_and_no_known_modification(synthetic_fragments, synthetic_modifications):
    fragment = synthetic_fragments[0]
    charge = 2
    mz = mz_from_neutral_mass(fragment.unmodified_mass, charge, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0, rt=5.0, scan_id="scan_1")]

    rows = build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, _config())

    assert len(rows) == 1
    row = rows[0]
    assert row.fragment_id == "FRAG_A"
    assert row.charge == charge
    assert row.delta_da == pytest.approx(0.0, abs=1e-6)
    assert row.recommended_formula == "0"
    assert row.known_modification is None
    assert row.modification_candidates == ""
    assert row.peak_id == "PK0001"


def test_modified_match_reports_formula_and_modification_independently(synthetic_fragments, synthetic_modifications):
    fragment = synthetic_fragments[0]
    charge = 2
    modified_mass = fragment.unmodified_mass + 14.0156  # m1A-like shift
    mz = mz_from_neutral_mass(modified_mass, charge, "negative")
    peaks = [Peak(mz=mz, intensity=2000.0, rt=6.0, scan_id="scan_2")]

    rows = build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, _config())

    assert len(rows) == 1
    row = rows[0]
    assert row.known_modification == "1-methyladenosine"
    assert "m1A" in row.modification_candidates
    assert row.recommended_formula is not None
    assert row.formula_candidates != ""
    # Formula and Modification candidates are computed independently: one
    # being non-empty does not depend on the other (原則3).
    assert row.recommended_formula != row.known_modification


def test_formula_candidates_independent_of_zero_modifications(synthetic_fragments):
    fragment = synthetic_fragments[0]
    charge = 1
    mz = mz_from_neutral_mass(fragment.unmodified_mass + 15.9949, charge, "negative")
    peaks = [Peak(mz=mz, intensity=500.0)]

    rows = build_mass_comparison_rows(synthetic_fragments, peaks, [], _config())

    assert len(rows) == 1
    row = rows[0]
    assert row.known_modification is None
    assert row.modification_candidates == ""
    assert row.recommended_formula == "O1"


def test_ms2_columns_default_when_not_provided(synthetic_fragments, synthetic_modifications):
    fragment = synthetic_fragments[0]
    mz = mz_from_neutral_mass(fragment.unmodified_mass, 2, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]

    rows = build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, _config())

    assert len(rows) == 1
    row = rows[0]
    assert row.ms2_spectrum_id is None
    assert row.ms2_matched_ion_count == 0
    assert row.ms2_matched_ions == ""
    # Recommended Formula / Known Modification must not depend on MS2 at all.
    assert row.recommended_formula == "0"


def test_no_peaks_returns_empty_list(synthetic_fragments, synthetic_modifications):
    assert build_mass_comparison_rows(synthetic_fragments, [], synthetic_modifications, _config()) == []


def test_disabled_mass_comparison_returns_empty_list(synthetic_fragments, synthetic_modifications):
    fragment = synthetic_fragments[0]
    mz = mz_from_neutral_mass(fragment.unmodified_mass, 2, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]
    config = _config(mass_comparison={"enabled": False, "mz_tolerance_ppm": 10})

    assert build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, config) == []


def test_multiple_fragments_share_peak_id_numbering(synthetic_fragments, synthetic_modifications):
    fragment_a, fragment_b = synthetic_fragments
    mz_a = mz_from_neutral_mass(fragment_a.unmodified_mass, 2, "negative")
    mz_b = mz_from_neutral_mass(fragment_b.unmodified_mass, 2, "negative")
    peaks = [Peak(mz=mz_a, intensity=1000.0), Peak(mz=mz_b, intensity=1500.0)]

    rows = build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, _config())

    peak_ids = {row.fragment_id: row.peak_id for row in rows}
    assert peak_ids["FRAG_A"] == "PK0001"
    assert peak_ids["FRAG_B"] == "PK0002"


def test_include_zero_delta_as_match_false_drops_exact_matches(synthetic_fragments, synthetic_modifications):
    fragment = synthetic_fragments[0]
    mz = mz_from_neutral_mass(fragment.unmodified_mass, 2, "negative")
    peaks = [Peak(mz=mz, intensity=1000.0)]
    config = _config()
    config.formula_candidate = dict(config.formula_candidate, include_zero_delta_as_match=False)

    assert build_mass_comparison_rows(synthetic_fragments, peaks, synthetic_modifications, config) == []
