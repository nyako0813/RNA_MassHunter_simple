"""Smoke tests for the modules vendored from nyako0813/RNA_MassHunter.

These exist to catch import-chain regressions (the whole reason
ms2_extraction.py / cca_processing.py / models.py / config.py were
hand-trimmed instead of imported as-is — see their module docstrings and
docs/design/RNA_MassHunter_再設計_実装仕様書.md §3.4, §3.9). If any of
these fail, something about the vendoring broke; they are not meant to
test the (not-yet-implemented) new logic — formula_candidate.py,
mass_comparison.py, ms2_support.py, simple_pipeline.py.
"""
from pathlib import Path

from rna_masshunter.cca_processing import (
    CCAMaturationState,
    RegisteredSequenceCCAMode,
    process_cca_tail,
)
from rna_masshunter.config import DEFAULT_CONFIG, load_config, resolve_paths, validate_config
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses
from rna_masshunter.models import RunConfig
from rna_masshunter.modifications import find_modifications_by_mass_shift, load_modifications
from rna_masshunter.ms2_extraction import extract_ms2_spectra, generate_theoretical_ms2_ions

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_cca_processing_appends_missing_tail():
    result = process_cca_tail("ACGU", RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert result.processed_sequence == "ACGUCCA"
    assert result.original_tail_state == CCAMaturationState.NONE


def test_cca_processing_no_op_when_already_cca():
    result = process_cca_tail("ACGCCA", RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert result.processed_sequence == "ACGCCA"
    assert result.added_nucleotide_count == 0


def test_config_round_trip(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("sequence:\n  sequence: ACGU\n", encoding="utf-8")
    config = load_config(config_path)
    assert isinstance(config, RunConfig)
    assert config.sequence["sequence"] == "ACGU"
    # every DEFAULT_CONFIG section should be present even though the YAML
    # above only set one key in one section
    for section in DEFAULT_CONFIG:
        assert getattr(config, section) is not None
    validate_config(config)  # should not raise on defaults
    resolve_paths(config, REPO_ROOT)  # should not raise


def test_base_masses_and_unmodified_mass():
    base_masses = load_base_masses(REPO_ROOT / "data" / "base_masses.yaml")
    mass = calculate_unmodified_rna_mass("ACGU", base_masses)
    assert mass is not None
    assert mass > 0


def test_digestion_produces_fragments():
    base_masses = load_base_masses(REPO_ROOT / "data" / "base_masses.yaml")
    config = load_config(REPO_ROOT / "config.yaml")
    config.sequence["sequence"] = "ACGUACGUACGU"
    fragments = digest_sequence(
        target_id="test",
        sequence="ACGUACGUACGU",
        position_map={i: i for i in range(1, 13)},
        config=config,
        base_masses=base_masses,
    )
    assert len(fragments) > 0
    assert all(f.unmodified_mass > 0 for f in fragments)


def test_modifications_yaml_loads_and_searches():
    modifications = load_modifications(REPO_ROOT / "data" / "modifications.yaml")
    assert len(modifications) > 50  # ~118 entries as of the 2026-09-02 vendoring
    # m1A has mass_shift_from_unmodified ~14.0156 per the design spec's YAML excerpt
    matches = find_modifications_by_mass_shift(modifications, 14.0156, tolerance_da=0.001)
    assert any(m.id == "m1A" for m in matches)


def test_ms2_extraction_functions_are_importable_and_return_empty_on_no_fragments():
    # generate_theoretical_ms2_ions with an empty fragment list should just
    # return an empty list rather than raising — confirms the trimmed
    # ms2_extraction.py wiring (masses/models/ms1_mapping) is intact.
    config = load_config(REPO_ROOT / "config.yaml")
    base_masses = load_base_masses(REPO_ROOT / "data" / "base_masses.yaml")
    ions = generate_theoretical_ms2_ions([], config, base_masses)
    assert ions == []
    # extract_ms2_spectra is exercised against a real mzML file in the
    # (not-yet-written) integration test once one is available; here we
    # only check the import + call shape against a nonexistent path fails
    # gracefully (via the warnings list) rather than raising.
    warnings: list[dict] = []
    spectra = extract_ms2_spectra("__does_not_exist__.mzML", config.ms2_annotation, warnings=warnings)
    assert spectra == []
    assert any(w["Level"] == "WARNING" for w in warnings)
