from pathlib import Path

import pytest

from rna_masshunter.masses import DEFAULT_PHOSPHATE_MASS, load_base_masses
from rna_masshunter.models import Modification
from rna_masshunter.nucleoside_targets import (
    build_modified_nucleoside_targets,
    build_nucleoside_target_universe,
    build_standard_nucleoside_targets,
    theoretical_mass_for_target,
)

REPO_ROOT = Path(__file__).resolve().parent.parent

# Published monoisotopic masses for the free ribonucleosides (Da).
_LITERATURE_MASSES = {"A": 267.0968, "C": 243.0855, "G": 283.0917, "U": 244.0695}


def test_standard_nucleoside_masses_match_literature_values():
    targets = build_standard_nucleoside_targets()
    by_label = {t.label: t for t in targets}

    assert set(by_label) == {"A", "C", "G", "U"}
    for label, expected_mass in _LITERATURE_MASSES.items():
        assert by_label[label].nucleoside_mass == pytest.approx(expected_mass, abs=1e-3)
        assert by_label[label].category == "standard"

    assert by_label["A"].name == "Adenosine"
    assert by_label["C"].name == "Cytidine"
    assert by_label["G"].name == "Guanosine"
    assert by_label["U"].name == "Uridine"


def test_build_modified_nucleoside_targets_from_real_modifications_yaml():
    from rna_masshunter.modifications import load_modifications

    modifications = load_modifications(REPO_ROOT / "data" / "modifications.yaml")
    targets = build_modified_nucleoside_targets(modifications)

    assert len(targets) == len(modifications)  # every real entry has modified_nucleoside_mass_mono
    by_label = {t.label: t for t in targets}
    assert by_label["m1A"].nucleoside_mass == pytest.approx(281.1124, abs=1e-4)
    assert by_label["m1A"].name == "1-methyladenosine"
    assert by_label["m1A"].category == "modified"


def test_missing_modified_nucleoside_mass_mono_is_skipped_with_warning():
    modifications = [
        Modification(id="hasmass", symbol="hm", mass_shift_from_unmodified=14.0, category="biological", target_bases=["A"], raw={"modified_nucleoside_mass_mono": 300.0}),
        Modification(id="nomass", symbol="nm", mass_shift_from_unmodified=14.0, category="biological", target_bases=["A"], raw={}),
    ]
    warnings: list[dict] = []

    targets = build_modified_nucleoside_targets(modifications, warnings)

    assert len(targets) == 1
    assert targets[0].label == "hm"
    assert any(w["Source"] == "nucleoside_targets" and "nomass" in str(w.get("Context")) for w in warnings)


def test_build_nucleoside_target_universe_combines_standard_and_modified():
    modifications = [
        Modification(id="m1A", symbol="m1A", mass_shift_from_unmodified=14.0156, category="biological", target_bases=["A"], raw={"modified_nucleoside_mass_mono": 281.1124}),
    ]

    universe = build_nucleoside_target_universe(modifications)

    assert len(universe) == 5  # 4 standard + 1 modified
    assert {t.category for t in universe} == {"standard", "modified"}


def test_theoretical_mass_for_target_dephosphorylated_vs_phosphorylated():
    target = build_standard_nucleoside_targets()[0]  # A, 267.0968
    base_masses = {"constants": {"phosphate": 79.966331}}

    free_mass = theoretical_mass_for_target(target, dephosphorylated=True, base_masses=base_masses)
    phospho_mass = theoretical_mass_for_target(target, dephosphorylated=False, base_masses=base_masses)

    assert free_mass == pytest.approx(target.nucleoside_mass, abs=1e-9)
    assert phospho_mass == pytest.approx(target.nucleoside_mass + 79.966331, abs=1e-9)


def test_theoretical_mass_for_target_falls_back_to_default_phosphate_mass():
    target = build_standard_nucleoside_targets()[0]

    phospho_mass = theoretical_mass_for_target(target, dephosphorylated=False, base_masses={})

    assert phospho_mass == pytest.approx(target.nucleoside_mass + DEFAULT_PHOSPHATE_MASS, abs=1e-9)


def test_theoretical_mass_for_target_uses_real_base_masses_yaml():
    base_masses = load_base_masses(REPO_ROOT / "data" / "base_masses.yaml")
    target = build_standard_nucleoside_targets()[0]

    phospho_mass = theoretical_mass_for_target(target, dephosphorylated=False, base_masses=base_masses)

    assert phospho_mass > target.nucleoside_mass
