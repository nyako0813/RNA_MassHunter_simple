"""Tests for rna_masshunter/hypothesis_check.py (hypothesis_mass_check_spec.md §8)."""
from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import openpyxl
import pytest

from rna_masshunter import simple_pipeline
from rna_masshunter.config import load_config
from rna_masshunter.hypothesis_check import (
    SOURCE_CONFIG_MANUAL,
    SOURCE_TRNA_LIBRARY_DEFAULT,
    Hypothesis,
    base_plus_elements_mass,
    build_label_mass_catalog,
    check_hypotheses,
    dinucleotide_mass,
    hypotheses_from_config_targets,
    hypotheses_from_conserved_modifications,
    hypothesis_from_target,
)
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Peak, RunConfig
from rna_masshunter.modifications import load_modifications
from rna_masshunter.trna_library import load_trna_library

REPO_ROOT = Path(__file__).resolve().parent.parent
REAL_MZML = REPO_ROOT / "data" / "input" / "05_Mix.mzML"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_simple_pipeline import _write_mzml  # noqa: E402  (shared synthetic-mzML writer)


@pytest.fixture(scope="module")
def modifications():
    return load_modifications(REPO_ROOT / "data" / "modifications.yaml")


@pytest.fixture(scope="module")
def catalog(modifications):
    return build_label_mass_catalog(modifications)


def _config(polarity: str = "positive", ppm: float = 15, max_charge: int = 4, targets=None) -> RunConfig:
    return RunConfig(
        instrument={"polarity": polarity},
        hypothesis_check={"enabled": True, "mass_tolerance_ppm": ppm, "max_charge": max_charge, "targets": targets or []},
    )


# --- §8-1: trna_library conserved_modifications --------------------------------

def test_all_58_trna_entries_have_conserved_modifications_including_archaeosine_pos15():
    library = load_trna_library(REPO_ROOT / "data" / "trna_library.yaml")
    assert len(library) == 58
    for trna_id, entry in library.items():
        mods = entry.get("conserved_modifications")
        assert mods, f"{trna_id} has no conserved_modifications"
        assert any(m["position"] == 15 and m.get("modification") == "G+" for m in mods), trna_id
        assert entry["sequence"][14] == "G", trna_id


# --- agmatidine (C+) at the tRNA-Ile(CAU) wobble (agmatidine_conserved_modification_spec.md) ----

ILE2_ID = "tRNA-Ile2-CAT-1-1"  # anticodon CAU = tRNA-Ile(CAU)


def test_only_ile2_cau_trna_carries_agmatidine_at_its_wobble_position():
    library = load_trna_library(REPO_ROOT / "data" / "trna_library.yaml")
    entry = library[ILE2_ID]
    mods = [(m["position"], m["modification"]) for m in entry["conserved_modifications"] if "modification" in m]
    assert mods == [(15, "G+"), (entry["wobble_position"], "C+")]
    assert entry["sequence"][entry["wobble_position"] - 1] == "C"
    assert entry["sequence"][entry["wobble_position"] - 1:][:3] == entry["anticodon"]
    # every other entry keeps just the archaeosine default
    for trna_id, other in library.items():
        if trna_id != ILE2_ID:
            assert [m["modification"] for m in other["conserved_modifications"] if "modification" in m] == ["G+"], trna_id


def test_agmatidine_catalog_mass_matches_cytidine_plus_agmatine_minus_water(catalog):
    from rna_masshunter.elemental_composition import ElementalComposition
    cytidine = ElementalComposition({"C": 9, "H": 13, "N": 3, "O": 5}).exact_mass
    agmatine = ElementalComposition({"C": 5, "H": 14, "N": 4}).exact_mass
    water = ElementalComposition({"H": 2, "O": 1}).exact_mass
    assert catalog["C+"] == pytest.approx(355.1968, abs=1e-4)
    assert catalog["C+"] == pytest.approx(cytidine + agmatine - water, abs=1e-3)


def test_ile2_archaeosine_and_agmatidine_are_checked_alongside_rule_candidates(catalog, modifications):
    library = load_trna_library(REPO_ROOT / "data" / "trna_library.yaml")
    warnings: list[dict] = []
    all_hypotheses = hypotheses_from_conserved_modifications(library[ILE2_ID], catalog, modifications, warnings)
    hypotheses = [h for h in all_hypotheses if h.name.endswith(("= G+", "= C+"))]
    assert [(h.name, round(h.theoretical_mass, 4), h.source) for h in hypotheses] == [
        (f"{ILE2_ID} position 15 = G+", 324.1182, SOURCE_TRNA_LIBRARY_DEFAULT),
        (f"{ILE2_ID} position 36 = C+", 355.1968, SOURCE_TRNA_LIBRARY_DEFAULT),
    ]
    assert [h.confidence for h in hypotheses] == ["high_probability", "confirmed"]
    assert len(all_hypotheses) > 2  # Ile(CAU) also gets the A37 (t6A family) candidates
    assert warnings == []  # every target matches the base at its position

    peaks = [
        Peak(mz=mz_from_neutral_mass(324.1182, 1, "positive"), intensity=100.0),
        Peak(mz=mz_from_neutral_mass(355.1968, 2, "positive"), intensity=200.0),
    ]
    rows = check_hypotheses(hypotheses, peaks, _config())
    assert [(r.hypothesis_name, r.match_found, r.charge) for r in rows] == [
        (f"{ILE2_ID} position 15 = G+", True, 1),
        (f"{ILE2_ID} position 36 = C+", True, 2),
    ]


# --- confidence tiers and position rules (conserved_modifications_confidence_tier_spec.md) ------

sys.path.insert(0, str(REPO_ROOT / "tools"))
import add_conserved_modifications as gen  # noqa: E402  (the reproducible generator for the library's rule entries)


@pytest.fixture(scope="module")
def library_entries():
    return list(load_trna_library(REPO_ROOT / "data" / "trna_library.yaml").values())


def test_discriminator_position_handles_with_and_without_cca():
    assert gen.discriminator_position("GGGCCA") == 3
    assert gen.discriminator_position("GGGCCC") == 6


def test_canonical_positions_reproduce_spec_validation(library_entries):
    # Anticodon 3rd base (standard 36) = U  ->  standard 37 is A, for all 16 such tRNAs.
    u36 = [e for e in library_entries if e["anticodon"][2] == "U"]
    assert len(u36) == 16
    assert all(gen.base_at(e, gen.canonical_positions(e)[37]) == "A" for e in u36)
    # tRNA-Phe(GAA): position 37 is G (wyosine substrate), not A.
    phe = [e for e in library_entries if e["id"].startswith("tRNA-Phe-GAA-")]
    assert len(phe) == 2
    assert all(gen.base_at(e, gen.canonical_positions(e)[37]) == "G" for e in phe)
    # Position 55 is U in 47 of 58 (81%).
    assert sum(gen.base_at(e, gen.canonical_positions(e)[55]) == "U" for e in library_entries) == 47


def test_wobble_position_matches_the_anticodon_in_every_entry(library_entries):
    mismatched = [e["id"] for e in library_entries if e["sequence"][e["wobble_position"] - 1:][:3] != e["anticodon"]]
    assert mismatched == []  # tRNA-Ala-TGC-1-1 was 37 (wrong) before ala_tgc_wobble_position_fix.md; now 34


def test_ala_tgc_wobble_fix_gives_it_the_u34_and_a37_rule_entries(library_entries):
    ala = next(e for e in library_entries if e["id"] == "tRNA-Ala-TGC-1-1")
    assert ala["wobble_position"] == 34
    assert ala["sequence"][33:36] == "UGC" == ala["anticodon"]
    rule_items = [(m["position"], m["modification_candidates"], m["confidence"]) for m in _candidate_entries(ala)]
    assert rule_items == [
        (34, ["cnm5U", "cmnm5U", "mnm5U", "mnm5s2U", "s2U", "s4U", "mcm5U", "mcm5s2U", "ncm5U", "ncm5Um"], "high_probability"),
        (37, ["t6A", "ms2t6A", "hn6A", "ms2hn6A"], "high_probability"),
    ]  # no U55 / A58: those count back from the discriminator, and this tRNA's T-loop bases are G55 / U58


def _candidate_entries(entry):
    return [m for m in entry["conserved_modifications"] if "modification_candidates" in m]


@pytest.mark.parametrize("first_candidate, expected_count", [
    ("Y", 47),        # U55            universal_U55_pseudouridine
    ("cnm5U", 15),    # wobble U34     ma_U34_main_target (14 + tRNA-Ala-TGC-1-1 after its wobble_position fix)
    ("t6A", 35),      # A37            ma_A37_t6A         (34 + tRNA-Ala-TGC-1-1)
    ("m1A", 47),      # A58            archaea_A58_methylation
    ("imG-14", 2),    # tRNA-Phe G37
])
def test_each_rule_applies_to_the_expected_number_of_entries(library_entries, first_candidate, expected_count):
    matching = [e for e in library_entries if any(m["modification_candidates"][0] == first_candidate for m in _candidate_entries(e))]
    assert len(matching) == expected_count


def test_rule_entries_carry_the_spec_candidate_lists_positions_and_tier(library_entries):
    expected_candidates = {
        "Y": ["Y"],
        "cnm5U": ["cnm5U", "cmnm5U", "mnm5U", "mnm5s2U", "s2U", "s4U", "mcm5U", "mcm5s2U", "ncm5U", "ncm5Um"],
        "t6A": ["t6A", "ms2t6A", "hn6A", "ms2hn6A"],
        "m1A": ["m1A", "m6A"],  # generic "methylation" (no defined mass) intentionally omitted
        "imG-14": ["imG-14", "imG", "imG2"],
    }
    for entry in library_entries:
        positions = gen.canonical_positions(entry)
        for item in _candidate_entries(entry):
            assert item["confidence"] == "high_probability"
            assert item["modification_candidates"] == expected_candidates[item["modification_candidates"][0]]
            assert item["position"] in positions.values()
    phe = next(e for e in library_entries if e["id"] == "tRNA-Phe-GAA-1-1")
    assert [(m["position"], m["modification_candidates"][0]) for m in _candidate_entries(phe) if m["modification_candidates"][0] == "imG-14"] == [(phe["wobble_position"] + 3, "imG-14")]


def test_no_rule_entries_where_the_base_does_not_match(library_entries):
    for entry in library_entries:
        positions = gen.canonical_positions(entry)
        for item in _candidate_entries(entry):
            base = gen.base_at(entry, item["position"])
            first = item["modification_candidates"][0]
            assert base == {"Y": "U", "cnm5U": "U", "t6A": "A", "m1A": "A", "imG-14": "G"}[first], entry["id"]
    # wobble C (Ile2) gets no U34 rule.
    ile2 = next(e for e in library_entries if e["id"] == ILE2_ID)
    assert not any(m["modification_candidates"][0] == "cnm5U" for m in _candidate_entries(ile2))


def test_a37_note_claims_the_anticodon_u_correlation_only_for_anticodon_3rd_base_u(library_entries):
    claim = "anticodon 3rd base = U correlation"
    with_claim = without_claim = 0
    for entry in library_entries:
        for item in _candidate_entries(entry):
            if item["modification_candidates"][0] != "t6A":
                continue
            assert ("ma_A37_t6A" in item["note"])
            if entry["anticodon"][2] == "U":
                assert claim in item["note"], entry["id"]
                with_claim += 1
            else:
                assert claim not in item["note"], entry["id"]
                without_claim += 1
    assert (with_claim, without_claim) == (16, 19)  # 16 anticodon-U tRNAs (all have A37) + 19 other A37 tRNAs
    ala = next(e for e in library_entries if e["id"] == "tRNA-Ala-TGC-1-1")
    assert ala["anticodon"][2] == "C"


def test_confidence_tiers_of_curated_entries(library_entries):
    for entry in library_entries:
        for item in entry["conserved_modifications"]:
            assert item["confidence"] in ("confirmed", "high_probability")
            if item.get("modification") == "G+":
                assert item["confidence"] == "high_probability" and "occasional absence" in item["note"]
            if item.get("modification") == "C+":
                assert item["confidence"] == "confirmed" and entry["id"] == ILE2_ID


def test_library_yaml_is_in_sync_with_the_generator():
    """Regenerating from the rules must reproduce data/trna_library.yaml exactly (hand edits belong in the generator)."""
    assert gen.main(["--check"]) == 0


def test_modification_candidates_expand_to_independent_hypotheses_sharing_position_and_confidence(catalog, modifications):
    entry = {
        "id": "t", "sequence": "A" * 40,
        "conserved_modifications": [{"position": 37, "modification_candidates": ["t6A", "ms2t6A", "hn6A", "ms2hn6A"], "confidence": "high_probability", "note": "n"}],
    }
    warnings: list[dict] = []
    hypotheses = hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings)
    assert [h.name for h in hypotheses] == [f"t position 37 = {label}" for label in ("t6A", "ms2t6A", "hn6A", "ms2hn6A")]
    assert [round(h.theoretical_mass, 4) for h in hypotheses] == [round(catalog[l], 4) for l in ("t6A", "ms2t6A", "hn6A", "ms2hn6A")]
    assert {h.confidence for h in hypotheses} == {"high_probability"}
    assert warnings == []

    # Each is matched independently: a peak at ms2t6A's mass only satisfies that candidate.
    peaks = [Peak(mz=mz_from_neutral_mass(catalog["ms2t6A"], 1, "positive"), intensity=10.0)]
    rows = check_hypotheses(hypotheses, peaks, _config())
    assert {r.hypothesis_name: r.match_found for r in rows} == {
        "t position 37 = t6A": False, "t position 37 = ms2t6A": True,
        "t position 37 = hn6A": False, "t position 37 = ms2hn6A": False,
    }
    assert {r.confidence for r in rows} == {"high_probability"}  # tier carried onto both Yes and No rows


def test_modification_and_candidates_in_one_entry_are_both_used_and_unknown_candidate_is_skipped(catalog, modifications):
    entry = {"id": "t", "sequence": "A" * 40, "conserved_modifications": [
        {"position": 37, "modification": "t6A", "modification_candidates": ["m1A", "methylation"], "confidence": "confirmed"},
    ]}
    warnings: list[dict] = []
    hypotheses = hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings)
    assert [h.name for h in hypotheses] == ["t position 37 = t6A", "t position 37 = m1A"]
    assert any(w["Level"] == "ERROR" and "methylation" in w["Message"] for w in warnings)


def test_unknown_confidence_value_warns(catalog, modifications):
    entry = {"id": "t", "sequence": "A" * 40, "conserved_modifications": [{"position": 37, "modification": "t6A", "confidence": "sure"}]}
    warnings: list[dict] = []
    hypotheses = hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings)
    assert len(hypotheses) == 1
    assert any("unknown confidence" in w["Message"] for w in warnings)


def test_every_rule_derived_label_exists_in_the_catalog_and_targets_the_right_base(library_entries, catalog, modifications):
    """Loading each tRNA's defaults must be warning-free: labels exist and target the base at their position."""
    for entry in library_entries:
        warnings: list[dict] = []
        hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings)
        assert warnings == [], (entry["id"], warnings)


# --- §8-2/§8-3: theoretical masses ---------------------------------------------

def test_dinucleotide_m22g_phosphorothioate_u_reproduces_real_data_mass(catalog):
    mass = dinucleotide_mass(catalog["m2,2G"], catalog["U"], "phosphorothioate")
    assert mass == pytest.approx(633.1254, abs=0.001)


def test_dinucleotide_m22g_phosphorothioate_c_mass(catalog):
    assert dinucleotide_mass(catalog["m2,2G"], catalog["C"], "phosphorothioate") == pytest.approx(632.1414, abs=0.001)


def test_dinucleotide_phosphodiester_matches_literature_ApA(catalog):
    # ApA (adenylyl-(3'→5')-adenosine), C20H26N10O10P: monoisotopic neutral 596.1493.
    assert dinucleotide_mass(catalog["A"], catalog["A"], "phosphodiester") == pytest.approx(596.1493, abs=0.001)


def test_phosphorothioate_is_S_minus_O_heavier_than_phosphodiester(catalog):
    normal = dinucleotide_mass(catalog["G"], catalog["U"], "phosphodiester")
    pt = dinucleotide_mass(catalog["G"], catalog["U"], "phosphorothioate")
    assert pt - normal == pytest.approx(31.97207 - 15.99491, abs=1e-4)


def test_dinucleotide_rejects_unknown_linkage(catalog):
    with pytest.raises(ValueError, match="linkage"):
        dinucleotide_mass(catalog["A"], catalog["U"], "phosphoramidate")


def test_base_plus_elements(catalog):
    assert base_plus_elements_mass(catalog["cnm5s2U"], {"O": 1, "S": 1}) == pytest.approx(347.0246, abs=0.001)
    assert base_plus_elements_mass(catalog["cnm5s2U"], {"S": 1}) == pytest.approx(331.0297, abs=0.001)


def test_base_plus_elements_rejects_unknown_element_and_non_integer_count(catalog):
    with pytest.raises(ValueError, match="Unsupported element"):
        base_plus_elements_mass(catalog["A"], {"Xx": 1})
    with pytest.raises(ValueError, match="integer"):
        base_plus_elements_mass(catalog["A"], {"O": 1.5})


def test_label_catalog_lookups(catalog):
    assert catalog["G+"] == pytest.approx(324.1182, abs=1e-4)
    assert catalog["G"] == pytest.approx(283.0917, abs=1e-4)
    assert catalog["G+"] - catalog["G"] == pytest.approx(41.0265, abs=1e-4)  # spec §2: archaeosine vs unmodified G
    assert catalog["ncm5s2U"] == pytest.approx(317.0682, abs=1e-4)


# --- target parsing --------------------------------------------------------------

def test_hypothesis_from_target_three_forms(catalog):
    label = hypothesis_from_target({"name": "G15 unmod", "label": "G"}, catalog)
    dinuc = hypothesis_from_target({"name": "d", "components": ["m2,2G", "U"], "linkage": "phosphorothioate"}, catalog)
    adduct = hypothesis_from_target({"name": "a", "base": "cnm5s2U", "add_elements": {"O": 1, "S": 1}}, catalog)
    assert (label.name, label.source) == ("G15 unmod", SOURCE_CONFIG_MANUAL)
    assert label.theoretical_mass == pytest.approx(283.0917, abs=1e-4)
    assert dinuc.theoretical_mass == pytest.approx(633.1254, abs=0.001)
    assert adduct.theoretical_mass == pytest.approx(347.0246, abs=0.001)


def test_dinucleotide_linkage_defaults_to_phosphodiester(catalog):
    hypothesis = hypothesis_from_target({"components": ["A", "A"]}, catalog)
    assert hypothesis.theoretical_mass == pytest.approx(596.1493, abs=0.001)


@pytest.mark.parametrize("target, message", [
    ({"name": "x", "label": "NOT_A_LABEL"}, "not found"),
    ({"name": "x", "label": "G", "base": "G", "add_elements": {"O": 1}}, "exactly one"),
    ({"name": "x"}, "exactly one"),
    ({"name": "x", "components": ["A"]}, "exactly 2"),
    ({"name": "x", "components": ["A", "NOT_A_LABEL"]}, "not found"),
    ({"name": "x", "components": ["A", "U"], "linkage": "bogus"}, "linkage"),
    ({"name": "x", "base": "A"}, "add_elements"),
    ({"name": "x", "base": "A", "add_elements": {"Xx": 1}}, "Unsupported element"),
])
def test_malformed_targets_raise_value_error(catalog, target, message):
    with pytest.raises(ValueError, match=message):
        hypothesis_from_target(target, catalog)


def test_hypotheses_from_config_targets_empty_by_default(catalog):
    assert hypotheses_from_config_targets(_config(), catalog) == []
    config = _config()
    config.hypothesis_check["targets"] = None  # `targets:` left blank in YAML
    assert hypotheses_from_config_targets(config, catalog) == []


# --- conserved_modifications -> hypotheses ---------------------------------------

def test_conserved_modifications_become_library_default_hypotheses(catalog, modifications):
    library = load_trna_library(REPO_ROOT / "data" / "trna_library.yaml")
    warnings: list[dict] = []
    hypotheses = hypotheses_from_conserved_modifications(library["tRNA-Gln-TTG-2-1"], catalog, modifications, warnings)
    assert all(h.source == SOURCE_TRNA_LIBRARY_DEFAULT for h in hypotheses)
    archaeosine = [h for h in hypotheses if h.name.endswith("= G+")]
    assert len(archaeosine) == 1
    assert archaeosine[0].theoretical_mass == pytest.approx(324.1182, abs=1e-4)
    assert "position 15" in archaeosine[0].name
    assert warnings == []  # G+ targets G and position 15 is G


def test_conserved_modification_unknown_label_is_skipped_with_error(catalog, modifications):
    entry = {"id": "t", "sequence": "G" * 20, "conserved_modifications": [{"position": 15, "modification": "NOPE"}]}
    warnings: list[dict] = []
    assert hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings) == []
    assert any(w["Level"] == "ERROR" and "NOPE" in w["Message"] for w in warnings)


def test_conserved_modification_base_mismatch_warns_but_still_checks(catalog, modifications):
    entry = {"id": "t", "sequence": "A" * 20, "conserved_modifications": [{"position": 15, "modification": "G+"}]}
    warnings: list[dict] = []
    hypotheses = hypotheses_from_conserved_modifications(entry, catalog, modifications, warnings)
    assert len(hypotheses) == 1
    assert any(w["Level"] == "WARNING" and "position 15 is A" in w["Message"] for w in warnings)


def test_no_trna_entry_gives_no_hypotheses(catalog, modifications):
    assert hypotheses_from_conserved_modifications(None, catalog, modifications) == []


# --- matching (§5) ------------------------------------------------------------------

def _hyp(mass: float) -> Hypothesis:
    return Hypothesis("h", SOURCE_CONFIG_MANUAL, "desc", mass)


@pytest.mark.parametrize("polarity", ["positive", "negative"])
@pytest.mark.parametrize("charge", [1, 2, 3, 4])
def test_exact_match_at_each_charge_and_polarity(polarity, charge):
    mass = 633.1254
    peaks = [Peak(mz=mz_from_neutral_mass(mass, charge, polarity), intensity=500.0, rt=3.0, scan_id="s")]
    rows = check_hypotheses([_hyp(mass)], peaks, _config(polarity=polarity))
    matched = [r for r in rows if r.match_found]
    assert [(r.charge, r.peak_id) for r in matched] == [(charge, "PK0001")]
    assert matched[0].delta_ppm == pytest.approx(0.0, abs=1e-6)
    assert matched[0].observed_mass == pytest.approx(mass, abs=1e-9)


def test_charge_above_max_charge_is_not_searched():
    mass = 633.1254
    peaks = [Peak(mz=mz_from_neutral_mass(mass, 3, "positive"), intensity=500.0)]
    assert not any(r.match_found for r in check_hypotheses([_hyp(mass)], peaks, _config(max_charge=2)))
    assert any(r.match_found for r in check_hypotheses([_hyp(mass)], peaks, _config(max_charge=3)))


def test_tolerance_boundary_is_on_neutral_mass_ppm():
    mass = 500.0
    inside = Peak(mz=mz_from_neutral_mass(mass * (1 + 14e-6), 1, "positive"), intensity=1.0)
    outside = Peak(mz=mz_from_neutral_mass(mass * (1 + 16e-6), 1, "positive"), intensity=1.0)
    rows = check_hypotheses([_hyp(mass)], [inside, outside], _config(ppm=15))
    assert [r.peak_id for r in rows if r.match_found] == ["PK0001"]
    assert rows[0].delta_ppm == pytest.approx(14.0, abs=1e-3)
    assert rows[0].delta_da == pytest.approx(mass * 14e-6, abs=1e-6)


def test_no_match_yields_single_no_row_with_blank_fields():
    rows = check_hypotheses([_hyp(633.1254)], [Peak(mz=100.0, intensity=1.0)], _config())
    assert len(rows) == 1
    row = rows[0]
    assert row.match_found is False
    assert (row.peak_id, row.charge, row.observed_mass, row.delta_da, row.delta_ppm, row.intensity) == (None,) * 6
    assert row.theoretical_mass == 633.1254


def test_multiple_peaks_expand_to_multiple_rows_and_peak_ids_follow_input_order():
    mass = 400.0
    mzs = [mz_from_neutral_mass(mass * (1 + ppm * 1e-6), 1, "positive") for ppm in (5, 0, -3)]
    peaks = [Peak(mz=mz, intensity=10.0) for mz in mzs]
    rows = check_hypotheses([_hyp(mass)], peaks, _config())
    assert sorted(r.peak_id for r in rows) == ["PK0001", "PK0002", "PK0003"]
    assert [r.peak_id for r in rows] == ["PK0002", "PK0003", "PK0001"]  # sorted by |Δppm|


def test_no_peaks_gives_no_rows_per_hypothesis():
    rows = check_hypotheses([_hyp(1.0e2), _hyp(2.0e2)], [], _config())
    assert [r.match_found for r in rows] == [False, False]


# --- config validation ----------------------------------------------------------------

def test_default_config_has_empty_targets():
    from rna_masshunter.config import DEFAULT_CONFIG
    hc = DEFAULT_CONFIG["hypothesis_check"]
    assert hc["targets"] == [] and hc["max_charge"] == 4 and hc["mass_tolerance_ppm"] == 15


@pytest.mark.parametrize("override, message", [
    ({"mass_tolerance_ppm": 0}, "mass_tolerance_ppm"),
    ({"max_charge": 0}, "max_charge"),
    ({"max_charge": 2.5}, "max_charge"),
    ({"targets": "oops"}, "targets"),
])
def test_validate_config_rejects_bad_hypothesis_check(tmp_path, override, message):
    from rna_masshunter.config import validate_config
    path = tmp_path / "c.yaml"
    path.write_text("instrument:\n  polarity: positive\ndigestion:\n  enzyme: Nuclease_P1\n", encoding="utf-8")
    config = load_config(path)
    config.hypothesis_check.update(override)
    with pytest.raises(ValueError, match=message):
        validate_config(config)


# --- §8-4: pipeline integration (synthetic mzML) ----------------------------------------------

def _pipeline_config(tmp_path: Path, mzml_path, hypothesis_yaml: str, trna_type: str = "", enzyme: str = "Nuclease_P1") -> Path:
    text = textwrap.dedent(f"""
        sequence:
          trna_type: "{trna_type}"
        instrument:
          polarity: positive
        digestion:
          enzyme: {enzyme}
        alkaline_phosphatase:
          enabled: true
        input:
          mzml_path: "{mzml_path}"
        project:
          output_dir: {tmp_path / "output"}
        ms1_peak_extraction:
          mz_min: 0
          mz_max: 5000
          intensity_threshold: 0
        reporting:
          excel_output: true
    """) + textwrap.dedent(hypothesis_yaml)
    path = tmp_path / "config.yaml"
    path.write_text(text, encoding="utf-8")
    return path


_HYPOTHESIS_TARGETS_YAML = """
    hypothesis_check:
      enabled: true
      mass_tolerance_ppm: 15
      max_charge: 4
      targets:
        - name: "m2,2G-PT-U dinucleotide"
          components: ["m2,2G", "U"]
          linkage: phosphorothioate
        - name: "m2,2G-PT-C dinucleotide"
          components: ["m2,2G", "C"]
          linkage: phosphorothioate
        - name: "cnm5s2U+O+S"
          base: cnm5s2U
          add_elements: {O: 1, S: 1}
"""


def test_pipeline_writes_07b_sheet_with_yes_and_no_rows_and_merges_trna_defaults(tmp_path, catalog):
    pt_u = dinucleotide_mass(catalog["m2,2G"], catalog["U"], "phosphorothioate")
    archaeosine = catalog["G+"]
    mzml_path = tmp_path / "hyp.mzML"
    _write_mzml(mzml_path, [
        {"id": "scan=1", "ms_level": 1, "rt": 1.0, "mzs": [mz_from_neutral_mass(pt_u, 1, "positive")], "intensities": [1000.0]},
        {"id": "scan=2", "ms_level": 1, "rt": 2.0, "mzs": [mz_from_neutral_mass(archaeosine, 2, "positive")], "intensities": [800.0]},
    ])
    config_path = _pipeline_config(tmp_path, mzml_path, _HYPOTHESIS_TARGETS_YAML, trna_type="tRNA-Gln-TTG-2-1")

    result = simple_pipeline.run(config_path, project_root=REPO_ROOT)

    by_name = {}
    for row in result["hypothesis_check_rows"]:
        by_name.setdefault(row.hypothesis_name, []).append(row)
    names = list(by_name)
    assert names[0] == "tRNA-Gln-TTG-2-1 position 15 = G+"  # library defaults first ...
    assert names[-3:] == ["m2,2G-PT-U dinucleotide", "m2,2G-PT-C dinucleotide", "cnm5s2U+O+S"]  # ... then config targets
    assert all(n.startswith("tRNA-Gln-TTG-2-1 position ") for n in names[:-3])
    assert {row.confidence for name in names[-3:] for row in by_name[name]} == {""}
    assert by_name["tRNA-Gln-TTG-2-1 position 15 = G+"][0].source == SOURCE_TRNA_LIBRARY_DEFAULT
    assert by_name["tRNA-Gln-TTG-2-1 position 15 = G+"][0].charge == 2
    assert by_name["tRNA-Gln-TTG-2-1 position 15 = G+"][0].confidence == "high_probability"
    # Peak IDs follow the pipeline's (m/z-sorted) peak list, same as every other sheet.
    pt_u_mz = mz_from_neutral_mass(pt_u, 1, "positive")
    pt_u_id = f"PK{next(i for i, p in enumerate(result['peaks']) if p.mz == pytest.approx(pt_u_mz)) + 1:04d}"
    assert by_name["m2,2G-PT-U dinucleotide"][0].match_found and by_name["m2,2G-PT-U dinucleotide"][0].peak_id == pt_u_id
    assert not by_name["m2,2G-PT-C dinucleotide"][0].match_found
    assert not by_name["cnm5s2U+O+S"][0].match_found

    wb = openpyxl.load_workbook(result["output_path"])
    assert wb.sheetnames == [
        "01_Index", "02_Input", "03_Nucleoside_Targets", "04_Observed_Mass", "05_Mass_Intensity",
        "06_Nucleoside_Comparison", "07_Modifications", "07b_Hypothesis_Check", "08_Visualization",
    ]
    ws = wb["07b_Hypothesis_Check"]
    assert [c.value for c in ws[3]] == [
        "Hypothesis_Name", "Source", "Confidence", "Formula_Description", "Theoretical_Mass", "Match_Found",
        "Matched_Peak_IDs", "Charge", "Observed_Mass", "ΔDa", "Δppm", "Intensity", "RT_Range",
    ]
    rows = {r[0].value: [c.value for c in r] for r in ws.iter_rows(min_row=4)}
    assert rows["m2,2G-PT-U dinucleotide"][5] == "Yes" and rows["m2,2G-PT-U dinucleotide"][6] == pt_u_id
    assert rows["m2,2G-PT-C dinucleotide"][5] == "No"
    assert rows["m2,2G-PT-U dinucleotide"][2] in (None, "")  # config targets carry no confidence
    assert rows["tRNA-Gln-TTG-2-1 position 15 = G+"][2] == "high_probability"
    assert wb["01_Index"]["A8"].value == "07b_Hypothesis_Check"


def test_pipeline_oligomer_mode_also_writes_07b(tmp_path, catalog):
    mzml_path = tmp_path / "hyp.mzML"
    mass = catalog["G"]
    _write_mzml(mzml_path, [{"id": "scan=1", "ms_level": 1, "rt": 1.0, "mzs": [mz_from_neutral_mass(mass, 1, "positive")], "intensities": [1000.0]}])
    yaml_text = """
        hypothesis_check:
          targets:
            - {name: "G", label: G}
    """
    config_path = _pipeline_config(tmp_path, mzml_path, yaml_text, enzyme="RNase_T1")

    result = simple_pipeline.run(config_path, project_root=REPO_ROOT)

    assert result["mode"] == "oligomer"
    assert [r.match_found for r in result["hypothesis_check_rows"]] == [True]
    assert "07b_Hypothesis_Check" in openpyxl.load_workbook(result["output_path"]).sheetnames


def test_pipeline_without_hypotheses_omits_07b_and_leaves_existing_output_unchanged(tmp_path):
    config_path = _pipeline_config(tmp_path, "", "")
    result = simple_pipeline.run(config_path, project_root=REPO_ROOT)
    assert result["hypotheses"] == [] and result["hypothesis_check_rows"] == []
    assert "07b_Hypothesis_Check" not in openpyxl.load_workbook(result["output_path"]).sheetnames


def test_pipeline_disabled_hypothesis_check_omits_07b(tmp_path):
    yaml_text = """
        hypothesis_check:
          enabled: false
          targets:
            - {name: "G", label: G}
    """
    result = simple_pipeline.run(_pipeline_config(tmp_path, "", yaml_text, trna_type="tRNA-Gln-TTG-2-1"), project_root=REPO_ROOT)
    assert result["hypotheses"] == []
    assert "07b_Hypothesis_Check" not in openpyxl.load_workbook(result["output_path"]).sheetnames


def test_pipeline_without_mzml_skips_check_with_warning(tmp_path):
    yaml_text = """
        hypothesis_check:
          targets:
            - {name: "G", label: G}
    """
    result = simple_pipeline.run(_pipeline_config(tmp_path, "", yaml_text), project_root=REPO_ROOT)
    assert len(result["hypotheses"]) == 1 and result["hypothesis_check_rows"] == []
    assert any("hypothesis_check was skipped" in w["Message"] for w in result["warnings"])


def test_pipeline_unknown_config_label_fails_fast(tmp_path):
    yaml_text = """
        hypothesis_check:
          targets:
            - {name: "typo", label: NOT_A_LABEL}
    """
    with pytest.raises(ValueError, match="NOT_A_LABEL"):
        simple_pipeline.run(_pipeline_config(tmp_path, "", yaml_text), project_root=REPO_ROOT)


def test_pipeline_ile2_trna_type_checks_archaeosine_and_agmatidine(tmp_path, catalog):
    mzml_path = tmp_path / "ile2.mzML"
    _write_mzml(mzml_path, [
        {"id": "scan=1", "ms_level": 1, "rt": 1.0, "mzs": [mz_from_neutral_mass(catalog["C+"], 1, "positive")], "intensities": [1000.0]},
    ])
    result = simple_pipeline.run(_pipeline_config(tmp_path, mzml_path, "", trna_type=ILE2_ID), project_root=REPO_ROOT)

    found = {r.hypothesis_name: (r.match_found, r.confidence) for r in result["hypothesis_check_rows"]}
    assert found[f"{ILE2_ID} position 15 = G+"] == (False, "high_probability")
    assert found[f"{ILE2_ID} position 36 = C+"] == (True, "confirmed")
    # Excel: Confidence column shows archaeosine as high_probability and agmatidine as confirmed.
    ws = openpyxl.load_workbook(result["output_path"])["07b_Hypothesis_Check"]
    header = [c.value for c in ws[3]]
    confidence_by_name = {r[header.index("Hypothesis_Name")].value: r[header.index("Confidence")].value for r in ws.iter_rows(min_row=4)}
    assert confidence_by_name[f"{ILE2_ID} position 15 = G+"] == "high_probability"
    assert confidence_by_name[f"{ILE2_ID} position 36 = C+"] == "confirmed"


# --- §8-4: real-data verification (needs the git-ignored raw mzML) ------------------------------

@pytest.mark.skipif(not REAL_MZML.exists(), reason="data/input/05_Mix.mzML is git-ignored raw data")
def test_real_data_05_mix_m22g_pt_u_matches_pk72358_to_72360(tmp_path):
    config_path = _pipeline_config(tmp_path, REAL_MZML, _HYPOTHESIS_TARGETS_YAML)
    text = config_path.read_text().replace("mz_min: 0", "mz_min: 120").replace("mz_max: 5000", "mz_max: 3000").replace("intensity_threshold: 0", "intensity_threshold: 1000")
    config_path.write_text(text, encoding="utf-8")

    result = simple_pipeline.run(config_path, project_root=REPO_ROOT)

    rows = result["hypothesis_check_rows"]
    pt_u = [r for r in rows if r.hypothesis_name == "m2,2G-PT-U dinucleotide"]
    assert all(r.match_found for r in pt_u)
    assert {"PK72358", "PK72359", "PK72360"} <= {r.peak_id for r in pt_u}
    assert all(r.charge == 1 and abs(r.delta_ppm) < 2 for r in pt_u)
    assert pt_u[0].theoretical_mass == pytest.approx(633.1254, abs=0.001)
    # §2: no match for the other hypotheses in this standards mix (C variant is ~26 ppm off).
    assert not any(r.match_found for r in rows if r.hypothesis_name in ("m2,2G-PT-C dinucleotide", "cnm5s2U+O+S"))

    ws = openpyxl.load_workbook(result["output_path"])["07b_Hypothesis_Check"]
    header = [c.value for c in ws[3]]
    name_col, match_col, peak_col = (header.index(n) for n in ("Hypothesis_Name", "Match_Found", "Matched_Peak_IDs"))
    yes_ids = {r[peak_col].value for r in ws.iter_rows(min_row=4) if r[name_col].value == "m2,2G-PT-U dinucleotide" and r[match_col].value == "Yes"}
    assert {"PK72358", "PK72359", "PK72360"} <= yes_ids
