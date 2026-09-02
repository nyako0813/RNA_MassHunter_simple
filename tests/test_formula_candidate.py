import pytest

from rna_masshunter.formula_candidate import (
    enumerate_formula_candidates,
    get_recommended_formula,
    rank_formula_candidates,
)

# Same provisional per-element caps as config.yaml's formula_candidate.element_limits
# (docs/design/RNA_MassHunter_再設計_実装仕様書.md §9.1). Without these, an
# unrestricted 7-atom search over C/H/N/O/P/S/Se can turn up a coincidental
# heavy-atom combination (e.g. involving P2/S1/Se-1) whose mass happens to be
# a hair closer to a small ΔDa than the "obvious" answer — a real property of
# an unbiased mass-only search, not a bug, but it makes the "obvious" answer
# an unstable choice for a targeted unit test unless the search space is
# bounded the way the pipeline actually configures it.
_PRACTICAL_LIMITS = {"C": 3, "H": 7, "N": 3, "O": 4, "P": 1, "S": 1, "Se": 1}


def test_single_oxygen_delta():
    candidates = enumerate_formula_candidates(15.9949, tolerance_da=0.001, element_limits=_PRACTICAL_LIMITS)
    assert get_recommended_formula(candidates).formula == "O1"


def test_ch2_delta():
    candidates = enumerate_formula_candidates(14.0157, tolerance_da=0.001, element_limits=_PRACTICAL_LIMITS)
    assert get_recommended_formula(candidates).formula in ("C1H2",)


@pytest.mark.parametrize("max_total_atoms,expected_atom_count_upper_bound", [(7, 7), (3, 3)])
def test_max_total_atoms_respected(max_total_atoms, expected_atom_count_upper_bound):
    candidates = enumerate_formula_candidates(50.0, tolerance_da=0.5, max_total_atoms=max_total_atoms)
    assert candidates
    assert all(c.atom_count <= expected_atom_count_upper_bound for c in candidates)


def test_ranking_not_biased_by_atom_count():
    # 3原子以上でも質量誤差が小さければ上位に来ることを明示的に検証（原則4）
    candidates = enumerate_formula_candidates(42.0106, tolerance_da=0.01)  # 例: acetylation相当
    assert candidates
    assert candidates == sorted(candidates, key=lambda c: (abs(c.mass_error_da), c.formula))
    assert not all(c.atom_count <= 2 for c in candidates)


def test_zero_delta_returns_empty_formula():
    candidates = enumerate_formula_candidates(0.0, tolerance_da=0.0001)
    assert get_recommended_formula(candidates).formula == "0"


def test_negative_delta_allows_loss_formula():
    candidates = enumerate_formula_candidates(-15.9949, tolerance_da=0.001, element_limits=_PRACTICAL_LIMITS)
    assert get_recommended_formula(candidates).formula == "O-1"


def test_no_candidates_within_tolerance_returns_empty_list():
    candidates = enumerate_formula_candidates(123.456, tolerance_da=0.0001, max_total_atoms=2)
    assert candidates == []
    assert get_recommended_formula(candidates) is None


def test_element_limits_restrict_search():
    # 3*O ~= 47.985 Da, within 0.5 Da tolerance of 48.0 and reachable only via O
    # under these element_limits (everything else capped at 0 atoms).
    candidates = enumerate_formula_candidates(
        48.0, tolerance_da=0.5, max_total_atoms=7, element_limits={"C": 0, "H": 0, "N": 0, "O": 4, "P": 0, "S": 0, "Se": 0},
    )
    assert candidates
    assert all(set(c.elements) <= {"O"} for c in candidates)


def test_max_candidates_truncates():
    candidates = enumerate_formula_candidates(50.0, tolerance_da=0.5, max_candidates=2)
    assert len(candidates) <= 2


def test_rank_formula_candidates_is_pure_function():
    candidates = enumerate_formula_candidates(42.0106, tolerance_da=0.01)
    shuffled = list(reversed(candidates))
    assert rank_formula_candidates(shuffled) == candidates
