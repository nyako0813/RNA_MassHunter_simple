import pytest
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES

def test_composition_addition(): assert (ElementalComposition({'C':1})+ElementalComposition({'H':2})).to_dict()=={'C':1,'H':2}
def test_composition_subtraction(): assert (ElementalComposition({'C':2})-ElementalComposition({'C':1})).to_dict()=={'C':1}
def test_negative_element_rejected():
    with pytest.raises(ValueError): ElementalComposition({'O':-1})
def test_subtraction_negative_rejected():
    with pytest.raises(ValueError): ElementalComposition({'O':1})-ElementalComposition({'O':2})
def test_exact_mass(): assert ElementalComposition({'C':1,'H':2}).exact_mass==pytest.approx(12+2*MONOISOTOPIC_ATOMIC_MASSES['H'])
def test_canonical_string(): assert ElementalComposition.delta({'S':1,'O':-1}).canonical_string()=='O-1S1'
def test_equality_hash_ignore_delta_construction_mode():
    a=ElementalComposition({'C':1}); b=ElementalComposition.delta({'C':1}); assert a==b and hash(a)==hash(b)
def test_o_to_s_delta(): assert ElementalComposition.delta({'O':-1,'S':1}).exact_mass==pytest.approx(15.97715658043)
def test_float_mass_comparison(): assert ElementalComposition({'H':2}).is_close_mass(ElementalComposition({'H':2}),1e-12)
def test_json_dict_order_is_deterministic(): assert list(ElementalComposition({'S':1,'C':2}).to_dict())==['C','S']
