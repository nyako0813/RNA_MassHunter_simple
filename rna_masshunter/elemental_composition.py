"""Deterministic elemental compositions for shadow chemical-state calculations."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Mapping
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES

_ELEMENT_ORDER=("C","H","N","O","P","S","Se")

@dataclass(frozen=True)
class ElementalComposition:
    _items: tuple[tuple[str,int], ...] = ()
    _signed: bool = field(default=False, compare=False, hash=False)

    def __init__(self, counts: Mapping[str,int] | None = None, *, allow_negative: bool = False):
        clean={}
        for element,value in (counts or {}).items():
            element=str(element); value=int(value)
            if element not in MONOISOTOPIC_ATOMIC_MASSES: raise ValueError(f"Unsupported element: {element}")
            if value < 0 and not allow_negative: raise ValueError(f"Negative element count: {element}={value}")
            if value: clean[element]=value
        order={e:i for i,e in enumerate(_ELEMENT_ORDER)}
        object.__setattr__(self,"_items",tuple(sorted(clean.items(),key=lambda x:(order.get(x[0],99),x[0]))))
        object.__setattr__(self,"_signed",bool(allow_negative or any(v<0 for v in clean.values())))

    @classmethod
    def delta(cls, counts: Mapping[str,int] | None = None): return cls(counts,allow_negative=True)
    def to_dict(self): return dict(self._items)
    def __add__(self,other):
        if not isinstance(other,ElementalComposition): return NotImplemented
        values=self.to_dict()
        for e,v in other._items: values[e]=values.get(e,0)+v
        return ElementalComposition(values,allow_negative=self._signed or other._signed)
    def __sub__(self,other):
        if not isinstance(other,ElementalComposition): return NotImplemented
        values=self.to_dict()
        for e,v in other._items: values[e]=values.get(e,0)-v
        return ElementalComposition(values,allow_negative=self._signed)
    @property
    def exact_mass(self): return sum(MONOISOTOPIC_ATOMIC_MASSES[e]*v for e,v in self._items)
    def canonical_string(self): return "".join(f"{e}{v}" for e,v in self._items) or "0"
    def is_close_mass(self,other,tolerance=1e-9): return abs(self.exact_mass-other.exact_mass)<=tolerance
    def __bool__(self): return bool(self._items)
