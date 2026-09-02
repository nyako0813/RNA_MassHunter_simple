"""Data structures for the RNA_MassHunter simple pipeline.

This is a trimmed-down copy of the dataclasses needed by the vendored
modules (masses/digestion/ms1_mapping/modifications/etc.) and by the new
simple-pipeline modules (formula_candidate/mass_comparison/ms2_support).

Origin: extracted from nyako0813/RNA_MassHunter `rna_masshunter/models.py`
(only the fields/classes actually used here are kept; many unrelated
SCIEX/shadow-audit dataclasses from the original file were intentionally
dropped). See docs/design/RNA_MassHunter_再設計_実装仕様書.md §8, §22.
"""
from dataclasses import dataclass, field
from typing import Any


@dataclass
class RunConfig:
    """Flat, dict-of-dicts config object, mirroring the pattern used by the
    original repository's config.py (each section is a plain dict, not a
    nested dataclass, so `.get()` with defaults works the same way)."""

    project: dict[str, Any] = field(default_factory=dict)
    input: dict[str, Any] = field(default_factory=dict)
    sequence: dict[str, Any] = field(default_factory=dict)
    instrument: dict[str, Any] = field(default_factory=dict)
    cca_processing: dict[str, Any] = field(default_factory=dict)
    digestion: dict[str, Any] = field(default_factory=dict)
    alkaline_phosphatase: dict[str, Any] = field(default_factory=dict)
    fragment_mapping: dict[str, Any] = field(default_factory=dict)
    ms1_peak_extraction: dict[str, Any] = field(default_factory=dict)
    ms2_annotation: dict[str, Any] = field(default_factory=dict)
    formula_candidate: dict[str, Any] = field(default_factory=dict)
    mass_comparison: dict[str, Any] = field(default_factory=dict)
    reporting: dict[str, Any] = field(default_factory=dict)
    visualization: dict[str, Any] = field(default_factory=dict)
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class Modification:
    id: str
    symbol: str | None
    mass_shift_from_unmodified: float
    category: str
    target_bases: list[str]
    detectability: Any = None
    curation: Any = None
    sources: Any = None
    source: Any = None
    source_priority: Any = None
    curation_status: str = ""
    candidate_policy: dict[str, Any] = field(default_factory=dict)
    chemical_group: str = ""
    near_isobaric_group: str = ""
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class Peak:
    mz: float
    intensity: float
    rt: float | None = None
    scan_id: str | None = None
    ms_level: int = 1
    tier: str | None = None


@dataclass
class Fragment:
    fragment_id: str
    target_id: str
    sequence: str
    start: int
    end: int
    standard_start: int | None
    standard_end: int | None
    enzyme: str
    missed_cleavages: int
    terminal_form: str
    unmodified_mass: float
    warnings: list[str] = field(default_factory=list)


@dataclass
class FragmentMS1Match:
    match_id: str
    fragment_id: str
    target_id: str
    sequence: str
    start: int
    end: int
    standard_start: int | None
    standard_end: int | None
    enzyme: str
    missed_cleavages: int
    terminal_form: str
    fragment_mass: float
    charge: int
    theoretical_mz: float
    observed_mz: float
    mass_error_da: float
    mass_error_ppm: float
    intensity: float
    rt: float | None
    scan_id: str | None


@dataclass
class MS2SpectrumInfo:
    spectrum_id: str
    scan_index: int
    rt: float | None
    precursor_mz: float | None
    precursor_charge: int | None
    precursor_intensity: float | None
    num_peaks: int
    base_peak_mz: float | None
    base_peak_intensity: float | None
    total_ion_current: float
    peaks: list[tuple[float, float]] = field(default_factory=list)
    raw_peaks: list[tuple[float, float]] | None = None
    scan_mz_min: float | None = None
    scan_mz_max: float | None = None
    effective_intensity_threshold: float | None = None
    threshold_information_available: bool = False


@dataclass
class TheoreticalMS2Ion:
    ion_id: str
    parent_fragment_id: str
    parent_sequence: str
    ion_type: str
    ion_sequence: str
    ion_start: int
    ion_end: int
    charge: int
    theoretical_mz: float
    theoretical_mass: float
    neutral_loss: str = ""
    modification_id: str = ""
    modification_name: str = ""
    comment: str = ""
    base_loss_position: int = 0
    base_loss_base: str = ""
