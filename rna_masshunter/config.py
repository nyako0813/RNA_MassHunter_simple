"""Config loading for the RNA_MassHunter simple pipeline.

Trimmed re-implementation of nyako0813/RNA_MassHunter `rna_masshunter/config.py`,
following the same pattern (`DEFAULT_CONFIG` dict merged over user YAML,
`RunConfig(**merged)`), but with only the sections this simplified pipeline
actually uses. The original file also has `sciex_profile`, `organism`,
`experiment`, `reconstruction` (used generically for MS1 peak extraction —
renamed here to `ms1_peak_extraction` since this project has no intact-mass
"reconstruction" step), `modification_search`, `unknown_modification_search`,
`p1_annotation`, `p1_sap_dinucleotide`, `modification_evidence_ranking`,
`biological_context`, `peak_filtering`, `performance`, `cca_tail` — none of
which the simple pipeline reads, so they were dropped. `formula_candidate`
and `mass_comparison` are new sections (see
docs/design/RNA_MassHunter_再設計_実装仕様書.md §9).
"""
from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.models import RunConfig
from rna_masshunter.warnings_manager import add_warning


DEFAULT_CONFIG: dict[str, dict[str, Any]] = {
    "project": {"name": "RNA_MassHunter_simple", "output_dir": "output", "log_dir": "logs"},
    "input": {"mzml_path": ""},
    "sequence": {
        "name": "target_tRNA",
        "type": "RNA",
        # Set to a data/trna_library.yaml id to auto-fill sequence/anticodon/
        # wobble_position instead of pasting them manually (see
        # rna_masshunter/trna_library.py).
        "trna_type": "",
        "sequence": "",
        "anticodon": "",
        "wobble_position": None,
    },
    "instrument": {"polarity": "negative"},

    "cca_processing": {
        "enabled": True,
        # "EXCLUDES_CCA" | "INCLUDES_COMPLETE_CCA" | "UNKNOWN".
        # Deliberately NOT auto-detected from the sequence string (see
        # rna_masshunter/cca_processing.py docstring) — always set explicitly.
        "registered_sequence_cca_mode": "EXCLUDES_CCA",
    },

    "digestion": {
        "enabled": True,
        "enzyme": "RNase_T1",
        "digestion_mode": "specific",
        "missed_cleavages": 1,
        "min_length": 2,
        "max_length": None,
        "include_terminal_forms": True,
        "allow_partial_digestion": True,
        "allow_nonspecific_cleavage": False,
    },
    "alkaline_phosphatase": {
        "enabled": False,
        "assume_complete": False,
        "allow_residual_phosphate": True,
        "allow_cyclic_phosphate": True,
    },

    "fragment_mapping": {
        "enabled": True,
        "mz_tolerance_ppm": 10,
        "min_charge": 1,
        "max_charge": 8,
        "polarity": "auto",
    },

    # Filters applied when pulling MS1 peaks out of the mzML file. Same
    # fields/semantics as the original repo's `reconstruction` section, but
    # renamed since this project has no intact-mass reconstruction step.
    "ms1_peak_extraction": {
        "rt_min": None,
        "rt_max": None,
        "mz_min": 500,
        "mz_max": 3000,
        "intensity_threshold": 1000,
        # §14B: collapse profile-mode split points (same scan, within
        # merge_tolerance_ppm of each other) down to one peak per ion before
        # anything downstream sees them. 10ppm comfortably covers the
        # ~0.00002-0.00005 Da split observed on real profile-mode data
        # without merging genuinely distinct nearby ions.
        "merge_profile_points": True,
        "merge_tolerance_ppm": 10,
        # §14C: after within-scan merging, additionally collapse the same
        # ion's repeated per-scan detections across a run of up to
        # merge_max_scan_gap consecutive scans (reuses merge_tolerance_ppm
        # for the m/z side). 2 is a provisional default — check it against
        # the real gap distribution between consecutive detections of the
        # same ion before trusting it unreviewed on a new dataset (§14C).
        "merge_across_scans": True,
        "merge_max_scan_gap": 2,
    },

    "ms2_annotation": {
        "enabled": True,
        "mz_tolerance_ppm": 20,
        "min_peak_intensity": 10,
        "min_relative_intensity_percent": 1.0,
        "max_peaks_per_spectrum": 500,
        "precursor_match_tolerance_ppm": 20,
        "use_theoretical_fragments": True,
        "ion_series": ["d", "w", "a", "z"],
        "charge_states": [1],
        "min_ion_length": 1,
        "max_ion_length": None,
    },

    "formula_candidate": {
        "enabled": True,
        "max_total_atoms": 7,
        "elements": ["C", "H", "N", "O", "P", "S", "Se"],
        # Provisional per-element caps — revisit once real ΔDa/modification
        # distributions from data/modifications.yaml have been reviewed
        # (see 仕様書 §9.1, §12, and 未確定事項リスト).
        "element_limits": {"C": 3, "H": 7, "N": 3, "O": 4, "P": 1, "S": 1, "Se": 1},
        "mass_tolerance_da": 0.01,
        "max_candidates_per_match": 10,
        "include_zero_delta_as_match": True,
    },

    # 2026-09-02訂正: ピーク探索は mz_tolerance_ppm ではなく max_delta_da
    # (絶対質量幅, Da単位) を主探索フィルタとする。狭いppm窓では実際の修飾
    # (14〜288 Da程度のシフト) を検出できないため。min_charge/max_charge/
    # polarity は fragment_mapping を流用する（仕様書 §9.2, §12.5参照）。
    "mass_comparison": {
        "enabled": True,
        "max_delta_da": 300,
        "max_matches_per_fragment": 50,
        # §24 (Phase 10): P1 complete-digestion mode (digestion.enzyme ==
        # "Nuclease_P1") matches peaks directly against a small, fixed set
        # of known nucleoside masses, so it uses a narrow ppm tolerance
        # (like fragment_mapping's, not the wide max_delta_da search) —
        # there's no unknown ΔDa to explore, just known masses to confirm.
        "nucleoside_mz_tolerance_ppm": 10,
    },

    # hypothesis_mass_check_spec.md §3.2: user-specified modification
    # hypotheses (in addition to a selected tRNA's conserved_modifications
    # from data/trna_library.yaml) whose theoretical mass is checked against
    # the raw MS1 peaks — a yes/no presence check only, never an automatic
    # identification of what an observed mass "is". Empty `targets` by
    # default; see rna_masshunter/hypothesis_check.py for the target forms
    # (label / components+linkage / base+add_elements).
    "hypothesis_check": {
        "enabled": True,
        "mass_tolerance_ppm": 15,
        "max_charge": 4,
        "targets": [],
    },

    "reporting": {
        "excel_output": True,
        "output_filename": "RNA_MassHunter_simple_report.xlsx",
        "max_excel_rows_per_sheet": 100000,
        "truncate_large_sheets": True,
    },

    # 08_Visualization scatter chart (Charge x Observed Neutral Mass, §17).
    # A dedicated section rather than overloading mass_comparison, since it
    # is purely about report presentation, not the matching logic itself.
    "visualization": {
        "enabled": True,
    },
}


def _merge_defaults(data: dict[str, Any], warnings: list[dict[str, Any]] | None) -> dict[str, Any]:
    merged = {}
    for section, defaults in DEFAULT_CONFIG.items():
        current = data.get(section, {})
        if not isinstance(current, dict):
            current = {}
            if warnings is not None:
                add_warning(warnings, "WARNING", "config", f"Config section '{section}' was not a mapping; defaults were used.")
        missing = sorted(set(defaults) - set(current))
        if missing and warnings is not None:
            add_warning(warnings, "WARNING", "config", f"Config section '{section}' missing keys filled by defaults.", missing)
        merged[section] = {**defaults, **current}
    merged["raw"] = data
    return merged


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def load_config(config_path: str | Path, warnings: list[dict[str, Any]] | None = None) -> RunConfig:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"config.yaml not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    merged = _merge_defaults(data, warnings)
    return RunConfig(**merged)


_VALID_CCA_MODES = {"EXCLUDES_CCA", "INCLUDES_COMPLETE_CCA", "UNKNOWN"}
_VALID_ELEMENTS = {"C", "H", "N", "O", "P", "S", "Se"}
_VALID_ENZYMES = {"RNase_T1", "RNase_A", "RNase_T2", "Nuclease_P1", "Benzonase", "U_specific_RNase"}


def validate_config(config: RunConfig, warnings: list[dict[str, Any]] | None = None) -> None:
    """Minimal structural validation. Extend as new config-driven checks are
    added (see 仕様書 §18 for the full list of error cases to cover)."""
    if config.instrument.get("polarity") not in ("negative", "positive"):
        raise ValueError("instrument.polarity must be 'negative' or 'positive'")

    if config.cca_processing.get("registered_sequence_cca_mode") not in _VALID_CCA_MODES:
        raise ValueError(f"cca_processing.registered_sequence_cca_mode must be one of {sorted(_VALID_CCA_MODES)}")

    enzyme = config.digestion.get("enzyme")
    if config.digestion.get("enabled") and enzyme not in _VALID_ENZYMES:
        raise ValueError(f"digestion.enzyme '{enzyme}' is not a known enzyme ({sorted(_VALID_ENZYMES)})")

    fc = config.formula_candidate
    if fc.get("enabled"):
        max_total_atoms = fc.get("max_total_atoms")
        if not isinstance(max_total_atoms, int) or max_total_atoms <= 0:
            raise ValueError("formula_candidate.max_total_atoms must be a positive integer")
        unknown_elements = set(fc.get("elements", [])) - _VALID_ELEMENTS
        if unknown_elements:
            raise ValueError(f"formula_candidate.elements contains unsupported elements: {sorted(unknown_elements)}")
        unknown_limit_elements = set(fc.get("element_limits", {})) - _VALID_ELEMENTS
        if unknown_limit_elements:
            raise ValueError(f"formula_candidate.element_limits contains unsupported elements: {sorted(unknown_limit_elements)}")
        tolerance = fc.get("mass_tolerance_da")
        if not isinstance(tolerance, (int, float)) or tolerance < 0:
            raise ValueError("formula_candidate.mass_tolerance_da must be a non-negative number")

    mc = config.mass_comparison
    if mc.get("enabled"):
        max_delta_da = mc.get("max_delta_da")
        if not isinstance(max_delta_da, (int, float)) or max_delta_da <= 0:
            raise ValueError("mass_comparison.max_delta_da must be a positive number")
        max_matches = mc.get("max_matches_per_fragment")
        if not isinstance(max_matches, int) or max_matches <= 0:
            raise ValueError("mass_comparison.max_matches_per_fragment must be a positive integer")

    hc = config.hypothesis_check
    if hc.get("enabled"):
        tolerance_ppm = hc.get("mass_tolerance_ppm")
        if isinstance(tolerance_ppm, bool) or not isinstance(tolerance_ppm, (int, float)) or tolerance_ppm <= 0:
            raise ValueError("hypothesis_check.mass_tolerance_ppm must be a positive number")
        max_charge = hc.get("max_charge")
        if isinstance(max_charge, bool) or not isinstance(max_charge, int) or max_charge <= 0:
            raise ValueError("hypothesis_check.max_charge must be a positive integer")
        if hc.get("targets") is not None and not isinstance(hc.get("targets"), list):
            raise ValueError("hypothesis_check.targets must be a list")


def resolve_paths(config: RunConfig, project_root: str | Path) -> RunConfig:
    root = Path(project_root)
    for section_name in ("project", "input"):
        section = getattr(config, section_name)
        for key, value in list(section.items()):
            if not value or not isinstance(value, str):
                continue
            if key.endswith("_dir") or key.endswith("_path"):
                path = Path(value)
                if not path.is_absolute():
                    section[key] = str(root / path)
    return config
