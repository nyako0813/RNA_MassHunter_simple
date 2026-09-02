"""MS2 spectrum extraction and unmodified theoretical ion generation.

Vendored (with edits) from nyako0813/RNA_MassHunter `rna_masshunter/ms2_annotation.py`.

WHY THIS IS A SEPARATE FILE AND NOT A DIRECT IMPORT OF ms2_annotation.py:
the original `ms2_annotation.py` imports, at module load time,
`base_loss_masses`, `base_loss_ions`, `modified_precursor`,
`modified_fragment_ions`, `ms2_unmatched_audit`, and
`ms2_zero_intensity_audit` — a cluster of modules for full MS2
localization/evidence scoring that this simplified project deliberately
does not use (see docs/design/RNA_MassHunter_再設計_実装仕様書.md §3.9,
§14A, and Claude Code rule #14). `from rna_masshunter.ms2_annotation import
extract_ms2_spectra` would therefore drag all of that in transitively.

Only two functions are needed here: `extract_ms2_spectra` (parses MS2
scans out of an mzML file) and `generate_theoretical_ms2_ions` (unmodified
d/w/a/z-series theoretical fragment ions). Both are copied below with one
substantive edit: the original `extract_ms2_spectra` unconditionally calls
`ms2_zero_intensity_audit.capture_source_spectrum` /
`record_parsed_spectrum` / `record_parser_error` for a diagnostics feature
this project does not use — those calls (and the now-dead
`source_record`/`parsed_mz_array`/`parsed_intensity_array`/
`annotation_indices` bookkeeping that only existed to feed them) have been
removed. No other logic was changed. If you re-vendor this file later from
a newer version of the original repo, re-apply the same removal.
"""
from typing import Any

import numpy as np

from rna_masshunter.masses import (
    C_TERMINAL_ION_TYPES,
    N_TERMINAL_ION_TYPES,
    calculate_unmodified_rna_mass,
    fragment_ion_series_offsets,
)
from rna_masshunter.models import Fragment, MS2SpectrumInfo, TheoreticalMS2Ion
from rna_masshunter.ms1_mapping import theoretical_mz_from_mass
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra
from rna_masshunter.warnings_manager import add_warning


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _optional_positive_int(value: Any, default: int | None = None) -> int | None:
    if value in (None, ""):
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _scan_window(spectrum: dict[str, Any]) -> tuple[float | None, float | None]:
    """Read acquisition scan bounds only when explicit mzML metadata is present."""
    scans = (spectrum.get("scanList") or {}).get("scan") or []
    for scan in scans:
        windows = ((scan or {}).get("scanWindowList") or {}).get("scanWindow") or []
        for window in windows:
            lower = _safe_float((window or {}).get("scan window lower limit"))
            upper = _safe_float((window or {}).get("scan window upper limit"))
            if lower is not None and upper is not None:
                return lower, upper
    return None, None


def _precursor_info(spectrum: dict[str, Any]) -> dict[str, Any]:
    info = {"mz": None, "charge": None, "intensity": None}
    precursor_list = spectrum.get("precursorList", {}).get("precursor", [])
    if not precursor_list:
        return info
    selected = precursor_list[0].get("selectedIonList", {}).get("selectedIon", [])
    if not selected:
        return info
    ion = selected[0]
    info["mz"] = _safe_float(ion.get("selected ion m/z") or ion.get("isolation window target m/z"))
    charge = ion.get("charge state")
    try:
        info["charge"] = int(charge) if charge is not None else None
    except (TypeError, ValueError):
        info["charge"] = None
    info["intensity"] = _safe_float(ion.get("peak intensity") or ion.get("ion intensity"))
    return info


def _ion_charges(ms2_config: dict[str, Any]) -> list[int]:
    raw = ms2_config.get("charge_states") or [1]
    values = raw if isinstance(raw, list) else [raw]
    charges = []
    for value in values:
        try:
            charge = abs(int(value))
        except (TypeError, ValueError):
            continue
        if charge > 0 and charge not in charges:
            charges.append(charge)
    return charges or [1]


def _ion_series(ms2_config: dict[str, Any]) -> list[str]:
    raw = ms2_config.get("ion_series") or ["d", "w", "a", "z"]
    values = raw if isinstance(raw, list) else [raw]
    series = []
    for value in values:
        ion_type = str(value).strip().lower()
        if ion_type in (N_TERMINAL_ION_TYPES | C_TERMINAL_ION_TYPES) and ion_type not in series:
            series.append(ion_type)
    return series or ["d", "w", "a", "z"]


def extract_ms2_spectra(
    mzml_path: str,
    ms2_config: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[MS2SpectrumInfo]:
    spectra: list[MS2SpectrumInfo] = []
    min_intensity = float(ms2_config.get("min_peak_intensity", 10) or 0)
    min_relative_percent = float(ms2_config.get("min_relative_intensity_percent", 1.0) or 0)
    max_peaks = _optional_positive_int(ms2_config.get("max_peaks_per_spectrum"), 500)

    try:
        for scan_index, spectrum in enumerate(iter_spectra(mzml_path), start=1):
            try:
                ms_level = int(spectrum.get("ms level", 0) or 0)
            except (TypeError, ValueError):
                ms_level = 0
            if ms_level != 2:
                continue

            mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
            intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
            if mz_array.size != intensity_array.size:
                if warnings is not None:
                    add_warning(warnings, "WARNING", "ms2_extraction", "MS2 spectrum m/z and intensity arrays had different lengths.", spectrum.get("id"))
                continue

            raw_base_peak_index = int(np.argmax(intensity_array)) if intensity_array.size else None
            base_peak_mz = float(mz_array[raw_base_peak_index]) if raw_base_peak_index is not None else None
            base_peak_intensity = float(intensity_array[raw_base_peak_index]) if raw_base_peak_index is not None else None
            total_ion_current = float(np.sum(intensity_array)) if intensity_array.size else 0.0
            raw_peaks = [(float(mz), float(intensity)) for mz, intensity in zip(mz_array, intensity_array, strict=False)]
            scan_mz_min, scan_mz_max = _scan_window(spectrum)
            effective_intensity_threshold = max(
                min_intensity,
                float(base_peak_intensity or 0.0) * min_relative_percent / 100.0,
            )

            if intensity_array.size:
                mask = intensity_array >= min_intensity
                if base_peak_intensity and min_relative_percent > 0:
                    mask = mask & (intensity_array >= base_peak_intensity * min_relative_percent / 100.0)
                mz_array = mz_array[mask]
                intensity_array = intensity_array[mask]

            if max_peaks and mz_array.size > max_peaks:
                order = np.argsort(intensity_array)[::-1][:max_peaks]
                order = order[np.argsort(mz_array[order])]
                mz_array = mz_array[order]
                intensity_array = intensity_array[order]

            precursor = _precursor_info(spectrum)
            spectrum_info = MS2SpectrumInfo(
                spectrum_id=str(spectrum.get("id") or f"scan_{scan_index}"),
                scan_index=scan_index,
                rt=_safe_float(_rt_minutes(spectrum)),
                precursor_mz=precursor.get("mz"),
                precursor_charge=precursor.get("charge"),
                precursor_intensity=precursor.get("intensity"),
                num_peaks=int(mz_array.size),
                base_peak_mz=base_peak_mz,
                base_peak_intensity=base_peak_intensity,
                total_ion_current=total_ion_current,
                peaks=[(float(mz), float(intensity)) for mz, intensity in zip(mz_array, intensity_array, strict=False)],
                raw_peaks=raw_peaks, scan_mz_min=scan_mz_min, scan_mz_max=scan_mz_max,
                effective_intensity_threshold=effective_intensity_threshold,
                threshold_information_available=True,
            )
            spectra.append(spectrum_info)
    except Exception as exc:
        if warnings is not None:
            add_warning(warnings, "WARNING", "ms2_extraction", "MS2 spectrum extraction failed; annotation skipped.", str(exc))
    return spectra


def generate_theoretical_ms2_ions(
    theoretical_fragments: list[Fragment],
    config: Any,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[TheoreticalMS2Ion]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2_config.get("use_theoretical_fragments"), True):
        return []
    min_length = max(1, int(ms2_config.get("min_ion_length", 1) or 1))
    max_length = _optional_positive_int(ms2_config.get("max_ion_length"), None)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    charges = _ion_charges(ms2_config)
    ion_series = _ion_series(ms2_config)
    offsets = fragment_ion_series_offsets(base_masses, warnings=warnings)
    n_terminal_series = [ion_type for ion_type in ion_series if ion_type in N_TERMINAL_ION_TYPES]
    c_terminal_series = [ion_type for ion_type in ion_series if ion_type in C_TERMINAL_ION_TYPES]

    ions: list[TheoreticalMS2Ion] = []
    for fragment in theoretical_fragments or []:
        sequence = (fragment.sequence or "").upper().replace("T", "U")
        if len(sequence) < 2:
            continue
        for cut in range(1, len(sequence)):
            termini = [
                (n_terminal_series, sequence[:cut], 1, cut),
                (c_terminal_series, sequence[cut:], cut + 1, len(sequence)),
            ]
            for series_types, ion_sequence, ion_start, ion_end in termini:
                if not series_types:
                    continue
                ion_length = len(ion_sequence)
                if ion_length < min_length:
                    continue
                if max_length is not None and ion_length > max_length:
                    continue
                base_mass = calculate_unmodified_rna_mass(ion_sequence, base_masses, warnings=warnings, terminal_form="default")
                if base_mass is None:
                    continue
                for ion_type in series_types:
                    mass = float(base_mass) + offsets.get(ion_type, 0.0)
                    for charge in charges:
                        ions.append(TheoreticalMS2Ion(
                            ion_id=f"MS2ION_{len(ions) + 1:06d}",
                            parent_fragment_id=fragment.fragment_id,
                            parent_sequence=sequence,
                            ion_type=ion_type,
                            ion_sequence=ion_sequence,
                            ion_start=ion_start,
                            ion_end=ion_end,
                            charge=charge,
                            theoretical_mass=mass,
                            theoretical_mz=theoretical_mz_from_mass(mass, charge, polarity),
                            modification_id="",
                            modification_name="",
                            comment="unmodified d/w/a/z series ion",
                        ))
    return ions
