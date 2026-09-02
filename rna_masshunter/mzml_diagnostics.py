from pathlib import Path
from typing import Any

import numpy as np

from rna_masshunter.mzml_reader import iter_spectra
from rna_masshunter.warnings_manager import add_warning


def _rt_minutes(spectrum: dict[str, Any]) -> float | None:
    scan_list = spectrum.get("scanList", {}).get("scan", [])
    if not scan_list:
        return None
    scan = scan_list[0]
    rt = scan.get("scan start time")
    if rt is None:
        return None
    unit = str(scan.get("unitName", "")).lower()
    return float(rt) / 60.0 if "second" in unit else float(rt)


def _has_precursor_key(spectrum: dict[str, Any], key: str) -> bool:
    precursor_list = spectrum.get("precursorList", {}).get("precursor", [])
    for precursor in precursor_list:
        selected = precursor.get("selectedIonList", {}).get("selectedIon", [])
        if selected and key in selected[0]:
            return True
    return False


def run_mzml_diagnostics(mzml_path: str | Path, logger, warnings: list[dict[str, Any]]) -> dict[str, Any]:
    ms1_count = 0
    ms2_count = 0
    rts = []
    mz_ranges = []
    intensity_ranges = []
    scan_ids = []
    has_precursor_mz = False
    has_precursor_charge = False
    array_lengths = []

    for spectrum in iter_spectra(mzml_path):
        ms_level = int(spectrum.get("ms level", 0))
        if ms_level == 1:
            ms1_count += 1
        elif ms_level == 2:
            ms2_count += 1
            has_precursor_mz = has_precursor_mz or _has_precursor_key(spectrum, "selected ion m/z")
            has_precursor_charge = has_precursor_charge or _has_precursor_key(spectrum, "charge state")
        rt = _rt_minutes(spectrum)
        if rt is not None:
            rts.append(rt)
        mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
        intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
        if mz_array.size:
            mz_ranges.extend([float(np.nanmin(mz_array)), float(np.nanmax(mz_array))])
            array_lengths.append(int(mz_array.size))
        if intensity_array.size:
            intensity_ranges.extend([float(np.nanmin(intensity_array)), float(np.nanmax(intensity_array))])
        if spectrum.get("id"):
            scan_ids.append(str(spectrum.get("id")))

    centroid_guess = "centroid-like" if array_lengths and float(np.median(array_lengths)) < 5000 else "profile-like or dense centroid"
    diag_warnings = []
    if ms1_count == 0:
        diag_warnings.append("No MS1 spectra found.")
        add_warning(warnings, "ERROR", "mzml_diagnostics", "No MS1 spectra found.")

    diagnostics = {
        "Number of MS1 spectra": ms1_count,
        "Number of MS2 spectra": ms2_count,
        "RT range": f"{min(rts):.4f} - {max(rts):.4f} min" if rts else "",
        "m/z range": f"{min(mz_ranges):.4f} - {max(mz_ranges):.4f}" if mz_ranges else "",
        "intensity range": f"{min(intensity_ranges):.4f} - {max(intensity_ranges):.4f}" if intensity_ranges else "",
        "Has precursor m/z": has_precursor_mz,
        "Has precursor charge": has_precursor_charge,
        "Centroid/profile guess": centroid_guess,
        "Scan ID pattern": scan_ids[0] if scan_ids else "",
        "Warnings": "; ".join(diag_warnings),
    }
    logger.info("mzML diagnostics: %s MS1, %s MS2", ms1_count, ms2_count)
    return diagnostics
