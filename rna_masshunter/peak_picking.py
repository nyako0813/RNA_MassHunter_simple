from pathlib import Path
from typing import Any

import numpy as np

from rna_masshunter.models import Peak
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra


def extract_ms1_peaks(mzml_path: str | Path, reconstruction_config: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> list[Peak]:
    rt_min = reconstruction_config.get("rt_min")
    rt_max = reconstruction_config.get("rt_max")
    mz_min = float(reconstruction_config.get("mz_min", 0))
    mz_max = float(reconstruction_config.get("mz_max", float("inf")))
    intensity_threshold = float(reconstruction_config.get("intensity_threshold", 0))
    peaks: list[Peak] = []

    for spectrum in iter_spectra(mzml_path):
        if int(spectrum.get("ms level", 0)) != 1:
            continue
        rt = _rt_minutes(spectrum)
        if rt is not None and rt_min is not None and rt < float(rt_min):
            continue
        if rt is not None and rt_max is not None and rt > float(rt_max):
            continue
        mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
        intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
        if mz_array.size != intensity_array.size:
            continue
        mask = (mz_array >= mz_min) & (mz_array <= mz_max) & (intensity_array >= intensity_threshold)
        for mz_value, intensity in zip(mz_array[mask], intensity_array[mask], strict=False):
            peaks.append(Peak(mz=float(mz_value), intensity=float(intensity), rt=rt, scan_id=str(spectrum.get("id", ""))))
    return peaks
