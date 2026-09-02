from pathlib import Path
from typing import Iterator

from pyteomics import mzml


def iter_spectra(mzml_path: str | Path) -> Iterator[dict]:
    with mzml.MzML(str(mzml_path)) as reader:
        for spectrum in reader:
            yield spectrum


def extract_ms1_spectra(mzml_path: str | Path) -> list[dict]:
    return [spectrum for spectrum in iter_spectra(mzml_path) if int(spectrum.get("ms level", 0)) == 1]


def extract_ms2_spectra(mzml_path: str | Path) -> list[dict]:
    return [spectrum for spectrum in iter_spectra(mzml_path) if int(spectrum.get("ms level", 0)) == 2]
