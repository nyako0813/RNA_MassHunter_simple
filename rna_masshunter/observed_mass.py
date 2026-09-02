"""Peak ID assignment and Observed_Mass/Mass_Intensity sheet row data.

New module (仕様書 §6, §8.3). Deliberately kept as a thin helper with no
real logic beyond ID bookkeeping and dataclass assembly — per 仕様書 §6
("既存機能で十分なら不要なモジュールは増やさない"), if it never grows
beyond this it can be folded into mass_comparison.py later.

`assign_peak_ids` is the single source of truth for Peak ID <-> Peak
mapping: mass_comparison.py and (later) excel_report.py both call it on the
same `peaks` list so that 04_Observed_Mass / 05_Mass_Intensity /
06_Mass_Comparison share identical Peak IDs (仕様書 §14, Peak IDによる
追跡性).
"""
from __future__ import annotations

from dataclasses import dataclass

from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import Peak


@dataclass
class ObservedMassRow:
    peak_id: str
    mz: float
    charge: int | None
    observed_mass: float | None
    intensity: float
    rt: float | None
    scan_id: str | None


def assign_peak_ids(peaks: list[Peak]) -> dict[int, str]:
    """Stable 'PK0001'-style IDs keyed by `peaks`' list index."""
    return {index: f"PK{index + 1:04d}" for index in range(len(peaks))}


def build_observed_mass_rows(
    peaks: list[Peak],
    charge_hint: dict[int, list[int]] | None = None,
    polarity: str = "negative",
) -> list[ObservedMassRow]:
    """One row per peak, or (仕様書 §15/§16.2, "1ピークにつき複数charge候補が
    あり得るため、charge単位で行を分ける" の推奨に従い) one row per
    (peak, charge) pair when `charge_hint` supplies more than one charge for
    a given peak index — e.g. the set of charges under which that peak
    showed up in Mass Comparison. Without a hint, charge/observed_mass stay
    None (Mass Comparison is what actually determines charge per peak)."""
    ids = assign_peak_ids(peaks)
    rows: list[ObservedMassRow] = []
    for index, peak in enumerate(peaks):
        charges = (charge_hint or {}).get(index) or [None]
        for charge in charges:
            observed_mass = neutral_mass_from_mz(peak.mz, charge, polarity) if charge else None
            rows.append(ObservedMassRow(
                peak_id=ids[index],
                mz=peak.mz,
                charge=charge,
                observed_mass=observed_mass,
                intensity=peak.intensity,
                rt=peak.rt,
                scan_id=peak.scan_id,
            ))
    return rows
