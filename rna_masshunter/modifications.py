from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.models import Modification
from rna_masshunter.warnings_manager import add_warning

MASS_SHIFT_KEYS = ("mass_shift_from_unmodified", "mass_shift", "delta_mass", "delta_mass_da", "mass_difference")


def _records_from_yaml(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict):
        for key in ("modifications", "entries", "data"):
            if isinstance(data.get(key), list):
                return [item for item in data[key] if isinstance(item, dict)]
        return [dict(value, id=value.get("id", name)) for name, value in data.items() if isinstance(value, dict)]
    return []


def _mass_shift(record: dict[str, Any]) -> float | None:
    for key in MASS_SHIFT_KEYS:
        if key in record:
            try:
                return float(record[key])
            except (TypeError, ValueError):
                return None
    return None


def load_modifications(path: str | Path, warnings: list[dict[str, Any]] | None = None) -> list[Modification]:
    yaml_path = Path(path)
    with yaml_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    modifications = []
    for record in _records_from_yaml(data):
        mass_shift = _mass_shift(record)
        mod_id = str(record.get("id") or record.get("symbol") or "")
        target_bases = record.get("target_bases", [])
        if isinstance(target_bases, str):
            target_bases = [target_bases]
        modifications.append(
            Modification(
                id=mod_id,
                symbol=record.get("symbol"),
                mass_shift_from_unmodified=float(mass_shift) if mass_shift is not None else float("nan"),
                category=str(record.get("category", "")),
                target_bases=list(target_bases) if isinstance(target_bases, list) else [],
                detectability=record.get("detectability"),
                curation=record.get("curation"),
                sources=record.get("sources"),
                source=record.get("source"),
                source_priority=record.get("source_priority") or (record.get("source", {}) or {}).get("source_priority"),
                curation_status=str(record.get("curation_status") or (record.get("source", {}) or {}).get("curation_status") or (record.get("curation", {}) or {}).get("status") or ""),
                candidate_policy=dict(record.get("candidate_policy") or {}),
                chemical_group=str(record.get("chemical_group") or ""),
                near_isobaric_group=str(record.get("near_isobaric_group") or ""),
                raw=record,
            )
        )
    if not modifications and warnings is not None:
        add_warning(warnings, "WARNING", "modifications", "No modifications were loaded.", str(yaml_path))
    return modifications


def validate_modifications(modifications: list[Modification], warnings: list[dict[str, Any]] | None = None) -> None:
    for item in modifications:
        if not item.id and not item.symbol and warnings is not None:
            add_warning(warnings, "ERROR", "modifications", "Modification lacks id and symbol.", item.raw)
        if item.mass_shift_from_unmodified != item.mass_shift_from_unmodified and warnings is not None:
            add_warning(warnings, "ERROR", "modifications", "Modification mass shift is not numeric.", item.id or item.symbol)
        if not item.category and warnings is not None:
            add_warning(warnings, "ERROR", "modifications", "Modification category is missing.", item.id or item.symbol)
        if not item.target_bases and warnings is not None:
            add_warning(warnings, "ERROR", "modifications", "Modification target_bases is missing.", item.id or item.symbol)


def find_modifications_by_mass_shift(modifications: list[Modification], mass_shift: float, tolerance_da: float = 0.01) -> list[Modification]:
    return [item for item in modifications if abs(item.mass_shift_from_unmodified - mass_shift) <= tolerance_da]
