from typing import Any

from rna_masshunter.masses import PROTON_MASS
from rna_masshunter.models import Fragment, FragmentMS1Match, Peak, RunConfig
from rna_masshunter.warnings_manager import add_warning


def map_ms1_fragments(*args, **kwargs):
    """Backward-compatible alias for MVP-3 MS1 fragment mapping."""
    return map_fragments_to_ms1_peaks(*args, **kwargs)


def theoretical_mz_from_mass(
    neutral_mass: float,
    charge: int,
    polarity: str,
    proton_mass: float = PROTON_MASS,
) -> float:
    z = abs(int(charge))
    mode = str(polarity).lower()
    if mode == "negative":
        return (float(neutral_mass) - z * proton_mass) / z
    if mode == "positive":
        return (float(neutral_mass) + z * proton_mass) / z
    raise ValueError(f"Unknown polarity '{polarity}'. Expected 'negative' or 'positive'.")


def ppm_error(observed: float, theoretical: float) -> float:
    if theoretical == 0:
        return 0.0
    return (float(observed) - float(theoretical)) / float(theoretical) * 1_000_000


def _confidence(error_ppm: float, tolerance_ppm: float, peak_tier: str | None) -> str:
    tier = str(peak_tier or "").lower()
    abs_error = abs(error_ppm)
    if abs_error <= tolerance_ppm / 3 and tier in {"major", "minor"}:
        return "High"
    if abs_error <= tolerance_ppm / 2:
        return "Medium"
    return "Low"


def _coerce_charge(value: Any, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def retention_sort_key(item: FragmentMS1Match) -> tuple[float, float]:
    """Return the unchanged formal truncation key: abs ppm asc, intensity desc."""
    return (abs(item.mass_error_ppm), -item.intensity)


def _physical_peak_id(peak: Peak, peak_index: int) -> str:
    scan_id = str(getattr(peak, "scan_id", None) or "no_scan")
    rt = getattr(peak, "rt", None)
    rt_text = "" if rt is None else f"{float(rt):.8f}"
    return f"PK_{peak_index:06d}_{scan_id}_{rt_text}_{float(peak.mz):.8f}"


def _eligible_peaks(peaks: list[Peak], mapping_config: dict[str, Any]) -> list[Peak]:
    if not mapping_config.get("use_peak_tiers", True):
        return list(peaks)

    include_trace = bool(mapping_config.get("include_trace_peaks", True))
    eligible = []
    for peak in peaks:
        tier = getattr(peak, "tier", None)
        tier_key = str(tier or "").lower()
        if tier_key in {"major", "minor"}:
            eligible.append(peak)
        elif include_trace and tier_key == "trace":
            eligible.append(peak)
        elif tier is None:
            eligible.append(peak)
    return eligible


def map_fragments_to_ms1_peaks(
    fragments: list[Fragment],
    peaks: list[Peak],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None = None,
    audit_context: dict[str, Any] | None = None,
) -> list[FragmentMS1Match]:
    mapping_config = config.fragment_mapping or {}
    if not mapping_config.get("enabled", True):
        return []

    missing_inputs = False
    if not fragments:
        missing_inputs = True
        if warnings is not None:
            add_warning(warnings, "WARNING", "ms1_mapping", "fragment_mapping.enabled is true but theoretical_fragments is empty.")
    if not peaks:
        missing_inputs = True
        if warnings is not None:
            add_warning(warnings, "WARNING", "ms1_mapping", "fragment_mapping.enabled is true but MS1 peaks are empty.")
    if missing_inputs:
        if audit_context is not None:
            configured_max = _coerce_charge(mapping_config.get("max_matches_per_fragment"), 20)
            if configured_max < 1:
                configured_max = 20
            audit_context["fragments"] = [
                {"fragment": fragment, "ranked_matches": [], "configured_max_matches": configured_max}
                for fragment in fragments
            ]
            audit_context["configured_max_matches"] = configured_max
            audit_context["eligible_peak_count"] = 0
            audit_context["sort_key"] = "(abs(mass_error_ppm) asc, intensity desc); stable ties use charge asc then eligible peak input order"
        return []

    polarity = str(mapping_config.get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str(config.instrument.get("polarity", "negative") or "negative").lower()
    if polarity not in {"negative", "positive"}:
        if warnings is not None:
            add_warning(warnings, "ERROR", "ms1_mapping", "fragment_mapping polarity is unknown.", polarity)
        return []

    min_charge = _coerce_charge(mapping_config.get("min_charge"), 1)
    max_charge = _coerce_charge(mapping_config.get("max_charge"), 8)
    if min_charge < 1 or max_charge < min_charge:
        if warnings is not None:
            add_warning(
                warnings,
                "ERROR",
                "ms1_mapping",
                "fragment_mapping min_charge / max_charge is invalid.",
                {"min_charge": min_charge, "max_charge": max_charge},
            )
        return []

    tolerance_ppm = float(mapping_config.get("mz_tolerance_ppm", 10) or 10)
    max_matches = _coerce_charge(mapping_config.get("max_matches_per_fragment"), 20)
    if max_matches < 1:
        max_matches = 20

    eligible_peaks = _eligible_peaks(peaks, mapping_config)
    if not eligible_peaks:
        if warnings is not None:
            add_warning(warnings, "WARNING", "ms1_mapping", "No MS1 peaks remained after fragment_mapping peak tier filtering.")
        if audit_context is not None:
            audit_context["fragments"] = [
                {"fragment": fragment, "ranked_matches": [], "configured_max_matches": max_matches}
                for fragment in fragments
            ]
            audit_context["configured_max_matches"] = max_matches
            audit_context["eligible_peak_count"] = 0
            audit_context["sort_key"] = "(abs(mass_error_ppm) asc, intensity desc); stable ties use charge asc then eligible peak input order"
        return []

    matches: list[FragmentMS1Match] = []
    for fragment in fragments:
        fragment_matches: list[FragmentMS1Match] = []
        for charge in range(min_charge, max_charge + 1):
            theoretical_mz = theoretical_mz_from_mass(fragment.unmodified_mass, charge, polarity)
            for peak_index, peak in enumerate(eligible_peaks, start=1):
                observed_mz = float(getattr(peak, "mz"))
                error_ppm = ppm_error(observed_mz, theoretical_mz)
                if abs(error_ppm) > tolerance_ppm:
                    continue
                error_da = observed_mz - theoretical_mz
                match = FragmentMS1Match(
                        match_id=f"{fragment.fragment_id}_z{charge}_{len(fragment_matches) + 1}",
                        fragment_id=fragment.fragment_id,
                        target_id=fragment.target_id,
                        sequence=fragment.sequence,
                        start=fragment.start,
                        end=fragment.end,
                        standard_start=fragment.standard_start,
                        standard_end=fragment.standard_end,
                        enzyme=fragment.enzyme,
                        missed_cleavages=fragment.missed_cleavages,
                        terminal_form=fragment.terminal_form,
                        fragment_mass=fragment.unmodified_mass,
                        charge=charge,
                        theoretical_mz=theoretical_mz,
                        observed_mz=observed_mz,
                        mass_error_da=error_da,
                        mass_error_ppm=error_ppm,
                        intensity=float(getattr(peak, "intensity", 0.0) or 0.0),
                        rt=getattr(peak, "rt", None),
                        scan_id=getattr(peak, "scan_id", None),
                        peak_tier=getattr(peak, "tier", None),
                        confidence=_confidence(error_ppm, tolerance_ppm, getattr(peak, "tier", None)),
                    )
                setattr(match, "_audit_peak_index", peak_index)
                setattr(match, "_audit_physical_peak_id", _physical_peak_id(peak, peak_index))
                setattr(match, "_audit_generation_order", len(fragment_matches) + 1)
                fragment_matches.append(match)

        ranked_for_audit = sorted(fragment_matches, key=retention_sort_key)
        if audit_context is not None:
            audit_context.setdefault("fragments", []).append(
                {
                    "fragment": fragment,
                    "ranked_matches": list(ranked_for_audit),
                    "configured_max_matches": max_matches,
                }
            )
            audit_context["sort_key"] = "(abs(mass_error_ppm) asc, intensity desc); stable ties use charge asc then eligible peak input order"
            audit_context["configured_max_matches"] = max_matches
            audit_context["eligible_peak_count"] = len(eligible_peaks)

        if len(fragment_matches) > max_matches:
            fragment_matches.sort(key=retention_sort_key)
            if warnings is not None:
                add_warning(
                    warnings,
                    "WARNING",
                    "ms1_mapping",
                    "Fragment MS1 matches were truncated by max_matches_per_fragment.",
                    {"fragment_id": fragment.fragment_id, "before": len(fragment_matches), "after": max_matches},
                )
            fragment_matches = fragment_matches[:max_matches]
        matches.extend(fragment_matches)

    for index, match in enumerate(matches, start=1):
        match.match_id = f"MS1F_{index:06d}"
    return matches
