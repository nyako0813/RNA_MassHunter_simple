from dataclasses import replace
from typing import Any

from rna_masshunter.enzymes import find_cleavage_sites, get_enzyme_rule, normalize_enzyme_name
from rna_masshunter.masses import calculate_unmodified_rna_mass
from rna_masshunter.models import Fragment, RunConfig
from rna_masshunter.warnings_manager import add_warning


NO_DIGESTION_ENZYMES = {"", "none", "no_digestion", "intact", "full_length"}


def digest_rna(*args, **kwargs):
    """Backward-compatible alias for MVP-2 digestion."""
    return digest_sequence(*args, **kwargs)


def _fragment_warning(fragment_warnings: list[str], message: str) -> None:
    if message not in fragment_warnings:
        fragment_warnings.append(message)

def _resolve_three_prime_tail_sequences(full_sequence: str) -> tuple[str, ...]:
    """Determine which 3'-CCA maturation states are still worth hypothesizing.

    The CCA-adding enzyme appends C, C, then A one residue at a time onto a
    genomic/precursor 3' end that does not itself encode CCA. If the target
    sequence supplied by the user already carries part (or all) of that tail,
    re-appending the full "C"/"CC"/"CCA" set on top would double count
    residues that are already present (e.g. an input ending in "...CCA"
    would otherwise become "...CCACCA").

    This inspects the *existing* 3' end of the full sequence and only
    returns the remaining maturation states that are still biologically
    possible from there.
    """
    if full_sequence.endswith("CCA"):
        # Already fully matured; no further tail growth to hypothesize.
        return ("",)
    if full_sequence.endswith("CC"):
        # Only the final A addition remains possible.
        return ("", "A")
    if full_sequence.endswith("C"):
        # Either stalled after the first C, or completed with "CA".
        return ("", "CA")
    # No CCA-related residues encoded yet; all maturation states are possible.
    return ("", "C", "CC", "CCA")


def _generate_three_prime_tail_candidates(
    fragment: Fragment,
    full_sequence: str,
    base_masses: dict,
    warnings: list[dict[str, Any]] | None = None,
) -> list[Fragment]:
    """Generate 3'-terminal C/CC/CCA sequence candidates.

    Only fragments reaching the original RNA 3' end are eligible.
    The candidates are ordinary Fragment objects so the existing MS1
    mapping can match their theoretical masses directly.
    """
    if fragment.end != len(full_sequence):
        return [fragment]

    tail_sequences = _resolve_three_prime_tail_sequences(full_sequence)
    candidates: list[Fragment] = []

    for tail in tail_sequences:
        sequence = fragment.sequence + tail
        end = fragment.end + len(tail)

        mass = calculate_unmodified_rna_mass(
            sequence,
            base_masses,
            warnings=warnings,
            terminal_form=fragment.terminal_form,
        )

        fragment_warnings = list(fragment.warnings)
        if mass is None:
            _fragment_warning(
                fragment_warnings,
                "unmodified mass could not be calculated",
            )
            mass = 0.0

        candidates.append(
            replace(
                fragment,
                fragment_id=(
                    f"{fragment.target_id}_{fragment.enzyme}_"
                    f"{fragment.start}_{end}_"
                    f"mc{fragment.missed_cleavages}_{fragment.terminal_form}"
                ),
                sequence=sequence,
                end=end,
                standard_end=(
                    fragment.standard_end
                    if not tail
                    else None
                ),
                unmodified_mass=float(mass),
                warnings=fragment_warnings,
            )
        )

    return candidates

def digest_sequence(
    target_id: str,
    sequence: str,
    position_map: dict[int, int | None],
    config: RunConfig,
    base_masses: dict,
    warnings: list[dict[str, Any]] | None = None,
) -> list[Fragment]:
    sequence = (sequence or "").upper().replace("T", "U")
    if not sequence:
        return []

    digestion_config = config.digestion or {}
    if not digestion_config.get("enabled", True):
        return []

    enzyme = normalize_enzyme_name(digestion_config.get("enzyme", "RNase_T1"))
    if enzyme.lower() in NO_DIGESTION_ENZYMES:
        if warnings is not None:
            add_warning(
                warnings,
                "INFO",
                "digestion",
                "No digestion enzyme was selected; theoretical digestion fragments were not generated.",
                {"target_id": target_id, "enzyme": enzyme},
            )
        return []

    digestion_mode = str(digestion_config.get("digestion_mode") or ("complete" if enzyme == "Nuclease_P1" else "specific")).lower()
    allow_nonspecific = bool(digestion_config.get("allow_nonspecific_cleavage", False)) or digestion_mode in {"complete", "nonspecific"}
    try:
        rule = get_enzyme_rule(enzyme)
    except ValueError as exc:
        if warnings is not None:
            add_warning(warnings, "ERROR", "digestion", str(exc), {"target_id": target_id})
        return []

    if not rule.get("specific", True) and not allow_nonspecific:
        if warnings is not None:
            add_warning(warnings, "WARNING", "digestion", "Nonspecific cleavage enzyme selected but allow_nonspecific_cleavage is false; no theoretical fragments were generated.", {"enzyme": enzyme})
        return []
    if not rule.get("specific", True) and allow_nonspecific and len(sequence) > 50 and warnings is not None:
        add_warning(warnings, "WARNING", "digestion", "Nonspecific cleavage is enabled and may generate many candidates.", {"enzyme": enzyme, "sequence_length": len(sequence)})

    ap_config = config.alkaline_phosphatase or {}
    if (
        warnings is not None
        and ap_config.get("enabled", False)
        and not ap_config.get("assume_complete", False)
        and digestion_config.get("include_terminal_forms", True)
    ):
        add_warning(
            warnings,
            "WARNING",
            "digestion",
            "alkaline_phosphatase.enabled is true but assume_complete is false; residual phosphate and cyclic phosphate forms may be included.",
            {"target_id": target_id},
        )

    try:
        cleavage_sites = find_cleavage_sites(sequence, enzyme, warnings=warnings)
    except ValueError:
        return []

    boundaries = [0] + sorted(set(cleavage_sites))
    if boundaries[-1] != len(sequence):
        boundaries.append(len(sequence))

    max_missed = max(0, int(digestion_config.get("missed_cleavages", 0) or 0))
    min_length = max(1, int(digestion_config.get("min_length", 1) or 1))
    max_length_value = digestion_config.get("max_length")
    max_length = int(max_length_value) if max_length_value not in (None, "") else None
    fragments: list[Fragment] = []

    fragment_ranges: list[tuple[int, int, int]] = []
    if digestion_mode in {"complete", "nonspecific"}:
        effective_max_length = max_length or len(sequence)
        for start in range(1, len(sequence) + 1):
            for end in range(start, min(len(sequence), start + effective_max_length - 1) + 1):
                fragment_ranges.append((start, end, 0))
    else:
        for start_boundary_index in range(len(boundaries) - 1):
            for missed in range(max_missed + 1):
                end_boundary_index = start_boundary_index + missed + 1
                if end_boundary_index >= len(boundaries):
                    continue
                fragment_ranges.append((boundaries[start_boundary_index] + 1, boundaries[end_boundary_index], missed))

    for start, end, missed in fragment_ranges:
        fragment_sequence = sequence[start - 1:end]
        if len(fragment_sequence) < min_length:
            continue
        if max_length is not None and len(fragment_sequence) > max_length:
            continue

        fragment_warnings: list[str] = []
        mass = calculate_unmodified_rna_mass(fragment_sequence, base_masses, warnings=warnings, terminal_form="default")
        if mass is None:
            _fragment_warning(fragment_warnings, "unmodified mass could not be calculated")
            mass = 0.0

        fragment = Fragment(
            fragment_id=f"{target_id}_{enzyme}_{start}_{end}_mc{missed}_default",
            target_id=target_id,
            sequence=fragment_sequence,
            start=start,
            end=end,
            standard_start=position_map.get(start),
            standard_end=position_map.get(end),
            enzyme=enzyme,
            missed_cleavages=missed,
            terminal_form="default",
            unmodified_mass=float(mass),
            warnings=fragment_warnings,
        )
        terminal_fragments = generate_terminal_forms(
            fragment,
            config,
            base_masses,
            warnings=warnings,
        )

        for terminal_fragment in terminal_fragments:
            fragments.extend(
                _generate_three_prime_tail_candidates(
                    terminal_fragment,
                    sequence,
                    base_masses,
                    warnings=warnings,
                )
            )
    return fragments


def generate_terminal_forms(
    fragment: Fragment,
    config: RunConfig,
    base_masses: dict,
    warnings: list[dict[str, Any]] | None = None,
) -> list[Fragment]:
    digestion_config = config.digestion or {}
    if not digestion_config.get("include_terminal_forms", True):
        return [fragment]

    ap_config = config.alkaline_phosphatase or {}
    if not ap_config.get("enabled", False):
        return [fragment]

    forms = ["dephosphorylated"]
    if not ap_config.get("assume_complete", False):
        if ap_config.get("allow_residual_phosphate", True):
            forms.append("residual_phosphate")
        if ap_config.get("allow_cyclic_phosphate", True):
            forms.append("cyclic_phosphate")

    terminal_fragments: list[Fragment] = []
    for terminal_form in forms:
        mass = calculate_unmodified_rna_mass(fragment.sequence, base_masses, warnings=warnings, terminal_form=terminal_form)
        fragment_warnings = list(fragment.warnings)
        if mass is None:
            _fragment_warning(fragment_warnings, "unmodified mass could not be calculated")
            mass = 0.0
        terminal_fragments.append(
            replace(
                fragment,
                fragment_id=f"{fragment.target_id}_{fragment.enzyme}_{fragment.start}_{fragment.end}_mc{fragment.missed_cleavages}_{terminal_form}",
                terminal_form=terminal_form,
                unmodified_mass=float(mass),
                warnings=fragment_warnings,
            )
        )
    return terminal_fragments