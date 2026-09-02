from typing import Any

from rna_masshunter.warnings_manager import add_warning


ENZYME_RULES: dict[str, dict[str, Any]] = {
    "RNase_T1": {"name": "RNase_T1", "cleaves_after": {"G"}, "specific": True},
    "RNase_A": {"name": "RNase_A", "cleaves_after": {"C", "U"}, "specific": True},
    "RNase_T2": {"name": "RNase_T2", "cleaves_after": {"A", "C", "G", "U"}, "specific": False},
    "Nuclease_P1": {"name": "Nuclease_P1", "cleaves_after": {"A", "C", "G", "U"}, "specific": False},
    "Benzonase": {"name": "Benzonase", "cleaves_after": {"A", "C", "G", "U"}, "specific": False},
    "U_specific_RNase": {"name": "U_specific_RNase", "cleaves_after": {"U"}, "specific": True},
}

ENZYME_ALIASES = {
    "p1": "Nuclease_P1",
    "nuclease_p1": "Nuclease_P1",
    "nuclease p1": "Nuclease_P1",
    "Nuclease_P1": "Nuclease_P1",
}


def load_enzyme_definitions(*args, **kwargs):
    """Return built-in MVP-2 enzyme definitions for backward compatibility."""
    return list(ENZYME_RULES.values())


def normalize_enzyme_name(enzyme_name: str) -> str:
    text = str(enzyme_name or "").strip()
    return ENZYME_ALIASES.get(text, ENZYME_ALIASES.get(text.lower(), text))


def get_enzyme_rule(enzyme_name: str) -> dict:
    enzyme_name = normalize_enzyme_name(enzyme_name)
    try:
        return ENZYME_RULES[str(enzyme_name)]
    except KeyError as exc:
        known = ", ".join(sorted(ENZYME_RULES))
        raise ValueError(f"Unknown enzyme '{enzyme_name}'. Known enzymes: {known}") from exc


def find_cleavage_sites(
    sequence: str,
    enzyme_name: str,
    warnings: list[dict[str, Any]] | None = None,
) -> list[int]:
    enzyme_name = normalize_enzyme_name(enzyme_name)
    try:
        rule = get_enzyme_rule(enzyme_name)
    except ValueError as exc:
        if warnings is not None:
            add_warning(warnings, "ERROR", "enzymes", str(exc), {"enzyme": enzyme_name})
        raise

    cleaves_after = set(rule.get("cleaves_after", set()))
    sites = [index for index, base in enumerate(sequence.upper(), start=1) if base in cleaves_after]
    if len(sequence) not in sites:
        sites.append(len(sequence))
    return sorted(set(site for site in sites if site > 0))
