from pathlib import Path

from rna_masshunter.models import RunConfig
from rna_masshunter.trna_library import apply_trna_type, load_trna_library

REPO_ROOT = Path(__file__).resolve().parent.parent
LIBRARY_PATH = REPO_ROOT / "data" / "trna_library.yaml"

_KNOWN_SEQUENCE = (
    "AGUCCUGUAGGGUAGUGGUCAAUCCUUCGGGCCUUUGGAGCCCGGGACAGCGGUUCGAAUCCGCUCAGGACUA"
)


def _make_config(trna_type: str = "", sequence: str = "", anticodon: str = "") -> RunConfig:
    return RunConfig(sequence={"name": "target_tRNA", "trna_type": trna_type, "sequence": sequence, "anticodon": anticodon})


def test_trna_library_has_58_unique_ids():
    library = load_trna_library(LIBRARY_PATH)
    assert len(library) == 58


def test_load_trna_library_returns_id_to_entry_mapping():
    library = load_trna_library(LIBRARY_PATH)
    entry = library["tRNA-Gln-TTG-2-1"]
    assert entry["amino_acid"] == "Gln"
    assert entry["anticodon"] == "UUG"
    assert entry["sequence"] == _KNOWN_SEQUENCE


def test_apply_trna_type_fills_sequence_anticodon_wobble_position():
    library = load_trna_library(LIBRARY_PATH)
    config = _make_config(trna_type="tRNA-Gln-TTG-2-1")
    warnings: list[dict] = []

    apply_trna_type(config, warnings, library)

    assert config.sequence["sequence"] == _KNOWN_SEQUENCE
    assert config.sequence["anticodon"] == "UUG"
    assert config.sequence["wobble_position"] == 35
    assert not any(w["Level"] == "ERROR" for w in warnings)


def test_apply_trna_type_unknown_id_records_error_and_leaves_sequence_unchanged():
    library = load_trna_library(LIBRARY_PATH)
    config = _make_config(trna_type="tRNA-Does-Not-Exist-1-1")
    warnings: list[dict] = []

    apply_trna_type(config, warnings, library)

    assert config.sequence["sequence"] == ""
    assert any(w["Level"] == "ERROR" and w["Source"] == "trna_library" for w in warnings)


def test_apply_trna_type_noop_when_trna_type_empty():
    library = load_trna_library(LIBRARY_PATH)
    config = _make_config(trna_type="", sequence="MANUAL_SEQUENCE")
    warnings: list[dict] = []

    apply_trna_type(config, warnings, library)

    assert config.sequence["sequence"] == "MANUAL_SEQUENCE"
    assert warnings == []


def test_apply_trna_type_warns_when_overwriting_conflicting_manual_values():
    library = load_trna_library(LIBRARY_PATH)
    config = _make_config(trna_type="tRNA-Gln-TTG-2-1", sequence="WRONG_SEQUENCE", anticodon="AAA")
    warnings: list[dict] = []

    apply_trna_type(config, warnings, library)

    assert config.sequence["sequence"] == _KNOWN_SEQUENCE
    assert config.sequence["anticodon"] == "UUG"
    warning_messages = [w["Message"] for w in warnings if w["Level"] == "WARNING"]
    assert any("sequence.sequence" in message for message in warning_messages)
    assert any("sequence.anticodon" in message for message in warning_messages)
