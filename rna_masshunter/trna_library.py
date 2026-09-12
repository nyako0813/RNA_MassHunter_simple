"""tRNA種類選択による配列自動入力機能 (tools/tRNA種類選択による配列自動入力機能_実装仕様書.md §4).

`config.sequence.trna_type` に data/trna_library.yaml の id を指定すると、
`sequence`/`anticodon`/`wobble_position` が自動入力される。空文字列のままなら
何もせず、従来どおり手入力の `sequence.sequence` が使われる（後方互換）。
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from rna_masshunter import config as config_module
from rna_masshunter.models import RunConfig
from rna_masshunter.warnings_manager import add_warning


def load_trna_library(path: str | Path) -> dict[str, dict[str, Any]]:
    """data/trna_library.yaml を読み込み、id -> entry の辞書を返す。"""
    source = Path(path)
    if not source.is_file():
        return {}
    data = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    entries = data.get("trna_library") or []
    return {str(entry["id"]): entry for entry in entries if entry.get("id")}


def apply_trna_type(config: RunConfig, warnings: list[dict[str, Any]] | None, library: dict[str, dict[str, Any]]) -> None:
    """config.sequence.trna_type が設定されていれば、配列・アンチコドン・wobble位を自動入力する。"""
    trna_type = str((config.sequence or {}).get("trna_type") or "").strip()
    if not trna_type:
        return

    entry = library.get(trna_type)
    if entry is None:
        available = ", ".join(sorted(library.keys())) or "(ライブラリが空です)"
        message = f"sequence.trna_type '{trna_type}' は data/trna_library.yaml に見つかりません。利用可能: {available}"
        if warnings is not None:
            add_warning(warnings, "ERROR", "trna_library", message)
        return

    manual_sequence = str((config.sequence or {}).get("sequence") or "").strip()
    manual_anticodon = str((config.sequence or {}).get("anticodon") or "").strip()

    if manual_sequence and manual_sequence != entry["sequence"]:
        if warnings is not None:
            add_warning(
                warnings, "WARNING", "trna_library",
                f"sequence.trna_type '{trna_type}' の配列で sequence.sequence の手入力値を上書きしました。",
            )
    if manual_anticodon and manual_anticodon != entry["anticodon"]:
        if warnings is not None:
            add_warning(
                warnings, "WARNING", "trna_library",
                f"sequence.trna_type '{trna_type}' のアンチコドンで sequence.anticodon の手入力値を上書きしました。",
            )

    config.sequence["sequence"] = entry["sequence"]
    config.sequence["anticodon"] = entry["anticodon"]
    config.sequence["wobble_position"] = entry["wobble_position"]
    default_name = config_module.DEFAULT_CONFIG["sequence"]["name"]
    if not str(config.sequence.get("name") or "").strip() or config.sequence.get("name") == default_name:
        config.sequence["name"] = entry["id"]
