"""End-to-end orchestration for the simple RNA_MassHunter pipeline.

New module (仕様書 §6, §8.4, §22 phases 3-4). This is the first place that
actually calls `cca_processing.process_cca_tail()` — it exists in the
vendored code but the original `main.py` never wires it in (see
cca_processing.py's module docstring and 仕様書 §3.4).

Deliberately independent of the original, much larger `main.py`: nothing
here imports it, and `main.py` itself is untouched (仕様書 §22, Claude Code
rule #1).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from rna_masshunter import config as config_module
from rna_masshunter.cca_processing import CCAProcessingResult, process_cca_tail
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.mass_comparison import MassComparisonRow, build_mass_comparison_rows
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses
from rna_masshunter.models import Fragment, Modification, Peak, RunConfig
from rna_masshunter.modifications import load_modifications, validate_modifications
from rna_masshunter.ms2_extraction import extract_ms2_spectra
from rna_masshunter.ms2_support import build_ms2_ion_index
from rna_masshunter.peak_picking import extract_ms1_peaks
from rna_masshunter.warnings_manager import add_warning

REPO_ROOT = Path(__file__).resolve().parent.parent


def run(config_path: str | Path, project_root: str | Path | None = None) -> dict[str, Any]:
    """仕様書 §8.4 の手順:

    1. config読込・パス解決・バリデーション
    2/3. 修飾YAML読込
    4. cca_processing.process_cca_tail() でCCA成熟化（配線: これが初めて
       process_cca_tail を呼び出す箇所）
    5/6. 理論質量計算
    7. digestion.digest_sequence() で理論断片生成
    8. peak_picking.extract_ms1_peaks() でMS1ピーク抽出（mzml_path未指定
       時はスキップし警告のみ記録、§18）
    9. config.ms2_annotation.enabled かつ mzMLがある場合のみMS2スペクトル
       抽出・理論イオンindex構築
    10. mass_comparison.build_mass_comparison_rows() で行生成
    11. config.reporting.excel_output が真の場合のみExcel出力
        （excel_report.py はPhase 5で追加されるため、ここでは遅延import —
        呼ばれるのは実際にExcel出力が有効な場合のみ）
    12. テスト・CLI双方から使えるよう dict で結果を返す
    """
    root = Path(project_root) if project_root is not None else REPO_ROOT
    warnings: list[dict[str, Any]] = []

    config: RunConfig = config_module.load_config(config_path, warnings)
    config = config_module.resolve_paths(config, root)
    config_module.validate_config(config, warnings)

    modifications: list[Modification] = load_modifications(root / "data" / "modifications.yaml", warnings)
    validate_modifications(modifications, warnings)
    base_masses = load_base_masses(root / "data" / "base_masses.yaml", warnings)

    raw_sequence = str(config.sequence.get("sequence") or "")
    fragments: list[Fragment] = []
    theoretical_mass: float | None = None
    cca_result: CCAProcessingResult | None = None
    if not raw_sequence:
        add_warning(warnings, "WARNING", "simple_pipeline", "sequence.sequence is empty; theoretical mass and fragments were not generated.")
    else:
        cca_result = process_cca_tail(
            raw_sequence,
            config.cca_processing.get("registered_sequence_cca_mode", "EXCLUDES_CCA"),
            enabled=bool(config.cca_processing.get("enabled", True)),
        )
        theoretical_mass = calculate_unmodified_rna_mass(cca_result.processed_sequence, base_masses, warnings)
        position_map = {i: i for i in range(1, len(cca_result.processed_sequence) + 1)}
        fragments = digest_sequence(
            target_id=str(config.sequence.get("name") or "target"),
            sequence=cca_result.processed_sequence,
            position_map=position_map,
            config=config,
            base_masses=base_masses,
            warnings=warnings,
        )

    mzml_path = str(config.input.get("mzml_path") or "")
    peaks: list[Peak] = []
    if not mzml_path:
        add_warning(warnings, "WARNING", "simple_pipeline", "input.mzml_path is empty; MS1 peaks were not extracted.")
    else:
        try:
            peaks = extract_ms1_peaks(mzml_path, config.ms1_peak_extraction, warnings)
        except Exception as exc:  # noqa: BLE001 - mirrors ms2_extraction's "skip with warning" pattern
            add_warning(warnings, "WARNING", "simple_pipeline", "MS1 peak extraction failed; Mass Comparison will be skipped.", str(exc))
            peaks = []

    ms2_spectra: list[Any] = []
    ms2_ion_index: dict[str, list[Any]] = {}
    if mzml_path and fragments and bool(config.ms2_annotation.get("enabled", True)):
        ms2_spectra = extract_ms2_spectra(mzml_path, config.ms2_annotation, warnings)
        if ms2_spectra:
            ms2_ion_index = build_ms2_ion_index(fragments, config, base_masses, warnings)

    rows: list[MassComparisonRow] = build_mass_comparison_rows(
        fragments, peaks, modifications, config, warnings=warnings,
        ms2_spectra=ms2_spectra or None, ms2_ion_index=ms2_ion_index or None,
    )

    output_path: Path | None = None
    if bool(config.reporting.get("excel_output", True)):
        from rna_masshunter.excel_report import write_simple_mass_hunter_report

        output_dir = Path(str(config.project.get("output_dir") or "output"))
        output_dir.mkdir(parents=True, exist_ok=True)
        output_filename = str(config.reporting.get("output_filename") or "RNA_MassHunter_simple_report.xlsx")
        output_path = output_dir / output_filename
        write_simple_mass_hunter_report(output_path, config, fragments, peaks, rows, modifications, warnings=warnings)

    return {
        "config": config,
        "warnings": warnings,
        "cca_result": cca_result,
        "theoretical_mass": theoretical_mass,
        "fragments": fragments,
        "peaks": peaks,
        "ms2_spectra": ms2_spectra,
        "mass_comparison_rows": rows,
        "output_path": str(output_path) if output_path else None,
    }
