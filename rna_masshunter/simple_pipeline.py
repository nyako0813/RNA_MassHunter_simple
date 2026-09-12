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
from rna_masshunter.nucleoside_comparison import NucleosideComparisonRow, build_nucleoside_comparison_rows
from rna_masshunter.nucleoside_targets import NucleosideTarget, build_nucleoside_target_universe
from rna_masshunter.peak_picking import extract_ms1_peaks, merge_adjacent_profile_points, merge_peaks_across_scans
from rna_masshunter.trna_library import apply_trna_type, load_trna_library
from rna_masshunter.warnings_manager import add_warning

# §24 (Phase 10): selecting this enzyme switches the whole pipeline into P1
# complete-digestion mode (known-nucleoside mass matching) instead of the
# oligomer fragment / ΔDa-exploration flow the rest of this module runs.
_P1_COMPLETE_DIGESTION_ENZYME = "Nuclease_P1"

REPO_ROOT = Path(__file__).resolve().parent.parent

# §24.5: P1 mode targets free nucleosides (~230-370 Da) up to their
# 5'-monophosphate form (~+80 Da, so up to roughly 450 Da) — a default
# oligomer-oriented mz_min like 500 would filter out essentially all of
# them before they ever reach nucleoside_comparison.py. This is just a
# heuristic trip-wire (not a hard validation error) to catch that
# misconfiguration early instead of silently returning zero matches.
_P1_MODE_MZ_MIN_WARNING_THRESHOLD = 400


def _warn_if_mz_min_too_high_for_p1_mode(config: RunConfig, warnings: list[dict[str, Any]]) -> None:
    mz_min = config.ms1_peak_extraction.get("mz_min")
    try:
        mz_min_value = float(mz_min) if mz_min is not None else None
    except (TypeError, ValueError):
        mz_min_value = None
    if mz_min_value is not None and mz_min_value >= _P1_MODE_MZ_MIN_WARNING_THRESHOLD:
        add_warning(
            warnings, "WARNING", "simple_pipeline",
            "ms1_peak_extraction.mz_min looks too high for Nuclease_P1 (P1 complete-digestion) mode; "
            "free nucleoside masses are typically ~230-370 Da (up to ~450 Da phosphorylated), so peaks "
            "may be filtered out before they can match. Consider lowering mz_min to ~100-150 for P1 mode.",
            {"mz_min": mz_min_value, "threshold": _P1_MODE_MZ_MIN_WARNING_THRESHOLD},
        )


def run(config_path: str | Path, project_root: str | Path | None = None) -> dict[str, Any]:
    """仕様書 §8.4 の手順:

    1. config読込・パス解決・バリデーション
    1.5. trna_library.apply_trna_type() で config.sequence.trna_type から
         sequence/anticodon/wobble_position を自動入力（tRNA種類選択機能、
         tools/tRNA種類選択による配列自動入力機能_実装仕様書.md §4.2）
    2/3. 修飾YAML読込
    4. cca_processing.process_cca_tail() でCCA成熟化（配線: これが初めて
       process_cca_tail を呼び出す箇所）
    5/6. 理論質量計算
    7. digestion.digest_sequence() で理論断片生成
    8. peak_picking.extract_ms1_peaks() でMS1ピーク抽出（mzml_path未指定
       時はスキップし警告のみ記録、§18）。
       config.ms1_peak_extraction.merge_profile_points が真の場合、続けて
       peak_picking.merge_adjacent_profile_points() でprofile-mode由来の
       分裂ピーク点を統合する（§14B）。以降の全ステップ（MS2前駆体マッチ、
       Mass Comparison、Excel出力）はこの統合後のpeaksを一貫して使う。
    9. config.ms2_annotation.enabled かつ mzMLがある場合のみMS2スペクトル
       抽出・理論イオンindex構築
    10. mass_comparison.build_mass_comparison_rows() で行生成
    11. config.reporting.excel_output が真の場合のみExcel出力
        （excel_report.py はPhase 5で追加されるため、ここでは遅延import —
        呼ばれるのは実際にExcel出力が有効な場合のみ）
    12. テスト・CLI双方から使えるよう dict で結果を返す

    §24（Phase 10）: `config.digestion.enzyme` が "Nuclease_P1" の場合は
    P1完全分解モードに分岐する。RNase A/T1等のオリゴマー断片フロー
    （配列・CCA処理・digest_sequence、いずれも位置依存のロジック）は
    完全にスキップし、代わりに既知ヌクレオシド質量ユニバース
    （nucleoside_targets.py）とMS1ピークを直接照合する
    （nucleoside_comparison.py、狭いppm許容差）。MS2参考情報・Formula
    Candidate探索はP1モードでは行わない（§24のスコープ外）。
    """
    root = Path(project_root) if project_root is not None else REPO_ROOT
    warnings: list[dict[str, Any]] = []

    config: RunConfig = config_module.load_config(config_path, warnings)
    config = config_module.resolve_paths(config, root)
    config_module.validate_config(config, warnings)

    trna_library = load_trna_library(root / "data" / "trna_library.yaml")
    apply_trna_type(config, warnings, trna_library)

    modifications: list[Modification] = load_modifications(root / "data" / "modifications.yaml", warnings)
    validate_modifications(modifications, warnings)
    base_masses = load_base_masses(root / "data" / "base_masses.yaml", warnings)

    is_p1_mode = str(config.digestion.get("enzyme") or "").strip() == _P1_COMPLETE_DIGESTION_ENZYME
    nucleoside_targets: list[NucleosideTarget] = []

    raw_sequence = str(config.sequence.get("sequence") or "")
    fragments: list[Fragment] = []
    theoretical_mass: float | None = None
    cca_result: CCAProcessingResult | None = None
    if is_p1_mode:
        nucleoside_targets = build_nucleoside_target_universe(modifications, warnings)
        if raw_sequence:
            add_warning(warnings, "INFO", "simple_pipeline", "digestion.enzyme is Nuclease_P1; sequence/CCA/fragment generation was skipped (P1 complete-digestion mode matches peaks directly against known nucleoside masses, §24).")
        if bool(config.ms2_annotation.get("enabled", True)):
            # §24.5: explicit, regardless of the config value — d/w/a/z
            # fragment ions require an oligomer backbone to cut, which does
            # not exist once P1 has hydrolyzed everything down to single
            # nucleosides, so MS2 annotation is not a meaningful concept in
            # this mode at all (not just "not yet implemented").
            add_warning(warnings, "INFO", "simple_pipeline", "digestion.enzyme is Nuclease_P1; MS2 annotation was explicitly skipped regardless of ms2_annotation.enabled — d/w/a/z fragment ions don't exist for single nucleosides (§24.5).")
        _warn_if_mz_min_too_high_for_p1_mode(config, warnings)
    elif not raw_sequence:
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
        merge_tolerance_ppm = float(config.ms1_peak_extraction.get("merge_tolerance_ppm", 10) or 10)
        if peaks and bool(config.ms1_peak_extraction.get("merge_profile_points", True)):
            before = len(peaks)
            peaks = merge_adjacent_profile_points(peaks, merge_tolerance_ppm)
            if len(peaks) != before:
                add_warning(warnings, "INFO", "simple_pipeline", "Merged adjacent profile-mode peak points (§14B).", {"before": before, "after": len(peaks)})
        if peaks and bool(config.ms1_peak_extraction.get("merge_across_scans", True)):
            before = len(peaks)
            max_scan_gap = int(config.ms1_peak_extraction.get("merge_max_scan_gap", 2) or 0)
            peaks = merge_peaks_across_scans(peaks, merge_tolerance_ppm, max_scan_gap)
            if len(peaks) != before:
                add_warning(warnings, "INFO", "simple_pipeline", "Merged the same ion's detections across nearby scans (§14C).", {"before": before, "after": len(peaks), "max_scan_gap": max_scan_gap})

    ms2_spectra: list[Any] = []
    ms2_ion_index: dict[str, list[Any]] = {}
    rows: list[MassComparisonRow] = []
    nucleoside_rows: list[NucleosideComparisonRow] = []
    if is_p1_mode:
        # §24: known-nucleoside matching only — no MS2, no Formula Candidate.
        nucleoside_rows = build_nucleoside_comparison_rows(peaks, nucleoside_targets, config, base_masses, warnings=warnings)
    else:
        if mzml_path and fragments and bool(config.ms2_annotation.get("enabled", True)):
            ms2_spectra = extract_ms2_spectra(mzml_path, config.ms2_annotation, warnings)
            if ms2_spectra:
                ms2_ion_index = build_ms2_ion_index(fragments, config, base_masses, warnings)

        rows = build_mass_comparison_rows(
            fragments, peaks, modifications, config, warnings=warnings,
            ms2_spectra=ms2_spectra or None, ms2_ion_index=ms2_ion_index or None,
        )

    output_path: Path | None = None
    if bool(config.reporting.get("excel_output", True)):
        output_dir = Path(str(config.project.get("output_dir") or "output"))
        output_dir.mkdir(parents=True, exist_ok=True)
        output_filename = str(config.reporting.get("output_filename") or "RNA_MassHunter_simple_report.xlsx")
        output_path = output_dir / output_filename

        if is_p1_mode:
            from rna_masshunter.excel_report import write_nucleoside_mass_hunter_report

            write_nucleoside_mass_hunter_report(output_path, config, nucleoside_targets, peaks, nucleoside_rows, modifications, base_masses, warnings=warnings)
        else:
            from rna_masshunter.excel_report import write_simple_mass_hunter_report

            write_simple_mass_hunter_report(output_path, config, fragments, peaks, rows, modifications, warnings=warnings)

    return {
        "config": config,
        "warnings": warnings,
        "mode": "p1_nucleoside" if is_p1_mode else "oligomer",
        "cca_result": cca_result,
        "theoretical_mass": theoretical_mass,
        "fragments": fragments,
        "peaks": peaks,
        "ms2_spectra": ms2_spectra,
        "mass_comparison_rows": rows,
        "nucleoside_targets": nucleoside_targets,
        "nucleoside_comparison_rows": nucleoside_rows,
        "output_path": str(output_path) if output_path else None,
    }
