"""Excel report writer for the simple RNA_MassHunter pipeline.

New file, NOT the same module as the original repository's much larger
`excel_report.py` (仕様書 §5 — a same-named-but-different file is
intentional: `write_excel_report()`'s 670-line body and its Intact/SCIEX/
shadow-audit sheet set target a different, unrelated output format). This
file only contains the new `write_simple_mass_hunter_report()` for the
01-08 sheet layout (仕様書 §16).

The handful of small formatting/index helpers below
(`_flatten_dict`, `_autosize_and_freeze`, `_sheet_link`, `_coerce_to_frame`,
`_json_safe`, `_excel_safe_cell`, `_add_index_and_backlinks`,
`_truncate_frame_if_needed`, `_append_excel_warning`, and the
`EXCEL_MAX_ROWS`/`DATA_START_ROW`/`EXCEL_DATA_ROW_LIMIT` constants) are
vendored from nyako0813/RNA_MassHunter `rna_masshunter/excel_report.py`
(仕様書 §4, §16.3 — reuse the existing Index/hyperlink/formatting
conventions rather than reimplementing them) — fetched read-only from the
original repo for reference, the same pattern already used for
models.py/config.py/ms2_extraction.py (see those files' docstrings). Only
these specific helpers were copied, not the surrounding 2314-line file or
any of the SCIEX/shadow-audit sheet-building logic.

One deliberate adaptation: the original hardcodes the index sheet's name as
"Index"; this project's sheet-naming convention numbers sheets ("01_Index",
"02_Input", ...), so `_add_index_and_backlinks` / `_autosize_and_freeze`
below take an `index_sheet_name` parameter instead of the literal string.
"""
from __future__ import annotations

import json
import math
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter

from rna_masshunter.mass_comparison import MassComparisonRow
from rna_masshunter.models import Fragment, Modification, Peak, RunConfig
from rna_masshunter.observed_mass import assign_peak_ids

# --- vendored helpers (see module docstring) ---------------------------------

EXCEL_MAX_ROWS = 1_048_576
DATA_START_ROW = 3
EXCEL_DATA_ROW_LIMIT = EXCEL_MAX_ROWS - DATA_START_ROW


def _flatten_dict(data: dict[str, Any], prefix: str = "") -> list[dict[str, Any]]:
    rows = []
    for key, value in data.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            rows.extend(_flatten_dict(value, full_key))
        else:
            rows.append({"Parameter": full_key, "Value": value})
    return rows


def _autosize_and_freeze(writer: pd.ExcelWriter, index_sheet_name: str = "01_Index") -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2" if worksheet.title == index_sheet_name else "A4"
        for column_cells in worksheet.columns:
            max_length = 0
            column = get_column_letter(column_cells[0].column)
            for cell in column_cells:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, min(len(value), 60))
            worksheet.column_dimensions[column].width = max(10, max_length + 2)


def _sheet_link(sheet_name: str, cell: str = "A1") -> str:
    return f"#'{sheet_name}'!{cell}"


def _coerce_to_frame(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list):
        return pd.DataFrame(value)
    if isinstance(value, dict):
        return pd.DataFrame([value])
    return pd.DataFrame([{"Value": value}])


def _json_safe(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, Enum):
        return _json_safe(value.value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "to_dict") and callable(value.to_dict):
        return _json_safe(value.to_dict())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        normalized = [_json_safe(item) for item in value]
        return sorted(normalized, key=lambda item: json.dumps(item, ensure_ascii=False, sort_keys=True, default=str))
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _excel_safe_cell(value: Any) -> Any:
    normalized = _json_safe(value)
    if normalized is None:
        return ""
    if isinstance(normalized, (list, dict)):
        return json.dumps(normalized, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return normalized


def _add_index_and_backlinks(writer: pd.ExcelWriter, sheet_names: list[str], index_sheet_name: str = "01_Index") -> None:
    workbook = writer.book
    index_sheet = workbook[index_sheet_name]
    for row_index, sheet_name in enumerate(sheet_names, start=2):
        link_cell = index_sheet.cell(row=row_index, column=1)
        link_cell.value = sheet_name
        link_cell.hyperlink = _sheet_link(sheet_name, "A1")
        link_cell.style = "Hyperlink"

    for sheet_name in sheet_names:
        worksheet = workbook[sheet_name]
        worksheet["A1"] = "← Back to Index"
        worksheet["A1"].hyperlink = _sheet_link(index_sheet_name, "A1")
        worksheet["A1"].style = "Hyperlink"


def _append_excel_warning(warnings: list[dict[str, Any]], sheet_name: str, original_rows: int, written_rows: int) -> None:
    warnings.append({
        "Timestamp": datetime.now().isoformat(timespec="seconds"),
        "Level": "WARNING",
        "Source": "excel_report",
        "Message": "Excel sheet was truncated because it exceeded max_excel_rows_per_sheet.",
        "Context": {"sheet": sheet_name, "original_rows": original_rows, "written_rows": written_rows},
    })


def _truncate_frame_if_needed(
    sheet_name: str,
    frame: pd.DataFrame,
    max_rows: int,
    truncate_large_sheets: bool,
    warnings: list[dict[str, Any]],
) -> pd.DataFrame:
    original_rows = len(frame)
    safe_limit = min(max_rows, EXCEL_DATA_ROW_LIMIT)
    if original_rows <= safe_limit:
        return frame
    written_rows = safe_limit if truncate_large_sheets else EXCEL_DATA_ROW_LIMIT
    written_rows = min(written_rows, original_rows)
    _append_excel_warning(warnings, sheet_name, original_rows, written_rows)
    return frame.head(written_rows).copy()


# --- sheet-specific helpers (new) --------------------------------------------

SHEET_DESCRIPTIONS = {
    "02_Input": "Run configuration (flattened) and any warnings raised during this run.",
    "03_Theoretical": "Theoretical RNase digestion fragments and their unmodified monoisotopic mass.",
    "04_Observed_Mass": "Observed MS1 peaks with the charge(s) under which each peak matched a fragment in 06_Mass_Comparison.",
    "05_Mass_Intensity": "Same as 04_Observed_Mass with Intensity as an explicit column.",
    "06_Mass_Comparison": "Fragment x charge x Peak matches: mass differences plus independent Formula/Modification Candidate columns.",
    "07_Modifications": "Known RNA modification database (data/modifications.yaml) used for Modification Candidate search.",
    "08_Visualization": "Scatter chart: Charge (x) vs Observed Neutral Mass (y), sourced from 04_Observed_Mass.",
}


def _modification_display_name(modification: Modification) -> str:
    raw = modification.raw or {}
    return str(raw.get("name") or modification.symbol or modification.id)


def _theoretical_frame(fragments: list[Fragment]) -> pd.DataFrame:
    rows = [
        {
            "Fragment ID": fragment.fragment_id,
            "Sequence": fragment.sequence,
            "Start": fragment.start,
            "End": fragment.end,
            "Length": fragment.end - fragment.start + 1,
            "Enzyme": fragment.enzyme,
            "Missed Cleavage": fragment.missed_cleavages,
            "Terminal Form": fragment.terminal_form,
            "Theoretical Mass": fragment.unmodified_mass,
        }
        for fragment in fragments
    ]
    return pd.DataFrame(rows, columns=["Fragment ID", "Sequence", "Start", "End", "Length", "Enzyme", "Missed Cleavage", "Terminal Form", "Theoretical Mass"])


def _format_rt_range(scan_count: int, rt_range: tuple[float, float] | None) -> str:
    """仕様書 §14C: scan_count<=1（cross-scanマージが起きていない通常の
    ピーク）なら空欄のまま。既存のRT列と重複する情報を増やさないため。"""
    if scan_count <= 1 or rt_range is None:
        return ""
    return f"{rt_range[0]:.3f}-{rt_range[1]:.3f}"


def _observed_rows(peaks: list[Peak], mass_comparison_rows: list[MassComparisonRow]) -> list[dict[str, Any]]:
    """仕様書 §15/§16.2 の推奨に従い、1ピークにつき複数chargeでマッチした
    場合はcharge単位で行を分ける。マッチが1件も無かったピークも
    (charge/observed massは空欄のまま) 1行として残し、抽出済みMS1ピーク
    全体が04/05シートから欠落しないようにする。Scan Count/RT Range列は
    peak_picking.merge_peaks_across_scans() が統合したピークについてのみ
    埋まる（仕様書 §14C）。"""
    ids_by_index = assign_peak_ids(peaks)
    charges_by_peak_id: dict[str, set[int]] = {}
    observed_mass_by_peak_charge: dict[tuple[str, int], float] = {}
    for row in mass_comparison_rows:
        charges_by_peak_id.setdefault(row.peak_id, set()).add(row.charge)
        observed_mass_by_peak_charge.setdefault((row.peak_id, row.charge), row.observed_mass)

    rows: list[dict[str, Any]] = []
    for index, peak in enumerate(peaks):
        peak_id = ids_by_index[index]
        scan_count = getattr(peak, "scan_count", 1)
        rt_range_str = _format_rt_range(scan_count, getattr(peak, "rt_range", None))
        charges = sorted(charges_by_peak_id.get(peak_id, set()))
        if not charges:
            rows.append({
                "Peak ID": peak_id, "m/z": peak.mz, "Charge": None, "Observed Mass": None,
                "Intensity": peak.intensity, "RT": peak.rt, "Scan ID": peak.scan_id,
                "Scan Count": scan_count, "RT Range": rt_range_str,
            })
            continue
        for charge in charges:
            rows.append({
                "Peak ID": peak_id, "m/z": peak.mz, "Charge": charge,
                "Observed Mass": observed_mass_by_peak_charge.get((peak_id, charge)),
                "Intensity": peak.intensity, "RT": peak.rt, "Scan ID": peak.scan_id,
                "Scan Count": scan_count, "RT Range": rt_range_str,
            })
    return rows


def _observed_mass_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    columns = ["Peak ID", "m/z", "Charge", "Observed Mass", "RT", "Scan ID", "Scan Count", "RT Range"]
    return pd.DataFrame(rows, columns=columns)


def _mass_intensity_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    columns = ["Peak ID", "m/z", "Charge", "Observed Mass", "Intensity", "RT", "Scan ID", "Scan Count", "RT Range"]
    return pd.DataFrame(rows, columns=columns)


_MASS_COMPARISON_COLUMNS = [
    "Peak ID", "Fragment ID", "Sequence", "Charge", "Intensity", "Observed Mass", "Theoretical Mass",
    "ΔDa", "Δppm", "Recommended Formula", "Formula Candidates", "Known Modification",
    "Modification Candidates", "MS2 Spectrum ID", "MS2 Matched Ions", "Scan Count", "RT Range",
]


def _mass_comparison_frame(rows: list[MassComparisonRow]) -> pd.DataFrame:
    data = [
        {
            "Peak ID": row.peak_id, "Fragment ID": row.fragment_id, "Sequence": row.sequence,
            "Charge": row.charge, "Intensity": row.intensity, "Observed Mass": row.observed_mass,
            "Theoretical Mass": row.theoretical_mass, "ΔDa": row.delta_da, "Δppm": row.delta_ppm,
            "Recommended Formula": row.recommended_formula or "", "Formula Candidates": row.formula_candidates,
            "Known Modification": row.known_modification or "", "Modification Candidates": row.modification_candidates,
            "MS2 Spectrum ID": row.ms2_spectrum_id or "", "MS2 Matched Ions": row.ms2_matched_ions,
            "Scan Count": row.scan_count, "RT Range": _format_rt_range(row.scan_count, row.rt_range),
        }
        for row in rows
    ]
    return pd.DataFrame(data, columns=_MASS_COMPARISON_COLUMNS)


def _modifications_frame(modifications: list[Modification]) -> pd.DataFrame:
    rows = [
        {
            "ID": modification.id, "Symbol": modification.symbol or "", "Name": _modification_display_name(modification),
            "Target Bases": ", ".join(modification.target_bases), "Category": modification.category,
            "Mass Shift (Da)": modification.mass_shift_from_unmodified, "Chemical Group": modification.chemical_group,
            "Near-Isobaric Group": modification.near_isobaric_group, "Detectability": modification.detectability,
            "Curation Status": modification.curation_status,
        }
        for modification in modifications
    ]
    columns = ["ID", "Symbol", "Name", "Target Bases", "Category", "Mass Shift (Da)", "Chemical Group", "Near-Isobaric Group", "Detectability", "Curation Status"]
    return pd.DataFrame(rows, columns=columns)


def _input_frame(config: RunConfig, warnings: list[dict[str, Any]]) -> pd.DataFrame:
    config_dict = {key: value for key, value in vars(config).items() if key != "raw"}
    rows = _flatten_dict(config_dict)
    if warnings:
        rows.append({"Parameter": "--- Warnings ---", "Value": f"{len(warnings)} warning(s)"})
        for warning in warnings:
            context = warning.get("Context")
            suffix = f" ({context!r})" if context not in (None, "") else ""
            rows.append({"Parameter": f"[{warning['Level']}] {warning['Source']}", "Value": f"{warning['Message']}{suffix}"})
    return pd.DataFrame(rows, columns=["Parameter", "Value"])


def _add_visualization_sheet(
    writer: pd.ExcelWriter,
    observed_mass_frame: pd.DataFrame,
    config: RunConfig,
    warnings: list[dict[str, Any]],
) -> None:
    """仕様書 §17: 08_Visualization に Charge(x) x Observed Neutral Mass(y)
    の散布図を1枚配置する。04_Observed_Mass シートのセルを直接参照するので
    (openpyxlの`Reference`は値のコピーではなくセル参照)、Excel上でデータ
    が更新されればチャートも追随する。チャート生成が失敗しても例外を
    握りつぶし、レポート全体の出力は継続する（§17, §18）。"""
    sheet = writer.book.create_sheet("08_Visualization")
    sheet["A1"] = "← Back to Index"
    sheet["A1"].hyperlink = _sheet_link("01_Index", "A1")
    sheet["A1"].style = "Hyperlink"

    viz_config = config.visualization or {}
    if not viz_config.get("enabled", True):
        sheet["A3"] = "Visualization is disabled (config.visualization.enabled = false)."
        return
    if observed_mass_frame.empty:
        sheet["A3"] = "No observed peaks to plot."
        return

    try:
        from openpyxl.chart import Reference, ScatterChart, Series

        source_sheet = writer.book["04_Observed_Mass"]
        columns = list(observed_mass_frame.columns)
        charge_col = columns.index("Charge") + 1
        mass_col = columns.index("Observed Mass") + 1
        data_first_row = DATA_START_ROW + 1  # header at DATA_START_ROW, data starts the row after
        data_last_row = DATA_START_ROW + len(observed_mass_frame)

        chart = ScatterChart()
        chart.title = "Observed Neutral Mass by Charge"
        chart.x_axis.title = "Charge"
        chart.y_axis.title = "Observed Neutral Mass (Da)"
        chart.style = 2

        x_values = Reference(source_sheet, min_col=charge_col, min_row=data_first_row, max_row=data_last_row)
        y_values = Reference(source_sheet, min_col=mass_col, min_row=data_first_row, max_row=data_last_row)
        series = Series(y_values, x_values, title="Observed peaks")
        series.marker.symbol = "circle"
        series.graphicalProperties.line.noFill = True
        chart.series.append(series)

        sheet.add_chart(chart, "A5")
    except Exception as exc:  # noqa: BLE001 - chart generation must never abort the whole report
        from rna_masshunter.warnings_manager import add_warning

        add_warning(warnings, "WARNING", "excel_report", "08_Visualization chart generation failed; report was written without it.", str(exc))
        sheet["A3"] = "Chart generation failed; see Warnings in 02_Input."


def _apply_number_formats(worksheet: Any, frame: pd.DataFrame, formats: dict[str, str], header_row: int = 3) -> None:
    if frame.empty:
        return
    for column_name, number_format in formats.items():
        if column_name not in frame.columns:
            continue
        column_index = list(frame.columns).index(column_name) + 1
        column_letter = get_column_letter(column_index)
        for row_offset in range(len(frame)):
            cell = worksheet[f"{column_letter}{header_row + 1 + row_offset}"]
            cell.number_format = number_format


def write_simple_mass_hunter_report(
    output_path: str | Path,
    config: RunConfig,
    fragments: list[Fragment],
    peaks: list[Peak],
    mass_comparison_rows: list[MassComparisonRow],
    modifications: list[Modification],
    warnings: list[dict[str, Any]] | None = None,
) -> Path:
    """仕様書 §16: 01_Index〜08_Visualization を出力する。既存の
    _add_index_and_backlinks/_autosize_and_freeze のIndex/ハイパーリンク/
    書式の流儀をそのまま踏襲する（§16.3）。"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    warnings = warnings if warnings is not None else []

    reporting = config.reporting or {}
    max_rows = int(reporting.get("max_excel_rows_per_sheet", 100000) or 100000)
    truncate_large_sheets = bool(reporting.get("truncate_large_sheets", True))

    observed_rows = _observed_rows(peaks, mass_comparison_rows)

    sheets: dict[str, pd.DataFrame] = {
        "02_Input": _input_frame(config, warnings),
        "03_Theoretical": _theoretical_frame(fragments),
        "04_Observed_Mass": _observed_mass_frame(observed_rows),
        "05_Mass_Intensity": _mass_intensity_frame(observed_rows),
        "06_Mass_Comparison": _mass_comparison_frame(mass_comparison_rows),
        "07_Modifications": _modifications_frame(modifications),
    }
    sheets = {name: _truncate_frame_if_needed(name, frame, max_rows, truncate_large_sheets, warnings) for name, frame in sheets.items()}
    # Cells may hold None / dict (e.g. Detectability) — normalize them so
    # openpyxl never chokes on an unsupported type (仕様書 §18 error handling
    # spirit: prefer a readable blank/JSON string over a crash).
    sheets = {name: frame.map(_excel_safe_cell) if not frame.empty else frame for name, frame in sheets.items()}

    all_sheet_names = list(sheets) + ["08_Visualization"]
    index_rows = [{"Sheet": name, "Description": SHEET_DESCRIPTIONS.get(name, ""), "Notes": "Data starts at A3."} for name in all_sheet_names]

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        pd.DataFrame(index_rows, columns=["Sheet", "Description", "Notes"]).to_excel(writer, sheet_name="01_Index", index=False)
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
        _add_visualization_sheet(writer, sheets["04_Observed_Mass"], config, warnings)
        _add_index_and_backlinks(writer, all_sheet_names, index_sheet_name="01_Index")

        # 仕様書 §16.4: mass系列は"0.000"、ΔDa系列は符号付き"+0.00000;-0.00000"、Δppmは"0.0"。
        _apply_number_formats(writer.sheets["03_Theoretical"], sheets["03_Theoretical"], {"Theoretical Mass": "0.000"})
        _apply_number_formats(writer.sheets["04_Observed_Mass"], sheets["04_Observed_Mass"], {"m/z": "0.000", "Observed Mass": "0.000"})
        _apply_number_formats(writer.sheets["05_Mass_Intensity"], sheets["05_Mass_Intensity"], {"m/z": "0.000", "Observed Mass": "0.000"})
        _apply_number_formats(writer.sheets["06_Mass_Comparison"], sheets["06_Mass_Comparison"], {
            "Observed Mass": "0.000", "Theoretical Mass": "0.000",
            "ΔDa": "+0.00000;-0.00000", "Δppm": "0.0",
        })

        _autosize_and_freeze(writer, index_sheet_name="01_Index")

    return output_path
