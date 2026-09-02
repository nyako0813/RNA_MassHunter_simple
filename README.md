# RNA_MassHunter (simple)

RNA/tRNA LC-MS/MSデータから、理論断片・観測MS1質量・元素組成候補(Formula Candidate)・既知修飾候補(Modification Candidate)をExcelに整理する簡易版パイプライン。

自動での修飾同定は目的にせず、MSデータを正確に整理し人間が判断できる情報を提供することを目的とする。詳細は `docs/design/RNA_MassHunter_再設計_実装仕様書.md` を参照。

## この段階でできること（ベンダリング済み）

`nyako0813/RNA_MassHunter`（大規模な既存リポジトリ）から、このプロジェクトが必要とする機能だけを移植・検証済み:

- 配列処理: CCA末端処理 (`rna_masshunter/cca_processing.py`)
- 消化: RNase A/T1・missed cleavage・Alkaline Phosphatase (`rna_masshunter/digestion.py`, `enzymes.py`)
- 質量計算: 理論質量・観測中性質量・元素質量 (`rna_masshunter/masses.py`, `elemental_composition.py`)
- MS1: mzML読み込み・ピーク抽出 (`rna_masshunter/mzml_reader.py`, `peak_picking.py`, `mzml_diagnostics.py`)
- MS1マッチング: 二分探索エンジン、ppm誤差計算 (`rna_masshunter/mass_shift_ms1_search.py`, `ms1_mapping.py`)
- 既知修飾: `modifications.yaml` ロード・質量差検索 (`rna_masshunter/modifications.py`)
- MS2参考情報: スペクトル抽出・理論d/w/a/zイオン生成 (`rna_masshunter/ms2_extraction.py`)
- 設定: `config.py` / `config.yaml`（このプロジェクト用に絞った新規セクション構成）

`rna_masshunter/ms2_extraction.py` と `rna_masshunter/cca_processing.py` は元ファイルをそのまま import せず、手を入れてコピーしている（元の `ms2_annotation.py` / `cca_tail_state.py` は本プロジェクトが使わない巨大な依存チェーンを引き込むため）。詳細は各ファイル冒頭のdocstringを参照。

## まだ実装されていないこと（次のステップ）

`formula_candidate.py`, `mass_comparison.py`, `observed_mass.py`, `ms2_support.py`, `simple_pipeline.py`, `simple_main.py`, Excel出力(`excel_report.py`)は未実装。実装手順は `docs/design/Claude_Code_実装引き継ぎ指示書.md` を参照。

## セットアップ

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
pytest
```

`config.yaml` を編集して `sequence.sequence` と `input.mzml_path` を設定してから実行する。
