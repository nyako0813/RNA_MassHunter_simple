# RNA_MassHunter (simple)

RNA/tRNA LC-MS/MSデータから、理論断片・観測MS1質量・元素組成候補(Formula Candidate)・既知修飾候補(Modification Candidate)をExcelに整理する簡易版パイプライン。

自動での修飾同定は目的にせず、MSデータを正確に整理し人間が判断できる情報を提供することを目的とする。詳細は `docs/design/RNA_MassHunter_再設計_実装仕様書.md` を参照。

## 実装状況

設計仕様書（`docs/design/RNA_MassHunter_再設計_実装仕様書.md`）§22のPhase 0〜9まで実装・検証済み（`pytest` 67件通過、実データでのエンドツーエンド動作確認済み）。

`nyako0813/RNA_MassHunter`（大規模な既存リポジトリ）から移植した機能:

- 配列処理: CCA末端処理 (`rna_masshunter/cca_processing.py`)
- 消化: RNase A/T1・missed cleavage・Alkaline Phosphatase (`rna_masshunter/digestion.py`, `enzymes.py`)
- 質量計算: 理論質量・観測中性質量・元素質量 (`rna_masshunter/masses.py`, `elemental_composition.py`)
- MS1: mzML読み込み・ピーク抽出 (`rna_masshunter/mzml_reader.py`, `peak_picking.py`, `mzml_diagnostics.py`)
- MS1マッチング: 二分探索エンジン、ppm誤差計算 (`rna_masshunter/mass_shift_ms1_search.py`, `ms1_mapping.py`)
- 既知修飾: `modifications.yaml` ロード・質量差検索 (`rna_masshunter/modifications.py`)
- MS2参考情報: スペクトル抽出・理論d/w/a/zイオン生成 (`rna_masshunter/ms2_extraction.py`)
- 設定: `config.py` / `config.yaml`（このプロジェクト用に絞った新規セクション構成）

`rna_masshunter/ms2_extraction.py` と `rna_masshunter/cca_processing.py` は元ファイルをそのまま import せず、手を入れてコピーしている（元の `ms2_annotation.py` / `cca_tail_state.py` は本プロジェクトが使わない巨大な依存チェーンを引き込むため）。詳細は各ファイル冒頭のdocstringを参照。

このプロジェクトで新規実装したもの:

- `rna_masshunter/formula_candidate.py`: ΔDaから元素組成候補(Formula Candidate)を網羅探索・ランキング（原子数バイアス無し）。内部実装は、`max_total_atoms`/`element_limits`/`elements`の組み合わせごとに「到達可能な組成の全体集合（ユニバース）」を初回のみDFSで構築してキャッシュし、以降の各ΔDa問い合わせはその質量ソート済みユニバースへの二分探索（bisect）で応答する。実データで1クエリあたり約5マイクロ秒（旧・毎回DFS方式の約3000倍高速）。
- `rna_masshunter/mass_comparison.py`, `observed_mass.py`: Fragment×charge×PeakマッチングとMass Comparison行生成。
  - ピーク探索は理論m/zの狭いppm窓ではなく `mass_comparison.max_delta_da`（絶対質量幅、既定300Da）を主フィルタとする（実際の修飾質量シフトを見逃さないため）。
  - `max_matches_per_fragment`による切り詰めの優先順位は「Formula/Modification候補が1件以上ある行を優先 → 同順位内は|ΔDa|昇順」（探索窓が広いため、候補の有無を無視した|ΔDa|のみの切り詰めだと、質量的に近いだけで情報の無い行が枠を占有してしまう問題への対応）。
  - `has_any_candidate()`で候補有無を判定する公開ヘルパーを提供。
- `rna_masshunter/ms2_support.py`: MS2参考情報（未修飾d/w/a/zイオンとの一致件数・一覧）の付加。Recommended Formula/Known Modificationの判定には一切影響しない。
- `rna_masshunter/peak_picking.py`（既存ファイルへの追加関数、§14B/§14C）:
  - `merge_adjacent_profile_points()`: mzMLがcentroid化されていないprofile modeの場合に、同一scan内でppm近接する分裂点をintensity最大の点に統合する。
  - `merge_peaks_across_scans()`: 同一イオンがLC溶出中に複数の連続scanにまたがって検出される場合（`merge_max_scan_gap`以内のscan順位差）に、1行へ統合する。統合による情報欠落を避けるため`Peak.scan_count`/`Peak.rt_range`に集約情報を保持する。
  - いずれも`config.ms1_peak_extraction`の`merge_profile_points`/`merge_across_scans`で個別に有効・無効を切り替え可能（既定は両方true）。
- `rna_masshunter/simple_pipeline.py`, `simple_main.py`: 全体オーケストレーション（`main.py` は変更していない、独立した新規エントリポイント）。ピーク抽出直後・Mass Comparisonより前に上記の統合処理を適用する。
- `rna_masshunter/excel_report.py`: 01_Index〜08_Visualizationのシート出力（元リポジトリの同名だが無関係な `excel_report.py` とは別物）。
  - 04_Observed_Mass/05_Mass_Intensity/06_Mass_Comparisonに「Scan Count」「RT Range」列（§14Cの統合情報）。
  - 08_Visualizationは06_Mass_Comparisonを情報源とするCharge×ΔDaの散布図で、Formula/Modification候補の有無で2系列に色分けする。

## 実データでの検証結果（Phase 8）

ユーザーの実mzML（tRNA-Gln2、profile-mode、2,167 MS1スペクトル、MS2スペクトルなし）で`simple_main.py --config config.yaml`を実行し、以下を確認済み:

- 実行時間: 約1分17秒
- profile-mode分裂点の統合（§14B）でMS1ピーク数 75,686 → 18,738件（-75%）
- 複数scanにまたがる同一イオン検出の統合（§14C、`merge_max_scan_gap=2`、実データの検出間隔分布から妥当性を確認済み）
- `06_Mass_Comparison`でFormula/Modification候補を持つ行の割合が切り詰め優先順位変更により55.6% → 82.5%に改善
- 01_Indexの全リンク・各シートの逆リンクが機能
- 08_Visualizationのチャートが正しく埋め込まれる

MS2列については、このmzMLがMS1専用データ（MS2スペクトルを含まない）だったため常に空欄——実装の不具合ではなく、データの性質を正しく反映した結果。

## セットアップ

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
pytest
```

`config.yaml` を編集して `sequence.sequence` と `input.mzml_path` を設定してから実行する。

```bash
python simple_main.py --config config.yaml
```

`output/RNA_MassHunter_simple_report.xlsx` が生成される。`data/input/`（実mzML等の生データ）と`output/`（生成物）は`.gitignore`でリポジトリから除外されている。
