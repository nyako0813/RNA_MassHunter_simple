# RNA_MassHunter (simple)

RNA/tRNA LC-MS/MSデータから、理論断片・観測MS1質量・元素組成候補(Formula Candidate)・既知修飾候補(Modification Candidate)をExcelに整理する簡易版パイプライン。

自動での修飾同定は目的にせず、MSデータを正確に整理し人間が判断できる情報を提供することを目的とする。詳細は `docs/design/RNA_MassHunter_再設計_実装仕様書.md` を参照。

## 実装状況

設計仕様書（`docs/design/RNA_MassHunter_再設計_実装仕様書.md`）§22のPhase 0〜9、および§24（Phase 10、P1完全分解モード）まで実装・検証済み（`pytest` 167件通過、実データでのエンドツーエンド動作確認済み）。

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
- `rna_masshunter/simple_pipeline.py`, `main.py`: 全体オーケストレーション（`main.py`は旧名`simple_main.py`から改名した、本体RNA_MassHunterの`main.py`とは独立した新規エントリポイント）。ピーク抽出直後・Mass Comparisonより前に上記の統合処理を適用する。
- `rna_masshunter/excel_report.py`: 01_Index〜08_Visualizationのシート出力（元リポジトリの同名だが無関係な `excel_report.py` とは別物）。
  - 04_Observed_Mass/05_Mass_Intensity/06_Mass_Comparisonに「Scan Count」「RT Range」列（§14Cの統合情報）。
  - 08_Visualizationは06_Mass_Comparisonを情報源とするCharge×ΔDaの散布図で、Formula/Modification候補の有無で2系列に色分けする。

### P1完全分解モード（§24, Phase 10）

`config.digestion.enzyme` に `"Nuclease_P1"` を指定すると、パイプライン全体が別モードに切り替わる。Nuclease P1はRNAを個々のリボヌクレオシド（5'-一リン酸）まで完全に加水分解するため、オリゴマー断片×ΔDa探索ではなく、既知ヌクレオシド質量との直接照合を行う。

- `rna_masshunter/nucleoside_targets.py`: 標準4塩基（`ElementalComposition`から算出、文献値と照合するテスト付き）+ `data/modifications.yaml`収載の修飾ヌクレオシド（`modified_nucleoside_mass_mono`フィールドを使用、118件全てに値あり）から成る既知ヌクレオシド質量ユニバースを構築する。
- `rna_masshunter/nucleoside_comparison.py`: 既存の`mass_shift_ms1_search.find_peaks_near_mz`（狭いppm許容差）でMS1ピークをこのユニバースと直接照合する。`mass_comparison.py`の広い`max_delta_da`探索とは異なり、既知の個別質量への確認照合のため。
- リン酸の有無は既存の`config.alkaline_phosphatase.enabled`をそのまま流用（P1モード専用の新規configキーは増やしていない）。有効なら遊離ヌクレオシド質量、無効なら5'-一リン酸化分を加えた質量で照合する。
- Formula Candidateによる未知修飾ヌクレオシドの推定は行わない（既知ヌクレオシドとの質量照合のみが§24のスコープ）。MS2参考情報もP1モードでは付加しない——`config.ms2_annotation.enabled`の値に関わらず明示的にスキップする（d/w/a/zイオンは断片化する主鎖があって初めて意味を持つ概念であり、単一ヌクレオシドには成立しないため、§24.5）。
- P1モード時は`sequence.sequence`・CCA処理・断片生成（いずれも配列上の位置に依存するロジック）を完全にスキップする。
- Excel出力は既存の01_Index〜08_Visualizationの枠組み（Index/ハイパーリンク/切り詰め/書式）をそのまま流用しつつ、03/06シートをNucleosideTargetベースの内容（`03_Nucleoside_Targets`, `06_Nucleoside_Comparison`）に差し替える（04/05/07は共通）。
- **`ms1_peak_extraction.mz_min` に注意**（§24.5）: 既定値`500`はオリゴマー断片を想定した値。P1モードで検出対象となる遊離ヌクレオシドは概ね230〜370 Da（リン酸付加でも+80 Da程度）と大幅に軽いため、`mz_min`を100〜150程度まで下げないと検出漏れが起きる。`mz_min`が400以上のままP1モードを実行するとWARNINGが記録される。
- `modified_nucleoside_mass_mono`が無い`data/modifications.yaml`エントリに備えたフォールバック（§24.2）: target_baseの遊離ヌクレオシド質量 + `mass_shift_from_unmodified`で計算する（target_basesが単一塩基でない場合はスキップし警告する）。現状の118件は全て`modified_nucleoside_mass_mono`を直接持っており（フォールバックは未使用）、フォールバック計算値との差は最大でも約0.05mDa（4桁丸め相当、系統的なズレ無し）であることをテストで確認済み。

### 修飾仮説の理論質量チェック（`07b_Hypothesis_Check`）

「この位置にこの修飾があるはず」という仮説の理論質量が、生のMS1ピークに実在するかを確認する機能（仕様: `claude_code/hypothesis_mass_check_spec.md`）。**観測質量から修飾を自動同定することはしない**——ユーザーが指定した仮説について一致・不一致（Yes/No）を返すだけで、不一致の場合の候補生成もしない。P1モード・オリゴマーモードのどちらでも動作する。

- 仮説は2系統をマージして同じ照合ロジックにかける。
  - `sequence.trna_type`で選んだtRNAの`data/trna_library.yaml` `conserved_modifications`（Source = `trna_library_default`）。各エントリは`position`（配列自身の通し番号）・`modification`（単一ラベル）または`modification_candidates`（ラベルのリスト。各候補が独立した仮説に展開される）・`confidence`・`note`を持つ。`Confidence`列（`07b`）にそのまま出る。
    - `confirmed`: この生物種・株で直接検出済み。現状は **C+（agmatidine）を`tRNA-Ile2-CAT-1-1`のwobble位置に**（*M. acetivorans* C2A Δhpt株のtRNAでLC-MS確認、Gregorova et al. 2020, RNA Biol, [DOI:10.1080/15476286.2020.1853385](https://doi.org/10.1080/15476286.2020.1853385)）のみ。
    - `high_probability`: 一般則に基づく推定で、個々のtRNAでの直接確認ではない。全58件の **position 15 = G+（archaeosine）**（全長リコンストラクトでarchaeosine欠損が疑われるピークがあり、100%存在するとは言えない）と、下記の位置・塩基ルール由来の候補。
    - 位置・塩基ルール（標準tRNA番号から配列内位置を推定: 34 = `wobble_position`、37 = wobble+3、55 / 58 = discriminator（標準73）から逆算。塩基が一致する場合のみ追加）: U55 → Y（47件、`universal_U55_pseudouridine`。ΨはUと同質量なので質量だけではUと区別できない）／wobble U34 → cnm5U等10候補（15件、`ma_U34_main_target`）／A37 → t6A系4候補（35件、`ma_A37_t6A`）／A58 → m1A・m6A（47件、`archaea_A58_methylation`。汎用ラベル`methylation`は質量を持たずm1A/m6Aと同質量なので除外）／`tRNA-Phe-GAA-1-1/1-2`のG37 → imG-14・imG・imG2（ワイオシン経路。実データ`05_Mix`でも321.11 / 335.12 Daが複数ピーク検出）。imG-14をPhe以外に付けていないのはtRNA割り当てが未確認のため。
    - これらは`tools/add_conserved_modifications.py`で再生成できる（冪等、`--check`で差分確認）。手で直す場合はスクリプト側を直すこと（`tests/test_hypothesis_check.py`が生成結果とYAMLの一致を検査する）。
    - 過去の不整合の修正: `tRNA-Ala-TGC-1-1`の`wobble_position`が`37`と誤っていた（アンチコドン`UGC`は配列の34〜36番目）ため`34`に修正済み。これにより、このエントリにもU34系10候補・A37系4候補が付く。他の57件は`wobble_position`と配列上のアンチコドンが一致している。
  - `config.yaml`の`hypothesis_check.targets`（Source = `config_manual`、既定は空）。書式は`config.yaml`のコメント参照。
- ターゲット形式は3種: `label`（カタログのヌクレオシド質量そのまま）／`components`+`linkage`（ジヌクレオチド。5'/3'は区別しない）／`base`+`add_elements`（カタログ質量+元素の単同位体質量）。ジヌクレオチドは N1 + N2 + HPO3 − H2O（`phosphorothioate`はさらに +S −O）で、本体RNA_MassHunterの`p1_sap_dinucleotide_candidates.py`と同じ元素組成モデル。
- 各chargeを1〜`max_charge`と仮定して観測中性質量を逆算し、理論質量から`mass_tolerance_ppm`以内のピークを全て収集する。1仮説に複数ピークが一致すれば複数行に展開、一致無しは`Match_Found = No`の1行。Peak IDは04/05/06シートと共通。
- 理論質量は遊離ヌクレオシド基準（リン酸付加は考慮しない）。AP無し（5'-リン酸残存）のデータでは`add_elements`等で補正するか、P1+AP前提の仮説として解釈すること。
- 仮説が1件も無い、または`hypothesis_check.enabled: false`の場合はシートを作らない（既存の出力は変わらない）。未知のlabel・不正なtargetsは起動時に`ValueError`（`trna_library`側の不明labelはERROR警告でスキップ）。
- 実データ検証（`05_Mix.mzML`、P1+AP、positive）: m2,2G-PT-U（633.1254 Da）が PK72358〜PK72360（charge 1、Δppm −1.48〜−0.56）で一致。m2,2G-PT-C（632.1414）は最近傍が約26 ppm離れており、15 ppm許容では不一致。cnm5s2U+O+S / +S / ncm5s2Uは不一致（標準品ミックスとして妥当）。この結果は`tests/test_hypothesis_check.py`に統合テストとして組み込んである（`data/input/05_Mix.mzML`が無い環境ではskip）。

## 実データでの検証結果（Phase 8）

ユーザーの実mzML（tRNA-Gln2、profile-mode、2,167 MS1スペクトル、MS2スペクトルなし）で`main.py --config config.yaml`を実行し、以下を確認済み:

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
python main.py                      # リポジトリ直下の config.yaml を使う
python main.py --config config.yaml # 設定ファイルを明示する場合
python main.py --list-trna          # data/trna_library.yaml の tRNA 種類一覧
```

（エントリスクリプトは旧名 `simple_main.py`。`main.py` に改名した。）

`output/RNA_MassHunter_simple_report.xlsx` が生成される。`data/input/`（実mzML等の生データ）と`output/`（生成物）は`.gitignore`でリポジトリから除外されている。
