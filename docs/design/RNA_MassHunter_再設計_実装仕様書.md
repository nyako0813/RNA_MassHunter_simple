# RNA_MassHunter 再設計・実装仕様書（簡易版パイプライン）

作成日: 2026-09-02（2026-09-02にMS2ピーク対応を追記・改訂）
対象リポジトリ: `nyako0813/RNA_MassHunter`（main branch、調査時点コミット `6be9ced`）
入力資料: ChatGPTによる設計引き継ぎメモ（本仕様書に統合済み。原文はプロジェクトの `claude/rna-masshunter-design-memo.md` に保存済み）

> **改訂メモ**: 初版はMS1のみのフローだったが、ユーザーの指示により「MS2ピークを参考情報としてMass Comparisonに付加する」機能を追加した（自動修飾同定・スコアリングは行わず、あくまで観測MS2ピークが理論d/w/a/zイオンと一致した件数・一覧を示す裏付け情報として提示する）。変更箇所は §3.9, §4, §6, §7, §8.2, §8.5, §9.5, §14A, §15, §16.2, §18, §19, §22, §23, 未確定事項リストを参照。

この仕様書は、既存の大規模パイプライン（SCIEX監査群・MS2アノテーション・複合修飾探索など多数の高度機能を持つ）とは別に、**「配列→理論断片→観測質量→質量差→元素組成候補／既知修飾候補→Excel」という簡易的なフロー**を新設するための実装仕様書です。Claude Code（またはCodex）がこのままコーディングに入れる粒度まで具体化しています。

---

## 目次
1. プロジェクト目的
2. 現在のコード構成
3. 現在存在する関数（再利用対象の実シグネチャ）
4. 再利用する関数（確定版）
5. 修正するファイル
6. 新規作成するファイル
7. データフロー
8. クラス/関数仕様
9. config仕様
10. YAML仕様
11. 質量計算仕様
12. 元素組成候補アルゴリズム
13. 修飾候補検索
14. MS1ピーク処理
15. Mass Comparison
16. Excel仕様
17. Visualization仕様
18. エラー処理
19. 単体テスト
20. 統合テスト
21. 実mzMLによる動作確認
22. 実装順序
23. Claude Codeへの実装ルール

---

## 1. プロジェクト目的

RNA_MassHunter（簡易版）は、tRNA/RNAのLC-MS/MSデータを解析し、

- RNA配列から理論的なRNase断片を生成
- 各断片の理論質量を計算
- mzMLから観測されたMS1ピークを取得
- m/zと電荷から観測中性質量を計算
- 理論質量と観測質量を比較し、質量差（ΔDa / Δppm）を計算
- 質量差から**元素組成候補**（Formula Candidate）を列挙
- 既知RNA修飾候補（Modification Candidate）を列挙
- それらをExcelに整理する

ことを目的とする。

**基本思想**: 本ソフト自身が修飾を「確定同定」することを目的にしない。MSデータを正確に整理し、人間が修飾を判断できる情報を提供することを主目的とする。複雑な生物学的スコアリングや自動推論（既存リポジトリの `biological_context.py`, `evidence_ranking.py`, `modification_hypothesis_audit.py` 等が担う機能）は主解析フローから明確に除外する。

---

## 2. 現在のコード構成

既存リポジトリは非常に大規模（`rna_masshunter/` 配下に120以上のモジュール、`main.py` 単体で1494行、`excel_report.py` が2314行）で、SCIEX機器プロファイル監査（`sciex_*.py`、約30ファイル）、MS2アノテーション、修飾位置仮説監査、複合修飾探索など、tRNA修飾の自動同定を志向した高度な機能群を含む。

一方、今回新設する簡易版パイプラインが必要とする機能（配列処理・消化・質量計算・MS1ピーク抽出・既知修飾検索）は、この巨大な依存関係グラフから**独立した小規模モジュール群**として既に存在している（詳細は §3, §8 の依存関係調査を参照）。したがって、

- 既存の巨大な `main()` 関数やSCIEX監査系には一切手を入れない・依存を追加しない
- 簡易版は独立した新規オーケストレーション層（新規モジュール）として実装し、既存の小規模・独立モジュールのみをimportして再利用する

という方針を確定する（詳細は §5, §22 参照）。

---

## 3. 現在存在する関数（再利用対象の実シグネチャ）

以下は実コードを確認済みの正確なシグネチャ。設計メモの想定とほぼ一致しているが、差異がある箇所は「メモとの差異」として明記する。

### 3.1 `rna_masshunter/masses.py`（完全独立、依存は `warnings_manager` のみ）

```python
PROTON_MASS = 1.007276466812
DEFAULT_WATER_MASS = 18.010564684
DEFAULT_PHOSPHATE_MASS = 79.966331

MONOISOTOPIC_ATOMIC_MASSES = {
    "C": 12.0, "H": 1.00782503223, "N": 14.00307400443,
    "O": 15.99491461957, "P": 30.97376199842,
    "S": 31.9720711744, "Se": 79.9165218,
}

def load_base_masses(path, warnings=None) -> dict
def calculate_unmodified_rna_mass(sequence, base_masses, warnings=None, terminal_form="default") -> float | None
def neutral_mass_from_mz(mz, charge, polarity="negative", proton_mass=PROTON_MASS) -> float
def mz_from_neutral_mass(neutral_mass, charge, polarity="negative", proton_mass=PROTON_MASS) -> float
def elemental_delta_mass(elements: dict[str, int]) -> float
def fragment_ion_series_offsets(base_masses, warnings=None) -> dict[str, float]
```

`terminal_form` は `"default" | "dephosphorylated" | "residual_phosphate" | "cyclic_phosphate"` の4種（`dephosphorylated` と `default` は同じ0調整、命名上は別）。`water`/`phosphate` の実値は `data/base_masses.yaml` の `constants` セクションから読み込まれ、欠損時のみ `DEFAULT_*` 定数にフォールバックする。**メモとの差異なし。そのまま利用可。**

### 3.2 `rna_masshunter/elemental_composition.py`（完全独立、依存は `masses.py` の定数のみ）

```python
_ELEMENT_ORDER = ("C", "H", "N", "O", "P", "S", "Se")

@dataclass(frozen=True)
class ElementalComposition:
    def __init__(self, counts: Mapping[str, int] | None = None, *, allow_negative: bool = False)
    @classmethod
    def delta(cls, counts=None)              # allow_negative=True のショートカット
    def to_dict(self) -> dict[str, int]
    def __add__(self, other) -> ElementalComposition
    def __sub__(self, other) -> ElementalComposition
    @property
    def exact_mass -> float
    def canonical_string(self) -> str          # 例: "C2H2O1"、負数混在は "O-1S1" のような表記
    def is_close_mass(self, other, tolerance=1e-9) -> bool
    def __bool__(self) -> bool
```

**元素順序・対象元素（C,H,N,O,P,S,Se）ともに設計メモの想定と完全一致。** ただし「ΔDaから元素組成候補を列挙する」網羅探索ロジックはこのファイルに存在しない。新規実装が必要（§12）。

### 3.3 `rna_masshunter/enzymes.py` / `rna_masshunter/digestion.py`（`models.RunConfig` にのみ依存、SCIEX非依存）

```python
ENZYME_RULES = {
    "RNase_T1": {"cleaves_after": {"G"}, "specific": True},
    "RNase_A":  {"cleaves_after": {"C", "U"}, "specific": True},
    "RNase_T2": {"cleaves_after": {"A","C","G","U"}, "specific": False},
    "Nuclease_P1": {...}, "Benzonase": {...}, "U_specific_RNase": {...},
}
def get_enzyme_rule(enzyme_name: str) -> dict
def find_cleavage_sites(sequence: str, enzyme_name: str, warnings=None) -> list[int]

def digest_sequence(
    target_id: str, sequence: str, position_map: dict[int, int | None],
    config: RunConfig, base_masses: dict, warnings: list[dict] | None = None,
) -> list[Fragment]
def generate_terminal_forms(fragment: Fragment, config: RunConfig, base_masses: dict, warnings=None) -> list[Fragment]
```

`config.digestion`（enabled/enzyme/digestion_mode/missed_cleavages/min_length/max_length/include_terminal_forms/allow_partial_digestion/allow_nonspecific_cleavage）と `config.alkaline_phosphatase`（enabled/assume_complete/allow_residual_phosphate/allow_cyclic_phosphate）は**既に独立したセクション**であり、メモが要求する「RNase処理とAP処理を混同しない」設計に既に合致している。**そのまま再利用可。**

### 3.4 `rna_masshunter/cca_processing.py` / `cca_tail_state.py`

```python
class CCAMaturationState(str, Enum):  # NONE / C / CC / CCA

@dataclass(frozen=True)
class CCAProcessingResult:
    original_sequence: str
    processed_sequence: str
    original_tail_state: CCAMaturationState
    added_suffix: str
    added_nucleotide_count: int
    cca_processing_applied: bool

def process_cca_tail(
    sequence: str,
    registered_sequence_cca_mode: RegisteredSequenceCCAMode | str,
    *, enabled: bool = True,
) -> CCAProcessingResult
```

**メモとの差異（重要）**: `process_cca_tail()` は実装済みだが、`main.py` から一度も呼び出されていない（デッドコード）。新設計で初めて配線する必要がある。`RegisteredSequenceCCAMode`（`EXCLUDES_CCA` / `INCLUDES_COMPLETE_CCA` / `UNKNOWN`）は `cca_tail_state.py` にあるEnumだが、`cca_processing.py` はこのEnumのみを利用し、`cca_tail_state.py` の重量級クラス（`CCATailVariant` 等、複数CCA仮説の網羅生成用）には依存しない。**新設計では `cca_processing.py` と `RegisteredSequenceCCAMode` Enumのみを取り込み、`cca_tail_state.py` の他の重い機能は使わない。**

### 3.5 `rna_masshunter/mzml_reader.py` / `peak_picking.py`

```python
def iter_spectra(mzml_path) -> Iterator[dict]
def extract_ms1_spectra(mzml_path) -> list[dict]
def extract_ms2_spectra(mzml_path) -> list[dict]

def extract_ms1_peaks(mzml_path, reconstruction_config: dict, warnings=None) -> list[Peak]
```

`extract_ms1_peaks()` の第2引数名は `reconstruction_config` だが、中身は `rt_min/rt_max/mz_min/mz_max/intensity_threshold` の汎用フィルタであり、既存 `main.py` でも `config.reconstruction` をそのまま渡している（呼び出し例は §8.3 参照）。**メモの想定通り。そのまま再利用可。**

`Peak` データクラス（`models.py`）:
```python
@dataclass
class Peak:
    mz: float
    intensity: float
    rt: float | None = None
    scan_id: str | None = None
    ms_level: int = 1
    tier: str | None = None
```

### 3.6 `rna_masshunter/modifications.py`

```python
def load_modifications(path, warnings=None) -> list[Modification]
def validate_modifications(modifications, warnings=None) -> None
def find_modifications_by_mass_shift(
    modifications: list[Modification], mass_shift: float, tolerance_da: float = 0.01,
) -> list[Modification]
```

**メモとの差異**: 許容差指定は絶対値Da固定（ppm指定は無い）。線形スキャンで絶対誤差フィルタのみ行い、**結果は一致度順にソートされない**（呼び出し側でソートが必要）。`Modification` データクラスのフィールドは修飾名(`symbol`/`raw["name"]`)・mass_shift_from_unmodified・category・target_bases・detectability・curation・sources・candidate_policy・chemical_group・near_isobaric_group すべて実在し、メモの想定と完全一致。

### 3.7 `rna_masshunter/mass_shift_ms1_search.py`（新設計の探索エンジンとして最重要）

```python
@dataclass
class SortedPeakIndex:
    peaks: list[Peak]
    mzs: list[float]

def build_sorted_peak_index(peaks: list[Peak]) -> SortedPeakIndex
def find_peaks_near_mz(index: SortedPeakIndex, theoretical_mz: float, tolerance_ppm: float) -> list[PeakMatch]
```

`bisect` による二分探索でtolerance窓内のピークのみ抽出する軽量な探索基盤。`ms1_mapping.ppm_error` / `_confidence` を内部で再利用している。**新設計の `mass_comparison.py` の探索エンジンとして、既存 `map_fragments_to_ms1_peaks`（全ピーク線形走査）より高速かつ簡潔なため、こちらを優先的に再利用する。**

### 3.8 `rna_masshunter/ms1_mapping.py`

```python
def theoretical_mz_from_mass(neutral_mass, charge, polarity, proton_mass=PROTON_MASS) -> float
def ppm_error(observed, theoretical) -> float
def _confidence(error_ppm, tolerance_ppm, peak_tier) -> str
def map_fragments_to_ms1_peaks(fragments, peaks, config, warnings=None, audit_context=None) -> list[FragmentMS1Match]
```

`ppm_error` と `_confidence`（High/Medium/Low判定）はそのまま流用する。`map_fragments_to_ms1_peaks` 自体は新設計では直接使わず（全ピーク線形走査のため §3.7 の二分探索版を優先）、代わりにロジックの一部（ppm誤差計算・信頼度判定）のみ再利用する。

### 3.9 `rna_masshunter/ms2_annotation.py`（MS2参考情報向け、2026-09-02追記）

**背景**: ユーザーの指示により、Mass ComparisonにMS2ピークを「参考情報」として付加する（自動同定・スコアリングは行わない）。既存の `ms2_annotation.py`（53KB、約1170行）は本格的なMS2アノテーション機能（前駆体候補探索・修飾レスキュー・エビデンスレベル判定・配列カバレッジ計算等）を持つが、そのすべてを取り込むのはメモの基本思想（自動同定を目的にしない）に反する。**以下の3関数のみを再利用し、それ以外（`find_parent_candidates`, `match_ms2_spectra`, `build_fragment_evidence`, `_evidence_level`, `_sequence_coverage` 等の前駆体レスキュー・エビデンス判定ロジック）は使わない。**

```python
def extract_ms2_spectra(
    mzml_path: str, ms2_config: dict, warnings: list[dict] | None = None,
    zero_intensity_context: dict | None = None,
) -> list[MS2SpectrumInfo]
```
mzMLからMS2スペクトルのみ抽出し、強度閾値・相対強度閾値・最大ピーク数でフィルタする。

**重要な訂正（2026-09-02、リポジトリ精査により判明）**: `ms2_annotation.py` はモジュール冒頭で `base_loss_masses.py`, `base_loss_ions.py`, `modified_precursor.py`, `modified_fragment_ions.py`, `ms2_unmatched_audit.py`, `ms2_zero_intensity_audit.py` を import しており、`extract_ms2_spectra`/`generate_theoretical_ms2_ions` の2関数だけが目的でも `from rna_masshunter.ms2_annotation import ...` を書くとこれら重量級モジュール一式がPythonのimport実行時に連鎖的に読み込まれる。さらに `extract_ms2_spectra` の内部は `ms2_zero_intensity_audit.capture_source_spectrum/record_parsed_spectrum/record_parser_error`（ゼロ強度監査用の記録処理）を**無条件に**呼び出しており、単純コピーでもこの依存を切れない。したがって**「そのまま再利用」はできない**。正しい対応は、`extract_ms2_spectra` と `generate_theoretical_ms2_ions` の関数本体をコピーしたうえで、(a) `capture_source_spectrum`/`record_parsed_spectrum`/`record_parser_error` の呼び出し行を削除する（ゼロ強度監査は今回のスコープ外のため実害なし）、(b) 両関数が実際に必要とする小さな非公開ヘルパー（`_scan_window`, `_precursor_info`, `_ion_charges`, `_ion_series`, `_as_bool`, `_optional_positive_int`, `_safe_float`）も一緒にコピーする、という**手を入れたうえでの再利用**とする。具体的な移植手順は実装引き継ぎ指示書（`claude/RNA_MassHunter_Claude_Code_実装引き継ぎ指示書.md`）に記載する。

```python
def generate_theoretical_ms2_ions(
    theoretical_fragments: list[Fragment], config: Any, base_masses: dict, warnings=None,
) -> list[TheoreticalMS2Ion]
```
各理論断片(`Fragment`)の全カット位置について、未修飾の d/w/a/z 系列イオン（`config.ms2_annotation.ion_series` で指定、既定 `["d","w","a","z"]`）の理論m/zを生成する。**コード内コメントに明記の通り「unmodified d/w/a/z series ion」であり、README記載の「c/yイオン」という表現は実装と不一致（実際に使われているのはd/w/a/z系列）。仕様書・ドキュメントではd/w/a/zで統一する。** この関数自体は `ms2_zero_intensity_audit` 等への依存を持たないが、`ms2_annotation.py` からのimportはモジュール冒頭の重量級import連鎖を引き込むため、上記 `extract_ms2_spectra` と同様に**関数本体をコピーして再利用**する（必要な非公開ヘルパー `_ion_charges`, `_ion_series`, `_as_bool` も同様にコピー）。

```python
def _best_ion_match(observed_mz: float, ions: list[TheoreticalMS2Ion], tolerance_ppm: float)
```
観測m/zに対しtolerance内で最も近い理論イオンを返す軽量関数（線形走査、`ppm_error`利用）。**モジュール非公開（アンダースコア始まり）のため、そのままimportするのではなく、新規 `ms2_support.py` 内に同等の5行程度のロジックとしてコピー実装することを推奨**（`ms2_annotation.py` 側の将来のリファクタで壊れるリスクを避けるため）。

`MS2SpectrumInfo` / `TheoreticalMS2Ion`（`models.py`）の主要フィールド:
```python
@dataclass
class MS2SpectrumInfo:
    spectrum_id: str; scan_index: int; rt: float | None
    precursor_mz: float | None; precursor_charge: int | None; precursor_intensity: float | None
    num_peaks: int; peaks: list[tuple[float, float]]  # (mz, intensity)のフィルタ後リスト
    # 他: base_peak_mz, base_peak_intensity, total_ion_current 等

@dataclass
class TheoreticalMS2Ion:
    ion_id: str; parent_fragment_id: str; parent_sequence: str
    ion_type: str          # "d" | "w" | "a" | "z"
    ion_sequence: str; ion_start: int; ion_end: int; charge: int
    theoretical_mz: float; theoretical_mass: float
```

`config.ms2_annotation`（`config.py` に既存、**新設計で新規キー追加は不要**）:
```yaml
ms2_annotation:
  enabled: true
  mz_tolerance_ppm: 20
  min_peak_intensity: 10
  min_relative_intensity_percent: 1.0
  max_peaks_per_spectrum: 500
  precursor_match_tolerance_ppm: 20
  use_theoretical_fragments: true
  ion_series: ["d", "w", "a", "z"]
  include_neutral_loss: false
  include_base_loss: false
  min_ion_length: 1
```
このセクションは既に必要なフィールドを全て備えているため、**新設計は既存の `ms2_annotation` セクションをそのまま読み、config.pyへの変更は不要**（§9.5参照）。

---

## 4. 再利用する関数（確定版）

| モジュール | 再利用する関数/クラス | 用途 |
|---|---|---|
| `masses.py` | `calculate_unmodified_rna_mass`, `neutral_mass_from_mz`, `mz_from_neutral_mass`, `elemental_delta_mass`, `load_base_masses`, `PROTON_MASS` 等定数 | 理論質量・観測中性質量の計算 |
| `elemental_composition.py` | `ElementalComposition`（`delta()`, `exact_mass`, `canonical_string()`） | Formula Candidateの質量計算・表記生成 |
| `enzymes.py` | `get_enzyme_rule`, `find_cleavage_sites`, `ENZYME_RULES` | RNase A/T1 定義（`digestion.py` 経由で間接利用） |
| `digestion.py` | `digest_sequence`, `generate_terminal_forms` | 理論断片生成（missed cleavage・AP対応込み） |
| `cca_processing.py` | `process_cca_tail`, `CCAMaturationState`, `CCAProcessingResult` | CCA末端処理 |
| `cca_tail_state.py` | `RegisteredSequenceCCAMode`（Enumのみ） | `process_cca_tail` の引数型 |
| `mzml_reader.py` | （直接は使わず `peak_picking.py` 経由） | — |
| `peak_picking.py` | `extract_ms1_peaks` | mzMLからのMS1ピーク抽出 |
| `modifications.py` | `load_modifications`, `validate_modifications`, `find_modifications_by_mass_shift` | Modification Candidate列挙 |
| `mass_shift_ms1_search.py` | `build_sorted_peak_index`, `find_peaks_near_mz`, `SortedPeakIndex` | Fragment×charge×Peakの高速マッチング |
| `ms1_mapping.py` | `ppm_error`, `_confidence` | Δppm計算・信頼度判定 |
| `models.py` | `Fragment`, `Peak`, `Modification`, `MS2SpectrumInfo`, `TheoreticalMS2Ion` | データ構造の基盤（流用、拡張はしない） |
| `ms2_annotation.py` | `extract_ms2_spectra`, `generate_theoretical_ms2_ions`（`_best_ion_match` はコピー実装、§3.9参照） | MS2参考情報（2026-09-02追記） |
| `warnings_manager.py` | `add_warning` | エラー/警告記録の統一 |
| `config.py` | `load_config`, `validate_config`, `resolve_paths`, `DEFAULT_CONFIG`, `_merge_defaults` | 設定ロード（新セクション追加のうえ流用） |
| `excel_report.py`（ヘルパーのみ） | `_add_index_and_backlinks`, `_autosize_and_freeze`, `_sheet_link`, `_coerce_to_frame`, `_excel_safe_cell`, `_flatten_dict`, `_truncate_frame_if_needed` | Index/ハイパーリンク/書式まわりの共通処理 |

**再利用しないもの（設計思想が異なるため）**:
- `unknown_modification.py`（`generate_unknown_modification_candidates`）: 固定6種のdelta限定探索であり、新設計の「7原子までの網羅探索」とは方式が根本的に異なる。骨格（fragment×charge走査）のみ参考にし、候補生成ロジックは新規実装する。
- `modification_search.py` の `search_known_modifications_by_mass_shift`: ロジックはModification Candidateに極めて近いが、出力列（`KnownModificationCandidate`）がSCIEX/shadow audit向けの余剰フィールドを多数含み、新設計の06シート構成に対して過剰。**探索エンジン部分（`build_sorted_peak_index`, `find_peaks_near_mz`）のみ再利用し、候補生成・出力整形は `mass_comparison.py` で新規に行う。**
- `excel_report.py` の `write_excel_report()` 本体（670行の巨大関数）: 既存の出力シート構成（Intact/SCIEX/shadow audit系）と新設計の01〜08シート構成は目的が異なるため、本体ロジックは再利用せず、新規関数として実装する（§5, §16参照）。
- `main.py` の `main()` 本体: 850行超で分岐が複雑に絡み合っており部分抽出は困難。新設計は薄い新規オーケストレーション層として実装する。
- `ms2_annotation.py` の `match_ms2_spectra` / `find_parent_candidates` / `build_fragment_evidence`、および `evidence_ranking.py`, `biological_context.py`, `modification_hypothesis_audit.py`, `modification_hypothesis_schema.py`: 前駆体候補レスキュー・エビデンスレベル判定・配列カバレッジ・生物学的コンテキスト評価など、いずれも「自動同定・スコアリング」寄りの機能であり、メモの基本思想（人間が判断するための情報整理に徹する）に反するため一切使わない（2026-09-02追記）。

---

## 5. 修正するファイル

| ファイル | 変更内容 | 理由 |
|---|---|---|
| `rna_masshunter/config.py` | `DEFAULT_CONFIG` に `formula_candidate` セクションを新規追加。`cca_processing` セクションに `registered_sequence_cca_mode` フィールドを追加。`validate_config()` に上記の簡易バリデーションを追加 | Formula CandidateとCCA処理を設定可能にするため（詳細は§9） |
| `rna_masshunter/models.py` | `RunConfig` に `formula_candidate: dict = field(default_factory=dict)` を追加。新規dataclass `FormulaCandidate`, `MassComparisonRow` を追加 | 新設計の中心データ構造を既存パターン（dataclass集約）に合わせて配置するため |
| `rna_masshunter/excel_report.py` | 既存関数は変更せず、新規関数 `write_simple_mass_hunter_report(...)` を末尾に追加。共通ヘルパー（`_add_index_and_backlinks` 等）はそのままimportして再利用 | 01〜08シート構成の新規Excel出力を既存の巨大関数に混ぜずに追加するため |
| `config.yaml` | `formula_candidate` セクションのサンプル値を追記。ついでに既存の `alkaline_phosphatase:` セクションの重複定義（62行目・67行目）を1つに統合（バグ修正、副次対応） | 新機能の設定サンプルを提供しつつ、調査で発見した軽微なバグを解消するため |
| `main.py` | **変更しない**（新規エントリポイントを別ファイルとして追加する。理由は§22参照） | 既存の巨大な `main()` への影響を避けるため |

---

## 6. 新規作成するファイル

| ファイル | 役割 |
|---|---|
| `rna_masshunter/formula_candidate.py` | ΔDaから元素組成候補（Formula Candidate）を列挙・ランキング |
| `rna_masshunter/mass_comparison.py` | Fragment×Peakマッチング、ΔDa/Δppm計算、Formula/Modification両候補の統合、Mass Comparison行の生成 |
| `rna_masshunter/observed_mass.py` | MS1ピークからObserved_Mass/Mass_Intensityシート用の行データ（Peak IDの安定採番含む）を生成する薄いヘルパー |
| `rna_masshunter/ms2_support.py`（2026-09-02追記） | MS2スペクトルを前駆体(fragment×charge)に紐付け、理論d/w/a/zイオンとの一致件数・一覧を返す薄いヘルパー（§8.5） |
| `rna_masshunter/simple_pipeline.py` | 新設計フローのオーケストレーション本体（config読込→CCA→消化→AP→MS1ピーク抽出→MS2スペクトル抽出→Mass Comparison→Excel出力） |
| `simple_main.py`（リポジトリ直下） | CLIエントリポイント（`main.py` とは独立。`argparse` で `--config` を受け取り `simple_pipeline.run()` を呼ぶだけの薄いラッパー） |
| `test_formula_candidate.py` | `formula_candidate.py` の単体テスト |
| `test_mass_comparison.py` | `mass_comparison.py` の単体テスト |
| `test_ms2_support.py`（2026-09-02追記） | `ms2_support.py` の単体テスト |
| `test_simple_pipeline.py` | 統合テスト（合成データによるエンドツーエンド確認） |

「既存機能で十分なら不要なモジュールは増やさない」という原則に従い、`observed_mass.py` はロジックを持たせず薄いヘルパーに留める（本格的なロジックが増えなければ `mass_comparison.py` に統合してもよい。実装時に判断）。

---

## 7. データフロー

```text
config.yaml
   │
   ▼
config.py: load_config / validate_config / resolve_paths   [既存, 再利用]
   │
   ▼
sequence (config.sequence.sequence)
   │
   ▼
cca_processing.process_cca_tail()                            [既存だが未配線→新規配線]
   │  processed_sequence
   ▼
digestion.digest_sequence()                                  [既存, 再利用]
   │  (RNase A/T1・missed cleavage・AP end formをconfig経由で内包)
   ▼
theoretical_fragments: list[Fragment]  (各Fragment.unmodified_massを保持)
   │
   ├───────────────────────────────────────┬───────────────────────────────┐
   ▼                                       │                               │
peak_picking.extract_ms1_peaks(mzml)       │   ms2_annotation.extract_ms2_spectra(mzml)   [既存, 再利用, 2026-09-02追記]
   │  observed_peaks: list[Peak]           │       │  ms2_spectra: list[MS2SpectrumInfo]
   ▼                                       │       ▼
mass_comparison.build_mass_comparison_rows(theoretical_fragments, observed_peaks, modifications, config, ms2_spectra=ms2_spectra)
   │   内部で:
   │     - mass_shift_ms1_search.build_sorted_peak_index / find_peaks_near_mz  [既存, 再利用]
   │     - masses.neutral_mass_from_mz                                        [既存, 再利用]
   │     - ΔDa = observed_mass - fragment.unmodified_mass
   │     - Δppm = ms1_mapping.ppm_error(observed, theoretical)                [既存, 再利用]
   │     - formula_candidate.enumerate_formula_candidates(ΔDa, tolerance)     [新規]
   │     - modifications.find_modifications_by_mass_shift(mods, ΔDa, tol)     [既存, 再利用]
   │     - ms2_support.find_ms2_support(fragment_id, charge, observed_mz, ms2_spectra, ion_index, config)  [新規, 2026-09-02追記]
   ▼
mass_comparison_rows: list[MassComparisonRow]
   │
   ▼
excel_report.write_simple_mass_hunter_report(...)             [新規、既存ヘルパー再利用]
   │
   ▼
output/*.xlsx  (01_Index 〜 08_Visualization)
```

memoの概念図（§7 補足として再掲、実装が完全対応することを確認済み）:

```text
Observed Mass
      ↓
Theoretical Massとの差
      ↓
     ΔDa
      │
      ├───────────────┐
      ↓               ↓
Formula Candidate   Modification Candidate
      ↓               ↓
Recommended         Known
Formula             Modification
      ↓               ↓
Formula Candidates  Modification Candidates
```

---

## 8. クラス/関数仕様

### 8.1 `rna_masshunter/formula_candidate.py`

```python
from dataclasses import dataclass
from rna_masshunter.elemental_composition import ElementalComposition, _ELEMENT_ORDER
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES


@dataclass(frozen=True)
class FormulaCandidateResult:
    formula: str            # ElementalComposition.canonical_string()
    exact_mass: float        # ElementalComposition.exact_mass
    mass_error_da: float     # exact_mass - mass_difference
    mass_error_ppm: float | None   # 参考値。理論質量が別途分かる場合のみ計算、無ければNone
    atom_count: int          # sum(abs(v) for v in composition.to_dict().values())
    elements: dict[str, int]  # composition.to_dict()


def calculate_delta_mass(observed_mass: float, theoretical_mass: float) -> float:
    """ΔDa = observed - theoretical。単なる引き算だが呼び出し側の意図を明確にするための薄いラッパー。"""
    return observed_mass - theoretical_mass


def enumerate_formula_candidates(
    mass_difference: float,
    tolerance_da: float,
    *,
    max_total_atoms: int = 7,
    element_limits: dict[str, int] | None = None,
    elements: tuple[str, ...] = _ELEMENT_ORDER,
) -> list[FormulaCandidateResult]:
    """
    質量差 mass_difference (Da) を tolerance_da 以内で説明できる元素組成を、
    C/H/N/O/P/S/Se の符号付き個数の組み合わせとして max_total_atoms 以下の
    範囲で全探索する。

    element_limits: {"C": 3, "H": 7, ...} のように元素ごとの絶対値上限を指定可能。
    未指定の要素は max_total_atoms が実質的な上限になる。

    戻り値は mass_error_da の絶対値昇順（合致度順）。ties は formula文字列の
    辞書順で安定ソートする。
    """


def rank_formula_candidates(
    candidates: list[FormulaCandidateResult],
) -> list[FormulaCandidateResult]:
    """abs(mass_error_da) 昇順にソートするだけの明示的なヘルパー
    （enumerate_formula_candidates が既にソート済みでも、外部から個別に呼べるようにする）。"""


def get_recommended_formula(candidates: list[FormulaCandidateResult]) -> FormulaCandidateResult | None:
    """ランキング済みリストの先頭を返す。空リストならNone。"""
    return candidates[0] if candidates else None
```

**探索アルゴリズムの詳細は §12 を参照。**

### 8.2 `rna_masshunter/mass_comparison.py`

```python
from dataclasses import dataclass, field


@dataclass
class MassComparisonRow:
    peak_id: str
    fragment_id: str
    sequence: str
    charge: int
    intensity: float
    observed_mz: float               # 2026-09-02追記: MS2前駆体マッチングにも使うため明示的に保持
    observed_mass: float
    theoretical_mass: float
    delta_da: float
    delta_ppm: float
    recommended_formula: str | None
    formula_candidates: str          # "C2H2O1 [+0.00003 Da]; C3H6 [+0.00121 Da]; ..." 形式
    known_modification: str | None
    modification_candidates: str     # "m1A [Δ=+0.0002 Da]; ..." 形式
    ms2_spectrum_id: str | None = None        # 2026-09-02追記
    ms2_matched_ion_count: int = 0            # 2026-09-02追記
    ms2_matched_ions: str = ""                # 2026-09-02追記: "d3 [-0.5 ppm]; w5 [+1.2 ppm]" 形式
    rt: float | None = None
    scan_id: str | None = None
    warnings: list[str] = field(default_factory=list)


def build_mass_comparison_rows(
    fragments: list["Fragment"],
    peaks: list["Peak"],
    modifications: list["Modification"],
    config: "RunConfig",
    warnings: list[dict] | None = None,
    ms2_spectra: list["MS2SpectrumInfo"] | None = None,   # 2026-09-02追記
    ms2_ion_index: dict[str, list["TheoreticalMS2Ion"]] | None = None,  # 2026-09-02追記
) -> list[MassComparisonRow]:
    """
    fragments × charge(config.fragment_mapping.min_charge..max_charge) × peaks を
    mass_shift_ms1_search の二分探索で総当たりし、tolerance内のマッチごとに
    MassComparisonRow を1行生成する。

    手順:
      1. build_sorted_peak_index(peaks) でピーク集合を1回だけソート
      2. 各fragment・各chargeについて理論m/z (ms1_mapping.theoretical_mz_from_mass) を計算
      3. find_peaks_near_mz(index, theoretical_mz, config.mass_comparison.mz_tolerance_ppm) でマッチ抽出
      4. マッチごとに observed_mass = masses.neutral_mass_from_mz(peak.mz, charge, polarity)
      5. delta_da = observed_mass - fragment.unmodified_mass
      6. delta_ppm = ms1_mapping.ppm_error(observed_mass, fragment.unmodified_mass)
      7. formula_candidates = formula_candidate.enumerate_formula_candidates(delta_da, tolerance_da, ...)
      8. modification_candidates = sorted(
             modifications.find_modifications_by_mass_shift(modifications, delta_da, tolerance_da),
             key=lambda m: abs(m.mass_shift_from_unmodified - delta_da),
         )
      9. peak_id は observed_mass.assign_peak_ids() で採番済みのIDをpeakオブジェクトから引く
                     （Peak自体にIDフィールドが無いため、peaksをリストとして受け取った時点の
                       インデックス順で採番し、Observed_Mass/Mass_Intensityシートと共有する）
      10. (2026-09-02追記) config.ms2_annotation.enabled かつ ms2_spectra が渡されていれば、
          ms2_support.find_ms2_support(fragment.fragment_id, charge, peak.mz, ms2_spectra,
          ms2_ion_index, config) を呼び、結果があれば ms2_spectrum_id / ms2_matched_ion_count /
          ms2_matched_ions を埋める。ms2_spectra が None、または該当スペクトルが無い場合は
          既定値（None / 0 / ""）のまま、エラーにはしない。
    """
```

**ΔDaが±0の場合の扱い**: `formula_candidates` に空文字列（元素差分なし＝一致）を含めてよいかは実装時にconfigで選択可能にする（`mass_comparison.include_zero_delta_as_match: bool`）。既定値は `true`（理論値と完全一致した観測は「差分なし」として明示的に1行残す）。

### 8.3 `rna_masshunter/observed_mass.py`

```python
@dataclass
class ObservedMassRow:
    peak_id: str
    mz: float
    charge: int | None
    observed_mass: float | None
    intensity: float
    rt: float | None
    scan_id: str | None


def assign_peak_ids(peaks: list["Peak"]) -> dict[int, str]:
    """peaksのリストインデックスに基づき 'PK0001' 形式の安定IDを採番して返す
    （{index: peak_id} の辞書）。mass_comparison.py / excel_report.py 双方から
    同じ採番結果を参照させることで Peak ID による追跡性（原則8）を担保する。"""


def build_observed_mass_rows(peaks: list["Peak"], charge_hint: dict[int, int] | None = None) -> list[ObservedMassRow]:
    """Observed_Mass / Mass_Intensity シート用の行データを生成する薄いヘルパー。
    charge_hint が無い場合、chargeはNone（後段のMass Comparisonで判明した
    charge値のみ埋める2パス処理も検討可、詳細は実装時に決定）。"""
```

### 8.4 `rna_masshunter/simple_pipeline.py`

```python
def run(config_path: str | Path) -> dict:
    """
    1. config = config.load_config(config_path); config = config.resolve_paths(config, project_root)
    2. modifications = modifications.load_modifications(...)
    3. sequence = config.sequence["sequence"].upper().replace("T", "U")
    4. cca_result = cca_processing.process_cca_tail(
           sequence,
           config.cca_processing.get("registered_sequence_cca_mode", "EXCLUDES_CCA"),
           enabled=config.cca_processing.get("enabled", True),
       )
    5. base_masses = masses.load_base_masses(...)
    6. theoretical_mass = masses.calculate_unmodified_rna_mass(cca_result.processed_sequence, base_masses)
    7. fragments = digestion.digest_sequence(target_id=..., sequence=cca_result.processed_sequence, ...)
    8. peaks = peak_picking.extract_ms1_peaks(mzml_path, config.reconstruction)  # 既存の流用パターンを踏襲
    9. (2026-09-02追記) ms2_spectra = ms2_annotation.extract_ms2_spectra(mzml_path, config.ms2_annotation, warnings) if config.ms2_annotation.get("enabled") else []
       ms2_ion_index = ms2_support.build_ms2_ion_index(fragments, config, base_masses, warnings) if ms2_spectra else {}
    10. rows = mass_comparison.build_mass_comparison_rows(fragments, peaks, modifications, config, ms2_spectra=ms2_spectra, ms2_ion_index=ms2_ion_index)
    11. excel_report.write_simple_mass_hunter_report(output_path, config, fragments, peaks, rows, modifications)
    12. 戻り値として集計結果(dict)を返す（テスト・CLI双方から使えるように）
    """
```

`simple_main.py`（リポジトリ直下）はこの `run()` を呼ぶだけの数行のCLIラッパーとする（`argparse` で `--config` のみ受け取る、既存 `main.py` の `parse_args`/`resolve_config_path` パターンを踏襲してよい）。

### 8.5 `rna_masshunter/ms2_support.py`（2026-09-02追記）

```python
from dataclasses import dataclass
from rna_masshunter.ms2_annotation import extract_ms2_spectra, generate_theoretical_ms2_ions
from rna_masshunter.ms1_mapping import ppm_error


@dataclass(frozen=True)
class MS2SupportResult:
    spectrum_id: str
    matched_ion_count: int
    matched_ions: str   # "d3 [-0.5 ppm]; w5 [+1.2 ppm]; a2 [+0.3 ppm]" 形式


def build_ms2_ion_index(
    theoretical_fragments: list["Fragment"], config: "RunConfig", base_masses: dict, warnings=None,
) -> dict[str, list["TheoreticalMS2Ion"]]:
    """generate_theoretical_ms2_ions() の結果を parent_fragment_id ごとにグルーピングして返す。
    Mass Comparisonの行数分（fragment×charge×peakマッチ数）呼び出されても再計算しないよう、
    simple_pipeline.py 側で1回だけ計算して使い回す。"""
    ions = generate_theoretical_ms2_ions(theoretical_fragments, config, base_masses, warnings)
    index: dict[str, list] = {}
    for ion in ions:
        index.setdefault(ion.parent_fragment_id, []).append(ion)
    return index


def find_ms2_support(
    fragment_id: str,
    charge: int,
    observed_precursor_mz: float,
    ms2_spectra: list["MS2SpectrumInfo"],
    ms2_ion_index: dict[str, list["TheoreticalMS2Ion"]],
    config: "RunConfig",
) -> MS2SupportResult | None:
    """
    1. config.ms2_annotation.precursor_match_tolerance_ppm 以内で、precursor_charge が一致する
       （またはNoneで未記録の）スペクトルを ms2_spectra から抽出する
       （ppm_error(observed_precursor_mz, spectrum.precursor_mz) で判定）。
    2. 複数一致する場合は base_peak_intensity が最大のスペクトルを採用する
       （タイブレークの単純化、原則2「自動同定が目的ではない」を踏まえ、複雑な優先順位付けはしない）。
    3. 一致するスペクトルが無ければ None を返す（MS2証拠なし。ΔDa比較自体はMS1のみで完結するため
       Mass Comparison行自体は生成され、MS2列が空欄になるだけ）。
    4. ms2_ion_index.get(fragment_id, []) の理論d/w/a/zイオンに対し、採用したスペクトルの
       各観測ピークを config.ms2_annotation.mz_tolerance_ppm 以内でマッチさせる
       （§3.9 で述べた _best_ion_match と同等のロジックをここに実装、線形走査で十分な規模）。
    5. マッチしたイオンをppm誤差昇順で並べ、MS2SupportResult(spectrum_id, matched_ion_count, matched_ions) を返す。
    """
```

**スコープ外にする点（原則2の徹底、2026-09-02追記）**:
- 修飾済みイオン（ΔDaを反映したd/w/a/zイオン）の探索はしない。あくまで未修飾理論イオンとの一致件数を示すだけであり、修飾位置の特定（localization）は行わない。
- `config.ms2_annotation.include_neutral_loss` / `include_base_loss` は既定 `false` のまま踏襲し、有効化しない（スコープ拡大を避けるため）。
- 前駆体候補の「レスキュー」（`find_modified_parent_candidates` 相当、修飾質量を仮定した前駆体再探索）は行わない。前駆体との対応づけは、Mass Comparisonの行が既に確定させた `(fragment_id, charge, observed_mz)` をそのまま使う。

---

## 9. config仕様

### 9.1 新規セクション: `formula_candidate`

```yaml
formula_candidate:
  enabled: true
  max_total_atoms: 7
  elements: ["C", "H", "N", "O", "P", "S", "Se"]
  element_limits:
    C: 3
    H: 7
    N: 3
    O: 4
    P: 1
    S: 1
    Se: 1
  mass_tolerance_da: 0.01
  max_candidates_per_match: 10
  include_zero_delta_as_match: true
```

`element_limits` の初期値は暫定値であり、実装時に `modifications.yaml` の既知修飾の元素組成分布（例えば最大の修飾がどの程度の原子数を要するか）を確認したうえで最終決定する（メモ§16の指示どおり）。`DEFAULT_CONFIG` への追加は既存の `modification_search` / `unknown_modification_search` セクションと同じ形式（`enabled`, `mz_tolerance_ppm` 等の慣例）に揃える。

### 9.2 新規セクション: `mass_comparison`

```yaml
mass_comparison:
  enabled: true
  mz_tolerance_ppm: 10          # 既存 fragment_mapping.mz_tolerance_ppm と同じ既定値を踏襲
  min_charge: 1                 # 既存 fragment_mapping.min_charge を踏襲（重複を避けるため参照のみでも可、実装時に決定）
  max_charge: 8
  polarity: "auto"
```

`fragment_mapping` セクションと値が重複するため、実装時に「`mass_comparison` は独自セクションを持たず `fragment_mapping` の該当フィールドを流用する」か「明示的に別セクションとして持つか」を選択する（**推奨: 既存 `fragment_mapping` セクションをそのまま流用し、`mass_comparison` セクションは `enabled` と `mz_tolerance_ppm` のみの最小構成にする**。重複設定によるユーザーの混乱を避けるため）。

### 9.3 `cca_processing` セクションの拡張

現状:
```yaml
cca_processing:
  enabled: true
```

拡張後:
```yaml
cca_processing:
  enabled: true
  registered_sequence_cca_mode: "EXCLUDES_CCA"   # "EXCLUDES_CCA" | "INCLUDES_COMPLETE_CCA" | "UNKNOWN"
```

これはユーザー自身が過去に確定した設計判断（配列末尾の機械的自動判定は採用せず、config側で明示指定する）と整合する。

### 9.4 出力パス

新規セクションは設けず、既存の `project.output_dir` をそのまま利用する（例: `output/RNA_MassHunter_simple_report.xlsx`）。ファイル名は `reporting` セクションに `simple_report_filename: "RNA_MassHunter_simple_report.xlsx"` を追加してもよい（実装時に既存の `reporting` セクションとの整合性を見て判断）。

### 9.5 `ms2_annotation` セクションの再利用（2026-09-02追記、新規config変更なし）

既存の `config.ms2_annotation` セクション（`enabled`, `mz_tolerance_ppm`, `precursor_match_tolerance_ppm`, `ion_series`, `min_peak_intensity`, `min_relative_intensity_percent`, `max_peaks_per_spectrum`, `use_theoretical_fragments`, `include_neutral_loss`, `include_base_loss`, `min_ion_length` 等）は、新設計のMS2参考情報機能に必要なフィールドを既に全て持っている。**`DEFAULT_CONFIG`/`RunConfig`への新規追加は不要**で、`simple_pipeline.py` はこのセクションをそのまま読む。

`config.ms2_annotation.enabled = false` の場合、簡易版パイプラインはMS2スペクトル抽出自体をスキップし、Mass Comparisonの `ms2_*` 列は全て空欄のまま出力する（MS1のみの初版フローと完全互換）。

### 9.6 `alkaline_phosphatase` 重複定義バグの修正

現行 `config.yaml` は `alkaline_phosphatase:` セクションが62行目・67行目の2箇所に重複記述されている（YAML仕様上は後者のみが有効で実害はないが、今回の変更のついでに1つに統合し、コメントで意図を明記する）。

---

## 10. YAML仕様

`data/modifications.yaml`（既存、変更不要）:

```yaml
schema_version: RNA_MassHunter_modifications_v0.1
modifications:
  - id: m1A
    symbol: m1A
    name: 1-methyladenosine
    target_bases: [A]
    category: biological
    modified_nucleoside_mass_mono: 281.1124
    mass_shift_from_unmodified: 14.0156
    mass_basis: nucleoside_delta
    detectability: {ms1: true, ms2: true}
    isobaric_group: methylation_group
    chemical_group: methylation
    near_isobaric_group: ''
    candidate_policy:
      include_by_mass_search: true
      include_if_position_rule_exists: true
      include_if_literature_supported: true
      include_if_user_specified: true
    source_priority: {...}
    sources: [...]
    curation: {status: manually_checked, notes: "..."}
```

118件収録済み、`modifications.py: load_modifications()` でそのままロード可能。**新設計はこのYAMLを変更しない。** `find_modifications_by_mass_shift()` が参照するのは主に `mass_shift_from_unmodified` フィールドであり、Modification Candidate列挙にそのまま使える。

`data/base_masses.yaml`（既存、変更不要）: `rna_residue_masses`, `constants.water`, `constants.phosphate` を保持。`calculate_unmodified_rna_mass()` が参照する。

---

## 11. 質量計算仕様

- 理論質量: `masses.calculate_unmodified_rna_mass(sequence, base_masses, terminal_form=...)`。`digest_sequence()` が断片ごとに内部で呼び出し、`Fragment.unmodified_mass` に格納済み。
- 観測中性質量: `masses.neutral_mass_from_mz(mz, charge, polarity)`。RNAは基本的に負イオンモード（`polarity="negative"`、`config.instrument.polarity` から取得）。
- ΔDa: `observed_mass - theoretical_mass`（`formula_candidate.calculate_delta_mass()` に集約）。
- Δppm: `ms1_mapping.ppm_error(observed, theoretical)` = `(observed - theoretical) / theoretical * 1e6`。
- 内部計算精度: Python `float`（倍精度、約15〜17桁）をそのまま維持し、途中で丸めない。
- Excel表示精度: mass列は小数点以下3桁（≈0.001 Da）、ΔDa列は小数点以下5桁（Formula Candidateのppmオーダーの一致確認に耐える精度、例 `+0.00003 Da` 表記のため）、Δppm列は小数点以下1〜2桁を基本とする（`number_format` をopenpyxl側で設定）。
- 既定の質量許容差: `±0.01 Da`（`formula_candidate.mass_tolerance_da` / `mass_comparison.mz_tolerance_ppm` として設定可能）。

---

## 12. 元素組成候補アルゴリズム

### 12.1 目的とアルゴリズム方針

ΔDaを説明できる元素組成 `{C: c, H: h, N: n, O: o, P: p, S: s, Se: se}`（各整数、正負両方あり得る＝付加も欠損も対象）を、`sum(|count|) <= max_total_atoms`（既定7）の制約下で網羅探索し、`|exact_mass - ΔDa| <= tolerance_da` を満たす候補を全て列挙する。

### 12.2 探索方式（DFS + 枝刈り）

```python
def enumerate_formula_candidates(mass_difference, tolerance_da, *, max_total_atoms=7,
                                  element_limits=None, elements=_ELEMENT_ORDER):
    element_limits = element_limits or {e: max_total_atoms for e in elements}
    results = []
    counts = {}

    def dfs(idx, remaining_atoms, running_mass):
        if idx == len(elements):
            error = running_mass - mass_difference
            if abs(error) <= tolerance_da:
                comp = ElementalComposition.delta(dict(counts))
                results.append(FormulaCandidateResult(
                    formula=comp.canonical_string(),
                    exact_mass=comp.exact_mass,
                    mass_error_da=error,
                    mass_error_ppm=None,
                    atom_count=sum(abs(v) for v in counts.values()),
                    elements=comp.to_dict(),
                ))
            return
        element = elements[idx]
        limit = min(element_limits.get(element, max_total_atoms), remaining_atoms)
        atomic_mass = MONOISOTOPIC_ATOMIC_MASSES[element]
        for count in range(-limit, limit + 1):
            if count == 0:
                counts.pop(element, None)
            else:
                counts[element] = count
            dfs(idx + 1, remaining_atoms - abs(count), running_mass + atomic_mass * count)
        counts.pop(element, None)

    dfs(0, max_total_atoms, 0.0)
    results.sort(key=lambda r: (abs(r.mass_error_da), r.formula))
    return results[:max_candidates]  # config.formula_candidate.max_candidates_per_match で切り詰め
```

**計算量の見積もり**: 7元素・各元素の探索幅が `2*limit+1`（既定 limit=7 なら最大15）だとしても、`sum(|count|) <= 7` の制約により実効的な組み合わせ数は Stars and Bars 的に大きく絞られる（符号付き7個以下の分割数のオーダーであり、厳密には数千〜数万通り程度。7原子・7元素の整数構成数は高々 `C(14,7) * 2^7` 程度の桁に収まり、実測でのベンチマークをテストに含めることを推奨）。Mass Comparisonの各行（fragment×charge×peakマッチごと）でこの探索を毎回行うとコストが積み上がる可能性があるため、**ΔDaを一定の粒度（例: 0.0001 Da）に丸めたうえで `functools.lru_cache` によりキャッシュする**ことを推奨する（同一ΔDa値が複数行で繰り返し出現するケースが多いため）。

### 12.3 順位付け

- 唯一のランキング基準は `abs(mass_error_da)`（絶対質量誤差）の昇順。
- 原則4「3原子以上でも順位を下げない」を明示的に満たすため、`atom_count` はランキングキーに一切使わない（ソートキーに含めない）。
- 同点（誤差が完全一致するケース、例えば構造的に等価な組成が生じるはずはないが数値誤差で偶然一致する場合）は `formula` 文字列の辞書順で安定ソートし、テストで再現性を保証する。
- `Recommended Formula` = `get_recommended_formula(candidates)`（先頭1件）。
- `Formula Candidates` = 全候補を `"C2H2O1 [+0.00003 Da]; C3H6 [+0.00121 Da]; N2O2 [+0.00342 Da]"` の形式で `; ` 区切り結合。符号は誤差の符号を明示（`+`/`-`）。

### 12.4 Formula CandidateとModification Candidateの分離（原則3の実装保証）

`formula_candidate.py` は `modifications.py` を一切importしない。`mass_comparison.py` 側で両者を並行して呼び出し、`MassComparisonRow` の別々のフィールド（`recommended_formula`/`formula_candidates` と `known_modification`/`modification_candidates`）に格納する。**Excel上でも列を隣接させるが混在させない**（§15, §16）。

---

## 13. 修飾候補検索

```python
def build_modification_candidates(delta_da: float, modifications: list[Modification], tolerance_da: float) -> list[Modification]:
    matches = find_modifications_by_mass_shift(modifications, delta_da, tolerance_da)
    return sorted(matches, key=lambda m: abs(m.mass_shift_from_unmodified - delta_da))
```

- 既存 `find_modifications_by_mass_shift()` をそのまま呼び出し、呼び出し側（`mass_comparison.py`）でソートを追加するだけで完結する（新規アルゴリズムは不要）。
- `Known Modification` = `matches[0]` の表示名（`raw.get("name") or symbol or id`）。
- `Modification Candidates` = 全候補を `"m1A [Δ=+0.0002 Da]; ho5C [Δ=-0.0011 Da]"` 形式で結合（Δは `mass_shift_from_unmodified - delta_da` の符号付き値）。
- `tolerance_da` は `formula_candidate.mass_tolerance_da` と共通の値を使う（Formula/Modification両方で同じ許容差を使うことで、ユーザーがExcel上で両者を同条件で比較できるようにする）。

---

## 14. MS1ピーク処理

- `peak_picking.extract_ms1_peaks(mzml_path, config.reconstruction, warnings)` をそのまま利用（既存の `main.py` と同じ呼び出しパターン、`config.reconstruction` の `rt_min/rt_max/mz_min/mz_max/intensity_threshold` フィールドを流用）。
- 簡易版では `classify_peak_tiers`（Major/Minor/Trace階層化）は**必須としない**（メモに明記が無く、シンプルフローの範囲外と判断）。ただし `config.mass_comparison` に `use_peak_tiers: bool`（既定 `false`）を持たせ、将来ノイズの多いデータで階層フィルタが必要になった場合に既存の `peak_filtering.py` / `classify_peak_tiers` を再利用できる拡張余地を残す。
- Peak IDの採番は `observed_mass.assign_peak_ids()`（§8.3）で行い、`04_Observed_Mass`・`05_Mass_Intensity`・`06_Mass_Comparison` の3シートで同一IDを共有する（原則8）。

### 14A. MS2ピーク処理（参考情報、2026-09-02追記）

- `config.ms2_annotation.enabled` が `true` かつ mzMLにMS2スペクトルが含まれる場合のみ実行する。
- `ms2_annotation.extract_ms2_spectra(mzml_path, config.ms2_annotation, warnings)` でMS2スペクトル一覧を取得（§3.9）。
- `ms2_support.build_ms2_ion_index(fragments, config, base_masses, warnings)` で理論d/w/a/zイオンをfragment_idごとに事前グルーピング（§8.5）。
- Mass Comparisonの各行生成時に `ms2_support.find_ms2_support(...)` を呼び、一致するMS2スペクトルがあれば `ms2_spectrum_id` / `ms2_matched_ion_count` / `ms2_matched_ions` を埋める（§8.2, §15）。
- **この処理はあくまで「参考情報の付加」であり、Recommended Formula / Known Modification の決定には一切影響しない**（原則3の精神を踏襲: Formula/Modification Candidateのランキングはmass_comparison.pyの既存ロジックのみで完結し、MS2の有無で変わらない）。
- MS2スペクトルが1つも無い、または一致するスペクトルが見つからない場合もエラーにはせず、該当列を空欄のまま出力する（§18）。

---

## 15. Mass Comparison

`06_Mass_Comparison` シートの列構成（確定版）:

| 列 | 型 | 由来 |
|---|---|---|
| Peak ID | str | `observed_mass.assign_peak_ids()` |
| Fragment ID | str | `Fragment.fragment_id`（`digest_sequence` が採番） |
| Sequence | str | `Fragment.sequence` |
| Charge | int | マッチ時のcharge |
| Intensity | float | `Peak.intensity` |
| Observed Mass | float | `masses.neutral_mass_from_mz(...)` |
| Theoretical Mass | float | `Fragment.unmodified_mass` |
| ΔDa | float | `observed_mass - theoretical_mass` |
| Δppm | float | `ms1_mapping.ppm_error(...)` |
| Recommended Formula | str | `formula_candidate.get_recommended_formula(...).formula` |
| Formula Candidates | str | `; ` 区切り文字列（§12.3） |
| Known Modification | str | 最有力 `Modification` の表示名 |
| Modification Candidates | str | `; ` 区切り文字列（§13） |
| MS2 Spectrum ID | str | 一致したMS2スペクトルのID（2026-09-02追記、§14A、§8.5） |
| MS2 Matched Ions | str | `; ` 区切りの一致理論イオン一覧、例 `"d3 [-0.5 ppm]; w5 [+1.2 ppm]"`（2026-09-02追記） |

`MS2 Spectrum ID` / `MS2 Matched Ions` は参考情報列であり、`Recommended Formula` や `Known Modification` の判定には影響しない（§14A）。MS2証拠が無い行では両列とも空欄になる。

生成関数 `mass_comparison.build_mass_comparison_rows()` の内部処理順序は §8.2, §7 のデータフロー図を参照。行数が多くなる場合（`reporting.max_excel_rows_per_sheet` 超過）は既存の `_truncate_frame_if_needed` を再利用して切り詰め、警告に記録する。

---

## 16. Excel仕様

### 16.1 シート構成（確定）

```text
01_Index
02_Input
03_Theoretical
04_Observed_Mass
05_Mass_Intensity
06_Mass_Comparison
07_Modifications
08_Visualization   (必要に応じて)
```

### 16.2 各シートの列定義

**02_Input**: 既存 `_flatten_dict(config_dict)` をそのまま利用し、`Parameter` / `Value` の2列でconfig全体をフラット化して出力（既存 `Input_parameters` シートと同じ形式）。

**03_Theoretical**:

| Fragment ID | Sequence | Start | End | Length | Enzyme | Missed Cleavage | Terminal Form | Theoretical Mass |
|---|---|---|---|---|---|---|---|---|

`Fragment` dataclassのフィールドとほぼ1:1対応（`Length` のみ `end - start + 1` で算出）。

**04_Observed_Mass**:

| Peak ID | m/z | Charge | Observed Mass | RT | Scan ID |
|---|---|---|---|---|---|

Chargeは「そのピークがMass Comparisonで何らかのFragmentとマッチした際のcharge」を採用（複数chargeでマッチした場合は最初にマッチした値、または全て別行に展開するかは実装時に決定。**推奨: 1ピークにつき複数charge候補があり得るため、charge単位で行を分ける**）。

**05_Mass_Intensity**:

| Peak ID | m/z | Charge | Observed Mass | Intensity | RT | Scan ID |
|---|---|---|---|---|---|---|

Intensityは独立列として保持（メモ§24の明示的要求。コメントには入れない）。

**06_Mass_Comparison**: §15参照。

**07_Modifications**: `data/modifications.yaml` の内容をそのまま整形出力。

| ID | Symbol | Name | Target Bases | Category | Mass Shift (Da) | Chemical Group | Near-Isobaric Group | Detectability | Curation Status |
|---|---|---|---|---|---|---|---|---|---|

**08_Visualization**: §17参照。

### 16.3 Index / ハイパーリンク

既存 `excel_report.py` の `_add_index_and_backlinks(writer, sheet_names)` と `_autosize_and_freeze(writer)` をそのまま再利用する。パターン:

1. `01_Index` シートに `{"Sheet": name, "Description": ..., "Notes": "Data starts at A3."}` の行を並べる
2. 各データシートは `startrow=2`（3行目からヘッダー開始、1〜2行目を注釈用に空ける既存流儀を踏襲）
3. `_add_index_and_backlinks` でIndex→各シートのハイパーリンク、各シートA1に「← Back to Index」の逆リンクを設定
4. 最後に `_autosize_and_freeze` で列幅自動調整・ウィンドウ枠固定

新規関数 `write_simple_mass_hunter_report()` はこの一連の流儀を完全に踏襲する（既存の `write_excel_report()` 内の該当箇所をそのままコピー呼び出しする形でよい。ロジックの再実装はしない）。

### 16.4 数値表示フォーマット

openpyxlの `cell.number_format` を用いて、mass系列は `"0.000"`、ΔDa系列は `"+0.00000;-0.00000"`、Δppm系列は `"0.0"` を設定する（§11参照）。

---

## 17. Visualization仕様

- 3D visualizationは不採用。
- `08_Visualization` シートに、X軸=Charge、Y軸=Observed Neutral Mass の2Dスキャッタープロットを1枚配置する。
- データソースは `04_Observed_Mass` シートの `Charge` 列・`Observed Mass` 列を直接参照する openpyxl の `ScatterChart` + `Reference` + `Series` を用いる（既存コードに前例が無いため新規実装、リスクは§18参照）。
- Heatmapは必須機能ではなく、実装コストと必要性を見て省略可（メモ§27の指示通り）。
- `config.mass_comparison`（または新設の `visualization` セクション）に `enabled: bool`（既定 `true`）を設け、チャート生成に失敗した場合でも例外を握りつぶしてExcel全体の出力自体は成功させる（`try/except` + warning記録。08シート自体は「データなし」として空でもよい）。

---

## 18. エラー処理

既存の `warnings_manager.add_warning(warnings, level, source, message, context=None)` パターン（`level` は `"INFO"|"WARNING"|"ERROR"`）を踏襲し、以下のケースを警告として記録する。

| ケース | レベル | 対応 |
|---|---|---|
| `config.sequence.sequence` が空 | WARNING | 理論質量・断片生成をスキップ（既存パターンを踏襲） |
| mzML未指定 / 読み込み失敗 | WARNING | MS1ピーク取得・Mass Comparisonをスキップし、Excelは他シートのみで出力 |
| `formula_candidate` 探索でtolerance内に候補0件 | INFO | `Recommended Formula`/`Formula Candidates` を空文字列にし、行自体は残す |
| `modifications` 探索でtolerance内に候補0件 | INFO | 同上（Modification側） |
| MS2スペクトル抽出失敗、またはmzMLにMS2スペクトルが存在しない（2026-09-02追記） | INFO | `ms2_annotation.extract_ms2_spectra` 自体が例外を握りつぶしてwarningを積む既存挙動をそのまま利用。`ms2_*` 列は空欄のまま出力を継続 |
| 一致するMS2前駆体スペクトルが複数あり自動選択した（2026-09-02追記） | INFO | 採用したスペクトルの `spectrum_id` のみ出力し、選択根拠（強度最大）は警告に記録 |
| `formula_candidate.element_limits` の合計が `max_total_atoms` と矛盾（例: 全要素0） | ERROR | config起動時バリデーションで検出、`validate_config()` に追加 |
| Excel出力時に `08_Visualization` のチャート生成が例外を送出 | WARNING | チャートなしで出力を継続（§17参照） |
| 未知元素（`element_limits` に `C/H/N/O/P/S/Se` 以外のキー） | ERROR | config起動時バリデーションで検出 |
| Peak ID採番の重複（理論上発生しないはずだが防御的に） | ERROR | assertまたは例外、テストでカバー |

すべての警告は既存の `Warnings` シート相当（新設計では `02_Input` シート末尾、または専用の軽量シートを追加するかは実装時に判断。**推奨: 既存の慣例に合わせて `02_Input` の下部に警告テーブルを併記**）に出力する。

---

## 19. 単体テスト

既存のテスト文化（pytest、関数ベース、`@pytest.mark.parametrize` 多用、`pytest.approx` で浮動小数比較、フィクスチャは最小限）を踏襲する。

### `test_formula_candidate.py`

```python
import pytest
from rna_masshunter.formula_candidate import enumerate_formula_candidates, get_recommended_formula

def test_single_oxygen_delta():
    candidates = enumerate_formula_candidates(15.9949, tolerance_da=0.001)
    assert get_recommended_formula(candidates).formula == "O1"

def test_ch2_delta():
    candidates = enumerate_formula_candidates(14.0157, tolerance_da=0.001)
    assert get_recommended_formula(candidates).formula in ("C1H2",)

@pytest.mark.parametrize("max_total_atoms,expected_atom_count_upper_bound", [(7, 7), (3, 3)])
def test_max_total_atoms_respected(max_total_atoms, expected_atom_count_upper_bound):
    candidates = enumerate_formula_candidates(50.0, tolerance_da=0.5, max_total_atoms=max_total_atoms)
    assert all(c.atom_count <= expected_atom_count_upper_bound for c in candidates)

def test_ranking_not_biased_by_atom_count():
    # 3原子以上でも質量誤差が小さければ上位に来ることを明示的に検証（原則4）
    candidates = enumerate_formula_candidates(42.0106, tolerance_da=0.01)  # 例: acetylation相当
    assert candidates == sorted(candidates, key=lambda c: abs(c.mass_error_da))

def test_zero_delta_returns_empty_formula():
    candidates = enumerate_formula_candidates(0.0, tolerance_da=0.0001)
    assert get_recommended_formula(candidates).formula == "0"

def test_negative_delta_allows_loss_formula():
    candidates = enumerate_formula_candidates(-15.9949, tolerance_da=0.001)
    assert get_recommended_formula(candidates).formula == "O-1"

def test_no_candidates_within_tolerance_returns_empty_list():
    candidates = enumerate_formula_candidates(123.456, tolerance_da=0.0001, max_total_atoms=2)
    assert candidates == []
    assert get_recommended_formula(candidates) is None
```

### `test_mass_comparison.py`

- 小規模フィクスチャ（3〜5個の合成 `Fragment` / `Peak` / `Modification`）を用意し、`build_mass_comparison_rows()` の出力行数・ΔDa/Δppm計算値・Recommended列の一致度順序を検証。
- Formula CandidateとModification Candidateが互いに独立して計算され、片方が0件でももう片方は影響を受けないことを検証（原則3）。
- （2026-09-02追記）`ms2_spectra`/`ms2_ion_index` を渡さない場合、`ms2_*` 列が常に既定値（None/0/""）になり、Recommended Formula/Known Modificationの値がMS2の有無で変化しないことを検証。

### `test_ms2_support.py`（2026-09-02追記）

```python
import pytest
from rna_masshunter.ms2_support import build_ms2_ion_index, find_ms2_support

def test_no_matching_precursor_returns_none(synthetic_fragments, synthetic_ms2_spectra, config, base_masses):
    ion_index = build_ms2_ion_index(synthetic_fragments, config, base_masses)
    result = find_ms2_support("FRAG_NOT_PRESENT", 3, 999.0, synthetic_ms2_spectra, ion_index, config)
    assert result is None

def test_matched_ions_sorted_by_ppm_error(synthetic_fragments, synthetic_ms2_spectra, config, base_masses):
    ion_index = build_ms2_ion_index(synthetic_fragments, config, base_masses)
    result = find_ms2_support(synthetic_fragments[0].fragment_id, 2, <precursor_mz>, synthetic_ms2_spectra, ion_index, config)
    assert result is not None
    assert result.matched_ion_count > 0

def test_ms2_disabled_short_circuits(synthetic_fragments, synthetic_ms2_spectra, config, base_masses):
    config.ms2_annotation["enabled"] = False
    # simple_pipeline.run() 側でms2_spectra=[]になることを想定したテスト
    # (ms2_support自体はconfig.ms2_annotation.enabledを見ないため、
    #  呼び出し側=simple_pipeline.pyでスキップされることを別途統合テストで検証)
```

---

## 20. 統合テスト

`test_simple_pipeline.py`:

- 合成RNA配列（数十塩基程度）と、`pyteomics` または手書きXMLで生成した最小限の合成mzML（3〜5ピーク、既知の修飾質量シフトを意図的に含める）を使い、`simple_pipeline.run()` をエンドツーエンドで実行する。
- 出力Excelファイルが生成されること、各シートが存在すること（`openpyxl.load_workbook` で確認）、Mass ComparisonシートのRecommended Formula/Known Modificationが期待値と一致することを検証。
- 既存 `test_main_config_path.py` のconfigパス解決パターンを参考にしてよい。
- mzML生成が手間な場合は、`peak_picking.extract_ms1_peaks` をモック（`unittest.mock.patch`）して `Peak` オブジェクトのリストを直接注入するテストでも可（既存テストのスタイルに合わせて判断）。
- （2026-09-02追記）合成mzMLにMS2スペクトル（ms level=2、既知のd/w/a/zイオンm/zを含む）も含め、出力Excelの `06_Mass_Comparison` の `MS2 Matched Ions` 列が期待通り埋まることを検証する。`config.ms2_annotation.enabled = false` のケースも別途実行し、MS1のみの初版フローと同じ出力になる（回帰しない）ことを確認する。

---

## 21. 実mzMLによる動作確認

コードでは自動化できない、人手での最終確認手順:

1. ユーザーの実際のmzMLファイルと、実際に使っている `config.yaml`（相当）を用意し、`simple_main.py --config <path>` を実行。
2. 生成されたExcelを開き、`01_Index` の全リンクが正しく遷移すること、各シートの「← Back to Index」リンクが機能することを目視確認。
3. `06_Mass_Comparison` で、ユーザーが既に知っている（wet lab等で確認済みの）修飾を含む断片について、`Known Modification` 列が期待通りの修飾名を返しているかを確認。
4. 同じ行の `Formula Candidates` 列が、その既知修飾の元素組成差分と整合するか（例: メチル化なら `CH2` が候補に含まれるか）を確認。
4.5.（2026-09-02追記）MS2スペクトルを含むmzMLの場合、`MS2 Matched Ions` 列がユーザーの既知の知見（例えば特定断片でd3イオンが強く観測されるはず等）と大きく矛盾していないかを確認。あくまで参考情報である旨を再確認し、`Recommended Formula`/`Known Modification` の判定がMS2の有無で変わっていないことも確認する。
5. `08_Visualization` の散布図が正しく charge×observed neutral mass をプロットしているか確認。
6. 大きなmzML（ピーク数が多い場合）で処理時間が実用的な範囲か計測し、必要なら `formula_candidate` のキャッシュ設定（§12.2）を調整。
7. Excel（Microsoft Excel / LibreOffice Calc 双方、ユーザーの利用環境に応じて）でファイルが正しく開けるかを確認（既存コードは `openpyxl` エンジンを使用しており通常は問題ないが、ハイパーリンクやチャートは念のため確認）。

---

## 22. 実装順序

既存の巨大な `main.py` の `main()` は分岐が複雑に絡み合っており部分抽出が困難なため、**新設計は独立した新規オーケストレーション層（`simple_pipeline.py` + `simple_main.py`）として実装し、`main.py` 自体には一切手を入れない**方針を確定する（既存パイプラインとの並行運用・段階的検証がしやすく、既存の巨大な回帰リスクを避けられるため）。

| フェーズ | 内容 | 完了条件 |
|---|---|---|
| 0 | `config.py` / `models.py` に `formula_candidate` セクション・`cca_processing.registered_sequence_cca_mode` を追加。`config.yaml` にサンプル値を追記、重複バグ修正 | config読込テストが通る |
| 1 | `formula_candidate.py` 実装 | `test_formula_candidate.py` 全通過（原則4の検証含む） |
| 2 | `mass_comparison.py` 実装（`observed_mass.py` のPeak ID採番含む、MS1のみで完結する版） | `test_mass_comparison.py` 全通過 |
| 2.5（2026-09-02追記） | `ms2_support.py` 実装（`build_ms2_ion_index`, `find_ms2_support`）、`mass_comparison.py` に `ms2_spectrum_id`/`ms2_matched_ion_count`/`ms2_matched_ions` 列を追加 | `test_ms2_support.py` 全通過、既存 `test_mass_comparison.py` が回帰しないこと（MS2未指定時は空欄のまま） |
| 3 | CCA処理の配線（`cca_processing.process_cca_tail` を `simple_pipeline.py` から呼び出す） | 単体で処理結果を確認 |
| 4 | `simple_pipeline.py` / `simple_main.py` 実装（既存の消化・MS1ピーク抽出・MS2スペクトル抽出・修飾検索を接続） | configのみでドライラン成功（mzMLなしでも警告付きで完走） |
| 5 | `excel_report.py` に `write_simple_mass_hunter_report()` を追加（01〜07シート、MS2列含む。08は次フェーズ） | 生成Excelを目視確認、Index/リンクが機能、MS2列が正しく埋まる |
| 6 | `08_Visualization`（散布図）を追加 | チャートが正しく表示される、失敗時もExcel全体は出力される |
| 7 | `test_simple_pipeline.py`（統合テスト、MS2ありなし両ケース） | 合成データでエンドツーエンド成功 |
| 8 | 実mzMLによる動作確認（§21、MS2確認含む） | ユーザーによる目視確認完了 |
| 9 | README更新・コミット/プッシュ | ユーザーの既存運用（GitHubへの都度保存）に沿って完了 |

各フェーズはそれぞれ独立してコミット可能な単位とし、フェーズ完了ごとにテストを通してから次に進む。

---

## 23. Claude Codeへの実装ルール

1. **既存の巨大な `main()` およびSCIEX監査系（`sciex_*.py`）・shadow audit系・MS2アノテーション系・evidence ranking系には一切手を入れない。依存も追加しない。**
2. Formula CandidateとModification Candidateのロジック・出力列は完全に分離すること（`formula_candidate.py` は `modifications.py` をimportしない）。
3. 元素組成候補のランキングで原子数（`atom_count`）を順位付けキーに使わないこと（3原子以上でも順位を下げない、原則4）。この制約はテスト（`test_ranking_not_biased_by_atom_count`）で明示的に担保する。
4. 新規モジュールは既存モジュールをimportして再利用し、ロジックを複製しないこと（特に `masses.py`, `elemental_composition.py`, `mass_shift_ms1_search.py`, `modifications.py` の関数をコピーせず直接呼び出す）。
5. config変更は `DEFAULT_CONFIG`（`config.py`）と `RunConfig`（`models.py`）の両方に反映し、`validate_config()` にも対応するバリデーションを追加すること。既存の `if not isinstance(...): raise ValueError(...)` パターンに合わせる。
6. 既存の `warnings_manager.add_warning()` パターン、`dataclass` + 型ヒントのスタイル、pytest関数ベース・`parametrize`・`pytest.approx` のテストスタイルに従うこと。
7. 内部計算はPythonの`float`（倍精度）を維持し、途中で丸めないこと。丸めはExcel表示（`number_format`）でのみ行う。
8. 数値以外の一致度順ソート（Formula Candidates / Modification Candidates の文字列生成）は、常に「質量差の絶対値昇順」を基準にすること。
9. `excel_report.py` の変更は新規関数の追加に限定し、既存の `write_excel_report()` 本体・既存シート生成ロジックには触れないこと。ただし `_add_index_and_backlinks` 等の共通ヘルパーは積極的に再利用すること。
10. 各フェーズ（§22）ごとに小さくコミットし、フェーズ完了時に対応するテストが通っていることを確認してから次のフェーズに進むこと。テストが落ちている状態でコミットしないこと。
11. `config.yaml` の `alkaline_phosphatase` 重複キーバグの修正は今回のスコープに含めてよいが、それ以外の既存機能（SCIEX系・MS2系・intact reconstruction系等）には触れないこと。
12. 実装中に本仕様書との齟齬や未確定事項（例: §8.2 の Peak ID とcharge展開の扱い、§16.2 の 04シートのcharge重複行の扱い）に気づいた場合は、独断で仕様を変更せず、いったん立ち止まってユーザーに確認すること。
13. コミットメッセージは既存リポジトリの慣例（日本語の簡潔な要約）に合わせること。
14. （2026-09-02追記）MS2関連の実装は §3.9・§8.5・§14A に明記した範囲（`extract_ms2_spectra`, `generate_theoretical_ms2_ions` の再利用、未修飾d/w/a/zイオンとの一致件数表示のみ）に厳密に留めること。`match_ms2_spectra`, `find_parent_candidates`, `build_fragment_evidence` や `evidence_ranking.py`, `biological_context.py`, `modification_hypothesis_audit.py` など、前駆体レスキュー・エビデンスレベル判定・生物学的コンテキスト評価に関わるモジュールは一切importしないこと。MS2列はあくまで参考情報であり、Formula Candidate/Modification Candidateのランキングや推奨値の決定ロジックに一切影響を与えないこと（テストで担保、§19）。

---

## 未確定事項・要確認リスト（実装着手前にユーザーへ確認推奨）

- §9.2: `mass_comparison` セクションを独立させるか、`fragment_mapping` を流用するか
- §15/16.2: 1つのPeakが複数chargeでマッチした場合の `04_Observed_Mass` の行展開方法
- §9.4: 出力Excelのファイル名・出力先の命名規則
- §12.2: `element_limits` の初期値（既存 `modifications.yaml` の分布を見て決定するとメモ§16に明記あり、実装フェーズ0で要確認）
- digestion.py内部の `_generate_three_prime_tail_candidates()`（3'末端フラグメントの末端形候補展開、AP関連）と、新規配線する `cca_processing.process_cca_tail()`（配列全体へのCCA成熟化）が役割上重複しないか、実装フェーズ3で実コード上の相互作用を確認すること
- （2026-09-02追記）§8.5 `find_ms2_support` のタイブレーク方針（同一tolerance内に複数のMS2前駆体スペクトルが一致した場合、強度最大のものを機械的に採用する設計）で実運用上問題ないか。RTの近さも考慮すべきか等はユーザーの実データを見て最終判断する
- （2026-09-02追記）`ion_series` に将来 b/c/x/y 系列を追加したいニーズがあるか（現状は `masses.fragment_ion_series_offsets` 内で「将来のCID/HCD refinement用」とコメントされているのみで、d/w/a/z系列のみが実運用対象）。今回のMS2参考情報機能では既存の既定値 `["d","w","a","z"]` をそのまま踏襲する前提としている
