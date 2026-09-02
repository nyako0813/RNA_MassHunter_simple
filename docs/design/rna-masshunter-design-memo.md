# RNA_MassHunter 設計引き継ぎメモ（ChatGPTによる設計、2026-09-02受領）

> このドキュメントは nyako さんがChatGPTとの設計セッションで作成し、Claude Codeへの実装引き継ぎ用に共有したメモの原文です。今後この内容を基にリポジトリ調査・実装仕様書作成・実装を進めます。

## 1. プロジェクトの目的
RNA_MassHunterは、tRNA/RNAのLC-MS/MSデータを解析し、
- RNA配列から理論的なRNase断片を生成
- 各断片の理論質量を計算
- mzMLから観測されたMS1ピークを取得
- m/zと電荷から観測中性質量を計算
- 理論質量と観測質量を比較
- 質量差（ΔDa / Δppm）を計算
- 質量差から元素組成候補を列挙
- 既知RNA修飾候補を列挙
- それらをExcelに整理
することを目的とする。

### 基本思想
RNA_MassHunter自身が修飾を「確定同定」することを主目的にはしない。
**MSデータを正確に整理し、人間が修飾を判断できる情報を提供すること**を主目的とする。
複雑な生物学的スコアリングや自動推論は主解析フローから外す。

---

## 2. 今回確定した解析フロー
```text
RNA sequence
    ↓
CCA processing
    ↓
RNase A / RNase T1
    ↓
missed cleavage
    ↓
Alkaline phosphatase (AP)
    ↓
theoretical fragments
    ↓
unmodified theoretical mass
    ↓
mzML MS1
    ↓
observed m/z
    ↓
charge
    ↓
observed neutral mass
    ↓
mass comparison
    ↓
    ┌─────────────────────┐
    │                     │
    ↓                     ↓
Formula Candidate    Modification Candidate
    │                     │
    ↓                     ↓
Recommended Formula  Known Modification
    │                     │
    └──────────┬──────────┘
               ↓
            Excel
```

---

## 3. CCA処理
既存のCCA処理機能を再利用する。
現在の仕様：
- ゲノム上のtRNA末端が `C` → `CA` を付加して `CCA`
- ゲノム上の末端が `CC` → `A` を付加して `CCA`
- すでに `CCA` → 変更しない

既存の `cca_processing.py` の処理を基本的に再利用する。CCA機能を作り直さない。

---

## 4. RNase処理
RNase A / RNase T1の既存機能を利用する。
ユーザー設定として、
- RNase A：ON/OFF
- RNase T1：ON/OFF
を設定できるようにする。

既存の `digestion.py` / `enzymes.py` をできるだけ再利用する。

---

## 5. Missed cleavage
missed cleavageは維持する。ユーザーが設定したmissed cleavage数に応じて理論断片を生成する。既存のdigestion機能を再利用する。

---

## 6. Alkaline Phosphatase
AP処理は独立した処理オプションとする。ON/OFFで制御できるようにする。RNase処理とAP処理を混同しない。

---

## 7. 理論質量
既存の質量計算機能を再利用する。特に `masses.py` の既存機能を利用する。

既存の主な機能：
- `calculate_unmodified_rna_mass()`
- `neutral_mass_from_mz()`
- `mz_from_neutral_mass()`
- `elemental_delta_mass()`

既存のmonoisotopic mass定義も再利用する。

既知の定数：
```text
PROTON_MASS = 1.007276466812
DEFAULT_WATER_MASS = 18.010564684
DEFAULT_PHOSPHATE_MASS = 79.966331
```

元素質量として既に `C, H, N, O, P, S, Se` が扱える。

---

## 8. 観測MS1
mzMLからMS1スペクトルを取得する。既存の `mzml_reader.py` / `peak_picking.py` を再利用する。

既存機能：
- `iter_spectra()`
- `extract_ms1_spectra()`
- `extract_ms2_spectra()`
- `extract_ms1_peaks()`

MS1ピークから、m/z / intensity / retention time / scan ID などを取得する。

---

## 9. 電荷と観測中性質量
ユーザーが希望している可視化では、横軸：charge、縦軸：observed neutral mass とする。そのため、MS1の観測m/zとchargeからneutral massを計算する。既存の `neutral_mass_from_mz()` を再利用する。

---

## 10. Mass Comparison
解析の中心となるシート。`06_Mass_Comparison` を予定。

基本的な列：

| 列 | 内容 |
| --- | --- |
| Peak ID | 観測ピーク識別子 |
| Fragment ID | 理論断片識別子 |
| Sequence | RNA断片配列 |
| Charge | 電荷 |
| Intensity | 強度 |
| Observed Mass | 観測中性質量 |
| Theoretical Mass | 理論質量 |
| ΔDa | Observed - Theoretical |
| Δppm | 質量差ppm |
| Recommended Formula | 最も合致する元素組成 |
| Formula Candidates | 元素組成候補一覧 |
| Known Modification | 最も合致する既知修飾 |
| Modification Candidates | 既知修飾候補一覧 |

---

## 11. 元素組成候補
今回新たに追加する重要機能。

### 目的
質量差から、「この質量差はどのような元素組成の追加で説明できるか」を列挙する。**修飾名を推測する機能ではない。**

例えば、
```text
ΔDa = +15.9949 → O
ΔDa = +14.0157 → CH2
```
のようにする。

---

## 12. 元素組成候補と修飾候補は完全に分離
重要。

### Formula Candidate
質量差から計算した元素組成。
```text
ΔDa
 ↓
C2H2O
```

### Modification Candidate
`modifications.yaml` に登録されている既知修飾。
```text
ΔDa
 ↓
既知修飾候補
```

両者を同じ候補リストに混ぜない。

---

## 13. Recommended Formula
元素組成候補の中から、**質量差に最も合致するものを1つ** `Recommended Formula` に入れる。3原子以上だからといって優先度を下げない。

例えば、`C2H2O` `C3H6` `N2O2` が候補なら、原子数ではなく質量差への一致度を基本として最も合致するものをRecommended Formulaとする。

---

## 14. Formula Candidates
候補をセル内に合致度順で並べる。

例えば、
```text
C2H2O [+0.00003 Da]; C3H6 [+0.00121 Da]; N2O2 [+0.00342 Da]
```
のような形式。

つまり、`Recommended Formula` → 最有力候補1つ、`Formula Candidates` → 候補を合致度順に列挙 とする。

---

## 15. 元素数制限
元素組成候補探索は無制限にしない。**合計7原子以下** を基本制限とする。
```text
C + H + N + O + P + S + Se <= 7
```
これは「7原子を超える化合物は物理的に存在しない」という意味ではない。あくまで、**質量差だけからRNA修飾候補として探索する範囲を7原子までに限定する** という解析上の制限。

---

## 16. 元素ごとの上限
合計7原子だけでは候補数が増えすぎる可能性があるため、元素ごとの上限も設定可能にする。

概念例：
```yaml
formula_candidate:
  enabled: true
  max_total_atoms: 7
  elements:
    C: 3
    H: 7
    N: 3
    O: 4
    P: 1
    S: 1
    Se: 1
```
ただし、この具体的な数値は実装時に既存コードと修飾データを確認して最終決定する。

---

## 17. 元素組成候補の順位付け
3原子以上だから優先度を下げる処理は行わない。基本的には、
1. ΔDaとの絶対質量誤差
2. 必要に応じて追加の合理性条件
を使う。ただし、これは**化学構造の同定スコアではない**。質量差に対する元素組成候補のランキングである。

---

## 18. Known Modification
既存の `modifications.yaml` を使用する。既存の `modifications.py` に `find_modifications_by_mass_shift()` が存在する。これを再利用する。

既知修飾については、修飾名 / mass shift / category / target bases / detectability / curation / sources / candidate policy / chemical group / near-isobaric group など既存データを利用する。

---

## 19. Modification Candidates
`modifications.yaml` の既知修飾のうち、ΔDaに合致するものを列挙する。

例えば、
```text
Modification Candidates
modification A [Δ=+0.0002 Da]; modification B [Δ=-0.0011 Da]; modification C [Δ=+0.0034 Da]
```
のように、質量差への合致度順に並べる。

---

## 20. Known Modification
既知修飾候補のうち、最も合致するものを1つ別セルに表示する。ただし、
**Known Modification = 自動的な修飾同定結果**
とは扱わない。あくまで、
**質量差が最も近い既知修飾候補**
である。

---

## 21. Excel構成
基本構成：
```text
01_Index
02_Input
03_Theoretical
04_Observed_Mass
05_Mass_Intensity
06_Mass_Comparison
07_Modifications
08_Visualization
```
`08_Visualization`は必要に応じて使用。

---

## 22. Index
最初のシートをIndexにする。各シートへのリンクを作成する。また各シートからIndexへ戻れるリンクを設置する。

---

## 23. Observed_Mass
`04_Observed_Mass`
観測されたneutral massを整理する。可視化では、X axis = Charge、Y axis = Observed neutral mass とする。Peak IDを表示して、元データとの追跡性を確保する。

---

## 24. Mass_Intensity
`05_Mass_Intensity`
Intensityはコメントではなく、**独立した列**として保存する。

例：
| Peak ID | m/z | Charge | Observed Mass | Intensity | RT | Scan ID |
| ------- | --: | -----: | ------------: | --------: | -: | ------- |

---

## 25. Theoretical
`03_Theoretical`
理論断片情報を保存。

想定列：
| Fragment ID | Sequence | Start | End | Length | Enzyme | Missed Cleavage | Terminal Form | Theoretical Mass |
| ----------- | -------- | ----: | --: | -----: | ------ | --------------: | ------------- | ---------------: |

---

## 26. Modifications
`07_Modifications`
既知修飾情報を整理する。`modifications.yaml` のデータを必要に応じて出力。

---

## 27. Visualization
3D visualizationは採用しない。必要なら2D visualizationとして、`Charge × Observed Neutral Mass` を使用する。Heatmapについては必須機能ではなく、必要性が低ければ削除可能。

---

## 28. 既存コードで再利用するもの
以下は基本的に再利用する。
```text
rna_masshunter/cca_processing.py
rna_masshunter/cca_tail_state.py
rna_masshunter/digestion.py
rna_masshunter/enzymes.py
rna_masshunter/masses.py
rna_masshunter/elemental_composition.py
rna_masshunter/mzml_reader.py
rna_masshunter/peak_picking.py
rna_masshunter/modifications.py
```
特に、`elemental_composition.py` は元素組成候補生成に利用できる可能性が高い。既存の `ElementalComposition` を再利用し、同じ元素質量定義を二重管理しない。

---

## 29. 部分的に修正する可能性があるもの
```text
config.py
models.py
ms1_mapping.py
excel_report.py
main.py
```
ただし、実際の現在コードを確認してから変更箇所を確定する。

---

## 30. 新規作成候補
主に、
```text
mass_comparison.py
formula_candidate.py
```
を新規作成する。必要に応じて、`observed_mass.py` などを追加する。ただし、既存機能で十分なら不要なモジュールは増やさない。

---

## 31. formula_candidate.py の役割
概念的には、`calculate_delta_mass()` / `enumerate_formula_candidates()` / `rank_formula_candidates()` / `get_recommended_formula()` などを担当。

入力：`mass_difference` / `tolerance` / `element constraints` / `max_total_atoms`

出力：`formula` / `exact_mass` / `mass_error` / `atom_count` など。

---

## 32. 数値精度
内部計算は高精度を維持する。少なくとも6桁以上の精度を保持する。Excel表示では、mass → 約0.001 Da、ΔDa → 約0.001 Da、Δppm → 適切な小数精度 程度を基本とする。

デフォルトの質量許容差は概念上 ±0.01 Da 程度を候補とするが、最終的には設定可能にする。

---

## 33. 重要な設計原則
1. 既存の動いている機能をできるだけ再利用する。
2. 自動同定を目的にしない。
3. 元素組成候補と修飾候補を分離する。
4. 3原子以上だからといって候補順位を下げない。
5. 元素組成候補は質量差との一致度を基本に並べる。
6. 候補探索は最大7原子に制限する。
7. 候補が複数あることを前提とする。
8. Peak ID / Fragment IDでデータを追跡可能にする。
9. Excelだけを見ても解析結果を理解できるようにする。

---

## 34. 最終的にClaude Codeへ渡すもの
最終成果物として、**「RNA_MassHunter 再設計・実装仕様書」** を作成する。Claude Codeがそのまま実装に入れるレベルまで具体化する。

設計書には、
1. プロジェクト目的
2. 現在のコード構成
3. 現在存在する関数
4. 再利用する関数
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
を含める。

---

## 35. 現時点でコード変更はしない
現在は設計段階。まだコードを変更しない。次のステップとして、GitHub上の現在のRNA_MassHunterについて、
```text
models.py
config.py
main.py
digestion.py
cca_processing.py
elemental_composition.py
ms1_mapping.py
excel_report.py
tests/
```
などを詳細確認する。その上で、**「現在のコードのどこを使い、どこを変更し、どこに新規コードを追加するか」** を関数単位で確定する。その結果を基に、最終的なClaude Code用設計書を作成する。

---

## 現時点での最重要確定事項

### Formula Candidate
**ΔDa → 元素組成候補**
- 修飾名ではない
- C/H/N/O/P/S/Seを基本対象
- 合計7原子以下
- 3原子以上でも順位を下げない
- 質量差への一致度順
- 最有力をRecommended Formula
- 全候補をFormula Candidates

### Modification Candidate
**ΔDa → modifications.yamlの既知修飾候補**
- 元素組成候補とは別管理
- 質量差への一致度順
- 最有力をKnown Modification
- 全候補をModification Candidates

### 最終的な比較
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

この設計を今後の基準仕様とする。
