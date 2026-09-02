# Claude Code 実装引き継ぎ指示書 — RNA_MassHunter簡易版

作成日: 2026-09-02
このファイルは、設計仕様書（`RNA_MassHunter_再設計_実装仕様書.md`、同ディレクトリ）を実装に移すために、ローカル環境で動くClaude Code（またはCodex）セッションへ渡すための引き継ぎ指示書です。**まずこのファイルを読み、次に設計仕様書を読んでください。** 設計仕様書が「何を作るか」の仕様、この指示書は「どこから・どういう手順で始めるか」の手順書です。

もしこのメッセージを読んでいるのが人間（nyakoさん）であれば: このファイル一式（`rna_masshunter_simple/` ディレクトリ全体）をローカル環境に展開し、Claude Codeセッションを開始してこのファイルのパスを最初のメッセージで渡してください。

---

## 0. 現状（すでに完了していること）

Claude（Cowork、クラウド環境）側で、`nyako0813/RNA_MassHunter`（既存の大規模リポジトリ）から今回の簡易版が必要とする部分だけを移植し、**動作確認済み**の状態でこの `rna_masshunter_simple/` ディレクトリ一式を用意しました。

- `rna_masshunter/` パッケージに、再利用対象モジュール（`masses.py`, `elemental_composition.py`, `enzymes.py`, `digestion.py`, `cca_processing.py`, `mzml_reader.py`, `mzml_diagnostics.py`, `peak_picking.py`, `modifications.py`, `mass_shift_ms1_search.py`, `ms1_mapping.py`, `warnings_manager.py`）をそのままコピー済み。
- `rna_masshunter/models.py` と `rna_masshunter/config.py` は、元の巨大なファイル（SCIEX監査系など今回不要なセクションを多数含む）から**必要なフィールドだけを抜き出して新規に書き直した**もの。
- `rna_masshunter/cca_processing.py` は、元は `cca_tail_state.py` から `RegisteredSequenceCCAMode` をimportしていましたが、**`cca_tail_state.py` は `intact_rna_mass.py` → `structure_fragment.py` → `backbone_state.py`/`modification_composer.py`等、今回使わない重量級モジュール一式を連鎖的に読み込んでしまう**ことが判明したため、このEnum（3値のみ）をファイル内に直接定義する形に書き換えてあります。
- `rna_masshunter/ms2_extraction.py` は新規ファイルで、元の `ms2_annotation.py` から `extract_ms2_spectra` / `generate_theoretical_ms2_ions` の2関数だけを抜き出したものです。**元のファイルをそのままimportすると `base_loss_masses`, `base_loss_ions`, `modified_precursor`, `modified_fragment_ions`, `ms2_unmatched_audit`, `ms2_zero_intensity_audit` まで芋づる式に読み込まれてしまう**ため、これも避けています。またゼロ強度監査用の記録処理（`capture_source_spectrum` 等）の呼び出しも、今回使わない機能なので削除してあります。
- `data/modifications.yaml`（118件収録の既知修飾データベース）と `data/base_masses.yaml`（元素質量定数）をコピー済み。
- `config.yaml` サンプルと、対応する `rna_masshunter/config.py`（`DEFAULT_CONFIG` を今回必要なセクションだけに絞って新規作成）を用意済み。
- `tests/test_elemental_composition.py`（元リポジトリからそのまま移植）と `tests/test_vendored_smoke.py`（今回新規、移植した各モジュールが実際に動くことを確認するスモークテスト）を用意し、**`pytest` で17件全て通過することをクラウド環境で確認済み**です。

**つまりPhase 0（設計仕様書§22の「移植」部分）は完了しています。** ここから先（`formula_candidate.py` 以降の新規ロジック）が実装対象です。

---

## 1. 実行環境: WSL2 (Ubuntu) を推奨

Windows/WSLどちらでも動かせますが、**WSL2上のUbuntu**を推奨します。理由:

- 依存パッケージ（`pyteomics`, `lxml`, `numpy`, `openpyxl` など）はLinux環境での実績が長く、Windowsネイティブより安定します。
- 元の大規模リポジトリの開発履歴に「Windows非互換のimport resourceを修正」というコミットがあり、**このコードベースは過去にWindowsネイティブ実行で一度つまずいている**実績があります（`resource` はUnix専用の標準ライブラリで、Windowsには存在しません）。今回移植したモジュール群はその問題箇所を含んでいないはずですが、今後のコード追加で同種の問題を踏むリスクを避けるため、最初からLinux環境で統一するのが安全です。
- Claude Codeやgit、pyenv/venvなど、Python開発でよく使われるツール群がLinux前提で書かれていることが多く、初心者の方がトラブルシューティング情報を探しやすい環境です。

セットアップ手順（WSL未導入の場合）:
```powershell
# Windows PowerShell（管理者権限）で
wsl --install -d Ubuntu
```
インストール後、WSL内のUbuntuで:
```bash
sudo apt update && sudo apt install -y python3-venv python3-pip git
```

どうしてもWindowsネイティブで進めたい場合は、`requirements.txt` の各パッケージがWindows用wheelを持っているか（`pip install` 時にビルドエラーが出ないか）を早い段階で確認してください。問題が出た場合はWSLへの切り替えを検討してください。

---

## 2. GitHubリポジトリの作成

このプロジェクト専用の**新規リポジトリ**を作成します（既存の `nyako0813/RNA_MassHunter` とは別リポジトリ。設計方針としてユーザーの承認済み — 設計仕様書冒頭の改訂メモおよび本指示書の性質を参照）。

リポジトリ名の案: `RNA_MassHunter_simple`（この指示書ではこの名前を使いますが、お好みで変更してかまいません）。公開範囲（public/private）はユーザーに確認してください（未公開の研究用ツールであることが多いため、迷う場合は `--private` を既定にすることを推奨します）。

### `gh` CLIが使える場合（推奨）
```bash
cd ~/  # 作業ディレクトリへ
# このディレクトリ(rna_masshunter_simple/)がまだ展開されていなければ、
# 送付されたzipを展開してから作業する
gh auth status   # 未認証なら `gh auth login` を先に実行
gh repo create RNA_MassHunter_simple --private --source=rna_masshunter_simple --remote=origin --push
```
（`gh repo create ... --source --push` は、ローカルディレクトリをそのまま新規GitHubリポジトリの初期コミットとしてpushします。事前に `git init` / `git add` / `git commit` を済ませておくか、`gh repo create` に任せる場合は先に `cd rna_masshunter_simple && git init && git add -A && git commit -m "..."` を済ませてから実行してください。）

### `gh` CLIが無い/使えない場合
1. https://github.com/new でリポジトリを手動作成（Owner: `nyako0813`、Repository name: `RNA_MassHunter_simple`、Public/Privateを選択、READMEなどの自動生成チェックは全て外す）。
2. ローカルで:
```bash
cd rna_masshunter_simple
git init
git add -A
git commit -m "初期コミット: 既存リポジトリからの移植分＋設計ドキュメント"
git branch -M main
git remote add origin https://github.com/nyako0813/RNA_MassHunter_simple.git
git push -u origin main
```

**注意**: このクラウドセッション（Claude Cowork）は nyakoさんのGitHub認証情報を持っていないため、リポジトリの作成・pushはローカル環境（WSL上のClaude Code）で実行してください。

---

## 3. ディレクトリ構成（現状）

```text
rna_masshunter_simple/
├── README.md
├── requirements.txt
├── requirements-dev.txt
├── .gitignore
├── conftest.py
├── config.yaml                  # サンプル設定（sequence/mzml_pathは空、埋めてから実行）
├── docs/design/
│   ├── RNA_MassHunter_再設計_実装仕様書.md   # 詳細仕様（必読）
│   ├── rna-masshunter-design-memo.md          # ChatGPTによる原設計メモ
│   └── Claude_Code_実装引き継ぎ指示書.md      # このファイル
├── rna_masshunter/
│   ├── __init__.py
│   ├── masses.py                # 移植そのまま
│   ├── elemental_composition.py # 移植そのまま
│   ├── enzymes.py                # 移植そのまま
│   ├── digestion.py             # 移植そのまま
│   ├── cca_processing.py        # 移植+編集（Enum inline化）
│   ├── mzml_reader.py           # 移植そのまま
│   ├── mzml_diagnostics.py      # 移植そのまま
│   ├── peak_picking.py          # 移植そのまま
│   ├── modifications.py         # 移植そのまま
│   ├── mass_shift_ms1_search.py # 移植そのまま
│   ├── ms1_mapping.py           # 移植そのまま
│   ├── ms2_extraction.py        # 移植+編集（不要依存の除去）新規ファイル
│   ├── warnings_manager.py      # 移植そのまま
│   ├── models.py                # 新規作成（必要フィールドのみ抜粋）
│   └── config.py                # 新規作成（必要セクションのみ）
├── data/
│   ├── modifications.yaml       # 移植そのまま（118件）
│   └── base_masses.yaml         # 移植そのまま
└── tests/
    ├── test_elemental_composition.py   # 移植そのまま
    └── test_vendored_smoke.py          # 新規（移植の動作確認用）
```

---

## 4. ここから実装すること（設計仕様書§6, §22の要約）

設計仕様書の「§8 クラス/関数仕様」に各モジュールの詳細な関数シグネチャ・処理内容が書かれています。以下は実装順序の要約です（詳細と完了条件は設計仕様書§22の表を参照）。

| Phase | 追加するファイル | 内容 |
|---|---|---|
| 1 | `rna_masshunter/formula_candidate.py` | ΔDaから元素組成候補を列挙するDFS探索（§8.1, §12） |
| 2 | `rna_masshunter/mass_comparison.py`, `rna_masshunter/observed_mass.py` | Fragment×charge×PeakマッチングとMassComparisonRow生成（§8.2, §8.3, §15） |
| 2.5 | `rna_masshunter/ms2_support.py` | MS2参考情報の付加（§8.5, §14A） |
| 3 | `rna_masshunter/simple_pipeline.py`, `simple_main.py` | 全体オーケストレーション（§8.4） |
| 4 | `rna_masshunter/excel_report.py`（新規ファイル、元の同名ファイルとは別物） | 01〜08シートのExcel出力（§16） |
| 5 | `tests/test_formula_candidate.py`, `tests/test_mass_comparison.py`, `tests/test_ms2_support.py`, `tests/test_simple_pipeline.py` | 各フェーズのテスト（§19, §20） |
| 6 | — | 実mzMLによる動作確認（§21、ユーザー立ち会い） |

**各フェーズはそれぞれ独立してコミットし、フェーズ完了時に対応するテストが通っていることを確認してから次のフェーズに進んでください。** 既存の `tests/test_vendored_smoke.py` を壊していないことも都度確認してください（`pytest` を実行するだけです）。

---

## 5. 遵守事項（設計仕様書§23と同一、重要なので再掲）

1. 既存の巨大な `nyako0813/RNA_MassHunter` リポジトリには一切手を入れない。参照専用（vendoring元の確認）以外でcloneする必要すらない。
2. Formula CandidateとModification Candidateのロジック・出力列を完全に分離する。
3. 元素組成候補のランキングで原子数を順位付けキーに使わない（3原子以上でも順位を下げない）。テストで担保する。
4. MS2関連の実装は「参考情報の付加」の範囲に厳密に留める。自動同定・スコアリング・修飾位置推定は行わない。
5. 内部計算はPythonの`float`（倍精度）を維持し、丸めはExcel表示側でのみ行う。
6. 各フェーズごとに小さくコミットし、テストが通った状態を保つ。テストが落ちている状態でコミットしない。
7. コミットメッセージは日本語の簡潔な要約とする（元リポジトリの慣例に合わせる）。
8. 実装中に設計仕様書との齟齬や未確定事項（仕様書末尾の「未確定事項・要確認リスト」を参照）に気づいた場合は、**独断で仕様を変更せず、いったん立ち止まってユーザー（nyakoさん）に確認する。**

---

## 6. 完了の定義（MVP）

以下が揃った時点でこの引き継ぎのゴールとします。

- `config.yaml` に実際のtRNA配列とmzMLパスを設定し、`python simple_main.py --config config.yaml` が例外なく完走する。
- 出力Excel（`output/RNA_MassHunter_simple_report.xlsx`）に `01_Index` 〜 `07_Modifications` の全シートが存在し、`01_Index` の全リンクが機能する。
- `06_Mass_Comparison` に、ΔDa/Δppm・Recommended Formula・Formula Candidates・Known Modification・Modification Candidatesの列が正しく埋まっている。
- `formula_candidate.py` の原子数バイアス無しランキング（原則3項目）がテストで明示的に検証されている。
- ユーザーの実データ（実mzML）で一度動作確認が完了している（設計仕様書§21）。
- `08_Visualization` は実装できていればなお良いが、MVPの必須条件ではない（設計仕様書§27の通り）。

---

## 7. 困ったとき

- 設計仕様書の内容と実際のコードが食い違う場合は、**実際のコードを優先しつつ、仕様書側の記述も更新して矛盾を残さない**（仕様書は「生きたドキュメント」として扱う）。
- `rna_masshunter/models.py` や `config.py` に新しいフィールドが必要になった場合、元の `nyako0813/RNA_MassHunter` リポジトリの該当箇所（`models.py`, `config.py`）を参照専用でcloneして確認して構いませんが、**このリポジトリにその依存を追加しない**（必要な定義だけをこのリポジトリにコピーする、これまでと同じ方針）。
- 大きな設計判断（新しいシートを増やす、既存の列構成を変える等）が必要になった場合はユーザーに確認してください。
