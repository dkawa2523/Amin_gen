# Semiconductor Gas Candidate Screening MVP

半導体プロセス向けガス候補を、候補生成、構造正規化、物性取得、PFAS/残留性/反応性スクリーニング、反応速度データ coverage、Excel 出力まで一通り流す Python MVP です。

この実装は初期スクリーニング用です。法規制判断、GHS 承認、装置適合性、プロセス保証の最終判断には使わず、Evidence と Review Required を専門家レビューの入口として扱ってください。

## できること

- 入力 CSV または内蔵 seed から候補分子を生成
- SMILES 正規化、重複除去、構造 prefilter
- `Tm`, `Tb`, `Tc`, `Pc`, `Pvap(25/40/60C)`, `GWP100_AR6` の候補値収集と採用
- `phase_25C_1atm` と `supply_class` の導出
- PFAS 構造ルール、PFAS リスト照合、残留性 screening
- ローカル表と PUG-View による反応性/GHS physical H-code 整理
- HF/O/F/OH/Cl/e- などの反応速度データ coverage 表示
- `Summary`, `Evidence`, `Coverage`, `Review Required`, `Rejected` の Excel 出力

## セットアップ

```bash
cd gas_screening_mvp_project
py -m venv .venv
.venv\Scripts\activate
py -m pip install -e ".[dev]"
```

macOS/Linux の場合:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

熱物性・構造正規化を強化する場合は optional dependencies を入れます。

```bash
py -m pip install -e ".[chem,dev]"
```

RDKit が無い環境でも軽量 SMILES parser でサンプル pipeline は動きます。ただし、本番の構造正規化、InChIKey 生成、複雑な SMILES 判定には RDKit を使ってください。

## 実行例

```bash
py -m gas_screening_mvp.cli run ^
  --config examples/demo_config.yml ^
  --input examples/sample_input.csv ^
  --output outputs/demo_screening.xlsx
```

インストール後は console script でも実行できます。

```bash
gas-screening run --config examples/demo_config.yml --input examples/sample_input.csv --output outputs/demo_screening.xlsx
```

### Run modes

- `enrichment`: デフォルト。ローカル生成とローカル表を中心に、安全側で remote は config/CLI で明示した場合のみ使います。
- `exploration`: 候補探索用。formula-only、LocalMutation、PubChem similarity expansion などの探索機能を許可します。ただし各 generator/API は個別 enable が必要です。
- `refresh`: 既存候補の再取得・再評価用のモード名です。remote 利用可否は `--remote` / `--no-remote` または config に従います。

CLI override:

```bash
py -m gas_screening_mvp.cli run ^
  --config examples/exploration_config.yml ^
  --input examples/sample_input.csv ^
  --output outputs/exploration_screening.xlsx ^
  --mode exploration ^
  --dry-run
```

`--dry-run` は外部 API を呼ばず、生成、正規化、重複排除、prefilter、予定 API 件数を `Run Stats` と `Planned API` sheet に出します。

`--remote` は PubChem/PUG-View を有効化し、`--no-remote` は無効化します。デフォルト config では remote は無効です。

## 入力 CSV

最小列:

```csv
name,smiles,cas,family
triethylamine,CCN(CC)CC,121-44-8,amine
carbon tetrafluoride,FC(F)(F)F,75-73-0,fluorocarbon
```

## ローカル表

`examples/demo_config.yml` から以下の CSV を指定できます。

- `curated_properties_csv`: 手動確認済み物性値
- `gwp_csv`: IPCC/EPA などから作成した GWP 表
- `pfas_list_csv`: PFAS リストの InChIKey 表
- `reactivity_csv`: CAMEO/GHS/社内ルール由来の反応性表
- `kinetics_csv`: NIST/LXCat/QDB などの coverage 表

`examples/*_sample.csv` は pipeline 確認用の小さな例です。本番では公式データまたは社内確認済みデータに置き換えてください。

## 外部 API

デフォルトでは外部 API は無効です。PubChem を使う場合:

```yaml
providers:
  pubchem_enabled: true
  pugview_enabled: true
  pubchem_rps: 3.0
  pugview_rps: 1.0

fetch_planning:
  min_identity_api_score: 0.45
  min_enrichment_api_score: 0.85
```

API 呼び出しは SQLite cache と negative cache を使います。PUG-View は shortlist enrichment として扱い、高頻度の物性抽出には使いません。

## テスト

```bash
py -m pytest -q
```

## コード構成

```text
src/gas_screening_mvp/
  domain/          dataclass model と API request signature
  generation/      seed/template/formula 候補生成
  normalization/   RDKit 正規化と軽量 SMILES fallback
  prefilter/       構造 filter と API priority score
  providers/       local/remote DB provider
  planning/        fetch planner と rate limiter
  selection/       単位変換、採用値選択、競合検出
  classification/  PFAS rule
  derivation/      phase/supply_class/summary
  storage/         SQLite API cache
  export/          Excel exporter
```
