# Variable and Graph Explanation

This file explains the main variables in the screening CSV files and the meaning of each generated graph.

この結果は semiconductor gas candidate screening MVP の探索結果です。EHS、法規制、SDS、工程安全の確定判定ではありません。候補比較とレビュー優先度付けのためのスクリーニング結果として扱ってください。

## Output Files

| File | Meaning |
|---|---|
| `summary.csv` | 全候補の最終サマリー。採用された物性値、判定値、候補由来、レビュー要否を含む。 |
| `amine_summary.csv` | アミン系候補だけを抽出したサマリー。アミン分類、置換基、半導体用途タグを含む。 |
| `amine_class_summary.csv` | アミン分類ごとの集計。1級、2級、3級、環状、金属アミドなどの件数を比較する。 |
| `candidate_breakdown.csv` | 候補ファミリー、生成ルール、品質、PFAS判定などの内訳集計。 |
| `evidence.csv` | 物性値や判定値の根拠を保存するシート。推定値、ローカル表、候補値などの生データ確認用。 |
| `coverage.csv` | 各候補でどの物性や判定が取得できたかのカバレッジ確認用。 |
| `review_required.csv` | 手動レビューが必要な候補の一覧。 |
| `rejected.csv` | prefilter で落ちた候補の一覧。 |
| `planned_api.csv` | API取得予定候補。現在の実行では remote 無効のため基本的に空。 |
| `run_stats.csv` | 生成数、重複排除数、prefilter通過数、API要求数などの実行統計。 |

## Key Variable Groups

### Identity Variables

| Variable | Meaning |
|---|---|
| `input_name` | 入力名または生成時の候補名。 |
| `preferred_name` | 表示用に採用した候補名。 |
| `cas` | CAS番号。主キーではなく参考情報。 |
| `pubchem_cid` | PubChem CID。remote無効時は空が多い。 |
| `formula` | 分子式。 |
| `canonical_smiles` | 正規化されたSMILES。構造同定に使う。 |
| `inchikey` | InChIKey。利用可能な場合の重複排除や構造同定キー。 |
| `identity_status` | 同定状態。`resolved` は同定済み、`manual_review_required` は確認が必要。 |
| `candidate_id` | パイプライン内での候補ID。 |

### Basic Molecular Variables

| Variable | Unit | Meaning |
|---|---:|---|
| `molecular_weight` | g/mol | 分子量。候補の揮発性、供給形態、分子サイズの目安。 |
| `heavy_atom_count` | count | H以外の原子数。候補サイズや探索制限に使う。 |
| `element_symbols` | text | 含有元素の一覧。例: `C;N;F`。 |
| `structure_status` | text | 構造の扱い。`valid` は処理可能、formula-onlyや未対応構造は注意。 |
| `structure_status_reason` | text | 構造状態の理由。RDKitなしの場合の軽量正規化理由などを含む。 |

### Temperature and Pressure Variables

| Variable | Unit | Meaning |
|---|---:|---|
| `tm_C` | degC | normal melting point。融点。25Cで固体かどうかの判断に使う。 |
| `tb_C` | degC | normal boiling point。沸点。常温で気体/液体か、供給しやすさの判断に使う。 |
| `tc_C` | degC | critical temperature。臨界温度。蒸気圧計算や超臨界的挙動の判定に使う。 |
| `pc_MPa` | MPa | critical pressure。臨界圧力。流体物性の目安。 |
| `pvap_25C_kPa` | kPa | 25Cでの蒸気圧。室温供給性の重要指標。 |
| `pvap_40C_kPa` | kPa | 40Cでの蒸気圧。軽加温時の供給性の指標。 |
| `pvap_60C_kPa` | kPa | 60Cでの蒸気圧。加温ソースでの供給性の指標。 |
| `pvap_25C_status` | text | 25C蒸気圧の状態。値あり、欠損、臨界温度超過など。 |
| `pvap_40C_status` | text | 40C蒸気圧の状態。 |
| `pvap_60C_status` | text | 60C蒸気圧の状態。 |

### Phase and Supply Variables

| Variable | Meaning |
|---|---|
| `phase_25C_1atm` | 25C、1 atmでの相推定。`gas`, `liquid`, `solid`, `gas_or_supercritical`, `unknown`。 |
| `supply_class` | 半導体装置向け供給形態の目安。圧縮/液化ガス、液体バブラー、加温ソース、固体昇華レビューなど。 |

### Environmental Screening Variables

| Variable | Meaning |
|---|---|
| `gwp100_ar6` | IPCC AR6系の100年GWP値。単位は kg CO2e / kg。ローカル表に一致した候補のみ値が入る。 |
| `gwp100_ar6_status` | GWP値の採用状態。`selected` は採用、`missing` はローカル表に値なし。 |
| `gwp100_ar6_source` | GWP値の由来。例: `GWP_local`。 |
| `pfas_flag` | PFASスクリーニング結果。`yes`, `possible`, `no`, `unknown`。確定判定ではなくレビュー用フラグ。 |
| `pfas_basis` | PFAS判定の根拠。構造ルール、ローカルリストなど。 |
| `pfas_list_hits` | ローカルPFASリストで一致した場合のヒット情報。 |
| `persistence_screen` | 残留性の簡易スクリーニング。`likely_persistent`, `unknown` など。 |
| `persistence_basis` | 残留性判定の根拠。 |

### Reactivity and Coverage Variables

| Variable | Meaning |
|---|---|
| `reactive_groups` | 反応性グループ。例: amine/basic nitrogen, silane/reducing gas。 |
| `reactivity_flags` | レビュー用反応性フラグ。例: `acid_base_reactive`, `combustible_or_flammable_review`。 |
| `ghs_physical_h_codes` | ローカル表にある物理危険性Hコード。確定GHS分類ではない。 |
| `kinetics_coverage` | kinetics情報のカバレッジ。`available`, `partial`, `not_checked` など。 |
| `kinetics_sources` | kinetics候補情報の由来。 |
| `data_quality` | 選択値の品質ランク。`A/B/C/D`, `Partial`, `Conflict`, `Missing`。 |
| `review_required` | 手動レビューが必要かどうか。未解決同定、PFAS可能性、品質不足などでTrueになる。 |

### Generation Variables

| Variable | Meaning |
|---|---|
| `candidate_source` | 候補の由来。例: `seed`, `template_generation`, `semiconductor_amine`, `formula_generation`。 |
| `candidate_family` | 候補ファミリー。例: `amine`, `fluorocarbon`, `inorganic`。 |
| `generation_rule` | 生成ルール。例: `primary_amine`, `tertiary_amine`, `curated_fluorocarbon_formula`。 |
| `parent_candidate_id` | LocalMutationなど親候補がある場合の親ID。 |
| `candidate_scope` | 候補の探索スコープ。例: curated semiconductor amine/process gas。 |

### Amine-Specific Variables

| Variable | Meaning |
|---|---|
| `amine_class` | アミン分類の内部ラベル。例: `primary`, `secondary`, `tertiary`, `cyclic`。 |
| `amine_class_label` | 表記ゆれを避ける英語ラベル。例: `primary_amine`, `tertiary_amine`。 |
| `amine_detail` | フッ素化なども含めた詳細ラベル。 |
| `amine_substituents` | N原子まわりの置換基。例: `methyl; ethyl`。 |
| `amine_substituent_profile` | 置換基の種類。例: `alkyl`, `fluorinated_alkyl`, `unsaturated`, `cyclic_ring`。 |
| `amine_substituent_count` | 置換基数。1級/2級/3級の理解に使う。 |
| `fluorinated_amine` | フッ素を含むアミン候補かどうか。 |
| `fluorinated_substituent_count` | フッ素化置換基の数。 |
| `amine_fluorination_level` | フッ素化の程度。`none`, `single_fluorinated_substituent`, `multiple_fluorinated_substituents`。 |
| `unsaturated_amine` | vinyl/allylなど不飽和置換基を含むか。 |
| `cyclic_amine_or_substituent` | 環状アミンまたは環状置換基を含むか。 |
| `ring_name` | 環状候補の場合のリング名。 |
| `precursor_family` | 半導体前駆体としてのファミリー。例: aminosilane, metal amideなど。 |
| `semiconductor_process_roles` | 半導体プロセスで想定される役割タグ。例: ALD/CVD precursor, etch/clean candidate。 |
| `semiconductor_relevance_basis` | 半導体用途タグの根拠。 |

## Run Stats Meaning

| Metric | Meaning |
|---|---|
| `generated_candidates` | generatorが出した候補数。 |
| `normalized_candidates` | 正規化処理後の候補数。 |
| `deduplicated_candidates` | 重複排除後の候補数。 |
| `prefilter_passed` | prefilterを通過した候補数。 |
| `prefilter_rejected` | prefilterで除外された候補数。 |
| `pubchem_requests` | PubChem PUG-RESTへの実リクエスト数。 |
| `pugview_requests` | PubChem PUG-Viewへの実リクエスト数。 |
| `cache_hits` | API cacheで回答できた件数。 |
| `negative_cache_hits` | negative cacheでAPI不要と判断できた件数。 |
| `planned_api_requests` | dry-runまたはplanner上のAPI予定件数。 |
| `planned_api_candidates` | API取得予定候補数。 |
| `remote_enrichment_candidates` | remote enrichment対象候補数。 |
| `final_summary_rows` | Summaryに残った最終候補数。 |

## Graph Reading Guide

### General Rules

| Visual Element | Meaning |
|---|---|
| Point color by `pfas_flag` | 緑: `pfas_no`, 黄: `pfas_possible`, 赤: `pfas_yes`。PFAS懸念の有無を比較する。 |
| Point color by `amine_class_label` | 1級/2級/3級/環状など、アミン分類ごとの差を見る。 |
| Log scale | 蒸気圧やGWPのように桁差が大きい値を見やすくするための対数軸。距離は比率を表す。 |
| Boxplot | 中央値、四分位、分布幅を見る。外れ値は抑えて全体傾向を見やすくしている。 |
| Heatmap | 行と列の組み合わせ件数を色の濃さで示す。 |

## Graph-by-Graph Explanation

| No. | File | Meaning |
|---:|---|---|
| 01 | `01_amine_class_fluorination.png` | アミン分類ごとに、非フッ素化候補とフッ素化候補の数を比較する。どのクラスでフッ素化候補が多いかを見る。 |
| 02 | `02_amine_class_pfas.png` | アミン分類ごとのPFAS判定内訳。`possible` や `yes` が多いクラスはレビュー優先度が高い。 |
| 03 | `03_property_coverage.png` | アミン候補で各物性値がどれだけ埋まっているかを見る。データ欠損の把握用。 |
| 04 | `04_boiling_point_by_class.png` | アミン分類ごとの沸点分布。低沸点側はガス/高揮発、 高沸点側は液体/加温供給寄り。 |
| 05 | `05_vapor_pressure_25C_by_class.png` | 25C蒸気圧の分類別分布。高いほど室温供給しやすい。縦軸はlog。 |
| 06 | `06_phase_supply_matrix.png` | 25C相と供給クラスの組み合わせ。装置供給の観点で候補群を俯瞰する。 |
| 07 | `07_mw_vs_boiling_pfas.png` | MWと沸点の関係。一般にMWが上がると沸点も上がりやすい。色でPFAS懸念も同時確認する。 |
| 08 | `08_precursor_family_counts.png` | curated半導体アミン前駆体ファミリーの件数。ALD/CVD関連候補の偏りを見る。 |
| 09 | `09_fluorination_pfas_matrix.png` | フッ素化レベルとPFAS判定の関係。フッ素化置換基がPFAS懸念にどうつながるかを見る。 |
| 10 | `10_process_role_counts.png` | 半導体プロセス用途タグの件数。候補がどの用途に多いかを見る。 |
| 11 | `11_property_scatter_matrix_by_class.png` | MW、融点、沸点、臨界温度、臨界圧、蒸気圧を総当たりで比較する。色はアミン分類。 |
| 12 | `12_tb_vs_pvap25_by_class.png` | 沸点と25C蒸気圧の関係。沸点が低い候補ほど蒸気圧が高くなりやすい。 |
| 13 | `13_mw_vs_pvap25_by_class.png` | MWと25C蒸気圧の関係。分類ごとの揮発性の違いを見る。 |
| 14 | `14_tc_vs_pc_by_class.png` | 臨界温度と臨界圧の関係。流体物性のまとまりや外れ候補を見る。 |
| 15 | `15_tm_vs_tb_by_class.png` | 融点と沸点の関係。25C近辺で固体/液体/気体になりやすい候補を把握する。 |
| 16 | `16_pvap25_vs_pvap60_by_class.png` | 25C蒸気圧と60C蒸気圧の関係。加温でどれだけ供給性が上がるかを見る。 |
| 17 | `17_gwp_availability_by_family.png` | 候補ファミリーごとのGWP値取得率。現状はフッ素系ガスにGWP値が多い。 |
| 18 | `18_gwp100_ar6_known_values.png` | GWP既知値ランキング。SF6、NF3、PFC/HFC系など高GWP候補をすぐ確認する。 |
| 19 | `19_gwp_by_pfas_flag.png` | PFAS判定別にGWPを表示する。PFAS構造判定と温暖化影響が必ず一致するわけではない点を見る。 |
| 20 | `20_gwp_vs_molecular_weight.png` | GWPとMWの関係。GWPはMWだけでは決まらないため、外れ候補の確認に使う。 |
| 21 | `21_mw_vs_melting_pfas.png` | MWと融点の関係。重い候補やPFAS懸念候補が固体化しやすいかを見る。 |
| 22 | `22_mw_vs_critical_temperature_pfas.png` | MWと臨界温度の関係。候補サイズと流体特性の傾向を見る。 |
| 23 | `23_mw_vs_critical_pressure_pfas.png` | MWと臨界圧の関係。高圧条件や物性モデルの外れ候補を見る。 |
| 24 | `24_mw_vs_pvap25_pfas.png` | MWと25C蒸気圧の関係。室温供給しやすい低MW候補と、低蒸気圧候補を比較する。 |
| 25 | `25_mw_vs_pvap40_pfas.png` | MWと40C蒸気圧の関係。軽加温で供給可能性が上がる候補を見る。 |
| 26 | `26_mw_vs_pvap60_pfas.png` | MWと60C蒸気圧の関係。加温ソースでの供給候補を比較する。 |

## Practical Interpretation for Semiconductor Gas Screening

| Question | Useful Variables / Graphs |
|---|---|
| 室温ガス候補を探したい | `phase_25C_1atm`, `pvap_25C_kPa`, graphs 05, 06, 24 |
| 加温供給で使えそうな候補を探したい | `pvap_40C_kPa`, `pvap_60C_kPa`, `supply_class`, graphs 16, 25, 26 |
| アミン分類ごとの候補分布を見たい | `amine_class_label`, `amine_substituent_profile`, graphs 01, 02, 04, 05 |
| PFAS懸念を確認したい | `pfas_flag`, `pfas_basis`, graphs 02, 07, 09, 19, 21-26 |
| GWPが高い候補を確認したい | `gwp100_ar6`, `gwp100_ar6_status`, graphs 17-20 |
| レビューが必要な候補を抽出したい | `review_required`, `data_quality`, `identity_status`, `structure_status` |
| 半導体前駆体として有用そうなアミンを見たい | `precursor_family`, `semiconductor_process_roles`, graphs 08, 10 |

## Current Dataset Notes

- `run_stats.csv` shows `pubchem_requests = 0` and `pugview_requests = 0`, so the current output was generated without live remote API calls.
- `final_summary_rows = 687` and `amine_summary.csv` contains 672 amine candidates.
- GWP values are only filled when the candidate matched the local `examples/gwp_sample.csv` table. Missing GWP is not the same as low GWP.
- Many amine thermodynamic values are local estimates or local/curated values. Check `evidence.csv` before using a value for engineering design.
