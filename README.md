# tWhEEL - トランスクリプトーム機能解析パイプライン

**Transcriptome Workflow for Hierarchical Exploration and Expression Levels**

tWhEELは、トランスクリプトアセンブリデータから遺伝子予測、機能アノテーション、CRISPRガイド設計、GO機能解析を行い、ノックアウト実験のためのターゲット遺伝子を特定するための包括的なパイプラインです。

## 特徴

- **遺伝子予測**: TransDecoderを使用したORF予測とBUSCO品質評価
- **機能アノテーション**: UniProtデータベースに対するFASTA検索と遺伝子分類
- **CRISPRガイド設計**: CrisFlashを使用したSpCas9ガイドRNA設計
- **GO機能解析**: GO Knock-out Feasibility Index (GO-KF)の計算
- **インタラクティブレポート**: HTMLベースのウェブレポート生成
- **柔軟な設定**: YAMLベースの設定ファイル

## 必要な環境

- Python 3.8以上
- Docker
- 必要なPythonパッケージ（requirements.txtを参照）

## インストール

```bash
# リポジトリをクローン
git clone https://github.com/KazukiNakamae/tWhEEL.git
cd tWhEEL

# 必要なパッケージをインストール
pip install -r requirements.txt

# 必要に応じてDockerイメージをビルド
docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/
docker build -t kazukinakamae/compute_go_kfi:v1.0 src/docker/compute_go_kfi/
docker build -t kazukinakamae/count_go_parent_terms:v1.0 src/docker/count_go_parent_terms/
```

## 使用方法

### 1. 設定ファイルの準備

`config/default.yaml`をコピーして、入力ファイルのパスを設定します：

```yaml
input:
  genome_assembly: "/path/to/genome.fna"
  transcriptome: "/path/to/transcriptome.fsa_nt"
  uniprot_db: "/path/to/uniprot_sprot.fasta"
```

### 2. 実行

```bash
# 全ステップを実行
python src/twheel.py --config config/my_config.yaml

# 特定のステップのみ実行
python src/twheel.py --config config/my_config.yaml --step annotation

# 利用可能なステップを確認
python src/twheel.py --list-steps
```

### 3. 結果の確認

実行後、以下のファイルが生成されます：

- `outputs/summary.txt`: テキストサマリー
- `outputs/report.html`: HTMLインタラクティブレポート
- `outputs/web_data/`: Webアプリケーション用データ
- `outputs/intermediate/`: 中間ファイル

## パイプライン構成

### 1. 遺伝子予測 (gene_prediction)
- **TransDecoder**: ORF予測とタンパク質配列生成
- **BUSCO**: アセンブリ品質評価
- **seqkit**: 配列統計情報の取得

### 2. アノテーション (annotation)
- **FASTA検索**: UniProtデータベースに対するggsearch36
- **遺伝子分類**: Multi/Single遺伝子の分類
- **配列抽出**: 分類された遺伝子配列の抽出

### 3. ガイド設計 (guide_design)
- **CrisFlash**: SpCas9用ガイドRNA設計
- **サイト特異性**: 特異的標的の抽出とフィルタリング
- **アノテーション統合**: ガイド情報とアノテーションの結合

### 4. GO解析 (go_analysis)
- **GOアノテーション**: UniProt IDからのGO term取得
- **GO-KFI計算**: ノックアウト実行可能性指数の計算
- **機能分類**: Fully/Partially/Not targetableの分類

### 5. レポート生成 (report_generation)
- **テキストサマリー**: 実行統計とサマリー
- **HTMLレポート**: インタラクティブなウェブレポート
- **データベース**: SQLiteとJSON形式でのデータ保存

## 出力ファイル

### 主要出力
- `summary.txt`: 実行サマリー
- `report.html`: HTMLレポート
- `web_data/twheel_data.db`: SQLiteデータベース
- `web_data/summary.json`: JSON形式の統計情報

### 中間ファイル
- `intermediate/query.pep`: 予測タンパク質配列
- `intermediate/annotation_classified.txt`: 分類されたアノテーション
- `intermediate/annotation_with_guides_final.tsv`: ガイド情報統合結果
- `intermediate/go_kfi_table.csv`: GO-KFIテーブル
- `intermediate/go_parent_counts.png`: GO階層可視化

## 設定オプション

### 入力設定
```yaml
input:
  genome_assembly: ""  # ゲノムアセンブリファイル
  transcriptome: ""    # トランスクリプトームファイル
  uniprot_db: ""       # UniProtデータベース
```

### CRISPR設定
```yaml
crispr:
  pam_pattern: "NGG"   # PAMパターン
  max_mismatch: 2      # 最大ミスマッチ数
  threads: 10          # スレッド数
```

### GO解析設定
```yaml
go_analysis:
  add_go_image: "kazukinakamae/add-go-columns:v1.0"
  compute_kfi_image: "kazukinakamae/compute_go_kfi:v1.0"
  count_parent_image: "kazukinakamae/count_go_parent_terms:v1.0"
```

## HTMLレポート機能

生成されるHTMLレポートには以下の機能が含まれます：

### ダッシュボード
- 処理統計のサマリー
- GO Knock-out Feasibility分布グラフ

### 遺伝子ブラウザ
- Multi/Single遺伝子の検索・フィルタリング
- UniProt ID、同一性、E値の表示
- ガイドRNA数の表示

### ガイドRNA
- 設計されたガイドRNAの一覧
- 染色体位置、配列、ストランド情報
- 検索・フィルタリング機能

### GO解析
- GO termの検索・フィルタリング
- KFIスコアによる分類
- ターゲタブル分類の表示

## トラブルシューティング

### Docker関連の問題
```bash
# Dockerイメージが見つからない場合
docker build -t kazukinakamae/add-go-columns:v1.0 src/docker/add_go_columns/

# 権限の問題
sudo usermod -aG docker $USER
```

### メモリ不足の問題
```yaml
# 設定ファイルでメモリ使用量を調整
crispr:
  memory: "20g"  # メモリ使用量を削減
```

### パフォーマンス調整
```yaml
# CPU使用量の調整
crispr:
  threads: 4     # スレッド数を削減
busco:
  cpu: 4         # CPU数を削減
```

## 参考文献

- **TransDecoder**: Haas et al. (2013) Nature Protocols
- **BUSCO**: Simão et al. (2015) Bioinformatics
- **CrisFlash**: Nakamae et al. (2024) *準備中*
- **UniProt**: UniProt Consortium (2023) Nucleic Acids Research

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。詳細は[LICENSE](LICENSE)ファイルを参照してください。

## 貢献

プロジェクトへの貢献を歓迎します。バグ報告や機能要求はIssuesページから、コードの改善はPull Requestからお願いします。

## 作成者

- **Kazuki Nakamae** - [GitHub](https://github.com/KazukiNakamae)

## 更新履歴

- **v1.0.0** (2025-01-18): 初回リリース
  - 基本パイプライン機能
  - HTMLレポート生成
  - GO-KFI計算
