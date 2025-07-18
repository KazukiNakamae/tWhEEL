#!/usr/bin/env python3
"""
GO解析モジュール
GO Knock-out Feasibility Index (GO-KF)の計算とGO機能解析
"""

import os
import subprocess
import logging
import csv
from pathlib import Path


class GOAnalysis:
    """GO解析処理クラス"""
    
    def __init__(self, config):
        """
        初期化
        
        Args:
            config (dict): 設定辞書
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.docker_dir = os.path.join('src', 'docker')
        
    def build_docker_images(self):
        """カスタムDockerイメージのビルド"""
        self.logger.info("カスタムDockerイメージをビルド中...")
        
        # 必要なDockerイメージのリスト
        docker_images = [
            ('add_go_columns', self.config['go_analysis']['add_go_image']),
            ('compute_go_kfi', self.config['go_analysis']['compute_kfi_image']),
            ('count_go_parent_terms', self.config['go_analysis']['count_parent_image'])
        ]
        
        for folder_name, image_name in docker_images:
            docker_path = os.path.join(self.docker_dir, folder_name)
            
            if os.path.exists(docker_path):
                self.logger.info(f"Dockerイメージをビルド: {image_name}")
                
                cmd = [
                    "docker", "build", 
                    "-t", image_name,
                    docker_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    self.logger.error(f"Dockerビルドエラー ({image_name}): {result.stderr}")
                    raise RuntimeError(f"Dockerイメージのビルドに失敗しました: {image_name}")
                
                self.logger.info(f"Dockerイメージビルド完了: {image_name}")
            else:
                self.logger.warning(f"Dockerディレクトリが見つかりません: {docker_path}")
        
        self.logger.info("カスタムDockerイメージのビルドが完了しました")
    
    def add_go_annotations(self, input_path):
        """GOアノテーションの追加"""
        self.logger.info("GOアノテーションを追加中...")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'annotation_with_guides_go.tsv')
        
        # DockerでGOアノテーション追加
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.getcwd()}:/DATA",
            self.config['go_analysis']['add_go_image'],
            f"/DATA/{input_path}",
            f"/DATA/{output_path}"
        ]
        
        self.logger.debug(f"GOアノテーション追加コマンド: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GOアノテーション追加エラー: {result.stderr}")
            raise RuntimeError("GOアノテーション追加に失敗しました")
        
        self.logger.info(f"GOアノテーション追加完了: {output_path}")
        return output_path
    
    def compute_go_kfi(self, input_path):
        """GO Knock-out Feasibility Index (GO-KF)の計算"""
        self.logger.info("GO-KFIを計算中...")
        
        # 出力ファイルパス
        kfi_table_path = os.path.join(self.intermediate_dir, 'go_kfi_table.csv')
        kfi_stat_path = os.path.join(self.intermediate_dir, 'go_kfi_stat.csv')
        
        # DockerでGO-KFI計算
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.getcwd()}:/DATA",
            self.config['go_analysis']['compute_kfi_image'],
            "-i", f"/DATA/{input_path}",
            "-o", f"/DATA/{kfi_table_path}",
            "-s", f"/DATA/{kfi_stat_path}"
        ]
        
        self.logger.debug(f"GO-KFI計算コマンド: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GO-KFI計算エラー: {result.stderr}")
            raise RuntimeError("GO-KFI計算に失敗しました")
        
        # 結果をログに出力
        self.logger.info("GO-KFI計算結果:")
        self.logger.info(result.stdout)
        
        self.logger.info(f"GO-KFI計算完了: {kfi_table_path}, {kfi_stat_path}")
        return kfi_table_path, kfi_stat_path
    
    def extract_fully_targetable_gos(self, kfi_table_path):
        """Fully targetable GOsの抽出"""
        self.logger.info("Fully targetable GOsを抽出中...")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'go_kfi_table_fully_targetable.csv')
        
        # CSV読み込みと処理
        with open(kfi_table_path, 'r') as f_in:
            with open(output_path, 'w') as f_out:
                reader = csv.reader(f_in)
                writer = csv.writer(f_out)
                
                # ヘッダーを書き込み
                header = next(reader)
                writer.writerow(header)
                
                # データを処理
                for row in reader:
                    if len(row) >= 5:  # 4列目（インデックス4）をチェック
                        try:
                            # 4列目が1（Fully targetable）の行のみを抽出
                            if float(row[4]) == 1.0:
                                writer.writerow(row)
                        except (ValueError, IndexError):
                            continue
        
        # 抽出されたGO数を数える
        with open(output_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # ヘッダーをスキップ
            fully_targetable_count = sum(1 for _ in reader)
        
        self.logger.info(f"Fully targetable GOs数: {fully_targetable_count}")
        self.logger.info(f"Fully targetable GOs抽出完了: {output_path}")
        
        return output_path
    
    def count_go_parent_terms(self, fully_targetable_path):
        """GO親termの計算とvisualization"""
        self.logger.info("GO親termの計算とvisualizationを実行中...")
        
        # 出力ファイルパス
        parent_csv_path = os.path.join(self.intermediate_dir, 'go_parent_counts.csv')
        parent_png_path = os.path.join(self.intermediate_dir, 'go_parent_counts.png')
        
        # DockerでGO親term分析
        cmd = [
            "docker", "run", "-it", "--rm",
            "-v", f"{os.getcwd()}:/DATA",
            self.config['go_analysis']['count_parent_image'],
            "--go_table", f"/DATA/{fully_targetable_path}",
            "--go_col", "0",
            "--treemap",
            "--out_csv", f"/DATA/{parent_csv_path}",
            "--out_png", f"/DATA/{parent_png_path}"
        ]
        
        self.logger.debug(f"GO親term分析コマンド: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GO親term分析エラー: {result.stderr}")
            raise RuntimeError("GO親term分析に失敗しました")
        
        self.logger.info(f"GO親term分析完了: {parent_csv_path}, {parent_png_path}")
        return parent_csv_path, parent_png_path
    
    def parse_go_kfi_statistics(self, kfi_stat_path):
        """GO-KFI統計の解析"""
        self.logger.info("GO-KFI統計を解析中...")
        
        stats = {}
        
        try:
            with open(kfi_stat_path, 'r') as f:
                content = f.read()
                
                # 統計情報を抽出
                import re
                
                # 基本統計を抽出
                pattern = r'F:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,M:(\d+\.?\d*)%\],P:(\d+\.?\d*)%,N:(\d+\.?\d*)%,n:(\d+)'
                match = re.search(pattern, content)
                
                if match:
                    stats['fully_targetable_percent'] = float(match.group(1))
                    stats['simply_targetable_percent'] = float(match.group(2))
                    stats['multiply_targetable_percent'] = float(match.group(3))
                    stats['partially_targetable_percent'] = float(match.group(4))
                    stats['not_targetable_percent'] = float(match.group(5))
                    stats['total_gos'] = int(match.group(6))
                
                # 個別カウントを抽出
                fully_match = re.search(r'(\d+)\s+Fully targetable GOs', content)
                if fully_match:
                    stats['fully_targetable_count'] = int(fully_match.group(1))
                
                simply_match = re.search(r'(\d+)\s+Fully and simply targetable GOs', content)
                if simply_match:
                    stats['simply_targetable_count'] = int(simply_match.group(1))
                
                multiply_match = re.search(r'(\d+)\s+Fully and multiply targetable GOs', content)
                if multiply_match:
                    stats['multiply_targetable_count'] = int(multiply_match.group(1))
                
                partially_match = re.search(r'(\d+)\s+Partially targetable GOs', content)
                if partially_match:
                    stats['partially_targetable_count'] = int(partially_match.group(1))
                
                not_match = re.search(r'(\d+)\s+Not targetable.*GOs', content)
                if not_match:
                    stats['not_targetable_count'] = int(not_match.group(1))
                
        except Exception as e:
            self.logger.error(f"GO-KFI統計解析エラー: {e}")
            stats = {}
        
        self.logger.info(f"GO-KFI統計解析完了: {stats}")
        return stats
    
    def create_go_reverse_lookup(self, go_annotated_path):
        """GO termから逆引き検索可能なデータ構造を作成"""
        self.logger.info("GO逆引き検索データを作成中...")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'go_reverse_lookup.tsv')
        
        # GOアノテーションファイルを読み込み
        go_to_genes = {}
        
        with open(go_annotated_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # GO termの列を特定
            go_columns = [i for i, col in enumerate(header) if 'GO' in col.upper()]
            
            for row in reader:
                if len(row) > max(go_columns) if go_columns else False:
                    # 遺伝子情報を取得
                    gene_info = {
                        'query_id': row[0] if len(row) > 0 else '',
                        'subject_id': row[1] if len(row) > 1 else '',
                        'description': row[6] if len(row) > 6 else '',
                        'unique_id': row[9] if len(row) > 9 else '',
                        'guide_targets': row[-1] if len(row) > 10 else ''
                    }
                    
                    # 各GO termに対して遺伝子情報を登録
                    for col_idx in go_columns:
                        if col_idx < len(row) and row[col_idx]:
                            go_terms = row[col_idx].split(';')
                            for go_term in go_terms:
                                go_term = go_term.strip()
                                if go_term:
                                    if go_term not in go_to_genes:
                                        go_to_genes[go_term] = []
                                    go_to_genes[go_term].append(gene_info)
        
        # 逆引きデータをファイルに保存
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("go_term\tgene_count\tgenes\n")
            
            for go_term, genes in go_to_genes.items():
                gene_count = len(genes)
                gene_list = []
                
                for gene in genes:
                    gene_str = f"{gene['unique_id']}|{gene['description']}|{gene['guide_targets']}"
                    gene_list.append(gene_str)
                
                f.write(f"{go_term}\t{gene_count}\t{';'.join(gene_list)}\n")
        
        self.logger.info(f"GO逆引き検索データ作成完了: {output_path}")
        return output_path
    
    def run(self):
        """GO解析の全体実行"""
        try:
            # 前段階の結果を確認
            input_path = os.path.join(self.intermediate_dir, 'annotation_with_guides_final.tsv')
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"アノテーション-ガイド結合ファイルが見つかりません: {input_path}")
            
            # カスタムDockerイメージのビルド
            self.build_docker_images()
            
            # GOアノテーションの追加
            go_annotated_path = self.add_go_annotations(input_path)
            
            # GO-KFI計算
            kfi_table_path, kfi_stat_path = self.compute_go_kfi(go_annotated_path)
            
            # Fully targetable GOsの抽出
            fully_targetable_path = self.extract_fully_targetable_gos(kfi_table_path)
            
            # GO親termの計算とvisualization
            parent_csv_path, parent_png_path = self.count_go_parent_terms(fully_targetable_path)
            
            # GO-KFI統計の解析
            kfi_stats = self.parse_go_kfi_statistics(kfi_stat_path)
            
            # GO逆引き検索データの作成
            reverse_lookup_path = self.create_go_reverse_lookup(go_annotated_path)
            
            # 結果の保存
            results = {
                'input_path': input_path,
                'go_annotated_path': go_annotated_path,
                'kfi_table_path': kfi_table_path,
                'kfi_stat_path': kfi_stat_path,
                'fully_targetable_path': fully_targetable_path,
                'parent_csv_path': parent_csv_path,
                'parent_png_path': parent_png_path,
                'reverse_lookup_path': reverse_lookup_path,
                'kfi_statistics': kfi_stats
            }
            
            # 結果をファイルに保存
            results_file = os.path.join(self.intermediate_dir, 'go_analysis_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== GO解析結果 ===\n\n")
                f.write(f"入力ファイル: {input_path}\n")
                f.write(f"GOアノテーション付きファイル: {go_annotated_path}\n")
                f.write(f"GO-KFIテーブル: {kfi_table_path}\n")
                f.write(f"GO-KFI統計: {kfi_stat_path}\n")
                f.write(f"Fully targetable GOs: {fully_targetable_path}\n")
                f.write(f"GO親term分析: {parent_csv_path}\n")
                f.write(f"GO親term可視化: {parent_png_path}\n")
                f.write(f"GO逆引き検索: {reverse_lookup_path}\n\n")
                
                f.write("=== GO-KFI統計 ===\n")
                for key, value in kfi_stats.items():
                    f.write(f"{key}: {value}\n")
            
            self.logger.info(f"GO解析結果を保存: {results_file}")
            self.logger.info("GO解析処理が完了しました")
            
            return results
            
        except Exception as e:
            self.logger.error(f"GO解析処理でエラーが発生しました: {e}")
            raise
