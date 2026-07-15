#!/usr/bin/env python3
"""
Gene Ontology Analysis Module
GO Knock-out Feasibility Index (GO-KF)の計算とGO機能解析
"""

import os
import subprocess
import logging
import csv
from pathlib import Path


class GOAnalysis:
    """Gene Ontology analysis processing class"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.docker_dir = os.path.join('src', 'docker')
    
    def _check_docker_image_exists(self, image_name):
        """Check if a Docker image exists."""
        try:
            cmd = ["docker", "images", "-q", image_name]
            result = subprocess.run(cmd, capture_output=True, text=True)
            return bool(result.stdout.strip())
        except Exception as e:
            self.logger.warning(f"Error during Docker image existence check: {e}")
            return False
        
    def build_docker_images(self):
        """Build custom Docker images (with existing image check)"""
        self.logger.info("Start preparing your custom Docker image...")
        
        # 必要なDockerイメージのリスト
        docker_images = [
            ('add_go_columns', self.config['go_analysis']['add_go_image']),
            ('compute_go_kfi', self.config['go_analysis']['compute_kfi_image']),
            ('count_go_parent_terms', self.config['go_analysis']['count_parent_image'])
        ]
        
        build_needed = False
        
        for folder_name, image_name in docker_images:
            # 既存イメージのチェック
            if self._check_docker_image_exists(image_name):
                self.logger.info(f"A Docker image already exists (skip build).: {image_name}")
                continue
            
            # イメージが存在しない場合、ビルドを実行
            docker_path = os.path.join(self.docker_dir, folder_name)
            
            if os.path.exists(docker_path):
                self.logger.info(f"Building Docker image: {image_name}")
                build_needed = True
                
                cmd = [
                    "sudo", "docker", "build", 
                    "-t", image_name,
                    docker_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    self.logger.error(f"Docker build error ({image_name}): {result.stderr}")
                    
                    # 権限エラーの場合は追加情報を提供
                    if "permission denied" in result.stderr.lower():
                        self.logger.error("A Docker permission error occurred. Please check the following::")
                        self.logger.error("1. Verify that the Docker daemon is running.")
                        self.logger.error("2. Verify that the user belongs to the docker group.")
                        self.logger.error("3. Consider using 'sudo docker build'.")
                    
                    raise RuntimeError(f"The Docker image build failed.: {image_name}")
                
                self.logger.info(f"Docker image build complete: {image_name}")
            else:
                self.logger.warning(f"Docker directory not found.: {docker_path}")
                self.logger.warning(f"Processing cannot continue because the image '{image_name}' does not exist.")
                raise FileNotFoundError(f"Docker directory not found.: {docker_path}")
        
        if build_needed:
            self.logger.info("The custom Docker image build is complete.")
        else:
            self.logger.info("All necessary Docker images are already prepared.")
    
    def add_go_annotations(self, input_path):
        """Adding GO annotations"""
        self.logger.info("Adding GO annotations...")
        
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
        
        self.logger.debug(f"GO annotation addition command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GO annotation addition error: {result.stderr}")
            raise RuntimeError("Failed to add GO annotation.")
        
        self.logger.info(f"GO annotation added.: {output_path}")
        return output_path
    
    def compute_go_kfi(self, input_path):
        """Calculation of the GO Knock-out Feasibility Index (GO-KF)"""
        self.logger.info("Calculating GO-KFI...")
        
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
        
        self.logger.debug(f"GO-KFI calculation command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GO-KFI calculation error: {result.stderr}")
            raise RuntimeError("GO-KFI calculation failed.")
        
        # 結果をログに出力
        self.logger.info("GO-KFI calculation results:")
        self.logger.info(result.stdout)
        
        self.logger.info(f"GO-KFI calculation complete: {kfi_table_path}, {kfi_stat_path}")
        return kfi_table_path, kfi_stat_path
    
    def extract_fully_targetable_gos(self, kfi_table_path):
        """Extraction of Fully Targetable Goals"""
        self.logger.info("Extracting fully targetable GOs....")
        
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
        
        self.logger.info(f"Number of fully targetable GOs: {fully_targetable_count}")
        self.logger.info(f"Fully targetable GOs extraction completed: {output_path}")
        
        return output_path
    
    def count_go_parent_terms(self, fully_targetable_path):
        """GO親termの計算とvisualization"""
        self.logger.info("Performing calculation and visualization of GO parent term....")
        
        # 出力ファイルパス
        parent_csv_path = os.path.join(self.intermediate_dir, 'go_parent_counts.csv')
        parent_png_path = os.path.join(self.intermediate_dir, 'go_parent_counts.png')
        
        # Fully targetable GOsファイルの内容をチェック
        go_count = 0
        try:
            with open(fully_targetable_path, 'r') as f:
                # ヘッダー行をスキップ
                next(f)
                go_count = sum(1 for _ in f)
        except Exception as e:
            self.logger.warning(f"Fully targetable GOs file loading error: {e}")
        
        if go_count == 0:
            self.logger.warning("Since there are 0 fully targetable GOs, the GO parent term analysis will be skipped.")
            self.logger.info("No fully targetable GO terms were found in this dataset.")
            
            # 空の結果ファイルを作成
            with open(parent_csv_path, 'w') as f:
                f.write("go_term,count,parent_term\n")
            
            # 空のPNGファイルの代替として説明ファイルを作成
            no_data_message = os.path.join(self.intermediate_dir, 'go_parent_counts_no_data.txt')
            with open(no_data_message, 'w') as f:
                f.write("GO parent term visualization results\n")
                f.write("===================\n\n")
                f.write("Fully targetable GO terms: 0 pieces\n")
                f.write("The graph was not generated because there was no data to visualize.\n\n")
                f.write("Possible causes:\n")
                f.write("- CRISPR guide targets not designed for GO terms\n")
                f.write("- Threshold setting in GO-KFI calculation\n")
                f.write("- Filtering conditions for input data\n")
            
            self.logger.info(f"GO Parent term analysis skipped: Results explanation file {no_data_message}")
            return parent_csv_path, no_data_message
        
        self.logger.info(f"Fully targetable GOs to be processed: {go_count}個")
        
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
        
        self.logger.debug(f"GO parent term analysis command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"GO parent term analysis error: {result.stderr}")
            raise RuntimeError("GO parent term analysis failed.")
        
        self.logger.info(f"GO parent term analysis complete: {parent_csv_path}, {parent_png_path}")
        return parent_csv_path, parent_png_path
    
    def parse_go_kfi_statistics(self, kfi_stat_path):
        """GO-KFI統計の解析"""
        self.logger.info("Analyzing GO-KFI statistics...")
        
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
            self.logger.error(f"GO-KFI statistical analysis error: {e}")
            stats = {}
        
        self.logger.info(f"GO-KFI statistical analysis completed: {stats}")
        return stats
    
    def create_go_reverse_lookup(self, go_annotated_path):
        """GO termから逆引き検索可能なデータ構造を作成"""
        self.logger.info("Creating GO reverse lookup data....")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'go_reverse_lookup.tsv')
        
        # CSVフィールドサイズ制限を増加
        max_field_size = csv.field_size_limit()
        try:
            csv.field_size_limit(1000000)  # 1MBに増加
        except OverflowError:
            csv.field_size_limit(max_field_size)  # システムの最大値を使用
            self.logger.warning(f"Set the CSV field size limit to the maximum value {max_field_size}.")
        
        # GOアノテーションファイルを読み込み
        go_to_genes = {}
        
        try:
            with open(go_annotated_path, 'r', encoding='utf-8') as f:
                # 手動でTSVを解析（csvモジュールを使わない方法）
                header_line = f.readline().strip()
                header = header_line.split('\t')
                
                # GO termの列を特定
                go_columns = [i for i, col in enumerate(header) if 'GO' in col.upper()]
                self.logger.debug(f"GO annotation column: {[header[i] for i in go_columns]}")
                
                line_count = 0
                for line in f:
                    line_count += 1
                    if line_count % 10000 == 0:
                        self.logger.debug(f"処理済み行数: {line_count}")
                    
                    try:
                        row = line.strip().split('\t')
                        
                        if len(row) > max(go_columns) if go_columns else False:
                            # 遺伝子情報を取得
                            gene_info = {
                                'query_id': row[0] if len(row) > 0 else '',
                                'subject_id': row[1] if len(row) > 1 else '',
                                'description': row[6] if len(row) > 6 else '',
                                'unique_id': row[9] if len(row) > 9 else '',
                                'guide_targets': row[-1] if len(row) > 10 else ''
                            }
                            
                            # guide_targetsが非常に長い場合は切り詰め
                            if len(gene_info['guide_targets']) > 50000:
                                gene_info['guide_targets'] = gene_info['guide_targets'][:50000] + "...[truncated]"
                            
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
                    
                    except Exception as e:
                        self.logger.warning(f"Error processing line {line_count}: {e}")
                        continue
                
                self.logger.info(f"Total number of rows processed: {line_count}")
        
        except Exception as e:
            self.logger.error(f"GO annotation file loading error: {e}")
            raise
        finally:
            # CSVフィールドサイズ制限を元に戻す
            csv.field_size_limit(max_field_size)
        
        # 逆引きデータをファイルに保存
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("go_term\tgene_count\tgenes\n")
            
            for go_term, genes in go_to_genes.items():
                gene_count = len(genes)
                gene_list = []
                
                for gene in genes:
                    # 長すぎる場合は切り詰める
                    description = gene['description'][:500] + "..." if len(gene['description']) > 500 else gene['description']
                    guide_targets = gene['guide_targets'][:1000] + "..." if len(gene['guide_targets']) > 1000 else gene['guide_targets']
                    
                    gene_str = f"{gene['unique_id']}|{description}|{guide_targets}"
                    gene_list.append(gene_str)
                
                f.write(f"{go_term}\t{gene_count}\t{';'.join(gene_list)}\n")
        
        self.logger.info(f"GO reverse lookup data creation complete.: {output_path}")
        self.logger.info(f"Total number of GO terms: {len(go_to_genes)}")
        return output_path
    
    def run(self):
        """Execute entire GO analysis"""
        try:
            # 前段階の結果を確認
            input_path = os.path.join(self.intermediate_dir, 'annotation_with_guides_final.tsv')
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Annotation-guide combined file not found.: {input_path}")
            
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
                f.write("=== GO analysis results ===\n\n")
                f.write(f"input file: {input_path}\n")
                f.write(f"GO annotated file: {go_annotated_path}\n")
                f.write(f"GO-KFI Table: {kfi_table_path}\n")
                f.write(f"GO-KFI Statistics: {kfi_stat_path}\n")
                f.write(f"Fully targetable GOs: {fully_targetable_path}\n")
                f.write(f"GO Parent Term Analysis: {parent_csv_path}\n")
                f.write(f"GO parent term visualization: {parent_png_path}\n")
                f.write(f"GO Reverse Lookup: {reverse_lookup_path}\n\n")
                
                f.write("=== GO-KFI Statistics ===\n")
                for key, value in kfi_stats.items():
                    f.write(f"{key}: {value}\n")
            
            self.logger.info(f"Save GO analysis results: {results_file}")
            self.logger.info("GO analysis process completed.")
            
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during the GO parsing process.: {e}")
            raise
