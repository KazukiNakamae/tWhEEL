#!/usr/bin/env python3
"""
アノテーションモジュール
UniProtデータベースに対するFASTA検索と機能アノテーション
"""

import os
import subprocess
import logging
import re
from pathlib import Path
import shutil


class Annotation:
    """アノテーション処理クラス"""
    
    def __init__(self, config):
        """
        初期化
        
        Args:
            config (dict): 設定辞書
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        
    def prepare_uniprot_db(self):
        """UniProtデータベースの準備"""
        uniprot_path = self.config['input']['uniprot_db']
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'uniprot_sprot.fasta')
        
        # ファイルの解凍とコピー
        if uniprot_path.endswith('.gz'):
            self.logger.info("UniProtデータベースファイルを解凍中...")
            import gzip
            with gzip.open(uniprot_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("UniProtデータベースファイルをコピー中...")
            shutil.copy2(uniprot_path, output_path)
        
        self.logger.info(f"UniProtデータベースの準備が完了: {output_path}")
        return output_path
    
    def run_fasta_search(self, query_path, db_path):
        """FASTA検索の実行"""
        self.logger.info("UniProtデータベースに対するFASTA検索を開始...")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'ggsearch_results.txt')
        
        # ggsearch36実行
        cmd = [
            "docker", "run",
            "-m", "20G",
            "-v", f"{os.getcwd()}:/data",
            "-w", "/data",
            "--rm",
            self.config['fasta_search']['docker_image'],
            "ggsearch36",
            "-Q", "-T", str(self.config['fasta_search']['threads']),
            "-d", "1", "-m10",
            "-E", str(self.config['fasta_search']['evalue']),
            query_path, db_path
        ]
        
        self.logger.debug(f"FASTA検索実行コマンド: {' '.join(cmd)}")
        
        with open(output_path, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"FASTA検索実行エラー: {result.stderr}")
            raise RuntimeError("FASTA検索に失敗しました")
        
        self.logger.info(f"FASTA検索結果を保存: {output_path}")
        return output_path
    
    def parse_fasta_results(self, search_results_path):
        """FASTA検索結果のパース"""
        self.logger.info("FASTA検索結果を解析中...")
        
        parsed_path = os.path.join(self.intermediate_dir, 'annotation_table.txt')
        
        # 結果パース用のPerlスクリプト相当の処理をPythonで実装
        self._parse_ggsearch_output(search_results_path, parsed_path)
        
        self.logger.info(f"アノテーション結果を保存: {parsed_path}")
        return parsed_path
    
    def _parse_ggsearch_output(self, input_path, output_path):
        """ggsearch出力をパースしてアノテーションテーブルを作成"""
        annotations = []
        
        with open(input_path, 'r') as f:
            content = f.read()
        
        # 各クエリのセクションを分割
        query_sections = re.split(r'^Query: (.+)$', content, flags=re.MULTILINE)
        
        for i in range(1, len(query_sections), 2):
            query_id = query_sections[i].strip()
            query_content = query_sections[i+1] if i+1 < len(query_sections) else ""
            
            # ヒットを抽出
            hits = self._extract_hits(query_content)
            
            if hits:
                # 最良のヒットを選択
                best_hit = hits[0]
                
                # アノテーション情報を構築
                annotation = {
                    'query_id': query_id,
                    'subject_id': best_hit.get('subject_id', ''),
                    'identity': best_hit.get('identity', ''),
                    'length': best_hit.get('length', ''),
                    'evalue': best_hit.get('evalue', ''),
                    'bit_score': best_hit.get('bit_score', ''),
                    'description': best_hit.get('description', '')
                }
                
                annotations.append(annotation)
        
        # 結果をファイルに保存
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\n")
            for ann in annotations:
                f.write(f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                       f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t{ann['description']}\n")
    
    def _extract_hits(self, query_content):
        """クエリ内容からヒット情報を抽出"""
        hits = []
        
        # ヒットセクションを検索
        hit_pattern = r'>>(.+?)\s+\((\d+)\s+aa\)'
        matches = re.finditer(hit_pattern, query_content)
        
        for match in matches:
            subject_info = match.group(1).strip()
            subject_length = match.group(2)
            
            # スコア情報を抽出
            hit_start = match.start()
            next_match = re.search(r'>>', query_content[hit_start+1:])
            hit_end = next_match.start() + hit_start + 1 if next_match else len(query_content)
            
            hit_section = query_content[hit_start:hit_end]
            
            # E値、スコア、同一性を抽出
            score_match = re.search(r'E\(\d+\):\s+([\d.e-]+).*?(\d+\.?\d*)\s+bits', hit_section)
            identity_match = re.search(r'identity:\s+(\d+\.?\d*)%', hit_section)
            
            hit = {
                'subject_id': subject_info.split()[0] if subject_info else '',
                'description': subject_info,
                'length': subject_length,
                'evalue': score_match.group(1) if score_match else '',
                'bit_score': score_match.group(2) if score_match else '',
                'identity': identity_match.group(1) if identity_match else ''
            }
            
            hits.append(hit)
        
        return hits
    
    def classify_genes(self, annotation_path):
        """遺伝子の分類（Multi/Single）"""
        self.logger.info("遺伝子の分類を開始...")
        
        # アノテーション結果を読み込み
        genes = {}
        
        with open(annotation_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    query_id = fields[0]
                    subject_id = fields[1]
                    description = fields[6]
                    
                    # 遺伝子名を抽出
                    gene_name = self._extract_gene_name(description)
                    
                    if gene_name not in genes:
                        genes[gene_name] = []
                    
                    genes[gene_name].append({
                        'query_id': query_id,
                        'subject_id': subject_id,
                        'description': description
                    })
        
        # Multi/Single分類
        classified_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
        
        with open(classified_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\tgene_name\tclassification\tunique_id\n")
            
            # 元のファイルを再読み込みして分類情報を追加
            with open(annotation_path, 'r', encoding='utf-8') as f_orig:
                header = f_orig.readline().strip()
                for line in f_orig:
                    fields = line.strip().split('\t')
                    if len(fields) >= 7:
                        query_id = fields[0]
                        description = fields[6]
                        gene_name = self._extract_gene_name(description)
                        
                        # 分類を決定
                        if len(genes[gene_name]) > 1:
                            classification = "Multi"
                            unique_id = f"Multi_{gene_name}_{hash(query_id) % 100000:05d}"
                        else:
                            classification = "Single"
                            unique_id = f"Single_{gene_name}_{hash(query_id) % 100000:05d}"
                        
                        # 分類情報を追加して出力
                        f.write(f"{line.strip()}\t{gene_name}\t{classification}\t{unique_id}\n")
        
        self.logger.info(f"遺伝子分類結果を保存: {classified_path}")
        return classified_path
    
    def _extract_gene_name(self, description):
        """遺伝子名をdescriptionから抽出"""
        # UniProtのdescriptionから遺伝子名を抽出
        # 例: "sp|P12345|GENE_HUMAN Gene name GN=GENE PE=1 SV=1"
        
        # GN=パターンを検索
        gn_match = re.search(r'GN=([^\s]+)', description)
        if gn_match:
            return gn_match.group(1)
        
        # フォールバック: 最初の単語を使用
        parts = description.split()
        if len(parts) > 0:
            return parts[0].split('|')[-1] if '|' in parts[0] else parts[0]
        
        return "unknown"
    
    def extract_gene_sequences(self, classified_path):
        """Multi/Single遺伝子の配列を抽出"""
        self.logger.info("遺伝子配列の抽出を開始...")
        
        # 分類されたファイルを読み込み
        multi_genes = {}
        single_genes = {}
        
        with open(classified_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 10:
                    query_id = fields[0]
                    classification = fields[8]
                    unique_id = fields[9]
                    
                    if classification == "Multi":
                        if unique_id not in multi_genes:
                            multi_genes[unique_id] = []
                        multi_genes[unique_id].append(query_id)
                    else:
                        single_genes[unique_id] = query_id
        
        # 元のトランスクリプトームファイルから配列を抽出
        transcriptome_path = os.path.join(self.intermediate_dir, 'transcriptome.fa')
        
        # Multi遺伝子の配列抽出
        if multi_genes:
            multi_seq_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            self._extract_sequences_by_ids(transcriptome_path, multi_seq_path, multi_genes)
        
        # Single遺伝子の配列抽出
        if single_genes:
            single_seq_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
            self._extract_sequences_by_ids(transcriptome_path, single_seq_path, single_genes)
        
        self.logger.info("遺伝子配列の抽出が完了しました")
        
        return {
            'multi_genes': multi_genes,
            'single_genes': single_genes
        }
    
    def _extract_sequences_by_ids(self, fasta_path, output_path, gene_dict):
        """指定されたIDの配列を抽出"""
        # 全配列を読み込み
        sequences = {}
        
        with open(fasta_path, 'r') as f:
            current_id = None
            current_seq = ""
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = current_seq
                    current_id = line[1:].split()[0]  # 最初の部分のみを取得
                    current_seq = ""
                else:
                    current_seq += line
            
            # 最後の配列を処理
            if current_id and current_seq:
                sequences[current_id] = current_seq
        
        # 指定されたIDの配列を出力
        with open(output_path, 'w') as f:
            for unique_id, query_ids in gene_dict.items():
                if isinstance(query_ids, list):
                    # Multi遺伝子の場合
                    for query_id in query_ids:
                        if query_id in sequences:
                            f.write(f">{unique_id}_{query_id}\n{sequences[query_id]}\n")
                else:
                    # Single遺伝子の場合
                    query_id = query_ids
                    if query_id in sequences:
                        f.write(f">{unique_id}\n{sequences[query_id]}\n")
    
    def run(self):
        """アノテーション処理の全体実行"""
        try:
            # 前段階の結果を確認
            query_path = os.path.join(self.intermediate_dir, 'query.pep')
            if not os.path.exists(query_path):
                raise FileNotFoundError(f"予測タンパク質配列ファイルが見つかりません: {query_path}")
            
            # UniProtデータベースの準備
            db_path = self.prepare_uniprot_db()
            
            # FASTA検索
            search_results = self.run_fasta_search(query_path, db_path)
            
            # 結果のパース
            annotation_table = self.parse_fasta_results(search_results)
            
            # 遺伝子分類
            classified_table = self.classify_genes(annotation_table)
            
            # 遺伝子配列の抽出
            gene_sequences = self.extract_gene_sequences(classified_table)
            
            # 結果の保存
            results = {
                'db_path': db_path,
                'search_results': search_results,
                'annotation_table': annotation_table,
                'classified_table': classified_table,
                'gene_sequences': gene_sequences
            }
            
            # 結果をファイルに保存
            results_file = os.path.join(self.intermediate_dir, 'annotation_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== アノテーション結果 ===\n\n")
                f.write(f"UniProtデータベース: {db_path}\n")
                f.write(f"FASTA検索結果: {search_results}\n")
                f.write(f"アノテーションテーブル: {annotation_table}\n")
                f.write(f"分類テーブル: {classified_table}\n")
                f.write(f"Multi遺伝子数: {len(gene_sequences['multi_genes'])}\n")
                f.write(f"Single遺伝子数: {len(gene_sequences['single_genes'])}\n")
            
            self.logger.info(f"アノテーション結果を保存: {results_file}")
            self.logger.info("アノテーション処理が完了しました")
            
            return results
            
        except Exception as e:
            self.logger.error(f"アノテーション処理でエラーが発生しました: {e}")
            raise
