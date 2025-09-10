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
        self.protein_transcript_mapping = None
        
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
    
    def load_protein_transcript_mapping(self):
        """タンパク質ID→トランスクリプトIDマッピングファイルの読み込み"""
        mapping_path = self.config['input'].get('protein_transcript_mapping', '')
        
        if not mapping_path or not os.path.exists(mapping_path):
            self.logger.info("タンパク質-トランスクリプトマッピングファイルが指定されていません")
            return None
        
        self.logger.info(f"タンパク質-トランスクリプトマッピングファイルを読み込み中: {mapping_path}")
        
        mapping_dict = {}
        try:
            with open(mapping_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):  # 空行やコメント行をスキップ
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        transcript_id = parts[0].strip()
                        protein_id = parts[1].strip()
                        # protein_id -> transcript_id のマッピング
                        mapping_dict[protein_id] = transcript_id
                    else:
                        self.logger.warning(f"マッピングファイル行 {line_num}: 形式が不正です - {line}")
            
            self.logger.info(f"マッピング辞書を構築: {len(mapping_dict)} エントリ")
            self.protein_transcript_mapping = mapping_dict
            return mapping_dict
            
        except Exception as e:
            self.logger.error(f"マッピングファイル読み込みエラー: {e}")
            raise RuntimeError(f"マッピングファイルの読み込みに失敗しました: {mapping_path}")
    
    def run_fasta_search(self, query_path, db_path):
        """FASTA検索の実行"""
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'ggsearch_results.txt')
        
        # 既存ファイルのチェック
        if os.path.exists(output_path):
            # ファイルサイズをチェックして空でないことを確認
            if os.path.getsize(output_path) > 0:
                self.logger.info(f"既存のFASTA検索結果が見つかりました: {output_path}")
                self.logger.info("FASTA検索をスキップします")
                return output_path
            else:
                self.logger.warning(f"既存のFASTA検索結果ファイルが空です。再実行します: {output_path}")
        
        self.logger.info("UniProtデータベースに対するFASTA検索を開始...")
        
        # ggsearch36実行（intermediate_dirを作業ディレクトリとして使用）
        query_filename = os.path.basename(query_path)
        db_filename = os.path.basename(db_path)
        
        cmd = [
            "docker", "run",
            "-m", "20G",
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/fasta_work",
            "-w", "/fasta_work",
            "--rm",
            self.config['fasta_search']['docker_image'],
            "ggsearch36",
            "-Q", "-T", str(self.config['fasta_search']['threads']),
            "-d", "1", "-m10",
            "-E", str(self.config['fasta_search']['evalue']),
            query_filename, db_filename
        ]
        
        self.logger.debug(f"FASTA検索実行コマンド: {' '.join(cmd)}")
        
        # 標準出力をファイルに保存
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
        
        self.logger.debug(f"ggsearch output first 500 chars: {content[:500]}")
        
        # ggsearch36の実際のフォーマットに基づくパース
        lines = content.split('\n')
        current_query = None
        processing_hits = False
        
        for line in lines:
            line = line.strip()
            
            # クエリ開始行を検索: "1>>>query_id ... - xxx aa"
            query_match = re.match(r'^\d+>>>(.+?)\s+.+?\s+-\s+\d+\s+aa', line)
            if query_match:
                current_query = query_match.group(1).strip()
                processing_hits = False
                self.logger.debug(f"Found query: {current_query}")
                continue
            
            # ヒットセクション開始を検出: "The best scores are:"
            if "The best scores are:" in line:
                processing_hits = True
                continue
            
            # ヒット行をパース: "sp|ID|NAME description (length) score bits E-value"
            if processing_hits and current_query and line.startswith('sp|'):
                # より柔軟な正規表現でヒット行をパース
                hit_match = re.match(r'sp\|([^|]+)\|([^\s]+)\s+(.+?)\s+\(\s*(\d+)\)\s+(\d+)\s+([\d.]+)\s+([\d.e-]+)', line)
                if hit_match:
                    subject_acc = hit_match.group(1)
                    subject_name = hit_match.group(2)
                    description = hit_match.group(3).strip()
                    length = hit_match.group(4)
                    score = hit_match.group(5)
                    bits = hit_match.group(6)
                    evalue = hit_match.group(7)
                    
                    annotation = {
                        'query_id': current_query,
                        'subject_id': f"sp|{subject_acc}|{subject_name}",
                        'identity': '',  # ggsearch36では直接的な同一性情報がない
                        'length': length,
                        'evalue': evalue,
                        'bit_score': bits,
                        'description': description
                    }
                    
                    annotations.append(annotation)
                    self.logger.debug(f"Parsed hit: {subject_acc} for query {current_query}")
                    # 最初のヒットのみを取得して次のクエリへ
                    current_query = None
                    processing_hits = False
                    continue
                else:
                    self.logger.debug(f"Failed to parse hit line: {line}")
            
            # 次のクエリセクションまたはファイル終了でヒット処理を停止
            if processing_hits and (line.startswith('>>>') or not line):
                processing_hits = False
        
        # 結果をファイルに保存
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\n")
            for ann in annotations:
                f.write(f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                       f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t{ann['description']}\n")
        
        self.logger.info(f"Parsed {len(annotations)} annotations from ggsearch output")
    
    def _parse_ggsearch_line_based(self, content, output_path):
        """行ベースでggsearch結果をパース（フォールバック）"""
        annotations = []
        lines = content.split('\n')
        
        current_query = None
        in_hits_section = False
        
        for line in lines:
            line = line.strip()
            
            # クエリIDを検索
            if line.startswith('>>>') and ',' in line:
                current_query = line[3:].split(',')[0].strip()
                in_hits_section = False
                continue
            
            # ヒット開始を検出
            if current_query and '>>' in line and 'aa)' in line:
                in_hits_section = True
                
                # ヒット情報を抽出
                hit_match = re.search(r'>>(.+?)\s+\((\d+)\s+aa\)', line)
                if hit_match:
                    subject_info = hit_match.group(1).strip()
                    subject_length = hit_match.group(2)
                    
                    # 簡単なアノテーション情報を構築
                    annotation = {
                        'query_id': current_query,
                        'subject_id': subject_info.split()[0] if subject_info else '',
                        'identity': '',
                        'length': subject_length,
                        'evalue': '',
                        'bit_score': '',
                        'description': subject_info
                    }
                    
                    annotations.append(annotation)
                    # 1つ目のヒットのみを取得
                    current_query = None
                    in_hits_section = False
        
        # 結果をファイルに保存
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\n")
            for ann in annotations:
                f.write(f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                       f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t{ann['description']}\n")
        
        self.logger.info(f"Line-based parsing found {len(annotations)} annotations")
    
    def _extract_hits(self, query_content):
        """クエリ内容からヒット情報を抽出（改善版）"""
        hits = []
        
        # 複数のヒットパターンを試行
        hit_patterns = [
            r'>>(.+?)\s+\((\d+)\s+aa\)',  # 元のパターン
            r'>>(.+?)\s+\((\d+)\s+nt\)',  # 核酸配列の場合
            r'>(.+?)\s+Length:\s*(\d+)'    # 別のフォーマット
        ]
        
        for pattern in hit_patterns:
            matches = list(re.finditer(pattern, query_content))
            if matches:
                self.logger.debug(f"Found {len(matches)} hits with pattern: {pattern}")
                break
        
        if not matches:
            return hits
        
        for match in matches:
            subject_info = match.group(1).strip()
            subject_length = match.group(2)
            
            # スコア情報を抽出
            hit_start = match.start()
            next_match = re.search(r'>>', query_content[hit_start+1:])
            hit_end = next_match.start() + hit_start + 1 if next_match else len(query_content)
            
            hit_section = query_content[hit_start:hit_end]
            
            # E値、スコア、同一性を抽出（複数のパターンを試行）
            evalue = ''
            bit_score = ''
            identity = ''
            
            # E値の抽出
            evalue_patterns = [
                r'E\(\d+\):\s+([\d.e-]+)',
                r'E-value:\s*([\d.e-]+)',
                r'Expect\s*=\s*([\d.e-]+)'
            ]
            for e_pattern in evalue_patterns:
                e_match = re.search(e_pattern, hit_section)
                if e_match:
                    evalue = e_match.group(1)
                    break
            
            # スコアの抽出
            score_patterns = [
                r'(\d+\.?\d*)\s+bits',
                r'Score\s*=\s*(\d+\.?\d*)',
                r'bits:\s*(\d+\.?\d*)'
            ]
            for s_pattern in score_patterns:
                s_match = re.search(s_pattern, hit_section)
                if s_match:
                    bit_score = s_match.group(1)
                    break
            
            # 同一性の抽出
            identity_patterns = [
                r'identity:\s+(\d+\.?\d*)%',
                r'Identity\s*=\s*(\d+\.?\d*)%',
                r'(\d+\.?\d*)%\s+identity'
            ]
            for i_pattern in identity_patterns:
                i_match = re.search(i_pattern, hit_section)
                if i_match:
                    identity = i_match.group(1)
                    break
            
            hit = {
                'subject_id': subject_info.split()[0] if subject_info else '',
                'description': subject_info,
                'length': subject_length,
                'evalue': evalue,
                'bit_score': bit_score,
                'identity': identity
            }
            
            hits.append(hit)
        
        return hits
    
    def classify_genes(self, annotation_path):
        """遺伝子の分類（Multi/Single）"""
        self.logger.info("遺伝子の分類を開始...")
        
        # アノテーション結果の統計用カウンタ
        total_annotations = 0
        excluded_count = 0
        included_count = 0
        
        # subject_idをキーとしてアノテーション結果を読み込み
        subjects = {}
        
        with open(annotation_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    query_id = fields[0]
                    subject_id = fields[1]
                    description = fields[6]
                    
                    total_annotations += 1
                    
                    # existing protein modeでマッピングファイルが読み込まれている場合のフィルタリング
                    if self.protein_transcript_mapping is not None:
                        if query_id not in self.protein_transcript_mapping:
                            # マッピングファイルに存在しないIDは除外
                            excluded_count += 1
                            self.logger.debug(f"除外されたID（マッピングファイルに未登録）: {query_id}")
                            continue
                    
                    # 解析対象として処理
                    included_count += 1
                    
                    if subject_id not in subjects:
                        subjects[subject_id] = []
                    
                    subjects[subject_id].append({
                        'query_id': query_id,
                        'subject_id': subject_id,
                        'description': description,
                        'line_data': line.strip()
                    })
        
        # subject_idベースでMulti/Single判定と連番でunique_ID割り当て
        multi_counter = 1
        single_counter = 1
        subject_to_unique_id = {}
        
        for subject_id, entries in subjects.items():
            if len(entries) > 1:
                # Multi: 複数のquery_idが同じsubject_idにマッチ
                unique_id = f"Multi_{multi_counter:010d}"
                subject_to_unique_id[subject_id] = unique_id
                multi_counter += 1
            else:
                # Single: 1つのquery_idが1つのsubject_idにマッチ
                unique_id = f"Single_{single_counter:010d}"
                subject_to_unique_id[subject_id] = unique_id
                single_counter += 1
        
        # 分類結果を出力
        classified_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
        
        with open(classified_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\tgene_name\tclassification\tunique_id\n")
            
            # subject_idごとに結果を出力
            for subject_id, entries in subjects.items():
                unique_id = subject_to_unique_id[subject_id]
                classification = "Multi" if unique_id.startswith("Multi") else "Single"
                
                for entry in entries:
                    # 遺伝子名を抽出
                    gene_name = self._extract_gene_name(entry['description'])
                    
                    # 分類情報を追加して出力
                    f.write(f"{entry['line_data']}\t{gene_name}\t{classification}\t{unique_id}\n")
        
        # フィルタリング結果と統計情報のログ出力
        if self.protein_transcript_mapping is not None:
            self.logger.info(f"=== マッピングベースフィルタリング結果 ===")
            self.logger.info(f"総アノテーション数: {total_annotations}")
            self.logger.info(f"解析対象数: {included_count}")
            self.logger.info(f"除外数: {excluded_count}")
            if total_annotations > 0:
                inclusion_rate = (included_count / total_annotations) * 100
                self.logger.info(f"解析対象率: {inclusion_rate:.2f}%")
            
            if excluded_count > 0:
                self.logger.warning(f"{excluded_count} 個のIDがマッピングファイルに存在しないため除外されました")
        else:
            self.logger.info(f"TransDecoder mode: 全 {total_annotations} アノテーションを解析対象とします")
        
        self.logger.info(f"遺伝子分類結果を保存: {classified_path}")
        self.logger.info(f"Multi遺伝子数: {multi_counter - 1}, Single遺伝子数: {single_counter - 1}")
        self.logger.info(f"総解析対象遺伝子数: {(multi_counter - 1) + (single_counter - 1)}")
        
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
        
        # Multi遺伝子の配列抽出と個別ファイル作成
        if multi_genes:
            multi_seq_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            self._extract_sequences_by_ids(transcriptome_path, multi_seq_path, multi_genes)
            
            # Multi遺伝子を個別ファイルに分割
            self._create_individual_multi_gene_files(transcriptome_path, multi_genes)
            
            # アライメント処理を実行
            self._run_alignment_process()
        
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
                        sequence_id = self._resolve_sequence_id(query_id)
                        if sequence_id and sequence_id in sequences:
                            f.write(f">{unique_id}_{query_id}\n{sequences[sequence_id]}\n")
                        else:
                            self.logger.warning(f"配列が見つかりません (Multi): {query_id} -> {sequence_id}")
                else:
                    # Single遺伝子の場合
                    query_id = query_ids
                    sequence_id = self._resolve_sequence_id(query_id)
                    if sequence_id and sequence_id in sequences:
                        f.write(f">{unique_id}\n{sequences[sequence_id]}\n")
                    else:
                        self.logger.warning(f"配列が見つかりません (Single): {query_id} -> {sequence_id}")
    
    def _resolve_sequence_id(self, query_id):
        """
        クエリIDから実際のシーケンスIDを解決
        
        1. タンパク質-トランスクリプトマッピングがある場合、それを使用
        2. マッピングがない場合、従来の処理（.p1などを除去）を実行
        
        Args:
            query_id (str): クエリID（通常はprotein_id）
            
        Returns:
            str: 実際にFASTAファイル内で使用されるID、またはNone
        """
        # マッピングファイルが読み込まれている場合
        if self.protein_transcript_mapping:
            # protein_id -> transcript_id マッピングを試行
            if query_id in self.protein_transcript_mapping:
                transcript_id = self.protein_transcript_mapping[query_id]
                self.logger.debug(f"マッピング適用: {query_id} -> {transcript_id}")
                return transcript_id
            else:
                self.logger.debug(f"マッピングが見つかりません: {query_id}")
        
        # マッピングがない場合の従来の処理
        # query_idから .p1 などを除去してトランスクリプトIDとして使用
        base_id = query_id.split('.p')[0] if '.p' in query_id else query_id
        self.logger.debug(f"従来の処理適用: {query_id} -> {base_id}")
        return base_id
    
    def _validate_existing_protein_mode_requirements(self):
        """
        existing protein modeでの入力要件を検証
        
        既存のタンパク質配列を使用する場合、マッピングファイルが必要であることを確認
        """
        # existing protein modeかどうかを判定
        protein_sequences = self.config['input'].get('protein_sequences', '')
        is_existing_protein_mode = bool(protein_sequences and protein_sequences.strip())
        
        if is_existing_protein_mode:
            self.logger.info("existing protein modeを検出しました")
            
            # マッピングファイルの存在確認
            mapping_path = self.config['input'].get('protein_transcript_mapping', '')
            
            if not mapping_path or not mapping_path.strip():
                raise ValueError(
                    "existing protein modeではprotein_transcript_mappingパラメータが必須です。"
                    "タンパク質ID→トランスクリプトIDマッピングファイルを指定してください。"
                )
            
            if not os.path.exists(mapping_path):
                raise FileNotFoundError(
                    f"指定されたマッピングファイルが見つかりません: {mapping_path}\n"
                    f"existing protein modeではprotein_transcript_mappingファイルが必要です。"
                )
            
            # ファイル形式の基本的な確認
            try:
                with open(mapping_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()[:10]  # 最初の10行をチェック
                    
                valid_lines = 0
                for line_num, line in enumerate(lines, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        valid_lines += 1
                    else:
                        self.logger.warning(f"マッピングファイル行 {line_num}: TSV形式ではない可能性があります - {line}")
                
                if valid_lines == 0:
                    raise ValueError(
                        f"マッピングファイルに有効なデータが見つかりません: {mapping_path}\n"
                        f"ファイル形式がタブ区切り（TSV）であることを確認してください。"
                    )
                
                self.logger.info(f"マッピングファイル形式確認完了: {valid_lines} 有効行を検出")
                    
            except Exception as e:
                raise RuntimeError(f"マッピングファイルの読み込み/検証に失敗しました: {mapping_path}\nエラー: {e}")
            
        else:
            self.logger.info("TransDecoder modeを使用します（マッピングファイルは不要）")
    
    def _create_individual_multi_gene_files(self, transcriptome_path, multi_genes):
        """Multi遺伝子を個別ファイルに分割"""
        self.logger.info("Multi遺伝子を個別ファイルに分割中...")
        
        # multi_genesディレクトリを作成
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        Path(multi_genes_dir).mkdir(parents=True, exist_ok=True)
        
        # 全配列を読み込み
        sequences = {}
        with open(transcriptome_path, 'r') as f:
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
        
        # Multi遺伝子ごとに個別ファイルを作成
        for unique_id, query_ids in multi_genes.items():
            output_file = os.path.join(multi_genes_dir, f"{unique_id}.fa")
            
            with open(output_file, 'w') as f:
                for query_id in query_ids:
                    sequence_id = self._resolve_sequence_id(query_id)
                    if sequence_id and sequence_id in sequences:
                        f.write(f">{unique_id}_{query_id}\n{sequences[sequence_id]}\n")
                    else:
                        self.logger.warning(f"配列が見つかりません (Multi individual): {query_id} -> {sequence_id}")
            
            self.logger.debug(f"Multi遺伝子個別ファイル作成: {output_file}")
        
        self.logger.info(f"Multi遺伝子個別ファイル作成完了: {len(multi_genes)} ファイル")
    
    def _run_alignment_process(self):
        """アライメント処理を実行"""
        self.logger.info("Multi遺伝子アライメント処理を開始...")
        
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        multi_genes_aligned_dir = os.path.join(self.intermediate_dir, 'multi_genes_aligned')
        
        # 出力ディレクトリを作成
        Path(multi_genes_aligned_dir).mkdir(parents=True, exist_ok=True)
        
        # Dockerイメージをビルド
        docker_build_dir = "src/docker/run_alignment_and_mask"
        self.logger.info("run_alignment_and_maskのDockerイメージをビルド中...")
        
        build_cmd = [
            "sudo", "docker", "build",
            "-t", "run_alignment_and_mask/v1.0",
            docker_build_dir
        ]
        
        build_result = subprocess.run(build_cmd, capture_output=True, text=True)
        if build_result.returncode != 0:
            self.logger.error(f"Dockerビルドエラー: {build_result.stderr}")
            raise RuntimeError("Dockerイメージのビルドに失敗しました")
        
        self.logger.info("Dockerイメージビルド完了")
        
        # アライメント処理を実行（intermediate_dirを基準に相対パスで指定）
        multi_genes_rel_dir = os.path.relpath(multi_genes_dir, self.intermediate_dir)
        multi_genes_aligned_rel_dir = os.path.relpath(multi_genes_aligned_dir, self.intermediate_dir)
        
        alignment_cmd = [
            "docker", "run",
            "-m", "4G",
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/data",
            "-w", "/data",
            "--rm",
            "run_alignment_and_mask/v1.0",
            "python", "/app/run_alignment_and_mask.py",
            "-i", multi_genes_rel_dir,
            "-o", multi_genes_aligned_rel_dir
        ]
        
        self.logger.debug(f"アライメント実行コマンド: {' '.join(alignment_cmd)}")
        
        alignment_result = subprocess.run(alignment_cmd, capture_output=True, text=True)
        if alignment_result.returncode != 0:
            self.logger.error(f"アライメント処理エラー: {alignment_result.stderr}")
            raise RuntimeError("アライメント処理に失敗しました")
        
        # 実行結果をログ出力
        if alignment_result.stdout:
            self.logger.info(f"アライメント処理出力:\n{alignment_result.stdout}")
        
        self.logger.info("Multi遺伝子アライメント処理完了")
        
        # 結果ファイル数を確認
        aligned_files = list(Path(multi_genes_aligned_dir).glob("*_aligned_consensus.fa"))
        self.logger.info(f"生成されたアライメント済みファイル数: {len(aligned_files)}")
    
    def run(self):
        """アノテーション処理の全体実行"""
        try:
            # 前段階の結果を確認
            query_path = os.path.join(self.intermediate_dir, 'query.pep')
            if not os.path.exists(query_path):
                raise FileNotFoundError(f"予測タンパク質配列ファイルが見つかりません: {query_path}")
            
            # 入力検証: existing protein modeでのマッピングファイル要件チェック
            self._validate_existing_protein_mode_requirements()
            
            # タンパク質-トランスクリプトマッピングファイルの読み込み
            self.load_protein_transcript_mapping()
            
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
