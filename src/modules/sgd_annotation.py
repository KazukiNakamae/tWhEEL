#!/usr/bin/env python3
"""
Annotation Module
FASTA search against SGD database and functional annotation
"""

import os
import subprocess
import logging
import re
import csv
from pathlib import Path
import shutil

from .fasta_search_cache import FastaSearchCache
from .input_mode import EXISTING_PROTEIN_MODE, determine_protein_input_mode
from .multi_gene_alignment import (
    MultiGeneAlignmentRunner,
    annotation_workspace_lock,
    clear_multi_gene_artifacts,
    clear_single_gene_artifacts,
)
from .transcript_mapping import (
    TranscriptMappingError,
    load_mapping_file,
    read_fasta_sequences,
    validate_mapping_against_fastas,
)


class Sgd_annotation:
    """Annotation processing class"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.protein_transcript_mapping = None
        
    def prepare_sgd_db(self):
        """sgdデータベースの準備"""
        sgd_path = self.config['input']['sgd_db']
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'sgd.fasta')
        
        # ファイルの解凍とコピー
        if sgd_path.endswith('.gz'):
            self.logger.info("Unpacking SGD database file....")
            import gzip
            with gzip.open(sgd_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("Copying SGD database files....")
            shutil.copy2(sgd_path, output_path)
        
        self.logger.info(f"SGD database preparation complete.: {output_path}")
        return output_path
    
    def load_protein_transcript_mapping(self):
        """Load protein to transcript ID mapping file"""
        mapping_path = self.config['input'].get('protein_transcript_mapping', '')
        
        if not mapping_path or not os.path.exists(mapping_path):
            self.logger.info("No protein-transcript mapping file has been specified.")
            return None
        
        self.logger.info(f"Loading protein-transcript mapping file.: {mapping_path}")

        mapping_dict = load_mapping_file(mapping_path)
        self.logger.info(
            f"Built the mapping dictionary: {len(mapping_dict)} entries"
        )
        self.protein_transcript_mapping = mapping_dict
        return mapping_dict

    def validate_protein_transcript_mapping(
        self,
        protein_fasta_path,
        transcript_fasta_path,
    ):
        """Validate mapping coverage against the prepared pipeline FASTAs."""
        if self.protein_transcript_mapping is None:
            raise TranscriptMappingError(
                "Protein-transcript mapping has not been loaded"
            )

        validation = validate_mapping_against_fastas(
            self.protein_transcript_mapping,
            protein_fasta_path,
            transcript_fasta_path,
        )
        self.logger.info(
            "Protein-transcript mapping validation completed: "
            f"{validation['protein_count']} proteins, "
            f"{validation['mapped_transcript_count']} mapped transcripts"
        )
        return validation
    
    def run_fasta_search(self, query_path, db_path):
        """Execute FASTA search"""
        self.logger.info("Started a FASTA search on the SGD database...")
        
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
        
        self.logger.debug(f"FASTA search execution command: {' '.join(cmd)}")

        cache = FastaSearchCache(
            self.intermediate_dir,
            namespace='sgd',
            logger=self.logger,
        )
        return cache.run(
            cmd,
            query_path,
            db_path,
            search_settings={
                'docker_image': self.config['fasta_search']['docker_image'],
                'threads': self.config['fasta_search']['threads'],
                'evalue': self.config['fasta_search']['evalue'],
                'max_hits': 1,
                'output_format': 'm10',
            },
        )
    
    def parse_fasta_results(self, search_results_path):
        """Parse FASTA search results"""
        self.logger.info("Analyzing FASTA search results...")
        
        parsed_path = os.path.join(self.intermediate_dir, 'annotation_table.txt')
        
        # 結果パース用のPerlスクリプト相当の処理をPythonで実装
        self._parse_ggsearch_output(search_results_path, parsed_path)
        
        self.logger.info(f"Save annotation results: {parsed_path}")
        return parsed_path
    
    def _parse_ggsearch_output(self, input_path, output_path):
        """Parse ggsearch output and create annotation table"""
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
            
            # The FASTA description is optional, so capture only the first
            # whitespace-delimited identifier before the final length field.
            query_match = re.match(
                r'^\d+>>>(?P<query_id>\S+)(?:\s+.*?)?\s+-\s+\d+\s+aa(?:\s|$)',
                line,
            )
            if query_match:
                current_query = query_match.group('query_id')
                processing_hits = False
                self.logger.debug(f"Found query: {current_query}")
                continue
            
            # ヒットセクション開始を検出: "The best scores are:"
            if "The best scores are:" in line:
                processing_hits = True
                continue
            
            # Parse the first SGD hit and retain the gene symbol from the FASTA identifier.
            if processing_hits and current_query and 'SGDID:S' in line:
                hit_match = re.match(
                    r'^(?P<systematic_id>[^|\s]+)\|'
                    r'(?P<gene_symbol>[^|\s]+)\|'
                    r'(?P<sgd_id>SGDID:S\d{9})'
                    r'(?:\s+(?P<description>.*?))?\s+'
                    r'\(\s*(?P<length>\d+)\)\s+'
                    r'(?P<raw_score>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+'
                    r'(?P<bits>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+'
                    r'(?P<evalue>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$',
                    line,
                )
                if hit_match:
                    systematic_id = hit_match.group('systematic_id')
                    subject_acc = hit_match.group('sgd_id')
                    gene_symbol = hit_match.group('gene_symbol')
                    description = (hit_match.group('description') or '').strip()
                    length = hit_match.group('length')
                    bits = hit_match.group('bits')
                    evalue = hit_match.group('evalue')
                    
                    annotation = {
                        'query_id': current_query,
                        'subject_id': f"{subject_acc}",
                        'identity': '',  # ggsearch36では直接的な同一性情報がない
                        'length': length,
                        'evalue': evalue,
                        'bit_score': bits,
                        'description': description,
                        'gene_name': gene_symbol,
                        'annotation_source_database': 'SGD',
                        'sgd_id': subject_acc[len('SGDID:'):],
                        'ensembl_id': systematic_id,
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
            f.write(
                "query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t"
                "description\tgene_name\tannotation_source_database\tsgd_id\t"
                "ensembl_id\n"
            )
            for ann in annotations:
                f.write(f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                       f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t"
                       f"{ann['description']}\t{ann['gene_name']}\t"
                       f"{ann['annotation_source_database']}\t{ann['sgd_id']}\t"
                       f"{ann['ensembl_id']}\n")
        
        self.logger.info(f"Parsed {len(annotations)} annotations from ggsearch output")
    
    def _parse_ggsearch_line_based(self, content, output_path):
        """Parse ggsearch results line by line (fallback)"""
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
        """Extract hit information from query content (improved)"""
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
        """Classify genes (Multi/Single)"""
        self.logger.info("Started classifying genes...")
        
        # アノテーション結果の統計用カウンタ
        total_annotations = 0
        included_count = 0
        
        # subject_idをキーとしてアノテーション結果を読み込み
        subjects = {}
        
        with open(annotation_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            source_columns = (
                'query_id',
                'subject_id',
                'identity',
                'length',
                'evalue',
                'bit_score',
                'description',
            )
            for row_number, row in enumerate(reader, 2):
                if (
                    self.protein_transcript_mapping is not None
                    and any(row.get(column) is None for column in source_columns)
                ):
                    raise TranscriptMappingError(
                        "Malformed annotation row cannot be reconciled with "
                        "the protein-transcript mapping at "
                        f"{annotation_path}:{row_number}"
                    )

                query_id = (row.get('query_id') or '').strip()
                subject_id = (row.get('subject_id') or '').strip()
                description = (row.get('description') or '').strip()
                gene_name = (
                    row.get('gene_name', '').strip()
                    or self._extract_gene_name(description)
                )

                if not query_id or not subject_id:
                    if self.protein_transcript_mapping is not None:
                        raise TranscriptMappingError(
                            "Annotation row has an empty query_id or "
                            "subject_id and cannot be reconciled with the "
                            "protein-transcript mapping at "
                            f"{annotation_path}:{row_number}"
                        )
                    self.logger.warning(
                        "Skipping annotation row without query_id or subject_id."
                    )
                    continue

                annotation_source_database = row.get(
                    'annotation_source_database', ''
                ).strip()
                sgd_id = row.get('sgd_id', '').strip()
                ensembl_id = row.get('ensembl_id', '').strip()

                # Legacy tables can recover only the exact SGD identifier encoded
                # in subject_id. An Ensembl identifier must never be inferred from
                # a gene symbol or free-text description.
                if subject_id.startswith('SGDID:S'):
                    if not annotation_source_database:
                        annotation_source_database = 'SGD'
                    if not sgd_id:
                        sgd_id = subject_id[len('SGDID:'):]

                line_data = '\t'.join(
                    row.get(column) or '' for column in source_columns
                )
                    
                total_annotations += 1
                    
                # existing protein modeでマッピングファイルが読み込まれている場合のフィルタリング
                if self.protein_transcript_mapping is not None:
                    if query_id not in self.protein_transcript_mapping:
                        raise TranscriptMappingError(
                            "Annotation query ID is absent from the "
                            "validated protein-transcript mapping: "
                            f"{query_id!r}. Refusing to classify this "
                            "entry as untargetable. Rerun annotation "
                            "with matching query and mapping inputs."
                        )
                    
                # 解析対象として処理
                included_count += 1
                    
                if subject_id not in subjects:
                    subjects[subject_id] = []
                    
                subjects[subject_id].append({
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'description': description,
                    'gene_name': gene_name,
                    'line_data': line_data,
                    'annotation_source_database': annotation_source_database,
                    'sgd_id': sgd_id,
                    'ensembl_id': ensembl_id,
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
            f.write(
                "query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t"
                "description\tgene_name\tclassification\tunique_id\t"
                "annotation_source_database\tsgd_id\tensembl_id\n"
            )
            
            # subject_idごとに結果を出力
            for subject_id, entries in subjects.items():
                unique_id = subject_to_unique_id[subject_id]
                classification = "Multi" if unique_id.startswith("Multi") else "Single"
                
                for entry in entries:
                    gene_name = entry['gene_name']
                    
                    # 分類情報を追加して出力
                    f.write(
                        f"{entry['line_data']}\t{gene_name}\t{classification}\t"
                        f"{unique_id}\t{entry['annotation_source_database']}\t"
                        f"{entry['sgd_id']}\t{entry['ensembl_id']}\n"
                    )
        
        # フィルタリング結果と統計情報のログ出力
        if self.protein_transcript_mapping is not None:
            self.logger.info("=== Mapping-validated annotation results ===")
            self.logger.info(f"Total number of annotations: {total_annotations}")
            self.logger.info(f"Number of subjects analyzed: {included_count}")
        else:
            self.logger.info(f"TransDecoder mode: All {total_annotations} annotations will be analyzed.")
        
        self.logger.info(f"Save the gene classification results.: {classified_path}")
        self.logger.info(f"Multi-gene count: {multi_counter - 1}, Single gene count: {single_counter - 1}")
        self.logger.info(f"Total number of genes analyzed: {(multi_counter - 1) + (single_counter - 1)}")
        
        return classified_path
    
    def _extract_gene_name(self, description):
        """Extract gene name from description"""
        # sgdのdescriptionから遺伝子名を抽出
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
        """Extract one locked generation of Multi/Single gene sequences."""
        with annotation_workspace_lock(self.intermediate_dir):
            return self._extract_gene_sequences_locked(classified_path)

    def _extract_gene_sequences_locked(self, classified_path):
        """Extract gene sequences for Multi/Single genes"""
        self.logger.info("Gene sequence extraction started....")
        
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
        transcript_sequences = self._load_transcript_sequences(
            transcriptome_path
        )

        # Validate every Single and Multi member before creating any artifact.
        self._resolve_required_sequences(multi_genes, transcript_sequences)
        self._resolve_required_sequences(single_genes, transcript_sequences)

        # Multi遺伝子の配列抽出と個別ファイル作成
        if multi_genes:
            multi_seq_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            self._extract_sequences_by_ids(
                transcriptome_path,
                multi_seq_path,
                multi_genes,
                sequences=transcript_sequences,
            )
            
            # Split multi genes into individual files
            self._create_individual_multi_gene_files(
                transcriptome_path,
                multi_genes,
                sequences=transcript_sequences,
            )
            
            # Execute alignment process
            self._run_alignment_process(multi_genes)
        else:
            clear_multi_gene_artifacts(self.intermediate_dir)
        
        # Single遺伝子の配列抽出
        if single_genes:
            single_seq_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
            self._extract_sequences_by_ids(
                transcriptome_path,
                single_seq_path,
                single_genes,
                sequences=transcript_sequences,
            )
        else:
            clear_single_gene_artifacts(self.intermediate_dir)
        
        self.logger.info("Gene sequence extraction is complete.")
        
        return {
            'multi_genes': multi_genes,
            'single_genes': single_genes
        }
    
    @staticmethod
    def _load_transcript_sequences(fasta_path):
        """Load transcript sequences keyed by the first FASTA header token."""
        return read_fasta_sequences(fasta_path, 'transcript')

    def _extract_sequences_by_ids(
        self,
        fasta_path,
        output_path,
        gene_dict,
        sequences=None,
    ):
        """Extract sequences for specified IDs"""
        if sequences is None:
            sequences = self._load_transcript_sequences(fasta_path)

        resolved_sequences = self._resolve_required_sequences(
            gene_dict, sequences
        )

        # Output only after every classified query resolves.
        with open(output_path, 'w') as f:
            for unique_id, query_ids in gene_dict.items():
                if isinstance(query_ids, list):
                    # Multi遺伝子の場合
                    for query_id, sequence_id in resolved_sequences[unique_id]:
                        f.write(
                            f">{unique_id}_{query_id}\n"
                            f"{sequences[sequence_id]}\n"
                        )
                else:
                    # Single遺伝子の場合
                    _, sequence_id = resolved_sequences[unique_id][0]
                    f.write(f">{unique_id}\n{sequences[sequence_id]}\n")

    def _resolve_required_sequences(self, gene_dict, sequences):
        """Resolve every classified query before publishing sequence FASTAs."""
        resolved = {}
        missing = []

        for unique_id, query_ids in gene_dict.items():
            query_id_list = (
                query_ids if isinstance(query_ids, list) else [query_ids]
            )
            resolved[unique_id] = []
            for query_id in query_id_list:
                sequence_id = self._resolve_sequence_id(query_id)
                if sequence_id not in sequences:
                    missing.append(f"{query_id}->{sequence_id}")
                    continue
                resolved[unique_id].append((query_id, sequence_id))

        if missing:
            examples = ', '.join(missing[:5])
            if len(missing) > 5:
                examples += f", ... (+{len(missing) - 5})"
            if self.protein_transcript_mapping is not None:
                corrective_action = (
                    "correct the protein-transcript mapping and rerun "
                    "annotation"
                )
            else:
                corrective_action = (
                    "verify that query.pep uses TransDecoder-style "
                    "<transcript_id>.p<number> identifiers matching "
                    "transcriptome.fa, then rerun gene prediction and annotation"
                )
            raise TranscriptMappingError(
                "Classified protein query ID(s) could not be resolved to "
                "transcriptome FASTA records: "
                f"{examples}. Refusing to create incomplete gene FASTAs; "
                f"{corrective_action}."
            )

        return resolved
    
    def _resolve_sequence_id(self, query_id):
        """
        Resolve actual sequence ID from query ID
        
        1. タンパク質-トランスクリプトマッピングがある場合、それを使用
        2. マッピングがない場合、従来の処理（.p1などを除去）を実行
        
        Args:
            query_id (str): クエリID（通常はprotein_id）
            
        Returns:
            str: 実際にFASTAファイル内で使用されるID、またはNone
        """
        # マッピングファイルが読み込まれている場合
        if self.protein_transcript_mapping is not None:
            # protein_id -> transcript_id マッピングを試行
            if query_id in self.protein_transcript_mapping:
                transcript_id = self.protein_transcript_mapping[query_id]
                self.logger.debug(f"Apply mapping: {query_id} -> {transcript_id}")
                return transcript_id
            raise TranscriptMappingError(
                "Protein query ID is absent from the loaded "
                f"protein-transcript mapping: {query_id!r}"
            )
        
        # マッピングがない場合の従来の処理
        # query_idから .p1 などを除去してトランスクリプトIDとして使用
        base_id = re.sub(r'\.p\d+$', '', query_id)
        self.logger.debug(f"Conventional processing: {query_id} -> {base_id}")
        return base_id
    
    def _validate_existing_protein_mode_requirements(self):
        """
        Validate input requirements for existing protein mode
        
        既存のタンパク質配列を使用する場合、マッピングファイルが必要であることを確認
        """
        input_mode = determine_protein_input_mode(self.config)

        if input_mode == EXISTING_PROTEIN_MODE:
            self.logger.info("existing protein mode is detected")
            
            # マッピングファイルの存在確認
            mapping_path = self.config['input'].get('protein_transcript_mapping', '')
            
            if not mapping_path or not mapping_path.strip():
                raise ValueError(
                    "In existing protein mode, input.protein_transcript_mapping "
                    "is required. Specify a headerless "
                    "transcript_id<TAB>protein_id mapping file."
                )
            
            if not os.path.exists(mapping_path):
                raise FileNotFoundError(
                    f"The specified mapping file was not found: {mapping_path}\n"
                    "Existing protein mode requires a readable "
                    "input.protein_transcript_mapping file."
                )
            
        else:
            self.logger.info("Use TransDecoder mode (no mapping file required).")

        return input_mode
    
    def _create_individual_multi_gene_files(
        self,
        transcriptome_path,
        multi_genes,
        sequences=None,
    ):
        """Split multi genes into individual files"""
        self.logger.info("Split multi genes into individual files中...")
        
        # multi_genesディレクトリを作成
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        multi_genes_dir_path = Path(multi_genes_dir)
        multi_genes_dir_path.mkdir(parents=True, exist_ok=True)

        # Do not let groups from an older classification enter this run.
        expected_group_ids = set(multi_genes)
        for previous_group in multi_genes_dir_path.glob('*.fa'):
            if previous_group.stem not in expected_group_ids:
                previous_group.unlink()
        
        if sequences is None:
            sequences = self._load_transcript_sequences(transcriptome_path)

        resolved_sequences = self._resolve_required_sequences(
            multi_genes, sequences
        )

        # Multi遺伝子ごとに個別ファイルを作成
        for unique_id, query_ids in multi_genes.items():
            output_file = os.path.join(multi_genes_dir, f"{unique_id}.fa")
            
            with open(output_file, 'w') as f:
                for query_id, sequence_id in resolved_sequences[unique_id]:
                    f.write(
                        f">{unique_id}_{query_id}\n"
                        f"{sequences[sequence_id]}\n"
                    )
            
            self.logger.debug(f"Create individual files for multiple genes.: {output_file}")
        
        self.logger.info(f"Multi-gene individual file creation complete.: {len(multi_genes)} files")
    
    def _run_alignment_process(self, expected_multi_ids=None):
        """Run cached, bounded-parallel alignment for all Multi genes."""
        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            self.config,
            self.logger,
            run_command=subprocess.run,
        )
        return runner.run(expected_multi_ids)
    
    def run(self):
        """Execute one serialized annotation generation."""
        with annotation_workspace_lock(self.intermediate_dir):
            return self._run_locked()

    def _run_locked(self):
        """Execute entire annotation processing"""
        try:
            # 前段階の結果を確認
            query_path = os.path.join(self.intermediate_dir, 'query.pep')
            if not os.path.exists(query_path):
                raise FileNotFoundError(f"Predicted protein sequence file not found.: {query_path}")
            
            # 入力検証: existing protein modeでのマッピングファイル要件チェック
            input_mode = self._validate_existing_protein_mode_requirements()

            transcriptome_path = os.path.join(
                self.intermediate_dir, 'transcriptome.fa'
            )
            if input_mode == EXISTING_PROTEIN_MODE:
                # Validate the exact staged inputs before the expensive search.
                self.load_protein_transcript_mapping()
                self.validate_protein_transcript_mapping(
                    query_path,
                    transcriptome_path,
                )
            else:
                # A configured mapping belongs only to existing protein mode.
                self.protein_transcript_mapping = None
            
            # sgdデータベースの準備
            db_path = self.prepare_sgd_db()
            
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
                f.write("=== Annotation results ===\n\n")
                f.write(f"SGD database: {db_path}\n")
                f.write(f"FASTA search results: {search_results}\n")
                f.write(f"Annotation table: {annotation_table}\n")
                f.write(f"Classification table: {classified_table}\n")
                f.write(f"Multi-gene count: {len(gene_sequences['multi_genes'])}\n")
                f.write(f"Single gene count: {len(gene_sequences['single_genes'])}\n")
            
            self.logger.info(f"Save annotation results: {results_file}")
            self.logger.info("Annotation processing is complete.")
            
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during annotation processing.: {e}")
            raise
