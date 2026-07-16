#!/usr/bin/env python3
"""
CRISPRガイド設計モジュール
CrisFlashを使用したSpCas9ガイドRNA設計
"""

import os
import subprocess
import logging
import gzip
import shutil

from .docker_workspace import DockerWorkspace
import csv
from pathlib import Path


class GuideDesign:
    """CRISPRガイド設計処理クラス"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        
    def prepare_genome_assembly(self):
        """Prepare genome assembly file (using transcriptome as fallback)"""
        genome_path = self.config['input'].get('genome_assembly', '')
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'reference_for_offtarget.fa')
        
        # genome_assemblyが提供されているかチェック
        if genome_path and os.path.exists(genome_path):
            self.logger.info("Genome assembly provided - using for comprehensive off-target search")
            source_path = genome_path
            source_type = "genome assembly"
        else:
            # transcriptomeをフォールバックとして使用
            transcriptome_path = self.config['input']['transcriptome']
            self.logger.warning("Genome assembly not provided - using transcriptome for off-target search")
            self.logger.warning("WARNING: Off-target search using transcriptome may be less comprehensive than using genome assembly")
            source_path = transcriptome_path
            source_type = "transcriptome"
        
        # ファイルの解凍とコピー
        if source_path.endswith('.gz'):
            self.logger.info(f"Extracting {source_type} file...")
            with gzip.open(source_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info(f"Copying {source_type} file...")
            shutil.copy2(source_path, output_path)
        
        # 配列を1行にフォーマット
        self._format_fasta_oneline(output_path)
        
        # Convert degenerate bases to N (for CrisFlash compatibility)
        self._clean_degenerate_bases(output_path)
        
        self.logger.info(f"Reference sequence prepared for off-target search: {output_path}")
        return output_path
    
    def _format_fasta_oneline(self, fasta_path):
        """FASTAファイルを1行形式に変換"""
        temp_path = fasta_path + '.temp'
        
        with open(fasta_path, 'r') as f_in:
            with open(temp_path, 'w') as f_out:
                sequence = ""
                header = ""
                
                for line in f_in:
                    line = line.strip()
                    if line.startswith('>'):
                        if header and sequence:
                            f_out.write(f"{header}\n{sequence}\n")
                        header = line
                        sequence = ""
                    else:
                        sequence += line
                
                # 最後のシーケンスを出力
                if header and sequence:
                    f_out.write(f"{header}\n{sequence}\n")
        
        # 元のファイルを置換
        shutil.move(temp_path, fasta_path)
        self.logger.debug("The FASTA file has been converted to single-line format.")
    
    def _clean_degenerate_bases(self, fasta_path):
        """Convert degenerate bases to N (for CrisFlash compatibility)"""
        self.logger.info("Converting degenerate base to N...")
        
        # 縮重塩基の定義（大文字・小文字両方）
        degenerate_bases = {
            'R': 'N', 'r': 'N',  # A or G
            'Y': 'N', 'y': 'N',  # C or T
            'S': 'N', 's': 'N',  # G or C
            'W': 'N', 'w': 'N',  # A or T
            'K': 'N', 'k': 'N',  # G or T
            'M': 'N', 'm': 'N',  # A or C
            'B': 'N', 'b': 'N',  # C or G or T
            'D': 'N', 'd': 'N',  # A or G or T
            'H': 'N', 'h': 'N',  # A or C or T
            'V': 'N', 'v': 'N'   # A or C or G
        }
        
        temp_path = fasta_path + '.clean_temp'
        conversion_count = 0
        total_bases = 0
        
        with open(fasta_path, 'r') as f_in:
            with open(temp_path, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        # ヘッダー行はそのまま出力
                        f_out.write(line)
                    else:
                        # 配列行の縮重塩基を変換
                        cleaned_sequence = ""
                        for base in line.strip():
                            total_bases += 1
                            if base in degenerate_bases:
                                cleaned_sequence += degenerate_bases[base]
                                conversion_count += 1
                            else:
                                cleaned_sequence += base
                        f_out.write(cleaned_sequence + '\n')
        
        # 元のファイルを置換
        shutil.move(temp_path, fasta_path)
        
        # 統計情報をログ出力
        if conversion_count > 0:
            conversion_rate = (conversion_count / total_bases) * 100 if total_bases > 0 else 0
            self.logger.info(f"Degenerate base conversion complete: {conversion_count}/{total_bases} bases ({conversion_rate:.2f}%) converted to N")
        else:
            self.logger.info("No degenerate bases were found.")
    
    @staticmethod
    def _multi_gene_ids_from_fasta(multi_genes_path):
        """Return exact Multi group IDs represented in a grouped FASTA file."""
        path = Path(multi_genes_path)
        if not path.is_file():
            return set()

        multi_ids = set()
        with path.open('r', encoding='utf-8') as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith('>'):
                    continue
                header_text = line[1:].strip()
                identifier = (
                    header_text.split(maxsplit=1)[0] if header_text else ''
                )
                identifier_parts = identifier.split('_')
                if (
                    len(identifier_parts) < 2
                    or identifier_parts[0] != 'Multi'
                    or not identifier_parts[1]
                ):
                    raise ValueError(
                        "Invalid Multi-gene FASTA identifier at "
                        f"{path}:{line_number}: {identifier!r}"
                    )
                multi_ids.add('_'.join(identifier_parts[:2]))
        return multi_ids

    def _expected_multi_gene_ids(self, multi_genes_path):
        """Reconcile declared Multi groups with the extracted sequence FASTA."""
        fasta_multi_ids = self._multi_gene_ids_from_fasta(multi_genes_path)
        classified_path = Path(self.intermediate_dir) / 'annotation_classified.txt'
        if not classified_path.is_file():
            return fasta_multi_ids

        with classified_path.open(
            'r', encoding='utf-8', newline=''
        ) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required_columns = {'classification', 'unique_id'}
            available_columns = set(reader.fieldnames or [])
            missing_columns = sorted(required_columns - available_columns)
            if missing_columns:
                raise ValueError(
                    "Annotation classification file is missing required "
                    "columns: " + ','.join(missing_columns)
                )

            declared_multi_ids = set()
            for row_number, row in enumerate(reader, 2):
                if (row.get('classification') or '').strip() != 'Multi':
                    continue
                unique_id = (row.get('unique_id') or '').strip()
                if not unique_id:
                    raise ValueError(
                        "Multi-gene annotation has no unique_id at "
                        f"{classified_path}:{row_number}"
                    )
                declared_multi_ids.add(unique_id)

        missing_multi_ids = sorted(declared_multi_ids - fasta_multi_ids)
        unexpected_multi_ids = sorted(fasta_multi_ids - declared_multi_ids)
        if missing_multi_ids or unexpected_multi_ids:
            details = []
            if missing_multi_ids:
                details.append('missing=' + ','.join(missing_multi_ids))
            if unexpected_multi_ids:
                details.append(
                    'unexpected=' + ','.join(unexpected_multi_ids)
                )
            raise RuntimeError(
                "Multi-gene sequence FASTA does not match annotation "
                "declarations: " + '; '.join(details)
            )

        return declared_multi_ids

    @staticmethod
    def _single_gene_ids_from_fasta(single_genes_path):
        """Return validated Single IDs represented in an extracted FASTA file."""
        path = Path(single_genes_path)
        if not path.is_file():
            return set()

        single_ids = set()
        current_id = None
        current_sequence_length = 0
        with path.open('r', encoding='utf-8') as handle:
            for line_number, line in enumerate(handle, 1):
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                if stripped_line.startswith('>'):
                    if current_id is not None and current_sequence_length == 0:
                        raise RuntimeError(
                            "Single-gene FASTA contains an empty sequence for "
                            f"{current_id}: {path}"
                        )
                    header_text = stripped_line[1:].strip()
                    identifier = (
                        header_text.split(maxsplit=1)[0]
                        if header_text else ''
                    )
                    if not identifier:
                        raise ValueError(
                            "Invalid Single-gene FASTA identifier at "
                            f"{path}:{line_number}"
                        )
                    if identifier in single_ids:
                        raise RuntimeError(
                            "Duplicate Single-gene FASTA identifier at "
                            f"{path}:{line_number}: {identifier}"
                        )
                    single_ids.add(identifier)
                    current_id = identifier
                    current_sequence_length = 0
                else:
                    if current_id is None:
                        raise RuntimeError(
                            f"Invalid Single-gene FASTA sequence at "
                            f"{path}:{line_number}"
                        )
                    current_sequence_length += len(stripped_line)

        if current_id is not None and current_sequence_length == 0:
            raise RuntimeError(
                "Single-gene FASTA contains an empty sequence for "
                f"{current_id}: {path}"
            )
        return single_ids

    def _expected_single_gene_ids(self, single_genes_path):
        """Reconcile declared Single genes with the extracted sequence FASTA."""
        fasta_single_ids = self._single_gene_ids_from_fasta(single_genes_path)
        classified_path = Path(self.intermediate_dir) / 'annotation_classified.txt'
        if not classified_path.is_file():
            return fasta_single_ids

        with classified_path.open(
            'r', encoding='utf-8', newline=''
        ) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required_columns = {'classification', 'unique_id'}
            available_columns = set(reader.fieldnames or [])
            missing_columns = sorted(required_columns - available_columns)
            if missing_columns:
                raise ValueError(
                    "Annotation classification file is missing required "
                    "columns: " + ','.join(missing_columns)
                )

            declared_single_ids = set()
            for row_number, row in enumerate(reader, 2):
                if (row.get('classification') or '').strip() != 'Single':
                    continue
                unique_id = (row.get('unique_id') or '').strip()
                if not unique_id:
                    raise ValueError(
                        "Single-gene annotation has no unique_id at "
                        f"{classified_path}:{row_number}"
                    )
                declared_single_ids.add(unique_id)

        missing_single_ids = sorted(declared_single_ids - fasta_single_ids)
        unexpected_single_ids = sorted(fasta_single_ids - declared_single_ids)
        if missing_single_ids or unexpected_single_ids:
            details = []
            if missing_single_ids:
                details.append('missing=' + ','.join(missing_single_ids))
            if unexpected_single_ids:
                details.append(
                    'unexpected=' + ','.join(unexpected_single_ids)
                )
            raise RuntimeError(
                "Single-gene sequence FASTA does not match annotation "
                "declarations: " + '; '.join(details)
            )

        return declared_single_ids

    @staticmethod
    def _validate_multi_consensus_file(consensus_path, expected_multi_id):
        """Require one non-empty FASTA record matching its Multi group."""
        path = Path(consensus_path)
        record_ids = []
        sequence_lengths = []
        current_sequence_length = None
        with path.open('r', encoding='utf-8') as handle:
            for line in handle:
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                if stripped_line.startswith('>'):
                    if current_sequence_length is not None:
                        sequence_lengths.append(current_sequence_length)
                    header_text = stripped_line[1:].strip()
                    record_ids.append(
                        header_text.split(maxsplit=1)[0]
                        if header_text else ''
                    )
                    current_sequence_length = 0
                else:
                    if current_sequence_length is None:
                        raise RuntimeError(
                            f"Invalid aligned consensus FASTA: {path}"
                        )
                    current_sequence_length += len(stripped_line)
        if current_sequence_length is not None:
            sequence_lengths.append(current_sequence_length)

        if (
            record_ids != [expected_multi_id]
            or len(sequence_lengths) != 1
            or sequence_lengths[0] == 0
        ):
            raise RuntimeError(
                f"Invalid aligned consensus FASTA: {path}"
            )

    def _prepare_multi_genes_aligned(self, aligned_dir, expected_multi_ids):
        """Validate and integrate all expected multi-gene consensus files."""
        self.logger.info("Integrating multi-gene alignment files....")
        
        aligned_dir_path = Path(aligned_dir)
        
        # アライメント済みファイルを探す
        aligned_files = sorted(
            aligned_dir_path.glob("*_aligned_consensus.fa")
        )
        observed_multi_ids = {
            aligned_file.name[:-len('_aligned_consensus.fa')]
            for aligned_file in aligned_files
        }
        missing_multi_ids = sorted(
            set(expected_multi_ids) - observed_multi_ids
        )
        unexpected_multi_ids = sorted(
            observed_multi_ids - set(expected_multi_ids)
        )
        if missing_multi_ids or unexpected_multi_ids:
            details = []
            if missing_multi_ids:
                details.append('missing=' + ','.join(missing_multi_ids))
            if unexpected_multi_ids:
                details.append('unexpected=' + ','.join(unexpected_multi_ids))
            raise RuntimeError(
                "Multi-gene alignment outputs do not match the expected "
                "groups: " + '; '.join(details)
            )
        if not aligned_files:
            raise RuntimeError(
                "No aligned consensus files were produced for expected "
                "Multi-gene groups"
            )
        
        self.logger.info(f"Number of aligned files found: {len(aligned_files)}")
        
        # 一時的なoneliner形式ファイルを作成するディレクトリ
        oneline_dir = os.path.join(self.intermediate_dir, 'multi_genes_aligned_oneline')
        Path(oneline_dir).mkdir(parents=True, exist_ok=True)
        
        # 各ファイルを1行形式に変換
        oneline_files = []
        for aligned_file in aligned_files:
            multi_id = aligned_file.name[:-len('_aligned_consensus.fa')]
            self._validate_multi_consensus_file(aligned_file, multi_id)
            # 出力ファイル名
            oneline_filename = aligned_file.name
            oneline_path = os.path.join(oneline_dir, oneline_filename)
            
            # 1行形式に変換
            self._convert_to_oneline_format(aligned_file, oneline_path)
            oneline_files.append(oneline_path)
            
            self.logger.debug(f"Single line format conversion completed: {oneline_filename}")
        
        # 全てのファイルを統合
        combined_path = os.path.join(self.intermediate_dir, 'multi_genes_aligned_oneline.fa')
        
        with open(combined_path, 'w') as f_out:
            for oneline_file in oneline_files:
                # ファイル名からベース名を抽出（拡張子を除去）
                base_name = Path(oneline_file).stem
                
                # ヘッダーとWindowsのCR対策を適用して統合
                self._merge_file_with_header_replacement(oneline_file, base_name, f_out)
        
        self.logger.info(f"Multi-gene alignment file integration complete.: {combined_path}")
        self.logger.info(f"Number of integrated files: {len(oneline_files)}")
        
        return combined_path
    
    def _convert_to_oneline_format(self, input_path, output_path):
        """Convert FASTA files to oneliner format (processing equivalent to Perl in FanFlow)."""
        with open(input_path, 'r') as f_in:
            with open(output_path, 'w') as f_out:
                sequence_buffer = ""
                header_line = ""
                
                for line in f_in:
                    line = line.rstrip('\r\n')  # 改行とCRを除去
                    
                    if line.startswith('>'):
                        # 前のシーケンスを出力
                        if header_line and sequence_buffer:
                            f_out.write(f"{header_line}\n{sequence_buffer}\n")
                        
                        header_line = line
                        sequence_buffer = ""
                    else:
                        # 配列行を連結
                        sequence_buffer += line
                
                # 最後のシーケンスを出力
                if header_line and sequence_buffer:
                    f_out.write(f"{header_line}\n{sequence_buffer}\n")
    
    def _merge_file_with_header_replacement(self, file_path, header_name, output_file):
        """Merge file with header replacement to output file"""
        with open(file_path, 'r') as f_in:
            for line in f_in:
                # 行末のCR(\r)を削除（Windowsの改行対策）
                line = line.rstrip('\r\n')
                
                if line.startswith('>'):
                    # ヘッダー行を差し替え
                    output_file.write(f">{header_name}\n")
                else:
                    # 配列行はそのまま出力
                    if line.strip():  # 空行でない場合のみ
                        output_file.write(f"{line}\n")
    
    def prepare_gene_sequences(self):
        """Prepare gene sequence files"""
        prepared_files = {}
        
        # Multi遺伝子の処理（アライメント済みファイルを統合）
        multi_genes_aligned_dir = os.path.join(self.intermediate_dir, 'multi_genes_aligned')
        multi_genes_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
        expected_multi_ids = self._expected_multi_gene_ids(multi_genes_path)
        aligned_dir_path = Path(multi_genes_aligned_dir)
        aligned_files_exist = bool(list(
            aligned_dir_path.glob('*_aligned_consensus.fa')
        )) if aligned_dir_path.is_dir() else False

        if expected_multi_ids:
            if not aligned_dir_path.is_dir():
                raise FileNotFoundError(
                    "Multi-gene alignment directory is missing for expected "
                    "groups: " + ','.join(sorted(expected_multi_ids))
                )
            multi_genes_oneline_path = self._prepare_multi_genes_aligned(
                multi_genes_aligned_dir, expected_multi_ids
            )
            prepared_files['multi_genes'] = multi_genes_oneline_path
            self.logger.info(
                "Prepare a file with multi-gene alignment.: "
                f"{multi_genes_oneline_path}"
            )
        elif aligned_files_exist:
            raise RuntimeError(
                "Aligned Multi-gene consensus files exist, but no expected "
                "Multi-gene groups were declared"
            )
        else:
            self.logger.info("No Multi-gene groups require alignment.")
        
        # Single遺伝子の処理
        single_genes_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
        expected_single_ids = self._expected_single_gene_ids(single_genes_path)
        if expected_single_ids:
            # Single遺伝子を1行形式に変換
            self._format_fasta_oneline(single_genes_path)
            
            # Convert degenerate bases to N (for CrisFlash compatibility)
            self._clean_degenerate_bases(single_genes_path)
            
            prepared_files['single_genes'] = single_genes_path
            self.logger.info(f"Prepare a single gene file.: {single_genes_path}")
        
        if not prepared_files:
            raise FileNotFoundError("No multi-gene or single-gene file was found.")
        
        return prepared_files
    
    def run_crisflash(self, genome_path, gene_path, gene_type):
        """CrisFlashによるガイド標的設計"""
        self.logger.info(f"Running CrisFlash on the {gene_type} genes...")
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, f'{gene_type}_targets.bed')
        workspace = DockerWorkspace(self.intermediate_dir, '/data')
        
        # CrisFlashコマンド
        cmd = [
            "docker", "run", "--rm",
            "-v", workspace.volume_spec,
            "-m", self.config['crispr']['memory'],
            "--cpus", str(self.config['crispr']['threads']),
            self.config['crispr']['docker_image'],
            "crisflash",
            "-g", workspace.container_path(genome_path),
            "-s", workspace.container_path(gene_path),
            "-p", self.config['crispr']['pam_pattern'],
            "-m", str(self.config['crispr']['max_mismatch']),
            "-t", str(self.config['crispr']['threads']),
            "-o", workspace.container_path(output_path)
        ]
        
        self.logger.debug(f"CrisFlash execution command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"CrisFlash execution error: {result.stderr}")
            raise RuntimeError(f"CrisFlash execution failed. ({gene_type})")
        
        self.logger.info(f"Guide target design for {gene_type} is complete.: {output_path}")
        
        # 結果の統計を取得
        if os.path.exists(output_path):
            line_count = self._count_lines(output_path)
            self.logger.info(f"Number of guide targets for {gene_type}: {line_count}")
        
        return output_path
    
    def _count_lines(self, file_path):
        """Count lines in file"""
        try:
            with open(file_path, 'r') as f:
                return sum(1 for _ in f)
        except Exception:
            return 0

    @staticmethod
    def _is_site_specific_target(target_field):
        """Return whether a CrisFlash target has MM0=1 and MM1-2=0."""
        parts = target_field.strip().split('/')
        if len(parts) != 3 or not parts[0]:
            return False

        try:
            mm0 = int(parts[1])
            mm1_2 = int(parts[2])
        except ValueError:
            return False

        return mm0 == 1 and mm1_2 == 0
    
    def filter_site_specific_targets(self, bed_path, gene_type):
        """Extract site-specific targets"""
        self.logger.info(f"Extracting site-specific targets for {gene_type}...")
        
        # Output file paths
        specific_path = os.path.join(self.intermediate_dir, f'{gene_type}_targets_specific.bed')
        cleaned_path = os.path.join(self.intermediate_dir, f'{gene_type}_targets_specific_clean.bed')
        
        # The fourth column uses the format <sequence>/<MM0>/<MM1-2>.
        with open(bed_path, 'r') as f_in:
            with open(specific_path, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    if len(fields) >= 4 and self._is_site_specific_target(fields[3]):
                        f_out.write(line)
        
        # Remove trailing whitespace from every field.
        with open(specific_path, 'r') as f_in:
            with open(cleaned_path, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    cleaned_fields = [field.rstrip() for field in fields]
                    f_out.write('\t'.join(cleaned_fields) + '\n')
        
        # Log the number of retained site-specific targets.
        specific_count = self._count_lines(cleaned_path)
        self.logger.info(f"Number of site-specific targets for {gene_type}: {specific_count}")
        
        return cleaned_path
    
    def merge_with_annotation(self, targets_path, gene_type):
        """Merge guide target information with annotation"""
        self.logger.info(f"{gene_type}のMerge guide target information with annotation中...")
        
        # アノテーション分類ファイル
        annotation_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
        
        if not os.path.exists(annotation_path):
            self.logger.warning("Annotation classification file not found.")
            return None
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, f'{gene_type}_annotation_with_guides.tsv')
        
        # BEDファイルからターゲット情報を読み込み
        targets_map = {}
        
        with open(targets_path, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6:
                    # 1列目からキーを抽出
                    gene_key = fields[0]
                    
                    # Multi遺伝子の場合は"_aligned_consensus"を除去
                    if gene_type == 'multi_genes':
                        gene_key = gene_key.replace('_aligned_consensus', '')
                    
                    # BEDの1-6列目を結合
                    bed_info = '/'.join(fields[:6])
                    
                    if gene_key not in targets_map:
                        targets_map[gene_key] = []
                    targets_map[gene_key].append(bed_info)
        
        # アノテーションファイルと結合
        with open(annotation_path, 'r') as f_in:
            with open(output_path, 'w') as f_out:
                header = f_in.readline().strip()
                f_out.write(f"{header}\tguide_targets\n")
                
                for line in f_in:
                    annotation_line = line.rstrip('\r\n')
                    fields = annotation_line.split('\t')
                    if len(fields) >= 10:
                        unique_id = fields[9]  # 10列目がunique_id
                        classification = fields[8]  # 9列目がclassification
                        
                        # 対象のgene_typeかチェック
                        target_type = "Multi" if gene_type == 'multi_genes' else "Single"
                        
                        if classification == target_type:
                            # ガイド標的情報を取得
                            guide_targets = targets_map.get(unique_id, [])
                            guide_targets_str = ','.join(guide_targets) if guide_targets else ''
                            
                            f_out.write(f"{annotation_line}\t{guide_targets_str}\n")
        
        self.logger.info(f"Save the annotation-guide binding results for {gene_type}.: {output_path}")
        return output_path
    
    def create_final_combined_table(self, multi_path, single_path):
        """Create final combined table"""
        self.logger.info("Create final combined table中...")
        
        output_path = os.path.join(self.intermediate_dir, 'annotation_with_guides_final.tsv')
        
        # ヘッダーを書き込み
        with open(output_path, 'w') as f_out:
            header_written = False
            
            # Multi遺伝子の結果を追加
            if multi_path and os.path.exists(multi_path):
                with open(multi_path, 'r') as f_in:
                    header = f_in.readline().strip()
                    if not header_written:
                        f_out.write(header + '\n')
                        header_written = True
                    
                    for line in f_in:
                        f_out.write(line)
            
            # Single遺伝子の結果を追加
            if single_path and os.path.exists(single_path):
                with open(single_path, 'r') as f_in:
                    header = f_in.readline().strip()
                    if not header_written:
                        f_out.write(header + '\n')
                        header_written = True
                    
                    for line in f_in:
                        f_out.write(line)
        
        self.logger.info(f"Save the final joined table: {output_path}")
        return output_path
    
    def run(self):
        """CRISPRガイド設計の全体実行"""
        try:
            # ゲノムアセンブリの準備
            genome_path = self.prepare_genome_assembly()
            
            # Prepare gene sequence files
            gene_files = self.prepare_gene_sequences()
            
            results = {
                'genome_path': genome_path,
                'gene_files': gene_files,
                'target_results': {},
                'annotation_results': {}
            }
            
            # 各遺伝子タイプに対してガイド設計を実行
            for gene_type, gene_path in gene_files.items():
                # CrisFlashによるガイド標的設計
                targets_path = self.run_crisflash(genome_path, gene_path, gene_type)
                
                # サイト特異的標的の抽出
                specific_path = self.filter_site_specific_targets(targets_path, gene_type)
                
                # アノテーションとの結合
                annotation_with_guides = self.merge_with_annotation(specific_path, gene_type)
                
                results['target_results'][gene_type] = {
                    'all_targets': targets_path,
                    'specific_targets': specific_path,
                    'annotation_with_guides': annotation_with_guides
                }
            
            # 最終的な結合テーブルの作成
            multi_annotation = results['target_results'].get('multi_genes', {}).get('annotation_with_guides')
            single_annotation = results['target_results'].get('single_genes', {}).get('annotation_with_guides')
            
            final_table = self.create_final_combined_table(multi_annotation, single_annotation)
            results['final_table'] = final_table
            
            # 結果をファイルに保存
            results_file = os.path.join(self.intermediate_dir, 'guide_design_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== CRISPR guided design results ===\n\n")
                f.write(f"Genome assembly: {genome_path}\n")
                f.write(f"Gene file: {gene_files}\n\n")
                
                for gene_type, result in results['target_results'].items():
                    f.write(f"=== {gene_type} ===\n")
                    f.write(f"All guide targets: {result['all_targets']}\n")
                    f.write(f"Specific target: {result['specific_targets']}\n")
                    f.write(f"Annotation related: {result['annotation_with_guides']}\n\n")
                
                f.write(f"Final Join Table: {final_table}\n")
            
            self.logger.info(f"Save CRISPR guided design results: {results_file}")
            self.logger.info("The CRISPR guided design process has been completed.")
            
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during the CRISPR guided design process.: {e}")
            raise
