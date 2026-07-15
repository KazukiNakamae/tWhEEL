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
        genome_path = self.config['input']['genome_assembly']
        
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
    
    def _prepare_multi_genes_aligned(self, aligned_dir):
        """Integrate multi-gene alignment files"""
        self.logger.info("Integrating multi-gene alignment files....")
        
        aligned_dir_path = Path(aligned_dir)
        
        # アライメント済みファイルを探す
        aligned_files = list(aligned_dir_path.glob("*_aligned_consensus.fa"))
        
        if not aligned_files:
            self.logger.warning("Aligned files not found.")
            return None
        
        self.logger.info(f"Number of aligned files found: {len(aligned_files)}")
        
        # 一時的なoneliner形式ファイルを作成するディレクトリ
        oneline_dir = os.path.join(self.intermediate_dir, 'multi_genes_aligned_oneline')
        Path(oneline_dir).mkdir(parents=True, exist_ok=True)
        
        # 各ファイルを1行形式に変換
        oneline_files = []
        for aligned_file in aligned_files:
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
        if os.path.exists(multi_genes_aligned_dir):
            multi_genes_oneline_path = self._prepare_multi_genes_aligned(multi_genes_aligned_dir)
            if multi_genes_oneline_path:
                prepared_files['multi_genes'] = multi_genes_oneline_path
                self.logger.info(f"Prepare a file with multi-gene alignment.: {multi_genes_oneline_path}")
        else:
            # フォールバック: 従来のmulti_genes.faを使用
            multi_genes_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            if os.path.exists(multi_genes_path):
                self._format_fasta_oneline(multi_genes_path)
                prepared_files['multi_genes'] = multi_genes_path
                self.logger.info(f"Prepare a multi-gene file (non-aligned).: {multi_genes_path}")
        
        # Single遺伝子の処理
        single_genes_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
        if os.path.exists(single_genes_path):
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
        
        # CrisFlashコマンド
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.getcwd()}:/data",
            "-m", self.config['crispr']['memory'],
            "--cpus", str(self.config['crispr']['threads']),
            self.config['crispr']['docker_image'],
            "crisflash",
            "-g", f"/data/{genome_path}",
            "-s", f"/data/{gene_path}",
            "-p", self.config['crispr']['pam_pattern'],
            "-m", str(self.config['crispr']['max_mismatch']),
            "-t", str(self.config['crispr']['threads']),
            "-o", f"/data/{output_path}"
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
    
    def filter_site_specific_targets(self, bed_path, gene_type):
        """Extract site-specific targets"""
        self.logger.info(f"Extracting site-specific targets for {gene_type}...")
        
        # 出力ファイルパス
        specific_path = os.path.join(self.intermediate_dir, f'{gene_type}_targets_specific.bed')
        cleaned_path = os.path.join(self.intermediate_dir, f'{gene_type}_targets_specific_clean.bed')
        
        # Extract site-specific targets（4列目が/1/パターン）
        with open(bed_path, 'r') as f_in:
            with open(specific_path, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    if len(fields) >= 4:
                        # 4列目が/1/<数値>パターンかチェック
                        if '/1/' in fields[3]:
                            f_out.write(line)
        
        # 末尾の空白を除去
        with open(specific_path, 'r') as f_in:
            with open(cleaned_path, 'w') as f_out:
                for line in f_in:
                    # 各フィールドの末尾空白を除去
                    fields = line.strip().split('\t')
                    cleaned_fields = [field.rstrip() for field in fields]
                    f_out.write('\t'.join(cleaned_fields) + '\n')
        
        # 統計情報をログ出力
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
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        unique_id = fields[9]  # 10列目がunique_id
                        classification = fields[8]  # 9列目がclassification
                        
                        # 対象のgene_typeかチェック
                        target_type = "Multi" if gene_type == 'multi_genes' else "Single"
                        
                        if classification == target_type:
                            # ガイド標的情報を取得
                            guide_targets = targets_map.get(unique_id, [])
                            guide_targets_str = ','.join(guide_targets) if guide_targets else ''
                            
                            f_out.write(f"{line.strip()}\t{guide_targets_str}\n")
        
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
