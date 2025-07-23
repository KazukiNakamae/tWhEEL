#!/usr/bin/env python3
"""
Report Generation Module
Text summary and HTML interactive report generation
"""

import os
import json
import logging
import sqlite3
import csv
import requests
import time
import itertools
from datetime import datetime
from pathlib import Path
import shutil

# CSVフィールドサイズ制限を増加
csv.field_size_limit(1 << 30)

class ReportGenerator:
    """Report generation processing class"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.output_dir = self.config['output']['base_dir']
        self.web_data_dir = os.path.join(self.output_dir, 'web_data')
        
    def collect_all_results(self):
        """Collect results from all steps"""
        self.logger.info("Collecting results from all steps...")
        
        results = {
            'timestamp': datetime.now().isoformat(),
            'config': self.config,
            'files': {},
            'statistics': {},
            'data': {}
        }
        
        # Collect result files from each step
        result_files = {
            'gene_prediction': 'gene_prediction_results.txt',
            'annotation': 'annotation_results.txt',
            'guide_design': 'guide_design_results.txt',
            'go_analysis': 'go_analysis_results.txt'
        }
        
        for step, filename in result_files.items():
            filepath = os.path.join(self.intermediate_dir, filename)
            if os.path.exists(filepath):
                results['files'][step] = filepath
        
        # Collect main data files
        data_files = {
            'transcriptome': 'transcriptome.fa',
            'protein_sequences': 'query.pep',
            'annotation_table': 'annotation_classified.txt',
            'final_annotation_guides': 'annotation_with_guides_final.tsv',
            'go_annotated': 'annotation_with_guides_go.tsv',
            'go_kfi_table': 'go_kfi_table.csv',
            'go_kfi_stat': 'go_kfi_stat.csv',
            'fully_targetable_gos': 'go_kfi_table_fully_targetable.csv',
            'go_reverse_lookup': 'go_reverse_lookup.tsv'
        }
        
        for name, filename in data_files.items():
            filepath = os.path.join(self.intermediate_dir, filename)
            if os.path.exists(filepath):
                results['data'][name] = filepath
        
        # Collect statistics (pass results for new guide RNA statistics)
        results['statistics'] = self._collect_statistics(results)
        
        self.logger.info("Results collection completed")
        return results
    
    def _collect_statistics(self, results):
        """Collect statistics information"""
        stats = {
            'gene_prediction': {},
            'annotation': {},
            'guide_design': {},
            'go_analysis': {}
        }
        
        # Gene prediction statistics
        try:
            # Count total transcripts from transcriptome.fa
            transcriptome_path = self.config['input'].get('transcriptome')
            if transcriptome_path and os.path.exists(transcriptome_path):
                stats['gene_prediction']['transcript_count'] = self._count_fasta_sequences(transcriptome_path)
            else:
                stats['gene_prediction']['transcript_count'] = 0
            
            protein_path = os.path.join(self.intermediate_dir, 'query.pep')
            if os.path.exists(protein_path):
                protein_data = self._parse_protein_file(protein_path)
                stats['gene_prediction']['protein_count'] = len(protein_data['proteins'])
                stats['gene_prediction']['protein_stats'] = protein_data['stats']
                stats['gene_prediction']['proteins'] = protein_data['proteins']
        except Exception as e:
            self.logger.warning(f"Gene prediction statistics collection error: {e}")
            stats['gene_prediction']['transcript_count'] = 0
        
        # Annotation statistics
        try:
            annotation_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
            if os.path.exists(annotation_path):
                multi_count, single_count = self._count_gene_classifications(annotation_path)
                stats['annotation']['multi_genes'] = multi_count
                stats['annotation']['single_genes'] = single_count
                stats['annotation']['total_genes'] = multi_count + single_count
        except Exception as e:
            self.logger.warning(f"Annotation statistics collection error: {e}")
        
        # Guide design statistics
        try:
            multi_targets = os.path.join(self.intermediate_dir, 'multi_genes_targets_specific_clean.bed')
            single_targets = os.path.join(self.intermediate_dir, 'single_genes_targets_specific_clean.bed')
            
            if os.path.exists(multi_targets):
                stats['guide_design']['multi_targets'] = self._count_lines(multi_targets)
            if os.path.exists(single_targets):
                stats['guide_design']['single_targets'] = self._count_lines(single_targets)
            
            # 新しい詳細なガイドRNA統計を追加
            guide_rna_stats = self._calculate_guide_rna_statistics(results)
            stats['guide_design'].update(guide_rna_stats)
            
        except Exception as e:
            self.logger.warning(f"Guide design statistics collection error: {e}")
        
        # GO analysis statistics
        try:
            kfi_stat_path = os.path.join(self.intermediate_dir, 'go_kfi_stat.csv')
            if os.path.exists(kfi_stat_path):
                stats['go_analysis'] = self._parse_go_kfi_stats(kfi_stat_path)
        except Exception as e:
            self.logger.warning(f"GO analysis statistics collection error: {e}")
        
        return stats
    
    def _count_fasta_sequences(self, fasta_path):
        """Count number of sequences in FASTA file"""
        count = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    
    def _count_gene_classifications(self, annotation_path):
        """Count gene classifications"""
        multi_count = 0
        single_count = 0
        
        with open(annotation_path, 'r') as f:
            header = f.readline()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    classification = fields[8]
                    if classification == 'Multi':
                        multi_count += 1
                    elif classification == 'Single':
                        single_count += 1
        
        return multi_count, single_count
    
    def _count_lines(self, filepath):
        """Count lines in file"""
        try:
            with open(filepath, 'r') as f:
                return sum(1 for _ in f)
        except Exception:
            return 0
    
    def _calculate_guide_rna_statistics(self, results):
        """Calculate detailed guide RNA statistics"""
        stats = {
            'total_guides': 0,
            'site_specific_guides': 0,
            'site_specific_multi_target': 0,
            'site_specific_single_target': 0
        }
        
        annotation_path = results['data'].get('final_annotation_guides')
        if not annotation_path or not os.path.exists(annotation_path):
            return stats
        
        try:
            with open(annotation_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)
                
                for row in reader:
                    if len(row) >= 11 and row[-1]:  # guide_targets列が存在
                        guide_targets = row[-1]
                        
                        if guide_targets and guide_targets.strip():
                            targets = guide_targets.split(',')
                            for target in targets:
                                if target.strip():
                                    guide_parts = target.strip().split('/')
                                    if len(guide_parts) >= 8:
                                        # Total Guides
                                        stats['total_guides'] += 1
                                        
                                        # MM0とMM1-2の値を取得
                                        mm0 = int(guide_parts[4]) if guide_parts[4].isdigit() else 0
                                        mm1_2 = int(guide_parts[5]) if guide_parts[5].isdigit() else 0
                                        reference = guide_parts[0]
                                        
                                        # Site-Specific判定 (MM0=1 and MM1-2=0 → off_target_score=1)
                                        if mm0 == 1 and mm1_2 == 0:
                                            stats['site_specific_guides'] += 1
                                            
                                            # Multi/Single判定
                                            if reference.startswith('Multi_'):
                                                stats['site_specific_multi_target'] += 1
                                            elif reference.startswith('Single_'):
                                                stats['site_specific_single_target'] += 1
            
            return stats
            
        except Exception as e:
            self.logger.warning(f"Guide RNA statistics calculation error: {e}")
            return stats
    
    def _parse_go_kfi_stats(self, kfi_stat_path):
        """Parse GO-KFI statistics"""
        stats = {}
        
        try:
            with open(kfi_stat_path, 'r') as f:
                content = f.read()
                
            import re
            
            # Extract basic statistics
            pattern = r'F:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,M:(\d+\.?\d*)%\],P:(\d+\.?\d*)%,N:(\d+\.?\d*)%,n:(\d+)'
            match = re.search(pattern, content)
            
            if match:
                stats['fully_targetable_percent'] = float(match.group(1))
                stats['simply_targetable_percent'] = float(match.group(2))
                stats['multiply_targetable_percent'] = float(match.group(3))
                stats['partially_targetable_percent'] = float(match.group(4))
                stats['not_targetable_percent'] = float(match.group(5))
                stats['total_gos'] = int(match.group(6))
                
                # Extract detailed counts for chart visualization
                count_patterns = [
                    (r'(\d+)\s+Fully targetable GOs \(F\)', 'fully_targetable_count'),
                    (r'(\d+)\s+Fully and simply targetable GOs \(S\)', 'simply_targetable_count'),
                    (r'(\d+)\s+Fully and multiply targetable GOs \(M\)', 'multiply_targetable_count'),
                    (r'(\d+)\s+Partially targetable GOs \(P\)', 'partially_targetable_count'),
                    (r'(\d+)\s+Not targetable GOs \(N\)', 'not_targetable_count')
                ]
                
                for pattern, key in count_patterns:
                    count_match = re.search(pattern, content)
                    if count_match:
                        stats[key] = int(count_match.group(1))
                
                # Create chart data for JavaScript
                stats['chart_data'] = {
                    'fully_targetable': {
                        'percent': stats['fully_targetable_percent'],
                        'count': stats.get('fully_targetable_count', 0),
                        'simply_percent': stats['simply_targetable_percent'],
                        'multiply_percent': stats['multiply_targetable_percent']
                    },
                    'partially_targetable': {
                        'percent': stats['partially_targetable_percent'],
                        'count': stats.get('partially_targetable_count', 0)
                    },
                    'not_targetable': {
                        'percent': stats['not_targetable_percent'],
                        'count': stats.get('not_targetable_count', 0)
                    },
                    'total': stats['total_gos']
                }
                
        except Exception as e:
            self.logger.warning(f"GO-KFI statistics parsing error: {e}")
        
        return stats
    
    def _parse_protein_file(self, protein_path):
        """Parse protein FASTA file and extract summary information"""
        self.logger.info(f"Parsing protein file: {protein_path}")
        
        proteins = []
        
        try:
            with open(protein_path, 'r', encoding='utf-8') as f:
                current_header = None
                current_sequence = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Process previous protein if exists
                        if current_header is not None:
                            sequence = ''.join(current_sequence)
                            protein_info = self._process_protein(current_header, sequence)
                            proteins.append(protein_info)
                        
                        # Start new protein
                        current_header = line[1:]  # Remove '>' character
                        current_sequence = []
                    elif line and current_header is not None:
                        # Add sequence line
                        current_sequence.append(line)
                
                # Process last protein
                if current_header is not None:
                    sequence = ''.join(current_sequence)
                    protein_info = self._process_protein(current_header, sequence)
                    proteins.append(protein_info)
        
        except Exception as e:
            self.logger.error(f"Error parsing protein file: {e}")
            return {'proteins': [], 'stats': {}}
        
        # Calculate summary statistics
        stats = self._calculate_protein_summary_stats(proteins)
        
        self.logger.info(f"Parsed {len(proteins)} proteins from {protein_path}")
        
        return {
            'proteins': proteins,
            'stats': stats
        }
    
    def _process_protein(self, header, sequence):
        """Process individual protein and extract information"""
        length = len(sequence)
        category = self._categorize_protein_length(length)
        mass_kda = self._calculate_protein_mass(sequence)
        
        return {
            'id': header,
            'length': length,
            'category': category,
            'mass_kda': round(mass_kda, 2)
        }
    
    def _categorize_protein_length(self, length):
        """Categorize protein by length according to user specifications"""
        if length <= 100:
            return 'short'
        elif length >= 1000:
            return 'long'
        else:
            return 'medium'
    
    def _calculate_protein_mass(self, sequence):
        """Calculate estimated atomic mass of protein in kDa"""
        # Average molecular weights of amino acids (Da)
        aa_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
            'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
        }
        
        total_mass = 0.0
        valid_residues = 0
        
        for aa in sequence.upper():
            if aa in aa_weights:
                total_mass += aa_weights[aa]
                valid_residues += 1
        
        if valid_residues == 0:
            return 0.0
        
        # Subtract water molecules for peptide bonds (n-1 peptide bonds for n amino acids)
        if valid_residues > 1:
            total_mass -= (valid_residues - 1) * 18.015  # H2O molecular weight
        
        # Convert to kDa
        return total_mass / 1000.0
    
    def _calculate_protein_summary_stats(self, proteins):
        """Calculate summary statistics for proteins"""
        if not proteins:
            return {
                'total_count': 0,
                'avg_length': 0,
                'max_length': 0,
                'avg_mass_kda': 0,
                'max_mass_kda': 0,
                'categories': {'short': 0, 'medium': 0, 'long': 0}
            }
        
        lengths = [p['length'] for p in proteins]
        masses = [p['mass_kda'] for p in proteins]
        
        # Count categories
        categories = {'short': 0, 'medium': 0, 'long': 0}
        for protein in proteins:
            categories[protein['category']] += 1
        
        return {
            'total_count': len(proteins),
            'avg_length': round(sum(lengths) / len(lengths), 1),
            'max_length': max(lengths),
            'avg_mass_kda': round(sum(masses) / len(masses), 2),
            'max_mass_kda': round(max(masses), 2),
            'categories': categories
        }
    
    def generate_text_summary(self, results):
        """Generate text summary"""
        self.logger.info("Generating text summary...")
        
        summary_path = os.path.join(self.output_dir, 'summary.txt')
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("tWhEEL - Transcriptome Functional Analysis Pipeline\n")
            f.write("Execution Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Execution Date: {results['timestamp']}\n\n")
            
            # Processing statistics
            f.write("=== Processing Statistics ===\n")
            f.write("-" * 30 + "\n")
            
            stats = results['statistics']
            
            # Gene prediction statistics
            if 'gene_prediction' in stats:
                f.write(f"Protein Annotated Contigs: {stats['gene_prediction'].get('protein_count', 'N/A')}\n")
            
            # Annotation statistics
            if 'annotation' in stats:
                f.write(f"Multi (Redundant) Contigs</div>: {stats['annotation'].get('multi_genes', 'N/A')}\n")
                f.write(f"Single (Non-Redundant) Contigs: {stats['annotation'].get('single_genes', 'N/A')}\n")
                f.write(f"Total Annotated Contigs: {stats['annotation'].get('total_genes', 'N/A')}\n")
            
            # Guide design statistics
            if 'guide_design' in stats:
                f.write(f"Multi (Redundant) Contig Targets: {stats['guide_design'].get('multi_targets', 'N/A')}\n")
                f.write(f"Single (Non-Redundant) Contig Targets: {stats['guide_design'].get('single_targets', 'N/A')}\n")
            
            # GO analysis statistics
            if 'go_analysis' in stats:
                go_stats = stats['go_analysis']
                f.write(f"Total GO terms: {go_stats.get('total_gos', 'N/A')}\n")
                f.write(f"Fully targetable GO percentage: {go_stats.get('fully_targetable_percent', 'N/A')}%\n")
                f.write(f"Partially targetable GO percentage: {go_stats.get('partially_targetable_percent', 'N/A')}%\n")
                f.write(f"Not targetable GO percentage: {go_stats.get('not_targetable_percent', 'N/A')}%\n")
            
            f.write("\n")
            
            # Configuration information
            f.write("=== Configuration Settings ===\n")
            f.write("-" * 30 + "\n")
            f.write(f"Input files:\n")
            f.write(f"  - Genome assembly: {results['config']['input'].get('genome_assembly', 'N/A')}\n")
            f.write(f"  - Transcriptome: {results['config']['input'].get('transcriptome', 'N/A')}\n")
            f.write(f"  - UniProt database: {results['config']['input'].get('uniprot_db', 'N/A')}\n")
            f.write(f"PAM pattern: {results['config']['crispr'].get('pam_pattern', 'N/A')}\n")
            f.write(f"Max mismatch: {results['config']['crispr'].get('max_mismatch', 'N/A')}\n")
            f.write(f"BUSCO lineage: {results['config']['busco'].get('lineage', 'N/A')}\n")
            
            f.write("\n")
            
            # Output files
            f.write("=== Output Files ===\n")
            f.write("-" * 30 + "\n")
            f.write(f"HTML report: {os.path.join(self.output_dir, 'report.html')}\n")
            f.write(f"Web data: {self.web_data_dir}\n")
            f.write(f"Intermediate files: {self.intermediate_dir}\n")
            
            f.write("\n")
            f.write("For detailed analysis results, please refer to the HTML report (report.html).\n")
        
        self.logger.info(f"Text summary generated: {summary_path}")
        return summary_path
    
    def _check_database_exists(self, db_path):
        """Check if database exists and is valid"""
        if not os.path.exists(db_path):
            self.logger.info("Database file does not exist.")
            return False
        
        # Check file size
        if os.path.getsize(db_path) == 0:
            self.logger.warning("Database file exists but is empty.")
            return False
        
        try:
            # Try to connect to database
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Check database structure
            if not self._verify_database_structure(cursor):
                conn.close()
                return False
            
            # Check if database has data
            if not self._check_database_data(cursor):
                conn.close()
                return False
            
            conn.close()
            self.logger.info(f"既存の有効なデータベースを検出しました: {db_path}")
            return True
            
        except sqlite3.Error as e:
            self.logger.warning(f"Database connection error: {e}")
            return False
    
    def _verify_database_structure(self, cursor):
        """Verify that all required tables exist"""
        required_tables = ['genes', 'guide_rnas', 'go_terms', 'gene_go_relations']
        
        try:
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            existing_tables = [row[0] for row in cursor.fetchall()]
            
            missing_tables = [table for table in required_tables if table not in existing_tables]
            
            if missing_tables:
                self.logger.warning(f"Missing tables: {missing_tables}")
                return False
            
            self.logger.info("データベーステーブル構造を確認しました。")
            return True
            
        except sqlite3.Error as e:
            self.logger.warning(f"Database structure verification error: {e}")
            return False
    
    def _check_database_data(self, cursor):
        """Check if database contains data"""
        try:
            # Check if genes table has data
            cursor.execute("SELECT COUNT(*) FROM genes")
            gene_count = cursor.fetchone()[0]
            
            if gene_count == 0:
                self.logger.warning("Database tables exist but genes table is empty.")
                return False
            
            self.logger.info(f"データベースに {gene_count} 件の遺伝子データが存在します。")
            return True
            
        except sqlite3.Error as e:
            self.logger.warning(f"Database data check error: {e}")
            return False

    def prepare_web_data(self, results):
        """Prepare data for web application in JSON and SQLite formats"""
        self.logger.info("Preparing web data...")
        
        # Create web data directory if it doesn't exist
        Path(self.web_data_dir).mkdir(parents=True, exist_ok=True)
        
        # Create SQLite database
        db_path = os.path.join(self.web_data_dir, 'twheel_data.db')
        
        # Check for force rebuild flag in config
        force_rebuild = self.config.get('database', {}).get('force_rebuild', False)
        
        if not force_rebuild and self._check_database_exists(db_path):
            self.logger.info("既存のデータベースを使用します。データベース作成をスキップします。")
        else:
            if force_rebuild:
                self.logger.info("強制再生成フラグが設定されています。データベースを新規作成します。")
            else:
                self.logger.info("データベースが存在しないか無効です。新規作成します。")
            
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create tables
            self._create_database_tables(cursor)
            
            # Insert data
            self._insert_gene_data(cursor, results)
            self._insert_guide_data(cursor, results)
            self._insert_go_data(cursor, results)
            
            conn.commit()
            conn.close()
            
            self.logger.info(f"データベース作成完了: {db_path}")
        
        # Save basic statistics in JSON format
        json_data = {
            'timestamp': results['timestamp'],
            'statistics': results['statistics'],
            'config': {
                'input': results['config']['input'],
                'crispr': results['config']['crispr'],
                'busco': results['config']['busco']
            }
        }
        
        json_path = os.path.join(self.web_data_dir, 'summary.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, ensure_ascii=False, indent=2)
        
        self.logger.info(f"Web data prepared: {db_path}, {json_path}")
        return db_path, json_path
    
    def _create_database_tables(self, cursor):
        """Create database tables"""
        
        # Gene table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS genes (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                query_id TEXT,
                subject_id TEXT,
                gene_name TEXT,
                classification TEXT,
                unique_id TEXT,
                identity REAL,
                evalue REAL,
                description TEXT,
                guide_count INTEGER DEFAULT 0
            )
        """)
        
        # Guide RNA table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS guide_rnas (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_unique_id TEXT,
                chromosome TEXT,
                start INTEGER,
                end INTEGER,
                sequence TEXT,
                strand TEXT,
                score REAL,
                pam_sequence TEXT,
                gc_content REAL,
                FOREIGN KEY (gene_unique_id) REFERENCES genes(unique_id)
            )
        """)
        
        # GO term table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS go_terms (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                go_id TEXT,
                go_name TEXT,
                go_definition TEXT,
                go_namespace TEXT,
                total_genes INTEGER,
                genes_with_guides INTEGER,
                kfi_score REAL,
                targetable_type TEXT,
                annotated_genes TEXT,
                annotated_genes_full TEXT
            )
        """)
        
        # Gene-GO relation table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS gene_go_relations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_unique_id TEXT,
                go_id TEXT,
                FOREIGN KEY (gene_unique_id) REFERENCES genes(unique_id),
                FOREIGN KEY (go_id) REFERENCES go_terms(go_id)
            )
        """)
    
    def _insert_gene_data(self, cursor, results):
        """Insert gene data"""
        
        annotation_path = results['data'].get('annotation_table')
        if not annotation_path or not os.path.exists(annotation_path):
            return
        
        with open(annotation_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            for row in reader:
                if len(row) >= 10:
                    cursor.execute("""
                        INSERT INTO genes (
                            query_id, subject_id, gene_name, classification,
                            unique_id, identity, evalue, description
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        row[0], row[1], row[7], row[8], row[9],
                        float(row[2]) if row[2] else 0.0,
                        float(row[4]) if row[4] else 0.0,
                        row[6]
                    ))
    
    def _insert_guide_data(self, cursor, results):
        """Insert guide RNA data and update guide counts"""
        
        final_table_path = results['data'].get('final_annotation_guides')
        if not final_table_path or not os.path.exists(final_table_path):
            return
        
        with open(final_table_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # Track guide counts for each gene
            guide_counts = {}
            
            for row in reader:
                if len(row) >= 11 and row[-1]:  # guide_targets column exists
                    query_id = row[0] if len(row) > 0 else None  # Get query_id
                    unique_id = row[9]
                    guide_targets = row[-1]
                    
                    # Initialize guide count for this gene
                    if query_id not in guide_counts:
                        guide_counts[query_id] = 0
                    
                    # Parse guide targets
                    if guide_targets and guide_targets.strip():
                        guides = guide_targets.split(',')
                        for guide in guides:
                            if guide.strip():
                                guide_parts = guide.strip().split('/')
                                if len(guide_parts) >= 6:
                                    cursor.execute("""
                                        INSERT INTO guide_rnas (
                                            gene_unique_id, chromosome, start, end,
                                            sequence, strand
                                        ) VALUES (?, ?, ?, ?, ?, ?)
                                    """, (
                                        unique_id,
                                        guide_parts[0],
                                        int(guide_parts[1]) if guide_parts[1].isdigit() else 0,
                                        int(guide_parts[2]) if guide_parts[2].isdigit() else 0,
                                        guide_parts[3],
                                        guide_parts[5]
                                    ))
                                    
                                    # Increment guide count
                                    guide_counts[query_id] += 1
            
            # Update guide counts in genes table
            for query_id, count in guide_counts.items():
                if query_id:
                    cursor.execute("""
                        UPDATE genes SET guide_count = ? WHERE query_id = ?
                    """, (count, query_id))
    
    def _insert_go_data(self, cursor, results):
        """Insert GO term data with enhanced information"""
        
        kfi_table_path = results['data'].get('go_kfi_table')
        if not kfi_table_path or not os.path.exists(kfi_table_path):
            return
        
        # Get GO IDs from kfi_table
        go_ids = []
        go_data_map = {}
        
        with open(kfi_table_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            
            for row in reader:
                if len(row) >= 6:
                    go_id = row[0]
                    total_genes = int(row[1]) if row[1] else 0
                    genes_with_guides = int(row[2]) if row[2] else 0
                    kfi_score = float(row[4]) if row[4] else 0.0
                    targetable_type = row[5] if len(row) > 5 else 'Not'
                    
                    go_ids.append(go_id)
                    go_data_map[go_id] = {
                        'total_genes': total_genes,
                        'genes_with_guides': genes_with_guides,
                        'kfi_score': kfi_score,
                        'targetable_type': targetable_type
                    }
        
        # Fetch GO terms and definitions from API
        try:
            go_terms_info = self._fetch_go_terms_from_api(go_ids)
        except Exception as e:
            self.logger.warning(f"Failed to fetch GO terms from API: {e}")
            go_terms_info = {}
        
        # Parse GO-annotated genes
        try:
            go_genes_map = self._parse_go_annotated_genes(results)
        except Exception as e:
            self.logger.warning(f"Failed to parse GO-annotated genes: {e}")
            go_genes_map = {}
        
        # Insert data into database
        for go_id in go_ids:
            go_info = go_terms_info.get(go_id, {})
            go_data = go_data_map.get(go_id, {})
            annotated_genes = go_genes_map.get(go_id, [])
            
            # Format annotated genes (first 5 + "...and X more")
            formatted_genes = self._format_annotated_genes_list(annotated_genes)
            
            # Create full genes list for expandable functionality
            full_genes_list = ", ".join(annotated_genes) if annotated_genes else ""
            
            cursor.execute("""
                INSERT INTO go_terms (
                    go_id, go_name, go_definition, total_genes, 
                    genes_with_guides, kfi_score, targetable_type, annotated_genes, annotated_genes_full
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                go_id,
                go_info.get('name', ''),
                go_info.get('definition', ''),
                go_data.get('total_genes', 0),
                go_data.get('genes_with_guides', 0),
                go_data.get('kfi_score', 0.0),
                go_data.get('targetable_type', 'Not'),
                formatted_genes,
                full_genes_list
            ))
    
    def _fetch_go_terms_from_api(self, go_ids):
        """Fetch GO terms and definitions from QuickGO API"""
        
        def chunks(iterable, size):
            it = iter(iterable)
            while (batch := list(itertools.islice(it, size))):
                yield batch
        
        def fetch(ids, tries=3, timeout=30):
            url_template = ("https://www.ebi.ac.uk/QuickGO/services/"
                           "ontology/go/terms/{ids}?fields=id,name,definition")
            url = url_template.format(ids=",".join(ids))
            
            for attempt in range(1, tries + 1):
                try:
                    self.logger.info(f"Fetching GO terms from API (attempt {attempt}/{tries}): {len(ids)} terms")
                    response = requests.get(url, timeout=timeout)
                    if response.ok:
                        return response.json()["results"]
                    else:
                        self.logger.warning(f"API request failed with status {response.status_code}")
                except Exception as e:
                    self.logger.warning(f"API request error (attempt {attempt}/{tries}): {e}")
                
                if attempt < tries:
                    time.sleep(2**attempt)
            
            raise RuntimeError(f"Failed after {tries} tries: {url}")
        
        go_terms_info = {}
        chunk_size = 200
        
        self.logger.info(f"Fetching GO terms for {len(go_ids)} GO IDs")
        
        try:
            for batch in chunks(go_ids, chunk_size):
                results = fetch(batch)
                for term in results:
                    go_terms_info[term["id"]] = {
                        'name': term["name"],
                        'definition': term["definition"]["text"] if term.get("definition") else ""
                    }
                
                # Add small delay between requests to be respectful to API
                time.sleep(0.5)
            
            self.logger.info(f"Successfully fetched {len(go_terms_info)} GO terms from API")
            
        except Exception as e:
            self.logger.error(f"Failed to fetch GO terms from API: {e}")
            # Return empty dict, calling function will handle this
            
        return go_terms_info
    
    def _parse_go_annotated_genes(self, results):
        """Parse GO-annotated genes from go_reverse_lookup.tsv"""
        
        go_reverse_lookup_path = results['data'].get('go_reverse_lookup')
        if not go_reverse_lookup_path or not os.path.exists(go_reverse_lookup_path):
            self.logger.warning("go_reverse_lookup.tsv not found")
            return {}
        
        go_genes_map = {}
        
        try:
            self.logger.info(f"Parsing GO-annotated genes from: {go_reverse_lookup_path}")
            
            with open(go_reverse_lookup_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f):
                    line = line.strip()
                    if not line:
                        continue
                    
                    # Skip header line if present
                    if line_num == 0 and line.startswith('go_term'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        go_id = parts[0]
                        gene_count = parts[1]
                        genes_string = parts[2]
                        
                        # Parse genes - they appear to be separated by semicolons
                        # and contain gene names with GO annotations
                        gene_names = []
                        if genes_string:
                            # Split by semicolons and extract gene names
                            gene_entries = genes_string.split(';')
                            for entry in gene_entries:
                                if '|' in entry:
                                    # Extract gene name before the first '|'
                                    gene_name = entry.split('|')[0].strip()
                                    if gene_name and gene_name not in gene_names:
                                        gene_names.append(gene_name)
                        
                        go_genes_map[go_id] = gene_names
                    
                    # Process only first 10000 lines to avoid memory issues
                    if line_num > 10000:
                        break
            
            self.logger.info(f"Parsed GO-annotated genes for {len(go_genes_map)} GO terms")
            
        except Exception as e:
            self.logger.error(f"Error parsing GO-annotated genes: {e}")
            return {}
        
        return go_genes_map
    
    def _fetch_go_terms_from_database(self, go_ids):
        """Fetch GO terms and definitions from existing database to avoid duplicate API calls"""
        
        db_path = os.path.join(self.web_data_dir, 'twheel_data.db')
        
        if not os.path.exists(db_path):
            self.logger.warning("Database not found, falling back to API")
            return self._fetch_go_terms_from_api(go_ids)
        
        go_terms_info = {}
        
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create a placeholder string for the IN clause
            placeholders = ','.join(['?' for _ in go_ids])
            
            cursor.execute(f"""
                SELECT go_id, go_name, go_definition 
                FROM go_terms 
                WHERE go_id IN ({placeholders})
            """, go_ids)
            
            results = cursor.fetchall()
            
            for go_id, go_name, go_definition in results:
                go_terms_info[go_id] = {
                    'name': go_name or '',
                    'definition': go_definition or ''
                }
            
            conn.close()
            
            self.logger.info(f"Successfully loaded {len(go_terms_info)} GO terms from database")
            
            # Check for missing GO terms
            missing_go_ids = [go_id for go_id in go_ids if go_id not in go_terms_info]
            if missing_go_ids:
                self.logger.info(f"Found {len(missing_go_ids)} GO terms missing from database, fetching from API")
                try:
                    missing_go_info = self._fetch_go_terms_from_api(missing_go_ids)
                    go_terms_info.update(missing_go_info)
                except Exception as e:
                    self.logger.warning(f"Failed to fetch missing GO terms from API: {e}")
            
        except sqlite3.Error as e:
            self.logger.error(f"Database error when fetching GO terms: {e}")
            # Fall back to API if database fails
            self.logger.info("Falling back to API due to database error")
            go_terms_info = self._fetch_go_terms_from_api(go_ids)
        
        return go_terms_info
    
    def _format_annotated_genes_list(self, gene_list):
        """Format gene list as 'first 5 genes...and X more' format"""
        
        if not gene_list:
            return ""
        
        # Remove duplicates while preserving order
        unique_genes = []
        seen = set()
        for gene in gene_list:
            if gene not in seen:
                unique_genes.append(gene)
                seen.add(gene)
        
        if len(unique_genes) <= 5:
            return ", ".join(unique_genes)
        else:
            first_five = ", ".join(unique_genes[:5])
            remaining_count = len(unique_genes) - 5
            return f"{first_five}...and {remaining_count} more"
    
    def generate_html_report(self, results):
        """Generate HTML report"""
        self.logger.info("Generating HTML report...")
        
        # Create template files
        self._create_html_template()
        self._create_css_files()
        self._create_js_files()
        self._copy_static_assets()
        
        # Generate main HTML file
        html_path = os.path.join(self.output_dir, 'report.html')
        
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(self._generate_html_content(results))
        
        self.logger.info(f"HTML report generated: {html_path}")
        return html_path
    
    def _copy_static_assets(self):
        """Copy static assets like images to output directory"""
        # Create img directory in output folder
        img_dir = os.path.join(self.output_dir, 'static', 'img')
        Path(img_dir).mkdir(parents=True, exist_ok=True)
        
        # Copy logo file if it exists
        src_logo = os.path.join('src', 'web', 'static', 'img', 'logo.png')
        dst_logo = os.path.join(img_dir, 'logo.png')
        
        if os.path.exists(src_logo):
            shutil.copy2(src_logo, dst_logo)
            self.logger.info(f"Logo copied: {src_logo} -> {dst_logo}")
        
        # Copy GO-related images from intermediate directory
        go_images = [
            'go_parent_counts.png',
            'treemap_go_parent_counts.png'
        ]
        
        for image_name in go_images:
            src_image = os.path.join(self.intermediate_dir, image_name)
            dst_image = os.path.join(img_dir, image_name)
            
            if os.path.exists(src_image):
                shutil.copy2(src_image, dst_image)
                self.logger.info(f"GO image copied: {src_image} -> {dst_image}")
            else:
                self.logger.warning(f"GO image not found: {src_image}")
    
    def _create_html_template(self):
        """Create HTML template"""
        template_dir = os.path.join('src', 'web', 'templates')
        Path(template_dir).mkdir(parents=True, exist_ok=True)
        
        # Basic template files will be created later
        pass
    
    def _create_css_files(self):
        """Create CSS files in output directory"""
        # Create CSS directory in output folder
        css_dir = os.path.join(self.output_dir, 'static', 'css')
        Path(css_dir).mkdir(parents=True, exist_ok=True)
        
        # Main CSS file
        main_css_path = os.path.join(css_dir, 'main.css')
        with open(main_css_path, 'w', encoding='utf-8') as f:
            f.write("""
/* tWhEEL Report CSS - MultiQC Style */
* {
    box-sizing: border-box;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    line-height: 1.5;
    margin: 0;
    padding: 0;
    background-color: #f8f9fa;
    color: #212529;
}

/* Layout */
.report-layout {
    display: flex;
    min-height: 100vh;
}

/* Sidebar */
.sidebar {
    width: 280px;
    background-color: #fff;
    border-right: 1px solid #dee2e6;
    padding: 0;
    position: fixed;
    height: 100vh;
    overflow-y: auto;
    z-index: 1000;
}

.sidebar-header {
    padding: 1.5rem;
    border-bottom: 1px solid #dee2e6;
    text-align: center;
}

.sidebar-header img {
    max-width: 150px;
    height: auto;
    margin-bottom: 0.5rem;
}

.sidebar-header h1 {
    font-size: 1.75rem;
    font-weight: 600;
    color: #5cb85c;
    margin: 0;
    text-decoration: none;
}

.sidebar-header p {
    color: #6c757d;
    font-size: 0.875rem;
    margin: 0.25rem 0 0 0;
    font-weight: 400;
}

.sidebar-nav {
    padding: 0;
}

.sidebar-nav ul {
    list-style: none;
    padding: 0;
    margin: 0;
}

.sidebar-nav li {
    margin: 0;
}

.sidebar-nav a {
    display: block;
    padding: 0.75rem 1.5rem;
    color: #495057;
    text-decoration: none;
    font-size: 0.9rem;
    border-left: 3px solid transparent;
    transition: all 0.2s ease;
}

.sidebar-nav a:hover {
    background-color: #f8f9fa;
    color: #007bff;
    border-left-color: #dee2e6;
}

.sidebar-nav a.active {
    background-color: #e7f3ff;
    color: #007bff;
    border-left-color: #007bff;
    font-weight: 500;
}

/* Main content */
.main-content {
    margin-left: 280px;
    flex: 1;
    padding: 0;
}

.content-header {
    background-color: #fff;
    border-bottom: 1px solid #dee2e6;
    padding: 1.5rem 2rem;
}

.content-header h1 {
    font-size: 2rem;
    font-weight: 600;
    color: #212529;
    margin: 0;
}

.content-header .subtitle {
    color: #6c757d;
    font-size: 1rem;
    margin: 0.25rem 0 0 0;
}

.content-body {
    padding: 2rem;
}

/* Sections */
.content-section {
    background-color: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    margin-bottom: 2rem;
    overflow: hidden;
}

.section-header {
    background-color: #f8f9fa;
    border-bottom: 1px solid #dee2e6;
    padding: 1rem 1.5rem;
}

.section-header h2 {
    font-size: 1.25rem;
    font-weight: 600;
    color: #212529;
    margin: 0;
}

.section-body {
    padding: 1.5rem;
}

/* Statistics Grid */
.stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
}

.stat-card {
    background: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    padding: 1.25rem;
    text-align: center;
    transition: box-shadow 0.15s ease-in-out;
}

.stat-card:hover {
    box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
}

.stat-value {
    font-size: 2rem;
    font-weight: 700;
    color: #28a745;
    line-height: 1.2;
}

.stat-label {
    color: #6c757d;
    font-size: 0.875rem;
    margin-top: 0.5rem;
    font-weight: 500;
}

/* Tables */
.table-responsive {
    overflow-x: auto;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
}

.data-table th,
.data-table td {
    padding: 0.75rem;
    text-align: left;
    border-bottom: 1px solid #dee2e6;
    vertical-align: top;
}

.data-table th {
    background-color: #f8f9fa;
    font-weight: 600;
    color: #495057;
    border-top: 1px solid #dee2e6;
    position: sticky;
    top: 0;
}

.data-table tbody tr:hover {
    background-color: #f8f9fa;
}

/* Form controls */
.filter-group {
    display: flex;
    gap: 1rem;
    margin-bottom: 1.5rem;
    flex-wrap: wrap;
    align-items: center;
}

.search-box {
    flex: 1;
    min-width: 250px;
    padding: 0.5rem 0.75rem;
    border: 1px solid #ced4da;
    border-radius: 0.375rem;
    font-size: 0.9rem;
}

.filter-select {
    padding: 0.5rem 0.75rem;
    border: 1px solid #ced4da;
    border-radius: 0.375rem;
    font-size: 0.9rem;
    min-width: 150px;
}

/* Charts */
.chart-container {
    height: 400px;
    margin: 1.5rem 0;
    background-color: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    display: flex;
    align-items: center;
    justify-content: center;
}

/* Special elements */
.guide-sequence {
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    background-color: #f8f9fa;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    border: 1px solid #dee2e6;
}

.badge {
    display: inline-block;
    padding: 0.25em 0.6em;
    font-size: 0.75em;
    font-weight: 700;
    line-height: 1;
    text-align: center;
    white-space: nowrap;
    vertical-align: baseline;
    border-radius: 0.375rem;
}

.badge-success {
    background-color: #28a745;
    color: #fff;
}

.badge-warning {
    background-color: #ffc107;
    color: #212529;
}

.badge-danger {
    background-color: #dc3545;
    color: #fff;
}

.targetable-fully {
    background-color: #d1ecf1;
    color: #0c5460;
    border: 1px solid #bee5eb;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}

.targetable-partially {
    background-color: #fff3cd;
    color: #856404;
    border: 1px solid #ffeaa7;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}

.targetable-not {
    background-color: #f8d7da;
    color: #721c24;
    border: 1px solid #f5c6cb;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}

/* Loading and error states */
.loading {
    text-align: center;
    padding: 3rem;
    color: #6c757d;
    font-size: 0.9rem;
}

.error {
    background-color: #f8d7da;
    border: 1px solid #f5c6cb;
    color: #721c24;
    padding: 1rem;
    border-radius: 0.375rem;
    margin: 1rem 0;
}

.alert {
    padding: 0.75rem 1rem;
    margin-bottom: 1rem;
    border: 1px solid transparent;
    border-radius: 0.375rem;
}

.alert-info {
    color: #0c5460;
    background-color: #d1ecf1;
    border-color: #bee5eb;
}

/* Responsive design */
@media (max-width: 768px) {
    .sidebar {
        transform: translateX(-100%);
        transition: transform 0.3s ease;
    }
    
    .sidebar.show {
        transform: translateX(0);
    }
    
    .main-content {
        margin-left: 0;
    }
    
    .filter-group {
        flex-direction: column;
        align-items: stretch;
    }
    
    .search-box {
        min-width: auto;
    }
    
    .stats-grid {
        grid-template-columns: 1fr;
    }
}

/* Smooth scrolling */
html {
    scroll-behavior: smooth;
}

/* Protein Sequence Viewer Styles */
.protein-item {
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    margin-bottom: 0.5rem;
    background-color: #fff;
}

.protein-header {
    padding: 0.75rem 1rem;
    background-color: #f8f9fa;
    border-bottom: 1px solid #dee2e6;
    cursor: pointer;
    display: flex;
    justify-content: space-between;
    align-items: center;
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.9rem;
    transition: background-color 0.2s ease;
}

.protein-header:hover {
    background-color: #e9ecef;
}

.protein-header.expanded {
    border-bottom-color: transparent;
}

.protein-id {
    font-weight: 600;
    color: #495057;
    flex: 1;
}

.protein-stats {
    display: flex;
    gap: 1rem;
    font-size: 0.8rem;
    color: #6c757d;
}

.expand-toggle {
    font-size: 0.7rem;
    color: #6c757d;
    font-weight: bold;
    min-width: 20px;
    text-align: center;
}

.protein-content {
    display: none;
    padding: 1rem;
    background-color: #fff;
    border-top: 1px solid #dee2e6;
}

.protein-content.expanded {
    display: block;
}

.sequence-container {
    background-color: #f8f9fa;
    border: 1px solid #dee2e6;
    border-radius: 0.25rem;
    padding: 1rem;
    margin-bottom: 1rem;
}

.protein-sequence {
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.85rem;
    line-height: 1.6;
    word-break: break-all;
    word-wrap: break-word;
    color: #495057;
    white-space: pre-wrap;
}

.sequence-controls {
    margin-bottom: 1rem;
    display: flex;
    gap: 0.5rem;
    flex-wrap: wrap;
}

.btn {
    padding: 0.375rem 0.75rem;
    border: 1px solid #ced4da;
    border-radius: 0.375rem;
    background-color: #fff;
    color: #495057;
    cursor: pointer;
    font-size: 0.875rem;
    text-decoration: none;
    display: inline-block;
    transition: all 0.2s ease;
}

.btn:hover {
    background-color: #e9ecef;
    border-color: #adb5bd;
}

.btn-primary {
    background-color: #007bff;
    border-color: #007bff;
    color: #fff;
}

.btn-primary:hover {
    background-color: #0056b3;
    border-color: #004085;
}

.btn-secondary {
    background-color: #6c757d;
    border-color: #6c757d;
    color: #fff;
}

.btn-secondary:hover {
    background-color: #545b62;
    border-color: #4e555b;
}

.btn-outline {
    background-color: transparent;
    border-color: #ced4da;
    color: #495057;
}

.btn-outline:hover {
    background-color: #f8f9fa;
    border-color: #adb5bd;
}

.sequence-info {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    margin-bottom: 1rem;
    font-size: 0.875rem;
}

.info-item {
    display: flex;
    justify-content: space-between;
    padding: 0.5rem;
    background-color: #f8f9fa;
    border-radius: 0.25rem;
}

.info-label {
    font-weight: 600;
    color: #495057;
}

.info-value {
    color: #6c757d;
}

.loading-spinner {
    text-align: center;
    padding: 2rem;
    color: #6c757d;
}

.pagination-controls {
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 1rem;
    margin-top: 2rem;
    padding: 1rem;
    border-top: 1px solid #dee2e6;
}

.no-results {
    text-align: center;
    padding: 2rem;
    color: #6c757d;
    font-style: italic;
}

/* Length-based coloring */
.protein-header[data-length="short"] .protein-id::before {
    content: "●";
    color: #28a745;
    margin-right: 0.5rem;
}

.protein-header[data-length="medium"] .protein-id::before {
    content: "●";
    color: #ffc107;
    margin-right: 0.5rem;
}

.protein-header[data-length="long"] .protein-id::before {
    content: "●";
    color: #dc3545;
    margin-right: 0.5rem;
}

/* Search highlighting */
.highlight {
    background-color: #fff3cd;
    padding: 0.125rem 0.25rem;
    border-radius: 0.25rem;
    font-weight: 600;
}

/* Animation for expand/collapse */
.protein-content {
    transition: all 0.3s ease;
    overflow: hidden;
}

/* Expandable genes styles */
.expandable-genes {
    cursor: pointer;
    position: relative;
    transition: background-color 0.2s ease;
}

.expandable-genes:hover {
    background-color: #f8f9fa;
}

.expandable-genes.expanded {
    background-color: #e7f3ff;
}

.expand-icon {
    margin-left: 0.5rem;
    font-size: 0.75rem;
    opacity: 0.7;
    transition: opacity 0.2s ease;
}

.expandable-genes:hover .expand-icon {
    opacity: 1.0;
}

.genes-content {
    display: inline;
}

/* GO Terms Stat Card Colors */
.stat-value-fully-targetable { 
    background-color: #1B420A; 
    color: white; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}
.stat-value-simply-targetable { 
    background-color: #59BF10; 
    color: white; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}
.stat-value-multiply-targetable { 
    background-color: #D7DB01; 
    color: #333; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}
.stat-value-partially-targetable { 
    background-color: #F0E579; 
    color: #333; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}
.stat-value-not-targetable { 
    background-color: #D92000; 
    color: white; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}
.stat-value-total-go { 
    background-color: #777777; 
    color: white; 
    padding: 0.5rem;
    border-radius: 0.375rem;
    font-weight: bold;
}

/* Print styles */
@media print {
    .sidebar {
        display: none;
    }
    
    .main-content {
        margin-left: 0;
    }
    
    .protein-content {
        display: block !important;
    }
    
    .sequence-controls {
        display: none;
    }
    
    .expandable-genes {
        cursor: default;
    }
    
    .expand-icon {
        display: none;
    }
}
""")
    
    def _create_js_files(self):
        """Create JavaScript files in output directory"""
        # Create JS directory in output folder
        js_dir = os.path.join(self.output_dir, 'static', 'js')
        Path(js_dir).mkdir(parents=True, exist_ok=True)
        
        # Main JS file
        main_js_path = os.path.join(js_dir, 'main.js')
        with open(main_js_path, 'w', encoding='utf-8') as f:
            f.write("""
// tWhEEL Report JavaScript - MultiQC Style

class tWhEELReport {
    constructor() {
        this.currentSection = 'general-stats';
        this.data = null;
        this.proteins = [];
        this.filteredProteins = [];
        this.currentPage = 1;
        this.pageSize = 10;
        
        // セクションの読み込み状態を管理
        this.sectionsLoaded = {
            geneAnnotation: false,
            guideRna: false,
            goAnalysis: false,
            geneSequences: false
        };
        
        this.init();
    }
    
    async init() {
        await this.loadData();
        this.setupEventListeners();
        this.initializeNavigation();
        this.scrollToActiveSection();
    }
    
    async loadData() {
        try {
            // Load embedded data instead of fetching from external file
            const dataElement = document.getElementById('twheel-data');
            if (dataElement) {
                this.data = JSON.parse(dataElement.textContent);
                this.populateGeneralStats();
            } else {
                throw new Error('Embedded data not found');
            }
        } catch (error) {
            console.error('Failed to load embedded data:', error);
            this.showError('Failed to load data. Please check the report file.');
        }
    }
    
    setupEventListeners() {
        // Sidebar navigation
        document.querySelectorAll('.nav-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = e.target.getAttribute('href').substring(1);
                this.scrollToSection(targetId);
            });
        });
        
        // Search functionality
        document.addEventListener('input', (e) => {
            if (e.target.classList.contains('search-box')) {
                this.handleSearch(e.target.value, e.target.dataset.target);
            }
        });
        
        // Filter functionality
        document.addEventListener('change', (e) => {
            if (e.target.classList.contains('filter-select')) {
                this.handleFilter(e.target.value, e.target.dataset.filter);
            }
        });
        
        // Scroll spy functionality
        window.addEventListener('scroll', () => {
            this.updateActiveNavigation();
        });
        
        // Protein sequence viewer event listeners
        this.setupProteinViewerListeners();
        
        // Mobile sidebar toggle (future implementation)
        const toggleBtn = document.getElementById('sidebar-toggle');
        if (toggleBtn) {
            toggleBtn.addEventListener('click', () => {
                document.querySelector('.sidebar').classList.toggle('show');
            });
        }
    }
    
    setupProteinViewerListeners() {
        // Protein search
        const proteinSearch = document.getElementById('protein-search');
        if (proteinSearch) {
            proteinSearch.addEventListener('input', debounce((e) => {
                this.searchProteins(e.target.value);
            }, 300));
        }
        
        // Length filter
        const lengthFilter = document.getElementById('length-filter');
        if (lengthFilter) {
            lengthFilter.addEventListener('change', (e) => {
                this.filterProteinsByLength(e.target.value);
            });
        }
        
        // Expand/Collapse all buttons
        const expandAllBtn = document.getElementById('expand-all');
        const collapseAllBtn = document.getElementById('collapse-all');
        
        if (expandAllBtn) {
            expandAllBtn.addEventListener('click', () => this.expandAllProteins());
        }
        
        if (collapseAllBtn) {
            collapseAllBtn.addEventListener('click', () => this.collapseAllProteins());
        }
        
        // Pagination controls
        const prevPageBtn = document.getElementById('prev-page');
        const nextPageBtn = document.getElementById('next-page');
        const pageSizeSelect = document.getElementById('page-size');
        
        if (prevPageBtn) {
            prevPageBtn.addEventListener('click', () => this.goToPreviousPage());
        }
        
        if (nextPageBtn) {
            nextPageBtn.addEventListener('click', () => this.goToNextPage());
        }
        
        if (pageSizeSelect) {
            pageSizeSelect.addEventListener('change', (e) => {
                this.pageSize = parseInt(e.target.value);
                this.currentPage = 1;
                this.renderProteins();
            });
        }
    }
    
    initializeNavigation() {
        // Set up smooth scrolling and initial active states
        this.updateActiveNavigation();
    }
    
    scrollToSection(sectionId) {
        const element = document.getElementById(sectionId);
        if (element) {
            const offsetTop = element.offsetTop - 20; // Account for fixed elements
            window.scrollTo({
                top: offsetTop,
                behavior: 'smooth'
            });
        }
    }
    
    scrollToActiveSection() {
        // Scroll to the first section on load
        if (this.currentSection) {
            this.scrollToSection(this.currentSection);
        }
    }
    
    updateActiveNavigation() {
        const sections = document.querySelectorAll('.content-section');
        const navLinks = document.querySelectorAll('.nav-link');
        
        let currentActive = null;
        const scrollPos = window.scrollY + 100; // Offset for better UX
        
        sections.forEach(section => {
            const sectionTop = section.offsetTop;
            const sectionHeight = section.offsetHeight;
            
            if (scrollPos >= sectionTop && scrollPos < sectionTop + sectionHeight) {
                currentActive = section.id;
            }
        });
        
        // Update navigation active states
        navLinks.forEach(link => {
            link.classList.remove('active');
            if (currentActive && link.getAttribute('href') === `#${currentActive}`) {
                link.classList.add('active');
                
                // Load section-specific data when section becomes active
                if (this.currentSection !== currentActive) {
                    this.currentSection = currentActive;
                    this.onSectionChange(currentActive);
                }
            }
        });
    }
    
    onSectionChange(sectionId) {
        // Handle section-specific loading
        switch(sectionId) {
            case 'gene-prediction':
                if (!this.sectionsLoaded.geneSequences && this.proteins.length === 0) {
                    this.loadProteinSequences();
                }
                break;
            case 'gene-annotation':
                if (!this.sectionsLoaded.geneAnnotation) {
                    this.loadGeneData();
                }
                break;
            case 'guide-rna':
                if (!this.sectionsLoaded.guideRna) {
                    this.loadGuideData();
                }
                break;
            case 'go-analysis':
                if (!this.sectionsLoaded.goAnalysis) {
                    this.loadGOData();
                    this.loadGOCharts();
                }
                break;
        }
    }
    
    populateGeneralStats() {
        if (!this.data) return;
        
        const stats = this.data.statistics;
        
        // Update statistics cards
        this.updateStatCard('protein-count', stats.gene_prediction?.protein_count || 0);
        this.updateStatCard('total-genes', stats.annotation?.total_genes || 0);
        this.updateStatCard('multi-genes', stats.annotation?.multi_genes || 0);
        this.updateStatCard('single-genes', stats.annotation?.single_genes || 0);
        this.updateStatCard('total-gos', stats.go_analysis?.total_gos || 0);
        this.updateStatCard('fully-targetable', (stats.go_analysis?.fully_targetable_percent || 0).toFixed(1) + '%');
    }
    
    updateStatCard(id, value) {
        const element = document.getElementById(id);
        if (element) {
            if (typeof value === 'number') {
                element.textContent = value.toLocaleString();
            } else {
                element.textContent = value;
            }
        }
    }
    
    showError(message) {
        const alertDiv = document.createElement('div');
        alertDiv.className = 'alert alert-danger';
        alertDiv.innerHTML = `<strong>Error:</strong> ${message}`;
        
        const contentBody = document.querySelector('.content-body');
        if (contentBody) {
            contentBody.insertBefore(alertDiv, contentBody.firstChild);
        }
    }
    
    showSuccess(message) {
        const alertDiv = document.createElement('div');
        alertDiv.className = 'alert alert-success';
        alertDiv.innerHTML = `<strong>Success:</strong> ${message}`;
        
        const contentBody = document.querySelector('.content-body');
        if (contentBody) {
            contentBody.insertBefore(alertDiv, contentBody.firstChild);
        }
        
        // Auto-remove after 3 seconds
        setTimeout(() => {
            alertDiv.remove();
        }, 3000);
    }
    
    async loadGeneData() {
        console.log('Loading gene data...');
        try {
            const geneTable = document.querySelector('#genes-table tbody');
            if (!geneTable) return;
            
            // Load embedded gene annotation data
            if (this.data && this.data.gene_annotation && Array.isArray(this.data.gene_annotation)) {
                const geneData = this.data.gene_annotation;
                
                if (geneData.length === 0) {
                    geneTable.innerHTML = '<tr><td colspan="6" class="loading">No gene annotation data available</td></tr>';
                    return;
                }
                
                // Generate table rows
                const rows = geneData.map(gene => {
                    const classificationClass = gene.classification === 'Multi' ? 'classification-multi' : 'classification-single';
                    const guideCountClass = gene.target_count > 0 ? 'guide-count' : 'guide-count-zero';
                    
                    // Truncate long descriptions for display
                    const description = gene.description.length > 100 
                        ? gene.description.substring(0, 100) + '...' 
                        : gene.description;
                    
                    return `
                        <tr>
                            <td><span class="query-id">${gene.query_id}</span></td>
                            <td><span class="uniprot-id">${gene.subject_id}</span></td>
                            <td><strong>${gene.gene_name}</strong></td>
                            <td class="description-cell" title="${gene.description}">${description}</td>
                            <td>
                                <span class="${classificationClass}">${gene.classification}</span><br>
                                <small>${gene.unique_id}</small>
                            </td>
                            <td><span class="${guideCountClass}">${gene.target_count}</span></td>
                        </tr>
                    `;
                }).join('');
                
                geneTable.innerHTML = rows;
                console.log(`Loaded ${geneData.length} gene annotation entries`);
                
                // Setup search and filter functionality for gene data
                this.setupGeneTableFilters(geneData);
                
                // Mark section as loaded
                this.sectionsLoaded.geneAnnotation = true;
                
            } else {
                geneTable.innerHTML = '<tr><td colspan="6" class="error">Gene annotation data not found in embedded data</td></tr>';
            }
            
        } catch (error) {
            console.error('Failed to load gene data:', error);
            const geneTable = document.querySelector('#genes-table tbody');
            if (geneTable) {
                geneTable.innerHTML = '<tr><td colspan="6" class="error">Error loading gene data</td></tr>';
            }
        }
    }
    
    setupGeneTableFilters(geneData) {
        // Store original data for filtering and pagination
        this.originalGeneData = [...geneData];
        this.filteredGeneData = [...geneData];
        this.currentGenePage = 1;
        this.genePageSize = 10;
        
        // Search functionality
        const searchBox = document.querySelector('#gene-annotation .search-box');
        if (searchBox) {
            searchBox.addEventListener('input', debounce((e) => {
                this.currentGenePage = 1; // Reset to first page on search
                this.filterGeneData();
            }, 300));
        }
        
        // Classification filter
        const classificationFilter = document.querySelector('#gene-annotation .filter-select[data-filter="classification"]');
        if (classificationFilter) {
            classificationFilter.addEventListener('change', () => {
                this.currentGenePage = 1; // Reset to first page on filter change
                this.filterGeneData();
            });
        }
        
        // Guide RNA filter
        const guidesFilter = document.querySelector('#gene-annotation .filter-select[data-filter="guides"]');
        if (guidesFilter) {
            guidesFilter.addEventListener('change', () => {
                this.currentGenePage = 1; // Reset to first page on filter change
                this.filterGeneData();
            });
        }
        
        // Pagination controls
        const prevPageBtn = document.getElementById('genes-prev-page');
        const nextPageBtn = document.getElementById('genes-next-page');
        const pageSizeSelect = document.getElementById('genes-page-size');
        
        if (prevPageBtn) {
            prevPageBtn.addEventListener('click', () => this.goToGenesPreviousPage());
        }
        
        if (nextPageBtn) {
            nextPageBtn.addEventListener('click', () => this.goToGenesNextPage());
        }
        
        if (pageSizeSelect) {
            pageSizeSelect.addEventListener('change', (e) => {
                this.genePageSize = parseInt(e.target.value);
                this.currentGenePage = 1;
                this.renderFilteredGeneData();
            });
        }
        
        // Initial render with pagination
        this.renderFilteredGeneData();
    }
    
    filterGeneData() {
        const searchQuery = document.querySelector('#gene-annotation .search-box').value.toLowerCase();
        const classificationFilter = document.querySelector('#gene-annotation .filter-select[data-filter="classification"]').value;
        const guidesFilter = document.querySelector('#gene-annotation .filter-select[data-filter="guides"]').value;
        
            this.filteredGeneData = this.originalGeneData.filter(gene => {
                // Search filter
                const matchesSearch = !searchQuery || 
                    gene.query_id.toLowerCase().includes(searchQuery) ||
                    gene.gene_name.toLowerCase().includes(searchQuery) ||
                    gene.description.toLowerCase().includes(searchQuery) ||
                    gene.classification.toLowerCase().includes(searchQuery) ||
                    gene.unique_id.toLowerCase().includes(searchQuery);
            
            // Classification filter
            const matchesClassification = !classificationFilter || gene.classification === classificationFilter;
            
            // Guides filter
            let matchesGuides = true;
            if (guidesFilter === 'with-guides') {
                matchesGuides = gene.target_count > 0;
            } else if (guidesFilter === 'without-guides') {
                matchesGuides = gene.target_count === 0;
            }
            
            return matchesSearch && matchesClassification && matchesGuides;
        });
        
        // Update table with filtered data
        this.renderFilteredGeneData();
    }
    
    renderFilteredGeneData() {
        const geneTable = document.querySelector('#genes-table tbody');
        if (!geneTable) return;
        
        if (this.filteredGeneData.length === 0) {
            geneTable.innerHTML = '<tr><td colspan="6" class="no-results">No genes found matching the current filters.</td></tr>';
            this.updateGenePaginationControls();
            return;
        }
        
        // Calculate pagination for gene data
        const totalPages = Math.ceil(this.filteredGeneData.length / this.genePageSize);
        const startIndex = (this.currentGenePage - 1) * this.genePageSize;
        const endIndex = Math.min(startIndex + this.genePageSize, this.filteredGeneData.length);
        const currentGeneData = this.filteredGeneData.slice(startIndex, endIndex);
        
        const rows = currentGeneData.map(gene => {
            const classificationClass = gene.classification === 'Multi' ? 'classification-multi' : 'classification-single';
            const guideCountClass = gene.target_count > 0 ? 'guide-count' : 'guide-count-zero';
            
            // Truncate long descriptions for display
            const description = gene.description.length > 100 
                ? gene.description.substring(0, 100) + '...' 
                : gene.description;
            
            return `
                <tr>
                    <td><span class="query-id">${gene.query_id}</span></td>
                    <td><span class="uniprot-id">${gene.subject_id}</span></td>
                    <td><strong>${gene.gene_name}</strong></td>
                    <td class="description-cell" title="${gene.description}">${description}</td>
                    <td>
                        <span class="${classificationClass}">${gene.classification}</span><br>
                        <small>${gene.unique_id}</small>
                    </td>
                    <td><span class="${guideCountClass}">${gene.target_count}</span></td>
                </tr>
            `;
        }).join('');
        
        geneTable.innerHTML = rows;
        
        // Update pagination controls
        this.updateGenePaginationControls();
    }
    
    // Gene Annotation pagination methods
    goToGenesPreviousPage() {
        if (this.currentGenePage > 1) {
            this.currentGenePage--;
            this.renderFilteredGeneData();
        }
    }
    
    goToGenesNextPage() {
        const totalPages = Math.ceil(this.filteredGeneData.length / this.genePageSize);
        if (this.currentGenePage < totalPages) {
            this.currentGenePage++;
            this.renderFilteredGeneData();
        }
    }
    
    updateGenePaginationControls() {
        const paginationControls = document.getElementById('genes-pagination-controls');
        const prevBtn = document.getElementById('genes-prev-page');
        const nextBtn = document.getElementById('genes-next-page');
        const pageInfo = document.getElementById('genes-page-info');
        const showingInfo = document.getElementById('genes-showing-info');
        
        if (!paginationControls) return;
        
        const totalPages = Math.ceil(this.filteredGeneData.length / this.genePageSize);
        const startIndex = (this.currentGenePage - 1) * this.genePageSize;
        const endIndex = Math.min(startIndex + this.genePageSize, this.filteredGeneData.length);
        
        // Show/hide pagination controls
        paginationControls.style.display = totalPages > 1 ? 'flex' : 'none';
        
        // Update button states
        if (prevBtn) prevBtn.disabled = this.currentGenePage <= 1;
        if (nextBtn) nextBtn.disabled = this.currentGenePage >= totalPages;
        
        // Update page info
        if (pageInfo) pageInfo.textContent = `Page ${this.currentGenePage} of ${totalPages}`;
        
        // Update showing info
        if (showingInfo && this.filteredGeneData.length > 0) {
            showingInfo.textContent = `Showing ${startIndex + 1}-${endIndex} of ${this.filteredGeneData.length} entries`;
        } else if (showingInfo) {
            showingInfo.textContent = 'Showing 0-0 of 0 entries';
        }
    }
    
    async loadGuideData() {
        console.log('Loading guide data...');
        try {
            const guideTable = document.querySelector('#guides-table tbody');
            if (!guideTable) return;
            
            // Load embedded guide RNA data
            if (this.data && this.data.guide_rna && Array.isArray(this.data.guide_rna)) {
                const guideData = this.data.guide_rna;
                
                if (guideData.length === 0) {
                    guideTable.innerHTML = '<tr><td colspan="9" class="loading">No guide RNA data available</td></tr>';
                    return;
                }
                
                // Generate table rows
                const rows = guideData.map(guide => {
                    const strandClass = guide.strand === '+' ? 'strand-plus' : 'strand-minus';
                    const mm0Class = guide.mm0 === 0 ? 'mm-perfect' : 'mm-mismatch';
                    const mm12Class = guide.mm1_2 === 0 ? 'mm-perfect' : 'mm-mismatch';
                    
                    return `
                        <tr>
                            <td><span class="query-id">${guide.query_id}</span></td>
                            <td><span class="reference-id">${guide.reference}</span></td>
                            <td>${guide.start}</td>
                            <td>${guide.end}</td>
                            <td><span class="guide-sequence">${guide.target_sequence}</span></td>
                            <td><span class="${strandClass}">${guide.strand}</span></td>
                            <td><span class="${mm0Class}">${guide.mm0}</span></td>
                            <td><span class="${mm12Class}">${guide.mm1_2}</span></td>
                            <td>${guide.off_target_score.toFixed(6)}</td>
                        </tr>
                    `;
                }).join('');
                
                guideTable.innerHTML = rows;
                console.log(`Loaded ${guideData.length} guide RNA entries`);
                
                // Setup search and filter functionality for guide data
                this.setupGuideTableFilters(guideData);
                
                // Mark section as loaded
                this.sectionsLoaded.guideRna = true;
                
            } else {
                guideTable.innerHTML = '<tr><td colspan="9" class="error">Guide RNA data not found in embedded data</td></tr>';
            }
            
        } catch (error) {
            console.error('Failed to load guide data:', error);
            const guideTable = document.querySelector('#guides-table tbody');
            if (guideTable) {
                guideTable.innerHTML = '<tr><td colspan="9" class="error">Error loading guide data</td></tr>';
            }
        }
    }
    
    setupGuideTableFilters(guideData) {
        // Store original data for filtering and pagination
        this.originalGuideData = [...guideData];
        this.filteredGuideData = [...guideData];
        this.currentGuidePage = 1;
        this.guidePageSize = 10;
        
        // Search functionality
        const searchBox = document.querySelector('#guide-rna .search-box');
        if (searchBox) {
            searchBox.addEventListener('input', debounce((e) => {
                this.currentGuidePage = 1;
                this.filterGuideData();
            }, 300));
        }
        
        // MM0 filter
        const mm0Filter = document.querySelector('#guide-rna .filter-select[data-filter="mm0"]');
        if (mm0Filter) {
            mm0Filter.addEventListener('change', () => {
                this.currentGuidePage = 1;
                this.filterGuideData();
            });
        }
        
        // MM1-2 filter
        const mm12Filter = document.querySelector('#guide-rna .filter-select[data-filter="mm1-2"]');
        if (mm12Filter) {
            mm12Filter.addEventListener('change', () => {
                this.currentGuidePage = 1;
                this.filterGuideData();
            });
        }
        
        // Pagination controls
        const prevPageBtn = document.getElementById('guides-prev-page');
        const nextPageBtn = document.getElementById('guides-next-page');
        const pageSizeSelect = document.getElementById('guides-page-size');
        
        if (prevPageBtn) {
            prevPageBtn.addEventListener('click', () => this.goToGuidesPreviousPage());
        }
        
        if (nextPageBtn) {
            nextPageBtn.addEventListener('click', () => this.goToGuidesNextPage());
        }
        
        if (pageSizeSelect) {
            pageSizeSelect.addEventListener('change', (e) => {
                this.guidePageSize = parseInt(e.target.value);
                this.currentGuidePage = 1;
                this.renderFilteredGuideData();
            });
        }
        
        // Initial render with pagination
        this.renderFilteredGuideData();
    }
    
    filterGuideData() {
        const searchQuery = document.querySelector('#guide-rna .search-box').value.toLowerCase();
        const mm0Filter = document.querySelector('#guide-rna .filter-select[data-filter="mm0"]').value;
        const mm12Filter = document.querySelector('#guide-rna .filter-select[data-filter="mm1-2"]').value;
        
        this.filteredGuideData = this.originalGuideData.filter(guide => {
            // Search filter
            const matchesSearch = !searchQuery || 
                guide.query_id.toLowerCase().includes(searchQuery) ||
                guide.reference.toLowerCase().includes(searchQuery) ||
                guide.target_sequence.toLowerCase().includes(searchQuery);
            
            // MM0 filter
            let matchesMM0 = true;
            if (mm0Filter) {
                if (mm0Filter === '2') {
                    matchesMM0 = guide.mm0 >= 2;
                } else {
                    matchesMM0 = guide.mm0 === parseInt(mm0Filter);
                }
            }
            
            // MM1-2 filter
            let matchesMM12 = true;
            if (mm12Filter) {
                if (mm12Filter === '2') {
                    matchesMM12 = guide.mm1_2 >= 2;
                } else {
                    matchesMM12 = guide.mm1_2 === parseInt(mm12Filter);
                }
            }
            
            return matchesSearch && matchesMM0 && matchesMM12;
        });
        
        // Update table with filtered data
        this.renderFilteredGuideData();
    }
    
    renderFilteredGuideData() {
        const guideTable = document.querySelector('#guides-table tbody');
        if (!guideTable) return;
        
        if (this.filteredGuideData.length === 0) {
            guideTable.innerHTML = '<tr><td colspan="9" class="no-results">No guide RNAs found matching the current filters.</td></tr>';
            this.updateGuidePaginationControls();
            return;
        }
        
        // Calculate pagination for guide data
        const totalPages = Math.ceil(this.filteredGuideData.length / this.guidePageSize);
        const startIndex = (this.currentGuidePage - 1) * this.guidePageSize;
        const endIndex = Math.min(startIndex + this.guidePageSize, this.filteredGuideData.length);
        const currentGuideData = this.filteredGuideData.slice(startIndex, endIndex);
        
        const rows = currentGuideData.map(guide => {
            const strandClass = guide.strand === '+' ? 'strand-plus' : 'strand-minus';
            const mm0Class = guide.mm0 === 0 ? 'mm-perfect' : 'mm-mismatch';
            const mm12Class = guide.mm1_2 === 0 ? 'mm-perfect' : 'mm-mismatch';
            
            return `
                <tr>
                    <td><span class="query-id">${guide.query_id}</span></td>
                    <td><span class="reference-id">${guide.reference}</span></td>
                    <td>${guide.start}</td>
                    <td>${guide.end}</td>
                    <td><span class="guide-sequence">${guide.target_sequence}</span></td>
                    <td><span class="${strandClass}">${guide.strand}</span></td>
                    <td><span class="${mm0Class}">${guide.mm0}</span></td>
                    <td><span class="${mm12Class}">${guide.mm1_2}</span></td>
                    <td>${guide.off_target_score.toFixed(6)}</td>
                </tr>
            `;
        }).join('');
        
        guideTable.innerHTML = rows;
        
        // Update pagination controls
        this.updateGuidePaginationControls();
    }
    
    // Guide RNA pagination methods
    goToGuidesPreviousPage() {
        if (this.currentGuidePage > 1) {
            this.currentGuidePage--;
            this.renderFilteredGuideData();
        }
    }
    
    goToGuidesNextPage() {
        const totalPages = Math.ceil(this.filteredGuideData.length / this.guidePageSize);
        if (this.currentGuidePage < totalPages) {
            this.currentGuidePage++;
            this.renderFilteredGuideData();
        }
    }
    
    updateGuidePaginationControls() {
        const paginationControls = document.getElementById('guides-pagination-controls');
        const prevBtn = document.getElementById('guides-prev-page');
        const nextBtn = document.getElementById('guides-next-page');
        const pageInfo = document.getElementById('guides-page-info');
        const showingInfo = document.getElementById('guides-showing-info');
        
        if (!paginationControls) return;
        
        const totalPages = Math.ceil(this.filteredGuideData.length / this.guidePageSize);
        const startIndex = (this.currentGuidePage - 1) * this.guidePageSize;
        const endIndex = Math.min(startIndex + this.guidePageSize, this.filteredGuideData.length);
        
        // Show/hide pagination controls
        paginationControls.style.display = totalPages > 1 ? 'flex' : 'none';
        
        // Update button states
        if (prevBtn) prevBtn.disabled = this.currentGuidePage <= 1;
        if (nextBtn) nextBtn.disabled = this.currentGuidePage >= totalPages;
        
        // Update page info
        if (pageInfo) pageInfo.textContent = `Page ${this.currentGuidePage} of ${totalPages}`;
        
        // Update showing info
        if (showingInfo && this.filteredGuideData.length > 0) {
            showingInfo.textContent = `Showing ${startIndex + 1}-${endIndex} of ${this.filteredGuideData.length} entries`;
        } else if (showingInfo) {
            showingInfo.textContent = 'Showing 0-0 of 0 entries';
        }
    }
    
    async loadGOData() {
        console.log('Loading GO data...');
        try {
            const goTable = document.querySelector('#go-table tbody');
            if (!goTable) return;
            
            // Load embedded GO analysis data
            if (this.data && this.data.go_analysis && Array.isArray(this.data.go_analysis)) {
                const goData = this.data.go_analysis;
                
                if (goData.length === 0) {
                    goTable.innerHTML = '<tr><td colspan="10" class="loading">No GO analysis data available</td></tr>';
                    return;
                }
                
                // Generate table rows
                const rows = goData.map(go => {
                    const targetableClass = this.getTargetableClass(go.targetable_type);
                    const quickgoLink = `https://www.ebi.ac.uk/QuickGO/term/${go.go_id}`;
                    
                    return `
                        <tr>
                            <td><span class="go-id">${go.go_id}</span></td>
                            <td><strong>${go.go_name || 'Loading...'}</strong></td>
                            <td>${go.go_definition || 'Loading...'}</td>
                            <td><a href="${quickgoLink}" target="_blank" class="btn btn-outline">QuickGO</a></td>
                            <td>${go.total_genes}</td>
                            <td>${go.genes_with_guides}</td>
                            <td>${go.kfi_score.toFixed(4)}</td>
                            <td><span class="${targetableClass}">${go.targetable_type}</span></td>
                            <td>${go.annotated_genes || 'Loading...'}</td>
                        </tr>
                    `;
                }).join('');
                
                goTable.innerHTML = rows;
                console.log(`Loaded ${goData.length} GO analysis entries`);
                
                // Setup search and filter functionality for GO data
                this.setupGOTableFilters(goData);
                
                // Mark section as loaded
                this.sectionsLoaded.goAnalysis = true;
                
            } else {
                goTable.innerHTML = '<tr><td colspan="9" class="error">GO analysis data not found in embedded data</td></tr>';
            }
            
        } catch (error) {
            console.error('Failed to load GO data:', error);
            const goTable = document.querySelector('#go-table tbody');
            if (goTable) {
                goTable.innerHTML = '<tr><td colspan="10" class="error">Error loading GO data</td></tr>';
            }
        }
    }
    
    getTargetableClass(targetableType) {
        if (targetableType === 'Fully') {
            return 'targetable-fully';
        } else if (targetableType === 'Partially') {
            return 'targetable-partially';
        } else {
            return 'targetable-not';
        }
    }
    
    setupGOTableFilters(goData) {
        // Store original data for filtering and pagination
        this.originalGOData = [...goData];
        this.filteredGOData = [...goData];
        this.currentGOPage = 1;
        this.goPageSize = 25;
        
        // Search functionality
        const searchBox = document.querySelector('#go-analysis .search-box');
        if (searchBox) {
            searchBox.addEventListener('input', debounce((e) => {
                this.currentGOPage = 1;
                this.filterGOData();
            }, 300));
        }
        
        // Targetable filter
        const targetableFilter = document.querySelector('#go-analysis .filter-select[data-filter="targetable"]');
        if (targetableFilter) {
            targetableFilter.addEventListener('change', () => {
                this.currentGOPage = 1;
                this.filterGOData();
            });
        }
        
        // Pagination controls
        const prevPageBtn = document.getElementById('go-prev-page');
        const nextPageBtn = document.getElementById('go-next-page');
        const pageSizeSelect = document.getElementById('go-page-size');
        
        if (prevPageBtn) {
            prevPageBtn.addEventListener('click', () => this.goToGOPreviousPage());
        }
        
        if (nextPageBtn) {
            nextPageBtn.addEventListener('click', () => this.goToGONextPage());
        }
        
        if (pageSizeSelect) {
            pageSizeSelect.addEventListener('change', (e) => {
                this.goPageSize = parseInt(e.target.value);
                this.currentGOPage = 1;
                this.renderFilteredGOData();
            });
        }
        
        // Initial render with pagination
        this.renderFilteredGOData();
    }
    
    filterGOData() {
        const searchQuery = document.querySelector('#go-analysis .search-box').value.toLowerCase();
        const targetableFilter = document.querySelector('#go-analysis .filter-select[data-filter="targetable"]').value;
        
        this.filteredGOData = this.originalGOData.filter(go => {
            // Search filter
            const matchesSearch = !searchQuery || 
                go.go_id.toLowerCase().includes(searchQuery) ||
                go.go_name.toLowerCase().includes(searchQuery) ||
                go.go_definition.toLowerCase().includes(searchQuery) ||
                go.annotated_genes.toLowerCase().includes(searchQuery);
            
            // Targetable filter
            const matchesTargetable = !targetableFilter || go.targetable_type === targetableFilter;
            
            return matchesSearch && matchesTargetable;
        });
        
        // Update table with filtered data
        this.renderFilteredGOData();
    }
    
    renderFilteredGOData() {
        const goTable = document.querySelector('#go-table tbody');
        if (!goTable) return;
        
        if (this.filteredGOData.length === 0) {
            goTable.innerHTML = '<tr><td colspan="10" class="no-results">No GO terms found matching the current filters.</td></tr>';
            this.updateGOPaginationControls();
            return;
        }
        
        // Calculate pagination for GO data
        const totalPages = Math.ceil(this.filteredGOData.length / this.goPageSize);
        const startIndex = (this.currentGOPage - 1) * this.goPageSize;
        const endIndex = Math.min(startIndex + this.goPageSize, this.filteredGOData.length);
        const currentGOData = this.filteredGOData.slice(startIndex, endIndex);
        
        const rows = currentGOData.map(go => {
            const targetableClass = this.getTargetableClass(go.targetable_type);
            const quickgoLink = `https://www.ebi.ac.uk/QuickGO/term/${go.go_id}`;
            
            // Create expandable genes cell
            const hasMoreGenes = go.annotated_genes && go.annotated_genes.includes('...and');
            const genesCell = hasMoreGenes ? 
                `<td class="genes-cell expandable-genes" data-go-id="${go.go_id}" data-full-genes="${go.annotated_genes_full || go.annotated_genes}">
                    <span class="genes-content">${go.annotated_genes}</span>
                    <span class="expand-icon">📋</span>
                </td>` :
                `<td class="genes-cell">${go.annotated_genes || 'Loading...'}</td>`;
            
            // Create download button for Guide RNA CSV
            const downloadButton = go.genes_with_guides > 0 ? 
                `<button class="btn btn-primary btn-download" onclick="window.tWheelReport.downloadGOGuideRNAs('${go.go_id}', '${go.go_name}')" title="Download Guide RNA CSV for ${go.go_id}">
                    📥 CSV
                </button>` :
                `<span class="text-muted">No guides</span>`;
            
            return `
                <tr>
                    <td><span class="go-id">${go.go_id}</span></td>
                    <td><strong>${go.go_name || 'Loading...'}</strong></td>
                    <td>${go.go_definition || 'Loading...'}</td>
                    <td><a href="${quickgoLink}" target="_blank" class="btn btn-outline">QuickGO</a></td>
                    <td>${go.total_genes}</td>
                    <td>${go.genes_with_guides}</td>
                    <td>${go.kfi_score.toFixed(4)}</td>
                    <td><span class="${targetableClass}">${go.targetable_type}</span></td>
                    ${genesCell}
                    <td>${downloadButton}</td>
                </tr>
            `;
        }).join('');
        
        goTable.innerHTML = rows;
        
        // Add click event listeners for expandable genes
        this.setupExpandableGenes();
        
        // Update pagination controls
        this.updateGOPaginationControls();
    }
    
    // GO Analysis pagination methods
    goToGOPreviousPage() {
        if (this.currentGOPage > 1) {
            this.currentGOPage--;
            this.renderFilteredGOData();
        }
    }
    
    goToGONextPage() {
        const totalPages = Math.ceil(this.filteredGOData.length / this.goPageSize);
        if (this.currentGOPage < totalPages) {
            this.currentGOPage++;
            this.renderFilteredGOData();
        }
    }
    
    updateGOPaginationControls() {
        const paginationControls = document.getElementById('go-pagination-controls');
        const prevBtn = document.getElementById('go-prev-page');
        const nextBtn = document.getElementById('go-next-page');
        const pageInfo = document.getElementById('go-page-info');
        const showingInfo = document.getElementById('go-showing-info');
        
        if (!paginationControls) return;
        
        const totalPages = Math.ceil(this.filteredGOData.length / this.goPageSize);
        const startIndex = (this.currentGOPage - 1) * this.goPageSize;
        const endIndex = Math.min(startIndex + this.goPageSize, this.filteredGOData.length);
        
        // Show/hide pagination controls
        paginationControls.style.display = totalPages > 1 ? 'flex' : 'none';
        
        // Update button states
        if (prevBtn) prevBtn.disabled = this.currentGOPage <= 1;
        if (nextBtn) nextBtn.disabled = this.currentGOPage >= totalPages;
        
        // Update page info
        if (pageInfo) pageInfo.textContent = `Page ${this.currentGOPage} of ${totalPages}`;
        
        // Update showing info
        if (showingInfo && this.filteredGOData.length > 0) {
            showingInfo.textContent = `Showing ${startIndex + 1}-${endIndex} of ${this.filteredGOData.length} entries`;
        } else if (showingInfo) {
            showingInfo.textContent = 'Showing 0-0 of 0 entries';
        }
    }
    
    setupExpandableGenes() {
        console.log('Setting up expandable genes...');
        // Add click event listeners for expandable gene cells
        const expandableCells = document.querySelectorAll('.expandable-genes');
        console.log(`Found ${expandableCells.length} expandable gene cells`);
        
        expandableCells.forEach((cell, index) => {
            console.log(`Setting up cell ${index}:`, {
                goId: cell.dataset.goId,
                hasFullGenes: !!cell.dataset.fullGenes,
                fullGenesLength: cell.dataset.fullGenes ? cell.dataset.fullGenes.length : 0
            });
            
            cell.addEventListener('click', (e) => {
                console.log('Expandable gene cell clicked:', cell.dataset.goId);
                
                try {
                    const goId = cell.dataset.goId;
                    const fullGenes = cell.dataset.fullGenes;
                    const genesContent = cell.querySelector('.genes-content');
                    
                    console.log('Click data:', {
                        goId,
                        hasFullGenes: !!fullGenes,
                        fullGenesPreview: fullGenes ? fullGenes.substring(0, 50) + '...' : 'null',
                        currentText: genesContent ? genesContent.textContent.substring(0, 50) + '...' : 'no content element'
                    });
                    
                    if (!goId || !fullGenes || !genesContent) {
                        console.error('Missing required elements:', { goId, fullGenes: !!fullGenes, genesContent: !!genesContent });
                        return;
                    }
                    
                    if (cell.classList.contains('expanded')) {
                        console.log('Collapsing genes list...');
                        // Collapse - show truncated version
                        const originalGO = this.filteredGOData.find(go => go.go_id === goId);
                        if (originalGO) {
                            genesContent.textContent = originalGO.annotated_genes;
                            cell.classList.remove('expanded');
                            cell.querySelector('.expand-icon').textContent = '📋';
                            cell.style.cursor = 'pointer';
                            console.log('Collapsed successfully');
                        } else {
                            console.error('Original GO data not found for:', goId);
                        }
                    } else {
                        console.log('Expanding genes list...');
                        // Expand - show full version
                        genesContent.textContent = fullGenes;
                        cell.classList.add('expanded');
                        cell.querySelector('.expand-icon').textContent = '📄';
                        cell.style.cursor = 'pointer';
                        console.log('Expanded successfully, new text length:', genesContent.textContent.length);
                    }
                } catch (error) {
                    console.error('Error in expandable genes click handler:', error);
                }
                
                e.stopPropagation();
            });
            
            // Add hover effects
            cell.style.cursor = 'pointer';
            cell.title = 'クリックして全ての遺伝子を表示/非表示';
        });
        
        console.log('Expandable genes setup completed');
    }
    
    loadGOCharts() {
        console.log('Loading GO charts...');
        try {
            // Draw GO-KFI bar chart
            this.drawGOKFIChart();
            
            // Load GO parent count images
            this.loadGOParentImages();
            
            // Load GO parent counts table
            this.loadGOParentTable();
            
        } catch (error) {
            console.error('Failed to load GO charts:', error);
        }
    }
    
    loadGOParentTable() {
        console.log('Loading GO parent counts table...');
        const parentTable = document.querySelector('#go-parent-table tbody');
        if (!parentTable) return;
        
        // Load embedded GO parent counts data
        if (this.data && this.data.go_parent_counts && Array.isArray(this.data.go_parent_counts)) {
            const parentData = this.data.go_parent_counts;
            
            if (parentData.length === 0) {
                parentTable.innerHTML = '<tr><td colspan="3" class="loading">No GO parent terms data available</td></tr>';
                return;
            }
            
            // Generate table rows
            const rows = parentData.map(parent => {
                return `
                    <tr>
                        <td><strong>${parent.parent_go}</strong></td>
                        <td>${parent.description}</td>
                        <td>${parent.count}</td>
                    </tr>
                `;
            }).join('');
            
            parentTable.innerHTML = rows;
            console.log(`Loaded ${parentData.length} GO parent terms entries`);
            
        } else {
            parentTable.innerHTML = '<tr><td colspan="3" class="error">GO parent terms data not found in embedded data</td></tr>';
        }
    }
    
    drawGOKFIChart() {
        const chartContainer = document.getElementById('go-kfi-chart');
        if (!chartContainer) return;
        
        // Get GO analysis statistics
        const goStats = this.data?.statistics?.go_analysis;
        if (!goStats || !goStats.chart_data) {
            chartContainer.innerHTML = '<p class="error">GO-KFI chart data not available</p>';
            return;
        }
        
        const chartData = goStats.chart_data;
        
        // Create horizontal bar chart using SVG
        const chartSVG = this.createGOKFIBarChart(chartData);
        
        chartContainer.innerHTML = '';
        chartContainer.appendChild(chartSVG);
        
        console.log('GO-KFI chart rendered');
    }
    
    createGOKFIBarChart(data) {
        const width = 900;
        const height = 300;
        const margin = { top: 80, right: 50, bottom: 70, left: 80 };
        
        // Create SVG element
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.setAttribute('width', width);
        svg.setAttribute('height', height);
        svg.style.fontFamily = 'Arial, sans-serif';
        svg.style.fontSize = '12px';
        
        // Chart dimensions
        const chartWidth = width - margin.left - margin.right;
        const chartHeight = height - margin.top - margin.bottom;
        
        // Bar position and height
        const barY = margin.top + (chartHeight - 60) / 2;
        const barHeight = 60;
        
        // Calculate segment widths based on percentages
        const simplyPercent = data.fully_targetable.simply_percent || 0;
        const multiplyPercent = data.fully_targetable.multiply_percent || 0;
        const partiallyPercent = data.partially_targetable.percent || 0;
        const notPercent = data.not_targetable.percent || 0;
        
        // Segment data with proper positioning
        const segments = [
            {
                label: 'Fully (F) and Simply Targetable (S)',
                percent: simplyPercent,
                count: Math.round(data.total * simplyPercent / 100),
                color: '#59BF10',
                x: 0
            },
            {
                label: 'Fully (F) and Multiply Targetable (M)',
                percent: multiplyPercent,
                count: Math.round(data.total * multiplyPercent / 100),
                color: '#D7DB01',
                x: (simplyPercent / 100) * chartWidth
            },
            {
                label: 'Partially Targetable (P)',
                percent: partiallyPercent,
                count: data.partially_targetable.count,
                color: '#F0E579',
                x: ((simplyPercent + multiplyPercent) / 100) * chartWidth
            },
            {
                label: 'Not Targetable (N)',
                percent: notPercent,
                count: data.not_targetable.count,
                color: '#D92000',
                x: ((simplyPercent + multiplyPercent + partiallyPercent) / 100) * chartWidth
            }
        ];
        
        // Title
        const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        title.setAttribute('x', width / 2);
        title.setAttribute('y', 20);
        title.setAttribute('text-anchor', 'middle');
        title.setAttribute('fill', '#333');
        title.setAttribute('font-weight', 'bold');
        title.setAttribute('font-size', '16px');
        title.textContent = 'GO-KFI Assessment Results';
        svg.appendChild(title);
        
        // Legend - arranged in 2 rows, 2 columns
        const legendStartX = 200;
        const legendY1 = 40;
        const legendY2 = 55;
        const legendSquareSize = 12;
        const legendSpacingX = 300;
        
        // First row - Simply and Multiply
        [segments[0], segments[1]].forEach((segment, i) => {
            const legendX = legendStartX + (i * legendSpacingX);
            
            // Legend square
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', legendX);
            rect.setAttribute('y', legendY1 - legendSquareSize / 2);
            rect.setAttribute('width', legendSquareSize);
            rect.setAttribute('height', legendSquareSize);
            rect.setAttribute('fill', segment.color);
            svg.appendChild(rect);
            
            // Legend text
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', legendX + legendSquareSize + 5);
            text.setAttribute('y', legendY1 + 3);
            text.setAttribute('fill', '#333');
            text.setAttribute('font-size', '11px');
            text.textContent = segment.label;
            svg.appendChild(text);
        });
        
        // Second row - Partially and Not
        [segments[2], segments[3]].forEach((segment, i) => {
            const legendX = legendStartX + (i * legendSpacingX);
            
            // Legend square
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', legendX);
            rect.setAttribute('y', legendY2 - legendSquareSize / 2);
            rect.setAttribute('width', legendSquareSize);
            rect.setAttribute('height', legendSquareSize);
            rect.setAttribute('fill', segment.color);
            svg.appendChild(rect);
            
            // Legend text
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', legendX + legendSquareSize + 5);
            text.setAttribute('y', legendY2 + 3);
            text.setAttribute('fill', '#333');
            text.setAttribute('font-size', '11px');
            text.textContent = segment.label;
            svg.appendChild(text);
        });
        
        // Y-axis label
        const yLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        yLabel.setAttribute('x', 25);
        yLabel.setAttribute('y', barY + barHeight / 2);
        yLabel.setAttribute('text-anchor', 'middle');
        yLabel.setAttribute('fill', '#333');
        yLabel.setAttribute('font-weight', 'bold');
        yLabel.textContent = 'Contigs';
        svg.appendChild(yLabel);
        
        // Draw stacked bar segments
        segments.forEach(segment => {
            if (segment.percent > 0) {
                const segmentWidth = (segment.percent / 100) * chartWidth;
                
                // Segment rectangle
                const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                rect.setAttribute('x', margin.left + segment.x);
                rect.setAttribute('y', barY);
                rect.setAttribute('width', segmentWidth);
                rect.setAttribute('height', barHeight);
                rect.setAttribute('fill', segment.color);
                rect.setAttribute('stroke', '#fff');
                rect.setAttribute('stroke-width', '1');
                svg.appendChild(rect);
            }
        });
        
        // Combined statistics text for first segment
        const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        text.setAttribute('x', 225);
        text.setAttribute('y', barY + barHeight / 2);
        text.setAttribute('text-anchor', 'middle');
        text.setAttribute('fill', '#333');
        text.setAttribute('font-weight', 'bold');
        text.setAttribute('font-size', '12px');
        text.textContent = `F:${data.fully_targetable.count} [S:${segments[0].count}, M:${segments[1].count}, P:${segments[2].count}, N:${segments[3].count}, n:${data.total}]`;
        svg.appendChild(text);
        
        // X-axis scale
        const scaleY = barY + barHeight + 20;
        const scaleMarks = [0, 20, 40, 60, 80, 100];
        
        // X-axis line
        const xAxisLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        xAxisLine.setAttribute('x1', margin.left);
        xAxisLine.setAttribute('y1', scaleY);
        xAxisLine.setAttribute('x2', margin.left + chartWidth);
        xAxisLine.setAttribute('y2', scaleY);
        xAxisLine.setAttribute('stroke', '#333');
        xAxisLine.setAttribute('stroke-width', '1');
        svg.appendChild(xAxisLine);
        
        // X-axis scale marks
        scaleMarks.forEach(mark => {
            const x = margin.left + (mark / 100) * chartWidth;
            
            // Scale mark
            const tick = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            tick.setAttribute('x1', x);
            tick.setAttribute('y1', scaleY);
            tick.setAttribute('x2', x);
            tick.setAttribute('y2', scaleY + 5);
            tick.setAttribute('stroke', '#333');
            tick.setAttribute('stroke-width', '1');
            svg.appendChild(tick);
            
            // Scale label
            const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            label.setAttribute('x', x);
            label.setAttribute('y', scaleY + 18);
            label.setAttribute('text-anchor', 'middle');
            label.setAttribute('fill', '#333');
            label.setAttribute('font-size', '10px');
            label.textContent = mark;
            svg.appendChild(label);
        });
        
        // X-axis label
        const xLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        xLabel.setAttribute('x', margin.left + chartWidth / 2);
        xLabel.setAttribute('y', height - 10);
        xLabel.setAttribute('text-anchor', 'middle');
        xLabel.setAttribute('fill', '#333');
        xLabel.setAttribute('font-weight', 'bold');
        xLabel.textContent = '% GO-KFI';
        svg.appendChild(xLabel);
        
        return svg;
    }
    
    loadGOParentImages() {
        const imagesContainer = document.getElementById('go-parent-images');
        if (!imagesContainer) return;
        
        // Clear loading message
        imagesContainer.innerHTML = '';
        
        // GO parent counts image
        const goParentImg = document.createElement('img');
        goParentImg.src = 'static/img/go_parent_counts.png';
        goParentImg.alt = 'GO Parent Counts Distribution';
        goParentImg.style.maxWidth = '70%';
        goParentImg.style.height = 'auto';
        goParentImg.style.border = '1px solid #dee2e6';
        goParentImg.style.borderRadius = '0.375rem';
        goParentImg.style.marginBottom = '1rem';
        goParentImg.style.display = 'block';
        goParentImg.style.margin = '0 auto 1rem auto';
        
        goParentImg.onload = () => {
            console.log('GO parent counts image loaded successfully');
        };
        
        goParentImg.onerror = () => {
            console.log('GO parent counts image not found, trying treemap image');
            goParentImg.style.display = 'none';
            
            // Try treemap image instead
            const treemapImg = document.createElement('img');
            treemapImg.src = 'static/img/treemap_go_parent_counts.png';
            treemapImg.alt = 'GO Parent Counts Treemap';
            treemapImg.style.maxWidth = '70%';
            treemapImg.style.height = 'auto';
            treemapImg.style.border = '1px solid #dee2e6';
            treemapImg.style.borderRadius = '0.375rem';
            treemapImg.style.display = 'block';
            treemapImg.style.margin = '0 auto';
            
            treemapImg.onload = () => {
                console.log('GO treemap image loaded successfully');
            };
            
            treemapImg.onerror = () => {
                console.log('GO treemap image not found');
                const noImageMsg = document.createElement('p');
                noImageMsg.textContent = 'GO parent distribution images not available';
                noImageMsg.style.color = '#6c757d';
                noImageMsg.style.fontStyle = 'italic';
                noImageMsg.style.textAlign = 'center';
                imagesContainer.appendChild(noImageMsg);
            };
            
            imagesContainer.appendChild(treemapImg);
        };
        
        imagesContainer.appendChild(goParentImg);
    }
    
    // Protein sequence viewer methods
    async loadProteinSequences() {
        try {
            // Use embedded protein data instead of fetching external file
            if (this.data && this.data.statistics && this.data.statistics.gene_prediction && this.data.statistics.gene_prediction.proteins) {
                this.proteins = this.data.statistics.gene_prediction.proteins;
                this.filteredProteins = [...this.proteins];
                
                // Update protein stats from embedded data
                if (this.data.statistics.gene_prediction.protein_stats) {
                    const proteinStats = this.data.statistics.gene_prediction.protein_stats;
                    this.updateStatCard('avg-protein-length', proteinStats.avg_length);
                    this.updateStatCard('longest-protein', proteinStats.max_length);
                    this.updateStatCard('total-proteins', proteinStats.total_count);
                }
                
                this.renderProteins();
                
                // Mark section as loaded
                this.sectionsLoaded.geneSequences = true;
                
            } else {
                throw new Error('No embedded protein data found');
            }
            
        } catch (error) {
            console.error('Error loading protein sequences:', error);
            this.showProteinError(`
                <h4>プロテイン情報の表示</h4>
                <p>軽量版のプロテインサマリーが実装されました。ヘッダー名、配列長、予測分子量が表示されます。</p>
                <p><strong>含まれる情報:</strong></p>
                <ul>
                    <li>プロテインヘッダー名</li>
                    <li>配列長（アミノ酸数）</li>
                    <li>長さによる分類（Short: ≤100 AA, Medium: 101-999 AA, Long: ≥1000 AA）</li>
                    <li>予測分子量（kDa）</li>
                </ul>
                <p>このバージョンでは配列表示はありませんが、CORSエラーなしに即座に利用できます。</p>
            `);
        }
    }
    
    parseFastaFile(fastaText) {
        const proteins = [];
        const lines = fastaText.split('\\n');
        let currentProtein = null;
        
        for (const line of lines) {
            const trimmed = line.trim();
            if (trimmed.startsWith('>')) {
                // Header line
                if (currentProtein) {
                    proteins.push(currentProtein);
                }
                currentProtein = {
                    id: trimmed.substring(1),
                    sequence: '',
                    length: 0
                };
            } else if (currentProtein && trimmed) {
                // Sequence line
                currentProtein.sequence += trimmed;
            }
        }
        
        // Add last protein
        if (currentProtein) {
            proteins.push(currentProtein);
        }
        
        // Calculate lengths and add metadata
        proteins.forEach(protein => {
            protein.length = protein.sequence.length;
            protein.lengthCategory = this.categorizeProteinLength(protein.length);
            protein.gcContent = this.calculateGCContent(protein.sequence);
        });
        
        return proteins;
    }
    
    categorizeProteinLength(length) {
        if (length < 100) return 'short';
        if (length <= 1000) return 'medium';
        return 'long';
    }
    
    calculateGCContent(sequence) {
        // For protein sequences, we'll calculate hydrophobic vs hydrophilic ratio instead
        const hydrophobic = 'AILMFPWYVQ';
        let hydrophobicCount = 0;
        
        for (const aa of sequence) {
            if (hydrophobic.includes(aa.toUpperCase())) {
                hydrophobicCount++;
            }
        }
        
        return sequence.length > 0 ? (hydrophobicCount / sequence.length * 100).toFixed(1) : 0;
    }
    
    calculateProteinStats() {
        if (this.proteins.length === 0) return;
        
        const lengths = this.proteins.map(p => p.length);
        const avgLength = lengths.reduce((a, b) => a + b, 0) / lengths.length;
        const maxLength = Math.max(...lengths);
        
        this.updateStatCard('avg-protein-length', Math.round(avgLength));
        this.updateStatCard('longest-protein', maxLength);
        this.updateStatCard('total-proteins', this.proteins.length);
    }
    
    renderProteins() {
        const container = document.getElementById('protein-container');
        const paginationControls = document.getElementById('pagination-controls');
        
        if (!container) return;
        
        if (this.filteredProteins.length === 0) {
            container.innerHTML = '<div class="no-results">No proteins found matching the current filters.</div>';
            paginationControls.style.display = 'none';
            return;
        }
        
        // Calculate pagination
        const totalPages = Math.ceil(this.filteredProteins.length / this.pageSize);
        const startIndex = (this.currentPage - 1) * this.pageSize;
        const endIndex = Math.min(startIndex + this.pageSize, this.filteredProteins.length);
        const currentProteins = this.filteredProteins.slice(startIndex, endIndex);
        
        // Render proteins
        const proteinElements = currentProteins.map(protein => this.createProteinElement(protein));
        container.innerHTML = proteinElements.join('');
        
        // Update pagination controls
        this.updatePaginationControls(totalPages);
        
        // Add event listeners to protein headers
        this.addProteinEventListeners();
    }
    
    createProteinElement(protein) {
        return `
            <div class="protein-item">
                <div class="protein-header" data-length="${protein.category}" data-protein-id="${protein.id}">
                    <div class="protein-id">${protein.id}</div>
                    <div class="protein-stats">
                        <span>${protein.length} AA</span>
                        <span>${protein.mass_kda} kDa</span>
                        <span class="category-badge category-${protein.category}">${protein.category.toUpperCase()}</span>
                    </div>
                </div>
            </div>
        `;
    }
    
    formatSequence(sequence, lineLength = 60) {
        const formatted = [];
        for (let i = 0; i < sequence.length; i += lineLength) {
            formatted.push(sequence.substring(i, i + lineLength));
        }
        return formatted.join('\\n');
    }
    
    addProteinEventListeners() {
        document.querySelectorAll('.protein-header').forEach(header => {
            header.addEventListener('click', (e) => {
                const proteinId = header.dataset.proteinId;
                const content = document.querySelector(`[data-protein-id="${proteinId}"].protein-content`);
                const toggle = header.querySelector('.expand-toggle');
                
                if (content.classList.contains('expanded')) {
                    content.classList.remove('expanded');
                    header.classList.remove('expanded');
                    toggle.textContent = '+';
                } else {
                    content.classList.add('expanded');
                    header.classList.add('expanded');
                    toggle.textContent = '−';
                }
            });
        });
    }
    
    updatePaginationControls(totalPages) {
        const paginationControls = document.getElementById('pagination-controls');
        const prevBtn = document.getElementById('prev-page');
        const nextBtn = document.getElementById('next-page');
        const pageInfo = document.getElementById('page-info');
        
        if (!paginationControls) return;
        
        paginationControls.style.display = totalPages > 1 ? 'flex' : 'none';
        
        if (prevBtn) prevBtn.disabled = this.currentPage <= 1;
        if (nextBtn) nextBtn.disabled = this.currentPage >= totalPages;
        if (pageInfo) pageInfo.textContent = `Page ${this.currentPage} of ${totalPages}`;
    }
    
    // Navigation methods
    goToPreviousPage() {
        if (this.currentPage > 1) {
            this.currentPage--;
            this.renderProteins();
        }
    }
    
    goToNextPage() {
        const totalPages = Math.ceil(this.filteredProteins.length / this.pageSize);
        if (this.currentPage < totalPages) {
            this.currentPage++;
            this.renderProteins();
        }
    }
    
    
    // Expand/Collapse methods
    expandAllProteins() {
        document.querySelectorAll('.protein-content').forEach(content => {
            content.classList.add('expanded');
        });
        document.querySelectorAll('.protein-header').forEach(header => {
            header.classList.add('expanded');
        });
        document.querySelectorAll('.expand-toggle').forEach(toggle => {
            toggle.textContent = '−';
        });
    }
    
    collapseAllProteins() {
        document.querySelectorAll('.protein-content').forEach(content => {
            content.classList.remove('expanded');
        });
        document.querySelectorAll('.protein-header').forEach(header => {
            header.classList.remove('expanded');
        });
        document.querySelectorAll('.expand-toggle').forEach(toggle => {
            toggle.textContent = '+';
        });
    }
    
    // Utility methods for protein viewer
    copySequence(proteinId) {
        const protein = this.proteins.find(p => p.id === proteinId);
        if (protein) {
            this.copyToClipboard(protein.sequence);
        }
    }
    
    copyHeader(proteinId) {
        this.copyToClipboard(proteinId);
    }
    
    showProteinError(message) {
        const container = document.getElementById('protein-container');
        if (container) {
            container.innerHTML = `<div class="error">${message}</div>`;
        }
    }
    
    handleSearch(query, target) {
        console.log('handleSearch called:', query, target);
        if (target === 'proteins') {
            this.searchProteins(query);
            return;
        }
        
        // Handle other search targets
        console.log(`Searching for "${query}" in ${target}`);
    }
    
    handleFilter(value, filterType) {
        console.log('handleFilter called:', value, filterType);
        if (filterType === 'length') {
            this.filterProteinsByLength(value);
            return;
        }
        
        console.log(`Filtering by ${filterType}: ${value}`);
    }
    
    searchProteins(query) {
        console.log('searchProteins called with:', query);
        console.log('Total proteins:', this.proteins.length);
        
        if (!query.trim()) {
            this.filteredProteins = [...this.proteins];
        } else {
            const lowerQuery = query.toLowerCase();
            this.filteredProteins = this.proteins.filter(protein => 
                protein.id.toLowerCase().includes(lowerQuery)
            );
        }
        
        console.log('Filtered proteins:', this.filteredProteins.length);
        this.currentPage = 1;
        this.renderProteins();
    }
    
    filterProteinsByLength(category) {
        console.log('filterProteinsByLength called with:', category);
        console.log('Total proteins:', this.proteins.length);
        
        if (!category) {
            this.filteredProteins = [...this.proteins];
        } else {
            this.filteredProteins = this.proteins.filter(protein => {
                console.log('Protein category:', protein.category, 'Filter:', category);
                return protein.category === category;
            });
        }
        
        console.log('Filtered proteins:', this.filteredProteins.length);
        this.currentPage = 1;
        this.renderProteins();
    }
    
    // Utility methods
    formatNumber(num) {
        if (typeof num !== 'number') return num;
        return num.toLocaleString();
    }
    
    formatPercentage(num, decimals = 1) {
        if (typeof num !== 'number') return '0.0%';
        return num.toFixed(decimals) + '%';
    }
    
    // Guide RNA CSV download functionality
    async downloadGOGuideRNAs(goId, goName) {
        console.log('Downloading Guide RNAs for GO:', goId, goName);
        
        try {
            // Find GO term data
            const goTerm = this.filteredGOData.find(go => go.go_id === goId);
            if (!goTerm) {
                console.error('GO term not found:', goId);
                return;
            }
            
            // Get genes for this GO term
            const goGenes = goTerm.annotated_genes_full || goTerm.annotated_genes;
            if (!goGenes) {
                console.error('No genes found for GO term:', goId);
                this.showError('No genes found for this GO term');
                return;
            }
            
            // Parse gene names from the genes string
            const geneNames = goGenes.split(',').map(gene => gene.trim()).filter(gene => gene.length > 0);
            console.log('Found genes for GO term:', geneNames);
            
            // Find guide RNAs for these genes
            const relevantGuides = [];
            
            if (this.data.guide_rna && Array.isArray(this.data.guide_rna)) {
                for (const guide of this.data.guide_rna) {
                    // Check if this guide belongs to any of the GO term's genes
                    // We need to match by gene name or query_id
                    const matchingGene = geneNames.find(geneName => {
                        return guide.reference.includes(geneName);
                    });
                    
                    if (matchingGene) {
                        relevantGuides.push({
                            ...guide,
                            matched_gene: matchingGene,
                            go_id: goId,
                            go_term: goName
                        });
                    }
                }
            }
            
            console.log(`Found ${relevantGuides.length} guide RNAs for GO term ${goId}`);
            
            if (relevantGuides.length === 0) {
                this.showError('No guide RNAs found for genes in this GO term');
                return;
            }
            
            // Generate CSV content
            const csvContent = this.generateGOGuideRNACSV(relevantGuides, goId, goName);
            
            // Download CSV file
            const filename = `GO_${goId}_GuideRNAs_${new Date().toISOString().slice(0, 10)}.csv`;
            this.downloadCSV(csvContent, filename);
            
            this.showSuccess(`Downloaded ${relevantGuides.length} guide RNAs for GO term ${goId}`);
            
        } catch (error) {
            console.error('Error downloading GO guide RNAs:', error);
            this.showError('Failed to download guide RNA data');
        }
    }
    
    generateGOGuideRNACSV(guides, goId, goName) {
        // CSV header
        const headers = [
            'GO_ID',
            'GO_Term',
            'Matched_Gene_Group',
            'Query_ID',
            'Reference',
            'Start',
            'End',
            'Guide_Sequence',
            'Strand',
            'MM0',
            'MM1_2',
            'Off_Target_Score'
        ];
        
        // Create CSV rows
        const rows = guides.map(guide => [
            goId,
            `"${goName}"`, // Quote to handle commas in GO names
            guide.matched_gene,
            guide.query_id,
            guide.reference,
            guide.start,
            guide.end,
            guide.target_sequence,
            guide.strand,
            guide.mm0,
            guide.mm1_2,
            guide.off_target_score
        ]);
        
        // Combine header and rows
        const csvLines = [headers.join(',')].concat(rows.map(row => row.join(',')));
        
        return csvLines.join('\\n');
    }
    
    downloadCSV(csvContent, filename) {
        // Create blob and download
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        
        if (link.download !== undefined) {
            // Modern browsers
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', filename);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            URL.revokeObjectURL(url);
        } else {
            // Fallback for older browsers
            console.error('CSV download not supported in this browser');
            this.showError('CSV download not supported in this browser');
        }
    }

    // Gene Annotation CSV download functionality
    downloadGeneAnnotationCSV() {
        console.log('Downloading Gene Annotation CSV...');
        
        try {
            if (!this.filteredGeneData || this.filteredGeneData.length === 0) {
                this.showError('No gene annotation data available for download');
                return;
            }
            
            // CSV header
            const headers = [
                'Query_ID',
                'UniProt_Name',
                'Gene_Name',
                'Description',
                'tWhEEL_Classification',
                'Gene_Group',
                'Number_of_Targets'
            ];
            
            // Create CSV rows
            const rows = this.filteredGeneData.map(gene => [
                gene.query_id,
                gene.subject_id,
                `"${gene.gene_name}"`, // Quote to handle commas
                `"${gene.description}"`, // Quote to handle commas
                gene.classification,
                gene.unique_id,
                gene.target_count
            ]);
            
            // Combine header and rows
            const csvLines = [headers.join(',')].concat(rows.map(row => row.join(',')));
            const csvContent = csvLines.join('\\n');
            
            // Download CSV file
            const filename = `Gene_Annotation_${new Date().toISOString().slice(0, 10)}.csv`;
            this.downloadCSV(csvContent, filename);
            
            this.showSuccess(`Downloaded ${this.filteredGeneData.length} gene annotation entries`);
            
        } catch (error) {
            console.error('Error downloading gene annotation CSV:', error);
            this.showError('Failed to download gene annotation data');
        }
    }

    // Guide RNA CSV download functionality
    downloadGuideRNACSV() {
        console.log('Downloading Guide RNA CSV...');
        
        try {
            if (!this.filteredGuideData || this.filteredGuideData.length === 0) {
                this.showError('No guide RNA data available for download');
                return;
            }
            
            // CSV header
            const headers = [
                'Query_ID',
                'Reference',
                'Start',
                'End',
                'Target_Sequence',
                'Strand',
                'MM0',
                'MM1_2',
                'Off_Target_Score'
            ];
            
            // Create CSV rows
            const rows = this.filteredGuideData.map(guide => [
                guide.query_id,
                guide.reference,
                guide.start,
                guide.end,
                guide.target_sequence,
                guide.strand,
                guide.mm0,
                guide.mm1_2,
                guide.off_target_score
            ]);
            
            // Combine header and rows
            const csvLines = [headers.join(',')].concat(rows.map(row => row.join(',')));
            const csvContent = csvLines.join('\\n');
            
            // Download CSV file
            const filename = `Guide_RNA_${new Date().toISOString().slice(0, 10)}.csv`;
            this.downloadCSV(csvContent, filename);
            
            this.showSuccess(`Downloaded ${this.filteredGuideData.length} guide RNA entries`);
            
        } catch (error) {
            console.error('Error downloading guide RNA CSV:', error);
            this.showError('Failed to download guide RNA data');
        }
    }

    // GO Analysis CSV download functionality
    downloadGOAnalysisCSV() {
        console.log('Downloading GO Analysis CSV...');
        
        try {
            if (!this.filteredGOData || this.filteredGOData.length === 0) {
                this.showError('No GO analysis data available for download');
                return;
            }
            
            // CSV header
            const headers = [
                'GO_ID',
                'GO_Term',
                'GO_Definition',
                'GO_annotated_Gene_Count',
                'Targetable_GO_annotated_Gene_Count',
                'GO_KFI',
                'GO_KFI_Class',
                'GO_annotated_Gene_Groups'
            ];
            
            // Create CSV rows
            const rows = this.filteredGOData.map(go => [
                go.go_id,
                `"${go.go_name}"`, // Quote to handle commas
                `"${go.go_definition}"`, // Quote to handle commas  
                go.total_genes,
                go.genes_with_guides,
                go.kfi_score.toFixed(4),
                go.targetable_type,
                `"${go.annotated_genes_full || go.annotated_genes}"` // Quote to handle commas
            ]);
            
            // Combine header and rows
            const csvLines = [headers.join(',')].concat(rows.map(row => row.join(',')));
            const csvContent = csvLines.join('\\n');
            
            // Download CSV file
            const filename = `GO_Analysis_${new Date().toISOString().slice(0, 10)}.csv`;
            this.downloadCSV(csvContent, filename);
            
            this.showSuccess(`Downloaded ${this.filteredGOData.length} GO analysis entries`);
            
        } catch (error) {
            console.error('Error downloading GO analysis CSV:', error);
            this.showError('Failed to download GO analysis data');
        }
    }

    // GO Parent Terms CSV download functionality
    downloadGOParentTermsCSV() {
        console.log('Downloading GO Parent Terms CSV...');
        
        try {
            if (!this.data || !this.data.go_parent_counts || this.data.go_parent_counts.length === 0) {
                this.showError('No GO parent terms data available for download');
                return;
            }
            
            // CSV header
            const headers = [
                'Parent_GO_Term',
                'Description',
                'Count'
            ];
            
            // Create CSV rows
            const rows = this.data.go_parent_counts.map(parent => [
                parent.parent_go,
                `"${parent.description}"`, // Quote to handle commas
                parent.count
            ]);
            
            // Combine header and rows
            const csvLines = [headers.join(',')].concat(rows.map(row => row.join(',')));
            const csvContent = csvLines.join('\\n');
            
            // Download CSV file
            const filename = `GO_Parent_Terms_${new Date().toISOString().slice(0, 10)}.csv`;
            this.downloadCSV(csvContent, filename);
            
            this.showSuccess(`Downloaded ${this.data.go_parent_counts.length} GO parent terms entries`);
            
        } catch (error) {
            console.error('Error downloading GO parent terms CSV:', error);
            this.showError('Failed to download GO parent terms data');
        }
    }

    // Copy to clipboard functionality (for sequences, IDs, etc.)
    async copyToClipboard(text) {
        try {
            await navigator.clipboard.writeText(text);
            this.showSuccess('Copied to clipboard');
        } catch (error) {
            console.error('Failed to copy to clipboard:', error);
        }
    }
}

// Utility functions
function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

// Application initialization
document.addEventListener('DOMContentLoaded', () => {
    window.tWheelReport = new tWhEELReport();
});

// Export for potential module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = tWhEELReport;
}
""")
    
    def _generate_html_content(self, results):
        """Generate HTML content"""
        stats = results['statistics']
        
        # Load Gene Annotation data for embedding
        gene_annotation_data = self._load_gene_annotation_data_for_js(results)
        # Load Guide RNA data for embedding
        guide_rna_data = self._load_guide_rna_data_for_js(results)
        # Load GO analysis data for embedding
        go_analysis_data = self._load_go_data_for_js(results)
        # Load GO parent counts data for embedding
        go_parent_counts_data = self._load_go_parent_counts_data_for_js(results)
        
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Transcriptome Functional Analysis Report</title>
    <link rel="stylesheet" href="static/css/main.css">
</head>
<body>
    <!-- Embedded Data -->
    <script type="application/json" id="twheel-data">{json.dumps({
        'timestamp': results['timestamp'],
        'statistics': results['statistics'],
        'config': {
            'input': results['config']['input'],
            'crispr': results['config']['crispr'],
            'busco': results['config']['busco']
        },
        'gene_annotation': gene_annotation_data,
        'guide_rna': guide_rna_data,
        'go_analysis': go_analysis_data,
        'go_parent_counts': go_parent_counts_data
    }, ensure_ascii=False)}</script>

    <div class="report-layout">
        <!-- Sidebar -->
        <aside class="sidebar">
            <div class="sidebar-header">
                <img src="static/img/logo.png" alt="tWhEEL Logo">
                <p class="subtitle">A pipeline to aggregate results from gene functional analysis, guide RNA design, and gene ontology (GO)-based knock-out feasibility assessment into a single report.</p>
                <p class="meta">Report generated on {results['timestamp'][:19].replace('T', ' ')} UTC based on data in: intermediate/</p>
            </div>
            <nav class="sidebar-nav">
                <ul>
                    <li><a href="#general-stats" class="nav-link active">General Statistics</a></li>
                    <li><a href="#gene-prediction" class="nav-link">Gene Prediction</a></li>
                    <li><a href="#gene-annotation" class="nav-link">Gene Annotation</a></li>
                    <li><a href="#guide-rna" class="nav-link">Guide RNA Design</a></li>
                    <li><a href="#go-analysis" class="nav-link">GO-KFI Assessment</a></li>
                    <li><a href="#configuration" class="nav-link">Configuration</a></li>
                </ul>
            </nav>
        </aside>

        <!-- Main Content -->
        <main class="main-content">
            <header class="content-header">
                <h1>tWhEEL</h1>
                <p class="subtitle">Transcriptome Workflow for Hierarchical Exploration and Expression Levels</p>
            </header>

            <div class="content-body">
                <!-- General Statistics Section -->
                <section id="general-stats" class="content-section">
                    <div class="section-header">
                        <h2>General Statistics</h2>
                    </div>
                    <div class="section-body">
                        <div class="alert alert-info">
                            <strong>Welcome!</strong> Not sure where to start? 
                            <a href="#go-analysis">Jump to GO-KFI Assessment</a> section to see detailed results.
                        </div>

                        <div class="table-responsive">
                            <table class="data-table">
                                <thead>
                                    <tr>
                                        <th></th>
                                        <th>Total transcripts</th>
                                        <th>Protein Annotated</th>
                                        <th>Total Annotated Genes</th>
                                        <th>Multi (Redundant)</th>
                                        <th>Single (Non-Redundant)</th>
                                        <th>GO Terms</th>
                                        <th>Targetable GO Terms [%]</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td>Contigs</td>
                                        <td>{stats.get('gene_prediction', {}).get('transcript_count', 0):,}</td>
                                        <td>{stats.get('gene_prediction', {}).get('protein_count', 0):,}</td>
                                        <td>{stats.get('annotation', {}).get('total_genes', 0):,}</td>
                                        <td>{stats.get('annotation', {}).get('multi_genes', 0):,}</td>
                                        <td>{stats.get('annotation', {}).get('single_genes', 0):,}</td>
                                        <td>{stats.get('go_analysis', {}).get('total_gos', 0):,}</td>
                                        <td>{stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>
                    </div>
                </section>

                <!-- Gene Prediction Section -->
                <section id="gene-prediction" class="content-section">
                    <div class="section-header">
                        <h2>Gene Prediction</h2>
                    </div>
                    <div class="section-body">
                        <p>TransDecoder was used to predict protein-coding genes from the transcriptome assembly.</p>
                        
                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value" id="total-proteins">
                                    {stats.get('gene_prediction', {}).get('protein_count', 0):,}
                                </div>
                                <div class="stat-label">Total Predicted Proteins</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value" id="avg-protein-length">
                                    <span class="loading">...</span>
                                </div>
                                <div class="stat-label">Average Length (AA)</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value" id="longest-protein">
                                    <span class="loading">...</span>
                                </div>
                                <div class="stat-label">Longest Protein (AA)</div>
                            </div>
                        </div>

                        <div class="filter-group">
                            <input type="text" class="search-box" id="protein-search" placeholder="Search protein headers..." data-target="proteins">
                            <select class="filter-select" id="length-filter" data-filter="length">
                                <option value="">All Lengths</option>
                                <option value="short">Short (≤100 AA)</option>
                                <option value="medium">Medium (101-999 AA)</option>
                                <option value="long">Long (≥1000 AA)</option>
                            </select>
                        </div>

                        <div id="protein-container">
                            <div class="loading-spinner">
                                <p>Loading protein sequences...</p>
                            </div>
                        </div>

                        <div id="pagination-controls" class="pagination-controls" style="display: none;">
                            <button id="prev-page" class="btn btn-outline">Previous</button>
                            <span id="page-info">Page 1 of 1</span>
                            <button id="next-page" class="btn btn-outline">Next</button>
                            <select id="page-size" class="filter-select">
                                <option value="10">10 per page</option>
                                <option value="50">50 per page</option>
                                <option value="100">100 per page</option>
                                <option value="200">200 per page</option>
                            </select>
                        </div>
                    </div>
                </section>

                <!-- Gene Annotation Section -->
                <section id="gene-annotation" class="content-section">
                    <div class="section-header">
                        <h2>Gene Annotation</h2>
                    </div>
                    <div class="section-body">
                        <p>Protein sequences were searched against UniProt Knowledgebase (UniProtKB) database using FASTA for functional annotation.</p>     
                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('annotation', {}).get('total_genes', 0):,}</div>
                                <div class="stat-label">Total Annotated Contigs</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('annotation', {}).get('multi_genes', 0):,}</div>
                                <div class="stat-label">Multi (Redundant) Contigs</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('annotation', {}).get('single_genes', 0):,}</div>
                                <div class="stat-label">Single (Non-Redundant) Contigs</div>
                            </div>
                        </div>
                        <div class="filter-group">
                            <input type="text" class="search-box" placeholder="Search by query ID, gene name, classification, or gene group ID..." data-target="genes">
                            <select class="filter-select" data-filter="classification">
                                <option value="">All Classifications</option>
                                <option value="Multi">Multi Genes</option>
                                <option value="Single">Single Genes</option>
                            </select>
                            <select class="filter-select" data-filter="guides">
                                <option value="">Guide RNA Status</option>
                                <option value="with-guides">With Guide RNA</option>
                                <option value="without-guides">Without Guide RNA</option>
                            </select>
                            <button class="btn btn-primary" onclick="window.tWheelReport.downloadGeneAnnotationCSV()" title="Download Gene Annotation CSV">
                                📥 Download CSV
                            </button>
                        </div>

                        <div class="table-responsive">
                            <table class="data-table" id="genes-table">
                                <thead>
                                    <tr>
                                        <th>Query ID</th>
                                        <th>UniProt Name</th>
                                        <th>Gene Name</th>
                                        <th>Description</th>
                                        <th>tWhEEL Classification<br>Gene Group</th>
                                        <th>Number of Targets</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td colspan="6" class="loading">Loading gene data...</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>

                        <div id="genes-pagination-controls" class="pagination-controls" style="display: none;">
                            <div class="pagination-info">
                                <span id="genes-showing-info">Showing 0-0 of 0 entries</span>
                            </div>
                            <div class="pagination-buttons">
                                <button id="genes-prev-page" class="btn btn-outline">Previous</button>
                                <span id="genes-page-info">Page 1 of 1</span>
                                <button id="genes-next-page" class="btn btn-outline">Next</button>
                            </div>
                            <div class="pagination-size">
                                <select id="genes-page-size" class="filter-select">
                                    <option value="10" selected>10 per page</option>
                                    <option value="25">25 per page</option>
                                    <option value="50">50 per page</option>
                                    <option value="100">100 per page</option>
                                </select>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- Guide RNA Design Section -->
                <section id="guide-rna" class="content-section">
                    <div class="section-header">
                        <h2>Guide RNA Design</h2>
                    </div>
                    <div class="section-body">
                        <p>CRISPR guide RNAs were designed for gene knockout using CRISPRFlash with {results['config']['crispr'].get('pam_pattern', 'NGG')}-3′ PAM motif.</p>

                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('guide_design', {}).get('total_guides', 0):,}</div>
                                <div class="stat-label">Total Guides</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('guide_design', {}).get('site_specific_guides', 0):,}</div>
                                <div class="stat-label">Site-Specific Guides</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('guide_design', {}).get('site_specific_multi_target', 0):,}</div>
                                <div class="stat-label">Site-Specific Multi-Target Guides</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value">{stats.get('guide_design', {}).get('site_specific_single_target', 0):,}</div>
                                <div class="stat-label">Site-Specific Single-Target Guides</div>
                            </div>
                        </div>

                        <div class="filter-group">
                            <input type="text" class="search-box" placeholder="Search by query ID, reference, or target sequence..." data-target="guides">
                            <select class="filter-select" data-filter="mm0">
                                <option value="">All MM0</option>
                                <option value="0">MM0 = 0</option>
                                <option value="1">MM0 = 1</option>
                                <option value="2">MM0 ≥ 2</option>
                            </select>
                            <select class="filter-select" data-filter="mm1-2">
                                <option value="">All MM1-2</option>
                                <option value="0">MM1-2 = 0</option>
                                <option value="1">MM1-2 = 1</option>
                                <option value="2">MM1-2 ≥ 2</option>
                            </select>
                            <button class="btn btn-primary" onclick="window.tWheelReport.downloadGuideRNACSV()" title="Download Guide RNA CSV">
                                📥 Download CSV
                            </button>
                        </div>

                        <div class="table-responsive">
                            <table class="data-table" id="guides-table">
                                <thead>
                                    <tr>
                                        <th>Query ID</th>
                                        <th>Reference</th>
                                        <th>Start</th>
                                        <th>End</th>
                                        <th>Target Sequence</th>
                                        <th>Strand</th>
                                        <th>MM0</th>
                                        <th>MM1-2</th>
                                        <th>Off-target Match Score</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td colspan="9" class="loading">Loading guide RNA data...</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>

                        <div id="guides-pagination-controls" class="pagination-controls" style="display: none;">
                            <div class="pagination-info">
                                <span id="guides-showing-info">Showing 0-0 of 0 entries</span>
                            </div>
                            <div class="pagination-buttons">
                                <button id="guides-prev-page" class="btn btn-outline">Previous</button>
                                <span id="guides-page-info">Page 1 of 1</span>
                                <button id="guides-next-page" class="btn btn-outline">Next</button>
                            </div>
                            <div class="pagination-size">
                                <select id="guides-page-size" class="filter-select">
                                    <option value="10" selected>10 per page</option>
                                    <option value="25">25 per page</option>
                                    <option value="50">50 per page</option>
                                    <option value="100">100 per page</option>
                                </select>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- GO-KFI Assessment Section -->
                <section id="go-analysis" class="content-section">
                    <div class="section-header">
                        <h2>GO-KFI Assessment</h2>
                    </div>
                    <div class="section-body">
                        <p>Gene Ontology functional analysis with knock-out feasibility assessment (GO-KFI assessment).</p>

                        <!-- GO Terms Statistics Cards -->
                        <div class="stats-grid">
                            <div class="stat-card">
                                <div class="stat-value stat-value-fully-targetable">{stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%</div>
                                <div class="stat-label">Fully Targetable GO terms</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value stat-value-simply-targetable">{stats.get('go_analysis', {}).get('simply_targetable_percent', 0):.1f}%</div>
                                <div class="stat-label">Fully and Simply Targetable GO terms</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value stat-value-multiply-targetable">{stats.get('go_analysis', {}).get('multiply_targetable_percent', 0):.1f}%</div>
                                <div class="stat-label">Fully and Multiply Targetable GO terms</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value stat-value-partially-targetable">{stats.get('go_analysis', {}).get('partially_targetable_percent', 0):.1f}%</div>
                                <div class="stat-label">Partially Targetable GO terms</div>
                            </div>
                            
                            <div class="stat-card">
                                <div class="stat-value stat-value-not-targetable">{stats.get('go_analysis', {}).get('not_targetable_percent', 0):.1f}%</div>
                                <div class="stat-label">Not Targetable GO terms</div>
                            </div>           
                        </div>

                        <!-- GO-KFI Distribution Chart Section -->
                        <div class="chart-container" id="go-kfi-chart-container" style="height: 350px; margin-bottom: 2rem;">
                            <div id="go-kfi-chart">
                                <p class="loading">Loading GO-KFI distribution chart...</p>
                            </div>
                        </div>

                        <div class="filter-group">
                            <input type="text" class="search-box" placeholder="Search by GO ID, term, or genes..." data-target="go-terms">
                            <select class="filter-select" data-filter="targetable">
                                <option value="">All GO-KFI Classes</option>
                                <option value="Fully">Fully targetable</option>
                                <option value="Partially">Partially targetable</option>
                                <option value="Not">Not targetable</option>
                            </select>
                            <button class="btn btn-primary" onclick="window.tWheelReport.downloadGOAnalysisCSV()" title="Download GO Analysis CSV">
                                📥 Download CSV
                            </button>
                        </div>

                        <div class="table-responsive">
                            <table class="data-table" id="go-table">
                                <thead>
                                    <tr>
                                        <th>GO ID</th>
                                        <th>GO Term</th>
                                        <th>GO Definition</th>
                                        <th>QuickGO</th>
                                        <th>GO-annotated Gene Count</th>
                                        <th>Targetable GO-annotated Gene Count</th>
                                        <th>GO-KFI</th>
                                        <th>GO-KFI Class</th>
                                        <th>GO-annotated Gene Groups</th>
                                        <th>Guide RNA CSV</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td colspan="10" class="loading">Loading GO-KFI assessment data...</td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>

                        <div id="go-pagination-controls" class="pagination-controls" style="display: none;">
                            <div class="pagination-info">
                                <span id="go-showing-info">Showing 0-0 of 0 entries</span>
                            </div>
                            <div class="pagination-buttons">
                                <button id="go-prev-page" class="btn btn-outline">Previous</button>
                                <span id="go-page-info">Page 1 of 1</span>
                                <button id="go-next-page" class="btn btn-outline">Next</button>
                            </div>
                            <div class="pagination-size">
                                <select id="go-page-size" class="filter-select">
                                    <option value="25" selected>25 per page</option>
                                    <option value="50">50 per page</option>
                                    <option value="100">100 per page</option>
                                    <option value="200">200 per page</option>
                                </select>
                            </div>
                        </div>

                        <!-- GO Parent Terms Distribution Section -->
                        <div style="margin-top: 2rem;">
                            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
                                <h4 style="margin: 0; color: #333;">GO Parent Terms Distribution</h4>
                                <p style="margin: 0; color: #333;">Search and count of Parent GOs for GOs that are considered Fully Targetable GO terms</p>
                                <button class="btn btn-primary" onclick="window.tWheelReport.downloadGOParentTermsCSV()" title="Download GO Parent Terms CSV">
                                    📥 Download CSV
                                </button>
                            </div>
                            
                            <!-- GO Parent Terms Images -->
                            <div class="chart-container" id="go-parent-images-container" style="height: auto; flex-direction: column; text-align: center; margin-bottom: 2rem;">
                                <div id="go-parent-images">
                                    <p class="loading">Loading GO parent terms distribution images...</p>
                                </div>
                            </div>
                            
                            <!-- GO Parent Terms Table -->
                            <div class="table-responsive">
                                <table class="data-table" id="go-parent-table">
                                    <thead>
                                        <tr>
                                            <th>Parent GO Term</th>
                                            <th>Description</th>
                                            <th>Count</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td colspan="3" class="loading">Loading GO parent terms distribution data...</td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </section>

                <!-- Configuration Section -->
                <section id="configuration" class="content-section">
                    <div class="section-header">
                        <h2>Analysis Configuration</h2>
                    </div>
                    <div class="section-body">
                        <div class="stats-grid">
                            <div class="stat-card">
                                <h3>Input Files</h3>
                                <ul style="text-align: left; font-size: 0.85rem;">
                                    <li><strong>Genome Assembly:</strong><br>{results['config']['input'].get('genome_assembly', 'N/A')}</li>
                                    <li><strong>Transcriptome:</strong><br>{results['config']['input'].get('transcriptome', 'N/A')}</li>
                                    <li><strong>UniProt Database:</strong><br>{results['config']['input'].get('uniprot_db', 'N/A')}</li>
                                </ul>
                            </div>
                            
                            <div class="stat-card">
                                <h3>CRISPR Settings</h3>
                                <ul style="text-align: left; font-size: 0.85rem;">
                                    <li><strong>PAM Pattern:</strong> {results['config']['crispr'].get('pam_pattern', 'N/A')}</li>
                                    <li><strong>Max Mismatch:</strong> {results['config']['crispr'].get('max_mismatch', 'N/A')}</li>
                                    <li><strong>Threads:</strong> {results['config']['crispr'].get('threads', 'N/A')}</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </section>
            </div>
        </main>
    </div>
    
    <script src="static/js/main.js"></script>
</body>
</html>
"""
        
        return html_content
    
    def _generate_static_html_content(self, results):
        """Generate static HTML content (serverless-friendly)"""
        stats = results['statistics']
        
        # Load Gene Annotation data
        gene_annotation_rows = self._load_gene_annotation_data(results)
        
        static_html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Static Report</title>
    <style>
/* Inline CSS for serverless compatibility */
* {{
    box-sizing: border-box;
}}

body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    line-height: 1.5;
    margin: 0;
    padding: 2rem;
    background-color: #f8f9fa;
    color: #212529;
}}

.container {{
    max-width: 1200px;
    margin: 0 auto;
    background-color: #fff;
    border-radius: 0.375rem;
    box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
    overflow: hidden;
}}

.header {{
    background: linear-gradient(135deg, #5cb85c, #4cae4c);
    color: #fff;
    padding: 2rem;
    text-align: center;
}}

.header h1 {{
    font-size: 2.5rem;
    margin: 0 0 0.5rem 0;
    font-weight: 700;
}}

.header .subtitle {{
    font-size: 1.1rem;
    opacity: 0.9;
    margin: 0;
}}

.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    padding: 2rem;
    background-color: #f8f9fa;
}}

.stat-card {{
    background: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    padding: 1.25rem;
    text-align: center;
    box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
}}

.stat-value {{
    font-size: 2rem;
    font-weight: 700;
    color: #28a745;
    line-height: 1.2;
}}

.stat-label {{
    color: #6c757d;
    font-size: 0.875rem;
    margin-top: 0.5rem;
    font-weight: 500;
}}

.section {{
    margin: 0;
    border-top: 1px solid #dee2e6;
}}

.section-header {{
    background-color: #f8f9fa;
    border-bottom: 1px solid #dee2e6;
    padding: 1.5rem 2rem;
}}

.section-header h2 {{
    font-size: 1.5rem;
    font-weight: 600;
    color: #212529;
    margin: 0;
}}

.section-body {{
    padding: 2rem;
}}

.table-container {{
    overflow-x: auto;
    margin-top: 1.5rem;
}}

.data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
    background-color: #fff;
}}

.data-table th,
.data-table td {{
    padding: 0.75rem;
    text-align: left;
    border-bottom: 1px solid #dee2e6;
    vertical-align: top;
}}

.data-table th {{
    background-color: #f8f9fa;
    font-weight: 600;
    color: #495057;
    border-top: 1px solid #dee2e6;
    position: sticky;
    top: 0;
    z-index: 10;
}}

.data-table tbody tr:nth-child(even) {{
    background-color: #f8f9fa;
}}

.data-table tbody tr:hover {{
    background-color: #e9ecef;
}}

.classification-multi {{
    background-color: #d1ecf1;
    color: #0c5460;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}}

.classification-single {{
    background-color: #fff3cd;
    color: #856404;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}}

.guide-count {{
    font-weight: 600;
    color: #28a745;
}}

.guide-count-zero {{
    color: #dc3545;
}}

.description-cell {{
    max-width: 300px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}}

.query-id {{
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.8rem;
    background-color: #f8f9fa;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
}}

.uniprot-id {{
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.8rem;
    color: #0066cc;
}}

.info-box {{
    background-color: #d1ecf1;
    border: 1px solid #bee5eb;
    border-radius: 0.375rem;
    padding: 1rem;
    margin-bottom: 1.5rem;
    color: #0c5460;
}}

.footer {{
    background-color: #6c757d;
    color: #fff;
    text-align: center;
    padding: 1rem;
    font-size: 0.875rem;
}}

@media (max-width: 768px) {{
    body {{
        padding: 1rem;
    }}
    
    .stats-grid {{
        grid-template-columns: 1fr;
        padding: 1rem;
    }}
    
    .section-header,
    .section-body {{
        padding: 1rem;
    }}
    
    .header {{
        padding: 1.5rem 1rem;
    }}
    
    .header h1 {{
        font-size: 2rem;
    }}
}}

@media print {{
    body {{
        background-color: #fff;
    }}
    
    .container {{
        box-shadow: none;
    }}
}}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>tWhEEL</h1>
            <p class="subtitle">transcript-aware Whole-genomic Evaluation for Editing Locus pipeline</p>
            <p>Static Report - Generated on {results['timestamp'][:19].replace('T', ' ')}</p>
        </div>

        <!-- Statistics -->
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{stats.get('gene_prediction', {}).get('protein_count', 0):,}</div>
                <div class="stat-label">Predicted Proteins</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('total_genes', 0):,}</div>
                <div class="stat-label">Total Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('multi_genes', 0):,}</div>
                <div class="stat-label">Multi Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('single_genes', 0):,}</div>
                <div class="stat-label">Single Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('go_analysis', {}).get('total_gos', 0):,}</div>
                <div class="stat-label">Total GO Terms</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%</div>
                <div class="stat-label">Fully Targetable GO</div>
            </div>
        </div>

        <!-- Gene Annotation Section -->
        <div class="section">
            <div class="section-header">
                <h2>Gene Annotation Results</h2>
            </div>
            <div class="section-body">
                <div class="info-box">
                    <strong>Static Report:</strong> This is a serverless-friendly version that displays gene annotation results in a static table format.
                    Query ID-based display with guide RNA target counts is shown below.
                </div>
                
                <p><strong>Total Entries:</strong> {len(gene_annotation_rows):,} genes</p>
                
                <div class="table-container">
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>Query ID</th>
                                <th>UniProt Name</th>
                                <th>Gene Name</th>
                                <th>Description</th>
                                <th>tWhEEL Classification<br>Gene Group</th>
                                <th>Number of Targets</th>
                            </tr>
                        </thead>
                        <tbody>
                            {self._generate_gene_annotation_table_rows(gene_annotation_rows)}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Configuration Section -->
        <div class="section">
            <div class="section-header">
                <h2>Analysis Configuration</h2>
            </div>
            <div class="section-body">
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">Input Files</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>Genome Assembly:</strong><br>{results['config']['input'].get('genome_assembly', 'N/A')}</li>
                            <li><strong>Transcriptome:</strong><br>{results['config']['input'].get('transcriptome', 'N/A')}</li>
                            <li><strong>UniProt Database:</strong><br>{results['config']['input'].get('uniprot_db', 'N/A')}</li>
                        </ul>
                    </div>
                    
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">CRISPR Settings</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>PAM Pattern:</strong> {results['config']['crispr'].get('pam_pattern', 'N/A')}</li>
                            <li><strong>Max Mismatch:</strong> {results['config']['crispr'].get('max_mismatch', 'N/A')}</li>
                            <li><strong>Threads:</strong> {results['config']['crispr'].get('threads', 'N/A')}</li>
                        </ul>
                    </div>
                    
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">BUSCO Settings</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>Lineage:</strong> {results['config']['busco'].get('lineage', 'N/A')}</li>
                            <li><strong>CPU Count:</strong> {results['config']['busco'].get('cpu', 'N/A')}</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>

        <!-- Footer -->
        <div class="footer">
            <p>tWhEEL Static Report | Generated from intermediate data files</p>
        </div>
    </div>
</body>
</html>
"""
        
        return static_html_content
    
    def _load_gene_annotation_data(self, results):
        """Load gene annotation data from TSV file"""
        annotation_path = results['data'].get('final_annotation_guides')
        if not annotation_path or not os.path.exists(annotation_path):
            return []
        
        gene_data = []
        
        try:
            with open(annotation_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)
                
                for row in reader:
                    if len(row) >= 11:
                        query_id = row[0]
                        subject_id = row[1]
                        gene_name = row[7] if len(row) > 7 else ''
                        classification = row[8] if len(row) > 8 else ''
                        unique_id = row[9] if len(row) > 9 else ''
                        description = row[6] if len(row) > 6 else ''
                        guide_targets = row[-1] if len(row) > 10 else ''
                        
                        # Count guide targets
                        target_count = 0
                        if guide_targets and guide_targets.strip():
                            targets = guide_targets.split(',')
                            target_count = len([t for t in targets if t.strip()])
                        
                        gene_data.append({
                            'query_id': query_id,
                            'subject_id': subject_id,
                            'gene_name': gene_name,
                            'classification': classification,
                            'unique_id': unique_id,
                            'description': description,
                            'target_count': target_count
                        })
                        
        except Exception as e:
            self.logger.error(f"Error loading gene annotation data: {e}")
            return []
        
        return gene_data
    
    def _load_gene_annotation_data_for_js(self, results):
        """Load gene annotation data for JavaScript embedding"""
        annotation_path = results['data'].get('final_annotation_guides')
        if not annotation_path or not os.path.exists(annotation_path):
            return []
        
        gene_data = []
        
        try:
            with open(annotation_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)
                
                for i, row in enumerate(reader):
                        
                    if len(row) >= 11:
                        query_id = row[0]
                        subject_id = row[1]
                        gene_name = row[7] if len(row) > 7 else ''
                        classification = row[8] if len(row) > 8 else ''
                        unique_id = row[9] if len(row) > 9 else ''
                        description = row[6] if len(row) > 6 else ''
                        guide_targets = row[-1] if len(row) > 10 else ''
                        
                        # Count guide targets
                        target_count = 0
                        if guide_targets and guide_targets.strip():
                            targets = guide_targets.split(',')
                            target_count = len([t for t in targets if t.strip()])
                        
                        # Truncate long descriptions for JSON embedding
                        description_short = description[:200] if len(description) > 200 else description
                        
                        gene_data.append({
                            'query_id': query_id,
                            'subject_id': subject_id,
                            'gene_name': gene_name,
                            'classification': classification,
                            'unique_id': unique_id,
                            'description': description_short,
                            'target_count': target_count
                        })
                        
        except Exception as e:
            self.logger.error(f"Error loading gene annotation data for JS: {e}")
            return []
        
        return gene_data
    
    def _load_guide_rna_data_for_js(self, results):
        """Load guide RNA data for JavaScript embedding from annotation_with_guides_final.tsv"""
        annotation_path = results['data'].get('final_annotation_guides')
        if not annotation_path or not os.path.exists(annotation_path):
            return []
        
        guide_data = []
        
        try:
            with open(annotation_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)
                
                for i, row in enumerate(reader):
                        
                    if len(row) >= 11 and row[-1]:  # guide_targets column exists
                        query_id = row[0]
                        guide_targets = row[-1]
                        
                        # Parse guide targets - each guide becomes a separate row
                        if guide_targets and guide_targets.strip():
                            targets = guide_targets.split(',')
                            for target in targets:
                                if target.strip():
                                    guide_parts = target.strip().split('/')
                                    if len(guide_parts) >= 8:  # Ensure we have all required fields
                                        guide_data.append({
                                            'query_id': query_id,
                                            'reference': guide_parts[0],
                                            'start': int(guide_parts[1]) if guide_parts[1].isdigit() else 0,
                                            'end': int(guide_parts[2]) if guide_parts[2].isdigit() else 0,
                                            'target_sequence': guide_parts[3],
                                            'mm0': int(guide_parts[4]) if guide_parts[4].isdigit() else 0,
                                            'mm1_2': int(guide_parts[5]) if guide_parts[5].isdigit() else 0,
                                            'off_target_score': float(guide_parts[6]) if guide_parts[6].replace('.', '').isdigit() else 0.0,
                                            'strand': guide_parts[7]
                                        })
                                        
                                        # Stop if we reach 1000 guide RNAs
                                        if len(guide_data) >= 1000:
                                            break
                                
                                if len(guide_data) >= 1000:
                                    break
                        
        except Exception as e:
            self.logger.error(f"Error loading guide RNA data for JS: {e}")
            return []
        
        self.logger.info(f"Loaded {len(guide_data)} guide RNA entries for JavaScript embedding")
        return guide_data
    
    def _load_go_data_for_js(self, results):
        """Load GO analysis data for JavaScript embedding from go_kfi_table.csv with API data"""
        kfi_table_path = results['data'].get('go_kfi_table')
        if not kfi_table_path or not os.path.exists(kfi_table_path):
            return []
        
        go_data = []
        go_ids = []
        
        try:
            # First, load basic GO data from go_kfi_table.csv
            with open(kfi_table_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)
                
                for i, row in enumerate(reader):
                        
                    if len(row) >= 6:
                        go_id = row[0]
                        total_genes = int(row[1]) if row[1] else 0
                        genes_with_guides = int(row[2]) if row[2] else 0
                        kfi_score = float(row[4]) if row[4] else 0.0
                        targetable_type = row[5] if len(row) > 5 else 'Not'
                        
                        go_data.append({
                            'go_id': go_id,
                            'go_name': '',  # Will be filled by API
                            'go_definition': '',  # Will be filled by API
                            'total_genes': total_genes,
                            'genes_with_guides': genes_with_guides,
                            'kfi_score': kfi_score,
                            'targetable_type': targetable_type,
                            'annotated_genes': ''  # Will be filled from reverse lookup
                        })
                        go_ids.append(go_id)
            
            # Load GO terms and definitions from database instead of API to avoid duplicate fetching
            self.logger.info(f"Loading GO terms information from database for {len(go_ids)} GO IDs for JavaScript embedding")
            try:
                go_terms_info = self._fetch_go_terms_from_database(go_ids)
            except Exception as e:
                self.logger.warning(f"Failed to load GO terms from database for JS: {e}")
                go_terms_info = {}
            
            # Parse GO-annotated genes from reverse lookup
            self.logger.info("Parsing GO-annotated genes for JavaScript embedding")
            try:
                go_genes_map = self._parse_go_annotated_genes(results)
            except Exception as e:
                self.logger.warning(f"Failed to parse GO-annotated genes for JS: {e}")
                go_genes_map = {}
            
            # Merge API data and annotated genes data into go_data
            for go_entry in go_data:
                go_id = go_entry['go_id']
                
                # Add GO term information from API
                if go_id in go_terms_info:
                    go_entry['go_name'] = go_terms_info[go_id].get('name', '')
                    go_entry['go_definition'] = go_terms_info[go_id].get('definition', '')
                
                # Add annotated genes information
                if go_id in go_genes_map:
                    annotated_genes = go_genes_map[go_id]
                    go_entry['annotated_genes'] = self._format_annotated_genes_list(annotated_genes)
                    # Always include the full genes list for expandable functionality
                    go_entry['annotated_genes_full'] = ", ".join(annotated_genes) if annotated_genes else ""
                else:
                    # If no genes found, set both to empty
                    go_entry['annotated_genes'] = ""
                    go_entry['annotated_genes_full'] = ""
                        
        except Exception as e:
            self.logger.error(f"Error loading GO data for JS: {e}")
            return []
        
        self.logger.info(f"Loaded {len(go_data)} GO analysis entries with enriched data for JavaScript embedding")
        return go_data
    
    def _load_go_parent_counts_data_for_js(self, results):
        """Load GO parent counts data for JavaScript embedding from go_parent_counts.csv"""
        go_parent_counts_path = os.path.join(self.intermediate_dir, 'go_parent_counts.csv')
        if not os.path.exists(go_parent_counts_path):
            self.logger.warning(f"GO parent counts file not found: {go_parent_counts_path}")
            return []
        
        go_parent_data = []
        
        try:
            with open(go_parent_counts_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)  # Skip header: Parent GO,Description,count
                
                for i, row in enumerate(reader):
                    if len(row) >= 3:
                        parent_go = row[0]
                        description = row[1]
                        count = int(row[2]) if row[2].isdigit() else 0
                        
                        go_parent_data.append({
                            'parent_go': parent_go,
                            'description': description,
                            'count': count
                        })
                        
        except Exception as e:
            self.logger.error(f"Error loading GO parent counts data for JS: {e}")
            return []
        
        self.logger.info(f"Loaded {len(go_parent_data)} GO parent counts entries for JavaScript embedding")
        return go_parent_data
    
    def _generate_gene_annotation_table_rows(self, gene_data):
        """Generate HTML table rows for gene annotation data"""
        if not gene_data:
            return '<tr><td colspan="6">No gene annotation data available</td></tr>'
        
        rows = []
        for i, gene in enumerate(gene_data):
            # Limit to first 1000 rows to prevent excessive file size
            if i >= 1000:
                rows.append(f'''
                    <tr>
                        <td colspan="6" style="text-align: center; padding: 1rem; font-style: italic; color: #6c757d;">
                            ... and {len(gene_data) - 1000:,} more entries (showing first 1,000 for performance)
                        </td>
                    </tr>
                ''')
                break
            
            classification_class = 'classification-multi' if gene['classification'] == 'Multi' else 'classification-single'
            guide_count_class = 'guide-count' if gene['target_count'] > 0 else 'guide-count-zero'
            
            # Truncate long descriptions
            description = gene['description'][:100] + '...' if len(gene['description']) > 100 else gene['description']
            
            row = f'''
                <tr>
                    <td><span class="query-id">{gene['query_id']}</span></td>
                    <td><span class="uniprot-id">{gene['subject_id']}</span></td>
                    <td><strong>{gene['gene_name']}</strong></td>
                    <td class="description-cell" title="{gene['description']}">{description}</td>
                    <td><span class="{classification_class}">{gene['classification']}</span><br><small>{gene['unique_id']}</small></td>
                    <td><span class="{guide_count_class}">{gene['target_count']}</span></td>
                </tr>
            '''
            rows.append(row)
        
        return ''.join(rows)
    
    def run(self):
        return html_content
    
    def _generate_static_html_content(self, results):
        """Generate static HTML content (serverless-friendly)"""
        stats = results['statistics']
        
        # Load Gene Annotation data
        gene_annotation_rows = self._load_gene_annotation_data(results)
        
        static_html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Static Report</title>
    <style>
/* Inline CSS for serverless compatibility */
* {{
    box-sizing: border-box;
}}

body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    line-height: 1.5;
    margin: 0;
    padding: 2rem;
    background-color: #f8f9fa;
    color: #212529;
}}

.container {{
    max-width: 1200px;
    margin: 0 auto;
    background-color: #fff;
    border-radius: 0.375rem;
    box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
    overflow: hidden;
}}

.header {{
    background: linear-gradient(135deg, #5cb85c, #4cae4c);
    color: #fff;
    padding: 2rem;
    text-align: center;
}}

.header h1 {{
    font-size: 2.5rem;
    margin: 0 0 0.5rem 0;
    font-weight: 700;
}}

.header .subtitle {{
    font-size: 1.1rem;
    opacity: 0.9;
    margin: 0;
}}

.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    padding: 2rem;
    background-color: #f8f9fa;
}}

.stat-card {{
    background: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    padding: 1.25rem;
    text-align: center;
    box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
}}

.stat-value {{
    font-size: 2rem;
    font-weight: 700;
    color: #28a745;
    line-height: 1.2;
}}

.stat-label {{
    color: #6c757d;
    font-size: 0.875rem;
    margin-top: 0.5rem;
    font-weight: 500;
}}

.section {{
    margin: 0;
    border-top: 1px solid #dee2e6;
}}

.section-header {{
    background-color: #f8f9fa;
    border-bottom: 1px solid #dee2e6;
    padding: 1.5rem 2rem;
}}

.section-header h2 {{
    font-size: 1.5rem;
    font-weight: 600;
    color: #212529;
    margin: 0;
}}

.section-body {{
    padding: 2rem;
}}

.table-container {{
    overflow-x: auto;
    margin-top: 1.5rem;
}}

.data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.875rem;
    background-color: #fff;
}}

.data-table th,
.data-table td {{
    padding: 0.75rem;
    text-align: left;
    border-bottom: 1px solid #dee2e6;
    vertical-align: top;
}}

.data-table th {{
    background-color: #f8f9fa;
    font-weight: 600;
    color: #495057;
    border-top: 1px solid #dee2e6;
    position: sticky;
    top: 0;
    z-index: 10;
}}

.data-table tbody tr:nth-child(even) {{
    background-color: #f8f9fa;
}}

.data-table tbody tr:hover {{
    background-color: #e9ecef;
}}

.classification-multi {{
    background-color: #d1ecf1;
    color: #0c5460;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}}

.classification-single {{
    background-color: #fff3cd;
    color: #856404;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    font-size: 0.8rem;
    font-weight: 500;
}}

.guide-count {{
    font-weight: 600;
    color: #28a745;
}}

.guide-count-zero {{
    color: #dc3545;
}}

.description-cell {{
    max-width: 300px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}}

.query-id {{
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.8rem;
    background-color: #f8f9fa;
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
}}

.uniprot-id {{
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
    font-size: 0.8rem;
    color: #0066cc;
}}

.info-box {{
    background-color: #d1ecf1;
    border: 1px solid #bee5eb;
    border-radius: 0.375rem;
    padding: 1rem;
    margin-bottom: 1.5rem;
    color: #0c5460;
}}

.footer {{
    background-color: #6c757d;
    color: #fff;
    text-align: center;
    padding: 1rem;
    font-size: 0.875rem;
}}

@media (max-width: 768px) {{
    body {{
        padding: 1rem;
    }}
    
    .stats-grid {{
        grid-template-columns: 1fr;
        padding: 1rem;
    }}
    
    .section-header,
    .section-body {{
        padding: 1rem;
    }}
    
    .header {{
        padding: 1.5rem 1rem;
    }}
    
    .header h1 {{
        font-size: 2rem;
    }}
}}

@media print {{
    body {{
        background-color: #fff;
    }}
    
    .container {{
        box-shadow: none;
    }}
}}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>tWhEEL</h1>
            <p class="subtitle">Transcriptome Workflow for Hierarchical Exploration and Expression Levels</p>
            <p>Static Report - Generated on {results['timestamp'][:19].replace('T', ' ')}</p>
        </div>

        <!-- Statistics -->
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{stats.get('gene_prediction', {}).get('protein_count', 0):,}</div>
                <div class="stat-label">Predicted Proteins</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('total_genes', 0):,}</div>
                <div class="stat-label">Total Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('multi_genes', 0):,}</div>
                <div class="stat-label">Multi Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('annotation', {}).get('single_genes', 0):,}</div>
                <div class="stat-label">Single Genes</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('go_analysis', {}).get('total_gos', 0):,}</div>
                <div class="stat-label">Total GO Terms</div>
            </div>
            
            <div class="stat-card">
                <div class="stat-value">{stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%</div>
                <div class="stat-label">Fully Targetable GO</div>
            </div>
        </div>

        <!-- Gene Annotation Section -->
        <div class="section">
            <div class="section-header">
                <h2>Gene Annotation Results</h2>
            </div>
            <div class="section-body">
                <div class="info-box">
                    <strong>Static Report:</strong> This is a serverless-friendly version that displays gene annotation results in a static table format.
                    Query ID-based display with guide RNA target counts is shown below.
                </div>
                
                <p><strong>Total Entries:</strong> {len(gene_annotation_rows):,} genes</p>
                
                <div class="table-container">
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>Query ID</th>
                                <th>UniProt Name</th>
                                <th>Gene Name</th>
                                <th>Description</th>
                                <th>tWhEEL Classification<br>Gene Group</th>
                                <th>Number of Targets</th>
                            </tr>
                        </thead>
                        <tbody>
                            {self._generate_gene_annotation_table_rows(gene_annotation_rows)}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Configuration Section -->
        <div class="section">
            <div class="section-header">
                <h2>Analysis Configuration</h2>
            </div>
            <div class="section-body">
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">Input Files</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>Genome Assembly:</strong><br>{results['config']['input'].get('genome_assembly', 'N/A')}</li>
                            <li><strong>Transcriptome:</strong><br>{results['config']['input'].get('transcriptome', 'N/A')}</li>
                            <li><strong>UniProt Database:</strong><br>{results['config']['input'].get('uniprot_db', 'N/A')}</li>
                        </ul>
                    </div>
                    
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">CRISPR Settings</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>PAM Pattern:</strong> {results['config']['crispr'].get('pam_pattern', 'N/A')}</li>
                            <li><strong>Max Mismatch:</strong> {results['config']['crispr'].get('max_mismatch', 'N/A')}</li>
                            <li><strong>Threads:</strong> {results['config']['crispr'].get('threads', 'N/A')}</li>
                        </ul>
                    </div>
                    
                    <div class="stat-card">
                        <h3 style="font-size: 1.1rem; margin: 0 0 1rem 0;">BUSCO Settings</h3>
                        <ul style="text-align: left; font-size: 0.85rem; margin: 0; padding-left: 1rem;">
                            <li><strong>Lineage:</strong> {results['config']['busco'].get('lineage', 'N/A')}</li>
                            <li><strong>CPU Count:</strong> {results['config']['busco'].get('cpu', 'N/A')}</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>

        <!-- Footer -->
        <div class="footer">
            <p>tWhEEL Static Report | Generated from intermediate data files</p>
        </div>
    </div>
</body>
</html>
"""
        
        return static_html_content
    
    def run(self):
        """Execute overall report generation"""
        try:
            # Collect all results
            results = self.collect_all_results()
            
            # Generate text summary
            if self.config['report']['generate_summary']:
                summary_path = self.generate_text_summary(results)
                results['summary_path'] = summary_path
            
            # Prepare web data
            db_path, json_path = self.prepare_web_data(results)
            results['db_path'] = db_path
            results['json_path'] = json_path
            
            # Generate HTML report
            if self.config['report']['generate_html']:
                html_path = self.generate_html_report(results)
                results['html_path'] = html_path
            
            # Save results to file
            results_file = os.path.join(self.intermediate_dir, 'report_generation_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== Report Generation Results ===\n\n")
                f.write(f"Execution date: {results['timestamp']}\n\n")
                
                if 'summary_path' in results:
                    f.write(f"Text summary: {results['summary_path']}\n")
                
                if 'html_path' in results:
                    f.write(f"HTML report: {results['html_path']}\n")
                
                f.write(f"SQLite database: {results['db_path']}\n")
                f.write(f"JSON summary: {results['json_path']}\n\n")
                
                f.write("=== Statistics ===\n")
                for category, stats in results['statistics'].items():
                    f.write(f"{category}:\n")
                    for key, value in stats.items():
                        f.write(f"  {key}: {value}\n")
                    f.write("\n")
            
            self.logger.info(f"Report generation results saved: {results_file}")
            self.logger.info("Report generation process completed")
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error occurred during report generation: {e}")
            raise
