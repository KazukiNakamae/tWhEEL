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
from datetime import datetime
from pathlib import Path
import shutil


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
        
        # Collect statistics
        results['statistics'] = self._collect_statistics()
        
        self.logger.info("Results collection completed")
        return results
    
    def _collect_statistics(self):
        """Collect statistics information"""
        stats = {
            'gene_prediction': {},
            'annotation': {},
            'guide_design': {},
            'go_analysis': {}
        }
        
        # Gene prediction statistics
        try:
            protein_path = os.path.join(self.intermediate_dir, 'query.pep')
            if os.path.exists(protein_path):
                stats['gene_prediction']['protein_count'] = self._count_fasta_sequences(protein_path)
        except Exception as e:
            self.logger.warning(f"Gene prediction statistics collection error: {e}")
        
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
                
        except Exception as e:
            self.logger.warning(f"GO-KFI statistics parsing error: {e}")
        
        return stats
    
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
                f.write(f"Predicted proteins: {stats['gene_prediction'].get('protein_count', 'N/A')}\n")
            
            # Annotation statistics
            if 'annotation' in stats:
                f.write(f"Multi genes: {stats['annotation'].get('multi_genes', 'N/A')}\n")
                f.write(f"Single genes: {stats['annotation'].get('single_genes', 'N/A')}\n")
                f.write(f"Total genes: {stats['annotation'].get('total_genes', 'N/A')}\n")
            
            # Guide design statistics
            if 'guide_design' in stats:
                f.write(f"Multi gene targets: {stats['guide_design'].get('multi_targets', 'N/A')}\n")
                f.write(f"Single gene targets: {stats['guide_design'].get('single_targets', 'N/A')}\n")
            
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
    
    def prepare_web_data(self, results):
        """Prepare data for web application in JSON and SQLite formats"""
        self.logger.info("Preparing web data...")
        
        # Create SQLite database
        db_path = os.path.join(self.web_data_dir, 'twheel_data.db')
        
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
                go_namespace TEXT,
                kfi_score REAL,
                targetable_type TEXT,
                gene_count INTEGER
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
        """Insert guide RNA data"""
        
        final_table_path = results['data'].get('final_annotation_guides')
        if not final_table_path or not os.path.exists(final_table_path):
            return
        
        with open(final_table_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            for row in reader:
                if len(row) >= 11 and row[-1]:  # guide_targets column exists
                    unique_id = row[9]
                    guide_targets = row[-1]
                    
                    # Parse guide targets
                    if guide_targets:
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
                    
                    # Update gene guide count
                    cursor.execute("""
                        UPDATE genes SET guide_count = (
                            SELECT COUNT(*) FROM guide_rnas WHERE gene_unique_id = ?
                        ) WHERE unique_id = ?
                    """, (unique_id, unique_id))
    
    def _insert_go_data(self, cursor, results):
        """Insert GO term data"""
        
        kfi_table_path = results['data'].get('go_kfi_table')
        if not kfi_table_path or not os.path.exists(kfi_table_path):
            return
        
        with open(kfi_table_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            
            for row in reader:
                if len(row) >= 5:
                    go_id = row[0]
                    kfi_score = float(row[4]) if row[4] else 0.0
                    
                    # Determine targetable_type
                    if kfi_score == 1.0:
                        targetable_type = 'Fully'
                    elif kfi_score > 0.0:
                        targetable_type = 'Partially'
                    else:
                        targetable_type = 'Not'
                    
                    cursor.execute("""
                        INSERT INTO go_terms (
                            go_id, kfi_score, targetable_type, gene_count
                        ) VALUES (?, ?, ?, ?)
                    """, (go_id, kfi_score, targetable_type, 0))
    
    def generate_html_report(self, results):
        """Generate HTML report"""
        self.logger.info("Generating HTML report...")
        
        # Create template files
        self._create_html_template()
        self._create_css_files()
        self._create_js_files()
        
        # Generate main HTML file
        html_path = os.path.join(self.output_dir, 'report.html')
        
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(self._generate_html_content(results))
        
        self.logger.info(f"HTML report generated: {html_path}")
        return html_path
    
    def _create_html_template(self):
        """Create HTML template"""
        template_dir = os.path.join('src', 'web', 'templates')
        Path(template_dir).mkdir(parents=True, exist_ok=True)
        
        # Basic template files will be created later
        pass
    
    def _create_css_files(self):
        """Create CSS files"""
        css_dir = os.path.join('src', 'web', 'static', 'css')
        
        # Main CSS file
        main_css_path = os.path.join(css_dir, 'main.css')
        with open(main_css_path, 'w', encoding='utf-8') as f:
            f.write("""
/* tWhEEL Report CSS */
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    line-height: 1.6;
    margin: 0;
    padding: 0;
    background-color: #f5f5f5;
}

.container {
    max-width: 1200px;
    margin: 0 auto;
    padding: 20px;
}

.header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 2rem 0;
    text-align: center;
}

.nav-tabs {
    display: flex;
    background: white;
    border-bottom: 2px solid #e0e0e0;
    margin-bottom: 2rem;
}

.nav-tab {
    padding: 1rem 2rem;
    cursor: pointer;
    border: none;
    background: none;
    font-size: 1rem;
    transition: all 0.3s ease;
}

.nav-tab:hover {
    background: #f0f0f0;
}

.nav-tab.active {
    background: #667eea;
    color: white;
}

.tab-content {
    display: none;
    background: white;
    padding: 2rem;
    border-radius: 8px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
}

.tab-content.active {
    display: block;
}

.stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 1.5rem;
    margin-bottom: 2rem;
}

.stat-card {
    background: white;
    padding: 1.5rem;
    border-radius: 8px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    text-align: center;
}

.stat-value {
    font-size: 2rem;
    font-weight: bold;
    color: #667eea;
}

.stat-label {
    color: #666;
    margin-top: 0.5rem;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
    margin-top: 1rem;
}

.data-table th,
.data-table td {
    padding: 12px;
    text-align: left;
    border-bottom: 1px solid #e0e0e0;
}

.data-table th {
    background: #f8f9fa;
    font-weight: 600;
}

.search-box {
    width: 100%;
    padding: 12px;
    margin-bottom: 1rem;
    border: 1px solid #ddd;
    border-radius: 4px;
    font-size: 1rem;
}

.filter-group {
    display: flex;
    gap: 1rem;
    margin-bottom: 1rem;
    flex-wrap: wrap;
}

.filter-select {
    padding: 8px 12px;
    border: 1px solid #ddd;
    border-radius: 4px;
}

.chart-container {
    height: 400px;
    margin: 2rem 0;
}

.guide-sequence {
    font-family: 'Courier New', monospace;
    background: #f8f9fa;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 0.9rem;
}

.targetable-fully {
    background: #d4edda;
    color: #155724;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 0.8rem;
}

.targetable-partially {
    background: #fff3cd;
    color: #856404;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 0.8rem;
}

.targetable-not {
    background: #f8d7da;
    color: #721c24;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 0.8rem;
}

.loading {
    text-align: center;
    padding: 3rem;
    color: #666;
}

.error {
    background: #f8d7da;
    color: #721c24;
    padding: 1rem;
    border-radius: 4px;
    margin: 1rem 0;
}

@media (max-width: 768px) {
    .nav-tabs {
        flex-direction: column;
    }
    
    .filter-group {
        flex-direction: column;
    }
    
    .stats-grid {
        grid-template-columns: 1fr;
    }
}
""")
    
    def _create_js_files(self):
        """Create JavaScript files"""
        js_dir = os.path.join('src', 'web', 'static', 'js')
        
        # Main JS file
        main_js_path = os.path.join(js_dir, 'main.js')
        with open(main_js_path, 'w', encoding='utf-8') as f:
            f.write("""
// tWhEEL Report JavaScript

class tWhEELReport {
    constructor() {
        this.currentTab = 'dashboard';
        this.data = null;
        this.init();
    }
    
    async init() {
        await this.loadData();
        this.setupEventListeners();
        this.showTab('dashboard');
    }
    
    async loadData() {
        try {
            const response = await fetch('web_data/summary.json');
            this.data = await response.json();
            this.populateDashboard();
        } catch (error) {
            console.error('Failed to load data:', error);
            this.showError('Failed to load data');
        }
    }
    
    setupEventListeners() {
        // Tab switching
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.addEventListener('click', (e) => {
                const tabId = e.target.dataset.tab;
                this.showTab(tabId);
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
    }
    
    showTab(tabId) {
        // Tab switching
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.classList.remove('active');
        });
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.remove('active');
        });
        
        document.querySelector(`[data-tab="${tabId}"]`).classList.add('active');
        document.querySelector(`#${tabId}`).classList.add('active');
        
        this.currentTab = tabId;
        
        // Tab-specific initialization
        switch(tabId) {
            case 'genes':
                this.loadGeneData();
                break;
            case 'guides':
                this.loadGuideData();
                break;
            case 'go-analysis':
                this.loadGOData();
                break;
        }
    }
    
    populateDashboard() {
        if (!this.data) return;
        
        const stats = this.data.statistics;
        
        // Update statistics cards
        this.updateStatCard('protein-count', stats.gene_prediction?.protein_count || 0);
        this.updateStatCard('total-genes', stats.annotation?.total_genes || 0);
        this.updateStatCard('multi-genes', stats.annotation?.multi_genes || 0);
        this.updateStatCard('single-genes', stats.annotation?.single_genes || 0);
        this.updateStatCard('total-gos', stats.go_analysis?.total_gos || 0);
        this.updateStatCard('fully-targetable', stats.go_analysis?.fully_targetable_percent || 0);
    }
    
    updateStatCard(id, value) {
        const element = document.getElementById(id);
        if (element) {
            element.textContent = typeof value === 'number' ? value.toLocaleString() : value;
        }
    }
    
    showError(message) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'error';
        errorDiv.textContent = message;
        document.body.insertBefore(errorDiv, document.body.firstChild);
    }
    
    async loadGeneData() {
        // Gene data loading (to be implemented)
        console.log('Loading gene data...');
    }
    
    async loadGuideData() {
        // Guide RNA data loading (to be implemented)
        console.log('Loading guide data...');
    }
    
    async loadGOData() {
        // GO data loading (to be implemented)
        console.log('Loading GO data...');
    }
    
    handleSearch(query, target) {
        // Search functionality (to be implemented)
        console.log('Search:', query, 'Target:', target);
    }
    
    handleFilter(value, filter) {
        // Filter functionality (to be implemented)
        console.log('Filter:', filter, 'Value:', value);
    }
}

// Application initialization
document.addEventListener('DOMContentLoaded', () => {
    new tWhEELReport();
});
""")
    
    def _generate_html_content(self, results):
        """Generate HTML content"""
        stats = results['statistics']
        
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Transcriptome Functional Analysis Report</title>
    <link rel="stylesheet" href="src/web/static/css/main.css">
</head>
<body>
    <div class="header">
        <div class="container">
            <h1>tWhEEL</h1>
            <p>Transcriptome Workflow for Hierarchical Exploration and Expression Levels</p>
            <p>Execution Date: {results['timestamp']}</p>
        </div>
    </div>
    
    <div class="container">
        <div class="nav-tabs">
            <button class="nav-tab active" data-tab="dashboard">Dashboard</button>
            <button class="nav-tab" data-tab="genes">Gene Browser</button>
            <button class="nav-tab" data-tab="guides">Guide RNA</button>
            <button class="nav-tab" data-tab="go-analysis">GO Analysis</button>
            <button class="nav-tab" data-tab="settings">Settings</button>
        </div>
        
        <!-- Dashboard -->
        <div id="dashboard" class="tab-content active">
            <h2>Analysis Results Summary</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value" id="protein-count">
                        {stats.get('gene_prediction', {}).get('protein_count', 0):,}
                    </div>
                    <div class="stat-label">Predicted Proteins</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="total-genes">
                        {stats.get('annotation', {}).get('total_genes', 0):,}
                    </div>
                    <div class="stat-label">Total Genes</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="multi-genes">
                        {stats.get('annotation', {}).get('multi_genes', 0):,}
                    </div>
                    <div class="stat-label">Multi Genes</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="single-genes">
                        {stats.get('annotation', {}).get('single_genes', 0):,}
                    </div>
                    <div class="stat-label">Single Genes</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="total-gos">
                        {stats.get('go_analysis', {}).get('total_gos', 0):,}
                    </div>
                    <div class="stat-label">Total GO Terms</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="fully-targetable">
                        {stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%
                    </div>
                    <div class="stat-label">Fully Targetable GO</div>
                </div>
            </div>
            
            <div class="chart-container">
                <h3>GO Knock-out Feasibility Distribution</h3>
                <div id="go-kfi-chart">
                    <p class="loading">Loading chart...</p>
                </div>
            </div>
        </div>
        
        <!-- Gene Browser -->
        <div id="genes" class="tab-content">
            <h2>Gene Browser</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="Search by gene name..." data-target="genes">
                <select class="filter-select" data-filter="classification">
                    <option value="">All</option>
                    <option value="Multi">Multi Genes</option>
                    <option value="Single">Single Genes</option>
                </select>
                <select class="filter-select" data-filter="guides">
                    <option value="">All</option>
                    <option value="with-guides">With Guide RNA</option>
                    <option value="without-guides">Without Guide RNA</option>
                </select>
            </div>
            
            <table class="data-table" id="genes-table">
                <thead>
                    <tr>
                        <th>Gene Name</th>
                        <th>Classification</th>
                        <th>UniProt ID</th>
                        <th>Identity</th>
                        <th>E-value</th>
                        <th>Guide Count</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="7" class="loading">Loading data...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- Guide RNA -->
        <div id="guides" class="tab-content">
            <h2>Guide RNA Design Results</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="Search by guide sequence..." data-target="guides">
                <select class="filter-select" data-filter="chromosome">
                    <option value="">All Chromosomes</option>
                </select>
                <select class="filter-select" data-filter="strand">
                    <option value="">All Strands</option>
                    <option value="+">+</option>
                    <option value="-">-</option>
                </select>
            </div>
            
            <table class="data-table" id="guides-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Chromosome</th>
                        <th>Start</th>
                        <th>End</th>
                        <th>Guide Sequence</th>
                        <th>Strand</th>
                        <th>Score</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="7" class="loading">Loading data...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- GO Analysis -->
        <div id="go-analysis" class="tab-content">
            <h2>GO Functional Analysis Results</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="Search by GO term..." data-target="go-terms">
                <select class="filter-select" data-filter="targetable">
                    <option value="">All</option>
                    <option value="Fully">Fully targetable</option>
                    <option value="Partially">Partially targetable</option>
                    <option value="Not">Not targetable</option>
                </select>
            </div>
            
            <table class="data-table" id="go-table">
                <thead>
                    <tr>
                        <th>GO ID</th>
                        <th>GO Name</th>
                        <th>Namespace</th>
                        <th>KFI Score</th>
                        <th>Targetable</th>
                        <th>Gene Count</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="6" class="loading">Loading data...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- Settings -->
        <div id="settings" class="tab-content">
            <h2>Configuration Settings</h2>
            
            <div class="setting-group">
                <h3>Input Files</h3>
                <ul>
                    <li><strong>Genome Assembly:</strong> {results['config']['input'].get('genome_assembly', 'N/A')}</li>
                    <li><strong>Transcriptome:</strong> {results['config']['input'].get('transcriptome', 'N/A')}</li>
                    <li><strong>UniProt Database:</strong> {results['config']['input'].get('uniprot_db', 'N/A')}</li>
                </ul>
            </div>
            
            <div class="setting-group">
                <h3>CRISPR Settings</h3>
                <ul>
                    <li><strong>PAM Pattern:</strong> {results['config']['crispr'].get('pam_pattern', 'N/A')}</li>
                    <li><strong>Max Mismatch:</strong> {results['config']['crispr'].get('max_mismatch', 'N/A')}</li>
                    <li><strong>Threads:</strong> {results['config']['crispr'].get('threads', 'N/A')}</li>
                </ul>
            </div>
            
            <div class="setting-group">
                <h3>BUSCO Settings</h3>
                <ul>
                    <li><strong>Lineage:</strong> {results['config']['busco'].get('lineage', 'N/A')}</li>
                    <li><strong>CPU Count:</strong> {results['config']['busco'].get('cpu', 'N/A')}</li>
                </ul>
            </div>
        </div>
    </div>
    
    <script src="src/web/static/js/main.js"></script>
</body>
</html>
"""
        
        return html_content
    
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
