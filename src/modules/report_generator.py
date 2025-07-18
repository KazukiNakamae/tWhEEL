#!/usr/bin/env python3
"""
レポート生成モジュール
テキストサマリーとHTMLインタラクティブレポートの生成
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
    """レポート生成処理クラス"""
    
    def __init__(self, config):
        """
        初期化
        
        Args:
            config (dict): 設定辞書
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.output_dir = self.config['output']['base_dir']
        self.web_data_dir = os.path.join(self.output_dir, 'web_data')
        
    def collect_all_results(self):
        """すべてのステップの結果を収集"""
        self.logger.info("すべてのステップの結果を収集中...")
        
        results = {
            'timestamp': datetime.now().isoformat(),
            'config': self.config,
            'files': {},
            'statistics': {},
            'data': {}
        }
        
        # 各ステップの結果ファイルを収集
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
        
        # 主要データファイルの収集
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
        
        # 統計情報の収集
        results['statistics'] = self._collect_statistics()
        
        self.logger.info("結果の収集が完了しました")
        return results
    
    def _collect_statistics(self):
        """統計情報の収集"""
        stats = {
            'gene_prediction': {},
            'annotation': {},
            'guide_design': {},
            'go_analysis': {}
        }
        
        # 遺伝子予測統計
        try:
            protein_path = os.path.join(self.intermediate_dir, 'query.pep')
            if os.path.exists(protein_path):
                stats['gene_prediction']['protein_count'] = self._count_fasta_sequences(protein_path)
        except Exception as e:
            self.logger.warning(f"遺伝子予測統計の収集エラー: {e}")
        
        # アノテーション統計
        try:
            annotation_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
            if os.path.exists(annotation_path):
                multi_count, single_count = self._count_gene_classifications(annotation_path)
                stats['annotation']['multi_genes'] = multi_count
                stats['annotation']['single_genes'] = single_count
                stats['annotation']['total_genes'] = multi_count + single_count
        except Exception as e:
            self.logger.warning(f"アノテーション統計の収集エラー: {e}")
        
        # ガイド設計統計
        try:
            multi_targets = os.path.join(self.intermediate_dir, 'multi_genes_targets_specific_clean.bed')
            single_targets = os.path.join(self.intermediate_dir, 'single_genes_targets_specific_clean.bed')
            
            if os.path.exists(multi_targets):
                stats['guide_design']['multi_targets'] = self._count_lines(multi_targets)
            if os.path.exists(single_targets):
                stats['guide_design']['single_targets'] = self._count_lines(single_targets)
        except Exception as e:
            self.logger.warning(f"ガイド設計統計の収集エラー: {e}")
        
        # GO解析統計
        try:
            kfi_stat_path = os.path.join(self.intermediate_dir, 'go_kfi_stat.csv')
            if os.path.exists(kfi_stat_path):
                stats['go_analysis'] = self._parse_go_kfi_stats(kfi_stat_path)
        except Exception as e:
            self.logger.warning(f"GO解析統計の収集エラー: {e}")
        
        return stats
    
    def _count_fasta_sequences(self, fasta_path):
        """FASTAファイルの配列数をカウント"""
        count = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    
    def _count_gene_classifications(self, annotation_path):
        """遺伝子分類のカウント"""
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
        """ファイルの行数をカウント"""
        try:
            with open(filepath, 'r') as f:
                return sum(1 for _ in f)
        except Exception:
            return 0
    
    def _parse_go_kfi_stats(self, kfi_stat_path):
        """GO-KFI統計の解析"""
        stats = {}
        
        try:
            with open(kfi_stat_path, 'r') as f:
                content = f.read()
                
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
                
        except Exception as e:
            self.logger.warning(f"GO-KFI統計解析エラー: {e}")
        
        return stats
    
    def generate_text_summary(self, results):
        """テキストサマリーの生成"""
        self.logger.info("テキストサマリーを生成中...")
        
        summary_path = os.path.join(self.output_dir, 'summary.txt')
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("tWhEEL - トランスクリプトーム機能解析パイプライン\n")
            f.write("実行サマリーレポート\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"実行日時: {results['timestamp']}\n\n")
            
            # 処理統計
            f.write("【処理統計】\n")
            f.write("-" * 30 + "\n")
            
            stats = results['statistics']
            
            # 遺伝子予測統計
            if 'gene_prediction' in stats:
                f.write(f"予測タンパク質数: {stats['gene_prediction'].get('protein_count', 'N/A')}\n")
            
            # アノテーション統計
            if 'annotation' in stats:
                f.write(f"Multi遺伝子数: {stats['annotation'].get('multi_genes', 'N/A')}\n")
                f.write(f"Single遺伝子数: {stats['annotation'].get('single_genes', 'N/A')}\n")
                f.write(f"総遺伝子数: {stats['annotation'].get('total_genes', 'N/A')}\n")
            
            # ガイド設計統計
            if 'guide_design' in stats:
                f.write(f"Multi遺伝子標的数: {stats['guide_design'].get('multi_targets', 'N/A')}\n")
                f.write(f"Single遺伝子標的数: {stats['guide_design'].get('single_targets', 'N/A')}\n")
            
            # GO解析統計
            if 'go_analysis' in stats:
                go_stats = stats['go_analysis']
                f.write(f"総GO term数: {go_stats.get('total_gos', 'N/A')}\n")
                f.write(f"Fully targetable GO割合: {go_stats.get('fully_targetable_percent', 'N/A')}%\n")
                f.write(f"Partially targetable GO割合: {go_stats.get('partially_targetable_percent', 'N/A')}%\n")
                f.write(f"Not targetable GO割合: {go_stats.get('not_targetable_percent', 'N/A')}%\n")
            
            f.write("\n")
            
            # 設定情報
            f.write("【設定情報】\n")
            f.write("-" * 30 + "\n")
            f.write(f"入力ファイル:\n")
            f.write(f"  - ゲノムアセンブリ: {results['config']['input'].get('genome_assembly', 'N/A')}\n")
            f.write(f"  - トランスクリプトーム: {results['config']['input'].get('transcriptome', 'N/A')}\n")
            f.write(f"  - UniProtデータベース: {results['config']['input'].get('uniprot_db', 'N/A')}\n")
            f.write(f"PAMパターン: {results['config']['crispr'].get('pam_pattern', 'N/A')}\n")
            f.write(f"最大ミスマッチ数: {results['config']['crispr'].get('max_mismatch', 'N/A')}\n")
            f.write(f"BUSCOリネージ: {results['config']['busco'].get('lineage', 'N/A')}\n")
            
            f.write("\n")
            
            # 出力ファイル
            f.write("【出力ファイル】\n")
            f.write("-" * 30 + "\n")
            f.write(f"HTMLレポート: {os.path.join(self.output_dir, 'report.html')}\n")
            f.write(f"Webデータ: {self.web_data_dir}\n")
            f.write(f"中間ファイル: {self.intermediate_dir}\n")
            
            f.write("\n")
            f.write("詳細な解析結果は HTMLレポート (report.html) をご覧ください。\n")
        
        self.logger.info(f"テキストサマリーを生成: {summary_path}")
        return summary_path
    
    def prepare_web_data(self, results):
        """WebアプリケーションのためのデータをJSONとSQLiteで準備"""
        self.logger.info("Webデータを準備中...")
        
        # SQLiteデータベースの作成
        db_path = os.path.join(self.web_data_dir, 'twheel_data.db')
        
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # テーブルの作成
        self._create_database_tables(cursor)
        
        # データの挿入
        self._insert_gene_data(cursor, results)
        self._insert_guide_data(cursor, results)
        self._insert_go_data(cursor, results)
        
        conn.commit()
        conn.close()
        
        # JSON形式でも基本統計を保存
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
        
        self.logger.info(f"Webデータを準備: {db_path}, {json_path}")
        return db_path, json_path
    
    def _create_database_tables(self, cursor):
        """データベーステーブルの作成"""
        
        # 遺伝子テーブル
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
        
        # ガイドRNAテーブル
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
        
        # GO termテーブル
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
        
        # 遺伝子-GO関連テーブル
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
        """遺伝子データの挿入"""
        
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
        """ガイドRNAデータの挿入"""
        
        final_table_path = results['data'].get('final_annotation_guides')
        if not final_table_path or not os.path.exists(final_table_path):
            return
        
        with open(final_table_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            for row in reader:
                if len(row) >= 11 and row[-1]:  # guide_targets列が存在
                    unique_id = row[9]
                    guide_targets = row[-1]
                    
                    # ガイド標的をパース
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
                    
                    # 遺伝子のガイド数を更新
                    cursor.execute("""
                        UPDATE genes SET guide_count = (
                            SELECT COUNT(*) FROM guide_rnas WHERE gene_unique_id = ?
                        ) WHERE unique_id = ?
                    """, (unique_id, unique_id))
    
    def _insert_go_data(self, cursor, results):
        """GO termデータの挿入"""
        
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
                    
                    # targetable_typeを判定
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
        """HTMLレポートの生成"""
        self.logger.info("HTMLレポートを生成中...")
        
        # テンプレートファイルの作成
        self._create_html_template()
        self._create_css_files()
        self._create_js_files()
        
        # メインHTMLファイルの生成
        html_path = os.path.join(self.output_dir, 'report.html')
        
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(self._generate_html_content(results))
        
        self.logger.info(f"HTMLレポートを生成: {html_path}")
        return html_path
    
    def _create_html_template(self):
        """HTMLテンプレートの作成"""
        template_dir = os.path.join('src', 'web', 'templates')
        Path(template_dir).mkdir(parents=True, exist_ok=True)
        
        # 基本テンプレートファイルは後で作成
        pass
    
    def _create_css_files(self):
        """CSSファイルの作成"""
        css_dir = os.path.join('src', 'web', 'static', 'css')
        
        # メインCSSファイル
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
        """JavaScriptファイルの作成"""
        js_dir = os.path.join('src', 'web', 'static', 'js')
        
        # メインJSファイル
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
            console.error('データの読み込みに失敗しました:', error);
            this.showError('データの読み込みに失敗しました');
        }
    }
    
    setupEventListeners() {
        // タブ切り替え
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.addEventListener('click', (e) => {
                const tabId = e.target.dataset.tab;
                this.showTab(tabId);
            });
        });
        
        // 検索機能
        document.addEventListener('input', (e) => {
            if (e.target.classList.contains('search-box')) {
                this.handleSearch(e.target.value, e.target.dataset.target);
            }
        });
        
        // フィルター機能
        document.addEventListener('change', (e) => {
            if (e.target.classList.contains('filter-select')) {
                this.handleFilter(e.target.value, e.target.dataset.filter);
            }
        });
    }
    
    showTab(tabId) {
        // タブの切り替え
        document.querySelectorAll('.nav-tab').forEach(tab => {
            tab.classList.remove('active');
        });
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.remove('active');
        });
        
        document.querySelector(`[data-tab="${tabId}"]`).classList.add('active');
        document.querySelector(`#${tabId}`).classList.add('active');
        
        this.currentTab = tabId;
        
        // タブ固有の初期化
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
        
        // 統計カードの更新
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
        // 遺伝子データの読み込み（実装予定）
        console.log('Loading gene data...');
    }
    
    async loadGuideData() {
        // ガイドRNAデータの読み込み（実装予定）
        console.log('Loading guide data...');
    }
    
    async loadGOData() {
        // GOデータの読み込み（実装予定）
        console.log('Loading GO data...');
    }
    
    handleSearch(query, target) {
        // 検索機能（実装予定）
        console.log('Search:', query, 'Target:', target);
    }
    
    handleFilter(value, filter) {
        // フィルター機能（実装予定）
        console.log('Filter:', filter, 'Value:', value);
    }
}

// アプリケーションの初期化
document.addEventListener('DOMContentLoaded', () => {
    new tWhEELReport();
});
""")
    
    def _generate_html_content(self, results):
        """HTMLコンテンツの生成"""
        stats = results['statistics']
        
        html_content = f"""
<!DOCTYPE html>
<html lang="ja">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - トランスクリプトーム機能解析レポート</title>
    <link rel="stylesheet" href="src/web/static/css/main.css">
</head>
<body>
    <div class="header">
        <div class="container">
            <h1>tWhEEL</h1>
            <p>トランスクリプトーム機能解析パイプライン</p>
            <p>実行日時: {results['timestamp']}</p>
        </div>
    </div>
    
    <div class="container">
        <div class="nav-tabs">
            <button class="nav-tab active" data-tab="dashboard">ダッシュボード</button>
            <button class="nav-tab" data-tab="genes">遺伝子ブラウザ</button>
            <button class="nav-tab" data-tab="guides">ガイドRNA</button>
            <button class="nav-tab" data-tab="go-analysis">GO解析</button>
            <button class="nav-tab" data-tab="settings">設定情報</button>
        </div>
        
        <!-- ダッシュボード -->
        <div id="dashboard" class="tab-content active">
            <h2>解析結果サマリー</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value" id="protein-count">
                        {stats.get('gene_prediction', {}).get('protein_count', 0):,}
                    </div>
                    <div class="stat-label">予測タンパク質数</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="total-genes">
                        {stats.get('annotation', {}).get('total_genes', 0):,}
                    </div>
                    <div class="stat-label">総遺伝子数</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="multi-genes">
                        {stats.get('annotation', {}).get('multi_genes', 0):,}
                    </div>
                    <div class="stat-label">Multi遺伝子数</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="single-genes">
                        {stats.get('annotation', {}).get('single_genes', 0):,}
                    </div>
                    <div class="stat-label">Single遺伝子数</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="total-gos">
                        {stats.get('go_analysis', {}).get('total_gos', 0):,}
                    </div>
                    <div class="stat-label">総GO term数</div>
                </div>
                
                <div class="stat-card">
                    <div class="stat-value" id="fully-targetable">
                        {stats.get('go_analysis', {}).get('fully_targetable_percent', 0):.1f}%
                    </div>
                    <div class="stat-label">Fully targetable GO</div>
                </div>
            </div>
            
            <div class="chart-container">
                <h3>GO Knock-out Feasibility分布</h3>
                <div id="go-kfi-chart">
                    <p class="loading">チャートを読み込み中...</p>
                </div>
            </div>
        </div>
        
        <!-- 遺伝子ブラウザ -->
        <div id="genes" class="tab-content">
            <h2>遺伝子ブラウザ</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="遺伝子名で検索..." data-target="genes">
                <select class="filter-select" data-filter="classification">
                    <option value="">全て</option>
                    <option value="Multi">Multi遺伝子</option>
                    <option value="Single">Single遺伝子</option>
                </select>
                <select class="filter-select" data-filter="guides">
                    <option value="">全て</option>
                    <option value="with-guides">ガイドRNAあり</option>
                    <option value="without-guides">ガイドRNAなし</option>
                </select>
            </div>
            
            <table class="data-table" id="genes-table">
                <thead>
                    <tr>
                        <th>遺伝子名</th>
                        <th>分類</th>
                        <th>UniProt ID</th>
                        <th>同一性</th>
                        <th>E値</th>
                        <th>ガイド数</th>
                        <th>説明</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="7" class="loading">データを読み込み中...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- ガイドRNA -->
        <div id="guides" class="tab-content">
            <h2>ガイドRNA設計結果</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="ガイド配列で検索..." data-target="guides">
                <select class="filter-select" data-filter="chromosome">
                    <option value="">全染色体</option>
                </select>
                <select class="filter-select" data-filter="strand">
                    <option value="">全ストランド</option>
                    <option value="+">+</option>
                    <option value="-">-</option>
                </select>
            </div>
            
            <table class="data-table" id="guides-table">
                <thead>
                    <tr>
                        <th>遺伝子</th>
                        <th>染色体</th>
                        <th>開始位置</th>
                        <th>終了位置</th>
                        <th>ガイド配列</th>
                        <th>ストランド</th>
                        <th>スコア</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="7" class="loading">データを読み込み中...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- GO解析 -->
        <div id="go-analysis" class="tab-content">
            <h2>GO機能解析結果</h2>
            
            <div class="filter-group">
                <input type="text" class="search-box" placeholder="GO termで検索..." data-target="go-terms">
                <select class="filter-select" data-filter="targetable">
                    <option value="">全て</option>
                    <option value="Fully">Fully targetable</option>
                    <option value="Partially">Partially targetable</option>
                    <option value="Not">Not targetable</option>
                </select>
            </div>
            
            <table class="data-table" id="go-table">
                <thead>
                    <tr>
                        <th>GO ID</th>
                        <th>GO名</th>
                        <th>ネームスペース</th>
                        <th>KFIスコア</th>
                        <th>ターゲタブル</th>
                        <th>遺伝子数</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td colspan="6" class="loading">データを読み込み中...</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <!-- 設定情報 -->
        <div id="settings" class="tab-content">
            <h2>設定情報</h2>
            
            <div class="setting-group">
                <h3>入力ファイル</h3>
                <ul>
                    <li><strong>ゲノムアセンブリ:</strong> {results['config']['input'].get('genome_assembly', 'N/A')}</li>
                    <li><strong>トランスクリプトーム:</strong> {results['config']['input'].get('transcriptome', 'N/A')}</li>
                    <li><strong>UniProtデータベース:</strong> {results['config']['input'].get('uniprot_db', 'N/A')}</li>
                </ul>
            </div>
            
            <div class="setting-group">
                <h3>CRISPR設定</h3>
                <ul>
                    <li><strong>PAMパターン:</strong> {results['config']['crispr'].get('pam_pattern', 'N/A')}</li>
                    <li><strong>最大ミスマッチ数:</strong> {results['config']['crispr'].get('max_mismatch', 'N/A')}</li>
                    <li><strong>スレッド数:</strong> {results['config']['crispr'].get('threads', 'N/A')}</li>
                </ul>
            </div>
            
            <div class="setting-group">
                <h3>BUSCO設定</h3>
                <ul>
                    <li><strong>リネージ:</strong> {results['config']['busco'].get('lineage', 'N/A')}</li>
                    <li><strong>CPU数:</strong> {results['config']['busco'].get('cpu', 'N/A')}</li>
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
        """レポート生成の全体実行"""
        try:
            # 全結果の収集
            results = self.collect_all_results()
            
            # テキストサマリーの生成
            if self.config['report']['generate_summary']:
                summary_path = self.generate_text_summary(results)
                results['summary_path'] = summary_path
            
            # Webデータの準備
            db_path, json_path = self.prepare_web_data(results)
            results['db_path'] = db_path
            results['json_path'] = json_path
            
            # HTMLレポートの生成
            if self.config['report']['generate_html']:
                html_path = self.generate_html_report(results)
                results['html_path'] = html_path
            
            # 結果をファイルに保存
            results_file = os.path.join(self.intermediate_dir, 'report_generation_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== レポート生成結果 ===\n\n")
                f.write(f"実行日時: {results['timestamp']}\n\n")
                
                if 'summary_path' in results:
                    f.write(f"テキストサマリー: {results['summary_path']}\n")
                
                if 'html_path' in results:
                    f.write(f"HTMLレポート: {results['html_path']}\n")
                
                f.write(f"SQLiteデータベース: {results['db_path']}\n")
                f.write(f"JSONサマリー: {results['json_path']}\n\n")
                
                f.write("=== 統計情報 ===\n")
                for category, stats in results['statistics'].items():
                    f.write(f"{category}:\n")
                    for key, value in stats.items():
                        f.write(f"  {key}: {value}\n")
                    f.write("\n")
            
            self.logger.info(f"レポート生成結果を保存: {results_file}")
            self.logger.info("レポート生成処理が完了しました")
            
            return results
            
        except Exception as e:
            self.logger.error(f"レポート生成処理でエラーが発生しました: {e}")
            raise
