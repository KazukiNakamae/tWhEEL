#!/usr/bin/env python3
"""
tWhEEL - トランスクリプトーム機能解析パイプライン
Transcriptome Workflow for Hierarchical Exploration and Expression Levels

このスクリプトは、トランスクリプトアセンブリデータから遺伝子予測、機能アノテーション、
CRISPRガイド設計、GO機能解析を行い、ノックアウト実験のためのターゲット遺伝子を
特定するためのパイプラインです。
"""

import argparse
import sys
import os
import yaml
import logging
from pathlib import Path
from datetime import datetime

# モジュールのインポート
from modules.gene_prediction import GenePrediction
from modules.annotation import Annotation
from modules.guide_design import GuideDesign
from modules.go_analysis import GOAnalysis
from modules.report_generator import ReportGenerator


class tWhEELPipeline:
    """tWhEELパイプラインのメインクラス"""
    
    def __init__(self, config_path):
        """
        パイプラインの初期化
        
        Args:
            config_path (str): 設定ファイルのパス
        """
        self.config = self._load_config(config_path)
        self.setup_logging()
        self.setup_directories()
        
    def _load_config(self, config_path):
        """設定ファイルの読み込み"""
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            return config
        except Exception as e:
            print(f"設定ファイルの読み込みに失敗しました: {e}")
            sys.exit(1)
    
    def setup_logging(self):
        """ログ設定"""
        log_level = logging.DEBUG if self.config['execution']['verbose'] else logging.INFO
        
        # ログフォーマット
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # コンソールハンドラ
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        
        # ファイルハンドラ
        log_file = os.path.join(self.config['output']['base_dir'], 'twheel.log')
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        
        # ルートロガーの設定
        logging.basicConfig(
            level=log_level,
            handlers=[console_handler, file_handler]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info("tWhEELパイプラインを開始します")
    
    def setup_directories(self):
        """出力ディレクトリの作成"""
        directories = [
            self.config['output']['base_dir'],
            self.config['output']['intermediate'],
            os.path.join(self.config['output']['base_dir'], 'web_data')
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
            self.logger.debug(f"ディレクトリを作成: {directory}")
    
    def validate_inputs(self):
        """入力ファイルの検証"""
        required_files = [
            ('genome_assembly', self.config['input']['genome_assembly']),
            ('transcriptome', self.config['input']['transcriptome']),
            ('uniprot_db', self.config['input']['uniprot_db'])
        ]
        
        for name, filepath in required_files:
            if not filepath:
                self.logger.error(f"必須パラメータが設定されていません: {name}")
                return False
            if not os.path.exists(filepath):
                self.logger.error(f"ファイルが見つかりません: {filepath}")
                return False
                
        self.logger.info("入力ファイルの検証が完了しました")
        return True
    
    def run_step(self, step_name):
        """個別ステップの実行"""
        self.logger.info(f"ステップを開始: {step_name}")
        start_time = datetime.now()
        
        try:
            if step_name == "gene_prediction":
                predictor = GenePrediction(self.config)
                predictor.run()
                
            elif step_name == "annotation":
                annotator = Annotation(self.config)
                annotator.run()
                
            elif step_name == "guide_design":
                guide_designer = GuideDesign(self.config)
                guide_designer.run()
                
            elif step_name == "go_analysis":
                go_analyzer = GOAnalysis(self.config)
                go_analyzer.run()
                
            elif step_name == "report_generation":
                report_gen = ReportGenerator(self.config)
                report_gen.run()
                
            else:
                self.logger.error(f"不明なステップ名: {step_name}")
                return False
                
        except Exception as e:
            self.logger.error(f"ステップ{step_name}でエラーが発生しました: {e}")
            return False
        
        end_time = datetime.now()
        duration = end_time - start_time
        self.logger.info(f"ステップ{step_name}が完了しました (実行時間: {duration})")
        return True
    
    def run(self):
        """パイプライン全体の実行"""
        if not self.validate_inputs():
            return False
            
        start_time = datetime.now()
        self.logger.info("tWhEELパイプラインの実行を開始します")
        
        steps = self.config['execution']['steps']
        
        for step in steps:
            if not self.run_step(step):
                self.logger.error(f"パイプラインがステップ{step}で停止しました")
                return False
        
        end_time = datetime.now()
        total_duration = end_time - start_time
        self.logger.info(f"tWhEELパイプラインが正常に完了しました (総実行時間: {total_duration})")
        
        return True


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(
        description='tWhEEL - トランスクリプトーム機能解析パイプライン',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  python twheel.py --config config/default.yaml
  python twheel.py --config my_config.yaml --step annotation
  python twheel.py --list-steps
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        default='config/default.yaml',
        help='設定ファイルのパス (デフォルト: config/default.yaml)'
    )
    
    parser.add_argument(
        '--step', '-s',
        help='実行する特定のステップ (全ステップを実行する場合は指定不要)'
    )
    
    parser.add_argument(
        '--list-steps',
        action='store_true',
        help='利用可能なステップ一覧を表示'
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='tWhEEL version 1.0.0'
    )
    
    args = parser.parse_args()
    
    if args.list_steps:
        print("利用可能なステップ:")
        steps = ["gene_prediction", "annotation", "guide_design", "go_analysis", "report_generation"]
        for step in steps:
            print(f"  - {step}")
        return 0
    
    # パイプラインの実行
    try:
        pipeline = tWhEELPipeline(args.config)
        
        if args.step:
            # 特定のステップのみ実行
            success = pipeline.run_step(args.step)
        else:
            # 全ステップ実行
            success = pipeline.run()
            
        return 0 if success else 1
        
    except KeyboardInterrupt:
        print("\n実行がユーザーによって中断されました")
        return 1
    except Exception as e:
        print(f"予期しないエラーが発生しました: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
