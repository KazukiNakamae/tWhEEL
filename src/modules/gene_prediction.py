#!/usr/bin/env python3
"""
遺伝子予測モジュール
TransDecoderを使用したORF予測とBUSCO評価
"""

import os
import subprocess
import logging
from pathlib import Path
import gzip
import shutil


class GenePrediction:
    """遺伝子予測処理クラス"""
    
    def __init__(self, config):
        """
        初期化
        
        Args:
            config (dict): 設定辞書
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        
    def prepare_input(self):
        """入力ファイルの準備"""
        transcriptome_path = self.config['input']['transcriptome']
        
        # 出力ファイルパス
        output_path = os.path.join(self.intermediate_dir, 'transcriptome.fa')
        
        # ファイルの解凍とコピー
        if transcriptome_path.endswith('.gz'):
            self.logger.info("トランスクリプトームファイルを解凍中...")
            with gzip.open(transcriptome_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("トランスクリプトームファイルをコピー中...")
            shutil.copy2(transcriptome_path, output_path)
        
        # 配列を1行にフォーマット
        self._format_fasta_oneline(output_path)
        
        self.logger.info(f"入力ファイルの準備が完了: {output_path}")
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
        self.logger.debug("FASTAファイルを1行形式に変換しました")
    
    def run_transdecoder(self, input_path):
        """TransDecoderによるORF予測"""
        self.logger.info("TransDecoderによるORF予測を開始...")
        
        # Docker設定
        docker_image = self.config['transdecoder']['docker_image']
        
        # TransDecoder.LongOrfs実行
        cmd_longorfs = [
            "docker", "run", 
            "-v", f"{os.getcwd()}:/data",
            "-w", "/data",
            "--rm",
            docker_image,
            "TransDecoder.LongOrfs",
            "-t", input_path
        ]
        
        self.logger.debug(f"実行コマンド: {' '.join(cmd_longorfs)}")
        result = subprocess.run(cmd_longorfs, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"TransDecoder.LongOrfs実行エラー: {result.stderr}")
            raise RuntimeError("TransDecoder.LongOrfs実行に失敗しました")
        
        # TransDecoder.Predict実行
        cmd_predict = [
            "docker", "run",
            "-v", f"{os.getcwd()}:/data",
            "-w", "/data", 
            "--rm",
            docker_image,
            "TransDecoder.Predict",
            "-t", input_path
        ]
        
        self.logger.debug(f"実行コマンド: {' '.join(cmd_predict)}")
        result = subprocess.run(cmd_predict, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"TransDecoder.Predict実行エラー: {result.stderr}")
            raise RuntimeError("TransDecoder.Predict実行に失敗しました")
        
        # 予測されたタンパク質配列ファイルをコピー
        # TransDecoderは入力ファイルのベースネームを使って出力ファイルを作成する
        input_filename = os.path.basename(input_path)
        protein_path = f"{input_filename}.transdecoder.pep"
        query_path = os.path.join(self.intermediate_dir, 'query.pep')
        
        if os.path.exists(protein_path):
            shutil.copy2(protein_path, query_path)
            self.logger.info(f"予測タンパク質配列を保存: {query_path}")
        else:
            raise FileNotFoundError(f"TransDecoder出力ファイルが見つかりません: {protein_path}")
        
        return query_path
    
    def get_sequence_stats(self, protein_path):
        """タンパク質配列の統計情報を取得"""
        self.logger.info("タンパク質配列の統計を取得中...")
        
        cmd = [
            "docker", "run",
            "-v", f"{os.getcwd()}:{os.getcwd()}",
            "-w", os.getcwd(),
            "--rm",
            "quay.io/biocontainers/seqkit:2.8.0--h9ee0642_0",
            "seqkit", "stats", protein_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.warning(f"seqkit統計取得エラー: {result.stderr}")
            return None
        
        self.logger.info("タンパク質配列統計:")
        self.logger.info(result.stdout)
        
        return result.stdout
    
    def run_busco(self, protein_path):
        """BUSCO評価の実行"""
        self.logger.info("BUSCO評価を開始...")
        
        # BUSCO出力ディレクトリ名（フルパスではなく名前のみ）
        busco_output_name = 'busco_output'
        busco_output_path = os.path.join(self.intermediate_dir, busco_output_name)
        
        # protein_pathの相対パスを取得
        protein_filename = os.path.basename(protein_path)
        
        cmd = [
            "docker", "run",
            "--platform", "linux/amd64",  # Apple Silicon対応
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/busco_work",
            "-w", "/busco_work",
            "--rm",
            self.config['busco']['docker_image'],
            "busco",
            "-m", "protein",
            "-i", protein_filename,  # 相対パス
            "-o", busco_output_name,  # ディレクトリ名のみ
            "-l", self.config['busco']['lineage'],
            "-c", str(self.config['busco']['cpu'])
        ]
        
        self.logger.debug(f"BUSCO実行コマンド: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"BUSCO実行エラー: {result.stderr}")
            raise RuntimeError("BUSCO評価に失敗しました")
        
        # BUSCO結果の読み取り
        busco_summary = self._parse_busco_results(busco_output_path)
        
        self.logger.info("BUSCO評価結果:")
        self.logger.info(busco_summary)
        
        return busco_summary
    
    def _parse_busco_results(self, busco_output):
        """BUSCO結果の解析"""
        summary_files = list(Path(busco_output).glob("**/short_summary*.txt"))
        
        if not summary_files:
            self.logger.warning("BUSCO短縮サマリーファイルが見つかりません")
            return "BUSCO結果を解析できませんでした"
        
        summary_file = summary_files[0]
        
        try:
            with open(summary_file, 'r') as f:
                content = f.read()
            return content
        except Exception as e:
            self.logger.error(f"BUSCO結果読み取りエラー: {e}")
            return "BUSCO結果を読み取れませんでした"
    
    def run(self):
        """遺伝子予測の全体実行"""
        try:
            # 入力ファイルの準備
            input_path = self.prepare_input()
            
            # TransDecoder実行
            protein_path = self.run_transdecoder(input_path)
            
            # 配列統計取得
            stats = self.get_sequence_stats(protein_path)
            
            # BUSCO評価
            busco_result = self.run_busco(protein_path)
            
            # 結果の保存
            results = {
                'input_path': input_path,
                'protein_path': protein_path,
                'sequence_stats': stats,
                'busco_result': busco_result
            }
            
            # 結果をファイルに保存
            results_file = os.path.join(self.intermediate_dir, 'gene_prediction_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== 遺伝子予測結果 ===\n\n")
                f.write(f"入力ファイル: {input_path}\n")
                f.write(f"予測タンパク質配列: {protein_path}\n\n")
                f.write("=== 配列統計 ===\n")
                f.write(stats if stats else "統計情報を取得できませんでした\n")
                f.write("\n=== BUSCO評価結果 ===\n")
                f.write(busco_result)
            
            self.logger.info(f"遺伝子予測結果を保存: {results_file}")
            self.logger.info("遺伝子予測処理が完了しました")
            
            return results
            
        except Exception as e:
            self.logger.error(f"遺伝子予測処理でエラーが発生しました: {e}")
            raise
