#!/usr/bin/env python3
"""
GO Knock-out Feasibility Index (GO-KFI) 計算スクリプト
各GO termに対してノックアウトの実行可能性を評価
"""

import sys
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict


def parse_guide_targets(guide_targets_str):
    """ガイド標的文字列を解析してガイド数を返す"""
    if pd.isna(guide_targets_str) or guide_targets_str == '':
        return 0
    
    # カンマ区切りでガイド標的を分割
    guides = guide_targets_str.split(',')
    return len([g for g in guides if g.strip()])


def calculate_go_kfi(df):
    """GO-KFIを計算"""
    
    # GO termごとの遺伝子情報を集計
    go_stats = defaultdict(lambda: {
        'total_genes': 0,
        'genes_with_guides': 0,
        'total_guides': 0,
        'gene_list': []
    })
    
    # GO列を特定
    go_columns = [col for col in df.columns if col.startswith('GO_')]
    
    for index, row in df.iterrows():
        unique_id = row.get('unique_id', '')
        guide_targets = row.get('guide_targets', '')
        guide_count = parse_guide_targets(guide_targets)
        
        # 各GO列を処理
        for go_col in go_columns:
            go_terms_str = row.get(go_col, '')
            
            if pd.notna(go_terms_str) and go_terms_str != '':
                # セミコロン区切りでGO termを分割
                go_terms = [go.strip() for go in go_terms_str.split(';') if go.strip()]
                
                for go_term in go_terms:
                    go_stats[go_term]['total_genes'] += 1
                    go_stats[go_term]['total_guides'] += guide_count
                    go_stats[go_term]['gene_list'].append(unique_id)
                    
                    if guide_count > 0:
                        go_stats[go_term]['genes_with_guides'] += 1
    
    # GO-KFIを計算
    go_kfi_results = []
    
    for go_term, stats in go_stats.items():
        total_genes = stats['total_genes']
        genes_with_guides = stats['genes_with_guides']
        total_guides = stats['total_guides']
        
        # KFIスコアを計算
        if total_genes > 0:
            kfi_score = genes_with_guides / total_genes
        else:
            kfi_score = 0.0
        
        # ターゲタブル分類
        if kfi_score == 1.0:
            targetable_type = 'Fully'
        elif kfi_score > 0.0:
            targetable_type = 'Partially'
        else:
            targetable_type = 'Not'
        
        go_kfi_results.append({
            'GO_ID': go_term,
            'total_genes': total_genes,
            'genes_with_guides': genes_with_guides,
            'total_guides': total_guides,
            'kfi_score': kfi_score,
            'targetable_type': targetable_type
        })
    
    return pd.DataFrame(go_kfi_results)


def generate_statistics(kfi_df):
    """統計情報を生成"""
    
    total_gos = len(kfi_df)
    
    if total_gos == 0:
        return ""
    
    # 分類別のカウント
    fully_count = len(kfi_df[kfi_df['targetable_type'] == 'Fully'])
    partially_count = len(kfi_df[kfi_df['targetable_type'] == 'Partially'])
    not_count = len(kfi_df[kfi_df['targetable_type'] == 'Not'])
    
    # さらに詳細な分類
    fully_simple = len(kfi_df[(kfi_df['targetable_type'] == 'Fully') & (kfi_df['total_genes'] == 1)])
    fully_multi = fully_count - fully_simple
    
    # パーセンテージ計算
    fully_pct = (fully_count / total_gos) * 100
    fully_simple_pct = (fully_simple / total_gos) * 100
    fully_multi_pct = (fully_multi / total_gos) * 100
    partially_pct = (partially_count / total_gos) * 100
    not_pct = (not_count / total_gos) * 100
    
    # 統計レポート生成
    stats_report = f"""[DONE] {total_gos} GO terms → /DATA/go_kfi_table.csv
    ----------------------------------------------------
    |GO Knock-out Feasibility\t\t\t\t|
    ----------------------------------------------------
    |F:{fully_pct:.1f}%[S:{fully_simple_pct:.1f}%,M:{fully_multi_pct:.1f}%],P:{partially_pct:.1f}%,N:{not_pct:.1f}%,n:{total_gos}\t|
    |{fully_count}\tFully targetable GOs (F)\t\t|
    |{fully_simple}\tFully and simply targetable GOs (S)\t|
    |{fully_multi}\tFully and multiply targetable GOs (M)\t|
    |{partially_count}\tPartially targetable GOs (P)\t\t|
    |{not_count}\tNot targetable GOs (N)\t\t\t|
    |{total_gos}\tTotal GOs annotated\t\t\t|
    ---------------------------------------------------

[DONE] wrote /DATA/go_kfi_stat.csv
"""
    
    return stats_report


def main():
    parser = argparse.ArgumentParser(description='GO Knock-out Feasibility Index (GO-KFI) 計算')
    parser.add_argument('-i', '--input', required=True, help='入力TSVファイル')
    parser.add_argument('-o', '--output', required=True, help='出力KFIテーブルファイル')
    parser.add_argument('-s', '--stats', required=True, help='統計ファイル')
    
    args = parser.parse_args()
    
    print(f"Reading input file: {args.input}")
    
    # 入力ファイルを読み込み
    try:
        df = pd.read_csv(args.input, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(df)} rows for GO-KFI calculation...")
    
    # GO-KFIを計算
    kfi_df = calculate_go_kfi(df)
    
    # 結果を保存
    try:
        kfi_df.to_csv(args.output, index=False)
        print(f"GO-KFI table saved to: {args.output}")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 統計情報を生成
    stats_report = generate_statistics(kfi_df)
    
    # 統計情報を保存
    try:
        with open(args.stats, 'w') as f:
            f.write(stats_report)
        print(f"Statistics saved to: {args.stats}")
    except Exception as e:
        print(f"Error writing statistics file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 統計情報を標準出力にも表示
    print(stats_report)


if __name__ == "__main__":
    main()
