#!/usr/bin/env python3
"""
GO親term分析スクリプト
GO termの親分類を分析し、TreeMapで可視化
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import json


def get_go_parent_info():
    """GO親term情報を取得（簡易版）"""
    # 実際の実装では、GO OBOファイルを解析する必要があります
    # ここでは、主要なGO親termの例を示します
    
    go_parent_mapping = {
        # Biological Process
        'GO:0008150': 'biological_process',
        'GO:0009987': 'cellular_process',
        'GO:0050896': 'response_to_stimulus',
        'GO:0065007': 'biological_regulation',
        'GO:0008152': 'metabolic_process',
        'GO:0032501': 'multicellular_organismal_process',
        'GO:0032502': 'developmental_process',
        'GO:0040011': 'locomotion',
        'GO:0051179': 'localization',
        'GO:0000003': 'reproduction',
        
        # Molecular Function
        'GO:0003674': 'molecular_function',
        'GO:0003824': 'catalytic_activity',
        'GO:0005215': 'transporter_activity',
        'GO:0005488': 'binding',
        'GO:0038023': 'signaling_receptor_activity',
        'GO:0005198': 'structural_molecule_activity',
        'GO:0045182': 'translation_regulator_activity',
        'GO:0140110': 'transcription_regulator_activity',
        
        # Cellular Component
        'GO:0005575': 'cellular_component',
        'GO:0005623': 'cell',
        'GO:0032991': 'protein-containing_complex',
        'GO:0005576': 'extracellular_region',
        'GO:0005737': 'cytoplasm',
        'GO:0005634': 'nucleus',
        'GO:0005694': 'chromosome',
        'GO:0016020': 'membrane',
        'GO:0005856': 'cytoskeleton',
        'GO:0005773': 'vacuole'
    }
    
    return go_parent_mapping


def classify_go_terms(go_list):
    """GO termを親分類に分類"""
    parent_mapping = get_go_parent_info()
    
    # GO termを分類
    classified = defaultdict(list)
    unclassified = []
    
    for go_term in go_list:
        if not go_term.startswith('GO:'):
            continue
            
        # 直接マッピングがある場合
        if go_term in parent_mapping:
            classified[parent_mapping[go_term]].append(go_term)
        else:
            # 簡易的な分類（実際にはOBOファイルを使用）
            if go_term.startswith('GO:0008'):
                classified['biological_process'].append(go_term)
            elif go_term.startswith('GO:0003'):
                classified['molecular_function'].append(go_term)
            elif go_term.startswith('GO:0005'):
                classified['cellular_component'].append(go_term)
            else:
                unclassified.append(go_term)
    
    return classified, unclassified


def create_treemap(parent_counts, output_path):
    """TreeMapを作成"""
    try:
        import squarify
        
        # データの準備
        labels = []
        sizes = []
        
        for parent, count in parent_counts.items():
            labels.append(f"{parent}\n({count})")
            sizes.append(count)
        
        # カラーマップの設定
        colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
        
        # TreeMapの作成
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        squarify.plot(sizes=sizes, label=labels, color=colors, alpha=0.8, ax=ax)
        
        ax.set_title('GO Parent Terms Distribution', fontsize=16)
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"TreeMap saved to: {output_path}")
        
    except ImportError:
        print("Warning: squarify not available, creating simple bar chart instead")
        create_bar_chart(parent_counts, output_path)


def create_bar_chart(parent_counts, output_path):
    """バーチャートを作成（squarifyが使用できない場合）"""
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    parents = list(parent_counts.keys())
    counts = list(parent_counts.values())
    
    bars = ax.bar(parents, counts, color='skyblue', alpha=0.8)
    
    ax.set_title('GO Parent Terms Distribution', fontsize=16)
    ax.set_xlabel('GO Parent Terms', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    
    # バーの上に数値を表示
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}', ha='center', va='bottom')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Bar chart saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='GO親term分析とTreeMap可視化')
    parser.add_argument('--go_table', required=True, help='GO termテーブルファイル')
    parser.add_argument('--go_col', type=int, default=0, help='GO termの列番号')
    parser.add_argument('--treemap', action='store_true', help='TreeMapを作成')
    parser.add_argument('--out_csv', required=True, help='出力CSVファイル')
    parser.add_argument('--out_png', required=True, help='出力PNGファイル')
    
    args = parser.parse_args()
    
    print(f"Reading GO table: {args.go_table}")
    
    # GO tableを読み込み
    try:
        df = pd.read_csv(args.go_table, low_memory=False)
    except Exception as e:
        print(f"Error reading GO table: {e}", file=sys.stderr)
        sys.exit(1)
    
    # GO termの列を取得
    go_column = df.columns[args.go_col]
    go_terms = df[go_column].dropna().tolist()
    
    print(f"Processing {len(go_terms)} GO terms...")
    
    # GO termを分類
    classified, unclassified = classify_go_terms(go_terms)
    
    # 親termごとのカウント
    parent_counts = {}
    for parent, terms in classified.items():
        parent_counts[parent] = len(terms)
    
    # 未分類のGO termもカウント
    if unclassified:
        parent_counts['unclassified'] = len(unclassified)
    
    # 結果をDataFrameに変換
    result_df = pd.DataFrame([
        {'parent_term': parent, 'count': count}
        for parent, count in parent_counts.items()
    ])
    
    # CSVファイルに保存
    try:
        result_df.to_csv(args.out_csv, index=False)
        print(f"Parent counts saved to: {args.out_csv}")
    except Exception as e:
        print(f"Error writing CSV file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 可視化
    if args.treemap:
        create_treemap(parent_counts, args.out_png)
    else:
        create_bar_chart(parent_counts, args.out_png)
    
    # 結果をコンソールに出力
    print("\nGO Parent Term Distribution:")
    print("-" * 40)
    for parent, count in sorted(parent_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"{parent}: {count}")


if __name__ == "__main__":
    main()
