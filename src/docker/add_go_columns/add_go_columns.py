#!/usr/bin/env python3
"""
GO列追加スクリプト
UniProt IDからGO annotationを取得してTSVファイルに追加
"""

import sys
import pandas as pd
import requests
import time
from urllib.parse import quote


def fetch_go_annotations(uniprot_id):
    """UniProt IDからGO annotationを取得"""
    try:
        # UniProt REST APIを使用してGO情報を取得
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            
            # GO termを抽出
            go_terms = []
            if 'dbReferences' in data:
                for ref in data['dbReferences']:
                    if ref.get('type') == 'GO':
                        go_id = ref.get('id', '')
                        properties = ref.get('properties', {})
                        go_aspect = properties.get('GoAspect', '')
                        go_term = properties.get('GoTerm', '')
                        
                        if go_id:
                            go_terms.append({
                                'id': go_id,
                                'aspect': go_aspect,
                                'term': go_term
                            })
            
            # GO termを分類別に整理
            go_bp = []  # Biological Process
            go_mf = []  # Molecular Function
            go_cc = []  # Cellular Component
            
            for go in go_terms:
                if 'P:' in go['aspect']:
                    go_bp.append(go['id'])
                elif 'F:' in go['aspect']:
                    go_mf.append(go['id'])
                elif 'C:' in go['aspect']:
                    go_cc.append(go['id'])
            
            return {
                'GO_BP': ';'.join(go_bp) if go_bp else '',
                'GO_MF': ';'.join(go_mf) if go_mf else '',
                'GO_CC': ';'.join(go_cc) if go_cc else '',
                'GO_ALL': ';'.join([go['id'] for go in go_terms])
            }
            
    except Exception as e:
        print(f"Error fetching GO for {uniprot_id}: {e}", file=sys.stderr)
    
    return {
        'GO_BP': '',
        'GO_MF': '',
        'GO_CC': '',
        'GO_ALL': ''
    }


def extract_uniprot_id(subject_id):
    """subject_idからUniProt IDを抽出"""
    # UniProt IDの形式: sp|P12345|GENE_HUMAN
    if '|' in subject_id:
        parts = subject_id.split('|')
        if len(parts) >= 2:
            return parts[1]
    return subject_id


def main():
    if len(sys.argv) != 3:
        print("Usage: python add_go_columns.py input.tsv output.tsv", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"Reading input file: {input_file}")
    
    # TSVファイルを読み込み
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(df)} rows...")
    
    # GO列を初期化
    df['GO_BP'] = ''
    df['GO_MF'] = ''
    df['GO_CC'] = ''
    df['GO_ALL'] = ''
    
    # 処理済みUniProt IDのキャッシュ
    go_cache = {}
    
    # 各行を処理
    for index, row in df.iterrows():
        if index % 100 == 0:
            print(f"Processed {index}/{len(df)} rows...")
        
        # UniProt IDを取得
        subject_id = row.get('subject_id', '')
        uniprot_id = extract_uniprot_id(subject_id)
        
        if uniprot_id and uniprot_id not in go_cache:
            # GO annotationを取得
            go_annotations = fetch_go_annotations(uniprot_id)
            go_cache[uniprot_id] = go_annotations
            
            # APIリクエスト間隔を調整
            time.sleep(0.1)
        
        # GO情報を設定
        if uniprot_id in go_cache:
            go_data = go_cache[uniprot_id]
            df.at[index, 'GO_BP'] = go_data['GO_BP']
            df.at[index, 'GO_MF'] = go_data['GO_MF']
            df.at[index, 'GO_CC'] = go_data['GO_CC']
            df.at[index, 'GO_ALL'] = go_data['GO_ALL']
    
    print(f"Writing output file: {output_file}")
    
    # 結果を保存
    try:
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Successfully wrote {len(df)} rows to {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
