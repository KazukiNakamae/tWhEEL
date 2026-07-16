#!/usr/bin/env python3
"""
GO列追加スクリプト
UniProt IDからGO annotationを取得してTSVファイルに追加
"""

import sys
import csv
import requests
import time
import itertools
import re

# CSVフィールドサイズ制限を増加
csv.field_size_limit(1 << 30)

BATCH_SIZE = 100  # バッチサイズ
APIURL = "https://rest.uniprot.org/uniprotkb/stream"
GO_ID_PATTERN = re.compile(r'GO:\d{7}')


def extract_go_ids(go_field):
    """Extract unique GO IDs from a UniProt aspect field in source order."""
    if not go_field:
        return []
    return list(dict.fromkeys(GO_ID_PATTERN.findall(go_field)))


def classify_go_terms(go_terms_by_aspect):
    """Normalize GO IDs already classified by UniProt into BP, MF, and CC."""
    classified = {}
    all_go_ids = []
    seen = set()

    for column in ('GO_BP', 'GO_MF', 'GO_CC'):
        go_ids = list(dict.fromkeys(go_terms_by_aspect.get(column, [])))
        classified[column] = ';'.join(go_ids)
        for go_id in go_ids:
            if go_id not in seen:
                seen.add(go_id)
                all_go_ids.append(go_id)

    classified['GO_ALL'] = ';'.join(all_go_ids)
    return classified


def fetch_go_batch(uniprot_ids):
    """複数のUniProt IDからGO annotationをバッチで取得"""
    # 空文字と重複を除去
    uniprot_ids = sorted({uid for uid in uniprot_ids if uid})
    if not uniprot_ids:
        return {}
    
    # クエリを構築: accession:ID1 OR accession:ID2 ...
    query = " OR ".join(f"accession:{uid}" for uid in uniprot_ids)
    
    params = {
        "query": query,
        "fields": "accession,go_p,go_f,go_c",
        "format": "tsv"
    }
    
    headers = {
        'User-Agent': 'tWhEEL/1.0 (https://github.com/KazukiNakamae/tWhEEL) Python/3.x'
    }
    
    try:
        response = requests.get(APIURL, params=params, headers=headers, timeout=30)
        response.raise_for_status()
        
        result = {}
        reader = csv.reader(response.text.splitlines(), delimiter='\t')
        next(reader, None)
        for row in reader:
            if not row:
                continue

            padded_row = row + [''] * (4 - len(row))
            accession, go_bp, go_mf, go_cc = padded_row[:4]
            if not accession:
                continue

            result[accession] = classify_go_terms({
                'GO_BP': extract_go_ids(go_bp),
                'GO_MF': extract_go_ids(go_mf),
                'GO_CC': extract_go_ids(go_cc),
            })

        # Cache accessions with no returned annotations to avoid repeated requests.
        empty_annotation = classify_go_terms({})
        for accession in uniprot_ids:
            result.setdefault(accession, empty_annotation.copy())
        
        print(f"  Retrieved GO annotations for {len(result)}/{len(uniprot_ids)} UniProt IDs")
        return result
        
    except requests.exceptions.RequestException as e:
        print(f"  API error: {e}")
        return {}


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
    
    try:
        with open(input_file, 'r', newline='') as fin, \
             open(output_file, 'w', newline='') as fout:
            
            reader = csv.reader(fin, delimiter='\t')
            writer = csv.writer(fout, delimiter='\t')
            
            # ヘッダー行を処理
            try:
                header = next(reader)
                # GO列をヘッダーに追加
                new_header = header + ['GO_BP', 'GO_MF', 'GO_CC', 'GO_ALL']
                writer.writerow(new_header)
            except StopIteration:
                print("Error: Empty input file", file=sys.stderr)
                sys.exit(1)
            
            # グローバルキャッシュ
            go_cache = {}
            processed_rows = 0
            
            while True:
                # バッチサイズ分の行を読み込み
                chunk = list(itertools.islice(reader, BATCH_SIZE))
                if not chunk:
                    break
                
                # UniProt IDを抽出
                uniprot_ids = []
                for row in chunk:
                    if len(row) > 1:  # subject_idがある場合
                        uniprot_id = extract_uniprot_id(row[1])
                        uniprot_ids.append(uniprot_id)
                    else:
                        uniprot_ids.append('')
                
                # 新しいUniProt IDsのみAPIクエリ
                new_ids = [uid for uid in uniprot_ids if uid and uid not in go_cache]
                if new_ids:
                    batch_results = fetch_go_batch(new_ids)
                    # キャッシュに追加
                    for uid, go_data in batch_results.items():
                        go_cache[uid] = go_data
                    
                    # APIリクエスト間隔を調整
                    time.sleep(1)
                
                # 各行にGO情報を追加して出力
                for row, uniprot_id in zip(chunk, uniprot_ids):
                    if uniprot_id in go_cache:
                        go_data = go_cache[uniprot_id]
                        new_row = row + [
                            go_data['GO_BP'],
                            go_data['GO_MF'],
                            go_data['GO_CC'],
                            go_data['GO_ALL']
                        ]
                    else:
                        # GO情報がない場合は空文字
                        new_row = row + ['', '', '', '']
                    
                    writer.writerow(new_row)
                    processed_rows += 1
                
                if processed_rows % 100 == 0:
                    print(f"Processed {processed_rows} rows...")
            
            print(f"Successfully processed {processed_rows} rows and wrote to {output_file}")
            
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
